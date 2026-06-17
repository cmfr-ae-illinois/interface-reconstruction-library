
#include <AMReX_DistributionMapping.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <fstream>
#include <sstream>

#include "examples/amrex_advector/cases.h"

using namespace amrex;

struct AmrCoreFill {
  AMREX_GPU_DEVICE
  void operator()(const amrex::IntVect& /*iv*/,
                  amrex::Array4<amrex::Real> const& /*data*/,
                  const int /*dcomp*/, const int /*numcomp*/,
                  amrex::GeometryData const& /*geom*/,
                  const amrex::Real /*time*/, const amrex::BCRec* /*bcr*/,
                  const int /*bcomp*/, const int /*orig_comp*/) const {
    // do something for external Dirichlet (BCType::ext_dir)
  }
};

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
AmrCoreAdv::AmrCoreAdv() {
  ReadParameters();

  // Geometry on all levels has been defined already.

  // No valid BoxArray and DistributionMapping have been defined.
  // But the arrays for them have been resized.

  int nlevs_max = max_level + 1;

  istep.resize(nlevs_max, 0);

  t_new.resize(nlevs_max, 0.0);
  t_old.resize(nlevs_max, -1.e100);
  dt.resize(nlevs_max, 1.e100);

  moments_new.resize(nlevs_max);
  moments_old.resize(nlevs_max);
  band_id.resize(nlevs_max);
  interface.resize(nlevs_max);

  facevel.resize(nlevs_max);

  int bc_lo[AMREX_SPACEDIM];
  int bc_hi[AMREX_SPACEDIM];

  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    if (Geom(0).isPeriodic()[idim] == 1) {
      bc_lo[idim] = bc_hi[idim] = BCType::int_dir;  // periodic
    } else {
      bc_lo[idim] = bc_hi[idim] = BCType::foextrap;  // walls (Neumann)
    }
  }

  bcs.resize(1);  // Setup 1-component
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    bcs[0].setLo(idim, BCType::int_dir);
    bcs[0].setHi(idim, BCType::int_dir);
  }

  // stores fluxes at coarse-fine interface for synchronization
  // this will be sized "nlevs_max+1"
  // NOTE: the flux register associated with flux_reg[lev] is associated
  // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
  // therefore flux_reg[0] is never actually used in the reflux operation
  flux_reg.resize(nlevs_max + 1);

  ParmParse ppcase("case");
  ppcase.get("name", case_name);

  ParmParse pprec("reconstruction");
  pprec.get("name", reconstruction_name);

  ParmParse ppadv("advection");
  ppadv.get("name", advection_name);
}

AmrCoreAdv::~AmrCoreAdv() {}

// advance solution to final time
void AmrCoreAdv::Evolve() {
  Real cur_time = t_new[0];
  int last_plot_file_step = 0;

  for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step) {
    amrex::Print() << "\nCoarse STEP " << step + 1 << " starts ..."
                   << std::endl;

    ComputeDt();

    int lev = 0;
    int iteration = 1;
    timeStepNoSubcycling(cur_time, iteration);

    cur_time += dt[0];

    // sum phi to check conservation
    Real sum_phi = moments_new[0].sum();

    amrex::Print() << "Coarse STEP " << step + 1 << " ends."
                   << " TIME = " << cur_time << " DT = " << dt[0]
                   << " Sum(Phi) = " << sum_phi << std::endl;

    // sync up time
    for (lev = 0; lev <= finest_level; ++lev) {
      t_new[lev] = cur_time;
    }

    if (plot_int > 0 && (step + 1) % plot_int == 0) {
      last_plot_file_step = step + 1;
      WritePlotFile();
    }

    if (chk_int > 0 && (step + 1) % chk_int == 0) {
      WriteCheckpointFile();
    }

    if (cur_time >= stop_time - 1.e-6 * dt[0]) break;
  }

  {
    // amrex::MultiFab uniform_final;
    // BuildUniformFinestMoments(uniform_final);
    // const std::string final_filename = UniformMomentsBinaryFileName("final");
    // WriteUniformMomentsBinary(uniform_final, final_filename);
    // amrex::Print() << "Initial uniform domain = "
    //                << uniform_initial_moments.boxArray().minimalBox() <<
    //                "\n";

    // amrex::Print() << "Final uniform domain   = "
    //                << uniform_final.boxArray().minimalBox() << "\n";

    // amrex::Print() << "Initial nComp = " << uniform_initial_moments.nComp()
    //                << "\n";

    // amrex::Print() << "Final nComp   = " << uniform_final.nComp() << "\n";
    // ComputeUniformMomentL1Errors(uniform_initial_moments, uniform_final);
  }

  {
    const amrex::Real final_liquid_mass = ComputeCompositeM0();
    const amrex::Real change_liquid_mass =
        final_liquid_mass - initial_liquid_mass;

    const amrex::Real relative_change_liquid_mass =
        change_liquid_mass / initial_liquid_mass;

    amrex::Print() << "\nLiquid mass summary from AMR hierarchy\n";
    amrex::Print() << "  Initial M0        = " << std::setprecision(17)
                   << initial_liquid_mass << "\n";
    amrex::Print() << "  Final M0          = " << std::setprecision(17)
                   << final_liquid_mass << "\n";
    amrex::Print() << "  Change M0         = " << std::setprecision(17)
                   << change_liquid_mass << "\n";
    amrex::Print() << "  Relative change   = " << std::setprecision(17)
                   << relative_change_liquid_mass << "\n";
  }

  // new error norm
  const amrex::Real l1_error_m0 = ComputeL1ErrorM0();
  amrex::Print() << std::scientific << std::setprecision(16)
                 << "L1   M0 = " << l1_error_m0 << std::endl;

  if (plot_int > 0 && istep[0] > last_plot_file_step) {
    GetReconstruction(finest_level);
    WritePlotFile();
  }
}

// initializes multilevel data
void AmrCoreAdv::InitData() {
  if (restart_chkfile == "") {
    // start simulation from the beginning
    const Real time = 0.0;
    InitFromScratch(time);
    AverageDown();

    // amrex::Print() << "  Number of levels = " << finest_level + 1 << "\n";
    // for (int lev = 0; lev <= finest_level; ++lev) {
    //   amrex::Print() << "  Level " << lev
    //                  << " : number of boxes = " << grids[lev].size() << "\n";
    // }

    initial_liquid_mass = ComputeCompositeM0();
    initial_moments.define(moments_new[0].boxArray(),
                           moments_new[0].DistributionMap(),
                           moments_new[0].nComp(), 0);
    MultiFab::Copy(initial_moments, moments_new[0], 0, 0,
                   moments_new[0].nComp(), 0);
    {
      // amrex::MultiFab uniform_initial;
      // BuildUniformFinestMoments(uniform_initial_moments);
      // const std::string initial_filename =
      //     UniformMomentsBinaryFileName("initial");
      // WriteUniformMomentsBinary(uniform_initial_moments, initial_filename);
    }

    if (chk_int > 0) {
      WriteCheckpointFile();
    }

  } else {
    // restart from a checkpoint
    ReadCheckpointFile();
  }

  if (plot_int > 0) {
    WritePlotFile();
  }
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void AmrCoreAdv::MakeNewLevelFromCoarse(int lev, Real time, const BoxArray& ba,
                                        const DistributionMapping& dm) {
  const int ncomp = moments_new[lev - 1].nComp();
  const int nghost = moments_new[lev - 1].nGrow();

  moments_new[lev].define(ba, dm, ncomp, nghost);
  moments_old[lev].define(ba, dm, ncomp, nghost);
  band_id[lev].define(ba, dm, 1, 1);
  interface[lev].define(ba, dm, 1, nghost);

  t_new[lev] = time;
  t_old[lev] = time - 1.e200;

  // This clears the old MultiFab and allocates the new one
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    facevel[lev][idim] = MultiFab(
        amrex::convert(ba, IntVect::TheDimensionVector(idim)), dm, 1, num_grow);
  }

  if (lev > 0 && do_reflux) {
    flux_reg[lev].reset(
        new FluxRegister(ba, dm, refRatio(lev - 1), lev, ncomp));
  }

  FillCoarsePatch(lev, time, moments_new[lev], 0, ncomp);
}

DistributionMapping AmrCoreAdv::MakeDistributionMapWithWeights(
    int lev, const BoxArray& ba) {
  // Build weights from band_id on the OLD layout
  // one weight per box in the OLD grids[lev]
  int nboxes = ba.size();
  Vector<Real> weights(nboxes, 1.0);

  for (MFIter mfi(band_id[lev]); mfi.isValid(); ++mfi) {
    // sum band_id over the box to get a scalar cost per box
    auto const& arr = band_id[lev].const_array(mfi);
    const Box& bx = mfi.validbox();
    Real local_sum = 0;
    Loop(bx, [&](int i, int j, int k) { local_sum += arr(i, j, k); });
    weights[mfi.index()] = std::max(1.0, local_sum);
  }

  // Reduce across all MPI ranks so every rank has the full weight vector
  ParallelDescriptor::ReduceRealSum(weights.data(), nboxes);

  return DistributionMapping::makeKnapSack(weights);
}

void AmrCoreAdv::RedistributeLevel(int lev) {
  // Build new DM from current weights
  const BoxArray& ba = grids[lev];
  DistributionMapping new_dm = MakeDistributionMapWithWeights(lev, ba);

  // If DM hasn't changed, nothing to do
  if (new_dm == dmap[lev]) return;

  // --- Rebuild new_moments_new ---
  {
    MultiFab new_moments_new(ba, new_dm, moments_new[lev].nComp(),
                             moments_new[lev].nGrow());
    new_moments_new.ParallelCopy(
        moments_new[lev], 0, 0, moments_new[lev].nComp(),
        moments_new[lev].nGrow(), moments_new[lev].nGrow());
    std::swap(moments_new[lev], new_moments_new);
  }

  // --- Rebuild phi_old ---
  {
    MultiFab new_moments_old(ba, new_dm, moments_old[lev].nComp(),
                             moments_old[lev].nGrow());
    new_moments_old.ParallelCopy(
        moments_old[lev], 0, 0, moments_old[lev].nComp(),
        moments_old[lev].nGrow(), moments_old[lev].nGrow());
    std::swap(moments_old[lev], new_moments_old);
  }

  // --- Rebuild band_id ---
  {
    MultiFab new_band_id(ba, new_dm, 1, band_id[lev].nGrow());
    new_band_id.ParallelCopy(band_id[lev], 0, 0, 1, band_id[lev].nGrow(),
                             band_id[lev].nGrow());
    std::swap(band_id[lev], new_band_id);
  }

  // --- Rebuild interface ---
  {
    SepUnionMultiFab new_interface(ba, new_dm, 1, interface[lev].nGrow());
    new_interface.ParallelCopy(interface[lev], 0, 0, 1, interface[lev].nGrow(),
                               interface[lev].nGrow());
    std::swap(interface[lev], new_interface);
  }

  // --- Rebuild flux registers if needed ---
  // This clears the old MultiFab and allocates the new one
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    facevel[lev][idim] =
        MultiFab(amrex::convert(ba, IntVect::TheDimensionVector(idim)), new_dm,
                 1, num_grow);
  }

  if (lev > 0 && do_reflux) {
    flux_reg[lev].reset(new FluxRegister(ba, new_dm, refRatio(lev - 1), lev,
                                         moments_new[lev].nComp()));
  }

  // Tell AMReX about the new DM
  SetDistributionMap(lev, new_dm);
}

void AmrCoreAdv::regrid(int lbase, Real time, bool) {
  if (lbase >= max_level) {
    return;
  }

  int new_finest;
  Vector<BoxArray> new_grids(finest_level + 2);
  MakeNewGrids(lbase, time, new_finest, new_grids);

  BL_ASSERT(new_finest <= finest_level + 1);

  bool coarse_ba_changed = false;
  for (int lev = lbase + 1; lev <= new_finest; ++lev) {
    if (lev <= finest_level)  // an old level
    {
      bool ba_changed = (new_grids[lev] != grids[lev]);
      if (ba_changed || coarse_ba_changed) {
        BoxArray level_grids = grids[lev];
        DistributionMapping level_dmap = dmap[lev];
        if (ba_changed) {
          level_grids = new_grids[lev];
          level_dmap = MakeDistributionMap(lev, level_grids);
          // level_dmap = MakeDistributionMapWithWeights(lev, level_grids);
        }
        const auto old_num_setdm = num_setdm;
        RemakeLevel(lev, time, level_grids, level_dmap);
        SetBoxArray(lev, level_grids);
        if (old_num_setdm == num_setdm) {
          SetDistributionMap(lev, level_dmap);
        }
      }
      coarse_ba_changed = ba_changed;
      ;
    } else  // a new level
    {
      DistributionMapping new_dmap = MakeDistributionMap(lev, new_grids[lev]);
      // DistributionMapping new_dmap =
      //     MakeDistributionMapWithWeights(lev, new_grids[lev]);
      const auto old_num_setdm = num_setdm;
      MakeNewLevelFromCoarse(lev, time, new_grids[lev], new_dmap);
      SetBoxArray(lev, new_grids[lev]);
      if (old_num_setdm == num_setdm) {
        SetDistributionMap(lev, new_dmap);
      }
    }
  }

  for (int lev = new_finest + 1; lev <= finest_level; ++lev) {
    ClearLevel(lev);
    ClearBoxArray(lev);
    ClearDistributionMap(lev);
  }

  finest_level = new_finest;
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void AmrCoreAdv::RemakeLevel(int lev, Real time, const BoxArray& ba,
                             const DistributionMapping& dm) {
  const int ncomp = moments_new[lev].nComp();
  const int nghost = moments_new[lev].nGrow();

  MultiFab new_state(ba, dm, ncomp, nghost);
  MultiFab old_state(ba, dm, ncomp, nghost);
  MultiFab old_band_id(ba, dm, 1, 1);
  SepUnionMultiFab new_interface(ba, dm, 1, nghost);

  FillPatch(lev, time, new_state, 0, ncomp);
  new_interface.ParallelCopyWithPeriodicShift(interface[lev], 0, 0, 1, nghost,
                                              nghost, geom[lev]);

  std::swap(new_state, moments_new[lev]);
  std::swap(old_state, moments_old[lev]);
  std::swap(old_band_id, band_id[lev]);
  std::swap(new_interface, interface[lev]);

  t_new[lev] = time;
  t_old[lev] = time - 1.e200;

  // This clears the old MultiFab and allocates the new one
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    facevel[lev][idim] = MultiFab(
        amrex::convert(ba, IntVect::TheDimensionVector(idim)), dm, 1, num_grow);
  }

  if (lev > 0 && do_reflux) {
    flux_reg[lev].reset(
        new FluxRegister(ba, dm, refRatio(lev - 1), lev, ncomp));
  }
}

// Delete level data
// overrides the pure virtual function in AmrCore
void AmrCoreAdv::ClearLevel(int lev) {
  moments_new[lev].clear();
  moments_old[lev].clear();
  band_id[lev].clear();
  interface[lev].clear();
  flux_reg[lev].reset(nullptr);
}

// Make a new level from scratch using provided BoxArray and
// DistributionMapping. Only used during initialization. overrides the pure
// virtual function in AmrCore
void AmrCoreAdv::MakeNewLevelFromScratch(int lev, Real time, const BoxArray& ba,
                                         const DistributionMapping& dm) {
  const int ncomp = ncomp_moments;
  const int nghost = 0;

  moments_new[lev].define(ba, dm, ncomp, nghost);
  moments_old[lev].define(ba, dm, ncomp, nghost);
  band_id[lev].define(ba, dm, 1, 1);
  band_id[lev].setVal(0.0);
  interface[lev].define(ba, dm, 1, nghost);

  t_new[lev] = time;
  t_old[lev] = time - 1.e200;

  // This clears the old MultiFab and allocates the new one
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    facevel[lev][idim] = MultiFab(
        amrex::convert(ba, IntVect::TheDimensionVector(idim)), dm, 1, num_grow);
  }

  if (lev > 0 && do_reflux) {
    flux_reg[lev].reset(
        new FluxRegister(ba, dm, refRatio(lev - 1), lev, ncomp));
  }

  MultiFab& moments_lev = moments_new[lev];
  SepUnionMultiFab& interface_lev = interface[lev];

  const auto problo = Geom(lev).ProbLoArray();
  const auto dx = Geom(lev).CellSizeArray();

  // int num_tiles_local = 0;

  // for (MFIter mfi(moments_lev, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
  //   ++num_tiles_local;
  // }

  // int num_tiles_global = num_tiles_local;
  // amrex::ParallelDescriptor::ReduceIntSum(
  //     num_tiles_global, amrex::ParallelDescriptor::IOProcessorNumber());

  // amrex::Print() << "Level " << lev << "\n";
  // amrex::Print() << "  number of boxes = " << moments_lev.boxArray().size()
  //                << "\n";
  // amrex::Print() << "  number of tiles = " << num_tiles_global << "\n";

  for (MFIter mfi(moments_lev, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    Array4<Real> moments_fab = moments_lev[mfi].array();
    Array4<IRL::SeparatorUnion> interface_fab = interface_lev[mfi].array();
    const Box& box = mfi.tilebox();

    if (case_name == "deformation3d") {
      amrex::launch(box, [=] AMREX_GPU_DEVICE(const Box& tbx) {
        Deformation3D::initialize_case(tbx, moments_fab, interface_fab, problo,
                                       dx);
      });
    } else if (case_name == "translation3d" || case_name == "default") {
      amrex::launch(box, [=] AMREX_GPU_DEVICE(const Box& tbx) {
        Translation3D::initialize_case(tbx, moments_fab, interface_fab, problo,
                                       dx);
      });
    } else if (case_name == "rotation3d") {
      amrex::launch(box, [=] AMREX_GPU_DEVICE(const Box& tbx) {
        Rotation3D::initialize_case(tbx, moments_fab, interface_fab, problo,
                                    dx);
      });
    } else {
      throw std::runtime_error("Unknown case");
    }
  }

  UpdateBand();
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void AmrCoreAdv::ErrorEst(int lev, TagBoxArray& tags, Real /*time*/,
                          int /*ngrow*/) {
  static bool first = true;

  // only do this during the first call to ErrorEst
  if (first) {
    first = false;
    // read in an array of "phierr", which is the tagging threshold
    // in this example, we tag values of "phi" which are greater than phierr
    // for that particular level
    // in subroutine state_error, you could use more elaborate tagging, such
    // as more advanced logical expressions, or gradients, etc.
    // ParmParse pp("adv");
    // int n = pp.countval("phierr");
    // if (n > 0) {
    //   pp.getarr("phierr", phierr, 0, n);
    // }
  }

  //   if (lev >= phierr.size()) return;

  const int clearval = TagBox::CLEAR;
  const int tagval = TagBox::SET;

  const MultiFab& state = band_id[lev];
  GeometryData geomdata = geom[lev].data();
  const Real* AMREX_RESTRICT dx = geomdata.CellSize();
  const Real vol_inv = 1.0 / (dx[0] * dx[1] * dx[2]);

  for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const auto statefab = state.array(mfi);
    const auto tagfab = tags.array(mfi);
    //   Real phierror = phierr[lev];

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // state_error(i, j, k, tagfab, statefab, phierror,
      // tagval);
      if (statefab(i, j, k) != 0.0) {
        tagfab(i, j, k) = tagval;
        // } else {
        //   tagfab(i, j, k) = clearval;
      }
    });
  }
}

// read in some parameters from inputs file
void AmrCoreAdv::ReadParameters() {
  {
    ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
    pp.query("max_step", max_step);
    pp.query("stop_time", stop_time);
  }

  {
    ParmParse pp("amr");  // Traditionally, these have prefix, amr.

    pp.query("regrid_int", regrid_int);
    pp.query("plot_file", plot_file);
    pp.query("plot_int", plot_int);
    pp.query("chk_file", chk_file);
    pp.query("chk_int", chk_int);
    pp.query("restart", restart_chkfile);
  }

  {
    ParmParse pp("adv");

    pp.query("cfl", cfl);
    num_grow = std::max(2, 1 + static_cast<int>(std::ceil(cfl)));
    amrex::Print() << "Target CFL = " << cfl << ", requiring " << num_grow
                   << " ghost layers\n";
    pp.query("do_reflux", do_reflux);
  }
}

// set covered coarse cells to be the average of overlying fine cells
void AmrCoreAdv::AverageDown() {
  for (int lev = finest_level - 1; lev >= 0; --lev) {
    amrex::average_down(moments_new[lev + 1], moments_new[lev], geom[lev + 1],
                        geom[lev], 0, moments_new[lev].nComp(), refRatio(lev));
  }
}

// more flexible version of AverageDown() that lets you average down across
// multiple levels
void AmrCoreAdv::AverageDownTo(int crse_lev) {
  amrex::average_down(moments_new[crse_lev + 1], moments_new[crse_lev],
                      geom[crse_lev + 1], geom[crse_lev], 0,
                      moments_new[crse_lev].nComp(), refRatio(crse_lev));
}

// compute a new multifab by coping in phi from valid region and filling ghost
// cells works for single level and 2-level cases (fill fine grid ghost by
// interpolating from coarse)
template <typename MF>
void AmrCoreAdv::FillPatch(int lev, Real time, MF& mf, int icomp, int ncomp) {
  if (lev == 0) {
    Vector<MF*> smf;
    Vector<Real> stime;
    GetData(0, time, smf, stime);

    if (Gpu::inLaunchRegion()) {
      GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
      PhysBCFunct<GpuBndryFuncFab<AmrCoreFill>> physbc(geom[lev], bcs,
                                                       gpu_bndry_func);
      amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                  geom[lev], physbc, 0);
    } else {
      CpuBndryFuncFab bndry_func(
          nullptr);  // Without EXT_DIR, we can pass a nullptr.
      PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bndry_func);
      amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                  geom[lev], physbc, 0);
    }
  } else {
    Vector<MF*> cmf, fmf;
    Vector<Real> ctime, ftime;
    GetData(lev - 1, time, cmf, ctime);
    GetData(lev, time, fmf, ftime);

    Interpolater* mapper = &cell_cons_interp;

    if (Gpu::inLaunchRegion()) {
      GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
      PhysBCFunct<GpuBndryFuncFab<AmrCoreFill>> cphysbc(geom[lev - 1], bcs,
                                                        gpu_bndry_func);
      PhysBCFunct<GpuBndryFuncFab<AmrCoreFill>> fphysbc(geom[lev], bcs,
                                                        gpu_bndry_func);

      amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime, 0, icomp,
                                ncomp, geom[lev - 1], geom[lev], cphysbc, 0,
                                fphysbc, 0, refRatio(lev - 1), mapper, bcs, 0);
    } else {
      CpuBndryFuncFab bndry_func(
          nullptr);  // Without EXT_DIR, we can pass a nullptr.
      PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev - 1], bcs, bndry_func);
      PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev], bcs, bndry_func);

      amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime, 0, icomp,
                                ncomp, geom[lev - 1], geom[lev], cphysbc, 0,
                                fphysbc, 0, refRatio(lev - 1), mapper, bcs, 0);
    }
  }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
template <typename MF>
void AmrCoreAdv::FillCoarsePatch(int lev, Real time, MF& mf, int icomp,
                                 int ncomp) {
  BL_ASSERT(lev > 0);

  Vector<MultiFab*> cmf;
  Vector<Real> ctime;
  GetData(lev - 1, time, cmf, ctime);
  Interpolater* mapper = &cell_cons_interp;

  if (cmf.size() != 1) {
    amrex::Abort("FillCoarsePatch: how did this happen?");
  }

  if (Gpu::inLaunchRegion()) {
    GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
    PhysBCFunct<GpuBndryFuncFab<AmrCoreFill>> cphysbc(geom[lev - 1], bcs,
                                                      gpu_bndry_func);
    PhysBCFunct<GpuBndryFuncFab<AmrCoreFill>> fphysbc(geom[lev], bcs,
                                                      gpu_bndry_func);

    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp,
                                 geom[lev - 1], geom[lev], cphysbc, 0, fphysbc,
                                 0, refRatio(lev - 1), mapper, bcs, 0);
  } else {
    CpuBndryFuncFab bndry_func(
        nullptr);  // Without EXT_DIR, we can pass a nullptr.
    PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev - 1], bcs, bndry_func);
    PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev], bcs, bndry_func);

    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp,
                                 geom[lev - 1], geom[lev], cphysbc, 0, fphysbc,
                                 0, refRatio(lev - 1), mapper, bcs, 0);
  }
}

// utility to copy in data from moments_old and/or moments_new into another
// multifab
template <>
void AmrCoreAdv::GetData(int lev, Real time, Vector<MultiFab*>& data,
                         Vector<Real>& datatime) {
  data.clear();
  datatime.clear();

  const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

  if (time > t_new[lev] - teps && time < t_new[lev] + teps) {
    data.push_back(&moments_new[lev]);
    datatime.push_back(t_new[lev]);
  } else if (time > t_old[lev] - teps && time < t_old[lev] + teps) {
    data.push_back(&moments_old[lev]);
    datatime.push_back(t_old[lev]);
  } else {
    data.push_back(&moments_old[lev]);
    data.push_back(&moments_new[lev]);
    datatime.push_back(t_old[lev]);
    datatime.push_back(t_new[lev]);
  }
}

// utility to copy in data from interfaces into another
// multifab
template <>
void AmrCoreAdv::GetData(int lev, Real time, Vector<SepUnionMultiFab*>& data,
                         Vector<Real>& datatime) {
  data.clear();
  datatime.clear();
  data.push_back(&interface[lev]);
  datatime.push_back(t_new[lev]);
}

// Advance all the levels with the same dt
void AmrCoreAdv::timeStepNoSubcycling(Real time, int iteration) {
  if (max_level > 0 && regrid_int > 0)  // We may need to regrid
  {
    if (istep[0] % regrid_int == 0) {
      UpdateBand();
      regrid(0, time);
      // RedistributeLevel(finest_level);
    }
  }

  if (Verbose()) {
    for (int lev = 0; lev <= finest_level; lev++) {
      amrex::Print() << "[Level " << lev << " step " << istep[lev] + 1 << "] ";
      amrex::Print() << "ADVANCE with time = " << t_new[lev]
                     << " dt = " << dt[0] << std::endl;
    }
  }

  DefineVelocityAllLevels(time + 0.5_rt * dt[0]);
  AdvanceAllLevels(time, dt[0], iteration);

  // Make sure the coarser levels are consistent with the finer levels
  AverageDown();

  for (int lev = 0; lev <= finest_level; lev++) ++istep[lev];

  if (Verbose()) {
    for (int lev = 0; lev <= finest_level; lev++) {
      amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
      amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
    }
  }
}

// a wrapper for EstTimeStep
void AmrCoreAdv::ComputeDt() {
  Vector<Real> dt_tmp(finest_level + 1);

  for (int lev = 0; lev <= finest_level; ++lev) {
    dt_tmp[lev] = EstTimeStep(lev, t_new[lev]);
  }
  ParallelDescriptor::ReduceRealMin(&dt_tmp[0], dt_tmp.size());

  constexpr Real change_max = 1.1;
  Real dt_0 = dt_tmp[0];

  for (int lev = 0; lev <= finest_level; ++lev) {
    dt_tmp[lev] = std::min(dt_tmp[lev], change_max * dt[lev]);
    dt_0 = std::min(dt_0, dt_tmp[lev]);
  }

  // Limit dt's by the value of stop_time.
  const Real eps = 1.e-3 * dt_0;

  if (t_new[0] + dt_0 > stop_time - eps) {
    dt_0 = stop_time - t_new[0];
  }

  dt[0] = dt_0;

  for (int lev = 1; lev <= finest_level; ++lev) {
    dt[lev] = dt[lev - 1];
  }
}

// compute dt from CFL considerations
Real AmrCoreAdv::EstTimeStep(int lev, Real time) {
  BL_PROFILE("AmrCoreAdv::EstTimeStep()");

  Real dt_est = std::numeric_limits<Real>::max();

  const Real* dx = geom[lev].CellSize();

  Real max_vel = 0.0;
  if (case_name == "deformation3d") {
    max_vel = Deformation3D::get_max_vel();
  } else if (case_name == "translation3d") {
    max_vel = Translation3D::get_max_vel();
  } else if (case_name == "rotation3d") {
    max_vel = Rotation3D::get_max_vel();
  }

  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    dt_est = amrex::min(dt_est, dx[idim] / max_vel);
  }

  dt_est *= cfl;

  return dt_est;
}

// get plotfile name
std::string AmrCoreAdv::PlotFileName(int lev) const {
  return amrex::Concatenate(plot_file, lev, 5);
}

// put together an array of multifabs for writing
Vector<const MultiFab*> AmrCoreAdv::PlotFileMF() const {
  amrex::Print() << "finest_level = " << finest_level << "\n";
  Vector<const MultiFab*> r(finest_level + 1);
  Vector<MultiFab> plotmf(finest_level + 1);
  const int ncomp_output = 5;
  for (int lev = 0; lev <= finest_level; ++lev) {
    plotmf[lev].define(grids[lev], dmap[lev], ncomp_output, 0);
    int comp = 0;
    MultiFab::Copy(plotmf[lev], moments_new[lev], 0, comp, 4, 0);
    comp += 4;
    MultiFab::Copy(plotmf[lev], band_id[lev], 0, comp, 1, 0);
    r[lev] = &plotmf[lev];
  }
  return r;
}

// set plotfile variable names
Vector<std::string> AmrCoreAdv::PlotFileVarNames() const {
  return {"m0", "m1x", "m1y", "m1z", "band_id"};
}

// write plotfile to disk
void AmrCoreAdv::WritePlotFile() {
  const std::string& plotfilename = PlotFileName(istep[0]);
  // const auto& mf = PlotFileMF();

  Vector<const MultiFab*> mf(finest_level + 1);
  Vector<MultiFab> plotmf(finest_level + 1);
  const int ncomp_output = 5;
  for (int lev = 0; lev <= finest_level; ++lev) {
    plotmf[lev].define(grids[lev], dmap[lev], ncomp_output, 0);
    int comp = 0;
    MultiFab::Copy(plotmf[lev], moments_new[lev], 0, comp, 4, 0);
    comp += 4;
    MultiFab::Copy(plotmf[lev], band_id[lev], 0, comp, 1, 0);
    mf[lev] = &plotmf[lev];
  }

  // const auto& varnames = PlotFileVarNames();
  const auto& varnames =
      Vector<std::string>({"m0", "m1x", "m1y", "m1z", "band_id"});

  amrex::Print() << "Writing plotfile " << plotfilename << "\n";

  amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, mf, varnames,
                                 Geom(), t_new[0], istep, refRatio());

  // Printing interface at finest level
  {
    const int rank = amrex::ParallelDescriptor::MyProc();
    const int size = amrex::ParallelDescriptor::NProcs();
    // Print PVD file
    if (rank == 0) {
      // Write file header
      const std::string pvdfile = std::string("interface.pvd");
      if (istep[0] == 0) {
        std::ofstream outFile(pvdfile);
        outFile << "<?xml version=\"1.0\"?>\n";
        outFile << "<VTKFile type=\"Collection\" version=\"0.1\" "
                   "byte_order=\"LittleEndian\">\n";
        outFile << "  <Collection nsteps=\"1\">\n";
        outFile << "    <DataSet timestep=\"" << std::scientific
                << std::setprecision(6) << t_new[0] << "\" file=\""
                << std::string(PlotFileName(0) + "/Interface/interface.vtu")
                << "\"/>\n";
        outFile << "  </Collection>\n</VTKFile>";
        outFile.close();
      } else {
        std::ifstream inFile(pvdfile);
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(inFile, line)) {
          lines.push_back(line);
        }
        inFile.close();
        std::ostringstream oss;
        oss << "  <Collection nsteps=\"" << istep[0] + 1 << "\">";
        const std::string linensteps = oss.str();
        oss.str("");
        oss.clear();
        oss << "    <DataSet timestep=\"" << std::setprecision(17) << t_new[0]
            << "\" file=\""
            << std::string(PlotFileName(istep[0]) + "/Interface/interface.vtu")
            << "\"/>";
        const std::string linedata = oss.str();
        lines[2] = linensteps;
        const int insertPos = lines.size() - 2;
        lines.insert(lines.begin() + insertPos, linedata);
        std::ofstream outFile(pvdfile);
        for (int i = 0; i < lines.size(); i++) {
          outFile << lines[i] << "\n";
        }
        outFile.close();
      }
    }

    // Triangulate interface
    IRL::MixedPolygonBezierSurface surface;

    MultiFab& moments_lev = moments_new[finest_level];
    SepUnionMultiFab& interface_lev = interface[finest_level];
    const auto problo = Geom(finest_level).ProbLoArray();
    const auto dx = Geom(finest_level).CellSizeArray();

    for (MFIter mfi(moments_lev, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<Real> moments_fab = moments_lev[mfi].array();
      Array4<IRL::SeparatorUnion> interface_fab = interface_lev[mfi].array();
      const Box& box = mfi.tilebox();
      const auto lo = lbound(box), hi = ubound(box);
      for (int k = lo.z; k <= hi.z; ++k) {
        const double z = problo[2] + k * dx[2];
        for (int j = lo.y; j <= hi.y; ++j) {
          const double y = problo[1] + j * dx[1];
          for (int i = lo.x; i <= hi.x; ++i) {
            const double x = problo[0] + i * dx[0];
            const IRL::Pt lower_cell_pt(x, y, z);
            const IRL::Pt upper_cell_pt(x + dx[0], y + dx[1], z + dx[2]);
            const auto cell = IRL::RectangularCuboid::fromBoundingPts(
                lower_cell_pt, upper_cell_pt);

            if (interface_fab(i, j, k).type() ==
                IRL::SeparatorUnion::SeparatorType::OnePlane) {
              const auto planar_sep = IRL::PlanarSeparator::fromOnePlane(
                  interface_fab(i, j, k).getPlane());
              const auto polygon =
                  IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
                      cell, planar_sep, planar_sep[0]);
              surface.addPolygon(polygon);
            } else if (interface_fab(i, j, k).type() ==
                       IRL::SeparatorUnion::SeparatorType::Paraboloid) {
              using VolumeAndSuface = IRL::AddSurfaceOutput<
                  IRL::Volume, IRL::ParaboloidParametrizedSurfaceOutput>;
              const auto paraboloid = interface_fab(i, j, k).getParaboloid();
              auto volume_and_surface =
                  IRL::getVolumeMoments<VolumeAndSuface>(cell, paraboloid);
              surface.addSurface(volume_and_surface.getSurface()
                                     .getQuadraticBezierTriangleApprox());
            }
          }
        }
      }
    }

    const auto points = surface.getPointList();
    const auto polygons = surface.getPolygonList();
    const auto bezier_triangles = surface.getBezierTriangleList();
    const int number_of_points = points.size();
    const int number_of_polygons = polygons.first.size();
    const int number_of_triangles = bezier_triangles.size();

    // Now write
    std::string interfacefoldername = plotfilename;
    interfacefoldername.append("/").append("Interface");

    if (rank == 0) {
      // Create folder to store interface files
      amrex::UtilCreateCleanDirectory(interfacefoldername,
                                      false);  // ---- dont call barrier
      // Write file header
      FILE* file;
      file = fopen(std::string(interfacefoldername + "/interface.vtu").c_str(),
                   "w");
      fprintf(file, "<?xml version=\"1.0\"?>\n");
      fprintf(file,
              "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
              "byte_order=\"LittleEndian\">\n");
      fprintf(file, "  <UnstructuredGrid>\n");
      fclose(file);
    }

    for (int r = 0; r < size; r++) {
      if (rank == r) {
        FILE* file;
        // Write ASCII header for local piece
        file = fopen(
            std::string(interfacefoldername + "/interface.vtu").c_str(), "a");
        fprintf(file, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
                number_of_points, number_of_triangles + number_of_polygons);
        fprintf(file,
                "<Points>\n<DataArray type=\"Float64\" "
                "NumberOfComponents=\"3\">\n");
        for (IRL::UnsignedIndex_t i = 0; i < number_of_points; ++i) {
          fprintf(file, "%15.8E %15.8E %15.8E ", std::get<0>(points[i])[0],
                  std::get<0>(points[i])[1], std::get<0>(points[i])[2]);
        }
        fprintf(file, "\n</DataArray>\n</Points>\n");
        fprintf(file,
                "<PointData RationalWeights=\"RationalWeights\">\n<DataArray "
                "type=\"Float64\" Name=\"RationalWeights\" "
                "format=\"ascii\">\n");
        for (IRL::UnsignedIndex_t i = 0; i < number_of_points; ++i) {
          fprintf(file, "%15.8E ", std::get<1>(points[i]));
        }
        fprintf(file, "\n</DataArray>\n</PointData>\n");

        fprintf(file, "<Cells>\n");
        fprintf(file,
                "<DataArray type=\"Int64\" Name=\"connectivity\" "
                "format=\"ascii\">\n");
        for (IRL::UnsignedIndex_t i = 0; i < number_of_triangles; ++i) {
          for (IRL::UnsignedIndex_t j = 0; j < bezier_triangles[i].size();
               j++) {
            fprintf(file, "%d ", bezier_triangles[i][j]);
          }
        }
        for (IRL::UnsignedIndex_t i = 0; i < polygons.second.size(); ++i) {
          fprintf(file, "%d ", polygons.second[i]);
        }
        fprintf(file, "\n</DataArray>\n");

        fprintf(
            file,
            "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
        IRL::UnsignedIndex_t count = 0;
        for (IRL::UnsignedIndex_t i = 0; i < number_of_triangles; ++i) {
          count += bezier_triangles[i].size();
          fprintf(file, "%d ", count);
        }
        for (IRL::UnsignedIndex_t i = 0; i < number_of_polygons; ++i) {
          count += polygons.first[i];
          fprintf(file, "%d ", count);
        }
        fprintf(file, "\n</DataArray>\n");
        fprintf(file,
                "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
        for (IRL::UnsignedIndex_t i = 0; i < number_of_triangles; ++i) {
          fprintf(file, "76 ");
        }
        for (IRL::UnsignedIndex_t i = 0; i < number_of_polygons; ++i) {
          fprintf(file, "7 ");
        }
        fprintf(file, "\n</DataArray>\n");
        fprintf(file, "</Cells>\n");
        fprintf(file, "</Piece>\n");
        fclose(file);
      }
      ParallelDescriptor::Barrier();
    }

    if (rank == 0) {
      // Write file footer
      FILE* file;
      file = fopen(std::string(interfacefoldername + "/interface.vtu").c_str(),
                   "a");
      fprintf(file, "</UnstructuredGrid>\n</VTKFile>");
      fclose(file);
    }
  }
}

void AmrCoreAdv::WriteCheckpointFile() const {
  // chk00010            write a checkpoint file with this root directory
  // chk00010/Header     this contains information you need to save (e.g.,
  // finest_level, t_new, etc.) and also
  //                     the BoxArrays at each level
  // chk00010/Level_0/
  // chk00010/Level_1/
  // etc.                these subdirectories will hold the MultiFab data at
  // each level of refinement

  // checkpoint file name, e.g., chk00010
  const std::string& checkpointname = amrex::Concatenate(chk_file, istep[0]);

  amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

  const int nlevels = finest_level + 1;

  // ---- prebuild a hierarchy of directories
  // ---- dirName is built first.  if dirName exists, it is renamed.  then
  // build
  // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
  // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
  // ---- after all directories are built
  // ---- ParallelDescriptor::IOProcessor() creates the directories
  amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

  // write Header file
  if (ParallelDescriptor::IOProcessor()) {
    std::string HeaderFileName(checkpointname + "/Header");
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    std::ofstream HeaderFile;
    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
    if (!HeaderFile.good()) {
      amrex::FileOpenFailed(HeaderFileName);
    }

    HeaderFile.precision(17);

    // write out title line
    HeaderFile << "Checkpoint file for AmrCoreAdv\n";

    // write out finest_level
    HeaderFile << finest_level << "\n";

    // write out array of istep
    for (int i = 0; i < istep.size(); ++i) {
      HeaderFile << istep[i] << " ";
    }
    HeaderFile << "\n";

    // write out array of dt
    for (int i = 0; i < dt.size(); ++i) {
      HeaderFile << dt[i] << " ";
    }
    HeaderFile << "\n";

    // write out array of t_new
    for (int i = 0; i < t_new.size(); ++i) {
      HeaderFile << t_new[i] << " ";
    }
    HeaderFile << "\n";

    // write the BoxArray at each level
    for (int lev = 0; lev <= finest_level; ++lev) {
      boxArray(lev).writeOn(HeaderFile);
      HeaderFile << '\n';
    }
  }

  // write the MultiFab data to, e.g., chk00010/Level_0/
  for (int lev = 0; lev <= finest_level; ++lev) {
    VisMF::Write(moments_new[lev],
                 amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_",
                                               "moments"));
  }
}

namespace {
// utility to skip to next line in Header
void GotoNextLine(std::istream& is) {
  constexpr std::streamsize bl_ignore_max{100000};
  is.ignore(bl_ignore_max, '\n');
}
}  // namespace

void AmrCoreAdv::ReadCheckpointFile() {
  amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

  // Header
  std::string File(restart_chkfile + "/Header");

  VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

  Vector<char> fileCharPtr;
  ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
  std::string fileCharPtrString(fileCharPtr.dataPtr());
  std::istringstream is(fileCharPtrString, std::istringstream::in);

  std::string line, word;

  // read in title line
  std::getline(is, line);

  // read in finest_level
  is >> finest_level;
  GotoNextLine(is);

  // read in array of istep
  std::getline(is, line);
  {
    std::istringstream lis(line);
    int i = 0;
    while (lis >> word) {
      istep[i++] = std::stoi(word);
    }
  }

  // read in array of dt
  std::getline(is, line);
  {
    std::istringstream lis(line);
    int i = 0;
    while (lis >> word) {
      dt[i++] = std::stod(word);
    }
  }

  // read in array of t_new
  std::getline(is, line);
  {
    std::istringstream lis(line);
    int i = 0;
    while (lis >> word) {
      t_new[i++] = std::stod(word);
    }
  }

  for (int lev = 0; lev <= finest_level; ++lev) {
    // read in level 'lev' BoxArray from Header
    BoxArray ba;
    ba.readFrom(is);
    GotoNextLine(is);

    // create a distribution mapping
    DistributionMapping dm{ba, ParallelDescriptor::NProcs()};

    // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H
    // class
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    // build MultiFab and FluxRegister data
    int ncomp = 4;
    int nghost = 0;
    moments_old[lev].define(grids[lev], dmap[lev], ncomp, nghost);
    moments_new[lev].define(grids[lev], dmap[lev], ncomp, nghost);

    if (lev > 0 && do_reflux) {
      flux_reg[lev].reset(new FluxRegister(grids[lev], dmap[lev],
                                           refRatio(lev - 1), lev, ncomp));
    }

    // build face velocity MultiFabs
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      facevel[lev][idim] =
          MultiFab(amrex::convert(ba, IntVect::TheDimensionVector(idim)), dm, 1,
                   num_grow);
    }
  }

  // read in the MultiFab data
  for (int lev = 0; lev <= finest_level; ++lev) {
    VisMF::Read(moments_new[lev],
                amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_",
                                              "moments"));
  }
}

// inline IRL::Vec3<double> AmrCoreAdv::getVelocity(const IRL::Pt& pt,
//                                                  Array4<Real const> const&
//                                                  vx, Array4<Real const>
//                                                  const& vy, Array4<Real
//                                                  const> const& vz, const Box&
//                                                  bx, const int lev) {
//   const auto& dx = geom[lev].CellSizeArray();
//   const auto& prob_lo = geom[lev].ProbLoArray();
//   const auto& lo = lbound(bx);
//   const auto& hi = ubound(bx);
//   // Find which cell the point falls in
//   int i = static_cast<int>(amrex::Math::floor((pt[0] - prob_lo[0]) / dx[0]));
//   int j = static_cast<int>(amrex::Math::floor((pt[1] - prob_lo[1]) / dx[1]));
//   int k = static_cast<int>(amrex::Math::floor((pt[2] - prob_lo[2]) / dx[2]));
//   if (!bx.contains(i, j, k)) {
//     std::ostringstream oss;
//     oss << "Position (" << pt[0] << "," << pt[1] << "," << pt[2]
//         << ") leads to index (" << i << "," << j << "," << k
//         << ") outside of box (" << lo.x << "," << lo.y << "," << lo.z <<
//         ")x("
//         << hi.x << "," << hi.y << "," << hi.z << ")";
//     throw std::runtime_error(oss.str());
//   }
//   // Cell lo face positions
//   Real xlo = prob_lo[0] + i * dx[0];
//   Real ylo = prob_lo[1] + j * dx[1];
//   Real zlo = prob_lo[2] + k * dx[2];
//   // Normalized position within cell [0,1]
//   Real tx = (pt[0] - xlo) / dx[0];
//   Real ty = (pt[1] - ylo) / dx[1];
//   Real tz = (pt[2] - zlo) / dx[2];
//   return IRL::Vec3<double>(vx(i, j, k) * (1.0 - tx) + vx(i + 1, j, k) * tx,
//                            vy(i, j, k) * (1.0 - ty) + vy(i, j + 1, k) * ty,
//                            vz(i, j, k) * (1.0 - tz) + vz(i, j, k + 1) * tz);
// }

// inline IRL::Pt AmrCoreAdv::project_vertex(const IRL::Pt& pt, const double dt,
//                                           Array4<Real const> const& vx,
//                                           Array4<Real const> const& vy,
//                                           Array4<Real const> const& vz,
//                                           const Box& bx, const int lev) {
//   using Pt = IRL::Pt;
//   auto v1 = getVelocity(pt, vx, vy, vz, bx, lev);
//   auto v2 = getVelocity(pt + Pt::fromVec3(0.5 * dt * v1), vx, vy, vz, bx,
//   lev); auto v3 = getVelocity(pt + Pt::fromVec3(0.5 * dt * v2), vx, vy, vz,
//   bx, lev); auto v4 = getVelocity(pt + Pt::fromVec3(dt * v3), vx, vy, vz, bx,
//   lev); return pt + Pt::fromVec3(dt * (v1 + 2.0 * v2 + 2.0 * v3 + v4) / 6.0);
// }

void AmrCoreAdv::UpdateBand() {
  for (int lev = 0; lev <= finest_level; lev++) {
    // Get geometrical information
    const auto dx = geom[lev].CellSizeArray();
    const auto problo = geom[lev].ProbLoArray();
    const Real vol = dx[0] * dx[1] * dx[2];
    const Real vol_inv = 1.0 / vol;

    moments_new[lev].FillBoundary(geom[lev].periodicity());
    band_id[lev].setVal(0.0);

    for (MFIter mfi(moments_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<Real const> moments = moments_new[lev].const_array(mfi);
      Array4<Real> band = band_id[lev].array(mfi);
      const Box& bx = mfi.tilebox();
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (moments(i, j, k, 0) * vol_inv > IRL::global_constants::VF_LOW &&
            moments(i, j, k, 0) * vol_inv < IRL::global_constants::VF_HIGH) {
          band(i, j, k) = 1.0;
        }
      });
    }
    band_id[lev].FillBoundary(geom[lev].periodicity());

    const int nlayers = static_cast<int>(static_cast<double>(num_grow) /
                                         std::pow(2.0, finest_level - lev));
    for (int n = 0; n < nlayers; ++n) {
      for (MFIter mfi(band_id[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<Real> band = band_id[lev].array(mfi);
        const Box& bx = mfi.tilebox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          if (std::fabs(band(i, j, k)) < 0.5) {
            for (int ii = -1; ii <= 1; ++ii) {
              for (int jj = -1; jj <= 1; ++jj) {
                for (int kk = -1; kk <= 1; ++kk) {
                  if (std::fabs(band(i + ii, j + jj, k + kk) - (n + 1)) < 0.5) {
                    band(i, j, k) = n + 2;
                  }
                }
              }
            }
          }
        });
      }
      band_id[lev].FillBoundary(geom[lev].periodicity());
    }
  }
}

void AmrCoreAdv::AdvanceAllLevels(Real time, Real dt_lev, int /*iteration*/) {
  // Update band_id identifying region close to interface
  UpdateBand();

  // for (int lev = 0; lev <= finest_level; lev++) {
  const int lev = finest_level;
  {
    // Get geometrical information
    const auto dx = geom[lev].CellSizeArray();
    const auto problo = geom[lev].ProbLoArray();
    const Real vol = dx[0] * dx[1] * dx[2];
    const Real vol_inv = 1.0 / vol;

    // Store old moments
    MultiFab::Copy(moments_old[lev], moments_new[lev], 0, 0,
                   moments_new[lev].nComp(), moments_new[lev].nGrow());
    // std::swap(moments_old[lev], moments_new[lev]);
    t_old[lev] = t_new[lev];
    t_new[lev] += dt_lev;

    // Build tmp moment multifab with ghost layers
    MultiFab moments_with_ghost(grids[lev], dmap[lev], ncomp_moments, num_grow);
    MultiFab::Copy(moments_with_ghost, moments_old[lev], 0, 0,
                   moments_old[lev].nComp(), moments_old[lev].nGrow());
    moments_with_ghost.FillBoundary(geom[lev].periodicity());

    // Build tmp interface multifab with ghost layers
    SepUnionMultiFab interface_with_ghost(grids[lev], dmap[lev],
                                          interface[lev].nComp(), num_grow);

    // Reconstruct interface and update ghosts
    amrex::ParallelDescriptor::Barrier();
    const auto rec_start = amrex::second();
    GetReconstruction(interface[lev], interface_with_ghost, moments_with_ghost,
                      geom[lev]);
    amrex::ParallelDescriptor::Barrier();
    reconstruction_time += amrex::second() - rec_start;

    // Advect moments using reconstructed interface
    amrex::ParallelDescriptor::Barrier();
    const auto adv_start = amrex::second();
    TransportMoments(interface_with_ghost, facevel[lev], band_id[lev],
                     moments_new[lev], geom[lev], dt[lev]);
    amrex::ParallelDescriptor::Barrier();
    advection_time += amrex::second() - adv_start;
  }  // end lev
}

void AmrCoreAdv::DefineVelocityAllLevels(Real time) {
  for (int lev = 0; lev <= finest_level; ++lev)
    DefineVelocityAtLevel(lev, time);
}

void AmrCoreAdv::DefineVelocityAtLevel(int lev, Real time) {
  for (MFIter mfi(moments_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    // ======== GET FACE VELOCITY =========
    GpuArray<Box, AMREX_SPACEDIM> nbx;
    AMREX_D_TERM(nbx[0] = mfi.nodaltilebox(0);, nbx[1] = mfi.nodaltilebox(1);
                 , nbx[2] = mfi.nodaltilebox(2););

    AMREX_D_TERM(
        const Box& ngbxx = amrex::grow(mfi.nodaltilebox(0), num_grow);
        , const Box& ngbxy = amrex::grow(mfi.nodaltilebox(1), num_grow);
        , const Box& ngbxz = amrex::grow(mfi.nodaltilebox(2), num_grow););

    GpuArray<Array4<Real>, AMREX_SPACEDIM> vel{
        AMREX_D_DECL(facevel[lev][0].array(mfi), facevel[lev][1].array(mfi),
                     facevel[lev][2].array(mfi))};

    GeometryData geomdata = geom[lev].data();
    const Real* AMREX_RESTRICT dx = geomdata.CellSize();
    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();

    if (case_name == "deformation3d") {
      amrex::ParallelFor(
          AMREX_D_DECL(ngbxx, ngbxy, ngbxz),
          AMREX_D_DECL(
              // X-faces
              [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Real x = prob_lo[0] + i * dx[0];
                Real y = prob_lo[1] + (j + 0.5) * dx[1];
                Real z = prob_lo[2] + (k + 0.5) * dx[2];
                vel[0](i, j, k) =
                    Deformation3D::get_face_velocity_x(x, y, z, time);
              },
              // Y-faces
              [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Real x = prob_lo[0] + (i + 0.5) * dx[0];
                Real y = prob_lo[1] + j * dx[1];
                Real z = prob_lo[2] + (k + 0.5) * dx[2];
                vel[1](i, j, k) =
                    Deformation3D::get_face_velocity_y(x, y, z, time);
              },
              // Z-faces
              [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Real x = prob_lo[0] + (i + 0.5) * dx[0];
                Real y = prob_lo[1] + (j + 0.5) * dx[1];
                Real z = prob_lo[2] + k * dx[2];
                vel[2](i, j, k) =
                    Deformation3D::get_face_velocity_z(x, y, z, time);
              }));
    } else if (case_name == "translation3d" || case_name == "default") {
      amrex::ParallelFor(
          AMREX_D_DECL(ngbxx, ngbxy, ngbxz),
          AMREX_D_DECL(
              // X-faces
              [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                vel[0](i, j, k) = Translation3D::get_face_velocity_x();
              },
              // Y-faces
              [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                vel[1](i, j, k) = Translation3D::get_face_velocity_y();
              },
              // Z-faces
              [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                vel[2](i, j, k) = Translation3D::get_face_velocity_z();
              }));
    } else if (case_name == "rotation3d") {
      amrex::ParallelFor(
          AMREX_D_DECL(ngbxx, ngbxy, ngbxz),
          AMREX_D_DECL(
              // X-faces
              [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Real x = prob_lo[0] + i * dx[0];
                Real y = prob_lo[1] + (j + 0.5) * dx[1];
                Real z = prob_lo[2] + (k + 0.5) * dx[2];
                vel[0](i, j, k) =
                    Rotation3D::get_face_velocity_x(x, y, z, time);
              },
              // Y-faces
              [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Real x = prob_lo[0] + (i + 0.5) * dx[0];
                Real y = prob_lo[1] + j * dx[1];
                Real z = prob_lo[2] + (k + 0.5) * dx[2];
                vel[1](i, j, k) =
                    Rotation3D::get_face_velocity_y(x, y, z, time);
              },
              // Z-faces
              [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Real x = prob_lo[0] + (i + 0.5) * dx[0];
                Real y = prob_lo[1] + (j + 0.5) * dx[1];
                Real z = prob_lo[2] + k * dx[2];
                vel[2](i, j, k) =
                    Rotation3D::get_face_velocity_z(x, y, z, time);
              }));
    } else {
      throw std::runtime_error("Unknown case");
    }
  }
}

Real AmrCoreAdv::RecTime() { return reconstruction_time; }

Real AmrCoreAdv::AdvTime() { return advection_time; }

// void AmrCoreAdv::BuildUniformFinestMoments(
//     amrex::MultiFab& a_uniform_moments) const {
//   using namespace amrex;

//   const int finest = finest_level;
//   const int ncomp = ncomp_moments;
//   const int ngrow = 0;

//   //
//   ---------------------------------------------------------------------------
//   // 1. Build full-domain uniform grid at finest AMR resolution.
//   //
//   ---------------------------------------------------------------------------
//   BoxArray uniform_ba(Geom(finest).Domain());

//   // Optional: controls box decomposition only, not grid spacing.
//   uniform_ba.maxSize(32);

//   DistributionMapping uniform_dm(uniform_ba);

//   a_uniform_moments.define(uniform_ba, uniform_dm, ncomp, ngrow);
//   a_uniform_moments.setVal(0.0);

//   // Finest-level geometry.
//   const auto fine_dx = Geom(finest).CellSizeArray();
//   const auto problo = Geom(finest).ProbLoArray();

//   const Real fine_vol = fine_dx[0] * fine_dx[1] * fine_dx[2];

//   const Box& fine_domain = Geom(finest).Domain();
//   const IntVect fine_dom_lo = fine_domain.smallEnd();

//   const int fine_lo_x = fine_dom_lo[0];
//   const int fine_lo_y = fine_dom_lo[1];
//   const int fine_lo_z = fine_dom_lo[2];

//   //
//   ---------------------------------------------------------------------------
//   // 2. Fill from coarse to fine.
//   //
//   // Each AMR level is expanded into finest-level index space.
//   // Then ParallelCopy transfers that data into the full uniform diagnostic
//   // grid.
//   //
//   // Coarse levels fill first. Finer levels overwrite refined regions later.
//   //
//   ---------------------------------------------------------------------------
//   for (int lev = 0; lev <= finest; ++lev) {
//     // Ratio from level lev to finest level.
//     IntVect ratio_vect(AMREX_D_DECL(1, 1, 1));

//     for (int l = lev; l < finest; ++l) {
//       ratio_vect *= refRatio(l);
//     }

//     const int rx = ratio_vect[0];
//     const int ry = ratio_vect[1];
//     const int rz = ratio_vect[2];

//     const auto lev_dx = Geom(lev).CellSizeArray();
//     const Real lev_vol = lev_dx[0] * lev_dx[1] * lev_dx[2];

//     const Box& lev_domain = Geom(lev).Domain();
//     const IntVect lev_dom_lo = lev_domain.smallEnd();

//     const int lev_lo_x = lev_dom_lo[0];
//     const int lev_lo_y = lev_dom_lo[1];
//     const int lev_lo_z = lev_dom_lo[2];

//     //
//     -------------------------------------------------------------------------
//     // Build a temporary MultiFab whose boxes are the AMR boxes on this
//     level,
//     // but refined into finest-level index space.
//     //
//     // Use the same DistributionMapping as moments_new[lev], so each rank
//     // expands the source boxes it already owns.
//     //
//     -------------------------------------------------------------------------
//     BoxArray lev_on_fine_ba = moments_new[lev].boxArray();
//     lev_on_fine_ba.refine(ratio_vect);

//     MultiFab lev_on_fine(lev_on_fine_ba, moments_new[lev].DistributionMap(),
//                          ncomp, 0);

//     lev_on_fine.setVal(0.0);

//     //
//     -------------------------------------------------------------------------
//     // Locally fill lev_on_fine from moments_new[lev].
//     //
//     -------------------------------------------------------------------------
//     for (MFIter mfi(lev_on_fine, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
//       const Box& fine_box = mfi.tilebox();

//       auto fine_arr = lev_on_fine.array(mfi);

//       // Important:
//       // mfi iterates over lev_on_fine, whose BoxArray is refined compared
//       with
//       // moments_new[lev]. Use mfi.index() to access the corresponding source
//       // FAB.
//       auto lev_arr = moments_new[lev].const_array(mfi.index());

//       amrex::ParallelFor(
//           fine_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
//             // Physical centroid of the uniform finest cell.
//             const Real x =
//                 problo[0] + (static_cast<Real>(i) + Real(0.5)) * fine_dx[0];
//             const Real y =
//                 problo[1] + (static_cast<Real>(j) + Real(0.5)) * fine_dx[1];
//             const Real z =
//                 problo[2] + (static_cast<Real>(k) + Real(0.5)) * fine_dx[2];

//             if (lev == finest) {
//               //
//               ----------------------------------------------------------------
//               // Finest AMR level:
//               // Copy stored moments directly. These cells may be mixed,
//               liquid,
//               // or gas, but the stored moments are already at finest
//               // resolution.
//               //
//               ----------------------------------------------------------------
//               fine_arr(i, j, k, 0) = lev_arr(i, j, k, 0);
//               fine_arr(i, j, k, 1) = lev_arr(i, j, k, 1);
//               fine_arr(i, j, k, 2) = lev_arr(i, j, k, 2);
//               fine_arr(i, j, k, 3) = lev_arr(i, j, k, 3);

//             } else {
//               //
//               ----------------------------------------------------------------
//               // Coarser AMR level:
//               // Map finest-level index back to this AMR level.
//               //
//               ----------------------------------------------------------------
//               const int ic = lev_lo_x + (i - fine_lo_x) / rx;
//               const int jc = lev_lo_y + (j - fine_lo_y) / ry;
//               const int kc = lev_lo_z + (k - fine_lo_z) / rz;

//               const Real vf = lev_arr(ic, jc, kc, 0) / lev_vol;

//               if (vf > IRL::global_constants::VF_HIGH) {
//                 // Pure liquid coarse cell.
//                 // Every generated finest cell inside it is full liquid.
//                 fine_arr(i, j, k, 0) = fine_vol;
//                 fine_arr(i, j, k, 1) = fine_vol * x;
//                 fine_arr(i, j, k, 2) = fine_vol * y;
//                 fine_arr(i, j, k, 3) = fine_vol * z;

//               } else if (vf < IRL::global_constants::VF_LOW) {
//                 // Pure gas coarse cell.
//                 // Every generated finest cell inside it is empty.
//                 fine_arr(i, j, k, 0) = Real(0.0);
//                 fine_arr(i, j, k, 1) = Real(0.0);
//                 fine_arr(i, j, k, 2) = Real(0.0);
//                 fine_arr(i, j, k, 3) = Real(0.0);
//               }
//               // else {
//               //   // Unexpected case:
//               //   // A mixed cell exists below the finest level.
//               //   //
//               //   // This should not happen if AMR refinement is doing what
//               you
//               //   // want. But do not leave it unwritten. Use a simple
//               //   uniform-vf
//               //   // fallback.
//               //   fine_arr(i, j, k, 0) = vf * fine_vol;
//               //   fine_arr(i, j, k, 1) = vf * fine_vol * x;
//               //   fine_arr(i, j, k, 2) = vf * fine_vol * y;
//               //   fine_arr(i, j, k, 3) = vf * fine_vol * z;
//               // }
//             }
//           });
//     }

//     //
//     -------------------------------------------------------------------------
//     // MPI-safe copy into the full uniform diagnostic grid.
//     //
//     // lev_on_fine and a_uniform_moments are both in finest-level index
//     space.
//     // ParallelCopy handles communication between MPI ranks.
//     //
//     -------------------------------------------------------------------------
//     a_uniform_moments.ParallelCopy(lev_on_fine, 0, 0, ncomp, 0, 0,
//                                    Geom(finest).periodicity());
//   }
// }

// void AmrCoreAdv::WriteUniformMomentsBinary(
//     const amrex::MultiFab& a_uniform_moments,
//     const std::string& a_filename) const {
//   using namespace amrex;

//   const int finest = finest_level;
//   const int ncomp = a_uniform_moments.nComp();

//   const Box& domain = Geom(finest).Domain();
//   const IntVect lo = domain.smallEnd();
//   const IntVect hi = domain.bigEnd();

//   const int nx = domain.length(0);
//   const int ny = domain.length(1);
//   const int nz = domain.length(2);

//   // Gather the distributed uniform MultiFab onto the IO processor
//   // as one full-domain box, so the binary file has a simple ordering.
//   BoxArray single_ba(domain);

//   Vector<int> pmap(1, ParallelDescriptor::IOProcessorNumber());
//   DistributionMapping single_dm(pmap);

//   MultiFab uniform_single(single_ba, single_dm, ncomp, 0);
//   uniform_single.setVal(0.0);

//   uniform_single.ParallelCopy(a_uniform_moments, 0, 0, ncomp, 0, 0,
//                               Geom(finest).periodicity());

//   if (ParallelDescriptor::IOProcessor()) {
//     std::ofstream ofs(a_filename, std::ios::out | std::ios::binary);

//     if (!ofs.good()) {
//       amrex::Abort("Could not open binary output file: " + a_filename);
//     }

//     //
//     -----------------------------------------------------------------------
//     // Header
//     //
//     -----------------------------------------------------------------------

//     // Magic number and version.
//     // These help MATLAB check that it is reading the right file type.
//     const std::int32_t magic = 20260504;
//     const std::int32_t version = 1;

//     ofs.write(reinterpret_cast<const char*>(&magic), sizeof(std::int32_t));
//     ofs.write(reinterpret_cast<const char*>(&version), sizeof(std::int32_t));

//     // Write case.name
//     {
//       const std::int32_t len = static_cast<std::int32_t>(case_name.size());

//       ofs.write(reinterpret_cast<const char*>(&len), sizeof(std::int32_t));
//       ofs.write(case_name.data(), len);
//     }

//     // Write reconstruction.name
//     {
//       const std::int32_t len =
//           static_cast<std::int32_t>(reconstruction_name.size());

//       ofs.write(reinterpret_cast<const char*>(&len), sizeof(std::int32_t));
//       ofs.write(reconstruction_name.data(), len);
//     }

//     // Integer mesh/header data.
//     // 13 int32 values:
//     // lo_x lo_y lo_z
//     // hi_x hi_y hi_z
//     // nx ny nz
//     // ncomp
//     // finest_level
//     // max_level
//     // reserved
//     std::int32_t header[13];

//     header[0] = static_cast<std::int32_t>(lo[0]);
//     header[1] = static_cast<std::int32_t>(lo[1]);
//     header[2] = static_cast<std::int32_t>(lo[2]);

//     header[3] = static_cast<std::int32_t>(hi[0]);
//     header[4] = static_cast<std::int32_t>(hi[1]);
//     header[5] = static_cast<std::int32_t>(hi[2]);

//     header[6] = static_cast<std::int32_t>(nx);
//     header[7] = static_cast<std::int32_t>(ny);
//     header[8] = static_cast<std::int32_t>(nz);

//     header[9] = static_cast<std::int32_t>(ncomp);

//     header[10] = static_cast<std::int32_t>(finest_level);
//     header[11] = static_cast<std::int32_t>(max_level);
//     header[12] = static_cast<std::int32_t>(0);  // reserved

//     ofs.write(reinterpret_cast<const char*>(header), sizeof(header));

//     // Optional printout
//     amrex::Print() << "Writing uniform moments binary file: " << a_filename
//                    << "\n";
//     amrex::Print() << "  case.name = " << case_name << "\n";
//     amrex::Print() << "  reconstruction.name = " << reconstruction_name <<
//     "\n"; amrex::Print() << "  uniform finest mesh = " << nx << " x " << ny
//     << " x "
//                    << nz << "\n";
//     amrex::Print() << "  ncomp = " << ncomp << "\n";

//     //
//     -----------------------------------------------------------------------
//     // Data layout:
//     //
//     // for k = lo_z : hi_z
//     //   for j = lo_y : hi_y
//     //     for i = lo_x : hi_x
//     //       for n = 0 : ncomp-1
//     //         write moment(i,j,k,n)
//     //
//     // Components:
//     //   n = 0 : m0
//     //   n = 1 : m1x
//     //   n = 2 : m1y
//     //   n = 3 : m1z
//     //
//     -----------------------------------------------------------------------

//     for (MFIter mfi(uniform_single); mfi.isValid(); ++mfi) {
//       auto const arr = uniform_single.const_array(mfi);

//       for (int k = lo[2]; k <= hi[2]; ++k) {
//         for (int j = lo[1]; j <= hi[1]; ++j) {
//           for (int i = lo[0]; i <= hi[0]; ++i) {
//             for (int n = 0; n < ncomp; ++n) {
//               const double value = static_cast<double>(arr(i, j, k, n));
//               ofs.write(reinterpret_cast<const char*>(&value),
//               sizeof(double));
//             }
//           }
//         }
//       }
//     }

//     ofs.close();
//   }
// }

// std::string AmrCoreAdv::UniformMomentsBinaryFileName(
//     const std::string& a_label) const {
//   const int finest = finest_level;

//   const amrex::Box& domain = Geom(finest).Domain();

//   const int nx = domain.length(0);
//   const int ny = domain.length(1);
//   const int nz = domain.length(2);

//   std::ostringstream oss;
//   oss << "uniform_moments_" << a_label << "_" << case_name << "_"
//       << reconstruction_name << "_" << nx << ".bin";

//   return oss.str();
// }

// amrex::Real AmrCoreAdv::ComputeCompositeM0() const {
//   Real local_sum = 0.0;

//   for (int lev = 0; lev <= finest_level; ++lev) {
//     // Build a mask: 1 means include this cell, 0 means covered by finer
//     level. iMultiFab mask(moments_new[lev].boxArray(),
//                    moments_new[lev].DistributionMap(), 1, 0);

//     mask.setVal(1);

//     if (lev < finest_level) {
//       // Mark cells covered by level lev+1 as 0.
//       BoxArray fine_ba = moments_new[lev + 1].boxArray();

//       // Convert fine boxes down to this level's index space.
//       fine_ba.coarsen(refRatio(lev));

//       for (MFIter mfi(mask); mfi.isValid(); ++mfi) {
//         const Box& bx = mfi.validbox();
//         auto mask_arr = mask.array(mfi);

//         for (int ibox = 0; ibox < fine_ba.size(); ++ibox) {
//           Box covered = bx & fine_ba[ibox];

//           if (!covered.ok()) continue;

//           amrex::Loop(covered, [=](int i, int j, int k) noexcept {
//             mask_arr(i, j, k) = 0;
//           });
//         }
//       }
//     }

//     // Sum m0 only where mask == 1.
//     for (MFIter mfi(moments_new[lev]); mfi.isValid(); ++mfi) {
//       const Box& bx = mfi.validbox();

//       auto m_arr = moments_new[lev].const_array(mfi);
//       auto mask_arr = mask.const_array(mfi);

//       Real fab_sum = 0.0;

//       amrex::Loop(bx, [=, &fab_sum](int i, int j, int k) noexcept {
//         if (mask_arr(i, j, k) == 1) {
//           fab_sum += m_arr(i, j, k, 0);
//         }
//       });

//       local_sum += fab_sum;
//     }
//   }

//   ParallelDescriptor::ReduceRealSum(local_sum);

//   return local_sum;
// }

// amrex::Real AmrCoreAdv::ComputeCompositeM0() const {
//   amrex::Real local_sum = 0.0;

//   const int lev = 0;

//   for (MFIter mfi(moments_new[lev]); mfi.isValid(); ++mfi) {
//     const Box& bx = mfi.validbox();

//     auto m_arr = moments_new[lev].const_array(mfi);

//     amrex::Real fab_sum = 0.0;

//     amrex::Loop(bx, [=, &fab_sum](int i, int j, int k) noexcept {
//       fab_sum += m_arr(i, j, k, 0);
//     });

//     local_sum += fab_sum;
//   }

//   amrex::ParallelDescriptor::ReduceRealSum(local_sum);

//   return local_sum;
// }

// summing zeroth moment on a mesh with no refinement
amrex::Real AmrCoreAdv::ComputeCompositeM0() const {
  amrex::Real local_sum = 0.0;

  for (MFIter mfi(moments_new[0]); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.validbox();

    auto m_arr = moments_new[0].const_array(mfi);

    amrex::Real fab_sum = 0.0;

    for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
      for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
          fab_sum += m_arr(i, j, k, 0);
        }
      }
    }

    local_sum += fab_sum;
  }

  amrex::ParallelDescriptor::ReduceRealSum(local_sum);

  return local_sum;
}

// computing error norms for a mesh with
// void AmrCoreAdv::ComputeUniformMomentL1Errors(
//     const amrex::MultiFab& a_initial, const amrex::MultiFab& a_final) const {
//   using namespace amrex;

//   if (a_initial.nComp() < 4 || a_final.nComp() < 4) {
//     amrex::Abort(
//         "ComputeUniformMomentL1Errors: expected at least 4 components.");
//   }

//   if (a_initial.boxArray().ixType() != a_final.boxArray().ixType()) {
//     amrex::Abort(
//         "ComputeUniformMomentL1Errors: initial/final index types differ.");
//   }

//   // The initial and final uniform grids must represent the same
//   // resolution/domain.
//   if (a_initial.boxArray().minimalBox() != a_final.boxArray().minimalBox()) {
//     amrex::Print() << "Initial uniform domain = "
//                    << a_initial.boxArray().minimalBox() << "\n";
//     amrex::Print() << "Final uniform domain   = "
//                    << a_final.boxArray().minimalBox() << "\n";
//     amrex::Abort(
//         "ComputeUniformMomentL1Errors: initial/final uniform domains
//         differ.");
//   }

//   // Put final data on the exact same BoxArray/DistributionMapping as
//   initial.
//   // This makes MFIter indexing safe.
//   MultiFab final_on_initial_layout(
//       a_initial.boxArray(), a_initial.DistributionMap(), a_initial.nComp(),
//       0);

//   final_on_initial_layout.setVal(0.0);

//   final_on_initial_layout.ParallelCopy(a_final, 0, 0, 4, 0, 0,
//                                        Geom(finest_level).periodicity());

//   Real local_L1_M0 = 0.0;
//   Real local_L1_M1mag = 0.0;

//   for (MFIter mfi(a_initial, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
//     const Box& bx = mfi.tilebox();

//     auto init = a_initial.const_array(mfi);
//     auto fin = final_on_initial_layout.const_array(mfi);

//     Real fab_L1_M0 = 0.0;
//     Real fab_L1_M1mag = 0.0;

//     amrex::Loop(bx,
//                 [=, &fab_L1_M0, &fab_L1_M1mag](int i, int j, int k) noexcept
//                 {
//                   const Real M0_i = init(i, j, k, 0);
//                   const Real M0_f = fin(i, j, k, 0);

//                   const Real M1x_i = init(i, j, k, 1);
//                   const Real M1y_i = init(i, j, k, 2);
//                   const Real M1z_i = init(i, j, k, 3);

//                   const Real M1x_f = fin(i, j, k, 1);
//                   const Real M1y_f = fin(i, j, k, 2);
//                   const Real M1z_f = fin(i, j, k, 3);

//                   const Real M1mag_i =
//                       std::sqrt(M1x_i * M1x_i + M1y_i * M1y_i + M1z_i *
//                       M1z_i);

//                   const Real M1mag_f =
//                       std::sqrt(M1x_f * M1x_f + M1y_f * M1y_f + M1z_f *
//                       M1z_f);

//                   fab_L1_M0 += std::abs(M0_f - M0_i);
//                   fab_L1_M1mag += std::abs(M1mag_f - M1mag_i);
//                 });

//     local_L1_M0 += fab_L1_M0;
//     local_L1_M1mag += fab_L1_M1mag;
//   }

//   ParallelDescriptor::ReduceRealSum(local_L1_M0);
//   ParallelDescriptor::ReduceRealSum(local_L1_M1mag);

//   amrex::Print() << "\nUniform fine-grid L1 error norms\n";
//   amrex::Print() << "  L1(M0)       = " << std::setprecision(17) <<
//   local_L1_M0
//                  << "\n";
//   amrex::Print() << "  L1(|M1|)     = " << std::setprecision(17)
//                  << local_L1_M1mag << "\n";
// }

// computing error norm for a uniform mesh (no refinements)
amrex::Real AmrCoreAdv::ComputeL1ErrorM0() const {
  Real local_l1 = 0.0;

  for (MFIter mfi(moments_new[0]); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.validbox();

    auto final_arr = moments_new[0].const_array(mfi);
    auto init_arr = initial_moments.const_array(mfi);

    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    Real fab_l1 = 0.0;

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
          const Real mom_err = final_arr(i, j, k, 0) - init_arr(i, j, k, 0);
          fab_l1 += std::abs(mom_err);
        }
      }
    }

    local_l1 += fab_l1;
  }

  ParallelDescriptor::ReduceRealSum(local_l1);

  return local_l1;
}