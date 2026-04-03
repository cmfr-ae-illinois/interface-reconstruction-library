
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <fstream>
#include <sstream>

// #include <Kernels.H>
#include "examples/amrex_advector/deformation_3d.h"

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

  ParmParse pp("case");
  pp.get("name", case_name);
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

  if (plot_int > 0 && istep[0] > last_plot_file_step) {
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
  new_interface.ParallelCopy(interface[lev], 0, 0, 1, nghost, nghost,
                             geom[lev].periodicity());

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

  for (MFIter mfi(moments_lev, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    Array4<Real> moments_fab = moments_lev[mfi].array();
    Array4<IRL::SeparatorUnion> interface_fab = interface_lev[mfi].array();
    const Box& box = mfi.tilebox();

    if (case_name == "deformation3d" || case_name == "default") {
      amrex::launch(box, [=] AMREX_GPU_DEVICE(const Box& tbx) {
        Deformation3D::initialize_case(tbx, moments_fab, interface_fab, problo,
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
  if (case_name == "deformation3d" || case_name == "default") {
    max_vel = Deformation3D::get_max_vel();
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

inline IRL::Vec3<double> AmrCoreAdv::getVelocity(const IRL::Pt& pt,
                                                 Array4<Real const> const& vx,
                                                 Array4<Real const> const& vy,
                                                 Array4<Real const> const& vz,
                                                 const Box& bx, const int lev) {
  const auto& dx = geom[lev].CellSizeArray();
  const auto& prob_lo = geom[lev].ProbLoArray();
  const auto& lo = lbound(bx);
  const auto& hi = ubound(bx);
  // Find which cell the point falls in
  int i = static_cast<int>(amrex::Math::floor((pt[0] - prob_lo[0]) / dx[0]));
  int j = static_cast<int>(amrex::Math::floor((pt[1] - prob_lo[1]) / dx[1]));
  int k = static_cast<int>(amrex::Math::floor((pt[2] - prob_lo[2]) / dx[2]));
  if (!bx.contains(i, j, k)) {
    std::ostringstream oss;
    oss << "Position (" << pt[0] << "," << pt[1] << "," << pt[2]
        << ") leads to index (" << i << "," << j << "," << k
        << ") outside of box (" << lo.x << "," << lo.y << "," << lo.z << ")x("
        << hi.x << "," << hi.y << "," << hi.z << ")";
    throw std::runtime_error(oss.str());
  }
  // Cell lo face positions
  Real xlo = prob_lo[0] + i * dx[0];
  Real ylo = prob_lo[1] + j * dx[1];
  Real zlo = prob_lo[2] + k * dx[2];
  // Normalized position within cell [0,1]
  Real tx = (pt[0] - xlo) / dx[0];
  Real ty = (pt[1] - ylo) / dx[1];
  Real tz = (pt[2] - zlo) / dx[2];
  return IRL::Vec3<double>(vx(i, j, k) * (1.0 - tx) + vx(i + 1, j, k) * tx,
                           vy(i, j, k) * (1.0 - ty) + vy(i, j + 1, k) * ty,
                           vz(i, j, k) * (1.0 - tz) + vz(i, j, k + 1) * tz);
}

inline IRL::Pt AmrCoreAdv::project_vertex(const IRL::Pt& pt, const double dt,
                                          Array4<Real const> const& vx,
                                          Array4<Real const> const& vy,
                                          Array4<Real const> const& vz,
                                          const Box& bx, const int lev) {
  using Pt = IRL::Pt;
  auto v1 = getVelocity(pt, vx, vy, vz, bx, lev);
  auto v2 = getVelocity(pt + Pt::fromVec3(0.5 * dt * v1), vx, vy, vz, bx, lev);
  auto v3 = getVelocity(pt + Pt::fromVec3(0.5 * dt * v2), vx, vy, vz, bx, lev);
  auto v4 = getVelocity(pt + Pt::fromVec3(dt * v3), vx, vy, vz, bx, lev);
  return pt + Pt::fromVec3(dt * (v1 + 2.0 * v2 + 2.0 * v3 + v4) / 6.0);
}

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
  int lev = finest_level;
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
    Real dtdx = dt_lev / dx[0];
    Real dtdy = dt_lev / dx[1];
    Real dtdz = dt_lev / dx[2];

    // State with ghost cells
    MultiFab moments_with_ghost(grids[lev], dmap[lev], ncomp_moments, num_grow);
    MultiFab::Copy(moments_with_ghost, moments_old[lev], 0, 0,
                   moments_old[lev].nComp(), moments_old[lev].nGrow());
    moments_with_ghost.FillBoundary(geom[lev].periodicity());
    // FillPatch(lev, time, moments_with_ghost, 0,
    // moments_with_ghost.nComp());

    // Reconstruct interface
    for (MFIter mfi(interface[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> interface_lev = interface[lev].array(mfi);
      Array4<Real const> moments = moments_with_ghost.const_array(mfi);
      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (moments(i, j, k) * vol_inv <= IRL::global_constants::VF_LOW) {
          interface_lev(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
        } else if (moments(i, j, k) * vol_inv >=
                   IRL::global_constants::VF_HIGH) {
          interface_lev(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
        } else {
          //////////////////////////////////////// MOF
          // const double x = problo[0] + i * dx[0];
          // const double y = problo[1] + j * dx[1];
          // const double z = problo[2] + k * dx[2];
          // const IRL::Pt lower_cell_pt(x, y, z);
          // const IRL::Pt upper_cell_pt(x + dx[0], y + dx[1], z + dx[2]);
          // const IRL::Pt cell_center = 0.5 * (lower_cell_pt + upper_cell_pt);
          // const auto cell = IRL::RectangularCuboid::fromBoundingPts(
          //     lower_cell_pt, upper_cell_pt);
          // const double liq_m0 = moments(i, j, k, 0);
          // const double gas_m0 = vol - liq_m0;
          // const IRL::Pt liq_m1 =
          //     (1.0 / liq_m0) * IRL::Pt(moments(i, j, k, 1), moments(i, j, k,
          //     2),
          //                              moments(i, j, k, 3));
          // const IRL::Pt gas_m1 =
          //     (1.0 / gas_m0) *
          //     (vol * cell_center - IRL::Pt(moments(i, j, k, 1),
          //                                  moments(i, j, k, 2),
          //                                  moments(i, j, k, 3)));
          // const IRL::SeparatedMoments<IRL::VolumeMoments> svm(
          //     IRL::VolumeMoments(liq_m0, liq_m1),
          //     IRL::VolumeMoments(gas_m0, gas_m1));
          // interface_lev(i, j, k) =
          //     IRL::reconstructionWithMOF3D(cell, svm, 0.5, 0.5);
          ///////////////////////////////// LVIRA
          IRL::ELVIRANeighborhood elvira_neighborhood;
          IRL::LVIRANeighborhood<IRL::RectangularCuboid> lvira_neighborhood;
          elvira_neighborhood.resize(27);
          lvira_neighborhood.resize(27);
          IRL::RectangularCuboid cells[27];
          std::array<double, 27> cells_vfrac;
          for (int kk = k - 1; kk < k + 2; ++kk) {
            for (int jj = j - 1; jj < j + 2; ++jj) {
              for (int ii = i - 1; ii < i + 2; ++ii) {
                // Reversed order, bad for cache locality but thats okay..
                const int neigh_id =
                    (kk - k + 1) * 9 + (jj - j + 1) * 3 + (ii - i + 1);
                const double xx = problo[0] + ii * dx[0];
                const double yy = problo[1] + jj * dx[1];
                const double zz = problo[2] + kk * dx[2];
                cells[neigh_id] = IRL::RectangularCuboid::fromBoundingPts(
                    IRL::Pt(xx, yy, zz),
                    IRL::Pt(xx + dx[0], yy + dx[1], zz + dx[2]));
                cells_vfrac[neigh_id] = moments(ii, jj, kk, 0) * vol_inv;
                elvira_neighborhood.setMember(&cells[neigh_id],
                                              &cells_vfrac[neigh_id], ii - i,
                                              jj - j, kk - k);
                lvira_neighborhood.setMember(
                    static_cast<IRL::UnsignedIndex_t>(neigh_id),
                    &cells[neigh_id], &cells_vfrac[neigh_id]);
                // Trap center cell
                if (ii == i && jj == j && kk == k) {
                  lvira_neighborhood.setCenterOfStencil(neigh_id);
                }
              }
            }
          }
          // Now perform actual LVIRA and obtain interface PlanarSeparator
          interface_lev(i, j, k) =
              IRL::reconstructionWithELVIRA3D(elvira_neighborhood);
          const auto planar_separator = IRL::PlanarSeparator::fromOnePlane(
              interface_lev(i, j, k).getPlane());
          interface_lev(i, j, k) = IRL::reconstructionWithLVIRA3D(
              lvira_neighborhood, planar_separator);
        }
      });
    }  // end mfi

    SepUnionMultiFab interface_with_ghost(grids[lev], dmap[lev],
                                          interface[lev].nComp(), num_grow);
    for (MFIter mfi(interface[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion const> interface_lev =
          interface[lev].const_array(mfi);
      Array4<IRL::SeparatorUnion> tmp_interface =
          interface_with_ghost.array(mfi);
      amrex::ParallelFor(mfi.tilebox(),
                         [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                           tmp_interface(i, j, k) = interface_lev(i, j, k);
                         });
    }
    interface_with_ghost.FillBoundary(geom[lev].periodicity());

    for (MFIter mfi(moments_with_ghost, TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      Array4<Real const> velx = facevel[lev][0].const_array(mfi);
      Array4<Real const> vely = facevel[lev][1].const_array(mfi);
      Array4<Real const> velz = facevel[lev][2].const_array(mfi);
      Array4<IRL::SeparatorUnion const> interface_lev =
          interface_with_ghost.const_array(mfi);

      const Box& bx = mfi.tilebox();
      const Box& grown_bx = mfi.growntilebox();
      // const Box& gbx = amrex::grow(bx, 1);

      Array4<Real const> old_moments_lev = moments_old[lev].const_array(mfi);
      Array4<Real const> band_id_lev = band_id[lev].const_array(mfi);
      Array4<Real> new_moments_lev = moments_new[lev].array(mfi);

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (band_id_lev(i, j, k) != 0.0) {
          std::array<double, 6> flux_volumes;
          std::array<IRL::Pt, 8> cell;
          std::array<IRL::Pt, 14> preimage;
          IRL::CappedDodecahedron flux;
          const double x = problo[0] + i * dx[0];
          const double y = problo[1] + j * dx[1];
          const double z = problo[2] + k * dx[2];
          // Initialize cell corners
          cell[0] = IRL::Pt(x + dx[0], y, z + dx[2]);
          cell[1] = IRL::Pt(x + dx[0], y, z);
          cell[2] = IRL::Pt(x + dx[0], y + dx[1], z);
          cell[3] = IRL::Pt(x + dx[0], y + dx[1], z + dx[2]);
          cell[4] = IRL::Pt(x, y, z + dx[2]);
          cell[5] = IRL::Pt(x, y, z);
          cell[6] = IRL::Pt(x, y + dx[1], z);
          cell[7] = IRL::Pt(x, y + dx[1], z + dx[2]);
          // Compute the back projected cell corners
          for (int n = 0; n < 8; ++n) {
            preimage[n] = project_vertex(cell[n], -dt[lev], velx, vely, velz,
                                         grown_bx, lev);
          }
          // Compute soleinoidal flux volumes
          flux_volumes[0] = dt[lev] * velx(i, j, k) * dx[1] * dx[2];
          flux_volumes[1] = dt[lev] * velx(i + 1, j, k) * dx[1] * dx[2];
          flux_volumes[2] = dt[lev] * vely(i, j, k) * dx[0] * dx[2];
          flux_volumes[3] = dt[lev] * vely(i, j + 1, k) * dx[0] * dx[2];
          flux_volumes[4] = dt[lev] * velz(i, j, k) * dx[0] * dx[1];
          flux_volumes[5] = dt[lev] * velz(i, j, k + 1) * dx[0] * dx[1];
          // Create face flux hexahedra to compute correction
          for (int f = 0; f < 6; f++) {
            for (int i = 0; i < 4; i++) {
              flux[i] = cell[flux_id_table[f][i]];
              flux[i + 4] = preimage[flux_id_table[f][i]];
            }
            flux[8] =
                project_vertex(0.25 * (flux[0] + flux[1] + flux[2] + flux[3]),
                               -dt[lev], velx, vely, velz, grown_bx, lev);
            flux.adjustCapToMatchVolume(flux_volumes[f]);
            preimage[face_center_id_table[f]] = flux[8];
          }
          const auto preimage_cell =
              IRL::Polyhedron24::fromRawPtPointer(14, preimage.data());

          // Compute bounding box of preimage
          Real xlo = preimage[0][0];
          Real ylo = preimage[0][1];
          Real zlo = preimage[0][2];
          Real xhi = preimage[0][0];
          Real yhi = preimage[0][1];
          Real zhi = preimage[0][2];
          for (int i = 1; i < 14; i++) {
            xlo = preimage[i][0] < xlo ? preimage[i][0] : xlo;
            ylo = preimage[i][1] < ylo ? preimage[i][1] : ylo;
            zlo = preimage[i][2] < zlo ? preimage[i][2] : zlo;
            xhi = preimage[i][0] > xhi ? preimage[i][0] : xhi;
            yhi = preimage[i][1] > yhi ? preimage[i][1] : yhi;
            zhi = preimage[i][2] > zhi ? preimage[i][2] : zhi;
          }
          int ilo =
              static_cast<int>(amrex::Math::floor((xlo - problo[0]) / dx[0]));
          int jlo =
              static_cast<int>(amrex::Math::floor((ylo - problo[1]) / dx[1]));
          int klo =
              static_cast<int>(amrex::Math::floor((zlo - problo[2]) / dx[2]));
          int ihi =
              static_cast<int>(amrex::Math::floor((xhi - problo[0]) / dx[0]));
          int jhi =
              static_cast<int>(amrex::Math::floor((yhi - problo[1]) / dx[1]));
          int khi =
              static_cast<int>(amrex::Math::floor((zhi - problo[2]) / dx[2]));

          if (!grown_bx.contains(ilo, jlo, klo)) {
            std::cout << "Cell " << i << " " << j << " " << k << " has lo "
                      << ilo << " " << jlo << " " << klo << std::endl;
          }
          if (!grown_bx.contains(ihi, jhi, khi)) {
            std::cout << "Cell " << i << " " << j << " " << k << " has hi "
                      << ihi << " " << jhi << " " << khi << std::endl;
          }

          // Intersect preimage
          for (int n = 0; n < 4; n++) {
            new_moments_lev(i, j, k, n) = 0.0;
          }
          for (int kk = klo; kk <= khi; kk++) {
            for (int jj = jlo; jj <= jhi; jj++) {
              for (int ii = ilo; ii <= ihi; ii++) {
                const double xloc = problo[0] + ii * dx[0];
                const double yloc = problo[1] + jj * dx[1];
                const double zloc = problo[2] + kk * dx[2];
                const auto cell_loc = IRL::RectangularCuboid::fromBoundingPts(
                    IRL::Pt(xloc, yloc, zloc),
                    IRL::Pt(xloc + dx[0], yloc + dx[1], zloc + dx[2]));
                IRL::PlanarLocalizer localizer = cell_loc.getLocalizer();
                if (!(interface_lev(ii, jj, kk).type() ==
                          IRL::SeparatorUnion::SeparatorType::OnePlane ||
                      interface_lev(ii, jj, kk).type() ==
                          IRL::SeparatorUnion::SeparatorType::Paraboloid)) {
                  std::cout
                      << "Neigh " << ii << " " << jj << " " << kk
                      << " has type "
                      << static_cast<int>(interface_lev(ii, jj, kk).type())
                      << std::endl;
                }
                IRL::LocalizedSeparatorUnion local_sep(
                    &localizer, &interface_lev(ii, jj, kk));
                auto cut_moments = IRL::getVolumeMoments<IRL::VolumeMoments>(
                    preimage_cell, local_sep);
                const double m0 = cut_moments.volume();
                IRL::Pt centroid =
                    (1.0 / IRL::safelyEpsilon(m0)) * cut_moments.centroid();
                centroid[0] = std::max(centroid[0], xloc);
                centroid[0] = std::min(centroid[0], xloc + dx[0]);
                centroid[1] = std::max(centroid[1], yloc);
                centroid[1] = std::min(centroid[1], yloc + dx[1]);
                centroid[2] = std::max(centroid[2], zloc);
                centroid[2] = std::min(centroid[2], zloc + dx[2]);
                centroid = project_vertex(centroid, dt[lev], velx, vely, velz,
                                          grown_bx, lev);
                new_moments_lev(i, j, k, 0) += m0;
                new_moments_lev(i, j, k, 1) += m0 * centroid[0];
                new_moments_lev(i, j, k, 2) += m0 * centroid[1];
                new_moments_lev(i, j, k, 3) += m0 * centroid[2];
              }
            }
          }

          // Update moments
          // new_moments_lev(i, j, k) = old_moments_lev(i, j, k);
        }
      });
    }  // end mfi
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

    if (case_name == "deformation3d" || case_name == "default") {
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
    } else {
      throw std::runtime_error("Unknown case");
    }
  }
}