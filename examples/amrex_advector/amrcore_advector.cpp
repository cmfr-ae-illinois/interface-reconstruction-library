
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

// #include <Kernels.H>
#include "examples/amrex_advector/amrcore_advector.h"

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_psi(amrex::Box const& bx, const amrex::Real time,
                           amrex::Array4<amrex::Real> const& psi,
                           amrex::GeometryData const& geomdata) {
  using namespace amrex;
  constexpr Real PI = 3.1415926535897932384626;

  const auto lo = lbound(bx);
  const auto hi = ubound(bx);

  const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
  const Real* AMREX_RESTRICT dx = geomdata.CellSize();

  for (int j = lo.y; j <= hi.y; ++j) {
    Real y = dx[1] * (0.5 + j) + prob_lo[1];
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
      Real x = dx[0] * (0.5 + i) + prob_lo[0];
      psi(i, j, 0) = std::pow(std::sin(PI * x), 2) *
                     std::pow(std::sin(PI * y), 2) * std::cos(PI * time / 2.0) *
                     1.0 / PI;
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_x(int i, int j, int k,
                         amrex::Array4<amrex::Real> const& vx,
                         amrex::Array4<amrex::Real> const& psi,
                         amrex::Real dy) {
  vx(i, j, k) = -((psi(i, j + 1, 0) + psi(i - 1, j + 1, 0)) -
                  (psi(i, j - 1, 0) + psi(i - 1, j - 1, 0))) *
                (0.25 / dy);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_y(int i, int j, int k,
                         amrex::Array4<amrex::Real> const& vy,
                         amrex::Array4<amrex::Real> const& psi,
                         amrex::Real dx) {
  vy(i, j, k) = ((psi(i + 1, j, 0) + psi(i + 1, j - 1, 0)) -
                 (psi(i - 1, j, 0) + psi(i - 1, j - 1, 0))) *
                (0.25 / dx);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_z(int i, int j, int k,
                         amrex::Array4<amrex::Real> const& vz) {
  vz(i, j, k) = 0.0;
}

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
  nsubsteps.resize(nlevs_max, 1);
  if (do_subcycle) {
    for (int lev = 1; lev <= max_level; ++lev) {
      nsubsteps[lev] = MaxRefRatio(lev - 1);
    }
  }

  t_new.resize(nlevs_max, 0.0);
  t_old.resize(nlevs_max, -1.e100);
  dt.resize(nlevs_max, 1.e100);

  moments_new.resize(nlevs_max);
  moments_old.resize(nlevs_max);

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
    // lo-side BCs
    if (bc_lo[idim] == BCType::int_dir ||  // periodic uses "internal Dirichlet"
        bc_lo[idim] == BCType::foextrap ||  // first-order extrapolation
        bc_lo[idim] == BCType::ext_dir) {   // external Dirichlet
      bcs[0].setLo(idim, bc_lo[idim]);
    } else {
      amrex::Abort("Invalid bc_lo");
    }

    // hi-side BCSs
    if (bc_hi[idim] == BCType::int_dir ||  // periodic uses "internal Dirichlet"
        bc_hi[idim] == BCType::foextrap ||  // first-order extrapolation
        bc_hi[idim] == BCType::ext_dir) {   // external Dirichlet
      bcs[0].setHi(idim, bc_hi[idim]);
    } else {
      amrex::Abort("Invalid bc_hi");
    }
  }

  // stores fluxes at coarse-fine interface for synchronization
  // this will be sized "nlevs_max+1"
  // NOTE: the flux register associated with flux_reg[lev] is associated
  // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
  // therefore flux_reg[0] is never actually used in the reflux operation
  flux_reg.resize(nlevs_max + 1);
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

#ifdef AMREX_MEM_PROFILING
    {
      std::ostringstream ss;
      ss << "[STEP " << step + 1 << "]";
      MemProfiler::report(ss.str());
    }
#endif

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

  t_new[lev] = time;
  t_old[lev] = time - 1.e200;

  // This clears the old MultiFab and allocates the new one
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    facevel[lev][idim] = MultiFab(
        amrex::convert(ba, IntVect::TheDimensionVector(idim)), dm, 1, 1);
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

  FillPatch(lev, time, new_state, 0, ncomp);

  std::swap(new_state, moments_new[lev]);
  std::swap(old_state, moments_old[lev]);

  t_new[lev] = time;
  t_old[lev] = time - 1.e200;

  // This clears the old MultiFab and allocates the new one
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    facevel[lev][idim] = MultiFab(
        amrex::convert(ba, IntVect::TheDimensionVector(idim)), dm, 1, 1);
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
  flux_reg[lev].reset(nullptr);
}

// Make a new level from scratch using provided BoxArray and
// DistributionMapping. Only used during initialization. overrides the pure
// virtual function in AmrCore
void AmrCoreAdv::MakeNewLevelFromScratch(int lev, Real time, const BoxArray& ba,
                                         const DistributionMapping& dm) {
  const int ncomp = 1;
  const int nghost = 0;

  moments_new[lev].define(ba, dm, ncomp, nghost);
  moments_old[lev].define(ba, dm, ncomp, nghost);

  t_new[lev] = time;
  t_old[lev] = time - 1.e200;

  // This clears the old MultiFab and allocates the new one
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    facevel[lev][idim] = MultiFab(
        amrex::convert(ba, IntVect::TheDimensionVector(idim)), dm, 1, 1);
  }

  if (lev > 0 && do_reflux) {
    flux_reg[lev].reset(
        new FluxRegister(ba, dm, refRatio(lev - 1), lev, ncomp));
  }

  MultiFab& state = moments_new[lev];

  const auto problo = Geom(lev).ProbLoArray();
  const auto dx = Geom(lev).CellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    Array4<Real> fab = state[mfi].array();
    const Box& box = mfi.tilebox();

    amrex::launch(box, [=] AMREX_GPU_DEVICE(Box const& tbx) {
      //   initdata(tbx, fab, problo, dx);
      const auto lo = lbound(tbx);
      const auto hi = ubound(tbx);

      for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
          Real y = problo[1] + (0.5 + j) * dx[1];
          for (int i = lo.x; i <= hi.x; ++i) {
            Real x = problo[0] + (0.5 + i) * dx[0];
            Real r2 = (std::pow(x - 0.5, 2) + std::pow((y - 0.75), 2)) / 0.01;
            fab(i, j, k) = 1.0 + std::exp(-r2);
          }
        }
      }
    });
  }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void AmrCoreAdv::ErrorEst(int lev, TagBoxArray& tags, Real /*time*/,
                          int /*ngrow*/) {
  static bool first = true;
  static Vector<Real> phierr;

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

  //    const int clearval = TagBox::CLEAR;
  const int tagval = TagBox::SET;

  const MultiFab& state = moments_new[lev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      const auto statefab = state.array(mfi);
      const auto tagfab = tags.array(mfi);
      //   Real phierror = phierr[lev];

      amrex::ParallelFor(bx,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                           // state_error(i, j, k, tagfab, statefab, phierror,
                           // tagval);
                           if (statefab(i, j, k) > 1.05) {
                             tagfab(i, j, k) = tagval;
                           }
                         });
    }
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
    pp.query("do_reflux", do_reflux);
    pp.query("do_subcycle", do_subcycle);
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
void AmrCoreAdv::FillPatch(int lev, Real time, MultiFab& mf, int icomp,
                           int ncomp) {
  if (lev == 0) {
    Vector<MultiFab*> smf;
    Vector<Real> stime;
    GetData(0, time, smf, stime);

    if (Gpu::inLaunchRegion()) {
      GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
      PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > physbc(geom[lev], bcs,
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
    Vector<MultiFab*> cmf, fmf;
    Vector<Real> ctime, ftime;
    GetData(lev - 1, time, cmf, ctime);
    GetData(lev, time, fmf, ftime);

    Interpolater* mapper = &cell_cons_interp;

    if (Gpu::inLaunchRegion()) {
      GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
      PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[lev - 1], bcs,
                                                         gpu_bndry_func);
      PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[lev], bcs,
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
void AmrCoreAdv::FillCoarsePatch(int lev, Real time, MultiFab& mf, int icomp,
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
    PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[lev - 1], bcs,
                                                       gpu_bndry_func);
    PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[lev], bcs,
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

// Advance all the levels with the same dt
void AmrCoreAdv::timeStepNoSubcycling(Real time, int iteration) {
  if (max_level > 0 && regrid_int > 0)  // We may need to regrid
  {
    if (istep[0] % regrid_int == 0) {
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
  int n_factor = 1;

  for (int lev = 0; lev <= finest_level; ++lev) {
    dt_tmp[lev] = std::min(dt_tmp[lev], change_max * dt[lev]);
    n_factor *= nsubsteps[lev];
    dt_0 = std::min(dt_0, n_factor * dt_tmp[lev]);
  }

  // Limit dt's by the value of stop_time.
  const Real eps = 1.e-3 * dt_0;

  if (t_new[0] + dt_0 > stop_time - eps) {
    dt_0 = stop_time - t_new[0];
  }

  dt[0] = dt_0;

  for (int lev = 1; lev <= finest_level; ++lev) {
    dt[lev] = dt[lev - 1] / nsubsteps[lev];
  }
}

// compute dt from CFL considerations
Real AmrCoreAdv::EstTimeStep(int lev, Real time) {
  BL_PROFILE("AmrCoreAdv::EstTimeStep()");

  Real dt_est = std::numeric_limits<Real>::max();

  const Real* dx = geom[lev].CellSize();

  if (time == 0.0) {
    DefineVelocityAtLevel(lev, time);
  } else {
    Real t_nph_predicted = time + 0.5 * dt[lev];
    DefineVelocityAtLevel(lev, t_nph_predicted);
  }

  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    Real est = facevel[lev][idim].norminf(0, 0, true);
    dt_est = amrex::min(dt_est, dx[idim] / est);
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
  Vector<const MultiFab*> r;
  for (int i = 0; i <= finest_level; ++i) {
    r.push_back(&moments_new[i]);
  }
  return r;
}

// set plotfile variable names
Vector<std::string> AmrCoreAdv::PlotFileVarNames() const { return {"phi"}; }

// write plotfile to disk
void AmrCoreAdv::WritePlotFile() const {
  const std::string& plotfilename = PlotFileName(istep[0]);
  const auto& mf = PlotFileMF();
  const auto& varnames = PlotFileVarNames();

  amrex::Print() << "Writing plotfile " << plotfilename << "\n";

  amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, mf, varnames,
                                 Geom(), t_new[0], istep, refRatio());
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
  // ---- dirName is built first.  if dirName exists, it is renamed.  then build
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
    VisMF::Write(moments_new[lev], amrex::MultiFabFileFullPrefix(
                                       lev, checkpointname, "Level_", "phi"));
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

    // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    // build MultiFab and FluxRegister data
    int ncomp = 1;
    int nghost = 0;
    moments_old[lev].define(grids[lev], dmap[lev], ncomp, nghost);
    moments_new[lev].define(grids[lev], dmap[lev], ncomp, nghost);

    if (lev > 0 && do_reflux) {
      flux_reg[lev].reset(new FluxRegister(grids[lev], dmap[lev],
                                           refRatio(lev - 1), lev, ncomp));
    }

    // build face velocity MultiFabs
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      facevel[lev][idim] = MultiFab(
          amrex::convert(ba, IntVect::TheDimensionVector(idim)), dm, 1, 1);
    }
  }

  // read in the MultiFab data
  for (int lev = 0; lev <= finest_level; ++lev) {
    VisMF::Read(moments_new[lev], amrex::MultiFabFileFullPrefix(
                                      lev, restart_chkfile, "Level_", "phi"));
  }
}

void AmrCoreAdv::AdvanceAllLevels(Real time, Real dt_lev, int /*iteration*/) {
  constexpr int num_grow = 3;

  Vector<Array<MultiFab, AMREX_SPACEDIM> > fluxes(finest_level + 1);
  for (int lev = 0; lev <= finest_level; lev++) {
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      BoxArray ba = grids[lev];
      ba.surroundingNodes(idim);
      fluxes[lev][idim] = MultiFab(ba, dmap[lev], 1, 0);
    }
  }

  for (int lev = 0; lev <= finest_level; lev++) {
    std::swap(moments_old[lev], moments_new[lev]);
    t_old[lev] = t_new[lev];
    t_new[lev] += dt_lev;

    const auto dx = geom[lev].CellSizeArray();
    AMREX_D_TERM(Real dtdx = dt_lev / dx[0];, Real dtdy = dt_lev / dx[1];
                 , Real dtdz = dt_lev / dx[2]);

    // State with ghost cells
    MultiFab Sborder(grids[lev], dmap[lev], moments_new[lev].nComp(), num_grow);
    FillPatch(lev, time, Sborder, 0, Sborder.nComp());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      FArrayBox tmpfab;
      for (MFIter mfi(moments_new[lev], TilingIfNotGPU()); mfi.isValid();
           ++mfi) {
        AMREX_D_TERM(
            Array4<Real const> velx = facevel[lev][0].const_array(mfi);
            , Array4<Real const> vely = facevel[lev][1].const_array(mfi);
            , Array4<Real const> velz = facevel[lev][2].const_array(mfi));

        const Box& bx = mfi.tilebox();
        const Box& gbx = amrex::grow(bx, 1);

        Array4<Real const> statein = Sborder.const_array(mfi);

        AMREX_D_TERM(Array4<Real> fluxx = fluxes[lev][0].array(mfi);
                     , Array4<Real> fluxy = fluxes[lev][1].array(mfi);
                     , Array4<Real> fluxz = fluxes[lev][2].array(mfi));

        int ntmpcomps = (AMREX_SPACEDIM == 2) ? 2 : 3;
        tmpfab.resize(amrex::grow(bx, 2), ntmpcomps);
        Elixir tmpeli = tmpfab.elixir();
        int itmp = 0;

        // compute longitudinal fluxes
        // ===========================

        // x -------------------------
        Array4<Real> phix = tmpfab.array(itmp);
        itmp += 1;
        Box b = gbx;
        Array4<Real const> phix_c = phix;
        amrex::ParallelFor(
            b.grow(Direction::x, -1).surroundingNodes(Direction::x),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) {
              phix(i, j, k) = ((velx(i, j, k) < 0) ? statein(i, j, k)
                                                   : statein(i - 1, j, k));
            });

        // y -------------------------
        Array4<Real> phiy = tmpfab.array(itmp);
        itmp += 1;
        b = gbx;
        Array4<Real const> phiy_c = phiy;
        amrex::ParallelFor(
            b.grow(Direction::y, -1).surroundingNodes(Direction::y),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) {
              phiy(i, j, k) = ((vely(i, j, k) < 0) ? statein(i, j, k)
                                                   : statein(i, j - 1, k));
            });

#if (AMREX_SPACEDIM > 2)
        // z -------------------------
        Array4<Real> phiz = tmpfab.array(itmp);
        itmp += 1;
        b = gbx;
        Array4<Real const> phiz_c = phiz;
        amrex::ParallelFor(
            b.grow(Direction::z, -1).surroundingNodes(Direction::z),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) {
              phiz(i, j, k) = ((velz(i, j, k) < 0) ? statein(i, j, k)
                                                   : statein(i, j, k - 1));
            });
#endif

        // final edge states
        // ===========================
        amrex::ParallelFor(mfi.nodaltilebox(0),
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                             fluxx(i, j, k) = velx(i, j, k) * phix(i, j, k);
                           });

        amrex::ParallelFor(mfi.nodaltilebox(1),
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                             fluxy(i, j, k) = vely(i, j, k) * phiy(i, j, k);
                           });
#if (AMREX_SPACEDIM == 3)
        amrex::ParallelFor(mfi.nodaltilebox(2),
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                             fluxz(i, j, k) = velz(i, j, k) * phiz(i, j, k);
                           });
#endif
        AMREX_ASSERT(itmp == ntmpcomps);
      }  // end mfi
    }  // end omp
  }  // end lev

  // =======================================================
  // Average down the fluxes before using them to update phi
  // =======================================================
  for (int lev = finest_level; lev > 0; lev--) {
    average_down_faces(amrex::GetArrOfConstPtrs(fluxes[lev]),
                       amrex::GetArrOfPtrs(fluxes[lev - 1]), refRatio(lev - 1),
                       Geom(lev - 1));
  }

  for (int lev = 0; lev <= finest_level; lev++) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      const auto dx = geom[lev].CellSizeArray();
      AMREX_D_TERM(Real dtdx = dt_lev / dx[0];, Real dtdy = dt_lev / dx[1];
                   , Real dtdz = dt_lev / dx[2]);

      // ===========================================
      // Compute moments_new using a conservative update
      // ===========================================
      for (MFIter mfi(moments_new[lev], TilingIfNotGPU()); mfi.isValid();
           ++mfi) {
        Array4<Real const> statein = moments_old[lev].const_array(mfi);
        Array4<Real> stateout = moments_new[lev].array(mfi);

        AMREX_D_TERM(
            Array4<Real const> fluxx = fluxes[lev][0].const_array(mfi);
            , Array4<Real const> fluxy = fluxes[lev][1].const_array(mfi);
            , Array4<Real const> fluxz = fluxes[lev][2].const_array(mfi));

        const Box& bx = mfi.tilebox();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          stateout(i, j, k) =
              statein(i, j, k) +
              (AMREX_D_TERM((fluxx(i, j, k) - fluxx(i + 1, j, k)) * dtdx,
                            +(fluxy(i, j, k) - fluxy(i, j + 1, k)) * dtdy,
                            +(fluxz(i, j, k) - fluxz(i, j, k + 1)) * dtdz));
        });
      }  // end mfi
    }  // end omp
  }  // end lev
}

void AmrCoreAdv::DefineVelocityAllLevels(Real time) {
  for (int lev = 0; lev <= finest_level; ++lev)
    DefineVelocityAtLevel(lev, time);
}

void AmrCoreAdv::DefineVelocityAtLevel(int lev, Real time) {
  const auto dx = geom[lev].CellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    for (MFIter mfi(moments_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      // ======== GET FACE VELOCITY =========
      GpuArray<Box, AMREX_SPACEDIM> nbx;
      AMREX_D_TERM(nbx[0] = mfi.nodaltilebox(0);, nbx[1] = mfi.nodaltilebox(1);
                   , nbx[2] = mfi.nodaltilebox(2););

      AMREX_D_TERM(const Box& ngbxx = amrex::grow(mfi.nodaltilebox(0), 1);
                   , const Box& ngbxy = amrex::grow(mfi.nodaltilebox(1), 1);
                   , const Box& ngbxz = amrex::grow(mfi.nodaltilebox(2), 1););

      GpuArray<Array4<Real>, AMREX_SPACEDIM> vel{
          AMREX_D_DECL(facevel[lev][0].array(mfi), facevel[lev][1].array(mfi),
                       facevel[lev][2].array(mfi))};

      const Box& psibox =
          Box(IntVect(AMREX_D_DECL(
                  std::min(ngbxx.smallEnd(0) - 1, ngbxy.smallEnd(0) - 1),
                  std::min(ngbxx.smallEnd(1) - 1, ngbxy.smallEnd(1) - 1), 0)),
              IntVect(AMREX_D_DECL(
                  std::max(ngbxx.bigEnd(0), ngbxy.bigEnd(0) + 1),
                  std::max(ngbxx.bigEnd(1) + 1, ngbxy.bigEnd(1)), 0)));

      FArrayBox psifab(psibox, 1);
      Elixir psieli = psifab.elixir();
      Array4<Real> psi = psifab.array();
      GeometryData geomdata = geom[lev].data();

      amrex::launch(psibox, [=] AMREX_GPU_DEVICE(const Box& tbx) {
        get_face_velocity_psi(tbx, time, psi, geomdata);
      });

      amrex::ParallelFor(AMREX_D_DECL(ngbxx, ngbxy, ngbxz),
                         AMREX_D_DECL(
                             [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                               get_face_velocity_x(i, j, k, vel[0], psi, dx[1]);
                             },
                             [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                               get_face_velocity_y(i, j, k, vel[1], psi, dx[0]);
                             },
                             [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                               get_face_velocity_z(i, j, k, vel[2]);
                             }));
    }
  }
}