#ifndef EXAMPLES_AMREX_ADVECTOR_AMRCORE_H_
#define EXAMPLES_AMREX_ADVECTOR_AMRCORE_H_

#include <limits>
#include <memory>
#include <string>

#include <AMReX_AmrCore.H>
#include <AMReX_BCRec.H>
#include <AMReX_FluxRegister.H>

#include "irl/amrex/sepunion_multifab.h"
#include "irl/generic_cutting/cut_polygon.h"

#include "examples/amrex_advector/advection_helpers.h"

struct InterfaceScalarField {
  std::string name;

  amrex::MultiFab polygon_scalar_data;
  amrex::MultiFab paraboloid_scalar_data;

  std::vector<double> flattened_polygon_scalar_data;
  std::vector<double> flattened_paraboloid_scalar_data;

  InterfaceScalarField() = default;

  InterfaceScalarField(const std::string& a_name, const amrex::BoxArray& a_ba,
                       const amrex::DistributionMapping& a_dm,
                       const int a_ngrow = 0)
      : name(a_name),
        polygon_scalar_data(a_ba, a_dm, 1, a_ngrow),
        paraboloid_scalar_data(a_ba, a_dm, 1, a_ngrow) {
    polygon_scalar_data.setVal(0.0);
    paraboloid_scalar_data.setVal(0.0);
  }

  void clearFlattenedData() {
    flattened_polygon_scalar_data.clear();
    flattened_paraboloid_scalar_data.clear();
  }
};

class AmrCoreAdv : public amrex::AmrCore {
 public:
  ////////////////
  // public member functions

  // constructor - reads in parameters from inputs file
  //             - sizes multilevel arrays and data structures
  AmrCoreAdv();
  virtual ~AmrCoreAdv();

  // advance solution to final time
  void Evolve();

  // initializes multilevel data
  void InitData();

  amrex::DistributionMapping MakeDistributionMapWithWeights(
      int lev, const amrex::BoxArray& ba);

  void RedistributeLevel(int lev);

  virtual void regrid(int lbase, amrex::Real time,
                      bool initial = false) override;

  // Make a new level using provided BoxArray and DistributionMapping and
  // fill with interpolated coarse level data.
  // overrides the pure virtual function in AmrCore
  virtual void MakeNewLevelFromCoarse(
      int lev, amrex::Real time, const amrex::BoxArray& ba,
      const amrex::DistributionMapping& dm) override;

  // Remake an existing level using provided BoxArray and DistributionMapping
  // and fill with existing fine and coarse data. overrides the pure virtual
  // function in AmrCore
  virtual void RemakeLevel(int lev, amrex::Real time, const amrex::BoxArray& ba,
                           const amrex::DistributionMapping& dm) override;

  // Delete level data
  // overrides the pure virtual function in AmrCore
  virtual void ClearLevel(int lev) override;

  // Make a new level from scratch using provided BoxArray and
  // DistributionMapping. Only used during initialization. overrides the pure
  // virtual function in AmrCore
  virtual void MakeNewLevelFromScratch(
      int lev, amrex::Real time, const amrex::BoxArray& ba,
      const amrex::DistributionMapping& dm) override;

  // tag all cells for refinement
  // overrides the pure virtual function in AmrCore
  virtual void ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time,
                        int ngrow) override;

  void UpdateBand();

  // Advance phi at all levels for a single time step
  void AdvanceAllLevels(amrex::Real time, amrex::Real dt_lev, int iteration);

  // Define the advection velocity
  void DefineVelocityAtLevel(int lev, amrex::Real time);

  void DefineVelocityAllLevels(amrex::Real time);

  // compute dt from CFL considerations
  amrex::Real EstTimeStep(int lev, amrex::Real time);

  amrex::Real RecTime();
  amrex::Real AdvTime();

  void PostProcessCheckpointPair(const std::string& initial_checkpoint,
                                 const std::string& final_checkpoint);

 private:
  ////////////////
  // private member functions

  // read in some parameters from inputs file
  void ReadParameters();

  // set covered coarse cells to be the average of overlying fine cells
  void AverageDown();

  // more flexible version of AverageDown() that lets you average down across
  // multiple levels
  void AverageDownTo(int crse_lev);

  // compute a new multifab by coping in phi from valid region and filling ghost
  // cells works for single level and 2-level cases (fill fine grid ghost by
  // interpolating from coarse)
  template <typename MF>
  void FillPatch(int lev, amrex::Real time, MF& mf, int icomp, int ncomp);

  // fill an entire multifab by interpolating from the coarser level
  // this comes into play when a new level of refinement appears
  template <typename MF>
  void FillCoarsePatch(int lev, amrex::Real time, MF& mf, int icomp, int ncomp);

  // utility to copy in data from phi_old and/or phi_new into another multifab
  template <typename MF>
  void GetData(int lev, amrex::Real time, amrex::Vector<MF*>& data,
               amrex::Vector<amrex::Real>& datatime);

  // Advance a level by dt - includes a recursive call for finer levels
  void timeStepWithSubcycling(int lev, amrex::Real time, int iteration);

  // Advance all levels by the same dt
  void timeStepNoSubcycling(amrex::Real time, int iteration);

  // a wrapper for EstTimeStep
  void ComputeDt();

  // get plotfile name
  std::string PlotFileName(int lev) const;

  // put together an array of multifabs for writing
  amrex::Vector<const amrex::MultiFab*> PlotFileMF() const;

  // set plotfile variables names
  amrex::Vector<std::string> PlotFileVarNames() const;

  // write plotfile to disk
  void WritePlotFile();

  // write checkpoint file to disk
  void WriteCheckpointFile() const;

  // read checkpoint file from disk
  void ReadCheckpointFile();

  void GetReconstruction(const int lev);
  void GetReconstruction(
      amrex::SepUnionMultiFab& interface,
      amrex::SepUnionMultiFab& interface_with_ghost,
      const amrex::MultiFab& moments, const amrex::Geometry& a_geom,
      std::vector<InterfaceScalarField>* scalar_fields = nullptr);

  void TransportMoments(
      const amrex::SepUnionMultiFab& a_interface_with_ghost,
      const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& a_facevel,
      const amrex::MultiFab& a_band_id, amrex::MultiFab& a_moments,
      const amrex::Geometry& a_geom, const double a_dt, const double a_time);

  void BuildUniformFinestMoments(amrex::MultiFab& a_uniform_moments) const;

  amrex::Real ComputeCompositeM0() const;

  void ComputeUniformMomentL1Errors(const amrex::MultiFab& a_initial,
                                    const amrex::MultiFab& a_final) const;

  amrex::Real ComputeL1ErrorM0() const;

  ////////////////
  // private data members

  amrex::Vector<int> istep;  // which step?

  // keep track of old time, new time, and time step at each level
  amrex::Vector<amrex::Real> t_new;
  amrex::Vector<amrex::Real> t_old;
  amrex::Vector<amrex::Real> dt;

  // array of multifabs to store the solution at each level of refinement
  // after advancing a level we use "swap".
  int ncomp_moments = 4;
  amrex::Vector<amrex::MultiFab> moments_new;
  amrex::Vector<amrex::MultiFab> moments_old;
  amrex::Vector<amrex::MultiFab> band_id;
  amrex::Vector<amrex::SepUnionMultiFab> interface;

  // Scalar data attached to reconstructed interface surfaces.
  amrex::Vector<std::vector<InterfaceScalarField>> interface_scalar_fields;

  // this is essentially a 2*DIM integer array storing the physical boundary
  // condition types at the lo/hi walls in each direction
  amrex::Vector<amrex::BCRec> bcs;  // 1-component

  // stores fluxes at coarse-fine interface for synchronization
  // this will be sized "nlevs_max+1"
  // NOTE: the flux register associated with flux_reg[lev] is associated
  // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
  // therefore flux_reg[0] and flux_reg[nlevs_max] are never actually
  // used in the reflux operation
  amrex::Vector<std::unique_ptr<amrex::FluxRegister>> flux_reg;

  // Velocity on all faces at all levels
  amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> facevel;

  ////////////////
  // runtime parameters

  // maximum number of steps and stop time
  int max_step = std::numeric_limits<int>::max();
  amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

  // if >= 0 we restart from a checkpoint
  std::string restart_chkfile = "";
  std::string checkpoint_path = "";
  std::string interface_output_path = "";

  // case name
  std::string case_name = "default";

  VelocityFieldType velocity_field_type = VelocityFieldType::Interpolated;

  // reconstruction name
  std::string reconstruction_name = "default";

  // advection name
  std::string advection_name = "default";

  // advective cfl number - dt = cfl*dx/umax
  amrex::Real cfl = 0.7;

  // timers
  amrex::Real reconstruction_time = 0.0;
  amrex::Real advection_time = 0.0;

  // how often each level regrids the higher levels of refinement
  // (after a level advances that many time steps)
  int regrid_int = 2;

  // hyperbolic refluxing as part of multilevel synchronization
  int do_reflux = 1;

  // plotfile prefix and frequency
  std::string plot_file{"plt"};
  int plot_int = -1;

  // checkpoint prefix and frequency
  std::string chk_file{"chk"};
  int chk_int = -1;

  // Number of ghost layers needed for advection
  int num_grow = 1;

  amrex::Real initial_liquid_mass = 0.0;
  amrex::MultiFab uniform_initial_moments;
  amrex::MultiFab initial_moments;
};

#include "examples/amrex_advector/amrcore_advector.tpp"

#include "examples/amrex_advector/advection.tpp"
#include "examples/amrex_advector/reconstruction.tpp"

#endif  // EXAMPLES_AMREX_ADVECTOR_AMRCORE_H_
