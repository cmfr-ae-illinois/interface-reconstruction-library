
#include <AMReX_Geometry.H>
#include "irl/amrex/sepunion_multifab.h"

using namespace amrex;

extern "C" {

void amrex_fi_new_sepunionmfab(SepUnionMultiFab*& mf, const BoxArray*& ba,
                               const DistributionMapping*& dm, int nc,
                               const int* ng, const int* nodal) {
  mf = new SepUnionMultiFab(amrex::convert(*ba, IntVect(nodal)), *dm, nc,
                            IntVect(ng));
  mf->setVal(IRL::SeparatorUnion());
  ba = &(mf->boxArray());
  dm = &(mf->DistributionMap());
}

void amrex_fi_new_sepunionmfab_alias(SepUnionMultiFab*& mf,
                                     const SepUnionMultiFab* srcmf, int comp,
                                     int ncomp) {
  mf = new SepUnionMultiFab(*srcmf, amrex::make_alias, comp, ncomp);
  mf->setVal(IRL::SeparatorUnion());
}

void amrex_fi_delete_sepunionmfab(SepUnionMultiFab* mf) { delete mf; }

int amrex_fi_sepunionmfab_ncomp(const SepUnionMultiFab* mf) {
  return mf->nComp();
}

void amrex_fi_sepunionmfab_ngrow(const SepUnionMultiFab* mf, int* ngv) {
  IntVect const& ng = mf->nGrowVect();
  AMREX_D_TERM(ngv[0] = ng[0];, ngv[1] = ng[1];, ngv[2] = ng[2]);
}

const BoxArray* amrex_fi_sepunionmfab_boxarray(const SepUnionMultiFab* mf) {
  return &(mf->boxArray());
}

const DistributionMapping* amrex_fi_sepunionmfab_distromap(
    const SepUnionMultiFab* mf) {
  return &(mf->DistributionMap());
}

void amrex_fi_sepunionmfab_dataptr_iter(SepUnionMultiFab* mf, MFIter* mfi,
                                        IRL::SeparatorUnion*& dp, int lo[3],
                                        int hi[3]) {
  SepUnionArrayBox& fab = (*mf)[*mfi];
  dp = fab.dataPtr();
  const Box& bx = fab.box();
  const int* lov = bx.loVect();
  const int* hiv = bx.hiVect();
  for (int i = 0; i < BL_SPACEDIM; ++i) {
    lo[i] = lov[i];
    hi[i] = hiv[i];
  }
}

void amrex_fi_sepunionmfab_dataptr_int(SepUnionMultiFab* mf, int igrd,
                                       IRL::SeparatorUnion*& dp, int lo[3],
                                       int hi[3]) {
  SepUnionArrayBox& fab = (*mf)[igrd];
  dp = fab.dataPtr();
  const Box& bx = fab.box();
  const int* lov = bx.loVect();
  const int* hiv = bx.hiVect();
  for (int i = 0; i < BL_SPACEDIM; ++i) {
    lo[i] = lov[i];
    hi[i] = hiv[i];
  }
}

void amrex_fi_sepunionmfab_copy(SepUnionMultiFab* dstmf,
                                const SepUnionMultiFab* srcmf, int srccomp,
                                int dstcomp, int nc, const int* ng) {
  SepUnionMultiFab::Copy(*dstmf, *srcmf, srccomp, dstcomp, nc, IntVect(ng));
}

void amrex_fi_sepunionmfab_parallelcopy(SepUnionMultiFab* dstmf,
                                        const SepUnionMultiFab* srcmf,
                                        int srccomp, int dstcomp, int nc,
                                        int srcng, int dstng,
                                        const Geometry* geom) {
  dstmf->ParallelCopyWithPeriodicShift(*srcmf, srccomp, dstcomp, nc, srcng,
                                       dstng, *geom);
}

void amrex_fi_sepunionmfab_parallelcopy_gv(SepUnionMultiFab* dstmf,
                                           const SepUnionMultiFab* srcmf,
                                           int srccomp, int dstcomp, int nc,
                                           const int* srcng, const int* dstng,
                                           const Geometry* geom) {
  IntVect sg(AMREX_D_DECL(srcng[0], srcng[1], srcng[2]));
  IntVect dg(AMREX_D_DECL(dstng[0], dstng[1], dstng[2]));
  dstmf->ParallelCopyWithPeriodicShift(*srcmf, srccomp, dstcomp, nc, sg, dg,
                                       *geom);
}

void amrex_fi_sepunionmfab_fill_boundary(SepUnionMultiFab* mf,
                                         const Geometry* geom, int c, int nc,
                                         int cross) {
  mf->FillBoundaryWithPeriodicShift(c, nc, *geom, cross);
}

void amrex_fi_write_sepunionmfab(const SepUnionMultiFab* mf, const char* name) {
  SepUnionMultiFab_Write(*mf, std::string(name));
}

void amrex_fi_read_sepunionmfab(SepUnionMultiFab* mf, const char* name) {
  BL_ASSERT(mf != nullptr);
  SepUnionMultiFab_Read(*mf, std::string(name));
}

}  // extern "C"