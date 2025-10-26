
#ifndef EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_RECONSTRUCTION_METRICS_TPP_
#define EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_RECONSTRUCTION_METRICS_TPP_

#include "examples/implicit_surface_reconstruction/reconstruction_metrics.h"

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
MomentDiffNorms compute_moment_diff_norms(
    const Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                         IRL::GeneralSurfaceMoments3D<SM_ORDER>>>& A,
    const Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                         IRL::GeneralSurfaceMoments3D<SM_ORDER>>>& B) {
  const BasicMesh& mA = A.getMesh();
  const BasicMesh& mB = B.getMesh();
  if (mA.imin() != mB.imin() || mA.imax() != mB.imax() ||
      mA.jmin() != mB.jmin() || mA.jmax() != mB.jmax() ||
      mA.kmin() != mB.kmin() || mA.kmax() != mB.kmax()) {
    throw std::runtime_error("Meshes do not match for A and B.");
  }

  // Accumulators
  MomentDiffNorms out{};
  out.vol_M0_Linf = out.vol_M1_Linf = out.vol_M2_Linf = 0.0;
  out.surf_M0_Linf = out.surf_M1_Linf = out.surf_M2_Linf = 0.0;

  long long Ncells = 0;
  long double vol_M0_L2_acc = 0.0L, vol_M1_L2_acc = 0.0L, vol_M2_L2_acc = 0.0L;
  long double surf_M0_L2_acc = 0.0L, surf_M1_L2_acc = 0.0L,
              surf_M2_L2_acc = 0.0L;

  for (int i = mA.imin(); i <= mA.imax(); ++i) {
    for (int j = mA.jmin(); j <= mA.jmax(); ++j) {
      for (int k = mA.kmin(); k <= mA.kmax(); ++k) {
        const auto& a = A(i, j, k);
        const auto& b = B(i, j, k);

        // M0
        const double dM0v = std::abs(a.first[M0] - b.first[M0]);
        out.vol_M0_Linf = std::max(out.vol_M0_Linf, dM0v);
        vol_M0_L2_acc += static_cast<long double>(dM0v) * dM0v;

        // M1
        const double dM1xv = a.first[M1x] - b.first[M1x];
        const double dM1yv = a.first[M1y] - b.first[M1y];
        const double dM1zv = a.first[M1z] - b.first[M1z];
        const double nM1v =
            std::sqrt(dM1xv * dM1xv + dM1yv * dM1yv + dM1zv * dM1zv);
        out.vol_M1_Linf = std::max(out.vol_M1_Linf, nM1v);
        vol_M1_L2_acc += static_cast<long double>(nM1v) * nM1v;

        // M2
        const double dMxxv = a.first[Mxx] - b.first[Mxx];
        const double dMxyv = a.first[Mxy] - b.first[Mxy];
        const double dMxzv = a.first[Mxz] - b.first[Mxz];
        const double dMyyv = a.first[Myy] - b.first[Myy];
        const double dMyzv = a.first[Myz] - b.first[Myz];
        const double dMzzv = a.first[Mzz] - b.first[Mzz];
        const double nM2v_sq =
            dMxxv * dMxxv + dMyyv * dMyyv + dMzzv * dMzzv +
            2.0 * (dMxyv * dMxyv + dMxzv * dMxzv + dMyzv * dMyzv);
        const double nM2v = std::sqrt(nM2v_sq);
        out.vol_M2_Linf = std::max(out.vol_M2_Linf, nM2v);
        vol_M2_L2_acc += static_cast<long double>(nM2v_sq);

        // surface
        const double dM0s = std::abs(a.second[M0] - b.second[M0]);
        out.surf_M0_Linf = std::max(out.surf_M0_Linf, dM0s);
        surf_M0_L2_acc += static_cast<long double>(dM0s) * dM0s;

        const double dM1xs = a.second[M1x] - b.second[M1x];
        const double dM1ys = a.second[M1y] - b.second[M1y];
        const double dM1zs = a.second[M1z] - b.second[M1z];
        const double nM1s =
            std::sqrt(dM1xs * dM1xs + dM1ys * dM1ys + dM1zs * dM1zs);
        out.surf_M1_Linf = std::max(out.surf_M1_Linf, nM1s);
        surf_M1_L2_acc += static_cast<long double>(nM1s) * nM1s;

        const double dMxxs = a.second[Mxx] - b.second[Mxx];
        const double dMxys = a.second[Mxy] - b.second[Mxy];
        const double dMxzs = a.second[Mxz] - b.second[Mxz];
        const double dMyys = a.second[Myy] - b.second[Myy];
        const double dMyzs = a.second[Myz] - b.second[Myz];
        const double dMzzs = a.second[Mzz] - b.second[Mzz];
        const double nM2s_sq =
            dMxxs * dMxxs + dMyys * dMyys + dMzzs * dMzzs +
            2.0 * (dMxys * dMxys + dMxzs * dMxzs + dMyzs * dMyzs);
        const double nM2s = std::sqrt(nM2s_sq);
        out.surf_M2_Linf = std::max(out.surf_M2_Linf, nM2s);
        surf_M2_L2_acc += static_cast<long double>(nM2s_sq);

        ++Ncells;
      }
    }
  }

  // L2
  const long double invN =
      (Ncells > 0) ? 1.0L / static_cast<long double>(Ncells) : 0.0L;

  out.vol_M0_L2 = std::sqrt(static_cast<double>(vol_M0_L2_acc * invN));
  out.vol_M1_L2 = std::sqrt(static_cast<double>(vol_M1_L2_acc * invN));
  out.vol_M2_L2 = std::sqrt(static_cast<double>(vol_M2_L2_acc * invN));
  out.surf_M0_L2 = std::sqrt(static_cast<double>(surf_M0_L2_acc * invN));
  out.surf_M1_L2 = std::sqrt(static_cast<double>(surf_M1_L2_acc * invN));
  out.surf_M2_L2 = std::sqrt(static_cast<double>(surf_M2_L2_acc * invN));

  return out;
}

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_RECONSTRUCTION_METRICS_TPP_
