#ifndef EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_BINARY_TPP_
#define EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_BINARY_TPP_

#include "examples/implicit_surface_reconstruction/binary.h"

// writing moments to binary file
template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void writeMomentsToBinary(
    const Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                         IRL::GeneralSurfaceMoments3D<SM_ORDER>>>& moments,
    const std::string& filename) {
  constexpr int NV = (VM_ORDER + 1) * (VM_ORDER + 2) * (VM_ORDER + 3) / 6;
  constexpr int NS = (SM_ORDER + 1) * (SM_ORDER + 2) * (SM_ORDER + 3) / 6;

  const BasicMesh& mesh = moments.getMesh();

  struct Record {
    int i, j, k;
    double vol[NV];
    double surf[NS];
  };

  std::ofstream ofs(filename, std::ios::binary);
  if (!ofs)
    throw std::runtime_error("Cannot open file for writing: " + filename);

  // writing grid bounds
  int imin = mesh.imin(), imax = mesh.imax();
  int jmin = mesh.jmin(), jmax = mesh.jmax();
  int kmin = mesh.kmin(), kmax = mesh.kmax();
  ofs.write(reinterpret_cast<char*>(&imin), sizeof(int));
  ofs.write(reinterpret_cast<char*>(&imax), sizeof(int));
  ofs.write(reinterpret_cast<char*>(&jmin), sizeof(int));
  ofs.write(reinterpret_cast<char*>(&jmax), sizeof(int));
  ofs.write(reinterpret_cast<char*>(&kmin), sizeof(int));
  ofs.write(reinterpret_cast<char*>(&kmax), sizeof(int));

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const auto& pair = moments(i, j, k);
        Record rec{};
        rec.i = i;
        rec.j = j;
        rec.k = k;
        for (int t = 0; t < NV; ++t) rec.vol[t] = pair.first[t];
        for (int t = 0; t < NS; ++t) rec.surf[t] = pair.second[t];
        ofs.write(reinterpret_cast<const char*>(&rec), sizeof(rec));
      }
    }
  }

  ofs.close();
  std::cout << "✅ Moments written to binary file: " << filename << std::endl;
}

// reading moments from binary file
template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void readMomentsFromBinary(
    const std::string& filename,
    Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                   IRL::GeneralSurfaceMoments3D<SM_ORDER>>>* moments) {
  constexpr int NV = (VM_ORDER + 1) * (VM_ORDER + 2) * (VM_ORDER + 3) / 6;
  constexpr int NS = (SM_ORDER + 1) * (SM_ORDER + 2) * (SM_ORDER + 3) / 6;

  struct Record {
    int i, j, k;
    double vol[NV];
    double surf[NS];
  };

  std::ifstream ifs(filename, std::ios::binary);
  if (!ifs)
    throw std::runtime_error("Cannot open file for reading: " + filename);

  int imin, imax, jmin, jmax, kmin, kmax;
  ifs.read(reinterpret_cast<char*>(&imin), sizeof(int));
  ifs.read(reinterpret_cast<char*>(&imax), sizeof(int));
  ifs.read(reinterpret_cast<char*>(&jmin), sizeof(int));
  ifs.read(reinterpret_cast<char*>(&jmax), sizeof(int));
  ifs.read(reinterpret_cast<char*>(&kmin), sizeof(int));
  ifs.read(reinterpret_cast<char*>(&kmax), sizeof(int));

  const BasicMesh& mesh = moments->getMesh();

  Record rec{};
  while (ifs.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
    IRL::GeneralMoments3D<VM_ORDER> vol{};
    IRL::GeneralSurfaceMoments3D<SM_ORDER> surf{};
    for (int t = 0; t < NV; ++t) vol[t] = rec.vol[t];
    for (int t = 0; t < NS; ++t) surf[t] = rec.surf[t];
    (*moments)(rec.i, rec.j, rec.k).first = vol;
    (*moments)(rec.i, rec.j, rec.k).second = surf;
  }

  ifs.close();
  std::cout << "✅ Moments loaded from binary file: " << filename << std::endl;
}

inline void kahan_add(double& sum, double& c, const double x) {
  const double y = x - c;
  const double t = sum + y;
  c = (t - sum) - y;
  sum = t;
}

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void coarsenMomentsFromBinary(
    const std::string& fine_filename, int factor,
    Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                   IRL::GeneralSurfaceMoments3D<SM_ORDER>>>* coarse_moments) {
  if (factor <= 0) throw std::invalid_argument("factor must be > 0");
  if (!coarse_moments) throw std::invalid_argument("coarse_moments is null");

  constexpr int NV = (VM_ORDER + 1) * (VM_ORDER + 2) * (VM_ORDER + 3) / 6;
  constexpr int NS = (SM_ORDER + 1) * (SM_ORDER + 2) * (SM_ORDER + 3) / 6;

  struct Record {
    int i, j, k;
    double vol[NV];
    double surf[NS];
  };

  std::ifstream ifs(fine_filename, std::ios::binary);
  if (!ifs) throw std::runtime_error("Cannot open fine file: " + fine_filename);

  // reading fine mesh index bounds from file header
  int fi_min, fi_max, fj_min, fj_max, fk_min, fk_max;
  ifs.read(reinterpret_cast<char*>(&fi_min), sizeof(int));
  ifs.read(reinterpret_cast<char*>(&fi_max), sizeof(int));
  ifs.read(reinterpret_cast<char*>(&fj_min), sizeof(int));
  ifs.read(reinterpret_cast<char*>(&fj_max), sizeof(int));
  ifs.read(reinterpret_cast<char*>(&fk_min), sizeof(int));
  ifs.read(reinterpret_cast<char*>(&fk_max), sizeof(int));

  const int nx_f = fi_max - fi_min + 1;
  const int ny_f = fj_max - fj_min + 1;
  const int nz_f = fk_max - fk_min + 1;

  // coarse mesh info
  const BasicMesh& cmesh = coarse_moments->getMesh();
  const int ci_min = cmesh.imin(), ci_max = cmesh.imax();
  const int cj_min = cmesh.jmin(), cj_max = cmesh.jmax();
  const int ck_min = cmesh.kmin(), ck_max = cmesh.kmax();

  const int nx_c = ci_max - ci_min + 1;
  const int ny_c = cj_max - cj_min + 1;
  const int nz_c = ck_max - ck_min + 1;

  // checking factor value
  if (nx_f != factor * nx_c || ny_f != factor * ny_c || nz_f != factor * nz_c) {
    throw std::runtime_error(
        "Grid size mismatch: fine dims must equal factor * coarse dims in each "
        "direction.");
  }

  // coarse sums
  for (int I = ci_min; I <= ci_max; ++I)
    for (int J = cj_min; J <= cj_max; ++J)
      for (int K = ck_min; K <= ck_max; ++K) {
        auto& p = (*coarse_moments)(I, J, K);
        for (int t = 0; t < NV; ++t) p.first[t] = 0.0;
        for (int t = 0; t < NS; ++t) p.second[t] = 0.0;
      }

  // compensation fields
  Data<IRL::GeneralMoments3D<VM_ORDER>> cVol(&cmesh);
  Data<IRL::GeneralSurfaceMoments3D<SM_ORDER>> cSurf(&cmesh);

  for (int I = ci_min; I <= ci_max; ++I)
    for (int J = cj_min; J <= cj_max; ++J)
      for (int K = ck_min; K <= ck_max; ++K) {
        auto& cv = cVol(I, J, K);
        auto& cs = cSurf(I, J, K);
        for (int t = 0; t < NV; ++t) cv[t] = 0.0;
        for (int t = 0; t < NS; ++t) cs[t] = 0.0;
      }

  // accumulate moments with kahan summation
  Record rec{};
  while (ifs.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
    // Map fine indices to coarse indices
    const int di = rec.i - fi_min;
    const int dj = rec.j - fj_min;
    const int dk = rec.k - fk_min;

    const int bi = di / factor;
    const int bj = dj / factor;
    const int bk = dk / factor;

    const int I = ci_min + bi;
    const int J = cj_min + bj;
    const int K = ck_min + bk;

    // check
    if (I < ci_min || I > ci_max || J < cj_min || J > cj_max || K < ck_min ||
        K > ck_max) {
      throw std::runtime_error(
          "Computed coarse index out of bounds while coarsening.");
    }

    auto& sumPair = (*coarse_moments)(I, J, K);
    auto& cV = cVol(I, J, K);
    auto& cS = cSurf(I, J, K);

    // Compensated accumulation per component
    for (int t = 0; t < NV; ++t) kahan_add(sumPair.first[t], cV[t], rec.vol[t]);
    for (int t = 0; t < NS; ++t)
      kahan_add(sumPair.second[t], cS[t], rec.surf[t]);
  }
}

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_BINARY_TPP_
