#ifndef EXAMPLES_AMREX_ADVECTOR_CONNECTED_COMPONENTS_LABELING_H_
#define EXAMPLES_AMREX_ADVECTOR_CONNECTED_COMPONENTS_LABELING_H_

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_iMultiFab.H>

#ifdef AMREX_USE_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace CCL {

class UnionFind {
 public:
  explicit UnionFind(int n) : parent_(n), size_(n, 1) {
    std::iota(parent_.begin(), parent_.end(),
              0);  // initially each element is its parent and size of each
                   // component is 1
  }

  // determines the root of the set containing element x.
  int find(int x) {
    while (parent_[x] != x) {
      parent_[x] = parent_[parent_[x]];
      x = parent_[x];
    }
    return x;
  }

  bool unite(int a, int b) {
    // find the roots
    a = find(a);
    b = find(b);

    if (a == b) return false;  // skip if they belong to the same component

    if (size_[a] < size_[b])
      std::swap(a, b);  // smaller component is attached to the larger one

    parent_[b] = a;
    size_[a] += size_[b];

    return true;
  }

 private:
  std::vector<int> parent_;
  std::vector<int> size_;
};

inline int CountConnectedComponents(
    const amrex::MultiFab& vf, const amrex::Geometry& geom,
    amrex::Real threshold = 1.0e-10, int vf_comp = 0,
    amrex::MultiFab* component_labels = nullptr) {
  using namespace amrex;

  if (vf_comp < 0 || vf_comp >= vf.nComp()) {
    amrex::Abort("CCL: invalid volume-fraction component");
  }

  const int rank = ParallelDescriptor::MyProc();
  const int nprocs = ParallelDescriptor::NProcs();

  // component labels
  iMultiFab labels(vf.boxArray(), vf.DistributionMap(), 1, 1);
  labels.setVal(0);

  // Number of distinct components found locally
  // this is the number of components found in boxes owned by this MPI rank
  int nlocal = 0;

  // labeling components per box
  for (MFIter mfi(vf); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.validbox();

    const auto arr = vf.const_array(mfi);
    const auto lab = labels.array(mfi);

    const int ilo = bx.smallEnd(0);
    const int jlo = bx.smallEnd(1);
    const int klo = bx.smallEnd(2);

    const int nx = bx.length(0);
    const int ny = bx.length(1);
    const int nz = bx.length(2);

    const int int_max = std::numeric_limits<int>::max();

    if (nx <= 0 || ny <= 0 || nz <= 0 || nx > int_max / ny ||
        nx * ny > int_max / nz) {
      amrex::Abort("CCL: AMReX box is too large");
    }

    const int ncells = nx * ny * nz;

    UnionFind uf(ncells);

    // linear indexing
    auto index = [=](int i, int j, int k) {
      return ((k - klo) * ny + (j - jlo)) * nx + (i - ilo);
    };

    // merging face connected liquid cells
    for (int k = klo; k <= bx.bigEnd(2); ++k) {
      for (int j = jlo; j <= bx.bigEnd(1); ++j) {
        for (int i = ilo; i <= bx.bigEnd(0); ++i) {
          if (arr(i, j, k, vf_comp) <= threshold) continue;

          const int current = index(i, j, k);

          // Negative x neighbor.
          if (i > ilo && arr(i - 1, j, k, vf_comp) > threshold) {
            uf.unite(current, index(i - 1, j, k));
          }

          // Negative y neighbor.
          if (j > jlo && arr(i, j - 1, k, vf_comp) > threshold) {
            uf.unite(current, index(i, j - 1, k));
          }

          // Negative z neighbor.
          if (k > klo && arr(i, j, k - 1, vf_comp) > threshold) {
            uf.unite(current, index(i, j, k - 1));
          }
        }
      }
    }

    // assigning one ID to each local component
    std::vector<int> root_to_id(ncells, 0);

    for (int k = klo; k <= bx.bigEnd(2); ++k) {
      for (int j = jlo; j <= bx.bigEnd(1); ++j) {
        for (int i = ilo; i <= bx.bigEnd(0); ++i) {
          if (arr(i, j, k, vf_comp) <= threshold) continue;

          const int root = uf.find(index(i, j, k));

          if (root_to_id[root] == 0) {
            if (nlocal == std::numeric_limits<int>::max()) {
              amrex::Abort("CCL: too many local components");
            }

            root_to_id[root] = ++nlocal;
          }
          // all cells with the same root get the same label
          lab(i, j, k) = root_to_id[root];
        }
      }
    }
  }

  // making component IDs globally unique
  std::vector<int> counts(nprocs, 0);

  // gather the number of local components from all MPI ranks and store them in
  // the vector
#ifdef AMREX_USE_MPI
  MPI_Allgather(&nlocal, 1, MPI_INT, counts.data(), 1, MPI_INT,
                ParallelDescriptor::Communicator());
#else
  counts[0] = nlocal;
#endif

  int total = 0;   // total number of components across all MPI ranks
  int offset = 0;  // offset of this rank's components in the global ID space

  for (int p = 0; p < nprocs; ++p) {
    if (counts[p] < 0 ||
        counts[p] > std::numeric_limits<int>::max() - 1 - total) {
      amrex::Abort("CCL: total provisional labels exceed INT_MAX");
    }

    if (p < rank) offset += counts[p];

    total += counts[p];
  }

  // The label MultiFab uses 32-bit integers.
  if (total >= std::numeric_limits<int>::max()) {
    amrex::Abort("CCL: total provisional labels exceed INT_MAX");
  }

  // Shift rank-local component IDs into global ID space.
  for (MFIter mfi(labels); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.validbox();
    const auto lab = labels.array(mfi);

    for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
      for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
          if (lab(i, j, k) > 0) lab(i, j, k) += offset;
        }
      }
    }
  }

  // exchange ghost cell labels across boxes on same rank and different ranks
  labels.FillBoundary(geom.periodicity());

  // Identify component equivalences across boxes.

  std::vector<std::pair<int, int>> edges;

  for (MFIter mfi(labels); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.validbox();
    const auto lab = labels.const_array(mfi);

    for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
      for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
          const int a = lab(i, j, k);

          if (a <= 0) continue;

          // Positive x face.
          {
            const int b = lab(i + 1, j, k);
            if (b > 0 && a != b) {
              edges.emplace_back(std::min(a, b), std::max(a, b));
            }
          }

          // Positive y face.
          {
            const int b = lab(i, j + 1, k);
            if (b > 0 && a != b) {
              edges.emplace_back(std::min(a, b), std::max(a, b));
            }
          }

          // Positive z face.
          {
            const int b = lab(i, j, k + 1);
            if (b > 0 && a != b) {
              edges.emplace_back(std::min(a, b), std::max(a, b));
            }
          }
        }
      }
    }
  }

  // Remove duplicate equivalences locally
  std::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

  // Gather component equivalences onto MPI root.

  // Each pair (a,b) represents an edge in the component graph.

  // edge pairs to flattended array
  std::vector<int> sendbuf;
  sendbuf.reserve(2 * edges.size());
  for (const auto& edge : edges) {
    sendbuf.push_back(edge.first);
    sendbuf.push_back(edge.second);
  }

  if (sendbuf.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    amrex::Abort("CCL: too many local equivalences");
  }

  const int sendcount = static_cast<int>(sendbuf.size());

  std::vector<int> recvcounts(
      nprocs, 0);  // number of elements to receive from each rank
  std::vector<int> displs(
      nprocs, 0);  // displacement of each rank's data in the receive buffer
  std::vector<int> recvbuf;  // buffer to store received data on the root rank
                             // from all other ranks

#ifdef AMREX_USE_MPI

  // each rank sends 2 times the number of local edges
  MPI_Gather(&sendcount, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0,
             ParallelDescriptor::Communicator());

  if (rank == 0) {
    int ntotal = 0;

    for (int p = 0; p < nprocs; ++p) {
      if (recvcounts[p] < 0 ||
          recvcounts[p] > std::numeric_limits<int>::max() - ntotal) {
        amrex::Abort("CCL: equivalence buffer too large");
      }

      displs[p] = ntotal;
      ntotal += recvcounts[p];
    }

    if (ntotal > std::numeric_limits<int>::max()) {
      amrex::Abort("CCL: equivalence buffer too large");
    }

    recvbuf.resize(static_cast<std::size_t>(ntotal));
  }

  // gather the flattened edge arrays from all ranks to the root rank
  MPI_Gatherv(sendbuf.data(), sendcount, MPI_INT,
              rank == 0 ? recvbuf.data() : nullptr,
              rank == 0 ? recvcounts.data() : nullptr,
              rank == 0 ? displs.data() : nullptr, MPI_INT, 0,
              ParallelDescriptor::Communicator());

#else

  recvbuf = sendbuf;

#endif

  // Global Union-Find and final component count
  int ncomponents = 0;
  std::vector<int> final_labels(total + 1, 0);

  if (rank == 0) {
    // IDs are 1,...,total.
    UnionFind global_uf(total + 1);

    ncomponents = total;

    for (std::size_t e = 0; e + 1 < recvbuf.size(); e += 2) {
      const int a = recvbuf[e];
      const int b = recvbuf[e + 1];
      // A successful union removes one connected component
      if (global_uf.unite(a, b)) {
        --ncomponents;
      }
    }

    std::vector<int> root_to_component(total + 1, 0);
    int next_component = 0;
    for (int id = 1; id <= total; ++id) {
      const int root = global_uf.find(id);
      if (root_to_component[root] == 0) {
        root_to_component[root] = ++next_component;
      }
      final_labels[id] = root_to_component[root];
    }
  }

#ifdef AMREX_USE_MPI
  MPI_Bcast(&ncomponents, 1, MPI_INT, 0, ParallelDescriptor::Communicator());
  MPI_Bcast(final_labels.data(), total + 1, MPI_INT, 0,
            ParallelDescriptor::Communicator());
#else
  for (int id = 1; id <= total; ++id) {
    final_labels[id] = id;
  }
#endif

  for (MFIter mfi(labels); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.validbox();
    const auto lab = labels.array(mfi);
    for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
      for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
          const int provisional_label = lab(i, j, k);
          if (provisional_label > 0) {
            lab(i, j, k) = final_labels[provisional_label];
          }
        }
      }
    }
  }

  if (component_labels != nullptr) {
    component_labels->define(vf.boxArray(), vf.DistributionMap(), 1, 0);
    component_labels->setVal(0.0);
    for (MFIter mfi(*component_labels); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.validbox();
      const auto lab = labels.const_array(mfi);
      const auto output = component_labels->array(mfi);
      for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
          for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
            output(i, j, k, 0) = static_cast<amrex::Real>(lab(i, j, k));
          }
        }
      }
    }
  }

  return ncomponents;
}

}  // namespace CCL

#endif  // EXAMPLES_AMREX_ADVECTOR_CONNECTED_COMPONENTS_LABELING_H_