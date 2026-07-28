#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_NEIGHBORHOOD_TPP_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_NEIGHBORHOOD_TPP_

namespace IRL {

template <class CellType>
PUNeighborhood<CellType>::PUNeighborhood(void) {}

template <class CellType>
void PUNeighborhood<CellType>::addMember(const Pt* a_centroid,
                                         const SeparatorVariant* a_separator,
                                         const double a_weight,
                                         const double a_scalar) {
  assert(a_centroid != nullptr);
  assert(a_separator != nullptr);
  centroids_m.push_back(*a_centroid);
  separators_m.push_back(*a_separator);
  weights_m.push_back(a_weight);
  a_scalar_m.push_back(a_scalar);
}

template <class CellType>
void PUNeighborhood<CellType>::emptyNeighborhood(void) {
  centroids_m.resize(0);
  separators_m.resize(0);
  weights_m.resize(0);
  a_scalar_m.resize(0);
}

template <class CellType>
void PUNeighborhood<CellType>::setMember(const UnsignedIndex_t a_index,
                                         const Pt* a_centroid,
                                         const SeparatorVariant* a_separator,
                                         const double a_weight,
                                         const double a_scalar) {
  assert(a_separator != nullptr);
  this->checkIndex(a_index);
  centroids_m[a_index] = *a_centroid;
  separators_m[a_index] = *a_separator;
  weights_m[a_index] = a_weight;
  a_scalar_m[a_index] = a_scalar;
}

template <class CellType>
const CellType& PUNeighborhood<CellType>::getCenterCell(void) const {
  return center_cell_m;
}
template <class CellType>
double PUNeighborhood<CellType>::getScalar(const Pt& a_pt) const {
  // Implementation for scalar interpolation
  double scalar_value = 0.0;
  // Use inverse distance weighting for interpolation
  for (int i = 0; i < a_scalar_m.size(); ++i) {
    double distance = std::sqrt(std::pow(a_pt[0] - centroids_m[i][0], 2) +
                                std::pow(a_pt[1] - centroids_m[i][1], 2) +
                                std::pow(a_pt[2] - centroids_m[i][2], 2));
    if (distance > 0) {
      scalar_value += a_scalar_m[i] / distance;
    }
  }
  return scalar_value;
}

template <class CellType>
void PUNeighborhood<CellType>::setCenterCell(const CellType* a_cell) {
  center_cell_m = *a_cell;
}

template <class CellType>
const SeparatorVariant& PUNeighborhood<CellType>::getSeparator(
    const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return separators_m[a_index];
}

template <class CellType>
const double PUNeighborhood<CellType>::getWeight(
    const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return weights_m[a_index];
}

template <class CellType>
const Pt& PUNeighborhood<CellType>::getCentroid(
    const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return centroids_m[a_index];
}

template <class CellType>
void PUNeighborhood<CellType>::resize(const UnsignedIndex_t a_size) {
  centroids_m.resize(a_size);
  separators_m.resize(a_size);
  weights_m.resize(a_size);
  a_scalar_m.resize(a_size);
}

template <class CellType>
void PUNeighborhood<CellType>::reserve(const UnsignedIndex_t a_size) {
  centroids_m.reserve(a_size);
  separators_m.reserve(a_size);
  weights_m.reserve(a_size);
  a_scalar_m.reserve(a_size);
}

template <class CellType>
UnsignedIndex_t PUNeighborhood<CellType>::size(void) const {
  return static_cast<UnsignedIndex_t>(centroids_m.size());
}

template <class CellType>
const std::vector<SeparatorVariant> PUNeighborhood<CellType>::getSeparators()
    const {
  return separators_m;
}

template <class CellType>
const std::vector<Pt> PUNeighborhood<CellType>::getCentroids() const {
  return centroids_m;
}

template <class CellType>
const std::vector<double> PUNeighborhood<CellType>::getWeights() const {
  return weights_m;
}

template <class CellType>
void PUNeighborhood<CellType>::checkIndex(UnsignedIndex_t a_index) const {
  assert(a_index < separators_m.size());
}

}  // End Namespace IRL

#endif