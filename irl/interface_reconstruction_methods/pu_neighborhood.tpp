#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_NEIGHBORHOOD_TPP_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_NEIGHBORHOOD_TPP_

namespace IRL {

template <class CellType>
PUSTNeighborhood<CellType>::PUSTNeighborhood(void) {}

template <class CellType>
void PUSTNeighborhood<CellType>::addMember(
    const Pt* a_centroid, const SeparatorVariant* a_separator) {
  assert(a_centroid != nullptr);
  assert(a_separator != nullptr);
  centroids_m.push_back(*a_centroid);
  separators_m.push_back(*a_separator);
}

template <class CellType>
void PUSTNeighborhood<CellType>::emptyNeighborhood(void) {
  centroids_m.resize(0);
  separators_m.resize(0);
}

template <class CellType>
void PUSTNeighborhood<CellType>::setMember(
    const UnsignedIndex_t a_index, const Pt* a_centroid,
    const SeparatorVariant* a_separator) {
  assert(a_cell != nullptr);
  assert(a_plane != nullptr);
  this->checkIndex(a_index);
  centroids_m[a_index] = *a_centroid;
  separators_m[a_index] = *a_separator;
}

// template <class CellType>
// void PUSTNeighborhood<CellType>::setCenterOfStencil(const UnsignedIndex_t
// a_index) {
//     this->checkIndex(a_index);
//     center_cell_index_m = a_index;
// }

template <class CellType>
const CellType& PUSTNeighborhood<CellType>::getCenterCell(void) const {
  return center_cell_m;
}

template <class CellType>
void PUSTNeighborhood<CellType>::setCenterCell(const CellType* a_cell) {
  center_cell_m = *a_cell;
}

// template <class CellType>
// UnsignedIndex_t PUSTNeighborhood<CellType>::getCenterOfStencilIndex(void)
// const {
//     this->checkCenterStencilSet();
//     return center_cell_index_m;
// }

// template <class CellType>
// const SeparatorType&
// PUSTNeighborhood<CellType>::getCenterCellStoredMoments(void) const {
//     this->checkCenterStencilSet();
//     return this->getStoredMoments(center_cell_index_m);
// }

// template <class CellType>
// const typename PUSTNeighborhood<CellType>::CGP::cell_type&
// PUSTNeighborhood<CellType>::getCell(const UnsignedIndex_t a_index) const {
//     this->checkIndex(a_index);
//     return collection_m.getCell(a_index);
// }

// template <class CellType>
// const SeparatorVariant& PUSTNeighborhood<CellType>::getStoredMoments(const
// UnsignedIndex_t a_index) const {
//     this->checkIndex(a_index);
//     return collection_m.getStoredMoments(a_index);
// }

template <class CellType>
const SeparatorVariant& PUSTNeighborhood<CellType>::getSeparator(
    const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return separators_m[a_index];
}

template <class CellType>
const Pt& PUSTNeighborhood<CellType>::getCentroid(
    const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return centroids_m[a_index];
}

template <class CellType>
void PUSTNeighborhood<CellType>::resize(const UnsignedIndex_t a_size) {
  centroids_m.resize(a_size);
  separators_m.resize(a_size);
}

template <class CellType>
void PUSTNeighborhood<CellType>::reserve(const UnsignedIndex_t a_size) {
  centroids_m.reserve(a_size);
  separators_m.reserve(a_size);
}

template <class CellType>
UnsignedIndex_t PUSTNeighborhood<CellType>::size(void) const {
  return static_cast<UnsignedIndex_t>(centroids_m.size());
}

template <class CellType>
const std::vector<SeparatorVariant> PUSTNeighborhood<CellType>::getSeparators()
    const {
  return separators_m;
}

template <class CellType>
const std::vector<Pt> PUSTNeighborhood<CellType>::getCentroids() const {
  return centroids_m;
}
// template <class CellType>
// typename PUSTNeighborhood<CellType>::iterator
// PUSTNeighborhood<CellType>::begin(void) noexcept {
//     return collection_m.begin();
// }

// template <class CellType>
// typename PUSTNeighborhood<CellType>::const_iterator
// PUSTNeighborhood<CellType>::begin(void) const noexcept {
//     return this->cbegin();
// }

// template <class CellType>
// typename PUSTNeighborhood<CellType>::const_iterator
// PUSTNeighborhood<CellType>::end(void) const noexcept {
//     return this->cend();
// }

// template <class CellType>
// typename PUSTNeighborhood<CellType>::const_iterator
// PUSTNeighborhood<CellType>::cbegin(void) const noexcept {
//     return collection_m.cbegin();
// }

// template <class CellType>
// typename PUSTNeighborhood<CellType>::iterator
// PUSTNeighborhood<CellType>::end(
//     void) noexcept {
//     return collection_m.end();
// }

// template <class CellType>
// typename PUSTNeighborhood<CellType>::const_iterator
// PUSTNeighborhood<CellType>::cend(void) const noexcept {
//     return collection_m.cend();
// }

template <class CellType>
void PUSTNeighborhood<CellType>::checkIndex(UnsignedIndex_t a_index) const {
  assert(a_index < collection_m.size());
}

// template <class CellType>
// void PUSTNeighborhood<CellType>::checkCenterStencilSet(void) const {
//     assert(center_cell_index_m != static_cast<UnsignedIndex_t>(-1));
// }

}  // End Namespace IRL

#endif