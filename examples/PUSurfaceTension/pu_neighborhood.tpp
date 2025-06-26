#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_NEIGHBORHOOD_TPP_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_NEIGHBORHOOD_TPP_

namespace IRL {

    template <class CellType>
    PUSTNeighborhood<CellType>::PUSTNeighborhood(void)
        : center_cell_index_m(static_cast<IRL::UnsignedIndex_t>(-1)) {}
    
    template <class CellType>
    void PUSTNeighborhood<CellType>::addMember(const CellType* a_cell,
                                                const PlanarSeparator* a_plane) {
        assert(a_cell != nullptr);
        assert(a_plane != nullptr);
        collection_m.push_back(CGP(a_cell, a_plane));
    }

    template <class CellType>
    void PUSTNeighborhood<CellType>::emptyNeighborhood(void) {
        collection_m.resize(0);
        center_cell_index_m = static_cast<IRL::UnsignedIndex_t>(-1);
    }

    template <class CellType>
    void PUSTNeighborhood<CellType>::setMember(const UnsignedIndex_t a_index,
                                                const CellType* a_cell,
                                                const PlanarSeparator* a_plane) {
        assert(a_cell != nullptr);
        assert(a_plane != nullptr);
        this->checkIndex(a_index);
        collection_m[a_index] = CGP(a_cell,a_plane);
    }

    template <class CellType>
    void PUSTNeighborhood<CellType>::setCenterOfStencil(const UnsignedIndex_t a_index) {
        this->checkIndex(a_index);
        center_cell_index_m = a_index;
    }

    template <class CellType>
    const CellType& PUSTNeighborhood<CellType>::getCenterCell(void) const {
        this->checkCenterStencilSet();
        return this->getCell(center_cell_index_m);
    }

    template <class CellType>
    UnsignedIndex_t PUSTNeighborhood<CellType>::getCenterOfStencilIndex(void) const {
        this->checkCenterStencilSet();
        return center_cell_index_m;
    }

    template <class CellType>
    const PlanarSeparator& PUSTNeighborhood<CellType>::getCenterCellStoredMoments(void) const {
        this->checkCenterStencilSet();
        return this->getStoredMoments(center_cell_index_m);
    }

    template <class CellType>
    const typename PUSTNeighborhood<CellType>::CGP::cell_type&
    PUSTNeighborhood<CellType>::getCell(const UnsignedIndex_t a_index) const {
        this->checkIndex(a_index);
        return collection_m.getCell(a_index);
    }

    template <class CellType>
    const PlanarSeparator& PUSTNeighborhood<CellType>::getStoredMoments(const UnsignedIndex_t a_index) const {
        this->checkIndex(a_index);
        return collection_m.getStoredMoments(a_index);
    }

    template <class CellType>
    void PUSTNeighborhood<CellType>::resize(const UnsignedIndex_t a_size) {
        collection_m.resize(a_size);    
    }

    template <class CellType>
    UnsignedIndex_t PUSTNeighborhood<CellType>::size(void) const {
        return static_cast<UnsignedIndex_t>(collection_m.size());   
    }

    template <class CellType>
    typename PUSTNeighborhood<CellType>::iterator
    PUSTNeighborhood<CellType>::begin(void) noexcept {
        return collection_m.begin();
    }
        
    template <class CellType>
    typename PUSTNeighborhood<CellType>::const_iterator
    PUSTNeighborhood<CellType>::begin(void) const noexcept {
        return this->cbegin();
    }

    template <class CellType>
    typename PUSTNeighborhood<CellType>::const_iterator
    PUSTNeighborhood<CellType>::end(void) const noexcept {
        return this->cend();
    }

    template <class CellType>
    typename PUSTNeighborhood<CellType>::const_iterator
    PUSTNeighborhood<CellType>::cbegin(void) const noexcept {
        return collection_m.cbegin();
    }

    template <class CellType>
    typename PUSTNeighborhood<CellType>::iterator PUSTNeighborhood<CellType>::end(
        void) noexcept {
        return collection_m.end();
    }

    template <class CellType>
    typename PUSTNeighborhood<CellType>::const_iterator
    PUSTNeighborhood<CellType>::cend(void) const noexcept {
        return collection_m.cend();
    }
    
    template <class CellType>
    void PUSTNeighborhood<CellType>::checkIndex(UnsignedIndex_t a_index) const {
        assert(a_index < collection_m.size());
    }

    template <class CellType>
    void PUSTNeighborhood<CellType>::checkCenterStencilSet(void) const {
        assert(center_cell_index_m != static_cast<UnsignedIndex_t>(-1));
    }

} // End Namespace IRL

#endif