#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_NEIGHBORHOOD_H_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_NEIGHBORHOOD_H_

#include "irl/moments/cell_collection.h"
#include "irl/moments/cell_grouped_moments.h"

#include "irl/variant_reconstruction/separator_variant.h"

namespace IRL {
    /// \brief Neighborhood storaged used in the surface tension 
    /// using the partition of unity method. This stores the 
    /// CellGroupedMoments of the cell and planar Separator.
    template <class CellType>
    class PUSTNeighborhood {
        using CGP = CellGroupedMoments<CellType, PlanarSeparator>;
        using iterator = typename CellCollection<CGP>::iterator;
        using const_iterator = typename CellCollection<CGP>::const_iterator;

        public:
            using cell_type = CellType;

            // Default Constructor Tested
            PUSTNeighborhood(void);

            /// \brief Construct CellGroupedMoments and add to end of collection. Tested
            void addMember(const CellType* a_cell, const PlanarSeparator* a_plane);

            /// \brief Construct CellGroupedMoments and place into collection. Tested
            void setMember(const UnsignedIndex_t a_index, const CellType* a_cell,
                            const PlanarSeparator* a_plane);

            /// \brief Reset neighborhood size to 0. Tested
            void emptyNeighborhood(void);

            /// \brief Set size of the neighborhood. Tested
            void resize(const UnsignedIndex_t a_size);

            /// \brief Set the index for the center cell in the collection. Tested
            void setCenterOfStencil(const UnsignedIndex_t a_index);

            /// \brief Return the index for the center stencil Tested
            UnsignedIndex_t getCenterOfStencilIndex(void) const;

            /// \brief Return the center cell. Tested
            const CellType& getCenterCell(void) const;

            /// \brief Return the center cell moments Tested
            const PlanarSeparator& getCenterCellStoredMoments(void) const;

            /// \brief Return the cell stored at the index Tested
            const typename CGP::cell_type& getCell(const UnsignedIndex_t a_index) const;

            /// \brief Return moments stored at the index Tested
            const PlanarSeparator& getStoredMoments(const UnsignedIndex_t a_index) const;

            /// \brief Get size of the vector. Tested
            UnsignedIndex_t size(void) const;

            iterator begin(void) noexcept;
            const_iterator begin(void) const noexcept;
            const_iterator end(void) const noexcept;
            const_iterator cbegin(void) const noexcept;
            iterator end(void) noexcept;
            const_iterator cend(void) const noexcept;

            /// \brief Default destructor.
            ~PUSTNeighborhood(void) = default;

        private:
            /// \brief Make sure index is not larger than current collection size.
            void checkIndex(UnsignedIndex_t a_index) const;
            void checkCenterStencilSet(void) const;

            /// \brief Collection of cells and planes.
            CellCollection<CGP> collection_m;
            
            /// \brief Center stencil cell index in the list of added cells.
            UnsignedIndex_t center_cell_index_m;
    };
} // End Namespace IRL

#include "examples/PUSurfaceTension/pu_neighborhood.tpp"

#endif