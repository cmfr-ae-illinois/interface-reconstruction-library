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
class PUNeighborhood {
 public:
  using cell_type = CellType;

  // Default Constructor Tested
  PUNeighborhood(void);

  /// \brief Construct CellGroupedMoments and add to end of collection. Tested
  void addMember(const Pt* a_centroid, const SeparatorVariant* a_separator,
                 const double a_weight = 1.0, const double a_scalar = 0.0);

  /// \brief Construct CellGroupedMoments and place into collection. Tested
  void setMember(const UnsignedIndex_t a_index, const Pt* a_centroid,
                 const SeparatorVariant* a_separator,
                 const double a_weight = 1.0, const double a_scalar = 0.0);

  /// \brief Reset neighborhood size to 0. Tested
  void emptyNeighborhood(void);

  /// \brief Set size of the neighborhood. Tested
  void resize(const UnsignedIndex_t a_size);
  void reserve(const UnsignedIndex_t a_size);

  // \brief Return an interpolated scalar
  double getScalar(const Pt& a_pt) const;

  // / \brief Return the center cell. Tested
  const CellType& getCenterCell(void) const;

  void setCenterOfStencil(const UnsignedIndex_t a_index);

  const UnsignedIndex_t& getCenterOfStencil(void) const;

  const SeparatorVariant& getSeparator(const UnsignedIndex_t a_index) const;

  const Pt& getCentroid(const UnsignedIndex_t a_index) const;

  const double getWeight(const UnsignedIndex_t a_index) const;

  const std::vector<SeparatorVariant> getSeparators() const;

  const std::vector<Pt> getCentroids() const;

  const std::vector<double> getWeights() const;

  void setCenterCell(const CellType* a_cell);

  /// \brief Get size of the vector. Tested
  UnsignedIndex_t size(void) const;

  /// \brief Default destructor.
  ~PUNeighborhood(void) = default;

 private:
  /// \brief Make sure index is not larger than current collection size.
  void checkIndex(UnsignedIndex_t a_index) const;

  std::vector<IRL::Pt> centroids_m;
  std::vector<IRL::SeparatorVariant> separators_m;
  std::vector<double> weights_m;
  std::vector<double> a_scalar_m;
  CellType center_cell_m;
  UnsignedIndex_t center_cell_index_m;
};
}  // End Namespace IRL

#include "irl/interface_reconstruction_methods/pu_neighborhood.tpp"

#endif