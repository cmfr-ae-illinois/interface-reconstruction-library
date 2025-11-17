#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_H_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_H_

#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "irl/interface_reconstruction_methods/pu_neighborhood.h"

#include "irl/moments/cell_collection.h"
#include "irl/moments/cell_grouped_moments.h"

#include "irl/geometry/general/normal.h"
#include "irl/variant_reconstruction/separator_variant.h"
namespace IRL {
/// This file contains all the functions
/// needed to calculate the surface tension
/// Forces.

// template<class SeparatorType>
class PUImplicitSurface {
 private:
  const std::vector<Pt> centroids;
  const std::vector<SeparatorVariant> separators;
  const double kernel_size;

 public:
  // Constructor
  PUImplicitSurface(const std::vector<Pt>& centroids_,
                    const std::vector<SeparatorVariant>& separators_,
                    const double& kernel_size_)
      : centroids(centroids_),
        separators(separators_),
        kernel_size(kernel_size_) {}

  // Get Value
  inline void evaluate(const Pt& x, double* retVal);
  // Get Value, Grad
  inline void evaluate(const Pt& x, std::pair<double, Eigen::Vector3d>* retVal);
  // Get Value, Grad, Hessian
  void evaluate(Pt& x,
                std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal);
  // Get Value
  void getTotalWeight(Pt& x, double* retVal);
  // Find intersection between the implicit curve and a provided line.
  std::vector<Pt> intersectEdge(const Pt& x0, const Pt& x1,
                                const int& Npartitions, const double& tresh);

  std::vector<Pt> intersectEdgeCylinder(const Pt& x0, const Pt& x1,
                                        double radius, Pt& center,
                                        const int& Npartitions);
  inline void evaluate(
      const Pt& x,
      std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal);

  // Find intersection between the implicit curve and a provided line.
  inline std::vector<Pt> intersectEdge(const Pt& x0, const Pt& x1,
                                       const int& Npartitions);
  // Signed Distance of Separator
  static inline void implicitSeparator(const Pt& a_pt, const Pt& a_centroid,
                                       const SeparatorVariant* a_sepPtr,
                                       double* retVal);
  // Signed Distance and Gradient of Separator
  static inline void implicitSeparator(
      const Pt& a_pt, const Pt& a_centroid, const SeparatorVariant* a_sepPtr,
      std::pair<double, Eigen::Vector3d>* retVal);
  // Signed Distance, Gradient, and Hessian of Separator
  static inline void implicitSeparator(
      const Pt& a_pt, const Pt& a_centroid, const SeparatorVariant* a_sepPtr,
      std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal);

  void evaluateCylinder(
      Pt& x, double radius, Pt& center,
      std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal);

  void evaluateCylinder(Pt& x, double radius, Pt& center,
                        std::pair<double, Eigen::Vector3d>* retVal);

  void evaluateCylinder(Pt& x, double radius, Pt& center, double* retVal);

  //
  // Find the Tangent and Curvature at the point
  inline Normal getTangent(Pt& x);
  Normal getTangentCylinder(Pt& x, double radius, Pt& center);
  inline double getCurvature(Pt& x);
  // Print the Implicit Surface Properties
  void printSurface();
  inline Normal getTangent(Pt& x);
  inline double getCurvature(Pt& x);

  // Project point onto implicit surface
  inline const Pt projectOntoPU(const Pt& a_pt);
};

template <class CellType>
class PUST {
 private:
  PUSTNeighborhood<CellType> stencil_m;
  // PUImplicitSurface surface_m;
  double intersection_threshold_m;

 public:
  // Constructor
  PUST(PUSTNeighborhood<CellType> stencil_);
  // Default Constructor;
  PUST(void);
  // Neighborhood Setter
  void setNeighborhood(PUSTNeighborhood<CellType> stencil_);
  // Threshold Setter
  void setThreshold(double tresh);
  // Takes Neighborhood and Returns the Implicit Surface
  PUImplicitSurface neighborhoodToImplicitSurface(double delta);
  // Edge Solve Method - Returns the surface tension force vector
  // for edge
  Normal solveEdge(const double STIn, const Pt& P0, const Pt& P1,
                   const double delta, const double Pressure,
                   const Normal& Marangoni);

  Normal solveEdgeCylinder(double STCoeff, Pt& P0, Pt& P1, double radius,
                           Pt& center, double delta);
  // Get value function for neighborhood
  double getValue(double x, double y, double z, double delta);
  double getValueCylinder(double x, double y, double z, double radius,
                          Pt center);
  // Get Tangent for Neighborhood
  Normal getTangent(double x, double y, double z, double delta);
  Normal getTangentCylinder(double x, double y, double z, double radius,
                            Pt center);
  // Get Total weight
  double getWeight(double x, double y, double z, double delta);
  double getWeightCylinder(double x, double y, double z, double radius,
                           Pt center);
  // Print
  void printSolver();
  Normal solveEdge(double STCoeff, Pt& P0, Pt& P1);
  /// \brief Solve the system for the reconstruction
  Paraboloid solve(const PUSTNeighborhood<CellType>* a_neighborhood_pointer,
                   const Pt& a_centroid, const double a_delta = -1.0);
};

}  // End Namespace IRL

#include "irl/interface_reconstruction_methods/pu_solve.tpp"

#endif