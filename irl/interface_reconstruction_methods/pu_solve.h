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
  const std::vector<double> weights;
  const double kernel_size;

 public:
  // Constructor
  PUImplicitSurface(const std::vector<Pt>& centroids_,
                    const std::vector<SeparatorVariant>& separators_,
                    const std::vector<double>& weights_,
                    const double& kernel_size_)
      : centroids(centroids_),
        separators(separators_),
        weights(weights_),
        kernel_size(kernel_size_) {}

  // Get Value
  inline void evaluate(const Pt& x, double* retVal);
  // Get Value, Grad
  inline void evaluate(const Pt& x, std::pair<double, Eigen::Vector3d>* retVal);
  // Get Value, Grad, Hessian
  inline void evaluate(
      const Pt& x,
      std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal);
  // Get Value
  inline void getTotalWeight(Pt& x, double* retVal);
  // Find intersection between the implicit curve and a provided line.
  inline std::vector<Pt> intersectEdge(const Pt& x0, const Pt& x1,
                                       const int& Npartitions,
                                       const double& tresh, bool& blocked);

  inline std::vector<Pt> intersectEdgeCylinder(const Pt& x0, const Pt& x1,
                                               double radius, Pt& center,
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

  inline void evaluateCylinder(
      Pt& x, double radius, Pt& center,
      std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal);

  inline void evaluateCylinder(Pt& x, double radius, Pt& center,
                               std::pair<double, Eigen::Vector3d>* retVal);

  inline void evaluateCylinder(Pt& x, double radius, Pt& center,
                               double* retVal);

  // Ellipsoid Values
  inline void evaluateEllipsoid(
      Pt& x, const Normal& column1, const Normal& column2,
      const Normal& column3, const Pt& center,
      std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal);

  inline void evaluateEllipsoid(Pt& x, const Normal& column1,
                                const Normal& column2, const Normal& column3,
                                const Pt& center,
                                std::pair<double, Eigen::Vector3d>* retVal);

  inline void evaluateEllipsoid(Pt& x, const Normal& column1,
                                const Normal& column2, const Normal& column3,
                                const Pt& center, double* retVal);
  inline Normal getNormalEllipsoid(Pt& x, const Normal& column1,
                                   const Normal& column2, const Normal& column3,
                                   const Pt& center);
  inline double getPlaneCurvatureEllipsoid(Pt& x, const Normal& column1,
                                           const Normal& column2,
                                           const Normal& column3,
                                           const Pt& center,
                                           Normal& planeNormal);
  inline double getPlaneCurvature(Pt& x, Normal& planeNormal);
  inline double getMeanCurvatureEllipsoid(Pt& x, const Normal& column1,
                                          const Normal& column2,
                                          const Normal& column3,
                                          const Pt& center);
  // Find the Tangent and Curvature at the point
  inline Normal getTangentCylinder(Pt& x, double radius, Pt& center);
  // Print the Implicit Surface Properties
  inline void printSurface();
  inline Normal getTangent(Pt& x);
  inline Normal getNormal(Pt& x);
  inline double getCurvature(Pt& x);

  // Project point onto implicit surface
  inline const Pt projectOntoPU(const Pt& a_pt);

  inline const Pt projectOntoEllipsoid(const Pt& a_pt, const Normal& column1,
                                       const Normal& column2,
                                       const Normal& column3, const Pt& center);
};

template <class CellType>
class PU {
 private:
  PUNeighborhood<CellType> stencil_m;
  // PUImplicitSurface surface_m;
  double intersection_threshold_m;

 public:
  using RationalBezierArc = RationalBezierArcBase<double>;
  // Constructor
  PU(PUNeighborhood<CellType> stencil_);
  // Default Constructor;
  PU(void);
  // Neighborhood Setter
  void setNeighborhood(PUNeighborhood<CellType> stencil_);
  // Threshold Setter
  void setThreshold(double tresh);
  // Takes Neighborhood and Returns the Implicit Surface
  PUImplicitSurface neighborhoodToImplicitSurface(double delta);
  // Edge Solve Method - Returns the surface tension force vector
  // for edge
  Normal solveEdge(const double STIn, const Pt& P0, const Pt& P1,
                   const double delta, const double Pressure,
                   const Normal& Marangoni);

  Normal solveFace(const double STin, const Pt& P0, const Pt& P1, const Pt& P2,
                   const Pt& P3, const double delta, const double Pressure,
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
  Normal getNormal(double x, double y, double z, double delta);
  // Get Total weight
  double getWeight(double x, double y, double z, double delta);
  double getWeight(Pt& in, double delta);
  double getWeightCylinder(double x, double y, double z, double radius,
                           Pt center);
  // Curvature
  double getCurvature(double x, double y, double z, double delta);
  double getCurvature(Pt& in, double delta);
  double getMeanCurvatureEllipsoid(double x, double y, double z,
                                   const Normal& column1, const Normal& column2,
                                   const Normal& column3, const Pt& center);
  double getMeanCurvatureEllipsoid(Pt& in, const Normal& column1,
                                   const Normal& column2, const Normal& column3,
                                   const Pt& center);

  // Exact Values for Ellipsoid Stress
  Normal solveFaceEllipsoid(const double STin, const Pt& P0, const Pt& P1,
                            const Pt& P2, const Pt& P3, const Normal& column1,
                            const Normal& column2, const Normal& column3,
                            const Pt& center, const double Pressure,
                            const Normal& Marangoni);
  Normal getNormalEllipsoid(double x, double y, double z, const Normal& column1,
                            const Normal& column2, const Normal& column3,
                            const Pt& center);
  // Print
  void printSolver();
  /// \brief Solve the system for the reconstruction
  Paraboloid solve(const PUNeighborhood<CellType>* a_neighborhood_pointer,
                   const Pt& a_centroid, const double a_delta = -1.0);
  /// \brief Project onto the PU surface
  IRL::Pt projectOntoPU(const Pt& a_pt, const double a_delta);
  /// \brief Project onto the PU surface
  IRL::Pt projectOntoEllipsoid(const Pt& a_pt, const Normal& column1,
                               const Normal& column2, const Normal& column3,
                               const Pt& center);
};

}  // End Namespace IRL

#include "irl/interface_reconstruction_methods/pu_solve.tpp"

#endif