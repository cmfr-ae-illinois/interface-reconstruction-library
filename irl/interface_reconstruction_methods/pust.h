#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_H_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_H_

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
#include "irl/interface_reconstruction_methods/pu.h"
#include "irl/variant_reconstruction/separator_variant.h"
namespace IRL {
// This file contains all functions for surface tension forces using partition
// of unity methods.

template <class CellType>
class PUST : public PU<CellType> {
 public:
  using RationalBezierArc = RationalBezierArcBase<double>;
  // Constructor
  explicit PUST(const PUNeighborhood<CellType>& stencil_);

  PUST(const PUNeighborhood<CellType>& stencil_, const double kernel_size);
  // Default Constructor;
  PUST(void);
  // Solve Edge
  Normal solveEdge(const double STCoeff, const Pt& P0, const Pt& P1,
                   const double delta, const double Pressure,
                   const Normal& Marangoni);
  // Solve Face
  Normal solveFace(const double STCoeff, const Pt& P0, const Pt& P1,
                   const Pt& P2, const Pt& P3, const double delta,
                   const double Pressure, const Normal& Marangoni);
  // Debug
  void printSolver(void);

  // Ellipsoid Exact Values
  IRL::Pt projectOntoEllipsoid(const Pt& a_pt, const Normal& column1,
                               const Normal& column2, const Normal& column3,
                               const Pt& center);
  Normal getNormalEllipsoid(double x, double y, double z, const Normal& column1,
                            const Normal& column2, const Normal& column3,
                            const Pt& center);
  double getMeanCurvatureEllipsoid(Pt& x, const Normal& column1,
                                   const Normal& column2, const Normal& column3,
                                   const Pt& center);
  Normal solveFaceEllipsoid(const double STin, const Pt& P0, const Pt& P1,
                            const Pt& P2, const Pt& P3, const Normal& column1,
                            const Normal& column2, const Normal& column3,
                            const Pt& center, const double Pressure,
                            const Normal& Marangoni);
  double getPlaneCurvatureEllipsoid(Pt& x, const Normal& column1,
                                    const Normal& column2,
                                    const Normal& column3, const Pt& center,
                                    Normal& planeNormal);

  void evaluateEllipsoid(
      Pt& x, const Normal& column1, const Normal& column2,
      const Normal& column3, const Pt& center,
      std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal);

  void evaluateEllipsoid(Pt& x, const Normal& column1, const Normal& column2,
                         const Normal& column3, const Pt& center,
                         std::pair<double, Eigen::Vector3d>* retVal);

  void evaluateEllipsoid(Pt& x, const Normal& column1, const Normal& column2,
                         const Normal& column3, const Pt& center,
                         double* retVal);
};

}  // namespace IRL

#include "irl/interface_reconstruction_methods/pust.tpp"
#endif