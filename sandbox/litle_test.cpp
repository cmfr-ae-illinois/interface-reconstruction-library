//#define VALDEBUG2

#include "sandbox/float128_stream.h"

#include "irl/geometry/general/new_pt_calculation_functors.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/pt_with_data.h"
#include "irl/geometry/general/rotations.h"
// #include "irl/paraboloid_reconstruction/gradient_paraboloid.h"
// #include "irl/paraboloid_reconstruction/hessian_paraboloid.h"

#include <sys/time.h>
#include <cmath>
#include <random>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/half_edge_structures/half_edge_polyhedron_quadratic.h"
#include "irl/geometry/half_edge_structures/segmented_half_edge_polyhedron_quadratic.h"
#include "irl/geometry/polyhedrons/general_polyhedron.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/planar_reconstruction/planar_separator.h"

using namespace IRL;

int main() {

  using quadPt = IRL::PtBase<__float128>;

  auto cube = StoredRectangularCuboid<quadPt>::fromBoundingPts(quadPt(0.0, 0.0, 0.0), quadPt(1.0, 1.0, 1.0));

  HalfEdgePolyhedronQuadratic<quadPt> half_edge;
  cube.setHalfEdgeVersion(&half_edge);
  auto seg_half_edge = half_edge.generateSegmentedPolyhedron();

  IRL::SegmentedHalfEdgePolyhedronQuadratic<IRL::FaceQuadratic<IRL::HalfEdgeQuadratic<IRL::VertexQuadratic<quadPt>>>, IRL::VertexQuadratic<quadPt>> P2;

  NoSurfaceOutput* surf = nullptr;
  nudgePolyhedron(&seg_half_edge, &half_edge,
                  0, surf);

  splitHalfEdgePolytope(
      &seg_half_edge, &P2, &half_edge,
      PlaneBase<__float128>(NormalBase<__float128>(0.0, 0.0, -1.0), 0.0));

  std:: cout << "here are some information about P2\n";
  std:: cout << "nb of faces " << P2.getNumberOfFaces() << "\n";
  std:: cout << "nb of vertex " << P2.getNumberOfVertices() << "\n";
  auto face0 = P2[0];
  std:: cout << "first face : \n " << face0 << "\n";
  auto first_half_edge = face0->getStartingHalfEdge();
  std:: cout << "first half edge : \n " << first_half_edge << "\n";
  auto previous_point = first_half_edge->getPreviousVertex();
  std:: cout << "vertex : \n " << previous_point << "\n";
  std:: cout << "vertex position pointer : \n " << &(previous_point->getLocation()) << "\n";
  std:: cout << "vertex position : \n " << previous_point->getLocation() << "\n";
  return 0;
}
