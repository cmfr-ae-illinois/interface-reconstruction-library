#ifndef EXAMPLES_2D_ADVECTOR_FUNCTIONS_H_
#define EXAMPLES_2D_ADVECTOR_FUNCTIONS_H_

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace ASHISH {

// coordinates
struct Point {
  double x, y;
};

struct PointAndControl {
  Point P, C;
};

struct Parabola {
  double x, y, curv, angl;
};

using BezierList = std::vector<PointAndControl>;

// vector
struct Vector {
  double u, v;
};

// Point with tangential direction
// Used for edge transports
struct Point_N {
  Point P;
  Vector N;

  bool operator!=(const Point_N& other) const {
    return (P.x != other.P.x || P.y != other.P.y || N.u != other.N.u ||
            N.v != other.N.v);
  }

  void flip_dir(void) {
    N.u = -N.u;
    N.v = -N.v;
  }
};

// Point with t value struct
// Used for AnalyticIntersections
struct Point_T {
  Point P;
  double t;
};

struct Point_T_sides {
  Point P;
  double t;
  int side_num;

  // Constructor
  Point_T_sides(Point P_val, double t_val, int side_num_val)
      : P(P_val), t(t_val), side_num(side_num_val) {}
};

struct HalfEdge {  // Half edge Linked list
  Point_T_sides data;
  // Point control_point; - TO IMPLEMENT
  HalfEdge* next;

  // Constructor
  HalfEdge(Point_T_sides pts) : data(pts), next(nullptr) {}
};

// Define the LinkedList class
class PTS_LinkedList {
 public:
  PTS_LinkedList() = default;

  // Method to add a new node at the end of the list
  void append(Point_T_sides value);

  // Setters and getters
  void setHead(Point_T_sides value);
  void setTail(Point_T_sides value);
  HalfEdge* getHead();
  HalfEdge* getTail();

  void deleteLastNode();
  void shiftNodesBack();

  // Method to print the list of t values
  void printList_t() const;

  // Destructor to free the allocated memory
  ~PTS_LinkedList();

 private:
  HalfEdge* head;
  HalfEdge* tail;
};

// nodes for Point Linked list
struct CoordArray {
  Point data;
  CoordArray* next;

  // Constructor
  CoordArray(Point pts) : data(pts), next(nullptr) {}
};

// Point Linked list
class P_LinkedList {
 public:
  P_LinkedList() = default;

  // Method to add a new node at the end of the list
  void append(Point value);

  // Setters and getters
  void setHead(Point value);
  void setTail(Point value);
  CoordArray* getHead();
  CoordArray* getTail();

  void deleteLastNode();

  // Method to print the list of t values
  void printList_t() const;

  // Destructor to free the allocated memory
  // ~P_LinkedList();

 private:
  CoordArray* head;
  CoordArray* tail;
};

// define the volume fraction solver class
class VolumeFractionSolver {
 public:
  // General parabola function
  double static GeneralFunc(double x, double a, double b, double c);

  /*
   * FUNCTION: Analytically solves for the intersection points between a
   * parabola and the curved cell edge INPUTS: (a, b, c) - intersecting
   * parabola coefficients (point_one, point_two) - cell boundary vertices
   * (control) - control point for parabola OUTPUT: vector of intersection
   * points with their corresponding t value: (Point (x, y)) - intersection
   * point (double t) - corresponding t value
   */
  std::vector<Point_T> static AnalyticIntersections(double a, double b,
                                                    double c, Point point_one,
                                                    Point point_two,
                                                    Point control);

  // FUNCTIONS: Employs divergence theorem given 2 points, a control point,
  // and their corresponding t values
  double static TotalDivergenceThm(Point current, Point next, Point control);
  double static DivergenceThm(Point current, Point next, Point control,
                              double ct, double nt);

  /*
   *   FUNCTION: Removes vertices that fall outside of the bounds defined by
   * the parabola INPUTS: head of list of vertices and 'a, b, c' coefficients
   * of a parabola OUTPUTS: head of list of new vertices that lie within the
   * bounds defined by the parabola
   */
  void VertexRemover(PTS_LinkedList& vertices, double a, double b,
                     double c);  // INCOMPLETE

  // FUNCTION: Analytically solves for the volume fraction of the region
  // bounded by a parabola within a cell boundary
  std::pair<double, P_LinkedList> static VolumeFraction(
      P_LinkedList vertices, P_LinkedList controls, double a, double b,
      double c);  // INCOMPLETE

  Point static PointFinder(CoordArray* array, int point_num);

  /*
   *   FUNCTIONS: Analytically solves for the roots of a quartic equation
   * (sasamil solver) INPUTS: (a, b, c, d) - coefficients of the equation x^4
   * + ax^3 + bx^2 + cx + d = 0
   */
  std::vector<double> static solve_quartic(double a, double b, double c,
                                           double d);
  unsigned int static solveP3(double* x, double a, double b, double c);

  std::vector<double> static solve_cubic(double a, double b, double c,
                                         double d);

 private:
  static constexpr double tol = std::numeric_limits<double>::epsilon();
  // Sasamil Solver - https://github.com/sasamil/Quartic
  static constexpr double PI = 3.141592653589793238463L;
  static constexpr double M_2PI = 2 * PI;
  static constexpr double eps = 1e-12;
};

}  // namespace ASHISH

// #include "examples/2d_advector/linked_list.tpp"
// //
// #include "examples/2d_advector/functions.tpp"

#endif  // EXAMPLES_2D_ADVECTOR_FUNCTIONS_H_
