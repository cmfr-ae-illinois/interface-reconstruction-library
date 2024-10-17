
#include "examples/2d_advector/functions.h"

namespace ASHISH {

/*---------------------------------------------------------------------------
 *	solve cubic equation x^3 + a*x^2 + b*x + c
 *	x - array of size 3
 *	In case 3 real roots: => x[0], x[1], x[2], return 3
 *	        2 real roots: x[0], x[1],          return 2
 *	        1 real root : x[0], x[1] ± i*x[2], return 1
 */
unsigned int VolumeFractionSolver::solveP3(double* x, double a, double b,
                                           double c) {
  double a2 = a * a;
  double q = (a2 - 3 * b) / 9;
  double r = (a * (2 * a2 - 9 * b) + 27 * c) / 54;
  double r2 = r * r;
  double q3 = q * q * q;
  double A, B;
  if (r2 < q3) {
    double t = r / sqrt(q3);
    if (t < -1) t = -1;
    if (t > 1) t = 1;
    t = acos(t);
    a /= 3;
    q = -2 * sqrt(q);
    x[0] = q * cos(t / 3) - a;
    x[1] = q * cos((t + M_2PI) / 3) - a;
    x[2] = q * cos((t - M_2PI) / 3) - a;
    return 3;
  } else {
    A = -pow(fabs(r) + sqrt(r2 - q3), 1. / 3);
    if (r < 0) A = -A;
    B = (0 == A ? 0 : q / A);

    a /= 3;
    x[0] = (A + B) - a;
    x[1] = -0.5 * (A + B) - a;
    x[2] = 0.5 * sqrt(3.) * (A - B);
    if (fabs(x[2]) < eps) {
      x[2] = x[1];
      return 2;
    }

    return 1;
  }
}

// Solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d
std::vector<double> VolumeFractionSolver::solve_quartic(double a, double b,
                                                        double c, double d) {
  double a3 = -b;
  double b3 = a * c - 4. * d;
  double c3 = -a * a * d - c * c + 4. * b * d;

  double x3[3];
  unsigned int iZeroes = solveP3(x3, a3, b3, c3);

  double q1, q2, p1, p2, D, sqD, y;

  y = x3[0];
  // THE ESSENCE - choosing Y with maximal absolute value !
  if (iZeroes != 1) {
    if (std::fabs(x3[1]) > std::fabs(y)) y = x3[1];
    if (std::fabs(x3[2]) > std::fabs(y)) y = x3[2];
  }

  // h1+h2 = y && h1*h2 = d  <=>  h^2 -y*h + d = 0    (h === q)

  D = y * y - 4 * d;
  if (std::fabs(D) < eps) {  // in other words - D==0
    q1 = q2 = y * 0.5;
    // g1+g2 = a && g1+g2 = b-y   <=>   g^2 - a*g + b-y = 0    (p === g)
    D = a * a - 4 * (b - y);
    if (std::fabs(D) < eps)  // in other words - D==0
      p1 = p2 = a * 0.5;

    else {
      sqD = std::sqrt(D);
      p1 = (a + sqD) * 0.5;
      p2 = (a - sqD) * 0.5;
    }
  } else {
    sqD = std::sqrt(D);
    q1 = (y + sqD) * 0.5;
    q2 = (y - sqD) * 0.5;
    // g1+g2 = a && g1*h2 + g2*h1 = c       ( && g === p )  Krammer
    p1 = (a * q1 - c) / (q1 - q2);
    p2 = (c - a * q2) / (q1 - q2);
  }

  std::vector<double> retval;

  // solving quadratic eq. - x^2 + p1*x + q1 = 0
  D = p1 * p1 - 4 * q1;
  if (D < 0.0 &&
      std::abs(std::sqrt(-D) * 0.5) < std::numeric_limits<double>::epsilon()) {
    retval.push_back(-p1 * 0.5);
  } else {
    sqD = sqrt(D);
    retval.push_back((-p1 + sqD) * 0.5);
    retval.push_back((-p1 - sqD) * 0.5);
  }

  // solving quadratic eq. - x^2 + p2*x + q2 = 0
  D = p2 * p2 - 4. * q2;
  if (D < 0.0 && abs(sqrt(-D) * 0.5) < std::numeric_limits<double>::epsilon()) {
    retval.push_back(-p2 * 0.5);
  } else {
    sqD = std::sqrt(D);
    retval.push_back((-p2 + sqD) * 0.5);
    retval.push_back((-p2 - sqD) * 0.5);
  }

  return retval;
}
// end of 4th order solver

std::vector<double> VolumeFractionSolver::solve_cubic(double a, double b,
                                                      double c, double d) {
  double delta_zero = b * b - 3 * a * c;
  double delta_one = 2 * b * b * b - 9 * a * b * c + 27 * a * a * d;
  double tol = std::numeric_limits<double>::epsilon();
  std::vector<double> retval;

  if (delta_zero <= 0 && delta_one <= 0) {
    retval.push_back(-b / (3 * a));
    retval.push_back(-b / (3 * a));
    retval.push_back(-b / (3 * a));
  } else {
    double C_part = std::pow(
        delta_one * delta_one - 4 * delta_zero * delta_zero * delta_zero,
        1 / 2.);
    double C1 = std::pow((delta_one + C_part) / 2., 1 / 3.);
    double C2 = std::pow((delta_one - C_part) / 2., 1 / 3.);
    double C;

    if (C1 != 0) {
      C = C1;
    } else {
      C = C2;
    }

    double epsilon = (-1. + std::pow(-3, 1 / 2.)) / 2.;
    std::vector<std::complex<double>> roots;
    for (int i = 0; i < 3; i++) {
      double roots_part = C * std::pow(epsilon, i);
      roots.push_back(-(b + roots_part + delta_zero / roots_part) / (3. * a));
    }
    for (const std::complex<double>& root : roots) {
      if (root.imag() < tol) {
        retval.push_back(root.real());
      }
    }
  }

  return retval;
}

double VolumeFractionSolver::GeneralFunc(double x, double a, double b,
                                         double c) {
  return a * x * x + b * x + c;
}

std::vector<Point_T> VolumeFractionSolver::AnalyticIntersections(
    double a, double b, double c, Point point_one, Point point_two,
    Point control) {
  double x1 = point_one.x;
  double y1 = point_one.y;
  double x2 = point_two.x;
  double y2 = point_two.y;
  double xc = control.x;
  double yc = control.y;
  std::vector<double> t_vals;

  // parabola-parabola case
  double A = a * (x1 * x1 + 2 * x1 * x2 - 4 * x1 * xc + x2 * x2 - 4 * x2 * xc +
                  4 * xc * xc);
  double B = a * (-4 * x1 * x1 - 4 * x1 * x2 + 12 * x1 * xc + 4 * x2 * xc -
                  8 * xc * xc);
  double C = a * (6 * x1 * x1 + 2 * x1 * x2 - 12 * x1 * xc + 4 * xc * xc) +
             b * (x1 + x2 - 2 * xc) - y1 - y2 + 2 * yc;
  double D = 4 * a * x1 * (xc - x1) + 2 * b * (xc - x1) + 2 * y1 - 2 * yc;
  double E = a * x1 * x1 + b * x1 + c - y1;
  int eqn_order = 0;

  std::vector<double> t_solutions;
  if (abs(A) > tol) {
    // if A != 0, solve quartic equation
    t_solutions = solve_quartic(B / A, C / A, D / A, E / A);
    eqn_order = 4;
  } else if (abs(B) > tol) {
    // solve cubic equation
    t_solutions = solve_cubic(B, C, D, E);
    eqn_order = 3;
  } else if (abs(C) > tol) {
    // solve quadratic equation
    t_solutions.push_back((-D - sqrt(D * D - 4 * C * E)) / (2 * C));
    t_solutions.push_back((-D + sqrt(D * D - 4 * C * E)) / (2 * C));
    eqn_order = 2;
  } else if (abs(D) > tol) {
    // solve linear equation
    t_solutions.push_back(-E / D);
    eqn_order = 1;
  } else {
    t_solutions = {};
    eqn_order = 0;
  }

  for (size_t i = 0; i < eqn_order; ++i) {
    if (t_solutions[i] > tol && t_solutions[i] < (1. - tol)) {
      t_vals.push_back(t_solutions[i]);
    }
  }

  std::sort(t_vals.begin(), t_vals.end());

  std::vector<Point_T> intersections;

  for (auto t : t_vals) {
    double B0 = (1 - t) * (1 - t);
    double B1 = 2 * (1 - t) * t;
    double B2 = t * t;

    double x = B0 * x1 + B1 * xc + B2 * x2;
    double y = B0 * y1 + B1 * yc + B2 * y2;

    intersections.push_back({{x, y}, t});
  }

  return intersections;
}

double VolumeFractionSolver::TotalDivergenceThm(Point current, Point next,
                                                Point control) {
  double x1 = current.x;
  double y1 = current.y;
  double x2 = next.x;
  double y2 = next.y;
  double xc = control.x;
  double yc = control.y;

  return 1. / 6. *
         (-x1 * (2. * yc + y2) + 2. * xc * (y1 - y2) + x2 * (y1 + 2. * yc));
}

double VolumeFractionSolver::DivergenceThm(Point current, Point next,
                                           Point control, double ct,
                                           double nt) {
  double x1 = current.x;
  double y1 = current.y;
  double x2 = next.x;
  double y2 = next.y;
  double xc = control.x;
  double yc = control.y;
  double t1 = ct;
  double t2 = nt;

  double part1 =
      t1 * t1 * t1 - 3. * t1 * t1 + 3. * t1 - t2 * (t2 * t2 - 3. * t2 + 3.);
  double part2 = 2. * t1 * t1 * t1 - 3. * t1 * t1 - t2 * t2 * (2. * t2 - 3.);
  double part3 = t2 * t2 * t2 - t1 * t1 * t1;

  return 1. / 6. *
         (2. * x1 * yc * part1 - x1 * y2 * part2 -
          2. * xc * (y1 * part1 + y2 * part3) + x2 * y1 * part2 +
          2. * x2 * yc * part3);
}

Point VolumeFractionSolver::PointFinder(CoordArray* array, int point_num) {
  CoordArray* current_array = array;
  for (int i = 0; i < point_num; i++) {
    current_array = current_array->next;
  }
  return current_array->data;
}

void VolumeFractionSolver::VertexRemover(PTS_LinkedList& vertices, double a,
                                         double b, double c) {
  HalfEdge* head = vertices.getHead();
  HalfEdge* tail = vertices.getTail();
  HalfEdge* current = head;

  while (current->next != head) {
    // If intersection (t != 2), skip remover
    if (current->next->data.t != 2) {
      current = current->next;
      continue;
    }

    double x = current->next->data.P.x;
    double y = current->next->data.P.y;
    if (y > GeneralFunc(x, a, b, c)) {
      if (current->next == tail) {
        vertices.deleteLastNode();
      } else {
        HalfEdge* temp = current->next;
        current->next = current->next->next;
        delete temp;
      }
    } else {
      current = current->next;
    }
  }

  vertices.printList_t();
}

std::pair<double, P_LinkedList> VolumeFractionSolver::VolumeFraction(
    P_LinkedList vertices, P_LinkedList controls, double a, double b,
    double c) {
  // Defining PTS_LinkedList 'new_vertices'
  PTS_LinkedList new_vertices;
  CoordArray* vertices_head = vertices.getHead();  // head of 'vertices' pointer
  CoordArray* vertices_current =
      vertices_head;  // steps through 'vertices' pointer
  CoordArray* controls_head = controls.getHead();  // head of 'controls' pointer
  CoordArray* control_current =
      controls_head;  // steps through 'controls' pointer
  int side_num = 0;
  bool shift_vertices = true;
  bool last = false;

  while (last == false) {
    if (vertices_current->data.y <=
        GeneralFunc(vertices_current->data.x, a, b, c)) {
      Point_T_sides to_add = {vertices_current->data, 2, side_num};
      new_vertices.setTail(to_add);
      if (side_num == 0) {
        shift_vertices = false;
      }
    }
    std::vector<Point_T> intersections = AnalyticIntersections(
        a, b, c, vertices_current->data, vertices_current->next->data,
        control_current->data);

    for (int i = 0; i < intersections.size(); i++) {
      Point_T_sides to_add = {intersections[i].P, intersections[i].t, side_num};
      new_vertices.setTail(to_add);
    }

    // breaks the loop
    if (vertices_current->next == vertices_head) {
      last = true;
    }

    // moves to next point
    control_current = control_current->next;
    vertices_current = vertices_current->next;
    side_num++;
  }

  // Shift vertices if necessary
  if (shift_vertices) {
    new_vertices.shiftNodesBack();
  }

  // VertexRemover(new_vertices, a, b, c);

  // SKIPPING VF = 0 case
  // SKIPPING starting vertex outside CB case

  // Calculating area within bounds
  double area_within_bounds = 0;
  HalfEdge* new_vertices_head = new_vertices.getHead();
  HalfEdge* new_vertices_current = new_vertices_head;

  // Extract data from current node in 'new_vertices'
  Point current_coords = new_vertices_current->data.P;
  double current_t = new_vertices_current->data.t;
  int current_side = new_vertices_current->data.side_num;

  Point point_current_vertex = PointFinder(vertices_head, current_side);
  Point point_control_current = PointFinder(controls_head, current_side);
  bool use_parabola = (current_t == 2) ? false : true;
  last = false;

  Point next_coords;
  Point next_vertex;
  double next_t;

  double x0;
  double x1;
  double to_add;

  while (last == false) {
    if (new_vertices_current->next == new_vertices_head) {
      last = true;
      next_coords = new_vertices_head->data.P;
      next_vertex = PointFinder(vertices_head, 0);
      next_t = new_vertices_head->data.t;
    } else {
      next_coords = new_vertices_current->next->data.P;
      next_vertex = PointFinder(vertices_head, current_side + 1);
      next_t = new_vertices_current->next->data.t;
    }

    if (use_parabola) {
      use_parabola = false;
      x0 = current_coords.x;
      x1 = next_coords.x;
      to_add =
          1. / 6. *
          (-3. * c * x0 + a * x0 * x0 * x0 + 3. * c * x1 - a * x1 * x1 * x1);
      area_within_bounds += to_add;
      // std::cout << "Entered Parabola " << to_add << std::endl;
    } else {
      if (next_t != 2) {
        use_parabola = true;
      }
      if (current_t != 2) {
        if (next_t != 2) {
          to_add = DivergenceThm(point_current_vertex, next_vertex,
                                 point_control_current, current_t, next_t);
        } else {
          to_add = DivergenceThm(point_current_vertex, next_vertex,
                                 point_control_current, current_t, 1.);
        }
      } else if (next_t != 2) {
        to_add = DivergenceThm(point_current_vertex, next_vertex,
                               point_control_current, 0., next_t);
      } else {
        to_add = TotalDivergenceThm(point_current_vertex, next_vertex,
                                    point_control_current);
      }
      area_within_bounds += to_add;
      // std::cout << "Entered Divergence " << to_add << std::endl;
    }

    new_vertices_current = new_vertices_current->next;
    current_coords = new_vertices_current->data.P;
    current_t = new_vertices_current->data.t;
    current_side = new_vertices_current->data.side_num;
    point_current_vertex = PointFinder(vertices_head, current_side);
    point_control_current = PointFinder(controls_head, current_side);
  }

  // Calculating area
  double area = 0;
  vertices_current = vertices_head;
  control_current = controls_head;
  point_current_vertex = vertices_current->data;
  point_control_current = control_current->data;
  last = false;

  while (last == false) {
    if (vertices_current->next == vertices_head) {
      last = true;
      next_vertex = vertices_head->data;
    } else {
      next_vertex = vertices_current->next->data;
    }

    to_add = TotalDivergenceThm(point_current_vertex, next_vertex,
                                point_control_current);
    area += to_add;

    vertices_current = vertices_current->next;
    control_current = control_current->next;
    point_current_vertex = vertices_current->data;
    point_control_current = control_current->data;
  }

  std::cout << "P Linked List: ";
  new_vertices.printList_t();

  // std::cout << "area: " << area << std::endl;

  P_LinkedList final_new_vertices;
  final_new_vertices.setTail({2, 3});
  double volume_frac = area_within_bounds / area;

  return {volume_frac, final_new_vertices};
}

}  // namespace ASHISH