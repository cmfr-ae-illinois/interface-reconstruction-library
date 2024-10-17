#ifndef EXAMPLES_2D_ADVECTOR_EDGE_TRANSPORT_H_
#define EXAMPLES_2D_ADVECTOR_EDGE_TRANSPORT_H_

#include <cmath>
#include <iostream>
#include <string>
#include "examples/2d_advector/functions.h"

namespace ASHISH {

struct transport_edge_return {
  std::vector<Point> start_pathline;
  std::vector<Point> transported_edge;
  std::vector<Point> end_pathline;
};

struct transport_edge_return_tangents {
  std::vector<Point> start_pathline;
  // std::vector<Point_N> transported_edge;
  std::vector<Point> end_pathline;
};

struct Two_By_Two_Matrix {
  std::pair<std::pair<double, double>, std::pair<double, double>> matrix;
};

class EdgeTransport {
 public:
  // Test velocity fields
  Vector static vertical_func_order_1(double t, Point P);
  Vector static vertical_func_order_2(double t, Point P);
  Vector static vertical_func_order_3(double t, Point P);
  Two_By_Two_Matrix static vertical_func_1_gradient(double t, Point u);
  Two_By_Two_Matrix static vertical_func_2_gradient(double t, Point u);
  Two_By_Two_Matrix static vertical_func_3_gradient(double t, Point u);
  Vector static owkes_func(double t, Point P);
  Vector static tangent_rate_func(double t, Vector n,
                                  Two_By_Two_Matrix gradient);

  // Helper functions
  std::vector<double> static linspace(double start, double end);
  // True: rays do intersect, recursion not necessary
  // False: rays do NOT intersect, recursion necessary
  std::pair<bool, double> static rays_intersect(Point_N P0, Point_N P1);

  // RK4 functions
  Point static rk4(Point P, double dt, double total_time,
                   Vector (*func)(double t, Point P));
  std::vector<Point> static rk4Vertices(Point P, double dt, double total_time,
                                        Vector (*func)(double t, Point P));
  Point_N static rk4_with_tangents(Point_N Pn, double dt, double total_time,
                                   Vector (*func)(double t, Point P),
                                   Two_By_Two_Matrix (*func_2)(double t,
                                                               Point P));

  // Main functions
  transport_edge_return static transport_edge(Point P0, Point P1, double dt,
                                              double total_time,
                                              Vector (*func)(double t,
                                                             Point P));
  transport_edge_return_tangents static transport_edge_with_tangents(
      Point_N P0, Point_N P1, double dt, double total_time,
      Vector (*func)(double t, Point P),
      Two_By_Two_Matrix (*func_2)(double t, Point P), int recursion_num,
      std::vector<Point_N>& transported_edge, bool add_start = true,
      bool add_end = true);

 private:
  static const int Nx = 5;  // number of spatial interpolation points,
                            // including start and end points
  // static double t_steps = 5.; // number of temporal interpolation points,
  // including start and end points int total_time = 0.5;
};

}  // namespace ASHISH

// #include "examples/2d_advector/edge_transport.tpp"

#endif  // EXAMPLES_2D_ADVECTOR_EDGE_TRANSPORT_H_