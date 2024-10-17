#include "examples/2d_advector/edge_transport.h"

namespace ASHISH {

Vector EdgeTransport::vertical_func_order_1(double t, Point P) {
  return {0., P.x};
}

Vector EdgeTransport::vertical_func_order_2(double t, Point P) {
  return {0., P.x * P.x};
}

Vector EdgeTransport::vertical_func_order_3(double t, Point P) {
  return {0., P.x * P.x * P.x};
}

Two_By_Two_Matrix EdgeTransport::vertical_func_3_gradient(double t, Point u) {
  return {{{0., 0.}, {3 * u.x * u.x, 0.}}};
}

Vector EdgeTransport::tangent_rate_func(double t, Vector n,
                                        Two_By_Two_Matrix gradient) {
  double partials_x =
      gradient.matrix.first.first * n.u + gradient.matrix.first.second * n.v;
  double partials_y =
      gradient.matrix.second.first * n.u + gradient.matrix.second.second * n.v;

  return {partials_x, partials_y};
}

Vector EdgeTransport::owkes_func(double t, Point P) {
  double pi = 2 * asin(1.0);
  double v_x =
      -2. * sin(pi * P.x) * sin(pi * P.x) * sin(pi * P.y) * cos(pi * P.y);
  double v_y =
      2. * sin(pi * P.y) * sin(pi * P.y) * sin(pi * P.x) * cos(pi * P.x);

  return {v_x, v_y};
}

std::vector<double> EdgeTransport::linspace(double start, double end) {
  double diff = end - start;
  std::vector<double> points_along_edge;

  for (int i = 0; i < Nx; i++) {
    points_along_edge.push_back(start + diff * i / (Nx - 1));
  }

  return points_along_edge;
}

std::pair<bool, double> EdgeTransport::rays_intersect(Point_N P0, Point_N P1) {
  // likely not problematic
  // checked over theory with Evrard
  // honestly should probably check again by hand
  double tol = std::numeric_limits<double>::epsilon() * 1000.;

  double dx = P1.P.x - P0.P.x;
  double dy = P1.P.y - P0.P.y;

  double det = P1.N.u * P0.N.v - P1.N.v * P0.N.u;
  if (std::abs(det) < tol) {
    double det_2 = P0.N.u * dy - P0.N.v * dx;
    if (std::abs(det_2) < tol) {
      double dot_product = P0.N.u * dx + P0.N.v * dy;
      if (dot_product > tol) {
        return std::pair<bool, double>(true,
                                       0.5 * std::sqrt(dx * dx + dy * dy));
      }
    }
    return std::pair<bool, double>(false, 0.0);
  }

  double t0 = (dy * P1.N.u - dx * P1.N.v) / det;
  double t1 = (dy * P0.N.u - dx * P0.N.v) / det;

  // std::cout << "t0: " << t0 << std::endl;
  // std::cout << "t1: " << t1 << std::endl;
  // std::cout << std::endl;

  return std::pair<bool, double>(t0 > 0.0 && t1 > 0.0, t0);
}

Point EdgeTransport::rk4(Point P, double dt, double total_time,
                         Vector (*func)(double t, Point P)) {
  Point P_new = P;

  Vector k1, k2, k3, k4;
  Point pk2, pk3, pk4;
  k1 = func(total_time, P_new);

  pk2 = {P_new.x + (dt / 2.) * k1.u, P_new.y + (dt / 2.) * k1.v};
  k2 = func(total_time + dt / 2., pk2);

  pk3 = {P_new.x + (dt / 2.) * k2.u, P_new.y + (dt / 2.) * k2.v};
  k3 = func(total_time + dt / 2., pk3);

  pk4 = {P_new.x + dt * k3.u, P_new.y + dt * k3.v};
  k4 = func(total_time + dt, pk4);

  P_new.x += (dt / 6.) * (k1.u + 2 * k2.u + 2 * k3.u + k4.u);
  P_new.y += (dt / 6.) * (k1.v + 2 * k2.v + 2 * k3.v + k4.v);

  return P_new;
}

std::vector<Point> EdgeTransport::rk4Vertices(
    Point P, double dt, double total_time, Vector (*func)(double t, Point P)) {
  Point P_new = P;
  std::vector<Point> P_vec;
  P_vec.push_back(P);

  Vector k1, k2, k3, k4;
  Point pk2, pk3, pk4;
  k1 = func(total_time, P_new);
  pk2 = {P_new.x + (dt / 2.) * k1.u, P_new.y + (dt / 2.) * k1.v};
  k2 = func(total_time + dt / 2., pk2);
  pk3 = {P_new.x + (dt / 2.) * k2.u, P_new.y + (dt / 2.) * k2.v};
  k3 = func(total_time + dt / 2., pk3);
  pk4 = {P_new.x + dt * k3.u, P_new.y + dt * k3.v};
  k4 = func(total_time + dt, pk4);

  P_new.x += (dt / 6.) * (k1.u + 2 * k2.u + 2 * k3.u + k4.u);
  P_new.y += (dt / 6.) * (k1.v + 2 * k2.v + 2 * k3.v + k4.v);

  P_vec.push_back(P_new);

  return P_vec;
}

Point_N EdgeTransport::rk4_with_tangents(
    Point_N Pn, double dt, double total_time, Vector (*func)(double t, Point P),
    Two_By_Two_Matrix (*func_2)(double t, Point P)) {
  Point P_new = Pn.P;
  Vector n_new = Pn.N;

  // Vectors or points?
  Vector k1, k2, k3, k4;
  Point pk2, pk3, pk4;
  Vector q1, q2, q3, q4;
  Point pq2, pq3, pq4;
  // Updating n_new
  q1 = tangent_rate_func(total_time, Pn.N, func_2(total_time, P_new));

  pq2 = {P_new.x + (dt / 2.) * q1.u, P_new.y + (dt / 2.) * q1.v};
  q2 = tangent_rate_func(total_time + dt / 2., Pn.N,
                         func_2(total_time + dt / 2., pq2));

  pq3 = {P_new.x + (dt / 2.) * q2.u, P_new.y + (dt / 2.) * q2.v};
  q3 = tangent_rate_func(total_time + dt / 2., Pn.N,
                         func_2(total_time + dt / 2., pq3));

  pq4 = {P_new.x + dt * q3.u, P_new.y + dt * q3.v};
  q4 = tangent_rate_func(total_time + dt, Pn.N, func_2(total_time + dt, pq3));

  n_new.u += (dt / 6.) * (q1.u + 2 * q2.u + 2 * q3.u + q4.u);
  n_new.v += (dt / 6.) * (q1.v + 2 * q2.v + 2 * q3.v + q4.v);

  // Updating u_new
  k1 = func(total_time, P_new);

  pk2 = {P_new.x + (dt / 2.) * k1.u, P_new.y + (dt / 2.) * k1.v};
  k2 = func(total_time + dt / 2., pk2);

  pk3 = {P_new.x + (dt / 2.) * k2.u, P_new.y + (dt / 2.) * k2.v};
  k3 = func(total_time + dt / 2., pk3);

  pk4 = {P_new.x + dt * k3.u, P_new.y + dt * k3.v};
  k4 = func(total_time + dt, pk4);

  P_new.x += (dt / 6.) * (k1.u + 2 * k2.u + 2 * k3.u + k4.u);
  P_new.y += (dt / 6.) * (k1.v + 2 * k2.v + 2 * k3.v + k4.v);

  return {P_new, n_new};
}

transport_edge_return EdgeTransport::transport_edge(
    Point P0, Point P1, double dt, double total_time,
    Vector (*func)(double t, Point P)) {
  double x_start = P0.x;
  double y_start = P0.y;
  double x_end = P1.x;
  double y_end = P1.y;
  std::vector<Point> start_pathline;
  std::vector<Point> transported_edge;
  std::vector<Point> end_pathline;
  Point point_to_add;

  // Uses linspace, not necessary when including tangents
  std::vector<double> x_range = linspace(x_start, x_end);
  std::vector<double> y_range = linspace(y_start, y_end);

  for (int i = 0; i < Nx; i++) {
    Point curr_point = {x_range[i], y_range[i]};
    if (i == 0) {
      start_pathline = rk4Vertices(curr_point, dt, total_time, func);
      transported_edge.push_back(start_pathline.back());
    } else if (i == Nx - 1) {
      end_pathline = rk4Vertices(curr_point, dt, total_time, func);
      transported_edge.push_back(end_pathline.back());
    } else {
      point_to_add = rk4(curr_point, dt, total_time, func);
      transported_edge.push_back(point_to_add);
    }
  }

  return {start_pathline, transported_edge, end_pathline};
}

transport_edge_return_tangents EdgeTransport::transport_edge_with_tangents(
    Point_N P0, Point_N P1, double dt, double total_time,
    Vector (*func)(double t, Point P),
    Two_By_Two_Matrix (*func_2)(double t, Point P), int recursion_num,
    std::vector<Point_N>& transported_edge, bool add_start, bool add_end) {
  // Note: include control points in 'Point_N' struct as part of output

  // Overwriting tangent vectors with normalized tangent vectors from P0 to
  // P1, and vice versa
  P0.N.u = P1.P.x - P0.P.x;
  P0.N.v = P1.P.y - P0.P.y;
  double norm = std::sqrt(P0.N.u * P0.N.u + P0.N.v * P0.N.v);
  P0.N.u /= norm;
  P0.N.v /= norm;
  P1.N.u = -P0.N.u;
  P1.N.v = -P0.N.v;

  std::vector<Point> start_pathline;
  std::vector<Point> end_pathline;
  Point point_to_add;

  /*
   * Note: Turn rk4Vertices into recursive function to well approximate
   * pathline Pass by reference 'start_pathline'
   */
  start_pathline = rk4Vertices(P0.P, dt, total_time, func);
  Point_N start_new_vertex =
      rk4_with_tangents(P0, dt, total_time, func, func_2);
  // std::cout << "Start New Vertex: P0 [(" << start_new_vertex.P.x << ", "
  //           << start_new_vertex.P.y << "), (" << start_new_vertex.N.u << ",
  //           "
  //           << start_new_vertex.N.v << ")]" << std::endl;

  // Add transported point P0 to transported edge
  if (add_start == true) {
    transported_edge.push_back(start_new_vertex);
  }

  /*
   * Note: Turn rk4Vertices into recursive function to well approximate
   * pathline Pass by reference 'end_pathline'
   */
  end_pathline = rk4Vertices(P1.P, dt, total_time, func);
  Point_N end_new_vertex = rk4_with_tangents(P1, dt, total_time, func, func_2);
  // std::cout << "End New Vertex: P1 [(" << end_new_vertex.P.x << ", "
  //           << end_new_vertex.P.y << "), (" << end_new_vertex.N.u << ", "
  //           << end_new_vertex.N.v << ")]" << std::endl;

  // Recursion
  auto not_need_vertex = rays_intersect(start_new_vertex, end_new_vertex);
  if (!not_need_vertex.first && recursion_num < 3) {
    Point_N P05;
    P05.P.x = (P0.P.x + P1.P.x) / 2.;
    P05.P.y = (P0.P.y + P1.P.y) / 2.;
    P05.N = P1.N;
    // std::cout << "Within recursion: P0 [(" << P0.P.x << ", " << P0.P.y <<
    // "),
    // ("
    //           << P0.N.u << ", " << P0.N.v << ")], P05 [(" << P05.P.x << ",
    //           "
    //           << P05.P.y << "), (" << P05.N.u << ", " << P05.N.v << ")]"
    //           << std::endl;
    transport_edge_return_tangents data_first_split =
        transport_edge_with_tangents(P0, P05, dt, total_time, func, func_2,
                                     recursion_num + 1, transported_edge, false,
                                     true);
    P05.N = P0.N;
    // std::cout << "P05 [(" << P05.P.x << ", " << P05.P.y << "), (" <<
    // P05.N.u
    //           << ", " << P05.N.v << ")], P1 [(" << P1.P.x << ", " << P1.P.y
    //           << "), (" << P1.N.u << ", " << P1.N.v << ")]" << std::endl;
    transport_edge_return_tangents data_second_split =
        transport_edge_with_tangents(P05, P1, dt, total_time, func, func_2,
                                     recursion_num + 1, transported_edge, false,
                                     false);
  }

  if (add_end == true) {
    // Always point tangents from P0 to P1
    end_new_vertex.flip_dir();
    transported_edge.push_back(end_new_vertex);
  }

  return {start_pathline, end_pathline};
}

}  // namespace ASHISH