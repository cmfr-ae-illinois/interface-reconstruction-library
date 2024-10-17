#include "examples/2d_advector/irl2d.h"

namespace IRL2D {

double operator*(const Vec& a_vec_0, const Vec& a_vec_1) {
  return a_vec_0[0] * a_vec_1[0] + a_vec_0[1] * a_vec_1[1];
}
Vec operator+(const Vec& a_vec_0, const Vec& a_vec_1) {
  return Vec(a_vec_0[0] + a_vec_1[0], a_vec_0[1] + a_vec_1[1]);
}
Vec operator-(const Vec& a_vec_0, const Vec& a_vec_1) {
  return Vec(a_vec_0[0] - a_vec_1[0], a_vec_0[1] - a_vec_1[1]);
}
Vec operator*(const double a_scalar, const Vec& a_vec) {
  auto vec_to_return = a_vec;
  vec_to_return *= a_scalar;
  return vec_to_return;
}
Vec operator*(const Vec& a_vec, const double a_scalar) {
  return a_scalar * a_vec;
}
Vec operator/(const Vec& a_vec, const double a_scalar) {
  Vec vec_to_return = a_vec;
  vec_to_return /= a_scalar;
  return vec_to_return;
}
Vec operator*(const Mat& a_mat, const Vec& a_vec) {
  return Vec(a_mat[0] * a_vec, a_mat[1] * a_vec);
}

void Print(const BezierList& list) {
  std::cout << "Cell:";
  for (int i = 0; i < list.size(); i++) {
    std::cout << "\n  Vec" << i << "  = " << list[i].first;
    std::cout << "\n  Ctrl" << i << " = " << list[i].second;
  }
  std::cout << std::endl;
}

void ToVTK(const std::vector<BezierList>& list, const std::string& filename) {
  const int nsamples = 100;
  const double w = 1. / static_cast<double>(nsamples);
  int npoints = 0;
  for (int i = 0; i < list.size(); i++) {
    npoints += nsamples * list[i].size();
  }
  std::ofstream file;
  file.open(filename + std::string(".vtu"));
  file << "<VTKFile type=\"UnstructuredGrid\">\n<UnstructuredGrid>\n";
  file << "<Piece NumberOfPoints=\"" << npoints << "\" NumberOfCells=\""
       << list.size() << "\">\n";
  file << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n";
  for (int i = 0; i < list.size(); i++) {
    for (int n = 0; n < list[i].size(); n++) {
      const auto P0 = list[i][n].first;
      const auto P1 = list[i][n].second;
      const auto P2 = list[i][(n + 1) % list[i].size()].first;
      for (int m = 0; m < nsamples; m++) {
        const double t = static_cast<double>(m) * w;
        const auto P =
            P0 * (1. - t) * (1. - t) + P1 * 2. * (1 - t) * t + P2 * t * t;
        file << P.x() << " " << P.y() << " 0. \n";
      }
    }
  }
  file << "</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int32\" "
          "Name=\"connectivity\" format=\"ascii\">\n";
  int offset = 0;
  for (int i = 0; i < list.size(); i++) {
    int count = 0;
    for (int n = 0; n < list[i].size(); n++) {
      for (int m = 0; m < nsamples; m++) {
        file << offset + count++ << " ";
      }
    }
    offset = count;
  }
  file << "\n</DataArray>\n";
  file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  int count = 0;
  for (int i = 0; i < list.size(); i++) {
    count += nsamples * list[i].size();
    file << count << " ";
  }
  file << "\n</DataArray>\n";
  file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (int i = 0; i < list.size(); i++) {
    file << "7 ";
  }
  file << "\n</DataArray>\n";
  file << "</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
  file.close();

  npoints = 0;
  for (int i = 0; i < list.size(); i++) {
    npoints += 2 * list[i].size();
  }
  file.open(filename + std::string("_skeleton.vtu"));
  file << "<VTKFile type=\"UnstructuredGrid\">\n<UnstructuredGrid>\n";
  file << "<Piece NumberOfPoints=\"" << npoints << "\" NumberOfCells=\""
       << list.size() << "\">\n";
  file << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n";
  for (int i = 0; i < list.size(); i++) {
    for (int n = 0; n < list[i].size(); n++) {
      const auto P0 = list[i][n].first;
      const auto P1 = list[i][n].second;
      file << P0.x() << " " << P0.y() << " 0. \n";
      file << P1.x() << " " << P1.y() << " 0. \n";
    }
  }
  file << "</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int32\" "
          "Name=\"connectivity\" format=\"ascii\">\n";
  offset = 0;
  for (int i = 0; i < list.size(); i++) {
    int count = 0;
    for (int n = 0; n < list[i].size(); n++) {
      file << offset + count++ << " ";
      file << offset + count++ << " ";
    }
    offset = count;
  }
  file << "\n</DataArray>\n";
  file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  count = 0;
  for (int i = 0; i < list.size(); i++) {
    count += 2 * list[i].size();
    file << count << " ";
  }
  file << "\n</DataArray>\n";
  file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (int i = 0; i < list.size(); i++) {
    file << "7 ";
  }
  file << "\n</DataArray>\n";
  file << "</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
  file.close();
}

void ToVTK(const BezierList& list, const std::string& filename) {
  ToVTK(std::vector<BezierList>{list}, filename);
}

BezierList RectangleListFromBounds(const double x0, const double x1,
                                   const double y0, const double y1) {
  Vec p00(x0, y0), p10(x1, y0), p11(x1, y1), p01(x0, y1);
  PtAndControl pc0{p00, .5 * (p00 + p10)};
  PtAndControl pc1{p10, .5 * (p10 + p11)};
  PtAndControl pc2{p11, .5 * (p11 + p01)};
  PtAndControl pc3{p01, .5 * (p01 + p00)};
  return BezierList{pc0, pc1, pc2, pc3};
}

unsigned int solveP3(double* x, const double a, const double b,
                     const double c) {
  /*---------------------------------------------------------------------------
   *	solve cubic equation x^3 + a*x^2 + b*x + c
   *	x - array of size 3
   *	In case 3 real roots: => x[0], x[1], x[2], return 3
   *	        2 real roots: x[0], x[1],          return 2
   *	        1 real root : x[0], x[1] ± i*x[2], return 1
   */
  double a2 = a * a;
  double q = (a2 - 3 * b) / 9.;
  double r = (a * (2. * a2 - 9. * b) + 27. * c) / 54.;
  double r2 = r * r;
  double q3 = q * q * q;
  const double tol = std::numeric_limits<double>::epsilon();
  double A, B;
  if (r2 < q3) {
    double t = r / std::sqrt(q3);
    if (t < -1.) t = -1.;
    if (t > 1.) t = 1.;
    t = std::acos(t);
    q = -2. * std::sqrt(q);
    x[0] = q * std::cos(t / 3.) - a / 3.;
    x[1] = q * std::cos((t + 2. * M_PI) / 3.) - a / 3.;
    x[2] = q * std::cos((t - 2. * M_PI) / 3.) - a / 3.;
    return 3;
  } else {
    A = -std::pow(std::fabs(r) + std::sqrt(r2 - q3), 1. / 3.);
    if (r < 0.) A = -A;
    B = (0. == A ? 0. : q / A);
    x[0] = (A + B) - a / 3.;
    x[1] = -0.5 * (A + B) - a / 3.;
    x[2] = 0.5 * std::sqrt(3.) * (A - B);
    if (std::fabs(x[2]) < tol) {
      x[2] = x[1];
      return 2;
    }
    return 1;
  }
}

std::vector<double> solve_quartic(const double a, const double b,
                                  const double c, const double d) {
  // Solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d
  const double a3 = -b;
  const double b3 = a * c - 4. * d;
  const double c3 = -a * a * d - c * c + 4. * b * d;
  const double tol = std::numeric_limits<double>::epsilon();
  double x3[3];
  const unsigned int iZeroes = solveP3(x3, a3, b3, c3);
  double q1, q2, p1, p2, D, sqD, y;

  y = x3[0];
  // THE ESSENCE - choosing Y with maximal absolute value !
  if (iZeroes != 1) {
    if (std::fabs(x3[1]) > std::fabs(y)) y = x3[1];
    if (std::fabs(x3[2]) > std::fabs(y)) y = x3[2];
  }

  // h1+h2 = y && h1*h2 = d  <=>  h^2 -y*h + d = 0    (h === q)

  D = y * y - 4. * d;
  if (std::fabs(D) < tol) {  // in other words - D==0
    q1 = q2 = y * 0.5;
    // g1+g2 = a && g1+g2 = b-y   <=>   g^2 - a*g + b-y = 0    (p === g)
    D = a * a - 4. * (b - y);
    if (std::fabs(D) < tol)  // in other words - D==0
    {
      p1 = p2 = a * 0.5;
    } else {
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
  D = p1 * p1 - 4. * q1;
  if (D < 0.0 && std::abs(std::sqrt(-D) * 0.5) < tol) {
    retval.push_back(-p1 * 0.5);
  } else {
    sqD = sqrt(D);
    retval.push_back((-p1 + sqD) * 0.5);
    retval.push_back((-p1 - sqD) * 0.5);
  }

  // solving quadratic eq. - x^2 + p2*x + q2 = 0
  D = p2 * p2 - 4. * q2;
  if (D < 0.0 && abs(sqrt(-D) * 0.5) < tol) {
    retval.push_back(-p2 * 0.5);
  } else {
    sqD = std::sqrt(D);
    retval.push_back((-p2 + sqD) * 0.5);
    retval.push_back((-p2 - sqD) * 0.5);
  }

  return retval;
}

std::vector<double> solve_cubic(const double a, const double b, const double c,
                                const double d) {
  double delta_zero = b * b - 3. * a * c;
  double delta_one = 2. * b * b * b - 9. * a * b * c + 27. * a * a * d;
  double tol = std::numeric_limits<double>::epsilon();
  std::vector<double> retval;

  if (delta_zero <= 0 && delta_one <= 0) {
    retval.push_back(-b / (3. * a));
    retval.push_back(-b / (3. * a));
    retval.push_back(-b / (3. * a));
  } else {
    double C_part = std::pow(
        delta_one * delta_one - 4. * delta_zero * delta_zero * delta_zero,
        1. / 2.);
    double C1 = std::pow((delta_one + C_part) / 2., 1. / 3.);
    double C2 = std::pow((delta_one - C_part) / 2., 1. / 3.);
    double C;

    if (C1 != 0.) {
      C = C1;
    } else {
      C = C2;
    }

    double epsilon = (-1. + std::pow(-3., 1. / 2.)) / 2.;
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

std::vector<double> AnalyticIntersections(const Parabola& parabola,
                                          const Vec& p0, const Vec& p1,
                                          const Vec& p2) {
  double tol = 100. * std::numeric_limits<double>::epsilon();
  const double a = -parabola.coeff();
  const double x1 = p0.x(), y1 = p0.y(), x2 = p2.x(), y2 = p2.y(), xc = p1.x(),
               yc = p1.y();
  double A = a * (x1 * x1 + 2. * x1 * x2 - 4. * x1 * xc + x2 * x2 -
                  4. * x2 * xc + 4. * xc * xc);
  double B = a * (-4. * x1 * x1 - 4. * x1 * x2 + 12. * x1 * xc + 4. * x2 * xc -
                  8. * xc * xc);
  double C = a * (6. * x1 * x1 + 2. * x1 * x2 - 12. * x1 * xc + 4. * xc * xc) -
             y1 - y2 + 2. * yc;
  double D = 4. * a * x1 * (xc - x1) + 2. * y1 - 2. * yc;
  double E = a * x1 * x1 - y1;
  int eqn_order = 0;

  std::vector<double> t_vals, t_solutions;
  if (std::abs(A) > tol) {
    // if A != 0, solve quartic equation
    t_solutions = solve_quartic(B / A, C / A, D / A, E / A);
    // std::cout << "Quartic solve (ABCDE = " << A << ", " << B << ", " << C
    //           << ", " << D << ", " << E << std::endl;
    eqn_order = 4;
  } else if (std::abs(B) > tol) {
    // solve cubic equation
    t_solutions = solve_cubic(B, C, D, E);
    eqn_order = 3;
    // std::cout << "Cubic solve (ABCDE = " << A << ", " << B << ", " << C << ",
    // "
    //           << D << ", " << E << std::endl;
  } else if (std::abs(C) > tol) {
    // solve quadratic equation
    t_solutions.push_back((-D - std::sqrt(D * D - 4. * C * E)) / (2. * C));
    t_solutions.push_back((-D + std::sqrt(D * D - 4. * C * E)) / (2. * C));
    eqn_order = 2;
    // std::cout << "Quadratic solve (ABCDE = " << A << ", " << B << ", " << C
    //           << ", " << D << ", " << E << std::endl;
  } else if (std::abs(D) > tol) {
    // solve linear equation
    t_solutions.push_back(-E / D);
    eqn_order = 1;
    // std::cout << "Linear solve (ABCDE = " << A << ", " << B << ", " << C <<
    // ", "
    //           << D << ", " << E << std::endl;
  } else {
    t_solutions = {};
    eqn_order = 0;
    // std::cout << "No solve (ABCDE = " << A << ", " << B << ", " << C << ", "
    //           << D << ", " << E << std::endl;
  }

  for (size_t i = 0; i < eqn_order; ++i) {
    if (t_solutions[i] >= 0. && t_solutions[i] <= 1.) {
      t_vals.push_back(t_solutions[i]);
    }
  }

  std::sort(t_vals.begin(), t_vals.end());
  return t_vals;

  // std::vector<std::pair<Vec, double>> intersections;
  // for (auto t : t_vals) {
  //   double B0 = (1. - t) * (1. - t);
  //   double B1 = 2. * (1. - t) * t;
  //   double B2 = t * t;
  //   double x = B0 * x1 + B1 * xc + B2 * x2;
  //   double y = B0 * y1 + B1 * yc + B2 * y2;
  //   intersections.push_back(std::pair<Vec, double>{{x, y}, t});
  // }
  // return intersections;
}

const Vec BezierPoint(const Vec& p0, const Vec& p1, const Vec& p2,
                      const double t) {
  return p0 * (1. - t) * (1. - t) + p1 * 2. * (1 - t) * t + p2 * t * t;
}

bool IsBelow(const Parabola& parabola, const Vec& pt) {
  return pt.y() < -parabola.coeff() * pt.x() * pt.x();
}

std::vector<BezierList> ParabolaClip(const BezierList& original_cell,
                                     const Parabola& parabola) {
  // If empty cell, return empty cell
  if (original_cell.size() == 0 || parabola.isAlwaysBelow())
    return std::vector<BezierList>(0);
  if (parabola.isAlwaysAbove()) return std::vector<BezierList>{original_cell};

  // Specify constants
  double tol = std::numeric_limits<double>::epsilon();
  bool parabola_contains_vertex = true;
  const int itmax = 100;

  // Initialize list of clipped regions
  std::vector<BezierList> clipped_cell;
  // Reserve 5 times original size (i.e., 4 intersections per arc)
  clipped_cell.reserve(5 * original_cell.size());

  // Store parabola properties
  Vec datum = parabola.datum();
  ReferenceFrame frame = parabola.frame();
  double coeff = parabola.coeff();

  // This while loop is needed in case we might want to nudge the parabola
  int it = 0;
  while (parabola_contains_vertex && it < itmax) {
    parabola_contains_vertex = false;
    it++;
    // Copy cell in frame of reference of parabola
    const int nvert_init = original_cell.size();
    BezierList cell, cell_with_intersections, intersections;
    cell.assign(original_cell.begin(), original_cell.end());
    for (int i = 0; i < nvert_init; i++) {
      cell[i].first = frame * (cell[i].first - datum);
      cell[i].second = frame * (cell[i].second - datum);
    }

    // Compute all intersections and insert them in cell
    int nintersections = 0;
    std::vector<int>
        vertex_nature;  // 1 = below, 0 = above, 2 = entry, 3 = exit
    std::vector<std::pair<double, int>> sorted_list_for_mapping;
    for (int i = 0; i < nvert_init; i++) {
      // Get points of current bezier arc
      const auto p0 = cell[i].first;
      const auto p1 = cell[i].second;
      const auto p2 = cell[(i + 1) % cell.size()].first;
      // Add p0 and p1 to new list
      cell_with_intersections.push_back(PtAndControl{p0, p1});
      // Is p0 below or above parabola?
      vertex_nature.push_back(IsBelow(parabola, p0) ? 1 : 0);
      // Compute t-values of potential intersections
      auto t_vals = AnalyticIntersections(parabola, p0, p1, p2);
      // Tmp variables needs for arc splitting
      double t0 = 0.;
      auto p0new = p0, p1new = p1;
      // Loop over all detected intersections
      for (int j = 0; j < t_vals.size(); j++) {
        const double t = t_vals[j];
        // If t = 0 or t = 1, we must nudge! Go back to beginning of function
        if (std::fabs(t) < tol || std::fabs(1. - t) < tol) {
          parabola_contains_vertex = true;
          break;
        }
        //// Spilling of the arc (using de Casteljau algorithm)
        // First: Update previous control point
        cell_with_intersections.back().second =
            p0new + (t - t0) / (1. - t0) * (p1new - p0new);
        // Second: Add intersection and new control point
        p0new = BezierPoint(p0, p1, p2, t);
        p1new = p1new + (t - t0) / (1. - t0) * (p2 - p1new);
        const int inter_id = cell_with_intersections.size();
        sorted_list_for_mapping.push_back(std::make_pair(-p0new.x(), inter_id));
        cell_with_intersections.push_back(PtAndControl{p0new, p1new});
        nintersections++;
        // Check if intersection is entry of exit
        // If: previous vertex is below or exit, then current vertex is entry
        // Otherwise: the current vertex is exit
        vertex_nature.push_back(
            (vertex_nature.back() == 1 || vertex_nature.back() == 2) ? 3 : 2);
        t0 = t;
      }
      if (parabola_contains_vertex) break;
    }

    // If non-even # of intersections or parabola intersects with vertex, nudge!
    if (nintersections % 2 != 0 || parabola_contains_vertex) {
      datum += 10. * tol * frame[1];
      continue;
    }

    // If no intersections at all; leave loop
    if (nintersections == 0) {
      if (vertex_nature[0] == 1) {
        clipped_cell.push_back(cell_with_intersections);
      }
      break;
    }

    // Now create mappings to/from intersections
    const int nvertices = cell_with_intersections.size();
    std::sort(sorted_list_for_mapping.begin(), sorted_list_for_mapping.end());
    std::vector<int> mapping_inter_to_cell(nintersections);
    std::vector<int> mapping_cell_to_inter(nvertices);
    for (int i = 0; i < nvertices; i++) {
      mapping_cell_to_inter[i] = -1;
    }
    for (int i = 0; i < nintersections; i++) {
      mapping_inter_to_cell[i] = sorted_list_for_mapping[i].second;
      mapping_cell_to_inter[mapping_inter_to_cell[i]] = i;
    }

    ///////// Weiler-Atherton clipping algorithm
    int count_entries = 0;
    while (count_entries < nintersections / 2) {
      // Create new bezier list
      BezierList closed_clipped_region;
      // First, find entry
      int start = -1;
      for (int i = 0; i < nvertices; i++) {
        if (vertex_nature[i] == 2) {
          start = i;
          break;
        }
      }
      assert(start >= 0);
      // Add entry to new bezier list and set as "above"
      closed_clipped_region.push_back(cell_with_intersections[start]);
      vertex_nature[start] = 0;
      count_entries++;
      // Store start id
      const int check_start = start;
      // Loop over cell until next exit or start point itself;
      for (int i = 0; i < nvertices; i++) {
        // Fint next vertex on cell
        const int next_id = (start + 1 + i) % nvertices;
        // If next vertex is start, we are done
        if (next_id == check_start) {
          break;
        }
        // Add next vertex to list
        closed_clipped_region.push_back(cell_with_intersections[next_id]);
        // If next vertex is exit, we switch to parabola and find next entry
        if (vertex_nature[next_id] == 3) {
          // Find next entry in cell
          const int next_entry_id =
              mapping_inter_to_cell[mapping_cell_to_inter[next_id] + 1];
          // Set entry flag to zero
          vertex_nature[next_entry_id] = 0;
          if (next_entry_id != check_start) count_entries++;
          // Move i until before the next entry
          while ((start + 2 + i) % nvertices != next_entry_id) i++;
          // Create new control point for bezier arc extracted from parabola
          const auto p0 = cell_with_intersections[next_id].first;
          const auto p2 = cell_with_intersections[next_entry_id].first;
          closed_clipped_region.back().second =
              Vec(0.5 * (p0.x() + p2.x()),
                  p0.y() + coeff * (p0.x() - p2.x()) * p0.x());
        }
        vertex_nature[next_id] = 0;
      }
      clipped_cell.push_back(closed_clipped_region);
    }
  }

  if (parabola_contains_vertex) {
    std::cout << "WARNING: Parabola contains vertex!" << itmax
              << " nudges were not enough." << std::endl;
  }

  // Move clipped cell back to canonical frame of reference
  for (int i = 0; i < clipped_cell.size(); i++) {
    for (int j = 0; j < clipped_cell[i].size(); j++) {
      clipped_cell[i][j].first =
          parabola.datum() + frame.transpose() * clipped_cell[i][j].first;
      clipped_cell[i][j].second =
          parabola.datum() + frame.transpose() * clipped_cell[i][j].second;
    }
  }

  return clipped_cell;
}

std::vector<BezierList> ParabolaClip(
    const std::vector<BezierList>& original_cell, const Parabola& parabola) {
  // If empty cell, return empty cell
  if (original_cell.size() == 0 || parabola.isAlwaysBelow())
    return std::vector<BezierList>(0);
  if (parabola.isAlwaysAbove()) return original_cell;

  auto clipped_cell = ParabolaClip(original_cell[0], parabola);
  for (int i = 1; i < original_cell.size(); i++) {
    const auto tmp_clipped_cell = ParabolaClip(original_cell[i], parabola);
    clipped_cell.insert(clipped_cell.end(), tmp_clipped_cell.begin(),
                        tmp_clipped_cell.end());
  }
  return clipped_cell;
}

std::vector<BezierList> ClipByRectangleAndParabola(
    const BezierList& original_cell, const Vec& x0, const Vec& x1,
    const Parabola& parabola) {
  std::array<Parabola, 4> localizers;
  localizers[0] = Parabola(x1, ReferenceFrame(0.), 0.);
  localizers[1] = Parabola(x0, ReferenceFrame(M_PI / 2.), 0.);
  localizers[2] = Parabola(x0, ReferenceFrame(M_PI), 0.);
  localizers[3] = Parabola(x1, ReferenceFrame(3. * M_PI / 2.), 0.);
  auto clipped_cell = ParabolaClip(original_cell, parabola);
  for (int i = 0; i < 4; i++) {
    clipped_cell = ParabolaClip(clipped_cell, localizers[i]);
  }
  return clipped_cell;
}

double CellVolume(const BezierList& cell) {
  if (cell.size() == 0) return 0.0;

  double area = 0.0;
  for (int i = 0; i < cell.size(); i++) {
    const auto p0 = cell[i].first;
    const auto p1 = cell[i].second;
    const auto p2 = cell[(i + 1) % cell.size()].first;
    area -= (2. * p1.x() * p0.y() + p2.x() * p0.y() - 2. * p0.x() * p1.y() +
             2. * p2.x() * p1.y() - p0.x() * p2.y() - 2. * p1.x() * p2.y()) /
            6.;
  }
  return area;
}

double CellVolume(const std::vector<BezierList>& cell) {
  if (cell.size() == 0) return 0.0;

  double area = 0.0;
  for (int i = 0; i < cell.size(); i++) {
    area += CellVolume(cell[i]);
  }

  return area;
}

}  // namespace IRL2D
