#include <iostream>
#include <string>

#include "examples/pu_analysis/functions.h"
#include "examples/pu_analysis/surfaces.h"

// selecting implicit surface
void selectSurface(const std::string& implicit_surface, ImplicitF& F,
                   GradientF& gradF, HessianF& hessF) {
  if (implicit_surface == "circle") {
    F = F_circle;
    gradF = gradF_circle;
    hessF = hessF_circle;
  } else if (implicit_surface == "ellipse") {
    F = F_ellipse;
    gradF = gradF_ellipse;
    hessF = hessF_ellipse;
  } else if (implicit_surface == "sine") {
    F = F_sine;
    gradF = gradF_sine;
    hessF = hessF_sine;
  } else {
    std::cerr << "Unknown input: " << implicit_surface << "\n";
    std::exit(1);
  }
}

// setting mesh
BasicMesh setMesh(const int& Nx, const IRL2D::Vec& lower_domain,
                  const IRL2D::Vec& upper_domain) {
  const int GC = 5;
  BasicMesh mesh(Nx, Nx, GC);
  mesh.setCellBoundaries(lower_domain, upper_domain);
  return mesh;
}

// projecting point onto implicit surface
Eigen::Vector2d project_onto_surface(const Eigen::Vector2d& x0, ImplicitF F,
                                     GradientF gradF) {
  Eigen::Vector2d x_proj = x0;
  const int max_iter = 50;
  const double tol = 1e-10;
  for (int i = 0; i < max_iter; i++) {
    double f = F(x_proj.x(), x_proj.y());
    Eigen::Vector2d g = gradF(x_proj.x(), x_proj.y());
    double g_norm2 = g.squaredNorm();
    if (g_norm2 < 1e-14) break;
    x_proj -= (f / g_norm2) * g;
    if (std::abs(f) < tol) break;
  }

  return x_proj;
}

// build parabola
IRL2D::Parabola build_parabola(const Eigen::Vector2d& x0, GradientF gradF,
                               HessianF hessF) {
  Eigen::Vector2d g = gradF(x0.x(), x0.y());
  IRL2D::Vec normal = {g.x(), g.y()};
  normal.normalize();
  IRL2D::Vec tangent = {normal.y(), -normal.x()};
  Eigen::Matrix2d h = hessF(x0.x(), x0.y());
  double curvature = (g(0) * g(0) * h(1, 1) - 2.0 * g(0) * g(1) * h(0, 1) +
                      g(1) * g(1) * h(0, 0)) /
                     (std::pow(g(0) * g(0) + g(1) * g(1), 1.5));

  IRL2D::Parabola parabola;
  parabola.datum() = IRL2D::Vec(x0.x(), x0.y());
  parabola.frame() = {tangent, normal};
  parabola.coeff() = 0.5 * curvature;
  return parabola;
}

// Finding F with caching
double evaluate_or_cache(double x, double y, ImplicitF F,
                         std::unordered_map<Point, double, PointHash>& cache) {
  Point pt = {x, y};
  auto it = cache.find(pt);
  if (it != cache.end()) {
    return it->second;
  }
  double val = F(x, y);
  cache[pt] = val;
  return val;
}

// mixed cell?
CellStatus get_cell_status(
    const Cell& cell, ImplicitF F, GradientF gradF, HessianF hessF,
    std::unordered_map<Point, double, PointHash>& F_cache) {
  // target cell points
  double x = cell.x, y = cell.y, dx = cell.dx;
  std::vector<Point> pts = {{x, y},
                            {x + dx, y},
                            {x, y + dx},
                            {x + dx, y + dx},
                            {x + dx / 2, y},
                            {x, y + dx / 2},
                            {x + dx, y + dx / 2},
                            {x + dx / 2, y + dx}};
  std::vector<double> Fvals;
  for (const auto& pt : pts) {
    Fvals.push_back(evaluate_or_cache(pt.first, pt.second, F, F_cache));
  }

  bool all_pos = std::all_of(Fvals.begin(), Fvals.end(),
                             [](double Fval) { return Fval > 0; });
  bool all_neg = std::all_of(Fvals.begin(), Fvals.end(),
                             [](double Fval) { return Fval < 0; });
  if (!all_pos && !all_neg) return CellStatus::Mixed;

  // 3x3 stencil
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      double cx = x + i * dx;
      double cy = y + j * dx;
      std::vector<Point> stencil_pts = {{cx, cy},
                                        {cx + dx, cy},
                                        {cx, cy + dx},
                                        {cx + dx, cy + dx},
                                        {cx + dx / 2, cy},
                                        {cx, cy + dx / 2},
                                        {cx + dx, cy + dx / 2},
                                        {cx + dx / 2, cy + dx}};
      for (const auto& pt : stencil_pts) {
        double val = evaluate_or_cache(pt.first, pt.second, F, F_cache);
        if ((val > 0 && !all_pos) || (val < 0 && !all_neg)) {  // mismatch?
          Eigen::Vector2d x_center(x + 0.5 * dx, y + 0.5 * dx);
          Eigen::Vector2d x_proj = project_onto_surface(x_center, F, gradF);
          IRL2D::Parabola parabola = build_parabola(x_proj, gradF, hessF);
          IRL2D::BezierList rectangle = IRL2D::RectangleFromBounds(
              IRL2D::Vec(x, y), IRL2D::Vec(x + dx, y + dx));
          IRL2D::BezierList clipped_interface =
              IRL2D::ParabolaClip(rectangle, parabola, true);
          return clipped_interface.empty()
                     ? (all_pos ? CellStatus::Above : CellStatus::Below)
                     : CellStatus::Mixed;
        }
      }
    }
  }
  return all_pos ? CellStatus::Above : CellStatus::Below;
}

// refine cell
void refine_cell(std::unique_ptr<Cell>& cell, ImplicitF F, const int& max_level,
                 GradientF gradF, HessianF hessF,
                 std::unordered_map<Point, double, PointHash>& F_cache) {
  cell->status = get_cell_status(*cell, F, gradF, hessF, F_cache);
  if (cell->status != CellStatus::Mixed || cell->level >= max_level) return;

  double x = cell->x, y = cell->y, dx = cell->dx / 2.0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      auto child =
          std::make_unique<Cell>(x + i * dx, y + j * dx, dx, cell->level + 1);
      refine_cell(child, F, max_level, gradF, hessF, F_cache);
      cell->children.push_back(std::move(child));
    }
  }
}

// refining all cells
std::vector<std::unique_ptr<Cell>> refine_grid(
    const BasicMesh& mesh, ImplicitF F, const int& max_level, GradientF gradF,
    HessianF hessF, std::unordered_map<Point, double, PointHash>& F_cache) {
  std::vector<std::unique_ptr<Cell>> grid;

  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      double x = mesh.x(i);
      double y = mesh.y(j);
      double dx = mesh.dx();
      auto cell = std::make_unique<Cell>(x, y, dx, 0);
      refine_cell(cell, F, max_level, gradF, hessF, F_cache);
      grid.push_back(std::move(cell));
    }
  }
  return grid;
}

// area above and below interface in refined cell
std::pair<double, double> get_mixed_area(const double& x, const double& y,
                                         const double& dx, ImplicitF F,
                                         GradientF gradF, HessianF hessF) {
  // projecting cell center onto implicit surface
  Eigen::Vector2d x_center(x + 0.5 * dx, y + 0.5 * dx);
  Eigen::Vector2d x_proj = project_onto_surface(x_center, F, gradF);

  // generate parabola
  IRL2D::Parabola parabola = build_parabola(x_proj, gradF, hessF);

  // computing areas
  IRL2D::BezierList cell =
      IRL2D::RectangleFromBounds(IRL2D::Vec(x, y), IRL2D::Vec(x + dx, y + dx));
  IRL2D::Moments moments = IRL2D::ComputeMoments(cell, parabola);
  double area_inside = moments.m0();
  double area_outside = (dx * dx) - area_inside;

  return {area_inside, area_outside};
}

// traverse through tree to get area
void get_cell_area(const Cell& cell, AreaResult& result, ImplicitF F,
                   GradientF gradF, HessianF hessF) {
  if (!cell.children.empty()) {
    for (const auto& child : cell.children) {
      get_cell_area(*child, result, F, gradF, hessF);
    }
    return;
  }

  double area = cell.area();
  switch (cell.status) {
    case CellStatus::Below:
      result.inside += area;
      break;
    case CellStatus::Above:
      result.outside += area;
      break;
    case CellStatus::Mixed:
      std::pair<double, double> mixed_area =
          get_mixed_area(cell.x, cell.y, cell.dx, F, gradF, hessF);
      result.inside += mixed_area.first;
      result.outside += mixed_area.second;
      break;
  }
}

// getting leaf cells
void collect_leaf_cells(const Cell& cell, std::vector<const Cell*>& output,
                        std::vector<const Cell*>& output_mixed) {
  if (cell.children.empty()) {
    output.push_back(&cell);
    if (cell.status == CellStatus::Mixed) {
      output_mixed.push_back(&cell);
    }
    return;
  }
  for (const auto& child : cell.children) {
    collect_leaf_cells(*child, output, output_mixed);
  }
}
