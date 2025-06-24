#include <iostream>
#include <string>

#include "examples/initialize_vf/functions.h"

// selecting implicit surface
void selectSurface(const std::string& implicit_surface, ImplicitF& F,
                   GradientF& gradF, HessianF& hessF){

  if (implicit_surface == "circle"){
    F = F_circle;
    gradF = gradF_circle;
    hessF = hessF_circle;
  } else if (implicit_surface == "ellipse"){
    F = F_ellipse;
    gradF = gradF_ellipse;
    hessF = hessF_ellipse;
  } else {
    std::cerr << "Unknown input: " << implicit_surface << "\n";
    std::exit(1);
  }

}

// setting mesh
BasicMesh setMesh(const int& Nx, const IRL2D::Vec& lower_domain,
                  const IRL2D::Vec& upper_domain){  
  const int GC = 3;
  BasicMesh mesh(Nx, Nx, GC);
  mesh.setCellBoundaries(lower_domain, upper_domain);
  return mesh;
}

// mixed cell?
CellStatus get_cell_status(const Cell& cell, ImplicitF F){

  double x = cell.x, y = cell.y, dx = cell.dx;

  std::vector<double> Fvals = {
    F(x, y), F(x + dx, y), F(x, y + dx), F(x + dx, y + dx)
  };

  bool all_pos = std::all_of(Fvals.begin(), Fvals.end(), [](double Fval){return Fval > 0;});  
  bool all_neg = std::all_of(Fvals.begin(), Fvals.end(), [](double Fval){return Fval < 0;});

  if (all_pos) return CellStatus::Above;
  if (all_neg) return CellStatus::Below;
  return CellStatus::Mixed;
}

// refine cell
void refine_cell(std::unique_ptr<Cell>& cell, ImplicitF F, const int& max_level){

  cell->status = get_cell_status(*cell, F);
  if (cell->status != CellStatus::Mixed || cell->level >= max_level) return;

  double x = cell->x, y = cell->y, dx = cell->dx/2.0;
  for (int i = 0; i < 2; i++){
    for (int j = 0; j < 2; j++){
      auto child = std::make_unique<Cell>(x + i*dx, y + j*dx, dx, cell->level + 1);
            refine_cell(child, F, max_level);
            cell->children.push_back(std::move(child));
    }
  }

}

// refining all cells
std::vector<std::unique_ptr<Cell>> refine_grid(const BasicMesh& mesh, ImplicitF F,
                                               const int& max_level){
                                                
  std::vector<std::unique_ptr<Cell>> grid;

  for (int i = mesh.imin(); i <= mesh.imax(); i++){
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++){
      double x = mesh.x(i);
      double y = mesh.y(j);
      double dx = mesh.dx();
      auto cell = std::make_unique<Cell>(x, y, dx, 0);
      refine_cell(cell, F, max_level);
      grid.push_back(std::move(cell));
    }
  }

  return grid;
}

// area above and below interface in refined cell
std::pair<double, double> get_mixed_area(const double& x, const double& y,
                                         const double& dx, ImplicitF F,
                                         GradientF gradF, HessianF hessF){
  
  // projecting cell center onto surface
  Eigen::Vector2d x_center(x + 0.5*dx, y + 0.5*dx);
  Eigen::Vector2d x_proj = x_center;
  const int max_iter = 50;
  const double tol = 1e-10;
  for (int i = 0; i < max_iter; i++){
    double f = F(x_proj.x(), x_proj.y());
    Eigen::Vector2d g = gradF(x_proj.x(), x_proj.y());
    double g_norm2 = g.squaredNorm();
    if (g_norm2 < 1e-14) break;
    Eigen::Vector2d dx = (f / g_norm2) * g;
    x_proj -= dx;
    if (std::abs(f) < tol) break;
  }

  // generate parabola
  IRL2D::Parabola parabola;
  parabola.datum() = IRL2D::Vec(x_proj.x(), x_proj.y());
  Eigen::Vector2d g = gradF(x_proj.x(), x_proj.y());
  IRL2D::Vec normal = {g.x(), g.y()};
  normal.normalize();
  IRL2D::Vec tangent = {normal.y(), -normal.x()};
  parabola.frame() = {tangent, normal};
  Eigen::Matrix2d h = hessF(x_proj.x(), x_proj.y());
  double curvature = (g(0) * g(0) * h(1,1) - 2.0 * g(0) * g(1) * h(0,1) + g(1) * g(1) * h(0,0)) / 
                     (std::pow(g(0) * g(0) + g(1) * g(1), 1.5));
  parabola.coeff() = 0.5 * curvature;
  //std::cout << parabola.coeff() << std::endl;

  // computing areas
  IRL2D::BezierList cell = IRL2D::RectangleFromBounds( IRL2D::Vec(x, y), 
                                                       IRL2D::Vec(x + dx, y + dx) );
  IRL2D::Moments moments = IRL2D::ComputeMoments(cell, parabola);
  double area_inside = moments.m0();
  double area_outside = (dx * dx) - area_inside;

  return {area_inside, area_outside};

}

// traverse through tree to get area
void get_cell_area(const Cell& cell, AreaResult& result, ImplicitF F,
                   GradientF gradF, HessianF hessF){
  
  if (!cell.children.empty()){
    for (const auto& child : cell.children){
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
      std::pair<double, double> mixed_area = get_mixed_area(cell.x, cell.y, cell.dx,
                                                            F, gradF, hessF);
      result.inside += mixed_area.first;
      result.outside += mixed_area.second;
      break;
  }
}

