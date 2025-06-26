#include <iostream>
#include <string>

#include "examples/initialize_vf/functions.h"

// selecting implicit surface
void selectSurface(const std::string& implicit_surface, ImplicitF& F,
                   GradientF& gradF, HessianF& hessF){

  if (implicit_surface == "sphere"){
    F = F_sphere;
    gradF = gradF_sphere;
    hessF = hessF_sphere;
  } else if (implicit_surface == "ellipsoid"){
    F = F_ellipsoid;
    gradF = gradF_ellipsoid;
    hessF = hessF_ellipsoid;
  } else {
    std::cerr << "Unknown input: " << implicit_surface << "\n";
    std::exit(1);
  }

}

// projecting point onto implicit surface
Eigen::Vector3d project_onto_surface(const Eigen::Vector3d& x0, ImplicitF F,
                                     GradientF gradF){
  
  Eigen::Vector3d x_proj = x0;
  const int max_iter = 100;
  const double tol = 1e-10;
  for (int i = 0; i < max_iter; i++){
    double f = F(x_proj(0), x_proj(1), x_proj(2));
    Eigen::Vector3d g = gradF(x_proj(0), x_proj(1), x_proj(2));
    double g_norm2 = g.squaredNorm();
    if (g_norm2 < 1e-14) break;
    x_proj -= (f / g_norm2) * g;
    if (std::abs(f) < tol) break;
  }

  return x_proj;
}

// generate paraboloid
// IRL::Paraboloid generateLocalParaboloid(const IRL::Pt& a_pt,
//                                         const Eigen::Vector3d& a_gradF,
//                                         const Eigen::Matrix3d& a_hessF) {
//   const Eigen::Matrix3d hessF = 0.5 * (a_hessF + a_hessF.transpose());
//   // This uses the method described in
//   // https://www.geometrictools.com/Documentation/PrincipalCurvature.pdf
//   const double inv_gradF_norm = 1. / a_gradF.norm();
//   const IRL::Normal normal =
//       IRL::Normal(a_gradF(0), a_gradF(1), a_gradF(2)) * inv_gradF_norm;
//   IRL::ReferenceFrame frame = referenceFrameFromNormal(normal);
//   Eigen::MatrixXd J(3, 2);
//   for (int i = 0; i < 3; i++) {
//     for (int j = 0; j < 2; j++) {
//       J(i, j) = frame[j][i];
//     }
//   }
//   const Eigen::Matrix2d A = (J.transpose() * hessF * J) * inv_gradF_norm;
//   Eigen::EigenSolver<Eigen::Matrix2d> eigensolver(A);
//   const double eval1 = eigensolver.eigenvalues()(0).real();
//   const double eval2 = eigensolver.eigenvalues()(1).real();
//   const Eigen::Vector2d evec1 =
//       Eigen::Vector2d(eigensolver.eigenvectors()(0, 0).real(),
//                       eigensolver.eigenvectors()(1, 0).real());
//   const Eigen::Vector3d T1 = J * evec1;
//   frame[0] = IRL::Normal(T1(0), T1(1), T1(2));
//   frame[0].normalize();
//   frame[1] = IRL::crossProduct(frame[2], frame[0]);
//   return IRL::Paraboloid(a_pt, frame, 0.5 * eval1, 0.5 * eval2);
// }

// Finding F with caching
double evaluate_or_cache(const double& x, const double& y, const double& z, ImplicitF F,
                         std::unordered_map<Eigen::Vector3d, double, Vector3dHash>& cache) {
  Eigen::Vector3d pt(x, y, z);
  auto it = cache.find(pt);
  if (it != cache.end()) {
      return it->second;
  }
  double val = F(pt[0], pt[1], pt[2]);
  cache[pt] = val;
  return val;
}

// sample points
std::vector<Eigen::Vector3d> get_sample_points(const double& x, const double& y,
                                               const double& z, const double& dx) {

  std::vector<Eigen::Vector3d> pts;

  // corners
  for (int i = 0; i <= 1; ++i)
    for (int j = 0; j <= 1; ++j)
      for (int k = 0; k <= 1; ++k)
        pts.emplace_back(x + i * dx, y + j * dx, z + k * dx);

  // edge centers
  pts.emplace_back(x + dx / 2, y, z);
  pts.emplace_back(x + dx / 2, y + dx, z);
  pts.emplace_back(x, y + dx / 2, z);
  pts.emplace_back(x + dx, y + dx / 2, z);
  pts.emplace_back(x, y, z + dx / 2);
  pts.emplace_back(x + dx, y, z + dx / 2);

  // face centers 
  pts.emplace_back(x + dx / 2, y + dx / 2, z);
  pts.emplace_back(x + dx / 2, y, z + dx / 2);
  pts.emplace_back(x, y + dx / 2, z + dx / 2);
  pts.emplace_back(x + dx, y + dx / 2, z + dx / 2);
  pts.emplace_back(x + dx / 2, y + dx, z + dx / 2);
  pts.emplace_back(x + dx / 2, y + dx / 2, z + dx);

  return pts;
}


// mixed cell?
CellStatus get_cell_status(const Cell& cell, ImplicitF F,
                           GradientF gradF, HessianF hessF,
                           std::unordered_map<Eigen::Vector3d, double, Vector3dHash>& F_cache){
  
  // target cell points
  double x = cell.x, y = cell.y, z = cell.z, dx = cell.dx;

  std::vector<Eigen::Vector3d> pts = get_sample_points(x, y, z, dx);

  std::vector<double> Fvals;
  for (const auto& pt : pts){
    Fvals.push_back(evaluate_or_cache(pt[0], pt[1], pt[2], F, F_cache));
  }
  
  bool all_pos = std::all_of(Fvals.begin(), Fvals.end(), [](double Fval){return Fval > 0;});  
  bool all_neg = std::all_of(Fvals.begin(), Fvals.end(), [](double Fval){return Fval < 0;});
  if (!all_pos && !all_neg) return CellStatus::Mixed;

  // 3x3x3 stencil
  // for (int i = -1; i <= 1; i++){
  //   for (int j = -1; j <= 1; j++){
  //     for (int k = -1; k <= 1; k++){
  //       double cx = x + i * dx;
  //       double cy = y + j * dx;
  //       double cz = z + k * dx;
  //       std::vector<Eigen::Vector3d> stencil_pts = get_sample_points(cx, cy, cz, dx);
  //       for (const auto& pt : stencil_pts){
  //         double val = evaluate_or_cache(pt[0], pt[1], pt[2], F, F_cache);
  //         if ((val > 0 && !all_pos) || (val < 0 && !all_neg)){
  //           // check clipping
  //           Eigen::Vector3d x_center(x + 0.5*dx, y + 0.5*dx, z + 0.5*dx);
  //           Eigen::Vector3d x_proj = project_onto_surface(x_center, F, gradF);
  //           IRL::Pt x_proj_pt(x_proj(0), x_proj(1), x_proj(2));
  //           IRL::Paraboloid paraboloid  = generateLocalParaboloid(x_proj_pt, gradF, hessF);
  //           // IRL2D::BezierList rectangle = IRL2D::RectangleFromBounds( IRL2D::Vec(x, y), 
  //           //                                                           IRL2D::Vec(x + dx, y + dx) );
  //           // IRL2D::BezierList clipped_interface = IRL2D::ParabolaClip(rectangle, parabola, true);
  //           // return clipped_interface.empty() ? (all_pos ? CellStatus::Above : CellStatus::Below) : CellStatus::Mixed;
  //         }
  //       }
  //     }
  //   }
  // } 

  return all_pos ? CellStatus::Above : CellStatus::Below;
}

// refine cell
void refine_cell(std::unique_ptr<Cell>& cell, ImplicitF F,
                 GradientF gradF, HessianF hessF, const int& max_level,
                 std::unordered_map<Eigen::Vector3d, double, Vector3dHash>& F_cache){

  cell->status = get_cell_status(*cell, F, gradF, hessF, F_cache);
  if (cell->status != CellStatus::Mixed || cell->level >= max_level) return;

  double x = cell->x, y = cell->y, z = cell->z, dx = cell->dx/2.0;

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        auto child = std::make_unique<Cell>(
            x + i * dx,
            y + j * dx,
            z + k * dx,
            dx,
            cell->level + 1
        );
        refine_cell(child, F, gradF, hessF, max_level, F_cache);
        cell->children.push_back(std::move(child));
      }
    }
  }

}


// refining all cells
std::vector<std::unique_ptr<Cell>> refine_grid(const BasicMesh& mesh, ImplicitF F,
                                               GradientF gradF, HessianF hessF, const int& max_level,
                                               std::unordered_map<Eigen::Vector3d, double, Vector3dHash>& F_cache){
                                                
  std::vector<std::unique_ptr<Cell>> grid;

  for (int i = mesh.imin(); i <= mesh.imax(); i++){
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++){
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++){
        double x = mesh.x(i), y = mesh.y(j), z = mesh.z(k);
        double dx = mesh.dx();
        auto cell = std::make_unique<Cell>(x, y, z, dx, 0);
        refine_cell(cell, F, gradF, hessF, max_level, F_cache);
        grid.push_back(std::move(cell));
      }
      
    }
  }
  return grid;
}

// area above and below interface in refined cell
std::pair<double, double> get_mixed_area(const Eigen::Vector3d x0,
                                         const double& dx, ImplicitF F,
                                         GradientF gradF, HessianF hessF){
  
  // projecting cell center onto implicit surface
  Eigen::Vector3d x_center(x0(0) + 0.5*dx, x0(1) + 0.5*dx, x0(2) + 0.5 * dx);
  Eigen::Vector3d x_proj = project_onto_surface(x_center, F, gradF);

  // // generate parabola
  // IRL2D::Parabola parabola = build_parabola(x_proj, gradF, hessF);

  // // computing areas
  // IRL2D::BezierList cell = IRL2D::RectangleFromBounds( IRL2D::Vec(x, y), 
  //                                                      IRL2D::Vec(x + dx, y + dx) );
  // IRL2D::Moments moments = IRL2D::ComputeMoments(cell, parabola);
  // double area_inside = moments.m0();
  // double area_outside = (dx * dx) - area_inside;

  double area_inside = 0.0;
  double area_outside = 0.0;

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
      Eigen::Vector3d x0(cell.x, cell.y, cell.z);
      std::pair<double, double> mixed_area = get_mixed_area(x0, cell.dx,
                                                            F, gradF, hessF);
      result.inside += mixed_area.first;
      result.outside += mixed_area.second;
      break;
  }
}

// getting leaf cells
void collect_leaf_cells(const Cell& cell, std::vector<const Cell*>& output,
                        std::vector<const Cell*>& output_mixed){
  if (cell.children.empty()){
    output.push_back(&cell);
    if (cell.status == CellStatus::Mixed){
      output_mixed.push_back(&cell);
    }
    return;
  }
  for (const auto& child : cell.children){
    collect_leaf_cells(*child, output, output_mixed);
  }
}