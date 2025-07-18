#include <iostream>
#include <string>

#include "examples/initialize_vf/functions.h"

// select implicit surface
std::unique_ptr<Surface> makeSurface(const std::string& name) {
  if (name == "sphere") return std::make_unique<Sphere>();
  if (name == "ellipsoid") return std::make_unique<Ellipsoid>();
  if (name == "genus") return std::make_unique<Genus>();
  if (name == "orthocircle") return std::make_unique<Orthocircle>();
  std::cerr << "Unknown input: " << name << "\n";
  std::exit(1);
}

// initialize
void initialize(std::unique_ptr<Surface>& surface,
                ImplicitF& F, GradientF& gradF,
                HessianF& hessF, BasicMesh& mesh) {
  F = [surface_ptr = surface.get()](const double& x, const double& y, const double& z) {
    return surface_ptr->F(x, y, z);
  };
  gradF = [surface_ptr = surface.get()](const double& x, const double& y, const double& z) {
    return surface_ptr->gradF(x, y, z);
  };
  hessF = [surface_ptr = surface.get()](const double& x, const double& y, const double& z) {
    return surface_ptr->hessF(x, y, z);
  };

  auto [lower, upper] = surface->domainBounds();
  mesh.setCellBoundaries(lower, upper);
}

// projecting point onto implicit surface
Eigen::Vector3d project_onto_surface(const Eigen::Vector3d& x0, ImplicitF F,
                                     GradientF gradF){
  
  Eigen::Vector3d x_proj = x0;
  const int max_iter = 500;
  const double tol = 1e-10;
  for (int i = 0; i < max_iter; i++){
    double f = F(x_proj(0), x_proj(1), x_proj(2));
    Eigen::Vector3d g = gradF(x_proj(0), x_proj(1), x_proj(2));
    double g_norm2 = g.squaredNorm();
    if (g_norm2 < 1e-14) break;
    x_proj -= (f / g_norm2) * g;
    if (std::abs(f) < tol) break;
    if (i == (max_iter - 1)){
      std::cout << "Max iterations reached. " << "f = " << std::abs(f) << std::endl;
    }
  }
  return x_proj;
}

// Finding F with caching
double evaluate_or_cache(const Eigen::Vector3d& pt, ImplicitF F,
                         std::unordered_map<Eigen::Vector3d, double, Vector3dHash>& cache) {
  auto it = cache.find(pt);
  if (it != cache.end()) {
      return it->second;
  }
  double val = F(pt[0], pt[1], pt[2]);
  cache[pt] = val;
  return val;
}

// sample points
// std::vector<Eigen::Vector3d> get_sample_points(const double& x, const double& y,
//                                                const double& z, const double& dx) {
//   std::vector<Eigen::Vector3d> pts;
//   // corners
//   for (int i = 0; i <= 1; ++i)
//     for (int j = 0; j <= 1; ++j)
//       for (int k = 0; k <= 1; ++k)
//         pts.emplace_back(x + i * dx, y + j * dx, z + k * dx);
//   // edge centers
//   // pts.emplace_back(x + dx / 2, y, z);
//   // pts.emplace_back(x + dx / 2, y + dx, z);
//   // pts.emplace_back(x, y + dx / 2, z);
//   // pts.emplace_back(x + dx, y + dx / 2, z);
//   // pts.emplace_back(x, y, z + dx / 2);
//   // pts.emplace_back(x + dx, y, z + dx / 2);
//   // pts.emplace_back(x + dx / 2, y, z + dx);
//   // pts.emplace_back(x + dx / 2, y + dx, z + dx);
//   // pts.emplace_back(x, y + dx / 2, z + dx);
//   // pts.emplace_back(x + dx/2, y + dx / 2, z + dx);
//   // pts.emplace_back(x, y + dx, z + dx / 2);
//   // pts.emplace_back(x + dx, y + dx, z + dx / 2);
//   // face centers 
//   pts.emplace_back(x + dx / 2, y + dx / 2, z);
//   pts.emplace_back(x + dx / 2, y, z + dx / 2);
//   pts.emplace_back(x, y + dx / 2, z + dx / 2);
//   pts.emplace_back(x + dx, y + dx / 2, z + dx / 2);
//   pts.emplace_back(x + dx / 2, y + dx, z + dx / 2);
//   pts.emplace_back(x + dx / 2, y + dx / 2, z + dx);
//   return pts;
// }

std::vector<Eigen::Vector3d> get_sample_points(const double& x, const double& y,
                                               const double& z, const double& dx,
                                               const bool& use_stencil) {
  std::vector<Eigen::Vector3d> pts;
  double sx, sy, sz, h;
  if (use_stencil) {
    sx = x - dx;
    sy = y - dx;
    sz = z - dx;
    h = 3.0 * dx;
  } else {
    sx = x;
    sy = y;
    sz = z;
    h = dx;
  }

  // corners
  for (int i = 0; i <= 1; ++i)
    for (int j = 0; j <= 1; ++j)
      for (int k = 0; k <= 1; ++k)
        pts.emplace_back(sx + i * h, sy + j * h, sz + k * h);

  // face centers
  pts.emplace_back(sx + h / 2, sy + h / 2, sz);
  pts.emplace_back(sx + h / 2, sy, sz + h / 2);
  pts.emplace_back(sx, sy + h / 2, sz + h / 2);
  pts.emplace_back(sx + h, sy + h / 2, sz + h / 2);
  pts.emplace_back(sx + h / 2, sy + h, sz + h / 2);
  pts.emplace_back(sx + h / 2, sy + h / 2, sz + h);

  return pts;
}


// mixed cell?
CellStatus get_cell_status(const Cell& cell, ImplicitF F,
                           GradientF gradF, HessianF hessF,
                           std::unordered_map<Eigen::Vector3d, double, Vector3dHash>& F_cache){
  
  // target cell points
  double x = cell.x, y = cell.y, z = cell.z, dx = cell.dx;

  std::vector<Eigen::Vector3d> pts = get_sample_points(x, y, z, dx, false);

  // std::vector<double> Fvals;
  // for (const auto& pt : pts){
  //   Fvals.push_back(evaluate_or_cache(pt, F, F_cache));
  // }
  // bool all_pos = std::all_of(Fvals.begin(), Fvals.end(), [](double Fval){return Fval > 0;});  
  // bool all_neg = std::all_of(Fvals.begin(), Fvals.end(), [](double Fval){return Fval < 0;});
  // if (!all_pos && !all_neg) return CellStatus::Mixed;

   bool all_pos = false, all_neg = false;
   for (const auto& pt : pts){
    double Fval = evaluate_or_cache(pt, F, F_cache);
    if (Fval > 0) all_pos = true;
    else if  (Fval < 0) all_neg = true;
    if (all_pos && all_neg) return CellStatus::Mixed;
   }

  // 3x3x3 stencil
  std::vector<Eigen::Vector3d> stencil_pts = get_sample_points(x, y, z, dx, true);
  for (const auto& pt : stencil_pts){
    double val = evaluate_or_cache(pt, F, F_cache);
    if ((val > 0 && !all_pos) || (val < 0 && !all_neg)){
      Eigen::Vector3d x_center(x + 0.5*dx, y + 0.5*dx, z + 0.5*dx);
      Eigen::Vector3d x_proj = project_onto_surface(x_center, F, gradF);
      IRL::Pt x_proj_pt(x_proj(0), x_proj(1), x_proj(2));
      Eigen::Vector3d gradF_val = gradF(x_proj(0), x_proj(1), x_proj(2));
      Eigen::Matrix3d hessF_val = hessF(x_proj(0), x_proj(1), x_proj(2));
      IRL::Paraboloid paraboloid = IRL::Paraboloid::fromDerivatives(x_proj_pt, gradF_val, hessF_val);
      auto rectangle = IRL::RectangularCuboid::fromBoundingPts(IRL::Pt(x, y, z), 
                                                               IRL::Pt(x + dx, y + dx, z + dx));
      auto moments = IRL::getNormalizedVolumeMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(rectangle, paraboloid);
      // auto moments = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(rectangle, paraboloid);
      double liq_vf = moments[0].volume();
      bool not_mixed = true;
      if (liq_vf >= IRL::global_constants::VF_LOW && liq_vf <= IRL::global_constants::VF_HIGH) not_mixed = false;
      return not_mixed ? (all_pos ? CellStatus::Above : CellStatus::Below) : CellStatus::Mixed;
    }
  }
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

  // generate parabola
  IRL::Pt x_proj_pt(x_proj(0), x_proj(1), x_proj(2));
  Eigen::Vector3d gradF_val = gradF(x_proj(0), x_proj(1), x_proj(2));
  Eigen::Matrix3d hessF_val = hessF(x_proj(0), x_proj(1), x_proj(2));
  IRL::Paraboloid paraboloid = IRL::Paraboloid::fromDerivatives(x_proj_pt, gradF_val, hessF_val);

  // computing areas
  auto cell = IRL::RectangularCuboid::fromBoundingPts(IRL::Pt(x0(0), x0(1), x0(2)), 
                                                      IRL::Pt(x0(0) + dx, x0(1) + dx, x0(2) + dx));
  auto moments = IRL::getNormalizedVolumeMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(cell, paraboloid);
  double area_inside = moments[0].volume();
  double area_outside = moments[1].volume();

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

// vf initializer class definitions
vfInitializer::vfInitializer(const int& Nx, const int& max_refine_level,
                                     const std::string& surface_name,
                                     const std::string& vtk_output_path)
    : Nx_(Nx),
      max_refine_level_(max_refine_level),
      surface_name_(surface_name),
      vtk_output_path_(vtk_output_path),
      mesh_(Nx, Nx, Nx, 5),
      vf_(&mesh_)  {
    surface_ = makeSurface(surface_name_);
}

void vfInitializer::run() {
    ImplicitF F;
    GradientF gradF;
    HessianF hessF;
    initialize(surface_, F, gradF, hessF, mesh_);

    auto start = std::chrono::high_resolution_clock::now();

    std::unordered_map<Eigen::Vector3d, double, Vector3dHash> F_cache;
    grid_ = refine_grid(mesh_, F, gradF, hessF, max_refine_level_, F_cache);

    int idx = 0;
    for (int i = mesh_.imin(); i <= mesh_.imax(); i++) {
      for (int j = mesh_.jmin(); j <= mesh_.jmax(); j++) {
        for (int k = mesh_.kmin(); k <= mesh_.kmax(); k++) {
          AreaResult result;
          get_cell_area(*grid_[idx], result, F, gradF, hessF);
          vf_(i,j,k) = result.volume_fraction();
          volume_amr_ += result.inside;
          ++idx;
        }
      }
    }

    auto end = std::chrono::high_resolution_clock::now();
    elapsed_ = end - start;

    for (const auto& root : grid_) {
        collect_leaf_cells(*root, leaf_cells_, mixed_leaf_cells_);
    }

    write_vtr(vtk_output_path_ + "vf.vtr", vf_, mesh_);
    write_vtu(vtk_output_path_ + "mixed_amr.vtu", mixed_leaf_cells_);
}

// triangle class
Triangle::Triangle(const Eigen::Vector3d& v0,
                   const Eigen::Vector3d& v1,
                   const Eigen::Vector3d& v2)
    : v0_(v0), v1_(v1), v2_(v2) {}

double Triangle::computeM0() const {
  return 0.5 * ((v1_ - v0_).cross(v2_ - v0_)).norm();
}

Eigen::Vector3d Triangle::computeM1() const {
  double A = computeM0();
  return A * (v0_ + v1_ + v2_) / 3.0;
}

Eigen::Matrix3d Triangle::computeM2() const {
  double A = computeM0();
  Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
  const std::vector<Eigen::Vector3d> v = {v0_, v1_, v2_};
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      M += v[i] * v[j].transpose();
    }
  }
  return A / 12.0 * M;
}

// surface metric class function definitions
SurfaceMetrics::SurfaceMetrics(const Cell& base_cell, ImplicitF F, const int num_divisions)
    : base_cell_(base_cell),
      F_([F](const IRL::Pt& p) { return F(p[0], p[1], p[2]); }),
      num_divisions_(num_divisions) {}

std::vector<const Cell*> SurfaceMetrics::collectMixedLeafCells() const {
  std::vector<const Cell*> mixed_leaf_cells;
  std::function<void(const Cell&)> recurse = [&](const Cell& c){
    if (c.children.empty()){
        if (c.status == CellStatus::Mixed){
          mixed_leaf_cells.push_back(&c);
        }
    } else {
        for (const auto& child : c.children){
          recurse(*child);
        }
    }
  };
  recurse(base_cell_);
  return mixed_leaf_cells;
}

// metric for surface within a base cell
Metrics SurfaceMetrics::computeSurfaceMetrics() {
  
  auto mixed_leaves = collectMixedLeafCells();

  // loop over all leaves
  for (const auto& c : mixed_leaves){
    auto domain = IRL::RectangularCuboid::fromBoundingPts(IRL::Pt(c->x, c->y, c->z), 
                                                          IRL::Pt(c->x + c->dx, c->y + c->dx, c->z + c->dx));
    IRL::MarchingCubes mc(domain, F_);
    auto surface = mc.triangulate(num_divisions_);
    auto triangles = surface.getTriangleList();

    // looping over all triangles within a leaf
    for (const auto& tri : triangles){
      Eigen::Vector3d v0(tri[0][0], tri[0][1], tri[0][2]);
      Eigen::Vector3d v1(tri[1][0], tri[1][1], tri[1][2]);
      Eigen::Vector3d v2(tri[2][0], tri[2][1], tri[2][2]);
      Triangle triangle(v0, v1, v2);
      metrics_.sM0 += triangle.computeM0();
      metrics_.sM1 += triangle.computeM1();
      metrics_.sM2 += triangle.computeM2();
    } 
  }

  return metrics_;
}


// for surface moment calculations
// double peskin_delta(const double& phi, const double& epsilon) {
//   if (std::abs(phi) > epsilon) return 0.0;
//   return 0.5 / epsilon * (1.0 + std::cos(M_PI * phi / epsilon));
// }
// double gaussian_delta_truncated(double phi, double sigma) {
//   constexpr double sqrt2pi = std::sqrt(2.0 * M_PI);
//   if (std::abs(phi) > 3.0 * sigma)
//       return 0.0;
//   return (1.0 / (sqrt2pi * sigma)) * std::exp(- (phi * phi) / (2.0 * sigma * sigma));
// }
// // within each child mixed cell
// SurfaceMoments compute_surface_moments(const Cell& cell, const int& N, 
//                                        ImplicitF F, GradientF gradF){
//   SurfaceMoments moments;
//   double x0 = cell.x, y0 = cell.y, z0 = cell.z;
//   double eps = 1.0 * cell.dx / N; // for peskin delta
//   double dv = std::pow(cell.dx / N, 3.0);
//   for (int i = 0; i < N; i++){
//     double x = x0 + (i + 0.5) * cell.dx / N;
//     for (int j = 0; j < N; j++){
//       double y = y0 + (j + 0.5) * cell.dx / N;
//       for (int k = 0; k < N; k++){
//         double z = z0 + (k + 0.5) * cell.dx / N;
//         double phi = F(x, y, z);
//         double delta = peskin_delta(phi, eps);
//         if (delta == 0.0) continue;
//         Eigen::Vector3d gradPhi = gradF(x, y, z);
//         double gradNorm = gradPhi.norm();
//         double w = delta * gradNorm * dv;
//         Eigen::Vector3d pos(x, y, z);
//         moments.M0 += w;
//         moments.M1 += pos * w;
//         moments.M2 += (pos * pos.transpose()) * w;
//       }
//     }
//   }
//   return moments;
// }
// // for base cells
// void get_shape_moments(const Cell& cell, SurfaceMoments& sm, const int& N,
//                        ImplicitF F, GradientF gradF){
//   if (cell.status != CellStatus::Mixed){
//     return;
//   }
//   if (!cell.children.empty() && cell.status == CellStatus::Mixed){
//     for (const auto& child : cell.children){
//       get_shape_moments(*child, sm, N, F, gradF);
//     }
//     return;
//   }
//   if (cell.children.empty() && cell.status == CellStatus::Mixed){
//     sm += compute_surface_moments(cell, N, F, gradF);
//   }
// }
