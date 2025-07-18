#include <iostream>
#include <string>
#include <functional>

#include "examples/initialize_vf/surfaces.h"
#include "examples/variant_advector/data.h"
#include "examples/variant_advector/basic_mesh.h"
#include "examples/variant_advector/reconstruction_types.h"

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/surface_mesher/marching_cubes.h"

#ifndef EXAMPLES_INITIALIZE_VF_FUNCTIONS_H_
#define EXAMPLES_INITIALIZE_VF_FUNCTIONS_H_

using ImplicitF = std::function<double(double,double,double)>;
using GradientF = std::function<Eigen::Vector3d(double, double, double)>;
using HessianF  = std::function<Eigen::Matrix3d(double, double, double)>;


std::unique_ptr<Surface> makeSurface(const std::string& name);

void initialize(std::unique_ptr<Surface>& surface,
                ImplicitF& F, GradientF& gradF,
                HessianF& hessF, BasicMesh& mesh);

Eigen::Vector3d project_onto_surface(const Eigen::Vector3d& x0, ImplicitF F,
                                     GradientF gradF);

enum class CellStatus { Above, Below, Mixed};

struct Cell {
  double x, y, z, dx;
  int level;
  CellStatus status;
  std::vector<std::unique_ptr<Cell>> children;

  Cell(double x_, double y_, double z_, double dx_, int level_)
    : x(x_), y(y_), z(z_), dx(dx_), level(level_), status(CellStatus::Mixed) {}
  
  bool is_leaf() const {
    return children.empty();
  }

  double area() const {
    return dx * dx * dx;
  }
};

// struct Vector3dHash {
//     std::size_t operator()(const Eigen::Vector3d& v) const {
//         std::hash<double> hasher;
//         std::size_t h1 = hasher(v[0]);
//         std::size_t h2 = hasher(v[1]);
//         std::size_t h3 = hasher(v[2]);
//         return h1 ^ (h2 << 1) ^ (h3 << 2);
//     }
// };

struct Vector3dHash {
    std::size_t operator()(const Eigen::Vector3d& v) const {
        std::hash<double> hasher;
        std::size_t seed = 0;
        hash_combine(seed, hasher(v[0]));
        hash_combine(seed, hasher(v[1]));
        hash_combine(seed, hasher(v[2]));
        return seed;
    }

private:
    static void hash_combine(std::size_t& seed, std::size_t hash) {
        seed ^= hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
};

double evaluate_or_cache(const Eigen::Vector3d& pt, ImplicitF F,
                         std::unordered_map<Eigen::Vector3d, double, Vector3dHash>& cache);

// std::vector<Eigen::Vector3d> get_sample_points(const double& x, const double& y,
//                                                const double& z, const double& dx);

std::vector<Eigen::Vector3d> get_sample_points(const double& x, const double& y,
                                               const double& z, const double& dx,
                                               const bool& use_stencil);

CellStatus get_cell_status(const Cell& cell, ImplicitF F,
                           GradientF gradF, HessianF hessF,
                           std::unordered_map<Eigen::Vector3d, double, Vector3dHash>& F_cache);
                                            
void refine_cell(std::unique_ptr<Cell>& cell, ImplicitF F,
                 GradientF gradF, HessianF hessF, const int& max_level,
                 std::unordered_map<Eigen::Vector3d, double, Vector3dHash>& F_cache);

std::vector<std::unique_ptr<Cell>> refine_grid(const BasicMesh& mesh, ImplicitF F,
                                               GradientF gradF, HessianF hessF, const int& max_level,
                                               std::unordered_map<Eigen::Vector3d, double, Vector3dHash>& F_cache);

std::pair<double, double> get_mixed_area(const Eigen::Vector3d x0,
                                         const double& dx, ImplicitF F,
                                         GradientF gradF, HessianF hessF);

struct AreaResult{
  double inside = 0;
  double outside = 0;

  double total() const {return inside + outside;}
  double volume_fraction() const {
    double denom = total();
    return (denom > 0.0) ? inside / denom : 0.0;
  }
};

void get_cell_area(const Cell& cell, AreaResult& result, ImplicitF F,
                   GradientF gradF, HessianF hessF);

void collect_leaf_cells(const Cell& cell, std::vector<const Cell*>& output,
                        std::vector<const Cell*>& output_mixed);

void write_vtr(const std::string& filepath, const Data<double>& vf, const BasicMesh& mesh);

void write_vtu(const std::string& filepath, const std::vector<const Cell*>& cells);

class vfInitializer{
  public:
    vfInitializer(const int& Nx, const int& max_refine_level,
                  const std::string& surface_name,
                  const std::string& vtk_output_path);
    
    void run();

    const Data<double>& getVolumeFractions() const { return vf_; }
    const std::vector<const Cell*>& getLeafCells() const { return leaf_cells_; }
    const std::vector<const Cell*>& getMixedLeafCells() const { return mixed_leaf_cells_; }
    double getTotalVolume() const { return volume_amr_; }
    const Surface* getSurface() const { return surface_.get();}
    std::chrono::duration<double> getElapsedTime() const { return elapsed_; }
    const std::vector<std::unique_ptr<Cell>>& getRefinedGrid() { return grid_; }

  private:
    int Nx_, max_refine_level_;
    std::string surface_name_, vtk_output_path_;
    BasicMesh mesh_;
    std::vector<std::unique_ptr<Cell>> grid_;
    Data<double> vf_;
    double volume_amr_ = 0.0;
    std::chrono::duration<double> elapsed_;
    std::vector<const Cell*> leaf_cells_, mixed_leaf_cells_;
    std::unique_ptr<Surface> surface_;
};

void performConvergence(const std::string& csv_path, const int& max_refine_level,
                        const vfInitializer& vfi);

class Triangle{
  public:
    Triangle(const Eigen::Vector3d& v0,
             const Eigen::Vector3d& v1,
             const Eigen::Vector3d& v2);
            
    double computeM0() const;
    Eigen::Vector3d computeM1() const;
    Eigen::Matrix3d computeM2() const;
  
  private:
    Eigen::Vector3d v0_, v1_, v2_;
};

struct Metrics{
  // geometric surface moments
  double sM0 = 0.0;
  Eigen::Vector3d sM1 = Eigen::Vector3d::Zero();
  Eigen::Matrix3d sM2 = Eigen::Matrix3d::Zero();

  // geometric volume moments
  double vM0 = 0.0;
  Eigen::Vector3d vM1 = Eigen::Vector3d::Zero();
  Eigen::Matrix3d vM2 = Eigen::Matrix3d::Zero();

  // other metrics

  // operator definitions
  Metrics& operator+=(const Metrics& other){
    sM0 += other.sM0;
    sM1 += other.sM1;
    sM2 += other.sM2;
    vM0 += other.vM0;
    vM1 += other.vM1;
    vM2 += other.vM2;
    return *this;
  }
  
};

class SurfaceMetrics{
  public:
    using ImplicitF_pt = std::function<double(IRL::Pt)>;

    SurfaceMetrics(const Cell& base_cell, ImplicitF F, const int num_divisions);

    Metrics computeSurfaceMetrics();

  private:
    const Cell& base_cell_;
    ImplicitF_pt F_; 
    int num_divisions_; 
    Metrics metrics_;
    std::vector<const Cell*> collectMixedLeafCells() const;
};

// surface moment calculations
// double peskin_delta(const double& phi, const double& epsilon);
// double gaussian_delta_truncated(double phi, double sigma);
// struct SurfaceMoments {
//     double M0 = 0;            
//     Eigen::Vector3d M1 = Eigen::Vector3d::Zero(); 
//     Eigen::Matrix3d M2 = Eigen::Matrix3d::Zero(); 
//     SurfaceMoments& operator+=(const SurfaceMoments& other) {
//       M0 += other.M0;
//       M1 += other.M1;
//       M2 += other.M2;
//       return *this;
//     }
//     SurfaceMoments operator+(const SurfaceMoments& other) const {
//       SurfaceMoments result = *this;
//       result += other;
//       return result;
//     }
// };
// SurfaceMoments compute_surface_moments(const Cell& cell, const int& N, 
//                                        ImplicitF F, GradientF gradF);                                
// void get_shape_moments(const Cell& cell, SurfaceMoments& sm, const int& N,
//                        ImplicitF F, GradientF gradF);

#endif