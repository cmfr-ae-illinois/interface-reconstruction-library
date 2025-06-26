#include <iostream>
#include <string>
#include <functional>

#include "examples/initialize_vf/surfaces.h"
#include "examples/variant_advector/data.h"
#include "examples/variant_advector/basic_mesh.h"
#include "examples/variant_advector/reconstruction_types.h"

#ifndef EXAMPLES_INITIALIZE_VF_FUNCTIONS_H_
#define EXAMPLES_INITIALIZE_VF_FUNCTIONS_H_

using ImplicitF = std::function<double(double,double,double)>;
using GradientF = std::function<Eigen::Vector3d(double, double, double)>;
using HessianF  = std::function<Eigen::Matrix3d(double, double, double)>;

void selectSurface(const std::string& implicit_surface, ImplicitF& F,
                   GradientF& gradF, HessianF& hessF);

Eigen::Vector3d project_onto_surface(const Eigen::Vector3d& x0, ImplicitF F,
                                     GradientF gradF);

// IRL::Paraboloid generateLocalParaboloid(const IRL::Pt& a_pt,
//                                         const Eigen::Vector3d& a_gradF,
//                                         const Eigen::Matrix3d& a_hessF);

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

struct Vector3dHash {
    std::size_t operator()(const Eigen::Vector3d& v) const {
        std::hash<double> hasher;
        std::size_t h1 = hasher(v[0]);
        std::size_t h2 = hasher(v[1]);
        std::size_t h3 = hasher(v[2]);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

double evaluate_or_cache(const double& x, const double& y, const double& z, ImplicitF F,
                         std::unordered_map<Eigen::Vector3d, double, Vector3dHash>& cache);

std::vector<Eigen::Vector3d> get_sample_points(const double& x, const double& y,
                                               const double& z, const double& dx);

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

#endif