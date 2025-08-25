#include <functional>
#include <iostream>
#include <string>

#include "examples/2d_advector/basic_mesh.h"
#include "examples/2d_advector/data.h"
#include "examples/2d_advector/irl2d.h"
#include "examples/pu_analysis/surfaces.h"

#ifndef EXAMPLES_INITIALIZE_VF_FUNCTIONS_H_
#define EXAMPLES_INITIALIZE_VF_FUNCTIONS_H_

using ImplicitF = std::function<double(double, double)>;
using GradientF = std::function<Eigen::Vector2d(double, double)>;
using HessianF = std::function<Eigen::Matrix2d(double, double)>;
using Point = std::pair<double, double>;

void selectSurface(const std::string& implicit_surface, ImplicitF& F,
                   GradientF& gradF, HessianF& hessF);

BasicMesh setMesh(const int& Nx, const IRL2D::Vec& lower_domain,
                  const IRL2D::Vec& upper_domain);

Eigen::Vector2d project_onto_surface(const Eigen::Vector2d& x0, ImplicitF F,
                                     GradientF gradF);

IRL2D::Parabola build_parabola(const Eigen::Vector2d& x0, GradientF gradF,
                               HessianF hessF);

enum class CellStatus { Above, Below, Mixed };

struct Cell {
  double x, y, dx;
  int level;
  CellStatus status;
  std::vector<std::unique_ptr<Cell>> children;

  Cell(double x_, double y_, double dx_, int level_)
      : x(x_), y(y_), dx(dx_), level(level_), status(CellStatus::Mixed) {}

  bool is_leaf() const { return children.empty(); }

  double area() const { return dx * dx; }
};

struct PointHash {
  std::size_t operator()(const Point& p) const {
    std::hash<double> hasher;
    std::size_t h1 = hasher(p.first);
    std::size_t h2 = hasher(p.second);
    return h1 ^ (h2 << 1);
  }
};

double evaluate_or_cache(double x, double y, ImplicitF F,
                         std::unordered_map<Point, double, PointHash>& cache);

CellStatus get_cell_status(
    const Cell& cell, ImplicitF F, GradientF gradF, HessianF hessF,
    std::unordered_map<Point, double, PointHash>& F_cache);

void refine_cell(std::unique_ptr<Cell>& cell, ImplicitF F, const int& max_level,
                 GradientF gradF, HessianF hessF,
                 std::unordered_map<Point, double, PointHash>& F_cache);

std::vector<std::unique_ptr<Cell>> refine_grid(
    const BasicMesh& mesh, ImplicitF F, const int& max_level, GradientF gradF,
    HessianF hessF, std::unordered_map<Point, double, PointHash>& F_cache);

std::pair<double, double> get_mixed_area(const double& x, const double& y,
                                         const double& dx, ImplicitF F,
                                         GradientF gradF, HessianF hessF);

struct AreaResult {
  double inside = 0;
  double outside = 0;

  double total() const { return inside + outside; }
  double volume_fraction() const {
    double denom = total();
    return (denom > 0.0) ? inside / denom : 0.0;
  }
};

void get_cell_area(const Cell& cell, AreaResult& result, ImplicitF F,
                   GradientF gradF, HessianF hessF);

void collect_leaf_cells(const Cell& cell, std::vector<const Cell*>& output,
                        std::vector<const Cell*>& output_mixed);

void collect_mixed_leaf_cells(const Cell& cell,
                              std::vector<const Cell*>& output);

void write_vtr(const std::string& filepath, const Data<double>& vf,
               const BasicMesh& mesh);

void write_vtu(const std::string& filepath,
               const std::vector<const Cell*>& cells);

#endif