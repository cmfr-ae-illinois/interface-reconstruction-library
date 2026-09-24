#ifndef EXAMPLES_LEVEL_SET_RECONSTRUCTION_VISUALIZATION_H_
#define EXAMPLES_LEVEL_SET_RECONSTRUCTION_VISUALIZATION_H_

#include "examples/level_set_reconstruction/level_set.h"
#include "examples/level_set_reconstruction/pu_field.h"
#include "examples/variant_advector/vtk.h"

namespace LevelSetVisualization {
// Save actual clipped LVIRA polygon centroids before optional Jibben fitting.
void storePolygonCentroids(const Data<double>& fractions,
                           const Data<IRL::SeparatorVariant>& planes,
                           Data<IRL::Pt>* centroids);

// Replace diagnostics with mean curvature and signed reference error.
void addInterfaceDiagnostics(const Data<double>& fractions,
                             const Data<IRL::SeparatorVariant>& interfaces,
                             const LevelSet& reference,
                             const std::string& label,
                             std::vector<InterfaceScalarField>* fields);

// Supported sampling cells only; point data is ready for ParaView Contour.
void writePUField(const ReconstructedPU& pu, const BasicMesh& mesh,
                  const LevelSet& reference, int sample_nx,
                  const std::string& filename);

void reconstructPUPPIC(const ReconstructedPU& pu, const Data<double>& fractions,
                       const Data<IRL::SeparatorVariant>& source,
                       const Data<IRL::Pt>& representative_points,
                       Data<IRL::SeparatorVariant>* fitted);
}  // namespace LevelSetVisualization
#endif
