#include <vtkCellCenters.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkDataArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/ml_classification/ml_classifier_notorch.h"
#include "irl/ml_classification/ml_classifier_weights_and_biases.h"

#include <iostream>
#include <vector>

static const int N = 5;
static const int cid = N / 2;

using scalar_array = std::array<std::array<std::array<double, N>, N>, N>;
using vector_array =
    std::array<std::array<std::array<std::array<double, 3>, N>, N>, N>;

inline void reflect(scalar_array& vfrac, vector_array& liq_bary,
                    const int dir) {
  const scalar_array tmp_vfrac = vfrac;
  const vector_array tmp_liq_bary = liq_bary;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      for (int k = 0; k < N; ++k) {
        const int mirror_i = (dir == 0) ? (N - 1 - i) : i;
        const int mirror_j = (dir == 1) ? (N - 1 - j) : j;
        const int mirror_k = (dir == 2) ? (N - 1 - k) : k;
        vfrac[i][j][k] = tmp_vfrac[mirror_i][mirror_j][mirror_k];
        liq_bary[i][j][k][0] =
            (dir == 0) ? -tmp_liq_bary[mirror_i][mirror_j][mirror_k][0]
                       : tmp_liq_bary[mirror_i][mirror_j][mirror_k][0];
        liq_bary[i][j][k][1] =
            (dir == 1) ? -tmp_liq_bary[mirror_i][mirror_j][mirror_k][1]
                       : tmp_liq_bary[mirror_i][mirror_j][mirror_k][1];
        liq_bary[i][j][k][2] =
            (dir == 2) ? -tmp_liq_bary[mirror_i][mirror_j][mirror_k][2]
                       : tmp_liq_bary[mirror_i][mirror_j][mirror_k][2];
      }
    }
  }
}

inline void permute(scalar_array& vfrac, vector_array& liq_bary,
                    const int dir) {
  if (dir == 1) {
    const scalar_array tmp_vfrac = vfrac;
    const vector_array tmp_liq_bary = liq_bary;
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        for (int k = 0; k < N; ++k) {
          vfrac[i][j][k] = tmp_vfrac[k][i][j];
          liq_bary[i][j][k][0] = tmp_liq_bary[k][i][j][1];
          liq_bary[i][j][k][1] = tmp_liq_bary[k][i][j][2];
          liq_bary[i][j][k][2] = tmp_liq_bary[k][i][j][0];
        }
      }
    }
  } else if (dir == 2) {
    const scalar_array tmp_vfrac = vfrac;
    const vector_array tmp_liq_bary = liq_bary;
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        for (int k = 0; k < N; ++k) {
          vfrac[i][j][k] = tmp_vfrac[j][k][i];
          liq_bary[i][j][k][0] = tmp_liq_bary[j][k][i][2];
          liq_bary[i][j][k][1] = tmp_liq_bary[j][k][i][0];
          liq_bary[i][j][k][2] = tmp_liq_bary[j][k][i][1];
        }
      }
    }
  }
}

inline void swap_xy(scalar_array& vfrac, vector_array& liq_bary) {
  const scalar_array tmp_vfrac = vfrac;
  const vector_array tmp_liq_bary = liq_bary;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      for (int k = 0; k < N; ++k) {
        vfrac[i][j][k] = tmp_vfrac[j][i][k];
        liq_bary[i][j][k][0] = tmp_liq_bary[j][i][k][1];
        liq_bary[i][j][k][1] = tmp_liq_bary[j][i][k][0];
        liq_bary[i][j][k][2] = tmp_liq_bary[j][i][k][2];
      }
    }
  }
}

inline int largest_value_index(const double v0, const double v1,
                               const double v2) {
  double largest = v0;
  int idx = 0;
  if (v1 > largest) {
    largest = v1;
    idx = 1;
  }
  if (v2 > largest) {
    largest = v2;
    idx = 2;
  }
  return idx;
}

inline int torchfree_classifier(
    scalar_array& vfrac,
    vector_array& liq_bary) {  // Clean up nonconnected liquid regions

  bool visited[N][N][N] = {false};
  std::vector<std::array<int, 3>> q;
  q.reserve(N * N * N);

  assert(N % 2 == 1);  // ensure we have a central cell
  const double epsilon_connect = 1.0e-12;

  auto in_bounds = [&](int i, int j, int k) {
    return (i >= 0 && i < N && j >= 0 && j < N && k >= 0 && k < N);
  };
  static const int n6[6][3] = {{+1, 0, 0}, {-1, 0, 0}, {0, +1, 0},
                               {0, -1, 0}, {0, 0, +1}, {0, 0, -1}};

  // Seed with central cell
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      for (int k = 0; k < N; ++k) {
        visited[i][j][k] = false;
      }
    }
  }
  visited[cid][cid][cid] = true;
  q.push_back({cid, cid, cid});

  // Flood fill
  for (int qi = 0; qi < q.size(); ++qi) {
    const int i = q[qi][0];
    const int j = q[qi][1];
    const int k = q[qi][2];

    for (int n = 0; n < 6; ++n) {
      const int ni = i + n6[n][0];
      const int nj = j + n6[n][1];
      const int nk = k + n6[n][2];
      if (!in_bounds(ni, nj, nk)) continue;
      if (visited[ni][nj][nk]) continue;
      if (vfrac[ni][nj][nk] > epsilon_connect) {
        visited[ni][nj][nk] = true;
        q.push_back({ni, nj, nk});
      }
    }
  }

  // Zero everything not connected
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      for (int k = 0; k < N; ++k) {
        if (visited[i][j][k]) continue;
        if (vfrac[i][j][k] > 0.0) {
          vfrac[i][j][k] = 0.0;
          liq_bary[i][j][k][0] = 0.0;
          liq_bary[i][j][k][1] = 0.0;
          liq_bary[i][j][k][2] = 0.0;
        }
      }
    }
  }

  // Compute global liquid centroid in stencil
  auto global_liq_bary = IRL::Pt(0.0, 0.0, 0.0);
  double total_volume = 0.0;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      for (int k = 0; k < N; ++k) {
        total_volume += vfrac[i][j][k];
        global_liq_bary += IRL::Pt(liq_bary[i][j][k][0], liq_bary[i][j][k][1],
                                   liq_bary[i][j][k][2]);
      }
    }
  }
  global_liq_bary *= 1.0 / IRL::safelyTiny(total_volume);

  // Reflect into first octant: ensure cx,cy,cz >= 0 using x/y/z reflections.
  if (global_liq_bary[0] < 0.0) {
    reflect(vfrac, liq_bary, 0);
    global_liq_bary[0] = -global_liq_bary[0];
  }
  if (global_liq_bary[1] < 0.0) {
    reflect(vfrac, liq_bary, 1);
    global_liq_bary[1] = -global_liq_bary[1];
  }
  if (global_liq_bary[2] < 0.0) {
    reflect(vfrac, liq_bary, 2);
    global_liq_bary[2] = -global_liq_bary[2];
  }

  const int idx = largest_value_index(global_liq_bary[2], global_liq_bary[0],
                                      global_liq_bary[1]);
  if (idx == 1) {
    permute(vfrac, liq_bary, 1);
    // (cx,cy,cz) -> (cy,cz,cx)
    const double tmp = global_liq_bary[0];
    global_liq_bary[0] = global_liq_bary[1];
    global_liq_bary[1] = global_liq_bary[2];
    global_liq_bary[2] = tmp;
  } else if (idx == 2) {
    permute(vfrac, liq_bary, 2);
    // (cx,cy,cz) -> (cz,cx,cy)
    const double tmp = global_liq_bary[0];
    global_liq_bary[0] = global_liq_bary[2];
    global_liq_bary[2] = global_liq_bary[1];
    global_liq_bary[1] = tmp;
  }

  if (global_liq_bary[1] > global_liq_bary[0] /*+ 1e-8*/) {
    // count small differences
    swap_xy(vfrac, liq_bary);
    const double tmp = global_liq_bary[0];
    global_liq_bary[0] = global_liq_bary[1];
    global_liq_bary[1] = tmp;
  }

  // Flatten stencil into 1D vector
  // size = vfrac + (mx, my, mz) for each cell + 3 eigenvalues
  const int n_inputs = 4 * 5 * 5 * 5;
  std::vector<float> flattened_state;
  flattened_state.reserve(n_inputs);
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      for (int k = 0; k < 5; ++k) {
        flattened_state.push_back(static_cast<float>(vfrac[i][j][k]));
        flattened_state.push_back(static_cast<float>(liq_bary[i][j][k][0]));
        flattened_state.push_back(static_cast<float>(liq_bary[i][j][k][1]));
        flattened_state.push_back(static_cast<float>(liq_bary[i][j][k][2]));
      }
    }
  }

  // // Add eigenvalues of reconstructed inertia tensor
  // Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
  // for (int i = 0; i < 5; ++i) {
  //   for (int j = 0; j < 5; ++j) {
  //     for (int k = 0; k < 5; ++k) {
  //       const auto local_liq_bary =
  //           IRL::Pt(liq_bary[i][j][k][0], liq_bary[i][j][k][1],
  //                   liq_bary[i][j][k][2]) /
  //           IRL::safelyEpsilon(vfrac[i][j][k]);
  //       const auto local_rel_liq_bary = local_liq_bary - global_liq_bary;
  //       for (int a = 0; a < 3; ++a) {
  //         for (int b = 0; b < 3; ++b) {
  //           I(a, b) +=
  //               vfrac[i][j][k] * local_rel_liq_bary[a] *
  //               local_rel_liq_bary[b];
  //         }
  //       }
  //     }
  //   }
  // }
  // I = I.trace() * Eigen::Matrix3d::Identity() - I;
  // Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(I);
  // Eigen::Vector3d evals = solver.eigenvalues();
  // std::sort(evals.data(), evals.data() + 3, std::greater<double>());
  // double I1 = evals[0], I2 = evals[1], I3 = evals[2];
  // flattened_state.push_back(static_cast<float>(I1));
  // flattened_state.push_back(static_cast<float>(I2));
  // flattened_state.push_back(static_cast<float>(I3));

  // Classify
  // std::vector<float> out_probs;
  // IRL::MLClassifierNoTorch classifier(N);
  // int predicted_class = classifier.classify(flattened_state, &out_probs);

  // Classify
  std::vector<double> logits;
  // forwardLogits(flattened_state, logits);
  std::array<double, IRL::mlclassifier::hidden_size1> h1{};
  std::array<double, IRL::mlclassifier::hidden_size2> h2{};
  std::array<double, IRL::mlclassifier::hidden_size3> h3{};

  for (int j = 0; j < IRL::mlclassifier::hidden_size1; ++j) {
    h1[j] = IRL::mlclassifier::fc1_bias[j];

    for (int i = 0; i < IRL::mlclassifier::input_size; ++i) {
      h1[j] += static_cast<double>(flattened_state[i]) *
               IRL::mlclassifier::fc1_weight[j][i];
    }

    if (h1[j] < 0.0) {
      h1[j] = 0.0;
    }
  }

  for (int j = 0; j < IRL::mlclassifier::hidden_size2; ++j) {
    h2[j] = IRL::mlclassifier::fc2_bias[j];

    for (int i = 0; i < IRL::mlclassifier::hidden_size1; ++i) {
      h2[j] += h1[i] * IRL::mlclassifier::fc2_weight[j][i];
    }

    if (h2[j] < 0.0) {
      h2[j] = 0.0;
    }
  }

  for (int j = 0; j < IRL::mlclassifier::hidden_size3; ++j) {
    h3[j] = IRL::mlclassifier::fc3_bias[j];

    for (int i = 0; i < IRL::mlclassifier::hidden_size2; ++i) {
      h3[j] += h2[i] * IRL::mlclassifier::fc3_weight[j][i];
    }

    if (h3[j] < 0.0) {
      h3[j] = 0.0;
    }
  }

  logits.assign(IRL::mlclassifier::output_size, 0.0);

  for (int j = 0; j < IRL::mlclassifier::output_size; ++j) {
    logits[j] = IRL::mlclassifier::fc4_bias[j];

    for (int i = 0; i < IRL::mlclassifier::hidden_size3; ++i) {
      logits[j] += h3[i] * IRL::mlclassifier::fc4_weight[j][i];
    }
  }

  int predicted_class = std::distance(
      logits.begin(), std::max_element(logits.begin(), logits.end()));

  return predicted_class;
}

inline std::string classid_to_string(const int id) {
  if (id == 1) {
    return "ligament";
  } else if (id == 2) {
    return "droplet";
  } else if (id == 3) {
    return "sheet/film";
  } else if (id == 4) {
    return "ligament end";
  } else if (id == 5) {
    return "sheet end";
  } else if (id == 0) {
    return "well-resolved";
  }
  return "unknown";
}

int main(int argc, char* argv[]) {
  scalar_array vfrac;
  vector_array liq_bary;

  std::string datafile =
      "../examples/ml_minimal_classifier_notorch/test_data/data.vtr";
  std::string plicfile =
      "../examples/ml_minimal_classifier_notorch/test_data/plic.vtu";

  auto data_reader = vtkSmartPointer<vtkXMLRectilinearGridReader>::New();
  data_reader->SetFileName(datafile.c_str());
  data_reader->Update();

  vtkRectilinearGrid* grid = data_reader->GetOutput();
  vtkCellData* cellData = grid->GetCellData();
  vtkDataArray* xCoords = grid->GetXCoordinates();
  vtkDataArray* yCoords = grid->GetYCoordinates();
  vtkDataArray* zCoords = grid->GetZCoordinates();
  vtkDataArray* vfrac_array =
      vtkDataArray::SafeDownCast(cellData->GetArray("VOF"));
  vtkDataArray* liq_bary_array =
      vtkDataArray::SafeDownCast(cellData->GetArray("liquid_bary"));
  int dims[3];
  grid->GetDimensions(dims);  // dims are for points, so cells are dims-1
  int cellDims[3] = {dims[0] - 1, dims[1] - 1, dims[2] - 1};
  auto interface_type = vtkSmartPointer<vtkIntArray>::New();
  interface_type->SetName("class");
  interface_type->SetNumberOfComponents(1);
  interface_type->SetNumberOfTuples(grid->GetNumberOfCells());
  for (int i = 0; i < grid->GetNumberOfCells(); i++) {
    interface_type->SetValue(i, -1);
  }

  std::cout << "Files read; processing data...\n";
  for (int i = 2; i < cellDims[0] - 2; ++i) {
    for (int j = 2; j < cellDims[1] - 2; ++j) {
      for (int k = 2; k < cellDims[2] - 2; ++k) {
        const int idx = i + j * cellDims[0] + k * cellDims[0] * cellDims[1];
        const double cell_vfrac = vfrac_array->GetComponent(idx, 0);
        if (cell_vfrac < IRL::global_constants::VF_LOW ||
            cell_vfrac > IRL::global_constants::VF_HIGH) {
          continue;
        }
        std::cout << "Processing cell "
                  << k + j * cellDims[2] + i * cellDims[2] * cellDims[1]
                  << " / " << cellDims[0] * cellDims[1] * cellDims[2] << "\r"
                  << std::flush;
        const double x0 = xCoords->GetComponent(i, 0);
        const double y0 = yCoords->GetComponent(j, 0);
        const double z0 = zCoords->GetComponent(k, 0);
        const double x1 = xCoords->GetComponent(i + 1, 0);
        const double y1 = yCoords->GetComponent(j + 1, 0);
        const double z1 = zCoords->GetComponent(k + 1, 0);
        const double dx = x1 - x0;
        const double dy = y1 - y0;
        const double dz = z1 - z0;
        const double cell_volume = std::abs(dx * dy * dz);
        IRL::Pt cell_center(0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * (z0 + z1));
        for (int ii = 0; ii < 5; ++ii) {
          for (int jj = 0; jj < 5; ++jj) {
            for (int kk = 0; kk < 5; ++kk) {
              const int local_idx = (i + ii - 2) + (j + jj - 2) * cellDims[0] +
                                    (k + kk - 2) * cellDims[0] * cellDims[1];
              vfrac[ii][jj][kk] = vfrac_array->GetComponent(local_idx, 0);
              IRL::Pt cell_liq_bary(liq_bary_array->GetComponent(local_idx, 0),
                                    liq_bary_array->GetComponent(local_idx, 1),
                                    liq_bary_array->GetComponent(local_idx, 2));
              cell_liq_bary -= cell_center;
              cell_liq_bary[0] /= dx;
              cell_liq_bary[1] /= dy;
              cell_liq_bary[2] /= dz;
              cell_liq_bary *= vfrac[ii][jj][kk];
              liq_bary[ii][jj][kk][0] = cell_liq_bary[0];
              liq_bary[ii][jj][kk][1] = cell_liq_bary[1];
              liq_bary[ii][jj][kk][2] = cell_liq_bary[2];
            }
          }
        }
        interface_type->SetValue(idx, torchfree_classifier(vfrac, liq_bary));
      }
    }
  }

  int count_classes[6] = {0};
  for (int i = 0; i < grid->GetNumberOfCells(); i++) {
    int class_id = interface_type->GetValue(i);
    if (class_id >= 0 && class_id < 6) {
      count_classes[class_id]++;
    }
  }

  std::cout << "\nClass counts:\n";
  for (int i = 0; i < 6; ++i) {
    std::cout << "  " << classid_to_string(i) << "\t: " << count_classes[i]
              << "\n";
  }

  auto plic_reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  plic_reader->SetFileName(plicfile.c_str());
  plic_reader->Update();
  vtkUnstructuredGrid* plic_grid = plic_reader->GetOutput();

  auto plic_interface_type = vtkSmartPointer<vtkIntArray>::New();
  plic_interface_type->SetName("class");
  plic_interface_type->SetNumberOfComponents(1);
  plic_interface_type->SetNumberOfTuples(plic_grid->GetNumberOfCells());
  for (int i = 0; i < plic_grid->GetNumberOfCells(); i++) {
    plic_interface_type->SetValue(i, -1);
  }

  // Cell locator
  auto locator = vtkSmartPointer<vtkCellLocator>::New();
  locator->SetDataSet(grid);
  locator->BuildLocator();
  auto plic_cell_center_filter = vtkSmartPointer<vtkCellCenters>::New();
  plic_cell_center_filter->SetInputData(plic_grid);
  plic_cell_center_filter->Update();
  auto plic_cell_centers = plic_cell_center_filter->GetOutput()->GetPoints();

  // Convert interface_type to PLIC array
  for (int i = 0; i < plic_grid->GetNumberOfCells(); i++) {
    auto type = interface_type->GetValue(
        locator->FindCell(plic_cell_centers->GetPoint(i)));
    plic_interface_type->SetValue(i, type);
  }

  // Write PLIC file
  plic_grid->GetCellData()->AddArray(plic_interface_type);
  auto plic_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  plic_writer->SetFileName("plic_with_classes.vtu");
  plic_writer->SetInputData(plic_grid);
  plic_writer->Write();

  return 0;
}