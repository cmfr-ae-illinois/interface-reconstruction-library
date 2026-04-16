#include <chrono>
#include <ctime>
#include <iostream>
#include <string>
#include <random>
#include <Eigen/Dense>
#include <vector>
#include <numeric>
#include <iomanip>
#include <cmath>
#include <vtkCellCenters.h>

#include "irl/ml_classification/data_gen.h"
#include "irl/ml_classification/inertia_classifier.h"
#include "irl/ml_classification/vtk_in.h"
#include "irl/ml_classification/ml_classifier.h"
#include "irl/ml_classification/e3nn.h"

// Utility: apply 90° rotation about Z axis to barycenter vectors and positions
torch::Tensor rotate90_z(const torch::Tensor& v) {
    // v: [B, N, 1, 3] or [N, 1, 3]
    auto R = torch::tensor({
        {0.0f, -1.0f, 0.0f},
        {1.0f,  0.0f, 0.0f},
        {0.0f,  0.0f, 1.0f}
    }, v.options());
    if (v.dim() == 4)
        return torch::einsum("bnvc,cd->bnvd", {v, R});
    else
        return torch::einsum("nvc,cd->nvd", {v, R});
}

// Helper: print full stencil data
void print_stencil(const torch::Tensor& vof, const torch::Tensor& bary, int stencil_size) {
    int64_t N = int64_t(stencil_size) * stencil_size * stencil_size;

    // Create local lvalue tensors before accessing
    torch::Tensor vof_flat  = vof.squeeze(0).squeeze(-1).contiguous();   // [N]
    torch::Tensor bary_flat = bary.squeeze(0).squeeze(1).contiguous();   // [N,3]

    auto vof_acc  = vof_flat.accessor<float, 1>();
    auto bary_acc = bary_flat.accessor<float, 2>();

    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < stencil_size; ++i) {
        for (int j = 0; j < stencil_size; ++j) {
            for (int k = 0; k < stencil_size; ++k) {
                int n = (i * stencil_size + j) * stencil_size + k;
                std::cout << "[" << i << ", " << j << ", " << k << "] "
                          << vof_acc[n] << ", "
                          << bary_acc[n][0] << ", "
                          << bary_acc[n][1] << ", "
                          << bary_acc[n][2] << "\n";
            }
        }
    }
}

/*

void evaluateMLClassifier(int stencil_size, int n_samples_per_class, IRL::MLClassifier& ml, IRL::Data_gen& data_gen) {
    const std::vector<std::string> classNames = {"Paraboloid", "Ligaments", "Spheres", "Sheets"};
    const int n_classes = classNames.size();

    for (int true_class = 0; true_class < n_classes; ++true_class) {
        std::vector<float> sum_probs(n_classes, 0.0f);

        for (int i = 0; i < n_samples_per_class; ++i) {
            // generate a synthetic state
            std::vector<double> flattened_state = data_gen.generateState(true_class, stencil_size, false);

            // classify and get full probability vector
            std::vector<float> probs;
            int detectedClass = ml.classify(flattened_state, &probs);

            // accumulate class probabilities
            for (int c = 0; c < n_classes; ++c)
                sum_probs[c] += probs[c];
        }

        // normalize to average over samples and convert to %
        for (float& p : sum_probs) p = (p / n_samples_per_class) * 100.0f;

        // print results
        std::cout << "\nTrue class: " << classNames[true_class] << std::endl;
        for (int c = 0; c < n_classes; ++c) {
            std::cout << "  " << classNames[c] << ": " << sum_probs[c] << " %" << std::endl;
        }
    }
}

*/


static std::string classNameFromType(int t) {
    switch (t) {
        case 1: return "Cylinder";
        case 2: return "Sphere";
        case 3: return "Sheet";
        default: return "Paraboloid";
    }
}

// Per-cell data: volume fraction and first moments
struct CellData {
    double vfrac;
    double mx, my, mz;  // first moments in x, y, z
};

// Read one stencil from a flat [vfrac, mx, my, mz] array into CellData array
static void unpackStencil(const std::vector<double>& flat,
                          std::vector<CellData>& stencil, int N)
{
    const int nCells = N * N * N;
    stencil.resize(nCells);

    // Expect include_Moments >= 1, so layout is [vfrac, mx, my, mz] per cell
    for (int idx = 0; idx < nCells; ++idx) {
        stencil[idx].vfrac = flat[4 * idx + 0];
        stencil[idx].mx    = flat[4 * idx + 1];
        stencil[idx].my    = flat[4 * idx + 2];
        stencil[idx].mz    = flat[4 * idx + 3];
    }
}

void simulation_stability_analysis() {
    int stencil_size = 5;

    // Net Parameters
    int input_size = 4 * stencil_size * stencil_size * stencil_size; // For vof only stencil_size * stencil_size * stencil_size; // 27 if stencil_size=3 and only vof
    int hidden_size1 = 256;
    int hidden_size2 = 64;
    int hidden_size3 = 32;
    int output_size = 4;

    //Training parameters
    double learning_rate = 0.001; //was 0.01 for SGD optimizer
    int batch_size = 64;
    int epochs = 5;

    //Data parameters
    int no_batches = 4096/2;
    int include_Moments = 1;
    double paraboloid_coeff_stddev = 0.1;
    double sheet_coeff_stddev = 0.1;
    double max_sheet_thickness = 0.5;
    double sheet_thickness_stddev = 0.0;
    double max_cylinder_radius = 0.5;
    double cylinder_radius_stddev = 0.0;
    bool include_truncated_cylinder = true;
    double max_sphere_radius = 0.5;
    double sphere_radius_stddev = 0.0;

    //IRL::InertiaClassifier inertia_classifier(stencil_size, 1, 0.85, 1.5);
    //IRL::MLClassifier_E3NN ml(stencil_size, hidden_size1, hidden_size2, hidden_size3, output_size);
    IRL::MLClassifier ml(stencil_size, input_size, hidden_size1, hidden_size2, hidden_size3, output_size);
    
    ml.updateDataParameters(no_batches, include_Moments,
                            paraboloid_coeff_stddev,
                            sheet_coeff_stddev,
                            max_sheet_thickness, sheet_thickness_stddev,
                            max_cylinder_radius, cylinder_radius_stddev, include_truncated_cylinder,
                            max_sphere_radius, sphere_radius_stddev, false);                    
    
    ml.loadDataset("/home/quirin/mlcfd/Datasets/float/From1/s5_3M/data/data.bin");
    //ml.generateDataset();
    int canonicalize_symmetries = 48;
    bool preProcess = true;
    ml.canonicalize_data(canonicalize_symmetries);
    ml.updateTrainingParameters(learning_rate, batch_size, epochs);
    std::string filename = "/home/quirin/mlcfd/Repositories/jet/nga.case";


    const int no_runs = 10;
    std::vector<std::vector<int>> preds;  // preds[m][t] = class

    for (int r = 0; r < no_runs; ++r) {
        ml.trainModel();

        std::vector<int> savedClasses;
        //IRL::classify_simulation(ml, filename, canonicalize_symmetries, preProcess, &savedClasses);

        if (r > 0 && savedClasses.size() != preds[0].size()) {
            throw std::runtime_error("Saved class vector size differs between runs. "
                                     "Use cell-id alignment (recommended) or ensure deterministic selection.");
        }

        preds.push_back(std::move(savedClasses));
    }

    // Compute mean instability u = 1 - (max vote fraction)
    const size_t N = preds[0].size();
    double sum_instability = 0.0;

    for (size_t t = 0; t < N; ++t) {
        std::array<int,4> counts{0,0,0,0};
        for (int run = 0; run < no_runs; ++run) {
            int c = preds[run][t];
            if (c < 0 || c >= 4) continue; // should not happen
            counts[c]++;
        }
        int max_count = std::max(std::max(counts[0], counts[1]), std::max(counts[2], counts[3]));
        double instability = 1.0 - static_cast<double>(max_count) / static_cast<double>(no_runs);
        sum_instability += instability;
    }

    double mean_instability = (N > 0) ? (sum_instability / static_cast<double>(N)) : 0.0;
    std::cout << "Mean per-cell instability over " << N << " cells = " << mean_instability << "\n";
}

static void print_flat_state(const std::vector<float>& state, const std::string& title) {
    std::cout << "\n==== " << title << " (flat) ====\n";
    for (size_t i = 0; i < state.size(); ++i) {
        std::cout << std::setw(4) << i << ": "
                  << std::setw(14) << std::setprecision(8) << state[i] << '\n';
    }
}

static void print_decoded_state(const std::vector<float>& flat,
                                int N,
                                int include_moments,
                                bool include_Eigenvalues,
                                const std::string& title,
                                bool stored_as_local_centroid)
{
    std::vector<IRL::CellData> stencil;
    IRL::SecondMoments I{};
    IRL::SecondMoments* Ip = (include_moments >= 2) ? &I : nullptr;
    IRL::Eigenvalues eig{};
    IRL::Eigenvalues* eigp = include_Eigenvalues ? &eig : nullptr;

    IRL::unpackStencil(flat, stencil, Ip, eigp, N, include_moments, include_Eigenvalues);

    const float c0 = 0.5f * (static_cast<float>(N) - 1.0f);

    std::cout << "\n==== " << title << " ====\n";

    for (int i = 0; i < N; ++i) {
        std::cout << "\n----- i = " << i << " -----\n";
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                const auto& c = stencil[IRL::cellIndex(i, j, k, N)];

                const float x_center = static_cast<float>(i) - c0;
                const float y_center = static_cast<float>(j) - c0;
                const float z_center = static_cast<float>(k) - c0;

                std::cout << "[" << i << "," << j << "," << k << "] "
                          << "v=" << c.vfrac;

                if (include_moments >= 1) {
                    if (!stored_as_local_centroid) {
                        std::cout << "  m1=(" << c.mx << ", " << c.my << ", " << c.mz << ")";
                        if (std::fabs(c.vfrac) > 1e-12f) {
                            const float cx_global = c.mx / c.vfrac;
                            const float cy_global = c.my / c.vfrac;
                            const float cz_global = c.mz / c.vfrac;

                            const float cx_local = cx_global - x_center;
                            const float cy_local = cy_global - y_center;
                            const float cz_local = cz_global - z_center;

                            std::cout << "  cg=(" << cx_global << ", " << cy_global << ", " << cz_global << ")"
                                      << "  cl=(" << cx_local << ", " << cy_local << ", " << cz_local << ")";
                        }
                    } else {
                        // after conversion: mx,my,mz ARE already local centroids
                        const float cx_local  = c.mx;
                        const float cy_local  = c.my;
                        const float cz_local  = c.mz;

                        const float cx_global = cx_local + x_center;
                        const float cy_global = cy_local + y_center;
                        const float cz_global = cz_local + z_center;

                        std::cout << "  cl=(" << cx_local << ", " << cy_local << ", " << cz_local << ")"
                                  << "  cg=(" << cx_global << ", " << cy_global << ", " << cz_global << ")";
                    }
                }

                std::cout << "\n";
            }
        }
    }
}


int main(int argc, char* argv[]) {
    
    IRL::Data_gen gen;

    constexpr int stencil_size = 5;
    constexpr int include_moments = 1;
    constexpr bool include_Eigenvalues = false;

    for (int t = 0; t < 6; ++t) {
        // do it 10 times per state
        for (int i = 0; i < 10; ++i) {
            // print current state for debugging
            std::cout << "State " << t << std::endl;
            std::vector<float> state = gen.generateState(t,5,1,false,0.1,0.1,0.5,0.0,0.5,0.0,0.5,0.0,false);

            IRL::preprocess_stencil(state,
                                    stencil_size,
                                    1,                  // no_symmetries
                                    include_moments,
                                    include_Eigenvalues,
                                    0.0f,               // noise_stddev
                                    1e-12f);            // epsilon_connect
        }
    }

    

    /*

    constexpr int stencil_size = 5;
    constexpr int include_moments = 1;
    constexpr bool include_Eigenvalues = false;

    print_flat_state(state, "ORIGINAL STATE");
    print_decoded_state(state, stencil_size, include_moments, include_Eigenvalues, "ORIGINAL STATE", false);

    // For debugging, keep this as clean as possible:
    // - no noise
    // - minimal symmetry handling
    // - standard connectivity threshold
    IRL::preprocess_stencil(state,
                            stencil_size,
                            1,                  // no_symmetries
                            include_moments,
                            include_Eigenvalues,
                            0.0f,               // noise_stddev
                            1e-12f);            // epsilon_connect

    print_flat_state(state, "PREPROCESSED STATE");
    print_decoded_state(state, stencil_size, include_moments, include_Eigenvalues, "PREPROCESSED STATE", true);

    */

    //simulation_stability_analysis();

    /*
    int num_tests = 100;
    for (int i=0; i<4; i++) {
        float error = 0.0f;
        for (int j=0; j<num_tests; j++) {
            std::vector<float> errorVectorOut = gen.generateState(i,5,2,false,0.1,0.1,0.5,0.0,0.5,0.0,0.5,0.0,false,true);
            error += errorVectorOut[0];  // The first element is the error
        }
        error = error / static_cast<float>(num_tests);
        std::cout << "Class " << classNameFromType(i) << " average frobenius error: " << error << std::endl;
    }
    */

    /*
    For generating barycenter data
    const int stencil_size = 5;
    const int samples_per_class = 100;

    IRL::Data_gen gen;

    std::ofstream out("barycenters.csv");
    if (!out) {
        std::cerr << "Error: could not open barycenters.csv for writing\n";
        return 1;
    }

    // Optional header (Excel will happily read this)
    // For "two-line per sample" format, header is less meaningful; keep it simple:
    out << "x,y,z\n";

    // datapoint_type mapping based on your switch:
    // 0 -> default -> Paraboloid
    // 1 -> Cylinder
    // 2 -> Sphere
    // 3 -> Sheet
    const std::vector<int> types = {0, 1, 2, 3};

    for (int t : types) {
        const std::string cname = classNameFromType(t);

        for (int s = 0; s < samples_per_class; ++s) {
            // include_Moments must be >= 1 so we get vfrac + (mx,my,mz) per cell
            std::vector<double> flat = gen.generateState(
                t,
                stencil_size,
                1,
                false,
            );

            std::vector<CellData> stencil;
            unpackStencil(flat, stencil, stencil_size);

            // Compute global centroid of first moments (exactly like your snippet)
            double cx = 0.0, cy = 0.0, cz = 0.0;
            double Vsum = 0.0;

            for (const auto& c : stencil) {
                cx += c.mx;
                cy += c.my;
                cz += c.mz;
                Vsum += c.vfrac;
            }

            if (Vsum > 0.0) {
                cx /= Vsum;
                cy /= Vsum;
                cz /= Vsum;
            } else {
                cx = cy = cz = 0.0;
            }

            // --- Write "two lines per sample" (as you requested) ---
            out << cx << "," << cy << "," << cz << "," << cname << "\n";

            // --- Alternative (cleaner Excel): one row per sample ---
            // out << cx << "," << cy << "," << cz << "," << cname << "\n";
        }
    }

    std::cout << "Wrote barycenters.csv\n";
    return 0;
    */
    
    
    /*
    // Test e3nn

    torch::manual_seed(42);

    // ---- Parameters ----
    int stencil_size = 5;                // try 3, 5, or 7
    int hidden1 = 64, hidden2 = 32, hidden3 = 32, output = 4;
    int64_t N = int64_t(stencil_size) * stencil_size * stencil_size;
    torch::Device device(torch::kCPU);

    // ---- Create model ----
    IRL::e3nn model(stencil_size, hidden1, hidden2, hidden3, output, device);
    model->eval(); // inference mode

    // ---- Generate random stencil ----
    int B = 1;
    auto vof = torch::rand({B, N, 1});              // volume fractions [0,1]
    auto bary = torch::randn({B, N, 1, 3}) * 10.0;  // barycenters (random coords)

    std::cout << "Original stencil (before rotation):\n";
    print_stencil(vof, bary, stencil_size);
    std::cout << "----------------------------------------\n";

    // ---- Create rotated copy (90° around z) ----
    auto bary_rot = rotate90_z(bary).clone();

    std::cout << "Rotated stencil (90° about Z):\n";
    print_stencil(vof, bary_rot, stencil_size);
    std::cout << "----------------------------------------\n";

    // ---- Flatten to model input [B, 4*N] ----
    auto flatten_stencil = [&](const torch::Tensor& vof_, const torch::Tensor& mom_){
        auto vflat = vof_.reshape({B, N, 1});
        auto mflat = mom_.reshape({B, N, 3});
        auto combo = torch::cat({vflat, mflat}, -1);   // [B,N,4]
        return combo.reshape({B, 4*N});                // [B,4N]
    };

    auto x1 = flatten_stencil(vof, bary);
    auto x2 = flatten_stencil(vof, bary_rot);

    // ---- Forward both through model ----
    auto y1 = model->forward(x1);
    auto y2 = model->forward(x2);

    // ---- Print results ----
    std::cout << std::setprecision(6);
    std::cout << "Input 1 logits:\n" << y1 << "\n";
    std::cout << "Input 2 (rotated) logits:\n" << y2 << "\n";
    std::cout << "Absolute difference:\n" << (y1 - y2).abs() << "\n";
    std::cout << "Mean abs diff: " << (y1 - y2).abs().mean().item<double>() << std::endl;
    */



    

    return 0;
    
}