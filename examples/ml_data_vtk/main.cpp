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
#include <limits>

#include <limits>
#include <stdexcept>

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
    
    /*
    ml.updateDataParameters(no_batches, include_Moments,
                            paraboloid_coeff_stddev,
                            sheet_coeff_stddev,
                            max_sheet_thickness, sheet_thickness_stddev,
                            max_cylinder_radius, cylinder_radius_stddev, include_truncated_cylinder,
                            max_sphere_radius, sphere_radius_stddev, false);     
    */               
    
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

// Below for calculating surface areas of generated shapes:

double computeTotalSurfaceArea(const std::vector<float>& state,
                               int stencil_size,
                               int include_moments,
                               bool include_surface_area)
{
    if (!include_surface_area) {
        return 0.0;
    }

    const std::size_t n_cells =
        static_cast<std::size_t>(stencil_size) *
        static_cast<std::size_t>(stencil_size) *
        static_cast<std::size_t>(stencil_size);

    const std::size_t comps_per_cell =
        1 +                               // vfrac
        (include_moments >= 1 ? 3 : 0) + // mx, my, mz
        1;                               // area

    const std::size_t area_offset =
        1 + (include_moments >= 1 ? 3 : 0);

    const std::size_t required_size = n_cells * comps_per_cell;
    if (state.size() < required_size) {
        throw std::runtime_error("State vector is smaller than expected.");
    }

    double total_area = 0.0;
    for (std::size_t c = 0; c < n_cells; ++c) {
        total_area += static_cast<double>(state[c * comps_per_cell + area_offset]);
    }

    return total_area;
}

struct Stats {
    double mean = 0.0;
    double stddev = 0.0;
    double min = 0.0;
    double max = 0.0;
};

Stats computeStats(const std::vector<double>& values)
{
    Stats s;
    if (values.empty()) {
        return s;
    }

    double sum = 0.0;
    s.min = std::numeric_limits<double>::max();
    s.max = std::numeric_limits<double>::lowest();

    for (double v : values) {
        sum += v;
        if (v < s.min) s.min = v;
        if (v > s.max) s.max = v;
    }

    s.mean = sum / static_cast<double>(values.size());

    double sq_sum = 0.0;
    for (double v : values) {
        const double d = v - s.mean;
        sq_sum += d * d;
    }

    // sample standard deviation
    if (values.size() > 1) {
        s.stddev = std::sqrt(sq_sum / static_cast<double>(values.size() - 1));
    } else {
        s.stddev = 0.0;
    }

    return s;
}


int main(int argc, char* argv[]) {


int stencil_size = 5;

    //Data parameters
    int no_batches;
    int include_Moments = 1;
    bool include_Surface_Area = false;
    bool include_Eigenvalues = true;
    double paraboloid_coeff_stddev = 0.1;
    double hyperbolic_cylinder_opening_angle_stddev = 20; //degrees
    double sheet_coeff_stddev = 0.1;
    double sheet_thickness_stddev = 0.0;
    double cylinder_radius_stddev = 0.0;
    double radius_circle_min = 2.5;
    double radius_circle_max = 10.0;
    double sphere_radius_stddev = 0.0;
    double ellipsoid_subgrid_stddev = 0.7;
    double min_long_ellipsoid_axis = 3.0;
    double max_long_ellipsoid_axis = 5.0;
    bool exact_2nd_moment = false;  // enable calculation of exact 2nd moments for data generation
    bool visualize = false; // if true, print centroids and / or write surfaces
    double machineZero = 1e-12;
    double lower_limit_subgrid = machineZero;
    double upper_limit_subgrid = std::sqrt(3.0);
    double class0_max_characteristic = 2.5;
    float epsilon_connectivity = 1e-12f;

    // Net Parameters
    int input_size = stencil_size * stencil_size * stencil_size 
    * (include_Moments >= 1 ? /*(include_Surface_Area ? 5 : 4) ---> remove the 4 again if uncommenting ->*/ 4 : 1)  // 4 if include_Moments >= 1 because we have vfrac + (mx,my,mz) per cell, otherwise just vfrac
    + (include_Moments >= 2 ? 6 : 0)  // +6 if include_Moments >= 2 because we have (xx, yy, zz, xy, xz, yz) components of the 2nd moment tensor; otherwise none
    + (include_Eigenvalues ? 3 : 0);
    int hidden_size1 = 256;
    int hidden_size2 = 64;
    int hidden_size3 = 32;
    int output_size = 6; //CHANGED 4 to 6

    //Training parameters
    double learning_rate = 0.001; //was 0.01 for SGD optimizer
    int batch_size = 64;
    int max_epochs = 50;
    int reduce_lr_patience = 4;
    int early_stop_patience = 8;

    IRL::Data_gen gen; 
    
    gen.updateDataParameters(
            no_batches*batch_size,
            include_Moments,
            include_Surface_Area,
            include_Eigenvalues,
            paraboloid_coeff_stddev,
            hyperbolic_cylinder_opening_angle_stddev,
            sheet_coeff_stddev,
            sheet_thickness_stddev,
            cylinder_radius_stddev,
            radius_circle_min,
            radius_circle_max,
            sphere_radius_stddev,
            ellipsoid_subgrid_stddev,
            min_long_ellipsoid_axis,
            max_long_ellipsoid_axis,
            exact_2nd_moment,
            visualize,
            machineZero,
            lower_limit_subgrid,
            upper_limit_subgrid,
            class0_max_characteristic
        );   

    gen.setVisualize();

    while (true) {
        int input_class;

        std::cout << "\nEnter class index to generate, or -1 to quit: ";
        std::cin >> input_class;

        if (!std::cin) {
            std::cout << "Invalid input. Please enter an integer.\n";
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }

        if (input_class < 0) {
            std::cout << "Exiting.\n";
            break;
        }

        if (input_class > 5) {
            std::cout << "Invalid class. Please enter a class from 0 to 5.\n";
            continue;
        }

        std::vector<float> state = gen.generateState(false, input_class);

        std::cout << "Generated state for class " << input_class
                << " with size " << state.size() << ".\n";
    }

    // bool well_resolved = false;



    // for (int i = 1; i < 6; i++) {
    //     std::cout<<"Generating subtype "<< i << std::endl;
    //     for (int j = 0; j < 10; j++) {
    //         gen.generateState(
    //                 well_resolved,
    //                 i
    //             );
    //     }
    // }
    // well_resolved = true;

    // for (int i = 0; i < 6; i++) {
    //     std::cout<<"Generating subtype "<< i << std::endl;
    //     for (int j = 0; j < 100; j++) {
    //         gen.generateState(
    //                 well_resolved,
    //                 i
    //             );
    //     }
    // }
    //IRL::Data_gen gen; 
    //gen.setVisualize();
    //std::vector<float> state = gen.generateState(false, 14);  
    /*
    IRL::preprocess_stencil(state,
                            stencil_size,
                            1,                  // no_symmetries
                            include_moments,
                            true,               // include_Surface_Area
                            include_Eigenvalues,
                            0.0f,               // noise_stddev
                            1e-12f);            // epsilon_connect

    float total_surface_area = 0.0f; 
    for (size_t i = 0; i < state.size(); i += 5) { 
        total_surface_area += state[i + 4]; // 5th component is surface area if include_Surface_Area is true, otherwise it would be the 1st moment mx. 
    } 
    std::cout << "Total surface area in stencil: " << total_surface_area << std::endl;
    
    // For debugging, print the entire stencil 
    for (size_t i = 0; i < stencil_size*stencil_size*stencil_size*5; i += 5) { 
        std::cout << "Cell " << (i / 5) << ": vfrac=" << state[i] << ", mx=" << state[i + 1] << ", my=" << state[i + 2] << ", mz=" << state[i + 3] << ", area=" << state[i + 4] << "\n"; 
    }
    for (size_t i = stencil_size*stencil_size*stencil_size*5; i < state.size(); i++)
    {
        std::cout << state[i] << "\n";
    }
    */
    

    /*
    // Below for getting surface areas for all classes and comparing statistics:
    IRL::Data_gen gen;

    constexpr int stencil_size = 5;
    constexpr int include_moments = 1;
    constexpr bool include_surface_area = true;
    constexpr bool include_eigenvalues = false;

    constexpr int first_datapoint_type = 0;
    constexpr int last_datapoint_type = 11;
    constexpr int n_examples_per_class = 100;

    // generation parameters
    constexpr double paraboloid_coeff_stddev = 0.1;
    constexpr double sheet_coeff_stddev = 0.1;
    constexpr double max_sheet_thickness = 0.5;
    constexpr double sheet_thickness_stddev = 0.0;
    constexpr double max_cylinder_radius = 0.5;
    constexpr double cylinder_radius_stddev = 0.0;
    constexpr double max_sphere_radius = 0.5;
    constexpr double sphere_radius_stddev = 0.0;
    constexpr bool visualize = false;

    std::cout << std::fixed << std::setprecision(8);

    for (int datapoint_type = first_datapoint_type;
         datapoint_type <= last_datapoint_type;
         ++datapoint_type)
    {
        std::vector<double> total_areas;
        total_areas.reserve(n_examples_per_class);

        for (int sample_idx = 0; sample_idx < n_examples_per_class; ++sample_idx) {
            std::vector<float> state = gen.generateState(
                datapoint_type,
                stencil_size,
                include_moments,
                include_surface_area,
                include_eigenvalues,
                paraboloid_coeff_stddev,
                sheet_coeff_stddev,
                max_sheet_thickness,
                sheet_thickness_stddev,
                max_cylinder_radius,
                cylinder_radius_stddev,
                max_sphere_radius,
                sphere_radius_stddev,
                visualize
            );

            const double total_surface_area =
                computeTotalSurfaceArea(state, stencil_size, include_moments, include_surface_area);

            total_areas.push_back(total_surface_area);
        }

        const Stats stats = computeStats(total_areas);

        std::cout << "Datapoint type " << datapoint_type << '\n';
        std::cout << "  average total surface area = " << stats.mean << '\n';
        std::cout << "  stddev                     = " << stats.stddev << '\n';
        std::cout << "  min                        = " << stats.min << '\n';
        std::cout << "  max                        = " << stats.max << "\n\n";
    }
    */

    
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