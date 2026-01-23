#include <chrono>
#include <ctime>
#include <iostream>
#include <string>
#include <random>
#include <Eigen/Dense>
#include <vector>
#include <numeric>
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


int main(int argc, char* argv[]) {
    
    IRL::Data_gen gen;

    gen.generateState(4,5,1,false,0.1,0.1,0.5,0.0,0.5,0.0,0.5,0.0,true);

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