#include "irl/ml_classification/stencil_rotator.h"
#include "irl/ml_classification/data_gen.h"
#include "irl/ml_classification/ml_classifier.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace IRL {

using ScalarField3D = std::vector<std::vector<std::vector<double>>>;
using VectorField3D = std::vector<std::vector<std::vector<Eigen::Vector3d>>>;

static ScalarField3D makeScalarField3D(int N, double value = 0.0) {
    return ScalarField3D(
        N,
        std::vector<std::vector<double>>(
            N,
            std::vector<double>(N, value)
        )
    );
}

static VectorField3D makeVectorField3D(int N) {
    return VectorField3D(
        N,
        std::vector<std::vector<Eigen::Vector3d>>(
            N,
            std::vector<Eigen::Vector3d>(N, Eigen::Vector3d::Zero())
        )
    );
}

// For a domain of size N, coarse cells span [-N/2, N/2] with unit cell size.
// Cell i therefore has center at -N/2 + i + 0.5.
static double cellCenterCoordFromIndex(int idx, int N) {
    return -0.5 * static_cast<double>(N) + static_cast<double>(idx) + 0.5;
}

// Convert full 3D domain fields into the provided CellData layout.
static std::vector<CellData> deconstructDomainToCells(
    const ScalarField3D& vfrac,
    const VectorField3D& firstMoment,
    int N)
{
    std::vector<CellData> cells(N * N * N);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                const int idx = cellIndex(i, j, k, N);
                cells[idx].vfrac = static_cast<float>(vfrac[i][j][k]);
                cells[idx].mx    = static_cast<float>(firstMoment[i][j][k].x());
                cells[idx].my    = static_cast<float>(firstMoment[i][j][k].y());
                cells[idx].mz    = static_cast<float>(firstMoment[i][j][k].z());
            }
        }
    }

    return cells;
}

// Extract one local stencil of size stencil_size x stencil_size x stencil_size,
// flatten it, and convert the raw global first moments to the same local,
// center-cell-relative representation that your classifier pipeline expects.
static std::vector<float> extractFlattenedSubstencil(
    const std::vector<CellData>& domain_cells,
    int domain_size,
    int ic,
    int jc,
    int kc,
    int stencil_size,
    int include_moments)
{
    const int half = stencil_size / 2;

    const double cx_c = cellCenterCoordFromIndex(ic, domain_size);
    const double cy_c = cellCenterCoordFromIndex(jc, domain_size);
    const double cz_c = cellCenterCoordFromIndex(kc, domain_size);

    std::vector<double> flattened;
    flattened.reserve(
        stencil_size * stencil_size * stencil_size * perCellStride(include_moments) +
        globalTailStride(include_moments)
    );

    for (int di = -half; di <= half; ++di) {
        for (int dj = -half; dj <= half; ++dj) {
            for (int dk = -half; dk <= half; ++dk) {
                const int i = ic + di;
                const int j = jc + dj;
                const int k = kc + dk;

                const CellData& c = domain_cells[cellIndex(i, j, k, domain_size)];
                const double alpha = static_cast<double>(c.vfrac);

                // vfrac
                flattened.push_back(alpha);

                if (include_moments >= 1) {
                    // Domain generator stores raw global first moments:
                    //   m_global = alpha * bary_global
                    //
                    // Your classifier input uses first moments in the local
                    // center-cell-relative frame:
                    //   m_local = alpha * (bary_global - center_cell_center)
                    //
                    // Therefore:
                    //   m_local = m_global - alpha * center_cell_center
                    const double mx_local = static_cast<double>(c.mx) - alpha * cx_c;
                    const double my_local = static_cast<double>(c.my) - alpha * cy_c;
                    const double mz_local = static_cast<double>(c.mz) - alpha * cz_c;

                    flattened.push_back(mx_local);
                    flattened.push_back(my_local);
                    flattened.push_back(mz_local);
                }
            }
        }
    }

    if (include_moments >= 2) {
        Eigen::Matrix3d secondMoment =
            IRL::compute2ndMoment(flattened, stencil_size, /*from_ith_moment=*/1);

        flattened.push_back(secondMoment(0, 0)); // xx
        flattened.push_back(secondMoment(1, 1)); // yy
        flattened.push_back(secondMoment(2, 2)); // zz
        flattened.push_back(secondMoment(0, 1)); // xy
        flattened.push_back(secondMoment(0, 2)); // xz
        flattened.push_back(secondMoment(1, 2)); // yz
    }

    return std::vector<float>(flattened.begin(), flattened.end());
}

static const char* classNameFromId(int cls) {
    switch (cls) {
        case 0: return "Well-resolved";
        case 1: return "Ligament";
        case 2: return "Drop";
        case 3: return "Sheet";
        case 4: return "Ligament tip";
        case 5: return "Sheet edge";
        default: return "Unknown";
    }
}

/*
void shell_testcase(IRL::Classifier& classifier)
{
    constexpr int domain_size   = 36;
    constexpr int stencil_size  = 5;
    constexpr int half          = stencil_size / 2;

    // Match this to your trained model input:
    // 0 -> vfrac only
    // 1 -> vfrac + first moments
    // 2 -> vfrac + first moments + 2nd moments
    constexpr int include_moments = 1;

    constexpr bool cannonicalize_symmetries = true;
    constexpr float noise_stddev = 0.0f;
    constexpr double epsilon = 1.0e-10;

    // Bubble definition:
    // outer_diameter = bubble_diameter + thickness
    // inner_diameter = bubble_diameter - thickness
    constexpr double bubble_diameter = 30.0;
    constexpr double bubble_radius   = 0.5 * bubble_diameter;

    // Thickness sweep:
    // step_size, 2*step_size, ..., thickness_multiplier*step_size
    constexpr double step_size = 0.2;
    constexpr int thickness_multiplier = 10;

    using ScalarField3D = std::vector<std::vector<std::vector<double>>>;
    using VectorField3D = std::vector<std::vector<std::vector<Eigen::Vector3d>>>;

    auto makeScalarField3D = [](int N, double value = 0.0) {
        return ScalarField3D(
            N,
            std::vector<std::vector<double>>(
                N,
                std::vector<double>(N, value)
            )
        );
    };

    auto makeVectorField3D = [](int N) {
        return VectorField3D(
            N,
            std::vector<std::vector<Eigen::Vector3d>>(
                N,
                std::vector<Eigen::Vector3d>(N, Eigen::Vector3d::Zero())
            )
        );
    };

    auto cellCenterCoordFromIndex = [](int idx, int N) {
        return -0.5 * static_cast<double>(N) + static_cast<double>(idx) + 0.5;
    };

    auto deconstructDomainToCells =
        [&](const ScalarField3D& vfrac,
            const VectorField3D& firstMoment,
            int N)
    {
        std::vector<CellData> cells(N * N * N);

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k) {
                    const int idx = cellIndex(i, j, k, N);
                    cells[idx].vfrac = static_cast<float>(vfrac[i][j][k]);
                    cells[idx].mx    = static_cast<float>(firstMoment[i][j][k].x());
                    cells[idx].my    = static_cast<float>(firstMoment[i][j][k].y());
                    cells[idx].mz    = static_cast<float>(firstMoment[i][j][k].z());
                }
            }
        }

        return cells;
    };

    auto extractFlattenedSubstencil =
        [&](const std::vector<CellData>& domain_cells,
            int domain_size_local,
            int ic, int jc, int kc,
            int stencil_size_local,
            int include_moments_local)
    {
        const int half_local = stencil_size_local / 2;

        const double cx_c = cellCenterCoordFromIndex(ic, domain_size_local);
        const double cy_c = cellCenterCoordFromIndex(jc, domain_size_local);
        const double cz_c = cellCenterCoordFromIndex(kc, domain_size_local);

        std::vector<double> flattened;
        flattened.reserve(
            stencil_size_local * stencil_size_local * stencil_size_local *
            perCellStride(include_moments_local) +
            globalTailStride(include_moments_local)
        );

        for (int di = -half_local; di <= half_local; ++di) {
            for (int dj = -half_local; dj <= half_local; ++dj) {
                for (int dk = -half_local; dk <= half_local; ++dk) {
                    const int i = ic + di;
                    const int j = jc + dj;
                    const int k = kc + dk;

                    const CellData& c = domain_cells[cellIndex(i, j, k, domain_size_local)];
                    const double alpha = static_cast<double>(c.vfrac);

                    flattened.push_back(alpha);

                    if (include_moments_local >= 1) {
                        // Convert raw global first moments to center-cell-relative first moments
                        const double mx_local = static_cast<double>(c.mx) - alpha * cx_c;
                        const double my_local = static_cast<double>(c.my) - alpha * cy_c;
                        const double mz_local = static_cast<double>(c.mz) - alpha * cz_c;

                        flattened.push_back(mx_local);
                        flattened.push_back(my_local);
                        flattened.push_back(mz_local);
                    }
                }
            }
        }

        if (include_moments_local >= 2) {
            Eigen::Matrix3d secondMoment =
                IRL::compute2ndMoment(flattened, stencil_size_local, 1);

            flattened.push_back(secondMoment(0, 0)); // xx
            flattened.push_back(secondMoment(1, 1)); // yy
            flattened.push_back(secondMoment(2, 2)); // zz
            flattened.push_back(secondMoment(0, 1)); // xy
            flattened.push_back(secondMoment(0, 2)); // xz
            flattened.push_back(secondMoment(1, 2)); // yz
        }

        return std::vector<float>(flattened.begin(), flattened.end());
    };

    IRL::Data_gen data_gen;
    const Eigen::Vector3d origin(0.0, 0.0, 0.0);

    // Store results over all thicknesses
    std::vector<double> thicknesses;
    std::vector<int> well_resolved_counts;
    std::vector<int> ligament_counts;
    std::vector<int> drop_counts;
    std::vector<int> sheet_counts;
    std::vector<int> ligament_tip_counts;
    std::vector<int> sheet_edge_counts;
    std::vector<int> skipped_empty_or_full_counts;
    std::vector<int> classified_interfacial_counts;

    for (int m = 1; m <= thickness_multiplier; ++m) {

        const double thickness    = m * step_size;
        const double outer_radius = bubble_radius + 0.5 * thickness;
        const double inner_radius = bubble_radius - 0.5 * thickness;

        if (inner_radius <= 0.0) {
            std::cerr << "Skipping thickness " << thickness
                      << " because inner_radius <= 0.\n";
            continue;
        }

        if (outer_radius >= 0.5 * static_cast<double>(domain_size)) {
            std::cerr << "Skipping thickness " << thickness
                      << " because outer sphere reaches domain boundary.\n";
            continue;
        }

        auto outer_vfrac       = makeScalarField3D(domain_size, 0.0);
        auto inner_vfrac       = makeScalarField3D(domain_size, 0.0);
        auto shell_vfrac       = makeScalarField3D(domain_size, 0.0);

        auto outer_firstMoment = makeVectorField3D(domain_size);
        auto inner_firstMoment = makeVectorField3D(domain_size);
        auto shell_firstMoment = makeVectorField3D(domain_size);

        auto dummy_centroid = makeVectorField3D(domain_size);
        std::vector<IRL::ParaboloidParametrizedSurfaceOutput> dummy_surfaces;
        std::vector<double> coarse_coords(domain_size + 1, 0.0);

        data_gen.generateSpecificSphere(
            outer_vfrac,
            outer_firstMoment,
            dummy_centroid,
            dummy_surfaces,
            domain_size,
            origin,
            outer_radius,
            coarse_coords,
            false,
            nullptr
        );

        dummy_surfaces.clear();

        data_gen.generateSpecificSphere(
            inner_vfrac,
            inner_firstMoment,
            dummy_centroid,
            dummy_surfaces,
            domain_size,
            origin,
            inner_radius,
            coarse_coords,
            false,
            nullptr
        );

        // shell = outer sphere - inner sphere
        for (int i = 0; i < domain_size; ++i) {
            for (int j = 0; j < domain_size; ++j) {
                for (int k = 0; k < domain_size; ++k) {

                    double alpha = outer_vfrac[i][j][k] - inner_vfrac[i][j][k];
                    Eigen::Vector3d m1 =
                        outer_firstMoment[i][j][k] - inner_firstMoment[i][j][k];

                    if (alpha < epsilon) {
                        alpha = 0.0;
                        m1.setZero();
                    } else if (alpha > 1.0) {
                        alpha = 1.0;
                    }

                    shell_vfrac[i][j][k] = alpha;
                    shell_firstMoment[i][j][k] = m1;
                }
            }
        }

        std::vector<CellData> domain_cells =
            deconstructDomainToCells(shell_vfrac, shell_firstMoment, domain_size);

        std::array<int, 6> counts_interface_classes{};
        int skipped_empty_or_full = 0;
        int classified_interfacial = 0;

        for (int i = half; i < domain_size - half; ++i) {
            for (int j = half; j < domain_size - half; ++j) {
                for (int k = half; k < domain_size - half; ++k) {

                    const int center_id = cellIndex(i, j, k, domain_size);
                    const float alpha_center = domain_cells[center_id].vfrac;

                    const bool center_is_interface =
                        (alpha_center > epsilon && alpha_center < 1.0f - epsilon);

                    // Skip fully empty and fully full center cells
                    if (!center_is_interface) {
                        ++skipped_empty_or_full;
                        continue;
                    }

                    std::vector<float> flattened_state =
                        extractFlattenedSubstencil(
                            domain_cells,
                            domain_size,
                            i, j, k,
                            stencil_size,
                            include_moments
                        );

                    IRL::preprocess_stencil(
                        flattened_state,
                        stencil_size,
                        cannonicalize_symmetries,
                        include_moments,

                        noise_stddev,
                        1e-2f
                    );

                    std::vector<float> out_probs;
                    const int cls = classifier.classify(flattened_state, &out_probs);

                    ++classified_interfacial;

                    if (cls >= 0 && cls < 6) {
                        counts_interface_classes[cls]++;
                    }
                }
            }
        }

        thicknesses.push_back(thickness);
        well_resolved_counts.push_back(counts_interface_classes[0]);
        ligament_counts.push_back(counts_interface_classes[1]);
        drop_counts.push_back(counts_interface_classes[2]);
        sheet_counts.push_back(counts_interface_classes[3]);
        ligament_tip_counts.push_back(counts_interface_classes[4]);
        sheet_edge_counts.push_back(counts_interface_classes[5]);
        skipped_empty_or_full_counts.push_back(skipped_empty_or_full);
        classified_interfacial_counts.push_back(classified_interfacial);
    }

    // Print in Excel-friendly column format
    std::cout << "Shell testcase\n";
    std::cout << "domain_size              = " << domain_size << "\n";
    std::cout << "stencil_size             = " << stencil_size << "\n";
    std::cout << "step_size                = " << step_size << "\n";
    std::cout << "thickness_multiplier     = " << thickness_multiplier << "\n";
    std::cout << "bubble_diameter          = " << bubble_diameter << "\n";

    std::cout << "Thicknesses:\n";
    for (double x : thicknesses) {
        std::cout << x << "\n";
    }

    std::cout << "Well-resolved:\n";
    for (int x : well_resolved_counts) {
        std::cout << x << "\n";
    }

    std::cout << "Ligament:\n";
    for (int x : ligament_counts) {
        std::cout << x << "\n";
    }

    std::cout << "Drop:\n";
    for (int x : drop_counts) {
        std::cout << x << "\n";
    }

    std::cout << "Sheet:\n";
    for (int x : sheet_counts) {
        std::cout << x << "\n";
    }

    std::cout << "Ligament tip:\n";
    for (int x : ligament_tip_counts) {
        std::cout << x << "\n";
    }

    std::cout << "Sheet edge:\n";
    for (int x : sheet_edge_counts) {
        std::cout << x << "\n";
    }

    std::cout << "skipped empty/full cells:\n";
    for (int x : skipped_empty_or_full_counts) {
        std::cout << x << "\n";
    }

    std::cout << "classified interface:\n";
    for (int x : classified_interfacial_counts) {
        std::cout << x << "\n";
    }
}
*/

void shell_testcase(IRL::Classifier& classifier)
{
    constexpr int domain_size   = 36;
    constexpr int stencil_size  = 5;
    constexpr int half          = stencil_size / 2;

    // 0 -> vfrac only
    // 1 -> vfrac + first moments
    // 2 -> vfrac + first moments + global 2nd moments
    constexpr int include_moments = 1;

    // If true and include_moments >= 2, append 3 eigenvalues after the 6 second moments
    constexpr bool include_Eigenvalues = false;

    constexpr int no_symmetries = 48;
    constexpr float noise_stddev = 0.0f;
    constexpr float epsilon_connect = 1e-12f;
    constexpr double epsilon = 1.0e-10;

    constexpr double bubble_diameter = 30.0;
    constexpr double bubble_radius   = 0.5 * bubble_diameter;

    constexpr double step_size = 0.2;
    constexpr int thickness_multiplier = 10;

    using ScalarField3D = std::vector<std::vector<std::vector<double>>>;
    using VectorField3D = std::vector<std::vector<std::vector<Eigen::Vector3d>>>;

    auto makeScalarField3D = [](int N, double value = 0.0) {
        return ScalarField3D(
            N,
            std::vector<std::vector<double>>(
                N,
                std::vector<double>(N, value)
            )
        );
    };

    auto makeVectorField3D = [](int N) {
        return VectorField3D(
            N,
            std::vector<std::vector<Eigen::Vector3d>>(
                N,
                std::vector<Eigen::Vector3d>(N, Eigen::Vector3d::Zero())
            )
        );
    };

    auto cellCenterCoordFromIndex = [](int idx, int N) {
        return -0.5 * static_cast<double>(N) + static_cast<double>(idx) + 0.5;
    };

    auto deconstructDomainToCells =
        [&](const ScalarField3D& vfrac,
            const VectorField3D& firstMoment,
            int N)
    {
        std::vector<CellData> cells(N * N * N);

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k) {
                    const int idx = cellIndex(i, j, k, N);
                    cells[idx].vfrac = static_cast<float>(vfrac[i][j][k]);
                    cells[idx].mx    = static_cast<float>(firstMoment[i][j][k].x());
                    cells[idx].my    = static_cast<float>(firstMoment[i][j][k].y());
                    cells[idx].mz    = static_cast<float>(firstMoment[i][j][k].z());
                }
            }
        }

        return cells;
    };

    auto computeEigenvaluesFromFlattenedState =
        [](const std::vector<double>& flattened_state_local,
        int stencil_size_local,
        int include_moments_local,
        double machineZero = 1.0e-12) -> Eigenvalues
    {
        Eigen::Matrix3d I_tensor = IRL::computeInertiaTensor(
            flattened_state_local,
            stencil_size_local,
            include_moments_local,
            machineZero
        );

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(I_tensor);
        Eigenvalues eig{0.0f, 0.0f, 0.0f};

        if (solver.info() == Eigen::Success) {
            Eigen::Vector3d evals = solver.eigenvalues();

            // Sort descending: I1 >= I2 >= I3
            std::sort(evals.data(), evals.data() + 3, std::greater<double>());

            eig.lambda1 = static_cast<float>(evals[0]);
            eig.lambda2 = static_cast<float>(evals[1]);
            eig.lambda3 = static_cast<float>(evals[2]);
        }

        return eig;
    };

    auto extractFlattenedSubstencilWithHelpers =
        [&](const std::vector<CellData>& domain_cells,
            int domain_size_local,
            int ic, int jc, int kc,
            int stencil_size_local,
            int include_moments_local,
            bool include_Eigenvalues_local) -> std::vector<float>
    {
        const int half_local = stencil_size_local / 2;

        const double cx_c = cellCenterCoordFromIndex(ic, domain_size_local);
        const double cy_c = cellCenterCoordFromIndex(jc, domain_size_local);
        const double cz_c = cellCenterCoordFromIndex(kc, domain_size_local);

        std::vector<CellData> local_stencil(stencil_size_local * stencil_size_local * stencil_size_local);

        for (int di = -half_local; di <= half_local; ++di) {
            for (int dj = -half_local; dj <= half_local; ++dj) {
                for (int dk = -half_local; dk <= half_local; ++dk) {

                    const int gi = ic + di;
                    const int gj = jc + dj;
                    const int gk = kc + dk;

                    const CellData& src = domain_cells[cellIndex(gi, gj, gk, domain_size_local)];

                    const int si = di + half_local;
                    const int sj = dj + half_local;
                    const int sk = dk + half_local;

                    CellData dst{};
                    dst.vfrac = src.vfrac;

                    if (include_moments_local >= 1) {
                        const double alpha = static_cast<double>(src.vfrac);

                        // convert raw global first moments to center-cell-relative stencil moments
                        dst.mx = static_cast<float>(static_cast<double>(src.mx) - alpha * cx_c);
                        dst.my = static_cast<float>(static_cast<double>(src.my) - alpha * cy_c);
                        dst.mz = static_cast<float>(static_cast<double>(src.mz) - alpha * cz_c);
                    } else {
                        dst.mx = dst.my = dst.mz = 0.0f;
                    }

                    local_stencil[cellIndex(si, sj, sk, stencil_size_local)] = dst;
                }
            }
        }

        SecondMoments I{};
        SecondMoments* Ip = nullptr;
        
        Eigenvalues eig{};
        Eigenvalues* eigp = nullptr;

        // Build tmp_flat whenever we need either inertia or eigenvalues
        std::vector<double> tmp_flat;
        if (include_moments_local >= 2 || include_Eigenvalues_local) {
            tmp_flat.reserve(stencil_size_local * stencil_size_local * stencil_size_local * 4);

            for (int si = 0; si < stencil_size_local; ++si) {
                for (int sj = 0; sj < stencil_size_local; ++sj) {
                    for (int sk = 0; sk < stencil_size_local; ++sk) {
                        const CellData& c = local_stencil[cellIndex(si, sj, sk, stencil_size_local)];
                        tmp_flat.push_back(static_cast<double>(c.vfrac));
                        tmp_flat.push_back(static_cast<double>(c.mx));
                        tmp_flat.push_back(static_cast<double>(c.my));
                        tmp_flat.push_back(static_cast<double>(c.mz));
                    }
                }
            }
        }

        if (include_moments_local >= 2) {
            Eigen::Matrix3d M = IRL::compute2ndMoment(tmp_flat, stencil_size_local, /*from_ith_moment=*/1);

            I.Ixx = static_cast<float>(M(0,0));
            I.Iyy = static_cast<float>(M(1,1));
            I.Izz = static_cast<float>(M(2,2));
            I.Ixy = static_cast<float>(M(0,1));
            I.Ixz = static_cast<float>(M(0,2));
            I.Iyz = static_cast<float>(M(1,2));
            Ip = &I;
        }

        if (include_Eigenvalues_local) {
            eig = computeEigenvaluesFromFlattenedState(
                tmp_flat,
                stencil_size_local,
                include_moments_local,
                1.0e-12
            );
            eigp = &eig;
        }

        std::vector<float> flat;
        repackStencil(
            flat,
            local_stencil,
            Ip,
            eigp,
            stencil_size_local,
            include_moments_local,
            include_Eigenvalues_local
        );

        return flat;
    };

    auto printStencilSlices =
        [&](const std::vector<float>& flat_stencil,
            int stencil_size_local,
            int include_moments_local,
            bool include_Eigenvalues_local,
            int center_i, int center_j, int center_k,
            double thickness_local,
            int drop_counter)
    {
        std::vector<CellData> stencil_cells;
        SecondMoments I{};
        SecondMoments* Ip = (include_moments_local >= 2) ? &I : nullptr;
        Eigenvalues eig{};
        Eigenvalues* eigp = include_Eigenvalues_local ? &eig : nullptr;

        unpackStencil(
            flat_stencil,
            stencil_cells,
            Ip,
            eigp,
            stencil_size_local,
            include_moments_local,
            include_Eigenvalues_local
        );

        std::cout << "\n============================================================\n";
        std::cout << "DROP STENCIL DEBUG #" << drop_counter << "\n";
        std::cout << "thickness = " << thickness_local << "\n";
        std::cout << "domain center cell = [" << center_i << "," << center_j << "," << center_k << "]\n";
        std::cout << "stencil_size = " << stencil_size_local << "\n";

        for (int k = 0; k < stencil_size_local; ++k) {
            std::cout << "\n----- z-slice k = " << k << " -----\n";

            for (int j = 0; j < stencil_size_local; ++j) {
                for (int i = 0; i < stencil_size_local; ++i) {
                    const CellData& c = stencil_cells[cellIndex(i, j, k, stencil_size_local)];

                    float cx = 0.0f, cy = 0.0f, cz = 0.0f;
                    if (include_moments_local >= 1 && c.vfrac > 1.0e-12f) {
                        cx = c.mx / c.vfrac;
                        cy = c.my / c.vfrac;
                        cz = c.mz / c.vfrac;
                    }

                    std::cout << "[" << i << "," << j << "," << k << "] "
                            << c.vfrac << " "
                            << "(" << cx << ", " << cy << ", " << cz << ")";

                    if (i < stencil_size_local - 1) {
                        std::cout << "    ";
                    }
                }
                std::cout << "\n";
            }
        }

        std::cout << "Eigenvalues (if computed): "
                  << ((eigp) ? (std::to_string(eig.lambda1) + ", " +
                                std::to_string(eig.lambda2) + ", " +
                                std::to_string(eig.lambda3)) : "N/A")
                  << "\n";

        std::cout << "============================================================\n";
    };

    IRL::Data_gen data_gen;
    const Eigen::Vector3d origin(0.0, 0.0, 0.0);

    std::vector<double> thicknesses;
    std::vector<int> well_resolved_counts;
    std::vector<int> ligament_counts;
    std::vector<int> drop_counts;
    std::vector<int> sheet_counts;
    std::vector<int> ligament_tip_counts;
    std::vector<int> sheet_edge_counts;
    std::vector<int> skipped_empty_or_full_counts;
    std::vector<int> classified_interfacial_counts;

    for (int m = 1; m <= thickness_multiplier; ++m) {

        const double thickness    = m * step_size;
        const double outer_radius = bubble_radius + 0.5 * thickness;
        const double inner_radius = bubble_radius - 0.5 * thickness;

        if (inner_radius <= 0.0) {
            std::cerr << "Skipping thickness " << thickness
                      << " because inner_radius <= 0.\n";
            continue;
        }

        if (outer_radius >= 0.5 * static_cast<double>(domain_size)) {
            std::cerr << "Skipping thickness " << thickness
                      << " because outer sphere reaches domain boundary.\n";
            continue;
        }

        auto outer_vfrac       = makeScalarField3D(domain_size, 0.0);
        auto inner_vfrac       = makeScalarField3D(domain_size, 0.0);
        auto shell_vfrac       = makeScalarField3D(domain_size, 0.0);

        auto outer_firstMoment = makeVectorField3D(domain_size);
        auto inner_firstMoment = makeVectorField3D(domain_size);
        auto shell_firstMoment = makeVectorField3D(domain_size);

        auto dummy_centroid = makeVectorField3D(domain_size);
        std::vector<IRL::ParaboloidParametrizedSurfaceOutput> dummy_surfaces;
        std::vector<double> coarse_coords(domain_size + 1, 0.0);

        data_gen.generateSpecificSphere(
            outer_vfrac,
            outer_firstMoment,
            dummy_centroid,
            dummy_surfaces,
            domain_size,
            origin,
            outer_radius,
            coarse_coords,
            false,
            nullptr
        );

        dummy_surfaces.clear();

        data_gen.generateSpecificSphere(
            inner_vfrac,
            inner_firstMoment,
            dummy_centroid,
            dummy_surfaces,
            domain_size,
            origin,
            inner_radius,
            coarse_coords,
            false,
            nullptr
        );

        for (int i = 0; i < domain_size; ++i) {
            for (int j = 0; j < domain_size; ++j) {
                for (int k = 0; k < domain_size; ++k) {

                    double alpha = outer_vfrac[i][j][k] - inner_vfrac[i][j][k];
                    Eigen::Vector3d m1 =
                        outer_firstMoment[i][j][k] - inner_firstMoment[i][j][k];

                    if (alpha < epsilon) {
                        alpha = 0.0;
                        m1.setZero();
                    } else if (alpha > 1.0) {
                        alpha = 1.0;
                    }

                    shell_vfrac[i][j][k] = alpha;
                    shell_firstMoment[i][j][k] = m1;
                }
            }
        }

        std::vector<CellData> domain_cells =
            deconstructDomainToCells(shell_vfrac, shell_firstMoment, domain_size);

        std::array<int, 6> counts_interface_classes{};
        int skipped_empty_or_full = 0;
        int classified_interfacial = 0;
        int printed_drop_stencils = 0;

        for (int i = half; i < domain_size - half; ++i) {
            for (int j = half; j < domain_size - half; ++j) {
                for (int k = half; k < domain_size - half; ++k) {

                    const int center_id = cellIndex(i, j, k, domain_size);
                    const float alpha_center = domain_cells[center_id].vfrac;

                    bool center_is_interface =
                        (alpha_center > epsilon && alpha_center < 1.0f - epsilon);

                    if (!center_is_interface) {
                        ++skipped_empty_or_full;
                        continue;
                    }

                    std::vector<float> flattened_state =
                        extractFlattenedSubstencilWithHelpers(
                            domain_cells,
                            domain_size,
                            i, j, k,
                            stencil_size,
                            include_moments,
                            include_Eigenvalues
                        );

                    //std::cout << "Now in thickness " << thickness << ", cell (" << i << ", " << j << ", " << k << ")\n";
                    //std::cout << "Length of flattened state before: " << flattened_state.size() << "\n";

                    IRL::preprocess_stencil(
                        flattened_state,
                        stencil_size,
                        no_symmetries,
                        include_moments,
                        include_Eigenvalues,
                        noise_stddev,
                        epsilon_connect
                    );

                    //std::cout << "Length of flattened state: " << flattened_state.size() << "\n";

                    std::vector<float> out_probs;
                    const int cls = classifier.classify(flattened_state, &out_probs);

                    // Debug: for thickness == 0.4, print first 5 stencils classified as drops
                    if (cls == 1 && printed_drop_stencils < 5) {
                        std::cout << "Alpha center: " << alpha_center << std::endl;
                        ++printed_drop_stencils;

                        printStencilSlices(
                            flattened_state,
                            stencil_size,
                            include_moments,
                            include_Eigenvalues,
                            i, j, k,
                            thickness,
                            printed_drop_stencils
                        );
                    }

                    ++classified_interfacial;

                    if (cls >= 0 && cls < 6) {
                        counts_interface_classes[cls]++;
                    }
                }
            }
        }

        thicknesses.push_back(thickness);
        well_resolved_counts.push_back(counts_interface_classes[0]);
        ligament_counts.push_back(counts_interface_classes[1]);
        drop_counts.push_back(counts_interface_classes[2]);
        sheet_counts.push_back(counts_interface_classes[3]);
        ligament_tip_counts.push_back(counts_interface_classes[4]);
        sheet_edge_counts.push_back(counts_interface_classes[5]);
        skipped_empty_or_full_counts.push_back(skipped_empty_or_full);
        classified_interfacial_counts.push_back(classified_interfacial);
    }

    std::cout << "Shell testcase\n";
    std::cout << "domain_size              = " << domain_size << "\n";
    std::cout << "stencil_size             = " << stencil_size << "\n";
    std::cout << "step_size                = " << step_size << "\n";
    std::cout << "thickness_multiplier     = " << thickness_multiplier << "\n";
    std::cout << "bubble_diameter          = " << bubble_diameter << "\n";

    std::cout << "Thicknesses:\n";
    for (double x : thicknesses) std::cout << x << "\n";

    std::cout << "Well-resolved:\n";
    for (int x : well_resolved_counts) std::cout << x << "\n";

    std::cout << "Ligament:\n";
    for (int x : ligament_counts) std::cout << x << "\n";

    std::cout << "Drop:\n";
    for (int x : drop_counts) std::cout << x << "\n";

    std::cout << "Sheet:\n";
    for (int x : sheet_counts) std::cout << x << "\n";

    std::cout << "Ligament tip:\n";
    for (int x : ligament_tip_counts) std::cout << x << "\n";

    std::cout << "Sheet edge:\n";
    for (int x : sheet_edge_counts) std::cout << x << "\n";

    std::cout << "skipped empty/full cells:\n";
    for (int x : skipped_empty_or_full_counts) std::cout << x << "\n";

    std::cout << "classified interface:\n";
    for (int x : classified_interfacial_counts) std::cout << x << "\n";
}


} // namespace IRL