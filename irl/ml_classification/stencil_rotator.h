#pragma once
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <array>
#include <csignal> //debugging
#include <Eigen/Dense>

namespace IRL {

// Fowrard declarations:
inline Eigen::Matrix3d compute2ndMoment(const std::vector<float>& flat_stencil,
                                        int stencil_size,
                                        int from_ith_moment,
                                        bool include_Surface_Area,
                                        double machineZero,
                                        double cell_volume);

inline void appendInertiaEigenvalues(std::vector<float>& flattened_state,
                                     int stencil_size,
                                     int include_Moments,
                                     int from_ith_moment,
                                     bool include_Surface_Area,
                                     double machineZero);

// Per-cell data: volume fraction and first moments
struct CellData {
    float vfrac;
    float mx, my, mz;  // first moments in x, y, z
    float area;
};

// Global symmetric 2nd moment / inertia tensor packed as 6 unique components
struct SecondMoments {
    float Ixx, Iyy, Izz;
    float Ixy, Ixz, Iyz; // symmetric tensor => Iyx=Ixy, etc.
};

struct Eigenvalues {
    float lambda1, lambda2, lambda3;
};

// Convert (i,j,k) to flat cell index in the 1D array of cells
inline int cellIndex(int i, int j, int k, int N) {
    return (i * N * N + j * N + k);
}

// Strides in the flat vector depending on include_moments
inline int perCellStride(int include_moments, bool include_surface_area) {
    // include_moments:
    // 0 -> [vfrac]
    // 1 -> [vfrac, mx, my, mz]
    // include_surface_area: +1
    int stride = (include_moments >= 1) ? 4 : 1;
    stride += (include_surface_area) ? 1 : 0;
    return stride;
}

inline int globalTailStride(int include_moments, bool include_Eigenvalues = false) {
    // only include_moments==2 has the +6 tensor at the end, include_Eigenvalues adds +3 eigenvalues after that
    int tail = 0;
    if (include_moments == 2) tail += 6;
    if (include_Eigenvalues) tail += 3;
    return tail;
}

// Read one stencil from a flat array into CellData array (and optional global inertia)
static void unpackStencil(const std::vector<float>& flat,
                          std::vector<CellData>& stencil,
                          SecondMoments* I,
                          Eigenvalues* eigenvalues,
                          int stencil_size,
                          int include_moments,
                          bool include_Surface_Area = false,
                          bool include_Eigenvalues = false)
{
    const int nCells = stencil_size * stencil_size * stencil_size;
    const int stride = perCellStride(include_moments, include_Surface_Area);
    stencil.resize(nCells);

    for (int idx = 0; idx < nCells; ++idx) {
        stencil[idx].vfrac = flat[stride * idx + 0];
        if (stencil[idx].vfrac <= 1e-12f) {
            stencil[idx].vfrac = 0.0f;
        }

        if (include_moments >= 1) {
            stencil[idx].mx = flat[stride * idx + 1];
            stencil[idx].my = flat[stride * idx + 2];
            stencil[idx].mz = flat[stride * idx + 3];
        } else {
            stencil[idx].mx = stencil[idx].my = stencil[idx].mz = 0.0f;
        }

        if (include_Surface_Area) {
            stencil[idx].area = flat[stride * idx + 4];
        } else {
            stencil[idx].area = 0.0f;
        }
    }

    int base = stride * nCells;
    if (include_moments >= 2 && I) {
        int space_for_eigenvaues_in_back_of_vector = 0;
        if (include_Eigenvalues) {
            space_for_eigenvaues_in_back_of_vector = 3;
        }
        I->Ixx = flat[flat.size()-6-space_for_eigenvaues_in_back_of_vector];
        I->Iyy = flat[flat.size()-5-space_for_eigenvaues_in_back_of_vector];
        I->Izz = flat[flat.size()-4-space_for_eigenvaues_in_back_of_vector];
        I->Ixy = flat[flat.size()-3-space_for_eigenvaues_in_back_of_vector];
        I->Ixz = flat[flat.size()-2-space_for_eigenvaues_in_back_of_vector];
        I->Iyz = flat[flat.size()-1-space_for_eigenvaues_in_back_of_vector];
    }

    if (include_Eigenvalues) {
        // read eigenvalues into struct
        if (eigenvalues) {
            eigenvalues->lambda1 = flat[flat.size()-3];
            eigenvalues->lambda2 = flat[flat.size()-2];
            eigenvalues->lambda3 = flat[flat.size()-1];
        }
    }
}

inline void convert_to_local_centroids(std::vector<CellData>& stencil,
                                       int N,
                                       int include_moments)
{
    if (include_moments < 1) return;

    const float c0 = 0.5f * (static_cast<float>(N) - 1.0f);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                CellData& c = stencil[cellIndex(i, j, k, N)];

                const float x_center = static_cast<float>(i) - c0;
                const float y_center = static_cast<float>(j) - c0;
                const float z_center = static_cast<float>(k) - c0;

                if (c.vfrac <= 1e-12f) {
                    c.vfrac = 0.0f;
                    c.mx = 0.0f;
                    c.my = 0.0f;
                    c.mz = 0.0f;
                    c.area = 0.0f;
                    continue;
                }

                //const float cx_global = c.mx / IRL::safelyEpsilon(c.vfrac);
                //const float cy_global = c.my / IRL::safelyEpsilon(c.vfrac);
                //const float cz_global = c.mz / IRL::safelyEpsilon(c.vfrac);

                const float cx_global = c.mx / c.vfrac;
                const float cy_global = c.my / c.vfrac;
                const float cz_global = c.mz / c.vfrac;

                c.mx = cx_global;// - x_center;
                c.my = cy_global;// - y_center;
                c.mz = cz_global;// - z_center;
                /*
                // Check if conversion happened correctly (debugging), is it in -0.5, 0.5?
                if (std::abs(c.mx) > 0.5f || std::abs(c.my) > 0.5f || std::abs(c.mz) > 0.5f) {
                    std::cerr << "Warning: converted local centroid out of expected range [-0.5, 0.5] for cell (" << i << "," << j << "," << k << "): "
                              << "local centroid = (" << c.mx << ", " << c.my << ", " << c.mz << "), "
                              << "global centroid = (" << cx_global << ", " << cy_global << ", " << cz_global << ")\n";
                }
                */
            }
        }
    }
}


// Write CellData stencil back into flat array (and optional global inertia)
static void repackStencil(std::vector<float>& flat,
                          const std::vector<CellData>& stencil,
                          const SecondMoments* I,
                          const Eigenvalues* eigenvalues,
                          int stencil_size,
                          int include_moments,
                          bool include_Surface_Area = false,
                          //bool include_Total_Surface_Area = false,
                          bool include_Eigenvalues = false)
{
    std::vector<CellData> packed_stencil = stencil;
    //convert_to_local_centroids(packed_stencil, N, include_moments);

    // Normalize surface areas by the liquid cell volume ^2/3
    /*
    for (auto& c : packed_stencil) {
        //c.area /= std::pow(c.vfrac, 2.0f / 3.0f) + 1e-12f; // add small epsilon to avoid division by zero
        // get thickness parameter
        c.area = c.vfrac / (c.area + 1e-12f);
    }
    */    
    

    const int nCells = stencil_size * stencil_size * stencil_size;
    int stride = perCellStride(include_moments, include_Surface_Area);
    const int tail   = globalTailStride(include_moments, include_Eigenvalues);

    flat.resize(stride * nCells); //This is the size not including stencil wide 2nd moments or eigenvalues, those get appended later
    /* omit surface area if total surface area is wanted
    if (include_Surface_Area && include_Total_Surface_Area) {
        flat.resize((stride-1) * nCells + tail + 1); // +1 for total surface area at the end
    }
    */

    for (int idx = 0; idx < nCells; ++idx) {
        flat[stride * idx + 0] = packed_stencil[idx].vfrac;

        if (include_moments >= 1) {
            flat[stride * idx + 1] = packed_stencil[idx].mx;
            flat[stride * idx + 2] = packed_stencil[idx].my;
            flat[stride * idx + 3] = packed_stencil[idx].mz;
        }

        if (include_Surface_Area) {
            flat[stride * idx + 4] = packed_stencil[idx].area;
        }
    }
    /* Old: Use previouly calculated 2nd moments and eigenvalues
    int base = stride * nCells;
    if (include_moments >= 2 && I) {
        flat[base + 0] = I->Ixx;
        flat[base + 1] = I->Iyy;
        flat[base + 2] = I->Izz;
        flat[base + 3] = I->Ixy;
        flat[base + 4] = I->Ixz;
        flat[base + 5] = I->Iyz;
    }

    if (include_Eigenvalues) {
        if (include_moments >= 2 && I) {
            base += 6;
        }
        if (eigenvalues) {
            flat[base + 0] = eigenvalues->lambda1;
            flat[base + 1] = eigenvalues->lambda2;
            flat[base + 2] = eigenvalues->lambda3;
        }
    }
    */

    //New: Re-calculate 2nd moments and eigenvalues
    if (include_moments >= 2) {
        // Append approximate 2nd moments to flattened_state
        Eigen::Matrix3d approxSecondMoment = IRL::compute2ndMoment(flat, stencil_size, 1, include_Surface_Area, 1e-12, 1.0);
        flat.push_back(static_cast<float>(approxSecondMoment(0, 0))); // Ixx
        flat.push_back(static_cast<float>(approxSecondMoment(1, 1))); // Iyy
        flat.push_back(static_cast<float>(approxSecondMoment(2, 2))); // Izz
        flat.push_back(static_cast<float>(approxSecondMoment(0, 1))); // Ixy
        flat.push_back(static_cast<float>(approxSecondMoment(0, 2))); // Ixz
        flat.push_back(static_cast<float>(approxSecondMoment(1, 2))); // Iyz
    }
    if (include_Eigenvalues) {
        IRL::appendInertiaEigenvalues(flat, stencil_size, include_moments, 1, include_Surface_Area, 1e-12);
    }

    /* Below was a quick test to see if total surface area is benefitial, did not seem to help, so no better embedding into the code has happened
    if (include_Total_Surface_Area) {
        float total_area = 0.0f;
        for (const auto& c : packed_stencil) {
            total_area += c.area;
        }
        flat.back() = total_area;
    }

    std::cout << "Repacked stencil length = " << flat.size() << "\n";
    */
}

// Symmetry helpers


// Reflections into the first octant

// Reflect across plane x -> -x  (mirror in x, i <-> N-1-i; mx -> -mx)
inline void reflect_x(std::vector<CellData>& S, int N, int include_moments) {
    const std::vector<CellData> C = S;  // copy source
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k) {
                CellData d = C[cellIndex(N - 1 - i, j, k, N)];
                // First moment transforms as a vector under reflection
                if (include_moments >= 1) {
                    d.mx = -d.mx;       // x component flips sign
                    // y,z unchanged
                }
                S[cellIndex(i, j, k, N)] = d;
            }
}

// Reflect across plane y -> -y  (mirror in y, j <-> N-1-j; my -> -my)
inline void reflect_y(std::vector<CellData>& S, int N, int include_moments) {
    const std::vector<CellData> C = S;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k) {
                CellData d = C[cellIndex(i, N - 1 - j, k, N)];
                if (include_moments >= 1) {
                    d.my = -d.my;       // y component flips sign
                }
                S[cellIndex(i, j, k, N)] = d;
            }
}

// Reflect across plane z -> -z  (mirror in z, k <-> N-1-k; mz -> -mz)
inline void reflect_z(std::vector<CellData>& S, int N, int include_moments) {
    const std::vector<CellData> C = S;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k) {
                CellData d = C[cellIndex(i, j, N - 1 - k, N)];
                if (include_moments >= 1) {
                    d.mz = -d.mz;       // z component flips sign
                }
                S[cellIndex(i, j, k, N)] = d;
            }
}

// Transform global 2nd moments under reflections (I' = R I R^T)
inline void reflect_x(SecondMoments& I) {
    I.Ixy = -I.Ixy;
    I.Ixz = -I.Ixz;
    // Ixx,Iyy,Izz,Iyz unchanged
}
inline void reflect_y(SecondMoments& I) {
    I.Ixy = -I.Ixy;
    I.Iyz = -I.Iyz;
}
inline void reflect_z(SecondMoments& I) {
    I.Ixz = -I.Ixz;
    I.Iyz = -I.Iyz;
}

// 120° rotations about the (1,1,1) diagonal

inline void permute_xyz(std::vector<CellData>& S, int N, int type, int include_moments) {
    const std::vector<CellData> C = S;
    if (type == 1) {
        // (x,y,z)->(y,z,x)
        for (int i=0;i<N;++i)
        for (int j=0;j<N;++j)
            for (int k=0;k<N;++k) {
                // new(i,j,k) = old(k,i,j)
                CellData d = C[cellIndex(k, i, j, N)];

                if (include_moments >= 1) {
                    const float mx_old=d.mx, my_old=d.my, mz_old=d.mz;
                    d.mx = my_old;
                    d.my = mz_old;
                    d.mz = mx_old;
                }

                S[cellIndex(i, j, k, N)] = d;
            }
    } else if (type == 2) {
        // (x,y,z)->(z,x,y)
        for (int i=0;i<N;++i)
        for (int j=0;j<N;++j)
            for (int k=0;k<N;++k) {
                // new(i,j,k) = old(j,k,i)
                CellData d = C[cellIndex(j, k, i, N)];

                if (include_moments >= 1) {
                    const float mx_old=d.mx, my_old=d.my, mz_old=d.mz;
                    d.mx = mz_old;
                    d.my = mx_old;
                    d.mz = my_old;
                }

                S[cellIndex(i, j, k, N)] = d;
            }
    }
}

// Transform global 2nd moments under axis permutations (I' = P I P^T)
inline void permute_xyz(SecondMoments& I, int type) {
    if (type == 1) {
        // (x,y,z)->(y,z,x); p(x)=y, p(y)=z, p(z)=x
        SecondMoments O = I;
        I.Ixx = O.Iyy;
        I.Iyy = O.Izz;
        I.Izz = O.Ixx;
        I.Ixy = O.Iyz;
        I.Ixz = O.Ixy;
        I.Iyz = O.Ixz;
    } else if (type == 2) {
        // (x,y,z)->(z,x,y); p(x)=z, p(y)=x, p(z)=y
        SecondMoments O = I;
        I.Ixx = O.Izz;
        I.Iyy = O.Ixx;
        I.Izz = O.Iyy;
        I.Ixy = O.Ixz;
        I.Ixz = O.Iyz;
        I.Iyz = O.Ixy;
    }
}

// Reflection across the diagonal plane x = y

inline void swap_xy(std::vector<CellData>& S, int N, int include_moments) {
    const std::vector<CellData> C = S;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k) {
                CellData d = C[cellIndex(j, i, k, N)];
                if (include_moments >= 1) {
                    const float mx_old = d.mx, my_old = d.my;
                    d.mx = my_old;  // x' component = old y
                    d.my = mx_old;  // y' component = old x
                    // z unchanged
                }
                S[cellIndex(i, j, k, N)] = d;
            }
}

// Transform global 2nd moments under swap x<->y
inline void swap_xy(SecondMoments& I) {
    SecondMoments O = I;
    I.Ixx = O.Iyy;
    I.Iyy = O.Ixx;
    // Izz unchanged
    I.Ixy = O.Ixy;
    I.Ixz = O.Iyz;
    I.Iyz = O.Ixz;
}

// Approximate centroid from volume fractions assuming per-cell centroid at cell center
// This is only used when include_moments == 0.
static void approximateCentroidFromVfrac(const std::vector<CellData>& stencil,
                                         int N,
                                         float& cx, float& cy, float& cz,
                                         float& Vsum)
{
    cx = cy = cz = 0.0f;
    Vsum = 0.0f;

    const float c0 = 0.5f * (static_cast<float>(N) - 1.0f); // center index

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k) {
                const CellData& c = stencil[cellIndex(i, j, k, N)];
                const float v = c.vfrac;
                if (v == 0.0f) continue;

                // cell-center coordinates relative to stencil center
                const float x = static_cast<float>(i) - c0;
                const float y = static_cast<float>(j) - c0;
                const float z = static_cast<float>(k) - c0;

                cx += v * x;
                cy += v * y;
                cz += v * z;
                Vsum += v;
            }

    if (Vsum > 0.0f) {
        cx /= Vsum;
        cy /= Vsum;
        cz /= Vsum;
    } else {
        cx = cy = cz = 0.0f;
    }
}

inline float clampf(float x, float lo, float hi) {
    return std::max(lo, std::min(x, hi));
}

inline std::mt19937& globalNoiseRng() {
    static thread_local std::mt19937 rng(std::random_device{}());
    return rng;
}

inline int sampleNoiseCategory() {
    // 0: none (50%)
    // 1: mild (30%)
    // 2: moderate (15%)
    // 3: strong (5%)
    static thread_local std::discrete_distribution<int> dist{25, 25, 25, 25};
    return dist(globalNoiseRng());
}

inline float sampleGaussian(float sigma) {
    if (sigma <= 0.0f) return 0.0f;
    std::normal_distribution<float> dist(0.0f, sigma);
    return dist(globalNoiseRng());
}

inline void preprocess_stencil(std::vector<float>& flat_stencil,
                           int stencil_size, int no_symmetries, int include_moments = 1, bool include_Surface_Area = false, 
                           bool include_Eigenvalues = false, float noise_stddev = 0.0f, float epsilon_connect = 1e-12f)
{
    // Unpack
    std::vector<CellData> stencil;
    SecondMoments I{};
    SecondMoments* Ip = (include_moments >= 2) ? &I : nullptr;
    Eigenvalues eig{};
    Eigenvalues* eigp = include_Eigenvalues ? &eig : nullptr;
    unpackStencil(flat_stencil, stencil, Ip, eigp, stencil_size, include_moments, include_Surface_Area, include_Eigenvalues);

    // Noise
    if (noise_stddev > 1e-10f) {
        const float epsilon_noise = 1.0e-6f;
        const float c0 = 0.5f * (static_cast<float>(stencil_size) - 1.0f);

        // Probability to inject small noise into an empty cell
        const float empty_cell_probability = 1.0f / 6.0f;

        // Keep synthetic empty-cell noise small
        const float empty_vfrac_scale    = 0.25f;
        const float empty_centroid_scale = 0.25f;

        std::bernoulli_distribution add_empty_noise(empty_cell_probability);

        auto shapeFactor = [](float v) -> float {
            // max at v=0.5, zero at v=0 and v=1
            return 4.0f * v * (1.0f - v);
        };

        for (int i = 0; i < stencil_size; ++i) {
            for (int j = 0; j < stencil_size; ++j) {
                for (int k = 0; k < stencil_size; ++k) {

                    auto& c = stencil[cellIndex(i, j, k, stencil_size)];
                    const float v_old = c.vfrac;

                    // cell-center coordinates in global stencil frame
                    const float x_center = static_cast<float>(i) - c0;
                    const float y_center = static_cast<float>(j) - c0;
                    const float z_center = static_cast<float>(k) - c0;

                    // Case 1: cell already has liquid
                    if (v_old > epsilon_noise) {

                        // Multiply the sampled Gaussian by f(v)
                        const float v_new = clampf(
                            v_old + sampleGaussian(noise_stddev) * shapeFactor(v_old),
                            1e-12f,
                            1.0f - 1e-12f
                        );

                        c.vfrac = v_new;

                        if (include_moments >= 1) {
                            // recover centroid in global stencil coordinates
                            float cx_cell = c.mx / v_old;
                            float cy_cell = c.my / v_old;
                            float cz_cell = c.mz / v_old;

                            // add centroid noise directly, no clipping to cell bounds
                            cx_cell += sampleGaussian(noise_stddev) * (1-v_new);
                            cy_cell += sampleGaussian(noise_stddev) * (1-v_new);
                            cz_cell += sampleGaussian(noise_stddev) * (1-v_new);

                            c.mx = v_new * cx_cell;
                            c.my = v_new * cy_cell;
                            c.mz = v_new * cz_cell;
                        }

                        continue;
                    }
                    
                    // Case 2: currently empty cell
                    if (v_old <= 1e-12f && add_empty_noise(globalNoiseRng())) {

                        // choose vfrac uniformly in [1e-12, 1-1e-12]
                        static thread_local std::uniform_real_distribution<float> uniform_vfrac(
                            1e-12f, 1.0f - 1e-12f
                        );
                        const float v_new = uniform_vfrac(globalNoiseRng());

                        c.vfrac = v_new;

                        if (include_moments >= 1) {
                            // assume old centroid is at the cell center
                            float cx_cell = x_center;
                            float cy_cell = y_center;
                            float cz_cell = z_center;

                            // add centroid noise exactly as for non-empty cells
                            // add centroid noise directly, no clipping to cell bounds
                            cx_cell += sampleGaussian(noise_stddev) * (1-v_new);
                            cy_cell += sampleGaussian(noise_stddev) * (1-v_new);
                            cz_cell += sampleGaussian(noise_stddev) * (1-v_new);

                            c.mx = v_new * cx_cell;
                            c.my = v_new * cy_cell;
                            c.mz = v_new * cz_cell;
                        } else {
                            c.mx = c.my = c.mz = 0.0f;
                        }
                    }
                    else {
                        c.vfrac = 0.0f;
                        c.mx = c.my = c.mz = 0.0f;
                    }
                    
                }
            }
        }
    }

    const int N = stencil_size;
    const int c = N / 2;
    
    // Perform BFS Search
    //const float eps = 1e-2f; // "positive" threshold for connectivity, value here for simulation purposes but can be made smaller for sythetic data 

    // helper: bounds check
    auto in_bounds = [&](int i, int j, int k) {
        return (i >= 0 && i < N &&
                j >= 0 && j < N &&
                k >= 0 && k < N);
    };

    // 6-neighborhood
    static const int n6[6][3] = {
        {+1, 0, 0}, {-1, 0, 0},
        { 0,+1, 0}, { 0,-1, 0},
        { 0, 0,+1}, { 0, 0,-1}
    };

    // If center is not positive → zero everything and return
    if (stencil[cellIndex(c,c,c,N)].vfrac <= epsilon_connect) {
        /*
        for (auto& cell : stencil) {
            cell.vfrac = 0.0f;
            cell.mx = cell.my = cell.mz = 0.0f;
        }
        */
        repackStencil(flat_stencil, stencil, Ip, eigp, N, include_moments, include_Surface_Area, include_Eigenvalues);
        return;
    }

    // visited mask
    std::vector<uint8_t> visited(N * N * N, 0);

    // BFS queue
    std::vector<std::array<int,3>> q;
    q.reserve(N * N * N);

    // seed
    visited[cellIndex(c,c,c,N)] = 1;
    q.push_back({c,c,c});

    // flood fill
    for (size_t qi = 0; qi < q.size(); ++qi) {
        int i = q[qi][0];
        int j = q[qi][1];
        int k = q[qi][2];

        for (int n = 0; n < 6; ++n) {
            int ni = i + n6[n][0];
            int nj = j + n6[n][1];
            int nk = k + n6[n][2];

            if (!in_bounds(ni,nj,nk)) continue;

            int id = cellIndex(ni,nj,nk,N);
            if (visited[id]) continue;

            if (stencil[id].vfrac > epsilon_connect) {
                visited[id] = 1;
                q.push_back({ni,nj,nk});
            }
        }
    }

    // zero everything not connected
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {

                int id = cellIndex(i,j,k,N);

                if (visited[id]) continue;

                if (stencil[id].vfrac > 0.0f) {
                    stencil[id].vfrac = 0.0f;
                    stencil[id].mx = 0.0f;
                    stencil[id].my = 0.0f;
                    stencil[id].mz = 0.0f;
                    stencil[id].area = 0.0f;
                }
            }
        }
    }
    

    if (no_symmetries < 8) {
        repackStencil(flat_stencil, stencil, Ip, eigp, stencil_size, include_moments, include_Surface_Area, include_Eigenvalues);
        return;
    }

    // Compute global centroid of first moments
    float cx = 0.0f, cy = 0.0f, cz = 0.0f;
    float Vsum = 0.0f;

    if (include_moments >= 1) {
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
    } else {
        // include_moments == 0: approximate centroid from vfrac + cell centers
        approximateCentroidFromVfrac(stencil, stencil_size, cx, cy, cz, Vsum);
    }

    // Reflect into first octant: ensure cx,cy,cz >= 0 using x/y/z reflections.
    if (cx < 0.0) {
        reflect_x(stencil, stencil_size, include_moments);
        if (include_moments >= 2) reflect_x(I);
        cx = -cx;   // centroid x-component flips sign as well
    }
    if (cy < 0.0) {
        reflect_y(stencil, stencil_size, include_moments);
        if (include_moments >= 2) reflect_y(I);
        cy = -cy;
    }
    if (cz < 0.0) {
        reflect_z(stencil, stencil_size, include_moments);
        if (include_moments >= 2) reflect_z(I);
        cz = -cz;
    }

    if (no_symmetries <= 8) {
        repackStencil(flat_stencil, stencil, Ip, eigp, stencil_size, include_moments, include_Surface_Area, include_Eigenvalues);
        return;
    }

    // 120° rotations about (1,1,1) diagonal.

    const float z0 = cz;
    const float z1 = cx;
    const float z2 = cy;

    int best = 0;
    float bestz = z0;
    if (z1 > bestz) { best = 1; bestz = z1; }
    if (z2 > bestz) { best = 2; bestz = z2; }

    if (best == 1) {
        permute_xyz(stencil, stencil_size, 1, include_moments);
        if (include_moments >= 2) permute_xyz(I, 1);
        // (cx,cy,cz) -> (cy,cz,cx)
        const float tmp = cx;
        cx = cy;
        cy = cz;
        cz = tmp;
    } else if (best == 2) {
        permute_xyz(stencil, stencil_size, 2, include_moments);
        if (include_moments >= 2) permute_xyz(I, 2);
        // (cx,cy,cz) -> (cz,cx,cy)
        const float tmp = cx;
        cx = cz;
        cz = cy;
        cy = tmp;
    }

    // If we only want 24 symmetries (rotations + octant choice), stop here.
    if (no_symmetries <= 24) {
        repackStencil(flat_stencil, stencil, Ip, eigp, stencil_size, include_moments, include_Surface_Area, include_Eigenvalues);
        return;
    }

    // Mirror symmetry around diagonal x = y.
    if (cy > cx /*+ 1e-8*/) {
        // count small differences
        swap_xy(stencil, stencil_size, include_moments);
        if (include_moments >= 2) swap_xy(I);
    }
    // Pack back into flat storage
    repackStencil(flat_stencil, stencil, Ip, eigp, stencil_size, include_moments, include_Surface_Area, include_Eigenvalues);

    //m1.x / IRL::safelyEpsilon(vfrac);, recheck if -0.5,0.5 after
    // check in sheet if it is well resolved by comparing V1 and V2
}

} // namespace IRL
#include "irl/ml_classification/inertia_calc.h"