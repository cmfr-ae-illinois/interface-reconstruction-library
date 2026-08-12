#pragma once
#include <vector>
#include "common_functions.h"

namespace IRL {

// Abstract base classifier
class Classifier {
protected:
    int stencil_size;
    // Parameters for what information the classifier works with
    int include_Moments = 0;
    bool include_Surface_Area = false;
    bool include_Eigenvalues = false;
    // Parameters for how a stencil gets preprocessed before classification
    float epsilon_connect;
    int no_symmetries;
    float noise_stddev;

public:
    explicit Classifier(
        int stencil,
        float epsilon_connect_ = 1e-12f,
        int no_symmetries_ = 0,
        float noise_stddev_ = 0.0f)
        : stencil_size(stencil),
          epsilon_connect(epsilon_connect_),
          no_symmetries(no_symmetries_),
          noise_stddev(noise_stddev_)
    {}

    virtual ~Classifier() = default;

    int getStencilSize() const {
        return stencil_size;
    }

    int getIncludeMoments() const {
        return include_Moments;
    }

    bool getIncludeSurfaceArea() const {
        return include_Surface_Area;
    }

    bool getIncludeEigenvalues() const {
        return include_Eigenvalues;
    }

    float getEpsilonConnect() const {
        return epsilon_connect;
    }

    int getNoSymmetries() const {
        return no_symmetries;
    }

    float getNoiseStddev() const {
        return noise_stddev;
    }

    // Every classifier must implement classify
    virtual int classify(
        const std::vector<float>& flattened_state,
        std::vector<float>* out_probs = nullptr) = 0;

    inline void preprocess_stencil(std::vector<float>& flat_stencil){
        // Unpack
        std::vector<CellData> stencil;
        SecondMoments I{};
        SecondMoments* Ip = (include_Moments >= 2) ? &I : nullptr;
        Eigenvalues eig{};
        Eigenvalues* eigp = include_Eigenvalues ? &eig : nullptr;
        unpackStencil(flat_stencil, stencil, Ip, eigp, stencil_size, include_Moments, include_Surface_Area, include_Eigenvalues);

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

                            if (include_Moments >= 1) {
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

                            if (include_Moments >= 1) {
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
            repackStencil(flat_stencil, stencil, Ip, eigp, N, include_Moments, include_Surface_Area, include_Eigenvalues);
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
            repackStencil(flat_stencil, stencil, Ip, eigp, stencil_size, include_Moments, include_Surface_Area, include_Eigenvalues);
            return;
        }

        // Compute global centroid of first moments
        float cx = 0.0f, cy = 0.0f, cz = 0.0f;
        float Vsum = 0.0f;

        if (include_Moments >= 1) {
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
            // include_Moments == 0: approximate centroid from vfrac + cell centers
            approximateCentroidFromVfrac(stencil, stencil_size, cx, cy, cz, Vsum);
        }

        // Reflect into first octant: ensure cx,cy,cz >= 0 using x/y/z reflections.
        if (cx < 0.0) {
            reflect_x(stencil, stencil_size, include_Moments);
            if (include_Moments >= 2) reflect_x(I);
            cx = -cx;   // centroid x-component flips sign as well
        }
        if (cy < 0.0) {
            reflect_y(stencil, stencil_size, include_Moments);
            if (include_Moments >= 2) reflect_y(I);
            cy = -cy;
        }
        if (cz < 0.0) {
            reflect_z(stencil, stencil_size, include_Moments);
            if (include_Moments >= 2) reflect_z(I);
            cz = -cz;
        }

        if (no_symmetries <= 8) {
            repackStencil(flat_stencil, stencil, Ip, eigp, stencil_size, include_Moments, include_Surface_Area, include_Eigenvalues);
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
            permute_xyz(stencil, stencil_size, 1, include_Moments);
            if (include_Moments >= 2) permute_xyz(I, 1);
            // (cx,cy,cz) -> (cy,cz,cx)
            const float tmp = cx;
            cx = cy;
            cy = cz;
            cz = tmp;
        } else if (best == 2) {
            permute_xyz(stencil, stencil_size, 2, include_Moments);
            if (include_Moments >= 2) permute_xyz(I, 2);
            // (cx,cy,cz) -> (cz,cx,cy)
            const float tmp = cx;
            cx = cz;
            cz = cy;
            cy = tmp;
        }

        // If we only want 24 symmetries (rotations + octant choice), stop here.
        if (no_symmetries <= 24) {
            repackStencil(flat_stencil, stencil, Ip, eigp, stencil_size, include_Moments, include_Surface_Area, include_Eigenvalues);
            return;
        }

        // Mirror symmetry around diagonal x = y.
        if (cy > cx /*+ 1e-8*/) {
            // count small differences
            swap_xy(stencil, stencil_size, include_Moments);
            if (include_Moments >= 2) swap_xy(I);
        }
        // Pack back into flat storage
        repackStencil(flat_stencil, stencil, Ip, eigp, stencil_size, include_Moments, include_Surface_Area, include_Eigenvalues);
    }
};

} // namespace IRL