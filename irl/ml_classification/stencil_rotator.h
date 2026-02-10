#pragma once
#include <vector>
#include <cmath>

namespace IRL {

// Per-cell data: volume fraction and first moments
struct CellData {
    float vfrac;
    float mx, my, mz;  // first moments in x, y, z
};

// Global symmetric 2nd moment / inertia tensor packed as 6 unique components
struct SecondMoments {
    float Ixx, Iyy, Izz;
    float Ixy, Ixz, Iyz; // symmetric tensor => Iyx=Ixy, etc.
};

// Convert (i,j,k) to flat cell index in the 1D array of cells
inline int cellIndex(int i, int j, int k, int N) {
    return (i * N * N + j * N + k);
}

// Strides in the flat vector depending on include_moments
inline int perCellStride(int include_moments) {
    // include_moments:
    // 0 -> [vfrac]
    // 1 -> [vfrac, mx, my, mz]
    // 2 -> [vfrac, mx, my, mz] + 6 global inertia at end
    return (include_moments >= 1) ? 4 : 1;
}

inline int globalTailStride(int include_moments) {
    // only include_moments==2 has the +6 tensor at the end
    return (include_moments >= 2) ? 6 : 0;
}

// Read one stencil from a flat array into CellData array (and optional global inertia)
static void unpackStencil(const std::vector<float>& flat,
                          std::vector<CellData>& stencil,
                          SecondMoments* I,
                          int N,
                          int include_moments)
{
    const int nCells = N * N * N;
    const int stride = perCellStride(include_moments);
    stencil.resize(nCells);

    for (int idx = 0; idx < nCells; ++idx) {
        stencil[idx].vfrac = flat[stride * idx + 0];

        if (include_moments >= 1) {
            stencil[idx].mx = flat[stride * idx + 1];
            stencil[idx].my = flat[stride * idx + 2];
            stencil[idx].mz = flat[stride * idx + 3];
        } else {
            stencil[idx].mx = stencil[idx].my = stencil[idx].mz = 0.0f;
        }
    }

    if (include_moments >= 2 && I) {
        const int base = stride * nCells;
        I->Ixx = flat[base + 0];
        I->Iyy = flat[base + 1];
        I->Izz = flat[base + 2];
        I->Ixy = flat[base + 3];
        I->Ixz = flat[base + 4];
        I->Iyz = flat[base + 5];
    }
}

// Write CellData stencil back into flat array (and optional global inertia)
static void repackStencil(std::vector<float>& flat,
                          const std::vector<CellData>& stencil,
                          const SecondMoments* I,
                          int N,
                          int include_moments)
{
    const int nCells = N * N * N;
    const int stride = perCellStride(include_moments);
    const int tail   = globalTailStride(include_moments);

    flat.resize(stride * nCells + tail);

    for (int idx = 0; idx < nCells; ++idx) {
        flat[stride * idx + 0] = stencil[idx].vfrac;

        if (include_moments >= 1) {
            flat[stride * idx + 1] = stencil[idx].mx;
            flat[stride * idx + 2] = stencil[idx].my;
            flat[stride * idx + 3] = stencil[idx].mz;
        }
    }

    if (include_moments >= 2 && I) {
        const int base = stride * nCells;
        flat[base + 0] = I->Ixx;
        flat[base + 1] = I->Iyy;
        flat[base + 2] = I->Izz;
        flat[base + 3] = I->Ixy;
        flat[base + 4] = I->Ixz;
        flat[base + 5] = I->Iyz;
    }
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

inline void rotate_stencil(std::vector<float>& flat_stencil,
                           int stencil_size, int no_symmetries, int include_moments = 1)
{
    if (no_symmetries < 8) {
        return; // no rotation
    }

    // Unpack the flat vector into structured cells
    std::vector<CellData> stencil;
    SecondMoments I{};
    SecondMoments* Ip = (include_moments >= 2) ? &I : nullptr;
    unpackStencil(flat_stencil, stencil, Ip, stencil_size, include_moments);

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
        repackStencil(flat_stencil, stencil, Ip, stencil_size, include_moments);
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
        repackStencil(flat_stencil, stencil, Ip, stencil_size, include_moments);
        return;
    }

    // Mirror symmetry around diagonal x = y.
    if (cy > cx /*+ 1e-8*/) {
        // count small differences
        swap_xy(stencil, stencil_size, include_moments);
        if (include_moments >= 2) swap_xy(I);
    }
    // Pack back into flat storage
    repackStencil(flat_stencil, stencil, Ip, stencil_size, include_moments);
}

} // namespace IRL
