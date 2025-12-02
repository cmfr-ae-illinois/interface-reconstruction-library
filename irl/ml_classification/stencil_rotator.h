#pragma once
#include <vector>
#include <cmath>

namespace IRL {

// Per-cell data: volume fraction and first moments
struct CellData {
    double vfrac;
    double mx, my, mz;  // first moments in x, y, z
};

// Convert (i,j,k) to flat cell index in the 1D array of cells
inline int cellIndex(int i, int j, int k, int N) {
    return (i * N * N + j * N + k);
}

// Read one stencil from a flat [vfrac, mx, my, mz] array into CellData array
static void unpackStencil(const std::vector<double>& flat,
                          std::vector<CellData>& stencil, int N)
{
    const int nCells = N * N * N;
    stencil.resize(nCells);
    for (int idx = 0; idx < nCells; ++idx) {
        stencil[idx].vfrac = flat[4 * idx + 0];
        stencil[idx].mx    = flat[4 * idx + 1];
        stencil[idx].my    = flat[4 * idx + 2];
        stencil[idx].mz    = flat[4 * idx + 3];
    }
}

// Write CellData stencil back into flat [vfrac, mx, my, mz] array
static void repackStencil(std::vector<double>& flat,
                          const std::vector<CellData>& stencil, int N)
{
    const int nCells = N * N * N;
    flat.resize(4 * nCells);
    for (int idx = 0; idx < nCells; ++idx) {
        flat[4 * idx + 0] = stencil[idx].vfrac;
        flat[4 * idx + 1] = stencil[idx].mx;
        flat[4 * idx + 2] = stencil[idx].my;
        flat[4 * idx + 3] = stencil[idx].mz;
    }
}

// Symmetry helpers


// Reflections into the first octant

// Reflect across plane x -> -x  (mirror in x, i <-> N-1-i; mx -> -mx)
inline void reflect_x(std::vector<CellData>& S, int N) {
    std::vector<CellData> C = S;  // copy source
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k) {
                CellData d = C[cellIndex(N - 1 - i, j, k, N)];
                // First moment transforms as a vector under reflection
                d.mx = -d.mx;       // x component flips sign
                // y,z unchanged
                S[cellIndex(i, j, k, N)] = d;
            }
}

// Reflect across plane y -> -y  (mirror in y, j <-> N-1-j; my -> -my)
inline void reflect_y(std::vector<CellData>& S, int N) {
    std::vector<CellData> C = S;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k) {
                CellData d = C[cellIndex(i, N - 1 - j, k, N)];
                d.my = -d.my;       // y component flips sign
                S[cellIndex(i, j, k, N)] = d;
            }
}

// Reflect across plane z -> -z  (mirror in z, k <-> N-1-k; mz -> -mz)
inline void reflect_z(std::vector<CellData>& S, int N) {
    std::vector<CellData> C = S;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k) {
                CellData d = C[cellIndex(i, j, N - 1 - k, N)];
                d.mz = -d.mz;       // z component flips sign
                S[cellIndex(i, j, k, N)] = d;
            }
}

// 120° rotations about the (1,1,1) diagonal

inline void permute_xyz(std::vector<CellData>& S, int N, int type) {
    std::vector<CellData> C = S;
    if (type == 1) {
        // (x,y,z)->(y,z,x)
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                for (int k = 0; k < N; ++k) {
                    CellData d = C[cellIndex(j, k, i, N)];
                    double mx_old = d.mx, my_old = d.my, mz_old = d.mz;
                    d.mx = my_old;  // x' component = old y
                    d.my = mz_old;  // y' component = old z
                    d.mz = mx_old;  // z' component = old x
                    S[cellIndex(i, j, k, N)] = d;
                }
    } else if (type == 2) {
        // (x,y,z)->(z,x,y)
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                for (int k = 0; k < N; ++k) {
                    CellData d = C[cellIndex(k, i, j, N)];
                    double mx_old = d.mx, my_old = d.my, mz_old = d.mz;
                    d.mx = mz_old;  // x' component = old z
                    d.my = mx_old;  // y' component = old x
                    d.mz = my_old;  // z' component = old y
                    S[cellIndex(i, j, k, N)] = d;
                }
    }
}

// Reflection across the diagonal plane x = y

inline void swap_xy(std::vector<CellData>& S, int N) {
    std::vector<CellData> C = S;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k) {
                CellData d = C[cellIndex(j, i, k, N)];
                double mx_old = d.mx, my_old = d.my;
                d.mx = my_old;  // x' component = old y
                d.my = mx_old;  // y' component = old x
                // z unchanged
                S[cellIndex(i, j, k, N)] = d;
            }
}


inline void rotate_stencil(std::vector<double>& flat_stencil,
                           int stencil_size, int no_symmetries)
{
    if (no_symmetries < 8) {
        return; // no rotation
    }

    // Unpack the flat vector into structured cells
    std::vector<CellData> stencil;
    unpackStencil(flat_stencil, stencil, stencil_size);

    // Compute global centroid of first moments
    double cx = 0.0, cy = 0.0, cz = 0.0;
    const int stencil_volume = stencil_size * stencil_size * stencil_size;
    for (const auto& c : stencil) {
        cx += c.mx;
        cy += c.my;
        cz += c.mz;
    }
    cx /= stencil_volume;
    cy /= stencil_volume;
    cz /= stencil_volume;

    // Reflect into first octant: ensure cx,cy,cz >= 0 using x/y/z reflections.
    if (cx < 0.0) {
        reflect_x(stencil, stencil_size);
        cx = -cx;   // centroid x-component flips sign as well
    }
    if (cy < 0.0) {
        reflect_y(stencil, stencil_size);
        cy = -cy;
    }
    if (cz < 0.0) {
        reflect_z(stencil, stencil_size);
        cz = -cz;
    }

    if (no_symmetries <= 8) {
        repackStencil(flat_stencil, stencil, stencil_size);
        return;
    }

    // 120° rotations about (1,1,1) diagonal.

    double z0 = cz;
    double z1 = cx;
    double z2 = cy;

    int best = 0;
    double bestz = z0;
    if (z1 > bestz) { best = 1; bestz = z1; }
    if (z2 > bestz) { best = 2; bestz = z2; }

    if (best == 1) {
        permute_xyz(stencil, stencil_size, 1);
        // (cx,cy,cz) -> (cy,cz,cx)
        double tmp = cx;
        cx = cy;
        cy = cz;
        cz = tmp;
    } else if (best == 2) {
        permute_xyz(stencil, stencil_size, 2);
        // (cx,cy,cz) -> (cz,cx,cy)
        double tmp = cx;
        cx = cz;
        cz = cy;
        cy = tmp;
    }

    // If we only want 24 symmetries (rotations + octant choice), stop here.
    if (no_symmetries <= 24) {
        repackStencil(flat_stencil, stencil, stencil_size);
        return;
    }

    // Mirror symmetry around diagonal x = y.

    swap_xy(stencil, stencil_size);

    // Pack back into flat storage
    repackStencil(flat_stencil, stencil, stencil_size);
}

} // namespace IRL
