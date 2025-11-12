#pragma once
#include <torch/torch.h>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>

namespace IRL {

// ===================== Helpers =====================
inline int64_t cube_index_3d_to_1d(int s, int i,int j,int k){
    // row-major: i fastest or slowest? We'll keep i fastest -> (i + s*j + s*s*k)
    // But we must be consistent everywhere. We'll pick: (i*s + j)*s + k (i slowest).
    return (static_cast<int64_t>(i)*s + j)*s + k;
}

// ===================== Cube group (48 rotations incl. mirrors) =====================
// We generate all 3x3 signed permutation matrices (entries in {0,±1}); that's 3!*2^3 = 48.
// Each acts on grid coordinates (centered) exactly; we precompute the induced node permutation.
struct CubeElement {
    torch::Tensor R;      // [3,3] float32 (0,±1)
    torch::Tensor P;      // [N,N] float32 permutation matrix (optional for projection)
    torch::Tensor perm;   // [N]   int64 permutation indices for fast index_select
};

struct CubeGroup {
    int stencil_size;
    int64_t N;        // N = s^3
    int c;            // center index (s/2)
    std::vector<CubeElement> G;

    CubeGroup(int stencil_size_, torch::Device device, c10::ScalarType dtype)
        : stencil_size(stencil_size_), N(static_cast<int64_t>(stencil_size_)*stencil_size_*stencil_size_), c(stencil_size_/2)
    {
        if (stencil_size < 1) throw std::invalid_argument("stencil_size must be >=1");
        // Generate all permutations of axes
        std::array<int,3> base = {0,1,2};
        std::vector<std::array<int,3>> perms;
        std::sort(base.begin(), base.end());
        do { perms.push_back(base); } while(std::next_permutation(base.begin(), base.end()));
        // All sign combinations (2^3)
        const std::array<std::array<int,3>,8> signs = {{
            {+1,+1,+1},{+1,+1,-1},{+1,-1,+1},{+1,-1,-1},
            {-1,+1,+1},{-1,+1,-1},{-1,-1,+1},{-1,-1,-1}
        }};

        for (const auto& p : perms){
            for (const auto& s : signs){
                // Build R whose columns are signed basis vectors e_{p[c]} * s[c]
                auto R = torch::zeros({3,3}, torch::TensorOptions().dtype(dtype).device(device));
                for (int col=0; col<3; ++col){
                    R.index_put_({p[col], col}, static_cast<float>(s[col]));
                }
                auto perm_idx = make_perm_index(R).to(torch::kLong).to(device);
                auto P = make_perm_matrix(perm_idx, device, dtype);
                // Avoid duplicates (shouldn't occur)
                bool dup=false;
                for (auto& g : G){
                    if (torch::allclose(R, g.R)) { dup=true; break; }
                }
                if (!dup) G.push_back({R, P, perm_idx});
            }
        }
        // We expect 48 elements
        // fprintf(stderr, "O_h size = %zu (expected 48)\n", G.size());
    }

    torch::Tensor make_perm_matrix(const torch::Tensor& perm_idx, torch::Device device, c10::ScalarType dtype){
        auto P = torch::zeros({N,N}, torch::TensorOptions().dtype(dtype).device(device));
        auto pacc = perm_idx.accessor<int64_t,1>();
        for (int64_t i=0;i<N;++i){
            P.index_put_({i, pacc[i]}, 1.0f);
        }
        return P;
    }

    torch::Tensor make_perm_index(const torch::Tensor& R){
        // perm[n] = n' such that node n maps to n' under R
        std::vector<int64_t> v(N);
        for (int i=0;i<stencil_size;++i)
        for (int j=0;j<stencil_size;++j)
        for (int k=0;k<stencil_size;++k){
            // position in centered integer coordinates
            auto r = torch::tensor({float(i-c), float(j-c), float(k-c)}, R.options());
            auto rp = torch::matmul(R, r);
            int ip = int(std::llround(rp[0].item<double>())) + c;
            int jp = int(std::llround(rp[1].item<double>())) + c;
            int kp = int(std::llround(rp[2].item<double>())) + c;
            ip = std::clamp(ip, 0, stencil_size-1);
            jp = std::clamp(jp, 0, stencil_size-1);
            kp = std::clamp(kp, 0, stencil_size-1);
            int64_t n  = cube_index_3d_to_1d(stencil_size, i,j,k);
            int64_t np = cube_index_3d_to_1d(stencil_size, ip,jp,kp);
            v[n] = np;
        }
        return torch::from_blob(v.data(), {N}, torch::TensorOptions().dtype(torch::kLong)).clone();
    }
};

// ===================== Equivariant dense layer (no conv) =====================
//
// Channels are split into scalars (l=0) and vectors (l=1 with 3 components).
// Node mixing matrices are projected into the commutant of the group: Π(A) = (1/|G|) Σ_g P_g A P_g^T
//
struct EquivariantDenseImpl : torch::nn::Module {
    int64_t N;                   // number of nodes (stencil_size^3)
    int S_in, V_in, S_out, V_out;
    torch::Tensor W_ss;          // [S_in, S_out] scalar channel mixer
    torch::Tensor W_vv;          // [V_in, V_out] vector channel mixer (does not touch the 3 spatial comps)
    torch::Tensor A_s;           // [N, N] base node mixing (scalars)
    torch::Tensor A_v;           // [N, N] base node mixing (vectors)
    CubeGroup group;             // carries P_g and R_g

    EquivariantDenseImpl(int stencil_size, int S_in_, int V_in_, int S_out_, int V_out_,
                         torch::Device device, c10::ScalarType dtype)
        : N(static_cast<int64_t>(stencil_size)*stencil_size*stencil_size),
          S_in(S_in_), V_in(V_in_), S_out(S_out_), V_out(V_out_),
          group(stencil_size, device, dtype)
    {
        W_ss = register_parameter("W_ss", torch::randn({S_in, S_out}, torch::TensorOptions().dtype(dtype).device(device))*0.1);
        W_vv = register_parameter("W_vv", torch::randn({V_in, V_out}, torch::TensorOptions().dtype(dtype).device(device))*0.1);
        A_s  = register_parameter("A_s",  torch::randn({N, N}, torch::TensorOptions().dtype(dtype).device(device))*0.01);
        A_v  = register_parameter("A_v",  torch::randn({N, N}, torch::TensorOptions().dtype(dtype).device(device))*0.01);
    }

    torch::Tensor project_commutant(const torch::Tensor& A) {
        // Π(A) = (1/|G|) Σ_g P_g A P_g^T
        auto out = torch::zeros_like(A);
        for (const auto& g : group.G) {
            out.addmm_(g.P, torch::mm(A, g.P.transpose(0,1)));
        }
        out /= static_cast<float>(group.G.size());
        return out;
    }

    std::pair<torch::Tensor, torch::Tensor>
    forward(const torch::Tensor& s_in, const torch::Tensor& v_in) {
        // s_in: [B,N,S_in]
        // v_in: [B,N,V_in,3]
        auto s_ch = torch::matmul(s_in, W_ss);                        // [B,N,S_out]
        auto v_ch = torch::einsum("bnfc,fg->bngc", {v_in, W_vv});     // [B,N,V_out,3]

        // Node mixing (projected)
        auto As = project_commutant(A_s);                              // [N,N]
        auto Av = project_commutant(A_v);                              // [N,N]
        auto s_out = torch::einsum("bns,nm->bms", {s_ch, As});         // [B,N,S_out]
        auto v_out = torch::einsum("bnvc,nm->bmvc", {v_ch, Av});       // [B,N,V_out,3]
        return {s_out, v_out};
    }
};
TORCH_MODULE(EquivariantDense);

// ===================== Gated nonlinearity =====================
//
// Scalars: GELU. Vectors: scale magnitudes by sigmoid(gate(s)), keep directions.
//
struct GatedNonlinImpl : torch::nn::Module {
    torch::nn::Linear gate{nullptr}; // S_out → V_out
    GatedNonlinImpl(int S_out, int V_out, torch::Device device, c10::ScalarType dtype) {
        gate = register_module("gate", torch::nn::Linear(torch::nn::LinearOptions(S_out, V_out).bias(true)));
        gate->to(device, dtype);
    }

    std::pair<torch::Tensor, torch::Tensor>
    forward(const torch::Tensor& s, const torch::Tensor& v) {
        // s: [B,N,S], v: [B,N,V,3]
        auto s_act = torch::gelu(s);
        auto g = torch::sigmoid(gate->forward(s_act)); // [B,N,V]
        g = g.unsqueeze(-1);                           // [B,N,V,1]
        auto v_norm = (v.pow(2).sum(-1, true) + 1e-8).sqrt(); // [B,N,V,1]
        auto v_dir  = v / v_norm;
        auto v_act  = g * v_dir * v_norm; // scale magnitudes only
        return {s_act, v_act};
    }
};
TORCH_MODULE(GatedNonlin);

// ===================== e3nn model (non-convolutional) =====================
//
// Public name: e3nn   (replaces your previous Net)
// Two constructors:
//   (A) e3nn(int stencil_size, int hidden1, int hidden2, int hidden3, int output)
//   (B) e3nn(int input_size,  int hidden1, int hidden2, int hidden3, int output)  // compat
//
// Input format: [B, 4*N] flattened per-cell (vof, mx, my, mz)
// Output: logits [B, output_size]
//
struct e3nnImpl : torch::nn::Module {
    int stencil_size;
    int64_t N;
    EquivariantDense l1{nullptr}, l2{nullptr};
    GatedNonlin nl1{nullptr}, nl2{nullptr};
    torch::nn::Linear head{nullptr};

    torch::Device device;
    c10::ScalarType dtype;

    // Preferred constructor
    e3nnImpl(int stencil_size_,
             int hidden_size1, int hidden_size2, int /*hidden_size3*/, int output_size,
             torch::Device dev = torch::kCPU, c10::ScalarType dt = torch::kFloat32)
        : stencil_size(stencil_size_),
          N(static_cast<int64_t>(stencil_size_)*stencil_size_*stencil_size_),
          device(dev), dtype(dt)
    {
        const int S_in = 1, V_in = 1;
        const int S_h1 = hidden_size1, V_h1 = hidden_size1;
        const int S_h2 = hidden_size2, V_h2 = hidden_size2;

        l1  = register_module("l1",  EquivariantDense(stencil_size, S_in, V_in, S_h1, V_h1, device, dtype));
        nl1 = register_module("nl1", GatedNonlin(S_h1, V_h1, device, dtype));
        l2  = register_module("l2",  EquivariantDense(stencil_size, S_h1, V_h1, S_h2, V_h2, device, dtype));
        nl2 = register_module("nl2", GatedNonlin(S_h2, V_h2, device, dtype));
        head = register_module("head", torch::nn::Linear(torch::nn::LinearOptions(S_h2 + V_h2, output_size)));
        head->to(device, dtype);
    }

    // Tag type for backward-compat constructor
    struct from_input_size_t { explicit from_input_size_t() = default; };
    static constexpr from_input_size_t from_input_size{};

    // Backward-compat constructor: derive stencil_size from input_size = 4*N
    e3nnImpl(from_input_size_t, int input_size,
            int hidden_size1, int hidden_size2, int hidden_size3, int output_size)
        : e3nnImpl(derive_stencil_size_from_input(input_size), hidden_size1, hidden_size2, hidden_size3, output_size)
    {}

    static int derive_stencil_size_from_input(int input_size){
        if (input_size % 4 != 0) throw std::invalid_argument("input_size must be divisible by 4 (vof+3 moments per cell).");
        double Nf = static_cast<double>(input_size) / 4.0;
        double s = std::round(std::cbrt(Nf));
        if (std::abs(s*s*s - Nf) > 1e-6) {
            throw std::invalid_argument("input_size/4 must be a perfect cube.");
        }
        int si = static_cast<int>(s);
        if (si <= 0) throw std::invalid_argument("Derived stencil_size must be > 0.");
        return si;
    }

    torch::Tensor forward(torch::Tensor x) {
        // x: [B, 4*N] flattened (vof, mx, my, mz) for N = stencil_size^3 cells
        x = x.to(device, dtype);
        const int64_t B = x.size(0);
        const int64_t expected = 4 * N;
        if (x.size(1) != expected) {
            throw std::invalid_argument("e3nn::forward - input second dimension must be 4*N where N = stencil_size^3.");
        }

        // Unflatten into scalar + vector node features
        auto s_in = torch::zeros({B, N, 1}, x.options());     // scalars
        auto v_in = torch::zeros({B, N, 1, 3}, x.options());  // vectors
        for (int64_t n=0; n<N; ++n){
            s_in.index_put_({torch::indexing::Slice(), n, 0},   x.index({torch::indexing::Slice(), 4*n + 0}));
            v_in.index_put_({torch::indexing::Slice(), n, 0, 0}, x.index({torch::indexing::Slice(), 4*n + 1}));
            v_in.index_put_({torch::indexing::Slice(), n, 0, 1}, x.index({torch::indexing::Slice(), 4*n + 2}));
            v_in.index_put_({torch::indexing::Slice(), n, 0, 2}, x.index({torch::indexing::Slice(), 4*n + 3}));
        }

        // Two equivariant blocks
        auto [s1, v1] = l1->forward(s_in, v_in);
        std::tie(s1, v1) = nl1->forward(s1, v1);
        auto [s2, v2] = l2->forward(s1, v1);
        std::tie(s2, v2) = nl2->forward(s2, v2);

        // Invariant pooling to logits
        auto s_pool = s2.mean(1);                     // [B, S_h2]
        auto v_pool = v2.norm(2, /*dim=*/-1).mean(1); // [B, V_h2]
        auto pooled = torch::cat({s_pool, v_pool}, -1);
        auto logits = head->forward(pooled);          // [B, output_size]
        return logits;
    }
};
TORCH_MODULE(e3nn);

} // namespace IRL
