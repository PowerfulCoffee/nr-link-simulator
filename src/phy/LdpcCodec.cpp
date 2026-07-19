#include "phy/PhyInterfaces.h"
#include "phy/LdpcTables.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <utility>
#include <bitset>
#include <unordered_map>
#include <memory>

namespace nr {
namespace phy {

using namespace nr::ldpc;

namespace {

int find_zc_set(int zc) {
    for (int s = 0; s < 8; s++) {
        int len = ZC_SET_LENGTHS[s];
        for (int i = 0; i < len; i++) {
            if (ZC_SETS[s][i] == zc) return s;
        }
    }
    return 0;
}

int get_bg_shift(int bgn, int set_idx, int row, int col) {
    if (bgn == 1) {
        return BG1_SHIFTS[set_idx][row][col];
    } else {
        return BG2_SHIFTS[set_idx][row][col];
    }
}

struct LdpcDims {
    int rows;
    int cols;
    int n_sys;
    int N;
    int K;
};

LdpcDims get_dims(int bgn, int zc) {
    LdpcDims d;
    if (bgn == 1) {
        d.rows = 46;
        d.cols = 68;
    } else {
        d.rows = 42;
        d.cols = 52;
    }
    d.n_sys = d.cols - d.rows;
    d.N = d.cols * zc;
    d.K = d.n_sys * zc;
    return d;
}

struct Edge {
    int col;
    int shift;
};

std::vector<std::vector<Edge>> build_edges(int bgn, int zc) {
    auto d = get_dims(bgn, zc);
    int set_idx = find_zc_set(zc);
    std::vector<std::vector<Edge>> edges(d.rows);
    for (int i = 0; i < d.rows; i++) {
        for (int j = 0; j < d.cols; j++) {
            int sh = get_bg_shift(bgn, set_idx, i, j);
            if (sh >= 0) {
                int sh_mod = sh % zc;
                if (sh_mod < 0) sh_mod += zc;
                edges[i].push_back({j, sh_mod});
            }
        }
    }
    return edges;
}

struct CheckToVarEdge {
    int var_node;
};

struct VarToCheckEdge {
    int check_node;
    int edge_idx;
};

struct LdpcGraph {
    int M;
    int N;
    int K;
    int zc;
    int bgn;
    int total_edges;
    
    std::vector<int> check_edge_offset;
    std::vector<std::vector<CheckToVarEdge>> check_edges;
    std::vector<std::vector<VarToCheckEdge>> var_edges;
    std::vector<std::vector<int>> parity_check_edges;
};

std::shared_ptr<LdpcGraph> build_graph(int bgn, int zc) {
    static std::unordered_map<uint64_t, std::shared_ptr<LdpcGraph>> cache;
    uint64_t key = (static_cast<uint64_t>(bgn) << 32) | static_cast<uint32_t>(zc);
    auto it = cache.find(key);
    if (it != cache.end()) {
        return it->second;
    }
    
    auto d = get_dims(bgn, zc);
    auto bg_edges = build_edges(bgn, zc);
    
    auto g = std::make_shared<LdpcGraph>();
    g->M = d.rows * zc;
    g->N = d.N;
    g->K = d.K;
    g->zc = zc;
    g->bgn = bgn;
    g->check_edge_offset.resize(g->M + 1, 0);
    g->check_edges.resize(g->M);
    g->var_edges.resize(g->N);
    g->parity_check_edges.resize(g->M);
    
    for (int i = 0; i < d.rows; i++) {
        for (int l = 0; l < zc; l++) {
            int cnode = i * zc + l;
            for (size_t e = 0; e < bg_edges[i].size(); e++) {
                int j = bg_edges[i][e].col;
                int sh = bg_edges[i][e].shift;
                int vnode = j * zc + ((l + sh) % zc);
                g->check_edges[cnode].push_back({vnode});
                g->parity_check_edges[cnode].push_back(vnode);
            }
        }
    }
    
    g->total_edges = 0;
    for (int i = 0; i < g->M; i++) {
        g->check_edge_offset[i] = g->total_edges;
        int deg = (int)g->check_edges[i].size();
        for (int e = 0; e < deg; e++) {
            int vnode = g->check_edges[i][e].var_node;
            g->var_edges[vnode].push_back({i, g->total_edges + e});
        }
        g->total_edges += deg;
    }
    g->check_edge_offset[g->M] = g->total_edges;
    
    cache[key] = g;
    return g;
}

}

class LdpcEncoderImpl : public ILdpcEncoder {
public:
    BitVec encode(const BitVec& info_bits, int bgn, int zc) override {
        auto d = get_dims(bgn, zc);
        auto edges = build_edges(bgn, zc);

        if ((int)info_bits.size() > d.K) {
            throw std::runtime_error("info_bits too long for given zc/bgn");
        }

        BitVec cw(d.N, 0);
        int n_info = (int)info_bits.size();
        for (int i = 0; i < n_info; i++) {
            cw[i] = info_bits[i];
        }

        int Np = d.rows * zc;
        int Ncore = 4 * zc;
        std::vector<uint8_t> p(Np, 0);

        std::vector<uint8_t> phi(Np, 0);
        for (int i = 0; i < d.rows; i++) {
            for (const auto& e : edges[i]) {
                if (e.col < d.n_sys) {
                    for (int l = 0; l < zc; l++) {
                        phi[i * zc + l] ^= cw[e.col * zc + ((l + e.shift) % zc)];
                    }
                }
            }
        }

        {
            int nvar = Ncore;
            std::vector<std::vector<int>> mat(nvar);
            std::vector<uint8_t> rhs(nvar, 0);

            for (int i = 0; i < 4; i++) {
                for (int l = 0; l < zc; l++) {
                    int row = i * zc + l;
                    rhs[row] = phi[row];
                    for (const auto& e : edges[i]) {
                        if (e.col >= d.n_sys) {
                            int pt = e.col - d.n_sys;
                            if (pt < 4) {
                                int col = pt * zc + ((l + e.shift) % zc);
                                mat[row].push_back(col);
                            }
                        }
                    }
                }
            }

            std::vector<int> pivot_col(nvar, -1);
            int r = 0;
            for (int c = 0; c < nvar && r < nvar; c++) {
                int pr = -1;
                for (int rr = r; rr < nvar; rr++) {
                    for (int cc : mat[rr]) {
                        if (cc == c) { pr = rr; break; }
                    }
                    if (pr >= 0) break;
                }
                if (pr == -1) continue;
                std::swap(mat[r], mat[pr]);
                std::swap(rhs[r], rhs[pr]);
                pivot_col[r] = c;

                for (int rr = 0; rr < nvar; rr++) {
                    if (rr == r) continue;
                    bool has = false;
                    for (int cc : mat[rr]) {
                        if (cc == c) { has = true; break; }
                    }
                    if (has) {
                        int i1 = 0, i2 = 0;
                        std::vector<int> new_mat;
                        while (i1 < (int)mat[r].size() && i2 < (int)mat[rr].size()) {
                            if (mat[r][i1] < mat[rr][i2]) { new_mat.push_back(mat[r][i1++]); }
                            else if (mat[r][i1] > mat[rr][i2]) { new_mat.push_back(mat[rr][i2++]); }
                            else { i1++; i2++; }
                        }
                        while (i1 < (int)mat[r].size()) new_mat.push_back(mat[r][i1++]);
                        while (i2 < (int)mat[rr].size()) new_mat.push_back(mat[rr][i2++]);
                        mat[rr] = std::move(new_mat);
                        rhs[rr] ^= rhs[r];
                    }
                }
                r++;
            }

            for (int rr = 0; rr < nvar; rr++) {
                if (pivot_col[rr] >= 0) {
                    int c = pivot_col[rr];
                    p[c] = rhs[rr];
                }
            }
        }

        for (int i = 0; i < 4; i++) {
            for (int l = 0; l < zc; l++) {
                cw[d.K + i * zc + l] = p[i * zc + l];
            }
        }

        for (int i = 4; i < d.rows; i++) {
            for (int l = 0; l < zc; l++) {
                uint8_t val = phi[i * zc + l];
                for (const auto& e : edges[i]) {
                    if (e.col >= d.n_sys) {
                        int pt = e.col - d.n_sys;
                        if (pt != i) {
                            val ^= cw[d.K + pt * zc + ((l + e.shift) % zc)];
                        }
                    }
                }
                cw[d.K + i * zc + l] = val;
            }
        }

        return cw;
    }
};

namespace {

inline double phi(double x) {
    if (x < 1e-8) return 1e8;
    if (x > 40.0) return 0.0;
    return -std::log(std::tanh(x / 2.0));
}

}

class LdpcDecoderImpl : public ILdpcDecoder {
    mutable std::vector<double> msg_c2v_;
    mutable std::vector<double> msg_v2c_;
    mutable std::vector<double> ch_llr_;
    mutable std::vector<double> post_;
    
public:
    std::pair<BitVec, bool> decode(const SoftVec& llr_in, int bgn, int zc,
                                   int n_iter, bool early_term) override {
        auto graph = build_graph(bgn, zc);
        int M = graph->M;
        int N = graph->N;
        int K = graph->K;
        int E = graph->total_edges;
        
        const double llr_max = 20.0;
        
        ch_llr_.assign(N, 0.0);
        post_.assign(N, 0.0);
        msg_c2v_.assign(E, 0.0);
        msg_v2c_.assign(E, 0.0);
        
        int n_pcb = 2 * zc;
        int in_len = (int)llr_in.size();

        if (in_len == N - n_pcb) {
            for (int i = 0; i < in_len; i++) {
                double v = std::max(-llr_max, std::min(llr_max, llr_in[i]));
                ch_llr_[n_pcb + i] = v;
                post_[n_pcb + i] = v;
            }
        } else if (in_len >= N) {
            for (int i = 0; i < N; i++) {
                double v = std::max(-llr_max, std::min(llr_max, llr_in[i]));
                ch_llr_[i] = v;
                post_[i] = v;
            }
        } else {
            for (int i = 0; i < in_len; i++) {
                double v = std::max(-llr_max, std::min(llr_max, llr_in[i]));
                ch_llr_[i] = v;
                post_[i] = v;
            }
        }

        bool converged = false;
        std::vector<double> phi_buf;
        std::vector<int> sgn_buf;

        for (int iter = 0; iter < n_iter; iter++) {
            for (int v = 0; v < N; v++) {
                double post_v = post_[v];
                for (const auto& edge : graph->var_edges[v]) {
                    int eidx = edge.edge_idx;
                    msg_v2c_[eidx] = post_v - msg_c2v_[eidx];
                    msg_v2c_[eidx] = std::max(-llr_max, std::min(llr_max, msg_v2c_[eidx]));
                }
            }

            for (int i = 0; i < M; i++) {
                int deg = (int)graph->check_edges[i].size();
                if (deg == 0) continue;
                
                int base = graph->check_edge_offset[i];
                
                phi_buf.resize(deg);
                sgn_buf.resize(deg);
                int sgn_prod = 1;
                double phi_sum = 0.0;
                
                for (int e = 0; e < deg; e++) {
                    double val = msg_v2c_[base + e];
                    double a = std::fabs(val);
                    int s = (val >= 0) ? 1 : -1;
                    sgn_buf[e] = s;
                    sgn_prod *= s;
                    double p = phi(a);
                    phi_buf[e] = p;
                    phi_sum += p;
                }
                
                for (int e = 0; e < deg; e++) {
                    double a = std::fabs(msg_v2c_[base + e]);
                    double phi_sum_e = phi_sum - phi_buf[e];
                    double res_mag = phi(phi_sum_e);
                    if (!std::isfinite(res_mag) || res_mag < 0) res_mag = 0.0;
                    if (res_mag > llr_max) res_mag = llr_max;
                    int res_sgn = sgn_prod * sgn_buf[e];
                    msg_c2v_[base + e] = res_sgn * res_mag;
                }
            }

            for (int v = 0; v < N; v++) {
                double sum = ch_llr_[v];
                for (const auto& edge : graph->var_edges[v]) {
                    sum += msg_c2v_[edge.edge_idx];
                }
                post_[v] = std::max(-llr_max, std::min(llr_max, sum));
            }

            if (early_term) {
                bool ok = true;
                for (int i = 0; i < M; i++) {
                    uint8_t par = 0;
                    for (int v : graph->parity_check_edges[i]) {
                        par ^= (post_[v] < 0) ? 1 : 0;
                    }
                    if (par != 0) { ok = false; break; }
                }
                if (ok) { converged = true; break; }
            }
        }

        BitVec info_out(K);
        for (int i = 0; i < K; i++) {
            info_out[i] = (post_[i] < 0) ? 1 : 0;
        }

        return {info_out, converged};
    }
};

std::unique_ptr<ILdpcEncoder> create_ldpc_encoder() {
    return std::make_unique<LdpcEncoderImpl>();
}

std::unique_ptr<ILdpcDecoder> create_ldpc_decoder() {
    return std::make_unique<LdpcDecoderImpl>();
}

}
}
