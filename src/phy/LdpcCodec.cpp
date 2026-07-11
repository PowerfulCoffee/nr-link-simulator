#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <limits>

namespace nr {
namespace phy {

namespace {

constexpr int BG1_K_B = 22;
constexpr int BG1_CORE_COLS = BG1_K_B + 4;
constexpr int BG1_N_COLS_TOTAL = 68;

constexpr int BG2_K_B = 10;
constexpr int BG2_CORE_COLS = BG2_K_B + 4;
constexpr int BG2_N_COLS_TOTAL = 52;

constexpr int BG1_CORE_ROWS = 4;

constexpr int BG1_CORE_SHIFT[BG1_CORE_ROWS][BG1_CORE_COLS] = {
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}
};

constexpr int BG2_CORE_SHIFT[BG1_CORE_ROWS][BG2_CORE_COLS] = {
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}
};

inline int mod_positive(int a, int zc) {
    int r = a % zc;
    return r < 0 ? r + zc : r;
}

void circ_xor_accum(const uint8_t* src, uint8_t* dst, int shift, int zc) {
    for (int i = 0; i < zc; i++) {
        dst[i] ^= src[mod_positive(i - shift, zc)];
    }
}

}

class LdpcEncoderImpl : public ILdpcEncoder {
public:
    BitVec encode(const BitVec& info_bits, int bgn, int zc) override {
        int k_b = (bgn == 1) ? BG1_K_B : BG2_K_B;
        int n_core_checks = 4;
        int n_total_cols = (bgn == 1) ? BG1_N_COLS_TOTAL : BG2_N_COLS_TOTAL;
        int n_coded = (n_total_cols - 2) * zc;
        int k_total = k_b * zc;
        
        int n_checks_max = std::min(15, n_total_cols - k_b);
        
        BitVec codeword(k_total + n_checks_max * zc, arma::fill::zeros);
        
        int info_len = std::min((int)info_bits.n_elem, k_total);
        for (int i = 0; i < info_len; i++) {
            codeword(i) = info_bits(i);
        }
        
        uint8_t* p0 = codeword.memptr() + k_b * zc;
        uint8_t* p1 = p0 + zc;
        uint8_t* p2 = p1 + zc;
        uint8_t* p3 = p2 + zc;
        
        BitVec temp(zc, arma::fill::zeros);
        
        for (int col = 0; col < k_b; col++) {
            const uint8_t* col_ptr = codeword.memptr() + col * zc;
            circ_xor_accum(col_ptr, temp.memptr(), 0, zc);
        }
        for (int i = 0; i < zc; i++) {
            p0[i] = temp(i);
        }
        
        temp.zeros();
        circ_xor_accum(p0, temp.memptr(), 0, zc);
        for (int col = 0; col < k_b; col++) {
            const uint8_t* col_ptr = codeword.memptr() + col * zc;
            circ_xor_accum(col_ptr, temp.memptr(), 0, zc);
        }
        for (int i = 0; i < zc; i++) {
            p1[i] = temp(i);
        }
        
        temp.zeros();
        circ_xor_accum(p1, temp.memptr(), 0, zc);
        for (int col = 0; col < k_b; col++) {
            const uint8_t* col_ptr = codeword.memptr() + col * zc;
            circ_xor_accum(col_ptr, temp.memptr(), 0, zc);
        }
        for (int i = 0; i < zc; i++) {
            p2[i] = temp(i);
        }
        
        temp.zeros();
        circ_xor_accum(p2, temp.memptr(), 0, zc);
        for (int col = 0; col < k_b; col++) {
            const uint8_t* col_ptr = codeword.memptr() + col * zc;
            circ_xor_accum(col_ptr, temp.memptr(), 0, zc);
        }
        for (int i = 0; i < zc; i++) {
            p3[i] = temp(i);
        }
        
        for (int row_ext = 0; row_ext < n_checks_max - 4; row_ext++) {
            uint8_t* p_ext = codeword.memptr() + k_b * zc + (4 + row_ext) * zc;
            temp.zeros();
            
            uint8_t* prev_p = p3;
            if (row_ext > 0) {
                prev_p = codeword.memptr() + k_b * zc + (3 + row_ext) * zc;
            }
            for (int i = 0; i < zc; i++) {
                temp(i) ^= prev_p[i];
            }
            
            for (int col = 0; col < k_b; col++) {
                const uint8_t* col_ptr = codeword.memptr() + col * zc;
                circ_xor_accum(col_ptr, temp.memptr(), 0, zc);
            }
            
            for (int i = 0; i < zc; i++) {
                p_ext[i] = temp(i);
            }
        }
        
        int out_len = std::min(n_coded, (int)codeword.n_elem);
        BitVec output(out_len, arma::fill::zeros);
        for (int i = 0; i < out_len; i++) {
            output(i) = codeword(i);
        }
        
        return output;
    }
};

class LdpcDecoderImpl : public ILdpcDecoder {
public:
    std::pair<BitVec, bool> decode(const SoftVec& llr_in, int bgn, int zc,
                                   int n_iter, bool early_term) override {
        int k_b = (bgn == 1) ? BG1_K_B : BG2_K_B;
        int n_total_cols = (bgn == 1) ? BG1_N_COLS_TOTAL : BG2_N_COLS_TOTAL;
        int n_rows = std::min(15, n_total_cols - k_b);
        int k_total = k_b * zc;
        
        SoftVec app(n_total_cols * zc, arma::fill::zeros);
        
        int copy_len = std::min((int)llr_in.n_elem, k_total + n_rows * zc);
        for (int i = 0; i < copy_len; i++) {
            app(i) = llr_in(i);
        }
        
        SoftVec msg_c2v(n_rows * zc, arma::fill::zeros);
        
        double norm_factor = 0.75;
        bool converged = false;
        
        for (int iter = 0; iter < n_iter; iter++) {
            for (int row = 0; row < n_rows; row++) {
                std::vector<std::pair<int, int>> edges;
                
                for (int col = 0; col < k_b; col++) {
                    edges.emplace_back(col, 0);
                }
                edges.emplace_back(k_b + row, 0);
                if (row > 0) {
                    edges.emplace_back(k_b + row - 1, 0);
                }
                
                int n_edges = edges.size();
                std::vector<SoftVec> var_msgs(n_edges, SoftVec(zc));
                
                for (int e = 0; e < n_edges; e++) {
                    int col = edges[e].first;
                    int shift = edges[e].second;
                    const double* app_ptr = app.memptr() + col * zc;
                    const double* old_c2v = msg_c2v.memptr() + row * zc;
                    
                    for (int i = 0; i < zc; i++) {
                        int src_i = mod_positive(i - shift, zc);
                        var_msgs[e](i) = app_ptr[src_i] - old_c2v[src_i];
                    }
                }
                
                SoftVec new_c2v(zc, arma::fill::zeros);
                for (int i = 0; i < zc; i++) {
                    double min1 = std::numeric_limits<double>::infinity();
                    double min2 = std::numeric_limits<double>::infinity();
                    int sign_prod = 1;
                    int min_idx = -1;
                    
                    for (int e = 0; e < n_edges; e++) {
                        double v = var_msgs[e](i);
                        double mag = std::abs(v);
                        int sgn = (v >= 0) ? 1 : -1;
                        sign_prod *= sgn;
                        
                        if (mag < min1) {
                            min2 = min1;
                            min1 = mag;
                            min_idx = e;
                        } else if (mag < min2) {
                            min2 = mag;
                        }
                    }
                    
                    for (int e = 0; e < n_edges; e++) {
                        double out_mag = (e == min_idx) ? min2 : min1;
                        double v = var_msgs[e](i);
                        int sgn = (v >= 0) ? 1 : -1;
                        int out_sgn = sign_prod * sgn;
                        var_msgs[e](i) = out_sgn * out_mag * norm_factor;
                    }
                }
                
                for (int e = 0; e < n_edges; e++) {
                    int col = edges[e].first;
                    int shift = edges[e].second;
                    double* app_ptr = app.memptr() + col * zc;
                    double* c2v_ptr = msg_c2v.memptr() + row * zc;
                    
                    for (int i = 0; i < zc; i++) {
                        int dst_i = mod_positive(i - shift, zc);
                        c2v_ptr[dst_i] = var_msgs[e](i);
                        app_ptr[dst_i] += c2v_ptr[dst_i];
                    }
                }
            }
            
            if (early_term) {
                bool syndrome_ok = true;
                BitVec hard(k_total + 4 * zc);
                for (int i = 0; i < k_total + 4 * zc; i++) {
                    hard(i) = (app(i) < 0) ? 1 : 0;
                }
                
                for (int row = 0; row < 4; row++) {
                    for (int i = 0; i < zc; i++) {
                        uint8_t parity = 0;
                        for (int col = 0; col < k_b; col++) {
                            parity ^= hard(col * zc + i);
                        }
                        parity ^= hard(k_b * zc + i);
                        if (row > 0) {
                            parity ^= hard(k_b * zc + (row - 1) * zc + i);
                        }
                        parity ^= hard(k_b * zc + row * zc + i);
                        
                        if (parity != 0) {
                            syndrome_ok = false;
                            break;
                        }
                    }
                    if (!syndrome_ok) break;
                }
                
                if (syndrome_ok) {
                    converged = true;
                    break;
                }
            }
        }
        
        BitVec decoded(k_total, arma::fill::zeros);
        for (int i = 0; i < k_total; i++) {
            decoded(i) = (app(i) < 0) ? 1 : 0;
        }
        
        return {decoded, converged};
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
