#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>

namespace nr {
namespace phy {

namespace {

constexpr int TYPE1_RE_POS[] = {0, 2, 4, 6, 8, 10};
constexpr int TYPE1_N_RE = 6;
constexpr int TYPE1_N_CDM_GROUPS = 2;

constexpr int TYPE2_RE_POS[] = {0, 1, 6, 7};
constexpr int TYPE2_N_RE = 4;
constexpr int TYPE2_N_CDM_GROUPS = 3;

constexpr double INV_SQRT2 = 0.7071067811865475;

Complex get_dmrs_symbol(const std::vector<uint8_t>& c_seq, int m) {
    int c0 = c_seq[2 * m];
    int c1 = c_seq[2 * m + 1];
    double re = (1 - 2 * c0) * INV_SQRT2;
    double im = (1 - 2 * c1) * INV_SQRT2;
    return Complex(re, im);
}

// TS 38.211 §7.4.1.1.1: DMRS sequence generation
// c_init = (2^17*(N_symb_slot*n_s,f + l + 1)*(2*N_ID + 1)
//           + 2^17*floor(lambda/2) + 2*N_ID + n_SCID) mod 2^31
// lambda = CDM group index
uint32_t compute_cinit(int n_symb_slot, int n_slot, int l, uint16_t N_id,
                       int cdm_group, int n_scid) {
    uint64_t t = (1ULL << 17) * (uint64_t)(n_symb_slot * n_slot + l + 1)
                 * ((uint64_t)(2 * N_id) + 1);
    t += ((1ULL << 17) * (uint64_t)(cdm_group / 2));
    t += (uint64_t)(2 * N_id) + (uint64_t)n_scid;
    return (uint32_t)(t % (1ULL << 31));
}

int get_dmrs_symbol_position(int additional_pos, int duration, int idx) {
    int dmrs_symbols[4];
    int n_dmrs = 1;
    dmrs_symbols[0] = 2;
    
    if (additional_pos >= 1) {
        dmrs_symbols[n_dmrs++] = (duration == 1) ? 11 : 7;
    }
    if (additional_pos >= 2) {
        dmrs_symbols[n_dmrs++] = (duration == 1) ? 7 : 11;
    }
    if (additional_pos >= 3) {
        dmrs_symbols[n_dmrs++] = 9;
    }
    
    if (idx >= 0 && idx < n_dmrs) {
        return dmrs_symbols[idx];
    }
    return -1;
}

int get_num_dmrs_symbols(int additional_pos) {
    int n_dmrs = 1;
    if (additional_pos >= 1) n_dmrs++;
    if (additional_pos >= 2) n_dmrs++;
    if (additional_pos >= 3) n_dmrs++;
    return n_dmrs;
}

}

class DmrsGeneratorImpl : public IDmrsGenerator {
public:
    void generate_dmrs(const SimulationConfig& config, ResourceGrid& grid,
                       int slot_idx, int symbol_start) override {
        int n_rb = config.n_rb;
        int n_sc = grid.n_subcarriers;
        int n_sym = grid.n_symbols;
        int n_ant = grid.n_ant;
        bool is_type1 = (config.dmrs_type == DmrsType::TYPE1);
        // TS 38.211 Table 7.4.1.1.2-1/-2: number of CDM groups and max ports
        // Type1: 2 CDM groups (delta=0/1), ports 1000-1003 (4 layers)
        // Type2: 3 CDM groups (delta=0/1/2), ports 1000-1005 (6 layers)
        int n_cdm_groups = is_type1 ? 2 : 3;
        int max_ports = is_type1 ? 4 : 6;
        int n_layers = std::min(config.n_layers, max_ports);

        int n_dmrs_sym = get_num_dmrs_symbols(config.dmrs_additional_pos);

        uint16_t N_id = static_cast<uint16_t>(config.n_id_dmrs);
        int n_symb_slot = 14;

        for (int dmrs_idx = 0; dmrs_idx < n_dmrs_sym; dmrs_idx++) {
            int l = get_dmrs_symbol_position(config.dmrs_additional_pos, config.dmrs_duration, dmrs_idx);
            if (l < 0 || l >= n_sym) continue;
            if (symbol_start >= 0 && l < symbol_start) continue;

            int sym = l;

            // One PN sequence per CDM group per symbol
            // seq covers r(2n+k') with n up to n_rb * (is_type1 ? 3 : 2), k' in {0,1}
            int seq_len = is_type1 ? (n_rb * 3 + 1) * 2 : (n_rb * 2 + 1) * 2;
            if (seq_len < 2) seq_len = 2;

            for (int layer = 0; layer < n_layers; layer++) {
                // Port 1000+p: CDM group = p % n_cdm_groups, freq OCC index = p / n_cdm_groups
                int cdm = layer % n_cdm_groups;
                int occ_idx = layer / n_cdm_groups;
                int delta = cdm;

                int ant = layer % n_ant;
                uint32_t cinit = compute_cinit(n_symb_slot, slot_idx, sym, N_id, cdm, 0);
                std::vector<uint8_t> c_seq;
                generate_pn_sequence(cinit, 2 * seq_len, c_seq);

                // Mapping per 38.211 7.4.1.1.2:
                //   Type1: k = 4n + 2k' + delta, k' = 0,1  (n = 0,1,2 within each PRB)
                //   Type2: k = 6n + k' + delta,  k' = 0,1  (n = 0,1 within each PRB)
                for (int prb = 0; prb < n_rb; prb++) {
                    int n_max = is_type1 ? 3 : 2;
                    for (int n = 0; n < n_max; n++) {
                        for (int kp = 0; kp < 2; kp++) {
                            int sc = prb * 12 + delta + (is_type1 ? 4 * n + 2 * kp : 6 * n + kp);
                            if (sc >= n_sc) continue;
                            // r(2n + k') : sequence index, continuous across PRBs
                            int r_idx = 2 * n_max * prb + 2 * n + kp;
                            if (2 * r_idx + 1 >= 2 * seq_len) continue;
                            // TS 38.211 7.4.1.1.2 freq OCC w_f(k') = [1, +/-1] applied
                            // across the RE PAIR (kp=0, kp=1): the second port of a
                            // CDM group keeps the same sign on the pair's first RE
                            // and flips only on the second RE. This makes the pair
                            // (y0 = (h_a+h_b)*r0, y1 = (h_a-h_b)*r1) separable at
                            // the receiver. A constant -1 over the whole comb would
                            // collapse every pair to a single inseparable sum.
                            double w_f = (occ_idx == 0 || kp == 0) ? 1.0 : -1.0;
                            Complex dmrs_val = get_dmrs_symbol(c_seq, r_idx) * w_f;
                            grid.set_re(ant, sym, sc, dmrs_val);
                        }
                    }
                }
            }
        }
    }
    
    ComplexCube extract_dmrs(const ResourceGrid& rx_grid, const SimulationConfig& config,
                             int /*slot_idx*/, int symbol_start) override {
        int n_rb = config.n_rb;
        int n_sc = rx_grid.n_subcarriers;
        int n_sym = rx_grid.n_symbols;
        int n_rx_ant = rx_grid.n_ant;
        int n_re_per_prb = (config.dmrs_type == DmrsType::TYPE1) ? TYPE1_N_RE : TYPE2_N_RE;
        const int* re_pos = (config.dmrs_type == DmrsType::TYPE1) ? TYPE1_RE_POS : TYPE2_RE_POS;
        
        int n_dmrs_sym = get_num_dmrs_symbols(config.dmrs_additional_pos);
        
        ComplexCube dmrs_grid(n_sc, n_sym, n_rx_ant, arma::fill::zeros);
        
        for (int dmrs_idx = 0; dmrs_idx < n_dmrs_sym; dmrs_idx++) {
            int l = get_dmrs_symbol_position(config.dmrs_additional_pos, config.dmrs_duration, dmrs_idx);
            if (l < 0 || l >= n_sym) continue;
            if (symbol_start >= 0 && l < symbol_start) continue;
            
            int sym = l;
            
            for (int prb = 0; prb < n_rb; prb++) {
                for (int re_idx = 0; re_idx < n_re_per_prb; re_idx++) {
                    int sc = prb * 12 + re_pos[re_idx];
                    if (sc >= n_sc) continue;
                    
                    for (int rx_ant = 0; rx_ant < n_rx_ant; rx_ant++) {
                        Complex val = rx_grid.get_re(rx_ant, sym, sc);
                        dmrs_grid(sc, sym, rx_ant) = val;
                    }
                }
            }
        }
        
        return dmrs_grid;
    }
};

std::unique_ptr<IDmrsGenerator> create_dmrs_generator() {
    return std::make_unique<DmrsGeneratorImpl>();
}

}
}
