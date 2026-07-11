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

uint32_t compute_cinit(int n_symb_slot, int n_slot, int l, uint16_t N_id) {
    uint64_t term1 = (1ULL << 17);
    uint64_t term2 = (uint64_t)(n_symb_slot * n_slot + l + 1);
    uint64_t term3 = 2ULL * N_id + 1ULL;
    uint32_t cinit = (uint32_t)(term1 * term2 * term3 + 2ULL * N_id);
    return cinit;
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
        int n_layers = std::min(config.n_layers, (config.dmrs_type == DmrsType::TYPE1) ? 4 : 6);
        
        int n_re_per_prb = (config.dmrs_type == DmrsType::TYPE1) ? TYPE1_N_RE : TYPE2_N_RE;
        int n_cdm_groups = (config.dmrs_type == DmrsType::TYPE1) ? TYPE1_N_CDM_GROUPS : TYPE2_N_CDM_GROUPS;
        const int* re_pos = (config.dmrs_type == DmrsType::TYPE1) ? TYPE1_RE_POS : TYPE2_RE_POS;
        
        int n_dmrs_sym = get_num_dmrs_symbols(config.dmrs_additional_pos);
        
        uint16_t N_id = 0;
        int n_symb_slot = 14;
        
        for (int dmrs_idx = 0; dmrs_idx < n_dmrs_sym; dmrs_idx++) {
            int l = get_dmrs_symbol_position(config.dmrs_additional_pos, config.dmrs_duration, dmrs_idx);
            if (l < 0 || l >= n_sym) continue;
            if (symbol_start >= 0 && l < symbol_start) continue;
            
            int sym = l;
            
            int seq_len = n_rb * n_re_per_prb;
            uint32_t cinit = compute_cinit(n_symb_slot, slot_idx, sym, N_id);
            std::vector<uint8_t> c_seq;
            generate_pn_sequence(cinit, 2 * seq_len, c_seq);
            
            for (int prb = 0; prb < n_rb; prb++) {
                for (int re_idx = 0; re_idx < n_re_per_prb; re_idx++) {
                    int sc = prb * 12 + re_pos[re_idx];
                    if (sc >= n_sc) continue;
                    
                    int m = prb * n_re_per_prb + re_idx;
                    if (m >= seq_len) continue;
                    
                    Complex base_sym = get_dmrs_symbol(c_seq, m);
                    
                    int cdm_group;
                    if (config.dmrs_type == DmrsType::TYPE1) {
                        cdm_group = re_idx / (n_re_per_prb / n_cdm_groups);
                        cdm_group = std::min(cdm_group, n_cdm_groups - 1);
                    } else {
                        if (re_idx < 2) cdm_group = 0;
                        else cdm_group = 1;
                    }
                    
                    for (int layer = 0; layer < n_layers; layer++) {
                        double w_f = 1.0;
                        
                        if (config.dmrs_type == DmrsType::TYPE1) {
                            int port_offset = layer;
                            if (port_offset < 2) {
                                w_f = (port_offset % 2 == 0) ? 1.0 : -1.0;
                                cdm_group = 0;
                            } else {
                                w_f = (port_offset % 2 == 0) ? 1.0 : -1.0;
                                cdm_group = 1;
                            }
                        } else {
                            int cdm_g = layer / 2;
                            int port_in_g = layer % 2;
                            cdm_group = std::min(cdm_g, n_cdm_groups - 1);
                            w_f = (port_in_g == 0) ? 1.0 : -1.0;
                        }
                        
                        Complex dmrs_val = base_sym * Complex(w_f, 0.0);
                        
                        if (config.n_layers <= n_ant) {
                            grid.set_re(layer % n_ant, sym, sc, dmrs_val);
                        } else {
                            double scale = 1.0 / std::sqrt((double)n_layers);
                            for (int ant = 0; ant < n_ant; ant++) {
                                grid.set_re(ant, sym, sc, dmrs_val * scale);
                            }
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
