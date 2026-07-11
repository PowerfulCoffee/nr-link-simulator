#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <vector>
#include <cmath>

namespace nr {
namespace phy {

class DmrsGeneratorImpl : public IDmrsGenerator {
public:
    void generate_dmrs(const SimulationConfig& config, ResourceGrid& grid,
                       int slot_idx, int symbol_start) override {
        DmrsPattern pattern = get_dmrs_pattern(config.dmrs_type, config.dmrs_additional_pos, 
                                               config.dmrs_duration);
        
        int n_sc = grid.n_subcarriers;
        int n_sym = grid.n_symbols;
        int n_ant = grid.n_ant;
        
        std::vector<uint8_t> pn_seq;
        uint32_t cinit = (slot_idx << 14) | (symbol_start << 12) | (config.n_rb & 0x3FF);
        generate_pn_sequence(cinit, n_sc * 2, pn_seq);
        
        for (int sym = 0; sym < n_sym; sym++) {
            if (pattern.re_per_prb[sym] == 0) continue;
            
            for (int prb = 0; prb < config.n_rb; prb++) {
                for (int re_in_prb = 0; re_in_prb < 12; re_in_prb++) {
                    if (!is_dmrs_re(config.dmrs_type, prb, re_in_prb, sym, symbol_start)) continue;
                    
                    int sc = prb * 12 + re_in_prb;
                    if (sc >= n_sc) continue;
                    
                    int seq_idx = (sc * 2) % pn_seq.size();
                    double re = (pn_seq[seq_idx] == 0) ? 1.0 / std::sqrt(2.0) : -1.0 / std::sqrt(2.0);
                    double im = (pn_seq[seq_idx + 1] == 0) ? 1.0 / std::sqrt(2.0) : -1.0 / std::sqrt(2.0);
                    Complex dmrs_sym(re, im);
                    
                    for (int ant = 0; ant < std::min(n_ant, 2); ant++) {
                        double phase = ant * M_PI * sc / n_sc;
                        grid.set_re(ant, sym, sc, dmrs_sym * Complex(std::cos(phase), std::sin(phase)));
                    }
                }
            }
        }
    }
    
    ComplexCube extract_dmrs(const ResourceGrid& rx_grid, const SimulationConfig& config,
                             int slot_idx, int symbol_start) override {
        int n_sc = rx_grid.n_subcarriers;
        int n_sym = rx_grid.n_symbols;
        int n_rx_ant = rx_grid.n_ant;
        int n_layers = config.n_layers;
        
        ComplexCube dmrs(n_sc, n_sym, n_rx_ant * n_layers, arma::fill::zeros);
        
        for (int sym = 0; sym < n_sym; sym++) {
            for (int sc = 0; sc < n_sc; sc++) {
                int prb = sc / 12;
                int re_in_prb = sc % 12;
                
                if (!is_dmrs_re(config.dmrs_type, prb, re_in_prb, sym, symbol_start)) continue;
                
                for (int rx = 0; rx < n_rx_ant; rx++) {
                    Complex val = rx_grid.get_re(rx, sym, sc);
                    for (int l = 0; l < n_layers; l++) {
                        dmrs(sc, sym, rx * n_layers + l) = val;
                    }
                }
            }
        }
        
        return dmrs;
    }

private:
    bool is_dmrs_re(DmrsType type, int prb, int re_in_prb, int sym, int sym_start) {
        if (type == DmrsType::TYPE1) {
            int dmrs_sym = sym_start;
            if (sym != dmrs_sym && sym != dmrs_sym + 1) return false;
            
            if (re_in_prb == 0 || re_in_prb == 2 || re_in_prb == 4 || re_in_prb == 6 ||
                re_in_prb == 8 || re_in_prb == 10) {
                return true;
            }
        } else {
            int dmrs_sym = sym_start;
            if (sym != dmrs_sym && sym != dmrs_sym + 1) return false;
            
            if (re_in_prb == 0 || re_in_prb == 1 || re_in_prb == 6 || re_in_prb == 7) {
                return true;
            }
        }
        return false;
    }
};

std::unique_ptr<IDmrsGenerator> create_dmrs_generator() {
    return std::make_unique<DmrsGeneratorImpl>();
}

} // namespace phy
} // namespace nr
