#include "phy/PhyInterfaces.h"
#include <cmath>
#include <vector>
#include <algorithm>

namespace nr {
namespace phy {

namespace {

std::vector<int> find_dmrs_symbols(const ResourceGrid& dmrs_grid, int n_tx_ant) {
    std::vector<int> dmrs_symbols;
    int n_sym = dmrs_grid.n_symbols;
    int n_sc = dmrs_grid.n_subcarriers;
    
    for (int sym = 0; sym < n_sym; sym++) {
        bool has_dmrs = false;
        for (int sc = 0; sc < n_sc; sc++) {
            for (int tx = 0; tx < n_tx_ant; tx++) {
                if (std::abs(dmrs_grid.get_re(tx, sym, sc)) > 1e-10) {
                    has_dmrs = true;
                    break;
                }
            }
            if (has_dmrs) break;
        }
        if (has_dmrs) {
            dmrs_symbols.push_back(sym);
        }
    }
    return dmrs_symbols;
}

std::vector<int> find_dmrs_subcarriers(const ResourceGrid& dmrs_grid, int sym, int tx) {
    std::vector<int> sc_list;
    int n_sc = dmrs_grid.n_subcarriers;
    for (int sc = 0; sc < n_sc; sc++) {
        if (std::abs(dmrs_grid.get_re(tx, sym, sc)) > 1e-10) {
            sc_list.push_back(sc);
        }
    }
    return sc_list;
}

void interpolate_frequency(ComplexCube& h_est, const std::vector<int>& dmrs_sc,
                           int sym, int ch_idx, int n_sc) {
    if (dmrs_sc.empty()) return;
    
    std::vector<Complex> h_at_dmrs(dmrs_sc.size());
    for (size_t i = 0; i < dmrs_sc.size(); i++) {
        h_at_dmrs[i] = h_est(dmrs_sc[i], sym, ch_idx);
    }
    
    int first_dmrs_sc = dmrs_sc.front();
    int last_dmrs_sc = dmrs_sc.back();
    
    for (int sc = 0; sc < first_dmrs_sc; sc++) {
        h_est(sc, sym, ch_idx) = h_at_dmrs.front();
    }
    
    for (size_t i = 0; i + 1 < dmrs_sc.size(); i++) {
        int sc_a = dmrs_sc[i];
        int sc_b = dmrs_sc[i + 1];
        Complex h_a = h_at_dmrs[i];
        Complex h_b = h_at_dmrs[i + 1];
        int dist = sc_b - sc_a;
        if (dist <= 0) continue;
        
        for (int sc = sc_a; sc <= sc_b; sc++) {
            double alpha = static_cast<double>(sc - sc_a) / dist;
            h_est(sc, sym, ch_idx) = h_a * (1.0 - alpha) + h_b * alpha;
        }
    }
    
    for (int sc = last_dmrs_sc; sc < n_sc; sc++) {
        h_est(sc, sym, ch_idx) = h_at_dmrs.back();
    }
}

} // anonymous namespace

class LsChannelEstimator : public IChannelEstimator {
public:
    ComplexCube estimate(const ResourceGrid& rx_grid, const ResourceGrid& dmrs_grid,
                         const SimulationConfig& config) override {
        int n_sc = rx_grid.n_subcarriers;
        int n_sym = rx_grid.n_symbols;
        int n_rx_ant = rx_grid.n_ant;
        int n_tx_ant = dmrs_grid.n_ant;
        int n_layers = config.n_layers;
        
        int n_ch = n_rx_ant * n_layers;
        ComplexCube h_est(n_sc, n_sym, n_ch, arma::fill::zeros);
        
        std::vector<int> dmrs_symbols = find_dmrs_symbols(dmrs_grid, n_tx_ant);
        if (dmrs_symbols.empty()) {
            for (int rx = 0; rx < n_rx_ant; rx++) {
                for (int l = 0; l < n_layers; l++) {
                    int ch_idx = rx * n_layers + l;
                    for (int sc = 0; sc < n_sc; sc++) {
                        for (int sym = 0; sym < n_sym; sym++) {
                            h_est(sc, sym, ch_idx) = Complex(1.0, 0.0);
                        }
                    }
                }
            }
            return h_est;
        }
        
        int ref_tx_ant = std::min(n_tx_ant, n_layers);
        for (int dmrs_sym : dmrs_symbols) {
            for (int rx = 0; rx < n_rx_ant; rx++) {
                for (int l = 0; l < n_layers; l++) {
                    int ch_idx = rx * n_layers + l;
                    int tx = l % ref_tx_ant;
                    
                    std::vector<int> dmrs_sc = find_dmrs_subcarriers(dmrs_grid, dmrs_sym, tx);
                    
                    for (int sc : dmrs_sc) {
                        Complex rx_val = rx_grid.get_re(rx, dmrs_sym, sc);
                        Complex ref_val = dmrs_grid.get_re(tx, dmrs_sym, sc);
                        
                        if (std::abs(ref_val) > 1e-10) {
                            h_est(sc, dmrs_sym, ch_idx) = rx_val / ref_val;
                        } else {
                            h_est(sc, dmrs_sym, ch_idx) = Complex(1.0, 0.0);
                        }
                    }
                    
                    interpolate_frequency(h_est, dmrs_sc, dmrs_sym, ch_idx, n_sc);
                }
            }
        }
        
        for (int rx = 0; rx < n_rx_ant; rx++) {
            for (int l = 0; l < n_layers; l++) {
                int ch_idx = rx * n_layers + l;
                
                int first_dmrs_sym = dmrs_symbols.front();
                int last_dmrs_sym = dmrs_symbols.back();
                
                for (int sym = 0; sym < first_dmrs_sym; sym++) {
                    for (int sc = 0; sc < n_sc; sc++) {
                        h_est(sc, sym, ch_idx) = h_est(sc, first_dmrs_sym, ch_idx);
                    }
                }
                
                for (size_t i = 0; i + 1 < dmrs_symbols.size(); i++) {
                    int sym_a = dmrs_symbols[i];
                    int sym_b = dmrs_symbols[i + 1];
                    int dist = sym_b - sym_a;
                    if (dist <= 1) continue;
                    
                    for (int sym = sym_a + 1; sym < sym_b; sym++) {
                        double alpha = static_cast<double>(sym - sym_a) / dist;
                        for (int sc = 0; sc < n_sc; sc++) {
                            Complex h_a = h_est(sc, sym_a, ch_idx);
                            Complex h_b = h_est(sc, sym_b, ch_idx);
                            h_est(sc, sym, ch_idx) = h_a * (1.0 - alpha) + h_b * alpha;
                        }
                    }
                }
                
                for (int sym = last_dmrs_sym + 1; sym < n_sym; sym++) {
                    for (int sc = 0; sc < n_sc; sc++) {
                        h_est(sc, sym, ch_idx) = h_est(sc, last_dmrs_sym, ch_idx);
                    }
                }
            }
        }
        
        return h_est;
    }
    
    std::string get_name() const override {
        return "LS";
    }
};

std::unique_ptr<IChannelEstimator> create_ls_channel_estimator() {
    return std::make_unique<LsChannelEstimator>();
}

} // namespace phy
} // namespace nr
