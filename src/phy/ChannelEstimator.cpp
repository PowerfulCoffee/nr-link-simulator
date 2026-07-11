#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <cmath>
#include <vector>

namespace nr {
namespace phy {

class LsChannelEstimator : public IChannelEstimator {
public:
    ComplexCube estimate(const ResourceGrid& rx_grid, const ResourceGrid& dmrs_grid,
                         const SimulationConfig& config) override {
        int n_sc = rx_grid.n_subcarriers;
        int n_sym = rx_grid.n_symbols;
        int n_rx_ant = rx_grid.n_ant;
        int n_tx_ant = config.n_tx_ant;
        int n_layers = config.n_layers;
        
        ComplexCube h_est(n_sc, n_sym, n_rx_ant * n_layers, arma::fill::zeros);
        
        std::vector<int> dmrs_symbols;
        for (int sym = 0; sym < n_sym; sym++) {
            bool has_dmrs = false;
            for (int sc = 0; sc < n_sc; sc++) {
                for (int tx = 0; tx < std::min(2, n_tx_ant); tx++) {
                    if (std::abs(dmrs_grid.get_re(tx, sym, sc)) > 0.1) {
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
        
        if (dmrs_symbols.empty()) {
            dmrs_symbols.push_back(0);
        }
        
        for (int dmrs_sym : dmrs_symbols) {
            for (int sc = 0; sc < n_sc; sc++) {
                for (int rx = 0; rx < n_rx_ant; rx++) {
                    Complex rx_val = rx_grid.get_re(rx, dmrs_sym, sc);
                    
                    for (int l = 0; l < n_layers; l++) {
                        Complex ref_val(0, 0);
                        int ref_ant = l % std::min(2, n_tx_ant);
                        ref_val = dmrs_grid.get_re(ref_ant, dmrs_sym, sc);
                        
                        Complex h_ls(0, 0);
                        if (std::abs(ref_val) > 1e-10) {
                            h_ls = rx_val / ref_val;
                        }
                        
                        h_est(sc, dmrs_sym, rx * n_layers + l) = h_ls;
                    }
                }
            }
        }
        
        for (int l = 0; l < n_layers; l++) {
            for (int rx = 0; rx < n_rx_ant; rx++) {
                int ch_idx = rx * n_layers + l;
                
                for (int dmrs_i = 0; dmrs_i < (int)dmrs_symbols.size(); dmrs_i++) {
                    int sym_a = dmrs_symbols[dmrs_i];
                    int sym_b = (dmrs_i + 1 < (int)dmrs_symbols.size()) ? 
                                dmrs_symbols[dmrs_i + 1] : n_sym;
                    
                    for (int sym = sym_a; sym < sym_b; sym++) {
                        double alpha = 0.0;
                        if (sym_b < n_sym) {
                            alpha = static_cast<double>(sym - sym_a) / (sym_b - sym_a);
                        }
                        
                        for (int sc = 0; sc < n_sc; sc++) {
                            if (sym == sym_a) {
                                h_est(sc, sym, ch_idx) = h_est(sc, sym_a, ch_idx);
                            } else {
                                Complex h_a = h_est(sc, sym_a, ch_idx);
                                Complex h_b = (sym_b < n_sym) ? h_est(sc, sym_b, ch_idx) : h_a;
                                h_est(sc, sym, ch_idx) = h_a * (1.0 - alpha) + h_b * alpha;
                            }
                        }
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

std::unique_ptr<IChannelEstimator> create_channel_estimator(const std::string& name) {
    if (name == "LS" || name == "ls") {
        return create_ls_channel_estimator();
    }
    return create_ls_channel_estimator();
}

} // namespace phy
} // namespace nr
