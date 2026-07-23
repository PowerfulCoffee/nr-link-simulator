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
            estimated_noise_var_ = 1e-6;
            return h_est;
        }

        double noise_power_sum = 0.0;
        int noise_sample_count = 0;
        constexpr double MIN_NOISE_VAR = 1e-10;
        // DMRS Type1 power boost: beta=sqrt(2), so |x|^2=2, h_ls noise variance = sigma^2/2
        // Second-order difference (internal SC): E[|noise|^2] = 3/2 * (sigma^2/2) = 3/4 sigma^2 -> scale by 4/3
        // First-order difference (edge SC): E[|noise|^2] = 2 * (sigma^2/2) = sigma^2 -> scale by 1.0
        constexpr double SCALE_INTERNAL = 4.0 / 3.0;
        constexpr double SCALE_EDGE = 1.0;

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

                    if (dmrs_sc.size() >= 2) {
                        {
                            int sc_first = dmrs_sc[0];
                            int sc_next = dmrs_sc[1];
                            Complex h_first = h_est(sc_first, dmrs_sym, ch_idx);
                            Complex h_next = h_est(sc_next, dmrs_sym, ch_idx);
                            Complex noise_sample = h_next - h_first;
                            noise_power_sum += SCALE_EDGE * std::norm(noise_sample);
                            noise_sample_count++;
                        }

                        for (size_t i = 1; i + 1 < dmrs_sc.size(); i++) {
                            int sc_prev = dmrs_sc[i - 1];
                            int sc_curr = dmrs_sc[i];
                            int sc_next = dmrs_sc[i + 1];
                            Complex h_prev = h_est(sc_prev, dmrs_sym, ch_idx);
                            Complex h_curr = h_est(sc_curr, dmrs_sym, ch_idx);
                            Complex h_next = h_est(sc_next, dmrs_sym, ch_idx);
                            Complex noise_sample = (h_prev + h_next) / 2.0 - h_curr;
                            noise_power_sum += SCALE_INTERNAL * std::norm(noise_sample);
                            noise_sample_count++;
                        }

                        {
                            int sc_last = dmrs_sc.back();
                            int sc_prev = dmrs_sc[dmrs_sc.size() - 2];
                            Complex h_last = h_est(sc_last, dmrs_sym, ch_idx);
                            Complex h_prev = h_est(sc_prev, dmrs_sym, ch_idx);
                            Complex noise_sample = h_prev - h_last;
                            noise_power_sum += SCALE_EDGE * std::norm(noise_sample);
                            noise_sample_count++;
                        }
                    }

                    int cdm_group_size = 2;
                    for (size_t i = 0; i + 1 < dmrs_sc.size(); i += cdm_group_size) {
                        size_t j = i + 1;
                        if (j >= dmrs_sc.size()) break;
                        int sc_a = dmrs_sc[i];
                        int sc_b = dmrs_sc[j];
                        if (sc_b - sc_a > 2) continue;
                        Complex avg = (h_est(sc_a, dmrs_sym, ch_idx) + h_est(sc_b, dmrs_sym, ch_idx)) / 2.0;
                        h_est(sc_a, dmrs_sym, ch_idx) = avg;
                        h_est(sc_b, dmrs_sym, ch_idx) = avg;
                    }
                    
                    interpolate_frequency(h_est, dmrs_sc, dmrs_sym, ch_idx, n_sc);
                }
            }
        }

        if (noise_sample_count > 0) {
            double avg_noise_power = noise_power_sum / noise_sample_count;
            estimated_noise_var_ = std::max(avg_noise_power, MIN_NOISE_VAR);
        } else {
            estimated_noise_var_ = MIN_NOISE_VAR;
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

    double get_estimated_noise_var() const override {
        return estimated_noise_var_;
    }

private:
    double estimated_noise_var_ = 0.0;
};

class PerfectChannelEstimator : public IChannelEstimator {
public:
    ComplexCube estimate(const ResourceGrid& rx_grid, const ResourceGrid& /*dmrs_grid*/,
                         const SimulationConfig& config) override {
        int n_sc = rx_grid.n_subcarriers;
        int n_sym = rx_grid.n_symbols;
        int n_rx_ant = rx_grid.n_ant;
        int n_layers = config.n_layers;
        int n_ch = n_rx_ant * n_layers;

        ComplexCube h_est(n_sc, n_sym, n_ch, arma::fill::zeros);

        if (perfect_h_.n_rows == static_cast<arma::uword>(n_sc) &&
            perfect_h_.n_cols == static_cast<arma::uword>(n_sym) &&
            perfect_h_.n_slices == static_cast<arma::uword>(n_ch)) {
            h_est = perfect_h_;
        } else {
            for (int rx = 0; rx < n_rx_ant; rx++) {
                for (int l = 0; l < n_layers; l++) {
                    int ch_idx = rx * n_layers + l;
                    for (int sc = 0; sc < n_sc; sc++) {
                        for (int sym = 0; sym < n_sym; sym++) {
                            if (rx == l) {
                                h_est(sc, sym, ch_idx) = Complex(1.0, 0.0);
                            }
                        }
                    }
                }
            }
        }

        noise_var_ = perfect_noise_var_;
        return h_est;
    }

    std::string get_name() const override {
        return "Perfect";
    }

    double get_estimated_noise_var() const override {
        return noise_var_;
    }

    void set_perfect_channel(const ComplexCube& h, double noise_var) override {
        perfect_h_ = h;
        perfect_noise_var_ = noise_var;
    }

private:
    ComplexCube perfect_h_;
    double perfect_noise_var_ = 0.0;
    double noise_var_ = 0.0;
};

std::unique_ptr<IChannelEstimator> create_ls_channel_estimator() {
    return std::make_unique<LsChannelEstimator>();
}

std::unique_ptr<IChannelEstimator> create_perfect_channel_estimator() {
    return std::make_unique<PerfectChannelEstimator>();
}

} // namespace phy
} // namespace nr
