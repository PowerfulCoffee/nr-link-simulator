#include "phy/PhyInterfaces.h"
#include <cmath>
#include <vector>
#include <algorithm>

namespace nr {
namespace phy {

namespace {

int get_fft_size_from_sc(int n_sc) {
    int fft_size = 1;
    while (fft_size < n_sc + 1) fft_size <<= 1;
    return fft_size;
}

int get_cp_len(int fft_size, bool is_first) {
    int cp = 144 * fft_size / 2048;
    return is_first ? cp + 16 * fft_size / 2048 : cp;
}

double get_ofdm_symbol_time(double scs, int fft_size, int sym_idx) {
    double tu = 1.0 / scs;
    double t = 0.0;
    for (int s = 0; s < sym_idx; s++) {
        bool is_first = (s == 0);
        int cp = get_cp_len(fft_size, is_first);
        t += tu * (fft_size + cp) / fft_size;
    }
    bool is_first = (sym_idx == 0);
    int cp = get_cp_len(fft_size, is_first);
    t += tu * (fft_size + cp / 2) / fft_size;
    return t;
}

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
        
        int n_ch = n_rx_ant * n_tx_ant;  // per-port channels (rx x tx)
        ComplexCube h_est(n_sc, n_sym, n_ch, arma::fill::zeros);

        std::vector<int> dmrs_symbols = find_dmrs_symbols(dmrs_grid, n_tx_ant);
        if (dmrs_symbols.empty()) {
            for (int rx = 0; rx < n_rx_ant; rx++) {
                for (int tx = 0; tx < n_tx_ant; tx++) {
                    int ch_idx = rx * n_tx_ant + tx;
                    for (int sc = 0; sc < n_sc; sc++) {
                        for (int sym = 0; sym < n_sym; sym++) {
                            h_est(sc, sym, ch_idx) = Complex(1.0, 0.0);
                        }
                    }
                }
            }
            estimated_noise_var_ = 1e-6;
            estimated_doppler_ = 0.0;
            return h_est;
        }

        double noise_power_sum = 0.0;
        int noise_sample_count = 0;
        constexpr double MIN_NOISE_VAR = 1e-10;
        // Edge two-point difference: s = h_next - h_prev, Var(s) = 2*sigma^2
        constexpr double SCALE_EDGE = 0.5;

        // ---- Type1 CDM-group pairing (TS 38.211 7.4.1.1.2 freq OCC) ----
        // Ports p and p+2 of a CDM group share the same DMRS RE pairs; the
        // OCC w_f = [1, +/-1] across the pair yields, per RE pair,
        //   y(first)/r0  = h_p + h_{p+2}
        //   y(second)/r1 = h_p - h_{p+2}
        // which is jointly solvable for both ports. Without the joint solve
        // a per-port LS returns inseparable sums and MIMO-LS collapses.
        // Enabled only when the paired port exists (>=3 layers); SISO and
        // 2-layer runs keep the legacy per-port LS path below unchanged.
        const bool is_type1_dmrs = (config.dmrs_type == DmrsType::TYPE1);
        const int max_ports = is_type1_dmrs ? 4 : 6;
        const int n_gen_layers = std::min(config.n_layers, max_ports);
        const bool pair_mode = is_type1_dmrs && n_gen_layers >= 3;

        for (int dmrs_sym : dmrs_symbols) {
            for (int rx = 0; rx < n_rx_ant; rx++) {
                if (pair_mode) {
                    auto accumulate_noise = [&](const std::vector<int>& scs, int ch, double scale) {
                        if (scs.size() < 2) return;
                        {
                            Complex ns = h_est(scs[1], dmrs_sym, ch) - h_est(scs[0], dmrs_sym, ch);
                            noise_power_sum += SCALE_EDGE * std::norm(ns) * scale;
                            noise_sample_count++;
                        }
                        for (size_t i = 1; i + 1 < scs.size(); i++) {
                            double d1 = static_cast<double>(scs[i] - scs[i - 1]);
                            double d2 = static_cast<double>(scs[i + 1] - scs[i]);
                            if (d1 <= 0.0 || d2 <= 0.0) continue;
                            double d12 = d1 + d2;
                            Complex ns = (h_est(scs[i - 1], dmrs_sym, ch) * (d2 / d12)
                                        + h_est(scs[i + 1], dmrs_sym, ch) * (d1 / d12))
                                       - h_est(scs[i], dmrs_sym, ch);
                            double w2 = d2 / d12, w1 = d1 / d12;
                            noise_power_sum += std::norm(ns) * scale / (1.0 + w2 * w2 + w1 * w1);
                            noise_sample_count++;
                        }
                        {
                            size_t n = scs.size();
                            Complex ns = h_est(scs[n - 2], dmrs_sym, ch) - h_est(scs[n - 1], dmrs_sym, ch);
                            noise_power_sum += SCALE_EDGE * std::norm(ns) * scale;
                            noise_sample_count++;
                        }
                    };

                    for (int g = 0; g < 2; g++) {
                        int tx_base = g;
                        if (tx_base >= n_gen_layers) continue;
                        int tx_pair = g + 2;
                        const bool has_pair = (tx_pair < n_gen_layers);
                        const int ch_base = rx * n_tx_ant + tx_base;
                        const int ch_pair = rx * n_tx_ant + tx_pair;

                        std::vector<int> dmrs_sc = find_dmrs_subcarriers(dmrs_grid, dmrs_sym, tx_base);

                        for (size_t i = 0; i + 1 < dmrs_sc.size(); i += 2) {
                            int s0 = dmrs_sc[i], s1 = dmrs_sc[i + 1];
                            Complex ref0 = dmrs_grid.get_re(tx_base, dmrs_sym, s0);
                            Complex ref1 = dmrs_grid.get_re(tx_base, dmrs_sym, s1);
                            Complex y0 = rx_grid.get_re(rx, dmrs_sym, s0);
                            Complex y1 = rx_grid.get_re(rx, dmrs_sym, s1);
                            Complex sv = (std::abs(ref0) > 1e-10) ? y0 / ref0 : Complex(0.0, 0.0);
                            Complex dv = (std::abs(ref1) > 1e-10) ? y1 / ref1 : Complex(0.0, 0.0);
                            if (has_pair) {
                                Complex h_b = 0.5 * (sv + dv);
                                Complex h_p = 0.5 * (sv - dv);
                                h_est(s0, dmrs_sym, ch_base) = h_b;
                                h_est(s1, dmrs_sym, ch_base) = h_b;
                                h_est(s0, dmrs_sym, ch_pair) = h_p;
                                h_est(s1, dmrs_sym, ch_pair) = h_p;
                            } else {
                                h_est(s0, dmrs_sym, ch_base) = sv;
                                h_est(s1, dmrs_sym, ch_base) = dv;
                            }
                        }
                        if (dmrs_sc.size() % 2 == 1) { // defensive odd tail
                            int s0 = dmrs_sc.back();
                            Complex ref0 = dmrs_grid.get_re(tx_base, dmrs_sym, s0);
                            Complex y0 = rx_grid.get_re(rx, dmrs_sym, s0);
                            if (std::abs(ref0) > 1e-10) h_est(s0, dmrs_sym, ch_base) = y0 / ref0;
                        }

                        // joint LS writes the same value to both REs of a pair,
                        // so thin to one sample per pair (value is N(h, sigma^2/2))
                        // for the differencing statistic; otherwise the repeated
                        // values shrink the variance estimate by 3x.
                        std::vector<int> scs_stat;
                        if (has_pair) {
                            scs_stat.reserve((dmrs_sc.size() + 1) / 2);
                            for (size_t i = 0; i < dmrs_sc.size(); i += 2) {
                                scs_stat.push_back(dmrs_sc[i]);
                            }
                        } else {
                            scs_stat = dmrs_sc;
                        }
                        accumulate_noise(scs_stat, ch_base, has_pair ? 2.0 : 1.0);
                        if (has_pair) accumulate_noise(scs_stat, ch_pair, 2.0);

                        interpolate_frequency(h_est, dmrs_sc, dmrs_sym, ch_base, n_sc);
                        if (has_pair) interpolate_frequency(h_est, dmrs_sc, dmrs_sym, ch_pair, n_sc);
                    }
                    continue; // pair path done for this (dmrs_sym, rx)
                }
                for (int tx = 0; tx < n_tx_ant; tx++) {
                    int ch_idx = rx * n_tx_ant + tx;

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

                    // Noise estimation based on the actual DMRS subcarrier
                    // spacing (comb-2 for Type 1: d = 2; Type 2 has uneven
                    // gaps). Each LS estimate h = h + n has Var(n) = sigma^2
                    // (unit-power DMRS), independent of spacing; only the
                    // smoothing weights depend on the gaps.
                    if (dmrs_sc.size() >= 2) {
                        // Edge first difference
                        {
                            int sc_first = dmrs_sc[0];
                            int sc_next = dmrs_sc[1];
                            Complex noise_sample = h_est(sc_next, dmrs_sym, ch_idx)
                                                 - h_est(sc_first, dmrs_sym, ch_idx);
                            noise_power_sum += SCALE_EDGE * std::norm(noise_sample);
                            noise_sample_count++;
                        }

                        // Interior three-point linear smoothing with spacing
                        // (d1, d2): s = (d2*h_prev + d1*h_next)/(d1+d2) - h_curr
                        // Var(s) = sigma^2 * (1 + (d2/(d1+d2))^2 + (d1/(d1+d2))^2)
                        // For comb-2 (d1 = d2 = 2): Var = 3/2*sigma^2 -> scale 2/3
                        for (size_t i = 1; i + 1 < dmrs_sc.size(); i++) {
                            int sc_prev = dmrs_sc[i - 1];
                            int sc_curr = dmrs_sc[i];
                            int sc_next = dmrs_sc[i + 1];
                            double d1 = static_cast<double>(sc_curr - sc_prev);
                            double d2 = static_cast<double>(sc_next - sc_curr);
                            double d12 = d1 + d2;
                            if (d1 <= 0.0 || d2 <= 0.0) continue;
                            Complex h_prev = h_est(sc_prev, dmrs_sym, ch_idx);
                            Complex h_curr = h_est(sc_curr, dmrs_sym, ch_idx);
                            Complex h_next = h_est(sc_next, dmrs_sym, ch_idx);
                            Complex noise_sample = (h_prev * (d2 / d12) + h_next * (d1 / d12)) - h_curr;
                            double w2 = d2 / d12, w1 = d1 / d12;
                            double var_factor = 1.0 + w2 * w2 + w1 * w1;
                            noise_power_sum += std::norm(noise_sample) / var_factor;
                            noise_sample_count++;
                        }

                        // Edge last difference
                        {
                            int sc_last = dmrs_sc.back();
                            int sc_prev = dmrs_sc[dmrs_sc.size() - 2];
                            Complex noise_sample = h_est(sc_prev, dmrs_sym, ch_idx)
                                                 - h_est(sc_last, dmrs_sym, ch_idx);
                            noise_power_sum += SCALE_EDGE * std::norm(noise_sample);
                            noise_sample_count++;
                        }
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

        int fft_size = get_fft_size_from_sc(n_sc);
        double scs = config.scs;
        estimated_doppler_ = 0.0;

        if (dmrs_symbols.size() >= 2) {
            Complex cross_corr(0.0, 0.0);
            double ref_power = 0.0;
            int sym_a = dmrs_symbols[0];
            int sym_b = dmrs_symbols[1];
            double t_a = get_ofdm_symbol_time(scs, fft_size, sym_a);
            double t_b = get_ofdm_symbol_time(scs, fft_size, sym_b);
            double delta_t = t_b - t_a;

            for (int rx = 0; rx < n_rx_ant; rx++) {
                for (int tx = 0; tx < n_tx_ant; tx++) {
                    int ch_idx = rx * n_tx_ant + tx;
                    for (int sc = 0; sc < n_sc; sc++) {
                        Complex h_a = h_est(sc, sym_a, ch_idx);
                        Complex h_b = h_est(sc, sym_b, ch_idx);
                        cross_corr += h_b * std::conj(h_a);
                        ref_power += std::norm(h_a);
                    }
                }
            }

            if (std::abs(cross_corr) > 1e-10 && delta_t > 1e-10) {
                double phase_diff = std::arg(cross_corr);
                estimated_doppler_ = phase_diff / (2.0 * M_PI * delta_t);
            }
        }

        for (int rx = 0; rx < n_rx_ant; rx++) {
            for (int tx = 0; tx < n_tx_ant; tx++) {
                int ch_idx = rx * n_tx_ant + tx;
                
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
                    if (dist <= 0) continue;
                    
                    for (int sym = sym_a + 1; sym < sym_b; sym++) {
                        double alpha = static_cast<double>(sym - sym_a) / dist;
                        for (int sc = 0; sc < n_sc; sc++) {
                            h_est(sc, sym, ch_idx) = h_est(sc, sym_a, ch_idx) * (1.0 - alpha) 
                                                   + h_est(sc, sym_b, ch_idx) * alpha;
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
        return "LS-Doppler";
    }

    double get_estimated_noise_var() const override {
        return estimated_noise_var_;
    }
    
    double get_estimated_doppler() const {
        return estimated_doppler_;
    }

private:
    double estimated_noise_var_ = 0.0;
    double estimated_doppler_ = 0.0;
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
