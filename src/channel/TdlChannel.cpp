#include "channel/ChannelModels.h"
#include <cmath>
#include <random>
#include <algorithm>

namespace nr {
namespace channel {

TdlChannel::TdlChannel(ChannelType type) : type_(type), seed_(12345) {
    init_tdl_params();
}

void TdlChannel::set_seed(uint64_t seed) {
    seed_ = seed;
}

void TdlChannel::set_config(const SimulationConfig& config) {
    config_ = config;
}

void TdlChannel::init_tdl_params() {
    tdl_params_.n_taps = 0;
    tdl_params_.los = false;
    tdl_params_.k_factor_dB = 0.0;
    
    switch (type_) {
        case ChannelType::TDL_A:
            tdl_params_.n_taps = 23;
            tdl_params_.los = false;
            {
                double delays[] = {0.0, 0.3819, 0.4025, 0.5868, 0.4610, 0.5375, 0.6708, 0.5750,
                                   0.7618, 1.5375, 1.8978, 2.2242, 2.1718, 2.4942, 2.5119,
                                   3.0582, 4.0810, 4.4579, 4.5695, 4.7966, 5.0066, 5.3043, 9.6586};
                double powers[] = {-13.4, 0.0, -2.2, -4.0, -6.0, -8.2, -9.9, -10.5,
                                   -7.5, -15.9, -6.6, -16.7, -12.4, -15.2, -10.8,
                                   -11.3, -12.7, -16.2, -18.3, -18.9, -16.6, -19.9, -29.7};
                for (int i = 0; i < 23; i++) {
                    tdl_params_.norm_delay[i] = delays[i];
                    tdl_params_.norm_pwr[i] = std::pow(10.0, powers[i] / 10.0);
                }
            }
            break;
        case ChannelType::TDL_B:
            tdl_params_.n_taps = 23;
            tdl_params_.los = false;
            {
                double delays[] = {0.0, 0.1072, 0.2155, 0.2095, 0.2870, 0.2986, 0.3752, 0.5055,
                                   0.3681, 0.3697, 0.5700, 0.5283, 1.1021, 1.2756, 1.5474,
                                   1.7842, 2.0169, 2.8294, 3.0219, 3.6187, 4.1067, 4.2790, 4.7834};
                double powers[] = {0.0, -2.2, -4.0, -3.2, -9.8, -3.2, -3.4, -7.2,
                                   -3.4, -3.0, -5.0, -4.8, -10.4, -13.6, -16.0,
                                   -18.2, -17.3, -19.2, -16.7, -21.2, -23.8, -23.2, -28.4};
                for (int i = 0; i < 23; i++) {
                    tdl_params_.norm_delay[i] = delays[i];
                    tdl_params_.norm_pwr[i] = std::pow(10.0, powers[i] / 10.0);
                }
            }
            break;
        case ChannelType::TDL_C:
            tdl_params_.n_taps = 24;
            tdl_params_.los = false;
            {
                double delays[] = {0.0, 0.2099, 0.2219, 0.2329, 0.2176, 0.6366, 0.6448, 0.6560,
                                   0.6584, 0.7935, 0.8213, 0.9336, 1.2285, 1.3083, 2.1704,
                                   2.7104, 4.2589, 4.5998, 4.6001, 4.7966, 5.0066, 5.3043, 9.6586, 10.0000};
                double powers[] = {-4.4, -1.2, -3.5, -5.2, -2.5, 0.0, -2.2, -3.9,
                                   -7.4, -7.1, -10.7, -11.1, -5.1, -6.8, -8.7,
                                   -13.2, -20.1, -19.0, -17.6, -21.2, -22.9, -19.9, -29.7, -28.0};
                for (int i = 0; i < 24; i++) {
                    tdl_params_.norm_delay[i] = delays[i];
                    tdl_params_.norm_pwr[i] = std::pow(10.0, powers[i] / 10.0);
                }
            }
            break;
        case ChannelType::TDL_D:
        case ChannelType::TDL_E:
            tdl_params_.n_taps = 13;
            tdl_params_.los = true;
            tdl_params_.k_factor_dB = (type_ == ChannelType::TDL_D) ? 13.3 : 22.0;
            {
                double delays[] = {0.0, 0.035, 0.612, 1.363, 1.405, 1.804, 2.596,
                                   1.775, 4.042, 7.937, 9.424, 9.708, 12.525};
                double powers[] = {-13.4, 0.0, -2.2, -4.0, -6.0, -8.2, -9.9,
                                   -10.5, -15.9, -16.7, -12.4, -15.2, -18.3};
                for (int i = 0; i < 13; i++) {
                    tdl_params_.norm_delay[i] = delays[i];
                    tdl_params_.norm_pwr[i] = std::pow(10.0, powers[i] / 10.0);
                }
            }
            break;
        default:
            tdl_params_.n_taps = 1;
            tdl_params_.norm_delay[0] = 0.0;
            tdl_params_.norm_pwr[0] = 1.0;
            break;
    }
}

ComplexCube TdlChannel::get_channel(int n_prbs, int n_symbols, double sample_rate) {
    int n_sc = n_prbs * 12;
    int n_rx_ant = config_.n_rx_ant;
    int n_tx_ant = std::min(config_.n_tx_ant, config_.n_layers);
    int n_layers = config_.n_layers;
    
    ComplexCube h(n_sc, n_symbols, n_rx_ant * n_layers, arma::fill::zeros);
    
    std::mt19937 rng(seed_);
    std::normal_distribution<double> dist(0.0, 1.0);
    
    double max_delay = 0.0;
    for (int i = 0; i < tdl_params_.n_taps; i++) {
        max_delay = std::max(max_delay, tdl_params_.norm_delay[i] * config_.delay_spread);
    }
    
    for (int sym = 0; sym < n_symbols; sym++) {
        for (int rx = 0; rx < n_rx_ant; rx++) {
            for (int l = 0; l < n_layers; l++) {
                ComplexMat tap_coeffs(tdl_params_.n_taps, 1);
                double pwr_sum = 0.0;
                for (int t = 0; t < tdl_params_.n_taps; t++) {
                    pwr_sum += tdl_params_.norm_pwr[t];
                }
                
                for (int t = 0; t < tdl_params_.n_taps; t++) {
                    double re = dist(rng) * std::sqrt(tdl_params_.norm_pwr[t] / (2.0 * pwr_sum));
                    double im = dist(rng) * std::sqrt(tdl_params_.norm_pwr[t] / (2.0 * pwr_sum));
                    tap_coeffs(t, 0) = Complex(re, im);
                }
                
                if (tdl_params_.los && rx == l) {
                    double k_lin = std::pow(10.0, tdl_params_.k_factor_dB / 10.0);
                    double los_scale = std::sqrt(k_lin / (k_lin + 1.0));
                    double nlos_scale = std::sqrt(1.0 / (k_lin + 1.0));
                    tap_coeffs(0, 0) = tap_coeffs(0, 0) * nlos_scale + Complex(los_scale, 0.0);
                }
                
                for (int sc = 0; sc < n_sc; sc++) {
                    double freq = sc * config_.scs;
                    Complex h_freq(0, 0);
                    for (int t = 0; t < tdl_params_.n_taps; t++) {
                        double delay = tdl_params_.norm_delay[t] * config_.delay_spread;
                        double phase = -2.0 * M_PI * freq * delay;
                        h_freq += tap_coeffs(t, 0) * Complex(std::cos(phase), std::sin(phase));
                    }
                    h(sc, sym, rx * n_layers + l) = h_freq;
                }
            }
        }
    }
    
    return h;
}

ComplexVec TdlChannel::apply_channel(const ComplexVec& tx_signal, const ComplexCube& h) {
    return tx_signal;
}

double TdlChannel::get_thermal_noise(double sinr_db, double bandwidth, double noise_figure) {
    double sinr_linear = std::pow(10.0, sinr_db / 10.0);
    return 1.0 / sinr_linear;
}

ComplexVec TdlChannel::add_noise(const ComplexVec& signal, double noise_var) {
    std::mt19937 rng(seed_++);
    std::normal_distribution<double> dist(0.0, std::sqrt(noise_var / 2.0));
    
    ComplexVec noisy(signal.n_elem);
    for (int i = 0; i < signal.n_elem; i++) {
        double noise_re = dist(rng);
        double noise_im = dist(rng);
        noisy(i) = signal(i) + Complex(noise_re, noise_im);
    }
    
    return noisy;
}

std::string TdlChannel::get_name() const {
    switch (type_) {
        case ChannelType::TDL_A: return "TDL-A";
        case ChannelType::TDL_B: return "TDL-B";
        case ChannelType::TDL_C: return "TDL-C";
        case ChannelType::TDL_D: return "TDL-D";
        case ChannelType::TDL_E: return "TDL-E";
        default: return "TDL";
    }
}

} // namespace channel
} // namespace nr
