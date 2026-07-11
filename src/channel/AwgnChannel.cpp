#include "channel/ChannelModels.h"
#include <cmath>
#include <random>

namespace nr {
namespace channel {

AwgnChannel::AwgnChannel() : seed_(12345) {}

void AwgnChannel::set_seed(uint64_t seed) {
    seed_ = seed;
}

void AwgnChannel::set_config(const SimulationConfig& config) {
    config_ = config;
}

ComplexCube AwgnChannel::get_channel(int n_prbs, int n_symbols, double /*sample_rate*/) {
    int n_sc = n_prbs * 12;
    int n_rx = config_.n_rx_ant;
    int n_tx_layers = config_.n_layers;
    
    ComplexCube h(n_sc, n_symbols, n_rx * n_tx_layers, arma::fill::zeros);
    
    for (int sc = 0; sc < n_sc; sc++) {
        for (int sym = 0; sym < n_symbols; sym++) {
            for (int rx = 0; rx < n_rx; rx++) {
                for (int l = 0; l < n_tx_layers; l++) {
                    if (rx == l) {
                        h(sc, sym, rx * n_tx_layers + l) = Complex(1.0, 0.0);
                    }
                }
            }
        }
    }
    
    return h;
}

ComplexVec AwgnChannel::apply_channel(const ComplexVec& tx_signal, const ComplexCube& /*h*/) {
    return tx_signal;
}

double AwgnChannel::get_thermal_noise(double sinr_db, double /*bandwidth*/, double /*noise_figure*/) {
    double sinr_linear = std::pow(10.0, sinr_db / 10.0);
    double signal_power = 1.0;
    double noise_var_per_complex = signal_power / sinr_linear;
    return noise_var_per_complex;
}

ComplexVec AwgnChannel::add_noise(const ComplexVec& signal, double noise_var) {
    std::mt19937 rng(seed_++);
    double sigma_per_dim = std::sqrt(noise_var / 2.0);
    std::normal_distribution<double> dist(0.0, sigma_per_dim);
    
    ComplexVec noisy(signal.n_elem);
    for (arma::uword i = 0; i < signal.n_elem; i++) {
        double noise_re = dist(rng);
        double noise_im = dist(rng);
        noisy(i) = signal(i) + Complex(noise_re, noise_im);
    }
    
    return noisy;
}

} // namespace channel
} // namespace nr
