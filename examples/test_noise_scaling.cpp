#include "phy/PhyInterfaces.h"
#include "phy/ModuleFactory.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>

using namespace nr;
using namespace nr::phy;

int main() {
    const int n_prbs = 25;
    const int n_sc = n_prbs * 12;
    const int n_sym = 14;
    const int n_ant = 1;
    const double scs = 15e3;
    
    auto ofdm = create_ofdm_modulator();
    
    // Create a resource grid with unit-power QPSK symbols
    ResourceGrid grid(n_ant, n_sym, n_sc);
    std::mt19937 rng(42);
    std::normal_distribution<double> nd(0.0, 1.0);
    
    for (int ant = 0; ant < n_ant; ant++) {
        for (int sym = 0; sym < n_sym; sym++) {
            for (int sc = 0; sc < n_sc; sc++) {
                double i = (nd(rng) > 0) ? 1.0/M_SQRT2 : -1.0/M_SQRT2;
                double q = (nd(rng) > 0) ? 1.0/M_SQRT2 : -1.0/M_SQRT2;
                grid.set_re(ant, sym, sc, Complex(i, q));
            }
        }
    }
    
    // Measure frequency domain power
    double freq_power = 0.0;
    int n_freq_re = 0;
    for (int ant = 0; ant < n_ant; ant++) {
        for (int sym = 0; sym < n_sym; sym++) {
            for (int sc = 0; sc < n_sc; sc++) {
                Complex v = grid.get_re(ant, sym, sc);
                freq_power += std::norm(v);
                n_freq_re++;
            }
        }
    }
    freq_power /= n_freq_re;
    std::cout << "Freq-domain RE power (active): " << freq_power << std::endl;
    
    ComplexVec tx_signal = ofdm->modulate(grid, scs);
    int n_samp = tx_signal.n_elem;
    
    double time_power = 0.0;
    for (int i = 0; i < n_samp; i++) {
        time_power += std::norm(tx_signal(i));
    }
    time_power /= n_samp;
    std::cout << "Time-domain signal power: " << time_power << " (n_sc/fft_size = " << (double)n_sc/512.0 << ")\n";
    
    const double snr_db = 10.0;
    const double snr_lin = std::pow(10.0, snr_db/10.0);
    
    // Test 1: OLD method (noise_var_time = sig_pwr/SNR)
    std::cout << "\n=== OLD method (noise_var = sig_pwr/SNR = " << time_power/snr_lin << ") ===" << std::endl;
    {
        std::mt19937 test_rng(123);
        std::normal_distribution<double> dist(0.0, 1.0);
        double noise_var_time = time_power / snr_lin;
        ComplexVec rx = tx_signal;
        for (int i = 0; i < n_samp; i++) {
            double ni = dist(test_rng) * std::sqrt(noise_var_time / 2.0);
            double nq = dist(test_rng) * std::sqrt(noise_var_time / 2.0);
            rx(i) += Complex(ni, nq);
        }
        ResourceGrid demod = ofdm->demodulate(rx, n_ant, scs, n_sym);
        int n_sc_rx = demod.n_subcarriers;
        int offset = (n_sc_rx - n_sc) / 2;
        double rx_sig = 0, rx_noise = 0;
        int cnt = 0;
        for (int sym = 0; sym < n_sym; sym++) {
            for (int sc = 0; sc < n_sc; sc++) {
                Complex tx_v = grid.get_re(0, sym, sc);
                Complex rx_v = demod.get_re(0, sym, sc + offset);
                rx_sig += std::norm(tx_v);
                rx_noise += std::norm(rx_v - tx_v);
                cnt++;
            }
        }
        rx_sig /= cnt; rx_noise /= cnt;
        std::cout << "  After demod: sig_pwr=" << rx_sig << ", noise_pwr=" << rx_noise << std::endl;
        std::cout << "  Actual SNR: " << 10*std::log10(rx_sig/rx_noise) << " dB (expected 10 dB)" << std::endl;
        std::cout << "  Expected noise_var for receiver (1/SNR): " << 1.0/snr_lin << std::endl;
        std::cout << "  Mismatch factor: " << rx_noise/(1.0/snr_lin) << "x (causes LLR error)\n";
    }
    
    // Test 2: NEW method (noise_var_time = 1/SNR)
    std::cout << "\n=== NEW method (noise_var = 1/SNR = " << 1.0/snr_lin << ") ===" << std::endl;
    {
        std::mt19937 test_rng(123);
        std::normal_distribution<double> dist(0.0, 1.0);
        double noise_var_time = 1.0 / snr_lin;
        ComplexVec rx = tx_signal;
        for (int i = 0; i < n_samp; i++) {
            double ni = dist(test_rng) * std::sqrt(noise_var_time / 2.0);
            double nq = dist(test_rng) * std::sqrt(noise_var_time / 2.0);
            rx(i) += Complex(ni, nq);
        }
        ResourceGrid demod = ofdm->demodulate(rx, n_ant, scs, n_sym);
        int n_sc_rx = demod.n_subcarriers;
        int offset = (n_sc_rx - n_sc) / 2;
        double rx_sig = 0, rx_noise = 0;
        int cnt = 0;
        for (int sym = 0; sym < n_sym; sym++) {
            for (int sc = 0; sc < n_sc; sc++) {
                Complex tx_v = grid.get_re(0, sym, sc);
                Complex rx_v = demod.get_re(0, sym, sc + offset);
                rx_sig += std::norm(tx_v);
                rx_noise += std::norm(rx_v - tx_v);
                cnt++;
            }
        }
        rx_sig /= cnt; rx_noise /= cnt;
        std::cout << "  After demod: sig_pwr=" << rx_sig << ", noise_pwr=" << rx_noise << std::endl;
        std::cout << "  Actual SNR: " << 10*std::log10(rx_sig/rx_noise) << " dB (expected 10 dB)" << std::endl;
        std::cout << "  Matches receiver's expected noise_var = 1/SNR: " 
                  << (std::abs(rx_noise - 1.0/snr_lin) < 0.002 ? "YES" : "NO") << std::endl;
    }
    
    return 0;
}
