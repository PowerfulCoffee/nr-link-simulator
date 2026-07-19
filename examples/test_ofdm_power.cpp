#include "phy/PhyInterfaces.h"
#include "phy/ModuleFactory.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <random>
#include <armadillo>

using namespace nr;
using namespace nr::phy;

int main() {
    const int n_prb = 25;
    const int n_sc = n_prb * 12;  // 300
    const int n_sym = 14;
    const int n_ant = 1;
    
    // Compute FFT size (same as OfdmModulator)
    int fft_size = 1;
    while (fft_size < n_sc + 1) fft_size <<= 1;
    
    std::cout << "n_sc=" << n_sc << ", fft_size=" << fft_size << std::endl;
    std::cout << "n_active/n_total = " << (double)n_sc / fft_size << std::endl;
    std::cout << "Expected time domain power (n_active/N): " << (double)n_sc / fft_size << std::endl;
    std::cout << "Expected power ratio (linear): " << (double)fft_size / n_sc << " = " 
              << 10*log10((double)fft_size / n_sc) << " dB" << std::endl;
    
    // Create resource grid with random QPSK symbols on all subcarriers
    ResourceGrid grid(n_ant, n_sym, n_sc);
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(0.0, 1.0);
    
    for (int sym = 0; sym < n_sym; sym++) {
        for (int sc = 0; sc < n_sc; sc++) {
            // QPSK symbol: (±1 ± j)/sqrt(2) -> unit power
            int b0 = (rng() % 2) * 2 - 1;
            int b1 = (rng() % 2) * 2 - 1;
            Complex val(b0 * 0.7071067811865475, b1 * 0.7071067811865475);
            grid.set_re(0, sym, sc, val);
        }
    }
    
    // Set DC to zero
    int dc_pos = fft_size / 2;
    int half_sc = n_sc / 2;
    // DC is at position half_sc in our grid (sc=150 for n_sc=300)
    // But our grid doesn't include DC - let me check...
    // In OfdmModulator::modulate(), DC is explicitly set to 0 at freq[dc_pos]
    // So the n_sc=300 subcarriers are mapped around DC with a gap at DC.
    // That means 300 subcarriers + 1 DC null = 301 used subcarriers out of 512 FFT bins.
    
    // Measure grid power
    double grid_pwr = 0.0;
    int n_re = 0;
    for (int sym = 0; sym < n_sym; sym++) {
        for (int sc = 0; sc < n_sc; sc++) {
            grid_pwr += std::norm(grid.get_re(0, sym, sc));
            n_re++;
        }
    }
    grid_pwr /= n_re;
    std::cout << "\nFreq domain avg RE power: " << grid_pwr << std::endl;
    
    // OFDM modulate
    auto ofdm = create_ofdm_modulator();
    ComplexVec tx_signal = ofdm->modulate(grid, 15000);
    
    double time_pwr = arma::mean(arma::real(tx_signal % arma::conj(tx_signal)));
    std::cout << "Time domain avg power: " << time_pwr << " (expected ~" << (double)n_sc/fft_size << ")" << std::endl;
    std::cout << "Time signal length: " << tx_signal.n_elem << std::endl;
    
    // OFDM demodulate
    ResourceGrid demod = ofdm->demodulate(tx_signal, n_ant, 15000, n_sym);
    std::cout << "Demod grid n_subcarriers: " << demod.n_subcarriers << std::endl;
    
    // Check reconstruction error
    double max_err = 0;
    double mse = 0;
    int n_compared = 0;
    for (int sym = 0; sym < n_sym; sym++) {
        for (int sc = 0; sc < n_sc; sc++) {
            if (sc < (int)demod.n_subcarriers) {
                Complex orig = grid.get_re(0, sym, sc);
                Complex rec = demod.get_re(0, sym, sc);
                double e = std::norm(orig - rec);
                max_err = std::max(max_err, std::abs(orig - rec));
                mse += e;
                n_compared++;
            }
        }
    }
    mse /= n_compared;
    std::cout << "\nReconstruction MSE: " << mse << std::endl;
    std::cout << "Max |error|: " << max_err << std::endl;
    
    // Add noise in time domain and measure noise power after demod
    double sinr_db = 10.0;
    double sinr_lin = std::pow(10.0, sinr_db/10.0);
    double noise_var_slow = time_pwr / sinr_lin;
    
    ComplexVec noisy(tx_signal.n_elem);
    double sigma_dim = std::sqrt(noise_var_slow / 2.0);
    for (arma::uword i = 0; i < tx_signal.n_elem; i++) {
        double nr = dist(rng) * sigma_dim;
        double ni = dist(rng) * sigma_dim;
        noisy(i) = tx_signal(i) + Complex(nr, ni);
    }
    
    ResourceGrid noisy_demod = ofdm->demodulate(noisy, n_ant, 15000, n_sym);
    
    // Measure noise power: for each active RE, noise = rx - tx
    double measured_noise_pwr = 0.0;
    int n_noise_re = 0;
    for (int sym = 0; sym < n_sym; sym++) {
        for (int sc = 0; sc < n_sc; sc++) {
            if (sc < (int)noisy_demod.n_subcarriers) {
                Complex orig = grid.get_re(0, sym, sc);
                Complex rec = noisy_demod.get_re(0, sym, sc);
                measured_noise_pwr += std::norm(rec - orig);
                n_noise_re++;
            }
        }
    }
    measured_noise_pwr /= n_noise_re;
    
    std::cout << "\n=== Noise analysis at SNR=" << sinr_db << " dB ===" << std::endl;
    std::cout << "noise_var_slow (time domain): " << noise_var_slow << std::endl;
    std::cout << "Measured freq-domain noise per active RE: " << measured_noise_pwr << std::endl;
    std::cout << "Expected: noise_var_slow = " << noise_var_slow << std::endl;
    std::cout << "noise_var used in receive() = 1/SNR = " << 1.0/sinr_lin << std::endl;
    std::cout << "Ratio (measured/expected_in_receive): " << measured_noise_pwr / (1.0/sinr_lin) 
              << " = " << 10*log10(measured_noise_pwr / (1.0/sinr_lin)) << " dB" << std::endl;
    
    return 0;
}
