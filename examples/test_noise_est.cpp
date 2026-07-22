#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include "phy/ModuleFactory.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>

using namespace nr;
using namespace nr::phy;
using namespace nr::channel;

int main() {
    std::cout << "=== Noise Estimation Accuracy Test (AWGN channel, LS CE) ===\n\n";

    SimulationConfig config;
    config.mcs_index = 7;
    config.n_rb = 25;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.channel_type = ChannelType::AWGN;
    config.perfect_csi = false;
    config.n_ldpc_iterations = 25;
    config.random_seed = 42;
    config.scs = 15000;
    config.mod_scheme = mcs_to_modulation(config.mcs_index);
    config.code_rate = mcs_to_code_rate(config.mcs_index);
    config.dmrs_type = DmrsType::TYPE1;
    config.dmrs_additional_pos = 0;
    config.dmrs_duration = 1;

    PdschProcessor proc(config);
    int n_test_blocks = 100;

    for (int snr_db = 0; snr_db <= 25; snr_db += 5) {
        double ideal_noise_var = std::pow(10.0, -snr_db / 10.0);
        double est_sum = 0.0;
        double est_min = 1e9;
        double est_max = 0.0;
        int n_ok = 0;

        proc.set_seed(12345 + snr_db);
        for (int b = 0; b < n_test_blocks; b++) {
            TransportBlock tb = proc.generate_transport_block();
            PdschTxResult tx_res = proc.transmit(tb, b);

            ResourceGrid rx_grid = tx_res.tx_grid;
            double sigma = std::sqrt(ideal_noise_var / 2.0);
            std::mt19937 rng(42 + b * 1000 + snr_db);
            std::normal_distribution<double> nd(0.0, sigma);
            for (int ant = 0; ant < rx_grid.n_ant; ant++) {
                for (int sym = 0; sym < rx_grid.n_symbols; sym++) {
                    for (int sc = 0; sc < rx_grid.n_subcarriers; sc++) {
                        Complex val = rx_grid.get_re(ant, sym, sc);
                        double nr = nd(rng);
                        double ni = nd(rng);
                        rx_grid.set_re(ant, sym, sc, val + Complex(nr, ni));
                    }
                }
            }

            PdschRxResult rx_res = proc.receive(rx_grid, tx_res, snr_db, b);
            double est_noise_var = rx_res.noise_var_est;

            if (est_noise_var > 1e-12) {
                est_sum += est_noise_var;
                est_min = std::min(est_min, est_noise_var);
                est_max = std::max(est_max, est_noise_var);
                n_ok++;
            }
        }

        double avg_est = est_sum / n_ok;
        double est_snr_db = -10.0 * std::log10(avg_est);
        double error_pct = (avg_est - ideal_noise_var) / ideal_noise_var * 100.0;

        std::cout << "SNR = " << std::setw(2) << snr_db << " dB: "
                  << "ideal_var = " << std::scientific << std::setprecision(4) << ideal_noise_var
                  << ", avg_est = " << avg_est
                  << " (est_SNR = " << std::fixed << std::setprecision(2) << est_snr_db << " dB)"
                  << ", bias = " << std::showpos << std::setprecision(2) << error_pct << "%"
                  << std::noshowpos
                  << ", range [" << std::scientific << est_min << ", " << est_max << "]\n";
    }

    std::cout << "\n=== Test complete ===\n";
    std::cout << "Notes:\n";
    std::cout << "  - Noise is estimated from second-order differences of adjacent DMRS LS estimates\n";
    std::cout << "  - Internal DMRS subcarriers only (edge points excluded to avoid channel slope bias)\n";
    std::cout << "  - Variance scaling factor 2/3 applied for unbiased estimate\n";
    return 0;
}
