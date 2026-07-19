#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include <iostream>
#include <iomanip>

using namespace nr;
using namespace nr::phy;

int main() {
    SimulationConfig config;
    config.mcs_index = 27;
    config.n_rb = 3;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.max_blocks_per_sinr = 5;
    config.target_block_errors = 100;
    config.sinr_start = 25.0;
    config.sinr_end = 25.0;
    config.sinr_step = 1.0;
    config.channel_type = ChannelType::AWGN;
    config.mod_scheme = mcs_to_modulation(config.mcs_index);
    config.code_rate = mcs_to_code_rate(config.mcs_index);
    config.n_ldpc_iterations = 20;
    config.early_termination = true;
    config.random_seed = 42;
    config.scs = 15000;

    int qm = mcs_to_bits_per_symbol(config.mcs_index);
    double R = mcs_to_code_rate(config.mcs_index);
    int n_re_per_prb = 13*12;
    int tbs = calculate_tbs(config.n_rb, n_re_per_prb, qm, config.n_layers, R);
    std::cout << "TBS = " << tbs << " bits, Qm=" << qm << ", R=" << R << "\n";

    auto results = run_bler_simulation(config, nullptr, "LS");
    for (auto& r : results) {
        std::cout << "SINR=" << r.sinr_db << " dB: " << r.n_errors << "/" << r.n_blocks
                  << " errors, BLER=" << std::fixed << std::setprecision(4) << r.bler << "\n";
    }
    return 0;
}
