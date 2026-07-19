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
    config.max_blocks_per_sinr = 200;
    config.target_block_errors = 30;
    config.sinr_start = 20.0;
    config.sinr_end = 32.0;
    config.sinr_step = 2.0;
    config.n_sinr_points = 7;
    config.channel_type = ChannelType::TDL_A;
    config.delay_spread = 30e-9;
    config.max_doppler = 0.0;
    config.perfect_csi = true;
    config.mod_scheme = mcs_to_modulation(config.mcs_index);
    config.code_rate = mcs_to_code_rate(config.mcs_index);
    config.n_ldpc_iterations = 20;
    config.early_termination = true;
    config.random_seed = 42;
    config.scs = 15000;

    std::cout << "Testing TDL-A (perfect CSI) MCS27, 3PRB, 15kHz SCS, DS=30ns\n";
    auto results = run_bler_simulation(config, nullptr, "perfect");
    for (auto& res : results) {
        std::cout << "SNR=" << std::fixed << std::setprecision(1) << res.sinr_db
                  << "dB: BLER=" << res.bler << " (" << res.n_errors << "/" << res.n_blocks << ")\n";
    }
    return 0;
}
