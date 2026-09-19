#include "common/Types.h"
#include "phy/PdschProcessor.h"
#include "phy/ModuleFactory.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <memory>

using namespace nr;
using namespace nr::phy;
using namespace nr::channel;

int main() {
    SimulationConfig config;
    config.mcs_index = 4;
    config.mcs_table = 2;
    config.n_rb = 30;
    config.n_tx_ant = 4;
    config.n_rx_ant = 4;
    config.n_layers = 4;
    config.max_blocks_per_sinr = 100;
    config.target_block_errors = 30;
    config.sinr_start = -4.0;
    config.sinr_end = 12.0;
    config.sinr_step = 2.0;
    config.n_sinr_points = 9;
    config.channel_type = ChannelType::AWGN;  // Use AWGN but we'll generate flat Rayleigh manually? No, use single-tap via TDL?
    config.delay_spread = 0.0;  // DS=0 means single tap flat fading!
    config.max_doppler = 5.0;
    config.enable_los = false;
    config.dmrs_type = DmrsType::TYPE1;
    config.dmrs_additional_pos = 1;
    config.dmrs_duration = 0;
    config.perfect_csi = true;
    config.n_ldpc_iterations = 20;
    config.mod_scheme = mcs_to_modulation(config.mcs_index, config.mcs_table);
    config.code_rate = mcs_to_code_rate(config.mcs_index, config.mcs_table);

    std::cout << "4x4 Flat Rayleigh (DS=0) QPSK MCS4\n";
    auto channel = ModuleFactory::create_channel(ChannelType::TDL_A);
    channel->set_config(config);

    std::vector<BlerResult> results = run_bler_simulation(config, channel.get(), "Perfect");

    std::cout << "\nResults:\n";
    std::cout << std::setw(10) << "SINR(dB)" << std::setw(12) << "Blocks" << std::setw(12) << "Errors" << std::setw(12) << "BLER" << "\n";
    for (auto& r : results) {
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(10) << r.sinr_db
                  << std::setw(12) << r.n_blocks
                  << std::setw(12) << r.n_errors
                  << std::setw(12) << std::setprecision(6) << r.bler
                  << "\n";
    }
    return 0;
}
