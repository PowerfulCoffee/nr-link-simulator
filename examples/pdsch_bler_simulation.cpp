#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include "phy/ModuleFactory.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <memory>

using namespace nr;
using namespace nr::phy;
using namespace nr::channel;

int main() {
    SimulationConfig config;
    config.mcs_index = 10;
    config.n_rb = 25;
    config.n_tx_ant = 2;
    config.n_rx_ant = 2;
    config.n_layers = 2;
    config.max_blocks_per_sinr = 20;
    config.target_block_errors = 5;
    config.sinr_start = 0.0;
    config.sinr_end = 10.0;
    config.sinr_step = 2.0;
    config.channel_type = ChannelType::AWGN;
    config.mod_scheme = mcs_to_modulation(config.mcs_index);
    config.code_rate = mcs_to_code_rate(config.mcs_index);

    std::cout << "========================================\n";
    std::cout << "  NR PDSCH BLER Simulation (Quick Test)\n";
    std::cout << "========================================\n";
    std::cout << "Configuration:\n";
    std::cout << "  MCS Index:        " << config.mcs_index << "\n";
    std::cout << "  Bandwidth:        " << config.n_rb << " PRBs\n";
    std::cout << "  TX Antennas:      " << config.n_tx_ant << "\n";
    std::cout << "  RX Antennas:      " << config.n_rx_ant << "\n";
    std::cout << "  Layers:           " << config.n_layers << "\n";
    std::cout << "  SINR Range:       [" << config.sinr_start << ", " << config.sinr_end << "] dB, step " << config.sinr_step << " dB\n";
    std::cout << "  Max blocks/SINR:  " << config.max_blocks_per_sinr << "\n";
    std::cout << "  Target errors:    " << config.target_block_errors << "\n";
    std::cout << "========================================\n\n";

    std::cout << "Starting BLER simulation...\n\n";

    std::vector<BlerResult> results = run_bler_simulation(config, nullptr, "LS");

    std::cout << "\n========================================\n";
    std::cout << "  Results Summary\n";
    std::cout << "========================================\n";
    std::cout << std::setw(10) << "SINR(dB)" << std::setw(12) << "Blocks" << std::setw(12) << "Errors" << std::setw(12) << "BLER" << "\n";
    std::cout << "----------------------------------------\n";

    std::ofstream outfile("bler_results.csv");
    outfile << "SINR_dB,Blocks,Errors,BLER\n";

    for (size_t i = 0; i < results.size(); i++) {
        const BlerResult& res = results[i];
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(10) << res.sinr_db
                  << std::setw(12) << res.n_blocks
                  << std::setw(12) << res.n_errors
                  << std::setw(12) << std::setprecision(6) << res.bler
                  << "\n";
        outfile << std::fixed << std::setprecision(2)
                << res.sinr_db << ","
                << res.n_blocks << ","
                << res.n_errors << ","
                << std::setprecision(6) << res.bler << "\n";
    }

    outfile.close();
    std::cout << "========================================\n";
    std::cout << "Results saved to: bler_results.csv\n";

    return 0;
}
