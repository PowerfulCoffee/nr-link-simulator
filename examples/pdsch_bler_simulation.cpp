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

int main(int argc, char* argv[]) {
    int mcs = 5;
    double sinr_start = -2.0;
    double sinr_end = 8.0;
    double sinr_step = 1.0;
    int max_blocks = 200;
    int target_errors = 30;
    bool perfect_csi = false;
    
    if (argc > 1) mcs = atoi(argv[1]);
    if (argc > 2) sinr_start = atof(argv[2]);
    if (argc > 3) sinr_end = atof(argv[3]);
    if (argc > 4) sinr_step = atof(argv[4]);
    if (argc > 5) perfect_csi = atoi(argv[5]) != 0;
    if (argc > 6) max_blocks = atoi(argv[6]);
    if (argc > 7) target_errors = atoi(argv[7]);
    
    SimulationConfig config;
    config.mcs_index = mcs;
    config.n_rb = 25;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.max_blocks_per_sinr = max_blocks;
    config.target_block_errors = target_errors;
    config.sinr_start = sinr_start;
    config.sinr_end = sinr_end;
    config.sinr_step = sinr_step;
    config.n_sinr_points = static_cast<int>((sinr_end - sinr_start) / sinr_step) + 1;
    config.channel_type = ChannelType::AWGN;
    config.mod_scheme = mcs_to_modulation(config.mcs_index);
    config.code_rate = mcs_to_code_rate(config.mcs_index);
    config.dmrs_type = DmrsType::TYPE1;
    config.dmrs_additional_pos = 0;
    config.dmrs_duration = 1;
    config.perfect_csi = perfect_csi;
    
    std::cout << "========================================\n";
    std::cout << "  NR PDSCH BLER Simulation (" << (perfect_csi ? "Perfect CE" : "LS-CE") << ")\n";
    std::cout << "========================================\n";
    std::cout << "Configuration:\n";
    std::cout << "  MCS Index:        " << config.mcs_index << "\n";
    std::cout << "  Modulation:       " << (int)mcs_to_bits_per_symbol(config.mcs_index) << " QAM\n";
    std::cout << "  Code rate:        " << config.code_rate << "\n";
    std::cout << "  Bandwidth:        " << config.n_rb << " PRBs\n";
    std::cout << "  SINR Range:       [" << config.sinr_start << ", " << config.sinr_end << "] dB, step " << config.sinr_step << " dB\n";
    std::cout << "  Max blocks/SINR:  " << config.max_blocks_per_sinr << "\n";
    std::cout << "  Target errors:    " << config.target_block_errors << "\n";
    std::cout << "  DMRS:             Type1, single symbol, pos=2\n";
    std::cout << "  Channel Est:      " << (perfect_csi ? "Ideal/Perfect" : "LS + noise est") << "\n";
    std::cout << "========================================\n\n";
    
    std::vector<BlerResult> results = run_bler_simulation(config, nullptr, perfect_csi ? "Perfect" : "LS");
    
    std::cout << "\n========================================\n";
    std::cout << "  Results Summary\n";
    std::cout << "========================================\n";
    std::cout << std::setw(10) << "SINR(dB)" << std::setw(12) << "Blocks" << std::setw(12) << "Errors" << std::setw(12) << "BLER" << "\n";
    std::cout << "----------------------------------------\n";
    
    for (size_t i = 0; i < results.size(); i++) {
        const BlerResult& res = results[i];
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(10) << res.sinr_db
                  << std::setw(12) << res.n_blocks
                  << std::setw(12) << res.n_errors
                  << std::setw(12) << std::setprecision(6) << res.bler
                  << "\n";
    }
    
    return 0;
}
