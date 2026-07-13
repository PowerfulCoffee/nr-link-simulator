#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include <iostream>
#include <iomanip>

using namespace nr;
using namespace nr::phy;

int main(int argc, char** argv) {
    int mcs = 9;
    if (argc > 1) mcs = atoi(argv[1]);
    double sinr_test = -5.0;
    if (argc > 2) sinr_test = atof(argv[2]);
    
    SimulationConfig config;
    config.mcs_index = mcs;
    config.n_rb = 25;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.max_blocks_per_sinr = 50;
    config.target_block_errors = 100;
    config.sinr_start = sinr_test;
    config.sinr_end = sinr_test;
    config.sinr_step = 1.0;
    config.channel_type = ChannelType::AWGN;
    config.mod_scheme = mcs_to_modulation(config.mcs_index);
    config.code_rate = mcs_to_code_rate(config.mcs_index);
    config.n_ldpc_iterations = 20;
    
    std::cout << "Testing MCS " << mcs << " at SINR=" << sinr_test << "dB\n";
    std::cout << "Modulation: ";
    if (config.mod_scheme == ModulationScheme::QPSK) std::cout << "QPSK";
    else if (config.mod_scheme == ModulationScheme::QAM16) std::cout << "16QAM";
    else std::cout << "Other";
    std::cout << ", R=" << config.code_rate << "\n";
    
    auto results = run_bler_simulation(config, nullptr, "LS");
    for (auto& r : results) {
        std::cout << "SINR " << r.sinr_db << "dB: " << r.n_errors << "/" << r.n_blocks 
                  << " errors, BLER=" << std::fixed << std::setprecision(4) << r.bler << "\n";
    }
    return 0;
}
