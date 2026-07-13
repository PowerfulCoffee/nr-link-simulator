#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace nr;
using namespace nr::phy;

int main() {
    SimulationConfig config;
    config.mcs_index = 5;
    config.n_rb = 6;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.max_blocks_per_sinr = 200;
    config.target_block_errors = 30;
    config.sinr_start = -3.0;
    config.sinr_end = 5.0;
    config.sinr_step = 1.0;
    config.channel_type = ChannelType::AWGN;
    config.mod_scheme = mcs_to_modulation(config.mcs_index);
    config.code_rate = mcs_to_code_rate(config.mcs_index);
    config.n_ldpc_iterations = 30;
    config.early_termination = true;
    
    const char* mod_name = "QPSK";
    if (config.mod_scheme == ModulationScheme::QAM16) mod_name = "16QAM";
    else if (config.mod_scheme == ModulationScheme::QAM64) mod_name = "64QAM";
    
    std::cout << "MCS " << config.mcs_index << ": " << mod_name << ", R=" << config.code_rate << "\n";
    
    auto results = run_bler_simulation(config, nullptr, "LS");
    
    std::ofstream out("bler_qpsk.csv");
    out << "SINR_dB,Blocks,Errors,BLER\n";
    std::cout << std::setw(10) << "SINR(dB)" << std::setw(10) << "BLER\n";
    for (auto& r : results) {
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(10) << r.sinr_db
                  << std::setw(12) << std::setprecision(4) << r.bler << "\n";
        out << std::setprecision(2) << r.sinr_db << "," << r.n_blocks << "," << r.n_errors << "," << std::setprecision(6) << r.bler << "\n";
    }
    out.close();
    return 0;
}
