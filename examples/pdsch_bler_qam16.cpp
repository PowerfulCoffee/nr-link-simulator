#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace nr;
using namespace nr::phy;

int main(int argc, char** argv) {
    int mcs = 10;
    int n_rb = 25;
    double sinr_start = 0.0;
    double sinr_end = 12.0;
    double sinr_step = 0.5;
    int max_blocks = 200;
    int target_errs = 50;
    
    if (argc > 1) mcs = atoi(argv[1]);
    if (argc > 2) n_rb = atoi(argv[2]);
    
    SimulationConfig config;
    config.mcs_index = mcs;
    config.n_rb = n_rb;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.max_blocks_per_sinr = max_blocks;
    config.target_block_errors = target_errs;
    config.sinr_start = sinr_start;
    config.sinr_end = sinr_end;
    config.sinr_step = sinr_step;
    config.channel_type = ChannelType::AWGN;
    config.mod_scheme = mcs_to_modulation(config.mcs_index);
    config.code_rate = mcs_to_code_rate(config.mcs_index);
    config.n_ldpc_iterations = 20;
    
    const char* mod_name = "BPSK";
    if (config.mod_scheme == ModulationScheme::QPSK) mod_name = "QPSK";
    else if (config.mod_scheme == ModulationScheme::QAM16) mod_name = "16QAM";
    else if (config.mod_scheme == ModulationScheme::QAM64) mod_name = "64QAM";
    
    std::cout << "========================================\n";
    std::cout << "  NR PDSCH BLER Simulation\n";
    std::cout << "========================================\n";
    std::cout << "  MCS = " << mcs << " (" << mod_name << ", R=" << config.code_rate << ")\n";
    std::cout << "  PRBs = " << n_rb << "\n";
    std::cout << "  SISO AWGN\n";
    std::cout << "========================================\n\n";
    
    auto results = run_bler_simulation(config, nullptr, "LS");
    
    std::ofstream out("bler_results.csv");
    out << "SINR_dB,Blocks,Errors,BLER\n";
    std::cout << std::setw(10) << "SINR(dB)" << std::setw(10) << "Blocks" << std::setw(10) << "Errors" << std::setw(12) << "BLER\n";
    for (auto& r : results) {
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(10) << r.sinr_db
                  << std::setw(10) << r.n_blocks
                  << std::setw(10) << r.n_errors
                  << std::setw(12) << std::setprecision(4) << r.bler << "\n";
        out << std::setprecision(2) << r.sinr_db << "," << r.n_blocks << "," << r.n_errors << "," << std::setprecision(6) << r.bler << "\n";
    }
    out.close();
    std::cout << "\nSaved to bler_results.csv\n";
    return 0;
}
