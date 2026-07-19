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
    SimulationConfig config;
    config.mcs_index = 27;
    config.n_rb = 3;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.max_blocks_per_sinr = 800;
    config.target_block_errors = 60;
    config.sinr_start = 17.0;
    config.sinr_end = 21.0;
    config.sinr_step = 0.5;
    config.n_sinr_points = static_cast<int>((config.sinr_end - config.sinr_start) / config.sinr_step) + 1;
    config.channel_type = ChannelType::AWGN;
    config.mod_scheme = mcs_to_modulation(config.mcs_index);
    config.code_rate = mcs_to_code_rate(config.mcs_index);
    config.n_ldpc_iterations = 20;
    config.early_termination = true;
    config.random_seed = 123;
    config.scs = 15000;

    int qm_int = mcs_to_bits_per_symbol(config.mcs_index);
    double R = mcs_to_code_rate(config.mcs_index);
    const char* mod_name = "64QAM";

    int n_re_per_prb = 13 * 12;
    int tbs = calculate_tbs(config.n_rb, n_re_per_prb, qm_int, config.n_layers, R);

    std::cout << "=== NR PDSCH BLER Simulation (C++ nr-link-simulator) ===\n";
    std::cout << "MCS " << config.mcs_index << ": " << mod_name
              << ", Qm=" << qm_int << ", R=" << R << " (" << (int)(R*1024) << "/1024)\n";
    std::cout << "PRBs=" << config.n_rb << ", TBS=" << tbs << " bits, SISO AWGN\n";
    std::cout << "SINR: [" << config.sinr_start << ", " << config.sinr_end
              << "] dB, step " << config.sinr_step << " dB\n";
    std::cout << "Max blocks/SINR=" << config.max_blocks_per_sinr
              << ", target errors=" << config.target_block_errors << "\n\n";

    auto results = run_bler_simulation(config, nullptr, "LS");

    std::string out_file = "bler_mcs27_awgn_fine.csv";
    std::ofstream out(out_file);
    out << "SINR_dB,Blocks,Errors,BLER\n";
    std::cout << std::setw(10) << "SINR(dB)" << std::setw(10) << "Blocks"
              << std::setw(10) << "Errors" << std::setw(12) << "BLER\n";
    std::cout << std::string(46, '-') << "\n";
    for (auto& res : results) {
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(10) << res.sinr_db
                  << std::setw(10) << res.n_blocks
                  << std::setw(10) << res.n_errors
                  << std::setw(12) << std::setprecision(4) << res.bler << "\n";
        out << std::setprecision(2) << res.sinr_db << "," << res.n_blocks << ","
            << res.n_errors << "," << std::setprecision(6) << res.bler << "\n";
    }
    out.close();
    std::cout << "\nResults saved to " << out_file << "\n";
    return 0;
}
