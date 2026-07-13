#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace nr;
using namespace nr::phy;

int main() {
    SimulationConfig config;
    config.mcs_index = 5;
    config.n_rb = 6;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.max_blocks_per_sinr = 500;
    config.target_block_errors = 50;
    config.sinr_start = -4.0;
    config.sinr_end = 4.0;
    config.sinr_step = 0.5;
    config.channel_type = ChannelType::AWGN;
    config.mod_scheme = mcs_to_modulation(config.mcs_index);
    config.code_rate = mcs_to_code_rate(config.mcs_index);
    config.n_ldpc_iterations = 30;
    config.early_termination = true;
    
    std::cout << "=== NR Link Simulator - PDSCH BLER (AWGN, Ideal CSI) ===\n";
    std::cout << "MCS " << config.mcs_index << ": QPSK, R=" << config.code_rate << "\n";
    std::cout << "PRBs = " << config.n_rb << ", 1x1 SISO, LDPC iter=" << config.n_ldpc_iterations << "\n\n";
    
    std::ofstream out("../bler_awgn_ideal.csv");
    out << "SINR_dB,Blocks,Errors,BLER\n";
    std::cout << std::setw(10) << "SINR(dB)" << std::setw(10) << "Blocks" << std::setw(10) << "Errors" << std::setw(12) << "BLER\n";
    
    for (double sinr = config.sinr_start; sinr <= config.sinr_end + 0.01; sinr += config.sinr_step) {
        SimulationConfig cfg = config;
        PdschProcessor proc(cfg);
        auto channel = channel::create_channel(cfg.channel_type);
        proc.set_channel_model(std::move(channel));
        BlerResult res;
        proc.process_single_snr_point(sinr, res);
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(10) << res.sinr_db
                  << std::setw(10) << res.n_blocks
                  << std::setw(10) << res.n_errors
                  << std::setw(12) << std::setprecision(4) << res.bler << "\n";
        out << std::setprecision(2) << res.sinr_db << "," << res.n_blocks << "," << res.n_errors << "," << std::setprecision(6) << res.bler << "\n";
        out.flush();
    }
    out.close();
    std::cout << "\nResults saved to bler_awgn_ideal.csv\n";
    return 0;
}
