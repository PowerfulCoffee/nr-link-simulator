#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include "phy/ModuleFactory.h"
#include <iostream>
#include <iomanip>

using namespace nr;
using namespace nr::phy;
using namespace nr::channel;

int main() {
    SimulationConfig config;
    config.mcs_index = 5;
    config.n_rb = 25;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.max_blocks_per_sinr = 10;
    config.target_block_errors = 5;
    config.sinr_start = 10.0;
    config.sinr_end = 10.0;
    config.sinr_step = 2.0;
    config.channel_type = ChannelType::AWGN;
    config.mod_scheme = mcs_to_modulation(config.mcs_index);
    config.code_rate = mcs_to_code_rate(config.mcs_index);

    std::cout << "SISO Test: 1TX, 1RX, 1 layer, " << config.n_rb << " RBs, MCS " << config.mcs_index << "\n";
    std::cout << "Modulation: ";
    switch (config.mod_scheme) {
        case ModulationScheme::QPSK: std::cout << "QPSK"; break;
        case ModulationScheme::QAM16: std::cout << "16QAM"; break;
        case ModulationScheme::QAM64: std::cout << "64QAM"; break;
        case ModulationScheme::QAM256: std::cout << "256QAM"; break;
    }
    std::cout << ", Code rate: " << config.code_rate << "\n\n";

    std::vector<BlerResult> results = run_bler_simulation(config, nullptr, "LS");

    std::cout << "Results:\n";
    std::cout << std::setw(10) << "SINR(dB)" << std::setw(12) << "Blocks" << std::setw(12) << "Errors" << std::setw(12) << "BLER" << "\n";
    for (const auto& res : results) {
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(10) << res.sinr_db
                  << std::setw(12) << res.n_blocks
                  << std::setw(12) << res.n_errors
                  << std::setw(12) << std::setprecision(4) << res.bler
                  << "\n";
    }

    bool any_ok = false;
    for (const auto& res : results) {
        if (res.n_blocks > 0) {
            any_ok = true;
            if (res.bler < 1.0) {
                std::cout << "\nSISO link works! BLER=" << res.bler << " at SINR=" << res.sinr_db << " dB\n";
            }
        }
    }

    if (!any_ok) {
        std::cout << "\nNo blocks processed!\n";
        return 1;
    }

    std::cout << "Test completed.\n";
    return 0;
}
