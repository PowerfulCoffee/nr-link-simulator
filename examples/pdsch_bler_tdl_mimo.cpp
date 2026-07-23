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
    int mcs = 10;
    double sinr_start = 0.0;
    double sinr_end = 20.0;
    double sinr_step = 2.0;
    int max_blocks = 100;
    int target_errors = 30;
    bool perfect_csi = false;
    int n_rb = 25;
    int n_rx_ant = 4;
    int n_layers = 4;
    int channel_type_int = 1;
    double delay_spread = 100e-9;
    double max_doppler = 70.0;
    int dmrs_add_pos = 1;
    std::string out_file = "";
    
    if (argc > 1) mcs = atoi(argv[1]);
    if (argc > 2) sinr_start = atof(argv[2]);
    if (argc > 3) sinr_end = atof(argv[3]);
    if (argc > 4) sinr_step = atof(argv[4]);
    if (argc > 5) perfect_csi = atoi(argv[5]) != 0;
    if (argc > 6) max_blocks = atoi(argv[6]);
    if (argc > 7) target_errors = atoi(argv[7]);
    if (argc > 8) n_rb = atoi(argv[8]);
    if (argc > 9) n_rx_ant = atoi(argv[9]);
    if (argc > 10) n_layers = atoi(argv[10]);
    if (argc > 11) channel_type_int = atoi(argv[11]);
    if (argc > 12) delay_spread = atof(argv[12]);
    if (argc > 13) max_doppler = atof(argv[13]);
    if (argc > 14) dmrs_add_pos = atoi(argv[14]);
    if (argc > 15) out_file = argv[15];
    
    ChannelType ch_type;
    std::string ch_name;
    switch (channel_type_int) {
        case 0: ch_type = ChannelType::AWGN; ch_name = "AWGN"; break;
        case 1: ch_type = ChannelType::TDL_A; ch_name = "TDL-A"; break;
        case 2: ch_type = ChannelType::TDL_B; ch_name = "TDL-B"; break;
        case 3: ch_type = ChannelType::TDL_C; ch_name = "TDL-C"; break;
        case 4: ch_type = ChannelType::TDL_D; ch_name = "TDL-D"; break;
        case 5: ch_type = ChannelType::TDL_E; ch_name = "TDL-E"; break;
        default: ch_type = ChannelType::TDL_A; ch_name = "TDL-A"; break;
    }
    
    int n_tx_ant = n_layers;
    
    SimulationConfig config;
    config.mcs_index = mcs;
    config.n_rb = n_rb;
    config.n_tx_ant = n_tx_ant;
    config.n_rx_ant = n_rx_ant;
    config.n_layers = n_layers;
    config.max_blocks_per_sinr = max_blocks;
    config.target_block_errors = target_errors;
    config.sinr_start = sinr_start;
    config.sinr_end = sinr_end;
    config.sinr_step = sinr_step;
    config.n_sinr_points = static_cast<int>((sinr_end - sinr_start) / sinr_step) + 1;
    config.channel_type = ch_type;
    config.delay_spread = delay_spread;
    config.max_doppler = max_doppler;
    config.enable_los = (ch_type == ChannelType::TDL_D || ch_type == ChannelType::TDL_E);
    config.mod_scheme = mcs_to_modulation(config.mcs_index);
    config.code_rate = mcs_to_code_rate(config.mcs_index);
    config.dmrs_type = DmrsType::TYPE1;
    config.dmrs_additional_pos = dmrs_add_pos;
    config.dmrs_duration = 1;
    config.perfect_csi = perfect_csi;
    
    std::cout << "========================================\n";
    std::cout << "  NR PDSCH BLER TDL MIMO Simulation\n";
    std::cout << "========================================\n";
    std::cout << "Configuration:\n";
    std::cout << "  MCS Index:        " << config.mcs_index << " (" << (int)mcs_to_bits_per_symbol(config.mcs_index) << "QAM, R=" << config.code_rate << ")\n";
    std::cout << "  Bandwidth:        " << config.n_rb << " PRBs\n";
    std::cout << "  MIMO:             " << config.n_tx_ant << "Tx x " << config.n_rx_ant << "Rx, " << config.n_layers << " layers\n";
    std::cout << "  Channel:          " << ch_name << ", DS=" << delay_spread*1e9 << "ns, fd=" << max_doppler << "Hz\n";
    std::cout << "  SINR Range:       [" << config.sinr_start << ", " << config.sinr_end << "] dB, step " << config.sinr_step << " dB\n";
    std::cout << "  Max blocks/SINR:  " << config.max_blocks_per_sinr << "\n";
    std::cout << "  Target errors:    " << config.target_block_errors << "\n";
    std::cout << "  DMRS:             Type1, add_pos=" << dmrs_add_pos << "\n";
    std::cout << "  Channel Est:      " << (perfect_csi ? "Ideal/Perfect" : "LS + Doppler est/comp") << "\n";
    std::cout << "========================================\n\n";
    
    if (!out_file.empty()) {
        std::ofstream f(out_file);
        f << "SINR_dB,Blocks,Errors,BLER\n";
        f.close();
    }
    
    std::vector<BlerResult> results = run_bler_simulation(config, nullptr, perfect_csi ? "Perfect" : "LS-Doppler");
    
    std::cout << "\n========================================\n";
    std::cout << "  Results Summary\n";
    std::cout << "========================================\n";
    std::cout << std::setw(10) << "SINR(dB)" << std::setw(12) << "Blocks" << std::setw(12) << "Errors" << std::setw(12) << "BLER" << "\n";
    std::cout << "----------------------------------------\n";
    
    if (!out_file.empty()) {
        std::ofstream f(out_file, std::ios::app);
        for (size_t i = 0; i < results.size(); i++) {
            const BlerResult& res = results[i];
            f << res.sinr_db << "," << res.n_blocks << "," << res.n_errors << "," << res.bler << "\n";
        }
        f.close();
    }
    
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
