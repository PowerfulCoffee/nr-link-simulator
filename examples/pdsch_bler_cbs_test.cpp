#include "phy/PdschProcessor.h"
#include "common/NrTables.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace nr;
using namespace nr::phy;

int main() {
    std::cout << "=== NR PDSCH Code Block Segmentation Test ===" << std::endl;

    SimulationConfig config;
    config.mcs_index = 27;
    config.mod_scheme = ModulationScheme::QAM64;
    config.code_rate = 910.0 / 1024.0;
    config.n_rb = 25;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.scs = 15e3;
    config.channel_type = ChannelType::AWGN;
    config.perfect_csi = false;
    config.random_seed = 42;
    config.n_ldpc_iterations = 20;
    config.early_termination = true;
    config.max_blocks_per_sinr = 100;
    config.target_block_errors = 50;
    config.sinr_start = 30.0;
    config.sinr_end = 40.0;
    config.sinr_step = 2.0;
    config.n_sinr_points = 6;

    int qm = mcs_to_bits_per_symbol(config.mcs_index);
    double target_coderate = mcs_to_code_rate(config.mcs_index);
    int n_re_per_prb = 13 * 12;

    int tbs = calculate_tbs(config.n_rb, n_re_per_prb, qm, config.n_layers, target_coderate);
    std::cout << "Configuration: " << config.n_rb << " PRBs, MCS" << config.mcs_index
              << ", Qm=" << qm << ", R=" << target_coderate << std::endl;
    std::cout << "Calculated TBS: " << tbs << " bits" << std::endl;

    int tb_crc_len = get_crc_length(tbs);
    int b_after_tb_crc = tbs + tb_crc_len;
    std::cout << "After TB CRC (" << tb_crc_len << " bits): " << b_after_tb_crc << " bits" << std::endl;

    int G = config.n_rb * n_re_per_prb * qm * config.n_layers;
    CodeBlockSegParams cbs = compute_cb_segmentation(tbs, b_after_tb_crc,
                                                      n_re_per_prb, config.n_rb,
                                                      qm, config.n_layers, target_coderate);
    std::cout << "Code Block Segmentation:" << std::endl;
    std::cout << "  num_cb = " << cbs.num_cb << std::endl;
    std::cout << "  tb_crc_len = " << cbs.tb_crc_len << std::endl;
    std::cout << "  cb_crc_len = " << cbs.cb_crc_len << std::endl;
    std::cout << "  cb_info_bits = " << cbs.cb_info_bits << std::endl;
    std::cout << "  cb_size_with_crc = " << cbs.cb_size_with_crc << std::endl;
    std::cout << "  bgn = " << cbs.bgn << std::endl;
    std::cout << "  zc = " << cbs.zc << std::endl;
    std::cout << "  k_b = " << cbs.k_b << std::endl;
    std::cout << "  cb_k = " << cbs.cb_k << std::endl;
    std::cout << "  total cw_length = " << cbs.cw_length << " (G=" << G << ")" << std::endl;

    std::cout << "\nRunning BLER test (high SNR to verify zero-error decoding):" << std::endl;
    std::cout << std::setw(12) << "SINR(dB)" << std::setw(10) << "Blocks"
              << std::setw(10) << "Errors" << std::setw(12) << "BLER" << std::endl;
    std::cout << "----------------------------------------------" << std::endl;

    PdschProcessor processor(config);

    std::ofstream csv("bler_cbs_test.csv");
    csv << "sinr_db,n_blocks,n_errors,bler,num_cb" << std::endl;

    for (int i = 0; i < config.n_sinr_points; i++) {
        double sinr_db = config.sinr_start + i * config.sinr_step;
        BlerResult res;
        processor.process_single_snr_point(sinr_db, res);

        std::cout << std::setw(12) << std::fixed << std::setprecision(2) << res.sinr_db
                  << std::setw(10) << res.n_blocks
                  << std::setw(10) << res.n_errors
                  << std::setw(12) << std::setprecision(4) << res.bler << std::endl;
        csv << res.sinr_db << "," << res.n_blocks << "," << res.n_errors << "," << res.bler << "," << cbs.num_cb << std::endl;

        if (res.n_errors == 0 && res.n_blocks >= 50) {
            std::cout << "Zero BLER achieved at " << sinr_db << " dB - CBS works!" << std::endl;
            break;
        }
    }

    csv.close();
    std::cout << "\nResults saved to bler_cbs_test.csv" << std::endl;

    return 0;
}
