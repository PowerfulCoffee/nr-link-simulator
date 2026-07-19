/*
 * PDSCH BLER vs SNR simulation for ALL MCS indices (0-27) on AWGN channel.
 * Uses run_bler_simulation() for each MCS with per-MCS SNR ranges.
 * Uses PDSCH MCS Table 2 (256QAM, TS 38.214 Table 5.1.3.1-2).
 *
 * Per-MCS SNR range (1dB steps, adjusted for modulation/code rate):
 *   Low-rate QPSK   (MCS 0-1):   -10 to 6 dB
 *   Mid-rate QPSK   (MCS 2-3):   -6 to 10 dB
 *   High-rate QPSK  (MCS 4):     -4 to 12 dB
 *   Low-rate 16QAM  (MCS 5-7):    0 to 16 dB
 *   High-rate 16QAM (MCS 8-10):   2 to 18 dB
 *   Low-rate 64QAM  (MCS 11-14):  4 to 22 dB
 *   Mid-rate 64QAM  (MCS 15-17):  8 to 24 dB
 *   High-rate 64QAM (MCS 18-19): 10 to 26 dB
 *   Low-rate 256QAM (MCS 20-22): 12 to 28 dB
 *   Mid-rate 256QAM (MCS 23-25): 16 to 32 dB
 *   High-rate256QAM (MCS 26-27): 18 to 34 dB
 */

#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <chrono>
#include <filesystem>

using namespace nr;
using namespace nr::phy;
using namespace std::chrono;
namespace fs = std::filesystem;

struct MCSConfig {
    int mcs;
    double snr_start;
    double snr_end;
};

int main() {
    std::vector<MCSConfig> mcs_configs;
    for (int mcs = 0; mcs < 28; mcs++) {
        double snr_start, snr_end;
        int qm = mcs_to_bits_per_symbol(mcs);
        double R = mcs_to_code_rate(mcs);
        if (qm == 2) {
            if (R < 0.20)       { snr_start = -10.0; snr_end = 6.0; }
            else if (R < 0.45)  { snr_start = -6.0;  snr_end = 10.0; }
            else                { snr_start = -4.0;  snr_end = 12.0; }
        } else if (qm == 4) {
            if (R < 0.50)       { snr_start = 0.0;   snr_end = 16.0; }
            else                { snr_start = 2.0;   snr_end = 18.0; }
        } else if (qm == 6) {
            if (R < 0.60)       { snr_start = 4.0;   snr_end = 22.0; }
            else if (R < 0.76)  { snr_start = 8.0;   snr_end = 24.0; }
            else                { snr_start = 10.0;  snr_end = 26.0; }
        } else { // qm == 8 (256QAM)
            if (R < 0.75)       { snr_start = 12.0;  snr_end = 28.0; }
            else if (R < 0.87)  { snr_start = 16.0;  snr_end = 32.0; }
            else                { snr_start = 18.0;  snr_end = 34.0; }
        }
        mcs_configs.push_back({mcs, snr_start, snr_end});
    }

    const int n_prb = 3;
    const int n_re_per_prb = 156;
    const double sinr_step = 1.0;
    const int max_blocks_per_snr = 5000;
    const int target_block_errors = 100;
    const int n_ldpc_iters = 25;

    std::cout << "PDSCH BLER vs SNR - All MCS (0-27), AWGN, nPRB=" << n_prb
              << " (" << n_re_per_prb << " RE/PRB), SISO, Table 2 (256QAM)" << std::endl;
    std::cout << "Step=" << sinr_step << "dB, max_blocks=" << max_blocks_per_snr
              << ", target_errors=" << target_block_errors
              << ", ldpc_iters=" << n_ldpc_iters << std::endl << std::endl;

    fs::create_directories("results");
    std::ofstream summary("results/bler_all_mcs_awgn.csv");
    summary << "mcs,qm,rate,tbs,snr_db,bler,n_blocks,n_errors\n";

    auto total_start = high_resolution_clock::now();
    int total_snr_points = 0;

    for (const auto& mc : mcs_configs) {
        int mcs = mc.mcs;
        int qm = mcs_to_bits_per_symbol(mcs);
        double R = mcs_to_code_rate(mcs);
        int tbs = calculate_tbs(n_prb, n_re_per_prb, qm, 1, R);
        int n_sinr = static_cast<int>((mc.snr_end - mc.snr_start) / sinr_step) + 1;

        std::cout << "=== MCS " << std::setw(2) << mcs
                  << " (Qm=" << qm << ", R=" << std::fixed << std::setprecision(3) << R
                  << ", TBS=" << tbs << ") SNR "
                  << std::setprecision(0) << mc.snr_start << "~" << mc.snr_end << "dB ==="
                  << std::endl;

        SimulationConfig config;
        config.mcs_index = mcs;
        config.n_rb = n_prb;
        config.n_tx_ant = 1;
        config.n_rx_ant = 1;
        config.n_layers = 1;
        config.mod_scheme = mcs_to_modulation(mcs);
        config.code_rate = R;
        config.channel_type = ChannelType::AWGN;
        config.perfect_csi = true;
        config.max_blocks_per_sinr = max_blocks_per_snr;
        config.target_block_errors = target_block_errors;
        config.sinr_start = mc.snr_start;
        config.sinr_end = mc.snr_end;
        config.sinr_step = sinr_step;
        config.n_sinr_points = n_sinr;
        config.n_ldpc_iterations = n_ldpc_iters;
        config.early_termination = true;
        config.random_seed = 123;
        config.scs = 15000;

        auto results = run_bler_simulation(config, nullptr, "perfect");

        for (const auto& res : results) {
            std::cout << "  SNR " << std::setw(5) << std::fixed << std::setprecision(1)
                      << res.sinr_db << " dB: BLER=" << std::setprecision(4) << res.bler
                      << " (blocks=" << res.n_blocks << ", errors=" << res.n_errors << ")"
                      << std::endl;
            summary << mcs << "," << qm << "," << std::fixed << std::setprecision(4) << R
                    << "," << tbs << "," << std::setprecision(1) << res.sinr_db
                    << "," << std::setprecision(6) << res.bler << ","
                    << res.n_blocks << "," << res.n_errors << "\n";
            total_snr_points++;
        }
        summary.flush();
        std::cout << std::endl;
    }

    auto total_end = high_resolution_clock::now();
    double total_sec = duration_cast<seconds>(total_end - total_start).count();
    std::cout << "Total " << total_snr_points << " SNR points completed in "
              << total_sec << "s (" << (int)(total_sec/60) << " min)" << std::endl;
    std::cout << "Results saved to results/bler_all_mcs_awgn.csv" << std::endl;

    return 0;
}
