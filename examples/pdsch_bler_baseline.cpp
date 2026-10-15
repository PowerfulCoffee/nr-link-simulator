/*
 * Baseline PDSCH link curve sweep (BEFORE any code fixes).
 * Full TX chain (TB->CRC->LDPC->RateMatch->Scramble->Modulate) + AWGN
 * + RX chain (Equalize->Demodulate->Descramble->RateRecover->LDPC->CRC),
 * AWGN channel with perfect CSI (fast path, SISO).
 * Output: results/baseline/baseline_bler_awgn.csv
 */
#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <filesystem>

using namespace nr;
using namespace nr::phy;
using namespace std::chrono;
namespace fs = std::filesystem;

int main(int argc, char** argv) {
    // Representative MCS subset spanning QPSK/16QAM/64QAM/256QAM
    const std::vector<int> mcs_list = {0, 2, 4, 6, 9, 11, 13, 16, 19, 21, 24, 26, 27};

    // argv[1]: csi mode "perfect" (default) or "ls" (LS channel estimation,
    //          exercises DMRS generation + LS noise estimation)
    // argv[2]: output csv path (default results/baseline/baseline_bler_awgn.csv)
    bool perfect_csi = true;
    if (argc > 1 && std::string(argv[1]) == "ls") perfect_csi = false;
    std::string out_csv = (argc > 2) ? argv[2]
                                     : std::string("results/baseline/baseline_bler_awgn.csv");

    const int n_prb = 3;
    const int n_re_per_prb = 156;
    const double sinr_step = 1.0;
    const int max_blocks_per_snr = 2000;
    const int target_block_errors = 60;
    const int n_ldpc_iters = 25;

    std::cout << "BASELINE PDSCH sweep: " << mcs_list.size() << " MCS, nPRB=" << n_prb
              << ", AWGN, " << (perfect_csi ? "perfect CSI" : "LS CSI") << ", SISO" << std::endl;

    fs::create_directories("results/baseline");
    std::ofstream csv(out_csv);
    csv << "mcs,qm,rate,tbs,snr_db,bler,n_blocks,n_errors\n";

    auto t0 = high_resolution_clock::now();
    int total_pts = 0;

    for (int mcs : mcs_list) {
        int qm = mcs_to_bits_per_symbol(mcs);
        double R = mcs_to_code_rate(mcs);
        int tbs = calculate_tbs(n_prb, n_re_per_prb, qm, 1, R);

        // per-MCS SNR range (same logic as pdsch_bler_all_mcs)
        double snr_start, snr_end;
        if (qm == 2) {
            if (R < 0.20)      { snr_start = -10.0; snr_end = 6.0; }
            else if (R < 0.45) { snr_start = -6.0;  snr_end = 10.0; }
            else               { snr_start = -4.0;  snr_end = 12.0; }
        } else if (qm == 4) {
            if (R < 0.50)      { snr_start = 0.0;   snr_end = 16.0; }
            else               { snr_start = 2.0;   snr_end = 18.0; }
        } else if (qm == 6) {
            if (R < 0.60)      { snr_start = 4.0;   snr_end = 22.0; }
            else if (R < 0.76) { snr_start = 8.0;   snr_end = 24.0; }
            else               { snr_start = 10.0;  snr_end = 26.0; }
        } else {
            if (R < 0.75)      { snr_start = 12.0;  snr_end = 28.0; }
            else if (R < 0.87) { snr_start = 16.0;  snr_end = 32.0; }
            else               { snr_start = 18.0;  snr_end = 34.0; }
        }
        int n_pts = static_cast<int>((snr_end - snr_start) / sinr_step) + 1;

        SimulationConfig config;
        config.mcs_index = mcs;
        config.n_rb = n_prb;
        config.n_tx_ant = 1;
        config.n_rx_ant = 1;
        config.n_layers = 1;
        config.mod_scheme = mcs_to_modulation(mcs);
        config.code_rate = R;
        config.channel_type = ChannelType::AWGN;
        config.perfect_csi = perfect_csi;
        config.max_blocks_per_sinr = max_blocks_per_snr;
        config.target_block_errors = target_block_errors;
        config.sinr_start = snr_start;
        config.sinr_end = snr_end;
        config.sinr_step = sinr_step;
        config.n_sinr_points = n_pts;
        config.n_ldpc_iterations = n_ldpc_iters;
        config.early_termination = true;
        config.random_seed = 123;
        config.scs = 15000;

        auto results = run_bler_simulation(config, nullptr, "perfect");

        std::cout << "MCS " << std::setw(2) << mcs << " (Qm=" << qm << ", R="
                  << std::fixed << std::setprecision(3) << R << ", TBS=" << tbs << "):";
        for (const auto& res : results) {
            std::cout << " " << std::setprecision(1) << res.sinr_db << ":"
                      << std::setprecision(3) << res.bler;
            csv << mcs << "," << qm << "," << std::fixed << std::setprecision(4) << R
                << "," << tbs << "," << std::setprecision(1) << res.sinr_db
                << "," << std::setprecision(6) << res.bler << ","
                << res.n_blocks << "," << res.n_errors << "\n";
            total_pts++;
        }
        csv.flush();
        std::cout << std::endl;
    }

    double sec = duration_cast<seconds>(high_resolution_clock::now() - t0).count();
    std::cout << "Done: " << total_pts << " SNR points in " << sec << " s" << std::endl;
    std::cout << "Saved: " << out_csv << std::endl;
    return 0;
}
