#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include <iostream>
#include <iomanip>

using namespace nr;
using namespace nr::phy;

int main() {
    const int n_prb = 3;
    std::cout << std::setw(4) << "MCS" << std::setw(6) << "Qm"
              << std::setw(8) << "R" << std::setw(8) << "TBS"
              << std::setw(8) << "G" << std::setw(6) << "C"
              << std::setw(6) << "bgn" << std::setw(6) << "Zc"
              << std::setw(8) << "K_cb" << std::setw(8) << "cb_k"
              << std::setw(8) << "E" << "\n";
    std::cout << std::string(80, '-') << "\n";

    for (int mcs = 0; mcs <= 27; mcs++) {
        SimulationConfig config;
        config.mcs_index = mcs;
        config.n_rb = n_prb;
        config.n_tx_ant = 1; config.n_rx_ant = 1; config.n_layers = 1;
        config.dmrs_type = DmrsType::TYPE1;
        config.dmrs_additional_pos = 0;
        config.dmrs_duration = 1;
        config.tdd_enabled = false;
        config.mod_scheme = mcs_to_modulation(mcs);
        config.code_rate = mcs_to_code_rate(mcs);
        config.scs = 15000;
        config.max_doppler = 0;
        config.delay_spread = 30e-9;
        config.perfect_csi = true;
        config.channel_type = ChannelType::AWGN;
        config.n_ldpc_iterations = 20;
        config.early_termination = true;
        config.random_seed = 42;

        PdschProcessor proc(config);

        int qm = mcs_to_bits_per_symbol(mcs);
        double R = mcs_to_code_rate(mcs);
        int n_re_per_prb = 13 * 12;
        int G = n_prb * n_re_per_prb * qm;
        int tbs = calculate_tbs(n_prb, n_re_per_prb, qm, 1, R);
        int k_info = tbs;
        int crc_len = get_crc_length(tbs);
        int k_crc = k_info + crc_len;
        auto cbs = compute_cb_segmentation(tbs, k_crc, n_re_per_prb, n_prb, qm, 1, R);

        int bgn = cbs.bgn;
        int zc = cbs.zc;
        int k_b = cbs.k_b;
        int C = cbs.num_cb;
        int cb_k = (C == 1) ? k_crc : cbs.cb_size_with_crc;
        int E_r = (C == 1) ? G : (((G + C - 1) / C + qm - 1) / qm) * qm;
        int n_filler = k_b * zc - cb_k;

        std::cout << std::fixed << std::setprecision(3)
                  << std::setw(4) << mcs
                  << std::setw(6) << qm
                  << std::setw(8) << R
                  << std::setw(8) << tbs
                  << std::setw(8) << G
                  << std::setw(6) << C
                  << std::setw(6) << bgn
                  << std::setw(6) << zc
                  << std::setw(8) << k_b*zc
                  << std::setw(8) << cb_k
                  << std::setw(8) << E_r;
        if (n_filler > 0) std::cout << "  filler=" << n_filler;
        std::cout << "\n";
    }
    return 0;
}
