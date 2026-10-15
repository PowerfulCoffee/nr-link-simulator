// Functional check for the Type I CSI precoder (TS 38.214 5.2.1.2):
// 2Tx / 1Rx / rank-1 AWGN link with W = [1, e^{j*pi*m/2}]^T / sqrt2.
// Expected: BLER waterfall similar to the identity (1Tx) case, since the
// precoded power is conserved (unit-norm columns).
#include "phy/PdschProcessor.h"
#include <iostream>
#include <cmath>

using namespace nr;
using namespace nr::phy;

int main() {
    SimulationConfig config;
    config.n_rb = 3;
    config.n_tx_ant = 2;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.mcs_index = 9;
    config.channel_type = ChannelType::AWGN;
    config.max_blocks_per_sinr = 50;
    config.target_block_errors = 50;
    config.random_seed = 42;
    config.perfect_csi = true;
    config.precoder_mode = PrecoderMode::TYPE1_CSI;
    config.csi_n1 = 1;
    config.csi_n2 = 1;
    config.csi_i2 = 1; // phase e^{j*pi/2} -> [1, j]/sqrt2

    PdschProcessor proc(config);

    // Print the precoding matrix
    ComplexMat w = build_type1_csi_precoder(2, 1, config.csi_n1, config.csi_n2,
                                            config.csi_i1_phi, config.csi_i11,
                                            config.csi_i12, config.csi_i2);
    std::cout << "Type1 W (2x1): ";
    for (int p = 0; p < 2; p++) {
        std::cout << "(" << w(p, 0).real() << ", " << w(p, 0).imag() << ") ";
    }
    std::cout << std::endl;

    BlerResult r;
    double snrs[] = {0.0, 3.0, 4.0, 5.0, 6.0, 9.0};
    for (double snr : snrs) {
        proc.process_single_snr_point(snr, r);
        std::cout << "SNR=" << snr << " dB: BLER=" << r.bler
                  << " (" << r.n_errors << "/" << r.n_blocks << ")" << std::endl;
    }
    return 0;
}
