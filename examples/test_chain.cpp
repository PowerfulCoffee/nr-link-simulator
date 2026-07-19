#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace nr;
using namespace nr::phy;

int main() {
    SimulationConfig config;
    config.mcs_index = 27;
    config.n_rb = 3;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.channel_type = ChannelType::AWGN;
    config.n_ldpc_iterations = 20;
    config.early_termination = true;

    std::cout << "MCS " << (int)config.mcs_index << std::endl;

    PdschProcessor proc(config);

    TransportBlock tb = proc.generate_transport_block();
    std::cout << "TB bits: " << tb.tb_size << std::endl;

    PdschTxResult tx_res = proc.transmit(tb, 0);
    std::cout << "Tx n_info_bits=" << tx_res.n_info_bits
              << " bgn=" << (int)tx_res.bgn
              << " zc=" << tx_res.zc
              << " k_b=" << tx_res.k_b
              << " n_coded_bits=" << tx_res.n_coded_bits
              << std::endl;

    for (int snr_db : {40, 60, 80, 100}) {
        int n_ok = 0;
        int N = 10;
        for (int b = 0; b < N; b++) {
            TransportBlock tb2 = proc.generate_transport_block();
            PdschTxResult tx2 = proc.transmit(tb2, b+1);
            PdschRxResult rx2 = proc.receive(tx2.tx_grid, tx2, snr_db, b+1);
            if (rx2.crc_ok) n_ok++;
        }
        std::cout << "SNR=" << snr_db << "dB: " << n_ok << "/" << N
                  << " passed" << std::endl;
    }

    return 0;
}
