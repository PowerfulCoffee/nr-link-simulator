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
    config.max_blocks_per_sinr = 1;
    config.target_block_errors = 5;
    config.sinr_start = 100.0;
    config.sinr_end = 100.0;
    config.sinr_step = 2.0;
    config.channel_type = ChannelType::AWGN;
    config.mod_scheme = mcs_to_modulation(config.mcs_index);
    config.code_rate = mcs_to_code_rate(config.mcs_index);
    config.n_ldpc_iterations = 50;

    std::cout << "Debug SISO Test: Direct loopback (no channel noise)\n";
    std::cout << "MCS " << config.mcs_index << ", " << config.n_rb << " RBs\n";
    std::cout << "Modulation: QPSK, Code rate: " << config.code_rate << "\n\n";

    SimulationConfig sim_config = config;
    sim_config.n_tx_ant = 1;
    sim_config.n_layers = 1;
    sim_config.n_rx_ant = 1;

    PdschProcessor processor(sim_config);

    TransportBlock tb = processor.generate_transport_block();
    std::cout << "TB size: " << tb.tb_size << " bits\n";
    std::cout << "First 20 TB bits: ";
    for (int i = 0; i < std::min(20, tb.tb_size); i++) {
        std::cout << (int)tb.bits(i);
    }
    std::cout << "\n";

    PdschTxResult tx_res = processor.transmit(tb, 0);
    std::cout << "LDPC: bgn=" << tx_res.bgn << ", zc=" << tx_res.zc << "\n";
    std::cout << "E (rate matched bits): " << tx_res.n_coded_bits << "\n";
    std::cout << "Scrambling seed: " << tx_res.scrambling_seed << "\n";

    PdschRxResult rx_res = processor.receive(tx_res.tx_grid, tx_res, 20.0, 0);

    std::cout << "\nCRC ok: " << (rx_res.crc_ok ? "YES" : "NO") << "\n";
    if (rx_res.crc_ok) {
        std::cout << "First 20 decoded bits: ";
        for (int i = 0; i < std::min(20, (int)rx_res.decoded_bits.n_elem); i++) {
            std::cout << (int)rx_res.decoded_bits(i);
        }
        std::cout << "\n";

        bool bits_match = true;
        int min_len = std::min(tb.tb_size, (int)rx_res.decoded_bits.n_elem);
        for (int i = 0; i < min_len; i++) {
            if (tb.bits(i) != rx_res.decoded_bits(i)) {
                bits_match = false;
                std::cout << "Mismatch at bit " << i << ": tx=" << (int)tb.bits(i)
                          << ", rx=" << (int)rx_res.decoded_bits(i) << "\n";
                if (i > 10) break;
            }
        }
        if (bits_match) {
            std::cout << "All bits match!\n";
        }
    }

    return rx_res.crc_ok ? 0 : 1;
}
