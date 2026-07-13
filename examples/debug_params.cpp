#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include <iostream>
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
    config.channel_type = ChannelType::AWGN;
    config.mod_scheme = mcs_to_modulation(config.mcs_index);
    config.code_rate = mcs_to_code_rate(config.mcs_index);
    config.n_ldpc_iterations = 20;
    
    std::cout << "Configuration:\n";
    std::cout << "  MCS = " << config.mcs_index << "\n";
    std::cout << "  Modulation = ";
    switch(config.mod_scheme) {
        case ModulationScheme::QPSK: std::cout << "QPSK"; break;
        case ModulationScheme::QAM16: std::cout << "16QAM"; break;
        case ModulationScheme::QAM64: std::cout << "64QAM"; break;
        default: std::cout << "Other";
    }
    std::cout << "\n";
    std::cout << "  Target code rate = " << config.code_rate << "\n";
    std::cout << "  PRBs = " << config.n_rb << "\n";
    
    int n_re_per_prb = 13 * 12;
    int qm = mcs_to_bits_per_symbol(config.mcs_index);
    int E = calculate_num_coded_bits(config.n_rb, n_re_per_prb, qm, config.n_layers);
    int tbs = calculate_tbs(config.n_rb, n_re_per_prb, qm, config.n_layers, config.code_rate);
    
    std::cout << "\nPDSCH Resource Calculation:\n";
    std::cout << "  RE per PRB (13 data symbols) = " << n_re_per_prb << "\n";
    std::cout << "  Bits per symbol (Qm) = " << qm << "\n";
    std::cout << "  Total coded bits (E) = " << E << "\n";
    std::cout << "  TBS (3GPP standard) = " << tbs << " bits\n";
    
    PdschProcessor proc(config);
    
    TransportBlock tb = proc.generate_transport_block();
    std::cout << "  Actual TB size from generate = " << tb.tb_size << " bits\n";
    
    PdschTxResult tx_res = proc.transmit(tb, 0);
    
    std::cout << "\nLDPC Parameters:\n";
    std::cout << "  TB + CRC bits (K_info) = " << tx_res.tb_bits_after_crc.size() << "\n";
    std::cout << "  BGN = " << tx_res.bgn << "\n";
    std::cout << "  Zc = " << tx_res.zc << "\n";
    
    LdpcParams ldpc_info = select_ldpc_params(tx_res.tb_bits_after_crc.size(), config.code_rate);
    std::cout << "  K (systematic bits) = " << ldpc_info.k << "\n";
    std::cout << "  N (codeword size) = " << ldpc_info.n << "\n";
    std::cout << "  N_cb (after 2*Zc puncturing) = " << ldpc_info.n - 2*ldpc_info.zc << "\n";
    
    double actual_code_rate = static_cast<double>(tb.tb_size) / tx_res.n_coded_bits;
    std::cout << "\nActual code rate = " << actual_code_rate << " (info bits / transmitted bits)\n";
    
    std::cout << "\nRunning BLER simulation with this configuration...\n";
    std::cout << "QPSK at rate ~0.37 typically works around 0-2 dB in AWGN\n";
    
    return 0;
}
