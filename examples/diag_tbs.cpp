#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include <iostream>
#include <iomanip>

using namespace nr;
using namespace nr::phy;

int main(int argc, char** argv) {
    int mcs = 5;
    if (argc > 1) mcs = atoi(argv[1]);
    
    SimulationConfig config;
    config.mcs_index = mcs;
    config.n_rb = 25;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.max_blocks_per_sinr = 1;
    config.target_block_errors = 1;
    config.sinr_start = 10;
    config.sinr_end = 10;
    config.channel_type = ChannelType::AWGN;
    config.mod_scheme = mcs_to_modulation(mcs);
    config.code_rate = mcs_to_code_rate(mcs);
    config.n_ldpc_iterations = 20;
    config.perfect_csi = false;
    
    DmrsPattern pat = get_dmrs_pattern(DmrsType::TYPE1, 0, 1);
    int n_dmrs_syms = 0;
    int dmrs_re_per_prb = 0;
    for (int s = 0; s < 14; s++) {
        if (pat.re_per_prb[s] > 0) {
            n_dmrs_syms++;
            dmrs_re_per_prb = pat.re_per_prb[s];
        }
    }
    std::cout << "DMRS: " << n_dmrs_syms << " symbols, " << dmrs_re_per_prb << " RE/PRB/symbol\n";
    
    int n_re_per_prb = 0;
    for (int s = 0; s < 14; s++) {
        if (pat.re_per_prb[s] == 0) {
            n_re_per_prb += 12;
        } else {
            n_re_per_prb += (12 - pat.re_per_prb[s]);
        }
    }
    std::cout << "PDSCH RE per PRB (including non-DMRS RE on DMRS symbols): " << n_re_per_prb << "\n";
    
    int n_re_per_prb_current = 0;
    for (int s = 0; s < 14; s++) {
        if (pat.re_per_prb[s] == 0) n_re_per_prb_current += 12;
    }
    std::cout << "PDSCH RE per PRB (current, excluding DMRS symbols entirely): " << n_re_per_prb_current << "\n";
    
    PdschProcessor proc(config);
    auto tx = proc.transmit(0);
    
    int qm = tx.qm;
    int n_coded = tx.n_coded_bits;
    int n_info = tx.n_info_bits;
    std::cout << "\n=== C++ (current code) ===\n";
    std::cout << "MCS " << mcs << ": Qm=" << qm << ", R=" << config.code_rate << "\n";
    std::cout << "n_info_bits (TBS+CRC)=" << n_info << "\n";
    std::cout << "N_coded_bits=" << n_coded << "\n";
    std::cout << "N_data_RE=" << n_coded/qm << "\n";
    std::cout << "N_re_per_prb=" << n_coded/qm/config.n_rb << "\n";
    
    return 0;
}
