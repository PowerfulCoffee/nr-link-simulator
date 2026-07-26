#include "common/Types.h"
#include "common/NrTables.h"
#include <iostream>
#include <iomanip>

using namespace nr;

int main(int argc, char* argv[]) {
    int mcs = 3;
    int n_prb = 25;
    int n_layers = 1;
    int dmrs_add_pos = 1;
    
    if (argc > 1) mcs = atoi(argv[1]);
    if (argc > 2) n_prb = atoi(argv[2]);
    if (argc > 3) n_layers = atoi(argv[3]);
    if (argc > 4) dmrs_add_pos = atoi(argv[4]);
    
    int qm = mcs_to_bits_per_symbol(mcs);
    double R = mcs_to_code_rate(mcs);
    
    int n_dmrs_sym = (dmrs_add_pos == 0) ? 1 : (dmrs_add_pos == 1 ? 2 : (dmrs_add_pos == 2 ? 3 : 4));
    int n_pdsch_sym = 14 - n_dmrs_sym;
    int n_re_per_prb = n_pdsch_sym * 12;
    
    int tbs = calculate_tbs(n_prb, n_re_per_prb, qm, n_layers, R);
    
    int n_re_total = n_pdsch_sym * n_prb * 12;
    int n_coded = (n_re_total * qm * n_layers / qm) * qm;
    if (n_coded < (tbs + 24)) n_coded = ((tbs + 24 + qm - 1) / qm) * qm;
    
    int n_info = tbs;
    
    std::cout << "MCS=" << mcs << " (" << qm << "QAM), n_prb=" << n_prb 
              << ", layers=" << n_layers << ", dmrs_add_pos=" << dmrs_add_pos << std::endl;
    std::cout << "qm=" << qm << std::endl;
    std::cout << "target_code_rate=" << std::fixed << std::setprecision(6) << R << std::endl;
    std::cout << "n_pdsch_sym=" << n_pdsch_sym << std::endl;
    std::cout << "n_re_per_prb=" << n_re_per_prb << std::endl;
    std::cout << "tbs(k_info)=" << tbs << std::endl;
    std::cout << "k_crc=" << tbs + 24 << std::endl;
    std::cout << "n_coded=" << n_coded << std::endl;
    std::cout << "actual_rate=" << std::fixed << std::setprecision(4) << (double)(tbs + 24)/n_coded << std::endl;
    
    return 0;
}
