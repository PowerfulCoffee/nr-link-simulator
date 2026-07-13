#include "common/NrTables.h"
#include "phy/LdpcTables.h"
#include <cmath>
#include <algorithm>

namespace nr {

using namespace nr::ldpc;

const std::array<McsTableEntry, MAX_MCS_INDEX + 1> MCS_TABLE_1 = {{
    {2, 120},   {2, 157},   {2, 193},   {2, 251},   {2, 308},
    {2, 379},   {2, 449},   {2, 526},   {2, 602},   {2, 679},
    {4, 340},   {4, 378},   {4, 434},   {4, 490},   {4, 553},
    {4, 616},   {4, 658},   {6, 438},   {6, 466},   {6, 517},
    {6, 567},   {6, 616},   {6, 666},   {6, 719},   {6, 772},
    {6, 822},   {6, 873},   {6, 910},   {6, 948}
}};

ModulationScheme mcs_to_modulation(int mcs) {
    if (mcs < 0 || mcs > MAX_MCS_INDEX) {
        return ModulationScheme::QPSK;
    }
    int qm = MCS_TABLE_1[mcs].mod_order;
    switch (qm) {
        case 1: return ModulationScheme::BPSK;
        case 2: return ModulationScheme::QPSK;
        case 4: return ModulationScheme::QAM16;
        case 6: return ModulationScheme::QAM64;
        case 8: return ModulationScheme::QAM256;
        default: return ModulationScheme::QPSK;
    }
}

double mcs_to_code_rate(int mcs) {
    if (mcs < 0 || mcs > MAX_MCS_INDEX) {
        return 0.3;
    }
    return MCS_TABLE_1[mcs].code_rate_x1024 / 1024.0;
}

int mcs_to_bits_per_symbol(int mcs) {
    if (mcs < 0 || mcs > MAX_MCS_INDEX) {
        return 2;
    }
    return MCS_TABLE_1[mcs].mod_order;
}

int calculate_num_coded_bits(int n_prb, int n_re_per_prb, int qm, int n_layers) {
    return n_prb * n_re_per_prb * qm * n_layers;
}

int calculate_tbs(int n_prb, int n_re_per_prb, int qm, int n_layers, double target_coderate) {
    int n_re_per_prb_tbs = std::min(156, n_re_per_prb);
    int n_re = n_re_per_prb_tbs * n_prb;
    
    double target_tb_size = target_coderate * n_re * qm * n_layers;
    
    double n_info_q;
    if (target_tb_size <= 3824) {
        double n_small = std::max(3.0, floor(log2(target_tb_size)) - 6);
        n_info_q = std::max(24.0, pow(2, n_small) * floor(target_tb_size / pow(2, n_small)));
    } else {
        double n_large = floor(log2(target_tb_size - 24)) - 5;
        n_info_q = std::max(3840.0, pow(2, n_large) * round((target_tb_size - 24) / pow(2, n_large)));
    }
    
    int num_cb;
    if (n_info_q <= 3824) {
        num_cb = 1;
    } else if (target_coderate <= 0.25) {
        num_cb = static_cast<int>(ceil((n_info_q + 24) / 3816.0));
    } else if (n_info_q > 8424) {
        num_cb = static_cast<int>(ceil((n_info_q + 24) / 8424.0));
    } else {
        num_cb = 1;
    }
    
    static const int tbs_table[] = {
        -1, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128,
        136, 144, 152, 160, 168, 176, 184, 192, 208, 224, 240, 256,
        272, 288, 304, 320, 336, 352, 368, 384, 408, 432, 456, 480,
        504, 528, 552, 576, 608, 640, 672, 704, 736, 768, 808, 848,
        888, 928, 984, 1032, 1064, 1128, 1160, 1192, 1224, 1256,
        1288, 1320, 1352, 1416, 1480, 1544, 1608, 1672, 1736, 1800,
        1864, 1928, 2024, 2088, 2152, 2216, 2280, 2408, 2472, 2536,
        2600, 2664, 2728, 2792, 2856, 2976, 3104, 3240, 3368, 3496,
        3624, 3752, 3824
    };
    const int tbs_table_len = sizeof(tbs_table) / sizeof(tbs_table[0]);
    
    int tb_size;
    if (n_info_q <= 3824) {
        int idx = 0;
        while (idx < tbs_table_len && tbs_table[idx] < static_cast<int>(n_info_q)) {
            idx++;
        }
        if (idx >= tbs_table_len) idx = tbs_table_len - 1;
        tb_size = tbs_table[idx];
    } else {
        tb_size = static_cast<int>(8 * num_cb * ceil((n_info_q + 24) / (8.0 * num_cb)) - 24);
    }
    
    return tb_size;
}

int mod_to_bits_per_symbol(ModulationScheme mod) {
    switch (mod) {
        case ModulationScheme::BPSK: return 1;
        case ModulationScheme::QPSK: return 2;
        case ModulationScheme::QAM16: return 4;
        case ModulationScheme::QAM64: return 6;
        case ModulationScheme::QAM256: return 8;
        default: return 2;
    }
}

int get_crc_length(int tb_size) {
    if (tb_size > 3824) {
        return 24;
    } else {
        return 16;
    }
}

LdpcParams select_ldpc_params(int k_info, double target_coderate) {
    LdpcParams p{};
    p.n_info_bits = k_info;
    
    if (k_info <= 292) {
        p.bgn = 2;
    } else if (k_info <= 3824 && target_coderate <= 0.67) {
        p.bgn = 2;
    } else if (target_coderate <= 0.25) {
        p.bgn = 2;
    } else {
        p.bgn = 1;
    }
    
    if (p.bgn == 1) {
        p.k_b = 22;
    } else {
        if (k_info > 640) p.k_b = 10;
        else if (k_info > 560) p.k_b = 9;
        else if (k_info > 192) p.k_b = 8;
        else p.k_b = 6;
    }
    
    int zc_min = (k_info + p.k_b - 1) / p.k_b;
    
    p.zc = 384;
    for (int s = 0; s < 8; s++) {
        int len = ZC_SET_LENGTHS[s];
        for (int i = 0; i < len; i++) {
            int z = ZC_SETS[s][i];
            if (z >= zc_min && z < p.zc) {
                p.zc = z;
            }
        }
    }
    
    p.k = p.k_b * p.zc;
    if (p.bgn == 1) {
        p.n = 68 * p.zc;
    } else {
        p.n = 52 * p.zc;
    }
    
    return p;
}

void generate_pn_sequence(uint32_t cinit, int length, std::vector<uint8_t>& seq) {
    seq.resize(length);
    const int Nc = 1600;
    
    int N = Nc + length;
    std::vector<uint8_t> x1(N + 31, 0);
    std::vector<uint8_t> x2(N + 31, 0);
    x1[0] = 1;
    for (int i = 0; i < 31; i++) {
        x2[i] = (cinit >> i) & 1;
    }
    
    for (int n = 0; n < N; n++) {
        x1[n + 31] = (x1[n + 3] + x1[n]) % 2;
        x2[n + 31] = (x2[n + 3] + x2[n + 2] + x2[n + 1] + x2[n]) % 2;
    }
    
    for (int n = 0; n < length; n++) {
        seq[n] = (x1[n + Nc] + x2[n + Nc]) % 2;
    }
}

DmrsPattern get_dmrs_pattern(DmrsType type, int additional_pos, int duration) {
    DmrsPattern pattern;
    pattern.type = (type == DmrsType::TYPE1) ? 1 : 2;
    
    for (int i = 0; i < 14; i++) {
        pattern.re_per_prb[i] = 0;
        pattern.cdm_groups_per_prb[i] = 0;
    }
    
    if (type == DmrsType::TYPE1) {
        pattern.n_cdm_groups = 2;
    } else {
        pattern.n_cdm_groups = 3;
    }
    
    int dmrs_symbols[4];
    int n_dmrs = 1;
    dmrs_symbols[0] = 2;
    
    if (additional_pos >= 1) {
        dmrs_symbols[n_dmrs++] = (duration == 1) ? 11 : 7;
    }
    if (additional_pos >= 2) {
        dmrs_symbols[n_dmrs++] = (duration == 1) ? 7 : 11;
    }
    if (additional_pos >= 3) {
        dmrs_symbols[n_dmrs++] = 9;
    }
    
    int re_per_symbol = (type == DmrsType::TYPE1) ? 6 : 4;
    for (int i = 0; i < n_dmrs; i++) {
        int sym = dmrs_symbols[i];
        if (sym >= 0 && sym < 14) {
            pattern.re_per_prb[sym] = re_per_symbol;
            pattern.cdm_groups_per_prb[sym] = pattern.n_cdm_groups;
        }
    }
    
    return pattern;
}

}
