#include "common/NrTables.h"
#include <cmath>
#include <algorithm>

namespace nr {

const std::array<McsTableEntry, MAX_MCS_INDEX + 1> MCS_TABLE_1 = {{
    {2, 120},   {2, 157},   {2, 193},   {2, 251},   {2, 308},
    {2, 363},   {2, 449},   {2, 532},   {2, 615},   {2, 658},
    {4, 340},   {4, 378},   {4, 434},   {4, 490},   {4, 553},
    {4, 616},   {4, 658},   {4, 700},   {4, 747},   {4, 797},
    {4, 841},   {4, 895},   {4, 948},   {4, 996},   {6, 806},
    {6, 869},   {6, 932},   {6, 988},   {8, 948}
}};

const int ZC_VALUES[] = {
    2, 4, 8, 16, 32, 64, 128, 256,
    3, 6, 12, 24, 48, 96, 192, 384,
    5, 10, 20, 40, 80, 160, 320,
    7, 14, 28, 56, 112, 224,
    9, 18, 36, 72, 144, 288,
    11, 22, 44, 88, 176, 352,
    13, 26, 52, 104, 208,
    15, 30, 60, 120, 240
};
const int NUM_ZC_VALUES = sizeof(ZC_VALUES) / sizeof(ZC_VALUES[0]);

ModulationScheme mcs_to_modulation(int mcs) {
    if (mcs < 0 || mcs > MAX_MCS_INDEX) {
        return ModulationScheme::QPSK;
    }
    int qm = MCS_TABLE_1[mcs].mod_order;
    switch (qm) {
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

int mod_to_bits_per_symbol(ModulationScheme mod) {
    switch (mod) {
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

LdpcBaseGraphInfo select_ldpc_params(int tb_size, double code_rate) {
    LdpcBaseGraphInfo info;
    info.k_cb = tb_size;
    
    if (tb_size <= 292 || (tb_size <= 3824 && code_rate <= 0.67) || code_rate <= 0.25) {
        info.bgn = 2;
    } else {
        info.bgn = 1;
    }
    
    int k_b;
    if (info.bgn == 1) {
        k_b = 22;
    } else {
        if (tb_size > 640) k_b = 10;
        else if (tb_size > 560) k_b = 9;
        else if (tb_size > 440) k_b = 8;
        else if (tb_size > 380) k_b = 7;
        else k_b = 6;
    }
    
    int zc_min = (tb_size + k_b - 1) / k_b;
    info.zc = ZC_VALUES[NUM_ZC_VALUES - 1];
    for (int i = 0; i < NUM_ZC_VALUES; i++) {
        if (ZC_VALUES[i] >= zc_min) {
            info.zc = ZC_VALUES[i];
            break;
        }
    }
    
    return info;
}

int get_tbs(int mcs, int n_prb, int n_symbols, int n_layers) {
    if (mcs < 0 || mcs > MAX_MCS_INDEX) return 0;
    
    int qm = MCS_TABLE_1[mcs].mod_order;
    double r = MCS_TABLE_1[mcs].code_rate_x1024 / 1024.0;
    
    int n_re_prb = 12 * n_symbols - 8;
    int n_re = std::min(n_re_prb * n_prb, 152 * n_prb);
    
    int n_info = (int)(n_re * r * qm * n_layers);
    
    int tbs;
    if (n_info <= 3824) {
        int n = std::max(3, (int)std::ceil(std::log2(n_info)) - 6);
        int n_info_quant = std::max(24, ((n_info + (1 << n) - 1) >> n) << n);
        
        static const int tbs_table[] = {
            24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152,
            160, 168, 176, 184, 192, 208, 224, 240, 256, 272, 288, 304, 320, 336, 352,
            368, 384, 408, 432, 456, 480, 504, 528, 552, 576, 608, 640, 672, 704, 736,
            768, 808, 848, 888, 928, 984, 1032, 1064, 1128, 1160, 1192, 1224, 1256, 1288,
            1320, 1352, 1416, 1480, 1544, 1608, 1672, 1736, 1800, 1864, 1928, 2024, 2088,
            2152, 2216, 2280, 2408, 2472, 2536, 2600, 2664, 2728, 2792, 2856, 2976, 3104,
            3240, 3368, 3496, 3624, 3752, 3824
        };
        int n_tbs = sizeof(tbs_table) / sizeof(tbs_table[0]);
        tbs = tbs_table[n_tbs - 1];
        for (int i = 0; i < n_tbs; i++) {
            if (tbs_table[i] >= n_info_quant) {
                tbs = tbs_table[i];
                break;
            }
        }
    } else {
        int n = (int)std::ceil(std::log2(n_info - 24)) - 5;
        int n_info_quant = ((n_info - 24 + (1 << n) - 1) >> n) << n + 24;
        
        int C;
        if (n_info_quant > 8424) {
            C = (n_info_quant + 8424 - 1) / 8424;
            tbs = 8 * C * ((n_info_quant + 8 * C - 1) / (8 * C));
        } else {
            C = 1;
            tbs = 8 * ((n_info_quant + 7) / 8);
        }
    }
    
    return tbs;
}

uint32_t gold_sequence(uint32_t state, int n) {
    uint32_t x1 = state & 0x7FFFFFFF;
    uint32_t x2 = 0x00000001;
    
    for (int i = 0; i < n; i++) {
        uint32_t new_x1 = ((x1 >> 24) ^ (x1 >> 20)) & 1;
        x1 = (x1 << 1) | new_x1;
        
        uint32_t new_x2 = ((x2 >> 24) ^ (x2 >> 23) ^ (x2 >> 21) ^ (x2 >> 20)) & 1;
        x2 = (x2 << 1) | new_x2;
    }
    
    return (x1 ^ x2) & 1;
}

void generate_pn_sequence(uint32_t cinit, int length, std::vector<uint8_t>& seq) {
    seq.resize(length);
    const int Nc = 1600;
    
    uint32_t x1 = 1;
    uint32_t x2 = cinit;
    
    for (int n = 0; n < Nc + length; n++) {
        uint32_t x1_new = ((x1 >> 24) ^ (x1 >> 20)) & 1;
        x1 = ((x1 << 1) | x1_new) & 0x7FFFFFFF;
        
        uint32_t x2_new = ((x2 >> 24) ^ (x2 >> 23) ^ (x2 >> 21) ^ (x2 >> 20)) & 1;
        x2 = ((x2 << 1) | x2_new) & 0x7FFFFFFF;
        
        if (n >= Nc) {
            seq[n - Nc] = (x1 ^ x2) & 1;
        }
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
    
    int dmrs_symbols[] = {2};
    int n_dmrs = 1;
    
    if (additional_pos >= 1) {
        dmrs_symbols[n_dmrs++] = (duration == 1) ? 11 : 7;
    }
    if (additional_pos >= 2) {
        dmrs_symbols[n_dmrs++] = (duration == 1) ? 7 : 11;
    }
    if (additional_pos == 3) {
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

} // namespace nr
