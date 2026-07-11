#pragma once

#include "common/Types.h"
#include <array>

namespace nr {

struct McsTableEntry {
    int mod_order;
    double code_rate_x1024;
};

struct LdpcBaseGraphInfo {
    int bgn;
    int zc;
    int k_cb;
};

constexpr int MAX_MCS_INDEX = 28;

extern const std::array<McsTableEntry, MAX_MCS_INDEX + 1> MCS_TABLE_1;

constexpr int MAX_ZC = 384;
extern const int ZC_VALUES[];
extern const int NUM_ZC_VALUES;

LdpcBaseGraphInfo select_ldpc_params(int tb_size, double code_rate);
int get_tbs(int mcs, int n_prb, int n_symbols, int n_layers);
int get_crc_length(int tb_size);
ModulationScheme mcs_to_modulation(int mcs);
double mcs_to_code_rate(int mcs);
int mod_to_bits_per_symbol(ModulationScheme mod);

constexpr int get_rb_size(int scs) {
    return 12;
}

constexpr double get_slot_duration(int scs) {
    return 1.0e-3 / (scs / 15e3);
}

constexpr int get_symbols_per_slot() {
    return 14;
}

constexpr int get_fft_size(int n_rb, int scs) {
    int n_sc = n_rb * 12;
    int fft_size = 1;
    while (fft_size < n_sc) fft_size <<= 1;
    return fft_size;
}

struct DmrsPattern {
    int type;
    int n_cdm_groups;
    int re_per_prb[14];
    int cdm_groups_per_prb[14];
};

DmrsPattern get_dmrs_pattern(DmrsType type, int additional_pos, int duration);

uint32_t gold_sequence(uint32_t state, int n);
void generate_pn_sequence(uint32_t cinit, int length, std::vector<uint8_t>& seq);

struct BlerResult {
    double sinr_db;
    int n_blocks;
    int n_errors;
    double bler;
    
    BlerResult() : sinr_db(0), n_blocks(0), n_errors(0), bler(0) {}
};

} // namespace nr
