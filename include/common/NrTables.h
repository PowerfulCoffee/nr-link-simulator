#pragma once

#include "common/Types.h"
#include <array>

namespace nr {

struct McsTableEntry {
    int mod_order;
    double code_rate_x1024;
};

struct LdpcParams {
    int bgn;
    int zc;
    int k_b;
    int k;
    int n;
    int n_info_bits;
};

struct CodeBlockSegParams {
    int num_cb;
    int tb_crc_len;
    int cb_crc_len;
    int cb_info_bits;
    int cb_size_with_crc;
    int cb_k;
    int bgn;
    int zc;
    int k_b;
    int cw_length;
};

constexpr int MAX_MCS_INDEX = 28;

extern const std::array<McsTableEntry, MAX_MCS_INDEX + 1> MCS_TABLE_1;

LdpcParams select_ldpc_params(int k_info, double target_coderate);
int get_crc_length(int tb_size);
ModulationScheme mcs_to_modulation(int mcs);
double mcs_to_code_rate(int mcs);
int mcs_to_bits_per_symbol(int mcs);
int mod_to_bits_per_symbol(ModulationScheme mod);

int calculate_tbs(int n_prb, int n_re_per_prb, int qm, int n_layers, double target_coderate);
int calculate_num_coded_bits(int n_prb, int n_re_per_prb, int qm, int n_layers);
CodeBlockSegParams compute_cb_segmentation(int tb_size, int n_info_bits_after_tb_crc,
                                           int n_re_per_prb, int n_prb,
                                           int qm, int n_layers, double target_coderate);

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

void generate_pn_sequence(uint32_t cinit, int length, std::vector<uint8_t>& seq);

struct BlerResult {
    double sinr_db;
    int n_blocks;
    int n_errors;
    double bler;
    
    BlerResult() : sinr_db(0), n_blocks(0), n_errors(0), bler(0) {}
};

} // namespace nr
