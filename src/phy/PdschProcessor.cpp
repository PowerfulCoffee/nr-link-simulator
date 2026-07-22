#include "phy/PdschProcessor.h"
#include "phy/ModuleFactory.h"
#include "common/NrTables.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <armadillo>
#include <memory>
#include <random>

namespace nr {
namespace phy {

namespace {

constexpr int NUM_SYMBOLS_PER_SLOT = 14;
constexpr double LLR_CLIP = 20.0;

DmrsPattern get_active_dmrs_pattern(const SimulationConfig& config) {
    return get_dmrs_pattern(config.dmrs_type, config.dmrs_additional_pos, config.dmrs_duration);
}

bool is_dmrs_symbol(int sym, const DmrsPattern& pattern) {
    return pattern.re_per_prb[sym] > 0;
}

bool is_pdsch_symbol(int sym, const DmrsPattern& pattern, const SimulationConfig& config, int slot_idx) {
    if (is_dmrs_symbol(sym, pattern)) return false;
    if (!config.tdd_enabled) return true;

    int pat_idx = slot_idx % TddConfig::slots_per_frame;
    SlotType st = config.tdd_config.slot_pattern[pat_idx];

    if (st == SlotType::DOWNLINK) return true;
    if (st == SlotType::UPLINK) return false;

    if (st == SlotType::SPECIAL) {
        return sym < config.tdd_config.s_slot_dl_symbols;
    }
    return true;
}

int count_pdsch_re(int n_rb, const DmrsPattern& pattern, const SimulationConfig& config, int slot_idx) {
    int count = 0;
    for (int sym = 0; sym < NUM_SYMBOLS_PER_SLOT; sym++) {
        if (is_pdsch_symbol(sym, pattern, config, slot_idx)) {
            count += n_rb * 12;
        }
    }
    return count;
}

int count_pdsch_re_per_prb(const DmrsPattern& pattern, const SimulationConfig& config, int slot_idx) {
    int count = 0;
    for (int sym = 0; sym < NUM_SYMBOLS_PER_SLOT; sym++) {
        if (is_pdsch_symbol(sym, pattern, config, slot_idx)) {
            count += 12;
        }
    }
    return count;
}

void map_pdsch_to_grid(ResourceGrid& grid, const ComplexMat& precoded,
                       int rb_start, int rb_len, const DmrsPattern& pattern,
                       const SimulationConfig& config, int slot_idx) {
    int n_ant = grid.n_ant;
    int re_idx = 0;
    int total_re = static_cast<int>(precoded.n_rows);

    for (int sym = 0; sym < NUM_SYMBOLS_PER_SLOT; sym++) {
        if (!is_pdsch_symbol(sym, pattern, config, slot_idx)) {
            continue;
        }

        for (int rb = rb_start; rb < rb_start + rb_len; rb++) {
            for (int sc_in_rb = 0; sc_in_rb < 12; sc_in_rb++) {
                int sc = rb * 12 + sc_in_rb;
                if (sc >= grid.n_subcarriers) {
                    continue;
                }
                if (re_idx >= total_re) {
                    return;
                }

                for (int ant = 0; ant < n_ant; ant++) {
                    if (ant < static_cast<int>(precoded.n_cols)) {
                        grid.set_re(ant, sym, sc, precoded(re_idx, ant));
                    }
                }
                re_idx++;
            }
        }
    }
}

ComplexMat extract_pdsch_from_grid(const ResourceGrid& grid,
                                    int rb_start, int rb_len, const DmrsPattern& pattern,
                                    const SimulationConfig& config, int slot_idx) {
    int n_ant = grid.n_ant;
    int total_re = count_pdsch_re(rb_len, pattern, config, slot_idx);

    ComplexMat extracted(total_re, n_ant, arma::fill::zeros);
    int re_idx = 0;

    for (int sym = 0; sym < NUM_SYMBOLS_PER_SLOT; sym++) {
        if (!is_pdsch_symbol(sym, pattern, config, slot_idx)) {
            continue;
        }

        for (int rb = rb_start; rb < rb_start + rb_len; rb++) {
            for (int sc_in_rb = 0; sc_in_rb < 12; sc_in_rb++) {
                int sc = rb * 12 + sc_in_rb;
                if (sc >= grid.n_subcarriers) {
                    continue;
                }
                if (re_idx >= total_re) {
                    continue;
                }

                for (int ant = 0; ant < n_ant; ant++) {
                    extracted(re_idx, ant) = grid.get_re(ant, sym, sc);
                }
                re_idx++;
            }
        }
    }

    return extracted;
}

ComplexCube extract_channel_at_pdsch(const ComplexCube& h_est,
                                      int n_rb,
                                      int n_sc_total,
                                      const DmrsPattern& pattern,
                                      const SimulationConfig& config, int slot_idx) {
    int n_sc = static_cast<int>(h_est.n_rows);
    int n_sym = static_cast<int>(h_est.n_cols);
    int n_ch = static_cast<int>(h_est.n_slices);
    int total_re = count_pdsch_re(n_rb, pattern, config, slot_idx);

    ComplexCube h_at_pdsch(total_re, 1, n_ch, arma::fill::zeros);
    int re_idx = 0;

    for (int sym = 0; sym < NUM_SYMBOLS_PER_SLOT; sym++) {
        if (!is_pdsch_symbol(sym, pattern, config, slot_idx)) {
            continue;
        }
        int sym_use = sym;
        if (sym_use >= n_sym) {
            sym_use = n_sym - 1;
        }

        for (int rb = 0; rb < n_rb; rb++) {
            for (int sc_in_rb = 0; sc_in_rb < 12; sc_in_rb++) {
                int sc = rb * 12 + sc_in_rb;
                if (sc >= n_sc || sc >= n_sc_total) {
                    continue;
                }
                if (re_idx >= total_re) {
                    continue;
                }

                for (int c = 0; c < n_ch; c++) {
                    h_at_pdsch(re_idx, 0, c) = h_est(sc, sym_use, c);
                }
                re_idx++;
            }
        }
    }

    return h_at_pdsch;
}

} // anonymous namespace

PdschProcessor::PdschProcessor(const SimulationConfig& config)
    : config_(config), seed_(config.random_seed), rng_(config.random_seed) {
    n_pdsch_rbs_ = config_.n_rb;
    dmrs_pattern_ = get_active_dmrs_pattern(config_);
    n_dmrs_symbols_ = 0;
    for (int s = 0; s < NUM_SYMBOLS_PER_SLOT; s++) {
        if (dmrs_pattern_.re_per_prb[s] > 0) n_dmrs_symbols_++;
    }

    config_.n_tx_ant = std::min(config_.n_tx_ant, 4);
    config_.n_layers = std::min(config_.n_layers, config_.n_tx_ant);
    config_.n_rx_ant = std::min(config_.n_rx_ant, config_.n_tx_ant > 1 ? 4 : 1);

    if (config_.mcs_index >= 0 && config_.mcs_index <= MAX_MCS_INDEX) {
        config_.mod_scheme = mcs_to_modulation(config_.mcs_index);
        config_.code_rate = mcs_to_code_rate(config_.mcs_index);
    }

    pdsch_sym_alloc_.start = 0;
    pdsch_sym_alloc_.length = NUM_SYMBOLS_PER_SLOT;

    init_default_modules();

    int ref_slot = 0;
    n_pdsch_symbols_ = 0;
    for (int s = 0; s < NUM_SYMBOLS_PER_SLOT; s++) {
        if (is_pdsch_symbol(s, dmrs_pattern_, config_, ref_slot)) n_pdsch_symbols_++;
    }
}

void PdschProcessor::set_channel_estimator(std::unique_ptr<IChannelEstimator> estimator) {
    channel_estimator_ = std::move(estimator);
}

void PdschProcessor::set_channel_model(std::unique_ptr<channel::IChannelModel> channel) {
    channel_model_ = std::move(channel);
}

void PdschProcessor::set_seed(uint64_t seed) {
    seed_ = seed;
    rng_.seed(seed);
}

void PdschProcessor::init_default_modules() {
    crc_encoder_ = create_crc_encoder();
    ldpc_encoder_ = create_ldpc_encoder();
    ldpc_decoder_ = create_ldpc_decoder();
    rate_matcher_ = create_rate_matcher();
    scrambler_ = create_scrambler();
    modulator_ = create_modulator();
    layer_mapper_ = create_layer_mapper();
    precoder_ = create_precoder();
    dmrs_generator_ = create_dmrs_generator();
    resource_mapper_ = create_resource_mapper();
    ofdm_modulator_ = create_ofdm_modulator();
    equalizer_ = create_mmse_equalizer();
    channel_estimator_ = create_ls_channel_estimator();
    channel_model_ = channel::create_channel(config_.channel_type);
    channel_model_->set_config(config_);
}

int PdschProcessor::calculate_pdsch_capacity() {
    return count_pdsch_re(n_pdsch_rbs_, dmrs_pattern_, config_, 0);
}

ComplexMat PdschProcessor::get_identity_precoding_matrix() {
    int n_ports = config_.n_tx_ant;
    int n_layers = config_.n_layers;
    ComplexMat w(n_ports, n_layers, arma::fill::zeros);
    for (int l = 0; l < n_layers; l++) {
        w(l, l) = Complex(1.0, 0.0);
    }
    return w;
}

TransportBlock PdschProcessor::generate_transport_block() {
    int n_re_per_prb = count_pdsch_re_per_prb(dmrs_pattern_, config_, 0);
    int qm = mcs_to_bits_per_symbol(config_.mcs_index);
    double target_coderate = mcs_to_code_rate(config_.mcs_index);
    int n_layers = config_.n_layers;
    
    int tbs = calculate_tbs(n_pdsch_rbs_, n_re_per_prb, qm, n_layers, target_coderate);

    TransportBlock tb(tbs);
    tb.mcs = config_.mcs_index;
    tb.rv = 0;
    tb.generate_random_bits(seed_++);
    return tb;
}

PdschTxResult PdschProcessor::transmit(const TransportBlock& tb, int slot_idx) {
    PdschTxResult result;
    result.n_info_bits = tb.tb_size;

    int n_layers = config_.n_layers;
    int n_ports = config_.n_tx_ant;
    int n_sc = n_pdsch_rbs_ * 12;

    int crc_len = get_crc_length(tb.tb_size);
    BitVec bits_with_crc = crc_encoder_->encode(tb.bits, crc_len);
    result.tb_bits_after_crc = bits_with_crc;

    int k_info = static_cast<int>(bits_with_crc.size());
    double target_coderate = mcs_to_code_rate(config_.mcs_index);
    int qm = mcs_to_bits_per_symbol(config_.mcs_index);
    int n_re_per_prb = count_pdsch_re_per_prb(dmrs_pattern_, config_, slot_idx);
    int G = calculate_num_coded_bits(n_pdsch_rbs_, n_re_per_prb, qm, n_layers);

    CodeBlockSegParams cbs = compute_cb_segmentation(tb.tb_size, k_info,
                                                      n_re_per_prb, n_pdsch_rbs_,
                                                      qm, n_layers, target_coderate);
    result.bgn = cbs.bgn;
    result.zc = cbs.zc;
    result.k_b = cbs.k_b;
    result.qm = qm;
    result.num_cb = cbs.num_cb;
    result.cb_crc_len = cbs.cb_crc_len;
    result.cb_info_bits = cbs.cb_info_bits;
    result.cb_size_with_crc = cbs.cb_size_with_crc;

    BitVec rm_bits_concat;
    result.cb_info.clear();

    int cb_offset = 0;
    int E_r;
    if (cbs.num_cb == 1) {
        E_r = G;
    } else {
        E_r = (G + cbs.num_cb - 1) / cbs.num_cb;
        E_r = ((E_r + qm - 1) / qm) * qm;
    }
    result.cb_e_bits = E_r;

    for (int c = 0; c < cbs.num_cb; c++) {
        BitVec cb_bits;
        int k_cb;

        if (cbs.num_cb == 1) {
            cb_bits = bits_with_crc;
            k_cb = k_info;
        } else {
            int payload_per_cb = cbs.cb_info_bits;
            cb_bits.resize(payload_per_cb, 0);
            int start_bit = c * payload_per_cb;
            int end_bit = std::min((c + 1) * payload_per_cb, k_info);
            int copy_len = end_bit - start_bit;
            for (int i = 0; i < copy_len; i++) {
                cb_bits[i] = bits_with_crc[start_bit + i];
            }
            for (int i = copy_len; i < payload_per_cb; i++) {
                cb_bits[i] = 0;
            }
            BitVec cb_with_crc = crc_encoder_->encode(cb_bits, cbs.cb_crc_len);
            cb_bits = cb_with_crc;
            k_cb = static_cast<int>(cb_bits.size());
            cb_bits.resize(cbs.cb_k, 0);
        }

        BitVec coded_bits_full = ldpc_encoder_->encode(cb_bits, cbs.bgn, cbs.zc);

        int n_punctured = 2 * cbs.zc;
        BitVec coded_bits(coded_bits_full.begin() + n_punctured, coded_bits_full.end());

        int n_filler = cbs.cb_k - k_cb;
        int filler_start = k_cb - n_punctured;
        BitVec coded_bits_comp;
        coded_bits_comp.reserve(coded_bits.size() - n_filler);
        if (n_filler > 0 && filler_start >= 0 && filler_start + n_filler <= (int)coded_bits.size()) {
            coded_bits_comp.insert(coded_bits_comp.end(), coded_bits.begin(), coded_bits.begin() + filler_start);
            coded_bits_comp.insert(coded_bits_comp.end(), coded_bits.begin() + filler_start + n_filler, coded_bits.end());
        } else {
            coded_bits_comp = coded_bits;
        }

        int this_E = E_r;
        if (c == cbs.num_cb - 1 && cbs.num_cb > 1) {
            this_E = G - (cbs.num_cb - 1) * E_r;
            this_E = std::max(this_E, qm);
        }

        BitVec rm_bits_cb = rate_matcher_->rate_match(coded_bits_comp, this_E, tb.rv, cbs.bgn, cbs.zc, qm, 0);

        CodeBlockInfo cbi;
        cbi.offset = cb_offset;
        cbi.length = this_E;
        cbi.e_bits = this_E;
        result.cb_info.push_back(cbi);

        rm_bits_concat.insert(rm_bits_concat.end(), rm_bits_cb.begin(), rm_bits_cb.end());
        cb_offset += this_E;
    }

    result.n_coded_bits = static_cast<int>(rm_bits_concat.size());
    if (result.n_coded_bits > G) {
        rm_bits_concat.resize(G);
        result.n_coded_bits = G;
    }

    result.scrambling_seed = static_cast<uint32_t>(slot_idx + 1);
    BitVec scrambled = scrambler_->scramble(rm_bits_concat, result.scrambling_seed);

    ComplexVec modulated = modulator_->modulate(scrambled, config_.mod_scheme);

    ComplexMat layered = layer_mapper_->map(modulated, n_layers);

    ComplexMat w = get_identity_precoding_matrix();
    ComplexMat precoded = precoder_->precode(layered, w);

    int n_sym = get_symbols_per_slot();
    ResourceGrid tx_grid(n_ports, n_sym, n_sc);
    tx_grid.slot_idx = slot_idx;
    tx_grid.reset();

    dmrs_generator_->generate_dmrs(config_, tx_grid, slot_idx, 0);

    result.dmrs_grid = tx_grid;

    map_pdsch_to_grid(tx_grid, precoded, 0, n_pdsch_rbs_, dmrs_pattern_, config_, slot_idx);

    result.tx_grid = tx_grid;

    ComplexVec tx_signal = ofdm_modulator_->modulate(tx_grid, config_.scs);
    result.tx_signal = tx_signal;

    return result;
}

bool PdschProcessor::decode_transport_block(const SoftVec& descrambled_llr, const PdschTxResult& tx_info,
                                             BitVec& decoded_info_bits) {
    decoded_info_bits.clear();
    int C = tx_info.num_cb;
    int K = tx_info.k_b * tx_info.zc;
    int zc = tx_info.zc;
    int bgn = tx_info.bgn;
    int qm = tx_info.qm;
    int n_punctured = 2 * zc;
    int N = (bgn == 1 ? 68 : 52) * zc;
    int N_cb = N - n_punctured;

    std::vector<BitVec> cb_decoded(C);
    bool all_cb_ok = true;

    for (int c = 0; c < C; c++) {
        int cb_llr_start = tx_info.cb_info[c].offset;
        int E_cb = tx_info.cb_info[c].e_bits;

        SoftVec cb_llr(descrambled_llr.begin() + cb_llr_start,
                       descrambled_llr.begin() + cb_llr_start + E_cb);

        int cb_k_info;
        int cb_crc_len = tx_info.cb_crc_len;

        if (C == 1) {
            cb_k_info = tx_info.n_info_bits + get_crc_length(tx_info.n_info_bits);
        } else {
            cb_k_info = tx_info.cb_size_with_crc;
        }

        int n_filler = K - cb_k_info;
        int filler_start_comp = cb_k_info - n_punctured;
        int n_cb_comp = N_cb - n_filler;

        SoftVec recovered_comp = rate_matcher_->rate_recover(cb_llr, n_cb_comp, 0, bgn, zc, qm, 0);

        SoftVec full_llr(N, 0.0);
        if (n_filler > 0 && filler_start_comp >= 0 && filler_start_comp <= n_cb_comp) {
            for (int i = 0; i < filler_start_comp; i++) {
                full_llr[n_punctured + i] = recovered_comp[i];
            }
            for (int i = 0; i < n_filler; i++) {
                full_llr[cb_k_info + i] = LLR_CLIP;
            }
            for (int i = filler_start_comp; i < n_cb_comp; i++) {
                int dst_idx = cb_k_info + n_filler + (i - filler_start_comp);
                if (dst_idx < N) {
                    full_llr[dst_idx] = recovered_comp[i];
                }
            }
        } else {
            for (int i = 0; i < n_cb_comp && i < N_cb; i++) {
                full_llr[n_punctured + i] = recovered_comp[i];
            }
        }

        for (int i = 0; i < N; i++) {
            if (full_llr[i] > LLR_CLIP) full_llr[i] = LLR_CLIP;
            else if (full_llr[i] < -LLR_CLIP) full_llr[i] = -LLR_CLIP;
        }

        auto [decoded_cb_bits, conv_ok] = ldpc_decoder_->decode(full_llr, bgn, zc,
                                                                 config_.n_ldpc_iterations,
                                                                 config_.early_termination);
        (void)conv_ok;

        if ((int)decoded_cb_bits.size() < cb_k_info) {
            all_cb_ok = false;
            cb_decoded[c] = BitVec(cb_k_info, 0);
            continue;
        }

        BitVec cb_bits(decoded_cb_bits.begin(), decoded_cb_bits.begin() + cb_k_info);

        if (C > 1 && cb_crc_len > 0) {
            auto [cb_info, cb_check_ok] = crc_encoder_->decode(cb_bits, cb_crc_len);
            if (!cb_check_ok) {
                all_cb_ok = false;
            }
            cb_decoded[c] = cb_info;
        } else {
            cb_decoded[c] = cb_bits;
        }
    }

    if (C == 1) {
        int crc_len = get_crc_length(tx_info.n_info_bits);
        int total_len = tx_info.n_info_bits + crc_len;
        if ((int)cb_decoded[0].size() < total_len) {
            return false;
        }
        BitVec rx_bits(cb_decoded[0].begin(), cb_decoded[0].begin() + total_len);
        auto [info_bits, tb_check_ok] = crc_encoder_->decode(rx_bits, crc_len);
        decoded_info_bits = info_bits;
        return tb_check_ok;
    } else {
        BitVec tb_bits_with_crc;
        for (int c = 0; c < C; c++) {
            tb_bits_with_crc.insert(tb_bits_with_crc.end(), cb_decoded[c].begin(), cb_decoded[c].end());
        }
        int tb_crc_len = 24;
        int total_tb_len = tx_info.n_info_bits + tb_crc_len;
        if ((int)tb_bits_with_crc.size() < total_tb_len) {
            return false;
        }
        BitVec rx_bits(tb_bits_with_crc.begin(), tb_bits_with_crc.begin() + total_tb_len);
        auto [info_bits, tb_check_ok] = crc_encoder_->decode(rx_bits, tb_crc_len);
        decoded_info_bits = info_bits;
        return tb_check_ok && all_cb_ok;
    }
}

static ResourceGrid fix_ofdm_demod_grid(const ResourceGrid& demod_grid, int n_ant, int n_sym, int n_sc) {
    int n_sc_rx = static_cast<int>(demod_grid.n_subcarriers);

    if (n_sc_rx == n_sc) {
        ResourceGrid fixed(n_ant, n_sym, n_sc);
        for (int ant = 0; ant < n_ant; ant++) {
            for (int sym = 0; sym < n_sym; sym++) {
                for (int sc = 0; sc < n_sc; sc++) {
                    fixed.set_re(ant, sym, sc, demod_grid.get_re(ant, sym, sc));
                }
            }
        }
        return fixed;
    }

    ResourceGrid fixed(n_ant, n_sym, n_sc);

    int fft_size_rx = 1;
    while (fft_size_rx < n_sc_rx + 1) {
        fft_size_rx <<= 1;
    }
    int fft_size_target = 1;
    while (fft_size_target < n_sc + 1) {
        fft_size_target <<= 1;
    }

    if (fft_size_rx == fft_size_target) {
        int offset = (n_sc_rx - n_sc) / 2;
        for (int ant = 0; ant < n_ant; ant++) {
            for (int sym = 0; sym < n_sym; sym++) {
                for (int sc = 0; sc < n_sc; sc++) {
                    int sc_rx = sc + offset;
                    if (sc_rx >= 0 && sc_rx < n_sc_rx) {
                        fixed.set_re(ant, sym, sc, demod_grid.get_re(ant, sym, sc_rx));
                    } else {
                        fixed.set_re(ant, sym, sc, Complex(0.0, 0.0));
                    }
                }
            }
        }
    } else {
        for (int ant = 0; ant < n_ant; ant++) {
            for (int sym = 0; sym < n_sym; sym++) {
                for (int sc = 0; sc < n_sc; sc++) {
                    fixed.set_re(ant, sym, sc, Complex(0.0, 0.0));
                }
            }
        }
    }

    return fixed;
}

PdschRxResult PdschProcessor::receive(const ResourceGrid& rx_grid_in, const PdschTxResult& tx_info,
                                      double sinr_db, int slot_idx) {
    PdschRxResult result;
    result.crc_ok = false;

    int n_sc = n_pdsch_rbs_ * 12;
    int n_sym = get_symbols_per_slot();
    int n_rx_ant = rx_grid_in.n_ant;
    int n_layers = config_.n_layers;

    ResourceGrid rx_grid = (static_cast<int>(rx_grid_in.n_subcarriers) == n_sc)
                               ? rx_grid_in
                               : fix_ofdm_demod_grid(rx_grid_in, n_rx_ant, n_sym, n_sc);

    ResourceGrid dmrs_grid = tx_info.dmrs_grid;

    SimulationConfig ch_config = config_;
    ch_config.n_layers = n_layers;
    ComplexCube channel_est = channel_estimator_->estimate(rx_grid, dmrs_grid, ch_config);

    double noise_var_ideal = std::pow(10.0, -sinr_db / 10.0);
    double noise_var_est = channel_estimator_->get_estimated_noise_var();
    double noise_var = (noise_var_est > 1e-12) ? noise_var_est : noise_var_ideal;

    ComplexMat rx_pdsch = extract_pdsch_from_grid(rx_grid, 0, n_pdsch_rbs_, dmrs_pattern_, config_, slot_idx);

    ComplexCube h_at_pdsch = extract_channel_at_pdsch(channel_est, n_pdsch_rbs_, n_sc, dmrs_pattern_, config_, slot_idx);

    ComplexMat equalized = equalizer_->equalize(rx_pdsch, h_at_pdsch, noise_var, n_layers);

    ComplexVec demapped = layer_mapper_->demap(equalized, n_layers);

    std::vector<double> eff_noise = equalizer_->get_eff_noise_var();
    int n_pdsch_re = equalized.n_rows;
    std::vector<double> demapped_noise(demapped.n_elem, noise_var);
    for (int i = 0; i < n_pdsch_re; i++) {
        for (int l = 0; l < n_layers; l++) {
            int out_idx = i * n_layers + l;
            int eff_idx = i * n_layers + l;
            if (out_idx < (int)demapped.n_elem && eff_idx < (int)eff_noise.size()) {
                demapped_noise[out_idx] = (eff_noise[eff_idx] > 1e-12) ? eff_noise[eff_idx] : noise_var;
            }
        }
    }

    SoftVec llr = modulator_->demodulate(demapped, config_.mod_scheme, demapped_noise);
    for (int i = 0; i < (int)llr.size(); i++) {
        if (llr[i] > LLR_CLIP) llr[i] = LLR_CLIP;
        else if (llr[i] < -LLR_CLIP) llr[i] = -LLR_CLIP;
    }

    SoftVec descrambled = scrambler_->descramble(llr, tx_info.scrambling_seed);

    BitVec decoded_info;
    result.crc_ok = decode_transport_block(descrambled, tx_info, decoded_info);
    result.decoded_bits = decoded_info;

    result.sinr_est = -10.0 * std::log10(std::max(noise_var, 1e-12));
    result.noise_var_est = noise_var;

    return result;
}

bool PdschProcessor::process_single_snr_point(double sinr_db, BlerResult& result) {
    result.sinr_db = sinr_db;
    result.n_blocks = 0;
    result.n_errors = 0;
    result.bler = 0.0;

    SimulationConfig ch_config = config_;
    ch_config.n_rx_ant = config_.n_rx_ant;
    channel_model_->set_seed(seed_);
    channel_model_->set_config(ch_config);

    int n_sc = n_pdsch_rbs_ * 12;
    int n_sym = get_symbols_per_slot();

    bool fast_awgn = (config_.perfect_csi &&
                      config_.channel_type == ChannelType::AWGN &&
                      config_.n_tx_ant == 1 && config_.n_rx_ant == 1 && config_.n_layers == 1);

    bool fast_fading = (config_.perfect_csi &&
                        (config_.channel_type == ChannelType::TDL_A ||
                         config_.channel_type == ChannelType::TDL_B ||
                         config_.channel_type == ChannelType::TDL_C ||
                         config_.channel_type == ChannelType::TDL_D ||
                         config_.channel_type == ChannelType::TDL_E ||
                         config_.channel_type == ChannelType::CDL_A ||
                         config_.channel_type == ChannelType::CDL_B ||
                         config_.channel_type == ChannelType::CDL_C ||
                         config_.channel_type == ChannelType::CDL_D ||
                         config_.channel_type == ChannelType::CDL_E) &&
                        config_.n_tx_ant == 1 && config_.n_rx_ant == 1 && config_.n_layers == 1);

    double sinr_lin = std::pow(10.0, sinr_db / 10.0);
    double noise_var = 1.0 / sinr_lin;
    double sigma_dim = std::sqrt(noise_var / 2.0);
    std::mt19937 local_rng(seed_);

    for (int block = 0; block < config_.max_blocks_per_sinr; block++) {
        int slot_idx = block;
        TransportBlock tb = generate_transport_block();
        PdschTxResult tx_res = transmit(tb, slot_idx);

        PdschRxResult rx_res;

        if (fast_awgn) {
            ResourceGrid rx_grid = tx_res.tx_grid;
            std::normal_distribution<double> nd(0.0, sigma_dim);
            for (int ant = 0; ant < rx_grid.n_ant; ant++) {
                for (int sym = 0; sym < rx_grid.n_symbols; sym++) {
                    for (int sc = 0; sc < rx_grid.n_subcarriers; sc++) {
                        Complex val = rx_grid.get_re(ant, sym, sc);
                        double nr = nd(local_rng);
                        double ni = nd(local_rng);
                        rx_grid.set_re(ant, sym, sc, val + Complex(nr, ni));
                    }
                }
            }

            rx_res.crc_ok = false;

            ComplexMat rx_pdsch = extract_pdsch_from_grid(rx_grid, 0, n_pdsch_rbs_, dmrs_pattern_, config_, slot_idx);

            ComplexVec rx_symbols(rx_pdsch.n_rows);
            for (int i = 0; i < (int)rx_pdsch.n_rows; i++) {
                rx_symbols(i) = rx_pdsch(i, 0);
            }

            SoftVec llr = modulator_->demodulate(rx_symbols, config_.mod_scheme, noise_var);
            for (int i = 0; i < (int)llr.size(); i++) {
                if (llr[i] > LLR_CLIP) llr[i] = LLR_CLIP;
                else if (llr[i] < -LLR_CLIP) llr[i] = -LLR_CLIP;
            }
            SoftVec descrambled = scrambler_->descramble(llr, tx_res.scrambling_seed);

            BitVec decoded_info;
            rx_res.crc_ok = decode_transport_block(descrambled, tx_res, decoded_info);
            rx_res.decoded_bits = decoded_info;
            rx_res.sinr_est = sinr_db;
            rx_res.noise_var_est = noise_var;
        } else if (fast_fading) {
            int fft_size = get_fft_size(n_pdsch_rbs_, config_.scs);
            double sample_rate = config_.scs * static_cast<double>(fft_size);
            ComplexCube h = channel_model_->get_channel(n_pdsch_rbs_, n_sym, sample_rate);

            std::normal_distribution<double> nd(0.0, sigma_dim);

            int qam_m;
            switch (config_.mod_scheme) {
                case ModulationScheme::QPSK:   qam_m = 1; break;
                case ModulationScheme::QAM16:  qam_m = 2; break;
                case ModulationScheme::QAM256: qam_m = 4; break;
                case ModulationScheme::QAM64:
                default:                       qam_m = 3; break;
            }
            int M_pam = 1 << qam_m;
            int bits_per_sym = qam_m * 2;

            std::vector<Complex> rx_eq;
            std::vector<double> rx_amp;
            rx_eq.reserve(n_sc * n_sym);
            rx_amp.reserve(n_sc * n_sym);

            for (int sym = 0; sym < n_sym; sym++) {
                if (!is_pdsch_symbol(sym, dmrs_pattern_, config_, slot_idx)) continue;
                for (int sc = 0; sc < n_sc; sc++) {
                    Complex tx_val = tx_res.tx_grid.get_re(0, sym, sc);
                    Complex h_val = h(sc, sym, 0);
                    double nr = nd(local_rng);
                    double ni = nd(local_rng);
                    Complex y = tx_val * h_val + Complex(nr, ni);
                    double h_amp = std::abs(h_val);
                    if (h_amp < 1e-9) h_amp = 1e-9;
                    Complex y_mrc = y * std::conj(h_val) / h_amp;
                    rx_eq.push_back(y_mrc);
                    rx_amp.push_back(h_amp);
                }
            }

            int n_pdsch_res = static_cast<int>(rx_eq.size());
            SoftVec llr(n_pdsch_res * bits_per_sym);

            double ms_pam = 0.0;
            for (int k = 0; k < M_pam; k++) { double v = (2.0*k - M_pam + 1.0); ms_pam += v*v; }
            ms_pam /= M_pam;
            double a_pam = 1.0 / std::sqrt(2.0 * ms_pam);
            double var_per_dim = std::max(noise_var / 2.0, 1e-12);

            auto gray_encode = [](int k) { return k ^ (k >> 1); };

            for (int i = 0; i < n_pdsch_res; i++) {
                Complex y = rx_eq[i];
                double amp = rx_amp[i];
                double bits_i[4], bits_q[4];
                double yiq[2] = { y.real(), y.imag() };
                double* llrs_iq[2] = { bits_i, bits_q };

                for (int dim = 0; dim < 2; dim++) {
                    double y_pam = yiq[dim];
                    double* out_l = llrs_iq[dim];
                    for (int b = 0; b < qam_m; b++) {
                        double d0 = 1e30, d1 = 1e30;
                        for (int k = 0; k < M_pam; k++) {
                            double s = (2.0*k - M_pam + 1.0) * a_pam * amp;
                            double dist = (y_pam - s) * (y_pam - s);
                            int g = gray_encode(k);
                            int bit_val = (g >> (qam_m - 1 - b)) & 1;
                            if (bit_val == 0) { if (dist < d0) d0 = dist; }
                            else              { if (dist < d1) d1 = dist; }
                        }
                        out_l[b] = (d1 - d0) / (2.0 * var_per_dim);
                    }
                }
                for (int b = 0; b < qam_m; b++) {
                    llr[i*bits_per_sym + b] = bits_i[b];
                    llr[i*bits_per_sym + qam_m + b] = bits_q[b];
                }
            }

            for (int i = 0; i < (int)llr.size(); i++) {
                if (llr[i] > LLR_CLIP) llr[i] = LLR_CLIP;
                else if (llr[i] < -LLR_CLIP) llr[i] = -LLR_CLIP;
            }

            SoftVec descrambled = scrambler_->descramble(llr, tx_res.scrambling_seed);

            BitVec decoded_info;
            rx_res.crc_ok = decode_transport_block(descrambled, tx_res, decoded_info);
            rx_res.decoded_bits = decoded_info;
            rx_res.sinr_est = sinr_db;
            rx_res.noise_var_est = noise_var;
        } else {
            int fft_size = get_fft_size(n_pdsch_rbs_, config_.scs);
            double sample_rate = config_.scs * static_cast<double>(fft_size);

            ComplexCube h = channel_model_->get_channel(n_pdsch_rbs_, n_sym, sample_rate);
            ComplexVec rx_signal = channel_model_->apply_channel(tx_res.tx_signal, h);

            double noise_var_slow = 1.0 / sinr_lin;
            rx_signal = channel_model_->add_noise(rx_signal, noise_var_slow);

            int n_rx_ant = ch_config.n_rx_ant;
            ResourceGrid demod_grid = ofdm_modulator_->demodulate(rx_signal, n_rx_ant,
                                                                  config_.scs, n_sym);
            ResourceGrid rx_grid = fix_ofdm_demod_grid(demod_grid, n_rx_ant, n_sym, n_sc);
            rx_res = receive(rx_grid, tx_res, sinr_db, block);
        }

        result.n_blocks++;
        if (!rx_res.crc_ok) {
            result.n_errors++;
        }

        result.bler = static_cast<double>(result.n_errors) / result.n_blocks;

        if (result.n_errors >= config_.target_block_errors && result.n_blocks >= 100) {
            break;
        }
        if (result.n_errors == 0 && result.n_blocks >= 1000) {
            break;
        }
    }

    seed_ = local_rng();
    return result.n_errors >= config_.target_block_errors;
}

std::vector<BlerResult> run_bler_simulation(const SimulationConfig& config,
                                             std::unique_ptr<IChannelEstimator> estimator,
                                             const std::string& estimator_name) {
    (void)estimator_name;
    std::vector<BlerResult> results;

    SimulationConfig sim_config = config;
    sim_config.n_tx_ant = std::min(sim_config.n_tx_ant, 2);
    sim_config.n_layers = std::min(sim_config.n_layers, sim_config.n_tx_ant);
    sim_config.n_rx_ant = std::max(1, std::min(sim_config.n_rx_ant, sim_config.n_tx_ant > 1 ? 4 : 1));

    if (sim_config.n_tx_ant < 1) sim_config.n_tx_ant = 1;
    if (sim_config.n_layers < 1) sim_config.n_layers = 1;
    if (sim_config.n_rx_ant < 1) sim_config.n_rx_ant = 1;

    PdschProcessor processor(sim_config);
    if (estimator) {
        processor.set_channel_estimator(std::move(estimator));
    }

    auto channel = channel::create_channel(sim_config.channel_type);
    processor.set_channel_model(std::move(channel));

    for (int i = 0; i < sim_config.n_sinr_points; i++) {
        double sinr_db = sim_config.sinr_start + i * sim_config.sinr_step;
        if (sinr_db > sim_config.sinr_end + 1e-9) {
            break;
        }

        BlerResult res;
        processor.process_single_snr_point(sinr_db, res);
        results.push_back(res);
    }

    return results;
}

}
}
