#include "phy/PdschProcessor.h"
#include "phy/ModuleFactory.h"
#include "common/NrTables.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace nr {
namespace phy {

extern std::unique_ptr<ICrcEncoder> create_crc_encoder();
extern std::unique_ptr<ILdpcEncoder> create_ldpc_encoder();
extern std::unique_ptr<ILdpcDecoder> create_ldpc_decoder();
extern std::unique_ptr<IRateMatcher> create_rate_matcher();
extern std::unique_ptr<IScrambler> create_scrambler();
extern std::unique_ptr<IModulator> create_modulator();
extern std::unique_ptr<ILayerMapper> create_layer_mapper();
extern std::unique_ptr<IPrecoder> create_precoder();
extern std::unique_ptr<IDmrsGenerator> create_dmrs_generator();
extern std::unique_ptr<IResourceMapper> create_resource_mapper();
extern std::unique_ptr<IOfdmModulator> create_ofdm_modulator();
extern std::unique_ptr<IEqualizer> create_mmse_equalizer();
extern std::unique_ptr<IChannelEstimator> create_ls_channel_estimator();

namespace {

constexpr int DMRS_SYMBOL = 2;
constexpr int NUM_SYMBOLS_PER_SLOT = 14;

int get_ldpc_encoded_length(int bgn, int zc) {
    constexpr int BG1_K_B = 22;
    constexpr int BG2_K_B = 10;
    constexpr int N_CHECKS = 15;
    int k_b = (bgn == 1) ? BG1_K_B : BG2_K_B;
    return (k_b + N_CHECKS) * zc;
}

static bool is_dmrs_symbol(int sym) {
    return sym == DMRS_SYMBOL;
}

int count_pdsch_re(int n_rb) {
    int count = 0;
    for (int sym = 0; sym < NUM_SYMBOLS_PER_SLOT; sym++) {
        if (!is_dmrs_symbol(sym)) {
            count += n_rb * 12;
        }
    }
    return count;
}

void map_pdsch_to_grid(ResourceGrid& grid, const ComplexMat& precoded,
                       int rb_start, int rb_len) {
    int n_ant = grid.n_ant;
    int re_idx = 0;
    int total_re = static_cast<int>(precoded.n_rows);

    for (int sym = 0; sym < NUM_SYMBOLS_PER_SLOT; sym++) {
        if (is_dmrs_symbol(sym)) {
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
                                    int rb_start, int rb_len) {
    int n_ant = grid.n_ant;
    int total_re = count_pdsch_re(rb_len);

    ComplexMat extracted(total_re, n_ant, arma::fill::zeros);
    int re_idx = 0;

    for (int sym = 0; sym < NUM_SYMBOLS_PER_SLOT; sym++) {
        if (is_dmrs_symbol(sym)) {
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
                                      int n_sc_total) {
    int n_sc = static_cast<int>(h_est.n_rows);
    int n_sym = static_cast<int>(h_est.n_cols);
    int n_ch = static_cast<int>(h_est.n_slices);
    int total_re = count_pdsch_re(n_rb);

    ComplexCube h_at_pdsch(total_re, 1, n_ch, arma::fill::zeros);
    int re_idx = 0;

    for (int sym = 0; sym < NUM_SYMBOLS_PER_SLOT; sym++) {
        if (is_dmrs_symbol(sym)) {
            continue;
        }
        if (sym >= n_sym) {
            sym = n_sym - 1;
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
                    h_at_pdsch(re_idx, 0, c) = h_est(sc, sym, c);
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
    n_dmrs_symbols_ = 1;

    config_.n_tx_ant = std::min(config_.n_tx_ant, 2);
    config_.n_layers = std::min(config_.n_layers, config_.n_tx_ant);
    config_.n_rx_ant = std::min(config_.n_rx_ant, config_.n_tx_ant > 1 ? 4 : 1);

    pdsch_sym_alloc_.start = 0;
    pdsch_sym_alloc_.length = NUM_SYMBOLS_PER_SLOT;

    init_default_modules();

    n_pdsch_symbols_ = NUM_SYMBOLS_PER_SLOT - n_dmrs_symbols_;
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
    return count_pdsch_re(n_pdsch_rbs_);
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
    int pdsch_re_count = calculate_pdsch_capacity();
    int bits_per_sym = mod_to_bits_per_symbol(config_.mod_scheme);
    int n_layers = config_.n_layers;
    int tbs = static_cast<int>((pdsch_re_count * n_layers * bits_per_sym * config_.code_rate) / 8) * 8;
    tbs = std::max(tbs, 24);
    tbs = std::min(tbs, pdsch_re_count * n_layers * bits_per_sym - 32);

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

    LdpcBaseGraphInfo ldpc_info = select_ldpc_params(tb.tb_size, config_.code_rate);
    result.bgn = ldpc_info.bgn;
    result.zc = ldpc_info.zc;

    BitVec coded_bits = ldpc_encoder_->encode(bits_with_crc, ldpc_info.bgn, ldpc_info.zc);

    int bits_per_sym = mod_to_bits_per_symbol(config_.mod_scheme);
    int pdsch_re_count = calculate_pdsch_capacity();
    int E = pdsch_re_count * n_layers * bits_per_sym;
    result.n_coded_bits = E;

    BitVec rm_bits = rate_matcher_->rate_match(coded_bits, E, tb.rv, ldpc_info.bgn, ldpc_info.zc);

    result.scrambling_seed = static_cast<uint32_t>(slot_idx + 1);
    BitVec scrambled = scrambler_->scramble(rm_bits, result.scrambling_seed);

    ComplexVec modulated = modulator_->modulate(scrambled, config_.mod_scheme);

    ComplexMat layered = layer_mapper_->map(modulated, n_layers);

    ComplexMat w = get_identity_precoding_matrix();
    ComplexMat precoded = precoder_->precode(layered, w);

    int n_sym = get_symbols_per_slot();
    ResourceGrid tx_grid(n_ports, n_sym, n_sc);
    tx_grid.slot_idx = slot_idx;
    tx_grid.reset();

    dmrs_generator_->generate_dmrs(config_, tx_grid, slot_idx, 0);

    map_pdsch_to_grid(tx_grid, precoded, 0, n_pdsch_rbs_);

    result.tx_grid = tx_grid;

    ComplexVec tx_signal = ofdm_modulator_->modulate(tx_grid, config_.scs);
    result.tx_signal = tx_signal;

    return result;
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
                                      double sinr_db, int /*slot_idx*/) {
    PdschRxResult result;
    result.crc_ok = false;

    int n_sc = n_pdsch_rbs_ * 12;
    int n_sym = get_symbols_per_slot();
    int n_rx_ant = rx_grid_in.n_ant;
    int n_layers = config_.n_layers;

    ResourceGrid rx_grid = (static_cast<int>(rx_grid_in.n_subcarriers) == n_sc)
                               ? rx_grid_in
                               : fix_ofdm_demod_grid(rx_grid_in, n_rx_ant, n_sym, n_sc);

    ResourceGrid dmrs_grid = tx_info.tx_grid;

    SimulationConfig ch_config = config_;
    ch_config.n_layers = n_layers;
    ComplexCube channel_est = channel_estimator_->estimate(rx_grid, dmrs_grid, ch_config);

    double noise_var = std::pow(10.0, -sinr_db / 10.0);

    ComplexMat rx_pdsch = extract_pdsch_from_grid(rx_grid, 0, n_pdsch_rbs_);

    ComplexCube h_at_pdsch = extract_channel_at_pdsch(channel_est, n_pdsch_rbs_, n_sc);

    ComplexMat equalized = equalizer_->equalize(rx_pdsch, h_at_pdsch, noise_var, n_layers);

    ComplexVec demapped = layer_mapper_->demap(equalized, n_layers);

    SoftVec llr = modulator_->demodulate(demapped, config_.mod_scheme, noise_var);

    SoftVec descrambled = scrambler_->descramble(llr, tx_info.scrambling_seed);

    int N = get_ldpc_encoded_length(tx_info.bgn, tx_info.zc);
    SoftVec recovered = rate_matcher_->rate_recover(descrambled, N, 0, tx_info.bgn, tx_info.zc);

    auto [decoded_bits, _] = ldpc_decoder_->decode(recovered, tx_info.bgn, tx_info.zc,
                                                      config_.n_ldpc_iterations,
                                                      config_.early_termination);
    (void)_;

    result.crc_ok = false;
    result.decoded_bits = BitVec();

    int crc_len = get_crc_length(tx_info.n_info_bits);
    int total_len = tx_info.n_info_bits + crc_len;

    if (static_cast<int>(decoded_bits.n_elem) >= total_len) {
        BitVec rx_bits(total_len);
        for (int i = 0; i < total_len; i++) {
            rx_bits(i) = decoded_bits(i);
        }
        auto [info_bits, check_ok] = crc_encoder_->decode(rx_bits, crc_len);
        result.decoded_bits = info_bits;
        result.crc_ok = check_ok;
    }

    result.sinr_est = sinr_db;
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
    int fft_size = get_fft_size(n_pdsch_rbs_, config_.scs);
    double sample_rate = config_.scs * static_cast<double>(fft_size);

    for (int block = 0; block < config_.max_blocks_per_sinr; block++) {
        TransportBlock tb = generate_transport_block();

        PdschTxResult tx_res = transmit(tb, block);

        ComplexCube h = channel_model_->get_channel(n_pdsch_rbs_, n_sym, sample_rate);

        ComplexVec rx_signal = channel_model_->apply_channel(tx_res.tx_signal, h);
        double bandwidth = static_cast<double>(n_pdsch_rbs_) * 12.0 * config_.scs;
        double noise_var = channel_model_->get_thermal_noise(sinr_db, bandwidth, 5.0);
        rx_signal = channel_model_->add_noise(rx_signal, noise_var);

        int n_rx_ant = ch_config.n_rx_ant;
        ResourceGrid demod_grid = ofdm_modulator_->demodulate(rx_signal, n_rx_ant,
                                                              config_.scs, n_sym);

        ResourceGrid rx_grid = fix_ofdm_demod_grid(demod_grid, n_rx_ant, n_sym, n_sc);

        PdschRxResult rx_res = receive(rx_grid, tx_res, sinr_db, block);

        result.n_blocks++;
        if (!rx_res.crc_ok) {
            result.n_errors++;
        }

        result.bler = static_cast<double>(result.n_errors) / result.n_blocks;

        if (result.n_errors >= config_.target_block_errors && result.n_blocks >= 10) {
            break;
        }
    }

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

} // namespace phy
} // namespace nr
