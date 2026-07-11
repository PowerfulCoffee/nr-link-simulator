#include "phy/PdschProcessor.h"
#include "common/NrTables.h"
#include <cmath>
#include <limits>

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
extern std::unique_ptr<IEqualizer> create_equalizer();
extern std::unique_ptr<IChannelEstimator> create_ls_channel_estimator();

PdschProcessor::PdschProcessor(const SimulationConfig& config)
    : config_(config), seed_(config.random_seed), rng_(config.random_seed) {
    init_default_modules();
    
    n_pdsch_rbs_ = config_.n_rb;
    n_dmrs_symbols_ = config_.dmrs_duration;
    pdsch_sym_alloc_.start = 0;
    pdsch_sym_alloc_.length = 14;
    n_pdsch_symbols_ = calculate_pdsch_capacity();
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
    equalizer_ = create_equalizer();
    channel_estimator_ = create_ls_channel_estimator();
    channel_model_ = channel::create_channel(config_.channel_type);
    channel_model_->set_config(config_);
}

int PdschProcessor::calculate_pdsch_capacity() {
    DmrsPattern pattern = get_dmrs_pattern(config_.dmrs_type, config_.dmrs_additional_pos,
                                           config_.dmrs_duration);
    int n_pdsch_sym = 0;
    int dmrs_re_per_prb = 0;
    for (int sym = 0; sym < 14; sym++) {
        if (pattern.re_per_prb[sym] > 0) {
            dmrs_re_per_prb += pattern.re_per_prb[sym];
        } else {
            n_pdsch_sym++;
        }
    }
    return n_pdsch_sym;
}

ComplexMat PdschProcessor::get_identity_precoding_matrix() {
    int n_layers = config_.n_layers;
    int n_ant = std::min(config_.n_tx_ant, 2);
    ComplexMat w(n_ant, n_layers, arma::fill::zeros);
    for (int l = 0; l < std::min(n_layers, n_ant); l++) {
        w(l, l) = Complex(1.0, 0.0);
    }
    return w;
}

TransportBlock PdschProcessor::generate_transport_block() {
    int tbs = get_tbs(config_.mcs_index, n_pdsch_rbs_, n_pdsch_symbols_, config_.n_layers);
    TransportBlock tb(tbs);
    tb.mcs = config_.mcs_index;
    tb.rv = 0;
    tb.generate_random_bits(seed_++);
    return tb;
}

PdschTxResult PdschProcessor::transmit(const TransportBlock& tb, int slot_idx) {
    PdschTxResult result;
    result.n_info_bits = tb.tb_size;
    
    int crc_len = get_crc_length(tb.tb_size);
    BitVec bits_with_crc = crc_encoder_->encode(tb.bits, crc_len);
    result.tb_bits_after_crc = bits_with_crc;
    
    LdpcBaseGraphInfo ldpc_info = select_ldpc_params(tb.tb_size, config_.code_rate);
    result.bgn = ldpc_info.bgn;
    result.zc = ldpc_info.zc;
    
    BitVec coded_bits = ldpc_encoder_->encode(bits_with_crc, ldpc_info.bgn, ldpc_info.zc);
    
    int bits_per_sym = mod_to_bits_per_symbol(config_.mod_scheme);
    int total_re = n_pdsch_symbols_ * n_pdsch_rbs_ * 12;
    int n_coded_symbols = total_re * config_.n_layers;
    int E = n_coded_symbols * bits_per_sym;
    result.n_coded_bits = E;
    
    BitVec rm_bits = rate_matcher_->rate_match(coded_bits, E, tb.rv, ldpc_info.bgn, ldpc_info.zc);
    
    result.scrambling_seed = static_cast<uint32_t>((slot_idx << 16) | (config_.n_rb & 0xFFFF));
    BitVec scrambled = scrambler_->scramble(rm_bits, result.scrambling_seed);
    
    ComplexVec modulated = modulator_->modulate(scrambled, config_.mod_scheme);
    
    ComplexMat layered = layer_mapper_->map(modulated, config_.n_layers);
    
    ComplexMat w = get_identity_precoding_matrix();
    ComplexMat precoded = precoder_->precode(layered, w);
    
    int n_sc = n_pdsch_rbs_ * 12;
    int n_sym = get_symbols_per_slot();
    int n_tx_ant = std::min(config_.n_tx_ant, 2);
    ResourceGrid tx_grid(n_tx_ant, n_sym, n_sc);
    tx_grid.slot_idx = slot_idx;
    tx_grid.reset();
    
    dmrs_generator_->generate_dmrs(config_, tx_grid, slot_idx, 2);
    
    resource_mapper_->map_pdsch(tx_grid, precoded, 0, n_pdsch_rbs_, pdsch_sym_alloc_);
    result.tx_grid = tx_grid;
    
    ComplexVec tx_signal = ofdm_modulator_->modulate(tx_grid, config_.scs);
    result.tx_signal = tx_signal;
    
    return result;
}

PdschRxResult PdschProcessor::receive(const ResourceGrid& rx_grid, const PdschTxResult& tx_info,
                                      double sinr_db, int slot_idx) {
    PdschRxResult result;
    
    ResourceGrid dmrs_grid = tx_info.tx_grid;
    
    ComplexCube channel_est = channel_estimator_->estimate(rx_grid, dmrs_grid, config_);
    
    double noise_var = std::pow(10.0, -sinr_db / 10.0);
    
    ComplexMat rx_pdsch = resource_mapper_->extract_pdsch(rx_grid, 0, n_pdsch_rbs_, pdsch_sym_alloc_);
    
    int n_rx_ant = rx_grid.n_ant;
    int n_layers = config_.n_layers;
    int n_pdsch_re = rx_pdsch.n_rows;
    
    ComplexMat sym_for_eq(n_pdsch_re, n_rx_ant);
    for (int i = 0; i < n_pdsch_re; i++) {
        for (int rx = 0; rx < n_rx_ant; rx++) {
            sym_for_eq(i, rx) = rx_pdsch(i, rx);
        }
    }
    
    ComplexMat equalized = equalizer_->equalize(sym_for_eq, channel_est, noise_var, n_layers);
    
    ComplexVec demapped = layer_mapper_->demap(equalized, n_layers);
    
    SoftVec llr = modulator_->demodulate(demapped, config_.mod_scheme, noise_var);
    
    SoftVec descrambled = scrambler_->descramble(llr, tx_info.scrambling_seed);
    
    int N_cb;
    if (tx_info.bgn == 1) {
        N_cb = 66 * tx_info.zc;
    } else {
        N_cb = 50 * tx_info.zc;
    }
    SoftVec recovered = rate_matcher_->rate_recover(descrambled, N_cb, 0, tx_info.bgn, tx_info.zc);
    
    auto [decoded_bits, crc_ok] = ldpc_decoder_->decode(recovered, tx_info.bgn, tx_info.zc,
                                                        config_.n_ldpc_iterations,
                                                        config_.early_termination);
    
    if (crc_ok) {
        int crc_len = get_crc_length(tx_info.n_info_bits);
        if ((int)decoded_bits.n_elem >= tx_info.n_info_bits + crc_len) {
            BitVec rx_bits(tx_info.n_info_bits + crc_len);
            for (int i = 0; i < tx_info.n_info_bits + crc_len && i < (int)decoded_bits.n_elem; i++) {
                rx_bits(i) = decoded_bits(i);
            }
            auto [info_bits, check_ok] = crc_encoder_->decode(rx_bits, crc_len);
            result.decoded_bits = info_bits;
            result.crc_ok = check_ok;
        } else {
            result.decoded_bits = BitVec();
            result.crc_ok = false;
        }
    } else {
        result.decoded_bits = BitVec();
        result.crc_ok = false;
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
    
    channel_model_->set_seed(seed_);
    channel_model_->set_config(config_);
    
    double bandwidth = config_.n_rb * 12 * config_.scs;
    double noise_var = channel_model_->get_thermal_noise(sinr_db, bandwidth, 5.0);
    
    for (int block = 0; block < config_.max_blocks_per_sinr; block++) {
        TransportBlock tb = generate_transport_block();
        
        PdschTxResult tx_res = transmit(tb, 0);
        
        ComplexCube h = channel_model_->get_channel(config_.n_rb, get_symbols_per_slot(),
                                                    config_.scs * 4096);
        
        ComplexVec rx_signal = channel_model_->apply_channel(tx_res.tx_signal, h);
        rx_signal = channel_model_->add_noise(rx_signal, noise_var);
        
        ResourceGrid rx_grid = ofdm_modulator_->demodulate(rx_signal, config_.n_rx_ant,
                                                           config_.scs, get_symbols_per_slot());
        
        PdschRxResult rx_res = receive(rx_grid, tx_res, sinr_db, 0);
        
        result.n_blocks++;
        if (!rx_res.crc_ok) {
            result.n_errors++;
        }
        
        result.bler = static_cast<double>(result.n_errors) / result.n_blocks;
        
        if (result.n_errors >= config_.target_block_errors && result.n_blocks >= 100) {
            break;
        }
    }
    
    return result.n_errors >= config_.target_block_errors;
}

std::vector<BlerResult> run_bler_simulation(const SimulationConfig& config,
                                             std::unique_ptr<IChannelEstimator> estimator,
                                             const std::string& estimator_name) {
    std::vector<BlerResult> results;
    
    PdschProcessor processor(config);
    if (estimator) {
        processor.set_channel_estimator(std::move(estimator));
    }
    
    auto channel = channel::create_channel(config.channel_type);
    processor.set_channel_model(std::move(channel));
    
    for (int i = 0; i < config.n_sinr_points; i++) {
        double sinr_db = config.sinr_start + i * config.sinr_step;
        if (sinr_db > config.sinr_end + 1e-9) break;
        
        BlerResult res;
        processor.process_single_snr_point(sinr_db, res);
        results.push_back(res);
    }
    
    return results;
}

} // namespace phy
} // namespace nr
