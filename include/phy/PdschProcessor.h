#pragma once

#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PhyInterfaces.h"
#include "channel/ChannelModels.h"
#include <memory>
#include <random>

namespace nr {
namespace phy {

struct PdschTxResult {
    ResourceGrid tx_grid;
    ComplexVec tx_signal;
    BitVec tb_bits_after_crc;
    int n_info_bits;
    int n_coded_bits;
    int bgn;
    int zc;
    uint32_t scrambling_seed;
};

struct PdschRxResult {
    BitVec decoded_bits;
    bool crc_ok;
    double sinr_est;
    double noise_var_est;
};

class PdschProcessor {
public:
    PdschProcessor(const SimulationConfig& config);
    ~PdschProcessor() = default;
    
    void set_channel_estimator(std::unique_ptr<IChannelEstimator> estimator);
    void set_channel_model(std::unique_ptr<channel::IChannelModel> channel);
    
    void set_seed(uint64_t seed);
    
    TransportBlock generate_transport_block();
    
    PdschTxResult transmit(const TransportBlock& tb, int slot_idx);
    PdschRxResult receive(const ResourceGrid& rx_grid, const PdschTxResult& tx_info,
                          double sinr_db, int slot_idx);
    
    bool process_single_snr_point(double sinr_db, BlerResult& result);
    
private:
    SimulationConfig config_;
    uint64_t seed_;
    std::mt19937 rng_;
    
    std::unique_ptr<ICrcEncoder> crc_encoder_;
    std::unique_ptr<ILdpcEncoder> ldpc_encoder_;
    std::unique_ptr<ILdpcDecoder> ldpc_decoder_;
    std::unique_ptr<IRateMatcher> rate_matcher_;
    std::unique_ptr<IScrambler> scrambler_;
    std::unique_ptr<IModulator> modulator_;
    std::unique_ptr<ILayerMapper> layer_mapper_;
    std::unique_ptr<IPrecoder> precoder_;
    std::unique_ptr<IDmrsGenerator> dmrs_generator_;
    std::unique_ptr<IResourceMapper> resource_mapper_;
    std::unique_ptr<IOfdmModulator> ofdm_modulator_;
    std::unique_ptr<IEqualizer> equalizer_;
    std::unique_ptr<IChannelEstimator> channel_estimator_;
    std::unique_ptr<channel::IChannelModel> channel_model_;
    
    SymbolAllocation pdsch_sym_alloc_;
    int n_pdsch_symbols_;
    int n_pdsch_rbs_;
    int n_dmrs_symbols_;
    
    void init_default_modules();
    int calculate_pdsch_capacity();
    ComplexMat get_identity_precoding_matrix();
};

std::vector<BlerResult> run_bler_simulation(const SimulationConfig& config,
                                             std::unique_ptr<IChannelEstimator> estimator,
                                             const std::string& estimator_name = "LS");

} // namespace phy
} // namespace nr
