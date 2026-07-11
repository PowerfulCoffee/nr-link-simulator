#pragma once

#include "phy/PhyInterfaces.h"
#include <memory>

namespace nr {
namespace phy {

std::unique_ptr<ICrcEncoder> create_crc_encoder();
std::unique_ptr<ILdpcEncoder> create_ldpc_encoder();
std::unique_ptr<ILdpcDecoder> create_ldpc_decoder();
std::unique_ptr<IRateMatcher> create_rate_matcher();
std::unique_ptr<IScrambler> create_scrambler();
std::unique_ptr<IModulator> create_modulator();
std::unique_ptr<ILayerMapper> create_layer_mapper();
std::unique_ptr<IPrecoder> create_precoder();
std::unique_ptr<IDmrsGenerator> create_dmrs_generator();
std::unique_ptr<IResourceMapper> create_resource_mapper();
std::unique_ptr<IOfdmModulator> create_ofdm_modulator();
std::unique_ptr<IChannelEstimator> create_ls_channel_estimator();
std::unique_ptr<IEqualizer> create_mmse_equalizer();

}
}
