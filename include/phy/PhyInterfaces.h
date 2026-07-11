#pragma once

#include "common/Types.h"

namespace nr {
namespace phy {

class ICrcEncoder {
public:
    virtual ~ICrcEncoder() = default;
    virtual BitVec encode(const BitVec& bits, int crc_poly) = 0;
    virtual std::pair<BitVec, bool> decode(const BitVec& bits, int crc_poly) = 0;
};

class ILdpcEncoder {
public:
    virtual ~ILdpcEncoder() = default;
    virtual BitVec encode(const BitVec& info_bits, int bgn, int zc) = 0;
};

class ILdpcDecoder {
public:
    virtual ~ILdpcDecoder() = default;
    virtual std::pair<BitVec, bool> decode(const SoftVec& llr, int bgn, int zc,
                                           int n_iter, bool early_term) = 0;
};

class IRateMatcher {
public:
    virtual ~IRateMatcher() = default;
    virtual BitVec rate_match(const BitVec& coded_bits, int E, int rv, int bgn, int zc) = 0;
    virtual SoftVec rate_recover(const SoftVec& llr, int N, int rv, int bgn, int zc) = 0;
};

class IScrambler {
public:
    virtual ~IScrambler() = default;
    virtual BitVec scramble(const BitVec& bits, uint32_t cinit) = 0;
    virtual SoftVec descramble(const SoftVec& llr, uint32_t cinit) = 0;
};

class IModulator {
public:
    virtual ~IModulator() = default;
    virtual ComplexVec modulate(const BitVec& bits, ModulationScheme scheme) = 0;
    virtual SoftVec demodulate(const ComplexVec& symbols, ModulationScheme scheme, double noise_var) = 0;
};

class ILayerMapper {
public:
    virtual ~ILayerMapper() = default;
    virtual ComplexMat map(const ComplexVec& symbols, int n_layers) = 0;
    virtual ComplexVec demap(const ComplexMat& layered_symbols, int n_layers) = 0;
};

class IPrecoder {
public:
    virtual ~IPrecoder() = default;
    virtual ComplexMat precode(const ComplexMat& layered_symbols, const ComplexMat& w) = 0;
    virtual ComplexMat deprecode(const ComplexMat& rx_symbols, const ComplexMat& channel_est, int n_layers) = 0;
};

class IDmrsGenerator {
public:
    virtual ~IDmrsGenerator() = default;
    virtual void generate_dmrs(const SimulationConfig& config, ResourceGrid& grid,
                               int slot_idx, int symbol_start) = 0;
    virtual ComplexCube extract_dmrs(const ResourceGrid& rx_grid, const SimulationConfig& config,
                                     int slot_idx, int symbol_start) = 0;
};

class IChannelEstimator {
public:
    virtual ~IChannelEstimator() = default;
    virtual ComplexCube estimate(const ResourceGrid& rx_grid, const ResourceGrid& dmrs_grid,
                                 const SimulationConfig& config) = 0;
    virtual std::string get_name() const = 0;
};

class IEqualizer {
public:
    virtual ~IEqualizer() = default;
    virtual ComplexMat equalize(const ComplexMat& rx_symbols, const ComplexCube& channel_est,
                                double noise_var, int n_layers) = 0;
};

class IResourceMapper {
public:
    virtual ~IResourceMapper() = default;
    virtual void map_pdsch(ResourceGrid& grid, const ComplexMat& precoded_symbols,
                           int rb_start, int rb_len, const SymbolAllocation& sym_alloc) = 0;
    virtual ComplexMat extract_pdsch(const ResourceGrid& grid, int rb_start, int rb_len,
                                     const SymbolAllocation& sym_alloc) = 0;
};

class IOfdmModulator {
public:
    virtual ~IOfdmModulator() = default;
    virtual ComplexVec modulate(const ResourceGrid& grid, double scs) = 0;
    virtual ResourceGrid demodulate(const ComplexVec& signal, int n_ant, double scs, int n_symbols) = 0;
};

} // namespace phy
} // namespace nr
