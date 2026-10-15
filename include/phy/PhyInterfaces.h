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
    virtual BitVec rate_match(const BitVec& coded_bits, int E, int rv, int bgn, int zc,
                              int qm = 0, int n_filler = 0) = 0;
    virtual SoftVec rate_recover(const SoftVec& llr, int N, int rv, int bgn, int zc,
                                  int qm = 0, int n_filler = 0) = 0;
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
    virtual SoftVec demodulate(const ComplexVec& symbols, ModulationScheme scheme,
                               const std::vector<double>& noise_vars) = 0;
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
    virtual double get_estimated_noise_var() const { return 0.0; }
    virtual double get_estimated_doppler() const { return 0.0; }
    virtual void set_perfect_channel(const ComplexCube& /*h*/, double /*noise_var*/) {}
};

class IEqualizer {
public:
    virtual ~IEqualizer() = default;
    virtual ComplexMat equalize(const ComplexMat& rx_symbols, const ComplexCube& channel_est,
                                double noise_var, int n_layers) = 0;
    virtual std::vector<double> get_eff_noise_var() const { return {}; }
};

class IResourceMapper {
public:
    virtual ~IResourceMapper() = default;
    // dmrs_symbols: OFDM symbol indices carrying DMRS; PDSCH data must not be
    // mapped on them (caller passes empty for non-DMRS-aware legacy use)
    virtual void map_pdsch(ResourceGrid& grid, const ComplexMat& precoded_symbols,
                           int rb_start, int rb_len, const SymbolAllocation& sym_alloc,
                           const std::vector<int>& dmrs_symbols = {}) = 0;
    virtual ComplexMat extract_pdsch(const ResourceGrid& grid, int rb_start, int rb_len,
                                     const SymbolAllocation& sym_alloc,
                                     const std::vector<int>& dmrs_symbols = {}) = 0;
};

// TS 38.214 5.2.1.2 Type I single-panel CSI codebook precoding matrix:
// W (n_tx_ant x n_layers), unit-norm columns; rank 1/2 exact, rank > 2 falls
// back to identity. n_tx_ant must equal 2*N1*N2 (dual polarization).
ComplexMat build_type1_csi_precoder(int n_tx_ant, int n_layers,
                                    int n1, int n2, int i1_phi,
                                    int i11, int i12, int i2);

class IOfdmModulator {
public:
    virtual ~IOfdmModulator() = default;
    virtual ComplexVec modulate(const ResourceGrid& grid, double scs) = 0;
    virtual ResourceGrid demodulate(const ComplexVec& signal, int n_ant, double scs, int n_symbols) = 0;
};

} // namespace phy
} // namespace nr
