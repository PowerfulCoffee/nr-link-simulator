#pragma once

#include <complex>
#include <vector>
#include <memory>
#include <cstdint>
#include <cstddef>
#include <armadillo>

namespace nr {

using Complex = std::complex<double>;
using Bit = uint8_t;
using SoftBit = double;

template <typename T>
using Vec = arma::Col<T>;

template <typename T>
using Mat = arma::Mat<T>;

template <typename T>
using Cube = arma::Cube<T>;

using ComplexVec = Vec<Complex>;
using ComplexMat = Mat<Complex>;
using ComplexCube = Cube<Complex>;
using BitVec = std::vector<Bit>;
using SoftVec = std::vector<SoftBit>;

enum class ChannelType {
    AWGN,
    TDL_A,
    TDL_B,
    TDL_C,
    TDL_D,
    TDL_E,
    CDL_A,
    CDL_B,
    CDL_C,
    CDL_D,
    CDL_E
};

enum class ModulationScheme {
    BPSK,
    QPSK,
    QAM16,
    QAM64,
    QAM256
};

enum class DmrsType {
    TYPE1,
    TYPE2
};

enum class SlotType {
    DOWNLINK,
    UPLINK,
    SPECIAL
};

struct SymbolAllocation {
    int start = 0;
    int length = 14;
};

struct TddConfig {
    static constexpr int slots_per_frame = 20;
    static constexpr int symbols_per_slot = 14;
    
    SlotType slot_pattern[20] = {
        SlotType::DOWNLINK, SlotType::DOWNLINK, SlotType::DOWNLINK,
        SlotType::SPECIAL, SlotType::UPLINK, SlotType::UPLINK,
        SlotType::DOWNLINK, SlotType::DOWNLINK, SlotType::DOWNLINK, SlotType::DOWNLINK,
        SlotType::DOWNLINK, SlotType::DOWNLINK, SlotType::DOWNLINK,
        SlotType::SPECIAL, SlotType::UPLINK, SlotType::UPLINK,
        SlotType::DOWNLINK, SlotType::DOWNLINK, SlotType::DOWNLINK, SlotType::DOWNLINK
    };
    
    int s_slot_dl_symbols = 6;
    int s_slot_gp_symbols = 4;
    int s_slot_ul_symbols = 4;
};

struct SimulationConfig {
    double scs = 30e3;
    int n_rb = 273;
    int n_tx_ant = 64;
    int n_rx_ant = 4;
    int n_layers = 4;
    
    ChannelType channel_type = ChannelType::AWGN;
    double delay_spread = 30e-9;
    double max_doppler = 5.0;
    bool enable_los = false;
    
    DmrsType dmrs_type = DmrsType::TYPE1;
    int dmrs_additional_pos = 0;
    int dmrs_duration = 1;
    bool tdd_enabled = false;
    
    int mcs_index = 15;
    ModulationScheme mod_scheme = ModulationScheme::QPSK;
    double code_rate = 0.3;
    
    TddConfig tdd_config;
    
    int n_sinr_points = 21;
    double sinr_start = -2.0;
    double sinr_end = 18.0;
    double sinr_step = 1.0;
    
    int max_blocks_per_sinr = 10000;
    int target_block_errors = 100;
    uint64_t random_seed = 12345;
    
    int n_ldpc_iterations = 25;
    bool early_termination = true;
    bool perfect_csi = false;
};

struct ResourceGrid {
    int n_ant = 0;
    int n_symbols = 0;
    int n_subcarriers = 0;
    int slot_idx = 0;
    ComplexCube data;
    
    ResourceGrid() = default;
    ResourceGrid(int ant, int sym, int sc)
        : n_ant(ant), n_symbols(sym), n_subcarriers(sc),
          data(sc, sym, ant, arma::fill::zeros) {}
    
    void set_re(int ant, int sym, int sc, Complex val) {
        data(sc, sym, ant) = val;
    }
    
    Complex get_re(int ant, int sym, int sc) const {
        return data(sc, sym, ant);
    }
    
    void reset() {
        data.zeros();
    }
};

struct TransportBlock {
    int tb_size = 0;
    int mcs = 0;
    int rv = 0;
    BitVec bits;
    bool crc_ok = false;
    
    TransportBlock() = default;
    TransportBlock(int size) : tb_size(size), bits(size, 0) {}
    
    void generate_random_bits(uint64_t seed);
};

} // namespace nr
