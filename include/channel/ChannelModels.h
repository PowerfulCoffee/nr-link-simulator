#pragma once

#include "common/Types.h"
#include <random>

namespace nr {
namespace channel {

class IChannelModel {
public:
    virtual ~IChannelModel() = default;
    
    virtual void set_seed(uint64_t seed) = 0;
    virtual void set_config(const SimulationConfig& config) = 0;
    
    virtual ComplexCube get_channel(int n_prbs, int n_symbols, double sample_rate) = 0;
    virtual ComplexVec apply_channel(const ComplexVec& tx_signal, const ComplexCube& h) = 0;
    virtual double get_thermal_noise(double sinr_db, double bandwidth, double noise_figure) = 0;
    virtual ComplexVec add_noise(const ComplexVec& signal, double noise_var) = 0;
    virtual std::string get_name() const = 0;
};

struct TdlParams {
    double norm_delay[24];
    double norm_pwr[24];
    int n_taps;
    bool los;
    double k_factor_dB;
};

struct CdlParams {
    double norm_delay[24];
    double norm_pwr[24];
    Complex aod[24];
    Complex aoa[24];
    Complex zod[24];
    Complex zoa[24];
    double xpr[24];
    int n_clusters;
    bool los;
    double k_factor_dB;
};

class AwgnChannel : public IChannelModel {
public:
    AwgnChannel();
    void set_seed(uint64_t seed) override;
    void set_config(const SimulationConfig& config) override;
    ComplexCube get_channel(int n_prbs, int n_symbols, double sample_rate) override;
    ComplexVec apply_channel(const ComplexVec& tx_signal, const ComplexCube& h) override;
    double get_thermal_noise(double sinr_db, double bandwidth, double noise_figure) override;
    ComplexVec add_noise(const ComplexVec& signal, double noise_var) override;
    std::string get_name() const override { return "AWGN"; }
    
private:
    SimulationConfig config_;
    uint64_t seed_;
};

class TdlChannel : public IChannelModel {
public:
    TdlChannel(ChannelType type);
    void set_seed(uint64_t seed) override;
    void set_config(const SimulationConfig& config) override;
    ComplexCube get_channel(int n_prbs, int n_symbols, double sample_rate) override;
    ComplexVec apply_channel(const ComplexVec& tx_signal, const ComplexCube& h) override;
    double get_thermal_noise(double sinr_db, double bandwidth, double noise_figure) override;
    ComplexVec add_noise(const ComplexVec& signal, double noise_var) override;
    std::string get_name() const override;
    
private:
    ChannelType type_;
    SimulationConfig config_;
    uint64_t seed_;
    std::mt19937 rng_;
    TdlParams tdl_params_;
    
    void init_tdl_params();
};

class CdlChannel : public IChannelModel {
public:
    CdlChannel(ChannelType type);
    void set_seed(uint64_t seed) override;
    void set_config(const SimulationConfig& config) override;
    ComplexCube get_channel(int n_prbs, int n_symbols, double sample_rate) override;
    ComplexVec apply_channel(const ComplexVec& tx_signal, const ComplexCube& h) override;
    double get_thermal_noise(double sinr_db, double bandwidth, double noise_figure) override;
    ComplexVec add_noise(const ComplexVec& signal, double noise_var) override;
    std::string get_name() const override;
    
private:
    ChannelType type_;
    SimulationConfig config_;
    uint64_t seed_;
    std::mt19937 rng_;
    CdlParams cdl_params_;
    
    void init_cdl_params();
    ComplexMat generate_ula_response(int n_ant, double angle_spacing);
};

std::unique_ptr<IChannelModel> create_channel(ChannelType type);

} // namespace channel
} // namespace nr
