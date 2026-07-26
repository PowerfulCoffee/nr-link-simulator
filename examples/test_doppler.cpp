#include "common/Types.h"
#include "phy/PhyInterfaces.h"
#include "phy/ModuleFactory.h"
#include "channel/ChannelModels.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace nr;
using namespace nr::phy;
using namespace nr::channel;

std::vector<int> find_dmrs_positions(const ResourceGrid& grid, int n_tx_ant) {
    std::vector<int> pos;
    for (int sym = 0; sym < grid.n_symbols; sym++) {
        bool has = false;
        for (int sc = 0; sc < grid.n_subcarriers && !has; sc++) {
            for (int tx = 0; tx < n_tx_ant && !has; tx++) {
                if (std::abs(grid.get_re(tx, sym, sc)) > 1e-10) has = true;
            }
        }
        if (has) pos.push_back(sym);
    }
    return pos;
}

int main() {
    std::cout << "=== Doppler Shift Estimation Test ===" << std::endl;
    std::cout << std::endl;

    SimulationConfig config;
    config.scs = 30e3;
    config.n_rb = 25;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.channel_type = ChannelType::TDL_D;
    config.delay_spread = 30e-9;
    config.max_doppler = 100.0;
    config.enable_los = true;
    config.dmrs_additional_pos = 1;

    const int n_sc = config.n_rb * 12;
    const int n_sym = 14;
    const double sample_rate = config.scs * 2048;

    auto channel = create_channel(config.channel_type);
    channel->set_config(config);
    channel->set_seed(42);

    auto dmrs_gen = create_dmrs_generator();
    auto chan_est = create_ls_channel_estimator();

    ResourceGrid tx_grid(config.n_tx_ant, n_sym, n_sc);
    dmrs_gen->generate_dmrs(config, tx_grid, 0, 0);
    
    std::vector<int> dmrs_pos = find_dmrs_positions(tx_grid, config.n_tx_ant);
    std::cout << "DMRS symbol positions (dmrs_additional_pos=" << config.dmrs_additional_pos << "):";
    for (int p : dmrs_pos) std::cout << " " << p;
    std::cout << std::endl;
    
    int fft_size = 2048;
    auto get_ofdm_time = [&](int sym_idx) {
        double tu = 1.0 / config.scs;
        double t = 0.0;
        for (int s = 0; s < sym_idx; s++) {
            bool is_first = (s == 0);
            int cp = 144 * fft_size / 2048;
            if (is_first) cp += 16 * fft_size / 2048;
            t += tu * (fft_size + cp) / fft_size;
        }
        bool is_first = (sym_idx == 0);
        int cp = 144 * fft_size / 2048;
        if (is_first) cp += 16 * fft_size / 2048;
        t += tu * (fft_size + cp/2) / fft_size;
        return t;
    };
    
    if (dmrs_pos.size() >= 2) {
        double t0 = get_ofdm_time(dmrs_pos[0]);
        double t1 = get_ofdm_time(dmrs_pos[1]);
        double delta_t = t1 - t0;
        std::cout << "Time between DMRS symbols: " << delta_t*1e6 << " us (" << delta_t*1e3 << " ms)" << std::endl;
        std::cout << "Max unambiguous Doppler: ±" << 1.0/(2*delta_t) << " Hz" << std::endl;
    }
    std::cout << std::endl;

    std::vector<double> test_dopplers = {0.0, 50.0, 100.0, 200.0, 300.0};

    std::cout << "--- TDL-D (LoS, K-factor=13.3dB) Channel Tests ---" << std::endl;
    for (double fd_test : test_dopplers) {
        config.max_doppler = fd_test;
        channel->set_config(config);
        channel->set_seed(42);

        ComplexCube h = channel->get_channel(config.n_rb, n_sym, sample_rate);

        ResourceGrid rx_grid(config.n_rx_ant, n_sym, n_sc);
        for (int sym = 0; sym < n_sym; sym++) {
            for (int sc = 0; sc < n_sc; sc++) {
                Complex tx_val = tx_grid.get_re(0, sym, sc);
                Complex h_val = h(sc, sym, 0);
                rx_grid.set_re(0, sym, sc, tx_val * h_val);
            }
        }

        ComplexCube h_est = chan_est->estimate(rx_grid, tx_grid, config);
        double est_fd = chan_est->get_estimated_doppler();

        std::cout << "True fd: " << std::setw(8) << fd_test << " Hz | "
                  << "Est. fd: " << std::setw(10) << std::fixed << std::setprecision(2) << est_fd << " Hz";
        if (fd_test > 0) {
            double err = std::abs(est_fd - fd_test);
            std::cout << " | Abs Error: " << std::setw(8) << err << " Hz";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "--- Verify with noiseless true channel (fd=200Hz) ---" << std::endl;
    config.max_doppler = 200.0;
    channel->set_config(config);
    channel->set_seed(123);
    ComplexCube h = channel->get_channel(config.n_rb, n_sym, sample_rate);

    Complex cross(0,0);
    for (int sc = 0; sc < n_sc; sc++) {
        Complex h0 = h(sc, dmrs_pos[0], 0);
        Complex h1 = h(sc, dmrs_pos[1], 0);
        cross += h1 * std::conj(h0);
    }
    double delta_t = get_ofdm_time(dmrs_pos[1]) - get_ofdm_time(dmrs_pos[0]);
    double phase = std::arg(cross);
    double est_fd_true = phase / (2*M_PI*delta_t);
    std::cout << "True H cross-corr phase = " << phase << " rad" << std::endl;
    std::cout << "Doppler from true H: " << est_fd_true << " Hz (expected ~200 Hz)" << std::endl;

    return 0;
}
