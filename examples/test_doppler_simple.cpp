#include "common/Types.h"
#include "phy/PhyInterfaces.h"
#include "phy/ModuleFactory.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace nr;
using namespace nr::phy;

int main() {
    std::cout << "=== Simple Doppler Estimation Verification ===" << std::endl;
    
    SimulationConfig config;
    config.scs = 30e3;
    config.n_rb = 25;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.n_layers = 1;
    config.dmrs_type = DmrsType::TYPE1;
    config.dmrs_additional_pos = 1;
    config.dmrs_duration = 1;
    
    const int n_sc = config.n_rb * 12;
    const int n_sym = 14;
    
    auto dmrs_gen = create_dmrs_generator();
    auto chan_est = create_ls_channel_estimator();
    
    ResourceGrid tx_grid(config.n_tx_ant, n_sym, n_sc);
    dmrs_gen->generate_dmrs(config, tx_grid, 0, 0);
    
    std::vector<int> dmrs_pos;
    for (int sym = 0; sym < n_sym; sym++) {
        bool has = false;
        for (int sc = 0; sc < n_sc && !has; sc++) {
            if (std::abs(tx_grid.get_re(0, sym, sc)) > 1e-10) has = true;
        }
        if (has) dmrs_pos.push_back(sym);
    }
    std::cout << "DMRS positions:";
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
    
    double t0 = get_ofdm_time(dmrs_pos[0]);
    double t1 = get_ofdm_time(dmrs_pos[1]);
    double delta_t = t1 - t0;
    std::cout << "delta_t = " << delta_t*1e6 << " us" << std::endl;
    
    double test_fd = 200.0;
    std::cout << "\n--- Test with known Doppler: fd = " << test_fd << " Hz ---" << std::endl;
    
    ResourceGrid rx_grid(config.n_rx_ant, n_sym, n_sc);
    Complex h_const(0.707, 0.707);
    
    for (int sym = 0; sym < n_sym; sym++) {
        double t_sym = get_ofdm_time(sym);
        double phase = 2.0 * M_PI * test_fd * t_sym;
        Complex doppler_rot(std::cos(phase), std::sin(phase));
        Complex h_now = h_const * doppler_rot;
        
        for (int sc = 0; sc < n_sc; sc++) {
            Complex tx_val = tx_grid.get_re(0, sym, sc);
            rx_grid.set_re(0, sym, sc, tx_val * h_now);
        }
    }
    
    ComplexCube h_est = chan_est->estimate(rx_grid, tx_grid, config);
    double est_fd = chan_est->get_estimated_doppler();
    std::cout << "Estimated Doppler: " << est_fd << " Hz (expected: " << test_fd << " Hz)" << std::endl;
    std::cout << "Estimation error: " << std::abs(est_fd - test_fd) << " Hz" << std::endl;
    
    std::cout << "\nCheck h_est at DMRS positions:" << std::endl;
    for (int s : dmrs_pos) {
        Complex h_check = h_est(n_sc/2, s, 0);
        double t_sym = get_ofdm_time(s);
        double expected_phase = 2.0*M_PI*test_fd*t_sym + std::arg(h_const);
        std::cout << "  Symbol " << s << ": h_est = " << h_check 
                  << ", angle = " << std::arg(h_check) 
                  << " rad (expected " << expected_phase << " rad)" << std::endl;
    }
    
    std::cout << "\n--- Test with multiple Doppler values ---" << std::endl;
    for (double fd : {0.0, 50.0, 100.0, 200.0, 500.0}) {
        for (int sym = 0; sym < n_sym; sym++) {
            double t_sym = get_ofdm_time(sym);
            double phase = 2.0 * M_PI * fd * t_sym;
            Complex doppler_rot(std::cos(phase), std::sin(phase));
            Complex h_now = h_const * doppler_rot;
            for (int sc = 0; sc < n_sc; sc++) {
                Complex tx_val = tx_grid.get_re(0, sym, sc);
                rx_grid.set_re(0, sym, sc, tx_val * h_now);
            }
        }
        h_est = chan_est->estimate(rx_grid, tx_grid, config);
        est_fd = chan_est->get_estimated_doppler();
        std::cout << "True fd = " << std::setw(8) << fd << " Hz | Est fd = " 
                  << std::setw(10) << std::fixed << std::setprecision(2) << est_fd << " Hz | Error = " 
                  << std::setw(8) << std::abs(est_fd - fd) << " Hz" << std::endl;
    }
    
    return 0;
}
