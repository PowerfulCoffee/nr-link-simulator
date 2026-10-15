/*
 * Diagnostic: 4x4 AWGN + LS channel estimation + identity precoder.
 * Prints intermediate quantities (DMRS grid layout, LS estimate diagonal
 * vs off-diagonal, effective noise, CRC result) to localize the MIMO-LS
 * failure observed via the web runner.
 */
#include "common/Types.h"
#include "common/NrTables.h"
#include "channel/ChannelModels.h"
#include "phy/PdschProcessor.h"
#include "phy/ModuleFactory.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace nr;
using namespace nr::phy;

static void dump_dmrs_grid(const ResourceGrid& g) {
    std::cout << "[DMRS] grid n_ant=" << g.n_ant << " n_sym=" << g.n_symbols
              << " n_sc=" << g.n_subcarriers << "\n";
    for (int ant = 0; ant < g.n_ant; ant++) {
        int count = 0;
        int kmod4_hist[4] = {0, 0, 0, 0};
        Complex first = g.get_re(ant, 2, 0);
        bool have_first = false;
        for (int sym = 0; sym < g.n_symbols; sym++) {
            for (int sc = 0; sc < g.n_subcarriers; sc++) {
                if (std::abs(g.get_re(ant, sym, sc)) > 1e-10) {
                    count++;
                    if (sym == 2) kmod4_hist[sc % 4]++;
                    if (!have_first) { first = g.get_re(ant, sym, sc); have_first = true; }
                }
            }
        }
        std::cout << "  ant" << ant << ": nonzero=" << count
                  << " sym2 k%4={" << kmod4_hist[0] << "," << kmod4_hist[1]
                  << "," << kmod4_hist[2] << "," << kmod4_hist[3] << "}"
                  << " first=" << std::fixed << std::setprecision(3)
                  << first.real() << (first.imag() >= 0 ? "+" : "") << first.imag() << "j\n";
    }
}

static void dump_h_est(const ComplexCube& h, int n_rx, int n_layers) {
    int n_sc = (int)h.n_rows, n_sym = (int)h.n_cols, n_ch = (int)h.n_slices;
    std::cout << "[H_EST] cube (n_sc=" << n_sc << ", n_sym=" << n_sym
              << ", n_ch=" << n_ch << ")\n";
    int sym = 2;
    for (int rx = 0; rx < n_rx; rx++) {
        for (int l = 0; l < n_layers; l++) {
            int ch = rx * n_layers + l;
            // sample a few subcarriers of this channel slice
            Complex h00 = h(10, sym, ch), h01 = h(11, sym, ch);
            std::cout << "  h[rx" << rx << ",l" << l << "] sc10="
                      << std::setprecision(3) << h00.real()
                      << (h00.imag() >= 0 ? "+" : "") << h00.imag() << "j"
                      << " sc11=" << h01.real()
                      << (h01.imag() >= 0 ? "+" : "") << h01.imag() << "j\n";
        }
    }
}

int main() {
    SimulationConfig cfg;
    cfg.mcs_index = 9;
    cfg.mcs_table = 1;
    cfg.mod_scheme = mcs_to_modulation(9, 1);
    cfg.code_rate = mcs_to_code_rate(9, 1);
    cfg.n_rb = 10;
    cfg.n_tx_ant = 4;
    cfg.n_rx_ant = 4;
    cfg.n_layers = 4;
    cfg.channel_type = ChannelType::AWGN;
    cfg.perfect_csi = false;
    cfg.precoder_mode = PrecoderMode::IDENTITY;
    cfg.dmrs_type = DmrsType::TYPE1;
    cfg.dmrs_additional_pos = 2;
    cfg.dmrs_duration = 1;
    cfg.n_ldpc_iterations = 25;
    cfg.early_termination = true;
    cfg.max_blocks_per_sinr = 1;
    cfg.target_block_errors = 1;
    cfg.random_seed = 123;
    cfg.scs = 15000;

    PdschProcessor proc(cfg);
    proc.set_seed(123);

    TransportBlock tb = proc.generate_transport_block();
    std::cout << "[TB] size=" << tb.tb_size << "\n";
    PdschTxResult tx = proc.transmit(tb, 0);
    dump_dmrs_grid(tx.dmrs_grid);
    std::cout << "[TX] tx_signal len=" << tx.tx_signal.n_elem << "\n";

    // AWGN channel: apply + noise
    auto ch = channel::create_channel(ChannelType::AWGN);
    ch->set_config(cfg);
    int n_sc = cfg.n_rb * 12;
    int fft_size = 1;
    while (fft_size < n_sc + 1) fft_size <<= 1;
    double sample_rate = cfg.scs * (double)fft_size;
    ComplexCube h_true = ch->get_channel(cfg.n_rb, 14, sample_rate);
    ComplexVec rx_sig = ch->apply_channel(tx.tx_signal, h_true);
    double sinr_db = 12.0;
    double noise_var = 1.0 / std::pow(10.0, sinr_db / 10.0);
    rx_sig = ch->add_noise(rx_sig, noise_var);
    std::cout << "[CH] noise_var=" << std::scientific << noise_var << std::fixed << "\n";

    // OFDM demodulate
    auto ofdm = create_ofdm_modulator();
    int n_sym = 14;
    ResourceGrid demod = ofdm->demodulate(rx_sig, cfg.n_rx_ant, cfg.scs, n_sym);
    std::cout << "[DEMOD] n_ant=" << demod.n_ant << " n_sym=" << demod.n_symbols
              << " n_sc=" << demod.n_subcarriers << "\n";

    // LS estimate
    auto est = create_ls_channel_estimator();
    ComplexCube h_est = est->estimate(demod, tx.dmrs_grid, cfg);
    std::cout << "[NOISE_EST] " << std::scientific
              << est->get_estimated_noise_var() << std::fixed << "\n";
    dump_h_est(h_est, cfg.n_rx_ant, cfg.n_layers);

    // full RX
    PdschRxResult rx = proc.receive(demod, tx, sinr_db, 0);
    std::cout << "[RX] crc_ok=" << rx.crc_ok
              << " sinr_est=" << std::setprecision(2) << rx.sinr_est << "\n";
    return 0;
}
