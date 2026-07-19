#include "phy/PhyInterfaces.h"
#include <cmath>
#include <vector>

namespace nr {
namespace phy {

class MmseEqualizer : public IEqualizer {
public:
    ComplexMat equalize(const ComplexMat& rx_symbols, const ComplexCube& channel_est,
                        double noise_var, int n_layers) override {
        int n_re = rx_symbols.n_rows;
        int n_rx_ant = rx_symbols.n_cols;
        int n_ch_sc = channel_est.n_rows;
        int n_ch_sym = channel_est.n_cols;

        ComplexMat equalized(n_re, n_layers, arma::fill::zeros);
        eff_noise_var_.resize(n_re * n_layers);

        bool flat_channel = (n_ch_sym == 1);

        for (int i = 0; i < n_re; i++) {
            int sc, sym;
            if (flat_channel) {
                sc = i;
                sym = 0;
            } else {
                sc = i % n_ch_sc;
                int sym_idx = i / n_ch_sc;
                sym = (sym_idx < n_ch_sym) ? sym_idx : (n_ch_sym - 1);
            }

            ComplexMat H(n_rx_ant, n_layers);
            for (int rx = 0; rx < n_rx_ant; rx++) {
                for (int l = 0; l < n_layers; l++) {
                    int ch_idx = rx * n_layers + l;
                    if (ch_idx < (int)channel_est.n_slices) {
                        H(rx, l) = channel_est(sc, sym, ch_idx);
                    } else {
                        H(rx, l) = (rx == l) ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
                    }
                }
            }

            ComplexMat y(n_rx_ant, 1);
            for (int rx = 0; rx < n_rx_ant; rx++) {
                y(rx, 0) = rx_symbols(i, rx);
            }

            if (n_rx_ant == 1 && n_layers == 1) {
                Complex h = H(0, 0);
                double h_mag_sq = std::norm(h);

                Complex gy = std::conj(h) * y(0, 0) / (h_mag_sq + noise_var);
                double d = h_mag_sq / (h_mag_sq + noise_var);
                equalized(i, 0) = gy / d;
                eff_noise_var_[i] = (1.0 / d - 1.0);
            } else {
                ComplexMat Hh = H.st();
                ComplexMat HHh = H * Hh;
                for (int rx = 0; rx < n_rx_ant; rx++) {
                    HHh(rx, rx) += noise_var;
                }

                ComplexMat G = Hh * HHh.i();
                ComplexMat gy = G * y;
                ComplexMat GH = G * H;

                for (int l = 0; l < n_layers; l++) {
                    double d = GH(l, l).real();
                    if (std::abs(d) < 1e-12) d = 1e-12;
                    equalized(i, l) = gy(l, 0) / d;
                    eff_noise_var_[i * n_layers + l] = (1.0 / d - 1.0);
                }
            }
        }

        return equalized;
    }

    std::vector<double> get_eff_noise_var() const override {
        return eff_noise_var_;
    }

private:
    std::vector<double> eff_noise_var_;
};

std::unique_ptr<IEqualizer> create_mmse_equalizer() {
    return std::make_unique<MmseEqualizer>();
}

} // namespace phy
} // namespace nr
