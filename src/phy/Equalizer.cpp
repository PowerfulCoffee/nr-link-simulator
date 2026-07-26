#include "phy/PhyInterfaces.h"
#include <cmath>
#include <vector>
#include <stdexcept>

namespace nr {
namespace phy {

namespace {

inline ComplexMat solve_via_cholesky(const ComplexMat& A, const ComplexMat& B) {
    arma::cx_mat L;
    bool success = arma::chol(L, A, "lower");
    if (!success) {
        ComplexMat A_reg = A;
        for (arma::uword i = 0; i < A.n_rows; i++) {
            A_reg(i, i) += Complex(1e-8, 0.0);
        }
        success = arma::chol(L, A_reg, "lower");
        if (!success) {
            return A.i() * B;
        }
    }
    arma::cx_mat Y = arma::solve(arma::trimatl(L), B, arma::solve_opts::fast);
    arma::cx_mat X = arma::solve(arma::trimatu(L.t()), Y, arma::solve_opts::fast);
    return ComplexMat(X.memptr(), X.n_rows, X.n_cols);
}

}

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
            } else if (n_layers <= n_rx_ant) {
                ComplexMat Hh = H.st();
                ComplexMat HhH = Hh * H;
                ComplexMat A = HhH;
                for (int l = 0; l < n_layers; l++) {
                    A(l, l) += Complex(noise_var, 0.0);
                }

                ComplexMat Hhy = Hh * y;
                ComplexMat gy = solve_via_cholesky(A, Hhy);

                ComplexMat eye_nl(n_layers, n_layers, arma::fill::eye);
                ComplexMat GH = eye_nl - noise_var * solve_via_cholesky(A, eye_nl);

                for (int l = 0; l < n_layers; l++) {
                    double d = GH(l, l).real();
                    if (std::abs(d) < 1e-12) d = 1e-12;
                    equalized(i, l) = gy(l, 0) / d;
                    eff_noise_var_[i * n_layers + l] = (1.0 / d - 1.0);
                }
            } else {
                ComplexMat Hh = H.st();
                ComplexMat HHh = H * Hh;
                for (int rx = 0; rx < n_rx_ant; rx++) {
                    HHh(rx, rx) += Complex(noise_var, 0.0);
                }

                ComplexMat G = Hh * solve_via_cholesky(HHh, ComplexMat(n_rx_ant, n_rx_ant, arma::fill::eye));
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
