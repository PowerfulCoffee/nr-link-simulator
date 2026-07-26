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
    arma::cx_mat X = arma::solve(arma::trimatu(L.ht()), Y, arma::solve_opts::fast);
    return ComplexMat(X.memptr(), X.n_rows, X.n_cols);
}

inline ComplexMat invert_via_cholesky(const ComplexMat& A) {
    arma::cx_mat L;
    bool success = arma::chol(L, A, "lower");
    if (!success) {
        ComplexMat A_reg = A;
        for (arma::uword i = 0; i < A.n_rows; i++) {
            A_reg(i, i) += Complex(1e-8, 0.0);
        }
        success = arma::chol(L, A_reg, "lower");
        if (!success) {
            return A.i();
        }
    }
    arma::cx_mat L_inv = arma::inv(arma::trimatl(L));
    return L_inv.ht() * L_inv;
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

        double sqrt_no = std::sqrt(std::max(noise_var, 1e-12));

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

            ComplexMat H_tilde = H / sqrt_no;
            ComplexMat y_tilde = y / sqrt_no;

            if (n_layers <= n_rx_ant) {
                ComplexMat Hh = H_tilde.ht();
                ComplexMat HhH = Hh * H_tilde;
                ComplexMat A = HhH;
                for (int l = 0; l < n_layers; l++) {
                    A(l, l) += Complex(1.0, 0.0);
                }

                ComplexMat Hhy = Hh * y_tilde;
                ComplexMat x_hat = solve_via_cholesky(A, Hhy);

                ComplexMat A_inv = invert_via_cholesky(A);

                for (int l = 0; l < n_layers; l++) {
                    double err_var = std::real(A_inv(l, l));
                    err_var = std::max(err_var, 1e-12);
                    double d = 1.0 - err_var;
                    if (std::abs(d) < 1e-12) d = (d >= 0) ? 1e-12 : -1e-12;
                    equalized(i, l) = x_hat(l, 0) / d;
                    eff_noise_var_[i * n_layers + l] = err_var / d;
                }
            } else {
                ComplexMat Hh = H_tilde.ht();
                ComplexMat HHh = H_tilde * Hh;
                for (int rx = 0; rx < n_rx_ant; rx++) {
                    HHh(rx, rx) += Complex(1.0, 0.0);
                }

                ComplexMat HHh_inv = invert_via_cholesky(HHh);
                ComplexMat G = Hh * HHh_inv;
                ComplexMat x_hat = G * y_tilde;

                ComplexMat eye_nl(n_layers, n_layers, arma::fill::eye);
                ComplexMat GH = eye_nl - G * H_tilde;

                for (int l = 0; l < n_layers; l++) {
                    double d = GH(l, l).real();
                    if (std::abs(d) < 1e-12) d = (d >= 0) ? 1e-12 : -1e-12;
                    equalized(i, l) = x_hat(l, 0) / d;
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
