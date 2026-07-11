#include "phy/PhyInterfaces.h"
#include <cmath>

namespace nr {
namespace phy {

class MmseEqualizer : public IEqualizer {
public:
    ComplexMat equalize(const ComplexMat& rx_symbols, const ComplexCube& channel_est,
                        double noise_var, int n_layers) override {
        int n_re = rx_symbols.n_rows;
        int n_rx_ant = rx_symbols.n_cols;
        int n_sc = channel_est.n_rows;
        int n_sym = channel_est.n_cols;
        
        ComplexMat equalized(n_re, n_layers, arma::fill::zeros);
        
        for (int i = 0; i < n_re; i++) {
            int sc = i % n_sc;
            int sym_idx = i / n_sc;
            int sym = (sym_idx < n_sym) ? sym_idx : (n_sym - 1);
            
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
                Complex w = std::conj(h) / (h_mag_sq + noise_var);
                equalized(i, 0) = w * y(0, 0);
            } else {
                ComplexMat Hh = H.st();
                ComplexMat HHh = H * Hh;
                for (int rx = 0; rx < n_rx_ant; rx++) {
                    HHh(rx, rx) += noise_var;
                }
                
                ComplexMat W = Hh * HHh.i();
                ComplexMat x_hat = W * y;
                
                for (int l = 0; l < n_layers; l++) {
                    equalized(i, l) = x_hat(l, 0);
                }
            }
        }
        
        return equalized;
    }
};

std::unique_ptr<IEqualizer> create_mmse_equalizer() {
    return std::make_unique<MmseEqualizer>();
}

} // namespace phy
} // namespace nr
