#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <cmath>

namespace nr {
namespace phy {

class MmseEqualizer : public IEqualizer {
public:
    ComplexMat equalize(const ComplexMat& rx_symbols, const ComplexCube& channel_est,
                        double noise_var, int n_layers) override {
        int n_symbols = rx_symbols.n_rows;
        int n_rx_ant = rx_symbols.n_cols;
        int n_sc = channel_est.n_rows;
        
        ComplexMat equalized(n_symbols, n_layers, arma::fill::zeros);
        
        for (int i = 0; i < n_symbols; i++) {
            int sc = i % n_sc;
            
            ComplexMat H(n_rx_ant, n_layers);
            for (int rx = 0; rx < n_rx_ant; rx++) {
                for (int l = 0; l < n_layers; l++) {
                    H(rx, l) = channel_est(sc, 0, rx * n_layers + l);
                }
            }
            
            ComplexMat y(n_rx_ant, 1);
            for (int rx = 0; rx < n_rx_ant; rx++) {
                y(rx, 0) = rx_symbols(i, rx);
            }
            
            ComplexMat Hh = H.t();
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
        
        return equalized;
    }
};

std::unique_ptr<IEqualizer> create_mmse_equalizer() {
    return std::make_unique<MmseEqualizer>();
}

std::unique_ptr<IEqualizer> create_equalizer() {
    return create_mmse_equalizer();
}

} // namespace phy
} // namespace nr
