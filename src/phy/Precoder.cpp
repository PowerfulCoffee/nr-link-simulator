#include "phy/PhyInterfaces.h"
#include <armadillo>
#include <cmath>

namespace nr {
namespace phy {

class PrecoderImpl : public IPrecoder {
public:
    ComplexMat precode(const ComplexMat& layered_symbols, const ComplexMat& w) override {
        int n_symbols = static_cast<int>(layered_symbols.n_rows);
        int n_layers = static_cast<int>(layered_symbols.n_cols);
        int n_ant = static_cast<int>(w.n_rows);

        ComplexMat precoded(n_symbols, n_ant, arma::fill::zeros);

        for (int i = 0; i < n_symbols; i++) {
            for (int ant = 0; ant < n_ant; ant++) {
                Complex sum(0, 0);
                for (int l = 0; l < n_layers; l++) {
                    sum += w(ant, l) * layered_symbols(i, l);
                }
                precoded(i, ant) = sum;
            }
        }

        return precoded;
    }

    ComplexMat deprecode(const ComplexMat& rx_symbols, const ComplexMat& channel_est, int n_layers) override {
        int n_symbols = static_cast<int>(rx_symbols.n_rows);
        int n_rx_ant = static_cast<int>(rx_symbols.n_cols);

        ComplexMat deprecoded(n_symbols, n_layers, arma::fill::zeros);

        for (int i = 0; i < n_symbols; i++) {
            for (int l = 0; l < n_layers; l++) {
                if (l < n_rx_ant) {
                    deprecoded(i, l) = rx_symbols(i, l);
                } else {
                    deprecoded(i, l) = rx_symbols(i, 0);
                }
            }
        }
        (void)channel_est;

        return deprecoded;
    }
};

std::unique_ptr<IPrecoder> create_precoder() {
    return std::make_unique<PrecoderImpl>();
}

}
}
