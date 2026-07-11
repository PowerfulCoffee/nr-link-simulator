#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <cmath>

namespace nr {
namespace phy {

class PrecoderImpl : public IPrecoder {
public:
    ComplexMat precode(const ComplexMat& layered_symbols, const ComplexMat& w) override {
        int n_symbols = layered_symbols.n_rows;
        int n_layers = layered_symbols.n_cols;
        int n_ant = w.n_rows;
        
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
    
    ComplexMat deprecode(const ComplexMat& rx_symbols, const ComplexMat& w, int n_layers) override {
        int n_symbols = rx_symbols.n_rows;
        int n_ant = rx_symbols.n_cols;
        
        ComplexMat wh = w.t();
        for (int ant = 0; ant < w.n_rows; ant++) {
            for (int l = 0; l < w.n_cols; l++) {
                wh(l, ant) = std::conj(w(ant, l));
            }
        }
        
        ComplexMat deprecoded(n_symbols, n_layers, arma::fill::zeros);
        for (int i = 0; i < n_symbols; i++) {
            for (int l = 0; l < n_layers; l++) {
                Complex sum(0, 0);
                for (int ant = 0; ant < n_ant; ant++) {
                    sum += wh(l, ant) * rx_symbols(i, ant);
                }
                deprecoded(i, l) = sum;
            }
        }
        
        return deprecoded;
    }
};

std::unique_ptr<IPrecoder> create_precoder() {
    return std::make_unique<PrecoderImpl>();
}

} // namespace phy
} // namespace nr
