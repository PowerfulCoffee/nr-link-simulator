#include "phy/PhyInterfaces.h"
#include <armadillo>
#include <cmath>

namespace nr {
namespace phy {

namespace {

// DFT beam vector for the Type I single-panel codebook (TS 38.214 5.2.1.2):
//   v_l = exp( j*2*pi*( m1*(l mod N1)/N1 + m2*floor(l/N1)/N2 ) ),  l = 0..N1*N2-1
ComplexVec type1_beam(int n1, int n2, int m1, int m2) {
    int n_elem = n1 * n2;
    ComplexVec v(n_elem);
    for (int l = 0; l < n_elem; l++) {
        int l1 = l % n1;
        int l2 = l / n1;
        double phase = 2.0 * M_PI * (static_cast<double>(m1) * l1 / n1
                                     + static_cast<double>(m2) * l2 / std::max(n2, 1));
        v[l] = Complex(std::cos(phase), std::sin(phase));
    }
    return v;
}

} // anonymous namespace

ComplexMat build_type1_csi_precoder(int n_tx_ant, int n_layers,
                                    int n1, int n2, int i1_phi,
                                    int i11, int i12, int i2) {
    // TS 38.214 5.2.1.2 Type I single-panel:
    //   W = W1 * W2,  W1 = blkdiag(V, alpha*V),  alpha = e^{j*pi*phi/2}
    //   rank 1: V = [v_m1,m2],              W2 = 1
    //   rank 2: V = [v_m, v_{m+pi/2}],      W2 = (1/sqrt2)*[[1, 1], [1, -1]]
    // rank > 2 falls back to identity (single-layer-per-port), flagged by
    // returning an identity-like matrix sized n_tx_ant x n_layers.
    int n_elem = n1 * n2;
    int n_ports_expected = 2 * n_elem;

    ComplexMat w(n_tx_ant, n_layers, arma::fill::zeros);

    if (n_layers > 2 || n_tx_ant != n_ports_expected) {
        // Fallback: identity mapping (layer l -> port l)
        for (int l = 0; l < n_layers && l < n_tx_ant; l++) {
            w(l, l) = Complex(1.0, 0.0);
        }
        return w;
    }

    double alpha = std::cos(M_PI * i1_phi / 2.0);
    double alpha_im = std::sin(M_PI * i1_phi / 2.0);
    Complex pol_phase(alpha, alpha_im);

    // Beam 1 from (i11, i12); beam 2 is the adjacent DFT beam (pi/2-shifted pair)
    int m1 = ((i11 % n1) + n1) % n1;
    int m2 = ((i12 % n2) + n2) % n2;
    int m1_b = (m1 + 1) % n1;

    ComplexVec beam1 = type1_beam(n1, n2, m1, m2);
    ComplexVec beam2 = type1_beam(n1, n2, m1_b, m2);

    auto put_beam = [&](const ComplexVec& beam, int col, double sign) {
        // W1 columns: [beam; alpha*beam], then W2 Hadamard folding is merged
        // analytically: col = (b + alpha*b)/sqrt2 or (b - alpha*b)/sqrt2
        for (int p = 0; p < n_elem; p++) {
            w(p, col) = beam[p] / std::sqrt(2.0);
        }
        for (int p = 0; p < n_elem; p++) {
            w(n_elem + p, col) = beam[p] * pol_phase * sign / std::sqrt(2.0);
        }
    };

    if (n_layers == 1) {
        // 2Tx single-pol codebook is covered by n_elem == 1 with alpha phase;
        // i2 extends the polarization phase to 4 codewords for n_elem == 1.
        Complex phase = pol_phase;
        if (n_elem == 1) {
            double ph = M_PI * i2 / 2.0;
            phase = Complex(std::cos(ph), std::sin(ph));
        }
        for (int p = 0; p < n_elem; p++) {
            w(p, 0) = beam1[p] / std::sqrt(2.0);
        }
        for (int p = 0; p < n_elem; p++) {
            w(n_elem + p, 0) = beam1[p] * phase / std::sqrt(2.0);
        }
    } else { // n_layers == 2
        put_beam(beam1, 0, +1.0);
        put_beam(beam2, 1, -1.0);
    }

    // Normalize each column to unit power
    for (int l = 0; l < n_layers; l++) {
        double norm = 0.0;
        for (int p = 0; p < n_tx_ant; p++) {
            norm += std::norm(w(p, l));
        }
        if (norm > 1e-12) {
            double s = 1.0 / std::sqrt(norm);
            for (int p = 0; p < n_tx_ant; p++) {
                w(p, l) *= s;
            }
        }
    }

    return w;
}

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
