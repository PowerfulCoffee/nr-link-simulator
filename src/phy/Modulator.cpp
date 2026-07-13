#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <cmath>
#include <vector>
#include <complex>
#include <armadillo>

namespace nr {
namespace phy {

class ModulatorImpl : public IModulator {
public:
    ComplexVec modulate(const BitVec& bits, ModulationScheme scheme) override {
        int bits_per_sym = mod_to_bits_per_symbol(scheme);
        int n_symbols = static_cast<int>(bits.size()) / bits_per_sym;
        ComplexVec symbols(n_symbols);
        
        for (int i = 0; i < n_symbols; i++) {
            if (scheme == ModulationScheme::BPSK) {
                double b = bits[i] ? -1.0 : 1.0;
                symbols(i) = Complex(b, 0.0);
            } else if (scheme == ModulationScheme::QPSK) {
                double inv_sqrt2 = 1.0 / std::sqrt(2.0);
                double re = (bits[2*i] == 0) ? inv_sqrt2 : -inv_sqrt2;
                double im = (bits[2*i+1] == 0) ? inv_sqrt2 : -inv_sqrt2;
                symbols(i) = Complex(re, im);
            } else if (scheme == ModulationScheme::QAM16) {
                double inv_sqrt10 = 1.0 / std::sqrt(10.0);
                double re = ((bits[4*i] ? -1.0 : 1.0) * (1.0 + 2.0 * bits[4*i+1])) * inv_sqrt10;
                double im = ((bits[4*i+2] ? -1.0 : 1.0) * (1.0 + 2.0 * bits[4*i+3])) * inv_sqrt10;
                symbols(i) = Complex(re, im);
            } else if (scheme == ModulationScheme::QAM64) {
                double inv_sqrt42 = 1.0 / std::sqrt(42.0);
                double re = ((bits[6*i] ? -1.0 : 1.0) * (1.0 + 2.0*bits[6*i+1] + 4.0*(bits[6*i+2]^bits[6*i+1]))) * inv_sqrt42;
                double im = ((bits[6*i+3] ? -1.0 : 1.0) * (1.0 + 2.0*bits[6*i+4] + 4.0*(bits[6*i+5]^bits[6*i+4]))) * inv_sqrt42;
                symbols(i) = Complex(re, im);
            } else {
                double inv_sqrt2 = 1.0 / std::sqrt(2.0);
                double re = (bits[2*i] == 0) ? inv_sqrt2 : -inv_sqrt2;
                double im = (bits[2*i+1] == 0) ? inv_sqrt2 : -inv_sqrt2;
                symbols(i) = Complex(re, im);
            }
        }
        
        return symbols;
    }
    
    SoftVec demodulate(const ComplexVec& symbols, ModulationScheme scheme, double noise_var) override {
        int bits_per_sym = mod_to_bits_per_symbol(scheme);
        int n_symbols = static_cast<int>(symbols.n_elem);
        SoftVec llr(n_symbols * bits_per_sym);
        
        double sig2 = std::max(noise_var, 1e-10);
        
        for (int i = 0; i < n_symbols; i++) {
            if (scheme == ModulationScheme::BPSK) {
                llr[i] = 2.0 * symbols(i).real() / sig2;
            } else if (scheme == ModulationScheme::QPSK) {
                double re = symbols(i).real();
                double im = symbols(i).imag();
                llr[2*i]   = 2.0 * std::sqrt(2.0) * re / sig2;
                llr[2*i+1] = 2.0 * std::sqrt(2.0) * im / sig2;
            } else if (scheme == ModulationScheme::QAM16) {
                double inv_sqrt10 = 1.0 / std::sqrt(10.0);
                double re = symbols(i).real();
                double im = symbols(i).imag();
                llr[4*i]   = 4.0 * re / (sig2 * std::sqrt(10.0));
                llr[4*i+1] = 4.0 * (2.0 * inv_sqrt10 - std::fabs(re)) / sig2;
                llr[4*i+2] = 4.0 * im / (sig2 * std::sqrt(10.0));
                llr[4*i+3] = 4.0 * (2.0 * inv_sqrt10 - std::fabs(im)) / sig2;
            } else if (scheme == ModulationScheme::QAM64) {
                double inv_sqrt42 = 1.0 / std::sqrt(42.0);
                double re = symbols(i).real();
                double im = symbols(i).imag();
                llr[6*i]   = 4.0 * re / (sig2 * std::sqrt(42.0));
                llr[6*i+1] = 4.0 * (2.0 * inv_sqrt42 - std::fabs(re)) / sig2;
                llr[6*i+2] = 4.0 * (2.0 * inv_sqrt42 - std::fabs(std::fabs(re) - 4.0 * inv_sqrt42)) / sig2;
                llr[6*i+3] = 4.0 * im / (sig2 * std::sqrt(42.0));
                llr[6*i+4] = 4.0 * (2.0 * inv_sqrt42 - std::fabs(im)) / sig2;
                llr[6*i+5] = 4.0 * (2.0 * inv_sqrt42 - std::fabs(std::fabs(im) - 4.0 * inv_sqrt42)) / sig2;
            } else {
                llr[2*i]   = 2.0 * std::sqrt(2.0) * symbols(i).real() / sig2;
                llr[2*i+1] = 2.0 * std::sqrt(2.0) * symbols(i).imag() / sig2;
            }
        }
        
        return llr;
    }
};

std::unique_ptr<IModulator> create_modulator() {
    return std::make_unique<ModulatorImpl>();
}

}
}
