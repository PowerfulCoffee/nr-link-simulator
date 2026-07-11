#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <cmath>

namespace nr {
namespace phy {

class ModulatorImpl : public IModulator {
public:
    ComplexVec modulate(const BitVec& bits, ModulationScheme scheme) override {
        int bits_per_sym = mod_to_bits_per_symbol(scheme);
        int n_symbols = bits.n_elem / bits_per_sym;
        ComplexVec symbols(n_symbols);
        
        for (int i = 0; i < n_symbols; i++) {
            BitVec sym_bits = bits(arma::span(i * bits_per_sym, (i + 1) * bits_per_sym - 1));
            symbols(i) = modulate_symbol(sym_bits, scheme);
        }
        
        return symbols;
    }
    
    SoftVec demodulate(const ComplexVec& symbols, ModulationScheme scheme, double noise_var) override {
        int bits_per_sym = mod_to_bits_per_symbol(scheme);
        int n_symbols = symbols.n_elem;
        SoftVec llr(n_symbols * bits_per_sym);
        
        for (int i = 0; i < n_symbols; i++) {
            SoftVec sym_llr = demodulate_symbol(symbols(i), scheme, noise_var);
            for (int j = 0; j < bits_per_sym; j++) {
                llr(i * bits_per_sym + j) = sym_llr(j);
            }
        }
        
        return llr;
    }

private:
    Complex modulate_symbol(const BitVec& bits, ModulationScheme scheme) {
        switch (scheme) {
            case ModulationScheme::QPSK:
                return modulate_qpsk(bits);
            case ModulationScheme::QAM16:
                return modulate_qam16(bits);
            case ModulationScheme::QAM64:
                return modulate_qam64(bits);
            case ModulationScheme::QAM256:
                return modulate_qam256(bits);
            default:
                return modulate_qpsk(bits);
        }
    }
    
    Complex modulate_qpsk(const BitVec& bits) {
        double inv_sqrt2 = 1.0 / std::sqrt(2.0);
        double re = (bits(0) == 0) ? inv_sqrt2 : -inv_sqrt2;
        double im = (bits(1) == 0) ? inv_sqrt2 : -inv_sqrt2;
        return Complex(re, im);
    }
    
    Complex modulate_qam16(const BitVec& bits) {
        double inv_sqrt10 = 1.0 / std::sqrt(10.0);
        double levels[4] = {-3.0, -1.0, 1.0, 3.0};
        
        int re_idx = bits(0) * 2 + bits(1);
        int im_idx = bits(2) * 2 + bits(3);
        
        double re = levels[re_idx ^ (bits(0) ? 0 : 0)] * inv_sqrt10;
        double im = levels[im_idx ^ (bits(2) ? 0 : 0)] * inv_sqrt10;
        
        re = ((bits(0) ? -1.0 : 1.0) * (1.0 + 2.0 * bits(1))) * inv_sqrt10;
        im = ((bits(2) ? -1.0 : 1.0) * (1.0 + 2.0 * bits(3))) * inv_sqrt10;
        
        return Complex(re, im);
    }
    
    Complex modulate_qam64(const BitVec& bits) {
        double inv_sqrt42 = 1.0 / std::sqrt(42.0);
        
        double re = ((bits(0) ? -1.0 : 1.0) * (1.0 + 2.0 * bits(1) + 4.0 * (bits(2) ^ bits(1)))) * inv_sqrt42;
        double im = ((bits(3) ? -1.0 : 1.0) * (1.0 + 2.0 * bits(4) + 4.0 * (bits(5) ^ bits(4)))) * inv_sqrt42;
        
        return Complex(re, im);
    }
    
    Complex modulate_qam256(const BitVec& bits) {
        double inv_sqrt170 = 1.0 / std::sqrt(170.0);
        
        double re = ((bits(0) ? -1.0 : 1.0) * (1.0 + 2.0 * bits(1) + 4.0 * (bits(2) ^ bits(1)) + 8.0 * (bits(3) ^ bits(2) ^ bits(1)))) * inv_sqrt170;
        double im = ((bits(4) ? -1.0 : 1.0) * (1.0 + 2.0 * bits(5) + 4.0 * (bits(6) ^ bits(5)) + 8.0 * (bits(7) ^ bits(6) ^ bits(5)))) * inv_sqrt170;
        
        return Complex(re, im);
    }
    
    SoftVec demodulate_symbol(Complex sym, ModulationScheme scheme, double noise_var) {
        switch (scheme) {
            case ModulationScheme::QPSK:
                return demodulate_qpsk(sym, noise_var);
            case ModulationScheme::QAM16:
                return demodulate_qam16(sym, noise_var);
            case ModulationScheme::QAM64:
                return demodulate_qam64(sym, noise_var);
            case ModulationScheme::QAM256:
                return demodulate_qam256(sym, noise_var);
            default:
                return demodulate_qpsk(sym, noise_var);
        }
    }
    
    SoftVec demodulate_qpsk(Complex sym, double noise_var) {
        double scale = 2.0 / noise_var;
        SoftVec llr(2);
        llr(0) = scale * sym.real();
        llr(1) = scale * sym.imag();
        return llr;
    }
    
    SoftVec demodulate_qam16(Complex sym, double noise_var) {
        double inv_sqrt10 = 1.0 / std::sqrt(10.0);
        double scale = 2.0 * inv_sqrt10 / noise_var;
        SoftVec llr(4);
        
        llr(0) = scale * sym.real();
        llr(1) = -scale * (std::abs(sym.real()) - 2.0 * inv_sqrt10);
        llr(2) = scale * sym.imag();
        llr(3) = -scale * (std::abs(sym.imag()) - 2.0 * inv_sqrt10);
        
        return llr;
    }
    
    SoftVec demodulate_qam64(Complex sym, double noise_var) {
        double inv_sqrt42 = 1.0 / std::sqrt(42.0);
        double scale = 2.0 * inv_sqrt42 / noise_var;
        SoftVec llr(6);
        
        double re = sym.real();
        double im = sym.imag();
        
        llr(0) = scale * re;
        llr(1) = -scale * (std::abs(re) - 2.0 * inv_sqrt42);
        llr(2) = -scale * (std::abs(std::abs(re) - 4.0 * inv_sqrt42) - 2.0 * inv_sqrt42);
        llr(3) = scale * im;
        llr(4) = -scale * (std::abs(im) - 2.0 * inv_sqrt42);
        llr(5) = -scale * (std::abs(std::abs(im) - 4.0 * inv_sqrt42) - 2.0 * inv_sqrt42);
        
        return llr;
    }
    
    SoftVec demodulate_qam256(Complex sym, double noise_var) {
        double inv_sqrt170 = 1.0 / std::sqrt(170.0);
        double scale = 2.0 * inv_sqrt170 / noise_var;
        SoftVec llr(8);
        
        double re = sym.real();
        double im = sym.imag();
        
        llr(0) = scale * re;
        llr(1) = -scale * (std::abs(re) - 2.0 * inv_sqrt170);
        llr(2) = -scale * (std::abs(std::abs(re) - 4.0 * inv_sqrt170) - 2.0 * inv_sqrt170);
        llr(3) = -scale * (std::abs(std::abs(std::abs(re) - 4.0 * inv_sqrt170) - 2.0 * inv_sqrt170) - 2.0 * inv_sqrt170);
        llr(4) = scale * im;
        llr(5) = -scale * (std::abs(im) - 2.0 * inv_sqrt170);
        llr(6) = -scale * (std::abs(std::abs(im) - 4.0 * inv_sqrt170) - 2.0 * inv_sqrt170);
        llr(7) = -scale * (std::abs(std::abs(std::abs(im) - 4.0 * inv_sqrt170) - 2.0 * inv_sqrt170) - 2.0 * inv_sqrt170);
        
        return llr;
    }
};

std::unique_ptr<IModulator> create_modulator() {
    return std::make_unique<ModulatorImpl>();
}

} // namespace phy
} // namespace nr
