#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <cmath>
#include <complex>
#include <armadillo>

namespace nr {
namespace phy {

namespace {

inline int gray_encode(int k) { return k ^ (k >> 1); }
inline int gray_decode(int g) {
    int b = g;
    b ^= (b >> 1); b ^= (b >> 2); b ^= (b >> 4); b ^= (b >> 8);
    return b;
}

void pam_mod(const uint8_t* bits, int m, double a, double& out) {
    int M = 1 << m;
    int g = 0;
    for (int i = 0; i < m; i++) g = (g << 1) | bits[i];
    int k = gray_decode(g);
    out = (2.0*k - M + 1.0) * a;
}

void pam_demod(double y, double var_per_dim, double a, int m, double* llr_out) {
    int M = 1 << m;
    for (int b = 0; b < m; b++) {
        double d0 = 1e30, d1 = 1e30;
        for (int k = 0; k < M; k++) {
            double s = (2.0*k - M + 1.0) * a;
            double dist = (y - s) * (y - s);
            int g = gray_encode(k);
            int bit_val = (g >> (m - 1 - b)) & 1;
            if (bit_val == 0) { if (dist < d0) d0 = dist; }
            else              { if (dist < d1) d1 = dist; }
        }
        llr_out[b] = (d1 - d0) / (2.0 * var_per_dim);
    }
}

double pam_norm(int m) {
    int M = 1 << m;
    double sum = 0.0;
    for (int k = 0; k < M; k++) {
        double v = (2.0*k - M + 1.0);
        sum += v*v;
    }
    double ms = sum / M;
    return 1.0 / std::sqrt(2.0 * ms);
}

}

class ModulatorImpl : public IModulator {
public:
    ComplexVec modulate(const BitVec& bits, ModulationScheme scheme) override {
        int bps = mod_to_bits_per_symbol(scheme);
        int n = static_cast<int>(bits.size()) / bps;
        ComplexVec sym(n);

        for (int i = 0; i < n; i++) {
            if (scheme == ModulationScheme::BPSK) {
                sym(i) = Complex(bits[i] ? -1.0 : 1.0, 0.0);
            } else if (scheme == ModulationScheme::QPSK) {
                double a = 1.0 / std::sqrt(2.0);
                double re = bits[2*i] ? -a : a;
                double im = bits[2*i+1] ? -a : a;
                sym(i) = Complex(re, im);
            } else if (scheme == ModulationScheme::QAM16) {
                double a = pam_norm(2);
                double re, im;
                uint8_t rb[2] = {bits[4*i], bits[4*i+1]};
                uint8_t ib[2] = {bits[4*i+2], bits[4*i+3]};
                pam_mod(rb, 2, a, re); pam_mod(ib, 2, a, im);
                sym(i) = Complex(re, im);
            } else if (scheme == ModulationScheme::QAM64) {
                double a = pam_norm(3);
                double re, im;
                uint8_t rb[3] = {bits[6*i], bits[6*i+1], bits[6*i+2]};
                uint8_t ib[3] = {bits[6*i+3], bits[6*i+4], bits[6*i+5]};
                pam_mod(rb, 3, a, re); pam_mod(ib, 3, a, im);
                sym(i) = Complex(re, im);
            } else if (scheme == ModulationScheme::QAM256) {
                double a = pam_norm(4);
                double re, im;
                uint8_t rb[4] = {bits[8*i], bits[8*i+1], bits[8*i+2], bits[8*i+3]};
                uint8_t ib[4] = {bits[8*i+4], bits[8*i+5], bits[8*i+6], bits[8*i+7]};
                pam_mod(rb, 4, a, re); pam_mod(ib, 4, a, im);
                sym(i) = Complex(re, im);
            } else {
                double a = 1.0 / std::sqrt(2.0);
                sym(i) = Complex(bits[2*i] ? -a : a, bits[2*i+1] ? -a : a);
            }
        }
        return sym;
    }

    SoftVec demodulate(const ComplexVec& sym, ModulationScheme scheme, double noise_var) override {
        int bps = mod_to_bits_per_symbol(scheme);
        int n = static_cast<int>(sym.n_elem);
        SoftVec llr(n * bps);
        double var_per_dim = std::max(noise_var / 2.0, 1e-12);

        for (int i = 0; i < n; i++) {
            if (scheme == ModulationScheme::BPSK) {
                llr[i] = 2.0 * sym(i).real() / var_per_dim;
            } else if (scheme == ModulationScheme::QPSK) {
                double a = 1.0 / std::sqrt(2.0);
                llr[2*i]   = 2.0 * a * sym(i).real() / var_per_dim;
                llr[2*i+1] = 2.0 * a * sym(i).imag() / var_per_dim;
            } else if (scheme == ModulationScheme::QAM16) {
                double a = pam_norm(2);
                double lr[2], li[2];
                pam_demod(sym(i).real(), var_per_dim, a, 2, lr);
                pam_demod(sym(i).imag(), var_per_dim, a, 2, li);
                for (int b = 0; b < 2; b++) { llr[4*i+b]=lr[b]; llr[4*i+2+b]=li[b]; }
            } else if (scheme == ModulationScheme::QAM64) {
                double a = pam_norm(3);
                double lr[3], li[3];
                pam_demod(sym(i).real(), var_per_dim, a, 3, lr);
                pam_demod(sym(i).imag(), var_per_dim, a, 3, li);
                for (int b = 0; b < 3; b++) { llr[6*i+b]=lr[b]; llr[6*i+3+b]=li[b]; }
            } else if (scheme == ModulationScheme::QAM256) {
                double a = pam_norm(4);
                double lr[4], li[4];
                pam_demod(sym(i).real(), var_per_dim, a, 4, lr);
                pam_demod(sym(i).imag(), var_per_dim, a, 4, li);
                for (int b = 0; b < 4; b++) { llr[8*i+b]=lr[b]; llr[8*i+4+b]=li[b]; }
            } else {
                double a = 1.0 / std::sqrt(2.0);
                llr[2*i]   = 2.0 * a * sym(i).real() / var_per_dim;
                llr[2*i+1] = 2.0 * a * sym(i).imag() / var_per_dim;
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
