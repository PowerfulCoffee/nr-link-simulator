#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <cmath>
#include <complex>
#include <armadillo>
#include <vector>

namespace nr {
namespace phy {

namespace {

inline int gray_encode(int k) { return k ^ (k >> 1); }
inline int gray_decode(int g) {
    int b = g;
    b ^= (b >> 1); b ^= (b >> 2); b ^= (b >> 4); b ^= (b >> 8);
    return b;
}

double pam_gray_3gpp(const uint8_t* bits, int m) {
    if (m == 1) {
        return 1.0 - 2.0 * bits[0];
    }
    double rest = pam_gray_3gpp(bits + 1, m - 1);
    return (1.0 - 2.0 * bits[0]) * ((1 << (m - 1)) - rest);
}

void pam_mod_3gpp(const uint8_t* bits, int m, double a, double& out) {
    double val = pam_gray_3gpp(bits, m);
    out = val * a;
}

double qam_norm(int m_per_dim) {
    int M = 1 << m_per_dim;
    std::vector<uint8_t> dummy(m_per_dim, 0);
    double sum = 0.0;
    for (int k = 0; k < M; k++) {
        for (int i = 0; i < m_per_dim; i++) {
            dummy[i] = (k >> (m_per_dim - 1 - i)) & 1;
        }
        double v = pam_gray_3gpp(dummy.data(), m_per_dim);
        sum += v * v;
    }
    double ms = sum / M;
    return 1.0 / std::sqrt(2.0 * ms);
}

inline double logsumexp(const double* vals, int n) {
    if (n <= 0) return -1e30;
    double max_val = vals[0];
    for (int i = 1; i < n; i++) {
        if (vals[i] > max_val) max_val = vals[i];
    }
    if (max_val < -500.0) return max_val;
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double d = vals[i] - max_val;
        if (d > -50.0) {
            sum += std::exp(d);
        }
    }
    return max_val + std::log(sum);
}

void pam_demod_app(double y, double var_per_dim, double a, int m, double* llr_out) {
    int M = 1 << m;
    int half_M = M / 2;
    std::vector<double> s_vals(M);
    uint8_t bits[8];
    for (int k = 0; k < M; k++) {
        for (int i = 0; i < m; i++) {
            bits[i] = (k >> (m - 1 - i)) & 1;
        }
        s_vals[k] = pam_gray_3gpp(bits, m) * a;
    }

    double two_var = 2.0 * var_per_dim;
    std::vector<double> exps0(half_M);
    std::vector<double> exps1(half_M);

    for (int b = 0; b < m; b++) {
        int cnt0 = 0, cnt1 = 0;
        for (int k = 0; k < M; k++) {
            for (int i = 0; i < m; i++) {
                bits[i] = (k >> (m - 1 - i)) & 1;
            }
            double dist = (y - s_vals[k]) * (y - s_vals[k]);
            double exponent = -dist / two_var;
            if (bits[b] == 0) {
                exps0[cnt0++] = exponent;
            } else {
                exps1[cnt1++] = exponent;
            }
        }
        double lse0 = logsumexp(exps0.data(), cnt0);
        double lse1 = logsumexp(exps1.data(), cnt1);
        llr_out[b] = lse1 - lse0;
    }
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
                uint8_t b = bits[i];
                double v = 1.0 - 2.0 * b;
                sym(i) = Complex(v, 0.0);
            } else if (scheme == ModulationScheme::QPSK) {
                double a = 1.0 / std::sqrt(2.0);
                uint8_t bits_i[1] = {bits[2*i]};
                uint8_t bits_q[1] = {bits[2*i+1]};
                double re, im;
                pam_mod_3gpp(bits_i, 1, a, re);
                pam_mod_3gpp(bits_q, 1, a, im);
                sym(i) = Complex(re, im);
            } else if (scheme == ModulationScheme::QAM16) {
                int m = 2;
                double a = qam_norm(m);
                double re, im;
                uint8_t rb[2] = {bits[4*i], bits[4*i+2]};
                uint8_t ib[2] = {bits[4*i+1], bits[4*i+3]};
                pam_mod_3gpp(rb, m, a, re);
                pam_mod_3gpp(ib, m, a, im);
                sym(i) = Complex(re, im);
            } else if (scheme == ModulationScheme::QAM64) {
                int m = 3;
                double a = qam_norm(m);
                double re, im;
                uint8_t rb[3] = {bits[6*i], bits[6*i+2], bits[6*i+4]};
                uint8_t ib[3] = {bits[6*i+1], bits[6*i+3], bits[6*i+5]};
                pam_mod_3gpp(rb, m, a, re);
                pam_mod_3gpp(ib, m, a, im);
                sym(i) = Complex(re, im);
            } else if (scheme == ModulationScheme::QAM256) {
                int m = 4;
                double a = qam_norm(m);
                double re, im;
                uint8_t rb[4] = {bits[8*i], bits[8*i+2], bits[8*i+4], bits[8*i+6]};
                uint8_t ib[4] = {bits[8*i+1], bits[8*i+3], bits[8*i+5], bits[8*i+7]};
                pam_mod_3gpp(rb, m, a, re);
                pam_mod_3gpp(ib, m, a, im);
                sym(i) = Complex(re, im);
            } else {
                double a = 1.0 / std::sqrt(2.0);
                uint8_t bits_i[1] = {bits[2*i]};
                uint8_t bits_q[1] = {bits[2*i+1]};
                double re, im;
                pam_mod_3gpp(bits_i, 1, a, re);
                pam_mod_3gpp(bits_q, 1, a, im);
                sym(i) = Complex(re, im);
            }
        }
        return sym;
    }

    SoftVec demodulate(const ComplexVec& sym, ModulationScheme scheme, double noise_var) override {
        std::vector<double> noise_vars(sym.n_elem, noise_var);
        return demodulate(sym, scheme, noise_vars);
    }

    SoftVec demodulate(const ComplexVec& sym, ModulationScheme scheme,
                       const std::vector<double>& noise_vars) override {
        int bps = mod_to_bits_per_symbol(scheme);
        int n = static_cast<int>(sym.n_elem);
        SoftVec llr(n * bps);

        for (int i = 0; i < n; i++) {
            double nv = (i < (int)noise_vars.size()) ? noise_vars[i] : noise_vars.back();
            double var_per_dim = std::max(nv / 2.0, 1e-12);
            if (scheme == ModulationScheme::BPSK) {
                llr[i] = -2.0 * sym(i).real() / var_per_dim;
            } else if (scheme == ModulationScheme::QPSK) {
                double a = 1.0 / std::sqrt(2.0);
                double lr[1], li[1];
                pam_demod_app(sym(i).real(), var_per_dim, a, 1, lr);
                pam_demod_app(sym(i).imag(), var_per_dim, a, 1, li);
                llr[2*i]   = lr[0];
                llr[2*i+1] = li[0];
            } else if (scheme == ModulationScheme::QAM16) {
                int m = 2;
                double a = qam_norm(m);
                double lr[2], li[2];
                pam_demod_app(sym(i).real(), var_per_dim, a, m, lr);
                pam_demod_app(sym(i).imag(), var_per_dim, a, m, li);
                llr[4*i]   = lr[0];
                llr[4*i+2] = lr[1];
                llr[4*i+1] = li[0];
                llr[4*i+3] = li[1];
            } else if (scheme == ModulationScheme::QAM64) {
                int m = 3;
                double a = qam_norm(m);
                double lr[3], li[3];
                pam_demod_app(sym(i).real(), var_per_dim, a, m, lr);
                pam_demod_app(sym(i).imag(), var_per_dim, a, m, li);
                llr[6*i]   = lr[0];
                llr[6*i+2] = lr[1];
                llr[6*i+4] = lr[2];
                llr[6*i+1] = li[0];
                llr[6*i+3] = li[1];
                llr[6*i+5] = li[2];
            } else if (scheme == ModulationScheme::QAM256) {
                int m = 4;
                double a = qam_norm(m);
                double lr[4], li[4];
                pam_demod_app(sym(i).real(), var_per_dim, a, m, lr);
                pam_demod_app(sym(i).imag(), var_per_dim, a, m, li);
                llr[8*i]   = lr[0];
                llr[8*i+2] = lr[1];
                llr[8*i+4] = lr[2];
                llr[8*i+6] = lr[3];
                llr[8*i+1] = li[0];
                llr[8*i+3] = li[1];
                llr[8*i+5] = li[2];
                llr[8*i+7] = li[3];
            } else {
                double a = 1.0 / std::sqrt(2.0);
                double lr[1], li[1];
                pam_demod_app(sym(i).real(), var_per_dim, a, 1, lr);
                pam_demod_app(sym(i).imag(), var_per_dim, a, 1, li);
                llr[2*i]   = lr[0];
                llr[2*i+1] = li[0];
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
