#include "phy/PhyInterfaces.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdexcept>

namespace nr {
namespace phy {

namespace {

int calculate_k0(int rv, int bgn, int zc, int n_cb) {
    int n_cols = (bgn == 1) ? 66 : 50;
    static const int coeffs_bg1[] = {0, 17, 33, 56};
    static const int coeffs_bg2[] = {0, 13, 25, 43};
    const int* coeffs = (bgn == 1) ? coeffs_bg1 : coeffs_bg2;
    int c = coeffs[rv % 4];
    int k0 = (c * n_cb) / (n_cols * zc);
    k0 = k0 * zc;
    return k0;
}

void apply_bicm_interleave(BitVec& bits, int qm) {
    if (qm <= 1) return;
    int E = static_cast<int>(bits.size());
    if (E % qm != 0) return;
    int n_rows = E / qm;
    BitVec out(E);
    for (int j = 0; j < n_rows; j++) {
        for (int i = 0; i < qm; i++) {
            out[i + j * qm] = bits[i * n_rows + j];
        }
    }
    bits = out;
}

void apply_bicm_deinterleave(SoftVec& llr, int qm) {
    if (qm <= 1) return;
    int E = static_cast<int>(llr.size());
    if (E % qm != 0) return;
    int n_rows = E / qm;
    SoftVec out(E, 0.0);
    for (int j = 0; j < n_rows; j++) {
        for (int i = 0; i < qm; i++) {
            out[i * n_rows + j] = llr[i + j * qm];
        }
    }
    llr = out;
}

}

class RateMatcherImpl : public IRateMatcher {
public:
    BitVec rate_match(const BitVec& coded_bits, int E, int rv, int bgn, int zc,
                      int qm, int /*n_filler*/) override {
        int n_cb = static_cast<int>(coded_bits.size());

        int k0 = calculate_k0(rv, bgn, zc, n_cb);

        BitVec matched(E);
        for (int i = 0; i < E; i++) {
            int idx = (k0 + i) % n_cb;
            if (idx < 0) idx += n_cb;
            matched[i] = coded_bits[idx];
        }

        if (qm > 1) {
            apply_bicm_interleave(matched, qm);
        }

        return matched;
    }

    SoftVec rate_recover(const SoftVec& llr, int n_cb, int rv, int bgn, int zc,
                          int qm, int /*n_filler*/) override {
        SoftVec llr_deint = llr;
        if (qm > 1) {
            apply_bicm_deinterleave(llr_deint, qm);
        }

        int E = static_cast<int>(llr_deint.size());

        int k0 = calculate_k0(rv, bgn, zc, n_cb);

        SoftVec recovered(n_cb, 0.0);
        for (int i = 0; i < E; i++) {
            int idx = (k0 + i) % n_cb;
            if (idx < 0) idx += n_cb;
            if (idx >= 0 && idx < n_cb) {
                recovered[idx] += llr_deint[i];
            }
        }

        return recovered;
    }
};

std::unique_ptr<IRateMatcher> create_rate_matcher() {
    return std::make_unique<RateMatcherImpl>();
}

}
}
