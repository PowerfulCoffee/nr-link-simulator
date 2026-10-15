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

// Rebuild the circular buffer exactly as TS 38.212 5.4.2 defines it:
// the buffer holds the systematic part (filler bits kept in place as
// <NULL>) followed by the parity part; the first 2*Zc systematic bits are
// already punctured by the encoder and are not part of the buffer.
// The compacted input has the NULL (filler) bits removed; they are
// re-inserted here so that k0 is computed on the full Ncb and bit
// selection skips the filler positions as required by the spec.
struct NullSkipBuffer {
    std::vector<int> pos;   // buffer position -> compacted index, -1 = NULL
    int n_cb;

    NullSkipBuffer(int n_compact, int bgn, int zc, int n_filler) {
        n_cb = n_compact + n_filler;
        pos.assign(n_cb, -1);
        if (n_filler <= 0) {
            for (int i = 0; i < n_cb; i++) pos[i] = i;
            return;
        }
        int n_sys_p = ((bgn == 1) ? 22 : 10) * zc - 2 * zc; // systematic minus punctured
        int filler_start = n_sys_p - n_filler;
        if (filler_start < 0) filler_start = 0;
        int c = 0;
        for (int i = 0; i < n_cb; i++) {
            if (i >= filler_start && i < filler_start + n_filler) continue; // NULL
            pos[i] = c++;
        }
    }
};

// TS 38.212 5.4.2.2 bit interleaving: f_{i + j*Qm} = e_{i*E/Qm + j}
void apply_bicm_interleave(BitVec& bits, int qm) {
    if (qm <= 1) return;
    int E = static_cast<int>(bits.size());
    if (E % qm != 0) {
        throw std::runtime_error("rate matching: E must be a multiple of Qm per TS 38.212 5.4.2.2");
    }
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
                      int qm, int n_filler) override {
        int n_compact = static_cast<int>(coded_bits.size());
        NullSkipBuffer buf(n_compact, bgn, zc, n_filler);

        int k0 = calculate_k0(rv, bgn, zc, buf.n_cb);

        // Bit selection per TS 38.212 5.4.2.1: start at k0, walk the circular
        // buffer, skip <NULL> (filler) positions, collect E bits. E may exceed
        // one buffer lap (repetition is legal); only a full lap with zero
        // collected bits (all-NULL buffer) is an error.
        BitVec e;
        e.reserve(E);
        int p = k0;
        int n_null_skipped = 0; // NULL positions seen in the current lap
        while (static_cast<int>(e.size()) < E) {
            if (buf.pos[p] >= 0) {
                e.push_back(coded_bits[buf.pos[p]]);
            } else {
                n_null_skipped++;
            }
            p = (p + 1) % buf.n_cb;
            if (p == k0) { // completed one full lap from k0
                if (n_null_skipped >= buf.n_cb) {
                    throw std::runtime_error("rate matching: buffer exhausted (E too large)");
                }
                n_null_skipped = 0;
            }
        }

        if (qm > 1) {
            apply_bicm_interleave(e, qm);
        }
        return e;
    }

    SoftVec rate_recover(const SoftVec& llr, int n_cb_compact, int rv, int bgn, int zc,
                          int qm, int n_filler) override {
        SoftVec llr_e = llr;
        if (qm > 1) {
            apply_bicm_deinterleave(llr_e, qm);
        }

        int E = static_cast<int>(llr_e.size());
        NullSkipBuffer buf(n_cb_compact, bgn, zc, n_filler);

        int k0 = calculate_k0(rv, bgn, zc, buf.n_cb);

        // Walk the same selection order, accumulating LLRs into the compacted
        // buffer (NULL positions contribute nothing). Multiple laps are legal
        // and repeat selections add up (soft combining), mirroring rate_match.
        SoftVec recovered(n_cb_compact, 0.0);
        int p = k0;
        int k = 0;
        int n_null_skipped = 0;
        while (k < E) {
            if (buf.pos[p] >= 0) {
                recovered[buf.pos[p]] += llr_e[k];
                k++;
            } else {
                n_null_skipped++;
            }
            p = (p + 1) % buf.n_cb;
            if (p == k0) { // completed one full lap from k0
                if (n_null_skipped >= buf.n_cb) {
                    break; // degenerate all-NULL buffer: nothing recoverable
                }
                n_null_skipped = 0;
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
