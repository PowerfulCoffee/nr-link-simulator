#include "phy/PhyInterfaces.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace nr {
namespace phy {

class RateMatcherImpl : public IRateMatcher {
public:
    BitVec rate_match(const BitVec& coded_bits, int E, int rv, int bgn, int zc) override {
        (void)bgn;
        int N_cb = static_cast<int>(coded_bits.size());
        int k0 = calculate_k0(rv, N_cb, zc);

        BitVec matched(E);
        for (int i = 0; i < E; i++) {
            int idx = (k0 + i) % N_cb;
            if (idx < 0) idx += N_cb;
            matched[i] = coded_bits[idx];
        }

        return matched;
    }

    SoftVec rate_recover(const SoftVec& llr, int N, int rv, int bgn, int zc) override {
        (void)bgn;
        int E = static_cast<int>(llr.size());
        int N_cb = N;

        SoftVec recovered(N_cb, 0.0);
        int k0 = calculate_k0(rv, N_cb, zc);

        for (int i = 0; i < E; i++) {
            int idx = (k0 + i) % N_cb;
            if (idx < 0) idx += N_cb;
            if (idx >= 0 && idx < N_cb) {
                recovered[idx] += llr[i];
            }
        }

        return recovered;
    }

private:
    int calculate_k0(int rv, int N_cb, int zc) const {
        int Z = zc;
        if (Z <= 0) Z = 1;
        switch (rv % 4) {
            case 0: return 0;
            case 1: return N_cb / 4;
            case 2: return N_cb / 2;
            case 3: return 3 * N_cb / 4;
            default: return 0;
        }
    }
};

std::unique_ptr<IRateMatcher> create_rate_matcher() {
    return std::make_unique<RateMatcherImpl>();
}

}
}
