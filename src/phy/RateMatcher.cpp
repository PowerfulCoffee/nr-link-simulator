#include "phy/PhyInterfaces.h"
#include <cmath>
#include <algorithm>

namespace nr {
namespace phy {

class RateMatcherImpl : public IRateMatcher {
public:
    BitVec rate_match(const BitVec& coded_bits, int E, int rv, int bgn, int zc) override {
        (void)bgn;
        (void)zc;
        int N = static_cast<int>(coded_bits.n_elem);
        int N_cb = N;
        int k0 = calculate_k0(rv, N_cb);

        BitVec matched(E);
        for (int i = 0; i < E; i++) {
            int idx = (k0 + i) % N_cb;
            if (idx < 0) idx += N_cb;
            matched(i) = coded_bits(idx);
        }

        return matched;
    }

    SoftVec rate_recover(const SoftVec& llr, int N, int rv, int bgn, int zc) override {
        (void)bgn;
        (void)zc;
        int E = static_cast<int>(llr.n_elem);
        int N_cb = N;

        SoftVec recovered(N, arma::fill::zeros);
        int k0 = calculate_k0(rv, N_cb);

        for (int i = 0; i < E; i++) {
            int idx = (k0 + i) % N_cb;
            if (idx < 0) idx += N_cb;
            if (idx >= 0 && idx < N) {
                recovered(idx) += llr(i);
            }
        }

        return recovered;
    }

private:
    int calculate_k0(int rv, int N_cb) const {
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
