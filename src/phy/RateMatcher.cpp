#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <cmath>
#include <algorithm>

namespace nr {
namespace phy {

class RateMatcherImpl : public IRateMatcher {
public:
    BitVec rate_match(const BitVec& coded_bits, int E, int rv, int bgn, int zc) override {
        int N = coded_bits.n_elem;
        BitVec circular_buffer(N);
        
        for (int i = 0; i < N; i++) {
            circular_buffer(i) = coded_bits(i);
        }
        
        int k0 = calculate_k0(rv, bgn, zc, N);
        
        BitVec matched(E);
        for (int i = 0; i < E; i++) {
            int idx = (k0 + i) % N;
            matched(i) = circular_buffer(idx);
        }
        
        return matched;
    }
    
    SoftVec rate_recover(const SoftVec& llr, int N, int rv, int bgn, int zc) override {
        int E = llr.n_elem;
        SoftVec recovered(N, arma::fill::zeros);
        
        int k0 = calculate_k0(rv, bgn, zc, N);
        
        for (int i = 0; i < E; i++) {
            int idx = (k0 + i) % N;
            recovered(idx) += llr(i);
        }
        
        return recovered;
    }

private:
    int calculate_k0(int rv, int bgn, int zc, int N) const {
        int Z = zc;
        int N_cb;
        
        if (bgn == 1) {
            N_cb = 66 * Z;
        } else {
            N_cb = 50 * Z;
        }
        
        if (N_cb > N) {
            N_cb = N;
        }
        
        int rv_offset;
        switch (rv % 4) {
            case 0: rv_offset = 0; break;
            case 1: rv_offset = N_cb / 4; break;
            case 2: rv_offset = N_cb / 2; break;
            case 3: rv_offset = 3 * N_cb / 4; break;
            default: rv_offset = 0;
        }
        
        return rv_offset;
    }
};

std::unique_ptr<IRateMatcher> create_rate_matcher() {
    return std::make_unique<RateMatcherImpl>();
}

} // namespace phy
} // namespace nr
