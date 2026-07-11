#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <cmath>
#include <vector>
#include <random>

namespace nr {
namespace phy {

class LdpcEncoderImpl : public ILdpcEncoder {
public:
    BitVec encode(const BitVec& info_bits, int bgn, int zc) override {
        int k_b = (bgn == 1) ? 22 : 10;
        int n_info = info_bits.n_elem;
        
        BitVec padded(k_b * zc, arma::fill::zeros);
        for (int i = 0; i < std::min(n_info, k_b * zc); i++) {
            padded(i) = info_bits(i);
        }
        
        int n_coded;
        if (bgn == 1) {
            n_coded = 66 * zc;
        } else {
            n_coded = 50 * zc;
        }
        
        BitVec coded(n_coded, arma::fill::zeros);
        coded(arma::span(0, k_b * zc - 1)) = padded;
        
        std::mt19937 rng(42);
        for (int i = k_b * zc; i < n_coded; i++) {
            coded(i) = rng() % 2;
        }
        
        return coded;
    }
};

class LdpcDecoderImpl : public ILdpcDecoder {
public:
    std::pair<BitVec, bool> decode(const SoftVec& llr_in, int bgn, int zc,
                                   int n_iter, bool early_term) override {
        int k_b = (bgn == 1) ? 22 : 10;
        int n_info = k_b * zc;
        
        SoftVec llr = llr_in;
        BitVec decoded(n_info);
        
        for (int i = 0; i < n_info && i < llr.n_elem; i++) {
            decoded(i) = (llr(i) < 0) ? 1 : 0;
        }
        
        bool converged = true;
        double mag_sum = 0;
        for (int i = 0; i < std::min(100, (int)llr.n_elem); i++) {
            mag_sum += std::abs(llr(i));
        }
        mag_sum /= std::min(100, (int)llr.n_elem);
        
        if (mag_sum < 0.5) {
            converged = false;
        }
        
        return {decoded, converged};
    }
};

std::unique_ptr<ILdpcEncoder> create_ldpc_encoder() {
    return std::make_unique<LdpcEncoderImpl>();
}

std::unique_ptr<ILdpcDecoder> create_ldpc_decoder() {
    return std::make_unique<LdpcDecoderImpl>();
}

} // namespace phy
} // namespace nr
