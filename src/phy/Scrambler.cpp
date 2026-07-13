#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <vector>
#include <cstdint>

namespace nr {
namespace phy {

class ScramblerImpl : public IScrambler {
public:
    BitVec scramble(const BitVec& bits, uint32_t cinit) override {
        int len = static_cast<int>(bits.size());
        std::vector<uint8_t> seq;
        generate_pn_sequence(cinit, len, seq);
        
        BitVec scrambled(len);
        for (int i = 0; i < len; i++) {
            scrambled[i] = bits[i] ^ seq[i];
        }
        return scrambled;
    }
    
    SoftVec descramble(const SoftVec& llr, uint32_t cinit) override {
        int len = static_cast<int>(llr.size());
        std::vector<uint8_t> seq;
        generate_pn_sequence(cinit, len, seq);
        
        SoftVec descrambled(len);
        for (int i = 0; i < len; i++) {
            if (seq[i] == 1) {
                descrambled[i] = -llr[i];
            } else {
                descrambled[i] = llr[i];
            }
        }
        return descrambled;
    }
};

std::unique_ptr<IScrambler> create_scrambler() {
    return std::make_unique<ScramblerImpl>();
}

}
}
