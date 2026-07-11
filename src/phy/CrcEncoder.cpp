#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"

namespace nr {
namespace phy {

class CrcEncoderImpl : public ICrcEncoder {
public:
    BitVec encode(const BitVec& bits, int crc_len) override {
        uint32_t poly;
        if (crc_len == 24) {
            poly = 0x1864CFB;
        } else {
            poly = 0x11021;
            crc_len = 16;
        }
        
        int n = bits.n_elem;
        BitVec result(n + crc_len);
        result(arma::span(0, n - 1)) = bits;
        
        uint32_t crc = 0;
        for (int i = 0; i < n; i++) {
            uint32_t bit = ((crc >> (crc_len - 1)) ^ bits(i)) & 1;
            crc <<= 1;
            if (bit) {
                crc ^= poly;
            }
        }
        
        for (int i = 0; i < crc_len; i++) {
            result(n + i) = (crc >> (crc_len - 1 - i)) & 1;
        }
        
        return result;
    }
    
    std::pair<BitVec, bool> decode(const BitVec& bits, int crc_len) override {
        if (crc_len == 0) {
            crc_len = (bits.n_elem > 3824 + 24) ? 24 : 16;
        }
        
        if (crc_len != 16 && crc_len != 24) {
            crc_len = 24;
        }
        
        int n = bits.n_elem - crc_len;
        if (n <= 0) {
            return {BitVec(), false};
        }
        
        BitVec info_bits = bits(arma::span(0, n - 1));
        BitVec encoded = encode(info_bits, crc_len);
        
        bool crc_ok = true;
        for (int i = 0; i < bits.n_elem; i++) {
            if (bits(i) != encoded(i)) {
                crc_ok = false;
                break;
            }
        }
        
        return {info_bits, crc_ok};
    }
};

std::unique_ptr<ICrcEncoder> create_crc_encoder() {
    return std::make_unique<CrcEncoderImpl>();
}

} // namespace phy
} // namespace nr
