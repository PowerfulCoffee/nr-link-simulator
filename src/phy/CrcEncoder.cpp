#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <vector>
#include <cstdint>
#include <utility>

namespace nr {
namespace phy {

class CrcEncoderImpl : public ICrcEncoder {
public:
    BitVec encode(const BitVec& bits, int crc_len) override {
        uint32_t poly;
        if (crc_len == 24) {
            poly = 0x1864CFB;
        } else if (crc_len == 16) {
            poly = 0x11021;
        } else {
            crc_len = 24;
            poly = 0x1864CFB;
        }
        
        int n = static_cast<int>(bits.size());
        BitVec result(n + crc_len, 0);
        for (int i = 0; i < n; i++) {
            result[i] = bits[i];
        }
        
        uint32_t crc = 0;
        for (int i = 0; i < n; i++) {
            uint32_t bit = ((crc >> (crc_len - 1)) ^ bits[i]) & 1;
            crc <<= 1;
            crc &= (1U << crc_len) - 1;
            if (bit) {
                crc ^= poly;
            }
        }
        
        for (int i = 0; i < crc_len; i++) {
            result[n + i] = (crc >> (crc_len - 1 - i)) & 1;
            crc <<= 1;
            crc &= (1U << crc_len) - 1;
        }
        
        return result;
    }
    
    std::pair<BitVec, bool> decode(const BitVec& bits, int crc_len) override {
        if (crc_len <= 0) {
            crc_len = (static_cast<int>(bits.size()) > 3824 + 24) ? 24 : 16;
        }
        if (crc_len != 16 && crc_len != 24) {
            crc_len = 24;
        }
        
        int n = static_cast<int>(bits.size()) - crc_len;
        if (n <= 0) {
            return {BitVec(), false};
        }
        
        BitVec info_bits(bits.begin(), bits.begin() + n);
        BitVec encoded = encode(info_bits, crc_len);
        
        bool crc_ok = true;
        for (int i = 0; i < static_cast<int>(bits.size()); i++) {
            if (bits[i] != encoded[i]) {
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

}
}
