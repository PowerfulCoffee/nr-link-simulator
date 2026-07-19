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
        uint32_t mask;
        if (crc_len == 24) {
            // CRC24A: x^24 + x^23 + x^18 + x^17 + x^14 + x^11 + x^10 + x^7 + x^6 + x^5 + x^4 + x^3 + x + 1
            poly = 0x864CFB;
            mask = 0xFFFFFF;
        } else if (crc_len == 16) {
            // CRC16: x^16 + x^12 + x^5 + 1
            poly = 0x1021;
            mask = 0xFFFF;
        } else {
            crc_len = 24;
            poly = 0x864CFB;
            mask = 0xFFFFFF;
        }
        
        int n = static_cast<int>(bits.size());
        BitVec result(n + crc_len, 0);
        for (int i = 0; i < n; i++) {
            result[i] = bits[i];
        }
        
        uint32_t crc = 0;
        // Process all information bits
        for (int i = 0; i < n; i++) {
            uint32_t msb = (crc >> (crc_len - 1)) & 1;
            crc = ((crc << 1) | bits[i]) & mask;
            if (msb) {
                crc ^= poly;
            }
        }
        // Process crc_len zeros to flush the register
        for (int i = 0; i < crc_len; i++) {
            uint32_t msb = (crc >> (crc_len - 1)) & 1;
            crc = (crc << 1) & mask;
            if (msb) {
                crc ^= poly;
            }
        }
        
        // Extract CRC bits MSB first
        for (int i = 0; i < crc_len; i++) {
            result[n + i] = (crc >> (crc_len - 1 - i)) & 1;
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
