#include "phy/PhyInterfaces.h"
#include "phy/ModuleFactory.h"
#include <iostream>

using namespace nr;
using namespace nr::phy;

int main() {
    std::cout << "Testing CRC Encoder/Decoder...\n";

    auto crc = create_crc_encoder();

    const int num_bits = 100;
    BitVec bits(num_bits);
    for (int i = 0; i < num_bits; i++) {
        bits(i) = i % 2;
    }

    BitVec encoded = crc->encode(bits, 16);
    auto decoded_pair = crc->decode(encoded, 16);
    BitVec decoded = decoded_pair.first;
    bool crc_ok = decoded_pair.second;

    if (!crc_ok) {
        std::cerr << "FAIL: CRC check failed for clean data!\n";
        return 1;
    }

    if (decoded.n_elem != bits.n_elem) {
        std::cerr << "FAIL: Length mismatch! Expected " << bits.n_elem << ", got " << decoded.n_elem << "\n";
        return 1;
    }

    for (int i = 0; i < (int)bits.n_elem; i++) {
        if (decoded(i) != bits(i)) {
            std::cerr << "FAIL: Bit mismatch at position " << i << "!\n";
            return 1;
        }
    }

    std::cout << "PASS: CRC test passed!\n";
    return 0;
}
