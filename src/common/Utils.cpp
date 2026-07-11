#include "common/Types.h"
#include <random>
#include <cmath>

namespace nr {

void TransportBlock::generate_random_bits(uint64_t seed) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> dist(0, 1);
    bits.set_size(tb_size);
    for (int i = 0; i < tb_size; i++) {
        bits(i) = dist(rng);
    }
}

namespace utils {

double db_to_linear(double db) {
    return std::pow(10.0, db / 10.0);
}

double linear_to_db(double linear) {
    return 10.0 * std::log10(linear);
}

Complex qpsk_modulate(Bit b0, Bit b1) {
    double real = (b0 == 0) ? 1.0 / std::sqrt(2.0) : -1.0 / std::sqrt(2.0);
    double imag = (b1 == 0) ? 1.0 / std::sqrt(2.0) : -1.0 / std::sqrt(2.0);
    return Complex(real, imag);
}

std::pair<SoftBit, SoftBit> qpsk_demodulate(Complex sym, double noise_var) {
    double scale = 1.0 / (noise_var * std::sqrt(2.0));
    return {sym.real() * scale * -2.0, sym.imag() * scale * -2.0};
}

} // namespace utils
} // namespace nr
