#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/ModuleFactory.h"
#include <iostream>
#include <iomanip>
#include <random>

using namespace nr;
using namespace nr::phy;

int main() {
    auto mod = create_modulator();
    auto scram = create_scrambler();

    int mcs = 27;
    int qm = mcs_to_bits_per_symbol(mcs);
    ModulationScheme scheme = mcs_to_modulation(mcs);
    std::cout << "Testing QAM" << (qm==6?"64":(qm==8?"256":"?")) << " demod LLR sign...\n";

    // Generate some bits, modulate, demodulate with tiny noise, check hard decisions match
    int n_bits = 6000;
    BitVec bits(n_bits);
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> bit_dist(0, 1);
    for (int i = 0; i < n_bits; i++) bits[i] = bit_dist(rng);

    uint32_t seed = 123;
    BitVec scrambled = scram->scramble(bits, seed);

    ComplexVec syms = mod->modulate(scrambled, scheme);
    std::cout << "n_symbols = " << syms.n_elem << ", expected = " << n_bits/qm << "\n";

    // Check power normalization
    double pwr = 0;
    for (int i = 0; i < (int)syms.n_elem; i++) {
        pwr += std::norm(syms(i));
    }
    pwr /= syms.n_elem;
    std::cout << "Average symbol power = " << pwr << " (should be ~1.0)\n";

    // Demodulate with very small noise
    double noise_var = 1e-10;
    SoftVec llr = mod->demodulate(syms, scheme, noise_var);

    // Hard decisions: LLR > 0 -> bit 0, LLR < 0 -> bit 1
    int bit_errors = 0;
    for (int i = 0; i < n_bits; i++) {
        int hard_bit = (llr[i] < 0) ? 1 : 0;
        if (hard_bit != scrambled[i]) {
            bit_errors++;
            if (bit_errors <= 5) {
                std::cout << "  Bit error at " << i << ": llr=" << llr[i]
                          << " hard=" << hard_bit << " expected=" << (int)scrambled[i] << "\n";
            }
        }
    }
    std::cout << "Demod hard decision errors (before descramble): " << bit_errors << "/" << n_bits << "\n";

    // Descramble
    SoftVec descrambled = scram->descramble(llr, seed);

    int bit_errors2 = 0;
    for (int i = 0; i < n_bits; i++) {
        int hard_bit = (descrambled[i] < 0) ? 1 : 0;
        if (hard_bit != bits[i]) {
            bit_errors2++;
        }
    }
    std::cout << "After descramble hard errors: " << bit_errors2 << "/" << n_bits << "\n";

    return 0;
}
