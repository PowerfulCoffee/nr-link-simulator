#include "phy/PhyInterfaces.h"
#include "phy/ModuleFactory.h"
#include <iostream>
#include <cmath>

using namespace nr;
using namespace nr::phy;

int main() {
    std::cout << "Testing QPSK Modulator...\n";

    auto mod = create_modulator();

    const int num_bits = 200;
    BitVec bits(num_bits);
    for (int i = 0; i < num_bits; i++) {
        bits(i) = i % 2;
    }

    ComplexVec syms = mod->modulate(bits, ModulationScheme::QPSK);
    if (syms.n_elem != 100) {
        std::cerr << "FAIL: QPSK symbol count wrong: expected 100, got " << syms.n_elem << "\n";
        return 1;
    }

    for (int i = 0; i < (int)syms.n_elem; i++) {
        double pwr = std::norm(syms(i));
        if (std::abs(pwr - 1.0) > 0.1) {
            std::cerr << "FAIL: QPSK symbol " << i << " power not 1: " << pwr << "\n";
            return 1;
        }
    }

    std::cout << "PASS: QPSK modulation symbol power is 1!\n";
    std::cout << "PASS: Modulator test passed!\n";
    return 0;
}
