#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/ModuleFactory.h"
#include <iostream>
#include <iomanip>
#include <memory>

using namespace nr;
using namespace nr::phy;

int main() {
    auto mod = create_modulator();
    auto layer_map = create_layer_mapper();
    auto precode = create_precoder();
    auto scrambler = create_scrambler();

    const int n_bits = 7800;
    BitVec bits(n_bits, arma::fill::zeros);
    for (int i = 0; i < n_bits; i++) {
        bits(i) = i % 2;
    }

    uint32_t cinit = 1;
    BitVec scrambled = scrambler->scramble(bits, cinit);

    ComplexVec modulated = mod->modulate(scrambled, ModulationScheme::QPSK);
    std::cout << "Modulated symbols: " << modulated.n_elem << "\n";

    int n_layers = 1;
    ComplexMat layered = layer_map->map(modulated, n_layers);
    std::cout << "Layered: " << layered.n_rows << " x " << layered.n_cols << "\n";

    ComplexMat w(1, 1, arma::fill::zeros);
    w(0,0) = Complex(1.0, 0.0);
    ComplexMat precoded = precode->precode(layered, w);
    std::cout << "Precoded: " << precoded.n_rows << " x " << precoded.n_cols << "\n";

    double noise_var = 0.01;
    ComplexMat rx_symbols = precoded;

    ComplexMat equalized = rx_symbols;

    ComplexVec demapped = layer_map->demap(equalized, n_layers);
    std::cout << "Demapped: " << demapped.n_elem << "\n";

    SoftVec llr = mod->demodulate(demapped, ModulationScheme::QPSK, noise_var);
    std::cout << "LLR: " << llr.n_elem << "\n";

    SoftVec descrambled = scrambler->descramble(llr, cinit);

    BitVec hard_bits(n_bits);
    int errors = 0;
    for (int i = 0; i < n_bits; i++) {
        hard_bits(i) = (descrambled(i) < 0) ? 1 : 0;
        if (hard_bits(i) != bits(i)) {
            errors++;
            if (errors <= 10) {
                std::cout << "Bit error at " << i << ": tx=" << (int)bits(i)
                          << ", rx=" << (int)hard_bits(i) << ", llr=" << descrambled(i) << "\n";
            }
        }
    }
    std::cout << "Total bit errors: " << errors << " / " << n_bits << "\n";

    return (errors == 0) ? 0 : 1;
}
