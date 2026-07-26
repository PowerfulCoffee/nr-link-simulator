#include "phy/PhyInterfaces.h"
#include "phy/ModuleFactory.h"
#include "common/NrTables.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>

using namespace nr;
using namespace nr::phy;
using namespace std;

int main() {
    auto mod = create_modulator();
    mt19937 rng(42);
    normal_distribution<double> nd(0.0, 1.0);

    cout << "=== Simple QPSK LLR test ===" << endl;
    const int N = 10000;
    BitVec bits(2*N);
    for (int i = 0; i < 2*N; i++) bits[i] = rand() % 2;

    ComplexVec sym = mod->modulate(bits, ModulationScheme::QPSK);
    cout << "First 5 QPSK symbols:" << endl;
    for (int i = 0; i < 5; i++) {
        cout << "  bits[" << i << "]=(" << (int)bits[2*i] << "," << (int)bits[2*i+1] << ")"
             << "  sym=(" << sym(i).real() << "," << sym(i).imag() << ")" << endl;
    }

    double snr_db = 10.0;
    double noise_var = 1.0 / pow(10.0, snr_db/10.0);
    cout << "noise_var = " << noise_var << " (Es/N0=" << snr_db << "dB)" << endl;
    double sigma = sqrt(noise_var / 2.0);

    ComplexVec noisy(N);
    for (int i = 0; i < N; i++) {
        noisy(i) = Complex(sym(i).real() + sigma*nd(rng), sym(i).imag() + sigma*nd(rng));
    }

    vector<double> nvars(N, noise_var);
    SoftVec llr = mod->demodulate(noisy, ModulationScheme::QPSK, nvars);
    cout << "First 10 LLRs (should be positive for bit=1, negative for bit=0):" << endl;
    int err = 0;
    for (int i = 0; i < 2*N; i++) {
        int hard = (llr[i] > 0) ? 1 : 0;
        if (hard != bits[i]) err++;
        if (i < 20) {
            cout << "  bit=" << (int)bits[i] << " llr=" << llr[i] << " hard=" << hard << endl;
        }
    }
    cout << "BER = " << (double)err/(2*N) << endl;

    cout << "\n=== Test noiseless case ===" << endl;
    SoftVec llr_clean = mod->demodulate(sym, ModulationScheme::QPSK, vector<double>(N, 1e-10));
    int err_clean = 0;
    for (int i = 0; i < 2*N; i++) {
        int hard = (llr_clean[i] > 0) ? 1 : 0;
        if (hard != bits[i]) err_clean++;
        if (i < 20) {
            cout << "  bit=" << (int)bits[i] << " llr=" << llr_clean[i] << " hard=" << hard << endl;
        }
    }
    cout << "Noiseless BER = " << (double)err_clean/(2*N) << endl;

    return 0;
}
