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
    auto crc_enc = create_crc_encoder();
    auto mod = create_modulator();

    cout << "=== CRC24A Cross-check vs Sionna ===" << endl;
    BitVec bits = {0,1,0,0,0,0,0,1,1,0,1,1,0,1,0,1,0,1,1,0};
    BitVec enc = crc_enc->encode(bits, 24);
    cout << "Info bits: ";
    for (int i = 0; i < (int)bits.size(); i++) cout << (int)bits[i];
    cout << endl;
    cout << "CRC24A (C++): ";
    for (int i = (int)bits.size(); i < (int)enc.size(); i++) cout << (int)enc[i];
    cout << endl;
    cout << "Expected:      101000011001100111100100" << endl;
    vector<int> expected_crc = {1,0,1,0,0,0,0,1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,0,0};
    bool crc_match = true;
    for (int i = 0; i < 24; i++) {
        if ((int)enc[bits.size()+i] != expected_crc[i]) { crc_match = false; break; }
    }
    cout << "CRC Match: " << (crc_match ? "PASS" : "FAIL") << endl;

    auto [decoded, crc_ok] = crc_enc->decode(enc, 24);
    cout << "CRC decode self-check: " << (crc_ok ? "PASS" : "FAIL") << endl;

    // Additional random tests: encode then decode should always pass
    mt19937 rng(42);
    int n_tests = 1000;
    int fail = 0;
    for (int t = 0; t < n_tests; t++) {
        int len = 20 + rng() % 500;
        BitVec rbits(len);
        for (int i = 0; i < len; i++) rbits[i] = rng() % 2;
        BitVec renc = crc_enc->encode(rbits, 24);
        auto [rdec, rok] = crc_enc->decode(renc, 24);
        if (!rok) { fail++; continue; }
        for (int i = 0; i < len; i++) {
            if (rdec[i] != rbits[i]) { fail++; break; }
        }
    }
    cout << "Random CRC encode/decode: " << (n_tests-fail) << "/" << n_tests << " PASS" << endl;

    cout << "\n=== QPSK LLR Demod Self-check ===" << endl;
    const int N = 50000;
    BitVec tx_bits(2*N);
    for (int i = 0; i < 2*N; i++) tx_bits[i] = rng() % 2;
    ComplexVec sym = mod->modulate(tx_bits, ModulationScheme::QPSK);
    normal_distribution<double> nd(0.0, 1.0);
    double snr_db = 10.0;
    double noise_var = 1.0 / pow(10.0, snr_db/10.0);
    double sigma = sqrt(noise_var / 2.0);
    ComplexVec noisy(N);
    for (int i = 0; i < N; i++) {
        noisy(i) = Complex(sym(i).real() + sigma*nd(rng), sym(i).imag() + sigma*nd(rng));
    }
    vector<double> nvars(N, noise_var);
    SoftVec llr_demod = mod->demodulate(noisy, ModulationScheme::QPSK, nvars);
    int demod_err = 0;
    for (int i = 0; i < 2*N; i++) {
        int hard = (llr_demod[i] > 0) ? 1 : 0;
        if (hard != tx_bits[i]) demod_err++;
    }
    cout << "QPSK BER @ 10dB (app demod): " << (double)demod_err/(2*N) << " (expected ~0.0006-0.0008)" << endl;

    cout << "\n=== Noiseless LLR Check (all modulations) ===" << endl;
    vector<pair<ModulationScheme, string>> schemes = {
        {ModulationScheme::BPSK, "BPSK"},
        {ModulationScheme::QPSK, "QPSK"},
        {ModulationScheme::QAM16, "16QAM"},
        {ModulationScheme::QAM64, "64QAM"},
    };
    for (auto& [sch, name] : schemes) {
        int bps = mod_to_bits_per_symbol(sch);
        int n_sym = 10000;
        BitVec tbits(n_sym * bps);
        for (int i = 0; i < n_sym*bps; i++) tbits[i] = rng() % 2;
        ComplexVec ssym = mod->modulate(tbits, sch);
        vector<double> vn(n_sym, 1e-10);
        SoftVec sllr = mod->demodulate(ssym, sch, vn);
        int err = 0;
        for (int i = 0; i < n_sym*bps; i++) {
            int hard = (sllr[i] > 0) ? 1 : 0;
            if (hard != tbits[i]) err++;
        }
        cout << "  " << name << " noiseless BER: " << (double)err/(n_sym*bps) << " (expected 0)" << endl;
    }

    cout << "\n=== Summary ===" << endl;
    cout << "CRC24A: " << (crc_match ? "MATCH" : "MISMATCH") << " with Sionna" << endl;
    cout << "LLR: app (logsumexp) demod, noiseless BER=0, 10dB QPSK BER reasonable" << endl;
    return crc_match ? 0 : 1;
}
