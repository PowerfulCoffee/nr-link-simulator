#include "phy/PhyInterfaces.h"
#include "phy/ModuleFactory.h"
#include "common/NrTables.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <vector>
#include <fstream>

using namespace nr;
using namespace nr::phy;
using namespace std;

int main() {
    cout << "=== CRC Cross-Check ===" << endl;
    auto crc_enc = create_crc_encoder();

    vector<vector<uint8_t>> test_cases = {
        {1},
        {0},
        {1,0,1,0,1,1,0,0,1,1},
    };
    for (int i = 0; i < 100; i++) {
        vector<uint8_t> bits(100);
        for (int j = 0; j < 100; j++) bits[j] = rand() % 2;
        test_cases.push_back(bits);
    }

    bool crc_ok = true;
    for (size_t tc = 0; tc < test_cases.size(); tc++) {
        auto& bits = test_cases[tc];
        auto encoded = crc_enc->encode(BitVec(bits.begin(), bits.end()), 24);
        auto [decoded, valid] = crc_enc->decode(encoded, 24);
        if (!valid) {
            cout << "  CRC test " << tc << " FAILED: self-check not ok" << endl;
            crc_ok = false;
            break;
        }
        if (decoded.size() != bits.size()) {
            cout << "  CRC test " << tc << " FAILED: size mismatch" << endl;
            crc_ok = false;
            break;
        }
        for (size_t i = 0; i < bits.size(); i++) {
            if (decoded[i] != bits[i]) {
                cout << "  CRC test " << tc << " FAILED: bit mismatch at " << i << endl;
                crc_ok = false;
                break;
            }
        }
        if (!crc_ok) break;
    }
    if (crc_ok) {
        cout << "  CRC self-check PASSED for all test cases" << endl;
    }

    ofstream crc_out("/workspace/crc_test_vectors.txt");
    for (int tc = 0; tc < 10; tc++) {
        vector<uint8_t> bits(50);
        for (int j = 0; j < 50; j++) bits[j] = rand() % 2;
        auto encoded = crc_enc->encode(BitVec(bits.begin(), bits.end()), 24);
        for (auto b : encoded) crc_out << (int)b;
        crc_out << endl;
    }
    crc_out.close();
    cout << "  CRC test vectors written to /workspace/crc_test_vectors.txt" << endl;

    cout << "\n=== LLR Demod Cross-Check ===" << endl;
    auto mod = create_modulator();

    const int NSYM = 10000;
    mt19937 rng(42);
    normal_distribution<double> nd(0.0, 1.0);
    uniform_int_distribution<int> uid(0, 1);

    vector<ModulationScheme> schemes = {
        ModulationScheme::QPSK, ModulationScheme::QAM16,
        ModulationScheme::QAM64, ModulationScheme::QAM256
    };
    vector<string> scheme_names = {"QPSK", "16QAM", "64QAM", "256QAM"};

    ofstream llr_out("/workspace/llr_test_vectors.txt");
    llr_out << scientific << setprecision(15);
    for (size_t s = 0; s < schemes.size(); s++) {
        int bps = mod_to_bits_per_symbol(schemes[s]);
        int n_sym = NSYM / bps;
        int total_bits = n_sym * bps;
        BitVec bits(total_bits);
        for (int i = 0; i < total_bits; i++) bits[i] = uid(rng);
        ComplexVec sym = mod->modulate(bits, schemes[s]);

        vector<double> noise_vars(n_sym);
        for (int i = 0; i < n_sym; i++) {
            noise_vars[i] = 0.1 + 0.9 * (i / (double)n_sym);
        }

        ComplexVec noisy_sym(n_sym);
        for (int i = 0; i < n_sym; i++) {
            double std = sqrt(noise_vars[i] / 2.0);
            noisy_sym(i) = Complex(sym(i).real() + std*nd(rng), sym(i).imag() + std*nd(rng));
        }

        SoftVec llr = mod->demodulate(noisy_sym, schemes[s], noise_vars);

        double err_count = 0;
        BitVec hard_bits(total_bits);
        for (int i = 0; i < total_bits; i++) {
            hard_bits[i] = (llr[i] < 0) ? 1 : 0;
            if (hard_bits[i] != bits[i]) err_count++;
        }
        double ber = err_count / total_bits;
        cout << "  " << scheme_names[s] << ": LLR BER = " << ber << " (at SNR ~5-10dB)" << endl;

        llr_out << scheme_names[s] << endl;
        for (int i = 0; i < n_sym; i++) {
            llr_out << noisy_sym(i).real() << " " << noisy_sym(i).imag() << " " << noise_vars[i];
            for (int b = 0; b < bps; b++) {
                llr_out << " " << llr[i*bps + b];
            }
            llr_out << endl;
        }
    }
    llr_out.close();
    cout << "  LLR test vectors written to /workspace/llr_test_vectors.txt" << endl;

    cout << "\n=== LDPC Decode Cross-Check ===" << endl;
    auto ldpc_enc = create_ldpc_encoder();
    auto ldpc_dec = create_ldpc_decoder();

    const int zc = 24;
    const int bgn = 2;
    struct LDPCDims { int rows, cols, n_sys, N, K; };
    auto get_ldpc_dims = [] (int bgn, int zc) -> LDPCDims {
        LDPCDims d;
        if (bgn == 1) { d.rows=46; d.cols=68; } else { d.rows=42; d.cols=52; }
        d.n_sys = d.cols - d.rows; d.N = d.cols*zc; d.K = d.n_sys*zc;
        return d;
    };
    auto ldpc_d = get_ldpc_dims(bgn, zc);
    int K_ldpc = ldpc_d.K;
    int N_ldpc = ldpc_d.N;

    BitVec info_bits(K_ldpc - 24);
    for (int i = 0; i < (int)info_bits.size(); i++) info_bits[i] = uid(rng);
    auto info_with_crc = crc_enc->encode(info_bits, 24);

    BitVec enc_in(K_ldpc, 0);
    for (int i = 0; i < (int)info_with_crc.size() && i < K_ldpc; i++) enc_in[i] = info_with_crc[i];

    BitVec cw = ldpc_enc->encode(enc_in, bgn, zc);

    ofstream ldpc_out("/workspace/ldpc_test_vectors.txt");
    ldpc_out << scientific << setprecision(15);
    ldpc_out << bgn << " " << zc << endl;

    for (int snr_idx = 0; snr_idx < 5; snr_idx++) {
        double es_n0_db = -2 + snr_idx * 2;
        double es_n0_lin = pow(10.0, es_n0_db / 10.0);
        double noise_var = 1.0 / es_n0_lin;
        double std = sqrt(noise_var / 2.0);

        SoftVec llr_full(N_ldpc);
        for (int i = 0; i < N_ldpc; i++) {
            double tx_bpsk = 1.0 - 2.0 * cw[i];
            double rx_re = tx_bpsk + std * nd(rng);
            llr_full[i] = 2.0 * rx_re / (noise_var/2.0);
        }

        auto [decoded, conv_ok] = ldpc_dec->decode(llr_full, bgn, zc, 20, true);

        int bit_errs = 0;
        int min_len = min((int)decoded.size(), K_ldpc);
        for (int i = 0; i < min_len; i++) {
            if (decoded[i] != enc_in[i]) bit_errs++;
        }
        bool crc_pass = false;
        if (bit_errs == 0) {
            BitVec dec_info(decoded.begin(), decoded.begin() + info_with_crc.size());
            auto [_, crc_ok_dec] = crc_enc->decode(dec_info, 24);
            crc_pass = crc_ok_dec;
        }
        cout << "  SNR=" << es_n0_db << "dB: conv=" << conv_ok << " bit_errs=" << bit_errs << " crc_pass=" << crc_pass << endl;

        ldpc_out << es_n0_db << " " << noise_var << " " << conv_ok << " " << bit_errs;
        for (int i = 0; i < N_ldpc; i++) ldpc_out << " " << llr_full[i];
        ldpc_out << endl;
    }
    ldpc_out.close();
    cout << "  LDPC test vectors written to /workspace/ldpc_test_vectors.txt" << endl;

    cout << "\nAll C++ test vectors written. Please run Python verification." << endl;
    return 0;
}
