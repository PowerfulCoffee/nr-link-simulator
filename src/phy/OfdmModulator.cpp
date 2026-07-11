#include "phy/PhyInterfaces.h"
#include <cmath>
#include <algorithm>

namespace nr {
namespace phy {

namespace {

int get_fft_size(int n_subcarriers) {
    int fft_size = 1;
    while (fft_size < n_subcarriers + 1) {
        fft_size <<= 1;
    }
    return fft_size;
}

int get_cp_len(int fft_size, bool is_first_symbol) {
    int cp_short = 144 * fft_size / 2048;
    if (is_first_symbol) {
        return cp_short + 16 * fft_size / 2048;
    }
    return cp_short;
}

int estimate_fft_size(int total_samples, int n_ant, int n_symbols) {
    double ratio = static_cast<double>(total_samples) * 2048.0 /
                   static_cast<double>(n_ant * (n_symbols * 2192 + 16));
    int approx = static_cast<int>(std::round(ratio));
    int fft_size = 1;
    while (fft_size < approx) {
        fft_size <<= 1;
    }
    int fft_size_lower = fft_size >> 1;
    auto expected_len = [&](int fft) {
        int cp_l = get_cp_len(fft, true);
        int cp_s = get_cp_len(fft, false);
        return n_ant * (n_symbols * fft + cp_l + (n_symbols - 1) * cp_s);
    };
    if (fft_size_lower >= 128) {
        int diff_upper = std::abs(expected_len(fft_size) - total_samples);
        int diff_lower = std::abs(expected_len(fft_size_lower) - total_samples);
        if (diff_lower < diff_upper) {
            fft_size = fft_size_lower;
        }
    }
    return fft_size;
}

}

class OfdmModulatorImpl : public IOfdmModulator {
public:
    ComplexVec modulate(const ResourceGrid& grid, double scs) override {
        int n_ant = grid.n_ant;
        int n_sym = grid.n_symbols;
        int n_sc = grid.n_subcarriers;
        int fft_size = get_fft_size(n_sc);

        std::vector<int> cp_lens(n_sym);
        int total_per_ant = 0;
        for (int sym = 0; sym < n_sym; sym++) {
            cp_lens[sym] = get_cp_len(fft_size, sym == 0);
            total_per_ant += fft_size + cp_lens[sym];
        }
        int total_samples = n_ant * total_per_ant;

        ComplexVec signal(total_samples, arma::fill::zeros);

        for (int ant = 0; ant < n_ant; ant++) {
            int sample_offset = ant * total_per_ant;
            for (int sym = 0; sym < n_sym; sym++) {
                int cp_len = cp_lens[sym];
                ComplexVec freq(fft_size, arma::fill::zeros);

                int dc_pos = fft_size / 2;
                int half_sc = n_sc / 2;
                for (int sc = 0; sc < half_sc; sc++) {
                    freq(dc_pos - half_sc + sc) = grid.get_re(ant, sym, sc);
                }
                for (int sc = half_sc; sc < n_sc; sc++) {
                    freq(dc_pos + 1 + (sc - half_sc)) = grid.get_re(ant, sym, sc);
                }
                freq(dc_pos) = Complex(0.0, 0.0);

                ComplexVec time = arma::ifft(freq) * std::sqrt(static_cast<double>(fft_size));

                for (int i = 0; i < cp_len; i++) {
                    signal(sample_offset + i) = time(fft_size - cp_len + i);
                }
                for (int i = 0; i < fft_size; i++) {
                    signal(sample_offset + cp_len + i) = time(i);
                }

                sample_offset += fft_size + cp_len;
            }
        }

        return signal;
    }

    ResourceGrid demodulate(const ComplexVec& signal, int n_ant, double scs, int n_symbols) override {
        int fft_size = estimate_fft_size(static_cast<int>(signal.n_elem), n_ant, n_symbols);

        std::vector<int> cp_lens(n_symbols);
        int total_per_ant = 0;
        for (int sym = 0; sym < n_symbols; sym++) {
            cp_lens[sym] = get_cp_len(fft_size, sym == 0);
            total_per_ant += fft_size + cp_lens[sym];
        }

        int n_sc = 0;
        int max_rb = (fft_size - 1) / 12;
        for (int rb = max_rb; rb >= 1; rb--) {
            int test_sc = rb * 12;
            if (get_fft_size(test_sc) == fft_size) {
                n_sc = test_sc;
                break;
            }
        }
        if (n_sc == 0) {
            n_sc = fft_size - 1;
            if (n_sc % 12 != 0) {
                n_sc = (n_sc / 12) * 12;
            }
        }

        ResourceGrid grid(n_ant, n_symbols, n_sc);

        for (int ant = 0; ant < n_ant; ant++) {
            int sample_offset = ant * total_per_ant;
            for (int sym = 0; sym < n_symbols; sym++) {
                int cp_len = cp_lens[sym];

                ComplexVec time(fft_size);
                for (int i = 0; i < fft_size; i++) {
                    int sig_idx = sample_offset + cp_len + i;
                    if (sig_idx < static_cast<int>(signal.n_elem)) {
                        time(i) = signal(sig_idx);
                    } else {
                        time(i) = Complex(0.0, 0.0);
                    }
                }

                ComplexVec freq = arma::fft(time) / std::sqrt(static_cast<double>(fft_size));

                int dc_pos = fft_size / 2;
                int half_sc = n_sc / 2;
                for (int sc = 0; sc < half_sc; sc++) {
                    grid.set_re(ant, sym, sc, freq(dc_pos - half_sc + sc));
                }
                for (int sc = half_sc; sc < n_sc; sc++) {
                    grid.set_re(ant, sym, sc, freq(dc_pos + 1 + (sc - half_sc)));
                }

                sample_offset += fft_size + cp_len;
            }
        }

        return grid;
    }
};

std::unique_ptr<IOfdmModulator> create_ofdm_modulator() {
    return std::make_unique<OfdmModulatorImpl>();
}

}
}
