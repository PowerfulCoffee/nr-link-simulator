#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <cmath>

namespace nr {
namespace phy {

class OfdmModulatorImpl : public IOfdmModulator {
public:
    ComplexVec modulate(const ResourceGrid& grid, double scs) override {
        int n_ant = grid.n_ant;
        int n_sym = grid.n_symbols;
        int n_sc = grid.n_subcarriers;
        int fft_size = get_fft_size(grid.n_subcarriers / 12, scs);
        int cp_len = fft_size / 8;
        
        int samples_per_symbol = fft_size + cp_len;
        int total_samples = n_sym * samples_per_symbol * n_ant;
        
        ComplexVec signal(total_samples, arma::fill::zeros);
        
        for (int ant = 0; ant < n_ant; ant++) {
            for (int sym = 0; sym < n_sym; sym++) {
                ComplexVec freq(fft_size, arma::fill::zeros);
                
                int dc_offset = fft_size / 2;
                int sc_start = dc_offset - n_sc / 2;
                for (int sc = 0; sc < n_sc; sc++) {
                    int fft_idx = (sc_start + sc) % fft_size;
                    freq(fft_idx) = grid.get_re(ant, sym, sc);
                }
                
                ComplexVec time = arma::ifft(freq) * std::sqrt(static_cast<double>(fft_size));
                
                int sample_offset = (ant * n_sym + sym) * samples_per_symbol;
                
                for (int i = 0; i < cp_len; i++) {
                    signal(sample_offset + i) = time(fft_size - cp_len + i);
                }
                for (int i = 0; i < fft_size; i++) {
                    signal(sample_offset + cp_len + i) = time(i);
                }
            }
        }
        
        return signal;
    }
    
    ResourceGrid demodulate(const ComplexVec& signal, int n_ant, double scs, int n_symbols) override {
        int fft_size = 0;
        int n_rb_est = signal.n_elem / (n_symbols * n_ant * 14);
        int n_sc_est = n_rb_est * 12;
        while (fft_size < n_sc_est) fft_size = (fft_size == 0) ? 128 : fft_size * 2;
        if (fft_size == 0) fft_size = 1024;
        
        int n_sc = fft_size;
        for (int rb = 273; rb >= 1; rb--) {
            if (rb * 12 <= fft_size) {
                n_sc = rb * 12;
                break;
            }
        }
        
        int cp_len = fft_size / 8;
        int samples_per_symbol = fft_size + cp_len;
        
        ResourceGrid grid(n_ant, n_symbols, n_sc);
        
        for (int ant = 0; ant < n_ant; ant++) {
            for (int sym = 0; sym < n_symbols; sym++) {
                int sample_offset = (ant * n_symbols + sym) * samples_per_symbol;
                
                ComplexVec time(fft_size);
                for (int i = 0; i < fft_size; i++) {
                    int sig_idx = sample_offset + cp_len + i;
                    if (sig_idx < signal.n_elem) {
                        time(i) = signal(sig_idx);
                    }
                }
                
                ComplexVec freq = arma::fft(time) / std::sqrt(static_cast<double>(fft_size));
                
                int dc_offset = fft_size / 2;
                int sc_start = dc_offset - n_sc / 2;
                for (int sc = 0; sc < n_sc; sc++) {
                    int fft_idx = (sc_start + sc) % fft_size;
                    grid.set_re(ant, sym, sc, freq(fft_idx));
                }
            }
        }
        
        return grid;
    }
};

std::unique_ptr<IOfdmModulator> create_ofdm_modulator() {
    return std::make_unique<OfdmModulatorImpl>();
}

} // namespace phy
} // namespace nr
