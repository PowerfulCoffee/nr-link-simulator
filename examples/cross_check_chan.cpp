#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include "phy/ModuleFactory.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include <random>
#include <vector>
#include <complex>
#include <cstdint>
#include <cstring>

using namespace nr;
using namespace nr::phy;

namespace fs = std::filesystem;

template <typename T>
void write_npy_file(const std::string& filename, const T* data, size_t n_elements,
                    const std::string& dtype_str, const std::vector<size_t>& shape) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    const char magic[] = "\x93NUMPY";
    file.write(magic, 6);

    uint8_t version[2] = {1, 0};
    file.write(reinterpret_cast<const char*>(version), 2);

    std::string shape_str = "(";
    for (size_t i = 0; i < shape.size(); i++) {
        shape_str += std::to_string(shape[i]);
        if (i < shape.size() - 1) shape_str += ", ";
    }
    if (shape.size() == 1) shape_str += ",";
    shape_str += ")";

    bool fortran_order = false;
    std::string header = "{'descr': '" + dtype_str + "', 'fortran_order': " +
                         (fortran_order ? "True" : "False") + ", 'shape': " + shape_str + ", }";

    size_t header_len = header.size() + 1;
    size_t total_header_len = 6 + 2 + 2 + header_len;
    size_t pad = (16 - (total_header_len % 16)) % 16;
    header_len += pad;
    header += std::string(pad, ' ');
    header += '\n';

    uint16_t header_len_le = static_cast<uint16_t>(header.size());
    file.write(reinterpret_cast<const char*>(&header_len_le), 2);
    file.write(header.c_str(), header.size());

    file.write(reinterpret_cast<const char*>(data), n_elements * sizeof(T));
    file.close();
}

void write_uint8_npy(const std::string& filename, const std::vector<uint8_t>& data,
                     const std::vector<size_t>& shape) {
    write_npy_file<uint8_t>(filename, data.data(), data.size(), "|u1", shape);
}

void write_double_npy(const std::string& filename, const std::vector<double>& data,
                      const std::vector<size_t>& shape) {
    write_npy_file<double>(filename, data.data(), data.size(), "<f8", shape);
}

void write_complex_npy(const std::string& filename, const std::vector<std::complex<double>>& data,
                       const std::vector<size_t>& shape) {
    std::vector<double> interleaved;
    interleaved.reserve(data.size() * 2);
    for (auto& c : data) {
        interleaved.push_back(c.real());
        interleaved.push_back(c.imag());
    }
    write_npy_file<double>(filename, interleaved.data(), interleaved.size(), "<c16", shape);
}

void write_grid_npy(const std::string& filename, const ResourceGrid& grid) {
    int n_ant = grid.n_ant;
    int n_sym = grid.n_symbols;
    int n_sc = grid.n_subcarriers;
    std::vector<std::complex<double>> data(n_ant * n_sym * n_sc);
    for (int ant = 0; ant < n_ant; ant++) {
        for (int sym = 0; sym < n_sym; sym++) {
            for (int sc = 0; sc < n_sc; sc++) {
                Complex v = grid.get_re(ant, sym, sc);
                size_t idx = static_cast<size_t>(ant) * n_sym * n_sc +
                             static_cast<size_t>(sym) * n_sc + sc;
                data[idx] = std::complex<double>(v.real(), v.imag());
            }
        }
    }
    write_complex_npy(filename, data, {(size_t)n_ant, (size_t)n_sym, (size_t)n_sc});
}

void write_cube_npy(const std::string& filename, const ComplexCube& cube) {
    int n_sc = cube.n_rows;
    int n_sym = cube.n_cols;
    int n_ch = cube.n_slices;
    std::vector<std::complex<double>> data(n_sc * n_sym * n_ch);
    for (int ch = 0; ch < n_ch; ch++) {
        for (int sym = 0; sym < n_sym; sym++) {
            for (int sc = 0; sc < n_sc; sc++) {
                Complex v = cube(sc, sym, ch);
                size_t idx = static_cast<size_t>(ch) * n_sym * n_sc +
                             static_cast<size_t>(sym) * n_sc + sc;
                data[idx] = std::complex<double>(v.real(), v.imag());
            }
        }
    }
    write_complex_npy(filename, data, {(size_t)n_ch, (size_t)n_sym, (size_t)n_sc});
}

void write_mat_npy(const std::string& filename, const ComplexMat& mat) {
    int n_re = mat.n_rows;
    int n_cols = mat.n_cols;
    std::vector<std::complex<double>> data(n_re * n_cols);
    for (int i = 0; i < n_re; i++) {
        for (int j = 0; j < n_cols; j++) {
            Complex v = mat(i, j);
            data[i * n_cols + j] = std::complex<double>(v.real(), v.imag());
        }
    }
    write_complex_npy(filename, data, {(size_t)n_re, (size_t)n_cols});
}

int main() {
    const int mcs_index = 11;
    const int n_rb = 3;
    const int n_layers = 1;
    const uint64_t random_seed = 42;
    const double es_n0_db = 20.0;
    const int rv = 0;
    const int n_ldpc_iter = 20;

    std::string output_dir = "/workspace/nr-link-simulator/results/cross_check_chan/";
    fs::create_directories(output_dir);

    SimulationConfig config;
    config.mcs_index = mcs_index;
    config.n_rb = n_rb;
    config.n_layers = n_layers;
    config.n_tx_ant = 1;
    config.n_rx_ant = 1;
    config.random_seed = random_seed;
    config.channel_type = ChannelType::AWGN;
    config.dmrs_type = DmrsType::TYPE1;
    config.dmrs_additional_pos = 0;
    config.dmrs_duration = 1;
    config.scs = 15000;
    config.n_ldpc_iterations = n_ldpc_iter;

    int qm = mcs_to_bits_per_symbol(mcs_index);
    double target_coderate = mcs_to_code_rate(mcs_index);
    ModulationScheme mod_scheme = mcs_to_modulation(mcs_index);

    DmrsPattern dmrs_pat = get_dmrs_pattern(DmrsType::TYPE1, 0, 1);
    int n_re_per_prb = 0;
    for (int sym = 0; sym < 14; sym++) {
        if (dmrs_pat.re_per_prb[sym] == 0) {
            n_re_per_prb += 12;
        }
    }

    std::cout << "MCS " << mcs_index << ": Qm=" << qm << ", R=" << target_coderate
              << ", n_re_per_prb=" << n_re_per_prb << std::endl;

    int tbs = calculate_tbs(n_rb, n_re_per_prb, qm, n_layers, target_coderate);
    int G = calculate_num_coded_bits(n_rb, n_re_per_prb, qm, n_layers);

    std::cout << "TBS = " << tbs << " bits, G = " << G << " coded bits" << std::endl;

    TransportBlock tb(tbs);
    tb.mcs = mcs_index;
    tb.rv = rv;
    std::mt19937 rng_tb(random_seed);
    std::uniform_int_distribution<int> bit_dist(0, 1);
    for (int i = 0; i < tbs; i++) {
        tb.bits[i] = static_cast<uint8_t>(bit_dist(rng_tb));
    }

    write_uint8_npy(output_dir + "tb_bits.npy", tb.bits, {static_cast<size_t>(tbs)});

    auto crc_encoder = create_crc_encoder();
    auto ldpc_encoder = create_ldpc_encoder();
    auto ldpc_decoder = create_ldpc_decoder();
    auto rate_matcher = create_rate_matcher();
    auto scrambler = create_scrambler();
    auto modulator = create_modulator();
    auto layer_mapper = create_layer_mapper();
    auto precoder = create_precoder();
    auto dmrs_generator = create_dmrs_generator();
    auto channel_estimator = create_ls_channel_estimator();
    auto equalizer = create_mmse_equalizer();

    int n_sc = n_rb * 12;
    int n_sym = 14;
    int n_ports = 1;

    int tb_crc_len = get_crc_length(tbs);
    BitVec tb_bits_with_crc = crc_encoder->encode(tb.bits, tb_crc_len);

    int k_info = static_cast<int>(tb_bits_with_crc.size());
    CodeBlockSegParams cbs = compute_cb_segmentation(tbs, k_info, n_re_per_prb, n_rb,
                                                      qm, n_layers, target_coderate);

    int E_r;
    if (cbs.num_cb == 1) {
        E_r = G;
    } else {
        E_r = (G + cbs.num_cb - 1) / cbs.num_cb;
        E_r = ((E_r + qm - 1) / qm) * qm;
    }

    BitVec rm_bits_all;
    for (int c = 0; c < cbs.num_cb; c++) {
        BitVec cb_bits;
        int k_cb;
        if (cbs.num_cb == 1) {
            cb_bits = tb_bits_with_crc;
            k_cb = k_info;
        } else {
            int payload_per_cb = cbs.cb_info_bits;
            cb_bits.resize(payload_per_cb, 0);
            int start_bit = c * payload_per_cb;
            int end_bit = std::min((c + 1) * payload_per_cb, k_info);
            int copy_len = end_bit - start_bit;
            for (int i = 0; i < copy_len; i++) {
                cb_bits[i] = tb_bits_with_crc[start_bit + i];
            }
            for (int i = copy_len; i < payload_per_cb; i++) {
                cb_bits[i] = 0;
            }
            BitVec cb_with_crc = crc_encoder->encode(cb_bits, cbs.cb_crc_len);
            cb_bits = cb_with_crc;
            k_cb = static_cast<int>(cb_bits.size());
            cb_bits.resize(cbs.cb_k, 0);
        }
        BitVec cb_bits_for_encoder = cb_bits;
        if (cbs.num_cb == 1) cb_bits_for_encoder.resize(cbs.cb_k, 0);
        BitVec coded_bits_full = ldpc_encoder->encode(cb_bits_for_encoder, cbs.bgn, cbs.zc);
        int n_punctured = 2 * cbs.zc;
        BitVec coded_bits(coded_bits_full.begin() + n_punctured, coded_bits_full.end());
        int n_filler = cbs.cb_k - k_cb;
        int filler_start = k_cb - n_punctured;
        BitVec coded_bits_comp;
        coded_bits_comp.reserve(coded_bits.size() - n_filler);
        if (n_filler > 0 && filler_start >= 0 && filler_start + n_filler <= (int)coded_bits.size()) {
            coded_bits_comp.insert(coded_bits_comp.end(), coded_bits.begin(), coded_bits.begin() + filler_start);
            coded_bits_comp.insert(coded_bits_comp.end(), coded_bits.begin() + filler_start + n_filler, coded_bits.end());
        } else {
            coded_bits_comp = coded_bits;
        }
        int this_E = E_r;
        if (c == cbs.num_cb - 1 && cbs.num_cb > 1) {
            this_E = G - (cbs.num_cb - 1) * E_r;
            this_E = std::max(this_E, qm);
        }
        BitVec rm_bits_cb = rate_matcher->rate_match(coded_bits_comp, this_E, rv, cbs.bgn, cbs.zc, qm, 0);
        rm_bits_all.insert(rm_bits_all.end(), rm_bits_cb.begin(), rm_bits_cb.end());
    }
    if ((int)rm_bits_all.size() > G) rm_bits_all.resize(G);

    uint32_t scrambling_seed = 1;
    BitVec scrambled_bits = scrambler->scramble(rm_bits_all, scrambling_seed);
    ComplexVec modulated = modulator->modulate(scrambled_bits, mod_scheme);
    ComplexMat layered = layer_mapper->map(modulated, n_layers);
    ComplexMat w(n_ports, n_layers, arma::fill::zeros);
    for (int l = 0; l < n_layers; l++) w(l, l) = Complex(1.0, 0.0);
    ComplexMat precoded = precoder->precode(layered, w);

    ResourceGrid tx_grid(n_ports, n_sym, n_sc);
    tx_grid.reset();
    dmrs_generator->generate_dmrs(config, tx_grid, 0, 0);

    ResourceGrid dmrs_only_grid = tx_grid;

    int re_idx = 0;
    for (int sym = 0; sym < n_sym; sym++) {
        if (dmrs_pat.re_per_prb[sym] > 0) continue;
        for (int rb = 0; rb < n_rb; rb++) {
            for (int sc_in_rb = 0; sc_in_rb < 12; sc_in_rb++) {
                int sc = rb * 12 + sc_in_rb;
                if (sc >= n_sc || re_idx >= (int)precoded.n_rows) continue;
                for (int ant = 0; ant < n_ports; ant++) {
                    tx_grid.set_re(ant, sym, sc, precoded(re_idx, ant));
                }
                re_idx++;
            }
        }
    }

    write_grid_npy(output_dir + "tx_grid.npy", tx_grid);

    double sinr_lin = std::pow(10.0, es_n0_db / 10.0);
    double noise_var = 1.0 / sinr_lin;
    double sigma_dim = std::sqrt(noise_var / 2.0);

    std::mt19937 noise_rng(random_seed + 1000);
    std::normal_distribution<double> nd(0.0, sigma_dim);

    ResourceGrid rx_grid(n_ports, n_sym, n_sc);
    rx_grid.reset();
    for (int ant = 0; ant < n_ports; ant++) {
        for (int sym = 0; sym < n_sym; sym++) {
            for (int sc = 0; sc < n_sc; sc++) {
                Complex val = tx_grid.get_re(ant, sym, sc);
                double nr = nd(noise_rng);
                double ni = nd(noise_rng);
                rx_grid.set_re(ant, sym, sc, val + Complex(nr, ni));
            }
        }
    }

    write_grid_npy(output_dir + "rx_grid.npy", rx_grid);

    ComplexCube h_est = channel_estimator->estimate(rx_grid, dmrs_only_grid, config);
    write_cube_npy(output_dir + "h_est.npy", h_est);

    int n_pdsch_re = 0;
    for (int sym = 0; sym < n_sym; sym++) {
        if (dmrs_pat.re_per_prb[sym] == 0) n_pdsch_re += n_sc;
    }

    ComplexMat rx_pdsch(n_pdsch_re, n_ports, arma::fill::zeros);
    ComplexCube h_at_pdsch(n_pdsch_re, 1, n_ports * n_layers, arma::fill::zeros);
    int pdsch_idx = 0;
    for (int sym = 0; sym < n_sym; sym++) {
        if (dmrs_pat.re_per_prb[sym] > 0) continue;
        for (int sc = 0; sc < n_sc; sc++) {
            for (int ant = 0; ant < n_ports; ant++) {
                rx_pdsch(pdsch_idx, ant) = rx_grid.get_re(ant, sym, sc);
            }
            for (int rx = 0; rx < n_ports; rx++) {
                for (int l = 0; l < n_layers; l++) {
                    int ch_idx = rx * n_layers + l;
                    h_at_pdsch(pdsch_idx, 0, ch_idx) = h_est(sc, sym, ch_idx);
                }
            }
            pdsch_idx++;
        }
    }

    write_mat_npy(output_dir + "rx_pdsch.npy", rx_pdsch);
    write_cube_npy(output_dir + "h_at_pdsch.npy", h_at_pdsch);

    ComplexMat equalized = equalizer->equalize(rx_pdsch, h_at_pdsch, noise_var, n_layers);
    write_mat_npy(output_dir + "equalized.npy", equalized);

    std::vector<double> eff_noise = equalizer->get_eff_noise_var();
    write_double_npy(output_dir + "eff_noise_var.npy", eff_noise, {eff_noise.size()});

    ComplexVec demapped(n_pdsch_re * n_layers);
    for (int i = 0; i < n_pdsch_re; i++) {
        for (int l = 0; l < n_layers; l++) {
            demapped(i * n_layers + l) = equalized(i, l);
        }
    }

    std::vector<double> demapped_noise(demapped.n_elem, noise_var);
    for (int i = 0; i < n_pdsch_re; i++) {
        for (int l = 0; l < n_layers; l++) {
            int out_idx = i * n_layers + l;
            int eff_idx = i * n_layers + l;
            if (eff_idx < (int)eff_noise.size()) {
                demapped_noise[out_idx] = (eff_noise[eff_idx] > 1e-12) ? eff_noise[eff_idx] : noise_var;
            }
        }
    }

    SoftVec llr = modulator->demodulate(demapped, mod_scheme, demapped_noise);
    write_double_npy(output_dir + "llr_values.npy", llr, {llr.size()});

    SoftVec descrambled = scrambler->descramble(llr, scrambling_seed);

    int C = cbs.num_cb;
    int K = cbs.cb_k;
    int zc = cbs.zc;
    int bgn = cbs.bgn;
    int N = (bgn == 1) ? 68 * zc : 52 * zc;
    int n_punctured = 2 * zc;
    int N_cb = N - n_punctured;

    bool crc_ok = true;
    BitVec decoded_info_bits;
    std::vector<BitVec> cb_decoded(C);
    bool all_cb_ok = true;
    std::vector<CodeBlockInfo> cb_info_list;
    int cb_offset = 0;
    for (int c = 0; c < C; c++) {
        int this_E = E_r;
        if (c == C - 1 && C > 1) this_E = G - (C - 1) * E_r;
        CodeBlockInfo cbi;
        cbi.offset = cb_offset;
        cbi.length = this_E;
        cbi.e_bits = this_E;
        cb_info_list.push_back(cbi);
        cb_offset += this_E;
    }

    for (int c = 0; c < C; c++) {
        int cb_llr_start = cb_info_list[c].offset;
        int E_cb = cb_info_list[c].e_bits;
        SoftVec cb_llr(descrambled.begin() + cb_llr_start,
                       descrambled.begin() + cb_llr_start + E_cb);
        int cb_k_info = (C == 1) ? tbs + tb_crc_len : cbs.cb_size_with_crc;
        int cb_crc_len = cbs.cb_crc_len;
        int n_filler = K - cb_k_info;
        int filler_start_comp = cb_k_info - n_punctured;
        int n_cb_comp = N_cb - n_filler;
        SoftVec recovered_comp = rate_matcher->rate_recover(cb_llr, n_cb_comp, rv, bgn, zc, qm, 0);
        SoftVec full_llr(N, 0.0);
        if (n_filler > 0 && filler_start_comp >= 0 && filler_start_comp <= n_cb_comp) {
            for (int i = 0; i < filler_start_comp; i++) full_llr[n_punctured + i] = recovered_comp[i];
            for (int i = 0; i < n_filler; i++) full_llr[cb_k_info + i] = 20.0;
            for (int i = filler_start_comp; i < n_cb_comp; i++) {
                int dst_idx = cb_k_info + n_filler + (i - filler_start_comp);
                if (dst_idx < N) full_llr[dst_idx] = recovered_comp[i];
            }
        } else {
            for (int i = 0; i < n_cb_comp && i < N_cb; i++) full_llr[n_punctured + i] = recovered_comp[i];
        }
        for (int i = 0; i < N; i++) {
            if (full_llr[i] > 20.0) full_llr[i] = 20.0;
            else if (full_llr[i] < -20.0) full_llr[i] = -20.0;
        }
        auto [decoded_cb_bits, conv_ok] = ldpc_decoder->decode(full_llr, bgn, zc, n_ldpc_iter, true);
        (void)conv_ok;
        if ((int)decoded_cb_bits.size() < cb_k_info) { all_cb_ok = false; cb_decoded[c] = BitVec(cb_k_info, 0); continue; }
        BitVec cb_bits_out(decoded_cb_bits.begin(), decoded_cb_bits.begin() + cb_k_info);
        if (C > 1 && cb_crc_len > 0) {
            auto [cb_info_bits_out, cb_check_ok] = crc_encoder->decode(cb_bits_out, cb_crc_len);
            if (!cb_check_ok) all_cb_ok = false;
            cb_decoded[c] = cb_info_bits_out;
        } else {
            cb_decoded[c] = cb_bits_out;
        }
    }

    if (C == 1) {
        int total_len = tbs + tb_crc_len;
        if ((int)cb_decoded[0].size() < total_len) crc_ok = false;
        else {
            BitVec rx_bits(cb_decoded[0].begin(), cb_decoded[0].begin() + total_len);
            auto [info_bits, tb_check_ok] = crc_encoder->decode(rx_bits, tb_crc_len);
            decoded_info_bits = info_bits;
            crc_ok = tb_check_ok;
        }
    } else {
        BitVec tb_bits_with_crc_decoded;
        for (int c = 0; c < C; c++) tb_bits_with_crc_decoded.insert(tb_bits_with_crc_decoded.end(), cb_decoded[c].begin(), cb_decoded[c].end());
        int total_tb_len = tbs + 24;
        if ((int)tb_bits_with_crc_decoded.size() < total_tb_len) crc_ok = false;
        else {
            BitVec rx_bits(tb_bits_with_crc_decoded.begin(), tb_bits_with_crc_decoded.begin() + total_tb_len);
            auto [info_bits, tb_check_ok] = crc_encoder->decode(rx_bits, 24);
            decoded_info_bits = info_bits;
            crc_ok = tb_check_ok && all_cb_ok;
        }
    }

    write_uint8_npy(output_dir + "decoded_bits.npy", decoded_info_bits, {static_cast<size_t>(tbs)});

    {
        std::ofstream info_file(output_dir + "config.txt");
        info_file << "mcs_index=" << mcs_index << "\n";
        info_file << "n_rb=" << n_rb << "\n";
        info_file << "n_layers=" << n_layers << "\n";
        info_file << "n_sc=" << n_sc << "\n";
        info_file << "n_sym=" << n_sym << "\n";
        info_file << "qm=" << qm << "\n";
        info_file << "tbs=" << tbs << "\n";
        info_file << "G=" << G << "\n";
        info_file << "es_n0_db=" << es_n0_db << "\n";
        info_file << "noise_var=" << noise_var << "\n";
        info_file << "n_pdsch_re=" << n_pdsch_re << "\n";
        info_file << "crc_ok=" << (crc_ok ? "PASS" : "FAIL") << "\n";
        info_file.close();
    }

    int bit_errors = 0;
    if ((int)decoded_info_bits.size() == tbs) {
        for (int i = 0; i < tbs; i++) {
            if (decoded_info_bits[i] != tb.bits[i]) bit_errors++;
        }
    } else {
        bit_errors = tbs;
    }

    std::cout << "Processing complete." << std::endl;
    std::cout << "CRC check: " << (crc_ok ? "PASS" : "FAIL") << std::endl;
    std::cout << "Bit errors: " << bit_errors << "/" << tbs << std::endl;
    std::cout << "Output files written to: " << output_dir << std::endl;

    return 0;
}
