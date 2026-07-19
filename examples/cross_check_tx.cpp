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

constexpr double LLR_CLIP = 20.0;

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
    std::vector<size_t> cshape = shape;
    if (cshape.size() >= 1) {
    }
    write_npy_file<double>(filename, interleaved.data(), interleaved.size(), "<c16", shape);
}

int main() {
    const int mcs_index = 11;
    const int n_rb = 3;
    const int n_layers = 1;
    const uint64_t random_seed = 42;
    const double es_n0_db = 10.0;
    const int rv = 0;
    const int n_ldpc_iter = 20;

    std::string output_dir = "/workspace/nr-link-simulator/results/cross_check/";
    fs::create_directories(output_dir);

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

    int tb_crc_len = get_crc_length(tbs);
    BitVec tb_bits_with_crc = crc_encoder->encode(tb.bits, tb_crc_len);

    write_uint8_npy(output_dir + "tb_crc_bits.npy", tb_bits_with_crc, {tb_bits_with_crc.size()});

    int k_info = static_cast<int>(tb_bits_with_crc.size());
    CodeBlockSegParams cbs = compute_cb_segmentation(tbs, k_info, n_re_per_prb, n_rb,
                                                      qm, n_layers, target_coderate);

    int N = (cbs.bgn == 1) ? 68 * cbs.zc : 52 * cbs.zc;
    int n_punctured = 2 * cbs.zc;

    int E_r;
    if (cbs.num_cb == 1) {
        E_r = G;
    } else {
        E_r = (G + cbs.num_cb - 1) / cbs.num_cb;
        E_r = ((E_r + qm - 1) / qm) * qm;
    }

    int cb0_e_bits = E_r;
    if (cbs.num_cb > 1) {
        cb0_e_bits = E_r;
    }

    std::cout << "CB segmentation: bgn=" << cbs.bgn << ", Zc=" << cbs.zc
              << ", Kb=" << cbs.k_b << ", C=" << cbs.num_cb
              << ", cb_info_bits=" << cbs.cb_info_bits
              << ", cb_size_with_crc=" << cbs.cb_size_with_crc
              << ", K=" << cbs.cb_k << ", E=" << E_r << std::endl;

    {
        std::ofstream info_file(output_dir + "cb_info.txt");
        info_file << "bgn=" << cbs.bgn << "\n";
        info_file << "zc=" << cbs.zc << "\n";
        info_file << "k_b=" << cbs.k_b << "\n";
        info_file << "num_cb=" << cbs.num_cb << "\n";
        info_file << "cb_info_bits=" << cbs.cb_info_bits << "\n";
        info_file << "cb_size_with_crc=" << cbs.cb_size_with_crc << "\n";
        info_file << "cb_k=" << cbs.cb_k << "\n";
        info_file << "cb_e_bits=" << E_r << "\n";
        info_file << "G=" << G << "\n";
        info_file << "qm=" << qm << "\n";
        info_file << "crc_len=" << tb_crc_len << "\n";
        info_file << "N=" << N << "\n";
        info_file << "n_punctured=" << n_punctured << "\n";
        info_file << "tb_size=" << tbs << "\n";
        info_file << "mcs_index=" << mcs_index << "\n";
        info_file << "n_rb=" << n_rb << "\n";
        info_file << "n_re_per_prb=" << n_re_per_prb << "\n";
        info_file << "n_layers=" << n_layers << "\n";
        info_file << "es_n0_db=" << es_n0_db << "\n";
        info_file.close();
    }

    BitVec cb0_bits_with_crc;
    BitVec cb0_ldpc_input;
    BitVec cb0_ldpc_output;
    BitVec cb0_punctured;
    BitVec cb0_no_filler;
    BitVec cb0_rate_matched;
    BitVec rm_bits_all;
    std::vector<CodeBlockInfo> cb_info_list;

    int cb_offset = 0;
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

        if (c == 0) {
            cb0_bits_with_crc.resize(k_cb);
            for (int i = 0; i < k_cb; i++) {
                cb0_bits_with_crc[i] = cb_bits[i];
            }
            cb0_ldpc_input.assign(cbs.cb_k, 0);
            for (int i = 0; i < k_cb; i++) {
                cb0_ldpc_input[i] = cb_bits[i];
            }
        }

        BitVec cb_bits_for_encoder = cb_bits;
        if (cbs.num_cb == 1) {
            cb_bits_for_encoder.resize(cbs.cb_k, 0);
        }

        BitVec coded_bits_full = ldpc_encoder->encode(cb_bits_for_encoder, cbs.bgn, cbs.zc);

        if (c == 0) {
            cb0_ldpc_output = coded_bits_full;
        }

        BitVec coded_bits(coded_bits_full.begin() + n_punctured, coded_bits_full.end());

        if (c == 0) {
            cb0_punctured = coded_bits;
        }

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

        if (c == 0) {
            cb0_no_filler = coded_bits_comp;
        }

        int this_E = E_r;
        if (c == cbs.num_cb - 1 && cbs.num_cb > 1) {
            this_E = G - (cbs.num_cb - 1) * E_r;
            this_E = std::max(this_E, qm);
        }

        BitVec rm_bits_cb = rate_matcher->rate_match(coded_bits_comp, this_E, rv, cbs.bgn, cbs.zc, qm, 0);

        if (c == 0) {
            cb0_rate_matched = rm_bits_cb;
        }

        CodeBlockInfo cbi;
        cbi.offset = cb_offset;
        cbi.length = this_E;
        cbi.e_bits = this_E;
        cb_info_list.push_back(cbi);

        rm_bits_all.insert(rm_bits_all.end(), rm_bits_cb.begin(), rm_bits_cb.end());
        cb_offset += this_E;
    }

    if ((int)rm_bits_all.size() > G) {
        rm_bits_all.resize(G);
    }

    write_uint8_npy(output_dir + "cb0_bits_with_crc.npy", cb0_bits_with_crc, {cb0_bits_with_crc.size()});
    write_uint8_npy(output_dir + "cb0_ldpc_input.npy", cb0_ldpc_input, {static_cast<size_t>(cbs.cb_k)});
    write_uint8_npy(output_dir + "cb0_ldpc_output.npy", cb0_ldpc_output, {static_cast<size_t>(N)});
    write_uint8_npy(output_dir + "cb0_punctured.npy", cb0_punctured, {cb0_punctured.size()});
    write_uint8_npy(output_dir + "cb0_no_filler.npy", cb0_no_filler, {cb0_no_filler.size()});
    write_uint8_npy(output_dir + "cb0_rate_matched.npy", cb0_rate_matched, {static_cast<size_t>(cb0_e_bits)});
    write_uint8_npy(output_dir + "rm_bits_all.npy", rm_bits_all, {static_cast<size_t>(G)});

    uint32_t scrambling_seed = 1;
    BitVec scrambled_bits = scrambler->scramble(rm_bits_all, scrambling_seed);
    write_uint8_npy(output_dir + "scrambled_bits.npy", scrambled_bits, {static_cast<size_t>(G)});

    ComplexVec modulated = modulator->modulate(scrambled_bits, mod_scheme);

    std::vector<std::complex<double>> mod_sym_std(modulated.n_elem);
    for (size_t i = 0; i < modulated.n_elem; i++) {
        mod_sym_std[i] = std::complex<double>(modulated(i).real(), modulated(i).imag());
    }
    write_complex_npy(output_dir + "modulated_symbols.npy", mod_sym_std, {modulated.n_elem});

    double sinr_lin = std::pow(10.0, es_n0_db / 10.0);
    double noise_var = 1.0 / sinr_lin;
    double sigma_dim = std::sqrt(noise_var / 2.0);

    std::mt19937 noise_rng(random_seed + 1000);
    std::normal_distribution<double> nd(0.0, sigma_dim);

    ComplexVec rx_symbols(modulated.n_elem);
    for (size_t i = 0; i < modulated.n_elem; i++) {
        double nr = nd(noise_rng);
        double ni = nd(noise_rng);
        rx_symbols(i) = modulated(i) + Complex(nr, ni);
    }

    std::vector<std::complex<double>> rx_sym_std(rx_symbols.n_elem);
    for (size_t i = 0; i < rx_symbols.n_elem; i++) {
        rx_sym_std[i] = std::complex<double>(rx_symbols(i).real(), rx_symbols(i).imag());
    }
    write_complex_npy(output_dir + "rx_symbols.npy", rx_sym_std, {rx_symbols.n_elem});

    SoftVec llr_values = modulator->demodulate(rx_symbols, mod_scheme, noise_var);
    for (size_t i = 0; i < llr_values.size(); i++) {
        if (llr_values[i] > LLR_CLIP) llr_values[i] = LLR_CLIP;
        else if (llr_values[i] < -LLR_CLIP) llr_values[i] = -LLR_CLIP;
    }

    write_double_npy(output_dir + "llr_values.npy", llr_values, {llr_values.size()});

    SoftVec descrambled_llr = scrambler->descramble(llr_values, scrambling_seed);
    write_double_npy(output_dir + "descrambled_llr.npy", descrambled_llr, {descrambled_llr.size()});

    int C = cbs.num_cb;
    int K = cbs.cb_k;
    int zc = cbs.zc;
    int bgn = cbs.bgn;
    int N_cb = N - n_punctured;

    std::vector<double> llr_for_decoder_all;
    BitVec decoded_tb_bits;
    bool crc_ok = true;

    std::vector<BitVec> cb_decoded(C);
    bool all_cb_ok = true;

    for (int c = 0; c < C; c++) {
        int cb_llr_start = cb_info_list[c].offset;
        int E_cb = cb_info_list[c].e_bits;

        SoftVec cb_llr(descrambled_llr.begin() + cb_llr_start,
                       descrambled_llr.begin() + cb_llr_start + E_cb);

        int cb_k_info;
        int cb_crc_len = cbs.cb_crc_len;

        if (C == 1) {
            cb_k_info = tbs + tb_crc_len;
        } else {
            cb_k_info = cbs.cb_size_with_crc;
        }

        int n_filler = K - cb_k_info;
        int filler_start_comp = cb_k_info - n_punctured;
        int n_cb_comp = N_cb - n_filler;

        SoftVec recovered_comp = rate_matcher->rate_recover(cb_llr, n_cb_comp, rv, bgn, zc, qm, 0);

        SoftVec full_llr(N, 0.0);
        if (n_filler > 0 && filler_start_comp >= 0 && filler_start_comp <= n_cb_comp) {
            for (int i = 0; i < filler_start_comp; i++) {
                full_llr[n_punctured + i] = recovered_comp[i];
            }
            for (int i = 0; i < n_filler; i++) {
                full_llr[cb_k_info + i] = LLR_CLIP;
            }
            for (int i = filler_start_comp; i < n_cb_comp; i++) {
                int dst_idx = cb_k_info + n_filler + (i - filler_start_comp);
                if (dst_idx < N) {
                    full_llr[dst_idx] = recovered_comp[i];
                }
            }
        } else {
            for (int i = 0; i < n_cb_comp && i < N_cb; i++) {
                full_llr[n_punctured + i] = recovered_comp[i];
            }
        }

        for (int i = 0; i < N; i++) {
            if (full_llr[i] > LLR_CLIP) full_llr[i] = LLR_CLIP;
            else if (full_llr[i] < -LLR_CLIP) full_llr[i] = -LLR_CLIP;
        }

        if (c == 0) {
            llr_for_decoder_all.assign(full_llr.begin(), full_llr.end());
            write_double_npy(output_dir + "llr_for_decoder.npy", llr_for_decoder_all, {static_cast<size_t>(N)});
        }

        auto [decoded_cb_bits, conv_ok] = ldpc_decoder->decode(full_llr, bgn, zc,
                                                               n_ldpc_iter, true);
        (void)conv_ok;

        if ((int)decoded_cb_bits.size() < cb_k_info) {
            all_cb_ok = false;
            cb_decoded[c] = BitVec(cb_k_info, 0);
            continue;
        }

        BitVec cb_bits_out(decoded_cb_bits.begin(), decoded_cb_bits.begin() + cb_k_info);

        if (C > 1 && cb_crc_len > 0) {
            auto [cb_info_bits_out, cb_check_ok] = crc_encoder->decode(cb_bits_out, cb_crc_len);
            if (!cb_check_ok) {
                all_cb_ok = false;
            }
            cb_decoded[c] = cb_info_bits_out;
        } else {
            cb_decoded[c] = cb_bits_out;
        }
    }

    BitVec decoded_info_bits;
    if (C == 1) {
        int total_len = tbs + tb_crc_len;
        if ((int)cb_decoded[0].size() < total_len) {
            crc_ok = false;
        } else {
            BitVec rx_bits(cb_decoded[0].begin(), cb_decoded[0].begin() + total_len);
            auto [info_bits, tb_check_ok] = crc_encoder->decode(rx_bits, tb_crc_len);
            decoded_info_bits = info_bits;
            crc_ok = tb_check_ok;
        }
    } else {
        BitVec tb_bits_with_crc_decoded;
        for (int c = 0; c < C; c++) {
            tb_bits_with_crc_decoded.insert(tb_bits_with_crc_decoded.end(),
                                            cb_decoded[c].begin(), cb_decoded[c].end());
        }
        int total_tb_len = tbs + 24;
        if ((int)tb_bits_with_crc_decoded.size() < total_tb_len) {
            crc_ok = false;
        } else {
            BitVec rx_bits(tb_bits_with_crc_decoded.begin(),
                           tb_bits_with_crc_decoded.begin() + total_tb_len);
            auto [info_bits, tb_check_ok] = crc_encoder->decode(rx_bits, 24);
            decoded_info_bits = info_bits;
            crc_ok = tb_check_ok && all_cb_ok;
        }
    }

    write_uint8_npy(output_dir + "decoded_bits.npy", decoded_info_bits, {static_cast<size_t>(tbs)});

    {
        std::ofstream crc_file(output_dir + "crc_ok.txt");
        crc_file << (crc_ok ? "PASS" : "FAIL") << "\n";
        crc_file.close();
    }

    int bit_errors = 0;
    if ((int)decoded_info_bits.size() == tbs) {
        for (int i = 0; i < tbs; i++) {
            if (decoded_info_bits[i] != tb.bits[i]) {
                bit_errors++;
            }
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
