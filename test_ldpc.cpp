#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/ModuleFactory.h"
#include <iostream>
#include <iomanip>
#include <memory>

using namespace nr;
using namespace nr::phy;

int main() {
    auto crc_enc = create_crc_encoder();
    auto ldpc_enc = create_ldpc_encoder();
    auto ldpc_dec = create_ldpc_decoder();
    auto rate_match = create_rate_matcher();
    auto scrambler = create_scrambler();

    const int tb_size = 1000;
    BitVec tb_bits(tb_size, arma::fill::zeros);
    for (int i = 0; i < tb_size; i++) {
        tb_bits(i) = (i * 7 + 3) % 2;
    }

    int crc_len = get_crc_length(tb_size);
    BitVec bits_with_crc = crc_enc->encode(tb_bits, crc_len);
    std::cout << "TB size: " << tb_size << ", CRC len: " << crc_len << "\n";
    std::cout << "Bits with CRC: " << bits_with_crc.n_elem << " bits\n";

    LdpcBaseGraphInfo ldpc_info = select_ldpc_params(tb_size, 0.5);
    std::cout << "LDPC: bgn=" << ldpc_info.bgn << ", zc=" << ldpc_info.zc << "\n";

    BitVec coded_bits = ldpc_enc->encode(bits_with_crc, ldpc_info.bgn, ldpc_info.zc);
    std::cout << "LDPC encoded length: " << coded_bits.n_elem << "\n";

    int E = 3000;
    int rv = 0;
    BitVec rm_bits = rate_match->rate_match(coded_bits, E, rv, ldpc_info.bgn, ldpc_info.zc);
    std::cout << "Rate matched to E=" << E << ": " << rm_bits.n_elem << " bits\n";

    uint32_t cinit = 12345;
    BitVec scrambled = scrambler->scramble(rm_bits, cinit);

    BitVec scrambled_copy = scrambled;
    SoftVec llr_hard(E);
    for (int i = 0; i < E; i++) {
        llr_hard(i) = (scrambled_copy(i) == 0) ? 10.0 : -10.0;
    }

    SoftVec descrambled = scrambler->descramble(llr_hard, cinit);

    int N = static_cast<int>(coded_bits.n_elem);
    std::cout << "Rate recover N=" << N << "\n";
    SoftVec recovered = rate_match->rate_recover(descrambled, N, rv, ldpc_info.bgn, ldpc_info.zc);
    std::cout << "Recovered length: " << recovered.n_elem << "\n";

    auto [decoded_bits, conv_ok] = ldpc_dec->decode(recovered, ldpc_info.bgn, ldpc_info.zc, 50, true);
    std::cout << "Decoded length: " << decoded_bits.n_elem << ", conv=" << conv_ok << "\n";

    int total_len = tb_size + crc_len;
    bool crc_ok = false;
    BitVec info_bits;
    if (static_cast<int>(decoded_bits.n_elem) >= total_len) {
        BitVec rx_bits(total_len);
        for (int i = 0; i < total_len; i++) {
            rx_bits(i) = decoded_bits(i);
        }
        auto [ib, check] = crc_enc->decode(rx_bits, crc_len);
        info_bits = ib;
        crc_ok = check;
    }

    std::cout << "CRC ok: " << (crc_ok ? "YES" : "NO") << "\n";

    if (crc_ok) {
        int errors = 0;
        int min_len = std::min(tb_size, (int)info_bits.n_elem);
        for (int i = 0; i < min_len; i++) {
            if (tb_bits(i) != info_bits(i)) {
                errors++;
                if (errors <= 10) {
                    std::cout << "Bit error at " << i << ": tx=" << (int)tb_bits(i)
                              << ", rx=" << (int)info_bits(i) << "\n";
                }
            }
        }
        std::cout << "Total bit errors: " << errors << "\n";
    }

    return crc_ok ? 0 : 1;
}
