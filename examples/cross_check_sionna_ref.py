"""Cross-check script for CRC, LLR demodulation, and LDPC against Sionna."""
import sys
sys.path.insert(0, '/workspace/sionna/src')

import torch
import numpy as np
from sionna.phy.mapping import Mapper, Demapper, Constellation, pam_gray
from sionna.phy.fec.crc import CRCEncoder, CRCDecoder
from sionna.phy.fec.ldpc import LDPC5GEncoder, LDPC5GDecoder

torch.manual_seed(123)
np.random.seed(123)

print("="*70)
print("CROSS-CHECK: QPSK/16QAM/64QAM Constellation")
print("="*70)

for order, name in [(2,"QPSK"), (4,"16QAM"), (6,"64QAM")]:
    const = Constellation("qam", order)
    points = const().detach().numpy()
    print(f"\n{name} constellation (normalized power=1):")
    m = order // 2
    for i, p in enumerate(points):
        bits_str = ''.join(str((i >> (order-1-j)) & 1) for j in range(order))
        print(f"  bits=[{','.join(bits_str)}] -> ({p.real:+.6f},{p.imag:+.6f})")

print()
print("="*70)
print("CROSS-CHECK: PAM Gray mapping (1D)")
print("="*70)
for m in [1, 2, 3]:
    M = 2**m
    print(f"\nPAM m={m}:")
    for k in range(M):
        bits = np.array([(k >> (m-1-i)) & 1 for i in range(m)])
        val = pam_gray(bits)
        bits_str = ','.join(str(b) for b in bits)
        print(f"  bits=[{bits_str}] -> {val:.0f}")

print()
print("="*70)
print("CROSS-CHECK: CRC24A")
print("="*70)

crc_enc = CRCEncoder("CRC24A")
crc_dec = CRCDecoder(crc_enc)

info_bits = torch.tensor([[0,1,0,0,0,0,0,1,1,0,1,1,0,1,0,1,0,1,1,0]], dtype=torch.int32)
crc_out = crc_enc(info_bits)
crc_bits = crc_out[0, info_bits.shape[1]:].tolist()
print(f"Info bits: {info_bits[0].tolist()}")
print(f"CRC24A: {''.join(str(b) for b in crc_bits)} = {crc_bits}")
crc_decoded, crc_valid = crc_dec(crc_out)
print(f"CRC verify: {crc_valid.item()}")

print()
print("="*70)
print("CROSS-CHECK: LDPC5G Encode/Decode (perfect LLR)")
print("="*70)

k_info = 200
crc_len = 24
k_total = k_info + crc_len
n_ldpc = 480

crc_enc2 = CRCEncoder("CRC24A")
ldpc_enc2 = LDPC5GEncoder(k_total, n_ldpc)
ldpc_dec2 = LDPC5GDecoder(ldpc_enc2, num_iter=20, hard_out=True, return_infobits=True, cn_update="boxplus-phi")
crc_dec2 = CRCDecoder(crc_enc2)

print(f"k_info={k_info}, k_total={k_total}, n_ldpc={n_ldpc}")
info2 = torch.randint(0, 2, (1, k_info), dtype=torch.int32)
info_crc = crc_enc2(info2)
cword2 = ldpc_enc2(info_crc.to(torch.float32))
print(f"Codeword length n={cword2.shape[-1]}")

# Perfect LLR: ln(P(b=1)/P(b=0)), bit=0 -> -inf, bit=1 -> +inf
# Sionna LDPC5GDecoder clips to [-20, 20]
llr_ch2 = (2.0 * cword2 - 1.0) * 20.0
decoded_llr = ldpc_dec2(llr_ch2)
decoded_info, crc_ok = crc_dec2(decoded_llr.to(torch.int32))
print(f"Perfect decode: info matches={torch.equal(decoded_info.int(), info2.int())}, crc_valid={crc_ok.item()}")

# Save reference data
np.save('/workspace/nr-link-simulator/examples/ref_info.npy', info2[0].numpy().astype(np.int32))
np.save('/workspace/nr-link-simulator/examples/ref_info_crc.npy', info_crc[0].numpy().astype(np.int32))
np.save('/workspace/nr-link-simulator/examples/ref_codeword.npy', cword2[0].numpy().astype(np.int32))
print(f"\nReference saved: k_info={k_info}, k_total={k_total}, n={cword2.shape[-1]}")
print(f"Info bits (first 40): {info2[0,:40].tolist()}")
print(f"Codeword first 40 bits: {cword2[0,:40].to(torch.int32).tolist()}")
