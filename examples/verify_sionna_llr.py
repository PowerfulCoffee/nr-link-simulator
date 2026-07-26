"""Verify Sionna LLR sign convention end-to-end."""
import sys
sys.path.insert(0, '/workspace/sionna/src')

import torch
import numpy as np
from sionna.phy.mapping import Mapper, Demapper, Constellation
from sionna.phy.fec.crc import CRCEncoder, CRCDecoder
from sionna.phy.fec.ldpc import LDPC5GEncoder, LDPC5GDecoder

torch.manual_seed(42)
np.random.seed(42)

print("="*60)
print("1. Verify LLR sign convention (QPSK, no noise)")
print("="*60)

mapper = Mapper("qam", 2)
demapper = Demapper("app", "qam", 2)

# All-zero bits (QPSK: I bits first, then Q bits)
bits = torch.zeros(1, 8, dtype=torch.int32)
x = mapper(bits)
print(f"TX bits: {bits[0]}")
print(f"TX symbols: {x[0]}")

# No noise: LLRs should strongly indicate bit=0
no = torch.tensor(0.01)
llr = demapper(x, no)
print(f"LLRs (no noise, all-zero bits): {llr[0]}")
print(f"LLR > 0 means bit=1 more likely (Sionna convention)")
print(f"Hard decision: {(llr < 0).int()[0]} (1 where LLR<0, i.e. bit=0)")
print()

# All-one bits
bits1 = torch.ones(1, 8, dtype=torch.int32)
x1 = mapper(bits1)
llr1 = demapper(x1, no)
print(f"TX bits: {bits1[0]}")
print(f"TX symbols: {x1[0]}")
print(f"LLRs (no noise, all-one bits): {llr1[0]}")
print()

print("="*60)
print("2. Verify CRC consistency")
print("="*60)

crc_enc = CRCEncoder("CRC24A")
crc_dec = CRCDecoder(crc_enc)

# Random info bits
k = 100
info_bits = torch.randint(0, 2, (1, k), dtype=torch.int32)
crc_bits = crc_enc(info_bits)
crc_checksum = crc_bits[0, k:]
print(f"Info bits (first 20): {info_bits[0,:20]}")
print(f"CRC24A checksum (24 bits): {crc_checksum.tolist()}")
print(f"Checksum length: {len(crc_checksum)}")

# Verify decode
crc_decode, crc_valid = crc_dec(crc_bits)
print(f"CRC decode matches info: {torch.equal(crc_decode, info_bits)}")
print(f"CRC valid: {crc_valid}")
print()

print("="*60)
print("3. Verify LDPC encode/decode (BGN=1, Z=4, G=1*Z+22=26, k=G*Z=104)")
print("="*60)

bgn = 1
z = 4
k = bgn * z + 22  # 4 + 22 = 26? Wait let me check
# Actually for 5G LDPC: k = 22*Z for BGN1 (BGN2 is 10*Z)
# Let me use a standard configuration
k = 256
n_ldpc = 512

ldpc_enc = LDPC5GEncoder(k, n_ldpc)
ldpc_dec = LDPC5GDecoder(ldpc_enc, num_iter=20, hard_out=True, return_infobits=True)
print(f"LDPC encoder: k={ldpc_enc.k}, n={ldpc_enc.n}, z={ldpc_enc.z}, bgn={ldpc_enc.bgn}")
print(f"Filler bits k_filler: {ldpc_enc.k_filler}")

# Random info bits
info = torch.randint(0, 2, (1, k), dtype=torch.int32)
cword = ldpc_enc(info)
print(f"Codeword length: {cword.shape[-1]}")

# Check codeword is binary
print(f"Codeword unique values: {torch.unique(cword).tolist()}")

# BPSK modulate: 0->+1, 1->-1 (Sionna convention for AWGN channel LLR test)
x = 1 - 2*cword
no_awgn = torch.tensor(0.5)
noise = torch.randn_like(x) * torch.sqrt(no_awgn/2)
y = x + noise

# Compute LLRs for BPSK/AWGN: llr = 2*y/no_awgn (for 0->+1, 1->-1: P(b=1)/P(b=0))
# Actually: P(y|b=0) ~ exp(-(y-1)^2/no), P(y|b=1) ~ exp(-(y+1)^2/no)
# LLR = ln(P(b=1)/P(b=0)) = [-(y+1)^2 + (y-1)^2]/no = [-4y]/no = -4y/no
# Wait let me compute correctly:
# squared_dist for constellation 0->+1: (y-1)^2, for 1->-1: (y+1)^2
# exponents: -dist/no
# LLR = (exp1 - exp0 via logsumexp) for b=1 vs b=0
# LLR = [-(y+1)^2/no] - [-(y-1)^2/no] = (-(y+1)^2 + (y-1)^2)/no = (-4y)/no = -4y/no
# Wait that gives: if y=+1 (b=0), LLR = -4/no < 0 → bit=0? But LLR>0=bit=1
# Hmm let me be careful: constellation points are at +1 and -1 for PAM/BPSK
# bit 0 maps to which point? Let me check:

pam_const = Constellation("pam", 1)
print(f"\nBPSK/PAM constellation points: {pam_const.points}")
print(f"Bit labels: {pam_const.bit_labels}")

# Use Sionna demapper directly for BPSK
llr_ch = demapper = None
# Let me use a simple BPSK demapper approach:
# For BPSK (1 bit per symbol), constellation has 2 points
# point 0: bit label [0] → value?
print(f"\nPAM bit_labels[0]: {pam_const.bit_labels[0].tolist()}, points[0]: {pam_const.points[0]}")
print(f"PAM bit_labels[1]: {pam_const.bit_labels[1].tolist()}, points[1]: {pam_const.points[1]}")

# Let's just use the Demapper for BPSK
demapper_bpsk = Demapper("app", "pam", 1)
y_bpsk = torch.complex(real=y, imag=torch.zeros_like(y))
llr_bpsk = demapper_bpsk(y_bpsk, no_awgn)
print(f"\nBPSK LLRs (first 10): {llr_bpsk[0,:10].tolist()}")
print(f"TX codeword (first 10): {cword[0,:10].tolist()}")
# Hard decision: LLR < 0 → bit=0? Let's check:
hard_llr = (llr_bpsk < 0).int()[0]
ber = (hard_llr != cword[0,:n_ldpc].int()).float().mean()
print(f"Hard decision BER (from demapper LLR): {ber.item():.4f}")

# Now decode with LDPC
info_hat = ldpc_dec(llr_bpsk)
ber_info = (info_hat != info).float().mean()
print(f"Decoded info BER: {ber_info.item():.4f}")
print(f"Info matches: {torch.equal(info_hat, info)}")
print()

print("="*60)
print("4. Check filler bits handling")
print("="*60)
# Filler bits are set to -llr_max at input of LDPC5GDecoder (before *-1)
# meaning after internal *-1 they become +llr_max → hard decision to 0
# Let's verify k_filler
print(f"k_filler = {ldpc_enc.k_filler}")
print(f"k = {ldpc_enc.k} (total k including filler bits)")
print(f"k - k_filler = {ldpc_enc.k - ldpc_enc.k_filler} (information bits)")
