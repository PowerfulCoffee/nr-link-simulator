"""Quick Sionna AWGN BLER sweep to find SNR range"""
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import torch
import numpy as np
from sionna.phy.fec.ldpc import LDPC5GEncoder, LDPC5GDecoder
from sionna.phy.fec.crc import CRCEncoder, CRCDecoder
from sionna.phy.mapping import Constellation, Mapper, Demapper

crc_enc = CRCEncoder("CRC24A")
crc_dec = CRCDecoder(crc_enc)

MCS_LIST = {
    3:  {"Qm": 2, "k_info": 1800, "n_coded": 7200},
    10: {"Qm": 4, "k_info": 4736, "n_coded": 14400},
    17: {"Qm": 6, "k_info": 9224, "n_coded": 21600},
    24: {"Qm": 6, "k_info": 16392, "n_coded": 21600},
}

n_prb = 25
n_sc = n_prb * 12
n_ofdm_sym = 14
dmrs_pos = [2, 11]
dmrs_scs = list(range(0, n_sc, 2))
n_pdsch_sym = n_ofdm_sym - len(dmrs_pos)
n_pdsch_re = n_pdsch_sym * n_sc

for mcs, p in MCS_LIST.items():
    Qm = p["Qm"]
    k_info = p["k_info"]
    n_coded = p["n_coded"]
    k_crc = k_info + 24
    n_symbols = n_coded // Qm
    
    print(f"\n=== MCS{mcs}: Qm={Qm}, k_info={k_info}, n_coded={n_coded}, rate={k_crc/n_coded:.3f} ===")
    
    ldpc_enc = LDPC5GEncoder(k_crc, n_coded, num_bits_per_symbol=Qm)
    ldpc_dec = LDPC5GDecoder(ldpc_enc, num_iter=20, hard_out=True, return_info=False)
    constellation = Constellation("qam", Qm)
    mapper = Mapper(constellation=constellation)
    demapper = Demapper("app", constellation=constellation)
    
    if mcs == 3:
        snrs = range(-6, 12, 2)
    elif mcs == 10:
        snrs = range(0, 22, 2)
    elif mcs == 17:
        snrs = range(6, 28, 2)
    else:
        snrs = range(12, 34, 2)
    
    for esn0_db in snrs:
        no = 10.0 ** (-esn0_db / 10.0)
        n_errors = 0
        n_blocks = 0
        bs = 50
        n_batches = 4
        for _ in range(n_batches):
            b = torch.randint(0, 2, (bs, k_info), dtype=torch.float32)
            b_crc = crc_enc(b)
            c = ldpc_enc(b_crc)
            x = mapper(c)[:, :n_symbols]
            
            rx = torch.zeros(bs, n_ofdm_sym, n_sc, dtype=torch.complex64)
            tx = torch.zeros(bs, n_ofdm_sym, n_sc, dtype=torch.complex64)
            
            idx = 0
            for sym in range(n_ofdm_sym):
                if sym in dmrs_pos:
                    for sc in dmrs_scs:
                        tx[:, sym, sc] = 1.0
                else:
                    for sc in range(n_sc):
                        if idx < n_symbols:
                            tx[:, sym, sc] = x[:, idx]
                            idx += 1
            
            h = torch.ones(bs, n_ofdm_sym, n_sc, dtype=torch.complex64)
            noise_std = np.sqrt(no / 2.0)
            noise = (torch.randn_like(tx, dtype=torch.float32) + 1j*torch.randn_like(tx, dtype=torch.float32)) * noise_std
            rx = tx * h + noise
            
            llr_list = []
            idx = 0
            for sym in range(n_ofdm_sym):
                if sym in dmrs_pos:
                    continue
                for sc in range(n_sc):
                    if idx < n_symbols:
                        y_eq = rx[:, sym, sc]
                        llr = demapper(y_eq, torch.tensor(no, dtype=torch.float32))
                        llr_list.append(llr)
                        idx += 1
            
            llr_cat = torch.stack(llr_list, dim=1).reshape(bs, -1)[:, :n_coded]
            c_hat = ldpc_dec(llr_cat)
            _, crc_valid = crc_dec(c_hat)
            n_errors += int((crc_valid < 0.5).float().sum())
            n_blocks += bs
        
        bler = n_errors / n_blocks
        print(f"  EsN0={esn0_db:2d}dB: BLER={bler:.3f}")
        if bler < 0.01:
            break
