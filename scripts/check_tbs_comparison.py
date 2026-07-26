#!/usr/bin/env python3
"""Compare TBS calculation between C++ and Sionna for our target config"""
import numpy as np

MCS_TABLE2_PDSCH = [
    (2, 120), (2, 193), (2, 308), (2, 449), (2, 602),
    (4, 378), (4, 434), (4, 490), (4, 553), (4, 616), (4, 658),
    (6, 466), (6, 517), (6, 567), (6, 616), (6, 666), (6, 719), (6, 772), (6, 822), (6, 873),
    (8, 682.5), (8, 711), (8, 754), (8, 797), (8, 841), (8, 885), (8, 916.5), (8, 948),
]

def calculate_tbs_py(n_prb, n_re_per_prb, qm, n_layers, target_coderate):
    n_re_per_prb_tbs = min(156, n_re_per_prb)
    n_re = n_re_per_prb_tbs * n_prb
    target_tb_size = target_coderate * n_re * qm * n_layers
    
    if target_tb_size <= 3824:
        n_small = max(3, int(np.floor(np.log2(target_tb_size)) - 6))
        n_info_q = max(24, int(2**n_small * np.floor(target_tb_size / 2**n_small)))
    else:
        n_large = int(np.floor(np.log2(target_tb_size - 24)) - 5)
        n_info_q = max(3840, int(2**n_large * round((target_tb_size - 24) / 2**n_large)))
    
    if n_info_q <= 3824:
        num_cb = 1
    elif target_coderate <= 0.25:
        num_cb = int(np.ceil((n_info_q + 24) / 3816.0))
    elif n_info_q > 8424:
        num_cb = int(np.ceil((n_info_q + 24) / 8424.0))
    else:
        num_cb = 1
    
    crc_len = 24 if n_info_q > 3824 else 16
    
    tbs_table = [
        -1, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128,
        136, 144, 152, 160, 168, 176, 184, 192, 208, 224, 240, 256,
        272, 288, 304, 320, 336, 352, 368, 384, 408, 432, 456, 480,
        504, 528, 552, 576, 608, 640, 672, 704, 736, 768, 808, 848,
        888, 928, 984, 1032, 1064, 1128, 1160, 1192, 1224, 1256,
        1288, 1320, 1352, 1416, 1480, 1544, 1608, 1672, 1736, 1800,
        1864, 1928, 2024, 2088, 2152, 2216, 2280, 2408, 2472, 2536,
        2600, 2664, 2728, 2792, 2856, 2976, 3104, 3240, 3368, 3496,
        3624, 3752, 3824
    ]
    
    if n_info_q <= 3824:
        tbs = None
        for candidate in tbs_table:
            if candidate >= n_info_q:
                tbs = candidate
                break
        if tbs is None:
            tbs = 3824
    else:
        C = num_cb
        n_info_p = n_info_q
        if C > 1:
            n_info_p = n_info_q - 24
            n_cb_info = int(np.ceil(n_info_p / C / 8) * 8)
            tbs = C * n_cb_info - 24
        else:
            n_cb_info = int(np.ceil(n_info_p / 8) * 8)
            tbs = n_cb_info
    
    return tbs, crc_len, n_info_q

n_prb = 30
n_layers = 4
dmrs_symbols = [2, 7]
n_ofdm_sym = 14
n_dmrs_sym = len(dmrs_symbols)
n_pdsch_sym = n_ofdm_sym - n_dmrs_sym
n_re_per_prb = n_pdsch_sym * 12

print(f"Config: {n_prb}PRB, {n_layers}layers, DMRS symbols={dmrs_symbols}")
print(f"n_pdsch_sym={n_pdsch_sym}, n_re_per_prb={n_re_per_prb}")
print()

# Test representative MCS for each modulation
# QPSK: MCS 4 (Qm=2, R=602/1024)
# 16QAM: MCS 10 (Qm=4, R=658/1024)
# 64QAM: MCS 18 (Qm=6, R=772/1024)
# 256QAM: MCS 26 (Qm=8, R=916.5/1024)
test_mcs = [
    ("QPSK", 4),
    ("16QAM", 10),
    ("64QAM", 18),
    ("256QAM", 26),
]

for mod_name, mcs in test_mcs:
    Qm, R_x1024 = MCS_TABLE2_PDSCH[mcs]
    R = R_x1024 / 1024.0
    tbs, crc_len, n_info_q = calculate_tbs_py(n_prb, n_re_per_prb, Qm, n_layers, R)
    
    n_sc = n_prb * 12
    n_pdsch_re_total = 0
    for sym in range(n_ofdm_sym):
        if sym not in dmrs_symbols:
            n_pdsch_re_total += n_sc
    n_coded_total = n_pdsch_re_total * Qm * n_layers
    G = n_coded_total
    G = ((G // (Qm * n_layers)) * (Qm * n_layers))
    n_symbols_per_layer = G // (Qm * n_layers)
    
    print(f"{mod_name} MCS{mcs}: Qm={Qm}, R={R:.4f}")
    print(f"  n_info_q={n_info_q}, tbs(k_info)={tbs}, crc_len={crc_len}")
    print(f"  n_pdsch_re_total={n_pdsch_re_total}, G={G}, n_symbols_per_layer={n_symbols_per_layer}")
    print()
