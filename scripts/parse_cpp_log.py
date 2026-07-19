import re
import csv
import os

log_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'bler_all_mcs_cpp.log')
out_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'cpp_all_mcs_parsed.csv')

MCS_TABLE = [
    (2, 120), (2, 193), (2, 308), (2, 449), (2, 602),
    (4, 378), (4, 434), (4, 490), (4, 553), (4, 616), (4, 658),
    (6, 466), (6, 517), (6, 567), (6, 616), (6, 666),
    (6, 719), (6, 772), (6, 822), (6, 873),
    (8, 682.5), (8, 711), (8, 754), (8, 797), (8, 841), (8, 885), (8, 916.5), (8, 948),
]
CPP_TBS = {
    0:104, 1:176, 2:288, 3:408, 4:552,
    5:704, 6:808, 7:888, 8:1032, 9:1128, 10:1224,
    11:1288, 12:1416, 13:1608, 14:1736, 15:1864,
    16:2024, 17:2152, 18:2280, 19:2408,
    20:2472, 21:2600, 22:2792, 23:2976, 24:3104, 25:3240, 26:3368, 27:3496
}

mcs_re = re.compile(r'=== MCS\s+(\d+)\s+\(Qm=(\d+),\s+R=([\d.]+)')
snr_re = re.compile(r'SNR\s+([-\d.]+)\s+dB:\s+BLER=([\d.]+)\s+\(blocks=(\d+),\s+errors=(\d+)\)')

current_mcs = None
results = []

with open(log_path, 'r') as f:
    for line in f:
        m = mcs_re.search(line)
        if m:
            current_mcs = int(m.group(1))
            continue
        s = snr_re.search(line)
        if s and current_mcs is not None:
            snr = float(s.group(1))
            bler = float(s.group(2))
            blocks = int(s.group(3))
            errors = int(s.group(4))
            qm, r1024 = MCS_TABLE[current_mcs]
            R = r1024 / 1024.0
            tbs = CPP_TBS[current_mcs]
            results.append((current_mcs, qm, R, tbs, snr, bler, blocks, errors))

with open(out_csv, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['mcs','qm','rate','tbs','snr_db','bler','n_blocks','n_errors'])
    for r in results:
        w.writerow([r[0], r[1], f"{r[2]:.4f}", r[3], f"{r[4]:.1f}", f"{r[5]:.6f}", r[6], r[7]])

print(f"Parsed {len(results)} data points from {len(set(r[0] for r in results))} MCS indices")
print(f"Saved to {out_csv}")

bler10_cpp = {}
for mcs in sorted(set(r[0] for r in results)):
    mcs_data = [(r[4], r[5]) for r in results if r[0] == mcs]
    mcs_data.sort()
    for i in range(len(mcs_data)-1):
        s1, b1 = mcs_data[i]
        s2, b2 = mcs_data[i+1]
        if (b1 - 0.1) * (b2 - 0.1) <= 0 and b1 != b2:
            s_bler10 = s1 + (0.1 - b1) * (s2 - s1) / (b2 - b1)
            bler10_cpp[mcs] = s_bler10
            break
    if mcs not in bler10_cpp:
        for s, b in mcs_data:
            if b < 0.1:
                bler10_cpp[mcs] = s
                break

print("\nBLER=0.1 SNR thresholds (C++):")
mod_names = {2:'QPSK', 4:'16QAM', 6:'64QAM', 8:'256QAM'}
for mcs in sorted(bler10_cpp.keys()):
    qm, r1024 = MCS_TABLE[mcs]
    print(f"  MCS{mcs:2d} ({mod_names[qm]:6s}, R={r1024/1024:.3f}): SNR = {bler10_cpp[mcs]:5.2f} dB")
