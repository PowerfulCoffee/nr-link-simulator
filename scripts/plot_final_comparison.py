import os, csv, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

results_dir = '/workspace/nr-link-simulator/results'
out_dir = '/workspace/nr-link-simulator/output'

def load_csv(path):
    data = {}
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            mcs = int(row['mcs'])
            snr = float(row['snr_db'])
            bler = float(row['bler'])
            if mcs not in data: data[mcs] = []
            data[mcs].append((snr, bler))
    for mcs in data: data[mcs].sort()
    return data

def bler10(sbl):
    sbl = sorted(sbl)
    for i in range(len(sbl)-1):
        s1, b1 = sbl[i]; s2, b2 = sbl[i+1]
        if b1 >= 0.1 and b2 <= 0.1:
            return s1 + (0.1-b1)*(s2-s1)/(b2-b1)
    return None

cpp = load_csv(os.path.join(results_dir, 'cpp_all_mcs_parsed.csv'))
sionna_fast = load_csv(os.path.join(results_dir, 'sionna_all_mcs_awgn_fast.csv'))
refine_path = os.path.join(results_dir, 'sionna_refined_key_mcs.csv')
sionna_refine = load_csv(refine_path) if os.path.exists(refine_path) else {}

sionna = dict(sionna_fast)
for mcs, pts in sionna_refine.items():
    rs = set(round(s,1) for s,b in pts)
    existing = [(s,b) for s,b in sionna.get(mcs,[]) if round(s,1) not in rs]
    sionna[mcs] = sorted(existing + pts)

MCS = [(2,120),(2,193),(2,308),(2,449),(2,602),(4,378),(4,434),(4,490),(4,553),(4,616),(4,658),
       (6,466),(6,517),(6,567),(6,616),(6,666),(6,719),(6,772),(6,822),(6,873),
       (8,682.5),(8,711),(8,754),(8,797),(8,841),(8,885),(8,916.5),(8,948)]
modn = {2:'QPSK',4:'16QAM',6:'64QAM',8:'256QAM'}
modc = {2:'#1f77b4',4:'#ff7f0e',6:'#2ca02c',8:'#d62728'}

fig, axes = plt.subplots(2, 2, figsize=(18, 14))
fig.suptitle('NR PDSCH BLER vs Es/N0: C++ (Offset Min-Sum) vs Sionna (BP) — AWGN SISO, 3PRB',
             fontsize=14, fontweight='bold', y=0.98)

groups = {'QPSK (MCS 0-4)':range(0,5),'16QAM (MCS 5-10)':range(5,11),
          '64QAM (MCS 11-19)':range(11,20),'256QAM (MCS 20-27)':range(20,28)}

for ai,(title,mcs_list) in enumerate(groups.items()):
    ax = axes[ai//2][ai%2]
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel('Es/N0 (dB)', fontsize=10)
    ax.set_ylabel('BLER', fontsize=10)
    ax.set_yscale('log'); ax.set_ylim(5e-4, 1.0); ax.grid(True, alpha=0.3, which='both')
    ax.axhline(y=0.1, color='gray', ls='--', alpha=0.5, lw=0.8)
    for mcs in mcs_list:
        qm,r = MCS[mcs]; R = r/1024.0
        c = modc[qm]; a = 0.35 + 0.65*(mcs-mcs_list[0])/max(1,len(mcs_list)-1)
        if mcs in cpp:
            s = [x[0] for x in cpp[mcs]]; b = [x[1] for x in cpp[mcs]]
            ax.plot(s,b,'-',c=c,alpha=a,lw=1.6,label=f'MCS{mcs} C++ (R={R:.2f})')
            ax.plot(s,b,'o',c=c,alpha=a,ms=3)
        if mcs in sionna:
            s = [x[0] for x in sionna[mcs]]; b = [x[1] for x in sionna[mcs]]
            ax.plot(s,b,'--',c=c,alpha=a,lw=1.6,label=f'MCS{mcs} Sionna (R={R:.2f})')
            ax.plot(s,b,'s',c=c,alpha=a,ms=3)
    ax.legend(fontsize=6.5, loc='lower left', ncol=2)
plt.tight_layout(rect=[0,0,1,0.96])
plt.savefig(os.path.join(out_dir,'bler_all_mcs_waterfall.png'), dpi=150, bbox_inches='tight')
plt.close()
print("Saved waterfall plot")

fig2, (ax1,ax2) = plt.subplots(1,2,figsize=(17,6))
b10c, b10s = {}, {}
for mcs in range(28):
    v = bler10(cpp.get(mcs,[]))
    if v is not None: b10c[mcs]=v
    v = bler10(sionna.get(mcs,[]))
    if v is not None: b10s[mcs]=v

mc = sorted(b10c.keys())
for qm in [2,4,6,8]:
    idx = [i for i,m in enumerate(mc) if MCS[m][0]==qm]
    if idx:
        ax1.scatter([MCS[mc[i]][1]/1024 for i in idx],[b10c[mc[i]] for i in idx],
                    c=modc[qm],label=f'{modn[qm]} C++',marker='o',s=60,zorder=5,edgecolors='k',linewidth=0.5)
ms = sorted(b10s.keys())
for qm in [2,4,6,8]:
    idx = [i for i,m in enumerate(ms) if MCS[m][0]==qm]
    if idx:
        ax1.scatter([MCS[ms[i]][1]/1024 for i in idx],[b10s[ms[i]] for i in idx],
                    c=modc[qm],label=f'{modn[qm]} Sionna',marker='s',s=50,zorder=4,alpha=0.6)

ax1.set_xlabel('Code Rate R', fontsize=11)
ax1.set_ylabel('Es/N0 @ BLER=0.1 (dB)', fontsize=11)
ax1.set_title('Required SNR for BLER=0.1 vs Code Rate', fontsize=12, fontweight='bold')
ax1.grid(True, alpha=0.3); ax1.legend(fontsize=8, ncol=2, loc='upper left')

common = sorted(set(b10c.keys())&set(b10s.keys()))
diffs = [b10c[m]-b10s[m] for m in common]
for m in common:
    qm = MCS[m][0]
    ax2.bar(m, b10c[m]-b10s[m], color=modc[qm], alpha=0.8, width=0.7, edgecolor='k', linewidth=0.3)
if diffs:
    md = np.mean(diffs)
    ax2.axhline(md, color='black', ls='--', lw=1.5, label=f'Mean = {md:+.2f} dB')
ax2.set_title('C++ − Sionna SNR Gap @ BLER=0.1 (positive: C++ worse)', fontsize=12, fontweight='bold')
ax2.set_xlabel('MCS Index', fontsize=11); ax2.set_ylabel('SNR Gap (dB)', fontsize=11)
ax2.grid(True, alpha=0.3, axis='y'); ax2.legend(fontsize=10)
ax2.set_xticks(range(0,28,2))
ax2.axhspan(-0.3, 0.3, alpha=0.1, color='green')
plt.tight_layout()
plt.savefig(os.path.join(out_dir,'bler_all_mcs_summary.png'), dpi=150, bbox_inches='tight')
plt.close()
print("Saved summary plot")

print("\n=== C++ vs Sionna BLER=0.1 Gap Summary ===")
print(f"MCS compared: {len(common)} (MCS0-1 skipped: rate<0.2 needs repetition)")
print(f"Mean gap: {np.mean(diffs):+.2f} dB")
print(f"Std dev:  {np.std(diffs):.2f} dB")
print(f"Within ±0.5dB: {sum(1 for d in diffs if abs(d)<=0.5)}/{len(diffs)} MCS")
