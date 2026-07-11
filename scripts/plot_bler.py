#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import glob
import os

def main():
    parser = argparse.ArgumentParser(description='Plot BLER vs SINR curves')
    parser.add_argument('files', nargs='*', help='CSV result files')
    parser.add_argument('-o', '--output', default='bler_curve.png', help='Output image')
    parser.add_argument('-t', '--title', default='NR PDSCH BLER vs SINR', help='Plot title')
    args = parser.parse_args()
    
    if not args.files:
        args.files = glob.glob('*.csv')
    
    if not args.files:
        print("No CSV files found")
        return
    
    plt.figure(figsize=(10, 6))
    
    markers = ['o', 's', '^', 'D', 'v', 'p', '*']
    colors = ['b', 'r', 'g', 'm', 'c', 'y', 'k']
    
    for i, f in enumerate(args.files):
        try:
            df = pd.read_csv(f)
            label = os.path.splitext(os.path.basename(f))[0]
            plt.semilogy(df['SINR_dB'], df['BLER'], 
                        marker=markers[i % len(markers)],
                        color=colors[i % len(colors)],
                        linewidth=1.5, markersize=6,
                        label=label)
        except Exception as e:
            print("Error reading %s: %s" % (f, e))
    
    plt.xlabel('SINR (dB)', fontsize=12)
    plt.ylabel('BLER', fontsize=12)
    plt.title(args.title, fontsize=14)
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.legend(fontsize=10)
    plt.ylim([1e-3, 1])
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print("Plot saved to %s" % args.output)

if __name__ == '__main__':
    main()
