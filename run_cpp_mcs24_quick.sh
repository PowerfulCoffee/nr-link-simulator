#!/bin/bash
set -e
cd /workspace/nr-link-simulator
RESULTS=results/tdl_crosscheck_4x4

echo "Running C++ 256QAM MCS24 key SNR points..."

echo "SNR_dB,Blocks,Errors,BLER" > $RESULTS/cpp_ideal_mcs24.csv

for snr in 20 22 24 26 28 30; do
    echo "  SNR $snr dB..."
    ./build/examples/pdsch_bler_tdl_mimo 24 $snr $snr 2 1 30 15 30 4 4 1 100 5 1 0 2 /tmp/tmp_mcs24.csv > /dev/null 2>&1
    if [ -f /tmp/tmp_mcs24.csv ]; then
        tail -n 1 /tmp/tmp_mcs24.csv >> $RESULTS/cpp_ideal_mcs24.csv
    fi
done

echo "C++ 256QAM done!"
cat $RESULTS/cpp_ideal_mcs24.csv
