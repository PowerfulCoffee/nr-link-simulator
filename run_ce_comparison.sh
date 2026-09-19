#!/bin/bash
cd /workspace/nr-link-simulator
RESULTS=results/ce_compare
mkdir -p $RESULTS

MAX=30
TERR=15
BATCH=4

run_cpp() {
    local label=$1
    local mcs=$2
    local snr_start=$3
    local snr_end=$4
    local ce_type=$5
    local outfile=$6
    echo "=== C++ $label MCS=$mcs ($ce_type) SNR=[$snr_start,$snr_end] ==="
    ./build/examples/pdsch_bler_tdl_mimo $mcs $snr_start $snr_end 1 $ce_type $MAX $TERR 30 4 4 1 15 5 1 0 2 $outfile 2>&1 | tail -20
}

run_sionna() {
    local label=$1
    local mcs=$2
    local snr_start=$3
    local snr_end=$4
    local ce_type=$5
    local outfile=$6
    echo "=== Sionna $label MCS=$mcs ($ce_type) SNR=[$snr_start,$snr_end] ==="
    python scripts/sionna_tdl_crosscheck.py --mcs $mcs --snr-start $snr_start --snr-end $snr_end --snr-step 1 \
        --max-blocks $MAX --target-errors $TERR --batch-size $BATCH --ce-type $ce_type \
        --delay-spread 15 --doppler 5 --output $outfile 2>&1 | tail -20
}

# QPSK MCS4 (Ideal CE: already have from earlier run, but re-run for consistency)
run_cpp "QPSK" 4 -2 6 1 $RESULTS/cpp_ideal_mcs4.csv
run_cpp "QPSK" 4 4 12 0 $RESULTS/cpp_ls_mcs4.csv

# 16QAM MCS8 (R=0.540)
run_cpp "16QAM" 8 6 16 1 $RESULTS/cpp_ideal_mcs8.csv
run_cpp "16QAM" 8 10 20 0 $RESULTS/cpp_ls_mcs8.csv

# 64QAM MCS15 (R=0.650)
run_cpp "64QAM" 15 14 24 1 $RESULTS/cpp_ideal_mcs15.csv
run_cpp "64QAM" 15 18 28 0 $RESULTS/cpp_ls_mcs15.csv

# Sionna QPSK MCS4
run_sionna "QPSK" 4 -2 6 ideal $RESULTS/sionna_ideal_mcs4.csv
run_sionna "QPSK" 4 4 12 ls $RESULTS/sionna_ls_mcs4.csv

# Sionna 16QAM MCS8
run_sionna "16QAM" 8 6 16 ideal $RESULTS/sionna_ideal_mcs8.csv
run_sionna "16QAM" 8 10 20 ls $RESULTS/sionna_ls_mcs8.csv

# Sionna 64QAM MCS15
run_sionna "64QAM" 15 14 24 ideal $RESULTS/sionna_ideal_mcs15.csv
run_sionna "64QAM" 15 18 28 ls $RESULTS/sionna_ls_mcs15.csv

echo "=== ALL SIMULATIONS COMPLETE ==="
ls -la $RESULTS/*.csv
