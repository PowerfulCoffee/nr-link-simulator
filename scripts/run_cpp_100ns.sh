#!/bin/bash
set -e

BUILD_DIR="/workspace/nr-link-simulator/build/examples"
RESULTS_DIR="/workspace/nr-link-simulator/results/ce_compare_100ns"
CPP_BIN="${BUILD_DIR}/pdsch_bler_tdl_mimo"

# Common parameters
N_RB=30
N_RX=4
N_LAYERS=4
CHANNEL_TYPE=1  # TDL-A
DELAY_SPREAD=100  # ns
MAX_DOPPLER=10  # Hz
DMRS_ADD_POS=1
DMRS_DURATION=0
MCS_TABLE=2
MAX_BLOCKS=100
TARGET_ERRORS=30
SNR_STEP=1

echo "========================================"
echo "Running C++ BLER simulations (100ns, 10Hz)"
echo "========================================"

run_cpp_sim() {
    local mcs=$1
    local ce_type=$2
    local snr_start=$3
    local snr_end=$4
    local perfect_csi=0
    local label="ls"
    
    if [ "$ce_type" = "ideal" ]; then
        perfect_csi=1
        label="ideal"
    fi
    
    local outfile="${RESULTS_DIR}/cpp_${label}_mcs${mcs}.csv"
    
    if [ -f "$outfile" ]; then
        echo "Skipping MCS${mcs} ${ce_type} (output exists)"
        return
    fi
    
    echo ""
    echo "=== Running C++ MCS${mcs} ${ce_type} CE: SNR ${snr_start}~${snr_end}dB ==="
    
    $CPP_BIN $mcs $snr_start $snr_end $SNR_STEP $perfect_csi $MAX_BLOCKS $TARGET_ERRORS \
              $N_RB $N_RX $N_LAYERS $CHANNEL_TYPE $DELAY_SPREAD $MAX_DOPPLER \
              $DMRS_ADD_POS $DMRS_DURATION $MCS_TABLE $outfile
}

# MCS4: QPSK
run_cpp_sim 4 ideal -2 10
run_cpp_sim 4 ls 4 20

# MCS8: 16QAM
run_cpp_sim 8 ideal 6 20
run_cpp_sim 8 ls 12 26

# MCS15: 64QAM
run_cpp_sim 15 ideal 14 28
run_cpp_sim 15 ls 20 35

# MCS27: 256QAM
run_cpp_sim 27 ideal 22 35
run_cpp_sim 27 ls 28 42

echo ""
echo "========================================"
echo "C++ simulations complete!"
echo "========================================"
ls -la ${RESULTS_DIR}/cpp_*.csv
