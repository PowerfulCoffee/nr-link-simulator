#!/bin/bash
set -e
cd /workspace/nr-link-simulator
RESULTS=results/tdl_crosscheck_4x4
mkdir -p $RESULTS

echo "Running Sionna QPSK MCS4..."
python scripts/sionna_tdl_crosscheck.py \
    --mcs 4 --snr-start -4 --snr-end 12 --snr-step 2 \
    --max-blocks 50 --target-errors 20 --batch-size 10 \
    --channel tdl --delay-spread 100 --max-doppler 5 \
    --dmrs-add-pos 1 --dmrs-duration 0 --ce-type ideal \
    --ldpc-iters 20 --mcs-table 2 \
    --output $RESULTS/sionna_ideal_mcs4.csv \
    > $RESULTS/sionna_ideal_mcs4.log 2>&1

echo "Running Sionna 16QAM MCS10..."
python scripts/sionna_tdl_crosscheck.py \
    --mcs 10 --snr-start 4 --snr-end 24 --snr-step 2 \
    --max-blocks 50 --target-errors 20 --batch-size 10 \
    --channel tdl --delay-spread 100 --max-doppler 5 \
    --dmrs-add-pos 1 --dmrs-duration 0 --ce-type ideal \
    --ldpc-iters 20 --mcs-table 2 \
    --output $RESULTS/sionna_ideal_mcs10.csv \
    > $RESULTS/sionna_ideal_mcs10.log 2>&1

echo "Running Sionna 64QAM MCS16..."
python scripts/sionna_tdl_crosscheck.py \
    --mcs 16 --snr-start 10 --snr-end 30 --snr-step 2 \
    --max-blocks 50 --target-errors 20 --batch-size 5 \
    --channel tdl --delay-spread 100 --max-doppler 5 \
    --dmrs-add-pos 1 --dmrs-duration 0 --ce-type ideal \
    --ldpc-iters 20 --mcs-table 2 \
    --output $RESULTS/sionna_ideal_mcs16.csv \
    > $RESULTS/sionna_ideal_mcs16.log 2>&1

echo "Running Sionna 256QAM MCS24..."
python scripts/sionna_tdl_crosscheck.py \
    --mcs 24 --snr-start 16 --snr-end 36 --snr-step 2 \
    --max-blocks 50 --target-errors 20 --batch-size 5 \
    --channel tdl --delay-spread 100 --max-doppler 5 \
    --dmrs-add-pos 1 --dmrs-duration 0 --ce-type ideal \
    --ldpc-iters 20 --mcs-table 2 \
    --output $RESULTS/sionna_ideal_mcs24.csv \
    > $RESULTS/sionna_ideal_mcs24.log 2>&1

echo "All Sionna simulations done!"
