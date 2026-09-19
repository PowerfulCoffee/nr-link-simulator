#!/bin/bash
set -e
cd /workspace/nr-link-simulator
RESULTS=results/tdl_crosscheck_4x4
mkdir -p $RESULTS

BLOCKS=30
ERRORS=15
BATCH=3

echo "=== Sionna Ideal CE (1dB step) ==="

echo "QPSK MCS4 Ideal CE (-2 to 6 dB)..."
python scripts/sionna_tdl_crosscheck.py --mcs 4 --snr-start -2 --snr-end 6 --snr-step 1 \
    --max-blocks $BLOCKS --target-errors $ERRORS --batch-size $BATCH \
    --channel tdl --delay-spread 100 --max-doppler 5 \
    --dmrs-add-pos 1 --dmrs-duration 0 --ce-type ideal \
    --ldpc-iters 20 --mcs-table 2 \
    --output $RESULTS/sionna_ideal_mcs4_1db.csv

echo "16QAM MCS10 Ideal CE (6 to 16 dB)..."
python scripts/sionna_tdl_crosscheck.py --mcs 10 --snr-start 6 --snr-end 16 --snr-step 1 \
    --max-blocks $BLOCKS --target-errors $ERRORS --batch-size $BATCH \
    --channel tdl --delay-spread 100 --max-doppler 5 \
    --dmrs-add-pos 1 --dmrs-duration 0 --ce-type ideal \
    --ldpc-iters 20 --mcs-table 2 \
    --output $RESULTS/sionna_ideal_mcs10_1db.csv

echo "64QAM MCS16 Ideal CE (14 to 26 dB)..."
python scripts/sionna_tdl_crosscheck.py --mcs 16 --snr-start 14 --snr-end 26 --snr-step 1 \
    --max-blocks $BLOCKS --target-errors $ERRORS --batch-size 2 \
    --channel tdl --delay-spread 100 --max-doppler 5 \
    --dmrs-add-pos 1 --dmrs-duration 0 --ce-type ideal \
    --ldpc-iters 20 --mcs-table 2 \
    --output $RESULTS/sionna_ideal_mcs16_1db.csv

echo "256QAM MCS24 Ideal CE (22 to 32 dB)..."
python scripts/sionna_tdl_crosscheck.py --mcs 24 --snr-start 22 --snr-end 32 --snr-step 1 \
    --max-blocks $BLOCKS --target-errors $ERRORS --batch-size 2 \
    --channel tdl --delay-spread 100 --max-doppler 5 \
    --dmrs-add-pos 1 --dmrs-duration 0 --ce-type ideal \
    --ldpc-iters 20 --mcs-table 2 \
    --output $RESULTS/sionna_ideal_mcs24_1db.csv

echo ""
echo "=== Sionna LS CE (1dB step) ==="

echo "QPSK MCS4 LS CE (0 to 10 dB)..."
python scripts/sionna_tdl_crosscheck.py --mcs 4 --snr-start 0 --snr-end 10 --snr-step 1 \
    --max-blocks $BLOCKS --target-errors $ERRORS --batch-size $BATCH \
    --channel tdl --delay-spread 100 --max-doppler 5 \
    --dmrs-add-pos 1 --dmrs-duration 0 --ce-type ls \
    --ldpc-iters 20 --mcs-table 2 \
    --output $RESULTS/sionna_ls_mcs4.csv

echo "16QAM MCS10 LS CE (10 to 22 dB)..."
python scripts/sionna_tdl_crosscheck.py --mcs 10 --snr-start 10 --snr-end 22 --snr-step 1 \
    --max-blocks $BLOCKS --target-errors $ERRORS --batch-size $BATCH \
    --channel tdl --delay-spread 100 --max-doppler 5 \
    --dmrs-add-pos 1 --dmrs-duration 0 --ce-type ls \
    --ldpc-iters 20 --mcs-table 2 \
    --output $RESULTS/sionna_ls_mcs10.csv

echo "64QAM MCS16 LS CE (18 to 32 dB)..."
python scripts/sionna_tdl_crosscheck.py --mcs 16 --snr-start 18 --snr-end 32 --snr-step 1 \
    --max-blocks $BLOCKS --target-errors $ERRORS --batch-size 2 \
    --channel tdl --delay-spread 100 --max-doppler 5 \
    --dmrs-add-pos 1 --dmrs-duration 0 --ce-type ls \
    --ldpc-iters 20 --mcs-table 2 \
    --output $RESULTS/sionna_ls_mcs16.csv

echo "256QAM MCS24 LS CE (26 to 40 dB)..."
python scripts/sionna_tdl_crosscheck.py --mcs 24 --snr-start 26 --snr-end 40 --snr-step 1 \
    --max-blocks $BLOCKS --target-errors $ERRORS --batch-size 2 \
    --channel tdl --delay-spread 100 --max-doppler 5 \
    --dmrs-add-pos 1 --dmrs-duration 0 --ce-type ls \
    --ldpc-iters 20 --mcs-table 2 \
    --output $RESULTS/sionna_ls_mcs24.csv

echo ""
echo "All Sionna simulations done!"
