#!/bin/bash
set -e
cd /workspace/nr-link-simulator
mkdir -p results
RESULTS=results

MAX_BLOCKS=50
TARGET_ERRORS=20
BATCH=2
DELAY_SPREAD=15
DOPPLER=5
N_PRB=30
N_RX=4
N_LAYERS=4
CH_TYPE=1
DMRS_ADD_POS=1
DMRS_DUR=0
MCS_TABLE=2

echo "=== Running C++ Ideal CE ==="
# QPSK MCS6
./build/examples/pdsch_bler_tdl_mimo 6 0 10 1 1 $MAX_BLOCKS $TARGET_ERRORS $N_PRB $N_RX $N_LAYERS $CH_TYPE $DELAY_SPREAD $DOPPLER $DMRS_ADD_POS $DMRS_DUR $MCS_TABLE $RESULTS/cpp_ideal_mcs6.csv
# 16QAM MCS12
./build/examples/pdsch_bler_tdl_mimo 12 8 18 1 1 $MAX_BLOCKS $TARGET_ERRORS $N_PRB $N_RX $N_LAYERS $CH_TYPE $DELAY_SPREAD $DOPPLER $DMRS_ADD_POS $DMRS_DUR $MCS_TABLE $RESULTS/cpp_ideal_mcs12.csv
# 64QAM MCS23
./build/examples/pdsch_bler_tdl_mimo 23 15 27 1 1 $MAX_BLOCKS $TARGET_ERRORS $N_PRB $N_RX $N_LAYERS $CH_TYPE $DELAY_SPREAD $DOPPLER $DMRS_ADD_POS $DMRS_DUR $MCS_TABLE $RESULTS/cpp_ideal_mcs23.csv

echo "=== Running C++ LS CE ==="
./build/examples/pdsch_bler_tdl_mimo 6 3 13 1 0 $MAX_BLOCKS $TARGET_ERRORS $N_PRB $N_RX $N_LAYERS $CH_TYPE $DELAY_SPREAD $DOPPLER $DMRS_ADD_POS $DMRS_DUR $MCS_TABLE $RESULTS/cpp_ls_mcs6.csv
./build/examples/pdsch_bler_tdl_mimo 12 11 21 1 0 $MAX_BLOCKS $TARGET_ERRORS $N_PRB $N_RX $N_LAYERS $CH_TYPE $DELAY_SPREAD $DOPPLER $DMRS_ADD_POS $DMRS_DUR $MCS_TABLE $RESULTS/cpp_ls_mcs12.csv
./build/examples/pdsch_bler_tdl_mimo 23 18 30 1 0 $MAX_BLOCKS $TARGET_ERRORS $N_PRB $N_RX $N_LAYERS $CH_TYPE $DELAY_SPREAD $DOPPLER $DMRS_ADD_POS $DMRS_DUR $MCS_TABLE $RESULTS/cpp_ls_mcs23.csv

echo "=== Running Sionna Ideal CE ==="
python scripts/sionna_tdl_crosscheck.py --mcs 6 --snr-start 0 --snr-end 10 --snr-step 1 \
    --max-blocks $MAX_BLOCKS --target-errors $TARGET_ERRORS --batch-size $BATCH \
    --ce-type ideal --output $RESULTS/sionna_ideal_mcs6.csv
python scripts/sionna_tdl_crosscheck.py --mcs 12 --snr-start 8 --snr-end 18 --snr-step 1 \
    --max-blocks $MAX_BLOCKS --target-errors $TARGET_ERRORS --batch-size $BATCH \
    --ce-type ideal --output $RESULTS/sionna_ideal_mcs12.csv
python scripts/sionna_tdl_crosscheck.py --mcs 23 --snr-start 15 --snr-end 27 --snr-step 1 \
    --max-blocks $MAX_BLOCKS --target-errors $TARGET_ERRORS --batch-size $BATCH \
    --ce-type ideal --output $RESULTS/sionna_ideal_mcs23.csv

echo "=== Running Sionna LS CE ==="
python scripts/sionna_tdl_crosscheck.py --mcs 6 --snr-start 3 --snr-end 13 --snr-step 1 \
    --max-blocks $MAX_BLOCKS --target-errors $TARGET_ERRORS --batch-size $BATCH \
    --ce-type ls --output $RESULTS/sionna_ls_mcs6.csv
python scripts/sionna_tdl_crosscheck.py --mcs 12 --snr-start 11 --snr-end 21 --snr-step 1 \
    --max-blocks $MAX_BLOCKS --target-errors $TARGET_ERRORS --batch-size $BATCH \
    --ce-type ls --output $RESULTS/sionna_ls_mcs12.csv
python scripts/sionna_tdl_crosscheck.py --mcs 23 --snr-start 18 --snr-end 30 --snr-step 1 \
    --max-blocks $MAX_BLOCKS --target-errors $TARGET_ERRORS --batch-size $BATCH \
    --ce-type ls --output $RESULTS/sionna_ls_mcs23.csv

echo "=== All simulations complete! ==="
ls -la results/
