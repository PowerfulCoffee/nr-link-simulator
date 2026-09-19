#!/bin/bash
set -e
cd /workspace/nr-link-simulator
RESULTS=results/tdl_crosscheck_4x4
mkdir -p $RESULTS
BIN=/workspace/nr-link-simulator/build/examples/pdsch_bler_tdl_mimo

BLOCKS=30
ERRORS=15
STEP=1

echo "=== Running C++ LS CE simulations ==="

echo "QPSK MCS4 LS CE (-2 to 8 dB)..."
$BIN 4 -2 8 $STEP 0 $BLOCKS $ERRORS 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ls_mcs4.csv

echo "16QAM MCS10 LS CE (6 to 18 dB)..."
$BIN 10 6 18 $STEP 0 $BLOCKS $ERRORS 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ls_mcs10.csv

echo "64QAM MCS16 LS CE (14 to 26 dB)..."
$BIN 16 14 26 $STEP 0 $BLOCKS $ERRORS 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ls_mcs16.csv

echo "256QAM MCS24 LS CE (22 to 34 dB)..."
$BIN 24 22 34 $STEP 0 $BLOCKS $ERRORS 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ls_mcs24.csv

echo "C++ LS CE done!"
