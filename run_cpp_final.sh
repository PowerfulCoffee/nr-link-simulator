#!/bin/bash
set -e
cd /workspace/nr-link-simulator
RESULTS=results/tdl_crosscheck_4x4
mkdir -p $RESULTS
BIN=./build/examples/pdsch_bler_tdl_mimo

BLOCKS=30
ERRORS=15

echo "=== C++ Ideal CE (1dB step, waterfall region) ==="

echo "QPSK MCS4 Ideal CE (-2 to 6 dB)..."
$BIN 4 -2 6 1 1 $BLOCKS $ERRORS 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ideal_mcs4_1db.csv

echo "16QAM MCS10 Ideal CE (6 to 16 dB)..."
$BIN 10 6 16 1 1 $BLOCKS $ERRORS 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ideal_mcs10_1db.csv

echo "64QAM MCS16 Ideal CE (14 to 26 dB)..."
$BIN 16 14 26 1 1 $BLOCKS $ERRORS 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ideal_mcs16_1db.csv

echo "256QAM MCS24 Ideal CE (22 to 32 dB)..."
$BIN 24 22 32 1 1 $BLOCKS $ERRORS 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ideal_mcs24_1db.csv

echo ""
echo "=== C++ LS CE (1dB step) ==="

echo "QPSK MCS4 LS CE (0 to 10 dB)..."
$BIN 4 0 10 1 0 $BLOCKS $ERRORS 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ls_mcs4.csv

echo "16QAM MCS10 LS CE (10 to 22 dB)..."
$BIN 10 10 22 1 0 $BLOCKS $ERRORS 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ls_mcs10.csv

echo "64QAM MCS16 LS CE (18 to 32 dB)..."
$BIN 16 18 32 1 0 $BLOCKS $ERRORS 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ls_mcs16.csv

echo "256QAM MCS24 LS CE (26 to 40 dB)..."
$BIN 24 26 40 1 0 $BLOCKS $ERRORS 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ls_mcs24.csv

echo ""
echo "All C++ simulations done!"
