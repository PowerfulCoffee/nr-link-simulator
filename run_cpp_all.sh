#!/bin/bash
set -e
cd /workspace/nr-link-simulator
RESULTS=results/tdl_crosscheck_4x4
mkdir -p $RESULTS

echo "Running C++ QPSK MCS4..."
./build/examples/pdsch_bler_tdl_mimo 4 -4 12 2 1 50 20 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ideal_mcs4.csv > $RESULTS/cpp_ideal_mcs4.log 2>&1

echo "Running C++ 16QAM MCS10..."
./build/examples/pdsch_bler_tdl_mimo 10 4 24 2 1 50 20 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ideal_mcs10.csv > $RESULTS/cpp_ideal_mcs10.log 2>&1

echo "Running C++ 64QAM MCS16..."
./build/examples/pdsch_bler_tdl_mimo 16 10 30 2 1 50 20 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ideal_mcs16.csv > $RESULTS/cpp_ideal_mcs16.log 2>&1

echo "Running C++ 256QAM MCS24..."
./build/examples/pdsch_bler_tdl_mimo 24 16 36 2 1 50 20 30 4 4 1 100 5 1 0 2 $RESULTS/cpp_ideal_mcs24.csv > $RESULTS/cpp_ideal_mcs24.log 2>&1

echo "All C++ simulations done!"
