#!/bin/bash
set -e
cd /workspace/nr-link-simulator/build

CPP=./examples/pdsch_bler_tdl_mimo
RESULTS=../results/tdl_crosscheck

# Config: 30PRB, 4Rx, 4 layers, TDL-A DS=100ns, fd=5Hz, DMRS Type1 add-pos=1 (sym 2&7, duration=0), MCS table 2
# Args: mcs sinr_start sinr_end sinr_step perfect_csi max_blocks target_errors n_rb n_rx_ant n_layers ch_type ds fd dmrs_add_pos dmrs_duration mcs_table out_file

N_RB=30
N_RX=4
N_LAYERS=4
CH_TYPE=1
DS=100e-9
FD=5
DMRS_ADD=1
DMRS_DUR=0
MCS_TBL=2
MAX_BLOCKS=100
TARGET_ERR=20

echo "=== Running C++ QPSK (MCS4) Ideal CE ==="
$CPP 4 -6 18 2 1 $MAX_BLOCKS $TARGET_ERR $N_RB $N_RX $N_LAYERS $CH_TYPE $DS $FD $DMRS_ADD $DMRS_DUR $MCS_TBL $RESULTS/cpp_ideal_mcs4_table2_30prb_4x4.csv

echo "=== Running C++ QPSK (MCS4) LS CE ==="
$CPP 4 -6 18 2 0 $MAX_BLOCKS $TARGET_ERR $N_RB $N_RX $N_LAYERS $CH_TYPE $DS $FD $DMRS_ADD $DMRS_DUR $MCS_TBL $RESULTS/cpp_ls_mcs4_table2_30prb_4x4.csv
