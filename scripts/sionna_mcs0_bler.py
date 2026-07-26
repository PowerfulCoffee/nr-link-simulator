#!/usr/bin/env python3
"""Quick Sionna PDSCH BLER test for MCS 0 comparison with C++"""
import sys
sys.path.insert(0, '/workspace/sionna/src')
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import numpy as np
import torch
from sionna.phy.nr import PDSCHConfig, PDSCHTransmitter, PDSCHReceiver
from sionna.phy.channel import AWGN, OFDMModulator, OFDMDemodulator
from sionna.phy.mimo import StreamManagement
from sionna.phy.ofdm import ResourceGrid, ResourceGridMapper, LSChannelEstimator, LMMSEEqualizer, RemoveNulledScs
from sionna.phy import Block

# PDSCH parameters matching C++ test_bler
mcs_index = 0
n_rb = 25
n_tx = 1
n_rx = 1
n_layers = 1
num_iterations = 20
max_blocks = 100
target_errors = 50

# MCS 0: QPSK, rate 0.117188
modulation = 'qpsk'
target_rate = 0.117188

print(f"=== Sionna PDSCH BLER: MCS {mcs_index} ({modulation}), R={target_rate:.4f} ===")
print(f"PRBs: {n_rb}, TX: {n_tx}, RX: {n_rx}, Layers: {n_layers}")

# PDSCH configuration
config = PDSCHConfig()
config.num_layers = n_layers
config.mcs_index = mcs_index
config.mcs_table_index = 1
config.n_rb = n_rb
config.scs = 15
config.dmrs.config_type = 1
config.dmrs.length = 1
config.dmrs.additional_position = 1
config.dmrs.num_cdm_groups_without_data = 1

# Create transmitter and receiver
pusch_transmitter = PDSCHTransmitter(config)
rx_config = PDSCHConfig()
rx_config.num_layers = n_layers
rx_config.mcs_index = mcs_index
rx_config.mcs_table_index = 1
rx_config.n_rb = n_rb
rx_config.scs = 15
rx_config.dmrs.config_type = 1
rx_config.dmrs.length = 1
rx_config.dmrs.additional_position = 1
rx_config.dmrs.num_cdm_groups_without_data = 1

pusch_receiver = PDSCHReceiver(rx_config, num_iter=num_iterations)

# SNR points to test
snr_points = [-4, -3.5, -3, -2.5, -2]

channel = AWGN()

for snr_db in snr_points:
    snr_lin = 10.0**(snr_db/10.0)
    n_errors = 0
    n_blocks = 0
    
    for block in range(max_blocks):
        bsz = 1
        # Generate random transport block
        # First get TB size from transmitter
        tb, cb_len = pusch_transmitter(batch_size=bsz)
        
        # The transmitter returns (b, c) where b is TB bits, c is codeword bits?
        # Let's transmit properly
        x, b, cw, tb_info = pusch_transmitter(batch_size=bsz, return_tb=True)
        
        # x is the resource grid? Let's check
        # Actually let me just check the shape
        
        # Actually let's do this differently - use the full PUSCHReceiver
        # Let's set a fixed noise level
        pass
    
    # Let me just check how to properly use the API
    break

# Let me look at how to properly use PDSCHTransmitter/Receiver
print("\n--- Checking API ---")
try:
    import inspect
    sig = inspect.signature(PDSCHTransmitter.forward)
    print(f"PDSCHTransmitter.forward signature: {sig}")
except Exception as e:
    print(f"Error: {e}")
