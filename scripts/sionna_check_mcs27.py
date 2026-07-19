import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import torch
import numpy as np
from sionna.nr import PDSCHConfig, PDSCHTransmitter, PDSCHReceiver
from sionna.channel import AWGN

mcs_index = 27
n_rb = 3
n_layers = 1
n_iter = 20

pdsch_cfg = PDSCHConfig()
pdsch_cfg.bandwidth = n_rb
pdsch_cfg.mcs_index = mcs_index
pdsch_cfg.num_layers = n_layers

print(f"MCS {mcs_index}:")
print(f"  num_bits_per_symbol (Qm) = {pdsch_cfg.num_bits_per_symbol}")
print(f"  target_coderate = {pdsch_cfg.target_coderate}")
print(f"  tb_size = {pdsch_cfg.tb_size}")
print(f"  n_re/PRB = {pdsch_cfg.num_resource_elements_per_prb}")
print(f"  dmrs: type={pdsch_cfg.dmrs.config_type}, add_pos={pdsch_cfg.dmrs.additional_positions}")
