"""
Sionna PUSCH BLER simulation for comparison with C++ nr-link-simulator.
Single-port, SISO, AWGN channel (freq domain), LS channel estimation.
"""
import sys
sys.path.insert(0, '/workspace/sionna/src')

import numpy as np
import torch

from sionna.phy.nr import PUSCHConfig, PUSCHTransmitter, PUSCHReceiver

def run_sionna_bler(mcs_index=5, n_size_bwp=25,
                    sinr_start=-4, sinr_end=6, sinr_step=1.0,
                    max_blocks=200, target_errors=30,
                    batch_size=20):
    
    pusch_config = PUSCHConfig()
    pusch_config.n_size_bwp = n_size_bwp
    pusch_config.num_antenna_ports = 1
    pusch_config.num_layers = 1
    
    pusch_config.dmrs.config_type = 1
    pusch_config.dmrs.length = 1
    pusch_config.dmrs.additional_position = 0
    pusch_config.dmrs.num_cdm_groups_without_data = 1
    pusch_config.dmrs.dmrs_port_set = [0]
    
    pusch_config.tb.mcs_index = mcs_index
    pusch_config.tb.mcs_table = 1
    
    pusch_transmitter = PUSCHTransmitter(pusch_config, output_domain="freq")
    
    rg = pusch_transmitter.resource_grid
    print(f"Resource grid: {rg.num_ofdm_symbols} symbols, {rg.fft_size} FFT size, {rg.num_effective_subcarriers} effective SCs")
    print(f"DMRS symbols: {rg.pilot_pattern.num_pilot_symbols} pilot symbols")
    sys.stdout.flush()
    
    pusch_receiver = PUSCHReceiver(
        pusch_transmitter,
        channel_estimator=None,  # None = LS estimator (default)
        input_domain="freq",
        return_tb_crc_status=True,
        disable_inter_stream_sic=True
    )
    
    results = []
    
    sinr_db = sinr_start
    while sinr_db <= sinr_end + 1e-9:
        sinr_lin = 10**(sinr_db / 10.0)
        no = 1.0 / sinr_lin
        
        n_blocks = 0
        n_errors = 0
        
        while n_blocks < max_blocks:
            bsz = min(batch_size, max_blocks - n_blocks)
            
            x, b = pusch_transmitter(bsz)
            # x: [batch_size, num_tx=1, num_tx_ant=1, num_ofdm_symbols, fft_size]
            
            noise = (torch.randn_like(x) + 1j*torch.randn_like(x)) * np.sqrt(no/2.0)
            y = x + noise
            
            no_t = torch.tensor(no, dtype=torch.float32)
            
            b_hat, crc_status = pusch_receiver(y, no_t)
            
            for i in range(bsz):
                n_blocks += 1
                if not crc_status[i, 0]:
                    n_errors += 1
            
            if n_errors >= target_errors and n_blocks >= 30:
                break
            if n_errors == 0 and n_blocks >= 100:
                break
        
        bler = n_errors / n_blocks if n_blocks > 0 else 1.0
        results.append((sinr_db, n_blocks, n_errors, bler))
        print(f"  SNR={sinr_db:5.1f} dB:  blocks={n_blocks:5d}, errors={n_errors:5d}, BLER={bler:.4f}")
        sys.stdout.flush()
        
        sinr_db += sinr_step
    
    return results


if __name__ == "__main__":
    print("="*60)
    print("Sionna PUSCH BLER Simulation (AWGN SISO, LS CE)")
    print("="*60)
    print()
    
    results = run_sionna_bler(
        mcs_index=5,
        n_size_bwp=25,
        sinr_start=-4,
        sinr_end=6,
        sinr_step=1.0,
        max_blocks=200,
        target_errors=30,
        batch_size=20
    )
    
    print()
    print("="*60)
    print("Sionna BLER Summary")
    print("="*60)
    print(f"{'SNR(dB)':>10} {'Blocks':>10} {'Errors':>10} {'BLER':>12}")
    print("-"*44)
    for snr, nb, ne, bler in results:
        print(f"{snr:10.1f} {nb:10d} {ne:10d} {bler:12.6f}")
