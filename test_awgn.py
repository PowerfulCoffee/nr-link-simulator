import torch
import numpy as np
from sionna.phy.nr import TBEncoder, TBDecoder, calculate_tb_size
from sionna.phy.mapping import Constellation, Mapper, Demapper

torch.manual_seed(42); np.random.seed(42)
qm=2; R=363/1024; n_prb=30; n_sc=360; n_ofdm_sym=14
n_layers=4; n_rx=4; n_tx=4
dmrs_symbols=[2,7]
pdsch_symbols=[s for s in range(n_ofdm_sym) if s not in dmrs_symbols]
n_pdsch=len(pdsch_symbols)

tb = calculate_tb_size(modulation_order=qm, target_coderate=R, num_prbs=n_prb,
    num_ofdm_symbols=n_ofdm_sym, num_dmrs_per_prb=24, num_layers=n_layers, verbose=False)
tb_size=int(tb[0]); cw_length=int(torch.sum(tb[5]))
G = n_sc*n_pdsch*qm*n_layers

const = Constellation("qam", num_bits_per_symbol=qm, normalize=True)
map_ = Mapper(constellation=const)
demap = Demapper("app", constellation=const)
enc = TBEncoder(tb_size, G, R, qm, num_layers=n_layers, channel_type="PDSCH", verbose=False)
dec = TBDecoder(enc, num_bp_iter=20)

for snr_db in [-2, 0, 2, 4, 6, 8]:
    noise_var = 1.0/(10**(snr_db/10))
    n_err = 0; n_blk = 10
    for _ in range(n_blk):
        bits = torch.randint(0,2,(1,tb_size))
        cw = enc(bits)
        syms = map_(cw)
        syms_grid = syms.reshape(1, n_layers, n_pdsch, n_sc)
        tx_grid = torch.zeros(1,n_layers,n_ofdm_sym,n_sc,dtype=torch.complex64)
        for i,s in enumerate(pdsch_symbols):
            tx_grid[:,:,s,:] = syms_grid[:,:,i,:]
        dmrs_val=1/np.sqrt(2)
        for l in range(n_layers):
            for sym in dmrs_symbols:
                for sc in range(l,n_sc,4):
                    tx_grid[:,l,sym,sc] = dmrs_val
        h = (torch.randn(1,n_ofdm_sym,n_sc,n_rx,n_tx)+1j*torch.randn(1,n_ofdm_sym,n_sc,n_rx,n_tx))/np.sqrt(2)
        rx = torch.zeros(1,n_rx,n_ofdm_sym,n_sc,dtype=torch.complex64)
        for sym in range(n_ofdm_sym):
            for sc in range(n_sc):
                rx[0,:,sym,sc] = h[0,sym,sc,:,:] @ tx_grid[0,:,sym,sc]
        rx = rx + (torch.randn_like(rx)+1j*torch.randn_like(rx))*np.sqrt(noise_var/2)
        # Ideal CE
        h_est = h.permute(0,3,1,2,4)
        rx_d = torch.zeros(1,n_rx,n_pdsch,n_sc,dtype=torch.complex64)
        h_d = torch.zeros(1,n_rx,n_pdsch,n_sc,n_layers,dtype=torch.complex64)
        for i,s in enumerate(pdsch_symbols):
            rx_d[:,:,i,:]=rx[:,:,s,:]
            h_d[:,:,i,:,:]=h_est[:,:,s,:,:]
        Ht = h_d.permute(0,2,3,1,4).reshape(n_pdsch*n_sc,n_rx,n_layers)
        yt = rx_d.permute(0,2,3,1).reshape(n_pdsch*n_sc,n_rx)
        xh = torch.zeros_like(yt[:,:n_layers])
        for k in range(Ht.shape[0]):
            xh[k] = torch.linalg.lstsq(Ht[k], yt[k].unsqueeze(1)).solution.squeeze(1)
        xh = xh.reshape(1,n_pdsch,n_sc,n_layers).permute(0,3,1,2)
        rx_syms = xh.reshape(1,n_layers*n_pdsch*n_sc)[:,:cw_length]
        llr = demap(rx_syms, no=noise_var)
        d, crc = dec(llr)
        if not crc[0] or not torch.equal(bits[0], d[0].long()):
            n_err+=1
    print(f"SNR={snr_db}dB: BLER={n_err/n_blk:.3f} ({n_err}/{n_blk})")
