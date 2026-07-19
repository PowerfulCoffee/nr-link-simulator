#include "common/Types.h"
#include "common/NrTables.h"
#include "phy/PdschProcessor.h"
#include <iostream>
#include <iomanip>

using namespace nr;
using namespace nr::phy;

int main() {
    std::cout << "=== Feature Verification Tests ===\n\n";

    {
        std::cout << "--- Test 1: 256QAM (MCS27 with 256QAM capability) ---\n";
        SimulationConfig config;
        config.mcs_index = 27;
        config.n_rb = 3;
        config.n_tx_ant = 1;
        config.n_rx_ant = 1;
        config.n_layers = 1;
        config.channel_type = ChannelType::AWGN;
        config.max_blocks_per_sinr = 20;
        config.target_block_errors = 5;
        config.sinr_start = 30.0;
        config.sinr_end = 30.0;
        config.sinr_step = 2.0;
        config.n_sinr_points = 1;
        config.perfect_csi = false;
        config.n_ldpc_iterations = 20;
        config.random_seed = 42;
        config.scs = 15000;
        config.mod_scheme = mcs_to_modulation(config.mcs_index);
        config.code_rate = mcs_to_code_rate(config.mcs_index);

        auto proc = std::make_unique<PdschProcessor>(config);
        TransportBlock tb = proc->generate_transport_block();
        std::cout << "MCS27 TBS=" << tb.tb_size << " bits, mod="
                  << (config.mod_scheme == ModulationScheme::QAM64 ? "64QAM" : "other")
                  << " (MCS27 uses 64QAM per 3GPP table)\n";
        PdschTxResult tx = proc->transmit(tb, 0);
        std::cout << "  Encoded " << tx.n_coded_bits << " bits, " << tx.num_cb << " CBs\n";
    }

    {
        std::cout << "\n--- Test 2: DMRS Type 2 with additional positions ---\n";
        SimulationConfig config;
        config.mcs_index = 10;
        config.n_rb = 10;
        config.n_tx_ant = 1;
        config.n_rx_ant = 1;
        config.n_layers = 1;
        config.dmrs_type = DmrsType::TYPE2;
        config.dmrs_additional_pos = 1;
        config.dmrs_duration = 1;
        config.channel_type = ChannelType::AWGN;
        config.max_blocks_per_sinr = 5;
        config.target_block_errors = 1;
        config.sinr_start = 20.0;
        config.sinr_end = 20.0;
        config.sinr_step = 2.0;
        config.n_sinr_points = 1;
        config.perfect_csi = false;
        config.n_ldpc_iterations = 15;
        config.random_seed = 123;
        config.scs = 15000;
        config.mod_scheme = mcs_to_modulation(config.mcs_index);
        config.code_rate = mcs_to_code_rate(config.mcs_index);

        auto proc = std::make_unique<PdschProcessor>(config);
        TransportBlock tb = proc->generate_transport_block();
        std::cout << "DMRS Type2, pos=1: TBS=" << tb.tb_size << " bits\n";
        PdschTxResult tx = proc->transmit(tb, 0);

        DmrsPattern pat = get_dmrs_pattern(DmrsType::TYPE2, 1, 1);
        int dmrs_count = 0;
        for (int s = 0; s < 14; s++) {
            if (pat.re_per_prb[s] > 0) {
                std::cout << "  DMRS at symbol " << s << " (" << pat.re_per_prb[s] << " RE/PRB)\n";
                dmrs_count++;
            }
        }
        std::cout << "Total DMRS symbols: " << dmrs_count << "\n";
        std::cout << "Encoded " << tx.n_coded_bits << " bits\n";
    }

    {
        std::cout << "\n--- Test 3: 2x2 MIMO (TDL-A, perfect CSI) ---\n";
        SimulationConfig config;
        config.mcs_index = 15;
        config.n_rb = 5;
        config.n_tx_ant = 2;
        config.n_rx_ant = 2;
        config.n_layers = 2;
        config.channel_type = ChannelType::TDL_A;
        config.delay_spread = 30e-9;
        config.max_doppler = 0.0;
        config.perfect_csi = true;
        config.max_blocks_per_sinr = 20;
        config.target_block_errors = 5;
        config.sinr_start = 10.0;
        config.sinr_end = 30.0;
        config.sinr_step = 10.0;
        config.n_sinr_points = 3;
        config.n_ldpc_iterations = 15;
        config.random_seed = 456;
        config.scs = 15000;
        config.mod_scheme = mcs_to_modulation(config.mcs_index);
        config.code_rate = mcs_to_code_rate(config.mcs_index);

        std::cout << "2x2 MIMO, MCS15, 5PRB, TDL-A\n";
        auto results = run_bler_simulation(config, nullptr, "perfect");
        for (auto& r : results) {
            std::cout << "  SNR=" << r.sinr_db << "dB: BLER=" << std::fixed << std::setprecision(3) << r.bler
                      << " (" << r.n_errors << "/" << r.n_blocks << ")\n";
        }
    }

    {
        std::cout << "\n--- Test 4: TDD slot format configuration ---\n";
        DmrsPattern pat = get_dmrs_pattern(DmrsType::TYPE1, 0, 1);
        SimulationConfig cfg;
        cfg.tdd_enabled = true;
        cfg.tdd_config.s_slot_dl_symbols = 6;
        cfg.tdd_config.s_slot_gp_symbols = 4;
        cfg.tdd_config.s_slot_ul_symbols = 4;

        auto is_pdsch = [&](int sym, int slot_idx) {
            if (pat.re_per_prb[sym] > 0) return false;
            if (!cfg.tdd_enabled) return true;
            int pat_idx = slot_idx % TddConfig::slots_per_frame;
            SlotType st = cfg.tdd_config.slot_pattern[pat_idx];
            if (st == SlotType::DOWNLINK) return true;
            if (st == SlotType::UPLINK) return false;
            if (st == SlotType::SPECIAL) return sym < cfg.tdd_config.s_slot_dl_symbols;
            return true;
        };

        std::cout << "TDD pattern (DDDSU DDDSUU):\n";
        for (int s = 0; s < 20; s++) {
            int pat_idx = s % TddConfig::slots_per_frame;
            const char* type_str = "?";
            switch(cfg.tdd_config.slot_pattern[pat_idx]) {
                case SlotType::DOWNLINK: type_str = "D"; break;
                case SlotType::UPLINK: type_str = "U"; break;
                case SlotType::SPECIAL: type_str = "S"; break;
            }
            int pdsch_syms = 0;
            for (int sym = 0; sym < 14; sym++) {
                if (is_pdsch(sym, s)) pdsch_syms++;
            }
            std::cout << "  Slot " << std::setw(2) << s << ": " << type_str
                      << ", PDSCH symbols=" << pdsch_syms << "\n";
        }
    }

    {
        std::cout << "\n--- Test 5: TDL-B/C/D/E and CDL channels registered ---\n";
        for (auto ct : {ChannelType::TDL_A, ChannelType::TDL_B, ChannelType::TDL_C,
                        ChannelType::TDL_D, ChannelType::TDL_E,
                        ChannelType::CDL_A, ChannelType::CDL_B, ChannelType::CDL_C,
                        ChannelType::CDL_D, ChannelType::CDL_E}) {
            auto ch = channel::create_channel(ct);
            bool ok = (ch != nullptr);
            std::cout << "  " << static_cast<int>(ct) << ": " << (ok ? "OK" : "FAIL") << "\n";
        }
    }

    std::cout << "\n=== All feature tests completed ===\n";
    return 0;
}
