/*
 * Scenario-driven PDSCH E2E BLER sweep for the web platform.
 * Reads simulation parameters from CLI args (mapped by web/server.py from
 * a scenario JSON), runs the full TX->channel->RX chain per SNR point,
 * prints [PROGRESS] lines for live log streaming and writes a CSV.
 */
#include "common/Types.h"
#include "common/NrTables.h"
#include "channel/ChannelModels.h"
#include "phy/PdschProcessor.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <filesystem>
#include <algorithm>

using namespace nr;
using namespace nr::phy;
using namespace std::chrono;
namespace fs = std::filesystem;

namespace {

struct WebScenario {
    std::vector<int> mcs_list = {0, 9, 16, 27};
    int mcs_table = 1;
    int n_rb = 3;
    int n_tx = 1;
    int n_rx = 1;
    int n_layers = 1;
    double scs = 15000.0;
    bool perfect_csi = false;
    ChannelType channel = ChannelType::AWGN;
    double delay_spread_ns = 30.0;
    double doppler_hz = 5.0;
    bool los = false;
    int max_blocks = 2000;
    int target_errors = 60;
    int ldpc_iters = 25;
    uint64_t seed = 123;
    // DMRS
    int dmrs_type = 1;      // 1 or 2
    int dmrs_add_pos = 0;
    int dmrs_duration = 1;
    // scrambling (TS 38.211 7.3.1.1)
    int n_rnti = 1;
    int n_id = 0;
    int n_id_dmrs = 0;
    // precoder
    std::string precoder = "identity";   // identity | type1
    int csi_n1 = 1, csi_n2 = 1;
    int csi_i1_phi = 0, csi_i11 = 0, csi_i12 = 0, csi_i2 = 0;
    // SNR sweep; both snr_start/snr_end <= -900 -> auto range per MCS
    double snr_start = -999.0, snr_end = -999.0, snr_step = 1.0;
    std::string out_csv = "results/runs/unnamed/bler.csv";
};

std::vector<int> parse_int_list(const std::string& s) {
    std::vector<int> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (!tok.empty()) out.push_back(std::stoi(tok));
    }
    return out;
}

ChannelType parse_channel(const std::string& s) {
    if (s == "AWGN") return ChannelType::AWGN;
    if (s == "TDL_A") return ChannelType::TDL_A;
    if (s == "TDL_B") return ChannelType::TDL_B;
    if (s == "TDL_C") return ChannelType::TDL_C;
    if (s == "TDL_D") return ChannelType::TDL_D;
    if (s == "TDL_E") return ChannelType::TDL_E;
    if (s == "CDL_A") return ChannelType::CDL_A;
    if (s == "CDL_B") return ChannelType::CDL_B;
    if (s == "CDL_C") return ChannelType::CDL_C;
    if (s == "CDL_D") return ChannelType::CDL_D;
    if (s == "CDL_E") return ChannelType::CDL_E;
    std::cerr << "[WARN] unknown channel '" << s << "', using AWGN" << std::endl;
    return ChannelType::AWGN;
}

WebScenario parse_args(int argc, char** argv) {
    WebScenario sc;
    auto get = [&](int& i) -> std::string {
        if (i + 1 >= argc) {
            std::cerr << "[ERROR] missing value for " << argv[i] << std::endl;
            exit(2);
        }
        return argv[++i];
    };
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--mcs") sc.mcs_list = parse_int_list(get(i));
        else if (a == "--mcs_table") {
            std::string v = get(i);
            // tolerate symbolic names from web UI: qam64/qam256 (see Types.h mcs_table)
            if (v == "qam64" || v == "1") sc.mcs_table = 1;
            else if (v == "qam256" || v == "2") sc.mcs_table = 2;
            else sc.mcs_table = std::stoi(v);
        }
        else if (a == "--n_rb") sc.n_rb = std::stoi(get(i));
        else if (a == "--n_tx") sc.n_tx = std::stoi(get(i));
        else if (a == "--n_rx") sc.n_rx = std::stoi(get(i));
        else if (a == "--n_layers") sc.n_layers = std::stoi(get(i));
        else if (a == "--scs") sc.scs = std::stod(get(i));
        else if (a == "--csi") sc.perfect_csi = (get(i) == "perfect");
        else if (a == "--channel") sc.channel = parse_channel(get(i));
        else if (a == "--delay_spread_ns") sc.delay_spread_ns = std::stod(get(i));
        else if (a == "--doppler") sc.doppler_hz = std::stod(get(i));
        else if (a == "--los") {
            // NOTE: must call get(i) exactly once -- "get(i)=="1" || get(i)=="true"
            // would consume the NEXT argument too (|| evaluates both sides
            // when the first fails), swallowing the following flag's value.
            std::string v = get(i);
            sc.los = (v == "1" || v == "true");
        }
        else if (a == "--max_blocks") sc.max_blocks = std::stoi(get(i));
        else if (a == "--target_errors") sc.target_errors = std::stoi(get(i));
        else if (a == "--ldpc_iters") sc.ldpc_iters = std::stoi(get(i));
        else if (a == "--seed") sc.seed = static_cast<uint64_t>(std::stoull(get(i)));
        else if (a == "--dmrs_type") sc.dmrs_type = std::stoi(get(i));
        else if (a == "--dmrs_addpos") sc.dmrs_add_pos = std::stoi(get(i));
        else if (a == "--dmrs_duration") sc.dmrs_duration = std::stoi(get(i));
        else if (a == "--rnti") sc.n_rnti = std::stoi(get(i));
        else if (a == "--n_id") sc.n_id = std::stoi(get(i));
        else if (a == "--n_id_dmrs") sc.n_id_dmrs = std::stoi(get(i));
        else if (a == "--precoder") sc.precoder = get(i);
        else if (a == "--csi_n1") sc.csi_n1 = std::stoi(get(i));
        else if (a == "--csi_n2") sc.csi_n2 = std::stoi(get(i));
        else if (a == "--csi_i1_phi") sc.csi_i1_phi = std::stoi(get(i));
        else if (a == "--csi_i11") sc.csi_i11 = std::stoi(get(i));
        else if (a == "--csi_i12") sc.csi_i12 = std::stoi(get(i));
        else if (a == "--csi_i2") sc.csi_i2 = std::stoi(get(i));
        else if (a == "--snr_start") sc.snr_start = std::stod(get(i));
        else if (a == "--snr_end") sc.snr_end = std::stod(get(i));
        else if (a == "--snr_step") sc.snr_step = std::stod(get(i));
        else if (a == "--out") sc.out_csv = get(i);
        else std::cerr << "[WARN] unknown arg " << a << std::endl;
    }
    return sc;
}

// Per-MCS SNR range heuristic (same logic as pdsch_bler_baseline)
void auto_snr_range(int qm, double R, double& start, double& end) {
    if (qm == 2) {
        if (R < 0.20)      { start = -10.0; end = 6.0; }
        else if (R < 0.45) { start = -6.0;  end = 10.0; }
        else               { start = -4.0;  end = 12.0; }
    } else if (qm == 4) {
        if (R < 0.50)      { start = 0.0;   end = 16.0; }
        else               { start = 2.0;   end = 18.0; }
    } else if (qm == 6) {
        if (R < 0.60)      { start = 4.0;   end = 22.0; }
        else if (R < 0.76) { start = 8.0;   end = 24.0; }
        else               { start = 10.0;  end = 26.0; }
    } else {
        if (R < 0.75)      { start = 12.0;  end = 28.0; }
        else if (R < 0.87) { start = 16.0;  end = 32.0; }
        else               { start = 18.0;  end = 34.0; }
    }
}

} // namespace


int main(int argc, char** argv) {
    WebScenario sc = parse_args(argc, argv);

    std::cout << "[INFO] === PDSCH E2E scenario run (web runner) ===" << std::endl;
    std::cout << "[INFO] channel=" << (sc.channel == ChannelType::AWGN ? "AWGN" : "TDL/CDL")
              << " csi=" << (sc.perfect_csi ? "perfect" : "LS")
              << " n_rb=" << sc.n_rb << " tx=" << sc.n_tx << " rx=" << sc.n_rx
              << " layers=" << sc.n_layers
              << " precoder=" << sc.precoder
              << " dmrs_type=" << sc.dmrs_type
              << " scs=" << sc.scs << "Hz" << std::endl;
    std::cout << "[INFO] mcs_list=";
    for (size_t i = 0; i < sc.mcs_list.size(); i++) {
        std::cout << (i ? "," : "") << sc.mcs_list[i];
    }
    std::cout << " max_blocks=" << sc.max_blocks
              << " target_errors=" << sc.target_errors
              << " ldpc_iters=" << sc.ldpc_iters
              << " seed=" << sc.seed << std::endl;

    fs::path out_path(sc.out_csv);
    if (out_path.has_parent_path()) {
        fs::create_directories(out_path.parent_path());
    }
    std::ofstream csv(sc.out_csv);
    if (!csv.is_open()) {
        std::cerr << "[ERROR] cannot open output CSV: " << sc.out_csv << std::endl;
        return 3;
    }
    csv << "mcs,qm,rate,tbs,snr_db,bler,n_blocks,n_errors\n";

    // Pre-compute total SNR points for progress reporting
    const int n_re_per_prb = 156;
    std::vector<int> tbs_list(sc.mcs_list.size());
    std::vector<double> snr_starts(sc.mcs_list.size()), snr_ends(sc.mcs_list.size());
    int total_pts = 0;
    for (size_t m = 0; m < sc.mcs_list.size(); m++) {
        int mcs = sc.mcs_list[m];
        int qm = mcs_to_bits_per_symbol(mcs, sc.mcs_table);
        double R = mcs_to_code_rate(mcs, sc.mcs_table);
        tbs_list[m] = calculate_tbs(sc.n_rb, n_re_per_prb, qm, sc.n_layers, R);
        double s0, s1;
        if (sc.snr_start > -900.0 && sc.snr_end > -900.0) {
            s0 = sc.snr_start;
            s1 = sc.snr_end;
        } else {
            auto_snr_range(qm, R, s0, s1);
        }
        snr_starts[m] = s0;
        snr_ends[m] = s1;
        total_pts += static_cast<int>((s1 - s0) / sc.snr_step) + 1;
    }
    std::cout << "[INFO] total_snr_points=" << total_pts
              << " out=" << sc.out_csv << std::endl;

    auto t0 = high_resolution_clock::now();
    int done = 0;

    for (size_t m = 0; m < sc.mcs_list.size(); m++) {
        int mcs = sc.mcs_list[m];
        int qm = mcs_to_bits_per_symbol(mcs, sc.mcs_table);
        double R = mcs_to_code_rate(mcs, sc.mcs_table);
        int tbs = tbs_list[m];

        SimulationConfig cfg;
        cfg.mcs_index = mcs;
        cfg.mcs_table = sc.mcs_table;
        cfg.mod_scheme = mcs_to_modulation(mcs, sc.mcs_table);
        cfg.code_rate = R;
        cfg.n_rb = sc.n_rb;
        cfg.n_tx_ant = sc.n_tx;
        cfg.n_rx_ant = sc.n_rx;
        cfg.n_layers = sc.n_layers;
        cfg.scs = sc.scs;
        cfg.channel_type = sc.channel;
        cfg.delay_spread = sc.delay_spread_ns * 1e-9;
        cfg.max_doppler = sc.doppler_hz;
        cfg.enable_los = sc.los;
        cfg.perfect_csi = sc.perfect_csi;
        cfg.dmrs_type = (sc.dmrs_type == 2) ? DmrsType::TYPE2 : DmrsType::TYPE1;
        cfg.dmrs_additional_pos = sc.dmrs_add_pos;
        cfg.dmrs_duration = sc.dmrs_duration;
        cfg.n_rnti = sc.n_rnti;
        cfg.n_id_pdsch = sc.n_id;
        cfg.n_id_dmrs = sc.n_id_dmrs;
        cfg.precoder_mode = (sc.precoder == "type1") ? PrecoderMode::TYPE1_CSI
                                                     : PrecoderMode::IDENTITY;
        cfg.csi_n1 = sc.csi_n1;
        cfg.csi_n2 = sc.csi_n2;
        cfg.csi_i1_phi = sc.csi_i1_phi;
        cfg.csi_i11 = sc.csi_i11;
        cfg.csi_i12 = sc.csi_i12;
        cfg.csi_i2 = sc.csi_i2;
        cfg.max_blocks_per_sinr = sc.max_blocks;
        cfg.target_block_errors = sc.target_errors;
        cfg.n_ldpc_iterations = sc.ldpc_iters;
        cfg.early_termination = true;
        cfg.random_seed = sc.seed;
        cfg.sinr_start = snr_starts[m];
        cfg.sinr_end = snr_ends[m];
        cfg.sinr_step = sc.snr_step;
        cfg.n_sinr_points =
            static_cast<int>((snr_ends[m] - snr_starts[m]) / sc.snr_step) + 1;

        std::cout << "[INFO] --- MCS " << mcs << " (Qm=" << qm << ", R=" << std::fixed
                  << std::setprecision(4) << R << ", TBS=" << tbs
                  << "), SNR " << snr_starts[m] << ".." << snr_ends[m]
                  << " dB step " << sc.snr_step << " ---" << std::endl;

        PdschProcessor processor(cfg);
        processor.set_seed(sc.seed + static_cast<uint64_t>(mcs));

        for (int i = 0; i < cfg.n_sinr_points; i++) {
            double snr = snr_starts[m] + i * sc.snr_step;
            if (snr > snr_ends[m] + 1e-9) break;

            BlerResult res;
            processor.process_single_snr_point(snr, res);

            csv << mcs << "," << qm << "," << std::fixed << std::setprecision(4) << R
                << "," << tbs << "," << std::setprecision(2) << res.sinr_db
                << "," << std::setprecision(6) << res.bler << ","
                << res.n_blocks << "," << res.n_errors << "\n";
            csv.flush();

            done++;
            std::cout << "[PROGRESS] " << done << "/" << total_pts
                      << " mcs=" << mcs << " snr=" << std::setprecision(2) << res.sinr_db
                      << " bler=" << std::setprecision(6) << res.bler
                      << " errors=" << res.n_errors << "/" << res.n_blocks
                      << std::endl;
        }
    }

    double sec = duration_cast<milliseconds>(high_resolution_clock::now() - t0).count() / 1000.0;
    std::cout << "[INFO] DONE points=" << done << " elapsed_s=" << std::fixed
              << std::setprecision(1) << sec << " csv=" << sc.out_csv << std::endl;
    return 0;
}

