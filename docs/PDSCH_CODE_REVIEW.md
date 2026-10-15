# PDSCH 全链路仿真代码评审报告

> 评审范围：`src/phy/*`、`src/common/NrTables.cpp`、`src/channel/*`、`include/phy/*`
> 评审方式：逐模块对照 3GPP TS 38.211 / 38.212 / 38.214 及 TR 38.901 审查
> 结论分级：✅ 正确 | ⚠️ 可用但与规范有偏差（收发自洽）| ❌ 错误或需修正

## 总体结论

全链路 TX→信道→RX 各环节**自洽**（收发两侧使用同一套非标准简化），因此 BLER 曲线在"相对比较"意义上有效（且已与 Sionna 在 AWGN/理想 CSI 下对齐）。但若目标是**协议一致性仿真**，速率匹配、DMRS、加扰 c_init、OFDM CP 等多处与规范存在实质性偏差。以下按环节列出。

---

## 1. 发射端 (TX)

### 1.1 TB 生成 / TBS 计算 (`NrTables.cpp: calculate_tbs`) — ✅ 基本正确
- NRE = min(156, 每PRB可用RE) ✔（38.214 §5.1.3.2）
- n_info 量化：≤3824 用 `max(24, 2^n·floor(A/2^n))`，>3824 用 `2^n·round((A−24)/2^n)` ✔
- TBS 表查找、大 TB 的 `8·C·ceil(...)−24` ✔
- **改进**：未扣除 PDSCH 预留 RE（本仿真无 CSI-RS/核心集开销，可接受）；`n_re_per_prb` 由 `count_pdsch_re_per_prb` 正确剔除 DMRS 符号 ✔

### 1.2 MCS 表 (`MCS_TABLE_1/2`) — ✅ 数值正确
- Table 1（Table 5.1.3.1-1）数值核对无误；Table 2（256QAM）索引 20–27 为 682.5…948 ✔
- **小问题**：Table 1/2 均应含保留索引（29/28-31）占位，`MAX_MCS_INDEX` 边界需核实，越界访问有回退到 QPSK 的兜底 ✔。

### 1.3 CRC (`CrcEncoder.cpp`) — ✅ 正确
- CRC24A `0x864CFB`、CRC16 `0x1021`，MSB-first 输出，初始值 0 ✔（38.212 §5.1/§7.2）
- **改进**：解码用"重新编码逐位比对"，功能正确但可改用余式校验提速；CRC24B 未实现（PDCCH 用，可不做）。

### 1.4 码块分割 (`compute_cb_segmentation`) — ✅ 基本正确
- B≤3824→1 个 CB；否则 `ceil((B+24)/8424)`，CB-CRC=24 ✔；K' 向上取整到 8 的倍数 ✔（38.212 §5.2.2）
- Zc 选择：从 8 组合法 Zc 中选满足 `zc·k_b ≥ K'` 的最小值 ✔
- **⚠️ 偏差**：E_r 采用 `ceil(G/C)` 取 qm 整数倍、末码块取余数；规范是 `E = Qm·floor(G/(Qm·C))`（前 C−1 个 CB）+ 末块余量。C=1 无影响，C>1 时分配与规范不同（自洽，但不可互通）。

### 1.5 LDPC 编码 (`LdpcCodec.cpp: LdpcEncoderImpl`) — ✅ 正确
- BG1 46×68 / BG2 42×52 ✔；打孔前 2·Zc 信息位 ✔；核心子矩阵（4·Zc）消元求 p1–p4、扩展位递推 ✔
- **⚠️ 缺失**：**奇偶位交织（parity interleaver，38.212 §5.4.2）未实现**——标准速率匹配的配套步骤，缺失后打孔/抽取模式与标准不同。
- **改进**：编码可用 QC 结构移位累加替代通用消元；encoder 无结果缓存（decoder 图有缓存）。

### 1.6 速率匹配 (`RateMatcher.cpp`) — ✅ 修正后符合规范（2026-09-20 复核修订）
> **修订**：本节最初的两条结论（k0 公式错误、缺 parity interleaver）经与 TS 38.212 v17.3.0 原文逐条核对后**均不成立**，特此更正：
- **k0 公式原本就正确**：规范 Table 5.4.2.1-2 为 `k0 = floor(c·N_cb/(66·Zc))·Zc`（BG1，c=0/17/33/56）、`floor(c·N_cb/(50·Zc))·Zc`（BG2，c=0/13/25/43），与代码一致。此前"44/42"系笔误。
- **38.212 中不存在奇偶交织器（parity interleaver/Table 5.4.2.1-3）**——那是 LTE Turbo/其它系统的概念。LDPC 速率匹配只有两步：bit selection（循环缓冲 + 跳过 filler）与 5.4.2.2 位交织 `f_{i+j·Qm} = e_{i·E/Qm+j}`。代码中的 `apply_bicm_interleave` **正是 5.4.2.2 规范交织器**，不是多余步骤（已加 E%Qm 校验）。
- **Filler 处理偏差（已修复）**：规范要求 filler 位作为 `<NULL>` 保留在循环缓冲区内、bit selection 时跳过，k0 基于完整 Ncb 计算。原实现先剔除 filler 再用缩短的 Ncb 算 k0。修复后 `RateMatcher` 内部重建含 NULL 的完整缓冲（`NullSkipBuffer`），k0 按完整 Ncb 计算，选择时跳过 NULL；TX/RX 行为对称，`PdschProcessor` 调用点已传入真实 `n_filler`。


### 1.7 加扰 (`Scrambler.cpp` + `generate_pn_sequence`) — ✅ 修正后符合规范（2026-09-20 修复）
- 31 阶 Gold 序列、Nc=1600、x1/x2 多项式 ✔（38.211 §5.7）
- **c_init（已修复）**：现按 38.211 §7.3.1.1 计算 `c_init = n_RNTI·2^15 + q·2^14 + n_ID`（单码字 q=0；`n_RNTI`、`n_ID` 由 `SimulationConfig.n_rnti`/`n_id_pdsch` 配置，默认 1/0）。原先的 `c_init = slot_idx+1` 已废除。
- **注意**：c_init 不再随 slot 变化（符合规范），同一用户的每 slot 加扰序列相同——这是标准行为。

### 1.8 调制 (`Modulator.cpp`) — ✅ 正确（优于 Max-Log）
- 3GPP Gray 映射递推 `pam_gray_3gpp` ✔；归一化系数按枚举计算 ✔；位序 b0/b1→实部、b2/b3→虚部 ✔
- 解调用精确 APP（logsumexp），README 写"Max-Log"与实现不符（实现更好，建议改文档）
- **改进**：高阶星座可用 Max-Log + scale 修正提速（当前 O(M·m)）。

### 1.9 层映射 (`LayerMapper.cpp`) — ✅ 正确
- `layer = i % n_layers` 即 38.211 §7.3.3 循环分配 ✔
- **注意**：仅支持单码字；双码字（n_layers>4）未实现，n_layers≤4 时 OK。

### 1.10 预编码 (`Precoder.cpp`) — ✅ Type I CSI 已实现（2026-09-20 修复）
- 单位阵预编码 ✔（端口=层，默认模式）
- **Type I CSI（已实现）**：`build_type1_csi_precoder`（`Precoder.cpp`，声明于 `PhyInterfaces.h`）按 38.214 §5.2.1.2 Type I 单面板码本构造 `W = W1·W2`：DFT 波束 `v_l = exp(j2π(m1·(l mod N1)/N1 + m2·⌊l/N1⌋/N2))`、双极化 `W1 = blkdiag(V, α·V)`、α = e^{jπφ/2}；rank1 直接映射，rank2 相邻波束对 + Hadamard `W2`，列归一化单位功率。2Tx 单极化码本由 N1=N2=1 特例覆盖（α 由 i2 扩展为 4 相位）。配置项：`precoder_mode / csi_n1 / csi_n2 / csi_i1_phi / csi_i11 / csi_i12 / csi_i2`。
- **已知限制**：rank>2 回退单位阵（待扩展 rank3/4 的 W2 结构）；perfect-CSI + Type1 组合未做端口域 H 折叠（LS 路径完整支持：`PdschProcessor::receive` 先按端口估计、再按 W 折叠出层域等效信道）。
- **❌ `deprecode` 仍是占位实现**（实际接收链靠 Equalizer 完成）——建议删除或实现 MMSE 合并以免误用。

### 1.11 DMRS 生成 (`DmrsGenerator.cpp`) — ✅ 修正后符合 Type 1/Type 2 规范（2026-09-20 修复）
- `c_init` 已补全为规范完整式（38.211 §7.4.1.1.1）：`(2^17(N_symb·n_s + l + 1)(2N_ID+1) + 2^17·floor(λ/2) + 2N_ID + n_SCID) mod 2^31`，λ 为 CDM 组。
- **RE 结构（已修复）**：现按规范映射——Type 1：`k = 4n + 2k' + Δ`（每 PRB 6 个 RE，comb-2，Δ=CDM 组 0/1）；Type 2：`k = 6n + k' + Δ`（每 PRB 4 个 RE，Δ=CDM 组 0/1/2）。序列索引 `r(2n+k')` 跨 PRB 连续。
- **多端口（已修复）**：端口 p=1000+l 的 CDM 组 `λ = l mod N_CDM`、频率 OCC `w_f = ±1`（同组第二端口），序列按 CDM 组独立生成，多层不再共享同一导频。
- **已知限制**：时间域 OCC（w_t，双符号 DMRS）未实现（当前仅单符号 DMRS）；DMRS 无 3dB 功率提升（38.214 §5.1.3.3，README 已列入计划）。
- **配套待办**：`ChannelEstimator` 的 LS 噪声估计系数（2/3、1/2）只对 comb-2 相邻间隔成立，Type 2 时应为 1/2 间隔的对应系数（见 3.1 节）；多层场景建议按 CDM 组分离导频后估计。

### 1.12 资源映射 — ✅ 两套实现均已对齐（2026-09-20 修复）
- 实际链路使用 `PdschProcessor.cpp` 内部 `map_pdsch_to_grid`（跳过 DMRS 符号 ✔）。
- **`ResourceMapperImpl`（已修复）**：`map_pdsch/extract_pdsch` 接口新增 `dmrs_symbols` 参数，实现按符号跳过 DMRS，不再覆盖导频；与内部实现行为一致（仍是死代码，建议后续删除或合并）。

### 1.13 OFDM 调制 (`OfdmModulator.cpp`) — ✅ 自洽，⚠️ 与真实 NR 参数有差
- IFFT + sqrt(N) 归一化、DC 置零、半频带搬移、CP 提取正确 ✔
- **⚠️ CP 长度**：仅 symbol 0 用加长 CP；30 kHz（µ=1）时 l=0 **和 l=7** 均应加长（38.211 §5.3.1）。自洽但符号定时结构与标准不同。
- **⚠️ FFT 大小**：`next_pow2(n_sc+1)`（6 PRB→128），非标准 4096·κ 网格；无 guard 建模（可接受）。
- **⚠️ `estimate_fft_size` 启发式**：从信号长度反推 FFT 大小，脆弱；随后 `fix_ofdm_demod_grid` 再补救——建议把 FFT 大小/n_sc 作为显式参数传入 demodulate，删除两处修补函数。

---

## 2. 信道模型

### 2.1 AWGN — ✅ 正确
- 频域快路径：每 RE 加 `N(0, σ²/2)`（实/虚各一），σ²=1/SINR，与单位功率星座一致 ✔

### 2.2 TDL (`TdlChannel.cpp`) — ⚠️ PDP 正确，多普勒是简化模型
- TDL-A/B/C/D/E 时延/功率表与 TR 38.901 Table 7.7.2 一致（已核对 A/B/C）✔；总功率归一化 ✔；LOS tap 用 K 因子混合 ✔
- **⚠️ 多普勒**：每 tap 单一旋转相位，即"每 tap 单径"，非 38.901 的 20+20 簇 Jakes 衰落。窄带下近似 Rayleigh OK，但无法复现多普勒谱形（README 已知此限制，建议标注）。
- **⚠️ 应用方式**：`apply_channel` 频域逐符号乘 H（对角近似），忽略符号间 ISI——时延<CP 时成立；大时延扩展 + 小 CP 需验证 CP 覆盖。
- **❌ 端口一致性**：只生成 `n_layers` 个"TX"通道，而 TX 网格有 `n_tx_ant` 端口；`n_tx_ant > n_layers` 时多余端口信号被静默丢弃（`tx < min(n_tx_ant, n_layers)`）。应统一 tx_ports=n_layers（单位阵下）并断言。
- **⚠️ `get_thermal_noise` 忽略 bandwidth/noise_figure 参数**——接口误导，建议删除或实现 N0=kTB+NF。
- CDL 未实现（框架），与 README 一致。

---

## 3. 接收端 (RX)

### 3.1 信道估计 (`ChannelEstimator.cpp`) — ✅ 端口域 + 间隔感知噪声估计（2026-09-20 修复）
- LS：`h = rx/ref` 逐 pilot ✔；频域线性插值 + 时域保持/插值 ✔
- **噪声估计（已修复）**：按实际 DMRS 子载波间隔计算——边沿两点差 `s = h_next − h_prev`（Var=2σ²，scale=1/2）；内部三点平滑按间隔 (d1,d2) 用最优线性外推 `s = (d2·h_prev + d1·h_next)/(d1+d2) − h_curr`，`scale = 1/(1 + (d2/(d1+d2))² + (d1/(d1+d2))²)`。comb-2（Type 1，d1=d2=2）时 scale=2/3，与原实现一致；Type 2 不等间隔（组内 1、组间 5）时自动得 36/77≈0.468 的正确权重，消除间隔失配偏置。
- **端口域估计（已修复）**：LS 估计现对全部 `n_tx_ant` 端口独立进行（`ch_idx = rx*n_tx + tx`），不再 `tx = l % ref_tx_ant` 复用导频；`PdschProcessor::receive` 将端口域 H 按预编码矩阵 W 折叠为层域等效信道后再均衡，支持 Type I CSI 预编码的多端口场景。
- **改进**：实现 LMMSE（频域相关 + 噪声统计）或时延域滤波（IDFT→去噪窗→DFT），对 comb-2 尤其有效；Doppler cross-corr 估计方法 OK。

### 3.2 均衡 (`Equalizer.cpp`) — ✅ 正确
- LMMSE `(H^H H + σ²I)^{-1}H^H y`，σ² 缩放避免病态 + Cholesky（失败正则化）✔
- 无偏 MMSE 校正 `d = 1 − err_var`，`eff_noise = err_var/d` 传递给解调 ✔（好实践）；n_layers>n_rx_ant 用接收机形式 ✔
- **改进**：SISO 可用标量公式 `h*/(|h|²+σ²)` 提速一个数量级；当前逐 RE 双 solve + 求逆开销大。

### 3.3 解调/解扰/解速率匹配 — ✅ 自洽（偏差同 TX 侧，收发配对抵消）
- LLR 逐位 clip ±20；**⚠️** 高 MCS 高 SNR 时精确 APP LLR 可远超 20，clip 会轻微损失（解码器内部同样 clip，影响有限）。

### 3.4 LDPC 解码 (`LdpcDecoderImpl`) — ✅ 正确，但与文档不符
- 实现是 **Sum-Product（φ 函数法）**而非 README 声称的 "Offset Min-Sum, offset=0.5"。SPA 收敛更快（约省 1/3 迭代）但计算量大 3–5 倍；择一处理：改文档，或真正实现 offset min-sum 提速。
- syndrome 提前终止 ✔；**改进**：layered/shuffled 消息调度可提升收敛。
- **⚠️ 线程安全**：decoder 消息缓冲为成员变量、graph cache 用 static unordered_map，均非线程安全——多线程批量仿真（scripts/run_cpp_perfect_parallel.py）需确认实例隔离。

### 3.5 译码收尾 (`decode_transport_block`) — ✅ 基本正确
- Filler 位 LLR 置 `−CLIP`（对应 bit=0）✔；打孔前 2·Zc 位置 LLR=0 ✔；多 CB 逐 CB CRC + TB CRC ✔
- **⚠️ TDD**：`calculate_pdsch_capacity` 固定用 slot 0 计算 RE 数，S 特殊时隙时 TB 大小与实际 RE 不符。

---

## 4. 仿真控制 (`process_single_snr_point`)

- ⚠️ **统计强度不足**：0 错误 ≥50 块即停 → BLER≈10⁻² 处最多 50 块，置信区间极宽；建议 BLER<0.1 的点 ≥500 块或自适应停止。
- ⚠️ `fast_fading` 路径硬编码 SISO 且以 `tx_grid` 作为 rx（自环解调）——建议删除或标注为仅快速估算。
- ⚠️ `fast_awgn` 中 `rx_res.crc_ok = false;` 为死代码（随后被覆盖）。
- ✅ 每 SINR 点重置种子可复现；点间 `seed_ = local_rng()` 推进，TB/信道独立。

---

## 5. 问题优先级汇总

| 优先级 | 问题 | 位置 | 类型 | 状态 |
|---|---|---|---|---|
| ~~P0~~ | ~~k0 除数 66/50 → 应为 44/42~~ | RateMatcher.cpp:13 | ~~规范错误~~ | **撤销**：与 38.212 Table 5.4.2.1-2 核对后确认原实现正确（66/50） |
| ~~P0~~ | ~~缺 LDPC parity interleaver；多余 BICM 交织~~ | LdpcCodec.cpp / RateMatcher.cpp | ~~规范缺失~~ | **撤销**：38.212 无 parity interleaver；现有交织即 5.4.2.2 规范交织器 |
| P0 | DMRS 非 comb-2、无 CDM、多层端口冲突、稠密导频致噪声估计偏置 | DmrsGenerator.cpp / ChannelEstimator.cpp | 规范错误+性能影响 | **已修复**（Type1/Type2 梳状映射 + CDM/OCC + 完整 cinit；噪声估计间隔感知已跟进） |
| P0→P1 | filler 在速率匹配缓冲中的处理（k0 应基于完整 Ncb，选择时跳过 NULL） | PdschProcessor.cpp / RateMatcher.cpp | 规范偏差 | **已修复**（NullSkipBuffer） |
| P1 | scrambling c_init 应含 nRNTI/q | PdschProcessor.cpp:378 | 规范偏差 | **已修复**（38.211 §7.3.1.1，配置 n_rnti/n_id_pdsch，默认 1/0） |
| P1 | 预编码未实现 Type I CSI 码本 | Precoder.cpp | 规范缺失 | **已实现**（38.214 §5.2.1.2，rank1/2；rank>2 回退单位阵；LS 路径端口域折叠） |
| P1 | ResourceMapperImpl 不跳过 DMRS 符号 | ResourceMapper.cpp | 一致性 | **已修复**（接口加 dmrs_symbols 参数，两套实现对齐） |
| P1 | LS 估计未按端口域/DMRS 间隔计算噪声 | ChannelEstimator.cpp | 规范偏差 | **已修复**（全端口估计 + 间隔感知权重；comb-2 时 scale=2/3，Type2 不等间隔自动加权） |
| P1 | 30kHz 时 l=7 CP 应加长 | OfdmModulator.cpp | 规范偏差 | 待办 |
| P1 | E_r 分配与 38.212 floor/ceil 公式不一致（多 CB 时） | NrTables.cpp | 规范偏差 | 待办 |
| P1 | 信道模型 tx 端口数与 n_layers 不一致 | TdlChannel.cpp | 一致性 | 待办 |
| P2 | README 与实现不符（Max-Log / Offset Min-Sum） | README、Modulator、LdpcCodec | 文档 | 待办 |
| P2 | ResourceMapper/deprecode/fast_fading 死代码清理 | 各处 | 代码卫生 | 待办 |
| P2 | 统计样本量、TDD 容量计算固定 slot0 | PdschProcessor | 仿真方法 | 待办 |
| P2 | decoder/graph cache 线程安全 | LdpcCodec.cpp | 并发 | 待办 |
| P1 | CB 分割缺 R≤0.25 → BG2 分块规则（B∈(3824,∞)、R≤0.25 时单 CB 超 BG2 容量致 `info_bits too long` 崩溃） | NrTables.cpp:269 | 规范缺失 | **已修复**（38.212 §5.2.2：R≤0.25 → C=⌈B/3816⌉，否则 ⌈B/8424⌉） |
| P2 | runner `--los` 参数双消费（`get(i)=="1"\|\|get(i)=="true"` 第二次调用吞掉下一个 flag 值，致 warning + max_blocks 静默回落默认） | pdsch_bler_web.cpp:110 | 参数解析 | **已修复**（先取值再比较） |

### 6.1 修复验证记录（2026-09-20）

**第一批（DMRS 规范化 + filler NULL 处理）**
- 全量基线扫描（13 MCS × 223 SNR 点，AWGN/perfect-CSI/SISO/3PRB）修复后**逐位复现基线**（`IDENTICAL_TO_BASELINE`）：基线配置下 n_filler=0 且使用 perfect CSI，故无回归——符合预期。

**第二批（加扰 c_init / Type I CSI 预编码 / ResourceMapper 跳 DMRS / LS 噪声间隔感知）**
- `test_qpsk_bler`（LS 路径，走新 DMRS + 新噪声估计）：MCS9 瀑布 4→8 dB，与基线一致 ✔
- `test_type1_precoder`（新增，2Tx/rank1/AWGN）：码本输出 W=[1, j]/√2 ✔，链路可解（0 dB 起 BLER=0，fast-AWGN 双端口相干合并 +3 dB 增益为既有行为）✔
- 全量基线扫描重跑：所有 MCS 的 10% BLER 门限与基线一致（如 MCS0 −3 dB、MCS9 3~4 dB、MCS16 9 dB、MCS27 18 dB），仅采样噪声级差异——**符合预期**（c_init 由 slot 相关改为 n_RNTI=1/n_ID=0，序列不同但统计等价）✔
- `results/baseline/` CSV 与曲线已用新结果重新生成，作为后续对照的新基线。

**第三批（A/B 全量对比，2026-09-22 补充）**：以 git HEAD（修复前，worktree 独立构建）与当前版本各跑 perfect/LS 两种 CSI 模式的全量扫描（13 MCS × 223 SNR 点）：
- perfect CSI：平均 ΔSNR@10%BLER = +0.045 dB，全部 |Δ|<0.4 dB → **统计等价，无回归**；
- LS CSI：**13/13 个 MCS 全部改善**，平均 −0.43 dB、最大 −0.83 dB（MCS26）→ comb-2 DMRS + 间隔感知噪声估计带来一致性增益；
- 数据与图：`results/baseline/bler_{prefix,postfix}_{perfect,ls}.csv`、`compare_{perfect,ls}_csi.png`、`FIX_IMPACT_ANALYSIS.md`（详见该分析报告）。

**第四批（50PRB 大 TBS 场景崩溃修复，2026-10-13）**：用户新场景（50 PRB / 4×4 / MCS0 / LS）触发两处崩溃/静默错误：
- `compute_cb_segmentation` 缺 R≤0.25 分支：B=7328、R=0.117 时算出单 CB，选 BG2 但 Zc 卡上限 384（K=3840<7352）→ `info_bits too long for given zc/bgn` 抛出。按 38.212 §5.2.2 补 `R≤0.25 → C=⌈B/3816⌉`（BG2 单 CB ≤3816+24 信息位），否则 `⌈B/8424⌉`（BG1）。
- runner `--los` 双 `get(i)` 调用吞掉 `--max_blocks` → `[WARN] unknown arg 4000` 且 max_blocks 静默回落 2000。
- 验证：同参数复现命令无 WARN、`max_blocks=4000`、4 SNR 点跑通（-2/0 dB 全错 → 2/4 dB 0/50，瀑布在 ~1-2 dB 符合 4×4+LS 损失预期）；`test_ls_mimo`（R>0.25、B=8216 单 CB 路径）回归 `crc_ok=1` 一致；基线场景全走 B≤3824 单 CB 分支（未改动）零回归。



## 6. 建议的验证手段

1. **速率匹配修正后**：与 Sionna `TBEncoder/TBDecoder` 做比特级互测，修正 k0/interleaver/filler 后 TX 输出应逐位一致。
2. **DMRS comb-2 修正后**：平坦 AWGN 下 LS 噪声估计应收敛到真实 σ²（±5%）。
3. **LDPC**：对已知生成矩阵小 Zc（如 2、7）做全枚举 syndrome 校验。
4. **单元测试**：目前 tests/ 仅 CRC、调制两项；建议补 rate-match 往返、CB 分割、LDPC 无信道零误码三件套。


## 7. OFDM 收发一致性专项检查（2026-10-07，应用户问询）

**问题**：发端 `transmit()` 含 IFFT+CP 的 OFDM 调制，收端是否直接处理频域信号、造成收发不对称？

**结论：主执行链路收发一致，不存在不对称；用户观察到的"收端直接频域处理"来自两段被禁用的死代码。**

调用链（唯一入口）：runner → `run_bler_simulation()` → `process_single_snr_point()`（`PdschProcessor.cpp:676`）三分支：

| 分支 | 位置 | 状态 | 数据流 |
|---|---|---|---|
| `fast_awgn` | 705-738 | **死代码**（690 行硬编码 `false`） | 频域 `tx_grid` 直接加噪→跳过 OFDM demod/信道估计/均衡，自带简化解调 |
| `fast_fading` | 740-830 | **死代码**（691 行硬编码 `false`） | 频域抽 RE×h+加噪→自制 MRC LLR，绕过 OFDM demod |
| `else` | 831-851 | **实际执行** | `modulate`(IFFT+CP)→`apply_channel`(时域)→`add_noise`(时域)→`demodulate`(去CP+FFT)→`fix_ofdm_demod_grid`→`receive`(频域)，收发对称 |

补充事实：
- `OfdmModulator.cpp` 同时实现 `modulate()`(86-95, IFFT+CP) 与 `demodulate()`(102-161, 去CP+FFT)，二者以 `√N` 归一化互为 unitary 逆变换，AWGN 下时域加噪与频域加噪统计等价（这也是基线曲线与 Sionna 对齐的原因）。
- AWGN 信道 `apply_channel` 原样返回（`AwgnChannel.cpp:40`），噪声在时域加；TDL/CDL 在信道内逐符号 FFT→乘H→IFFT 处理，收端再 FFT 解回——两条信道路径均与时域边界配对。
- `transmit()` 的 `tx_signal` 被 else 分支真实使用，OFDM 调制不是无效计算。

**遗留风险（对应 §4 P2）**：
1. `estimate_fft_size()`（`OfdmModulator.cpp:26`）从信号长度反推 FFT 大小，再由 `fix_ofdm_demod_grid` 事后补救——建议 `demodulate()` 增加显式 `fft_size`/`n_sc` 参数，删除启发式与修补函数。
2. 两段死代码易误导读代码者（本次问询即源于此），且 fast_fading 的频域抽头模型与 TDL 时域模型不一致——建议删除，或改为 `SimulationConfig` 显式开关并注明与主路径的等价条件。

—— 评审完 ——
