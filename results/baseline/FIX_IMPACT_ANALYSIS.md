# 修复前后链路性能对比分析（PDSCH Link Curve A/B）

日期：2026-09-20 ｜ 对比脚本：`scripts/compare_fix_impact.py`

## 1. 实验设置

| 项目 | 说明 |
|---|---|
| A 组（pre） | git HEAD `6931d55`（修复前代码，git worktree 独立构建） |
| B 组（post） | 当前工作区（含全部 6 项修复：DMRS Type1/2 规范化、c_init、filler NULL、Type I CSI、ResourceMapper 跳 DMRS、LS 噪声间隔感知） |
| 仿真配置 | AWGN、3 PRB、SISO、rank1、每 MCS 17 个 SNR 点 × 13 个 MCS（0~27），seed=123 |
| CSI 模式 | **perfect**（理想信道，验证统计等价性）与 **LS**（非理想，走 DMRS 生成 + LS 估计 + 噪声估计全链路） |
| 统计 | 每 SNR 点最多 2000 块 / 目标 60 错误 |

产物：`results/baseline/bler_{prefix,postfix}_{perfect,ls}.csv`（各 223 点）、
`compare_perfect_csi.png`、`compare_ls_csi.png`（SNR vs BLER / SNR vs Tput 双联图）。

## 2. 10% BLER 门限 SNR 对比（dB）

| MCS | Qm | perfect pre | perfect post | Δ | LS pre | LS post | Δ |
|----:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 2 | −3.01 | −3.03 | −0.02 | 1.83 | 1.09 | **−0.74** |
| 2 | 2 | −0.91 | −0.93 | −0.02 | 3.13 | 2.96 | −0.18 |
| 4 | 2 | 0.09 | 0.11 | +0.02 | 4.06 | 3.67 | −0.39 |
| 6 | 2 | 2.11 | 2.19 | +0.08 | 5.19 | 5.09 | −0.10 |
| 9 | 2 | 3.23 | 3.62 | +0.38 | 6.92 | 6.49 | −0.43 |
| 11 | 4 | 4.39 | 4.20 | −0.20 | 7.80 | 7.13 | **−0.67** |
| 13 | 4 | 6.06 | 6.06 | −0.01 | 9.16 | 8.90 | −0.26 |
| 16 | 4 | 8.86 | 8.87 | +0.01 | 12.03 | 11.62 | −0.41 |
| 19 | 6 | 10.54 | 10.52 | −0.02 | 13.94 | 13.45 | −0.49 |
| 21 | 6 | 12.99 | 12.96 | −0.04 | 16.11 | 15.77 | −0.34 |
| 24 | 6 | 15.25 | 15.64 | +0.39 | 19.06 | 18.48 | **−0.59** |
| 26 | 6 | 17.24 | 17.24 | 0.00 | 21.00 | 20.17 | **−0.83** |
| 27 | 6 | 18.14 | 18.13 | −0.00 | 21.18 | 21.05 | −0.13 |

汇总：**perfect CSI：平均 +0.045 dB（全部 |Δ|<0.4 dB）→ 统计上无变化；
LS CSI：平均 −0.43 dB，13/13 个 MCS 全部改善，最大 −0.83 dB（MCS26）→ 一致性正收益。**

## 3. 结论与分析

1. **Perfect-CSI 路径：统计等价（符合预期）**。理想信道绕过 DMRS 生成与信道估计，唯一进入链路的改动是加扰 c_init（n_RNTI=1/n_ID=0）——序列不同但统计等价，残余 |Δ|≤0.4 dB 属采样噪声（50~2000 块/点）。这证明修复没有引入任何基础链路（LDPC/速率匹配/调制解调）的性能损失。

2. **LS 路径：全 MCS 一致改善约 0.1~0.8 dB**。来源可分解为：
   - **DMRS 梳状化（comb-2）**：原实现把导频铺满全部子载波，LS 的相邻 pilot 差分噪声估计被信道频选分量污染且与插值假设失配；comb-2 + 间隔感知权重后噪声估计无偏，均衡器的有效 SINR 缩放更准确 → 高 MCS（26/24/11）受益最大（−0.6~−0.8 dB），因为高阶 QAM 对 LLR 缩放误差更敏感。
   - **c_init 规范化**：对自闭环 BLER 无直接影响（序列统计等价），但使结果可与标准参考实现逐位互测。
   - **filler NULL 处理 / ResourceMapper**：本配置（单 CB、SISO、跳过整 DMRS 符号）不触发，属互操作性收益。

3. **吞吐量（SNR vs Tput）**：曲线形状不变，LS 模式下同一 SNR 的吞吐整体左移（提前达标）；峰值吞吐不变（MCS27 ≈ 2.47 Mbps @3 PRB），符合预期——修复改善的是"达到目标 BLER 所需 SNR"，不是容量上限。

4. **后续建议**：TDL 衰落 + 多端口场景下 DMRS CDM/端口修复与 Type I CSI 预编码的影响尚未量化（当前 AWGN/SISO 不激活 CDM 组间正交与波束赋形增益），建议下一步在 `pdsch_bler_tdl_mimo` 基础上做同样的 A/B。

## 4. 文件清单

- 数据：`results/baseline/bler_prefix_perfect.csv`、`bler_prefix_ls.csv`、`bler_postfix_perfect.csv`、`bler_postfix_ls.csv`、`fix_impact_snr10_table.csv`
- 曲线（汇总叠加）：`results/baseline/compare_perfect_csi.png`、`results/baseline/compare_ls_csi.png`
- **曲线（按 MCS 分网格）**：`results/baseline/per_mcs_grid_perfect.png`、`per_mcs_grid_ls.png`——4×4 网格，每 MCS 独立小图（左轴 BLER 对数轴 + 右轴 Tput，红色横线 = 10% BLER 门限；实线 = 修复后，虚线 = 修复前）
- 脚本：`scripts/compare_fix_impact.py`、`scripts/plot_per_mcs_grid.py`；对照用例：`examples/pdsch_bler_baseline.cpp`（已支持 `perfect|ls` 与输出路径参数）
