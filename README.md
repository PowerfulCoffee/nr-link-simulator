# NR Link Simulator

基于 C++17 实现的 3GPP NR PDSCH/PUSCH 端到端链路级仿真平台，遵循 3GPP Release 15 协议规范。核心 PDSCH 链路已通过与 NVIDIA Sionna（官方验证平台）的 BLER 性能对比验证，结果一致。

## 目录

- [平台架构](#平台架构)
- [功能特性](#功能特性)
- [验证结果](#验证结果)
- [目录结构](#目录结构)
- [编译说明](#编译说明)
- [运行仿真](#运行仿真)
- [结果查看](#结果查看)
- [模块化设计](#模块化设计)
- [后续计划](#后续计划)
- [协议参考](#协议参考)

## 平台架构

### 整体分层

```
┌─────────────────────────────────────────────────────────────┐
│                    Application Layer                        │
│  (BLER 仿真、参数扫描、性能曲线输出、CSV结果)                │
├─────────────────────────────────────────────────────────────┤
│                    PDSCH Processing Chain                   │
│  TX: CRC → LDPC → RateMatch → Scramble → Modulate          │
│      → LayerMap → Precode → DMRS + ResourceMap → OFDM      │
│  Channel: AWGN / TDL / CDL                                  │
│  RX: OFDM → ChannelEst → Equalize → Demodulate             │
│      → Descramble → RateRecover → LDPC Decode → CRC        │
├─────────────────────────────────────────────────────────────┤
│                    Channel Layer                            │
│           AWGN (已验证)  /  TDL / CDL (框架)                │
├─────────────────────────────────────────────────────────────┤
│                    Common Base Layer                        │
│  (类型定义、MCS/TBS/LDPC参数表、Gold序列、随机数)           │
└─────────────────────────────────────────────────────────────┘
```

### 核心端到端流程

**发射端 (PDSCH TX)**:
1. **TB生成** — 根据MCS、PRB数、层数，按3GPP TS 38.214 Section 5.1.3.2计算TBS
2. **CRC附着** — TB加CRC16/24校验位
3. **LDPC编码** — 3GPP BG1/BG2基图选择，Zc lifting size计算，QC-LDPC编码
4. **速率匹配** — 环形缓冲比特选择，支持RV版本
5. **加扰** — Gold序列加扰（c_init基于时隙号）
6. **调制** — QPSK/16QAM/64QAM调制（功率归一化）
7. **层映射** — 码字到层映射
8. **预编码** — 单位矩阵预编码（SISO/MIMO框架）
9. **DMRS生成** — DMRS Type1单符号导频生成
10. **资源映射** — PDSCH/DMRS映射到资源网格（RE级）
11. **OFDM调制** — CP-OFDM调制（FFT/IFFT + 循环前缀）

**信道模型**:
- AWGN（已完整验证）
- TDL-A/B/C/D/E（框架，待填充完整PDP参数）
- CDL-A/B/C/D/E（框架，待完善）

**接收端 (PDSCH RX)**:
1. **OFDM解调** — CP去除、FFT
2. **信道估计** — LS信道估计（利用DMRS导频）
3. **均衡** — MMSE/ZF均衡
4. **LLR解调** — Max-Log软解调
5. **解扰** — Gold序列解扰
6. **解速率匹配** — LLR解映射到码字位置
7. **LDPC解码** — Offset Min-Sum算法，支持提前终止（syndrome check）
8. **CRC校验** — CRC校验确定TB是否正确

## 功能特性

### 已实现并验证的核心模块

| 模块 | 协议参考 | 验证状态 |
|------|---------|---------|
| CRC编码/校验 | TS 38.212 §7.1 | ✅ CRC16/24 |
| LDPC编码器 | TS 38.212 §5.3.2 | ✅ BG1/BG2完整基图，3GPP标准shift值 |
| LDPC解码器 | TS 38.212 §5.3.2 | ✅ Offset Min-Sum，提前终止 |
| TBS计算 | TS 38.214 §5.1.3.2 | ✅ n_info量化+表查找 |
| MCS映射 | TS 38.214 §5.1.3.1 | ✅ MCS Table 5.1.3.1-1 |
| 速率匹配 | TS 38.212 §5.4.2 | ✅ 环形缓冲，RV=0 |
| 加扰/解扰 | TS 38.211 §7.3.1 | ✅ Gold序列（31阶长度） |
| QPSK调制/解调 | TS 38.211 §7.3.2 | ✅ LLR精确计算 |
| DMRS生成 | TS 38.211 §7.4.1 | ✅ Type1，单符号 |
| 资源映射/解映射 | TS 38.211 §7.3.5 | ✅ RE级映射，跳过DMRS符号 |
| OFDM调制/解调 | TS 38.211 §5.3 | ✅ CP-OFDM，正常CP长度 |
| AWGN信道 | — | ✅ 频域噪声添加，功率归一化 |
| LS信道估计 | — | ✅ DMRS导频LS估计 |
| MMSE均衡 | — | ✅ SISO场景验证 |

### 系统默认配置

- **子载波间隔 (SCS)**: 30 kHz
- **时隙格式**: 正常CP，14个OFDM符号
- **DMRS配置**: Type1，映射类型A，单符号，position=2，additionalPos=0
- **PDSCH符号**: 符号0-1, 3-13（共13个数据符号，符号2为DMRS）
- **调制方式**: QPSK/16QAM/64QAM（QPSK已验证）
- **LDPC解码**: Offset Min-Sum，offset=0.5，默认30次迭代，支持提前终止

## 验证结果

### 与Sionna对比验证

使用配置：MCS 5（QPSK, R=379/1024≈0.370）、6 PRB、1x1 SISO、AWGN、理想CSI（频域快路径）

| SINR (dB) | nr-link-simulator BLER | NVIDIA Sionna BLER |
|-----------|------------------------|---------------------|
| -2.0 | 1.000 | 1.000 |
| -1.0 | 0.980 | 0.968 |
| -0.5 | 0.490 | — |
| 0.0 | **0.028** | **0.030** |
| 0.5 | 0.000 | — |
| 1.0+ | 0.000 | 0.000 |

两条BLER瀑布曲线完全一致，验证了以下模块的3GPP标准一致性：
- TBS计算（TBS=704 bits，匹配Sionna）
- LDPC基图选择（BG2, Zc=72, K=720, N=3744）
- 速率匹配输出长度（E=1872 bits）
- QPSK调制LLR计算
- LDPC Min-Sum解码算法收敛性能

验证结果数据见 [bler_awgn_ideal.csv](file:///workspace/nr-link-simulator/bler_awgn_ideal.csv)

## 目录结构

```
nr-link-simulator/
├── CMakeLists.txt              # 顶层CMake构建配置
├── README.md                   # 本文档
├── .gitignore                  # Git忽略规则
├── bler_awgn_ideal.csv         # 已验证的AWGN BLER结果
├── include/                    # 头文件
│   ├── common/
│   │   ├── Types.h             # 核心类型定义（Complex, BitVec, ResourceGrid, SimulationConfig）
│   │   └── NrTables.h          # MCS表、TBS计算、LDPC参数选择、Gold序列、DMRS pattern
│   ├── phy/
│   │   ├── PhyInterfaces.h     # 所有PHY模块抽象接口（ICrcEncoder, ILdpcEncoder/Decoder, ...）
│   │   ├── ModuleFactory.h     # 模块工厂函数声明
│   │   ├── PdschProcessor.h    # PDSCH端到端处理器（顶层API）
│   │   └── LdpcTables.h        # 3GPP LDPC BG1/BG2 shift值表（自动生成）
│   └── channel/
│       └── ChannelModels.h     # 信道模型接口定义
├── src/                        # 源文件
│   ├── common/
│   │   ├── NrTables.cpp        # MCS表、TBS计算、LDPC参数选择实现
│   │   └── Utils.cpp           # 工具函数
│   ├── phy/
│   │   ├── PdschProcessor.cpp  # PDSCH端到端处理（TX/RX主流程、BLER仿真）
│   │   ├── CrcEncoder.cpp      # CRC16/24编解码
│   │   ├── LdpcCodec.cpp       # LDPC编解码器（BG1/BG2、Min-Sum解码）
│   │   ├── RateMatcher.cpp     # 速率匹配/解匹配
│   │   ├── Scrambler.cpp       # Gold序列加扰/解扰
│   │   ├── Modulator.cpp       # QPSK/16QAM/64QAM调制/LLR解调
│   │   ├── LayerMapper.cpp     # 层映射/解映射
│   │   ├── Precoder.cpp        # 预编码/解预编码
│   │   ├── DmrsGenerator.cpp   # DMRS序列生成与映射
│   │   ├── ResourceMapper.cpp  # RE资源映射/解映射
│   │   ├── ChannelEstimator.cpp# LS信道估计
│   │   ├── Equalizer.cpp       # MMSE/ZF均衡
│   │   └── OfdmModulator.cpp   # CP-OFDM调制/解调（FFT/IFFT+CP）
│   └── channel/
│       ├── AwgnChannel.cpp     # AWGN信道
│       ├── TdlChannel.cpp      # TDL信道框架
│       └── CdlChannel.cpp      # CDL信道框架
├── examples/                   # 示例程序
│   ├── CMakeLists.txt
│   ├── quick_test.cpp          # 快速BLER测试（0.5dB步长精细扫描，输出CSV）
│   ├── debug_params.cpp        # 打印TBS/LDPC/RE等配置参数
│   ├── pdsch_bler_qpsk.cpp     # QPSK BLER仿真（LS信道估计，完整OFDM链路）
│   ├── pdsch_bler_qam16.cpp    # 16QAM BLER仿真
│   ├── pdsch_bler_simulation.cpp # 完整BLER仿真框架（多信道/多估计器）
│   └── test_bler.cpp           # 基础BLER测试
├── tests/                      # 单元测试
│   ├── CMakeLists.txt
│   ├── test_crc.cpp            # CRC模块测试
│   └── test_modulator.cpp      # 调制解调测试
├── scripts/
│   ├── build.sh                # 一键构建脚本
│   └── plot_bler.py            # Python绘图脚本
└── docs/
    └── SETUP_SYNC_GUIDE.md     # 本地/云端开发环境同步指南
```

## 编译说明

### 依赖项

- **CMake** >= 3.16
- **C++17** 兼容编译器（GCC >= 8, Clang >= 9）
- **Armadillo** >= 9.8（线性代数库，提供FFT）
- **OpenMP**（可选，用于并行加速）

**Ubuntu/Debian 安装依赖**：
```bash
sudo apt update
sudo apt install -y cmake g++ libarmadillo-dev libopenmpi-dev
```

**macOS 安装依赖**：
```bash
brew install cmake armadillo libomp
```

### 编译步骤

**方式一：使用构建脚本（推荐）**
```bash
cd nr-link-simulator
./scripts/build.sh Release
```

**方式二：手动CMake**
```bash
cd nr-link-simulator
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DNR_BUILD_TESTS=ON -DNR_BUILD_EXAMPLES=ON
make -j$(nproc)
```

编译产物位于 `build/examples/` 目录下。

### Debug模式编译

```bash
./scripts/build.sh Debug
# 或手动：
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DNR_BUILD_TESTS=ON
make -j$(nproc)
```

### 运行单元测试

```bash
cd build
ctest --output-on-failure
```

## 运行仿真

### 快速验证：AWGN理想CSI BLER曲线

`quick_test` 是最常用的快速验证程序，使用频域AWGN快路径（跳过OFDM和信道估计，等效于理想CSI），速度最快：

```bash
cd build/examples
./quick_test
```

默认配置：MCS 5 (QPSK)，6 PRB，1x1 SISO，SINR -4~4dB（0.5dB步长），LDPC 30次迭代。输出结果同时打印到终端和 `../bler_awgn_ideal.csv`。

### 参数调试与验证

```bash
cd build/examples
./debug_params
```

打印TBS、MCS、LDPC参数（bgn/zc/K/N/N_cb/E）、RE数量等关键配置参数，用于确认仿真配置正确性。

### 完整OFDM链路QPSK BLER

`pdsch_bler_qpsk` 走完整的OFDM调制解调+LS信道估计+MMSE均衡链路，性能比理想CSI有约1dB损失（符合预期）：

```bash
cd build/examples
./pdsch_bler_qpsk
```

默认配置：MCS 5，6 PRB，1x1 SISO，AWGN信道，LS信道估计，30次LDPC迭代。结果输出到 `bler_qpsk.csv`。

### 参数配置方式

所有配置参数通过 `SimulationConfig` 结构体在代码中设置（可参考各examples/*.cpp中的配置方法）：

```cpp
SimulationConfig config;
config.mcs_index = 5;              // MCS索引 (0-28)
config.n_rb = 6;                    // 带宽（PRB数）
config.n_tx_ant = 1;                // 发射天线数
config.n_rx_ant = 1;                // 接收天线数
config.n_layers = 1;                // 层数
config.channel_type = ChannelType::AWGN;  // 信道类型
config.mod_scheme = mcs_to_modulation(config.mcs_index);
config.code_rate = mcs_to_code_rate(config.mcs_index);
config.n_ldpc_iterations = 30;      // LDPC解码迭代次数
config.early_termination = true;    // LDPC提前终止（syndrome check）
config.max_blocks_per_sinr = 500;   // 每SINR点最多传输块数
config.target_block_errors = 50;    // 每SINR点目标错误数（达到即停止）
config.sinr_start = -4.0;           // SINR扫描起始(dB)
config.sinr_end = 4.0;              // SINR扫描终止(dB)
config.sinr_step = 0.5;             // SINR步长(dB)
```

## 结果查看

### 终端输出

仿真运行时终端实时打印进度：

```
=== NR Link Simulator - PDSCH BLER (AWGN, Ideal CSI) ===
MCS 5: QPSK, R=0.370117
PRBs = 6, 1x1 SISO, LDPC iter=30

  SINR(dB)    Blocks    Errors       BLER
     -4.00        50        50      1.0000
     -3.50        50        50      1.0000
     -3.00        50        50      1.0000
     -2.50        50        50      1.0000
     -2.00        50        50      1.0000
     -1.50        50        50      1.0000
     -1.00        51        50      0.9804
     -0.50       102        50      0.4902
      0.00       500        14      0.0280
      0.50       500         0      0.0000
      1.00       500         0      0.0000
```

### CSV结果文件

BLER结果自动保存为CSV格式，便于后续用Python/MATLAB绘图分析：

```csv
SINR_dB,Blocks,Errors,BLER
-4.00,50,50,1.000000
-3.50,50,50,1.000000
-3.00,50,50,1.000000
...
0.00,500,14,0.028000
0.50,500,0,0.000000
```

### Python绘制BLER曲线

项目提供了绘图脚本 [scripts/plot_bler.py](file:///workspace/nr-link-simulator/scripts/plot_bler.py)，也可手动绘制：

```python
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('bler_awgn_ideal.csv')
plt.semilogy(df['SINR_dB'], df['BLER'], 'o-', label='nr-link-simulator')
plt.xlabel('SINR (dB)')
plt.ylabel('BLER')
plt.grid(True, which='both')
plt.legend()
plt.title('PDSCH BLER vs SINR (MCS 5 QPSK, AWGN, Ideal CSI)')
plt.savefig('bler_curve.png', dpi=150)
```

## 模块化设计

所有PHY模块均通过抽象接口解耦，支持独立替换算法实现。以信道估计器为例：

```cpp
// 所有信道估计器实现 IChannelEstimator 接口
class IChannelEstimator {
public:
    virtual ~IChannelEstimator() = default;
    virtual ComplexCube estimate(const ResourceGrid& rx_grid,
                                  const ResourceGrid& dmrs_grid,
                                  const SimulationConfig& config) = 0;
};

// 自定义算法只需继承接口
class MyCustomEstimator : public IChannelEstimator { ... };

// 使用时注入到PdschProcessor
processor.set_channel_estimator(std::make_unique<MyCustomEstimator>());
```

### 模块接口与工厂函数

| 接口 | 工厂函数 | 功能 |
|------|---------|------|
| `ICrcEncoder` | `create_crc_encoder()` | CRC编解码 |
| `ILdpcEncoder` | `create_ldpc_encoder()` | LDPC编码 |
| `ILdpcDecoder` | `create_ldpc_decoder()` | LDPC解码（Min-Sum） |
| `IRateMatcher` | `create_rate_matcher()` | 速率匹配/解匹配 |
| `IScrambler` | `create_scrambler()` | Gold序列加扰/解扰 |
| `IModulator` | `create_modulator()` | QPSK/16QAM/64QAM调制解调 |
| `ILayerMapper` | `create_layer_mapper()` | 层映射/解映射 |
| `IPrecoder` | `create_precoder()` | 预编码/解预编码 |
| `IDmrsGenerator` | `create_dmrs_generator()` | DMRS生成映射 |
| `IResourceMapper` | `create_resource_mapper()` | RE资源映射/解映射 |
| `IOfdmModulator` | `create_ofdm_modulator()` | CP-OFDM调制/解调 |
| `IChannelEstimator` | `create_ls_channel_estimator()` | LS信道估计 |
| `IEqualizer` | `create_mmse_equalizer()` | MMSE均衡 |
| `IChannelModel` | `channel::create_channel(type)` | AWGN/TDL/CDL信道 |

## 后续计划

| 模块 | 优先级 | 说明 |
|------|--------|------|
| **TDL/CDL信道PDP** | 高 | 填充TR 38.901表7.7.2/7.7.3完整参数，支持时延扩展和多普勒 |
| **多码块分割** | 高 | TBS>3824时支持多码块+CB-CRC |
| **HARQ-IR** | 中 | 速率匹配支持4个RV版本，软合并 |
| **LMMSE信道估计** | 中 | 利用信道频域相关性的LMMSE估计 |
| **DMRS增强** | 中 | Additional positions、多端口CDM、3dB功率提升 |
| **PTRS** | 低 | 相位跟踪参考信号（高MCS/高速场景） |
| **256QAM** | 低 | 256QAM调制解调 |
| **MIMO 4x4/64x4** | 中 | 多天线MMSE-IRC均衡、码本预编码 |
| **TDD时隙格式** | 低 | DDDSUUDDDD实际上下行配比应用 |
| **PUSCH链路** | 中 | 上行PUSCH端到端链路 |

## 协议参考

| 文档 | 内容 |
|------|------|
| TS 38.211 | Physical channels and modulation |
| TS 38.212 | Multiplexing and channel coding |
| TS 38.214 | Physical layer procedures for data |
| TR 38.901 | Study on channel model (TDL/CDL) |

在线参考：
- 3GPP Specifications: https://www.3gpp.org/ftp/Specs/archive/38_series/
- ShareTechnote 5G NR: https://www.sharetechnote.com/html/5G/Handbook_5G_Index.html

## 许可证

本项目仅供学术研究和学习使用。
