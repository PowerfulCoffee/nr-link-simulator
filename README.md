# NR PDSCH 链路级仿真平台

基于 C++17 实现的 3GPP NR PDSCH 端到端链路级仿真平台，遵循 3GPP Release 15/16 协议规范。

## 目录

- [平台架构](#平台架构)
- [功能特性](#功能特性)
- [模块设计](#模块设计)
- [需求补充说明](#需求补充说明)
- [编译与运行](#编译与运行)
- [开发环境配置](#开发环境配置)
- [GitHub 同步方案](#github-同步方案)
- [扩展开发指南](#扩展开发指南)

## 平台架构

### 整体架构分层

```
┌─────────────────────────────────────────────────────────────┐
│                    Application Layer                        │
│  (BLER 仿真、性能对比、结果可视化、参数扫描)                  │
├─────────────────────────────────────────────────────────────┤
│                    PHY Processing Layer                     │
│  ┌──────┐ ┌──────┐ ┌──────┐ ┌──────┐ ┌──────┐ ┌──────┐   │
│  │ CRC  │ │LDPC  │ │Rate  │ │Scram-│ │Modu- │ │Layer │   │
│  │Encode│ │Codec │ │Match │ │bler  │ │lator │ │Map   │   │
│  └──────┘ └──────┘ └──────┘ └──────┘ └──────┘ └──────┘   │
│  ┌──────┐ ┌──────┐ ┌──────┐ ┌──────┐ ┌──────┐ ┌──────┐   │
│  │Pre-  │ │DMRS  │ │Re-   │ │OFDM  │ │ChEst │ │Equal-│   │
│  │coder │ │Gen   │ │Map   │ │Mod   │ │      │ │izer  │   │
│  └──────┘ └──────┘ └──────┘ └──────┘ └──────┘ └──────┘   │
├─────────────────────────────────────────────────────────────┤
│                    Channel Layer                            │
│  ┌────────────┐  ┌────────────┐  ┌────────────┐           │
│  │    AWGN    │  │    TDL     │  │    CDL     │           │
│  └────────────┘  └────────────┘  └────────────┘           │
├─────────────────────────────────────────────────────────────┤
│                    Common Base Layer                        │
│  (类型定义、参数表、工具函数、随机数、日志统计)               │
└─────────────────────────────────────────────────────────────┘
```

### 目录结构

```
nr-link-simulator/
├── include/                     # 头文件目录
│   ├── common/                  # 公共基础模块
│   │   ├── Types.h             # 基础类型定义
│   │   └── NrTables.h          # NR 参数表 (MCS, ZC, DMRS等)
│   ├── phy/                     # PHY层模块
│   │   ├── PhyInterfaces.h     # 所有PHY模块抽象接口
│   │   ├── ModuleFactory.h     # 模块工厂函数声明
│   │   └── PdschProcessor.h    # PDSCH端到端处理器
│   └── channel/                 # 信道模型
│       └── ChannelModels.h     # 信道模型接口与类定义
├── src/                         # 源文件目录
│   ├── common/                  # 公共模块实现
│   ├── phy/                     # PHY层实现
│   └── channel/                 # 信道模型实现
├── tests/                       # 单元测试
├── examples/                    # 示例程序
├── scripts/                     # 构建/运行脚本
├── cmake/                       # CMake模块
├── docs/                        # 文档
└── third_party/                 # 第三方依赖
```

## 功能特性

### 发射端 (TX)

| 模块 | 3GPP协议参考 | 功能说明 |
|------|-------------|---------|
| CRC编码 | TS 38.212 Section 7.1 | CRC16/CRC24A/CRC24B |
| LDPC编码 | TS 38.212 Section 7.2 | BG1/BG2 基图选择、 lifting size Zc |
| 速率匹配 | TS 38.212 Section 7.3 | 码块级联、比特选择、交织、RV版本 |
| 码块级联 | TS 38.212 Section 7.4 | 码块拼接 |
| 加扰 | TS 38.211 Section 7.3.1 | Gold序列生成 |
| 调制 | TS 38.211 Section 7.3.2 | QPSK/16QAM/64QAM/256QAM |
| 层映射 | TS 38.211 Section 7.3.3 | 1-4层映射 |
| 天线端口映射 | TS 38.211 Section 7.3.4 | 预编码 |
| DMRS生成 | TS 38.211 Section 7.4.1 | Type1/Type2, single DMRS |
| 资源映射 | TS 38.211 Section 7.3.5 | RE映射 |
| OFDM调制 | TS 38.211 Section 5.3 | CP-OFDM 调制/解调 |

### 信道模型

| 模型 | 3GPP协议参考 | 说明 |
|------|-------------|------|
| AWGN | - | 加性高斯白噪声信道 |
| TDL | TR 38.901 Section 7.7 | TDL-A/B/C/D/E |
| CDL | TR 38.901 Section 7.7 | CDL-A/B/C/D/E (支持ULA阵列) |

### 接收端 (RX)

| 模块 | 可选算法 | 说明 |
|------|---------|------|
| 信道估计 | LS (默认), MMSE, LMMSE | DMRS导频信道估计+插值 |
| 均衡 | MMSE, ZF | MIMO均衡 |
| 层解映射 | - | 层到码字映射 |
| LLR解调 | 精确LLR, Max-Log | 软判决解映射 |
| 解扰 | - | Gold序列解扰 |
| 解速率匹配 | - | RV解映射、HARQ合并 |
| LDPC解码 | BP, Min-Sum, Layered BP | 支持提前终止 |
| CRC校验 | - | 错误检测 |

### 系统配置

- **子载波间隔**: 30kHz
- **TDD配比**: DDDSUUDDDD 重复 (8DL:2UL per 10 slots)
- **特殊时隙**: 6D:4GP:4U symbols
- **天线配置**: 64 TX ports, 4 RX ports, 4 layers
- **带宽**: 可配置 (默认273 PRB = 100MHz)
- **DMRS**: Type1/Type2, single-symbol DMRS, additionalPos=0

## 模块设计

### 模块化设计原则

每个模块均通过抽象接口解耦，支持独立替换：

```cpp
// 示例：替换信道估计算法
class IChannelEstimator {
public:
    virtual ~IChannelEstimator() = default;
    virtual ComplexCube estimate(...) = 0;
    virtual std::string get_name() const = 0;
};

// 自定义LS估计器
class MyLsEstimator : public IChannelEstimator { ... };

// 使用时注入
processor-&gt;set_channel_estimator(std::make_unique&lt;MyLsEstimator&gt;());
```

### 关键接口一览

| 接口类 | 工厂函数 | 说明 |
|--------|---------|------|
| `ICrcEncoder` | `create_crc_encoder()` | CRC编解码 |
| `ILdpcEncoder` | `create_ldpc_encoder()` | LDPC编码 |
| `ILdpcDecoder` | `create_ldpc_decoder()` | LDPC解码 |
| `IRateMatcher` | `create_rate_matcher()` | 速率匹配 |
| `IScrambler` | `create_scrambler()` | 加扰/解扰 |
| `IModulator` | `create_modulator()` | 调制/解调 |
| `ILayerMapper` | `create_layer_mapper()` | 层映射 |
| `IPrecoder` | `create_precoder()` | 预编码 |
| `IDmrsGenerator` | `create_dmrs_generator()` | DMRS生成 |
| `IResourceMapper` | `create_resource_mapper()` | 资源映射 |
| `IOfdmModulator` | `create_ofdm_modulator()` | OFDM调制 |
| `IChannelEstimator` | `create_ls_channel_estimator()` | 信道估计 |
| `IEqualizer` | `create_mmse_equalizer()` | 均衡 |
| `IChannelModel` | `create_channel(type)` | 信道模型 |

## 需求补充说明

基于资深NR协议仿真经验，对原始需求补充如下设计要点：

### 1. 已补充的关键模块

| 补充项 | 说明 | 必要性 |
|--------|------|--------|
| **OFDM调制/解调** | FFT/IFFT、CP插入/去除 | 链路级仿真必备，原需求遗漏 |
| **TB大小计算** | 按TS 38.214 Section 5.1.3.2计算 | 协议一致性必须 |
| **MCS选择映射** | MCS索引到Qm和R的映射表 | 仿真配置需要 |
| **LDPC基图选择** | BG1/BG2选择逻辑 | 协议规定 |
| **Gold序列生成** | PN序列c_init计算 | 加扰、DMRS必需 |
| **预编码矩阵** | 默认单位矩阵，支持码本 | MIMO传输必需 |

### 2. 建议后续增强

| 增强项 | 优先级 | 说明 |
|--------|--------|------|
| **PTRS支持** | 中 | 相位跟踪参考信号，高MCS/高速场景 |
| **CSI-RS支持** | 低 | 信道状态信息参考信号 |
| **HARQ-ACK反馈** | 中 | 上行控制信息 |
| **完整TDL/CDL功率延迟分布** | 高 | 需填充TR 38.901表7.7.2/7.7.3完整参数 |
| **真正的LDPC BG矩阵** | 高 | 当前为框架，需替换为协议标准校验矩阵 |
| **MMSE/MMSE-IRC均衡** | 中 | 多用户/干扰场景 |
| **LMMSE信道估计** | 中 | 利用信道频域相关性 |
| **多普勒/频偏模型** | 中 | 高速移动场景 |
| **相位噪声** | 低 | 高频段(FR2)需要 |
| **定点仿真支持** | 低 | 硬件验证 |

### 3. 并行化设计

仿真平台使用OpenMP支持SINR点级和TB级并行：
- 不同SINR点间可并行
- 同一SINR下不同TB可并行（需独立RNG）
- OpenMP可通过CMake开关启用

## 编译与运行

### 依赖项

- CMake &gt;= 3.16
- C++17 兼容编译器 (GCC &gt;= 8, Clang &gt;= 9)
- Armadillo &gt;= 9.8 (线性代数库，自带FFT)
- OpenMP (可选，用于并行加速)

**Ubuntu/Debian 安装依赖**:
```bash
sudo apt install cmake g++ libarmadillo-dev libopenmp-dev
```

**macOS 安装依赖**:
```bash
brew install cmake armadillo libomp
```

### 编译

```bash
# 使用构建脚本
cd nr-link-simulator
./scripts/build.sh Release

# 或手动编译
mkdir build &amp;&amp; cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### 运行仿真

```bash
cd build/examples

# 查看帮助
./pdsch_bler_simulation --help

# AWGN信道, MCS10, SINR -2~12dB
./pdsch_bler_simulation --channel AWGN --mcs 10 --sinr-start -2 --sinr-end 12 --output awgn_mcs10.csv

# TDL-A信道, DMRS Type2, MMSE估计器
./pdsch_bler_simulation --channel TDLA --dmrs-type 2 --estimator MMSE --output tdl_mmse.csv

# 快速测试 (减少块数)
./pdsch_bler_simulation --max-blocks 100 --target-errors 10
```

### 命令行参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--mcs &lt;index&gt;` | MCS索引(0-28) | 15 |
| `--sinr-start &lt;dB&gt;` | 起始SINR | -2 |
| `--sinr-end &lt;dB&gt;` | 终止SINR | 18 |
| `--sinr-step &lt;dB&gt;` | SINR步长 | 1 |
| `--channel &lt;type&gt;` | 信道类型(AWGN/TDLA/.../CDLE) | AWGN |
| `--dmrs-type &lt;1/2&gt;` | DMRS类型 | 1 |
| `--max-blocks &lt;N&gt;` | 每SINR最大TB数 | 10000 |
| `--target-errors &lt;N&gt;` | 每SINR目标错误数 | 100 |
| `--estimator &lt;name&gt;` | 信道估计器(LS/MMSE) | LS |
| `--output &lt;file&gt;` | 输出CSV文件 | bler_results.csv |

### 运行测试

```bash
cd build
ctest --output-on-failure
```

## 开发环境配置

### 本地PC (TRAE IDE)

#### VS Code/Trae 配置

在 `.vscode/` 目录下建议创建：

**c_cpp_properties.json**:
```json
{
    "configurations": [
        {
            "name": "Linux",
            "includePath": [
                "${workspaceFolder}/include",
                "/usr/include",
                "/usr/local/include"
            ],
            "defines": [],
            "compilerPath": "/usr/bin/g++",
            "cStandard": "c17",
            "cppStandard": "c++17",
            "intelliSenseMode": "linux-gcc-x64",
            "compileCommands": "${workspaceFolder}/build/compile_commands.json"
        }
    ],
    "version": 4
}
```

**settings.json**:
```json
{
    "cmake.buildDirectory": "${workspaceFolder}/build",
    "cmake.configureSettings": {
        "CMAKE_BUILD_TYPE": "Debug",
        "NR_BUILD_TESTS": "ON"
    },
    "C_Cpp.default.cppStandard": "c++17"
}
```

#### 本地调试

```bash
# Debug模式编译
./scripts/build.sh Debug

# GDB调试
gdb --args ./build/examples/pdsch_bler_simulation --max-blocks 10
```

### 云端环境 (推荐配置)

#### Docker 开发环境

创建 `Dockerfile.dev`:
```dockerfile
FROM ubuntu:22.04

RUN apt update &amp;&amp; apt install -y \
    build-essential cmake gdb git \
    libarmadillo-dev libopenmpi-dev \
    python3 python3-pip \
    &amp;&amp; rm -rf /var/lib/apt/lists/*

RUN pip3 install numpy matplotlib pandas jupyter

WORKDIR /workspace
```

构建并运行：
```bash
docker build -t nr-sim-dev -f Dockerfile.dev .
docker run -it -v $(pwd):/workspace -p 8888:8888 nr-sim-dev
```

#### 云端CI (GitHub Actions)

在 `.github/workflows/` 创建CI配置（可选）用于自动化测试。

## GitHub 同步方案

### 1. 初始化仓库

```bash
cd nr-link-simulator
git init
git add .
git commit -m "Initial commit: NR PDSCH link-level simulator framework"
git branch -M main
```

### 2. 关联GitHub远程仓库

```bash
# 在GitHub上创建空仓库 nr-link-simulator 后:
git remote add origin https://github.com/&lt;your-username&gt;/nr-link-simulator.git
git push -u origin main
```

### 3. 日常工作流

**本地开发**:
```bash
# 创建功能分支
git checkout -b feature/new-channel-estimator

# 开发...修改代码...
# 本地编译测试
./scripts/build.sh Debug
cd build &amp;&amp; ctest

# 提交
git add .
git commit -m "Add LMMSE channel estimator"
git push origin feature/new-channel-estimator
```

**云端同步**:
```bash
# 在云端环境拉取最新代码
git pull origin main

# 云端编译运行
./scripts/build.sh Release
cd build/examples &amp;&amp; ./pdsch_bler_simulation ...
```

### 4. 推荐分支策略

| 分支 | 用途 |
|------|------|
| `main` | 稳定可运行版本 |
| `develop` | 开发集成分支 |
| `feature/*` | 新功能开发 |
| `bugfix/*` | 问题修复 |
| `perf/*` | 算法性能优化 |

### 5. 推荐 .gitattributes (可选)

```
*.h text eol=lf
*.cpp text eol=lf
*.cmake text eol=lf
*.sh text eol=lf
*.csv binary
```

## 扩展开发指南

### 添加新的信道估计算法

```cpp
// 1. 在 ChannelEstimator.cpp 或新文件中创建新类
class LmmseEstimator : public IChannelEstimator {
public:
    ComplexCube estimate(const ResourceGrid&amp; rx_grid, 
                         const ResourceGrid&amp; dmrs_grid,
                         const SimulationConfig&amp; config) override {
        // 实现你的LMMSE估计算法
        ComplexCube h_est(...);
        // ...
        return h_est;
    }
    std::string get_name() const override { return "LMMSE"; }
};

// 2. 在主程序中使用
processor-&gt;set_channel_estimator(std::make_unique&lt;LmmseEstimator&gt;());
```

### 添加新的LDPC解码器

```cpp
// 实现 ILdpcDecoder 接口
class BeliefPropagationDecoder : public ILdpcDecoder {
public:
    std::pair&lt;BitVec, bool&gt; decode(const SoftVec&amp; llr, int bgn, int zc,
                                   int n_iter, bool early_term) override {
        // 实现BP解码算法
    }
};

// 替换默认解码器
// 在 PdschProcessor::init_default_modules() 中修改
```

### 性能曲线绘制 (Python)

仿真输出CSV后，可用Python绘图：

```python
import matplotlib.pyplot as plt
import pandas as pd

# 读取结果
df_ls = pd.read_csv('bler_ls.csv')
df_mmse = pd.read_csv('bler_mmse.csv')

# 绘制BLER曲线
plt.semilogy(df_ls['SINR_dB'], df_ls['BLER'], 'o-', label='LS Estimator')
plt.semilogy(df_mmse['SINR_dB'], df_mmse['BLER'], 's-', label='MMSE Estimator')
plt.xlabel('SINR (dB)')
plt.ylabel('BLER')
plt.grid(True)
plt.legend()
plt.title('NR PDSCH BLER vs SINR (MCS 15, TDL-A 30ns)')
plt.savefig('bler_curve.png', dpi=150)
```

## 协议参考文档

| 文档 | 内容 |
|------|------|
| TS 38.211 | Physical channels and modulation |
| TS 38.212 | Multiplexing and channel coding |
| TS 38.214 | Physical layer procedures for data |
| TR 38.901 | Study on channel model (TDL/CDL) |
| TS 38.101 | User Equipment (UE) radio transmission/reception |
| TS 38.104 | Base Station (BS) radio transmission/reception |

在线参考：
- ShareTechnote 5G NR Handbook: https://www.sharetechnote.com/html/5G/Handbook_5G_Index.html
- 3GPP Specifications: https://www.3gpp.org/ftp/Specs/archive/38_series/

## 性能目标与验证

- **AWGN信道**: BLER曲线应与3GPP约定性能匹配，MCS10 QPSK在~2dB达到BLER=0.1
- **TDL-A 30ns**: 较AWGN有2-4dB损失（依信道估计质量）
- **算法对比**: 相同信道下对比LS vs MMSE信道估计，观察BLER性能差异

## 许可证

本项目仅供学术研究和学习使用。
