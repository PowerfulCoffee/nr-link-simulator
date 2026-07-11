# GitHub上传与本地/云端同步指南

## 概述

本项目代码已在 TRAE Work（云端）环境中完成 Git 初始提交。本文档说明：
1. 如何将代码上传到 GitHub
2. 如何在本地 TRAE IDE 中配置开发环境
3. 本地 TRAE IDE 与 TRAE Work（云端）之间的代码同步方案

---

## 第一步：在 GitHub 创建空仓库

1. 打开浏览器访问：https://github.com/new
2. 填写仓库信息：
   - **Repository name**: `nr-link-simulator`
   - **Description**: `3GPP NR PDSCH Link-Level Simulator (C++)`
   - 选择 **Private** 或 **Public**
   - **不要勾选** "Add a README file"、"Add .gitignore"、"Choose a license"（因为代码已有这些文件）
3. 点击 **Create repository**
4. 创建后页面会显示仓库地址，类似：
   - HTTPS: `https://github.com/你的用户名/nr-link-simulator.git`
   - SSH: `git@github.com:你的用户名/nr-link-simulator.git`

---

## 第二步：获取代码到本地 TRAE IDE

由于 TRAE Work 是云端环境，代码需要先到你的本地电脑。有三种方式：

### 方式 A：通过 TRAE IDE 直接连接 TRAE Work（推荐）

如果你的 TRAE IDE 支持直接打开 TRAE Work 空间：
1. 在 TRAE IDE 中选择 "Open Workspace" 或 "连接云端工作空间"
2. 选择当前的 Work 空间，IDE 会自动同步文件
3. 跳过下方的方式 B/C，直接进入第三步

### 方式 B：下载代码压缩包

1. 在 TRAE Work 的文件浏览器中，找到 `/workspace/nr-link-simulator.tar.gz`
2. 右键下载到本地电脑
3. 解压到你希望的本地开发目录，例如：
   ```
   Windows: D:\Projects\nr-link-simulator
   macOS/Linux: ~/Projects/nr-link-simulator
   ```

### 方式 C：在本地 clone 空仓库

如果你已经在 GitHub 上创建了仓库（第一步），可以先 clone 空仓库，再把下载的代码文件复制进去。

---

## 第三步：在本地 TRAE IDE 配置开发环境

### 3.1 安装依赖

**Windows (推荐使用 WSL2 Ubuntu)**:
```bash
# WSL2 Ubuntu 下：
sudo apt update
sudo apt install -y cmake g++ gdb libarmadillo-dev git

# Python 绘图依赖
pip install matplotlib pandas numpy
```

**macOS**:
```bash
brew install cmake armadillo libomp
pip3 install matplotlib pandas numpy
```

**Linux (Ubuntu/Debian)**:
```bash
sudo apt install -y cmake g++ gdb libarmadillo-dev git
pip3 install matplotlib pandas numpy
```

### 3.2 配置 GitHub 认证（二选一）

#### 方案一：使用 SSH Key（推荐，一次配置长期使用）

```bash
# 1. 生成 SSH Key（如果还没有）
ssh-keygen -t ed25519 -C "你的邮箱@example.com"
# 直接回车使用默认路径，建议设置 passphrase

# 2. 查看公钥
cat ~/.ssh/id_ed25519.pub
# 复制输出的全部内容

# 3. 在 GitHub 网页添加：
# 打开 https://github.com/settings/keys
# 点击 "New SSH key"
# Title 填 "TRAE IDE Local"
# Key 粘贴刚才复制的公钥
# 点击 "Add SSH key"

# 4. 测试连接
ssh -T git@github.com
# 看到 "Hi 用户名! You've successfully authenticated" 表示成功
```

#### 方案二：使用 Personal Access Token (HTTPS)

如果不想用 SSH：
1. 打开 https://github.com/settings/tokens
2. 生成一个新的 Token（勾选 `repo` 权限）
3. 保存好 Token（只显示一次）
4. 第一次 push 时 Git 会提示输入用户名和密码，密码处粘贴 Token

### 3.3 首次推送代码到 GitHub

进入本地项目目录，执行：

```bash
cd /path/to/nr-link-simulator  # 进入解压后的代码目录

# 检查是否在git仓库中
git status
```

如果 `.git` 目录不存在（比如解压 tar.gz 后）：
```bash
git init
git add .
git commit -m "Initial commit: NR PDSCH link-level simulator framework"
git branch -M main
```

添加远程仓库并推送（**把"你的用户名"替换为你实际的GitHub用户名**）：
```bash
# 使用 SSH（如果配置了SSH key，推荐）：
git remote add origin git@github.com:你的用户名/nr-link-simulator.git

# 或者使用 HTTPS：
# git remote add origin https://github.com/你的用户名/nr-link-simulator.git

# 首次推送
git push -u origin main
```

推送成功后，刷新 GitHub 仓库页面，就能看到所有代码文件。

### 3.4 验证本地编译

```bash
# 编译
./scripts/build.sh Debug

# 运行单元测试
cd build
ctest --output-on-failure

# 运行一个简单仿真
cd examples
./pdsch_bler_simulation --max-blocks 10
```

---

## 第四步：本地 TRAE IDE 与 TRAE Work（云端）同步方案

GitHub 作为"中央仓库"，本地和云端都通过 push/pull 同步：

```
  本地 TRAE IDE          GitHub.com         TRAE Work (云端)
      │                    │                    │
      ├── push ──────────&gt; │                    │
      │                    │ &lt;──── pull ────────┤
      │                    │                    │
      ├────── pull &lt;───────┤                    │
      │                    │ ────────── push ──&gt;│
      │                    │                    │
```

### 4.1 日常开发工作流

#### 场景 1：在本地开发，同步到云端和 GitHub

**本地电脑（TRAE IDE）**:
```bash
# 开始工作前先拉取最新
git pull origin main

# 开发、修改代码、编译测试...

# 测试通过后提交
git add -A
git commit -m "描述你的修改，比如: Add TDL-A channel model taps"
git push origin main
```

**TRAE Work（云端）获取最新代码**:
```bash
cd /workspace/nr-link-simulator
git pull origin main
```

#### 场景 2：在云端调试/跑仿真后，同步回本地

如果在云端做了修改需要保留：
```bash
# 云端（TRAE Work）:
cd /workspace/nr-link-simulator
git add -A
git commit -m "Fix channel estimation interpolation bug"
git push origin main
```
注意：云端 push 需要先配置 GitHub 认证（见 4.3 节）。如果无法配置，建议只在本地做提交，云端只做编译运行。

**本地获取云端提交**:
```bash
git pull origin main
```

### 4.2 推荐的最佳实践

1. **以本地开发为主，GitHub 为真相源**
   - 主要编码、单步调试在本地 TRAE IDE 中进行（调试更方便，不受网络影响）
   - 每完成一个功能模块就 commit + push
   - 云端 TRAE Work 用于需要服务器资源的场景（大规模BLER仿真跑数、长时间运算）

2. **功能分支工作流（推荐用于新功能开发）**
   ```bash
   # 本地创建功能分支
   git checkout -b feature/ldpc-min-sum-decoder
   
   # 在分支上开发、测试...
   git add -A
   git commit -m "Implement min-sum LDPC decoder"
   git push origin feature/ldpc-min-sum-decoder
   
   # 功能验证完成后合并到 main
   git checkout main
   git pull origin main
   git merge feature/ldpc-min-sum-decoder
   git push origin main
   ```

3. **避免冲突**
   - 不要同时在本地和云端修改同一文件
   - 每次开始工作前先 `git pull` 拉取最新代码
   - 提交前先 `git pull --rebase` 保持提交历史线性

4. **云端运行大规模仿真的典型流程**
   - 本地调试通过小 case（`--max-blocks 50`）确认代码正确
   - commit & push 到 GitHub
   - 在 TRAE Work 中 `git pull` 拉取最新代码
   - 在云端编译 Release 版本（`./scripts/build.sh Release`）
   - 跑大规模仿真（`--max-blocks 10000 --target-errors 200`）
   - 仿真结果 CSV 可以下载到本地用 Python 绘图

### 4.3 在 TRAE Work 配置 GitHub 认证（支持云端 push）

如果希望能直接在云端提交并 push 到 GitHub，需要配置 SSH key：

```bash
# 在 TRAE Work 终端执行
ssh-keygen -t ed25519 -C "trae-work" -f ~/.ssh/id_ed25519 -N ""
cat ~/.ssh/id_ed25519.pub
```

复制输出的公钥内容，打开 https://github.com/settings/keys 添加新 SSH key（Title 填 "TRAE Work Cloud"），之后就可以在云端正常执行 `git push` 了。

---

## 第五步：本地 TRAE IDE 调试配置

### .vscode/launch.json（GDB 调试配置）

在项目下创建 `.vscode/launch.json`：
```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "Debug BLER Simulation",
            "type": "cppdbg",
            "request": "launch",
            "program": "${workspaceFolder}/build/examples/pdsch_bler_simulation",
            "args": ["--max-blocks", "10", "--sinr-start", "0", "--sinr-end", "4"],
            "stopAtEntry": false,
            "cwd": "${workspaceFolder}/build/examples",
            "environment": [],
            "externalConsole": false,
            "MIMode": "gdb",
            "setupCommands": [
                {
                    "description": "Enable pretty-printing for gdb",
                    "text": "-enable-pretty-printing",
                    "ignoreFailures": true
                }
            ]
        },
        {
            "name": "Debug CRC Test",
            "type": "cppdbg",
            "request": "launch",
            "program": "${workspaceFolder}/build/tests/test_crc",
            "args": [],
            "stopAtEntry": false,
            "cwd": "${workspaceFolder}/build/tests",
            "MIMode": "gdb"
        }
    ]
}
```

### .vscode/c_cpp_properties.json（IntelliSense 配置）
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

---

## 快速命令参考

| 操作 | 命令 |
|------|------|
| 查看当前状态 | `git status` |
| 拉取最新代码 | `git pull origin main` |
| 添加所有修改 | `git add -A` |
| 提交 | `git commit -m "说明"` |
| 推送到 GitHub | `git push origin main` |
| 创建功能分支 | `git checkout -b feature/xxx` |
| 切换到主分支 | `git checkout main` |
| 查看提交历史 | `git log --oneline -10` |
| 编译 Debug 版本 | `./scripts/build.sh Debug` |
| 运行单元测试 | `cd build &amp;&amp; ctest --output-on-failure` |
| 查看仿真参数帮助 | `./build/examples/pdsch_bler_simulation --help` |

---

## 常见问题

1. **Push 被拒绝 (rejected)**：先执行 `git pull --rebase origin main`，解决冲突后再 push
2. **合并冲突**：打开带有冲突标记的文件，手动解决 `&lt;&lt;&lt;&lt;&lt;&lt;&lt;` / `=======` / `&gt;&gt;&gt;&gt;&gt;&gt;&gt;` 标记的部分，然后 `git add 文件名` 再 `git rebase --continue`
3. **编译找不到 armadillo**：Ubuntu/Debian 执行 `sudo apt install libarmadillo-dev`；macOS 执行 `brew install armadillo`
4. **云端无法 push**：按 4.3 节配置 SSH key，或者把代码文件下载到本地后从本地 push
