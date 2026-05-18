# 链间二硫键预测工具 v1.4
### 交流，调试请联系 jmwang24@outlook.com 
## 简介
基于 AlphaFold-Multimer 预测两条蛋白链的复合物结构，通过分析跨链半胱氨酸硫原子距离自动识别链间二硫键。支持单条 FASTA、批量块格式、多序列成对模式。

## 系统要求
- Windows 10/11 专业版/企业版 + WSL2 (Ubuntu 22.04/24.04)
- 磁盘 ≥ 20 GB（模型权重 + 临时文件）
- 内存 ≥ 8 GB
- 网络（首次下载模型权重约 3.8 GB，后续预测需访问公共 MSA 服务器）

## 安装步骤

### 1. 安装 WSL2 和 Ubuntu
以管理员身份打开 PowerShell：
```powershell
wsl --install -d Ubuntu-24.04
```
重启电脑，启动 Ubuntu 并设置用户名/密码。

### 2. 在 WSL 中安装 Miniconda 和 ColabFold
进入 WSL 终端（`wsl`），执行：
```bash
# 下载 Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh   # 按提示操作，最后输入 yes 初始化

source ~/.bashrc

# 创建环境并安装 ColabFold
conda create -n colabfold python=3.10 -y
conda activate colabfold
pip install colabfold[alphafold] biopython -i https://pypi.tuna.tsinghua.edu.cn/simple
```

### 3. 在 Windows 中安装 Biopython（用于脚本）
打开 **Anaconda Prompt**（或 PowerShell 并激活你的 Python 环境，如 `scanpy_analysis`）：
```bash
pip install biopython
```

### 4. 模型权重自动下载
首次运行 `colabfold_batch` 时会自动下载 AlphaFold 模型权重（约 3.8 GB）到 `~/.cache/colabfold/`，请保持网络畅通。

## 使用方法

### 输入文件格式

#### 方式一：单条 FASTA（两条序列）
```fasta
>chain_A
GIVEQCCTSICSLYQLENYCN
>chain_B
FVNQHLCGSHLVEALYLVCGERGFFYTPKT
```
运行：
```bash
python predict_interchain.py input.fasta --threshold 2.5 --recycle 3 --outdir results
```

#### 方式二：批量块格式（带标识符）
```text
pair_001
>chain_A
...
>chain_B
...

pair_002
>chain_A
...
>chain_B
...
```
运行：
```bash
python predict_interchain.py batch.txt --batch --threshold 2.5 --recycle 3 --outdir results
```

#### 方式三：多序列成对模式（自动两两配对）
一个 FASTA 文件包含多条序列，按顺序每两条配对为一个蛋白对。标识符自动提取 `pair_xxxx` 或生成序号。
```bash
python predict_interchain.py all_sequences.fasta --paired --threshold 2.5 --recycle 3 --outdir results
```

### 参数说明
| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--threshold` | 2.5 | 二硫键距离阈值（Å），低于该值视为候选 |
| `--recycle` | 3 | AlphaFold 迭代次数（增加提高精度但耗时） |
| `--outdir` | `results` | 输出目录 |
| `--batch` | - | 输入为批量块格式 |
| `--paired` | - | 输入为多序列 FASTA，自动两两配对 |

## 输出结果
- 屏幕打印每个蛋白对的处理进度、所有跨链 Cys 对距离、预测的二硫键列表。
- 输出目录下生成两个文件：
  - `标识符.json`：包含二硫键列表、序列等。
  - `标识符.pdb`：预测的复合物结构（永久保存）。

## 常见问题

**Q1: 提示 `ModuleNotFoundError: No module named 'Bio'`**  
A: 确保运行脚本的 Python 环境安装了 `biopython`（`pip install biopython`）。

**Q2: 提示 `wsl: command not found`**  
A: 在 **Windows PowerShell** 或 **Anaconda Prompt** 中运行脚本，不要在 WSL 内部运行。

**Q3: 预测卡在 MSA 搜索或超时**  
A: 网络问题。可改用单序列模式（修改脚本中的 `--msa-mode mmseqs2_uniref` 为 `--msa-mode single_sequence`），或配置本地 MSA 数据库（需额外 70 GB 磁盘空间）。

**Q4: 如何断点续跑？**  
A: 脚本已自动实现。若输出 JSON 已存在，该蛋白对会被跳过。

**Q5: 如何修改 WSL 用户名？**  
A: 脚本中 `wsl -u jmwang24` 的 `jmwang24` 请改为你的实际 WSL 用户名。

## 版本历史
- **v1.4** (2026-04-02)：增加 `--paired` 模式，自动保存 PDB，断点续跑。

## 依赖清单（完整）
- **Windows 环境**：Python 3.8+，Biopython
- **WSL 环境**：Miniconda，ColabFold（包含 AlphaFold-Multimer v3），MMseqs2（自动安装）
- **自动下载**：AlphaFold 模型权重（~3.8 GB）
- **可选**：UniRef30 数据库（~70 GB，用于本地 MSA）

## 许可证
本工具仅供Tsinghua-IGEM2026研究使用。AlphaFold 模型权重受 DeepMind 许可限制，用户需自行接受条款。
