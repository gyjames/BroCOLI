# BroCOLI：基于 Bron–Kerbosch 算法的长读长转录本异构体校准工具

**语言：** [English](README.md) | **简体中文**

[![生命周期：稳定](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![持续维护](https://img.shields.io/badge/Maintained%3F-Yes-red)](https://github.com/gyjames/BroCOLI/graphs/contributors)
[![安装](https://img.shields.io/badge/Install-Github-blue)](#installation)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C.svg)](https://en.cppreference.com/w/cpp/17)
[![许可证：GPL v3](https://img.shields.io/badge/License-GPLv3-ff6b35.svg)](https://www.gnu.org/licenses/gpl-3.0)

## 项目简介

BroCOLI（**Bron-Kerbosch calibrator of Long-read Isoform**）利用高效算法，从长读长 RNA-Seq 数据中识别和定量转录本，支持常规（bulk）、单细胞及空间转录组应用，同时在处理大规模数据集时保持较低的内存占用和较高的运行速度。

该工具包包含**两个可执行程序**：

| 可执行程序 | 功能 |
| :--- | :--- |
| **`prebrocoli`** | 从单细胞及空间转录组文库的原始长读长 FASTQ 数据中提取细胞条形码（barcode）和 UMI |
| **`brocoli`** | 从已比对的长读长数据中识别并定量转录本异构体；一个程序即可处理常规、单细胞和空间转录组数据 |

```mermaid
%%{init: {'themeVariables': { 'fontSize': '20px' }}}%%
graph LR
  classDef input fill:#e3f2fd,stroke:#1565c0,stroke-width:2px,color:black;
  classDef process fill:#fff9c4,stroke:#fbc02d,stroke-width:2px,color:black;
  classDef result fill:#e8f5e9,stroke:#2e7d32,stroke-width:2px,color:black;

  A[原始 fastq.gz]:::input
  C[匹配后的 FASTQ]:::process
  B{排序后的 BAM}:::process
  E[基因与转录本计数]:::result

  A --> |"常规转录组：minimap2"| B
  A --> |"单细胞 / 空间：prebrocoli"| C
  C --> |"minimap2"| B
  B --> |"brocoli"| E

  linkStyle 0 stroke:#ff9800,stroke-width:4px;
  linkStyle 1,2 stroke:#9c27b0,stroke-width:4px;
  linkStyle 3 stroke:#2e7d32,stroke-width:4px;
```

典型的单细胞或空间转录组分析流程为 `prebrocoli` → 比对工具（例如 minimap2）→ `brocoli`。常规转录组数据完成比对后，可直接使用 `brocoli`。

## 目录

- [环境要求](#requirements)
- [安装](#installation)
- [快速开始](#quick-start)
- [支持的文库类型](#supported-chemistries)
- [输出文件](#outputs)
- [仓库结构](#repository-layout)
- [文档](#documentation)
- [参考文献](#reference)
- [联系方式](#contact)

<a id="requirements"></a>
## 🌏 环境要求

| 依赖项 | 说明 |
| :--- | :--- |
| **C++17 编译器** | g++ 17 或更高版本 |
| **htslib** | 用于解析 BAM/CRAM 文件，通过 `pkg-config` 定位 |
| **WFA2-lib** | 用于仿射缺口罚分比对，需要 C++ 绑定（`bindings/cpp/WFAligner.hpp`） |
| **edlib** | 以第三方源码形式引入；将 `edlib.h` 和 `edlib.cpp` 放入 `third_party/edlib/` |
| **zlib** | 用于读取 gzip 压缩的 FASTQ 文件 |

上述依赖均已在 `environment.yml` 中声明，推荐使用 conda 安装和管理。

<a id="installation"></a>
## 🛠️ 安装

```console
git clone https://github.com/gyjames/BroCOLI.git
cd BroCOLI
```

### 🥦 方式 A：使用 conda（推荐）

`setup_conda.sh` 可以一步完成环境的创建（或更新）、单元测试以及两个可执行程序的编译：

```console
bash setup_conda.sh
```

```console
conda activate brocoli
./brocoli -h
./prebrocoli -h
```

### 方式 B：使用 make

如果当前环境已安装所需依赖，可以执行：

```console
make check-env      # 检查 htslib / WFAligner.hpp / edlib，均应显示 "found"
make                # 编译 brocoli 和 prebrocoli
make test           # 仅运行单元测试
```

两种安装方式都会将编译生成的可执行文件放在仓库根目录。

<a id="quick-start"></a>
## 🚀 快速开始

完整的命令行参数说明、参数调优建议和使用示例，请参阅 [BroCOLI 官方文档](https://weiwei4396.github.io/BroCOLI/)。

```console
conda activate brocoli

# 1. 提取 barcode 和 UMI（仅适用于单细胞 / 空间转录组数据）
./prebrocoli -q visium -w whitelist.txt -p 16 -o demux/ reads.fastq.gz

# 2. 使用你偏好的长读长比对工具对拆分后的 reads 进行比对
minimap2 -ax splice -uf ... reference.fa demux/preBroCOLI_matched.fastq | samtools sort -o aligned.bam

# 3. 识别并定量转录本异构体
./brocoli -h
```

<a id="supported-chemistries"></a>
## 🐋 支持的文库类型

使用 `prebrocoli -q <chemistry>` 指定文库类型：

| 文库类型 | 文库结构 | 支持状态 |
| :--- | :--- | :---: |
| `magicseq` | 两段或三段 barcode（BC8 + BC8 [+ BC8]）+ UMI12 | ✅ |
| `visium` | 10x Visium 空间转录组，BC16 + UMI12 | ✅ |
| `10x3v3` | 10x Genomics 3′ v3 | 🚧 |

🚧 表示命令行能够识别该类型，但功能尚未实现。运行时程序会输出说明性错误并终止，而不会在没有提示的情况下返回零个 barcode。

对于 MAGIC-seq，可通过 `--type 2` 或 `--type 3` 选择文库结构。白名单（whitelist）文件第一列的 barcode 长度必须与所选结构对应：

- **16 bp（XY）：** primer → BCX8 → LinkerXY → BCY8 → UMI12 → polyT10。
- **24 bp（XYZ）：** primer → BCX8 → LinkerXY → BCY8 → LinkerYZ → BCZ8 → UMI12 → polyT10。

每个白名单文件只能使用一种 barcode 长度。可以在第二列提供可选的、已知的 12 bp UMI，列之间用制表符或空格分隔。如果一个 barcode 对应多个 UMI，可以在不同行中重复填写该 barcode。空行及以 `#` 开头的注释行会被忽略。

如果首行第一列为 `barcode`、`barcodes`、`CellBarcode`、`CellBarcodes`、`cell_barcode`、`cell_barcodes`、`CB` 或 `BC`（不区分大小写），程序会将该行识别为可选表头。也支持带括号标注的表头，例如 `Barcode (16bp)` 或 `Barcode(16bp)`。其他格式错误的行仍会触发错误，并显示对应行号。

白名单文件也可以仅包含 barcode；此时会保留从 reads 中提取的 UMI，不执行基于参考 UMI 的校正。如果白名单中的 barcode 长度与 `--type` 冲突，程序将终止，并指出有问题的文件及行号。

如果未指定 `--type`，且标准输入（stdin）和标准错误（stderr）都连接到终端，程序会通过 stderr 显示英文提示，要求选择 `2` 或 `3`；直接按 Enter 将选择 `2`。仅重定向标准输出（stdout）不会隐藏这一提示。

如果 stdin 或 stderr 被重定向（包括通过 `nohup` 在后台运行），程序会默认选择 `2`，并在日志中记录该选择。XYZ 文库应显式添加 `--type 3`。输入 reads 必须以文件形式提供，可以作为位置参数传入，也可以使用 `-I`；不支持从标准输入（`-`）读取。已有的 XYZ 分析脚本也必须补充 `--type 3`，即使使用了 `-x/-y/-z` 参数也是如此。

程序会检查所选 `--type` 是否与白名单相匹配，但**不会自动检测或验证 reads 的文库结构**。

MAGIC-seq 的 XY 和 XYZ 文库默认均使用 **corr+top32** 策略：

1. 先按照现有逻辑分别校正各 barcode 片段，再将其组合。如果校正后的组合存在于完整白名单中，则直接接受，不进行 WFA 比对。
2. 否则，使用现有的 9-mer 索引重新召回候选 barcode，根据各片段候选组合的倒排表投票数进行排序，最多选取 32 个候选使用 WFA 评分。如果投票数和 WFA 分数均相同，则按 barcode 字典序确定顺序。

`--topk N` 可修改进入 WFA 评分阶段的候选数量上限；设置为 `0` 表示保留全部候选。`--no-corr` 可以禁用“校正后组合直接命中”的快速路径。

因此，`--no-corr --topk 32` 表示仅使用 top32 候选策略；`--no-corr --topk 0` 表示对全部候选进行 WFA 评分，并采用确定性的并列结果处理规则。后一种比较模式的输出不一定与旧版本逐字节一致，因为旧版本通过哈希表的迭代顺序解决并列情况。

上述选项**仅影响 MAGIC-seq**。运行摘要会报告校正组合命中次数（corr hits）、WFA 搜索次数、截断前的候选数量以及实际 WFA 调用次数。这些计数均以**每次搜索尝试**为单位，而非按输入 read 数量统计。

这种策略在准确率与速度之间作出了权衡：尤其当白名单包含所有可能的组合时，校正后的片段组合可能匹配到错误的 barcode。外部实验报告的速度与准确率数据尚未在此处复现。

UMI 的提取、校正和序列修剪（trimming）仍沿用原有逻辑；不过，如果 barcode 分配发生变化，用于 UMI 校正的参考集合也可能随之变化。

```console
./prebrocoli -q magicseq --type 2 -w whitelist.txt -p 16 -o demux/ reads.fastq.gz

# 三片段（XYZ）文库，使用 24 bp barcode 的白名单
./prebrocoli -q magicseq --type 3 -w whitelist_xyz.txt -p 16 -o demux_xyz/ reads.fastq.gz

# 后台运行：默认采用 XY；如果是 XYZ，请添加 --type 3
nohup ./prebrocoli -q magicseq -w whitelist.txt -p 16 -o demux/ reads.fastq.gz > run.log 2>&1 &
```

`make test` 包含无需 WFA2 即可运行的白名单及文库布局测试。`make test-pre` 会运行完整的 XY/XYZ 提取回归测试，需要 edlib 和 WFA2。

编译完成后，在 Linux 上运行以下命令，可以检查交互式类型选择、后台运行的默认值、白名单冲突处理，以及仅接受文件输入的行为：

```console
python3 tests/test_pre_cli.py ./prebrocoli
```

<a id="outputs"></a>
## 🔥 输出文件

`prebrocoli` 会将数据文件写入 `--outdir` 指定的目录。日志输出到 stdout/stderr，交互式提示输出到 stderr。只有需要保存日志时，才需要使用 shell 重定向。

| 文件 | 内容 |
| :--- | :--- |
| `<prefix>_matched.fastq` | 成功分配 barcode 的 reads；序列 ID 重写为 `barcode_umi#readid`，并附加 `CB:Z:` / `UB:Z:` 标签 |
| `<prefix>_unmatched.fastq` | 未匹配到 barcode 的 reads；仅在设置 `--keep_unmatched` 时生成 |
| `<prefix>_reads_barcodes.tsv` | 每条 read 的 barcode、UMI 和编辑距离 |
| `<prefix>_barcode_counts.tsv` | barcode 对应的 read 数量，按丰度排序 |
| `<prefix>_summary.txt` | 与终端输出相同的运行摘要 |

<a id="repository-layout"></a>
## 👉🏻 仓库结构

```text
BroCOLI/
├── src/                  BroCOLI：转录本异构体发现与定量
├── pre/                  preBroCOLI：barcode 和 UMI 提取
├── tests/                单元测试
├── third_party/edlib/    将 edlib.h 和 edlib.cpp 放在这里
├── environment.yml       conda 依赖配置
├── setup_conda.sh        一键配置环境并编译
├── Makefile
└── CMakeLists.txt
```

`src/` 和 `pre/` 是两套独立的源码目录，分别编译。它们各自都有 `main.cpp` 和 `common.hpp`，因此**不要混用这两个目录中的源码文件**。

<a id="documentation"></a>
## 📘 文档

详细文档和教程请参阅 [BroCOLI 官方文档](https://weiwei4396.github.io/BroCOLI/)。

<a id="reference"></a>
## 🏆 参考文献

1. Li H. Minimap2: pairwise alignment for nucleotide sequences[J]. Bioinformatics, 2018, 34(18): 3094-3100. [Minimap2](https://github.com/lh3/minimap2)
2. Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa. "Fast gap-affine pairwise alignment using the wavefront algorithm.". Bioinformatics, 2020. [WFA2](https://github.com/smarco/WFA2-lib)
3. Martin Šošić, Mile Šikić; Edlib: a C/C++ library for fast, exact sequence alignment using edit distance. Bioinformatics 2017 btw753. doi: 1093/bioinformatics/btw753. [Edlib](https://github.com/Martinsos/edlib)
4. O.Cheng, Flexiplex: a versatile demultiplexer and search tool for omics data, Bioinformatics, Volume 40, Issue 3, 2024. [flexiplex](https://github.com/DavidsonGroup/flexiplex)

<a id="contact"></a>
## ✉️ 联系方式

如果在使用过程中遇到问题，或者有改进建议，欢迎联系 Wei Pan（weipan4396@gmail.com）。如果发现程序错误，也可以在 GitHub 仓库中提交 Issue。
