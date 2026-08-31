# BroCOLI : Bron-Kerbosch calibrator of Long-read Isoform

[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Maintained?](https://img.shields.io/badge/Maintained%3F-Yes-red)](https://github.com/gyjames/BroCOLI/graphs/contributors)
[![Install](https://img.shields.io/badge/Install-Github-blue)](#installation)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C.svg)](https://en.cppreference.com/w/cpp/17)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-ff6b35.svg)](https://www.gnu.org/licenses/gpl-3.0)

## About

BroCOLI (Bron-Kerbosch calibrator of Long-read Isoform) leverages efficient algorithms for transcript identification and quantification from long-read RNA-Seq data, supporting both bulk, single-cell and spatial applications, while maintaining low memory usage and fast performance for large-scale datasets.

The toolkit ships as **two executables**:

| Executable | Role |
| :--- | :--- |
| **`prebrocoli`** | Cell barcode and UMI extraction from raw long-read FASTQ, for single-cell and spatial libraries |
| **`brocoli`** | Isoform detection and quantification from aligned long reads — bulk, single-cell and spatial in one binary |

A typical single-cell/spatial data run is `prebrocoli` → aligner (e.g. minimap2) → `brocoli`. Bulk data goes straight to `brocoli`.

## Table of contents

- [Requirements](#requirements)
- [Installation](#🛠️Installation)
- [Quick start](#quick-start)
- [Supported chemistries](#supported-chemistries)
- [Outputs](#outputs)
- [Repository layout](#repository-layout)
- [Documentation](#documentation)
- [Reference](#reference)
- [Contact](#contact)

## 🌏Requirements

| Dependency | Notes |
| :--- | :--- |
| **C++17 compiler** | g++ 17 or later |
| **htslib** | BAM/CRAM parsing, located via `pkg-config` |
| **WFA2-lib** | gap-affine alignment, C++ bindings (`bindings/cpp/WFAligner.hpp`) |
| **edlib** | vendored — drop `edlib.h` and `edlib.cpp` into `third_party/edlib/` |
| **zlib** | gzipped FASTQ input |

All of the above are declared in `environment.yml`; the recommended route is to let conda provide them.

## 🛠️Installation

```console
git clone https://github.com/gyjames/BroCOLI.git
cd BroCOLI
```

### 🥦Option A — conda (recommended)

`setup_conda.sh` creates (or updates) the environment, runs the unit tests, and builds both executables in one shot:

```console
bash setup_conda.sh
```

```console
conda activate brocoli
./brocoli -h
./prebrocoli -h
```

### Option B — make

Inside an environment that already provides the dependencies:

```console
make check-env      # htslib / WFAligner.hpp / edlib must all report "found"
make                # builds brocoli and prebrocoli
make test           # unit tests only
```

Both routes leave the executables in the repository root.

## 🚀Quick start

Full argument lists, parameter tuning and worked examples live in the [documentation](#documentation).

```console
conda activate brocoli

# 1. barcode + UMI extraction (single-cell / spatial only)
./prebrocoli -q visium -w whitelist.txt -p 16 -o demux/ reads.fastq.gz

# 2. align the demultiplexed reads with your favourite long-read aligner
minimap2 -ax splice -uf ... reference.fa demux/preBroCOLI_matched.fastq | samtools sort -o aligned.bam

# 3. isoform detection and quantification
./brocoli -h
```


## 🐋Supported chemistries

`prebrocoli -q <chemistry>`:

| Chemistry | Library | Status |
| :--- | :--- | :---: |
| `magicseq` | three-segment barcode (BC8 + BC8 + BC8) + UMI | ✅ |
| `visium` | 10x Visium spatial, BC16 + UMI12 | ✅ |
| `10x3v3` | 10x Genomics 3′ v3 | 🚧 |

🚧 = recognised on the command line but not yet implemented; the run stops with an explanatory error rather than silently returning zero barcodes.

## 🔥Outputs

`prebrocoli` writes everything under `--outdir`; nothing goes to std::out, so no shell redirection is needed.

| File | Contents |
| :--- | :--- |
| `<prefix>_matched.fastq` | reads with an assigned barcode, ID rewritten as `barcode_umi#readid` plus `CB:Z:` / `UB:Z:` tags |
| `<prefix>_unmatched.fastq` | reads with no barcode — only with `--keep_unmatched` |
| `<prefix>_reads_barcodes.tsv` | per-read barcode, UMI and edit distances |
| `<prefix>_barcode_counts.tsv` | barcode → read count, sorted by abundance |
| `<prefix>_summary.txt` | the same summary printed to the terminal |

## 👉🏻Repository layout

```
BroCOLI/
├── src/                  BroCOLI — discovery and quantification
├── pre/                  preBroCOLI — barcode and UMI extraction
├── tests/                unit tests
├── third_party/edlib/    edlib.h and edlib.cpp here
├── environment.yml       conda dependencies
├── setup_conda.sh        one-shot environment setup and build
├── Makefile
└── CMakeLists.txt
```

`src/` and `pre/` are two independent source trees and are compiled separately — they each have their own `main.cpp` and their own `common.hpp`, so their files must not be mixed.

## 📘Documentation

Please check out the documentation and tutorials at [BroCOLI Documentation](https://weiwei4396.github.io/BroCOLI/).

## 🏆Reference

1. Li H. Minimap2: pairwise alignment for nucleotide sequences[J]. Bioinformatics, 2018, 34(18): 3094-3100. [Minimap2](https://github.com/lh3/minimap2)
2. Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa. "Fast gap-affine pairwise alignment using the wavefront algorithm.". Bioinformatics, 2020. [WFA2](https://github.com/smarco/WFA2-lib)
3. Martin Šošić, Mile Šikić; Edlib: a C/C++ library for fast, exact sequence alignment using edit distance. Bioinformatics 2017 btw753. doi: 1093/bioinformatics/btw753. [Edlib](https://github.com/Martinsos/edlib)
4. O.Cheng, Flexiplex: a versatile demultiplexer and search tool for omics data, Bioinformatics, Volume 40, Issue 3, 2024. [flexiplex](https://github.com/DavidsonGroup/flexiplex)

## ✉️Contact

If you come across any issues or have suggestions, please feel free to contact Wei Pan (weipan4396@gmail.com), or open an issue if you find bugs.
