# BroCOLI : Bron-Kerbosch calibrator of Long-read Isoform
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Maintained?](https://img.shields.io/badge/Maintained%3F-Yes-red)](https://github.com/gyjames/BroCOLI/graphs/contributors)
[![Install](https://img.shields.io/badge/Install-Github-blue)](#installation)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv2-ff6b35.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)

## About
BroCOLI (Bron-Kerbosch calibrator of Long-read Isoform) leverages efficient algorithms for transcript identification and quantification from long-read RNA-Seq data, supporting both bulk and single-cell applications, while maintaining low memory usage and fast performance for large-scale datasets. 

## Table of contents
- [Requirements](#Requirements)
- [Installation](#Installation)
- [Documentation](#Documentation)
- [Reference](#Reference)
- [Contact](#Contact)


## Requirements
BroCOLI requires a **C++11-compatible compiler** (e.g., g++ 4.8 or later).

## 🛠️Installation
In order to compile the BroCOLI source in this GitHub repository the following steps can be taken:
```console
git clone https://github.com/gyjames/BroCOLI.git
```
```console
cd BroCOLI
sh build.sh
```
Once compiled, two executables ( **BroCOLI_bulk** and **BroCOLI_sc** ) will appear in the **BroCOLI folder**. You can use either the -h (--help) argument or the demo data to verify if the program runs successfully.
```console
./BroCOLI_bulk -h
./BroCOLI_sc -h
```

## 📘Documentation
Please check out the documentation and tutorials at [BroCOLI Documentation](https://weiwei4396.github.io/BroCOLI/).

## Reference

1. Li H. Minimap2: pairwise alignment for nucleotide sequences[J]. Bioinformatics, 2018, 34(18): 3094-3100. [Minimap2](https://github.com/lh3/minimap2)
2. Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa. "Fast gap-affine pairwise alignment using the wavefront algorithm.". Bioinformatics, 2020. [WFA2](https://github.com/smarco/WFA2-lib)
3. Martin Šošić, Mile Šikić; Edlib: a C/C ++ library for fast, exact sequence alignment using edit distance. Bioinformatics 2017 btw753. doi: 1093/bioinformatics/btw753. [Edlib](https://github.com/Martinsos/edlib)
4. O.Cheng, Flexiplex: a versatile demultiplexer and search tool for omics data, Bioinformatics, Volume 40, Issue 3, 2024. [flexiplex](https://github.com/DavidsonGroup/flexiplex)
5. [C++11 ThreadPool](https://github.com/progschj/ThreadPool)


## Contact
If you come across any issues or have suggestions, please feel free to contact Wei Pan (weipan4396@gmail.com), or open an issue if you find bugs.








