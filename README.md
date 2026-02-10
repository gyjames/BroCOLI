# BroCOLI : Bron-Kerbosch calibrator of Long-read Isoform
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Maintained?](https://img.shields.io/badge/Maintained%3F-Yes-orange)](https://github.com/gyjames/BroCOLI/graphs/contributors)
[![Install](https://img.shields.io/badge/Install-Github-brightgreen)](#installation)
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

Full documentation can be found [here](https://weiwei4396.github.io/BroCOLI/).

## Reference

1. [C++11 ThreadPool](https://github.com/progschj/ThreadPool)
2. Li H. Minimap2: pairwise alignment for nucleotide sequences[J]. Bioinformatics, 2018, 34(18): 3094-3100.


## Contact
If you come across any issues or have suggestions, please feel free to contact Wei Pan (weipan4396@gmail.com), or open an issue if you find bugs.








