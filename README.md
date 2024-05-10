[![Python 3.8](https://img.shields.io/badge/python-3.8%20%7C%203.9%20%7C%203.10-blue)](https://www.python.org/downloads/release/python-3100/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<img src="/permGWAS_logo.png" data-canonical-src="/permGWAS_logo.png" height="80" />  

## permGWAS2

This is an improved version of permGWAS. The original version can be found at [permGWAS Version1](https://github.com/grimmlab/permGWAS/releases/tag/permGWAS)

permGWAS2 is an open source software tool written in python to efficiently perform genome-wide association studies (GWAS)
with permutation-based thresholds. It uses a batch-wise Linear Mixed Model to compute several univariate tests simultaneously. 
permGWAS2 provides support for multiple CPUs as well as for GPUs. 

In contrast to the original version, permGWAS2 allows for two different permutation strategies:

x (default): permute the fixed effects matrix including covariates and the SNP of interest (equivalent to permuting y and the covariance matrix)

y: permute only the phenotype vector (same method as in the original permGWAS)

Details on the architecture of permGWAS and permGWAS2, benchmarking results of the framework and on permutation-based thresholds can be found in our publications.

## Publications & Citation

John, M., Korte, A. & Grimm, D. G. (2023). 
**permGWAS2: Enhanced and Accelerated Permutation-based Genome-wide Association Studies**. bioRxiv.

DOI: https://doi.org/10.1101/2023.11.28.569016

John, M., Ankenbrand, M. J., Artmann, C., Freudenthal, J. A., Korte, A., & Grimm, D. G. (2022).  
**Efficient Permutation-based Genome-wide Association Studies for Normal and Skewed Phenotypic Distributions**.  
Bioinformatics, 2022.   

DOI: https://doi.org/10.1093/bioinformatics/btac455


## How to run permGWAS2
1. [Requirements & Installation](./docs/INSTALLATION.md)
2. [Quickstart Guide](./docs/QUICKSTART.md)
3. [Data Guide](./docs/DATAGUIDE.md)
4. [permGWAS2 with permutations](./docs/PERMUTATIONS.md)
5. [Create plots](./docs/PLOTS.md)
6. [Optional settings](./docs/OPTIONS.md)
