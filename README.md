[![Python 3.8](https://img.shields.io/badge/python-3.8%20%7C%203.9%20%7C%203.10-blue)](https://www.python.org/downloads/release/python-3100/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<img src="/permGWAS_logo.png" data-canonical-src="/permGWAS_logo.png" height="80" />  

## permGWAS2

This is an improved version of permGWAS. The original version can be found at [permGWAS Version1](https://github.com/grimmlab/permGWAS/releases/tag/permGWAS)

permGWAS2 is an open source software tool written in python to efficiently perform genome-wide association studies (GWAS)
with permutation-based thresholds. permGWAS2 provides support for multiple CPUs as well as for GPUs. 

In contrast to the original version, permGWAS2 allows for two different permutation strategies:

x (default): permute the fixed effects matrix including covariates and the SNP of interest (equivalent to permuting y and the covariance matrix)

y: permute only the phenotype vector (same method as in the original permGWAS)

Details on the architecture of permGWAS, benchmarking results of the framework and on permutation-based thresholds can be found in our publication.

## Publication & Citation

John, M., Korte, A. & Grimm, D. G. (2023). 
**permGWAS2: Enhanced and Accelerated Permutation-based Genome-wide Association Studies**. bioRxiv.

DOI: https://doi.org/10.1101/2023.11.28.569016

John, M., Ankenbrand, M. J., Artmann, C., Freudenthal, J. A., Korte, A., & Grimm, D. G. (2022).  
**Efficient Permutation-based Genome-wide Association Studies for Normal and Skewed Phenotypic Distributions**.  
Bioinformatics, 2022.   

DOI: https://doi.org/10.1093/bioinformatics/btac455


## Requirements

To ensure a stable working environment, we recommend using [docker](https://www.docker.com). To follow this recommendation, 
docker needs to be installed and running on your machine. We provide a Dockerfile based on CUDA 11.5 and Ubuntu 20.4.

## Installation

Clone the repository into the directory where you want to set up the project

```shell
git clone https://github.com/grimmlab/permGWAS.git
```

Navigate to `config` and build a Docker image using the provided Dockerfile

```shell
cd config
docker build -t IMAGENAME .
```

Run an interactive Docker container based on the created image.\
You have to mount the directory where the repository is located on your machine in the Docker container. 
If you want to work on GPU, specify the GPUs to mount.

```shell
docker run -it -v PATH/TO/REPO/FOLDER:NAME/OF/DIRECTORY/IN/CONTAINER --gpus=all --name CONTAINERNAME IMAGENAME
```


## Execution

In the Docker container, navigate to the repository in the mounted directory.

Run the script with the test data provided

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv 
```

## Config file

permGWAS2 accepts yaml config files where you can specify all flags and options instead of passing them all separately.

```shell
python3 permGWAS.py -config ./data/config.yaml 
```

## Input data
The minimal requirement is to provide a genotype and a phenotype file. We provide test data in the folder `data`.
permGWAS2 is designed to work with several genotype file formats:

### HDF5/H5/H5PY
The file has to contain the following keys:

- snps: genotype matrix, additively encoded (012)
- sample_ids: vector containing corresponding sample ids
- position_index: vector containing the positions of all SNPs
- chr_index: vector containing the corresponding chromosome number

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv 
```

We recommend to use permGWAS2 with HDF5/H5/H5PY files. For this we provide a function to create an H5 file which 
satisfies our requirements and takes CSV, PLINK and binary PLINK genotype files as an input. For more info on 
how to use this function, see the section create H5 file below.

### CSV
The first column should be the sample ids. The column names should be the SNP identifiers in the form "CHR_POSITION"
(e.g. Chr1_657). The values should be the genotype matrix in additive encoding. 

```shell
python3 permGWAS.py -x ./data/x_matrix.csv -y ./data/y_matrix.csv 
```

### PLINK
To use PLINK data, a .map and .ped file with the same prefix need to be in the same folder. 
To run permGWAS2 with PLINK files, you can use PREFIX.map or PREFIX.ped as option for the genotype file.

```shell
python3 permGWAS.py -x ./data/x_matrix.map -y ./data/y_matrix.pheno 
```

### binary PLINK
To use binary PLINK data, a .bed, .bim and .fam file with the same prefix need to be in the same folder. 
To run permGWAS2 with binary PLINK files, you can use PREFIX.bed, PREFIX.bim or PREFIX.fam as option for the genotype file.


### phenotype file 
permGWAS2 currently only accepts CSV, PHENO and TXT files for the phenotype. Here the first column should contain 
the sample ids. The remaining columns should contain the phenotype values with the phenotype name as column name. 
For TXT and PHENO files it is assumed that the values are separated by a single space. 
It is possible to run permGWAS with several traits one after another as long as they are stored in the same 
phenotype file.

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -trait phenotype_value phenotype_2
```
You can also run permGWAS2 for all available phenotypes in your phenotype file:

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -trait all
```

### kinship file
Per default permGWAS2 computes the realized relationship kernel as kinship matrix. 
It is also possible to provide a kinship matrix. Currently, permGWAS only accepts CSV, H5, HDF5, H5PY files as 
kinship file. For CSV files the first column should contain the sample ids. For H5, HDF5, H5PY files the kinship 
matrix should have the key 'kinship' and the corresponding sample ids the key 'sample_ids'
The sample ids need to match those of the genotype matrix.

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -k ./data/k_matrix.csv
```

### covariates file
It is possible to run permGWAS2 with covariates. If no covariates file is provided, only the intercept will be used as 
fixed effect. Currently, permGWAS2 only accepts CSV files for covariates. Here the first column should contain the 
sample ids. The sample ids must match those of the phenotype file.

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -cov ./data/cov_matrix.csv
```

## create H5 file
We provide a function to create an H5 file which satisfies our requirements. It is possible to create the H5 based on a 
CSV, PLINK or binary PLINK files which have to fulfil the same requirements as above. The function takes the genotype 
file path via the option -x and additionally one can specify a new directory to save the H5 file via -sd if the save 
directory is not specified, the new file will be stored in the same directory as the input file.

```shell
python3 create_h5_file.py -x ./data/x_matrix.map -sd ./data/test
```

## output files
Per default permGWAS2 creates a CSV output file and saves it in a directory called `results`. One can also specify a 
different directory for the output files via the flag `--out_dir`. The output file will be saved under the name
`p_values_NAME.csv`, where NAME will be the phenotype name per default, but can also be changed via `--out_file`.

The result file contains for each analyzed SNP:
- CHR: chromosome
- POS: position
- p_value: p-value
- test_stat: test statistic
- maf: minor allele frequency
- SE: standard error
- effect_size: coefficient beta

Additionally, a TXT file with summary statistics will be saved. 
This file contains the estimates of the variance components of the null model,
 the narrow-sense heritability, the Bonferroni threshold and, 
if activated, the permutation-based threshold.

## further options
### minor allele frequency (MAF)
It is possible to filter the markers for minor allele frequency. For this use the flag `-maf a`, where a should be any 
integer value between 0 and 30. Per default permGWAS2 does not filter for MAF (i.e. `-maf 0`).

### permutations
By including the flag `-perm q` for any integer number q>0, permGWAS2 first performs a normal GWAS and afterwards 
creates q permutations of the phenotype and computes permutation-based p-values as well as minimal p-values that can be 
used to compute a permutation-based threshold via the maxT/minP method. For more details on permutation-based p-values
and thresholds see our paper published alongside this software.
When performing permGWAS2 with permutations, the `p_values_NAME.csv` output file contains an additional column 
'adjusted_p_val', which contains the permutation-based p values. Additionally, a second file `min_p_values_NAME.csv` 
will be created, containing for each permutation the seed and the minimal p-value.

### GPU usage
For faster computations, permGWAS2 supports GPU usage. If one or several GPUs are available permGWAS2 will per default use 
the GPU device 0 for its computations. Via the flag `-device` the GPU device can be changed. 

### batch size
It is possible to adjust the batch size for the simultaneous computation of univariate tests via `-batch`. Here the 
default is set to 50000. If you run into memory errors while using permGWAS2 we suggest reducing the batch size.

When using permGWAS2 with permutations, several univariate tests will be computed for all permutations at once. 
To prevent running into memory errors, one can adjust the batch size for permutations separately via ``-batch_perm``.
Here the default value is set to 1000. We suggest adjusting this parameter depending on the number of samples and number
of permutations.

### batch-wise loading of genotype
As memory is a limiting factor, permGWAS2 is also capable to load the genotype matrix batch-wise from file under certain 
conditions. For this you have to provide a precomputed kinship matrix (see input data above) and the genotype matrix 
must be provided via an HDF5 file (see above for a function to create an HDF5 file). Currently, it is not possible to 
use MAF filtering when loading batch-wise. 

However, if memory is not an issue, we recommend loading the genotype file completely to improve the speed of permGWAS2.
When no precomputed kinship is provided or MAF filtering is used, the genotype matrix will be loaded completely per 
default. It is also possible to force permGWAS2 to load the genotype matrix completely even if a kinship is provided
via the flag `-load_genotype`. 

### Manhattan plot
When using the flag `-mplot`, permGWAS2 will generate and save a Manhattan plot with Bonferroni significance threshold. 
If permGWAS2 is run with permutations, additionally the permutation-based threshold will be plotted.

### QQ-plot
When using the flag `-qqplot`, permGWAS2 will generate and save a QQ-plot including the inflation factor lambda.


### Flags and options
|||
|---|---|
|-x (--genotype_file) | absolute or relative path to genotype file |
|-y (--phenotype_file) | absolute or relative path to phenotype file |
|-trait (--y_name)| name of phenotype (column) to be used in phenotype file, optional, default is "phenotype_value"|
|-k (-kinship_file) | absolute or relative path to kinship file, optional|
|-cov (--covariate_file) | absolute or relative path to covariates file, optional|
|-cov_list (--covariate_list) | names of covariates to use from covariate_file, optional |
|-maf (--maf_threshold) | minor allele frequency threshold as percentage value, optional, default is 0|
|-load_genotype | choose whether to load full genotype from file or batch-wise during computations, optional, default is False|
|-config (--config_file) | full path to yaml config file|
|-model | specify model name, only relevant if you define your own models, currently only lmm is available|
|-out_dir | name of the directory result-files should be stored in, optional, if not provided, files will be stored in folder "results" in current directory|
|-out_file | NAME of result files, will be stored as NAME_p_values and NAME_min_p_values, optional, if not provided name of phenotype will be used|
|-disable_gpu | use if you want to perform computations on CPU only though GPU would be available| 
|-device | GPU device to be used, optional, default is 0|
|-perm | number of permutations to be performed, optional, default is 0|
|-perm_method | method to use for permutations: y - permute only y, x - permute y and kinship matrix, default is x|
|-adj_p_value | additionally compute permutation-based adjusted p-values and store them in the p-value file, optional default is False|
|-batch (--batch_size) | number of SNPs to work on simultaneously, optional, default is 50000|
|-batch_perm (--perm_batch_size) | number of SNPs to work on simultaneously while using permutations, optional, default is 1000|
|-mplot (--plot, --manhattan)| creates Manhattan plot, optional|
|-qqplot | creates QQ-plot, optional|
|-not_add | use when genotype is not in additive encoding| 