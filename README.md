# permGWAS
Permutation-based Linear Mixed Models for GWAS

## Introduction
- GWAS allg - LMM
- permutationen
- 3D/4D architektur
- perm p-val, perm-thr
- output file



## Requirements

To ensure a stable working environment, we recommend using [docker](https://www.docker.com). To follow this recommendation, 
docker needs to be installed and running on your machine. We provide a Dockerfile based on CUDA 11.1 and Ubuntu 20.4.


## Installation

Clone the repository into the directory where you want to set up the project

```shell
git clone https://github.com/grimmlab/permGWAS.git
```

Navigate to `config` and build a Docker image using the provided Dockerfile

```shell
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


## Input data
The minimal requirement is to provide a genotype and a phenotype file. We provide test data in the folder `data`.
permGWAS is designed to work with several genotype file formats:

### .hdf5/.h5/.h5py
The file has to contain the following keys:

- snps: genotype matrix, 012 encoded
- sample_ids: vector containing corresponding sample ids
- position_index: vector containing the positions of all SNPs
- chr_index: vector containing the corresponding chromosome number

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv 
```

### .csv
The first column should be the sample ids. The column names should be the SNP identifiers in the form "CHR_POSITION"
(e.g. Chr1_657). The values should be the genotype matrix in additive encoding. 

```shell
python3 permGWAS.py -x ./data/x_matrix.csv -y ./data/y_matrix.csv 
```

### PLINK
To use PLINK data, a .map and .ped file with the same prefix need to be in the same folder. 
To run permGWAS with PLINK files, you can use PREFIX.map or PREFIX.ped as option for the genotype file.

```shell
python3 permGWAS.py -x ./data/x_matrix.map -y ./data/y_matrix.pheno 
```

### binary PLINK
To use binary PLINK data, a .bed, .bim and .fam file with the same prefix need to be in the same folder. 
To run permGWAS with binary PLINK files, you can use PREFIX.bed, PREFIX.bim or PREFIX.fam as option for the genotype file.


### phenotype file 
permGWAS currently only accepts .csv, .pheno and .txt files for the phenotype. Here the first column should contain the 
sample ids. The remaining columns should contain the phenotype values with the phenotype name as column name. 
For .txt files it is assumed that the values are separated by a single space.


### kinship file
Per default permGWAS computes the realized relationship kernel as kinship matrix. 
It is also possible to provide a kinship matrix. Currently, permGWAS only accepts .csv, .h5, .hdf5, .h5py files as 
kinship file. For .csv files the first column should contain the sample ids. For .h5, .hdf5, .h5py files the kinship 
matrix should have the key 'kinship' and the corresponding sample ids the key 'sample_ids'
The sample ids need to match those of the genotype matrix.

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv --k ./data/k_matrix.csv
```

### covariates file
It is possible to run permGWAS with covariates. If no covariates file is provided, only the intercept will be used as 
fixed effect. Currently, permGWAS only accepts .csv files for covariates. Here the first column should contain the 
sample ids. The sample ids must match those of the phenotype file.

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv --cov ./data/cov_matrix.csv
```

### Options
|||
|---|---|
|-x (-x_file, -genotype) | absolute or relative path to genotype file |
|-y (-y_file, -phenotype) | absolute or relative path to phenotype file |
|--y_name (--pt_name)| name of phenotype (column) to be used in phenotype file, optional, default is "phenotype_value"|
|--k (--k_file, --kinship) | absolute or relative path to kinship file, optional|
|--cov (--cov_file) | absolute or relative path to covariates file, optional|
|--load_genotype (--load_x) | choose whether to load full genotype file or not, optional, default is False|
|--maf | minor allele frequency threshold as percentage value, optional, default is 0|
|--perm | number of permutations to be performed, optional, default is 0|
|--out_dir | name of the directory result-files should be stored in, optional, if not provided, files will be stored in folder "results" in current directory|
|--out_file | NAME of result files, will be stored as NAME_p_values and NAME_min_p_values, optional, if not provided name of phenotype will be used|
|--device | GPU device to be used, optional, default is 0|
|--batch | number of SNPs to work on simultaneously, optional, default is 50000|
|--batch_perm | number of SNPs to work on simultaneously while using permutations, optional, default is 1000|
|--plot | creates manhattan plot, optional|
    


## Citation

tbd
