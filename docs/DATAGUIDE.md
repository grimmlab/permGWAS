# Data Guide

The minimal requirement is to provide a genotype and a phenotype file. We provide test data in the folder `data`.
permGWAS2 is designed to work with several genotype file formats:

## Genotype file
permGWAS needs **fully imputed** genotypes. We support our custom HDF5/H5/H5PY file, CSV PLINK and binary PLINK files. 
We recommend to use permGWAS2 with HDF5/H5/H5PY files. For this we provide a function to create an H5 file which satisfies 
our requirements and takes CSV, PLINK and binary PLINK genotype files as an input. For more info on how to use this function, 
see the section **Create H5 file** below.

### HDF5/H5/H5PY
The file has to contain the following keys:

- snps: genotype matrix, additively encoded (012)
- sample_ids: vector containing corresponding sample ids
- position_index: vector containing the positions of all SNPs
- chr_index: vector containing the corresponding chromosome number

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv 
```

### CSV
The **first column** should be the **sample ids**. The **column names** should be the **SNP identifiers** in the form 
"CHR_POSITION" (e.g. Chr1_657). The values should be the genotype matrix in **additive encoding**. 

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


## Phenotype file 
permGWAS2 currently only accepts CSV, PHENO and TXT files for the phenotype. Here the **first column** should contain 
the **sample ids**. The remaining columns should contain the phenotype values with the phenotype name as column name. 
For TXT and PHENO files it is assumed that the values are separated by a **single space**. The samples need not be in 
the same order as in the genotype file. permGWAS2 automatically matched genotype and phenotype and discards all samples 
where only one of both is available.
It is possible to run permGWAS with several traits one after another as long as they are stored in the same 
phenotype file.

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -trait phenotype_value phenotype_2
```
You can also run permGWAS2 for all available phenotypes in your phenotype file:

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -trait all
```

## Kinship file
Per default permGWAS2 computes the realized relationship kernel as kinship matrix. 
It is also possible to provide a kinship matrix. Currently, permGWAS only accepts CSV, H5, HDF5, H5PY files as 
kinship file. For CSV files the first column should contain the sample ids. For H5, HDF5, H5PY files the kinship 
matrix should have the key 'kinship' and the corresponding sample ids the key 'sample_ids'
The sample ids need to match those of the genotype matrix.

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -k ./data/k_matrix.csv
```

## Covariates file
It is possible to run permGWAS2 with covariates. If no covariates file is provided, only the intercept will be used as 
fixed effect. Currently, permGWAS2 only accepts CSV files for covariates. Here the first column should contain the 
sample ids. The sample ids must match those of the phenotype file.

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -cov ./data/cov_matrix.csv
```

## create H5 file
We provide a function to create an H5 file which satisfies our requirements. It is possible to create the H5 based on a 
CSV, PLINK or binary PLINK files which have to fulfil the same requirements as above. The function takes the genotype 
file path via the option `-x` and additionally one can specify a new directory to save the H5 file via `-sd` if the save 
directory is not specified, the new file will be stored in the same directory as the input file.

```shell
python3 create_h5_file.py -x ./data/x_matrix.map -sd ./data/test
```