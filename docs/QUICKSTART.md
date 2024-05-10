# Quickstart Guide

## Simple workflow using Docker

1. Create a new Docker container using our [Installation Guide](./INSTALLATION.md) or start an existing container with:

```shell
docker start -i CONTAINERNAME
```

2. Navigate to the directory where the permGWAS2 repository is located:

```shell
cd /REPO_DIRECTORY/permGWAS
```

3. Run the script with the test data provided in the `./data` folder:

```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv 
```

To use permGWAS2 without Docker, simply omit the first step. 


## Basic settings
### 1. Input Data
Details on the supported data types can be found in the [Data Guide](./DATAGUIDE.md). 
###### Genotype & Phenotype
- The minimal requirement is to provide a genotype and a phenotype file (as relative or absolute paths) via the 
flags `-x` and `-y`, respectively. 
- By default, permGWAS assumes that the phenotype in the phenotype file is called `phenotype_value`. You can specify a 
different name via the flag `-trait`:
```shell
python3 permGWAS.py -x PATH_TO_GENOTYPE -y PATH_TO_PHENOTYPE -trait PHENO_NAME
```
- It is possible to run permGWAS2 for several phenotypes located in the same phenotype file one after another. You can 
either specify a list of phenotypes or run permGWAS2 for all available phenotypes in the file by using the key word `all`:
```shell
python3 permGWAS.py -x PATH_TO_GENOTYPE -y PATH_TO_PHENOTYPE -trait PHENO_1 PHENO_2 PHENO_3

python3 permGWAS.py -x PATH_TO_GENOTYPE -y PATH_TO_PHENOTYPE -trait all
```

###### Kinship
By default, permGWAS2 computes the realized relationship kernel as kinship matrix. You can use a pre-computed genomic 
relationship matrix via the flag `-k`:
```shell
python3 permGWAS.py -x PATH_TO_GENOTYPE -y PATH_TO_PHENOTYPE -k PATH_TO_KINSHIP
```

###### Covariates
It is possible to run permGWAS2 with additional covariates. To specify the covariate file, use the flag `cov`. 
By default, this uses all available covariates in the file. If you only want to use certain columns/covariates, you 
have to use the flag `-cov_list` and specify the covariate names as a list:
```shell
python3 permGWAS.py -x PATH_TO_GENOTYPE -y PATH_TO_PHENOTYPE -cov PATH_TO_COVARIATE_FILE

python3 permGWAS.py -x PATH_TO_GENOTYPE -y PATH_TO_PHENOTYPE -cov PATH_TO_COVARIATE_FILE -cov_list COV_1 COV_2 COV_3
```

### 2. Config file
permGWAS2 accepts yaml config files where you can specify all flags and options instead of passing them all separately:
```shell
python3 permGWAS.py -config ./data/config.yaml 
```
The config file should have the following structure:
```YAML
---
genotype_file: "PATH_TO_GENOTYPE"
phenotype_file: "PATH_TO_PHENOTYPE"
trait: "PHENO_NAME"
kinship_file: "PATH_TO_KINSHIP"
covariate_file: "PATH_TO_COVARIATE_FILE"
covariate_list:
    - "COV_1"
    - "COV_2"
    - "COV_3"
```

### 3. Output files
Per default permGWAS2 creates a CSV output file and saves it in a directory called `results`. You can also specify a 
different directory for the output files via the flag `-out_dir`. The output file will be saved under the name
`p_values_NAME.csv`, where NAME will be the phenotype name by default, but can also be changed via `-out_file`.
```shell
python3 permGWAS.py -x PATH_TO_GENOTYPE -y PATH_TO_PHENOTYPE -out_dir RESULT_FILE_DIR -out_file RESULT_FILE_NAME
```
The result file contains for each analyzed SNP:
- CHR: chromosome number
- POS: position within chromosome
- p_value: computed p-value
- test_stat: computed test statistic
- maf: minor allele frequency of SNP
- SE: standard error
- effect_size: coefficient beta

Additionally, a TXT file with summary statistics will be saved. 
This file contains the estimates of the variance components of the null model,
 the narrow-sense heritability, the Bonferroni threshold and, 
if activated, the permutation-based threshold.


## Further options
The table below shows all available flags. For detailed explanations of further flags and options go to 
[permGWAS2 with permutations](./PERMUTATIONS.md), [Create plots](./PLOTS.md) and [Optional settings](./OPTIONS.md). 

|**flag**|**description**|
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