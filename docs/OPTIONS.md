# Optional settings
## Minor allele frequency (MAF)
It is possible to filter the markers for minor allele frequency. For this use the flag `-maf` and specify an integer 
value between 0 and 30. For example to remove all SNPs with MAF<10%:
```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -maf 10
```
Per default permGWAS2 does not filter for MAF.

## GPU usage
For faster computations, permGWAS2 supports GPU usage. If one or several GPUs are available permGWAS2 will per default use 
the GPU device 0 for its computations. If no GPUs are available, permGWAS will perform all computations on CPUs only. 
To change the GPU you can use the flag `-device` and specify the number of the GPU to use. If you do NOT want to use 
GPUs, although they are available, you can use the flag `disable_gpu`:
```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -device 1

python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -dasable_gpu
```

## Batch size
It is possible to adjust the batch size for the simultaneous computation of univariate tests via `-batch`. Here the 
default is set to 50000. If you run into memory errors while using permGWAS2 we suggest reducing the batch size.
```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -batch 10000
```
When using permGWAS2 with permutations, several univariate tests will be computed for all permutations at once. 
To prevent running into memory errors, one can adjust the batch size for permutations separately via `-batch_perm`.
Here the default value is set to 1000. We suggest adjusting this parameter depending on the number of samples and number
of permutations. For more information about permutations see [permGWAS2 with permutations](./PERMUTATIONS.md)
```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -perm 100 -batch_perm 500
```

## Batch-wise loading of genotype
As memory is a limiting factor, permGWAS2 is also capable to load the genotype matrix batch-wise from file under certain 
conditions. For this you have to provide a precomputed kinship matrix (see [DataGuide](./DATAGUIDE.md)) and the genotype matrix 
must be provided via an HDF5 file (see [DataGuide](./DATAGUIDE.md) for a function to create an HDF5 file). 

However, if memory is not an issue, we recommend loading the genotype file completely to improve the speed of permGWAS2.
When no precomputed kinship is provided, the genotype matrix will be loaded completely per default. It is also possible 
to force permGWAS2 to load the genotype matrix completely even if a kinship is provided via the flag `-load_genotype`. 
```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -load_genotype
```

## Model (coming soon)
permGWAS computes test statistics and p-values based on a Linear Mixed Model (LMM). In the future there will be other 
models available. The model can be chosen via `-model`. Currently, only `lmm` is available. 

## Non-additive encoding
permGWAS assumes that the genotypes are in additive encoding (i.e. number of minor alleles) and produces an error if the genotypes 
are encoded differently. If your data is **not additively encoded**, you can use the flag `-not_add`. For example if you 
are working with other data than SNP data. However, our framework was developed for SNP data, and we give no guarantee that it 
works for other purposes. 


See [Quickstart](./QUICKSTART.md), [permGWAS2 with permutations](./PERMUTATIONS.md) and [Create plots](./PLOTS.md) for 
detailed explanations of other flags and options.

## Overview of all flags and options
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