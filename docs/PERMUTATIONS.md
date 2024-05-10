# permGWAS2 with permutations

The main purpose of permGWAS2 is to perform GWAS with permutation-based thresholds. To use permGWAS2 with permutations, 
you have to specify the number of permutations *q* via the flag `-perm`:
```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -perm 100
```
This creates an additional result file `min_p_values_NAME.csv` containing for each permutation the seed and the minimal 
p-value. Additionally, the `summary_statistics_NAME.txt` output file now contains permutation-based significance 
thresholds for common significance levels $\alpha$.

### General workflow of permGWAS2 with permutations
1. Compute p-values for all available SNPs during normal GWAS run
2. Create *q* permutations
3. Compute the test statistic for each permutation and SNP in batches
4. For each permutation find the maximal test statistic over all SNPs and compute the corresponding minimal p-value
5. The permutation-based threshold is given as the ($1-\alpha$)th percentile for a significance level $\alpha$ 
(*maxT/minP method*)

### Additional settings
- permGWAS2 supports two different permutation strategies which can be selected via the flag `-perm_method`:
  1. `x`(default): permutes the fixed effects matrix including SNP of interest and covariates (equivalent to permuting 
  the phenotype and covariance matrix). This method considers the population structure while permuting.
  2. `y`: only permute the phenotype vector. This method is faster but breaks the population structure between the 
  samples 
```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -perm 100 -perm_method x

python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -perm 100 -perm_method y
```
- permGWAS2 supports computations on GPUs. If GPUs are available, it will automatically use the 0th GPU. If no GPUs are 
available, permGWAS will perform all computations on CPUs only. To change the GPU you can use the flag `-device` and 
specify the number of the GPU to use. If you do NOT want to use GPUs, although they are available, you can use the flag 
`disable_gpu`:
```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -perm 100 -device 1

python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -perm 100 -dasable_gpu
```
- Since permGWAS2 computes the test statistics for different SNPs and permutations simultaneously in batches, the 
available VRAM poses a limitation. To avoid running into memory errors (when using GPUs), you can manually adjust the 
batch-size, i.e. the number of SNPs to be processed simultaneously for all permutations, via the flag `-batch_perm` 
(The default are 1000 SNPs):
```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -perm 100 -batch_perm 500
```
- permGWAS is also able to compute permutation-based adjusted p-values and save them in the p_value output file via the 
flag `adj_p_value`. However, it should be noted that in order to get meaningful adjusted p-values, millions of 
permutations are needed. 
```shell
python3 permGWAS.py -x ./data/x_matrix.h5 -y ./data/y_matrix.csv -perm 100 -adj_p_value
```
