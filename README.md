# Simpy: Calculate imputation accuracy by comparing imputation results to WGS data.

Documentation and scripts provided for calculating imputation accuracy. Should work with any imputed VCF file that has GT or GT+DS format fields.

## Requirements

- bcftools: used to calculate MAFs
- python v3: source code was implemented and tested on python 3.6
  - `pandas`
  - `numpy`
  - `cyvcf2`

## Required command line arguments are:

The following inputs, in vcf.gz format, including its respective tabix .tbi file, are required to run.

- imputed: imputation results
- wgs: ground truth file, containing experimentally determined genotypes (i.e. Whole Genome Sequencing data)
- bwgs: same wgs file but in BCF format to speed up the process and .csi index file associated.

**NB PREREQUISITE:**
- All files provided must be in vcf.gz format (compressed, tabixed). 
- Alleles must match, in other words: no Swap, no flips, same build.
- It is not necessary to provide allele frequencies, since the tool will calculate it internally using bcftools.
- the tools works also if wgs data in BCF format are not provided, but will be slower.

## Usage:

The help shows all the required arguments listed above, plus optional arguments.

```
python3 Simpy.py -h
usage: Simpy.py --ga <input_genotype_array.vcf.gz> --imputed <imputed_file.vcf.gz> --wgs <whole_genome_file.vcf.gz>
Use -h or --help to display help.

arguments:
  -h, --help            show this help message and exit
  --ga GA               optional, path to genotype array file in vcf.gz format, with tbi
  --wgs WGS             required, path to whole genome file in vcf.gz format, with tbi
  --imputed IMPUTED     required, path to imputed file in vcf.gz format, with tbi
  --ref REF             optional, path to reference panel file in vcf.gz
                        format, with tbi. Used for MAF calculation. WGS file
                        will be used if no reference file is provided.
  --max_total_rows MAX_TOTAL_ROWS
                        optional, maximun number of rows or variants to be loaded
                        simultaneously, summing all chunks loaded by all cores
  --max_per_core MAX_PER_CORE
                        optional, maximun number of variants per chunk per core, lower
                        it to avoid RAM overload
  --min_per_core MIN_PER_CORE
                        optional, minimun number of variants per chunk per core,
                        increase to avoid interprocess communication overload
  --sout SOUT           optional, output file path/name per sample, default is
                        the same as the imputed file with
                        _per_sample_results.txt suffix
  --vout VOUT           optional, output file path/name per variant, default is
                        the same as the imputed file with
                        _per_variant_results.txt suffix
```

A detailed report with accuracy ratio, F1 score, Pearson correlation (r2) is generated and wrote to the output file (i.e accuracy_result.txt)

## How to run example:
```
python3.6 Simpy.py --imputed imputed.vcf.gz --wgs WGS.vcf.gz --bwgs WGS.bcf.gz --bimputed imputed.bcf.gz
```

## Results:

```
Processing  239 imputed samples
Processing chunk: 1 Max rows per chunk: 10000
1 Read imputed file time:  0.0004201233386993408
2 Chunking time:  8.239876478910446e-06
3 Calculation time:  0.16841873712837696
4 Merging calculations per sample time:  0.0037479093298316
Results per sample at: imputed.vcf.chr1_per_sample_results.txt
Results per variant at: imputed.vcf.chr1_per_variant_results.txt
Total run time (sec): 0.17389075085520744

```

The results will be displayed as the example bellow (per variant):
```
position        SNP     REF_MAF IMPUTED_MAF     WGS_MAF F-score concordance_P0  IQS     r2      precision       recall  TP      TN      FP      FN TP_ratio TN_ratio        FP_ratio        FN_ratio        RMSE
22:48742072     22:48742072_G_A 0.0003  0.0013  0.0093  0.992   0.984   0.246   0.143   1.0     0.984   754     738     0       12      1.0     0.984       0.0     0.016   0.126
22:48742498     22:48742498_G_A 0.0028  0.0     0.002   0.998   0.996   0.0     0.014   1.0     0.996   752     749     0       3       1.0     0.996       0.0     0.004   0.063
22:48743012     22:48743012_C_A 0.061   0.0578  0.0665  0.947   0.882   0.479   0.3     0.946   0.948   804     610     46      44      0.946   0.933       0.054   0.067   0.307
22:48742525     22:48742525_G_A 0.0041  0.0013  0.0033  0.998   0.996   0.57    0.36    1.0     0.996   754     747     0       3       1.0     0.996       0.0     0.004   0.066
22:48741511     22:48741511_C_T 0.0009  0.0007  0.0007  1.0     1.0     1.0     0.999   1.0     1.0     753     751     0       0       1.0     1.00.0      0.0     0.002
22:48741402     22:48741402_A_G 0.2096  0.0439  0.2646  0.836   0.621   0.233   0.185   0.879   0.796   814     370     112     208     0.879   0.640.121   0.36    0.591
22:48742178     22:48742178_G_A 0.0089  0.0007  0.0153  0.986   0.971   0.081   0.062   1.0     0.972   753     729     0       22      1.0     0.971       0.0     0.029   0.167
22:48741385     22:48741385_T_C 0.3718  0.113   0.4149  0.863   0.547   0.192   0.168   0.768   0.984   1078    83      325     18      0.768   0.822       0.232   0.178   0.656
22:48741702     22:48741702_C_T 0.001   0.0047  0.0293  0.976   0.956   0.287   0.165   0.995   0.958   755     712     4       33      0.995   0.956       0.005   0.044   0.238
```

Results per sample:
```
imputed_ids     WGS_ids F-score concordance_P0  r2      precision       recall  TP      TN      FP      FN      TP_ratio        TN_ratio        FP_ratio    FN_ratio        RMSE
A00003_A00003   A00003_A00003   0.976   0.95    0.497   0.972   0.981   103.0   94.0    3.0     2.0     0.972   0.979   0.028   0.021   0.152
A00018_A00018   A00018_A00018   0.932   0.851   0.696   0.972   0.896   103.0   84.0    3.0     12.0    0.972   0.875   0.028   0.125   0.296
A00056_A00056   A00056_A00056   0.99    0.98    0.613   0.981   1.0     104.0   96.0    2.0     0.0     0.981   1.0     0.019   0.0     0.162
A00080_A00080   A00080_A00080   0.986   0.97    0.671   0.981   0.991   106.0   93.0    2.0     1.0     0.981   0.989   0.019   0.011   0.141
A00083_A00083   A00083_A00083   0.99    0.98    0.56    0.981   1.0     104.0   96.0    2.0     0.0     0.981   1.0     0.019   0.0     0.128
A00099_A00099   A00099_A00099   0.99    0.98    0.139   0.99    0.99    102.0   98.0    1.0     1.0     0.99    0.99    0.01    0.01    0.138
A00120_A00120   A00120_A00120   0.976   0.96    0.5     0.962   0.99    101.0   96.0    4.0     1.0     0.962   0.99    0.038   0.01    0.21
A00146_A00146   A00146_A00146   0.972   0.941   0.64    0.963   0.981   103.0   93.0    4.0     2.0     0.963   0.979   0.037   0.021   0.176
A00152_A00152   A00152_A00152   0.986   0.97    0.403   0.981   0.99    103.0   96.0    2.0     1.0     0.981   0.99    0.019   0.01    0.17
```

The results can be interpreted as follows.

Metrics per variant:
- REF_MAF: Reference Panel MAF (if reference panel is provided)
- IMPUTED_MAF: Imputed MAF
- WGS_MAF: Whole Genome MAF
- F-score: macro F-score,
- concordance_P0: accuracy ratio ,
- IQS: imputation quality score
- precision: precision
- recall: recall
- TP: true positives
- TN: true negatives
- FP: false positives
- FN: false negatives
- TP_ratio: true positives ratio
- TN_ratio: true negatives ratio
- FP_ratio: false positives ratio
- FN_ratio: false negatives ratio
- RMSE: root mean squared error

Metrics per sample:
- F-core: F-score per sample
- concordance_P0: accuracy ratio
- r2: r-squared
- precision: precision
- recall: recall
- TP: true positives
- TN: true negatives
- FP: false positives
- FN: false negatives
- TP_ratio: true positives ratio
- TN_ratio: true negatives ratio
- FP_ratio: false positives ratio
- FN_ratio: false negatives ratio
- RMSE: root mean squared error


This message will show up if everything worked:

```
Processing chunk: 1 Max rows per chunk: 10000
...
Processing chunk: 52 Max rows per chunk: 10000
1 Read imputed file time:  5.211546601727605
2 Chunking time:  0.0284710843116045
3 Calculation time:  90.96355835348368
4 Merging calculations per sample time:  20.620813813060522
Results per sample at: imputed_file_100samples.vcf_per_sample_results.txt
Results per variant at: imputed_file_100samples.vcf_per_variant_results.txt
Total run time (sec): 117.97938374243677

```

After running this example, you can visualize the results for the test sample data in:
- imputed_file_per_sample_results.txt
- imputed_file_per_variant_results.txt

For example:

```
head imputed.vcf_per_variant_results.txt

position        SNP     IMPUTED_MAF     WGS_MAF F-score concordance_P0  IQS     r2      precision       recall  TP      TN      FP      FN TP_ratio TN_ratio        FP_ratio        FN_ratio        RMSE
22:16050783     22:16050783_A_G 0.0     0.0     1.0     1.0     0.0     0.0     1.0     1.0     100     100     0       0       1.0     1.00.0      0.0     0.001
22:16050922     22:16050922_T_G 0.0     0.0     1.0     1.0     0.0     0.0     1.0     1.0     100     100     0       0       1.0     1.00.0      0.0     0.0
22:16050984     22:16050984_C_G 0.0     0.0     1.0     1.0     0.0     0.0     1.0     1.0     100     100     0       0       1.0     1.00.0      0.0     0.0
22:16051269     22:16051269_G_T 0.0     0.0     1.0     1.0     0.0     0.0     1.0     1.0     100     100     0       0       1.0     1.00.0      0.0     0.0
22:16051477     22:16051477_C_A 0.0     0.0     1.0     1.0     0.0     0.0     1.0     1.0     100     100     0       0       1.0     1.00.0      0.0     0.006
22:16052240     22:16052240_C_G 0.0     0.0     1.0     1.0     0.0     0.0     1.0     1.0     100     100     0       0       1.0     1.00.0      0.0     0.0
22:16052271     22:16052271_G_A 0.0     0.0     1.0     1.0     0.0     0.0     1.0     1.0     100     100     0       0       1.0     1.00.0      0.0     0.001
22:16052428     22:16052428_G_A 0.0     0.0     1.0     1.0     0.0     0.0     1.0     1.0     100     100     0       0       1.0     1.00.0      0.0     0.0
22:16052639     22:16052639_C_T 0.0     0.0     1.0     1.0     0.0     0.0     1.0     1.0     100     100     0       0       1.0     1.00.0      0.0     0.0

```

## How to run example with Big Input file to avoid Memory error:
```
./imputation_accuracy.sh -i imputed.vcf.gz -w WGS.vcf.gz -bw WGS.bcf.gz -t 4
```

This will create multiple chunks of the VCF imputed in input and analyze each one with the Simpy.py software and at the end will use the rebuild_metrics.py script to re-calculate the results. 

This message will show up if everything worked:

```
Imputed VCF   : imputed_HG00479.vcf.gz
WGS VCF       : wgs_HG00479.vcf.gz
WGS BCF       : wgs_HG00479.bcf.gz
Threads       : 4
Skip chunks   : False
Skip analysis : False
Splitted Imputed file in chuncks of [100k]
BCF Imputed files Created
Joining files...
Deleting tmp files...
Process Completed.
${sample_name}_per_sample_results.tsv.gz
${sample_name}_per_variant_results.tsv.gz"
```

## References:

[1] Ramnarine S, Zhang J, Chen LS, Culverhouse R, Duan W, Hancock DB, Hartz SM, Johnson EO, Olfson E, Schwantes-An TH, Saccone NL. When does choice of accuracy measure alter imputation accuracy assessments?. PloS one. 2015 Oct 12;10(10):e0137601.
