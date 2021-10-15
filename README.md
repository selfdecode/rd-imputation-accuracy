# SelfDecode pipeline to analyze multiple phasing and imputation softwares simultaneously .

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
Usage: ./PIscore.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20>
Use -h or --help to display help.

author: SELFDECODE
contact email : adriano@selfdecode.com

SelfDecode pipeline to analyze multiple phasing and imputation softwares simultaneously 

Parameters:

      -h|--help
              show help
      -i|--input
              input file
      -r|--reference
              full path reference file without extension.
      -t|--threads
              number of cpus to use.
      -o|--output
              output prefix. No extension.
      -c|--chr
              chomosome to analyze allowed [1-22 X]. NO chr prefix.
      -ibeagle|--imp_beagle
              skip Beagle imputation
      -pbeagle|--phase_beagle
              skip Beagle haplotype estimation
      -impute5|--impute5
              skip Impute imputation
      -shapeit|--shapeit
              skip ShapeIT Phasing
      --minimac|--minimac
              skip Minimac imputation
      -eagle|--eagle
              skip Eagle phasing
      -bigref|--BIGREF
              use this option if you get memory allocate error during accuracy evaluation

[base] Usage: ./PIscore.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20>
[skip] Usage: ./PIscore.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20> -ibeagle x -impute5 x
[memo] Usage: ./PIscore.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20> -bigref x
```


## How to run example:
```
Usage: ./PIscore.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20>```
```

## Results:

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

## How to run Accuracy evaluation with Big Input file to avoid Memory error:
```
./imputation_accuracy.sh -i <imputed_phasing_imputation_CombinationSoftware.vcf.gz> -w WGS.vcf.gz -bw WGS.bcf.gz -t 4
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

[1] A. De Marino, A. Mahmoud, M. Bose, B. Novkovic, P. Yazdi. Which software should I use for phasing or imputation? A comparative analysis of current phasing and imputation software availables based on HMM https://docs.google.com/document/d/1i7rZzvK3ksIXI4wcfLpyqsfpHrlToC0buZx4d_-d6us/edit?pli=1#heading=h.af80tl7prv5v
