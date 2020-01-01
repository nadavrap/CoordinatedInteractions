# Coordinated interaction
This is a collection of scripts in python, R and shell scripts that were used in the analysis and data generation of 
the manuscript "Coordinated Interaction: A model and test for globally signed epistasis in complex traits", Sheppard et al. 2020.

## Software
### Required:
* Plink1.9
* Plink2
* Python3.6
* R3.6.1

Scripts were developed under the assumption of an available SGE cluster as many jobs are been submitted. <> (In case this is not the case, please contact the authors.)

### Optional:
* SAGE
* PRsice
* BOLT-LMM

## Data
We used UK biobank imputed data. You have to apply for, and it cannot be shared.
* Imputed genomic data
* Phenotipic data

## Generating PRS and tissue-specific PRS (contact: Nadav Rappoport)
1. Extract covariates from UKBB data `src/covariates_generate.sh`
1. Estimate effect size
    ```shell script
    python3 CoordinatedInteractions/src/even_odd.py --nfolds 10 --imputed --normalized --AllPCs 
    --phenotype T2D --betagenerator Plink --predictor None
    ```
1. Summarize effect sizes from file per chromosome to file by fold: 
    ```shell script
    f=T2D
    for i in `seq 0 9`; do
        # For each fold, first copy chromosome 1 with header, then append other chromosomes w/o headers
        cp ${f}/plink_glm_clumped_chr1_fold_not${i}_covariates_standarized_sex_MAF0.001_*.${f}.glm.logistic \
            ${f}/norm_plink_glm_imp_fold_${i}.tsv
        for c in `seq 2 22`; do
            cat ${f}/plink_glm_clumped_chr${c}_fold_not${i}_covariates_standarized_sex_MAF0.001_*.${f}.glm.logistic | \
                tail -n+2 >> ${f}/norm_plink_glm_imp_fold_${i}.tsv 
        done
    done
    ```
1. Estimate PRS
    ```shell script
    PVALUES="0.01 0.001
    python3 src/plink_score.py -f $f -p $PVALUES
    ```
2. Summarize over folds
    ```shell script
    src/plink_score.py -f $f -p $PVALUES --summarize
    ```
3. Test for coordinated interactions

When external summary statistics (variants' effect sizes) is taken from an esxternal source, and there is no need to 
use cross validation in order to estimate it, you can skip on the relevant step.

When analyzing for a subset of variants (e.g, tissue-specific), same scripts and similar commands can be used. 
See example in [Tissue_specific_example.md](Tissue_specific_example.md). 

## Estimating Coordinated Interaction (contact: Andy Dahl)

### Generating Figure 2 and Table 2: Even/Odd estimator for CI in UKB

This is performed by the R scripts in test_ci/ukb_eo. These scripts also generate Supp Figs 1 and 2, and Supplementary Table 1 (called table2.csv), a more complete version of Table 2.


### Generating Figure 3 and Table 3: Tissue-specific estimator for CI in UKB

This is performed by the R scripts in test_ci/ukb_tiss. These analyses depend on tissue/cell type-specific genomic annotations. These were obtained from https://github.com/bulik/ldsc, which were derived by (Finucane et al., 2018) using previously generated data (Fehrmann et al. and 2015; Pers et al., 2015).


### Generating Figure 4: Tissue pair-specific estimator for CI in UKB

This is performed by the R scripts in test_ci/ukb_tissXtiss.


## Dependencies
### R packages:
* ukbtools, argparse, preprocessCore
### Python
* math, argparse, glob, numpy, pandas, shlex, shutil, socket, subprocess, tempfile
