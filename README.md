# Coordinated interaction
This is a collection of scripts in python, R and shell scripts that were used in the analysis and data generation of 
the manuscript "PRS Coordinated interactions".

## Software
### Required:
* Plink1.9
* Plink2
* Python3.6
* R3.6.1

Scripts were developed under the assumption of an available SGE cluster as many jobs are been submitted. In case this 
is not the case, please contact the authors.

### Optional:
* SAGE
* PRsice
* BOLT-LMM

## Data
We used UK biobank imputed data. You have to apply for, and it cannot be shared.
* Imputed genomic data
* Phenotipic data

## Analysis is build from several steps with many details
1. Extract covariates from UKBB data `src/covariates_generate.sh`
1. Estimate effect size `python3 CoordinatedInteractions/src/even_odd.py --nfolds 10 --imputed --normalized --AllPCs 
--phenotype T2D --betagenerator Plink --predictor None`
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
2. Estimate PRS
```shell script
PVALUES="0.01 0.001
python3 src/plink_score.py -f $f -p $PVALUES
```
2. Summarize over folds
```shell script
src/plink_score.py -f $f -p $PVALUES --summarize
```
3. Test for coordinated interactions

## Dpendencies
### R packages:
* ukbtools
### Python
* 