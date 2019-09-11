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
1. Estimate effect size
2. Estimate PRS
3. Test for coordinated interactions

## Dpendencies
### R packages:
* ukbtools
### Python
* 