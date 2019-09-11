#!/bin/bash
# UKBB Data fields:
# 22000 - Genotyping batch
# 22009 - Genetic principal components
# 21001 - BMI
# 21003 - Age
# 31    - Sex
# 22011 - Related samples
# 22001 - Mismatching sex

# Extract basic covariates
RSCRIPT=`python parameters.py RSCRIPT`
WD=`python parameters.py WD`
mkdir -p $WD/CSV
$RSCRIPT src/ukb_extract_phenotypes.R 31 22009 22000 21001 21003 \
  -o=${WD}/CSV/covariates.cov --remove_non_caucasian --remove_related --remove_sex_mismatch

# Make dummy covariates (for centers and genotyping batch)

# Merge covariates

# Normalize covariates
# Will end with the next files:
# 1. covariates2.cov
# 2. covariates_standarized_sex