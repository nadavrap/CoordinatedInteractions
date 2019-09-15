DATA_DIR <- './'
BINARY_TRAITS <- c('disease_ALLERGY_ECZEMA_DIAGNOSED', 'disease_ASTHMA_DIAGNOSED', 'T2D', 'disease_HI_CHOL_SELF_REP',
                   'disease_CARDIOVASCULAR')
PHENOTYPES <- c('blood_EOSINOPHIL_COUNT', 'blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT', 
                'blood_LYMPHOCYTE_COUNT', 'blood_MEAN_CORPUSCULAR_HEMOGLOBIN', 'blood_MEAN_PLATELET_VOL',
                'blood_MEAN_SPHERED_CELL_VOL', 'blood_MONOCYTE_COUNT', 'blood_PLATELET_COUNT',
                'blood_PLATELET_DISTRIB_WIDTH', 'blood_RBC_DISTRIB_WIDTH', 'blood_RED_COUNT',
                'blood_WHITE_COUNT', 'bmd_HEEL_TSCOREz', 'BMI', 
                'cov_EDU_YEARS', 'impedance_BASAL_METABOLIC_RATEz', 'lung_FEV1FVCzSMOKE',
                'Height',
                'Triglycerides', 'LDL_direct', 'Glucose',
                BINARY_TRAITS)
METHODS <- c('Plink', 'BOLT')
TISSUES <- c('Adipocytes', 'Blood.Cells', 'A08.186.211.464.405.Hippocampus',
             'A03.620.Liver', 'A03.734.Pancreas', 'A08.186.211.Brain',
             'A10.690.Muscles')
DEFAULT_COV_FILE <- file.path(DATA_DIR, 'CSV/covariates.cov')

library(grid)

if (!requireNamespace("preprocessCore")) {
  install.packages("BiocManager")
  BiocManager::install()
  BiocManager::install('preprocessCore')
}

#' Get lm model's pvalue
#' 
#' Based on the next post:
#' http://www.gettinggeneticsdone.com/2011/01/rstats-function-for-extracting-f-test-p.html
#' 
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


# Read PRS values
getPRS <- function(phenotype, method, tissue=NA) {
  # Tissue-specific files have different name convention.
  if (is.na(tissue)) {
    #basename <- paste0('Summary_', tolower(method), '_prsice_', phenotype, '.txt')
    basename <- paste0(phenotype, '_PRS.txt')
  } else {
    basename <- paste0(phenotype, '_', tissue, '_', tolower(method), '_geno_PRS_summary.txt')
  }
  fname <- file.path(DATA_DIR, phenotype, basename)
  prs <- read.table(fname, header = TRUE, sep='', fill=TRUE,
                    stringsAsFactors = FALSE)
  # prs[,3] <- as.numeric(prs[,3])
  # Replace phenotype name with a generic column name
  # colnames(prs)[3] <- 'Pheno'
  prs
}

# Get covariates
# Need to modify some variables, so cache it in .rds file for next time.
getCovars <- function(fname) {
  cachefname <- paste0(fname, '.rds')
  if(file.exists(cachefname)) {
    return(readRDS(cachefname))
  }
  covs <- read.table(fname, header = TRUE, sep='')
  # Sex as logical
  covs$Male <- as.logical(covs$sex_f31_0_0 == 'Male')
  covs$sex_f31_0_0 <- NULL
  # Set genotype batch and assessment centers as factor
  covs$genotype_measurement_batch_f22000_0_0 <- as.factor(covs$genotype_measurement_batch_f22000_0_0)
  covs$uk_biobank_assessment_centre_f54_0_0 <- as.factor(covs$uk_biobank_assessment_centre_f54_0_0)
  # Remove PCs > 10
  # covs[,paste0('genetic_principal_components_f22009_0_', 11:40)] <- NULL
  # Take median of BMI in case of multiple measurements
  covs$BMI <- apply(covs[,grep('body_mass_index_bmi_f21001_', colnames(covs))], 
                    1, function(x) median(x, na.rm = TRUE))
  covs[,grep('body_mass_index_bmi_f21001_', colnames(covs))] <- NULL
  saveRDS(covs, file=cachefname)
  covs
}

#' Regress
#' 
#' Take always 'Pheno' as outcome, and first 10 PCs, age and Sex as variables. Other should
#' be explicitly given.
#' @examples 
#' regress(m, c('PRS_ukbb_even_chrs', 'PRS_ukbb_odd_chrs'))
regress <- function(m, vars=c(), interaction=c(), withBMI=TRUE, nPCs=10, outcome='Pheno') {
    pcs=paste("genetic_principal_components_f22009_0_",1:nPCs,sep="")
  if(withBMI) {
    bmi <- 'BMI'
  } else {
    bmi <- c()
  }
  ageAndSex <- c('age_when_attended_assessment_centre_f21003_0_0', 'Male')
  f <- paste(outcome, "~", paste(c(interaction, vars, ageAndSex, bmi, pcs), collapse="+"))
  #print(f)
  fit <- lm(as.formula(f), data=m)
}

#' Given a lm, returns R^2, Adjusted R^2 and AIC
performance <- function(y, fit) {
  s <- summary(fit)
  data.frame(Rsqrd=s$r.squared, AdjustedRsqrd=s$adj.r.squared,
             AIC=AIC(fit))
}


# Given a phenotype name and a method, returns PRS and covariates data.frame
get_data <- function(phenotype, method, tissue=NA) {
  prs <- getPRS(phenotype, method, tissue = tissue)
  if(!is.na(tissue)) {
    prs_all <- getPRS(phenotype, method)
    prs <- merge(prs, prs_all)
    colnames(prs)[colnames(prs)=='PRS_even'] <- 'PRS_even_tissue'
    colnames(prs)[colnames(prs)=='PRS_odd'] <- 'PRS_odd_tissue'
  }
  covs <- getCovars(DEFAULT_COV_FILE)
  m <- merge(prs, covs)
  m <- m[complete.cases(m),]
  m[,c('FID', 'IID')] <- NULL
  m
}

#' likelihood ratio test
logLikPVal <- function(f1, f2) {
  ll1 = logLik(f1)
  ll2 = logLik(f2)
  n = 2*(ll2-ll1)
  df = attr(ll2,'df') - attr(ll1,'df')
  1-pchisq(n,df)
}


#' Get variance explained
#' 
#' @param m data.frame with data
#' @param feature column name (e.g., 'PRS_ukbb_odd_chrs')
#' @return Percent varianced explained
single_var_exp <- function(m, phenotype, feature, indices=NULL, withBMI=TRUE) {
  if (! all(is.null(indices))) {
    m <- m[indices,]
  }
  v <- var(m[,feature])
  r <- regress(m, vars=c(feature), withBMI=withBMI)
  s <- summary(r)
  if(! feature %in% rownames(s$coefficients)) {return(NA)}
  b <- s$coefficients[feature,1]
  b^2*v / var(m$Pheno)*100
}
