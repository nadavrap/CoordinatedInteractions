#!Rscript
source('regression_functions.R')
BOOTSRAP_NUM=2
PVALUE <- 0.01
#PHENOTYPES <- c('BMI', 'Height', 'LDL_direct', 'Glucose', 'Triglycerides')
TISSUES <- c(paste0('Franke.', c(11, 14, 46, 69, 84, 89, 120)),
             c('Liver', 'Adipose', 'Brain', 'Brain_Hippocampus', 'Pancreas', 'Skeletal_Muscle'))
OSX_DATA_DIR <- "/Users/rappoportn/UKBB/data/"
DATA_DIR <- ifelse(dir.exists(OSX_DATA_DIR), OSX_DATA_DIR, '/zaitlen/netapp/group/UKBIOBANK/')
REMOVE_OUTLIERS <- TRUE
library(foreach)
library(doParallel)


get_pheno_and_covar <- function(phenotype, pvalue) {
  normalize_phenotype <- TRUE # Quantile normalization
  withBMI <- phenotype != 'BMI' # Include BMI as covariates
  prs <- read.csv(paste0(DATA_DIR, phenotype, '/tissue_specific_PRS_P', pvalue,'.csv'))
  covs <- getCovars(DEFAULT_COV_FILE)
  m <- merge(prs, covs)
  m <- m[complete.cases(m),]
  m[,c('FID', 'IID')] <- NULL
  # If no data available, raise warning and return a data.frame with NA
  if(nrow(m) == 0) {
    #warning(paste(phenotype, 'has not PRS results using', method))
    return(data.frame(N=NA))
  }
  if (normalize_phenotype) {
      m$Pheno <- preprocessCore::normalize.quantiles.use.target(as.matrix(m$Pheno), target=rnorm(nrow(m)))[,1]
  }
  m
}

table2_row_tissue <- function(phenotype, pvalue=0.01, tissue='Franke.11',m=c(), remove_outliers=FALSE) {
  withBMI=phenotype != 'BMI'
  #for each phenotype using blood_EOSINOPHIL_COUNT as an example
  if (length(m) == 0) {
    m <- get_pheno_and_covar(phenotype, pvalue)
  }
  if(nrow(m) == 0) {
    return(rep(NA, 6))
  }
  if (remove_outliers) {
      outlier_values <- boxplot.stats(m$Pheno)$ou
      m <- m[!m$Pheno %in% outlier_values,]
  }

  into <- paste0('Odd_', tissue)
  inte <- paste0('Even_', tissue)
  inta <- paste0('All_', tissue)
  # Normalize phenotype
  m$Pheno <- preprocessCore::normalize.quantiles.use.target(as.matrix(m$Pheno), target=rnorm(nrow(m)))[,1]
  #m$Pheno <- quantile_normalisation(as.matrix(m$Pheno))
  
  #for column 1 on Y~PRS (or now PRS) do:
   
  # var_exp <- function(m, indices=NA, int1, int2) {
  #   m <- m[indices,]
  #   # Covariance matrix
  #   v <- var(cbind(m[,int1], m[,int2]))
  #   v1 = v[1,1] # PRS_ukbb_odd_chrs variance 
  #   v2 = v[2,2] # PRS_ukbb_even_chrs variance
  #   cv = v[1,2] # Covariance
  #   # Regress Pheno ~ PRS_ukbb_even_chrs + PRS_ukbb_odd_chrs + covariates
  #   res1 <- regress(m, vars=c(int1, int2), withBMI=withBMI)
  #   s <- summary(res1)
  #   b1 = s$coefficients[int1,1] # PRS_ukbb_even_chrs' beta
  #   b2 = s$coefficients[int2,1] # PRS_ukbb_odd_chrs' beta
  #   #this is the nubmer to report in %VE
  #   pve1 = (b1^2*v1 + b2^2*v2 +2*cv*b1*b2)/var(m$Pheno)*100
  # }
  
  # bootstrapping with BOOTSRAP_NUM replications
  pve_all <- single_var_exp(m, phenotype, 'All')
  pve_tissue_all <- single_var_exp(m, phenotype, inta)
  pve_1 <- single_var_exp(m, phenotype, into)
  pve_2 <- single_var_exp(m, phenotype, inte)

  # results <- boot::boot(data=m, var_exp, R=BOOTSRAP_NUM, int1=into, int2=inte)
  # pve <- results$t0
  # pve_stderr <- sd(results$t)
  
  PRS_cor <- function(m, indices, int1, int2) {
    m <- m[indices,]
    res1 <- summary(regress(m, outcome=int1, withBMI=withBMI))$residuals
    res2 <- summary(regress(m, outcome=int2, withBMI=withBMI))$residuals
    cor(res1,res2)
  }
  # results <- boot::boot(data=m, PRS_cor, R=BOOTSRAP_NUM, int1=into, int2=inte)
  # cor2 <- results$t0
  # cor2_stderr <- sd(results$t)
  vars_cor <- PRS_cor(m, 1:nrow(m), into, inte)
  
  #Column three is the same, just put the value into the word doc. 
  # interaction <- regress(m, vars=c(int1, int2),
  #                          interaction=c(paste(int1, '*', int2)), withBMI=withBMI)
  # Pheno~PRS_odd + PRS_tissue_even + PRS_tissue_even*PRS_odd + covariates
  interaction <- regress(m, vars=c(inte, 'Odd'),
                           interaction=c(paste(inte, '* Odd')),
                         withBMI=withBMI)
  s <- summary(interaction)
  beta_interTeAo <- s$coefficients[paste(inte, 'Odd', sep = ':'),1]
  pval_interTeAo <- s$coefficients[paste(inte, 'Odd', sep = ':'),4]
  # Pheno~ PRS_even + PRS_tissue_odd + PRS_tissue_odd*PRS_even+ covariates
  interaction <- regress(m, vars=c(into, 'Even'),
                           interaction=c(paste(into, '* Even')),
                         withBMI=withBMI)
  s <- summary(interaction)
  beta_interToAe <- s$coefficients[paste(into, 'Even', sep = ':'),1]
  pval_interToAe <- s$coefficients[paste(into, 'Even', sep = ':'),4]
  
  #return(c(pve_all, pve_tissue_all, pve_1, pve_2, pve, pve_stderr, cor2, cor2_stderr, beta_interTeAo, -log10(pval_interTeAo), beta_interToAe, -log10(pval_interToAe)))
  return(c(pve_all, pve_tissue_all, pve_1, pve_2, vars_cor, beta_interTeAo, -log10(pval_interTeAo), beta_interToAe, -log10(pval_interToAe)))
}

remove_outliers <- REMOVE_OUTLIERS
registerDoParallel(min(length(PHENOTYPES)*length(TISSUES), detectCores()))
d <- foreach (phenotype=PHENOTYPES, .combine=rbind) %dopar% {
  m <- get_pheno_and_covar(phenotype, PVALUE)
  foreach (i=TISSUES, .combine = rbind) %dopar% {
    table2_row_tissue(phenotype, PVALUE, tissue = i, m=m, remove_outliers=remove_outliers)
  }
}
stopImplicitCluster()

#rownames(d) <- PHENOTYPES
colnames(d) <- c('%VE_all', '%VE_tissue_all', '%VE_1', '%VE_2', 'Cor', 'BTeAo', 'mLog10PvalTeAo', 'BToAe', 'mLog10PvalToAe')
d <- cbind(Pheno=rep(TISSUES, length(PHENOTYPES)), Tissue=rep(PHENOTYPES, each=length(TISSUES)), d)
is_outliers_remove <- ifelse(remove_outliers, 'outliers_removed_','')
write.csv(d, file=paste0('./Even_odd_Tissue_PRS_', is_outliers_remove, format(Sys.Date(), "%Y_%m_%d"),'.csv'), row.names=FALSE, quote = FALSE)
# 
