#!Rscript

# qsub -cwd -N PRS -j y -l h_rt=10:00:00 -pe smp 2- ~/bin/R_wrapper.sh ./scripts/prs_summary_even_odd.R
# qsub -cwd -N PRS -j y -l h_rt=10:00:00 -pe shared 2- ~/bin/R_wrapper.sh ./scripts/prs_summary_even_odd.R

#pvalue <- 0.01
PVALUES <- c('0.0001', 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, '1.0', '1e-5', '1e-6', '1e-07', '1e-08')
library(foreach)
library(doParallel)
source('regression_functions.R')

#BOOTSRAP_NUM=2
NSLOTS<-Sys.getenv('NSLOTS')
if(NSLOTS == '') {
    NSLOTS = detectCores()
} else {
  NSLOTS = as.integer(NSLOTS)
}

EXTERNAL_PHENOTYPES <- c('LDL_direct', 'cov_EDU_YEARS', 'T2D', 'Triglycerides', 'BMI', 'Height',
                            'disease_ASTHMA_DIAGNOSED', 'disease_CARDIOVASCULAR')

# Generated like that:
# set.seed(2019)
# L=do.call(c, lapply(0:9, function(i) {l = sample(1:22); l = list(sort(l[1:11]),sort(l[12:22])); names(l) <- paste(paste0('rand', i, '.'), 1:2, sep=''); l}))
# SPLITS = c(list(Odd=seq(2, 22, by=2), Even=seq(1, 22, by=2)), L)
# dput(SPLITS)
SPLITS=list(Even = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22), Odd = c(1,3, 5, 7, 9, 11, 13, 15, 17, 19, 21),
 	rand0.1 = c(1, 2, 7,8, 10, 12, 14, 15, 17, 18, 22),  rand0.2 = c(3, 4, 5, 6, 9, 11, 13, 16, 19, 20, 21),
    	rand1.1 = c(1, 3, 4, 7, 10, 11, 14, 16, 17, 19, 21), rand1.2 = c(2, 5, 6, 8, 9, 12, 13, 15, 18, 20, 22),
	rand2.1 = c(1, 2, 4, 5, 7, 11, 16, 18, 19, 20, 22),  rand2.2 = c(3, 6, 8, 9, 10, 12, 13, 14, 15, 17, 21),
    	rand3.1 = c(2, 3, 6, 7, 9, 12, 13, 17, 20, 21, 22),  rand3.2 = c(1, 4, 5, 8, 10, 11, 14, 15, 16, 18, 19),
    	rand4.1 = c(4, 5, 6, 7, 9, 12, 13, 14, 15, 17, 18),  rand4.2 = c(1, 2, 3, 8, 10, 11, 16, 19, 20, 21, 22),
    	rand5.1 = c(1, 2, 5, 6, 7, 9, 11, 12, 15, 17, 21),   rand5.2 = c(3, 4, 8, 10, 13, 14, 16, 18, 19, 20, 22),
    	rand6.1 = c(3, 6, 7, 8, 9, 11, 13, 19, 20, 21, 22),  rand6.2 = c(1, 2, 4, 5, 10, 12, 14, 15, 16, 17, 18),
    	rand7.1 = c(3, 6, 7, 8, 9, 10, 11, 14, 17, 19, 21),  rand7.2 = c(1, 2, 4, 5, 12, 13, 15, 16, 18, 20, 22),
    	rand8.1 = c(3, 4, 5, 8, 9, 12, 13, 14, 16, 18, 21),  rand8.2 = c(1, 2, 6, 7, 10, 11, 15, 17, 19, 20, 22),
    	rand9.1 = c(1, 4, 5, 6, 7, 10, 13, 15, 16, 19, 20),  rand9.2 = c(2, 3, 8, 9, 11, 12, 14, 17, 18, 21, 22))

#' Get the PRS from file
getPRS <- function(phenotype, method, pvalue, external=FALSE, impORpgen='pgen', nPCs=10) {
    #basename <- paste0('Summary_', tolower(method), '_prsice_', phenotype, '.txt')
    #basename <- paste0(phenotype, '_PRS.txt')
    if (external) {
        basename <- paste0('external_norm_score_', impORpgen, '_P', format(pvalue, scientific=FALSE), '.csv')
    } else {
        npcs <- ifelse(nPCs == 40, '_40PCs', '')
        basename <- paste0('norm_score_P', format(pvalue, scientific=FALSE), npcs, '.csv')
        #print(paste('Reading', basename, 'for phenotype', phenotype))
    }
    fname <- file.path(DATA_DIR, phenotype, basename)
    if(!file.exists(fname)) {return(NA)}
  prs <- read.table(fname, header = TRUE, sep=',', fill=TRUE, stringsAsFactors = FALSE)
  #prs[,3] <- as.numeric(prs[,3])
  # Replace phenotype name with a generic column name
  # colnames(prs)[3] <- 'Pheno'
  prs
}

#' Given a phenotype name and a method, returns PRS and covariates data.frame.
#' Cache the results for next run.
get_data <- function(phenotype, method, pvalue, external=FALSE, impORpgen='pgen', nPCs=10) {
  # use format to avoid scientific notation
    cachefname = paste0('Cache/', phenotype, method, format(pvalue, scientific=FALSE),
                        ifelse(external, paste0('_ext_', impORpgen), ''), '_', nPCs, 'PCs.rds')
  if(file.exists(cachefname)) {
      return(readRDS(cachefname))
  }
  prs <- getPRS(phenotype, method, pvalue=pvalue, external, impORpgen, nPCs)
  covs <- getCovars(DEFAULT_COV_FILE)
  m <- merge(prs, covs)
  m <- m[complete.cases(m),]
  m[,c('FID', 'IID')] <- NULL
  if (sum(sapply(m, class)=='factor')!=2) {
     stop('Check classes of columns, cannot find 2 factors: uk_biobank_assessment_centre_f54_0_0 and genotype_measurement_batch_f22000_0_0')
  }
  saveRDS(m, cachefname)
  m
}

sum_prs <- function(m, chromosomes) {
    rows <- paste0('X', chromosomes)
    apply(m[,rows], 1, sum)
}

auc <- function(phenos, prediction) {
    if(nlevels(as.factor(phenos)) != 2) {
        return(NA)
    } else {
        as.numeric(pROC::auc(pROC::roc(phenos, prediction)))
    }
}

table2_row <- function(phenotype, method='Plink', int1='Odd', int2='Even', pvalue, normalize_PRS=FALSE, remove_outliers=FALSE, external=FALSE, impORpgen='', controlForBMI=TRUE, nPCs=20, adjustForNpcs=40) {
  withBMI=phenotype != 'BMI' & controlForBMI
  #for each phenotype using blood_EOSINOPHIL_COUNT as an example
  m <- get_data(phenotype, method, pvalue, external, impORpgen, nPCs=nPCs)
  if(nrow(m) == 0) {
    return(rep(NA, 6))
  }
  if (phenotype %in% BINARY_TRAITS){
      is_binary <- TRUE
      remove_outliers <- FALSE
      if (!all(levels(as.factor(m$Pheno)) %in% c('-9', '1', '2')) & !all(levels(as.factor(m$Pheno)) %in% c('1', '2'))) {
          warning('Does not looks like a binary trait')
          return(NA)
      }
      m <- droplevels(m[m$Pheno >= 0,])
      m$OrgPheno <- m$Pheno - 1
  }
  if (remove_outliers) {
      outlier_values <- boxplot.stats(m$Pheno)$ou
      m <- m[!m$Pheno %in% outlier_values,]
  }
  # Normalize phenotype
  m$Pheno <- preprocessCore::normalize.quantiles.use.target(as.matrix(m$Pheno), target=rnorm(nrow(m)))[,1]
  #m$Pheno <- quantile_normalisation(as.matrix(m$Pheno))

  if(! 'All' %in% colnames(m)) {
         m[,'All'] <- sum_prs(m, 1:22)
  }
  m[,int1] <- sum_prs(m, SPLITS[[int1]])
  m[,int2] <- sum_prs(m, SPLITS[[int2]])

  if (normalize_PRS){
      m[,'All'] <- scale(m[,'All'], center=TRUE, scale=TRUE)[,1]
      m[,int1] <- scale(m[,int1], center=TRUE, scale=TRUE)[,1]
      m[,int2] <- scale(m[,int2], center=TRUE, scale=TRUE)[,1]
  }
  var_exp <- function(m, indices=NA, int1, int2) {
    m <- m[indices,]
    # Covariance matrix
    v <- var(cbind(m[,int1], m[,int2]))
    v1 = v[1,1] # PRS_ukbb_odd_chrs variance 
    v2 = v[2,2] # PRS_ukbb_even_chrs variance
    cv = v[1,2] # Covariance
    # Regress Pheno ~ PRS_ukbb_even_chrs + PRS_ukbb_odd_chrs + covariates
    res1 <- regress(m, vars=c(int1, int2), withBMI=withBMI, nPCs=adjustForNpcs)
    s <- summary(res1)
    b1 = s$coefficients[int1,1] # PRS_ukbb_even_chrs' beta
    b2 = s$coefficients[int2,1] # PRS_ukbb_odd_chrs' beta
    #this is the nubmer to report in %VE
    pve1 = (b1^2*v1 + b2^2*v2 +2*cv*b1*b2)/var(m$Pheno)*100
  }
  
  # bootstrapping with BOOTSRAP_NUM replications
  pve_all <- single_var_exp(m, phenotype, 'All', withBMI=withBMI)
  pve_1 <- single_var_exp(m, phenotype, int1, withBMI=withBMI)
  pve_2 <- single_var_exp(m, phenotype, int2, withBMI=withBMI)
  auc_all <- auc(m$Pheno, m$All)
  #results <- boot::boot(data=m, var_exp, R=BOOTSRAP_NUM, int1=int1, int2=int2)
  #pve <- results$t0
  #pve_stderr <- sd(results$t)
  
  PRS_cor <- function(m, indices, int1, int2, adjustForNpcs) {
    m <- m[indices,]
    res1 <- summary(regress(m, outcome=int1, withBMI=withBMI, nPCs=adjustForNpcs))$residuals
    res2 <- summary(regress(m, outcome=int2, withBMI=withBMI, nPCs=adjustForNpcs))$residuals
    cor(res1,res2)
  }
  #results <- boot::boot(data=m, PRS_cor, R=BOOTSRAP_NUM, int1=int1, int2=int2)
  #cor2 <- results$t0
  #cor2_stderr <- sd(results$t)
  vars_cor <- PRS_cor(m, 1:nrow(m), int1, int2, adjustForNpcs=adjustForNpcs)

  #Column three is the same, just put the value into the word doc. 
  interaction <- regress(m, vars=c(int1, int2),
                           interaction=c(paste(int1, '*', int2)), withBMI=withBMI, nPCs=adjustForNpcs)
  s <- summary(interaction)
  beta_inter <- s$coefficients[paste(int1, int2, sep = ':'),1]
  #p3 <- lmp(interaction)
  pval_inter <- s$coefficients[paste(int1, int2, sep = ':'),4]
  
  #return(c(pve_all, pve_1, pve_2, pve, pve_stderr, cor2, cor2_stderr, beta_inter, -log10(pval_inter)))
  return(c(pve_all, pve_1, pve_2, auc_all, vars_cor, beta_inter, -log10(pval_inter)))
}
analyze_even_odd <- function(pvalue, normalize_PRS=FALSE, remove_outliers=FALSE, phenotypes=PHENOTYPES, external=FALSE, impORpgen='', controlForBMI=TRUE, nPCs=10, adjustForNpcs=40) {
    registerDoParallel(min(length(phenotypes), NSLOTS)) 
    d_even_odd <- foreach (phenotype=phenotypes, .combine=rbind) %dopar% {
      table2_row(phenotype, pvalue=pvalue, normalize_PRS=normalize_PRS, remove_outliers=remove_outliers, external=external, impORpgen=impORpgen, controlForBMI=controlForBMI, nPCs=nPCs, adjustForNpcs=adjustForNpcs)
    }
    d_rands <- lapply(0:9, function(i) {
        d_rand <- foreach (phenotype=phenotypes, .combine=rbind) %dopar% {
            table2_row(phenotype, int1 = paste0('rand', i, '.1'), int2=paste0('rand', i, '.2'), pvalue=pvalue,
                       normalize_PRS=normalize_PRS, remove_outliers=remove_outliers, external=external, impORpgen=impORpgen, controlForBMI=controlForBMI, nPCs=nPCs, adjustForNpcs=adjustForNpcs)
        }
    })
    
    stopImplicitCluster()
    d_rands <- do.call(rbind, d_rands)
    d <- rbind(d_even_odd, d_rands)
    #rownames(d) <- phenotypes
    #d <- rbind(d_even_odd, d_rand)
    colnames(d) <- c('%VE_all', '%VE_1', '%VE_2', 'AUC', 'Cor', 'Beta', 'mLog10PVal')
    if (external) {
        d <- cbind(Pheno=phenotypes, Set=c(rep(paste0('EvenOdd_ext', impORpgen), nrow(d_even_odd)),
                                            rep(paste0('RandomChr_ext', impORpgen), nrow(d_rands))), d)
    } else {
        d <- cbind(Pheno=phenotypes, Set=c(rep('EvenOdd', nrow(d_even_odd)), rep('RandomChr', nrow(d_rands))), d)
    }
    #d
    is_normalized <- ifelse(normalize_PRS, 'normalized_', '')
    is_outliers_remove <- ifelse(remove_outliers, 'outliers_removed_','')
    is_external <- ifelse(external, paste0('external_', impORpgen, '_'), '')
    is_BMI <- ifelse(controlForBMI, '', 'woCovBMI_')
    is_PCs <- paste0(nPCs, 'PCs_AdjustFor', adjustForNpcs)
    write.csv(d, file=paste0('Interactions/Even_odd_rand_P', format(pvalue, scientific=FALSE),'_PRS_',is_normalized,
                                is_outliers_remove, is_external, is_BMI, is_PCs, '_', format(Sys.Date(), "%Y_%m_%d"),
                                '.csv'), row.names=FALSE)
    d
}

examles <- function() {
# for(pvalue in c('0.0001', 0.001, 0.01, 0.05, 0.1, '1.0')) {print(Sys.time()); analyze_even_odd(pvalue, normalize_PRS=TRUE)}
# for(pvalue in c('0.0001', 0.001, 0.01, 0.05, 0.1, '1.0')) {print(Sys.time()); analyze_even_odd(pvalue)}
# for(pvalue in c('0.0001', 0.001, 0.01, 0.05, 0.1, '1.0')) {print(Sys.time()); analyze_even_odd(pvalue, remove_outliers=TRUE)}
# pvalue = '0.01'; analyze_even_odd(pvalue, remove_outliers=TRUE, normalize_PRS=TRUE, phenotypes=EXTERNAL_PHENOTYPES, external=TRUE, impORpgen='pgen')
# pvalue = '0.01'; analyze_even_odd(pvalue, remove_outliers=TRUE, normalize_PRS=TRUE, controlForBMI=FALSE)
# for (nPCs in c(10,20, 40)) {for(pvalue in PVALUES) {print(Sys.time()); analyze_even_odd(pvalue, normalize_PRS=TRUE, remove_outliers=TRUE, nPCs=nPCs)}}
# nPCs=40; for(pvalue in PVALUES) {print(Sys.time()); analyze_even_odd(pvalue, normalize_PRS=TRUE, remove_outliers=TRUE, nPCs=nPCs)}

    pvals <- c('1e-07','1e-06','1e-05','0.0001','0.001','0.005','0.01','0.05','0.1','0.2','0.5','1.0')
    d <- data.frame(do.call(rbind, lapply(pvals, function(pvalue)
                                            table2_row('T2D', pvalue=pvalue, normalize_PRS=TRUE, remove_outliers=FALSE,
                                                       external=FALSE, impORpgen='', controlForBMI=TRUE, nPCs=10))))
    d$P = as.numeric(pvals); colnames(d) <- c('%VE_all', '%VE_1', '%VE_2', 'AUC', 'Cor', 'Beta', 'mLog10PVal', 'P-threshold')

    # Plot case ratio by PRS quantile
    quantiles <- quantile(m$All, probs = seq(0, 1, 0.1))
    do.call(rbind, lapply(2:length(quantiles), function(i) {
        M <- m[m$All >= quantiles[i-1] & m$All < quantiles[i],]; c(names(quantiles)[i], sum(M$Pheno == 2)/nrow(M))}))
    caseRatio <- function(pvalue, phenotype='T2D', nquantiles=10) {
        m <- get_data(phenotype, method, pvalue, external, impORpgen)
        m <- droplevels(m[m$Pheno >= 0,]); m[,'All'] <- sum_prs(m, 1:22)
        quantiles <- quantile(m$All, probs = seq(0, 1, 1/nquantiles));
        M1 <- m[m$All >= quantiles[1] & m$All < quantiles[2],]
        M10 <- m[m$All > quantiles[nquantiles] & m$All <= quantiles[nquantiles+1],]
        sum(M10$Pheno == 2)/sum(M1$Pheno == 2)
    }
    do.call(rbind, lapply(c('1e-07', '1e-06', '1e-05', PVALUES), function(p) c(p, caseRatio(pvalue=p))))
    do.call(rbind, lapply(PVALUES, function(p) c(p, caseRatio(phenotype='T2D', pvalue=p, nquantiles=40))))

    # External datasets
    for(pvalue in PVALUES){
        analyze_even_odd(pvalue, remove_outliers=TRUE, normalize_PRS=TRUE, phenotypes=EXTERNAL_PHENOTYPES,
                            external=TRUE, impORpgen='pgen')
        analyze_even_odd(pvalue, remove_outliers=TRUE, normalize_PRS=TRUE, phenotypes=EXTERNAL_PHENOTYPES,
                            external=TRUE, impORpgen='imp')
    }
}

#' Print for each phenotype the results for which the %Variance explained is maximized.
get_best_pve <- function(nPCs=10, date='2019_08_19') {
    #d <- lapply(PVALUES, function(p) {d <- read.csv(paste0('Interactions/', 'Even_odd_rand_P',p,'_PRS_normalized_outliers_removed_external_imp_',nPCs,'PCs_2019_08_12.csv')); d$Pval <- p; d})
    d <- lapply(PVALUES, function(p) {d <- read.csv(paste0('Interactions/', 'Even_odd_rand_P',p,
                                                            '_PRS_normalized_outliers_removed_',nPCs,'PCs_AdjustFor40',
                                                            date,'.csv'));
                                        d$Pval <- p;d})
    d <- data.frame(do.call(rbind,d))
    d <- d[rev(order(d$X.VE_all)),]
    pheno_pval <- d[!duplicated(d$Pheno),c('Pheno', 'Pval')]
    d <- merge(d, pheno_pval) # Keep pval that maximize PVE_all
    write.csv(d, file=paste0('Interactions/Even_odd_rand_maximizing_PVE', '_', format(Sys.Date(), "%Y_%m_%d"),'.csv'))
}
