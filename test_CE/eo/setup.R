#source( '../universal_setup.R' )
source('~/coordination_fresh/universal_setup.R')
nperm	<- 100 

get_pheno_and_covar <- function(phen, pval,remove_outliers=TRUE,qnorm=TRUE,internal,seed=0,return_ID=FALSE,allchrs=FALSE){

	### load PRS and Pheno
	if( internal ){
		prsfile	<- paste0(ukbdir, phen, '/norm_score_P'							, pval,'_40PCs.csv')
	} else {
		prsfile	<- paste0(ukbdir, phen, '/external_norm_score_imp_P', pval,'_40PCs.csv')
	} 
	if( !file.exists( prsfile ) ) stop(paste0( 'bad prsfile: ', prsfile )) 
	prs		<- read.csv(prsfile	)

	### Split EO -- not used in final version with more powerful test
	seed	<- as.numeric( seed )
	if( seed == 0 ){
		chrsE	<- 1:11*2
		chrsO	<- 1:11*2-1
	} else {
		set.seed( seed )
		chrsE	<- sort.list( runif(22) )[1:11]
		chrsO	<- setdiff( 1:22, chrsE )
	}

	nasub	<- which( prs[,'Pheno'] == -9 ) 
	if( phen %in% binphens ){
		if( length( nasub ) > 0 )
			prs	<- prs[-nasub,]
	} else {
		if( length( nasub ) > 0 )
			stop('-9 in qphen')
	}

	All			<- rowSums( prs[,paste0('X'	,c(chrsE,chrsO)	)], na.rm=T ) 
	if( allchrs ){
		prs	<- cbind( prs[,c('IID','FID','Pheno')], All, prs[,paste0('X',1:22)] )
	} else {
		Even		<- rowSums( prs[,paste0('X'	,chrsE					)], na.rm=T )
		Odd			<- rowSums( prs[,paste0('X'	,chrsO					)], na.rm=T )
		prs	<- cbind( prs[,c('IID','FID','Pheno')], All, Even, Odd )
	}

	### merge with covariates
	covs<- getCovars('~/coordination_fresh/final_covars.txt') #based on: paste0(ukbdir,'CSV/covariates_norel3.cov') 
	m		<- merge(prs, covs)
	m		<- m[complete.cases(m),]
	if( ! return_ID )
	m[,c('FID', 'IID')] <- NULL

	if( !phen %in% binphens ){
		if (remove_outliers) {
			outlier_values <- boxplot.stats(m$Pheno)$ou
			m <- m[!m$Pheno %in% outlier_values,]
		}
		if( qnorm )
			m$Pheno <- preprocessCore::normalize.quantiles.use.target(as.matrix(m$Pheno), target=rnorm(nrow(m)))[,1]
		m$Pheno <- scale(m$Pheno)
	}
	m
}
