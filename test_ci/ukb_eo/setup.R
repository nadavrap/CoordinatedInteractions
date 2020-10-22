library(foreach)
library(doParallel)
source( '../universal_setup.R' )

ukbdir	<- '/zaitlen/netapp/group/UKBIOBANK/'

get_pheno_and_covar <- function(phenotype, pvalue,remove_outliers=TRUE,qnorm=TRUE,eo_pc=FALSE,internal_prs_40,seed=0,return_ID=FALSE){
	if( internal_prs_40 ){
		prsfile		<- paste0(ukbdir, phenotype, '/norm_score_P'							, pvalue,'_40PCs.csv')
	} else {
		prsfile		<- paste0(ukbdir, phenotype, '/external_norm_score_imp_P'	, pvalue,'.csv')
	}
	if( !file.exists( prsfile  ) ) stop(paste0( 'bad prsfile: ', prsfile ))

	prs		<- read.csv(prsfile	)

	seed	<- as.numeric( seed )
	if( seed == 0 ){
		print( 'Running E/O' )
		chrsE	<- 1:11*2
		chrsO	<- 1:11*2-1
	} else {
		print( 'Running random Chr.' )
		set.seed( seed )
		chrsE	<- sort.list( runif(22) )[1:11]
		chrsO	<- setdiff( 1:22, chrsE )
	}

	nasub	<- which( prs[,'Pheno'] == -9 ) 
	if( phenotype %in% binphens ){
		if( length( nasub ) > 0 )
			prs	<- prs[-nasub,]
	} else {
		if( length( nasub ) > 0 )
			stop('-9 in qphen')
	}

	All			<- apply( prs[,paste0('X'	,c(chrsE,chrsO)	)], 1, sum )
	Even		<- apply( prs[,paste0('X'	,chrsE					)], 1, sum )
	Odd			<- apply( prs[,paste0('X'	,chrsO					)], 1, sum )

	prs	<- cbind( prs[,c('IID','FID','Pheno')], All, Even, Odd )
	covs<- getCovars(paste0(ukbdir,'CSV/covariates.cov'))
	m		<- merge(prs, covs)
	if( eo_pc ){
		epcs	<- read.table( 'flashpca_ukb_even/eigenvectors.txt' , head=T )
		opcs	<- read.table( 'flashpca_ukb_odd/eigenvectors.txt' , head=T )
		colnames(epcs)	<- c( 'IID', 'FID', paste0( 'EvenPC_', 1:40 ) )
		colnames(opcs)	<- c( 'IID', 'FID', paste0( 'OddPC_', 1:40 ) )
		m			<- merge(m,epcs)
		m			<- merge(m,opcs)
		rm( epcs, opcs )
	}
	m			<- m[complete.cases(m),]
	if( ! return_ID )
	m[,c('FID', 'IID')] <- NULL

	if( !phenotype %in% binphens ){
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
