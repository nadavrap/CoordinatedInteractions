source( '~/coordination_fresh/universal_setup.R' )
source( '~/coordination_fresh/eo/setup.R' )

get_pheno_covar_wtiss <- function(phen, tiss, pval,remove_outliers=TRUE,qnorm=TRUE,internal,seed=0,allchrs=FALSE,return_ID=FALSE,return_bg=FALSE){ 
	print( tiss ) 
	if( !internal ){
		prsfile0	<- paste0(ukbdir, phen, '/external_norm_score_imp_P'	, pval, '_', tiss,'_40PCs.csv')
	} else {
		prsfile0	<- paste0(ukbdir, phen, '/norm_score_P'	, pval, '_', tiss,'_40PCs.csv')
	}
	if( !file.exists( prsfile0 ) ) stop(paste0( 'bad prsfile0:', prsfile0 ))
	prs0	<- read.csv(prsfile0 )
	prs0	<- prs0[,-3]
	colnames(prs0)	<- c( 'FID', 'IID', paste0( tiss, '_', 1:22 ) ) 

	All_tiss<- rowSums( prs0[,paste0(tiss,'_',1:22)], na.rm=T )
	if( allchrs ){
		prs0	<- cbind( prs0, All_tiss )
		###new, not included in >99% of tiss runs but needed in tissXtiss
		##no longer needed in tissXtiss #colnames(prs0)[ncol(prs0)]	<- paste0( tiss, '_All' )
	} else {
		stop() 
	} 

	prs	<- get_pheno_and_covar(phen, pval, remove_outliers, qnorm, internal, seed=seed, return_ID=TRUE, allchrs=allchrs)
	prs	<- merge( prs, prs0 )
	rm( prs0 ) 

	if( return_bg ){ 
		prs0	<- get_bg_tprs( phen, pval, internal )
		prs	<- merge( prs, prs0 )
		rm( prs0 ) 
	}

	prs		<- prs[complete.cases(prs),]
	if( ! return_ID )
		prs[,c('FID', 'IID')] <- NULL

	prs
}

get_bg_tprs	<- function( phen, pval, internal ){
 	bgfile	<- paste0( '~/coordination_fresh/tiss/bg_tprs/', phen, '_'	, pval, '_40PCs.Rdata')
  if(file.exists(bgfile)) return(readRDS(bgfile)) 

	for( tiss in tissues ){
		tprs	<- get_pheno_covar_wtiss(phen, tiss, pval, internal=internal,seed=0,allchrs=TRUE,return_ID=TRUE,return_bg=FALSE)
		tprs	<- tprs[,c('FID','IID',paste0(tiss,'_',1:22))]
		stopifnot( mean( is.na( tprs ) ) == 0 )
		if( tiss == tissues[1]  ){
			tprs0	<- tprs
		} else {
			tprs0	<- merge( tprs, tprs0 )
		}
		rm( tprs )
	}
	stopifnot( mean( is.na( tprs0 ) ) == 0 )
	stopifnot( ncol( tprs0 ) == length(tissues)*22 + 2 )
	tprs0		<- cbind( 
		tprs0[,c('FID','IID')], 
		#sapply( tissues, function(tiss) rowSums( tprs0[,paste0(tiss, '_', 1:22	)] ) )
		sapply( 1:22, function(j) rowSums( tprs0[,paste0(tissues, '_', j	)] ) )
	)
	colnames( tprs0 )	<- c( 'FID', 'IID', paste0( 'BG_', 1:22) )
  saveRDS(tprs0, file=bgfile)
	tprs0
}
