library(foreach)
library(doParallel)
source( '../universal_setup.R' )
source( '../ukb_eo/setup.R' )

get_pheno_covar_alltiss <- function(phen, TISSUES, pvalue,remove_outliers=TRUE,qnorm=TRUE,eo_pc=FALSE,internal_prs_40,seed=0){
	prs0	<- get_pheno_and_covar(phen, pvalue, internal_prs_40=internal_prs_40, seed=seed, return_ID=TRUE, remove_outliers=remove_outliers,qnorm=qnorm,eo_pc=eo_pc )

	seed	<- as.numeric( seed )
	if( seed == 0 ){
		chrsE	<- 1:11*2
		chrsO	<- 1:11*2-1
	} else {
		set.seed( seed )
		chrsE	<- sort.list( runif(22) )[1:11]
		chrsO	<- setdiff( 1:22, chrsE )
	}

	#### avg ts-prs
	load( paste0('../ukb_tiss/Rdata/bg_TS_PRS/', pvalue, '_', phen, '_', type,  '.Rdata') )
	if( !phen %in% binphens ){
		prs$Pheno	<- NULL
	} else {
		nasub	<- which( prs[,'Pheno'] == -9 ) 
		if( length( nasub ) > 0 )
			prs	<- prs[-nasub,]
	}

	BG_All	<- rowSums( prs[,paste0( rep(TISSUES,each=22), '_', rep( c(chrsE,chrsO)	, ntiss ))] )/ntiss
	BG_Even	<- rowSums( prs[,paste0( rep(TISSUES,each=22), '_', rep( chrsE					, ntiss ))] )/ntiss
	BG_Odd	<- rowSums( prs[,paste0( rep(TISSUES,each=22), '_', rep( chrsO					, ntiss ))] )/ntiss


	TS_Even	<- sapply( TISSUES, function(tiss) rowSums( prs[,paste0(tiss, '_'	,chrsE	)] ) )
	TS_Odd	<- sapply( TISSUES, function(tiss) rowSums( prs[,paste0(tiss, '_'	,chrsO	)] ) )

	prs						<- cbind( prs[,c('IID','FID')], TS_Even, TS_Odd, BG_All, BG_Even, BG_Odd )
	colnames(prs)	<- c( 'IID', 'FID', paste0( TISSUES, '_Even' ), paste0( TISSUES, '_Odd' ), paste0( 'BG_', c( 'All', 'Even', 'Odd' ) ) )


	prs	<- merge( prs, prs0 )
	rm( prs0 )

	prs			<- prs[complete.cases(prs),]
	prs[,c('FID', 'IID')] <- NULL
	if( !phen %in% binphens )
		prs$Pheno <- scale(prs$Pheno)
	prs
}
