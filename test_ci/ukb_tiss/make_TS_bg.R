rm( list=ls() )
source('../functions.R')
source( '../universal_setup.R' )
#source('setup.R')

load( '../ukb_eo/Rdata/bestps.Rdata' )
bestps['ext','disease_ASTHMA_DIAGNOSED']	<- '0.0001'

ukbdir	<- '/wynton/scratch/UKBIOBANK/'

for( type in sample(c( 'int40', 'bin10', 'ext' )) )
	for( phen in sample(phens) )
{
	if( (type == 'ext'	& !phen %in% extphens) |
			(type == 'bin10'& !phen %in% binphens) |
			(type == 'int40'&  phen %in% binphens) ) next

	savefile	<- paste0('Rdata/bg_TS_PRS/', bestps[type,phen], '_', phen, '_', type, '.Rdata')
	if( file.exists( savefile	) ) next

	print( phen )
	for( tissue in TISSUES ){
		cat( tissue )
		if( type == 'bin10' ){
			prsfile0	<- paste0(ukbdir, phen, '/norm_score_P'							, bestps[type,phen], '_', tissue, '.csv')
		} else if( type == 'int40' ){
			prsfile0	<- paste0(ukbdir, phen, '/norm_score_P'							, bestps[type,phen], '_', tissue, '_40PCs.csv')
		} else if( type == 'ext' ){
			prsfile0	<- paste0(ukbdir, phen, '/external_norm_score_imp_P', bestps[type,phen], '_', tissue, '.csv')
		}
		if( !file.exists( prsfile0 ) ) stop(paste0( 'bad prsfile0:', prsfile0 ))
		prs0	<- read.csv(prsfile0 )
		colnames(prs0)	<- c( 'FID', 'IID', 'Pheno', paste0( tissue, '_', 1:22 ) )

		if( tissue == TISSUES[1] ){
			prs	<- prs0
		} else {
			prs	<- merge( prs, prs0 )
		}
		rm( prs0 )
	}
	cat( '\n\n' )
	save( prs, file=savefile )
}

for( type in c( 'int40', 'bin10', 'ext' ) )
	for( phen in sample(phens) )
{

	if( (type == 'ext'	& !phen %in% extphens) |
			(type == 'bin10'& !phen %in% binphens) |
			(type == 'int40'&  phen %in% binphens) ) next

	loadfile	<- paste0('Rdata/bg_TS_PRS/', bestps[type,phen], '_', phen, '_', type, '.Rdata')
	savefile	<- paste0('Rdata/bg_TS_PRS/', bestps[type,phen], '_', phen, '_', type, '_avgs.Rdata')
	if( ! file.exists( loadfile	) | file.exists( savefile	) ) next

	load( loadfile )
	BG_All	<- rowMeans( sapply( TISSUES, function(tiss) rowSums( prs[,paste0(tiss, '_', 1:22	)] ) ) )
	TS_All	<-					 sapply( TISSUES, function(tiss) rowSums( prs[,paste0(tiss, '_'	,1:22	)] ) )
	colnames(TS_All)	<- paste0( 'PRS_', TISSUES )

	chrcols	<- which( colnames(prs) %in% c(sapply( TISSUES, function(tiss) paste0(tiss, '_'	,1:22 ) ) ) )
	prs			<- prs[,-chrcols]

	prs	<- cbind( prs, BG_All, TS_All )
	colnames(prs)	<- c( 'FID', 'IID', 'Pheno', 'PRS_BG', paste0( 'PRS_', TISSUES )  )
	save( prs, file=savefile )
}
