rm( list=ls() )
load( '../ukb_eo/Rdata/bestps.Rdata' )
bestps['ext','disease_ASTHMA_DIAGNOSED']	<- '0.0001'

source('../functions.R')
source('setup.R')

#seeds	<- c( 0, round( 7 + (2:19)^3/42 ), 192 )
S			<- length(seeds)

for( type in c( 'bin10', 'int40', 'ext' ) )
	for( seed in seeds )
			for( phen in sample(phens) )
{
	if( (type == 'ext'	& !phen %in% extphens) |
			(type == 'bin10'& !phen %in% binphens) |
			(type == 'int40'&  phen %in% binphens) ) next

	savefile	<- paste0('Rdata/Tissue_evenodd_internalPRS_2019_08_21_', bestps[type,phen], '_', phen, '_', type, '_', seed, '.Rdata')
	sinkfile	<- paste0('Rout/Tissue_evenodd_internalPRS_2019_08_21_',  bestps[type,phen], '_', phen, '_', type, '_', seed, '.Rout')
	if( file.exists( savefile ) | file.exists( sinkfile ) ) next
	sink( sinkfile )
	#withBMI	<- (phen != 'BMI')
	withBMI	<- (! phen %in% c( 'BMI', 'T2D', 'disease_CARDIOVASCULAR','disease_HI_CHOL_SELF_REP'))
	nPC			<- ifelse( type == 'int40', 40, 10 )

	m	<- get_pheno_covar_alltiss(phen, TISSUES, bestps[type,phen], internal_prs_40=( type %in% c( 'int40', 'bin10' ) ), seed=seed )

	dat	<- lapply( TISSUES, function(tiss1){
		registerDoParallel(5)
		datloc	<- foreach (tiss2=TISSUES, .combine = rbind) %dopar% {
			inte		<- paste0( tiss1, '_Even' )
			into		<- paste0( tiss2, '_Odd'  )
			intobg	<- 'BG_Odd'
			intebg	<- 'BG_Even'
			s <- summary(regress(m, withBMI=withBMI, nPC=nPC, binary=( phen %in% binphens),
				vars=c(
					'Odd', 'Even', into, inte, intobg, intebg,
					'Odd * Even',  paste0( 'Odd*', inte ), paste0( 'Odd*', intebg ), paste0( 'Even*', into ), paste0( 'Even*', intobg ), paste0( intebg, '*', intobg )
				), interaction=paste0(inte, '*', into)))
			beta_inter <- s$coefficients[paste0(inte, ':', into),1]
			pval_inter <- s$coefficients[paste0(inte, ':', into),4]
			c(beta_inter, -log10(pval_inter))
		}
		stopImplicitCluster()
		rownames(datloc) <- NULL
		colnames(datloc) <- c('beta', 'mLog10Pval')
		cbind(Tissue=TISSUES, datloc)
	})
	names( dat )	<- TISSUES
	save( dat, file=savefile )
	sink()
}
