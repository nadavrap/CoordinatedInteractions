rm( list=ls() )
load( '../ukb_eo/Rdata/bestps.Rdata' )
bestps['ext','disease_ASTHMA_DIAGNOSED']	<- '0.0001'
source('../functions.R')
source('setup.R')

for( type in c( 'ext', 'bin10', 'int40' ) )
	for( seed in seeds )
		for( phen in sample(phens) )
{
	if( (type == 'ext'	& !phen %in% extphens) |
			(type == 'bin10'& !phen %in% binphens) |
			(type == 'int40'&  phen %in% binphens) ) next

	savefile	<- paste0('Rdata/Tissue_', bestps[type,phen], '_', phen, '_', type, '_', seed, '_nobg.Rdata')
	sinkfile	<- paste0( 'Rout/Tissue_', bestps[type,phen], '_', phen, '_', type, '_', seed, '_nobg.Rout' )
	if( file.exists( savefile ) | file.exists( sinkfile ) ) next
	sink( sinkfile )
	withBMI	<- (! phen %in% c( 'BMI', 'T2D', 'disease_CARDIOVASCULAR','disease_HI_CHOL_SELF_REP'))
	nPC			<- ifelse( type == 'int40', 40, 10 )
	binary	<- phen %in% binphens

try({
	registerDoParallel(4)
	dat <- foreach (i=TISSUES, .combine = rbind) %dopar% {
		m <- get_pheno_covar_tiss(phen, tissue=i, bestps[type,phen], internal_prs_40=( type %in% c( 'int40', 'bin10' ) ), seed=seed )
		print( i )

		inta <- paste0( i, '_All'  )
		into <- paste0( i, '_Odd'  )
		inte <- paste0( i, '_Even' )

		pve_all		<- single_var_exp(m, withBMI=withBMI, nPC=nPC, phenotype='Pheno', binary=binary, feature='All'		)
		pve_tiss	<- single_var_exp(m, withBMI=withBMI, nPC=nPC, phenotype='Pheno', binary=binary, feature=inta			)
		pve_1			<- single_var_exp(m, withBMI=withBMI, nPC=nPC, phenotype='Pheno', binary=binary, feature=into			)
		pve_2			<- single_var_exp(m, withBMI=withBMI, nPC=nPC, phenotype='Pheno', binary=binary, feature=inte			)

		res1 <- summary(regress(		m, withBMI=withBMI, nPC=nPC, outcome=into			, binary=FALSE ))$residuals
		res2 <- summary(regress(		m, withBMI=withBMI, nPC=nPC, outcome=inte			, binary=FALSE ))$residuals
		vars_cor <- cor(res1,res2)
		
		s1		<- summary(regress(		m, withBMI=withBMI, nPC=nPC, outcome='Pheno'	, binary=binary,
			vars=c( inte,	'Odd', 'Even','Odd * Even'), interaction=paste(inte, '* Odd')	))

		s2		<- summary(regress(		m, withBMI=withBMI, nPC=nPC, outcome='Pheno'	, binary=binary,
			vars=c( into,	'Odd', 'Even','Odd * Even'), interaction=paste(into, '* Even')	))
		
		return(c( i, pve_all, pve_tiss, pve_1, pve_2, vars_cor,
							s1$coefficients[paste(inte, 'Odd' , sep = ':'),1] ,
			-log10(	s1$coefficients[paste(inte, 'Odd' , sep = ':'),4]),
							s2$coefficients[paste(into, 'Even', sep = ':'),1] ,
			-log10(	s2$coefficients[paste(into, 'Even', sep = ':'),4]) ))
	}
	stopImplicitCluster()
	rownames(dat) <- NULL
	colnames(dat) <- c('Tissue','%VE', '%VE_tissue', '%VE_tissue_1', '%VE_tissue_2', 'Cor', 'BTeAo', 'mLog10PvalTeAo', 'BToAe', 'mLog10PvalToAe')
	save( dat, file=savefile )
})
	sink()
}
