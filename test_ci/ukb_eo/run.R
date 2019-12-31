rm( list=ls() )
source('../functions.R')
source('setup.R')

for( type in c( 'bin10', 'ext', 'int40' ) )
	for( phen in sample(phens) )
		for( seed in seeds )
			for( pv in pvs )
{
	if( (type == 'ext'	& !phen %in% extphens) |
			(type == 'bin10'& !phen %in% binphens) |
			(type == 'int40'&  phen %in% binphens) ) next

	if( seed > 0 ){
		if( !file.exists( 'Rdata/bestps.Rdata' ) ) next
		load( 'Rdata/bestps.Rdata' )
		#if( pv != bestps[type,phen] ) next
		if( pv != bestps[type,phen] & pv != pvs[which( pvs == bestps[type,phen] ) - 1] ) next
	}

	savefile	<- paste0('Rdata/Evenodd_', pv, '_', phen, '_', type, '_', seed, '.Rdata')
	sinkfile	<- paste0( 'Rout/Evenodd_', pv, '_', phen, '_', type, '_', seed, '.Rout')
	print( sinkfile )
	if( file.exists( savefile ) | file.exists( sinkfile ) ) next
	print( sinkfile )
	sink( sinkfile )

	withBMI	<- (! phen %in% c( 'BMI', 'T2D', 'disease_CARDIOVASCULAR','disease_HI_CHOL_SELF_REP'))
	binary	<- phen %in% binphens
	nPC			<- ifelse( type == 'int40', 40, 10 )

	into	<- 'Odd' 
	inte	<- 'Even'
	m			<- get_pheno_and_covar(phen, pv, internal_prs_40=( type %in% c( 'int40', 'bin10' ) ), seed=seed )

	pve_all	<- single_var_exp(m, withBMI=withBMI, nPC=nPC, phen, 'All'		, binary=binary )
	pve_1		<- single_var_exp(m, withBMI=withBMI, nPC=nPC, phen, into			, binary=binary	)
	pve_2		<- single_var_exp(m, withBMI=withBMI, nPC=nPC, phen, inte			, binary=binary	)

	res1	<- summary(regress(	m, withBMI=withBMI, nPC=nPC, outcome=into		, binary=FALSE ))$residuals
	res2	<- summary(regress(	m, withBMI=withBMI, nPC=nPC, outcome=inte		, binary=FALSE ))$residuals
	vars_cor <- cor(res1,res2)
	
	s			<- summary(regress(	m, withBMI=withBMI, nPC=nPC, outcome='Pheno', binary=binary, vars=c('Odd','Even'), interaction='Odd * Even'))
	beta_inter <- s$coefficients['Odd:Even',1]
	pval_inter <- s$coefficients['Odd:Even',4]

	dat	<- c(pve_all, pve_1, pve_2, vars_cor, beta_inter, -log10(pval_inter))
	print(dat)
	save( dat, file=savefile )

	sink()
}
