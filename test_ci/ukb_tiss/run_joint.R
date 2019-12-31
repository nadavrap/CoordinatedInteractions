rm( list=ls() )
load( '../ukb_eo/Rdata/bestps.Rdata' )
bestps['ext','disease_ASTHMA_DIAGNOSED']	<- '0.0001'
source('../functions.R')
source('../ukb_tissXtiss/setup.R')

for( type in c( 'ext', 'bin10', 'int40' ) )
	for( seed in seeds )
		for( phen in sample(phens) )
{
	if( (type == 'ext'	& !phen %in% extphens) |
			(type == 'bin10'& !phen %in% binphens) |
			(type == 'int40'&  phen %in% binphens) ) next

	savefile	<- paste0('Rdata/Tissue_', bestps[type,phen], '_', phen, '_', type, '_', seed, '_joint.Rdata')
	sinkfile	<- paste0( 'Rout/Tissue_', bestps[type,phen], '_', phen, '_', type, '_', seed, '_joint.Rout' )
	if( file.exists( savefile ) | file.exists( sinkfile ) ) next
	sink( sinkfile )
	withBMI	<- (! phen %in% c( 'BMI', 'T2D', 'disease_CARDIOVASCULAR','disease_HI_CHOL_SELF_REP'))
	nPC			<- ifelse( type == 'int40', 40, 10 )
	binary	<- phen %in% binphens

try({
	m	<- get_pheno_covar_alltiss(phen, TISSUES, bestps[type,phen], internal_prs_40=( type %in% c( 'int40', 'bin10' ) ), seed=seed )

	into <- paste0( TISSUES, '_Odd'  )
	inte <- paste0( TISSUES, '_Even' )

	s1		<- summary(regress(		m, withBMI=withBMI, nPC=nPC, outcome='Pheno'	, binary=binary,
		vars=c( inte,	'Odd', 'Even','Odd * Even'), interaction=paste('Odd :' ,inte)	))

	s2		<- summary(regress(		m, withBMI=withBMI, nPC=nPC, outcome='Pheno'	, binary=binary,
		vars=c( into,	'Odd', 'Even','Odd * Even'), interaction=paste('Even :',into)	))

	betas1	<- s1$coefficients[paste('Odd' , inte, sep = ':'),1]
	betas2	<- s2$coefficients[paste('Even', into, sep = ':'),1]

	ps1			<- s1$coefficients[paste('Odd' , inte, sep = ':'),4]
	ps2			<- s2$coefficients[paste('Even', into, sep = ':'),4]
	
	dat	<- cbind( betas1, betas2, ps1, ps2 )
	print( dim( dat ) )
	colnames( dat )	<- c( 'b1', 'b2', 'p1', 'p2' )
	rownames( dat )	<- TISSUES
	print( dat )

	save( dat, file=savefile )
})
sink()
}
