rm( list=ls() )

savefile	<- 'Rdata/allX.Rdata'
if( ! file.exists( savefile ) ){

	load( 'Rdata/bestps.Rdata' )
	source('../functions.R')
	source('setup.R')

	phen	<- 'Height'
	type	<- 'int40'

	pv		<- bestps[type,phen]
	allX	<- NULL
	seeds	<- 1:100
	for( seed in seeds ){
		print( seed )
		m			<- get_pheno_and_covar(phen, pv, internal_prs_40=( type %in% c( 'int40', 'bin10' ) ), seed=seed )
		allX	<- cbind( allX, scale(m[,'Odd']*m[,'Even']) )
		rm( m )
		print( cov(allX) )
	}
	save( allX, file=savefile )
}
load( savefile )

quantile( cor( allX )[ lower.tri( diag( 100 ), diag=F ) ] )
