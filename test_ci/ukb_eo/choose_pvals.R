rm( list=ls() )
source('../functions.R')
source('setup.R')

bestps	<- array( NA, dim=c(3,P), dimnames=list( c( 'ext', 'int40', 'bin10' ), phens ) )
for( type in c( 'ext', 'int40', 'bin10' ) )
	for( phen in phens )
try({
	if( (type == 'ext'	& !phen %in% extphens) |
			(type == 'bin10'& !phen %in% binphens) |
			(type == 'int40'&  phen %in% binphens) ) next

	scores				<- rep( NA, npv )
	names(scores)	<- pvs
	for( pv in pvs ){
		load( paste0( 'Rdata/Evenodd_', pv, '_', phen, '_', type, '_', 0, '.Rdata') )
		scores[pv]	<- dat[1]
		rm( dat )
	}
	if( phen %in% binphens ){
	print( phen )
	print( scores )
	}
	if( any( is.na( scores ) ) ) stop()
	bestps[type,phen]	<- pvs[which.max(scores)]
},silent=T)

t(bestps)
bestps['ext','disease_ASTHMA_DIAGNOSED']	<- '0.0001'

save( bestps, file='Rdata/bestps.Rdata' )
