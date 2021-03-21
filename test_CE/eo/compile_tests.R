rm( list=ls() )
source('setup.R') 

homts	<- pves	<- array( NA,
	dim=c(					3			, P			, npv ),
	dimnames=list(	types	, phens	, as.character(pvs) )
) 
for( type in types )
	for( phen in phens )
		for( pv in pvs )
{ 
	if( (type == 'ext' & !phen %in% extphens) |
			(type == 'bin' & !phen %in% binphens) |
			(type == 'int' &  phen %in% binphens) ) next 

	if( type == 'ext' & phen %in% extphens1 & pv != '1.0' ) next

	load( paste0('Rdata/hom_', pv, '_', phen, '_', type, '.Rdata') )

	if( phen %in% binphens ){
		pves [type,phen,as.character(pv)]	<- auc
	} else {
		pves [type,phen,as.character(pv)]	<- fit_summ$coef['All',1]^2
		homts[type,phen,as.character(pv)]	<- fit_summ$coef['All',3]
	}

	rm( fit_summ )
}
hom_taus	<- apply( homts, 1:2, function(x) ifelse( all(is.na(x)), NA, pvs[which.max(x)] ) ) 

output	<- array( NA,
	dim=c(					3			, P			, npv, 1+nperm ),
	dimnames=list(	types	, phens	, as.character(pvs), as.character( 0:100 )  )
) 
for( perm.i in 0:100 )
	for( type in types )
		for( pv in pvs )
			for( phen in phens )
{

	print( perm.i )
	if( (type == 'ext' & !phen %in% extphens) |
			(type == 'bin' & !phen %in% binphens) |
			(type == 'int' &  phen %in% binphens) ) next 
	if( type == 'ext' & phen %in% extphens1 & pv != '1.0' ) next

	load( paste0('Rdata/allpairs_', pv, '_', phen, '_', type, '_', perm.i, '.Rdata') )
	# confirmed this is due to non-convergence of glm in 2/2 evaluated cases: if( is.na( pval_inter ) ){ print( loadfile ); print( fit_summ ) }
	output[type,phen,as.character(pv),perm.i+1]	<- pval_inter
	#allbetas<- summary(lm1)$coef[subs,1]
	#allsds	<- summary(lm1)$coef[subs,2]
	#allps		<- summary(lm1)$coef[subs,4]
	rm( pval_inter )
}
#fdrs		<- apply( output, c(1,2,4), function(x) min(p.adjust(x,'fdr')) ) 
fdrs		<- apply( output, c(1,2,4), function(x) min(p.adjust(x,'fdr'),na.rm=T) ) 
#minps		<- apply( output, c(1,2,4), function(x) min(x) )
#emp_pvs	<- apply( minps	, 1:2			, function(x) mean( x[1] >= x ) ) 

apply( is.na(output), 1:2, sum )
minps		<- apply( output, c(1,2,4), function(x) ifelse( all( is.na(x) ), NA, min(x,na.rm=T) ) )
emp_pvs	<- apply( minps	, 1:2			, function(x) mean( x[1] >= x,na.rm=T ) )

save( homts, hom_taus, fdrs, emp_pvs, minps, output, pves, file='Rdata/compiled_output.Rdata' )
