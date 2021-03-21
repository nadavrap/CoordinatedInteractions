rm( list=ls() )
source('~/coordination_fresh/functions.R')
source('setup.R') 
load( 'Rdata/compiled_output.Rdata' )

corr_file	<- 'Rdata/corrs.Rdata'
if(file.exists( corr_file )){
	load( corr_file )
} else { 
	all_cors	<- array( NA,
		dim=c(					3			, P			, 22, 22 ),
		dimnames=list(	types	, phens	, 1:22, 1:22  )
	) 
	for( type in c( 'ext', 'int', 'bin' ) )
		for( phen in phens )
	try({ 
		if( (type == 'ext' & !phen %in% extphens) |
				(type == 'bin' & !phen %in% binphens) |
				(type == 'int' &  phen %in% binphens) ) next
		m		<- get_pheno_and_covar(phen, hom_taus[type,phen], internal=( type %in% c( 'int', 'bin' ) ), seed=0, allchrs=TRUE ) 
		all_cors[type,phen,,]	<- cor( m[,paste0('X',1:22)] )
		rm( m ) 
	}) 
	save( all_cors, file=corr_file ) 
} 

allb_file	<- 'Rdata/allbs.Rdata'
if(file.exists( allb_file )){
	load( allb_file )
} else { 
	allbs	<- array( NA,
		dim=c(					3			, P			, 22, 22 ),
		dimnames=list(	types	, phens	, 1:22, 1:22  )
	) 
	for( type in c( 'ext', 'int', 'bin' ) )
		for( phen in phens )
	{ 
		if( (type == 'ext' & !phen %in% extphens) |
				(type == 'bin' & !phen %in% binphens) |
				(type == 'int' &  phen %in% binphens) ) next 
		load( paste0('Rdata/allpairs_', hom_taus[type,phen], '_', phen, '_', type, '_', 0, '.Rdata') )
		allbs[type,phen,,]	<- sapply( 1:22, function(i) sapply( 1:22, function(j){
			x	<- NA
			try( x	<- fit_summ[paste0('X',i,':','X',j),1], silent=T )
			x
		})) 
	}
	save( allbs, file=allb_file ) 
}

allp_file	<- 'Rdata/allps.Rdata'
if(file.exists( allp_file )){
	load( allp_file )
} else { 
	allps	<- array( NA,
		dim=c(					3			, P			, 22, 22 ),
		dimnames=list(	types	, phens	, 1:22, 1:22  )
	) 
	for( type in c( 'ext', 'int', 'bin' ) )
		for( phen in phens )
	{ 
		if( (type == 'ext' & !phen %in% extphens) |
				(type == 'bin' & !phen %in% binphens) |
				(type == 'int' &  phen %in% binphens) ) next 
		load( paste0('Rdata/allpairs_', hom_taus[type,phen], '_', phen, '_', type, '_', 0, '.Rdata') )
		allps[type,phen,,]	<- sapply( 1:22, function(i) sapply( 1:22, function(j){
			x	<- NA
			try( x	<- fit_summ[paste0('X',i,':','X',j),4], silent=T )
			x
		})) 
	}
	save( allps, file=allp_file ) 
}

mains		<-  c( 'Internal, Quant Traits', 'Internal, Bin Traits', 'External' )
names(mains)	<- types

pdf( paste0( '~/coordination_fresh/figs/FigS2.pdf' ), width=10, height=10 )
par( mfcol=c(3,3) )
par( mar=c(4.3,4.3,1,1) )

for( type in types )
	for( jj in 1:3 )
{
	if( type == 'ext' ){
		phens	<- extphens
		ylim	<- c(-.004,.0008)
	} else if( type == 'bin' ){
		phens	<- binphens
		ylim	<- c(-2.0,1.0)
	} else {
		phens	<- qphens
		ylim	<- c(-1.1,.15)
	}
	P_bonf<- ifelse( type == 'ext', length(extphens), length(binphens)+length(qphens) ) 

	x	<- fdrs[type,phens,1]
	#y	<- -log10(sapply( phens, function(p){
	y	<- sapply( phens, function(p){
		thetas	<- c( all_cors[type,p,,] )
		if( jj == 1 ){
			gamps		<- c( allbs		[type,p,,] )
		} else if( jj == 2 ){
			gamps		<- c( allps		[type,p,,] )
		} else {
			gamps		<- -log10( c( allps		[type,p,,] ) )
		}
		summary( lm( gamps ~ thetas ) )$coef[2,4]
	})

	ylab	<- c(
		'beta_{ij} ~ theta_{ij}',
		'pval_{ij} ~ theta_{ij}',
		'-log10(pval_{ij}) ~ theta_{ij}'
		#'Regression p for:\nEO pval_{ij} ~ theta_{ij}',
		#'Regression p for:\nEO pval_{ij} ~ theta_{ij}',
		#'Regression p for:\nEO -log10(pval_{ij}) ~  theta_{ij}'
	)[jj]

	plot(x,y,ylim=0:1, type='n',xlab='FDR for EO Test of Coordination', ylab=ylab, main='' )

	abline( h=0, col='grey', lwd=1.4 )
	abline( h=1, col='grey', lwd=1.4 )

	points(x,y,pch=16)
	legend( 'bottomright', bty='n', leg=mains[type], cex=1.2 )

}
dev.off()
