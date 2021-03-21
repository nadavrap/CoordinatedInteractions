rm( list=ls() )
source('setup.R') 

eo_res					<- read.csv( '~/coordination_fresh/figs/table_s3.csv', head=T )
opt_pvs					<- as.character( eo_res[,3] ) 
names(opt_pvs)	<- paste0( as.character( eo_res[,1] ), '_', as.character( eo_res[,2] ) ) 

hitphens	<- as.character( eo_res[ eo_res[,'FDR'] < .1, 1 ] )
hitphens	<- phens[ nicephens[phens] %in% hitphens ]
phens			<- hitphens
P					<- length( phens )

nperm	<- 5

for( opt in c( '_bg', '_nondiag' ) ){

output	<- array( NA,
	dim=c(				3			, P			, ntiss		, 1+nperm			),
	dimnames=list(types	, phens	, tissues , 1:(nperm+1) )
) 

for( perm.i in 0:nperm )
	for( type in c( 'int', 'bin', 'ext' ) )
		for( phen in sample(phens) )
			for( tiss in tissues )
try({

	if( (type == 'ext' & !phen %in% extphens) |
			(type == 'bin' & !phen %in% binphens) |
			(type == 'int' &  phen %in% binphens) ) next 

	if( type == 'ext' & phen %in% extphens1 ){
		pv <- '1.0'
	} else {
		pv	<- opt_pvs[paste0( nicephens[phen], '_', nicetypes[type] )]
	}
	if( pv == '---')	pv	<- '1.0'
	if( pv == 1 )	pv	<- '1.0'
	#if( pv == 1e-4 )	pv	<- '0.0001'
	pv	<- as.character(pv) 

		tryCatch({
	load( paste0('Rdata/allpairs_', pv, '_', phen, '_', type, '_', tiss, '_', perm.i, opt, '.Rdata') )
		},error=function(e) print( paste0('Rdata/allpairs_', pv, '_', phen, '_', type, '_', tiss, '_', perm.i, opt, '.Rdata') ) )
	output[type,phen,tiss,perm.i+1]	<- pval_inter
	rm( pval_inter )
})

xypoints	<- function(ys,ymax,...)try({
	ys	<- sort( -log10( ys ) )
	n		<- length( ys )
	xs	<- sort(-log10( 1:n/(1+n) ) )
	pch	<- rep( 16, length(ys) )
	if( any(ys>ymax) ){
		pch[ ys > ymax ]	<- 17
		ys [ ys > ymax ]	<- ymax
	} 
	points( xs, ys, pch=pch, ... ) 
})

pdf( paste0( '~/coordination_fresh/figs/FigS5', opt, '.pdf' ), width=12, height=6 )
par( mfrow=c(1,2) ) 
for( otype in c( 'UKBB', 'External' ) ){

	if( otype == 'UKBB' ){
		loctypes	<- c( 'bin', 'int' )
	} else {
		loctypes	<- c( 'ext' )
	}

	plot( c(0,7), c(0,7), type='n',xlab=expression( -log[10](p['Expected'])), ylab=expression( -log[10](p['CE']) ), main=paste0( otype, ' Summ Stats' ) )
	abline( a=0, b=1, col='grey' )
	xypoints( c(output[loctypes,phens,,-1])  , ymax=7, col=2 ) 
	ys	<-	output[loctypes,phens,,1]
	xypoints( c(ys)  , ymax=7, col=1 ) 

	legend( cex=1.3, 'bottomright', fill=1:2, leg=c( 'Real Data', 'Permutated PRS' ), bty='n' )

}
dev.off()

pdf( paste0( '~/coordination_fresh/figs/qqplot_bytiss', opt, '.pdf' ), width=12, height=6 )
par( mfrow=c(3,6) )
for( tiss in tissues ){ 
		plot( c(0,7), c(0,7), type='n',xlab=expression( -log[10](p['Expected'])), ylab=expression( -log[10](p['CE']) ), main=tiss )
		abline( a=0, b=1, col='grey' )
		xypoints( c(output[,,tiss,-1])  , ymax=7, col=2 ) 
		xypoints( c(output[,,tiss, 1])  , ymax=7, col=1 ) 
		#legend( cex=1.3, 'bottomright', fill=1:2, leg=c( 'Real Data', 'Permutated PRS' ), bty='n' ) 
}
dev.off()

pdf( paste0( '~/coordination_fresh/figs/qqplot_byphen', opt, '.pdf' ), width=12, height=6 )
par( mfrow=c(3,6) )
for( phen in phens ){ 
		plot( c(0,7), c(0,7), type='n',xlab=expression( -log[10](p['Expected'])), ylab=expression( -log[10](p['CE']) ), main=phen )
		abline( a=0, b=1, col='grey' )
		xypoints( c(output[,phen,,-1])  , ymax=7, col=2 ) 
		xypoints( c(output[,phen,, 1])  , ymax=7, col=1 ) 
		#legend( cex=1.3, 'bottomright', fill=1:2, leg=c( 'Real Data', 'Permutated PRS' ), bty='n' ) 
}
dev.off()

}
