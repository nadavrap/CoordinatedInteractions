rm( list=ls() )
source('setup.R')
load( 'Rdata/bestps.Rdata' )

pdf( '~/figs/AM/Fig2_gammas.pdf', width=10, height=10 )
layout( rbind( 1:2, 3:4, 5:6, c(7,6) ), width=c(10,3.2), height=c(10,10,10,1.3) )

for( type in types ){
	print( type )
	source('setup.R')
	if( type == 'ext' ){
		phens	<- extphens
		ylim	<- c(-.004,.0008)
	} else if( type == 'bin10' ){
		phens	<- binphens
		ylim	<- c(-2.0,1.0)
	} else {
		phens	<- qphens
		ylim	<- c(-1.1,.15)
	}
	P_bonf<- ifelse( type == 'ext', length(extphens), length(binphens)+length(qphens) )

	tmp	<- array( NA,
		dim=c( length(phens), length(seeds), 6 ),
		dimnames=list( phens, as.character(seeds), c('%VE_global', '%VE_1', '%VE_2', 'Cor', 'Beta', 'mLog10Pval' ) )
	)
	for( seed in seeds )
		for( phen in phens )
	{
		loadfile	<- paste0('Rdata/Evenodd_', bestps[type,phen], '_', phen, '_', type, '_', seed, '.Rdata')
		if( ! file.exists( loadfile ) ) next
		load( loadfile )
		tmp[phen,as.character(seed),]	<- dat
		rm( dat )
	}
	print( mean ( is.na( tmp ) ) )
	tmp[is.na(tmp)]	<- 0

	eocor   <- tmp[,,'Cor']       
	gamma   <- tmp[,,'Beta'] * 1e2
	pval    <- tmp[,,'mLog10Pval']

	#print( rowMeans( gamma > 0 ) )


	pvlist	<- sort.list( rowMeans( pval,na.rm=T ), dec=T )
	pvlist1	<- sort.list( rowSums( pval > -log10( .05/P_bonf/S ) ), dec=T )

	#ylim	<- range(gamma,na.rm=T)
	par( mar=c(2.8,5,.5,1) )
	plot( c(-.04,.13), ylim*1.05, type='n', xlab='', ylab=expression( gamma['coord'] %*% 100 ), cex.lab=1.5, cex.lab=1.5 )
	legend( 'topleft', leg=letters[which( type == types )], cex=1, bty='n' )
	abline( v=0, col='grey', lty=3, lwd=2 )
	abline( h=0, col='grey', lty=3, lwd=2 )

	upsub <- which( c(gamma) > ylim[2] )
	if( (uu <- length( upsub )) > 0 )
		points( eocor[upsub], 1.08*rep( ylim[2], uu ), col=cols[rep( phens, S )[upsub]], pch=17, cex=0.8 ) 

	lsub  <- which( c(gamma) < ylim[1] ) 
	if( (ll <- length( lsub )) > 0 ) 
		points( eocor[lsub] , 1.08*rep( ylim[1], ll ), col=cols[rep( phens, S )[lsub]]	, pch=17, cex=0.8 ) 

	for( phen in phens[pvlist] )
		points( eocor	[phen,], gamma[phen,], col=cols[phen], pch=16, cex=0.5 )

	for( phen in phens[pvlist1] )
		points( eocor	[phen,], gamma[phen,], col=cols[phen],
			pch=15 +		 ( pval[phen,] < -log10( .05/P_bonf/S ) ),
			cex=1.1 - .7*( pval[phen,] < -log10( .05/P_bonf/S ) ) )

	for( phen in phens[pvlist] )
		points( mean(	eocor[phen,]), mean(gamma[phen,]), col=cols[phen], pch=16, cex=2.3 )

	par( mar=c(.4,0,.2,0) )
	plot.new()
	topphens	<- phens[ pvlist ]
	legend( 'top'  , bty='n', cex=0.95, leg=as.expression(sapply( topphens, function(phen) bquote( paste( .(nicephens[phen]), ',  ', tau == .(bestps[type,phen]) ) ) )), fill=cols[topphens]	 )

}
legend( 'bottom'	, bty='n', cex=1.1, leg=c( 'Each chromosome splits','Bonferroni significant','Average over splits','Truncated to plot'), pt.cex=c(1.0,1.2,2.0,1.0), pch=c(16,15,16,17), y.int=1.6 )

plot.new()
mtext( side=1, line=-1.5, 'Even/Odd PRS Partial Correlation' )

dev.off()
