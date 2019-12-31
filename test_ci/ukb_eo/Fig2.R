rm( list=ls() )
source('setup.R')
for( pvaltype in c( '', '_stepdown', '_allpts' ) ){

load( 'Rdata/bestps.Rdata' )
if( pvaltype == '_stepdown' )
	bestps	<- apply( bestps, 1:2, function( pv ){
		if( is.na( pv ) ) return( NA )
		pvs[ which( pvs == pv )-1 ]
	})

pdf( paste0( '~/figs/AM/Fig2', pvaltype, '.pdf' ), width=10, height=6 )
layout( rbind( 1:2, 3:4, 5:6, c(7,6) ), width=c(10,2.5), height=c(10,10,10,1.8) )

for( type in types ){
	print( type )
	source('setup.R')
	if( type == 'ext' ){
		phens	<- extphens
	} else if( type == 'bin10' ){
		phens	<- binphens
	} else {
		phens	<- qphens
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
	pval    <- tmp[,,'mLog10Pval']

	print( c( type,pvaltype ) )
	xx	<- apply( pval, 1, max )
	xx	<- xx[ xx > -log10( .5/S ) ]
	print( round(xx,1) )

	pvlist	<- sort.list( rowMeans( pval,na.rm=T ), dec=T )
	pvlist1	<- sort.list( rowSums( pval > -log10( .05/S ) ), dec=T )

	print( range( eocor ) )
	print( range( pval ) )

	par( mar=c(2.8,5,.5,1) )
	if( pvaltype == '_allpts' ){
		xlim	<- c(-.04,.125)
	} else {
		xlim	<- c(-.02,.125)
	}
	plot( xlim, range(c(pval,0,8)), type='n', xlab='', ylab=expression(-log[10](p)), cex.lab=1.3, cex.axis=.9 )



	legend( 'topleft', leg=letters[which( type == types )], cex=1.4, bty='n' )
	abline( v=0, col='grey', lty=3, lwd=2 )
	cuts	<- -log10(c(.05/S,.05/P_bonf/S))
	for( k in 1:2 ){
		col	<- ifelse( k == 1, 'red', 'grey' )
		text( ifelse( pvaltype == '_stepdown', .11, .08 ), cuts[k]+.5, lab=c('.05/(#Splits)','.05/(#Splits*#Phen)')[k], col=col, cex=.85 )
		abline( h=cuts[k], col=col )
	}

	if( pvaltype != '_allpts' ){
		for( phen in phens[pvlist] ){
			i	<- which.max( pval[phen,] )
			points( eocor	[phen,i], pval[phen,i], col=cols[phen], pch=16, cex=1.6 )
		}
	} else {
		for( phen in phens[pvlist] )
			points( eocor	[phen,], pval		[phen,], col=cols[phen], pch=16, cex=0.6 )
		for( phen in phens[pvlist1] )
			points( eocor	[phen,], pval		[phen,], col=cols[phen],
				pch=15 +		 ( pval[phen,] < -log10( .05/P_bonf/S ) ),
				cex=1.1 - .7*( pval[phen,] < -log10( .05/P_bonf/S ) ) )
		for( phen in phens[pvlist] )
			points( mean(	eocor[phen,]), mean(pval[phen,]), col=cols[phen], pch=16, cex=2.3 )
	}

	par( mar=c(.4,0,.2,0) )
	plot.new()
	topphens	<- unique( c( phens[ pvlist1[ 1:min( length(pvlist1), 14 ) ] ], phens[ pvlist[ 1:min( length(pvlist), 14 ) ] ] ) )[1:min( length(pvlist), 14 )]
	legend( 'top'  , bty='n', cex=0.9, leg=as.expression(sapply( topphens, function(phen) bquote( paste( .(nicephens[phen]), ',  ', tau == .(bestps[type,phen]) ) ) )), col=cols[topphens], pch=16, pt.cex=1.5 )

}
if( pvaltype == '_allpts' )
legend( 'bottom'	, bty='n', cex=1.0	, leg=c( 'Each chromosome splits','Bonferroni significant','Average over splits'), pt.cex=c(1.0,1.2,2.0), pch=c(16,15,16), y.int=1.6 )

plot.new()
mtext( side=1, line=-1.5, 'Even/Odd PRS Partial Correlation', cex=.9 )

dev.off()
}
