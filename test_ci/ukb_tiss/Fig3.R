rm( list=ls() )
source('setup.R')

nsels	<-  c( 26, 11, 20 )
names(nsels)	<- types

for( pvaltype in c( '', '_joint', '_nobg' ) )try({ #, '_stepdown'

load( '../ukb_eo/Rdata/bestps.Rdata' )
bestps['ext','disease_ASTHMA_DIAGNOSED']	<- '0.0001'

if( pvaltype == '_stepdown' )
	bestps	<- apply( bestps, 1:2, function( pv ){
		if( is.na( pv ) ) return( NA )
		pvs[ which( pvs == pv )-1 ]
	})

pdf( paste0( '~/figs/AM/Fig3', pvaltype, '.pdf' ), width=13, height=8 )
layout( rbind( 1:2, 3:4, 5:6, c(7,6) ), width=c(10,1.9), height=c(10,10,10,1.8) )

for( type in c( 'int40', 'bin10', 'ext' ) ){
	source('setup.R')
	if( type == 'ext' ){
		phens	<- extphens
	} else if( type == 'bin10' ){
		phens	<- binphens
	} else {
		phens	<- qphens
	}

	P_bonf<- ifelse( type == 'ext', length(extphens), length(binphens)+length(qphens) )
	P			<- length(phens)
	matphens<- matrix( phens	, P, ntiss )
	mattiss	<- matrix( TISSUES, P, ntiss, byrow=T )

	########################
	if( pvaltype == '_nobg' ){
		suff	<- '_nobg'
		tmp	<- array( NA,
			dim=c( length(phens), length(seeds), ntiss, 9 ),
			dimnames=list( phens, seeds, TISSUES, c('%VE_global', '%VE_tissue', '%VE_tissue_1', '%VE_tissue_2', 'Cor', 'BTeAo', 'mLog10PvalTeAo', 'BToAe', 'mLog10PvalToAe') )
		)
	} else if( pvaltype == '_joint' ){
		suff	<- '_joint'
		tmp	<- array( NA,
			dim=c( length(phens), length(seeds), ntiss, 4 ),
			dimnames=list( phens, seeds, TISSUES, c('BTeAo', 'BToAe', 'mLog10PvalTeAo', 'mLog10PvalToAe') )
		)
	} else {
		suff	<- ''
		tmp	<- array( NA,
			dim=c( length(phens), length(seeds), ntiss, 10 ),
			dimnames=list( phens, seeds, TISSUES, c('%VE_global', '%VE_tissue', '%VE_bg', '%VE_tissue_1', '%VE_tissue_2', 'Cor', 'BTeAo', 'mLog10PvalTeAo', 'BToAe', 'mLog10PvalToAe') )
		)
	}
	for( phen in phens )
		for( seed in seeds )
	{
		savefile	<- paste0('Rdata/Tissue_', bestps[type,phen], '_', phen, '_', type, '_', seed, suff, '.Rdata')
		if( ! file.exists( savefile ) ) next
		load( savefile )
		if( pvaltype=='_joint' )
			dat	<- cbind( rownames(dat), dat )
		stopifnot( all( dat[,1] == TISSUES ) )
		tmp[phen,as.character(seed),,]	<- as.numeric( dat[,-1] )
		rm( dat )
	}
	if( pvaltype=='_joint' )
		tmp[,,,3:4]	<- -log10( tmp[,,,3:4] )
	#print( apply ( is.na( tmp ), 1, mean ) )
	print( c( type, pvaltype ) )
	print( mean( is.na(tmp) ) )
	if( any( is.na( tmp ) ) )
	tmp[is.na(tmp)]	<- runif( sum(is.na(tmp)), 0, 1e-9 )

	########################
	par( mar=c(2.0,4.5,.5,1) )
	pvals   <- c( as.numeric(tmp[,,,'mLog10PvalToAe']), as.numeric(tmp[,,,'mLog10PvalTeAo']))
	metap	<- apply( tmp[,,,c('mLog10PvalToAe','mLog10PvalTeAo')], c(1,3), function(x) max(as.numeric(x),na.rm=T) )
	ss		<- sort.list(metap, dec=T)
	n			<- length( metap )
								
	plot( range(-log10( c(1,n)/(1+n) )) , range(c(0,7,pvals+.2),na.rm=T), type='n', xlab='', ylab=expression( -log[10](p['TS-CI']) ), axes=F )
	legend( 'topleft', leg=letters[which( type == types )], cex=1.4, bty='n' )
	axis(1)
	axis(2)

	cuts	<- -log10( .05/c( 2*length(seeds), 2*length(seeds)*ntiss, 2*length(seeds)*ntiss*P_bonf ) )
	for( k in 1:3 ){
		col	<- ifelse( k == 2, 'red', 'grey' )
		text( .3, cuts[k]+.2, lab=c('.05/(#Splits)','.05/(#Split*#Tiss)','.05/(#Splits*#Tiss*#Phen)')[k], col=col, cex=.7 )
		abline( h=cuts[k], col=col )
	}
	abline( a=0, b=1, col='grey', lwd=2, lty=3 )

	nsel	<- nsels[type]
	for( i in 1:4 )
		text( -log10( i/(1+n) ), metap[ss[i]]+.3, col=cols[matphens[ss[i]]], lab=nicetiss[mattiss[ss[i]]], cex=.6 )
	for( i in 5:nsel )
		text( -log10( i/(1+n) ), metap[ss[i]]+1.3, col=cols[matphens[ss[i]]], lab=nicetiss[mattiss[ss[i]]], cex=.6, srt=90 )

	sel <- ss[1:nsel]

	points( -log10( 1:n/(1+n) ), metap[ss], col=cols[matphens[ss]], pch=16 )
	#for( s.i in 1:length(ss) ){
	#	metap	<- apply( tmp[,,,c('mLog10PvalToAe','mLog10PvalTeAo')], c(1,3), function(x) as.numeric(x)[s.i] )
	#	points( -log10( 1:n/(1+n) ), metap[ss], col=cols[matphens[ss]], pch=16, cex=.4 )
	#}
	#	metap	<- apply( tmp[,,,c('mLog10PvalToAe','mLog10PvalTeAo')], c(1,3), function(x) median(as.numeric(x)) )
	#	points( -log10( 1:n/(1+n) ), metap[ss], col=cols[matphens[ss]], pch=17, cex=1.2 )

	topphens	<- unique(matphens[sel])

	par( mar=c(.4,0,.2,0) )
	plot.new()
	#legend( 'top', fill=cols[phenloc], leg=paste0( nicephens[phenloc], ', p=', bestps[type,phenloc] ), bty='n', cex=.8 )
	legend( 'top'  , bty='n', cex=0.95, leg=as.expression(sapply( topphens, function(phen) bquote( paste( .(nicephens[phen]), ',  ', tau == .(bestps[type,phen]) ) ) )), fill=cols[topphens]	 )
}
#legend( 'bottom'	, bty='n', cex=1.0	, leg=c( 'Each chromosome splits','Bonferroni significant','Average over splits'), pt.cex=c(1.0,1.2,2.0), pch=c(16,15,16), y.int=1.6 )

plot.new()
mtext( side=1, line=-1.5, expression( -log[10](p['Null']) ), cex=.8 )

dev.off()
})
