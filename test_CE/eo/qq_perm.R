rm( list=ls() )
source('setup.R') 
load( 'Rdata/compiled_output.Rdata' )

xypoints	<- function(ys,ymax,...){ 
	ys	<- sort( -log10( ys ) )
	n		<- length( ys )
	xs	<- sort(-log10( 1:n/(1+n) ) )
	pch	<- rep( 16, length(ys) )
	if( any(ys>ymax) ){
		pch[ ys > ymax ]	<- 17
		ys [ ys > ymax ]	<- ymax
	} 
	points( xs, ys, pch=pch, ... ) 
} 

modes	<- c( 'tau_hom', 'FDR', 'minp', 'empirical' )
mains	<- c( 'Top Hom Tau', 'FDR-adjusted', 'Top Het Tau', 'Permutation-Adjusted' )
pdfs	<- c( 'SFigX', 'FigS3', 'SFigXX', 'FigS4' )

names(mains)	<- names(pdfs)	<- modes

for( mode in modes ){

	#pdf( paste0( '~/coordination_fresh/figs/qqplot_eo_', mode, '.pdf' ), width=12, height=6 )
	pdf( paste0( '~/coordination_fresh/figs/', pdfs[mode], '.pdf' ), width=12, height=6 )
	par( mfrow=c(1,2) )

	for( otype in c( 'UKBB', 'External' ) ){

		if( otype == 'UKBB' ){
			loctypes	<- c( 'bin', 'int' )
		} else {
			loctypes	<- c( 'ext' )
		}

		ymax	<- ifelse( mode=='empirical', 2.2, 7 )
		plot( c(0,ymax), c(0,ymax), type='n',xlab=expression( -log[10](p['Expected'])), ylab=expression( -log[10](p['CE']) ), main=paste0( mains[mode], '; ', otype, ' Summ Stats' ) )
		abline( a=0, b=1, col='grey' )

		if( mode == 'tau_hom' ){ 
			ys	<- t(sapply( phens, function(phen){
				type	<- ifelse( loctypes[1] == 'ext', 'ext', c( 'bin', 'int' )[ 1 + phen %in% qphens ] )
				if( is.na(hom_taus[type,phen])) return( rep(NA,nperm+1) ) 
				output[ type, phen, hom_taus[type,phen],]
			}))
			xypoints( ys[,-1] , ymax=7, col=2 ) 
			xypoints( ys[,1]  , ymax=7, col=1 ) 
			legend( cex=1.3, 'bottomright', fill=1:2, leg=c( 'Real Data', 'Permutated PRS' ), bty='n' )

		} else if( mode == 'FDR' ){
			xypoints( fdrs[loctypes,,-1] , ymax=7, col=2 ) 
			xypoints( fdrs[loctypes,,1]  , ymax=7, col=1 ) 
			legend( cex=1.3, 'bottomright', fill=1:2, leg=c( 'Real Data', 'Permutated PRS' ), bty='n' )

		} else if( mode == 'empirical' ){
			abline( h=-log10( 1/101 )	, col=4, lwd=1.5, lty=2 )
			abline( h=-log10( .05 )		, col=5, lwd=1.5, lty=2 )
			xypoints( emp_pvs[loctypes,]  , ymax=7, col=1 ) 
			legend( cex=1.3, 'bottomright', col=c(4,5), lty=2, lwd=2.5, leg=c( 'Max Possible (p=1/101)', 'p=.05' ), bty='n' )
			
		} else if( mode == 'minp' ){ 
			xypoints( sapply( 1:1e5, function(i) min(runif(npv)) ) , ymax=7, col=3 )
			xypoints( minps[loctypes,,-1] , ymax=7, col=2 ) 
			xypoints( minps[loctypes,,1]  , ymax=7, col=1 ) 
			legend( cex=1.3, 'bottomright', fill=1:3, leg=c( 'Real Data', 'Permutated PRS', 'Min(runif(nthresholds))' ), bty='n' )

		}

	}
	dev.off()
}
