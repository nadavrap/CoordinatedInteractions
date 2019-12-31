rm( list=ls() )
load( '../ukb_eo/Rdata/bestps.Rdata' )
bestps['ext','disease_ASTHMA_DIAGNOSED']	<- '0.0001'
for( type in c( 'int40', 'ext' ) ){ #, 'bin10'
	source('setup.R')
	#seeds	<- c( 0, round( 7 + (2:19)^3/42 ), 192 )
	S			<- length(seeds)

	if( type == 'ext' ){
		xmax	<- 3
		phens	<- extphens
	} else if( type == 'bin10' ){
		xmax	<- 2.8
		phens	<- binphens
	} else {
		xmax	<- 3.5
		phens	<- qphens
	}
	nphen   <- length(phens)

	########################
	tmp	<- array( NA,
		dim=c(					nphen, 2*S																		, ntiss		, ntiss		, 2 ),
		dimnames=list(	phens, c(paste0(seeds,'1'),paste0(seeds,'2'))	, TISSUES	, TISSUES	, c('beta','pval') )
	)
	for( phen in phens )
		for( seed in as.character(seeds) )
	{
		savefile	<- paste0('Rdata/Tissue_evenodd_internalPRS_2019_08_21_', bestps[type,phen], '_', phen, '_', type, '_', seed, '.Rdata')
		if( ! file.exists( savefile ) ) next
		load( savefile )
		#for( tiss1 in TISSUES ){
		#	datloc	<- dat[[tiss1]]
		#	nameloc	<- datloc[,1]
		#	for( tiss2 in TISSUES )
		#		tmp[phen,paste0(seed,'1'),tiss1,tiss2,]	<- as.numeric( datloc[ which(nameloc==tiss2),-(1:2)] )
		#}
		for( tiss1 in TISSUES )
			for( tiss2 in TISSUES )
		tryCatch({
				tmp[phen,paste0(seed,'1'),tiss1,tiss2,]	<- as.numeric( (dat[[tiss1]])[ which( (dat[[tiss1]])[,1] == tiss2 ),2:3] )
			}, error=function(e){
		print( savefile )
			print( tiss1 )
			print( tiss2 )
			print( phen )
			print( as.numeric( (dat[[tiss1]])[ which( (dat[[tiss1]])[,1] == tiss2 ),2:3] ) )
			print( tmp[phen,paste0(seed,'1'),tiss1,tiss2,] )
			})

		for( k in 1:2 ){

			mat	<- t(tmp[phen,paste0(seed,'1'),,,k])
			mat[upper.tri(mat,diag=T)]		<- NA
			tmp[phen,paste0(seed,'2'),,,k]	<- mat


			mat	<- tmp[phen,paste0(seed,'1'),,,k]
			mat[upper.tri(mat,diag=T)]		<- NA
			tmp[phen,paste0(seed,'1'),,,k]	<- mat
		}
		stopifnot( mean( is.na( tmp[phen,paste0(seed,'1'),tiss1,,] ) ) == mean( is.na( tmp[phen,paste0(seed,'2'),tiss1,,] ) ))
		rm( dat )
	}
	print( mean ( is.na( tmp ) ) )
	print( 1-( ntiss*(ntiss-1)/2 ) / ntiss^2 )
	tmp[is.na(tmp)]	<- 0

	#########################
	matphens<- array( NA, dim=c(nphen,ntiss,ntiss), dimnames=list(phens,nicetiss,nicetiss))
	mattiss	<- array( NA, dim=c(nphen,ntiss,ntiss), dimnames=list(phens,nicetiss,nicetiss))
	for(p in phens)
		matphens[p,,]	<- p
	for(t1 in nicetiss)
		for(t2 in nicetiss)
			mattiss[,t1,t2]	<- paste( t1, t2, sep=':' )

	#########################
	pdf( paste0( '~/figs/AM/condition_TxT', type, '.pdf' ), width=12, height=6 )
	par( mar=c(4.1,4.1,.2,.2) )
	try({
		pvals   <- c(tmp[,,,,'pval'])
									
		plot( c(0,xmax), range(c(0,8,pvals+.4),na.rm=T), type='n', xlab=expression( -log[10](p['Null']) ), ylab=expression( -log[10](p['CI']) ), axes=F )
		axis(1)
		axis(2)

		cuts	<- -log10( .05/c( 2*length(seeds), 2*length(seeds)*( ntiss*(ntiss-1)/2 ), 2*length(seeds)*( ntiss*(ntiss-1)/2 )*nphen ) )
		for( k in 1:3 ){
			text( .3, cuts[k]+.1, lab=c('.05/(#Splits)','.05/(#Split*#TissPair)','.05/(#Splits*#TissPair*#Phen)')[k], col='grey', cex=.7 )
			abline( h=cuts[k], col='grey' )
		}
		abline( a=0, b=1, col='grey', lwd=2, lty=3 )

		metap	<- apply( tmp[,,,,'pval'], c(1,3,4), function(x) max(as.numeric(x),na.rm=T) )
		ss		<- sort.list(metap, dec=T)

		nz		<- (ntiss*(ntiss-1)*nphen/2)
	#	print( mean( metap == 0 ) )
	#	print( mean( metap[ss[  1:nz ]] == 0 ) )
	#	print( mean( metap[ss[-(1:nz)]] == 0 ) )

		nsel	<- 15
		for( i in 1:1 )
			text( -log10( i/(1+nz) ), metap[ss[i]]+0.3, col=cols[matphens[ss[i]]], lab=mattiss[ss[i]], cex=.5 )
		for( i in 2:nsel )
			text( -log10( i/(1+nz) ), metap[ss[i]]+1.1, col=cols[matphens[ss[i]]], lab=mattiss[ss[i]], cex=.5, srt=90 )

		points( -log10( 1:nz/(1+nz) ), metap[ss[1:nz]], col=cols[matphens[ss[1:nz]]], pch=16 )
		for( s.i in 1:nz )try({
			metap	<- apply( tmp[,,,,'pval'], c(1,3,4), function(x) as.numeric(x)[s.i] )
			points( -log10( 1:nz/(1+nz) ), metap[ss[1:nz]], col=cols[matphens[ss[1:nz]]], pch=16, cex=.4 )
		})
		#	metap	<- apply( tmp[,,,,'pval'], c(1,3), function(x) mean(as.numeric(x)) )
		#	points( -log10( 1:n/(1+n) ), metap[ss], col=cols[matphens[ss]], pch=17, cex=1.2 )


		phenloc	<- unique(matphens[ ss[1:nsel] ])
		legend( 'topleft', fill=cols[phenloc], leg=paste0( nicephens[phenloc], ', p=', bestps[type,phenloc] ), bty='n', cex=.8 )
	})
	dev.off()
}
