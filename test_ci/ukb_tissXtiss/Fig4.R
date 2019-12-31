rm( list=ls() )
load( '../ukb_eo/Rdata/bestps.Rdata' )
load( '../ukb_tiss/Rdata/table_output.Rdata' )
source('setup.R')

output
pchs	<- c( 16, 16, 17 )
names( pchs )	<- types

repeats	<- c( 26, 27, 37  )
output	<- output[-repeats,]


#########################
pdf( '~/figs/AM/Fig4.pdf', width=17, height=4.5 )
layout( matrix( 1:5, 1, 5 ), widths=c(0.8,10,10,10,5.55555) )
par( mar=c(0,0,0,0) )
plot.new()

par( mar=c(5.2,3,1.2,1) )
phen_vec	<- pv_vec	<- NULL
for( type in c(  'int40', 'bin10', 'ext' ) ){

	pmat		<- matrix( NA, 2*S, 0 )
	phenvec	<- tissvec	<- NULL
	for( phen in sample(phens) ){
		typelab	<- ifelse( type == 'ext', 'extern', 'intern' )

		sub			<- which( output[,1] == nicephens[phen] & output[,3] == typelab )
		if( length(sub) <= 1 ) next

		if( (type == 'ext'	& !phen %in% extphens) |
				(type == 'bin10'& !phen %in% binphens) |
				(type == 'int40'&  phen %in% binphens) ) next
		pv_vec		<- c( pv_vec	, output[sub[1],4] )
		phen_vec	<- c( phen_vec, output[sub[1],1] )

		tisses<- sort( output[sub,2] )
		nt		<- length( tisses )

		tmp	<- array( NA,
			dim=c(					2*S																		, ntiss		, ntiss		),
			dimnames=list(	c(paste0(seeds,'1'),paste0(seeds,'2')), TISSUES	, TISSUES )
		)
		for( seed in as.character(seeds) ){
			savefile	<- paste0('Rdata/Tissue_evenodd_internalPRS_2019_08_21_', bestps[type,phen], '_', phen, '_', type, '_', seed, '.Rdata')
			if( ! file.exists( savefile ) ) next
			load( savefile )
			for( tiss1 in sample(TISSUES) )
				for( tiss2 in sample(TISSUES) )
					tmp[paste0(seed,'1'),tiss1,tiss2]	<- as.numeric( (dat[[tiss1]])[ which( (dat[[tiss1]])[,1] == tiss2 ),3] )
			rm( dat )
		}
		
		for( seed in as.character(seeds) ){
			mat	<- t(tmp[paste0(seed,'1'),,])
			mat[upper.tri(mat,diag=T)]<- NA
			tmp[paste0(seed,'2'),,]		<- mat

			mat	<- tmp[paste0(seed,'1'),,]
			mat[upper.tri(mat,diag=T)]<- NA
			tmp[paste0(seed,'1'),,]		<- mat
		}
		for( tiss in TISSUES )
		stopifnot( mean( is.na( tmp[paste0(seed,'1'),tiss,] ) ) == mean( is.na( tmp[paste0(seed,'2'),tiss,] ) ))

		for( i in 1:nt )
			for( j in 1:nt )
		{
			if( i == j ) next
			tiss1	<- TISSUES[ which( nicetiss == tisses[i] ) ]
			tiss2	<- TISSUES[ which( nicetiss == tisses[j] ) ]
			if( all( is.na( tmp[,tiss1,tiss2] ) ) ) next

			pmat		<- cbind( pmat, tmp[,tiss1,tiss2] )
			phenvec	<- c( phenvec, phen )
			tissvec	<- c( tissvec, paste( nicetiss[tiss1], nicetiss[tiss2], sep=' : ' ) )
		}
		print( mean( is.na( pmat ) ) )
		pmat[ is.na( pmat ) ]	<- 0+1e-8*runif(sum(is.na(pmat)),0,1)
		dim( pmat )
	}

	nz		<- ncol( pmat )
	plot( range(-log10(1:nz/(1+nz))), range(c(0,6.5,pmat+.4),na.rm=T), type='n', xlab='', ylab='', axes=F )
	legend( 'topleft', leg=letters[which( type == types )], cex=1.4, bty='n' )
	axis(1)
	axis(2,at=0:3*2)
		mtext( side=1, line=4, expression( -log[10](p['Null']) )		, cex=1.1 )
	if( type == 'int40' ){
		mtext( side=2, line=4, expression( -log[10](p['TSxTS-CI']) ), cex=1.1 )
	}

	n_ext	<- sum( output == 'extern' )
	bonf	<- ifelse( type == 'ext', n_ext, nrow(output)-n_ext )
	cuts	<- -log10( .05/c( 2*length(seeds), 2*length(seeds) * bonf ) )
	for( k in 1:2 ){
		col	<- ifelse( k == 2, 'red', 'grey' )
		if( type == 'int40' )
		text( .2, cuts[k]+.2, lab=c('.05/(#Splits)','.05/(#Split*#TissPairs)')[k], col=col, cex=.7 )
		abline( h=cuts[k], col=col )
	}
	abline( a=0, b=1, col='grey', lwd=2, lty=3 )

	metap	<- apply( pmat, 2, max )
	ss		<- sort.list(metap, dec=T)

	for( i in 1:4 )
		text( -log10(i/(1+nz))		, metap[ss[i]]+1.3, col=cols[phenvec[ss[i]]], lab=tissvec[ss[i]], cex=.6, srt=90 )
	points( -log10(1:nz/(1+nz))	, metap[ss]				, col=cols[phenvec[ss]], pch=16 )

	#for( s.i in 1:nz )
	#	points( rep( -log10( s.i/(1+nz) ), 2*S ), pmat[,ss[s.i]], pch=16, cex=.4 , col=cols[phenvec[ss[s.i]]] ) #, col=cols[matphens[ss[1:nz]]]
}
par( mar=c(6.1,0.6,3.2,0) )
plot.new()

names(cols)		<- nicephens
phen_vec
legend( 'center'		, bty='n', cex=1.2, fill=cols[phen_vec], leg=as.expression(sapply( 1:length(phen_vec), function(i) bquote( paste( .(phen_vec[i]), ',  ', tau == .(pv_vec[i]) ) ) )) )
#legend( 'bottom', bty='n', cex=1.2, fill=cols[phen_vec], leg=as.expression(sapply( 1:length(phen_vec), function(i) bquote( paste( .(phen_vec[i]), ',  ', tau == .(pv_vec[i]) ) ) )) )

dev.off()
