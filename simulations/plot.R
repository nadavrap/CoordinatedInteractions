rm( list=ls() )
load('Rdata/setup.Rdata')

compilefile	<- 'Rdata/compiled.Rdata'
if( file.exists( compilefile ) ){
	load( compilefile )
} else {
	allout	<- array( NA, dim=c( nmode, 18, 1e4 ) )
	for( mode in 1:nmode )
		for( iter in 1:1e3 )
	tryCatch({
		if( iter == 1 ) print( mode )
		savefile	<- paste0( 'Rdata/', mode, '_', iter, '.Rdata' )
		if( !file.exists(savefile) ) next
		load( savefile )
		allout[mode,,iter]	<- c(t(out))
		rm( out )
	}, error=function(e) print(savefile))
	save( allout, file=compilefile )
}
round( apply( is.na(allout), 1:2, mean ), 2 )

######## plot vs N, h2, C
plotfxn	<- function( xs, ys, xlab, let='' ){
	plot( range(xs), 0:1, type='n', ylab='E-O Test Positive Rate', xlab=xlab, cex.lab=2 )
	abline( h=.05, col='grey', lwd=6, lty=2 )
	for( k in 1:3 ){
		lines( xs, rowMeans( ys[,6*k,] < .05, na.rm=T ), col=k, lwd=1.6 )
		points(xs, rowMeans( ys[,6*k,] < .05, na.rm=T ), col=k, pch=16, cex=1.8 )
	}
	legend( 'topleft', leg=let, cex=2.4, adj=c(2.1,1.1), bty='n' )
}

pdf( 'output/ci_supp_sims.pdf', width=13, height=4 )
layout( matrix(1:4,1,4), widths=c(6,6,6,2.9) )
par( mar=c(5,5,1,1) )
plotfxn( Ns		[Nseq]	, allout[Nseq,,]	, let='a', xlab='N' )
plotfxn( h2s	[h2seq]	, allout[h2seq,,]	, let='b', xlab=expression(h^2) )
plotfxn( gams [gamseq], allout[gamseq,,], let='c', xlab=expression(gamma) )

par( mar=c(1,1,1,1) )
plot.new()
legend( 'center', bty='n', leg=c( 'None', 'Uncoordinated', 'Coordinated' ), title='Epistasis Model', cex=1.7, fill=1:3 )
dev.off()


###### making tables
sim.table <- function(res, pop = FALSE, assort = FALSE){
	chi2s		<- apply(res[c(1,3,5),], 1, function(x) mean(x^2   ,na.rm=T))
	chi2sds	<- apply(res[c(1,3,5),], 1, function(x) sd  (x^2   ,na.rm=T))
	power		<- apply(res[c(2,4,6),], 1, function(x) mean(x<0.05,na.rm=T))
	# format means and sds
 	chi2s <- paste0(
 		formatC(chi2s  , digits = 1, format='f'), ' (', 
 		formatC(chi2sds, digits = 1, format='f'), ')'
 	)
 	power	<- format( power, nsmall=2, digits=1 )
	c( rbind( chi2s, power ) )
}

out			<- matrix( NA, 3*4+1, 2*3 )
out_ukb	<- 
out_trip<- matrix( '--', 3, 2*3 )
colnames(out)			<-
colnames(out_ukb)	<-
colnames(out_trip)<- c( 'prs', 'prs.pow', 'eocor', 'eocor.pow', 'eotest', 'eotest.pow' )

for( j in 1:3 ){ ### cycles over Null, Iso, CI
	for( k in 1:4 ) ### cycles over baseline, PopStruct, PopStruct+PC, assort
		out[(j-1)*4+k,]	<- sim.table(allout[c(1,3,17,4)[k],1:6+6*(j-1),])
	out_ukb [j,]			<- sim.table(allout[2             ,1:6+6*(j-1),])
	out_trip[j,]			<- sim.table(allout[5             ,1:6+6*(j-1),])
}
	out			[13,]			<- sim.table(allout[18            ,1:6+6*(j-1),])
assort		<- rep( c(0,0,0,1), 3 )
popstruct	<- rep( c(0,1,1,0), 3 )
pc_adjust	<- rep( c(0,0,1,0), 3 )
models		<- rep( c('null','iso','ci'), each=4 )
out				<- rbind( 
	cbind( models, popstruct, pc_adjust, assort, out[1:12,] ),
	c( 'ci-oracle', 0, 0, 0, out[13,] )
)

write.table(out			, file="output/sim_table_1.txt"       , quote=FALSE, row.names=FALSE, sep='\t')
write.table(out_trip, file="output/sim_table_trips.txt"   , quote=FALSE, row.names=FALSE, sep='\t')
write.table(out_ukb , file="output/sim_table_ukbbMAFs.txt", quote=FALSE, row.names=FALSE, sep='\t')
