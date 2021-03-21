rm( list=ls() )
source('setup.R')

eo_res	<- read.csv( '~/coordination_fresh/figs/table_s3.csv', head=T )
eo_res	<- eo_res[ eo_res[,2] != 'Extern', ]
hitphens	<- as.character( eo_res[ eo_res[,'FDR'] < .1, 1 ] )
hitphens	<- phens[ nicephens[phens] %in% hitphens ]
extphens	<- intersect(extphens,hitphens)

npair	<- 22*21/2
X		<- matrix( 1:22, 22, 22 )
is	<- t(X)[lower.tri(X,diag=F)]
js	<- X   [lower.tri(X,diag=F)]
sub0<- paste0('X',is,':','X',js)

tmp	<- array( NA,
	dim=c( 2, length(extphens), length(pvs), npair ),
	dimnames=list( c( 'int', 'ext' ), extphens, as.character(pvs), sub0 )
) 
for( phen in extphens )
	for( pv in pvs )
		for( type in c( 'ext', 'int' ) )
{

	if( type == 'int' ){
		loadtype	<- ifelse( phen %in% binphens, 'bin', 'int' )
	} else {
		loadtype	<- 'ext'
	}

	loadfile	<- paste0('Rdata/allpairs_', pv, '_', phen, '_', loadtype, '_', 0, '.Rdata') 
	if( type == 'ext' & phen %in% extphens1 & ! file.exists( loadfile ) ) next
	load( loadfile )

	subs<- intersect( sub0, rownames(fit_summ) ) 
	tmp[type,phen,as.character(pv),subs]	<- fit_summ[subs,1]

	rm( fit_summ )
}

eo_res					<- read.csv( '~/coordination_fresh/figs/table_s3.csv', head=T )
opt_pvs					<- eo_res[,3]
names(opt_pvs)	<- paste0( eo_res[,1], '_', eo_res[,2] )

pdf( paste0( '~/coordination_fresh/figs/Fig2.pdf' ), width=5*3.5, height=2*3.5 )
par( mar=c(5,5,1.2,1.2) )
layout( matrix( 1:10, 2, 5, byrow=T ) )
for( phen in extphens[ sort.list(nicephens[extphens]) ] ){

	pv	<- opt_pvs[paste0( nicephens[phen], '_Extern' )]
	if( pv == 1 )	pv	<- '1.0'
	if( pv == 1e-4 )	pv	<- '0.0001'
	if( pv == '---')	pv	<- '1.0'
	pv	<- as.character(pv)
	y	<- tmp['ext',phen,pv,]

	pv	<- opt_pvs[paste0( nicephens[phen], '_', ifelse( phen %in% binphens, nicetypes[2], nicetypes[1] ) )]
	if( pv == 1 )	pv	<- '1.0'
	if( pv == 1e-4 )	pv	<- '0.0001'
	pv	<- as.character(pv)
	x	<- tmp['int',phen,pv,]

	lims	<- range( c(x,y), na.rm=T ) 
	#plot( lims, lims, type='n', main='', xlab=expression( Internal~gamma[ij] ) , ylab=expression( External~gamma[ij] ), cex.lab=1.9 )
	plot( lims, lims, type='n', main='', xlab='Internal CE Estimates', ylab='External CE Estimates', cex.lab=1.5 )

	abline( a=0, b=1, lty=1, col='grey' )
	abline( h=0     , lty=1, col='grey' )
	abline( v=0     , lty=1, col='grey' )

	points( x, y, pch=16, cex=.7 )

	myfit	<- lm( y~x-1 )
	pval	<- summary( myfit )$coef[1,4] 
	abline( myfit, col=2, lwd=2 ) 
	pval	<- ifelse( pval > .01, round( pval, 2 ), format( pval, digit=2, scientific=T ) )
	#legend( bty='n', 'topleft', leg=paste0( 'P=', pval ), col=2, lty=1, lwd=2, title=nicephens[phen], cex=1.3 ) 
	legend( bty='n', 'topleft', leg=paste0( 'P=', pval ), title=nicephens[phen], cex=1.5 ) 

}
dev.off()
