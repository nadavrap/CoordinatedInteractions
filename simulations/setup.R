rm( list=ls() )

nmode		<- 18
Ns			<- rep( 1e4			, nmode )
h2s			<- rep( .5 			, nmode )
gams		<- rep( .2 			, nmode )
maftypes<- rep( 'unif'  , nmode )
simtypes<- rep( 'simple', nmode )
adj_pcs	<- rep( FALSE   , nmode )
trips		<- rep( FALSE   , nmode )
oracle	<- rep( FALSE   , nmode )

maftypes[2]	<- 'UKBB'
simtypes[3]	<- 'strat'
simtypes[4]	<- 'am'
trips   [5]	<- TRUE

Ns	[6:9]		<- c( 1e3, 2e3, 5e3, 2e4 )
h2s	[10:12]	<- c( .1, .3, .7 )
gams[13:16]	<- c( -.4, -.2, 0, .4 )


Nseq	<- c(6,7,8,1,9)
h2seq	<- c(10,11,1,12)
gamseq<- c(13,14,15,1,16)

simtypes[17]<- 'strat'
adj_pcs [17]<- TRUE

oracle	[18]<- TRUE

M	<- 1e3

save.image('Rdata/setup.Rdata')
