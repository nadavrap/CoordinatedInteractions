rm( list=ls() )
load( '../ukb_eo/Rdata/bestps.Rdata' )

output	<- matrix( NA, 0, 9 )
colnames(output)	<- c( 'Phenotype', 'Tissue', 'Int/Ext', 'PRS_p_thresh', 'Significance Level', 'Min-p', 'PRS-PVE', 'TissPRS-PVE', 'BGPRS-PVE' )

for( type in c( 'int40', 'bin10', 'ext' ) ){
	source('setup.R')
	if( type == 'ext' ){
		phens	<- extphens
	} else if( type == 'bin10' ){
		phens	<- binphens
	} else {
		phens	<- qphens
	}

	P_bonf<- length(qphens)+length(binphens)
	P			<- length(phens)
	matphens<- matrix( phens	, P, ntiss )
	mattiss	<- matrix( TISSUES, P, ntiss, byrow=T )

	########################
	tmp	<- array( NA,
		dim=c( length(phens), length(seeds), ntiss, 10 ),
		dimnames=list( phens, seeds, TISSUES, c('%VE_global', '%VE_tissue', '%VE_bg', '%VE_tissue_1', '%VE_tissue_2', 'Cor', 'BTeAo', 'mLog10PvalTeAo', 'BToAe', 'mLog10PvalToAe') )
	)
	for( phen in phens )
		for( seed in seeds )
	{
		savefile	<- paste0('Rdata/Tissue_', bestps[type,phen], '_', phen, '_', type, '_', seed, '.Rdata')
		if( ! file.exists( savefile ) ) next
		load( savefile )
		nameloc	<- dat[,1]
		for( tiss in TISSUES )
			tmp[phen,as.character(seed),tiss,]	<- as.numeric( dat[ which(nameloc==tiss),-1] )
		rm( dat )
	}

	pves		<- apply( tmp[,,,'%VE_global'	], c(1,3), function(x) mean(as.numeric(x),na.rm=T) )
	pves_bg	<- apply( tmp[,,,'%VE_bg'			], c(1,3), function(x) mean(as.numeric(x),na.rm=T) )
	pves_ts	<- apply( tmp[,,,'%VE_tissue'	], c(1,3), function(x) mean(as.numeric(x),na.rm=T) )
	pvals   <- c( as.numeric(tmp[,,,'mLog10PvalToAe']), as.numeric(tmp[,,,'mLog10PvalTeAo']))
	metap		<- apply( tmp[,,,c('mLog10PvalToAe','mLog10PvalTeAo')], c(1,3), function(x) max(as.numeric(x),na.rm=T) )

	cuts	<- -log10( .05/c( 2*length(seeds), 2*length(seeds)*ntiss, 2*length(seeds)*ntiss*P_bonf ) )
	sub		<- which( metap > cuts[1] )
	if( length( sub ) == 0 ) stop('shit')

	typelab	<- ifelse( type == 'ext', 'extern', 'intern' )
	print( typelab )
	for( i in sub )
		output	<- rbind( output, c(
			nicephens[matphens[i]],
			nicetiss [mattiss [i]],
			typelab,
			bestps[type,matphens[i]],
			paste0( rep( '*', sum( metap[i] > cuts ) ), collapse='' ),
			format( 10^-metap[i], digit=2, scientific=T ),
			format( pves     [i], digit=2, nsmall=1 ),
			format( pves_ts  [i], digit=2, nsmall=1 ),
			format( pves_bg  [i], digit=2, nsmall=1 )
		))

}
warnings()
output	<- output[ sort.list( output[,1] ), ]
write.table( output, file='~/figs/AM/table2.csv', row.names=FALSE, quote=FALSE, sep=',' )
save( output, file='Rdata/table_output.Rdata' )
