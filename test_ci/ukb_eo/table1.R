rm( list=ls() )
source('setup.R')
load( '../ukb_eo/Rdata/bestps.Rdata' )

output	<- matrix( NA, 0, 7 )
colnames(output)	<- c( 'Phenotype', 'PRS Type', 'PRS_p_thresh', 'Significance Level', 'Min-p', 'Gamma', 'PRS-PVE' )

for( type in types ){
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

	########################
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

	gammas	<- apply( tmp[,,'Beta'	], 1, function(x) mean(as.numeric(x),na.rm=T) )
	#gammafrac<- apply( tmp[,,'Beta'	], 1, function(x) sum(as.numeric(x)<0,na.rm=T) )
	pves		<- apply( tmp[,,'%VE_global'	], 1, function(x) mean(as.numeric(x),na.rm=T) )

	pvals   <- tmp[,,'mLog10Pval']
	metap		<- apply( pvals, 1, function(x) max(as.numeric(x),na.rm=T) )

	cuts		<- -log10( .05/c( length(seeds), length(seeds)*P_bonf ) )
	print( cuts )
	print( sapply( metap, function(pv) paste0( rep( '*', ( pv > cuts[1] ) + ( pv > cuts[2] ) ), collapse='' ) ) )
	ncuts		<- apply( pvals, 1, function(x) sum(as.numeric(x) > cuts[1],na.rm=T) )
	ncuts1	<- apply( pvals, 1, function(x) sum(as.numeric(x) > cuts[2],na.rm=T) )

	output_loc	<- cbind(
		nicephens[phens],
		nicetypes[type],
		bestps[type,phens],
		sapply( metap, function(pv) paste0( rep( '*', ( pv > cuts[1] ) + ( pv > cuts[2] ) ), collapse='' ) ),
		format( 10^-metap, digit=2, scientific=T ),
		format( gammas   , digit=2, scientific=T ),
		format( pves     , digit=2, nsmall=2 )#,
		#ncuts,
		#ncuts1
	)
	output_loc	<- output_loc[ sort.list( output_loc[,7], dec=T ), ]

	output	<- rbind( output, output_loc )
}
write.table( output, file='~/figs/AM/table1.csv', row.names=FALSE, quote=FALSE, sep=',' )
save( output, file='Rdata/table_output.Rdata' )
