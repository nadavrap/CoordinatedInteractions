rm( list=ls() )
source('~/coordination_fresh/functions.R')
source('setup.R')

load( 'Rdata/compiled_output.Rdata' )# fdrs, emp_pvs, maxes, output,  


tab2	<- matrix( NA, 0, 8 ) 
#colnames(tab2)	<- c( 'Phenotype', 'PRS Type', 'Opt P-value Thresh', '%VE/AUC', 'Min-p', 'Empirical p', 'FDR', 'Signif Level' )
colnames(tab2)	<- c( 'Phenotype', 'PRS Type', 'Opt P-value Thresh', '%VE', 'Min-p', 'Empirical p', 'FDR', 'Signif Level' )
for( type in c( 'int', 'bin', 'ext' ) ){

	if( type == 'ext' ){
		phens	<- extphens
	} else if( type == 'bin' ){
		phens	<- binphens
	} else {
		phens	<- qphens
	}

	output_loc	<- NULL
	for( phen in phens ){

		asts	<- paste0( rep( '*', ( fdrs[type,phen,1] < .1 ) + ( fdrs[type,phen,1] < .01  ) ), collapse='' )

		pve		<- round( pves[type,phen,which.min( output[type,phen,,1] )] * 100, 2 )

		tau	<- pvs[ which.min( output[type,phen,,1] ) ]
		if( type == 'ext' & phen %in% extphens1 ) tau	<- '---' 

		output_loc	<- rbind( output_loc, c(
			nicephens[phen],
			nicetypes[type],
			tau,
			paste0( pve, '%' ),
			tidy( minps[type,phen,1] ),
			tidy( emp_pvs[type,phen] ),
			tidy( fdrs[type,phen,1] ),
			asts
		))

	}
	output_loc<- output_loc[ sort.list( output_loc[,1], dec=F ), ]
	tab2		<- rbind( tab2, output_loc ) 
	rm( output_loc )
}

write.table( tab2																					, file='~/coordination_fresh/figs/table_s3.csv'	, row.names=FALSE, quote=FALSE, sep=',' )
write.table( tab2[as.numeric(tab2[,7]) < .1,c(1,2,3,4,7)]	, file='~/coordination_fresh/figs/table2.csv'		, row.names=FALSE, quote=FALSE, sep=',' )
