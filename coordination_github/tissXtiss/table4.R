rm( list=ls() )
source('~/coordination_fresh/functions.R')
source('~/coordination_fresh/tiss/setup.R')

eo_res					<- read.csv( '~/coordination_fresh/figs/table_s3.csv', head=T )
opt_pvs					<- as.character( eo_res[,3] ) 
names(opt_pvs)	<- paste0( as.character( eo_res[,1] ), '_', as.character( eo_res[,2] ) )

hitphens	<- as.character( eo_res[ eo_res[,'FDR'] < .1, 1 ] )
hitphens	<- phens[ nicephens[phens] %in% hitphens ]

opt	<- '_nondiag' 
tiss_res	<- read.csv( '../figs/table_s4.csv', head=T )[,1:5]

output	<- matrix( NA, 0, 6 )
colnames(output)	<- c( 'Phenotype', 'Tissue 1', 'Tissue 2', 'PRS Type', 'p-value', 'p-va-exact' ) #, 'Signif Level'

for( type in c( 'int', 'bin', 'ext' ) ){

	if( type == 'ext' ){
		phens	<- extphens
	} else if( type == 'bin' ){
		phens	<- binphens
	} else {
		phens	<- qphens
	}
	phens	<- intersect( phens, hitphens )

	output_loc	<- NULL
	for( phen in phens )
		for( tiss1.i in 1:(ntiss-1) )
			for( tiss2.i in (tiss1.i+1):ntiss )
	try({

		tiss1	<- tissues[ tiss1.i ]
		tiss2	<- tissues[ tiss2.i ]

		#phenrows	<- intersect( which( tiss_res[,'Phenotype'] == nicephens[phen] ), which( tiss_res[,'PRS.Type'] == nicetypes[type] ) )
		#i	<- intersect( phenrows, which( tiss_res[,'Tissue'] == nicetiss[tiss1] ) )
		#j	<- intersect( phenrows, which( tiss_res[,'Tissue'] == nicetiss[tiss2] ) ) 
		#if( any( tiss_res[c(i,j),5] > .05/ntiss ) ) next

		rows1	<- intersect( which( tiss_res[,'Phenotype'] == nicephens[phen] ), which( tiss_res[,'Tissue'] == nicetiss[tiss1] ) )
		minp1	<- min( tiss_res[rows1,5] )

		rows2	<- intersect( which( tiss_res[,'Phenotype'] == nicephens[phen] ), which( tiss_res[,'Tissue'] == nicetiss[tiss2] ) )
		minp2	<- min( tiss_res[rows2,5] ) 
		if( any( c(minp1,minp2) > .05/ntiss ) ) next 

		if( type == 'ext' & phen %in% extphens1 ){
			pv <- '1.0'
		} else {
			pv	<- opt_pvs[paste0( nicephens[phen], '_', nicetypes[type] )]
		}
		if( pv == 1 )	pv	<- '1.0'
		if( pv == 1e-4 )	pv	<- '0.0001'
		pv	<- as.character(pv) 

		load( paste0('Rdata/allpairs_', pv, '_', phen, '_', type, '_', 	tiss1, '_', tiss2, '_', 0, opt, '.Rdata') ) 
		if( type == 'ext' & phen %in% extphens1 ) pv	<- '---' 

		output_loc	<- rbind( output_loc, c(
			nicephens[phen],
			nicetiss[tiss1],
			nicetiss[tiss2],
			nicetypes[type],
			tidy( pval_inter ),
			##paste0( rep( '*', ( pval_inter < .05 ) + ( pval_inter < .05/ntiss  ) ), collapse='' ),
			pval_inter
		)) 
		rm( pval_inter )
	})

	output_loc<- output_loc[ sort.list( output_loc[,1], dec=F ), ]
	output		<- rbind( output, output_loc ) 
	rm( output_loc )
}

write.table( output[,-6]	, file='~/coordination_fresh/figs/table_s5.csv'	, row.names=FALSE, quote=FALSE, sep=',' )

ints	<- which( output[,4] != 'Extern' )
subi	<- as.numeric(output[ints,6]) < .05/nrow(output[ints,])
sube	<- as.numeric(output[-ints,6]) < .05/nrow(output[-ints,])
write.table( output[c(subi,sube), -6 ]	, file='~/coordination_fresh/figs/table4.csv'		, row.names=FALSE, quote=FALSE, sep=',' )
