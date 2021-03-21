rm( list=ls() )
source('~/coordination_fresh/functions.R')
source('setup.R')

eo_res					<- read.csv( '~/coordination_fresh/figs/table_s3.csv', head=T )
opt_pvs					<- as.character( eo_res[,3] ) 
names(opt_pvs)	<- paste0( as.character( eo_res[,1] ), '_', as.character( eo_res[,2] ) )

hitphens	<- as.character( eo_res[ eo_res[,'FDR'] < .1, 1 ] )
hitphens	<- phens[ nicephens[phens] %in% hitphens ]

opt	<- '_nondiag'

output	<- matrix( NA, 0, 8 )
colnames(output)	<- c( 'Phenotype', 'Tissue', 'PRS Type', 'TPRS %VE', 'p-value', 'Signif Level', 'p-va-exact', 'fdr' )

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
		for( tiss in tissues )
	{

		if( type == 'ext' & phen %in% extphens1 ){
			pv <- '1.0'
		} else {
			pv	<- opt_pvs[paste0( nicephens[phen], '_', nicetypes[type] )]
		}
		if( pv == 1 )	pv	<- '1.0'
		if( pv == 1e-4 )	pv	<- '0.0001'
		pv	<- as.character(pv) 

		load( paste0('Rdata/allpairs_', pv, '_', phen, '_', type, '_', tiss, '_', 0, opt, '.Rdata') ) 
		load( paste0('Rdata/hom_', pv, '_', phen, '_', type, '_', tiss, '.Rdata') ) 
		pve			<- round( fit_summ$coef['All_tiss',1]^2 * 100, 3 )

		output_loc	<- rbind( output_loc, c(
			nicephens[phen],
			nicetiss[tiss],
			nicetypes[type],
			paste0( pve, '%' ),
			tidy( pval_inter ),
			paste0( rep( '*', ( pval_inter < .05 ) + ( pval_inter < .05/ntiss  ) ), collapse='' ),
			pval_inter
		)) 
		rm( nicepv, fit_summ, pve, pval_inter )
	}

	fdrs	<- sapply( phens, function(phen){ 
		phensub	<- which( output_loc[,1] == nicephens[phen] )
		p.adjust(output_loc[phensub,7],'fdr')
	})
	fdrs	<- sapply( c( fdrs ), function(x) ifelse( x > .01, round( x, 2 ), format( x, digit=2, scientific=T ) ) )
	output_loc<- cbind( output_loc, fdrs )

	output_loc<- output_loc[ sort.list( output_loc[,1], dec=F ), ]
	output		<- rbind( output, output_loc ) 
	rm( output_loc )
}
write.table( output[																		, c(1:6,8)], file='~/coordination_fresh/figs/table_s4.csv', row.names=FALSE, quote=FALSE, sep=',' ) 
write.table( output[ as.numeric(output[,5]) < .05/ntiss	, 1:5			], file='~/coordination_fresh/figs/table3.csv'	, row.names=FALSE, quote=FALSE, sep=',' )
