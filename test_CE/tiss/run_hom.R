rm( list=ls() )
source('~/coordination_fresh/functions.R')
source('setup.R') 

eo_res					<- read.csv( '~/coordination_fresh/figs/table_s3.csv', head=T )
opt_pvs					<- eo_res[,3]
names(opt_pvs)	<- paste0( eo_res[,1], '_', eo_res[,2] )

for( type in c( 'ext', 'int', 'bin' ) )
	for( phen in sample(phens) )
		for( tiss in tissues )
try({
	if( (type == 'ext' & !phen %in% extphens) |
			(type == 'bin' & !phen %in% binphens) |
			(type == 'int' &  phen %in% binphens) ) next

	pv	<- opt_pvs[paste0( nicephens[phen], '_', nicetypes[type] )]
	if( pv == 1    )	pv	<- '1.0'
	if( pv == 1e-4 )	pv	<- '0.0001'
	if( pv == '---')	pv	<- '1.0'
	pv	<- as.character(pv) 


	savefile	<- paste0('Rdata/hom_', pv, '_', phen, '_', type, '_', tiss, '.Rdata')
	sinkfile	<- paste0( 'Rout/hom_', pv, '_', phen, '_', type, '_', tiss, '.Rout')
	if( file.exists( savefile ) | file.exists( sinkfile ) ) next
	print( sinkfile )
	sink( sinkfile )

	withBMI	<- ( phen != 'BMI' )
	binary	<- phen %in% binphens

	m			<- get_pheno_covar_wtiss( phen, tiss, pv, internal=( type %in% c( 'int', 'bin' ) ), allchrs=TRUE, return_bg=F )

	unused1	<- which( grepl('genotype_measurement_batch_f22000_0_0', colnames(m)) )
	unused2	<- which( grepl('uk_biobank_assessment_centre_f54_0_0' , colnames(m)) )
	m	<- m[,-c(unused1,unused2)]
	gc() 

	if( binary )
		m[,'Pheno']	<-  as.numeric(m[,'Pheno'])-1 

	m[,-1]	<- scale( m[,-1] ) 

	stopifnot( all( colMeans(is.na(m)) %in% 0:1 ) )
	nacols	<- which( colMeans(is.na(m)) > 0 )
	if( length( nacols ) > 0 )
		m[,nacols]	<- 0
	stopifnot( all( colMeans(is.na(m)) == 0 ) )

  f0 <- paste('Pheno', "~", paste(c(
		'All',
		'All_tiss',
		'age_when_attended_assessment_centre_f21003_0_0',
		'sex_f31_0_0',
		ifelse(withBMI,'BMI',''),
		paste("genetic_principal_components_f22009_0_",1:nPC,sep="")),
	collapse="+")) 

	if( binary ){
		lm0	<- glm(as.formula(f0), data=m,family=binomial)
	} else {
		lm0	<- lm( as.formula(f0), data=m)
	} 
	fit_summ	<- summary(lm0)
	save(  fit_summ, file=savefile )
	rm( m, fit_summ, lm0 ); gc()

	sink()
})
