rm( list=ls() )
source('~/coordination_fresh/functions.R')
source('~/coordination_fresh/tiss/setup.R')

eo_res					<- read.csv( '../figs/table_s3.csv', head=T )
opt_pvs					<- eo_res[,3]
names(opt_pvs)	<- paste0( eo_res[,1], '_', eo_res[,2] ) 

hitphens	<- as.character( eo_res[ eo_res[,'FDR'] < .1, 1 ] )
hitphens	<- phens[ nicephens[phens] %in% hitphens ]

opt	<- '_nondiag' 
tiss_res	<- read.csv( '../figs/table_s4.csv', head=T )[,1:5]

for( perm.i in 0:0 )
	for( type in c( 'int', 'bin', 'ext' ) )
		for( phen in sample(hitphens) )
			for( tiss1.i in 1:(ntiss-1) )
				for( tiss2.i in (tiss1.i+1):ntiss )
try({

	if( (type == 'ext' & !phen %in% extphens) |
			(type == 'bin' & !phen %in% binphens) |
			(type == 'int' &  phen %in% binphens) ) next 

	tiss1	<- tissues[ tiss1.i ]
	tiss2	<- tissues[ tiss2.i ]

	rows1	<- intersect( which( tiss_res[,'Phenotype'] == nicephens[phen] ), which( tiss_res[,'Tissue'] == nicetiss[tiss1] ) )
	minp1	<- min( tiss_res[rows1,5] )

	rows2	<- intersect( which( tiss_res[,'Phenotype'] == nicephens[phen] ), which( tiss_res[,'Tissue'] == nicetiss[tiss2] ) )
	minp2	<- min( tiss_res[rows2,5] ) 
	if( any( c(minp1,minp2) > .05/ntiss ) ) next

	pv			<- opt_pvs[		paste0( nicephens[phen], '_', nicetypes[type] )] 
	if( pv == 1    )	pv	<- '1.0'
	if( pv == 1e-4 )	pv	<- '0.0001'
	if( pv == '---')	pv	<- '1.0'
	pv	<- as.character(pv) 

	savefile	<- paste0('Rdata/allpairs_', pv, '_', phen, '_', type, '_', tiss1, '_', tiss2, '_', perm.i, opt, '.Rdata')
	sinkfile	<- paste0( 'Rout/allpairs_', pv, '_', phen, '_', type, '_', tiss1, '_', tiss2, '_', perm.i, opt, '.Rout')
	if( file.exists( savefile ) | file.exists( sinkfile ) ) next
	print( sinkfile )
	sink( sinkfile ) 
try({

	withBMI	<- ( phen != 'BMI' )
	binary	<- phen %in% binphens

		#get_pheno_covar_wtiss( phen, tiss1, pv, internal=( type %in% c( 'int', 'bin' ) ), allchrs=TRUE, return_ID=TRUE, return_bg=( opt %in% c( '_bg' ) ) ),
	m			<- merge(
		get_pheno_covar_wtiss( phen, tiss1, pv, internal=( type %in% c( 'int', 'bin' ) ), allchrs=TRUE, return_ID=TRUE, return_bg=F ),
		get_pheno_covar_wtiss( phen, tiss2, pv, internal=( type %in% c( 'int', 'bin' ) ), allchrs=TRUE, return_ID=TRUE, return_bg=F )[,c('FID','IID',paste0(tiss2,'_',1:22))]
	)
	m	<- m[,-which( colnames(m) %in% c( 'FID', 'IID' ) )]

	if( perm.i > 0 ){
		set.seed( perm.i )
		perm	<- sort.list( runif(nrow(m)) )
		m[,paste0(tiss1,'_',1:22)]	<- m[perm,paste0(tiss1,'_',1:22)]
		m[,paste0(tiss2,'_',1:22)]	<- m[perm,paste0(tiss2,'_',1:22)]
	} 

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
		paste0('X',1:22),
		paste0(tiss1,'_',1:22),
		paste0(tiss2,'_',1:22),
		'age_when_attended_assessment_centre_f21003_0_0',
		'sex_f31_0_0',
		ifelse(withBMI,'BMI',''),
		paste("genetic_principal_components_f22009_0_",1:nPC,sep="")),
	collapse="+")) 

	X	<- matrix( 1:22, 22, 22 )
	is	<- t(X)[lower.tri(X,diag=F)]
	js	<- X   [lower.tri(X,diag=F)]
  f0 <- paste0( c( f0, 
		paste0('X',is,':','X',js),
		paste0('X',is,':',tiss1, '_', js),
		paste0('X',js,':',tiss1, '_', is),
		paste0('X',is,':',tiss2, '_', js),
		paste0('X',js,':',tiss2, '_', is)
	), collapse="+") 

	#if( opt == '_bg' )
	#	stop()

	interxns	<- c(
		paste0( tiss1, '_', is,':',tiss2, '_', js),
		paste0( tiss1, '_', js,':',tiss2, '_', is)
	)
	f1		<- paste0( c( f0, interxns ), collapse="+") 

	if( binary ){
		lm0	<- glm(as.formula(f0), data=m,family=binomial)
		lm1	<- glm(as.formula(f1), data=m,family=binomial)
		pval_inter <- anova( lm1, lm0, test='Chisq' )$Pr[2] 
	} else {
		lm0	<- lm( as.formula(f0), data=m)
		lm1	<- lm( as.formula(f1), data=m)
		pval_inter <- anova( lm1, lm0 )$Pr[2] 
	} 
	fit_summ	<- summary(lm1)$coef 
	save(  fit_summ, pval_inter, file=savefile )
	rm( m, fit_summ, pval_inter, lm0, lm1 ); gc()

})
	sink()
})
