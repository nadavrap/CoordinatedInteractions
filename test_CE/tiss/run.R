rm( list=ls() )
source('~/coordination_fresh/functions.R')
source('setup.R') 

eo_res					<- read.csv( '~/coordination_fresh/figs/table_s3.csv', head=T )
opt_pvs					<- eo_res[,3]
names(opt_pvs)	<- paste0( eo_res[,1], '_', eo_res[,2] )

hitphens	<- as.character( eo_res[ eo_res[,'FDR'] < .1, 1 ] )
hitphens	<- phens[ nicephens[phens] %in% hitphens ]

for( perm.i in 0:5 )
	for( opt in c( '_nondiag', '_bg' ) )
		for( type in c( 'int', 'bin', 'ext' ) )
			for( phen in sample(hitphens) )
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

	savefile	<- paste0('Rdata/allpairs_', pv, '_', phen, '_', type, '_', tiss, '_', perm.i, opt, '.Rdata')
	sinkfile	<- paste0( 'Rout/allpairs_', pv, '_', phen, '_', type, '_', tiss, '_', perm.i, opt, '.Rout')
	if( file.exists( savefile ) | file.exists( sinkfile ) ) next
	print( sinkfile )
	sink( sinkfile ) 

	withBMI	<- ( phen != 'BMI' )
	binary	<- phen %in% binphens

	m			<- get_pheno_covar_wtiss( phen, tiss, pv, internal=( type %in% c( 'int', 'bin' ) ), allchrs=TRUE, return_bg=( opt %in% c( '_bg' ) ) )
	if( perm.i > 0 ){
		set.seed( perm.i )
		perm	<- sort.list( runif(nrow(m)) )
		mycols<- intersect( paste0(tiss,'_',1:22), colnames(m) )
		m[,mycols]	<- m[perm,mycols]
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
		paste0(tiss,'_',1:22),
		'age_when_attended_assessment_centre_f21003_0_0',
		'sex_f31_0_0',
		ifelse(withBMI,'BMI',''),
		paste("genetic_principal_components_f22009_0_",1:nPC,sep="")),
	collapse="+")) 

	X	<- matrix( 1:22, 22, 22 )
	is	<- t(X)[lower.tri(X,diag=F)]
	js	<- X   [lower.tri(X,diag=F)]
  f0 <- paste0( c( f0, paste0('X',is,':','X',js) ), collapse="+") 

	if( opt == '_bg' )
		f0		<- paste0( c( f0, 
			paste0('BG_',1:22),
			paste0( 'X', is,':BG_', js),
			paste0( 'X', js,':BG_', is)
		), collapse="+") 

	interxns	<- c(
		paste0( 'X', is,':',tiss, '_', js),
		paste0( 'X', js,':',tiss, '_', is)
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

	subs	<- c( 
		paste0('X',1:22),
		paste0(tiss,'_',1:22),
		paste0('X',is,':','X',js), 
		interxns
	) 
	subs		<- intersect( subs, rownames(summary(lm1)$coef) ) 
	covbet	<- vcov(lm1)[subs,subs]

	save(  fit_summ, pval_inter, covbet, file=savefile )
	rm( m, fit_summ, pval_inter, covbet, lm0, lm1 ); gc()

	sink()
})
