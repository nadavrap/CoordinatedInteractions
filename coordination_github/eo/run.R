rm( list=ls() )
source('~/coordination_fresh/functions.R')
source('setup.R') 

for( perm.i in 0:100 ) 
	for( type in c( 'ext', 'int', 'bin' ) )
		for( phen in sample(phens) )
			for( pv in sample(pvs) )
try({
	if( (type == 'ext' & !phen %in% extphens) |
			(type == 'bin' & !phen %in% binphens) |
			(type == 'int' &  phen %in% binphens) ) next 

	if( type == 'ext' & phen %in% extphens1 & pv != '1.0' ) next

	savefile	<- paste0('Rdata/allpairs_', pv, '_', phen, '_', type, '_', perm.i, '.Rdata')
	sinkfile	<- paste0( 'Rout/allpairs_', pv, '_', phen, '_', type, '_', perm.i, '.Rout')
	if( file.exists( savefile ) | file.exists( sinkfile ) ) next
	print( sinkfile )
	sink( sinkfile )

	withBMI	<- ( phen != 'BMI' )
	binary	<- phen %in% binphens

	m			<- get_pheno_and_covar(phen, pv, internal=( type %in% c( 'int', 'bin' ) ), seed=0, allchrs=TRUE ) 
	if( perm.i > 0 ){
		set.seed( perm.i )
		perm	<- sort.list( runif(nrow(m)) )
		m[,paste0('X',1:22)]	<- m[perm,paste0('X',1:22)]
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
		'age_when_attended_assessment_centre_f21003_0_0',
		'sex_f31_0_0',
		ifelse(withBMI,'BMI',''),
		paste("genetic_principal_components_f22009_0_",1:nPC,sep="")),
	collapse="+")) 

	X		<- matrix( 1:22, 22, 22 )
	is	<- t(X)[lower.tri(X,diag=F)]
	js	<- X   [lower.tri(X,diag=F)]
  f1	<- paste0( c( f0, paste0('X',is,':','X',js) ), collapse="+") 

	if( binary ){
		lm0	<- glm(as.formula(f0), data=m,family=binomial)
		lm1	<- glm(as.formula(f1), data=m,family=binomial)
		pval_inter <- anova( lm1, lm0, test='Chisq' )$Pr[2] 
	} else {
		lm0	<- lm( as.formula(f0), data=m)
		lm1	<- lm( as.formula(f1), data=m)
		pval_inter <- anova( lm1, lm0 )$Pr[2] 
	} 

	subs	<- intersect( paste0('X',is,':','X',js), rownames(summary(lm1)$coef) ) 
	fit_summ	<- summary(lm1)$coef 
	covbet	<- vcov(lm1)[subs,subs]

	save(  fit_summ, pval_inter, covbet, file=savefile )
	rm( m, fit_summ, pval_inter, covbet, lm0, lm1 ); gc()

	sink()
})
