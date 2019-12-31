library("ROCR")
regress <- function(m, vars=c(), interaction=c(), withBMI, outcome='Pheno', binary, nPC ) {
  f <- paste(outcome, "~", paste(c(
		interaction, vars,
		'age_when_attended_assessment_centre_f21003_0_0', 'Male', ifelse(withBMI,'BMI',''),
		paste("genetic_principal_components_f22009_0_",1:nPC,sep="")), collapse="+"))
		
		#print( 'form' )
		print( f )
		print( as.formula( f ) )
		print( 'xxxxx' )
	if( binary ){
		m[,outcome]	<-  as.numeric(m[,outcome])-1
		return( glm(as.formula(f), data=m,family=binomial) )
	} else {
		return( lm(	as.formula(f), data=m) )
	}
}

single_var_exp <- function(m, phenotype, feature, withBMI, nPC, binary){
  v <- var(m[,feature])
	s <- regress(m, vars=c(feature), withBMI=withBMI, nPC=nPC, binary=binary)
	if( binary ){
		pred	<- prediction(predict(s,type=c("response")), m$Pheno)
		return( as.numeric( performance(pred, "auc")@y.values ) )
	} else {
		b <- summary(s)$coefficients[feature,1]
		return( b^2*v / var(m$Pheno)*100 )
	}
}
