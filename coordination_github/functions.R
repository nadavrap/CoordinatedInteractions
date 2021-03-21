library("ROCR")

getCovars <- function(fname) {
  cachefname <- paste0(fname, '.rds')
  if(file.exists(cachefname)) return(readRDS(cachefname))
  covs <- read.table(fname, header = TRUE, sep='')
	#print( dim( covs ) ) #[1] 342815    171 
	#print( head( round( covs,2 ) ) )
	
	# recode sex
  # this was done with current eo, but SHOULD NOT ... now renamed to .rds.eo.old
  #covs$Male <- as.logical(covs$sex_f31_0_0 == 'Male')
  #covs$sex_f31_0_0 <- NULL

  # Take median of BMI in case of multiple measurments
  covs$BMI <- apply(covs[,grep('body_mass_index_bmi_f21001_', colnames(covs)),drop=F], 1, function(x) median(x, na.rm = TRUE))
  covs[,grep('body_mass_index_bmi_f21001_', colnames(covs))] <- NULL 

  # ## # Set genotype batch and assemsment ceters as factor
  # ## #covs$genotype_measurement_batch_f22000_0_0 <- as.factor(covs$genotype_measurement_batch_f22000_0_0)
  # ## #covs$uk_biobank_assessment_centre_f54_0_0 <- as.factor(covs$uk_biobank_assessment_centre_f54_0_0)

	#### check factors
	#print( 'cc' )
	#print( table(sample( unlist(covs[,grep('genotype_measurement_batch_f22000_0_0', colnames(covs))]), 1e5 ) ))
	#print( quantile( colSums( covs[,grep('genotype_measurement_batch_f22000_0_0', colnames(covs))], na.rm=T ) ) ) 
	#print( table(rowSums( covs[,grep('genotype_measurement_batch_f22000_0_0', colnames(covs))], na.rm=T ) ))
	#### TODO: is it ok some are 0s? :: it is the right size if one level was dropped. I see a similar thing below
	##       0      1
	##   3415 339400

	#print( 'dd' )
	#print( table(sample( unlist(covs[,grep('uk_biobank_assessment_centre_f54_0_0', colnames(covs))]), 1e5 ) ))
	#print( table(rowSums( covs[,grep('uk_biobank_assessment_centre_f54_0_0', colnames(covs))] ) ))
	#print(       colSums( covs[,grep('uk_biobank_assessment_centre_f54_0_0', colnames(covs))] ) )

	bads	<- which( apply( covs[,grep('genotype_measurement_batch_f22000_0_0', colnames(covs))], 1, function(x) any(x==-9) ) )
	print( length(bads) )
	if( length(bads) > 0 ) stop( length(bads) )
	#print( head( bads ) )
	#covs[bads,grep('genotype_measurement_batch_f22000_0_0', colnames(covs))]	<- NA

	#print( head( covs[,grep('genotype_measurement_batch_f22000_0_0', colnames(covs))] ) )
	#print( table( as.character(unlist(covs[,grep('genotype_measurement_batch_f22000_0_0', colnames(covs))])) %in% 0:1, na.rm=T ))
	#print( table( covs[,grep('genotype_measurement_batch_f22000_0_0', colnames(covs))], exclude=NULL ) )

  saveRDS(covs, file=cachefname)
  covs
}

tidy	<- function(x)
	ifelse( x > .01, round( x, 2 ), format( x, digit=2, scientific=T ) ) 
