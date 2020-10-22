types	<- c( 'int40', 'bin10', 'ext' )
nicetypes	<- c( 'Intern Cont.', 'Intern Bin', 'Extern' )
names(nicetypes)	<- types

pvs		<- c( '1e-08', '1e-07', '1e-06', '1e-05', '0.0001', '0.001', '0.01', '0.1', '1.0' ) #, '0.005', '0.05', '0.5'
seeds	<- as.character( c( 0, round( 7 + (2:99)^3/42 ), 192 ) )

extphens<- c( 'Triglycerides', 'cov_EDU_YEARS', 'BMI', 'Height', 'LDL_direct', 'T2D', 'disease_ASTHMA_DIAGNOSED', 'disease_CARDIOVASCULAR' )

binphens<- c('T2D', 'disease_ASTHMA_DIAGNOSED', 'disease_CARDIOVASCULAR', 'disease_HI_CHOL_SELF_REP', 'disease_ALLERGY_ECZEMA_DIAGNOSED')
qphens	<- c('blood_EOSINOPHIL_COUNT', 'blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT', 
			'blood_LYMPHOCYTE_COUNT', 'blood_MEAN_CORPUSCULAR_HEMOGLOBIN', 'blood_MEAN_PLATELET_VOL', 'blood_MEAN_SPHERED_CELL_VOL',
			'blood_MONOCYTE_COUNT', 'blood_PLATELET_COUNT', 'blood_PLATELET_DISTRIB_WIDTH', 'blood_RBC_DISTRIB_WIDTH', 'blood_RED_COUNT',
			'blood_WHITE_COUNT', 'bmd_HEEL_TSCOREz', 'BMI', 'cov_EDU_YEARS', 'impedance_BASAL_METABOLIC_RATEz', 'lung_FEV1FVCzSMOKE',
			'Height', 'Triglycerides', 'LDL_direct', 'Glucose' )
phens		<- c( qphens, binphens)
nicephens	<- c('Eosinophil #', 'Reticulocyte #', 'Lymphocyte #', 'Corp. Hemoglobin', 'Platelet Vol.', 'Sphered Cell Vol',
			'Monocyte #', 'Platelet #', 'Platelet Distn Width', 'RBC Distn Width', 'RBC #',
			'WBC #', 'Heel Bone Mineral Dens.', 'BMI', 'Edu Years', 'Basal Metabolic Rate', 'FEV1FVCzSMOKE',
			'Height', 'Triglycerides', 'LDL', 'Glucose', 'T2D', 'Asthma', '`CARDIOVASCULAR`','High Cholest.', 'Eczema')


cols0    <- c( '#e6beff','#fabebe', '#008080',  '#808080','#800000', '#911eb4', '#9a6324', '#ffd8b1', '#f58231', '#fffac8', '#4363d8', '#aaffc3',  '#46f0f0', '#bcf60c', '#ffe119',  '#e6194b', '#f032e6','#3cb44b','#808000','#000075', '#66ccff' )
cols		<- c( cols0, cols0[c(6,16,13,11,12)] )

nicetiss<- c(paste0('Franke.', c('Blood.Cells', 'Adipocytes', 'Brain', 'Hippocampus', 'Liver', 'Muscles', 'Pancreas')), c('Adipose', 'Brain', 'Hippocampus', 'Liver', 'Pancreas', 'Skeletal_Muscle'))
TISSUES <- c( 'Franke.11', 'Franke.120', 'Franke.14', 'Franke.46', 'Franke.69', 'Franke.84', 'Franke.89', 'Adipose', 'Brain', 'Brain_Hippocampus', 'Liver', 'Pancreas', 'Skeletal_Muscle')

S				<- length(seeds)
npv			<- length(pvs)
P				<- length(phens)
ntiss   <- length(TISSUES)


names(cols)			<- phens
names(nicephens)<- phens
names(nicetiss)	<- TISSUES


cbind( nicephens, nicephens )
cbind( TISSUES, nicetiss )

getCovars <- function(fname) {
  cachefname <- paste0(fname, '.rds')
  if(file.exists(cachefname)) return(readRDS(cachefname))
  covs <- read.table(fname, header = TRUE, sep='')
  covs$Male <- as.logical(covs$sex_f31_0_0 == 'Male')
  covs$sex_f31_0_0 <- NULL
  # Set genotype batch and assemsment ceters as factor
  covs$genotype_measurement_batch_f22000_0_0 <- as.factor(covs$genotype_measurement_batch_f22000_0_0)
  covs$uk_biobank_assessment_centre_f54_0_0 <- as.factor(covs$uk_biobank_assessment_centre_f54_0_0)
  # Take median of BMI in case of multiple measurments
  covs$BMI <- apply(covs[,grep('body_mass_index_bmi_f21001_', colnames(covs))], 1, function(x) median(x, na.rm = TRUE))
  covs[,grep('body_mass_index_bmi_f21001_', colnames(covs))] <- NULL
  saveRDS(covs, file=cachefname)
  covs
}
