ukbdir	<- '/u/project/sriram/nadavrap/UKBB/'

types	<- c( 'int', 'bin', 'ext' )
nicetypes	<- c( 'Intern Quant', 'Intern Bin', 'Extern' )
names(nicetypes)	<- types

pvs		<- c( '1e-08', '1e-07', '1e-06', '1e-05', '0.0001', '0.001', '0.01', '0.1', '1.0' ) #, '0.005', '0.05', '0.5'
#seeds	<- as.character( c( 0, round( 7 + (2:99)^3/42 ), 192 ) )

nperm	<- 2

extphens0	<- c( 'Triglycerides', 'cov_EDU_YEARS', 'BMI', 'Height', 'LDL_direct', 'T2D', 'disease_ASTHMA_DIAGNOSED', 'disease_CARDIOVASCULAR' )
extphens1	<- c( 'blood_PLATELET_COUNT', 'blood_MEAN_CORPUSCULAR_HEMOGLOBIN', 'blood_WHITE_COUNT', 'blood_MONOCYTE_COUNT', 'blood_MEAN_PLATELET_VOL', 'blood_RED_COUNT', 'blood_LYMPHOCYTE_COUNT', 'blood_EOSINOPHIL_COUNT', 'blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT' )
extphens	<- c( extphens0, extphens1 )

binphens<- c('T2D', 'disease_ASTHMA_DIAGNOSED', 'disease_CARDIOVASCULAR', 'disease_HI_CHOL_SELF_REP', 'disease_ALLERGY_ECZEMA_DIAGNOSED')
qphens	<- c('blood_EOSINOPHIL_COUNT', 'blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT', 
			'blood_LYMPHOCYTE_COUNT', 'blood_MEAN_CORPUSCULAR_HEMOGLOBIN', 'blood_MEAN_PLATELET_VOL', 'blood_MEAN_SPHERED_CELL_VOL',
			'blood_MONOCYTE_COUNT', 'blood_PLATELET_COUNT', 'blood_PLATELET_DISTRIB_WIDTH', 'blood_RBC_DISTRIB_WIDTH', 'blood_RED_COUNT',
			'blood_WHITE_COUNT', 'bmd_HEEL_TSCOREz', 'BMI', 'cov_EDU_YEARS', 'impedance_BASAL_METABOLIC_RATEz', 'lung_FEV1FVCzSMOKE',
			'Height', 'Triglycerides', 'LDL_direct', 'Glucose' )
phens		<- c( qphens, binphens)
nicephens	<- c('Eosinophil #', 'Reticulocyte #', 'Lymphocyte #', 'Corp. Hemoglobin', 'Platelet Vol.', 'Sphered Cell Vol',
			'Monocyte #', 'Platelet #', 'Platelet Distn Width', 'RBC Distn Width', 'RBC #',
			'WBC #', 'Bone Mineral Dens.', 'BMI', 'Edu Years', 'Basal Metabolic Rate', 'FEV/FVC',
			'Height', 'Triglycerides', 'LDL', 'Glucose', 'T2D', 'Asthma', 'Cardiovascular','High Cholest.', 'Eczema')



## #cols0    <- c( '#e6beff','#fabebe', '#008080',  '#808080','#800000', '#911eb4', '#9a6324', '#ffd8b1', '#f58231', '#fffac8', '#4363d8', '#aaffc3',  '#46f0f0', '#bcf60c', '#ffe119',  '#e6194b', '#f032e6','#3cb44b','#808000','#000075', '#66ccff' )
## cols0    <- c( '#e6beff','#fabebe', '#008080',  '#808080','#800000', '#911eb4', '#9a6324', '#ffd8b1', '#f58231', 1, '#4363d8', '#aaffc3',  '#46f0f0', '#bcf60c', '#ffe119',  '#e6194b', '#f032e6','#3cb44b','#808000','#000075', '#66ccff' )
## cols		<- c( cols0, cols0[c(11,16,13,7,12)] )

nicetiss<- c(paste0('Franke.', c('Blood.Cells', 'Adipocytes', 'Brain', 'Hippocampus', 'Liver', 'Muscles', 'Pancreas')), c('Adipose', 'Brain', 'Hippocampus', 'Liver', 'Pancreas', 'Skeletal_Muscle'))
tissues <- c( 'Franke.11', 'Franke.120', 'Franke.14', 'Franke.46', 'Franke.69', 'Franke.84', 'Franke.89', 'Adipose', 'Brain', 'Brain_Hippocampus', 'Liver', 'Pancreas', 'Skeletal_Muscle')

#S				<- length(seeds)
npv			<- length(pvs)
P				<- length(phens)
ntiss   <- length(tissues)

# #names(cols)			<- phens
names(nicephens)<- phens
names(nicetiss)	<- tissues

cbind( nicephens, nicephens )
cbind( tissues, nicetiss )

nPC			<- 40
