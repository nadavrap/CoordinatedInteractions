#!R
# For example:
# Rscript ukb_extract_phenotype.R --pheno_name BMI --out ./BMI/BMI.pheno
# Notice, not all phenotypes are implemented, as there are some that have some logic for extraction.
library(ukbtools)
suppressPackageStartupMessages(library("argparse"))

SCRIPTS_DIR <- '/u/project/sriram/nadavrap/UKBB/CoordinatedInteractions/src/'
PATH <- system(paste0('python3 ', SCRIPTS_DIR, 'parameters.py WD') ,intern=TRUE)
FNAME <- system(paste0('python3 ', SCRIPTS_DIR, 'parameters.py UKB_FNAME') ,intern=TRUE)
FNAMES <- c('ukb21608', 'ukb27646', 'ukb25904')

pheno2ID = list(T2D=2443,
                blood_EOSINOPHIL_COUNT=30150,
                blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT=30300,
                blood_LYMPHOCYTE_COUNT=30120,
                blood_MEAN_CORPUSCULAR_HEMOGLOBIN=30050,
                blood_MEAN_PLATELET_VOL=30100,
                blood_MEAN_SPHERED_CELL_VOL=30270,
                blood_MONOCYTE_COUNT=30130,
                blood_PLATELET_COUNT=30080,
                blood_PLATELET_DISTRIB_WIDTH=30110,
                blood_RBC_DISTRIB_WIDTH=30070,
                blood_RED_COUNT=30010,
                blood_WHITE_COUNT=30000,
                bmd_HEEL_TSCOREz=c(4106, 4125), # Take sum
                #body_BALDING1=2395,
                #body_BALDING4=2395,
                body_BMIz=21001,
                body_HEIGHTz=50,
                body_WHRadjBMIz=c(48,49), # Take ratio
                bp_DIASTOLICadjMEDz=4079,
                bp_SYSTOLICadjMEDz=4080,
                #cov_EDU_COLLEGE=6138~6138-0.5,
                #cov_EDU_YEARS=6138~6138-0.5,
                #cov_SMOKING_STATUS=20116,
                #disease_AID_ALL=20002~20002-0.28,
                #disease_AID_SURE=20002~20002-0.28,
                disease_ALLERGY_ECZEMA_DIAGNOSED=6152,
                disease_ASTHMA_DIAGNOSED=6152,
                disease_CVD=6150,
                disease_CARDIOVASCULAR=20002,
                #disease_DERMATOLOGY=20002~20002-0.28,
                #disease_HI_CHOL_SELF_REP=20002~20002-0.28,
                #disease_HYPOTHYROIDISM_SELF_REP=20002~20002-0.28,
                #disease_RESPIRATORY_ENT=20002~20002-0.28,
                #disease_T2D=2443 + 2976,
                impedance_BASAL_METABOLIC_RATEz=23105,
                lung_FEV1FVCzSMOKE=c(3062,3063), # Ratio
                lung_FVCzSMOKE=3062,
                mental_NEUROTICISM=20127,
                #other_MORNINGPERSON=1180,
                #pigment_HAIR=1747,
                pigment_HAIR_blonde=1747,
                pigment_HAIR_darkbrown=1747,
                #pigment_SKIN=1717,
                #pigment_SUNBURN=1737,
                #pigment_TANNING=1727,
                repro_MENARCHE_AGE=2714,
                repro_MENOPAUSE_AGE=3581,
                repro_NumberChildrenEverBorn_Males=2405,
                repro_NumberChildrenEverBorn_Females=2734,
                cov_EDU_YEARS=6138,
                disease_HI_CHOL_SELF_REP=20002, # 1473
		Glucose=30740,
		LDL_direct=30780,
		Triglycerides=30870
                )


# Takes about 25 minutes
get_data <- function(pref) {
	 fname=paste0(pref, "_data.rds")
	 if (!file.exists(fname)) {
	 my_ukb_data <- ukb_df(pref, path=PATH)
	 # Create a field-to-descriptive name key, as dataframe or named lookup vector.
	 my_ukb_key <- ukb_df_field(pref, path = PATH)
	# save(my_ukb_data, file = fname)
	saveRDS(my_ukb_data, file = fname)
	 } else {
	 my_ukb_data <- readRDS(fname)
	 }
	 my_ukb_data
}
var_map <- function(pref) {
	fname <- paste0(pref,'_variables_map.rds')
	if (!file.exists(fname)) {
	    my_ukb_key <- ukb_df_field(pref, path = PATH)
	    if (my_ukb_key[1,1] =='eid') {
	    	    my_ukb_key <- my_ukb_key[2:nrow(my_ukb_key),]
	    }
	    saveRDS(my_ukb_key, fname)
	} else {
	    my_ukb_key <- readRDS(fname)
	}
	my_ukb_key
}
get_cols_idxs <- function(my_ukb_key, data_id) {
    suppressWarnings({
	# Supress, as the first column is 'eid', but the others are integers.
	vars <- my_ukb_key$col.name[as.integer(my_ukb_key$field.showcase) %in% data_id]
    })
    if (length(vars) == 0) {
        warning('No columns found matching the given data id')
        return(NA)
    }
    vars
}

#' Take median of each data id, and then sum over data_ ids
median_and_agg <- function(data_ids, my_ukb_data, my_ukb_key, aggfunc='+') {
    print(paste('Aggregate using median_and_agg and', aggfunc))
    vals <- lapply(data_ids, function(data_id) {
        d  <- my_ukb_data[,get_cols_idxs(my_ukb_key, data_id)]
        agg_vals <- apply(sapply(d, as.numeric), 1, median, na.rm=TRUE)
    })
    if(length(vals) == 1) {
        return(vals[[1]])
    }
    res <- Reduce(function(x,y) mapply(aggfunc,x,y),vals)
    print('Summary of values:')
    print(summary(res))
    if(any(is.infinite(res))) {
        print('Found Inf values. Cannot continue')
        print(head(which(is.infinite(res))))
        for (v in vals) {
            print(v[head(which(is.infinite(res)))])
        }
        quit()
    }
    
    res
}

write_pheno <- function(data_id, agg_func=NA, field_out_name, fname, negativeAsNA=FALSE) {
    for (pref in FNAMES) {
        my_ukb_key <- var_map(pref)
        cols_idx <- get_cols_idxs(my_ukb_key, data_id)
        print(paste('Column index for', paste(data_id), 'is/are:', paste(cols_idx), 'in', pref))
        if (! all(is.na(cols_idx))) {
            break
        }
    }
    my_ukb_data <- get_data(pref)
    print(paste('Retrieved data for', nrow(my_ukb_data), 'samples.'))
    
    my_ukb_data <- my_ukb_data[,c('eid', cols_idx)]
    print(paste('Number of columns:', ncol(my_ukb_data)))
    
        if (negativeAsNA) {
            print(paste('Setting negative values as NA.',
                        'Useful as for some numeric data types',
                        '-1 represents "Do not know" and -3 represents "Prefer not to answer"'))
           my_ukb_data[my_ukb_data<0] <- NA
        }
        if (is.function(agg_func)) {
            agg_vals <- agg_func(data_id, my_ukb_data, my_ukb_key)
            my_ukb_data[field_out_name] <- agg_vals
            my_ukb_data <- my_ukb_data[, c('eid', field_out_name)]
            data_id <- field_out_name
        }
        print(paste('Print phenotype:', data_id))
	ukb_gen_write_plink(my_ukb_data, path=fname, ukb.variables=data_id)
        print(paste('Wrote', length(data_id), 'features to', fname))
}

#' Aggregation function that keep Yes/No values, and assign others to NA
#'
#' For some data fields in UKBB like T2D (2443), the values are represented as factor with names
#' Yes, No, 'Prefer not to answer' and 'Do not know'
#' This aggregation function, get a list of values, and make a decision for Yes, No or NA
#' @param l list of values
#' @default In case of conflict, this value take over.
case_control <- function(data_id, my_ukb_data, my_ukb_key, cases='Yes', missing=c('Prefer not to answer', 'Do not know'), all_missing_as_control=FALSE) {
    print(paste('Aggregate using case_control function where', paste(cases), 'values are set as cases'))
    my_ukb_data  <- my_ukb_data[,get_cols_idxs(my_ukb_key, data_id)]
    print('Values found:')
    print(summary(as.factor(my_ukb_data[,1])))
    res <- apply(my_ukb_data, 1, function(l) {
        if (any(cases %in% l)) return(2)
        if (any(missing %in% l)) return(-9)
        if (! all_missing_as_control) {
            if (all(is.na(l))) return(-9)
        }
        return(1)})
    print(summary(as.factor(res)))
    res
}

eductaion_mapping <- function(data_id, my_ukb_data, my_ukb_key) {
    # College or University degree = 20 years;
    # A levels/AS levels or equivalent = 13 years;
    # O levels/GCSEs or equivalent = 10 years;
    # CSEs or equivalent = 10 years;
    # NVQ or HND or HNC or equivalent = 19 years;
    # Other professional qualifications eg: nursing, teaching = 15 years;
    # None of the above = 7 years;
    # Prefer not to answer = NA
    # (Phenotype construction based on Okbay et al. 2016 Nature, Tables S1.2 and S1.14)
    my_ukb_data  <- my_ukb_data[,get_cols_idxs(my_ukb_key, data_id)]
    m <- c("College or University degree"=20, 'A levels/AS levels or equivalent' = 13, 'O levels/GCSEs or equivalent' = 10,
      'CSEs or equivalent' = 10, 'NVQ or HND or HNC or equivalent' = 19,
      'Other professional qualifications eg: nursing, teaching' = 15,
      'None of the above' = 7, 'Prefer not to answer' = NA)
    res <- apply(my_ukb_data, 1, function(x) max(m[x], na.rm=TRUE))
    res[is.infinite(res)] <- NA
    res
}

get_suitable_agg_function <- function(pheno_name) {
    agg_func <- median_and_agg
    if (tolower(pheno_name) == 't2d') {
        agg_func <- function(data_ids, my_ukb_data, my_ukb_key) case_control(data_ids, my_ukb_data, my_ukb_key, cases='Yes')
    } else if (pheno_name %in% c('body_WHRadjBMIz', 'lung_FEV1FVCzSMOKE')) { # Need to take the ratio
        agg_func <- function(data_ids, my_ukb_data, my_ukb_key) median_and_agg(data_ids, my_ukb_data, my_ukb_key, aggfunc='/')
    } else if (pheno_name == 'bmd_HEEL_TSCOREz') {
        agg_func <- median_and_agg
    } else if (pheno_name == 'disease_ALLERGY_ECZEMA_DIAGNOSED') {
        agg_func <- function(data_ids, my_ukb_data, my_ukb_key) case_control(data_ids, my_ukb_data, my_ukb_key, cases='Hayfever, allergic rhinitis or eczema')
    } else if (pheno_name == 'disease_ASTHMA_DIAGNOSED') {
        agg_func <- function(data_ids, my_ukb_data, my_ukb_key) case_control(data_ids, my_ukb_data, my_ukb_key, cases='Asthma')
    } else if (pheno_name == 'pigment_HAIR_blonde') {
        agg_func <- function(data_ids, my_ukb_data, my_ukb_key) case_control(data_ids, my_ukb_data, my_ukb_key, cases='Blonde')
    } else if (pheno_name == 'pigment_HAIR_darkbrown') {
        agg_func <- function(data_ids, my_ukb_data, my_ukb_key) case_control(data_ids, my_ukb_data, my_ukb_key, cases='Dark brown')
    } else if (pheno_name == 'disease_HI_CHOL_SELF_REP') {
        agg_func <- function(data_ids, my_ukb_data, my_ukb_key) case_control(data_ids, my_ukb_data, my_ukb_key, cases='1473', all_missing_as_control=TRUE)
    } else if (pheno_name == 'disease_CARDIOVASCULAR') {
        # 1065 (hypertension); 1066 (heart/cardiac problem); 1067 (peripheral vascular disease); 1068 (venous thromboembolic disease); 1081 (stroke); 1082 (transient ischaemic attack (tia)); 1083 (subdural haemorrhage/haematoma); 1425 (cerebral aneurysm); 1473 (high cholesterol); 1493 (other venous/lymphatic disease)
        agg_func <- function(data_ids, my_ukb_data, my_ukb_key) case_control(data_ids, my_ukb_data, my_ukb_key, cases=c('1065', '1066', '1067', '1068', '1081', '1082', '1083', '1425', '1473', '1493'), all_missing_as_control=TRUE)
    } else if (pheno_name == 'disease_CVD') {
        agg_func <- function(data_ids, my_ukb_data, my_ukb_key) case_control(data_ids, my_ukb_data, my_ukb_key, cases=c("Heart attack","Angina","Stroke","High blood pressure"))
    } else if (pheno_name == 'cov_EDU_YEARS') {
        agg_func <- eductaion_mapping
    } else if (pheno_name == 'cov_EDU_COLLEGE') {
        agg_func <- function(data_ids, my_ukb_data, my_ukb_key) case_control(data_ids, my_ukb_data, my_ukb_key, cases='College or University degree')
    }
    agg_func
}

# create parser object
parser <- ArgumentParser(description=paste('Extract a single phenotype. In case multi columns associated, ',
                            'decide whehter aggregate or not. Notice that not all phenotypes are implemented yet.'))
# Phenotype codes
parser$add_argument('-i', '--pheno_id', type="integer", help='Phenotype code (e.g., 50, 22001, 22006)')
parser$add_argument('-n', '--pheno_name', help=paste('Phenotype name - implemented only for', paste(names(pheno2ID), collapse=' ')), default=NULL)
parser$add_argument('-o', '--out', help='Output filename')
parser$add_argument('-f', '--data_file', help='UKB Data file prefix', default=FNAME)
#parser$add_argument('-a', '--agg_func', default=NULL, help='Aggregation function. Default no aggregation, but print all columns.')

args <- parser$parse_args()

agg_func <- NA

pheno_id = args$pheno_id
if (!is.null(args$pheno_name)) {
    pheno_id = pheno2ID[[args$pheno_name]]
    print(paste('Phenotype named', args$pheno_name, 'was mapped to id:', pheno_id))
}

agg_func <- get_suitable_agg_function(args$pheno_name)

if (args$pheno_name %in% c('repro_NumberChildrenEverBorn_Males', 'repro_NumberChildrenEverBorn_Females')) {
    negativeAsNA <- TRUE
} else {
    negativeAsNA <- FALSE
}

#write_pheno(my_ukb_data, pheno_id, agg_func, args$pheno_name, args$out, negativeAsNA=negativeAsNA)
write_pheno(pheno_id, agg_func, args$pheno_name, args$out, negativeAsNA=negativeAsNA)
