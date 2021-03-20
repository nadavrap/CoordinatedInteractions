#!R
# Example of execution:
# ukb_extract_phenotypes.R 20002 41142 -o=./data/ukbb_20002_41142_filtered.txt -p -r -s -c
# Takes about 20 minutes

library(ukbtools) # See examples here: https://kenhanscombe.github.io/ukbtools/
suppressPackageStartupMessages(library("argparse"))
SCRIPTS_DIR <- '/u/project/sriram/nadavrap/UKBB/CoordinatedInteractions/src/'
WD <- system(paste0('python3 ', SCRIPTS_DIR, 'parameters.py WD') ,intern=TRUE)
FNAME <- system(paste0('python3 ', SCRIPTS_DIR, 'parameters.py UKB_FNAME') ,intern=TRUE)
#FNAME='ukb27646'
#FNAME="ukb21608"
#FNAME='ukb25904'

# Takes about 25 minutes
get_data <- function() {
	 fname=paste0(WD, FNAME, ".rds")
	 if (!file.exists(fname)) {
	 my_ukb_data <- ukb_df(FNAME, path=WD)
	 # Create a field-to-descriptive name key, as dataframe or named lookup vector.
	 my_ukb_key <- ukb_df_field(FNAME, path = WD)
	# save(my_ukb_data, file = fname)
	saveRDS(my_ukb_data, file = fname)
	 } else {
	 my_ukb_data <- readRDS(fname)
	 }
	 my_ukb_data
}
var_map <- function() {
	fname <- paste0(WD, FNAME, '_variables_map.rds')
	if (!file.exists(fname)) {
	    my_ukb_key <- ukb_df_field(FNAME, path = WD)
	    if (my_ukb_key[1,1] =='eid') {
	    	    my_ukb_key <- my_ukb_key[2:nrow(my_ukb_key),]
	    }
	    saveRDS(my_ukb_key, fname)
	} else {
	    my_ukb_key <- readRDS(fname)
	}
	my_ukb_key
}

write_covar_and_pheno <- function(my_ukb_data, data_codes, fname) {
	my_ukb_key <- var_map()
	vars <- my_ukb_key$col.name[as.integer(my_ukb_key$field.showcase) %in% data_codes]

	# Age at Initial assessment visit
	more_cov <- my_ukb_key$col.name[my_ukb_key$field.html == '21003-0.0']
	# Sex
	more_cov <- c(more_cov, my_ukb_key$col.name[my_ukb_key$field.showcase == '31'])
	# Center of assessment
	# 'uk_biobank_assessment_centre_f54_0_0'
	more_cov <- c(more_cov, my_ukb_key$col.name[my_ukb_key$field.showcase == '54'])
	#pcs_col_names <- pcs_col_names[!is.na(pcs_col_names)]
	ukb_gen_write_plink(my_ukb_data, path=fname, ukb.variables=c(vars, more_cov))
        print(paste('Wrote', length(c(vars, more_cov)), 'features to', fname))
}

remove_poor_samples <- function(my_ukb_data) {
	my_ukb_key <- var_map()
	is_poor <- my_ukb_data[,my_ukb_key$col.name[my_ukb_key$field.showcase == 22010]]
	print(paste('Removing', sum(!is.na(is_poor)), 'samples out of', length(is_poor),
	'due to Field 22010: Recommended genomic analysis exclusions'))
	droplevels(my_ukb_data[is.na(is_poor),])
}
remove_all_related_samples <- function(my_ukb_data) {
	# Expecting to remove 19,381 samples
        my_ukb_key <- var_map()
	not_related <- apply(my_ukb_data[,my_ukb_key$col.name[my_ukb_key$field.showcase ==22011]],
						1, function(x) all(is.na(x)))
	print(paste('Removing', sum(!not_related), 'samples out of', length(not_related), 
			    'related samples. Notice that we remove here all samples from any pairs.',
				'A beter way would be to remove a single sample.'))
	droplevels(my_ukb_data[not_related,])
}
remove_non_caucasian <- function(my_ukb_data) {
     #Expecting to have 409,692 samples 
      my_ukb_key <- var_map()
      is_caucasian <- !is.na(my_ukb_data[,my_ukb_key$col.name[my_ukb_key$field.showcase ==22006]])
      print(paste('Removing', sum(!is_caucasian), 'samples out of', length(is_caucasian),
	  			'due to Field 22006: Caucasian'))
      droplevels(my_ukb_data[is_caucasian,])
}
#' Remove samples with mismattching sex
remove_sex_samples <- function(my_ukb_data) {
        my_ukb_key <- var_map()
	not_defined <- is.na(my_ukb_data[,my_ukb_key$col.name[my_ukb_key$field.showcase ==22001]])
	print(paste('Removed', sum(not_defined),'samples out of', length(not_defined),' due to mismatch sex'))
	droplevels(my_ukb_data[!not_defined,])
}


# create parser object
parser <- ArgumentParser(description='Phenotype codes to extract')
# Phenotype codes
parser$add_argument('integers', metavar='N', type="integer", nargs='+', help='Phenotype codes (data fields)')
parser$add_argument('-o', '--out', default='./ukbb_phenotypes.csv', help='Output filename')
parser$add_argument('-p', '--remove_poor', action='store_true', default=FALSE,
					help='Remove samples with poor genotyping')
parser$add_argument('-r', '--remove_related', action='store_true', default=FALSE,
					help='Remove all related samples')
parser$add_argument('-c', '--remove_non_caucasian', action='store_true', default=FALSE,
					help='Remove non-caucasian samples')
parser$add_argument('-s', '--remove_sex_mismatch', action='store_true', default=FALSE,
					help='Remove samples with mismatching sex')

args <- parser$parse_args()

my_ukb_data <- get_data()
print(paste('Retrieved data (', ncol(my_ukb_data), 'columns) for', nrow(my_ukb_data), 'samples.'))
if (args$remove_poor) {
    my_ukb_data <- remove_poor_samples(my_ukb_data)
}
if (args$remove_sex_mismatch) {
    my_ukb_data <- remove_sex_samples(my_ukb_data)
}
if (args$remove_related) {
    my_ukb_data <- remove_all_related_samples(my_ukb_data)
}
if (args$remove_non_caucasian) {
    my_ukb_data <- remove_non_caucasian(my_ukb_data)
}
if (args$remove_poor | args$remove_sex | args$remove_related | args$remove_non_caucasian) {
    print(paste('After filtration, left with', nrow(my_ukb_data), 'samples.'))
}

write_covar_and_pheno(my_ukb_data, args$integers, args$out)
