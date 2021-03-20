#!R
#' Given a phenotype and a list of patients, first remove outliers, then quantile normalize it.
#'
#'

get_pheno <- function(phenotype) {
    read.table(paste0(phenotype, '/', phenotype, '.pheno'), header=TRUE)
}

get_samples <- function(famfile) {
    read.table(famfile, header=TRUE)
}

remove_outliers <- function(d) {
    outlier_values <- boxplot.stats(d[,3])$ou
    print(paste('Removing', sum(d[,3] %in% outlier_values) , 'outliers out of',
                nrow(d), '(', round(sum(d[,3] %in% outlier_values)/nrow(d)*100),'%)'))
    d[!d[,3] %in% outlier_values, ]
}

quantile_normalization <- function(d) {
    d[,3] <- preprocessCore::normalize.quantiles.use.target(as.matrix(d[,3]), target=rnorm(nrow(d)))[,1]
    d
}

remove_and_normalize <- function(phenotype, fold_ix, phenos, not) {
    outfile <- paste0(phenotype, '/', phenotype, '_fold_', ifelse(not, 'not', ''), fold_ix, '.pheno')
    famfile <- paste0(phenotype, '/', 'fold_', ifelse(not, 'not', ''), fold_ix, '.fam')
    samples <- get_samples(famfile)
    phenos <- phenos[phenos$IID %in% samples$IID, ]
    write.table(phenos, file=outfile, row.names = FALSE, quote=FALSE)
    phenos
}

suppressPackageStartupMessages(library("argparse"))
# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--phenotype", help="Phenotype of interest.")
args <- parser$parse_args()
phenotype <-  args$phenotype
if(is.null(phenotype)) {
    print('Phenotype was not provided')
    quit()
} else {
    print(paste('Normalizing', phenotype))
}
phenos <- get_pheno(phenotype)
# Remove outliers
phenos <- remove_outliers(phenos)
# Normalize
phenos <- quantile_normalization(phenos)

res <- lapply(0:9, function(i) remove_and_normalize(phenotype, i, phenos, not=TRUE))
res <- lapply(0:9, function(i) remove_and_normalize(phenotype, i, phenos, not=FALSE))
res <- do.call(rbind, res)

write.table(res, file = paste0(phenotype, '/', phenotype, '_normalized.pheno'), row.names = FALSE, quote=FALSE)

