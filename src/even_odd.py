#!python3
"""
Split samples for cross validation, and estimate effect sizes.
Can be run with a specified fold number and explicit chromosome, or can submit all the needed jobs
(#folds * #chromosomes) as separated jobs.

Example of submissions:
# PHENO=BMI; qsub -pe shared 10- -N ${PHENO}_even_odd python3 scripts/even_odd.py --phenotype ${PHENO} \
--betagenerator Plink --predictor PRSice --evenodd --nfolds 10 --verbose
python3 CoordinatedInteractions/src/even_odd.py --nfolds 10 --imputed --normalized --AllPCs --phenotype Height --betagenerator Plink --predictor None --samples_to_include CSV/covariates_norel3.cov --verbose
"""
import argparse
import collections
import datetime
import hashlib
import socket
import time

import numpy
import pandas

from PRS_functions import *

# plink --bfile GEN/ukbb --freq --out snpQC/ukbb_freq
# awk -v MAF="${MAF}" '{if ($5 < MAF) {print $1}}' ukbb_freq.frq > maf_lt_${MAF}.txt
MAF_SNPS = WD + 'snpQC/maf_lt_%s.txt' % MAF
TABLES_DIR = BOLT_PATH + 'tables/'
EXTRACT_SCRIPT = SCRIPTS_PATH + 'ukb_extract_phenotype.R'
MERGED_GENOTYPE = GENOTYPE_DIR + 'ukbb'
GENOTYPE_BED_FILES = GENOTYPE_DIR + 'ukb_cal_chr{1:22}_v2.bed'
DEFAULT_FACTOR_COVS = 'uk_biobank_assessment_centre_f54_0_0,genotype_measurement_batch_f22000_0_0,sex_f31_0_0'
SAMPLES_TO_EXCLUDE = WD + 'sampleQC/samples_to_exclude.fam'
# DEFAULT_COVAR_FILE=WD+'CSV/covariates2.cov'
# Standarize covariates using plink:
# plink2 --bfile ukbb --covar covariates2.cov --covar-variance-standardize --write-covar --out covariates_standarized
# DEFAULT_COVAR_FILE=WD+'CSV/covariates_standarized.cov'
#DEFAULT_COVAR_FILE = WD + 'CSV/covariates_standarized_sex.cov'  # Where Sex replaced from Male/Femal/NONE to 1/0/NA
DEFAULT_COVAR_FILE = WD + 'CSV/covariates_stand_norel3.cov'  # Where Sex replaced from Male/Femal/NONE to 1/0/NA
BOLT_COVAR_FILE = WD + 'CSV/covariates.cov'
PRSICE_NO_COVAR = None
# Long string with all dummy covariates
DEFAULT_COVAR_NAMES = open(WD + 'CSV/default_covar_names.txt').readline().strip()
ALL_PCS_COVAR_NAMES = open(WD + 'CSV/covarnames_40PCs.txt').readline().strip()


def verboseprint(x):
    """function will be defined to print if verbose is one."""
    return None


def hash_list(*argv):
    s = ''.join(map(str, argv)).encode()
    hashed = hashlib.md5(s).hexdigest()
    # verboseprint(f'Hashing: {s} To: {hashed}')
    return hashed


def get_even_odd_chrs_files():
    """Need to implemnt generation of file. Currently just return filename"""
    '''plink2 --bfile ukbb --chr 1,3,5,7,9,11,13,15,17,19,21 --make-just-bim --out ukbb_odd_chrs
    plink2 --bfile ukbb --chr 2,4,6,8,10,12,14,16,18,20,22 --make-just-bim --out ukbb_even_chrs
    for pref in even odd; do cut -f2 ukbb_${pref}_chrs.bim > ukbb_${pref}_chrs.txt; done
    '''
    return GENOTYPE_DIR + 'ukbb_even_chrs.txt', GENOTYPE_DIR + 'ukbb_odd_chrs.txt'


def get_phenotypefile(phenotype, verbose=False):
    """Phenotype can be a name of a phenotype, or an ID of UKBB phenotype.
    Might also add implementation of reading phenotype from a file.
    Returns a name of the phenotype file.
    """
    phenofile = WD + phenotype + '/' + phenotype + '.pheno'
    if os.path.isfile(phenofile):
        verboseprint('Phenotype file already present for %s' % phenotype)
    else:
        cmd_str = f'''{RSCRIPT} {EXTRACT_SCRIPT} --pheno_name {phenotype} --out {phenofile}'''
        exec_str(cmd_str, verbose=verbose)
        verboseprint('End extracting phenotype')
    return phenofile


def run_bolt(fold_idx, phenofile, covarnames, phenotype, covarfile=BOLT_COVAR_FILE, extract=None, append='', **kwargs):
    """Run BOLT-LMM to estimate effect sizes
    BOLT runs 20 hours for Height with clumping
    TODO Need to implement additional covariates like BMI (if different from phenotype)
    For case/controls, BOLT output betas. To get Log OR:
    logOR = beta_bolt/(mu*(1-mu)) where mu being the prevalence: mu=case/(case+control)
    """
    keepsamples = WD + phenotype + '/fold_not%s.fam' % fold_idx
    removesamples = WD + phenotype + '/fold_%s.fam' % fold_idx
    verboseprint('Use as covar file: %s' % BOLT_COVAR_FILE)
    outfname = WD + phenotype + '/bolt_lmm_' + extractfilename(keepsamples) + '_' + extractfilename(covarfile) + '_' + (
        extractfilename(extract) if extract is not None else '') + '_' + hash_list(covarnames) + append
    statsfile = outfname + '.stats.gz'
    blupfile = outfname + '.betas'
    bmicovar = '--qCovarCol=body_mass_index_bmi_f21001_0_0' if 'body_mass_index_bmi_f21001_0_0' in covarnames else ''
    if not os.path.isfile(statsfile) or os.path.getsize(statsfile) == 0 or not os.path.isfile(
            blupfile) or os.path.getsize(blupfile) == 0:
        ldtable = TABLES_DIR + 'LDSCORE.1000G_EUR.tab.gz'
        maptable = TABLES_DIR + 'genetic_map_hg19_withX.txt.gz'
        cmd_str = f'''{BOLT_PATH} \
        --bed={GENOTYPE_DIR}ukb_cal_chr{1:22}_v2.bed \
        --bim={GENOTYPE_DIR}ukb_snp_chr{1:22}_v2.bim \
        --fam={GENOTYPE_DIR}ukbb_for_BOLT.fam \
        --remove={SAMPLES_TO_EXCLUDE} \
        --remove={removesamples} \
        --exclude={MAF_SNPS} \
        --phenoFile={phenofile} \
        --phenoCol={phenotype} \
        --covarFile={covarfile} \
        --covarCol=uk_biobank_assessment_centre_f54_0_0 \
        --covarCol=genotype_measurement_batch_f22000_0_0 \
        {bmicovar} \
        --covarCol=sex_f31_0_0 \
        --covarMaxLevels=110 \
        --qCovarCol=age_when_attended_assessment_centre_f21003_0_0 \
        --qCovarCol=genetic_principal_components_f22009_0_{1:10} \
        --LDscoresFile={ldtable} \
        --geneticMapFile={maptable} \
        --lmmForceNonInf \
        --numThreads={NSLOTS} \
        --statsFile={statsfile} \
        --modelSnps={extract} \
        --verboseStats \
        --predBetasFile={blupfile}'''
        exec_str(cmd_str)
        verboseprint('End BOLT-LMM')
    return blupfile


def run_plink_glm(fold_idx, phenofile, covarfile, covarnames, extract=None, phenotype=None, is_binary_trait=False,
                  append='', chrom=None, **kwargs):
    """Run Plink to estimate effect sizes.
    Time examples, Height with 10 cores took 20 hours for odd chromosomes
    Height with 28 cores took 20 hours
    Height with 28 cores, 9 hours odd chromosomes
    keepsamples - Filename with list of samples to keep.
    Do we need also samples to remove?
    extract - Filename with list of variants to keep.
    phenofile - Filename with phenotype. Assuem a single phenotype is the file. O/w need to change implementation to
    accept phenotype name.
    # phenotype - Name of phenotype - column name in phenofile.
    covarfile - Name of covariate file
    covarnames - List of covariate names to uses.
    """
    verboseprint('Lets train using Plink glm')
    verboseprint('Variants to extract:' + str(extract))
    keepsamples = WD + phenotype + '/fold_not%s.fam' % fold_idx
    outdir = WD + phenotype + '/'
    tmpdir = get_tmp_dir()
    fname = 'plink_glm_' + extractfilename(keepsamples) + '_' + extractfilename(covarfile) + '_' + (
        extractfilename(extract) if extract is not None else '') + '_' + hash_list(covarnames) + append
    outfname = outdir + fname
    tmpoutfname = tmpdir + fname
    resfilename = outfname + ('.%s.glm.' % phenotype) + ('logistic' if is_binary_trait else 'linear')
    logfilename = outfname + '.log'
    # Check if there is a log file and ended as expected
    if os.path.isfile(logfilename) and os.path.getsize(logfilename) > 0:
        # Read log file and make sure no Error or Failed, and ended in a good way
        with open(logfilename) as logfh:
            for line in logfh:
                if 'Error' in line or 'Failed' in line:
                    break
            if line.startswith('End time:'):
                verboseprint('Plink results already presents')
                return resfilename

    memory = NSLOTS * GIG * 750
    extract = '' if extract is None else ('--extract ' + extract)

    freq_file = f'{FREQ_DIR}ukb_imp_chr{chrom}_pgen_rmdup.afreq'
    freq_file_opt = f'--read-freq {freq_file} --extract {freq_file}' if file_not_empty(freq_file) else ''
    verboseprint('Frequencies files were %sfound' % ('' if freq_file else 'not '))

    glm_modifier = 'cc-residualize no-firth' if is_binary_trait else '' # firth-residualize

    verboseprint('Allocating %s MB (%s slots * %sGB per slots * 950 MB/GB)' % (memory, NSLOTS, GIG))
    args = f'''--bfile {MERGED_GENOTYPE} --keep {keepsamples} {extract} {freq_file_opt}\
    --threads {NSLOTS} --pheno {phenofile} --pheno-name {phenotype} --glm {glm_modifier} hide-covar --covar {covarfile} \
    --covar-name {covarnames} --out {tmpoutfname} --memory {memory}'''

    my_env = os.environ.copy()
    my_env["OMP_NUM_THREADS"] = "1"  # str(NSLOTS//2)
    exec_str(f'''{PLINK2} {args}''', env=my_env)
    verboseprint('End running plink glm')

    mv_files_from_tmp(tmpdir, outdir)

    return resfilename


def run_plink_glm_imputed(fold_idx, phenofile, covarfile, covarnames, extract=None, phenotype=None,
                          is_binary_trait=False, append='', chrom=None, maf=MAF, normalized=False, clumped=False):
    """Run Plink to estimate effect sizes.
    Time examples, Height with 10 cores took 20 hours for odd chromosomes
    Height with 28 cores took 20 hours
    Height with 28 cores, 9 hours odd chromosomes
    keepsamples - Filename with list of samples to keep.
    Do we need also samples to remove?
    extract - Filename with list of variants to keep.
    phenofile - Filename with phenotype. Assueme a single phenotype is the file. O/w need to change implementation to
    accept phenotype name.
    # phenotype - Name of phenotype - columname in phenofile.
    covarfile - Name of covariate file
    covarnames - List of covariate names to uses.
    normalized - Normalized phenotypes.
    clumped - Uses clumped data (default False), usefule for binary traits with super low run time.
    """
    verboseprint('Lets train using Plink glm')
    verboseprint('Variants to extract: %s' % ('All chromosome' if extract is None else extract))
    keepsamples = WD + phenotype + '/fold_not%s.fam' % fold_idx
    prefix = 'clumped' if clumped else 'imp'
    outdir = WD + phenotype + '/'
    tmpdir = get_tmp_dir()
    fname = f'/plink_glm_{prefix}_chr' + chrom + '_' + extractfilename(keepsamples) + '_' + \
            extractfilename(covarfile) + \
            (('_' + extractfilename(extract, keep_suffix=True)) if extract is not None else '') + '_' + 'MAF' + \
            str(maf) + ('_normalized' if normalized else '') + '_' + hash_list(covarnames) + append
    outfname = outdir + fname
    tmpoutfname = tmpdir + fname
    resfilename = outfname + ('.%s.glm.' % phenotype) + ('logistic' if is_binary_trait else 'linear')
    logfilename = outfname + '.log'
    # Check if there is a log file and ended as expected
    if file_not_empty(logfilename):
        with open(logfilename) as logfh:
            for line in logfh:
                if 'Error' in line:
                    break
            if line.startswith('End time:'):
                verboseprint('Plink results already presents')
                return resfilename
    if normalized:
        phenofile = f'{phenotype}/{phenotype}_fold_not{fold_idx}.pheno'

    memory = NSLOTS * GIG * 950
    to_extract = '' if extract is None else f'--extract {extract}'
    pfile = f'{CLUMP_DIR}clumped_chr{chrom}_pgen' if clumped else f'{IMP_DIR}ukb_imp_chr{chrom}_pgen'
    freq_file = f'{FREQ_DIR}ukb_imp_chr{chrom}_pgen_rmdup.afreq'
    freq_file_opt = f'--read-freq {freq_file} --extract {freq_file}' if file_not_empty(freq_file) else ''
    verboseprint('Frequencies files %sfound' % ('' if freq_file else 'not '))

    glm_modifier = 'cc-residualize no-firth' if is_binary_trait else '' # firth-residualize no-firth firth-residualize 
    
    verboseprint('Allocating %s MB (%s slots * %sGB per slots * 950 MB/GB)' % (memory, NSLOTS, GIG))
    exec_str(f'''{PLINK2} --pfile {pfile} --keep {keepsamples} {to_extract} --maf {MAF} {freq_file_opt} \
    --threads {NSLOTS} --pheno {phenofile} --pheno-name {phenotype} --glm {glm_modifier} hide-covar --covar {covarfile}\
    --covar-name {covarnames} --out {tmpoutfname} --memory {memory}''')

    verboseprint('End running plink glm imputed')

    # Move files from tmpdir to outdir
    mv_files_from_tmp(tmpdir, outdir)

    return resfilename


def add_dummy_pval(fname):
    """In case there is no P-values (like results of BOLT-LMM blups, then add dummy pvalue of 1."""
    outfname = '/scratch/' + os.path.basename(fname)
    verboseprint('Betas file %s have no P-values. Generating dummy p-values into %s' % (fname, outfname))
    with open(fname) as f:
        with open(outfname, 'w') as o:
            line = f.readline()
            o.write(line.strip() + '\tP\n')
            for line in f:
                o.write(line.strip() + '\t1\n')
    return outfname


def run_prsice(keepsamples, basefile, phenofile, phenotype, extract, covarfile=None, covarnames=None, to_clump=True,
               is_binary_trait=False, target=PRSICE_TARGET, factor_covs=DEFAULT_FACTOR_COVS):
    """Run PRSice to estimate PRS
    covarfile -
    covarnames -
    basefile - Filename with the effect sizes
    factor_covs -
    extract - Filename with list of variants to extract (to keep)
    Run time ~ 1 hour
    :param extract: Variants to extract
    :param phenotype: Name of phenotypes
    :param keepsamples: File with samples to keep
    :param to_clump: To clump the data
    :param covarnames: Name of covariates
    :param covarfile: Covariates file
    :param phenofile: File with phenotypes
    :param factor_covs: List of covariates that should be treated as factors
    :param target: Name of target (to compute PRS)
    :param is_binary_trait: Is it a binary trait?
    :rtype: basestring
    :type basefile: basestring
    """
    verboseprint('In run_prsice')
    basebasename = extractfilename(basefile)
    # basebasename = basebasename[:basebasename.rfind('_')]
    outfname = WD + phenotype + '/prsice_' + basebasename + '_' + (
        extractfilename(covarfile) if covarfile is not None else 'no_covar') + (
                   '_clump_' if to_clump else '_') + f'MAF{MAF}_' + hash_list(keepsamples, covarnames, extract,
                                                                              to_clump, is_binary_trait, target,
                                                                              factor_covs)
    resfilename = outfname + '.prsice'
    if not file_not_empty(resfilename) and not file_not_empty(outfname + '.best') and not file_not_empty(
            outfname + '.summary'):
        verboseprint('Will generate results to: ' + outfname)
        verboseprint('Will clump? : %s' % to_clump)
        with open(basefile) as f:
            values = f.readline().split()
            allel1_col = 'A1' if 'A1' in values else 'ALLELE1' if 'ALLELE1' in values else exit(
                'What is the column for affected allele?')
            snp_col = 'ID' if 'ID' in values else 'SNP' if 'SNP' in values else exit('What is the column for SNP id?')
            beta_col = 'BETA' if 'BETA' in values else 'OR' if 'OR' in values else exit(
                'What is the column with the betas?')
            if 'P' not in values:
                basefile = add_dummy_pval(basefile)
        pval_col = 'P'
        cov_str = ''
        to_extract = '' if (extract is None or extract == 'None') else f'--extract {extract}'
        if covarfile is not None and covarnames is not None:
            cov_str = ''' --cov-file %(covarfiles)s \
            --cov-factor %(factor_covs)s \
            --cov-col %(covarnames)s '''

        command_str = f'''%(prsice)s \
        --print-snp \
        --snp %(snp_col)s \
        --base %(basefile)s \
        --A1 %(allel1_col)s \
        --pvalue %(pval_col)s \
        --target %(target)s \
        --thread %(threads)s \
        --stat %(beta_col)s \
        %(beta)s \
        --binary-target %(is_binary)s \
        --keep %(keepsamples)s \
        --lower 0.0 \
        --pheno-file %(phenofile)s \
        --pheno-col %(phenotype)s \
        --score sum %(cov_str)s %(to_clump)s {to_extract} \
        --maf {MAF} \
        --ld-maf {MAF} \
        --out %(outfname)s ''' % {'prsice': PRSICE, 'cov_str': cov_str, 'basefile': basefile, 'allel1_col': allel1_col,
                                  'pval_col': pval_col,
                                  'target': target, 'beta_col': beta_col, 'is_binary': 'T' if is_binary_trait else 'F',
                                  'snp_col': snp_col,
                                  'keepsamples': keepsamples, 'threads': NSLOTS, 'phenofile': phenofile,
                                  'phenotype': phenotype,
                                  'to_clump': '' if to_clump else '--no-clump', 'extract': extract,
                                  'outfname': outfname, 'beta': '--beta' if beta_col == 'BETA' else ''}
        exec_str(command_str)
        verboseprint('End running PRSice')
    else:
        verboseprint('PRSice results %s already present' % resfilename)
    return outfname


def split_to_nfolds(phenotype, phenotype_fname, nfolds, samples_to_exclude=None, samples_to_include=None, normalized=False):
    """Need to remove samples:
    - with no phenotype
    - Related
    - Non caucasian
    - Poor genotypes
    """
    # Get phenotype filename
    if os.path.isfile(WD + phenotype + '/fold_0.fam'):
        verboseprint('Fam files already splited')
        return
    # Read phenotype file and shuffle
    df = pandas.read_csv(phenotype_fname, delim_whitespace=True).sample(frac=1)
    # Remove samples with no phenotype
    no_pheno = df[pandas.isna(df[phenotype])][['FID', 'IID']]
    df = df[pandas.notna(df[phenotype])][['FID', 'IID']]
    # Remove samples to exclude
    if samples_to_exclude:
        to_exclude = pandas.read_csv(samples_to_exclude, delim_whitespace=True, header=None, names=['FID', 'IID'])
        to_exclude = pandas.concat([to_exclude, no_pheno]).drop_duplicates().reset_index(drop=True)
        df_all = df.merge(to_exclude.drop_duplicates(), on=['FID', 'IID'], how='left', indicator=True)
        df = df_all[df_all['_merge'] == 'left_only'][['FID', 'IID']]
        to_exclude.to_csv(WD + phenotype + '/to_exclude.fam', sep='\t', index=False)
    if samples_to_include:
        to_include = pandas.read_csv(samples_to_include, delim_whitespace=True, usecols=[0])
        df = df[df.FID.isin(to_include.iloc[:,0].tolist())]
    # Split to N folds
    folds = numpy.array_split(df, nfolds)
    verboseprint('Writing cross val fam files')
    for i, fold in enumerate(folds):
        filename = WD + phenotype + '/fold_%s.fam' % i
        fold.to_csv(filename, sep='\t', index=False)
        filename = WD + phenotype + '/fold_not%s.fam' % i
        pandas.concat(folds[:i] + folds[i + 1:]).to_csv(filename, sep='\t', index=False)

    if normalized:
        # Normalze phenotype, and write phenotype by fold
        verboseprint('Normalized phenontype')
        cmd_str = f'''{RSCRIPT} {SCRIPTS_PATH}pheno_normalization.R --phenotype {phenotype}'''
        exec_str(cmd_str)


def is_binary_pheno(phenofile):
    with open(phenofile) as f:
        f.readline()
        c = collections.Counter([line.split()[-1] for line in f])
    if len(list(c)) == 3:
        if sorted(list(c)) == ['-9', '1', '2'] or sorted(list(c)) == ['0', '1', '2']:
            verboseprint('Phenotype looks binary. Found the next values:')
            verboseprint(c)
            return True
    elif len(list(c)) == 2 and sorted(list(c)) == ['1', '2']:
        verboseprint('Phenotype looks binary. Found the next values:')
        verboseprint(c)
        return True
    verboseprint('Phenotype seems to be continuous')
    return False


def merge_plink(pheno, fold_idx, chrs=range(1, 23), chrgroup='all'):
    dirpath = f'{WD}{pheno}/'
    outfile = f'{dirpath}plink_glm_imp_fold_{fold_idx}_{chrgroup}.tsv'
    if file_not_empty(outfile):
        return outfile
    f0 = glob.glob(
        f"{dirpath}plink_glm_imp_chr{chrs[0]}_fold_not{fold_idx}_covariates_standarized_sex_MAF{MAF}_*.{pheno}.glm"
        f".linear")[0]
    write_effect_sizes(f0, outfile, with_header=True)
    for c in chrs[1:]:
        f = glob.glob(f'{dirpath}plink_glm_imp_chr{c}_fold_not{fold_idx}_covariates_standarized_sex_MAF{MAF}_*'
                      '.{pheno}.glm.linear')[0]
        write_effect_sizes(f, outfile, with_header=False)
    return outfile


def summarize_prs(pattern, outfile, snps, phenotype, phenofile):
    """Given a pattern of PRSice results file (e.g., 'prsice_plink_glm_fold_*best'), summarize and write to outfile
    Summary includes mean, median and stdev of each sample.
    :param pattern: pattern of files
    :param outfile: Filename to output
    :param snps: SNPs set
    :param phenotype: Phenotype
    :param phenofile: Name of file with phenotype
    """
    files = glob.glob(pattern)
    prefixes = tuple(extractfilename(s) for s in snps)
    verboseprint('Found %s files' % len(files))

    with open(phenofile) as f:
        _ = f.readline()  # Read header
        phenotypes = {(line.split()[0], line.split()[1]): line.split()[2] for line in f}
    prss = {p: {} for p in prefixes}
    for filename in files:
        prefix = [p for p in prefixes if p in filename][0]
        verboseprint('Reading file: %s, matching SNPs set: %s' % (filename, prefix))
        with open(filename) as f:
            f.readline()  # Read header
            for line in f:
                vals = line.split()
                # prss[prefix].setdefault(tuple(vals[:2]), []).append(vals[-1])
                prss[prefix][tuple(vals[:2])] = vals[-1]
    with open(outfile, 'w') as o:
        o.write('\t'.join(('FID', 'IID', phenotype) + tuple('PRS_' + p for p in prefixes)) + '\n')
        for k in phenotypes.keys():
            o.write('\t'.join(k) + '\t' + phenotypes.get(k, '') + '\t' + '\t'.join(
                tuple(prss[p].get(k, '') for p in prefixes)) + '\n')


def summarize_plink_glm(phenotype, fold_idx):
    outfile = f'{phenotype}/norm_plink_glm_imp_fold_${fold_idx}.tsv'
    outfile = f'{phenotype}/norm_plink_glm_imp_fold_${fold_idx}_40pcs.tsv'
    # for i in `seq 0 9`; do rm -f ${f}/norm_plink_glm_imp_fold_${i}.tsv; cp ${f}/plink_glm_imp_chr1_fold_not${i}_covariates_standarized_sex_MAF0.001_normalized_*.${f}.glm.linear ${f}/norm_plink_glm_imp_fold_${i}.tsv; for c in `seq 2 22`; do grep -v '^#CHROM' ${f}/plink_glm_imp_chr${c}_fold_not${i}_covariates_standarized_sex_MAF0.001_normalized_*.${f}.glm.linear >> ${f}/norm_plink_glm_imp_fold_${i}.tsv; done; done
    # f=LDL_direct; for i in `seq 0 9`; do rm -f ${f}/norm_plink_glm_imp_fold_${i}_40pcs.tsv; cp ${f}/plink_glm_imp_chr1_fold_not${i}_covariates_standarized_sex_MAF0.001_normalized_a4e674b672f03793978b931770ca9426*.${f}.glm.linear ${f}/norm_plink_glm_imp_fold_${i}_40pcs.tsv; for c in `seq 2 22`; do grep -v '^#CHROM' ${f}/plink_glm_imp_chr${c}_fold_not${i}_covariates_standarized_sex_MAF0.001_normalized_a4e674b672f03793978b931770ca9426*.${f}.glm.linear >> ${f}/norm_plink_glm_imp_fold_${i}_40pcs.tsv; done; done
    # f=T2D; for i in `seq 0 9`; do rm -f ${f}/norm_plink_glm_imp_fold_${i}.tsv; cp ${f}/plink_glm_clumped_chr1_fold_not${i}_covariates_standarized_sex_MAF0.001_*.${f}.glm.logistic ${f}/norm_plink_glm_imp_fold_${i}.tsv; for c in `seq 2 22`; do grep -v '^#CHROM' ${f}/plink_glm_clumped_chr${c}_fold_not${i}_covariates_standarized_sex_MAF0.001_*.${f}.glm.logistic >> ${f}/norm_plink_glm_imp_fold_${i}.tsv; done; done


def main():
    parser = argparse.ArgumentParser(description='''A tool for Even-Odd PRS prediction. 
    Trying to be as modular as possible, that is why many arguments are given.''')
    parser.add_argument('-f', '--phenotype', help='UKBB Phenotype to use')
    parser.add_argument('-t', '--phenotypefile',
                        help='UKBB Phenotype file name. If not specified, will generate one based on the given '
                             'phenotype.',
                        default=None)
    # parser.add_argument('-c', '--covarfile', help='Covariates file', default=None)
    # parser.add_argument('-a', '--covarnames',
    #   help='Covriate names. Notice: Not implemented yet for prsice.', default=None)
    # parser.add_argument('-d', '--addinalcovarnames',
    #   help='Additional covarnames. A common one would be "body_mass_index_bmi_f21001_0_0"', default='')
    parser.add_argument('-b', '--betagenerator', help='Tool for effect size estimator', choices=('BOLT-LMM', 'Plink'),
                        default='Plink')
    parser.add_argument('-p', '--predictor',
                        help='Tool for PRS prediction. Notice that Plink , and Self are not yet implemented',
                        choices=('PRSice', 'Plink', 'None'), default='PRSice')
    group = parser.add_mutually_exclusive_group()  # List of snps is mutually exclusive with even/odd split
    group.add_argument('-s', '--snps', nargs='+', default=None,
                       help='List of files with list of variants to train seperately.')
    group.add_argument('-e', '--evenodd', action='store_true', default=False,
                       help='Split variants by even/odd chromosomes')
    parser.add_argument('-n', '--nfolds', type=int, default=10, help='Number of folds for cross validation')
    parser.add_argument('-v', '--verbose', action='store_true', default=True, help='Print verbose')
    parser.add_argument('-i', '--fold', type=int, default=None, help='Fold index to run.')
    parser.add_argument('-a', '--append', default='',
                        help='Suffix to append to betagenerator results file. Useful if want to test a slightly '
                             'different run.')
    parser.add_argument('-m', '--imputed', action='store_true', default=False, help='Use imputed data')
    parser.add_argument('-c', '--chromosome', default=None, help='Chromosome index.')
    parser.add_argument('-u', '--impsubset', default=None, type=int,
                        help='Number of files that bim/pvar file was splitted to.')
    parser.add_argument('--normalized', default=False, action='store_true',
                        help='Use normalized and outlier removed phenotypes.')
    parser.add_argument('--clumped', default=False, action='store_true',
                        help='Use pre-clumped data based on MAF (default False).')
    parser.add_argument('--AllPCs', default=False, action='store_true', help='Use 40 PCs instead of 10 as default.')
    parser.add_argument('--samples_to_include', default=None,
                        help='File name where first column are sample ids to keep. Is only used in split to folds.')

    args = parser.parse_args()

    _verboseprint = print if args.verbose else lambda *a, **k: None
    global verboseprint
    verboseprint = _verboseprint

    verboseprint(time.asctime(time.localtime(time.time())))
    verboseprint('Hostname: %s' % socket.gethostname())

    phenotype = args.phenotype
    verboseprint('Got the phenotype:' + phenotype)
    if not os.path.exists(WD + phenotype):
        os.makedirs(WD + phenotype)
        os.makedirs(WD + phenotype + '/log/')
    phenofile = get_phenotypefile(phenotype, verbose=args.verbose) if args.phenotypefile is None else args.phenotypefile
    evenodd = args.evenodd
    snps = args.snps
    if snps is None:
        snps = [None]
    nfolds = args.nfolds
    # covarfile = args.covarfile
    covarnames = ALL_PCS_COVAR_NAMES if args.AllPCs else DEFAULT_COVAR_NAMES
    if args.phenotype.lower() != 'bmi':
        verboseprint('Adding BMI as covariate to Plink')
        covarnames = covarnames + ',body_mass_index_bmi_f21001_0_0'

    to_clump = True
    split_to_nfolds(phenotype, phenofile, nfolds, samples_to_exclude=None, samples_to_include=args.samples_to_include, normalized=args.normalized)
    target = ''  # Will be set downstream
    betagenerator = args.betagenerator
    if betagenerator.lower() == 'plink':
        if args.imputed:
            betagenerator = run_plink_glm_imputed
        else:
            betagenerator = run_plink_glm
        target = PRSICE_TARGET
        verboseprint('Will train using Plink')
    elif betagenerator.lower() in ('bolt', 'bolt-lmm'):
        betagenerator = run_bolt
        to_clump = False
        target = '%sukb_cal_chr#_v2,%sukb30397_cal_chr1_v2_s488339.fam' % (GENOTYPE_DIR, GENOTYPE_DIR)
        verboseprint('Will train using BOLT-LMM')
    else:
        exit('Wrong training method. Should be Plink or BOLT-LMM')
    binary_trait = is_binary_pheno(phenofile)

    i = args.fold if args.fold is not None else SGE_TASK_ID
    if args.impsubset is not None and snps[0] is None:
        i = SGE_TASK_ID // args.impsubset
        subfile = (SGE_TASK_ID % args.impsubset)
        snps = [f'{IMP_DIR}TMP_CHR{args.chromosome}/ukb_imp_chr{args.chromosome}_pgen_{subfile}.pvar']
        verboseprint('Will run plink using a subset of variants for fold {i} and subfile {subfile}')
    if i is not None:
        train_fam_file = WD + phenotype + '/fold_not%s.fam' % i
        test_fam_file = WD + phenotype + '/fold_%s.fam' % i
        # Train
        resfile = betagenerator(fold_idx=i, phenofile=phenofile, covarfile=DEFAULT_COVAR_FILE, covarnames=covarnames,
                                phenotype=phenotype, extract=snps[0], is_binary_trait=binary_trait, append=args.append,
                                chrom=args.chromosome, normalized=args.normalized, clumped=args.clumped)
        # Test
        if args.predictor.lower() == 'prsice':
            run_prsice(keepsamples=test_fam_file, basefile=resfile, phenotype=phenotype, covarfile=PRSICE_NO_COVAR,
                       covarnames=PRSICE_NO_COVAR, phenofile=phenofile, extract=snps[0], to_clump=to_clump,
                       is_binary_trait=binary_trait, target=target, factor_covs=DEFAULT_FACTOR_COVS)
        else:
            verboseprint('Do not run predictor (PRSice)')
    else:
        append_str = ('-a %s' % args.append) if args.append != '' else ''
        all_PCs_opt = '--AllPCs' if args.AllPCs else ''
        normalized = '--normalized' if args.normalized else ''
        if evenodd:
            verboseprint('Split running to even/odd chromosomes')
            snps = get_even_odd_chrs_files()
        elif snps is None:  # Use all variants
            verboseprint('Does not subset variants, use all of them')
            snps = MERGED_GENOTYPE + '.bim'
        for i in range(nfolds):
            if args.imputed:
                for chrom in range(1, 23):
                    # resfile = run_plink_glm_imputed(fold_idx=i, phenofile=phenofile, covarfile=DEFAULT_COVAR_FILE,
                    # covarnames=covarnames, phenotype=phenotype, extract=None, is_binary_trait=binary_trait,
                    # append=args.append, chrom=str(i))
                    logfile = f'{WD}{phenotype}/log/fold_{i}_chr_{chrom}_%s.log' % (
                        datetime.datetime.today().strftime('%Y_%m_%d_%H_%M_%S'))
                    jobname = f'%s{i}_{chrom}' % phenotype.replace('blood_', '')
                    expected_run_time = '%s:00:00' % int(25-chrom) # Range for different chromosomes
                    cmd_str = f'''qsub -S /bin/bash -j y -wd {WD} -o {logfile} {MULTI_CORE} {QUEUE} -l h_rt={expected_run_time} \
                                       -N {jobname} {PYTHON3} "{SCRIPTS_PATH}/even_odd.py \
                                       --phenotype {phenotype} --betagenerator {args.betagenerator} \
                                       --predictor {args.predictor} \
                                       {all_PCs_opt} {normalized}\
                                       --nfolds {args.nfolds} -i {i} %(verbose)s {append_str} -t {phenofile} \
                                       %(imputed)s --chromosome {chrom}"''' % {
                        'verbose': '--verbose' if args.verbose else '', 'imputed': '--imputed' if args.imputed else ''}
                    exec_str(cmd_str, verbose=args.verbose, shell=False) #, env=os.environ.copy())
            else:
                for s in snps:
                    logfile = WD + phenotype + '/log/fold_%s_snps_%s_%s.log' % (
                        i, extractfilename(s), datetime.datetime.today().strftime('%Y_%m_%d_%H_%M_%S'))
                    try:
                        jobname = phenotype + str(i) + extractfilename(s)[5]  # Add o/e from odd/even
                    except:
                        jobname = phenotype + str(i) + os.path.basename(s)  # If no even/odd

                    # Submit a new job to the queue
                    cmd_str = f'''qsub -S /bin/bash {QUEUE} -j y -wd {WD} -o %(logfile)s {MULTI_CORE} -N %(jobname)s \
                    {PYTHON3} "{SCRIPTS_PATH}/even_odd.py --phenotype %(phenotype)s --betagenerator %(betagen)s --predictor \
                    %(predictor)s --snps %(snps)s --nfolds %(nfolds)s -i %(i)s %(verbose)s %(append)s -t %(phenofile)s \
                    {all_PCs_opt} {normalized} %(imputed)s"''' % {
                        'logfile': logfile, 'jobname': jobname, 'phenotype': phenotype, 'betagen': args.betagenerator,
                        'predictor': args.predictor, 'snps': s, 'nfolds': args.nfolds, 'i': i,
                        'verbose': '--verbose' if args.verbose else '', 'GIG': GIG, 'append': append_str,
                        'phenofile': phenofile, 'imputed': '--imputed' if args.imputed else ''}
                    exec_str(cmd_str, verbose=args.verbose, shell=True)
    # Print time at the end of script
    verboseprint(time.asctime(time.localtime(time.time())))


if __name__ == "__main__":
    # execute only if run as a script
    main()
