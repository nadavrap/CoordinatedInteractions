#!python3

"""
First need to summarize plink results:
for f in `cat PHENOTYPES`; do for i in `seq 0 9`; do cp ${f}/plink_glm_imp_chr1_fold_not${i}_covariates_standarized_sex_MAF0.001_normalized_9308c5520845776e5d54fef2900800f7.${f}.glm.linear ${f}/norm_plink_glm_imp_fold_${i}.tsv; for c in `seq 2 22`; do cat ${f}/plink_glm_imp_chr${c}_fold_not${i}_covariates_standarized_sex_MAF0.001_normalized_9308c5520845776e5d54fef2900800f7.${f}.glm.linear >> ${f}norm_plink_glm_imp_fold_${i}.tsv; done; done; done

For binary traits (phenotype not normalized):
for f in T2D; do for i in `seq 0 9`; do cp ${f}/plink_glm_imp_chr1_fold_not${i}_covariates_standarized_sex_MAF0.001_*.${f}.glm.logistic ${f}/norm_plink_glm_imp_fold_${i}.tsv; for c in `seq 2 22`; do cat ${f}/plink_glm_imp_chr${c}_fold_not${i}_covariates_standarized_sex_MAF0.001_*.${f}.glm.logistic | tail -n+2 >> ${f}/norm_plink_glm_imp_fold_${i}.tsv; done; done; done
for f in disease_ASTHMA_DIAGNOSED disease_ALLERGY_ECZEMA_DIAGNOSED disease_HI_CHOL_SELF_REP; do for i in `seq 0 9`; do rm -f ${f}/norm_plink_glm_imp_fold_${i}.tsv; cp ${f}/plink_glm_clumped_chr1_fold_not${i}_covariates_standarized_sex_MAF0.001_*.${f}.glm.logistic ${f}/norm_plink_glm_imp_fold_${i}.tsv; for c in `seq 2 22`; do cat ${f}/plink_glm_clumped_chr${c}_fold_not${i}_covariates_standarized_sex_MAF0.001_*.${f}.glm.logistic | tail -n+2 >> ${f}/norm_plink_glm_imp_fold_${i}.tsv; done; done; done

Plink score with p-value 0.001 and chromosome 22 (8500 variants) and bfile (./IMP/ukb_imp_chr22_MAF0.001) took ~10 minutes.
But using the subset clumped pgen file, it took 20 seconds.
Extract clumped variants:
for i in `seq 1 22`; do awk '{print $3}' clump_chr${i}_ex_dups.clumped > chr${i}_clumped.txt; done
Consider generating much smaller clumped pgen files:
~/bin/plink2_linux_x86_64/plink2 --bfile IMP/ukb_imp_chr22_MAF0.001 --extract ./CLUMP/chr22_clumped.txt --out CLUMP/clumped_chr22_pgen --make-pgen

for i in `seq 2 9` `seq 15 17` 19 20; do
awk '{print $3}' CLUMP/clump_chr${i}_ex_dups.clumped > CLUMP/chr${i}_clumped.txt
qsub -N pgen${i} -cwd -l mem_free=1G -pe shared 10 -j y -o log/pgen_${i}_`d`.log -b y -l h_rt=50:00:00 ~/bin/plink2_linux_x86_64/plink2 --bfile IMP/ukb_imp_chr${i}_MAF0.001 --extract CLUMP/chr${i}_clumped.txt --out CLUMP/clumped_chr${i}_pgen --make-pgen --threads 10; done

Run explicit chromosome:
for p in 0.01 0.001 0.0001; do for d in blood_EOSINOPHIL_COUNT; do for c in 21 22; do qsub -t 1-10 -N score${d}_${c} -cwd -l h_data=1G -j y -o ${d}/log/score${p}_\$TASK_ID_`d`.log ~/bin/py_wrapper.sh scripts/plink_score.py -f $d -p $p -c $c; done; done; done

# All chromosomes:
for p in 1.0 0.1 .0.05 0.01 0.001 0.0001; do for d in `cat PHENOTYPES`; do qsub -t 1-10 -N score${d}_${p} -cwd -l h_data=1G -j y -o ${d}/log/score${p}_\$TASK_ID_`d`.log ~/bin/py_wrapper.sh "scripts/plink_score.py -f $d -p $p"; done; done

Tissue specific
T=Liver; F=BMI; p=0.0001; i=0; python3 scripts/plink_score.py -f $F -p $p -i $i -e TissueSpecific/VariantsByTissue/${T}_chromatin_0.bim --suffix $T --imputed
for T in Liver Adipose Brain Brain_Hippocampus Pancreas Skeletal_Muscle; for F in BMI Height; for p in 1.0 0.1 0.01 0.05 0.001 0.0001; do qsub -t 1-10 -cwd -N ${F}${T} -j y -o ${d}/log/plink_score_tissue_${T}${p}_`d`.log -pe shared 1-10 -l h_rt=24:00:00 ~/bin/py_wrapper.sh scripts/plink_score.py -f $F -p $p -e TissueSpecific/VariantsByTissue/${T}_chromatin_0.bim --suffix $T

Use external summary statistics:
p=0.01; d=LDL_direct; for c in `seq 1 22`; do qsub -cwd -N ${d}${c}AUX -j y -o ${d}/log/score_aux_${c}_`d`.log  -pe shared 1-100 -l h_rt=23:00:00,h_data=1G ~/bin/py_wrapper.sh scripts/plink_score.py -f $d -p $p -c $c -b AUX_data/LDL/jointGwasMc_LDL.txt -m; sleep 10; done

# Summarize
for p in 1.0 0.1 0.05 0.01 0.001 0.0001; do for d in T2D; do qsub -N ${d}_${p} -cwd -l h_data=1G -j y -o ${d}/log/score${p}_`d`.log ~/bin/py_wrapper.sh "scripts/plink_score.py -f $d -p $p --summarize"; done; done
"""
import argparse
import time
from math import log  # ,exp

from PRS_functions import *


def run_plink_score(bfile, scorefile, outfile, maf=None, keep=None, extract=None, is_pgen=False, precomputed_freqs=''):
    """Given genotyping file, and estimated effect sizes, run plink2 --score which is simply multiplying \beta*G"""
    if file_not_empty(outfile + '.sscore') and file_not_empty(outfile + '.log'):
        with open(outfile + '.log') as f:
            for line in f:
                if 'Error' in line:
                    break
            if line.startswith('End time:'):
                print(f'Score results already presents: {outfile}')
                return outfile
    maf = '' if maf is None else f'--maf {maf}'
    keep = '' if keep is None else f'--keep {keep}'
    extract = '' if extract is None else f'--extract {extract}'
    target_type = 'p' if is_pgen else 'b'  # Prefix for pgen file or bgen file
    tmpdir = get_tmp_dir()
    outdir, baseoutfile = os.path.split(f'{outfile}')
    cmd = f'{PLINK2} {maf} --out {tmpdir}{baseoutfile} --{target_type}file {bfile} {keep} {extract} ' \
          f'--score {scorefile} header list-variants ignore-dup-ids ' \
          f'cols=maybefid,maybesid,phenos,nallele,dosagesum,scoreavgs,scoresums --threads {NSLOTS} ' \
          f'{precomputed_freqs} --rm-dup force-first '
    # Note, nallele was previosly named nmissallele
    try:
        exec_str(cmd)
    except Exception as e:
        print('Failed to execute, probably not enough variants left.')
        # Check for small chromosome, if no variants left, write a mock empty file
        if file_not_empty(f'{tmpdir}{baseoutfile}.log'):
            with open(f'{tmpdir}{baseoutfile}.log') as f:
                lines = f.readlines()
                if lines[-3].startswith('Error: No valid variants in --score file') and lines[-1].startswith(
                        'End time'):
                    print('Printing empty file' + outfile + '.sscore')
                    with open(outfile + '.sscore', 'w') as o:
                        o.write('\n')

    # Moving files from tmp dir to destination
    for f in glob.glob(tmpdir + '*'):
        destname = f'{outdir}/' + os.path.basename(f)
        if os.path.exists(destname):
            print(f'File {destname} exists and will be replaced')
            os.remove(destname)
        shutil.move(f, outdir)
    # Remove temporary directory
    shutil.rmtree(tmpdir)

    return outfile


def oddratio2beta(x):
    """Given odd ratio, return a matching \beta estimate"""
    if x == 'NA':
        return x
    x = float(x)
    # return exp(x)/(1+exp(x))
    return log(x)


def subset_by_pval(basefile, pvalue_thr):
    """Subset the effect size file to variants with p-value lower than the given threshold."""
    # print(f'In subset_by_pval. Basefile: {basefile}, Pvalue: {pvalue_thr}')
    outfname = os.path.splitext(basefile)[0] + f'_P{pvalue_thr}.tsv'
    waitfname = outfname + '.GENERATING'
    while os.path.isfile(waitfname):
        time.sleep(3)

    if file_not_empty(outfname):
        return outfname
    done_rsids = set()
    open(waitfname, 'w').close()
    with open(basefile) as f:
        with open(outfname, 'w') as w:
            values = f.readline().split()
            header = dict((j.upper(), i) for i, j in enumerate(values))
            ID = header['ID'] if 'ID' in header else header['SNP'] if 'SNP' in header else header[
                'SNPID'] if 'SNPID' in header else header['MARKERNAME'] if 'MARKERNAME' in header else header[
                'RSID'] if 'RSID' in header else header['SNPNAME'] if 'SNPNAME' in header else header['RS']
            A1 = header['A1'] if 'A1' in header else header['EFFECT_ALLELE'] if 'EFFECT_ALLELE' in header else header[
                'ALLELE1'] if 'ALLELE1' in header else header['ALLELE_1']
            BETA = header['BETA'] if 'BETA' in header else header['B'] if 'B' in header else header[
                'OR'] if 'OR' in header else header['OR_FIX']
            P = header['P'] if 'P' in header else header['PVAL'] if 'PVAL' in header else header[
                'P_VALUE'] if 'P_VALUE' in header else header['P-VALUE'] if 'P-VALUE' in header else header[
                'P_DGC'] if 'P_DGC' in header else header['P_FIX']
            or_not_beta = 'OR' in header

            # w.write('\t'.join([values[ID], values[A1], values[BETA]]) + '\n')
            pval = -1
            last = {'ID': 'ID', 'A1': 'A1', 'BETA': 'BETA', 'P': -1}
            for line in f:
                values = line.split()
                if or_not_beta:
                    values[BETA] = oddratio2beta(values[BETA])
                snpid, a1, beta, pval = (values[ID], values[A1], values[BETA], values[P])
                pval = 1.0 if pval == 'NA' else float(pval)
                if snpid == last['ID']:
                    if pval > last['P']:
                        continue
                else:
                    if last['P'] <= pvalue_thr:
                        lastbeta = str(last['BETA'])
                        if lastbeta != 'NA' and last['ID'] not in done_rsids:
                            w.write('\t'.join([last['ID'], last['A1'].upper(), lastbeta]) + '\n')
                            done_rsids.add(last['ID'])
                last = {'ID': values[ID], 'A1': values[A1], 'BETA': values[BETA], 'P': pval}

    os.remove(waitfname)
    return outfname


def main():
    parser = argparse.ArgumentParser(
        description=f'''Given a score file and a genomic file compute the score. Instead of providing chromosome and 
        fold, it cat be extract from environment variable SGE_TASK_ID.''')
    parser.add_argument('-f', '--phenotype', help='UKBB Phenotype to use')
    parser.add_argument('-p', '--pvalues', help='Maximal P-value(s)', type=float, nargs='+')
    parser.add_argument('-c', '--chromosome', help='Chromosome', choices=list(range(1, 23)), type=int)
    parser.add_argument('-i', '--fold_idx', help='Fold index (0-9)', choices=list(range(10)), type=int)
    parser.add_argument('-b', '--basefile', help='Base file with effect sizes')
    parser.add_argument('-s', '--summarize', help='Summarize results over chromosomes.', action='store_true')
    parser.add_argument('-e', '--extract', help='Subset to the variants in the given list.')
    parser.add_argument('-m', '--imputed', action='store_true', default=False,
                        help='Use imputed data (default False), o/w uses clumped data.')
    parser.add_argument('--suffix', help='Add suffix to output filenames.')
    parser.add_argument('--external', action='store_true', default=False,
                        help='Use external data set, and therefore no need to cross validation.')
    parser.add_argument('--AllPCs', action='store_true', default=False, help='Use gwas with 40PCs.')
    parser.add_argument('-t', '--tissuespecific', action='store_true', default=False,
                        help='Perform a tissue specific scoring. Conjugated with -b.')
    parser.add_argument('--remove_score_files', action='store_true', default=False, help='Remove score files (norm_score_imp_chr*_fold?_P*)')

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    phenotype = args.phenotype
    chromosome = args.chromosome
    if chromosome is None:
        chromosomes = range(1, 23)
    else:
        chromosomes = (chromosome,)
    fold_idx = args.fold_idx
    suffix = ('_' + args.suffix) if args.suffix else ''
    suffix = suffix + ('_40PCs' if args.AllPCs else '')
    for pvalue in args.pvalues:
        if args.tissuespecific:
            subdir = f'TissueSpecific/{pvalue}/'
        else:
            subdir = f'{pvalue}/'
        if not os.path.exists(f'{WD}{phenotype}/{subdir}'):
            os.makedirs(f'{WD}{phenotype}/{subdir}')
        impORpgen = 'imp' if args.imputed else 'pgen'

        outfilestr = '{WD}{phenotype}/{subdir}norm_score_{impORpgen}_chr{chromosome}_fold_{fold_idx}_P{pvalue}{suffix}'

        if args.summarize:
            summarize_scores(phenotype, pvalue, outfilestr, suffix, subdir, impORpgen, no_folds=args.external, to_remove_score_files=args.remove_score_files)
            continue

        if args.fold_idx is None:
            fold_idx = SGE_TASK_ID

        basefile = args.basefile
        keep = f'{WD}{phenotype}/fold_{fold_idx}.fam'
        if fold_idx is None:
            keep = f'{WD}{phenotype}/{phenotype}.pheno'
        if basefile is None:
            if args.AllPCs:
                basefile = f'{WD}{phenotype}/norm_plink_glm_imp_fold_{fold_idx}_40pcs.tsv'
            else:
                basefile = f'{WD}{phenotype}/norm_plink_glm_imp_fold_{fold_idx}.tsv'

        sumoutfile = f'{WD}{phenotype}/norm_score_P{pvalue}{suffix}.csv'
        if args.external:
            sumoutfile = f'{WD}{phenotype}/external_norm_score_{impORpgen}_P{pvalue}{suffix}.csv'
        if file_not_empty(sumoutfile) and len(open(sumoutfile).readlines()) > 140000:
            print(f'Summary file already exists: {sumoutfile}')
            continue

        for chromosome in chromosomes:
            if args.imputed:
                # bfile
                #bpfile = f'{IMP_DIR}ukb_imp_chr{chromosome}_MAF0.001'
                bpfile = f'{IMP_DIR}ukb_imp_chr{chromosome}_pgen'
                is_pgen = True
            else:
                # pgen
                bpfile = f'{CLUMP_DIR}clumped_chr{chromosome}_pgen'
                is_pgen = True
            outfile = eval('f"' + outfilestr + '"')
            # extract = f'{CLUMP_DIR}chr{chromosome}_clumped.txt'
            if not args.extract:
                extract = bpfile + '.pvar'
            else:
                extract = args.extract

            precomputed_freqs = ''
            # print(f'{FREQ_DIR}{phenotype}_imp_chr{chromosome}_freq_rm_dup.afreq')
            #if file_not_empty(f'{FREQ_DIR}{phenotype}_imp_chr{chromosome}_freq_rm_dup.afreq'):
            precomputed_freqs_file = f'{FREQ_DIR}ukb_imp_chr{chromosome}_pgen_rmdup.afreq' # f'{FREQ_DIR}ukb_imp_pgen_rmdup.afreq'
            if file_not_empty(precomputed_freqs_file):
                precomputed_freqs = f'--read-freq {precomputed_freqs_file}'

            scorefile = subset_by_pval(basefile, pvalue)
            # run_plink_score(bfile, scorefile, outfile, maf=None, keep=keep, extract=extract)
            try:
                run_plink_score(bpfile, scorefile, outfile, maf=None, keep=keep, extract=extract, is_pgen=is_pgen,
                                precomputed_freqs=precomputed_freqs)
            except:
                pass


def get_phenotype(phenotype):
    """Phenotype can be a name of a phenotype, or an ID of UKBB phenotype.
    Might also add implementation of reading phenotype from a file.
    Returns a name of the phenotype file.
    """
    phenofile = phenotype + '/' + phenotype + '.pheno'
    with open(phenofile) as f:
        _ = f.readline()  # Read header
        return {(line.split()[0], line.split()[1]): line.split()[2] for line in f}


def remove_score_files(phenotype, pvalue, outfilestr, suffix, subdir, impORpgen, no_folds, folds):
    """Remove output files from plink --score. They are not needed anymore once the results were summarized."""
    print('Removing files')
    for fold_idx in folds:
        for chromosome in range(1, 23):
            scorefile = eval('f"' + outfilestr + '"') + '.*'
            for f in glob.glob(scorefile, recursive=False):
                os.remove(f)


def summarize_scores(phenotype, pvalue, outfilestr, suffix, subdir, impORpgen, no_folds=False, to_remove_score_files=False):
    """Given you already have results (score) for each fold and for each chromosome, summarize it to a single file."""
    outfile = f'{WD}{phenotype}/norm_score_P{pvalue}{suffix}.csv'
    nfiles = 22 * 10  # 22 chromosomes times 10 folds
    files_counter = 0
    if no_folds:
        outfile = f'{WD}{phenotype}/external_norm_score_{impORpgen}_P{pvalue}{suffix}.csv'
        nfiles = 22  # No folds
    folds = ('None',) if no_folds else range(10)
    if file_not_empty(outfile):
        with open(outfile) as f:
            if len(f.readlines()) > 2:
                print(outfile + ' already presents.')
                if to_remove_score_files:
                    remove_score_files(phenotype, pvalue, outfilestr, suffix, subdir, impORpgen, no_folds, folds)
                return ()
    phenotypes = get_phenotype(phenotype)
    with open(outfile, 'w') as o:
        o.write('FID,IID,Pheno,' + ','.join(map(str, range(1, 23))) + '\n')
        for fold_idx in folds:
            prss = {}
            for chromosome in range(1, 23):
                # scorefile = f'{WD}{phenotype}/score_pgen_chr{chromosome}_fold_{fold_idx}_P{pvalue}.sscore'
                scorefile = eval('f"' + outfilestr + '"') + '.sscore'
                with open(scorefile) as f:
                    # read header
                    # FID    IID     NMISS_ALLELE_CT NAMED_ALLELE_DOSAGE_SUM SCORE1_AVG
                    _ = f.readline()
                    if len(_) < 2:  # If empty file
                        for k in prss:
                            prss[k].append('0')
                    for line in f:
                        fid, iid, nmiss, nad, prs, prs_sum = line.split()
                        prss.setdefault((fid, iid), list()).append(prs_sum)
                    files_counter = files_counter + 1
            for k, v in prss.items():
                v = v + ['0'] * (22 - len(v))  # Complete missing by zeros
                o.write(f'{k[0]},{k[1]},{phenotypes[k]},' + ','.join(v) + '\n')
    if files_counter == nfiles and to_remove_score_files:
        remove_score_files(phenotype, pvalue, outfilestr, suffix, subdir, impORpgen, no_folds, folds)


if __name__ == "__main__":
    # execute only if run as a script
    main()
