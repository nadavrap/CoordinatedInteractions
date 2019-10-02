#!python3
# For a given phenotype and a given p-value
# For each patient note:
# ID, Pheno, PRS all PRS odd, PRS even, PRS tissue_i all, PRS tissue_i odd, PRS tissue_i even
import argparse
from typing import List

from PRS_functions import *

WD = './'
TISSUES = ('', '_Liver', '_Adipose', '_Brain', '_Brain_Hippocampus', '_Pancreas', '_Skeletal_Muscle', '_Franke.11',
           '_Franke.120', '_Franke.14', '_Franke.46', '_Franke.69', '_Franke.84', '_Franke.89')


def prs_scores(scores_file):
    """
    :param scores_file: Filename with PRS scores
    :return: Dictionary: patient_id => tuple(phenotype, prs_all, prs_even)
    """
    prss = {}
    with open(scores_file) as f:
        _ = f.readline().split(',')
        for line in f:
            line_values: List[str] = line.split(',')
            fid, iid, pheno = line_values[:3]
            prs = list(map(float, line_values[3:]))
            odd = sum(prs[::2])
            even = sum(prs[1::2])
            all_prs = odd + even
            prss[fid] = (pheno, str(all_prs), str(odd), str(even))
    return prss


def summarize_prss(phenotype, pvalue):
    """
    Summarize PRS by tissue and write it to file.
    :param phenotype: Name of phenotype.
    :param pvalue: Pvalue to analyze.
    :return: None
    """
    outfile = f'{WD}/{phenotype}/tissue_specific_PRS_P{pvalue}.csv'
    if file_not_empty(outfile):
        return outfile
    prss = []
    for tissue in TISSUES:
        scores_file = f'{WD}/{phenotype}/norm_score_P{pvalue}{tissue}.csv'
        prss.append(prs_scores(scores_file))
    with open(outfile, 'w') as w:
        colnames = [['All' + t, 'Odd' + t, 'Even' + t] for t in TISSUES]
        colnames = [item for sublist in colnames for item in sublist]
        colnames = ['FID', 'Pheno'] + colnames
        w.write(','.join(colnames) + '\n')
        for pid in prss[0].keys():
            vals = [pid, prss[0][pid][0]]
            scores = [d[pid][1:] for d in prss]
            scores = [item for sublist in scores for item in sublist]
            w.write(','.join(vals + scores) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''A tool for summarize tissue specific PRSs''')
    parser.add_argument('-f', '--phenotype', help='UKBB Phenotype to use')
    parser.add_argument('-p', '--pvalue', help='Pvalue', default='0.01')
    args = parser.parse_args()
    summarize_prss(args.phenotype, args.pvalue)
