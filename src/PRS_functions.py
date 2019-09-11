#!python3
"""Set of functions needed in multiple scripts."""
import glob
import shlex
import shutil
import subprocess
import tempfile

from parameters import *


def mv_files_from_tmp(tmpdir, outdir):
    """Move files from tmpdir to outdir"""
    for f in glob.glob(tmpdir + '*'):
        destname = f'{outdir}/' + os.path.basename(f)
        if os.path.exists(destname):
            os.remove(destname)
        shutil.move(f, outdir)
    # Remove temporary directory
    shutil.rmtree(tmpdir)


def get_tmp_dir():
    """Create and return a temporary directory under TMP_DIR."""
    return tempfile.mkdtemp(dir=TMP_DIR) + '/'


def extractfilename(filename, keep_suffix=False):
    """Given a filename, return file name w/o path and w/o extension.
    :rtype: string
    """
    if keep_suffix:
        return os.path.basename(filename)
    return os.path.splitext(os.path.basename(filename))[0]


def file_not_empty(file_name):
    return os.path.isfile(file_name) and os.path.getsize(file_name) > 0


def exec_str(command, shell=False, stdout=None, verbose=False, env=None):
    """Execute the given string. If failed, exit and print status.
    """
    if verbose:
        print('executing: ' + command)
    if env is None:
        stats = subprocess.call(shlex.split(command), shell=shell, stdout=stdout)
    else:
        stats = subprocess.call(shlex.split(command), shell=shell, stdout=stdout, env=env)
    if stats != 0:
        sys.stdout.flush()
        sys.stderr.flush()
        raise (Exception(f'''Execution failed.\nExit status: {stats}'''))


# from typing import Set, Any, Tuple, Union
# def join_chromatin_ranges(tissue_name, prefixes, window_size=DEFAULT_WINDOW_SIZE):
#     """Given prefixes of .bed (range) files, join them and add window_size"""
#     ranges_file_name = f'{tissue_name}_{window_size}.ranges'
#     if file_not_empty(ranges_file_name):
#         return ranges_file_name
#     ranges: Set[Tuple[Any, int, Union[int, Any]]] = set()
#     for prefix in prefixes:
#         with open(prefix + 'bed') as f:
#             for line in f:
#                 c, b, e = map(int, line.replace('chr', '').split())
#                 b = min(0, b - window_size)
#                 e = e + window_size
#                 ranges.add((c, b, e))
#     with open(ranges_file_name, 'w') as o:
#         for val in ranges:
#             o.write('\t'.join(map(str, ranges)) + '\n')
#     return ranges_file_name


# def tissue_specifc_variants(tissue, source, join=True, window_size=50000):
#     """Use Hilary's data to extract tissue-specific variants.
#     - tissue Name of tissue. Can also be a prefix
#     - source should be GE for gene-expression or chromatin for chromatin markers.
#     - join - in case of multiple tissues (or same tissue with multiple markers) join them together. If join=False and
#         there is more than one match, raises an exception.
#     - window_size - window size around each variant (default 50Kb from each side).
#     """
#     data_dir = WD + 'TissueSpecific/'
#     outfile = f'{data_dir}/{tissue}_{source}_{window_size}'
#     bimoutfile = outfile + '.bim'
#     if file_not_empty(bimoutfile):
#         return bimoutfile
#     if source.lower() == 'ge':
#         mappingfile = data_dir + 'Multi_tissue_gene_expr.ldcts'
#         join_method = ''
#     elif source.lower() == 'chromatin':
#         mappingfile = data_dir + 'Multi_tissue_chromatin.ldcts'
#         join_method = join_chromatin_ranges
#
#     tissues = {}
#     with open(mappingfile) as f:
#         for line in f:
#             tissue_name, tissue_files = line.split()
#             if tissue_name.lower().startswith(tissue.lower()):
#                 tissues[tissue_name] = tissue_files.split(',')[0]
#     if len(tissues) == 0:
#         raise (Exception(f'Tissue {tissue} was not found'))
#     if len(tissues) > 1 and not join:
#         raise (Exception(f'Tissue {tissue} is ambigous. Found next possible matches: ' + (', '.join(tissues.keys()))))
#     fnames = list(tissues.values())[0]
#
#     ranges_file = join_chromatin_ranges(tissue, fnames, window_size)
#
#     cmd = f'''{PLINK2} --pfile IMP/ukb_imp_chr{chromosome}_pgen --extract range {tissue_set_range_file} \
#                 --make-just-bim --out {outfile} --threads {NSLOTS}'''
#     exec_str(cmd)
#     return (bimoutfile)


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def write_effect_sizes(infile, outfile, with_header):
    """Given an effect sizes file, write another file where A1 is the effect allel, and A2 the reference"""
    flag = 'w' if with_header else 'a'
    with open(outfile, flag) as o:
        if with_header:
            o.write('\t'.join(['CHROM', 'ID', 'POS', 'A1', 'A2', 'BETA', 'P']) + '\n')
        with open(infile) as f:
            _ = f.readline()
            for line in f:
                CHROM, POS, ID, REF, ALT, A1, TEST, OBS_CT, BETA, SE, T_STAT, P = line.split()
                A2 = REF if A1 != REF else ALT
                o.write('\t'.join([CHROM, ID, POS, A1, A2, BETA, P]) + '\n')
