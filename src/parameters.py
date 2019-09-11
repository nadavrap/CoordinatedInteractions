#!python3
"""
# This file is imported to every script, so make sure not to add variables that may conflict with others.
# This files include parameters and variables that are needed for the scripts.
# In order to load it to the global scope of a python script, add:
> from parameters import *
# To get a single parameter, run:
> python3 <param name>
so the value will be printed to the STDOUT.
"""
import os
import sys

#####################
# Paths to software #
#####################
'''Might be: '/zaitlen/netapp/group/UKBIOBANK/' or '/u/project/sriram/nadavrap/UKBB/ or '/coldstorage/nadav/UKBB/'''
WD = './'  # Working directory, make sure you have write permission
DATA_DIR = WD
TMP_DIR = os.environ["TMPDIR"] if "TMPDIR" in os.environ else '/tmp/'
'''On some systems, python2.7 is still the default, and python3 is not available. Therefore we need to use a wrapper
A wrapper shell script may look like one of the next two:
    #!/bin/bash
    source /u/local/Modules/default/init/modules.sh
    module load python/3.7.2 zlib &> /dev/null
    python3 "$@"
OR
    #!/bin/bash
    scl enable rh-python36 "python $@"
'''
HOME_DIR = os.path.expanduser('~') + '/'
PYTHON3 = HOME_DIR + 'bin/py_wrapper.sh'  # Replace with python3 if that works on your system.
RSCRIPT = 'Rscript'
PLINK1 = 'plink'
PLINK2 = 'plink2'
BOLT_PATH = 'BOTL-LMM/BOLT-LMM_v2.3.2/'
BOLT = BOLT_PATH + 'bolt'
'''If Rscript is not in the path, and need to be load, uses a R_wrapper.sh that looks like:
    #!/bin/bash
    module load CBI r
    Rscript $@
And replace value of next variable to 'R_wrapper.sh'
'''
PRSICE = RSCRIPT + f'''{HOME_DIR}/PRSice/PRSice.R --prsice {HOME_DIR}/PRSice/PRSice_linux --dir {HOME_DIR}/PRSice'''

###################
# UKBB data files #
###################
UKB_FNAME = 'ukb27646'
GENOTYPE_DIR = f'{WD}GEN/'
# PRSICE_TARGET='%sukb_cal_chr#_v2,%sukb30397_cal_chr1_v2_s488339.fam' % (GENOTYPE_DIR, GENOTYPE_DIR)
PRSICE_TARGET = '%sukbb' % GENOTYPE_DIR
IMP_DIR = DATA_DIR + 'IMP/'
CLUMP_DIR = DATA_DIR + 'CLUMP/'

##########################
# SGE cluster parameters #
##########################
PARALLEL_ENVIRONMENT = 'shared'  # Can be smp, parallel, shared
QUEUE = '-q long.q'  # Can be empty for default
GIG = 2  # GB per slots
MIN_CORES = 1
MULTI_CORE = f'-pe {PARALLEL_ENVIRONMENT} {MIN_CORES}-'

###############
# Parameters #
##############
# Minor allele frequency (fraction) in plink notation: 0.001 is 0.1%
MAF = 0.001
# Window size for variants around genes. This distance in bp will be taken before gene start and after gene's end
DEFAULT_WINDOW_SIZE = 50000

if __name__ == '__main__':
    '''This is used when need a variable not in python (e.g., R or shell script.)
    To get a variable in R:
    > PATH <- system('python parameters.py WD' ,intern=TRUE)
    To get a variable in shell:
    > WD=`python parameters.py WD`
    '''
    param_name = sys.argv[1]
    param_value = eval(param_name)
    print(param_value)
