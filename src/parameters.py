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
HOME_DIR = '/u/project/zaitlenlab/nadavrap/' #os.path.expanduser('~') + '/'
'''Might be: '/zaitlen/netapp/group/UKBIOBANK/' or '/u/project/sriram/nadavrap/UKBB/ or '/coldstorage/nadav/UKBB/'''
#WD = HOME_DIR + 'UKBB/'  # Working directory, make sure you have write permission
WD = '/u/project/sriram/nadavrap/UKBB/'
DATA_DIR = WD
SCRIPTS_PATH = f'{WD}CoordinatedInteractions/src/'
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
PYTHON3 = HOME_DIR + 'bin/py_wrapper.sh'  # Replace with python3 if that works on your system.
RSCRIPT = HOME_DIR + 'bin/R_wrapper.sh'
PLINK1 = 'plink'
PLINK2 = HOME_DIR + 'bin/plink2'
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
UKB_FNAME = 'ukb21608' #'ukb27646'
GENOTYPE_DIR = f'{WD}GEN/'
# PRSICE_TARGET='%sukb_cal_chr#_v2,%sukb30397_cal_chr1_v2_s488339.fam' % (GENOTYPE_DIR, GENOTYPE_DIR)
PRSICE_TARGET = '%sukbb' % GENOTYPE_DIR
IMP_DIR = DATA_DIR + 'IMP/'
CLUMP_DIR = DATA_DIR + 'CLUMP/'
FREQ_DIR = DATA_DIR + 'FREQ/'

##########################
# SGE cluster parameters #
##########################
PARALLEL_ENVIRONMENT = 'shared'  # Can be smp, parallel, shared (Hoffman2)
QUEUE = '' # '-q long.q'  # Can be empty for default
GIG = 2  # GB per slots
MIN_CORES = 4
MULTI_CORE = f'-pe {PARALLEL_ENVIRONMENT} {MIN_CORES}-'
NSLOTS = int(os.environ["NSLOTS"]) if "NSLOTS" in os.environ else 1
SGE_TASK_ID=os.environ["SGE_TASK_ID"] if "SGE_TASK_ID" in os.environ else None
if SGE_TASK_ID == 'undefined':
    SGE_TASK_ID = None
elif SGE_TASK_ID is not None:
    SGE_TASK_ID = int(SGE_TASK_ID) - 1


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
