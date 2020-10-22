#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o Rout/out_multi			#-- output directory (fill in)
#$ -e Rout/err_multi                     #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r n                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=18G                  #-- submits on nodes with enough free memory (required)
#$ -l h_rt=24:20:00                #-- runtime limit (see above; this requests 24 hours)
#$ -t 1:10

Rscript run.R
