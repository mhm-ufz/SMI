#! /bin/bash
#-----------------------------------------------------
#$ -S /bin/bash                       # Working Shell
#$ -N SMI_EU                        # Set job name
#$ -o /work/thober/StdOut/smi_$TASK_ID.txt  # Log file
#$ -e /work/thober/StdOut/smi_$TASK_ID.txt
#$ -j y                               # Write error in log file
# submission for all model ensembles, months, and years
# #$ -t 1-2591
# submission for all model members, months, and years
# #$ -t 1-32723
# submission for Grand Ensemble, or Ensemble for one model
# #$ -t 1-323
# submission for Statistical model
# #$ -t 1-3239
# submission for ESP 28 years X 12 months
#$ -t 1-325
#$ -l h_rt=0:10:00                    # Request resource: hard run time hours:minutes:seconds
#$ -l h_vmem=1G                       # Request resource: memory requirements/per slot
#$ -l centos6=false                   # run on eve idiv  
#$ -cwd                               # Change into directory where you wrote qsub
#$ -binding linear:1
#----------------------------------------------------

time ./launch_smi.sh -t $SGE_TASK_ID
