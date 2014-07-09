#! /bin/bash
#-----------------------------------------------------
#$ -S /bin/bash                       # Working Shell
#$ -N SMI_EU                          # Set job name
#$ -o /work/thober/StdOut/smi_$TASK_ID.txt  # Log file
#$ -e /work/thober/StdOut/smi_$TASK_ID.txt
#$ -j y                               # Write error in log file
# # submission for all model ensembles, months, and years
# #$ -t 1-2591
# # submission for all model members, months, and years
# #$ -t 1-32723
# submission for all model members, months, and years
#$ -t 1-323
#$ -l h_rt=2:00:00                    # Request resource: hard run time hours:minutes:seconds
#$ -l h_vmem=6G                       # Request resource: memory requirements/per slot
#$ -l centos6=true                    # run on eve idiv  
#$ -cwd                               # Change into directory where you wrote qsub
#----------------------------------------------------

time ./launch_smi.sh -t $SGE_TASK_ID
