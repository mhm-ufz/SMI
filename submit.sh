#! /bin/bash
#-----------------------------------------------------
#$ -S /bin/bash            # Working Shell
#$ -N smi_mHMv5_pan_de     # Set job name
#$ -o /data/stohyd/mHM_project/germany/pan_german/default_eve2/smi/$JOB_NAME.$JOB_ID # Log file
#$ -j y                    # Write error in log file
#$ -l h_rt=36:00:00         # Request resource: hard run time hours:minutes:seconds
#$ -l h_vmem=10G           # Request resource: memory requirement
#$ -l centos6=true         # OS
#$ -cwd                    # Change into directory where you wrote qsub
#$ -m ea                   # mail notification e - in case of job ended a - in case of abortion 
#$ -M matthias.zink@ufz.de # mail address
#----------------------------------------------------

export OMP_NUM_THREADS=1

cd /data/stohyd/mHM_project/germany/pan_german/default_eve2_2/smi_2/

time ./smi
