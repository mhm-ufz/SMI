#!/bin/bash
#
#set -e
#set -x
#
# bash script for executing the Program Eval (ME)
#
# author: Stephan Thober
#
# created: 5.10.2012
#
# ------------------------------------------------------------------------------
#
pid=$$

# get a task from option
task=-9999
while getopts "t:" Option ; do
    case ${Option} in
        t) task="${OPTARG}";;
        *) printf "Error: unimplemented option.\n\n";  exit 1;;
    esac
done

if [[ ${task} -eq -9999 ]]; then
    echo '***ERROR: negative task'
    exit 1
fi

# compile the program
# make > /dev/null && echo 'succesfully compiled!'

# JobThres of Jobs allowed to run
JobThres='1000'

# parameters
year_start=1983
Nyears=27
# year_start=1982
# Nyears=28
month_start=1
Nmonths=12
# !!! NMME models !!!
models="CMC1-CanCM3 COLA-RSMAS-CCSM3 GFDL-CM2p1 IRI-ECHAM4p5-AnomalyCoupled IRI-ECHAM4p5-DirectCoupled NASA-GMAO-062012 NCEP-CFSv1 NCEP-CFSv2"
Nmodels=8
cum_realizations="10 16 26 38 50 62 77 101"
Ndisag=25

# !!! ONE MODEL ONLY !!!
models="NCEP-CFSv2"
Nmodels=1
cum_realizations="24"
Ndisag=0
year_start=2009
Nyears=1
month_start=3
Nmonths=1

# !!! Statistical model !!!
# models="AR_forecast"
# Nmodels=1
# cum_realizations="10"

# !!! ESP !!!
# models=''
# Nmodels=28
# cum_realizations="0"

# extract model 
(( model_ID = ( ${task} / (${Nyears}*${Nmonths}) ) ))
(( task     =   ${task} - ${model_ID} * (${Nyears}*${Nmonths}) ))
(( month_ID =   ${task} / (${Nyears}) ))
(( year_ID  =   ${task} - ${month_ID} * ${Nyears} ))

# set model
(( model_ID += 1 ))
# cc=${model_ID}

# >>>>>>> for full model realizations >>>>>>>>>
# extract model and model realization from model id
cc=0
for c_id in ${cum_realizations}; do
    if [[ ! ${model_ID} -gt ${c_id} ]]; then
        if [[ cc -eq 0 ]]; then
            model_rea=${model_ID}
        else
            cc_id=$(echo ${cum_realizations} | cut -d ' ' -f ${cc})
            (( model_rea = ${model_ID} - ${cc_id} ))
            echo ${cc_id} ${cc}
        fi
        break
    fi
    (( cc++ ))
done
(( cc++ ))
# set member
mem=$(printf "%2.2i" ${model_rea})
# >>>>>>> for full model realizations >>>>>>>
# set model
model=$(echo ${models} | cut -d ' ' -f ${cc})
# model=${model_ID}
# set month
(( month_ID = ${month_ID} + ${month_start} ))
mm=$(printf "%2.2i" ${month_ID})
# set year
(( year_ID  = ${year_ID}  + ${year_start} ))
yy=$(printf "%2.2i" ${year_ID})

# work path
wpath='/work/thober/forecast_T+P/'
# wpath='/work/thober/esp/'

# set all variables that are going to be changed in mhm.nml
OutPath=${wpath}${yy}'_'${mm}'/subensembles/'         # submit for subensembles
InPath='/data/stohyd/mHM_project/europe/entire_domain/ST_RK_06_2014/forecast/pre_tavg/'
InPath=${InPath}${yy}'_'${mm}'/subensembles/'
#
# implementation of backward search algorithm
ens_config='11111111'
Nmodels=8 # Number of models equals the length of ens_config
#

# create outpath if not exists
if [[ ! -d ${wpath} ]]; then
    echo '***ERROR: working path does not exist:'
    echo ${wpath}
    exit 1
fi

# Loop over configurations (see average_subensembles.sh for details)
for (( cc = 0; cc < ${Nmodels}; cc++ )); do
    # if model has not been used, loop
    if [[ ${ens_config:${cc}:1} -eq 0 ]]; then continue; fi
    # copy configuration and set 0 at position cc
    config=${ens_config:0:${cc}}0${ens_config:((${cc}+1)):${Nmodels}}
    echo 'processing config: '${config}
    # set Paths for this configuration
    SubOutPath=${OutPath}${config}'/'
    SubInPath=${InPath}${config}'/'
    EvalFile=${SubInPath}'sm_anomaly.nc'
    OutFile=${SubInPath}mSMI.nc
    echo ${OutFile}
    
    # check whether outfile exists
    # if [[ -f ${OutFile} ]]; then 
    #     echo 'SKIPPING!'
    #     exit 0
    # fi

    # check if script is run on local machine or frontends of eve
    machine=$(uname -a | cut -d ' ' -f 2)

    if [[ ${machine%[0-9]} != 'ces220l' ]]; then
        
        echo 'processing InFile:  '${EvalFile}
        echo 'processing OutFile: '${OutFile}

        # create run dir in OutPath
        if [[ -d ${SubOutPath}run/ ]]; then
	    rm -r ${SubOutPath}run/
        fi

        mkdir -p ${SubOutPath}run/
        echo ${SubOutPath}run/
        # copy all required files in this directory
        cp main.dat ${SubOutPath}run
        ln -s ${PWD}/smi ${SubOutPath}run/smi
        # change directory
        cd ${SubOutPath}run/${Fname%.*}

        # modify main namelist
        sed -i -e "s# outpath        = .*# outpath        = \"${SubInPath}\"#" \
               -e "s# SM_eval_file   = .*# SM_eval_file   = \"${EvalFile}\"#" \
            main.dat
        
        # commit the job
        #time ./smi
        
        # go back
        cd - > /dev/null

    else
        
        echo '***ERROR: not implemented!'
        exit 1

    fi
done
#
#echo 'Done!'
