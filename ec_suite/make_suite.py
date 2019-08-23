# export ECF_HOST=datascience1
# export ECF_PORT=47345
# ecflow_server --port=47345
# ecflow_client --load=ecf_launch.def
# ecflow_client --begin=forward


import ecflow as ec
import os
from shutil import copyfile as cp

sm_path = '/work/thober/pgb/'
sm_files = [ff for ff in os.listdir(sm_path) if '.nc' in ff[-3:]] 

# start suite definition
suite_name = 'forward'
family_name = 'run_smi'
task_name = 'forward'

smi_exe = '/home/thober/lib/smi/smi'
ecf_home = '/work/thober/pgb/smi'
ecf_files = '/home/thober/lib/smi/ec_suite/' # not used in the sense that ECflow foresees

ecs = ec.Suite(suite_name,
               ec.Edit(ECF_HOME=ecf_home,
                       ECF_FILES=ecf_home,
                       ECF_INCLUDE=os.path.join(ecf_files, 'include'),
                       ECF_TRIES=1,
                       SMI_EXE=smi_exe))

# add family for each file
for ff in sm_files:
    dir_name = ff[:-3].replace('-', '_')
    dir_task = os.path.join(ecf_home, suite_name, dir_name)
    os.makedirs(dir_task, exist_ok=True)
    # copy task template
    file_name = task_name + '.ecf'
    template_file = os.path.join(ecf_files,suite_name,family_name,file_name)
    cp(template_file, os.path.join(dir_task, file_name))
    #
    ecs.add(ec.Family(dir_name,
            ec.Task(task_name),
            ec.Edit(SOILMOIST_FILE=os.path.join(sm_path, ff))))

defs = ec.Defs(ecs)
print(defs)

print("Checking ecflow setup")
# check job creation
print("    Checking job creation: .ecf -> .job0")
print(defs.check_job_creation(throw_on_error=True, verbose=True))

# check triggers
print("    Checking trigger expressions")
assert len(defs.check()) == 0,defs.check()

print("Saving definition to file 'ecf_launch.def'")
defs.save_as_defs("ecf_launch.def")

print('Done!')
