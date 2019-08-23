#!/bin/ksh
set -e          # stop the shell on first error
set -u          # fail when using an undefined variable
set -o pipefail # fail if last(rightmost) command exits with a non-zero status
set -x          # echo script lines as they are executed

# Defines the variables that are needed for any communication with ECF
export ECF_PORT=%ECF_PORT%    # The server port number
export ECF_HOST=%ECF_HOST%    # The name of ecf host that issued this task
export ECF_NAME=%ECF_NAME%    # The name of this current task
export ECF_PASS=%ECF_PASS%    # A unique password
export ECF_TRYNO=%ECF_TRYNO%  # Current try number of the task
export ECF_RID=$$
export ECF_TIMEOUT=300 # Only wait 5 minutes, if the server cannot be contacted (note default is 24 hours) before failing
#export ECF_DEBUG_CLIENT=1

# SANITY Check, typically only valid for new platforms. make sure hostname is resolvable to an IP address
# host %ECF_HOST%

# Tell ecFlow we have started
export PATH=${PATH}:/global/apps/thober/ecflow/4.14.0/bin
%ECF_CLIENT_EXE_PATH:ecflow_client% --init=$$

# Defined a error handler
ERROR() {
   echo 'ERROR() called'
   set +e                      # Clear -e flag, so we don't fail
   wait                        # wait for background process to stop
   %ECF_CLIENT_EXE_PATH:ecflow_client% --abort=trap   # Notify ecFlow that something went wrong, using 'trap' as the reason
   trap 0                      # Remove the trap
   exit 0                      # End the script
}

# Trap any calls to exit and errors caught by the -e flag
trap ERROR 0

# Trap any signal that may cause the script to fail
trap '{ echo "Killed by a signal"; ERROR ; }' 1 2 3 4 5 6 7 8 10 12 13 15
