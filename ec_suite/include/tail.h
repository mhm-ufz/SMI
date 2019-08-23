%ECF_CLIENT_EXE_PATH:ecflow_client% --complete    # Notify ecFlow of a normal end
echo "Done!"
trap 0                 # Remove all traps
exit 0                 # End the shell
