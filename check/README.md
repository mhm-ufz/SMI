# smi check cases

With the provided python script `smi_checks_run.py` you can run all checks in the `case_**` folders with one or more given smi executable(s).

It will be checked, if smi finishes normally and if the output matches the reference values.

A summary for all cases is given at the end.

## Usage

    python smi_checks_run.py [-h] [-e EXE [EXE ...]] [-v] [-l LOG_PATH] [-t OPENMP_THREADS] [-s [SKIP [SKIP ...]]] [-o [ONLY [ONLY ...]]]

Run the smi check cases with a given smi executable.

## Optional Arguments

    -h, --help            show this help message and exit
    -e EXE [EXE ...], --exe EXE [EXE ...]
                          Paths to smi exe[s]. (default: ['../smi'])
    -v, --verbose         Show the smi output. (default: no)
    -l LOG_PATH, --log_path LOG_PATH
                          Directory for smi-logs. (default: the resp. case dir)
    -t OPENMP_THREADS, --threads OPENMP_THREADS
                          Number of threads for openMP. No mpi allowed! (default: 0)
    -s [SKIP [SKIP ...]], --skip [SKIP [SKIP ...]]
                          skip cases (case_01 case_03 ..) (default: [])
    -o [ONLY [ONLY ...]], --only [ONLY [ONLY ...]]
                          only run cases (case_01 case_03 ..) (default: all)

## Examples

Run smi from parent directory in verbosity mode with mpi on 4 processes:

        python smi_checks_run.py -e ../smi -v -m 4

Silently run smi (given with an absolute path) with openmp on 4 threads:

        python smi_checks_run.py -e /abspath/smi_openmp -t 4

Silently run smi from parent directory:

        python smi_checks_run.py

Run with multiple smi exes:

        python smi_checks_run.py -e ../smi1 ../smi2

## Cleanup

To remove the created output of smi run:

        python smi_checks_clean.py

## Author

    Sebastian Mueller

## Contributors

    Stephan Thober, Robert Schweppe

Written Jan. 2022.
