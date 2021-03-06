# -*- coding: utf-8 -*-
"""
Run the SMI check cases with a given SMI executable.

Author
------
    Sebastian Mueller

Contributor
-----------
    Stephan Thober, Robert Schweppe

Created
-------
    Feb 2020

Examples
--------
    Run smi from parent directory in verbosity mode with mpi on 4 processes:
        python run_smi_checks.py -e ../smi -v -m 4

    Silently run smi (given with an absolute path) with openmp on 4 threads:
        python run_smi_checks.py -e /abspath/smi_openmp -t 4

    Silently run smi from parent directory:
        python run_smi_checks.py

    Run with multiple smi exes:
        python run_smi_checks.py -e ../smi1 ../smi2
"""
from hashlib import new
import sys
import os
import shutil
import glob
import argparse

# dependecies
import pexpect
from pexpect.popen_spawn import PopenSpawn
import numpy as np
import xarray as xr
import pandas as pd

# pexpect.spawn not present on windows
SPAWN = PopenSpawn if sys.platform == "win32" else pexpect.spawn

# Constants
RTOL = 1e-03  # allowed relative tolerance in output files
ATOL = 1e-04  # allowed relative tolerance in output files
# patterns for output files
PATTERNS = [
    "*.nc",
    "*.txt",
]
# output folder and reference folder
OUT = "output"
REF = "output_save"
# patterns for comparison of restart files (key: var name in output_save, value: var name in output_b1)
MATCH_VARS = {}
IGNORE_VARS = []
EXE = ["../smi"]
# those cases are currently ignored as there are some problems with init from restart, case11 needs to be reconsidered
SKIP = []
SKIP_CASES_MPI = []

# ARGUMENT PARSER #############################################################


def parse_args():
    """
    Parse the given arguments.

    Returns
    -------
    str
        Path to smi exe.
    bool
        Verbosity.
    int
        Number of mpi processes.
    int
        Number of openmp threads.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    parser.add_argument(
        "-e",
        "--exe",
        action="store",
        nargs="+",
        default=EXE,
        dest="exe",
        help="Paths to smi exe[s]. (default: {})".format(EXE),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        dest="verbose",
        help="Show the smi output. (default: no)",
    )
    parser.add_argument(
        "-l",
        "--log_path",
        action="store",
        default=None,
        dest="log_path",
        help="Directory for smi-logs. (default: the resp. case dir)",
    )
    parser.add_argument(
        "-m",
        "--mpi",
        action="store",
        default=0,
        type=int,
        dest="mpi_nop",
        help="Number of mpi processes. No openMP allowed! (default: 0)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        action="store",
        default=0,
        type=int,
        dest="openmp_threads",
        help="Number of threads for openMP. No mpi allowed! (default: 0)",
    )
    parser.add_argument(
        "-s",
        "--skip",
        action="store",
        nargs="*",
        default=SKIP,
        dest="skip",
        help="skip cases (case_01 case_03 ..) (default: {})".format(SKIP),
    )
    parser.add_argument(
        "-o",
        "--only",
        action="store",
        nargs="*",
        default=None,
        dest="only",
        help="only run cases (case_01 case_03 ..) (default: all)",
    )
    args = parser.parse_args()
    return (
        args.exe,
        args.verbose,
        args.log_path,
        args.mpi_nop,
        args.openmp_threads,
        args.skip,
        args.only,
    )


# HELPER FUNCTIONS ############################################################


def tab(n):
    """Tab for printing."""
    return " " * 2 * n


def sep_text(*texts, sep_n=70, tab_n=0):
    """Separated text for printing."""
    sep_line = os.linesep + tab(tab_n) + "#" * sep_n
    text_line = ""
    for text in texts:
        text_line += os.linesep + tab(tab_n) + "# " + text
    return sep_line + text_line + sep_line + os.linesep


def print_comparison(*result, tab_n=3):
    """Print the result of a file comparison."""
    diff_n, big_diff_n, total, miss, diff_vars = result
    print(tab(tab_n), "{} of {} records differ".format(diff_n, total))
    print(
        tab(tab_n),
        "{} of {} records differ more than {}".format(big_diff_n, total, ATOL),
    )
    if miss:
        print(tab(tab_n), "Missing: ", *miss)
    if diff_vars:
        print(tab(tab_n), "Differing: ", *diff_vars)


def create_out_dir(path):
    """Create the output directory. Erase exiting one."""
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


def get_files(patterns, folder):
    """Get list of files following given patterns in a given folder."""
    file_paths = []
    for pattern in patterns:
        file_paths += glob.glob(os.path.join(folder, pattern))
    base_names = [os.path.basename(path) for path in file_paths]
    return file_paths, base_names


# COMPARE ROUTINES ############################################################


def compare_patterns(patterns, new_dir, ref_dir, tab_n=2):
    """Compare files following given name patterns in two directories."""
    diff_sum = 0
    big_diff_sum = 0
    total_sum = 0
    miss_files = []
    info = {}
    new_files, new_names = get_files(patterns, new_dir)
    ref_files, ref_names = get_files(patterns, ref_dir)
    for ref_file, ref_name in zip(ref_files, ref_names):
        print(tab(tab_n), ref_name)
        if ref_name in new_names:
            diff_n, big_diff_n, total, miss, diff_vars = compare_files(
                new_files[new_names.index(ref_name)], ref_file
            )
            print_comparison(
                diff_n, big_diff_n, total, miss, diff_vars, tab_n=tab_n + 1
            )
            diff_sum += diff_n
            big_diff_sum += big_diff_n
            total_sum += total
            info[ref_name] = [diff_n, big_diff_n, total, miss, diff_vars]
        else:
            print(tab(tab_n + 1), "missing")
            miss_files.append(ref_name)
        print()
    return diff_sum, big_diff_sum, total_sum, miss_files, info


def compare_files(new_file, ref_file):
    """Compare two given files."""
    if ref_file.endswith(".txt"):
        return compare_csv_files(new_file, ref_file)
    return compare_nc_files(new_file, ref_file)


def compare_csv_files(new_file, ref_file):
    """Compare two given CSV files."""
    if "event_ids" in new_file:
      pd_new = pd.read_csv(new_file, sep=r"\s+").sort_values("<event>")
      pd_ref = pd.read_csv(ref_file, sep=r"\s+").sort_values("<event>")
    else:
      pd_new = pd.read_csv(new_file, sep=r"\s+")
      pd_ref = pd.read_csv(ref_file, sep=r"\s+")
    return compare_xarrays(pd_new.to_xarray(), pd_ref.to_xarray())


def compare_nc_files(new_file, ref_file):
    """Compare two given NC files."""
    ds_new = xr.load_dataset(new_file)
    ds_ref = xr.load_dataset(ref_file)
    return compare_xarrays(ds_new, ds_ref)


def compare_xarrays(ds_new, ds_ref, match=None, ignore=None):
    """Compare two given xarrays."""
    # renaming dictionary
    match = {} if match is None else match
    # variables to ignore
    ignore = {} if ignore is None else ignore
    diff_n = 0
    big_diff_n = 0
    total = len(ds_ref.data_vars)
    miss = []
    diff_vars = []
    for var_name in ds_ref.data_vars:
        if var_name in ignore:
            continue
        # rename the current var_name with the aid of the match-dictionary
        new_var_name = match.get(var_name, var_name)
        if new_var_name in ds_new.data_vars:
            if not ds_new[new_var_name].equals(ds_ref[var_name]):
                diff_n += 1
                if not np.allclose(
                    ds_new[new_var_name].values,
                    ds_ref[var_name].values,
                    equal_nan=True,
                    rtol=RTOL,
                    atol=ATOL,
                ):
                    diff_vars.append(str(var_name))
                    big_diff_n += 1
        else:
            miss.append(str(var_name))
            diff_n += 1
            big_diff_n += 1
    return diff_n, big_diff_n, total, miss, diff_vars


# CALL ROUTINES ###############################################################


class Output(object):
    """A class to duplicate an output stream to stdout.

    Parameters
    ----------
    file_or_name : filename or open filehandle (writable)
        File that will be duplicated
    print_log : bool, optional
        State if log should be printed. Default: True
    """

    def __init__(self, file_or_name, print_log=True):
        if hasattr(file_or_name, "write") and hasattr(file_or_name, "seek"):
            self.file = file_or_name
        else:
            self.file = open(file_or_name, "w")
        self._closed = False
        self.encoding = sys.stdout.encoding
        if not self.encoding:
            self.encoding = "utf-8"
        self.print_log = print_log
        self.last_line = ""

    def close(self):
        """Close the file and restore the channel."""
        self.flush()
        self.file.close()
        self._closed = True

    def write(self, data):
        """Write data to both channels."""
        try:
            self.last_line = data.decode(self.encoding)
        except AttributeError:
            self.last_line = data
        self.file.write(self.last_line)
        if self.print_log:
            sys.stdout.write(self.last_line)
            sys.stdout.flush()

    def flush(self):
        """Flush both channels."""
        self.file.flush()
        if self.print_log:
            sys.stdout.flush()

    def __del__(self):
        """Close and delete."""
        if not self._closed:
            self.close()


def run_model(
    exe,
    cwd=None,
    mpi_nop=0,
    openmp_threads=0,
    print_log=True,
    save_log=True,
    log_path=None,
    log_name=None,
    timeout=1800,
):
    """
    Run smi in a given directory.

    Parameters
    ----------
    exe : str
        path to the smi executable.
    cwd : str or None
        path to the working directory to run smi in. Default: os.get_cwd()
    mpi_nop : int
        Number of processes for mpi if wanted. Default: 0
    openmp_threads : int
        Number of openmp threads if wanted. Default: 0
    print_log : bool, optional
        state if the output should be displayed in the terminal.
        Default: True
    save_log : bool, optional
        Whether to save the smi output to a file.
        Default: True
    log_path : str or None, optional
        Path, where the log file should be saved. Default: None
        (the defined output directory or the task_root directory)
    log_name : str or None, optional
        Name of the log file. Default: None
        (exe+time+"_log.txt")
    timeout : int or None, optional
        Time to wait for smi to finish in seconds. Default: 1800

    Returns
    -------
    success : bool
        State if smi terminated 'normally'.
    """
    # use absolute path since we change the cwd in the smi call
    exe = os.path.abspath(exe)
    exe_name = os.path.basename(exe)
    cwd = os.getcwd() if cwd is None else os.path.abspath(cwd)
    # check that not mpi and openmp should be used at the same time
    mpi = int(mpi_nop) > 0
    if mpi and int(openmp_threads) > 0:
        raise ValueError("smi run: either use mpi or openmp, not both.")
    # set environment variable for openmp threads if wanted
    if openmp_threads > 0:
        os.environ["OMP_NUM_THREADS"] = str(int(openmp_threads))
    # define the call arguments
    args = ["mpirun", "-n", str(int(mpi_nop))] if mpi else []
    args += [exe]
    # prevent eraising files...
    if not save_log:
        log_name = None
    # set standard log_name
    if log_name is None:
        log_name = exe_name + "_" + os.path.basename(cwd) + "_log.txt"
    # put the logfile in the desired folder
    if log_path is None:
        log_path = cwd
    log = os.path.abspath(os.path.join(log_path, log_name))
    os.makedirs(os.path.dirname(log), exist_ok=True)
    # create a splitted output stream (to file and stdout)
    out = Output(log, print_log=print_log)
    # call smi with pexpect
    child = SPAWN(
        " ".join(args),
        timeout=timeout,
        logfile=out,
        encoding=out.encoding,
        cwd=cwd,
    )
    # wait for smi to finish
    child.expect(pexpect.EOF)
    if sys.platform != "win32":
        child.close()
        exitstatus = child.exitstatus
    else:
        exitstatus = child.wait()
    out.close()
    exit_ok = exitstatus == 0
    # check the last 4 lines of the log-file for "SMI: finished!"
    with open(log) as log_f:
        lines = log_f.readlines()
    smi_ok = any(["SMI: finished!" in line for line in lines[-4:]])
    success = exit_ok and smi_ok
    # give a hint about false negative exit status
    if exit_ok != smi_ok:
        print("...success confusion:", out.last_line)
    # return success status and log-file path
    if not save_log:
        os.remove(log)
        return success, None
    return success, log


# MAIN PART ###################################################################

if __name__ == "__main__":
    # get args
    exe_list, print_log, log_path, mpi_nop, omp_n, skip, only = parse_args()
    # checking path
    cases_path = os.path.dirname(os.path.realpath(__file__))
    # get all cases folders (in the cases_path)
    raw_paths = glob.glob(os.path.join(cases_path, "case*"))
    # sort out the mrm check cases
    cases = [path for path in raw_paths if "mrm" not in os.path.basename(path)]
    # sort the cases by name
    cases.sort()
    final_result = True
    final_exe_results = {}
    exe_case_results = {}
    # check "only" cases
    missing_cases = set()  # no missing cases by default
    if only is not None:
        missing_cases = set(only) - set(map(os.path.basename, cases))
        if missing_cases:
            final_result = False
    # skip some cases for mpi
    if int(mpi_nop) > 0:
        skip += SKIP_CASES_MPI
    # iterate of all smi exe-s given
    for exe in exe_list:
        # dict for checking results
        results = {}
        # result for the current exe
        exe_result = True
        # relative path to the current exe
        exe_name = os.path.relpath(exe, start=cases_path)
        print(
            sep_text(
                "checking exe: {}".format(exe_name),
                "print log: {}".format(print_log),
                "log path: {}".format(log_path),
                "mpi processes: {}".format(mpi_nop),
                "openMP threads: {}".format(omp_n),
                sep_n=70,
            )
        )
        # iterate over all cases
        for case in cases:
            # base name of the case
            case_base = os.path.basename(case)
            # check if this case should be run (defined by "only" arg)
            if only is not None and case_base not in only:
                continue
            # skip case if wanted
            if case_base in skip:
                print(sep_text("skip case: " + case_base, sep_n=60, tab_n=1))
                continue
            print(sep_text("checking case: " + case_base, sep_n=60, tab_n=1))
            # get the output and reference directories
            out_dir = os.path.join(case, OUT)
            ref_dir = os.path.join(case, REF)
            create_out_dir(out_dir)
            # run the current check case
            success, log = run_model(
                exe=exe,
                cwd=case,
                print_log=print_log,
                log_path=log_path,
                mpi_nop=mpi_nop,
                openmp_threads=omp_n,
            )
            print(tab(2), "run succeeded:", success)
            print(tab(2), "log path:", os.path.relpath(log, start=cases_path))
            print(sep_text("checking output...", sep_n=50, tab_n=2))
            # compare the output with the references
            diff_sum, big_diff_sum, total_sum, miss, info = compare_patterns(
                PATTERNS, out_dir, ref_dir, tab_n=3
            )
            print(sep_text("SUMMARY: " + case_base, sep_n=50, tab_n=2))
            print_comparison(
                diff_sum, big_diff_sum, total_sum, miss, [], tab_n=3
            )
            # define the result
            # "success" could be false-negative
            # ... = success and not big_diff_sum and not miss_files
            results[case_base] = not big_diff_sum and not miss
        # summary
        print(sep_text("CHECK SUMMARY: " + exe_name, sep_n=60, tab_n=1))
        for case_n in results:
            print(tab(2), case_n, "succeeded:", results[case_n])
            exe_result &= results[case_n]
        print(
            sep_text(
                "FINISHED: " + exe_name,
                "...succeeded: {}".format(exe_result),
                sep_n=70,
            )
        )
        # update the final result
        final_exe_results[exe_name] = exe_result
        final_result &= exe_result
        exe_case_results[exe_name] = results

    # final result
    print(sep_text("CHECK SUMMARY", sep_n=70))
    for exe_n in final_exe_results:
        print(tab(1), exe_n, "succeeded:", final_exe_results[exe_n])
        for case_n in exe_case_results[exe_n]:
            if not exe_case_results[exe_n][case_n]:
                print(tab(2), case_n, "failed")
    # print missing cases (specified by "only" arg)
    if missing_cases:
        miss_case_str = str(missing_cases)
        print(sep_text("CASES not found: " + miss_case_str, sep_n=60, tab_n=1))
    print(
        sep_text("FINISHED", "...succeeded: {}".format(final_result), sep_n=70)
    )
    # assert that the final result is positive
    if not final_result:
        raise AssertionError("smi checks were not successful.")
