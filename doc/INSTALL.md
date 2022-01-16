# Install

## Dependencies

For Windows, some Linux distributions and soon also MacOS
specific installation instructions for the following list
can be found below.

- a fortran compiler
- make (a tool to compile a program)
- cmake (version >= 3.14) (a tool to create a system dependent makefile)
- fitting netcdf-fortran libraries (libraries for the usage of the data format netcdf on which smi depends)
- (optional, but makes things much easier) git

Git is a version-control system. If you want to contribute to a project, it is highly recommended to
use Git. You can use Git to download (clone) the project to your local pc and have a look at the history or
synchronize it without copying the whole repository again. You can also download the project folder without
Git, but this would not allow you to pull updates from and push changes to our repository.

## System dependend installation instructions

### Windows
[Cygwin](https://cygwin.com/) is an environment with a terminal that allows to compile and
run programs of Unix-like systems. You can find further instructions to install cygwin on the webpage, as well as
instructions on how to install further dependencies after the installation.

After the installation of cygwin and its dependencies smi will be installed
using cygwin. All commands and the execution of smi only run in that environment.

Install cygwin by executing the cygwin setup and choose the following dependencies:

- [ ] gcc-fortran (the fortran compiler)
- [ ] make
- [ ] cmake (version >= 3.14)
- [ ] libnetcdf-fortran-devel
- [ ] Git *(optional, Git is also available outside of cygwin, [see the Git website](https://git-scm.com/downloads))*

While installing cygwin you will have to choose a mirror. A mirror is a server
on the internet where the files for the installation come from. Choose any server
located near your city and when in doubt, choose the first one in the list.
In the next step you can find all available packages provided by cygwin, set
the view to "full". In the search panel you can filter the packages
by the dependencies listed above (e.g. make). When you choose a
version, the newest one is usually a good choice if not marked as experimental.

*Note for UFZ members:* Install cygwin locally, do not choose a location on the
network for the installation.

Some cygwin versions create a new home directory for you. You may check e.g. here:

    C:\cygwin64\home\$username


### Ubuntu, WSL2, Mint and other apt-get based systems with matching repositories

    sudo apt-get install git # (optional)
    sudo apt-get install gfortran netcdf-bin libnetcdf-dev libnetcdff-dev cmake

### Archlinux

    sudo pacman -S git # (optional)
    sudo pacman -S gcc-libs netcdf-fortran cmake

### Module systems

If you are on a module system, load the modules gcc or intel depending on your
favorite compiler. Then, load the modules netcdf-fortran and cmake.

These modules will have system specific names, environments, etc.
You may use `module spider` to find the right packages and the
right dependencies, potentially use corresponding wiki pages.

#### On eve (the cluster at the UFZ)

From the source directory use a script provided in `moduleLoadScripts`,
for example for the GNU 7.3 compiler:

    source hpc-module-loads/eve.gfortran102

### Conda (on MacOS and other Unix systems supporting Conda)

Conda provides all necessary dependencies. Here we create a local conda environment
```bash
conda create -y --prefix ./smi_env
conda activate ./smi_env
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install -y cmake make fortran-compiler netcdf-fortran fypp
```

Then follow the compile instructions bellow.

## Specific setups

The following hints can replace the step `cmake` in the installation instruction.

You can skip this part and continue with "Installation", if you do not have a module system
setup (like on clusters) or if you have not installed all packages with a package manager,
such as cygwin or apt-get.

### Module systems

The executable can be build in a way that it runs independend of loaded modules in the end. The
module system, though, adds system paths in the backround the user should not care about too much, so
the setup is a workaround. (This would be the case with any other building tool aswell.)
It should be stable, anyway.

In case you want to have a module-independend build, instead of just executing `cmake`, either run

    cmake -DCMAKE_BUILD_MODULE_SYSTEM_INDEPENDEND:STRING=ON

or change the variable `CMAKE_BUILD_MODULE_SYSTEM_INDEPENDEND` with `ccmake` to `ON` after running `cmake`.

## Installation

1. Change to a directory where you want to store the source code.
2. Clone the corresponding smi repository into a folder using Git (if installed):

        git clone git@git.ufz.de:chs/progs/SMI.git
        cd SMI

3. Configure the build and generate a system dependent makefile

    Execute `cmake` with the path to the build folder (`-B`, folder will be created) and the Git source directory (`-S`) as parameter.

        cmake -B build -S .

    If everything worked well a Makefile was created with the corresponding paths.

    *Note: have a look at "Specific setups" above in case you are using module systems,
    or when the netcdf libraries are not located where the package manager usually installs libraries,
    or when they are not saved in environment variables (i.e., classical MacOS setups at CHS).*

4. Make the build:

    You can now use cmake to build the SMI program, which will use `make` internally:

        cmake --build build --parallel

    If this also worked fine, an executable was created, which has to be moved or copied to the Git source directory.

5. Execute the file:

        cp build/app/smi .

    Now you might execute smi:

        ./smi

6. Installation:

    In order to install the compiled `smi` program to access it system-wide, you can run the following:

        cmake --install build --prefix <your/install/prefix>

    `<your/install/prefix>` needs to be replaced with a location on you computer. For example:

    - within a conda-environment: `$CONDA_PREFIX`
    - local installation for current user: `~/.local`

## Building Realease or Debug versions

If you want to set up specific versions of the build, you can
create different folders for that. Assume a release and a debug
version. Then a good idea would be to create one folder named
`debug` and one folder named `release`

    cmake -DCMAKE_BUILD_TYPE=Release -B release -S .

    cmake -DCMAKE_BUILD_TYPE=Debug -B debug -S .

Executing

    cmake --build release --parallel

    cmake --build debug --parallel

would then always result in a release build and a debug build in the respective folder.
