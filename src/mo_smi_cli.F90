!> \file    mo_smi_cli.f90
!> \copydoc mo_smi_cli

!> \brief   Module to parse command line arguments of SMI.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Oct 2021
!> \details A simple parser for command line arguments for smi.
!!          You can pass the path to the config namelist (main.dat by default)
!!
!!          You can also pass the CWD as plain last argument and get a help or version text.
!!
!!          \code{.sh}
!!          smi -h
!!          smi -v
!!          \endcode
module mo_smi_cli

#ifdef NAG
  USE f90_unix_dir, ONLY : CHDIR
#endif

  implicit none

  private

  public :: parse_command_line

contains

  !> \brief parse the given command line arguments.
  subroutine parse_command_line()
    use mo_cli, only: cli_parser
    use mo_smi_info, only: version, file_namelist

    implicit none

    type(cli_parser) :: parser

    parser = cli_parser( &
      description="The Soil Moisture Index - SMI program", &
      add_version_option=.true., &
      version=version)

    call parser%add_option( &
      name="cwd", &
      blank=.true., &
      help="The desired working directory (optional).")

    call parser%add_option( &
      name="nml", &
      s_name="n", &
      has_value=.true., &
      value_name="path", &
      default="main.dat", &
      help="The SMI configuration namelist.")

    ! parse given command line arguments
    call parser%parse()

    ! change working directory first
    if (parser%option_was_read("cwd")) call chdir(parser%option_value("cwd"))

    ! set nml file path
    file_namelist = parser%option_value("nml")

  end subroutine parse_command_line

end module mo_smi_cli
