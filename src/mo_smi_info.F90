!> \file    mo_smi_info.f90
!> \copydoc mo_smi_info

#ifndef PROJVERSION
#define PROJVERSION "0.0.0-dev0"
#endif
#ifndef PROJDATE
#define PROJDATE __DATE__
#endif
#define set_version(x) character(len = *), parameter :: version = x
#define set_date(x) character(len = *), parameter :: version_date = x

!> \brief   module with SMI program information
!> \details Provides all information about the SMI program as parameter.
!!          The \p version parameter will be set during compilation
!!          to content of \a version.txt.
!!          The \p version_date parameter will be set during compilation
!!          to content of \a version_date.txt,
!!          if it is a release version, otherwise it will be the current date.
!> \author  Sebastian Mueller
!> \date    Jan 2022
module mo_smi_info

  implicit none

  set_version(PROJVERSION)
  !< Current program version

  set_date(PROJDATE)
  !< Time of current program version release

  !> Driver file
  character(len = *), parameter :: file_main = 'main.f90'
  !> Namelist file name
  character(:), allocatable :: file_namelist ! = 'main.dat'

end module mo_smi_info
