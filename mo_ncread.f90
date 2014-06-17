module mo_NcRead
  
  ! This module provides subroutines for reading arrays from nc file using the netcdf4 library
  
  ! Please make sure you include the netcdf4 library in your makefile in order to use this module

  ! numerical precision
  use mo_kind, only: i4, sp, dp
  !
  ! functions and constants of netcdf4 library
  use netcdf,  only: nf90_open, nf90_get_var, nf90_close, NF90_MAX_NAME , &
                     nf90_inq_varid, nf90_inquire_variable, &
                     nf90_inquire_dimension, NF90_NOWRITE, &
                     nf90_noerr, nf90_strerror
  !
  implicit none
  !
  private
  !
  public :: Get_NcDim ! get the dimensions of a Variable
  public :: Get_NcVar ! get the data of a Variable in a nc file
  public :: check
  !
  interface Get_NcVar
     module procedure Get_NcVar_1d_sp, Get_NcVar_1d_dp, Get_NcVar_2d_sp, Get_NcVar_2d_dp, &
                      Get_NcVar_3d_sp, Get_NcVar_3d_dp, Get_NcVar_4d_sp, Get_NcVar_4d_dp, &
                      Get_NcVar_5d_sp, Get_NcVar_5d_dp, Get_NcVar_1d_i4, Get_NcVar_2d_i4, &
                      Get_NcVar_3d_i4, Get_NcVar_4d_i4, Get_NcVar_5d_i4 
  end interface
  !
contains
  !
  ! ------------------------------------------------------------------------------
  ! 
  ! NAME
  !     Get_NcDim
  !
  ! PURPOSE
  !     gets the dimensions of variable in a netcdf file
  !
  ! CALLING SEQUENCE
  !     dim = Get_NcDim(Filename, Variable, PrintInfo=PrintInfo, ndims=ndims)
  !
  ! INTENT(IN)
  !     character(len=*) :: Filename - Filename of netcdf file
  !
  ! INTENT(IN)
  !     character(len=*) :: Variable - Variable name exactly as specified in the file
  !
  ! INTENT(IN), OPTIONAL
  !     logical       :: PrintInfo - if given and true, information about dimension
  !                                  and their lengths will be printed to standard output
  !
  ! INTENT(OUT)
  !     integer(i4), dimension(5) :: Get_NcDim - dimension length, 1 if dimension does not exist
  !
  ! INTENT(OUT), OPTIONAL
  !     integer(i4) :: ndims - # of dimensions
  !
  ! LITERATURE
  !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html
  !
  ! HISTORY
  !     Written,  Stephan Thober, Dec 2011
  !     Modified, Matthias Cuntz, Jan 2012 - ndims
  ! ------------------------------------------------------------------------------

  function Get_NcDim(Filename, Variable, PrintInfo, ndims)
    !
    implicit none
    !
    character(len=*),      intent(in)  :: Filename
    character(len=*),      intent(in)  :: Variable
    logical,     optional, intent(in)  :: PrintInfo
    integer(i4), optional, intent(out) :: ndims
    integer(i4), dimension(5)          :: Get_NcDim
    !
    logical     :: PrintFlag
    integer(i4) :: ncid    ! id of input stream
    integer(i4) :: varid   ! id of variable to be read
    integer(i4) :: vartype ! type of variable
    integer(i4) :: NumDims ! # of dimensions
    !
    ! Open NetCDF filename
    call check(nf90_open(Filename, NF90_NOWRITE, ncid))
    !
    PrintFlag = .false.
    if (present(PrintInfo)) PrintFlag = PrintInfo
    !
    ! Inquire file and check if VarName exists in the dataset,
    ! get number of dimensions and
    ! get the length of the dimensions
    call Get_Info(Variable, ncid, varid, vartype, Get_NcDim, Info=PrintFlag, ndims=NumDims)
    if (present(ndims)) ndims=NumDims
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end function Get_NcDim

  ! ------------------------------------------------------------------------------

  !    NAME
  !        Get_NcVar

  !    PURPOSE
  !        Reads a 2 - 5 dimensional array from a nc file given
  !        the variable name EXACTLY as specified in the file.
  !        When calling, the array has to be allocated correctly!
  !        If the dimension of the actual data is less than the ones
  !        of the array, then the dimension lengths of the array will
  !        be filled with ones.

  !    CALLING SEQUENCE
  !        call Get_NcVar(Filename, VarName, Dat, start=jdate, count=Nvalues)

  !    INTENT(IN)
  !        character(len=*) :: Filename - Name of the nc file

  !    INTENT(IN)
  !        character(len=*) :: VarName - Name of the Variable in the nc file
  
  !    INTENT(INOUT)
  !        real(sp/dp), dimension(:,:[,:[,:[,:]]]) :: array - array where data will be read
  
  !    INTENT(IN), OPTIONAL
  !        integer(i4), dimension(:) :: jdate ! starting indeces of first value to read
  !                                             len is the number of dimensions of array, default is 1, see example 
  
  !    INTENT(IN), OPTIONAL
  !        integer(i4), dimension(:) :: count ! same size as jdate, specifies how many values in each dimension
  !                                             is going to be read
  
  !    RESTRICTIONS
  !        Output array is a floating point of 2-5 dimensions.
  !        NOT yet tested for different compilers than intel11.1.075
  !        CANNOT read packed data

  !    EXAMPLE
  !        see test program in directory test_mo_NcRead

  !    LITERATURE
  !        http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html

  !    HISTORY
  !        Written,  Stephan Thober, Nov 2011
  !        Modified, Stephan Thober, Nov 2011 - added comments
  !        Modified, Matthias Cuntz, Jan 2012 - unified routines for different dimensions and data types
  !        Modified, Stephan Thober, Mar 2012 - corrected dynamical read of data
  ! ------------------------------------------------------------------------------

  subroutine Get_NcVar_1d_sp(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 1
    integer, parameter :: itype = 5 ! 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    real(sp),    dimension(:),               intent(inout) :: Dat    ! array where values should be stored
    integer(i4), dimension(:), optional, intent(in)    :: start
    integer(i4), dimension(:), optional, intent(in)    :: count
    !
    integer(i4), dimension(5) :: Rstart
    integer(i4), dimension(5) :: Rcount
    integer(i4)                   :: ncid    ! id of input stream
    integer(i4)                   :: varid   ! id of variable to be read
    integer(i4)                   :: vartype ! type of variable
    integer(i4)                   :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = 1
    Rcount(1:idims) = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) < size(shape(dat))) stop 'ERROR*** start has less values than data has dimensions. GetNcVar'
       if (size(start) > 5) stop 'ERROR*** start has dimension greater than 5. GetNcVar'
       Rstart(1:size(start)) = start
    end if
    !
    if (present(count)) then
       if (size(count) < size(shape(dat))) stop 'ERROR*** count has less values than data has dimensions. GetNcVar'
       if (size(count) > 5) stop 'ERROR*** count has dimension greater than 5. GetNcVar'
       Rcount(1:size(count)) = count
       do i=1, idims
          if (size(Dat,i) < Rcount(i)) stop 'ERROR*** try to read more data in dimension than there is. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_1d_sp


  subroutine Get_NcVar_1d_dp(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 1
    integer, parameter :: itype = 6 ! 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    real(dp),    dimension(:),               intent(inout) :: Dat    ! array where values should be stored
    integer(i4), dimension(:), optional, intent(in)    :: start
    integer(i4), dimension(:), optional, intent(in)    :: count
    !
    integer(i4), dimension(5) :: Rstart
    integer(i4), dimension(5) :: Rcount
    integer(i4)                   :: ncid    ! id of input stream
    integer(i4)                   :: varid   ! id of variable to be read
    integer(i4)                   :: vartype ! type of variable
    integer(i4)                   :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = 1
    Rcount(1:idims) = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) < size(shape(dat))) stop 'ERROR*** start has less values than data has dimensions. GetNcVar'
       if (size(start) > 5) stop 'ERROR*** start has dimension greater than 5. GetNcVar'
       Rstart(1:size(start)) = start
    end if
    !
    if (present(count)) then
       if (size(count) < size(shape(dat))) stop 'ERROR*** count has less values than data has dimensions. GetNcVar'
       if (size(count) > 5) stop 'ERROR*** count has dimension greater than 5. GetNcVar'
       Rcount(1:size(count)) = count
       do i=1, idims
          if (size(Dat,i) < Rcount(i)) stop 'ERROR*** try to read more data in dimension than there is. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_1d_dp


  subroutine Get_NcVar_2d_sp(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 2
    integer, parameter :: itype = 5 ! 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    real(sp),    dimension(:,:),             intent(inout) :: Dat    ! array where values should be stored
    integer(i4), dimension(:), optional, intent(in)    :: start
    integer(i4), dimension(:), optional, intent(in)    :: count
    !
    integer(i4), dimension(5) :: Rstart
    integer(i4), dimension(5) :: Rcount
    integer(i4)                   :: ncid    ! id of input stream
    integer(i4)                   :: varid   ! id of variable to be read
    integer(i4)                   :: vartype ! type of variable
    integer(i4)                   :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = 1
    Rcount(1:idims) = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) < size(shape(dat))) stop 'ERROR*** start has less values than data has dimensions. GetNcVar'
       if (size(start) > 5) stop 'ERROR*** start has dimension greater than 5. GetNcVar'
       Rstart(1:size(start)) = start
    end if
    !
    if (present(count)) then
       if (size(count) < size(shape(dat))) stop 'ERROR*** count has less values than data has dimensions. GetNcVar'
       if (size(count) > 5) stop 'ERROR*** count has dimension greater than 5. GetNcVar'
       Rcount(1:size(count)) = count
       do i=1, idims
          if (size(Dat,i) < Rcount(i)) stop 'ERROR*** try to read more data in dimension than there is. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_2d_sp


  subroutine Get_NcVar_2d_dp(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 2
    integer, parameter :: itype = 6 ! 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    real(dp),    dimension(:,:),             intent(inout) :: Dat    ! array where values should be stored
    integer(i4), dimension(:), optional, intent(in)    :: start
    integer(i4), dimension(:), optional, intent(in)    :: count
    !
    integer(i4), dimension(5) :: Rstart
    integer(i4), dimension(5) :: Rcount
    integer(i4)                   :: ncid    ! id of input stream
    integer(i4)                   :: varid   ! id of variable to be read
    integer(i4)                   :: vartype ! type of variable
    integer(i4)                   :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = 1
    Rcount(1:idims) = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) < size(shape(dat))) stop 'ERROR*** start has less values than data has dimensions. GetNcVar'
       if (size(start) > 5) stop 'ERROR*** start has dimension greater than 5. GetNcVar'
       Rstart(1:size(start)) = start
    end if
    !
    if (present(count)) then
       if (size(count) < size(shape(dat))) stop 'ERROR*** count has less values than data has dimensions. GetNcVar'
       if (size(count) > 5) stop 'ERROR*** count has dimension greater than 5. GetNcVar'
       Rcount(1:size(count)) = count
       do i=1, idims
          if (size(Dat,i) < Rcount(i)) stop 'ERROR*** try to read more data in dimension than there is. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_2d_dp


  subroutine Get_NcVar_3d_sp(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 3
    integer, parameter :: itype = 5 ! 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    real(sp),    dimension(:,:,:),           intent(inout) :: Dat    ! array where values should be stored
    integer(i4), dimension(:), optional, intent(in)    :: start
    integer(i4), dimension(:), optional, intent(in)    :: count
    !
    integer(i4), dimension(5) :: Rstart
    integer(i4), dimension(5) :: Rcount
    integer(i4)                   :: ncid    ! id of input stream
    integer(i4)                   :: varid   ! id of variable to be read
    integer(i4)                   :: vartype ! type of variable
    integer(i4)                   :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = 1
    Rcount(1:idims) = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) < size(shape(dat))) stop 'ERROR*** start has less values than data has dimensions. GetNcVar'
       if (size(start) > 5) stop 'ERROR*** start has dimension greater than 5. GetNcVar'
       Rstart(1:size(start)) = start
    end if
    !
    if (present(count)) then
       if (size(count) < size(shape(dat))) stop 'ERROR*** count has less values than data has dimensions. GetNcVar'
       if (size(count) > 5) stop 'ERROR*** count has dimension greater than 5. GetNcVar'
       Rcount(1:size(count)) = count
       do i=1, idims
          if (size(Dat,i) < Rcount(i)) stop 'ERROR*** try to read more data in dimension than there is. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_3d_sp


  subroutine Get_NcVar_3d_dp(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 3
    integer, parameter :: itype = 6 ! 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    real(dp),    dimension(:,:,:),           intent(inout) :: Dat    ! array where values should be stored
    integer(i4), dimension(:), optional, intent(in)    :: start
    integer(i4), dimension(:), optional, intent(in)    :: count
    !
    integer(i4), dimension(5) :: Rstart
    integer(i4), dimension(5) :: Rcount
    integer(i4)               :: ncid    ! id of input stream
    integer(i4)               :: varid   ! id of variable to be read
    integer(i4)               :: vartype ! type of variable
    integer(i4)               :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = 1
    Rcount(1:idims) = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) < size(shape(dat))) stop 'ERROR*** start has less values than data has dimensions. GetNcVar'
       if (size(start) > 5) stop 'ERROR*** start has dimension greater than 5. GetNcVar'
       Rstart(1:size(start)) = start
    end if
    !
    if (present(count)) then
       if (size(count) < size(shape(dat))) stop 'ERROR*** count has less values than data has dimensions. GetNcVar'
       if (size(count) > 5) stop 'ERROR*** count has dimension greater than 5. GetNcVar'
       Rcount(1:size(count)) = count
       do i=1, idims
          if (size(Dat,i) < Rcount(i)) stop 'ERROR*** try to read more data in dimension than there is. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_3d_dp


  subroutine Get_NcVar_4d_sp(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 4
    integer, parameter :: itype = 5 ! 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    real(sp),    dimension(:,:,:,:),         intent(inout) :: Dat    ! array where values should be stored
     integer(i4), dimension(:), optional, intent(in)    :: start
    integer(i4), dimension(:), optional, intent(in)    :: count
    !
    integer(i4), dimension(5) :: Rstart
    integer(i4), dimension(5) :: Rcount
    integer(i4)                   :: ncid    ! id of input stream
    integer(i4)                   :: varid   ! id of variable to be read
    integer(i4)                   :: vartype ! type of variable
    integer(i4)                   :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = 1
    Rcount(1:idims) = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) < size(shape(dat))) stop 'ERROR*** start has less values than data has dimensions. GetNcVar'
       if (size(start) > 5) stop 'ERROR*** start has dimension greater than 5. GetNcVar'
       Rstart(1:size(start)) = start
    end if
    !
    if (present(count)) then
       if (size(count) < size(shape(dat))) stop 'ERROR*** count has less values than data has dimensions. GetNcVar'
       if (size(count) > 5) stop 'ERROR*** count has dimension greater than 5. GetNcVar'
       Rcount(1:size(count)) = count
       do i=1, idims
          if (size(Dat,i) < Rcount(i)) stop 'ERROR*** try to read more data in dimension than there is. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_4d_sp


  subroutine Get_NcVar_4d_dp(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 4
    integer, parameter :: itype = 6 ! 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    real(dp),    dimension(:,:,:,:),         intent(inout) :: Dat    ! array where values should be stored
    integer(i4), dimension(:), optional, intent(in)    :: start
    integer(i4), dimension(:), optional, intent(in)    :: count
    !
    integer(i4), dimension(5) :: Rstart
    integer(i4), dimension(5) :: Rcount
    integer(i4)                   :: ncid    ! id of input stream
    integer(i4)                   :: varid   ! id of variable to be read
    integer(i4)                   :: vartype ! type of variable
    integer(i4)                   :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = 1
    Rcount(1:idims) = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) < size(shape(dat))) stop 'ERROR*** start has less values than data has dimensions. GetNcVar'
       if (size(start) > 5) stop 'ERROR*** start has dimension greater than 5. GetNcVar'
       Rstart(1:size(start)) = start
    end if
    !
    if (present(count)) then
       if (size(count) < size(shape(dat))) stop 'ERROR*** count has less values than data has dimensions. GetNcVar'
       if (size(count) > 5) stop 'ERROR*** count has dimension greater than 5. GetNcVar'
       Rcount(1:size(count)) = count
       do i=1, idims
          if (size(Dat,i) < Rcount(i)) stop 'ERROR*** try to read more data in dimension than there is. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_4d_dp


  subroutine Get_NcVar_5d_sp(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 5
    integer, parameter :: itype = 5 ! 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    real(sp),    dimension(:,:,:,:,:),       intent(inout) :: Dat    ! array where values should be stored
    integer(i4), dimension(:), optional, intent(in)    :: start
    integer(i4), dimension(:), optional, intent(in)    :: count
    !
    integer(i4), dimension(5) :: Rstart
    integer(i4), dimension(5) :: Rcount
    integer(i4)                   :: ncid    ! id of input stream
    integer(i4)                   :: varid   ! id of variable to be read
    integer(i4)                   :: vartype ! type of variable
    integer(i4)                   :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = 1
    Rcount(1:idims) = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) < size(shape(dat))) stop 'ERROR*** start has less values than data has dimensions. GetNcVar'
       if (size(start) > 5) stop 'ERROR*** start has dimension greater than 5. GetNcVar'
       Rstart(1:size(start)) = start
    end if
    !
    if (present(count)) then
       if (size(count) < size(shape(dat))) stop 'ERROR*** count has less values than data has dimensions. GetNcVar'
       if (size(count) > 5) stop 'ERROR*** count has dimension greater than 5. GetNcVar'
       Rcount(1:size(count)) = count
       do i=1, idims
          if (size(Dat,i) < Rcount(i)) stop 'ERROR*** try to read more data in dimension than there is. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_5d_sp


  subroutine Get_NcVar_5d_dp(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 5
    integer, parameter :: itype = 6 ! 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    real(dp),    dimension(:,:,:,:,:),       intent(inout) :: Dat    ! array where values should be stored
    integer(i4), dimension(:), optional, intent(in)    :: start
    integer(i4), dimension(:), optional, intent(in)    :: count
    !
    integer(i4), dimension(5) :: Rstart
    integer(i4), dimension(5) :: Rcount
    integer(i4)                   :: ncid    ! id of input stream
    integer(i4)                   :: varid   ! id of variable to be read
    integer(i4)                   :: vartype ! type of variable
    integer(i4)                   :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = 1
    Rcount(1:idims) = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) < size(shape(dat))) stop 'ERROR*** start has less values than data has dimensions. GetNcVar'
       if (size(start) > 5) stop 'ERROR*** start has dimension greater than 5. GetNcVar'
       Rstart(1:size(start)) = start
    end if
    !
    if (present(count)) then
       if (size(count) < size(shape(dat))) stop 'ERROR*** count has less values than data has dimensions. GetNcVar'
       if (size(count) > 5) stop 'ERROR*** count has dimension greater than 5. GetNcVar'
       Rcount(1:size(count)) = count
       do i=1, idims
          if (size(Dat,i) < Rcount(i)) stop 'ERROR*** try to read more data in dimension than there is. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_5d_dp
  
  subroutine Get_NcVar_1d_i4(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 1
    integer, parameter :: itype = 4 ! 3 = single 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    integer(i4), dimension(:),               intent(inout) :: Dat    ! array where values should be stored
    integer(i4), dimension(idims), optional, intent(in)    :: start
    integer(i4), dimension(idims), optional, intent(in)    :: count
    !
    integer(i4), dimension(idims) :: Rstart
    integer(i4), dimension(idims) :: Rcount
    integer(i4)                   :: ncid    ! id of input stream
    integer(i4)                   :: varid   ! id of variable to be read
    integer(i4)                   :: vartype ! type of variable
    integer(i4)                   :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) /= idims) stop 'ERROR*** size of start does not equal dimensions of data. GetNcVar'
       Rstart = start
    end if
    !
    if (present(count)) then
       if (size(count) /= idims) stop 'ERROR*** size of count does not equal dimensions of data. GetNcVar'
       Rcount = count
       do i=1, idims
          if (size(Dat,i) /= Rcount(i)) stop 'ERROR*** size mismatch. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_1d_i4

  subroutine Get_NcVar_2d_i4(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 2
    integer, parameter :: itype = 4 ! 3 = single, 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    integer(i4), dimension(:,:),             intent(inout) :: Dat    ! array where values should be stored
    integer(i4), dimension(idims), optional, intent(in)    :: start
    integer(i4), dimension(idims), optional, intent(in)    :: count
    !
    integer(i4), dimension(idims) :: Rstart
    integer(i4), dimension(idims) :: Rcount
    integer(i4)                   :: ncid    ! id of input stream
    integer(i4)                   :: varid   ! id of variable to be read
    integer(i4)                   :: vartype ! type of variable
    integer(i4)                   :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) /= idims) stop 'ERROR*** size of start does not equal dimensions of data. GetNcVar'
       Rstart = start
    end if
    !
    if (present(count)) then
       if (size(count) /= idims) stop 'ERROR*** size of count does not equal dimensions of data. GetNcVar'
       Rcount = count
       do i=1, idims
          if (size(Dat,i) /= Rcount(i)) stop 'ERROR*** size mismatch. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_2d_i4

  subroutine Get_NcVar_3d_i4(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 3
    integer, parameter :: itype = 4 ! 3 = single, 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    integer(i4), dimension(:,:,:),           intent(inout) :: Dat    ! array where values should be stored
    integer(i4), dimension(idims), optional, intent(in)    :: start
    integer(i4), dimension(idims), optional, intent(in)    :: count
    !
    integer(i4), dimension(idims) :: Rstart
    integer(i4), dimension(idims) :: Rcount
    integer(i4)                   :: ncid    ! id of input stream
    integer(i4)                   :: varid   ! id of variable to be read
    integer(i4)                   :: vartype ! type of variable
    integer(i4)                   :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) /= idims) stop 'ERROR*** size of start does not equal dimensions of data. GetNcVar'
       Rstart = start
    end if
    !
    if (present(count)) then
       if (size(count) /= idims) stop 'ERROR*** size of count does not equal dimensions of data. GetNcVar'
       Rcount = count
       do i=1, idims
          if (size(Dat,i) /= Rcount(i)) stop 'ERROR*** size mismatch. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_3d_i4

  subroutine Get_NcVar_4d_i4(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 4
    integer, parameter :: itype = 4 ! 3 = single, 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    integer(i4), dimension(:,:,:,:),         intent(inout) :: Dat    ! array where values should be stored
    integer(i4), dimension(idims), optional, intent(in)    :: start
    integer(i4), dimension(idims), optional, intent(in)    :: count
    !
    integer(i4), dimension(idims) :: Rstart
    integer(i4), dimension(idims) :: Rcount
    integer(i4)                   :: ncid    ! id of input stream
    integer(i4)                   :: varid   ! id of variable to be read
    integer(i4)                   :: vartype ! type of variable
    integer(i4)                   :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) /= idims) stop 'ERROR*** size of start does not equal dimensions of data. GetNcVar'
       Rstart = start
    end if
    !
    if (present(count)) then
       if (size(count) /= idims) stop 'ERROR*** size of count does not equal dimensions of data. GetNcVar'
       Rcount = count
       do i=1, idims
          if (size(Dat,i) /= Rcount(i)) stop 'ERROR*** size mismatch. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_4d_i4

  subroutine Get_NcVar_5d_i4(Filename, VarName, Dat, start, count)
    !
    implicit none
    !
    integer, parameter :: idims = 5
    integer, parameter :: itype = 4 ! 4 = single, 5 = float, 6 = double
    !
    character(len=*),                        intent(in)    :: Filename
    character(len=*),                        intent(in)    :: VarName ! Variable name
    integer(i4), dimension(:,:,:,:,:),       intent(inout) :: Dat    ! array where values should be stored
    integer(i4), dimension(idims), optional, intent(in)    :: start
    integer(i4), dimension(idims), optional, intent(in)    :: count
    !
    integer(i4), dimension(idims) :: Rstart
    integer(i4), dimension(idims) :: Rcount
    integer(i4)                   :: ncid    ! id of input stream
    integer(i4)                   :: varid   ! id of variable to be read
    integer(i4)                   :: vartype ! type of variable
    integer(i4)                   :: i
    !
    ! Defaults for Options Start and Count
    Rstart = 1
    Rcount = shape(Dat)
    !
    ! Assign options Start and Count if present
    if (present(start)) then
       if (size(start) /= idims) stop 'ERROR*** size of start does not equal dimensions of data. GetNcVar'
       Rstart = start
    end if
    !
    if (present(count)) then
       if (size(count) /= idims) stop 'ERROR*** size of count does not equal dimensions of data. GetNcVar'
       Rcount = count
       do i=1, idims
          if (size(Dat,i) /= Rcount(i)) stop 'ERROR*** size mismatch. Get_NcVar'
       end do
    end if
    !
    ! Open NetCDF filename
    call check(nf90_open(trim(Filename),NF90_NOWRITE, ncid))
    !
    ! Inquire file, check if VarName exists and get the id
    call Get_Info(Varname,ncid,varid,vartype)
    ! check variable type ( 5 equals float type, 6 equals double )
    if (vartype /= itype) stop 'ERROR*** type of variable does not match argument type. subroutine Get_NcVar'
    !
    ! get values by varid
    call check(nf90_get_var(ncid, varid, Dat, Rstart, Rcount))
    !
    ! close File
    call check(nf90_close(ncid))
    !
  end subroutine Get_NcVar_5d_i4
  ! ------------------------------------------------------------------------------
  !
  ! SUBROUTINE GET_INFO
  !
  ! This subroutine is PRIVATE and therefore does not exist outside of this module.
  !
  ! This subroutine inquires the nc file. Given the Variable name and the stream 
  ! of the nc file, this subroutine determines the variable id, the kind of the
  ! variable and the length of the dimensions in the file.
  !
  ! See: http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90/ -> NF90-INQUIRE
  ! for detailed information.
  !
  ! ------------------------------------------------------------------------------

  subroutine Get_Info(Varname, ncid, varid, xtype, dl, Info, ndims)
    !
    implicit none
    !
    character(len=*),                    intent(in)     :: Varname
    integer(i4),                         intent(in)     :: ncid 
    integer(i4),                         intent(out)    :: varid    ! variable id of data to be read
    integer(i4),                         intent(out)    :: xtype    ! type of the variable
    integer(i4), dimension(:), optional, intent(inout)  :: dl
    logical,                   optional, intent(in)     :: Info
    integer(i4),               optional, intent(out)    :: ndims    ! Number of Dimensions for specific variable
    !
    integer(i4), dimension(:), allocatable :: DimID   ! Id of dimension
    character(NF90_MAX_NAME)               :: name    ! name of Variables in the file
    integer(i4)                            :: NumDims ! Number of Dimensions for specific variable
    integer(i4)                            :: n       ! loop index
    character(256)                         :: form    ! output format
    integer(i4)                            :: itmp
    !
    call check(nf90_inq_varid(ncid, Varname, varid))
    call check(nf90_inquire_variable(ncid, varid, ndims=NumDims))
    if (present(ndims)) ndims=NumDims
    !
    ! get the dimension Ids and the data type of the variable
    allocate(DimId(NumDims))
    call check(nf90_inquire_variable(ncid, varid, xtype=xtype, dimids=DimId))
    !
    if ( present(dl) ) then
       ! check consistency of dimensions
       if ( NumDims > size(dl) ) &
            stop 'ERROR*** Dimension size of Variable is greater than dims of array. subroutine Get_Info'
       ! go through dimension ids and get its length
       dl(:) = 1 ! initialise
       dimloop: do n = 1, NumDims
          call check(nf90_inquire_dimension(ncid, DimId(n), name, itmp))
          dl(n) = itmp
          if (present(info)) then
             if (info) then
                write(form,'(a12,I03.3,a1)') "(a10,i1,a4,a", len(trim(name)), ")"
                write(*,form) 'Dimension ', n, ' is ', trim(name)
                write(*,'(a14,i5.5)') 'The Length is ', dl(n)
             end if
          end if
       end do dimloop
       !
    end if
    !
  end subroutine Get_Info

  ! -----------------------------------------------------------------------------
  !  private error checking routine
  subroutine check(status)
    !
    implicit none
    !
    integer(i4), intent(in) :: status
    !
    if (status /= nf90_noerr) then
       write(*,*) trim(nf90_strerror(status))
       stop
    end if
    !
  end subroutine check
  !
end module mo_NcRead
