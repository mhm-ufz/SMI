!> \file mo_smi_write.f90
!> \copydoc mo_smi_write

!> \brief Writing subroutines for SMI
!> \author Luis Samaniego
!> \date 09.02.2011
module mo_smi_write

  use mo_kind,                 only: i4, sp, dp
  use mo_smi_global_variables, only: period

  implicit none

  ! ! clusters
  integer(i4), dimension(:,:,:), allocatable      :: idCluster
  integer(i4)                                     :: nInterArea
  integer(i4)                                     :: nEvents           ! total number of drough events to estimate SAD
  integer(i4), dimension(:), allocatable          :: shortCnoList      ! consolidated cluster no. list
  integer(i4), dimension(:,:), allocatable        :: eventId           ! event identification number, cluster, month of occurence
  integer(i4)                                     :: nClusters         ! number od clusters
  integer(i4), dimension(:), allocatable          :: eIdPerm           ! permutation of the event Ids with ascending area
  !
  ! cluster statistics
  real(dp), dimension(:,:), allocatable          :: DAreaEvol    ! drought area evolution (fraction Germany)
  real(dp), dimension(:,:), allocatable          :: DTMagEvol    !         magnitud
  real(dp), dimension(:), allocatable            :: aDD          ! average (mean) drought duration space-time
  real(dp), dimension(:), allocatable            :: aDA          ! average (mean) drought area
  real(dp), dimension(:), allocatable            :: TDM          ! total drought magnitud
  real(dp), dimension(:,:,:), allocatable        :: SAD          ! severity area duration curves
  real(dp), dimension(:,:), allocatable          :: SADperc      ! percentilles of severity area duration curves
  integer(i4)                                    :: nDsteps      ! number of durations steps
  real(dp), dimension(:,:,:), allocatable        :: severity     ! severity for a given duration
  real(dp), dimension(:,:,:), allocatable        :: dASevol      ! evolution of the drought areas and severity

  integer(i4), parameter                         :: nDurations = 4         ! number of durations
  integer(i4), dimension(nDurations), parameter  :: durList = (/3,6,9,12/) !(/3,6,9,12/)   ! list of durations to evaluate
  integer(i4), parameter                         :: nQProp = 3             ! number of SAD percetiles for a given duration
  real(dp), dimension(nQProp), parameter         :: QProp = (/90._dp, 96._dp, 98._dp /)
                                                                           ! percentiles corresponding
                                                                           ! to return periods of 10,25,50 years

  integer(i4)                                    :: nLargerEvents = 600    ! n. largest drought events
  !
  ! Basin summary
  integer(i4), parameter                         :: nBasins = 6

contains

  !> \brief for writting the SMI to nc file
  !> \author Stephan Thober
  !> \date 5.7.2014
  !> \date 9.8.2019
  !!       - refactored to new period structure and leap day handling
  subroutine WriteSMI( outpath, SMI, mask, per, lats_1d, lons_1d, lats_2d, lons_2d, easting, northing)

    use mo_kind,          only: i4, sp
    use mo_string_utils,  only: num2str
    use mo_constants,     only: nodata_sp, nodata_dp
    use mo_netcdf,        only: NcDataset, NcVariable, NcDimension

    implicit none

    ! input variables
    character(len=*),                              intent(in) :: outpath     ! ouutput path for results

    logical,        dimension(:,:),                intent(in) :: mask
    type(period),                                  intent(in) :: per
    real(sp),       dimension(:,:),                intent(in) :: SMI
    real(dp),       dimension(:),   allocatable,   intent(in) :: lats_1d, lons_1d   ! latitude and longitude fields of input
    real(dp),       dimension(:,:),   allocatable,   intent(in) :: lats_2d, lons_2d   ! latitude and longitude fields of input
    real(dp), dimension(:), allocatable, intent(in) :: easting ! easting coordinates of input
    real(dp), dimension(:), allocatable, intent(in) :: northing ! easting coordinates of input
    ! local Variables
    type(NcDataset)                                           :: nc_out
    type(NcVariable)                                          :: nc_var
    type(NcDimension)                                         :: nc_y, nc_x, nc_tim
    character(256)                                            :: Fname
    integer(i4)                                               :: ii          ! month

    real(sp),       dimension(:,:,:), allocatable             :: dummy_D3_sp ! field real unpacked
    integer(i4)                                               :: y
    integer(i4)                                               :: x

    ! initialize dimension
    y = size(mask, 2)
    x = size(mask, 1)
    print *,'x: ', x
    print *,'y: ', y
    ! fname
    fName = trim(outpath) //'SMI.nc'
    nc_out = NcDataset(fName, 'w')
    if (allocated(lats_1d)) then
       nc_y = nc_out%setDimension("lat", y)
       nc_x = nc_out%setDimension("lon", x)
    else if (allocated(lats_2d)) then
       nc_y = nc_out%setDimension("northing", y)
       nc_x = nc_out%setDimension("easting", x)
    else
       nc_y = nc_out%setDimension("y", y)
       nc_x = nc_out%setDimension("x", x)
    end if
    nc_tim = nc_out%setDimension("time", -1)

    ! unpack estimated SMIp
    allocate( dummy_d3_sp( x, y, size( SMI, 2) ) )
    do ii = 1,  size( SMI, 2)
      dummy_d3_sp( :, :, ii) = unpack ( SMI(:,ii), mask, nodata_sp )
    end do

    ! save SMI (same sequence as in 1)
    nc_var = nc_out%setVariable('SMI', "f32", &
         (/ nc_x, nc_y, nc_tim /))
    call nc_var%setData(dummy_d3_sp)
    call nc_var%setAttribute('long_name', 'soil moisture index')
    call nc_var%setAttribute('missing_value', nodata_sp)
    call nc_var%setAttribute('units', '-')
    deallocate( dummy_d3_sp )

    ! add lat and lon
    if (allocated(lats_1d)) then
      print *, 'lat1d'
      nc_var = nc_out%setVariable('lat', "f64", (/ nc_y /))
      call nc_var%setData(lats_1d)
      call nc_var%setAttribute('long_name', 'latitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_north')
   else if (allocated(lats_2d)) then
      print *, 'lats2d'
      nc_var = nc_out%setVariable('lat', "f64", (/ nc_x, nc_y /))
      call nc_var%setData(lats_2d)
      call nc_var%setAttribute('long_name', 'latitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_north')

      ! nc_var = nc_out%setVariable('northing', "f64", (/ nc_y/))
      ! call nc_var%setData(northing)
      ! call nc_var%setAttribute('long_name', 'northing')
      ! call nc_var%setAttribute('missing_value', nodata_dp)
      ! call nc_var%setAttribute('units', 'meters')
   end if

    if (allocated(lons_1d)) then
      nc_var = nc_out%setVariable('lon', "f64", (/ nc_x /))
      call nc_var%setData(lons_1d)
      call nc_var%setAttribute('long_name', 'longitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_east')
    else if (allocated(lons_2d)) then
      nc_var = nc_out%setVariable('lon', "f64", (/ nc_x, nc_y /))
      call nc_var%setData(lons_2d)
      call nc_var%setAttribute('long_name', 'longitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_east')

      ! nc_var = nc_out%setVariable('easting', "f64", (/ nc_x /))
      ! call nc_var%setData(easting)
      ! call nc_var%setAttribute('long_name', 'easting')
      ! call nc_var%setAttribute('missing_value', nodata_dp)
      ! call nc_var%setAttribute('units', 'meters')


    end if

    ! add time
    nc_var = nc_out%setVariable('time', "i32", (/ nc_tim /))
    call nc_var%setData(per%time_points)
    call nc_var%setAttribute('long_name', 'time')
    call nc_var%setAttribute('units', 'days since '                              // &
                                      trim(num2str(per%y_start,   '(i4)')) // '-'     // &
                                      trim(num2str(per%m_start, '(i2.2)')) // '-'     // &
                                      trim(num2str(per%d_start, '(i2.2)')) // ' 00:00:00' )
    ! close file
    call nc_out%close()

  end subroutine WriteSMI

  !> \brief for writting the CDF information (that is the kernel and associated soil moisture values) to nc file
  !> \author Stephan Thober
  !> \date 8.8.2019
  subroutine WriteCDF( outpath, SM, hh, mask, per, nCalendarStepsYear, lats_1d, lons_1d, lats_2d, lons_2d, easting, northing)

    use mo_kind,          only: i4, sp
    use mo_message,       only: message
    use mo_string_utils,  only: num2str
    use mo_constants,     only: nodata_sp, nodata_dp
    use mo_netcdf,        only: NcDataset, NcVariable, NcDimension

    implicit none

    ! input variables
    character(len=*),                              intent(in) :: outpath     ! ouutput path for results

    logical,        dimension(:,:),                intent(in) :: mask
    type(period),                                  intent(in) :: per
    real(sp),       dimension(:,:),                intent(in) :: SM
    integer(i4),                                   intent(in) :: nCalendarStepsYear
    real(dp),       dimension(:),   allocatable,   intent(in) :: lats_1d, lons_1d   ! latitude and longitude fields of input
    real(dp),       dimension(:,:),   allocatable,   intent(in) :: lats_2d, lons_2d   ! latitude and longitude fields of input
    real(dp), dimension(:), allocatable, intent(in) :: easting ! easting coordinates of input
    real(dp), dimension(:), allocatable, intent(in) :: northing ! easting coordinates of input
    real(dp),       dimension(:,:),                intent(in) :: hh

    ! local Variables
    type(NcDataset)                                           :: nc_out
    type(NcVariable)                                          :: nc_var
    type(NcDimension)                                         :: nc_y, nc_x, nc_tim, nc_cal
    character(256)                                            :: Fname
    integer(i4)                                               :: mm          ! month

    real(sp),       dimension(:,:,:), allocatable             :: dummy_D3_sp ! field real unpacked
    real(dp),       dimension(:,:,:), allocatable             :: dummy_D3_dp ! field real unpacked
    integer(i4)                                               :: y
    integer(i4)                                               :: x

    ! initialize dimension
    y = size(mask, 2)
    x = size(mask, 1)

    ! fname
    fName = trim(outpath) //'cdf_info.nc'
    nc_out = NcDataset(fName, 'w')
    if (allocated(lats_1d)) then
       nc_y = nc_out%setDimension("lat", y)
       nc_x = nc_out%setDimension("lon", x)
    else if (allocated(lats_2d)) then
       nc_y = nc_out%setDimension("northing", y)
       nc_x = nc_out%setDimension("easting", x)
    else
       nc_y = nc_out%setDimension("y", y)
       nc_x = nc_out%setDimension("x", x)
    end if
    nc_tim = nc_out%setDimension("time", -1)

    ! unpack soil moisture estimated
    allocate( dummy_d3_sp( x, y, size( SM, 2) ) )
    do mm = 1,  size( SM, 2)
      dummy_d3_sp( :, :, mm) = unpack ( SM(:,mm), mask, nodata_sp )
    end do

    ! save SM (same sequence as in 1)
    nc_var = nc_out%setVariable('SM', "f32", &
         (/ nc_x, nc_y, nc_tim /))
    call nc_var%setData(dummy_d3_sp)
    call nc_var%setAttribute('long_name', 'soil moisture saturation')
    call nc_var%setAttribute('missing_value', nodata_sp)
    call nc_var%setAttribute('units', 'fraction of pore space')
    deallocate( dummy_d3_sp )

    ! write out kernel width if it has been optimised
    allocate( dummy_D3_dp( x, y, size( hh, 2 ) ) )
    do mm = 1, size( hh, 2 )
      dummy_D3_dp(:, :, mm) = unpack( hh(:,mm), mask, real( nodata_sp, dp ) )
    end do
    nc_cal = nc_out%setDimension("time_steps", size(dummy_d3_dp,3))
    nc_var = nc_out%setVariable('kernel_width', "f64", &
        (/ nc_x, nc_y, nc_cal /))
    call nc_var%setData(dummy_d3_dp)
    call nc_var%setAttribute('long_name', 'optimised kernel width')
    call nc_var%setAttribute('missing_value', nodata_dp)
    if (nCalendarStepsYear .eq. 12_i4) then
      call nc_var%setAttribute('units', 'months')
    else if (nCalendarStepsYear .eq. 365_i4) then
      call nc_var%setAttribute('units', 'days')
    else
      call message("***ERROR: nCalendarStepsYear has to be 12 or 365 in subroutine WriteCDF")
      stop 1
    end if

    deallocate( dummy_D3_dp )

    ! add lat and lon
    if (allocated(lats_1d)) then
      nc_var = nc_out%setVariable('lat', "f64", (/ nc_y /))
      call nc_var%setData(lats_1d)
      call nc_var%setAttribute('long_name', 'latitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_north')
    else if (allocated(lats_2d)) then
      nc_var = nc_out%setVariable('lat', "f64", (/ nc_x, nc_y /))
      call nc_var%setData(lats_2d)
      call nc_var%setAttribute('long_name', 'latitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_north')

      ! nc_var = nc_out%setVariable('northing', "f64", (/ nc_y /))
      ! call nc_var%setData(northing)
      ! call nc_var%setAttribute('long_name', 'northing')
      ! call nc_var%setAttribute('missing_value', nodata_dp)
      ! call nc_var%setAttribute('units', 'meters')
   end if

    if (allocated(lons_1d)) then
      nc_var = nc_out%setVariable('lon', "f64", (/ nc_x /))
      call nc_var%setData(lons_1d)
      call nc_var%setAttribute('long_name', 'longitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_east')

    else if (allocated(lons_2d)) then
      nc_var = nc_out%setVariable('lon', "f64", (/ nc_x, nc_y /))
      call nc_var%setData(lons_2d)
      call nc_var%setAttribute('long_name', 'longitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_east')

      ! nc_var = nc_out%setVariable('easting', "f64", (/ nc_x /))
      ! call nc_var%setData(easting)
      ! call nc_var%setAttribute('long_name', 'easting')
      ! call nc_var%setAttribute('missing_value', nodata_dp)
      ! call nc_var%setAttribute('units', 'meters')
    end if

    ! add time
    nc_var = nc_out%setVariable('time', "i32", (/ nc_tim /))
    call nc_var%setData(per%time_points)
    call nc_var%setAttribute('long_name', 'time')
    call nc_var%setAttribute('units', 'days since '                              // &
                                      trim(num2str(per%y_start,   '(i4)')) // '-'     // &
                                      trim(num2str(per%m_start, '(i2.2)')) // '-'     // &
                                      trim(num2str(per%d_start, '(i2.2)')) // ' 00:00:00' )
    ! close file
    call nc_out%close()

  end subroutine WriteCDF

  !> \brief WRITE netCDF files
  !!        FORMAT     netCDF
  !!        http://www.unidata.ucar.edu/software/netcdf/
  !> \author Luis E. Samaniego-Eguiguren, UFZ
  !> \date 16.02.2011
  subroutine WriteNetCDF(outpath, wFlag, per, lats_1d, lons_1d, lats_2d, lons_2d, easting, northing, &
        SMIc, SM_invert, duration)
    !
    use mo_kind,          only: i4
    use mo_string_utils,  only: num2str
    use mo_constants,     only: nodata_i4 , nodata_dp, nodata_sp
    use mo_netcdf,        only: NcDataset, NcVariable, NcDimension

    implicit none
    !
    ! input variables
    character(len=*),                            intent(in) :: outpath     ! ouutput path for results
    type(period),                                intent(in) :: per
    integer(i4),                                 intent(in) :: wFlag
    real(dp),       dimension(:),   allocatable,   intent(in) :: lats_1d, lons_1d   ! latitude and longitude fields of input
    real(dp),       dimension(:,:),   allocatable,   intent(in) :: lats_2d, lons_2d   ! latitude and longitude fields of input
    real(dp), dimension(:), allocatable, intent(in) :: easting ! easting coordinates of input
    real(dp), dimension(:), allocatable, intent(in) :: northing ! easting coordinates of input
    integer(i4),    dimension(:,:,:), optional,  intent(in) :: SMIc         ! Drought indicator
    real(sp),       dimension(:,:,:), optional,  intent(in) :: SM_invert
    integer(i4),                      optional,  intent(in) :: duration     ! optional, duration

    ! local Variables
    type(NcDataset)              :: nc_out
    type(NcVariable)             :: nc_var
    type(NcDimension)            :: nc_y, nc_x, nc_tim
    character(256)               :: Fname

    select case (wFlag)

    case (2)
       ! fname
       fName  = trim(outpath)//'SM_invert.nc'
       nc_out = NcDataset(fname, 'w')

       if (allocated(lats_1d)) then
          nc_y = nc_out%setDimension("lat", size(SM_invert, 2))
          nc_x = nc_out%setDimension("lon", size(SM_invert, 1))
       else if (allocated(lats_2d)) then
          nc_y = nc_out%setDimension("northing", size(SM_invert, 2))
          nc_x = nc_out%setDimension("easting", size(SM_invert, 1))
       else
          nc_y = nc_out%setDimension("y", size(SM_invert, 2))
          nc_x = nc_out%setDimension("x", size(SM_invert, 1))
       end if
       nc_tim = nc_out%setDimension("time", -1)

       nc_var = nc_out%setVariable('SM_Lall', "f32", &
            (/ nc_x, nc_y, nc_tim /))
       call nc_var%setData(SM_invert)
       call nc_var%setAttribute('long_name', 'SM according to inverse of SMI')
       call nc_var%setAttribute('missing_value', nodata_sp)
       call nc_var%setAttribute('units', 'mm')


    case (3)
       ! fname
       fName  = trim(outpath)//'SMIc.nc'
       nc_out = NcDataset(fname, 'w')
       print *, 'y SMIc: ',size(SMIc, 2)
       print *, 'x SMIc: ',size(SMIc, 1)
       if (allocated(lats_1d)) then
          nc_y = nc_out%setDimension("lat", size(SMIc, 2))
          nc_x = nc_out%setDimension("lon", size(SMIc, 1))
       else if (allocated(lats_2d)) then
          nc_y = nc_out%setDimension("northing", size(SMIc, 2))
          nc_x = nc_out%setDimension("easting", size(SMIc, 1))
       else
          nc_y = nc_out%setDimension("y", size(SMIc, 2))
          nc_x = nc_out%setDimension("x", size(SMIc, 1))
       end if
       nc_tim = nc_out%setDimension("time", -1)

       nc_var = nc_out%setVariable('mSMIc', "i32", &
            (/ nc_x, nc_y, nc_tim /))
       call nc_var%setData(SMIc)
       call nc_var%setAttribute('long_name', 'monthly SMI indicator SMI < th')
       call nc_var%setAttribute('missing_value', nodata_i4)
       call nc_var%setAttribute('units', '-')

    case (4)
       ! fname
       fName = trim(outpath)//'DCluster.nc'
       nc_out = NcDataset(fname, 'w')

       if (allocated(lats_1d)) then
          nc_y = nc_out%setDimension("lat", size(idCluster, 2))
          nc_x = nc_out%setDimension("lon", size(idCluster, 1))
       else if (allocated(lats_2d)) then
          nc_y = nc_out%setDimension("northing", size(idCluster, 2))
          nc_x = nc_out%setDimension("easting", size(idCluster, 1))
       else
          nc_y = nc_out%setDimension("y", size(idCluster, 2))
          nc_x = nc_out%setDimension("x", size(idCluster, 1))
       end if
       nc_tim = nc_out%setDimension("time", -1)

       nc_var = nc_out%setVariable('mDC', "i32", &
            (/ nc_x, nc_y, nc_tim /))
       call nc_var%setData(idCluster)
       call nc_var%setAttribute('long_name', 'consolidated cluster evolution')
       call nc_var%setAttribute('missing_value', nodata_i4)
       call nc_var%setAttribute('units', '-')

    case (5)
       ! fname
       write (fName, '(i2.2)') duration
       fName = trim(outpath)//'sev_'//trim(fName)//'.nc'
       nc_out = NcDataset(fname, 'w')

       if (allocated(lats_1d)) then
          nc_y = nc_out%setDimension("lat", size(severity, 2))
          nc_x = nc_out%setDimension("lon", size(severity, 1))
       else if (allocated(lats_2d)) then
          nc_y = nc_out%setDimension("northing", size(severity, 2))
          nc_x = nc_out%setDimension("easting", size(severity, 1))
       else
          nc_y = nc_out%setDimension("y", size(severity, 2))
          nc_x = nc_out%setDimension("x", size(severity, 1))
       end if
       nc_tim = nc_out%setDimension("time", -1)

       nc_var = nc_out%setVariable('Severity', "f64", &
            (/ nc_x, nc_y, nc_tim /))
       call nc_var%setData(idCluster)
       call nc_var%setAttribute('long_name', 'd-month severity')
       call nc_var%setAttribute('missing_value', nodata_dp)
       call nc_var%setAttribute('units', '-')

    end select

    ! add lat and lon
    if (allocated(lats_1d)) then
      nc_var = nc_out%setVariable('lat', "f64", (/ nc_y /))
      call nc_var%setData(lats_1d)
      call nc_var%setAttribute('long_name', 'latitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_north')
    else if (allocated(lats_2d)) then
      nc_var = nc_out%setVariable('lat', "f64", (/ nc_x, nc_y /))
      call nc_var%setData(lats_2d)
      call nc_var%setAttribute('long_name', 'latitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_north')

      ! nc_var = nc_out%setVariable('northing', "f64", (/ nc_y /))
      ! call nc_var%setData(northing)
      ! call nc_var%setAttribute('long_name', 'northing')
      ! call nc_var%setAttribute('missing_value', nodata_dp)
      ! call nc_var%setAttribute('units', 'meters')
   end if

    if (allocated(lons_1d)) then
       nc_var = nc_out%setVariable('lon', "f64", (/ nc_x /))
      call nc_var%setData(lons_1d)
      call nc_var%setAttribute('long_name', 'longitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_east')
    else if (allocated(lons_2d)) then
      nc_var = nc_out%setVariable('lon', "f64", (/ nc_x, nc_y /))
      call nc_var%setData(lons_2d)
      call nc_var%setAttribute('long_name', 'longitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_east')

      ! nc_var = nc_out%setVariable('easting', "f64", (/ nc_x /))
      ! call nc_var%setData(easting)
      ! call nc_var%setAttribute('long_name', 'easting')
      ! call nc_var%setAttribute('missing_value', nodata_dp)
      ! call nc_var%setAttribute('units', 'meters')
    end if

    ! add time
    if ((wflag .eq. 2) .or. &
        (wflag .eq. 3) .or. &
        (wflag .eq. 4)) then
      nc_var = nc_out%setVariable('time', "i32", (/ nc_tim /))
      call nc_var%setData(per%time_points)
      call nc_var%setAttribute('long_name', 'time')
      call nc_var%setAttribute('units', 'days since '                              // &
          trim(num2str(per%y_start,   '(i4)')) // '-'     // &
          trim(num2str(per%m_start, '(i2.2)')) // '-'     // &
          trim(num2str(per%d_start, '(i2.2)')) // ' 00:00:00' )
    end if
    ! close file
    call nc_out%close()
  end subroutine WriteNetCDF

!> \brief WRITE Results of the cluster analysis
!> \author Luis E. Samaniego-Eguiguren, UFZ
!> \date 17.03.2011
subroutine WriteResultsCluster(SMIc, outpath, wFlag, yStart, yEnd, nMonths, nCells, &
     deltaArea, cellsize, d)

  use mo_kind,          only : i4, dp
  use mo_constants,     only : YearMonths

  !
  implicit none

  ! input variables
  character(len=*),      intent(in)        :: outpath     ! ouutput path for results
  integer(i4), dimension(:,:,:),intent(in) :: SMIc        ! Drought indicator 1 - is under drought
  integer(i4),           intent(in)        :: wFlag
  integer(i4),           intent(in)        :: yStart
  integer(i4),           intent(in)        :: yEnd
  integer(i4),           intent(in)        :: nMonths
  integer(i4),           intent(in)        :: nCells
  integer(i4),           intent(in)        :: deltaArea
  real(sp),              intent(in)        :: cellsize
  integer(i4), optional, intent(in)        :: d
  real(dp), parameter                      :: eps = 1.0e-5_dp ! EPSILON(1.0_dp)


  ! local variables
  real(dp)                  :: pDArea
  !
  integer(i4)               :: i, j, t, k, y, m, mStart, mEnd
  character(len=256)        :: FMT
  character(len=256)        :: fName

  select case (wFlag)
    case (1)
      ! main statistics
      fName = trim(outpath) // 'results_ADM.txt'
      open  (10, file = fName, status='unknown')
      write (10, 100 ) 'i', 'c_Id', 'mStart', 'mEnd', 'aDD', 'aDA', 'TDM'
      do i = 1,nClusters
         ! find starting and ending months of every cluster
         mStart = 0
         mEnd = 0
         do t = 1, nMonths
            if ( DAreaEvol(t,i) .gt. 0.0_dp) then
               mEnd = t
               if (mStart .eq. 0) mStart = t
            end if
            if ( ( DAreaEvol(t,i) .lt. eps) .and. (mStart .gt. 0) ) exit
         end do
        write (10,110)  i, shortCnoList(i), mStart, mEnd, aDD(i), aDA(i), TDM(i)
      end do
      close (10)
      !
      fName = trim(outpath) // 'DArea_evol.txt'
      open  (11, file = fName, status='unknown')
      write (FMT, 120) '(a5,', nClusters, '(2x,a2,i7.7))'
      write (11, FMT ) 'm', ('c_', shortCnoList(i), i=1,nClusters)
      write (FMT, 130) '(i5,', nClusters, 'e11.5)'
      do t=1, nMonths
        write (11,FMT) t, (DAreaEvol(t,i), i=1,nClusters)
      end  do
      close (11)
      !
      fName = trim(outpath) // 'TDM_evol.txt'
      open  (12, file = fName, status='unknown')
      write (FMT, 120) '(a5,', nClusters, '(2x,a2,i7.7))'
      write (12, FMT ) 'm', ('c_', shortCnoList(i), i=1,nClusters)
      write (FMT, 130) '(i5,', nClusters, 'e11.5)'
      do t=1, nMonths
        write (12,FMT) t, (DTMagEvol(t,i), i=1,nClusters)
      end  do
      close (12)
      !
      fName = trim(outpath) // 'event_ids.txt'
      open  (13, file = fName, status='unknown')
      write (13, 140 ) '<event>', 'c_Id', 'month', 'nCells'
      ! sorted events in ascending order of areal extend
      do i=1, nEvents
        write (13,145) eIdPerm(i), (eventId(eIdPerm(i),j), j=1,3)
      end  do
      close (13)
      !
      ! NEW STATISTICS
      ! total drought area evolution (% w.r.t. whole domain)
      fName = trim(outpath) // 'DArea_evol_total.txt'
      open  (14, file = fName, status='unknown')
      write (14, 150 ) 'year', 'month', '%AreaEU'
      t = 0
      do y =yStart, yEnd
        do m = 1, int(YearMonths, i4)
          t = t + 1
          pdArea = real( count( SMIc(:,:,t) == 1 ), dp) / &
                   real( nCells, dp) * 1e2_dp
          write(14,160) y, m, pdArea
        end do
      end do
      close (14)
      !
      ! Time evolution of the cluster area (less than SMIc) and monthly severity
      fName = trim(outpath) // 'DcArea_sev_evol.txt'
      open  (15, file = fName, status='unknown')
      write (15, 170 ) 'year', 'month', '%cAreaEU','SevDE'
      t = 0
      do y =yStart, yEnd
        do m = 1, int(YearMonths, i4)
          t = t + 1
          write(15,180) y, m, dASevol(t,1,nBasins+1), dASevol(t,2,nBasins+1)
        end do
      end do
      close (15)
      !
    case(2)
      ! SAD curves for each event and duration
      do k = 1, nLargerEvents
         write (fName,200) 'SAD_e_',  eIdPerm( nEvents + 1 - k ) , '_', d, '.txt'
         fName = trim(outpath) // trim(fName)
         open (20, file=fName, status='unknown')
         write (20,210) 'Area[km2]', 'Severity'
         write (20,220) (( SAD(i, j, k ), j=1,2 ), i=1, nInterArea)
         close(20)
      end do

      ! SAD percentiles for a given duration
      write (fName,230) 'SAD_perc_', d, '.txt'
      fName = trim(outpath) // trim(fName)
      open (21, file=fName, status='unknown')
      write (FMT, 240) '(a15,', nQProp, '(9x,a2,f4.2))'
      write (21,FMT) 'Area[km2]', ('p_', QProp(i), i=1,nQProp)
      write (FMT, 250) '(f15.0,', nQProp, 'f15.5)'
      write (21,FMT) ( real(i * deltaArea, dp)*(real(cellsize,dp)**2.0_dp), (SADperc(i,j), j=1,nQProp), i=1, nInterArea)

      close (21)
  end select

  100 format (2a15, 2a10, 3a15)
  110 format (2i15, 2i10, 3f15.5)
  120 format (a4,i5,a13)
  130 format (a4,i5,a6)
  140 format (4a12)
  145 format (4i12)
  150 format (2a10,a10)
  160 format (2i10,f10.3)
  170 format (2a10,2a10)
  180 format (2i10,2f10.3)

  200 format (a6,i5.5,a,i2.2,a4)
  210 format (2a15)
  220 format (f15.0, f15.5)
  230 format (a9,i5.5,a,i2.2,a4)
  240 format (a5,i2,a13)
  250 format (a7,i2,a6)
end subroutine WriteResultsCluster


!> \brief WRITE Results of the cluster SMI basins
!> \author Luis E. Samaniego-Eguiguren, UFZ
!> \date 24.05.2011
subroutine WriteResultsBasins( outpath, SMI, mask, yStart, yEnd, nMonths, Basin_Id )

  use mo_kind,          only : i4
  use mo_constants,     only : nodata_dp, YearMonths

  implicit none

  ! input variables
  character(len=*),            intent(in) :: outpath     ! ouutput path for results
  real(sp),    dimension(:,:), intent(in) :: SMI
  logical,     dimension(:,:), intent(in) :: mask
  integer(i4),                 intent(in) :: yStart
  integer(i4),                 intent(in) :: yEnd
  integer(i4),                 intent(in) :: nMonths     ! number of simulated months
  integer(i4), dimension(:,:), intent(in) :: Basin_Id    ! IDs for basinwise drought analysis

  ! local variables
  real(dp),    dimension(size(mask,1), &
                        size(mask,2), nMonths) :: SMI_unpack
  character(len=256)                           :: fName
  integer(i4)                                  :: i, m, y, j
  real(dp),    dimension(:,:), allocatable     :: Basin_SMI

  !-------------------------
  ! basin wise
  !-------------------------
  allocate( Basin_SMI (nMonths, nBasins+1) )
  ! unpack SMI
  do i = 1, size(SMI,2)
     SMI_unpack(:,:,i) = unpack(real(SMI(:,i), dp), mask, nodata_dp)
  end do

  Basin_SMI = nodata_dp
  !
  do m = 1, nMonths
    do i = 1, nBasins
      !
      Basin_SMI(m,i) = sum( SMI_unpack(:,:,m), Basin_Id(:,:) == i ) / real(count(Basin_Id(:,:) == i), dp)
      !
    end do
    Basin_SMI(m,nBasins+1) = sum(SMI_unpack(:,:,m), mask ) /  &
         real(count(mask), dp)
  end do
  !
  ! Print Results
  fName = trim(outpath) // 'basin_avg_SMI.txt'
  open(20, file =fName , status = 'unknown')
  write(20, 1) 'Month_No', 'year', 'month', ('Basin_', i, i = 1, nBasins), 'Germany_Avg'
  !
  fName = trim(outpath) // 'basin_avg_dArea.txt'
  open(21, file =fName , status = 'unknown')
  write(21, 3) 'Month_No', 'year', 'month', ('Basin_', i, i = 1, nBasins)
  !
  fName = trim(outpath) // 'basin_avg_sev.txt'
  open(22, file =fName , status = 'unknown')
  write(22, 3) 'Month_No', 'year', 'month', ('Basin_', i, i = 1, nBasins)
  !! NEW !!
  m=0
  do y = yStart, yEnd
     do j = 1, int(YearMonths, i4)
        m = m + 1
        write(20, 2) m, y, j, (Basin_SMI(m,i), i = 1, nBasins), Basin_SMI(m,nBasins+1)
        write(21, 4) m, y, j, (dASevol(m,1,i), i = 1, nBasins)
        write(22, 4) m, y, j, (dASevol(m,2,i), i = 1, nBasins)
     end do
  end do
  close(20)
  close(21)
  close(22)
  !
  1 format(3a8, 2x, 6(a6, i2.2, 2x),  a13 )
  2 format(3i8, 2x, 6(f8.4,     2x),  f13.4)
  3 format(3a8, 2x, 6(a6, i2.2, 2x) )
  4 format(3i8, 2x, 6(f8.4,     2x) )

end subroutine WriteResultsBasins

end module mo_smi_write
