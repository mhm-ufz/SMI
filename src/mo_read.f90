!**********************************************************************
!  MODULE InputOutput                                                 *
!  PURPOSE     Input / output subroutines                             *
!  CREATED     Luis Samaniego, 09.02.2011                             *
!                                                                     *
!**********************************************************************
MODULE mo_read 

  USE mo_kind,   only : i4, sp, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ReadDataMain   ! read everything

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------


  !*********************************************************************
  !    SUBROUTINE Read Database
  !    PURPOSE    Reads main file
  !    AUTHOR:    Luis E. Samaniego-Eguiguren, UFZ 23.05.2007
  !    NOTES:
  !               packed fields are stored in dim1->dim2 sequence 
  !*********************************************************************
  subroutine ReadDataMain( SMI, do_cluster, cluster_ext_smi, eval_SMI, read_opt_h, silverman_h, opt_h, lats, lons,  &
                           basin_flag, mask, SM_est,  SM_eval, yStart, yEnd, mStart, dStart,  &
                           Basin_Id, timepoints, SMI_thld, outpath, cellsize, thCellClus, nCellInter, do_sad, deltaArea, &
                           nCalendarStepsYear )    !  tmask_eval, tmask_est,

    use mo_kind,          only: i4
    use mo_utils,         only: equal, notequal
    use mo_ncread,        only: Get_NcVar, Get_NcVarAtt

    use mo_smi_constants, only: nodata_dp, YearMonths

    USE mo_julian,        only: NDYIN, NDAYS

    
    implicit none

    ! input / output Variables
    real(sp),    dimension(:,:), allocatable, intent(out) :: SMI         ! SMI only read if cluster_ext_smi is TRUE
    logical,                                  intent(out) :: do_sad      ! do SAD analysis
    logical,                                  intent(out) :: do_cluster  ! do cluster calculation
    logical,                                  intent(out) :: cluster_ext_smi ! clsutering external data
    logical,                                  intent(out) :: eval_SMI    ! should SMI be calculated
    logical,                                  intent(out) :: read_opt_h  ! read kernel width
    logical,                                  intent(out) :: silverman_h ! optimize kernel width
    logical,                                  intent(out) :: basin_flag  ! basin flag
    logical,     dimension(:,:), allocatable, intent(out) :: mask        ! grid mask
    integer(i4),                              intent(out) :: yStart      ! starting year
    integer(i4),                              intent(out) :: yEnd        ! ending year
    integer(i4),                              intent(out) :: mStart      ! starting month
    integer(i4),                              intent(out) :: dStart      ! starting day
    integer(i4),                              intent(out) :: thCellClus  ! treshold  for cluster formation in space ~ 640 km2
    integer(i4),                              intent(out) :: nCellInter  ! number cells for joining clusters in time ~ 6400 km2
    integer(i4),                              intent(out) :: deltaArea   ! number of cells per area interval
    integer(i4),                              intent(out) :: nCalendarStepsYear  ! number of calendar time steps per
                                                                                 ! year (month=12, day=365)
    
    integer(i4), dimension(:,:), allocatable, intent(out) :: Basin_Id    ! IDs for basinwise drought analysis
    integer(i4), dimension(:),   allocatable, intent(out) :: timepoints  ! timestep as from input NetCDF
    real(sp),                                 intent(out) :: cellsize    ! cell edge lenght of input data
    real(sp),    dimension(:,:), allocatable, intent(out) :: SM_est      ! daily / monthly fields packed for estimation
!    logical,     dimension(:,:), allocatable, intent(out) :: tmask_est   ! temporal mask of estimated arr
    real(sp),    dimension(:,:), allocatable, intent(out) :: SM_eval     ! monthly fields packed for evaluation
!    logical,     dimension(:,:), allocatable, intent(out) :: tmask_eval  ! temporal mask of evaluated arr
    real(dp),    dimension(:,:), allocatable, intent(out) :: opt_h       ! optimized kernel width
    real(dp),    dimension(:,:), allocatable, intent(out) :: lats, lons  ! latitude and longitude fields of input
    real(sp),                                 intent(out) :: SMI_thld    ! SMI threshold for clustering
    character(len=256),                       intent(out) :: outpath     ! ouutput path for results

    ! local Variables
    integer(i4)                                     :: ii, iip
    integer(i4)                                     :: dd, mm, yy
    
    integer(i4)                                     :: datatype       ! datatype of attribute
    integer(i4)                                     :: nCells         ! number of effective cells
    integer(i4)                                     :: number_lag_days! number of lag days for daily files (0,7,31)
    integer(i4)                                     :: jDayStart
    integer(i4)                                     :: jDayEnd
    integer(i4)                                     :: jDay
    integer(i4)                                     :: nTimeSteps     ! number of time steps excluding leap days 

    
    ! directories, filenames, attributes
    character(256)                                  :: AttValues      ! netcdf attribute values
    character(256)                                  :: maskfName
    character(256)                                  :: opt_h_file
    character(256)                                  :: SM_eval_file
    character(256)                                  :: basinfName
    
    ! variable names in netcdf input files
    character(256)                                  :: soilmoist_file
    character(256)                                  :: mask_vname
    character(256)                                  :: SM_vname
    character(256)                                  :: SM_eval_vname
    character(256)                                  :: type_var_eval
    character(256)                                  :: type_time_eval
    character(256)                                  :: smi_file_clustering
    character(256)                                  :: basin_vname
    character(256)                                  :: opt_h_vname
    real(sp)                                        :: nodata_value ! local data nodata value in nc for mask creation

    real(sp),       dimension(:,:),     allocatable :: dummy_D2_sp
    real(sp),       dimension(:,:,:),   allocatable :: dummy_D3_sp
    real(dp),       dimension(:,:,:),   allocatable :: dummy_D3_dp

    ! read main config
    namelist/mainconfig/basin_flag, basinfName, basin_vname, maskfName, mask_vname, soilmoist_file, SM_vname,  &
         outpath, &
         nCalendarStepsYear, number_lag_days, &
         eval_SMI, silverman_h, read_opt_h, opt_h_vname, opt_h_file,     &
         SM_eval_file, SM_eval_vname, type_var_eval, type_time_eval, SMI_thld,    &
         do_cluster, cluster_ext_smi, smi_file_clustering, cellsize, thCellClus, nCellInter, &
         do_sad, deltaArea


    ! dummy init to avoid gnu compiler complaint
    outpath    = '-999'; cellsize=-999.0_sp; thCellClus=-999; nCellInter=-999; deltaArea=-999
    
    do_cluster = .FALSE.
    do_sad = .FALSE.
    SMI_thld   = 0.0_dp
    silverman_h = .FALSE.
    
    ! read namelist
    open (unit=10, file='main.dat', status='old')
      read(10, nml=mainconfig)
    close (10)
    print*, 'main.dat read ...ok'

    ! consistency check
    if ( ( trim( type_var_eval ) .ne. 'dp' ) .and. &
         ( trim( type_var_eval ) .ne. 'sp' ) ) &
         stop '***ERROR: mo_read: var_type_eval must be dp or sp'

    if ( .not.(( nCalendarStepsYear == 12) .or. ( nCalendarStepsYear == 365) ) ) &
         stop '***ERROR: number of time steps per year different from 12 or 365'

    if ( ( number_lag_days .gt. 0 ) .and. ( nCalendarStepsYear == 12)  ) &
         stop '***ERROR: number_lag_days = 0 for montly SM values'

    if ( do_sad .and. (.not. do_cluster) ) &
         stop '***ERROR: flag do_cluster must be set to .TRUE. to perform SAD analisys'
    
    ! read main mask
    call Get_ncVar( maskfName, trim(mask_vname), dummy_D2_sp )
    allocate( mask( size(dummy_D2_sp,1), size(dummy_D2_sp,2) ) )

    ! determine no data value
    call Get_NcVarAtt(maskfName, trim(mask_vname), 'missing_value', AttValues, dtype=datatype)
    ! convert to number 
    read(AttValues, *) nodata_value
    ! create mask
    mask = merge( .true., .false., notequal(dummy_D2_sp, nodata_value ) )
    deallocate( dummy_D2_sp )
    ! consistency check
    nCells   = count( mask )
    if ( nCells .eq. 0 ) then
       print *, '***ERROR: no cell selected in mask'
       stop
    end if 
    
    print*, 'mask read ...ok'

    ! read basin mask
    if ( basin_flag .AND. (.NOT. cluster_ext_smi)) then
       call Get_ncVar( basinfName, trim(basin_vname), Basin_Id )
       ! determine no data value
       call Get_NcVarAtt( maskfName, trim(mask_vname), 'missing_value', AttValues, dtype=datatype)
       ! convert to number
       read( AttValues, *) nodata_value 
       ! consistency check
       if ( ( size( Basin_Id, 1) .ne. size( mask, 1 ) ) .or. &
            ( size( Basin_Id, 2) .ne. size( mask, 2 ) ) ) then
          print *, '***ERROR: size mismatch between basin field and given mask file'
          stop
       end if
       ! intersect mask and basin_id
       mask     = merge( .true., .false., mask .and. ( Basin_Id .ne. int( nodata_value, i4) ) )
       Basin_Id = merge( Basin_Id, int( nodata_value, i4), mask )
       !Basin_Id = pack(  Basin_ID, mask)
       print*, 'basin ID read ...ok'
    end if

    ! read SM field
    if ( .NOT. cluster_ext_smi )  then
       call Get_ncVar( soilmoist_file, trim( SM_vname ), dummy_D3_dp )
       ! consistency check
       if ( ( size( dummy_D3_dp, 1) .ne. size( mask, 1 ) ) .or. &
            ( size( dummy_D3_dp, 2) .ne. size( mask, 2 ) ) ) then
          print *, '***ERROR: size mismatch between SM field and given mask file'
          stop
       end if

       ! find number of leap years in  SM data set                                    (moved from below)
       ! get timepoints in days and masks for climatologies at calendar day or month 
       call get_time( soilmoist_file,  size( dummy_D3_dp, 3 ), trim(type_time_eval), &
            yStart, mStart, dStart, &
            yEnd, timepoints) !, tmask_est)
      
       ! read lats and lon from file
       call Get_ncVar( soilmoist_file, 'lat', lats )
       call Get_ncVar( soilmoist_file, 'lon', lons )
       
       
       if ( nCalendarStepsYear .eq. YearMonths ) then
          allocate( SM_est( nCells, size( dummy_D3_dp, 3 ) ) )
          ! no lag average, use monthly data as read
          do ii = 1, size( dummy_D3_dp, 3 )
             SM_est(:,ii) = pack( real(dummy_D3_dp(:,:,ii),sp), mask )
          end do
       else
          ! make moving-lag average (365 days)
          ! count days removing leap days
          jDayStart =  NDAYS(dStart, mStart, yStart)
          jDayEnd = NDAYS(31, 12, yEnd)                             ! >>> to be improved
          nTimeSteps = 0
          do ii = 1, size( dummy_D3_dp, 3 )
             call NDYIN( (jDayStart+ii-1), dd, mm, yy)
             if ( (mm .eq. 2) .and. ( dd .eq. 29) ) cycle
             nTimeSteps = nTimeSteps + 1
          end do
          
          if ( abs( jDayEnd - jDayStart + 1 -  size( dummy_D3_dp, 3 ) ) .gt. 0    ) then
             print *, '***ERROR: SM field does not cover complete years'
             stop
          end if
          
          allocate( SM_est( nCells, nTimeSteps ) )
          ! exclude leap days
          ! running average from the end to the begin of the vector 
          jDay = nTimeSteps   
          do ii = 1, size( dummy_D3_dp, 3 )
             iip = size( dummy_D3_dp, 3 ) - ii + 1
             call NDYIN( (jDayEnd-ii+1), dd, mm, yy)
             if ( ( mm .eq. 2 ) .and. ( dd .eq. 29 ) ) cycle
             if ( jDay .lt. number_lag_days ) then
                SM_est(:,jDay) = pack( real(dummy_D3_dp(:,:,iip ),sp), mask )
             else
                SM_est(:,jDay) = pack( real( sum(dummy_D3_dp(:,:,iip-number_lag_days+1:iip), DIM=3) &
                     / real(number_lag_days,dp) ,sp), mask )
             end if
             jDay = jDay - 1
          end do
       end if
  
       deallocate( dummy_D3_dp )

       !>>>> moved up
       ! ! get timepoints in days and mask of months
       ! call get_time(soilmoist_file, size(SM_est, dim=2), trim(type_time_eval), &
       !      yStart, mStart, dStart, &
       !      yEnd, timepoints, tmask_est)
      
       ! ! read lats and lon from file
       ! call Get_ncVar( soilmoist_file, 'lat', lats )
       ! call Get_ncVar( soilmoist_file, 'lon', lons )
    end if

! *************************************************************************
! >>> EVAL - NOT YET DONE !!
! *************************************************************************
    
    ! check whether second SM field and optimized h should be read
    if ( eval_SMI  .AND. (.NOT. cluster_ext_smi)) then
       ! read second SM field that uses CDF of the first one
       if ( trim( type_var_eval ) .eq. 'dp' ) then
          ! var_type is double precision
          call Get_ncVar( trim( SM_eval_file ), trim( SM_eval_vname ), dummy_D3_dp )
          allocate( SM_eval( nCells, size( dummy_D3_dp, 3 ) ) )
          do ii = 1, size( dummy_D3_dp, 3 )
             SM_eval(:, ii) = pack( real(dummy_D3_dp(:,:,ii),sp), mask)
          end do
          deallocate( dummy_D3_dp )
       end if
       if ( trim( type_var_eval ) .eq. 'sp' ) then
          ! var_type is single precision
          call Get_ncVar( trim( SM_eval_file ), trim( SM_eval_vname ), dummy_D3_sp )
          allocate( SM_eval( nCells, size( dummy_D3_sp, 3 ) ) )
          do ii = 1, size( dummy_D3_sp, 3 )
             SM_eval(:, ii) = pack( dummy_D3_sp(:,:,ii), mask)
          end do
          deallocate( dummy_D3_sp )
       end if

       ! get timepoints in days and mask of months
       if ( allocated( timepoints ) ) deallocate( timepoints )
       call get_time( SM_eval_file, size(SM_eval, dim=2), trim(type_time_eval), &
            yStart, mStart, dStart, &
            yEnd, timepoints ) !, tmask_eval)
       
       ! if ( all( count( tmask_eval, dim = 1 ) .eq. 0_i4 ) ) &
       !      stop '***ERROR no data in eval given, check time axis'
       ! print*, 'read soil moisture field for evaluation... ok'
    end if
    
    ! initialize opt_h
    ! subroutine read_kernerl_width_h(opt_h_file, opt_h_vname, opt_h, mask, )

    allocate ( opt_h( ncells, nCalendarStepsYear ) )                         ! OLD: YearMonths 
    opt_h = nodata_dp
    !
    if ( read_opt_h  .AND. (.NOT. cluster_ext_smi)) then
       ! read optimized kernel width from file
       call Get_ncVar( trim( opt_h_file ), trim( opt_h_vname ), dummy_D3_dp )
       do ii = 1, size( dummy_D3_dp, 3 )
          opt_h( :, ii ) = pack( real( dummy_D3_dp(:,:,ii),sp ), mask )
       end do
       deallocate( dummy_D3_dp )
       if ( any( equal( opt_h, nodata_dp ) ) ) &
            stop '***ERROR kernel width contains nodata values'

       print *, 'read kernel width from file... ok'
    end if    

    if (do_cluster .AND. cluster_ext_smi) then
       call Get_ncVar( smi_file_clustering, trim('SMI'), dummy_D3_sp )
       ! consistency check
       if ( ( size( dummy_D3_sp, 1) .ne. size( mask, 1 ) ) .or. &
            ( size( dummy_D3_sp, 2) .ne. size( mask, 2 ) ) ) then
          print *, '***ERROR: size mismatch between SMI field and given mask file'
          stop
       end if
       allocate( SMI( nCells, size( dummy_D3_sp, 3 ) ) )
       do ii = 1, size( dummy_D3_sp, 3 )
          SMI(:,ii) = pack( real(dummy_D3_sp(:,:,ii),sp), mask )
       end do
       deallocate( dummy_D3_sp )

       ! get timepoints in days and mask of months
       call get_time( smi_file_clustering, size(SMI, dim=2), trim(type_time_eval), &
            yStart, mStart, dStart, &
            yEnd, timepoints)  !, tmask_est

       ! if ( any( count( tmask_est, dim = 1 ) .eq. 0_i4 ) ) &
       !      stop '***ERROR no data for estimation given for all calendar months, check time axis'

       ! read lats and lon from file
       call Get_ncVar( smi_file_clustering, 'lat', lats )
       call Get_ncVar( smi_file_clustering, 'lon', lons )
    end if

! <<<---------------------------------------------------------------------------------------------
  end subroutine ReadDataMain


  !
  !     PORPOSE
  !         Convert input time into months & check timesteps & determine mask for months
  
  !     CALLING SEQUENCE
  !         call get_time(fName, vname, timevector, maskvector)
  
  !     RESTRICTIONS
  !         None
  
  !     LITERATURE
  !         None
  
  !     HISTORY
  !         Written,  Matthias Zink, Oct 2012

  subroutine get_time(fName, nTimeSteps, dtype, yStart, mStart, dStart, yEnd, timepoints) !, mask)
    !
    use mo_julian,       only: date2dec, dec2date
    use mo_message,      only: message
    use mo_NcRead,       only: Get_NcVar, Get_NcVarAtt
    use mo_string_utils, only: DIVIDE_STRING

    use mo_smi_constants, only: nodata_i4, YearMonths, DayHours

    implicit none

    character(len=*),                         intent(in)  :: fName      ! name of NetCDF file
    integer(i4),                              intent(in)  :: nTimeSteps ! size of the time dimension
    character(len=*),                         intent(in)  :: dtype      ! data type of var time
    integer(i4),                              intent(out) :: yStart     ! start year  of the dataser
    integer(i4),                              intent(out) :: mStart     ! start month of the dataser
    integer(i4),                              intent(out) :: dStart     ! start day   of the dataser
    integer(i4),                              intent(out) :: yEnd       ! end year of the dataset
    integer(i4), dimension(:),   allocatable, intent(out) :: timepoints ! timestep in days or months
!    logical    , dimension(:,:), allocatable, intent(out) :: mask       ! masking days/months in timespace
    
    integer(i4)                               :: i                      ! loop variable
    integer(i4)                               :: yRef, dRef, mRef       ! reference time of NetCDF (unit attribute of
    integer(i4)                               :: month, year, d         ! 
    integer(i4)                               :: datatype               ! datatype of attribute
    !
    integer(i4),   dimension(:), allocatable  :: timesteps              ! time variable of NetCDF in input units
    real(sp),      dimension(:), allocatable  :: timesteps_sp           ! time variable of NetCDF in input units
    real(dp),      dimension(:), allocatable  :: timesteps_dp           ! time variable of NetCDF in input units
    !
    character(256)                            :: AttValues              ! netcdf attribute values
    character(256), dimension(:), allocatable :: strArr                 ! dummy for netcdf attribute handling
    character(256), dimension(:), allocatable :: date                   ! dummy for netcdf attribute handling
    real(dp)                                  :: ref_jday               ! refernce date of dateset as julian day 
    !
    ! get unit attribute of variable 'time'
    call Get_NcVarAtt(fName, 'time', 'units', AttValues, dtype=datatype)
    ! AttValues looks like "<unit> since YYYY-MM-DD HH:MM:SS"
    call DIVIDE_STRING(trim(AttValues), ' ', strArr)
    !
    call DIVIDE_STRING(trim(strArr(3)), '-', date) 
    read(date(1),*) yRef
    read(date(2),*) mRef
    read(date(3),*) dRef

    ! allocate and initialize
    allocate(timesteps  ( nTimeSteps ))
    allocate(timepoints ( nTimeSteps ))
!    allocate(mask ( nTimeSteps , YearMonths))
    timesteps  = nodata_i4
    timepoints = 0
!    mask       = .FALSE.

    !!
    select case ( trim(dtype) )
    case ( 'i4' ) 
       call Get_NcVar(fName, 'time', timesteps)
    case ( 'sp' )
       call Get_NcVar(fName, 'time', timesteps_sp )
       timesteps = int( timesteps_sp, i4 )
       deallocate( timesteps_sp )
    case ( 'dp' )
       call Get_NcVar(fName, 'time', timesteps_dp )
       timesteps = int( timesteps_dp, i4 )
       deallocate( timesteps_dp )
    case default
       print*, "***ERROR: Not a valid FORMAT for the time variable given in the namelist"
       stop
    end select
    
    ! strArr(1) is <unit>
    ref_jday = date2dec(dd=dRef, mm=mRef, yy=yRef)
    do i = 1, nTimeSteps
       select case (strArr(1))
       case('hours')
          call dec2date(real(timesteps(i), dp)/real(DayHours, dp) + ref_jday, dd=d, mm=month, yy=year)
          timepoints(i) = timepoints(i) + nint( real(timesteps(i) - timesteps(1), dp) / real(DayHours, dp) , i4)
       case('days')
          call dec2date(real(timesteps(i), dp) + ref_jday, dd=d, mm=month, yy=year) 
          timepoints(i) = timepoints(i) + timesteps(i) - timesteps(1)
       case('months')
          d       = dRef ! set to dref as default
          month   = mod( (timesteps(i) + mRef ), YearMonths )
          year    = yRef + floor(real( timesteps(i) + mRef, dp) / real(YearMonths, dp) )
          ! correct month and year for december
          if (month == 0 ) then
             month = 12
             year  = year - 1
          end if
          !
          if (i .EQ. 1) then
             timepoints(i) = 0
          else
             timepoints(i) = timepoints(i) +  nint(date2dec(dd=15 ,mm=month, yy=year) - date2dec(dd=15 ,mm=mStart, yy=yStart), i4)
          end if
       case DEFAULT
          call message('***ERROR: Time slicing in ' , trim(strArr(1)), ' not implemented!')
          stop
       end select

       if (i .EQ. 1) then
          yStart = year  ! save start year and month output timestamp
          mStart = month
          dStart = d
       end if
       if (i .EQ. nTimeSteps) yEnd   = year  ! save end year

!       ! save month in tim mask
!       mask(i,month) = .TRUE.
    end do

  end subroutine get_time

end module mo_read
