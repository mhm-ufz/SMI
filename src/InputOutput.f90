!**********************************************************************
!  MODULE InputOutput                                                 *
!  PURPOSE     Input / output subroutines                             *
!  CREATED     Luis Samaniego, 09.02.2011                             *
!                                                                     *
!**********************************************************************
module InputOutput

  use mo_kind,   only                              : i4, sp, dp

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

  ! ##################################################################
  ! subroutine for writting the SMI to nc file
  ! author: Stephan Thober
  ! created: 5.7.2014
  ! ##################################################################
  subroutine WriteSMI( outpath, SMI, mask, yStart, mStart, dStart,  yEnd, &
    timepoints, nCalendarStepsYear, lats, lons, hh ) 

    use mo_kind,          only: i4, sp
    use mo_string_utils,  only: num2str
    use mo_smi_constants, only: nodata_sp, nodata_dp
    use mo_julian,        only: NDAYS, NDYIN
    use mo_netcdf,        only: NcDataset, NcVariable, NcDimension

    implicit none

    ! input variables
    character(len=*),                              intent(in) :: outpath     ! ouutput path for results
 
    logical,        dimension(:,:),                intent(in) :: mask
    integer(i4),                                   intent(in) :: yStart
    integer(i4),                                   intent(in) :: mStart
    integer(i4),                                   intent(in) :: dStart
    integer(i4),                                   intent(in) :: yEnd
    real(sp),       dimension(:,:),                intent(in) :: SMI 
    integer(i4),    dimension(:),   allocatable,   intent(in) :: timepoints
    integer(i4),                                   intent(in) :: nCalendarStepsYear
    real(dp),       dimension(:,:), allocatable,   intent(in) :: lats, lons   ! latitude and longitude fields of input
    real(dp),       dimension(:,:), optional,      intent(in) :: hh

    ! local Variables
    type(NcDataset)                                           :: nc_out
    type(NcVariable)                                          :: nc_var
    type(NcDimension)                                         :: nc_row, nc_col, nc_tim, nc_cal
    character(256)                                            :: Fname
    integer(i4)                                               :: dd          ! day
    integer(i4)                                               :: mm          ! month
    integer(i4)                                               :: yy          ! year

    real(sp),       dimension(:,:,:), allocatable             :: dummy_D3_sp ! field real unpacked 
    real(dp),       dimension(:,:,:), allocatable             :: dummy_D3_dp ! field real unpacked 
    integer(i4)                                               :: nrows
    integer(i4)                                               :: ncols
    integer(i4)                                               :: jDayStart
    integer(i4)                                               :: jDayEnd
    integer(i4)                                               :: nTotalTimeSteps     ! number of time steps including leap days 
    integer(i4)                                               :: tt                  ! counter without leap days
    integer(i4)                                               :: jDay                ! jDay counting from jDayStart
    
    ! initialize dimension
    nrows = size(mask, 1)
    ncols = size(mask, 2)

    ! fname
    fName = trim(outpath) //'SMI.nc'
    nc_out = NcDataset(fName, 'w')
    nc_row = nc_out%setDimension("nrows", nrows)
    nc_col = nc_out%setDimension("ncols", ncols)
    nc_tim = nc_out%setDimension("time", -1)
    
    ! unpack estimated SMIp
    if (nCalendarStepsYear .eq. 12) then
      allocate( dummy_d3_sp( nrows, ncols, size( SMI, 2) ) )
      nTotalTimeSteps = size( SMI, 2)

      do mm = 1,  nTotalTimeSteps
        dummy_d3_sp( :, :, mm) = unpack ( SMI(:,mm), mask, nodata_sp )
      end do

    else
      ! store SMI including leap days, assumption => SMI (29.02.yyyy) ~ SMI (28.02.yyyy)
      jDayStart =  NDAYS(dStart, mStart, yStart)
      jDayEnd = NDAYS(31, 12, yEnd)                            
      nTotalTimeSteps = jDayEnd - jDayStart + 1
      allocate( dummy_d3_sp( nrows, ncols, nTotalTimeSteps  ) )

      tt = 0
      do jDay = 1,  nTotalTimeSteps
        call NDYIN( (jDayStart+jDay-1), dd, mm, yy)
        if ( ( mm .eq. 2 ) .and. ( dd .eq. 29 ) ) then
          dummy_d3_sp( :, :, jDay) = unpack ( SMI(:,tt), mask, nodata_sp )
        else
          tt = tt + 1
          dummy_d3_sp( :, :, jDay) = unpack ( SMI(:,tt), mask, nodata_sp )
        end if
      end do

    end if
    
    ! save SMI (same sequence as in 1)
    nc_var = nc_out%setVariable('SMI', "f32", &
         (/ nc_row, nc_col, nc_tim /))
    call nc_var%setData(dummy_d3_sp)
    call nc_var%setAttribute('long_name', 'soil moisture index')
    call nc_var%setAttribute('missing_value', nodata_sp)
    call nc_var%setAttribute('units', '-')
    deallocate( dummy_d3_sp )

    ! write out kernel width if it has been optimised
    if ( present( hh ) ) then
       allocate( dummy_D3_dp( nrows, ncols, size( hh, 2 ) ) )
       do mm = 1, size( hh, 2 )
          dummy_D3_dp(:, :, mm) = unpack( hh(:,mm), mask, real( nodata_sp, dp ) ) 
       end do
       nc_cal = nc_out%setDimension("calendar_steps", size(dummy_d3_dp,3))
       nc_var = nc_out%setVariable('kernel_width', "f64", &
            (/ nc_row, nc_col, nc_cal /))
       call nc_var%setData(dummy_d3_dp)
       call nc_var%setAttribute('long_name', 'optimised kernel width')
       call nc_var%setAttribute('missing_value', nodata_dp)
       call nc_var%setAttribute('units', '-')
       deallocate( dummy_D3_dp ) 
    end if

    ! add lat and lon
    if (allocated(lats)) then
      nc_var = nc_out%setVariable('lat', "f64", (/ nc_row, nc_col /))
      call nc_var%setData(lats)
      call nc_var%setAttribute('long_name', 'latitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_north')
    end if
    if (allocated(lons)) then
      nc_var = nc_out%setVariable('lon', "f64", (/ nc_row, nc_col /))
      call nc_var%setData(lons)
      call nc_var%setAttribute('long_name', 'longitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_east')
    end if

    ! add time
    nc_var = nc_out%setVariable('time', "i32", (/ nc_tim /))
    call nc_var%setData(timepoints)
    call nc_var%setAttribute('long_name', 'time')
    call nc_var%setAttribute('units', 'days since '                              // &
                                      trim(num2str(yStart,   '(i4)')) // '-'     // &
                                      trim(num2str(mStart, '(i2.2)')) // '-'     // &
                                      trim(num2str(dStart, '(i2.2)')) // ' 00:00:00' )
    ! close file
    call nc_out%close()
    
  end subroutine WriteSMI

  !*************************************************************************
  !    PURPOSE    WRITE netCDF files
  !    FORMAT     netCDF
  !               http://www.unidata.ucar.edu/software/netcdf/  
  !
  !    AUTHOR:    Luis E. Samaniego-Eguiguren, UFZ
  !    UPDATES
  !               Created        Sa   16.02.2011   
  !               Last Update    Sa   16.02.2011  
  !**************************************************************************
  subroutine WriteNetCDF(outpath, wFlag, yStart, dStart, mStart, timepoints, lats, lons, &
        SMIc, SM_invert, duration) 
    !
    use mo_kind,          only: i4
    use mo_string_utils,  only: num2str
    use mo_smi_constants, only: nodata_i4 , nodata_dp, nodata_sp
    use mo_netcdf,        only: NcDataset, NcVariable, NcDimension

    implicit none
    !
    ! input variables
    character(len=*),                            intent(in) :: outpath     ! ouutput path for results
    integer(i4),                                 intent(in) :: wFlag
    integer(i4),                                 intent(in) :: yStart
    integer(i4),                                 intent(in) :: mStart
    integer(i4),                                 intent(in) :: dStart
    integer(i4),    dimension(:),   allocatable, intent(in) :: timepoints
    real(dp),       dimension(:,:), allocatable, intent(in) :: lats, lons   ! latitude and longitude fields of input
    integer(i4),    dimension(:,:,:), optional,  intent(in) :: SMIc         ! Drought indicator
    real(sp),       dimension(:,:,:), optional,  intent(in) :: SM_invert
    integer(i4),                      optional,  intent(in) :: duration     ! optional, duration

    ! local Variables
    type(NcDataset)              :: nc_out
    type(NcVariable)             :: nc_var
    type(NcDimension)            :: nc_row, nc_col, nc_tim
    character(256)               :: Fname

    select case (wFlag)

    case (2)
       ! fname
       fName  = trim(outpath)//'SM_invert.nc'
       nc_out = NcDataset(fname, 'w')

       nc_row = nc_out%setDimension("nrows", size(SM_invert, 1))
       nc_col = nc_out%setDimension("ncols", size(SM_invert, 2))
       nc_tim = nc_out%setDimension("time", -1)

       nc_var = nc_out%setVariable('SM_Lall', "f32", &
            (/ nc_row, nc_col, nc_tim /))
       call nc_var%setData(SM_invert)
       call nc_var%setAttribute('long_name', 'SM according to inverse of SMI')
       call nc_var%setAttribute('missing_value', nodata_sp)
       call nc_var%setAttribute('units', 'mm')
       
       
    case (3)
       ! fname
       fName  = trim(outpath)//'SMIc.nc'
       nc_out = NcDataset(fname, 'w')

       nc_row = nc_out%setDimension("nrows", size(SMIc, 1))
       nc_col = nc_out%setDimension("ncols", size(SMIc, 2))
       nc_tim = nc_out%setDimension("time", -1)

       nc_var = nc_out%setVariable('mSMIc', "i32", &
            (/ nc_row, nc_col, nc_tim /))
       call nc_var%setData(SMIc)
       call nc_var%setAttribute('long_name', 'monthly SMI indicator SMI < th')
       call nc_var%setAttribute('missing_value', nodata_i4)
       call nc_var%setAttribute('units', '-')

    case (4)
       ! fname
       fName = trim(outpath)//'DCluster.nc'
       nc_out = NcDataset(fname, 'w')

       nc_row = nc_out%setDimension("nrows", size(idCluster, 1))
       nc_col = nc_out%setDimension("ncols", size(idCluster, 2))
       nc_tim = nc_out%setDimension("time", -1)

       nc_var = nc_out%setVariable('mDC', "i32", &
            (/ nc_row, nc_col, nc_tim /))
       call nc_var%setData(idCluster)
       call nc_var%setAttribute('long_name', 'consolidated cluster evolution')
       call nc_var%setAttribute('missing_value', nodata_i4)
       call nc_var%setAttribute('units', '-')
         
    case (5)
       ! fname
       write (fName, '(i2.2)') duration
       fName = trim(outpath)//'sev_'//trim(fName)//'.nc'
       nc_out = NcDataset(fname, 'w')

       nc_row = nc_out%setDimension("nrows", size(severity, 1))
       nc_col = nc_out%setDimension("ncols", size(severity, 2))
       nc_tim = nc_out%setDimension("time", -1)

       nc_var = nc_out%setVariable('Severity', "f64", &
            (/ nc_row, nc_col, nc_tim /))
       call nc_var%setData(idCluster)
       call nc_var%setAttribute('long_name', 'd-month severity')
       call nc_var%setAttribute('missing_value', nodata_dp)
       call nc_var%setAttribute('units', '-')

    end select

    ! add lat and lon
    if (allocated(lats)) then
      nc_var = nc_out%setVariable('lat', "f64", (/ nc_row, nc_col /))
      call nc_var%setData(lats)
      call nc_var%setAttribute('long_name', 'latitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_north')
    end if
    if (allocated(lons)) then
      nc_var = nc_out%setVariable('lon', "f64", (/ nc_row, nc_col /))
      call nc_var%setData(lons)
      call nc_var%setAttribute('long_name', 'longitude')
      call nc_var%setAttribute('missing_value', nodata_dp)
      call nc_var%setAttribute('units', 'degrees_east')
    end if
    
    ! add time
    if ((wflag .eq. 2) .or. &
        (wflag .eq. 3) .or. &
        (wflag .eq. 4)) then
       nc_var = nc_out%setVariable('time', "i32", (/ nc_tim /))
       call nc_var%setData(timepoints)
       call nc_var%setAttribute('long_name', 'time')
       call nc_var%setAttribute('units', 'days since '                              // &
            trim(num2str(yStart,   '(i4)')) // '-'     // &
            trim(num2str(mStart, '(i2.2)')) // '-'     // &
            trim(num2str(dStart, '(i2.2)')) // ' 00:00:00' )
    end if
    ! close file
    call nc_out%close()
  end subroutine WriteNetCDF

!**********************************************************************
!    PURPOSE    WRITE Results of the cluster analysis
!
!    AUTHOR:    Luis E. Samaniego-Eguiguren, UFZ
!    UPDATES
!               Created        Sa   17.03.2011
!               Last Update    Sa
!**********************************************************************
subroutine WriteResultsCluster(SMIc, outpath, wFlag, yStart, yEnd, nMonths, nCells, &
     deltaArea, cellsize, d)

  use mo_kind,          only : i4, dp
  use mo_smi_constants, only : YearMonths

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
        do m = 1, YearMonths
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
        do m = 1, YearMonths
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


!**********************************************************************
!    PURPOSE    WRITE Results of the cluster SMI basins
!
!    AUTHOR:    Luis E. Samaniego-Eguiguren, UFZ
!    UPDATES
!               Created        Sa   24.05.2011
!               Last Update    Sa
!**********************************************************************
subroutine WriteResultsBasins( outpath, SMI, mask, yStart, yEnd, nMonths, Basin_Id )

  use mo_kind,          only : i4
  use mo_smi_constants, only : nodata_dp, YearMonths

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
     do j = 1, YearMonths
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

end module InputOutput
