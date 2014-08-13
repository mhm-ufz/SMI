!**********************************************************************
!  MODULE InputOutput                                                 *
!  PURPOSE     Input / output subroutines                             *
!  CREATED     Luis Samaniego, 09.02.2011                             *
!                                                                     *
!**********************************************************************
module InputOutput

  use mo_kind,   only                              : i4, sp, dp
  use mo_ncread, only                              : Get_NcVar, Get_NcDim, Get_NcVarAtt

  implicit none

  character(256)                                  :: optPDFparfName    ! opt parameter PDF file 
  ! fields
  real(sp),    dimension(:,:), allocatable        :: h_opt             ! optimized kernel width h
  real(dp), dimension(:,:), allocatable           :: SMIp              ! SMI field packed
  real(dp), dimension(:,:,:), allocatable         :: SMI               ! SMI field unpacked
  integer(i4), dimension(:,:,:), allocatable      :: SMIc              ! SMI indicator
  !
  ! clusters
  integer(i4), dimension(:,:,:), allocatable      :: idCluster
  integer(i4), dimension(:,:), allocatable        :: cellCoor
  real(dp)                                        :: cellArea
  integer(i4)                                     :: nInterArea
  integer(i4)                                     :: nEvents           ! total number of drough events to estimate SAD
  integer(i4), dimension(:), allocatable          :: shortCnoList      ! consolidated cluster no. list
  integer(i4), dimension(:,:), allocatable        :: eventId           ! event identification number, cluster, month of occurence
  integer(i4)                                     :: nClusters
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

  !
  ! main parameters 
  integer(i4)                                    :: DAT_flag     ! daily of monthly inputs
  ! for mHM (4 x 4) km2
!!$  integer(i4), parameter                         :: thCellClus = 40        ! treshold  for cluster formation in space ~ 640 km2
!!$  integer(i4), parameter                         :: nCellInter = 400       ! number cells for joining clusters in time ~ 6400 km2
!!$  integer(i4), parameter                         :: deltaArea  = 20        ! number of cells per area interval
  ! for COSMO ~(7 x 7) km2
  integer(i4), parameter                         :: thCellClus = 13        ! treshold  for cluster formation in space ~ 640 km2
  integer(i4), parameter                         :: nCellInter = 130       ! number cells for joining clusters in time ~ 6400 km2
  integer(i4), parameter                         :: deltaArea  = 7         ! number of cells per area interval
  !
  integer(i4), parameter                         :: nDurations = 4         ! number of durations
  integer(i4), dimension(nDurations), parameter  :: durList = (/3,6,9,12/) !(/3,6,9,12/)   ! list of durations to evaluate
  integer(i4), parameter                         :: nQProp = 3             ! number of SAD percetiles for a given duration
  ! real(dp), dimension(nQProp), parameter         :: QProp = (/0.90_dp, 0.96_dp, 0.98_dp /)
  real(dp), dimension(nQProp), parameter         :: QProp = (/90._dp, 96._dp, 98._dp /)
                                                                           ! percentiles corresponding
                                                                           ! to return periods of 10,25,50 years
  integer(i4)                                    :: nLargerEvents = 600    ! n. largest drought events
  !
  ! Basin summary
  integer(i4), parameter                         :: nBasins = 6
   ! netCDF
  integer(i4), dimension(:), allocatable         :: eMask                   ! error field

contains

  ! ##################################################################
  ! subroutine for writting the SMI to nc file
  ! author: Stephan Thober
  ! created: 5.7.2014
  ! ##################################################################
  subroutine WriteSMI( outpath, SMI, mask, yStart, mStart, dStart, times, lats, lons, hh ) 

    use mo_kind,          only: i4, sp
    use mo_string_utils,  only: num2str
    use mo_ncwrite,       only: var2nc
    use mo_smi_constants, only: nodata_sp, YearMonths

    implicit none

    ! input variables
    character(len=*),                              intent(in) :: outpath     ! ouutput path for results
 
    logical,        dimension(:,:),                intent(in) :: mask
    integer(i4),                                   intent(in) :: yStart
    integer(i4),                                   intent(in) :: mStart
    integer(i4),                                   intent(in) :: dStart
    real(sp),       dimension(:,:),                intent(in) :: SMI 
    integer(i4),    dimension(:),   allocatable,   intent(in) :: times
    real(dp),       dimension(:,:),                intent(in) :: lats, lons   ! latitude and longitude fields of input
    real(dp),       dimension(:,:), optional,      intent(in) :: hh

    ! local Variables
    character(256)                                            :: Fname
    integer(i4)                                               :: mm          ! month
    real(sp),       dimension(:,:,:), allocatable             :: dummy_D3_sp ! field real unpacked 
    real(dp),       dimension(:,:,:), allocatable             :: dummy_D3_dp ! field real unpacked 
    ! dimension names
    character(256), dimension(3)                              :: dnames
    character(256), dimension(3)                              :: dims_hopt

    ! initialize dimension names
    dnames(1) = 'nrows'
    dnames(2) = 'ncols'
    dnames(3) = 'time'    
    dims_hopt(1) = 'nrows'
    dims_hopt(2) = 'ncols'
    dims_hopt(3) = 'months'

    ! fname
    fName = trim(outpath)//'mSMI.nc'

    ! unpack estimated SMIp
    allocate( dummy_d3_sp( size( mask, 1), size( mask, 2), size( SMI, 2) ) )
    do mm = 1, size( SMI, 2 )
       dummy_d3_sp( :, :, mm) = unpack ( SMI(:,mm), mask, nodata_sp )
    end do
    ! save SMI (same sequence as in 1)
    call var2nc( fName, dummy_d3_sp, dnames, v_name = 'SMI', &
         long_name = 'monthly soil moisture index', units = '-', &
         missing_value = nodata_sp, create = .true. ) ! &
    ! scale_factor = 1., coordinates = 'lon lat' )
    deallocate( dummy_d3_sp )

    ! write out kernel width if it has been optimised
    if ( present( hh ) ) then
       allocate( dummy_D3_dp( size( mask, 1 ), size( mask, 2), size( hh, 2 ) ) )
       do mm = 1, YearMonths
          dummy_D3_dp( :, :, mm ) = unpack( hh(:,mm), mask, real( nodata_sp, dp ) ) 
       end do
       call var2nc( fname, dummy_D3_dp, dims_hopt, v_name = 'kernel_width', &
            long_name = 'optimised kernel width', units = '-', &
            missing_value = real( nodata_sp, dp ) )
       deallocate( dummy_D3_dp ) 
    end if

    ! add lat and lon
    call var2nc( fname, lats, dnames(1:2), v_name = 'lat', &
         long_name = 'longitude', units = 'degrees_east' )
    call var2nc( fname, lons, dnames(1:2), v_name = 'lon', &
         long_name = 'latitude', units = 'degrees_north' )

    ! add time
    call var2nc( fname, times, dnames(3:3), v_name='time', &
         long_name = 'time', units = 'days since '                               // &
                                      trim(num2str(yStart,   '(i4)')) // '-'     // &
                                      trim(num2str(mStart, '(i2.2)')) // '-'     // &
                                      trim(num2str(dStart, '(i2.2)')) // ' 00:00:00' )
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
  subroutine WriteNetCDF(outpath, wFlag, opt_h, SM_est, mask, yStart, lats, lons, d, SMI_eval) 
    !
    use mo_kind,          only: i4, sp
    use mo_string_utils,  only: num2str
    use mo_ncwrite,       only: var2nc
    use mo_smi_constants, only: nodata_i4 , nodata_dp 

    implicit none
    !
    ! input variables
    character(len=*),                         intent(in) :: outpath     ! ouutput path for results 
    integer(i4),                              intent(in) :: wFlag
    logical,        dimension(:,:),           intent(in) :: mask
    real(sp),       dimension(:,:),           intent(in) :: SM_est
    real(dp),       dimension(:,:),           intent(in) :: opt_h
    integer(i4),                              intent(in) :: yStart
    integer(i4),                    optional, intent(in) :: d            ! optional, duration
    real(sp),       dimension(:,:), optional, intent(in) :: SMI_eval
    real(dp),       dimension(:,:),           intent(in) :: lats, lons   ! latitude and longitude fields of input

    ! local Variables
    character(256)                                :: Fname
    integer(i4), target                           :: m                  ! netCDF counter

    integer(i4)                                   :: nMonths       ! number of simulated months
    integer(i4),    dimension(:),     allocatable :: times
    integer(i4),    dimension(:,:,:), allocatable :: Ziu           ! field integer unpacked 

    !               dimension names
    character(256), dimension(3)                  :: dnames
    character(256), dimension(3)                  :: dims_hopt

    ! initialize time
    nMonths = size( SM_est, 2 )
    allocate( times( nMonths ) )
    forall( m = 1:nMonths ) times( m ) = m
    ! initialize dimension names
    dnames(1) = 'nrows'
    dnames(2) = 'ncols'
    dnames(3) = 'time'    
    dims_hopt(1) = 'nrows'
    dims_hopt(2) = 'ncols'
    dims_hopt(3) = 'months'
    !
    select case (wFlag)

    case (3)
       ! fname
       fName = trim(outpath)//'mSMIc.nc'
       ! save drought mSMIc (same sequence as in 1)
       call var2nc( fName, SMIc, dnames, v_name = 'mSMIc', &
            long_name = 'monthly SMI indicator SMI < th', units = '-', &
            missing_value = nodata_i4, create = .true. ) ! &
            ! scale_factor = 1., coordinates = 'lon lat' )
       ! add lat and lon
       call var2nc( fname, lats, dnames(1:2), v_name = 'lat', &
            long_name = 'longitude', units = 'degrees_east' )
       call var2nc( fname, lons, dnames(1:2), v_name = 'lon', &
            long_name = 'latitude', units = 'degrees_north' )
    case (4)
       ! fname
       fName = trim(outpath)//'mDCluster.nc'
       ! save drought clusters (same sequence as in 1)
       call var2nc( fName, idCluster, dnames, v_name = 'mDC', &
            long_name = 'consolidated cluster evolution', units = '-', &
            missing_value = nodata_i4, create = .true. ) ! &
            ! scale_factor = 1., coordinates = 'lon lat' )
       deallocate( Ziu )
       ! add lat and lon
       call var2nc( fname, lats, dnames(1:2), v_name = 'lat', &
            long_name = 'longitude', units = 'degrees_east' )
       call var2nc( fname, lons, dnames(1:2), v_name = 'lon', &
            long_name = 'latitude', units = 'degrees_north' )
    case (5)
       ! fname
       write (fName, '(i2.2)') d
       fName = trim(outpath)//'sev_'//trim(fName)//'.nc'
       ! save Severity for a given duration d (same sequence as in 1)
       call var2nc( fName, severity, dnames, v_name = 'Severity', &
            long_name = 'd-month severity', units = '-', &
            missing_value = nodata_dp, create = .true. ) ! &
            ! scale_factor = 1., coordinates = 'lon lat' )
       ! add lat and lon
       call var2nc( fname, lats, dnames(1:2), v_name = 'lat', &
            long_name = 'longitude', units = 'degrees_east' )
       call var2nc( fname, lons, dnames(1:2), v_name = 'lon', &
            long_name = 'latitude', units = 'degrees_north' )
    end select

  end subroutine WriteNetCDF

!**********************************************************************
!    PURPOSE    WRITE Results of the cluster analysis
!
!    AUTHOR:    Luis E. Samaniego-Eguiguren, UFZ
!    UPDATES
!               Created        Sa   17.03.2011
!               Last Update    Sa
!**********************************************************************
subroutine WriteResultsCluster(outpath, wFlag, yStart, yEnd, nMonths, nCells, d)

  use mo_kind,          only : i4, dp
  use mo_smi_constants, only : YearMonths

  !
  implicit none

  ! input variables
  character(len=*),      intent(in) :: outpath     ! ouutput path for results
  integer(i4),           intent(in) :: wFlag
  integer(i4),           intent(in) :: yStart
  integer(i4),           intent(in) :: yEnd
  integer(i4),           intent(in) :: nMonths
  integer(i4),           intent(in) :: nCells
  integer(i4), optional, intent(in) :: d

  ! local variables
  real(dp)                  :: pDArea
  !
  integer(i4)               :: i, j, t, k, y, m
  character(len=256)        :: FMT
  character(len=256)        :: fName

  select case (wFlag)
    case (1)
      ! main statistics
      fName = trim(outpath) // 'results_SAD.txt'
      open  (10, file = fName, status='unknown')
      write (10, 100 ) 'i', 'c_Id', 'aDD', 'aDA', 'TDM'
      do i=1,nClusters
        write (10,110)  i, shortCnoList(i), aDD(i), aDA(i), TDM(i)
      end do
      close (10)
      !
      fName = trim(outpath) // 'DArea_evol.txt'
      open  (11, file = fName, status='unknown')
      write (FMT, 120) '(a5,', nClusters, '(2x,a2,i6.6))'
      write (11, FMT ) 'm', ('c_', shortCnoList(i), i=1,nClusters)
      write (FMT, 130) '(i5,', nClusters, 'f10.3)'
      do t=1, nMonths
        write (11,FMT) t, (DAreaEvol(t,i), i=1,nClusters)
      end  do
      close (11)
      !
      fName = trim(outpath) // 'TDM_evol.txt'
      open  (12, file = fName, status='unknown')
      write (FMT, 120) '(a5,', nClusters, '(2x,a2,i6.6))'
      write (12, FMT ) 'm', ('c_', shortCnoList(i), i=1,nClusters)
      write (FMT, 130) '(i5,', nClusters, 'f10.3)'
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
      write (14, 150 ) 'year', 'month', '%AreaDE'
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
      write (15, 170 ) 'year', 'month', '%cAreaDE','SevDE' 
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
      write (21,FMT) ( real(i * deltaArea, dp)*cellArea, (SADperc(i,j), j=1,nQProp), i=1, nInterArea)

      close (21)
  end select

  100 format (2a10, 3a15)
  110 format (2i10, 4f15.3)
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
subroutine WriteResultsBasins( outpath, mask, yStart, yEnd, nMonths, Basin_Id )

  use mo_kind,          only : i4
  use mo_smi_constants, only : nodata_dp, YearMonths

  implicit none

  ! input variables
  character(len=*),            intent(in) :: outpath     ! ouutput path for results
  logical,     dimension(:,:), intent(in) :: mask
  integer(i4),                 intent(in) :: yStart
  integer(i4),                 intent(in) :: yEnd
  integer(i4),                 intent(in) :: nMonths     ! number of simulated months
  integer(i4), dimension(:,:), intent(in) :: Basin_Id    ! IDs for basinwise drought analysis

  ! local variables
  character(len=256)                         :: fName
  integer(i4)                                :: i, m, y, j
  real(dp),    dimension(:,:), allocatable   :: Basin_SMI        

  !-------------------------
  ! basin wise
  !-------------------------
  allocate( Basin_SMI (nMonths, nBasins+1) )
  Basin_SMI = nodata_dp
  !
  do m = 1, nMonths
    do i = 1, nBasins
      !
      Basin_SMI(m,i) = sum( SMI(:,:,m), Basin_Id(:,:) == i ) / real(count(Basin_Id(:,:) == i), dp)
      !
    end do
    Basin_SMI(m,nBasins+1) = sum(SMI(:,:,m), mask ) /  &
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
