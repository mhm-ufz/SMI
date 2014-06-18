!**********************************************************************
!  MODULE InputOutput                                                 *
!  PURPOSE     Input / output subroutines                             *
!  CREATED     Luis Samaniego, 09.02.2011                             *
!                                                                     *
!**********************************************************************
module InputOutput
  use mo_kind,   only                              : i4, sp, dp
  use mo_ncread, only                              : Get_NcVar, Get_NcDim
  implicit none
  character(256)                                  :: headerfName
  character(256)                                  :: basinfName
  character(256)                                  :: DataPathIn
  character(256)                                  :: DataPathOut
  character(256)                                  :: optPDFparfName    ! opt parameter PDF file 
  integer(i4)                                     :: yStart            ! starting year
  integer(i4)                                     :: yEnd              ! ending year
  integer(i4)                                     :: nYears            ! number of years
  integer(i4)                                     :: nMonths           ! number of simulated months
  integer(i4)                                     :: nCells            ! number of effective cells
  integer(i4), parameter                          :: nMy = 12          ! number of months per year
  ! fields
  real(sp),    dimension(:,:), allocatable        :: SM_est            ! monthly fields packed for estimation
  real(sp),    dimension(:,:), allocatable        :: SM_eval           ! monthly fields packed for evaluation
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
  integer(i4)                                    :: SMI_flag
  real(dp)                                       :: SMI_thld               ! SMI threshold
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
  integer(i4), dimension(:,:), allocatable       :: Basin_Id
  integer(i4), parameter                         :: nBasins = 6
  real(dp), dimension(:,:), allocatable          :: Basin_SMI        
  ! netCDF
  integer(i4), dimension(5)                      :: dl                      ! for displaying dimensions
  character(256)                                 :: varNameSM               ! variable 1
  character(256)                                 :: varNameBasinId          ! variable 2
  character(256)                                 :: varNameMask             ! variable 3
  real(dp),  dimension(:,:), allocatable, target :: lats, lons
  integer(i4), dimension(:), allocatable         :: eMask                   ! error field
contains
  !
  !*********************************************************************
  !    SUBROUTINE Read Database
  !    PURPOSE    Reads main file
  !    AUTHOR:    Luis E. Samaniego-Eguiguren, UFZ 23.05.2007
  !    NOTES:
  !               packed fields are stored in dim1->dim2 sequence 
  !*********************************************************************
  subroutine ReadDataMain( do_cluster, cal_SMI, opt_h, basin_flag, mask, nodata )
    use mo_kind,         only : i4
    use kernelSmoother,  only : flagKernelType, nInter, offSet, flagCDF, hOptDB  
    implicit none
    !
    ! input / output Variables
    logical,                              intent(out) :: do_cluster ! do cluster calculation
    logical,                              intent(out) :: cal_SMI    ! should SMI be calculated
    logical,                              intent(out) :: opt_h      ! optimize kernel width
    logical,                              intent(out) :: basin_flag ! basin flag
    logical, dimension(:,:), allocatable, intent(out) :: mask ! grid mask
    real(sp),                             intent(out) :: nodata

    !
    ! local Variables
    logical                   :: monthly_flag ! indicate whether monthly data is given
    integer(i4)               :: i, j, iu, ic, im, k, ios
    integer(i4)               :: ii
    character(256)            :: dummy
    ! directories and filenames
    character(256)            :: maskfName
    character(256)            :: opt_h_file
    character(256)            :: SM_eval_file
    
    ! variable names in netcdf input files
    character(256)            :: mask_vname
    character(256)            :: SM_vname
    character(256)            :: basin_vname

    real(sp), dimension(:,:),   allocatable    :: dummy_D2_sp
    real(dp), dimension(:,:,:), allocatable    :: dummy_D3_dp
    real(sp), dimension(:,:), allocatable      :: ZiT   ! field integer unpacked transposed
    real(sp), dimension(:,:,:), allocatable    :: ZrT   ! field float unpacked transposed
    real(sp), dimension(:,:,:,:), allocatable  :: Zr2T  ! field float unpacked transposed

    ! read main config
    namelist/mainconfig/basin_flag, basinfName, basin_vname, maskfName, mask_vname, DataPathIn, &
         SM_vname, yStart, yEnd, DataPathOut, monthly_flag, SMI_flag, cal_SMI, &
         opt_h, opt_h_file, SM_eval_file, &
         flagKerneltype, nInter, offSet, flagCDF, SMI_thld, nodata, do_cluster
    !   
    open (unit=10, file='main.txt', status='old')
    read(10, mainconfig)
    close (10)
    print*, 'Main.dat read ...ok'

    ! read Mask
    call Get_ncVar( maskfName, trim(mask_vname), dummy_D2_sp )
    allocate( mask( size(dummy_D2_sp,1), size(dummy_D2_sp,2) ) )
    mask = merge( .true., .false., (dummy_D2_sp > nodata ) )
    deallocate( dummy_D2_sp )
    print*, 'mask read ...ok'

    ! read basin mask
    if ( basin_flag ) then
       call Get_ncVar( basinfName, trim(basin_vname), Basin_Id )
       ! consistency check
       if ( ( size( Basin_Id, 1) .ne. size( mask, 1 ) ) .or. &
            ( size( Basin_Id, 2) .ne. size( mask, 2 ) ) ) then
          print *, '***ERROR: size mismatch between basin field and given mask file'
          stop
       end if
       ! intersect mask and basin_id
       mask     = merge( .true., .false., mask .and. ( Basin_Id .gt. nodata ) )
       Basin_Id = merge( Basin_Id, int(nodata,i4), mask )
    end if

    nCells   = count( mask )
    ! consistency check
    if ( nCells .eq. 0 ) then
       print *, '***ERROR: no cell selected in mask'
       stop
    end if
    print*, 'basin ID read ...ok'
    
    ! read SM field
    call Get_ncVar( DataPathIn, trim(SM_vname), dummy_D3_dp )
    ! consistency check
    if ( ( size( dummy_D3_dp, 1) .ne. size( mask, 1 ) ) .or. &
         ( size( dummy_D3_dp, 2) .ne. size( mask, 2 ) ) ) then
       print *, '***ERROR: size mismatch between SM field and given mask file'
       stop
    end if
    allocate( SM_est( nCells, size( dummy_D3_dp, 3 ) ) )
    do ii = 1, size( dummy_D3_dp, 3 )
       SM_est(:,ii) = pack( real(dummy_D3_dp(:,:,ii),sp), mask )
    end do
    
    ! read lats and lon from file
    call Get_ncVar( DataPathIn, 'lat', lats )
    call Get_ncVar( DataPathIn, 'lon', lons )

    ! check if monthly data has been given
    if ( .not. monthly_flag ) then
       print *, '***ERROR: daily averaging not implemented yet'
       stop
    end if
    print*, 'soil moisture field read ...ok'

    ! check whether second SM field and optimized h should be written
    if ( .not. cal_SMI ) then
       print *, '***ERROR: reading of optimized kernel width and SM evaluation not implemented!'
       stop
    end if

    print *, '***WARNING: yStart and yEnd are currently not considered'

    ! initialize further Variables
    nMonths = size( SM_est, 2 )
    nYears  = nMonths / nMy

  end subroutine ReadDataMain

  !**********************************************************************
  !    PURPOSE    WRITE Results of the kernel smoother
  !    FORMAT     ascii tables
  !
  !    AUTHOR:    Luis E. Samaniego-Eguiguren, UFZ
  !    UPDATES
  !               Created        Sa   21.03.2006
  !               Last Update    Sa   
  !**********************************************************************
  subroutine WriteResultsKernel(wFlag,iCell,iMonth) 
    use mo_kind,          only : i4
    use kernelSmoother,   only : edf, pdf, nInter, nObs, hOpt, flagKernelType, flagCDF, hOptDB
    !
    implicit none
    integer(i4), intent (in)  :: wFlag
    integer(i4), intent (in)  :: iCell
    integer(i4), intent (in)  :: iMonth
    !
    integer(i4)               :: i, j
    character(len=20)         :: dummy
    character(len=256)        :: fName
    !
    select case (wFlag)
    case (1)
       ! print *, 'Kernel density for cell   ', iCell, ' was estimated ... OK'
       write (dummy, '(i1,a,i5.5,a,i2.2,a4)')  flagKernelType, '_', icell, '_', iMonth, '.txt'
       fName = trim(dataPathOut) // '/cdf_opt_' // trim(dummy)
       open(10, file=fName, status='unknown')
       write(10, 100)  'x', 'pdf(x)'
       select case (flagCDF)
       case (.TRUE.)
          ! CDF
          do i = 1, nInter
             write (10, 101) (pdf(i,j), j=1,3,2)
          end do
       case (.FALSE.)
          ! PDF
          do i = 1, nInter
             write (10, 101) (pdf(i,j), j=1,2)
          end do
       end select
       close (10)  

    case (2)
       ! print *, 'Empirical density estimated   ... OK'
       write (dummy, '(i5.5,a,i2.2,a4)')  icell, '_', iMonth, '.txt'

       fName = trim(dataPathOut) // '/edf_' // trim(dummy)
       open(20, file=fName, status='unknown')
       write(20, 100)  'x', 'edf(x)'
       write (20, 101)  edf(0,1), edf(0,2)
       do i = 1, nObs
          write (20, 101)  edf(i,1), edf(i-1,2)
          write (20, 101) (edf(i,j), j=1,2)
       end do
       write (20, 101)  edf(nObs+1,1), edf(nObs+1,2)
       close (20)  

    case (3)
       ! print *, 'Optimised bandwidth   ... OK'
       fName = trim(dataPathOut)// '/opti.txt'
       open (30, file=fName, status='unknown', position='APPEND')
       write(30, 300) icell, iMonth, flagKernelType, hOPt
       close(30)  

    case (4)       
       ! print *, 'Optimised bandwidth   ... OK'
       fName = trim(dataPathOut)// '/opti.txt'
       open (30, file=fName, status='unknown')
 !      write(30, 300) icell, iMonth, flagKernelType, hOPt
       do i= 1, nCells
          do j = 1, nMy
             if ( hOptDB(i,j) < 0.0_dp  ) cycle
             write (30,300) i, j, flagKernelType, hOptDB(i,j)
          end do
       end do
       close(30) 
    end select
    !
    ! formats (variable)
100 format (  2a15 )  
101 format (  2e15.4 )
300 format (  i10, i5, i5, es15.5)
  end subroutine WriteResultsKernel

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
  subroutine WriteNetCDF(wFlag, mask, nodata, d) 
    use mo_kind, only              : i4, sp
    !
    !use netCDF_varDef
    use mo_ncwrite,     only: var2nc
    use kernelSmoother, only: hOptDB
    implicit none
    !
    !implicit none
    logical, dimension(:,:), intent(in)    :: mask
    real(sp),                intent(in)    :: nodata
    ! local Variables
    character(256)                         :: Fname
    integer(i4), intent (in), optional     :: d                ! optional, duration
    integer(i4), intent (in)               :: wFlag            !
    integer(i4), target                    :: m                ! netCDF counter
    integer(i4), save                      :: ncId             ! netCDF ID handler
    !
    integer(i4), dimension(:,:,:), allocatable  :: Ziu           ! field integer unpacked 
    real(sp),    dimension(:,:,:), allocatable  :: Zu            ! field real unpacked 
    real(dp),    dimension(:,:,:), allocatable  :: dummy_D3_dp   ! field real unpacked 

    integer(i4)                            :: iLoc(2)
    !
    ! dimension names
    character(256), dimension(3)           :: dnames
    character(256), dimension(3)           :: dims_hopt
    ! initialize dimension names
    dnames(1) = 'nrows'
    dnames(2) = 'ncols'
    dnames(3) = 'time'    
    dims_hopt(1) = 'nrows'
    dims_hopt(2) = 'ncols'
    dims_hopt(3) = 'months'
    !
    select case (wFlag)
    case (1)
       ! fname
       fName = trim(DataPathOut)//'mean_m_SM.nc'
       ! unpack estimated SM
       allocate( Zu( size( mask, 1), size( mask, 2), size( SM_est, 2) ) )
       do m = 1, nMonths
          Zu( :, :, m) = unpack ( real( SM_est(:,m), sp ), mask, nodata )
       end do
       ! save monthly soil moisture averages
       call var2nc( fName, Zu, dnames, v_name = 'mSMmean', &
            longname = 'monthly mean SM/SMs', units = '%', &
            fill_value = nodata, f_exists = .false. ) ! &
            ! scale_factor = 1., coordinates = 'lon lat' )
       deallocate( Zu )
       ! add lat and lon
       call var2nc( fname, lats, dnames(1:2), v_name = 'lat', &
            longname = 'longitude', units = 'degrees_east' )
       call var2nc( fname, lons, dnames(1:2), v_name = 'lon', &
            longname = 'latitude', units = 'degrees_north' )
    case (2)
       ! fname
       fName = trim(DataPathOut)//'mSMI.nc'
       ! unpack estimated SMIp
       allocate( Zu( size( mask, 1), size( mask, 2), size( SMIp, 2) ) )
       do m = 1, nMonths
          Zu( :, :, m) = unpack ( real( SMIp(:,m), sp ), mask, nodata )
       end do
       ! save SMIp (same sequence as in 1)
       call var2nc( fName, Zu, dnames, v_name = 'SMI', &
            longname = 'monthly soil moisture index', units = '-', &
            fill_value = nodata, f_exists = .false. ) ! &
            ! scale_factor = 1., coordinates = 'lon lat' )
       deallocate( Zu )
       ! write out kernel width if it has been optimised
       allocate( dummy_D3_dp( size( mask, 1 ), size( mask, 2), size( hOptDB, 2 ) ) )
       do m = 1, nMy
          dummy_D3_dp( :, :, m ) = unpack( hOptDB(:,m), mask, real( nodata, dp ) ) 
       end do
       call var2nc( fname, dummy_D3_dp, dims_hopt, v_name = 'h_opt', &
            longname = 'optimised kernel width', units = '-', &
            fill_value = real( nodata, dp ) )
       deallocate( dummy_D3_dp ) 
       ! add lat and lon
       call var2nc( fname, lats, dnames(1:2), v_name = 'lat', &
            longname = 'longitude', units = 'degrees_east' )
       call var2nc( fname, lons, dnames(1:2), v_name = 'lon', &
            longname = 'latitude', units = 'degrees_north' )
    case (3)
       ! fname
       fName = trim(DataPathOut)//'mSMIc.nc'
       ! save drought mSMIc (same sequence as in 1)
       call var2nc( fName, SMIc, dnames, v_name = 'mSMIc', &
            longname = 'monthly SMI indicator SMI < th', units = '-', &
            fill_value = int(nodata,i4), f_exists = .false. ) ! &
            ! scale_factor = 1., coordinates = 'lon lat' )
       ! add lat and lon
       call var2nc( fname, lats, dnames(1:2), v_name = 'lat', &
            longname = 'longitude', units = 'degrees_east' )
       call var2nc( fname, lons, dnames(1:2), v_name = 'lon', &
            longname = 'latitude', units = 'degrees_north' )
    case (4)
       ! fname
       fName = trim(DataPathOut)//'mDCluster.nc'
       ! save drought clusters (same sequence as in 1)
       call var2nc( fName, idCluster, dnames, v_name = 'mDC', &
            longname = 'consolidated cluster evolution', units = '-', &
            fill_value = int(nodata,i4), f_exists = .false. ) ! &
            ! scale_factor = 1., coordinates = 'lon lat' )
       deallocate( Ziu )
       ! add lat and lon
       call var2nc( fname, lats, dnames(1:2), v_name = 'lat', &
            longname = 'longitude', units = 'degrees_east' )
       call var2nc( fname, lons, dnames(1:2), v_name = 'lon', &
            longname = 'latitude', units = 'degrees_north' )
    case (5)
       ! fname
       write (fName, '(i2.2)') d
       fName = trim(DataPathOut)//'sev_'//trim(fName)//'.nc'
       ! save Severity for a given duration d (same sequence as in 1)
       call var2nc( fName, severity, dnames, v_name = 'Severity', &
            longname = 'd-month severity', units = '-', &
            fill_value = real(nodata,dp), f_exists = .false. ) ! &
            ! scale_factor = 1., coordinates = 'lon lat' )
       deallocate( Zu )
       ! add lat and lon
       call var2nc( fname, lats, dnames(1:2), v_name = 'lat', &
            longname = 'longitude', units = 'degrees_east' )
       call var2nc( fname, lons, dnames(1:2), v_name = 'lon', &
            longname = 'latitude', units = 'degrees_north' )
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
subroutine WriteResultsCluster(wFlag,d)
  use mo_kind, only                    : i4, dp
  !
  implicit none
  integer(i4), intent (in)            :: wFlag
  integer(i4), intent (in), optional  :: d
  real(dp)                            :: pDArea
  !
  integer(i4)               :: i, j, t, k, y, m
  character(len=256)        :: dummy, FMT
  character(len=256)        :: fName
  !
  select case (wFlag)
    case (1)
      ! main statistics
      fName = trim(dataPathOut) // 'results_SAD.txt'
      open  (10, file = fName, status='unknown')
      write (10, 100 ) 'i', 'c_Id', 'aDD', 'aDA', 'TDM'
      do i=1,nClusters
        write (10,110)  i, shortCnoList(i), aDD(i), aDA(i), TDM(i)
      end do
      close (10)
      !
      fName = trim(dataPathOut) // 'DArea_evol.txt'
      open  (11, file = fName, status='unknown')
      write (FMT, 120) '(a5,', nClusters, '(2x,a2,i6.6))'
      write (11, FMT ) 'm', ('c_', shortCnoList(i), i=1,nClusters)
      write (FMT, 130) '(i5,', nClusters, 'f10.3)'
      do t=1, nMonths
        write (11,FMT) t, (DAreaEvol(t,i), i=1,nClusters)
      end  do
      close (11)
      !
      fName = trim(dataPathOut) // 'TDM_evol.txt'
      open  (12, file = fName, status='unknown')
      write (FMT, 120) '(a5,', nClusters, '(2x,a2,i6.6))'
      write (12, FMT ) 'm', ('c_', shortCnoList(i), i=1,nClusters)
      write (FMT, 130) '(i5,', nClusters, 'f10.3)'
      do t=1, nMonths
        write (12,FMT) t, (DTMagEvol(t,i), i=1,nClusters)
      end  do
      close (12)
      !
      fName = trim(dataPathOut) // 'event_ids.txt'
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
      fName = trim(dataPathOut) // 'DArea_evol_total.txt'
      open  (14, file = fName, status='unknown')
      write (14, 150 ) 'year', 'month', '%AreaDE'
      t = 0
      do y =yStart, yEnd 
        do m = 1, nMy
          t = t + 1
          pdArea = real( count( SMIc(:,:,t) == 1 ), dp) / &
                   real( nCells, dp) * 1e2_dp
          write(14,160) y, m, pdArea
        end do
      end do
      close (14)
      !
      ! Time evolution of the cluster area (less than SMIc) and monthly severity
      fName = trim(dataPathOut) // 'DcArea_sev_evol.txt'
      open  (15, file = fName, status='unknown')
      write (15, 170 ) 'year', 'month', '%cAreaDE','SevDE' 
      t = 0
      do y =yStart, yEnd 
        do m = 1, nMy
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
         fName = trim(dataPathOut) // trim(fName)
         open (20, file=fName, status='unknown')
         write (20,210) 'Area[km2]', 'Severity'
         write (20,220) (( SAD(i, j, k ), j=1,2 ), i=1, nInterArea)
         close(20)
      end do

      ! SAD percentiles for a given duration
      write (fName,230) 'SAD_perc_', d, '.txt'
      fName = trim(dataPathOut) // trim(fName)
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
subroutine WriteResultsBasins( mask, nodata )
  use mo_kind, only          : i4
  !
  implicit none
  !
  logical, dimension(:,:), intent(in) :: mask
  integer(i4),             intent(in) :: nodata
  integer(i4)               :: i, m, k, y, j
  character(len=256)        :: fName
  !-------------------------
  ! basin wise
  !-------------------------
  allocate( Basin_SMI (nMonths, nBasins+1) )
  Basin_SMI = nodata
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
  fName = trim(dataPathOut) // 'basin_avg_SMI.txt'
  open(20, file =fName , status = 'unknown')
  write(20, 1) 'Month_No', 'year', 'month', ('Basin_', i, i = 1, nBasins), 'Germany_Avg'
  !
  fName = trim(dataPathOut) // 'basin_avg_dArea.txt'
  open(21, file =fName , status = 'unknown')
  write(21, 3) 'Month_No', 'year', 'month', ('Basin_', i, i = 1, nBasins)
  ! 
  fName = trim(dataPathOut) // 'basin_avg_sev.txt'
  open(22, file =fName , status = 'unknown')
  write(22, 3) 'Month_No', 'year', 'month', ('Basin_', i, i = 1, nBasins)
  !! NEW !!
  m=0
  do y = yStart, yEnd
     do j = 1, nMy
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
