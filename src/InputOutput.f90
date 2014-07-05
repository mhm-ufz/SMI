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
  character(256)                                  :: DataPathIn
  character(256)                                  :: DataPathOut
  character(256)                                  :: optPDFparfName    ! opt parameter PDF file 
  integer(i4)                                     :: nYears            ! number of years
  integer(i4)                                     :: nMonths           ! number of simulated months
  integer(i4)                                     :: nCells            ! number of effective cells
  integer(i4), parameter                          :: nMy = 12          ! number of months per year
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
  subroutine ReadDataMain( do_cluster, eval_SMI, read_opt_h, silverman_h, opt_h, basin_flag, mask, &
       SM_est, tmask_est, SM_eval, tmask_eval, yStart, yEnd, mstart, time_sp, nodata, offSet )
    use mo_kind,         only: i4
    use mo_utils,        only: equal
    use kernelSmoother,  only: flagKernelType, nInter, flagCDF   
    use mo_string_utils, only: DIVIDE_STRING

    implicit none
    !
    ! input / output Variables
    logical,                               intent(out) :: do_cluster ! do cluster calculation
    logical,                               intent(out) :: eval_SMI   ! should SMI be calculated
    logical,                               intent(out) :: read_opt_h ! read kernel width
    logical,                               intent(out) :: silverman_h ! optimize kernel width
    logical,                               intent(out) :: basin_flag ! basin flag
    logical,  dimension(:,:), allocatable, intent(out) :: mask       ! grid mask
    integer(i4),                           intent(out) :: yStart     ! starting year
    integer(i4),                           intent(out) :: yEnd       ! ending year
    integer(i4),                           intent(out) :: mStart     ! starting year
    real(sp), dimension(:),   allocatable, intent(out) :: time_sp
    real(sp),                              intent(out) :: nodata
    real(sp), dimension(:,:), allocatable, intent(out) :: SM_est     ! monthly fields packed for estimation
    logical,  dimension(:,:), allocatable, intent(out) :: tmask_est  ! temporal mask of estimated arr
    real(sp), dimension(:,:), allocatable, intent(out) :: SM_eval    ! monthly fields packed for evaluation
    logical,  dimension(:,:), allocatable, intent(out) :: tmask_eval ! temporal mask of evaluated arr
    real(dp), dimension(:,:), allocatable, intent(out) :: opt_h
    real(dp),                              intent(out) :: offSet ! shift starting/ending x

    !
    ! local Variables
    logical                   :: monthly_flag ! indicate whether monthly data is given
    integer(i4)               :: ii
    integer(i4)               :: mm
    ! directories and filenames
    character(256)            :: maskfName
    character(256)            :: opt_h_file
    character(256)            :: SM_eval_file
    character(256)            :: basinfName
    
    ! variable names in netcdf input files
    character(256)            :: mask_vname
    character(256)            :: SM_vname
    character(256)            :: SM_eval_vname
    character(256)            :: basin_vname
    character(256)            :: opt_h_vname
    character(256)            :: time_units

    character(256), dimension(:),       allocatable :: strArr  ! dummy for netcdf attribute handling
    real(sp),       dimension(:,:),     allocatable :: dummy_D2_sp
    real(sp),       dimension(:,:,:),   allocatable :: dummy_D3_sp
    real(dp),       dimension(:,:,:),   allocatable :: dummy_D3_dp

    ! read main config
    namelist/mainconfig/basin_flag, basinfName, basin_vname, maskfName, mask_vname, DataPathIn, &
         SM_vname, yStart, yEnd, DataPathOut, monthly_flag, SMI_flag, eval_SMI, &
         silverman_h, read_opt_h, opt_h_vname, opt_h_file, SM_eval_file, SM_eval_vname, &
         flagKerneltype, nInter, offSet, flagCDF, SMI_thld, nodata, do_cluster
    !   
    open (unit=10, file='main.dat', status='old')
    read(10, mainconfig)
    close (10)
    print*, 'main.dat read ...ok'

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
    deallocate( dummy_D3_dp )
    ! determine time mask 
    print *, '***CAUTION: time axis assumed to be starting with 1!'
    call Get_NcVarAtt( DataPathIn, 'time', 'units', time_units )
    call DIVIDE_STRING(trim(time_units), ' ', strArr)
    call DIVIDE_STRING(trim(strArr(3)), '-', strArr) 
    read( strArr(1), * ) yStart
    read( strArr(2), * ) mStart
    allocate( time_sp( size( SM_est, 2 ) ) )
    forall( mm = 1 : size( time_sp, 1 ) ) time_sp(mm) = real(mm,sp)
    allocate ( tmask_est( size(SM_est, 2), nMy ) )
    tmask_est = .false.
    do mm = 1, nMy
       tmask_est(:,mm) = ( mod( int(time_sp,i4) + mStart, nMy ) .eq. mod( mm, nMy ) )
    end do
    if ( any( count( tmask_est, dim = 1 ) .eq. 0_i4 ) ) &
         stop '***ERROR no data for estimation given for all calendar months, check time axis'
    
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
    if ( eval_SMI ) then
       ! read second SM field that uses CDF of the first one
       call Get_ncVar( trim( SM_eval_file ), trim( SM_eval_vname ), dummy_D3_sp )
       allocate( SM_eval( nCells, size( dummy_D3_sp, 3 ) ) )
       do ii = 1, size( dummy_D3_sp, 3 )
          SM_eval(:, ii) = pack( dummy_D3_sp(:,:,ii), mask)
       end do
       deallocate( dummy_D3_sp )
       ! read time axis
       if ( allocated( time_sp ) ) deallocate( time_sp )
       call Get_ncVar( trim( SM_eval_file ), 'time', time_sp )
       call Get_ncVarAtt( trim( SM_eval_file), 'time', 'units', time_units )
       call DIVIDE_STRING(trim(time_units), ' ', strArr)
       call DIVIDE_STRING(trim(strArr(3)), '-', strArr) 
       read(strArr(1),*) yStart
       read(strArr(2),*) mStart
       ! determine time mask
       allocate ( tmask_eval( size(SM_eval, 2), nMy ) )
       tmask_eval = .false.
       do mm = 1, nMy
          tmask_eval(:,mm) = ( mod( int(time_sp,i4) + mStart, nMy ) .eq. mod( mm, nMy ) )
       end do
       if ( all( count( tmask_eval, dim = 1 ) .eq. 0_i4 ) ) &
            stop '***ERROR no data in eval given, check time axis'
       print*, 'read soil moisture field for evaluation... ok'
    end if

    ! initialize opt_h
    allocate ( opt_h( ncells, nMy ) )
    opt_h = nodata
    !
    if ( read_opt_h ) then
       ! read optimized kernel width from file
       call Get_ncVar( trim( opt_h_file ), trim( opt_h_vname ), dummy_D3_dp )
       do ii = 1, size( dummy_D3_dp, 3 )
          opt_h( :, ii ) = pack( real( dummy_D3_dp(:,:,ii),sp ), mask )
       end do
       deallocate( dummy_D3_dp )
       if ( any( equal( opt_h, real(nodata, dp)) ) ) &
            stop '***ERROR kernel width contains nodata values'
       print *, 'read kernel width from file... ok'
    end if    

    print *, '***WARNING: yStart and yEnd are only considered for output writing'

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
  subroutine WriteResultsKernel(wFlag, iCell, iMonth, opt_h) 
    use mo_kind,          only : i4
    use kernelSmoother,   only : edf, pdf, nInter, nObs, hOpt, flagKernelType, flagCDF
    !
    implicit none
    integer(i4),              intent(in) :: wFlag
    integer(i4),              intent(in) :: iCell
    integer(i4),              intent(in) :: iMonth
    real(sp), dimension(:,:), intent(in) :: opt_h
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
             if ( opt_h(i,j) < 0.0_dp  ) cycle
             write (30,300) i, j, flagKernelType, opt_h(i,j)
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


  ! ##################################################################
  ! subroutine for writting the SMI to nc file
  ! author: Stephan Thober
  ! created: 5.7.2014
  ! ##################################################################
  subroutine WriteSMI( SMI, mask, nodata, yStart, mStart, time, hh ) 
    !
    use mo_kind,         only: i4, sp
    use mo_string_utils, only: num2str
    use mo_ncwrite,      only: var2nc
    !
    implicit none
    !

    ! input variables
    real(sp), dimension(:,:),           intent(in) :: SMI
    logical, dimension(:,:),            intent(in) :: mask
    real(sp),                           intent(in) :: nodata
    real(dp), dimension(:,:), optional, intent(in) :: hh
    integer(i4),                        intent(in) :: yStart
    integer(i4),                        intent(in) :: mStart
    real(sp), dimension(:), allocatable,intent(in) :: time

    ! local Variables
    character(256)                              :: Fname
    integer(i4)                                 :: mm            ! month
    integer(i4), dimension(:,:,:), allocatable  :: Ziu           ! field integer unpacked 
    real(sp),    dimension(:,:,:), allocatable  :: dummy_D3_sp   ! field real unpacked 
    real(dp),    dimension(:,:,:), allocatable  :: dummy_D3_dp   ! field real unpacked 
    ! dimension names
    character(256), dimension(3)                :: dnames
    character(256), dimension(3)                :: dims_hopt

    ! initialize dimension names
    dnames(1) = 'nrows'
    dnames(2) = 'ncols'
    dnames(3) = 'time'    
    dims_hopt(1) = 'nrows'
    dims_hopt(2) = 'ncols'
    dims_hopt(3) = 'months'

    ! fname
    fName = trim(DataPathOut)//'mSMI.nc'

    ! unpack estimated SMIp
    allocate( dummy_d3_sp( size( mask, 1), size( mask, 2), size( SMI, 2) ) )
    do mm = 1, size( SMI, 2 )
       dummy_d3_sp( :, :, mm) = unpack ( SMI(:,mm), mask, nodata )
    end do
    ! save SMI (same sequence as in 1)
    call var2nc( fName, dummy_d3_sp, dnames, v_name = 'SMI', &
         long_name = 'monthly soil moisture index', units = '-', &
         missing_value = nodata, create = .true. ) ! &
    ! scale_factor = 1., coordinates = 'lon lat' )
    deallocate( dummy_d3_sp )

    ! write out kernel width if it has been optimised
    if ( present( hh ) ) then
       allocate( dummy_D3_dp( size( mask, 1 ), size( mask, 2), size( hh, 2 ) ) )
       do mm = 1, nMy
          dummy_D3_dp( :, :, mm ) = unpack( hh(:,mm), mask, real( nodata, dp ) ) 
       end do
       call var2nc( fname, dummy_D3_dp, dims_hopt, v_name = 'kernel_width', &
            long_name = 'optimised kernel width', units = '-', &
            missing_value = real( nodata, dp ) )
       deallocate( dummy_D3_dp ) 
    end if

    ! add lat and lon
    call var2nc( fname, lats, dnames(1:2), v_name = 'lat', &
         long_name = 'longitude', units = 'degrees_east' )
    call var2nc( fname, lons, dnames(1:2), v_name = 'lon', &
         long_name = 'latitude', units = 'degrees_north' )

    ! add time
    call var2nc( fname, time, dnames(3:3), v_name='time', &
         long_name = 'time', units = 'months since ' // &
         trim(num2str(ystart, '(i4)')) // '-' // trim(num2str(mStart,'(i2.2)')) // '-01 00:00:00' )
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
  subroutine WriteNetCDF(wFlag, opt_h, SM_est, mask, nodata, yStart, d, SMI_eval) 
    !
    use mo_kind,         only: i4, sp
    use mo_string_utils, only: num2str
    use mo_ncwrite,      only: var2nc
    !
    implicit none
    !
    ! input variables
    integer(i4), intent (in)                       :: wFlag
    logical, dimension(:,:),            intent(in) :: mask
    real(sp),                           intent(in) :: nodata
    real(sp), dimension(:,:),           intent(in) :: SM_est
    real(dp), dimension(:,:),           intent(in) :: opt_h
    integer(i4),                        intent(in) :: yStart
    integer(i4),              optional, intent(in) :: d                ! optional, duration
    real(sp), dimension(:,:), optional, intent(in) :: SMI_eval
    ! local Variables
    character(256)                       :: Fname
    integer(i4), target                  :: m                ! netCDF counter
    !
    integer(i4), dimension(:),     allocatable  :: time
    integer(i4), dimension(:,:,:), allocatable  :: Ziu           ! field integer unpacked 
    real(sp),    dimension(:,:,:), allocatable  :: dummy_D3_sp   ! field real unpacked 
    real(dp),    dimension(:,:,:), allocatable  :: dummy_D3_dp   ! field real unpacked 

    ! dimension names
    character(256), dimension(3)           :: dnames
    character(256), dimension(3)           :: dims_hopt
    ! initialize time
    allocate( time( nMonths ) )
    forall( m = 1:nMonths ) time( m ) = m
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
       fName = trim(DataPathOut)//'mSMIc.nc'
       ! save drought mSMIc (same sequence as in 1)
       call var2nc( fName, SMIc, dnames, v_name = 'mSMIc', &
            long_name = 'monthly SMI indicator SMI < th', units = '-', &
            missing_value = int(nodata,i4), create = .true. ) ! &
            ! scale_factor = 1., coordinates = 'lon lat' )
       ! add lat and lon
       call var2nc( fname, lats, dnames(1:2), v_name = 'lat', &
            long_name = 'longitude', units = 'degrees_east' )
       call var2nc( fname, lons, dnames(1:2), v_name = 'lon', &
            long_name = 'latitude', units = 'degrees_north' )
    case (4)
       ! fname
       fName = trim(DataPathOut)//'mDCluster.nc'
       ! save drought clusters (same sequence as in 1)
       call var2nc( fName, idCluster, dnames, v_name = 'mDC', &
            long_name = 'consolidated cluster evolution', units = '-', &
            missing_value = int(nodata,i4), create = .true. ) ! &
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
       fName = trim(DataPathOut)//'sev_'//trim(fName)//'.nc'
       ! save Severity for a given duration d (same sequence as in 1)
       call var2nc( fName, severity, dnames, v_name = 'Severity', &
            long_name = 'd-month severity', units = '-', &
            missing_value = real(nodata,dp), create = .true. ) ! &
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
subroutine WriteResultsCluster(wFlag, yStart, yEnd, d)
  use mo_kind, only                    : i4, dp
  !
  implicit none
  ! input variables
  integer(i4), intent (in)            :: wFlag
  integer(i4), intent (in)            :: yStart
  integer(i4), intent (in)            :: yEnd
  integer(i4), intent (in), optional  :: d

  ! local variables
  real(dp)                            :: pDArea
  !
  integer(i4)               :: i, j, t, k, y, m
  character(len=256)        :: FMT
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
subroutine WriteResultsBasins( mask, nodata, yStart, yEnd )
  use mo_kind, only          : i4
  !
  implicit none
  ! input variables
  logical, dimension(:,:), intent(in) :: mask
  integer(i4),             intent(in) :: nodata
  integer(i4),             intent(in) :: yStart
  integer(i4),             intent(in) :: yEnd

  ! local variables
  integer(i4)               :: i, m, y, j
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
