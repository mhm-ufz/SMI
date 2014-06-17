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
  character(256)                                  :: errorfName
  character(256)                                  :: DataPathIn
  character(256)                                  :: DataPathOut
  character(256)                                  :: optPDFparfName    ! opt parameter PDF file 
  integer(i4)                                     :: yStart            ! starting year
  integer(i4)                                     :: yEnd              ! ending year
  integer(i4)                                     :: nYears            ! number of years
  integer(i4)                                     :: nMonths           ! number of simulated months
  integer(i4)                                     :: nCells            ! number of effective cells
  integer(i4), parameter                          :: nMy = 12          ! number of months per year
  !
  type gridGeoRef
     integer(i4)                                  :: ncols             ! number of columns
     integer(i4)                                  :: nrows             ! number of rows
     real(dp)                                     :: xllcorner         ! x coordinate of the lowerleft corner
     real(dp)                                     :: yllcorner         ! y coordinate of the lowerleft corner
     real(dp)                                     :: cellsize          ! cellsize x = cellsize y
     real(dp)                                     :: nodata_value      ! not belonging value
  end type gridGeoRef
  type (gridGeoRef)                               :: grid              ! grid definition for all fields
  !
  ! GRID mask
  integer(i4), dimension(:,:), allocatable        :: mask              ! Value of Mask
  ! fields
  real(dp),    dimension(:,:), allocatable        :: Z                 ! monthly fields packed
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
  real(dp), dimension(nQProp), parameter         :: QProp = (/0.90_dp, 0.96_dp, 0.98_dp /)
                                                                           ! percentiles corresponding
                                                                           ! to return periods of 10,25,50 years
  integer(i4), parameter                         :: nLargerEvents = 600    ! n. largest drought events
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
  real(dp),  dimension(:,:), allocatable         :: lats, lons
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
  subroutine ReadDataMain
    use mo_kind,         only : i4
    use kernelSmoother,  only : flagKernelType, nInter, offSet, flagCDF, hOptDB  
    implicit none
    !
    ! Variables
    integer(i4)               :: i, j, iu, ic, im, k, ios
    character(256)            :: dummy
    real(sp), dimension(:,:), allocatable      :: ZiT   ! field integer unpacked transposed
    real(sp), dimension(:,:,:), allocatable    :: ZrT   ! field float unpacked transposed
    real(sp), dimension(:,:,:,:), allocatable  :: Zr2T  ! field float unpacked transposed

    namelist/mainconfig/headerfName, basinfName, errorfName, DataPathIn, DataPathOut, &
         DAT_flag, SMI_flag, optPDFparfName, flagKerneltype, nInter, offSet, flagCDF, &
         yStart, yEnd, SMI_thld
    !   
    open (unit=10, file='main.txt', status='old')
    read(10, mainconfig)
    close (10)
    print*, 'Main.dat read ...'
    !
    ! ----------------------------------------------------------------------
    !	                               READ MASKs
    ! ----------------------------------------------------------------------
    if (SMI_flag == 4) then
       ! read NetCDF mask
       !
       varNameMask = 'Shape'
       dl = get_NcDim(headerfName, varNameMask)
       grid%ncols = dl(1)
       grid%nrows = dl(2) 
       allocate( mask (grid%nrows, grid%ncols), ZiT (grid%ncols, grid%nrows)  )
       !
       call Get_NcVar(headerfName, varNameMask, ZiT(:,:))
       mask = transpose(ZiT)
       ! some definitions
       grid%cellsize     = 7000._dp                      ! to be defined later
       grid%nodata_value = -9999._dp
       where (mask < 1) mask = int(grid%nodata_value, i4) 
       !
       ! read error mask (exclude cells with numerical problems)
       varNameMask = 'error'
       call Get_NcVar(errorfName, varNameMask, ZiT(:,:))
       where (transpose(ZiT) == 1) mask = int(grid%nodata_value, i4)

!!$       do i=1,grid%nrows
!!$          write (999, '(200i8)') (mask(i,j), j=1,grid%ncols)
!!$       end do
!!$ stop
       !
    else if (SMI_flag == 5) then
       ! read NetCDF mask
       !
       varNameMask = 'mask'
       dl = get_NcDim(headerfName, varNameMask)
       grid%ncols = dl(1)  ! 2
       grid%nrows = dl(2)  ! 3
       allocate( mask (grid%nrows, grid%ncols), ZrT (1,grid%ncols, grid%nrows),  &
                                                Zr2T(1,1,grid%ncols, grid%nrows)  )
       !
       call Get_NcVar(headerfName, varNameMask, ZrT(1,:,:))
       mask = int( transpose(ZrT(1,:,:)), i4)
       print*, trim(headerfName), ' ... read'
       !
       ! some definitions
       grid%cellsize     = 12220._dp                      ! to be defined later
       grid%nodata_value = -9999._dp
       where (mask < 1) mask = int(grid%nodata_value, i4) 
       !
       ! read error mask (exclude cells with numerical problems)
       varNameMask = 'mask'
       call Get_NcVar(errorfName, varNameMask, Zr2T(1,1,:,:))
       where (transpose(Zr2T(1,1,:,:)) > 0) mask = int(grid%nodata_value, i4)
       print*, trim(errorfName), ' ... read'
    else
       call ReadOneHeader(headerfName, iu)
       allocate (mask (grid%nrows , grid%ncols) )
       read (iu, *) (( mask(i,j), j=1,grid%ncols), i=1,grid%nrows )
       close(iu)
    end if
    
    ! ----------------------------------------------------------------------
    !	                  READ FIELDS and ESTIMATE TEMPORAL MEAN
    ! ----------------------------------------------------------------------
    ! STORAGE
    nYears  = yEnd - yStart + 1
    nMonths = nMy * nYears
    nCells  = count ( mask /= int(grid%nodata_value, i4) )
    cellArea   = grid%cellsize**2 / 1.0e6_dp                               ! in km2
    ! allways allocate to store hOpt 
    allocate (hOptDB(nCells,nMy))
    hOptDB = -9.0_dp    
    !
    select case (SMI_flag)
    case(1)
      allocate ( Z(nCells, nMonths) )                                      ! space_dim > time_dim
      call ReadFields_BinP
      print*, 'Soil moisture fields were read... '
      print*, 'Monthly soil moisture fields were estimated... '
    case(2)
      call ReadSMI
      print*, 'Monthly soil moisture index was read... '
    case(3)
      ! 
      allocate ( Z(nCells, nMonths) )                                      ! space_dim > time_dim
      call ReadFields_BinP
      print*, 'Soil moisture fields were read... '
      print*, 'Monthly soil moisture fields were estimated... '
      !
      open(10, file=optPDFparfName, status='old')
      ios = 0
      do while (.NOT. ios /= 0)
        read (10, *,iostat=ios) ic, im, k, hOptDB(ic,im)
      end do
      print*, 'Optimized bandwidth file: '
      print*, trim(optPDFparfName)
      print*, '... was read'
      print*, 'Number of cells-month without parameter hOpt =', count( hOptDB == -9.0_dp)
      print*, '  ... SMI will not be estimated for those grid cells'
    case(4:5)
      ! read NetCDF daily/montly SM fields
      allocate ( Z(nCells, nMonths) )                                      ! space_dim > time_dim
      call ReadFields_nc
      print*, 'Soil moisture fields were read... '
      if (SMI_flag==4) print*, 'Monthly soil moisture fields were estimated... '

    end select
    !
    ! ----------------------------------------------------------------------
    !	                                  BASINS
    ! ----------------------------------------------------------------------
    if (SMI_flag == 4 .or. SMI_flag == 5 ) then
       ! read NetCDF BasinId
       varNameBasinId = 'mask'
       print *, trim( basinfName )
       dl = get_NcDim(basinfName, varNameBasinId)
       if (dl(1) /= grid%ncols .or. dl(2) /= grid%nrows) then
          print*, 'Dimension of variable ', varNameBasinId, ' not correct'
          stop
       end if
       allocate( Basin_Id (grid%nrows, grid%ncols) )
       if ( allocated ( ZiT ) ) deallocate ( ZiT )
       allocate( ZiT (grid%ncols, grid%nrows)  )
       !
       ZiT = 0
       call Get_NcVar(basinfName,varNameBasinId, ZiT(:,:))
       Basin_Id = transpose(ZiT)
       ! some definitions
       where (Basin_Id < 1 .or. mask  < 1 ) Basin_Id = int(grid%nodata_value, i4)
       print*, trim(basinfName), ' ... read'
    else
       ! ascii 
       allocate( Basin_Id (grid%nrows , grid%ncols) )
       Basin_Id = int(grid%nodata_value,i4)
       !
       open(20, file=trim(basinfName), status = 'old')
       do i = 1, 6
          read(20, *) dummy
       end do
       !
       do i=1,grid%nrows
          read (20, *) (Basin_Id(i,j), j=1,grid%ncols)
       end do
       close(20)
    end if

    print*, 'Basin mask: '
    print*, trim(basinfName)
    print*, '... was read'
    if ( allocated (ZiT) ) deallocate (ZiT)
  end subroutine ReadDataMain

  !
  ! -----------------------------------------------------------------------
  !                     READ NetCDF DAILY/MONTHLY SM FIELDS
  ! -----------------------------------------------------------------------
  subroutine ReadFields_nc
    use mo_kind, only                       : i4, sp
    use imsl_libraries, only                : NDAYS, NDYIN
    implicit none

    character(256)                         :: dummy, fileNameSM
    integer(i4)                            :: d, dd, jd, js, je
    integer(i4)                            :: m
    integer(i4)                            :: y, yy
    integer(i4)                            :: tDays, mOld, mNew 
    real(sp), dimension(:), allocatable    :: V
    integer(i4), dimension(:), allocatable :: mDays
    ! netCDF
    integer(i4), parameter                 :: NDIMSc = 3   ! COSMO
    ! integer(i4), parameter               :: NDIMSw = 4   ! WRF-NOAH (original file, no dim reduction)
    integer(i4), parameter                 :: NDIMSw = 3   ! WRF-NOAH modified with nco
    integer(i4), dimension(NDIMSc)         :: startC, icountC
    integer(i4), dimension(NDIMSw)         :: startW, icountW
    character (len = *), parameter         :: LAT_NAME = 'lat'
    character (len = *), parameter         :: LON_NAME = 'lon'
    real(dp), dimension(:,:), allocatable  :: array
    integer(i4), dimension(5)              :: dl      ! for displaying dimensions
    integer(i4)                            :: iRec    ! reading record
    !
    allocate ( mDays(nMonths) )
    allocate ( V (nCells) )
    allocate ( lats  (grid%nCols, grid%nRows), lons(grid%nCols, grid%nRows) )
    allocate ( array (grid%nCols, grid%nRows) )
    !
    if (SMI_flag ==4 ) then
       m        = 1
       Z        = 0.0_dp
       mDays    = 0
       mOld     = 1 
       !
       FileNameSM = DataPathIn
       varNameSM  = 'dSMfr'
       dl = get_NcDim(FileNameSM, varNameSM)
       if (dl(1) /= grid%ncols .or. dl(2) /= grid%nrows) then
          print*, 'Dimension of variable ', varNameSM, ' not correct'
          stop
       end if
       !
       startC = (/1,1,1/)
       icountC = (/grid%ncols,grid%nrows,1/)
       !
       d  = 0
       do y = yStart, yEnd
          ! reading whole year
          js    = NDAYS (1,   1, y)
          je    = NDAYS (31,  12, y)
          jd    = js  
          tDays = je - js + 1
          do jd = js, je 
             d = d + 1
             call NDYIN (jd, dd, mNew, yy)
             if ( mNew /= mOld ) then
                ! reset counter        
                mOld = mNew
                m = m + 1
             end if
             mDays(m) = mDays(m) + 1
             ! >>> new
             startC(3) = d
             ! Open the file & pack transposed for consistency    
             call Get_NcVar(FileNameSM, varNameSM, array(:,:), startC, icountC)
             V =  pack ( transpose( array ) , (mask /= int(grid%nodata_value,i4) ) )
             ! accumulate for monthly values on the fly
             Z(:,m) = Z(:,m) + real(V,dp)
          end do
       end do
       ! estimate the montly SM means
       forall (m = 1 : nMonths ) Z(:,m) = Z(:,m) / real( mDays(m), dp)
!!$     print*,""
!!$     do m=1, nMonths
!!$       print*, m, mDays(m), Z(1,m)
!!$     end do
    else if(SMI_flag ==5 ) then
       Z          = 0.0_dp
       !
       FileNameSM = DataPathIn
       varNameSM  = 'SM_Lall'
       dl = get_NcDim(FileNameSM, varNameSM)
       if (dl(1) /= grid%ncols .or. dl(2) /= grid%nrows) then
          print*, 'Dimension of variable ', varNameSM, ' not correct'
          stop
       end if
       ! old counter for original sm file before dim reduction
       !   startW = (/1,1,1,1/)
       !   icountW = (/grid%ncols,grid%nrows,1,1/)
       ! new counter (3D)
       startW = (/1,1,1/)
       icountW = (/grid%ncols,grid%nrows,1/)
       !
       do m = 1, nMonths 
          ! reading every month
          !   startW(4) = m  
          startW(3) = m
          ! Open the file & pack transposed for consistency
          call Get_NcVar(FileNameSM, varNameSM, array(:,:), startW, icountW)
          V =  pack ( transpose( array ) , (mask /= int(grid%nodata_value,i4) ) )
          ! store monthly values
          Z(:,m) = real(V,dp)
       end do
    end if
    ! Read the latitude and longitude data.
    call Get_NcVar(FileNameSM, LAT_NAME, lats )
    call Get_NcVar(FileNameSM, LON_NAME, lons )
    deallocate(array, V)

  end subroutine ReadFields_nc
  !
  ! -----------------------------------------------------------------------
  !                         READ BINARY PACKED SMI
  ! -----------------------------------------------------------------------
  subroutine ReadSMI
    use mo_kind, only                    : i4, sp
    implicit none
    character(256)                      :: dummy, fileName
    integer(i4)                         :: m
    real(sp), dimension(:), allocatable :: V
    !
    ! Read Yearly Binary file
    fileName = trim(dataPathOut) //'mSMI.binP'
    ! open file
    open (unit=20, file=fileName, form='unformatted', access='direct', recl=4*nCells)
    allocate ( SMIp(nCells,nMonths) )
    allocate ( V(nCells) )
    !
    do m = 1, nMonths
      read  (20, rec=m) V
     ! store monthly values
      SMIp(:,m) = real(V,dp)
    end do
    close(20)
    deallocate ( V )
  end subroutine ReadSMI

  ! -----------------------------------------------------------------------
  !	                                 Header
  ! -----------------------------------------------------------------------
  subroutine ReadOneHeader(fileName,iu)
    implicit none
    character(len=*), intent(in)  :: fileName                             ! file name
    integer(i4),      intent(out) :: iu                                   ! input channel
    character(256)                :: dummy
    !
    iu = 10
    open (unit=iu, file=headerfName, status='old')
    read (iu, *) dummy, grid%ncols
    read (iu, *) dummy, grid%nrows
    read (iu, *) dummy, grid%xllcorner
    read (iu, *) dummy, grid%yllcorner
    read (iu, *) dummy, grid%cellsize
    read (iu, *) dummy, grid%nodata_value
    print *, 'Header of file: '
    print *, trim(fileName)
    print *, '... was read'
  end subroutine ReadOneHeader
  !
  ! -----------------------------------------------------------------------
  !	                             BINARY PACKED FILES
  ! -----------------------------------------------------------------------
  subroutine ReadFields_BinP
    use mo_kind,             only        : i4, sp
    use numerical_libraries, only        : NDAYS 
    implicit none
    character(256)                      :: dummy, fileName
    integer(i4)                         :: d, dd, jd, js, je
    integer(i4)                         :: m
    integer(i4)                         :: y, yy

    integer(i4)                         :: tDays, mOld, mNew 
    real(sp), dimension(:), allocatable :: V
    integer(i4), dimension(:), allocatable :: mDays
    !
    ! Read Yearly Binary file
    allocate ( V(nCells), mDays(nMonths) )
    m        = 1
    Z        = 0.0_dp
    mDays    = 0
    mOld     = 1 
    do y = yStart, yEnd
       ! open file
       write (dummy, '(i4,a5)') y, '.binP'
       fileName = trim(dataPathIn) // trim(dummy)
       open (unit=20, file=fileName, form='unformatted', access='direct', status='old', action='read', recl=4*nCells)
       ! reading whole year
       js    = NDAYS (1,   1, y)
       je    = NDAYS (31,  12, y)
       jd    = js  
       tDays = je - js + 1
       d     = 0
       do jd = js, je 
          d = d + 1
          call NDYIN (jd, dd, mNew, yy)
          if ( mNew /= mOld ) then
             ! reset counter        
             mOld = mNew
             m = m + 1
          end if
          mDays(m) = mDays(m) + 1
          read  (20, rec=d) V
          ! accumulate for monthly values on the fly
          Z(:,m) = Z(:,m) + real(V,dp)
       end do
       close(20)
    end do
    ! estimate the montly SM means
    forall (m = 1 : nMonths ) Z(:,m) = Z(:,m) / real( mDays(m), dp)

    ! do m=1, nMonths
    !   print*, m, mDays(m), Z(444,m)
    !  end do

    deallocate ( V, mDays )
  end subroutine ReadFields_BinP

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
  subroutine WriteNetCDF(wFlag, d) 
    use mo_kind, only              : i4, sp
    !
    use netCDF_varDef
    use netcdf
    implicit none
    !
    !implicit none
    integer(i4), intent (in), optional     :: d                ! optional, duration
    integer(i4), intent (in)               :: wFlag            !
    integer(i4), target                    :: m                ! netCDF counter
    integer(i4), save                      :: ncId             ! netCDF ID handler
    !
    real(sp), dimension(:,:), allocatable  :: Zu               ! field real unpacked 
    real(sp), dimension(:,:), allocatable, target:: ZuT        ! field real unpacked transposed
    integer(i4), dimension(:,:), allocatable  :: Ziu           ! field integer unpacked 
    integer(i4), dimension(:,:), allocatable, target:: ZiuT    ! field integer unpacked transposed

    integer(i4)                            :: iLoc(2)
    !
    select case (wFlag)
    case (1)
       ! save monthly soil moisture averages
       !
       ! set netCDF variables
       call set_netCDF_SM
       ! to create a new netCDF
       call create_netCDF(ncId) 
       ! write static variables  
       call write_static_netCDF(ncId)
       ! temporal storage
       allocate (Zu  (nLats , nLons))
       allocate (ZuT (nLons , nLats))
       ! save averages
       do m = 1, nMonths
          ! Unpack Z vector to their corresponding grids
          Zu = unpack ( real( Z(:,m), sp ), (mask /= int(grid%nodata_value,i4) ), real(grid%nodata_value,sp) )
          ZuT = transpose (Zu)
          ! put into the 4th variable
          V(4)%G2_f => ZuT 
          call write_dynamic_netCDF(ncId,m)
          !
          iLoc = 0
          iloc = maxloc(Zu, Zu > 1_dp)
          if (any(iLoc > 0)) then 
             write(6,'(a,2i5,a,i5)')  'WARNING: problem in cell (i,j): ',  iLoc, ' month:', m
          end if
       end do
       ! close file
       call close_netCDF(ncId) 

    case (2)
       ! save SMIp (same sequence as in 1)
       !
       call set_netCDF_SMI
       call create_netCDF(ncId)
       call write_static_netCDF(ncId)
       allocate (Zu  (nLats , nLons))
       allocate (ZuT (nLons , nLats))
       do m = 1, nMonths
          Zu = unpack ( real( SMIp(:,m), sp) , (mask /= int(grid%nodata_value,i4) ), real(grid%nodata_value,sp)  )
          ZuT = transpose (Zu)
          V(4)%G2_f => ZuT 
          call write_dynamic_netCDF(ncId,m)
       end do
       call close_netCDF(ncId)

    case (3)
       ! save drought mSMIc (same sequence as in 1)
       !
       call set_netCDF_mSMIc
       call create_netCDF(ncId)
       call write_static_netCDF(ncId)
       allocate (Ziu  (nLats , nLons))
       allocate (ZiuT (nLons , nLats))
       do m = 1, nMonths
          Ziu =  SMIc(:,:,m)
          ZiuT = transpose (Ziu)
          V(4)%G2_i => ZiuT 
          call write_dynamic_netCDF(ncId,m)
       end do
       call close_netCDF(ncId)

    case (4)
       ! save drought clusters (same sequence as in 1)
       !
       call set_netCDF_DC
       call create_netCDF(ncId)
       call write_static_netCDF(ncId)
       allocate (Ziu  (nLats , nLons))
       allocate (ZiuT (nLons , nLats))
       do m = 1, nMonths
          Ziu = idCluster(:,:,m)
          ZiuT = transpose (Ziu)
          V(4)%G2_i => ZiuT
          call write_dynamic_netCDF(ncId,m)
       end do
       call close_netCDF(ncId)

    case (5)
       ! save Severity for a given duration d (same sequence as in 1)
       !
       call set_netCDF_SEV(d)
       call create_netCDF(ncId)
       call write_static_netCDF(ncId)
       allocate (Zu  (nLats , nLons))
       allocate (ZuT (nLons , nLats))
       do m = 1, nDsteps
          Zu = severity(:,:,m)
          ZuT = transpose (Zu)
          V(4)%G2_f => ZuT 
          call write_dynamic_netCDF(ncId,m)
       end do
       call close_netCDF(ncId)

   case(6)
       ! save error file (cells with problems)
       call set_netCDF_eMask
       call create_netCDF(ncId)
       allocate (ZiuT (nLons , nLats))
       ZiuT = transpose( unpack ( eMask , (mask /= int(grid%nodata_value,i4) ), int(grid%nodata_value,i4) )   )
       V(5)%G2_i        => ZiuT
       call write_static_netCDF(ncId)
       call close_netCDF(ncId)

    end select

    if ( allocated (Zu) ) deallocate (Zu, ZuT)
    if ( allocated (Ziu) ) deallocate (Ziu, ZiuT)
  end subroutine WriteNetCDF

  !*************************************************************************
  !    PURPOSE    WRITE binP files
  !    FORMAT     binary packed for mHM
  !
  !    AUTHOR:    Luis E. Samaniego-Eguiguren, UFZ
  !    UPDATES
  !               Created        Sa   25.02.2011   
  !               Last Update    Sa     
  !**************************************************************************
  subroutine WritebinP(wFlag) 
    use mo_kind, only                       : i4, sp
    
    implicit none
    !
    integer(i4), intent (in)               :: wFlag            !
    integer(i4)                            :: m                ! 
    character(len=256)                     :: fName
    !
    select case (wFlag)
    case (1)
       ! save monthly soil moisture averages
       ! temporal storage
       fName = trim(dataPathOut) //'mean_m_SM.binP'
       open (unit=20, file=fName, form='unformatted', access='direct', recl=4*nCells)
       ! save averages
       do m = 1, nMonths
          write (20,rec=m) real( Z(:,m), sp ) 
      end do
      call writeHeader(fName, headerfName) 
      ! close file
      close(20)
    case (2)
       ! save SMIp
       !
       fName = trim(dataPathOut) //'mSMI.binP'
       open (unit=20, file=fName, form='unformatted', access='direct', recl=4*nCells)
       ! save averages
       do m = 1, nMonths
          write (20,rec=m) real( SMIp(:,m), sp)
      end do
      call writeHeader(fName, headerfName) 
      ! close file
      close(20)
    end select
    print *, 'binP file was created'
  end subroutine WritebinP

! *************************************************************************
! SUBROUTINE   WRITE HEADER binP
! *************************************************************************
  subroutine writeHeader(fName, fMask)
    implicit none
    character(256), intent(in)    :: fName, fMask
    character(256)                :: fNameH
    !
    fNameH = trim(fName)//'.header.txt'
    open (unit=100, file=fNameH, status='unknown')
    write (100, 1)  'ncols',        grid%ncols
    write (100, 1)  'nrows',        grid%nrows
    write (100, 2)  'xllcorner',    grid%xllcorner
    write (100, 2)  'yllcorner',    grid%yllcorner
    write (100, 2)  'cellsize',     grid%cellsize
    write (100, 2)  'NODATA_value', grid%nodata_value
    write (100, 3)  'mask_file',    trim(fMask)
    write (100, 1)  'nCells' ,      nCells
    write (100, 1)  'nRecords' ,    nMonths
    write (100, 1)  'recLenght' ,   4*nCells
    close(100)
    ! formats
1   format (a12, 2x, i10)
2   format (a12, 2x, f10.1)
3   format (a12, 2x, a)
  end subroutine writeHeader

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
subroutine WriteResultsBasins
  use mo_kind, only          : i4
  !
  implicit none
  !
  integer(i4)               :: i, m, k, y, j
  character(len=256)        :: fName
  !-------------------------
  ! basin wise
  !-------------------------
  allocate( Basin_SMI (nMonths, nBasins+1) )
  Basin_SMI = grid%nodata_value
  !
  do m = 1, nMonths
    do i = 1, nBasins
      !
      Basin_SMI(m,i) = sum( SMI(:,:,m), Basin_Id(:,:) == i ) / real(count(Basin_Id(:,:) == i), dp)
      !
    end do
    Basin_SMI(m,nBasins+1) = sum(SMI(:,:,m), mask /= int(grid%nodata_value,i4) ) / real(count(mask /= int(grid%nodata_value,i4) ), dp)
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
