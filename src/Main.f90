!**********************************************************************************
!  SOIL MOISTURE DROUGHT INDEX
!  PURPOSE:  
!            Estimate the soil moisture index using daily values generated by mHM
!  AUTHOR:   
!            Luis E. Samaniego-Eguiguren, UFZ, 02-2011
!  ALGORITHM: 
!            1) Read netcdf files daily/montly
!            2) Estimate monthly values of mHM-SM(sum over all layers) for each grid cell
!               - write them into netCDF
!            3) Estimate the empirical density function for each grid cell
!               using a non-parametric approach,  e.g. kernel smoother whose
!               bandwith is estimated with an unbiased cross-validation criterium
!            4) Estimate the mHM-drought index
!                    q(s) = \int_a^s f(s)ds 
!                    where:
!                         s = total soil moisture in mm over all soil layers in mHM
!                         a = lower limit of soil moisture a a given grid 
!                       q(s)= quantile corresponding to s
!            5) Estimate the following indices
!               - start
!               - duration
!               - magnitud
!               - severity
!               - affected area
!            6) Save results in netCDF format
!      
!  UPDATES:
!          Created   Sa  15.02.2011          main structure
!                    Sa  20.02.2011          debuging
!                    Sa  25.05.2011          v3. scenarios
!                    Sa  02.04.2012          v4. read COSMO SM
!                    Sa  22.06.2012          v4. read WRF-NOAH SM
!                    Zi  07.11.2016          modularized version
!                    Sa  20.03.2017          daily SMI, SAD flag, restructuring SAD 
!**********************************************************************************
program SM_Drought_Index

  use mo_kind,               only : i4, sp, dp
  use InputOutput,           only : WriteNetCDF, &
                                    nDurations, durList, nClusters, &
                                    writeResultsCluster, WriteResultsBasins, &
                                    WriteSMI
  use mo_read,               only : ReadDataMain
  use SMIndex,               only : optimize_width, calSMI
  use mo_drought_evaluation, only : droughtIndicator, ClusterEvolution, ClusterStats, calSAD
  use mo_smi_constants,      only : nodata_sp


  implicit none

  ! variables
  logical                                    :: do_cluster  ! flag indicating whether cluster should be calculated
  logical                                    :: do_sad      ! flag indicating whether SAD analysis should be done
  logical                                    :: cluster_ext_smi  ! flag indicating to read external data for clustering 
  logical                                    :: eval_SMI    ! flag indicating whether SMI should be
  !                                                         ! calculated or read from file
  logical                                    :: read_opt_h  ! read kernel width from file
  logical                                    :: silverman_h ! flag indicating whether kernel width 
  !                                                         ! should be optimized
  logical                                    :: do_basin    ! do_basin flag
  ! logical,     dimension(:,:), allocatable   :: tmask_est   ! monthly mask for estimated SM
  ! logical,     dimension(:,:), allocatable   :: tmask_eval  ! monthly mask for evaluated SM
  logical,     dimension(:,:), allocatable   :: mask

  integer(i4)                                :: yStart
  integer(i4)                                :: yEnd
  integer(i4)                                :: mStart
  integer(i4)                                :: dStart
  integer(i4)                                :: nCells     ! number of effective cells
  integer(i4)                                :: d

  integer(i4)                                :: nCalendarStepsYear  !  Number of calendar time steps per year (month=12, day=365)

  integer(i4)                                :: thCellClus ! treshold  for cluster formation in space ~ 640 km2
  integer(i4)                                :: nCellInter ! number cells for joining clusters in time ~ 6400 km2
  integer(i4)                                :: deltaArea  ! number of cells per area interval
  integer(i4), dimension(:,:), allocatable   :: Basin_Id   ! IDs for basinwise drought analysis
  integer(i4), dimension(:),   allocatable   :: timepoints ! vector containing timesteps for output writing
  integer(i4), dimension(:,:), allocatable   :: cellCoor   ! 

  real(sp)                                   :: SMI_thld   ! SMI threshold for clustering
  real(sp)                                   :: cellsize   ! cell edge lenght of input data
  real(sp),    dimension(:,:), allocatable   :: SM_est     ! monthly fields packed for estimation
  real(sp),    dimension(:,:), allocatable   :: SM_eval    ! monthly fields packed for evaluation
  real(sp),    dimension(:,:), allocatable   :: SMI        ! soil moisture index at evaluation array
  integer(i4), dimension(:,:,:), allocatable :: SMIc       ! Drought indicator 1 - is under drought
  !                                                        !                   0 - no drought
  real(dp),    dimension(:,:), allocatable   :: opt_h      ! optimized kernel width field
  real(dp),    dimension(:,:), allocatable   :: lats, lons ! latitude and longitude fields of input

  
  ! file handling 
  character(256)                             :: outpath    ! output path for results


  call ReadDataMain( SMI, do_cluster, cluster_ext_smi, eval_SMI, read_opt_h, silverman_h, opt_h, lats, lons, do_basin,    &
       mask, SM_est, SM_eval,  yStart, yEnd, mStart, dStart, Basin_Id, timepoints, &
       SMI_thld, outpath, cellsize, thCellClus, nCellInter, do_sad, deltaArea, &
       nCalendarStepsYear ) ! tmask_eval,  tmask_est,
  
  ! initialize some variables
  nCells = count( mask )     ! number of effective cells
  
  print*, 'FINISHED READING'

  ! optimize kernel width
  if ( (.NOT. read_opt_h) .AND. (.NOT. cluster_ext_smi)) then
     call optimize_width( opt_h, silverman_h, SM_est,  nCalendarStepsYear, yStart, yEnd )  ! tmask_est,
     print *, 'optimizing kernel width...ok'
  end if

  ! evaluate SMI at second data set SMI_eval
  if (.NOT. cluster_ext_smi) then 
     if ( eval_SMI ) then
        allocate( SMI( size( SM_eval, 1 ), size( SM_eval, 2 ) ) )
        SMI(:,:) = nodata_sp
        call calSMI( opt_h, SM_est, SM_eval, nCalendarStepsYear, SMI, yStart, yEnd ) ! tmask_est
     else
        allocate( SMI(size(SM_est,1),size(SM_est,2)) )
        SMI(:,:) = nodata_sp
        call calSMI( opt_h, SM_est, SM_est, nCalendarStepsYear, SMI, yStart, yEnd ) ! tmask_est
     end if
     print *, 'calculating SMI... ok'
  end if
     
  ! write output
  if (.NOT. cluster_ext_smi) then
     if ( read_opt_h ) then
        call WriteSMI( outpath, SMI, mask, yStart, mStart, dStart,  yEnd, &
             timepoints, nCalendarStepsYear, lats, lons )
     else
        call WriteSMI( outpath, SMI, mask, yStart, mStart, dStart,  yEnd, &
             timepoints, nCalendarStepsYear, lats, lons, hh = opt_h )
     end if
     print *, 'write SMI...ok'
  end if
     
  ! calculate drought cluster
  if ( do_cluster ) then
     ! drought indicator 
     call droughtIndicator( SMI, mask, SMI_thld, cellCoor, SMIc )
     call WriteNetCDF(outpath, SMIc, 3, yStart, dStart, mStart, timepoints, lats, lons)
     
     ! cluster indentification
     call ClusterEvolution( SMIc,  size( mask, 1), size( mask, 2 ), size(SMI, 2), nCells, cellCoor, nCellInter, thCellClus)
     call WriteNetCDF(outpath, SMIc, 4, yStart, dStart, mStart, timepoints, lats, lons)

     ! statistics  
     call ClusterStats(SMI, mask, size( mask, 1), size( mask, 2 ), size(SMI, 2), nCells, SMI_thld )

     ! write results
     if (nClusters > 0) call writeResultsCluster(SMIc, outpath, 1, yStart, yEnd, size(SMI, 2), nCells, deltaArea, cellsize)
     print *, 'Cluster evolution ...ok'
  end if

  ! SAD analysis
  if ( do_sad ) then
     do d = 1, nDurations
        call calSAD(SMI, mask, d, size( mask, 1), size( mask, 2 ), size(SMI, 2), nCells, deltaArea, cellsize)
        ! write SAD for a given duration + percentiles
        call writeResultsCluster(SMIc, outpath, 2, yStart, yEnd, size(SMI, 2), nCells, deltaArea, cellsize, durList(d))
        call WriteNetCDF(outpath, SMIc, 5, yStart, dStart, mStart, timepoints, lats, lons, durList(d))
     end do
  end if
  
     
  ! make basin averages
  if ( do_basin ) then
     ! write SMI average over major basins
     print *, 'calculate Basin Results ...'
     call WriteResultsBasins( outpath, SMI, mask, yStart, yEnd, size( SM_est, 2 ), Basin_Id )
  end if

  print *, 'DONE!'
  !
  !! print statement for check_cases
  print *, 'program: finished!'
  !
end program SM_Drought_Index
