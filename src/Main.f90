!**********************************************************************************
!  SOIL MOISTURE DROUGHT INDEX
!  PURPOSE:  
!            Estimate the soil moisture index using daily values generated by mHM
!  AUTHOR:   
!            Luis E. Samaniego-Eguiguren, UFZ, 02-2011
!  ALGORITHM: 
!            1) Read binp files daily
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
!          Created   Sa  15.02.2011            main structure
!                    Sa  20.02.2011            debuging
!                    Sa  25.05.2011            v3. scenarios
!                    Sa  02.04.2012            v4. read COSMO SM
!                    Sa  22.06.2012            v4. read WRF-NOAH SM 
!**********************************************************************************
program SM_Drought_Index

  use mo_kind,               only : i4, sp, dp
  use InputOutput,           only : WriteNetCDF, &
                                    nDurations, durList, &
                                    writeResultsCluster, WriteResultsBasins, &
                                    WriteSMI
  use mo_read,               only : ReadDataMain
  use SMIndex,               only : optimize_width, calSMI
  use mo_drought_evaluation, only: droughtIndicator, ClusterEvolution, ClusterStats, calSAD
  use mo_smi_constants,      only: nodata_sp


  implicit none

  ! variables
  logical                                  :: do_cluster  ! flag indicating whether cluster should
  !                                                       ! be calculated
  logical                                  :: eval_SMI    ! flag indicating whether SMI should be
  !                                                       ! calculated or read from file
  logical                                  :: read_opt_h  ! read kernel width from file
  logical                                  :: silverman_h ! flag indicating whether kernel width 
  !                                                       ! should be optimized
  logical                                  :: do_basin    ! do_basin flag
  logical,     dimension(:,:), allocatable :: tmask_est   ! monthly mask for estimated SM
  logical,     dimension(:,:), allocatable :: tmask_eval  ! monthly mask for evaluated SM
  logical,     dimension(:,:), allocatable :: mask

  integer(i4)                              :: yStart
  integer(i4)                              :: yEnd
  integer(i4)                              :: mStart
  integer(i4)                              :: nMonths    ! number of simulated months
  integer(i4)                              :: nCells     ! number of effective cells
  integer(i4)                              :: d
  integer(i4), dimension(:,:), allocatable :: Basin_Id   ! IDs for basinwise drought analysis

  real(sp),    dimension(:,:), allocatable :: SM_est     ! monthly fields packed for estimation
  real(sp),    dimension(:,:), allocatable :: SM_eval    ! monthly fields packed for evaluation
  real(sp),    dimension(:),   allocatable :: time
  real(sp),    dimension(:,:), allocatable :: SMI        ! soil moisture index at evaluation array

  real(dp)                                 :: SMI_thld   ! SMI threshold for clustering
  real(dp),    dimension(:,:), allocatable :: opt_h      ! optimized kernel width field
  real(dp),    dimension(:,:), allocatable :: lats, lons ! latitude and longitude fields of input
  
  ! file handling 
  character(256)                           :: outpath    ! ouutput path for results


  call ReadDataMain( do_cluster, eval_SMI, read_opt_h, silverman_h, opt_h, lats, lons, do_basin,    &
       mask, SM_est, tmask_est, SM_eval, tmask_eval, yStart, yEnd, mStart, Basin_Id, time,          &
       SMI_thld, outpath)
  !
  ! initialize some variables
  nMonths  = size( SM_est, 2 ) ! number of months in dataset
  nCells   = count( mask )     ! number of effective cells
  !
  print*, 'FINISHED READING'

  ! optimize kernel width
  if ( .not. read_opt_h ) then
     call optimize_width( opt_h, silverman_h, SM_est, tmask_est )
     print *, 'optimizing kernel width...ok'
  end if

  ! evaluate SMI at second data set SMI_eval
  if ( eval_SMI ) then
     allocate( SMI( size( SM_eval, 1 ), size( SM_eval, 2 ) ) )
     SMI = nodata_sp
     call calSMI( opt_h, SM_est, tmask_est, SM_eval, tmask_eval, SMI )
  else
     allocate( SMI( size( SM_est, 1 ), size( SM_est, 2 ) ) )
     SMI = nodata_sp
     call calSMI( opt_h, SM_est, tmask_est,  SM_est,  tmask_est, SMI )
  end if
  print *, 'calculating SMI... ok'

  ! write output
  if ( read_opt_h ) then
     call WriteSMI( outpath, SMI, mask, yStart, mStart, time, lats, lons )
  else
     call WriteSMI( outpath, SMI, mask, yStart, mStart, time, lats, lons, hh = opt_h )
  end if
  print *, 'write SMI...ok'

  ! calculate drought cluster
  if ( do_cluster ) then
     ! drought indicator 
     call droughtIndicator( mask, nMonths, SMI_thld )
     call WriteNetCDF(outpath, 3, opt_h, SM_est, mask, yStart, lats, lons)
     
     ! cluster indentification
     call ClusterEvolution( size( mask, 1), size( mask, 2 ), nMonths, nCells)
     call WriteNetCDF(outpath, 4, opt_h, SM_est, mask, yStart, lats, lons)
     ! statistics  
     call ClusterStats( size( mask, 1), size( mask, 2 ), nMonths, nCells, Basin_Id, SMI_thld )
     !
     ! SAD analysis
     do d = 1, nDurations
        call calSAD(d, size( mask, 1), size( mask, 2 ), nMonths, nCells)
        ! write SAD for a given duration + percentiles
        call writeResultsCluster(outpath, 2, yStart, yEnd, nMonths, nCells, durList(d) )
        call WriteNetCDF(outpath, 5, opt_h, SM_est, mask, yStart, lats, lons, durList(d))
     end do
     ! write results
     call writeResultsCluster(outpath, 1, yStart, yEnd, nMonths, nCells)
     print *, 'Cluster evolution ...ok'
  end if

  ! make basin averages
  if ( do_basin ) then
     ! write SMI average over major basins
     print *, 'calculate Basin Results ...'
     call WriteResultsBasins( outpath, mask, yStart, yEnd, nMonths, Basin_Id )
  end if

  print *, 'DONE!'
  !
end program SM_Drought_Index
