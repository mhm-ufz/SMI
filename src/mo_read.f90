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
  subroutine ReadDataMain( do_cluster, eval_SMI, read_opt_h, silverman_h, opt_h, lats, lons, basin_flag,   &
                           mask, SM_est, tmask_est, SM_eval, tmask_eval, yStart, yEnd, mstart,  Basin_Id,  &
                           time_sp, offSet, SMI_thld, outpath  )   

    use mo_kind,          only: i4
    use mo_utils,         only: equal
    use kernelSmoother,   only: flagKernelType, nInter, flagCDF   
    use mo_string_utils,  only: DIVIDE_STRING
    use mo_ncread,        only: Get_NcVar, Get_NcDim, Get_NcVarAtt

    use mo_smi_constants, only: nodata_dp, YearMonths
 
    implicit none

    ! input / output Variables
    logical,                                  intent(out) :: do_cluster  ! do cluster calculation
    logical,                                  intent(out) :: eval_SMI    ! should SMI be calculated
    logical,                                  intent(out) :: read_opt_h  ! read kernel width
    logical,                                  intent(out) :: silverman_h ! optimize kernel width
    logical,                                  intent(out) :: basin_flag  ! basin flag
    logical,     dimension(:,:), allocatable, intent(out) :: mask        ! grid mask
    integer(i4),                              intent(out) :: yStart      ! starting year
    integer(i4),                              intent(out) :: yEnd        ! ending year
    integer(i4),                              intent(out) :: mStart      ! starting year
    integer(i4), dimension(:,:), allocatable, intent(out) :: Basin_Id    ! IDs for basinwise drought analysis
    real(sp),    dimension(:),   allocatable, intent(out) :: time_sp
    real(sp),    dimension(:,:), allocatable, intent(out) :: SM_est      ! monthly fields packed for estimation
    logical,     dimension(:,:), allocatable, intent(out) :: tmask_est   ! temporal mask of estimated arr
    real(sp),    dimension(:,:), allocatable, intent(out) :: SM_eval     ! monthly fields packed for evaluation
    logical,     dimension(:,:), allocatable, intent(out) :: tmask_eval  ! temporal mask of evaluated arr
    real(dp),    dimension(:,:), allocatable, intent(out) :: opt_h
    real(dp),    dimension(:,:), allocatable, intent(out) :: lats, lons  ! latitude and longitude fields of input
    real(dp),                                 intent(out) :: offSet      ! shift starting/ending x
    real(dp),                                 intent(out) :: SMI_thld    ! SMI threshold for clustering
    character(len=*),                         intent(out) :: outpath     ! ouutput path for results

    ! local Variables
    integer(i4)                                     :: ii
    integer(i4)                                     :: mm
    integer(i4)                                     :: datatype  ! datatype of attribute
    integer(i4)                                     :: nCells    ! number of effective cells

    ! directories, filenames, attributes
    character(256)                                  :: AttValues    ! netcdf attribute values
    character(256)                                  :: maskfName
    character(256)                                  :: opt_h_file
    character(256)                                  :: SM_eval_file
    character(256)                                  :: basinfName
    
    ! variable names in netcdf input files
    character(256)                                  :: soilmoist_file
    character(256)                                  :: mask_vname
    character(256)                                  :: SM_vname
    character(256)                                  :: SM_eval_vname
    character(256)                                  :: basin_vname
    character(256)                                  :: opt_h_vname
    character(256)                                  :: time_units
    character(256), dimension(:),       allocatable :: strArr       ! dummy for netcdf attribute handling
    real(sp)                                        :: nodata_value ! local data nodata value in nc for mask creation

    real(sp),       dimension(:,:),     allocatable :: dummy_D2_sp
    real(sp),       dimension(:,:,:),   allocatable :: dummy_D3_sp
    real(dp),       dimension(:,:,:),   allocatable :: dummy_D3_dp

    ! read main config
    namelist/mainconfig/basin_flag, basinfName, basin_vname, maskfName, mask_vname, soilmoist_file, SM_vname,  &
                        yStart, yEnd, outpath, eval_SMI, silverman_h, read_opt_h, opt_h_vname, opt_h_file,     &
                        SM_eval_file, SM_eval_vname, flagKerneltype, nInter, offSet, flagCDF, SMI_thld, do_cluster

    ! read namelist
    open (unit=10, file='main.dat', status='old')
      read(10, nml=mainconfig)
    close (10)
    print*, 'main.dat read ...ok'




    ! read main mask
    call Get_ncVar( maskfName, trim(mask_vname), dummy_D2_sp )
    allocate( mask( size(dummy_D2_sp,1), size(dummy_D2_sp,2) ) )
    ! determine no data value
    call Get_NcVarAtt(maskfName, trim(mask_vname), 'missing_value', AttValues, dtype=datatype)
    ! convert to number 
    read(AttValues, *) nodata_value 
    ! create mask
    mask = merge( .true., .false., (dummy_D2_sp > nodata_value ) )
    deallocate( dummy_D2_sp )
    print*, 'mask read ...ok'



    ! read basin mask
    if ( basin_flag ) then
       call Get_ncVar( basinfName, trim(basin_vname), Basin_Id )
       ! determine no data value
       call Get_NcVarAtt(maskfName, trim(mask_vname), 'missing_value', AttValues, dtype=datatype)
       ! convert to number
       read(AttValues, *) nodata_value 
       ! consistency check
       if ( ( size( Basin_Id, 1) .ne. size( mask, 1 ) ) .or. &
            ( size( Basin_Id, 2) .ne. size( mask, 2 ) ) ) then
          print *, '***ERROR: size mismatch between basin field and given mask file'
          stop
       end if
       ! intersect mask and basin_id
       mask     = merge( .true., .false., mask .and. ( Basin_Id .gt. nodata_value ) )
       Basin_Id = merge( Basin_Id, int(nodata_value, i4), mask )
    end if

    nCells   = count( mask )
    ! consistency check
    if ( nCells .eq. 0 ) then
       print *, '***ERROR: no cell selected in mask'
       stop
    end if
    print*, 'basin ID read ...ok'
    
    ! read SM field
    call Get_ncVar( soilmoist_file, trim(SM_vname), dummy_D3_dp )
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
    print *, '***CAUTION: time axis assumed to be starting with 0!'
    call Get_NcVarAtt( soilmoist_file, 'time', 'units', time_units )
    call DIVIDE_STRING(trim(time_units), ' ', strArr)
    call DIVIDE_STRING(trim(strArr(3)), '-', strArr) 
    read( strArr(1), * ) yStart
    read( strArr(2), * ) mStart
    allocate( time_sp( size( SM_est, 2 ) ) )
    forall( mm = 1 : size( time_sp, 1 ) ) time_sp(mm) = real(mm-1,sp)
    allocate ( tmask_est( size(SM_est, 2), YearMonths ) )
    tmask_est = .false.
    do mm = 1, YearMonths
       tmask_est(:,mm) = ( mod( int(time_sp,i4) + mStart, YearMonths ) .eq. mod( mm, YearMonths ) )
    end do
    if ( any( count( tmask_est, dim = 1 ) .eq. 0_i4 ) ) &
         stop '***ERROR no data for estimation given for all calendar months, check time axis'
    
    ! read lats and lon from file
    call Get_ncVar( soilmoist_file, 'lat', lats )
    call Get_ncVar( soilmoist_file, 'lon', lons )

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
       allocate ( tmask_eval( size(SM_eval, 2), YearMonths ) )
       tmask_eval = .false.
       do mm = 1, YearMonths
          tmask_eval(:,mm) = ( mod( int(time_sp,i4) + mStart, YearMonths ) .eq. mod( mm, YearMonths ) )
       end do
       if ( all( count( tmask_eval, dim = 1 ) .eq. 0_i4 ) ) &
            stop '***ERROR no data in eval given, check time axis'
       print*, 'read soil moisture field for evaluation... ok'
    end if

    ! initialize opt_h
    allocate ( opt_h( ncells, YearMonths ) )
    opt_h = nodata_dp
    !
    if ( read_opt_h ) then
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

    print *, '***WARNING: yStart and yEnd are only considered for output writing'
    
  end subroutine ReadDataMain

end module mo_read
