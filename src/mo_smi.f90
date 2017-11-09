!**********************************************************************
!  SOIL MOISTURE INDEX                                                *
!  PURPOSE     Estimate SMI based on monthly values of SM_mHM         *
!  CREATED     Luis Samaniego, 15.02.2011                             *
!  HISTORY     Stephan Thober, 16.12.2014 - added evalSMI capability  *
!              Stephan Thober, 07.11.2017 - added inversion of CDFs   *
!**********************************************************************
module mo_smi
  !
  use mo_kind,          only: dp
  !
  implicit none
  ! public routines
  public :: optimize_width
  public :: calSMI
  public :: invSMI
  !
  ! global variables for inversion of cdf
  real(dp), public, dimension(:), allocatable :: xx_est
  real(dp), public                            :: hh_est
  real(dp), public                            :: y_val

  private

  
contains

  ! subroutine for estimating SMI for first array
  subroutine optimize_width( opt_h, silverman_h, SM, nCalendarStepsYear, yStart, yEnd) ! tmask,
    use mo_kind,          only : i4, sp
    use mo_kernel,        only : kernel_density_h
    use mo_smi_constants, only : YearMonths
    implicit none
    
    ! input variables
    logical,                  intent(in) :: silverman_h ! optimize kernel width
 !   logical,  dimension(:,:), intent(in) :: tmask
    real(sp), dimension(:,:), intent(in) :: SM
    integer(i4),              intent(in) :: nCalendarStepsYear
    integer(i4), intent(in)              :: yStart
    integer(i4), intent(in)              :: yEnd



    ! output variables
    real(dp), dimension(:,:), intent(inout) :: opt_h

    ! local variables
    integer(i4)                             :: ii     ! cell index
    integer(i4)                             :: mm     ! month index
    real(dp), dimension(:),     allocatable :: X
    integer(i4)                             :: nYears ! 

    nYears = yEnd - yStart + 1
    allocate ( X(nYears) )

    do ii = 1, size( SM, 1 )
       do mm = 1, nCalendarStepsYear                               !YearMonths    
          ! select values for month j (1...12)
          !if ( count( tmask(:,mm) ) .eq. 0_i4 ) cycle
          !allocate(X(count(tmask(:,mm))))
          !X(:) = pack( SM( ii, :), tmask(:,mm) )
          !
          X  = real(SM ( ii, mm:size(SM,2):nCalendarStepsYear ),  dp)
          ! determine kernel width if these have not been read from file
          ! call OPTI ! optimize with imsl
          opt_h(ii,mm) = kernel_density_h( X(:), silverman = silverman_h )
          !deallocate(X)
       end do
    end do
    deallocate(X)
    !! END
  end subroutine optimize_width


 
  !!======================================================
  !! ESTIMATE SMI 
  !!======================================================
  !! subroutine for calculating SMI for second array with pdf of first one
  subroutine calSMI( hh, sm_est, sm_eval,  nCalendarStepsYear, SMI, yStart, yEnd)  ! tmask_est,
    
    use mo_kind,          only: i4, sp, dp
    use mo_kernel,        only: kernel_cumdensity
    use mo_smi_constants, only: YearMonths, nodata_dp
    
    implicit none

    ! input variables
    real(dp), dimension(:,:), intent(in)  :: hh
    real(sp), dimension(:,:), intent(in)  :: sm_est
!    logical,  dimension(:,:), intent(in)  :: tmask_est
    real(sp), dimension(:,:), intent(in)  :: sm_eval
!    logical,  dimension(:,:), intent(in)  :: tmask_eval
    integer(i4),              intent(in)  :: nCalendarStepsYear

    integer(i4), intent(in)               :: yStart
    integer(i4), intent(in)               :: yEnd

  
    ! output variables
    real(sp), dimension(:,:), intent(out) :: SMI

    ! local variables
    integer(i4)                                        :: mm  ! loop index
    integer(i4)                                        :: ii  ! cell index
    real(dp), dimension(:), allocatable                :: cdf
    real(dp), dimension(:), allocatable                :: X_est
    real(dp), dimension(:), allocatable                :: X_eval
    integer(i4)                             :: nYears ! 


    nYears = yEnd - yStart + 1

    allocate ( X_est (nYears) )
    allocate ( X_eval(nYears) )
    allocate ( cdf   (nYears) )
    
    do ii = 1, size(SM_est,1)           ! cell loop
       do mm = 1,  nCalendarStepsYear   ! calendar time loop                                            !YearMonths

          !! cycle if month not present
          !if ( count( tmask_eval(:,mm) ) .eq. 0_i4 ) cycle

          ! allocate( X_est( count( tmask_est(:,mm)  ) ), &
          !          X_eval( count( tmask_eval(:,mm) ) ), &
          !             cdf( count( tmask_eval(:,mm) ) )   )

          X_est(:) =  nodata_dp
          X_eval(:) =  nodata_dp

          X_est(:)  = real(SM_est ( ii, mm:size(SM_est,2):nCalendarStepsYear ),  dp)
          X_eval(:) = real(SM_eval( ii, mm:size(SM_est,2):nCalendarStepsYear ),  dp)
          
          !X_est(:)  = pack( real(SM_est(ii,:),  dp), tmask_est(:,mm) )
          !X_eval(:) = pack( real(SM_eval(ii,:), dp), tmask_eval(:,mm))

          cdf(:)    = kernel_cumdensity(x_est, hh(ii,mm), xout=x_eval)
          SMI(ii, mm:size(SM_est,2):nCalendarStepsYear) =  real(cdf(:), sp)
          ! deallocate( X_est, X_eval, cdf )
       end do
    end do

    deallocate( X_est, X_eval, cdf )

  end subroutine calSMI

  
  ! create objective function of kernel_cumdensity and minimize it using
  ! nelmin because function is monotone
  subroutine invSMI(sm_est, hh, SMI_invert, nCalendarStepsYear, &
       SM_invert)
    use mo_kind,          only: i4, sp, dp
    use mo_smi_constants, only: nodata_sp, nodata_dp
    use mo_nelmin,        only: nelmin
    use mo_percentile,    only: percentile
    implicit none
    
    interface 
       function root_cdf(pp)
         use mo_kind, only: dp
         implicit none
         real(dp), intent (in) :: pp(:)
         real(dp) :: root_cdf
       end function root_cdf
    end interface

    ! input variables
    real(dp), dimension(:,:),              intent(in) :: hh
    real(sp), dimension(:,:),              intent(in) :: sm_est
    real(sp), dimension(:,:),              intent(in) :: SMI_invert
    integer(i4),                           intent(in) :: nCalendarStepsYear

    ! output variables
    real(sp), dimension(:,:), allocatable, intent(out) :: SM_invert

    ! local variables
    integer(i4)                         :: n_cells
    integer(i4)                         :: n_years_est, n_years_invert
    integer(i4)                         :: ii, yy, mm ! loop variables
    real(dp)                            :: xx_inv(1)
    real(dp)                            :: pstart(1)
    real(dp)                            :: funcbest
    real(dp), dimension(:), allocatable :: y_inv

    ! initialize extents
    n_cells        = size(SMI_invert, 1)
    n_years_est    = size(sm_est, 2)/nCalendarStepsYear
    n_years_invert = size(SMI_invert, 2)/nCalendarStepsYear
    
    ! initialize output array
    allocate(xx_est(n_years_est))
    allocate(y_inv(n_years_invert))
    allocate(SM_invert(n_cells, size(SMI_invert, 2)))

    y_inv     = nodata_dp
    SM_invert = nodata_sp
    xx_est    = nodata_dp
    
    print *, 'start inversion of CDF'
    print *, shape(SM_est)
    print *, shape(SMI_invert)
    print *, shape(hh)
    !$OMP parallel default(shared) &
    !$OMP private(mm, yy, xx_est, hh_est, y_inv, y_val, pstart, xx_inv)
    !$OMP do
    do ii = 1, n_cells
       if (modulo(ii, 1000) .eq. 0) print *, ii, n_cells
       do mm = 1, nCalendarStepsYear
          xx_est(:) = real(SM_est    ( ii, mm:size(sm_est, 2):nCalendarStepsYear    ),  dp)
          hh_est    = hh(ii, mm)
          y_inv(:)  = real(SMI_invert( ii, mm:size(SMI_invert, 2):nCalendarStepsYear),  dp)

          do yy = 1, n_years_invert
             y_val     = y_inv(yy)
             pstart(:) = percentile(xx_est, y_val * 100.)
             ! use nelmin for optimization because cdf is a monotonic function
             xx_inv = nelmin(root_cdf, pstart=pstart, funcmin=funcbest)
             SM_invert(ii, (yy - 1) * nCalendarStepsYear + mm) = xx_inv(1)
          end do
       end do
    end do
    !$OMP end do
    !$OMP end parallel

    ! free memory
    deallocate(y_inv, xx_est)
    
end subroutine invSMI

end module mo_smi
