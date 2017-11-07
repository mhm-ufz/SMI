!**********************************************************************
!  SOIL MOISTURE INDEX                                                *
!  PURPOSE     Estimate SMI based on monthly values of SM_mHM         *
!  CREATED     Luis Samaniego, 15.02.2011                             *
!**********************************************************************
module SMIndex
  !
  use mo_kind,          only: dp
  !
  implicit none
  ! public routines
  public :: optimize_width
  public :: calSMI
  !
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

  
  ! ! create objective function of kernel_cumdensity and minimize it using
  ! ! nelmin because function is monotone
  ! subroutine invSMI(sm_est, hh, SMI_invert, SM_invert)
  !   use mo_kind, only: i4, sp, dp
  !   ! input variables
  !   real(dp), dimension(:,:), intent(in)  :: hh
  !   real(sp), dimension(:,:), intent(in)  :: sm_est
  !   real(sp), dimension(:,:), intent(in)  :: SMI_invert
  !   ! output variables
  !   real(sp), dimension(:,:), intent(out) :: SM_invert

  !   ! local variables
  !   integer(i4) :: n_cells, n_time
  !   integer(i4) :: ii, tt ! loop variables

  !   ! shape of arrays
  !   n_cells = size(SMI_invert, 1)
  !   n_time  = size(SMI_invert, 2)
    
  !   do ii = 1, n_cells
  !      do tt = 1, n_time

  !         ! call optimizer for this inversion

  !         ! X_est(:)  =  nodata_dp
  !         ! X_eval(:) =  nodata_dp

  !         ! X_est(:)  = real(SM_est ( ii, mm:size(SM_est,2):nCalendarStepsYear ),  dp)
  !         ! X_eval(:) = real(SM_eval( ii, mm:size(SM_est,2):nCalendarStepsYear ),  dp)
          
  !         ! !X_est(:)  = pack( real(SM_est(ii,:),  dp), tmask_est(:,mm) )
  !         ! !X_eval(:) = pack( real(SM_eval(ii,:), dp), tmask_eval(:,mm))

  !         ! cdf(:)    = kernel_cumdensity(x_est, hh(ii,mm), xout=x_eval)
  !         ! SMI(ii, mm:size(SM_est,2):nCalendarStepsYear) =  real(cdf(:), sp)
  !         ! ! deallocate( X_est, X_eval, cdf )
  !      end do
  !   end do
   
  ! end subroutine invSMI

  ! ! define the root function for optimizer
  ! function root_cdf(pp, 
  !         root_cdf = SMI_invert - kernel_cumdensity(xx_ext, hh, xout=pp)
  ! end function root_cdf
end module SMIndex
