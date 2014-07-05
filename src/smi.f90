!**********************************************************************
!  SOIL MOISTURE INDEX                                                *
!  PURPOSE     Estimate SMI based on monthly values of SM_mHM         *
!  CREATED     Luis Samaniego, 15.02.2011                             *
!**********************************************************************
module SMIndex
  !
  use mo_kind,          only: i4, dp
  use InputOutput,      only: SMIp
  use InputOutput,      only: WriteNetCDF
  !
  implicit none
  ! public routines
  public :: optimize_width
  public :: calSMI
  !
  private

contains
  ! subroutine for estimating SMI for first array
  subroutine optimize_width( opt_h, silverman_h, SM, tmask, nodata, offSet )
    use mo_kind,             only : i4, sp
    use mo_utils,            only : equal
    ! use numerical_libraries, only : DGCDF
    ! use GCDF_INT
    use InputOutput,         only : nMy

    ! use setGRG2,             only : OPTI, checkLimits
    use mo_kernel,           only : kernel_density_h, kernel_cumdensity
    ! input variables
    logical,                  intent(in) :: silverman_h ! optimize kernel width
    logical,  dimension(:,:), intent(in) :: tmask
    real(sp),                 intent(in) :: nodata
    real(sp), dimension(:,:), intent(in) :: SM
    real(dp),                 intent(in) :: offSet

    ! output variables
    real(dp), dimension(:,:), intent(inout) :: opt_h

    ! local variables
    integer(i4)                             :: nYears
    integer(i4)                             :: nObs
    integer(i4)                             :: ii     ! cell index
    integer(i4)                             :: mm     ! month index
    real(dp), dimension(:),     allocatable :: X
    !
    ! Initialization
    nObs   = size( SM, 2 )
    nYears = nObs / nMy
    allocate ( X(nYears) )
    X    = -9999._dp

    do ii = 1, size( SM, 1 )
       do mm = 1, nMy
          ! select values for month j (1...12)
          if ( count( tmask(:,mm) ) .eq. 0_i4 ) cycle
          X(:) = pack( SM( ii, :), tmask(:,mm) )
          !
          ! determine kernel width if these have not been read from file
          ! call OPTI ! optimize with imsl
          opt_h(ii,mm) = kernel_density_h( X(:), silverman = silverman_h )
       end do
    end do

    ! clean up
    deallocate( X )

  end subroutine optimize_width

  ! subroutine for calculating SMI for second array with pdf of first one
  subroutine calSMI( hh, sm_est, tmask_est, sm_eval, tmask_eval, SMI )
    
    use mo_kind,     only: i4, sp, dp
    use InputOutput, only: nMy
    use mo_kernel,   only: kernel_cumdensity
    
    implicit none

    ! input variables
    real(dp), dimension(:,:), intent(in)  :: hh
    real(sp), dimension(:,:), intent(in)  :: sm_est
    logical,  dimension(:,:), intent(in)  :: tmask_est
    real(sp), dimension(:,:), intent(in)  :: sm_eval
    logical,  dimension(:,:), intent(in)  :: tmask_eval
    ! output variables
    real(sp), dimension(:,:), intent(out) :: SMI

    ! local variables
    integer(i4)                                        :: mm  ! loop index
    integer(i4)                                        :: ii  ! cell index
    integer(i4)                                        :: nObs_est
    integer(i4)                                        :: nObs_eval
    real(dp), dimension(:), allocatable                :: cdf
    real(dp), dimension(:), allocatable                :: X_est
    real(dp), dimension(:), allocatable                :: X_eval

    ! evaluate cumulative density
    do ii = 1, size( SM_est, 1 )
       do mm = 1, nMy
          ! cycle if month not present
          if ( count( tmask_eval(:,mm) ) .eq. 0_i4 ) cycle
          allocate( X_est( count( tmask_est(:,mm) ) ), &
               X_eval( count( tmask_eval(:,mm) ) ), cdf( count( tmask_eval(:,mm) ) ) )
          X_est(:)  = pack( real(SM_est(  ii, : ),dp), tmask_est( :,mm) )
          X_eval(:) = pack( real(SM_eval( ii, : ),dp), tmask_eval(:,mm) )
          cdf       = kernel_cumdensity( x_est, hh(ii,mm), xout = x_eval )
          SMI(ii,:) = unpack( real(cdf,sp), tmask_eval(:,mm), SMI(ii,:) )
          deallocate( X_est, X_eval, cdf )
       end do
    end do
  end subroutine calSMI
end module SMIndex
