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
  public :: calSMI
  public :: evalSMI
  !
  private

contains
  ! subroutine for estimating SMI for first array
  subroutine calSMI( opt_h, read_opt_h, silverman_h, SM_est, mask, nodata, offSet )
    use mo_kind,             only : i4, sp
    use mo_utils,            only : equal
   ! use numerical_libraries, only : DGCDF
   ! use GCDF_INT
    use InputOutput,         only : nMy, nYears, nMonths, nCells, &
                                    SMI_flag, eMask, SMI
    use kernelSmoother,      only : nObs, nInter, evalPDF, edf, pdf
    ! use setGRG2,             only : OPTI, checkLimits
    use mo_kernel,           only : kernel_density_h
    ! input variables
    logical,                  intent(in) :: read_opt_h ! reading kernel width
    logical,                  intent(in) :: silverman_h ! optimize kernel width
    logical,  dimension(:,:), intent(in) :: mask
    real(sp),                 intent(in) :: nodata
    real(sp), dimension(:,:), intent(in) :: SM_est
    real(dp),                 intent(in) :: offSet

    ! output variables
    real(sp), dimension(:,:), intent(inout):: opt_h

    ! local variables
    integer(i4)                              :: i, j, k, m, iStatus
    integer(i4), dimension(1)                :: minpos
    integer(i4)                              :: nrows
    integer(i4)                              :: ncols
    real(dp),    dimension(:),   allocatable :: X                   ! sample
    !
    ! Initialization
    nrows = size( mask, 1 )
    ncols = size( mask, 2 )
    nObs = nYears                                      ! monthly PDF
    if ( .not. allocated(SMIp) ) allocate ( SMIp(nCells,nMonths), eMask(nCells) )
    allocate ( X(nObs), edf(0:nObs+1,2), pdf(nInter,3) )
    ! check limits
    iStatus = 0
    eMask   = 0
    do i= 1, nCells
       do j = 1, nMy
          ! select values for month j (1...12)
          X(:) = SM_est(i, j : nMonths : nMy)
          ! find optimal bandwidth for PDF of cell i in month j
          ! if ( checkLimits(i) == .false. ) then
          !    iStatus = iStatus + 1
          !    eMask(i) = 1
          !    Z(i, j : nMonths : nMy ) = grid%nodata_value  !here a posible check
          ! end if
       end do
    end do
    !
    if (iStatus > 0) then
       ! save error locations
       ! call WriteNetCDF(6,mask,nodata)
       print*, 'WARNING: Input mistake... '
       print*, 'At least one cell has boundary problems      '
       print*, 'SMI will not be estimated in those cells  ...'

       print*, 'Program must stop.'
       stop
    end if
    !
    SMIp = nodata
    do i= 1, nCells
       print *, i, nCells
       do j = 1, nMy
          ! select values for month j (1...12)
          X(:) = SM_est( i, j : nMonths : nMy)
          if ( any( X .eq. nodata) ) cycle
          !
          ! determine kernel width if these have not been read from file
          if ( .not. read_opt_h ) then
             ! call OPTI ! optimize with imsl
             opt_h(i,j) = kernel_density_h( X(:), silverman = silverman_h )
          end if
          ! evaluate pdf
          call evalPDF(X(:), real(opt_h(i,j),dp), offSet)
          ! call evalEDF
          ! call WriteResultsKernel(1,i,j,opt_h)
          ! call WriteResultsKernel(2,i,j,opt_h)
          ! call WriteResultsKernel(3,i,j,opt_h)
          ! print*, 'PDF of cell: ', i, ' month: ', j
          !
          ! find quantiles for each month using the estimated PDF
          do k = 1, nYears
             m = (k-1) * nMy + j
             !SMIp(i,m) = DGCDF( Z(i,m), 4, nInter, pdf(:,1), pdf(:,2))
             ! SMIp(i,m) = GCDF( Z(i,m), pdf(:,1), pdf(:,2), IOPT=4 )
             ! calculate argmin
             minpos = minloc( abs(pdf(:,1) - SM_est(i,m)) )
             SMIp(i,m) = pdf( minpos(1), 3) ! sum( pdf(:minpos(1),2) )
          end do
          !SMIp(i,:) = interpol( pdf(:,2), pdf(:,1), Z(i,:) )
       end do
       !print*, 'SMIp of cell: ', i
    end do
    !
    ! put SMIp to SMIp
    allocate( SMI( nrows, ncols, nMonths ) )
    do m = 1, nMonths 
       SMI (:,:,m) = unpack ( SMIp(:,m),  mask, real(nodata,dp) )
       !
       ! filter for possible values
       where (mask .and. SMI(:,:,m) .lt. 0.0_dp )
          SMI(:,:,m) = 0.5_dp
       end where
    end do
  end subroutine calSMI

  ! subroutine for calculating SMI for second array with pdf of first one
  subroutine evalSMI( sm_est, hh, sm_eval, SMI_eval )
    
    use mo_kind,     only: i4, sp, dp
    use InputOutput, only: nMy
    use mo_kernel,   only: kernel_cumdensity
    
    implicit none

    ! input variables
    real(sp), dimension(:,:), intent(in)  :: sm_est
    real(sp), dimension(:,:), intent(in)  :: hh
    real(sp), dimension(:,:), intent(in)  :: sm_eval
    ! output variables
    real(sp), dimension(:,:), intent(out) :: SMI_eval

    ! local variables
    integer(i4)                                        :: mm  ! loop index
    integer(i4)                                        :: ii  ! cell index
    integer(i4)                                        :: nObs_est
    integer(i4)                                        :: nObs_eval
    real(dp), dimension(:), allocatable                :: cdf
    real(dp), dimension(:), allocatable                :: X_est
    real(dp), dimension(:), allocatable                :: X_eval

    ! initialize
    nObs_est  = size(  SM_est, 2 )
    allocate( X_est(  nObs_est / nMy ) )
    nObs_eval = size( SM_eval, 2 ) 
    allocate( X_eval( nObs_eval / nMy ) )
    allocate( cdf(    nObs_eval / nMy ) )

    ! evaluate cumulative density
    do ii = 1, size( SM_est, 1 )
       print *, ii, size( SM_est, 1 )
       do mm = 1, nMy
          X_est(:)  = real(SM_est(  ii, mm : nObs_est : nMy ),dp)
          X_eval(:) = real(SM_eval( ii, mm : nObs_eval: nMy ),dp)
          cdf       = kernel_cumdensity( x_est, real(hh(ii,mm),dp), xout = x_eval )
          cdf       = cdf / sum( cdf )
          SMI_eval(ii, mm : nObs_eval : nMy ) = real( cdf, sp )
       end do
    end do
          
    stop 'TESTING'

  end subroutine evalSMI
end module SMIndex
