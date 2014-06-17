module kernelSmoother
  use mo_kind, only                             : i4,dp
  implicit none
  ! parameters                                                         
  integer(i4)                                  :: flagKernelType      ! type 1...9
  logical                                      :: flagCDF             ! print CDF or PDF
  integer(i4)                                  :: nObs                ! sample size
  integer(i4), parameter                       :: nPar = 1            ! number of parameters
  integer(i4)                                  :: nInter              ! number of interval for the eval PDF
  real(dp), dimension(:), allocatable          :: X                   ! sample
  real(dp), dimension(:,:), allocatable        :: pdf                 ! x, frequency
  real(dp), dimension(:,:), allocatable        :: edf                 ! x, empirical cummulative distribution
  real(dp)                                     :: hOpt                ! smoothing parameter or bandwith
  real(dp)                                     :: xMax                ! sample max val
  real(dp)                                     :: xMin                ! sample min val
  real(dp)                                     :: hInt                ! integration interval
  real(dp)                                     :: offSet              ! shift starting/ending x
  real(dp), dimension(:,:), allocatable        :: hOptDB              ! data bank of opti parameter  
!
contains
  !
  !***********************************************************
  !  Kernel functions (see Wilks pg. 36)
  !  Luis Samaniego, 02.09.2011
  !***********************************************************
  function fKernel(h, x0, iOut)
    use mo_kind, only                             : i4, sp, dp
    implicit none
    real(dp)                 :: fKernel           ! f(x_o)
    real(dp), intent(in)     :: h                 ! given smoothing parameter
    real(dp), intent(in)     :: x0                ! given value
    integer(sp), intent(in)  :: iOut              ! exclude observation for crossvalidation
    !
    real(dp)              :: u, hInv
    !
    real(dp), parameter   :: c0 = 0.50000_dp   ! 1/2  
    real(dp), parameter   :: c1 = 0.75000_dp   ! 3/4
    real(dp), parameter   :: c2 = 0.93750_dp   ! 15/16
    real(dp), parameter   :: c3 = 1.09375_dp   ! 35/32
    !
    real(dp), parameter   :: ip   = 0.3183098861837906_dp ! 1/pi
    real(dp), parameter   :: p4   = 0.7853981633974483_dp ! pi/4 
    real(dp), parameter   :: p2   = 1.5707963267948966_dp ! pi/2
    real(dp), parameter   :: iS2p = 0.3989422804014326_dp ! 1/sqrt{ 2 pi }
    !  
    integer(i4)           :: i
    !  
    hInv    = 1.0_dp / h
    fKernel = 0.0_dp 
    select case (flagKernelType)
    case (1)
       ! uniform
       do i = 1, nObs
          if (i==iOut) cycle
          u = (x0 - X(i) ) * hInv
          if (dabs(u) > 1.0_dp) cycle
          fKernel = fKernel + c0
       end do
    case (2)  
       ! triangular	  
       do i = 1, nObs
          if (i==iOut) cycle
          u = (x0 - X(i) ) * hInv
          if (dabs(u) > 1.0_dp) cycle
          fKernel = fKernel + (1.0_dp - u)
       end do
    case (3)  
       ! quadratic (Epanechnikov)
       do i = 1, nObs
          if (i==iOut) cycle
          u = (x0 - X(i) ) * hInv
          if (dabs(u) > 1.0_dp) cycle
          fKernel = fKernel + c1*(1.0_dp - u*u)
       end do
    case (4)
       ! quartic (biweight)   
       do i = 1, nObs
          if (i==iOut) cycle
          u = (x0 - X(i) ) * hInv
          if (dabs(u) > 1.0_dp) cycle
          fKernel = fKernel + c2*(1.0_dp - u*u)**2
       end do
    case (5)  
       ! Triweight (tricube)
       do i = 1, nObs
          if (i==iOut) cycle
          u = (x0 - X(i) ) * hInv
          if (dabs(u) > 1.0_dp) cycle
          fKernel = fKernel + c3*(1.0_dp - u*u)**3
       end do
    case (6)  
       ! Cosine
       do i = 1, nObs
          if (i==iOut) cycle
          u = (x0 - X(i) ) * hInv
          if (dabs(u) > 1.0_dp) cycle
          fKernel = fKernel + p4*dcos(p2*u)
       end do
    case (7)  
       ! Gaussian
       do i = 1, nObs
          if (i==iOut) cycle
          u = (x0 - X(i) ) * hInv
          fKernel = fKernel + iS2p*dexp(-c0*u*u)
       end do
    case (8) 
       ! Cauchy
       do i = 1, nObs
          if (i==iOut) cycle
          u = (x0 - X(i) ) * hInv
          fKernel = fKernel + ip / (1.0_dp + u*u)
       end do
    case (9) 
       ! Picard
       do i = 1, nObs
          if (i==iOut) cycle
          u = (x0 - X(i) ) * hInv
          fKernel = fKernel + c0 * dexp( -dabs(u) )
       end do
    end select
    fKernel = fKernel / real(nObs,dp) / h
  end function fKernel

  !***********************************************************
  !  Kernel functions (see Wilks pg. 36)
  !  Luis Samaniego
  !  02.09.2011
  !***********************************************************
  subroutine evalPDF(h)
    use mo_kind, only        : i4, dp
    implicit none
    integer(i4)             :: i, ioacc
    integer(i4), parameter  :: iDigit = 2

    real(dp), intent(in)  :: h                                ! given smoothing parameter
    ! set nice boundaries
    ioacc = 10**iDigit
    xMin = real ( floor( (minval(X) - offSet) * real(ioacc, dp), dp ), dp   ) / real(ioacc, dp)
    xMax = real ( ceiling( (maxval(X) + offSet) * real(ioacc, dp), dp ), dp ) / real(ioacc, dp)
    hInt = (xMax - xMin) / real(nInter-1,dp)
    pdf  = 0.0_dp
    do i = 1, nInter
       ! set equally spaced abscissas
       pdf(i,1) = xMin + real((i-1),dp)*hInt
       ! evaluate the PDF
       pdf(i,2) = fKernel(h, pdf(i,1), -9) / hInt
       ! evaluate the CDF
       if (i > 1) pdf(i,3) = pdf(i-1,3) + pdf(i,2) * hInt 
    end do
    if ( pdf(1,3) > 0.0_dp .or. pdf(nInter,3) < 1.0_dp ) then
       print*, 'WARNING: CDF boundaries not reached'
       print*, 'Increase offSet value...'
    end if
    ! normalize PDF
    pdf(:,2) = pdf(:,2) / sum(pdf(:,2))
    ! normalize CDF
    pdf(:,3) = pdf(:,3) / pdf(nInter,3)
    ! print*,  pdf(nInter,3)
  end subroutine evalPDF

  !***********************************************************
  !  Empirical Distribution Function
  !  Luis Samaniego, 02.09.2011
  !
  !  NOTE: should be called after evalPDF
  !***********************************************************
  subroutine evalEDF
    use mo_kind, only                    : i4,dp
    use numerical_libraries, only        : DSVRGN
    implicit none
    integer(i4)                         :: i
    real(dp)                            :: freq      ! frequency
    !
    edf(1:nObs,2) = 1.0_dp / real(nObs,dp)
    ! sort x values
    call DSVRGN(nObs,X,edf(1:nObs,1))
    edf(0,1)      = pdf(1,1)
    edf(nObs+1,1) = pdf(nInter,1)
    edf(0,2)      = 0.0_dp
    edf(nObs+1,2) = 1.0_dp
    do i=1,nObs
       edf(i,2) = edf(i,2) + edf(i-1,2)
    end do
  end subroutine evalEDF

end module kernelSmoother

