!*****************************************************************************
!    Setting GRG2
!    AUTHOR
!         Luis E. Samaniego-Eguiguren, IREUS, 28.05.1999
!    DESCRIPTION
!         Initialization of the Nonlinear Optimization Subroutine GRG2 
!         The function to be optimized is suplied in subroutine GCOMP
!    DEFINITION OF INPUT/OUTPUT-CHANNELS
!         SCREEN OUTPUT:     *
!         OUTPUT FILE:       6
!******************************************************************************
module setGRG2
  use mo_kind,             only : i4, dp
  use kernelSmoother,      only : hOpt, nPar, X, nObs
  use numerical_libraries, only : DUVSTA, DEQTIL, DSVRGN
 ! use UVSTA_INT
 ! use EQTIL_INT
  ! inputs for GRG2
  IMPLICIT  DOUBLE PRECISION(A-H,O-Z), INTEGER(I,J,L,M,N)
  !
  ! Set problem size and objective
  !
  integer(4) , parameter :: NCORE = 5000          !DIMENSION OF THE {Z} ARRAY   
  integer(4) , parameter :: NNVARS = nPar         !NUMBER OF VARIABLES
  integer(4) , parameter :: NFUN = 1              !NUMBER OF FUNCTIONS INCLUDING OBJECTIVE
  integer(4) , parameter :: MAXBAS = NFUN         !UPPER LIMIT ON NUMBER OF BINDINGCONSTRAINTS
  integer(4) , parameter :: MAXHES = nPar         !MAXIMUM ALLOWABLE SIZE OF HESSIAN
  integer(4) , parameter :: NNOBJ  = 1            !INDEX OF COMPONENT OF VECTOR {G} CORR. TO OBJ. FUNCTION
  logical, parameter     :: INPRNT = .FALSE.	  !DO NOT PRINT ANY ECHO BACK OF INPUT DATA	 
  logical, parameter     :: OTPRNT = .FALSE.      !DO NOT PRINT ANY FINAL RESULTS	 
  ! setting output variables  	 
  real(8), dimension(NNVARS)   :: BLVAR, BUVAR
  real(8), dimension(NFUN)     :: BLCON, BUCON	
  real(8), dimension(NNVARS)   :: XX
  real(8), dimension(NFUN)     :: FCNS  		 !
  integer(4), dimension(NFUN)  :: INBIND
  real(8), dimension(NFUN)     :: RMULTS
  integer(4), dimension(NNVARS):: NONBAS
  real(8), dimension(NNVARS)   :: REDGR
  ! working space
  real(8), dimension(NCORE)    :: Z
  ! blanks
  real(8), dimension(NFUN)     :: RAMCON
  real(8), dimension(NNVARS)   :: RAMVAR
  ! arguments
  real(8), dimension(19)       :: DEFAUL, TITLE
  ! local variables
  real(dp), parameter                          ::  BIG = 1.0e31_dp
  real(dp), dimension(15,1)                    ::  stat
  integer(i4)                                  ::  nMiss
  integer(i4), parameter                       ::  nQprop = 2
  real(dp), parameter                          ::  c1 = 0.9_dp
  real(dp), parameter                          ::  c2 = 0.666666666666667_dp
  real(dp), parameter                          ::  c3 = 0.2_dp
  real(dp), dimension(nqProp)                  ::  Xemp
  real(dp), dimension(nQProp)                  ::  Xlo
  real(dp), dimension(nqProp)                  ::  Xhi
  real(dp), parameter, dimension(nqProp)       ::  QProp =(/0.25_dp, 0.75_dp/)
  real(dp)                                     ::  IQR
!
contains
!
subroutine OPTI
  ! Initialization of GRG2
  !
  ! Minimize the Obj. Function
  !
  DEFAUL(19) = 1.0d0
  !
  do i=1,18
    DEFAUL(I) = 1.0d0
  end do
  !
  ! Set tolerances
  !
	DEFAUL(1) = -1.0
	FPNEWT    =  1.0d-7
	DEFAUL(2) = -1.0
	FPINIT    =  1.0d-6
	DEFAUL(3) = -1.0
	FPSTOP    =  1.0d-6
	DEFAUL(6) = -1.0
	NNSTOP    = 10
  !
  ! Quadratic extrapolation
  !
  DEFAUL(15) = -1.0d0
  IIQUAD     = 1
  !
  ! Central difference approximation
  !
  DEFAUL(15) = -1.0d0
  LDERIV     = 1
  !
  ! Estimate statistics
  !
  call DUVSTA (0, nObs, 1, X, nObs, 0, 0, 1, 0.d0, 0.d0, 0, stat, 15, nMiss)
  call DEQTIL (nObs, X, nQprop, Qprop, Xemp, Xlo, Xhi, nMiss)
 ! call UVSTA (X, stat, MOPT=1, CONPRM=0.0_dp, CONPRV=0.0_dp)
 ! call EQTIL (X, nQprop, Qprop, Xemp, Xlo, Xhi, nMiss)
  IQR = Xemp(2) - Xemp(1)
  !
  ! Set initial values (Wilks, pg.38)
  XX(1) = min(c1*stat(3,1), c2*IQR) /dfloat(nObs)**c3
  !
  ! Set variable bounds
  !
  BLVAR(1) = (maxval(X) - minval(X))/dfloat(nObs-1)  ! average distance between observations
  BUVAR(1) = IQR
  !
  ! Set fuctions bounds
  !
  BLCON(1) = -BIG
  BUCON(1) = BIG
  !
  ! Other stuff
  !
  ! to write report USE:  OPEN(UNIT = 6,FILE = 'Report_OPT.sol') and do not CLOSE 6
  OPEN(UNIT = 6, FILE = 'Report_OPT.sol')
  !
  ! Call Nonlinear Optimization subroutine
  !
  CALL GRGSUB(INPRNT,OTPRNT,NCORE,NNVARS,NFUN,MAXBAS,&
      MAXHES,NNOBJ,TTITLE,BLVAR,BUVAR,BLCON,BUCON,DEFAUL,FPNEWT,FPINIT,&
      FPSTOP,FPSPIV,PPH1EP,NNSTOP,IITLIM,LLMSER,IIPR,IIPN4,IIPN5,&
      IIPN6,IIPER,IIDUMP,IIQUAD,LDERIV,MMODCG,&
      RAMCON,RAMVAR,XX,FCNS,INBIND,RMULTS,NONBAS,REDGR,&
      NBIND,NNONB,INFORM,Z)
  !
  ! delete report (partial results)
  !
    CLOSE (6, STATUS='DELETE')
  !
  ! keep parameters
  hOpt = real(XX(1),dp)
end subroutine OPTI

function checkLimits(iCell)
  use mo_kind,             only : i4,dp
  use kernelSmoother,      only : X, nObs
  !use numerical_libraries, only : DUVSTA, DEQTIL, DSVRGN
  use InputOutput,         only : DataPathOut

  !use UVSTA_INT
  !use EQTIL_INT
  implicit none
  integer(i4), intent(in) :: iCell
  real(dp), dimension(15,1)                    ::  stat
  integer(i4)                                  ::  nMiss
  integer(i4), parameter                       ::  nQprop = 2
  real(dp), parameter                          ::  c1 = 0.9_dp
  real(dp), parameter                          ::  c2 = 0.666666666666667_dp
  real(dp), parameter                          ::  c3 = 0.2_dp
  real(dp), dimension(nqProp)                  ::  Xemp
  real(dp), dimension(nQProp)                  ::  Xlo
  real(dp), dimension(nqProp)                  ::  Xhi
  real(dp), parameter, dimension(nqProp)       ::  QProp =(/0.25_dp, 0.75_dp/)
  real(dp)                                     ::  IQR
  real(dp)                                     ::  BLVAR1, BUVAR1
  logical                                      ::  checkLimits
  ! Estimate statistics
  !
  call DUVSTA (0, nObs,  1,    X, nObs, 0,    0,   1,    0.d0, 0.d0,     0, stat, 15, nMiss)
 ! call  UVSTA (X, stat, MOPT=1, CONPRM = 0.0_dp, CONPRV = 0.0_dp )
  call DEQTIL (nObs, X, nQprop, Qprop, Xemp, Xlo, Xhi, nMiss)
 ! call  EQTIL (X, nQprop, Qprop, Xemp, Xlo, Xhi, nMiss)
  IQR = Xemp(2) - Xemp(1)
  !
  ! Set variable bounds
  !
  BLVAR1 = (maxval(X) - minval(X)) / real((nObs-1),dp)  ! average distance between observations
  BUVAR1 = IQR
  checkLimits = .true.
  open (999, file=trim(dataPathOut) // 'errors.txt', status='unknown', access='APPEND' )
  if ( BUVAR1 <=  BLVAR1)  then
     checkLimits = .false.
     write (999, '(a25,i10)'  ) 'cell',  iCell
     write (999, '(a25,f10.5)') 'xMin',  minval(X)
     write (999, '(a25,f10.5)') 'xMax',  maxval(X)
     write (999, '(a25,f10.5)') 'p25',  Xemp(1)
     write (999, '(a25,f10.5)') 'p75',  Xemp(2)
     write (999, '(a40,f10.5)') 'lower boundary (xMax-xMin)/(nObs-1) =',  BLVAR1
     write (999, '(a40,f10.5)' ) 'upper boundary IQR = p75 - p25 =    ', BUVAR1
     write (999, '(f10.5)'    ) X
  end if
  close (999)
end function checkLimits

end module setGRG2


!******************************************************************************
! Leave-one-out log-likelihood Function to be maximized
! L. Samaniego, created 13.02.2010
!      modified         07.10.2011      intent()
!******************************************************************************
subroutine GCOMP(G,XU)
  use kernelSmoother, only     :  nObs, fKernel, X, nInter, pdf, evalPDF, hInt, nPar
  real(8), dimension(nPar)    :: XU
  real(8), dimension(1)       :: G(1)

!  real(8), dimension(:), intent(in)    ::  XU
!  real(8), dimension(:), intent(inout) ::  G
 
  real(8)                              ::  gCal, gCalk
  integer(4)                           ::  k, i
  !S
  G(1) = 0.0d0
  !
  ! Unbiased cross-validation criterium (to minimize)
  ! 1st term (see Scoot pg. 9)
  gCal = 0.0d0
  ! estimate density function
  call evalPDF(XU(1))
  do i = 1, nInter
    gCal = gCal + pdf(i,2)**2
  end do
  gCal = gCal * hInt
  ! 2nd term (umbiased estimator)
  gCalk = 0.0d0
  do k=1,nObs
     gCalk = gCalk + fKernel(XU(1), X(k), k)
  end do
  G(1)= gCal - 2.0d0 / dfloat(nObs) * gCalk
  !
  ! leave-one-out-likelihood estimator (Scoot, pg. 5)
  ! find kernel density at X(k) excluding observation k
  !do k=1,nObs
     !gCal = fKernel(XU(1), X(k), k)
     !if ( gCal < eps ) cycle
     ! minimize the function instead of maximizing
     !G(1)= G(1) - dlog (gCal)
  !end do
end subroutine GCOMP
