!**********************************************************************
!  SOIL MOISTURE INDEX                                                *
!  PURPOSE     Estimate SMI based on monthly values of SM_mHM         *
!  CREATED     Luis Samaniego, 15.02.2011                             *
!  HISTORY     Stephan Thober, 16.12.2014 - added evalSMI capability  *
!              Stephan Thober, 07.11.2017 - added inversion of CDFs   *
!              Stephan Thober, 02.08.2019 - added non-full year       *
!**********************************************************************
module mo_smi
  !
  use mo_global_variables, only: period
  use mo_kind,          only: dp
  !
  implicit none
  ! public routines
  public :: optimize_width
  public :: calSMI
  public :: invSMI
  !

  private

  
contains

  ! subroutine for estimating SMI for first array
  subroutine optimize_width( opt_h, silverman_h, SM, nCalendarStepsYear) ! tmask,
    use omp_lib
    use mo_kind,          only : i4, sp
    use mo_kernel,        only : kernel_density_h
    implicit none
    
    ! input variables
    logical,                  intent(in) :: silverman_h ! optimize kernel width
 !   logical,  dimension(:,:), intent(in) :: tmask
    real(sp), dimension(:,:), intent(in) :: SM
    integer(i4),              intent(in) :: nCalendarStepsYear

    ! output variables
    real(dp), dimension(:,:), intent(inout) :: opt_h

    ! local variables
    integer(i4)                             :: ii     ! cell index
    integer(i4)                             :: mm     ! month index
    real(dp), dimension(:),     allocatable :: X

    allocate( X(size(SM,2) / nCalendarStepsYear) )

    !$OMP parallel default(shared) &
    !$OMP private(ii,mm,X)
    !$OMP do COLLAPSE(2) SCHEDULE(STATIC)
    do ii = 1, size( SM, 1 )
       do mm = 1, nCalendarStepsYear    
         if ((mod(ii, 10) .eq. 0_i4) .and. (mm .eq. nCalendarStepsYear)) print *, ii, mm
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
    !$OMP end do
    !$OMP end parallel
    deallocate(X)
    !! END
  end subroutine optimize_width


 
  !!======================================================
  !! ESTIMATE SMI 
  !!======================================================
  !! subroutine for calculating SMI for second array with pdf of first one
  subroutine calSMI( hh, sm_est, sm_eval,  nCalendarStepsYear, SMI, per, per_eval)  ! tmask_est,
    
    use mo_kind,          only: i4, sp, dp
    use mo_kernel,        only: kernel_cumdensity
    use mo_smi_constants, only: nodata_dp
    
    implicit none

    ! input variables
    real(dp), dimension(:,:), intent(in)  :: hh
    real(sp), dimension(:,:), intent(in)  :: sm_est
    real(sp), dimension(:,:), intent(in)  :: sm_eval
    integer(i4),              intent(in)  :: nCalendarStepsYear
    type(period),             intent(in)  :: per
    type(period),             intent(in)  :: per_eval

    ! output variables
    real(sp), dimension(:,:), intent(out) :: SMI

    ! local variables
    integer(i4)                           :: mm  ! loop index
    integer(i4)                           :: ii  ! cell index
    real(dp), dimension(:), allocatable   :: cdf
    real(dp), dimension(:), allocatable   :: X_est
    real(dp), dimension(:), allocatable   :: X_eval
    integer(i4)                           :: n_time
    logical,                allocatable   :: t_mask(:)
    integer(i4),            allocatable   :: time(:)

    if (nCalendarStepsYear .eq. 12) then
      call get_time_indizes(time, per%j_start, per_eval%m_start, per_eval%n_months, nCalendarStepsYear)
    else
      call get_time_indizes(time, per%j_start, per_eval%d_start, per_eval%n_days, nCalendarStepsYear)
    end if
    n_time = size(time, 1)
    
    allocate ( X_est (size(SM_est,2) / nCalendarStepsYear) )
    X_est(:) = nodata_dp
    
    do ii = 1, size(SM_est,1)          ! cell loop
      do mm = 1,  nCalendarStepsYear   ! calendar time 

        t_mask = (time .eq. mm)
        n_time = count(t_mask)
        
        allocate ( X_eval(n_time) )
        allocate ( cdf   (n_time) )
        
        X_est(:)  = real(SM_est ( ii, mm:size(SM_est,2):nCalendarStepsYear ),  dp)
        ! X_eval(:) = real(SM_eval( ii, mm:size(SM_eval,2):nCalendarStepsYear ),  dp)
        
        !X_est(:)  = pack( real(SM_est(ii,:),  dp), tmask_est(:,mm) )
        X_eval(:) = pack(real(SM_eval(ii,:), dp), t_mask)
        
        cdf(:) = kernel_cumdensity(x_est, hh(ii,mm), xout=x_eval)
        SMI(ii, mm:size(SM_eval,2):nCalendarStepsYear) = real(cdf(:), sp)
        
        deallocate( X_eval, cdf )
      end do
    end do

    deallocate( X_est )

  end subroutine calSMI

  ! create objective function of kernel_cumdensity and minimize it using
  ! nelmin because function is monotone
  subroutine invSMI(sm_est, hh, SMI_invert, nCalendarStepsYear, &
       SM_invert)
    use mo_kind,          only: i4, sp, dp
    use mo_smi_constants, only: nodata_sp, nodata_dp
    use mo_kernel,        only: kernel_cumdensity
    
    implicit none
    
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
    integer(i4)                         :: xx_n_sample
    integer(i4), dimension(1)           :: idx_invert
    real(dp)                            :: hh_est
    real(dp)                            :: xx_min, xx_max, xx_h
    real(dp), dimension(:), allocatable :: y_inv
    real(dp), dimension(:), allocatable :: xx_cdf, yy_cdf
    real(dp), dimension(:), allocatable :: xx_est
   
    ! initialize extents
    n_cells        = size(SMI_invert, 1)
    n_years_est    = size(sm_est, 2) / nCalendarStepsYear
    n_years_invert = size(SMI_invert, 2) / nCalendarStepsYear
    xx_n_sample    = 2000_i4 ! gives precision of at least 0.0005 in SM
    allocate(xx_cdf(xx_n_sample))
    allocate(yy_cdf(xx_n_sample))
    xx_cdf = nodata_dp
    yy_cdf = nodata_dp

    ! initialize output array
    allocate(xx_est(n_years_est))
    allocate(y_inv(n_years_invert))
    allocate(SM_invert(n_cells, size(SMI_invert, 2)))
    xx_est    = nodata_dp
    y_inv     = nodata_dp
    SM_invert = nodata_sp

    print *, ''
    print *, '  start inversion of CDF'
    !$OMP parallel default(shared) &
    !$OMP private(mm, yy, xx_est, hh_est, y_inv, xx_min, xx_max, xx_h, xx_cdf, yy_cdf, idx_invert)
    !$OMP do
    do ii = 1, n_cells
       if (modulo(ii, 1000) .eq. 0) print *, ii, n_cells
       do mm = 1, nCalendarStepsYear
          xx_est(:) = real(SM_est    ( ii, mm:size(sm_est, 2):nCalendarStepsYear    ),  dp)
          hh_est    = hh(ii, mm)
          y_inv(:)  = real(SMI_invert( ii, mm:size(SMI_invert, 2):nCalendarStepsYear),  dp)

          ! sample cdf
          xx_min    = max(0._dp, minval(xx_est - 10._dp * hh_est))
          xx_max    = min(1._dp, maxval(xx_est + 10._dp * hh_est))
          xx_h      = (xx_max - xx_min) / real(xx_n_sample, dp)
          do yy = 1, xx_n_sample
             xx_cdf(yy) = xx_min + (yy - 1_i4) * xx_h
          end do
          xx_cdf = merge(1._dp, xx_cdf, xx_cdf .gt. 1._dp)
          yy_cdf = kernel_cumdensity(xx_est, hh_est, xout=xx_cdf)
          
          do yy = 1, n_years_invert
             idx_invert = minloc(abs(y_inv(yy) - yy_cdf))
             SM_invert(ii, (yy - 1) * nCalendarStepsYear + mm) = xx_cdf(idx_invert(1))
          end do
       end do
    end do
    !$OMP end do
    !$OMP end parallel
    print *, '  finish inversion of CDF... ok'
    print *, ''

    ! free memory
    deallocate(y_inv, xx_est, xx_cdf, yy_cdf)
    
  end subroutine invSMI

  subroutine get_time_indizes(time, j_start, start, n_time, nCalendarStepsYear)

    use mo_kind, only : i4, dp
    use mo_julian, only : dec2date
    
    implicit none

    integer(i4), intent(in) :: j_start, start, n_time, nCalendarStepsYear
    integer(i4), allocatable, intent(out) :: time(:)
    integer(i4) :: ii, jj, cc, dd, mm, yy

    ! remove leap days
    cc = 0
    do ii = 1, n_time
      call dec2date(real(j_start + ii - 1, dp), dd=dd, mm=mm, yy=yy)
      if ((dd .eq. 29) .and. (mm .eq. 2)) cc = cc + 1_i4
    end do
    
    allocate(time(n_time - cc))
    time(:) = 0_i4

    jj = 1
    do ii = 1, size(time)
      time(ii) = jj
      jj = jj + 1_i4
      if (jj .gt. nCalendarStepsYear) jj = 1_i4
    end do
    
  end subroutine get_time_indizes
  
end module mo_smi
