!**********************************************************************
!  SOIL MOISTURE INDEX                                                *
!  PURPOSE     Estimate SMI based on monthly values of SM_mHM         *
!  CREATED     Luis Samaniego, 15.02.2011                             *
!**********************************************************************
module SMIndex
  use mo_kind,          only: i4, dp
  use InputOutput,      only: SMIp
  use InputOutput,      only: WriteNetCDF
  implicit none
!
contains

  subroutine calSMI( opt_h, SM_est, mask, nodata )
    use mo_kind,             only  : i4, sp
   ! use numerical_libraries, only : DGCDF
   ! use GCDF_INT
    use mo_interpol,         only : interpol
    use InputOutput,         only : nMy, nYears, nMonths, nCells, WriteResultsKernel, &
                                    SMI_flag, eMask, SMI
    use kernelSmoother,      only : nObs, nInter, evalPDF, evalEDF, hOpt, X, edf, pdf, hOptDB
    ! use setGRG2,             only : OPTI, checkLimits
    use mo_kernel,           only : kernel_density_h
    ! input variables
    logical,                  intent(in)  :: opt_h ! optimize kernel width
    logical,  dimension(:,:), intent(in)  :: mask
    real(sp),                 intent(in)  :: nodata
    real(sp), dimension(:,:), intent (in) :: SM_est

    ! local variables
    integer(i4)                                  :: i, j, k, m, iStatus
    integer(i4), dimension(1)                    :: minpos
    integer(i4)                                  :: nrows
    integer(i4)                                  :: ncols
    !
    ! Initialization
    nrows = size( mask, 1 )
    ncols = size( mask, 2 )
    nObs = nYears                                      ! monthly PDF
    if ( .not. allocated(SMIp) ) allocate ( SMIp(nCells,nMonths), eMask(nCells) )
    allocate ( X(nObs), edf(0:nObs+1,2), pdf(nInter,3) )
    allocate ( hOptDB( ncells, nMy ) )

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
       do j = 1, nMy
          ! select values for month j (1...12)
          X(:) = SM_est( i, j : nMonths : nMy)
          if ( any( X .eq. nodata) ) cycle
          ! find optimal bandwidth for PDF of cell i in month j
          if (SMI_flag == 3) then
           if ( hOptDB(i,j) < 0.0_dp  ) cycle
           call evalPDF(hOptDB(i,j))
          else  
             ! call OPTI
             ! ! store to print after
             ! hOptDB(i,j)= hOpt
             if ( opt_h ) then
                hOptDB(i,j) = kernel_density_h( X(:), silverman = .false. )
             else
                hOptDB(i,j) = kernel_density_h( X(:), silverman = .true. )
             end if
             call evalPDF(hOptDB(i,j))
          end if
          ! call evalEDF
          ! call WriteResultsKernel(1,i,j)
          ! call WriteResultsKernel(2,i,j)
          ! call WriteResultsKernel(3,i,j)
          ! print*, 'PDF of cell: ', i, ' month: ', j
          !
          ! find quantiles for each month using the estimated PDF
          do k = 1, nYears
             m = (k-1) * nMy + j
             !SMIp(i,m) = DGCDF( Z(i,m), 4, nInter, pdf(:,1), pdf(:,2))
             ! SMIp(i,m) = GCDF( Z(i,m), pdf(:,1), pdf(:,2), IOPT=4 )
             ! calculate argmin
             minpos = minloc( abs(pdf(:,1) - SM_est(i,m)) )
             SMIp(i,m) = sum( pdf(:minpos(1),2) )
          end do
          !SMIp(i,:) = interpol( pdf(:,2), pdf(:,1), Z(i,:) )
       end do
       !print*, 'SMIp of cell: ', i
    end do
    !
    ! put SMIp to SMIp
    allocate( SMI( nrows, ncols, nMonths ) )
    do j = 1, nMonths 
       SMI (:,:,m) = unpack ( SMIp(:,m),  mask, real(nodata,dp) )
       !
       ! filter for possible values
       where (mask .and. SMI(:,:,m) .lt. 0.0_dp )
          SMI(:,:,m) = 0.5_dp
       end where
    end do
     print *, 'calculated SMI...ok'
  end subroutine calSMI
end module SMIndex
