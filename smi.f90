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

  subroutine calSMI
    use mo_kind,             only  : i4
   ! use numerical_libraries, only : DGCDF
    use GCDF_INT
    use InputOutput,         only : nMy, nYears, nMonths, nCells, Z, WriteResultsKernel, grid, &
                                    SMI_flag, eMask
    use kernelSmoother,      only : nObs, nInter, evalPDF, evalEDF, hOpt, X, edf, pdf, hOptDB
    use setGRG2,             only : OPTI, checkLimits
    ! local variables
    integer(i4)                                  :: i, j, k, m, iStatus
    !
    ! Initialization
    nObs = nYears                                                      ! monthly PDF
    if ( .not. allocated(SMIp) ) allocate ( SMIp(nCells,nMonths), eMask(nCells) )
    allocate ( X(nObs), edf(0:nObs+1,2), pdf(nInter,3) )

    ! check limits
    iStatus = 0
    eMask   = 0
    do i= 1, nCells
       do j = 1, nMy
          ! select values for month j (1...12)
          X(:) = Z(i, j : nMonths : nMy)
          ! find optimal bandwidth for PDF of cell i in month j
          if ( checkLimits(i) == .false. ) then
             iStatus = iStatus + 1
             eMask(i) = 1
             Z(i, j : nMonths : nMy ) = grid%nodata_value  !here a posible check
          end if
       end do
    end do
    !
    if (iStatus > 0) then
       ! save error locations
       call WriteNetCDF(6)
       print*, 'WARNING: Input mistake... '
       print*, 'At least one cell has boundary problems      '
       print*, 'SMI will not be estimated in those cells  ...'

       print*, 'Program must stop.'
       stop
    end if
    !
    SMIp = grid%nodata_value
    do i= 1, nCells
       do j = 1, nMy
          ! select values for month j (1...12)
          X(:) = Z( i, j : nMonths : nMy)
          if ( any( X == grid%nodata_value) ) cycle
          ! find optimal bandwidth for PDF of cell i in month j
          if (SMI_flag == 3) then
           if ( hOptDB(i,j) < 0.0_dp  ) cycle
           call evalPDF(hOptDB(i,j))
          else  
            call OPTI
            ! store to print after
            hOptDB(i,j)= hOpt
            call evalPDF(hOpt)
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
             SMIp(i,m) = GCDF( Z(i,m), pdf(:,1), pdf(:,2), IOPT=4 )
          end do
       end do
       !print*, 'SMIp of cell: ', i
    end do
    ! saving hOpt database
    print*, 'SMIp estimated ...'
    if (SMI_flag == 1 .or. SMI_flag == 4 .or. SMI_flag == 5) call WriteResultsKernel(4,0,0)   
    print*, 'hOpt saved...'
  end subroutine calSMI
end module SMIndex
