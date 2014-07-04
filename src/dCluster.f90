!*********************************************************************
!  PURPOSE: Find drought clusters and statistics
!           1) truncate SMIp < SMI_th
!
!  Luis Samaniego, created  28.02.2011
!                  updated  05.10.2011
!*********************************************************************
subroutine droughtIndicator( mask, nodata )
  use mo_kind,     only                          : i4, sp, dp
  use InputOutput, only                          : nMonths, nCells, &
                                                   SMIp, SMI_thld, &
                                                   SMIc, SMI, cellCoor
  ! input variable
  logical, dimension( :, : ), intent(in)        :: mask
  real(sp),                   intent(in)        :: nodata
  ! local variables
  integer(i4)                                   :: i, j, m, k
  integer(i4)                                   :: nrows
  integer(i4)                                   :: ncols
  ! initialize
  nrows = size( mask, 1 )
  ncols = size( mask, 2 )
  !
  print *, 'Drought  th: ', SMI_thld
  if ( .not. allocated (SMIc) ) allocate ( SMIc( nrows, ncols, nMonths) )
  do m = 1, nMonths
     !
     ! filter for possible error within domain
     where (SMI(:,:,m) .le. SMI_thld .and. SMI(:,:,m) .ne. nodata )
        SMIc(:,:,m) = 1
     elsewhere (SMI(:,:,m) .ne. nodata )
        SMIc(:,:,m) = 0
     elsewhere
        SMIc(:,:,m) = nodata
     end where
  end do

  !
  allocate (cellCoor(nCells,2))
  k = 0
  do j=1,ncols
     do i=1,nrows
        if ( SMIc(i,j,1) .ne. int( nodata, i4 ) ) then
           k = k + 1
           cellCoor(k,1) = i
           cellCoor(k,2) = j
           ! if (k==7645.or.k==2479) print*, i,j,k
        end if
     end do
  end do
  !
end subroutine droughtIndicator

!**********************************************************************
!  PROGRAM: dCluster
!
!  PURPOSE:  cluster drought clusters in space and time
!  DATE:     developed in Budapest, 10-11.03.2011
!**********************************************************************
subroutine ClusterEvolution( nrows, ncols, nodata )
  use mo_kind, only                                    : i4
  ! use numerical_libraries, only                        : SVIGN
  use mo_sort, only:                                     sort
  use InputOutput, only                                : SMIc,&
                                                         idCluster, &
                                                         nMonths, &
                                                         shortCnoList, &
                                                         nClusters, &
                                                         nCellInter,thCellClus
  implicit none
  ! input variables
  integer(i4),                            intent(in)  :: nrows
  integer(i4),                            intent(in)  :: ncols
  integer(i4),                            intent(in)  :: nodata
  ! local variables
  integer(i4)                                         :: i, j, t, k, l
  integer(i4), dimension(:), allocatable              :: nC
  integer(i4), parameter                              :: factor = 1e3
  integer(i4)                                         :: ncInter
  integer(i4), dimension(:,:), allocatable            :: cno
  integer(i4), dimension(:), allocatable              :: cnoList
  integer(i4), dimension(:), allocatable              :: vec
  integer(i4)                                         :: maxNc, nTotal, idRep
  character (256) ::dummy
  integer(i4), dimension(:,:), allocatable            :: idCnew
  !
  !
  if (.not. allocated(idCluster) )  allocate   (idCluster(nrows, ncols, nMonths))
  if (.not. allocated(nC) )         allocate   (nC(nMonths))
  if ( allocated(cno) )             deallocate (cno, vec, cnoList)
  if (.not. allocated(idCnew))  allocate   (idCnew(nrows, ncols))

  ! ---------------------------------------------
  ! CLUSTERING IN SPACE and remove small clusters
  ! ---------------------------------------------
  do t=1,nMonths
     if ( count(SMIc(:,:,t) == 1) < thCellClus) then
        nC(t) = 0
        idCluster(:,:,t) = nodata
        cycle
     end if
     call findClusters (t, idCluster(:,:,t), nC(t), nrows, ncols, nodata)
     !print*, 'Finding clusters time :', t, nC(t)
  end do
  print*, 'Clusters in space done ...'
  !
  ! maximum number of clusters at all time steps
  maxNc = maxval(nC(:))
  nTotal = nMonths*maxNc
  allocate ( cno(nMonths,maxNc), vec(nTotal), cnoList(nTotal) )
  cno = -9
  !
  ! set unique cluster numbers and store them
  ! NOTE 
  !      e.g. month*factor + running nr.,
  !           5*1000+1 = 5001              => fisrt cluster in 5th month
  ! 
  do t=1,nMonths
     if (nC(t) == 0) cycle
     do i=1, nC(t)
        cno(t,i) = t*factor + i
        where ( idCluster(:,:,t) == i ) idCluster(:,:,t) = cno(t,i)
     end do
  end do
  ! ---------------------------------
  ! FIND CLUSTERS IN TIME
  ! determine intersection sets, join
  ! ---------------------------------
  do t=2,nMonths
     do i=1,nC(t)
        if (cno(t,i) == -9) cycle
        do j=1,nC(t-1)
           if (cno(t-1,j) == -9) cycle
           ncInter = count ( idCluster(:,:,t) == cno(t,i) .and. idCluster(:,:,t-1) == cno(t-1,j) )
           if ( ncInter >= nCellInter ) then
              ! renumber all from 1 to t
              where ( idCluster(:,:,1:t) == cno(t-1,j) ) idCluster(:,:,1:t) =  cno(t,i)
              ! rename cluster id from cno
              where ( cno(1:t,:) == cno(t-1,j) ) cno(1:t,:) = cno(t,i)
           end if
        end do
     end do
  end do
  print*, 'Clustering in time done ...'
  ! --------------------------
  ! COMPILE CLUSTER IDS + list
  ! --------------------------
  vec = -9
  cnoList = pack (cno, mask = cno > 0, vector=vec)
  !
  ! 1. sort
  idRep = -9
  forall(i=1:nTotal) vec(i) = i
  ! call SVIGN(nTotal,cnoList,cnoList)
  call sort(cnoList)
  ! 2.  remove repeated values
  do i=nTotal, 2,-1
     if (cnoList(i) < 0 ) cycle
     if (.not. cnoList(i-1) == cnoList(i) ) then
        if (idRep > 0) idRep = -9
        cycle
     end if
     if (idRep == -9) idRep = cnoList(i)
     cnoList(i) = -9
  end do
  ! call SVIGN(nTotal,cnoList,cnoList)
  call sort(cnoList)
  !
  ! 3. final consolidated list
  nClusters = count(cnoList > 0)
  allocate ( shortCnoList(nClusters) )
  shortCnoList = pack(cnoList, mask = cnoList > 0)
  !
  print*, '# Consolidated Clusters :', nClusters
  !
  deallocate(vec, nC, cnoList, cno )
end subroutine ClusterEvolution

!-------------------------------------------------------
! SVAT statistics
!-------------------------------------------------------
subroutine ClusterStats( nrows, ncols )
  use mo_kind, only                                   : i4, dp
  !  use numerical_libraries, only                       : SVIGN , DSVRGN, DEQTIL
  use InputOutput, only                                : aDA, aDD, TDM, DTMagEvol, DAreaEvol, &
                                                         nCells, SMI, &
                                                         nClusters, nMonths, idCluster, &
                                                         shortCnoList, deltaArea, &
                                                         dASevol, nBasins, Basin_Id, SMI_thld
  implicit none
  ! input variables
  integer(i4), intent(in) :: nrows
  integer(i4), intent(in) :: ncols
  ! local variables
  integer(i4)                                         :: ic
  integer(i4)                                         :: t, i
  integer(i4) , dimension(:), allocatable             :: counterA
  integer(i4) , dimension(:,:), allocatable           :: aDDG
  real(dp) , dimension(:,:), allocatable              :: mSev
  !
  allocate ( aDDG     (nrows, ncols) )
  allocate ( DAreaEvol(nMonths,nClusters) )
  allocate ( DTMagEvol(nMonths,nClusters) )
  allocate ( aDD      (nClusters)   )
  allocate ( aDA      (nClusters)   )
  allocate ( counterA (nClusters)   )
  allocate ( TDM      (nClusters)   )
  allocate ( dASevol  (nMonths,2,nBasins+1)   )       ! (area, severity) whole Germany => nBasins+1
  allocate ( mSev     (nrows, ncols) )      ! mean  severity
  !
  DAreaEvol = 0.0_dp
  DTMagEvol = 0.0_dp
  counterA  = 0
  do ic = 1, nClusters
     !print*, 'Statistics of cluster : 'shortCnoList(ic)
     aDDG = 0
     do t = 1, nMonths
        where (idCluster(:,:,t) == shortCnoList(ic) )
           ! duration
           aDDG = aDDG + 1
        end where
        !
        !  drought area evolution
        DAreaEvol(t,ic) =  real( count(idCluster(:,:,t) == shortCnoList(ic) ), dp) / real(nCells, dp)
        if (DAreaEvol(t,ic) > 0) counterA(ic) = counterA(ic) + 1
        ! total magnitude (NEW  SM_tr -SMI) !!!
        DTMagEvol(t,ic) = sum( (SMI_thld - SMI(:,:,t)), mask = idCluster(:,:,t) == shortCnoList(ic) )
     end do
     ! AVERAGE (MEAN) DURATION
     aDD(ic) = real( sum(aDDG, mask = aDDG > 0),dp ) / real( count(aDDG > 0), dp)
     ! TOTAL MAGNITUD (sum over space and time of SMI over all
     !                 cells afected by the event ic )
     TDM(ic) = sum( DTMagEvol(:,ic) )
  end do
  !  AVERAGE (MEAN) DROUGHT AREA per cluster event
  aDA = sum(DAreaEvol, DIM=1) / real(counterA, dp)
  !
  ! (NEW) Time evolution of the area and monthly severity
  ! whole domain
  do t = 1, nMonths
     dASevol(t,1,nBasins+1) = real( count( idCluster(:,:,t) > 0 ), dp ) / real( nCells, dp) * 1e2_dp
     where (idCluster(:,:,t) > 0 )
         mSev(:,:) = 1.0_dp - SMI(:,:,t)
     end where
     dASevol(t,2,nBasins+1) = sum( mSev(:,:), mask = idCluster(:,:,t) > 0 ) / real( nCells, dp)  
  end do
  ! evolution at basin wise
  do t = 1, nMonths
    do i = 1, nBasins
      dASevol(t,1,i) = real ( count( idCluster(:,:,t) > 0 .and.  Basin_Id(:,:) == i  ), dp) &
                       / real(count(Basin_Id(:,:) == i), dp)  * 1e2_dp
      where (idCluster(:,:,t) > 0 .and.  Basin_Id(:,:) == i )
         mSev(:,:) = 1.0_dp - SMI(:,:,t)
      end where
      dASevol(t,2,i) = sum( mSev(:,:), mask = idCluster(:,:,t) > 0 .and. Basin_Id(:,:) == i ) &
                      / real(count(Basin_Id(:,:) == i), dp)
    end do
  end do
 print*, 'Cluster statistics were estimated ... '
  deallocate ( counterA, aDDG, mSev )
end subroutine ClusterStats

!-------------------------------------------------------
! SAD analysis
!-------------------------------------------------------
subroutine calSAD(iDur, nrows, ncols, nodata)
  use mo_kind, only                                    : i4, dp
  ! use numerical_libraries, only                        : DSVRGN, DEQTIL, SVIGP
  use mo_sort, only                                    : sort, sort_index
  use mo_percentile, only                              : percentile
  use InputOutput, only                                : shortCnoList, deltaArea, SMI, &
                                                         nInterArea, nEvents, nClusters,  idCluster, &
                                                         SAD, SADperc, DAreaEvol, severity, nDsteps, &
                                                         nDurations, durList, eventId, eIdPerm, &
                                                         writeResultsCluster, nLargerEvents, &
                                                         nCells, cellArea, nMonths, nQProp, QProp
  implicit none
  ! input variable
  integer(i4), intent(in)                             :: iDur
  integer(i4), intent(in)                             :: nrows
  integer(i4), intent(in)                             :: ncols
  integer(i4), intent(in)                             :: nodata
  ! local variable
  integer(i4)                                         :: t, i, k
  integer(i4)                                         :: ic
  integer(i4)                                         :: iDc
  integer(i4)                                         :: d
  integer(i4)                                         :: ms, me
  integer(i4)                                         :: ke
  integer(i4)                                         :: eCounter, eC
  integer(i4)                                         :: ncic
  integer(i4)                                         :: nIntA
  real(dp), dimension(:), allocatable                 :: sevP
  !
  real(dp), dimension(nQProp)                         :: Xlo
  real(dp), dimension(nqProp)                         :: Xhi
  integer(i4)                                         :: nObs
  integer(i4)                                         :: nMiss
  integer(i4), dimension(:), allocatable              :: vec
  !
  if (iDur == 1) then
     ! set up SAD curves
     nInterArea = int( ceiling( real(nCells,dp) / real(deltaArea, dp) ), i4 )
     !
     ! total number of events to be evaluated
     nEvents = count (DAreaEvol > 0.0_dp )
     nLargerEvents = nEvents
     allocate (  eventId   (nEvents, 3) )                !  dim1: running Nr.
                                                         !  dim2: 1 == cluster Id,
                                                         !        2 == month of ocurrence,
                                                         !        3 == number of cells
     allocate (  eIdPerm   (nEvents)    )
     allocate (  SAD       (nInterArea, 2, nLargerEvents ) )    !  dim2: 1 == area, 2 == Severity
     allocate (  SADperc   (nInterArea, nQProp ) )
     !
     ! keep event identification
     eCounter = 0
     do ic = 1,  nClusters
        do t = 1, nMonths
           if (DAreaEvol(t,ic) > 0.0_dp ) then
              eCounter = eCounter + 1
              eventId(eCounter,1) = shortCnoList(ic)
              eventId(eCounter,2) = t
              eventId(eCounter,3) = count( idCluster(:,:,t) == shortCnoList(ic) )
           end if
        end do
     end do
     ! sort event in order of magnitude
     forall(i=1:nEvents) eIdPerm(i) = i
     allocate (vec(nEvents))
     ! call SVIGP (nEvents, eventId(:,3), vec, eIdPerm)
     eIdPerm = sort_index( eventID(:,3) )
     deallocate (vec)
  end if
  !
  ! SAD
  print*, 'SAD for duration started: ', durList (iDur)
  SAD     = -9.0_dp
  SADperc = -9.0_dp
  nDsteps = ceiling (real(nMonths,dp) / real(durList (iDur),dp) )
  if ( allocated (severity) )  deallocate(severity)
  allocate ( severity (nrows, ncols, nDsteps ) )
  severity = nodata
  !
  ! estimate severities for a given duration over all time steps
  do d = 1, nDsteps
     ms = (d-1)*durList(iDur) + 1
     me = ms + durList(iDur) - 1
     if (me > nMonths) me = nMonths
     where ( SMI(:,:,1) .ne. nodata )
        severity(:,:,d) = 1.0_dp - sum( SMI(:,:,ms:me), DIM=3 ) / real( me-ms+1, dp )
     end where
  end do
  !
  ! estimate curves for larger events ONLY
  do eC = 1, nLargerEvents
     eCounter = eIdPerm( nEvents + 1 - eC )
     iDc = eventId(eCounter,1)
     t   = eventId(eCounter,2)
     d   = (t-1)/durList (iDur) + 1
     ! number of cells of the event in eCounter
     ncic = eventId(eCounter,3)
     if (allocated (sevP)) deallocate(sevP)
     allocate (sevP(ncic))
     sevP = pack (severity(:,:,d), mask = idCluster(:,:,t) == iDc )
     !
     ! ----------------------------------------------
     ! fast approach
     ! rank severities and estimate cummulative areas
     ! ----------------------------------------------
     ! call DSVRGN (ncic,sevP,sevP)
     call sort( sevP )
     nIntA = int( ceiling( real(ncic,dp) / real(deltaArea, dp) ), i4 )
     do k = 1, nIntA
        ke = k*deltaArea
        if (ke > ncic ) ke = ncic
        SAD(k, 1, eC ) = real(ke, dp) * cellArea                    ! cumulative area
        SAD(k, 2, eC ) = sum( sevP(ncic+1-ke:ncic) ) / real(ke, dp) ! sorted by incresing value
     end do
  end do
  !
  ! estimate percentilles for all intervals
  do k = 1, nInterArea
     nObs = count(SAD(k,2,:) > 0.0_dp)
     if (.not. nObs > 1) cycle
     if (allocated (sevP)) deallocate(sevP)
     allocate (sevP(nObs))
     sevP = pack (SAD(k,2,:), mask = SAD(k,2,:) > 0.0_dp )
     ! call DEQTIL (nObs, sevP, nQprop, Qprop, SADperc(k,:), Xlo, Xhi, nMiss)
     SADperc(k,:) =  percentile( sevP, Qprop )
  end do
  !
  if ( allocated(sevP) ) deallocate ( sevP )
end subroutine calSAD

!-------------------------------------------------------
! clustering in space
! assign a unique running number at time t
!-------------------------------------------------------
subroutine findClusters (t,iC,nCluster, nrows, ncols, nodata)
  use mo_kind, only                                    :  i4
  use InputOutput, only                                :  cellCoor, thCellClus, &
       SMIc, nCells

  implicit none
  integer(i4), intent(in)                               :: nrows
  integer(i4), intent(in)                               :: ncols
  integer(i4), intent(in)                               :: nodata
  integer(i4), intent(in)                               :: t                   !  time step
  integer(i4), intent(out)                              :: nCluster            !  number of clusters larger than a threshold
  integer(i4), dimension(nrows,ncols), intent(out)      :: iC
  integer(i4)                                           :: i, j, k, klu, klu1, kk
  integer(i4)                                           :: iul, idr, jul, jdr
  integer(i4)                                           :: krow, kcol
  integer(i4), dimension(:), allocatable                :: cno, vec
  integer(i4), dimension(:), allocatable                :: nCxCluster
  integer(i4)                                           :: nClusterR, ncc
  !
  iC = nodata
  nCluster = 0
  do k=1,nCells
     krow = cellCoor(k,1)
     kcol = cellCoor(k,2)
     iul  = max(krow-1,1)
     idr  = min(krow+1,nrows)
     jul  = kcol
     jdr  = min(kcol+1,ncols)
     ! SMIc of k
     klu  = SMIc(krow, kcol, t)
     if (klu /= 1) cycle
     if ( (klu == 1) .and. (iC(krow, kcol) == nodata) ) then
        nCluster=nCluster+1
        iC(krow, kcol) = nCluster
     end if
     do j=jul, jdr
        do i= iul, idr
           if (iC(i,j) == nodata .and. &
                SMIc(i,j,t) == klu                     ) then
              iC(i,j) = iC(krow, kcol)
           end if
        end do
     end do
  end do
  ! consolidate clusters
  if ( nCluster > 0 ) then
     allocate ( cno(nCluster), vec(nCluster) )
     cno = (/(i, i=1,nCluster)/)
     vec = -9
     nClusterR = nCluster
     do k=nCells,1,-1
        krow = cellCoor(k,1)
        kcol = cellCoor(k,2)
        klu  = iC(krow, kcol)
        iul  = max(krow-1,1)
        idr  = min(krow+1,nrows)
        jul  = kcol
        jdr  = min(kcol+1,ncols)
        do j=jul, jdr
           do i= iul, idr
              klu1 = iC(i, j)
              if ( klu  /= nodata .and. &
                   klu1 /= nodata .and. &
                   klu  /= klu1                       ) then
                 cno(klu1) = -9
                 nClusterR = nClusterR - 1
                 where ( iC == klu1 ) iC = klu
              end if
           end do
        end do
     end do
     !
     ! delete small clusters < thesh. area
     do i=1, nCluster
        if (cno(i) == -9 ) cycle
        ncc = count(iC == cno(i))
        if ( ncc <= thCellClus ) then
           where (iC == cno(i)) iC = nodata
           cno(i) = -9
           nClusterR = nClusterR -1
        end if
     end do
     !
     ! reordering
     cno = pack (cno, mask = cno > 0, vector=vec)
     where (cno <= 0) cno = 0
     !
     !
     allocate (nCxCluster(nClusterR))
     do i=1, nClusterR
        where (iC == cno(i)) iC = i
        nCxCluster(i) = count(iC == i)
     end do
     !
     deallocate (cno, nCxCluster)
  end if
  nCluster = nClusterR
end subroutine findClusters

