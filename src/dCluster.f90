!*********************************************************************
!  PURPOSE: Find drought clusters and statistics
!           1) truncate SMIp < SMI_th
!
!  Luis Samaniego, created  28.02.2011
!                  updated  05.10.2011
!*********************************************************************
subroutine droughtIndicator
  use mo_kind,     only                          : i4, dp
  use InputOutput, only                          : nMonths, nCells, &
                                                   SMIp, SMI_thld, &
                                                   SMIc, SMI, grid, &
                                                   mask, cellCoor
  ! local variables
  integer(i4)                                   :: i, j, m, k
  !
  print *, 'Drought  th: ', SMI_thld
  if ( .not. allocated (SMI) ) allocate  ( SMI (grid%nrows, grid%ncols, nMonths) )
  if ( .not. allocated (SMIc) ) allocate ( SMIc(grid%nrows, grid%ncols, nMonths) )
  do m = 1, nMonths
     SMI (:,:,m) = unpack ( SMIp(:,m),  (mask /= int(grid%nodata_value,i4) ), grid%nodata_value  )
     !
     ! filter for possible error within domain
     where (mask /= int(grid%nodata_value,i4) .and. SMI(:,:,m) < 0.0_dp )
        SMI(:,:,m) = 0.5_dp
     end where
     where (SMI(:,:,m) <= SMI_thld .and. SMI(:,:,m) /= grid%nodata_value )
        SMIc(:,:,m) = 1
     elsewhere (SMI(:,:,m) /= grid%nodata_value )
        SMIc(:,:,m) = 0
     elsewhere
        SMIc(:,:,m) = grid%nodata_value
     end where
  end do

  ! estimate probability of occurence
  !    do i = 1, nCells
  !      dghtPr(i) = real(count(dgthCluster(i,:) == 1),dp) / real(nMonths,dp)
  !      write (99,*) count(dgthCluster(i,:) == 1)
  !    end do
  ! unpack dgthCluster
  !
  allocate (cellCoor(nCells,2))
  k = 0
  do j=1,grid%ncols
     do i=1,grid%nrows
        if ( SMIc(i,j,1) /=  int(grid%nodata_value,i4) ) then
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
subroutine ClusterEvolution
  use mo_kind, only                                    : i4
  ! use numerical_libraries, only                        : SVIGN
  use mo_sort, only:                                     sort
  use InputOutput, only                                : grid, SMIc,&
                                                         idCluster, &
                                                         nMonths, &
                                                         shortCnoList, &
                                                         nClusters, &
                                                         nCellInter,thCellClus, mask
  implicit none
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
  if (.not. allocated(idCluster) )  allocate   (idCluster(grid%nrows, grid%ncols, nMonths))
  if (.not. allocated(nC) )         allocate   (nC(nMonths))
  if ( allocated(cno) )             deallocate (cno, vec, cnoList)
  if (.not. allocated(idCnew))  allocate   (idCnew(grid%nrows, grid%ncols))

  ! ---------------------------------------------
  ! CLUSTERING IN SPACE and remove small clusters
  ! ---------------------------------------------
  do t=1,nMonths
     if ( count(SMIc(:,:,t) == 1) < thCellClus) then
        nC(t) = 0
        idCluster(:,:,t) = int(grid%nodata_value, i4)
        cycle
     end if
     call findClusters (t, idCluster(:,:,t), nC(t))
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
  ! >>>
  ! to be deleted
  ! add mask + print grids
!!$  idRep = -9
!!$  cnoList = -9
!!$  forall(i=1:nTotal) vec(i) = i
!!$  do t=200,212
!!$      k = 0
!!$      do j=1,grid%ncols
!!$      do i=1,grid%nrows
!!$      
!!$         if (idCluster(i, j,t) > 0) then
!!$           k=k+1
!!$           cnoList(k) = idCluster(i, j,t)
!!$         end if   
!!$    
!!$      end do
!!$      end do
!!$  end do
!!$   call SVIGN(nTotal,cnoList,cnoList)
!!$  nClusters = count(cnoList > 0)
!!$  allocate ( shortCnoList(nClusters) )
!!$  shortCnoList = pack(cnoList, mask = cnoList > 0)

  !    last version to print clusters

!!$    do t=200,212
!!$     idCnew = mask-1
!!$    
!!$    do j=1,grid%ncols
!!$    do i=1,grid%nrows
!!$      if ( mask(i,j) /=  int(grid%nodata_value,i4) ) then
!!$         if (idCluster(i, j,t) > 0) then
!!$         do l=1, nClusters
!!$           if (idCluster(i, j,t) == shortCnoList(l)) exit
!!$           end do
!!$         idCnew(i,j)=l      
!!$         end if
!!$      end if
!!$    end do
!!$    end do
!!$    where (mask == int(grid%nodata_value,i4) ) idCnew = int(grid%nodata_value,i4)
!!$    write (dummy,'(i3.3)') t
!!$    open(999,file=trim(dummy)//'.asc', status='unknown')
!!$    write (999, 1)  'ncols',        grid%ncols
!!$    write (999, 1)  'nrows',        grid%nrows
!!$    write (999, 2)  'xllcorner',    grid%xllcorner
!!$    write (999, 2)  'yllcorner',    grid%yllcorner
!!$    write (999, 2)  'cellsize',     grid%cellsize
!!$    write (999, 2)  'NODATA_value', grid%nodata_value
!!$
!!$    
!!$1   format (a12, 2x, i10)
!!$2   format (a12, 2x, f10.1)
!!$3   format (a12, 2x, a)

  ! write(999,106) 'DSAA'
  ! write(999,108) grid%ncols+2, grid%nrows+2
  ! write(999,104) grid%xllcorner - 0.5d0*dfloat(grid%cellsize), &
  !      (grid%xllcorner - 0.5d0*dfloat(grid%cellsize)) + dfloat((grid%ncols + 1)*grid%cellsize)
  ! write(999,104) grid%yllcorner - 0.5d0*dfloat(grid%cellsize), &
  !      (grid%yllcorner - 0.5d0*dfloat(grid%cellsize)) + dfloat((grid%nrows + 1)*grid%cellsize)
  ! write (999,108) int(grid%nodata_value,i4), nClusters
  ! do k = grid%nrows, 1,-1
!!$      do k= 1, grid%nrows

!!$        if(k == grid%nrows) then
!!$           write (999, 114) (int(grid%nodata_value,i4), l=1, grid%ncols+2)
!!$        end if
  !     write (999, 114) int(grid%nodata_value,i4), ( idCnew(k, l), l=1, grid%ncols), int(grid%nodata_value,i4)
!!$   write (999, 114) ( idCnew(k, l), l=1, grid%ncols)
!!$     end do
!!$     !write (999, 114)(int(grid%nodata_value,i4), l=1, grid%ncols+2)
!!$     close (999)
!!$   end do

  !104  format(2f20.2)
  !106  format(a4)
  !108  format(2i10)
!!$114  format(175i10)
!!$     !
!!$
!!$
!!$  stop
  ! <<<
  !
  deallocate(vec, nC, cnoList, cno )
end subroutine ClusterEvolution

!-------------------------------------------------------
! SVAT statistics
!-------------------------------------------------------
subroutine ClusterStats
  use mo_kind, only                                   : i4, dp
  !  use numerical_libraries, only                       : SVIGN , DSVRGN, DEQTIL
  use InputOutput, only                                : aDA, aDD, TDM, DTMagEvol, DAreaEvol, &
                                                         nCells, grid, SMI, &
                                                         nClusters, nMonths, idCluster, &
                                                         shortCnoList, deltaArea, &
                                                         dASevol, nBasins, Basin_Id, SMI_thld
  implicit none
  !
  integer(i4)                                         :: ic
  integer(i4)                                         :: t, i
  integer(i4) , dimension(:), allocatable             :: counterA
  integer(i4) , dimension(:,:), allocatable           :: aDDG
  real(dp) , dimension(:,:), allocatable              :: mSev
  !
  allocate ( aDDG     (grid%nrows, grid%ncols) )
  allocate ( DAreaEvol(nMonths,nClusters) )
  allocate ( DTMagEvol(nMonths,nClusters) )
  allocate ( aDD      (nClusters)   )
  allocate ( aDA      (nClusters)   )
  allocate ( counterA (nClusters)   )
  allocate ( TDM      (nClusters)   )
  allocate ( dASevol  (nMonths,2,nBasins+1)   )       ! (area, severity) whole Germany => nBasins+1
  allocate ( mSev     (grid%nrows, grid%ncols) )      ! mean  severity
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
subroutine calSAD(iDur)
  use mo_kind, only                                    : i4, dp
  ! use numerical_libraries, only                        : DSVRGN, DEQTIL, SVIGP
  use mo_sort, only                                    : sort, sort_index
  use mo_percentile, only                              : percentile
  use InputOutput, only                                : grid, shortCnoList, deltaArea, SMI, &
                                                         nInterArea, nEvents, nClusters,  idCluster, &
                                                         SAD, SADperc, DAreaEvol, severity, nDsteps, &
                                                         nDurations, durList, eventId, eIdPerm, &
                                                         writeResultsCluster, nLargerEvents, &
                                                         nCells, cellArea, nMonths, nQProp, QProp
  implicit none
  !
  integer(i4), intent(in)                             :: iDur
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
     print*, nInterArea, 2, nLargerEvents
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
  allocate ( severity (grid%nrows, grid%ncols, nDsteps ) )
  severity = grid%nodata_value
  !
  ! estimate severities for a given duration over all time steps
  do d = 1, nDsteps
     ms = (d-1)*durList(iDur) + 1
     me = ms + durList(iDur) - 1
     if (me > nMonths) me = nMonths
     where ( SMI(:,:,1) /= grid%nodata_value )
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
  ! write SAD for a given duration + percentiles
  call writeResultsCluster(2, durList(iDur) )
  !
  if ( allocated(sevP) ) deallocate ( sevP )
end subroutine calSAD

!-------------------------------------------------------
! clustering in space
! assign a unique running number at time t
!-------------------------------------------------------
subroutine findClusters (t,iC,nCluster)
  use mo_kind, only                                    :  i4
  use InputOutput, only                                 :  grid, cellCoor, thCellClus, &
       SMIc, nCells

  implicit none
  integer(i4), intent(in)                               :: t                   !  time step
  integer(i4), intent(out)                              :: nCluster            !  number of clusters larger than a threshold
  integer(i4), dimension(grid%nrows, grid%ncols), intent(out)   :: iC
  integer(i4)                                           :: i, j, k, klu, klu1, kk
  integer(i4)                                           :: iul, idr, jul, jdr
  integer(i4)                                           :: krow, kcol
  integer(i4), dimension(:), allocatable                :: cno, vec
  integer(i4), dimension(:), allocatable                :: nCxCluster
  integer(i4)                                           :: nClusterR, ncc
  !
  iC = int(grid%nodata_value, i4)
  nCluster = 0
  do k=1,nCells
     krow = cellCoor(k,1)
     kcol = cellCoor(k,2)
     iul  = max(krow-1,1)
     idr  = min(krow+1,grid%nrows)
     jul  = kcol
     jdr  = min(kcol+1,grid%ncols)
     ! SMIc of k
     klu  = SMIc(krow, kcol, t)
     if (klu /= 1) cycle
     if ( (klu == 1) .and. (iC(krow, kcol) == int(grid%nodata_value,i4)) ) then
        nCluster=nCluster+1
        iC(krow, kcol) = nCluster
     end if
     do j=jul, jdr
        do i= iul, idr
           if (iC(i,j) == int(grid%nodata_value,i4)  .and. &
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
        idr  = min(krow+1,grid%nrows)
        jul  = kcol
        jdr  = min(kcol+1,grid%ncols)
        do j=jul, jdr
           do i= iul, idr
              klu1 = iC(i, j)
              if ( klu  /= int(grid%nodata_value,i4) .and. &
                   klu1 /= int(grid%nodata_value,i4) .and. &
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
           where (iC == cno(i)) iC = int(grid%nodata_value, i4)
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

