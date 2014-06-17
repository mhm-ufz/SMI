!******************************************************************************
!  SETTING netCDF  v1                                                     *
!                                                                             *
!  AUTHOR:  Luis Samaniego UFZ 2011                                           *
!  UPDATES:                                                                   *
!           created  Sa     22.01.2011     main structure v1                  *
!           update   Sa     30.01.2011                                        *
!-----------------------------------------------------------------------------*
!  LIMITATIONS:                                                               *
!           *  VARIABLE TYPES ALLOWED in v1 (most common in mHM, EDK, LCC)    *
!              FOR ARRAYS                                                     *
          !      NF90_INT                                                     *
!                NF90_FLOAT                                                   *
!                NF90_DOUBLE                                                  *
!           *  MAXIMUM ARRAY SIZE = 4                                         *
!******************************************************************************
subroutine set_netCDF_SMI
  use mo_kind, only                            : i4,dp
  use netCDF_varDef
  use InputOutput, only                        :  DataPathOut, grid, SMI_flag, lons, lats
  use netcdf
  implicit none
  !
  ! local variables  
  integer(i4)                                  :: i, j
  real(dp)                                     :: xc, yc
  ! 
  ! define output file
  fName_netCDF_Out = trim(DataPathOut)//'mSMI.nc'
  !
  ! define parameters               ! nr. dim types
  nDims  = 3                                               ! nr. dim types
  nVars  = 6                                               ! total nr. var to print
  nLats  = grid%nrows                                      ! latitude  => nRows
  nLons  = grid%ncols                                      ! longitude => nCols
  !
  ! allocate arrays
  allocate ( Dnc(nDims)  )
  allocate ( V(nVars)    )
  allocate ( yCoor(nLats) )
  allocate ( xCoor(nLons) )
  allocate ( rxCoor(nLons,nLats), ryCoor(nLons,nLats) )

  ! >>> NEW
  if (SMI_flag == 4 .or. SMI_flag == 5) then
     do i = 1, nLons
        xCoor(i) = real(i,sp)
     end do
     ! inverse for Panoply, ncview display
     do j = nLats,1,-1
        yCoor(j) =  real(j,sp) 
     end do
     rxCoor = lons
     ryCoor = lats
  else
     ! def northings and eastings arrays
     xCoor(1) =  grid%xllcorner + 0.5_sp*real(grid%cellsize,sp)
     do i = 2, nLons
        xCoor(i) =  xCoor(i-1) + real(grid%cellsize,sp)
     end do
     ! inverse for Panoply, ncview display
     yCoor(nLats) =  grid%yllcorner + 0.5_sp*real(grid%cellsize,sp)
     do j = nLats-1,1,-1
        yCoor(j) =  yCoor(j+1) + real(grid%cellsize,sp)
     end do
     ! find lat and lon (curvilinear orthogonal gridding, See Panoply ref)
     do i = 1, nLons
        do j = 1, nLats
           ! (x,y) -> (lon, lat)
           call CoorTransInv( real(xCoor(i),dp) , real(yCoor(j),dp), xc, yc)
           rxCoor(i,j) = real(xc,sp)
           ryCoor(i,j) = real(yc,sp)
        end do
     end do
  end if

  !
  ! define dimensions (coordinates) => corresponding variable must be created 
  i              = 1
  Dnc(i)%name      = 'easting'
  Dnc(i)%len       = nLons
  !
  i              = 2
  Dnc(i)%name      = 'northing'
  Dnc(i)%len       = nLats
  !
  i              = 3
  Dnc(i)%name      = 'time'
  Dnc(i)%len       = NF90_UNLIMITED
  !
  !
  ! DIMENSION VARIABLES
  i                =  1
  V(i)%name        =  Dnc(i)%name
  V(i)%xType       =  NF90_FLOAT
  V(i)%nLvls       =  0
  V(i)%nSubs       =  0
  V(i)%nDims       =  1
  V(i)%dimTypes    =  (/i,0,0,0,0/)
  ! printing
  V(i)%wFlag       =  .true.
  ! pointer      
  V(i)%G1_f        => xCoor
  ! attributes
  V(i)%nAtt          = 3
  !
  V(i)%att(1)%name   = 'units'
  V(i)%att(1)%xType  = NF90_CHAR
  V(i)%att(1)%nValues= 1
  V(i)%att(1)%values  = 'm'
  !
  V(i)%att(2)%name   = 'long_name'
  V(i)%att(2)%xType  = NF90_CHAR
  V(i)%att(2)%nValues= 1
  V(i)%att(2)%values = 'x-coordinate in cartesian coordinates GK4'
  !
  V(i)%att(3)%name   = 'valid_range'
  V(i)%att(3)%xType  = NF90_FLOAT
  V(i)%att(3)%nValues= 2
  write( V(i)%att(3)%values, '(2f15.2)')  xCoor(1), xCoor(nLons)
  !
  !
  i                =  2
  V(i)%name        =  Dnc(i)%name
  V(i)%xType       =  NF90_FLOAT
  V(i)%nLvls       =  0
  V(i)%nSubs       =  0
  V(i)%nDims       =  1
  V(i)%dimTypes    =  (/i,0,0,0,0/)
  ! printing
  V(i)%wFlag       =  .true.
  ! pointer      
  V(i)%G1_f        => yCoor
  ! attributes
  V(i)%nAtt          = 3
  !
  V(i)%att(1)%name   = 'units'
  V(i)%att(1)%xType  = NF90_CHAR
  V(i)%att(1)%nValues= 1
  V(i)%att(1)%values  = 'm'
  !
  V(i)%att(2)%name   = 'long_name'
  V(i)%att(2)%xType  = NF90_CHAR
  V(i)%att(2)%nValues= 1
  V(i)%att(2)%values = 'y-coordinate in cartesian coordinates GK4'
  !
  V(i)%att(3)%name   = 'valid_range'
  V(i)%att(3)%xType  = NF90_FLOAT
  V(i)%att(3)%nValues= 2
  write( V(i)%att(3)%values, '(2f15.2)')  yCoor(1), yCoor(nLats)
  !
  i                =  3
  V(i)%name        =  Dnc(i)%name
  V(i)%xType       =  NF90_INT
  V(i)%nLvls       =  0
  V(i)%nSubs       =  0
  V(i)%nDims       =  1
  V(i)%dimTypes    =  (/i,0,0,0,0/)
  ! printing
  V(i)%wFlag       =  .true.
  ! pointer      
  !                   during running time  
  ! attributes 
  V(i)%nAtt          = 2
  !
  V(i)%att(1)%name   = 'units'
  V(i)%att(1)%xType  = NF90_CHAR
  V(i)%att(1)%nValues= 1
  if (SMI_flag == 4) then
     V(i)%att(1)%values = 'months since 1970-12-31 0:0:0' ! COSMO
  else if (SMI_flag == 5) then
     V(i)%att(1)%values = 'months since 1990-12-31 0:0:0' ! WRF-NOAH
  else
     V(i)%att(1)%values = 'months since 1950-12-31 0:0:0' ! mHM
  end if
  !
  V(i)%att(2)%name   = 'long_name'
  V(i)%att(2)%xType  = NF90_CHAR
  V(i)%att(2)%nValues= 1
  V(i)%att(2)%values = 'time'
  !  
  ! FIELD VARIABLES
  i                =  4
  V(i)%name        =  'SMI'
  V(i)%xType       =  NF90_FLOAT
  V(i)%nLvls       =  1
  V(i)%nSubs       =  1
  V(i)%nDims       =  3
  V(i)%dimTypes    =  (/1,2,3,0,0/)
  ! printing
  V(i)%wFlag       =  .true.
  ! pointer      
  !                   during running time 
  ! attributes 
  V(i)%nAtt          = 7
  !
  V(i)%att(1)%name   = 'units'
  V(i)%att(1)%xType  = NF90_CHAR
  V(i)%att(1)%nValues= 1
  V(i)%att(1)%values = '-'
  !
  V(i)%att(2)%name   = 'long_name'
  V(i)%att(2)%xType  = NF90_CHAR
  V(i)%att(2)%nValues= 1
  V(i)%att(2)%values = 'montly soil moisture index'
  !
  V(i)%att(3)%name   = 'valid_range'
  V(i)%att(3)%xType  = NF90_FLOAT
  V(i)%att(3)%nValues= 2
  V(i)%att(3)%values = '0.   1.'
  !
  V(i)%att(4)%name   = 'scale_factor'
  V(i)%att(4)%xType  = NF90_FLOAT
  V(i)%att(4)%nValues= 1
  V(i)%att(4)%values = '1.'
  !
  V(i)%att(5)%name   = '_FillValue'
  V(i)%att(5)%xType  = NF90_FLOAT
  V(i)%att(5)%nValues= 1
  write( V(i)%att(5)%values, * ) grid%nodata_value
  !
  V(i)%att(6)%name   = 'missing_value'
  V(i)%att(6)%xType  = NF90_FLOAT
  V(i)%att(6)%nValues= 1
  write( V(i)%att(6)%values, * ) grid%nodata_value   
  !
  V(i)%att(7)%name   = 'coordinates'
  V(i)%att(7)%xType  = NF90_CHAR
  V(i)%att(7)%nValues= 1
  V(i)%att(7)%values = 'lon lat'
  !
  i                =  5
  V(i)%name        =  'lon'
  V(i)%xType       =  NF90_FLOAT
  V(i)%nLvls       =  1
  V(i)%nSubs       =  1
  V(i)%nDims       =  2
  V(i)%dimTypes    =  (/1,2,0,0,0/)
  ! printing
  V(i)%wFlag       =  .true.
  ! pointer
  V(i)%G2_f        => rxCoor(:,:)
  ! attributes
  V(i)%nAtt          = 3
  !
  V(i)%att(1)%name   = 'units'
  V(i)%att(1)%xType  = NF90_CHAR
  V(i)%att(1)%nValues= 1
  V(i)%att(1)%values = 'degrees_east'
  !
  V(i)%att(2)%name   = 'long_name'
  V(i)%att(2)%xType  = NF90_CHAR
  V(i)%att(2)%nValues= 1
  V(i)%att(2)%values = 'longitude'
  !
  V(i)%att(3)%name   = 'valid_range'
  V(i)%att(3)%xType  = NF90_FLOAT
  V(i)%att(3)%nValues= 2
  V(i)%att(3)%values = '4.5    16.'
  !
  i                  =  6
  V(i)%name        =  'lat'
  V(i)%xType       =  NF90_FLOAT
  V(i)%nLvls       =  1
  V(i)%nSubs       =  1
  V(i)%nDims       =  2
  V(i)%dimTypes    =  (/1,2,0,0,0/)
  ! printing
  V(i)%wFlag       =  .true.
  ! pointer
  V(i)%G2_f        => ryCoor(:,:)
  ! attributes
  V(i)%nAtt          = 3
  !
  V(i)%att(1)%name   = 'units'
  V(i)%att(1)%xType  = NF90_CHAR
  V(i)%att(1)%nValues= 1
  V(i)%att(1)%values = 'degrees_north'
  !
  V(i)%att(2)%name   = 'long_name'
  V(i)%att(2)%xType  = NF90_CHAR
  V(i)%att(2)%nValues= 1
  V(i)%att(2)%values = 'latitude'
  !
  V(i)%att(3)%name   = 'valid_range'
  V(i)%att(3)%xType  = NF90_FLOAT
  V(i)%att(3)%nValues= 2
  V(i)%att(3)%values = '46.5   55.5'
  !
  ! global attributes

  if (SMI_flag == 4) then
     globalAtt(1)%name    = 'title'
     globalAtt(1)%values   = 'Mean SMI over Germany 7 x 7 km'
     !
     globalAtt(2)%name    = 'history'
     globalAtt(2)%values   = 'COSMO_IMK, proceceed by L. Samaniego'
  else if (SMI_flag == 5) then
     globalAtt(1)%name    = 'title'
     globalAtt(1)%values   = 'Mean SMI over Germany 10 x 10 km'
     !
     globalAtt(2)%name    = 'history'
     globalAtt(2)%values   = 'CORDEX Europe WRF331 0.11 deg ERA-INTERIM, proceceed by L. Samaniego'
  else
     globalAtt(1)%name    = 'title'
     globalAtt(1)%values   = 'Mean SMI over Germany 4 x 4 km'
     !
     globalAtt(2)%name    = 'history'
     globalAtt(2)%values   = 'mHM v4.3_block, proceceed by L. Samaniego, R. Kumar & M. Zink'
end if
end subroutine  set_netCDF_SMI

