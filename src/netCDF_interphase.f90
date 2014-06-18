!******************************************************************************
!  NetCDF_Interface v1                                                        *
!  PURPOSE:                                                                   *
!           This subroutines read / writes outputs in netCDF format           *
!           To become the standard input/output of mHM                        *
!                                                                             *
!  GRID     lon-lat-lvl grid                                                  *
!           time                                                              *
!                                                                             *
!  AUTHOR:  Luis Samaniego UFZ 2011                                           *
!                                                                             *
!  UPDATES:                                                                   *
!           created  Sa     22.01.2011     main structure v1                  *
!           update   Sa     24.02.2011     attribute rot. coordinate          *
!-----------------------------------------------------------------------------*
!  ACKNOWLEDGMENT                                                             *
!  VERSION: NetCDF-3.6.1 F90 interface to NetCDF                              *
!           http://www.unidata.ucar.edu/software/netcdf/                      *
!           CVF6.*        --Compaq Visual Fortran compiler version 6.*        *
!           IVF9.*, 10.*  --Intel  Visual Fortran compiler version 9.*, 10.*  *
!           netcdf90 compiled by:   YAN Haoming  <wfllib@gmail.com>           *
!                                   http://www.whigg.ac.cn/yanhm/fortran.htm  *
!  REQUIREMENTS AND SUGESTIONS                                                *
!           Matthias Cuntz UFZ                                                *
!******************************************************************************
module netCDF_varDef
  use mo_kind, only                          : i4, sp, dp
  ! netCDF configuration for mHM
  ! see: http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html
  !
  ! definition of parameters
  integer(i4), parameter                    :: nMaxDim = 5         ! nr. max dimensions
  integer(i4), parameter                    :: nMaxAtt = 20        ! nr. max attributes
  integer(i4), parameter                    :: maxLen  = 256       ! nr. string lenght
  integer(i4), parameter                    :: nGAtt   = 2         ! nr. global attributes
  integer(i4), parameter                    :: nAttDim = 2         ! dim array of attribute values
  !
  ! definition of variables
  integer(i4)                               :: nVars               ! nr. variables
  integer(i4)                               :: nDims               ! nr. dimensions 
  integer(i4)                               :: nLats               ! latitude  => nRows
  integer(i4)                               :: nLons               ! longitude => nCols
  integer(i4)                               :: nLvls               ! nr. levels
  integer(i4)                               :: nSubs               ! nr. of sub units / tiles/ fractions
  integer(i4)                               :: nRecs               ! time intervals
  !
  type attribute
    character (len=maxLen)                  :: name                ! attribute name
    integer(i4)                             :: xType               ! attribute of the values
    integer(i4)                             :: nValues             ! number of attributes       
    character (len=maxLen)                  :: values               ! numbers or "characters" separed by spaces
  end type attribute
  !
  type variable
     character (len=maxLen)                 :: name                ! short name
     integer(i4)                            :: xType               ! NF90 var. type
     integer(i4)                            :: nLvls               ! number of levels  
     integer(i4)                            :: nSubs               ! number of subparts     
     logical                                :: unlimited           ! time limited
     integer(i4)                            :: variD               ! Id
     integer(i4)                            :: nDims               ! field dimension
     integer(i4), dimension(nMaxDim)        :: dimIds              ! passing var. dimensions
     integer(i4), dimension(nMaxDim)        :: dimTypes            ! type of dimensions
     integer(i4)                            :: nAtt                ! nr. attributes     
     type(attribute), dimension(nMaxAtt)    :: att                 ! var. attributes   
     integer(i4), dimension(nMaxDim)        :: start               ! starting indices for netCDF
     integer(i4), dimension(nMaxDim)        :: count               ! counter          for netCDF
     logical                                :: wFlag               ! write flag
     integer(i4), dimension(:      ), pointer :: G1_i              ! array pointing model variables
     integer(i4), dimension(:,:    ), pointer :: G2_i              ! array pointing model variables
     integer(i4), dimension(:,:,:  ), pointer :: G3_i              ! array pointing model variables
     integer(i4), dimension(:,:,:,:), pointer :: G4_i              ! array pointing model variables
     real(sp), dimension(:    ), pointer    :: G1_f                ! array pointing model variables
     real(sp), dimension(:,:    ), pointer  :: G2_f                ! array pointing model variables
     real(sp), dimension(:,:,:  ), pointer  :: G3_f                ! array pointing model variables
     real(sp), dimension(:,:,:,:), pointer  :: G4_f                ! array pointing model variables
     real(dp), dimension(:      ), pointer  :: G1_d                ! array pointing model variables
     real(dp), dimension(:,:    ), pointer  :: G2_d                ! array pointing model variables
     real(dp), dimension(:,:,:  ), pointer  :: G3_d                ! array pointing model variables
     real(dp), dimension(:,:,:,:), pointer  :: G4_d                ! array pointing model variables
  end type variable
  type(variable), dimension(:), allocatable :: V                   ! variable list
  !
  type(attribute), dimension(nGAtt)         :: globalAtt           ! global attributes for netCDF
  !
  type dims
     character (len=maxLen)                 :: name                ! dim. name
     integer(i4)                            :: len                 ! dim. lenght, undefined time => NF90_UNLIMITED
     integer(i4)                            :: dimId               ! dim. Id
  end type dims
  type (dims), dimension(:), allocatable    :: Dnc                 ! dimesions netCDF e.g. lon, lat, level, time
  real(sp), dimension(:), allocatable, target      :: yCoor        ! lats
  real(sp), dimension(:), allocatable, target      :: xCoor        ! lons
  real(sp), dimension(:,:), allocatable, target    :: rxCoor        ! tranformed (rotated) coordinates
  real(sp), dimension(:,:), allocatable, target    :: ryCoor        ! tranformed (rotated) coordinates
  !
  ! input/output  
  character(256)                            :: fName_netCDF_In     ! File Name in
  character(256)                            :: fName_netCDF_Out    ! File Name out
  !
contains

!******************************************************************************
!  CREATE netCDF                                                              *
!******************************************************************************
subroutine create_netCDF(ncid)
  use mo_kind, only                          : i4, sp
  use netcdf
  implicit none
  !
  ! Only required for CVF
!     use typeSizes   
!     use netcdf90
!     implicit none
!     include 'netcdf.inc'       
  
  ! netcdf related variables
  integer(i4), intent(out)                   :: ncId
  integer(i4)                                :: i, j, k
  integer(i4), dimension (nAttDim)           :: att_INT  
  real(sp), dimension (nAttDim)              :: att_FLOAT
  character(len=maxLen), dimension (nAttDim) :: att_CHAR
  !
  ! 1  Create netCDF dataset: enter define mode     ->  get ncId  
  call check( nf90_create ( trim(fName_netCDF_Out), NF90_CLOBBER, ncId ))
  !
  !
  ! 2  Define dimensions                                 -> get dimId
  do i=1, nDims
    call check( nf90_def_dim ( ncId, Dnc(i)%name , Dnc(i)%len, Dnc(i)%dimId ))
  end do
  !
  ! 3 Define dimids array, which is used to pass the dimids of the dimensions of
  ! the netCDF variables
  do i = 1, nVars
    V(i)%unlimited = .false.
    V(i)%dimids  = 0
    V(i)%start   = 1
    V(i)%count   = 1
    do k = 1, V(i)%nDims
      if ( Dnc( V(i)%dimTypes(k) )%len  == NF90_UNLIMITED ) V(i)%unlimited = .true.
      V(i)%dimids(k) = Dnc( V(i)%dimTypes(k) )%dimId
    end do
    if ( V(i)%unlimited ) then
      ! set counts for unlimited files (time is the always the last dimension
      if (V(i)%nDims == 1) cycle
      do k = 1, V(i)%nDims - 1
         V(i)%count(k)  = Dnc( V(i)%dimTypes(k) )%len
      end do
    end if
  end do
  !
  ! 4 Define the netCDF variables and atributes                            -> get varId
  do i=1, nVars
    if ( .not. V(i)%wFlag ) cycle
    call check( nf90_def_var ( ncId, V(i)%name, V(i)%xtype, V(i)%dimids(1:V(i)%nDims), V(i)%varId ))
    do k = 1, V(i)%nAtt
      select case ( V(i)%att(k)%xType )
        case (NF90_CHAR)
          ! read( V(i)%att(k)%values, *) ( att_CHAR(j), j =1, V(i)%att(k)%nValues )
          read( V(i)%att(k)%values, '(a)')  att_CHAR(1) 
          call check( nf90_put_att ( ncId, V(i)%varId, V(i)%att(k)%name, att_CHAR(1) ))     
        case (NF90_INT)
          read( V(i)%att(k)%values, *) ( att_INT(j), j =1, V(i)%att(k)%nValues )
          call check( nf90_put_att ( ncId, V(i)%varId, V(i)%att(k)%name, att_INT(1:V(i)%att(k)%nValues) ))
        case (NF90_FLOAT)
          read( V(i)%att(k)%values, *) ( att_FLOAT(j), j =1, V(i)%att(k)%nValues )
          call check( nf90_put_att ( ncId, V(i)%varId, V(i)%att(k)%name, att_FLOAT(1:V(i)%att(k)%nValues) ))
      end select
    end do
  end do
  !
  ! 5 Global attributes
  do k = 1, nGAtt
    call check( nf90_put_att ( ncId, NF90_GLOBAL, globalAtt(k)%name, globalAtt(k)%values ))
  end do
  !
  ! 6 end definitions: leave define mode
  call check ( nf90_enddef ( ncId ))
  print *, 'NetCDF file was created', ncId
  ! 
end subroutine create_netCDF

!******************************************************************************
!  WRITE netCDF                                                               *
!                                                                             *
!  NOTES: 1) netCDF file must be on *** data mode ***                         *
!         2) any numeric data type is allowed but NOT character               *
!                                                                             *
!******************************************************************************
subroutine write_static_netCDF(ncId)
  use mo_kind, only                  : i4
  use netcdf
  implicit none
  !
  ! netcdf related variables
  integer(i4), intent(in)           :: ncId
  integer(i4)                       :: i
 
  ! NOTES: 1) netCDF file must be on *** data mode ***
  !
  ! write all static variables
  do i = 1, nVars
    if ( V(i)%unlimited     ) cycle
    if ( .not. V(i)%wFlag   ) cycle
    print *, 'Var. ', trim(V(i)%name) ,' is  static'
    select case (V(i)%xtype)
      case (NF90_INT)
        select case (V(i)%nDims)
          case (1)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G1_i ))
          case (2)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G2_i ))
          case (3)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G3_i ))
          case (4)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G4_i ))
        end select
      case (NF90_FLOAT)
        select case (V(i)%nDims)
          case (1)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G1_f ))
          case (2)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G2_f ))
          case (3)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G3_f ))
          case (4)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G4_f ))
        end select
      case (NF90_DOUBLE)
        select case (V(i)%nDims)
          case (1)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G1_d ))
          case (2)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G2_d ))
          case (3)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G3_d ))
          case (4)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G4_d ))
        end select
    end select
  end do

end subroutine write_static_netCDF

!******************************************************************************
!  WRITE netCDF                                                               *
!                                                                             *
!  NOTES: 1) netCDF file must be on *** data mode ***                         *
!         2) any numeric data type is allowed but NOT character               *
!                                                                             *
!******************************************************************************
subroutine write_dynamic_netCDF(ncId, irec)
 ! use mainVar
  use mo_kind, only                  : i4
  use netcdf
  implicit none
  !
  ! Only required for CVF
!     use typeSizes   
!     use netcdf90
!     implicit none
!     include 'netcdf.inc'       
  
  ! netcdf related variables
  integer(i4), intent(in)           :: ncId
  integer(i4), intent(in), optional :: iRec
  integer(i4)                       :: i
  ! NOTES: 1) netCDF file must be on *** data mode ***
  !        2) start and end of the data chuck is controled by
  !           V(:)%start and  V(:)%count
  !
  ! set values for variables (one scalar or grid at a time)
  do i = 1, nVars
    if ( .not. V(i)%unlimited     ) cycle
    if ( .not. V(i)%wFlag   ) cycle
    if (iRec ==1) print *, 'Var. ',i , trim(V(i)%name) ,' is  dynamic' !, iRec
    V(i)%start ( V(i)%nDims ) = iRec
    select case (V(i)%xtype)
      case (NF90_INT)
        select case (V(i)%nDims-1)
          case (0)
            call check( nf90_put_var ( ncId,  V(i)%varId, iRec, V(i)%start )) 
          case (1)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G1_i, V(i)%start, V(i)%count ))
          case (2)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G2_i, V(i)%start, V(i)%count ))
          case (3)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G3_i, V(i)%start, V(i)%count ))
          case (4)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G4_i, V(i)%start, V(i)%count ))
        end select
      case (NF90_FLOAT)
        select case (V(i)%nDims-1)
          case (1)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G1_f, V(i)%start, V(i)%count ))
          case (2)         
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G2_f, V(i)%start, V(i)%count ))
          case (3)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G3_f, V(i)%start, V(i)%count ))
          case (4)
            call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G4_f, V(i)%start, V(i)%count ))
         end select
       case (NF90_DOUBLE)
         select case (V(i)%nDims-1)
           case (1)
             call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G1_d, V(i)%start, V(i)%count ))
           case (2)
             call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G2_d, V(i)%start, V(i)%count ))
           case (3)
             call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G3_d, V(i)%start, V(i)%count ))
           case (4)
             call check( nf90_put_var ( ncId,  V(i)%varId, V(i)%G4_d, V(i)%start, V(i)%count ))
         end select
    end select
  end do
      
end subroutine write_dynamic_netCDF

!******************************************************************************
!  CLOSE netCDF                                                               *
!******************************************************************************
subroutine close_netCDF(ncId)
  use mo_kind, only               : i4
  use netcdf
  implicit none
!
  !  Only required for CVF
!     use typeSizes   
!     use netcdf90
!     implicit none
!     include 'netcdf.inc'       
  ! 
  integer(i4), intent(in)        :: ncId
  !
  ! close: save new netCDF dataset
  call check( nf90_close( ncId ))
  deallocate ( Dnc )
  deallocate ( V   )
  deallocate ( xCoor )
  deallocate ( yCoor )
  deallocate ( rxCoor )
  deallocate ( ryCoor )
  print *, 'NetCDF file was saved'
  !
end subroutine close_netCDF

!******************************************************************************
!  ERROR STATUS                                                              *
!******************************************************************************
subroutine check(status)
  use netcdf
  implicit none
  ! Only required for CVF
!  use typeSizes  
!  use netcdf90
!  implicit none
!  include 'netcdf.inc'
  integer(i4), intent(in) :: status
  !
  if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     stop
  end if
end subroutine check  
!
end module netCDF_varDef
