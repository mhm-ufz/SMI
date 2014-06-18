!****************************************************************************
!  clean allcatable arrays
!****************************************************************************
subroutine clean
  use InputOutput
  use kernelSmoother
  use SMIndex
  use netCDF_varDef
  implicit none
  !
  if (allocated ( SMI )   )  deallocate ( SMI )
  if (allocated ( X )     )  deallocate ( X )
  if (allocated ( edf )   )  deallocate ( edf )
  if (allocated ( pdf )   )  deallocate ( pdf )
  if (allocated ( SAD )  )  deallocate (SAD )
  if (allocated ( SADperc) )  deallocate (SADperc )
  if (allocated ( eventId ) )  deallocate ( eventId )
  if (allocated ( severity ) )  deallocate ( severity )
  if (allocated ( eIdPerm  ) )  deallocate ( eIdPerm  )
  if (allocated ( dASevol  ) )  deallocate ( dASevol  )
  if (allocated ( eMask    ) )  deallocate ( eMask )
  !
  print*, 'Memory deallocated ...'
  print*, 'Ended successfully ...'
end subroutine clean
