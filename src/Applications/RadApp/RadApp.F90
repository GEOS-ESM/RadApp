#define I_AM_MAIN

#include "MAPL_Generic.h"

Program RadApp_Main
   use MPI
   use MAPL
!!$   use GEOS_GcsGridCompMod, only:  ROOT_SetServices => SetServices
   use RootRad_GridCompMod, only:  ROOT_SetServices => SetServices
   implicit none
!EOP

!EOC

   character(len=*), parameter :: Iam="GEOS5_Main"
   type (MAPL_Cap) :: cap
   type (MAPL_FlapCLI) :: cli
   type (MAPL_CapOptions) :: cap_options
   integer :: status

   cli = MAPL_FlapCLI(description = 'GEOS RADAPP', &
                                     authors     = 'GMAO')
   cap_options = MAPL_CapOptions(cli)
   cap = MAPL_Cap('RootRad', ROOT_SetServices, cap_options = cap_options)

   call cap%run(_RC)

end program RadApp_Main
   
