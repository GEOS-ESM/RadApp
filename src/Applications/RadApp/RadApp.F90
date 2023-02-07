#define I_AM_MAIN

#include "MAPL_Generic.h"

Program RadApp_Main

   use MPI
   use MAPL
   use RootRad_GridCompMod, only: ROOT_SetServices => SetServices
   implicit none

   character(len=*), parameter :: Iam="RadApp_Main"
   type (MAPL_Cap) :: cap
   type (MAPL_FlapCLI) :: cli
   type (MAPL_CapOptions) :: copts
   integer :: status

   cli = MAPL_FlapCLI (description = 'GEOS RadApp', authors = 'GMAO')
   copts = MAPL_CapOptions (cli)
   cap = MAPL_Cap ('RootRad', ROOT_SetServices, cap_options=copts)

   call cap%run (_RC)

end program RadApp_Main

