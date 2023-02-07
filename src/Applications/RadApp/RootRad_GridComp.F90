
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: RootRad_GridComp - Implements Interface to a minimalistic MAPL
!             GC that serves as a parent of radiation and pchem
!
! !INTERFACE:
!
#include "MAPL_Generic.h"

MODULE RootRad_GridCompMod
!
! !USES:
!
   USE ESMF
   USE MAPL_Mod

   USE GEOS_RadiationEnvGridCompMod, only : RadiationEnvSetServices => SetServices
   USE GEOS_RadiationGridCompMod,    only : RadiationSetServices    => SetServices
   use GEOS_ChemGridCompMod,         only : AChemSetServices        => SetServices

   IMPLICIT NONE
   PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC SetServices

   INTEGER :: RAD
   INTEGER :: RADENV
   INTEGER :: CHEM

contains

!=============================================================================
!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

   subroutine SetServices (GC, RC)

! !ARGUMENTS:

      implicit none
      type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
      integer,             intent(  OUT) :: RC  ! return code

! !DESCRIPTION: The SetServices registers the Radiation component.

!EOP
!=============================================================================

      type(ESMF_Config)          :: cf
      integer                    :: STATUS, aeroProvider
      CHARACTER(LEN=ESMF_MAXSTR) :: providerName, ratsName, RadCouple

!=============================================================================

      call ESMF_GridCompGet (GC, CONFIG=CF, _RC)

! Set the Initialize, Run, Finalize entry points
! ----------------------------------------------
      call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run, _RC)

! Determine RATS (Radiatively Active Tracers) provider ...
! For now PCHEM is the only option. It uses parameterized prediction from tabulated
! zonally-symmetric production and loss data, or, alternatively, specification from
! zonally-symmetric values. This data-based method avoids running the whole GCM. 
! ---------------------------------------------------------------------------------
      call ESMF_ConfigGetAttribute (CF, &
         ratsName, Default="PCHEM", Label="RATS_PROVIDER:", _RC)
      _ASSERT (ratsName == "PCHEM", "Only PCHEM allowed as RATs provider currently!")
      CHEM = MAPL_AddChild (GC, NAME='CHEMISTRY', SS=AChemSetServices, _RC)

! Register Radiation Child and its supporter RADENV
! -------------------------------------------------
      RADENV = MAPL_AddChild (GC, NAME='RADENV',    SS=RadiationEnvSetServices, _RC)
      RAD    = MAPL_AddChild (GC, NAME='RADIATION', SS=RadiationSetServices,    _RC)

! Get RATs via CHEM
! -----------------
      call MAPL_AddConnectivity (GC, &
         SHORT_NAME = [character(len=6):: 'OX','O3','CH4','N2O','CFC11','CFC12','HCFC22'], &
         DST_ID     = RAD,     &
         SRC_ID     = CHEM, _RC)

! Determine AERO provider
! -----------------------
      call ESMF_ConfigGetAttribute (CF, &
         providerName, Label="AERO_PROVIDER:", Default="GOCART2G.data", _RC)
      _ASSERT (providerName == "GOCART2G" .or. providerName == "GOCART2G.data", &
         "Unsupported Aerosol Provider!")

      ! GOCART AERO bundle also comes through Chemistry GC
      aeroProvider = CHEM

! Source AERO bundle
! ------------------
      call MAPL_AddConnectivity (GC,   &
         SHORT_NAME = (/'AERO'/),      &
         DST_ID     = RAD,             &
         SRC_ID     = aeroProvider, _RC)

! Terminate imports that are actually not required by GOCART2G[.data]
! -------------------------------------------------------------------
      call MAPL_TerminateImport (GC, &
         SHORT_NAME = [character(len=9):: 'FRLAND','FRLANDICE','FROCEAN','FRACI'],  &
         CHILD = aeroProvider, &
         RC = STATUS)

! Imports of Radiation calculated in RadEnv rather than thru ExtData
! ------------------------------------------------------------------
      call MAPL_AddConnectivity (GC, &
         SHORT_NAME  = [character(len=4):: 'PREF','RI','RL'], &
         DST_ID      =  RAD,         &
         SRC_ID      =  RADENV,   _RC)

! Alternative Radiation coupling
! ------------------------------
      call ESMF_ConfigGetAttribute (CF, &
         RadCouple, Default="USE_MODEL_RADCOUPLE_EXPORTS", Label="RADCOUPLE_METHOD:", _RC)
      if (trim(RadCouple) == "USE_MODEL_RADCOUPLE_EXPORTS") then
         if (MAPL_AM_I_ROOT()) then
            write(*,*) "Using default radiation coupling,", &
               "QV, QL, QI, FCLD will be provided by ExtData"
         end if
      else if (trim(RadCouple) == "NEW_METHOD_1") then
         if (MAPL_AM_I_ROOT()) then
            write(*,*) "Using new_method_1: QV,QL,QI,FCLD will be produced by RadEnv", &
               "from Q, QILS, QICN, QLLS, QLCN, CLLS, and CLCN"
         end if
         call MAPL_AddConnectivity (GC, &
            SHORT_NAME = [character(len=4):: 'QL','QI','FCLD','QV'], &
            DST_ID     = RAD,           &
            SRC_ID     = RADENV,     _RC)
      end if

! Generic Set Services
! --------------------
      call MAPL_GenericSetServices (GC, _RC)
      _RETURN(_SUCCESS)

   end subroutine SetServices

!=============================================================================
!BOP

! !IROUTINE: Run

! !INTERFACE:

   subroutine Run (GC, IMPORT, EXPORT, CLOCK, RC)

! !INPUT PARAMETERS:

      implicit none
      type (ESMF_Clock),    intent(INOUT) :: CLOCK   ! The clock

! !OUTPUT PARAMETERS:

      type (ESMF_GridComp), intent(INOUT) :: GC      ! Grid Component
      type (ESMF_State),    intent(INOUT) :: IMPORT  ! Import State
      type (ESMF_State),    intent(INOUT) :: EXPORT  ! Export State
      integer, optional,    intent(  OUT) :: RC      ! Return code

! !DESCRIPTION: Run radiation and its supporters.
!
! !REVISION HISTORY:
!
!  29 Sep 2011  Anton Darmenov   Cloned from the ExtData GC code
!  -- Feb 2023  Peter Norris     Tidy up, get working again with GOCART2G
!
!EOP
!=============================================================================

      type (MAPL_MetaComp),      pointer :: STATE
      type (ESMF_GridComp),      pointer :: GCS(:)
      type (ESMF_State),         pointer :: GIM(:)
      type (ESMF_State),         pointer :: GEX(:)
      character(len=ESMF_MAXSTR),pointer :: GCNames(:)
      integer :: STATUS, I

!=============================================================================

! Get my MAPL_Generic state
!--------------------------

      call MAPL_GetObjectFromGC (GC, STATE, _RC)

! Get the children`s states from the generic state
!-------------------------------------------------

      call MAPL_Get (STATE,             &
          GCS=GCS, GIM=GIM, GEX=GEX,    &
          GCNames=GCNames,           _RC)

! Run the Radiation Code and its supporters
! -----------------------------------------

      I=CHEM

      call ESMF_GridCompRun (GCS(I), &
         importState=GIM(I), exportState=GEX(I), clock=CLOCK, &
         phase=1, userRC=STATUS); _VERIFY(STATUS)
      call ESMF_GridCompRun (GCS(I), &
         importState=GIM(I), exportState=GEX(I), clock=CLOCK, &
         phase=2, userRC=STATUS); _VERIFY(STATUS)
      call MAPL_GenericRunCouplers (STATE, I, CLOCK, _RC)

      I=RADENV

      call ESMF_GridCompRun (GCS(I), &
         importState=GIM(I), exportState=GEX(I), clock=CLOCK, &
         userRC=STATUS); _VERIFY(STATUS)
      call MAPL_GenericRunCouplers (STATE, I, CLOCK, _RC)

      I=RAD

      call MAPL_TimerOn (STATE, GCNames(I))
      call ESMF_GridCompRun (GCS(I), &
         importState=GIM(I), exportState=GEX(I), clock=CLOCK, &
         userRC=STATUS); _VERIFY(STATUS)
      call MAPL_GenericRunCouplers (STATE, I, CLOCK, _RC)
      call MAPL_TimerOff(STATE, GCNames(I))

      _RETURN(_SUCCESS)

   end subroutine Run

end module RootRad_GridCompMod
