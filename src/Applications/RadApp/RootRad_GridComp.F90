
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

      type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
      integer,             intent(  OUT) :: RC  ! return code

! !DESCRIPTION: The SetServices registers the Radiation component

!EOP
!=============================================================================

      character(len=ESMF_MAXSTR) :: IAm
      integer                    :: STATUS
      character(len=ESMF_MAXSTR) :: COMP_NAME

      type(ESMF_Config)          :: cf
      integer                    :: aeroProvider
      CHARACTER(LEN=ESMF_MAXSTR) :: providerName, ratsName, RadCouple

!=============================================================================

! Get my name and set-up traceback handle
! ---------------------------------------

      Iam = 'SetServices'
      call ESMF_GridCompGet (GC, NAME=COMP_NAME, CONFIG=CF, _RC)
      Iam = trim(COMP_NAME) // "::" // Iam

! Set the Initialize, Run, Finalize entry points
! ----------------------------------------------
      call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run_, _RC)

! Register children
! -----------------
      RADENV = MAPL_AddChild (GC, NAME='RADENV',    SS=RadiationEnvSetServices, _RC)
      RAD    = MAPL_AddChild (GC, NAME='RADIATION', SS=RadiationSetServices,    _RC)

! Determine AERO and RATS provider
! --------------------------------
      call ESMF_ConfigGetAttribute (CF, &
         ratsName, Default="PCHEM", Label="RATS_PROVIDER:", _RC)
      if (trim(ratsName) == "PCHEM") then
         CHEM = MAPL_AddChild (GC, NAME='CHEMISTRY', SS=AChemSetServices, _RC)
      end if

! Imports of Radiation calculated in RadEnv not to come thru ExtData
! ------------------------------------------------------------------
      call MAPL_AddConnectivity (GC, &
         SHORT_NAME  = [character(len=4):: 'PREF','RI','RL'], &
         DST_ID      =  RAD,         &
         SRC_ID      =  RADENV,   _RC)

      call ESMF_ConfigGetAttribute (CF, &
         providerName, Label="AERO_PROVIDER:", Default="GOCART.data", _RC)

!   if provider of aerobundle is pchem then we must terminate imports for its contents
!?   so radenv doesn't expect them
      if (providerName == "GOCART.data" .or. providerName == "GOCART") then
         ASSERT_(trim(ratsName) == "PCHEM")
         aeroProvider = CHEM
         call MAPL_TerminateImport (GC, &
            SHORT_NAME = [character(len=9):: 'FRLAND','FRLANDICE','FROCEAN','FRACI'],  &
            CHILD = aeroProvider, &
            RC = STATUS)
      else
         ASSERT_(.false.)
      end if

      call MAPL_AddConnectivity (GC,   &
         SHORT_NAME = (/'AERO'/),      &
         DST_ID     = RAD,             &
         SRC_ID     = aeroProvider, _RC)

!    for now pchem is the rats provider, no alternative
!    --------------------------------------------------
      if (trim(ratsName) == "PCHEM") then
         call MAPL_AddConnectivity (GC, &
            SHORT_NAME = [character(len=6):: 'OX','O3','CH4','N2O','CFC11','CFC12','HCFC22'], &
            DST_ID     = RAD,           &
            SRC_ID     = CHEM,       _RC)
      end if

!    Now have alternative Radiation coupling
      call ESMF_ConfigGetAttribute (CF, &
         RadCouple, Default="USE_MODEL_RADCOUPLE_EXPORTS", Label="RADCOUPLE_METHOD:", _RC)
      if (trim(RadCouple) == "USE_MODEL_RADCOUPLE_EXPORTS") then
         if (MAPL_AM_I_ROOT()) then
            write(*,*) "Using default radiation coupling,", &
               "QV, QL, QI, FCLD will be provided by ExtData"
         end if
      else if (trim(RadCouple) == "NEW_METHOD_1") then
         if (MAPL_AM_I_ROOT()) then
            write(*,*) "Using new_method_1, QV,QL, QI, FCLD will be produced by RadEnv", &
               "from Q,QILS,QICN, QLLS,QLCN,CLLS, and CLCN"
         end if
         call MAPL_AddConnectivity (GC, &
            SHORT_NAME = [character(len=4):: 'QL','QI','FCLD','QV'], &
            DST_ID     = RAD,           &
            SRC_ID     = RADENV,     _RC)
      end if

! Generic Set Services
! --------------------
      call MAPL_GenericSetServices (GC, _RC)

      RETURN_(ESMF_SUCCESS)

   end subroutine SetServices

!=============================================================================
!BOP

! !IROUTINE: Run_

! !INTERFACE:

   subroutine Run_ (GC, IMPORT, EXPORT, CLOCK, RC)

! !USES:

     implicit NONE

! !INPUT PARAMETERS:

      type (ESMF_Clock),    intent(INOUT) :: CLOCK   ! The clock

! !OUTPUT PARAMETERS:

      type (ESMF_GridComp), intent(INOUT) :: GC      ! Grid Component
      type (ESMF_State),    intent(INOUT) :: IMPORT  ! Import State
      type (ESMF_State),    intent(INOUT) :: EXPORT  ! Export State
      integer,              intent(  OUT) :: RC      ! Return code

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  29 Sep 2011  Anton Darmenov   Cloned from the ExtData GC code
!
!EOP
!=============================================================================

      character(len=ESMF_MAXSTR)    :: Iam
      integer                       :: STATUS
      character(len=ESMF_MAXSTR)    :: comp_name

! Local derived type aliases

      type (MAPL_MetaComp),      pointer  :: STATE
      type (ESMF_GridComp),      pointer  :: GCS(:)
      type (ESMF_State),         pointer  :: GIM(:)
      type (ESMF_State),         pointer  :: GEX(:)
      character(len=ESMF_MAXSTR),pointer  :: GCNames(:)
      integer                             :: I

!=============================================================================

!  Get my name and set-up traceback handle
!  ---------------------------------------
      Iam = "Run_"
      call ESMF_GridCompGet( GC, name=comp_name, _RC)
      Iam = trim(comp_name) // '::' // trim(Iam)

! Get my MAPL_Generic state
!-----------------------------------

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
         phase=1, userRC=STATUS); VERIFY_(STATUS)
      call ESMF_GridCompRun (GCS(I), &
         importState=GIM(I), exportState=GEX(I), clock=CLOCK, &
         phase=2, userRC=STATUS); VERIFY_(STATUS)
      call MAPL_GenericRunCouplers (STATE, I, CLOCK, _RC)

      I=RADENV

      call ESMF_GridCompRun (GCS(I), &
         importState=GIM(I), exportState=GEX(I), clock=CLOCK, &
         userRC=STATUS); VERIFY_(STATUS)
      call MAPL_GenericRunCouplers (STATE, I, CLOCK, _RC)

      I=RAD

      call MAPL_TimerOn (STATE, GCNames(I))
      call ESMF_GridCompRun (GCS(I), &
         importState=GIM(I), exportState=GEX(I), clock=CLOCK, &
         userRC=STATUS ); VERIFY_(STATUS)
      call MAPL_GenericRunCouplers (STATE, I, CLOCK, _RC)
      call MAPL_TimerOff(STATE, GCNames(I))

      RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Run_

end module RootRad_GridCompMod
