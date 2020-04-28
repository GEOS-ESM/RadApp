
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: RootRad_GridComp - Implements Interface to a minimalistic MAPL GC that
!                   serves as a parent of radiation and pchem
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
   USE GEOS_RadiationGridCompMod, only : RadiationSetServices => SetServices
   use GEOS_ChemGridCompMod,       only : AChemSetServices     => SetServices

   IMPLICIT NONE
   PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC SetServices

   INTEGER :: RAD
   INTEGER :: RADENV
   INTEGER :: CHEM

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer,             intent(  OUT) :: RC  ! return code

! !DESCRIPTION:  The SetServices registers the Radiation component

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

    type(ESMF_Config)          :: cf
    integer                    :: aeroProvider
    CHARACTER(LEN=ESMF_MAXSTR) :: providerName, ratsName, RadCouple

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "::" // Iam

!   Set the Initialize, Run, Finalize entry points
!   ----------------------------------------------
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  Initialize_, RC=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,   Run_, RC=status)
    VERIFY_(STATUS)

!   Register children
!   -----------------
    RADENV = MAPL_AddChild(GC, NAME='RADENV', SS=RadiationEnvSetServices, RC=STATUS) 
    VERIFY_(STATUS)
    RAD = MAPL_AddChild(GC, NAME='RADIATION', SS=RadiationSetServices, RC=STATUS)
    VERIFY_(STATUS)

!   Determine AERO and RATS provider
!   --------------------------------
    call ESMF_ConfigGetAttribute(CF, ratsName, Default="PCHEM", &
                               Label="RATS_PROVIDER:", __RC__ )
    if (trim(ratsName) == "PCHEM") then
       chem = MAPL_AddChild(GC, NAME=    'CHEMISTRY', SS=AChemSetServices, RC=STATUS)
       VERIFY_(STATUS)
    end if

!   Imports of Radiation calculated in RadEnv not to come thru ExtData
!   ------------------------------------------------------------------
    call MAPL_AddConnectivity ( GC,                               &
        SHORT_NAME  = (/'PREF','RI','RL'/),                       &
        DST_ID      =  RAD,                                       &
        SRC_ID      =  RADENV,                                    &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(CF, providerName, Default="GOCART.data", &
                               Label="AERO_PROVIDER:", __RC__ )

!   if provider of aerobundle is pchem then we must terminate imports for its contents
!   so radenv doesn't expect them
    if (providerName == "GOCART.data" .or. providerName == "GOCART") then
       ASSERT_(trim(ratsName) == "PCHEM")
       aeroProvider = CHEM
       call MAPL_TerminateImport(GC, &
            SHORT_NAME = (/'FRLAND','FRLANDICE','FROCEAN','FRACI'/),  &
            CHILD = aeroProvider, &
            RC = STATUS)
    else
       ASSERT_(.false.)
    end if
        

     call MAPL_AddConnectivity ( GC,                               &
         SHORT_NAME  = (/'AERO'/),                                 &
         DST_ID      =  RAD,                                       &
         SRC_ID      =  aeroProvider,                              &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

!    for now pchem is the rats provider, no alternative
!    --------------------------------------------------
     if (trim(ratsName) == "PCHEM") then
        CALL MAPL_AddConnectivity( GC, &
            SHORT_NAME  = (/'OX    ','O3    ','CH4   ','N2O   ', &
                            'CFC11 ','CFC12 ','HCFC22'       /), &
                           DST_ID=RAD, SRC_ID=CHEM, RC=STATUS    )
        VERIFY_(STATUS)
     end if

!    Now have alternative Radiation coupling
     call ESMF_ConfigGetAttribute(CF, RadCouple, Default="USE_MODEL_RADCOUPLE_EXPORTS", &
                      Label="RADCOUPLE_METHOD:", __RC__ )
     if (trim(RadCouple) == "USE_MODEL_RADCOUPLE_EXPORTS") then
        if (MAPL_AM_I_ROOT()) then
           write(*,*)"Using default radiation coupling, QV, QL, QI, FCLD will be provided by ExtData"
        end if
     else if (trim(RadCouple) == "NEW_METHOD_1") then
        if (MAPL_AM_I_ROOT()) then
           write(*,*)"Using new_method_1, QV,QL, QI, FCLD will be produced by RadEnv from Q,QILS,QICN, QLLS,QLCN,CLLS, and CLCN"
        end if
        call MAPL_AddConnectivity(GC, &
             SHORT_NAME = (/'QL  ','QI  ','FCLD','QV  '/), &
             DST_ID=RAD, SRC_ID=RADENV, RC=STATUS)
        VERIFY_(STATUS)
     end if


!   Generic Set Services
!   --------------------
    call MAPL_GenericSetServices ( GC, rc=status)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

    end subroutine SetServices

!BOP
!
! !IROUTINE:  Initialize_ --- Initialize RUT
!
! !INTERFACE:
!

   SUBROUTINE Initialize_ ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:

   implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: CLOCK     ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout) :: GC      ! Grid Component
   type(ESMF_State), intent(inout) :: IMPORT     ! Import State
   type(ESMF_State), intent(inout) :: EXPORT     ! Export State
   integer, intent(out)            :: rc         ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  29 Sep 2011  Anton Darmenov   Cloned from the ExtData GC code
!
!EOP
!-------------------------------------------------------------------------


   type(ESMF_Config)           :: CF          ! Universal Config 
   character(len=ESMF_MAXSTR)  :: Iam
   integer                     :: status
   character(len=ESMF_MAXSTR)  :: comp_name
   integer :: im_world,jm_world,nx,ny,comm
   type(MAPL_MetaComp),pointer :: MAPL
   type(ESMF_VM) :: vm

!  Get my name and set-up traceback handle
!  ---------------------------------------
   Iam = "Initialize_"
   call ESMF_GridCompGet( GC, name=comp_name, config=CF, vm=vm, rc=status )
   VERIFY_(STATUS)
   call ESMF_VMGet(vm,mpiCommunicator=comm,rc=status)
   VERIFY_(status)
   call MAPL_GetObjectFromGC(GC,mapl,rc=status)
   VERIFY_(status)
   Iam = trim(comp_name) // '::' // trim(Iam)

!  Create grid for this GC
!  ------------------------
   call MAPL_GridCreate  (GC, RC=status )
   VERIFY_(STATUS)
   call MAPL_GetResource(MAPL,im_world,label='IM:',rc=status)
   VERIFY_(STATUS)
   call MAPL_GetResource(MAPL,jm_world,label='JM:',rc=status)
   VERIFY_(STATUS)
   call MAPL_GetResource(MAPL,nx,label='NX:',rc=status)
   VERIFY_(STATUS)
   call MAPL_GetResource(MAPL,ny,label='NY:',rc=status)
   VERIFY_(STATUS)

!  Initialize MAPL Generic
!  -----------------------
   call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, clock, RC=status)
   VERIFY_(STATUS)

!  All done
!  --------
   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Initialize_

   SUBROUTINE Run_ ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: CLOCK     ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: GC     ! Grid Component
   type(ESMF_State), intent(inout) :: IMPORT     ! Import State
   type(ESMF_State), intent(inout) :: EXPORT     ! Export State
   integer, intent(out) ::  rc                   ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  29 Sep 2011  Anton Darmenov   Cloned from the ExtData GC code
!
!EOP
!-------------------------------------------------------------------------

   character(len=ESMF_MAXSTR)    :: Iam
   integer                       :: STATUS
   character(len=ESMF_MAXSTR)    :: comp_name

! Local derived type aliases

   type (MAPL_MetaComp),      pointer  :: STATE
   type (ESMF_GridComp),      pointer  :: GCS(:)
   type (ESMF_State),         pointer  :: GIM(:)
   type (ESMF_State),         pointer  :: GEX(:)
   type (ESMF_State)                   :: INTERNAL
   type (ESMF_Config)                  :: CF
   character(len=ESMF_MAXSTR),pointer  :: GCNames(:)
   integer                             :: I
   integer                             :: IM, JM, LM

!  Get my name and set-up traceback handle
!  ---------------------------------------
   Iam = "Run_"
   call ESMF_GridCompGet( GC, name=comp_name, rc=status )
   VERIFY_(STATUS)
   Iam = trim(comp_name) // '::' // trim(Iam)

! Get my internal MAPL_Generic state
!-----------------------------------

   call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
   VERIFY_(STATUS)

! Get the children`s states from the generic state
!-------------------------------------------------

   call MAPL_Get ( STATE,   &
       GCS=GCS, GIM=GIM, GEX=GEX,       &
       IM = IM, JM = JM, LM = LM,       &
       GCNames = GCNames,               &
       INTERNAL_ESMF_STATE = INTERNAL,  &
                               RC=STATUS )
   VERIFY_(STATUS)

   I=CHEM

   call ESMF_GridCompRun (GCS(I), importState=GIM(I), exportState=GEX(I), clock=CLOCK, userRC=STATUS ); VERIFY_(STATUS)
   call MAPL_GenericRunCouplers (STATE, I,        CLOCK,    RC=STATUS ); VERIFY_(STATUS)

   I=RADENV

   call ESMF_GridCompRun (GCS(I), importState=GIM(I), exportState=GEX(I), clock=CLOCK, userRC=STATUS ); VERIFY_(STATUS)
   call MAPL_GenericRunCouplers (STATE, I,        CLOCK,    RC=STATUS ); VERIFY_(STATUS)

   I=RAD

   call MAPL_TimerOn (STATE,GCNames(I))
   call ESMF_GridCompRun (GCS(I), importState=GIM(I), exportState=GEX(I), clock=CLOCK, userRC=STATUS ); VERIFY_(STATUS)
   call MAPL_GenericRunCouplers (STATE, I,        CLOCK,    RC=STATUS ); VERIFY_(STATUS)
   call MAPL_TimerOff(STATE,GCNames(I))

!  All done
!  --------
   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Run_

   end module RootRad_GridCompMod
