
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: RootRad_GridComp - Implements Interface to a minimalistic MAPL GC that
!                   serves as a parent of the MAPL_ExtData GC
!
! !INTERFACE:
!
#include "MAPL_Generic.h"

   MODULE GEOS_RadiationEnvGridCompMod
!
! !USES:
!
   USE ESMF
   USE MAPL

   IMPLICIT NONE
   PRIVATE

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

  real, parameter :: P00          = MAPL_P00

! !DESCRIPTION: This is a Cinderella gridded component (GC)
!EOP

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


    type (ESMF_Config          )            :: CF
    integer                    :: i
    ! variables for alternative radiative coupling
    ! option 1
    integer, parameter         :: numRadCouple1_import = 4
    integer, parameter         :: numRadCouple1_export = 4
    character(len=ESMF_MAXSTR) :: RadCouple1_import(numRadCouple1_import) = (/"QLLS","QLCN","CLLS","CLCN"/)
    character(len=ESMF_MAXSTR) :: RadCouple1_export(numRadCouple1_export) = (/"QV  ","QL  ","QI  ","FCLD"/)
    character(len=ESMF_MAXSTR) :: RadCouple
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

! Set the state variable specs.
! -----------------------------

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'air_temperature',                                       &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'T',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
                                                        RC=STATUS  )

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME ='QILS',                                       &
        LONG_NAME  ='mass_fraction_of_large_scale_cloud_ice_water', &
        UNITS      ='kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME ='QICN',                                       &
        LONG_NAME  ='mass_fraction_of_convective_cloud_ice_water', &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'PL',                                        &
         LONG_NAME  = 'mid_level_pressure',                        &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PREF',                                      &
         LONG_NAME  = 'reference_air_pressure',                    &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='RL',                                          &
         LONG_NAME ='liquid_cloud_particle_effective_radius',      &
         UNITS     ='m',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RI',                                          &
         LONG_NAME = 'ice_phase_cloud_particle_effective_radius',   &
         UNITS     = 'm',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

!    Now have alternative Radiation coupling
     call ESMF_ConfigGetAttribute(CF, RadCouple, Default="USE_MODEL_RADCOUPLE_EXPORTS", &
                      Label="RADCOUPLE_METHOD:", __RC__ )

     if (trim(RadCouple) == "NEW_METHOD_1") then
        do i=1,numRadCouple1_import
           if (mapl_am_i_root()) write(*,*)'adding import ',trim(RadCouple1_import(i))
           call MAPL_AddImportSpec(GC,             &
                LONG_NAME  = 'UNKNOWN',            &
                UNITS      = 'UNKNOWN',            &
                SHORT_NAME = trim(RadCouple1_import(i)),   &
                DIMS       = MAPL_DimsHorzVert,    &
                VLOCATION  = MAPL_VLocationCenter, &
                RC         = STATUS )
        end do
        do i=1,numRadCouple1_export
           if (mapl_am_i_root()) write(*,*)'adding export ',trim(RadCouple1_export(i))
           call MAPL_AddExportSpec(GC,             &
                LONG_NAME  = 'UNKNOWN',            &
                UNITS      = 'UNKNOWN',            &
                SHORT_NAME = trim(RadCouple1_export(i)),   &
                DIMS       = MAPL_DimsHorzVert,    &
                VLOCATION  = MAPL_VLocationCenter, &
                RC         = STATUS )
        end do
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

! Local derived type aliases

    type (MAPL_MetaComp), pointer           :: MAPL
    type(ESMF_State)                        :: INTERNAL
    type(ESMF_FieldBundle)                  :: AEROSOLS
    type(ESMF_Grid)                         :: grid
    type(ESMF_Field)                        :: field
    integer                                 :: i
    integer                                 :: IM, JM, LM
    real, pointer                           :: ptr3(:,:,:) => null()

    real, pointer :: pref(:) => null()
    real :: a72(73)
    real :: b72(73)

      data a72 / &
       1.0000000,       2.0000002,       3.2700005,       4.7585009,       6.6000011, &
       8.9345014,       11.970302,       15.949503,       21.134903,       27.852606, &
       36.504108,       47.580610,       61.677911,       79.513413,       101.94402, &
       130.05102,       165.07903,       208.49704,       262.02105,       327.64307, &
       407.65710,       504.68010,       621.68012,       761.98417,       929.29420, &
       1127.6902,       1364.3402,       1645.7103,       1979.1604,       2373.0405, &
       2836.7806,       3381.0007,       4017.5409,       4764.3911,       5638.7912, &
       6660.3412,       7851.2316,       9236.5722,       10866.302,       12783.703, &
       15039.303,       17693.003,       20119.201,       21686.501,       22436.301, &
       22389.800,       21877.598,       21214.998,       20325.898,       19309.696, &
       18161.897,       16960.896,       15625.996,       14290.995,       12869.594, &
       11895.862,       10918.171,       9936.5219,       8909.9925,       7883.4220, &
       7062.1982,       6436.2637,       5805.3211,       5169.6110,       4533.9010, &
       3898.2009,       3257.0809,       2609.2006,       1961.3106,       1313.4804, &
       659.37527,       4.8048257,       0.0000000 /

      data b72 / &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,   8.1754130e-09,    0.0069600246,     0.028010041,     0.063720063, &
      0.11360208,      0.15622409,      0.20035011,      0.24674112,      0.29440312, &
      0.34338113,      0.39289115,      0.44374018,      0.49459020,      0.54630418, &
      0.58104151,      0.61581843,      0.65063492,      0.68589990,      0.72116594, &
      0.74937819,      0.77063753,      0.79194696,      0.81330397,      0.83466097, &
      0.85601798,      0.87742898,      0.89890800,      0.92038701,      0.94186501, &
      0.96340602,      0.98495195,       1.0000000 /

!  Get my name and set-up traceback handle
!  ---------------------------------------
   Iam = "Initialize_"
   call ESMF_GridCompGet( GC, name=comp_name, config=CF, rc=status )
   VERIFY_(STATUS)
   Iam = trim(comp_name) // '::' // trim(Iam)

!  Initialize MAPL Generic
!  -----------------------
   call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, clock,  __RC__ )

!  Retrieve the pointer to the state
!  ---------------------------------

   call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
   VERIFY_(STATUS)

   call MAPL_Get(MAPL, IM=IM, JM=JM, LM=LM, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS)
   VERIFY_(STATUS)

!  create PREF, ony needs to be done once
!  --------------------------------------
   call MAPL_GetPointer(EXPORT,PREF,'PREF',ALLOC=.true.,RC=STATUS)
   VERIFY_(STATUS)
   PREF = a72 + b72 * P00

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

   type (MAPL_MetaComp),      pointer  :: MAPL
   type (ESMF_GridComp),      pointer  :: GCS(:)
   type (ESMF_State),         pointer  :: GIM(:)
   type (ESMF_State),         pointer  :: GEX(:)
   type (ESMF_State)                   :: INTERNAL
   type (ESMF_Config)                  :: CF
   type (ESMF_FieldBundle)             :: AEROSOLS
   character(len=ESMF_MAXSTR),pointer  :: GCNames(:) => null()
   integer                             :: I,J,L
   integer                             :: IM, JM, LM
   real, pointer                       :: pref(:) => null()
   real, pointer                       :: U(:,:,:)    => null()
   real, pointer                       :: QILS(:,:,:) => null()
   real, pointer                       :: QICN(:,:,:) => null()
   real, pointer                       :: RL(:,:,:)   => null()
   real, pointer                       :: RI(:,:,:)   => null()
   real, pointer                       :: PL(:,:,:)   => null()
   real, pointer                       :: T(:,:,:)    => null()
   real, pointer                       :: aeroImport(:,:,:)   => null()
   real, pointer                       :: aeroInternal(:,:,:) => null()
   real                                :: tempor, plt
   integer, parameter :: numAEROs = 15
   character(len=ESMF_MAXSTR) :: aeroNames(numAEROs) = [character(len=ESMF_MAXSTR) :: &
                                                         "du001", "du002", "du003", "du004", &
                                                         "du005", "ss001", "ss002", "ss003", &
                                                         "ss004", "ss005", "SO4"  , "BCphobic", &
                                                         "BCphilic", "OCphobic", "OCphilic"]
   character(len=ESMF_MAXSTR) :: RadCouple
   ! pointers for alternative radiation coupling
   real, pointer                       ::   QV(:,:,:) => null()
   real, pointer                       ::   QL(:,:,:) => null()
   real, pointer                       ::   QI(:,:,:) => null()
   real, pointer                       :: FCLD(:,:,:) => null()
   real, pointer                       :: QLLS(:,:,:) => null()
   real, pointer                       :: QLCN(:,:,:) => null()
   real, pointer                       :: CLLS(:,:,:) => null()
   real, pointer                       :: CLCN(:,:,:) => null()

!  Get my name and set-up traceback handle
!  ---------------------------------------
   Iam = "Run_"
   call ESMF_GridCompGet( GC, name=comp_name, rc=status )
   VERIFY_(STATUS)
   call ESMF_GridCompGet( GC, name=comp_name, config=CF, rc=status )
   VERIFY_(STATUS)
   Iam = trim(comp_name) // '::' // trim(Iam)

! Get my internal MAPL_Generic state
!-----------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_Get(MAPL, IM=IM, JM=JM, LM=LM, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS)
   VERIFY_(STATUS)

   call MAPL_GetPointer(IMPORT, U, 'U' ,rc=status)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, T, 'T' ,rc=status)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, PL, 'PL' ,rc=status)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, QICN, 'QICN' ,rc=status)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, QILS, 'QILS' ,rc=status)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT, PREF, 'PREF' ,rc=status)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT, RI, 'RI' ,rc=status)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT, RL, 'RL' ,rc=status)
   VERIFY_(STATUS)

   do i=1,IM
      do j=1,JM
         call set_tempor(lm,pref,U(i,j,:),tempor)
         do l=1,LM
            plt = pl(i,j,l)/100.0
            call reff(T(i,j,l),plt,tempor,qils(i,j,l),qicn(i,j,l),rl(i,j,l),ri(i,j,l))
         end do
      end do
   end do

   ! finally do alternative radiation coupling
   ! -----------------------------------------
   call ESMF_ConfigGetAttribute(CF, RadCouple, Default="USE_MODEL_RADCOUPLE_EXPORTS", &
                    Label="RADCOUPLE_METHOD:", __RC__ )

   if (trim(RadCouple) == "NEW_METHOD_1") then
      call MAPL_GetPointer(IMPORT, QLLS, 'QLLS' ,rc=status)
      VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, QLCN, 'QLCN' ,rc=status)
      VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, CLLS, 'CLLS' ,rc=status)
      VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, CLCN, 'CLCN' ,rc=status)
      VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QV, 'QV' ,rc=status)
      VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QL, 'QL' ,rc=status)
      VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QI, 'QI' ,rc=status)
      VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, FCLD, 'FCLD' ,rc=status)
      VERIFY_(STATUS)
      ! now do some work on these
   end if
!  All done
!  --------
   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Run_

   subroutine SET_TEMPOR(LM,PREF,U,TEMPOR)

      integer, intent(in) :: LM
      real, intent(in ) :: PREF(0:LM)  ! [Pa]
      real, intent(in ) :: U(LM)       ! [m/s]
      real, intent(out) :: TEMPOR

      integer :: levs925, l

      levs925  = max(1,count(PREF < 92500.))

      TEMPOR = 0.
      do l = levs925, LM
        if (U(l).gt.4.) TEMPOR = 1.
      end do

   end subroutine SET_TEMPOR

   subroutine REFF(       &
         TE, PL, TEMPOR,  &
         QCiLS, QCiAN,    &
         RAD_RL, RAD_RI)

      real, intent(in ) :: TE              ! [K]
      real, intent(in ) :: PL              ! [hPa]
      real, intent(in ) :: QCiAN, QCiLS      ! [same units for both]
      real, intent(in ) :: TEMPOR
      real, intent(out) :: RAD_RL, RAD_RI  ! [m]

      ! parameters
      real, parameter :: MIN_RI = 20.e-6
      real, parameter :: MAX_RI = 40.e-6
      real, parameter :: RI_ANV = 30.e-6

      ! locals
      real :: RAD_RI_AN
      !!! calculate Effective Radii !!!

      if (PL < 150. ) then
         RAD_RI = MAX_RI
      end if
      if (PL >= 150. ) then
         RAD_RI = MAX_RI*150./PL
      end if

      !! weigh in a separate R_ice for Anvil Ice according to
      !
      !       R_net_eff = (q_anv + q_ls) / ( q_anv/R_ice_anv + q_ls/R_ice_ls )
      !-------------------------------------------------------------------------
      RAD_RI_AN  =  RAD_RI

      if ( ( QCiLS + QCiAN ) > 0.0 ) then
         RAD_RI_AN  = ( QCiLS + QCiAN ) / ( (QCiLS/RAD_RI) + (QCiAN/RI_ANV) )
      end if

      RAD_RI = MIN( RAD_RI, RAD_RI_AN )

      RAD_RI = MAX( RAD_RI, MIN_RI )

      ! Implement ramps for gradual change in effective radius
      if (PL < 300. ) then
         RAD_RL = 21.e-6
      end if
      if (PL >= 300. ) then
         RAD_RL = 21.e-6*300./PL
      end if
      RAD_RL = MAX( RAD_RL, 10.e-6 )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Thicken low high lat clouds
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! NOTE: Due to how tempor is calculated, it is now calculated in the
      ! GridComp and passed into progno_cloud

      if ( PL .GE. 775.  .AND. TE .LE.  275. .AND. (tempor.eq.1.) ) then
         RAD_RL = max(min(-0.1 * PL + 87.5, 10.),5.)*1.e-6
      end if
      if ( PL .GE. 825.  .AND. TE .LE.  282. .AND. (tempor.eq.1.) ) then
         RAD_RL = max(0.71 * TE - 190.25, 5.)*1.e-6
      end if
      if ( PL .GE. 775.  .AND. PL .LT. 825. .AND. TE .LE.  282. .AND. TE .GT. 275. .AND. (tempor.eq.1.) ) then
         RAD_RL = min(-0.1*PL + 0.71 * TE - 107.75, 10.)*1.e-6
      end if
      if ( PL .GE. 825.  .AND. TE .LE.  275. .AND. (tempor.eq.1.) ) then
         RAD_RL = 5.*1.e-6
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Thin low tropical clouds
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( PL .GE. 950.  .AND. TE .GE.  285. ) then
         RAD_RL = min(2.2 * TE - 617., 21.)*1.e-6
      end if
      if ( PL .GE. 925.  .AND. TE .GE.  290. ) then
         RAD_RL = min(0.44 * PL - 397., 21.)*1.e-6
      end if
      if ( PL .GE. 925.  .AND. PL .LT. 950. .AND. TE .GT.  285. .AND. TE .LT. 290.) then
         RAD_RL = max(min(0.44*PL + 2.2 * TE - 1035., 21.),10.)*1.e-6
      end if
      if ( PL .GE. 950.  .AND. TE .GE.  290. ) then
         RAD_RL = 21.*1.e-6
      end if

   end subroutine REFF

   end module GEOS_RadiationEnvGridCompMod
