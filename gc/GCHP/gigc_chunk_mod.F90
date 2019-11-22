!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               __          _______  ______       _____  _____                 !
!               \ \        / /  __ \|  ____|     / ____|/ ____|                !
!                \ \  /\  / /| |__) | |__ ______| |  __| |                     !
!                 \ \/  \/ / |  _  /|  __|______| | |_ | |                     !
!                  \  /\  /  | | \ \| |         | |__| | |____                 !
!                   \/  \/   |_|  \_\_|          \_____|\_____|                !
!                                                                              !
!--------- v1.0 FOR PAPER SUBMITTION - PERMANENT ARCHIVAL (20191121) ----------!
!
! WRF-GC: GEOS-Chem High Performance-powered Chemistry Add-On for WRF Model
! Developed by Haipeng Lin <hplin@g.harvard.edu>, Xu Feng <fengx7@pku.edu.cn>
!    January 2018, Peking University, Dept of Atmospheric and Oceanic Sciences
!    Correspondence to: Tzung-May Fu <fuzm@sustech.edu.cn>
!
! VERSION INFORMATION:
!    This copy of WRF-GC is permanently archived for submission of the WRF-GC, Pt. 1
!    one-way coupling paper.
!    It may not be the latest version. We always recommend using the latest version
!    of WRF-GC, as important fixes may have been incorporated.
!    Please visit http://wrf.geos-chem.org for information.
!
! COPYRIGHT STATEMENT:
!    Permission is hereby granted, free of charge, to any person obtaining a copy
!   of this software and associated documentation files (the "Software"), to 
!   use, copy, modify the Software, and to permit persons to whom the Software is
!   furnished to do so, subject to the following conditions:
!
!   - The above copyright notice and this permission notice shall be included in all
!   copies or substantial portions of the Software.
! 
!   - The Software, modified in part or in full may not be redistributed without
!   express permission from the copyright holder(s).
! 
!   Except as contained in this notice or in attribution, the name of the WRF-GC model
!   shall not be used as an endorsement for distributing modified copies of the
!   Software without prior written permission from the copyright holder(s).
! 
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!   SOFTWARE.
! 
!  WRF and the GEOS-Chem model, GCHP are (c) their original authors.
!
!------------------------------------------------------------------------------- 

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_chunk_mod
!
! !DESCRIPTION: Module GC\_CHUNK\_MOD is the module that contains the init,
!  run, and finalize methods for the ESMF interface to GEOS-Chem.
!\\
!\\
! !INTERFACE: 
!      
MODULE GIGC_Chunk_Mod
  USE GCHP_Utils
  USE Input_Opt_Mod,           ONLY : OptInput
  USE State_Chm_Mod,           ONLY : ChmState
  USE State_Diag_Mod,          ONLY : DgnState
  USE State_Met_Mod,           ONLY : MetState
  USE HCO_TYPES_MOD,           ONLY : ConfigObj

  USE GC_Stateful_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GIGC_Chunk_Init
  PUBLIC :: GIGC_Chunk_Run
  PUBLIC :: GIGC_Chunk_Final

!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: GIGC_Switch_Dims

!
! !PRIVATE TYPES:
!
  ! Derived type for chunk diagnostic output (for code validation)
  TYPE GC_DIAG
     LOGICAL                    :: DO_PRINT     ! Should we print out?
     INTEGER                    :: N_DIAG       ! # of diag quantities
     INTEGER                    :: COUNT        ! Counter for averaging
     CHARACTER(LEN=10), POINTER :: NAME(:)      ! Tracer names
     REAL*8,            POINTER :: TRACER(:,:)  ! Tracer concentrations
     CHARACTER(LEN=40)          :: FILENAME     ! File name for output
     INTEGER                    :: LUN          ! File unit # for output
  END TYPE GC_DIAG

  ! Derived type object for saving concentration diagnostics
  TYPE(GC_DIAG)                 :: DIAG_COL

!
! !PUBLIC TYPES:
!
  ! Derived type for chunk operator options (for GIGC_Chunk_Run), replacing previous
  ! phases.
  TYPE GIGC_Chunk_Operators
    logical                     :: Conv
    logical                     :: DryDep
    logical                     :: Emis
    logical                     :: Tend
    logical                     :: Turb
    logical                     :: Chem
    logical                     :: WetDep

    logical                     :: GCDiagn      ! Use GEOS-Chem History-based GCC/NetCDF Diagns?
  END TYPE GIGC_Chunk_Operators
  PUBLIC :: GIGC_Chunk_Operators

  INTEGER, PARAMETER :: KIND_R4 = selected_real_kind(3,25)
!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca & P. Le Sager - Chunkized & cleaned up.
!  09 Oct 2012 - R. Yantosca - Now pass am_I_Root to all routines
!  09 Oct 2012 - R. Yantosca - Added comments, cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_chunk_mod.F90
!  01 Nov 2012 - R. Yantosca - Now pass Input Options object to routines
!  15 Mar 2013 - R. Yantosca - Add routine GIGC_Cap_Tropopause_Prs
!  08 Mar 2018 - E. Lundgren - Move gigc_initialization_mod contents to 
!                              gigc_chunk_init now that LOC is much reduced
!  03 Aug 2018 - H.P. Lin    - Merge (temporarily) changes from GCHP working tree, v12.0.0-dev (as of today).
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_Switch_Dims
!
! !DESCRIPTION: Subroutine GIGC\_Switch\_Dims contains a set of very, very hackish
!  operations to attempt to "manage" the in-module "state" variables that should've
!  been stored inside a derived type object. This is a temporary shim until we can
!  coordinate with GCST regarding moving everything into derived type objects.
!  (This is also a good reference of which variables are stateful) (hplin, 6/4/18)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Switch_Dims( am_I_Root,                                   &
                               ID,         IM,          JM,         LM,     &
                               lonCtr,     latCtr,      lonEdge,    latEdge,&
                               Input_Opt,                                   &
                               State_Met,  State_Chm,   State_Diag,         &
                               CMN_Alloc,  GC_Alloc,                        &
                               RC )
!
! !USES:
!
    USE ErrCode_Mod

    USE CMN_SIZE_MOD
    USE CMN_O3_MOD
    USE CMN_FJX_MOD,    ONLY : ZPJ, ODMDUST, ODAER, ISOPOD, IRHARR, JVN_, NWVAA, NDUST, NAER
    USE CMN_FJX_MOD,    ONLY : L_, L1_, L2_, JVL_, JXL_, JXL1_, JXL2_, JTAUMX, N_

    USE VDIFF_PRE_Mod,  ONLY : Init_Vdiff_Pre
    USE GC_GRID_MOD
    USE GC_Environment_Mod, ONLY : GC_Init_Grid, GC_Init_Regridding
    USE PRESSURE_MOD,   ONLY : Cleanup_Pressure, Init_Pressure, GET_AP, GET_BP
    USE GRID_REGISTRY_MOD, ONLY : Init_Grid_Registry, Cleanup_Grid_Registry
    USE CARBON_MOD,     ONLY : INIT_CARBON, CLEANUP_CARBON
    USE AEROSOL_MOD,    ONLY : INIT_AEROSOL, CLEANUP_AEROSOL
    USE DIAG_OH_MOD,    ONLY : INIT_DIAG_OH, CLEANUP_DIAG_OH
    USE UCX_MOD,        ONLY : INIT_UCX, CLEANUP_UCX
    USE SULFATE_MOD,    ONLY : INIT_SULFATE, CLEANUP_SULFATE
    USE C2H6_MOD,       ONLY : INIT_C2H6, CLEANUP_C2H6
    USE HCOI_GC_Main_Mod, ONLY : HCOI_GC_Init, HCOI_GC_Final
    USE PBL_MIX_MOD,    ONLY : INIT_PBL_MIX, CLEANUP_PBL_MIX
    USE Regrid_A2A_Mod, ONLY : Cleanup_Map_A2A
    USE TOMS_MOD,       ONLY : Init_Toms, Cleanup_Toms

    USE HCO_Error_Mod,  ONLY : hp ! For hemco precision
    USE HCO_VertGrid_Mod, ONLY : HCO_VertGrid_Define, HCO_VertGrid_Cleanup

    USE HCO_INTERFACE_MOD, ONLY : HcoState, ExtState

    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState

    USE GC_Stateful_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN)    :: am_I_Root    ! Are we on the root CPU?
    INTEGER,            INTENT(IN)    :: ID           ! Domain identifier within this CPU
    INTEGER,            INTENT(IN)    :: IM           ! # lons, this CPU
    INTEGER,            INTENT(IN)    :: JM           ! # lats, this CPU
    INTEGER,            INTENT(IN)    :: LM           ! # levs, this CPU
    REAL(KIND_R4),      INTENT(IN)    :: lonCtr(:,:)  ! Lon centers [radians]
    REAL(KIND_R4),      INTENT(IN)    :: latCtr(:,:)  ! Lat centers [radians]
    REAL(KIND_R4),      INTENT(IN)    :: lonEdge(:,:) ! Lat centers [radians]
    REAL(KIND_R4),      INTENT(IN)    :: latEdge(:,:) ! Lat centers [radians]
    LOGICAL,            INTENT(IN)    :: CMN_Alloc    ! Do CMN allocations? If you skip GC_Allocate_All, use this
    LOGICAL,            INTENT(IN)    :: GC_Alloc     ! Do other allocations? Not in init

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),     INTENT(INOUT) :: Input_Opt
    TYPE(MetState),     INTENT(INOUT) :: State_Met
    TYPE(ChmState),     INTENT(INOUT) :: State_Chm
    TYPE(DgnState),     INTENT(INOUT) :: State_Diag
    ! Note that actually these have been mostly checked to be IN only -- especially Input_Opt.
    ! GIGC_Switch_Dims theoretically should NOT touch Input_Opt as it is shared between dims.

!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT)   :: RC
!
! !REMARKS:
!  This code is basically a set of hacks originally called GIGC\_Do\_Nested\_Hacks.
!  The name says it all - it attempts to "patch" the lack of derived type objects at a 
!  higher level than you're supposed to, but it enables in-CPU switching of array dimensions
!  without resorting to changing core GEOS-CHEM code. This should be later coordinated with
!  GCST to eventually remove all this using the "hacks" here as a reference to where should be
!  changed.
!
! !REVISION HISTORY:
!  4 Jun 2018 - H.P. Lin   - First crack.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, SAVE                     :: PREVIOUS_ID = -1
    INTEGER, SAVE                     :: PREVIOUS_GC_ID = -1
    LOGICAL                           :: HCO_ERROR
    LOGICAL, SAVE                     :: FIRST = .TRUE.
    LOGICAL                           :: INIT_ID

    REAL(hp)                          :: Ap(LM+1), Bp(LM+1)
    INTEGER                           :: L

    CALL GIGC_State_Get_Status(am_I_Root, ID, INIT_ID)

    ! If we are already in this domain we do not need to perform switching...
    if(PREVIOUS_GC_ID .eq. ID) then
        write(6, *) "Switching is not necessary as we are already"
        write(6, *) "in the right domain!"

        return
    endif

    write(6, *) "%%%%%%%%% GIGC IN-PET GRID SWITCHING HACKS %%%%%%%%%"
    write(6, *) "    (Not very) Proudly brought to you by hplin      "
    write(6, *) "Some debug information:"
    write(6, *) "ID:", ID
    write(6, *) "IM, JM, LM:", IM, JM, LM, " will be set as PAR"
    write(6, *) "Allocs: CMN", CMN_Alloc, "GC", GC_Alloc
    write(6, *) "RC at entry-point:", RC
    write(6, *) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

    !=================================================================
    ! A set of ugly hacks follow below (hplin, 1/6/18)
    !=================================================================

    ! -------- Headers/CMN_SIZE_MOD -------- !

    ! This needs to be performed before GC_Allocate_All & GIGC_Get_Options,
    ! since the latter touches DLON/DLAT and the former allocates them.
    call Cleanup_CMN_SIZE(am_I_Root, RC)

    ! The following is only necessary and replicated from INIT_CMN_SIZE,
    ! to allow in-PET switching. Run this in GIGC_Chunk_Run only.
    ! WARNING: Due to namespace conflicts, I_LO, J_LO/HI are NOT updated from here.
    ! Please start your indeces with 1... this works OK for WRF as we do conversion
    ! but this may not apply generically for other models, if your array indeces
    ! don't get unified to 1.
    !
    ! When you pass in State_Met, State_Chm, ... state variables, you should retrieve
    ! them from the respective PET using GC_Stateful_Mod, which retrieves the objects
    ! from the domain ID cache.
    IIPAR = IM
    JJPAR = JM
    LLPAR = LM

    ! Backwards compatibility
    IGLOB = IIPAR
    JGLOB = JJPAR
    LGLOB = LLPAR

    ! It is important to diagnose the size of chemistry grid.
    ! This will be re-diagnosed by GIGC_Chunk_Run later, but for now you must set
    ! it to a value that is > 0 or you won't get O2/N2 in the atmosphere...
    ! ...and that makes Fast-JX /really/ angry. (hplin, 8/21/18)
    LLTROP  = LM
    LLSTRAT = LM

    IF ( Input_Opt%LUCX ) THEN
      LLCHEM     = LLSTRAT
      LLCHEM_FIX = LLSTRAT
    ELSE
      LLCHEM     = LLTROP
      LLCHEM_FIX = LLTROP_FIX
    ENDIF

    ! You don't really need to carry about world size, do you G-C?
    IM_WORLD = IM
    JM_WORLD = JM
    LM_WORLD = LM

    ! Alllocate DLON, DLAT arrays here
    ! Do not use INIT_CMN_SIZE because some values (LO, HI, WORLD) are not present in this shim
    if(CMN_Alloc) then
      ALLOCATE( DLON( IIPAR, JJPAR, LLPAR ), STAT=RC )
      DLON = 0e+0_fpp

      ! Delta-latitudes [degrees]
      ALLOCATE( DLAT( IIPAR, JJPAR, LLPAR ), STAT=RC )
      DLAT = 0e+0_fpp
    endif

    ! -------- Headers/CMN_O3_MOD -------- !
    ! Deallocate arrays in CMN_O3_MOD, and allocate them using new dimensional sizes.
    ! This should really NOT be called in a higher-level function like this... but
    ! we are still coordinating with GCST regarding removing module-level state vars.
    call Cleanup_CMN_O3(am_I_Root, RC)
    if(CMN_Alloc) call Init_CMN_O3(am_I_Root, RC)

    ! -------- Headers/CMN_FJX_MOD -------- !
    ! Do not deallocate all CMN_FJX_MOD arrays as not all are IM,JM,LM dep'd.
    ! Some are even only read once: WVAA,RHAA,NRLAA,NCMAA,RDAA...
    if(allocated(ZPJ)) deallocate(ZPJ)
    if(allocated(ODMDUST)) deallocate(ODMDUST)
    if(allocated(ODAER)) deallocate(ODAER)
    if(allocated(ISOPOD)) deallocate(ISOPOD)
    if(allocated(IRHARR)) deallocate(IRHARR)
    if(CMN_Alloc) then
      allocate(ZPJ(LLPAR, JVN_, IIPAR, JJPAR))
      allocate(ODMDUST(IIPAR, JJPAR, LLPAR, NWVAA, NDUST))
      allocate(ODAER(IIPAR, JJPAR, LLPAR, NWVAA, NAER))
      allocate(ISOPOD(IIPAR, JJPAR, LLPAR, NWVAA))
      allocate(IRHARR(IIPAR, JJPAR, LLPAR))

      ! Also update quantities in CMN_FJX_MOD...
      ! L_, L1_, L2_, JVL_, JXL_, JXL1_, JXL2_, JTAUMX
      L_     = LLPAR    ! Number of CTM layers
      L1_    = L_+1     ! Number of CTM layer edges
      L2_    = L1_*2    ! Number of levels in FJX grid that
                        ! inc. both edges and mid-points
      JVL_   = LLPAR    ! Vertical levels for J-values

      JXL_   = LLPAR    ! Vertical levels for J-values computed within Fast-JX
      JXL1_  = JXL_+1   ! Vertical levels edges for J-values
      JXL2_  = 2*JXL_+2 ! Max # levels in the basic Fast-JX grid (mid-level)

      JTAUMX = ( N_ - 4*JXL_ ) / 2  ! Maximum number of divisions ( i.e., may
                                    ! not get to ATAUMN)
    endif

    !! ------------- HEMCO -------------- !!
    ! HEMCO's state is pointed by HcoState, ExtState. Depending on the domain identifier we need to check
    ! if we need to switch out HEMCO.
    if(PREVIOUS_ID .ne. -1) then
      ! Now also update HEMCO post-initialization (hplin, 2/26/19)
      call GIGC_State_Set_HCO(am_I_Root, PREVIOUS_ID, HcoState, RC)
      call GIGC_State_Set_HCOX(am_I_Root, PREVIOUS_ID, ExtState, RC)
    endif

    if(PREVIOUS_ID .ne. -1 .and. PREVIOUS_ID .ne. ID) then
      ! Switch out HEMCO for re-initialization, if not initialized previously
      if(.not. INIT_ID) then
        nullify(HcoState)
        nullify(ExtState)
        write(6,*) "HcoState for ID not previously init'd, assigning", ID
      else
      ! Or just retrieve it from the Stateful container...
        CALL GIGC_State_Get_HCO(am_I_Root, ID, HcoState, RC)
        CALL GIGC_State_Get_HCOX(am_I_Root, ID, ExtState, RC)
        write(6,*) "Switch HcoState from PREV->ID", PREVIOUS_ID, ID
        write(6,*) "Sanity check assoc", associated(HcoState), associated(ExtState)
        write(6,*) "Locations of Hco, Ext", loc(HcoState), loc(ExtState)
      endif
    else
      write(6,*) "No need to switch HcoState for ID", ID
    endif
    if(associated(HcoState)) then
      write(6,*) "Debug HcoState nSpc, NX, NY, NZ", HcoState%nSpc, HcoState%NX, HcoState%NY, HcoState%NZ
    endif

    ! -------- GeosUtil/GC_GRID_MOD -------- !
    call Cleanup_Grid(am_I_Root, RC)
    if(CMN_Alloc) then
      ! call Init_Grid(am_I_Root, Input_Opt, IM, JM, LM, RC)
      ! Now uses GC_Init_Grid (GC_Environment_Mod) in v12. We will replace this
      ! with a less "grid-dependent" version later after we talk to GCST...
      ! (hplin, 8/6/18)
      call GC_Init_Grid( am_I_Root, Input_Opt, RC )

      ! Set grid parameters...
      call SetGridFromCtrEdges(am_I_Root, IM, JM, lonCtr, latCtr, lonEdge, latEdge, RC)
      if(RC /= GC_SUCCESS) then
        write(6, *) "GIGC_Switch_Dims: Failed to SetGridFromCtrEdges! Abort"
        write(6, *) "lonCtr(s):"
        write(6, *) lonCtr
        write(6, *) "latCtr(s):"
        write(6, *) latCtr
        write(6, *) "lonEdge(s):"
        write(6, *) lonEdge
        write(6, *) "latEdge(s):"
        write(6, *) latEdge
        return
      endif

      ! The areas need to be updated on every call, because for integration with
      ! WRF nested domains in-PET switching of AREA_M2 and Met fields, even
      ! CMN_SIZE_Mod resolutions are necessary. (hplin, 5/12/18)
      !
      ! However, this should only execute if met fields are already available,
      ! which may not be the case -
      if(CMN_Alloc .and. GC_Alloc .and. associated(State_Met%AREA_M2)) then
        AREA_M2 = State_Met%AREA_M2
      endif

      ! Update the HEMCO Grid Status -- this is replicated from
      ! HCOI_GC_MAIN_MOD and should be merged into GEOS-Chem mainline, but right now
      ! we will patch this at a higher level because we are not sure if we are
      ! allowed to patch HEMCO (todo: email) (hplin, 8/20/18)
      if(GC_Alloc) then
        ! Convert Ap from [hPa] to [Pa]
        DO L = 1, LLPAR + 1
          Ap(L) = GET_AP(L) * 100_hp
          Bp(L) = GET_BP(L)
        ENDDO

        ! Update the HCO Vertical grid...
        call HCO_VertGrid_Cleanup(HcoState%Grid%zGrid)
        call HCO_VertGrid_Define(am_I_Root, HcoState%Config, zGrid = HcoState%Grid%zGrid, &
                                 nz = LM, Ap = Ap, Bp = Bp, RC = RC)

        ! Update HEMCO grid variables...
        HcoState%Grid%XMID%Val       => XMID   (:,:,1)
        HcoState%Grid%YMID%Val       => YMID   (:,:,1)
        HcoState%Grid%XEDGE%Val      => XEDGE  (:,:,1)
        HcoState%Grid%YEDGE%Val      => YEDGE  (:,:,1)
        HcoState%Grid%YSIN%Val       => YSIN   (:,:,1)
        HcoState%Grid%AREA_M2%Val    => AREA_M2(:,:,1)
      endif
    endif

    ! -------- GeosCore/gc_environment_mod -------- !
    ! only replicates features GC_Allocate_All because we want to shim away Set_Input_Opt !
    if(CMN_Alloc) then
      ! Fix xmap_r4r4 regridding issue in HEMCO (hplin, 8/10/18)
      call Cleanup_Map_A2A()
      call GC_Init_Regridding(am_I_Root, Input_Opt, RC)
    endif

    ! -------- GeosUtil/input_mod -------- !
    ! CT1 is a private variable so you cannot just "reset" it. This causes serious issues
    ! as we need to replicate "input_mod"'s remaining functionality in GC_Stateful_Mod,
    ! replacing the original GIGC_Get_Options.

    ! -------- GeosCore/pbl_mix_mod -------- !
    call CLEANUP_PBL_MIX()
    if(GC_Alloc) call INIT_PBL_MIX(am_I_Root, RC)

    ! -------- GeosCore/pressure_mod -------- !
    call Cleanup_Pressure()
    if(GC_Alloc) call INIT_PRESSURE(am_I_Root)

    ! -------- GeosUtil/Grid_Registry_mod -------- !
    ! Register the horizontal and vertical grid information so that
    ! the History component can use it for netCDF metadata
    ! call Cleanup_Grid_Registry( am_I_Root, RC )
    ! if(GC_Alloc) call Init_Grid_Registry( am_I_Root, RC )

    ! -------- GeosCore/toms_mod -------- !
    call Cleanup_Toms(am_I_Root, RC)
    if(GC_Alloc) call Init_Toms( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )

    ! -------- GeosCore/Carbon_Mod -------- !
    call CLEANUP_CARBON()
    if(GC_Alloc) call INIT_CARBON( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )

    ! -------- GeosCore/Aerosol_Mod -------- !
    call CLEANUP_AEROSOL()
    if(GC_Alloc) call INIT_AEROSOL( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )

    ! -------- GeosCore/diag_oh_mod -------- !
    call CLEANUP_DIAG_OH()
    if(GC_Alloc) call INIT_DIAG_OH(am_I_Root, Input_Opt, RC)

    ! -------- GeosCore/ucx_mod -------- !
    ! *** GeosCore Code Change Warning: We added a reset of NSFCMR = 0 in CLEANUP_UCX
    call CLEANUP_UCX(am_I_Root)
    if(GC_Alloc) call INIT_UCX(am_I_Root, Input_Opt, State_Chm, State_Diag)

    ! -------- GeosCore/sulfate_mod -------- !
    call CLEANUP_SULFATE()
    if(GC_Alloc) call INIT_SULFATE(am_I_Root, Input_Opt, State_Chm, State_Diag, RC)

    ! -------- GeosCore/c2h6_mod -------- !
    call CLEANUP_C2H6()
    if(GC_Alloc) call INIT_C2H6(am_I_Root, Input_Opt, RC)

    !
    ! Oh my, here it goes -
    !  We found love on an empty page
    !  Kill the stars above, trying to fight the fade
    !                 - Oh Wonder "Bigger than Love"
    !
    FIRST = .FALSE.
    PREVIOUS_ID = ID

    if(GC_Alloc) then
        PREVIOUS_GC_ID = ID
    endif

    write(6, *) "Exit RC: ", RC
    write(6, *) "%%%%%%%%%    FINISHED GIGC_Switch_Dims    %%%%%%%%%"
    
  END SUBROUTINE GIGC_Switch_Dims
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_chunk_init
!
! !DESCRIPTION: Subroutine GIGC\_CHUNK\_INIT is the init method for
!  the Grid-Independent GEOS-Chem (aka "GIGC").  This routine calls lower-
!  level routines to allocate arrays and read input files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Chunk_Init( am_I_Root, I_LO,      J_LO,       I_HI,      &
                              J_HI,      IM,        JM,         LM,        &
                              ID,                                          &
                              IM_WORLD,  JM_WORLD,  LM_WORLD,   nymdB,     &
                              nhmsB,     nymdE,     nhmsE,                 &
                              tsChem,    tsDyn,                            &
                              lonCtr,    latCtr,    lonEdge,    latEdge,   & 
                              myPET,                                       &
                              Input_Opt, State_Chm, State_Diag, State_Met, &
                              HcoConfig, RC,        MPI_COMM)
!
! !USES:
!
    USE Chemistry_Mod,           ONLY : Init_Chemistry
    USE CMN_Size_Mod,            ONLY : IIPAR, JJPAR, LLPAR, dLon, dLat
    USE DiagList_Mod
    USE Emissions_Mod,           ONLY : Emissions_Init
    USE Fast_JX_Mod,             ONLY : Init_FJX
    USE GC_Environment_Mod
    USE GC_Grid_Mod,             ONLY : SetGridFromCtrEdges
    USE HCO_Types_Mod,           ONLY : ConfigObj
    USE Input_Mod,               ONLY : Read_Input_File
    USE Input_Opt_Mod,           ONLY : OptInput, Set_Input_Opt
    USE Linoz_Mod,               ONLY : Linoz_Read
    USE PBL_Mix_Mod,             ONLY : Init_PBL_Mix
    USE PhysConstants,           ONLY : PI_180
    USE Pressure_Mod,            ONLY : Init_Pressure
    USE Roundoff_Mod,            ONLY : RoundOff
    USE State_Chm_Mod,           ONLY : ChmState
    USE State_Diag_Mod,          ONLY : DgnState
    USE State_Met_Mod,           ONLY : MetState
    USE Time_Mod,                ONLY : Set_Timesteps
    USE UCX_MOD,                 ONLY : INIT_UCX
    USE UnitConv_Mod,            ONLY : Convert_Spc_Units

    ! Added hplin 8/3/18
    USE Species_Mod,             ONLY : Species
    USE State_Chm_Mod,           ONLY : IND_
    USE ErrCode_Mod
    USE Error_Mod,               ONLY : Debug_Msg
    USE HCO_INTERFACE_MOD,       ONLY : HcoState, ExtState

    USE History_Mod,             ONLY : History_Init
    ! USE Grid_Registry_Mod,       ONLY : Init_Grid_Registry

    ! Stateful Module (hplin)
    USE GC_Stateful_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,            INTENT(IN)    :: I_LO        ! Min lon index, this CPU
    INTEGER,            INTENT(IN)    :: J_LO        ! Min lat index, this CPU
    INTEGER,            INTENT(IN)    :: I_HI        ! Max lon index, this CPU
    INTEGER,            INTENT(IN)    :: J_HI        ! Max lat index, this CPU
    INTEGER,            INTENT(IN)    :: IM          ! # lons, this CPU
    INTEGER,            INTENT(IN)    :: JM          ! # lats, this CPU
    INTEGER,            INTENT(IN)    :: LM          ! # levs, this CPU
    INTEGER,            INTENT(IN)    :: ID          ! Domain identifier within this CPU
    INTEGER,            INTENT(IN)    :: IM_WORLD    ! # lons, global grid
    INTEGER,            INTENT(IN)    :: JM_WORLD    ! # lats, global grid
    INTEGER,            INTENT(IN)    :: LM_WORLD    ! # levs, global grid
    INTEGER,            INTENT(IN)    :: myPET       ! Local PET
    INTEGER,            INTENT(IN)    :: nymdB       ! YYYYMMDD @ start of run
    INTEGER,            INTENT(IN)    :: nhmsB       ! hhmmss   @ start of run
    INTEGER,            INTENT(IN)    :: nymdE       ! YYYYMMDD @ end of run
    INTEGER,            INTENT(IN)    :: nhmsE       ! hhmmss   @ end of run
    REAL,               INTENT(IN)    :: tsChem      ! Chemistry timestep [s]
    REAL,               INTENT(IN)    :: tsDyn       ! Chemistry timestep [s]
    REAL(KIND_R4), INTENT(IN)    :: lonCtr(:,:)      ! Lon centers [radians]
    REAL(KIND_R4), INTENT(IN)    :: latCtr(:,:)      ! Lat centers [radians]
    REAL(KIND_R4), INTENT(IN)    :: lonEdge(:,:)     ! Lon edges   [radians]
    REAL(KIND_R4), INTENT(IN)    :: latEdge(:,:)     ! Lat edges   [radians]

    INTEGER, INTENT(IN)  :: MPI_COMM                 ! MPI Communicator #
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),     INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(ChmState),     INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState),     INTENT(INOUT) :: State_Diag  ! Diagnostics State object
    TYPE(MetState),     INTENT(INOUT) :: State_Met   ! Meteorology State object
    TYPE(ConfigObj),    POINTER       :: HcoConfig   ! HEMCO config obj 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Need to add better error checking
!
! !REVISION HISTORY: 
!  18 Jul 2011 - M. Long     - Initial Version
!  28 Mar 2012 - M. Long     - Rewrite per structure of BCC init interface
!  09 Oct 2012 - R. Yantosca - Added comments, cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  16 Oct 2012 - R. Yantosca - Renamed GC_MET argument to State_Met
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Chunk_Init
!  01 Nov 2012 - R. Yantosca - Now reference gigc_input_opt_mod.F90
!  01 Nov 2012 - R. Yantosca - Reordered arguments for clarit
!  28 Nov 2012 - M. Long     - Now pass lonCtr, latCtr, latEdg as arguments
!                              to routine GIGC_Init_Simulation
!  03 Dec 2012 - R. Yantosca - Now call Init_CMN_SIZE (in CMN_SIZE_mod.F)
!                              instead of GIGC_Init_Dimensions to initialize
!                              the size parameters.
!  03 Dec 2012 - R. Yantosca - Rename NI, NJ, NL to IM, JM, LM for clarity
!  03 Dec 2012 - R. Yantosca - Now pass I_LO, J_LO, I_HI, J_HI, IM_WORLD, 
!                              JM_WORLD, LM_WORLD via the arg list
!  05 Dec 2012 - R. Yantosca - Remove latEdg argument
!  06 Dec 2012 - R. Yantosca - Add nymdB, nhmsB, nymdB, nhmsB arguments,
!                              and remove nymd, nhms
!  06 Mar 2018 - E. Lundgren - Remove Set_Initial_MixRatios
!  08 Mar 2018 - E. Lundgren - Move gigc_initialized_mod code here and move
!                              dlat/dlon calculation to gc_init_grid;
!                              GC timesteps are now seconds;
!                              Call set_input_opt to initialize input_opt vars;
!                              Add error handling using MAPL Assert_;
!                              Rename Initialize_Geos_Grid to GC_Init_Grid;
!                              Now call GC_Allocate_All after input.geos read;
!                              Restructure grid init based on gcbe v11-02e;
!                              Remove all unused code and simplify comments
!  14 Jan 2019 - E. Lundgren - Read input.geos on all threads and remove broadcast
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                           :: STATUS
    TYPE(Species),      POINTER       :: ThisSpc     ! Species pointer for init bg values
    INTEGER                           :: IND         ! Current species index
    INTEGER                           :: I           ! Loop idx...
    INTEGER                           :: II, JJ, LL
    LOGICAL                           :: prtDebug
    LOGICAL,            SAVE          :: FIRST = .TRUE.

    ! A diagnostics list object for State_Diag_Mod (v12) provided by
    ! GIGC_State_Boot
    TYPE(DgnList)                     :: DiagList

    !=======================================================================
    ! GIGC_CHUNK_INIT begins here 
    !=======================================================================
    ! Assume success
    RC = GC_SUCCESS

    ! Create a splash page
    IF ( am_I_Root ) THEN 
       WRITE(6, *) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       WRITE(6, *) "%               GEOS-CHEM HIGH PERFORMANCE               %"
       WRITE(6, *) "%                                                        %"
       WRITE(6, *) "% This GEOS-Chem column code was based on work by the    %"
       WRITE(6, *) "% original GCHP authors.                                 %"
       WRITE(6, *) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       WRITE(6, *) "Simulation Start Date/Time:", nymdB, nhmsB
       WRITE(6, *) "Simulation End Date/Time:", nymdE, nhmsE
       WRITE(6, *) "Chemical Timestep (s):", tsChem
       WRITE(6, *) "Dynamics Timestep (s):", tsDyn
       WRITE(6, *) "%%% CONFIGURATION FOR DOMAIN", ID, "%%%"
       WRITE(6, *) "In-PET LON IDX:", I_LO, I_HI
       WRITE(6, *) "In-PET LAT IDX:", J_LO, J_HI
       WRITE(6, *) "IM, JM, LM:", IM, JM, LM
       WRITE(6, *) "%%% WORLD CONFIGURATION %%%"
       WRITE(6, *) "World IM, JM, LM:", IM_WORLD, JM_WORLD, LM_WORLD
    ENDIF

    !=======================================================================
    ! Initialize GEOS-Chem Stateful Module & Allocate Arrays (hplin, 6/12/18)
    !=======================================================================

    !-----------------------------------------------------------------------
    ! Read info from the "input.geos" file into the Input_Opt object
    ! on the root CPU.  MPI broadcast Input_Opt to non-root CPUs.
    ! Continue with non-root CPU setup.
    ! This is performed using GC_Stateful_Mod now (hplin, 6/12/18)
    !-----------------------------------------------------------------------
    call GIGC_State_Boot(am_I_Root = Am_I_Root,        & ! Are we on the root PET?
                         MPI_COMM  = MPI_COMM,         & ! MPI Communicator
                         RC        = RC)

    ! If input options are not read successfully, abort the simulation
    IF ( RC /= GC_SUCCESS ) THEN
      WRITE(6, *) "### GIGC_CHUNK_INIT/INIT_SIMULATION: fatal error in GIGC_State_Boot. stop."
      stop
    ENDIF

    IF ( prtDebug ) THEN
      CALL DEBUG_MSG( '### GIGC_INIT_SIMULATION: after GIGC_State_Boot' )
    ENDIF

    ! Get the Input Options from GIGC_State_Boot
    call GIGC_State_Get_Opt(am_I_Root, Input_Opt)

    ! Get the Diagnostics List from GIGC_State_Boot
    call GIGC_State_Get_DiagList(am_I_Root, DiagList)

    !=======================================================================
    ! Temporary support for in-PET grid switching: verify if we need to
    ! pre-update CMN_SIZE_MOD variables and de/reallocate ahead of time.
    !=======================================================================
    ! Now reallocations are required every time, because GIGC_State_Boot initializes a dummy
    ! sized I/J/L. (hplin, 6/12/18)
    call GIGC_Switch_Dims( am_I_Root, ID, IM, JM, LM,        &
                           lonCtr, latCtr, lonEdge, latEdge, &
                           Input_Opt, State_Met, State_Chm, State_Diag, .true., .false., RC )
    write(6, *) "GIGC_Chunk_Init: Received new grid IM, JM, LM =", IM, JM, LM

    !=======================================================================
    ! Initialize key GEOS-Chem sections
    !=======================================================================
    ! Determine if we have to print debug output
    ! prtDebug = ( Input_Opt%LPRT .and. am_I_Root )
    prtDebug = .true.

    ! Update Input_Opt with timing fields
    Input_Opt%NYMDb   = nymdB
    Input_Opt%NHMSb   = nhmsB
    Input_Opt%NYMDe   = nymdE
    Input_Opt%NHMSe   = nhmsE
    Input_Opt%TS_CHEM = INT( tsChem )   ! Chemistry timestep [sec]
    Input_Opt%TS_EMIS = INT( tsChem )   ! Chemistry timestep [sec]
    Input_Opt%TS_DYN  = INT( tsDyn  )   ! Dynamic   timestep [sec]
    Input_Opt%TS_CONV = INT( tsDyn  )   ! Dynamic   timestep [sec]

    Input_Opt%myCPU = myPET

    IF ( prtDebug ) THEN
      CALL DEBUG_MSG( '### GIGC_INIT_SIMULATION: after Input_Opt% setups' )
    ENDIF

    !========================================================================
    ! Compute the DLON and DLAT values. In GC-GEOS-5 this is a kludge but
    ! in WRF-GC the coupler passes lonEdge, latEdge data directly to GIGC_Chunk
    ! routines, so the result is rigorous. (hplin, 11/11/18)
    !========================================================================
    DO LL = 1, LLPAR
    DO JJ = 1, JJPAR
    DO II = 1, IIPAR
       ! Compute Delta-Longitudes [degrees]
       dLon(II,JJ,LL)    = RoundOff( (lonEdge(II+1, JJ) - lonEdge(II, JJ)) / PI_180, 4 )

       ! Compute Delta-Latitudes  [degrees]
       dLat(II,JJ,LL)    = RoundOff( (latEdge(II, JJ+1) - latEdge(II, JJ)) / PI_180, 4 )
    ENDDO
    ENDDO
    ENDDO

    ! Initialize horizontal grid parameters --
    ! CALL GC_Init_Grid( am_I_Root, Input_Opt, RC )
    ! Set grid based on passed mid-points
    ! CALL SetGridFromCtr( am_I_Root, IM,          JM, &
    !                      lonCtr,    latCtr,      RC      )
    ! >>>>>>> these are currently done by GIGC_Switch_Dims for now

    ! The grid is always initialized by GIGC_Switch_Dims which will switch
    ! grid dimension and information on-the-fly (hplin, 6/13/18)

    ! Set timesteps
    CALL Set_Timesteps( am_I_Root  = am_I_Root,                          &
                        Chemistry  = Input_Opt%TS_CHEM,                  &
                        Convection = Input_Opt%TS_CONV,                  &
                        Dynamics   = Input_Opt%TS_DYN,                   &
                        Emission   = Input_Opt%TS_EMIS,                  &
                        Radiation  = Input_Opt%TS_RAD,                   &
                        Unit_Conv  = MAX( Input_Opt%TS_DYN,              &
                                          Input_Opt%TS_CONV ),           &
                        Diagnos    = Input_Opt%TS_DIAG         )

    ! Initialize derived-type objects for met, chem, and diag
    CALL GC_Init_StateObj( am_I_Root, DiagList,   Input_Opt, &
                           State_Chm, State_Diag, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    IF ( prtDebug ) THEN
      CALL DEBUG_MSG( '### GIGC_INIT_SIMULATION: after GC_Init_All' )
    ENDIF

   ! Initialize other GEOS-Chem modules
    CALL GC_Init_Extra( am_I_Root, DiagList,   Input_Opt,    &
                        State_Chm, State_Diag, RC ) 
    IF ( RC /= GC_SUCCESS ) RETURN

    IF ( prtDebug ) THEN
      CALL DEBUG_MSG( '### GIGC_INIT_SIMULATION: after GC_Init_Extra' )
    ENDIF

    ! Initialize photolysis
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .OR.                     &
         Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       CALL Init_FJX( am_I_Root, Input_Opt, State_Chm, State_Diag, RC ) 
       IF ( RC /= GC_SUCCESS ) RETURN
          
          !### Debug
          IF ( prtDebug ) THEN
             CALL DEBUG_MSG( '### GIGC_INIT_SIMULATION: after INIT_FJX' )        
          ENDIF
    ENDIF

    ! Set State_Chm units
    State_Chm%Spc_Units = 'kg/kg dry'

    ! Initialize pressure module (set Ap & Bp)
    CALL Init_Pressure( am_I_Root )

    if(am_I_Root) then
        write(*,*) '# GIGC_CHUNK_MOD: after Init_Pressure'
    endif

    ! Register the horizontal and vertical grid information so that
    ! the History component can use it for netCDF metadata
    ! call Init_Grid_Registry( am_I_Root, RC )

    ! if(am_I_Root) then
    !     write(*,*) '# GIGC_CHUNK_MOD: after Init_Grid_Registry'
    ! endif

    ! Initialize PBL mixing module
    CALL Init_PBL_Mix( am_I_Root, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    if(am_I_Root) then
        write(*,*) '# GIGC_CHUNK_MOD: after Init_PBL_Mix'
    endif
    
    ! Initialize chemistry mechanism
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .OR. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       CALL Init_Chemistry( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    if(am_I_Root) then
        write(*,*) '# GIGC_CHUNK_MOD: after Init_Chemistry'
    endif

    ! Initialize HEMCO
    CALL EMISSIONS_INIT ( am_I_Root, Input_Opt, State_Met, State_Chm, RC, &
                          HcoConfig=HcoConfig )

    if(am_I_Root) then
        write(*,*) '# GIGC_CHUNK_MOD: after EMISSIONS_INIT'
    endif

    ! Stratosphere - can't be initialized without HEMCO because of STATE_PSC
    IF ( Input_Opt%LUCX ) THEN

       ! Initialize stratospheric routines
       CALL INIT_UCX( am_I_Root, Input_Opt, State_Chm, State_Diag )

    ENDIF

    ! Convert species units to internal state units (v/v dry)
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, 'v/v dry', RC )

    !=======================================================================
    ! Initialize simulation with background values.
    !=======================================================================
    ! These values are eventually replaced in State_Chm% by whatever external model
    ! you use to handle chemistry data (in case of WRF, it's the chem array)
    ! every time you run GIGC_Chunk_Run. GIGC in particular is NOT aware of
    ! any location-specific background values.
    !
    ! Ported from Includes_Before_Run.H (GCHP) (hplin, 4/26/2018)
    do I = 1, State_Chm%nSpecies
      ThisSpc => State_Chm%SpcData(I)%Info

      if(trim(ThisSpc%Name) == '') cycle
      IND = IND_(trim(ThisSpc%Name))
      if(IND < 0) cycle

      ! Initialize using background values from species database.
      call SET_BACKGROUND_CONC(am_I_Root, ThisSpc, State_Chm, State_Met, Input_Opt, IND, RC)

      ! Fix negatives, from chem_gridcompmod
      where(State_Chm%Species < 0.0e0)
        State_Chm%Species = 1.0e-36
      endwhere

      ThisSpc => NULL()
    enddo

    !=======================================================================
    ! %%% Replicate GEOS-Chem Classic (non-ESMF) functionality in GIGC %%%
    !
    ! Initializes the History component in GEOS-Chem.
    !=======================================================================
    ! For now, just hardwire the input file for the History component
    Input_Opt%HistoryInputFile = './HISTORY.rc'

    ! Initialize the GEOS-Chem history component
    ! CALL History_Init( am_I_root, Input_Opt, State_Met, State_Chm, State_Diag, RC )
    ! if(am_I_Root) then
    !     write(*,*) '# GIGC_CHUNK_MOD: after HISTORY_INIT'
    ! endif

    !=======================================================================
    ! Save stateful information into GC_Stateful_Mod (hplin, 6/12/18)
    !=======================================================================
    CALL GIGC_State_Set_Opt(am_I_Root, Input_Opt, HcoConfig)
    CALL GIGC_State_Init(am_I_Root    = am_I_Root,   &
                         ID           = ID,          &
                         State_Met    = State_Met,   &
                         State_Chm    = State_Chm,   &
                         State_Diag   = State_Diag,  &
                         HcoState     = HcoState,    &
                         ExtState     = ExtState,    &
                         RC           = RC           )

    write(6, *) '# GIGC_CHUNK_MOD: after SV save, ID =', ID

    FIRST = .FALSE.

  END SUBROUTINE GIGC_Chunk_Init

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_chunk_run
!
! !DESCRIPTION: Subroutine GIGC\_CHUNK\_RUN is the run method for
!  the Grid-Independent GEOS-Chem (aka "GIGC").  This routine is the driver
!  for the following operations:
!
! \begin{itemize}
! \item Dry deposition
! \item Chemistry
! \end{itemize}
!
! !INTERFACE:
!
  SUBROUTINE GIGC_Chunk_Run( am_I_Root,                                   &
                             IM,        JM,        LM,        ID,         &
                             nymd,      nhms,      year,      month,      &
                             day,       dayOfYr,   hour,      minute,     &
                             second,    hElapsed,                         &
                             lonCtr,    latCtr,    lonEdge,   latEdge,    &
                             Input_Opt, State_Chm, State_Met, State_Diag, &
                             Operators, IsChemTime,                       &
                             RC                                           )
!
! !USES:
!
    ! GEOS-Chem state objects 
    USE HCO_Interface_Mod,  ONLY : HcoState
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod
    USE State_Met_Mod,      ONLY : MetState

    ! GEOS-Chem components
    USE Aerosol_Mod,        ONLY : Set_AerMass_Diagnostic
    USE Chemistry_Mod,      ONLY : Do_Chemistry, Recompute_OD
    USE Convection_Mod,     ONLY : Do_Convection
    USE DryDep_Mod,         ONLY : Do_DryDep
    USE Emissions_Mod,      ONLY : Emissions_Run
    USE Mixing_Mod,         ONLY : Do_Tend, Do_Mixing
    USE Strat_Chem_Mod,     ONLY : Init_Strat_Chem, Minit_is_Set
    USE WetScav_Mod,        ONLY : Setup_WetScav, Do_WetDep

    ! Specialized subroutines
    USE Dao_Mod,            ONLY : AirQnt, Set_Dry_Surface_Pressure
    USE Dao_Mod,            ONLY : GIGC_Cap_Tropopause_Prs
    USE Set_Global_CH4_Mod, ONLY : Set_CH4
    USE PBL_Mix_Mod,        ONLY : Compute_PBL_Height
    USE Pressure_Mod,       ONLY : Set_Floating_Pressures
    USE TOMS_Mod,           ONLY : Compute_Overhead_O3
    USE UCX_Mod,            ONLY : Set_H2O_Trac

    ! Utilities
    USE ErrCode_Mod
    USE GC_Grid_Mod,        ONLY : AREA_M2
    USE HCO_Error_Mod
    USE HCO_Interface_Mod,  ONLY : SetHcoTime
    USE Pressure_Mod,       ONLY : Accept_External_Pedge
    USE State_Chm_Mod,      ONLY : IND_
    USE Time_Mod,           ONLY : Accept_External_Date_Time
    Use UnitConv_Mod,       ONLY : Convert_Spc_Units

    ! Diagnostics
    USE Diagnostics_Mod,    ONLY : Set_Diagnostics_EndofTimestep

    ! Added hplin 8/3/18
    USE CMN_Size_Mod
    USE Precision_Mod
    USE Error_Mod,          ONLY : Debug_Msg
    USE GC_Grid_Mod,        ONLY : SetGridFromCtrEdges
    USE History_Mod,        ONLY : History_Write, History_Update, History_SetTime

    ! Added hplin 12/1/18
    USE Olson_Landmap_Mod,  ONLY : Compute_Olson_Landmap_GCHP

!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)                 :: am_I_Root   ! Are we on root CPU?
    INTEGER,        INTENT(IN)                 :: IM          ! # of lons on this CPU
    INTEGER,        INTENT(IN)                 :: JM          ! # of lats on this CPU
    INTEGER,        INTENT(IN)                 :: LM          ! # of levs on this CPU
    INTEGER,        INTENT(IN)                 :: ID          ! Domain identifier within this CPU
    INTEGER,        INTENT(IN)                 :: nymd        ! YYYY/MM/DD @ current time
    INTEGER,        INTENT(IN)                 :: nhms        ! hh:mm:ss   @ current time
    INTEGER,        INTENT(IN)                 :: year        ! UTC year 
    INTEGER,        INTENT(IN)                 :: month       ! UTC month
    INTEGER,        INTENT(IN)                 :: day         ! UTC day
    INTEGER,        INTENT(IN)                 :: dayOfYr     ! UTC day of year
    INTEGER,        INTENT(IN)                 :: hour        ! UTC hour
    INTEGER,        INTENT(IN)                 :: minute      ! UTC minute
    INTEGER,        INTENT(IN)                 :: second      ! UTC second
    REAL(KIND_R4),  INTENT(IN)                 :: hElapsed    ! Hours elapsed in run
    TYPE(GIGC_Chunk_Operators),  INTENT(IN)    :: Operators   ! Operators to run (derived type)
    LOGICAL,        INTENT(IN)                 :: IsChemTime  ! Time for chemistry? 
    REAL(KIND_R4),  INTENT(IN)                 :: lonCtr(:,:) ! Lon centers [radians]
    REAL(KIND_R4),  INTENT(IN)                 :: latCtr(:,:) ! Lat centers [radians]
    REAL(KIND_R4), INTENT(IN)                  :: lonEdge(:,:)! Lon edges   [radians]
    REAL(KIND_R4), INTENT(IN)                  :: latEdge(:,:)! Lat edges   [radians]
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt   ! Input Options obj
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm   ! Chemistry State obj
    TYPE(MetState),      INTENT(INOUT) :: State_Met   ! Meteorology State obj
    TYPE(DgnState),      INTENT(INOUT) :: State_Diag  ! Diagnostics State obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!
! !REVISION HISTORY:
!  18 Jul 2011 - M. Long     - Initial Version
!  09 Oct 2012 - R. Yantosca - Added extra comments & cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  16 Oct 2012 - R. Yantosca - Renamed GC_MET argument to State_Met
!  17 Oct 2012 - R. Yantosca - Need to call AIRQNT before chemistry
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Chunk_Run
!  25 Oct 2012 - R. Yantosca - Now pass RC to GIGC_DO_CHEM
!  01 Nov 2012 - R. Yantosca - Now reference gigc_input_opt_mod.F90
!  08 Nov 2012 - R. Yantosca - Now pass Input_Opt to GIGC_Do_Chem
!  13 Nov 2012 - M. Long     - Added Dry Deposition method
!  29 Nov 2012 - R. Yantosca - Now block off calls to GIGC_DO_DRYDEP and
!                              GIGC_DO_CHEM w/ the appropriate logical flags
!  04 Dec 2012 - R. Yantosca - Now convert units of State_Chm%TRACERS here
!                              instead of in lower-level routines
!  07 Dec 2012 - R. Yantosca - Now call Accept_Date_Time_From_ESMF to pass the
!                              date & time from ESMF to GeosUtil/time_mod.F
!  07 Dec 2012 - R. Yantosca - Now pass UTC via Accept_Date_Time_From_ESMF;
!                              this ensures proper localtime computation
!  11 Dec 2012 - R. Yantosca - Now call DO_DRYDEP directly; no longer call
!                              GIGC_DO_DRYDEP, this is moved to obsolete dir.
!  11 Dec 2012 - R. Yantosca - Now call routine ACCEPT_EXTERNAL_PEDGE to pass
!                              the pressure edges from ESMF to GEOS-Chem
!  15 Mar 2013 - R. Yantosca - Now call GIGC_CAP_TROPOPAUSE_PRS to cap the
!                              State_Met%TROPP field to 200 hPa polewards
!                              of 60S and 60N.  We do this in the std G-C.
!  05 Jun 2013 - R. Yantosca - Remove obsolete code
!  22 Sep 2014 - C. Keller   - Added run phase argument
!  14 Oct 2014 - C. Keller   - Various updates to include drydep and emissions
!                              to tracer arrays, etc.
!  26 Nov 2014 - C. Keller   - Added IsChemTime variable.
!  19 Oct 2016 - R. Yantosca - Now call Set_Init_Cond_Strat_Chem after the
!                              1st call to AIRQNT to save initial conditions
!  01 Dec 2016 - E. Lundgren - Calculate LAI using new routine for GCHP
!  13 Feb 2018 - E. Lundgren - Call Recompute_OD at end of chem dt for aer diags
!EOP
!------------------------------------------------------------------------------
!BOC
    REAL*8                         :: DT
    INTEGER                        :: STATUS

    ! Local logicals to turn on/off individual components
    ! The parts to be executed are based on the input options,
    ! the time step and the phase.
    LOGICAL                        :: DoConv 
    LOGICAL                        :: DoDryDep
    LOGICAL                        :: DoEmis
    LOGICAL                        :: DoTend 
    LOGICAL                        :: DoTurb 
    LOGICAL                        :: DoChem
    LOGICAL                        :: DoWetDep

    ! First call?
    LOGICAL, SAVE                  :: FIRST = .TRUE.

    ! Update mixing ratios during AirQnt due to pressure change?
    LOGICAL                        :: pUpdate

    ! # of times this routine has been called. Only temporary for printing 
    ! processes on the first 10 calls.
    INTEGER, SAVE                  :: NCALLS = 0

    ! Loop variable for going through levs
    INTEGER                        :: L = 1
    INTEGER                        :: N

    ! GMT time in decimal hours
    REAL*4                         :: utc

    !=======================================================================
    ! GIGC_CHUNK_RUN begins here 
    !=======================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! The areas need to be updated on every call, because for integration with
    ! WRF nested domains in-PET switching of AREA_M2 and Met fields, even
    ! CMN_SIZE_Mod resolutions are necessary. (hplin, 5/12/18)
    AREA_M2 = State_Met%AREA_M2

    !=======================================================================
    ! Define processes to be covered in this phase
    !
    ! In the standard GEOS-Chem, the following operator sequence is used:
    ! 1. DryDep (kg)
    ! 2. Emissions (kg)
    ! 3. Turbulence (v/v)
    ! 4. Convection (v/v)
    ! 5. Chemistry (kg)
    ! 6. Wetdep (kg)
    !
    ! The GEOS-5 operator sequence is:
    ! 1. Gravity wave drag
    ! 2. Moist (convection)
    ! 3. Chemistry 1 (drydep and emissions)
    ! 4. Surface 1
    ! 5. Turbulence 1
    ! 6. Surface 2
    ! 7. Turbulence 2
    ! 8. Chemistry 2 (chemistry and wet deposition)
    ! 9. Radiation 
    !
    ! Here, we use the following operator sequence:
    ! 
    ! 1.  Convection (v/v) --> Phase 1
    ! 2.  DryDep (kg)      --> Phase 1
    ! 3.  Emissions (kg)   --> Phase 1
    ! 4a. Tendencies (v/v) --> Phase 1
    ! -------------------------------
    ! 4b. Turbulence (v/v) --> Phase 2 
    ! 5.  Chemistry (kg)   --> Phase 2
    ! 6.  WetDep (kg)      --> Phase 2     
    ! 
    ! Any of the listed processes is only executed if the corresponding switch
    ! in the input.geos file is enabled. If the physics component already
    ! covers convection or turbulence, they should not be applied here!
    ! The tendencies are only applied if turbulence is not done within
    ! GEOS-Chem (ckeller, 10/14/14).
    !
    ! For integration with other models and use GIGC_Chunk as a column model
    ! entry point, the "Phase" component has been replaced with a derived type
    ! object of type GIGC_Chunk_Operators, which contain switches that dictate
    ! whether to run Conv, DryDep, Emis, Turb, Chem, WetDep.
    ! However, these are global switches and still have to be controlled by the
    ! respective dynamic/chemistry time-steps & input.geos. (hplin, 5/10/18)
    !=======================================================================

    ! By default, do processes as defined in input.geos. DoTend is defined
    ! below. 
    DoConv   = Input_Opt%LCONV                    ! dynamic time step
    DoDryDep = Input_Opt%LDRYD .AND. IsChemTime   ! chemistry time step
    DoEmis   = Input_Opt%LEMIS .AND. IsChemTime   ! chemistry time step
    DoTurb   = Input_Opt%LTURB                    ! dynamic time step
    DoChem   = Input_Opt%LCHEM .AND. IsChemTime   ! chemistry time step
    DoWetDep = Input_Opt%LWETD                    ! dynamic time step 

    ! Only do selected processes for given operator options.
    DoConv   = DoConv .AND. Operators%Conv
    DoDryDep = DoDryDep .AND. Operators%DryDep
    DoEmis   = DoEmis .AND. Operators%Emis
    DoTurb   = DoTurb .AND. Operators%Turb
    DoChem   = DoChem .AND. Operators%Chem
    DoWetDep = DoWetDep .AND. Operators%WetDep

    ! Check if tendencies need be applied. The drydep and emission calls
    ! only calculates the emission / drydep rates, but do not apply the
    ! tendencies to the tracer array yet. If turbulence is done as part of
    ! GEOS-5, we need to make sure that these tendencies are applied to the
    ! tracer array. If turbulence is explicitly covered by GEOS-Chem,
    ! however, the tendencies become automatically applied within the PBL
    ! mixing routines (DO_MIXING), so we should never apply the tendencies
    ! in this case.
    !
    ! DoTend = ( DoEmis .OR. DoDryDep ) .AND. .NOT. Input_Opt%LTURB
    !
    ! hplin fix 12/6/2018: If turbulence is done, then mixing_mod:do_mixing
    ! will perform tendencies in the PBL mixing routine (call do_tend)
    ! hence it should NOT be done twice. Add .and. .not. DoTurb check here
    ! to match the comments above.
    DoTend = ((DoEmis .OR. DoDryDep) .AND. Operators%Tend) .AND. .NOT. DoTurb

    ! testing only
    IF ( am_I_Root .and. NCALLS < 10 ) THEN 
       write(*,*) 'GEOS-Chem Column Code Operators:'
       write(*,*) 'DoConv   : ', DoConv
       write(*,*) 'DoDryDep : ', DoDryDep
       write(*,*) 'DoEmis   : ', DoEmis
       write(*,*) 'DoTend   : ', DoTend
       write(*,*) 'DoTurb   : ', DoTurb
       write(*,*) 'DoChem   : ', DoChem
       write(*,*) 'DoWetDep : ', DoWetDep
       write(*,*) ' '
       write(*,*) 'Write G-C Diagns : ', Operators%GCDiagn
    ENDIF

    !-------------------------------------------------------------------------
    ! Pre-Run assignments
    !-------------------------------------------------------------------------

    !=======================================================================
    ! Temporary support for in-PET grid switching: verify if we need to
    ! pre-update CMN_SIZE_MOD variables and de/reallocate ahead of time.
    !=======================================================================
    call GIGC_Switch_Dims( am_I_Root, ID, IM, JM, LM, &
                           lonCtr, latCtr, lonEdge, latEdge, &
                           Input_Opt, State_Met, State_Chm, State_Diag, .true., .true., RC )
    write(6, *) "GIGC_Chunk_Run: Running for in-PET grid #", ID
    write(6, *) "GIGC_Chunk_Run: Received new grid IM, JM, LM =", IM, JM, LM

    ! Set lat/lon centers and calculate grid box parameters for the GEOS-Chem Grid.
    ! This is to satisfy HEMCO geographical coordinate requirements
    ! while using GIGC_Chunk_Run as a pure 1-D model that
    ! is technically unaware of any surrounding grid boxes.
    !
    ! As a note, you must ensure that IM, LM, JM are consistent across your calls to
    ! GIGC_Chunk_Init and GIGC_Chunk_Run, as the G-C grid is initialized by
    ! GIGC_Chunk_Init -> GIGC_Init_Simulation and arrays are opened to the extent of
    ! the original IM, JM, LM (_WORLD) parameters.
    ! Failure to do so will result in array overflows and almost certainly SIGSEGV.
    !
    ! Note that we've directly inferred NX, NY (IM, JM in GCHP-speak) from lonCtr, latCtr.
    ! This serves as a consistency check by making use of GC_GRID_MOD::SetGridFromCtr's
    ! checks against SIZE(XMID, 1) & SIZE(XMID, 2) respectively. (hplin, 4/24/18)
    IF (SIZE(lonCtr, 1) .ne. SIZE(latCtr, 1)) THEN
      WRITE(6, *) "GIGC_Chunk_Mod Fatal: Size of lonCtr, latCtr arrays on I do not match"
      RC = GC_FAILURE
      RETURN
    ENDIF

    IF (SIZE(lonCtr, 2) .ne. SIZE(latCtr, 2)) THEN
      WRITE(6, *) "GIGC_Chunk_Mod Fatal: Size of lonCtr, latCtr arrays on J do not match"
      RC = GC_FAILURE
      RETURN
    ENDIF

    CALL SetGridFromCtrEdges( am_I_Root, SIZE(lonCtr, 1), SIZE(lonCtr, 2), lonCtr, latCtr, lonEdge, latEdge, RC )
    IF ( RC /= GC_SUCCESS ) THEN
      WRITE(6, *) "GIGC_Chunk_Mod Fatal: Could not set GC grid parameters from lat/lon ctr/edges"
      RETURN
    ENDIF

    ! Fix negatives, from chem_gridcompmod
    ! This can be seen as an artifact of convection.    
    ! This should ONLY be done once - do not do it between processes.
    ! GEOS-Chem marks certain species as negative so they are not used. Converting these negatives
    ! has detrimental effects. (hplin, 8/21/18) after 6 hours of debugging
    where(State_Chm%Species < 0.0e0)
      State_Chm%Species = 1.0e-36
    endwhere

    ! Calculate troposphere and stratosphere # of levs given State_Met% information.
    ! This will override the values as specified as defaults in ESMF_ mode of
    ! gc_environment_mod::GC_Allocate_All -> CMN_SIZE_Mod::Init_CMN_Size
    ! which sets value_LLTROP = 40, value_LLSTRAT = 59 as fixed defaults.
    !
    ! LLCHEM = # of levs incl. in chemistry grid.
    ! LLSTRAT = max # of levs BELOW stratopause (trop+strat)
    !
    ! The # of levs for these WILL affect LLCHEM in CMN_Size_Mod, used by aerosol
    ! chemistry. If you don't update these quantities things will crash unexpectedly
    ! even if you thought passing LM values are enough to make GIGC work.
    !
    ! In the future we will replace CMN_SIZE_MOD with a N-domain in-PET switching project
    ! developed for WRF integration. (hplin, 5/11/2018)
    !
    ! FIXME: Using levs data based solely on (1, 1), assuming it remains consistent across
    ! the whole column group.
    LLTROP  = LM
    LLSTRAT = LM
    DO L = 1, LM
      IF(State_Met%TROPP(1, 1) .ge. State_Met%PEDGE(1, 1, L)) THEN
        LLTROP  = LLTROP - 1
      ENDIF
    ENDDO

    if(Input_Opt%LUCX) then
      LLCHEM     = LLSTRAT
      LLCHEM_FIX = LLSTRAT
    else
      LLCHEM     = LLTROP
      LLCHEM_FIX = LLTROP_FIX
    endif

    if(am_I_Root.and.NCALLS<10) write(*,*) 'GIGC_CHUNK_RUN: Diagnosed LLTROP, LLSTRAT =', LLTROP, LLSTRAT 

    if(am_I_Root.and.NCALLS<10) write(*,*) "- Running GEOS-Chem Column Code, Call", NCALLS

    ! Eventually initialize/reset wetdep
    IF ( DoConv .OR. DoChem .OR. DoWetDep ) THEN
       CALL SETUP_WETSCAV( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
    ENDIF

    ! Compute Olson Landmap State_Met values (IREG, ILAND...)
    ! from State_Met%LandTypeFrac(I,J,N). Copied from chem_gridcompmod, hplin 12/1/18
    ! Compute State_Met variables IREG, ILAND, IUSE, and FRCLND
    CALL Compute_Olson_Landmap_GCHP( am_I_Root, State_Met, RC )

    ! Pass time values to GEOS-Chem
    ! Elapsed_sec is for diagnosing zenith angle and stuff in HEMCO.
    ! This makes G-C not very "column" as it depends on time, FIXME later (hplin)
    utc =           ( DBLE( hour ) ) + ( DBLE( minute ) / 60e+0_f8 ) + ( DBLE( second ) / 3600e+0_f8 )
    CALL Accept_External_Date_Time( am_I_Root      = am_I_Root,  &
                                    value_NYMD     = nymd,       &  
                                    value_NHMS     = nhms,       &  
                                    value_YEAR     = year,       &  
                                    value_MONTH    = month,      &  
                                    value_DAY      = day,        &  
                                    value_DAYOFYR  = dayOfYr,    &  
                                    value_HOUR     = hour,       &  
                                    value_MINUTE   = minute,     &  
                                    value_SECOND   = second,     &
                                    value_HELAPSED = hElapsed,   &
                                    value_UTC      = utc,        &
                                    RC             = RC         )

    ! Set HEMCO time
    CALL SetHcoTime ( am_I_Root, DoEmis, RC )

    ! Set the pressure at level edges [hPa] from caller
    CALL Accept_External_Pedge    ( am_I_Root      = am_I_Root,  &
                                    State_Met      = State_Met,  &
                                    RC             = RC         )

    ! Initialize surface pressures to match the post-advection pressures
    State_Met%PSC2_WET = State_Met%PS2_WET
    State_Met%PSC2_DRY = State_Met%PS2_DRY
    CALL SET_FLOATING_PRESSURES( am_I_Root, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Define airmass and related quantities, and determine whether to scale
    ! species mixing ratios to account for mass conservation across 
    ! pressure changes. Only scale mixing ratios if transport is OFF and
    ! it is after the first timestep. If transport is ON then species
    ! should be "insensitive" to changes in pressure. If it is the first 
    ! timestep then pressure history is not available for the scaling.
    pUpdate = ((.not.FIRST).and.(.not.Input_Opt%LTRAN))
    CALL AirQnt( am_I_Root, Input_Opt, State_Met, State_Chm, RC, pUpdate )

    ! Save the initial tracer concentrations in the MINIT variable of
    ! GeosCore/strat_chem_mod.F90. This has to be done here, after the
    ! very first call to AIRQNT, because we need State_Chm%AD to have been
    ! populated with non-zero values.  Otherwise the unit conversions will
    ! blow up and cause GCHP to crash. (bmy, 10/19/16)
    IF ( FIRST .and. Input_Opt%LSCHEM ) THEN
       CALL INIT_STRAT_CHEM( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
       Minit_is_set = .true.
    ENDIF

    ! Cap the polar tropopause pressures at 200 hPa, in order to avoid
    ! tropospheric chemistry from happening too high up (cf. J. Logan)
    CALL GIGC_Cap_Tropopause_Prs  ( am_I_Root      = am_I_Root,  &
                                    IM             = IM,         &
                                    JM             = JM,         &
                                    Input_Opt      = Input_Opt,  &
                                    State_Met      = State_Met,  &
                                    RC             = RC         )

    ! Compute PBL quantities
    CALL COMPUTE_PBL_HEIGHT( am_I_Root, State_Met, RC )

    ! Convert species conc units to kg/kg dry prior to Phase 1/2 calls
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, State_Chm, &
                            'kg/kg dry', RC )
    
    ! SDE 05/28/13: Set H2O to STT if relevant
    IF ( IND_('H2O','A') > 0 ) THEN
       CALL SET_H2O_TRAC( am_I_Root, ((.NOT. Input_Opt%LUCX) .OR.    &
                          Input_Opt%LSETH2O ), Input_Opt, State_Met, &
                          State_Chm, RC )
       ! Only force strat once if using UCX
       IF (Input_Opt%LSETH2O) Input_Opt%LSETH2O = .FALSE.
    ENDIF

    !=======================================================================
    ! EMISSIONS. Pass HEMCO Phase 1 which only updates the HEMCO clock
    ! and the HEMCO data list. Should be called every time to make sure 
    ! that the HEMCO clock and the HEMCO data list are up to date.
    !=======================================================================
    CALL EMISSIONS_RUN( am_I_Root, Input_Opt,  State_Met,         &
                        State_Chm, State_Diag, DoEmis, 1, RC       )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!                                PHASE 1                                 !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Debug chemistry state and species information on PET (1, 1, 1) lowest lev.
    ! WRITE(6, *) "%%%% REGISTERED STATE_CHM SPECIES DUMP AFTER CHUNK-SETUP %%%%"
    ! do N = 1, MIN(State_Chm%nSpecies, 8)
    !     WRITE(6, *) "N:", N, "Name:", State_Chm%SpcData(N)%Info%Name
    !     ! WRITE(6, *) "Full Name:", State_Chm%SpcData(N)%Info%FullName
    !     WRITE(6, *) "Value at PET (1,1,1):", State_Chm%Species(1, 1, 1, N)
    !     WRITE(6, *) "Value at PET (7,1,1):", State_Chm%Species(7, 1, 1, N)
    ! enddo

    !=======================================================================
    ! 1. Convection (in v/v)
    ! 
    ! Call GEOS-Chem internal convection routines if convection is enabled
    ! in input.geos. This should only be done if convection is not covered
    ! by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    !=======================================================================
    IF ( DoConv ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do convection now' 

       CALL DO_CONVECTION ( am_I_Root, Input_Opt, State_Met, State_Chm, &
                            State_Diag, RC )

       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Convection done!'
    ENDIF

    !=======================================================================
    ! 2. Dry deposition
    !
    ! Calculate the deposition rates in [s-1].
    !=======================================================================
    IF ( DoDryDep ) THEN
       if(am_I_Root.and.NCALLS<10) THEN
          write(*,*) ' --- Do drydep now'
          write(*,*) '     Use FULL PBL: ', Input_Opt%PBL_DRYDEP
       endif

       ! Do dry deposition
       CALL Do_DryDep( am_I_Root, Input_Opt=Input_Opt, State_Chm=State_Chm, &
                       State_Met=State_Met, State_Diag=State_Diag, RC=RC ) 

       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Drydep done!'
    ENDIF ! Do drydep

    !=======================================================================
    ! 3. Emissions (HEMCO)
    !
    ! HEMCO must be called on first time step to make sure that the HEMCO 
    ! data lists are all properly set up. 
    !=======================================================================
    IF ( DoEmis ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do emissions now'

       ! Do emissions. Pass HEMCO Phase 2 which performs the emissions 
       ! calculations.
       CALL EMISSIONS_RUN ( am_I_Root,  Input_Opt, State_Met, State_Chm, &
                            State_Diag, DoEmis, 2, RC )

       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Emissions done!'
    ENDIF

    !=======================================================================
    ! If physics covers turbulence, simply add the emission and dry 
    ! deposition fluxes calculated above to the tracer array, without caring
    ! about the vertical distribution. The tracer tendencies are only added
    ! to the tracers array after emissions, drydep. So we need to use the
    ! emissions time step here.
    ! Subroutine DO_TEND operates in mass units, e.g. the tracers must be 
    ! in [kg].
    ! SDE 2016-04-05: Input units should now be v/v dry. 
    !=======================================================================
    IF ( DoTend ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*)   &
                           ' --- Add emissions and drydep to tracers'

       ! Get emission time step [s].
       DT = HcoState%TS_EMIS 

       ! Apply tendencies over entire PBL using emission time step.
       CALL DO_TEND ( am_I_Root, Input_Opt, State_Met, State_Chm,  &
                      State_Diag, .FALSE., RC, DT=DT )

       if(am_I_Root.and.NCALLS<10) write(*,*)   &
                                 ' --- Fluxes applied to tracers!' 
    ENDIF ! Tendencies

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!                              PHASE 2                                !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !=======================================================================
    ! 4. Turbulence (v/v)
    !
    ! Call GEOS-Chem internal turbulence routines if turbulence is enabled
    ! in input.geos. This should only be done if turbulence is not covered
    ! by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    ! Subroutine DO_MIXING operates in mixing ratios, e.g. the tracers must 
    ! be in [v/v]. 
    !=======================================================================
    IF ( DoTurb ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do turbulence now'

       ! Do mixing and apply tendencies. This will use the dynamic time step,
       ! which is fine since this call will be executed on every time step. 
       CALL DO_MIXING ( am_I_Root, Input_Opt, State_Met, State_Chm, &
                        State_Diag, RC )

       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Turbulence done!'
    ENDIF

    ! Set tropospheric CH4 concentrations and fill species array with
    ! current values. 
    IF ( .NOT. DoTurb .AND. Input_Opt%ITS_A_FULLCHEM_SIM  &
         .AND. IND_('CH4','A') > 0 ) THEN
       CALL SET_CH4 ( am_I_Root,  Input_Opt, State_Met, State_Chm, &
                      State_Diag, RC )
    ENDIF

    !=======================================================================
    ! 5. Chemistry
    !=======================================================================
    IF ( DoChem ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do chemistry now'

       ! Calculate TOMS O3 overhead. For now, always use it from the
       ! Met field. State_Met%TO3 is imported from PCHEM (ckeller, 10/21/2014).
       CALL COMPUTE_OVERHEAD_O3( am_I_Root, DAY, .TRUE., State_Met%TO3 )

       ! Set H2O to species value if H2O is advected
       IF ( IND_('H2O','A') > 0 ) THEN
          CALL SET_H2O_TRAC( am_I_Root, (.not. Input_Opt%LUCX), Input_Opt, &
                             State_Met, State_Chm, RC )
       ENDIF

       ! Do chemistry
       CALL Do_Chemistry( am_I_Root  = am_I_Root,            & ! Root CPU?
                          Input_Opt  = Input_Opt,            & ! Input Options
                          State_Chm  = State_Chm,            & ! Chemistry State
                          State_Met  = State_Met,            & ! Met State
                          State_Diag = State_Diag,           & ! Diagn State
                          RC         = RC                   )  ! Success?

       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Chemistry done!'
    ENDIF

    !=======================================================================
    ! 6. Wet deposition
    !=======================================================================
    IF ( DoWetDep ) THEN
       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do wetdep now'

       ! Do wet deposition
       CALL DO_WETDEP( am_I_Root, Input_Opt, State_Met, State_Chm,  &
                       State_Diag, RC )

       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Wetdep done!'
    ENDIF

    !=======================================================================
    ! Diagnostics
    !=======================================================================

    IF ( DoChem ) THEN
       ! Recalculate the optical depth at the wavelength(s) specified
       ! in the Radiation Menu. This must be done before the call to any
       ! diagnostic.
       CALL Recompute_OD( am_I_Root, Input_Opt,  State_Met,  &
                          State_Chm, State_Diag, RC         )
    ENDIF
    ! Set certain diagnostics dependent on state at end of step. This
    ! includes species concentration and dry deposition flux.
    CALL Set_Diagnostics_EndofTimestep( am_I_Root,  Input_Opt, &
                                        State_Met,  State_Chm, &
                                        State_Diag, RC )

    ! Archive aerosol mass and PM2.5 diagnostics
    IF ( State_Diag%Archive_AerMass ) THEN
       CALL Set_AerMass_Diagnostic( am_I_Root, Input_Opt,  State_Met, &
                                    State_Chm, State_Diag, RC         )
    ENDIF

    ! History component is temporarily disabled due to errors with
    ! Grid_Registry_Mod, hphlin 11/1/18
    IF ( Operators%GCDiagn .and. .false. ) THEN
      ! %%% REPLICATING GEOS-CHEM CLASSIC (NON-ESMF) FUNCTIONALITY IN GIGC %%% !
      ! (hplin, 9/28/18)
      IF (am_I_Root .and. NCALLS < 10) write(*,*) ' --- Do history now'

      ! Update each HISTORY ITEM from its data source
      CALL History_Update( am_I_Root, RC )

      ! Increment the timestep values by the heartbeat time
      ! This is becasue we need to write out data with the timestamp
      ! at the end of the heartbeat timestep (i.e. at end of run)
      CALL History_SetTime( am_I_Root, RC )

      ! Write HISTORY ITEMS in each diagnostic collection to disk
      ! (or skip writing if it is not the proper output time.)
      CALL History_Write( am_I_Root, State_Chm%Spc_Units, RC )
      IF (am_I_Root .and. NCALLS < 10) write(*,*) ' --- History done!'
    ENDIF

    !=======================================================================
    ! Clean up
    !=======================================================================

    ! testing only
    IF ( NCALLS < 10 ) NCALLS = NCALLS + 1 

    ! First call is done
    FIRST = .FALSE.

    ! Convert units to units of the internal state
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, State_Chm, &
                            'v/v dry', RC )

    ! Return success
    RC = GC_SUCCESS

  END SUBROUTINE GIGC_Chunk_Run
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_chunk_final
!
! !DESCRIPTION: Subroutine GIGC\_CHUNK\_FINAL is the ESMF finalize method for
!  GEOS-Chem.  This routine deallocates pointers and arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Chunk_Final( am_I_Root, Input_Opt,  State_Chm,             &
                               State_Met, State_Diag, RC                    )
!
! !USES:
!
    USE Input_Opt_Mod,    ONLY : OptInput, Cleanup_Input_Opt
    USE State_Chm_Mod,    ONLY : ChmState, Cleanup_State_Chm
    USE State_Met_Mod,    ONLY : MetState, Cleanup_State_Met
    USE State_Diag_Mod,   ONLY : DgnState, Cleanup_State_Diag
    USE HCOI_GC_MAIN_MOD, ONLY : HCOI_GC_FINAL
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt     ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm     ! Chemistry State object
    TYPE(MetState), INTENT(INOUT) :: State_Met     ! Meteorology State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag    ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure
!
! !REVISION HISTORY: 
!  18 Jul 2011 - M. Long     - Initial Version
!  09 Oct 2012 - R. Yantosca - Added comments & cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Chunk_Final
!  01 Nov 2012 - R. Yantosca - Now reference gigc_input_opt_mod.F90
!  19 Sep 2017 - E. Lundgren - Move gigc_finalize content to within this routine
!  13 Feb 2019 - H.P. Lin    - (WRF-GC only) Do not deallocate Input_Opt as its needed
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume succes
    RC = GC_SUCCESS

    ! Finalize HEMCO
    CALL HCOI_GC_FINAL( am_I_Root, .FALSE., RC )
    IF ( am_I_Root ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'HEMCO::Finalize... OK.'
       ELSE
          write(*,'(a)') 'HEMCO::Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Diagnostics State object
    CALL Cleanup_State_Diag( am_I_Root, State_Diag, RC )
    IF ( am_I_Root ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Diag Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Diag Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Chemistry State object
    CALL Cleanup_State_Chm( am_I_Root, State_Chm, RC )
    IF ( am_I_Root ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Chm Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Chm Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Meteorology State object
    CALL Cleanup_State_Met( am_I_Root, State_Met, RC )
    IF ( am_I_Root ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Met Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Met Finalize... FAILURE.'
       ENDIF
    ENDIF
  END SUBROUTINE GIGC_Chunk_Final
!EOC
END MODULE GIGC_Chunk_Mod
