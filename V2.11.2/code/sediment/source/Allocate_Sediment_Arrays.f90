!
! Copyright 2013 RBINS-MUMM
!
! Licensed under the EUPL, Version 1.1 or - as soon they will be approved by
! the European Commission - subsequent versions of the EUPL (the "Licence");
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at:
!
! http://ec.europa.eu/idabc/eupl
!
! Unless required by the applicable law or agreed to in writing, software
! distributed under the Licence is distributed on an "AS IS" basis,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and
! limitations under the Licence.

!************************************************************************
!
! *Allocate_Sediment_Arrays* Allocate/deallocate arrays for the sediment
!                        transport, morphological and dredging/relocation models
!
! Author - Boudewijn Decrop, Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Allocate_Sediment_Arrays.f90  V2.11.1
!
! $Date: 2018-05-23 10:09:17 +0200 (Wed, 23 May 2018) $
!
! $Revision: 1138 $
!
! Description -
!
! Routines - allocate_dar_arrays, allocate_morph_arrays, allocate_sed_arrays, 
!            deallocate_dar_arrays, deallocate_morph_arrays,
!            deallocate_sed_arrays
!
!************************************************************************

!=======================================================================

SUBROUTINE allocate_dar_arrays
!************************************************************************
!
! *allocate_dar_arrays* Allocate arrays for the dredging/relocation models
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Allocate_Dar_Arrays.f90  V2.10.3
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - error_alloc, error_alloc_dar_struc
!
!************************************************************************
!
USE dararrays
USE darpars
USE darswitches
USE dartypes
USE gridpars
USE iopars
USE sedpars
USE syspars
USE dartypes_init
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: nolocs


procname(pglev+1) = 'allocate_dar_arrays'
CALL log_timer_in()

nolocs = MERGE(nocampaigns,nodredgingsites,iopt_dar_dredging_criterion.EQ.3)

!
!1. Time arrays
!--------------
!

IF (iopt_dar_coupling_time.EQ.2) THEN
   ALLOCATE (nt_lag(nolocs),STAT=errstat)
   CALL error_alloc('nt_lag',1,(/nolocs/),kndint)
   nt_lag = 0
ENDIF

IF (iopt_dar_duration.EQ.2) THEN
   ALLOCATE (nt_spread(nolocs),STAT=errstat)
   CALL error_alloc('nt_spread',1,(/nolocs/),kndint)
   nt_spread = 0
ENDIF

IF (iopt_dar_coupling_time.EQ.3) THEN
   
   ALLOCATE (nt_cycle(nolocs),STAT=errstat)
   CALL error_alloc('nt_cycle',1,(/nolocs/),kndint)
   nt_cycle = 0
   
   ALLOCATE (nt_dredge(nolocs),STAT=errstat)
   CALL error_alloc('nt_dredge',1,(/nolocs/),kndint)
   nt_dredge = 0
   
   ALLOCATE (nt_dredge_lastcycle(nolocs),STAT=errstat)
   CALL error_alloc('nt_dredge_lastcycle',1,(/nolocs/),kndint)
   nt_dredge_lastcycle = 0.0
   
   ALLOCATE (nt_fill(nolocs),STAT=errstat)
   CALL error_alloc('nt_fill',1,(/nolocs/),kndint)
   nt_fill = 0
   
   ALLOCATE (nt_idle(nolocs),STAT=errstat)
   CALL error_alloc('nt_idle',1,(/nolocs/),kndint)
   nt_idle = 0
   
   ALLOCATE (nt_passed(nolocs),STAT=errstat)
   CALL error_alloc('nt_passed',1,(/nolocs/),kndint)
   nt_passed = 0
   
   ALLOCATE (nt_relocate(nolocs),STAT=errstat)
   CALL error_alloc('nt_relocate',1,(/nolocs/),kndint)
   nt_relocate = 0

   ALLOCATE (nt_relocate_lastcycle(nolocs),STAT=errstat)
   CALL error_alloc('nt_relocate_lastcycle',1,(/nolocs/),kndint)
   nt_relocate_lastcycle = 0
   
   ALLOCATE (nt_sail_dredge(nolocs),STAT=errstat)
   CALL error_alloc('nt_sail_dredge',1,(/nolocs/),kndint)
   nt_sail_dredge = 0
   
   ALLOCATE (nt_sail_relocate(nolocs),STAT=errstat)
   CALL error_alloc('nt_sail_relocate',1,(/nolocs/),kndint)
   nt_sail_relocate = 0

   ALLOCATE (nt_wait(nolocs),STAT=errstat)
   CALL error_alloc('nt_wait',1,(/nolocs/),kndint)
   nt_wait = 0

   IF (iopt_dar_time_dredge.EQ.2) THEN
      ALLOCATE (nt_dredge_subzone(nolocs,MaxSubzones),STAT=errstat)
      CALL error_alloc('nt_dredge_subzone',2,(/nolocs,MaxSubzones/),kndint)
      nt_dredge_subzone = 0
   ENDIF

   IF (iopt_dar_time_dredge.EQ.3) THEN
      ALLOCATE (nt_track(maxtrackpoints,nolocs),STAT=errstat)
      CALL error_alloc('nt_track',2,(/maxtrackpoints,nolocs/),kndint)
      nt_track = 0
   ENDIF
   
ENDIF

IF (iopt_dar_coupling_space.EQ.3) THEN
   ALLOCATE (nocouplingtimes(nolocs),STAT=errstat)
      CALL error_alloc('nocouplingtimes',1,(/nolocs/),kndint)
      nocouplingtimes = 0
ENDIF

!
!2. Volume arrays
!----------------
!

ALLOCATE (V_dredge(ncloc,nrloc,nolocs),STAT=errstat)
CALL error_alloc('V_dredge',3,(/ncloc,nrloc,nolocs/),kndrtype)
V_dredge = 0.0

ALLOCATE (V_dredge_total(nolocs),STAT=errstat)
CALL error_alloc('V_dredge_total',1,(/nolocs/),kndrtype)
V_dredge_total = 0.0

ALLOCATE (V_lastcycle(nolocs),STAT=errstat)
CALL error_alloc('V_lastcycle',1,(/nolocs/),kndrtype)
V_lastcycle = 0.0

ALLOCATE (V_relocate_extra(norelocationsites),STAT=errstat)
CALL error_alloc('V_relocate_extra',1,(/norelocationsites/),kndrtype)
V_relocate_extra = 0.0

ALLOCATE (V_relocate_total(norelocationsites),STAT=errstat)
CALL error_alloc('V_relocate_total',1,(/norelocationsites/),kndrtype)
V_relocate_total = 0.0

IF (iopt_dar_time_dredge.EQ.3) THEN
   ALLOCATE (V_dredge_max(ncloc,nrloc,nolocs),STAT=errstat)
   CALL error_alloc('V_dredge_max',3,(/ncloc,nrloc,nolocs/),kndrtype)
   V_dredge_max = 0.0
   
   ALLOCATE (V_dredge_rest(ncloc,nrloc,nolocs),STAT=errstat)
   CALL error_alloc('V_dredge_rest',3,(/ncloc,nrloc,nolocs/),kndrtype)
   V_dredge_rest = 0.0

   ALLOCATE (V_ship_rest(nolocs),STAT=errstat)
   CALL error_alloc('V_ship_rest',1,(/nolocs/),kndrtype)
   V_ship_rest = 0.0
ENDIF

!
!3. Counters
!-----------
!

ALLOCATE (couple_counter(nolocs),STAT=errstat)
CALL error_alloc('couple_counter',1,(/nolocs/),kndint)
couple_counter = 1

ALLOCATE (cycle_counter(nolocs),STAT=errstat)
CALL error_alloc('cycle_counter',1,(/nolocs/),kndint)
cycle_counter = 1

ALLOCATE (dredge_counter(nolocs),STAT=errstat)
CALL error_alloc('dredge_counter',1,(/nolocs/),kndint)
dredge_counter = 1

ALLOCATE (priority_counter(nolocs),STAT=errstat)
CALL error_alloc('priority_counter',1,(/nolocs/),kndint)
priority_counter = 1

ALLOCATE (relocate_counter(nolocs),STAT=errstat)
CALL error_alloc('relocate_counter',1,(/nolocs/),kndint)
relocate_counter = 1

ALLOCATE (spatial_counter(nolocs),STAT=errstat)
CALL error_alloc('spatial_counter',1,(/nolocs/),kndint)
spatial_counter = 1

ALLOCATE (track_cell_counter(nolocs),STAT=errstat)
CALL error_alloc('track_cell_counter',1,(/nolocs/),kndint)
track_cell_counter = 1

ALLOCATE (track_nt_counter(nolocs),STAT=errstat)
CALL error_alloc('track_nt_counter',1,(/nolocs/),kndint)
track_nt_counter = 1

ALLOCATE (track_points_counter(nolocs),STAT=errstat)
CALL error_alloc('track_points_counter',1,(/nolocs/),kndint)
track_points_counter = 1

!
!4. Masks
!--------
!

ALLOCATE (mask_dredging(ncloc,nrloc,nodredgingsites),STAT=errstat)
CALL error_alloc('mask_dredging',3,(/ncloc,nrloc,nodredgingsites/),kndlog)
mask_dredging = .FALSE.

ALLOCATE (mask_relocation(ncloc,nrloc,norelocationsites),STAT=errstat)
CALL error_alloc('mask_relocation',3,(/ncloc,nrloc,norelocationsites/),kndlog)
mask_relocation = .FALSE.

ALLOCATE (mask_subzones_dredging(ncloc,nrloc,nodredgingsites,&
     & MaxSubzones),STAT=errstat)
CALL error_alloc('mask_subzones_dredging',4,(/ncloc,nrloc,nodredgingsites,&
                & MaxSubzones/),kndlog)
mask_subzones_dredging = .FALSE.

ALLOCATE (mask_subzones_relocation(ncloc,nrloc,norelocationsites,&
     & MaxSubzones),STAT=errstat)
CALL error_alloc('mask_subzones_relocation',4,(/ncloc,nrloc,norelocationsites,&
                & MaxSubzones/),kndlog)
mask_subzones_relocation = .FALSE.

!
!5. Logical arrays
!-----------------
!

ALLOCATE (campaign(nolocs),STAT=errstat)
CALL error_alloc('campaign',1,(/nolocs/),kndlog)
campaign = .FALSE.

ALLOCATE (dredging(nolocs),STAT=errstat)
CALL error_alloc('dredging',1,(/nolocs/),kndlog)
dredging = .FALSE.

ALLOCATE (eligibility(nolocs),STAT=errstat)
CALL error_alloc('eligibility',1,(/nolocs/),kndlog)
eligibility = .TRUE.

ALLOCATE (relocation(nolocs),STAT=errstat)
CALL error_alloc('relocation',1,(/nolocs/),kndlog)
relocation = .FALSE.

!
!6. Plume modeling
!-----------------
!

ALLOCATE (p_overflow(nolocs,nf),STAT=errstat)
CALL error_alloc('p_overflow',2,(/nolocs,nf/),kndrtype)
p_overflow = 0.0

ALLOCATE (p_passive(nolocs,nf),STAT=errstat)
CALL error_alloc('p_passive',2,(/nolocs,nf/),kndrtype)
p_passive = 0.0

ALLOCATE (p_relocation(nolocs,nf),STAT=errstat)
CALL error_alloc('p_relocation',2,(/nolocs,nf/),kndrtype)
p_relocation = 0.0
  
ALLOCATE (p_suctionhead(nolocs,nf),STAT=errstat)
CALL error_alloc('p_suctionhead',2,(/nolocs,nf/),kndrtype)
p_suctionhead = 0.0
   
IF (iopt_dar_effect.EQ.2.OR.iopt_dar_effect.EQ.3) THEN

   ALLOCATE (D_dredge_max(nolocs),STAT=errstat)
   CALL error_alloc('D_dredge_max',1,(/nolocs/),kndrtype)
   D_dredge_max = 0.0
   
   ALLOCATE (D_dredge_min(nolocs),STAT=errstat)
   CALL error_alloc('D_dredge_min',1,(/nolocs/),kndrtype)
   D_dredge_min = 0.0

   ALLOCATE (D_relocate_max(nolocs),STAT=errstat)
   CALL error_alloc('D_relocate_max',1,(/nolocs/),kndrtype)
   D_relocate_max = 0.0
   
   ALLOCATE (D_relocate_min(nolocs),STAT=errstat)
   CALL error_alloc('D_relocate_min',1,(/nolocs/),kndrtype)
   D_relocate_min = 0.0

   ALLOCATE (H_dredge_max(nolocs),STAT=errstat)
   CALL error_alloc('H_dredge_max',1,(/nolocs/),kndrtype)
   H_dredge_max = 0.0
   
   ALLOCATE (H_dredge_min(nolocs),STAT=errstat)
   CALL error_alloc('H_dredge_min',1,(/nolocs/),kndrtype)
   H_dredge_min = 0.0

   ALLOCATE (H_relocate_max(nolocs),STAT=errstat)
   CALL error_alloc('H_relocate_max',1,(/nolocs/),kndrtype)
   H_relocate_max = 0.0
   
   ALLOCATE (H_relocate_min(nolocs),STAT=errstat)
   CALL error_alloc('H_relocate_min',1,(/nolocs/),kndrtype)
   H_relocate_min = 0.0

ENDIF

!
!7. Other
!--------
!

ALLOCATE (dredging_depth(ncloc,nrloc,nolocs),STAT=errstat)
CALL error_alloc('dredging_depth',3,(/ncloc,nrloc,nolocs/),kndrtype)
dredging_depth = 0.0

ALLOCATE (fraction_ship(noships,nf),STAT=errstat)
CALL error_alloc('fraction_ship',2,(/noships,nf/),kndrtype)
fraction_ship = 0.0

ALLOCATE (volume_ship(noships),STAT=errstat)
CALL error_alloc('volume_ship',1,(/noships/),kndrtype)
volume_ship = 0.0

IF (iopt_dar_coupling_time.EQ.3) THEN
   ALLOCATE (nr_cycles(nolocs),STAT=errstat)
   CALL error_alloc('nr_cycles',1,(/nolocs/),kndint)
   nr_cycles = 0
ENDIF

IF (iopt_dar_dredging_criterion.EQ.3) THEN
   ALLOCATE (campaign_start(nocampaigns),STAT=errstat)
   CALL error_alloc('campaign_start',1,(/nolocs/),kndint)
   campaign_start = 0
ENDIF

IF (iopt_dar_time_dredge.EQ.3) THEN

   ALLOCATE (notrackpoints(nodredgingsites),STAT=errstat)
   CALL error_alloc('notrackpoints',1,(/nodredgingsites/),kndint)
   notrackpoints = 0

   ALLOCATE (itrack(maxtrackpoints,nolocs),STAT=errstat)
   CALL error_alloc('itrack',2,(/maxtrackpoints,nolocs/),kndint)
   itrack = 0

   ALLOCATE (jtrack(maxtrackpoints,nolocs),STAT=errstat)
   CALL error_alloc('jtrack',2,(/maxtrackpoints,nodredgingsites/),kndint)
   jtrack = 0

   ALLOCATE (xtrack(maxtrackpoints,nolocs),STAT=errstat)
   CALL error_alloc('xtrack',2,(/maxtrackpoints,nolocs/),kndrtype)
   xtrack= 0.0

   ALLOCATE (ytrack(maxtrackpoints,nolocs),STAT=errstat)
   CALL error_alloc('ytrack',2,(/maxtrackpoints,nolocs/),kndrtype)
   ytrack = 0.0
    
ENDIF

!
!8. Derived data types
!---------------------
!

ALLOCATE (Campaigns(nocampaigns),STAT=errstat)
CALL error_alloc_dar_struc('Campaigns',1,(/nocampaigns/),'CampaignParams')
CALL campaigns_init(Campaigns)

ALLOCATE (DredgingSites(nodredgingsites),STAT=errstat)
CALL error_alloc_dar_struc('DredgingSites',1,(/nodredgingsites/),&
                         & 'DredgingParams')
CALL dredgingsites_init(DredgingSites)

ALLOCATE (OperationTimes(nolocs),STAT=errstat)
CALL error_alloc_dar_struc('OperationTimes',1,(/nolocs/),'OperationTimesParams')
CALL operationtimes_init(OperationTimes)

ALLOCATE (RelocationSites(norelocationsites),STAT=errstat)
CALL error_alloc_dar_struc('RelocationSites',1,(/norelocationsites/),&
                         & 'RelocationParams')
CALL relocationsites_init(RelocationSites)

ALLOCATE (Ships(noships),STAT=errstat)
CALL error_alloc_dar_struc('Ships',1,(/noships/),'ShipParams')
CALL ships_init(Ships)

ALLOCATE (SubzonesDredging(nodredgingsites),STAT=errstat)
CALL error_alloc_dar_struc('SubzonesDredging',1,(/nodredgingsites/),&
                         & 'SubzonesDredgingParams')
CALL subzonesdredging_init(SubzonesDredging)

ALLOCATE (SubzonesRelocation(norelocationsites),STAT=errstat)
CALL error_alloc_dar_struc('SubzonesRelocation',1,(/norelocationsites/), &
                         & 'SubzonesRelocationParams')
CALL subzonesrelocation_init(SubzonesRelocation)

CALL log_timer_out()


RETURN

END SUBROUTINE allocate_dar_arrays

!=======================================================================

SUBROUTINE allocate_morph_arrays 
!************************************************************************
!
! *allocate_morphology_arrays* Allocate arrays for the morphological model
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Allocate_Sediment_Arrays.f90  V2.10
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - error_alloc, error_alloc_struc
!
!************************************************************************
!
USE gridpars
USE iopars
USE morpharrays
USE morphpars
USE morphswitches
USE sedpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'allocate_morph_arrays'
CALL log_timer_in()

!
!1. General arrays
!-----------------
!

ALLOCATE (active_layer_thickness(ncloc,nrloc),STAT=errstat)
CALL error_alloc('active_layer_thickness',2,(/ncloc,nrloc/),kndrtype)
active_layer_thickness = 0.0

ALLOCATE (bed_fraction_old(ncloc,nrloc,nb,nf),STAT=errstat)
CALL error_alloc('bed_fraction_old',4,(/ncloc,nrloc,nb,nf/),kndrtype)
bed_fraction_old = 0.0

ALLOCATE (bed_layer_thickness(ncloc,nrloc,nb),STAT=errstat)
CALL error_alloc('bed_layer_thickness',3,(/ncloc,nrloc,nb/),kndrtype)
bed_layer_thickness = 0.0

ALLOCATE (bed_layer_thickness_old(ncloc,nrloc,nb),STAT=errstat)
CALL error_alloc('bed_layer_thickness_old',3,(/ncloc,nrloc,nb/),kndrtype)
bed_layer_thickness_old = 0.0

ALLOCATE (bed_porosity(ncloc,nrloc,nb),STAT=errstat)
CALL error_alloc('bed_porosity',3,(/ncloc,nrloc,nb/),kndrtype)
bed_porosity = 0.0

ALLOCATE (bed_update(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bed_update',2,(/ncloc,nrloc/),kndrtype)
bed_update = 0.0

ALLOCATE (bed_update_dep(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('bed_update_dep',3,(/ncloc,nrloc,nf/),kndrtype)
bed_update_dep = 0.0

ALLOCATE (bed_update_ero(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('bed_update_ero',3,(/ncloc,nrloc,nf/),kndrtype)
bed_update_ero = 0.0

ALLOCATE (bed_update_int(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bed_update_int',2,(/ncloc,nrloc/),kndrtype)
bed_update_int = 0.0

ALLOCATE (dar_morph(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('dar_morph',3,(/ncloc,nrloc,nf/),kndrtype)
dar_morph = 0.0

!
!2. Fixed layer
!--------------
!

IF (iopt_morph_fixed_layer.EQ.1) THEN

   ALLOCATE (sed_avail(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('sed_avail',3,(/ncloc,nrloc,nf/),kndrtype)
   sed_avail = 0.0

   ALLOCATE (sed_avail_tot(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('sed_avail_tot',2,(/ncloc,nrloc/),kndrtype)
   sed_avail_tot = 0.0

ENDIF

!
!3. Time integration
!-------------------
!

IF (iopt_morph_time_int.GT.1) THEN

   ALLOCATE (bed_update_dep_old1(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('bed_update_dep_old1',3,(/ncloc,nrloc,nf/),kndrtype)
   bed_update_dep_old1 = 0.0

   ALLOCATE (bed_update_ero_old1(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('bed_update_ero_old1',3,(/ncloc,nrloc,nf/),kndrtype)
   bed_update_ero_old1 = 0.0

   IF (iopt_morph_time_int.EQ.3) THEN

      ALLOCATE (bed_update_dep_old2(ncloc,nrloc,nf),STAT=errstat)
      CALL error_alloc('bed_update_dep_old2',3,(/ncloc,nrloc,nf/),kndrtype)
      bed_update_dep_old2 = 0.0

      ALLOCATE (bed_update_ero_old2(ncloc,nrloc,nf),STAT=errstat)
      CALL error_alloc('bed_update_ero_old2',3,(/ncloc,nrloc,nf/),kndrtype)
      bed_update_ero_old2 = 0.0
      
   ENDIF

ENDIF

!
!4. Arrays for tidal averaging
!-----------------------------
!

IF (iopt_tidal_accel.EQ.1) THEN

   ALLOCATE (depth_guess(ncloc,nrloc,number_tidal_steps),STAT=errstat)
   CALL error_alloc('depth_guess',3,(/ncloc,nrloc,number_tidal_steps/),&
                  & kndrtype)
   depth_guess = 0.0

   ALLOCATE (tidalstep(number_tidal_steps),STAT=errstat)
   CALL error_alloc('tidalstep',1,(/number_tidal_steps/),kndint)
   tidalstep = 0

   ALLOCATE (umvel_guess(ncloc+1,nrloc,number_tidal_steps),STAT=errstat)
   CALL error_alloc('umvel_guess',3,(/ncloc+1,nrloc,number_tidal_steps/),&
                  & kndrtype)
   umvel_guess = 0.0

   ALLOCATE (vmvel_guess(ncloc,nrloc+1,number_tidal_steps),STAT=errstat)
   CALL error_alloc('vmvel_guess',3,(/ncloc,nrloc+1,number_tidal_steps/),&
                  & kndrtype)
   vmvel_guess = 0.0

   ALLOCATE (zeta_guess(0:ncloc+1,0:nrloc+1,number_tidal_steps),STAT=errstat)
   CALL error_alloc('zeta_guess',3,(/ncloc+2,nrloc+2,number_tidal_steps/),&
                  & kndrtype)
   zeta_guess = 0.0

ENDIF

CALL log_timer_out()


END SUBROUTINE allocate_morph_arrays

!========================================================================

SUBROUTINE allocate_sed_arrays
!************************************************************************
!
! *allocate_sed_arrays* Allocate arrays for the sediment transport model
!
! Author - Boudewijn Decrop, Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Allocate_Sediment_Arrays.f90  V2.11.1
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - error_alloc
!
!************************************************************************
!
USE gridpars
USE iopars
USE morphpars
USE morphswitches
USE sedarrays
USE sedpars
USE sedswitches
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'allocate_sed_arrays'
CALL log_timer_in()

!
!1. Particle attributes
!----------------------
!

ALLOCATE (bstres_cr_cst(nfp),STAT=errstat)
CALL error_alloc('bstres_cr_cst',1,(/nfp/),kndrtype)
bstres_cr_cst = 0.0

ALLOCATE (dp(nfp),STAT=errstat)
CALL error_alloc('dp',1,(/nfp/),kndrtype)
dp = 0.0

ALLOCATE (massp(nfp),STAT=errstat)
CALL error_alloc('massp',1,(/nfp/),kndrtype)
massp = 0.0

ALLOCATE (rhos(nfp),STAT=errstat)
CALL error_alloc('rhos',1,(/nfp/),kndrtype)
rhos = 0.0

ALLOCATE (volp(nfp),STAT=errstat)
CALL error_alloc('volp',1,(/nfp/),kndrtype)
volp = 0.0

ALLOCATE (ws_cst(nfp),STAT=errstat)
CALL error_alloc('ws_cst',1,(/nfp/),kndrtype)
ws_cst = 0.0

!
!2. Sediment concentrations
!--------------------------
!

ALLOCATE (ceq(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('ceq',3,(/ncloc,nrloc,nf/),kndrtype)
ceq = 0.0

ALLOCATE (cref(0:ncloc,0:nrloc,nf),STAT=errstat)
CALL error_alloc('cref',3,(/ncloc+1,nrloc+1,nf/),kndrtype)
cref = 0.0

ALLOCATE (ctot(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
CALL error_alloc('ctot',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
ctot = 0.0

ALLOCATE (cvol(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz,nf),STAT=errstat)
CALL error_alloc('cvol',4,(/ncloc+2*nhalo,nrloc+2*nhalo,nz,nf/),kndrtype)
cvol = 0.0

IF (iopt_sed.EQ.1) THEN
   ALLOCATE (dar_sediment(ncloc,nrloc,nz,nf),STAT=errstat)
   CALL error_alloc('dar_sediment',4,(/ncloc,nrloc,nz,nf/),kndrtype)
   dar_sediment = 0.0
ENDIF

!
!3. Flocculation model
!---------------------
!

IF (iopt_sed.EQ.2) THEN
   
   ALLOCATE (cnump(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz,nf),STAT=errstat)
   CALL error_alloc('cnump',4,(/ncloc+2*nhalo,nrloc+2*nhalo,nz,nf/),kndrtype)
   cnump = 0.0

   ALLOCATE (floc_dens(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
   CALL error_alloc('floc_dens',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
   floc_dens = 0.0

   ALLOCATE (floc_dia(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('floc_dia',3,(/ncloc,nrloc,nz/),kndrtype)
   floc_dia = 0.0

   ALLOCATE (floc_nc(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
   CALL error_alloc('floc_nc',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
   floc_nc = 0.0
   
   ALLOCATE (floc_vol(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('floc_vol',3,(/ncloc,nrloc,nz/),kndrtype)
   floc_vol = 0.0

ENDIF
   
!
!4. Sediment loads
!-----------------
!
!---bed load
IF (iopt_sed_mode.EQ.1.OR.iopt_sed_mode.EQ.3.OR.iopt_sed_toteq.GT.4) THEN

   ALLOCATE (qbedatu(0:ncloc+1,nrloc,nf),STAT=errstat)
   CALL error_alloc('qbedatu',3,(/ncloc+2,nrloc,nf/),kndrtype)
   qbedatu = 0.0

   ALLOCATE (qbedatv(ncloc,0:nrloc+1,nf),STAT=errstat)
   CALL error_alloc('qbedatv',3,(/ncloc,nrloc+2,nf/),kndrtype)
   qbedatv = 0.0

   IF (iopt_morph.EQ.1) THEN
      ALLOCATE (qbedatu_int(0:ncloc+1,nrloc,nf),STAT=errstat)
      CALL error_alloc('qbedatu_int',3,(/ncloc+2,nrloc,nf/),kndrtype)
      qbedatu_int = 0.0

      ALLOCATE (qbedatv_int(ncloc,0:nrloc+1,nf),STAT=errstat)
      CALL error_alloc('qbedatv_int',3,(/ncloc,nrloc+2,nf/),kndrtype)
      qbedatv_int = 0.0
   ENDIF

ENDIF

!---total load
IF (iopt_sed_mode.NE.1) THEN

   ALLOCATE (qtotatu(0:ncloc+1,nrloc,nf),STAT=errstat)
   CALL error_alloc('qtotatu',3,(/ncloc+2,nrloc,nf/),kndrtype)
   qtotatu = 0.0

   ALLOCATE (qtotatv(ncloc,0:nrloc+1,nf),STAT=errstat)
   CALL error_alloc('qtotatv',3,(/ncloc,nrloc+2,nf/),kndrtype)
   qtotatv = 0.0

   IF (iopt_morph.EQ.1) THEN

      ALLOCATE (qtotatu_int(0:ncloc+1,nrloc,nf),STAT=errstat)
      CALL error_alloc('qtotatu_int',3,(/ncloc+2,nrloc,nf/),kndrtype)
      qtotatu_int = 0.0
      ALLOCATE (qtotatv_int(ncloc,0:nrloc+1,nf),STAT=errstat)
      CALL error_alloc('qtotatv_int',3,(/ncloc,nrloc+2,nf/),kndrtype)
      qtotatv_int = 0.0

   ENDIF

ENDIF

!---suspended load
IF (iopt_sed_nodim.EQ.2.AND.iopt_sed_bbc_eq.GT.1) THEN
   ALLOCATE (qsusatc(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('qsusatc',3,(/ncloc,nrloc,nf/),kndrtype)
   qsusatc = 0.0
ENDIF

IF (iopt_sed_mode.EQ.2.OR.iopt_sed_mode.EQ.3) THEN
   ALLOCATE (qsusatu(0:ncloc+1,nrloc,nz,nf),STAT=errstat)
   CALL error_alloc('qsusatu',4,(/ncloc+2,nrloc,nz,nf/),kndrtype)
   qsusatu = 0.0

   ALLOCATE (qsusatv(ncloc,0:nrloc+1,nz,nf),STAT=errstat)
   CALL error_alloc('qsusatv',4,(/ncloc,nrloc+2,nz,nf/),kndrtype)
   qsusatv = 0.0

   IF (iopt_morph.EQ.1) THEN
      ALLOCATE (qsusatu_int(0:ncloc+1,nrloc,nf),STAT=errstat)
      CALL error_alloc('qsusatu_int',3,(/ncloc+2,nrloc,nf/),kndrtype)
      qsusatu_int = 0.0
      ALLOCATE (qsusatv_int(ncloc,0:nrloc+1,nf),STAT=errstat)
      CALL error_alloc('qsusatv_int',3,(/ncloc,nrloc+2,nf/),kndrtype)
      qsusatv_int = 0.0

   ENDIF

ENDIF

!
!5. Bottom fluxes
!----------------
!

IF (iopt_sed_mode.GE.2) THEN
   
  ALLOCATE (bottom_sed_dep(ncloc,nrloc,nf),STAT=errstat)
  CALL error_alloc('bottom_sed_dep',3,(/ncloc,nrloc,nf/),kndrtype)
  bottom_sed_dep = 0.0

  ALLOCATE (bottom_sed_ero(ncloc,nrloc,nf),STAT=errstat)
  CALL error_alloc('bottom_sed_ero',3,(/ncloc,nrloc,nf/),kndrtype)
  bottom_sed_ero = 0.0
  
  ALLOCATE (bottom_sed_flux(ncloc,nrloc,nf),STAT=errstat)
  CALL error_alloc('bottom_sed_flux',3,(/ncloc,nrloc,nf/),kndrtype)
  bottom_sed_flux = 0.0

ENDIF

ALLOCATE (height_c(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('height_c',3,(/ncloc,nrloc,nf/),kndrtype)
height_c = 0.0

IF (iopt_sed_nodim.EQ.2) THEN

   ALLOCATE (t_equil(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('t_equil',3,(/ncloc,nrloc,nf/),kndrtype)
   t_equil = 0.0

ENDIF

ALLOCATE (bed_fraction(0:ncloc,0:nrloc,nb,nf),STAT=errstat)
CALL error_alloc('bed_fraction',4,(/ncloc+1,nrloc+1,nb,nf/),kndrtype)
bed_fraction = 1.0

!
!6. Bottom stress related arrays
!-------------------------------
!

ALLOCATE (bdragcoefatc_sed(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bdragcoefatc_sed',2,(/ncloc,nrloc/),kndrtype)
bdragcoefatc_sed = 0.0

ALLOCATE (bdragcoefatu_sed(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bdragcoefatu_sed',2,(/ncloc,nrloc/),kndrtype)
bdragcoefatu_sed = 0.0

ALLOCATE (bdragcoefatv_sed(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bdragcoefatv_sed',2,(/ncloc,nrloc/),kndrtype)
bdragcoefatv_sed = 0.0

ALLOCATE (bstres_cr(0:ncloc,0:nrloc,nf),STAT=errstat)
CALL error_alloc('bstres_cr',3,(/ncloc+1,nrloc+1,nf/),kndrtype)
bstres_cr = 0.0

ALLOCATE (bstresatc_sed(0:ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('bstresatc_sed',2,(/ncloc+1,nrloc+1/),kndrtype)
bstresatc_sed = 0.0

ALLOCATE (bstresatu_sed(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bstresatu_sed',2,(/ncloc,nrloc/),kndrtype)
bstresatu_sed = 0.0

ALLOCATE (bstresatv_sed(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bstresatv_sed',2,(/ncloc,nrloc/),kndrtype)
bstresatv_sed = 0.0

ALLOCATE (d50_bed(0:ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('d50_bed',2,(/ncloc+1,nrloc+1/),kndrtype)
d50_bed = 0.0

ALLOCATE (ubstresatu_sed(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('ubstresatu_sed',2,(/ncloc+1,nrloc/),kndrtype)
ubstresatu_sed = 0.0

ALLOCATE (vbstresatv_sed(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('vbstresatv_sed',2,(/ncloc,nrloc+1/),kndrtype)
vbstresatv_sed = 0.0

ALLOCATE (zroughatc_sed(0:ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('zroughatc_sed',2,(/ncloc+1,nrloc+1/),kndrtype)
zroughatc_sed = 0.0

ALLOCATE (zroughatu_sed(ncloc,nrloc),STAT=errstat)
CALL error_alloc('zroughatc_sed',2,(/ncloc,nrloc/),kndrtype)
zroughatu_sed = 0.0

ALLOCATE (zroughatv_sed(ncloc,nrloc),STAT=errstat)
CALL error_alloc('zroughatv_sed',2,(/ncloc,nrloc/),kndrtype)
zroughatv_sed = 0.0

IF (iopt_waves.GT.0) THEN
   
   ALLOCATE (bstresatc_mean_sed(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bstresatc_mean_sed',2,(/ncloc,nrloc/),kndrtype)
   bstresatc_mean_sed = 0.0
   
   ALLOCATE (bstresatc_wav_sed(0:ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('bstresatc_wav_sed',2,(/ncloc+1,nrloc+1/),kndrtype)
   bstresatc_wav_sed = 0.0

   ALLOCATE (fwave_sed(0:ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('fwave_sed',2,(/ncloc+1,nrloc+1/),kndrtype)
   fwave_sed = 0.0

   ALLOCATE (wavethickatc_sed(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('wavethickatc_sed',2,(/ncloc,nrloc/),kndrtype)
   wavethickatc_sed = 0.0

   ALLOCATE (zaroughatc_sed(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('zaroughatc_sed',2,(/ncloc,nrloc/),kndrtype)
   zaroughatc_sed = 0.0

ENDIF

!
!7. Diffusion coefficients
!-------------------------
!

ALLOCATE (beta_sed(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('beta_sed',3,(/ncloc,nrloc,nf/),kndrtype)
beta_sed = 0.0

ALLOCATE (vdiffcoef_sed(ncloc,nrloc,nz+1,nf),STAT=errstat)
CALL error_alloc('vdiffcoef_sed',4,(/ncloc,nrloc,nz+1,nf/),kndrtype)
vdiffcoef_sed = 0.0

!
!8. Bed slope arrays
!-------------------
!

IF (iopt_sed_slope.GT.0.OR.iopt_morph_avalanching.GT.0) THEN

   ALLOCATE (bed_slope_x(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bed_slope_x',2,(/ncloc,nrloc/),kndrtype)
   bed_slope_x = 0.0

   ALLOCATE (bed_slope_x_atu(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bed_slope_x_atu',2,(/ncloc,nrloc/),kndrtype)
   bed_slope_x_atu = 0.0

   ALLOCATE (bed_slope_x_atv(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bed_slope_x_atv',2,(/ncloc,nrloc/),kndrtype)
   bed_slope_x_atv = 0.0

   ALLOCATE (bed_slope_y(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bed_slope_y',2,(/ncloc,nrloc/),kndrtype)
   bed_slope_y = 0.0

   ALLOCATE (bed_slope_y_atu(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bed_slope_y_atu',2,(/ncloc,nrloc/),kndrtype)
   bed_slope_y_atu = 0.0

   ALLOCATE (bed_slope_y_atv(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bed_slope_y_atv',2,(/ncloc,nrloc/),kndrtype)
   bed_slope_y_atv = 0.0

ENDIF

!
!9. Density
!----------
!

ALLOCATE (densw(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
CALL error_alloc('densw',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
densw = 0.0

ALLOCATE (beta_state_sed(0:ncloc+1,0:nrloc+1,nz,nf),STAT=errstat)
CALL error_alloc('beta_state_sed',4,(/ncloc+2,nrloc+2,nz,nf/),kndrtype)
beta_state_sed = 0.0

!
!10. Fall velocity
!-----------------
!

ALLOCATE (wfall(0:ncloc,0:nrloc,nz+1,nf),STAT=errstat)
CALL error_alloc('wfall',4,(/ncloc+1,nrloc+1,nz+1,nf/),kndrtype)
wfall = 0.0

!
!11. Open boundary conditions
!----------------------------
!

IF (iopt_obc_sed.EQ.1) THEN

   IF (iopt_sed.EQ.1) THEN
      
      ALLOCATE (obcsedatu(nobu,nz,0:2,nf),STAT=errstat)
      CALL error_alloc('obcsedatu',4,(/nobu,nz,3,nf/),kndrtype)
      IF (nobu.GT.0) obcsedatu = 0.0

      ALLOCATE (obcsedatv(nobv,nz,0:2,nf),STAT=errstat)
      CALL error_alloc('obcsedatv',4,(/nobv,nz,3,nf/),kndrtype)
      IF (nobv.GT.0) obcsedatv = 0.0

   ELSEIF (iopt_sed.EQ.2) THEN
      
      ALLOCATE (obccnumpatu(nobu,nz,0:2,nf),STAT=errstat)
      CALL error_alloc('obccnumpatu',4,(/nobu,nz,3,nf/),kndrtype)
      IF (nobu.GT.0) obccnumpatu = 0.0

      ALLOCATE (obccnumpatv(nobv,nz,0:2,nf),STAT=errstat)
      CALL error_alloc('obccnumpatv',4,(/nobv,nz,3,nf/),kndrtype)
      IF (nobv.GT.0) obccnumpatv = 0.0
      
   ENDIF
      
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE allocate_sed_arrays

!========================================================================

SUBROUTINE deallocate_dar_arrays
!************************************************************************
!
! *deallocate_dar_arrays* Deallocate arrays for dredging/relocation model
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Allocate_Dar_Arrays.f90  V2.10
!
! Description -
!
! Reference -
!
! Calling program - simulation_end
!
!************************************************************************
!
USE iopars
USE dararrays
USE darswitches
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'deallocate_dar_arrays'
CALL log_timer_in()

!---time arrays
IF (iopt_dar_coupling_time.EQ.2) DEALLOCATE (nt_lag)
IF (iopt_dar_duration.EQ.2) DEALLOCATE (nt_spread)
IF (iopt_dar_coupling_time.EQ.3) THEN
   DEALLOCATE (nt_cycle,nt_dredge,nt_dredge_lastcycle,nt_fill,nt_idle,&
             & nt_passed,nt_relocate,nt_relocate_lastcycle,nt_sail_dredge,&
             & nt_sail_relocate,nt_wait)
   IF (iopt_dar_time_dredge.EQ.2) DEALLOCATE (nt_dredge_subzone)
   IF (iopt_dar_time_dredge.EQ.3) DEALLOCATE (nt_track)
ENDIF
IF (iopt_dar_coupling_space.EQ.3) DEALLOCATE (nocouplingtimes)

!--- volume arrays
DEALLOCATE (V_dredge,V_dredge_total,V_lastcycle,V_relocate_extra,&
          & V_relocate_total)
IF (iopt_dar_time_dredge.EQ.3) THEN
   DEALLOCATE (V_dredge_max,V_dredge_rest,V_ship_rest)
ENDIF

!---counters
DEALLOCATE (couple_counter,cycle_counter,dredge_counter,priority_counter,&
          & relocate_counter,spatial_counter,track_cell_counter,&
          & track_nt_counter,track_points_counter)

!--- masks
DEALLOCATE (mask_dredging,mask_relocation,mask_subzones_dredging,&
          & mask_subzones_relocation)

!---logical arrays
DEALLOCATE (campaign,dredging,eligibility,relocation)

!-- plume modelling
DEALLOCATE (p_overflow,p_passive,p_relocation,p_suctionhead)
IF (iopt_dar_effect.EQ.2.OR.iopt_dar_effect.EQ.3) THEN
   DEALLOCATE (D_dredge_max,D_dredge_min,D_relocate_max,D_relocate_min,&
             & H_dredge_max,H_dredge_min,H_relocate_max,H_relocate_min)
ENDIF

!--- other
DEALLOCATE (dredging_depth,fraction_ship,volume_ship)
IF (iopt_dar_coupling_time.EQ.3) DEALLOCATE (nr_cycles)
IF (iopt_dar_dredging_criterion.EQ.3) DEALLOCATE (campaign_start)
IF(iopt_dar_time_dredge.EQ.3) THEN
   DEALLOCATE (notrackpoints,itrack,jtrack,xtrack,ytrack)
ENDIF

!---derived data types
DEALLOCATE (Campaigns,DredgingSites,OperationTimes,RelocationSites,Ships,&
          & SubzonesDredging,SubzonesRelocation)

CALL log_timer_out()


RETURN

END SUBROUTINE deallocate_dar_arrays

!========================================================================

SUBROUTINE deallocate_morph_arrays
!************************************************************************
!
! *deallocate_morph_arrays* Deallocate arrays for the morphological model
!
! Author - Alexander Breugem, IMDC
!
! Version - @(COHERENS)Allocate_Sediment_Arrays.f90  V2.10
!
! Description -
!
! Reference -
!
! Calling program - simulation_end
!
!************************************************************************
!
USE iopars
USE morpharrays
USE morphpars
USE morphswitches
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'deallocate_morph_arrays'
CALL log_timer_in

!---general arrays
DEALLOCATE (active_layer_thickness,bed_fraction_old,bed_layer_thickness,&
          & bed_layer_thickness_old,bed_porosity,bed_update,bed_update_dep,&
          & bed_update_ero,bed_update_int,dar_morph)

!---sorting or fixed layer
IF (iopt_morph_fixed_layer.EQ.1) DEALLOCATE (sed_avail,sed_avail_tot)

!---time integration
IF (iopt_morph_time_int.GT.1) THEN
   DEALLOCATE (bed_update_dep_old1,bed_update_ero_old1)
   IF (iopt_morph_time_int.EQ.3) THEN
      DEALLOCATE (bed_update_dep_old2,bed_update_ero_old2)
   ENDIF
ENDIF

!---arrays for tidal averaging
IF (iopt_tidal_accel.EQ.1) THEN
   DEALLOCATE(depth_guess,tidalstep,umvel_guess,vmvel_guess,zeta_guess)
ENDIF

CALL log_timer_out()


END SUBROUTINE deallocate_morph_arrays

!=======================================================================

SUBROUTINE deallocate_sed_arrays
!************************************************************************
!
! *deallocate_sed_arrays* Deallocate arrays for the sediment transport model
!
! Author - Boudewijn Decrop, Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Allocate_Sediment_Arrays.f90  V2.11.1
!
! Description -
!
! Reference -
!
! External calls - deallocate_morph_arrays
!
! Calling program - simulation_end
!
!************************************************************************
!
USE iopars
USE morphswitches
USE sedarrays
USE sedswitches
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'deallocate_sed_arrays'
CALL log_timer_in()

!---particle attributes
DEALLOCATE (bstres_cr_cst,dp,massp,rhos,volp,ws_cst)

!---sediment concentrations
DEALLOCATE (ceq,cref,ctot,cvol)
IF (iopt_sed.EQ.1) DEALLOCATE (dar_sediment)

!---flocculation model
IF (iopt_sed.EQ.2) DEALLOCATE (cnump,floc_dens,floc_dia,floc_nc,floc_vol)

!---bed load
IF (iopt_sed_mode.EQ.1.OR.iopt_sed_mode.EQ.3.OR.iopt_sed_toteq.GT.4) THEN
   DEALLOCATE (qbedatu,qbedatv)
   IF (iopt_morph.EQ.1) DEALLOCATE (qbedatu_int,qbedatv_int)
ENDIF

!---total load
IF (iopt_sed_mode.NE.1) THEN
   DEALLOCATE (qtotatu,qtotatv)
   IF (iopt_morph.EQ.1) DEALLOCATE (qtotatu_int,qtotatv_int)
ENDIF

!---suspended load
IF (iopt_sed_nodim.EQ.2.AND.iopt_sed_bbc_eq.GT.1) DEALLOCATE (qsusatc) 
IF (iopt_sed_mode.EQ.2.OR.iopt_sed_mode.EQ.3) THEN
   DEALLOCATE (qsusatu,qsusatv)
   IF (iopt_morph.EQ.1) DEALLOCATE (qsusatu_int,qsusatv_int)
ENDIF

!---bottom fluxes
IF (iopt_sed_mode.GE.2) THEN
   DEALLOCATE (bottom_sed_dep,bottom_sed_ero,bottom_sed_flux)
ENDIF
DEALLOCATE (height_c,bed_fraction)
IF (iopt_sed_nodim.EQ.2) DEALLOCATE (t_equil)

!---bottom stress related arrays
DEALLOCATE (bdragcoefatc_sed,bdragcoefatu_sed,bdragcoefatv_sed,bstres_cr,&
          & bstresatc_sed,bstresatu_sed,bstresatv_sed,d50_bed,ubstresatu_sed,&
          & vbstresatv_sed,zroughatc_sed,zroughatu_sed,zroughatv_sed)
IF (iopt_waves.GT.0) THEN
   DEALLOCATE (bstresatc_mean_sed,bstresatc_wav_sed,fwave_sed,wavethickatc_sed,&
             & zaroughatc_sed)
ENDIF
   
!---diffusion coefficients
DEALLOCATE (beta_sed,vdiffcoef_sed)

!---bed slope arrays
IF (iopt_sed_slope.GT.0.OR.iopt_morph_avalanching.GT.0) THEN
   DEALLOCATE (bed_slope_x,bed_slope_y,bed_slope_x_atu,bed_slope_y_atu,&
             & bed_slope_x_atv,bed_slope_y_atv)
ENDIF

!---density
DEALLOCATE (beta_state_sed,densw)

!---fall velocity
DEALLOCATE (wfall)

!---open boundary conditions
IF (iopt_obc_sed.EQ.1) THEN
   IF (iopt_sed.EQ.1) THEN
      DEALLOCATE (obcsedatu,obcsedatv)
   ELSEIF (iopt_sed.EQ.2) THEN
      DEALLOCATE (obccnumpatu,obccnumpatv)
   ENDIF
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE deallocate_sed_arrays
