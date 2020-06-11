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

MODULE dararrays
!************************************************************************
!
! *dararrays* Dredging and relocation arrays
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)morphology_arrays.f90  V2.8
!
! $Date: 2016-02-18 12:05:20 +0100 (Thu, 18 Feb 2016) $
!
! $Revision: 908 $
!
! Description - 
!
!************************************************************************
!
!
USE dartypes 

!---derived types
TYPE (CampaignParams), ALLOCATABLE, DIMENSION(:) :: Campaigns
TYPE (DredgingParams), ALLOCATABLE, DIMENSION(:) :: DredgingSites
TYPE (OperationTimesParams), ALLOCATABLE, DIMENSION(:) :: Operationtimes
TYPE (RelocationParams), ALLOCATABLE, DIMENSION(:) :: RelocationSites
TYPE (ShipParams), ALLOCATABLE, DIMENSION(:) :: Ships
TYPE (SubzonesDredgingParams), ALLOCATABLE, DIMENSION(:) :: SubzonesDredging
TYPE (SubzonesRelocationParams), ALLOCATABLE, DIMENSION(:) :: SubzonesRelocation

!---logical arrays
LOGICAL, ALLOCATABLE, DIMENSION(:) :: campaign, dredging, eligibility, &
                                    & relocation
LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: mask_dredging, mask_relocation 
LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: mask_subzones_dredging, &
                                          & mask_subzones_relocation

!---integer arrays
INTEGER, ALLOCATABLE, DIMENSION(:) :: campaign_start, couple_counter, &
     & cycle_counter, dredge_counter, nocouplingtimes, notrackpoints, &
     & nr_cycles, nt_cycle, nt_dredge, nt_dredge_lastcycle, nt_fill, nt_idle, &
     & nt_lag, nt_passed, nt_relocate, nt_relocate_lastcycle, nt_sail_dredge, &
     & nt_sail_relocate, nt_spread, nt_wait, priority_counter, &
     & relocate_counter, spatial_counter, track_cell_counter, track_nt_counter,&
     & track_points_counter
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: itrack, jtrack, nt_dredge_subzone, &
                                      & nt_track

!---real arrays
REAL, ALLOCATABLE, DIMENSION(:) :: D_dredge_max, D_dredge_min, D_relocate_max, &
     & D_relocate_min, H_dredge_max, H_dredge_min, H_relocate_max, &
     & H_relocate_min, volume_ship, V_dredge_total, V_lastcycle, &
     & V_relocate_extra, V_relocate_total, V_ship_rest
REAL, ALLOCATABLE, DIMENSION(:,:) :: fraction_ship, p_overflow, p_passive, &
                                   & p_relocation, p_suctionhead, xtrack, ytrack
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dredging_depth, V_dredge, &
                                     & V_dredge_max, V_dredge_rest


! Name                     Type    Purpose
!------------------------------------------------------------------------------
!*Campaigns*               DERIVED Attributes for campaigns
!*campaign*                LOGICAL Activates a campaign
!*campaign_start*          INTEGER Start time index for a campaign
!*couple_counter*          INTEGER Spatial counter giving the relocation site
!                                  in case of a user-defined spatial coupling 
!                                  through time series
!*cycle_counter            INTEGER Counter increased by one after each
!                                  dredging/relocation cycle
!*DredgingSites*           DERIVED Attributes for dredging sites
!*dredging*                LOGICAL Activates a dredging operation
!*dredge_counter*          INTEGER Counter increased by one after each dredging
!                                  activity
!*dredging_depth*          REAL    Depth of the sediment layer to be dredged [m]
!*D_dredge_max*            REAL    Maximum depth (with respect to the sea
!                                  surface) for overloss of dredged sediment
!                                  at the dredging site                      [m]
!*D_dredge_min*            REAL    Minimum depth (with respect to the sea
!                                  surface) for overloss of dredged sediment
!                                  at the dredging site                      [m]
!*D_relocate_max*          REAL    Maximum depth (with respect to the sea
!                                  surface) for overflow of sediment by the ship
!                                  during relocation                         [m]
!*D_relocate_min*          REAL    Minimum depth (with respect to the sea
!                                  surface) for overflow of sediment by the ship
!                                  during relocation                         [m]
!*eligibility*             LOGICAL .TRUE. if a relocation site is eligible for
!                                  relocation
!*fraction_ship*           REAL    Size fraction distribution of the relocated
!                                  sediment loaded on the ship
!*H_dredge_max*            REAL    Maximum height (above to the sea bed) for
!                                  suctionhead losses of dredged sediment at the
!                                  dredging site                             [m]
!*H_dredge_min*            REAL    Minimum height (above to the sea bed) for
!                                  suctionhead losses of dredged sediment at the
!                                  dredging site                             [m]
!*H_relocate_max*          REAL    Maximum height (above to the sea bed) for
!                                  passive losses of dredged sediment at the
!                                  relocation site                           [m]
!*H_relocate_min*          REAL    Minimum height (above to the sea bed) for
!                                  passive losses of dredged sediment at the
!                                  relocation site                           [m]
!*itrack*                  INTEGER X-indices of the grid cells used by track
!                                  dredging 
!*jtrack*                  INTEGER Y-indices of the grid cells used by track
!                                  dredging
!*mask_dredging*           LOGICAL .TRUE. at grid cells within a dredging area
!*mask_relocation*         LOGICAL .TRUE. at grid cells within a relocation area
!*mask_subzones_dredging*  LOGICAL .TRUE. at grid cells within a sub-zone
!                                  dredging area
!*mask_subzones_relocation*LOGICAL .TRUE. at grid cells within a sub-zone
!                                  relocation area
!*nocouplingtimes*         INTEGER Number of coupling times in case of a
!                                  user-defined spatial coupling through time
!                                  series
!*notrackpoints*           INTEGER Number of ship track positions for a
!                                  specific site
!*nr_cycles*               INTEGER Number of cycles in case dredging/relocation
!                                  is spread over cycles
!*nt_cycle*                INTEGER Duration (in time steps) a
!                                  dredging/relocation cycle
!*nt_dredge*               INTEGER Duration (in time steps) of a complete
!                                  dredging cycle at a dredging site
!*nt_dredge_lastcycle*     INTEGER Duration (in time steps) of dredging
!                                  during the last cycle
!*nt_dredge_subzone*       INTEGER Duration (in time indices) of dredging at a
!                                  dredging sub-zone
!*nt_fill*                 INTEGER Time (in time steps) needed to fill a ship
!*nt_idle*                 INTEGER Idle time (in time steps) of a ship
!*nt_lag*                  INTEGER Time lag (in time steps) between dredging
!                                  and relocation
!*nt_passed*               INTEGER Time counter increased by one after each
!                                  dredging time step within a sub-zone
!                                  dredging area 
!*nt_relocate*             INTEGER Time (in time steps) needed for relocation
!                                  during a dredging/relocation cycle 
!*nt_relocate_lastcycle*   INTEGER Time (in time steps) needed for relocation
!                                  during the last dredging/relocation cycle
!*nt_sail_dredge*          INTEGER Time (in time steps) needed to sail to the
!                                  dredging site
!*nt_sail_relocate*        INTEGER Time (in time steps) needed to sail to the
!                                  relocation site
!*nt_spread*               INTEGER Duration (in number of time steps) of the
!                                  dredging in case of time spread dredging
!*nt_track*                INTEGER Duration (in number of time steps) of the
!                                  dredging at a grid cell during tracked
!                                  dredging
!*nt_wait*                 INTEGER Waiting time (in number of time steps) before
!                                  the start of the dredging/relocation campaign
!*Operationtimes*          DERIVED Time parameters for dredging/relocation
!                                  operations
!*priority_counter*        INTEGER Counter increased by one when a dredging is
!                                  completed in a dredging sub-zone in case of
!                                  priority dredging 
!*p_overflow*              REAL    Fractional loss of dredged sediment by
!                                  overflow from the ship at the dredging site
!*p_passive*               REAL    Fractional loss of relocated sediment leading
!                                  to passive plume formation near the bed at
!                                  the relocation site
!*p_relocation*            REAL    Fractional loss of relocated sediment by the
!                                  ship at the relocation site
!*p_suctionhead*           REAL    Fractional suctionhead losses of dredged
!                                  sediment near the bed at the dredging site
!*RelocationSites*         DERIVED Attributes for relocation sites
!*relocation*              LOGICAL Activates a relocation operation
!*relocate_counter*        INTEGER Counter increased by one after each
!                                  relocation activity
!*Ships*                   DERIVED Attributes for ship operations during
!                                  dredging/relocation campaigns
!*SubzonesDredging*        DERIVED Attributes for subzones at a dredging site
!*SubzonesRelocation*      DERIVED Attributes for subzones at a relocation site
!*spatial_counter*         INTEGER Counter used during relocation to define the
!                                  sub-zone for relocation and increased by one
!                                  when relocation is completed at the previous
!                                  relocation sub-zone
!*track_cell_counter*      INTEGER Spatial counter increased by one each time
!                                  the ship moves to the next dredging along
!                                  its track 
!*track_nt_counter*        INTEGER Time counter increased by one after each
!                                  time step when the ship dredges a specific
!                                  cell along its track
!*track_points_counter*    INTEGER Number of cells passed by the ship during
!                                  tracked dredging
!*volume_ship*             REAL    Amount of sediments currently loaded on the
!                                  ship                                    [m^3]
!*V_dredge*                REAL    Amount of sediment to be dredged at a
!                                  dredging cell                           [m^3]
!*V_dredge_max*            REAL    Maximum amount of sediment to be dredged at
!                                  a specific cells during track dredging  [m^3]
!*V_dredge_rest*           REAL    Remaining amount of sediment to be dredged
!                                  at a specific cell during track dredging
!                                                                          [m^3]
!*V_dredge_total*          REAL    Total amount of sediment to be dredged at a
!                                  dredging site                           [m^3]
!*V_lastcycle*             REAL    Amount of sediment to be dredged at a
!                                  specific cell during the last cycle in case
!                                  of cycled dredging and relocation       [m^3]
!*V_relocate_extra*        REAL    Amount of non-relocated sediment at a
!                                  relocation site                         [m^3]
!*V_relocate_total*        REAL    Amount of relocated sediment at a
!                                  relocation site                         [m^3]
!*V_ship_rest*             REAL    Amount of sediment which can be loaded on
!                                  the ship during track dredging          [m^3]
!*xtrack*                  REAL    X-coordinates of the ship during tracked
!                                  dredging                       [m or degrees]
!*ytrack*                  REAL    Y-coordinates of the ship during tracked
!                                  dredging                       [m or degrees]
!
!******************************************************************************
!

SAVE

END MODULE dararrays
