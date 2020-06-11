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

MODULE dartypes
!************************************************************************
!
! *dartypes* Derived type definitions for the dredging/relocation model
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)dartypes.f90  V2.8
!
! $Date: 2015-04-10 11:30:34 +0200 (Fri, 10 Apr 2015) $
!
! $Revision: 841 $
!
! Description - 
!
!************************************************************************
!
USE darpars
USE syspars

IMPLICIT NONE

!---attributes for campaigns
TYPE :: CampaignParams
   LOGICAL :: dredge, relocate
   INTEGER :: nr_dredgingsite, nr_relocationsite, nr_ship, &
            & nr_zones_coupled, t_lag, t_spread
   REAL :: campdepth, camplevel, campvolume
   REAL, DIMENSION(MaxFractions) :: campfractions
   CHARACTER (LEN=lentime) :: CampaignStartDateTime
   CHARACTER (LEN=lentime), DIMENSION(MaxTimes) :: CoupleTime
   INTEGER, DIMENSION(MaxSites) :: relocationsites
END TYPE CampaignParams

!---attributes for dredging sites
TYPE :: DredgingParams
   INTEGER :: monitoring_frequency, nr_coordpol, nr_ship, &
            & nr_zones_coupled, t_lag, t_spread 
   REAL :: critical_volume, over_depth, target_depth, dredging_thickness
   INTEGER, DIMENSION(MaxTimes) :: nr_coupled
   INTEGER, DIMENSION(MaxSites) :: relocationsites
   INTEGER, DIMENSION(MaxSubzones) :: priority_order
   CHARACTER (LEN=lentime), DIMENSION(MaxTimes) :: CoupleTime
   REAL, DIMENSION(MaxCoordPol) :: Xpol, Ypol
END TYPE DredgingParams

!---attributes for dredging/relocation operations
TYPE :: OperationTimesParams
   INTEGER :: campaign_end
   INTEGER, DIMENSION(MaxTimes) :: coupling_times, dredging_times, &
                                 & relocation_times
END TYPE OperationTimesParams

!---attributes for relocation sites
TYPE :: RelocationParams
   INTEGER :: nr_coordpol 
   LOGICAL :: eligibility
   REAL :: max_volume, min_depth
   REAL, DIMENSION(MaxCoordPol) :: Xpol, Ypol
END TYPE RelocationParams

!---attributes for dredging/relocation ships
TYPE :: ShipParams
   INTEGER :: t_fill, t_idle, t_relocate,  t_sail_dredge, t_sail_relocate, &
        & t_wait
   REAL :: D_dredge_max, D_dredge_min, D_relocate_max, &
         & D_relocate_min, H_dredge_max, H_dredge_min, &
         & H_relocate_max, H_relocate_min, R_circle, U_ship, V_ship
   REAL, DIMENSION(MaxFractions) :: p_overflow, p_passive, p_relocation, &
                                  & p_suctionhead      
END TYPE ShipParams

!---attributes for dredging sub-zones
TYPE :: SubzonesDredgingParams
   INTEGER:: nr_subzones
   INTEGER, DIMENSION(MaxSubzones) :: nr_coordpol
   REAL, DIMENSION(MaxCoordPol,MaxSubzones) :: Xpol, Ypol
END TYPE SubzonesDredgingParams

!---attributes for relocation sub-zones
TYPE :: SubzonesRelocationParams
   INTEGER:: nr_subzones
   INTEGER, DIMENSION(MaxSubzones) :: nr_coordpol
   REAL, DIMENSION(MaxCoordPol,MaxSubzones) :: Xpol, Ypol
END TYPE SubzonesRelocationParams

SAVE

! Name               Type    Purpose
!------------------------------------------------------------------------------
!*Campaignparams*         DERIVED Attributes for campaigns
!*%dredge*                LOGICAL Disables/enables dredging during the campaign
!*%relocate*              LOGICAL Disables/enables relocation during the
!                                 campaign
!*%nr_dredging_site*      INTEGER Number of dredging sites coupled to the
!                                 campaign
!*%nr_relocation_site*    INTEGER Number of relocation sites coupled to the
!                                 campaign
!*%nr_ship*               INTEGER Ship id of the campaign
!*%nr_zones_coupled*      INTEGER Number of relocation sites coupled to the
!                                 campaign
!*%t_lag*                 INTEGER Time lag between dredging and relocation
!*%t_spread*              INTEGER Duration (in time steps) of the dredging
!                                 works in case of time spread dredging
!*%campdepth*             REAL    Dredging depth                             [m]
!*%camplevel*             REAL    Dredging level                             [m]
!*%campvolume*            REAL    Dredging volume of the campaign, taken evenly
!                                 across the dredging site                   [m]
!*%campfractions*         REAL    Size fraction distribution of the relocated
!                                 sediment
!*%CampaignStartDateTime* CHAR    Campain start date/time
!*%CoupleTime*            CHAR    Times at which coupling between dredging and
!                                 relocation sites is defined
!*%relocationsites*       INTEGER Relocation sites coupled to the campaign
!
!*DredgingParams*         DERIVED Attributes for dredging sites
!*%monitoring_frequency*  INTEGER Frequency (in number of time steps at which
!                                 bathymetry is monitored and the dredging
!                                 criterion is checked
!*%nr_coordpol*           INTEGER Number of points marking the dredging site
!                                 polygon
!*%nr_ship*               INTEGER Id of the ship connected to the dredging site
!*%nr_zones_coupled*      INTEGER Number of zones in the dredging site
!*%t_lag*                 INTEGER Time lag (in number of time steps) between
!                                 dredging and relocation
!*%t_spread*              INTEGER Duration (in number of time steps) of the
!                                 dredging in case of time spread dredging
!*%critical_volume*       REAL    Critical dredging volume that triggers a
!                                 dredging operation
!*%over_depth*            REAL    Dredging overdepth                         [m]
!*%target_depth*          REAL    Dredging taerget depth                     [m]
!*%dredging_thickness*    REAL    Thichness of the layer that a ship can dredge
!                                 during one pass (in case of track dredging)[m]
!*%nr_coupled*            INTEGER Number of coupling times (in case of track
!                                 dredging)
!*%relocationsites*       INTEGER Relocation sites coupled to the dredging site
!*%priority_order*        INTEGER Order in which the dredging subzones are
!                                 given priority
!*%CoupleTime*            CHAR    Times at which coupling between dredging and
!                                 relocation sites are defined (in case of
!                                 user-defined spatial coupling)
!*%Xpol*                  REAL    X-coordinates of the dredging site polygon
!                                                                 [m or degrees]
!*%Ypol*                  REAL    Y-coordinates of the dredging site polygon
!                                                                 [m or degrees]
!
!*OperationTimesParams*   DERIVED Time parameters for dredging/relocation
!                                 operations
!*%campaign_end*          INTEGER End (in number of time steps since the start
!                                 of the simulation) of a dredging/relocation
!                                 campaign (in absence of track dredging)
!*%coupling_times*        INTEGER Time (in number of time steps since the start
!                                 of the simulation) at which coupling btween
!                                 dredging and relocation sites is defined
!*%dredging_times*        INTEGER Times (in number of time steps since the start
!                                 of the simulation) for the start/end of
!                                 dredging during a dredging/relocation
!                                 campaign
!*%relocation_times*      INTEGER Times (in number of time steps since the start
!                                 of the simulation) for the start/end of
!                                 relocation during a dredging/relocation
!                                 campaign
!*RelocationParams*       DERIVED Attributes for relocation sites
!*%nr_coordpol*           INTEGER Number of points marking the relocation site
!                                 polygon
!*%eligibility*           LOGICAL Indicates whether any relocation site is
!                                 still eligible for relocation (in case of an
!                                 active relocation criterion)
!*%max_volume*            REAL    Maximum allowed volume of sediment which can
!                                 be relocated (in case of a maximum volume
!                                 relocation criterion)                    [m^3]
!*%min_depth*             REAL    Minimum depth that needs to be assured at a
!                                 relocation site (in case of a minimum depth
!                                 relocation criterion)                      [m]
!*%Xpol*                  REAL    X-coordinates of the relocation site polygon
!                                                                 [m or degrees]
!*%Ypol*                  REAL    Y-coordinates of the relocation site polygon
!                                                                 [m or degrees]
!*ShipParams*             DERIVED Attributes for ship operations during
!                                 dredging/relocation campaigns
!*%t_fill*                INTEGER Time (in number of time steps) needed to fill
!                                 the ship
!*%t_idle*                INTEGER Idle time (in number of time steps) of the
!                                 ship
!*%t_relocate*            INTEGER Time (in number of time steps) needed to
!                                 empty the ship at a relocation site
!*%t_sail_dredge*         INTEGER Time (in number of time steps) needed to
!                                 sail to the dredging site
!*%t_sail_relocate*       INTEGER Time (in number of time steps) needed to sail
!                                 from the dredging to the relocation site
!*%t_wait*                INTEGER Waiting time (in number of time steps) before
!                                 the start of the dredging/relocation campaign
!*%D_dredge_max*          REAL    Maximum depth (with respect to the sea
!                                 surface) for overloss of dredged sediment
!                                 at the dredging site                       [m]
!*%D_dredge_min*          REAL    Minimum depth (with respect to the sea
!                                 surface) for overloss of dredged sediment
!                                 at the dredging site                       [m]
!*%D_relocate_max*        REAL    Maximum depth (with respect to the sea
!                                 surface) for overflow of sediment by the ship
!                                 during relocation                          [m]
!*%D_relocate_min*        REAL    Minimum depth (with respect to the sea
!                                 surface) for overflow of sediment by the ship
!                                 during relocation                          [m]
!*%H_dredge_max*          REAL    Maximum height (above to the sea bed) for
!                                 suctionhead losses of dredged sediment at the
!                                 dredging site                              [m]
!*%H_dredge_min*          REAL    Minimum height (above to the sea bed) for
!                                 suctionhead losses of dredged sediment at the
!                                 dredging site                              [m]
!*%H_relocate_max*        REAL    Maximum height (above to the sea bed) for
!                                 passive losses of dredged sediment at the
!                                 relocation site                            [m]
!*%H_relocate_min*        REAL    Minimum height (above to the sea bed) for
!                                 passive losses of dredged sediment at the
!                                 relocation site                            [m]
!*%R_circle*              REAL    Radius of the circular area where relocation
!                                 takes place in case of deepest point
!                                 relocation                                 [m]
!*%U_ship*                REAL    The ship's velocity in case of track dredging
!                                                                          [m/s]
!*%V_ship*                REAL    Volume amount of sediment which can be
!                                 carried by the ship                      [m^3]
!*%p_overflow*            REAL    Fractional loss of dredged sediment by
!                                 overflow from the ship at the dredging site
!*%p_passive*             REAL    Fractional loss of relocated sediment leading
!                                 to passive plume formation near the bed at
!                                 the relocation site
!*%p_relocation*          REAL    Fractional loss of relocated sediment by the
!                                 ship at the relocation site
!*%p_suctionhead*         REAL    Fractional suctionhead losses of dredged
!                                 sediment near the bed at the dredging site
!
!*SubzonesDredgingParams* DERIVED Attributes for subzones at the dredgingsite
!*%nr_subzones*           INTEGER Number of subzones at the dredging site
!*%nr_coordpol*           INTEGER Number of points marking the sub-zone
!                                 dredging polygons
!*%Xpol*                  REAL    X-coordinates for the subzone dredging
!                                 polygons                        [m or degrees]
!*%Ypol*                  REAL    Y-coordinates for the subzone dredging
!                                 polygons                        [m or degrees]
!
!*SubzonesRelocationParams*DERIVED Attributes of subzones at the relocation
!                                  site
!*%nr_subzones*           INTEGER Number of subzones at the relocation site
!*%nr_coordpol*           INTEGER Number of points marking the sub-zone
!                                 relocation polygons
!*%Xpol*                  REAL    X-coordinates of the subzone relocation
!                                 polygons                        [m or degrees]
!*%Ypol*                  REAL    Y-coordinates of the subzone relocation
!                                 polygons                        [m or degrees]
!
!*******************************************************************************
!

END MODULE dartypes
