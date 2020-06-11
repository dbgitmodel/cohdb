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

MODULE darpars
!************************************************************************
!
! *morphpars* Dredging and relocation parameters
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)morphpars.f90  V2.8
!
! $Date: 2015-04-10 11:30:34 +0200 (Fri, 10 Apr 2015) $
!
! $Revision: 841 $
!
! Description - 
!
!************************************************************************
!
!

INTEGER, PARAMETER :: MaxCoordPol = 50, MaxFractions = 50, MaxSites = 50, &
                    & MaxSubzones = 50, MaxTimes = 50

INTEGER:: maxtrackpoints = 1000, nocampaigns = 1, &
        & nocirclepoints = 16, nodredgingsites = 1, &
        & norelocationsites = 1, noships = 1 

SAVE

! Name               Type     Purpose
!------------------------------------------------------------------------------
!*MaxCoordPol*       INTEGER Maximum allowed number of polygon points used to
!                            delineate a (sub)-dredging or (sub)-relocation site
!*MaxFractions*      INTEGER Maximum allowed number of fractions used in
!                            sediment loss terms at dredging or relocation sites
!*MaxSites*          INTEGER Maximum allowed of relocation sites linked to a
!                            campaign or dredging site
!*MaxSubzones*       INTEGER Maximum allowed num of dredging/relocation
!                            sub-zones
!*MaxTimes*          INTEGER Maximum allowed number of coupling times between
!                            dredging and relocation times
!*maxtrackpoints*    INTEGER The maximum number of cells on the dredging track
!                            after conversion
!*nocampaigns*       INTEGER Number of campaigns
!*nocirclepoints*    INTEGER Number of points that make up the circle polygon in
!                            case of deepest point relocation
!*nodredgingsites*   INTEGER Number of dredging sites
!*norelocationsites* INTEGER Number of relocation sites
!*noships            INTEGER Number of ships available for dredging and/or
!                            relocation
!************************************************************************


END MODULE darpars
