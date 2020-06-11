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

MODULE darswitches
!************************************************************************
!
! *morphswitches* Dredging and relocation switches
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)sediment.f90  V2.5
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
 
IMPLICIT NONE

INTEGER ::  iopt_dar_coupling_time = 1, iopt_dar_coupling_space = 1, &
          & iopt_dar_dredging_criterion = 1, &
          & iopt_dar_duration = 1, iopt_dar_effect= 1, &
          & iopt_dar_relocation_criterion = 1, iopt_dar_scenario = 1, &
          & iopt_dar_time_dredge = 1, iopt_dar_time_relocate = 1, &
          & iopt_dar_user_volume = 1

SAVE

! Name                        Type    Purpose
!------------------------------------------------------------------------------
!
!*iopt_dar_coupling_time*     INTEGER Selects the type of temporal coupling
!                                     between dredging and relocation
!                               = 1 => Simultaneous
!                               = 2 => Time lagged
!                               = 3 => Determined by full dredging cycle
!*iopt_dar_coupling_space*    INTEGER Selects the type of spatial coupling 
!                                     between dredging and relocation
!                               = 1 => Cyclical order, changed per dredging
!                                      cycle
!                               = 2 => Sequential order, changed per relocation
!                                      criterion
!                               = 3 => Time series
!*iopt_dar_dredging_criterion*INTEGER Selects the type of dredging criterion
!                               = 1 => Target depth
!                               = 2 => Critical volume
!                               = 3 => User-defined dredging and/or relocation
!*iopt_dar_duration*          INTEGER Selects the duration of the dredging 
!                                     and relocation activities
!                               = 1 => Instantaneous
!                               = 2 => Spread over time
!                               = 3 => Spread over cycles
!*iopt_dar_effect*            INTEGER Selects the type of effect being modeled
!                               = 1 => Morphology
!                               = 2 => Plume modelling
!                               = 3 => Both
!*iopt_dar_relocation_criterion*INTEGER Determines whether or not relocation 
!                                       is allowed
!                               = 0 => disabled
!                               = 1 => minimum water-depth (user-defined)
!                               = 2 => maximum volume of relocated sediment
!*iopt_dar_scenario*          INTEGER Selects the type of scenario being
!                                     modeled (=short hand switch)
!                               = 1 => LT (long term)
!                               = 2 => ST (short term)
!                               = 3 => PD (plume dispersal)
!*iopt_dar_time_dredge*       INTEGER Selects the kind of time integration
!                               = 1 => Equal spreading
!                               = 2 => Predefined order
!                               = 3 => Dredging track
!*iopt_dar_time_relocate*     INTEGER Selects the kind of time integration
!                               = 1 => Equal spreading
!                               = 2 => Deepest point
!                               = 3 => Random order
!*iopt_dar_user_volume*       INTEGER Determines the way volumes are defined
!                                     in user-defined campaigns
!                               = 1 => Dredging depth
!                               = 2 => Total volume
!                               = 3 => Reference level
!
!******************************************************************************
!

END MODULE darswitches
