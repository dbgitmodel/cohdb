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

MODULE morphswitches
!************************************************************************
!
! *morphswitches* Morphological switches
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)sediment.f90  V2.10
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
 
IMPLICIT NONE

INTEGER :: iopt_morph_active_layer = 1, iopt_morph_avalanching = 0, &
         & iopt_morph_corr = 0, iopt_morph_fixed_layer = 0, &
         & iopt_morph_hdiff = 0, iopt_morph_tidal_scheme = 1, &
         & iopt_morph_time_int = 1, iopt_morph_vert_fluxes = 0, &

! Name              Type    Purpose
!------------------------------------------------------------------------------
!
!*iopt_morph_active_layer* INTEGER Select the equation that determines the!
!                                  active layer thickness
!                            = 0 => Disabled
!                            = 1 => Flat bed concitions, grain-size dependent
!                                   (Armanini)
!                            = 2 => In the presence of bedforms (Ribberink)
!                            = 3 => Combination of sediment features and
!                                   hydraulic  conditions (Harris & Wiberg)
!*iopt_morph_avalanching*  INTEGER Disables/enables avalanching (0/1)
!*iopt_morph_corr*         INTEGER Selects type of action taken for mass
!                                  conservation when more sediment is taken
!                                  from the bed layer than available
!                            = 0 => No action is taken a monitoring message is
!                                   written (if enabled)
!                            = 1 => Excess mass is taken from neighbouring
!                                   cells and a monitoring message
!                                   (if enabled) is written if the scheme fails)
!                            = 2 => Excess mass is taken from neighbouring
!                                   cells and the program is aborted with an
!                                   error message
!*iopt_morph_fixed_layer*  INTEGER Disables/enables fixed layers in the bed!
!                                  (0/1)
!*iopt_morph_hdiff*        INTEGER Disables/enables the numerical diffusion term
!                                  in the Exner equation (0/1)
!*iopt_morph_tidal_scheme* INTEGER Select the sceheme for tidal averaging
!                                  velocity correction 
!                            = 0 => No correction
!                            = 1 => Continuity correction
!                            = 2 => Friction correction
!                            = 3 => Mixed continuity/friction correction
!*iopt_morph_time_int*     INTEGER Select the kind of time integration
!                            = 1 => Explicit Euler (first order)
!                            = 2 => Adams-Bashforth (second order)
!                            = 3 => Adams-Bashforth (third order)
!*iopt_morph_vert_fluxes*   INTEGER Select the kind of vertical exchange model
!                             = 0 => Disabled
!                             = 1 => two-layer model (exchange layer concept
!                                    Ribberink)
!
!*******************************************************************************
!

SAVE

END MODULE morphswitches
