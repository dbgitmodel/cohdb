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

MODULE partswitches
!************************************************************************
!
! *partswitches* Switches for the Lagrangian particle module
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)partswitches.f90  V2.11
!  
! $Date: 2017-04-12 16:52:46 +0200 (Wed, 12 Apr 2017) $
!
! $Revision: 1012 $
!
! Description - 
!
!************************************************************************
!
IMPLICIT NONE

!---advection
INTEGER :: iopt_part_hadv = 2, iopt_part_vadv = 1

!---sign change in leeway drift
INTEGER :: iopt_part_leeway = 0

!---buoyancy
INTEGER :: iopt_part_dens = 0

!---diffusion
INTEGER :: iopt_part_hdif = 0, iopt_part_vdif = 0

!---particle cloud generator
INTEGER :: iopt_part_cloud = 1

!---use wind current to drive particle drift
INTEGER :: iopt_part_wind = 1

!---volume concentrations
INTEGER :: iopt_part_conc = 0

!---output from particle module
INTEGER :: iopt_part_out = 1

SAVE

!
! Name                  Type    Purpose
!-------------------------------------------------------------------------------
!*iopt_part_cloud*      INTEGER Disable/enable the generation of particle
!                               clouds (0/1)
!*iopt_part_conc*       INTEGER Select type of Eulerian concentrations for the
!                               particle module
!                        = 0 => No continuous particle distributions
!                        = 1 => Number densities and age concentrations
!                        = 2 => Volume and age concentrations
!                        = 3 => Mass and age concentrations
!*iopt_part_dens*       INTEGER Buoyancy effects for the vertical particle
!                               velocity
!                        = 0 => Buoyancy effects excluded
!                        = 1 => Using reference density
!                        = 2 => Using density from equation of state (obtained
!                               from COHERENS)
!*iopt_part_hadv*       INTEGER Type of time integration scheme for horizontal
!                               advection
!                        = 0 => Horizontal advection disabled
!                        = 1 => First order forward Euler
!                        = 2 => Fourth order Runge-Kutta 
!*iopt_part_hdif*       INTEGER Disables/enables horizontal diffusion (0/1)
!*iopt_part_leeway*     INTEGER Disables/enables leeway drift to change sign
!                               randomnly (0/1)
!*iopt_part_out*        INTEGER Disables/enables particle trajectory output
!                               (0/1)
!*iopt_part_vadv*       INTEGER Disables/enables vertical advection (0/1)
!*iopt_part_vdif*       INTEGER Selects type of vertical diffusion coefficient
!                        = 0 => Vertical diffusion disabled
!                        = 1 => Using constant diffusion coefficient
!                        = 2 => Using diffusion coefficient from COHERENS 
!*iopt_part_wind*       INTEGER Disables/enables wind input (0/1)
!
!************************************************************************
!

END MODULE partswitches
