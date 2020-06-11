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

MODULE reset_particles
!************************************************************************
!
! *reset_particles* Reset model parameters and arrays for the particle model
!                   if needed
!
! Author - Valerie Duliere and Patrick Luyten
!
! Revision - @(COHERENS)reset_model.F90  V2.x
!
! $Date: 2016-08-09 17:11:24 +0200 (Tue, 09 Aug 2016) $
!
! $Revision: 951 $
!
! Description - routines are empty but need to be declared
!
! Reference -
!
! Routines - reset_out_ppars, reset_part_params, reset_partics
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS

!========================================================================

SUBROUTINE reset_out_ppars(end_reset)
!************************************************************************
!
! *reset_out_ppars* Reset parameters for partyicle output
!
! Author - Valerie Duliere and Patrick Luyten
!
! Revision - @(COHERENS)reset_particles.f90  V2.x
!
! Description - empty default routine
!
! Calling program - particle_trajects_init
!  
! Module calls -
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: end_reset
  
!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*end_reset* LOGICAL Resets end value of time resolution (%tlims(2)) to a
!                    value lower than or equal to nstep if .TRUE.
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE reset_out_ppars

!========================================================================

SUBROUTINE reset_part_params
!************************************************************************
!
! *reset_part_params* Reset particle model parameters
!
! Author - Valerie Duliere and Patrick Luyten
!
! Revision - @(COHERENS)reset_particles.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model, initialise_particle_model
!
! Module calls -
!
!************************************************************************
!


RETURN

END SUBROUTINE reset_part_params

!=====================================================================

SUBROUTINE reset_partics
!************************************************************************
!
! *reset_partics* Reset initial conditions for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)reset_particles.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE reset_partics


END MODULE reset_particles
