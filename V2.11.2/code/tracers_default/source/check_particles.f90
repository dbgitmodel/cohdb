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

MODULE check_particles
!************************************************************************
!
! *check_particles* Check parameters for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!  
! Version - @(COHERENS)check_particles.f90  V2.x
!
! $Date: 2016-07-12 15:51:17 +0200 (Tue, 12 Jul 2016) $
!
! $Revision: 943 $
!
! Description - routines are empty but need to be declared
!
! Reference -
!
! Routines - check_out_ppars, check_partics, check_part_params
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS

!=====================================================================

SUBROUTINE check_out_ppars
!************************************************************************
!
! *check_out_ppars* check parameters for trajectory output
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)check_particles.f90  V2.x
!
! Description - empty default routine
!
! Reference -
!
! Calling program - particle_trajects_init
!
! Module calls - 
!
!************************************************************************
!


RETURN

END SUBROUTINE check_out_ppars

!=====================================================================

SUBROUTINE check_partics
!************************************************************************
!
! *check_partics* Check initial conditions for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)check_particles.f90  V2.x
!
! Description - routines are empty but need to be declared
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
! Module calls -
!
!************************************************************************
!


RETURN

END SUBROUTINE check_partics

!=====================================================================

SUBROUTINE check_part_params
!************************************************************************
!
! *check_part_params* Check parameters for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)check_particles.f90  V2.x
!
! Description - routines are empty but need to be declared
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
! Module calls -
!
!************************************************************************
!


RETURN

END SUBROUTINE check_part_params

END MODULE check_particles
