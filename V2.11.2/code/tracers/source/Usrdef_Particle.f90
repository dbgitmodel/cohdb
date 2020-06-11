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
! *Usrdef_Particle* User-defined setup for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description - empty default file
!
! Reference -
!
! Routines - usrdef_part_params, usrdef_partics, usrdef_part_output,
!            usrdef_part_spec, usrdef_pout_params, usrdef_parcld_data
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_part_params
!************************************************************************
!
! *usrdef_part_params* Define parameters for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11
!
! Description - empty default file
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_part_params

!========================================================================

SUBROUTINE usrdef_partics
!************************************************************************
!
! *usrdef_partics* Define initial conditions for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11
!
! Description - empty default file
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_partics

!========================================================================

SUBROUTINE usrdef_part_output
!************************************************************************
!
! *usrdef_part_output* User-formatted output for the particle model 
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11
!
! Description - empty default file
!
! Reference -
!
! Calling program - coherens_main, initialise_model, initialise_particle_model,
!                   particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_part_output

!========================================================================

SUBROUTINE usrdef_part_spec
!************************************************************************
!
! *usrdef_part_spec* Define specifier arrays for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11
!
! Description - empty default file
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_part_spec

!========================================================================

SUBROUTINE usrdef_pout_params
!************************************************************************
!
! *usrdef_pout_params* Define parameters for particle trajectory output
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11
!
! Description - empty default file
!
! Reference -
!
! Calling program - particle_trajects_init
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_pout_params

!========================================================================

SUBROUTINE usrdef_parcld_data(ifil,ciodatetime,disdata)
!************************************************************************
!
! *usrdef_parcld_data* Define locations for cloud discharges
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11
!
! Description - empty default file
!
! Reference -
!
! Calling program - define_parcld_data
!
!************************************************************************
!
USE syspars
  
IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(OUT) :: ciodatetime
INTEGER, INTENT(IN) :: ifil
REAL (KIND=kndrlong), INTENT(INOUT), DIMENSION(3) :: disdata

!
! Name         Type Purpose
!------------------------------------------------------------------------------
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*disdata*     REAL    Input data
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_parcld_data
