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
! *Usrdef_Sediment* User-defined model setup
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.10
!
! $Date: 2016-10-27 15:53:28 +0200 (Thu, 27 Oct 2016) $
!
! Description - test case morphtidav
!
! Reference -
!
! Routines - usrdef_sed_params, usrdef_morph_params, usrdef_dar_params,
!            usrdef_sedics, usrdef_morphics, usrdef_sed_spec, usrdef_dar_spec
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_sed_params
!************************************************************************
!
! *usrdef_sed_params* Define parameters sediment transport
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.10
!
! Description - test case morphtidav
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE sedswitches
USE sedpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_sed_params'

CALL log_timer_in()

!
!1. Switches
!-----------
!
!---general
iopt_sed_type = 1   ! Sand

!---load formulae
iopt_sed_mode = 4   ! Total load
iopt_sed_toteq = 1

!---formulation for critical shear stress
iopt_sed_bstres_cr = 2

!---type of settling velocity
iopt_sed_ws = 2

!
!2. Sediment parameters
!----------------------
!
!---time counter
icsed = 60

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_sed_params

!========================================================================

SUBROUTINE usrdef_sed_spec
!************************************************************************
!
! *usrdef_sed_charac* Define characteristics of sediment fractions
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.10
!
! Description - test case morphtidav
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE sedarrays
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_sed_charac'
CALL log_timer_in()

!--diameter
dp = 150.0E-06

!---density
rhos = 2650.0

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_sed_spec

!========================================================================

SUBROUTINE usrdef_sedics
!************************************************************************
!
! *usrdef_sedics* Define arrays for initial conditions sediment transport
!
! Author -
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.5
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_sedics

!========================================================================

SUBROUTINE usrdef_morph_params
!************************************************************************
!
! *usrdef_morph_params* Define parameters sediment transport
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.10
!
! Description - test case morphtidav
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE timepars
USE morphpars
USE morphswitches
USE sedpars
USE sedswitches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_morph_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!

SELECT CASE (runtitle(12:12))
   CASE ('B','C','D')
      iopt_morph_tidal_scheme = 0
   CASE ('E')
      iopt_morph_tidal_scheme = 1
   CASE ('F')
      iopt_morph_tidal_scheme = 2
   CASE ('G')
      iopt_morph_tidal_scheme = 3
END SELECT

!
!2. Morphological parameters
!----------------------------
!
!---time parameters
SELECT CASE (runtitle(12:12))
   CASE ('B')
      nstep_hydro = 34560; morph_steps = 2; number_tidal_steps = 576
   CASE ('C')
      nstep_hydro = 17280; morph_steps = 4; number_tidal_steps = 288
   CASE ('D','E','F','G')
      nstep_hydro = 17280; morph_steps = 2; number_tidal_steps = 288
END SELECT

!---bed porosity
bed_porosity_cst = 0.6

!---time counter
icmorph = icsed

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_morph_params

!========================================================================

SUBROUTINE usrdef_morphics
!************************************************************************
!
! *usrdef_morphics* Define arrays for initial conditions sediment transport
!
! Author -
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_morphics

!========================================================================

SUBROUTINE usrdef_dar_params
!************************************************************************
!
! *usrdef_dar_params* Define parameters for dredging and relocation
!
! Author -
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_dar_params

!========================================================================

SUBROUTINE usrdef_dar_spec
!************************************************************************
!
! *usrdef_dar_spec* Define site and ship properties for dredging and relocation
!
! Author -
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.8
!
! Description - empty default file
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_dar_spec
