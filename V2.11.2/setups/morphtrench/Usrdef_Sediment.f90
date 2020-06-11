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
! Author -
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.10
!
! Description - test case morphtrench
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
! Description - test case morphtrench
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE sedpars
USE sedswitches
USE gridpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

procname(pglev+1) = 'usrdef_sed_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!
!---general
iopt_sed_type = 1
SELECT CASE (runtitle(13:13))
   CASE ('A','E')
      iopt_sed_mode = 1
   CASE DEFAULT
      iopt_sed_mode = 3
      iopt_sed_bbc_eq = 2
END SELECT

!---bed load
iopt_sed_bedeq = 2

!---bed slope
SELECT CASE (runtitle(13:13))
   CASE ('C','G')
     iopt_sed_slope = 1
   CASE DEFAULT
     iopt_sed_slope = 0
END SELECT

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
! Description - test case morphtrench
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
dp = 10.0E-06

!---density
rhos = 2650.0

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_sed_spec

!========================================================================

SUBROUTINE usrdef_sedics
!************************************************************************
!
! *usrdef_sed_ics* Define arrays for initial conditions sediment transport
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

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_sedics

!========================================================================

SUBROUTINE usrdef_morph_params
!************************************************************************
!
! *usrdef_morph_params* Define parameters for the morphological model
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.10
!
! Description - test case morphtrench
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE morphpars
USE morphswitches
USe sedpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_morph_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!
!---avalanching
SELECT CASE (runtitle(13:13))
   CASE ('D','H')
      iopt_morph_avalanching = 1
   CASE DEFAULT
      iopt_morph_avalanching = 0
END SELECT

!---fixed layer
iopt_morph_fixed_layer = 1

!---enable mass correction
iopt_morph_corr = 1

!
!2. Morphological parameters
!----------------------------
!
!---morphological acceleration
morph_factor  = 100.0

!---bed porosity
bed_porosity_cst = 0.4

!---critical slopes for avalanching
stat_angle_wet = 4.0
dyn_angle_wet = 3.0

!---time counters
icmorph = icsed
icaval = 120

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_morph_params

!========================================================================

SUBROUTINE usrdef_morphics
!************************************************************************
!
! *usrdef_morphics* Define initial conditions for the morphological model
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.10
!
! Description - test case morphtrench
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE morpharrays
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_morph_ics'
CALL log_timer_in()

!---layer thickness
bed_layer_thickness(:,:,1) = 5.0

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_morphics

!========================================================================

SUBROUTINE usrdef_dar_params
!************************************************************************
!
! *usrdef_dar_params* Define parameters for dredging and relocation
!
! Author - Alexander Breugem and  Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.5
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

END SUBROUTINE usrdef_dar_params

!========================================================================

SUBROUTINE usrdef_dar_spec
!************************************************************************
!
! *usrdef_dar_spec* Define site and ship properties for dredging and relocation
!
! Author - Alexnader Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.5
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


