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
! $Date: 2016-02-18 12:05:20 +0100 (Thu, 18 Feb 2016) $
!
! $Revision: 908 $
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.10
!
! Description - test case morphctb
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
! Description - test case morphctb
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
iopt_sed_type = 1

!---load formulae
SELECT CASE (runtitle(10:10))
   CASE ('A','B')
      iopt_sed_mode = 4
      iopt_sed_toteq = 1
   CASE ('C','D','E','F','G','H')
      iopt_sed_mode = 3
      iopt_sed_bedeq = 1
END SELECT

!---formulation for critical shear stress
iopt_sed_bstres_cr = 2

!---type of settling velocity
iopt_sed_ws = 2

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_sed_params

!========================================================================

SUBROUTINE usrdef_sed_spec
!************************************************************************
!
! *usrdef_sed_spec* Define characteristics of sediment fractions
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.10
!
! Description - test case morphctb
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


procname(pglev+1) = 'usrdef_sed_spec'
CALL log_timer_in()

!--diameter
dp = 50.0E-06

!---diameter
rhos = 1650.0

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
! Description - test case morphctb
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
USE sedpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_morph_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!
!---fixed layer
SELECT CASE (runtitle(10:10))
   CASE ('A','B','C','D')
      iopt_morph_fixed_layer = 0
      iopt_morph_active_layer = 0
   CASE ('E','F')
      iopt_morph_fixed_layer = 1
      iopt_morph_active_layer = 0
   CASE ('G','H')
      iopt_morph_fixed_layer = 1
      iopt_morph_active_layer = 1
END SELECT

!---mass correction
iopt_morph_corr = 1

!
!2. Morphological parameters
!---------------------------
!
!---morphological factor
SELECT CASE (runtitle(10:10))
  CASE('A','C','E','G')
    morph_factor  = 1.0
  CASE('B','D','F','H')
    morph_factor = 5.0
END SELECT

!---bed porosity
bed_porosity_cst = 0.6

!---time counters
icmorph = icsed
icsedbal = 1440

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_morph_params

!========================================================================

SUBROUTINE usrdef_morphics
!************************************************************************
!
! *usrdef_morphics* Define arrays for initial conditions sediment transport
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.10
!
! Description - test case morphctb
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
bed_layer_thickness(:,:,1) = 1.0

CALL log_timer_out()


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

END SUBROUTINE usrdef_dar_params

!========================================================================

SUBROUTINE usrdef_dar_spec
!************************************************************************
!
! *usrdef_dar_spec* Define site and ship properties for dredging and relocation
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

END SUBROUTINE usrdef_dar_spec
