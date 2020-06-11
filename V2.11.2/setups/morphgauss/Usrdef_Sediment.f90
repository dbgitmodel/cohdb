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
! $Date: 2016-10-26 15:58:08 +0200 (Wed, 26 Oct 2016) $
!
! $Revision: 983 $
!
! Description - test case morphgauss
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
! Description - test case morphgauss
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE sedswitches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_sed_params'
CALL log_timer_in()

!
!1. Sediment switches
!--------------------
!
!1.1 General
!-----------
!
!---transport mode
iopt_sed_mode = 4
!---type of sediment
iopt_sed_type = 1

!
!1.2 Load formula
!----------------
!
!---bed load
SELECT CASE (runtitle(11:11))
  CASE('B')
     iopt_sed_toteq = 3
  CASE('C')
     iopt_sed_toteq = 4
  CASE('D')
     iopt_sed_toteq = 5
  CASE('E')
     iopt_sed_toteq = 6
  CASE('F')
     iopt_sed_toteq = 7
  CASE DEFAULT
     iopt_sed_toteq = 1
END SELECT

!
!1.3 Bottom stress
!-----------------
!
!---formulation for critical shear stress (0/1/2/3)
iopt_sed_bstres_cr = 2

!
!1.4 Settling
!------------
!
!---type of settling velocity
iopt_sed_ws = 2

!
!1.6 Open boundaries
!-------------------
!

iopt_sed_obc_flux = 1

CALL log_timer_out()

RETURN


END SUBROUTINE usrdef_sed_params

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
! Description - test case morphgauss
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
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_morph_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!

SELECT CASE (runtitle(11:11))
   CASE ('G')
      iopt_morph_fixed_layer = 1
      iopt_morph_active_layer = 0
   CASE ('H')
      iopt_morph_fixed_layer = 1
      iopt_morph_active_layer = 1
   CASE DEFAULT
      iopt_morph_fixed_layer = 0
      iopt_morph_active_layer = 0
END SELECT

SELECT CASE (runtitle(11:11))
   CASE ('I')
      iopt_morph_time_int = 2
      iopt_morph_hdiff = 1
   CASE ('J')
      iopt_morph_time_int = 3
      iopt_morph_hdiff = 0
   CASE DEFAULT
      iopt_morph_time_int = 1
      iopt_morph_hdiff = 0
END SELECT

!
!2. Morphology parameters
!------------------------
!
!---morphological acceleration
morph_factor = 10.0

!---bed porosity
bed_porosity_cst = 0.4

CALL log_timer_out()

RETURN


END SUBROUTINE usrdef_morph_params

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

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_dar_params

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

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_sedics

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
! Description - test case morphgauss
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE gridpars
USE iopars
USE morpharrays
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_morphics'
CALL log_timer_in()


!---layer thickness
bed_layer_thickness(1:ncloc,1:nrloc,1) = 5.0

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_morphics

!========================================================================

SUBROUTINE usrdef_sed_spec
!************************************************************************
!
! *usrdef_sed_charac* Define characteristics of sediment fractions
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.8
!
! Description - test case morphgauss
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

!---particle properties
dp = 150.0E-06
rhos = 2650.0

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_sed_spec

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

END SUBROUTINE usrdef_dar_spec
