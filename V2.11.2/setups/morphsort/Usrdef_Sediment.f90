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
! $Date: 2016-02-18 12:05:20 +0100 (Thu, 18 Feb 2016) $
!
! $Revision: 908 $
!
! Description - test case morphsort
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
! Description - test case morphsort
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
USE timepars
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
iopt_sed_mode = 1

!---type of sediment
iopt_sed_type = 1

!
!1.2 Load formula
!----------------
!
!---bed load (1/2/3/4/5/6/7)
SELECT CASE (runtitle(10:10))
   CASE ('H')
      iopt_sed_bedeq = 2
   CASE ('I')
      iopt_sed_bedeq = 4
   CASE ('J')
      iopt_sed_bedeq = 6
   CASE DEFAULT
      iopt_sed_bedeq = 1
END SELECT

!
!1.3 Bottom stress
!-----------------
!
!---formulation for critical shear stress
iopt_sed_bstres_cr = 2

!
!1.4 Settling
!------------
!
!---type of settling velocity
iopt_sed_ws = 2

!
!1.5 Open boundaries
!-------------------
!

iopt_sed_obc_flux = 1

!
!2. Sediment parameters
!----------------------
!
!---number of fractions
nf = 2
!---number of bed layers
nb = 8
!---sediment time step
icsed = ic3d

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
! Description - test case morphsort
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
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_morph_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!

iopt_morph_active_layer = MERGE(1,0,runtitle(10:10).NE.'D')
SELECT CASE (runtitle(10:10))
   CASE ('E')
      iopt_morph_time_int = 2
   CASE ('F','G')
      iopt_morph_time_int = 3
   CASE DEFAULT
      iopt_morph_time_int = 1
END SELECT
iopt_morph_hdiff = MERGE(0,1,runtitle(10:10).NE.'G')
iopt_morph_corr = MERGE(1,0,runtitle(10:10).NE.'C')

!
!2. Morphology parameters
!------------------------
!
!---morphological acceleration
morph_factor = MERGE(10.0,1.0,runtitle(10:10).NE.'B')

!---time counters
icmorph = ic3d
icsedbal = 360

!---bed porosity
bed_porosity_cst = 0.6

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_morph_params

!========================================================================

SUBROUTINE usrdef_dar_params
!************************************************************************
!
! *usrdef_dar_params* Define parameters for dredging and relocation
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
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
! Author - Alexander Breugem and Kevin Delecluse (IMDC)
!
! Last update - 4 Nov 2008  @(COHERENS)Usrdef_Model.f90  V2.10
!
! Description - test case morphsort
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE gridpars
USE iopars
USE sedarrays
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_sedics'
CALL log_timer_in()

!---bed fraction 1 in layers 1 and 3, fraction 2 in layer 2
bed_fraction(1:ncloc,1:nrloc,1,1) = 1.0
bed_fraction(1:ncloc,1:nrloc,1,2) = 0.0
bed_fraction(1:ncloc,1:nrloc,2,1) = 0.0
bed_fraction(1:ncloc,1:nrloc,2,2) = 1.0
bed_fraction(1:ncloc,1:nrloc,3,1) = 1.0
bed_fraction(1:ncloc,1:nrloc,3,2) = 0.0
bed_fraction(1:ncloc,1:nrloc,4:8,:) = 0.0
IF (nc1loc.EQ.1) bed_fraction(1,1:nrloc,:,:) = 0.0

CALL log_timer_out()


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
! Description - test case morphsort
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE gridpars
USE morpharrays
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_morphics'
CALL log_timer_in()

!---layer thickness
bed_layer_thickness(1:ncloc,1:nrloc,1) = 0.1
bed_layer_thickness(1:ncloc,1:nrloc,2) = 0.2
bed_layer_thickness(1:ncloc,1:nrloc,3) = 0.1
bed_layer_thickness(1:ncloc,1:nrloc,4:8) = 0.0
IF (nc1loc.EQ.1) bed_layer_thickness(1,1:nrloc,:) = 0.0

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_morphics

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
! Description - test case morphsort
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

dp(1) = 50.0E-06
dp(2) = 1000.0E-06
rhos(1:2) = 2650.0

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
