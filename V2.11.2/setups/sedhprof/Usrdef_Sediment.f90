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
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.10.3
!
! $Date: 2016-10-26 15:58:08 +0200 (Wed, 26 Oct 2016) $
!
! $Revision: 983 $
!
! Description - test case sedhprof
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
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.10.3
!
! Description - test case sedhprof
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
USE switches
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
!---transport mode (1/2/3/4)
iopt_sed_mode = 2
!---grid dimension (2/3)
iopt_sed_nodim = MERGE(2,3,iopt_grid_nodim.EQ.2)
!---type of sediment (1/2)
iopt_sed_type = MERGE(1,2,LLT(runtitle(9:9),'G'))

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
!---type of settling velocity (0/1/2/3/4/5)
iopt_sed_ws = 2

!
!1.5 Bottom boundary condition
!-----------------------------
!
!---formulation for reference concentration (1/2/3/4/5)
SELECT CASE (runtitle(9:9))
   CASE ('A')
      iopt_sed_bbc_ref = 1
      iopt_sed_bbc_eq = 1
   CASE ('B')
      iopt_sed_bbc_ref = 2
      iopt_sed_bbc_eq = 3
   CASE ('C')
      iopt_sed_bbc_ref = 1
      iopt_sed_bbc_eq = 4
   CASE ('D')
      iopt_sed_bbc_ref = 2
      iopt_sed_bbc_eq = 5
   CASE ('E')
      iopt_sed_bbc_ref = 1
      iopt_sed_bbc_eq = 6
   CASE ('F')
      iopt_sed_bbc_ref = 1
   CASE ('G')
      iopt_sed_bbc_ref = 1
    CASE ('H')
      iopt_sed_bbc_ref = 2
    CASE ('I')
      iopt_sed_bbc_ref = 2
  END SELECT

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_sed_params

!========================================================================

SUBROUTINE usrdef_morph_params
!************************************************************************
!
! *usrdef_morph_params* Define parameters for the morphological model
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


RETURN

END SUBROUTINE usrdef_morph_params

!========================================================================

SUBROUTINE usrdef_dar_params
!************************************************************************
!
! *usrdef_dar_params* Define parameters for dredging and relocation
!
! Author - Alexander Breugem and  Kevin Delecluyse (IMDC)
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

!=======================================================================

SUBROUTINE usrdef_morphics
!************************************************************************
!
! *usrdef_morphics* Define initial conditions for the morphological model
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


RETURN

END SUBROUTINE usrdef_morphics

!========================================================================

SUBROUTINE usrdef_sed_spec
!************************************************************************
!
! *usrdef_sed_spec* Define characteristics of sediment fractions
!
! Author - Alexander Breugem
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.5
!
! Description - test case sedhprof
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

!---diameter [m]
dp = 89.4777E-06

!---density [m^3/s]
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


RETURN

END SUBROUTINE usrdef_dar_spec
