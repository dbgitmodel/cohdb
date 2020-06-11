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
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.10.3
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description - test case seddens
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
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.10.3
!
! Description - test case seddens
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
use sedpars
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
iopt_sed_mode = 2
!---type of sediment
iopt_sed_type = 1

!
!1.3 Bottom stress
!-----------------
!
!---formulation for critical shear stress
iopt_sed_bstres_cr = 1

!
!1.4 Settling
!------------
!
!---type of settling velocity
SELECT CASE (runtitle(8:8))
   CASE ('A','B')
      iopt_sed_ws = 1
   CASE ('C')
      iopt_sed_ws = 2
   CASE ('D')
      iopt_sed_ws = 4
   CASE ('E')
      iopt_sed_ws = 5
END SELECT

!
!1.7 Density effects
!-------------------
!
!---disable/enable density effects
iopt_sed_dens_grad = 1

cgel = 0.1

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_sed_params

!=========================================================

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
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.5
!
! Description - test case seddens
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

!
!*Local variables
!
INTEGER :: i, iglb


procname(pglev+1) = 'usrdef_sedics'
CALL log_timer_in()

i_110: DO i=1,ncloc
   iglb = nc1loc-1+i
   cvol(i,1,:,1) = MERGE(0.05,0.0,iglb.LE.5)
ENDDO i_110

CALL log_timer_out()


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
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.5
!
! Description - test case seddens
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
dp = 250.0E-06

!---density [m^3/s]
rhos = 2650.0

!---(uniform) critical stress [m^2/s^2]
!ensure that there is no sediment entrainment
bstres_cr_cst = 10.0  

!---(uniform) settling velocity [m/s]
SELECT CASE (runtitle(8:8))
   CASE ('A')
      ws_cst = 1.0E-08
   CASE ('B')
!Hogg et al. (1999): dimensionless settling velocity beta = 5E-03:
!Reference velocity Ur = sqrt(g0')=0.8093 m/s => ws=beta*Ur
      ws_cst = 4.498E-03
END SELECT

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
