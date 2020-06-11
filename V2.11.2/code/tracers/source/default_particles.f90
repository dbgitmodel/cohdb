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

MODULE default_particles
!************************************************************************
!
! *default_particles* Default settings for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)default_particles.f90  V2.11
!
! $Date: 2017-04-12 16:52:46 +0200 (Wed, 12 Apr 2017) $
!
! $Revision: 1012 $
!
! Description - 
!
! Reference -
!
! Routines - default_out_ppars, default_part_params
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS

!====================================================================

SUBROUTINE default_out_ppars
!************************************************************************
!
! *default_out_ppars* Default values for trajectory output parameters
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)default_particles.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE partpars
USE partvars
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
INTEGER :: iset

procname(pglev+1) = 'default_out_ppars'
CALL log_timer_in()

iset_110: DO iset=1,nosetspart
   outppars(iset)%aging = .FALSE.
   outppars(iset)%ktype = 0
   outppars(iset)%label = 1
   outppars(iset)%nodim = 0
   outppars(iset)%ntype = 0
   outppars(iset)%refdate = cdatetime_undef
ENDDO iset_110

CALL log_timer_out()


RETURN

END SUBROUTINE default_out_ppars

!====================================================================

SUBROUTINE default_part_params
!************************************************************************
!
! *default_part_params* Default settings of parameters for the particle  model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)default_particles.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE partpars
USE partswitches
USE switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'default_part_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!
!---advection
iopt_part_hadv = MERGE(0,2,iopt_grid_nodim.EQ.1)
iopt_part_vadv = MERGE(0,1,iopt_grid_nodim.EQ.2)

!---sign change in leeway drift
iopt_part_leeway = 0

!---buoyancy
iopt_part_dens = 0

!---diffusion
iopt_part_hdif = 0; iopt_part_vdif = 0

!---particle cloud generator
iopt_part_cloud = 0

!---use wind current to drive particle drift
iopt_part_wind = 1

!---volume concentrations
iopt_part_conc = 0

!---output from particle module
iopt_part_out = 1

!
!2. Parameters
!-------------
!
!---reference date/time
IF (iopt_part_model.LT.4) THEN
   RefDateTime_part = CStartDateTime
ELSE
   RefDateTime_part = CEndDateTime
ENDIF

!---number of particles
nopart = 1

!---number of clouds
noclouds = 1

!---number of labels
nolabels = 1

!---time unit
ptime_unit = 1

!---counter for physcial update
icpart = 1

!---number of output trajectories
nosetspart = 1

!---current drift factor
cdrift_fac = 1.0

!---factors for wind drift
wdrift_angle1 =  40.0; wdrift_angle2 = 8.0
wdrift_slope = 0.0315

!---(uniform) diffusion coefficients
xdifpart_cst = 0.0;  ydifpart_cst = 0.0; zdifpart_cst = 0.0

CALL log_timer_out()


RETURN

END SUBROUTINE default_part_params

END MODULE default_particles
