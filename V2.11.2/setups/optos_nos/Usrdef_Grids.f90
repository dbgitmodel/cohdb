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
! *Usrdef_Grids* User-defined parallel grid generator
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Grids.f90  V2.0
!
! $Date: 2015-10-06 17:42:25 +0200 (Tue, 06 Oct 2015) $
!
! $Revision: 886 $
!
! Description - optos nos
!
! Reference -
!
! Routines - usrdef_grid_params, usrdef_model_grid
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_grid_params
!************************************************************************
!
! *usrdef_grid_params* Setup parameters for grid generator
!
! Author - Patrick Luyten
!
! Description - optos nos
!
! Reference -
!
! Calling program - generate_grids
!
!************************************************************************
!
USE gridpars
USE grid_params
USE iopars
USE paralpars
USE physpars
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
REAL :: delxdat, delydat


procname(pglev+1) = 'usrdef_grid_params'
CALL log_timer_in()

!
!1. Process numbers
!------------------
!

npwork = 12

!
!2. Switches
!-----------
!
!---Cartesian or spherical (0/1)
iopt_grid_sph = 1
!---postprocessing file (0/1)
iopt_partit_post = 1

!
!3. Model grid
!-------------
!
!---grid resolution
nc = 158; nr = 206; nz = 20
nobu = 36; nobv = 127

!---number of grid levels
nomglevels = 3

!---regular grid
delxdat = 1/12.0; delydat =  1/24.0
surfacegrids(igrd_model,1)%delxdat = delxdat
surfacegrids(igrd_model,1)%delydat = delydat
surfacegrids(igrd_model,1)%x0dat = -4.0 - 0.5*delxdat
surfacegrids(igrd_model,1)%y0dat = 48.5 - 0.5*delydat

!
!4. Title and output files
!-------------------------
!
!---run title
runtitle = 'optnos'

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_grid_params

!========================================================================

SUBROUTINE usrdef_model_grid
!************************************************************************
!
! *usrdef_model_grid* Define model grid
!
! Author - Patrick Luyten
!
! Description - optos nos
!
! Reference -
!
! Calling program - generate_grids
!
! Module calls - close_file, open_file
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: ii, iunit, j, jj

procname(pglev+1) = 'usrdef_model_grid'
CALL log_timer_in()


!
!1. Open grid file
!-----------------
!

CALL open_file(iunit,'nos_grid.dat','IN','A')
READ (iunit,'(/)')

!
!2. Water depths [m]
!-------------------
!

j_210: DO j=1,nr-1
   READ (iunit,*) depmeanatc(1:nc-1,j)
ENDDO j_210

!
!3. Open boundary locations
!--------------------------
!
!---U-nodes
READ (iunit,*)
ii_310: DO ii=1,nobu
   READ (iunit,*) iobu(ii), jobu(ii)
ENDDO ii_310

!---V-nodes
READ (iunit,*)
jj_320: DO jj=1,nobv
   READ (iunit,*) iobv(jj), jobv(jj)
ENDDO jj_320

CALL close_file(iunit,'asc')

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_model_grid
