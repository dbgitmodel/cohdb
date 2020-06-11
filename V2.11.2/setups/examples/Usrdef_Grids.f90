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
! *Usrdef_Grids* User-defined parallel grid generator (example routines)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Grids.f90  V2.0
!
! $Date: 2013-04-16 12:20:30 +0200 (Tue, 16 Apr 2013) $
!
! $Revision: 548 $
!
! Description -
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
! *usrdef_grid_params* Setup parameters for grid generator (example routine)
!
! Author - Patrick Luyten
!
!
! Description -
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
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_grid_params'
CALL log_timer_in()

!
!1. Process numbers
!------------------
!

npwork = ?; nprocsx = 0; nprocsy = 0

!
!2. Switches
!-----------
!
!---horizontal grid type (1/2/3)
iopt_grid_htype = 1
!---Cartesian or spherical (0/1)
iopt_grid_sph = ?
!---reset initial partition (0/1)
iopt_partit_reset = ?
!---postprocessing file (0/1)
iopt_partit_post = ?

!
!3. Model grid
!-------------
!
!---grid resolution
nc = ?; nr = ?; nz = ?
nobu = ?; nobv = ?

!---regular grid
IF (iopt_grid_htype.EQ.1) THEN
   surfacegrids(igrd_model,1)%delxdat = ?
   surfacegrids(igrd_model,1)%delydat = ?
   surfacegrids(igrd_model,1)%x0dat = ?
   surfacegrids(igrd_model,1)%y0dat = ?
ENDIF

!
!4. Title and output files
!-------------------------
!
!---run title
runtitle = ?
!---info file
info_file = ?
!---decomposition files
modfiles(io_mppmod,1)%form = 'A'
modfiles(io_mppmod,1)%filename = ' '
decomp_file = ?
!---mumap data file
post_file = ?

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_grid_params

!========================================================================

SUBROUTINE usrdef_model_grid
!************************************************************************
!
! *usrdef_model_grid* Define model grid (example routine)
!
! Author - Patrick Luyten
!
!
! Description -
!
! Reference -
!
! Calling program - generate_grids
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_model_grid'
CALL log_timer_in()

!
!1. Coordinate arrays
!--------------------
!
!---horizontal [m or degrees]
IF (iopt_grid_htype.EQ.2) THEN
   gxcoordglb(1:nc,1) = ?
   gxcoordglb(1,1:nr) = ?
ELSEIF (iopt_grid_htype.EQ.3) THEN
   gxcoordglb(1:nc,1:nr) = ?
   gycoordglb(1:nc,1:nr) = ?
ENDIF

!
!2. Water depths [m]
!-------------------
!

depmeanglb(1:nc-1,1:nr-1) = ?

!
!3. Open boundary locations
!--------------------------
!
!---U-nodes
IF (nobu.GT.0) THEN
   iobu(1:nobu) = ?
   jobu(1:nobu) = ?
ENDIF

!---V-nodes
IF (nobv.GT.0) THEN
   iobv(1:nobv) = ?
   jobv(1:nobv) = ?
ENDIF

CALL log_timer_out()


RETURN

