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
! *Usrdef_Time_Averages* Define time averaged output (example routines)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.1.2
!
! $Date: 2013-05-28 15:26:47 +0200 (Tue, 28 May 2013) $
!
! $Revision: 574 $
!
! Description -
!
! Reference -
!
! Routines - usrdef_avr_params, usrdef_avr0d_vals, usrdef_avr2d_vals,
!            usrdef_avr3d_vals
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_avr_params
!************************************************************************
!
! *usrdef_avr_params* Specifiers for time averaged output (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.1.2
!
! Description - if %ivarid is defined (non-zero value), the %f90_name,
!               %long_name, %units and %vector_name attributes do not need to
!               be defined, unless the user wants to use non-default values
!             - variables need to be defined in the order of increasing ranks:
!               0, 2, 3
!
! Reference -
!
! Calling program - time_averages_init
!
! External calls -
!
! Module calls -
!
!************************************************************************
!
USE iopars
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iset, istat, ivar, nostats


procname(pglev+1) = 'usrdef_avr_params'
CALL log_timer_in()

!
!1. Variable attributes
!----------------------
!

ivar_110: DO ivar=1,novarsavr
!  ---variable id
   avrvars(ivar)%ivarid  = 0
!  ---variable rank (0/2/3)
   avrvars(ivar)%nrank = ?
!  ---variable number
   avrvars(ivar)%numvar = -1
!  ---output operator (oopt_null,oopt_max,oopt_min,oopt_mean,oopt_klev,oopt_dep)
   avrvars(ivar)%oopt = oopt_null
!  ---vertical level for 2-D output of 3-D variables
   avrvars(ivar)%klev = 0
!  ---vertical depth for 2-D output of 3-D variables
   avrvars(ivar)%dep = 0.0
   IF (tsrvars(ivar)%ivarid.GT.0) THEN
!     ---fortran name
      avrvars(ivar)%f90_name = ?
!     ---long description name
      avrvars(ivar)%long_name = ?
!     ---variable unit
      avrvars(ivar)%units = ?
!     ---vector name (vector quantities only)
      avrvars(ivar)%vector_name = ?
   ENDIF
ENDDO ivar_110

!
!2. Variable indices
!-------------------
!

iset_210: DO iset=1,nosetsavr
   ivarsavr(iset,:) = ?
ENDDO iset_210

!
!3. File attributes
!-------------------
!

iset_310: DO iset=1,nosetsavr

!  ---defined
   avr0d(iset)%defined = .FALSE.; avr2d(iset)%defined = .FALSE.
   avr3d(iset)%defined = .FALSE.; avrgrd(iset)%defined = .FALSE.

!  ---file format ('A','U','N')
   avr0d(iset)%form = 'U'; avr2d(iset)%form = 'U'
   avr3d(iset)%form = 'U'; avrgrd(iset)%form = 'U'

!  ---file name
   avr0d(iset)%filename = ''; avr2d(iset)%filename = ''
   avr3d(iset)%filename = ''; avrgrd(iset)%filename = ''

!  ---info file
   avr0d(iset)%info = .TRUE.; avr2d(iset)%info = .TRUE.
   avr3d(iset)%info = .TRUE.; avrgrd(iset)%info = .TRUE.

ENDDO iset_310

!
!4. Output grid (space/time) attributes
!--------------------------------------
!

iset_410: DO iset=1,nosetsavr
!  ---regular or irregular (stations) grid
   avrgpars(iset)%gridded = .TRUE.
!  ---separate grid file
   avrgpars(iset)%grid_file = .FALSE.
!  ---apply land mask (regular grid)
   avrgpars(iset)%land_mask = .FALSE.
!  ---time-dependent spatial grid
   avrgpars(iset)%time_grid = .FALSE.
!  ---time resolution (start/end/period)
   avrgpars(iset)%tlims(1:3) = ?
!  ---date/time of first output (optional)
   avrgpars(iset)%startdate = cdatetime_undef
!  ---date/time of last output (optional)
   avrgpars(iset)%enddate = cdatetime_undef
!  ---reference date in case of a numerical time coordinate
   avrgpars(iset)%refdate = cdatetime_undef
!  ---grid dimension (0/2/3)
   avrgpars(iset)%nodim = ?
!  ---number of stations (irregular grid)
   avrgpars(iset)%nostats = ?
!  ---time format (0/1/2/3/4/5/6/7)
   avrgpars(iset)%time_format = 0
!  ---resolution in X-direction (regular grid)
   avrgpars(iset)%xlims(1:3) = ?
!  ---resolution in Y-direction (regular grid)
   avrgpars(iset)%ylims(1:3) = ?
!  ---resolution in  vertical direction 
   avrgpars(iset)%zlims(1:3) = ?

ENDDO iset_410

!
!5. Station attributes
!---------------------
!

istat_710: DO istat=1,nostatsavr
!  ---X-index
   avrstatlocs(istat)%ipos = ?
!  ---Y-index
   avrstatlocs(istat)%jpos = ?
!  ---station name
   avrstatlocs(istat)%name = ?
ENDDO istat_710

!
!6. Station labels
!-----------------
!

iset_610: DO iset=1,nosetsavr
   nostats = avrgpars(iset)%nostats
   lstatsavr(iset,1:nostats) = ?
ENDDO iset_610

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_avr_params

!========================================================================

SUBROUTINE usrdef_avr0d_vals(out0ddat,n0vars)
!************************************************************************
!
! *usrdef_avr0d_vals* 0-D output data (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.1.2
!
! Description - define output only for variables whose key id is set to zero
!
! Reference -
!
! Calling program - time_averages_update
!
! External calls -
!
! Module calls -
!
!************************************************************************
!

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: n0vars
REAL, INTENT(OUT), DIMENSION(n0vars) :: out0ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out0ddat* REAL    output data
!*n0vars*   INTEGER number of 0-D variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar


ivar_110: DO ivar=1,n0vars
   out0ddat(ivar) = ?
ENDDO ivar_110

      
RETURN

END SUBROUTINE usrdef_avr0d_vals

!========================================================================

SUBROUTINE usrdef_avr2d_vals(out2ddat,i,j,n2vars)
!************************************************************************
!
! *usrdef_avr2d_vals* 2-D output data (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.1.2
!
! Description - define output only for variables whose key id is set to zero 
!
! Reference -
!
! Calling program - time_averages_update
!
! External calls -
!
! Module calls -
!
!************************************************************************
!

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, n2vars
REAL, INTENT(OUT), DIMENSION(n2vars) :: out2ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out2ddat* REAL    output data
!*i*        INTEGER X-index of output location
!*j*        INTEGER Y-index of output location
!*n2vars*   INTEGER number of 2-D variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar


ivar_110: DO ivar=1,n2vars
   out2ddat(ivar) = ?
ENDDO ivar_110


RETURN

END SUBROUTINE usrdef_avr2d_vals

!========================================================================

SUBROUTINE usrdef_avr3d_vals(out3ddat,i,j,k,n3vars)
!************************************************************************
!
! *usrdef_avr3d_vals* 3-D output data (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.1.2
!
! Description - define output only for variables whose key id is set to zero 
!
! Reference -
!
! Calling program - time_averages_update
!
! External calls -
!
! Module calls -
!
!************************************************************************
!
  
IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, k, n3vars
REAL, INTENT(OUT), DIMENSION(n3vars) :: out3ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out3ddat* REAL    output data
!*i*        INTEGER X-index of output location
!*j*        INTEGER Y-index of output location
!*k*        INTEGER Vertical index of output location
!*n3vars*   INTEGER number of 3-D variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar


ivar_110: DO ivar=1,n3vars
   out3ddat(ivar) = ?
ENDDO ivar_110


RETURN

END SUBROUTINE usrdef_avr3d_vals
