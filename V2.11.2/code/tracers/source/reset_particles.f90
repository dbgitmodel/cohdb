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

MODULE reset_particles
!************************************************************************
!
! *reset_particles* Reset model parameters and arrays for the particle model
!                   if needed
!
! Author - Valerie Duliere and Patrick Luyten
!
! Revision - @(COHERENS)reset_particles.f90  V2.11.2
!
! $Date: 2018-07-16 17:34:12 +0200 (Mon, 16 Jul 2018) $
!
! $Revision: 1164 $
!
! Description - 
!
! Reference -
!
! Routines - reset_out_files_part, reset_out_ppars, reset_part_params,
!            reset_partics
!
!************************************************************************
!
USE iopars


IMPLICIT NONE

CONTAINS

!========================================================================

SUBROUTINE reset_out_files_part(filepars)
!************************************************************************
!
! *reset_out_files_part* Reset attributes of a user particle output file
!
! Author - Valerie Duliere and Patrick Luyten
!
! Revision - @(COHERENS)reset_particles.f90  V2.11
!
! Description - 
!
! Module calls -
!
!************************************************************************
!
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
TYPE (FileParams), INTENT(INOUT), DIMENSION(:) :: filepars  
!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*filepars*  DERIVED Attributes of user-output files
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize


nsize = SIZE(filepars)
IF (nsize.EQ.0) RETURN

procname(pglev+1) = 'reset_out_files_part'
CALL log_timer_in()

!---reset attributes
n_110: DO n=1,nsize
   IF (filepars(n)%defined) THEN
      filepars(n)%status = 'W'
      filepars(n)%iunit = -1
   ENDIF
ENDDO n_110

CALL log_timer_out()


RETURN

END SUBROUTINE reset_out_files_part

!========================================================================

SUBROUTINE reset_out_ppars(end_reset)
!************************************************************************
!
! *reset_out_ppars* Reset parameters for particle output
!
! Author - Valerie Duliere and Patrick Luyten
!
! Revision - @(COHERENS)reset_particles.f90  V2.11
!
! Description - 
!
! Module calls - add_secs_to_date, warning_reset_arr_struc
!
!************************************************************************
!
USE partpars
USE partvars  
USE timepars
USE syspars
USE error_routines, ONLY: warning_reset_arr_struc
USE time_routines, ONLY: add_secs_to_date, log_timer_in, log_timer_out

!
!*Arguments
!
LOGICAL, INTENT(IN) :: end_reset
  
!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*end_reset* LOGICAL Resets end value of time resolution (%tlims(2)) to a
!                    value lower than or equal to nstep if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: iset, limend
INTEGER, DIMENSION(3) :: tlims(3)
REAL :: deltout


procname(pglev+1) = 'reset_out_ppars'
CALL log_timer_in()

!
!1. Set output dates
!-------------------
!

iset_110: DO iset=1,nosetspart

!  ---reset end limit (if necessary)
   tlims = outppars(iset)%tlims
   IF (end_reset.AND.tlims(2).GT.nstep) THEN
      limend = tlims(1) + ((nstep-tlims(1))/tlims(3))*tlims(3)
      CALL warning_reset_arr_struc(tlims(2),'outppars','%tlims(2)',limend,1,&
                                & (/iset/))
      outppars(iset)%tlims = tlims
   ENDIF
   tlims = outppars(iset)%tlims

!  ---start date
   CALL add_secs_to_date(CStartDateTime,outppars(iset)%startdate,tlims(1),&
                       & delt2d)

!  ---end date
   CALL add_secs_to_date(CStartDateTime,outppars(iset)%enddate,tlims(2),delt2d)

!  ---reference date
   IF (outppars(iset)%refdate.EQ.cdatetime_undef) THEN
       outppars(iset)%refdate = RefDateTime_part
   ENDIF

ENDDO iset_110

!
!2. Grid and time dimensions
!---------------------------
!

iset_210: DO iset=1,nosetspart
   tlims = outppars(iset)%tlims
   outppars(iset)%nstepout = (tlims(2)-tlims(1))/tlims(3) + 1
   deltout = tlims(3)*delt2d
   outppars(iset)%deltout = deltout/time_convert(ptime_unit)
ENDDO iset_210

CALL log_timer_out()


RETURN

END SUBROUTINE reset_out_ppars

!========================================================================

SUBROUTINE reset_part_params
!************************************************************************
!
! *reset_part_params* Reset particle model parameters
!
! Author - Valerie Duliere and Patrick Luyten
!
! Revision - @(COHERENS)reset_particles.f90  V2.11.2
!
! Description - 
!
! Module calls - error_lbound_var, error_lbound_var_date, error_ubound_var_date,
!                num_time_steps, warning_reset_var
!
!************************************************************************
!
USE gridpars  
USE iopars
USE paralpars
USE partpars
USE partswitches
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_abort, error_lbound_var, warning_reset_var
USE time_routines, ONLY: error_lbound_var_date, error_ubound_var_date,&
                         & log_timer_in, log_timer_out, num_time_steps

!
!*Local variables
!
CHARACTER (LEN=15) :: cval
CHARACTER (LEN=22) :: cstep
INTEGER :: n


procname(pglev+1) = 'reset_part_params'
CALL log_timer_in()

!
!1. Reset switches if necessary
!------------------------------
!
!---1-D case
IF (iopt_grid_nodim.EQ.1) THEN
   CALL warning_reset_var(iopt_part_hadv,'iopt_part_hadv',0)
   CALL warning_reset_var(iopt_part_vadv,'iopt_part_vadv',0)
   CALL warning_reset_var(iopt_part_hdif,'iopt_part_hdif',0)
   CALL warning_reset_var(iopt_part_cloud,'iopt_part_cloud',0)
   CALL warning_reset_var(iopt_part_leeway,'iopt_part_leeway',0)
ENDIF

!---2-D case
IF (iopt_grid_nodim.EQ.2) THEN
   CALL warning_reset_var(iopt_part_vadv,'iopt_part_vadv',0)
   CALL warning_reset_var(iopt_part_vdif,'iopt_part_vdif',0)
   CALL warning_reset_var(iopt_part_wind,'iopt_part_wind',0)
ENDIF

!---meteo
IF (iopt_meteo.EQ.0) THEN
   CALL warning_reset_var(iopt_part_wind,'iopt_part_wind',0)
   CALL warning_reset_var(iopt_part_leeway,'iopt_part_leeway',0)
ELSEIF (iopt_part_wind.EQ.1) THEN
   IF ((iopt_part_model.EQ.2.AND.modelid.EQ.modelidpart).OR.&
      & iopt_part_model.GT.2) THEN
      CALL warning_reset_var(iopt_meteo_data,'iopt_meteo_data',1)
      CALL warning_reset_var(iopt_meteo_stres,'iopt_meteo_stres',1)
      CALL warning_reset_var(iopt_meteo_heat,'iopt_meteo_heat',0)
      CALL warning_reset_var(iopt_meteo_pres,'iopt_meteo_pres',0)
      CALL warning_reset_var(iopt_meteo_precip,'iopt_meteo_precip',0)
      CALL warning_reset_var(iopt_sflux_pars,'iopt_sflux_pars',0)
   ENDIF
ENDIF

!---uniform density
IF (iopt_dens.EQ.0.AND.iopt_part_dens.EQ.2) THEN
   CALL warning_reset_var(iopt_part_dens,'iopt_part_dens',1)
ENDIF
   
!
!2. Parameters
!-------------
!
!---grid dimensions
IF (iopt_part_model.GT.1) THEN
   ncloc = nc; nrloc = nr
ENDIF

!---particle clouds
IF (iopt_part_cloud.EQ.0) THEN
   CALL warning_reset_var(noclouds,'noclouds',1)
ENDIF

!---backward/forward drifting
back_drift = iopt_part_model.EQ.4
for_drift = iopt_part_model.LT.4

!---time step
deltpart = MERGE(delt3d,-delt3d,for_drift)

!---date/time
IF (iopt_part_model.EQ.3) THEN
   CALL error_lbound_var_date(CEndDateTime,'CEndDateTime',CStartDateTime,&
                           & .FALSE.)
   CALL num_time_steps(CStartDateTime,CEndDateTime,deltpart,nstep)
ELSEIF (iopt_part_model.EQ.4) THEN
   CALL error_ubound_var_date(CEndDateTime,'CEndDateTime',CStartDateTime,&
                           & .FALSE.)
   CALL num_time_steps(CEndDateTime,CStartDateTime,-deltpart,nstep)
ENDIF
CALL error_lbound_var(nstep,'nstep',0,.FALSE.)

!---restart times
n_210: DO n=1,norestarts
   IF (ntrestart(n).EQ.int_fill) ntrestart(n) = nstep
ENDDO n_210

!---write log info
IF (loglev1.GT.0) THEN
   WRITE (cstep,'(I22)') nstep; cstep = ADJUSTL(cstep)
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'nstep = '//TRIM(cstep)
   WRITE (cval,'(G15.7)') deltpart; cval = ADJUSTL(cval)
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'deltpart = '//TRIM(cval)
ENDIF

!---log file counter
IF (iopt_part_model.GT.2.AND.runlog_count.EQ.int_fill) runlog_count = nstep

CALL error_abort('reset_part_params',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE reset_part_params

!=====================================================================

SUBROUTINE reset_partics
!************************************************************************
!
! *reset_partics* Reset initial conditions for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)reset_particles.f90  V2.9
!
! Description - 
!
! Calling program - initialise_model
!
! Module calls -
!  
!************************************************************************
!
USE iopars
USE partvars  
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'reset_partics'
CALL log_timer_in()

IF (iopt_grid_nodim.EQ.2) THEN
   WHERE (part%state.EQ.2)
      part%state = 1
   END WHERE
ELSEIF (iopt_grid_nodim.EQ.1) THEN
   WHERE (part%state.EQ.1)
      part%state = 2
   END WHERE
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE reset_partics


END MODULE reset_particles
