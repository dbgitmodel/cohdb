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

MODULE check_particles
!************************************************************************
!
! *check_particles* Check parameters for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)check_particles.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
! Reference -
!
! Routines - check_out_ppars, check_partics, check_part_params
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS

!=====================================================================

SUBROUTINE check_out_ppars
!************************************************************************
!
! *check_out_ppars* check parameters for trajectory output
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)check_particles.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - particle_trajects_init
!
! Module calls - check_time_limits_arr_struc, error_abort,
!                error_limits_arr_struc, error_vals_arr_struc 
!
!************************************************************************
!
USE iopars
USE partpars
USE partvars
USE timepars
USE error_routines, ONLY: check_time_limits_arr_struc, error_abort, &
                        & error_limits_arr_struc, error_vals_arr_struc_int
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!* Local variables
!
INTEGER :: iset
INTEGER, DIMENSION(3) :: tlims


procname(pglev+1) = 'check_out_ppars'
CALL log_timer_in()


iset_110: DO iset=1,nosetspart
   
!  ---output times   
   tlims = outppars(iset)%tlims
   CALL check_time_limits_arr_struc(tlims,'outppars','%tlims',0,nstep,1,&
                 & (/iset/))

!  ---type of vertical output
   CALL error_limits_arr_struc(outppars(iset)%ktype,'outppars','ktype',0,2,&
                             & 1,(/iset/))

!  ---labels
   CALL error_limits_arr_struc(outppars(iset)%label,'outppars','label',1,&
                             & nolabels,1,(/iset/))

!  ---dimension
   CALL error_vals_arr_struc_int(outppars(iset)%nodim,'outppars','nodim',3,1,&
                              & (/0,2,3/),(/iset/))
   
!  ---type of particle output
   CALL error_limits_arr_struc(outppars(iset)%ntype,'outppars','ntype',0,3,1,&
                            & (/iset/))
   
ENDDO iset_110

CALL error_abort('check_out_ppars',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_out_ppars

!=====================================================================

SUBROUTINE check_partics
!************************************************************************
!
! *check_partics* Check initial conditions for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)check_particles.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
! Module calls - error_abort, error_lbound_arr, error_lbounds_arr_struc,
!                error_limits_arr_struc, error_value_arr_struc
!
!************************************************************************
!
USE iopars
USE partpars
USE partswitches
USE partvars
USE syspars
USE switches
USE error_routines, ONLY: error_abort, error_lbound_arr, &
                        & error_lbound_arr_struc, error_limits_arr_struc, &
                        & error_value_arr_struc
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
CHARACTER (LEN=12) :: cfix, cval 
INTEGER :: icld, l, p, nopcloud


  
procname(pglev+1) = 'check_partics'
CALL log_timer_in()

!
!1. Particulate matter
!---------------------
!

l_110: DO l=1,nolabels
   IF (iopt_part_conc.GE.2) THEN
      CALL error_lbound_arr(diamconc(l),'diamconc',0.0,.FALSE.,1,indx=(/l/))
   ENDIF
   IF (iopt_part_conc.EQ.3) THEN 
      CALL error_lbound_arr(rhoconc(l),'rhoconc',0.0,.FALSE.,1,indx=(/l/))
   ENDIF
ENDDO l_110

!
!2. Particle clouds
!------------------
!

IF (iopt_part_cloud.GT.0) THEN

!  ---number of particles
   nopcloud = 0
   icld_210: DO icld=1,noclouds
      nopcloud = nopcloud + cloud(icld)%nopart*cloud(icld)%noreleases
   ENDDO icld_210
   IF (nopcloud.NE.nopart) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cval,'(I12)') nopcloud; WRITE(cfix,'(I12)') nopart
         cval = ADJUSTL(cval); cfix = ADJUSTL(cfix)
         WRITE (ioerr,'(A)') 'Invalid value for the total number of cloud'//&
                            & ' particles: '//TRIM(cval)
         WRITE (ioerr,'(A)') 'Must be equal to the value of nopart: '//&
                           & TRIM(cfix)
      ENDIF
   ENDIF

 ! ---cloud attributes  
   icld_220: DO icld=1,noclouds
      CALL error_limits_arr_struc(cloud(icld)%kdistype,'cloud','kdistype',0,3,&
                                & 1,(/icld/))
      CALL error_limits_arr_struc(cloud(icld)%label,'cloud','label',1,nolabels,&
                                & 1,(/icld/))
      CALL error_lbound_arr_struc(cloud(icld)%noreleases,'cloud','noreleases',&
                                & 0,.FALSE.,1,(/icld/))
      CALL error_limits_arr_struc(cloud(icld)%state,'cloud','state',1,2,1,&
                               & (/icld/))
      CALL error_lbound_arr_struc(cloud(icld)%length,'cloud','length',0.0,&
                               & .FALSE.,1,(/icld/))
      CALL error_limits_arr_struc(cloud(icld)%orientation,'cloud',&
           & 'orientation',0.0,360.0,1,(/icld/))
      IF (cloud(icld)%state.EQ.2) THEN
         CALL error_lbound_arr_struc(cloud(icld)%thick,'cloud','thick',0.0,&
                                  & .FALSE.,1,(/icld/))
      ENDIF
      CALL error_lbound_arr_struc(cloud(icld)%width,'cloud','width',0.0,&
                               & .FALSE.,1,(/icld/))
   ENDDO icld_220

ENDIF

!
!3. Particle initial conditions
!------------------------------
!

IF (modfiles(io_inicon,ics_part,1)%defined) THEN

   p_310: DO p=1,nopart
      CALL error_limits_arr_struc(part(p)%drift_state,'part','drift_state',1,5,&
                                & 1,(/p/))
      IF (iopt_grid_nodim.EQ.2) THEN
         CALL error_value_arr_struc(part(p)%state,'part','state',1,1,(/p/))
      ELSE
         CALL error_limits_arr_struc(part(p)%state,'part','state',1,2,1,(/p/))
      ENDIF
      CALL error_limits_arr_struc(part(p)%label,'part','label',1,nolabels,1,&
                               & (/p/))
      IF (iopt_grid_sph.EQ.1) THEN
         CALL error_limits_arr_struc(part(p)%xpos,'part','xpos',&
                                   & -180.0_kndrlong,180.0_kndrlong,1,(/p/))
         CALL error_limits_arr_struc(part(p)%ypos,'part','ypos',&
                                   & -90.0_kndrlong,90.0_kndrlong,1,(/p/))
      ENDIF

   ENDDO p_310

ENDIF

CALL error_abort('check_partics',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_partics

!=====================================================================

SUBROUTINE check_part_params
!************************************************************************
!
! *check_part_params* Check parameters for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)check_particles.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
! Module calls - add_secs_to_date, check_date, error_abort, error_lbound_var,
!                error_limits_var, error_vals_arr_struc_char, error_value_var,
!                mult_index, warning_reset_arr_struc
!
!************************************************************************
!
USE iopars
USE partpars
USE partswitches
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_abort, error_lbound_var, &
                        & error_limits_var, error_vals_arr_struc_char, &
                        & error_value_var, warning_reset_arr_struc
USE time_routines, ONLY: add_secs_to_date, check_date, log_timer_in, &
                       & log_timer_out
USE utility_routines, ONLY: mult_index

!
!*Local variables
!
CHARACTER (LEN=15) :: cval
CHARACTER (LEN=lentime) :: cdatetimex


procname(pglev+1) = 'check_part_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!
!1.1 Ranges
!----------
!
!---particle cloud generator
CALL error_limits_var(iopt_part_cloud,'iopt_part_cloud',0,1)

!---volume concentrations
CALL error_limits_var(iopt_part_conc,'iopt_part_conc',0,3)

!---buoyancy
CALL error_limits_var(iopt_part_dens,'iopt_part_dens',0,2)

!---horizontal advection
CALL error_limits_var(iopt_part_hadv,'iopt_part_hadv',0,2)

!---horizontal diffusion
CALL error_limits_var(iopt_part_hdif,'iopt_part_hdif',0,1)

!---trajectory output
CALL error_limits_var(iopt_part_out,'iopt_part_out',0,1)

!---leeway drift
CALL error_limits_var(iopt_part_leeway,'iopt_part_leeway',0,1)

!---vertical advection
CALL error_limits_var(iopt_part_vadv,'iopt_part_vadv',0,1)

!---vertical diffusion
CALL error_limits_var(iopt_part_vdif,'iopt_part_vdif',0,2)

!---wind drift
CALL error_limits_var(iopt_part_wind,'iopt_part_wind',0,1)

!
!1.2 Incompatibilities
!---------------------
!
!---meteo
IF (iopt_part_wind.EQ.1) THEN
   IF (iopt_meteo.NE.1) THEN
      nerrs = nerrs + 1
      IF (errchk) THEN
         WRITE (cval,'(I12)') iopt_meteo; cval = ADJUSTL(cval)
         WRITE (ioerr,'(A)') 'Invalid value for switch iopt_meteo: '//TRIM(cval)
         WRITE (ioerr,'(A)') &
              & 'Must be set to 1 to enable wind input in the particle model'
      ENDIF
   ENDIF
   IF (iopt_meteo_data.NE.1) THEN 
      nerrs = nerrs + 1
      IF (errchk) THEN
         WRITE (cval,'(I12)') iopt_meteo_data; cval = ADJUSTL(cval)
         WRITE (ioerr,'(A)') 'Invalid value for switch iopt_meteo_data: '//&
                           & TRIM(cval)
         WRITE (ioerr,'(A)') &
                & 'Must be set to 1 to enable wind input in the particle model'
      ENDIF
   ENDIF
   IF (iopt_meteo_stres.NE.1) THEN 
      nerrs = nerrs + 1
      IF (errchk) THEN
         WRITE (cval,'(I12)') iopt_meteo_stres; cval = ADJUSTL(cval)
         WRITE (ioerr,'(A)') 'Invalid value for switch iopt_meteo_stres: '//&
                           & TRIM(cval)
         WRITE (ioerr,'(A)') &
                & 'Must be set to 1 to enable wind input in the particle model'
      ENDIF
   ENDIF
ENDIF

!---density
IF (iopt_grid_nodim.EQ.2) THEN
   IF (iopt_part_dens.EQ.2) THEN
      nerrs = nerrs + 1
      IF (errchk) THEN
         WRITE (cval,'(I12)') iopt_part_dens; cval = ADJUSTL(cval)
         WRITE (ioerr,'(A)') 'Invalid value for switch iopt_part_dens: '//&
                           & TRIM(cval)
         WRITE (ioerr,'(A)') &
              & 'Allowed values are 0 and 1 in case of 2-D runs'
      ENDIF
   ENDIF
ENDIF

!
!2. Parameters
!-------------
!
!---integer number of external time steps within simulation period
IF (iopt_part_model.EQ.3) THEN
   CALL add_secs_to_date(CStartDateTime,cdatetimex,nstep,deltpart)
   IF (cdatetimex.NE.CEndDateTime) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cval,'(G15.7)') delt2d
         cval = ADJUSTL(cval)
         WRITE (ioerr,'(A)') 'Invalid external time step: '//TRIM(cval)
         WRITE (ioerr,'(A)') 'Simulation period must contain an integer '//&
                          &  'number of external time steps'
      ENDIF
   ENDIF
ELSEIF (iopt_part_model.EQ.4) THEN
   CALL add_secs_to_date(CEndDateTime,cdatetimex,nstep,-deltpart)
   IF (cdatetimex.NE.CStartDateTime) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cval,'(G15.7)') delt2d
         cval = ADJUSTL(cval)
         WRITE (ioerr,'(A)') 'Invalid external time step: '//TRIM(cval)
         WRITE (ioerr,'(A)') 'Simulation period must contain an integer '//&
                          &  'number of external time steps'
      ENDIF
   ENDIF
ENDIF

!---time counter for communication between COHERENS and particle model
IF (iopt_part_model.GE.3.AND.(.NOT.mult_index(icpart,ic3d))) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') delt2d*icpart
      cval = ADJUSTL(cval)
      WRITE (ioerr,'(A)') 'Invalid coupling time step between COHERENS and '//&
                        & ' the particle model: '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be a multiple of the 3-D time step'
   ENDIF
ENDIF

!---labels
CALL error_lbound_var(nolabels,'nolabels',1,.TRUE.)

!---reference date/time
CALL check_date(RefDateTime_part,'RefDateTime_part')

!---time unit
CALL error_limits_var(ptime_unit,'ptime_unit',1,4)

!---diffusion coefficients
CALL error_lbound_var(xdifpart_cst,'xdifpart_cst',0.0,.TRUE.)
CALL error_lbound_var(ydifpart_cst,'ydifpart_cst',0.0,.TRUE.)
CALL error_lbound_var(zdifpart_cst,'zdifpart_cst',0.0,.TRUE.)

!---particle clouds
IF (iopt_part_cloud.EQ.1) THEN
   CALL error_lbound_var(noclouds,'noclouds',1,.TRUE.)
ENDIF

!---counter for log file
IF (iopt_part_model.GT.2.AND.loglev1.GT.0) THEN
   CALL error_limits_var(runlog_count,'runlog_count',1,nstep)
ENDIF

!
!3. Attributes of forcing files
!------------------------------
!

IF (iopt_part_write.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_pargrd,1,2)%status,&
                           & 'modfiles','%status','"W"',3,(/io_pargrd,1,2/))
   CALL error_vals_arr_struc_char(modfiles(io_parphs,1,2)%status,&
                           & 'modfiles','%status','"W"',3,(/io_parphs,1,2/))
ELSE
   CALL warning_reset_arr_struc(modfiles(io_pargrd,1,2)%status,&
                             & 'modfiles','%status','0',3,(/io_pargrd,1,2/))
   CALL warning_reset_arr_struc(modfiles(io_parphs,1,2)%status,&
                             & 'modfiles','%status','0',3,(/io_parphs,1,2/))
ENDIF
IF (iopt_part_model.GT.2) THEN
   CALL error_vals_arr_struc_char(modfiles(io_pargrd,1,1)%status,&
                           & 'modfiles','%status','"R" "N"',3,(/io_pargrd,1,1/))
   CALL error_vals_arr_struc_char(modfiles(io_parphs,1,1)%status,&
        & 'modfiles','%status','"R" "N"',3,(/io_parphs,1,1/))
ENDIF

CALL error_abort('check_part_params',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_part_params

END MODULE check_particles
