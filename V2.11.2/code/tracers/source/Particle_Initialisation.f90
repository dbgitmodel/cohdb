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
! *Particle_Initialisation* Series of routines for initialisation of the
!                           particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Initialisation.f90  V2.11
!
! $Date: 2018-06-14 15:25:04 +0200 (Thu, 14 Jun 2018) $
!
! $Revision: 1155 $
!
! Description - 
!
! Reference -
!
! Routines - define_partics, initialise_particle_model, particle_releases,
!            read_partics, read_part_spc
!
!************************************************************************
!

!========================================================================

SUBROUTINE define_partics
!************************************************************************
!
! *define_partics* Define initial conditions for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Model_Initialisation.F90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
! External calls - initialise_particle_arrays, particle_cloud, read_partics,
!                  rng_open, usrdef_partics
!
! Module calls - rng_open
!
!************************************************************************
!
USE iopars
USE partpars
USE partswitches
USE rng_library, ONLY: rng_open
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'define_partics'
CALL log_timer_in()

!---using pre- or user-defined locations   
IF (iopt_part_cloud.EQ.0) THEN
   IF (modfiles(io_inicon,ics_part,1)%status.EQ.'R') THEN
      CALL read_partics
   ELSEIF (modfiles(io_inicon,ics_part,1)%status.EQ.'N') THEN
      CALL usrdef_partics
   ENDIF
   CALL initialise_particle_arrays
!---from cloud releases   
ELSE
   CALL particle_cloud
ENDIF

!---open random generators
IF (iopt_part_vdif.GT.0.OR.iopt_part_hdif.GT.0) CALL rng_open(prangenwalk)
IF (iopt_part_leeway.EQ.1) CALL rng_open(prangenlw)

CALL log_timer_out()


RETURN

END SUBROUTINE define_partics

!========================================================================

SUBROUTINE initialise_particle_arrays
!************************************************************************
!
! *initialise_particle_arrays* initialise particle arrays
!
! Author -  Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Initialisation.f90  V2.11
!
! $Date: 2018-06-14 15:25:04 +0200 (Thu, 14 Jun 2018) $
!
! $Revision: 1155 $
!
! Description -
!
! Reference -
!
! Calling program - particle_model
!
! External calls -
!
! Module calls - add_secs_to_date, noearlier, nolater
!
!************************************************************************
!
USE iopars
USE partpars
USE partswitches
USE partvars
USE syspars
USE timepars
USE time_routines

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=lentime) :: cdatetimex
INTEGER :: ntstart, numsteps, p
REAL (KIND=kndrlong) :: dsecs


procname(pglev+1) = 'initialise_particle_arrays'
CALL log_timer_in()

!
!1. Particle release time
!------------------------
!

p_110: DO p=1,nopart

   dsecs = part(p)%tstart*time_convert(ptime_unit)
   numsteps = NINT(dsecs/delt2d) 
   CALL add_secs_to_date(RefDateTime_part,cdatetimex,numsteps,delt2d)
   IF (for_drift) THEN
      IF (cdatetimex.nolater.CStartDateTime) THEN
         ntstart = 0
      ELSE
         CALL num_time_steps(CStartDateTime,cdatetimex,delt2d,ntstart)
      ENDIF
   ELSE
      IF (cdatetimex.noearlier.CStartDateTime) THEN
         ntstart = 0
      ELSE
         CALL num_time_steps(cdatetimex,CStartDateTime,delt2d,ntstart)
      ENDIF
   ENDIF
   part(p)%ntstart = ntstart
   part(p)%drift_state = 1
ENDDO p_110

!
!2. Location of particles released at initial time
!-------------------------------------------------
!

CALL particle_releases

!
!3. Contaminant properties
!-------------------------
!

IF (iopt_part_dens.GT.0) THEN
   volumeconc = (pi/6.0)*diamconc**3
   massconc = rhoconc*volumeconc
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE initialise_particle_arrays

!========================================================================

SUBROUTINE initialise_particle_model
!************************************************************************
!
! *initialise_particle_model* initialise particle model
!
! Author -  Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Initialisation.f90  V2.11
!
! $Date: 2018-06-14 15:25:04 +0200 (Thu, 14 Jun 2018) $
!
! $Revision: 1155 $
!
! Description -
!
! Reference -
!
! Calling program - particle_model
!
! External calls - allocate_part_arrays, combine_particle_grid,
!                  combine_particle_phys_data, define_partics, meteo_input,
!                  particle_concentrations, particle_trajects,
!                  read_cif_params, read_particle_grid, read_part_spec,
!                  time_series, update_particle_phys_data, usrdef_mod_params,
!                  usrdef_output, usrdef_part_output, usrdef_part_params,
!                  usrdef_part_spec, write_partics, write_part_spec,
!                  write_cif_vars_mod, write_cif_vars_mon, write_cif_vars_part,
!                  write_cif_vars_tspart
!
! Module calls - check_mod_params, check_part_params, check_partics,
!                close_filepars, comms_recv_int, default_mod_params,
!                default_part_params, initialise_time, open_filepars,
!                reset_mod_params, reset_part_params, reset_partics, rng_init
!
!************************************************************************
!
USE grid  
USE gridpars
USE iopars
USE modids
USE paralpars
USE partpars
USE partswitches
USE switches
USE check_model, ONLY: check_mod_params
USE check_particles, ONLY: check_part_params, check_partics 
USE comms_MPI, ONLY: comms_recv_int
USE default_model, ONLY: default_mod_params
USE default_particles, ONLY: default_part_params
USE inout_routines, ONLY: close_filepars, open_filepars
USE reset_model, ONLY: reset_mod_params
USE reset_particles, ONLY: reset_part_params, reset_partics
USE rng_library, ONLY: rng_init
USE time_routines, ONLY: initialise_time, log_timer_in, log_timer_out
  
IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc


procname(pglev+1) = 'initialise_particle_model'
CALL log_timer_in(npcc)

!
!1. Physical model parameters
!----------------------------
!
!1.1 Defaults
!------------
!

CALL default_mod_params

!
!1.2 Define
!----------
!

IF (ciffile%status.EQ.'R') THEN
   CALL read_cif_params(icif_mod)
ELSE
   CALL usrdef_mod_params
ENDIF

!
!1.3 Reset and check
!-------------------
!

CALL reset_mod_params
CALL check_mod_params

!
!2. Particle model
!-----------------
!
!---defaults
CALL default_part_params

!---define
IF (ciffile%status.EQ.'R') THEN
   CALL read_cif_params(icif_part)
ELSE 
   CALL usrdef_part_params
ENDIF

!---reset and check
CALL reset_part_params
CALL check_part_params

!
!3. Initialise
!-------------
!
!---time
CALL initialise_time

!---random generator
CALL rng_init

!
!4. Allocate arrays
!------------------
!

CALL allocate_part_arrays

!
!5. Domain decomposition
!-----------------------
!

IF (iopt_part_model.EQ.2) THEN
   CALL comms_recv_int(nc1procs,nprocscoh,idmaster,1,comm_world_MPI,&
                     & iarr_mg_nc1procs)
   CALL comms_recv_int(nc2procs,nprocscoh,idmaster,2,comm_world_MPI,&
                     & iarr_mg_nc2procs)
   CALL comms_recv_int(nr1procs,nprocscoh,idmaster,3,comm_world_MPI,&
                     & iarr_mg_nr1procs)
   CALL comms_recv_int(nr2procs,nprocscoh,idmaster,4,comm_world_MPI,&
                     & iarr_mg_nr2procs)
   ncprocs = nc2procs - nc1procs + 1
   nrprocs = nr2procs - nr1procs + 1
ENDIF
nc1loc = 1; nc2loc = nc
nr1loc = 1; nr2loc = nr


!
!6. Model grid
!-------------
!

IF (iopt_part_model.EQ.2) THEN
   CALL combine_particle_grid
ELSE
   CALL read_particle_grid
ENDIF

!
!7. Initial conditions
!---------------------
!
!---physical
IF (iopt_part_model.EQ.2) THEN
   CALL combine_particle_phys_data
ELSE
   CALL update_particle_phys_data
ENDIF

!---specifier arrays
IF (modfiles(io_parspc,1,1)%status.EQ.'R') THEN
   CALL read_part_spec
ELSEIF (modfiles(io_parspc,1,1)%status.EQ.'N') THEN
   CALL usrdef_part_spec
ENDIF
IF (modfiles(io_parspc,1,2)%status.EQ.'W') CALL write_part_spec

!---particles
CALL define_partics

!---reset and check   
CALL reset_partics
CALL check_partics

!
!8. Meteo data
!-------------
!

IF (iopt_part_wind.EQ.1) CALL meteo_input

!
!9. Initialise model output
!--------------------------
!
!---restart files
CALL write_partics

!---user-defined
CALL particle_concentrations
CALL particle_trajects
CALL usrdef_part_output
IF (iopt_part_model.GT.2) THEN
   IF (iopt_out_tsers.GT.0) CALL time_series
   CALL usrdef_output
ENDIF

!
!10. CIF
!-------
!

IF (iopt_part_model.GT.2) THEN

!
!10.1  Close CIF for reading
!---------------------------
!

   IF (ciffile%status.EQ.'R') CALL close_filepars(ciffile)

!
!10.2 Write CIF
!--------------
!

   IF (ciffile%status.EQ.'W') THEN
!     ---open file   
      CALL open_filepars(ciffile)
!     ---general/monitoring parameters 
      CALL write_cif_vars_mon
!     ---model setup
      CALL write_cif_vars_mod
!     ---particle model parameters 
      CALL write_cif_vars_part
!     ---time series output 
      IF (iopt_out_tsers.GT.0) CALL write_cif_vars_tsout
!     ---particle output
      CALL write_cif_vars_tspart
!     ---close CIF   
      CALL close_filepars(ciffile)
   ENDIF

ENDIF

CALL log_timer_out(npcc,itm_init)


RETURN

END SUBROUTINE initialise_particle_model

!========================================================================

SUBROUTINE particle_releases
!************************************************************************
!
! *particle_releases* initialise particle postion at release time
!
! Author -  Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Initialisation.f90  V2.11
!
! $Date: 2018-06-14 15:25:04 +0200 (Thu, 14 Jun 2018) $
!
! $Revision: 1155 $
!
! Description -
!
! Reference -
!
! Calling program - particle_model
!
! External calls -
!
! Module calls - Cint_at_p, error_abort, hrel_coords_curv_abs,
!                hrel_coords_rect_nonunif_abs, hrel_coords_rect_unif_abs
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE partpars
USE partvars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_abort
USE grid_interp, ONLY: hrel_coords_curv_abs, hrel_coords_rect_nonunif_abs, &
                     & hrel_coords_rect_unif_abs
USE particle_routines, ONLY: Cint_at_p
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, istat, j, k, p
REAL :: delxdat, delydat, depmeanp, deptotp, gsig, x, x0dat, y, y0dat
REAL (KIND=kndrlong) :: xpos, ypos, zpos


procname(pglev+1) = 'particle_releases'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!

x0dat = surfacegrids(igrd_model,1)%x0dat
y0dat = surfacegrids(igrd_model,1)%y0dat
delxdat = surfacegrids(igrd_model,1)%delxdat
delydat = surfacegrids(igrd_model,1)%delydat

!
!2. Particle location
!--------------------
!

p_210: DO p=1,nopart

   IF (part(p)%ntstart.EQ.nt) THEN

!     ---horizontal relative coordinates
      xpos = part(p)%xpos; ypos = part(p)%ypos
      SELECT CASE (iopt_grid_htype)
         CASE (1)
            CALL hrel_coords_rect_unif_abs(xpos,ypos,x0dat,y0dat,delxdat,&
                                      & delydat,(/1,1/),(/nc,nr/),istat,i,j,&
                                      & xcoord=x,ycoord=y)
         CASE (2)
            CALL hrel_coords_rect_nonunif_abs(xpos,ypos,gxcoord(1:nc,1),&
                                      & gycoord(1,1:nr),nc,nr,(/1,1/),&
                                      & (/nc,nr/),istat,i,j,xcoord=x,ycoord=y)
         CASE (3)
            CALL hrel_coords_curv_abs(xpos,ypos,gxcoord(1:nc,1:nr),&
                                    & gycoord(1:nc,1:nr),nc,nr,(/1,1/),&
                                    & (/nc,nr/),istat,2,i,j,&
                                    & xcoord=x,ycoord=y)
      END SELECT
      IF (istat.EQ.0) THEN
         part(p)%icoord = i; part(p)%jcoord = j
         part(p)%xcoord = x; part(p)%ycoord = y
      ENDIF

!     ---drift state      
      IF (istat.EQ.1) THEN
         part(p)%drift_state = 4
      ELSEIF (.NOT.maskatc(i,j)) THEN
         part(p)%drift_state = 3
      ELSE
         part(p)%drift_state = MERGE(2,1,part(p)%ntstart.EQ.0)
      ENDIF
      
!     ---vertical relative coordinates
      IF (part(p)%state.EQ.2.AND.part(p)%drift_state.LT.3) THEN
         SELECT CASE (part(p)%kdistype)
            CASE (1)
               k = NINT(part(p)%zpos)
               IF (iopt_grid_vtype.LT.3) THEN
                  gsig = gsigcoordatc(k)
               ELSE
                  gsig = gscoordatc(i,j,k)
               ENDIF
               part(p)%zpos = gsig*deptotatc(i,j) - depmeanatc(i,j)
            CASE (2)
               deptotp = Cint_at_p(deptotatc(i-1:i+1,j-1:j+1),i,j,x,y)
               depmeanp = Cint_at_p(depmeanatc(i-1:i+1,j-1:j+1),i,j,x,y)
               zpos = part(p)%zpos
               zpos = MAX(0.0,MIN(zpos,deptotp))
               part(p)%zpos = zpos - depmeanp
            CASE (3)
               deptotp = Cint_at_p(deptotatc(i-1:i+1,j-1:j+1),i,j,x,y)
               depmeanp = Cint_at_p(depmeanatc(i-1:i+1,j-1:j+1),i,j,x,y)
               zpos = deptotp - part(p)%zpos
               zpos = MAX(0.0,MIN(zpos,deptotp))
               part(p)%zpos = zpos - depmeanp
         END SELECT
      ELSE
         part(p)%zpos = 0.0
      ENDIF
      part(p)%kdistype = 0

   ENDIF
   
ENDDO p_210

CALL error_abort('particle_releases',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE particle_releases

!========================================================================

SUBROUTINE read_partics
!************************************************************************
!
! *read_partics* Read initial conditions for the particle model in standard
!                format
!
! Author -  Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Initialisation.f90  V2.11
!
! $Date: 2018-06-14 15:25:04 +0200 (Thu, 14 Jun 2018) $
!
! $Revision: 1155 $
!
! Description -
!
! Reference -
!
! Calling program - particle_model
!
! External calls -
!
! Module calls -  close_filepars, error_alloc_struc, open_filepars,
!                 read_glbatts_mod, read_time, read_varatts_mod, read_vars,
!                 varatts_init 
!
!************************************************************************
!
USE datatypes
USE iopars
USE partswitches
USE partvars
USE switches
USE syspars
USE timepars
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_time, read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out


IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: fcomm
CHARACTER (LEN=lentime) :: ciodatetime
INTEGER :: numvars, varid
TYPE(FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_partics'
CALL log_timer_in()

filepars = modfiles(io_inicon,ics_part,1)

!
!1. File header
!--------------
!

CALL open_filepars(filepars)
CALL read_glbatts_mod(filepars)
numvars = filepars%novars + 1
ALLOCATE (varatts(numvars),STAT=errstat)
CALL varatts_init(varatts)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL read_varatts_mod(filepars,varatts,numvars)

!
!2. Read particle data
!---------------------
!

fcomm = .FALSE.

DO WHILE (.NOT.fcomm)

!
!2.1 Date/time
!-------------
!

   CALL read_time(ciodatetime,filepars)
   fcomm = ciodatetime.EQ.CDateTime
   varid = 1

!
!2.2 Initial data
!-----------------
!

   varid = varid + 1
   CALL read_vars(part%state,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(part%drift_state,filepars,varid,varatts=varatts)
   varid = varid + 1
   IF (iopt_grid_nodim.NE.2) THEN
      CALL read_vars(part%kdistype,filepars,varid,varatts=varatts)
      varid = varid + 1
   ENDIF
   CALL read_vars(part%label,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(part%xpos,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(part%ypos,filepars,varid,varatts=varatts)
   varid = varid + 1
   IF (iopt_grid_nodim.NE.2) THEN
      CALL read_vars(part%zpos,filepars,varid,varatts=varatts)
      varid = varid + 1
   ENDIF
   CALL read_vars(part%tstart,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(part%age,filepars,varid,varatts=varatts)

ENDDO

!
!3. Close
!--------
!

CALL close_filepars(filepars)
DEALLOCATE (varatts)

modfiles(io_inicon,ics_part,1) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE read_partics

!========================================================================

SUBROUTINE read_part_spec
!************************************************************************
!
! *read_part_spec* Read specifiers arrays for particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Initialisation.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE partswitches
USE partvars
USE switches
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_par_spec'
CALL log_timer_in()

filepars = modfiles(io_parspc,1,1)

!
!1. File header
!--------------
!

CALL open_filepars(filepars)
CALL read_glbatts_mod(filepars)
numvars = filepars%novars
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL read_varatts_mod(filepars,varatts,numvars)

!
!2. Read data
!------------
!
!2.1 Particle properties
!-----------------------
!

varid = 0
IF (iopt_part_conc.GT.1.OR.iopt_part_dens.GT.0) THEN
   varid = varid + 1
   CALL read_vars(diamconc,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(rhoconc,filepars,varid,varatts=varatts)
ENDIF

!
!2.2 Particle cloud parameters
!-----------------------------
!

IF (iopt_part_cloud.EQ.1) THEN
   
   IF (iopt_grid_nodim.NE.2) THEN
      varid = varid + 1
      CALL read_vars(cloud%kdistype,filepars,varid,varatts=varatts)
   ENDIF
   varid = varid + 1
   CALL read_vars(cloud%label,filepars,varid,varatts=varatts)
   varid = varid + 1   
   CALL read_vars(cloud%length,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(cloud%nopart,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(cloud%noreleases,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(cloud%orientation,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(cloud%state,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(cloud%thick,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(cloud%width,filepars,varid,varatts=varatts)

ENDIF
   
!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_parspc,1,1) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE read_part_spec
