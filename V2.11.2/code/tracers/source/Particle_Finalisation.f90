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
! *Particle_Finalisation* Finalise sediment model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Finalisation.f90  V2.11
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description -
!
! Reference -
!
! Routines - write_partics, write_part_spec
!
!************************************************************************
!

!========================================================================

SUBROUTINE write_partics
!************************************************************************
!
! *write_partics* Write initial conditions for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Finalisation.f90  V2.10.1
!
! Description - 
!
! Reference -
!
! Calling program - coherens_main, initialise_model, initialise_particle_model,
!                   particle_model
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, write_atts_mod,
!                write_time, write_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE partvars
USE switches
USE timepars
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_time, write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out  

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
INTEGER :: numvars, varid
INTEGER, SAVE :: itrestart
TYPE(FileParams) :: filepars
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


IF (norestarts.EQ.0.OR.modfiles(io_fincon,ics_part,2)%status.EQ.'0') RETURN
IF (nt.EQ.0) itrestart = 0
IF (itrestart.GT.norestarts) RETURN

procname(pglev+1) = 'write_partics'
flag = .FALSE.
filepars = modfiles(io_fincon,ics_part,2)

!
!1. Write file header on first call
!----------------------------------
!

IF (nt.EQ.ntrestart(1)) THEN

   CALL log_timer_in()
   flag = .TRUE.

!  ---file attributes
   CALL set_modfiles_atts(io_fincon,ics_part,2)
   filepars = modfiles(io_fincon,ics_part,2)
   numvars = filepars%novars + 1

!  ---open file
   CALL open_filepars(filepars)

!  ---variable attributes
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
   CALL set_modvars_atts(io_fincon,ics_part,2,varatts,numvars)

!  ---write attributes
   CALL write_atts_mod(filepars,varatts,numvars)

   itrestart = 1

ENDIF

!
!2. Write physical data
!---------------------
!

IF (itrestart.GT.0.AND.itrestart.LE.norestarts) THEN

   IF (nt.EQ.ntrestart(itrestart)) THEN
      IF (itrestart.GT.1) THEN
         CALL log_timer_in()
         flag = .TRUE.
      ENDIF

!
!2.1 Date/time
!-------------
!

      CALL write_time(CDateTime,filepars)
      varid = 1

!
!2.2 Initial data
!-----------------
!

      varid = varid + 1
      CALL write_vars(part%state,filepars,varid,varatts=varatts)
      varid = varid + 1
      CALL write_vars(part%drift_state,filepars,varid,varatts=varatts)
      varid = varid + 1
      IF (iopt_grid_nodim.NE.2) THEN
         CALL write_vars(part%kdistype,filepars,varid,varatts=varatts)
         varid = varid + 1
      ENDIF
      CALL write_vars(part%label,filepars,varid,varatts=varatts)
      varid = varid + 1
      CALL write_vars(part%xpos,filepars,varid,varatts=varatts)
      varid = varid + 1
      CALL write_vars(part%ypos,filepars,varid,varatts=varatts)
      varid = varid + 1
      IF (iopt_grid_nodim.NE.2) THEN
         CALL write_vars(part%zpos,filepars,varid,varatts=varatts)
         varid = varid + 1
      ENDIF
      CALL write_vars(part%tstart,filepars,varid,varatts=varatts)
      varid = varid + 1
      CALL write_vars(part%age,filepars,varid,varatts=varatts)

!
!2.3 Next output index
!---------------------
!

      itrestart = itrestart + 1

   ENDIF

ENDIF

!
!3. Finalise
!-----------
!

IF ((cold_start.AND.ntrestart(1).EQ.0).OR.&
  & (itrestart.GT.norestarts)) THEN
   CALL close_filepars(filepars)
   DEALLOCATE (varatts)
ENDIF

modfiles(io_fincon,ics_part,2) = filepars

IF (flag) CALL log_timer_out() 


RETURN

END SUBROUTINE write_partics

!========================================================================

SUBROUTINE write_part_spec
!************************************************************************
!
! *write_part_spec* Write specifiers arrays for particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Finalisation.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, write_atts_mod,
!                write_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE partswitches
USE partvars
USE switches
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'write_part_spec'
CALL log_timer_in()

!
!1. Write file header
!--------------------
!
!---file attributes
CALL set_modfiles_atts(io_parspc,1,2)
filepars = modfiles(io_parspc,1,2)
numvars = filepars%novars

!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL set_modvars_atts(io_parspc,1,2,varatts,numvars)

!---write
CALL write_atts_mod(filepars,varatts,numvars)

!
!2.1 Particle properties
!-----------------------
!

varid = 0
IF (iopt_part_conc.GT.1.OR.iopt_part_dens.GT.0) THEN
   varid = varid + 1
   CALL write_vars(diamconc,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(rhoconc,filepars,varid,varatts=varatts)
ENDIF

!
!2.2 Particle cloud parameters
!-----------------------------
!

IF (iopt_part_cloud.EQ.1) THEN
   
   IF (iopt_grid_nodim.NE.2) THEN
      varid = varid + 1
      CALL write_vars(cloud%kdistype,filepars,varid,varatts=varatts)
   ENDIF
   varid = varid + 1
   CALL write_vars(cloud%label,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(cloud%length,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(cloud%nopart,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(cloud%noreleases,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(cloud%orientation,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(cloud%state,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(cloud%thick,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(cloud%width,filepars,varid,varatts=varatts)

ENDIF

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_parspc,1,2) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE write_part_spec
