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
! *Usrdef_Particle* User-defined setup for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11
!
! $Date: 2016-10-27 16:15:25 +0200 (Thu, 27 Oct 2016) $
!
! $Revision: 986 $
!
! Description - test case partrot
!
! Reference -
!
! Routines - usrdef_part_params, usrdef_partics, usrdef_part_output,
!            usrdef_part_spec, usrdef_pout_params, usrdef_parcld_data
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_part_params
!************************************************************************
!
! *usrdef_part_params* Define parameters for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11
!
! Description - test case partrot
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!
USE iopars  
USE partpars
USE partswitches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_part_params'
CALL log_timer_in()

!
!1. Particle switches
!--------------------
!
!---time integration
iopt_part_hadv = MERGE (1,2,runtitle(8:8).EQ.'B')

!---output
iopt_part_out = 1

!
!2. Partcicle parameters
!-----------------------
!

nopart = 8
nosetspart = 1
ptime_unit = 3
SELECT CASE (runtitle(8:8))
   CASE ('C','D')
      icpart = 480
END SELECT

RefDateTime_part = '2003/01/01;00:00:00'

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_part_params

!========================================================================

SUBROUTINE usrdef_partics
!************************************************************************
!
! *usrdef_partics* Define initial conditions for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11
!
! Description - test case partrot
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!
USE iopars  
USE partpars
USE partvars
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j, p
INTEGER, DIMENSION(nopart) :: istart, jstart
REAL :: delx, dely


procname(pglev+1) = 'usrdef_partics'
CALL log_timer_in()

!---particle state
part%state = 1
part%drift_state = 1
part%tstart = MERGE(0.0,24.0,iopt_part_model.LT.4)

!---initial position
istart= (/9,17,25,33,49,57,65,73/)
jstart = 41
delx = surfacegrids(igrd_model,1)%delxdat
dely = surfacegrids(igrd_model,1)%delydat
p_110: DO p=1,nopart
   i = istart(p); j = jstart(p)
   part(p)%xpos = (i-0.5)*delx
   part(p)%ypos = (j-0.5)*dely
ENDDO p_110

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_partics

!========================================================================

SUBROUTINE usrdef_part_output
!************************************************************************
!
! *usrdef_part_output* User-formatted output for the particle model 
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11
!
! Description - test case partrot
!
! Reference -
!
! Calling program - coherens_main, initialise_model, initialise_particle_model,
!                   particle_model
!  
! Module calls - close_file, mult_index, open_file  
!
!************************************************************************
!
USE iopars  
USE partpars  
USE partvars
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: close_file, open_file
USE math_library, ONLY: complex_polar
USE time_routines, ONLY: diff_dates, log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
CHARACTER (LEN=1) :: cp
INTEGER, SAVE :: icout, iunit
INTEGER :: p
REAL, SAVE :: omega, xc, yc 
REAL :: dphase, phase, pphase, prad, xanal, xpos, yanal, ypos
REAL (KIND=kndrlong) :: mtime
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: phase0, radius


!
!1. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN

!  ---open output file
   CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
   WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
   WRITE (iunit,'(A)') 'Output parameters for test case partrot: '//&
                     & ' simulation '//TRIM(runtitle)
   WRITE (iunit,*)

!  ---output ferequency
   icout = NINT(900.0/delt3d)

!  ---location of rotation center
   xc = 20.25; yc = 20.25

!  ---angular speed
   omega = MERGE(pi/3600.0,-pi/3600.0,for_drift)

!  ---allocate
   ALLOCATE(radius(nopart),STAT=errstat)
   CALL error_alloc('radius',1,(/nopart/),kndrtype)
   ALLOCATE(phase0(nopart),STAT=errstat)
   CALL error_alloc('phase0',1,(/nopart/),kndrtype)
   
!  ---initial position in polar coordinates
   p_110: DO p=1,nopart
      xpos = part(p)%xpos
      radius(p) = ABS(xpos-xc)
      phase0(p) = MERGE(0.0,-pi,xpos.GT.xc)
   ENDDO p_110
   
ENDIF

!
!2.Write output
!--------------
!


IF (mult_index(nt,icout)) THEN

   procname(pglev+1) = 'usrdef_part_output'
   CALL log_timer_in()

   CALL diff_dates(RefDateTime_part,CDateTime,ptime_unit,rtime=mtime)

   WRITE (iunit,9001) mtime
   p_210: DO p=1,nopart
      WRITE (cp,'(I1)') p
      phase = MOD(phase0(p)+omega*nosecsrun,twopi)
      xanal = radius(p)*COS(phase)
      yanal = radius(p)*SIN(phase)
      xpos = part(p)%xpos - xc
      ypos = part(p)%ypos - yc
      CALL complex_polar(xpos,ypos,xamp=prad,xpha=pphase)
      dphase = pphase - phase
      dphase = MOD(dphase,twopi)
      IF (dphase.GE.pi) dphase = dphase - twopi
      IF (dphase.LE.-pi) dphase = dphase + twopi
      WRITE (iunit,9002) 'xpos'//cp, xpos
      WRITE (iunit,9002) 'ypos'//cp, ypos
      WRITE (iunit,9002) 'drad'//cp, (prad-radius(p))/radius(p)
      WRITE (iunit,9003) 'dphi'//cp, radtodeg*dphase
   ENDDO p_210
   
   CALL log_timer_out()

ENDIF

!
!3. Close output file
!--------------------
!

IF (cold_start.OR.nt.EQ.nstep) THEN
   CALL close_file(iunit,'A')
   DEALLOCATE (radius,phase0)
ENDIF

RETURN

9001 FORMAT ('Time: ',F5.2,' hours')
9002 FORMAT (T3,A,T12,': ',G12.4E3)
9003 FORMAT (T3,A,T12,': ',3G12.4E3)

END SUBROUTINE usrdef_part_output

!========================================================================

SUBROUTINE usrdef_part_spec
!************************************************************************
!
! *usrdef_part_spec* Define specifier arrays for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_part_spec

!========================================================================

SUBROUTINE usrdef_pout_params
!************************************************************************
!
! *usrdef_pout_params* Define parameters for particle trajectory output
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11
!
! Description - test case partrot
!
! Reference -
!
! Calling program - particle_trajects_init
!
!************************************************************************
!
USE iopars  
USE partvars  
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_pout_params'
CALL log_timer_in()

!---file parameters  
part1d(1)%defined = .TRUE.

!---output attributes
outppars(1)%nodim = 2
outppars(1)%tlims = (/0,nstep,1/)

CALL log_timer_out()  


RETURN

END SUBROUTINE usrdef_pout_params

!========================================================================

SUBROUTINE usrdef_parcld_data(ifil,ciodatetime,disdata)
!************************************************************************
!
! *usrdef_parcld_data* Define locations for cloud discharges
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - define_parcld_data
!
!************************************************************************
!
USE syspars
  
IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(OUT) :: ciodatetime
INTEGER, INTENT(IN) :: ifil
REAL (KIND=kndrlong), INTENT(INOUT), DIMENSION(3) :: disdata

!
! Name         Type Purpose
!------------------------------------------------------------------------------
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*disdata*     REAL    Input data
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_parcld_data
