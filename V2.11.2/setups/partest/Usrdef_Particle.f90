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
! Description - test case partest
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
! Description - test case partest
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
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_part_params'
CALL log_timer_in()

!
!1. Particle switches
!--------------------
!
!---time integration
SELECT CASE (runtitle(8:8))
   CASE ('C','E'); iopt_part_vdif = 2
   CASE DEFAULT; iopt_part_vdif = 0
END SELECT

!---output
SELECT CASE (runtitle(8:8))
   CASE ('D','E'); iopt_part_out = 0
   CASE DEFAULT; iopt_part_out = 1
END SELECT

!---concentrations
iopt_part_conc = 1

!
!2. Partcicle parameters
!-----------------------
!

SELECT CASE (runtitle(8:8))
   CASE ('D','E'); nopart = 3000
   CASE DEFAULT; nopart = 30   
END SELECT

nolabels = 3

nosetspart = MERGE(1,0,iopt_part_out.EQ.1)
ptime_unit = 3

IF (iopt_part_model.EQ.3) THEN
   icpart = 4
ELSEIF (iopt_part_model.EQ.4) THEN
   icpart = 1
ENDIF

IF (iopt_part_model.LT.4) THEN
   RefDateTime_part = '2003/01/02;00:00:00:000'
ELSE
   RefDateTime_part = '2003/01/03;00:00:00:000'
ENDIF

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
! Description - test case partest
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!
USE depths
USE gridpars  
USE iopars
USE partvars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: immersed
INTEGER :: i, j, k, l
REAL :: delh, delta


procname(pglev+1) = 'usrdef_partics'
CALL log_timer_in()

!
!1. State
!--------
!

part%state = MERGE(1,2,runtitle(8:8).EQ.'A'.OR.runtitle(8:8).EQ.'G')
part%tstart = 0.0

!
!2. Initial position
!-------------------
!

delh = surfacegrids(igrd_model,1)%delxdat
immersed = part(1)%state.EQ.2

SELECT CASE (runtitle(8:8))

   CASE ('A','B','C','F','G')
      part(1:10)%xpos = delh*19.5
      part(11:20)%xpos = delh*44.5
      part(21:30)%xpos = delh*94.5
      part%ypos = delh*5.0
      part(1:10)%label = 1
      part(11:20)%label = 2
      part(21:30)%label = 3
      IF (immersed) THEN
         part%kdistype = 1
         part(1:10)%zpos = (/(k,k=1,nz-1,2)/)
         part(11:20)%zpos = part(1:10)%zpos
         part(21:30)%zpos = part(1:10)%zpos
      ENDIF

   CASE ('D','E')
      part(1:1000)%xpos = delh*19.5
      part(1001:2000)%xpos = delh*44.5
      part(2001:3000)%xpos = delh*94.5
      part%ypos = delh*5.0
      part(1:1000)%label = 1
      part(1001:2000)%label = 2
      part(2001:3000)%label = 3
      part%kdistype = 0
      i = 20; j = 5; delta = depmeanatc(i,j)/1000.0
      part(1)%zpos = -depmeanatc(i,j) + 0.5*delta
      l_211: DO l=2,1000
         part(l)%zpos = part(1)%zpos + (l-1)*delta
      ENDDO l_211
      i = 45; delta = depmeanatc(i,j)/1000.0
      part(1001)%zpos = -depmeanatc(i,j) + 0.5*delta
      l_212: DO l=1002,2000
         part(l)%zpos = part(1001)%zpos + (l-1001)*delta
      ENDDO l_212
      i = 95; delta = depmeanatc(i,j)/1000.0
      part(2001)%zpos = -depmeanatc(i,j) + 0.5*delta
      l_213: DO l=2002,3000
         part(l)%zpos = part(2001)%zpos + (l-2001)*delta
      ENDDO l_213

END SELECT

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
! Version - @(COHERENS)Usrdef_Particle.f90  V2.11.1
!
! Description - test case partest
!
! Reference -
!
! Calling program - coherens_main, initialise_model, initialise_particle_model,
!                   particle_model
!  
! Module calls - close_file, distance_pts, error_alloc, mult_index, open_file
!
!************************************************************************
!
USE iopars  
USE partpars
USE partswitches
USE partvars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: distance_pts
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL, SAVE :: flag, immersed
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
CHARACTER (LEN=1) :: cl
CHARACTER (LEN=2) :: cp
INTEGER, SAVE :: icout, iunit, notracks = 9, state
INTEGER :: l, p, pp
INTEGER (KIND=kndilong), SAVE :: nosecsini
REAL :: dist, rtime
INTEGER (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:) :: ptrack
REAL (KIND=kndrlong), DIMENSION(nolabels) :: xcen, xdev, ycen, ydev, zcen, zdev
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:) :: htrack, vtrack, &
                                          & xtrack_old, ytrack_old, ztrack_old



IF (ciffile%status.EQ.'W') RETURN

!
!1. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN
   
   procname(pglev+1) = 'usrdef_part_output'
   CALL log_timer_in()
   
!  ---open output file
   CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
   WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
   WRITE (iunit,'(A)') 'Output parameters for test case partest: '//&
                     & ' simulation '//TRIM(runtitle)
   WRITE (iunit,*)

!  ---particle state
   immersed = part(1)%state.EQ.2

!  ---flag for small release
   flag = nopart.EQ.30
   
!  ---output ferequency
   icout = MERGE(NINT(3600.0/delt2d),NINT(3600.0/delt3d),iopt_part_model.LE.2)

!  ---initial time   
   IF (iopt_part_model.LE.2) THEN
      nosecsini = INT(nstep/2,KIND=kndilong)*NINT(delt2d,KIND=kndilong)
   ELSEIF (iopt_part_model.EQ.3) THEN
      nosecsini = 0
   ELSE
      nosecsini = INT(nstep,KIND=kndilong)*NINT(delt3d,KIND=kndilong)
   ENDIF
   
!  ---allocate
   IF (flag) THEN
      ALLOCATE(htrack(notracks),STAT=errstat)
      CALL error_alloc('htrack',1,(/notracks/),kndrlong)
      ALLOCATE(ptrack(notracks),STAT=errstat)
      CALL error_alloc('ptrack',1,(/notracks/),kndint)
      ptrack = (/1,5,10,11,15,20,21,25,30/)
      ALLOCATE(xtrack_old(notracks),STAT=errstat)
      CALL error_alloc('xtrack_old',1,(/notracks/),kndrlong)
      ALLOCATE(ytrack_old(notracks),STAT=errstat)
      CALL error_alloc('ytrack_old',1,(/notracks/),kndrlong)
      IF (immersed) THEN
         ALLOCATE(vtrack(notracks),STAT=errstat)
         CALL error_alloc('vtrack',1,(/notracks/),kndrlong)
         ALLOCATE(ztrack_old(notracks),STAT=errstat)
         CALL error_alloc('ztrack_old',1,(/notracks/),kndrlong)
      ENDIF
   ENDIF

   CALL log_timer_out()
   
ENDIF

!
!2. Return if no tracking
!----------------------
!

IF (nt.LT.part(1)%ntstart) RETURN

procname(pglev+1) = 'usrdef_part_output'
CALL log_timer_in()

!
!3. Without diffusion
!--------------------
!

IF (flag) THEN

!  ---locations and travelled distance   
   p_310: DO p=1,notracks
      pp = ptrack(p)
      IF (nt.EQ.part(pp)%ntstart) THEN
         htrack(p) = 0.0
         xtrack_old(p) = part(pp)%xpos
         ytrack_old(p) = part(pp)%ypos
         IF (immersed) THEN
            vtrack(p) = 0.0
            ztrack_old(p) = part(pp)%zpos
         ENDIF
      ELSE
         dist = distance_pts(xtrack_old(p),part(pp)%xpos,ytrack_old(p),&
                           & part(pp)%ypos,2,3)
         htrack(p) = htrack(p) + 0.001*dist
         xtrack_old(p) = part(pp)%xpos
         ytrack_old(p) = part(pp)%ypos
         IF (immersed) THEN
            vtrack(p) = vtrack(p) + ABS(part(pp)%zpos-ztrack_old(p))
            ztrack_old(p) = part(pp)%zpos
         ENDIF
      ENDIF
   ENDDO p_310

!  ---write output
   IF (mult_index(nt,icout)) THEN
      rtime = (nosecsrun-nosecsini)/3600.0
      IF (iopt_part_model.EQ.4) rtime = -rtime
      WRITE (iunit,9001) rtime
      pp_320: DO p=1,notracks
         pp = ptrack(p)
         WRITE (cp,'(I2.2)') pp; cp = ADJUSTL(cp)
         WRITE (iunit,9002) 'xpos'//cp, part(pp)%xpos
         WRITE (iunit,9002) 'ypos'//cp, part(pp)%ypos
         IF (immersed) WRITE (iunit,9002) 'zpos'//cp, part(pp)%zpos
         WRITE (iunit,9002) 'htrack'//cp, htrack(p)
         IF (immersed) WRITE (iunit,9002) 'vtrack'//cp, vtrack(p)
      ENDDO pp_320
   ENDIF

!
!4. Tests C and D
!----------------
!   

ELSEIF (mult_index(nt,icout)) THEN

!  ---mean location for each distribution
   xcen = 0; ycen = 0.0; zcen = 0.0
   p_410: DO p=1,nopart
      l = part(p)%label
      xcen(l) = xcen(l) + part(p)%xpos
      ycen(l) = ycen(l) + part(p)%ypos
      IF (immersed) zcen(l) = zcen(l) + part(p)%zpos
   ENDDO p_410
   xcen = 0.001*xcen; ycen = 0.001*ycen
   IF (immersed) zcen = 0.001*zcen
   
!  ---mean dispersion
   xdev = 0.0; ydev = 0.0
   IF (immersed) zdev = 0.0
   p_420: DO p=1,nopart
      l = part(p)%label
      xdev(l) = xdev(l) + (xcen(l)-part(p)%xpos)**2
      ydev(l) = ydev(l) + (ycen(l)-part(p)%ypos)**2
      IF (immersed) zdev(l) = zdev(l) + (zcen(l)-part(p)%zpos)**2
   ENDDO p_420
   xdev = SQRT(xdev/999.0); ydev = SQRT(ydev/999.0)
   IF (immersed) zdev = SQRT(zdev/999.0)

!  ---write output
   rtime = (nosecsrun-nosecsini)/3600.0
   WRITE (iunit,9001) rtime
   l_430: DO l=1,nolabels
      WRITE (cl,'(I1)') l; cl = ADJUSTL(cl)
      WRITE (iunit,9002) 'xcen'//cl, 0.001*xcen(l)
      WRITE (iunit,9002) 'ycen'//cl, 0.001*ycen(l)
      IF (immersed) WRITE (iunit,9002) 'zcen'//cl, zcen(l)
      WRITE (iunit,9002) 'xdev'//cl, 0.001*xdev(l)
      WRITE (iunit,9002) 'ydev'//cl, 0.001*ydev(l)
      IF (immersed) WRITE (iunit,9002) 'zdev'//cl, zdev(l)
   ENDDO l_430

ENDIF
   
!
!5. Close output file
!--------------------
!

IF (nt.EQ.nstep) THEN
   CALL close_file(iunit,'A')

!
!6. Deallocate
!-------------
!
   
   IF (flag) THEN
      DEALLOCATE (htrack,ptrack,xtrack_old,ytrack_old)
      IF (immersed) DEALLOCATE (vtrack,ztrack_old)
   ENDIF
   
ENDIF

CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',F5.2,' hours')
9002 FORMAT (T3,A,T12,': ',G12.4E3)

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
! Description - test case partest
!
! Reference -
!
! Calling program - particle_trajects_init
!
!************************************************************************
!
USE iopars  
USE partvars
USe switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_pout_params'
CALL log_timer_in()

!---file parameters  
part1d(1)%defined = .TRUE.

!---output attributes
IF (runtitle(8:8).EQ.'A'.OR.runtitle(8:8).EQ.'G') THEN
   outppars(1)%nodim = 2
ELSE
   outppars(1)%nodim = 0
ENDIF
outppars(1)%refdate = '2003/01/02;00:00:00:000'
outppars(1)%tlims(1) = MERGE(nstep/2,0,iopt_part_model.LE.2)
outppars(1)%tlims(2:3) = (/nstep,ic3d/)
outppars(1)%aging = .TRUE.

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
