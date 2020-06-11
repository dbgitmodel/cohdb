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
! Description - test case partrhone
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
! Description - test case partrhone
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
!---particle cloud
iopt_part_cloud = 1   
   
!---output
iopt_part_out = 1

!---concentrations
iopt_part_conc = 1

!---wind forcing
iopt_part_wind = MERGE(0,1,runtitle(10:10).EQ.'A')

!---horizontal diffusion
iopt_part_hdif = MERGE(1,0,runtitle(10:10).EQ.'C')

!---leeway drift
iopt_part_leeway = MERGE(1,0,runtitle(10:10).EQ.'D')

!
!2. Partcicle parameters
!-----------------------
!
!---total number of particles
nopart = 9000

!---particle labels
nolabels = 3

!---particle clouds
noclouds = 3

!--output sets
nosetspart = 1

!---horizontal diffusion coefficient
IF (iopt_part_hdif.EQ.1) THEN
   xdifpart_cst = 10.0; ydifpart_cst = 10.0
ENDIF

!---time unit
ptime_unit = 3

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
! Description -
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!


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
! Description - test case partrhone
!
! Reference -
!
! Calling program - coherens_main, initialise_model, initialise_particle_model,
!                   particle_model
!  
! Module calls - close_file, distance_pts, error_alloc, loop_index, open_file
!
!************************************************************************
!
USE iopars
USE partpars
USE partvars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: distance_pts
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: loop_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
CHARACTER (LEN=1) :: crel
INTEGER, SAVE :: iunit, nopartrel, numrel
INTEGER :: icld, irel, l, p, pp
INTEGER, SAVE, DIMENSION(3) :: ntout
REAL, PARAMETER :: rflag = -999.9 
REAL :: dist, rtime
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: pcen
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:) :: htrack, xdev, xrel, &
                                           & xtrack_old, ydev, yrel, ytrack_old


IF (ciffile%status.EQ.'W') RETURN

!
!1. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN
   
   procname(pglev+1) = 'usrdef_part_output'
   CALL log_timer_in()

!  ---open output file
   IF (.NOT.cold_start) THEN
      CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case partrhone: '//&
                        & ' simulation '//TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF

!  ---number of releases, particles per release
   numrel = 0
   icld_110: DO icld=1,noclouds
      numrel = numrel + cloud(icld)%noreleases
   ENDDO icld_110
   nopartrel = cloud(1)%nopart
   
!  ---allocate
   ALLOCATE (pcen(numrel),STAT=errstat)
   CALL error_alloc('pcen',1,(/numrel/),kndint)
   ALLOCATE (xrel(numrel),STAT=errstat)
   CALL error_alloc('xrel',1,(/numrel/),kndrlong)
   ALLOCATE (yrel(numrel),STAT=errstat)
   CALL error_alloc('yrel',1,(/numrel/),kndrlong)
   ALLOCATE (xdev(numrel),STAT=errstat)
   CALL error_alloc('xdev',1,(/numrel/),kndrlong)
   ALLOCATE (ydev(numrel),STAT=errstat)
   CALL error_alloc('ydev',1,(/numrel/),kndrlong)
   ALLOCATE (htrack(numrel),STAT=errstat)
   CALL error_alloc('htrack',1,(/numrel/),kndrlong)
   htrack = 0.0
   ALLOCATE (xtrack_old(numrel),STAT=errstat)
   CALL error_alloc('xtrack_old',1,(/numrel/),kndrlong)
   ALLOCATE (ytrack_old(numrel),STAT=errstat)
   CALL error_alloc('ytrack_old',1,(/numrel/),kndrlong)
   
!  ---particle index at cloud center
   irel = 0
   p_120: DO p=1,nopart
      IF (part(p)%center) THEN
         irel = irel + 1
         pcen(irel) = p
      ENDIF
   ENDDO p_120

!  ---output times
   ntout(1) = part(pcen(1))%ntstart
   ntout(2) = nstep
   ntout(3) = MERGE(720,40,iopt_hydro_impl.EQ.0)
   
   CALL log_timer_out()
   
ENDIF

IF (nt.LT.ntout(1)) RETURN

procname(pglev+1) = 'usrdef_part_output'
CALL log_timer_in()

!
!2. Update travelled distance
!----------------------------
!

p_210: DO irel=1,numrel
   p = pcen(irel)
   IF (nt.EQ.part(p)%ntstart) THEN
      htrack(irel) = 0.0
      xtrack_old(irel) = part(p)%xpos
      ytrack_old(irel) = part(p)%ypos
   ELSEIF (part(p)%drift_state.EQ.2) THEN
      dist = distance_pts(xtrack_old(irel),part(p)%xpos,ytrack_old(irel),&
                        & part(p)%ypos,2,3)
      htrack(irel) = htrack(irel) + 0.001*dist
      xtrack_old(irel) = part(p)%xpos
      ytrack_old(irel) = part(p)%ypos
   ENDIF
ENDDO p_210

!
!3. Store/write output parameters
!--------------------------------
!

IF (loop_index(ntout,nt)) THEN

!  ---time
   rtime = nosecsrun/3600.0   

!  ---location of cloud centers
   irel_310: DO irel=1,numrel
      p = pcen(irel)
      xrel(irel) = part(p)%xpos
      yrel(irel) = part(p)%ypos
   ENDDO irel_310

!  ---cloud spreading
   irel_320: DO irel=1,numrel
      p = pcen(irel)
      xdev(irel) = 0.0; ydev(irel) = 0.0
      pp_321: DO pp=p+1,p+nopartrel
         xdev(irel) = xdev(irel) + (xrel(irel)-part(pp)%xpos)**2
         ydev(irel) = ydev(irel) + (yrel(irel)-part(pp)%ypos)**2
      ENDDO pp_321
      xdev(irel) = SQRT(xdev(irel)/(nopartrel-1))
      ydev(irel) = SQRT(ydev(irel)/(nopartrel-1))
   ENDDO irel_320

!  ---write
   WRITE (iunit,9001) rtime
   irel_330: DO irel=1,numrel
      p = pcen(irel)
!      IF (part(p)%drift_state.EQ.2) THEN
         WRITE (crel,'(I1)') irel; crel = ADJUSTL(crel)
         WRITE (iunit,9002) 'xpos'//crel, xrel(irel)
         WRITE (iunit,9002) 'ypos'//crel, yrel(irel)
         WRITE (iunit,9002) 'htrack'//crel, htrack(irel)
         WRITE (iunit,9002) 'xdev'//crel, xdev(irel)
         WRITE (iunit,9002) 'ydev'//crel, ydev(irel)
!      ENDIF
   ENDDO irel_330

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
   
   DEALLOCATE (pcen,xrel,yrel,xdev,ydev,htrack,xtrack_old,ytrack_old)
   
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
! Description - test case partrhone
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!
USE iopars
USE partvars  
USE time_routines, ONLY: log_timer_in, log_timer_out
  
  
procname(pglev+1) = 'usrdef_part_spec'
CALL log_timer_in()

cloud%kdistype = 3
cloud%label = (/1,2,3/)
cloud%length = 0.5*surfacegrids(igrd_model,1)%delydat
cloud%nopart = 1000
cloud%noreleases = 3
cloud%orientation = 0.0
cloud%state = 1
cloud%thick = 0.0
cloud%width = 0.5*surfacegrids(igrd_model,1)%delxdat

CALL log_timer_out()


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
! Description - test case partrhone
!
! Reference -
!
! Calling program - particle_trajects_init
!
!************************************************************************
!
USE iopars  
USE partvars
USE switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_pout_params'
CALL log_timer_in()

!---file parameters  
part1d(1)%defined = .TRUE.

!---output attributes
outppars(1)%nodim = 2
outppars(1)%tlims(1) = MERGE(12960,720,iopt_hydro_impl.EQ.0)
outppars(1)%tlims(2:3) = (/nstep,1/)
outppars(1)%aging = .TRUE.
outppars(1)%ntype = 2

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
! Description - test case partrhone
!
! Reference -
!
! Calling program - define_parcld_data
!
! Module calls -
!  
!************************************************************************
!
USE iopars
USE grid
USE partpars
USE partvars
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

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
INTEGER :: ipos, jpos
INTEGER, SAVE, DIMENSION(3) :: numrel


procname(pglev+1) = 'usrdef_parcld_data'
CALL log_timer_in()

IF (modfiles(io_parcld,ifil,1)%iostat.EQ.0) THEN
   modfiles(io_parcld,ifil,1)%iostat = 1
   numrel(ifil) = 0
   GOTO 1000
ENDIF

numrel(ifil) = numrel(ifil) + 1

SELECT CASE (ifil)
   CASE (1); ipos = 59; jpos = 24
   CASE (2); ipos = 51; jpos = 22
   CASE (3); ipos = 67; jpos = 22
END SELECT
disdata(1) = 0.25*(gxcoord(ipos,jpos)+gxcoord(ipos+1,jpos)+&
                 & gxcoord(ipos,jpos+1)+gxcoord(ipos+1,jpos+1))
disdata(2) = 0.25*(gycoord(ipos,jpos)+gycoord(ipos+1,jpos)+&
                 & gycoord(ipos,jpos+1)+gycoord(ipos+1,jpos+1))
disdata(3) = 1.0

SELECT CASE (numrel(ifil))
   CASE (1); ciodatetime(1:19) = '2003/01/02;12:00:00' 
   CASE (2); ciodatetime(1:19) = '2003/01/03;12:00:00' 
   CASE (3); ciodatetime(1:19) = '2003/01/04;12:00:00'
END SELECT
ciodatetime(20:23) = ':000'

IF (numrel(ifil).EQ.cloud(ifil)%noreleases) THEN
   modfiles(io_parcld,ifil,1)%iostat = 2
ENDIF

1000 CALL log_timer_out()  


RETURN

END SUBROUTINE usrdef_parcld_data
