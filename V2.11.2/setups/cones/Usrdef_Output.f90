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

SUBROUTINE usrdef_output
!************************************************************************
!
! *usrdef_output* User-formatted output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.1.1
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - output parameters for test case cones
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_mod, mult_index, open_file
!
!************************************************************************
!
USE density
USE gridpars
USE iopars
USE modids
USE paralpars
USE syspars
USE timepars
USE inout_routines, ONLY: close_file, open_file
USE paral_comms, ONLY: combine_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index


IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER :: i, icenter, icout = 0, iunit = 0, j, jcenter
REAL :: cmin, cmax, dh = 0.0, mtime, omega = 0.0, phase0 = 0.0, phase, &
      & radc = 0.0, rflag = -999.9, width = 0.0, xcenter, xmin, xplus, &
      & ycenter, ymin, yplus
REAL, DIMENSION(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo,nz) :: salglb


!
!1. Reset forcing attributes for forcing files
!---------------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
!  ---initial conditions
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = modfiles(io_fincon,ics_phys,2)%form
   modfiles(io_inicon,ics_phys,1)%filename = &
 & modfiles(io_fincon,ics_phys,2)%filename
   modfiles(io_inicon,ics_phys,2)%status = '0'
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN

!  ---open output file
   IF (.NOT.cold_start.AND.master) THEN
      CALL open_file(iunit,runtitle(1:6)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case cones: simulation '&
                         & //TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF

!  ---grid spacing
   dh = 1.0

!  ---basin width
   width = nc-1

!  ---initial phase
   phase0 = ATAN2(dh,0.25*width-dh) + pi

!  ---distance of cone center to basin center
   radc = dh*SQRT(82.0)

!  ---angular frequency
   omega = 0.05/60.0

!  ---output frequency
   icout = 50

ENDIF

!
!3. Store parameters
!-------------------
!

IF ((nt.GT.0).AND.(mult_index(nt,icout))) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

   mtime = nosecsrun/60.0

!
!3.1 Combine salinity at master
!------------------------------
!

   CALL combine_mod(salglb,sal,(/1-nhalo,1-nhalo,1/),iarr_sal,0.0)
   IF (.NOT.master) GOTO 1000

!
!3.2 Position of cone center
!---------------------------
!

   phase = (omega*nt)*delt2d - phase0
   xcenter = radc*COS(phase) + 0.5*(width-dh)
   ycenter = radc*SIN(phase) + 0.5*(width-dh)
   icenter = xcenter/dh + 1
   jcenter = ycenter/dh + 1

!
!3.3 Distribution radii in X-direction
!-------------------------------------
!
!  ---left
   flag = .FALSE.
   i_331: DO i=icenter,1,-1
      IF((.NOT.flag).AND.(salglb(i,jcenter,2).LT.0.01)) THEN
         flag = .TRUE.
         xmin = xcenter - dh*&
              & (i+(0.01-0.5*(salglb(i+1,jcenter,2)+salglb(i,jcenter,2))/&
              & (salglb(i+1,jcenter,2)-salglb(i,jcenter,2))))

      ENDIF
   ENDDO i_331
   IF (.NOT.flag) xmin = rflag
 
!  ---right
   flag = .FALSE.
   i_332: DO i=icenter,nc-1
      IF((.NOT.flag).AND.(salglb(i,jcenter,2).LT.0.01)) THEN
         flag = .TRUE.
         xplus = dh*&
               & (i+(0.01+(0.5*salglb(i-1,jcenter,2)-1.5*salglb(i,jcenter,2))/&
               & (salglb(i,jcenter,2)-salglb(i-1,jcenter,2)))) - xcenter
      ENDIF
   ENDDO i_332
   IF (.NOT.flag) xplus = rflag

!
!3.4 Distribution radii in Y-direction
!-------------------------------------
!
!  ---left
   flag = .FALSE.
   j_341: DO j=jcenter,1,-1
      IF((.NOT.flag).AND.(salglb(icenter,j,2).LT.0.01)) THEN
         flag = .TRUE.
         ymin = ycenter - dh*&
              & (j+(0.01-0.5*(salglb(icenter,j+1,2)+salglb(icenter,j,2))/&
              & (salglb(icenter,j+1,2)-salglb(icenter,j,2))))
      ENDIF
   ENDDO j_341
   IF (.NOT.flag) ymin = rflag

!  ---right
   flag = .FALSE.
   j_342: DO j=jcenter,nr-1
      IF((.NOT.flag).AND.(salglb(icenter,j,2).LT.0.01)) THEN
         flag = .TRUE.
         yplus = dh*&
               & (j+(0.01+(0.5*salglb(icenter,j-1,2)-1.5*salglb(icenter,j,2))/&
               & (salglb(icenter,j,2)-salglb(icenter,j-1,2)))) - ycenter
      ENDIF
   ENDDO j_342
   IF (.NOT.flag) yplus = rflag

!
!3.5 Max/min concentrations
!--------------------------
!

   cmin = MINVAL(salglb(1:nc-1,1:nr-1,2))
   cmax = MAXVAL(salglb(1:nc-1,1:nr-1,2))

!
!3.6 Write output
!----------------
!

   WRITE (iunit,9001) mtime
   WRITE (iunit,9002) 'xmin', xmin
   WRITE (iunit,9002) 'xplus', xplus
   WRITE (iunit,9002) 'ymin', ymin
   WRITE (iunit,9002) 'yplus', yplus
   WRITE (iunit,9002) 'cmin', cmin
   WRITE (iunit,9002) 'cmax', cmax

1000 CALL log_timer_out()

ENDIF

!
!4. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')


RETURN

9001 FORMAT ('Time: ',F5.1,' minutes')
9002 FORMAT (T3,A,T12,': ',G12.4E3)

END SUBROUTINE usrdef_output
