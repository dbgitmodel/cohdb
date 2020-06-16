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
! Version - @(COHERENS)Usrdef_Output.f90  V2.7.1
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - output parameters for test case bohai
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, close_filepars, combine_mod, error_alloc,
!                open_file, open_filepars, read_output_metadata, read_submod,
!                read_vars
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE paralpars
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: close_file, close_filepars, open_file, &
                        & open_filepars, read_output_metadata, read_submod, &
                        & read_vars
USE paral_comms, ONLY: combine_mod
USE time_routines, ONLY: log_timer_in, log_timer_out


IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: packing
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
CHARACTER (LEN=12) :: cl
INTEGER, PARAMETER :: nostatsvarsell = 4, nostatsvars2d = 3, noutvars0d = 4, &
                    & noutvars2d = 7
INTEGER, SAVE :: iunit
INTEGER :: i, ivar, j, l, n, nocoords2d, noutmax
INTEGER, ALLOCATABLE, DIMENSION(:) :: vecids
REAL :: delx, dely, x0, y0
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: maskglb
REAL, DIMENSION(noutvars0d) :: out0dres
REAL, DIMENSION(2,2) :: xpos, ypos, zetmin
REAL, DIMENSION(nostatsanal,nostatsvars2d,2) :: ampdat, phadat
REAL, DIMENSION(nostatsanal,nostatsvarsell,2) :: elldat
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: out2damp


!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,1)%status = 'R'
   modfiles(io_modgrd,1,1)%form = modfiles(io_modgrd,1,2)%form 
   modfiles(io_modgrd,1,1)%filename = modfiles(io_modgrd,1,2)%filename
   modfiles(io_modgrd,1,2)%status = '0'
!  ---open boundary conditions (2-D)
   modfiles(io_2uvobc,1,1)%status = 'R'
   modfiles(io_2uvobc,1,1)%form = modfiles(io_2uvobc,1,2)%form 
   modfiles(io_2uvobc,1,1)%filename = modfiles(io_2uvobc,1,2)%filename
   modfiles(io_2uvobc,1,2)%status = '0'
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN

!  ---open output file
   IF (.NOT.cold_start.AND.master) THEN
      CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case bohai: simulation '&
                        & //TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF

ENDIF

IF (nt.LT.nstep) RETURN

procname(pglev+1) = 'usrdef_output'
CALL log_timer_in()

!
!3. Allocate/initialise
!----------------------
!
!---allocate arrays
ALLOCATE (maskglb(nc,nr),STAT=errstat)
CALL error_alloc('maskglb',2,(/nc,nr/),log_type)
ALLOCATE (out2damp(nc-1,nr-1,noutvars2d,2),STAT=errstat)
CALL error_alloc('out2damp',4,(/nc-1,nr-1,noutvars2d,2/),real_type)
noutmax = MAX(noutvars0d,noutvars2d)
ALLOCATE (vecids(noutmax),STAT=errstat)
CALL error_alloc('vecids',1,(/noutmax/),int_type)

!---global mask arrays
packing = analgpars(1)%packing
CALL combine_mod(maskglb,maskatc_int,(/1,1/),0,.FALSE.)

!---initialise parameters
x0 = 117.5; y0 = 37.0
delx = 1.0/12.0; dely = 1.0/12.0
nocoords2d = MERGE(5,4,analgpars(1)%packing)

!---return except master
IF (.NOT.master) GOTO 1001

!
!4. Read data
!------------
!
!4.1 Global residuals
!--------------------
!

res0d(1)%status = 'R'
CALL open_filepars(res0d(1))
CALL read_output_metadata(res0d(1))
READ (res0d(1)%iunit,*)
vecids(1:noutvars0d) = (/(ivar+1,ivar=1,noutvars0d)/)
CALL read_vars(out0dres,res0d(1),0,vecids=vecids(1:noutvars0d))
CALL close_filepars(res0d(1))

!
!4.2 Global amplitudes
!---------------------
!

n_420: DO n=1,2
   amp2d(1,n)%status = 'R'
   CALL open_filepars(amp2d(1,n))
   SELECT CASE (amp2d(1,n)%form)
      CASE ('A')
         CALL read_output_metadata(amp2d(1,n))
         READ (amp2d(1,n)%iunit,*)
      CASE ('U')
         CALL read_output_metadata(amp2d(1,n))
         READ (amp2d(1,n)%iunit)
   END SELECT
   vecids(1:noutvars2d) = (/(nocoords2d+ivar,ivar=1,noutvars2d)/)
   CALL read_submod(out2damp(:,:,:,n),amp2d(1,n),0,&
                  & vecids=vecids(1:noutvars2d),&
                  & maskvals=maskglb(1:nc-1,1:nr-1))
   CALL close_filepars(amp2d(1,n))
ENDDO n_420

!
!4.3 Station data
!----------------
!

n_430: DO n=1,2

!  ---amplitudes
   amp2d(2,n)%status = 'R'
   CALL open_filepars(amp2d(2,n))
   CALL read_output_metadata(amp2d(2,n))
   READ (amp2d(2,n)%iunit,*)
   ivar_431: DO ivar=1,nostatsvars2d
      READ (amp2d(2,n)%iunit,*)
      READ (amp2d(2,n)%iunit,*) ampdat(:,ivar,n)
   ENDDO ivar_431
   CALL close_filepars(amp2d(2,n))

!  ---phases
   pha2d(2,n)%status = 'R'
   CALL open_filepars(pha2d(2,n))
   CALL read_output_metadata(pha2d(2,n))
   READ (pha2d(2,n)%iunit,*)
   ivar_432: DO ivar=1,nostatsvars2d
      READ (pha2d(2,n)%iunit,*)
      READ (pha2d(2,n)%iunit,*) phadat(:,ivar,n)
   ENDDO ivar_432
   CALL close_filepars(pha2d(2,n))

!  ---elliptic parameters
   ell2d(2,n)%status = 'R'
   CALL open_filepars(ell2d(2,n))
   CALL read_output_metadata(ell2d(2,n))
   READ (ell2d(2,n)%iunit,*)
   ivar_433: DO ivar=1,nostatsvarsell
      READ (ell2d(2,n)%iunit,*)
      READ (ell2d(2,n)%iunit,*) elldat(:,ivar,n)
   ENDDO ivar_433
   CALL close_filepars(ell2d(2,n))

ENDDO n_430

!
!5. Define/write output parameters
!---------------------------------
!
!5.1 Location of amphidromic points
!----------------------------------
!

n_510: DO n=1,2
   zetmin(1,n) = MAXVAL(out2damp(:,:,5,n))
   zetmin(2,n) = zetmin(1,n)
   i_511: DO i=1,nc-1
   j_511: DO j=1,nr-1
      IF (maskglb(i,j).AND.(i.GT.25.AND.i.LT.43).AND.&
       & (j.GT.1.AND.j.LT.19)) THEN
         IF (out2damp(i,j,5,n).LT.zetmin(1,n)) THEN
            zetmin(1,n) = out2damp(i,j,5,n)
            xpos(1,n) = x0 + (i-0.5)*delx
            ypos(1,n) = y0 + (j-0.5)*dely
         ENDIF
      ENDIF
      IF (maskglb(i,j).AND.(i.GT.61.AND.i.LT.79).AND.&
       & (j.GT.1.AND.j.LT.19)) THEN
         IF (out2damp(i,j,5,n).LT.zetmin(2,n)) THEN
            zetmin(2,n) = out2damp(i,j,5,n)
            xpos(2,n) = x0 + (i-0.5)*delx
            ypos(2,n) = y0 + (j-0.5)*dely
         ENDIF
      ENDIF
   ENDDO j_511
   ENDDO i_511
ENDDO n_510
zetmin = 100.0*zetmin

!
!5.2 Write output parameters
!---------------------------
!
!---global parameters
WRITE (iunit,9001) nosecsrun/86400
WRITE (iunit,'(A)') 'Global parameters'
WRITE (iunit,9002) 'ekin', out0dres(1)
WRITE (iunit,9002) 'epot', out0dres(2)
WRITE (iunit,9002) 'etot', out0dres(3)
WRITE (iunit,9002) 'bdissip', out0dres(4)
WRITE (iunit,*)

!---amphidromic points
n_521: DO n=1,2
   IF (n.EQ.1) THEN
      WRITE (iunit,'(A)') 'Location of M2-amphidromic points'
   ELSE
      WRITE (iunit,'(A)') 'Location of S2-amphidromic points'
   ENDIF
   WRITE (iunit,9002) 'zetmin', zetmin(1,n)
   WRITE (iunit,9003) 'xpos', xpos(1,n)
   WRITE (iunit,9003) 'ypos', ypos(1,n)
   WRITE (iunit,9002) 'zetmin', zetmin(2,n)
   WRITE (iunit,9003) 'xpos', xpos(2,n)
   WRITE (iunit,9003) 'ypos', ypos(2,n)
ENDDO n_521
WRITE (iunit,*)

!---station data
WRITE (iunit,'(A)') 'Harmonic parameters'
l_522: DO l=1,nostatsanal
   WRITE (iunit,9004) 'Station', l
   WRITE (cl,'(I12)') l; cl = ADJUSTL(cl)
   WRITE (iunit,9002) 'M2zetamp'//TRIM(cl), 100.0*ampdat(l,3,1) 
   WRITE (iunit,9002) 'M2zetpha'//TRIM(cl), phadat(l,3,1)
   WRITE (iunit,9002) 'M2ellmaj'//TRIM(cl), 100.0*elldat(l,1,1)
   WRITE (iunit,9002) 'M2ellipt'//TRIM(cl), elldat(l,2,1)
   WRITE (iunit,9002) 'M2ellinc'//TRIM(cl), elldat(l,3,1)
   WRITE (iunit,9002) 'M2ellpha'//TRIM(cl), elldat(l,4,1)
   WRITE (iunit,9002) 'S2zetamp'//TRIM(cl), 100.0*ampdat(l,3,2) 
   WRITE (iunit,9002) 'S2zetpha'//TRIM(cl), phadat(l,3,2)
   WRITE (iunit,9002) 'S2ellmaj'//TRIM(cl), 100.0*elldat(l,1,2)
   WRITE (iunit,9002) 'S2ellipt'//TRIM(cl), elldat(l,2,2)
   WRITE (iunit,9002) 'S2ellinc'//TRIM(cl), elldat(l,3,2)
   WRITE (iunit,9002) 'S2ellpha'//TRIM(cl), elldat(l,4,2)
ENDDO l_522

!
!5.3 Close output file
!---------------------
!

CALL close_file(iunit,'A')

!
!6. Deallocate arrays
!--------------------
!

1001 CONTINUE
DEALLOCATE (maskglb,out2damp,vecids)

CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',I2,' days') 
9002 FORMAT (T3,A,T12,': ',G12.4E3)
9003 FORMAT (T3,A,T12,': ',F6.2)
9004 FORMAT (A,': ',I3)

END SUBROUTINE usrdef_output
