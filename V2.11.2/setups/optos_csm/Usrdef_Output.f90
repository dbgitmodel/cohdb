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
! Version - @(COHERENS)Usrdef_Output.f90  V2.4
!
! $Date: 2018-07-23 10:21:02 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1167 $
!
! Description - output parameters for test case optos csm
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, close_filepars, open_file, open_filepars,
!                read_output_metadata, read_vars
!
!************************************************************************
!
USE iopars
USE paralpars
USE switches
USE syspars
USE timepars
USE inout_routines, ONLY: close_file, close_filepars, open_file, &
                        & open_filepars, read_output_metadata, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out


IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
CHARACTER (LEN=12) :: cl
CHARACTER (LEN=lenfreq) :: harmname
INTEGER, PARAMETER :: noutvarsell = 4, noutvars0d = 4, noutvars2d = 3
INTEGER, SAVE :: iunit
INTEGER :: ivar, l, n
INTEGER, DIMENSION(noutvars0d) :: vecids
REAL, DIMENSION(noutvars0d) :: out0dres
REAL, DIMENSION(nostatsanal,noutvars2d,nofreqsanal) :: ampdat, phadat
REAL, DIMENSION(nostatsanal,noutvarsell,nofreqsanal) :: elldat


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
!  ---initial conditions
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = modfiles(io_fincon,ics_phys,2)%form
   modfiles(io_inicon,ics_phys,1)%filename = &
 & modfiles(io_fincon,ics_phys,2)%filename
   modfiles(io_fincon,ics_phys,2)%status = '0'
!  ---open boundary conditions (2-D)
   modfiles(io_2uvobc,1,1)%status = 'R'
   modfiles(io_2uvobc,1,1)%form = modfiles(io_2uvobc,1,2)%form
   modfiles(io_2uvobc,1,1)%filename = modfiles(io_2uvobc,1,2)%filename
   modfiles(io_2uvobc,1,2)%status = '0'
!  ---nesting
   IF (iopt_nests.EQ.1) THEN
      modfiles(io_nstspc,1,1)%status = 'R'
      modfiles(io_nstspc,1,1)%form = modfiles(io_nstspc,1,2)%form
      modfiles(io_nstspc,1,1)%filename = modfiles(io_nstspc,1,2)%filename
      modfiles(io_nstspc,1,2)%status = '0'
      modfiles(io_nstrel,1,1)%status = 'R'
      modfiles(io_nstrel,1,1)%form = modfiles(io_nstrel,1,2)%form
      modfiles(io_nstrel,1,1)%filename = modfiles(io_nstrel,1,2)%filename
      modfiles(io_nstrel,1,2)%status = '0'
   ENDIF
ENDIF

IF (.NOT.master.OR.runtitle(7:7).EQ.'0') RETURN

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN

!  ---open output file
   IF (.NOT.cold_start.AND.master) THEN
      CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case optos_csm:'//&
                        & ' simulation '//TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF

ENDIF

IF (nt.LT.nstep) RETURN

procname(pglev+1) = 'usrdef_output'
CALL log_timer_in()

!
!3. Read data
!------------
!
!3.1 Global residuals
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
!3.2 Station data
!----------------
!

n_320: DO n=1,nofreqsanal

!  ---amplitudes
   amp2d(2,n)%status = 'R'
   CALL open_filepars(amp2d(2,n))
   CALL read_output_metadata(amp2d(2,n))
   READ (amp2d(2,n)%iunit,*)
   ivar_321: DO ivar=1,noutvars2d
      READ (amp2d(2,n)%iunit,*)
      READ (amp2d(2,n)%iunit,*) ampdat(:,ivar,n)
   ENDDO ivar_321
   ampdat(:,3,n) = 100.0*ampdat(:,3,n)
   CALL close_filepars(amp2d(2,n))

!  ---phases
   pha2d(2,n)%status = 'R'
   CALL open_filepars(pha2d(2,n))
   CALL read_output_metadata(pha2d(2,n))
   READ (pha2d(2,n)%iunit,*)
   ivar_322: DO ivar=1,noutvars2d
      READ (pha2d(2,n)%iunit,*)
      READ (pha2d(2,n)%iunit,*) phadat(:,ivar,n)
   ENDDO ivar_322
   CALL close_filepars(pha2d(2,n))

!  ---elliptic parameters
   ell2d(2,n)%status = 'R'
   CALL open_filepars(ell2d(2,n))
   CALL read_output_metadata(ell2d(2,n))
   READ (ell2d(2,n)%iunit,*)
   ivar_323: DO ivar=1,noutvarsell
      READ (ell2d(2,n)%iunit,*)
      READ (ell2d(2,n)%iunit,*) elldat(:,ivar,n)
   ENDDO ivar_323
   elldat(:,1,n) = 100.0*elldat(:,1,n)
   CALL close_filepars(ell2d(2,n))

ENDDO n_320

!
!4. Write output parameters
!--------------------------
!

IF (master) THEN

!  ---global parameters
   WRITE (iunit,9001) nosecsrun/86400
   WRITE (iunit,'(A)') 'Global parameters'
   WRITE (iunit,9002) 'ekin', out0dres(1)
   WRITE (iunit,9002) 'epot', out0dres(2)
   WRITE (iunit,9002) 'etot', out0dres(3)
   WRITE (iunit,9002) 'bdissip', out0dres(4)
   WRITE (iunit,*)

!  ---station data
   WRITE (iunit,'(A)') 'Harmonic parameters'
   n_410: DO n=1,nofreqsanal
      harmname = TRIM(harm_freq_names(n))
      l_411: DO l=1,nostatsanal
         WRITE (iunit,9003) TRIM(analstatlocs(l)%name), TRIM(harmname)
         WRITE (cl,'(I12)') l; cl = ADJUSTL(cl)
         WRITE (iunit,9002) TRIM(harmname)//'zetamp'//TRIM(cl), ampdat(l,3,n) 
         WRITE (iunit,9002) TRIM(harmname)//'zetpha'//TRIM(cl), phadat(l,3,n)
         WRITE (iunit,9002) TRIM(harmname)//'ellmaj'//TRIM(cl), elldat(l,1,n)
         WRITE (iunit,9002) TRIM(harmname)//'ellipt'//TRIM(cl), elldat(l,2,n)
         WRITE (iunit,9002) TRIM(harmname)//'ellinc'//TRIM(cl), elldat(l,3,n)
         WRITE (iunit,9002) TRIM(harmname)//'ellpha'//TRIM(cl), elldat(l,4,n)
      ENDDO l_411
   ENDDO n_410
   
!  ---close output file
   CALL close_file(iunit,'A')

ENDIF

CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',I2,' days')
9002 FORMAT (T3,A,T13,': ',G12.4E3)
9003 FORMAT (A,T22,A)

END SUBROUTINE usrdef_output
