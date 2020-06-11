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
! Version - @(COHERENS)Usrdef_Output.f90  V2.7
!
! $Date: 2018-06-05 16:03:05 +0200 (Tue, 05 Jun 2018) $
!
! $Revision: 1142 $
!
! Description - output parameters for test obctang
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
INTEGER, PARAMETER :: nostatsvarsell = 4, nostatsvars2d = 3, noutvars0d = 4
INTEGER, SAVE :: iunit
INTEGER :: ivar, l
INTEGER, DIMENSION(noutvars0d) :: vecids
REAL, DIMENSION(noutvars0d) :: out0dres
REAL, DIMENSION(nostatsanal,nostatsvars2d) :: ampdat, phadat
REAL, DIMENSION(nostatsanal,nostatsvarsell) :: elldat


!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN

!
!1.1 Model grid
!--------------
!

   modfiles(io_modgrd,1,1)%status = 'R'
   modfiles(io_modgrd,1,1)%form = modfiles(io_modgrd,1,2)%form
   modfiles(io_modgrd,1,1)%filename = modfiles(io_modgrd,1,2)%filename
   modfiles(io_modgrd,1,2)%status = '0'

!
!1.2 Open boundary conditions
!----------------------------
!
!  ---2-D normal
   modfiles(io_2uvobc,1,1)%status = 'R'
   modfiles(io_2uvobc,1,1)%form = modfiles(io_2uvobc,1,2)%form
   modfiles(io_2uvobc,1,1)%filename = modfiles(io_2uvobc,1,2)%filename
   modfiles(io_2uvobc,1,2)%status = '0'

!  ---2-D tangential
   IF (iopt_obc_2D_tang.EQ.1) THEN
      modfiles(io_2xyobc,1,1)%status = 'R'
      modfiles(io_2xyobc,1,1)%form = modfiles(io_2xyobc,1,2)%form
      modfiles(io_2xyobc,1,1)%filename = modfiles(io_2xyobc,1,2)%filename
      modfiles(io_2xyobc,1,2)%status = '0'
   ENDIF

!  ---3-D normal
   IF (iopt_obc_3D.EQ.1) THEN
      modfiles(io_3uvobc,1,1)%status = 'R'
      modfiles(io_3uvobc,1,1)%form = modfiles(io_3uvobc,1,2)%form
      modfiles(io_3uvobc,1,1)%filename = modfiles(io_3uvobc,1,2)%filename
      modfiles(io_3uvobc,1,2)%status = '0'
   ENDIF

!  ---3-D tangential
   IF (iopt_obc_3D_tang.EQ.1) THEN
      modfiles(io_3xyobc,1,1)%status = 'R'
      modfiles(io_3xyobc,1,1)%form = modfiles(io_3xyobc,1,2)%form
      modfiles(io_3xyobc,1,1)%filename = modfiles(io_3xyobc,1,2)%filename
      modfiles(io_3xyobc,1,2)%status = '0'
   ENDIF

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

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN

!  ---open output file
   IF (.NOT.cold_start.AND.master) THEN
      CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case obctang:'//&
                        & ' simulation '//TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF

ENDIF

IF (.NOT.master.OR.nt.LT.nstep) RETURN

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
vecids = (/(ivar+1,ivar=1,noutvars0d)/)
CALL read_vars(out0dres,res0d(1),0,vecids=vecids)
CALL close_filepars(res0d(1))

!
!3.2 Station data
!----------------
!
!---amplitudes
amp2d(2,1)%status = 'R'
CALL open_filepars(amp2d(2,1))
CALL read_output_metadata(amp2d(2,1))
READ (amp2d(2,1)%iunit,*)
ivar_321: DO ivar=1,nostatsvars2d
   READ (amp2d(2,1)%iunit,*)
   READ (amp2d(2,1)%iunit,*) ampdat(:,ivar)
ENDDO ivar_321
CALL close_filepars(amp2d(2,1))

!---phases
pha2d(2,1)%status = 'R'
CALL open_filepars(pha2d(2,1))
CALL read_output_metadata(pha2d(2,1))
READ (pha2d(2,1)%iunit,*)
ivar_322: DO ivar=1,nostatsvars2d
   READ (pha2d(2,1)%iunit,*)
   READ (pha2d(2,1)%iunit,*) phadat(:,ivar)
ENDDO ivar_322
CALL close_filepars(pha2d(2,1))

!---elliptic parameters
ell2d(2,1)%status = 'R'
CALL open_filepars(ell2d(2,1))
CALL read_output_metadata(ell2d(2,1))
READ (ell2d(2,1)%iunit,*)
ivar_323: DO ivar=1,nostatsvarsell
   READ (ell2d(2,1)%iunit,*)
   READ (ell2d(2,1)%iunit,*) elldat(:,ivar)
ENDDO ivar_323
CALL close_filepars(ell2d(2,1))

!
!4. Define/write output parameters
!---------------------------------
!
!---global parameters
WRITE (iunit,9001) nosecsrun/86400
WRITE (iunit,'(A)') 'Global parameters'
WRITE (iunit,9002) 'ekin', 1.0E-09*out0dres(1)
WRITE (iunit,9002) 'epot', 1.0E-09*out0dres(2)
WRITE (iunit,9002) 'etot', 1.0E-09*out0dres(3)
WRITE (iunit,9002) 'bdissip', 0.001*out0dres(4)
WRITE (iunit,*)

!---station data
WRITE (iunit,'(A)') 'Harmonic parameters'
l_420: DO l=1,nostatsanal
   WRITE (iunit,9003) 'Station', l
   WRITE (cl,'(I12)') l; cl = ADJUSTL(cl)
   WRITE (iunit,9002) 'zetamp'//TRIM(cl), 100.0*ampdat(l,1) 
   WRITE (iunit,9002) 'zetpha'//TRIM(cl), phadat(l,1)
   WRITE (iunit,9002) 'umamp'//TRIM(cl), 100.0*ampdat(l,2)
   WRITE (iunit,9002) 'umpha'//TRIM(cl), phadat(l,2)
   WRITE (iunit,9002) 'vmamp'//TRIM(cl), 100.0*ampdat(l,3)
   WRITE (iunit,9002) 'vmpha'//TRIM(cl), phadat(l,3)
   WRITE (iunit,9002) 'ellmaj'//TRIM(cl), 100.0*elldat(l,1)
   WRITE (iunit,9002) 'ellipt'//TRIM(cl), elldat(l,2)
   WRITE (iunit,9002) 'ellinc'//TRIM(cl), elldat(l,3)
   WRITE (iunit,9002) 'ellpha'//TRIM(cl), elldat(l,4)
ENDDO l_420

!
!4.3 Close output file
!---------------------
!

CALL close_file(iunit,'A')

CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',I2,' days') 
9002 FORMAT (T3,A,T12,': ',G12.4E3)
9003 FORMAT (A,': ',I3)

END SUBROUTINE usrdef_output
