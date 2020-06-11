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
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.11.1
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - output parameters for test case wcprof
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, open_file
!
!************************************************************************
!
USE currents
USE fluxes
USE grid
USE iopars
USE physpars
USE syspars
USE timepars
USE wavevars
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER, SAVE :: iunit
REAL :: rhour, xvar


!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
!  ---surface forcing
   modfiles(io_1uvsur,1:2,1)%status = 'R'
   modfiles(io_1uvsur,1:2,1)%form = modfiles(io_1uvsur,1:2,2)%form
   modfiles(io_1uvsur,1:2,1)%filename = modfiles(io_1uvsur,1:2,2)%filename
   modfiles(io_1uvsur,1:2,2)%status = '0'
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (cold_start) RETURN

IF (nt.EQ.0.AND.runtitle(7:7).EQ.'1') THEN
!  ---open output file
   CALL open_file(iunit,runtitle(1:6)//runtitle(8:8)//'.'//suffix,'OUT','A')
   WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
   WRITE (iunit,'(A)') 'Output parameters for test case wcprof: simulation '&
                     & //runtitle(1:6)//runtitle(8:8)
   WRITE (iunit,*)

ENDIF

!
!3. Write output parameters
!--------------------------
!

IF (nt.EQ.nstep) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

!  ---write output
   rhour = nosecsrun/3600.0
   WRITE (iunit,9001) rhour
   WRITE (iunit,9002) 'ubot', uvel(2,1,1)
   WRITE (iunit,9002) 'umbot', umvel(2,1)
   WRITE (iunit,9002) 'wavevel', wavevel(1,1)
   WRITE (iunit,9002) 'waveheight', waveheight(1,1)
   WRITE (iunit,9002) 'waveperiod', waveperiod(1,1)
   WRITE (iunit,9002) 'waveangle', radtodeg*wavedir(1,1)
   xvar = MAX(0.5*delzatc(1,1,1)/zroughatc(1,1),zbtoz0lim)
   WRITE (iunit,9002) 'bstres_cur', &
                     & density_ref*(bcuratc(1,1)*ckar/LOG(xvar))**2
   WRITE (iunit,9002) 'bstres_mean', density_ref*bstresatc(1,1)
   WRITE (iunit,9002) 'bstres_max', density_ref*bstresatc_max(1,1)
   WRITE (iunit,9002) 'bstres_wav', density_ref*bstresatc_wav(1,1)
   WRITE (iunit,9002) 'zrough', zroughatc(1,1)
   WRITE (iunit,9002) 'zarough', zaroughatc(1,1)
   WRITE (iunit,9002) 'bdragcoef', bdragcoefatc(1,1)
   WRITE (iunit,9002) 'fwave', fwave(1,1)
   WRITE (iunit,9002) 'wavethick', 100.0*wavethickatc(1,1)

!  ---close output file
   IF (runtitle(7:7).EQ.'5') CALL close_file(iunit,'A')

   CALL log_timer_out()

ENDIF


RETURN

9001 FORMAT ('Time: ',F4.1,' hours')
9002 FORMAT (T3,A,T15,': ',G12.4E3)

END SUBROUTINE usrdef_output
