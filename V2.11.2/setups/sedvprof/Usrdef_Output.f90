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
! Version - @(COHERENS)Usrdef_Output.f90  V2.5
!
! $Date: 2018-07-23 10:21:02 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1167 $
!
! Description - output parameters for test case sedvprof
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, mult_indesx, open_file
!
!************************************************************************
!
USE currents
USE grid
USE gridpars
USE iopars
USE sedarrays
USE syspars
USE timepars
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER, SAVE :: iunit
REAL :: rhour


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
   modfiles(io_inicon,ics_phys,1)%form = &
 & modfiles(io_fincon,ics_phys,2)%form
   modfiles(io_inicon,ics_phys,1)%filename = &
 & modfiles(io_fincon,ics_phys,2)%filename
   modfiles(io_fincon,ics_phys,2)%status = '0'
!  ---surface forcing
   modfiles(io_1uvsur,1:2,1)%status = 'R'
   modfiles(io_1uvsur,1:2,1)%form = &
 & modfiles(io_1uvsur,1:2,2)%form
   modfiles(io_1uvsur,1:2,1)%filename = modfiles(io_1uvsur,1:2,2)%filename
   modfiles(io_1uvsur,1:2,2)%status = '0'
!  ---sediment particle attributes
   modfiles(io_sedspc,1,1)%status = 'R'
   modfiles(io_sedspc,1,1)%form =   modfiles(io_sedspc,1,2)%form
   modfiles(io_sedspc,1,1)%filename = modfiles(io_sedspc,1,2)%filename
   modfiles(io_sedspc,1,2)%status = '0'
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (cold_start) RETURN

IF (nt.EQ.0) THEN
!  ---open output file
   CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
   WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
   WRITE (iunit,'(A)') 'Output parameters for test case sedvprof: '&
                        & //TRIM(runtitle)
   WRITE (iunit,*)
ENDIF

!
!3. Write output parameters
!--------------------------
!

IF (nt.GT.0.AND.mult_index(nt,150)) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

   rhour = nosecsrun/3600.0
   WRITE (iunit,9001) rhour
   WRITE (iunit,9002) 'umean', umvel(1,1)
   WRITE (iunit,9002) 'ubot', uvel(1,1,1)
   WRITE (iunit,9002) 'usur', uvel(1,1,nz)
   WRITE (iunit,9002) 'eqconc', ceq(1,1,1)
   WRITE (iunit,9002) 'heightc', height_c(1,1,1)
   WRITE (iunit,9002) 'bsed_dep', bottom_sed_dep(1,1,1)
   WRITE (iunit,9002) 'bsed_ero', bottom_sed_ero(1,1,1)
   WRITE (iunit,9002) 'bsed_flux', bottom_sed_flux(1,1,1)
   WRITE (iunit,9002) 'sedbot', cvol(1,1,1,1)
   WRITE (iunit,9002) 'sedmin', MINVAL(ctot(1,1,:))
   WRITE (iunit,9002) 'sedmax', MAXVAL(ctot(1,1,:))
   WRITE (iunit,9002) 'sedint', SUM(delzatc(1,1,:)*ctot(1,1,:))

   CALL log_timer_out()

ENDIF

!
!4. Close output file
!--------------------
!

IF (nt.EQ.nstep) CALL close_file(iunit,'A')




RETURN

9001 FORMAT ('Time: ',F4.1,' hours')
9002 FORMAT (T3,A,T19,': ',G12.4E3)

END SUBROUTINE usrdef_output





