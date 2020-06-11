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
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - output parameters for test case rhone
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - usrdef_tsr0d_vals
!
! Module calls - close_file, mult_index, open_file
!
!************************************************************************
!
USE iopars
USE paralpars
USE switches
USE timepars
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER, SAVE :: icout, iunit
INTEGER :: ihour
REAL, DIMENSION(11) :: out0ddat


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
   modfiles(io_2uvobc,1:2,1)%status = 'R'
   modfiles(io_2uvobc,1:2,1)%form = modfiles(io_2uvobc,1:2,2)%form
   modfiles(io_2uvobc,1,1)%filename = modfiles(io_2uvobc,1,2)%filename
   modfiles(io_2uvobc,2,1)%filename = modfiles(io_2uvobc,2,2)%filename
   modfiles(io_2uvobc,1:2,2)%status = '0'
!  ---open boundary conditions (3-D)
   modfiles(io_3uvobc,1:2,1)%status = 'R'
   modfiles(io_3uvobc,1:2,1)%form = modfiles(io_3uvobc,1:2,2)%form
   modfiles(io_3uvobc,1,1)%filename = modfiles(io_3uvobc,1,2)%filename
   modfiles(io_3uvobc,2,1)%filename = modfiles(io_3uvobc,2,2)%filename
   modfiles(io_3uvobc,1:2,2)%status = '0'
!  ---open boundary conditions (salinity)
   modfiles(io_salobc,1:2,1)%status = 'R'
   modfiles(io_salobc,1:2,1)%form = modfiles(io_salobc,1:2,2)%form
   modfiles(io_salobc,1,1)%filename = modfiles(io_salobc,1,2)%filename
   modfiles(io_salobc,2,1)%filename = modfiles(io_salobc,2,2)%filename
   modfiles(io_salobc,1:2,2)%status = '0'
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
      WRITE (iunit,'(A)') 'Output parameters for test case rhone: simulation '&
                        & //TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF
!  ---output frequency
   icout = MERGE(2160,120,iopt_hydro_impl.EQ.0)
ENDIF

!
!3. Evaluate/write parameters
!----------------------------
!

IF ((nt.GT.0).AND.(mult_index(nt,icout))) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

   ihour = nosecsrun/3600

!  ---evaluate
   CALL usrdef_tsr0d_vals(out0ddat,11)

!  ---write
   IF (master) THEN
      WRITE (iunit,9001) ihour
      WRITE (iunit,9002) 'area34', out0ddat(1)
      WRITE (iunit,9002) 'area22', out0ddat(2)
      WRITE (iunit,9002) 'ekin', out0ddat(3)
      WRITE (iunit,9002) 'epot', out0ddat(4)
      WRITE (iunit,9002) 'etot', out0ddat(5)
      WRITE (iunit,9002) 'bdissip', out0ddat(6)
      WRITE (iunit,9002) 'enstr', out0ddat(7)
      WRITE (iunit,9002) 'curmax', out0ddat(8)
      WRITE (iunit,9002) 'wmax', out0ddat(9)
      WRITE (iunit,9002) 'wmin', out0ddat(10)
      WRITE (iunit,9002) 'd34max', out0ddat(11)
   ENDIF

   CALL log_timer_out()

ENDIF

!
!4. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')


RETURN

9001 FORMAT ('Time: ',I2,' hours')
9002 FORMAT (T3,A,T12,': ',G12.4E3)

END SUBROUTINE usrdef_output
