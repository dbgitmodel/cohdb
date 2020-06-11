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
! Version - @(COHERENS)Usrdef_Output.f90  V2.11
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - output parameters for test case obcest
!
! Reference -
!
! Calling program - coherens_main
!
! External calls -
!
! Module calls -
!
!************************************************************************
!
USE iopars
USE syspars  
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER, SAVE :: icout, iunit
INTEGER :: ihour
REAL, DIMENSION(15) :: out0ddat


!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.GT.0) RETURN

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
      
   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

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
      
   CALL log_timer_out()
   
ENDIF


RETURN

END SUBROUTINE usrdef_output
