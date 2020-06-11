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
! Description - output parameters for test case partrot
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls -
!
!************************************************************************
!
USE iopars
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out


IMPLICIT NONE


!
!1. Reset forcing attributes for forcing files
!---------------------------------------------
!

IF (nt.GT.0) RETURN

IF (ciffile%status.EQ.'W') THEN
   
   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

!  ---initial conditions (physcis)
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = modfiles(io_fincon,ics_phys,2)%form
   modfiles(io_inicon,ics_phys,1)%filename = &
 & modfiles(io_fincon,ics_phys,2)%filename
   modfiles(io_inicon,ics_phys,2)%status = '0'
!  ---initial conditions (particles)
   modfiles(io_inicon,ics_part,1)%status = 'R'
   modfiles(io_inicon,ics_part,1)%form = modfiles(io_fincon,ics_part,2)%form
   modfiles(io_inicon,ics_part,1)%filename = &
 & modfiles(io_fincon,ics_part,2)%filename
   modfiles(io_inicon,ics_part,2)%status = '0'
   
   CALL log_timer_out()

ENDIF


RETURN

END SUBROUTINE usrdef_output
