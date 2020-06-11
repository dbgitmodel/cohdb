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
! *Usrdef_Structures* User-defined model setup for structures
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Structures.f90  V2.6
!
! $Date: 2014-03-21 16:29:34 +0100 (Fri, 21 Mar 2014) $
!
! $Revision: 684 $
!
! Description - test case drythin
!
! Reference -
!
! Subroutines - usrdef_dry_cells, usrdef_thin_dams, usrdef_weirs,
!               usrdef_dischr_spec, usrdef_dischr_data
!
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_dry_cells
!************************************************************************
!
! *usrdef_dry_cells* Define locations of dry cells
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Structures.f90  V2.6
!
! Description - test case drythin
!
! Reference -
!
! Calling program - define_global_grid
!
!************************************************************************
!
USE iopars
USE structures
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
INTEGER :: l


procname(pglev+1) = 'usrdef_dry_cells'
CALL log_timer_in()

l_110: DO l=1,numdry
   idry(l)  = 41 + MOD(l-1,20)
   jdry(l) = (l-1)/20 + 16
ENDDO l_110

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_dry_cells

!========================================================================

SUBROUTINE usrdef_thin_dams
!************************************************************************
!
! *usrdef_thin_dams* Define locations of thin dams
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Structures.f90  V2.6
!
! Description - test case drythin
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE structures
USE time_routines, ONLY: log_timer_in, log_timer_out


IMPLICIT NONE

INTEGER :: ii, jj

procname(pglev+1) = 'usrdef_thin_dams'
CALL log_timer_in()

!
!1. Thin dams at U-nodes
!-----------------------
!

ii_110: DO ii=1,numthinu
   ithinu(ii) = 10*((ii-1)/10) + 21
   jthinu(ii) = MOD(ii-1,10) + 1
ENDDO ii_110

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_thin_dams

!========================================================================

SUBROUTINE usrdef_weirs
!************************************************************************
!
! *usrdef_weirs* Define locations and arrays for weirs and barriers
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Structures.f90  V2.6
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_weirs


!========================================================================

SUBROUTINE usrdef_dischr_spec
!************************************************************************
!
! *usrdef_discharges* Define specifier arrays for discharge module
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Structures.f90  V2.6
!
! Description -
!
! Reference -
!
! Calling program - define_dischr_spec
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_dischr_spec

!========================================================================

SUBROUTINE usrdef_dischr_data(iddesc,ifil,ciodatetime,disdata,nodat,novars)
!************************************************************************
!
! *usrdef_discharges* Define discharge data
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Structures.f90  V2.6
!
! Description -
!
! Reference - define_dischr_data
!
! Calling program -
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(OUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, nodat, novars
REAL, INTENT(INOUT), DIMENSION(nodat,novars) :: disdata

!
! Name         Type Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*disdata*     REAL    Input data
!*nodat*       INTEGER Number of discharge locations
!*novars*      INTEGER Number of input variables
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_dischr_data
