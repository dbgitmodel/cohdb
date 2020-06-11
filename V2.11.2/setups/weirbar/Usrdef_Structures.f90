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
! $Date: 2014-09-19 18:01:49 +0200 (Fri, 19 Sep 2014) $
!
! $Revision: 741 $
!
! Description - test case weirbar
!
! Reference -
!
! Subroutines - usrdef_dry_cells, usrdef_thin_dams, usrdef_weirs,
!               usrdef_dischr_spec, usrdef_dischr_data
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
! Description - 
!
! Reference -
!
! Calling program - define_global_grid
!
!************************************************************************
!

IMPLICIT NONE


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
! Description - test case weirbar
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

!
!*Local variables
!
INTEGER :: ii

procname(pglev+1) = 'usrdef_weirs'
CALL log_timer_in()

!
!1. Weir coordinates
!-------------------
!
!1.1 Weirs at U-nodes
!--------------------
!

IF (LLT(runtitle(8:8),'E')) THEN
   iwbaru(1) = 101
   jwbaru(1) = 1
ELSE
   ii_110: DO ii=1,numwbaru
      iwbaru(ii) = 10*((ii-1)/5) + 71
      jwbaru(ii) = MOD(ii-1,5) + 1
   ENDDO ii_110
ENDIF

!
!2. Weir parameters
!------------------
!
!2.1 Calibration coefficients
!----------------------------
!
!---friction coefficient (weirs)
wbarcoefu = 0.9
!---friction coefficient (orifices)
oricoefu = 0.04
!---modular limit
wbarmodlu = 0.7
!---relaxation parameter
wbarrlxu = 0.5

!
!2.2 Crest/orifice heights
!-------------------------
!

SELECT CASE (runtitle(8:8))
!  ---2-D weir
   CASE ('A') 
      wbarcrestu(1) = 9.0
      oriheightu(1) = 0.0              
      orisillu(1) = 0.0
!  ---3-D weir
   CASE ('B')
      wbarcrestu(1) = 8.0
      oriheightu(1) = 0.0              
      orisillu(1) = 0.0
!  ---3-D weir + orifice
   CASE ('C')
      wbarcrestu(1) = 8.0
      oriheightu(1) = 2.0              
      orisillu(1) = 1.0   
!  ---surface orifice 
   CASE ('D')
      wbarcrestu(1) = 15.0
      oriheightu(1) = 1.0              
      orisillu(1) = 8.5   
!  ---series of aligned 3-D weirs
   CASE ('E')
      wbarcrestu(1:5)   = 9.0
      wbarcrestu(6:10)  = 8.0
      wbarcrestu(11:15) = 7.0
      wbarcrestu(16:20) = 6.0
      wbarcrestu(21:25) = 7.0
      wbarcrestu(26:30) = 8.0
      wbarcrestu(31:35) = 9.0
      oriheightu = 0.0
      orisillu = 0.0
END SELECT

CALL log_timer_out()


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
