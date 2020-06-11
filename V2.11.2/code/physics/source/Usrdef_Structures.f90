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
! *Usrdef_Structures* User-defined model setup for structures and discharges
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Structures.f90  V2.6
!
! $Date: 2015-04-10 11:30:34 +0200 (Fri, 10 Apr 2015) $
!
! $Revision: 841 $
!
! Description - routines are empty but need to be declared
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
! Description - empty default file
!
! Reference -
!
! Calling program - define_global_grid
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_dry_cells

!========================================================================

SUBROUTINE usrdef_thin_dams
!************************************************************************
!
! *usrdef_weirs* Define locations of thin dams
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Structures.f90  V2.6
!
! Description - empty default file
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!


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
! Description - empty default file
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!


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
! Description - empty default file
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
! Description - empty default file
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
