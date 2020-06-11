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
! *Usrdef_Nested_Grids* User-defined setup for grid nesting
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Nested_Grids.f90  V2.0
!
! $Date: 2015-04-10 11:30:34 +0200 (Fri, 10 Apr 2015) $
!
! $Revision: 841 $
!
! Description - routines are empty but need to be declared
!
! Reference -
!
! Routines - usrdef_nstgrd_spec, usrdef_nstgrd_abs, usrdef_nstgrd_rel
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_nstgrd_spec
!************************************************************************
!
! *usrdef_nstgrd_spec* Define number of nested open boundaries and specifier
!                      arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Nested_Grids.f90  V2.0
!
! Description - empty default routine
!
! Reference -
!
! Calling program - define_nstgrd_spec
!
!************************************************************************
!

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_nstgrd_spec

!========================================================================

SUBROUTINE usrdef_nstgrd_abs(iset,nhdat,nzdat,xcoord,ycoord,zcoord,cnode)
!************************************************************************
!
! *usrdef_nstgrd_abs* Define absolute geographical positions of nested open
!                     boundaries
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Nested_Grids.f90  V2.0
!
! Description - empty default routine
!
! Reference -
!
! Calling program - define_nstgrd_locs
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: iset, nhdat, nzdat
REAL, INTENT(OUT), DIMENSION(nhdat) :: xcoord, ycoord
REAL, INTENT(OUT), DIMENSION(nhdat,nzdat) :: zcoord

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iset*      INTEGER Output file number
!*nhdat*     INTEGER Number of horizontal data locations
!*nzdat*     INTEGER Number of vertical data points
!*xcoord*    REAL    X-coordinates of nest locations
!*ycoord*    REAL    Y-coordinates of nest locations
!*zcoord*    REAL    Z-coordinates of nest locations
!*cnode*     CHAR    Nodal type
!
!------------------------------------------------------------------------------
!


xcoord = 0.0; ycoord = 0.0; zcoord = 0.0


RETURN

END SUBROUTINE usrdef_nstgrd_abs

!========================================================================

SUBROUTINE usrdef_nstgrd_rel(iset,nhdat,nzdat,nonodes,hnests,zcoord,cnode)
!************************************************************************
!
! *usrdef_nstgrd_rel* Define relative geographical positions of nested open
!                     boundaries
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Nested_Grids.f90  V2.0
!
! Description - empty default routine
!
! Reference -
!
! Calling program - define_nstgrd_locs
!
! Module calls - hrelativecoords_init
!
!************************************************************************
!
USE datatypes
USE datatypes_init, ONLY: hrelativecoords_init
USE syspars

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: iset, nhdat, nonodes, nzdat
TYPE (HRelativeCoords), INTENT(OUT), DIMENSION(nhdat,nonodes) :: hnests
REAL, INTENT(OUT), DIMENSION(nhdat,nzdat) :: zcoord

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iset*      INTEGER Output file number
!*nhdat*     INTEGER Number of horizontal data locations
!*nzdat*     INTEGER Number of vertical data points
!*nonodes*   INTEGER Number of nodes for interpolation
!*hnests*    DERIVED Relative coordinates of nest locations
!*zcoord*    REAL    Z-coordinates of nest locations
!*cnode*     CHAR    Nodal type
!
!------------------------------------------------------------------------------
!


CALL hrelativecoords_init(hnests,.FALSE.)
zcoord = 0.0


RETURN

END SUBROUTINE usrdef_nstgrd_rel
