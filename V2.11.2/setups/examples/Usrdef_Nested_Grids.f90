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
!                       (example routines)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Nested_Grids.f90  V2.0
!
! $Date: 2013-04-16 12:20:30 +0200 (Tue, 16 Apr 2013) $
!
! $Revision: 548 $
!
! Description -
!
! Reference -
!
! Routines - usrdef_nstgrd_abs, usrdef_nstgrd_rel, usrdef_nstgrd_spec
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_nstgrd_spec
!************************************************************************
!
! *usrdef_nstgrd_spec* Define number of nested open boundaries and specifier
!                      arrays (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Nested_Grids.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_nstgrd_spec
!
! External calls -
!
! Module calls -
!
!************************************************************************
!
USE iopars
USE nestgrids
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_nstgrd_spec'
CALL log_timer_in()

!---type of grids
nestcoords(1:nonestsets) = ?

!---number of horizontal nest points
nohnstglbc(1:nonestnests) = ?
nohnstglbu(1:nonestnests) = ?
nohnstglbv(1:nonestnests) = ?

!---number of vertical nest points
novnst(1:nonestsets) = ?

!---type of 2-D output
inst2dtype(1:nonestsets) = ?

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_nstgrd_spec

!========================================================================

SUBROUTINE usrdef_nstgrd_abs(iset,nhdat,nzdat,xcoord,ycoord,zcoord,cnode)
!************************************************************************
!
! *usrdef_nstgrd_abs* Define absolute geographical positions of nested open
!                     boundaries (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Nested_Grids.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_nstgrd_locs
!
! External calls -
!
! Module calls -
!
!************************************************************************
!
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=1), INTENT(IN) :: cnode
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


procname(pglev+1) = 'usrdef_nstgrd_abs'
CALL log_timer_in()

xcoord(1:nhdat) = ?
ycoord(1:nhdat) = ?
zcoord(1:nhdat,1:nzdat) = ?

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_nstgrd_abs

!========================================================================

SUBROUTINE usrdef_nstgrd_rel(iset,nhdat,nzdat,nonodes,hnests,zcoord,cnode)
!************************************************************************
!
! *usrdef_nstgrd_real* Define relative geographical positions of nested open
!                      boundaries (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Nested_Grids.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_nstgrd_locs
!
! External calls -
!
! Module calls -
!
!************************************************************************
!
USE datatypes
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=1), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: iset, nhdat, nonodes, nzdat
REAL, INTENT(OUT), DIMENSION(nhdat,nzdat) :: zcoord
TYPE (HRelativeCoords), INTENT(OUT), DIMENSION(nhdat,nonodes) :: hnests

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
!*Local variables
!
INTEGER :: inode


procname(pglev+1) = 'usrdef_nstgrd_rel'
CALL log_timer_in()

inode_110: DO inode=1,nonodes
   hnests(:,inode)%icoord = ?
   hnests(:,inode)%jcoord = ?
   hnests(:,inode)%xcoord = ?
   hnests(:,inode)%ycoord = ?
ENDDO inode_110
zcoord(1:nhdat,1:nzdat) = ?

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_nstgrd_rel
