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
! Version - @(COHERENS)Usrdef_Nested_Grids.f90  V2.9
!
! $Date: 2015-10-07 16:02:06 +0200 (Wed, 07 Oct 2015) $
!
! $Revision: 889 $
!
! Description - optos csm
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
!                      arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Nested_Grids.f90  V2.0
!
! Description - optos csm
!
! Reference -
!
! Calling program - define_nstgrd_spec
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
nestcoords = 2

!---number of horizontal nest points
nohnstglbc = 0
nohnstglbu = 36
nohnstglbv = 127

!---number of vertical nest points
novnst = 0

!---type of 2-D output
inst2dtype = 1

CALL log_timer_out()


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
! Description -
!
! Reference -
!
! Calling program - define_nstgrd_locs
!

!************************************************************************
!

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


RETURN

END SUBROUTINE usrdef_nstgrd_abs

!========================================================================

SUBROUTINE usrdef_nstgrd_rel(iset,nhdat,nzdat,nonodes,hnests,zcoord,cnode)
!************************************************************************
!
! *usrdef_nstgrd_real* Define relative geographical positions of nested open
!                      boundaries
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Nested_Grids.f90  V2.9
!
! Description - optos csm
!
! Reference -
!
! Calling program - define_nstgrd_locs
!
! Module calls - close_file, open_file
!
!************************************************************************
!
USE datatypes
USE iopars
USE inout_routines, ONLY: close_file, open_file
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
INTEGER, SAVE :: iunit
INTEGER :: inode
REAL, DIMENSION(2,2,nhdat) :: weights


procname(pglev+1) = 'usrdef_nstgrd_rel'
CALL log_timer_in()

IF (cnode.EQ.'U') THEN
   CALL open_file(iunit,'csm_nest.dat','IN','A')
   READ (iunit,'(/)')
ENDIF

inode_110: DO inode=1,nonodes
   READ (iunit,*)
   READ (iunit,9001) hnests(:,inode)%icoord
   READ (iunit,*)
   READ (iunit,9001) hnests(:,inode)%jcoord
   READ (iunit,*)
   READ (iunit,9002) weights
   hnests(:,inode)%weights(1,1) = weights(1,1,:)
   hnests(:,inode)%weights(2,1) = weights(2,1,:)
   hnests(:,inode)%weights(1,2) = weights(1,2,:)
   hnests(:,inode)%weights(2,2) = weights(2,2,:)
ENDDO inode_110

IF (cnode.EQ.'V') CALL close_file(iunit,'A')

CALL log_timer_out()


RETURN

9001 FORMAT (500I4)
9002 FORMAT (500ES15.7)

END SUBROUTINE usrdef_nstgrd_rel
