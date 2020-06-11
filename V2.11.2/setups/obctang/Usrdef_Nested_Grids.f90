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
! Version - @(COHERENS)Usrdef_Nested_Grids.f90  V2.7
!
! $Date: 2015-10-06 17:42:25 +0200 (Tue, 06 Oct 2015) $
!
! $Revision: 886 $
!
! Description - test case obctang
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
! Version - @(COHERENS)Usrdef_Nested_Grids.f90  V2.7
!
! Description - test case obctang
!
! Reference -
!
! Calling program - define_nstgrd_spec
!
!************************************************************************
!
USE gridpars
USE iopars
USE nestgrids
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_nstgrd_spec'
CALL log_timer_in()

!---type of grids
nestcoords = 2

!---number of horizontal nest points
nohnstglbc = 0
nohnstglbu = 200  
nohnstglbv = 400 
nohnstglbx = 198 
nohnstglby = 398

!---number of vertical nest points
novnst = MERGE(0,nz,iopt_grid_nodim.EQ.2)

!---type of 2-D output
SELECT CASE (runtitle(8:8))
   CASE ('A','B','E','F'); inst2dtype = 2
   CASE ('C','D','G','H'); inst2dtype = 1
END SELECT

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
! Version - @(COHERENS)Usrdef_Nested_Grids.f90  V2.7
!
! Description - test case obctang
!
! Reference -
!
! Calling program - define_nstgrd_locs
!
!************************************************************************
!
USE iopars
USE physpars
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
INTEGER :: ii, jj, k
REAL, DIMENSION(nhdat) :: xcoord, ycoord


procname(pglev+1) = 'usrdef_nstgrd_rel'
CALL log_timer_in()

!
!1. Horizontal locations
!-----------------------
!

SELECT CASE (TRIM(cnode))

!
!1.1 U-nodes
!-----------
!

   CASE ('U')

!     ---UtoU-West
      hnests(1:nhdat/2,1)%icoord = 51
      hnests(1,1)%jcoord = 25; hnests(nhdat/2,1)%jcoord = 75
      hnests(2:nhdat/2-2:2,1)%jcoord = (/(jj,jj=26,74)/)
      hnests(3:nhdat/2-1:2,1)%jcoord = (/(jj,jj=26,74)/)
      xcoord(1:nhdat/2) = 0.0; ycoord(1) = 0.75; ycoord(nhdat/2) = 0.25
      ycoord(2:nhdat/2-2:2) = 0.25; ycoord(3:nhdat/2-1:2) = 0.75

!     ---UtoU-East
      hnests(nhdat/2+1:nhdat,1)%icoord = 151
      hnests(nhdat/2+1,1)%jcoord = 25; hnests(nhdat,1)%jcoord = 75
      hnests(nhdat/2+2:nhdat-2:2,1)%jcoord = (/(jj,jj=26,74)/)
      hnests(nhdat/2+3:nhdat-1:2,1)%jcoord = (/(jj,jj=26,74)/)
      xcoord(nhdat/2+1:nhdat) = 0.0
      ycoord(nhdat/2+1) = 0.7; ycoord(nhdat) = 0.25
      ycoord(nhdat/2+2:nhdat-2:2) = 0.25
      ycoord(nhdat/2+3:nhdat-1:2) = 0.75

!     ---UtoU weight factors
      hnests(:,1)%weights(1,1) = (1.0-xcoord)*(1.0-ycoord)
      hnests(:,1)%weights(2,1) = xcoord*(1.0-ycoord)
      hnests(:,1)%weights(1,2) = ycoord*(1.0-xcoord)
      hnests(:,1)%weights(2,2) = xcoord*ycoord

!     ---CtoU-West
      hnests(1:nhdat/2,2)%icoord = 50
      hnests(1,2)%jcoord = 25; hnests(nhdat/2,2)%jcoord = 75
      hnests(2:nhdat/2-2:2,2)%jcoord = (/(jj,jj=26,74)/)
      hnests(3:nhdat/2-1:2,2)%jcoord = (/(jj,jj=26,74)/)
      xcoord(1:nhdat/2) = 0.5; ycoord(1) = 0.75; ycoord(nhdat/2) = 0.25
      ycoord(2:nhdat/2-2:2) = 0.25; ycoord(3:nhdat/2-1:2) = 0.75

!     ---CtoU-East
      hnests(nhdat/2+1:nhdat,2)%icoord = 150
      hnests(nhdat/2+1,2)%jcoord = 25; hnests(nhdat,2)%jcoord = 75
      hnests(nhdat/2+2:nhdat-2:2,2)%jcoord = (/(jj,jj=26,74)/)
      hnests(nhdat/2+3:nhdat-1:2,2)%jcoord = (/(jj,jj=26,74)/)
      xcoord(nhdat/2+1:nhdat) = 0.5
      ycoord(nhdat/2+1) = 0.75; ycoord(nhdat) = 0.25
      ycoord(nhdat/2+2:nhdat-2:2) = 0.25; ycoord(nhdat/2+3:nhdat-1:2) = 0.75

!     ---CtoU weight factors
      hnests(:,2)%weights(1,1) = (1.0-xcoord)*(1.0-ycoord)
      hnests(:,2)%weights(2,1) = xcoord*(1.0-ycoord)
      hnests(:,2)%weights(1,2) = ycoord*(1.0-xcoord)
      hnests(:,2)%weights(2,2) = xcoord*ycoord

!
!1.2 V-nodes
!-----------
!

   CASE ('V')

!     ---VtoV-South
      hnests(1,1)%icoord = 50; hnests(nhdat/2,1)%icoord = 150
      hnests(2:nhdat/2-2:2,1)%icoord = (/(ii,ii=51,149)/)
      hnests(3:nhdat/2-1:2,1)%icoord = (/(ii,ii=51,149)/)
      hnests(1:nhdat/2,1)%jcoord = 26
      xcoord(1) = 0.75; xcoord(nhdat/2) = 0.25
      xcoord(2:nhdat/2-2:2) = 0.25; xcoord(3:nhdat/2-1:2) = 0.75
      ycoord(1:nhdat/2) = 0.0

!     ---VtoV-North
      hnests(nhdat/2+1,1)%icoord = 50; hnests(nhdat,1)%icoord = 150
      hnests(nhdat/2+2:nhdat-2:2,1)%icoord = (/(ii,ii=51,149)/)
      hnests(nhdat/2+3:nhdat-1:2,1)%icoord = (/(ii,ii=51,149)/)
      hnests(nhdat/2+1:nhdat,1)%jcoord = 76
      xcoord(nhdat/2+1) = 0.75; xcoord(nhdat) = 0.25
      xcoord(nhdat/2+2:nhdat-2:2) = 0.25; xcoord(nhdat/2+3:nhdat-1:2) = 0.75
      ycoord(nhdat/2+1:nhdat) = 0.0

!     ---VtoV weight factors
      hnests(:,1)%weights(1,1) = (1.0-xcoord)*(1.0-ycoord)
      hnests(:,1)%weights(2,1) = xcoord*(1.0-ycoord)
      hnests(:,1)%weights(1,2) = ycoord*(1.0-xcoord)
      hnests(:,1)%weights(2,2) = xcoord*ycoord

!     ---CtoV-South
      hnests(1,2)%icoord = 50; hnests(nhdat/2,2)%icoord = 150
      hnests(2:nhdat/2-2:2,2)%icoord = (/(ii,ii=51,149)/)
      hnests(3:nhdat/2-1:2,2)%icoord = (/(ii,ii=51,149)/)
      hnests(1:nhdat/2,2)%jcoord = 25
      xcoord(1) = 0.75; xcoord(nhdat/2) = 0.25
      xcoord(2:nhdat/2-2:2) = 0.25; xcoord(3:nhdat/2-1:2) = 0.75
      ycoord(1:nhdat/2) = 0.5

!     ---CtoV-North
      hnests(nhdat/2+1,2)%icoord = 50; hnests(nhdat,2)%icoord = 150
      hnests(nhdat/2+2:nhdat-2:2,2)%icoord = (/(ii,ii=51,149)/)
      hnests(nhdat/2+3:nhdat-1:2,2)%icoord = (/(ii,ii=51,149)/)
      hnests(nhdat/2+1:nhdat,2)%jcoord = 75
      xcoord(nhdat/2+1) = 0.75; xcoord(nhdat) = 0.25
      xcoord(nhdat/2+2:nhdat-2:2) = 0.25; xcoord(nhdat/2+3:nhdat-1:2) = 0.75
      ycoord(nhdat/2+1:nhdat) = 0.5

!     ---CtoV weight factors
      hnests(:,2)%weights(1,1) = (1.0-xcoord)*(1.0-ycoord)
      hnests(:,2)%weights(2,1) = xcoord*(1.0-ycoord)
      hnests(:,2)%weights(1,2) = ycoord*(1.0-xcoord)
      hnests(:,2)%weights(2,2) = xcoord*ycoord

!
!1.3 X-nodes
!-----------
!

   CASE ('X')

!     ---VtoX-West
      hnests(1:nhdat/2,1)%icoord = 50
      hnests(1,1)%jcoord = 26
      hnests(2:nhdat/2-1:2,1)%jcoord = (/(jj,jj=27,75)/)
      hnests(3:nhdat/2:2,1)%jcoord = (/(jj,jj=27,75)/)
      xcoord(1:nhdat/2) = 0.25; ycoord(1) = 0.5
      ycoord(2:nhdat/2-1:2) = 0.0; ycoord(3:nhdat/2:2) = 0.5

!     ---VtoX-East
      hnests(nhdat/2+1:nhdat,1)%icoord = 150
      hnests(nhdat/2+1,1)%jcoord = 26
      hnests(nhdat/2+2:nhdat-1:2,1)%jcoord = (/(jj,jj=27,75)/)
      hnests(nhdat/2+3:nhdat:2,1)%jcoord = (/(jj,jj=27,75)/)
      xcoord(nhdat/2+1:nhdat) = 0.25; ycoord(nhdat/2+1) = 0.5
      ycoord(nhdat/2+2:nhdat-1:2) = 0.0; ycoord(nhdat/2+3:nhdat:2) = 0.5

!     ---VtoX weight factors
      hnests(:,1)%weights(1,1) = (1.0-xcoord)*(1.0-ycoord)
      hnests(:,1)%weights(2,1) = xcoord*(1.0-ycoord)
      hnests(:,1)%weights(1,2) = ycoord*(1.0-xcoord)
      hnests(:,1)%weights(2,2) = xcoord*ycoord

!
!1.4 Y-nodes
!-----------
!

   CASE ('Y')

!     ---UtoY-South
      hnests(1,1)%icoord = 51
      hnests(2:nhdat/2-1:2,1)%icoord = (/(ii,ii=52,150)/)
      hnests(3:nhdat/2:2,1)%icoord = (/(ii,ii=52,150)/)
      hnests(1:nhdat/2,1)%jcoord = 25
      xcoord(1) = 0.5; xcoord(2:nhdat/2-1:2) = 0.0
      xcoord(3:nhdat/2:2) = 0.5; ycoord(1:nhdat/2) = 0.25

!     ---UtoY-North
      hnests(nhdat/2+1,1)%icoord = 51
      hnests(nhdat/2+2:nhdat-1:2,1)%icoord = (/(ii,ii=52,150)/)
      hnests(nhdat/2+3:nhdat:2,1)%icoord = (/(ii,ii=52,150)/)
      hnests(nhdat/2+1:nhdat,1)%jcoord = 75
      xcoord(nhdat/2+1) = 0.5; xcoord(nhdat/2+2:nhdat-1:2) = 0.0
      xcoord(nhdat/2+3:nhdat:2) = 0.5; ycoord(nhdat/2+1:nhdat) = 0.25

!     ---UtoY weight factors
      hnests(:,1)%weights(1,1) = (1.0-xcoord)*(1.0-ycoord)
      hnests(:,1)%weights(2,1) = xcoord*(1.0-ycoord)
      hnests(:,1)%weights(1,2) = ycoord*(1.0-xcoord)
      hnests(:,1)%weights(2,2) = xcoord*ycoord

END SELECT

!
!2. Vertical locations
!---------------------
!

IF (nzdat.GT.0) THEN
   k_210: DO k=1,nzdat
      zcoord(:,k) = depmean_cst*(1.0-(k-0.5)/REAL(nzdat))
   ENDDO k_210
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_nstgrd_rel
