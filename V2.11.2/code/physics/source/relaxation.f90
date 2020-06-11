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

MODULE relaxation
!************************************************************************
!
! *relaxation* Arrays for relaxation at open boundaries
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)relaxation.f90  V2.0
!
! $Date: 2015-04-10 11:30:34 +0200 (Fri, 10 Apr 2015) $
!
! $Revision: 841 $
!
! Description - 
!
!************************************************************************
!

IMPLICIT NONE

INTEGER :: norlxzones = 0
INTEGER, DIMENSION(2) :: inodesrlx = 0
INTEGER, ALLOCATABLE, DIMENSION(:) :: idirrlx, iposrlx, ityprlx, jposrlx, &
                                    & ncrlx, nrrlx
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: indexrlxatc, indexrlxatuv
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: rlxwghtatc, rlxwghtatuv

SAVE

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*norlxzones*   INTEGER Number of relaxation zones
!*inodesrlx*    INTEGER Enables/disables (1/0) 'C'-type, 'U'-node and 'V'-node
!                       quantities
!*idirrlx*      INTEGER Orientation of relaxation zone
!                       (1,2,3,4 = West,East,South,North)
!*iposrlx*      INTEGER X-index of soutwest corner for each relaxation zone
!*ityprlx*      INTEGER Type of relaxation scheme near open boundaries
!                  = 1 => linear
!                  = 2 => quadratic
!                  = 3 => hyperbolic
!*jposrlx*      INTEGER Y-index of soutwest corner for each relaxation zone
!*ncrlx*        INTEGER X-dimension of relaxation zones
!*nrrlx*        INTEGER Y-dimension of relaxation zones
!*indexrlxatc*  INTEGER Indices of U- and V-o.b. where external profile is
!                       defined ('C'-type quantities)
!*indexrlxatuv* INTEGER Indices of U- and V-o.b. where external profile is
!                       defined ('UV'-type quantities)
!*rlxwghtatc*   REAL    Relaxation factor in X- and Y-direction for 'C'-type
!                       quantities
!*rlxwghtatuv*  REAL    Relaxation factor in X- and Y-direction for 'UV'-type
!                       quantities
!
!************************************************************************
!

END MODULE relaxation
