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

MODULE currents
!************************************************************************
!
! *currents* Velocity arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)currents.f90  V2.10.2
!
! $Date: 2016-10-10 16:36:15 +0200 (Mon, 10 Oct 2016) $
!
! $Revision: 973 $
!
! Description - 
!
!************************************************************************
!

IMPLICIT NONE

REAL, ALLOCATABLE, DIMENSION(:,:) :: p2dbcgradatu, udevint, udfvel, udvel, &
                                   & udvel_old, umpred, umvel, umvel_old 
REAL, ALLOCATABLE, DIMENSION(:,:) :: p2dbcgradatv, vdevint, vdfvel, vdvel, &
                                   & vdvel_old, vmpred, vmvel, vmvel_old
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: p3dbcgradatu, ufvel, uvel, uvel_old
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: p3dbcgradatv, vfvel, vvel, vvel_old
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: wfvel, wvel, wphys

SAVE

!
! Name          Type  Purpose
!------------------------------------------------------------------------------
!*p2dbcgradatu* REAL  Depth-integrated X-component of baroclinic pressure
!                     gradient                                        [m^2/s^2]
!                     or X-component of surface pressure gradient (1-D case)
!                                                                       [m/s^2]
!*udevint*      REAL  Depth-integrated terms in 2-D U-equation        [m^2/s^2]
!*udfvel*       REAL  X-component of filtered depth-integrated current  [m^2/s]
!*udvel*        REAL  X-component of depth-integrated current           [m^2/s]
!*udvel_old*    REAL  X-component of depth-integrated current at old time
!                     level                                             [m^2/s]
!*umpred*       REAL  X-component of depth-mean predicted current       [m^2/s]
!*umvel*        REAL  X-component of depth-mean current                   [m/s]
!*umvel_old*    REAL  X-component of depth-mean current at old time       [m/s]

!*p2dbcgradatv* REAL  Depth-integrated Y-component of baroclinic pressure
!                     gradient                                        [m^2/s^2]
!                     or Y-component of surface pressure gradient (1-D case)
!                                                                       [m/s^2]
!*vdevint*      REAL  Depth-integrated terms in 2-D V-equation        [m^2/s^2]
!*vdfvel*       REAL  Y-component of filtered depth-integrated current  [m^2/s]
!*vdvel*        REAL  Y-component of depth-integrated current           [m^2/s]
!*vdvel_old*    REAL  Y-component of depth-integrated current at old time
!                     level                                             [m^2/s]
!*vmpred*       REAL  Y-component of depth-mean predicted current       [m^2/s]
!*vmvel*        REAL  Y-component of depth-mean current                   [m/s]
!*vmvel_old*    REAL  Y-component of depth-mean current at old time       [m/s]
!*p3dbcgradatu* REAL  X-component of baroclinic pressure gradient       [m/s^2]
!*ufvel*        REAL  X-component of advective "filtered" velocity        [m/s]
!*uvel*         REAL  X-component of velocity                             [m/s]
!*uvel_old*     REAL  X-component of velocity at old time level           [m/s]
!*p3dbcgradatv* REAL  Y-component of baroclinic pressure gradient       [m/s^2]
!*vfvel*        REAL  Y-component of advective "filtered" velocity        [m/s]
!*vvel          REAL  X-component of velocity                             [m/s]
!*vvel_old*     REAL  Y-component of velocity at old time level           [m/s]
!*wfvel*        REAL  Vertical (transformed) advective velocity           [m/s]
!*wvel*         REAL  Vertical (transformed) velocity                     [m/s]
!*wphys*        REAL  Vertical (physical) velocity                        [m/s]
!
!************************************************************************
!

END MODULE currents
