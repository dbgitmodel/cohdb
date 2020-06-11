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

MODULE depths
!************************************************************************
!
! *depths*  Water depths and elevations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)depths.f90  V2.4
!
! $Date: 2018-02-21 09:29:52 +0100 (Wed, 21 Feb 2018) $
!
! $Revision: 1097 $
!
! Description - 
!
!************************************************************************
!

IMPLICIT NONE

REAL, ALLOCATABLE, DIMENSION(:,:) :: depmeanatc, deptotatc, deptotatc_old, &
                                   & deptotatc_prev
REAL, ALLOCATABLE, DIMENSION(:,:) :: depmeanatu, deptotatu, deptotatu_old
REAL, ALLOCATABLE, DIMENSION(:,:) :: depmeanatv, deptotatv, deptotatv_old
REAL, ALLOCATABLE, DIMENSION(:,:) :: depmeanatuv, deptotatuv, depmeanglb
REAL, ALLOCATABLE, DIMENSION(:,:) :: deptotatc_err, dzeta, zeta, zeta_old

SAVE

!
! Name           Type Purpose
!-----------------------------------------------------------------------------
!*depmeanatc*    REAL Mean water depth (local) at C-nodes                   [m]
!*deptotatc*     REAL Total water depth at C-nodes                          [m]
!*deptotatc_old* REAL Total water depth at C-nodes and old baroclinic time  [m]
!*deptotatc_prev*REAL Total water depth at C-nodes and previous 2-D or 3-D time
!                                                                           [m]
!*depmeanatu*    REAL Mean water depths at U-nodes                          [m]
!*deptotatu*     REAL Total water depth at U-nodes                          [m]
!*deptotatu_old* REAL Total water depth at U-nodes and old
!                     (baroclinic/barotropic) time step                     [m]
!*depmeanatv*    REAL Mean water depth at V-nodes                           [m]
!*deptotatv*     REAL Total water depth at V-nodes                          [m]
!*deptotatv_old* REAL Total water depth at V-nodes and
!                     (baroclinic/barotropic) time step                     [m]
!*depmeanatuv*   REAL Mean water depth at UV-nodes                          [m]
!*deptotatuv*    REAL Total water depth at UV-nodes                         [m]
!*depmeanglb*    REAL Mean water depth (global) at C- or -UV-nodes          [m]
!*deptotatc_err* REAL Error in total water depth C-nodes due to violation
!                     of mass conservation                                  [m]
!*dzeta*         REAL Free surface correction defined as the difference
!                     between the surface elevation at the next and previous
!                     iteration (implicit scheme only)                      [m]
!*zeta*          REAL Surface elevation                                     [m]
!*zeta_old*      REAL surface elevation at the start of the first (outer)
!                     iteration with the implicit scheme previous time step [m]
!
!************************************************************************
!

END MODULE depths
