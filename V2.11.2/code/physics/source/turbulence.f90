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

MODULE turbulence
!************************************************************************
!
! *turbulence* Turbulence arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)turbulence.f90  V2.0
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

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dissip, tke, tke_old, zlmix
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: buofreq2, shearfreq2

SAVE

!
! Name         Type  Purpose
!------------------------------------------------------------------------------
!*dissip*      REAL  Dissipation of turbulent kinetic energy             [W/kg]
!*tke*         REAL  Turbulent kinetic energy                            [J/kg]
!*tke_old*     REAL  Turbulent kinetic energy at old time step           [J/kg]
!*zlmix*       REAL  Turbulent mixing length                                [m]
!*buofreq2*    REAL  Squared buoyancy frequency                         [1/s^2]
!*shearfreq2*  REAL  Squared shear frequency                            [1/s^2]
!
!************************************************************************
!

END MODULE turbulence
