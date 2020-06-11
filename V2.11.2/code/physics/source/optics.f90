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

MODULE optics
!************************************************************************
!
! *optics* Optical arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)optics.f90  V2.0
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

REAL, ALLOCATABLE, DIMENSION(:,:) :: optattcoef2, qrad 
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: radiance

SAVE

!
! Name         Type  Purpose
!-----------------------------------------------------------------------------
!*optattcoef2* REAL  Invserse optical attenuation depth                 [m^-1]
!*qrad*        REAL  Surface solar irradiance                          [W/m^2]
!*radiance*    REAL  Solar irradiance                                  [W/m^2]
!
!*****************************************************************************
!

END MODULE optics
