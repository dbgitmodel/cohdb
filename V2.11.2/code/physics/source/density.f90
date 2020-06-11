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

MODULE density
!************************************************************************
!
! *density* Density arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)density.f90  V2.0
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

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: beta_sal, beta_temp, dens, sal, temp

SAVE

!
! Name      Type  Purpose
!------------------------------------------------------------------------------
!*beta_sal* REAL  Salinity expansion coefficient                       [1/PSU] 
!*beta_temp*REAL  Temperature expansion coefficient                  [1/deg C]
!*dens*     REAL  Mass density                                        [kg/m^3]
!*sal*      REAL  Salinity                                               [PSU]
!*temp*     REAL  Temperature                                          [deg C]
!
!************************************************************************

END MODULE density
