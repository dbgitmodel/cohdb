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

MODULE meteo
!************************************************************************
!
! *meteo* Meteorological data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)meteo.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
!************************************************************************
!

IMPLICIT NONE

REAL, ALLOCATABLE, DIMENSION(:,:) :: atmpres, uwindatc, uwindatc_old, &
                                   & vwindatc, vwindatc_old, windatc
REAL, ALLOCATABLE, DIMENSION(:,:) :: airtemp, cloud_cover, evapminprec,&
                                   & evaporation, precipitation, qspecdif, &
                                   & relhum, sst, tempdif, vappres_air

SAVE

!
! Name             Type Purpose
!------------------------------------------------------------------------------
!*airtemp*         REAL Air temperature                                  [deg C]
!*atmpres*         REAL Atmospheric pressure                             [N/m^2]
!*cloud_cover*     REAL Cloud cover (between 0 and 1)                        [-]
!*evapminprec*     REAL Evaporation minus precipitation rate          [kg/m^2/s]
!*evaporation*     REAL Evaporation rate                              [kg/m^2/s]
!*precipitation*   REAL Precipitation rate                            [kg/m^2/s]
!*qspecdif*        REAL Air-sea difference of specific humidity              [-]
!*relhum*          REAL Relative humidity (between 0 and 1)                  [-]
!*sst*             REAL Sea surface temperature                          [deg C]
!*uwindatc*        REAL X-component of 10m wind at C-nodes                 [m/s]
!*tempdif*         REAL Air-sea temperature difference                   [deg C] 
!*uwindatc_old*    REAL X-component of 10m wind at C-nodes and old time step
!                                                                          [m/s]
!*vappres_air*     REAL Saturated vapour pressure at reference height     [mbar]
!*vwindatc*        REAL Y-component of 10m wind at C-nodes                 [m/s]
!*vwindatc_old*    REAL Y-component of 10m wind at C-nodes and old time step
!                                                                          [m/s]
!*windatc*         REAL Relative wind speed at reference level             [m/s]
!
!************************************************************************
!

END MODULE meteo
