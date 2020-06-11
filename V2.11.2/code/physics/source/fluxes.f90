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

MODULE fluxes
!************************************************************************
!
! *fluxes* Bottom/surface fluxes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)fluxes.f90  V2.11.1
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

!---bottom stress
REAL, ALLOCATABLE, DIMENSION(:,:) :: bdragcoefatc, bdragcoefatu, bdragcoefatv, &
                                   & bfricatu, bfricatv, bstresatc, &
                                   & bstresatc_max, bstresatc_wav, bstresatu, &
                                   & bstresatv, fwave, ubstresatu, vbstresatv, &
                                   & wavethickatc, zaroughatc, zroughatc, &
                                   & zroughatu, zroughatv

!---bottom currents
REAL, ALLOCATABLE, DIMENSION(:,:) :: bcuratc, bcuratu, bcuratv
!---surface stress
REAL, ALLOCATABLE, DIMENSION(:,:) :: sstresatc, usstresatc, usstresatu, &
                                   & vsstresatc, vsstresatv
!---drag and exchange coefficients
REAL, ALLOCATABLE, DIMENSION(:,:) :: cds, ces, chs
!---surface fluxes
REAL, ALLOCATABLE, DIMENSION(:,:) :: qlatent, qlwave, qnonsol, qsensible, &
                                   & ssalflux

!---zero work space array
REAL, ALLOCATABLE, DIMENSION(:,:) :: zeros2d

SAVE

!
! Name         Type Purpose
!------------------------------------------------------------------------------
!*bdragcoefatc* REAL Bottom drag coefficient at C-nodes
!*bdragcoefatu* REAL Bottom drag coefficient at U-nodes
!*bdragcoefatv* REAL Bottom drag coefficient at V-nodes
!*bcuratc*      REAL Bottom current at C-nodes                            [m/s]
!*bcuratu*      REAL Bottom current at U-nodes                            [m/s]
!*bcuratv*      REAL Bottom current at V-nodes                            [m/s]
!*bfricatu*     REAL Bottom friction velocity at U-nodes                  [m/s]
!*bfricatv*     REAL Bottom friction velocity at V-nodes                  [m/s]
!*bstresatc*    REAL (Mean) bottom stress at C-nodes                  [m^2/s^2]
!*bstresatc_max*REAL Maximum current-wave bottom stress at C-nodes    [m^2/s^2]
!*bstresatc_wav*REAL Bottom stress (wave only) at C-nodes             [m^2/s^2]
!*bstresatu*    REAL (Mean) bottom stress at U-nodes                  [m^2/s^2]
!*bstresatv*    REAL (Mean) bottom stress at V-nodes                  [m^2/s^2]
!*cds*          REAL Surface drag coefficient
!*ces*          REAL Exchange coefficient for latent heat flux (Dalton number)
!*chs*          REAL Exchange coefficient for sensible heat flux
!                    (Stanton number)
!*fwave*        REAL Wave friction factor                                   [-]
!*qlatent*      REAL Latent upward surface heat flux                    [W/m^2]
!*qlwave*       REAL Long-wave upward surface heat flux                 [W/m^2]
!*qnonsol*      REAL Non-solar upward surface heat flux                 [W/m^2]
!*qsensible*    REAL Sensible upward surface heat flux                  [W/m^2]
!*ssalflux*     REAL Surface salinity flux                            [PSU*m/s]
!*sstresatc*    REAL Surface stress at C-nodes                        [m^2/s^2]
!*ubstresatu*   REAL X-component of the (totla) bottom stress at the U-nodes
!                                                                     [m^2/s^2]
!*usstresatc*   REAL X-component of surface stress at C-nodes         [m^2/s^2]
!*usstresatu*   REAL X-component of surface stress at U-nodes         [m^2/s^2]
!*vbstresatv*   REAL Y-component of the (total) bottom stress at the V-nodes
!                                                                     [m^2/s^2]
!*vsstresatc*   REAL Y-component of surface stress at C-nodes         [m^2/s^2]
!*vsstresatv*   REAL Y-component of surface stress at V-nodes         [m^2/s^2]
!*wavethickatc* REAL Thickness of the wave boundary layer at C-nodes   [m]
!*zaroughatc*   REAL Apparent bottom roughness at C-nodes                   [m]
!*zroughatc*    REAL Bottom roughness at C-nodes                            [m]
!*zroughatu*    REAL Bottom roughness at U-nodes                            [m]
!*zroughatv*    REAL Bottom roughness at V-nodes                            [m]
!*zeros2d*      REAL Zero array used as argument in transport modules
!
!************************************************************************
!

END MODULE fluxes
