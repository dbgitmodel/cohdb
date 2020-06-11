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


MODULE wavevars
!************************************************************************
!
! *wavevars* Surface wave arrays
!
! Author - Alexabder Breugem and Patrick Luyten
!
! Version - @(COHERENS)wavevars.f90  V2.11.1
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description - 
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

!---grid arrays
LOGICAL, DIMENSION(:,:), ALLOCATABLE :: maskglbwav
REAL, DIMENSION(:,:), ALLOCATABLE :: gxcoordglbwav, gycoordglbwav

!---model arrays
REAL, ALLOCATABLE, DIMENSION(:,:) :: stokessource2du, stokessource2dv, &
      & udstokesatu, umbwdissipatc, umstokesatc, umstokesatu, umswdissipatc, &
      & vdstokesatv, vmbwdissipatc, vmstokesatc, vmstokesatv, vmswdissipatc, &
      & wavedir, waveexcurs, wavefreq, waveheight, wavenum, waveperiod, &
      & wavepres, wavevel
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ubwdissipatc, ustokesatc, ustokesatu, &
      & uswdissipatc, vbwdissipatc, vstokesatc, vstokesatv, vswdissipatc, &
      & wstokesatw

SAVE

!
! Name            Type    Purpose
!------------------------------------------------------------------------------
!*gxcoordglbwav*  REAL    Global X-coordinates of the wave model grid 
!*gycoordglbwav*  REAL    Global Y-coordinates of the wave model grid
!*maskglbwav*     LOGICAL Local mask array set to .FALSE. for land points on
!                         the wave model grid
!*stokessource2du*REAL    Depth-intergrated Stokes term in the 2-D U-equation 
!*stokessource2dv*REAL    Depth-intergrated Stokes term in the 2-D V-equation 
!*ubwdissipatc*   REAL    X-component of the 3-D bed wave dissipation force 
!                                                                         [m/s2]
!*udstokesatu*    REAL    X-component of the depth-integrated Stokes drift at
!                         U-nodes                                         [m2/s]
!*umbwdissipatc*  REAL    X-component of the 2-D bed wave dissipation force 
!                                                                        [m2/s2]
!*umstokesatc*    REAL    X-component of the depth-mean Stokes drift at C-nodes
!                                                                          [m/s]
!*umstokesatu*    REAL    X-component of the depth-mean Stokes drift at U-nodes
!                                                                          [m/s]
!*umswdissipatc*  REAL    X-component of the 2-D surface wave dissipation force
!                                                                        [m2/s2]
!*ustokesatc*     REAL    X-component of the Stokes drift at C-nodes       [m/s]
!*ustokesatu*     REAL    X-component of the Stokes drift at U-nodes       [m/s]
!*uswdissipatc*   REAL    X-component of the 3-D surface wave dissipation force 
!                                                                         [m/s2]
!*vbwdissipatc*   REAL    Y-component of the 3-D surface wave dissipation force
!                                                                         [m/s2]
!*vdstokesatv*    REAL    Y-component of the depth-integrated Stokes drift at
!                         the V-nodes                                     [m2/s]
!*vmbwdissipatc*  REAL    Y-component of the 2-D bottom wave dissipation force
!                                                                        [m2/s2]
!*vmstokesatc*    REAL    Y-component of the depth-mean Stokes drift at C-nodes
!                                                                          [m/s]
!*vmstokesatv*    REAL    Y-component of the depth-mean Stokes drift at V-nodes
!                                                                          [m/s]
!*vmswdissipatc*  REAL    Y-component of the 2-D surface wave dissipation force
!                                                                        [m2/s2]
!*vstokesatc*     REAL    Y-component of the Stokes drift at C-nodes       [m/s]
!*vstokesatv*     REAL    Y-component of the Stokes drift at V-nodes       [m/s]
!*vswdissipatc*   REAL    Y-component of the 3-D surface wave dissipation force
!                                                                         [m/s2]
!*wavedir*        REAL    Mean wave direction at C-nodes                   [rad]
!*wavexcurs*      REAL    Near-bottom wave excursion amplitude at C-nodes    [m]
!*wavefreq*       REAL    Peak wave frequency at C-nodes                 [rad/s]
!*waveheight*     REAL    Significant wave height at C-nodes                 [m]
!*wavenum*        REAL    Wave number at C-nodes                           [1/m]
!*waveperiod*     REAL    Peak wave period at C-nodes                        [s]
!*wavepres*       REAL    Wave induced pressure                          [m2/s2]
!*wavevel*        REAL    Near-bottom wave orbital velocity at C-nodes     [m/s]
!*wstokesatw*     REAL    Transformed vertical Stokes drift at W-nodes     [m/s]
!
!************************************************************************
!

END MODULE wavevars
