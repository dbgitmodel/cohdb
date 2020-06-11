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

MODULE partpars
!************************************************************************
!
! *partpars* Parameters for the Lagrangian particle module
!
! Author - Valerie Duliere
!
! Version - @(COHERENS)partpars.f90  V2.11
!
! $Date: 2017-04-12 16:52:46 +0200 (Wed, 12 Apr 2017) $
!
! $Revision: 1012 $
!
! Description - 
!
!************************************************************************
!
USE syspars
  
IMPLICIT NONE

LOGICAL :: back_drift, for_drift
CHARACTER (LEN=lentime) :: RefDateTime_part
INTEGER, PARAMETER :: nsize_sb = 2
INTEGER :: icpart, noclouds, nolabels, nopart, nosetspart, prangenlw, &
         & prangenwalk, ptime_unit
REAL :: cdrift_fac = 1.0, deltpart, wdrift_angle1 = 40.0, wdrift_angle2 = 8.0, &
      & wdrift_slope = 0.0315, xdifpart_cst = 0.0, ydifpart_cst = 0.0, &
      & zdifpart_cst = 0.0

SAVE

!
! Name                Type    Purpose
!------------------------------------------------------------------------------
!*RefDateTime_part*   CHAR    Reference date/time used to determine particle
!                             release times
!*back_drift*         LOGICAL .TRUE./.FALSE. for backward/forward drifting
!*cdrift_fac*         REAL    Current drift coefficient
!*delpart*            REAL    Time step for tracking (=delt3d for forward and
!                             -delt3d for backward tracking    
!*for_drift*          LOGICAL .TRUE./.FALSE. for forward/backward drifting
!*icpart*             INTEGER Counter for update of physical data by the
!                             particle module
!*noclouds*           INTEGER Number of clouds
!*nolabels*           INTEGER Number of concentrations (labels)
!*nopart*             INTEGER Total number of particles
!*nosetspart*         INTEGER Number of trajectory output files ("sets")
!*nsize_sb*           INTEGER Size (number of columns/rows) of the box
!                             around the current location where a search is
!                             performed to locate the partcicle's updated
!                             location
!*npbnd*              INTEGER Particle inside a domain less than npbnd grid
!                             cells of the physical domain edges are taken as
!                             out of domain
!*prangenlw*          INTEGER Random generator for crosswind drift factor 
!*prangenwalk*        INTEGER Random generator for calculating random walk
!                             diffusive velocities
!*ptime_unit*         INTEGER Time unit used for calculating release and aging
!                             times
!                      = 0 => seconds and milliseconds
!                      = 1 => seconds
!                      = 2 => minutes
!                      = 3 => hours
!                      = 4 => days
!*wdrift_angle1*      REAL    Coefficient for calculating the wind drift angle
!                                                                          [rad]
!*wdrift_angle2*      REAL    Coefficient for calculating the wind drift angle
!                                                                          [rad]
!*xdifpart_cst*       REAL    Constant horizontal diffusivity coefficient in the
!                             X-direction                                [m^2/s]
!*ydifpart_cst*       REAL    Constant horizontal diffusivity coefficient in the
!                             Y-direction                                [m^2/s]
!*zdifpart_cst*       REAL    Constant horizontal diffusivity coefficient in the
!                             vertical direction                         [m^2/s]
!
!------------------------------------------------------------------------------

END MODULE partpars
