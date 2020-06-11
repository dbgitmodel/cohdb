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

MODULE parttypes
!************************************************************************
!
! *parttypes* Derived type definitions for the Lagrangain particle module
!
! Author - Valerie Duliere
!  
! Version - @(COHERENS)parttypes.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description - 
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

!---attributes of the slick initial position
TYPE :: PartCloud
   INTEGER :: kdistype, label, nopart, noreleases, state
   REAL :: length, orientation, thick, width
END TYPE PartCloud

!---attributes of particles for the Lagrangian model
TYPE :: PartAtts
   LOGICAL :: center
   INTEGER :: drift_state, icoord, jcoord, kcoord, kdistype, label, ntstart, &
            & state
   REAL :: displx, disply, displz, xcoord, ycoord
   REAL (KIND=kndrlong) :: age, tstart, xpos, ypos, zpos
END TYPE PartAtts

!---attributes for trajectory output
TYPE :: OutPartParams
   LOGICAL :: aging
   CHARACTER (LEN=lentime) :: enddate, refdate, startdate
   INTEGER :: ktype, label, nodim, nopartout, nstepout, ntype, time_format
   INTEGER, DIMENSION(3) :: tlims
   REAL :: deltout
END TYPE OutPartParams

SAVE

! Name             Type    Purpose
!------------------------------------------------------------------------------
!*PartCloud*       DERIVED Attributes for initial particle clouds
!*%kdistype*       INTEGER Selects type vertical location of the cloud centre
!                    = 0 =>  at a specified z-level
!                    = 1 =>  at a specified vertical grid level at the C-node 
!                    = 2 =>  at a fixed distance from the sea bed
!                    = 3 =>  at a fixed distance from the sea surface
!*%label*          INTEGER Type of cloud concentration
!*%length          REAL    Cloud length in the Y-direction                   [m]
!*%nopart*         INTEGER Number of particles within the cloud at each release
!*%noreleases*     INTEGER Number of particle releases per cloud
!*%orientation     REAL    Cloud ellipse orientation
!                          (0-North; 90-East; 180-South)               [degrees]
!*%thick*          REAL    Cloud length in the vertical direction            [m]
!*%state*          INTEGER Type of cloud particles
!                    = 1 => surface
!                    = 2 => submerged within water column 
!*%width           REAL    Cloud length in the X-direction                   [m]
!
!*PartAtts*        DERIVED Particle attributes
!*%age*            LONG    Particle age
!*%center*         LOGICAL .TRUE. if particle is at a cloud center
!*displx*          REAL    Particle displacement in the X-direction during one
!                          time step                    [m or degrees longitude]
!*disply*          REAL    Particle displacement in the Y-direction during one
!                          time step                     [m or degrees latitude]
!*displz*          REAL    Particle displacement in the vertical direction
!                          during one time step                              [m]
!*%drift_state*    INTEGER Particle drift state
!                    = 1 => particle has not yet been released
!                    = 2 => particle is released and drifting
!                    = 3 => particle is located on a dry/land cell
!                    = 4 => particle is outside the domain
!                    = 5 => particle is deposited on the sea bed
!*%icoord*         INTEGER X-index coordinate with respect to the lower left
!                          corner of the grid cell where the particle is
!                          residing
!*%jcoord*         INTEGER Y-index coordinate with respect to the lower left
!                          corner of the grid celll where the particle is
!                          residing
!*%kcoord*         INTEGER Vertical grid index with respect to the lower left
!                          corner of the grid celll where the particle is
!                          residing
!*%kdistype*       INTEGER Selects type of initial location in the vertical
!                    = 0 =>  at a specified z-level
!                    = 1 =>  at a specified vertical grid level at the C-node 
!                    = 2 =>  at a fixed distance from the sea bed
!                    = 3 =>  at a fixed distance from the sea surface
!*%label*          INTEGER Particle label to allow for different concentrations
!*%ntstart*        INTEGER Time index (with respect to the simulation start
!                          time) when the particle is released
!*%state*          INTEGER Type of particle drift
!                    = 1 => surface drifting
!                    = 2 => immersed in the water column
!*%tstart*         LONG    Particle release time  with respect to the start
!                          time of the simulation (negative in case of backward
!                          tracking). Time unit is given by the parameter
!                          ptime_unit.
!*%xcoord*         REAL    Normalised X-coordinate with respect to the lower
!                          left corner of the grid cell where the particle is
!                          residing (between 0 and 1)
!*%xpos*           REAL    X-coordinate at the particle's position
!                                                       [m or degrees longitude]
!*%ycoord*         REAL    Normalised Y-coordinate with respect to the lower
!                          left corner of the grid cell where the particle is
!                          residing (between 0 and 1)
!*%ypos*           REAL    Y-coordinate at the particle's position
!                                                        [m or degrees latitude]
!*%zpos*           REAL    Z-coordinate at particle's location
!
!*OutpartParams*   DERIVED Parameters for particle traject output
!*%aging*          LOGICAL Include particle age if .TRUE.
!*%deltout*        INTEGER Output time step
!*%enddate*        CHAR    End output date/time
!*%ktype*          INTEGER Selects type of output for vertical particle location
!                    = 0 =>  z-coordinate
!                    = 1 =>  distance from sea bed
!                    = 2 =>  distance from sea surface
!*%label*          INTEGER Particle label selected for output (0 means all
!                          labels are selected)
!*%nodim*          INTEGER Selects type output coordinates
!                    = 0 => horizontal and vertical coordintates for all
!                           particles selected by ntype
!                    = 2 => horizontal coodinates for surface drifiting
!                           particles selected by ntype
!                    = 3 => horizontal and vertical coordinates for immersed
!                           particles selected by ntype
!*%nopartout*      INTEGER Number of output particles
!*%ntype*          INTEGER Type of particles selected for output
!                    = 0 => all particles
!                    = 1 => particles with a given label
!                    = 2 => particles at cloud denters
!                    = 3 => particles at cloud centers with a given label
!*%nstepout*       INTEGER Number of outpout times
!*%refdate*        CHAR    Refence date/time
!*%tlims*          INTEGER Start/End/Stride time index
!*%startdate*      CHAR    Start output date/time
!
!************************************************************************
!

END MODULE parttypes
