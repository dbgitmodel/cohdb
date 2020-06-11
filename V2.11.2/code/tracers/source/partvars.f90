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

MODULE partvars
!************************************************************************
!
! *partvars* Arrays for the Lagrangian particle module
!
! Author - Valerie Duliere
!
! Version - @(COHERENS)partvars.f90  V2.11
!  
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description - 
!
!************************************************************************
!
USE datatypes
USE parttypes  
  
IMPLICIT NONE

!---particle attributes
TYPE (PartAtts), ALLOCATABLE, DIMENSION(:) :: part

!---particle clouds
TYPE (PartCloud), ALLOCATABLE, DIMENSION(:) :: cloud 

!---particulate matter
REAL, ALLOCATABLE, DIMENSION(:) :: diamconc, massconc, rhoconc, volumeconc

!---volume distributions
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ageconc_tot, conc_tot
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ageconc, conc  

!---masks at velocity nodes
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: maskatc, maskatu, maskatv

!---vertical diffusion coefficient
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: vdifcoefpart

!---random numbers for particle diffusion
REAL, ALLOCATABLE, DIMENSION(:) :: randirpart, ranmagpart

!---particle output
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:,:) :: outpvars
TYPE (Fileparams), ALLOCATABLE, DIMENSION(:) :: part1d
TYPE (OutPartParams), ALLOCATABLE, DIMENSION(:) :: outppars

SAVE

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*ageconc*      REAL    Age concentrations for each contaminant. Unit is
!                       defined by the the parameter ptime_unit.
!*ageconc_tot*  REAL    Age concentrations for all contaminants. Unit is
!                       defined by the the parameter ptime_unit.
!*cloud*        DERIVED Parameters for cloud (ship) releases  
!*conc*         REAL    Volume distributions for each contaminant        [1/m^3]
!*conc_tot*     REAL    Volume distributions for all contaminants        [1/m^3]
!*diamconc*     REAL    Particle size                                        [m]
!*maskatc*      LOGICAL .FALSE. at land cells, .TRUE. otherwise
!*maskatu*      LOGICAL .FALSE. at closed U-nodes, .TRUE. otherwise
!*maskatv*      LOGICAL .FALSE. at closed V-nodes, .TRUE. otherwise
!*massconc*     REAL    Particle mass                                       [kg]
!*ouptpars*     DERIVED Parameters for particle trajectory output
!*outpvars*     DERIVED Variable attributes for particle trajectory output
!*part*         DERIVED Particle attributes
!*part1d*       DERIVED Attributes of particle trajectories output files
!*randirpart*   REAL    Random numbers for random walk (direction of horizontal
!                       diffusive current)
!*ranmagpart*   REAL    Random numbers for random walk (magnitude of diffusive
!                       current)
!*rhoconc*      REAL    Particle mass density                            kg/m^3]

!*volumeconc*   REAL    Particle volume                                    [m^3]
!
!------------------------------------------------------------------------------

END MODULE partvars
