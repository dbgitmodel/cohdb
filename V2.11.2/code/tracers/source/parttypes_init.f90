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

MODULE parttypes_init
!************************************************************************
!
! *parttypes_init* Initialise derived type scalars and arrays for the particle
!                  model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)parttypes_init.f90  V2.11
!
! $Date: 2017-04-12 16:52:46 +0200 (Wed, 12 Apr 2017) $
!
! $Revision: 1012 $
!
! Description - 
!
! Routines - cloud_init, outppart_init, parts_init
!
!************************************************************************
!

IMPLICIT NONE
  
CONTAINS

!========================================================================

SUBROUTINE cloud_init(cloudpart)
!************************************************************************
!
! *cloud_init* Initialise derived type variable of type 'PartClouds'
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)parttypes_init.f90  V2.11
!
! Description - argument of rank 1
!
!************************************************************************
!
USE parttypes

!
!*Arguments
!
TYPE (PartCloud), INTENT(OUT), DIMENSION(:) :: cloudpart

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cloudpart* DERIVED Attributes of particle clouds
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: icld, psize


psize = SIZE(cloudpart)
IF (psize.EQ.0) RETURN
p_110: DO icld=1,psize
   cloudpart(icld)%kdistype = 0
   cloudpart(icld)%label = 1
   cloudpart(icld)%nopart = 0
   cloudpart(icld)%state = 2
   cloudpart(icld)%length = 0.0
   cloudpart(icld)%orientation = 0.0
   cloudpart(icld)%thick = 0.0
   cloudpart(icld)%width = 0.0
ENDDO p_110


RETURN

END SUBROUTINE cloud_init

!========================================================================

SUBROUTINE outppart_init(partout)
!************************************************************************
!
! *outppart_init* Initialise derived type variable of type 'OutPartParams'
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)parttypes_init.f90  V2.11
!
! Description -
!
!************************************************************************
!
USE parttypes 
  
!
!*Arguments
!
TYPE (OutPartParams), INTENT(OUT), DIMENSION(:) :: partout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*partout*   DERIVED Attributes for particle output
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, psize


psize = SIZE(partout)
IF (psize.EQ.0) RETURN

n_110: DO n=1,psize
   partout(n)%aging = .FALSE.
   partout(n)%enddate = ''
   partout(n)%refdate = ''
   partout(n)%startdate = ''
   partout(n)%ktype = 0
   partout(n)%label = 1
   partout(n)%nodim = 0
   partout(n)%nopartout = 0
   partout(n)%nstepout = 0
   partout(n)%ntype = 0
   partout(n)%tlims = 0
   partout(n)%deltout = 0.0
ENDDO n_110


RETURN

END SUBROUTINE outppart_init

!========================================================================

SUBROUTINE parts_init(attspart)
!************************************************************************
!
! *parts_init* Initialise derived type variable of type 'PartAtts'
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)parttypes_init.f90  V2.11
!
! Description -
!
!************************************************************************
!
USE parttypes

!
!*Arguments
!
TYPE (PartAtts), INTENT(OUT), DIMENSION(:) :: attspart

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*attspart*  DERIVED Particle attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: p, psize


psize = SIZE(attspart)
IF (psize.EQ.0) RETURN

p_110: DO p=1,psize
   attspart(p)%center = .FALSE.
   attspart(p)%drift_state = 0
   attspart(p)%icoord = 0
   attspart(p)%jcoord = 0
   attspart(p)%kcoord = 1
   attspart(p)%kdistype = 0 
   attspart(p)%label = 1
   attspart(p)%ntstart = 0
   attspart(p)%state = 0
   attspart(p)%displx = 0.0
   attspart(p)%disply = 0.0
   attspart(p)%displz = 0.0
   attspart(p)%xcoord = 0.0
   attspart(p)%ycoord = 0.0
   attspart(p)%age = 0.0
   attspart(p)%tstart = 0.0
   attspart(p)%xpos = 0.0
   attspart(p)%ypos = 0.0
   attspart(p)%zpos = 0.0
ENDDO p_110


RETURN

END SUBROUTINE parts_init


END MODULE parttypes_init
