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


MODULE particle_routines
!************************************************************************
!
! *particle_routines* Utility routines for the  particle  model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_routines.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description -
!
! Reference -
!
! Generic routines - Cint_at_p, Uint_at_p, Vint_at_p
!  
! Routines - Cint_at_Z, Uint_at_Z, Vint_at_Z, Wint_at_p, Wint_at_Z,
!            vert_coord_index_part
!
!************************************************************************
!

IMPLICIT NONE
  
INTERFACE Cint_at_p
   MODULE PROCEDURE Cint_at_p_2d, Cint_at_p_3d
END INTERFACE

INTERFACE Uint_at_p
   MODULE PROCEDURE Uint_at_p_2d, Uint_at_p_3d
END INTERFACE

INTERFACE Vint_at_p
   MODULE PROCEDURE Vint_at_p_2d, Vint_at_p_3d
END INTERFACE

CONTAINS

!========================================================================

FUNCTION Cint_at_p_2d(xin,i,j,x,y)
!************************************************************************
!
! *Cint_at_p_2d* interpolate a 2-D C-node variable at a particle location with
!                relative coordinates (i,j,x,y)
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_routines.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description - only for uniform rectangular grids
!
! Module calls -
!
!************************************************************************
!
USE partvars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, j
REAL, INTENT(IN) :: x, y
REAL, INTENT(IN), DIMENSION(3,3) :: xin
REAL :: Cint_at_p_2d

!
! Name     Type    Purpose
!-----------------------------------------------------------------------------
!*xin*     REAL    Model data to be interpolated
!*i*       INTEGER X-index coordinate with respect to the lower corner of the
!                  grid cell where the particle resides 
!*j*       INTEGER Y-index coordinate with respect to the lower corner of the
!                  grid cell where the particle resides
!*x*       REAL    Normalised X-coordinate (between 0 and 1) with respect
!                  to the origin of the grid cell where the particle resides
!*y*       REAL    Normalised Y-coordinate (between 0 and 1) with respect
!                  to the origin of the grid cell where the particle resides
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: igrd, i1, jgrd, j1
REAL :: xnew, ynew, wsum
REAL, DIMENSION(2,2) :: weights


!---relative coordinates with respect to C-grid
i1 = MERGE(1,2,INT(x+0.5).EQ.0)
j1 = MERGE(1,2,INT(y+0.5).EQ.0)
xnew = MERGE(x+0.5,x-0.5,i1.EQ.1)
ynew = MERGE(y+0.5,y-0.5,j1.EQ.1)

!---weight factors
igrd = i1+i-2; jgrd = j1+j-2
weights(1,1) = MERGE((1-xnew)*(1-ynew),0.0,maskatc(igrd,jgrd))
weights(2,1) = MERGE(xnew*(1-ynew),0.0,maskatc(igrd+1,jgrd))
weights(1,2) = MERGE((1-xnew)*ynew,0.0,maskatc(igrd,jgrd+1))
weights(2,2) = MERGE(xnew*ynew,0.0,maskatc(igrd+1,jgrd+1))
wsum = SUM(weights)

!---interpolate
IF (wsum.GT.0.0) THEN
   weights = weights/wsum
   Cint_at_p_2d = SUM(weights*xin(i1:i1+1,j1:j1+1))/wsum
ELSE
   Cint_at_p_2d = 0.0
ENDIF


RETURN

END FUNCTION Cint_at_p_2d

!========================================================================

FUNCTION Cint_at_p_3d(xin,i,j,x,y,zpos)
!************************************************************************
!
! *Cint_at_p_3d* interpolate a 3-D C-node variable at a particle location with
!                relative coordinates (i,j,x,y) and vertical coordinate zpos
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_routines.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description - only for uniform rectangular grids
!
! Module calls - Cint_at_p_2d, Cint_at_Z
!
!************************************************************************
!
USE gridpars  
USE partvars
USE syspars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, j
REAL, INTENT(IN) :: x, y
REAL (KIND=kndrlong), INTENT(IN) :: zpos
REAL, INTENT(IN), DIMENSION(3,3,nz) :: xin
REAL :: Cint_at_p_3d

!
! Name     Type    Purpose
!-----------------------------------------------------------------------------
!*xin*     REAL    Model data to be interpolated
!*i*       INTEGER X-index coordinate with respect to the lower corner of the
!                  grid cell where the particle resides 
!*j*       INTEGER Y-index coordinate with respect to the lower corner of the
!                  grid cell where the particle resides
!*x*       REAL    Normalised X-coordinate (between 0 and 1) with respect
!                  to the origin of the grid cell where the particle resides
!*y*       REAL    Normalised Y-coordinate (between 0 and 1) with respect
!                  to the origin of the grid cell where the particle resides
!*zpos*    REAL    Vertical (Z-)coordinate at the particle location
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: igrd, ii, isub, i1, jgrd, jj, jsub, j1
REAL :: xnew, ynew, wsum
REAL, DIMENSION(2,2) :: xval, weights


!---2-D case
IF (SIZE(xin,DIM=3).EQ.1) THEN
   Cint_at_p_3d = Cint_at_p_2d(xin,i,j,x,y)
   RETURN
ENDIF

!---relative coordinates with respect to C-grid
i1 = MERGE(1,2,INT(x+0.5).EQ.0)
j1 = MERGE(1,2,INT(y+0.5).EQ.0)
xnew = MERGE(x+0.5,x-0.5,i1.EQ.1)
ynew = MERGE(y+0.5,y-0.5,j1.EQ.1)

!---weight factors
igrd = i1+i-2; jgrd = j1+j-2
weights(1,1) = MERGE((1-xnew)*(1-ynew),0.0,maskatc(igrd,jgrd))
weights(2,1) = MERGE(xnew*(1-ynew),0.0,maskatc(igrd+1,jgrd))
weights(1,2) = MERGE((1-xnew)*ynew,0.0,maskatc(igrd,jgrd+1))
weights(2,2) = MERGE(xnew*ynew,0.0,maskatc(igrd+1,jgrd+1))
wsum = SUM(weights)

!---interpolate model data at the particle's vertical location
ii_110: DO ii=i1,i1+1
jj_110: DO jj=j1,j1+1
   isub = ii-i1+1; jsub = jj-j1+1
   igrd = ii+i-2; jgrd = jj+j-2
   IF (weights(isub,jsub).GT.0.0) THEN
      xval(isub,jsub) = Cint_at_Z(xin(ii,jj,:),igrd,jgrd,zpos)
   ELSE
      xval(isub,jsub) = 0.0
   ENDIF
ENDDO jj_110
ENDDO ii_110

!---interpolate
IF (wsum.GT.0.0) THEN
   weights = weights/wsum
   Cint_at_p_3d = SUM(weights*xval)
ELSE
   Cint_at_p_3d = 0.0
ENDIF


RETURN

END FUNCTION Cint_at_p_3d

!========================================================================

FUNCTION Cint_at_Z(xin,i,j,zpos)
!************************************************************************
!
! *Cint_at_Z* interpolate a C-node variable vertically at the vertical
!             location zpos 
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_routines.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description - only for uniform rectangular grids
!
! Module calls -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE switches
USE syspars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, j
REAL (KIND=kndrlong), INTENT(IN) :: zpos
REAL, INTENT(IN), DIMENSION(nz) :: xin
REAL :: Cint_at_Z

!
! Name     Type    Purpose
!-----------------------------------------------------------------------------
!*xin*     REAL    C-node model vector to be interpolated
!*i*       INTEGER X-index coordinate of the C-node vertical vector
!*j*       INTEGER Y-index coordinate of the C-node vertical vector
!*zpos*    REAL    Vertical location in Z-coordinates
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k
REAL :: spos


!---convert location to s-coordinate
spos = MAX(0.0,MIN(1.0,(zpos+depmeanatc(i,j)/deptotatc(i,j))))

!---uniform sigma grid
IF (iopt_grid_vtype.EQ.1) THEN
   k = MIN(nz,INT(spos*nz+ 0.5))
   IF (k.LT.1) THEN
      Cint_at_Z = xin(1)
   ELSEIF (k.EQ.nz) THEN
      Cint_at_Z = xin(nz)
   ELSE
      Cint_at_Z = (k+0.5-spos)*xin(k)+(spos-k+0.5)*xin(k+1)
   ENDIF
   
!---non-uniform sigma grid   
ELSEIF (iopt_grid_vtype.EQ.2) THEN
   IF (spos.LE.gsigcoordatc(1)) THEN
      Cint_at_Z = xin(1)
   ELSEIF (spos.GE.gsigcoordatc(nz)) THEN
      Cint_at_Z = xin(nz)
   ELSE
      k = 1
      DO WHILE (k.LT.nz.AND.spos.GT.gsigcoordatc(k+1))
         k = k + 1
      ENDDO
      Cint_at_Z = ((gsigcoordatc(k+1)-spos)*xin(k)+&
                 & (spos-gsigcoordatc(k))*xin(k+1))&
                 & /(gsigcoordatc(k+1)-gsigcoordatc(k))
   ENDIF

!---non-uniform s-grid
ELSEIF (iopt_grid_vtype.EQ.3) THEN
   IF (spos.LE.gscoordatc(i,j,1)) THEN
      Cint_at_Z = xin(1)
   ELSEIF (spos.GE.gscoordatc(i,j,nz)) THEN
      Cint_at_Z = xin(nz)
   ELSE
      k = 1
      DO WHILE (k.LT.nz.AND.spos.GT.gscoordatc(i,j,k+1))
         k = k + 1
      ENDDO
      Cint_at_Z = ((gscoordatc(i,j,k+1)-spos)*xin(k)+&
                 & (spos-gscoordatc(i,j,k))*xin(k+1))&
                 & /(gscoordatc(i,j,k+1)-gscoordatc(i,j,k))
   ENDIF
ENDIF


RETURN

END FUNCTION Cint_at_Z

!========================================================================

FUNCTION Uint_at_p_2d(xin,i,j,x,y)
!************************************************************************
!
! *Uint_at_p_2d* interpolate a 2-D U-node variable at a particle location
!                with relative coordinates (i,j,x,y)
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_routines.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description - only for uniform rectangular grids
!
! Module calls -
!
!************************************************************************
!
USE partvars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, j
REAL, INTENT(IN) :: x, y
REAL, INTENT(IN), DIMENSION(2,3) :: xin
REAL :: Uint_at_p_2d

!
! Name     Type    Purpose
!-----------------------------------------------------------------------------
!*xin*     REAL    Model data to be interpolated
!*i*       INTEGER X-index coordinate with respect to the lower corner of the
!                  grid cell where the particle resides 
!*j*       INTEGER Y-index coordinate with respect to the lower corner of the
!                  grid cell where the particle resides
!*x*       REAL    Normalised X-coordinate (between 0 and 1) with respect
!                  to the origin of the grid cell where the particle resides
!*y*       REAL    Normalised Y-coordinate (between 0 and 1) with respect
!                  to the origin of the grid cell where the particle resides
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: jgrd, j1
REAL :: ynew, wsum
REAL, DIMENSION(2,2) :: weights


!---relative coordinates with respect to C-grid
j1 = MERGE(1,2,INT(y+0.5).EQ.0)
ynew = MERGE(y+0.5,y-0.5,j1.EQ.1)

!---weight factors
jgrd = j1+j-2
weights(1,1) = MERGE((1-x)*(1-ynew),0.0,maskatu(i,jgrd))
weights(2,1) = MERGE(x*(1-ynew),0.0,maskatu(i+1,jgrd))
weights(1,2) = MERGE((1-x)*ynew,0.0,maskatu(i,jgrd+1))
weights(2,2) = MERGE(x*ynew,0.0,maskatu(i+1,jgrd+1))
wsum = SUM(weights)

!---interpolate
IF (wsum.GT.0.0) THEN
   weights = weights/wsum
   Uint_at_p_2d = SUM(weights*xin(:,j1:j1+1))
ELSE
   Uint_at_p_2d = 0.0
ENDIF


RETURN

END FUNCTION Uint_at_p_2d

!========================================================================

FUNCTION Uint_at_p_3d(xin,i,j,x,y,zpos)
!************************************************************************
!
! *Uint_at_p_3d* interpolate a 3-D U-node variable at a particle location with
!                relative coordinates (i,j,x,y) and vertical coordinate zpos
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_routines.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description - only for uniforl rectangular grids
!
! Reference -
!
! Calling program -
!
! Module calls - Uint_at_p_2d, Uint_at_Z
!
!************************************************************************
!
USE gridpars
USE partvars
USE syspars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, j
REAL, INTENT(IN) :: x, y
REAL (KIND=kndrlong), INTENT(IN) :: zpos
REAL, INTENT(IN), DIMENSION(2,3,nz) :: xin
REAL :: Uint_at_p_3d

!
! Name     Type    Purpose
!-----------------------------------------------------------------------------
!*xin*     REAL    Model data to be interpolated
!*i*       INTEGER X-index coordinate with respect to the lower corner of the
!                  grid cell where the particle resides 
!*j*       INTEGER Y-index coordinate with respect to the lower corner of the
!                  grid cell where the particle resides
!*x*       REAL    Normalised X-coordinate (between 0 and 1) with respect
!                  to the origin of the grid cell where the particle resides
!*y*       REAL    Normalised Y-coordinate (between 0 and 1) with respect
!                  to the origin of the grid cell where the particle resides
!*zpos*    REAL    Vertical (Z-)coordinate at particle location
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ii, isub, jgrd, jj, jsub, j1
REAL :: ynew, wsum
REAL, DIMENSION(2,2) :: xval, weights


!---2-D case
IF (SIZE(xin,DIM=3).EQ.1) THEN
   Uint_at_p_3d = Uint_at_p_2d(xin,i,j,x,y)
   RETURN
ENDIF

!---relative coordinates with respect to C-grid
j1 = MERGE(1,2,INT(y+0.5).EQ.0)
ynew = MERGE(y+0.5,y-0.5,j1.EQ.1)

!---weight factors
jgrd = j1+j-2
weights(1,1) = MERGE((1-x)*(1-ynew),0.0,maskatu(i,jgrd))
weights(2,1) = MERGE(x*(1-ynew),0.0,maskatu(i+1,jgrd))
weights(1,2) = MERGE((1-x)*ynew,0.0,maskatu(i,jgrd+1))
weights(2,2) = MERGE(x*ynew,0.0,maskatu(i+1,jgrd+1))
wsum = SUM(weights)

!---interpolate input data at vertical location
ii_110: DO ii=i,i+1
jj_110: DO jj=j1,j1+1
   isub = ii-i+1; jsub = jj-j1+1
   jgrd = jj+j-2
   IF (weights(isub,jsub).GT.0.0) THEN
      xval(isub,jsub) = Uint_at_Z(xin(isub,jj,:),ii,jgrd,zpos)
   ELSE
      xval(isub,jsub) = 0.0
   ENDIF
ENDDO jj_110
ENDDO ii_110

!---interpolate
IF (wsum.GT.0.0) THEN
   weights = weights/wsum
   Uint_at_p_3d = SUM(weights*xval)
ELSE
   Uint_at_p_3d = 0.0
ENDIF


RETURN

END FUNCTION Uint_at_p_3d

!========================================================================

FUNCTION Uint_at_Z(xin,i,j,zpos)
!************************************************************************
!
! *Uint_at_Z* interpolate a U-node variable vertically at the vertical
!             location zpos 
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_routines.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description - only for uniform rectangular grids
!
! Module calls -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE switches
USE syspars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, j
REAL (KIND=kndrlong), INTENT(IN) :: zpos
REAL, INTENT(IN), DIMENSION(nz) :: xin
REAL :: Uint_at_Z

!
! Name     Type    Purpose
!-----------------------------------------------------------------------------
!*xin*     REAL    U-node model vector to be interpolated
!*i*       INTEGER X-index coordinate of the C-node vertical vector
!*j*       INTEGER Y-index coordinate of the C-node vertical vector
!*zpos*    REAL    Vertical location in Z-coordinates
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k
REAL :: spos


!---convert location to s-coordinate
spos = (zpos+depmeanatu(i,j))/deptotatu(i,j)

!---uniform sigma grid
IF (iopt_grid_vtype.EQ.1) THEN
   k = MIN(nz,INT(spos*nz+0.5))
   IF (k.LT.1) THEN
      Uint_at_Z = xin(1)
   ELSEIF (k.EQ.nz) THEN
      Uint_at_Z = xin(nz)
   ELSE
      Uint_at_Z = (k+0.5-spos)*xin(k)+(spos-k+0.5)*xin(k+1)
   ENDIF
   
!---non-uniform sigma grid   
ELSEIF (iopt_grid_vtype.EQ.2) THEN
   IF (spos.LE.gsigcoordatc(1)) THEN
      Uint_at_Z = xin(1)
   ELSEIF (spos.GE.gsigcoordatc(nz)) THEN
      Uint_at_Z = xin(nz)
   ELSE
      k = 1
      DO WHILE (k.LT.nz.AND.spos.GT.gsigcoordatc(k+1))
         k = k + 1
      ENDDO
      Uint_at_Z = ((gsigcoordatc(k+1)-spos)*xin(k)+&
                 & (spos-gsigcoordatc(k))*xin(k+1))&
                 & /(gsigcoordatc(k+1)-gsigcoordatc(k))
   ENDIF

!---non-uniform s-grid
ELSEIF (iopt_grid_vtype.EQ.3) THEN
   IF (spos.LE.gscoordatu(i,j,1)) THEN
      Uint_at_Z = xin(1)
   ELSEIF (spos.GE.gscoordatu(i,j,nz)) THEN
      Uint_at_Z = xin(nz)
   ELSE
      k = 1
      DO WHILE (k.LT.nz.AND.spos.GT.gscoordatu(i,j,k+1))
         k = k + 1
      ENDDO
      Uint_at_Z = ((gscoordatu(i,j,k+1)-spos)*xin(k)+&
                 & (spos-gscoordatu(i,j,k))*xin(k+1))&
                 & /(gscoordatu(i,j,k+1)-gscoordatu(i,j,k))
   ENDIF
ENDIF


RETURN

END FUNCTION Uint_at_Z

!========================================================================

FUNCTION Vint_at_p_2d(xin,i,j,x,y)
!************************************************************************
!
! *Vint_at_p_2d* interpolate a 2-D V-node variable at a particle location
!                with relative coordinates (i,j,x,y)
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_routines.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description - only for uniform rectangular grids
!
! Module calls -
!
!************************************************************************
!
USE partvars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, j
REAL, INTENT(IN) :: x, y
REAL, INTENT(IN), DIMENSION(3,2) :: xin
REAL :: Vint_at_p_2d

!
! Name     Type    Purpose
!-----------------------------------------------------------------------------
!*xin*     REAL    Model data to be interpolated
!*i*       INTEGER X-index coordinate with respect to the lower corner of the
!                  grid cell where the particle resides 
!*j*       INTEGER Y-index coordinate with respect to the lower corner of the
!                  grid cell where the particle resides
!*x*       REAL    Normalised X-coordinate (between 0 and 1) with respect
!                  to the origin of the grid cell where the particle resides
!*y*       REAL    Normalised Y-coordinate (between 0 and 1) with respect
!                  to the origin of the grid cell where the particle resides
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: igrd, i1
REAL :: xnew, wsum
REAL, DIMENSION(2,2) :: weights


!---relative coordinates with respect to C-grid
i1 = MERGE(1,2,INT(x+0.5).EQ.0)
xnew = MERGE(x+0.5,x-0.5,i1.EQ.1)

!---weight factors
igrd = MERGE(i-1,i,i1.EQ.1)
weights(1,1) = MERGE((1-xnew)*(1-y),0.0,maskatv(igrd,j))
weights(2,1) = MERGE(xnew*(1-y),0.0,maskatv(igrd+1,j))
weights(1,2) = MERGE((1-xnew)*y,0.0,maskatv(igrd,j+1))
weights(2,2) = MERGE(xnew*y,0.0,maskatv(igrd+1,j+1))
wsum = SUM(weights)

!---interpolate
IF (wsum.GT.0) THEN
   weights = weights/wsum
   Vint_at_p_2d = SUM(weights*xin(i1:i1+1,:))
ELSE
   Vint_at_p_2d = 0.0
ENDIF


RETURN

END FUNCTION Vint_at_p_2d

!========================================================================

FUNCTION Vint_at_p_3d(xin,i,j,x,y,zpos)
!************************************************************************
!
! *Vint_at_p_3d* interpolate a 3-D V-node variable at a particle location with
!                relative coordinates (i,j,x,y) and vertical coordinate zpos
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_routines.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description - only for uniform rectangular grids
!
! Reference -
!
! Calling program -
!
! Module calls - Vint_at_p_2d, Vint_at_Z
!
!************************************************************************
!
USE gridpars
USE partvars
USE syspars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, j
REAL, INTENT(IN) :: x, y
REAL (KIND=kndrlong), INTENT(IN) :: zpos
REAL, INTENT(IN), DIMENSION(3,2,nz) :: xin
REAL :: Vint_at_p_3d

!
! Name     Type    Purpose
!-----------------------------------------------------------------------------
!*xin*     REAL    Model data to be interpolated
!*i*       INTEGER X-index coordinate with respect to the lower corner of the
!                  grid cell where the particle resides 
!*j*       INTEGER Y-index coordinate with respect to the lower corner of the
!                  grid cell where the particle resides
!*x*       REAL    Normalised X-coordinate (between 0 and 1) with respect
!                  to the origin of the grid cell where the particle resides
!*y*       REAL    Normalised Y-coordinate (between 0 and 1) with respect
!                  to the origin of the grid cell where the particle resides
!*zpos*    REAL    Vertical (Z-)coordinate at particle location
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: igrd, ii, isub, i1, jj, jsub
REAL :: xnew, wsum
REAL, DIMENSION(2,2) :: xval, weights


!---2-D case
IF (SIZE(xin,DIM=3).EQ.1) THEN
   Vint_at_p_3d = Vint_at_p_2d(xin,i,j,x,y)
   RETURN
ENDIF

!---relative coordinates with respect to C-grid
i1 = MERGE(1,2,INT(x+0.5).EQ.0)
xnew = MERGE(x+0.5,x-0.5,i1.EQ.1)

!---weight factors
igrd = i1+i-2
weights(1,1) = MERGE((1-xnew)*(1-y),0.0,maskatv(igrd,j))
weights(2,1) = MERGE(xnew*(1-y),0.0,maskatv(igrd+1,j))
weights(1,2) = MERGE((1-xnew)*y,0.0,maskatv(igrd,j+1))
weights(2,2) = MERGE(xnew*y,0.0,maskatv(igrd+1,j+1))
wsum = SUM(weights)

!---interpolate input data at vertical location
ii_110: DO ii=i1,i1+1
jj_110: DO jj=j,j+1
   isub = ii-i1+1; jsub = jj-j+1
   igrd = ii+i-2
   IF (weights(isub,jsub).GT.0.0) THEN
      xval(isub,jsub) = Vint_at_Z(xin(ii,jsub,:),igrd,jj,zpos)
   ELSE
      xval(isub,jsub) = 0.0
   ENDIF
ENDDO jj_110
ENDDO ii_110

!---interpolate
IF (wsum.GT.0.0) THEN
   weights = weights/wsum
   Vint_at_p_3d = SUM(weights*xval)
ELSE
   Vint_at_p_3d = 0.0
ENDIF


RETURN

END FUNCTION Vint_at_p_3d

!========================================================================

FUNCTION Vint_at_Z(xin,i,j,zpos)
!************************************************************************
!
! *Vint_at_Z* interpolate a V-node variable vertically at the vertical
!             location zpos 
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_routines.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description - only for uniform rectangular grids
!
! Module calls -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE switches
USE syspars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, j
REAL (KIND=kndrlong), INTENT(IN) :: zpos
REAL, INTENT(IN), DIMENSION(nz) :: xin
REAL :: Vint_at_Z

!
! Name     Type    Purpose
!-----------------------------------------------------------------------------
!*xin*     REAL    V-node model vector to be interpolated
!*i*       INTEGER X-index coordinate of the C-node vertical vector
!*j*       INTEGER Y-index coordinate of the C-node vertical vector
!*zpos*    REAL    Vertical location in Z-coordinates
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k
REAL :: spos


!---convert location to s-coordinate
spos = (zpos+depmeanatv(i,j))/deptotatv(i,j)

!---uniform sigma grid
IF (iopt_grid_vtype.EQ.1) THEN
   k = MIN(nz,INT(spos*nz+0.5))
   IF (k.LT.1) THEN
      Vint_at_Z = xin(1)
   ELSEIF (k.EQ.nz) THEN
      Vint_at_Z = xin(nz)
   ELSE
      Vint_at_Z = (k+0.5-spos)*xin(k)+(spos-k+0.5)*xin(k+1)
   ENDIF
   
!---non-uniform sigma grid   
ELSEIF (iopt_grid_vtype.EQ.2) THEN
   IF (spos.LE.gsigcoordatc(1)) THEN
      Vint_at_Z = xin(1)
   ELSEIF (spos.GE.gsigcoordatc(nz)) THEN
      Vint_at_Z = xin(nz)
   ELSE
      k = 1
      DO WHILE (k.LT.nz.AND.spos.GT.gsigcoordatc(k+1))
         k = k + 1
      ENDDO
      Vint_at_Z = ((gsigcoordatc(k+1)-spos)*xin(k)+&
                & (spos-gsigcoordatc(k))*xin(k+1)) &
                & /(gsigcoordatc(k+1)-gsigcoordatc(k))
   ENDIF

!---non-uniform s-grid
ELSEIF (iopt_grid_vtype.EQ.3) THEN
   IF (spos.LE.gscoordatv(i,j,1)) THEN
      Vint_at_Z = xin(1)
   ELSEIF (spos.GE.gscoordatv(i,j,nz)) THEN
      Vint_at_Z = xin(nz)
   ELSE
      k = 1
      DO WHILE (k.LT.nz.AND.spos.GT.gscoordatv(i,j,k+1))
         k = k + 1
      ENDDO
      Vint_at_Z = ((gscoordatv(i,j,k+1)-spos)*xin(k)+&
                 & (spos-gscoordatv(i,j,k))*xin(k+1))&
                 & /(gscoordatv(i,j,k+1)-gscoordatv(i,j,k))
   ENDIF
ENDIF


RETURN

END FUNCTION Vint_at_Z

!========================================================================

FUNCTION Wint_at_p(xin,i,j,x,y,zpos)
!************************************************************************
!
! *Wint_at_p_3d* interpolate a 3-D W-node variable at a particle location with
!                relative coordinates (i,j,x,y) and vertical coordinate zpos
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_routines.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description - only for uniform rectangular grids
!
! Module calls - Wint_at_Z
!
!************************************************************************
!
USE gridpars
use iopars
USE partvars  
USE syspars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, j
REAL, INTENT(IN) :: x, y
REAL (KIND=kndrlong), INTENT(IN) :: zpos
REAL, INTENT(IN) :: xin(3,3,nz+1)
REAL :: Wint_at_p

!
! Name     Type    Purpose
!-----------------------------------------------------------------------------
!*xin*     REAL    Model data to be interpolated
!*i*       INTEGER X-index coordinate with respect to the lower corner of the
!                  grid cell where the particle resides 
!*j*       INTEGER Y-index coordinate with respect to the lower corner of the
!                  grid cell where the particle resides
!*x*       REAL    Normalised X-coordinate (between 0 and 1) with respect
!                  to the origin of the grid cell where the particle resides
!*y*       REAL    Normalised Y-coordinate (between 0 and 1) with respect
!                  to the origin of the grid cell where the particle resides
!*zpos*    REAL    Vertical (Z-)coordinate at the particle location
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: igrd, ii, isub, i1, jgrd, jj, jsub, j1
REAL :: xnew, ynew, wsum
REAL, DIMENSION(2,2) :: xval, weights


!---relative coordinates with respect to C-grid
i1 = MERGE(1,2,INT(x+0.5).EQ.0)
j1 = MERGE(1,2,INT(y+0.5).EQ.0)
xnew = MERGE(x+0.5,x-0.5,i1.EQ.1)
ynew = MERGE(y+0.5,y-0.5,j1.EQ.1)

!---weight factors
igrd = i1+i-2; jgrd = j1+j-2
weights(1,1) = MERGE((1-xnew)*(1-ynew),0.0,maskatc(igrd,jgrd))
weights(2,1) = MERGE(xnew*(1-ynew),0.0,maskatc(igrd+1,jgrd))
weights(1,2) = MERGE((1-xnew)*ynew,0.0,maskatc(igrd,jgrd+1))
weights(2,2) = MERGE(xnew*ynew,0.0,maskatc(igrd+1,jgrd+1))
wsum = SUM(weights)

!---interpolate model data at the particle's vertical location
ii_110: DO ii=i1,i1+1
jj_110: DO jj=j1,j1+1
   isub = ii-i1+1; jsub = jj-j1+1
   igrd = ii+i-2; jgrd = jj+j-2
   IF (weights(isub,jsub).GT.0.0) THEN
      xval(isub,jsub) = Wint_at_Z(xin(ii,jj,:),igrd,jgrd,zpos)
   ELSE
      xval(isub,jsub) = 0.0
   ENDIF
ENDDO jj_110
ENDDO ii_110

!---interpolate
IF (wsum.GT.0) THEN
   weights = weights/wsum
   Wint_at_p = SUM(weights*xval)
ELSE
   Wint_at_p = 0.0
ENDIF


RETURN

END FUNCTION Wint_at_p

!========================================================================

FUNCTION Wint_at_Z(xin,i,j,zpos)
!************************************************************************
!
! *Wint_at_Z* interpolate a W-node variable vertically at the vertical
!             location zpos 
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_routines.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description - only for uniform rectangular grids
!
! Module calls -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE switches
USE syspars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, j
REAL (KIND=kndrlong), INTENT(IN) :: zpos
REAL, INTENT(IN), DIMENSION(nz+1) :: xin
REAL :: Wint_at_Z

!
! Name     Type    Purpose
!-----------------------------------------------------------------------------
!*xin*     REAL    C-node model vector to be interpolated
!*i*       INTEGER X-index coordinate of the C-node vertical vector
!*j*       INTEGER Y-index coordinate of the C-node vertical vector
!*zpos*    REAL    Vertical location in Z-coordinates
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k
REAL :: s, spos


!---convert location to s-coordinate
spos = (zpos+depmeanatc(i,j))/deptotatc_old(i,j)
IF (spos.GE.1.0) THEN
   Wint_at_Z = xin(nz+1)
   RETURN
ELSEIF (spos.LE.0.0) THEN
   Wint_at_Z = xin(1)
   RETURN
ENDIF

!---uniform sigma grid
IF (iopt_grid_vtype.EQ.1) THEN
   s = spos*nz
   k = MIN(nz,MAX(1,INT(s+1.0)))
   Wint_at_Z = (k-s)*xin(k)+(s-k+1)*xin(k+1)
   
!---non-uniform sigma grid   
ELSEIF (iopt_grid_vtype.EQ.2) THEN
   k = 1
   DO WHILE (k.LE.nz.AND.spos.GT.gsigcoordatw(k+1))
      k = k + 1
   ENDDO
   Wint_at_Z = ((gsigcoordatw(k+1)-spos)*xin(k)+&
              & (spos-gsigcoordatw(k))*xin(k+1))&
              & /(gsigcoordatw(k+1)-gsigcoordatw(k))

!---non-uniform s-grid
ELSEIF (iopt_grid_vtype.EQ.3) THEN
   k = 1
   DO WHILE (k.LE.nz.AND.spos.GT.gscoordatw(i,j,k+1))
      k = k + 1
   ENDDO
   Wint_at_Z = ((gscoordatw(i,j,k+1)-spos)*xin(k)+&
              & (spos-gscoordatw(i,j,k))*xin(k+1))&
              & /(gscoordatw(i,j,k+1)-gscoordatw(i,j,k))
ENDIF


RETURN

END FUNCTION Wint_at_Z

!========================================================================

SUBROUTINE vert_coord_index_part(icoord,jcoord,kcoord,xcoord,ycoord,zpos)
!************************************************************************
!
! *vert_coord_index_part* vertical grid index at the particle's location 
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_routines.f90  V2.11
!
! $Date: 2017-05-03 14:09:15 +0200 (Wed, 03 May 2017) $
!
! $Revision: 1020 $
!
! Description -
!
! Reference -
!
! Calling program -
!
! Module calls - Cint_at_p
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE switches
USE syspars

!
!*Arguments
!
INTEGER, INTENT(IN) :: icoord, jcoord
INTEGER, INTENT(OUT) :: kcoord
REAL, INTENT(IN) :: xcoord, ycoord
REAL (KIND=kndrlong), INTENT(IN) :: zpos

!
! Name         Type    Purpose
!----------------------------------------------------------------------------
!*icoord*      INTEGER X-index coordinate with respect to lower left
!                      corner at the new location
!*jcoord*      INTEGER Y-index coordinate with respect to lower left
!                      corner at the new location
!*kcoord*      INTEGER Returned vertical grid index
!*xcoord*      REAL    Normalised X-coordinate (between 0 and 1) with respect
!                      to the origin of the grid cell at the new location
!*ycoord*      REAL    Normalised Y-coordinate (between 0 and 1) with respect
!                      to the origin of the grid cellat the new location
!*zpos*        REAL    Z-coordinate at the particle's location
!
!----------------------------------------------------------------------------
!
!*Local variables
!
REAL :: depmeanp, deptotp, s, sposp


!---convert z- to sigma coordinate
deptotp = Cint_at_p(deptotatc(icoord-1:icoord+1,jcoord-1:jcoord+1),&
                  & icoord,jcoord,xcoord,ycoord)
depmeanp = Cint_at_p(depmeanatc(icoord-1:icoord+1,jcoord-1:jcoord+1),&
                   & icoord,jcoord,xcoord,ycoord)
sposp = MAX(0.0,MIN(1.0,(zpos+depmeanp)/deptotp))

!---uniform sigma grid
IF (iopt_grid_vtype.EQ.1) THEN
   s = sposp*nz
   kcoord = MIN(nz,MAX(1,INT(s+1.0)))
   
!---non-uniform sigma grid   
ELSEIF (iopt_grid_vtype.EQ.2) THEN
   kcoord = 1
   DO WHILE (kcoord.LE.nz.AND.sposp.GT.gsigcoordatw(kcoord+1))
      kcoord = kcoord + 1
   ENDDO
   kcoord = MIN(nz,kcoord)
   
!---non-uniform s-grid
ELSEIF (iopt_grid_vtype.EQ.3) THEN
   kcoord = 1
   DO WHILE (kcoord.LE.nz.AND.sposp.GT.gscoordatw(icoord,jcoord,kcoord+1))
      kcoord = kcoord + 1
   ENDDO
   kcoord = MIN(nz,kcoord)
ENDIF


RETURN

END SUBROUTINE vert_coord_index_part

END MODULE particle_routines
