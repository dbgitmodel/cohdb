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

MODULE grid_routines
!************************************************************************
!
! *grid_routines* Utility routines related to the model grid 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.11
!
! $Date: 2018-05-31 16:41:05 +0200 (Thu, 31 May 2018) $
!
! $Revision: 1141 $
!
! Generic routines - rotate_vec
! 
! Routines - check_cell_location, construct_rectgrid, convert_loc_to_char,
!            distance_minloc, distance_pts, find_obc_index_glb, global_mask,
!            grid_rotation_params, grid_rotation_transf, grid_spacings_curv,
!            grid_spacings_rectang, local_proc, mask_array, num_proc,
!            Zcoord_arr, Zcoord_var
!
!************************************************************************
!
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

INTERFACE check_cell_location
   MODULE PROCEDURE check_cell_location_s, check_cell_location_d
END INTERFACE check_cell_location

INTERFACE construct_rectgrid
   MODULE PROCEDURE construct_rectgrid_nonunif, construct_rectgrid_unif
END INTERFACE construct_rectgrid

INTERFACE distance_minloc
   MODULE PROCEDURE distance_minloc_s, distance_minloc_d
END INTERFACE distance_minloc

INTERFACE distance_pts
   MODULE PROCEDURE distance_pts_s, distance_pts_d
END INTERFACE distance_pts

INTERFACE rotate_vec
   MODULE PROCEDURE rotate_vec_0d, rotate_vec_1d, rotate_vec_2d, &
                  & rotate_vec_3d, rotate_vec_4d
END INTERFACE rotate_vec

CONTAINS

!========================================================================

SUBROUTINE check_cell_location_s(xpos,ypos,xvals,yvals,flag)
!************************************************************************
!
! *check_cell_location_s* check whether a location is situated within a given
!                         quadrilateral cell 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.11
!
! Description - xpos and ypos given in single precision
!
!************************************************************************
!
USE syspars
!
!*Arguments
!
LOGICAL, INTENT(OUT) :: flag
REAL (KIND=kndreal), INTENT(IN) :: xpos, ypos
REAL, INTENT(IN), DIMENSION(4) :: xvals, yvals

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*xpos*     REAL    X-coordinate of the input location [m or degrees longitude]
!*ypos*     REAL    Y-coordinate of the input location [m or degrees latitude]
!*xvals*    REAL    X-coordinates at the corners (SW,SE,NE,NW) of the
!                   quadrilateral cell                  [m or degrees longitude]
!*yvals*    REAL    Y-coordinates at the corners (SW,SE,NE,NW) of the
!                   quadrilateral cell                   [m or degrees latitude]
!*flag*     LOGICAL .TRUE./.FALSE. if location is inside/outside the cell
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, next_n
REAL :: crossprod, crossprod_old, vec1x, vec1y, vec2x, vec2y


flag = .TRUE.
n_100: DO n=1,4
   IF (flag) THEN
      next_n = MOD(n,4) + 1
      vec1x = xvals(next_n) - xvals(n)
      vec1y = yvals(next_n) - yvals(n)
      vec2x = xpos - xvals(n)
      vec2y = ypos - yvals(n)
      crossprod = vec1x*vec2y - vec2x*vec1y
      IF (n.EQ.1) crossprod_old = crossprod
      IF (crossprod*crossprod_old.LT.0.0) flag = .FALSE.
      crossprod_old = crossprod
   ENDIF
ENDDO n_100


RETURN

END SUBROUTINE check_cell_location_s

!========================================================================

SUBROUTINE check_cell_location_d(xpos,ypos,xvals,yvals,flag)
!************************************************************************
!
! *check_cell_location_d* check whether a location is situated within a given
!                         quadrilateral cell 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.11
!
! Description - xpos and ypos given in double precision
!
!************************************************************************
!
USE syspars
!
!*Arguments
!
LOGICAL, INTENT(OUT) :: flag
REAL (KIND=kndrlong), INTENT(IN) :: xpos, ypos
REAL, INTENT(IN), DIMENSION(4) :: xvals, yvals

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*xpos*     REAL    X-coordinate of the input location [m or degrees longitude]
!*ypos*     REAL    Y-coordinate of the input location [m or degrees latitude]
!*xvals*    REAL    X-coordinates at the corners (SW,SE,NE,NW) of the
!                   quadrilateral cell                  [m or degrees longitude]
!*yvals*    REAL    Y-coordinates at the corners (SW,SE,NE,NW) of the
!                   quadrilateral cell                   [m or degrees latitude]
!*flag*     LOGICAL .TRUE./.FALSE. if location is inside/outside the cell
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, next_n
REAL :: crossprod, crossprod_old, vec1x, vec1y, vec2x, vec2y


flag = .TRUE.
n_100: DO n=1,4
   IF (flag) THEN
      next_n = MOD(n,4) + 1
      vec1x = xvals(next_n) - xvals(n)
      vec1y = yvals(next_n) - yvals(n)
      vec2x = xpos - xvals(n)
      vec2y = ypos - yvals(n)
      crossprod = vec1x*vec2y - vec2x*vec1y
      IF (n.EQ.1) crossprod_old = crossprod
      IF (crossprod*crossprod_old.LT.0.0) flag = .FALSE.
      crossprod_old = crossprod
   ENDIF
ENDDO n_100


RETURN

END SUBROUTINE check_cell_location_d

!========================================================================

SUBROUTINE construct_rectgrid_nonunif(xstart,ystart,delx,dely,xcoord,ycoord,&
                                    & nxgrd,nygrd)
!************************************************************************
!
! *construct_rectgrid_nonunif* Construct non-uniform rectangular grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.5.1
!
! Description - 
!
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: nxgrd, nygrd
REAL, INTENT(IN) :: xstart, ystart
REAL, INTENT(IN), DIMENSION(nxgrd) :: delx
REAL, INTENT(IN), DIMENSION(nygrd) :: dely
REAL, INTENT(OUT), DIMENSION(nxgrd,nygrd) :: xcoord, ycoord

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*xstart*   REAL    X-coordinate of grid location (1,1)
!*ystart*   REAL    Y-coordinate of grid location (1,1)
!*delx*     REAL    Grid spacings in the X-direction
!*dely*     REAL    Grid spacings in the Y-direction
!*xcoord*   REAL    X-coordinates of cell corners
!*ycoord*   REAL    Y-coordinates of cell corners
!*nxgrd*    INTEGER Number of cell corners in X-direction
!*nygrd*    INTEGER Number of cell corners in Y-direction
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j


procname(pglev+1) = 'construct_rectgrid_nonunif'
CALL log_timer_in()

xcoord(1,:) = xstart
i_110: DO i=2,nxgrd
   xcoord(i,:) = xstart + SUM(delx(1:i-1))
ENDDO i_110

ycoord(:,1) = ystart
j_120: DO j=2,nygrd
   ycoord(:,j) = ystart + SUM(dely(1:j-1))
ENDDO j_120

CALL log_timer_out()


RETURN

END SUBROUTINE construct_rectgrid_nonunif

!========================================================================

SUBROUTINE construct_rectgrid_unif(xstart,ystart,delx,dely,xcoord,ycoord,&
                                 & nxgrd,nygrd)
!************************************************************************
!
! *construct_rectgrid_unif* Construct grid uniform rectangular grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.5.1
!
! Description - 
!
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: nxgrd, nygrd
REAL, INTENT(IN) :: xstart, ystart
REAL, INTENT(IN) :: delx, dely
REAL, INTENT(OUT), DIMENSION(nxgrd,nygrd) :: xcoord, ycoord

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*xstart*   REAL    X-coordinate of grid location (1,1)
!*ystart*   REAL    Y-coordinate of grid location (1,1)
!*delx*     REAL    Grid spacings in the X-direction
!*dely*     REAL    Grid spacings in the Y-direction
!*xcoord*   REAL    X-coordinates of cell corners
!*ycoord*   REAL    Y-coordinates of cell corners
!*nxgrd*    INTEGER Number of cell corners in X-direction
!*nygrd*    INTEGER Number of cell corners in Y-direction
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j


procname(pglev+1) = 'construct_rectgrid_unif'
CALL log_timer_in()

i_110: DO i=1,nxgrd
   xcoord(i,:) = xstart + (i-1)*delx
ENDDO i_110

j_120: DO j=1,nygrd
   ycoord(:,j) = ystart + (j-1)*dely
ENDDO j_120

CALL log_timer_out()


RETURN

END SUBROUTINE construct_rectgrid_unif

!========================================================================

FUNCTION convert_loc_to_char(xcoord,ycoord) RESULT(charloc)
!************************************************************************
!
! *convert_loc_to_char* Write geographical location into a string
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.0
!
! Description -
!
!************************************************************************
!
!*Arguments
!
REAL, INTENT(IN) :: xcoord, ycoord
CHARACTER (LEN=21) :: charloc

!
! Name     Type  Purpose
!------------------------------------------------------------------------------
!*xcoord*  REAL  X-coordinate of location
!*ycoord*  REAL  Y-coordinate of location 
!*charloc* CHAR  Geographical position in string format
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=1) :: clat, clon
INTEGER :: latdeg, latmin, latsec, londeg, lonmin, lonsec
REAL :: rmin


clat = MERGE('N','S',ycoord.GE.0.0)
latdeg = ABS(ycoord)
rmin = (ABS(ycoord)-latdeg)*60.0
latmin = rmin
latsec = (rmin-latmin)*60

clon = MERGE('E','W',xcoord.GE.0.0)
londeg = ABS(xcoord)
rmin = (ABS(xcoord)-londeg)*60.0
lonmin = rmin
lonsec = (rmin-lonmin)*60

WRITE (charloc,9001) latdeg, latmin, latsec, clat, '; ',&
                   & londeg, lonmin, lonsec, clon


RETURN

9001 FORMAT (I2,2(1X,I2.2),2A,I3,2(1X,I2.2),A)

END FUNCTION convert_loc_to_char

!========================================================================

SUBROUTINE distance_minloc_s(x,y,xcoord,ycoord,nx,ny,inunit,outunit,distmin,&
                           & locmin,mask)
!************************************************************************
!
! *distance_minloc_s* Find closest distance to or closest location of an input
!                     location to a series of input locations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.10.2
!
! Description - xpos and ypos given in single precision
!
! Module calls - complex_polar
!
!************************************************************************
!
USE physpars
USE switches
USE syspars
USE math_library, ONLY: complex_polar

!
!*Arguments
!
INTEGER, INTENT(IN) :: inunit, nx, ny, outunit
INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(2) :: locmin
REAL (KIND=kndreal), INTENT(IN) :: x, y
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask 
REAL, INTENT(IN), DIMENSION(nx,ny)  :: xcoord, ycoord
REAL, INTENT(OUT), OPTIONAL :: distmin

!
! Name     Type   Purpose
!------------------------------------------------------------------------------
!*x*       REAL    X-coordinate of input location [m or degrees longitude]
!*y*       REAL    Y-coordinate of input location [m or degrees latitude]
!*xcoord*  REAL    X-coordinates of data grid [m or degrees longitude]
!*ycoord*  REAL    Y-coordinates of data grid [m or degrees latitude]
!*nx*      INTEGER X-dimension of coordinate arrays
!*ny*      INTEGER Y-dimension of coordinate arrays
!*inunit*  INTEGER unit of input locations (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!*outunit* INTEGER unit of result (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!                   = 3 => meters
!*distmin* REAL    If PRESENT, minimum distance
!*locmin*  INTEGER If PRESENT, location of the closest distance
!*mask*    LOGICAL Excludes masked (dry) grid points if PRESENT and .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(nx,ny) :: dist, dx, dy

!
!1. Distance to input location
!-----------------------------
!

IF (iopt_grid_sph.EQ.0) THEN
   dx = xcoord-x; dy = ycoord-y
   CALL complex_polar(dx,dy,xamp=dist)
ELSEIF (iopt_grid_sph.EQ.1) THEN
   IF (inunit.EQ.1) THEN
      dx = 0.5*(xcoord-x); dy = 0.5*(ycoord-y)
      dist = 2.0*ASIN(SQRT(COS(y)*COS(ycoord)*SIN(dx)**2+SIN(dy)**2))
   ELSEIF (inunit.EQ.2) THEN
      dx = 0.5*degtorad*(xcoord-x); dy = 0.5*degtorad*(ycoord-y)
      dist = 2.0*ASIN(SQRT(COS(degtorad*y)*COS(degtorad*ycoord)*SIN(dx)**2+&
                         & SIN(dy)**2))
   ENDIF
ENDIF

!
!2. Return minimum distance and/or location
!------------------------------------------
!

IF (.NOT.PRESENT(mask)) THEN
   IF (PRESENT(distmin)) distmin = MINVAL(dist)
   IF (PRESENT(locmin)) locmin = MINLOC(dist)
ELSE
   IF (PRESENT(distmin)) distmin = MINVAL(dist,MASK=mask)
   IF (PRESENT(locmin)) locmin = MINLOC(dist,MASK=mask)
ENDIF
IF (PRESENT(distmin).AND.iopt_grid_sph.EQ.1) THEN
   IF (outunit.EQ.2) THEN
      distmin = radtodeg*distmin
   ELSEIF (outunit.EQ.3) THEN
      distmin = Rearth*distmin
   ENDIF
ENDIF


RETURN

END SUBROUTINE distance_minloc_s

!========================================================================

SUBROUTINE distance_minloc_d(x,y,xcoord,ycoord,nx,ny,inunit,outunit,distmin,&
                           & locmin,mask)
!************************************************************************
!
! *distance_minloc_d* Find closest distance to or closest location of an input
!                     location to a series of input locations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.10.2
!
! Description - xpos and ypos given in double precision
!
! Module calls - complex_polar
!
!************************************************************************
!
USE physpars
USE switches
USE syspars
USE math_library, ONLY: complex_polar

!
!*Arguments
!
INTEGER, INTENT(IN) :: inunit, nx, ny, outunit
INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(2) :: locmin
REAL (KIND=kndrlong), INTENT(IN) :: x, y
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask 
REAL, INTENT(IN), DIMENSION(nx,ny)  :: xcoord, ycoord
REAL, INTENT(OUT), OPTIONAL :: distmin

!
! Name     Type   Purpose
!------------------------------------------------------------------------------
!*x*       REAL    X-coordinate of input location [m or degrees longitude]
!*y*       REAL    Y-coordinate of input location [m or degrees latitude]
!*xcoord*  REAL    X-coordinates of data grid [m or degrees longitude]
!*ycoord*  REAL    Y-coordinates of data grid [m or degrees latitude]
!*nx*      INTEGER X-dimension of coordinate arrays
!*ny*      INTEGER Y-dimension of coordinate arrays
!*inunit*  INTEGER unit of input locations (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!*outunit* INTEGER unit of result (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!                   = 3 => meters
!*distmin* REAL    If PRESENT, minimum distance
!*locmin*  INTEGER If PRESENT, location of the closest distance
!*mask*    LOGICAL Excludes masked (dry) grid points if PRESENT and .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(nx,ny) :: dist, dx, dy

!
!1. Distance to input location
!-----------------------------
!

IF (iopt_grid_sph.EQ.0) THEN
   dx = xcoord-x; dy = ycoord-y
   CALL complex_polar(dx,dy,xamp=dist)
ELSEIF (iopt_grid_sph.EQ.1) THEN
   IF (inunit.EQ.1) THEN
      dx = 0.5*(xcoord-x); dy = 0.5*(ycoord-y)
      dist = 2.0*ASIN(SQRT(COS(y)*COS(ycoord)*SIN(dx)**2+SIN(dy)**2))
   ELSEIF (inunit.EQ.2) THEN
      dx = 0.5*degtorad*(xcoord-x); dy = 0.5*degtorad*(ycoord-y)
      dist = 2.0*ASIN(SQRT(COS(degtorad*y)*COS(degtorad*ycoord)*SIN(dx)**2+&
                         & SIN(dy)**2))
   ENDIF
ENDIF

!
!2. Return minimum distance and/or location
!------------------------------------------
!

IF (.NOT.PRESENT(mask)) THEN
   IF (PRESENT(distmin)) distmin = MINVAL(dist)
   IF (PRESENT(locmin)) locmin = MINLOC(dist)
ELSE
   IF (PRESENT(distmin)) distmin = MINVAL(dist,MASK=mask)
   IF (PRESENT(locmin)) locmin = MINLOC(dist,MASK=mask)
ENDIF
IF (PRESENT(distmin).AND.iopt_grid_sph.EQ.1) THEN
   IF (outunit.EQ.2) THEN
      distmin = radtodeg*distmin
   ELSEIF (outunit.EQ.3) THEN
      distmin = Rearth*distmin
   ENDIF
ENDIF


RETURN

END SUBROUTINE distance_minloc_d

!========================================================================

FUNCTION distance_pts_s(x1,x2,y1,y2,inunit,outunit) RESULT(dist)
!************************************************************************
!
! *distance_pts* Distance between two points
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.10.2
!
! Description - in case of Cartesian grid result is in meters
!             - in case of spherical grid result depends on value of outunit
!             - x1, x2, y1, y2 given in single precision    
!
! Module calls - complex_polar
!
!************************************************************************
!
USE physpars
USE switches
USE syspars
USE math_library, ONLY: complex_polar

!
!*Arguments
!
INTEGER, INTENT(IN) :: inunit, outunit
REAL (KIND=kndreal), INTENT(IN) :: x1, x2, y1, y2
REAL :: dist

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*x1*      REAL    X-coordinate (longitude) of first point
!*x2*      REAL    X-coordinate (longitude) of second point
!*y1*      REAL    Y-coordinate (latitude) of first point
!*y2*      REAL    Y-coordinate (latitude) of second point
!*inunit*  INTEGER unit of input locations (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!*outunit* INTEGER unit of result (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!                   = 3 => meters
!*dist*    REAL    Distance
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL :: dx, dy


IF (iopt_grid_sph.EQ.0) THEN
   dx = x2-x1; dy = y2-y1
   CALL complex_polar(dx,dy,xamp=dist)
ELSEIF (iopt_grid_sph.EQ.1) THEN
   IF (inunit.EQ.1) THEN
      dx = 0.5*(x2-x1); dy = 0.5*(y2-y1)
      dist = 2.0*ASIN(SQRT(COS(y1)*COS(y2)*SIN(dx)**2+SIN(dy)**2))
   ELSEIF (inunit.EQ.2) THEN
      dx = 0.5*degtorad*(x2-x1); dy = 0.5*degtorad*(y2-y1)
      dist = 2.0*ASIN(SQRT(COS(degtorad*y1)*COS(degtorad*y2)*SIN(dx)**2+&
                         & SIN(dy)**2))
   ENDIF
   IF (outunit.EQ.2) THEN
      dist = radtodeg*dist
   ELSEIF (outunit.EQ.3) THEN
      dist = Rearth*dist
   ENDIF
ENDIF


RETURN

END FUNCTION distance_pts_s

!========================================================================

FUNCTION distance_pts_d(x1,x2,y1,y2,inunit,outunit) RESULT(dist)
!************************************************************************
!
! *distance_pts* Distance between two points
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.10.2
!
! Description - in case of Cartesian grid result is in meters
!             - in case of spherical grid result depends on value of outunit
!             - x1, x2, y1, y2 given in double precision    
!
! Module calls - complex_polar
!
!************************************************************************
!
USE physpars
USE switches
USE syspars
USE math_library, ONLY: complex_polar

!
!*Arguments
!
INTEGER, INTENT(IN) :: inunit, outunit
REAL (KIND=kndrlong), INTENT(IN) :: x1, x2, y1, y2
REAL :: dist

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*x1*      REAL    X-coordinate (longitude) of first point
!*x2*      REAL    X-coordinate (longitude) of second point
!*y1*      REAL    Y-coordinate (latitude) of first point
!*y2*      REAL    Y-coordinate (latitude) of second point
!*inunit*  INTEGER unit of input locations (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!*outunit* INTEGER unit of result (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!                   = 3 => meters
!*dist*    REAL    Distance
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL :: dx, dy


IF (iopt_grid_sph.EQ.0) THEN
   dx = x2-x1; dy = y2-y1
   CALL complex_polar(dx,dy,xamp=dist)
ELSEIF (iopt_grid_sph.EQ.1) THEN
   IF (inunit.EQ.1) THEN
      dx = 0.5*(x2-x1); dy = 0.5*(y2-y1)
      dist = 2.0*ASIN(SQRT(COS(y1)*COS(y2)*SIN(dx)**2+SIN(dy)**2))
   ELSEIF (inunit.EQ.2) THEN
      dx = 0.5*degtorad*(x2-x1); dy = 0.5*degtorad*(y2-y1)
      dist = 2.0*ASIN(SQRT(COS(degtorad*y1)*COS(degtorad*y2)*SIN(dx)**2+&
                         & SIN(dy)**2))
   ENDIF
   IF (outunit.EQ.2) THEN
      dist = radtodeg*dist
   ELSEIF (outunit.EQ.3) THEN
      dist = Rearth*dist
   ENDIF
ENDIF


RETURN

END FUNCTION distance_pts_d

!========================================================================

FUNCTION find_obc_index_glb(i,j,cnode,ierr)
!************************************************************************
!
! *find_obc_index_glb* Find the global open boundary index of a global model
!                      grid point (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE grid
USE gridpars
USE syspars

!
!*Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: i, j
INTEGER, INTENT(OUT), OPTIONAL :: ierr
INTEGER :: find_obc_index_glb

!
! Name                Type    Purpose
!------------------------------------------------------------------------------
!*i*                  INTEGER X-index of model point
!*j*                  INTEGER Y-index of model point
!*cnode*              CHAR    Type of open boundary node ('U', 'V', 'X', 'Y')
!*ierr*               INTEGER Error code on exit
!*find_obc_index_glb* INTEGER Global open boundary index
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ii, indexob, jj


SELECT CASE (TRIM(cnode))
   CASE ('U')
      ii_110: DO ii=1,nobu
         IF ((iobu(ii).EQ.i).AND.(jobu(ii).EQ.j)) THEN
            indexob = ii
            EXIT ii_110
         ENDIF
      ENDDO ii_110
   CASE ('V')
      jj_120: DO jj=1,nobv
         IF ((iobv(jj).EQ.i).AND.(jobv(jj).EQ.j)) THEN
            indexob = jj
            EXIT jj_120
         ENDIF
      ENDDO jj_120
   CASE ('X')
      ii_130: DO ii=1,nobx
         IF ((iobx(ii).EQ.i).AND.(jobx(ii).EQ.j)) THEN
            indexob = ii
            EXIT ii_130
         ENDIF
      ENDDO ii_130
   CASE ('Y')
      jj_140: DO jj=1,noby
         IF ((ioby(jj).EQ.i).AND.(joby(jj).EQ.j)) THEN
            indexob = jj
            EXIT jj_140
         ENDIF
      ENDDO jj_140
END SELECT

IF (indexob.GT.0) THEN
   find_obc_index_glb = indexob
ELSE
   find_obc_index_glb = 0
   IF (PRESENT(ierr)) ierr = 1
ENDIF


RETURN

END FUNCTION find_obc_index_glb

!========================================================================

SUBROUTINE global_mask(maskglb,cnode,itype)
!************************************************************************
!
! *global_mask* Define a global land mask
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.10.1
!
! Description -
!
! Module calls - combine_mod
!
!************************************************************************
!
USE grid
USE gridpars
USE syspars
USE paral_comms, ONLY: combine_mod

!
!*Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: itype
LOGICAL, INTENT(OUT), DIMENSION(nc,nr) :: maskglb

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*maskglb* LOGICAL Output mask
!*cnode*   CHAR    Nodal type of output mask ('C','U','V','W')
!*itype*   INTEGER Selects type of mask
!              = 1 => .TRUE. except at land points
!              = 2 => .TRUE. at coastal boundaries and wet points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL, DIMENSION(ncloc,nrloc) :: maskloc


procname(pglev+1) = 'global_mask'
CALL log_timer_in()

!
!1. Local mask
!-------------
!

SELECT CASE (TRIM(cnode))
   CASE ('C','W')
      IF (itype.EQ.1) THEN
         maskloc = seapoint(1:ncloc,1:nrloc)
      ELSE
         maskloc = maskatc_int
      ENDIF
   CASE ('U')
      IF (itype.EQ.1) THEN
         maskloc = node2du(1:ncloc,1:nrloc).GT.0
      ELSEIF (itype.EQ.2) THEN
         maskloc = node2du(1:ncloc,1:nrloc).GT.1
      ENDIF
   CASE ('V')
      IF (itype.EQ.1) THEN
         maskloc = node2dv(1:ncloc,1:nrloc).GT.0
      ELSEIF (itype.EQ.2) THEN
         maskloc = node2dv(1:ncloc,1:nrloc).GT.1
      ENDIF
   CASE ('UV')
         maskloc = node2duv(1:ncloc,1:nrloc).GT.0
END SELECT

!
!2. Global mask
!--------------
!

IF (TRIM(cnode).EQ.'C'.OR.TRIM(cnode).EQ.'W'.AND.itype.EQ.1) THEN
   maskglb = seapointglb(1:nc,1:nr)
ELSE
   CALL combine_mod(maskglb,maskloc,(/1,1/),0,.FALSE.,mglevel=0,commall=.TRUE.)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE global_mask

!========================================================================

SUBROUTINE grid_rotation_params(iddesc)
!************************************************************************
!
! *grid_rotation_params* Determine parameters in coordinate transforms between
!                        reference and rotated grids   
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.5.1
!
! Description -
!
!************************************************************************
!
USE switches
USE syspars

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER  Grid file id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL :: gcos, gridangle, gsign, gsin, longpole, xcos, x0ref, x0rot, y0ref, y0rot


IF (iopt_grid_sph.EQ.0) RETURN

procname(pglev+1) = 'grid_rotation_params'
CALL log_timer_in()

!---local variables
x0ref = degtorad*surfacegrids(iddesc,1)%x0dat
y0ref = degtorad*surfacegrids(iddesc,1)%y0dat
y0rot = degtorad*surfacegrids(iddesc,1)%y0rot
gridangle =  degtorad*surfacegrids(iddesc,1)%gridangle
gcos = COS(gridangle); gsin = SIN(gridangle)
gsign = SIGN(1.0,gcos)

!---longitude of rotated pole
xcos = (SIN(y0rot)-gcos*SIN(y0ref))/(gsin*COS(y0ref))
xcos = MIN(xcos,1.0); xcos = MAX(xcos,-1.0)
longpole = x0ref-gsign*ABS(ACOS(xcos))
surfacegrids(iddesc,1)%longpole = radtodeg*longpole

!---longitude of SW corner in rotated grid
xcos = (gsin*SIN(y0ref)-gcos*COS(x0ref-longpole)*COS(y0ref))/COS(y0rot)
xcos = MIN(xcos,1.0); xcos = MAX(xcos,-1.0)
x0rot = SIGN(1.0,SIN(longpole-x0ref))*ACOS(xcos)
IF (x0rot.GT.pi) x0rot = x0rot - twopi
IF (x0rot.LT.-pi) x0rot = x0rot + twopi
surfacegrids(iddesc,1)%x0rot = radtodeg*x0rot

CALL log_timer_out()


RETURN

END SUBROUTINE grid_rotation_params

!========================================================================

SUBROUTINE grid_rotation_transf(xcoord,ycoord,idir,inunit,iddesc)
!************************************************************************
!
! *grid_rotation_transf* Transform geographical coordinates to the ones
!                        in a rotated grid or vice versa
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.5.1
!
! Description -
!
!************************************************************************
!
USE switches
USE syspars
USE error_routines, ONLY: error_alloc

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, idir, inunit
REAL, INTENT(INOUT), DIMENSION(:,:) :: xcoord, ycoord

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*xcoord*  REAL     X-coordinates on old coordinate system on input,
!                   on new coordinate system on output
!*ycoord*  REAL     Y-coordinates on old coordinate system on input,
!                   on new coordinate system on output
!*idir*    INTEGER  Type of transformation
!               = 1  => from reference to rotated coordinates
!               =-1  => from rotated to reference coordinates
!*inunit*  INTEGER unit of input locations (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!*iddesc*  INTEGER  Grid file id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: nx, ny
REAL :: alpha, beta, betax, gcos, gsin, x0ref, y0ref
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d, x, y, ycos, ysin


procname(pglev+1) = 'grid_rotation_transf'
CALL log_timer_in()

nx = SIZE(xcoord,DIM=1); ny = SIZE(xcoord,DIM=2)
alpha = degtorad*surfacegrids(iddesc,1)%gridangle

!
!1. Allocate
!-----------
!

IF (iopt_grid_sph.EQ.0) THEN
   ALLOCATE (x(nx,ny),STAT=errstat)
   CALL error_alloc('x',2,(/nx,ny/),kndrtype)
   ALLOCATE (y(nx,ny),STAT=errstat)
   CALL error_alloc('y',2,(/nx,ny/),kndrtype)
ELSE
   ALLOCATE (array2d(nx,ny),STAT=errstat)
   CALL error_alloc('array2d',2,(/nx,ny/),kndrtype)
   ALLOCATE (ycos(nx,ny),STAT=errstat)
   CALL error_alloc('ycos',2,(/nx,ny/),kndrtype)
   ALLOCATE (ysin(nx,ny),STAT=errstat)
   CALL error_alloc('ysin',2,(/nx,ny/),kndrtype)
ENDIF

!
!2. Cartesian case
!-----------------
!

IF (iopt_grid_sph.EQ.0) THEN

   gcos = COS(alpha); gsin = SIN(alpha)
   x0ref = surfacegrids(iddesc,1)%x0dat
   y0ref = surfacegrids(iddesc,1)%y0dat

!  ---transform to new coordinates
   IF (idir.EQ.1) THEN
      x = xcoord - x0ref; y = ycoord - y0ref
      xcoord = gcos*x + gsin*y
      ycoord = gcos*y - gsin*x

!  ---transform to old coordinates
   ELSE
      x = xcoord; y = ycoord
      xcoord = x0ref + gcos*x - gsin*y
      ycoord = y0ref + gcos*y + gsin*x
   END IF

!
!3. Spherical case
!-----------------
!

ELSE

   beta = degtorad*surfacegrids(iddesc,1)%longpole
   gcos = COS(alpha); gsin = SIN(alpha)
   betax = beta - SIGN(1.0,beta)*pi
   IF (inunit.EQ.2) THEN
      xcoord = degtorad*xcoord; ycoord = degtorad*ycoord
   ENDIF

!  ---transform to new coordinates
   IF (idir.EQ.1) THEN
      array2d = COS(xcoord-beta)
      ycos = COS(ycoord); ysin = SIN(ycoord)
      ycoord = ASIN(gcos*SIN(ycoord)+gsin*COS(ycoord)*array2d)
      array2d = (gsin*ysin-gcos*ycos*array2d)/COS(ycoord)
      array2d = MIN(1.0,array2d); array2d = MAX(-1.0,array2d)
      array2d = ACOS(array2d)
      xcoord = SIGN(1.0,SIN(beta-xcoord))*array2d
      WHERE (xcoord.GT.pi) xcoord = xcoord - twopi
      WHERE (xcoord.LT.-pi) xcoord = xcoord + twopi

!  ---transform to old coordinates
   ELSE
      ycos = COS(ycoord); ysin = SIN(ycoord)
      ycoord = ASIN(gcos*SIN(ycoord)+gsin*ycos*COS(xcoord))
      array2d = (gcos*ycos*COS(xcoord)-gsin*ysin)/COS(ycoord)
      array2d = MIN(1.0,array2d); array2d = MAX(-1.0,array2d)
      xcoord = SIGN(1.0,SIN(xcoord))*ACOS(array2d) + betax
      WHERE (xcoord.GT.pi) xcoord = xcoord - twopi
      WHERE (xcoord.LT.-pi) xcoord = xcoord + twopi
   ENDIF

   IF (inunit.EQ.2) THEN
      xcoord = radtodeg*xcoord; ycoord = radtodeg*ycoord
   ENDIF

ENDIF

!
!4. Deallocate
!-------------
!

IF (iopt_grid_sph.EQ.0) THEN
   DEALLOCATE (x,y)
ELSE
   DEALLOCATE (array2d,ycos,ysin)
ENDIF

CALL log_timer_out()    


RETURN

END SUBROUTINE grid_rotation_transf

!========================================================================

SUBROUTINE grid_spacings_curv(dist,x1,x2,y1,y2,inunit,outunit)
!************************************************************************
!
! *grid_spacings_curv*  Horizontal grid spacings for a model/data curvilinear
!                       grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.10.2
!
! Description - in case of Cartesian grid result is in meters
!             - in case of spherical grid result depends on value of outunit
!
! Module calls - complex_polar
!
!************************************************************************
!
USE physpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE math_library, ONLY: complex_polar

!
!*Arguments
!
INTEGER, INTENT(IN) :: inunit, outunit
REAL, INTENT(IN), DIMENSION(:,:)  :: x1, x2, y1, y2
REAL, INTENT(OUT), DIMENSION(:,:) :: dist

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*dist*    REAL    Array of grid spacings
!*x1*      REAL    X-coordinate (longitude) of first point
!*x2*      REAL    X-coordinate (longitude) of second point
!*y1*      REAL    Y-coordinate (latitude) of first point
!*y2*      REAL    Y-coordinate (latitude) of second point
!*inunit*  INTEGER Unit of input locations (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!*outunit* INTEGER Unit of result (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!                   = 3 => meters
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: nx, ny
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: dx, dy


procname(pglev+1) = 'grid_spacings_curv'
CALL log_timer_in()

!
!1. Allocate
!-----------
!

nx = SIZE(dist,DIM=1); ny = SIZE(dist,DIM=2)
ALLOCATE (dx(nx,ny),STAT=errstat)
CALL error_alloc('dx',2,(/nx,ny/),kndrtype)
ALLOCATE (dy(nx,ny),STAT=errstat)
CALL error_alloc('dy',2,(/nx,ny/),kndrtype)

!
!2. Cartesian grid
!-----------------
!

IF (iopt_grid_sph.EQ.0) THEN

   dx = x2-x1; dy = y2-y1
   CALL complex_polar(dx,dy,xamp=dist)

!
!3. Spherical grid
!-----------------
!
   
ELSEIF (iopt_grid_sph.EQ.1) THEN

   IF (inunit.EQ.1) THEN
      dx = 0.5*(x2-x1); dy = 0.5*(y2-y1)
      dist = 2.0*ASIN(SQRT(COS(y1)*COS(y2)*SIN(dx)**2+SIN(dy)**2))
   ELSEIF (inunit.EQ.2) THEN
      dx = 0.5*degtorad*(x2-x1); dy = 0.5*degtorad*(y2-y1)
      dist = 2.0*ASIN(SQRT(COS(degtorad*y1)*COS(degtorad*y2)*SIN(dx)**2+&
                         & SIN(dy)**2))
   ENDIF
   IF (outunit.EQ.2) dist = radtodeg*dist
   IF (outunit.EQ.3) dist = Rearth*dist

ENDIF

!
!4. Deallocate
!-------------
!

DEALLOCATE (dx,dy)

CALL log_timer_out()


RETURN

END SUBROUTINE grid_spacings_curv

!========================================================================

SUBROUTINE grid_spacings_rectang(dist,hdel,ycoord,inunit,outunit,cdir)
!************************************************************************
!
! *grid_spacings_rectang* Horizontal grid spacings for a model/data
!                         rectangular grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.10.2
!
! Description - in case of Cartesian grid result is in meters
!             - in case of spherical grid result depends on value of outunit
!
! Module calls -
!
!************************************************************************
!
USE physpars
USE switches

!
!*Arguments
!
CHARACTER (LEN=1), INTENT(IN) :: cdir
INTEGER, INTENT(IN) :: inunit, outunit
REAL, INTENT(IN), DIMENSION(:) :: hdel, ycoord
REAL, INTENT(OUT), DIMENSION(:,:) :: dist

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*dist*    REAL    Array of grid spacings
!*hdel*    REAL    Coordinate spacings
!*ycoord*  REAL    Y-coordinates (latitude)
!*inunit*  INTEGER Unit of input locations (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!*outunit* INTEGER Unit of result (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!                   = 3 => meters
!*cdir*    CHAR    Coordinate direction
!                   = 'X' => X-direction
!                   = 'Y' => Y-direction 
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j, nx, ny


procname(pglev+1) = 'grid_spacings_rectang'
CALL log_timer_in()

nx = SIZE(dist,DIM=1); ny = SIZE(dist,DIM=2)

!
!1. X-direction
!--------------
!

IF (cdir.EQ.'X') THEN

   IF (iopt_grid_sph.EQ.0) THEN
      j_111: DO j=1,ny
         dist(:,j) = hdel
      ENDDO j_111
   ELSE
      IF (inunit.EQ.1) THEN
         j_112: DO j=1,ny
            dist(:,j) = 2.0*ASIN(COS(ycoord(j))*SIN(0.5*hdel))
         ENDDO j_112
      ELSEIF (inunit.EQ.2) THEN
         j_113: DO j=1,ny
            dist(:,j) = 2.0*ASIN(COS(degtorad*ycoord(j))*SIN(0.5*degtorad*hdel))
         ENDDO j_113
      ENDIF
      IF (outunit.EQ.2) dist = radtodeg*dist
      IF (outunit.EQ.3) dist = Rearth*dist
   ENDIF

!
!2. Y-direction
!--------------
!

ELSE

   i_210: DO i=1,nx
      dist(i,:) = hdel
   ENDDO i_210
   IF (iopt_grid_sph.EQ.1) THEN
      IF (inunit.EQ.1) THEN
         IF (outunit.EQ.2) THEN
            dist = radtodeg*dist
         ELSEIF (outunit.EQ.3) THEN
            dist = Rearth*dist
         ENDIF
      ELSE
         IF (outunit.EQ.1) THEN
            dist = degtorad*dist
         ELSEIF (outunit.EQ.3) THEN
            dist = (Rearth*degtorad)*dist
         ENDIF
      ENDIF
   ENDIF
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE grid_spacings_rectang

!========================================================================

FUNCTION local_proc(i,j,iproc,mglevel)
!************************************************************************
!
! *local_proc* Returns .TRUE. if the grid point with global indices (i,j)
!              belongs to the local domain
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.8
!
! Description - If iproc is present, returns .TRUE. if (i,j) belongs
!               to process iproc
!
!************************************************************************
!
USE multigrid
USE grid
USE gridpars

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j
INTEGER, INTENT(IN), OPTIONAL :: iproc, mglevel
LOGICAL :: local_proc

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*i*       INTEGER Global or local X-index of grid point
!*j*       INTEGER Blobal or local Y-index of grid point
!*iproc*   INTEGER Process number
!*mglevel* INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
! 
!*Local variables
!
INTEGER :: lev


IF (PRESENT(mglevel)) THEN
   lev = mglevel
ELSE
   lev = 0
ENDIF

IF (lev.EQ.0) THEN
   IF (PRESENT(iproc)) THEN
      IF (i.GE.nc1procs(iproc).AND.i.LE.nc2procs(iproc).AND.&
        & j.GE.nr1procs(iproc).AND.j.LE.nr2procs(iproc)) THEN
         local_proc = .TRUE.
      ELSE
         local_proc = .FALSE.
      ENDIF
   ELSE
      IF ((i.GE.1.AND.i.LE.ncloc).AND.(j.GE.1.AND.j.LE.nrloc)) THEN
         local_proc = .TRUE.
      ELSE
         local_proc = .FALSE.
      ENDIF
   ENDIF
ELSE
   IF (PRESENT(iproc)) THEN
      IF (i.GE.mgvars(lev)%nc1procs(iproc).AND.&
        & i.LE.mgvars(lev)%nc2procs(iproc).AND.&
        & j.GE.mgvars(lev)%nr1procs(iproc).AND.&
        & j.LE.mgvars(lev)%nr2procs(iproc)) THEN
         local_proc = .TRUE.
      ELSE
         local_proc = .FALSE.
      ENDIF
   ELSE
      IF ((i.GE.1.AND.i.LE.mgvars(lev)%ncloc).AND.&
        & (j.GE.1.AND.j.LE.mgvars(lev)%nrloc)) THEN
         local_proc = .TRUE.
      ELSE
         local_proc = .FALSE.
      ENDIF
   ENDIF
ENDIF

RETURN

END FUNCTION local_proc

!========================================================================

FUNCTION mask_array(idim,jdim,lims,node)
!************************************************************************
!
! *mask_array* Mask array for user-defined output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE grid
USE gridpars

!
!*Arguments
!
CHARACTER (LEN=1), INTENT(IN) :: node
INTEGER, INTENT(IN) :: idim, jdim
INTEGER, INTENT(IN), DIMENSION(3,2) :: lims
LOGICAL, DIMENSION(idim,jdim) :: mask_array

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*idim*    INTEGER  X-dimension of result array
!*jdim*    INTEGER  Y-dimension of result array
!*lims*    INTEGER  Start/end/increment values of space indices
!                   in X- and Y-direction
!*node*    CHAR     type of node where result array is evaluated 
!                   ('C','U','V','W')
!
!------------------------------------------------------------------------------
!

SELECT CASE(node)
CASE ('C','W')
   mask_array = nodeatc(lims(1,1):lims(2,1):lims(3,1),&
                      & lims(1,2):lims(2,2):lims(3,2)).EQ.1
CASE ('U')
   mask_array = node2du(lims(1,1):lims(2,1):lims(3,1),&
                      & lims(1,2):lims(2,2):lims(3,2)).GT.1
CASE ('V')
   mask_array = node2dv(lims(1,1):lims(2,1):lims(3,1),&
                      & lims(1,2):lims(2,2):lims(3,2)).GT.1
END SELECT


RETURN

END FUNCTION mask_array

!========================================================================

FUNCTION num_proc(i,j)
!************************************************************************
!
! *num_proc* Returns the process number containing the point with global
!            indices (i,j)
!
! Author - Patrick Luyten and Stephanie Ponsar
!
! Version - @(COHERENS)grid_routines.f90  V2.0
!
! Description -
!
!************************************************************************
!
USE grid
USE paralpars

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j
INTEGER :: num_proc

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*i*        INTEGER Global X-index of grid point
!*j*        INTEGER Global Y-index of grid point
!*num_proc* INTEGER Process number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: iproc


num_proc = -1

iproc_110: DO iproc=1,nprocs
   IF ((i.GE.nc1procs(iproc).AND.i.LE.nc2procs(iproc)).AND.&
     & (j.GE.nr1procs(iproc).AND.j.LE.nr2procs(iproc))) THEN
      num_proc = iproc
      EXIT iproc_110
   ENDIF
ENDDO iproc_110


RETURN

END FUNCTION num_proc

!========================================================================

SUBROUTINE rotate_vec_0d(vxin,vyin,vxout,vyout,idir,angle)
!************************************************************************
!
! *rotate_vec_0d* Rotate a scalar vector over an angle equal to 'angle'
!                 if idir = 1 or '-angle' if idir = -1 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.5.1
!
! Description -
!
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: idir
REAL, INTENT(IN) :: angle
REAL, INTENT(IN) :: vxin, vyin
REAL, INTENT(OUT) :: vxout, vyout

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*vxin*    REAL     X-component of input vector
!*vyin*    REAL     Y-component of input vector
!*vxout*   REAL     X-component of rotated output vector
!*vyout*   REAL     Y-component of rotated output vector
!*idir*    INTEGER  Direction of the rotation operation
!                    =  1 => from "reference" to "rotated" grid 
!                    = -1 => from "rotated" to "reference" grid
!*angle*   REAL     Angle of rotation                                      [rad]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL :: gcos, gsin


gcos = COS(angle); gsin = idir*SIN(angle)
vxout = gcos*vxin + gsin*vyin
vyout = gcos*vyin - gsin*vxin


RETURN

END SUBROUTINE rotate_vec_0d

!========================================================================

SUBROUTINE rotate_vec_1d(vxin,vyin,vxout,vyout,idir,angle)
!************************************************************************
!
! *rotate_vec_1d* Rotate a vector profile over an angle equal to 'angle'
!                 if idir = 1 or '-angle' if idir = -1 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.5.1
!
! Description -
!
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: idir
REAL, INTENT(IN) :: angle
REAL, INTENT(IN), DIMENSION(:) :: vxin, vyin
REAL, INTENT(OUT), DIMENSION(:) :: vxout, vyout

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*vxin*    REAL     X-component of input vector
!*vyin*    REAL     Y-component of input vector
!*vxout*   REAL     X-component of rotated output vector
!*vyout*   REAL     Y-component of rotated output vector
!*idir*    INTEGER  Direction of the rotation operation
!                    =  1 => from "reference" to "rotated" grid 
!                    = -1 => from "rotated" to "reference" grid
!*angle*   REAL     Angle of rotation                                      [rad]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL :: gcos, gsin


gcos = COS(angle); gsin = idir*SIN(angle)
vxout = gcos*vxin + gsin*vyin
vyout = gcos*vyin - gsin*vxin


RETURN

END SUBROUTINE rotate_vec_1d

!========================================================================

SUBROUTINE rotate_vec_2d(vxin,vyin,vxout,vyout,angle,idir)
!************************************************************************
!
! *rotate_vec_2d* Rotate a 2-D vector over an angle given by 'angle' if
!                 idir = 1 or '-angle' if idir = -1 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.5.1
!
! Description -
!
!************************************************************************
!
USE syspars
USE error_routines, ONLY: error_alloc

!
!*Arguments
!
INTEGER, INTENT(IN) :: idir
REAL, INTENT(IN), DIMENSION(:,:) :: angle
REAL, INTENT(IN), DIMENSION(:,:) :: vxin, vyin
REAL, INTENT(OUT), DIMENSION(:,:) :: vxout, vyout

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*vxin*    REAL     X-component of input vector
!*vyin*    REAL     Y-component of input vector
!*vxout*   REAL     X-component of rotated output vector
!*vyout*   REAL     Y-component of rotated output vector
!*angle*   REAL     Angle of rotation                                      [rad]
!*idir*    INTEGER  Direction of the rotation operation
!                    =  1 => from "reference" to "rotated" grid 
!                    = -1 => from "rotated" to "reference" grid
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(2) :: ndims
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: gcos, gsin


procname(pglev+1) = 'rotate_vec_2d'
CALL log_timer_in()

ndims = SHAPE(angle)
ALLOCATE (gcos(ndims(1),ndims(2)),STAT=errstat)
CALL error_alloc('gcos',2,(/ndims(1),ndims(2)/),kndrtype)
ALLOCATE (gsin(ndims(1),ndims(2)),STAT=errstat)
CALL error_alloc('gsin',2,(/ndims(1),ndims(2)/),kndrtype)

gcos = COS(angle); gsin = idir*SIN(angle)
vxout = gcos*vxin + gsin*vyin
vyout = gcos*vyin - gsin*vxin

DEALLOCATE (gcos,gsin)

CALL log_timer_out()


RETURN

END SUBROUTINE rotate_vec_2d

!========================================================================

SUBROUTINE rotate_vec_3d(vxin,vyin,vxout,vyout,angle,idir)
!************************************************************************
!
! *rotate_vec_3d* Rotate a 3-D vector over an angle given by 'angle' if
!                 idir = 1 or '-angle' if idir = -1 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.5.1
!
! Description -
!
!************************************************************************
!
USE syspars
USE error_routines, ONLY: error_alloc

!
!*Arguments
!
INTEGER, INTENT(IN) :: idir
REAL, INTENT(IN), DIMENSION(:,:) :: angle
REAL, INTENT(IN), DIMENSION(:,:,:) :: vxin, vyin
REAL, INTENT(OUT), DIMENSION(:,:,:) :: vxout, vyout

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*vxin*    REAL     X-component of input vector
!*vyin*    REAL     Y-component of input vector
!*vxout*   REAL     X-component of rotated output vector
!*vyout*   REAL     Y-component of rotated output vector
!*angle*   REAL     Angle of rotation                                      [rad]
!*idir*    INTEGER  Direction of the rotation operation
!                    =  1 => from "reference" to "rotated" grid 
!                    = -1 => from "rotated" to "reference" grid

!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k
INTEGER, DIMENSION(3) :: ndims
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: gcos, gsin


procname(pglev+1) = 'rotate_vec_3d'
CALL log_timer_in()

ndims = SHAPE(vxin)
ALLOCATE (gcos(ndims(1),ndims(2)),STAT=errstat)
CALL error_alloc('gcos',2,(/ndims(1),ndims(2)/),kndrtype)
ALLOCATE (gsin(ndims(1),ndims(2)),STAT=errstat)
CALL error_alloc('gsin',2,(/ndims(1),ndims(2)/),kndrtype)

gcos = COS(angle); gsin = idir*SIN(angle)

k_110: DO k=1,ndims(3)
   vxout(:,:,k) = gcos*vxin(:,:,k) + gsin*vyin(:,:,k)
   vyout(:,:,k) = gcos*vyin(:,:,k) - gsin*vxin(:,:,k)
ENDDO k_110

DEALLOCATE (gcos,gsin)

CALL log_timer_out()


RETURN

END SUBROUTINE rotate_vec_3d

!========================================================================

SUBROUTINE rotate_vec_4d(vxin,vyin,vxout,vyout,angle,idir)
!************************************************************************
!
! *rotate_vec_4d* Rotate a 4-D vector over an angle given by 'angle' if
!                 idir = 1 or '-angle' if idir = -1 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.5.1
!
! Description -
!
!************************************************************************
!
USE syspars
USE error_routines, ONLY: error_alloc

!
!*Arguments
!
INTEGER, INTENT(IN) :: idir
REAL, INTENT(IN), DIMENSION(:,:) :: angle
REAL, INTENT(IN), DIMENSION(:,:,:,:) :: vxin, vyin
REAL, INTENT(OUT), DIMENSION(:,:,:,:) :: vxout, vyout

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*vxin*    REAL     X-component of input vector
!*vyin*    REAL     Y-component of input vector
!*vxout*   REAL     X-component of rotated output vector
!*vyout*   REAL     Y-component of rotated output vector
!*angle*   REAL     Angle of rotation                                      [rad]
!*idir*    INTEGER  Direction of the rotation operation
!                    =  1 => from "reference" to "rotated" grid 
!                    = -1 => from "rotated" to "reference" grid

!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, l
INTEGER, DIMENSION(4) :: ndims
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: gcos, gsin


procname(pglev+1) = 'rotate_vec_4d'
CALL log_timer_in()


ndims = SHAPE(vxin)
ALLOCATE (gcos(ndims(1),ndims(2)),STAT=errstat)
CALL error_alloc('gcos',2,(/ndims(1),ndims(2)/),kndrtype)
ALLOCATE (gsin(ndims(1),ndims(2)),STAT=errstat)
CALL error_alloc('gsin',2,(/ndims(1),ndims(2)/),kndrtype)

gcos = COS(angle); gsin = idir*SIN(angle)

l_110: DO l=1,ndims(4)
k_110: DO k=1,ndims(3)
   vxout(:,:,k,l) = gcos*vxin(:,:,k,l) + gsin*vyin(:,:,k,l)
   vyout(:,:,k,l) = gcos*vyin(:,:,k,l) - gsin*vxin(:,:,k,l)
ENDDO k_110
ENDDO l_110

DEALLOCATE (gcos,gsin)

CALL log_timer_out()


RETURN

END SUBROUTINE rotate_vec_4d

!========================================================================

SUBROUTINE Zcoord_arr(zcoord,lbounds,ubounds,cnode,meanlevel,outflag)
!************************************************************************
!
! *Zcoord_arr* Array of Z-coordinates
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.2
!
! Description - returns zero at dry cells or nodes
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE syspars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: meanlevel
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                           & lbounds(3):ubounds(3)) :: zcoord

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*zcoord*    REAL    Array of Z-coordinates
!*lbounds*   INTEGER Lower bounds of array
!*ubounds*   INTEGER Upper bounds of array
!*cnode*     CHAR    Type of node ('C', 'U', 'V', 'W', 'UW', 'VW')
!*meanlevel* LOGICAL Uses mean/total water depths if .TRUE./.FALSE.
!*outflag*   REAL    Output flag for dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, l1, l2, l3, u1, u2, u3
REAL :: flagval


procname(pglev+1) = 'Zcoord_arr'
CALL log_timer_in()

!
!1. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!2. Optional argument
!--------------------
!

IF (PRESENT(outflag)) THEN
   flagval = outflag
ELSE
   flagval = 0.0
ENDIF

!
!2. Z-coordinates using mean water depths
!----------------------------------------
!

IF (meanlevel) THEN

   SELECT CASE (TRIM(cnode))

!  ---C-nodes
   CASE ('C')
      k_210: DO k=l3,u3
         WHERE (nodeatc(l1:u1,l2:u2).GT.0)
            zcoord(:,:,k) = depmeanatc(l1:u1,l2:u2)*&
                          & (gscoordatc(l1:u1,l2:u2,k)-1.0)
         ELSEWHERE
            zcoord(:,:,k) = flagval
         END WHERE
      ENDDO k_210

!  ---U-nodes
   CASE ('U')
      k_220: DO k=l3,u3
         WHERE (node2du(l1:u1,l2:u2).GT.1)
            zcoord(:,:,k) = depmeanatu(l1:u1,l2:u2)*&
                          & (gscoordatu(l1:u1,l2:u2,k)-1.0)
         ELSEWHERE
            zcoord(:,:,k) = flagval
         END WHERE
      ENDDO k_220

!  ---V-nodes
   CASE ('V')
      k_230: DO k=l3,u3
         WHERE (node2dv(l1:u1,l2:u2).GT.1)
            zcoord(:,:,k) = depmeanatv(l1:u1,l2:u2)*&
                          & (gscoordatv(l1:u1,l2:u2,k)-1.0)
         ELSEWHERE
            zcoord(:,:,k) = flagval
         END WHERE
      ENDDO k_230

!  ---W-nodes
   CASE ('W')
      k_240: DO k=l3,u3
         WHERE (nodeatc(l1:u1,l2:u2).GT.0)
            zcoord(:,:,k) = depmeanatc(l1:u1,l2:u2)*&
                          & (gscoordatw(l1:u1,l2:u2,k)-1.0)
         ELSEWHERE
            zcoord(:,:,k) = flagval
         END WHERE
      ENDDO k_240
      
!  ---UW-nodes
   CASE ('UW')
      k_250: DO k=l3,u3
         WHERE (node2du(l1:u1,l2:u2).GT.1)
            zcoord(:,:,k) = depmeanatu(l1:u1,l2:u2)*&
                          & (gscoordatuw(l1:u1,l2:u2,k)-1.0)
         ELSEWHERE
            zcoord(:,:,k) = flagval
         END WHERE
      ENDDO k_250

!  ---V-nodes
   CASE ('VW')
      k_260: DO k=l3,u3
         WHERE (node2dv(l1:u1,l2:u2).GT.1)
            zcoord(:,:,k) = depmeanatv(l1:u1,l2:u2)*&
                          & (gscoordatvw(l1:u1,l2:u2,k)-1.0)
         ELSEWHERE
            zcoord(:,:,k) = flagval
         END WHERE
      ENDDO k_260

   END SELECT

!
!3. Z-coordinates using total water depths
!-----------------------------------------
!

ELSE
   SELECT CASE (TRIM(cnode))

!  ---C-nodes
   CASE ('C')
      k_310: DO k=l3,u3
         WHERE (nodeatc(l1:u1,l2:u2).GT.0) 
            zcoord(:,:,k) = gscoordatc(l1:u1,l2:u2,k)*&
                          & deptotatc(l1:u1,l2:u2) - depmeanatc(l1:u1,l2:u2)
         ELSEWHERE
            zcoord(:,:,k) = flagval
         END WHERE
      ENDDO k_310

!  ---U-nodes
   CASE ('U')
      k_320: DO k=l3,u3
         WHERE (node2du(l1:u1,l2:u2).GT.1) 
            zcoord(:,:,k) = gscoordatu(l1:u1,l2:u2,k)*&
                          & deptotatu(l1:u1,l2:u2) - depmeanatu(l1:u1,l2:u2)
         ELSEWHERE
            zcoord(:,:,k) = flagval
         END WHERE
      ENDDO k_320

!  ---V-nodes
   CASE ('V')
      k_330: DO k=l3,u3
         WHERE (node2dv(l1:u1,l2:u2).GT.1) 
            zcoord(:,:,k) = gscoordatv(l1:u1,l2:u2,k)*&
                          & deptotatv(l1:u1,l2:u2) - depmeanatv(l1:u1,l2:u2)
         ELSEWHERE
            zcoord(:,:,k) = flagval
         END WHERE
      ENDDO k_330

!  ---W-nodes
   CASE ('W')
      k_340: DO k=l3,u3
         WHERE (nodeatc(l1:u1,l2:u2).GT.0) 
            zcoord(:,:,k) = gscoordatw(l1:u1,l2:u2,k)*&
                          & deptotatc(l1:u1,l2:u2) - depmeanatc(l1:u1,l2:u2)
         ELSEWHERE
            zcoord(:,:,k) = flagval
         END WHERE
      ENDDO k_340

!  ---UW-nodes
   CASE ('UW')
      k_350: DO k=l3,u3
         WHERE (node2du(l1:u1,l2:u2).GT.1) 
            zcoord(:,:,k) = gscoordatuw(l1:u1,l2:u2,k)*&
                          & deptotatu(l1:u1,l2:u2) - depmeanatu(l1:u1,l2:u2)
         ELSEWHERE
            zcoord(:,:,k) = flagval
         END WHERE
      ENDDO k_350

!  ---V-nodes
   CASE ('VW')
      k_360: DO k=l3,u3
         WHERE (node2dv(l1:u1,l2:u2).GT.1) 
            zcoord(:,:,k) = gscoordatvw(l1:u1,l2:u2,k)*&
                          & deptotatv(l1:u1,l2:u2) - depmeanatv(l1:u1,l2:u2)
         ELSEWHERE
            zcoord(:,:,k) = flagval
         END WHERE
      ENDDO k_360

   END SELECT

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE Zcoord_arr

!========================================================================

FUNCTION Zcoord_var(i,j,k,cnode,meanlevel)
!************************************************************************
!
! *Zcoord_var* Z-coordinate at grid point (i,j,k)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_routines.f90  V2.0
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE syspars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: meanlevel
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: i, j, k
REAL :: Zcoord_var

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*i*          INTEGER X-index of gri point
!*j*          INTEGER Y-index of gri point
!*k*          INTEGER Vertical index of grid point
!*cnode*      CHAR    Type of node ('C', 'U', 'V', 'W', UW, 'VW')
!*meanlevel*  LOGICAL Uses mean/total water depths if .TRUE./.FALSE.
!*Zcoord_var* REAL    Z-coordinate 
!
!------------------------------------------------------------------------------
!
!1. Z-coordinate using mean water depth
!--------------------------------------
!

IF (meanlevel) THEN
   SELECT CASE (TRIM(cnode))
      CASE ('C')
         Zcoord_var = depmeanatc(i,j)*(gscoordatc(i,j,k)-1.0)
      CASE ('U')
         Zcoord_var = depmeanatu(i,j)*(gscoordatu(i,j,k)-1.0)
      CASE ('V')
         Zcoord_var = depmeanatv(i,j)*(gscoordatv(i,j,k)-1.0)
      CASE ('W')
         Zcoord_var = depmeanatc(i,j)*(gscoordatw(i,j,k)-1.0)
      CASE ('UW')
         Zcoord_var = depmeanatu(i,j)*(gscoordatuw(i,j,k)-1.0)
      CASE ('VW')
         Zcoord_var = depmeanatv(i,j)*(gscoordatvw(i,j,k)-1.0)
      END SELECT

!
!2. Z-coordinate using total water depth
!---------------------------------------
!
 
ELSE
   SELECT CASE (TRIM(cnode))
      CASE ('C')
         Zcoord_var = gscoordatc(i,j,k)*deptotatc(i,j) - depmeanatc(i,j)
      CASE ('U')
         Zcoord_var = gscoordatu(i,j,k)*deptotatu(i,j) - depmeanatu(i,j)
      CASE ('V')
         Zcoord_var = gscoordatv(i,j,k)*deptotatv(i,j) - depmeanatv(i,j)
      CASE ('W')
         Zcoord_var = gscoordatw(i,j,k)*deptotatc(i,j) - depmeanatc(i,j)
      CASE ('UW')
         Zcoord_var = gscoordatuw(i,j,k)*deptotatu(i,j) - depmeanatu(i,j)
      CASE ('VW')
         Zcoord_var = gscoordatvw(i,j,k)*deptotatv(i,j) - depmeanatv(i,j)
   END SELECT
ENDIF


RETURN

END FUNCTION Zcoord_var


END MODULE grid_routines
