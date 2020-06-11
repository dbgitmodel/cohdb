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

!************************************************************************
!
! *Usrdef_Surface_Data* User-defined surface data setup (example routines)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.0
!
! $Date: 2013-04-16 12:20:30 +0200 (Tue, 16 Apr 2013) $
!
! $Revision: 548 $
!
! Description - 
!
! Reference -
!
! Routines - usrdef_surface_absgrid, usrdef_surface_relgrd,
!            usrdef_surface_data
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_surface_absgrd(iddesc,ifil,n1dat,n2dat,xcoord,ycoord)
!************************************************************************
!
! *usrdef_surface_absgrd* Define coordinate arrays of surface grid(s)
!                         with respect to model grid (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_surface_input_grid, define_surface_output_grid
!
! External calls -
!
! Module calls -
!
!************************************************************************
!
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, n1dat, n2dat
REAL, INTENT(INOUT), DIMENSION(n1dat,n2dat) :: xcoord, ycoord

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Grid file id
!*ifil*      INTEGER No. of grid file
!*n1dat*     INTEGER X-dimension of data grid
!*n2dat*     INTEGER Y-dimension of data grid
!*xcoord*    REAL    X-coordinates of data grid
!*ycoord*    REAL    Y-coordinates of data grid
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'usrdef_surface_absgrd'
CALL log_timer_in()

xcoord = ?
ycoord = ?

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_surface_absgrd

!========================================================================

SUBROUTINE usrdef_surface_relgrd(iddesc,ifil,surfgridglb,nx,ny,nonodes)
!************************************************************************
!
! *usrdef_surface_relgrd* Define relative coordinate arrays of surface grid(s)
!                         with respect to model grid (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_surface_input_grid, define_surface_output_grid
!
! External calls -
!
! Module calls -
!
!************************************************************************
!
USE datatypes
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE
 
!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, nonodes, nx, ny
TYPE (HRelativeCoords), INTENT(INOUT), DIMENSION(nx,ny,nonodes) :: surfgridglb

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER Grid file id
!*ifil*        INTEGER No. of grid file
!*surfgridglb* DERIVED Relative coordinates of model grid to data grid or data
!                      grid to model grid
!*nx*          INTEGER X-dimension of data grid
!*ny*          INTEGER Y-dimension of data grid 
!*nonodes*     INTEGER Number of nodes
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: inode


procname(pglev+1) = 'usrdef_surface_relgrd'
CALL log_timer_in()

inode_110: DO inode=1,nonodes
   surfgridglb(1:nx,1:ny,inode)%icoord = ?
   surfgridglb(1:nx,1:ny,inode)%jcoord = ?
   surfgridglb(1:nx,1:ny,inode)%xcoord = ?
   surfgridglb(1:nx,1:ny,inode)%ycoord = ?
ENDDO inode_110

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_surface_relgrd

!========================================================================

SUBROUTINE usrdef_surface_data(iddesc,ifil,ciodatetime,surdata,n1dat,n2dat,&
                             & novars)
!************************************************************************
!
! *usrdef_surface_data* Define surface input data (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.0
!
! Description - called either with iddesc = io_metsur, io_sstsur or io_biosur
!             - the example code considers two methods
!      1: data are obtained from an external input file
!      2: data are obtained from some arbitrary source
!
! Reference -
!
! Calling program - define_surface_data
!
! External calls - 
!
! Module calls - open_filepars
!
!************************************************************************
!
USE iopars
USE syspars
USE inout_routines, ONLY: open_filepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, novars, n1dat, n2dat
REAL, INTENT(INOUT), DIMENSION(n1dat,n2dat,novars) :: surdata

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER  No. of data file
!*ciodatetime* CHAR     Date/time in data file
!*surdata*     REAL     Data array
!*n1dat*       INTEGER  X-dimension of data array
!*n2dat*       INTEGER  Y-dimension of data array
!*novars*      INTEGER  Number of data parameters
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: method = ?
INTEGER :: idgrd, iunit, ivar


procname(pglev+1) = 'usrdef_surface_data'
CALL log_timer_in()

!
!1. Surface grid id
!------------------
!

SELECT CASE (iddesc)
   CASE (io_metsur); idgrd = igrd_meteo
   CASE (io_sstsur); idgrd = igrd_sst
   CASE (io_biosur); idgrd = igrd_bio
END SELECT

!
!2. Data from external data file
!-------------------------------
!

IF (method.EQ.1) THEN

!  ---open input file
   IF (modfiles(iddesc,ifil,1)%status.EQ.0) THEN
      CALL open_filepars(modfiles(io_iddesc,ifil,1))
      GOTO 1000
   ENDIF

!  ---read data
   iunit = modfiles(iddesc,ifil,1)%iunit
   READ (iunit,?,END=99) ciodatetime
   IF (surfacegrids(idgrd,ifil)%nhtype.EQ.0) THEN
      READ (iunit,?) surdata(1,1,1:novars)
   ELSE
      ivar_210: DO ivar=1,novars
         READ (iunit,?) surdata(:,:,ivar)
      ENDDO ivar_210
   ENDIF

!
!3. Data from any source
!-----------------------
!

ELSEIF (method.EQ.2) THEN

!  ---set open status on first call
   IF (modfiles(iddesc,ifil,1)%status.EQ.0) THEN
      modfiles(iddesc,ifil,1)%iostat = 1
      GOTO 1000
   ENDIF

!  ---obtain data
   ciodatetime = ?
   IF (surfacegrids(idgrd,ifil)%nhtype.EQ.0) THEN
      surdata(1,1,1:novars) = ?
   ELSE
      ivar_310: DO ivar=1,novars
         surdata(:,:,ivar) = ?
      ENDDO ivar_310
   ENDIF

!  ---reset status to end of file condition at last call
   modfiles(iddesc,ifil,1)%iostat = 2

ENDIF

1000 CALL log_timer_out()


RETURN

99 modfiles(iddesc,ifil,1)%iostat = 2

END SUBROUTINE usrdef_surface_data
