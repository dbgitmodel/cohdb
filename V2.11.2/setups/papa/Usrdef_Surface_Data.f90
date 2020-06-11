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
! *Usrdef_Surface_Data* User-defined surface data setup
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.11.2
!
! $Date: 2013-05-28 15:26:47 +0200 (Tue, 28 May 2013) $
!
! $Revision: 574 $
!
! Description - test case papa
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
!                         with respect to model grid
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
!************************************************************************
!

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


RETURN

END SUBROUTINE usrdef_surface_absgrd

!========================================================================

SUBROUTINE usrdef_surface_relgrd(iddesc,ifil,surfgridglb,nx,ny,nonodes)
!************************************************************************
!
! *usrdef_surface_relgrd* Define relative coordinate array of surface grid(s)
!                         with respect to model grid
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
!************************************************************************
!
USE datatypes

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


RETURN

END SUBROUTINE usrdef_surface_relgrd

!========================================================================

SUBROUTINE usrdef_surface_data(iddesc,ifil,ciodatetime,surdata,n1dat,n2dat,&
                             & novars)
!************************************************************************
!
! *usrdef_surface_data* Define surface input data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.11.2
!
! Description - test case papa
!
! Reference -
!
! Calling program - define_surface_data
!
! Module calls - convert_date, open_filepars 
!
!************************************************************************
!
USE iopars
USE syspars
USE inout_routines, ONLY: open_filepars
USE time_routines, ONLY: convert_date, log_timer_in, log_timer_out

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
INTEGER, SAVE :: iunit
REAL :: airtemp_cst, cloud_cover_cst, relhum_cst, uwindatu_cst, vwindatv_cst

procname(pglev+1) = 'usrdef_surface_data'
CALL log_timer_in()

!
!1. Open file
!------------
!

IF (modfiles(io_metsur,1,1)%iostat.EQ.0) THEN
   CALL open_filepars(modfiles(io_metsur,1,1))
   iunit = modfiles(io_metsur,1,1)%iunit
   modfiles(io_metsur,1,1)%iostat = 1
   GOTO 1000
ENDIF

!
!2. Read data
!------------
!

ciodatetime(20:23) = ':000'
READ (iunit,9001) ciodatetime(1:19), uwindatu_cst, vwindatv_cst, airtemp_cst, &
                & relhum_cst, cloud_cover_cst

surdata(1,1,1) = uwindatu_cst; surdata(1,1,2) = vwindatv_cst
surdata(1,1,3) = airtemp_cst; surdata(1,1,4) = relhum_cst
surdata(1,1,5) = cloud_cover_cst

1000 CALL log_timer_out()


RETURN

9001 FORMAT(A,2(1X,F6.2),1X,F5.1,2(1X,F5.3))

END SUBROUTINE usrdef_surface_data
