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
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.11
!
! $Date: 2013-04-16 12:20:30 +0200 (Tue, 16 Apr 2013) $
!
! $Revision: 548 $
!
! Description - test case partrhone
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
USE datatypes

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
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.11
!
! Description - test case partrhone
!
! Reference -
!
! Calling program - define_surface_data
!
!************************************************************************
!
USE iopars
USE syspars
USE timepars
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
REAL :: uwindatu_cst, vwindatv_cst, wdir, wind


procname(pglev+1) = 'usrdef_surface_data'
CALL log_timer_in()

!---set open condition on first call
IF (modfiles(io_metsur,1,1)%iostat.EQ.0) THEN
   modfiles(io_metsur,1,1)%iostat = 1
   GOTO 1000
ENDIF

!---define data
ciodatetime = CDateTime
IF (nt.EQ.0) THEN
   uwindatu_cst = 0.0; vwindatv_cst = 0.0
ELSE
   wind = 10.0; wdir = 90.0
   wdir = wdir*degtorad
   uwindatu_cst = wind*COS(wdir)
   vwindatv_cst = wind*SIN(wdir)
ENDIF
surdata(:,:,1) = uwindatu_cst
surdata(:,:,2) = vwindatv_cst

1000 CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_surface_data
