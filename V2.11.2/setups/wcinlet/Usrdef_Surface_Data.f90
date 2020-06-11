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
! Author - Boudewijn Decrop
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.11.1
!
! $Date: 2018-01-16 16:50:11 +0100 (Tue, 16 Jan 2018) $
!
! $Revision: 1071 $
!
! Description - test case wavload
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
! Description - test case tidal inlet with bottom current/wave interaction
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
! Author - Boudewijn Decrop
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.11.1
!
! Description - test case tidal inlet with bottom current/wave interaction
!
! Reference -
!
! Calling program - add_secs_to_date, define_surface_data
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE syspars
USE wavevars
USE cf90_routines, ONLY: cf90_get_var, cf90_get_var_chars, cf90_inq_dimid, &
                       & cf90_inquire, cf90_inquire_dimension, &
                       & cf90_inquire_variable
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: close_filepars, open_filepars
USE math_library, ONLY: complex_polar
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
INTEGER, SAVE :: irec, iunit, maxrec, nvars
INTEGER :: ivar
CHARACTER (LEN=lenname), SAVE, ALLOCATABLE, DIMENSION(:) :: varname
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: wavedircos, wavedirsin


procname(pglev+1) = 'usrdef_surface_data'
CALL log_timer_in()

!
!1. Initialisation on first calls
!---------------------------------
!

IF (modfiles(io_wavsur,1,1)%iostat.EQ.0) THEN

!  ---open data file
   CALL open_filepars(modfiles(io_wavsur,1,1))
   iunit = modfiles(io_wavsur,1,1)%iunit

!  ---record number
   irec = 1

!  ---number of variables
   CALL cf90_inquire(iunit,nvariables=nvars)

!  ---allocate arrays
   ALLOCATE(varname(nvars),STAT=errstat)
   CALL error_alloc('varname',1,(/nvars/),kndchar,lenname)
   ALLOCATE(wavedircos(nc,nr),STAT=errstat)
   CALL error_alloc('wavedircos',2,(/nc,nr/),kndrtype)
   ALLOCATE(wavedirsin(nc,nr),STAT=errstat)
   CALL error_alloc('wavedirsin',2,(/nc,nr/),kndrtype)

!  ---variable names
   ivar_110: DO ivar=1,nvars
      CALL cf90_inquire_variable(iunit,ivar,name=varname(ivar))
   ENDDO ivar_110

!  ---time records
   CALL cf90_inquire_dimension(iunit,1,len=maxrec)
   
   GOTO 1000

ENDIF

!
!2. Wave data
!------------
!

ivar_210: DO ivar=1,nvars
   SELECT CASE(TRIM(varname(ivar)))
      CASE ('time'); CALL cf90_get_var_chars(iunit,1,ciodatetime,irec)
      CASE ('waveheight'); CALL cf90_get_var(iunit,ivar,surdata(:,:,1),irec)
      CASE ('waveperiod'); CALL cf90_get_var(iunit,ivar,surdata(:,:,2),irec)
      CASE ('wavedircos'); CALL cf90_get_var(iunit,ivar,wavedircos,irec)
      CASE ('wavedirsin'); CALL cf90_get_var(iunit,ivar,wavedirsin,irec)
   END SELECT
ENDDO ivar_210

!----wave direction
CALL complex_polar(wavedircos,wavedirsin,xpha=surdata(:,:,3),&
                 & maskvals=maskatc_int)

!---update record number
irec = irec + 3

!
!3. Finalise
!-----------
!

IF (irec.GT.maxrec) THEN

!  ---close datafile
   CALL close_filepars(modfiles(io_wavsur,1,1))

!  ---deallocate
   DEALLOCATE (varname,wavedircos,wavedirsin)

ENDIF

1000 CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_surface_data
