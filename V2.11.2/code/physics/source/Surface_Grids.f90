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
! *Surface_Grids* Define surface grids
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Grids.f90  V2.11.2
!
! $Date: 2018-05-31 16:41:05 +0200 (Thu, 31 May 2018) $
!
! $Revision: 1141 $
!
! Description - 
!
! Reference -
!
! Subroutines - define_surface_input_grid, define_surface_output_grid,
!               read_surface_absgrd, read_surface_relgrd,
!               write_surface_absgrd, write_surface_relgrd
!
!************************************************************************
!

!========================================================================

SUBROUTINE define_surface_input_grid(idgrd,ifil,surfcoords)
!************************************************************************
!
! *define_surface_input_grid* Define relative coordinates of model grid with
!                             respect to an input surface data grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Grids.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - meteo_input, temperature_equation, wave_input
!
! External calls - read_surface_absgrd, read_surface_relgrd,
!                  usrdef_surface_absgrd, usrdef_surface_relgrd,
!                  write_surface_absgrd , write_surface_relgrd
!
! Module calls - combine_mod_hgrid_2d, construct_rectgrid,
!                distribute_mod_hgrid_2d, error_alloc, error_alloc_struc,
!                model_to_data_hcoords
!  
!************************************************************************
!
USE datatypes
USE gridpars
USE iopars
USE syspars
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE grid_interp, ONLY: model_to_data_hcoords
USE grid_routines, ONLY: construct_rectgrid
USE paral_comms, ONLY: combine_mod_hgrid_2d, distribute_mod_hgrid_2d
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: idgrd, ifil
TYPE (HRelativeCoords), INTENT(OUT), DIMENSION(ncloc,nrloc) :: surfcoords

!
! Name        Type     Purpose
!------------------------------------------------------------------------------
!*idgrd*      INTEGER  Grid id
!*ifil*       INTEGER  No. of grid file
!*surfcoords* DERIVED  Relative coordinates and weight factors for interpolation
!                      of the local model grid on the surface data grid
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, idabs, idrel, j, nhtype, n1dat, n2dat
INTEGER, DIMENSION(4) :: nhdist 
REAL :: delxdat, delydat, x0dat, y0dat
REAL, ALLOCATABLE, DIMENSION(:,:) :: xcoord, ycoord
TYPE (HRelativeCoords), SAVE, ALLOCATABLE, DIMENSION(:,:) :: surfcoordsglb


procname(pglev+1) = 'define_surface_input_grid'
CALL log_timer_in()

!
! Name           Type    Purpose
!------------------------------------------------------------------------------
!*idabs*         INTEGER Absolute grid file id
!*idrel*         INTEGER Relative grid file id
!*nhtype*        INTEGER Type of data grid (0/1/2)
!*n1dat*         INTEGER Number of grid points in X-direction on data grid
!*n2dat*         INTEGER Number of grid points in Y-direction on data grid
!*surfcoordsglb* DERIVED Relative coordinates and weight factors for
!                        interpolation of the local model grid on the surface
!                        data grid
!*xcoord*        REAL    X-coordinates of data grid
!*ycoord*        REAL    Y-coordinates of data grid
!*delxdat*       REAL    Grid spacing in X-direction in case of regular data
!                        grid
!*delydat*       REAL    Grid spacing in Y-direction in case of regular data
!                        grid
!*x0dat*         REAL    X-coordinate of lower left corner in case of regular
!                        data grid
!*y0dat*         REAL    Y-coordinate of lower left corner in case of regular
!                        data grid
!
!------------------------------------------------------------------------------
!
!1. Surface grid id
!------------------
!

SELECT CASE (idgrd)
   CASE (igrd_meteo)
      idabs = io_metabs; idrel = io_metrel
   CASE (igrd_sst)
      idabs = io_sstabs; idrel = io_sstrel
   CASE (igrd_waves)
      idabs = io_wavabs; idrel = io_wavrel
   CASE (igrd_bio)
      idabs = io_bioabs; idrel = io_biorel
END SELECT

!
!2. Return if data and model grid coincide
!------------------------------------------
!

IF (surfacegrids(idgrd,ifil)%nhtype.EQ.0.OR.&
  & surfacegrids(idgrd,ifil)%nhtype.EQ.4) RETURN 

!
!3. Define external grid using absolute coordinates
!--------------------------------------------------
!

IF (surfacegrids(idgrd,ifil)%surcoords.EQ.1) THEN

!
!3.1 Initialise
!--------------
!
!  ---parameters
   n1dat = surfacegrids(idgrd,ifil)%n1dat
   n2dat = surfacegrids(idgrd,ifil)%n2dat
   nhtype = surfacegrids(idgrd,ifil)%nhtype

!  ---allocate coordinate arrays
   ALLOCATE (xcoord(n1dat,n2dat),STAT=errstat)
   CALL error_alloc('xcoord',2,(/n1dat,n2dat/),kndrtype)
   xcoord = 0.0
   ALLOCATE (ycoord(n1dat,n2dat),STAT=errstat)
   CALL error_alloc('ycoord',2,(/n1dat,n2dat/),kndrtype)
   ycoord = 0.0

!
!3.2 Construct data grid
!-----------------------
!

   IF (nhtype.EQ.1) THEN

!
!3.2.1 Uniform rectangular grid
!------------------------------
!
!     ---grid parameters
      x0dat = surfacegrids(idgrd,ifil)%x0dat
      y0dat = surfacegrids(idgrd,ifil)%y0dat
      delxdat = surfacegrids(idgrd,ifil)%delxdat
      delydat = surfacegrids(idgrd,ifil)%delydat

!     ---construct regular grid
      CALL construct_rectgrid(x0dat,y0dat,delxdat,delydat,xcoord,ycoord,&
                            & n1dat,n2dat)

!     ---write to standard file
      IF (modfiles(idabs,ifil,2)%status.EQ.'W') THEN
         CALL write_surface_absgrd(idabs,ifil,n1dat,n2dat,xcoord,ycoord)
      ENDIF

!
!3.2.2 Non-uniform rectangular or curvilinear data grid
!------------------------------------------------------
!

   ELSE

!     ---read from standard file
      IF (modfiles(idabs,ifil,1)%status.EQ.'R') THEN
         CALL read_surface_absgrd(idabs,ifil,n1dat,n2dat,xcoord,ycoord)
!     ---user-defined
      ELSEIF (modfiles(idabs,ifil,1)%status.EQ.'N') THEN
         CALL usrdef_surface_absgrd(idabs,ifil,n1dat,n2dat,xcoord,ycoord)
      ENDIF

!     ---non-uniform rectangular case
      IF (nhtype.EQ.2) THEN
         j_3221: DO j=2,n2dat
            xcoord(:,j) = xcoord(:,1)
         ENDDO j_3221
         i_3222: DO i=2,n1dat
            ycoord(i,:) = ycoord(1,:)
         ENDDO i_3222
      ENDIF

!     ---write to standard file
      IF (modfiles(idabs,ifil,2)%status.EQ.'W') THEN
         CALL write_surface_absgrd(idabs,ifil,n1dat,n2dat,xcoord,ycoord)
      ENDIF

   ENDIF

!
!3.3 Relative coordinates and weight factors
!-------------------------------------------
!

   CALL model_to_data_hcoords(surfacegrids(idgrd,ifil),surfcoords,'C  ',n1dat,&
                            & n2dat,xcoord,ycoord)

!   
!3.4 Write relative coordinates to standard file
!-----------------------------------------------
!   

   IF (modfiles(idrel,ifil,2)%status.EQ.'W') THEN

      ALLOCATE (surfcoordsglb(nc,nr),STAT=errstat)
      CALL error_alloc_struc('surfcoordsglb',2,(/nc,nr/),'HRelativeCoords')
      CALL combine_mod_hgrid_2d(surfcoordsglb,surfcoords,(/1,1/),0,0,0.0)
      CALL write_surface_relgrd(idrel,ifil,surfcoordsglb)
      DEALLOCATE (surfcoordsglb)
      
   ENDIF
   
!
!3.4 Deallocate arrays
!---------------------
!

   DEALLOCATE (xcoord,ycoord)

!
!4. Define external grid using relative coordinates
!--------------------------------------------------
!

ELSEIF (surfacegrids(idgrd,ifil)%surcoords.EQ.2) THEN

!  ---allocate
   ALLOCATE (surfcoordsglb(nc,nr),STAT=errstat)
   CALL error_alloc_struc('surfcoordsglb',2,(/nc,nr/),'HRelativeCoords')
   
!  ---read
   IF (modfiles(idrel,ifil,1)%status.EQ.'R') THEN
      CALL read_surface_relgrd(idrel,ifil,surfcoordsglb,nc,nr,1)
!  ---user-defined
   ELSEIF (modfiles(idrel,ifil,1)%status.EQ.'N') THEN
      CALL usrdef_surface_relgrd(idrel,ifil,surfcoordsglb,nc,nr,1)
   ENDIF

!  ---write
   IF (modfiles(idrel,ifil,2)%defined) THEN
      CALL write_surface_relgrd(idrel,ifil,surfcoordsglb,nc,nr,1)
   ENDIF

!  ---distribute
   nhdist = 0   
   CALL distribute_mod_hgrid_2d(surfcoordsglb,surfcoords,(/1,1/),nhdist,0,0,0.0)

!  ---deallocate
   DEALLOCATE (surfcoordsglb)
   
ENDIF
   
CALL log_timer_out()


RETURN

END SUBROUTINE define_surface_input_grid

!========================================================================

SUBROUTINE define_surface_output_grid(iddesc,ifil,surfcoordsglb,nhdat)
!************************************************************************
!
! *define_surface_output_grid* Relative coordinates and weight factors used to
!                              interpolate external model data to the model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Grids.f90  V2.9
!
! Description - 
!
! Reference -
!
! Calling program -
!
! External calls - read_surface_absgrd, read_surface_relgrd,
!                  usrdef_surface_absgrd, usrdef_surface_relgrd,
!                  write_surface_absgrd , write_surface_relgrd
!
! Module calls - construct_rectgrid, data_to_model_hcoords, error_alloc
!
!************************************************************************
!
USE datatypes
USE iopars
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_interp, ONLY: data_to_model_hcoords
USE grid_routines, ONLY: construct_rectgrid
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, nhdat
TYPE (HRelativeCoords), INTENT(OUT), DIMENSION(nhdat,3) :: surfcoordsglb

!
! Name           Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*        INTEGER  Grid file id
!*ifil*          INTEGER  No. of grid file
!*surfcoordsglb* DERIVED  Relative coordinates and weight factors for
!                         interpolation of the external model grid on the
!                         global model grid
!*nhdat*         INTEGER  Size of data grid
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, idgrd, j, nhtype, n1dat, n2dat
REAL :: delxdat, delydat, x0dat, y0dat
REAL, ALLOCATABLE, DIMENSION(:,:) :: xcoord, ycoord

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*idgrd*       INTEGER Surface grid id
!*nhtype*      INTEGER Type of data grid (0/1/2)
!*n1dat*       INTEGER Number of grid points in X-direction on data grid
!*n2dat*       INTEGER Number of grid points in Y-direction on data grid
!*xcoord*      REAL    X-coordinates of data grid
!*ycoord*      REAL    Y-coordinates of data grid
!*delxdat*     REAL    Grid spacing in X-direction in case of regular data grid
!*delydat*     REAL    Grid spacing in Y-direction in case of regular data grid
!*x0dat*       REAL    X-coordinate of lower left corner in case of regular
!                      data grid
!*y0dat*       REAL    Y-coordinate of lower left corner in case of regular
!                      data grid
!
!------------------------------------------------------------------------------
!
!1. Surface grid id
!------------------
!

SELECT CASE (iddesc)
   CASE (io_metabs,io_metrel); idgrd = igrd_meteo
   CASE (io_sstabs,io_sstrel); idgrd = igrd_sst
   CASE (io_wavabs,io_wavrel); idgrd = igrd_waves
   CASE (io_bioabs,io_biorel); idgrd = igrd_bio
END SELECT

!
!2. No data grid
!---------------
!

nhtype = surfacegrids(idgrd,ifil)%nhtype
IF (nhtype.EQ.0) RETURN 

procname(pglev+1) = 'define_surface_output_grid'
CALL log_timer_in()

!
!3. Uniform data grid 
!--------------------
!

n1dat = surfacegrids(idgrd,ifil)%n1dat
n2dat = surfacegrids(idgrd,ifil)%n2dat

IF (nhtype.EQ.1) THEN

!
!3.1 Grid parameters
!-------------------
!

   x0dat = surfacegrids(idgrd,ifil)%x0dat
   y0dat = surfacegrids(idgrd,ifil)%y0dat
   delxdat = surfacegrids(idgrd,ifil)%delxdat
   delydat = surfacegrids(idgrd,ifil)%delydat

!
!3.2 Allocate arrays
!-------------------
!

   ALLOCATE (xcoord(n1dat,n2dat),STAT=errstat)
   CALL error_alloc('xcoord',2,(/n1dat,n2dat/),kndrtype)

   ALLOCATE (ycoord(n1dat,n2dat),STAT=errstat)
   CALL error_alloc('ycoord',2,(/n1dat,n2dat/),kndrtype)

!
!3.3 Construct regular grid
!--------------------------
!

   CALL construct_rectgrid(x0dat,y0dat,delxdat,delydat,xcoord,ycoord,&
                         & n1dat,n2dat)

!
!3.4 Curvilinear coordinates of data grid
!----------------------------------------
!

   CALL data_to_model_hcoords(surfcoordsglb(:,1),'C  ',n1dat,n2dat,&
                            & xcoord,ycoord,extrapol=.TRUE.,land=.TRUE.)
   CALL data_to_model_hcoords(surfcoordsglb(:,2),'U  ',n1dat,n2dat,&
                            & xcoord,ycoord,extrapol=.TRUE.,land=.TRUE.)
   CALL data_to_model_hcoords(surfcoordsglb(:,3),'V  ',n1dat,n2dat,&
                            & xcoord,ycoord,extrapol=.TRUE.,land=.TRUE.)

!
!3.5 Deallocate arrays
!---------------------
!

   DEALLOCATE (xcoord,ycoord)

!
!4. Rectangular non-uniform data grid
!------------------------------------
!

ELSEIF (nhtype.EQ.2) THEN

!
!4.1 Allocate arrays
!-------------------
!

   ALLOCATE (xcoord(n1dat,n2dat),STAT=errstat)
   CALL error_alloc('xcoord',2,(/n1dat,n2dat/),kndrtype)
   xcoord = 0.0

   ALLOCATE (ycoord(n1dat,n2dat),STAT=errstat)
   CALL error_alloc('ycoord',2,(/n1dat,n2dat/),kndrtype)
   ycoord = 0.0

!
!4.2 Read/define/write grid coordinates
!--------------------------------------
!
!  ---read
   IF (modfiles(iddesc,ifil,1)%status.EQ.'R') THEN
      CALL read_surface_absgrd(iddesc,ifil,n1dat,n2dat,xcoord(:,1),ycoord(1,:))
!  ---user-defined
   ELSEIF (modfiles(iddesc,ifil,1)%status.EQ.'N') THEN
      CALL usrdef_surface_absgrd(iddesc,ifil,n1dat,n2dat,xcoord(:,1),&
                               & ycoord(1,:))
   ENDIF
   j_421: DO j=2,n2dat
      xcoord(:,j) = xcoord(:,1)
   ENDDO j_421
   i_422: DO i=2,n1dat
      ycoord(i,:) = ycoord(1,:)
   ENDDO i_422

!  ---write
   IF (modfiles(iddesc,ifil,2)%defined) THEN
      CALL write_surface_absgrd(iddesc,ifil,n1dat,n2dat,xcoord(:,1),&
                              & ycoord(1,:))
   ENDIF

!
!4.3 Curvilinear coordinate arrays
!---------------------------------
!

   CALL data_to_model_hcoords(surfcoordsglb(:,1),'C  ',n1dat,n2dat,&
                            & xcoord,ycoord,extrapol=.TRUE.,land=.TRUE.)
   CALL data_to_model_hcoords(surfcoordsglb(:,2),'U  ',n1dat,n2dat,&
                            & xcoord,ycoord,extrapol=.TRUE.,land=.TRUE.)
   CALL data_to_model_hcoords(surfcoordsglb(:,3),'V  ',n1dat,n2dat,&
                            & xcoord,ycoord,extrapol=.TRUE.,land=.TRUE.)

!
!4.4 Deallocate arrays
!---------------------
!

   DEALLOCATE (xcoord,ycoord)

!
!5. Curvilinear coordinates from data file
!-----------------------------------------
!

ELSEIF (nhtype.EQ.3) THEN

!  ---read
   IF (modfiles(iddesc,ifil,1)%status.EQ.'R') THEN
      CALL read_surface_relgrd(iddesc,ifil,surfcoordsglb,nhdat,1,3)
!  ---user-defined
   ELSEIF (modfiles(iddesc,ifil,1)%status.EQ.'N') THEN
      CALL usrdef_surface_relgrd(iddesc,ifil,surfcoordsglb,nhdat,1,3)
   ENDIF

!  ---write
   IF (modfiles(iddesc,ifil,2)%defined) THEN
      CALL write_surface_relgrd(iddesc,ifil,surfcoordsglb,nhdat,1,3)
   ENDIF

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE define_surface_output_grid

!========================================================================

SUBROUTINE read_surface_absgrd(iddesc,ifil,n1dat,n2dat,xcoord,ycoord)
!************************************************************************
!
! *read_surface_absgrd* Read absolute coordinate arrays for a surface
!                       data grid in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Grids.f90  V2.1.0
!
! Description - data grid is of type 2
!
! Reference -
!
! Calling program - define_surface_input_grid, define_surface_output_grid
!
! External calls -
!
! Module calls - close_filepars, errror_alloc_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, n1dat, n2dat
REAL, INTENT(OUT), DIMENSION(n1dat,n2dat) :: xcoord
REAL, INTENT(OUT), DIMENSION(n1dat,n2dat) :: ycoord

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER  Grid file id
!*ifil*      INTEGER  No. of grid file
!*n1dat*     INTEGER  X-dimension of data grid
!*n2dat*     INTEGER  Y-dimension of data grid
!*xcoord*    REAL     X-coordinates of data grid
!*ycoord*    REAL     Y-coordinates of data grid
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: numvars
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_surface_absgrd'
CALL log_timer_in()

filepars = modfiles(iddesc,ifil,1)

!
!1. File header
!--------------
!

CALL open_filepars(filepars)
CALL read_glbatts_mod(filepars)
numvars = filepars%novars
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL read_varatts_mod(filepars,varatts,numvars)

!
!2. Read data
!------------
!

CALL read_vars(xcoord,filepars,1,varatts=varatts)
CALL read_vars(ycoord,filepars,2,varatts=varatts)

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(iddesc,ifil,1) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE read_surface_absgrd

!========================================================================

SUBROUTINE read_surface_relgrd(iddesc,ifil,surfcoordsglb)
!************************************************************************
!
! *read_surface_relgrd* Read relative coordinate arrays for a surface
!                       data grid in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Grids.f90  V2.9
!
! Description - 
!
! Reference -
!
! Calling program - define_surface_input_grid, define_surface_output_grid
!
! External calls -
!
! Module calls - close_filepars, error_alloc, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE gridpars
USE iopars
USE syspars
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil
TYPE (HRelativeCoords), INTENT(OUT), DIMENSION(nc,nr) :: surfcoordsglb

! Name           Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*        INTEGER  Grid file id
!*ifil*          INTEGER  No. of grid file
!*surfcoordsglb* DERIVED  Relative coordinates and weight factors for
!                         interpolation between the (global) model grtid and
!                         an external grid
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, inode, j, numvars
TYPE (FileParams) :: filepars
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: weightsdat
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_surface_relgrd'
CALL log_timer_in()

filepars = modfiles(iddesc,ifil,1)

!
!1. File header
!--------------
!

CALL open_filepars(filepars)
CALL read_glbatts_mod(filepars)
numvars = filepars%novars
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL read_varatts_mod(filepars,varatts,numvars)

!
!2. Allocate array with weights data
!-----------------------------------
!

ALLOCATE(weightsdat(2,2,nc,nr),STAT=errstat)
CALL error_alloc('weightsdat',4,(/2,2,nc,nr/),kndrtype)

!
!3. Read data
!------------
!

CALL read_vars(surfcoordsglb%icoord,filepars,1,varatts=varatts)
CALL read_vars(surfcoordsglb%jcoord,filepars,2,varatts=varatts)
CALL read_vars(weightsdat,filepars,3,varatts=varatts)   
j_310: DO j=1,nr
i_310: DO i=1,nc
   surfcoordsglb(i,j)%weights = weightsdat(:,:,i,j)
ENDDO i_310
ENDDO j_310

!
!4. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(iddesc,ifil,1) = filepars
DEALLOCATE (varatts,weightsdat)


CALL log_timer_out()


RETURN

END SUBROUTINE read_surface_relgrd

!========================================================================

SUBROUTINE write_surface_absgrd(iddesc,ifil,n1dat,n2dat,xcoord,ycoord)
!************************************************************************
!
! *write_surface_absgrd* Write absolute coordinate arrays for a surface
!                        data grid in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Grids.f90  V2.0
!
! Description - data grid is of type 2
!
! Reference -
!
! Calling program - define_surface_input_grid, define_surface_output_grod
!
! External calls -
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, write_atts_mod,
!                write_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE paralpars
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, n1dat, n2dat
REAL, INTENT(IN), DIMENSION(n1dat,n2dat) :: xcoord
REAL, INTENT(IN), DIMENSION(n1dat,n2dat) :: ycoord

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Grid file id
!*ifil*        INTEGER  No. of grid file
!*n1dat*       INTEGER  X-dimension of data grid
!*n2dat*       INTEGER  Y-dimension of data grid
!*xcoord*      REAL     X-coordinates of data grid
!*ycoord*      REAL     Y-coordinates of data grid
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: numvars
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_surface_absgrd'
CALL log_timer_in()

!
!1. Write file header
!--------------------
!
!---file attributes
CALL set_modfiles_atts(iddesc,ifil,2)
filepars = modfiles(iddesc,ifil,2)
numvars = filepars%novars

!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL set_modvars_atts(iddesc,ifil,2,varatts,numvars)

!---write
CALL write_atts_mod(filepars,varatts,numvars)

!
!2. Write data
!-------------
!

CALL write_vars(xcoord,filepars,1,varatts=varatts)
CALL write_vars(ycoord,filepars,2,varatts=varatts)

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(iddesc,ifil,2) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE write_surface_absgrd

!========================================================================

SUBROUTINE write_surface_relgrd(iddesc,ifil,surfcoordsglb)
!************************************************************************
!
! *write_surface_relgrd* Write relative coordinate arrays for a surface
!                        data grid in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Grids.f90  V2.9
!
! Description - 
!
! Reference -
!
! Calling program - define_surface_input_grid, define_surface_output_grid
!
! External calls -
!
! Module calls - close_filepars, error_alloc, error_alloc_struc, open_filepars
!                set_modfiles_atts, set_modvars_atts, write_atts_mod, write_vars
!
!************************************************************************
!
USE datatypes
USE gridpars
USE iopars
USE paralpars
USE syspars
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil
TYPE (HRelativeCoords), INTENT(IN), DIMENSION(nc,nr) :: surfcoordsglb

!
! Name           Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*        INTEGER  Grid file id
!*ifil*          INTEGER  No. of grid file
!*surfcoordsglb* DERIVED  Relative coordinates and weight factors for
!                         interpolation between the (global) model grtid and an
!                         external grid
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, inode, j, numvars
TYPE (FileParams) :: filepars
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: weightsdat
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_surface_relgrd'
CALL log_timer_in()

!
!1. Write file header
!--------------------
!
!---file attributes
CALL set_modfiles_atts(iddesc,ifil,2)
filepars = modfiles(iddesc,ifil,2)
numvars = filepars%novars

!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL set_modvars_atts(iddesc,ifil,2,varatts,numvars)

!---write
CALL write_atts_mod(filepars,varatts,numvars)

!
!2. Allocate array with weights data
!-----------------------------------
!

ALLOCATE(weightsdat(2,2,nc,nr),STAT=errstat)
CALL error_alloc('weightsdat',4,(/2,2,nc,nr/),kndrtype)

!
!3. Write data
!-------------
!

CALL write_vars(surfcoordsglb%icoord,filepars,1,varatts=varatts)
CALL write_vars(surfcoordsglb%jcoord,filepars,2,varatts=varatts)
j_310: DO j=1,nr
i_310: DO i=1,nc
   weightsdat(:,:,i,j) = surfcoordsglb(i,j)%weights
ENDDO i_310
ENDDO j_310
CALL write_vars(weightsdat,filepars,3,varatts=varatts)

!
!4. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(iddesc,ifil,2) = filepars
DEALLOCATE (varatts,weightsdat)

CALL log_timer_out()


RETURN

END SUBROUTINE write_surface_relgrd
