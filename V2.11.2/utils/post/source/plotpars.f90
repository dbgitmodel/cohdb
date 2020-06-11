MODULE plotpars
!************************************************************************
!
! *plotpars* Parameters and arrays for postprocessor
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)plotpars.f90  V2.11.2
!
! $Date: 2018-07-23 16:55:25 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1171 $
!
! Description -
!
!************************************************************************
!
USE datatypes
USE syspars

IMPLICIT NONE


INTEGER, PARAMETER :: MaxContLevs = 30, MaxTrack = 100
INTEGER, PARAMETER :: Maxpars = MAX(30,MaxTrack), MaxString = 120

TYPE :: PostParams
   LOGICAL :: anim, contreg
   CHARACTER (LEN=lentime) :: PDateTime
   INTEGER :: iend, ipos, istart, jend, jpos, jstart, horz_unit, icontstyle, &
            & iplotform, kpos, linepsyms, nolevels, nolists, notracks, &
            & novarskip1, novarskip2, novecskip1, novecskip2, numfig, numstat, &
            & numvar, numvec1, numvec2, numvec3, vert_unit
   INTEGER, DIMENSION(MaxTrack) :: linestyles, ptrack
   REAL :: contmax, contmin, depmax, depmin, depplot, gxend, gxpos, gxstart, &
         & gyend, gypos, gystart, hlen_unit, hrefvec, vrefvec
   REAL, DIMENSION(MaxContLevs) :: contlevels
END TYPE PostParams
TYPE (PostParams), ALLOCATABLE, DIMENSION(:) :: postatts
TYPE (FileParams) :: filepars
TYPE (OutGridParams) :: outgpars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: coordvars, outvars

LOGICAL :: coast_def, fullinfo = .true., gridded, packing, time_grid, traject, &
         & ztime
CHARACTER (LEN=lentime) :: PEndDateTime, PRefDateTime, PStartDateTime
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=2) :: plottype
CHARACTER (LEN=leniofile) :: coast_name, outfile, plotfile, postfile
CHARACTER (LEN=lendesc), DIMENSION(4) :: plottitles
CHARACTER (LEN=*), PARAMETER, DIMENSION(8) :: &
          & time_unit = (/'seconds','minutes','hours  ','days   ','weeks  ',&
                       &  'months ','years  ','year   '/)
INTEGER :: io_filevis, io_parvis, ncout, nocoords, nodim, nofigs, nolists, &
         & nopartout, noplots, norecskip, noseries, nostats, novars, nrout, &
         & nstepout, nzout, plinenum, ptime_format, gtime_format, vcoord
INTEGER, ALLOCATABLE, DIMENSION(:) :: lmaskout, vecids
REAL :: border = 0.1
REAL (KIND=kndrlong) :: deltout, fill_value, fill_value_eps, offset
REAL, DIMENSION(4) :: plotcorners

CHARACTER (LEN=lendesc), ALLOCATABLE, DIMENSION(:) :: station_names
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: gsout, gzcoord, gzcoord_ini, outvals3d
REAL, ALLOCATABLE, DIMENSION(:) :: outvals1d
REAL, ALLOCATABLE, DIMENSION(:,:) :: depout, outvals2d
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: outvals4d


SAVE

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*MaxContLevs*   INTEGER Maximum allowed levels for contouring
!*PostParams*    DERIVED Derived type definition for plotting attributes
!*%anim*         LOGICAL Enables/disables animation for 'HT' and 'VT' plots
!*%contreg*      LOGICAL Uniform spacings for contouring if .TRUE.
!*%PDateTime*    CHAR    Date/time for  'ZP' and 'HP' plots
!*%iend*         INTEGER X-index of the upper right location for the data plot   
!*%ipos*         INTEGER Index X-coordinate of plot location
!                        ('TS','TZ','ZP' plots)
!*%istart*       INTEGER X-index of the lower left location for the data plot
!*%jend*         INTEGER Y-index of the upper right location for the data plot   
!*%jpos*         INTEGER Index Y-coordinate of plot location
!                        ('TS','TZ','ZP' plots)
!*%jstart*       INTEGER Y-index of the lower left location for the data plot
!*%horz_unit*    INTEGER Defines how to interpret horizontal coordinates in
!                        parameter file
!               (parameters gxend, gxpos, gxstart, gyend, gypos, gystart)
!                 = 1 => Index coordinates (with respect to output grid) 
!                 = 2 => Geographical coordinates
!*%icontstyle*   INTEGER Style for contouring
!                 = 1 => isolines annotated and drawn, no filling
!                 = 2 => isolines not annotated, not drawn, with filling
!                 = 3 => isolines not annotated, drawn, with filling
!*%iplotform*    INTEGER Selects how to use vertical data in plots without a
!                        vertical component
!                 = 1 => at vertical index coordinate 'kpos'
!                 = 2 => at surface distance 'depplot'
!                 = 3 => depth-averaged profile value
!                 = 4 => depth-integrated profile value
!                 = 5 => maximum value over the vertical profile
!*%kpos*         INTEGER Vertical index coordinate used if iplotform=1
!*%linepsyms*    INTEGER Linestyles used for curves in 'muplot'
!*%nolevels*     INTEGER Number of levels for contouring
!*%nolists*      INTEGER Number of plot lists
!*%notracks*     INTEGER Number of output trajectories
!*%novarskip1*   INTEGER Sampling interval for scalar output data in the
!                        X-direction
!*%novarskip2*   INTEGER Sampling interval for scalar output data in the
!                        Y-direction
!*%novecskip1*   INTEGER Sampling interval for vector output data in the
!                        X-direction
!*%novecskip2*   INTEGER Sampling interval for vector output data in the
!                        Y-direction
!*%numfig*       INTEGER Plot number where a curve is drawn
!*numpart*       INTEGER 
!*%numstat*      INTEGER Station number for non-gridded output data
!*%numvar*       INTEGER Variable id selected for plotting
!*%numvec1*      INTEGER Variable id of the X-component of the vectors used in
!                        vector plotting
!*%numvec2*      INTEGER Variable id of the Y-component of the vectors used in
!                        vector plotting
!*%numvec3*      INTEGER Variable id of the vertical component of the vectors
!                        used in vector plotting
!*%ptrack*       INTEGER Indices of the output particles used for plotting 
!*%vert_unit*    INTEGER Defines how to interpret the parameters 'depmin' and
!                        'depmax' in parameter file
!                 = 1 => as vertical index coordinates
!                 = 2 => as distances from the (mean) surface
!*%contmax*      REAL    Maximum contour value
!*%contmin*      REAL    Minimum contour value
!                (if contmin and contmax are zero, they are taken as minimum,
!                 resp maximum of the data plot array) 
!*%depmax*       REAL    Sets maximum water depth in data plot. Value is
!                        interpreted as max. vertical index if vert_unit=1, or
!                        max. depth if vert_unit=2
!*%depmin*       REAL    Sets minimum water depth in data plot. Value is
!                        interpreted as min. vertical index if vert_unit=1, or
!                        min. depth if vert_unit=2
!*%depplot*      REAL    Surface distance (in meter) used if iplotform=2
!*%gxend*        REAL    End X-position for horizontal data plots
!*%gxpos*        REAL    X-position for non-horizontal data plots
!*%gxstart*      REAL    Start X-position for horizontal data plots
!*%gyend*        REAL    End Y-position for horizontal data plots
!*%gypos*        REAL    Y-position for non-horizontal data plots
!*%gystart*      REAL    Start Y-position for horizontal data plots
!*%hlen_unit*    REAL    Scaling factor for X- and Y-coordinates (Cartesian
!                        grid only)
!*%hrefvec*      REAL    Reference value for horizontal vectors used in plot
!                        legend 
!*%vrefvec*      REAL    Reference value for vertical vectors used in plot
!                        legend
!*%contlevels*   REAL    Contour levels
!*postatts*      DERIVED Attributes for plotting
!*filepars*      DERIVED Attributes of output data file
!*gridpars*      DERIVED Atrributes of output grid file
!*outgpars*      DERIVED Attributes of output grid
!*coordvars*     DERIVED Attributes of output grid coordinates
!*outvars*       DERIVED Attributes of output data variables
!*coast_def*     LOGICAL Coast file is determined internally (externally) if
!                        .TRUE. (.FALSE.)
!*fullinfo*      LOGICAL Write info all data reads if .TRUE.
!*gridded*       LOGICAL .TRUE. (.FALSE.) for regularly (irregularly) gridded
!                        output data
!*packing*       LOGICAL .TRUE. (.FALSE.) if output data are stored in packed
!                        format
!*plinenum*      INTEGER Line number on parameter input file
!*time_grid*     LOGICAL .TRUE. (.FALSE.) if output grid is dependent (not
!                        dependent) on time
!*traject*       LOGICAL Output file contains particle trajects if .TRUE.
!*ztime*         LOGICAL If .TRUE. and time_grid is .TRUE., data plots include
!                        time-dependent vertical coordinates
!*PEndDateTime*  CHAR    End date/time used for data plots
!*PRefDateTime*  CHAR    Reference date/time to define coordinates in time axis
!*PStartDateTime*CHAR    Start date/time for data plots
!*outform*       CHAR    Output file format ('A', 'U', 'N')
!*plottype*      CHAR    Plottype
!                  = 'HT' => horizontal transects
!                  = 'VT' => vertical transects
!                  = 'TS' => time series
!                  = 'TZ' => time-depth contours
!                  = 'TH' => transect-time contours
!                  = 'ZP' => vertical profiles
!                  = 'HP' => horizontal profiles 
!*coast_name*    CHAR    Name of coast file
!*outfile*       CHAR    Name of output data file
!*plotfile*      CHAR    String to construct names of plotting files
!*postfile*      CHAR    Name of parameter file
!*plottitles*    CHAR    Plot titles
!*time_unit*     CHAR    Time unit
!*io_filevis*    INTEGER File unit of 'files.vis' file
!*io_parvis*     INTEGER File unit of '.vis' file
!*ncout*         INTEGER X-dimension of output data
!*nocoords*      INTEGER Number of coordinate variables in data file
!*nodim*         INTEGER Spatial dimension of output grid
!*nofigs*        INTEGER Number of plots for each parameter listing and plot
!                        type
!*nolists*       INTEGER Number of plot specifiers for a given plot type
!                        (i.e. size of array postatts)
!*noplots*       INTEGER Total number of plots per plot type
!*norecskip*     INTEGER Sampling interval for output data in time
!*noseries*      INTEGER Number of time specifiers for each plot type
!*nostats*       INTEGER Number of output stations (irregular data grid)
!*novars*        INTEGER Number of output data variables
!*nrout*         INTEGER Y-dimension of output data
!*nstepout*      INTEGER Number of time records in data file
!*nzout*         INTEGER Vertical dimension of output data
!*ptime_format*  INTEGER Type of time coordinate for plotting
!                 = 1 => seconds
!                 = 2 => minutes
!                 = 3 => hours
!                 = 4 => days
!                 = 5 => months
!                 = 6 => years
!                 = 7 => date in years
!*gtime_format*  INTEGER Time format of output data
!                 = 0 => date/time in string format
!                 = 1 => seconds
!                 = 2 => minutes
!                 = 3 => hours
!                 = 4 => days
!                 = 5 => months
!                 = 6 => years
!                 = 7 => date in years
!*vecids*        INTEGER IDs of non-coordinate variables
!*border*        REAL    Fractional size of border around all 'mumap' plot
!*deltout*       REAL    Output time step in seconds
!*offset*        REAL    Offset (in seconds) of plot reference date with
!                        respect to output reference date
!*plotcorners*   REAL    Coordinates of plot borders
!*station_names* CHAR    Station names
!*depout*        REAL    Output bathymetry
!*gzcoord*       REAL    Mean Z-coordinates of output data
!*outvals1d*     REAL    Output data stored as a vector array
!*outvals2d*     REAL    Output data stored as a 2-D array
!*outvals3d*     REAL    Output data stored as a 3-D array
!*outvals4d*     REAL    Output data stored as a 4-D array
!
!------------------------------------------------------------------------------
!

END MODULE plotpars
