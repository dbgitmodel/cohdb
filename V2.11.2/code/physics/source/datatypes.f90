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

MODULE datatypes
!************************************************************************
!
! *datatypes* Derived type definitions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

!---parameters for halo communications in parallel mode
TYPE :: ExchComms
   LOGICAL :: sfirst
   INTEGER :: iddest, idsrce, tag
   INTEGER, DIMENSION(2) :: i1rcv, i2rcv, j1rcv, j2rcv
   INTEGER, DIMENSION(2) :: i1snd, i2snd, j1snd, j2snd
END TYPE ExchComms

!---file attributes for model I/O
TYPE :: FileParams
   LOGICAL :: defined, fill, info, opened, packing, time_regular
   CHARACTER (LEN=1) :: floattype, status
   CHARACTER (LEN=lenform) :: form
   CHARACTER (LEN=lencomm) :: comment
   CHARACTER (LEN=leniofile) :: filename, pathname
   CHARACTER (LEN=lendesc) :: title
   INTEGER :: endfile, iostat, iunit, maxrecs, nocoords, nodim, &
            & novars, nstepout, timeid, timerec, tskips, varid, zvarid
   INTEGER, DIMENSION(3) :: tlims
END TYPE FileParams

!---Global attributes for files in standard (CF-1.6) format
TYPE :: GlbFileAtts
   CHARACTER (LEN=lendesc) :: name, value
END TYPE GlbFileAtts

!---attributes of surface grids
TYPE :: GridParams
   LOGICAL :: extrapol, land, rotated
   INTEGER :: datflag, nhtype, n1dat, n2dat, surcoords
   REAL :: delxdat, delydat, gridangle, longpole, x0dat, y0dat, x0rot, y0rot
END TYPE GridParams

!---horizontal relative coordinates
TYPE :: HRelativeCoords
   INTEGER :: icoord, jcoord
   REAL, DIMENSION(2,2) :: weights
END TYPE HRelativeCoords

!---space/time locations for user-defined output
TYPE :: OutGridParams
   LOGICAL :: bedvars, dredgevars, gridded, packing, relocatevars, time_grid
   CHARACTER (LEN=1) :: status
   CHARACTER (LEN=lentime) :: enddate, refdate, startdate
   INTEGER :: ncout, nodim, nostats, nowetout, nrout, nstepout, nzout, &
            & time_format, vcoord
   INTEGER, DIMENSION(3) :: tlims, xlims, ylims, zlims
   REAL :: deltout
END TYPE OutGridParams

!---station locations for irregular user-defined output
TYPE :: StationLocs
   INTEGER :: ipos, jpos
   CHARACTER (LEN=lendesc) :: name
END TYPE StationLocs

!---variable attributes
TYPE :: VariableAtts
   LOGICAL :: fill
   CHARACTER (LEN=1) :: axis
   CHARACTER (LEN=lenname) :: f90_name
   CHARACTER (LEN=lendesc) :: comment, coordinates, long_name, standard_name, &
                            & vector_name
   CHARACTER (LEN=lenunit) :: units
   CHARACTER (LEN=lennode) :: node
   INTEGER :: data_type, ivarid, klev, nrank, numvar, oopt
   INTEGER, DIMENSION(4) :: dimids
   INTEGER, DIMENSION(4) :: halo_dims
   INTEGER, DIMENSION(5) :: global_dims, local_dims
   REAL :: dep
   REAL (KIND=kndrlong) :: fill_value
END TYPE VariableAtts

!---vertical relative coordinates
TYPE :: VRelativeCoords
   INTEGER :: kcoord
   REAL, DIMENSION(2) :: weights
END TYPE VRelativeCoords

!
! Name              Type    Purpose
!------------------------------------------------------------------------------
!*ExchComms*        DERIVED Parameters for halo communications in parallel mode
!*%sfirst*          LOGICAL Switch which determines if send precedes receive in
!                           non-blocking communications 
!*%iddest*          INTEGER Destination process of send operation
!*%idsrce*          INTEGER Source process of receive operation
!*%tag*             INTEGER Communication tag
!*%i1rcv*           INTEGER Start index in X-direction for receive operation
!*%i2rcv*           INTEGER End index in X-direction for receive operation
!*%j1rcv*           INTEGER Start index in Y-direction for receive operation
!*%j2rcv*           INTEGER End index in Y-direction for receive operation
!*%i1snd*           INTEGER Start index in X-direction for send operation
!*%i2snd*           INTEGER End index in X-direction for send operation
!*%j1snd*           INTEGER Start index in Y-direction for send operation
!*%j2snd*           INTEGER End index in Y-direction for send operation
!
!*FileParams*       DERIVED Attributes of model and user output files
!*%defined%         LOGICAL .TRUE. if %status is not equal to '0'
!*%fill*            LOGICAL Allows the use of fill values for data at dry
!                           points if .TRUE.
!*%info*            LOGICAL Header information is separately written to an
!                           information file if .TRUE.
!*%opened*          LOGICAL .TRUE. (.FALSE.) if file is connected (not
!                           connected)
!*%packing*         LOGICAL Data are stored in condensed "packing" format (i.e.
!                           vector format omitting land points) if .TRUE.
!*%time_regular*    LOGICAL File data or regular/irregular in time if 
!                           .TRUE./.FALSE.
!*%floattype*       CHAR    Type of real file data
!                     = 'S' -=> single precision
!                     = 'D' => double precision
!*%form*            CHAR    File format ('A', 'U', 'N')
!*%status*          CHAR    Determines whether file data are read,
!                           user_defined or written or temporary ('A','N','W')
!                     = '0' => undefined (all files)
!                     = 'R' => file data are read in COHERENS standard format
!                              (input only)
!                     = 'N' => file data are user-obtained (input only)
!                     = 'W' => file data are written to a newly created file in
!                              COHERENS standard format (output only)
!*%filename*        CHAR    File name
!*%pathname*        CHAR    Preappended file path
!*%title*           CHAR    A succinct description of what is in the file
!*%institution*     CHAR    Specifies where the original data was produced
!*%history*         CHAR    Date/time when the file was created
!*%source*          CHAR    Coherens version
!*%Conventions*     CHAR    netCDF convention
!*%comment*         CHAR    Miscellaneous information about the data or methods
!*%endfile*         INTEGER Switch to select action in case of an end of file
!                          condition
!                     = 0 => abort program (default)
!                     = 1 => program continues without reading data
!                     = 2 => program waits until data can be read
!*%iostat*          INTEGER File status
!                     =-1 => unable to open file
!                     = 0 => at the start of the file or file not opened
!                     = 1 => in the middle of the file
!                     = 2 => after I/O of last data
!                     = 3 => end of file
!*%iunit*           INTEGER File unit number
!*%maxrecs*         INTEGER Total amount of time records
!*%nocoords*        INTEGER Number of coordinate arrays
!*%nodim*           INTEGER Variable ranks (user output only)
!*%novars*          INTEGER Number of variables
!*%nstepout*        INTEGER Number of time records
!*%timeid*          INTEGER netcdf varid for time variable
!*%timerec*         INTEGER Time record number
!*%tskips*          INTEGER Selects data to be skipped in output mode
!                     = 1 => data are written for all times
!                     = 2 => t_mod =< t_obs
!                     = 3 => t_mod  < t_obs
!                     = 4 => t_obs =< t_end
!                     = 5 => t_mod =< t_obs =< t_end
!                     = 6 => t_mod  < t_obs =< t_end
!                     = 7 => t_obs  < t_end
!                     = 8 => t_mod =< t_obs < t_end
!                     = 9 => t_mod  < t_obs < t_end
!               (option only used in absence of other time selection mechanism:
!                initial condition files, nested output, user output) 
!*%varid*           INTEGER Variable record number
!*%zvarid*          INTEGER Variable file ID for the time-varying vertical
!                           coordinate
!*%nstepout*        INTEGER Number of time step for output
!*%tlims*           INTEGER Start/End/Stride time indices (model files only)
!
!*GlbFileAtts*      DERIVED Global attributes for files in COHERENS standard
!                           (CF-1.6) format
!*%name*            CHAR    Attribute's name
!*%value*           CHAR    Atribute's value
!                
!*GridParams*       DERIVED Attributes of surface data grids
!*%extrapol*        LOGICAL Allows extrapolation in case model grid points are outside
!                           the surface data grid
!*%land*            LOGICAL Allows extrapolation in case model grid points by land
!                           cells on the surface data grid
!*%rotated*         LOGICAL Selects rotated rectangular grid (nhtype=1,2)
!*%datflag*         INTEGER Selects type of interpolation when the surface data
!                           contains flagged values
!                     = 0 => no flags
!                     = 1 => data may contain flagged values, interpolation is
!                            performed if all surrounding data values are not
!                            flagged
!                     = 2 => data may contain flagged values, interpolation is
!                            performed if at least one surrounding data value
!                            is not flagged
!*%nhtype*          INTEGER Horizontal structure of model or data grid
!                     = 0 => single point data grid or data grid not defined
!                     = 1 => rectangular grid with uniform grid spacings in
!                            Cartesian or spherical coordinates (data grid
!                            constructed from x0dat, y0dat, delxdat, delydat)
!                     = 2 => rectangular grid coordinates with non-uniform
!                            spacings
!                     = 3 => non-rectangular (curvilinear or non-structured)
!                            grid
!                     = 4 => data grid coincides with model grid
!*%n1dat*           INTEGER Number of grid points in X-direction or number of
!                           surface points for resp. gridded, non-gridded data
!*%n2dat*           INTEGER Number of grid points in Y-direction or equal to 1
!                           for resp. gridded, non-gridded data
!*%surcoords*       INTEGER Type of input coordinates used for interpolation
!                           from data grid
!                      = 1 => using absolutes coordinates
!                      = 2 => using relative coordinates
!*%delxdat*         REAL    Spacing in X-direction in case of uniform
!                           rectangular grid
!*%delydat*         REAL    Spacing in Y-direction in case of uniform
!                           rectangularr grid
!*%gridangle*       REAL    Rotation angle in case of rotated grid     [degrees]
!*%longpole*        REAL    Longitude of rotated North pole            [degrees]
!*%x0dat*           REAL    X-coordinate of lower left Southwest corner in the
!                           (non-rotated) grid                    [m or degrees]
!*%x0rot*           REAL    Longitude of Southwest corner in the rotated grid
!                                                                      [degrees]
!*%y0dat*           REAL    Y-coordinate of lower left (Southwest) corner in
!                           the (non-rotated) grid                [m or degrees]
!*%y0rot*           REAL    Latitude of Southwest corner in the rotated grid
!                                                                      [degrees]
!
!*HRelativeCoords*  DERIVED Horizontal relative coordinates
!*%icoord*          INTEGER X-index coordinate with respect to lower left
!                           corner of reference grid
!*%jcoord*          INTEGER Y-index coordinate with respect to lower left
!                           corner of reference grid
!*%weights*         REAL    Weight factors for horizontal interpolation
!
!*OutGridParams*    DERIVED Attributes of user output grid
!*%bedvars          LOGICAL .TRUE. (.FALSE.) when output variables are (not) 
!                             characteristics of the bed
!*%gridded*         LOGICAL .TRUE. (.FALSE.) for horizontally regular
!                           (irregular) output
!*%packing*         LOGICAL Data are stored in condensed "packing" format (i.e.
!                           vector format omitting land points) if .TRUE.
!*%time_grid*       LOGICAL Create time-dependent grid if .TRUE. and nodim = 3
!*%status*          CHAR    Status of corresponding output files
!                     = 'W' => written to a newly created file
!*%enddate*         CHAR    End date/time
!*%refdate*         CHAR    Reference date/time
!*%startdate*       CHAR    Start date/time
!*%ncout*           INTEGER X-dimension of (regular) output grid
!*%nodim*           INTEGER (spatial) dimension of output grid (0,2,3)
!*%nostats*         INTEGER Number of output stations (irregular output)
!*%nowetout*        INTEGER Number of wet output data in case of land mask
!*%nrout*           INTEGER Y-dimension of (regular) output grid
!*%nstepout*        INTEGER Number of output time steps
!*%nzout*           INTEGER Z-dimension of output grid
!*%time_format*     INTEGER Format of time coordinate
!                     = 0 => date/time in string format
!                     = 1 => seconds
!                     = 2 => minutes
!                     = 3 => hours
!                     = 4 => days
!                     = 5 => weeks
!                     = 6 => months
!                     = 7 => years
!                     = 8 => date in years
!*%vcoord*          INTEGER Selects type of vertical coordinate variable
!                     = 1 => z-coordinate
!                     = 2 => sigma coordinate
!*%tlims*           INTEGER Start/End/Stride time index
!*%xlims*           INTEGER Start/End/Stride space indices in X-direction
!                           (regular output)
!*%ylims*           INTEGER Start/End/Stride space indices in Y-direction
!                           (regular output)
!*%zlims*           INTEGER Start/End/Stride space indices in Z-direction
!*%deltout*         REAL    Time step (unit determined by time_format)
!
!*Statlocs*         DERIVED Station locations for irregular user-defined output
!*%ipos*            INTEGER X-indices of output stations
!*%jpos*            INTEGER Y-indices of output stations
!*%name*            CHAR    Station name
!
!*VariableAtts*     DERIVED Variable attributes
!*%fill*            LOGICAL If .TRUE., a fill value is defined by the
!                           'fill_value' attribute. No fill value is defined if
!                           'fill' is .FALSE. or the 'fill' attribute of the
!                           asssociated file is .FALSE. 
!*%axis*            CHAR    CF axis attribute
!*%f90_name*        CHAR    FORTRAN name
!*%standard_name*   CHAR    netCDF CF standard name
!*%comment*         CHAR    CF comment attribute
!*%coordinates%     CHAR    CF coordinates attribute
!*%long_name*       CHAR    Long descriptive name
!*%vector_name*     CHAR    Associated vector name
!*%units*           CHAR    Variable unit
!*%node*            CHAR    Variable node on output grid
!*%data_type*       INTEGER FORTRAN type of variable
!*%ivarid*          INTEGER Variable key id
!*%klev*            INTEGER Vertical output grid level in case oopt = oopt_klev
!*%nrank*           INTEGER Variable rank
!*%numvar*          INTEGER Variable number (in case the last dimension of an
!                           output array represents a variable dimension)
!*%oopt*            INTEGER Applied operator id in case of a user output
!                           variable
!                     = oopt_null: no operator applied
!                     = oopt_min:  minimum value
!                     = oopt_max:  maximum value
!                     = oopt_mean: spatially mean value
!                     = oopt_int: spatially integrated value
!                     = oopt_klev: at vertical (C-node) grid level 'klev'
!                     = oopt_dep:  at vertical depth 'dep'

!*%global_dims*     INTEGER Global array shape (excluding halos)
!*%halo_dims*       INTEGER Halo sizes of array (WESN)
!*%local_dims*      INTEGER Local array shape (excluding halos)
!*%dep*             REAL    Vertical depth for output data in case
!                           oopt = oopt_dep
!*%fill_value*      REAL    Fill value for output data
!
!*VRelativeCoords*  DERIVED Vertical coordinates with respect to reference grid
!*%kcoord*          INTEGER Z-index coordinate with respect to bottom lower
!                           left corner of reference grid
!*%weights*         REAL    Weight factors for vertical interpolation
!
!************************************************************************
!

END MODULE datatypes
