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

MODULE syspars
!************************************************************************
!
! *syspars* System parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)syspars.f90  V2.11.1
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
!************************************************************************
!
IMPLICIT NONE

!--kind parameters
INTEGER, PARAMETER :: kndchar = KIND('A'), kndlog = KIND(.TRUE.), &
                    & kndint = 4, kndilong = 8, kndreal = 4, kndrlong = 8, &
                    & kndcmplx = 4, kndclong = 8

!---kind type parameter of real data in the simulation
INTEGER :: kndrtype

!---data types
INTEGER, PARAMETER :: char_type = 1, log_type = 2, int_type = 3, &
                    & ilong_type = 4, real_type = 5, rlong_type = 6, &
                    & complx_type = 7, clong_type = 8
INTEGER :: float_type

!---universal parameters
REAL :: pi, halfpi, twopi, enap, degtorad, radtodeg
REAL, PARAMETER :: pi_s = 3.14159265, halfpi_s = 1.57079633, &
                 & twopi_s = 6.28318531, enap_s = 2.718282
REAL, PARAMETER :: degtorad_s = pi_s/180.0, radtodeg_s = 180.0/pi_s
REAL (KIND=kndrlong), PARAMETER :: pi_d = 3.1415926535897932_kndrlong
REAL (KIND=kndrlong), PARAMETER :: halfpi_d = 1.5707963267948966_kndrlong
REAL (KIND=kndrlong), PARAMETER :: twopi_d = 6.2831853071795865_kndrlong
REAL (KIND=kndrlong), PARAMETER :: enap_d = 2.7182818284590452_kndrlong
REAL (KIND=kndrlong), PARAMETER :: degtorad_d = pi_d/180.0_kndrlong, &
                                 & radtodeg_d = 180.0_kndrlong/pi_d

!---tidal parameters
INTEGER, PARAMETER :: MaxAstroTides = 56, MaxConstituents = 77

!---random generators
INTEGER, PARAMETER :: MaxGenerators = 32

!---MPI communications
INTEGER, PARAMETER :: MaxHaloComms = 8

!---Thatcher-Harleman BC
INTEGER, PARAMETER :: MaxVLevels = 500

!---number of multigrids
INTEGER, PARAMETER :: MaxMGLevels = 5

!---model variables
INTEGER, PARAMETER :: &
    & MaxModArids = 562, MaxBioArids = 24, MaxSedArids = 196, &
    & MaxPartArids = 41, &
    & MaxTotArids = MaxModArids + MaxSedArids + MaxBioArids + MaxPartArids

!---model I/O
INTEGER, PARAMETER :: MaxCIFTypes = 12, MaxCIFVars = 200, MaxGlbAtts = 20, &
                    & MaxGridTypes = 5, MaxGridFiles = 2, MaxInitFiles = 5, &
                    & MaxIOFiles = 33, MaxIOTypes = 54, MaxProgLevels  = 20, &
                    & MaxRestarts = 10

!---monitoring files and error coding
INTEGER, PARAMETER :: MaxErrCodes = 15, MaxErrMesgs = 50, MaxTimers = 42

!---character string lengths
INTEGER, PARAMETER :: lencifblock = 60, lencifline = 2000, lencifvar = 120,  &
                    & lencomm = 120, lendesc = 120, lenerrcode = 120, &
                    & lenformat = 120, lenfreq = 7, leniofile = 120, &
                    & leninst = 120, lenname = 35, lennode = 3, lenform = 1, &
                    & lentime = 23, lentitle = 20, lenunit = 60, lenversion = 15

!---cif file
CHARACTER (LEN=1), PARAMETER :: cifcom = '!', cifend ='#', cifsep =','

!---user output
LOGICAL, PARAMETER :: DegreesOut = .TRUE.
CHARACTER (LEN=lenversion), PARAMETER :: CF_version = 'CF-1.6', &
                                       & model_version = 'V2.11.2'

!---output formats
CHARACTER (LEN=lenformat), PARAMETER :: IntegerFormat='(50I11)', &
                                      & LogicalFormat ='(50L1)', &
                                      & RealFormat = '(50G16.7)', &
                                      & DoubleFormat = '(50(E23.14E3))'

!---undefined and zero values
INTEGER, PARAMETER :: fill_int_def = -2147483647, izero = 0
INTEGER (KIND=kndilong), PARAMETER :: izero_d = 0_kndilong
REAL, PARAMETER :: fill_real_def = 9.9692099683868690E+36, rzero = 0.0
REAL (KIND=kndrlong), PARAMETER :: fill_double_def = 9.9692099683868690E+36, &
                                 & rzero_d = 0.0_kndrlong
CHARACTER (LEN=lentime), PARAMETER :: cdatetime_undef = &
                                    & 'xxxx/xx/xx;00:00:00:000'

!
! Name             Type    Purpose
!------------------------------------------------------------------------------
!*DegreesOut*      LOGICAL Determines output unit of phase angles
!                          .TRUE.  => degrees
!                          .FALSE. => radians
!*IntegerFormat*   CHAR    Format specification for formatted integer output
!*MaxMGLevsls*     INTEGER Maximum number of multigrid levels (including main
!                          grid)
!*MaxAstroTides*   INTEGER Maximum number of constituents for astronomical
!                          force
!*MaxBioArids*     INTEGER Maximum number of biological array ids
!*MaxCIFTypes*     INTEGER Maximum number of CIF files
!*MaxCIFVars*      INTEGER Maximum number of data variables on a CIF data line
!*MaxConstituents* INTEGER Maximum number of tidal constituents at open
!                          boundaries
!*MaxErrCodes*     INTEGER Maximum type error messages
!*MaxErrMesgs*     INTEGER Default maximum of error messages
!*MaxGlbAtts*      INTEGER Maximum number of (user defined) global attributes
!                          in user output files
!*MaxGenerators*   INTEGER Maximum number of random generators which can be
!                          defined within the program
!*MaxGridFiles*    INTEGER Maximum number of surface grid types
!*MaxGridTypes*    INTEGER Maximum number of surface grid files per type
!*MaxHaloComms*    INTEGER Maximum number of allowed halo communications
!                          (send or receive)
!*MaxModArids*     INTEGER Maximum number physical model array ids
!*MaxInitFiles*    INTEGER Maximum number of initial condition files
!*MaxIOFiles*      INTEGER Maximum number of I/O files per type
!*MaxIOTypes*      INTEGER Maximum number of model I/O types
!*MaxPartArids*    INTEGER Maximum number of particle key ids
!*MaxVLevels*      INTEGER Maximum (vertical) dimension of the 'return_time' 
!                          vector (used for Thatcher-Harleman open boundary
!                          condition)
!*MaxProgLevels*   INTEGER Maximum number of subprogram levels
!*MaxRestarts*     INTEGER Maximum number for writing restart conditions
!*MaxSedArids*     INTEGER Maximum number of sediment array ids
!*MaxTimers*       INTEGER Maximum number of timers
!*MaxTotArids*     INTEGER Maximum number all model array ids
!*RealFormat*      CHAR    Format specification for formatted real output
!*cdatetime_undef* CHAR    Flag for undefined date/times
!*cifcom*          CHAR    Comment character on a CIF line
!*cifend*          CHAR    Marks the end of a block in a CIF
!*cifsep*          CHAR    Data separator data within a CIF line
!*degtorad*        REAL/LONG   Factor to convert degrees to radians
!                              (default kind)
!*degtorad_d*      LONG        Factor to convert degrees to radians
!                              (double precision)
!*degtorad_s*      SINGLE REAL Factor to convert degrees to radians
!                              (single precision)
!*enap*            REAL/LONG   Euler's number e (default kind)
!*enap_d*          LONG        Euler's number e (double precision)
!*enap_s*          SINGLE REAL Euler's number e (single precision)
!*fill_double_def* DOUBLE  Default flag for undefined/invalid/missing double
!                          real data or model variables (except netCDF data)
!*fill_int_def*    INTEGER Default flag for undefined/invalid/missing integer
!                          data or model variables (except netCDF data)
!*fill_longint_def*LONGINT Default flag for undefined/invalid/missing long
!                          (64bits) integer data or model variables
!*fill_real_def*   REAL    Default flag for undefined/invalid/missing real data
!                          or model variables (except netCDF data)
!*halfpi*          REAL/LONG   Number pi divided by 2 (default kind)
!*halfpi_d*        LONG        Number pi divided by 2 (double precision)
!*halfpi_s*        SINGLE REAL Number pi divided by 2 (single precision)
!*izero*           INTEGER Zero (single precision)
!*izero_d*         LONGINT Zero (double precision)
!*kndchar*         INTEGER Kind parameter for character variables
!*kndcmplx*        INTEGER Kind parameter for complex variables
!*kndilong*        INTEGER Kind parameter for double precision integer variables
!*kndint*          INTEGER Kind parameter for single precision integer variables
!*kndlog*          INTEGER Kind parameter for logical variables
!*kndreal*         INTEGER Kind parameter for single precision real variables
!*kndrlong*        INTEGER Kind parameter for double precision real variables
!*kndrtype*        INTEGER Kind parameter for real variables in the simulation
!                          (=kndreal or kndrlong depending on compiler option)
!*lencifline*      INTEGER Maximum length of a data line in a CIF file
!*lencifvar*       INTEGER Maximum length of a CIF data variable in string
!                          format
!*lencomm*         INTEGER Length of the "comment" global attribute
!*lendesc*         INTEGER Length of "long_name" attribute
!*lenerrcode*      INTEGER Length of error code messages
!*lenformat*       INTEGER Length of a format specification
!*lenfreq*         INTEGER Length of the name of a frequency
!*leninst*         INTEGER Length of "Institute" global attribute
!*leniofile*       INTEGER Maximum length for I/O file names
!*lenname*         INTEGER Length of "f90_name" attribute
!*lennode*         INTEGER Length of a nodal type string
!*lenform*         INTEGER Length of the file format
!*lentime*         INTEGER Length of date/time string
!*lentitle*        INTEGER Length of simulation title
!*lenunit*         INTEGER Length of "units" attribute
!*lenversion*      INTEGER Length of model_version string
!*model_version*   CHAR    COHERENS version
!*pi*              REAL/LONG   Number pi (default kind)
!*pi_d*            LONG        Number pi (double precision)
!*pi_s*            SINGLE REAL Number pi (single precision)
!*radtodeg*        REAL/LONG   Factor to convert radians to degrees
!                              (default kind)
!*radtodeg_d*      LONG        Factor to convert radians to degrees
!                              (double precision)
!*radtodeg_s*      SINGLE REAL Factor to convert radians to degrees
!                              (single precision)
!*rzero*           REAL    Zero (single precision)
!*rzero_d*         DOUBLE  Zero (double precision)
!*twopi*           REAL/LONG   Number pi times 2 (default kind)
!*twopi_d*         LONG        Number pi times 2 (double precision)
!*twopi_s*         SINGLE REAL Number pi tilmes 2 (single precision)
!
!************************************************************************
! 

END MODULE syspars
