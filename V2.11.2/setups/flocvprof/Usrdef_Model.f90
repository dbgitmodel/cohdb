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
! *Usrdef_Model* User-defined model setup
!
! Author - IMDC
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.1
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description - test case flocvprof
!
! Reference -
!
! Routines - usrdef_init_params, usrdef_mod_params, usrdef_grid,
!            usrdef_partition, usrdef_phsics, usrdef_1dobc_spec,
!            usrdef_2dobc_spec, usrdef_physobc_spec, usrdef_1dobc_data,
!            usrdef_2dobc_data, usrdef_physobc_data, usrdef_rlxobc_spec
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_init_params
!************************************************************************
!
! *usrdef_init_params* Define parameters for monitoring
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.1
!
! Description - test case flocvprof
!
! Reference -
!
! Calling program - simulation_start
!
!************************************************************************
!
USE iopars
USE paralpars

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iproc


!
!1. Cold/warm start
!------------------
!

IF (ciffile%status.EQ.'W') cold_start = .TRUE.

!
!3. Log files
!------------
!
!---program leveling in log files
iproc_310: DO iproc=1,npworld
   levprocs_ini(iproc) = 7
   levprocs_run(iproc) = 7
ENDDO iproc_310
levprocs_err = 2

!
!6. Timing
!---------
!
!---timer level (0/1/2/3)
levtimer = 3


RETURN

END SUBROUTINE usrdef_init_params

!========================================================================

SUBROUTINE usrdef_mod_params
!************************************************************************
!
! *usrdef_mod_params*  Define parameters for physical model
!
! Author - IMDC
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.1
!
! Description - test case flocvprof
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE gridpars
USE iopars
USE physpars
USE switches
USE tide
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: sflag


procname(pglev+1) = 'usrdef_mod_params'
CALL log_timer_in()

!
!1. Setup flag
!-------------
!

sflag = MERGE(.TRUE.,.FALSE.,runtitle(10:10).EQ.'0')

!
!2. Switches
!-----------
!
!---model grid
iopt_grid_nodim = 1
!---vertical grid type (1/2/3)
!iopt_grid_vtype = 2

!---sediments
iopt_sed = MERGE(0,2,sflag)

!---bottom boundary conditions
iopt_bstres_drag = 1

!---turbulence
iopt_turb_ntrans = 2

!---1-D applications
iopt_sur_1D = 1

!---user output (0/1)
iopt_out_tsers = MERGE(0,1,sflag)

!
!3. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
IF (sflag) THEN
   CStartDateTime(1:23) = '2009/04/01;00:00:00:000'
   CEndDateTime(1:23)   = '2009/04/02;00:00:00:000'
ELSE
   CStartDateTime(1:23) = '2009/04/02;00:00:00:000'
   CEndDateTime(1:23)   = '2009/04/04;00:00:00:000'
ENDIF

!---time step (seconds)
delt2d = 12.0

!
!4. Physical model constants
!---------------------------
!
!---model grid
nz = 20

!---number of tidal constituents
nconobc = 1

!---acceleration of gravity [m/s^2]
gacc_ref = 9.81

!---uniform mean water depth
depmean_cst = 4.25

!---uniform bottom roughness length
bdragcoef_cst = 0.00464

!---tidal indices
index_obc(1) = icon_M2

!
!6. Model I/O file properties
!----------------------------
!
!6.1 Input
!---------
!
!---model grid
modfiles(io_modgrd,1,1)%status='0'

!---initial conditions
IF (sflag) THEN
   modfiles(io_inicon,ics_phys,1)%status = '0'
   modfiles(io_inicon,ics_sed,1)%status = '0'
ELSE
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = 'N'
   modfiles(io_inicon,ics_phys,1)%filename = &
                               & runtitle(1:9)//runtitle(11:11)//'.phsfin.nc'
   modfiles(io_inicon,ics_sed,1)%status = 'N'
ENDIF

!---surface forcing
modfiles(io_1uvsur,1:2,1)%status = 'N'

!---sediment particle attributes
modfiles(io_sedspc,1,1)%status = MERGE('0','N',sflag)

!
!6.2 Output
!----------
!

norestarts = 0

IF (ciffile%status.EQ.'W') THEN
!  ---restart times
   norestarts = 1
   ntrestart(1) = 0
!  ---model grid
   modfiles(io_modgrd,1,2)%status = 'W'
   modfiles(io_modgrd,1,2)%form = 'N'
!  ---surface forcing
   modfiles(io_1uvsur,1:2,2)%status = 'W'
   modfiles(io_1uvsur,1:2,2)%form = 'N'
!  ---sediment particle attributes
   modfiles(io_sedspc,1,2)%status = 'W'
   modfiles(io_sedspc,1,2)%form = 'N'
ENDIF

!---initial conditions
IF (sflag) THEN
   norestarts = 1
   ntrestart(1) = int_fill
   modfiles(io_fincon,ics_phys,2)%status = 'W'
   modfiles(io_fincon,ics_phys,2)%form = 'N'
   modfiles(io_fincon,ics_phys,2)%filename = &
                               & runtitle(1:9)//runtitle(11:11)//'.phsfin.nc'
ELSE
   norestarts = 0
   modfiles(io_fincon,ics_phys,2)%status = '0'
ENDIF

!
!9. User output
!--------------
!
!---title for user-output files
outtitle = runtitle(1:9)//runtitle(11:11)

!----time series output
nosetstsr = 1
novarstsr = 8

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_mod_params

!========================================================================

SUBROUTINE usrdef_grid
!************************************************************************
!
! *usrdef_grid* Define model grid arrays
!
! Author - IMDC
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - test case flocvprof
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_grid

!========================================================================

SUBROUTINE usrdef_partition
!************************************************************************
!
! *usrdef_partition* Define domain decomposition
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - domain_decomposition
!
!************************************************************************
!

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_partition

!========================================================================

SUBROUTINE usrdef_phsics
!************************************************************************
!
! *usrdef_phsics* Define initial conditions for physics
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_phsics
!
!************************************************************************
!

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_phsics

!========================================================================

SUBROUTINE usrdef_1dsur_spec
!************************************************************************
!
! *usrdef_1dsur_spec* Define specifier arrays for 1-D forcing
!
! Author - IMDC
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.1
!
! Description - test case flocvprof
!
! Reference -
!
! Calling program - define_1dobc_spec
!
!************************************************************************
!
USE obconds
USE syspars  

IMPLICIT NONE


!---amplitudes and phases of X-component of pressure gradient
gxslope_amp(1) = 3.0E-04
gxslope_pha(1) = 0.0

!---amplitudes and phases of Y-component of pressure gradient
gyslope_amp(1) = 0.0
gyslope_pha(1) = 0.0

!---type of data 
isur1dtype = 3


RETURN

END SUBROUTINE usrdef_1dsur_spec

!========================================================================

SUBROUTINE usrdef_2dobc_spec(iddesc)
!************************************************************************
!
! *usrdef_2dobc_spec* Define specifier arrays for 2-D open boundary conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description -
!
! Reference -
!
! Calling program - define_2dobc_spec
!
!************************************************************************
!

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id (io_2uvobc or io_2xyobc)
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_2dobc_spec

!========================================================================

SUBROUTINE usrdef_profobc_spec(iddesc,itypobux,itypobvy,iprofobux,iprofobvy,&
                             & iprofrlx,noprofsd,indexprof,indexvar,novars,&
                             & nofiles,nobux,nobvy)
!************************************************************************
!
! *usrdef_profobc_spec* Define specifier arrays for open boundary conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description -
!
! Reference -
!
! Calling program - define_profobc_spec
!
!************************************************************************
!
USE relaxation

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, nobux, nobvy, nofiles, novars
INTEGER, INTENT(INOUT), DIMENSION(2:nofiles) :: noprofsd
INTEGER, INTENT(OUT), DIMENSION(nobux) :: itypobux
INTEGER, INTENT(OUT), DIMENSION(nobvy) :: itypobvy
INTEGER, INTENT(INOUT), DIMENSION(nobux,novars) :: iprofobux
INTEGER, INTENT(INOUT), DIMENSION(nobvy,novars) :: iprofobvy
INTEGER, INTENT(INOUT), DIMENSION(nobux+nobvy,2:nofiles) :: indexprof
INTEGER, INTENT(INOUT), DIMENSION(nobux+nobvy,2:nofiles) :: indexvar
INTEGER, INTENT(INOUT), DIMENSION(norlxzones) :: iprofrlx

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*     INTEGER Data file id
!*itypobux*   INTEGER Type of U- or X-open boundary condition
!*itypobvy*   INTEGER Type of V- or Y-open boundary condition
!*iprofobux*  INTEGER Profile numbers at U or X--open boundaries
!*iprofobvy*  INTEGER Profile numbers at V- or Y-open boundaries
!*iprofrlx*   INTEGER Disables/enables relaxation at open boundary zones
!*noprofsd*   INTEGER Number of profiles per data file
!*indexprof*  INTEGER Mapping array of the profile numbers in the data files to
!                     the profile numbers assigned to the open boundaries. The
!                     physical size of the first dimension equals the number of
!                     profiles in a data file.
!*indexvar*   INTEGER Defines the variable number of the profiles in a data
!                     file. The physical size of the first dimension equals the
!                     number of profiles in a data file.
!*novars*     INTEGER Total number of variables
!*nofiles*    INTEGER Number of data files
!*nobux*      INTEGER Number of nodes at U- or X-open boundaries
!*nobvy*      INTEGER Number of nodes at V- or Y-open boundaries
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_profobc_spec

!========================================================================

SUBROUTINE usrdef_1dsur_data(ciodatetime,data1d,novars)
!************************************************************************
!
! *usrdef_1dsur_data* Define data for 1-D forcing
!
! Author - IMDC
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.1
!
! Description -
!
! Reference -
!
! Calling program - define_1dsur_data
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
INTEGER, INTENT(IN) :: novars
REAL, INTENT(INOUT), DIMENSION(novars) :: data1d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*ciodatetime* CHAR    Date/time in data file
!*data1d*      REAL    Surface data
!*novars*      INTEGER Number of surface data
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'usrdef_1dsur_data'
CALL log_timer_in()

IF (modfiles(io_1uvsur,2,1)%iostat.EQ.0) THEN
   modfiles(io_1uvsur,2,1)%iostat = 1
   GOTO 1000
ENDIF

ciodatetime = CStartDateTime

data1d(1) = -4.0E-04
data1d(2) = 0.0

1000 CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_1dsur_data

!========================================================================

SUBROUTINE usrdef_2dobc_data(iddesc,ifil,ciodatetime,data2d,nodat,novars)
!************************************************************************
!
! *usrdef_2dobc_data* Define open boundary data for 2-D mode
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description -
!
! Reference -
!
! Calling program - define_2dobc_data
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, nodat, novars
REAL, INTENT(INOUT), DIMENSION(nodat,novars) :: data2d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER Data file id (io_2uvobc or io_2xyobc)
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*data2d*      REAL    Input data
!*nodat*       INTEGER Number of input data
!*novars*      INTEGER Number of input parameters
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_2dobc_data

!========================================================================

SUBROUTINE usrdef_profobc_data(iddesc,ifil,ciodatetime,psiprofdat,numprofs)
!************************************************************************
!
! *usrdef_profobc_data* Define physical open boundary profiles
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.2
!
! Description -
!
! Reference -
!
! Calling program - define_profobc_data
!
!************************************************************************
!
USE gridpars
USE syspars

IMPLICIT NONE

!
!*  Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, numprofs
REAL, INTENT(INOUT), DIMENSION(numprofs,nz) :: psiprofdat

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER  No. of data file
!*ciodatetime* CHAR     Date/time in data file
!*psiprofdat*  REAL     Profile arrays
!*numprofs*    INTEGER  Number of profiles in data file
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_profobc_data

!========================================================================

SUBROUTINE usrdef_rlxobc_spec
!************************************************************************
!
! *usrdef_rlxobc_spec* Define specifier arrays for relaxation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_rlxobc_spec
!
!************************************************************************
!

IMPLICIT NONE


RETURN


END SUBROUTINE usrdef_rlxobc_spec
