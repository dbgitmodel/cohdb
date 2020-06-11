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
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.8
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - test case drythin
!
! Reference -
!
! Subroutines - usrdef_init_params, usrdef_mod_params, usrdef_grid,
!               usrdef_partition, usrdef_phsics, usrdef_1dsur_spec,
!               usrdef_2dobc_spec, usrdef_profobc_spec, usrdef_1dsur_data,
!               usrdef_2dobc_data, usrdef_profobc_data, usrdef_rlxobc_spec
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_init_params
!************************************************************************
!
! *usrdef_init_params* Define parameters for monitoring
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.6
!
! Description - test case drythin
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
!1. Cold/warm start
!------------------
!

IF (ciffile%status.EQ.'W') cold_start = .TRUE.

!
!3. Log files
!------------
!
!---program leveling in log files
levprocs_ini = 3
levprocs_run = 3

!
!6. Timing
!---------
!

levtimer = 3

!
!7. Parallel setup
!-----------------
!

IF (npworld.GT.1) nprocscoh = 4


RETURN

END SUBROUTINE usrdef_init_params

!========================================================================

SUBROUTINE usrdef_mod_params
!************************************************************************
!
! *usrdef_mod_params* Define parameters for physical model
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.8
!
! Description - test case drythin
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE gridpars
USE iopars
USE paralpars
USE physpars
USE structures
USE switches
USE syspars
USE tide
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: idesc, ifil, iotype


procname(pglev+1) = 'usrdef_mod_params'
CALL log_timer_in()

!
!2. Switches
!-----------
!
!---2-D mode o.b.c (0/1)
iopt_obc_2D = 1

!---structures module (1/0)
SELECT CASE (runtitle(8:8))
   CASE ('A')
      iopt_drycel = 1
      iopt_thndam = 0
   CASE ('B')
      iopt_drycel = 0
      iopt_thndam = 1
END SELECT

!
!4. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
CStartDateTime(1:23) = '2003/01/01;00:00:00:000'
CEndDateTime(1:23) = '2003/01/03;00:00:00:000'

!---time step (milliseconds)
delt2d = 5.0
ic3d = 6

!
!5. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 101; nr = 51; nz = 10

!---number of structures/open sea boundaries
nosbu = 100

!---number of tidal constituents
nconobc = 1

!---acceleration of gravity [m/s^2]
gacc_ref = 9.81

!---uniform mean water depth
depmean_cst = 10.0

!---mean depth flag
depmean_flag = -99.9 

!---uniform bottom roughness length
zrough_cst = 0.006

!---tidal indices
index_obc(1) = icon_S2

!---number of structures
IF (iopt_drycel.EQ.1) numdry = 400
IF (iopt_thndam.EQ.1) THEN
   numthinu = 70; numthinv = 0
ENDIF

!
!6. Model I/O file properties
!----------------------------
!
!6.1 Input
!---------
!
!---model grid
modfiles(io_modgrd,1,1)%status = 'N'
!---open boundary conditions (2-D)
modfiles(io_2uvobc,1,1)%status = 'N'
!---dry cells
IF (iopt_drycel.EQ.1) modfiles(io_drycel,1,1)%status = 'N'
!---thin dams
IF (iopt_thndam.EQ.1) modfiles(io_thndam,1,1)%status = 'N'

!
!6.2 Output
!----------
!

IF (ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,2)%status = 'W'
   modfiles(io_modgrd,1,2)%form = 'U'
!  ---open boundary conditions (2-D)
   modfiles(io_2uvobc,1,2)%status = 'W'
   modfiles(io_2uvobc,1,2)%form = 'U'
!  ---dry cells
   IF (iopt_drycel.EQ.1) THEN
      modfiles(io_drycel,1,2)%status = 'W'
      modfiles(io_drycel,1,2)%form = 'A'
   ENDIF
!  ---thin dams
   IF (iopt_thndam.EQ.1) THEN
      modfiles(io_thndam,1,2)%status = 'W'
      modfiles(io_thndam,1,2)%form = 'A'
   ENDIF
ENDIF

!
!7. Surface grid parameters
!--------------------------
!

surfacegrids(igrd_model,1)%delxdat = 100.0
surfacegrids(igrd_model,1)%delydat = 100.0

!
!8. User output
!--------------
!

nosetstsr = 1
novarstsr = 9

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_mod_params

!========================================================================

SUBROUTINE usrdef_grid
!************************************************************************
!
! *usrdef_grid* Define model grid arrays
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.6
!
! Description - test case drythin
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!*Local variables
!
INTEGER :: j

procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!1. Open boundary locations
!--------------------------
!

iobu(1:50) = 1
iobu(51:100) = nc

jobu(1:50) = (/(j,j=1,50)/)
jobu(51:100) = (/(j,j=1,50)/) 

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_grid

!========================================================================

SUBROUTINE usrdef_partition
!************************************************************************
!
! *usrdef_partition* Define domain decomposition
!
! Author - ANTEA
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
! *usrdef_phsics* Define arrays initial conditions for physics
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - read_phsics
!
! Module calls - exchange_mod
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
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_1dsur_spec
!
!************************************************************************
!

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_1dsur_spec

!========================================================================

SUBROUTINE usrdef_2dobc_spec(iddesc)
!************************************************************************
!
! *usrdef_2dobc_spec* Define specifier arrays for 2-D open boundary conditions
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - test case drythin
!
! Reference -
!
! Calling program - define_2dobc_spec
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE obconds
USE physpars
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

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


procname(pglev+1) = 'usrdef_2dobc_spec'
CALL log_timer_in()

!---type of conditions at open boundaries
ityp2dobu(1:50) = 3 
ityp2dobu(51:100) = 13

!---location elevation point
iloczobu = 1

!---amplitudes and phases
zdatobu_amp(1:50,1) = 2.0
zdatobu_pha(1:50,1) = 0.0

CALL log_timer_out()


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
! Author - ANTEA
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
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_1dsur_data
!
!************************************************************************
!
USE syspars

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


RETURN

END SUBROUTINE usrdef_1dsur_data

!========================================================================

SUBROUTINE usrdef_2dobc_data(iddesc,ifil,ciodatetime,data2d,nodat,novars)
!************************************************************************
!
! *usrdef_2dobc_data* Define open boundary data for 2-D mode
!
! Author - ANTEA
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

SUBROUTINE usrdef_profobc_data(iddesc,ifil,ciodatetime,psiprofdat,nzdat,&
                             & numvars,numprofs)
!************************************************************************
!
! *usrdef_profobc_data* Define physical open boundary profiles
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.2
!
! Description - called by reader processes only
!
! Reference -
!
! Calling program - define_profobc_data
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

!
!*  Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, numprofs, numvars, nzdat
REAL, INTENT(INOUT), DIMENSION(nzdat,numvars,numprofs) :: psiprofdat

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER  No. of data file
!*ciodatetime* CHAR     Date/time in data file
!*psiprofdat*  REAL     Profile arrays
!*nzdat*       INTEGER  Number of vertical levels in profiles (1 or nz)
!*numvars*     INTEGER  Number of variables in data file
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
! Author - ANTEA
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


RETURN

END SUBROUTINE usrdef_rlxobc_spec
