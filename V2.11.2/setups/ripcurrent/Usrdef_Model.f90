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
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90 V2.9 
!
! $Date: 2018-05-09 09:43:10 +0200 (Wed, 09 May 2018) $
!
! $Revision: 1130 $
!
! Description - test case rip current
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
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.9
!
! Description - test case rip current
!
! Reference -
!
! Calling program - simulation_start
!
!************************************************************************
!
USE iopars
USE paralpars
USE syspars
USe switches

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
   levprocs_ini(iproc) = 3
   levprocs_run(iproc) = 3
ENDDO iproc_310

!
!6. Timing
!---------
!

levtimer = 3

!
!7. Parallel setup
!-----------------
!

IF (npworld.GT.1) THEN
   nprocscoh = 2
   nprocswav = 2
   iopt_waves_model = 1
ELSE
   nprocscoh = 1
   nprocswav = 0
   iopt_waves_model = 0
ENDIF


RETURN

END SUBROUTINE usrdef_init_params

!========================================================================

SUBROUTINE usrdef_mod_params
!************************************************************************
!
! *usrdef_mod_params* Define parameters for physical model
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.9
!
! Description - test case rip current
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE gridpars
USE iopars
USE nestgrids
USE paralpars
USE physpars
USE switches
USE syspars
USE tide
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
REAL :: delxdat, delydat


procname(pglev+1) = 'usrdef_mod_params'
CALL log_timer_in()

!
!2. Switches
!-----------
!
!---grid dimension (1/2/3)
iopt_grid_nodim = 2

!---dimension of bottom stress formulation (2/3)
iopt_bstres_nodim = 2

!---advection scheme for 2-D currents (0/1/2/3)
iopt_adv_2D = 1
!---advection scheme for 3-D currents (0/1/2/3)
iopt_adv_3D = 3

!---2-D mode o.b.c (0/1)
iopt_obc_2D = 1 
 
!---wave current interaction
iopt_waves = 1; iopt_waves_curr = 1
iopt_bstres_waves_bfric = 1
iopt_waves_couple = MERGE(1,3,runtitle(11:11).EQ.'E')
SELECT CASE (runtitle(11:11))
   CASE ('A','E')
      iopt_waves_form = 2
      iopt_waves_dissip = 1; iopt_waves_pres = 1
  CASE ('B')
      iopt_waves_form = 2
      iopt_waves_dissip = 0; iopt_waves_pres = 1
   CASE ('C')
      iopt_waves_form = 2
      iopt_waves_dissip = 1; iopt_waves_pres = 0
   CASE ('D')
      iopt_waves_form = 1
      iopt_waves_dissip = 1; iopt_waves_pres = 1
END SELECT
IF (npworld.EQ.1) iopt_waves_couple = 0

!
!3. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
CStartDateTime(1:19) = '2004/12/24;00:00:00'
CEndDateTime(1:19)   = '2004/12/24;00:10:00'

!---time step (seconds)
delt2d = 0.04

!
!4. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 80; nr = 94

!---number of open sea boundaries
nosbv = 2*(nc-1); nosbu = nr-1

!---uniform bottom roughness length
zrough_cst = 0.0015

!
!6. Model I/O file properties
!----------------------------
!
!6.1 Input
!---------
!
!---model grid
modfiles(io_modgrd,1,1)%status = 'N'
modfiles(io_modgrd,1,1)%form = 'A'
modfiles(io_modgrd,1,1)%filename = 'PB_COHERENS_RC_org.dep'

!---open boundary conditions
modfiles(io_2uvobc,1:2,1)%status = 'N'

!
!6.2 Output
!----------
!
!---wave data
IF (iopt_waves_couple.EQ.0) THEN
   modfiles(io_wavsur,1,1)%status = 'R'
   modfiles(io_wavsur,1,1)%filename = TRIM(runtitle)//'.wavsur.nc'
   modfiles(io_wavsur,1,1)%form = 'N'
   modfiles(io_wavsur,1,1)%tlims = (/0,int_fill,1500/)
ELSE
   modfiles(io_wavsur,1,2)%status = 'W'
   modfiles(io_wavsur,1,2)%form = 'N'
ENDIF

IF (ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,2)%status = 'W'
   modfiles(io_modgrd,1,2)%form = 'N'
!  ---open boundary conditions
   modfiles(io_2uvobc,1:2,2)%status = 'W'
   modfiles(io_2uvobc,1:2,2)%form = 'A'
ENDIF

!
!7. Surface grid parameters
!--------------------------
!
!---model
delxdat = 0.2; delydat = 0.2
surfacegrids(igrd_model,1)%x0dat = 0.0
surfacegrids(igrd_model,1)%y0dat = 0.0
surfacegrids(igrd_model,1)%delxdat = delxdat
surfacegrids(igrd_model,1)%delydat = delydat

!---wave grid
IF (iopt_waves_couple.EQ.0) THEN
   surfacegrids(igrd_waves,1)%nhtype = 4
ENDIF

!
!8. User output
!--------------
!

nosetstsr = 2
SELECT CASE (runtitle(11:11))
   CASE ('A','D','E')
      novarstsr = 14
   CASE ('B')
      novarstsr = 12
   CASE ('C')
      novarstsr = 13
END SELECT
nostatstsr = 6

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_mod_params

!========================================================================

SUBROUTINE usrdef_grid
!************************************************************************
!
! *usrdef_grid* Define model grid arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.9
!
! Description - test case rip current
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - close_filepars, open_filepars
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE inout_routines, ONLY: close_filepars, open_filepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: ii, iunit, i,j, jj


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!1. Open grid file
!-----------------
!

CALL open_filepars(modfiles(io_modgrd,1,1))
iunit = modfiles(io_modgrd,1,1)%iunit

!
!2. Water depths
!---------------
!
j_210: DO j=1,nr-1
   READ (iunit,*) depmeanglb(1:nc-1,j)
ENDDO j_210

!
!3. Open boundary locations
!--------------------------
!

!---West
iobu(1:nr-1) = 1
jobu (1:nr-1) = (/(i,i=1,nr-1)/)

!---South
iobv(1:nc-1) = (/(i,i=1,nc-1)/)
jobv(1:nc-1) = nr

!---North
iobv(nc:2*nc-2) = (/(i,i=1,nc-1)/)
jobv(nc:2*nc-2) = 1

!---West
iobu(1:nr-1) = 1
jobu (1:nr-1) = (/(i,i=1,nr-1)/)

!
!4. Close data file
!------------------
!

CALL close_filepars(modfiles(io_modgrd,1,1))

CALL log_timer_out()


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
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.1
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
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
! Author - Patrick Luyten
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
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.9
!
! Description - test case rip current
!
! Reference -
!
! Calling program - define_2dobc_spec
!
!************************************************************************
!
USE gridpars
USE iopars
USE obconds
USE tide
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

!
!*Local variables
!
INTEGER :: l


procname(pglev+1) = 'usrdef_2dobc_spec'
CALL log_timer_in()

!
!1. Type of conditions
!---------------------
!

ityp2dobu(1:nobu) = 12
ityp2dobv(1:nc-1) = 5
ityp2dobv(nc:2*(nc-2)) = 5

!
!2. Data file specifiers
!-----------------------
!

no2dobuv = nobu
iobc2dtype = 2

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
INTEGER, INTENT(INOUT), DIMENSION(nobux) :: itypobux
INTEGER, INTENT(INOUT), DIMENSION(nobvy) :: itypobvy
INTEGER, INTENT(INOUT), DIMENSION(nobux,novars) :: iprofobux
INTEGER, INTENT(INOUT), DIMENSION(nobvy,novars) :: iprofobvy
INTEGER, INTENT(INOUT), DIMENSION(novars*(nobux+nobvy),2:nofiles) :: indexprof
INTEGER, INTENT(INOUT), DIMENSION(novars*(nobux+nobvy),2:nofiles) :: indexvar
INTEGER, INTENT(INOUT), DIMENSION(norlxzones) :: iprofrlx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id
!*itypobux*  INTEGER Type of U- or X-open boundary condition
!*itypobvy*  INTEGER Type of V- or Y-open boundary condition
!*iprofobux* INTEGER Profile numbers at U- or X-open boundaries
!*iprofobvy* INTEGER Profile numbers at V- or Y-open boundaries
!*iprofrlx*  INTEGER Disables/enables relaxation at open boundary zones
!*noprofsd*  INTEGER Number of profiles per data file
!*indexprof* INTEGER Mapping array of the profile numbers in the data files to
!                    the profile numbers assigned to the open boundaries. The
!                    physical size of the first dimension equals the number of
!                    profiles in a data file.
!*indexvar*  INTEGER Defines the variable number of the profiles in a data
!                    file. The physical size of the first dimension equals the
!                    number of profiles in a data file.
!*novars*    INTEGER Total number of variables
!*nofiles*   INTEGER Number of data files (+1)
!*nobux*     INTEGER Number of nodes at U- or X-open boundaries
!*nobvy*     INTEGER Number of nodes at V- or Y-open boundaries
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
! Author - Patrick Luyten
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
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.9
!
! Description - test case rip current
!
! Reference -
!
! Calling program - define_2dobc_data
!
!************************************************************************
!
USE iopars
USE gridpars
USE syspars
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out
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
INTEGER :: ii


CALL log_timer_in()

IF (modfiles(io_2uvobc,ifil,1)%iostat.EQ.0) THEN
   modfiles(io_2uvobc,ifil,1)%iostat = 1
   GOTO 1000
ENDIF

ciodatetime = CStartDateTime
data2d = 0.0

1000 CALL log_timer_out()


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
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.2
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
