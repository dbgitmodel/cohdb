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
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - test case rhone
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.1
!
! Description - test case rhone
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

!---names of log files
runlog_file = runtitle(1:5)//'.runlog'

!
!6. Timing
!---------
!

levtimer = 3

!
!7. Parallel setup
!-----------------
!

IF (npworld.GT.1) nprocscoh = 8


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
! Version - @(COHERENS)Usrdef_Model.f90  V2.5.1
!
! Description - test case rhone
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
USE switches
USE syspars
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_mod_params'
CALL log_timer_in()

!
!2. Switches
!-----------
!
!---Cartesian or spherical (0/1)
iopt_grid_sph = 1
!---vertical regular grid (0/1)
iopt_grid_vtype = 2
!---vertical grid transform(0/11/12/13/21)
iopt_grid_vtype_transf = 12

!---equation of state (0/1/2)
iopt_dens = 2
!---formulation for baroclinic pressure gradient (0/1/2)
iopt_dens_grad = 3
!---salinity equation (0/1)
iopt_sal = 2

!---advection scheme for 2-D currents (0/1/2/3/4)
iopt_adv_2D = 3
!---advection scheme for 3-D currents (0/1/2/3/4)
iopt_adv_3D = 3
!---advection scheme for scalars (0/1/2/3/4)
iopt_adv_scal = 3
!---type of vertical diffusion scheme (0/1/2/3)
IF (runtitle(6:6).EQ.'C') THEN
   iopt_vdif_coef = 2
ENDIF

!---implicit scheme
iopt_hydro_impl = 1
IF (iopt_hydro_impl.EQ.1) iopt_mg_smoother = 2

!---type of meteo input for surface stress and pressure (0/1/2)
IF (runtitle(6:6).GT.'C') THEN
   iopt_meteo = 1
   iopt_meteo_stres = 1
ENDIF

!---salinity o.b.c (0/1)
iopt_obc_sal = 1
!---2-D mode o.b.c (0/1)
iopt_obc_2D = 1
!---3-D mode o.b.c (0/1)
iopt_obc_3D = 1

!
!3. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
CStartDateTime(1:19) = '2003/01/01;00:00:00'
CEndDateTime(1:19) = '2003/01/05;00:00:00'

!---time step
delt2d = MERGE(10.0,180.0,iopt_hydro_impl.EQ.0)

!---counter for 3-D mode
ic3d = MERGE(18,1,iopt_hydro_impl.EQ.0)

!
!4. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 109; nr = 51; nz = 30

!---number of open sea boundaries
nosbu = 63; nosbv = 108

!---number of river open boundaries at V-nodes
nrvbv = 1

!---restart times
norestarts = 0

!---acceleration of gravity [m/s^2]
gacc_ref = 9.81

!---reference salinity [PSU]
sal_ref = 38.1

!---model grid
sigstar_DJ = 0.5
sig0_DJ = 0.1

!---depth flag
depmean_flag = -99.9

!---bottom drag coefficient
zrough_cst = 0.006

!---implicit settings
IF (iopt_hydro_impl.EQ.1) nomglevels = 2

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
modfiles(io_modgrd,1,1)%filename = 'rhonegrid.dat'
!---open boundary conditions (2-D)
modfiles(io_2uvobc,1:2,1)%status = 'N'
!---open boundary conditions (3-D)
modfiles(io_3uvobc,1:2,1)%status = 'N'
!---open boundary conditions (salinity)
modfiles(io_salobc,1:2,1)%status = 'N'
!---meteo
IF (iopt_meteo.EQ.1) THEN
   modfiles(io_metsur,1,1)%status = 'N'
   modfiles(io_metsur,1,1)%time_regular = .TRUE.
   IF (iopt_hydro_impl.EQ.0) THEN
      modfiles(io_metsur,1,1)%tlims = (/0,int_fill,-17280/)
   ELSE
      modfiles(io_metsur,1,1)%tlims = (/0,int_fill,-960/)
   ENDIF
ENDIF

!
!6.2 Output
!----------
!

IF (ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,2)%status = 'W'
   modfiles(io_modgrd,1,2)%form = 'U'
!  ---open boundary conditions (2-D)
   modfiles(io_2uvobc,1:2,2)%status = 'W'
   modfiles(io_2uvobc,1:2,2)%form = 'U'
!  ---open boundary conditions (3-D)
   modfiles(io_3uvobc,1:2,2)%status = 'W'
   modfiles(io_3uvobc,1:2,2)%form = 'A'
!  ---open boundary conditions (salinity)
   modfiles(io_salobc,1:2,2)%status = 'W'
   modfiles(io_salobc,1:2,2)%form = 'A'
ENDIF

!
!7. Surface grid parameters
!--------------------------
!

surfacegrids(igrd_model,1)%x0dat = 4.116667
surfacegrids(igrd_model,1)%y0dat = 43.08333
surfacegrids(igrd_model,1)%delxdat = 0.125E-01
surfacegrids(igrd_model,1)%delydat = 0.833334E-02

!
!8. User output
!--------------
!

nosetstsr = MERGE(2,1,iopt_verif.EQ.0)
novarstsr = 18

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.5.1
!
! Description - test case rhone
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
INTEGER :: i, ii, iunit, j, jj


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!3. Water depths [m]
!-------------------
!

CALL open_filepars(modfiles(io_modgrd,1,1))
iunit = modfiles(io_modgrd,1,1)%iunit

j_310: DO j=1,nr-1
   READ (iunit,*) depmeanglb(1:nc-1,j)
ENDDO j_310

CALL close_filepars(modfiles(io_modgrd,1,1))

!
!4. Open boundary locations
!--------------------------
!
!---West
ii = 0
j_441: DO j=1,nr-1
   IF (depmeanglb(1,j).GT.0.0) THEN
      ii = ii + 1
      iobu(ii) = 1
      jobu(ii) = j
   ENDIF
ENDDO j_441

!---East
j_442: DO j=1,nr-1
   IF (depmeanglb(nc-1,j).GT.0.0) THEN
      ii = ii + 1
      iobu(ii) = nc
      jobu(ii) = j
   ENDIF
ENDDO j_442

!---South
jj = 0
i_443: DO i=1,nc-1
   IF (depmeanglb(i,1).GT.0.0) THEN
      jj = jj + 1
      iobv(jj) = i
      jobv(jj) = 1
   ENDIF
ENDDO i_443

!---river mouth
iobv(nobv) = 59
jobv(nobv) = 37

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - test case rhone
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

!---type of conditions at U-open boundaries
ityp2dobu = 6; iloczobu = 1

!---type of conditions at V-open boundaries
ityp2dobv(1:nobv-1) = 6; ityp2dobv(nobv) = 4
iloczobv = 1

!---type and form of output data
no2dobuv = 1
index2dobuv(1,2) = nobu + nobv
iobc2dtype = 3

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
! Description - test case rhone
!
! Reference -
!
! Calling program - define_profobc_spec
!
!************************************************************************
!
USE gridpars
USE iopars
USE relaxation
USE time_routines, ONLY: log_timer_in, log_timer_out

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


procname(pglev+1) = 'usrdef_profobc_spec'
CALL log_timer_in()

SELECT CASE (iddesc)
   CASE (io_3uvobc,io_salobc)
      noprofsd = 1
      iprofobvy(nobv,1) = 1
      indexprof(1,2) = 1
END SELECT

CALL log_timer_out()


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
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - test case rhone
!
! Reference -
!
! Calling program - define_2dobc_data
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
!*Local variables
!
REAL :: qdis, width


procname(pglev+1) = 'usrdef_2dobc_data'
CALL log_timer_in()

IF (modfiles(iddesc,2,1)%iostat.EQ.0) THEN
   modfiles(iddesc,2,1)%iostat = 1
   GOTO 1000
ENDIF

ciodatetime = CStartDateTime

IF (runtitle(6:6).NE.'B') THEN
   qdis = 1500.0
ELSE
   qdis = 6000.0
ENDIF
width = 1010.225
data2d(1,1) = qdis/width

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
! Description - test case rhone
!
! Reference -
!
! Calling program - define_profobc_data
!
!************************************************************************
!
USE gridpars
USE iopars
USE syspars
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

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
!*numvars*     INTEGER  Number of variables in data file
!*nzdat*       INTEGER  Number of vertical levels in profiles (1 or nz)
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'usrdef_profobc_data'
CALL log_timer_in()

IF (modfiles(iddesc,2,1)%iostat.EQ.0) THEN
   modfiles(iddesc,2,1)%iostat = 1
   GOTO 1000
ENDIF

ciodatetime = CStartDateTime

SELECT CASE (iddesc)
CASE (io_3uvobc)
   psiprofdat(1,:) = 0.0
CASE (io_salobc)
   psiprofdat(1,:) = 10.0
END SELECT

1000 CALL log_timer_out()


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


RETURN

END SUBROUTINE usrdef_rlxobc_spec
