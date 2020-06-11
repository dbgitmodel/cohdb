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
! Version - @(COHERENS)Usrdef_Model.f90  V2.9
!
! $Date: 2018-06-05 16:03:05 +0200 (Tue, 05 Jun 2018) $
!
! $Revision: 1142 $
!
! Description - optos nos
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
! Description - optos nos
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

IF (npworld.GT.1) nprocscoh = 11


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
! Description - optos nos
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
LOGICAL :: test
INTEGER :: icount
!REAL :: delxdat, delydat


procname(pglev+1) = 'usrdef_mod_params'
CALL log_timer_in()

test = runtitle(7:7).EQ.'0'

!
!2. Switches
!-----------
!
!---Cartesian or spherical (0/1)
iopt_grid_sph = 1

!---implicit scheme
iopt_hydro_impl = 1

!---flooding/drying scheme
iopt_fld = 0

!---dimension of bottom stress formulation (2/3)
iopt_bstres_nodim = 3
iopt_bstres_drag = 3

!---advection scheme for 2-D currents (0/1/2/3)
iopt_adv_2D = 1
!---advection scheme for 3-D currents (0/1/2/3)
iopt_adv_3D = 3

!---2-D mode o.b.c (0/1)
iopt_obc_2D = 1
!---open boundary condition for advective fluxes
iopt_obc_advflux = 2
!---relaxation open boundary condition for momentum (0/1)
iopt_obc_advrlx = 1
!---nodal corrections for tides (0/1)
iopt_astro_pars = 2

!---nesting
iopt_nests = 1

!---communication switches
iopt_MPI_partit = 2

!---implicit switches
iopt_mg_smoother = 1
iopt_mg_prolong = 2
iopt_mg_cycle = 1

!---harmonic analysis
iopt_out_anal = MERGE(0,1,test)

!
!3. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
CStartDateTime(1:19) = '2006/04/01;00:00:00'
IF (test) THEN
   CEndDateTime(1:19) = '2006/04/05;00:00:00'
ELSE
   CEndDateTime(1:19) = '2006/05/01;00:00:00'
ENDIF

!---time step
delt2d = MERGE(20.0,300.0,iopt_hydro_impl.EQ.0)

!---counter for 3-D mode
ic3d = MERGE(15,1,iopt_hydro_impl.EQ.0)

!---counter for update nodal factors and tidal phases
icnodal = 0

!
!4. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 158; nr = 206; nz = 20

!---number of open sea boundaries
nosbu = 36; nosbv = 127

!---number of nested grids
nonestsets = 1

!---restart times
norestarts = 1
ntrestart(1) = 0

!---reference salinity [PSU]
sal_ref = 33.0
!---reference temperature [deg C]
temp_ref = 12.0

!---depth flag
depmean_flag = -99.9

!---uniform bottom roughness length
zrough_cst = 0.0035

!---drying/flooding depth parameters
dcrit_fld = 0.5; dthd_fld = 0.1; dmin_fld = 0.02

!---multigrid parameters
nomglevels = 3
ur_smooth = 0.8
mg_tol = 1.0E-06

!---relaxation distance for momentum advection
distrlx_obc = 20000.0

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
modfiles(io_modgrd,1,1)%filename = 'nos_grid.dat'
!---open boundary conditions
modfiles(io_2uvobc,1,1)%status = 'N'
modfiles(io_2uvobc,2,1)%status = 'R'
modfiles(io_2uvobc,2,1)%form = 'N'
modfiles(io_2uvobc,2,1)%filename = 'OUTDIR/'//'optcsm_obc.nc'
icount = MERGE(15,1,iopt_hydro_impl.EQ.0)
modfiles(io_2uvobc,2,1)%tlims = (/0,int_fill,icount/)
modfiles(io_2uvobc,2,1)%endfile = 2

!---domain decomposition
IF (iopt_MPI_partit.EQ.2) THEN
   modfiles(io_mppmod,1,1)%status = 'R'
   modfiles(io_mppmod,1,1)%form = 'A'
   modfiles(io_mppmod,1,1)%filename = 'optnos11.mppmod.txt'
ENDIF

!---nesting
IF (iopt_nests.EQ.1) THEN
   modfiles(io_nstspc,1,1)%status = 'N'
   modfiles(io_nstrel,1,1)%status = 'N'
   modfiles(io_nstrel,1,1)%form = 'A'
   modfiles(io_nstrel,1,1)%filename = 'nos_nest.dat'
ENDIF

!
!6.2 Output
!----------
!

IF (ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,2)%status = 'W'
   modfiles(io_modgrd,1,2)%form = 'U'
!  ---open boundary conditions
   modfiles(io_2uvobc,1,2)%status = 'W'
   modfiles(io_2uvobc,1,2)%form = 'A'
!  ---nesting
   IF (iopt_nests.EQ.1) THEN
      modfiles(io_nstspc,1,2)%status = 'W'
      modfiles(io_nstspc,1,2)%form = 'A'
      modfiles(io_nstrel,1,2)%status = 'W'
      modfiles(io_nstrel,1,2)%form = 'U'
   ENDIF
ENDIF

!---nesting
IF (iopt_nests.EQ.1) THEN
   modfiles(io_2uvnst,1,2)%status = 'W'
   modfiles(io_2uvnst,1,2)%form = 'N'
   modfiles(io_2uvnst,1,2)%filename = 'OUTDIR/'//runtitle(1:6)//'_obc.nc'
   icount = MERGE(15,1,iopt_hydro_impl.EQ.0)
   modfiles(io_2uvnst,1,2)%tlims = (/0,int_fill,icount/)
ENDIF

!---wait time
nowaitsecs = 5
maxwaitsecs = 600

!
!7. Surface grid parameters
!--------------------------
!

surfacegrids(igrd_model,1)%x0dat = -4.041667
surfacegrids(igrd_model,1)%y0dat = 48.47917
surfacegrids(igrd_model,1)%delxdat = 0.8333334E-01
surfacegrids(igrd_model,1)%delydat = 0.4166667E-01

!
!8. User output
!--------------
!
!---time series
novarstsr = 3
IF (iopt_verif.EQ.0) THEN
   nosetstsr = 2
   nostatstsr = 17
ELSE
   nosetstsr = 1
   nostatstsr = 0
ENDIF

!---harmonic analysis
IF (iopt_out_anal.EQ.1) THEN
   nosetsanal = 2
   nofreqsanal = 8
   novarsanal = 7
   nostatsanal = 17
ENDIF

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - optos nos
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
USE switches
USE inout_routines, ONLY: close_filepars, open_filepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: ii, i, iunit, j, jj
REAL :: depmin


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!1. Open grid file
!-----------------
!

CALL open_filepars(modfiles(io_modgrd,1,1))
iunit = modfiles(io_modgrd,1,1)%iunit
READ (iunit,'(/)')

!
!2. Water depths
!---------------
!

j_210: DO j=1,nr-1
   READ (iunit,*) depmeanglb(1:nc-1,j)
ENDDO j_210

depmin = MERGE(10.0,1.0,iopt_fld.EQ.0)
j_220: DO j=1,nr-1
i_220: DO i=1,nc-1
   IF (depmeanglb(i,j).GT.0.0.AND.depmeanglb(i,j).LT.depmin) THEN
      depmeanglb(i,j) = depmin
   ENDIF
ENDDO i_220
ENDDO j_220

!
!3. Open boundary locations
!--------------------------
!
!---U-nodes
READ (iunit,*)
ii_310: DO ii=1,nobu
   READ (iunit,*) iobu(ii), jobu(ii)
ENDDO ii_310

!---V-nodes
READ (iunit,*)
jj_320: DO jj=1,nobv
   READ (iunit,*) iobv(jj), jobv(jj)
ENDDO jj_320

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
! Description - optos nos
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
!*Local variables
!
INTEGER :: l


procname(pglev+1) = 'usrdef_2dobc_spec'
CALL log_timer_in()

!
!1. Type of o.b.c
!----------------
!

ityp2dobu = 11
ityp2dobv = 11

iloczobu = 1
iloczobv = 1

!
!2. Data file specifiers
!-----------------------
!

no2dobuv(2) = nobu+nobv
iobc2dtype(2) = 1
index2dobuv(:,2) =  (/(l,l=1,nobu+nobv)/)

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
