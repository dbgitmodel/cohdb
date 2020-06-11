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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.2
!
! $Date: 2018-01-22 14:06:38 +0100 (Mon, 22 Jan 2018) $
!
! $Revision: 1077 $
!
! Description - test case convtest
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.2
!
! Description - test case convtest
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

IF (npworld.GT.1) nprocscoh = 4


RETURN

END SUBROUTINE usrdef_init_params

!========================================================================

SUBROUTINE usrdef_mod_params
!************************************************************************
!
! *usrdef_mod_params* Define parameters for control file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.2
!
! Description - test case convtest
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
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
REAL :: dhun


procname(pglev+1) = 'usrdef_mod_params'
CALL log_timer_in()

!
!1. Process numbers
!------------------
!

IF (npworld.GT.1) THEN
!  ---number of processes in X-direction
   nprocsx = 4
!  ---number of processes in Y-direction
   nprocsy = 1
ENDIF

!
!2. Switches
!-----------
!
!---current
iopt_curr = MERGE(0,2,runtitle(9:9).EQ.'E')

!---equation of state
iopt_dens = 2
!---temperature equation (0/1/2)
iopt_temp = 2
!---salinity
iopt_sal = 2
!---formulation for baroclinic pressure gradient (0/1/2)
iopt_dens_grad = MERGE(0,3,iopt_curr.EQ.0)
!---convective adjustment
iopt_dens_convect = MERGE(1,0,runtitle(9:9).GT.'C')

!---advection scheme for 2-D currents (0/1/2/3/4)
iopt_adv_2D = 3
!---advection scheme for 3-D currents (0/1/2/3/4)
iopt_adv_3D = 3
!---advection scheme for scalars (0/1/2/3/4)
iopt_adv_scal = 3

!---horizontal diffusion
iopt_hdif_coef = MERGE(0,1,runtitle(9:9).EQ.'A')
iopt_hdif_scal = 3
iopt_hdif_lim = 1

!---vertical diffusion
iopt_vdif_coef = 3

!--meteo
IF (runtitle(6:6).GT.'B') THEN
   iopt_meteo = 1
   iopt_meteo_data = 2
   iopt_meteo_stres = 0
ENDIF

!
!3. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
CStartDateTime(1:19) = '2005/01/01;00:00:00'
CEndDateTime(1:19) = '2005/02/01;00:00:00'

!---time step
delt2d = 30.0

!---counters
ic3d = 5

!
!4. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 51; nr = 2; nz = 60

!---depth flag
depmean_flag = -99.0

!---acceleration of gravity [m/s^2]
gacc_ref = 9.81

!---reference salinity [PSU]
sal_ref = 33.0

!---reference temperature [deg C]
temp_ref = 22.0

!---diffusion coefficient
hdifscal_cst = 10.0
vdifscal_cst = 0.0

!---isoneutral diffusion
rho_crit_iso = 0.01
skewdiff_cst = 0.0

!---bottom drag coefficient
zrough_cst = 0.006

!
!6. Model I/O file properties
!----------------------------
!
!6.1 Input
!---------
!
!---model grid
modfiles(io_modgrd,1,1)%status = 'N'
!---initial conditions
modfiles(io_inicon,ics_phys,1)%status = 'N'
!---meteo
IF (iopt_meteo.EQ.1) THEN
   modfiles(io_metsur,1,1)%status = 'N'
ENDIF

!
!6.2 Output
!----------
!

IF (ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,2)%status = 'W'
   modfiles(io_modgrd,1,2)%form = 'N'
!---restart times
   norestarts = 1
   ntrestart(1) = 0
!  ---initial conditions
   modfiles(io_fincon,ics_phys,2)%status = 'W'
   modfiles(io_fincon,ics_phys,2)%form = 'N'
ENDIF

!
!7. Model grid
!-------------
!

dhun = 1000.0
surfacegrids(igrd_model,1)%delxdat = dhun
surfacegrids(igrd_model,1)%delydat = dhun
surfacegrids(igrd_model,1)%x0dat = 0.0
surfacegrids(igrd_model,1)%y0dat = 0.0

!
!8. User output
!--------------
!

nosetstsr = 1
novarstsr = 6

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.2
!
! Description - test case convtest
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE depths
USE gridpars
USE iopars
USE physpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i
REAL :: a_d, b_d, c_d, dhun, length, x


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!---water depths
dhun = surfacegrids(igrd_model,1)%delxdat
length = (nc-1)*dhun 
a_d = 50.0; b_d = -10.0; c_d = (nc-1)*dhun
depmeanglb = depmean_flag
i_110: DO i=1,nc-1
   x = (i-0.5)*dhun
   depmeanglb(i,1) = a_d + b_d*TANH(c_d*(0.5*length-x))
ENDDO i_110

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.2
!
! Description - test case convtest
!
! Reference -
!
! Calling program - define_phsics
!
!************************************************************************
!
USE density
USE grid
USE gridpars
USE iopars  
USE physpars
USE switches
USE grid_routines, ONLY: Zcoord_var
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

INTEGER :: i, iglb, j, k
REAL :: a_h, b_h, c_h, b_v, c_v, dhun, length, x, z, zc


procname(pglev+1) =  'usrdef_phsics'
CALL log_timer_in()

!
!1.Temperature
!-------------
!

IF (iopt_temp.EQ.2) THEN
   dhun = surfacegrids(igrd_model,1)%delxdat
   length = (nc-1)*dhun
   a_h = -25.0; b_h = 10.0; c_h = 25.0/length
   b_v = 10.0; c_v = 0.1
   j_110: DO i=1,ncloc
   i_110: DO j=1,nrloc
      IF (maskatc_int(i,j)) THEN
         iglb = i+nc1loc-1
         x = (iglb-0.5)*dhun
         zc = a_h + b_h*TANH(c_h*(0.5*length-x))
         k_111: DO k=nz,nz/2,-1
            z = Zcoord_var(i,1,k,'C  ',.TRUE.)
            temp(i,j,k) = 12.0 + b_v*TANH(c_v*(z-zc))
         ENDDO k_111
         k_112: DO k=nz/2-1,1,-1
            temp(i,j,k) = temp(i,j,k+1)
         ENDDO k_112
      ENDIF
   ENDDO i_110
   ENDDO j_110
ENDIF

!
!2.Salinity
!----------
!

IF (iopt_sal.EQ.2) THEN
   dhun = surfacegrids(igrd_model,1)%delxdat
   length = (nc-1)*dhun
   a_h = -30.0; b_h = 10.0; c_h = 5.0/length
   b_v = 5.0; c_v = 0.01
   j_210: DO i=1,ncloc
   i_210: DO j=1,nrloc
      IF (maskatc_int(i,j)) THEN
         iglb = i+nc1loc-1
         x = (iglb-0.5)*dhun
         zc = a_h + b_h*TANH(c_h*(0.75*length-x))
         k_211: DO k=nz,nz/2,-1
            z = Zcoord_var(i,1,k,'C  ',.TRUE.)
            sal(i,j,k) = 30.0 - b_v*TANH(c_v*(z-zc))
         ENDDO k_211
         k_212: DO k=nz/2-1,1,-1
            sal(i,j,k) = sal(i,j,k+1)
         ENDDO k_212
      ENDIF
   ENDDO i_210
   ENDDO j_210
ENDIF

CALL log_timer_out()


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
!*iddesc*    INTEGER Data file id (io_2uvobc or io_2xyobc)
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


RETURN

END SUBROUTINE usrdef_rlxobc_spec

