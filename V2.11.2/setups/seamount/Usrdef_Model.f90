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
! $Date: 2018-07-23 16:55:25 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1171 $
!
! Description - test case seamount 
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
! Description - test case seamount
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
! Description - test case seamount 
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
!2. Switches
!-----------
!
!---model grid
iopt_grid_vtype = 3
!iopt_grid_vtype_transf = 21

!---type of array interpolation (0/1)
iopt_arrint_vreg = 1

!---implicit scheme
iopt_hydro_impl = 1

!---equation of state
iopt_dens = 1

!---formulation for baroclinic pressure gradient (1/2/3)
SELECT CASE (runtitle(9:9))
   CASE ('A','D'); iopt_dens_grad = 1
   CASE ('B','E'); iopt_dens_grad = 2
   CASE ('C','F'); iopt_dens_grad = 3
END SELECT
   
!---temperature equation (0/1/2)
iopt_temp = 1
!---optical module (0/1)
iopt_temp_optic = 0

!---bottom stress formulation (0/1/2)
iopt_bstres_form = 0

!---advection scheme for scalars (0/1/2/3/4)
iopt_adv_scal = 3
!---advection scheme for 3-D currents (0/1/2/3/4)
iopt_adv_3D = 3

!---type of vertical diffusion scheme (0/1/2/3)
iopt_vdif_coef = 0

!---2-D mode o.b.c
iopt_obc_2D = 1

!
!3. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
CStartDateTime(1:19) = '2001/01/01;00:00:00'
CEndDateTime(1:19) = '2001/01/08;00:00:00'

!---time step and 3-D counter
IF (iopt_hydro_impl.EQ.1) THEN
   delt2d = 300.0
ELSE
   SELECT CASE (runtitle(9:9))
      CASE ('A','B','C')
         delt2d = 15.0; ic3d = 40
      CASE ('D','E','F')
         delt2d = 7.5; ic3d = 80
   END SELECT
ENDIF

!
!4. Physical model constants
!---------------------------
!
!---grid dimensions
SELECT CASE (runtitle(9:9))
   CASE ('A','B','C')
      nc = 51; nr = 51; nz = 25
   CASE ('D','E','F')
      nc = 101; nr = 101; nz = 50
END SELECT

!---number of o.b. points
nosbu = 2*(nr-1); nosbv = 2*(nc-1)

!---acceleration of gravity [m/s^2]
gacc_ref = 9.81

!---reference latitude
dlat_ref = 45.0

!---reference salinity [PSU]
sal_ref = 33.0

!---reference temperature [deg C]
temp_ref = 22.0

!---implicit settings
nomglevels = MERGE(2,3,nc.EQ.51)

!
!6. Model I/O file properties
!----------------------------
!
!---model grid
modfiles(io_modgrd,1,1)%status = 'N'
!---initial conditions
modfiles(io_inicon,ics_phys,1)%status = 'N'
!---open boundary conditions (2-D)
modfiles(io_2uvobc,1,1)%status = 'N'

!
!7. Model grid
!-------------
!

dhun = 320000.0/REAL(nc-1)
surfacegrids(igrd_model,1)%delxdat = dhun
surfacegrids(igrd_model,1)%delydat = dhun
surfacegrids(igrd_model,1)%x0dat = 0.0
surfacegrids(igrd_model,1)%y0dat = 0.0

!
!8. User output
!--------------
!

nosetstsr = 2
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.2
!
! Description - test case seamount
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

!
!*Local variables
!
INTEGER :: i, j, k
REAL :: depmax = 5000.0, hmin = 500.0, hmount = 4500.0, hscale = 40000.0, &
      & theta = 3.0
REAL :: a, b, dhun, dist, sk, x, xcen, y, ycen


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!2. Water depths [m]
!-------------------
!

dhun = surfacegrids(igrd_model,1)%delxdat
xcen = 0.5*(nc-1)*dhun; ycen = xcen
i_210: DO i=1,nc-1
j_210: DO j=1,nr-1
   x = (i-0.5)*dhun; y = (j-0.5)*dhun
   dist = (x-xcen)**2+(y-ycen)**2
   depmeanglb(i,j) = depmax - hmount*EXP(-dist/(hscale**2))
ENDDO j_210
ENDDO i_210

!
!3. Sigma coordinates
!--------------------
!

k_310: DO k=2,nz
   sk = (k-1)/REAL(nz) - 1
   a = sk*hmin
   b = SINH(theta*sk)/SINH(theta)
   gscoordglb(1:nc-1,1:nr-1,k) = (a+b*(depmeanglb(1:nc-1,1:nr-1)-hmin))&
                               & /depmeanglb(1:nc-1,1:nr-1) + 1.0
ENDDO k_310

!
!4. Open boundary locations
!--------------------------
!
!---U-nodes
iobu(1:nr-1) = 1; iobu(nr:2*nr-2) = nc
jobu(1:nr-1) = (/(j,j=1,nr-1)/); jobu(nr:2*nr-2) = jobu(1:nr-1)

!---V-nodes
iobv(1:nc-1) = (/(i,i=1,nc-1)/); iobv(nc:2*nc-2) = iobv(1:nc-1)
jobv(1:nc-1) = 1; jobv(nc:2*nc-2) = nr

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
! Description - test case seamount
!
! Reference -
!
! Calling program - define_phsics
!
! Externals - equation_of_state
!
!************************************************************************
!
USE density
USE grid  
USE gridpars
USE physpars
USE grid_routines, ONLY: Zcoord_var
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

INTEGER :: i, j, k
REAL :: delta = 500.0, rhoamp = 3.0
REAL :: tamp, tsur, z


procname(pglev+1) =  'usrdef_phsics'
CALL log_timer_in()

!
!3.Density arrays
!----------------
!

CALL equation_of_state

tamp = rhoamp/(density_ref*beta_temp_ref)
tsur = temp_ref
j_310: DO j=1,nrloc
i_310: DO i=1,ncloc
   IF (maskatc_int(i,j)) THEN
      k_311: DO k=1,nz
         z = Zcoord_var(i,j,k,'C  ',.TRUE.)
         temp(i,j,k) = tsur - tamp*EXP(z/delta)
      ENDDO k_311
   ENDIF
ENDDO i_310
ENDDO j_310

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.2
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.2
!
! Description - test case seamount
!
! Reference -
!
! Calling program - define_2dobc_spec
!
!************************************************************************
!
USE iopars  
USE gridpars
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


!---type of conditions at open boundaries (0-13)
ityp2dobu = 6
ityp2dobv = 6

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

