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
! *Usrdef_Model* User-defined model setup for Hjelmfelt & Lenau (1970) testcase
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.10.3
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description - test case sedhprof
!   sedhprofA: 2DH, equation of Engelund-Hansen (1967)
!   sedhprofB: 2DH, equation of Ackers-White (1973)
!   sedhprofC: 2DH, equation of van Rijn (2003)
!   sedhprofD: 2DH, equation of Wu et al (2000)
!   sedhprofE: 2DH, Rouse profile integrated with 7 points, and boundary
!                   condition according to van Rijn (1984)
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.1
!
! Description - test case sedhprof
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
!---timer level (0/1/2/3)
levtimer = 2

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
! *usrdef_mod_params*  Define parameters for physical model
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.10.3
!
! Description - test case sedhprof
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


procname(pglev+1) = 'usrdef_mod_params'
CALL log_timer_in()

!
!1. Process numbers
!------------------
!

IF (npworld.GT.1) THEN
!  ---number of processes in X-direction
   nprocsx = 1
!  ---number of processes in Y-direction
   nprocsy = 4
ENDIF

!
!2. Switches
!-----------
!
!2.1 Model grid
!--------------
!
!---grid dimension (1/2/3)
iopt_grid_nodim = MERGE(3,2,runtitle(9:9).EQ.'F'.OR.runtitle(9:9).EQ.'I')

!
!2.3 Currents
!------------
!
!---currents (0/1/2)
iopt_curr = 1

!
!2.5 Sediments
!-------------
!
!---disable/enable sediments (0/1)
iopt_sed = 1

!
!2.6 Bottom boundary conditions
!------------------------------
!
!---formulation for bottom drag coefficient (0/1/2/3/4)
iopt_bstres_drag = 3
!---dimension of bottom stress formulation (2/3)
iopt_bstres_nodim = 2

!
!2.8 Advection
!-------------
!
!---advection schemes
iopt_adv_2D = 0
iopt_adv_scal = 3

!
!2.9 Diffusion
!-------------
!
!---type of vertical diffusion scheme (0/1/2/3)
iopt_vdif_coef = MERGE(0,1,iopt_grid_nodim.EQ.2)

!
!2.10 Turbulence
!---------------
!
!---type of algebraic turbulence closure (1/2/3/4/5/6)
iopt_turb_alg = 6

!
!2.13 Open Boundary conditions
!-----------------------------
!
!---sediment o.b.c.
iopt_obc_sed = 1

!
!3. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
CStartDateTime(1:23) = '2009/04/01;00:00:00:000'
CEndDateTime(1:23)   = '2009/04/01;03:00:00:000'

!---time step (seconds)
delt2d = 10.0

!
!4. Physical model constants
!---------------------------
!
!---model grid
nc = 2; nr = 121
nz = MERGE(1,30,iopt_grid_nodim.EQ.2)

!---number of open sea boundaries
nosbu = 0; nosbv = 2

!---uniform mean water depth
depmean_cst = 10.0

!---uniform bottom roughness length,
zrough_cst = 0.001

!
!6. Model I/O file properties
!----------------------------
!6.1 Input
!---------
!
!---model grid
modfiles(io_modgrd,1,1)%status='N'
!---initial conditions
modfiles(io_inicon,ics_phys,1)%status = 'N'
!---sediment particle attributes
modfiles(io_sedspc,1,1)%status = 'N'
!---sediment open boundary conditions
modfiles(io_sedobc,1:2,1)%status = 'N'

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
   modfiles(io_modgrd,1,2)%form = 'U'
!  ---sediment open boundary conditions
   modfiles(io_sedobc,1:2,2)%status = 'W'
   modfiles(io_sedobc,1:2,2)%form = 'U'
!  ---sediment particle attributes
   modfiles(io_sedspc,1,2)%status = 'W'
   modfiles(io_sedspc,1,2)%form = 'U'
ENDIF

!
!7. Surface grid parameters
!--------------------------
!
!---model grid
surfacegrids(igrd_model,1)%delxdat = 40.0
surfacegrids(igrd_model,1)%delydat = 40.0
surfacegrids(igrd_model,1)%x0dat = 0.0
surfacegrids(igrd_model,1)%y0dat = 0.0

!
!9. User output
!--------------
!
!---time series output
nosetstsr = 1
novarstsr = 2

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_mod_params

!========================================================================

SUBROUTINE usrdef_grid
!************************************************************************
!
! *usrdef_grid* Define model grid arrays
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.5
!
! Description - test case sedhprof
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!3. Open boundary locations
!--------------------------
!
!---V-nodes
iobv = 1
jobv = (/1,nr/)

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
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.5
!
! Description - test case sedhprof
!
! Reference -
!
! Calling program - define_phsics
!
!************************************************************************
!
USE currents
USE diffusion
USE grid
USE gridpars
USE iopars
USE physpars
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: k
REAL :: ustar, vmean 


procname(pglev+1) = 'usrdef_phsics'
CALL log_timer_in()

!---currents
vmean = 0.6
ustar = 0.028178
WHERE (node2dv(1:ncloc,1:nrloc).GT.1)
   vdvel(1:ncloc,1:nrloc) = vmean*depmean_cst
END WHERE
IF (iopt_grid_nodim.EQ.3) THEN
   k_110: DO k=1,nz
      WHERE (node2dv(1:ncloc,1:nrloc).GT.1)
         vvel(1:ncloc,1:nrloc,k) = vmean
      END WHERE
   ENDDO k_110
ENDIF

!---diffusion coefficient
IF (iopt_grid_nodim.EQ.3) THEN
   k_120: DO k=1,nz+1
      WHERE (maskatc_int)
         vdifcoefmom(1:ncloc,1:nrloc,k) =  ustar*depmean_cst*ckar* &
          & gscoordatw(1:ncloc,1:nrloc,k)*(1.0 - gscoordatw(1:ncloc,1:nrloc,k))
      END WHERE
   ENDDO k_120
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
! Calling program - define_1dobc_spec
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
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - test case sedhprof
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

noprofsd(2:nofiles) = 1
iprofobvy(1:nobv,1) = 1
indexprof = 1
indexvar = 1

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
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.5
!
! Description - test case sedhprof
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
!
!------------------------------------------------------------------------------
!

procname(pglev+1) = 'usrdef_profobc_data'
CALL log_timer_in()

IF (modfiles(iddesc,ifil,1)%iostat.EQ.0) THEN
   modfiles(iddesc,ifil,1)%iostat = 1
   GOTO 1000
ENDIF

ciodatetime = CStartDateTime

SELECT CASE (iddesc)
   CASE (io_sedobc)
      psiprofdat(1,:) = 0.0
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

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_rlxobc_spec
