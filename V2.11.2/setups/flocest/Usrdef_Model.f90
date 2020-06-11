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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.1
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - test case flocest
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.1
!
! Description - test case flocest
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
! *usrdef_mod_params* Define parameters for physical model
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.1
!
! Description - test case flocest
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE gridpars
USE iopars
USE obconds
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
LOGICAL :: sflag


procname(pglev+1) = 'usrdef_mod_params'
CALL log_timer_in()

!
!1. Setup flag
!-------------
!

sflag = MERGE(.TRUE.,.FALSE.,runtitle(8:8).EQ.'0')

!
!2. Process numbers
!------------------
!

IF (npworld.GT.1) THEN
!  ---number of processes in X-direction
   nprocsx = 4
!  ---number of processes in Y-direction
   nprocsy = 1
ENDIF

!
!2.Switches
!----------
!
!---equation of state (0/1/2)
iopt_dens = 2
!---formulation for baroclinic pressure gradient (0/1/2)
iopt_dens_grad = 1
!---salinity equation (0/1)
iopt_sal = 2

!---sediments
iopt_sed = MERGE(0,2,sflag)

!---advection scheme for 2-D currents (0/1/2/3/4)
iopt_adv_2D = 1
!---advection scheme for 3-D currents (0/1/2/3/4)
iopt_adv_3D = 3
!---advection scheme for scalars (0/1/2/3/4)
iopt_adv_scal = 3

!---turbulence
iopt_turb_ntrans = 2

!---salinity o.b.c (0/1)
iopt_obc_sal = 1
!---2-D mode o.b.c (0/1)
iopt_obc_2D = 1
!---3-D mode o.b.c (0/1)
iopt_obc_3D = 1

!---user output (0/1)
iopt_out_tsers = MERGE(0,1,sflag)

!
!3. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
IF (sflag) THEN
   CStartDateTime(1:19) = '2003/01/01;00:00:00'
   CEndDateTime(1:19)   = '2003/01/02;00:00:00'
ELSE
   CStartDateTime(1:19) = '2003/01/02;00:00:00'
   CEndDateTime(1:19)   = '2003/01/04;00:00:00'
ENDIF

!---time step
delt2d = 3.0   

!---counter for 3-D mode
ic3d = 10

!
!4. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 101; nr = 11; nz = 20

!---number of open boundaries
nosbu = nr-1; nosbv = 0
nrvbu = 2; nrvbv = 0

!---number of tidal constituents
nconobc = 1

!---reference latitude
dlat_ref = 52.0

!---acceleration of gravity [m/s^2]
gacc_ref = 9.81

!---reference salinity [PSU]
sal_ref = 35.0

!---uniform bottom roughness length
zrough_cst = 0.006

!---tidal indices
index_obc(1) = icon_S2

!
!5. Model I/O file properties
!----------------------------
!
!5.1 Input
!---------
!
!---model grid
modfiles(io_modgrd,1,1)%status = 'N'
!---open boundary conditions (2-D)
modfiles(io_2uvobc,1:2,1)%status = 'N'
!---open boundary conditions (3-D)
modfiles(io_3uvobc,1:2,1)%status = 'N'
!---open boundary conditions (salinity)
modfiles(io_salobc,1:2,1)%status = 'N'
!---sediment particle attributes
modfiles(io_sedspc,1,1)%status = MERGE('0','N',sflag)
!---initial conditions
IF (sflag) THEN
   modfiles(io_inicon,ics_phys,1)%status = '0'
   modfiles(io_inicon,ics_sed,1)%status = '0'
ELSE
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = 'N'
   modfiles(io_inicon,ics_phys,1)%filename = runtitle(1:7)//'.phsfin.nc'
   modfiles(io_inicon,ics_sed,1)%status = MERGE('N','0',runtitle(8:8).EQ.'B')
ENDIF

!
!5.2 Output
!----------
!

IF (ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,2)%status = 'W'
   modfiles(io_modgrd,1,2)%form = 'N'
!  ---open boundary conditions (2-D)
   modfiles(io_2uvobc,1:2,2)%status = 'W'
   modfiles(io_2uvobc,1:2,2)%form = 'N'
!  ---open boundary conditions (3-D)
   modfiles(io_3uvobc,1:2,2)%status = 'W'
   modfiles(io_3uvobc,1:2,2)%form = 'N'
!  ---open boundary conditions (salinity)
   modfiles(io_salobc,1:2,2)%status = 'W'
   modfiles(io_salobc,1:2,2)%form = 'N'
!  ---sediment particle attributes
   modfiles(io_sedspc,1,2)%status = 'W'
   modfiles(io_sedspc,1,2)%form = 'N'
!  ---sediment initial conditions
   IF (runtitle(8:8).EQ.'B') THEN
      modfiles(io_fincon,ics_sed,2)%status = 'W'
      modfiles(io_fincon,ics_sed,2)%form = 'N'
   ENDIF
ENDIF

!---initial conditions
IF (sflag) THEN
   norestarts = 1
   ntrestart(1) = int_fill
   modfiles(io_fincon,ics_phys,2)%status = 'W'
   modfiles(io_fincon,ics_phys,2)%form = 'N'
   modfiles(io_fincon,ics_phys,2)%filename = runtitle(1:7)//'.phsfin.nc'
ELSE
   norestarts = 1
   ntrestart(1) = 0
   modfiles(io_fincon,ics_phys,2)%status = '0'
ENDIF

!
!6. Surface grid parameters
!--------------------------
!

surfacegrids(igrd_model,1)%delxdat = 100.0
surfacegrids(igrd_model,1)%delydat = 100.0
surfacegrids(igrd_model,1)%x0dat = 0.0
surfacegrids(igrd_model,1)%y0dat = 0.0

!
!7. User output
!--------------
!
!---time series
nosetstsr = 2
novarstsr = 17

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - test case flocest
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE depths
USE physpars
USE grid
USE gridpars
USE physpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j
REAL :: hleft, hright


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!1. Water depths
!---------------
!

hleft = 15.0; hright = 5.0
j_110: DO j=1,nr-1
i_110: DO i=1,nc-1
   IF (i.LE.50) THEN
      depmeanglb(i,j) = ((50-i)*hleft+(i-1)*hright)/49.0
   ELSEIF (j.GE.5.AND.j.LE.6) THEN
      depmeanglb(i,j) = hright
   ELSE
      depmeanglb(i,j) = depmean_flag
   ENDIF
ENDDO i_110
ENDDO j_110

!
!2. Open boundary locations
!--------------------------
!
!---west
iobu(1:nosbu) = 1
jobu(1:nosbu) = (/(i,i=1,nr-1)/)

!---east
iobu(nobu-1:nobu) = nc
jobu(nobu-1:nobu) =  (/5,6/)

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
! *usrdef_phsics* Define arrays initial conditions for physics
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE density
USE grid
USE gridpars
USE physpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.1
!
! Description - test case flocest
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
USE physpars
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc

!
! Name      Type      Purpose
!-------------------------------------------------------------------------
!*iddesc*   INTEGER   Data file id
!
!-------------------------------------------------------------------------
!
!*Local variables
!
REAL :: delydat, hleft, zamp

procname(pglev+1) = 'usrdef_2dobc_spec'
CALL log_timer_in()


!---type of condition at open boundaries
ityp2dobu(1:nosbu) = 11
ityp2dobu(nobu-1:nobu) = 15

!---contents of data file
no2dobuv(2) = 2
iobc2dtype(2) = 3
index2dobuv(1:2,2) = (/nosbu+1,nosbu+2/)

!---amplitudes, phases
hleft = 15.0
zamp = 1.0
zdatobu_amp(1:nosbu,1) = zamp
zdatobu_pha(1:nosbu,1) = halfpi
udatobu_amp(1:nosbu,1) = zamp*SQRT(hleft*gacc_mean)
udatobu_pha(1:nosbu,1) = halfpi

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.1
!
! Description - test case flocest
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
   CASE (io_3uvobc)
      noprofsd = 1
      iprofobux(nobu-1:nobu,1) = 1
      indexprof(1,2) = 1
   CASE (io_salobc)
      noprofsd = 2
      iprofobux(1:nobu-2,1) = 1
      iprofobux(nobu-1:nobu,1) = 2
      indexprof(1:2,2) = (/1,2/)
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.1
!
! Description -
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
!*iddesc*      INTEGER Data file id
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'usrdef_2dobc_data'
CALL log_timer_in()

IF (modfiles(iddesc,2,1)%iostat.EQ.0) THEN
   modfiles(iddesc,2,1)%iostat = 1
   GOTO 1000
ENDIF

ciodatetime = CStartDateTime

data2d(:,1) = 30.0

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11.1
!
! Description - test case flocest
!
! Reference -
!
! Calling program - define_profobc_data
!
!************************************************************************
!
USE gridpars
USE iopars
USE physpars
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

IF (modfiles(iddesc,2,1)%iostat.EQ.0) THEN
   modfiles(iddesc,2,1)%iostat = 1
   GOTO 1000
ENDIF

ciodatetime = CStartDateTime

SELECT CASE (iddesc)
   CASE (io_3uvobc)
      psiprofdat(1,:) = 0.0
   CASE (io_salobc)
      psiprofdat(1,:) = sal_ref
      psiprofdat(2,:) = 0.0
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
