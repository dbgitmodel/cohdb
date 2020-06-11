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
! Version - @(COHERENS)Usrdef_Model.f90  V2.x
!
! $Date: 2013-05-28 15:26:47 +0200 (Tue, 28 May 2013) $
!
! $Revision: 574 $
!
! Description - test case river
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
! Description - test case river
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
!3. Log files
!------------
!
!---program leveling in log files
iproc_310: DO iproc=1,npworld
   levprocs_ini(iproc) = 7
   levprocs_run(iproc) = 3
ENDDO iproc_310

!
!6. Timing
!---------
!

levtimer = 3


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
! Version - @(COHERENS)Usrdef_Model.f90  V2.4
!
! Description - test case river
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
USE tide
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: sflag
CHARACTER (LEN=1) :: form


procname(pglev+1) = 'usrdef_mod_params'
CALL log_timer_in()

!
!1. Setup flag
!-------------
!

sflag = MERGE(.TRUE.,.FALSE.,runtitle(6:6).EQ.'0')

!
!2. Process numbers
!------------------
!

IF (npworld.GT.1) THEN
!  ---total number of processes
   nprocs = 4
!  ---number of processes in X-direction
   nprocsx = 4
!  ---number of processes in Y-direction
   nprocsy = 1
ENDIF

!
!3. Switches
!-----------
!

IF (.NOT.sflag) THEN
!  ---equation of state (0/1/2)
   iopt_dens = 1
!  ---formulation for baroclinic pressure gradient (0/1/2)
   iopt_dens_grad = 1
ENDIF

!---salinity equation (0/1)
iopt_sal = MERGE(1,2,sflag)

!---advection scheme for scalars (0/1/2/3/4)
iopt_adv_scal = 3
!---advection scheme for 2-D/3-D currents (0,1,2,3,4)
iopt_adv_2D = 1
iopt_adv_3D = 1

!---type of background mixing scheme (0/1/2)
IF (.NOT.sflag) iopt_turb_iwlim = 1

!---2-D mode o.b.c (0/1)
iopt_obc_2D = 1

!---user output (0/1)
IF (sflag) iopt_out_tsers = 0

!
!4. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
IF (sflag) THEN
   CStartDateTime(1:19) = '2003/01/01;00:00:00'
   CEndDateTime(1:19) = '2003/01/03;00:00:00'
ELSE
   CStartDateTime(1:19) = '2003/01/03;00:00:00'
   CEndDateTime(1:19) = '2003/01/06;00:00:00'
ENDIF

!---time step
delt2d = 30.0

!---counter for 3-D mode
ic3d = 10

!
!5. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 141; nr = 2; nz = 20

!---number of river open sea boundaries
nosbu = 2

!---number of tidal constituents
nconobc = 1

!---restart times
IF (sflag) THEN
   norestarts = 1
   ntrestart(1) = int_fill
ENDIF

!---acceleration of gravity [m/s^2]
gacc_ref = 9.81

!---uniform mean water depth
depmean_cst = 20.0

!---uniform bottom roughness length
zrough_cst = 0.006

!---tidal indices
index_obc(1) = icon_S2

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
IF (sflag) THEN
   modfiles(io_inicon,ics_phys,1)%status = '0'
ELSE
   modfiles(io_inicon,ics_phys,1)%status = 'N'
   modfiles(io_inicon,ics_phys,1)%form = 'N'
   form = modfiles(io_inicon,ics_phys,1)%form
   modfiles(io_inicon,ics_phys,1)%filename = 'riverT.phsfin'//form
ENDIF
!---open boundary conditions (2-D)
modfiles(io_2uvobc,1,1)%status = 'N'

!
!6.2 Output
!----------
!
!---model grid
IF (sflag) THEN
   modfiles(io_modgrd,1,2)%status = 'W'
   modfiles(io_modgrd,1,2)%form = 'N'
ELSE
   modfiles(io_modgrd,1,2)%status = '0'
ENDIF
!---initial conditions
IF (sflag) THEN
   modfiles(io_fincon,ics_phys,2)%status = 'W'
   modfiles(io_fincon,ics_phys,2)%form = 'N'
   form = modfiles(io_fincon,ics_phys,1)%form
   modfiles(io_fincon,ics_phys,1)%filename = 'riverT.phsfin'//form
ENDIF
!---open boundary conditions (2-D)
IF (sflag) THEN
   modfiles(io_2uvobc,1,2)%status = 'W'
   modfiles(io_2uvobc,1,2)%form = 'A'
ELSE
   modfiles(io_2uvobc,1,2)%status = '0'
ENDIF

!
!7. Surface grid parameters
!--------------------------
!

surfacegrids(igrd_model,1)%delxdat = 1000.0
surfacegrids(igrd_model,1)%delydat = 1000.0
surfacegrids(igrd_model,1)%x0dat = 0.0
surfacegrids(igrd_model,1)%y0dat = 0.0

!
!8. User output
!--------------
!
!---title for forcing files
intitle = runtitle(1:5)//runtitle(7:7)
!---title for user-output files
outtitle = runtitle(1:5)//runtitle(7:7)
!---number of output files
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - test case river
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


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!3. Open boundary locations
!--------------------------
!
!---U-nodes
iobu = (/1,nc/)
jobu = 1

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.1
!
! Description - test case river
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
USE density
USE grid
USE gridpars
USE iopars
USE modids
USE switches
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j
REAL :: s, sl, sr, xc, xl, xr
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhexch


procname(pglev+1) =  'usrdef_phsics'
CALL log_timer_in()

!
!1. Read initial conditions from spin-up phase
!---------------------------------------------
!

!modfiles(io_inicon,1,1)%filename = TRIM(intitle)//'.phsfin'//&
!                                 & modfiles(io_inicon,1,1)%form
CALL read_phsics

!
!2. Salinity
!-----------
!

sl = 30.0; sr = 35.0
xl = 25000.0; xr = 45000.0

i_210: DO i=1,ncloc
j_210: DO j=1,nrloc
   IF (maskatc_int(i,j)) THEN
      xc = 0.5*(gxcoord(i,j)+gxcoord(i+1,j))
      IF (xc.LT.xl) THEN
         s = sl
      ELSEIF (xc.GT.xr) THEN
         s = sr
      ELSE
         s = (sr-sl)*(xc-xl)/(xr-xl)+sl
      ENDIF
      sal(i,j,:) = s
   ENDIF
ENDDO j_210
ENDDO i_210

!
!3. Exchange halos
!-----------------
!

IF (iopt_MPI.EQ.1) THEN
   lbounds = (/1-nhalo,1-nhalo,1/); nhexch = nhdens
   CALL exchange_mod(sal,lbounds,nhexch,iarr_sal)
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.x
!
! Description - test case river
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
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id (io_2uvobc or io_2xyobc)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL :: amp, crad, phase 


procname(pglev+1) = 'usrdef_2dobc_spec'
CALL log_timer_in()

!---type of conditions at open boundaries
ityp2dobu(1) = 11; ityp2dobu(2) = 13

!---location elevation point
iloczobu = 1; iloczobv = 1

!---amplitudes
crad = SQRT(gacc_mean*depmean_cst)
amp = 0.8; phase = -halfpi
udatobu_amp(1,1) = crad*amp
zdatobu_amp(1,1) = amp

!---phases
udatobu_pha(1,1) = phase
zdatobu_pha(1,1) = phase

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.x
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

SUBROUTINE usrdef_2dobc_data(ifil,ciodatetime,data2d,nodat,novars)
!************************************************************************
!
! *usrdef_2dobc_data* Define open boundary data for 2-D mode
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.x
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
INTEGER, INTENT(IN) :: ifil, nodat, novars
REAL, INTENT(INOUT), DIMENSION(nodat,novars) :: data2d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
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
