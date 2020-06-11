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
! Author - ANTEA and MUMM
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! $Date: 2018-02-21 09:29:52 +0100 (Wed, 21 Feb 2018) $
!
! $Revision: 1097 $
!
! Description - test case flood2d
!             - flooding/drying channel experiments
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
! Author - ANTEA and MUMM
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.3
!
! Description - test case flood2d
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
! Author - ANTEA and MUMM
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.4.2
!
! Description - test case flood2d
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
USE tide
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
REAL :: dhun, domsize, eps, href


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
!---inundation scheme
SELECT CASE (runtitle(8:8))
   CASE ('A','C'); iopt_fld = 1
   CASE ('B','D')
      iopt_fld = 2; fld_mask(1) = 1
END SELECT

!---formulation for bottom drag coefficient
iopt_bstres_drag = 3

!---2-D mode o.b.c (0/1)
iopt_obc_2D = 1

!
!3. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
CStartDateTime(1:23) = '2003/01/01;00:00:00:000'
CEndDateTime(1:23) = '2003/01/02;00:00:00:000'

!---time step and counter for 3-D mode
SELECT CASE (runtitle(8:8))
   CASE ('A','B')
      delt2d = 3.0; ic3d = 5
   CASE ('C','D')
      delt2d = 3.0; ic3d = 1
END SELECT

!
!4. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 101; nr = 2; nz = 20

!---number of open sea boundaries
nosbu = 1

!---number of tidal constituents
nconobc = 1

!---restart times
norestarts = 1
ntrestart(1) = 0

!---acceleration of gravity [m/s^2]
gacc_ref = 9.81

!---uniform roughness length
zrough_cst = 0.001

!---critical and minimum water depths
dcrit_fld = 0.1; dthd_fld = 0.1; dmin_fld = 0.02

!--bathymetry
depmean_flag = -99.0

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
modfiles(io_modgrd,1,1)%form = 'A'
IF (LLT(runtitle(8:8),'C')) THEN
   modfiles(io_modgrd,1,1)%filename='flood2dbath.dat'
ENDIF
!---initial conditions
modfiles(io_inicon,ics_phys,1)%status = '0'
!---open boundary conditions (2-D)
modfiles(io_2uvobc,1,1)%status = 'N'

!
!5.2 Output
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
!  ---initial conditions
   modfiles(io_fincon,ics_phys,2)%status = 'W'
   modfiles(io_fincon,ics_phys,2)%form = 'U'
!  ---open boundary conditions (2-D)
   modfiles(io_2uvobc,1,2)%status = 'W'
   modfiles(io_2uvobc,1,2)%form = 'U'
ENDIF

!
!6. Surface grid parameters
!--------------------------
!

dhun = 68.09864
surfacegrids(igrd_model,1)%delxdat = dhun
surfacegrids(igrd_model,1)%delydat = dhun
surfacegrids(igrd_model,1)%x0dat = 0.0
surfacegrids(igrd_model,1)%y0dat = 0.0

!
!7. User output
!--------------
!

nosetstsr = MERGE(2,1,iopt_verif.EQ.0)
novarstsr = 10

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_mod_params

!========================================================================

SUBROUTINE usrdef_grid
!************************************************************************
!
! *usrdef_grid* Define model grid arrays
!
! Author - ANTEA and MUMM
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.3
!
! Description - test case flood2d
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
USE physpars
USE inout_routines, ONLY : close_filepars, open_filepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!

INTEGER :: i, iunit, n1, n2, n3, n4
REAL :: hleft, hmin, hright


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!1. Bathymetry
!-------------
!

depmeanglb = depmean_flag


!1.1 Channel with sloping sea bed
!--------------------------------
!

IF (LLT(runtitle(8:8),'C')) THEN
   CALL open_filepars(modfiles(io_modgrd,1,1))
   iunit = modfiles(io_modgrd,1,1)%iunit
   READ (iunit,*) depmeanglb(1:nc-1,1)
   CALL close_filepars(modfiles(io_modgrd,1,1))

!
!1.2 Channel flow over an obstacle
!---------------------------------
!

ELSE
   hleft = 9.0; hmin = -1.0; hright = 3.0
   n1 = 30; n2 = 45; n3 = 55; n4 = 70
   i_120: DO i=1,nc-1
      IF (i.LE.n1) THEN
         depmeanglb(i,1) = hleft
      ELSEIF (i.LE.n2) THEN
         depmeanglb(i,1) = ((n2-i)*hleft+(i-n1)*hmin)/REAL(n2-n1)
      ELSEIF (i.LE.n3) THEN
         depmeanglb(i,1) = hmin
      ELSEIF (i.LE.n4) THEN
         depmeanglb(i,1) = ((n4-i)*hmin+(i-n3)*hright)/REAL(n4-n3)
      ELSE
         depmeanglb(i,1) = hright
      ENDIF
   ENDDO i_120
   
ENDIF

!
!2. Open boundary locations
!--------------------------
!

iobu = 1; jobu = 1

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_grid

!========================================================================

SUBROUTINE usrdef_partition
!************************************************************************
!
! *usrdef_partition* Define domain decomposition
!
! Author -
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.3
!
! Description - test case flood2d
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
! Author - ANTEA and MUMM
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.3
!
! Description - test case flood2d
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
! Author -
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.3
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
! Author - ANTEA and MUMM
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - test case flood2d
!
! Reference -
!
! Calling program - define_2dobc_spec
!
!************************************************************************
!
USE iopars
USE obconds
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


!---type of conditions
ityp2dobu = 3
iloczobu = 1; iloczobv = 1

!---amplitudes and phases
zdatobu_amp(:,1) = 3.0
zdatobu_pha(:,1) = -0.5*pi


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
! Author -
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
! Author -
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.3
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
! Author -
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
! Author -
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.3
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
! Author -
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.3
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
