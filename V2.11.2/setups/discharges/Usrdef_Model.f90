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
!

!************************************************************************
!
! *Usrdef_Model* User-defined model setup
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - test case discharges 
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
! Description - test case discharges
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

levprocs_err = 2

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.6
!
! Description - test case discharges
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
!3. Switches
!-----------
!

!---equation of state (0/1/2)
iopt_dens = 2
!---formulation for baroclinic pressure gradient (0/1/2)
iopt_dens_grad = 0

!---salinity equation (0/1/2)
iopt_sal = 2

!---advection scheme
iopt_adv_2D = 1
iopt_adv_3D = 1

!---structures module (1/0)
iopt_dischr = 1

!
!4. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
CStartDateTime(1:23) = '2003/01/01;00:00:00:000'
CEndDateTime(1:23) = '2003/01/01;00:05:00:000'

!---time step (milliseconds)
delt2d = 2.5

!---counter for 3-D mode
ic3d = 4

!
!5. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 201; nr = 201; nz = 20

!---number of structures/open sea boundaries
nosbu = 0; nosbv = 0

!---acceleration of gravity [m/s^2]
gacc_ref = 9.81

!---reference salinity [PSU]
sal_ref = 35.01

!---uniform mean water depth
depmean_cst = 20.0

!---uniform bottom roughness length
zrough_cst = 0.006

!
!6. Structures
!-------------
!
!---discharges
numdis = 4

!
!7. Model I/O file properties
!----------------------------
!
!
!7.1 Input
!---------
!
!---discharge module
!--specifiers
modfiles(io_disspc,1,1)%status = 'N'
!--locations
modfiles(io_disloc,1,1)%status = 'N'
IF (runtitle(11:11).EQ.'C'.OR.runtitle(11:11).EQ.'D') THEN
   modfiles(io_disloc,1,1)%tlims = (/0,120,4/)
ELSE   
   modfiles(io_disloc,1,1)%tlims = (/0,0,1/)
ENDIF
!--discharge rates
modfiles(io_disvol,1,1)%status = 'N'
IF (runtitle(11:11).EQ.'D') THEN
   modfiles(io_disvol,1,1)%tlims = (/0,120,4/)
ELSE
   modfiles(io_disvol,1,1)%tlims = (/0,0,1/)
ENDIF
!--momentum discharge
IF (runtitle(11:11).EQ.'A') THEN
   modfiles(io_discur,1,1)%status = '0'
ELSE
   modfiles(io_discur,1,1)%status = 'N'
ENDIF
modfiles(io_discur,1,1)%tlims = (/0,0,1/)
!--salinity discharge
modfiles(io_dissal,1,1)%status = 'N'
IF (runtitle(11:11).EQ.'D') THEN
   modfiles(io_dissal,1,1)%tlims = (/0,120,-60/)
ELSE
   modfiles(io_dissal,1,1)%tlims = (/0,0,1/)
ENDIF

!
!7.2 Output
!----------
!

IF (ciffile%status.EQ.'W') THEN

!  ---discharge specifiers
   modfiles(io_disspc,1,2)%status = 'W'
   modfiles(io_disspc,1,2)%form = 'U'
!  ---discharge locations
   IF (LLT(runtitle(11:11),'C')) THEN
      modfiles(io_disloc,1,2)%status = 'W'
      modfiles(io_disloc,1,2)%form = 'U'
   ENDIF
!  ---discharge rates
   IF (LLT(runtitle(11:11),'D')) THEN
      modfiles(io_disvol,1,2)%status = 'W'
      modfiles(io_disvol,1,2)%form = 'U'
   ENDIF
!  ---momentum discharge
   IF (LGT(runtitle(11:11),'A')) THEN
      modfiles(io_discur,1,2)%status = 'W'
      modfiles(io_discur,1,2)%form = 'U'
   ENDIF
!  ---salinity discharge
   IF (LLT(runtitle(11:11),'D')) THEN
      modfiles(io_dissal,1,2)%status = 'W'
      modfiles(io_dissal,1,2)%form = 'U'
   ENDIF

ENDIF

!
!8. Surface grid parameters
!--------------------------
!

surfacegrids(igrd_model,1)%delxdat = 100.0
surfacegrids(igrd_model,1)%delydat = 100.0
surfacegrids(igrd_model,1)%x0dat = 0.0
surfacegrids(igrd_model,1)%y0dat = 0.0

!
!9. User output
!--------------
!
!---number of output files
nosetstsr = 2
novarstsr = 11

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

SUBROUTINE usrdef_2dobc_data(ifil,ciodatetime,data2d,nodat,novars)
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
! Description -
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
