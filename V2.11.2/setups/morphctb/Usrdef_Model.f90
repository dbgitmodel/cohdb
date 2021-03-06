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
! *Usrdef_Model*
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.10
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description - test case morphctb
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
! Description -
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

!
!6. Timing
!---------
!
!---timer level (0/1/2/3)
levtimer = 2


RETURN

END SUBROUTINE usrdef_init_params

!========================================================================

SUBROUTINE usrdef_mod_params
!************************************************************************
!
! *usrdef_mod_params*  Define parameters for physical model
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Last update - 4 Nov 2008  @(COHERENS)Usrdef_Model.f90  V2.10
!
! Description - test case morphctb
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
USE turbpars
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

sflag = MERGE(.TRUE.,.FALSE.,runtitle(9:9).EQ.'0')

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
!2. Switches
!-----------
!
!---sediment and morphology
iopt_sed = 1
iopt_morph = MERGE(0,1,sflag)

!---formulation for bottom drag coefficient
iopt_bstres_drag = 3

!---2-D mode o.b.c
iopt_obc_2D = 1

!---user outoput (0/1)
iopt_out_tsers = MERGE(0,1,sflag)

!
!3. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
IF (sflag) THEN
   CStartDateTime(1:23) = '2009/04/01;00:00:00:000'
   CEndDateTime(1:23)   = '2009/04/03;00:00:00:000'
ELSE
   CStartDateTime(1:23) = '2009/04/03;00:00:00:000'
   CEndDateTime(1:23)   = '2009/04/06;00:00:00:000'
ENDIF

!---time step (seconds)
delt2d = 15.0

!---counter for 3-D mode
ic3d = 20

!
!4. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 61; nr = 2; nz = 10

!---number of open sea boundaries
nosbu = 2*(nr-1)

!---number of tidal constituents
nconobc = 1

!---tidal indices
index_obc(1) =  icon_S2

!---uniform mean water depth
depmean_cst = 10.0

!---uniform bottom roughness length,
zrough_cst = 2.5E-06

!
!6. Model I/O file properties
!----------------------------
!
!6.1 Model input
!---------------
!
!---model grid
modfiles(io_modgrd,1,1)%status='N'

!---initial conditions
modfiles(io_inicon,ics_phys,1)%status = MERGE('0','R',sflag)
modfiles(io_inicon,ics_phys,1)%form = 'N'
modfiles(io_inicon,ics_phys,1)%filename = MERGE('         ','phsics.nc',sflag)
modfiles(io_inicon,ics_sed,1)%status = MERGE('0','R',sflag)
modfiles(io_inicon,ics_sed,1)%filename = MERGE('         ','sedics.nc',sflag)
modfiles(io_inicon,ics_sed,1)%form = 'N'
modfiles(io_inicon,ics_morph,1)%status = MERGE('0','N',sflag)

!---open boundary conditions
modfiles(io_2uvobc,1,1)%status = 'N'

!---sediment particle attributes
modfiles(io_sedspc,1,1)%status = 'N'

!
!6.2 Output
!----------
!
!---initial conditions (physics/sediments)
IF (sflag) THEN
   norestarts = 1
   ntrestart(1) = int_fill
   modfiles(io_fincon,ics_phys,2)%status = 'W'
   modfiles(io_fincon,ics_phys,2)%form = 'N'
   modfiles(io_fincon,ics_phys,2)%filename = 'phsics.nc'
   modfiles(io_fincon,ics_sed,2)%status = 'W'
   modfiles(io_fincon,ics_sed,2)%form =  'N'
   modfiles(io_fincon,ics_sed,2)%filename = 'sedics.nc'
ENDIF

IF (ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,2)%status = 'W'
   modfiles(io_modgrd,1,2)%form = 'N'
!  ---open boundary conditions
   modfiles(io_2uvobc,1,2)%status = 'W'
   modfiles(io_2uvobc,1,2)%form = 'N'
!  ---sediment particle attributes
   modfiles(io_sedspc,1,2)%status = 'W'
   modfiles(io_sedspc,1,2)%form = 'N'
!  ---initial conditions (morphology)
   IF (.NOT.sflag) THEN
      norestarts = 1
      ntrestart(1) = 0
      modfiles(io_fincon,ics_morph,2)%status = 'W'
      modfiles(io_fincon,ics_morph,2)%form =  'N'
      modfiles(io_fincon,ics_morph,2)%filename = 'morphics.nc'
   ENDIF
ENDIF

!
!7. Surface grid parameters
!--------------------------
!

surfacegrids(igrd_model,1)%delxdat = 400.0
surfacegrids(igrd_model,1)%delydat = 400.0
surfacegrids(igrd_model,1)%x0dat = 0.0
surfacegrids(igrd_model,1)%y0dat = 0.0

!
!7. User output
!--------------
!
!---title for user-output files
outtitle = runtitle(1:8)//runtitle(10:10)
!---time series output
IF (iopt_out_tsers.EQ.1) THEN
   nosetstsr = 1
   SELECT CASE (runtitle(10:10))
      CASE ('A','B')
         novarstsr = 10
      CASE ('C','D')
         novarstsr = 12
      CASE ('E','F','G','H')
         novarstsr = 13
   END SELECT
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
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.10
!
! Description - test case morphctb
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

!
!*Local variables
!
INTEGER :: i


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!2. Open boundary locations
!--------------------------
!
!---west    
iobu(1:nr-1) = 1
jobu(1:nr-1) = (/(i,i=1,nr-1)/)

!---east
iobu(nr:2*nr-2) = nc
jobu(nr:2*nr-2) = (/(i,i=1,nr-1)/)

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



RETURN

END SUBROUTINE usrdef_1dsur_spec

!========================================================================

SUBROUTINE usrdef_2dobc_spec(iddesc)
!************************************************************************
!
! *usrdef_2dobc_spec* Define specifier arrays for 2-D open boundary conditions
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.10
!
! Description - test case morphctb
!
! Reference -
!
! Calling program - define_2dobc_spec
!
!
!************************************************************************
!
USE gridpars
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
!*Local variables
!
INTEGER :: i


procname(pglev+1) = 'usrdef_2dobc_spec'
CALL log_timer_in()


!---type of open boundary conditions
ityp2dobu(1:nr-1) = 12
ityp2dobu(nr:2*nr-2) = 5

!---amplitudes and phases
zdatobu_amp(1:nr-1,1) = 1.5
zdatobu_pha(1:nr-1,1) = -halfpi

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
! Author - Alexander Btreugem (IMDC)
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.10
!
! Description - test case morphctb
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

procname(pglev+1) = 'usrdef_2dobc_data'
CALL log_timer_in()


IF (modfiles(io_2uvobc,ifil,1)%iostat.EQ.0) THEN
   modfiles(io_2uvobc,ifil,1)%iostat = 1
   GOTO 1000
ENDIF

ciodatetime = CStartDateTime

data2d = 10.0

1000  CALL log_timer_out()


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
