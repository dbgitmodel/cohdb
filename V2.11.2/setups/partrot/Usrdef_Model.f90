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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11
!
! $Date: 2016-10-27 16:15:25 +0200 (Thu, 27 Oct 2016) $
!
! $Revision: 986 $
!
! Description - test case partrot
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11
!
! Description - test case partrot 
!
! Reference -
!
! Calling program - simulation_start
!
!************************************************************************
!
USE iopars
USE paralpars
USE switches

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
!7. Particle model
!-----------------
!

SELECT CASE (runtitle(8:8))
  CASE ('A','B'); iopt_part_model = MERGE(2,1,npworld.GT.1)
  CASE ('C'); iopt_part_model = 3
  CASE ('D'); iopt_part_model = 4
END SELECT

nprocscoh = 1

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11
!
! Description - test case partrot
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
!---grid dimension
iopt_grid_nodim = 2

!---currents
iopt_curr = 1

!---bottom stress formulation (0/1/2)
iopt_bstres_form = 0

!---particle model
iopt_part_write = MERGE(1,0,runtitle(8:8).EQ.'A')

!---output
iopt_out_tsers = 0

!
!3. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
IF (iopt_part_model.LT.4) THEN
   CStartDateTime(1:19) = '2003/01/01;00:00:00'
   CEndDateTime(1:19) = '2003/01/02;00:00:00'
ELSE
   CStartDateTime(1:19) = '2003/01/02;00:00:00'
   CEndDateTime(1:19) = '2003/01/01;00:00:00'
ENDIF
   
!---time step (milliseconds)
delt2d = 7.5

!
!4. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 82; nr = 82

!---uniform mean water depth
depmean_cst = 6.0

!
!6. Model I/O file properties
!----------------------------
!
!6.1 Input
!---------
!
!---model grid
modfiles(io_modgrd,1,1)%status = '0'
!---initial conditions
modfiles(io_inicon,ics_phys,1)%status = 'N'
modfiles(io_inicon,ics_part,1)%status = 'N'

!
!6.2 Output
!----------
!

norestarts = 0

IF (ciffile%status.EQ.'W') THEN
!  ---restart times
   norestarts = 1
   ntrestart(1) = 0
!  ---initial conditions (physics)
   modfiles(io_fincon,ics_phys,2)%status = 'W'
   modfiles(io_fincon,ics_phys,2)%form = 'N'
!  ---initial conditions (particles)
   modfiles(io_fincon,ics_part,2)%status = 'W'
   modfiles(io_fincon,ics_part,2)%form = 'N'
ENDIF

!---writing data for particle module
IF (iopt_part_write.EQ.1) THEN
   modfiles(io_pargrd,1,2)%status = 'W'
   modfiles(io_pargrd,1,2)%form = 'N'
   modfiles(io_parphs,1,2)%status = 'W'
   modfiles(io_parphs,1,2)%form = 'N'
   modfiles(io_parphs,1,2)%tlims = (/0,int_fill,480/)
ENDIF

!---reading physical data from COHERENS
IF (iopt_part_model.GE.3) THEN
   modfiles(io_pargrd,1,1)%status = 'R'
   modfiles(io_pargrd,1,1)%form = 'N'
   modfiles(io_pargrd,1,1)%filename = 'partrotA.pargrd.nc'
   modfiles(io_parphs,1,1)%status = 'R'
   modfiles(io_parphs,1,1)%form = 'N'
   modfiles(io_parphs,1,1)%filename = 'partrotA.parphs.nc'
ENDIF
   
!
!7. Surface grid parameters
!--------------------------
!

surfacegrids(igrd_model,1)%delxdat = 0.5
surfacegrids(igrd_model,1)%delydat = 0.5
surfacegrids(igrd_model,1)%x0dat = 0.0
surfacegrids(igrd_model,1)%y0dat = 0.0

!
!8. User output
!--------------
!


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
! Description - test case partrot
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
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - test case partrot
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.11
!
! Description - test case partrot
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
!************************************************************************
!
USE currents
USE density
USE depths
USE grid
USE gridpars
USE iopars
USE switches
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j, k
REAL :: dh, dist, omega, rads0, x, xs0, y, ys0, width


procname(pglev+1) = 'usrdef_phsics'
CALL log_timer_in()

!
!1. Currents
!-----------
!
!---angular speed
omega = pi/3600.0
!---grid spacing and basin width
dh = 0.5

!---X-current
j_110: DO j=1,nrloc
i_110: DO i=1,ncloc
   IF (node2du(i,j).GT.1) THEN
      umvel(i,j) = -omega*(j-41)*dh
      uvel(i,j,:) = umvel(i,j)
      udvel(i,j) = umvel(i,j)*deptotatu(i,j)
   ENDIF
ENDDO i_110
ENDDO j_110

!---Y-current
j_120: DO j=1,nrloc
i_120: DO i=1,ncloc
   IF (node2dv(i,j).GT.1) THEN
      vmvel(i,j) = omega*(i-41)*dh
      vvel(i,j,:) = vmvel(i,j)
      vdvel(i,j) = vmvel(i,j)*deptotatv(i,j)
   ENDIF
ENDDO i_120
ENDDO j_120

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
INTEGER, INTENT(INOUT), DIMENSION(nobux) :: itypobux
INTEGER, INTENT(INOUT), DIMENSION(nobvy) :: itypobvy
INTEGER, INTENT(INOUT), DIMENSION(nobux,novars) :: iprofobux
INTEGER, INTENT(INOUT), DIMENSION(nobvy,novars) :: iprofobvy
INTEGER, INTENT(INOUT), DIMENSION(novars*(nobux+nobvy),2:nofiles) :: indexprof
INTEGER, INTENT(INOUT), DIMENSION(novars*(nobux+nobvy),2:nofiles) :: indexvar
INTEGER, INTENT(INOUT), DIMENSION(norlxzones) :: iprofrlx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id
!*itypobux*  INTEGER Type of U- or X-open boundary condition
!*itypobvy*  INTEGER Type of V- or Y-open boundary condition
!*iprofobux* INTEGER Profile numbers at U- or X-open boundaries
!*iprofobvy* INTEGER Profile numbers at V- or Y-open boundaries
!*iprofrlx*  INTEGER Disables/enables relaxation at open boundary zones
!*noprofsd*  INTEGER Number of profiles per data file
!*indexprof* INTEGER Mapping array of the profile numbers in the data files to
!                    the profile numbers assigned to the open boundaries. The
!                    physical size of the first dimension equals the number of
!                    profiles in a data file.
!*indexvar*  INTEGER Defines the variable number of the profiles in a data
!                    file. The physical size of the first dimension equals the
!                    number of profiles in a data file.
!*novars*    INTEGER Total number of variables
!*nofiles*   INTEGER Number of data files (+1)
!*nobux*     INTEGER Number of nodes at U- or X-open boundaries
!*nobvy*     INTEGER Number of nodes at V- or Y-open boundaries
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.2
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
