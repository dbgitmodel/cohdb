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
! *Usrdef_Model* User-defined model setup for thacker test case
!
! Author - IMDC
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! $Date: 2017-10-10 13:57:00 +0200 (Tue, 10 Oct 2017) $
!
! $Revision: 1055 $
!
! Description - test case thacker
!   Flow and sediment transport in a cicular basin
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
! Author - IMDC
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.5
!
! Description - test case thacker
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
!2. Switches
!-----------
!
!---grid dimension (1/2/3)
SELECT CASE (runtitle(8:8))
  CASE ('A','B'); iopt_grid_nodim = 2
  CASE ('C','D'); iopt_grid_nodim = 3
END SELECT

!---disable/enable sediments (0/1)
iopt_sed = 1

!---dimension of bottom stress formulation
iopt_bstres_nodim = MERGE(2,3,iopt_grid_nodim.EQ.2)
!---formulation for bottom drag coefficient
iopt_bstres_drag = 3

!
!2.11 Flooding/drying
!--------------------
!
!---flooding/drying scheme (0/1/2)
iopt_fld = 1

!
!3. Date/time parameters
!-----------------------
!
!---Start/End date, time step
SELECT CASE (runtitle(8:8))
  CASE ('A','C')
     CStartDateTime(1:23) = '2010/04/01;00:00:00:000'
     CEndDateTime(1:23)   = '2010/04/01;01:00:00:000'
     delt2d = 1.0
  CASE ('B','D')
     CStartDateTime(1:23) = '2010/04/01;00:00:00:000'
     CEndDateTime(1:23)   = '2010/04/01;00:06:00:000'
     delt2d = 0.1
END SELECT

!
!4. Physical model constants
!---------------------------
!
!---model grid
nc = 111; nr = 111
nz = MERGE(1,10,iopt_grid_nodim.EQ.2)

!---drying/flooding depth parameters [m]
dmin_fld = 0.02

!---uniform bottom roughness length, 
zrough_cst = 0.001

!
!6. Model I/O file properties
!----------------------------
!
!6.1 Input
!---------
!
!---model grid
modfiles(io_modgrd,1,1)%status='N'
!---initial conditions
modfiles(io_inicon,ics_phys,1)%status = 'N'
!---sediment particle attributes
modfiles(io_sedspc,1,1)%status = 'N'

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
!  ---initial conditions
   modfiles(io_fincon,ics_phys,2)%status = 'W'
   modfiles(io_fincon,ics_phys,2)%form = 'U'
!  ---sediment particle attributes
   modfiles(io_sedspc,1,2)%status = 'W'
   modfiles(io_sedspc,1,2)%form = 'U'
ENDIF

!
!7. Surface grid parameters
!--------------------------
!

SELECT CASE (runtitle(8:8))
   CASE ('A','C')
      surfacegrids(igrd_model,1)%delxdat = 100.0
      surfacegrids(igrd_model,1)%delydat = 100.0
      surfacegrids(igrd_model,1)%x0dat = -500.0
      surfacegrids(igrd_model,1)%y0dat = -500.0
   CASE ('B','D')
      surfacegrids(igrd_model,1)%delxdat = 10.0
      surfacegrids(igrd_model,1)%delydat = 10.0
      surfacegrids(igrd_model,1)%x0dat = -50.0
      surfacegrids(igrd_model,1)%y0dat = -50.0
END SELECT

!
!9. User output
!--------------
!
!---time series output
nosetstsr = 1
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
! Author - IMDC
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.5
!
! Description - test case thacker
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
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!Local variables
!
INTEGER :: i, j
REAL :: dep, dx, dy, h0, x, radius, x0, y, y0


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!2. Water depths [m]
!-------------------
!

dx = surfacegrids(igrd_model,1)%delxdat
dy = surfacegrids(igrd_model,1)%delxdat
x0 = surfacegrids(igrd_model,1)%x0dat
y0 = surfacegrids(igrd_model,1)%y0dat
h0 = 10.0
radius = 0.5*(nc-11)*dx

i_210: DO i=1,nc-1
   x = x0 + (i-0.5)*dx - radius
   j_211: DO j=1,nr-1
      y = y0 + (j-0.5)*dy - radius
      dep = h0 + h0*(1-(x**2+y**2)/radius**2) 
      depmeanglb(i,j) = MERGE(dep,0.0,(dep.GT.0.0))   
   ENDDO j_211
ENDDO i_210

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
! Author - IMDC
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.5
!
! Description - test case thacker
!
! Reference -
!
! Calling program - define_phsics
!
! Module calls - distribute_mod, error_alloc
!
!************************************************************************
!
USE depths
USE gridpars
USE iopars
USE modids
USE syspars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: distribute_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!Local variables
!
INTEGER :: i, j
INTEGER, DIMENSION(4) :: nhdist
REAL :: dx, dy, h0, radius, x, x0, y, y0, zeta0
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: zetaglb


procname(pglev+1) = 'usrdef_phsics'
CALL log_timer_in()

!---allocate
ALLOCATE (zetaglb(0:nc+1,0:nr+1),STAT=errstat)
CALL error_alloc('zetaglb',2,(/nc+2,nr+2/),real_type)

!---initialise parameters
dx = surfacegrids(igrd_model,1)%delxdat
dy = surfacegrids(igrd_model,1)%delxdat
x0 = surfacegrids(igrd_model,1)%x0dat
y0 = surfacegrids(igrd_model,1)%y0dat
h0 = 10.0
radius = 0.5*(nc-11)*dx
SELECT CASE (runtitle(8:8))
   CASE ('A','C')
      zeta0 = 50.0
   CASE ('B','D')
      zeta0 = 20.0
END SELECT

!---initial surface elevation
i_110: DO i=1,nc-1
   x = x0 + (i-0.5)*dx - radius
   j_111: DO j=1,nr-1
      y = y0 + (j-0.5)*dy
      IF (depmeanglb(i,j).GT.0.0) THEN
         zetaglb(i,j)= 2.0*zeta0*(h0/(radius**2))*(x-0.5*zeta0) - 10.0
      ELSE
         zetaglb(i,j) = 0.0
      ENDIF
   ENDDO j_111
ENDDO i_110

!---distribute
nhdist = 1
CALL distribute_mod(zetaglb,zeta,(/0,0/),nhdist,iarr_zeta,0.0)

!---deallocate
DEALLOCATE (zetaglb)

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

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_rlxobc_spec
