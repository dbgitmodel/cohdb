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
! $Date: 2017-10-10 13:57:00 +0200 (Tue, 10 Oct 2017) $
!
! $Revision: 1055 $
!
! Description - test cases flood3d
!             - flooding/drying experiment over a 3-D obstacle
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

! Author - ANTEA and MUMM
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.3
!
! Description - test case flood3d
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
! Description - test case flood3d
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
REAL :: dhun

procname(pglev+1) = 'usrdef_mod_params'
CALL log_timer_in()

!
!2. Switches
!-----------
!
!---grid dimension (1/2/3)
iopt_grid_nodim = MERGE(3,2,LLT(runtitle(8:8),'D'))

!---bottom stres
iopt_bstres_drag = 3
iopt_bstres_nodim = MERGE(3,2,iopt_grid_nodim.EQ.3)

!---inundation scheme
iopt_fld = 2

!---2-D mode o.b.c (0/1)
iopt_obc_2D = 1

!
!3. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
CStartDateTime(1:23) = '2003/01/01;00:00:00:000'
CEndDateTime(1:23) = '2003/01/02;00:00:00:000'

!---time step
delt2d = 3.0
SELECT CASE (runtitle(8:8))
   CASE ('A')
      delt2d = 3.0; ic3d = 10
   CASE ('B')
      delt2d = 3.0; ic3d = 5
   CASE ('C','D')
      delt2d = 3.0; ic3d = 1
END SELECT

!
!4. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 101; nr = 101
nz = MERGE(10,1,iopt_grid_nodim.EQ.3)

!---number of river open sea boundaries
nosbu = nr-1; nosbv = 0

!---number of tidal constituents
nconobc = 1

!---acceleration of gravity [m/s^2]
gacc_ref = 9.81

!---critical and minimum water depths
dcrit_fld = 0.10; dthd_fld = 0.1; dmin_fld = 0.02

!--bathymetry
depmean_flag = -99.0

!---uniform roughness length
zrough_cst = 0.001

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
!---initial conditions
modfiles(io_inicon,ics_phys,1)%status = '0'
!---open boundary conditions (2D)
modfiles(io_2uvobc,1,1)%status = 'N'

!
!5.2 Output
!----------
!

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
novarstsr = MERGE(12,9,iopt_grid_nodim.EQ.3)

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
! Description - test case flood3d
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - construct_rectgrid
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE physpars
USE syspars
USE grid_routines, ONLY: construct_rectgrid
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j, jj
REAL :: delxdat, delydat, domsize, domsizex, domsizey, hin, hmin, hout, href, &
      & pfac, r, r1, r2, r3, r4, xcen, x0dat, ycen, y0dat
REAL, DIMENSION(nc-1,nr-1) :: cendist, gxcoordatc, gycoordatc, hdelta


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!1. Coordinates at cell corners
!------------------------------
!

delxdat = surfacegrids(igrd_model,1)%delxdat
delydat = surfacegrids(igrd_model,1)%delydat
x0dat = surfacegrids(igrd_model,1)%x0dat + 0.5*delxdat
y0dat = surfacegrids(igrd_model,1)%y0dat + 0.5*delydat
CALL construct_rectgrid(x0dat,y0dat,delxdat,delydat,gxcoordatc,gycoordatc,&
                      & nc-1,nr-1)
domsizex = (nc-1)*delxdat; domsizey = (nr-1)*delydat
xcen = surfacegrids(igrd_model,1)%x0dat + 0.5*domsizex
ycen = surfacegrids(igrd_model,1)%y0dat + 0.5*domsizey
cendist = 2.0*(((gxcoordatc-xcen)/domsizex)**2+((gycoordatc-ycen)/domsizey)**2)

!
!2. Bathymetry
!-------------
!

depmeanglb = depmean_flag

SELECT CASE (runtitle(8:8))

   CASE ('A')
      hmin = 0.3; href = 10.0
      depmeanglb(1:nc-1,1:nr-1) = hmin + (href-hmin)*cendist

   CASE ('B')
      hmin = 0.3; href = 10.0; pfac = 1.5
      domsize = SQRT(domsizex**2+domsizey**2)
      hdelta = twopi*SQRT((gxcoordatc-xcen)**2+(gycoordatc-ycen)**2)/domsize
      hdelta = pfac*COS(hdelta)**3
      depmeanglb(1:nc-1,1:nr-1) = hmin + (href-hmin)*cendist + hdelta

   CASE ('C','D')
      hmin = 0.0; hin = 4.0; hout = 10.0
      r1 = 0.25; r2 = 0.4; r3 = 0.5; r4 = 0.6
      j_220: DO j=1,nr-1
      i_220: DO i=1,nc-1
         r = SQRT(cendist(i,j))
         IF (r.LE.r1) THEN
            depmeanglb(i,j) = hin
         ELSEIF (r.LE.r2) THEN
            depmeanglb(i,j) = ((r2-r)*hin+(r-r1)*hmin)/(r2-r1)
         ELSEIF (r.LE.r3) THEN
            depmeanglb(i,j) = hmin
         ELSEIF (r.LE.r4) THEN
            depmeanglb(i,j) = ((r-r3)*hout+(r4-r)*hmin)/(r4-r3)
         ELSE
            depmeanglb(i,j) = hout
         ENDIF
      ENDDO i_220
      ENDDO j_220
      
END SELECT



!
!3. Open boundary locations
!--------------------------
!

iobu = 1
jobu = (/(jj,jj=1,nr-1)/)

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
! Author - ANTEA and MUMM
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.3
!
! Description - test case flood3d
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
! Description - test case flood3d
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

ityp2dobu = 3
iloczobu = 1

zdatobu_amp = 2.0
zdatobu_pha = -0.5*pi

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

SUBROUTINE usrdef_profobc_data(iddesc,ifil,ciodatetime,psiprofdat,nzdat,&
                             & numvars,numprofs)
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
