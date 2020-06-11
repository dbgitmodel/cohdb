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
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! $Date: 2018-06-05 16:03:05 +0200 (Tue, 05 Jun 2018) $
!
! $Revision: 1142 $
!
! Description - test case obctang
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - test case obctang
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
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - test case obctang
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
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=1) :: cexp


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
!---open boundary conditions
iopt_obc_2D = 1
SELECT CASE (runtitle(8:8))
   CASE('A','C')
      iopt_obc_3D = 0
      iopt_obc_2D_tang = 1
      iopt_obc_3D_tang = 0
   CASE('B','D')
      iopt_obc_3D = 1
      iopt_obc_2D_tang = 1
      iopt_obc_3D_tang = 1
   CASE('E','G')
      iopt_obc_3D = 0
      iopt_obc_2D_tang = 0
      iopt_obc_3D_tang = 0
   CASE('F','H')
      iopt_obc_3D = 1
      iopt_obc_2D_tang = 0
      iopt_obc_3D_tang = 0
END SELECT

!---nesting
SELECT CASE (runtitle(8:8))
  CASE ('A','C','E','G'); iopt_nests = 1
  CASE ('B','D','F','H'); iopt_nests = 0
END SELECT

!---harmonic analysis
iopt_out_anal = 1

!
!4. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
CStartDateTime(1:19) = '2003/01/01;00:00:00'
CEndDateTime(1:19) =   '2003/01/03;00:00:00'

!---time step (milliseconds)
delt2d = MERGE(30.0,15.0,iopt_nests.EQ.1)
ic3d = 10

!
!5. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 201; nr = 101; nz = 20

!---number of structures/open sea boundaries
nosbu = 200; nosbv = 400

!---number of tidal constituents
nconobc = 1

!---relaxation time for 2-D open boundary conditions
ntobcrlx = 1440

!---number of nested grids
IF (iopt_nests.EQ.1) nonestsets = 1

!---acceleration of gravity [m/s^2]
gacc_ref = 9.81

!---uniform mean water depth
depmean_cst = 20.0

!---uniform linear bottom drag
bdraglin = 0.001

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
!6.1.1 Model grid
!----------------
!

modfiles(io_modgrd,1,1)%status = 'N'

!
!6.1.2 Open boundary conditions
!------------------------------
!

SELECT CASE (runtitle(8:8))
   CASE ('B'); cexp = 'A'
   CASE ('D'); cexp = 'C'
   CASE ('F'); cexp = 'E'
   CASE ('H'); cexp = 'G'
END SELECT

!---2-D normal
modfiles(io_2uvobc,1,1)%status = 'N'
IF (iopt_nests.EQ.0) THEN
   modfiles(io_2uvobc,2,1)%status = 'R'
   modfiles(io_2uvobc,2,1)%filename = runtitle(1:7)//cexp//'_2uvobc.nc'
   modfiles(io_2uvobc,2,1)%form = 'N'
   modfiles(io_2uvobc,2,1)%tlims = (/0,int_fill,1/)
ENDIF

!---2-D tangential
IF (iopt_obc_2D_tang.EQ.1) THEN
   modfiles(io_2xyobc,1,1)%status = 'N'
   IF (iopt_nests.EQ.0) THEN
      modfiles(io_2xyobc,2,1)%status = 'R'
      modfiles(io_2xyobc,2,1)%filename = runtitle(1:7)//cexp//'_2xyobc.nc'
      modfiles(io_2xyobc,2,1)%form = 'N'
      modfiles(io_2xyobc,2,1)%tlims = (/0,int_fill,1/)
   ENDIF
ENDIF

!---3-D normal
IF (iopt_obc_3D.EQ.1) THEN
   modfiles(io_3uvobc,1,1)%status = 'N'
   modfiles(io_3uvobc,2,1)%status = 'R'
   modfiles(io_3uvobc,2,1)%filename = runtitle(1:7)//cexp//'_3uvobc.nc'
   modfiles(io_3uvobc,2,1)%form = 'N'
   modfiles(io_3uvobc,2,1)%tlims = (/0,int_fill,20/)
ENDIF

!---3-D tangential
IF (iopt_obc_3D_tang.EQ.1) THEN
   modfiles(io_3xyobc,1,1)%status = 'N'
   modfiles(io_3xyobc,2,1)%status = 'R'
   modfiles(io_3xyobc,2,1)%filename = runtitle(1:7)//cexp//'_3xyobc.nc'
   modfiles(io_3xyobc,2,1)%form = 'N'
   modfiles(io_3xyobc,2,1)%tlims = (/0,int_fill,20/)
ENDIF

!---nesting
IF (iopt_nests.EQ.1) THEN
!  ---specifiers
   modfiles(io_nstspc,1,1)%status = 'N'
!  ---nest locations
   modfiles(io_nstrel,1,1)%status = 'N'
ENDIF

!
!6.2 Output
!----------
!

IF (ciffile%status.EQ.'W') THEN

!
!6.2.1 Model grid
!----------------
!

   modfiles(io_modgrd,1,2)%status = 'W'
   modfiles(io_modgrd,1,2)%form = 'N'

!
!6.2.2 Open boundary conditions
!------------------------------
!
!  ---2-D normal
   modfiles(io_2uvobc,1,2)%status = 'W'
   modfiles(io_2uvobc,1,2)%form = 'N'
!  ---2-D tangential
   IF (iopt_obc_2D_tang.EQ.1) THEN
      modfiles(io_2xyobc,1,2)%status = 'W'
      modfiles(io_2xyobc,1,2)%form = 'N'
   ENDIF
!  ---3-D normal
   IF (iopt_obc_3D.EQ.1) THEN 
      modfiles(io_3uvobc,1,2)%status = 'W'
      modfiles(io_3uvobc,1,2)%form = 'N'
   ENDIF
!  ---3-D tangential
   IF (iopt_obc_3D_tang.EQ.1) THEN 
      modfiles(io_3xyobc,1,2)%status = 'W'
      modfiles(io_3xyobc,1,2)%form = 'N'
   ENDIF
!  ---nesting
   IF (iopt_nests.EQ.1) THEN
      modfiles(io_nstspc,1,2)%status = 'W'
      modfiles(io_nstspc,1,2)%form = 'A'
      modfiles(io_nstrel,1,2)%status = 'W'
      modfiles(io_nstrel,1,2)%form = 'N'
   ENDIF

ENDIF

!
!6.3 Nesting
!-----------
!

IF (iopt_nests.EQ.1) THEN

!  ---2-D normal
   modfiles(io_2uvnst,1,2)%status = 'W'
   modfiles(io_2uvnst,1,2)%form = 'N'
   modfiles(io_2uvnst,1,2)%filename = runtitle(1:8)//'_2uvobc.nc'
   modfiles(io_2uvnst,1,2)%tlims = (/0,int_fill,2/)

!  ---2-D tangential
   IF (runtitle(8:8).EQ.'A'.OR.runtitle(8:8).EQ.'C') THEN
      modfiles(io_2xynst,1,2)%status = 'W'
      modfiles(io_2xynst,1,2)%form = 'N'
      modfiles(io_2xynst,1,2)%filename = runtitle(1:8)//'_2xyobc.nc'
      modfiles(io_2xynst,1,2)%tlims = (/0,int_fill,2/)
   ENDIF

!  ---3-D normal
   modfiles(io_3uvnst,1,2)%status = 'W'
   modfiles(io_3uvnst,1,2)%form = 'N'
   modfiles(io_3uvnst,1,2)%filename = runtitle(1:8)//'_3uvobc.nc'
   modfiles(io_3uvnst,1,2)%tlims = (/0,int_fill,10/)

!  ---3-D tangential
   IF (runtitle(8:8).EQ.'A'.OR.runtitle(8:8).EQ.'C') THEN
      modfiles(io_3xynst,1,2)%status = 'W'
      modfiles(io_3xynst,1,2)%form = 'N'
      modfiles(io_3xynst,1,2)%filename = runtitle(1:8)//'_3xyobc.nc'
      modfiles(io_3xynst,1,2)%tlims = (/0,int_fill,10/)
   ENDIF

ENDIF

!
!7. Surface grid parameters
!--------------------------
!

IF (iopt_nests.EQ.1) THEN
   surfacegrids(igrd_model,1)%delxdat = 1000.0
   surfacegrids(igrd_model,1)%delydat = 1000.0
   surfacegrids(igrd_model,1)%x0dat = 0.0
   surfacegrids(igrd_model,1)%y0dat = 0.0
ELSE
   surfacegrids(igrd_model,1)%delxdat = 500.0
   surfacegrids(igrd_model,1)%delydat = 500.0
   surfacegrids(igrd_model,1)%x0dat = 50000.0
   surfacegrids(igrd_model,1)%y0dat = 25000.0
ENDIF

!
!8. User output
!--------------
!
!---time series
nosetstsr = 2
novarstsr = 11

!---harmonic output
nosetsanal = 2
nofreqsanal = 1
novarsanal = 10
nostatsanal = MERGE(9,5,iopt_nests.EQ.1)

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
! Description - test case obctang
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
INTEGER :: i, j


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!1. Open boundary locations
!--------------------------
!
!---West
iobu(1:nosbu/2) = 1
jobu(1:nosbu/2) = (/(j,j=1,nr-1)/)

!---East
iobu(nosbu/2+1:nosbu) = nc
jobu(nosbu/2+1:nosbu) = (/(j,j=1,nr-1)/)

!---South
iobv(1:nosbv/2) = (/(i,i=1,nc-1)/)
jobv(1:nosbv/2) = 1

!---North
iobv(nosbv/2+1:nosbv) = (/(i,i=1,nc-1)/)
jobv(nosbv/2+1:nosbv) = nr

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
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - test case obctang
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
USE switches
USE syspars
USE tide
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc

!
! Name      Type      Purpose
!-------------------------------------------------------------------------
!*iddesc*   INTEGER   Data file id (io_2uvobc or io_2xyobc)
!
!-------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: uflag
INTEGER :: ii, jj, l
REAL :: dphase, freq, hdel, uamp, zamp


procname(pglev+1) = 'usrdef_2dobc_spec'
CALL log_timer_in()

hdel = surfacegrids(igrd_model,1)%delxdat
freq = tidal_spectrum(icon_S2)
dphase = freq*hdel/SQRT(2.0*gacc_ref*depmean_cst)

SELECT CASE (runtitle(8:8))
   CASE ('A','B','E','F'); uflag = .FALSE.
   CASE ('C','D','G','H'); uflag = .TRUE.
END SELECT

SELECT CASE(iddesc)

!
!1. Normal open boundary conditions
!----------------------------------
!

   CASE (io_2uvobc)

!
!1.1 Type of open boundary conditions
!------------------------------------
!

      ityp2dobu = MERGE(11,12,uflag)
      ityp2dobv = MERGE(11,12,uflag)

!
!1.2 Amplitudes
!--------------
!
      
      IF (iopt_nests.EQ.1) THEN

         uamp = depmean_cst/2.0
         zamp = SQRT(depmean_cst/(2.0*gacc_ref))
         IF (uflag) THEN
            udatobu_amp(:,1) = uamp
            vdatobv_amp(:,1) = uamp
         ENDIF
         zdatobu_amp(:,1) = zamp
         zdatobv_amp(:,1) = zamp 

!
!1.3 Phases
!----------
!
!        ---West
         zdatobu_pha(1,1) = 0.0
         ii_131: DO ii=2,nobu/2
            zdatobu_pha(ii,1) = zdatobu_pha(ii-1,1) + dphase
         ENDDO ii_131
!        ---South
         zdatobv_pha(1,1) = zdatobu_pha(1,1)
         jj_132: DO jj=2,nobv/2
            zdatobv_pha(jj,1) = zdatobv_pha(jj-1,1) + dphase
         ENDDO jj_132
!        ---East
         zdatobu_pha(nobu/2+1,1) = zdatobv_pha(nobv/2,1) + dphase
         ii_133: DO ii=nobu/2+2,nobu
            zdatobu_pha(ii,1) = zdatobu_pha(ii-1,1) + dphase
         ENDDO ii_133
!        ---North
         zdatobv_pha(nobv/2+1,1) = zdatobu_pha(nobu/2,1) + dphase
         jj_134: DO jj=nobv/2+2,nobv
            zdatobv_pha(jj,1) = zdatobv_pha(jj-1,1) + dphase
         END DO jj_134

         zdatobu_pha = MOD(zdatobu_pha,twopi)
         zdatobv_pha = MOD(zdatobv_pha,twopi)

         IF (uflag) THEN
            udatobu_pha = zdatobu_pha
            vdatobv_pha = zdatobv_pha
         ENDIF

!
!1.4 Data file specifiers
!------------------------
!

      ELSE
         no2dobuv(2) = nobu+nobv
         iobc2dtype(2) = MERGE(1,2,uflag)
         index2dobuv(:,2) =  (/(l,l=1,nobu+nobv)/)
      ENDIF

!
!2. Tangential open boundary conditions
!--------------------------------------
!

   CASE (io_2xyobc)

!
!2.1 Type of open boundary conditions
!------------------------------------
!

      ityp2dobx = 2
      ityp2doby = 2

!
!2.2 Amplitudes
!--------------
!

      IF (iopt_nests.EQ.1) THEN

         vdatobx_amp(:,1) = depmean_cst/2.0
         udatoby_amp(:,1) = depmean_cst/2.0

!
!2.3 Phases
!----------
!
!        ---West
         ii_231: DO ii=1,nobx/2
            vdatobx_pha(ii,1) = zdatobu_pha(ii,1)
         ENDDO ii_231
!        ---South
         jj_232: DO jj=1,noby/2
            udatoby_pha(jj,1) = zdatobv_pha(jj,1)
         ENDDO jj_232
!        ---East
         ii_233: DO ii=1,nobx/2
            vdatobx_pha(nobx/2+ii,1) = zdatobu_pha(nobu/2+ii+1,1)
         ENDDO ii_233
!        ---North
         jj_234: DO jj=1,noby/2
            udatoby_pha(noby/2+jj,1) = zdatobv_pha(nobv/2+jj+1,1)
         ENDDO jj_234

!
!2.4 Data file specifiers
!------------------------
!

      ELSE
         no2dobxy(2) = nobx+noby
         index2dobxy(:,2) =  (/(l,l=1,nobx+noby)/)
      ENDIF

END SELECT

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
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - test case obctang
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
!*Local variables
!
INTEGER :: l


procname(pglev+1) = 'usrdef_profobc_spec'
CALL log_timer_in()

SELECT CASE (iddesc)

!  ---normal conditions
   CASE (io_3uvobc)
      iprofobux(:,1) = (/(l,l=1,nobu)/) 
      iprofobvy(:,1) = (/(l,l=nobu+1,nobu+nobv)/) 
      noprofsd(2) = nobu+nobv
      indexprof(:,2) = (/(l,l=1,nobu+nobv)/) 

!  ---tangential conditions
   CASE (io_3xyobc)
      itypobux = 2; itypobvy = 2
      iprofobux(:,1) = (/(l,l=1,nobx)/) 
      iprofobvy(:,1) = (/(l,l=nobx+1,nobx+noby)/) 
      noprofsd(2) = nobx+noby
      indexprof(:,2) = (/(l,l=1,nobx+noby)/) 
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
!*iddesc*      INTEGER Data file id
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
! Description - called by reader processes only
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
