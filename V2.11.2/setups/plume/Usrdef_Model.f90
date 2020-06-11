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
! Version - @(COHERENS)Usrdef_Model.f90  V2.8
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - test case plume
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
! Description - test case plume
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
! Version - @(COHERENS)Usrdef_Model.f90  V2.8
!
! Description - test case plume
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


procname(pglev+1) = 'usrdef_mod_params'
CALL log_timer_in()

!
!1. Setup flag
!-------------
!

sflag = MERGE(.TRUE.,.FALSE.,runtitle(6:6).EQ.'0')

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

!---implicit scheme
iopt_hydro_impl = 1

!---salinity equation (0/1)
iopt_sal = MERGE(1,2,sflag)

!---advection scheme for 2-D/3-D currents (0,1,2,3,4)
IF ((runtitle(7:7).EQ.'A').OR.(runtitle(7:7).EQ.'G')) THEN
   iopt_adv_2D = 3
   iopt_adv_3D = 3
ELSE
   iopt_adv_2D = 1
   iopt_adv_3D = 1
ENDIF
!---advection scheme for scalars (0,1,2,3,4)
iopt_adv_scal = 3

!---horizontal diffusion scheme for 2-D/3-D currents and scalars (0/1)
IF ((runtitle(7:7).EQ.'A').OR.(runtitle(7:7).EQ.'B').OR.&
  & (runtitle(7:7).EQ.'G')) THEN
   iopt_hdif_2D = 0
   iopt_hdif_3D = 0
   iopt_hdif_scal = 0
ELSE
   iopt_hdif_2D = 1
   iopt_hdif_3D = 1
   iopt_hdif_scal = 1
ENDIF
!---formulation for horizontal diffusion coefficients (0/1/2)
SELECT CASE (runtitle(7:7))
CASE ('A','B','G')
   iopt_hdif_coef = 0
CASE ('C')
   iopt_hdif_coef = 2
CASE DEFAULT
   iopt_hdif_coef = 1
END SELECT

!---type of background mixing scheme (0/1/2)
IF (.NOT.sflag) iopt_turb_iwlim = 1

!---salinity o.b.c (0/1)
iopt_obc_sal = 1
!---2-D mode o.b.c (0/1)
iopt_obc_2D = 1
!---3-D mode o.b.c (0/1)
iopt_obc_3D = 1

!---implicit scheme
iopt_mg_smoother = 2

!---user output (0/1)
IF (sflag) THEN
   iopt_out_tsers = 0
ELSE
   iopt_out_anal = 1
ENDIF

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
delt2d = MERGE(30.0,300.0,iopt_hydro_impl.EQ.0)

!---counter for 3-D mode
ic3d = MERGE(10,1,iopt_hydro_impl.EQ.0)

!
!5. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 121; nr = 41; nz = 20

!---number of open sea boundaries
nosbu = 80; nosbv = 120

!---number of river open boundaries
nrvbv = 1

!---number of tidal constituents
IF (runtitle(7:7).NE.'G') nconobc = 1

!---reference latitude
dlat_ref = 52.0

!---uniform mean water depth
depmean_cst = 20.0

!---uniform horizontal diffusion coefficient for momentum [m^2/s]
IF (iopt_hdif_coef.EQ.1) hdifmom_cst = 100.0
!---uniform horizontal diffusion coefficient for scalars [m^2/s]
IF (iopt_hdif_coef.EQ.1) hdifscal_cst = 100.0

!---uniform bottom roughness length
zrough_cst = 0.006

!---multigrid parameters
nomglevels = 3
mg_tol = 1.0E-07

!---tidal indices
index_obc(1) = icon_M2

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
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = 'U'
   modfiles(io_inicon,ics_phys,1)%filename = &
                                & runtitle(1:5)//runtitle(7:7)//'.phsfin.bin'
ENDIF

!---open boundary conditions (2-D)
modfiles(io_2uvobc,1:2,1)%status = 'N'
!---open boundary conditions (3-D)
modfiles(io_3uvobc,1:2,1)%status = 'N'
!---open boundary conditions (salinity)
modfiles(io_salobc,1:2,1)%status = 'N'

!
!6.2 Output
!----------
!
!---model grid
modfiles(io_modgrd,1,2)%status = 'W'
modfiles(io_modgrd,1,2)%form = 'A'
modfiles(io_modgrd,1,2)%filename = &
                        & runtitle(1:5)//runtitle(7:7)//'.modgrd.txt'

!---initial conditions
IF (sflag) THEN
   modfiles(io_fincon,ics_phys,2)%status = 'W'
   modfiles(io_fincon,ics_phys,2)%form = 'U'
   modfiles(io_fincon,ics_phys,2)%filename = &
                               & runtitle(1:5)//runtitle(7:7)//'.phsfin.bin'
ELSE
   modfiles(io_fincon,1,2)%status = '0'
ENDIF

!---restart times
norestarts = 1
ntrestart(1) = MERGE(int_fill,1,sflag)

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
!---title for user-output files
outtitle = runtitle(1:5)//runtitle(7:7)
!---time series output
IF (iopt_verif.EQ.0) THEN
   nosetstsr = 4
   novarstsr = 9
ELSE
   nosetstsr = 1
   novarstsr = 10
ENDIF
!---harmonic analysis
nosetsanal = 1
nofreqsanal = 1
novarsanal = 7

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
! Description - test case plume
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
!2. Open boundary locations
!--------------------------
!
!---U-nodes
iobu(1:40) = 1; iobu(41:80) = nc
jobu(1:40) = (/(j,j=1,40)/); jobu(41:80) = (/(j,j=1,40)/) 

!---V-nodes
iobv(1:120) = (/(i,i=1,120)/); jobv(1:120) = nr
iobv(121) = 26; jobv(121) = 1

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
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.1
!
! Description - test case plume
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
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
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - test case plume
!
! Reference -
!
! Calling program - define_2dobc_spec
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE obconds
USE physpars
USE syspars
USE tide
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
INTEGER :: ii, j
REAL :: crad0, dhun, freq, rl, yrel, zeta0


procname(pglev+1) = 'usrdef_2dobc_spec'
CALL log_timer_in()

!---tidal frequencies
IF (runtitle(7:7).NE.'G') THEN
   freq = tidal_spectrum(icon_M2)
ENDIF

!---type of conditions at U-open boundaries
IF (runtitle(7:7).EQ.'E') THEN
   ityp2dobu(1:nobu/2) = 12
ELSEIF (runtitle(7:7).EQ.'F') THEN
   ityp2dobu(1:nobu/2) = 4
ELSEIF (runtitle(7:7).EQ.'G') THEN
   ityp2dobu(1:nobu/2) = 13
ELSE
   ityp2dobu(1:nobu/2) = 11
ENDIF
ityp2dobu(nobu/2+1:nobu) = 13
iloczobu = 1

!---type of conditions at V-open boundaries
ityp2dobv(1:nobv-1) = 13
ityp2dobv(nobv) = 4
iloczobv = 1

!---contents of data file
no2dobuv(2) = 1
iobc2dtype(2) = 3
index2dobuv(1,2) = nobu+nobv

!---amplitudes/phases at U-nodes
IF (runtitle(7:7).NE.'G') THEN
   dhun = 1000.0; zeta0 = 0.8
   ii_130: DO ii=1,nobu/2
      j = jobu(ii)
      crad0 = SQRT(gacc_mean*depmean_cst)
      rl = crad0/coriolatu(1,1)
      yrel = dhun*(j-0.5)/rl
      zdatobu_amp(ii,1) = zeta0*EXP(-yrel)
      udatobu_amp(ii,1) = crad0*zeta0*EXP(-yrel)
      udatobu_pha(ii,1) = halfpi
      zdatobu_pha(ii,1) = halfpi
   ENDDO ii_130
ENDIF

!---amplitudes/phases at V-nodes
IF (runtitle(7:7).NE.'G') THEN
   vdatobv_amp(nobv,1) = 0.6*depmean_cst
   crad0 = SQRT(gacc_mean*depmean_cst)
   vdatobv_pha(nobv,1) = freq*25500.0/crad0
ENDIF

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
! Description - test case plume
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
   CASE (io_3uvobc,io_salobc)
      noprofsd = 1
      iprofobvy(nobv,1) = 1
      indexprof(1,2) = 1
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
! Description - test case plume
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
!*Local variables
!
REAL :: q, wr


procname(pglev+1) = 'usrdef_2dobc_data'
CALL log_timer_in()


IF (modfiles(iddesc,ifil,1)%iostat.EQ.0) THEN 
   modfiles(iddesc,ifil,1)%iostat = 1
   GOTO 1000
ENDIF

ciodatetime = CStartDateTime

q = 1500.0; wr = 1000.0
data2d(1,1) = q/wr

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
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.2
!
! Description - test case plume
!
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
!*Local variables
!
INTEGER :: k
REAL :: gzc


procname(pglev+1) = 'usrdef_profobc_data'
CALL log_timer_in()


IF (modfiles(iddesc,ifil,1)%iostat.EQ.0) THEN
   modfiles(iddesc,ifil,1)%iostat = 1
   GOTO 1000
ENDIF

ciodatetime = CStartDateTime

SELECT CASE (iddesc)

CASE (io_3uvobc)
   k_110: DO k=1,nz
      gzc = (k-0.5)/REAL(nz)
      IF (gzc.GT.0.75) THEN
         psiprofdat(1,k) = 0.6
      ELSE
         psiprofdat(1,k) = -0.2
      ENDIF
   ENDDO k_110

CASE (io_salobc)
   IF (runtitle(6:6).EQ.'0') THEN
      psiprofdat = float_fill
   ELSE
      k_120: DO k=1,nz
         gzc = (k-0.5)/REAL(nz)
         IF (gzc.GT.0.75) THEN
            psiprofdat(1,k) = 10.0
         ELSE
            psiprofdat(1,k) = float_fill
         ENDIF
      ENDDO k_120
   ENDIF

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
