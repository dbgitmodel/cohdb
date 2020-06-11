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
! *Usrdef_Model* User-defined model setup (example routines)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - test case csnsp
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
! *usrdef_init_params* Define file statuses for model I/O
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.1
!
! Description - test cse csnsp
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


RETURN

END SUBROUTINE usrdef_init_params

!========================================================================

SUBROUTINE usrdef_mod_params
!************************************************************************
!
! *usrdef_mod_params* Define parameters for control file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.1
!
! Description - test case csnsp
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
USE tide
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) ='usrdef_mod_params'
CALL log_timer_in()

!
!3. Switches
!-----------
!
!---grid dimension (1/2/3)
iopt_grid_nodim = 1

!---equation of state (0/1/2)
iopt_dens = 2

!---temperature equation (0/1)
iopt_temp = 2
!---optical module (0/1)
iopt_temp_optic = MERGE(0,1,runtitle(6:6).EQ.'I')

!---advection scheme for 2-D currents (0/1/2/3/4)
iopt_adv_2D = 0
!---advection scheme for 3-D currents (0/1/2/3/4)
iopt_adv_3D = 0
!---advection scheme for scalars (0/1/2/3/4)
iopt_adv_scal = 0

!---type of vertical diffusion scheme (0/1/2/3)
iopt_vdif_coef = MERGE(2,3,runtitle(6:6).EQ.'D')
!---type of algebraic turbulence closure (1/2/3/4/5)
IF (iopt_vdif_coef.EQ.2) iopt_turb_alg = 2
!---type of background mixing scheme (0/1/2)
IF ((runtitle(6:6).EQ.'B').OR.(runtitle(6:6).EQ.'C')) THEN
   iopt_turb_iwlim = 0
ELSE
   iopt_turb_iwlim = 1
ENDIF
!---form of stability functions (1/2/3)
iopt_turb_stab_form = MERGE(2,3,runtitle(6:6).EQ.'C')

!---type of implicity for Coriolis terms (0/1/2)
iopt_cor_impl = 2

!---type of surface forcing in 1-D case (0/1/2)
iopt_sur_1D = 1

!---meteo input (0/1)
iopt_meteo = 1

!---surface exchange coefficients
SELECT CASE (runtitle(6:6))
   CASE ('E'); iopt_sflux_pars = 5
   CASE ('G'); iopt_sflux_pars = 6
   CASE ('H'); iopt_sflux_pars = 1
   CASE DEFAULT; iopt_sflux_pars = 7
END SELECT
      
!---stratification effects for exchange coefficients
iopt_sflux_strat = 1
IF ((runtitle(6:6).EQ.'F').OR.(runtitle(6:6).EQ.'H')) THEN
   iopt_sflux_strat = 0
ELSE
   iopt_sflux_strat = 1
ENDIF

!---time-averaged output (0/1)
iopt_out_avrgd = 1

!
!4. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
CStartDateTime(1:19) = '1989/01/03;00:00:00'
CEndDateTime(1:19) = '1989/10/31;21:00:00'

!---time step (milliseconds)
delt2d = 1200.0

!
!5. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 3; nr = 3; nz = 50
nobu = 4; nobv = 4

!---number of tidal constituents
nconobc = 1

!---reference density [kg/m^3]
density_ref = 1026.8
!---reference latitude
dlat_ref = 55.515
!---reference longitude
dlon_ref = 0.9083333
!---reference salinity [PSU]
sal_ref = 34.8
!---reference temperature [deg C]
temp_ref = 10.0

!---uniform mean water depth
depmean_cst = 83.0

!---uniform bottom roughness length
zrough_cst = 0.0024

!---implicity factor for Coriolis terms (between 0 and 1)
theta_cor = 1.0

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
modfiles(io_inicon,ics_phys,1)%status = 'N'
modfiles(io_inicon,ics_phys,1)%filename = 'tin.cs'
modfiles(io_inicon,ics_phys,1)%form = 'A'
!---surface forcing
modfiles(io_1uvsur,1,1)%status = 'N'
!---meteo
modfiles(io_metsur,1,1)%status = 'N'
modfiles(io_metsur,1,1)%filename = 'metin.cs'
modfiles(io_metsur,1,1)%form = 'A'
modfiles(io_metsur,1,1)%tlims = (/0,int_fill,9/)

!
!6.2 Output
!----------
!

norestarts = 0

IF (ciffile%status.EQ.'W') THEN
!  ---restart times
   norestarts = 1
   ntrestart(1) = 0
!  ---initial conditions
   modfiles(io_fincon,ics_phys,2)%status = 'W'
   modfiles(io_fincon,ics_phys,2)%form = 'U'
!  ---surface forcing
   modfiles(io_1uvsur,1,2)%status = 'W'
   modfiles(io_1uvsur,1,2)%form = 'U'
ENDIF

!
!8. User output
!--------------
!

IF (iopt_verif.EQ.0) THEN
   nosetstsr = 2
   nosetsavr = 1
   novarstsr = 2
ELSE
   nosetstsr = 1
   nosetsavr = 1
   novarstsr = 7
ENDIF
novarsavr = 4

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_mod_params

!========================================================================

SUBROUTINE usrdef_grid
!************************************************************************
!
! *usrdef_grid* Define arrays for grid file
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

END SUBROUTINE usrdef_grid

!========================================================================

SUBROUTINE usrdef_partition
!************************************************************************
!
! *usrdef_partition* Define domain decomposition for parallel mode
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
! *usrdef_phsics* Define arrays for initial conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - test case csnsp
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - close_filepars, error_abort, open_filepars
!
!************************************************************************
!
USE density
USE gridpars
USE iopars
USE error_routines, ONLY: nerrs, error_abort
USE inout_routines, ONLY: close_filepars, open_filepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iunit, k, nnz
REAL, DIMENSION(nz) :: tin


procname(pglev+1) = 'usrdef_phsics'
CALL log_timer_in()

!
!3.Temperature
!-------------
!
!---open file with initial profile
CALL open_filepars(modfiles(io_inicon,ics_phys,1))
iunit = modfiles(io_inicon,ics_phys,1)%iunit

!---check vertical dimension
READ (iunit,*) nnz
IF (nnz.NE.nz) THEN
   nerrs = 1
   WRITE (0,'(A)')    'Invalid value for parameter nz'
   WRITE (0,'(A)')    '  must be the same as in file tin.cs'
   WRITE (0,'(A,I3)') '  Input value from file tin.cs : ',nnz
   WRITE (0,'(A,I3)') '  Present value                : ',nz
   CALL error_abort('usrdef_phsics',ierrno_input)
ENDIF

!---read initial conditions
READ (iunit,*) tin
k_310: DO k=1,nz
   temp(1:nc-1,1:nr-1,k) = tin(k)
ENDDO k_310

!---close file
CALL close_filepars(modfiles(io_inicon,ics_phys,1))

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
! Description - test case csnsp
!
! Reference -
!
! Calling program - define_1dsur_spec
!
!************************************************************************
!
USE iopars
USE obconds
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_1dsur_spec'
CALL log_timer_in()

!---amplitudes and phases of X-component of pressure gradient
gxslope_amp(1) = 4.2227563E-05
gxslope_pha(1) = 248.2057*degtorad

!---amplitudes and phases of Y-component of pressure gradient
gyslope_amp(1) = 4.6298766E-05
gyslope_pha(1) = 335.6429*degtorad

!---type of data in eventual data file (1/2/3)
isur1dtype = 3

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_1dsur_spec

!========================================================================

SUBROUTINE usrdef_2dobc_spec(iddesc)
!************************************************************************
!
! *usrdef_2dobc_spec* Define specifier arrays for 2-D open boundary conditions
!                     at U- and V-nodes
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
! *usrdef_2dobc_data* Define open boundary data for 2-D mode at U- and V-nodes
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
!*iddesc*    INTEGER Data file id (io_2uvobc or io_2xyobc)
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
! *usrdef_rlxobc_spec* Define specifier arrays for relaxation scheme
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
