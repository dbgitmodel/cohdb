77!
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
! $Date: 2014-09-19 18:01:49 +0200 (Fri, 19 Sep 2014) $
!
! $Revision: 741 $
!
! Description - 
!
! Reference -
!
! Routines - usrdef_init_params, usrdef_mod_params, usrdef_grid,
!            usrdef_partition, usrdef_phsics, usrdef_1dsur_spec,
!            usrdef_2dobc_spec, usrdef_physobc_spec, usrdef_1dsur_data,
!            usrdef_2dobc_data, usrdef_physobc_data, usrdef_rlxobc_spec
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_init_params
!************************************************************************
!
! *usrdef_init_params* Define parameters for monitoring (example routine)
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
! Module calls -
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

cold_start = .FALSE.

!
!2. Log files
!------------
!
!---program leveling in log files
iproc_210: DO iproc=1,npworld
   levprocs_ini(iproc) = 0
   levprocs_run(iproc) = 0
ENDDO iproc_210

!---(generic) names of log files
inilog_file = TRIM(runtitle)//'.inilogA'
runlog_file = TRIM(runtitle)//'.runlogA'

!---log parameters
exitlog = .TRUE.
runlog_count = int_fill

!
!3. Error coding
!---------------
!
!---maximum number of error messages
maxerrors = MaxErrMesgs

!---error trapping level
iproc_310: DO iproc=1,npworld
   levprocs_err(iproc) = 1
ENDDO iproc_310

!---(generic) name of error file
errlog_file = TRIM(runtitle)//'.errlogA'

!
!4. Warnings
!-----------
!

warning = .TRUE.
warlog_file = TRIM(runtitle)//'.warlogA'

!
!5. Timing
!---------
!
!---timer level (0/1/2/3)
levtimer = 0

!---name of timer report file
timing_file = TRIM(runtitle)//'.timingA'

!---time format (1/2/3/4)
timer_format = 1


RETURN

END SUBROUTINE usrdef_init_params

!========================================================================

SUBROUTINE usrdef_mod_params
!************************************************************************
!
! *usrdef_mod_params*  Define parameters for physical model (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.5.1
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - 
!
! Module calls -
!
!************************************************************************
!
USE gridpars
USE iopars
USE nestgrids
USE paralpars
USE physpars
USE relaxation
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
INTEGER :: iddesc, ifil, iset, nhtype


procname(pglev+1) = 'usrdef_mod_params'
CALL log_timer_in()

!
!1. Process numbers
!------------------
!

IF (npworld.GT.1) THEN
!  ---total number of processes
   nprocs = 1
!  ---number of processes in X-direction
   nprocsx = 0
!  ---number of processes in Y-direction
   nprocsy = 0
ENDIF

!
!2. Switches
!-----------
!
!2.1 Model grid
!--------------
!
!---horizontal grid type (1/2/3)
iopt_grid_htype = 1
!---grid dimension (1/2/3)
iopt_grid_nodim = 3
!---Cartesian or spherical (0/1)
iopt_grid_sph = 0
!---vertical grid type (1/2/3)
iopt_grid_vtype = 1
!---vertical grid transform(0/11/12/13/21)
iopt_grid_vtype_transf = 0

!
!2.2 Interpolation
!-----------------
!
!---horizontal interpolation (0/1)
iopt_arrint_hreg = 0
!---vertical interpolation (0/1)
iopt_arrint_vreg = 0
!---dimension of weight factors (0/1)
iopt_arrint_3D = 0

!
!2.3 Hydrodynamics
!-----------------
!
!---currents (0/1/2)
iopt_curr = 2
!---settling velocity (1/2)
iopt_curr_wfall = 1

!
!2.4 Density
!-----------
!
!---equation of state (0/1/2/3)
iopt_dens = 0
!---formulation for baroclinic pressure gradient (0/1/2/3)
iopt_dens_grad = 1
!---salinity equation (0/1/2)
iopt_sal = 0
!---salinity boundary condition (0/1)
iopt_sal_sbc = 0
!---temperature equation (0/1/2)
iopt_temp = 0
!---optical module (0/1)
iopt_temp_optic = 1
!---temperature surface boundary condition (1/2/3)
iopt_temp_sbc = 1

!
!2.5 Biology/sediments
!---------------------
!
!---disable/enable biology (0/1)
iopt_biolgy = 0
!---disable/enable sediments (0/1)
iopt_sed = 0

!
!2.6 Bottom boundary conditions
!------------------------------
!
!---formulation for bottom drag coefficient (0/1/2/3/4)
iopt_bstres_drag = 3
!---bottom stress formulation (0/1/2)
iopt_bstres_form = 2
!---dimension of bottom stress formulation (2/3)
iopt_bstres_nodim = 3

!
!2.7 Transport
!------------
!
!---3D/4D (0/1)
iopt_transp_full = 0

!
!2.8 Advection
!-------------
!
!---advection scheme for scalars (0/1/2/3)
iopt_adv_scal = 3
!---advection scheme for turbulence variables (0/1/2/3)
iopt_adv_turb = 0
!---type limiter for TVD scheme (1/2)
iopt_adv_tvd = 1
!---advection scheme for 2-D currents (0/1/2/3)
iopt_adv_2D = 1
!---advection scheme for 3-D currents (0/1/2/3)
iopt_adv_3D = 1

!
!2.9 Diffusion
!-------------
!
!---formulation for horizontal diffusion coefficients (0/1/2)
iopt_hdif_coef = 0
!---horizontal diffusion scheme for scalars (0/1)
iopt_hdif_scal = 0
!---horizontal diffusion scheme for turbulence variables (0/1)
iopt_hdif_turb = 0
!---horizontal diffusion scheme for 2-D currents (0/1)
iopt_hdif_2D = 0
!---horizontal diffusion scheme for 3-D currents (0/1)
iopt_hdif_3D = 0
!---type of kinematic viscosity (0/1)
iopt_kinvisc = 0
!---type of vertical diffusion scheme (0/1/2/3)
iopt_vdif_coef = 3

!
!2.10 Turbulence
!---------------
!
!---type of algebraic turbulence closure (1/2/3/4/5/6)
iopt_turb_alg = 1
!---type bottom b.c. for dissipation equation (1/2)
iopt_turb_dis_bbc = 2
!---type surface b.c. for dissipation equation (1/2)
iopt_turb_dis_sbc = 2
!---type of internal wave mixing (0/1/2)
iopt_turb_iwlim = 0
!---type of background mixing (0/1)
iopt_turb_kinvisc = 0
!---type of algebraic mixing length formulation (1/2/3/4)
iopt_turb_lmix = 4
!---number of transport equations for turbulence (0/1/2)
iopt_turb_ntrans = 1
!---k-l or k-eps type of turbulence scheme (1/2)
iopt_turb_param = 2
!---form of stability functions (1/2/3)
iopt_turb_stab_form = 3
!---quasi-equilibrium or non-equilibrium scheme (1/2)
iopt_turb_stab_lev = 1
!---model for stability functions if iopt_turb_stab_form=3 (1/2/3/4/5/6)
iopt_turb_stab_mod = 4
!---form of stability function for turbulence energy (1/2/3)
iopt_turb_stab_tke = 2
!---type bottom b.c. for turbulence energy (1/2)
iopt_turb_tke_bbc = 2
!---type surface b.c. for turbulence energy (1/2)
iopt_turb_tke_sbc = 2

!
!2.11 Flooding/drying
!--------------------
!
!---flooding/drying scheme (0/1/2)
iopt_fld = 0

!
!2.12 Structures
!---------------
!
!---dry cells (0/1)
iopt_drycel = 0
!---thin dams (0/1)
iopt_thndam = 0
!---weirs/barriers (0/1)
iopt_weibar = 0
!---discharges (0/1)
iopt_dischr = 0

!
!2.13 Time integration
!---------------------
!
!---type of implicity for Coriolis terms (0/1/2)
iopt_cor_impl = 1
!---explicit/implicit scheme (0/1)
iopt_hydro_simpl = 0
!---deposition flux formulation (0/1/2)
iopt_scal_depos = 1
!---type of implicity for vertical advection (0/1/2)
iopt_vadv_impl = 1
!---type of implicity for vertical diffusion (0/1/2)
iopt_vdif_impl = 2

!
!2.14 Open boundary conditions
!-----------------------------
!
!---open boundary condition for advective fluxes (1/2)
iopt_obc_advflux = 1
!---relaxation open boundary condition for momentum (0/1)
iopt_obc_advrlx = 0
!---biology o.b.c (0/1)
iopt_obc_bio = 0
!---inverse barometric effect (0/1)
iopt_obc_invbar = 0
!---open boundary relaxation (0/1)
iopt_obc_relax = 0
!---salinity o.b.c (0/1)
iopt_obc_sal = 0
!---sediment o.b.c (0/1)
iopt_obc_sed = 0
!---temperature o.b.c (0/1)
iopt_obc_temp = 0
!---2-D mode normal o.b.c (0/1)
iopt_obc_2D = 0
!---3-D mode normal o.b.c (0/1)
iopt_obc_3D = 0
!---2-D mode tangential o.b.c (0/1)
iopt_obc_2D_tang = 0
!---3-D mode tangential o.b.c (0/1)
iopt_obc_3D_tang = 0
!---thatcher-Harleman (0/1)
iopt_obc_th = 0

!
!2.15 Tides
!----------
!
!---astronomical phases for harmonic analysis (0/1)
iopt_astro_anal = 0
!---nodal factors and astronomical argument (0/1/2)
iopt_astro_pars = 0
!---astronomical tidal force (0/1)
iopt_astro_tide = 0

!
!2.16 1-D applications
!---------------------
!
!---1-D applications (0/1)
iopt_sur_1D = 0

!
!2.17 Meteo
!----------
!
!---meteo input (0/1)
iopt_meteo = 0
!---type of meteo input for heat fluxes (0/1/2/3/4/5)
iopt_meteo_heat = 0
!---type of meteo input for salinity fluxes (0/1/2)
iopt_meteo_salflx = 0
!---type of meteo input for surface stress and pressure (0/1/2)
iopt_meteo_stres = 0

!
!2.18 Surface boundary conditions
!--------------------------------
!
!---formulation for surface drag coefficient (0/1/2/3/4/5/6)
iopt_sflux_cds = 0
!---formulation for surface exchange coefficients (0/1/2/3/4)
iopt_sflux_cehs = 0
!---stratification effects for surface drag and exchange coefficients (0/1/2)
iopt_sflux_strat = 0

!
!2.19 Nesting
!------------
!
!---nested grids (0/1)
iopt_nests = 0

!
!2.19 Waves
!----------
!
!---wave module (0/1)
iopt_waves = 0

!
!2.20 Parallel mode (MPI)
!------------------------
!
!---MPI switches for parallel processing
IF (parallel_set) THEN
!  ---error coding (0/1)
   iopt_MPI_abort = 0
!  --all-to-all communication (1/2/3/4)
   iopt_MPI_comm_all = 2
!  ---collective calls (0/1)
   iopt_MPI_comm_coll = 0
!  ----exchange communications (1/2/3/4/5)
   iopt_MPI_comm_exch = 2
!  ---exchange of "full" arrays (0/1)
   iopt_MPI_comm_full = 0
!  ---combine communications (1/2/3/4)
   iopt_MPI_comm_gath = 2
!  ---distribute communications (1/2/3/4)
   iopt_MPI_comm_scat = 2
!  --type of domain decomposition (1/2)
   iopt_MPI_partit = 1
!  ---synchronisation call (0/1)
   iopt_MPI_sync = 0
ENDIF

!
!2.21 PetSC library options
!--------------------------
!
!---type of 
IF (iopt_petsc.EQ.1) THEN
!  ---type of preconditioner (1-12)
   iopt_petsc_precond = 5
!  ---type of solver (1-12)
   iopt_petsc_solver = 5
ENDIF

!
!2.22 User output
!----------------
!
!---harmonic analysis (0/1)
iopt_out_anal = 0
!---time-averaged output (0/1)
iopt_out_avrgd = 0
!---time series output (0/1)
iopt_out_tsers = 1

!
!2.23 netcdf
!-----------
!
!---error coding (0/1)
iopt_CDF_abort = 0
!---fill mode (0/1)
iopt_CDF_fill = 0
!---file format (0/1)
iopt_CDF_format = 0

!
!3. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
CStartDateTime(1:23) = 'xxxx/xx/xx;00:00:00,000'
CEndDateTime(1:23) = 'xxxx/xx/xx;00:00:00,000'

!---time step (seconds)
delt2d = ?

!---counter for 3-D mode
ic3d = 1

!---counter for update nodal factors and tidal phases
icnodal = 0

!---time zone
time_zone = 0.0

!
!4. Physical model constants
!---------------------------
!
!---model grid
IF (iopt_grid_nodim.GT.1) THEN
   nc = ?; nr = ?
ENDIF
IF (iopt_grid_nodim.NE.2) nz = ?

!---number of open sea boundaries
IF (iopt_grid_nodim.GT.1) THEN
   nosbu = 0; nosbv = 0
ENDIF

!---number of river open boundaries
IF (iopt_grid_nodim.GT.1) THEN
   nrvbu = 0; nrvbv = 0
ENDIF

!---number of relaxation zones
IF (iopt_obc_relax.EQ.1) norlxzones = 0

!---number of tidal constituents
IF (iopt_astro_tide.EQ.1) nconastro = 0
IF (iopt_obc_2D.EQ.1) nconobc = 0

!---number of nested grids
IF (iopt_nests.EQ.1) nonestsets = 0

!---wait time
nowaitsecs = 0
maxwaitsecs = 3600

!---structures
!---dry cells
IF (iopt_drycel.EQ.1) numdry = ?
!---thin dams
IF (iopt_thndam.EQ.1) THEN
   numthinu = 0
   numthinv = 0
ENDIF
!---weirs/barriers
IF (iopt_weibar.EQ.1) THEN
   numwbaru = 0
   numwbarv = 0
ENDIF

!---discharges
IF (iopt_dischr.EQ.1) numdis = ?

!---unit record size (in bytes)
nrecunit = 4

!---restart times
norestarts = 1
ntrestart(1) = int_fill

!---relaxation time for 2-D open boundary conditions
ntobcrlx = 0

!---return tiles (TH obc condition)
return_time = ?

!---master process id
idmaster = 0

!---Earth's radius [m]
Rearth = 6371000.0
!---air density [kg/m^3]
rho_air = 1.2
!---specific heat [J/kg/deg C]
specheat = 3987.5

!---reference atmospheric pressure [Pa]
atmpres_ref = 101325.0
!---reference latitude [degrees]
If (iopt_grid_sph.EQ.0) dlat_ref = 0.0
!---reference longitude [degrees]
IF (iopt_grid_sph.EQ.0) dlon_ref = 0.0
!---reference longitudes for astronomical phases [degrees]
IF (iopt_astro_anal.EQ.1) dlon_ref_anal = 0.0
IF (iopt_astro_pars.GT.0) dlon_ref_obc = 0.0
!---acceleration of gravity [m/s^2]
gacc_ref = real_fill
!---reference salinity [PSU]
sal_ref = 33.0
!---reference sst [degC]
sst_ref = 12.0
!---minimum temperature [degC]
IF (iopt_temp.GT.0) temp_min = 0.0
!---reference temperature [deg C]
temp_ref = 12.0

!---model grid
b_SH = 0.1; dl_BB = 1.5; du_BB = 1.5; hcrit_SH = 200.0; sigstar_DJ = 0.0
sig0_DJ = 0.1; theta_SH = 8.0

!---uniform mean water depth [m]
depmean_cst = 0.0

!---mean depth flag
depmean_flag = 0.0 

!---bottom to roughness height limit
zbtoz0lim = 2.0

!---drying/flooding depth parameters [m]
IF (iopt_fld.GT.0) THEN
   dcrit_fld = 0.1
   dmin_fld = 0.02
ENDIF
IF (iopt_fld.EQ.2) THEN
   dthd_fld = 0.1
   fld_mask(1) = 1; fld_mask(2:nofldmasks) = 0
ENDIF

!---uniform horizontal diffusion coefficient for momentum [m^2/s]
IF (iopt_hdif_coef.EQ.1) hdifmom_cst = 0.0
!---uniform horizontal diffusion coefficient for scalars [m^2/s]
IF (iopt_hdif_coef.EQ.1) hdifscal_cst = 0.0
!---uniform kinematic viscosity [m^2/s]
kinvisc_cst = 1.0E-06
!---Smagorinsky coefficient for momentum
IF (iopt_hdif_coef.EQ.2) smag_coef_mom = 0.1
!---Smagorinsky coefficient for scalars
IF (iopt_hdif_coef.EQ.2) smag_coef_scal = 0.1
!---uniform/molecular vertical diffusion coefficient for momentum [m^2/s]
IF (iopt_vdif_coef.GT.0) vdifmom_cst = 1.0E-06
!---uniform/molecular vertical diffusion coefficient for scalars [m^2/s]
IF (iopt_vdif_coef.GT.0) vdifscal_cst = 1.0E-06

!---bottom linear drag coefficient [m/s]
IF (iopt_bstres_form.EQ.1) bdraglin = 0.0
!---uniform bottom drag coefficient
IF (iopt_bstres_drag.EQ.1) bdragcoef_cst = 0.0
!---uniform bottom roughness length [m]
IF (iopt_bstres_drag.EQ.3) zrough_cst = 0.0
!---uniform surface drag coefficient
IF (iopt_meteo.EQ.1.AND.iopt_sflux_cds.EQ.0) cds_cst = 0.0013
!---uniform surface exchange coefficient for latent heat
IF (iopt_meteo.EQ.1.AND.iopt_sflux_ces.EQ.0) ces_cst = 0.0013 
!---uniform surface exchange coefficient for sensible heat
IF (iopt_meteo.EQ.1.AND.iopt_sflux_ces.EQ.0) chs_cst = 0.0013
!---Charnock's constant
ccharno = 0.014
!---von Karman's constant
ckar = 0.4
!---atmospheric reference height [m]
zref_atm  = 10.0

!---parameters for exchange coefficients in tabular form
IF (iopt_sflux_strat.EQ.2) THEN
!  ---relative humidity [-]
   drelhum = 0.05; relhummax = 1.0; relhummin = 0.5
!  ---air temperature [degC]
   dtempdif = 1.0; dtempmax = 5.0; dtempmin = -5.0
!  ---surface wind [m/s]
   dwind = 0.25; uwindmax = 50.0; uwindmin = 3.0 
ENDIF

!---optical coefficients
IF (iopt_temp_optic.EQ.1) THEN
!  ---attenuation coefficient for long waves [1/m]
   optattcoef1_cst = 10.0
!   ---attenuation coefficient for short waves [1/m]
   optattcoef2_cst = 2.06
!  ---infrared fraction of irradiance absorbed at sea surface
   opt_frac = 0.54
ENDIF

!---relaxation factor at weirs/barriers
IF (iopt_weibar.EQ.1) THEN
   wbarrlxu = 1.0
   wbarrlxv = 1.0
ENDIF

!---implicity factor for Coriolis terms (between 0 and 1)
theta_cor = 0.5
!---implicity factor for vertical advection (between 0 and 1)
theta_vadv = 0.501
!---implicity factor for vertical diffusion (between 0 and 1)
theta_vdif = 1.0

!---parameters for implicit scheme and solver
maxitsimp = 1
dzetaresid_conv = 1.0E-14
petsc_tol = 1.0E-07

!---ratio of internal/external wave speed (o.b.c.)
cgravratio = 0.03

!---relaxation distance for momentum advection
distrlx_obc = 0.0

!---tidal index IDs
IF (nconobc.GT.0) index_obc(1:nconobc) = ?
IF (nconastro.GT.0) index_astro(1:nconastro) = ?

!
!5. Turbulence model constants
!-----------------------------
!
!---numerical lower limits
dissipmin = 1.0E-012
tkemin = 1.0E-014
zlmixmin = 1.7E-010

!---diffusion term in t.k.e. equation
c_sk = 0.15
sigma_k = 1.0

!---roughness lengths
zrough_bot = 0.0
zrough_sur = 0.0

!---surface t.k.e. flux parameter
wfltke = 0.0

!---Mellor_Yamada parameters
e1_my = 1.8
e2_my = 1.33
e3_my = 1.0
sq_my = 0.2

!---k-eps model parameters
c1_eps = 1.44
c2_eps = 1.92
c31_eps = 0.2
c32_eps = 1.0
skeps = 0.09

!---mixing length formulations
alpha_Black = 0.2
beta_Xing = 2.0

!---background mixing schemes
tkelim = 1.0E-06
riccrit_iw = 0.7
vdifmom_iw = 1.0E-04
vdifscal_iw = 5.0E-05
vdifshear_iw = 0.005

!---Pacanowski-Philander
alpha_pp = 5.0
expmom_pp= 2.0
vbmom_pp = 1.0E-04
vbscal_pp = 1.0E-05
vmax_pp = 3.0
v0dif_pp = 0.01

!---Munk_anderson
alpha_ma = 10.0
beta_ma = 3.33
expmom_ma = 0.5
expscal_ma = 1.5
vmaxmom_ma = 3.0
vmaxscal_ma = 4.0
v0dif_ma = 0.06

!---flow-dependent schemes
cnu_ad = 2.0
delta1_ad = 0.0
delta2_ad = 0.0
k1_ad = 0.0025
k2_ad = 2.0E-05
lambda_ad = 0.0
omega1_ad = 1.0E-04
r1_ad = 1.0
r2_ad = 1.0

!
!6. Surface grid parameters
!--------------------------
!
!---model grid
IF (iopt_grid_nodim.GT.1.AND.iopt_grid_htype.EQ.1) THEN
   surfacegrids(igrd_model,1)%rotated = .FALSE.
   surfacegrids(igrd_model,1)%delxdat = 1.0
   surfacegrids(igrd_model,1)%delydat = 1.0
   surfacegrids(igrd_model,1)%gridangle = 0.0
   surfacegrids(igrd_model,1)%longpole = 0.0
   surfacegrids(igrd_model,1)%x0dat = 0.0
   surfacegrids(igrd_model,1)%y0dat = 0.0
ENDIF

!---meteo
IF (iopt_meteo.EQ.1) THEN
   surfacegrids(igrd_meteo,1)%nhtype = ?
   SELECT CASE (surfacegrids(igrd_meteo,1)%nhtype)
      CASE (1)
         surfacegrids(igrd_meteo,1)%n1dat = ?
         surfacegrids(igrd_meteo,1)%n2dat = ?
         surfacegrids(igrd_meteo,1)%delxdat = ?
         surfacegrids(igrd_meteo,1)%delydat = ?
         surfacegrids(igrd_meteo,1)%x0dat = ?
         surfacegrids(igrd_meteo,1)%y0dat = ?
      CASE (2:3)
         surfacegrids(igrd_meteo,1)%n1dat = ?
         surfacegrids(igrd_meteo,1)%n2dat = ?
    END SELECT
ENDIF

!---sst
IF (iopt_temp_sbc.GT.1) THEN
   surfacegrids(igrd_sst,1)%nhtype = ?
   SELECT CASE (surfacegrids(igrd_sst,1)%nhtype)
      CASE (1)
         surfacegrids(igrd_sst,1)%n1dat = ?
         surfacegrids(igrd_sst,1)%n2dat = ?
         surfacegrids(igrd_sst,1)%delxdat = ?
         surfacegrids(igrd_sst,1)%delydat = ?
         surfacegrids(igrd_sst,1)%x0dat = ?
         surfacegrids(igrd_sst,1)%y0dat = ?
      CASE (2:3)
         surfacegrids(igrd_sst,1)%n1dat = ?
         surfacegrids(igrd_sst,1)%n2dat = ?
   END SELECT
ENDIF

!---biology
IF (iopt_biolgy.GT.0) THEN
   ifil_610: DO ifil=1,MaxGridFiles
      surfacegrids(igrd_bio,ifil)%nhtype = ?
      SELECT CASE (surfacegrids(igrd_bio,ifil)%nhtype)
         CASE (1)
            surfacegrids(igrd_bio,ifil)%n1dat = ?
            surfacegrids(igrd_bio,ifil)%n2dat = ?
            surfacegrids(igrd_bio,ifil)%delxdat = ?
            surfacegrids(igrd_bio,ifil)%delydat = ?
            surfacegrids(igrd_bio,ifil)%x0dat = ?
            surfacegrids(igrd_bio,ifil)%y0dat = ?
         CASE (2:3)
            surfacegrids(igrd_bio,ifil)%n1dat = ?
            surfacegrids(igrd_bio,ifil)%n2dat = ?
      END SELECT
   ENDDO ifil_610
ENDIF

!
!7. Forcing input file properties
!--------------------------------
!
!7.1 Status attribute
!-------------------
!
!---domain decomposition ('N','R')
IF (parallel_set.AND.iopt_MPI_partit.EQ.2) THEN
   modfiles(io_mppmod,1,1)%status = ?
ENDIF

!---model grid ('N','R')
modfiles(io_modgrd,1,1)%status = '0'

!---initial conditions ('0,'N','R')
modfiles(io_inicon,ics_phys,1)%status = '0'

!---external surface grids ('N','R')
IF (surfacegrids(igrd_meteo,1)%nhtype.EQ.2.OR.&
  & surfacegrids(igrd_meteo,1)%nhtype.EQ.3) THEN
   modfiles(io_metgrd,1,1)%status = ?
ENDIF
IF (surfacegrids(igrd_sst,1)%nhtype.EQ.2.OR.&
  & surfacegrids(igrd_sst,1)%nhtype.EQ.3) THEN
   modfiles(io_sstgrd,1,1)%status = ?
ENDIF
IF (surfacegrids(igrd_bio,1)%nhtype.EQ.2.OR.&
  & surfacegrids(igrd_bio,1)%nhtype.EQ.3) THEN
   modfiles(io_biogrd,1,1)%status = ?
ENDIF

!---nested sub-grids ('N','R')
IF (iopt_nests.EQ.1) THEN
   iset_711: DO iset=1,nonestsets
      modfiles(io_nstgrd,iset,1)%status = ?
   ENDDO iset_711
ENDIF
!---specifiers for open boundary conditions ('N','R')
IF (iopt_sur_1D.EQ.1) modfiles(io_1uvsur,1,1)%status = ?
IF (iopt_obc_2D.EQ.1) modfiles(io_2uvobc,1,1)%status = ?
IF (iopt_obc_3D.EQ.1) modfiles(io_3uvobc,1,1)%status = ?
IF (iopt_obc_sal.EQ.1) modfiles(io_salobc,1,1)%status = ?
IF (iopt_obc_temp.EQ.1) modfiles(io_tmpobc,1,1)%status = ?
IF (iopt_obc_sed.EQ.1) modfiles(io_sedobc,1,1)%status = ?
IF (iopt_obc_bio.EQ.1) modfiles(io_bioobc,1,1)%status = ?

!---open boundary data ('0','N','R')
IF (iopt_sur_1D.EQ.1) modfiles(io_1uvsur,2,1)%status = ?
ifil_712: DO ifil=2,MaxIOFiles
   IF (iopt_obc_2D.EQ.1) modfiles(io_2uvobc,ifil,1)%status = ? 
   IF (iopt_obc_3D.EQ.1) modfiles(io_3uvobc,ifil,1)%status = ?
   IF (iopt_obc_sal.EQ.1) modfiles(io_salobc,ifil,1)%status = ? 
   IF (iopt_obc_temp.EQ.1) modfiles(io_tmpobc,ifil,1)%status = ?  
   IF (iopt_obc_sed.EQ.1) modfiles(io_sedobc,ifil,1)%status = ?  
   IF (iopt_obc_bio.EQ.1) modfiles(io_bioobc,ifil,1)%status = ?  
ENDDO ifil_712

!---structures
IF (iopt_drycell.EQ.1) modfiles(io_drycel,1,1)%status = ?
IF (iopt_weirbar.EQ.1) modfiles(io_weirbar,1,1)%status = ?
IF (iopt_dischr.EQ.1) modfiles(io_dischr,1,1)%status = ?


!---specifiers for nesting ('N','R')
IF (iopt_nests.EQ.1) modfiles(io_nstspc,1,1)%status = ?

!---surface data ('N','R')
IF (iopt_meteo.EQ.1) modfiles(io_metsur,1,1)%status = ?
IF (iopt_temp_sbc.GT.1) modfiles(io_sstsur,1,1)%status = ?
IF (iopt_waves.GT.0) modfiles(io_wavsur,1,1)%status = ?
IF (iopt_biolgy.GT.0) modfiles(io_biosur,1,1)%status = ?


!
!7.2 Other atrributes
!--------------------
!

iddesc_720: DO iddesc=1,MaxIOTypes
ifil_720: DO ifil=1,MaxIOFiles
   IF (modfiles(iddesc,ifil,1)%status.NE.'0') THEN
!     ---file format ('A','U','N')
      modfiles(iddesc,ifil,1)%form = ?
!     ---(non-default) filename (string of length 'leniofile')
      modfiles(iddesc,ifil,1)%filename = ''
!     ---update times (time series files only)
      modfiles(iddesc,ifil,1)%tlims(1:3) = ?
!     ---info file
      modfiles(iddesc,ifil,1)%info = .FALSE.
!     ---procedure in case of end of file condition and time series (0/1/2)
      modfiles(iddesc,ifil,1)%endfile = 0
   ENDIF
ENDDO ifil_720
ENDDO iddesc_720

!
!8. Forcing output file properties
!--------------------------------
!

iddesc_810: DO iddesc=1,MaxIOTypes
   SELECT CASE (iddesc)
      CASE (io_2uvnst,io_3uvnst,io_salnst,io_tmpnst,io_sednst,io_bionst)
         iset_811: DO iset=1,nonestsets
!           ---file format ('A','U','N')
            modfiles(iddesc,ifil,2)%form = ?
!           ---file name
            modfiles(iddesc,ifil,2)%filename = ?
!           ---info file
            modfiles(iddesc,ifil,2)%info = .FALSE.
         ENDDO iset_811
      CASE DEFAULT
         ifil_812: DO ifil=1,MaxIOFiles
            IF (modfiles(iddesc,ifil,1)%status.NE.'0') THEN
!              ---file format ('A','U','N')
               modfiles(iddesc,ifil,2)%form = ?
!              ---file name
               modfiles(iddesc,ifil,2)%filename = ''
!              ---info file
               modfiles(iddesc,ifil,2)%info = .FALSE.
            ENDIF
         ENDDO ifil_812
   END SELECT
ENDDO iddesc_810

!
!9. CIF files
!------------
!

ciffiles(icif_sed)%status = ?
ciffiles(icif_bio)%status = ?

!
!10. User output
!---------------
!
!---time series output
IF (iopt_out_tsers.EQ.1) THEN
   nosetstsr = ?
   nostatstsr = ?
   novarstsr = ?
ENDIF

!---time-averaged output
IF (iopt_out_avrgd.EQ.1)
   nosetsavr = ?
   nostatsavr = ?
   novarsavr = ?
ENDIF

!---harmonic analysis
IF (iopt_out_anal.EQ.1) THEN
   nosetsanal = ?
   nofreqsanal = ?
   nostatsanal = ?
   novarsanal = ?
ENDIF

!---titles
intitle = runtitle
outtitle = runtitle

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_mod_params

!========================================================================

SUBROUTINE usrdef_grid
!************************************************************************
!
! *usrdef_grid* Define model grid arrays (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.5.1
!
! Description - called if modfiles(io_modgrd,1,1)%status = 'N'
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - 
!
! Module calls -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!1. Coordinate arrays
!--------------------
!
!---horizontal [m or degrees]
IF (iopt_grid_htype.EQ.2) THEN
   gdelxglb(1:nc) = ?
   gdelyglb(1:nr) = ?
ELSEIF (iopt_grid_htype.EQ.3) THEN
   gxcoordglb(1:nc,1:nr) = ?
   gycoordglb(1:nc,1:nr) = ?
ENDIF

!---vertical coordinates [-]
IF (iopt_grid_vtype.EQ.2) THEN
   gsigcoord = ?
ELSEIF (iopt_grid_vtype.EQ.3) THEN
   gscoordglb(1:nc-1,1:nr-1,:) = ?
ENDIF

!
!2. Water depths [m]
!-------------------
!

depmeanglb(1:nc-1,1:nr-1) = ?

!
!3. Open boundary locations
!--------------------------
!
!---U-nodes
IF (nobu.GT.0) THEN
   iobu(1:nobu) = ?
   jobu(1:nobu) = ?
ENDIF

!---V-nodes
IF (nobv.GT.0) THEN
   iobv(1:nobv) = ?
   jobv(1:nobv) = ?
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_grid

!========================================================================

SUBROUTINE usrdef_partition
!************************************************************************
!
! *usrdef_partition* Define domain decomposition (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - if parallel mode is on, iopt_MPI_partit=2 and
!               modfiles(io_mppmod,1,1)%status = 'N'
!
! Reference -
!
! Calling program - domain_decomposition
!
! External calls - 
!
! Module calls -
!
!************************************************************************
!
USE iopars
USE paralpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iproc


IF (parallel_set) THEN

   procname(pglev+1) = 'usrdef_partition'
   CALL log_timer_in()

   iproc_110: DO iproc=1,nprocs
      nc1procs(iproc) = ?
      nc2procs(iproc) = ?
      nr1procs(iproc) = ?
      nr2procs(iproc) = ?
   ENDDO iproc_110

   CALL log_timer_out()

ENDIF


RETURN

END SUBROUTINE usrdef_partition

!========================================================================

SUBROUTINE usrdef_phsics
!************************************************************************
!
! *usrdef_phsics* Define initial conditions for physics (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - initial arrays are defined on local grid including halos
!             - if not defined, each initial array is set to zero by default
!             - the example code here considers two methods
!      1: arrays are defined locally, e.g. by reading from an external file
!         (always used in serial mode)
!      2: arrays are first defined globally and distributed
!         to the local domains afterwards. The global data in the example are
!         obtained from an external data file 
! Reference -
!
! Calling program - define_phsics
!
! External calls - exchange_physics
!
! Module calls - close_filepars, copy_vars, distribute_mod, error_alloc,
!                open_filepars
!
!************************************************************************
!
USE currents
USE density
USE depths
USE fluxes
USE grid
USE gridpars
USE iopars
USE modids
USE obconds
USE physpars
USE structures
USE switches
USE syspars
USE tide
USE turbulence
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: close_filepars, open_filepars
USE paral_comms, ONLY: copy_vars, distribute_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: method = ?
INTEGER :: iunit, k
INTEGER, DIMENSION(3) :: lbounds 
INTEGER, DIMENSION(4) :: nhdist = 0
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: realglb2d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realglb3d


procname(pglev+1) = 'usrdef_phsics'
CALL log_timer_in()


!
!1. "Local" method
!-----------------
!

IF (method.EQ.1) THEN

!
!1.1 2-D mode
!------------
!
!  ---depth-integrated currents
   udvel(1:ncloc,1:nrloc) = ?
   vdvel(1:ncloc,1:nrloc) = ?

!  ---surface elevations
   zeta(1:ncloc,1:nrloc) = ?

!
!1.2 3-D currents
!----------------
!

   IF (iopt_grid_nodim.NE.2) THEN
      uvel(1:ncloc,1:nrloc,1:nz) = ?
      vvel(1:ncloc,1:nrloc,1:nz) = ?
   ENDIF
   IF (iopt_grid_nodim.EQ.3) THEN
      wvel(1:ncloc,1:nrloc,1:nz+1) = ?
   ENDIF

!
!1.3 Density arrays
!------------------
!

   IF (iopt_temp.GT.0) temp(1:ncloc,1:nrloc,1:nz) = temp_ref
   IF (iopt_sal.GT.0) sal(1:ncloc,1:nrloc,1:nz) = sal_ref

!
!1.4 Turbulence arrays
!---------------------
!

   IF (iopt_vdif_coef.EQ.3) THEN
      tke(1:ncloc,1:nrloc,1:nz+1) = ?
      IF (iopt_turb_ntrans.EQ.2.AND.iopt_turb_param.EQ.1) THEN
         zlmix(1:ncloc,1:nrloc,1:nz+1) = ?
      ELSEIF (iopt_turb_ntrans.EQ.2.AND.iopt_turb_param.EQ.2) THEN
         dissip(1:ncloc,1:nrloc,1:nz+1) = ?
      ENDIF
   ENDIF

!
!1.5 Bottom stress arrays
!------------------------
!

   IF (iopt_bstres_form.EQ.2) THEN
      IF (iopt_bstres_drag.EQ.2) THEN
         bdragcoefatc(1:ncloc,1:nrloc) = ?
      ELSEIF (iopt_bstres_drag.EQ.4) THEN
         zroughatc(1:ncloc,1:nrloc) = ?
      ENDIF
   ENDIF
   
!
!1.6 Tidal arrays
!----------------
!

   IF (nconobc.GT.0) THEN
      IF (iopt_astro_pars.EQ.0) phase_obc(1:nconobc) = ?
      IF (iopt_astro_pars.GT.1) fnode_obc(1:nconobc) = ?

   ENDIF
   IF (iopt_astro_tide.EQ.1.AND.nconastro.GT.0) THEN
      IF (iopt_astro_pars.EQ.0) phase_astro(1:nconastro) = ?
      IF (iopt_astro_pars.GT.1) fnode_astro(1:nconastro) = ?
   ENDIF

!
!1.7 Energy loss at weirs
!------------------------
!

   IF (iopt_weirbar.EQ.1) THEN
      IF (numwbaru.GT.0) wbarelossu(1:numwbaru) = ?
      IF (numwbarv.GT.0) wbarelossv(1:numwbarv) = ?
   ENDIF

!
!1.8 Open boundary arrays
!------------------------
!
!  ---2-D mode
   IF (iopt_obc_2D.EQ.1) THEN
      obc2uvatu(1:nobu,1:2) = ?
      obc2uvatv(1:nobv,1:2) = ?
   ENDIF

!  ---3-D mode
   IF (iopt_obc_3D.EQ.1) THEN
      obc3uvatu(1:nobu,1:nz,1:2) = ?
      obc3uvatv(1:nobv,1:nz,1:2) = ?
   ENDIF

!  ---salinity
   IF (iopt_obc_sal.EQ.1) THEN
      obcsalatu(1:nobu,1:nz,0:2) = ?
      obcsalatv(1:nobv,1:nz,0:2) = ?
   ENDIF

!  ---temperature
   IF (iopt_obc_temp.EQ.1) THEN
      obctmpatu(1:nobu,1:nz,0:2) = ?
      obctmpatv(1:nobv,1:nz,0:2) = ?
   ENDIF

!
!2. "Global" method
!------------------
!

ELSEIF (method.EQ.2) THEN

!
!2.1 Allocate work space arrays
!------------------------------
!

   ALLOCATE (realglb2d(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo),STAT=errstat)
   CALL error_alloc('realglb2d',2,(/nc+2*nhalo,nr+2*nhalo/),real_type)
   realglb2d = 0.0
   ALLOCATE (realglb3d(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo,nz+1),STAT=errstat)
   CALL error_alloc('realglb3d',3,(/nc+2*nhalo,nr+2*nhalo,nz+1/),real_type)
   realglb3d = 0.0

!
!2.2 Open data file
!------------------
!
!  ---it is assumed that the name of the data file is defined in 
!     usrdef_mod_params
   CALL open_filepars(modfiles(io_inicon,ics_phys,1))
   iunit = modfiles(io_inicon,ics_phys,1)%iunit

!
!2.3 2-D mode
!------------
!
!  ---depth-integrated currents
   lbounds(1:2) = 1-nhalo
   realglb2d = 0.0
   READ (iunit,?) realglb2d(1:nc,1:nr)
   CALL distribute_mod(realglb2d,udvel,lbounds(1:2),nhdist,iarr_udvel,0.0)
   realglb2d = 0.0
   READ (iunit,?) realglb2d(1:nc,1:nr)
   CALL distribute_mod(realglb2d,vdvel,lbounds(1:2),nhdist,iarr_vdvel,0.0)

!  ---surface elevations
   lbounds(1:2) = 0
   READ (iunit,?) realglb2d(1:nc-1,1:nr-1)
   CALL distribute_mod(realglb2d,zeta,lbounds(1:2),nhdist,iarr_zeta,0.0)

!
!2.4 3-D currents
!----------------
!

   IF (iopt_grid_nodim.EQ.1.OR.iopt_grid_nodim.EQ.3) THEN
!     ---X-current
      lbounds = (/1-nhalo,1-nhalo,1/)
      realglb3d = 0.0
      READ (iunit,?) realglb3d(1:nc,1:nr,1:nz)
      CALL distribute_mod(realglb3d(:,:,1:nz),uvel,lbounds,nhdist,iarr_uvel,0.0)
!     ---Y-current
      realglb3d = 0.0
      READ (iunit) realglb3d(1:nc,1:nr,1:nz)
      CALL distribute_mod(realglb3d(:,:,1:nz),vvel,lbounds,nhdist,iarr_vvel,0.0)
   ENDIF

   IF (iopt_grid_nodim.EQ.3) THEN
!     ---vertical current
      realglb3d = 0.0
      lbounds  = (/0,0,1/)
      READ (iunit,?) realglb3d(1:nc-1,1:nr-1,:)
      CALL distribute_mod(realglb3d(0:nc,0:nr,:),wvel,lbounds,nhdist,iarr_wvel,&
                        & 0.0)
   ENDIF

!
!2.5 Density arrays
!------------------
!
!  ---temperature
   IF (iopt_temp.GT.0) THEN
      lbounds = (/1-nhalo,1-nhalo,1/)
      realglb3d = 0.0
      READ (iunit,?) realglb3d(1:nc-1,1:nr-1,1:nz)
      CALL distribute_mod(realglb3d(:,:,1:nz),temp,lbounds,nhdist,iarr_temp,0.0)
   ENDIF

!  ---salinity
   IF (iopt_sal.GT.0) THEN
      lbounds = (/1-nhalo,1-nhalo,1/)
      realglb3d = 0.0
      READ (iunit,?) realglb3d(1:nc-1,1:nr-1,1:nz)
      CALL distribute_mod(realglb3d(:,:,1:nz),sal,lbounds,nhdist,iarr_sal,0.0)
   ENDIF

!
!2.6 Turbulence arrays
!---------------------
!

   IF (iopt_vdif_coef.EQ.3) THEN

!     ---turbulence energy
      lbounds = (/1-nhalo,1-nhalo,1/)
      realglb3d = 0.0
      READ (iunit,?) realglb3d(1:nc-1,1:nr-1,2:nz)
      CALL distribute_mod(realglb3d,tke,lbounds,nhdist,iarr_tke,0.0)

!     ---mixing length
      IF (iopt_turb_ntrans.EQ.2.AND.iopt_turb_param.EQ.1) THEN
         lbounds = (/1-nhalo,1-nhalo,1/)
         realglb3d = 0.0
         READ (iunit,?) realglb3d(1:nc-1,1:nr-1,2:nz)
         CALL distribute_mod(realglb3d,zlmix,lbounds,nhdist,iarr_zlmix,0.0)

!     ---dissipation rate
      ELSEIF (iopt_turb_ntrans.EQ.2.AND.iopt_turb_param.EQ.2) THEN
         lbounds = (/1-nhalo,1-nhalo,1/)
         realglb3d = 0.0
         READ (iunit,?) realglb3d(1:nc-1,1:nr-1,2:nz)
         CALL distribute_mod(realglb3d,dissip,lbounds,nhdist,iarr_dissip,0.0)
      ENDIF

   ENDIF

!
!2.7 Bottom stress arrays
!------------------------
!

   IF (iopt_bstres_form.EQ.2) THEN
      lbounds(1:2) = 1
      realglb2d = 0.0
      READ (iunit,?) realglb2d(1:nc-1,1:nr-1)
      IF (iopt_bstres_drag.EQ.2) THEN
         CALL distribute_mod(realglb2d,bdragcoefatc,lbounds(1:2),nhdist,&
                           & iarr_bdragcoefatc,0.0)
      ELSEIF (iopt_bstres_drag.EQ.4) THEN
         CALL distribute_mod(realglb2d,zroughatc,lbounds(1:2),nhdist,&
                           & iarr_zroughatc,0.0)
      ENDIF
   ENDIF
   
!
!2.8 Tidal arrays
!----------------
!

   IF (nconobc.GT.0) THEN
      READ (iunit,?) phase_obc
      READ (iunit,?) fnode_obc
   ENDIF
   IF (iopt_astro_tide.EQ.1) THEN
      READ (iunit,?) phase_astro
      READ (iunit,?) fnode_astro
   ENDIF

!
!2.9 Energy losses terms at weirs
!--------------------------------
!
  IF (iopt_weibar.EQ.1) THEN
      IF (numwbaru.GT.0) READ (iunit,?) wbarelossu
      IF (numwbarv.GT.0) READ (iunit,?) wbarelossv
   ENDIF

!
!2.10 Open boundary arrays
!-------------------------
!
!  ---2-D mode
   IF (iopt_obc_2D.EQ.1) THEN
      IF (nobu.GT.0) THEN
         READ (iunit,?) obc2uvatu(1:nobu,1:2)
      ENDIF
      IF (nobv.GT.0) THEN
         READ (iunit,?) obc2uvatv(1:nobv,1:2)
      ENDIF
   ENDIF

!  ---3-D mode
   IF (iopt_obc_3D.EQ.1) THEN
      IF (nobu.GT.0) THEN
         READ (iunit,?) obc3uvatu(1:nobu,1:nz,1:2)
      ENDIF
      IF (nobv.GT.0) THEN
         READ (iunit,?) obc3uvatv(1:nobv,1:nz,1:2)
      ENDIF
   ENDIF

!  ---salinity
   IF (iopt_obc_sal.EQ.1) THEN
      IF (nobu.GT.0) THEN
         READ (iunit,?) obcsalatu(1:nobu,1:nz,1:2)
      ENDIF
      IF (nobv.GT.0) THEN
         READ (iunit,?) obcsalatv(1:nobv,1:nz,1:2)
      ENDIF
   ENDIF

!  ---temperature
   IF (iopt_obc_temp.EQ.1) THEN
      IF (nobu.GT.0) THEN
         READ (iunit,?) obctmpatu(1:nobu,1:nz,1:2)
      ENDIF
      IF (nobv.GT.0) THEN
         READ (iunit,?) obctmpatv(1:nobv,1:nz,1:2)
      ENDIF
   ENDIF

!
!2.11 Close data file
!--------------------
!

   CALL close_filepars(modfiles(io_inicon,ics_phys,1))

!
!2.12 Deallocate
!---------------
!

   DEALLOCATE (realglb2d,realgbl3d)

!
!2.13 Flag land areas (if needed)
!--------------------------------
!
!  ---2-D mode
   WHERE (nodeatu.EQ.0) udvel = 0.0
   WHERE (nodeatv.EQ.0) vdvel = 0.0
   WHERE (nodeatc(0:ncloc+1,0:nrloc+1).EQ.0) zeta = 0.0

!  ---3-D currents
   IF (iopt_grid_nodim.NE.2) THEN
      k_2131: DO k=1,nz
         WHERE (nodeatu.EQ.0) uvel(:,:,k) = 0.0
         WHERE (nodeatv.EQ.0) vvel(:,:,k) = 0.0
      ENDDO k_2131
   ENDIF
   IF (iopt_grid_nodim.EQ.3) THEN
      k_2312: DO k=1,nz+1
         WHERE (nodeatc(0:ncloc,0:nrloc.EQ.0).EQ.0) wvel(:,:,k) = 0.0
      ENDDO k_2312
   ENDIF

!  ---density arrays
   k_2133: DO k=1,nz
      IF (iopt_temp.GT.0) THEN
         WHERE (nodeatc.EQ.0) temp(:,:,k) = 0.0
      ENDIF
      IF (iopt_sal.GT.0) THEN
         WHERE (nodeatc.EQ.0) sal(:,:,k) = 0.0
      ENDIF
   ENDDO k_2133

!  ---turbulence arrays
   k_2134: DO k=1,nz+1
      WHERE (nodeatc.EQ.0) tke(:,:,k) = 0.0
      IF (iopt_turb_ntrans.EQ.2.AND.iopt_turb_param.EQ.1) THEN
         WHERE (nodeatc.EQ.0) zlmix(:,:,k) = 0.0
      ELSEIF (iopt_turb_ntrans.EQ.2.AND.iopt_turb_param.EQ.2) THEN
         WHERE (nodeatc.EQ.0) dissip(:,:,k) = 0.0
      ENDIF
   ENDDO k_2134

ENDIF
   
CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_phsics

!========================================================================

SUBROUTINE usrdef_1dsur_spec
!************************************************************************
!
! *usrdef_1dsur_spec* Define specifier arrays for 1-D forcing (example routine)
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
! External calls - 
!
! Module calls -
!
!************************************************************************
!
USE iopars
USE obconds
USE tide
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_1dsur_spec'
CALL log_timer_in()

!---amplitudes and phases of X-component of pressure gradient
gxslope_amp(1:nconobc) = ?
gyslope_pha(1:nconobc) = ?

!---amplitudes and phases of Y-component of pressure gradient 
gyslope_amp(1:nconobc) = ?
gyslope_pha(1:nconobc) = ?
     
!---surface elevation
zeta_amp(1:nconobc) = ?
zeta_pha(1:nconobc) = ?

!---type of data in (eventual) data file (1/2/3)
isur1dtype = ? 

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_1dsur_spec

!========================================================================

SUBROUTINE usrdef_2dobc_spec(iddesc)
!************************************************************************
!
! *usrdef_2dobc_spec* Define specifier arrays for 2-D open boundary conditions
!                     (example routine)
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
! External calls - 
!
! Module calls -
!
!************************************************************************
!
USE gridpars
USE iopars
USE obconds
USE switches
USE tide
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id (io_2uvobc or io_2xyobc)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ifil, l


procname(pglev+1) = 'usrdef_2dobc_spec'
CALL log_timer_in()

SELECT CASE (iddesc)

!
!1. Normal conditions
!--------------------
!

   CASE (io_2uvobc)

!
!1.1 Type of conditions
!----------------------
!
!     ---type of conditions (0-13)
      ityp2dobu(1:nobu) = ?
      ityp2dobv(1:nobv) = ?

!     ---location of elevation point (0/1/2)
      iloczobu(1:nobu) = ?
      iloczobv(1:nobv) = ?

!
!1.2 Specifiers for data files
!-----------------------------
!

      ifil_120: DO ifil=2,maxdatafiles(io_2uvobc,1)
!        ---number of data per file
         no2dobuv(ifil) = ?
!        ---type of data (1/2/3)
         iobc2dtype(ifil) = ?
!        ---map between data and open boundary index
         l_121: DO l=1,no2dobuv(ifil)
            index2dobuv(l,ifil) = ?
         ENDDO l_121
      ENDDO ifil_120

!
!1.3 Harmonic data
!------------------
!
!     ---amplitudes
      udatobu_amp(1:nobu,1:nconobc) = ?
      vdatobv_amp(1:nobv,1:nconobc) = ?
      zdatobu_amp(1:nobu,1:nconobc) = ?
      zdatobv_amp(1:nobv,1:nconobc) = ?
      
!     ---phases
      udatobu_pha(1:nobu,1:nconobc) = ?
      vdatobv_pha(1:nobv,1:nconobc) = ?
      zdatobu_pha(1:nobu,1:nconobc) = ?
      zdatobv_pha(1:nobv,1:nconobc) = ?

!
!1.4 Discharge along open boundary sections
!------------------------------------------
!

      IF (nqsecobu.GT.0) iqsecobu(1:noqsecobu,1:2) = ?
      IF (nqsecobv.GT.0) jqsecobv(1:noqsecobv,1:2) = ?

!
!2. Tangential conditions
!------------------------
!

   CASE (io_2xyobc)

!
!2.1 Type of conditions
!----------------------
!
!     ---type of conditions (0-13)
      ityp2dobx(1:nobx) = ?
      ityp2doby(1:noby) = ?

!
!2.2 Specifiers for data files
!-----------------------------
!

      ifil_220: DO ifil=2,maxdatafiles(io_2xyobc,1)
!        ---number of data per file
         no2dobxy(ifil) = ?
!        ---map between data and open boundary index
         l_221: DO l=1,no2dobxy(ifil)
            index2dobxy(l,ifil) = ?
         ENDDO l_121
      ENDDO ifil_120

!
!2.3 Harmonic data
!------------------
!
!     ---amplitudes
      vdatobx_amp(1:nobx,1:nconobc) = ?
      udatoby_amp(1:noby,1:nconobc) = ?
      
!     ---phases
      vdatobx_pha(1:nobx,1:nconobc) = ?
      udatoby_pha(1:noby,1:nconobc) = ?

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
!                       (example routine)
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
!*Local variables
!
INTEGER :: ifil, l


procname(pglev+1) = 'usrdef_profobc_spec'
CALL log_timer_in()

SELECT CASE (iddesc)

!
!1. 3-D normal baroclinic current
!--------------------------------
!

   CASE (io_3uvobc)

!     ---type of condition (0/1/2/3/4)
      itypobux(1:nobu) = ?
      itypobvy(1:nobv) = ?

!     ---profile number
      iprofobux(1:nobu,1) = ?
      iprofobvy(1:nobv,1) = ?

!     ---specifiers for data files
      ifil_110: DO ifil=2,nofiles
!        --number of profiles per file
         noprofsd(ifil) = ?
!        --maps between data index and global profile/index number
         l_111: DO l=1,noprofsd(ifil)
            indexprof(l,ifil) = ?
            indexvar(l,ifil) = ?
         ENDDO l_111
      ENDDO ifil_110
      iprofrlx(1:norlxzones) = ?

!
!2. 3-D tangential baroclinic current
!------------------------------------
!


!     ---type of condition (0/1/2)
      itypobux(1:nobu) = ?
      itypobvy(1:nobv) = ?

!     ---profile number
      ...


!      
!3. Salinity
!-----------
!

   CASE (io_salobc)

!     ---type of condition (0/1/2)
      itypobux(1:nobu) = ?
      itypobvy(1:nobv) = ?

!     ---profile number
      ...

!
!4. Temperature
!--------------
!

   CASE (io_tmpobc)
      ...

!
!5. Sediments
!------------
!

   CASE (io_sedobc)
      ...

!
!6. Biology
!----------
!

   CASE (io_bioobc)
      ...

END SELECT

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_profobc_spec

!========================================================================

SUBROUTINE usrdef_1dsur_data(ciodatetime,data1d,novars)
!************************************************************************
!
! *usrdef_1dsur_data* Define data for 1-D forcing (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - the example code considers two methods
!      1: data are obtained from an external input file
!      2: data are obtained from some arbitrary source
!
! Reference -
!
! Calling program - define_1dsur_data
!
! External calls - 
!
! Module calls - open_filepars
!
!************************************************************************
!
USE iopars
USE syspars
USE inout_routines, ONLY: open_filepars
USE time_routines, ONLY: log_timer_in, log_timer_out

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
!*Local variables
!
INTEGER :: method = ?
INTEGER :: iunit


procname(pglev+1) = 'usrdef_1dsur_data'
CALL log_timer_in()

!
!1. Data from external data file
!-------------------------------
!

IF (method.EQ.1) THEN

!  ---open input file
   IF (modfiles(io_1uvsur,2,1)%status.EQ.0) THEN
      CALL open_filepars(modfiles(io_1uvsur,2,1))
      GOTO 1000
   ENDIF

!  ---read data
   iunit = modfiles(io_1uvsur,2,1)%iunit
   READ (iunit,?,END=99) ciodatetime
   READ (iunit,?) data1d

!
!2. Data from any source
!-----------------------
!

ELSEIF (method.EQ.2) THEN

!  ---set open status on first call
   IF (modfiles(io_1uvsur,2,1)%status.EQ.0) THEN
      modfiles(io_1uvsur,2,1)%iostat = 1
      GOTO 1000
   ENDIF

!  ---obtain data
   ciodatetime = ?
   data1d = ?

!  ---reset status to end of file condition at last call
   modfiles(io_1uvsur,2,1)%iostat = 2

ENDIF

1000 CALL log_timer_out()


RETURN

99 modfiles(io_1uvsur,2,1)%iostat = 2

END SUBROUTINE usrdef_1dsur_data

!========================================================================

SUBROUTINE usrdef_2dobc_data(iddesc,ifil,ciodatetime,data2d,nodat,novars)
!************************************************************************
!
! *usrdef_2dobc_data* Define open boundary data for 2-D mode (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - the example code considers two methods
!      1: data are obtained from an external input file
!      2: data are obtained from some arbitrary source
!
! Reference -
!
! Calling program - define_2dobc_data
!
! External calls -
!
! Module calls - open_filepars
!
!************************************************************************
!
USE iopars
USE syspars
USE inout_routines, ONLY: open_filepars
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
INTEGER :: method = ?
INTEGER :: iunit, ivar


procname(pglev+1) = 'usrdef_2dobc_data'
CALL log_timer_in()

!
!1. Data from external data file
!-------------------------------
!

IF (method.EQ.1) THEN

!  ---open input file
   IF (modfiles(iddesc,ifil,1)%status.EQ.0) THEN
      CALL open_filepars(modfiles(iddesc,ifil,1))
      GOTO 1000
   ENDIF

!  ---read data
   iunit = modfiles(iddesc,ifil,1)%iunit
   READ (iunit,?,END=99) ciodatetime
   ivar_110: DO ivar=1,novars
      READ (iunit,?) data2d(:,ivar)
   ENDDO ivar_110

!
!2. Data from any source
!-----------------------
!

ELSEIF (method.EQ.2) THEN

!  ---set open status on first call
   IF (modfiles(iddesc,ifil,1)%status.EQ.0) THEN
      modfiles(iddesc,ifil,1)%iostat = 1
      GOTO 1000
   ENDIF

!  ---obtain data
   ciodatetime = ?
   ivar_210: DO ivar=1,novars
      data2d(:,ivar) = ?
   ENDDO ivar_210

!  ---reset status to end of file condition at last call
   modfiles(iddesc,ifil,1)%iostat = 2

ENDIF

1000 CALL log_timer_out()


RETURN

99 modfiles(iddesc,ifil,1)%iostat = 2

END SUBROUTINE usrdef_2dobc_data

!========================================================================

SUBROUTINE usrdef_profobc_data(iddesc,ifil,ciodatetime,psiprofdat,numprofs)
!************************************************************************
!
! *usrdef_profobc_data* Define physical open boundary profiles (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.2
!
! Description - called either with iddesc = io_3uvobc, io_salobc, io_tmpobc or
!               io_bioobc
!             - the example code considers two methods
!      1: data are obtained from an external input file
!      2: data are obtained from some arbitrary source
!
! Reference -
!
! Calling program - define_profobc_data
!
! Module calls - open_filepars
!
!************************************************************************
!
USE gridpars
USE iopars
USE syspars
USE inout_routines, ONLY: open_filepars
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
!* Local variables
!
INTEGER :: method = ?
INTEGER :: iunit


procname(pglev+1) = 'usrdef_profobc_data'
CALL log_timer_in()

!
!1. Data from external data file
!-------------------------------
!

IF (method.EQ.1) THEN

!  ---open input file
   IF (modfiles(iddesc,ifil,1)%status.EQ.0) THEN
      CALL open_filepars(modfiles(io_iddesc,ifil,1))
      GOTO 1000
   ENDIF

!  ---read data
   iunit = modfiles(iddesc,ifil,1)%iunit
   READ (iunit,?,END=99) ciodatetime
   READ (iunit,?) psiprofdat

!
!2. Data from any source
!-----------------------
!

ELSEIF (method.EQ.2) THEN

!  ---set open status on first call
   IF (modfiles(iddesc,ifil,1)%status.EQ.0) THEN
      modfiles(iddesc,ifil,1)%iostat = 1
      GOTO 1000
   ENDIF

!  ---obtain data
   ciodatetime = ?
   psiprofdat = ?

!  ---reset status to end of file condition at last call
   modfiles(iddesc,ifil,1)%iostat = 2

ENDIF

1000 CALL log_timer_out()


RETURN

99 modfiles(iddesc,ifil,1)%iostat = 2

END SUBROUTINE usrdef_profobc_data

!========================================================================

SUBROUTINE usrdef_rlxobc_spec
!************************************************************************
!
! *usrdef_rlxobc_spec* Define specifier arrays for relaxation (example routine)
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
! External calls - 
!
! Module calls -
!
!************************************************************************
!
USE gridpars
USE iopars
USE relaxation
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_rlxobc_spec'
CALL log_timer_in()

inodesrlx(1:2) = ?
idirrlx(1:norlxzones) = ?
ityprlx(1:norlxzones) = ?
iposrlx(1:norlxzones) = ?
jposrlx(1:norlxzones) = ?
ncrlx(1:norlxzones) = ?
nrrlx(1:norlxzones) = ?

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_rlxobc_spec
