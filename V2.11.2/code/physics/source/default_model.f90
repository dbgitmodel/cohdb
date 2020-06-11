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

MODULE default_model
!************************************************************************
!
! *default_model* Default model settings
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)default_model.f90  V2.11.2
!
! $Date: 2018-09-28 13:01:37 +0200 (Fri, 28 Sep 2018) $
!
! $Revision: 1188 $
!
! Description - 
!
! Reference -
!
! Generic routines - default_out_files
!
! Routines - default_ellvars, default_init_params, default_mod_params,
!            default_out_gpars
!
!************************************************************************
!
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

INTERFACE default_out_files
   MODULE PROCEDURE default_out_files_1d, default_out_files_2d
END INTERFACE


CONTAINS

!========================================================================

SUBROUTINE default_ellvars(ellvars)
!************************************************************************
!
! *default_ellvars* Default attributes for elliptic variables 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)default_model.f90  V2.0
!
! Description -
!
! Module calls - inquire_var
!
!************************************************************************
!
USE datatypes
USE modids
USE modvars_routines, ONLY: inquire_var

!
!*Arguments
!
TYPE (VariableAtts), INTENT(OUT), DIMENSION(14) :: ellvars

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ellvars* DERIVED Attributes of elliptic variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar, varid
INTEGER, DIMENSION(14) :: iellids


procname(pglev+1) = 'default_ellvars'
CALL log_timer_in()

iellids(1:7) = (/iarr_ellmaj2d,iarr_ellmin2d,iarr_ellip2d,iarr_ellinc2d,&
               & iarr_ellpha2d,iarr_ellcc2d,iarr_ellac2d/) 
iellids(8:14)= (/iarr_ellmaj3d,iarr_ellmin3d,iarr_ellip3d,iarr_ellinc3d,&
               & iarr_ellpha3d,iarr_ellcc3d,iarr_ellac3d/)

ivar_100: DO ivar=1,14
   varid = iellids(ivar)
   CALL inquire_var(varid,varatts=ellvars(ivar))
ENDDO ivar_100

CALL log_timer_out()
 

RETURN

END SUBROUTINE default_ellvars

!========================================================================

SUBROUTINE default_init_params
!************************************************************************
!
! *default_init_params* Default values for monitoring parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)default_model.f90  V2.11
!
! Description - parameters can be reset by the user in usrdef_init_params
!
! Reference -
!
! Module calls - error_alloc
!
!************************************************************************
!
USE paralpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc


!
!1. Cold/warm start
!------------------
!

cold_start = .FALSE.

!
!2. Log files
!------------
!
!---allocate
ALLOCATE (levprocs_ini(npworld),STAT=errstat)
CALL error_alloc('levprocs_ini',1,(/npworld/),kndint)
levprocs_ini = 0
ALLOCATE (levprocs_run(npworld),STAT=errstat)
CALL error_alloc('levprocs_run',1,(/npworld/),kndint)
levprocs_run = 0

!---default file names and resolution   
inilog_file = TRIM(runtitle)//'.inilog'
runlog_file = TRIM(runtitle)//'.runlog'

!---default parameters
exitlog = .TRUE.
runlog_count =  int_fill

!
!3. Error coding
!---------------
!
!---allocate
ALLOCATE (levprocs_err(npworld),STAT=errstat)
CALL error_alloc('levprocs_err',1,(/npworld/),kndint)
levprocs_err = 1

!---default file names and resolution   
errlog_file = TRIM(runtitle)//'.errlog'

!---default parameters
maxerrors = MaxErrMesgs

!
!4. Warnings
!-----------
!

warning = .TRUE.
warlog_file = TRIM(runtitle)//'.warlog'

!
!5. Monitoring
!-------------
!
!---implicit algorithm
monlog = 0
monlog_file = TRIM(runtitle)//'.monlog'

!---sediments
sedlog = .TRUE.
sedlog_file = TRIM(runtitle)//'.sedlog' 

!
!6. Timing
!---------
!

levtimer = 0
nopcc = 0
timing_file = TRIM(runtitle)//'.timing'
timer_format = 1
timer = .FALSE.

!
!7. Dredging and relocation log file
!-----------------------------------
!

darlog_file = TRIM(runtitle)//'.darlog'

!
!8. Parallel application
!-----------------------
!
!---numper of processes
nprocscoh = npworld; nprocspart = 0; nprocswav = 0

!---switches for model coupling
iopt_part_model = 0; iopt_waves_model = 0


RETURN

END SUBROUTINE default_init_params

!========================================================================

SUBROUTINE default_mod_params
!************************************************************************
!
! *default_mod_params* Default settings for model parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)default_model.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Module calls - filepars_init, gridpars_init
!
!************************************************************************
!
USE gridpars
USE nestgrids
USE obconds
USE paralpars
USE physpars
USE relaxation
USE structures
USE switches
USE syspars
USE tide
USE timepars
USE turbpars
USE datatypes_init, ONLY: filepars_init, gridpars_init

!
!* Local variables
!
INTEGER :: idesc, ifil, iglb, igrd, iotype


procname(pglev+1) = 'default_mod_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!
!---grid
iopt_grid_htype = 1; iopt_grid_nodim = 3; iopt_grid_sph = 0
iopt_grid_vtype = 1; iopt_grid_vtype_transf = 0

!---interpolation
iopt_arrint_depths = 1; iopt_arrint_hreg = 0; iopt_arrint_vreg = 0
iopt_arrint_3D = 0

!---hydrodynamics
iopt_curr = 2; iopt_curr_wfall = 1

!---density
iopt_dens = 0; iopt_dens_convect = 0; iopt_dens_grad = 1; iopt_sal = 0
iopt_temp = 0; iopt_temp_optic = 1; iopt_temp_sbc = 1

!---biology
iopt_biolgy = 0

!---sediments
iopt_dar = 0; iopt_morph = 0; iopt_sed = 0; iopt_tidal_accel = 0

!---particle model
iopt_part_write = 0

!---bottom stress
iopt_bstres_drag = 3; iopt_bstres_form = 2; iopt_bstres_nodim = 3
iopt_bstres_waves_bfric = 0

!---tranport
iopt_transp_full = 0

!---advection
iopt_adv_scal = 3; iopt_adv_turb = 0; iopt_adv_tvd = 1; iopt_adv_2D = 1
iopt_adv_3D = 1  

!---diffusion
iopt_hdif_coef = 0; iopt_hdif_lim = 1; iopt_hdif_scal = 0; iopt_hdif_turb = 0
iopt_hdif_2D = 0; iopt_hdif_3D = 0; iopt_kinvisc = 0; iopt_vdif_coef = 3
iopt_vdif_rot = 1

!---turbulence
iopt_turb_alg = 1; iopt_turb_dis_bbc = 2; iopt_turb_dis_sbc = 2
iopt_turb_iwlim = 0; iopt_turb_kinvisc = 0; iopt_turb_lmix = 4
iopt_turb_ntrans = 1; iopt_turb_param = 2; iopt_turb_stab_form = 3
iopt_turb_stab_lev = 1; iopt_turb_stab_mod = 4; iopt_turb_stab_tke = 2
iopt_turb_tke_bbc = 2; iopt_turb_tke_sbc = 2

!---drying/wetting
iopt_fld = 0; iopt_fld_alpha = 1

!---structures/discharges
iopt_dischr = 0; iopt_dischr_land = 1; iopt_drycel = 0; iopt_thndam = 0
iopt_weibar = 0

!---explicit/implicit integration
iopt_cor_impl = 1; iopt_hydro_impl = 0; iopt_scal_depos = 1
iopt_vadv_impl = 1; iopt_vdif_impl = 2

!---open boundary conditions
iopt_obc_advflux = 1; iopt_obc_advrlx = 0; iopt_obc_bio = 0
iopt_obc_invbar = 0; iopt_obc_relax = 0
iopt_obc_sal = 0; iopt_obc_sed = 0; iopt_obc_temp = 0; iopt_obc_th = 0
iopt_obc_2D = 0; iopt_obc_2D_tang = 0; iopt_obc_3D = 0; iopt_obc_3D_tang = 0

!---astronomical tide
iopt_astro_anal = 0; iopt_astro_pars = 2; iopt_astro_tide = 0

!---1-D applications
iopt_sur_1D = 0

!---meteo surface forcing
iopt_meteo = 0; iopt_meteo_data = 1; iopt_meteo_heat = 1; iopt_meteo_pres = 1
iopt_meteo_precip = 0; iopt_meteo_stres = 1

!---surface fluxes
iopt_sflux_pars = 1; iopt_sflux_precip = 0; iopt_sflux_qlong = 1
iopt_sflux_qshort = 1; iopt_sflux_strat = 0

!---surface waves
iopt_waves = 0; iopt_waves_couple = 0; iopt_waves_curr = 0
iopt_waves_dissip = 0; iopt_waves_extrapol = 1; iopt_waves_form = 1
iopt_waves_pres = 1

!---nesting
iopt_nests = 0

!---parallel processing (MPI)
iopt_MPI_abort = 0; iopt_MPI_comm_all = 2; iopt_MPI_comm_coll = 0
iopt_MPI_comm_exch = 2; iopt_MPI_comm_full = 0; iopt_MPI_comm_gath = 2
iopt_MPI_comm_scat = 2; iopt_MPI_partit = 1; iopt_MPI_sync = 0

!---output
iopt_out_avrgd = 0; iopt_out_anal = 0; iopt_out_tsers = 1

!---netcdf
iopt_CDF_abort = 0; iopt_CDF_fill = 0; iopt_CDF_format = 1
iopt_CDF_shared = 0; iopt_CDF_sync = 0; iopt_CDF_tlim = 1

!---multigrid
iopt_mg_cycle = 1; iopt_mg_prolong = 2; iopt_mg_smoother = 1

!---random generator
iopt_rng_seed = 1

!
!2. Model grid
!-------------
!

nc = 0; nr = 0; nz = 0; nosbu = 0; nosbv = 0; nrvbu = 0; nrvbv = 0

!
!3. Process numbers
!------------------
!

nprocsx = 0; nprocsy = 0

!
!4. Date/time parameters
!-----------------------
!

CEndDateTime = cdatetime_undef; CStartDateTime = cdatetime_undef
delt2d = 0.0; iccvt = 0; ic3d = 1; icnodal = 0
time_zone = 0.0; nt = 0; hydro_end = 0

!
!5. Physical parameters
!----------------------
!
!---number of tidal constituents
nconastro = 0; nconobc = 0

!---number of relaxation zones
norlxzones = 0

!---number of nested grids
nonestsets = 0

!---wait time
nowaitsecs = 0
maxwaitsecs = 3600

!---unit record size (in bytes)
nrecunit = 4

!---restart times
norestarts = 1
ntrestart(1) = int_fill

!---relaxation time for 2-D open boundary conditions
ntobcrlx = 0

!---general
rho_air = 1.2; specheat = 3987.5

!---reference values
atmpres_ref = 101325.0; dlat_ref = 0.0; dlon_ref = 0.0; dlon_ref_anal = 0
dlon_ref_obc = 0.0; gacc_ref = float_fill; sal_ref = 33.0; sst_ref = 12.0
temp_min = 0.0; temp_ref = 12.0

!---model grid
b_SH = 0.1; dl_BB = 1.5; du_BB = 1.5; hcrit_SH = 200.0
sigstar_DJ = 0.0; sig0_DJ = 0.1; theta_SH = 8.0

!---diffusion coefficients
hdifmom_cst = 0.0; hdifscal_cst = 0.0; kinvisc_cst = 1.0E-06
rho_crit_iso = 0.01; skewdiff_cst = 1000.0; slopemax_iso = 0.01
smag_coef_mom = 0.1; smag_coef_scal = 0.1
vdifmom_cst = 1.0E-06; vdifscal_cst = 1.0E-06

!---water depths
depmean_cst = 0.0; depmean_flag = 0.0

!---inundation schemes
dcrit_fld = 0.1; dmin_fld = 0.02; dthd_fld = 0.1; ndwexp = 1.0
fld_mask(1) = 1; fld_mask(2:nofldmasks) = 0

!---bottom/surface fluxes
bdragcoef_cst = 0.0; bdraglin = 0.0; ccharno = 0.014; ces_scst = 0.0013
ces_ucst = 0.00113; chs_scst = 0.00113; chs_ucst = 0.00113; ckar = 0.4
zbtoz0lim = 2.0; zref_temp = 10.0; zref_wind = 10.0; zrough_cst = 0.0
cdspars(1) = 0.0; cdspars(2) = 0.00113; cdspars(3:4) = 0.0

!---open boundary conditions
cgravratio = 0.03; distrlx_obc = 0.0
nqsecobu = 0; nqsecobv = 0
return_time = 0.0

!---optical parameters
optattcoef1_cst = 10.0; optattcoef2_cst = 0.067; opt_frac = 0.54

!---numerical
theta_cor = 0.5; theta_sur = 0.5; theta_vadv = 0.501; theta_vdif = 1.0

!---implicit scheme
maxitsimp = 1; nomgiterations = 100; nomglevels = 1;  nopostsweeps = 2
nopresweeps = 2; nosmoothsteps = 2
dzetaresid_conv = 1.0E-07; mg_tol = 1.0E-07; ur_mg = 1.0; ur_smooth = 0.8

!---tidal
index_obc = 0; index_astro = 0

!---structure/discharge parameters
numdis = 0; numdry = 0; numthinu = 0; numthinv = 0
numwbaru = 0; numwbarv = 0
wbarrlxu = 1.0; wbarrlxv = 1.0

!---wave parameters
wave_penetration_bed = 10.0; wave_penetration_surf = 2.0; wavethick_cst = 0.1

!
!6. Turbulence parameters
!------------------------
!
!---stability functions
skeps = 0.09; sq_my = 0.2; c_sk = 0.15

!---limiting conditions
tkelim = 1.0E-06; tkemin = 1.0E-014; dissipmin = 1.0E-012; zlmixmin = 1.7E-010

!---t.k.e.-equation
sigma_k = 1.0; wfltke = 0.0

!---kl-equation
e1_my = 1.8; e2_my = 1.33; e3_my = 1.0

!---eps-equation
c1_eps = 1.44; c2_eps = 1.92; c31_eps = 0.2; c32_eps = 1.0

!---roughness lengths
zrough_bot = 0.0; zrough_sur = 0.0

!---KPP scheme for internal waves
riccrit_iw = 0.7; vdifmom_iw = 1.0E-04; vdifscal_iw = 5.0E-05
vdifshear_iw = 0.005

!---Pacanowski-Philander relations
alpha_pp = 5.0; expmom_pp = 2.0; vbmom_pp = 1.0E-04; vbscal_pp = 1.0E-05
vmax_pp = 3.0; v0dif_pp = 0.01

!---Munk-Anderson relations
alpha_ma = 10.0; beta_ma = 3.33; expmom_ma = 0.5; expscal_ma = 1.5
vmaxmom_ma = 3.0; vmaxscal_ma = 4.0; v0dif_ma = 0.06

!---algebraic flow-dependent relations
delta1_ad = 0.0; delta2_ad = 0.0; r1_ad = 1.0; r2_ad = 1.0; k1_ad = 0.0025
k2_ad = 2.0E-05; lambda_ad = 0.0; cnu_ad = 2.0; omega1_ad = 1.0E-04

!---algebraic mixing length formulations
alpha_Black = 0.2; beta_Xing = 2.0

!
!7. I/O parameters
!-----------------
!

intitle = runtitle; outtitle = runtitle

nosetstsr = 0; novarstsr = 0; nostatstsr = 0
nosetsavr = 0; novarsavr = 0; nostatsavr = 0
nosetsanal = 0; nofreqsanal = 0; novarsanal = 0; nostatsanal = 0

!
!8. Attributes of model files
!----------------------------
!

CALL filepars_init(modfiles)

idesc_810: DO idesc=1,MaxIOTypes
   ifil_811: DO ifil=1,MaxIOFiles

!     ---time regular
      modfiles(idesc,ifil,:)%time_regular = .FALSE.

!     ---form
      SELECT CASE (idesc)
         CASE (io_mppmod,io_1uvsur,io_2uvobc,io_2xyobc,io_3uvobc,io_3xyobc,&
             & io_salobc,io_tmpobc,io_bioobc,io_nstspc,io_sedobc)
            IF (ifil.EQ.1) THEN
               modfiles(idesc,ifil,:)%form = 'A'
            ELSE
               modfiles(idesc,ifil,:)%form = 'N'
            ENDIF
         CASE DEFAULT
            modfiles(idesc,ifil,:)%form = 'N'
      END SELECT

!     ---info
      modfiles(idesc,ifil,:)%info = .FALSE.

!     ---status
      modfiles(idesc,ifil,:)%status = '0'

!     ---filling
      modfiles(idesc,ifil,:)%fill = .FALSE.

!     ---packing
      modfiles(idesc,ifil,:)%packing = .FALSE.      
      
!     ---end file
      modfiles(idesc,ifil,:)%endfile = 0

!     ---unit
      modfiles(idesc,ifil,:)%iunit = int_fill

!     ---time record number
      modfiles(idesc,ifil,:)%timerec = 0

!     ---time skips
      modfiles(idesc,ifil,:)%tskips = 1

   ENDDO ifil_811

!  ---time limits
   iotype_812: DO iotype=1,2
   ifil_812: DO ifil=1,MaxIOFiles
      SELECT CASE (idesc)
      CASE (io_mppmod,io_modgrd,io_metabs,io_sstabs,io_wavabs,io_bioabs,&
          & io_nstabs,io_metrel,io_sstrel,io_wavrel,io_biorel,io_nstrel,&
          & io_biospc,io_drycel,io_thndam,io_weibar,io_disspc)
            modfiles(idesc,ifil,iotype)%tlims = 0
         CASE (io_1uvsur,io_2uvobc,io_2xyobc,io_3uvobc,io_3xyobc,io_salobc,&
             & io_tmpobc,io_sedobc,io_bioobc)
            IF (ifil.EQ.1) THEN
               modfiles(idesc,1,iotype)%tlims = 0
            ELSE
               modfiles(idesc,ifil,iotype)%tlims = (/0,0,1/)
            ENDIF
         CASE DEFAULT
            modfiles(idesc,ifil,iotype)%tlims = (/0,0,1/)
      END SELECT
   ENDDO ifil_812
   ENDDO iotype_812

ENDDO idesc_810

!
!9. Surface grid attributes
!--------------------------
!

CALL gridpars_init(surfacegrids)

!---model grid
surfacegrids(igrd_model,1)%rotated = .FALSE.
surfacegrids(igrd_model,1)%nhtype = 1
surfacegrids(igrd_model,1)%n1dat = nc
surfacegrids(igrd_model,1)%n2dat = nr
surfacegrids(igrd_model,1)%delxdat = 1.0
surfacegrids(igrd_model,1)%delydat = 1.0
surfacegrids(igrd_model,1)%x0dat = 0.0
surfacegrids(igrd_model,1)%y0dat = 0.0

!---data flags
igrd_1010: DO igrd=1,MaxGridTypes
   IF (igrd.EQ.igrd_model.OR.igrd.EQ.igrd_meteo) THEN
      surfacegrids(igrd,1)%datflag = 0
   ELSE
      surfacegrids(igrd,1)%datflag = 2
   ENDIF
ENDDO igrd_1010

!
!11. CF global attributes
!------------------------
!

numglbatts = MERGE(7,8,iopt_CDF.EQ.0)

Conventions_CF = CF_version
history_CF = ''
institution_CF = 'RBINS – Operational Directorate Natural Environment'
source_CF =  'Coherens version '//TRIM(model_version)
comment_CF = ''
references_CF = ''

iglb_110: DO iglb=1,MaxGlbAtts
   glbatts(iglb)%name = ''
   glbatts(iglb)%value = ''
ENDDO iglb_110

CALL log_timer_out()


RETURN

END SUBROUTINE default_mod_params

!========================================================================

SUBROUTINE default_out_files_1d(filepars,file_type)
!************************************************************************
!
! *default_out_files_1d* Default values for output file attributes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)default_model.f90  V2.10.1
!
! Description -  argument has rank 1
!
! Reference -
!
! Calling program -
!
! Module calls - filepars_init
!
!************************************************************************
!
USE datatypes
USE switches
USE syspars
USE datatypes_init, ONLY: filepars_init

!
!*Arguments
!
CHARACTER (LEN=2), INTENT(IN) :: file_type
TYPE (FileParams), INTENT(OUT), DIMENSION(:) :: filepars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*filepars*  DERIVED Attributes of user-output files
!*file_type* CHAR    Type of user output
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize


nsize = SIZE(filepars)
IF (nsize.EQ.0) RETURN

procname(pglev+1) = 'default_out_files_1d'
CALL log_timer_in()

CALL filepars_init(filepars)

n_110: DO n=1,nsize
   filepars(n)%status = '0'
   filepars(n)%floattype = MERGE('S','D',float_type.EQ.real_type)
   filepars(n)%form = MERGE('U','N',iopt_CDF.EQ.0)
   filepars(n)%fill = .TRUE.
   filepars(n)%packing = .FALSE.
   filepars(n)%info = .FALSE.
   filepars(n)%comment = ''
   SELECT CASE (file_type)
      CASE ('TS')
         filepars(n)%title = 'Time series data'
      CASE ('TA')
         filepars(n)%title = 'Time averaged data'
      CASE ('HR')
         filepars(n)%title = 'Harmonic residual data'
      CASE ('HA')
         filepars(n)%title = 'Harmonic amplitude data'
      CASE ('HP')
         filepars(n)%title = 'Harmonic phase data'
      CASE ('HE')
         filepars(n)%title = 'Harmonic elliptic data'
      CASE ('PT')
         filepars(n)%title = 'Particle trajectories'
   END SELECT
   filepars(n)%title = TRIM(filepars(n)%title)//': '//TRIM(runtitle)
ENDDO n_110

CALL log_timer_out()


RETURN

END SUBROUTINE default_out_files_1d

!========================================================================

SUBROUTINE default_out_files_2d(filepars,file_type)
!************************************************************************
!
! *default_out_files_2d* Default values for output file attributes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)default_model.f90  V2.10.1
!
! Description -  argument has rank 2
!
! Reference -
!
! Calling program -
!
! Module calls - filepars_init
!
!************************************************************************
!
USE datatypes
USE switches
USE syspars
USE datatypes_init, ONLY: filepars_init

!
!*Arguments
!
CHARACTER (LEN=2), INTENT(IN) :: file_type
TYPE (FileParams), INTENT(OUT), DIMENSION(:,:) :: filepars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*filepars*  DERIVED Attributes of user-output files
!*file_type* CHAR    Type of user output
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n1, n2
INTEGER, DIMENSION(2) :: ndims


IF (SIZE(filepars).EQ.0) RETURN

procname(pglev+1) = 'default_out_files_2d'
CALL log_timer_in()

CALL filepars_init(filepars)

ndims = SHAPE(filepars)

n2_110: DO n2=1,ndims(2)
n1_110: DO n1=1,ndims(1)
   filepars(n1,n2)%status = '0'
   filepars(n1,n2)%floattype = MERGE('S','D',float_type.EQ.real_type)
   filepars(n1,n2)%form = MERGE('U','N',iopt_CDF.EQ.0)
   filepars(n1,n2)%fill = .TRUE.
   filepars(n1,n2)%packing = .FALSE.
   filepars(n1,n2)%info = .FALSE.
   filepars(n1,n2)%comment = ''
   SELECT CASE (file_type)
      CASE ('TS')
         filepars(n1,n2)%title = 'Time series'
      CASE ('TA')
         filepars(n1,n2)%title = 'Time averages'
      CASE ('HR')
         filepars(n1,n2)%title = 'Harmonic residuals'
      CASE ('HA')
         filepars(n1,n2)%title = 'Harmonic amplitudes'
      CASE ('HP')
         filepars(n1,n2)%title = 'Harmonic phases'
      CASE ('HE')
         filepars(n1,n2)%title = 'Elliptic parameters'
   END SELECT
ENDDO n1_110
ENDDO n2_110

CALL log_timer_out()


RETURN

END SUBROUTINE default_out_files_2d

!========================================================================

SUBROUTINE default_out_gpars(outgpars)
!************************************************************************
!
! *default_out_gpars* Default values for output grid parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)default_model.f90  V2.10.2
!
! Description -  argument has rank 1
!
! Reference -
!
! Calling program -
!
! Module calls - outgpars_init
!
!************************************************************************
!
USE datatypes
USE gridpars
USE switches
USE syspars
USE datatypes_init, ONLY: outgpars_init

!
!*Arguments
!
TYPE (OutGridParams), INTENT(OUT), DIMENSION(:) :: outgpars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*filepars*  DERIVED Attributes of output grid parameters
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: iset, nsize


nsize = SIZE(outgpars)
IF (nsize.EQ.0) RETURN

procname(pglev+1) = 'default_out_gpars'
CALL log_timer_in()

CALL outgpars_init(outgpars)

iset_110: DO iset=1,nsize
   outgpars(iset)%bedvars = .FALSE.
   outgpars(iset)%dredgevars = .FALSE.
   outgpars(iset)%relocatevars = .FALSE.
   outgpars(iset)%gridded = .TRUE.
   outgpars(iset)%time_grid = .FALSE.
   outgpars(iset)%status = 'W'
   outgpars(iset)%enddate = cdatetime_undef
   outgpars(iset)%refdate = cdatetime_undef
   outgpars(iset)%startdate = cdatetime_undef
   outgpars(iset)%nodim = MERGE(2,3,iopt_grid_nodim.EQ.2)
   outgpars(iset)%time_format = 1
   outgpars(iset)%vcoord = 2
   outgpars(iset)%xlims = (/1,nc-1,1/)
   outgpars(iset)%ylims = (/1,nr-1,1/)
   outgpars(iset)%zlims = (/1,nz,1/)
ENDDO iset_110

CALL log_timer_out()


RETURN

END SUBROUTINE default_out_gpars


END MODULE default_model
