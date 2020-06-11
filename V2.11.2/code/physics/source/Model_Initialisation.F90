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
! *Model_Initialisation* Series of routines for model initialisation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Initialisation.F90  V2.11.2
!
! $Date: 2018-07-16 17:34:12 +0200 (Mon, 16 Jul 2018) $
!
! $Revision: 1164 $
!
! Description - 
!
! Reference -
!
! Routines - coherens_start, default_grid, default_phsics, define_phsics,
!            exchange_phsics, initialise_model, initialise_physical_arrays,
!            real_phsics, simulation_start
!
!************************************************************************
!

!========================================================================

SUBROUTINE coherens_start
!************************************************************************
!
! *coherens_start* Procedures at start-up of the program
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Initialisation.F90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - coherens_main
!
! External calls -
!
! Module calls - comms_comm_rank, comms_comm_size, comms_errhandler,
!                comms_initialise, error_abort, error_alloc, mat_initialize
!
!************************************************************************
!
USE iopars
USE paralpars
USE switches
USE syspars
USE timepars
USE comms_MPI, ONLY: comms_comm_rank, comms_comm_size, comms_errhandler, &
                   & comms_initialise
USE error_routines, ONLY: error_abort, error_alloc

#ifdef CDF
   USE netcdf
#endif /*CDF*/

IMPLICIT NONE

#ifdef MPI
   INCLUDE 'mpif.h'
#endif /*MPI*/

!
!*Local variables
!
INTEGER :: iproc
REAL :: test


pglev = pglev + 1
procname(pglev) = 'coherens_start'

!
!1. Initial setup of parameters
!------------------------------
!
!---log files
exitlog = .FALSE.; loglev1 = 0; loglev2 = 0

!---error files
ioerr = 0; maxerrors = MaxErrMesgs

!---warning file
warnflag = .FALSE.

!---monitoring files
monflag = .FALSE.
sedflag = .FALSE.

!---system clock count rate
CALL SYSTEM_CLOCK(COUNT_RATE=npcc_rate,COUNT_MAX=npcc_max)

!---kind parameter for real data
kndrtype = KIND(test)

!---data type for real data
float_type = MERGE(real_type,rlong_type,kndrtype.EQ.kndreal)

!---universal parameters
IF (kndrtype.EQ.kndreal) THEN
   pi = pi_s; halfpi = halfpi_s; twopi = twopi_s
   enap = enap_s; degtorad = degtorad_s; radtodeg = radtodeg_s
ELSE
   pi = pi_d; halfpi = halfpi_d; twopi = twopi_d
   enap = enap_d; degtorad = degtorad_d; radtodeg = radtodeg_d
ENDIF

!
!2. Set parameters and switches from compiler options
!----------------------------------------------------
!
!---parallel model
#ifdef MPI
   parallel_set = .TRUE.
   iopt_MPI = 1
#else
   parallel_set = .FALSE.
   iopt_MPI = 0
   npworld =  1
#endif /*MPI*/

!---netCDF I/O
#ifdef CDF
   iopt_CDF = 1
#else
   iopt_CDF = 0
#endif /*CDF*/

!---MCTlibrary
#ifdef MCT
   iopt_MCT = 1
#else
   iopt_MCT = 0
#endif /*MCT*/

!---verification procedure
#ifdef VERIF
   iopt_verif = 1
#else
   iopt_verif = 0
#endif /*VERIF*/

!
!3. Initialise MPI
!-----------------
!

IF (iopt_MPI.EQ.1) THEN

!---MPI definitions
#ifdef MPI
#define MPI_any_source MPI_ANY_SOURCE
#define MPI_any_tag MPI_ANY_TAG
#define MPI_bsend_overhead MPI_BSEND_OVERHEAD
#define MPI_character MPI_CHARACTER
#define MPI_comm_null MPI_COMM_NULL
#define MPI_comm_world MPI_COMM_WORLD
#define MPI_complex MPI_COMPLEX
#define MPI_double_precision MPI_DOUBLE_PRECISION
#define MPI_integer MPI_INTEGER
#define MPI_logical MPI_LOGICAL   
#define MPI_proc_null MPI_PROC_NULL
#define MPI_real MPI_REAL
#define MPI_undefined MPI_UNDEFINED
   any_source_MPI = MPI_any_source
   any_tag_MPI = MPI_any_tag
   bsend_overhead_MPI = MPI_bsend_overhead
   character_MPI = MPI_character
   comm_null_MPI = MPI_comm_null
   comm_world_MPI = MPI_comm_world
   complex_MPI = MPI_complex
   double_precision_MPI = MPI_double_precision
   integer_MPI = MPI_integer
   logical_MPI = MPI_logical
   real_MPI = MPI_real
   real2_MPI = MPI_2real
   double_precision2_MPI = MPI_2double_precision
   proc_null_MPI = MPI_proc_null
   undefined_MPI = MPI_undefined
   float_MPI = MERGE(real_MPI,double_precision_MPI,kndrtype.EQ.kndreal)
   float2_MPI = MERGE(real2_MPI,double_precision2_MPI,kndrtype.EQ.kndreal)
#endif /*MPI*/
   
!  ---initialise
   CALL comms_initialise

!  ---(provisional) parameters for error output
   errchk = master

!  ---handles for communicators
   comm_model = comm_null_MPI

!  ---set error handler
   CALL comms_errhandler(comm_world_MPI)

!  ---model id
   modelid = 1; modelidcoh = 1

!  ---number of processes
   CALL comms_comm_size(comm_world_MPI,npworld)
   nprocscoh = npworld; nprocspart = 0; nprocswav = 0

!  ---process ranks
   CALL comms_comm_rank(comm_world_MPI,idlocglb)
   iprocglb = idlocglb + 1

!  ---master process
   master = idlocglb.EQ.idmaster
   mastermod = master
   idmastercoh = idmaster; idmastermod = idmaster
   
!  ---process ids
   IF (npworld.GT.1) THEN
      ALLOCATE(idprocsglb(npworld),STAT=errstat)
      IF (master) CALL error_alloc('idprocsglb',1,(/npworld/),kndint)
      idprocsglb = (/(iproc,iproc=0,npworld-1)/)
   ENDIF

ENDIF

!
!4. Serial processing
!--------------------
!

IF (npworld.EQ.1) THEN

!  ---(provisional) parameters for error/log output
   errchk = .TRUE.

!  ---handles for communicators
   comm_model = 0

!  ---model id
   modelid = 1; modelidcoh = 1; modelidwav = 0; modelidpart = 0

!  ---number of processes
   nprocs = 1; nprocscoh = 1; nprocspart = 0; nprocswav = 0
   
!  ---process id's
   idlocglb = 0; iprocglb = 1; idloc = 0; iprocloc = 1
   ALLOCATE (idprocsglb(1),STAT=errstat)
   CALL error_alloc('idprocsglb',1,(/1/),kndint)
   idprocsglb = 0
   ALLOCATE (idprocs(1),STAT=errstat)
   CALL error_alloc('idprocs',1,(/1/),kndint)
   idprocs = 0

!  ---master processes
   master = .TRUE.; mastermod = .TRUE.
   idmastercoh = 0; idmastermod = 0; idmasterpart = 0; idmasterwav = 0

!  ---reset switch
   iopt_MPI = 0

ENDIF

!
!5. Parameters for netCDF
!------------------------
!

#ifdef CDF
   char_NF90 = NF90_char; clobber_NF90 = NF90_clobber
   double_NF90 = NF90_double; fill_NF90 = NF90_fill
   fill_double_NF90 = NF90_fill_double
   fill_int_NF90 = NF90_fill_int; fill_real_NF90 = NF90_fill_real 
   global_NF90 = NF90_global; int_NF90 = NF90_int
   noerr_NF90 = NF90_noerr; nofill_NF90 = NF90_nofill
   nowrite_NF90 = NF90_nowrite; offset_64bit_NF90 = NF90_64bit_offset
   real_NF90 = NF90_real; share_NF90 = NF90_share
   sizehint_default_NF90 = NF90_sizehint_default
   unlimited_NF90 = NF90_unlimited; write_NF90 = NF90_write
#endif /*CDF*/
 
!
!6. Data flags
!-------------
!

#ifdef CDF
   int_fill = fill_int_NF90
   real_fill = fill_real_NF90
   double_fill = fill_double_NF90
#else
   int_fill = fill_int_def
   real_fill = fill_real_def
   double_fill = fill_double_def
#endif /*CDF*/
real_fill_eps = 0.00001*ABS(real_fill)
double_fill_eps = 0.00001*ABS(double_fill)
IF (float_type.EQ.real_type) THEN
   float_fill = real_fill; float_fill_eps = real_fill_eps
ELSE
   float_fill = double_fill; float_fill_eps = double_fill_eps
ENDIF

CALL error_abort('coherens_start',ierrno_comms)

pglev = pglev - 1


RETURN

END SUBROUTINE coherens_start

!========================================================================

SUBROUTINE default_grid
!************************************************************************
!
! *default_grid* Default values for grid arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Initialisation.F90  V2.0
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
USE depths
USE gridpars
USE iopars
USE physpars
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'default_grid'
CALL log_timer_in()

!---water depths
depmeanglb(1:nc-1,1:nr-1) = depmean_cst

CALL log_timer_out()


RETURN

END SUBROUTINE default_grid

!========================================================================

SUBROUTINE default_phsics
!************************************************************************
!
! *default_phsics* Default initial values for physical arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Initialisation.F90  V2.5
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - init_turbulence, mixing_length, water_depths, zlmix_to_dissip
!
! Module calls - 
!
!************************************************************************
!
USE density
USE diffusion
USE grid
USE gridpars
USE iopars
USE physpars
USE switches
USE turbulence
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
INTEGER :: iiopt_turb_lmix, k
REAL, PARAMETER :: tke_cst = 1.0E-06


procname(pglev+1) = 'default_phsics'
CALL log_timer_in()

!
!1. Total water depths
!---------------------
!

CALL water_depths

!
!2. Temperature and salinity
!---------------------------
!

k_210: DO k=1,nz
   WHERE (nodeatc.EQ.1)
      temp(:,:,k) = temp_ref
      sal(:,:,k) = sal_ref
   END WHERE
ENDDO k_210

!
!3. Uniform diffusion coefficients
!---------------------------------
!

IF (iopt_vdif_coef.EQ.1) THEN
   k_310: DO k=1,nz+1
      WHERE (nodeatc(0:ncloc,0:nrloc).GT.0) 
         vdifcoefmom(:,:,k) = vdifmom_cst
      END WHERE
      WHERE (maskatc_int)
         vdifcoefscal(:,:,k) = vdifscal_cst
         vdifcoeftke(:,:,k) = 0.0
      END WHERE
   ENDDO k_310
ENDIF

!
!4. Turbulence arrays
!--------------------
!

IF (iopt_vdif_coef.EQ.3) THEN

!
!4.1 Turbulence energy
!---------------------
!

   k_420: DO k=1,nz+1
      WHERE (nodeatc.EQ.1)
         tke(:,:,k) = tke_cst
      END WHERE
   ENDDO k_420

!
!4.2 Mixing length
!-----------------
!

   CALL init_turbulence
   iiopt_turb_lmix = iopt_turb_lmix
   iopt_turb_lmix = MERGE(4,iiopt_turb_lmix,iiopt_turb_lmix.EQ.0)
   CALL mixing_length
   iopt_turb_lmix = iiopt_turb_lmix

!
!4.3 Dissipation rate
!--------------------
!

   CALL zlmix_to_dissip

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE default_phsics

!========================================================================

SUBROUTINE define_phsics
!************************************************************************
!
! *define_phsics* Define or redefine physical initial conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Initialisation.F90  V2.0
!
! Description - 
!
! Reference -
!
! Calling program - coherens_main, initialise_model
!
! External calls - exchange_phsics, read_phsics, usrdef_phsics
!
! Module calls -
!
!************************************************************************
!
USE iopars
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'define_phsics'
CALL log_timer_in()

IF (modfiles(io_inicon,ics_phys,1)%status.EQ.'R') THEN
   CALL read_phsics
ELSEIF (modfiles(io_inicon,ics_phys,1)%status.EQ.'N') THEN
   CALL usrdef_phsics
ENDIF

IF (iopt_MPI.EQ.1) CALL exchange_phsics

CALL log_timer_out()


RETURN

END SUBROUTINE define_phsics

!========================================================================

SUBROUTINE exchange_phsics
!************************************************************************
!
! *exchange_phsics* Exchange initial conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Initialisation.F90  V2.0
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - exchange_mod
!
!************************************************************************
!
USE currents
USE density
USE depths
USE fluxes
USE gridpars
USE iopars
USE modids
USE switches
USE turbulence
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhexch


procname(pglev+1) = 'exchange_phsics'
CALL log_timer_in()

!
!1. 2-D mode
!-----------
!

lbounds(1:2) = 1-nhalo; nhexch = nh2vel
CALL exchange_mod(udvel,lbounds(1:2),nhexch,iarr_udvel)
CALL exchange_mod(vdvel,lbounds(1:2),nhexch,iarr_vdvel)
lbounds(1:2) = 0; nhexch = 1
CALL exchange_mod(zeta,lbounds(1:2),nhexch,iarr_zeta)

!
!2. 3-D currents
!---------------
!

IF ((iopt_grid_nodim.EQ.1).OR.(iopt_grid_nodim.EQ.3)) THEN
   lbounds = (/1-nhalo,1-nhalo,1/); nhexch = nh3vel
   CALL exchange_mod(uvel,lbounds,nhexch,iarr_uvel)
   CALL exchange_mod(vvel,lbounds,nhexch,iarr_vvel)
ENDIF

IF (iopt_grid_nodim.EQ.3) THEN
   lbounds = (/0,0,1/); nhexch = (/1,0,1,0/)
   CALL exchange_mod(wvel,lbounds,nhexch,iarr_wvel)
ENDIF

!
!3. Density arrays
!-----------------
!

lbounds = (/1-nhalo,1-nhalo,1/); nhexch = nhdens
IF (iopt_temp.GT.0) CALL exchange_mod(temp,lbounds,nhexch,iarr_temp)
IF (iopt_sal.GT.0) CALL exchange_mod(sal,lbounds,nhexch,iarr_sal)

!
!4. Turbulence arrays
!--------------------
!

IF (iopt_vdif_coef.EQ.3) THEN
   lbounds = (/1-nhalo,1-nhalo,1/); nhexch = nhturb
   CALL exchange_mod(tke,lbounds,nhexch,iarr_tke)
   IF ((iopt_turb_ntrans.EQ.2).AND.(iopt_turb_param.EQ.1)) THEN
      CALL exchange_mod(zlmix,lbounds,nhexch,iarr_zlmix)
   ELSEIF ((iopt_turb_ntrans.EQ.2).AND.(iopt_turb_param.EQ.2)) THEN
      CALL exchange_mod(dissip,lbounds,nhexch,iarr_dissip)
   ENDIF
ENDIF

!
!5. Bottom stress arrays
!-----------------------
!

IF (iopt_bstres_form.EQ.2) THEN
   lbounds(1:2) = 0; nhexch = (/1,0,1,0/)
   IF (iopt_bstres_drag.LE.2) THEN
      CALL exchange_mod(bdragcoefatc,lbounds(1:2),nhexch,iarr_bdragcoefatc)
   ELSE
      CALL exchange_mod(zroughatc,lbounds(1:2),nhexch,iarr_zroughatc)
   ENDIF
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE exchange_phsics

!========================================================================

SUBROUTINE initialise_model
!************************************************************************
!
! *initialise_model* Initialise model
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Initialisation.F90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - allocate_bio_arrays, allocate_dar_arrays,
!   allocate_global_grid, allocate_mod_arrays, allocate_morph_arrays,
!   allocate_part_arrays, *llocate_sed_arrays, allocate_struc_arrays,
!   astronomical_tides, baroclinic_gradient, biological_model, bottom_stress,
!   coh2wav_coupler, coh2wav_grid, coh2wav_params, coh2wav_switches,
!   combine_particle_grid, combine_particle_phys_data, correct_free_surface,
!   current_corr, current_pred, current_2d, dar_init, default_grid,
!   default_phsics, define_dischr_spc, define_global_grid, define_local_grid,
!   define_partics, define_phsics, define_rlxobc_spec, depth_at_nodes,
!   discharges, domain_decomposition, drying_factor, equation_of_state,
!   exchange_bioics, exchange_sedics, grid_arrays, harmonic_analysis, heat_flux,
!   horizontal_diff_coefs, initialise_biology_arrays, initialise_multigrid,
!   initialise_physical_arrays, initialise_sediment_arrays, local_struc_arrays,
!   meteo_input, morphology_equation, nest_locations, open_boundary_arrays,
!   pointer_arrays, particle_concentrations, particle_trajects,
!   particle_trajects_init, read_bioics, read_bio_spec, read_cif_params,
!   read_dar_spec, read_grid, read_morphics, read_part_spec, read_sedics,
!   read_sed_spec, read_thin_dams, read_weirs, relaxation_weights,
!   salinity_equation, sediment_equation, set_procnums, surface_stress,
!   temperature_equation, time_averages, time_series, update_nest_data_2d,
!   update_wave_params, update_1dsur_data, usrdef_bioics, usrdef_bio_params,
!   usrdef_bio_spec, usrdef_dar_params, usrdef_dar_spec, usrdef_grid,
!   usrdef_mod_params, usrdef_morphics, usrdef_morph_params, usrdef_output,
!   usrdef_part_output, usrdef_part_params, usrdef_part_spec,
!   usrdef_sed_params, usrdef_sedics, usrdef_sed_spec, usrdef_thin_dams,
!   usrdef_weirs, vertical_diff_coefs, vertical_diff_rot, water_depths,
!   wave_input, write_bioics, write_bio_spec, write_cif_vars_anal,
!   write_cif_vars_avrgd, write_cif_vars_bio, write_cif_vars_dar,
!   write_cif_vars_harm, write_cif_vars_mod, write_cif_vars_mon,
!   write_cif_vars_morph, write_cif_vars_part, write_cif_vars_sed,
!   write_cif_vars_tsout, write_cif_vars_tspart, write_dar_spec, write_grid,
!   write_morphics, write_particle_grid, write_particle_phys_data,
!   write_partics, write_part_spec, write_phsics, write_sedics, write_sed_spec,
!   write_thin_dams, write_weirs
!
! Module calls - check_bioics, check_bio_params, check_dar_params,
!   check_dar_spec, check_grid_arrays, check_initial_conditions,
!   check_mod_params, check_morphics, check_morph_params, check_partics,
!   check_partition, check_part_params, check_sed_params, check_sedics,
!   check_structs_locs, close_filepars, default_bio_params, default_bioics,
!   default_dar_params, default_mod_params, default_morphics,
!   default_morph_params, default_part_params, default_sed_params,
!   default_sed_spec, default_sedics, initialise_time, open_file,
!   open_filepars, reset_bio_params, reset_bioics, reset_dar_params,
!   reset_initial_conditions, reset_mod_params, reset_morph_params,
!   reset_part_params, reset_partics, reset_sed_params, read_sedics, rng_init,
!   sum_vars
!
!************************************************************************
!
USE iopars
USE grid
USE gridpars
USE obconds
USE paralpars
USE structures
USE switches
USE syspars
USE timepars
USE check_biology, ONLY: check_bioics, check_bio_params 
USE check_model, ONLY: check_grid_arrays, check_initial_conditions, &
     & check_mod_params, check_partition, check_struct_locs
USE check_particles, ONLY: check_part_params, check_partics 
USE check_sediments, ONLY: check_dar_params, check_dar_spec, check_morphics, &
                         & check_morph_params, check_sedics, check_sed_params
USE default_biology, ONLY: default_bioics, default_bio_params
USE default_model, ONLY: default_mod_params
USE default_particles, ONLY: default_part_params
USE default_sediments, ONLY: default_dar_params, default_morphics, &
                           & default_morph_params, default_sedics, &
                           & default_sed_params, default_sed_spec
USE inout_routines, ONLY: close_filepars, open_file, open_filepars
USE paral_utilities, ONLY: sum_vars
USE reset_biology, ONLY: reset_bioics, reset_bio_params
USE reset_model, ONLY: reset_initial_conditions, reset_mod_params
USE reset_particles, ONLY: reset_part_params, reset_partics
USE reset_sediments, ONLY: reset_dar_params, reset_morph_params, reset_sedics, &
                         & reset_sed_params
USE rng_library, ONLY: rng_init
USE time_routines, ONLY: initialise_time, log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc


procname(pglev+1) = 'initialise_model'
CALL log_timer_in(npcc)

!
!1. Physical model parameters
!----------------------------
!
!1.1 Defaults
!------------
!

CALL default_mod_params

!
!1.2 Define
!----------
!

IF (ciffile%status.EQ.'R') THEN
   CALL read_cif_params(icif_mod)
ELSE
   CALL usrdef_mod_params
ENDIF

!
!1.3 Reset and check
!-------------------
!

CALL reset_mod_params
CALL check_mod_params

!
!1.4 Monitoring files
!--------------------
!
!---implicit algorithm
monflag = monlog.GT.0.AND.master.AND.iopt_hydro_impl.EQ.1
IF (.NOT.monflag) monlog = 0
IF (monflag) CALL open_file(iomon,monlog_file,'OUT','A')

!---sediment monitoring
sedflag = sedlog.AND.master.AND.iopt_sed.GT.0
IF (sedflag) CALL open_file(iosed,sedlog_file,'OUT','A')

!
!1.5 Process numbers
!-------------------
!

CALL set_procnums

!
!1.6 Initialise
!--------------
!
!---time
CALL initialise_time
!---random generator
CALL rng_init

!
!2. Sediment parameters
!----------------------
!
!2.1 Water column sediment
!-------------------------
!

IF (iopt_sed.GT.0) THEN

!  ---defaults
   CALL default_sed_params

!  ---define
   IF (ciffile%status.EQ.'R') THEN
      CALL read_cif_params(icif_sed)
   ELSE 
      CALL usrdef_sed_params
   ENDIF

!  ---reset and check
   CALL reset_sed_params
   CALL check_sed_params

ENDIF

!
!2.2 Morphology
!--------------
!

IF (iopt_morph.GT.0) THEN

!  ---defaults
   CALL default_morph_params

!  ---define
   IF (ciffile%status.EQ.'R') THEN
      CALL read_cif_params(icif_morph)
   ELSE 
      CALL usrdef_morph_params
   ENDIF

!  ---reset and check
   CALL reset_morph_params
   CALL check_morph_params

ENDIF

!
!2.3. Dredging/relocation parameters
!------------------------
!

IF (iopt_dar.GT.0) THEN

!  ---defaults
   CALL default_dar_params

!  ---define
   IF  (ciffile%status.EQ.'R') THEN
      CALL read_cif_params(icif_dar)
   ELSE 
      CALL usrdef_dar_params
   ENDIF

!  ---reset and check
   CALL reset_dar_params
   CALL check_dar_params

ENDIF

!
!3. Biological parameters
!------------------------
!

IF (iopt_biolgy.GT.0) THEN

!  ---defaults
   CALL default_bio_params

!  ---define
   IF (ciffile%status.EQ.'R') THEN
      CALL read_cif_params(icif_bio)
   ELSE
      CALL usrdef_bio_params
   ENDIF

!  --- reset and check
   CALL reset_bio_params
   CALL check_bio_params

ENDIF

!
!4. Parameters for particle model
!--------------------------------
!

IF (iopt_part_model.GT.0.OR.iopt_part_write.EQ.1) THEN

!  ---defaults
   CALL default_part_params

!  ---define
   IF (ciffile%status.EQ.'R') THEN
      CALL read_cif_params(icif_part)
   ELSE
      CALL usrdef_part_params
   ENDIF

!  ---reset and check
   CALL reset_part_params
   CALL check_part_params
   
ENDIF
   
!
!5. Model grid
!-------------
!
!5.1 Grid coordinates and bathymetry
!-----------------------------------
!
!---allocate global arrays
CALL allocate_global_grid

!---allocate structures arrays
CALL allocate_struc_arrays

!---defaults
CALL default_grid

!---global grid
IF (modfiles(io_modgrd,1,1)%status.EQ.'R') THEN
   CALL read_grid
ELSEIF (modfiles(io_modgrd,1,1)%status.EQ.'N') THEN
   CALL usrdef_grid
ENDIF
CALL define_global_grid

!---check
CALL check_grid_arrays

!---write
IF (modfiles(io_modgrd,1,2)%status.EQ.'W') CALL write_grid

!---domain decomposition
CALL domain_decomposition

!---allocate model arrays
CALL allocate_mod_arrays
IF (iopt_sed.GT.0) CALL allocate_sed_arrays
IF (iopt_morph.GT.0) CALL allocate_morph_arrays
IF (iopt_dar.GT.0) CALL allocate_dar_arrays
IF (iopt_biolgy.GT.0) CALL allocate_bio_arrays
IF (iopt_part_model.GT.0) CALL allocate_part_arrays

!---distribute coordinate arrays and bathymetry onto local processes
CALL define_local_grid

!---check partition
IF (master.AND.iopt_MPI.EQ.1) THEN
   CALL check_partition
ENDIF

!
!5.2 Structure parameters
!------------------------
!
!---thin dams
IF (iopt_thndam.EQ.1) THEN
   IF (modfiles(io_thndam,1,1)%status.EQ.'R') THEN
      CALL read_thin_dams
   ELSEIF (modfiles(io_thndam,1,1)%status.EQ.'N') THEN
      CALL usrdef_thin_dams
   ENDIF
   IF (modfiles(io_thndam,1,2)%status.EQ.'W') CALL write_thin_dams
ENDIF

!---weirs/barriers
IF (iopt_weibar.EQ.1) THEN
   IF (modfiles(io_weibar,1,1)%status.EQ.'R') THEN
      CALL read_weirs
   ELSEIF (modfiles(io_weibar,1,1)%status.EQ.'N') THEN
      CALL usrdef_weirs
   ENDIF
   IF (modfiles(io_weibar,1,2)%status.EQ.'W') CALL write_weirs
ENDIF

!---check locations of structures
CALL check_struct_locs

!---definition of local arrays
CALL local_struc_arrays

!
!5.3 Define arrays related to model grid
!---------------------------------------
!
!---pointer arrays
CALL pointer_arrays

!---arrays for open boundary conditions
CALL open_boundary_arrays

!---grid arrays
CALL grid_arrays

!---mean water depths at nodes
CALL depth_at_nodes

!---send grid data to particle model
IF (iopt_part_write.EQ.1) CALL write_particle_grid
IF (iopt_part_model.EQ.2) CALL combine_particle_grid

!
!5.4 Number of active points
!---------------------------
!

IF (iopt_fld.GT.0) THEN

!  ---C-nodes
   nowetatcloc = COUNT(maskatc_int)
   CALL sum_vars(nowetatcloc,nowetatc,0,commall=.TRUE.)

!  ---U-nodes
   nowetatuloc =  COUNT(node2du(1:ncloc,1:nrloc).GT.1)
   CALL sum_vars(nowetatuloc,nowetatu,0,commall=.TRUE.)

!  ---V-nodes
   nowetatvloc =  COUNT(node2dv(1:ncloc,1:nrloc).GT.1)
   CALL sum_vars(nowetatvloc,nowetatv,0,commall=.TRUE.)

ENDIF

!
!6. Discharge locations
!----------------------
!

IF (iopt_dischr.EQ.1) THEN
   CALL define_dischr_spec
   CALL discharges
ENDIF

!
!7. Initial conditions
!---------------------
!
!7.1 Physics
!-----------
!
!---defaults
CALL default_phsics

!---define
IF (physinit) CALL define_phsics

!---reset and check
CALL reset_initial_conditions
CALL check_initial_conditions

!---total water depths
CALL water_depths

!---"alpha"-factor for flooding/drying scheme
IF (iopt_fld_alpha.GT.0) CALL drying_factor

!---initialise
CALL initialise_physical_arrays

!
!7.2 Sediment
!------------
!
!7.2.1 Water column sediment
!---------------------------
!

IF (iopt_sed.GT.0) THEN
!  ---set defaults
   CALL default_sedics

!  ---specifier arrays
   CALL default_sed_spec
   IF (modfiles(io_sedspc,1,1)%status.EQ.'R') THEN
      CALL read_sed_spec
   ELSEIF (modfiles(io_sedspc,1,1)%status.EQ.'N') THEN
      CALL usrdef_sed_spec
      IF (modfiles(io_sedspc,1,2)%status.EQ.'W') CALL write_sed_spec
   ENDIF

!  ---initial conditions
   IF (modfiles(io_inicon,ics_sed,1)%status.EQ.'R') THEN
      CALL read_sedics
   ELSEIF (modfiles(io_inicon,ics_sed,1)%status.EQ.'N') THEN
      CALL usrdef_sedics
      IF (iopt_MPI.EQ.1) CALL exchange_sedics
   ELSE
      IF (warnflag) THEN
         WRITE (iowarn,'(A)') &
                    & 'WARNING: using default initial conditions for sediments'
      ENDIF
   ENDIF

!  ---reset and check
   CALL reset_sedics
   CALL check_sedics

!  ---initialise sediment arrays
   CALL initialise_sediment_arrays

ENDIF

!
!7.2.2 Morphology
!----------------
!

IF (iopt_morph.GT.0) THEN

!  ---set defaults
   CALL default_morphics

!  ---initial conditions
   IF (modfiles(io_inicon,ics_morph,1)%status.EQ.'R') THEN
      CALL read_morphics
   ELSEIF (modfiles(io_inicon,ics_morph,1)%status.EQ.'N') THEN
      CALL usrdef_morphics
   ELSE
      IF (warnflag) THEN
         WRITE (iowarn,'(A)') &
                   & 'WARNING: using default initial conditions for morphology'
      ENDIF
   ENDIF

!  ---check
   CALL check_morphics

ENDIF

!
!7.2.3 Dredging and relocation
!-----------------------------
!

IF (iopt_dar.GT.0) THEN

!  ---specifier arrays
   IF (modfiles(io_darspc,1,1)%status.EQ.'R') THEN
      CALL read_dar_spec
   ELSEIF (modfiles(io_darspc,1,1)%status.EQ.'N') THEN
      CALL usrdef_dar_spec
      IF (modfiles(io_darspc,1,2)%status.EQ.'W') CALL write_dar_spec
   ENDIF

!  --- check
   CALL check_dar_spec

ENDIF

!
!7.3 Biology
!-----------
!

IF (iopt_biolgy.GT.0) THEN

!  ---set defaults
   CALL default_bioics

!  ---initial conditions
   IF (modfiles(io_inicon,ics_bio,1)%status.EQ.'R') THEN
      CALL read_bioics
   ELSEIF (modfiles(io_inicon,ics_bio,1)%status.EQ.'N') THEN
      CALL usrdef_bioics
      IF (iopt_MPI.EQ.1) CALL exchange_bioics
   ENDIF

!  ---specifier arrays
   IF (modfiles(io_biospc,1,1)%status.EQ.'R') THEN
      CALL read_bio_spec
   ELSEIF (modfiles(io_biospc,1,1)%status.EQ.'N') THEN
      CALL usrdef_bio_spec
      IF (modfiles(io_biospc,1,2)%status.EQ.'W') CALL write_bio_spec
   ENDIF

!  ---reset and check
   CALL reset_bioics
   CALL check_bioics

!  ---initialise biology arrays
   CALL initialise_biology_arrays

ENDIF

!
!7.4 Particle model
!------------------
!

IF (iopt_part_model.EQ.1) THEN

!  ---specifier arrays
   IF (modfiles(io_parspc,1,1)%status.EQ.'R') THEN
      CALL read_part_spec
   ELSEIF (modfiles(io_parspc,1,1)%status.EQ.'N') THEN
      CALL usrdef_part_spec
   ENDIF
   IF (modfiles(io_parspc,1,2)%status.EQ.'W') CALL write_part_spec
   
!  ---intialise physical arrays
   CALL define_particle_phys_data
   
!  ---particles attributes
   CALL define_partics

!  ---reset and chzeck   
   CALL reset_partics
   CALL check_partics

ENDIF
   
!
!8. Density arrays
!-----------------
! 
!---density and expansion coefficients
CALL equation_of_state

!---baroclinic pressure gradient
IF (iopt_mode_3D.EQ.2.AND.iopt_dens_grad.GT.0) CALL baroclinic_gradient

!
!9. Astronomical tide
!---------------------
!

IF (iopt_astro_tide.EQ.1) CALL astronomical_tides

!
!10. Surface/bottom forcing
!--------------------------
!
!---meteo data
IF (iopt_meteo.EQ.1) CALL meteo_input

!---initialise coupler with wave model 
IF (.NOT.cold_start.AND.iopt_waves_couple.GT.0)  THEN
!     --exchange model parameters
      CALL coh2wav_params
!     --exchange model grids
      CALL coh2wav_grid
!     --exchange switches  
      CALL coh2wav_switches
!     --initialise coupler and interpolation procedures
      CALL coh2wav_coupler
ENDIF

!---read or receive wave data
IF (iopt_waves.EQ.1.AND.(.NOT.cold_start.OR.iopt_waves_couple.EQ.0)) THEN
   CALL wave_input
ENDIF

!---surface/bottom stresses
IF (iopt_bstres_form.GT.0) CALL bottom_stress
IF (iopt_meteo_stres.EQ.1) CALL surface_stress

!---wave
IF (iopt_waves.EQ.1.AND.(.NOT.cold_start.OR.iopt_waves_couple.EQ.0)) THEN
   CALL update_wave_params
ENDIF

!
!11. Diffusion coefficients
!--------------------------
!
!---vertical diffusion coefficients
IF (iopt_vdif_coef.GT.0) CALL vertical_diff_coefs
!---horizontal diffusion coefficients
IF (iopt_hdif_coef.GT.0) CALL horizontal_diff_coefs
!---vertical rotated diffusion
IF (iopt_vdif_rot.EQ.1) CALL vertical_diff_rot

!
!12. Send physical data to the particle model
!--------------------------------------------
!

IF (iopt_part_write.EQ.1) CALL write_particle_phys_data(CStartDateTime)
IF (iopt_part_model.EQ.2) CALL combine_particle_phys_data

!
!13. Parameters for relaxation scheme
!------------------------------------
!

IF (iopt_obc_relax.EQ.1) THEN
   CALL define_rlxobc_spec
   CALL relaxation_weights
ENDIF

!
!14. Locations for nested output
!-------------------------------
!

IF (iopt_nests.EQ.1) CALL nest_locations

!
!15. Initialise forcing and open boundary conditions
!---------------------------------------------------
!
!---1-D mode
IF (iopt_sur_1D.EQ.1) CALL update_1dsur_data

!---3-D mode (predictor)
IF (iopt_grid_nodim.NE.2.AND.iopt_mode_3D.EQ.2) THEN
   CALL current_pred
ENDIF

!---2-D mode
IF (iopt_hydro_impl.EQ.0.AND.iopt_grid_nodim.NE.1) THEN
   CALL current_2d
ENDIF

!---initialise multigrid scheme
IF (iopt_hydro_impl.EQ.1) THEN
   CALL correct_free_surf
   CALL initialise_multigrid 
ENDIF

IF (iopt_nests.EQ.1) THEN
   CALL update_nest_data_2d(io_2uvnst)
   CALL update_nest_data_2d(io_2xynst)
ENDIF

!---3-D mode (corrector)
IF (iopt_grid_nodim.EQ.3.AND.iopt_mode_3D.EQ.2) CALL current_corr

!---density
IF (iopt_temp.EQ.2) CALL temperature_equation
IF (iopt_sal.EQ.2) CALL salinity_equation

!---sediment
IF (iopt_sed.GT.0) CALL sediment_equation

!---morphology
IF (iopt_morph.GT.0) CALL morphology_equation

!---dredging and relocation
IF (iopt_dar.GT.0) CALL dar_init

!---biology
IF (iopt_biolgy.GT.0) CALL biological_model

!
!16. Initialise model output
!---------------------------
!
!---restart files
CALL write_phsics
IF (iopt_sed.GT.0) CALL write_sedics
IF (iopt_morph.GT.0) CALL write_morphics
IF (iopt_biolgy.GT.0) CALL write_bioics
IF (iopt_part_model.EQ.1) CALL write_partics

!---user-defined
IF (iopt_out_tsers.GT.0) CALL time_series
IF (iopt_out_avrgd.GT.0) CALL time_averages
IF (iopt_out_anal.GT.0) CALL harmonic_analysis
IF (iopt_part_model.GT.0) CALL particle_concentrations
IF (iopt_part_model.EQ.1) CALL particle_trajects
CALL usrdef_output
IF (iopt_part_model.EQ.1) CALL usrdef_part_output

!
!17. CIF
!-------
!
!17.1  Close CIF for reading
!---------------------------
!

IF (ciffile%status.EQ.'R') CALL close_filepars(ciffile)

!
!17.2 Write CIF
!--------------
!

IF (master.AND.ciffile%status.EQ.'W') THEN
!  ---open file   
   CALL open_filepars(ciffile)
!  ---general/monitoring parameters 
   CALL write_cif_vars_mon
!  ---model setup
   CALL write_cif_vars_mod
!  ---biological model parameters 
   IF (iopt_biolgy.GT.0) CALL write_cif_vars_bio
!  ---sediment model parameters 
   IF (iopt_sed.GT.0) CALL write_cif_vars_sed
!  ---morphological model parameters 
   IF (iopt_morph.GT.0) CALL write_cif_vars_morph
!  ---dredging/relocation model parameters 
   IF (iopt_dar.GT.0) CALL write_cif_vars_dar
!  ---particle model parameters 
   IF (iopt_part_model.GT.0) CALL write_cif_vars_part
!  ---time series output 
   IF (iopt_out_tsers.GT.0) CALL write_cif_vars_tsout
!  ---time averaged output 
   IF (iopt_out_avrgd.GT.0) CALL write_cif_vars_avrgd
!  ---harmonic analysis
   IF (iopt_out_anal.GT.0) CALL write_cif_vars_harm
!  ---harmonic output
   IF (iopt_out_anal.GT.0) CALL write_cif_vars_anal
!  ---particle output
   IF (iopt_part_model.EQ.2) CALL particle_trajects_init
   IF (iopt_part_model.GT.0) CALL write_cif_vars_tspart
!  ---close CIF   
   CALL close_filepars(ciffile)
ENDIF

!
!18. Deallocate global grid arrays
!---------------------------------
!

DEALLOCATE (gscoordglb)

CALL log_timer_out(npcc,itm_init)


RETURN

END SUBROUTINE initialise_model

!========================================================================

SUBROUTINE initialise_physical_arrays
!************************************************************************
!
! *initialise_physical_arrays* Initialise arrays for main program
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Initialisation.F90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - physical_vertical_current, store_depths_old,
!                  transf_vertical_current
!
! Module calls - Carr_at_U, Carr_at_V, error_alloc, exchange_mod
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
USE meteo
USE modids
USE optics
USE physpars
USE switches
USE syspars
USE array_interp, ONLY: Carr_at_U, Carr_at_V
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: k
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu, maskatv


procname(pglev+1) = 'initialise_physical_arrays'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (maskatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc/),kndlog)

!
!2. Mask arrays
!--------------
!

maskatu = node2du(1:ncloc,1:nrloc).GT.1
maskatv = node2dv(1:ncloc,1:nrloc).GT.1

!
!3. Currents
!-----------
!
!3.1 Depth-mean
!--------------
!
!---U-nodes
WHERE (maskatu)
   umvel(1:ncloc,1:nrloc) = udvel(1:ncloc,1:nrloc)/deptotatu(1:ncloc,1:nrloc)
END WHERE
!---V-nodes
WHERE (maskatv)
   vmvel(1:ncloc,1:nrloc) = vdvel(1:ncloc,1:nrloc)/deptotatv(1:ncloc,1:nrloc)
END WHERE
!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   lbounds(1:2)  = (/1-nhalo,0/); nhexch = (/nh2vel,nh2vel,1,1/)
   CALL exchange_mod(umvel,lbounds(1:2),nhexch,iarr_umvel)
   lbounds(1:2)  = (/0,1-nhalo/); nhexch = (/1,1,nh2vel,nh2vel/)
   CALL exchange_mod(vmvel,lbounds(1:2),nhexch,iarr_vmvel)
ENDIF

!---predicted currents
umpred = umvel(0:ncloc+1,1:nrloc); vmpred = vmvel(1:ncloc,0:nrloc+1)

!
!3.2 3-D currents for 2-D case
!-----------------------------
!

IF (iopt_grid_nodim.EQ.2) THEN
!  ---U-nodes
   k_321: DO k=1,nz
      uvel(1:ncloc,1:nrloc,k) = umvel(1:ncloc,1:nrloc)
   ENDDO k_321
!  ---V-nodes
   k_322: DO k=1,nz
      vvel(1:ncloc,1:nrloc,k) = vmvel(1:ncloc,1:nrloc)
   ENDDO k_322
!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = (/1-nhalo,1-nhalo,1/); nhexch = nh3vel
      CALL exchange_mod(uvel,lbounds,nhexch,iarr_uvel)
      CALL exchange_mod(vvel,lbounds,nhexch,iarr_vvel)
   ENDIF
ENDIF

!
!3.3 Filtered currents
!--------------------
!

IF (iopt_curr.EQ.1.OR.iopt_grid_nodim.NE.3.OR.&
 & (iopt_grid_nodim.EQ.3.AND.iopt_hydro_impl.EQ.1)) THEN
!  ---U-nodes
   k_331: DO k=1,nz
      WHERE (maskatu)
         ufvel(1:ncloc,:,k) = uvel(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_331
!  ---V-nodes
   k_332: DO k=1,nz
      WHERE (maskatv)
         vfvel(:,1:nrloc,k) = vvel(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_332
!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = (/2-nhalo,1,1/); nhexch = (/nhfvel-1,nhfvel,0,0/)
      CALL exchange_mod(ufvel,lbounds,nhexch,iarr_ufvel)
      lbounds = (/1,2-nhalo,1/); nhexch = (/0,0,nhfvel-1,nhfvel/)
      CALL exchange_mod(vfvel,lbounds,nhexch,iarr_vfvel)
   ENDIF
ENDIF

!
!3.4 Store old values
!--------------------
!

uvel_old = uvel; vvel_old = vvel

!
!3.5 Vertical currents
!---------------------
!

IF (iopt_mode_3D.EQ.1.AND.iopt_grid_nodim.EQ.3) THEN
!  ---transformed
   CALL transf_vertical_current
!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = (/0,0,1/); nhexch = (/1,0,1,0/)
      CALL exchange_mod(wvel,lbounds,nhexch,iarr_wvel)
   ENDIF
!  ---physical
   CALL physical_vertical_current
ENDIF

!
!4. Roughness at velocity nodes 
!------------------------------
!

IF (iopt_bstres_drag.GT.2) THEN
   CALL Carr_at_U(zroughatc(:,1:nrloc),zroughatu,1,1,(/0,1,nz/),&
               & (/ncloc,nrloc,nz/),1,iarr_zroughatc,.TRUE.)
   CALL Carr_at_V(zroughatc(1:ncloc,:),zroughatv,1,1,(/1,0,nz/),&
               & (/ncloc,nrloc,nz/),1,iarr_zroughatc,.TRUE.)
ENDIF

!
!5. Store old depth arrays
!-------------------------
!

CALL store_depths_old

!
!6. Apply freezing point limit
!-----------------------------
!

IF (iopt_temp.GT.0) THEN

!  ---imposed lower limit
   IF (ABS(temp_min-float_fill).GT.float_fill_eps) THEN
      k_610: DO k=1,nz
         WHERE (maskatc_int)
            temp(1:ncloc,1:nrloc,k) = MAX(temp_min,temp(1:ncloc,1:nrloc,k))
         END WHERE
      ENDDO k_610
      
!  ---freezing point limit
   ELSE
      k_620: DO k=1,nz
         WHERE (maskatc_int)
            temp(1:ncloc,1:nrloc,k) = MAX(-0.0575*sal(1:ncloc,1:nrloc,k),&
                 & temp(1:ncloc,1:nrloc,k))
         END WHERE
      ENDDO k_620
   ENDIF

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = (/1-nhalo,1-nhalo,1/); nhexch = nhdens
      CALL exchange_mod(temp,lbounds,nhexch,iarr_temp)
   ENDIF

ENDIF

!
!7. Other arrays
!---------------
!

WHERE (maskatc_int)
   optattcoef2 = optattcoef2_cst
END WHERE
WHERE (nodeatc(0:ncloc,0:nrloc).EQ.1)
   atmpres = atmpres_ref
END WHERE

!
!8. Deallocate arrays
!---------------------
!

DEALLOCATE (maskatu,maskatv)

CALL log_timer_out()


RETURN

END SUBROUTINE initialise_physical_arrays

!========================================================================

SUBROUTINE read_phsics
!************************************************************************
!
! *read_phsics* Read physical initial conditions in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Initialisation.F90  V2.10.1
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - close_filepars, error_alloc_struc, open_filepars, read_time,
!                read_varatts_mod, read_vars, read_distribute_mod,
!                read_glbatts_mod, varatts_init
!
!************************************************************************
!
USE currents
USE datatypes
USE density
USE depths
USE fluxes
USE gridpars
USE iopars
USE obconds
USE switches
USE structures
USE syspars
USE tide
USE timepars
USE turbulence
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_alloc_struc
USE inout_paral, ONLY: read_distribute_mod
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_time, read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: fcomm
CHARACTER (LEN=lentime) :: ciodatetime
INTEGER :: numvars, varid
TYPE(FileParams) :: filepars
INTEGER, DIMENSION(3) :: lbounds, ubounds
INTEGER, DIMENSION(4) :: nhdims
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_phsics'
CALL log_timer_in()

filepars = modfiles(io_inicon,ics_phys,1)

!
!1. File header
!--------------
!

CALL open_filepars(filepars)
CALL read_glbatts_mod(filepars)
numvars = filepars%novars + 1
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL varatts_init(varatts)
CALL read_varatts_mod(filepars,varatts,numvars)

!
!2. Read physical data
!---------------------
!

fcomm = .FALSE.

DO WHILE (.NOT.fcomm)

!
!2.1 Date/time
!-------------
!

   CALL read_time(ciodatetime,filepars)
   fcomm = ciodatetime.EQ.CStartDateTime
   varid = 1

!
!2.2 2-D mode
!------------
!

   IF (iopt_mode_2D.GT.0) THEN

!     ---X-current
      varid = varid + 1; nhdims = nhalo
      lbounds(1:2) = (/1-nhalo,1-nhalo/); ubounds(1:2) = (/nc+nhalo,nr+nhalo/)
      CALL read_distribute_mod(udvel,filepars,varid,lbounds(1:2),ubounds(1:2),&
                             & nhdims,'U  ',fcomm,varatts=varatts)

!     ---Y-current
      varid = varid + 1
      CALL read_distribute_mod(vdvel,filepars,varid,lbounds(1:2),ubounds(1:2),&
                             & nhdims,'V  ',fcomm,varatts=varatts)

!     ---surface elevation
      varid = varid + 1; nhdims = 1
      lbounds(1:2) = 0; ubounds(1:2) = (/nc+1,nr+1/)
      CALL read_distribute_mod(zeta,filepars,varid,lbounds(1:2),ubounds(1:2),&
                             & nhdims,'C  ',fcomm,varatts=varatts)

   ENDIF

!
!2.3 3-D currents
!----------------
!

   IF (iopt_mode_3D.GT.0) THEN

!     ---X-current
      varid = varid + 1; nhdims = nhalo
      lbounds = (/1-nhalo,1-nhalo,1/); ubounds = (/nc+nhalo,nr+nhalo,nz/)
      CALL read_distribute_mod(uvel,filepars,varid,lbounds,ubounds,nhdims,&
                             & 'U  ',fcomm,varatts=varatts)

!     ---Y-current
      varid = varid + 1
      CALL read_distribute_mod(vvel,filepars,varid,lbounds,ubounds,nhdims,&
                             & 'V  ',fcomm,varatts=varatts)

   ENDIF

!  ---vertical current
   IF (iopt_mode_3D.GT.0.AND.iopt_grid_nodim.EQ.3) THEN
      varid = varid + 1; nhdims = (/1,0,1,0/)
      lbounds = (/0,0,1/); ubounds = (/nc,nr,nz+1/)
      CALL read_distribute_mod(wvel,filepars,varid,lbounds,ubounds,nhdims,&
                             & 'W  ',fcomm,varatts=varatts)
   ENDIF

!
!2.4 Density arrays
!------------------
!

   nhdims = nhalo
   lbounds = (/1-nhalo,1-nhalo,1/); ubounds = (/nc+nhalo,nr+nhalo,nz/)
   
!  ---temperature
   IF (iopt_temp.GT.0) THEN
      varid = varid + 1
      CALL read_distribute_mod(temp,filepars,varid,lbounds,ubounds,nhdims,&
                             & 'C  ',fcomm,varatts=varatts)
   ENDIF

!  ---salinity
   IF (iopt_sal.GT.0) THEN
      varid = varid + 1
      CALL read_distribute_mod(sal,filepars,varid,lbounds,ubounds,nhdims,&
                             & 'C  ',fcomm,varatts=varatts)
   ENDIF

!
!2.5 Turbulence arrays
!---------------------
!

   IF (iopt_vdif_coef.EQ.3) THEN

!     ---turbulence energy
      varid = varid + 1; nhdims = nhalo; ubounds = (/nc+nhalo,nr+nhalo,nz+1/)
      CALL read_distribute_mod(tke,filepars,varid,lbounds,ubounds,nhdims,'W  ',&
                             & fcomm,varatts=varatts)

!     ---mixing length
      IF (iopt_turb_ntrans.EQ.2.AND.iopt_turb_param.EQ.1) THEN
         varid = varid + 1
         CALL read_distribute_mod(zlmix,filepars,varid,lbounds,ubounds,&
                                & nhdims,'W  ',fcomm,varatts=varatts)

!     ---dissipation rate
      ELSEIF (iopt_turb_ntrans.EQ.2.AND.iopt_turb_param.EQ.2) THEN
         varid = varid + 1
         CALL read_distribute_mod(dissip,filepars,varid,lbounds,ubounds,&
                                & nhdims,'W  ',fcomm,varatts=varatts)
      ENDIF

   ENDIF

!
!2.6 Bottom stress arrays
!------------------------
!

   IF (iopt_bstres_form.EQ.2) THEN
      lbounds(1:2) = 0; ubounds(1:2) = (/nc,nr/)
      nhdims = (/1,0,1,0/)

      IF (iopt_bstres_drag.EQ.2) THEN
!        ---bottom drag coefficient
         varid = varid + 1
         CALL read_distribute_mod(bdragcoefatc,filepars,varid,lbounds(1:2),&
                                & ubounds(1:2),nhdims,'C  ',fcomm,&
                                & varatts=varatts)

      ELSEIF (iopt_bstres_drag.EQ.4) THEN
!        ---roughness length
         varid = varid + 1
         CALL read_distribute_mod(zroughatc,filepars,varid,lbounds(1:2),&
                                & ubounds(1:2),nhdims,'C  ',fcomm,&
                                & varatts=varatts)
      ENDIF

   ENDIF

!
!2.7 Tidal arrays
!----------------
!

   IF (nconobc.GT.0) THEN 
      varid = varid + 1
      CALL read_vars(fnode_obc,filepars,varid,varatts=varatts)
      varid = varid + 1
      CALL read_vars(phase_obc,filepars,varid,varatts=varatts)
   ENDIF
   IF (iopt_astro_tide.EQ.1) THEN
      varid = varid + 1
      CALL read_vars(fnode_astro,filepars,varid,varatts=varatts)
      varid = varid + 1
      CALL read_vars(phase_astro,filepars,varid,varatts=varatts)
   ENDIF

!
!2.8 Energy loss terms at weirs
!------------------------------
!

   IF (iopt_weibar.EQ.1) THEN
      IF (numwbaru.GT.0) THEN
         varid = varid + 1
         CALL read_vars(wbarelossu,filepars,varid,varatts=varatts)
      ENDIF
      IF (numwbarv.GT.0) THEN
         varid = varid + 1
         CALL read_vars(wbarelossv,filepars,varid,varatts=varatts)
      ENDIF
   ENDIF

!
!2.9 Open boundary arrays
!------------------------
!

   IF (iopt_obc_sal.EQ.1) THEN
      IF (nobu.GT.0) THEN
         varid = varid + 1
         CALL read_vars(obcsalatu,filepars,varid,varatts=varatts)
      ENDIF
      IF (nobv.GT.0) THEN
         varid = varid + 1
         CALL read_vars(obcsalatv,filepars,varid,varatts=varatts)
      ENDIF
   ENDIF
   IF (iopt_obc_temp.EQ.1) THEN
      IF (nobu.GT.0) THEN
         varid = varid + 1
         CALL read_vars(obctmpatu,filepars,varid,varatts=varatts)
      ENDIF
      IF (nobv.GT.0) THEN
         varid = varid + 1
         CALL read_vars(obctmpatv,filepars,varid,varatts=varatts)
      ENDIF
   ENDIF
   IF (iopt_obc_2D.EQ.1) THEN
      IF (nobu.GT.0) THEN
         varid = varid + 1
         CALL read_vars(obc2uvatu,filepars,varid,varatts=varatts)
      ENDIF
      IF (nobv.GT.0) THEN
         varid = varid + 1
         CALL read_vars(obc2uvatv,filepars,varid,varatts=varatts)
      ENDIF
   ENDIF
   IF (iopt_obc_3D.EQ.1) THEN
      IF (nobu.GT.0) THEN
         varid = varid + 1
         CALL read_vars(obc3uvatu,filepars,varid,varatts=varatts)
      ENDIF
      IF (nobv.GT.0) THEN
         varid = varid + 1
         CALL read_vars(obc3uvatv,filepars,varid,varatts=varatts)
      ENDIF
   ENDIF

ENDDO

!
!3. Close
!--------
!

CALL close_filepars(filepars)
DEALLOCATE (varatts)

modfiles(io_inicon,ics_phys,1) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE read_phsics

!========================================================================

SUBROUTINE simulation_start
!************************************************************************
!
! *simulation_start* Start next simulation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Initialisation.F90  V2.11
!
! Description - read title for next simulation
!             - return exit status
!             - open log and error files for intialisation
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - read_cif_params, usrdef_init_params
!
! Module calls - check_cif_lbounds_novars, check_init_params,
!                clock_date_time, conv_from_chars, default_init_params,
!                define_model_comms, error_abort, error_limits_var,
!                error_value_var, filepars_init, monitor_files, open_file,
!                open_filepars, read_cif_line
!
!************************************************************************
!

USE iopars
USE paralpars
USE switches
USE syspars
USE timepars
USE check_model, ONLY: check_init_params
USE cif_routines, ONLY: check_cif_lbound_novars, conv_from_chars, read_cif_line
USE datatypes_init, ONLY: filepars_init
USE default_model, ONLY: default_init_params
USE error_routines, ONLY: error_abort, error_limits_var,error_value_var
USE inout_routines, ONLY: monitor_files, open_file, open_filepars
USE time_routines, ONLY: clock_date_time

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: first_call = .TRUE.
CHARACTER (LEN=12) :: csimul
CHARACTER (LEN=lencifline) :: cline
CHARACTER (LEN=lencifvar), DIMENSION(MaxCIFVars) :: cvals
INTEGER :: i, iend, ierr,  lvar, numvars
INTEGER, SAVE :: iunit, numsimul


pglev = pglev + 1
procname(pglev) = 'simulation_start'

!
!1. Default parameters
!---------------------
!
!---log files
exitlog = .TRUE.; loglev1 = 0; loglev2 = 0

!---error files
errchk = .TRUE.; ioerr = 0; maxerrors = MaxErrMesgs; nerrs = 0

!---warning file
warning = .TRUE.; warnflag = .FALSE.

!---cif files
CALL filepars_init(ciffile)

!---timer
CALL clock_date_time(ClockTime)

!
!2. Number of simulations on first call
!--------------------------------------
!

IF (first_call) THEN
   ciflinenum = 0
   isimul = 0
   CALL open_file(iunit,'defruns','IN','A')
   ierr = 0; nosimul = 0; numsimul = 0
   DO WHILE (ierr.EQ.0)
      READ (iunit,*,IOSTAT=ierr) cline
      IF (ierr.EQ.0) THEN
         numsimul = numsimul + 1
         IF (cline(1:1).NE.cifcom) nosimul = nosimul + 1
      ENDIF
   END DO
   REWIND (iunit)
ENDIF

!
!3. Check for a  next simulation
!-------------------------------
!

next_simul = .FALSE.
i_310: DO i=isimul+1,numsimul
   READ (iunit,'(A)') cline
   ciflinenum = ciflinenum + 1
   iend = INDEX(cline,cifcom)
   IF (iend.EQ.1) THEN
      cline = ''
   ELSEIF (iend.GT.1) THEN
      cline = cline(1:iend-1)
   ENDIF
   IF (LEN_TRIM(cline).NE.0) THEN
      next_simul = .TRUE.
      isimul = i
      EXIT i_310
   ENDIF
ENDDO i_310

IF (.NOT.next_simul) RETURN

!
!4. Read parameters from 'defruns'
!---------------------------------
!
!---read input line
CALL read_cif_line(cline,cvals,numvars,filename='defruns')

!---run title
runtitle = ''
lvar = 1
CALL conv_from_chars(cvals(lvar),runtitle,lvar)
IF (LEN_TRIM(runtitle).EQ.0) THEN
   next_simul = .FALSE.
   RETURN
ENDIF
CALL check_cif_lbound_novars(numvars,3)

!---status of CIF file
ciffile%status = '0'
lvar = 2
CALL conv_from_chars(cvals(lvar),ciffile%status,lvar)

!---form of CIF file
ciffile%form = 'A'

!---name of CIF file
lvar = 3
CALL conv_from_chars(cvals(lvar),ciffile%filename,lvar)
IF (LEN_TRIM(ciffile%filename).EQ.0) THEN
   ciffile%filename = TRIM(runtitle)//'.cifmod'
ENDIF

!
!5. Open CIF file
!----------------
!

IF (ciffile%status.EQ.'R') THEN
   CALL open_filepars(ciffile)
   ciflinenum = 0
ENDIF

!
!6. Monitoring files
!-------------------
!
!---initialise by default
CALL default_init_params

!---user-defined
IF (ciffile%status.NE.'R') THEN
   CALL usrdef_init_params
ELSE
   CALL read_cif_params(icif_mon)
ENDIF

!---number of processes for the particle model
nprocspart = MERGE(1,0,iopt_part_model.EQ.2)

!---number of processes for flow model
IF (npworld.GT.1) nprocscoh = npworld - nprocspart - nprocswav

!---number of coupled models (MCT)
nocoupledmodels = MERGE(0,2,iopt_waves_model.EQ.0)

!---open files
CALL monitor_files

!---check initial parameters
CALL check_init_params
CALL error_abort('simulation_start',ierrno_comms)

!
!7. Start current-wave coupling
!------------------------------
!

IF (first_call.AND.iopt_MPI.EQ.1) CALL define_model_comms

!
!8. Write info
!-------------
!

IF (loglev1.GT.0) THEN
   WRITE (csimul,'(I12)') isimul; csimul = ADJUSTL(csimul)
   WRITE (iolog,'(A)') 'Simulation no.: '//TRIM(csimul)
   WRITE (iolog,'(A)') 'Execute: '//TRIM(runtitle)
ENDIF

first_call = .FALSE.
  
pglev = pglev - 1


RETURN

END SUBROUTINE simulation_start
