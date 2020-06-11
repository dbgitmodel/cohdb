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

MODULE reset_model
!************************************************************************
!
! *reset_model* Reset model parameters and arrays if needed
!
! Author - Patrick Luyten
!
! Revision - @(COHERENS)reset_model.f90  V2.11.1
!
! $Date: 2018-06-05 16:03:05 +0200 (Tue, 05 Jun 2018) $
!
! $Revision: 1142 $
!
! Description - 
!
! Reference -
!
! Generic routines - reset_out_files
!
! Routines - reset_initial_conditions, reset_mod_params, reset_mod_vars,
!            reset_out_gpars, reset_out_stats, reset_out_vars
!
!************************************************************************
!
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

INTERFACE reset_out_files
   MODULE PROCEDURE reset_out_files_1d, reset_out_files_2d
END INTERFACE


CONTAINS

!========================================================================

SUBROUTINE reset_initial_conditions
!************************************************************************
!
! *reset_initial_conditions* Reset initial conditions
!
! Author - Patrick Luyten
!
! Revision - @(COHERENS)reset_model.f90  V2.11.1
!
! Description - 
!
! Module calls -
!
!************************************************************************
!
USE density
USE fluxes
USE grid
USE gridpars
USE physpars
USE switches
USE wavevars

!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'reset_initial_conditions'
CALL log_timer_in()

!
!1. Density arrays
!-----------------
!
!---temperature
IF (iopt_temp.EQ.0) THEN
   k_110: DO k=1,nz
      WHERE (nodeatc.GT.0)
         temp(:,:,k) = temp_ref
      END WHERE
   ENDDO k_110
ENDIF

!---salinity
IF (iopt_sal.EQ.0) THEN
   k_120: DO k=1,nz
      WHERE (nodeatc.GT.0)
         sal(:,:,k) = sal_ref
      END WHERE
   ENDDO k_120
ENDIF

!
!2. Bottom drag coefficient/roughness length
!-------------------------------------------
!

IF (iopt_bstres_form.EQ.2) THEN
   IF (iopt_bstres_drag.EQ.1) THEN
      WHERE (nodeatc(0:ncloc,0:nrloc).GT.0)
         bdragcoefatc = bdragcoef_cst
      END WHERE
   ELSEIF  (iopt_bstres_drag.EQ.3) THEN
      WHERE (nodeatc(0:ncloc,0:nrloc).GT.0)
         zroughatc = zrough_cst
      END WHERE
   ENDIF
ENDIF

!
!3. Thicknes of the wave boundary layer
!--------------------------------------
!

IF (iopt_waves.EQ.1.AND.iopt_bstres_waves_bfric.EQ.0) THEN
   WHERE (maskatc_int)
      wavethickatc = wavethick_cst
   END WHERE
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE reset_initial_conditions

!========================================================================

SUBROUTINE reset_mod_params
!************************************************************************
!
! *reset_mod_params* Reset model parameters
!
! Author - Patrick Luyten
!
! Revision - @(COHERENS)reset_model.f90  V2.11
!
! Description - 
!
! Module calls - check_date, error_lbound_var, error_lbound_var_date,
!                num_time_steps, open_file, set_modfiles_name,
!                warning_reset_arr_struc, warning_reset_var
!
!************************************************************************
!
USE gridpars
USE nestgrids
USE paralpars
USE physpars
USE structures
USE switches
USE syspars
USE tide
USE timepars
USE cf90_routines, ONLY: cf90_inq_libvers
USE error_routines, ONLY: error_lbound_var, warning_reset_arr_struc, &
                        & warning_reset_var
USE modvars_routines, ONLY: set_modfiles_name
USE time_routines, ONLY: check_date, convert_date, error_lbound_var_date, &
                       & num_time_steps

!
!*Local variables
!
CHARACTER (LEN=lentime) :: cdatetimex
CHARACTER (LEN=12) :: cvar
CHARACTER (LEN=15) :: cval
CHARACTER (LEN=22) :: cstep
INTEGER :: idesc, ifil, iglb, igrd, iotype, msdelt, n
INTEGER, DIMENSION(8) :: intdate


procname(pglev+1) = 'reset_mod_params'
CALL log_timer_in()

!
!1. Reset switches if necessary
!------------------------------
!
!---hydrodynamics
IF (iopt_curr.EQ.0.OR.iopt_grid_nodim.EQ.1) THEN
   iopt_mode_2D = 0
ELSEIF (iopt_curr.EQ.1.OR.(iopt_grid_nodim.EQ.3.AND.iopt_hydro_impl.EQ.1)) THEN
   iopt_mode_2D = 1
ELSE
   iopt_mode_2D = 2
ENDIF
IF (iopt_curr.EQ.0.OR.iopt_grid_nodim.EQ.2) THEN
   iopt_mode_3D = 0
ELSEIF (iopt_curr.EQ.1) THEN
   iopt_mode_3D = 1
ELSE
   iopt_mode_3D = 2
ENDIF

!---no current case
IF (iopt_curr.EQ.0) THEN
   CALL warning_reset_var(iopt_hydro_impl,'iopt_hydro_impl',0)
   CALL warning_reset_var(iopt_adv_2D,'iopt_adv_2D',0)
   CALL warning_reset_var(iopt_adv_3D,'iopt_adv_3D',0)
   CALL warning_reset_var(iopt_adv_scal,'iopt_adv_scal',0)
   CALL warning_reset_var(iopt_hdif_2D,'iopt_hdif_2D',0)
   CALL warning_reset_var(iopt_hdif_3D,'iopt_hdif_3D',0)
   CALL warning_reset_var(iopt_bstres_form,'iopt_bstres_form',0)
   CALL warning_reset_var(iopt_astro_tide,'iopt_astro_tide',0)
   IF (iopt_hdif_scal.LT.2) THEN
      CALL warning_reset_var(iopt_vdif_coef,'iopt_vdif_coef',0)
   ENDIF
ENDIF

!---1-D application
IF (iopt_grid_nodim.EQ.1) THEN
   CALL warning_reset_var(iopt_grid_htype,'iopt_grid_htype',1)
   CALL warning_reset_var(iopt_grid_sph,'iopt_grid_sph',0)
   CALL warning_reset_var(iopt_dens_grad,'iopt_dens_grad',0)
   CALL warning_reset_var(iopt_hydro_impl,'iopt_hydro_impl',0)
   CALL warning_reset_var(iopt_astro_tide,'iopt_astro_tide',0)
   CALL warning_reset_var(iopt_adv_scal,'iopt_adv_scal',0)
   CALL warning_reset_var(iopt_adv_2D,'iopt_adv_2D',0)
   CALL warning_reset_var(iopt_adv_3D,'iopt_adv_3D',0)
   CALL warning_reset_var(iopt_adv_turb,'iopt_adv_turb',0)
   CALL warning_reset_var(iopt_hdif_coef,'iopt_hdif_coef',0)
   CALL warning_reset_var(iopt_nests,'iopt_nests',0)
   CALL warning_reset_var(iopt_obc_bio,'iopt_obc_bio',0)
   CALL warning_reset_var(iopt_obc_sal,'iopt_obc_sal',0)
   CALL warning_reset_var(iopt_obc_sed,'iopt_obc_sed',0)
   CALL warning_reset_var(iopt_obc_temp,'iopt_obc_temp',0)
   CALL warning_reset_var(iopt_obc_th,'iopt_obc_temp',0)
   CALL warning_reset_var(iopt_obc_2D,'iopt_obc_2D',0)
   CALL warning_reset_var(iopt_obc_2D_tang,'iopt_obc_2D_tang',0)
   CALL warning_reset_var(iopt_obc_3D,'iopt_obc_3D',0)
   CALL warning_reset_var(iopt_obc_3D_tang,'iopt_obc_3D_tang',0)
   CALL warning_reset_var(iopt_drycel,'iopt_drycel',0)
   CALL warning_reset_var(iopt_thndam,'iopt_thndam',0)
   CALL warning_reset_var(iopt_weibar,'iopt_weibar',0)
   IF (iopt_curr_wfall.EQ.2) THEN
      CALL warning_reset_var(iopt_curr_wfall,'iopt_curr_wfall',1)
   ENDIF
ELSE
   CALL warning_reset_var(iopt_sur_1D,'iopt_sur_1D',0)
ENDIF

!---2-D application
IF (iopt_grid_nodim.EQ.2) THEN
   CALL warning_reset_var(iopt_grid_vtype,'iopt_grid_vtype',1)
   CALL warning_reset_var(iopt_bstres_nodim,'iopt_bstres_nodim',2)
   CALL warning_reset_var(iopt_temp_optic,'iopt_temp_optic',0)
   CALL warning_reset_var(iopt_vdif_coef,'iopt_vdif_coef',0)
   CALL warning_reset_var(iopt_adv_3D,'iopt_adv_3D',0)
   CALL warning_reset_var(iopt_hdif_3D,'iopt_hdif_3D',0)
ENDIF

!---no 2-D mode
IF (iopt_mode_2D.LE.1) THEN
   CALL warning_reset_var(iopt_adv_2D,'iopt_avd_2D',0)
   CALL warning_reset_var(iopt_hdif_2D,'iopt_hdif_2D',0)
ENDIF

!---no 3-D mode
IF (iopt_mode_3D.LE.1) THEN
   CALL warning_reset_var(iopt_adv_3D,'iopt_avd_3D',0)
   CALL warning_reset_var(iopt_hdif_3D,'iopt_hdif_3D',0)
ENDIF

!---vertical grid
IF (iopt_grid_vtype.EQ.1) THEN
   CALL warning_reset_var(iopt_grid_vtype_transf,'iopt_grid_vtype_transf',0)
ENDIF

!---density
IF ((iopt_temp.EQ.0).AND.(iopt_sal.EQ.0)) THEN
   CALL warning_reset_var(iopt_dens,'iopt_dens',0)
ENDIF
IF (iopt_dens.EQ.0.AND.iopt_sed.EQ.0) THEN
   CALL warning_reset_var(iopt_dens_grad,'iopt_dens_grad',0)
ENDIF
IF (iopt_dens.LT.2.AND.iopt_dens_convect.EQ.1) THEN
   CALL warning_reset_var(iopt_dens_convect,'iopt_dens_convect',0)
ENDIF
IF (iopt_temp_sbc.GT.1) THEN
   CALL warning_reset_var(iopt_temp_optic,'iopt_temp_optic',0)
ENDIF

!---horizontal diffusion
IF (iopt_hdif_coef.EQ.0) THEN
   CALL warning_reset_var(iopt_hdif_2D,'iopt_hdif_2D',0)
   CALL warning_reset_var(iopt_hdif_3D,'iopt_hdif_3D',0)
   CALL warning_reset_var(iopt_hdif_scal,'iopt_hdif_scal',0)
   CALL warning_reset_var(iopt_hdif_turb,'iopt_hdif_turb',0)
ENDIF
IF ((iopt_hdif_2D.EQ.0).AND.(iopt_hdif_3D.EQ.0).AND.(iopt_hdif_scal.EQ.0).AND.&
  & (iopt_hdif_turb.EQ.0)) THEN
   CALL warning_reset_var(iopt_hdif_coef,'iopt_hdif_coef',0)
ENDIF
IF (iopt_hdif_scal.LT.2.OR.iopt_vdif_coef.EQ.0) THEN
   CALL warning_reset_var(iopt_vdif_rot,'iopt_vdif_rot',0)
ENDIF

!---bottom stress
IF (iopt_bstres_form.LT.2) THEN
   CALL warning_reset_var(iopt_bstres_drag,'iopt_bstres_drag',0)
ENDIF

!---weirs/barriers
IF (iopt_weibar.EQ.1) CALL warning_reset_var(iopt_arrint_3D,'iopt_arrint_3D',1)

!---waves
IF (iopt_waves.EQ.0) THEN
   CALL warning_reset_var(iopt_bstres_waves_bfric,'iopt_bstres_waves_bfric',0)
   CALL warning_reset_var(iopt_waves_couple,'iopt_waves_couple',0)
   CALL warning_reset_var(iopt_waves_curr,'iopt_waves_curr',0)
ENDIF
IF (iopt_waves_curr.EQ.0) THEN
   CALL warning_reset_var(iopt_waves_dissip,'iopt_waves_dissip',0)
   CALL warning_reset_var(iopt_waves_pres,'iopt_waves_pres',0)
ENDIF

!---particle model
IF (iopt_part_model.GT.2) THEN
   CALL warning_reset_var(iopt_part_write,'iopt_part_write',0)
ENDIF

!---morphology
IF (iopt_morph.EQ.0) THEN
   CALL warning_reset_var(iopt_tidal_accel,'iopt_tidal_accel',0)
ENDIF

!---number of open boundaries
nobu = nrvbu + nosbu; nobv = nrvbv + nosbv

!---no open boundary conditions if no open boundaries
IF (nobu.EQ.0.AND.nobv.EQ.0) THEN
   CALL warning_reset_var(iopt_obc_2D,'iopt_obc_2D',0)
   CALL warning_reset_var(iopt_obc_2D_tang,'iopt_obc_2D_tang',0)
   CALL warning_reset_var(iopt_obc_3D,'iopt_obc_3D',0)
   CALL warning_reset_var(iopt_obc_3D_tang,'iopt_obc_3D_tang',0)
ENDIF
IF ((nobu.EQ.0.AND.nobv.EQ.0).OR.iopt_adv_scal.EQ.0) THEN
   CALL warning_reset_var(iopt_obc_sal,'iopt_obc_sal',0)
   CALL warning_reset_var(iopt_obc_temp,'iopt_obc_temp',0)
   CALL warning_reset_var(iopt_obc_th,'iopt_obc_th',0)
   IF (iopt_sed.EQ.0) CALL warning_reset_var(iopt_obc_sed,'iopt_obc_sed',0)
   IF (iopt_biolgy.EQ.0) CALL warning_reset_var(iopt_obc_bio,'iopt_obc_bio',0)
ENDIF

!---open boundary switches
IF (iopt_sal.EQ.0) CALL warning_reset_var(iopt_obc_sal,'iopt_obc_sal',0)
IF (iopt_temp.EQ.0) CALL warning_reset_var(iopt_obc_temp,'iopt_obc_temp',0)
IF (iopt_sed.EQ.0) CALL warning_reset_var(iopt_obc_sed,'iopt_obc_sed',0)
IF (iopt_biolgy.EQ.0) CALL warning_reset_var(iopt_obc_bio,'iopt_obc_bio',0)

!---drying/wetting
IF (iopt_fld.EQ.0) CALL warning_reset_var(iopt_fld_alpha,'iopt_fld_alpha',0)

!---no tidal harmonics if no tides
IF (nconobc.EQ.0.AND.nconastro.EQ.0) THEN
   CALL warning_reset_var(iopt_astro_pars,'iopt_astro_pars',0)
ENDIF

!---astronomical tide
IF (iopt_astro_tide.EQ.1)  THEN
   CALL warning_reset_var(iopt_astro_pars,'iopt_astro_pars',2)
ENDIF

!---meteo
IF (iopt_meteo.EQ.1.AND.&
  &(surfacegrids(igrd_meteo,1)%nhtype.EQ.0.OR.iopt_grid_nodim.EQ.1)) THEN
   CALL warning_reset_var(iopt_meteo_precip,'iopt_meteo_precip',0)
   CALL warning_reset_var(iopt_meteo_pres,'iopt_meteo_pres',0)
ENDIF
IF (iopt_meteo.EQ.0) THEN
   CALL warning_reset_var(iopt_meteo_data,'iopt_meteo_data',0)
   CALL warning_reset_var(iopt_meteo_heat,'iopt_meteo_heat',0)
   CALL warning_reset_var(iopt_meteo_precip,'iopt_meteo_precip',0)
   CALL warning_reset_var(iopt_meteo_pres,'iopt_meteo_pres',0)
   CALL warning_reset_var(iopt_meteo_stres,'iopt_meteo_stres',0)
ENDIF
IF (iopt_temp.LT.2) THEN
   CALL warning_reset_var(iopt_meteo_heat,'iopt_meteo_heat',0)
ENDIF

!---surface fluxes
IF (iopt_meteo_data.EQ.2) THEN
   CALL warning_reset_var(iopt_sflux_pars,'iopt_sflux_pars',0)
   CALL warning_reset_var(iopt_sflux_strat,'iopt_sflux_strat',0)
ENDIF
IF (iopt_meteo_precip.EQ.0) THEN
   CALL warning_reset_var(iopt_sflux_precip,'iopt_sflux_precip',0)
ENDIF

!---turbulence
IF (iopt_vdif_coef.EQ.3) THEN
   IF( iopt_turb_ntrans.EQ.0) THEN
      CALL warning_reset_var(iopt_turb_stab_lev,'iopt_turb_stab_lev',1)
      CALL warning_reset_var(iopt_turb_iwlim,'iopt_turb_iwlim',0)
   ENDIF
ELSEIF (iopt_turb_ntrans.NE.0) THEN
   CALL warning_reset_var(iopt_turb_ntrans,'iopt_turb_ntrans',0)
ENDIF

!---tides
IF (iopt_grid_sph.EQ.0) THEN
   CALL warning_reset_var(iopt_astro_tide,'iopt_astro_tide',0)
ENDIF
IF (iopt_astro_tide.EQ.0) THEN
   CALL warning_reset_var(iopt_astro_anal,'iopt_astro_anal',0)
ENDIF

!---sediment transport and morphology
IF (iopt_sed.EQ.0) THEN
   CALL warning_reset_var(iopt_morph,'iopt_morph',0)
   CALL warning_reset_var(iopt_dar,'iopt_dar',0)
ENDIF

!---morphological acceleration
IF (iopt_tidal_accel.EQ.1) THEN
   CALL warning_reset_var(iopt_grid_nodim,'iopt_grid_nodim',2)
   CALL warning_reset_var(iopt_dens,'iopt_dens',0)
   CALL warning_reset_var(iopt_sal,'iopt_sal',0)
   CALL warning_reset_var(iopt_temp,'iopt_temp',0)
   CALL warning_reset_var(iopt_meteo,'iopt_meteo',0)
   CALL warning_reset_var(iopt_waves,'iopt_waves',0)
   CALL warning_reset_var(iopt_grid_vtype,'iopt_grid_vtype',1)
   CALL warning_reset_var(iopt_bstres_nodim,'iopt_bstres_nodim',2)
   CALL warning_reset_var(iopt_vdif_coef,'iopt_vdif_coef',0)
   CALL warning_reset_var(iopt_adv_3D,'iopt_adv_3D',0)
   CALL warning_reset_var(iopt_hdif_3D,'iopt_hdif_3D',0)
   CALL warning_reset_var(ic3d,'ic3d',1)
ENDIF

!
!2. Model parameters
!-------------------
!
!---time step
IF (delt2d.LT.1000.0) THEN
   msdelt = NINT(1000.0*delt2d)
   delt2d = msdelt/1000.0
ELSE
   delt2d = INT(delt2d)
ENDIF
CALL error_lbound_var(delt2d,'delt2d',0.0,.FALSE.)

!---time parameters
IF (iopt_grid_nodim.NE.3.OR.iopt_hydro_impl.EQ.1.OR.iopt_part_model.GT.2) THEN
   CALL warning_reset_var(ic3d,'ic3d',1)
ENDIF
IF (iopt_dens_convect.EQ.1.AND.iccvt.EQ.0) THEn
   CALL warning_reset_var(iccvt,'iccvt',ic3d)
ENDIF

!---3-D time step
delt3d = MERGE(delt2d,ic3d*delt2d,iopt_grid_nodim.EQ.1.OR.&
             & iopt_hydro_impl.EQ.1.OR.iopt_part_model.GE.3)
WRITE (cval,'(G15.7)') delt3d; cval = ADJUSTL(cval)
WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'delt3d = '//TRIM(cval)

!---date/time
IStartDateTime = convert_date(CStartDateTime)
CALL check_date(IStartDateTime,'IStartDateTime')
IEndDateTime = convert_date(CEndDateTime)
CALL check_date(IEndDateTime,'IEndDateTime')

!---number of time steps
IF (iopt_part_model.LE.2) THEN
   CALL error_lbound_var_date(CEndDateTime,'CEndDateTime',CStartDateTime,&
                           & .FALSE.)
   CALL num_time_steps(CStartDateTime,CEndDateTime,delt2d,nstep)
   CALL error_lbound_var(nstep,'nstep',0,.FALSE.)
   hydro_end = nstep
   IF (loglev1.GT.0) THEN
      WRITE (cstep,'(I22)') nstep; cstep = ADJUSTL(cstep)
      WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'nstep = '//TRIM(cstep)
   ENDIF
ENDIF

!---processor clock
IF (loglev1.GE.0) THEN
   WRITE (cvar,'(I12)') npcc_max/npcc_rate; cvar = ADJUSTL(cvar)
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//&
                     & 'Maximum time interval registered by processor '//&
                     & 'clock: '//TRIM(cvar)
ENDIF

!---optical parameters
IF (iopt_temp_optic.EQ.1) THEN
   IF (opt_frac.LT.0.0) CALL warning_reset_var(opt_frac,'opt_frac',0.0)
   IF (opt_frac.GT.1.0) CALL warning_reset_var(opt_frac,'opt_frac',1.0)
ENDIF

!---factors for horizontal diffusion
SELECT CASE (iopt_hdif_coef)
   CASE (0)
      hdifmom_fac = 0.0; hdifscal_fac = 0.0
   CASE (1)
      hdifmom_fac = hdifmom_cst; hdifscal_fac = hdifscal_cst
   CASE(2)
      hdifmom_fac = smag_coef_mom; hdifscal_fac = smag_coef_scal
END SELECT
   
!---bathymetry flagging
depmean_flag_eps = 0.00001*ABS(depmean_flag)

!---minimum water depth
IF (iopt_fld.EQ.0.AND.iopt_weibar.EQ.0) THEN
   dcrit_fld = 0.0; dmin_fld = 0.0
ENDIF

!---implicity factors
IF (iopt_hydro_impl.EQ.0) THEN
   CALL warning_reset_var(nomglevels,'nomglevels',1)
   CALL warning_reset_var(theta_sur,'theta_sur',0.0)
ENDIF
SELECT CASE (iopt_cor_impl)
   CASE(0); CALL warning_reset_var(theta_cor,'theta_cor',0.0)
   CASE(1)
      IF (theta_cor.LT.0.0) CALL warning_reset_var(theta_cor,'theta_cor',0.0)
      IF (theta_cor.GT.1.0) CALL warning_reset_var(theta_cor,'theta_cor',1.0)
   CASE(2); CALL warning_reset_var(theta_cor,'theta_cor',1.0)
END SELECT
SELECT CASE (iopt_vadv_impl)
   CASE(0); CALL warning_reset_var(theta_vadv,'theta_vadv',0.0)
   CASE(1)
      IF (theta_vadv.LT.0.0) &
         & CALL warning_reset_var(theta_vadv,'theta_vadv',0.0)
      IF (theta_vadv.GT.1.0) &
         & CALL warning_reset_var(theta_vadv,'theta_vadv',1.0)
   CASE(2); CALL warning_reset_var(theta_vadv,'theta_vadv',1.0)
END SELECT
SELECT CASE (iopt_vdif_impl)
   CASE(0); CALL warning_reset_var(theta_vdif,'theta_vdif',0.0)
   CASE(1)
      IF (theta_vdif.LT.0.0) CALL warning_reset_var(theta_vdif,'theta_vdif',0.0)
      IF (theta_vdif.GT.1.0) CALL warning_reset_var(theta_vdif,'theta_vdif',1.0)
   CASE(2); CALL warning_reset_var(theta_vdif,'theta_vdif',1.0)
END SELECT

!---number of outer iterations in the implicit scheme
IF (iopt_hydro_impl.EQ.0) CALL warning_reset_var(itsimp,'itsimp',1)

!---number of processes
IF (.NOT.parallel_set) CALL warning_reset_var(nprocs,'nprocs',1)

!---restart times
n_210: DO n=1,norestarts
   IF (ntrestart(n).EQ.int_fill) ntrestart(n) = nstep
ENDDO n_210

!---structures
IF (iopt_drycel.EQ.0) CALL warning_reset_var(numdry,'numdry',0)
IF (iopt_thndam.EQ.0) THEN
   CALL warning_reset_var(numthinu,'numthinu',0)
   CALL warning_reset_var(numthinv,'numthinv',0)
ENDIF
IF (iopt_weibar.EQ.0) THEN
   CALL warning_reset_var(numwbaru,'numwbaru',0)
   CALL warning_reset_var(numwbarv,'numwbarv',0)
ENDIF
IF (iopt_dischr.EQ.0) CALL warning_reset_var(numdis,'numdis',0)

!---log file counter
IF (iopt_part_model.LE.2.AND.runlog_count.EQ.int_fill) runlog_count = nstep

!---timer report
IF (iopt_verif.EQ.1) CALL warning_reset_var(levtimer,'levtimer',3)

!---number of time series output sets
!IF (iopt_verif.EQ.1) THEN
!   CALL warning_reset_var(nosetstsr,'nosetstsr',1)
!ENDIF

!
!3. Reset file attributes
!------------------------
!
!3.1 Status
!----------
!

idesc_311: DO idesc=1,MaxIOTypes
ifil_311: DO ifil=1,MaxIOFiles
   SELECT CASE (idesc)
      CASE (io_mppmod,io_modgrd,io_metabs,io_sstabs,io_wavabs,io_bioabs,&
          & io_metrel,io_sstrel,io_wavrel,io_biorel,io_sedspc,io_biospc,&
          & io_rlxobc,io_nstspc,io_metsur,io_sstsur,io_biosur,io_wavsur,&
          & io_drycel,io_thndam,io_weibar,io_disspc,io_disloc,io_disvol,&
          & io_discur,io_dissal,io_dissed,io_distmp)
         IF (ifil.GT.1) modfiles(idesc,ifil,:)%status = '0'
      CASE (io_inicon)
         IF (ifil.GT.MaxInitFiles) modfiles(idesc,ifil,:)%status = '0'
   END SELECT

   SELECT CASE (idesc)
      CASE (io_fincon)
         IF (norestarts.EQ.0) THEN
            IF (ifil.EQ.ics_phys) THEN
               CALL warning_reset_arr_struc(modfiles(idesc,ifil,2)%status,&
                                 & 'modfiles','%status','0',3,(/idesc,ifil,2/))
            ELSEIF (ifil.EQ.ics_bio.AND.iopt_biolgy.EQ.0) THEN
               CALL warning_reset_arr_struc(modfiles(idesc,ifil,2)%status,&
                                 & 'modfiles','%status','0',3,(/idesc,ifil,2/))
            ELSEIF (ifil.EQ.ics_sed.AND.iopt_sed.EQ.0) THEN
               CALL warning_reset_arr_struc(modfiles(idesc,ifil,2)%status,&
                                 & 'modfiles','%status','0',3,(/idesc,ifil,2/))
            ENDIF
         ENDIF
      CASE (io_2uvnst,io_2xynst,io_3uvnst,io_3xynst,io_salnst,io_tmpnst,&
          & io_sednst,io_bionst)
         IF (iopt_nests.EQ.0) THEN
            CALL warning_reset_arr_struc(modfiles(idesc,ifil,2)%status,&
                                 & 'modfiles','%status','0',3,(/idesc,ifil,2/))
         ENDIF
      CASE (io_bioabs,io_biorel)
         IF (ifil.LE.MaxGridFiles) THEN
            iotype_3115: DO iotype=1,2
               IF (iopt_biolgy.EQ.0.OR.&
                  & surfacegrids(igrd_bio,ifil)%nhtype.LE.1.OR.&
                  & surfacegrids(igrd_bio,ifil)%nhtype.EQ.4) THEN
                  CALL warning_reset_arr_struc(&
                              & modfiles(idesc,ifil,iotype)%status,'modfiles',&
                              & '%status','0',3,(/idesc,ifil,iotype/))
               ENDIF
            ENDDO iotype_3115
         ENDIF
   END SELECT
ENDDO ifil_311
ENDDO idesc_311

iotype_312: DO iotype=1,2
   IF (.NOT.parallel_set) THEN
      CALL warning_reset_arr_struc(modfiles(io_mppmod,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_mppmod,1,iotype/))
   ENDIF
   IF (iopt_sed.EQ.0) THEN
      modfiles(io_inicon,ics_sed,:)%status = '0'
      modfiles(io_fincon,ics_sed,:)%status = '0'
      modfiles(io_dissed,1,:)%status = '0'
   ENDIF
   IF (iopt_biolgy.EQ.0) THEN
      modfiles(io_inicon,ics_bio,:)%status = '0'
      modfiles(io_fincon,ics_bio,:)%status = '0'
   ENDIF
   IF (iopt_meteo.EQ.0.OR.surfacegrids(igrd_meteo,1)%nhtype.LE.1.OR.&
     & surfacegrids(igrd_meteo,1)%nhtype.EQ.4) THEN
      CALL warning_reset_arr_struc(modfiles(io_metabs,1,iotype)%status,&
           & 'modfiles','%status','0',3,(/io_metabs,1,iotype/))
      CALL warning_reset_arr_struc(modfiles(io_metrel,1,iotype)%status,&
           & 'modfiles','%status','0',3,(/io_metrel,1,iotype/))
   ENDIF
   IF (iopt_meteo.EQ.0) THEN
      CALL warning_reset_arr_struc(modfiles(io_metsur,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_metsur,1,iotype/))
   ENDIF
   IF (iopt_temp_sbc.EQ.1.OR.surfacegrids(igrd_sst,1)%nhtype.LE.1.OR.&
      & surfacegrids(igrd_sst,1)%nhtype.EQ.4) THEN
      CALL warning_reset_arr_struc(modfiles(io_sstabs,1,iotype)%status,&
           & 'modfiles','%status','0',3,(/io_sstabs,1,iotype/))
      CALL warning_reset_arr_struc(modfiles(io_sstrel,1,iotype)%status,&
           & 'modfiles','%status','0',3,(/io_sstrel,1,iotype/))
   ENDIF
   IF (iopt_temp_sbc.EQ.1) THEN
      CALL warning_reset_arr_struc(modfiles(io_sstsur,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_sstsur,1,iotype/))
   ENDIF
   IF (iopt_waves.EQ.0.OR.surfacegrids(igrd_waves,1)%nhtype.LE.1.OR.&
        & surfacegrids(igrd_waves,1)%nhtype.EQ.4) THEN
      CALL warning_reset_arr_struc(modfiles(io_wavabs,1,iotype)%status,&
           & 'modfiles','%status','0',3,(/io_wavabs,1,iotype/))
      CALL warning_reset_arr_struc(modfiles(io_wavrel,1,iotype)%status,&
           & 'modfiles','%status','0',3,(/io_wavrel,1,iotype/))
   ENDIF
   IF (iopt_waves.EQ.0) THEN
      CALL warning_reset_arr_struc(modfiles(io_wavsur,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_wavsur,1,iotype/))
   ENDIF
   IF (iopt_drycel.EQ.0) THEN
      CALL warning_reset_arr_struc(modfiles(io_drycel,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_drycel,1,iotype/))
   ENDIF
   IF (iopt_thndam.EQ.0) THEN
      CALL warning_reset_arr_struc(modfiles(io_thndam,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_thndam,1,iotype/))
   ENDIF
   IF (iopt_weibar.EQ.0) THEN
      CALL warning_reset_arr_struc(modfiles(io_weibar,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_weibar,1,iotype/))
   ENDIF
   IF (iopt_biolgy.EQ.0) THEN
      CALL warning_reset_arr_struc(modfiles(io_biosur,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_biosur,1,iotype/))
   ENDIF
   IF (iopt_nests.EQ.0) THEN
      CALL warning_reset_arr_struc(modfiles(io_nstspc,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_nstspc,1,iotype/))
   ENDIF
   IF (iopt_obc_relax.EQ.0) THEN
      CALL warning_reset_arr_struc(modfiles(io_rlxobc,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_rlxobc,1,iotype/))
   ENDIF

   IF (iopt_dischr.EQ.0) THEN
      CALL warning_reset_arr_struc(modfiles(io_disspc,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_disspc,1,iotype/))
      CALL warning_reset_arr_struc(modfiles(io_disloc,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_disloc,1,iotype/))
      CALL warning_reset_arr_struc(modfiles(io_disvol,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_disvol,1,iotype/))
      CALL warning_reset_arr_struc(modfiles(io_discur,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_discur,1,iotype/))
      CALL warning_reset_arr_struc(modfiles(io_dissal,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_dissal,1,iotype/))
      CALL warning_reset_arr_struc(modfiles(io_dissed,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_dissed,1,iotype/))
      CALL warning_reset_arr_struc(modfiles(io_distmp,1,iotype)%status,&
                            & 'modfiles','%status','0',3,(/io_distmp,1,iotype/))
   ENDIF

   ifil_3121: DO ifil=1,MaxIOFiles
      IF (iopt_nests.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_nstabs,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_nstabs,ifil,iotype/))
         CALL warning_reset_arr_struc(modfiles(io_nstrel,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_nstrel,ifil,iotype/))
         CALL warning_reset_arr_struc(modfiles(io_nstspc,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_nstspc,ifil,iotype/))
         CALL warning_reset_arr_struc(modfiles(io_2uvnst,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_2uvnst,ifil,iotype/))
         CALL warning_reset_arr_struc(modfiles(io_2xynst,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_2xynst,ifil,iotype/))
         CALL warning_reset_arr_struc(modfiles(io_3uvnst,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_3uvnst,ifil,iotype/))
         CALL warning_reset_arr_struc(modfiles(io_3xynst,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_3xynst,ifil,iotype/))
         CALL warning_reset_arr_struc(modfiles(io_salnst,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_salnst,ifil,iotype/))
         CALL warning_reset_arr_struc(modfiles(io_tmpnst,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_tmpnst,ifil,iotype/))
         CALL warning_reset_arr_struc(modfiles(io_sednst,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_sednst,ifil,iotype/))
         CALL warning_reset_arr_struc(modfiles(io_bionst,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_bionst,ifil,iotype/))
      ENDIF
      IF (iopt_sur_1D.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_1uvsur,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_1uvsur,ifil,iotype/))
      ENDIF
      IF (iopt_obc_2D.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_2uvobc,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_2uvobc,ifil,iotype/))
      ENDIF
      IF (iopt_obc_2D_tang.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_2xyobc,ifil,iotype)%status,&
                                   & 'modfiles','%status','0',3,&
                                   & (/io_2xyobc,ifil,iotype/))
      ENDIF
      IF (iopt_obc_3D.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_3uvobc,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_3uvobc,ifil,iotype/))
      ENDIF
      IF (iopt_obc_3D_tang.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_3xyobc,ifil,iotype)%status,&
                                   & 'modfiles','%status','0',3,&
                                   & (/io_3xyobc,ifil,iotype/))
      ENDIF
      IF (iopt_obc_sal.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_salobc,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_salobc,ifil,iotype/))
      ENDIF
      IF (iopt_sal.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_salnst,ifil,iotype)%status,&
                                   & 'modfiles','%status','0',3,&
                                   & (/io_salnst,ifil,iotype/))
      ENDIF
      IF (iopt_obc_temp.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_tmpobc,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_tmpobc,ifil,iotype/))
      ENDIF
      IF (iopt_temp.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_tmpnst,ifil,iotype)%status,&
                                   & 'modfiles','%status','0',3,&
                                   & (/io_tmpnst,ifil,iotype/))
      ENDIF
      IF (iopt_obc_sed.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_sedobc,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_sedobc,ifil,iotype/))
      ENDIF
      IF (iopt_sed.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_sednst,ifil,iotype)%status,&
                                   & 'modfiles','%status','0',3,&
                                   & (/io_sednst,ifil,iotype/))
      ENDIF
      IF (iopt_obc_bio.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_bioobc,ifil,iotype)%status,&
                                   & 'modfiles','%status','0',3,&
                                   & (/io_bioobc,ifil,iotype/))
      ENDIF
      IF (iopt_biolgy.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_bionst,ifil,iotype)%status,&
                                   & 'modfiles','%status','0',3,&
                                   & (/io_bionst,ifil,iotype/))
      ENDIF
      IF (iopt_dischr.EQ.0) THEN
         CALL warning_reset_arr_struc(modfiles(io_disloc,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_disloc,ifil,iotype/))
         CALL warning_reset_arr_struc(modfiles(io_disvol,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_disvol,ifil,iotype/))
         CALL warning_reset_arr_struc(modfiles(io_discur,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_discur,ifil,iotype/))
         CALL warning_reset_arr_struc(modfiles(io_dissal,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_dissal,ifil,iotype/))
         CALL warning_reset_arr_struc(modfiles(io_distmp,ifil,iotype)%status,&
                         & 'modfiles','%status','0',3,(/io_distmp,ifil,iotype/))
      ENDIF
   ENDDO ifil_3121
ENDDO iotype_312

iotype_313: DO iotype=1,2

   IF (surfacegrids(igrd_meteo,1)%surcoords.EQ.1) THEN
      CALL warning_reset_arr_struc(modfiles(io_metrel,1,iotype)%status,&
                                & 'modfiles','%status','0',3,&
                                & (/io_metrel,1,iotype/))
   ELSEIF (surfacegrids(igrd_meteo,1)%surcoords.EQ.2) THEN
      CALL warning_reset_arr_struc(modfiles(io_metabs,1,iotype)%status,&
                                & 'modfiles','%status','0',3,&
                                & (/io_metabs,1,iotype/))
   ENDIF
   IF (surfacegrids(igrd_sst,1)%surcoords.EQ.1) THEN
      CALL warning_reset_arr_struc(modfiles(io_sstrel,1,iotype)%status,&
                                & 'modfiles','%status','0',3,&
                                & (/io_sstrel,1,iotype/))
   ELSEIF (surfacegrids(igrd_sst,1)%surcoords.EQ.2) THEN
      CALL warning_reset_arr_struc(modfiles(io_sstabs,1,iotype)%status,&
                                & 'modfiles','%status','0',3,&
                                & (/io_sstabs,1,iotype/))
   ENDIF
   IF (surfacegrids(igrd_waves,1)%surcoords.EQ.1) THEN
      CALL warning_reset_arr_struc(modfiles(io_wavrel,1,iotype)%status,&
                                & 'modfiles','%status','0',3,&
                                & (/io_wavrel,1,iotype/))
   ELSEIF (surfacegrids(igrd_waves,1)%surcoords.EQ.2) THEN
      CALL warning_reset_arr_struc(modfiles(io_wavabs,1,iotype)%status,&
                                & 'modfiles','%status','0',3,&
                                & (/io_wavabs,1,iotype/))
   ENDIF
   IF (surfacegrids(igrd_bio,1)%surcoords.EQ.1) THEN
      CALL warning_reset_arr_struc(modfiles(io_biorel,1,iotype)%status,&
                                & 'modfiles','%status','0',3,&
                                & (/io_biorel,1,iotype/))
   ELSEIF (surfacegrids(igrd_bio,1)%surcoords.EQ.2) THEN
      CALL warning_reset_arr_struc(modfiles(io_bioabs,1,iotype)%status,&
                                & 'modfiles','%status','0',3,&
                                & (/io_bioabs,1,iotype/))
   ENDIF

ENDDO iotype_313

!
!3.2 Defined
!-----------
!

idesc_320: DO idesc=1,MaxIOTypes
ifil_320: DO ifil=1,MaxIOFiles
   modfiles(idesc,ifil,:)%defined = modfiles(idesc,ifil,:)%status.NE.'0'
ENDDO ifil_320
ENDDO idesc_320

!
!3.3 Filename
!------------
!

iotype_330: DO iotype=1,2
idesc_330: DO idesc=1,MaxIOTypes
ifil_330: DO ifil=1,MaxIOFiles
   IF (modfiles(idesc,ifil,iotype)%defined) THEN
      CALL set_modfiles_name(idesc,ifil,iotype)
   ENDIF
ENDDO ifil_330
ENDDO idesc_330
ENDDO iotype_330

!
!3.4 Largest index of active files
!---------------------------------
!

iotype_340: DO iotype=1,2
idesc_340: DO idesc=1,MaxIOTypes
   maxdatafiles(idesc,iotype) = 0
   ifil_341: DO ifil=1,MaxIOFiles
      IF (modfiles(idesc,ifil,iotype)%defined) THEN
         maxdatafiles(idesc,iotype) = ifil
      ENDIF
   ENDDO ifil_341
ENDDO idesc_340
ENDDO iotype_340

!
!3.5 File units
!--------------
!

idesc_350: DO idesc=1,MaxIOTypes
ifil_350: DO ifil=1,MaxIOFiles
   modfiles(idesc,ifil,:)%iunit = int_fill
ENDDO ifil_350
ENDDO idesc_350

!
!3.6 Float type
!--------------
!

iotype_360: DO iotype=1,2
idesc_360: DO idesc=1,MaxIOTypes
ifil_360: DO ifil=1,MaxIOFiles
   IF (TRIM(modfiles(idesc,ifil,iotype)%floattype).EQ.'') THEN
      modfiles(idesc,ifil,iotype)%floattype = &
                                       & MERGE('S','D',float_type.EQ.real_type)
   ENDIF
ENDDO ifil_360
ENDDO idesc_360
ENDDO iotype_360

!
!3.7 Disabling fill values
!-------------------------
!

idesc_370: DO idesc=1,MaxIOTypes
   SELECT CASE (idesc)
   CASE (io_mppmod,io_modgrd,io_metabs,io_sstabs,io_wavabs,io_bioabs,&
       & io_nstabs,io_metrel,io_sstrel,io_wavrel,io_biorel,io_nstrel,&
       & io_sedspc,io_biospc,io_darspc,io_1uvsur,io_2uvobc,io_2xyobc,&
       & io_3uvobc,io_3xyobc,io_salobc,io_tmpobc,io_sedobc,io_bioobc,&
       & io_rlxobc,io_nstspc,io_2uvnst,io_2xynst,io_3uvnst,io_3xynst,&
       & io_salnst,io_tmpnst,io_sednst,io_bionst,io_metsur,io_drycel,&
       & io_thndam,io_weibar,io_disspc,io_disloc,io_disvol,io_discur,&
       & io_dissal,io_dissed,io_distmp)
         modfiles(idesc,:,:)%fill = .FALSE.
   END SELECT
ENDDO idesc_370
      
!
!3.8 Time coordinate
!-------------------
!

idesc_380: DO idesc=1,MaxIOTypes
ifil_380: DO ifil=1,MaxIOFiles
   SELECT CASE (idesc)
   CASE (io_mppmod,io_modgrd,io_metabs,io_sstabs,io_wavabs,io_bioabs,&
        & io_nstabs,io_metrel,io_sstrel,io_wavrel,io_biorel,io_nstrel,&
        & io_sedspc,io_biospc,io_rlxobc,io_nstspc,io_drycel,io_thndam,&
        & io_weibar,io_disspc)
         modfiles(idesc,ifil,:)%nocoords = 0
      CASE (io_1uvsur,io_2uvobc,io_2xyobc,io_3uvobc,io_3xyobc,io_salobc,&
          & io_tmpobc,io_sedobc,io_bioobc)
         IF (ifil.EQ.1) THEN
            modfiles(idesc,ifil,:)%nocoords = 0
         ELSE
            modfiles(idesc,ifil,:)%nocoords = 1
         ENDIF
      CASE (io_inicon,io_fincon)
         modfiles(idesc,ifil,:)%nocoords = MERGE(0,1,ifil.GT.3)
      CASE (io_2uvnst,io_2xynst,io_3uvnst,io_3xynst,io_salnst,io_tmpnst,&
          & io_sednst,io_bionst,io_metsur,io_sstsur,io_wavsur,io_biosur,&
          & io_disloc,io_disvol,io_discur,io_dissal,io_dissed,io_distmp)
         modfiles(idesc,ifil,:)%nocoords = 1
   END SELECT
ENDDO ifil_380
ENDDO idesc_380

!
!3.9 Time limits
!---------------
!

idesc_390: DO idesc=1,MaxIOTypes

   SELECT CASE (idesc)
      CASE (io_mppmod,io_modgrd,io_metabs,io_sstabs,io_wavabs,io_bioabs,&
          & io_nstabs,io_metrel,io_sstrel,io_wavrel,io_biorel,io_nstrel,&
          & io_sedspc,io_biospc,io_rlxobc,io_nstspc,io_drycel,io_thndam,&
          & io_weibar,io_disspc)
         modfiles(idesc,1,1)%tlims = 0
      CASE (io_1uvsur,io_2uvobc,io_3uvobc,io_salobc,io_tmpobc,io_sedobc,&
          & io_bioobc)
         modfiles(idesc,1,1)%tlims = 0
   END SELECT

   ifil_391: DO ifil=1,MaxIOFiles
      SELECT CASE (idesc)
         CASE (io_2uvnst,io_2xynst,io_3uvnst,io_3xynst,io_salnst,io_tmpnst,&
             & io_sednst,io_bionst,io_parphs)
            n_3911: DO n=1,3
               modfiles(idesc,ifil,2)%tlims(n) = MERGE(nstep,&
                                  & modfiles(idesc,ifil,2)%tlims(n),&
                                  & modfiles(idesc,ifil,2)%tlims(n).EQ.int_fill)
            ENDDO n_3911
            modfiles(idesc,ifil,1)%tlims = modfiles(idesc,ifil,2)%tlims
         CASE DEFAULT
            n_3912: DO n=1,3
               modfiles(idesc,ifil,1)%tlims(n) = MERGE(nstep,&
                                  & modfiles(idesc,ifil,1)%tlims(n),&
                                  & modfiles(idesc,ifil,1)%tlims(n).EQ.int_fill)
            ENDDO n_3912
            modfiles(idesc,ifil,2)%tlims = modfiles(idesc,ifil,1)%tlims
      END SELECT
   ENDDO ifil_391

ENDDO idesc_390

!
!4. Reset data grid attributes
!------------------------------
!
!---1-D application
IF (iopt_grid_nodim.EQ.1) iopt_grid_htype = 1

!---model grid
surfacegrids(igrd_model,1)%nhtype = iopt_grid_htype

!---meteo
IF (iopt_meteo.EQ.1.AND.&
  &(surfacegrids(igrd_meteo,1)%nhtype.EQ.0.OR.iopt_grid_nodim.EQ.1)) THEN
   CALL warning_reset_var(iopt_meteo_pres,'iopt_meteo_pres',0)
ENDIF

!---SW corners
igrd_410: DO igrd=1,MaxGridTypes
   IF (igrd.NE.igrd_model.AND.surfacegrids(igrd,1)%nhtype.GE.3) THEN
      surfacegrids(igrd,1)%x0dat = 0.0
      surfacegrids(igrd,1)%y0dat = 0.0
   ENDIF
ENDDO igrd_410

!---resolution
WHERE (surfacegrids%nhtype.EQ.0)
   surfacegrids%n1dat = 1
   surfacegrids%n2dat = 1
END WHERE
WHERE (surfacegrids%nhtype.EQ.4)
   surfacegrids%n1dat = nc
   surfacegrids%n2dat = nr
END WHERE

!---rotated grids
WHERE (surfacegrids%nhtype.GE.3)
   surfacegrids%rotated = .FALSE.
END WHERE
WHERE (.NOT.surfacegrids%rotated)
   surfacegrids%x0rot = surfacegrids%x0dat
   surfacegrids%y0rot = surfacegrids%y0dat
   surfacegrids%gridangle = 0.0
   surfacegrids%longpole = 0.0
END WHERE

!---grid spacings
igrd_420: DO igrd=1,MaxGridTypes
   IF (igrd.NE.igrd_model.AND.surfacegrids(igrd,1)%nhtype.NE.1) THEN
      surfacegrids(igrd,1)%delxdat = 0.0
      surfacegrids(igrd,1)%delydat = 0.0
      surfacegrids(igrd,1)%x0dat = 0.0
      surfacegrids(igrd,1)%y0dat = 0.0
   ENDIF
ENDDO igrd_420

!---surface waves
IF (iopt_waves_couple.GT.0) THEN
   surfacegrids(igrd_waves,1)%extrapol = .TRUE.
   surfacegrids(igrd_waves,1)%land = .TRUE.
ENDIF

!
!5. 1-D/2-D application
!----------------------
!

IF (iopt_grid_nodim.EQ.1) THEN
   nc = 3; nr = 3
   nosbu = 4; nosbv = 4
   nrvbu = 0; nrvbv = 0 
   nobu = 4; nobv = 4
   surfacegrids(igrd_model,1)%delxdat = 1.0
   surfacegrids(igrd_model,1)%delydat = 1.0
   surfacegrids(igrd_model,1)%x0dat = 0.0
   surfacegrids(igrd_model,1)%y0dat = 0.0
   surfacegrids(igrd_model,1)%rotated = .FALSE.
ELSEIF (iopt_grid_nodim.EQ.2) THEN
   CALL warning_reset_var(nz,'nz',1)
ENDIF

!
!6. CF global attributes
!-----------------------
!

iglb_710: DO iglb=1,numglbatts
   SELECT CASE (iglb)
      CASE (1)
         glbatts(iglb)%name = 'Conventions'
         glbatts(iglb)%value = CF_version
      CASE (2)
         glbatts(iglb)%name = 'title'
      CASE (3)
         glbatts(iglb)%name = 'history'
         CALL DATE_AND_TIME(VALUES=intdate)
         intdate(4)= intdate(5); intdate(5) = intdate(6)
         intdate(6)= intdate(7); intdate(7) = intdate(8)
         cdatetimex = convert_date(intdate(1:7),seps='- ::')
         IF (TRIM(history_CF).EQ.'') THEN
            history_CF = 'created '//cdatetimex(1:19)
         ENDIF
         glbatts(iglb)%value = history_CF
      CASE (4)
         glbatts(iglb)%name = 'institution'
         glbatts(iglb)%value = institution_CF
      CASE (5)
         glbatts(iglb)%name = 'source'
         glbatts(iglb)%value = 'Coherens version '//TRIM(model_version)
      CASE (6)
         glbatts(iglb)%name = 'comment'
         glbatts(iglb)%value = comment_CF
      CASE (7)
         glbatts(iglb)%name = 'references'
         glbatts(iglb)%value = references_CF
      CASE (8)
         IF (iopt_CDF.EQ.1) THEN
            glbatts(iglb)%name = 'netcdf'
            glbatts(iglb)%value = cf90_inq_libvers()
         ENDIF
   END SELECT
ENDDO iglb_710


CALL log_timer_out()


RETURN

END SUBROUTINE reset_mod_params

!========================================================================

SUBROUTINE reset_mod_vars(varatts)
!************************************************************************
!
! *reset_mod_vars* Reset model variable attributes
!
! Author - Patrick Luyten
!
! Revision - @(COHERENS)reset_model.f90  V2.10.1
!
! Description - 
!
! Module calls - inquire_var
!
!************************************************************************
!
USE datatypes
USE modvars_routines, ONLY: inquire_var

!
!*Arguments
!
TYPE (VariableAtts), INTENT(INOUT) :: varatts

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*varatts* DERIVED Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cnum
INTEGER :: ivarid


procname(pglev+1) = 'reset_mod_vars'
CALL log_timer_in()

ivarid = varatts%ivarid
IF (ivarid.GT.0) THEN
   CALL inquire_var(ivarid,f90_name=varatts%f90_name,&
                  & standard_name=varatts%standard_name,&
                  & comment=varatts%comment,&
                  & long_name=varatts%long_name,units=varatts%units)
   IF (varatts%nrank.LT.0) CALL inquire_var(ivarid,nrank=varatts%nrank)
   IF (varatts%numvar.GT.0) THEN
      WRITE (cnum,'(I12)') varatts%numvar; cnum = ADJUSTL(cnum)
      varatts%f90_name = TRIM(varatts%f90_name)//'_'//TRIM(cnum)
      varatts%long_name = TRIM(varatts%long_name)//'_'//TRIM(cnum)
   ENDIF
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE reset_mod_vars

!========================================================================

SUBROUTINE reset_out_files_1d(filepars,arrname,outgpars,tseries)
!************************************************************************
!
! *reset_out_files_1d* Reset attributes of a user output file
!
! Author - Patrick Luyten
!
! Revision - @(COHERENS)reset_model.f90  V2.10.1
!
! Description -  argument has rank 1
!
! Reference -
!
! Calling program - harmonic_analysis_init, time_averages_init, time_series_init
!
! Module calls - warning_reset_arr_struc
!
!************************************************************************
!
USE datatypes
USE switches
USE error_routines, ONLY: warning_reset_arr_struc
!
!*Arguments
!
TYPE (FileParams), INTENT(INOUT), DIMENSION(:) :: filepars
CHARACTER (LEN=*), INTENT(IN) :: arrname
LOGICAL, INTENT(IN) :: tseries
TYPE (OutGridParams), INTENT(INOUT), DIMENSION(:) :: outgpars


!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*filepars*  DERIVED Attributes of user-output files
!*arrname*   CHAR    Name of the derived type array
!*outgpars*  CHAR    Attributes of corresponding output grid
!*tseries*   LOGICAL .TRUE. for time series, .FALSE. for time-averaged and
!                    harmonic output
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize


nsize = SIZE(filepars)
IF (nsize.EQ.0) RETURN

procname(pglev+1) = 'reset_out_files_1d'
CALL log_timer_in()

!---reset attributes
n_110: DO n=1,nsize
   IF (filepars(n)%defined) THEN
      filepars(n)%status = outgpars(n)%status
      filepars(n)%packing = outgpars(n)%packing
      filepars(n)%iunit = -1
   ENDIF
ENDDO n_110

!---reset attributes for verification procedure
IF (iopt_verif.EQ.1.AND.tseries) THEN
   CALL warning_reset_arr_struc(filepars(1)%form,arrname,'%form','N',1,(/1/))
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE reset_out_files_1d

!========================================================================

SUBROUTINE reset_out_files_2d(filepars,arrname,outgpars,tseries)
!************************************************************************
!
! *reset_out_files_2d* Reset attributes of a user output file
!
! Author - Patrick Luyten
!
! Revision - @(COHERENS)reset_model.f90  V2.4
!
! Description -  argument has rank 2
!
! Reference -
!
! Calling program - harmonic_analysis_init
!
! Module calls - warning_reset_arr_struc
!
!************************************************************************
!
USE datatypes
USE switches
USE error_routines, ONLY: warning_reset_arr_struc

!
!*Arguments
!
TYPE (FileParams), INTENT(INOUT), DIMENSION(:,:) :: filepars
CHARACTER (LEN=*), INTENT(IN) :: arrname
LOGICAL, INTENT(IN) :: tseries
TYPE (OutGridParams), INTENT(INOUT), DIMENSION(SIZE(filepars,DIM=1)) :: outgpars


!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*filepars*  DERIVED Attributes of user-output files
!*arrname*   CHAR    Name of the derived type array
!*outgpars*  CHAR    Attributes of corresponding output grid
!*tseries*   LOGICAL .TRUE. for time series, .FALSE. for time-averaged and
!                    harmonic output
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n1, n2
INTEGER, DIMENSION(2) :: ndims


IF (SIZE(filepars).EQ.0) RETURN

procname(pglev+1) = 'reset_out_files_2d'
CALL log_timer_in()

ndims = SHAPE(filepars)

!---reset attributes
n2_110: DO n2=1,ndims(2)
n1_110: DO n1=1,ndims(1)
   IF (filepars(n1,n2)%defined) THEN
      filepars(n1,n2)%status = outgpars(n1)%status
      filepars(n1,n2)%packing = outgpars(n1)%packing
      filepars(n1,n2)%iunit = -1
   ENDIF
ENDDO n1_110
ENDDO n2_110

!---reset attributes for verification procedure
IF (iopt_verif.EQ.1.AND.tseries) THEN
   n2_210: DO n2=1,ndims(2)
      CALL warning_reset_arr_struc(filepars(1,n2)%form,arrname,'%form','N',2,&
                                & (/1,n2/))
   ENDDO n2_210
ENDIF


CALL log_timer_out()


RETURN

END SUBROUTINE reset_out_files_2d

!========================================================================

SUBROUTINE reset_out_gpars(outgpars,arrname,tseries,end_reset)
!************************************************************************
!
! *reset_out_pars* Reset output parameters
!
! Author - Patrick Luyten
!
! Revision - @(COHERENS)reset_model.f90  V2.11
!
! Description - 
!
! Module calls - add_secs_to_date, convert_date, lim_dims,
!                warning_reset_arr_struc
!
!************************************************************************
!
USE datatypes
USE gridpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: warning_reset_arr_struc
USE time_routines, ONLY: add_secs_to_date, convert_date, num_time_steps
USE utility_routines, ONLY: lim_dims

!
!*Arguments
!
LOGICAL, INTENT(IN) :: end_reset, tseries
CHARACTER (LEN=*), INTENT(IN) :: arrname
TYPE (OutGridParams), INTENT(INOUT), DIMENSION(:) :: outgpars

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*outgpars* DERIVED  Array of output grid parameters
!*arrname*  CHAR     Array name
!*tseries*  LOGICAL  .TRUE. for time series, .FALSE. for time-averaged and
!                    harmonic output
!*end_reset* LOGICAL Resets end value of time resolution (%tlims(2)) to a
!                    value lower than or equal to nstep if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: gridded
INTEGER :: iset, limend, nodim, nosets, numsteps, time_format
INTEGER, DIMENSION(3) :: tlims
INTEGER, DIMENSION(7) :: intdate
REAL :: deltout


nosets = SIZE(outgpars)
IF (nosets.EQ.0) RETURN

procname(pglev+1) = 'reset_out_gpars'
CALL log_timer_in()

!
!1. Set output dates
!-------------------
!

iset_110: DO iset=1,nosets

!  ---reset end limit (if necessary)
   tlims = outgpars(iset)%tlims
   IF (end_reset.AND.tlims(2).GT.nstep) THEN
      limend = tlims(1) + ((nstep-tlims(1))/tlims(3))*tlims(3)
      CALL warning_reset_arr_struc(tlims(2),arrname,'%tlims(2)',limend,1,&
                                & (/iset/))
      outgpars(iset)%tlims = tlims
   ENDIF
   tlims = outgpars(iset)%tlims

!  ---start date
   IF (tseries) THEN
      numsteps = 2*tlims(1)
   ELSE
      numsteps = 2*tlims(1) + tlims(3)
   ENDIF
   CALL add_secs_to_date(CStartDateTime,outgpars(iset)%startdate,numsteps,&
                       & 0.5*delt2d)

!  ---reference date
   IF (outgpars(iset)%refdate.EQ.cdatetime_undef) THEN
      intdate = convert_date(outgpars(iset)%startdate)
   ELSE
      intdate = convert_date(outgpars(iset)%refdate)
   ENDIF
   SELECT CASE (outgpars(iset)%time_format)
      CASE (1); intdate(7) = 0
      CASE (2); intdate(6:7) = 0
      CASE (3); intdate(5:7) = 0
      CASE (4); intdate(4:7) = 0
      CASE (5); intdate(4:7) = 0
      CASE (6); intdate(3:7) = 0
      CASE (7); intdate(2:7) = 0
      CASE (8); intdate(2:7) = 0
   END SELECT
   outgpars(iset)%refdate = convert_date(intdate,seps='- ::')

!  ---end date
   IF (tseries) THEN
      numsteps = 2*tlims(2)
   ELSE
      numsteps = 2*tlims(2) - tlims(3)
   ENDIF
   CALL add_secs_to_date(CStartDateTime,outgpars(iset)%enddate,numsteps,&
                       & 0.5*delt2d)

ENDDO iset_110

!
!2. Grid and time dimensions
!---------------------------
!

iset_210: DO iset=1,nosets
   nodim = outgpars(iset)%nodim; gridded = outgpars(iset)%gridded
   IF (nodim.EQ.2) outgpars(iset)%zlims = 1
   SELECT CASE (nodim)
   CASE (0)
      outgpars(iset)%xlims = 0; outgpars(iset)%ylims = 0
      outgpars(iset)%zlims = 0; outgpars(iset)%nostats = 0
   CASE (2,3)   
      outgpars(iset)%ncout = MERGE(lim_dims(outgpars(iset)%xlims),0,gridded)
      outgpars(iset)%nrout = MERGE(lim_dims(outgpars(iset)%ylims),0,gridded)
   END SELECT
   outgpars(iset)%nzout = MERGE(lim_dims(outgpars(iset)%zlims),0,nodim.EQ.3)
   time_format = outgpars(iset)%time_format
   tlims = outgpars(iset)%tlims
   outgpars(iset)%nstepout = (tlims(2)-tlims(1))/tlims(3) + 1
   deltout = tlims(3)*delt2d
   CALL num_time_steps(outgpars(iset)%startdate,outgpars(iset)%enddate,&
                     & deltout,outgpars(iset)%nstepout)
   outgpars(iset)%nstepout = outgpars(iset)%nstepout + 1
   IF (time_format.GT.0) THEN
      outgpars(iset)%deltout = deltout/time_convert(time_format)
   ENDIF
ENDDO iset_210

!
!3. Packing
!----------
!

iset_310: DO iset=1,nosets
   IF (outgpars(iset)%nodim.EQ.0.OR.(.NOT.outgpars(iset)%gridded)) THEN
      CALL warning_reset_arr_struc(outgpars(iset)%packing,arrname,&
                               & '%packing',.FALSE.,1,(/iset/))
   ENDIF
ENDDO iset_310

!
!4. Stations
!-----------
!

iset_410: DO iset=1,nosets
   IF (outgpars(iset)%gridded) THEN
      CALL warning_reset_arr_struc(outgpars(iset)%nostats,arrname,'%nostats',&
                                 & 0,1,(/iset/))
   ENDIF
ENDDO iset_410

!
!5. Time-dependent grid
!----------------------
!

iset_510: DO iset=1,nosets
   IF (outgpars(iset)%nodim.LT.3.OR.(.NOT.tseries)) THEN
      CALL warning_reset_arr_struc(outgpars(iset)%time_grid,arrname,&
                                & '%time_grid',.FALSE.,1,(/iset/))
   ENDIF
ENDDO iset_510

!
!6. Reset attributes for verification procedure
!----------------------------------------------
!

IF (tseries.AND.iopt_verif.EQ.1.AND.iset.EQ.1) THEN

   CALL warning_reset_arr_struc(outgpars(iset)%gridded,arrname,&
                            & '%gridded',.TRUE.,1,(/iset/))
   CALL warning_reset_arr_struc(outgpars(iset)%packing,arrname,&
                            & '%packing',.FALSE.,1,(/iset/))
   CALL warning_reset_arr_struc(outgpars(iset)%tlims(1),arrname,&
                            & '%tlims(1)',0,1,(/iset/))
   CALL warning_reset_arr_struc(outgpars(iset)%tlims(2),arrname,&
                            & '%tlims(2)',nstep,1,(/iset/))

ENDIF


!
!7. Bedvars
!-----------
!

IF (iopt_morph.EQ.0) THEN
   iset_710: DO iset=1,nosets
      CALL warning_reset_arr_struc(outgpars(iset)%bedvars,arrname,&
                               & '%gridded',.FALSE.,1,(/iset/))
   ENDDO iset_710
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE reset_out_gpars

!========================================================================

SUBROUTINE reset_out_stats(statlocs)
!************************************************************************
!
! *reset_out_stats* Reset station data
!
! Author - Patrick Luyten
!
! Revision - @(COHERENS)reset_model.f90  V2.0
!
! Description - 
!
! Module calls - add_secs_to_date, convert_date, lim_dims,
!                warning_reset_arr_struc
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
TYPE (StationLocs), INTENT(INOUT), DIMENSION(:) :: statlocs

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*outstats* DERIVED Array of output station attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cstat
INTEGER :: istat, nostats


nostats = SIZE(statlocs)
IF (nostats.EQ.0) RETURN

procname(pglev+1) = 'reset_out_stats'
CALL log_timer_in()

istat_110: DO istat=1,nostats
   IF (LEN_TRIM(statlocs(istat)%name).EQ.0) THEN
      WRITE (cstat,'(I12)') istat
      statlocs(istat)%name = 'Station '//TRIM(ADJUSTL(cstat))
   ENDIF
ENDDO istat_110

CALL log_timer_out()


END SUBROUTINE reset_out_stats

!========================================================================

SUBROUTINE reset_out_vars(varatts)
!************************************************************************
!
! *reset_out_vars* Reset output variable attributes
!
! Author - Patrick Luyten
!
! Revision - @(COHERENS)reset_model.f90  V2.7.1
!
! Description - 
!
! Module calls - inquire_var
!
!************************************************************************
!
USE datatypes
USE syspars
USE modvars_routines, ONLY: inquire_var

!
!*Arguments
!
TYPE (VariableAtts), INTENT(INOUT), DIMENSION(:) :: varatts

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*varatts* DERIVED Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cvar
INTEGER :: ivar, ivarid, novars, numvar


novars = SIZE(varatts)
IF (novars.EQ.0) RETURN

procname(pglev+1) = 'reset_out_vars'
CALL log_timer_in()

ivar_100: DO ivar=1,novars
   ivarid = varatts(ivar)%ivarid
   IF (ivarid.GT.0) THEN
      numvar = varatts(ivar)%numvar
      IF (numvar.GT.0) THEN
         WRITE (cvar,'(I12)') numvar; cvar = ADJUSTL(cvar)
      ENDIF
      IF (TRIM(varatts(ivar)%f90_name).EQ.'') THEN
         CALL inquire_var(ivarid,f90_name=varatts(ivar)%f90_name)
         IF (numvar.GT.0) THEN
            varatts(ivar)%f90_name = TRIM(varatts(ivar)%f90_name)//'_'//&
                                   & TRIM(cvar)
         ENDIF
      ENDIF
      IF (TRIM(varatts(ivar)%standard_name).EQ.'') THEN
         CALL inquire_var(ivarid,standard_name=varatts(ivar)%standard_name)
         IF (numvar.GT.0) THEN
            varatts(ivar)%standard_name = TRIM(varatts(ivar)%standard_name)//&
                                        & '_'//TRIM(cvar)
         ENDIF
      ENDIF
      IF (TRIM(varatts(ivar)%long_name).EQ.'') THEN
         CALL inquire_var(ivarid,long_name=varatts(ivar)%long_name)
         IF (numvar.GT.0) THEN
            varatts(ivar)%long_name = TRIM(varatts(ivar)%long_name)//'_'//&
                                    & TRIM(cvar)
         ENDIF
      ENDIF
      IF (TRIM(varatts(ivar)%vector_name).EQ.'') THEN
         CALL inquire_var(ivarid,vector_name=varatts(ivar)%vector_name)
      ENDIF
      IF (TRIM(varatts(ivar)%units).EQ.'') THEN
         CALL inquire_var(ivarid,units=varatts(ivar)%units)
      ENDIF
      IF (varatts(ivar)%nrank.LT.0) THEN
         CALL inquire_var(ivarid,nrank=varatts(ivar)%nrank)
      ENDIF
!      IF (varatts(ivar)%fill_value.EQ.0.0_kndrlong) THEN
!         CALL inquire_var(ivarid,fill_value=varatts(ivar)%fill_value)
!      ENDIF
!      varatts(ivar)%fill = .TRUE.
!      varatts(ivar)%fill_value  = float_fill
   ENDIF
ENDDO ivar_100

CALL log_timer_out()


RETURN

END SUBROUTINE reset_out_vars


END MODULE reset_model
