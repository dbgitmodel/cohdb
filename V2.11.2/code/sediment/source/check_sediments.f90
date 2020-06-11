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

MODULE check_sediments
!************************************************************************
!
! *check_sed_params* Check sediment model parameters and arrays for errors
!
! Author - Boudewijn Decrop and Alexander Breugem (IMDC)
!
! Version - @(COHERENS)check_sediments.f90  V2.11.1
!
! $Date: 2018-05-23 10:09:17 +0200 (Wed, 23 May 2018) $
!
! $Revision: 1138 $
!
! Description - 
!
! Reference -
!
! Routines - check_dar_params, check_dar_spec, check_morphics,
!            check_morph_params, check_sedics, check_sed_params
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS

!=====================================================================

SUBROUTINE check_sed_params
!************************************************************************
!
! *check_sed_params* Check parameters for the sediment transport model
!
! Author - Boudewijn Decrop and Alexander Breugem (IMDC)
!
! Version - @(COHERENS)check_sediments.f90  V2.11.1
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - check_mult_counter, error_abort, error_lbound_var,
!                error_limits_var, error_value_var
!
!************************************************************************
!
USE iopars
USE paralpars
USE sedpars
USE sedswitches
USE switches
USE timepars
USE error_routines, ONLY: check_mult_counter, error_abort, error_lbound_var, &
                        & error_limits_var, error_ubound_var, error_value_var
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
REAL, PARAMETER :: real_small = 1.0E-5


IF (.NOT.master) RETURN

procname(pglev+1) = 'check_sed_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!
!1.1 Ranges
!----------
!
!---type of transport
CALL error_limits_var(iopt_sed_bedeq,'iopt_sed_bedeq',0,7)
CALL error_limits_var(iopt_sed_mode,'iopt_sed_mode',1,4)
CALL error_limits_var(iopt_sed_nodim,'iopt_sed_nodim',2,3)
CALL error_limits_var(iopt_sed_toteq,'iopt_sed_toteq',0,7)
CALL error_limits_var(iopt_sed_type,'iopt_sed_type',1,2)

!---bottom boundary condition
CALL error_limits_var(iopt_sed_bbc,'iopt_sed_bbc',0,2)
CALL error_limits_var(iopt_sed_bbc_eq,'iopt_sed_bbc_eq',1,6)
CALL error_limits_var(iopt_sed_bbc_ref,'iopt_sed_bbc_ref',1,2)

!---sediment density
CALL error_limits_var(iopt_sed_dens_grad,'iopt_sed_dens_grad',0,1)

!---sediment diffusion
CALL error_limits_var(iopt_sed_beta,'iopt_sed_beta',1,3)
CALL error_limits_var(iopt_sed_wave_diff,'iopt_sed_wave_diff',0,1)

!---bed stress
CALL error_limits_var(iopt_sed_bstres,'iopt_sed_bstres_cr',0,1)
CALL error_limits_var(iopt_sed_bstres_cr,'iopt_sed_bstres_cr',1,4)
CALL error_limits_var(iopt_sed_hiding,'iopt_sed_hiding',0,2)
CALL error_limits_var(iopt_sed_rough,'iopt_sed_rough',1,3)

!---settling velocity
CALL error_limits_var(iopt_sed_ws_floc,'iopt_sed_ws_floc',0,3)
CALL error_limits_var(iopt_sed_ws_hindset,'iopt_sed_ws_hindset',0,2)
CALL error_limits_var(iopt_sed_ws,'iopt_sed_ws',1,6)
CALL error_limits_var(iopt_sed_ws_lim,'iopt_sed_ws_lim',0,1)
CALL error_limits_var(iopt_sed_vadv,'iopt_sed_vadv',0,3)

!---slope effects
CALL error_limits_var(iopt_sed_slope,'iopt_sed_slope',0,1)

!---open boundary conditions
CALL error_limits_var(iopt_sed_obc_flux,'iopt_sed_obc_flux',0,1)

!---numerical
CALL error_limits_var(iopt_sed_filter,'iopt_sed_filter',0,1)
CALL error_limits_var(iopt_sed_median,'iopt_sed_median',1,2)
CALL error_limits_var(iopt_sed_source_impl,'iopt_sed_source_impl',0,3)

!
!1.2 Incompatibilities
!---------------------
!
!---bed and total load
IF (iopt_sed_mode.EQ.1.OR.iopt_sed_mode.EQ.3) THEN
   CALL error_limits_var(iopt_sed_bedeq,'iopt_sed_bedeq',1,7)
ELSEIF (iopt_sed_mode.EQ.4) THEN
   CALL error_limits_var(iopt_sed_toteq,'iopt_sed_toteq',1,7)   
ENDIF

!---only quadratic bottom friction is allowed for recalculation of bed stresses
IF (iopt_bstres_form.GT.0.AND.iopt_sed_rough.EQ.1) THEN
   CALL error_value_var(iopt_bstres_form,'iopt_bstres_form',2)
ENDIF

!---no morphology in absence of bed or total load
IF (iopt_sed.EQ.2) CALL error_value_var(iopt_morph,'iopt_morph',0)

!
!2. Model Parameters
!-------------------
!
!---time counter
SELECT CASE (iopt_sed_mode)
   CASE (1,4); CALL check_mult_counter(icsed,'icsed',ic3d)
   CASE (2,3); 
      IF (iopt_grid_nodim.NE.2.AND.iopt_hydro_impl.EQ.0) THEN
         CALL error_value_var(icsed,'icsed',ic3d)
      ENDIF
END SELECT

IF (iopt_sed_nodim.EQ.2.AND.iopt_sed_bbc_eq.EQ.1) THEN
   CALL error_lbound_var(nrquad_sed,'nrquad_sed',3,.TRUE.)
ENDIF

IF (iopt_waves.EQ.1) THEN
   CALL error_lbound_var(nrquad_wav,'nrquad_wav',3,.TRUE.)
ENDIF
CALL error_lbound_var(beta_sed_cst,'beta_sed_cst',0.0,.FALSE.)

CALL error_lbound_var(n_RichZaki,'n_RichZaki',0.0,.FALSE.)
CALL error_lbound_var(parth_exp,'parth_exp',0.0,.FALSE.)
CALL error_lbound_var(parth_coef,'parth_coef',0.0,.FALSE.)

CALL error_lbound_var(z0_coef,'z0_coef',1.0,.TRUE.)

CALL error_lbound_var(a_leussen,'a_leussen',0.0,.FALSE.)
CALL error_lbound_var(b_leussen,'b_leussen',0.0,.FALSE.)
CALL error_lbound_var(alpha_VR,'alpha_VR',0.0,.FALSE.)

CALL error_lbound_var(floc_VR_min,'floc_VR_min',0.0,.FALSE.)
CALL error_lbound_var(floc_VR_max,'floc_VR_max',floc_VR_min,.FALSE.)

CALL error_lbound_var(zrough_grain,'zrough_grain',0.0,.FALSE.)
CALL error_limits_var(height_c_cst,'height_c_cst',real_small,1.0)

CALL error_lbound_var(coef_bed_grad,'coef_bed_grad',0.0,.FALSE.)

IF (iopt_sed_ws_floc.GE.2.OR.iopt_sed_ws_hindset.EQ.2.OR.&
 & (iopt_sed_bbc.EQ.3.AND.iopt_sed_nodim.EQ.2)) THEN
   CALL error_limits_var(cgel,'cgel',real_small,1.0)
ENDIF

CALL error_limits_var(cmax,'cmax',real_small,1.0)

IF (iopt_sed_mode.EQ.2.OR.iopt_sed_mode.EQ.3) THEN
   IF (iopt_grid_nodim.EQ.2) THEN
      CALL  error_value_var(iopt_sed_nodim,'iopt_sed_nodim',2)
   ENDIF
   IF (iopt_grid_nodim.EQ.3) THEN
      CALL  error_value_var(iopt_sed_nodim,'iopt_sed_nodim',3)
   ENDIF
ENDIF

CALL error_abort('check_sed_params',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_sed_params

!=====================================================================

SUBROUTINE check_morph_params
!************************************************************************
!
! *check_morph_params* Check parameters for the morphological model
!
! Author - Boudewijn Decrop and Alexander Breugem (IMDC)
!
! Version - @(COHERENS)check_sediments.f90  V2.10
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - check_date, check_mult_counter, convert_date, error_abort,
!                error_lbound_var, error_lbound_var_date, error_limits_var,
!                error_value_var, num_time_steps
!
!************************************************************************
!
USE iopars
USE paralpars
USE morphpars
USE morphswitches
USE sedpars
USE sedswitches
USE switches
USE timepars
USE error_routines, ONLY: check_mult_counter, error_abort, error_lbound_var, &
                        & error_limits_var, error_value_var
USE time_routines, ONLY: check_date, convert_date, error_lbound_var_date, &
                       & log_timer_in, log_timer_out, num_time_steps

IMPLICIT NONE


IF (.NOT.master) RETURN

procname(pglev+1) = 'check_morph_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!
!1.1  Ranges
!-----------
!

CALL error_limits_var(iopt_morph_active_layer,'iopt_morph_active_layer',0,3)
CALL error_limits_var(iopt_morph_avalanching,'iopt_morph_avalanching',0,1)
CALL error_limits_var(iopt_morph_corr,'iopt_morph_corr',0,2)
CALL error_limits_var(iopt_morph_fixed_layer,'iopt_morph_fixed_layer',0,1)
CALL error_limits_var(iopt_morph_hdiff,'iopt_morph_hdiff',0,1)
CALL error_limits_var(iopt_morph_tidal_scheme,'iopt_morph_tidal_scheme',0,3)
CALL error_limits_var(iopt_morph_time_int,'iopt_morph_time_int',1,3)
CALL error_limits_var(iopt_morph_vert_fluxes,'iopt_morph_vert_fluxes',0,1)

!
!1.2 Incompatibilities
!---------------------
!
!---no morphology without sediment
CALL error_value_var(iopt_sed,'iopt_sed',1)

!
!2. Model Parameters
!-------------------
!
!---time counter
SELECT CASE (iopt_sed_mode)
   CASE (1,4); CALL error_value_var(icmorph,'icmorph',icsed)
   CASE (2,3); CALL check_mult_counter(icmorph,'icmorph',icsed)
END SELECT

IF (iopt_morph_avalanching.EQ.1) THEN
   CALL check_mult_counter(icaval,'icaval',icmorph)
ENDIF

!---parameters which must be positive 
CALL error_lbound_var(average_bedform_height,'average_bedform_height',&
                    & 0.0,.FALSE.)
CALL error_lbound_var(average_bedform_length,'average_bedform_length',&
                    & 0.0,.FALSE.)
CALL error_limits_var(bed_porosity_cst,'bed_porosity_cst',0.0,1.0)
CALL error_lbound_var(dune_celerity,'dune_celerity',0.0,.FALSE.)
CALL error_lbound_var(k1_harris,'k1_harris',0.0,.FALSE.)
CALL error_lbound_var(k2_harris,'k2_harris',0.0,.FALSE.)
CALL error_lbound_var(k2_jameson,'k2_jameson',0.0,.FALSE.)
CALL error_lbound_var(morph_factor,'morph_factor',0.0,.FALSE.)
CALL error_lbound_var(similarity_range,'similarity_range',0.0,.FALSE.)
CALL error_lbound_var(trough_probability,'trough_probability',0.0,.FALSE.)

!---critical slope angles for avalanching
CALL error_limits_var(dyn_angle_wet,'dyn_angle_wet',0.0,90.0)
CALL error_limits_var(dyn_angle_dry,'dyn_angle_dry',0.0,90.0)
CALL error_limits_var(stat_angle_wet,'stat_angle_wet',0.0,90.0)
CALL error_limits_var(stat_angle_dry,'stat_angle_dry',0.0,90.0)

!---date/time tidal averaging
IF (iopt_tidal_accel.EQ.1) THEN
   CALL error_lbound_var(morph_steps,'morph_steps',1,.TRUE.)
   CALL error_lbound_var(nstep_hydro,'nstep_hydro',1,.TRUE.)
   CALL error_lbound_var(number_tidal_steps,'number_tidal_steps',1,.TRUE.)
ENDIF

!---iteration parameters
CALL error_lbound_var(max_iterations_aval,'max_iterations_aval',1,.TRUE.)

CALL error_abort('check_morph_params',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_morph_params

!========================================================================

SUBROUTINE check_dar_params
!************************************************************************
!
! *check_dar_params* Check parameters for dredging/relocation
!
! Author - Boudewijn Decrop and Alexander Breugem (IMDC)
!
! Version - @(COHERENS)check_sediments.f90  V2.8
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - error_abort, error_lbound_var, error_limits_var,
!                error_value_var
!
!************************************************************************
!
USE darpars
USE darswitches
USE iopars
USE paralpars
USE switches
USE error_routines, ONLY: error_abort, error_lbound_var, error_limits_var, &
                        & error_value_var
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


IF (.NOT.master) RETURN

procname(pglev+1) = 'check_dar_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!
!---ranges
CALL error_limits_var(iopt_dar_coupling_time,'iopt_dar_coupling_time',1,3)
CALL error_limits_var(iopt_dar_coupling_space,'iopt_dar_coupling_space',1,3)
CALL error_limits_var(iopt_dar_dredging_criterion,&
                   & 'iopt_dar_dredging_criterion',1,3)
CALL error_limits_var(iopt_dar_duration,'iopt_dar_duration',1,3)
CALL error_limits_var(iopt_dar_effect,'iopt_dar_effect',1,3)
CALL error_limits_var(iopt_dar_relocation_criterion,&
                   & 'iopt_dar_relocation_criterion',0,2)
CALL error_limits_var(iopt_dar_scenario,'iopt_dar_scenario',1,3)
CALL error_limits_var(iopt_dar_time_dredge,'iopt_dar_time_dredge',1,3)
CALL error_limits_var(iopt_dar_time_relocate,'iopt_dar_time_relocate',1,3)
CALL error_limits_var(iopt_dar_user_volume,'iopt_dar_user_volume',1,3)

!---incompatibilities
SELECT CASE (iopt_dar_effect)
   CASE (1); CALL error_value_var(iopt_morph,'iopt_morph',1)
   CASE (2); CALL error_value_var(iopt_sed,'iopt_sed',1)
   CASE (3)
      CALL error_value_var(iopt_sed,'iopt_sed',1)
      CALL error_value_var(iopt_morph,'iopt_morph',1)
END SELECT

!
!2. Parameters
!-------------
!

CALL error_lbound_var(maxtrackpoints,'maxtrackpoints',1,.TRUE.)
CALL error_lbound_var(nocampaigns,'nocampaigns',1,.TRUE.)
CALL error_lbound_var(nocirclepoints,'nocirclepoints',1,.TRUE.)
CALL error_lbound_var(nodredgingsites,'nodredgingsites',1,.TRUE.)
CALL error_lbound_var(norelocationsites,'norelocationsites',1,.TRUE.)
CALL error_lbound_var(noships,'noships',1,.TRUE.)

CALL error_abort('check_dar_params',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_dar_params

!========================================================================

SUBROUTINE check_dar_spec
!************************************************************************
!
! *check_dar_spec* Check attribute arrays for dredging/relocation
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)check_sediments.f90  V2.8
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - error_abort, error_lbound_arr_struc, error_limits_arr_struc
!
!************************************************************************
!
USE dararrays
USE darpars
USE darswitches
USE iopars
USE paralpars
USE sedpars
USE error_routines, ONLY: error_abort, error_lbound_arr_struc, &
                       &  error_limits_arr_struc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: f, l, m, n, noerrs
REAL :: sumfrac
CHARACTER (LEN=12) :: cval, cval2
CHARACTER (LEN=27) :: cvar

IF (.NOT.master) RETURN

procname(pglev+1) = 'check_dar_spec'
CALL log_timer_in()

!
!1. Campaign attributes
!----------------------
!

IF (iopt_dar_dredging_criterion.EQ.3) THEN

   n_110: DO n=1,nocampaigns
      
      CALL error_limits_arr_struc(Campaigns(n)%nr_dredgingsite,'Campaigns',&
                             & '%nr_dredgingsite',1,nodredgingsites,1,(/n/))
      CALL error_limits_arr_struc(Campaigns(n)%nr_relocationsite,'Campaigns',&
                             & '%nr_relocationsite',1,norelocationsites,1,(/n/))
      CALL error_limits_arr_struc(Campaigns(n)%nr_ship,'Campaigns','%nr_ship',&
                             & 1,noships,1,(/n/))
      CALL error_lbound_arr_struc(Campaigns(n)%t_lag,'Campaigns','%t_lag',1,&
                             & .TRUE.,1,(/n/))
      CALL error_lbound_arr_struc(Campaigns(n)%t_spread,'Campaigns',&
                             & '%t_spread',1,.TRUE.,1,(/n/))
      CALL error_lbound_arr_struc(Campaigns(n)%campdepth,'Campaigns',&
                             & '%campdepth',0.0,.TRUE.,1,(/n/))
      CALL error_lbound_arr_struc(Campaigns(n)%camplevel,'Campaigns',&
                             & '%camplevel',0.0,.TRUE.,1,(/n/))
      CALL error_lbound_arr_struc(Campaigns(n)%campvolume,'Campaigns',&
                             & '%campvolume',0.0,.TRUE.,1,(/n/))

      IF (.NOT.Campaigns(n)%dredge) THEN
!        ---sum of the campaign volume fractions must be equal to one
         sumfrac  = SUM(Campaigns(n)%campfractions)
         IF (sumfrac.NE.1.0) THEN
            nerrs = nerrs + 1
            IF (errchk.AND.nerrs.LE.maxerrors) THEN
               WRITE (cval,'(G12.7)') sumfrac
               WRITE (ioerr,'(A)') 'Invalid sum of fractions in real array'//&
                                   ' Campaigns%campfractions'//': '//TRIM(cval)
               WRITE (ioerr,'(A)') 'Sum over all fractions must be equal to one'
            ENDIF
         ENDIF
      ENDIF
      
      noerrs = COUNT((Campaigns(n)%relocationsites(1:MaxSites).LT.1.OR.&
               & Campaigns(n)%relocationsites(1:MaxSites).GT.norelocationsites))
      IF (noerrs.GT.0) THEN
         m_111: DO m=1,MaxSites
         WRITE (cval,'(I12)') m; cvar = '('//TRIM(ADJUSTL(cval))//')'
            CALL error_limits_arr_struc(Campaigns(n)%relocationsites(m),&
                                 & 'Campaigns','%relocationsites'//TRIM(cvar), &
                                 & 1,norelocationsites,1,(/n/))
         ENDDO m_111
      ENDIF

   ENDDO n_110

ENDIF

!
!2. Dredging sites attributes
!----------------------------
!

n_210: DO n=1,nodredgingsites

   CALL error_lbound_arr_struc(DredgingSites(n)%monitoring_frequency,&
                     & 'DredgingSites','%monitoring_frequency',1,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(DredgingSites(n)%nr_coordpol,'DredgingSites',&
                     & '%nr_coordpol',1,.TRUE.,1,(/n/))
   CALL error_limits_arr_struc(DredgingSites(n)%nr_ship,'DredgingSites',&
                     & '%nr_ship',1,noships,1,(/n/))
   CALL error_lbound_arr_struc(DredgingSites(n)%t_lag,'DredgingSites',&
                     & '%t_lag',1,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(DredgingSites(n)%critical_volume,&
                     & 'DredgingSites','%critical_volume',0.0,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(DredgingSites(n)%over_depth,'DredgingSites',&
                     & '%over_depth',0.0,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(DredgingSites(n)%target_depth,'DredgingSites',&
                     & '%target_dept',0.0,.TRUE.,1,(/n/))

   noerrs = COUNT(DredgingSites(n)%relocationsites(1:MaxSites).LT.1.OR. &
           & DredgingSites(n)%relocationsites(1:MaxSites).GT.norelocationsites)
   IF (noerrs.GT.0) THEN
      m_211: DO m=1,MaxSites
         WRITE (cval,'(I12)') m; cvar = '('//TRIM(ADJUSTL(cval))//')'
         CALL error_limits_arr_struc(DredgingSites(n)%relocationsites(m),&
                             & 'DredgingSites','%relocationsites'//TRIM(cvar),&
                             & 1,norelocationsites,1,(/n/))
      ENDDO m_211
   ENDIF

   noerrs = COUNT(DredgingSites(n)%priority_order(1:MaxSubzones).LT.1.OR. &
             & DredgingSites(n)%relocationsites(1:MaxSites).GT.nodredgingsites)
   IF (noerrs.GT.0) THEN
      m_212: DO m=1,MaxSubzones
         WRITE (cval,'(I12)') m; cvar = '('//TRIM(ADJUSTL(cval))//')'
         CALL error_limits_arr_struc(DredgingSites(n)%priority_order(m),&
                            & 'DredgingSites','%priority_order'//TRIM(cvar),&
                            & 1,nodredgingsites,1,(/n/))
      ENDDO m_212
   ENDIF

   noerrs = COUNT(DredgingSites(n)%Xpol(1:MaxCoordPol).LT.0.0 )
   IF (noerrs.GT.0) THEN
      m_213: DO m=1,MaxCoordPol
         WRITE (cval,'(I12)') m; cvar = '('//TRIM(ADJUSTL(cval))//')'
         CALL error_lbound_arr_struc(DredgingSites(n)%Xpol(m),'DredgingSites',&
                                  & '%Xpol'//TRIM(cvar), 0.0,.TRUE.,1,(/n/))
      ENDDO m_213
   ENDIF

   noerrs = COUNT(DredgingSites(n)%Ypol(1:MaxCoordPol).LT.0.0 )
   IF (noerrs.GT.0) THEN
      m_214: DO m=1,MaxCoordPol
         WRITE (cval,'(I12)') m; cvar = '('//TRIM(ADJUSTL(cval))//')'
         CALL error_lbound_arr_struc(DredgingSites(n)%Ypol(m),'DredgingSites',&
                                  & '%Ypol'//TRIM(cvar),0.0,.TRUE.,1,(/n/))
      ENDDO m_214
   ENDIF

ENDDO n_210

!
!3. Relocation sites attributes
!------------------------------
!

n_310: DO n=1,norelocationsites

   CALL error_lbound_arr_struc(RelocationSites(n)%nr_coordpol,&
                            & 'RelocationSites','%nr_coordpol',1,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(RelocationSites(n)%max_volume,&
                           & 'RelocationSites','%max_volume',0.0,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(RelocationSites(n)%min_depth,'RelocationSites',&
                            & '%min_depth',0.0,.TRUE.,1,(/n/))

   noerrs = COUNT(RelocationSites(n)%Xpol(1:MaxCoordPol).LT.0.0 )
   IF (noerrs.GT.0) THEN
      m_311: DO m=1,MaxCoordPol
         WRITE (cval,'(I12)') m; cvar = '('//TRIM(ADJUSTL(cval))//')'
         CALL error_lbound_arr_struc(RelocationSites(n)%Xpol(m),&
                     & 'RelocationSites','%Xpol'//TRIM(cvar),0.0,.TRUE.,1,(/n/))
      ENDDO m_311
   ENDIF

   noerrs = COUNT(RelocationSites(n)%Ypol(1:MaxCoordPol).LT.0.0 )
   IF (noerrs.GT.0) THEN
      m_312: DO m=1,MaxCoordPol
         WRITE (cval,'(I12)') m; cvar = '('//TRIM(ADJUSTL(cval))//')'
         CALL error_lbound_arr_struc(RelocationSites(n)%Ypol(m),&
                     & 'RelocationSites','%Ypol'//TRIM(cvar),0.0,.TRUE.,1,(/n/))
      ENDDO m_312
   ENDIF

ENDDO n_310

!
!4. Ships attributes
!-------------------
!

n_410: DO n=1,noships

   CALL error_lbound_arr_struc(Ships(n)%t_fill,'Ships','%t_fill',0,&
                            & .TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%t_idle,'Ships','%t_idle',0,&
                            & .TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%t_relocate,'Ships','%t_relocate',0,&
                            & .TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%t_sail_dredge,'Ships','%t_sail_dredge',&
                             & 0,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%t_sail_relocate,'Ships',&
                            & '%t_sail_relocate',0,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%t_wait,'Ships','%t_wait',0,&
                            & .TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%D_dredge_max,'Ships','%D_dredge_max',&
                            & 0.0,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%D_dredge_min,'Ships','%D_dredge_min',&
                            & 0.0,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%D_relocate_max,'Ships',&
                            & '%D_relocate_max',0.0,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%D_relocate_min,'Ships',&
                            & '%D_relocate_min',0.0,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%H_dredge_max,'Ships','%H_dredge_max',&
                             & 0.0,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%H_dredge_min,'Ships','%H_dredge_min',&
                             & 0.0,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%H_relocate_max,'Ships',&
                            & '%H_relocate_max',0.0,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%H_relocate_min,'Ships',&
                            & '%H_relocate_min',0.0,.TRUE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%R_circle,'Ships','%R_circle',0.0,&
                            & .FALSE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%U_ship,'Ships','%U_ship',0.0,&
                            & .FALSE.,1,(/n/))
   CALL error_lbound_arr_struc(Ships(n)%V_ship,'Ships','%V_ship',0.0,&
                            & .FALSE.,1,(/n/))

   noerrs = COUNT(Ships(n)%p_overflow.LT.0.0.OR.Ships(n)%p_overflow.GT.1.0)
   IF (noerrs.GT.0) THEN
      f_411: DO f=1,nf
         WRITE (cval,'(I12)') f; cvar = '('//TRIM(ADJUSTL(cval))//')'
         CALL error_limits_arr_struc(Ships(n)%p_overflow(f),'Ships',&
                             & '%p_overflow'//TRIM(cvar),0.0, 1.0,1,(/n/))
      ENDDO f_411
   ENDIF

   noerrs = COUNT(Ships(n)%p_passive.LT.0.0.OR.Ships(n)%p_relocation.GT.1.0)
   IF (noerrs.GT.0) THEN
      f_412: DO f=1,nf
         WRITE (cval,'(I12)') f; cvar = '('//TRIM(ADJUSTL(cval))//')'
         CALL error_limits_arr_struc(Ships(n)%p_passive(f),'Ships',&
                                  & '%p_passive'//TRIM(cvar),0.0, 1.0,1,(/n/))
      ENDDO f_412
   ENDIF

   noerrs = COUNT(Ships(n)%p_relocation.LT.0.0.OR.Ships(n)%p_relocation.GT.1.0)
   IF (noerrs.GT.0) THEN
      f_413: DO f=1,nf
         WRITE (cval,'(I12)') f; cvar = '('//TRIM(ADJUSTL(cval))//')'
         CALL error_limits_arr_struc(Ships(n)%p_relocation(f),'Ships',&
                                 & '%p_relocation'//TRIM(cvar),0.0,1.0,1,(/n/))
      ENDDO f_413
   ENDIF

   noerrs = COUNT(Ships(n)%p_suctionhead.LT.0.0.OR.&
                & Ships(n)%p_suctionhead.GT.1.0)
   IF (noerrs.GT.0) THEN
      f_414: DO f=1,nf
         WRITE (cval,'(I12)') f; cvar = '('//TRIM(ADJUSTL(cval))//')'
         CALL error_limits_arr_struc(Ships(n)%p_suctionhead(f),'Ships',&
                            & '%p_suctionhead'//TRIM(cvar),0.0, 1.0,1,(/n/))
      ENDDO f_414
   ENDIF

ENDDO n_410

!
!5. Subzone dredging attributes
!-------------------------------
!

IF (iopt_dar_time_dredge.EQ.2) THEN

   n_510: DO n=1,nodredgingsites

      CALL error_lbound_arr_struc(SubzonesDredging(n)%nr_subzones,&
                           & 'SubzonesDredging','%nr_subzones',1,.TRUE.,1,(/n/))

      noerrs = COUNT(SubzonesDredging(n)%nr_coordpol(1:MaxSubzones).LT.1)
      IF (noerrs.GT.0) THEN
         m_511: DO m=1,MaxSubzones
            WRITE (cval,'(I12)') m; cvar = '('//TRIM(ADJUSTL(cval))//')'
            CALL error_lbound_arr_struc(SubzonesDredging(n)%nr_coordpol(m),&
               & 'SubzonesDredging','%nr_coordpol'//TRIM(cvar),0,.TRUE.,1,(/n/))
         ENDDO m_511
      ENDIF
      
      noerrs=COUNT(SubzonesDredging(n)%Xpol(1:MaxCoordPol,1:MaxSubzones).LT.0.0)
      IF (noerrs.GT.0) THEN
         m_512: DO m=1,MaxCoordPol
            WRITE (cval,'(I12)') m; cval = ADJUSTL(cval)
            l_5121: DO l=1,MaxSubzones
               WRITE (cval2,'(I12)') l; cval2 = ADJUSTL(cval2)
               cvar = '('//TRIM(cval)//','//TRIM(cval2)//')'
               CALL error_lbound_arr_struc(SubzonesDredging(n)%Xpol(m,l),&
                  & 'SubzonesDredging','%Xpol'//TRIM(cvar),0.0,.TRUE.,1,(/n/))
            ENDDO l_5121
         ENDDO m_512
      ENDIF
      
      noerrs=COUNT(SubzonesDredging(n)%Ypol(1:MaxCoordPol,1:MaxSubzones).LT.0.0)
      IF (noerrs.GT.0) THEN
         m_513: DO m=1,MaxCoordPol
            WRITE (cval,'(I12)') m; cval = ADJUSTL(cval)
            l_5131: DO l=1,MaxSubzones
               WRITE (cval2,'(I12)') l; cval2 = ADJUSTL(cval2)
               cvar = '('//TRIM(cval)//','//TRIM(cval2)//')'
               CALL error_lbound_arr_struc(SubzonesDredging(n)%Ypol(m,l),&
                  & 'SubzonesDredging','%Ypol'//TRIM(cvar),0.0,.TRUE.,2,(/m,l/))
            ENDDO l_5131
         ENDDO m_513
      ENDIF

   ENDDO n_510

ENDIF

!
!6. Subzones relocation attributes
!---------------------------------
!

IF (iopt_dar_time_relocate.EQ.3) THEN

   n_610: DO n=1,norelocationsites

      CALL error_lbound_arr_struc(SubzonesRelocation(n)%nr_subzones,&
                         & 'SubzonesRelocation','%nr_subzones',1,.TRUE.,1,(/n/))
      
      noerrs = COUNT(SubzonesRelocation(n)%nr_coordpol(1:MaxSubzones).LT.1)
      IF (noerrs.GT.0) THEN
         m_611: DO m=1,MaxSubzones
         WRITE (cval,'(I12)') m; cvar = '('//TRIM(ADJUSTL(cval))//')'
            CALL error_lbound_arr_struc(SubzonesRelocation(n)%nr_coordpol(m),&
                         & 'SubzonesRelocation','%nr_coordpol'//TRIM(cvar), &
                 & 0,.TRUE.,1,(/m/))
         ENDDO m_611
      ENDIF
      
      noerrs = COUNT(SubzonesRelocation(n)&
                   & %Xpol(1:MaxCoordPol,1:MaxSubzones).LT.0.0)
      IF (noerrs.GT.0) THEN
         m_612: DO m=1,MaxCoordPol
            WRITE (cval,'(I12)') m; cval = ADJUSTL(cval)
            l_6121: DO l=1,MaxSubzones
               WRITE (cval2,'(I12)') l; cval2 = ADJUSTL(cval2)
               cvar = '('//TRIM(cval)//','//TRIM(cval2)//')'
               CALL error_lbound_arr_struc(SubzonesRelocation(n)%Xpol(m,l),&
                        & 'SubzonesRelocation','%Xpol'//TRIM(cvar),0.0,.TRUE.,&
                        & 2,(/m,l/))
            ENDDO l_6121
         ENDDO m_612
      ENDIF
      
      noerrs = COUNT(SubzonesRelocation(n)&
                   & %Ypol(1:MaxCoordPol,1:MaxSubzones).LT.0.0)
      IF (noerrs.GT.0) THEN
         m_613: DO m=1,MaxCoordPol
            WRITE (cval,'(I12)') m; cval = ADJUSTL(cval)
            l_6131: DO l=1,MaxSubzones
               WRITE (cval2,'(I12)') l; cval2 = ADJUSTL(cval2)
               cvar = '('//TRIM(cval)//','//TRIM(cval2)//')'
               CALL error_lbound_arr_struc(SubzonesRelocation(n)%Ypol(m,l),&
                        & 'SubzonesRelocation','%Ypol'//TRIM(cvar),0.0,.TRUE.,&
                        & 2,(/m,l/))
            ENDDO l_6131
         ENDDO m_613
      ENDIF
      
   ENDDO n_610

ENDIF

! TODO: polygon error checking in check_dredging!!!

!!$  ! checking if polygon is closed, does not have double points or intersects
!!$   ! with itself
!!$   xpol = DredgingSites(i)%Xpol
!!$   ypol = DredgingSites(i)%Ypol
!!$   npol = SIZE(xpol)
!!$   CALL polycheck(DredgingSites(i)%Xpol,DredgingSites(i)%Ypol, SIZE(xpol), doublepoints)
!!$   CALL polyintersect(xpol,ypol,npol,intersect)
!!$   ! IF so, return an error
!!$   IF (ANY(isok)) THEN
!!$   ! error double points
!!$ELSEIF (ANY(intersect)) THEN
!!$   ! error intersect
!!$   ! If not, then determine points inside poygon

CALL error_abort('check_dar_spec',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_dar_spec

!=====================================================================

SUBROUTINE check_sedics
!************************************************************************
!
! *check_sedics* Check initial conditions for the sediment transport model
!
! Author - Boudewijn Decrop and Alexander Breugem (IMDC)
!
! Version - @(COHERENS)check_sediments.f90  V2.10
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - error_abort, error_lbound_arr, error_limits_arr
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphpars
USE morphswitches
USE sedarrays
USE sedpars
USE sedswitches
USE switches
USE error_routines, ONLY: error_abort, error_lbound_arr, error_limits_arr
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=12) :: ciglb, cjglb, cval
INTEGER :: f, i, iglb, j, jglb, k, noerrs
REAL, DIMENSION(ncloc,nrloc) :: array2dc
REAL, PARAMETER :: real_small = 1.0E-06


procname(pglev+1) = 'check_sedics'
CALL log_timer_in()

!
!1. Particle attributes
!----------------------
!

noerrs = COUNT(dp.LE.0.0)
IF (noerrs.GT.0) THEN
   f_110: DO f=1,nfp
      CALL error_lbound_arr(dp(f),'dp',0.0,.FALSE.,1,indx=(/f/))
   ENDDO f_110
ENDIF

noerrs = COUNT(rhos.LE.0.0)
IF (noerrs.GT.0) THEN
   f_120: DO f=1,nfp
      CALL error_lbound_arr(rhos(f),'rhos',0.0,.FALSE.,1,indx=(/f/))
   ENDDO f_120
ENDIF

IF (iopt_sed_bstres_cr.EQ.1) THEN
   noerrs = COUNT(bstres_cr_cst.LE.0.0)
   IF (noerrs.GT.0) THEN
      f_130: DO f=1,nfp
         CALL error_lbound_arr(bstres_cr_cst(f),'bstres_cr_cst',0.0,.FALSE.,1,&
                             & indx=(/f/))
      ENDDO f_130
   ENDIF
ENDIF

IF (iopt_sed.EQ.1.AND.iopt_sed_ws.EQ.1) THEN
   noerrs = COUNT(ws_cst.LE.0.0)
   IF (noerrs.GT.0) THEN
      f_140: DO f=1,nfp
         CALL error_lbound_arr(ws_cst(f),'ws_cst',0.0,.FALSE.,1,indx=(/f/))
      ENDDO f_140
   ENDIF
ENDIF

!
!2. Bed fractions
!----------------
!
!---bed fraction values
noerrs = COUNT((bed_fraction(1:ncloc,1:nrloc,1:nb,1:nf).LT.0.0.OR.&
              & bed_fraction(1:ncloc,1:nrloc,1:nb,1:nf).GT.1.0))
IF (noerrs.GT.0) THEN
   f_210: DO f=1,nf
   k_210: DO k=1,nb
   j_210: DO j=1,nrloc
   i_210: DO i=1,ncloc
      IF (maskatc_int(i,j)) THEN            
         iglb = nc1loc+i-1; jglb = nr1loc+j-1
         CALL error_limits_arr(bed_fraction(i,j,k,f),'bed_fraction',0.0,1.0,&
                             & 4,(/iglb,jglb,k,f/))
      ENDIF
   ENDDO i_210
   ENDDO j_210
   ENDDO k_210
   ENDDO f_210
ENDIF

!---sum of the bed fractions must be one
array2dc  = SUM(bed_fraction(1:ncloc,1:nrloc,1,:),DIM=3)
noerrs = COUNT(maskatc_int.AND.ANY(array2dc.NE.0.0).AND.&
             & ABS(array2dc-1.0).GT.real_small)
IF (noerrs.GT.0) THEN
   j_221: DO j=1,nrloc
   i_221: DO i=1,ncloc
      IF (maskatc_int(i,j).AND.array2dc(i,j).NE.0.0.AND.&
        & ABS(array2dc(i,j)-1.0).GT.real_small) THEN
         nerrs = nerrs + 1
         IF (errchk.AND.nerrs.LE.maxerrors) THEN
            iglb = nc1loc+i-1; jglb = nr1loc+j-1
            WRITE (cval,'(G12.7)') array2dc(i,j)
            WRITE (ciglb,'(I12)') iglb; ciglb = ADJUSTL(ciglb)
            WRITE (cjglb,'(I12)') jglb; cjglb = ADJUSTL(cjglb)
            WRITE (ioerr,'(A)') 'Invalid value for element ('&
                              & //TRIM(ciglb)//','//TRIM(cjglb)//',1)'//&
                              & ' of real array bed_fraction'//': '//TRIM(cval)
            WRITE (ioerr,'(A)') 'Sum over all fractions must be equal'//&
                              &' to one'
         ENDIF
      ENDIF
   ENDDO i_221
   ENDDO j_221
ENDIF

!
!3. Check volume concentrations
!------------------------------
!

noerrs = COUNT((cvol(1:ncloc,1:nrloc,1:nz,1:nf).LT.0.0.OR.&
              & cvol(1:ncloc,1:nrloc,1:nz,1:nf).GT.cmax))
IF (noerrs.GT.0) THEN
   j_310: DO j=1,nrloc
   i_310: DO i=1,ncloc
      iglb = nc1loc+i-1; jglb = nr1loc+j-1
      IF (nodeatc(i,j).EQ.1) THEN
         f_311: DO f=1,nf
         k_311: DO k=1,nz
            CALL error_limits_arr(cvol(i,j,k,f),'cvol',0.0,cmax,4,&
                               & (/iglb,jglb,k,f/))
         ENDDO k_311
         ENDDO f_311
      ENDIF
   ENDDO i_310
   ENDDO j_310
ENDIF

CALL error_abort('check_sedics',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_sedics

!=====================================================================

SUBROUTINE check_morphics
!************************************************************************
!
! *check_morphics* Check initial conditions for the morphological model
!
! Author - Boudewijn Decrop and Alexander Breugem (IMDC)
!
! Version - @(COHERENS)check_sediments.f90  V2.8
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - error_abort, error_lbound_arr, error_limits_arr
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphswitches
USE sedarrays
USE sedpars
USE error_routines, ONLY: error_abort, error_lbound_arr, error_limits_arr
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=12) :: ciglb, cjglb, ckglb, cval
INTEGER :: i, iglb, j, jglb, k, noerrs
REAL, PARAMETER :: real_small = 1.0E-06
REAL, DIMENSION(ncloc,nrloc,nb) :: array3dc


procname(pglev+1) = 'check_morphics'
CALL log_timer_in()

!---bed_porosity must be positive and lower than 1
noerrs = COUNT(bed_porosity.LT.0.0.OR.bed_porosity.GT.1.0)
IF (noerrs.GT.0) THEN
   j_110: DO j=1,nrloc
   i_110: DO i=1,ncloc
      IF (seapoint(i,j)) THEN
         k_111: DO k=1,nb
            iglb = nc1loc+i-1; jglb = nr1loc+j-1
            CALL error_limits_arr(bed_porosity(i,j,k),'bed_porosity', &
                                & 0.0,1.0,3,(/iglb,jglb/))
         ENDDO k_111
      ENDIF
   ENDDO i_110
   ENDDO j_110
ENDIF

!---bed layer thickness must be positive
IF (iopt_morph_fixed_layer.EQ.1) THEN
   noerrs = COUNT(bed_layer_thickness.LT.0.0)
   IF (noerrs.GT.0) THEN
      j_210: DO j=1,nrloc
      i_210: DO i=1,ncloc
         IF (seapoint(i,j)) THEN
            k_211: DO k=1,nb
               iglb = nc1loc+i-1; jglb = nr1loc+j-1
               CALL error_lbound_arr(bed_layer_thickness(i,j,k),&
                                  & 'bed_layer_thickness',0.0,.TRUE.,3,&
                                   & indx=(/iglb,jglb,k/))
            ENDDO k_211
         ENDIF
      ENDDO i_210
      ENDDO j_210
   ENDIF
ENDIF

!---sum of the bed fractions must be one
IF (nf.GT.1) THEN
   array3dc  = SUM(bed_fraction(1:ncloc,1:nrloc,:,:),DIM=4)
   k_321: DO k=1,nb
   j_321: DO j=1,nrloc
   i_321: DO i=1,ncloc
      IF (maskatc_int(i,j).AND.bed_layer_thickness(i,j,k).GT.0.0.AND.&
        & ABS(array3dc(i,j,k)-1.0).GT.real_small) THEN
         nerrs = nerrs + 1
         IF (errchk.AND.nerrs.LE.maxerrors) THEN
            iglb = nc1loc+i-1; jglb = nr1loc+j-1
            WRITE (cval,'(G12.7)') array3dc(i,j,k)
            WRITE (ciglb,'(I12)') iglb; ciglb = ADJUSTL(ciglb)
            WRITE (cjglb,'(I12)') jglb; cjglb = ADJUSTL(cjglb)
            WRITE (ckglb,'(I12)') k; ckglb = ADJUSTL(ckglb)
            WRITE (ioerr,'(A)') 'Invalid value for element ('&
                              & //TRIM(ciglb)//','//TRIM(cjglb)//','&
                              & //TRIM(ckglb)//') of real array bed_fraction'//&
                              & ': '//TRIM(cval)
            WRITE (ioerr,'(A)') 'Sum over all fractions must be equal'//&
                                  &' to one'
         ENDIF
      ENDIF
   ENDDO i_321
   ENDDO j_321
   ENDDO k_321
ENDIF

CALL error_abort('check_morphics',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_morphics

END MODULE check_sediments
