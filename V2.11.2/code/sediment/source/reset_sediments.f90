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

MODULE reset_sediments
!************************************************************************
!
! *reset_sediments* Reset sediment model parameters if needed
!
! Author - Boudewijn Decrop, Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)reset_sediments.f90  V2.11.1
!
! $Date: 2018-05-23 10:09:17 +0200 (Wed, 23 May 2018) $
!
! $Revision: 1138 $
!
! Description - 
!
! Reference -
!
! Routines - reset_dar_params, reset_morph_params, reset_sed_params,
!            reset_sedics
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS

!=====================================================================

SUBROUTINE reset_sed_params
!************************************************************************
!
! *reset_sed_params* Reset parameters for the sediment transport model
!
! Author - Boudewijn Decrop and Alexander Breugem (IMDC)
!
! Version - @(COHERENS)reset_sediments.f90  V2.11.1
!
! Description -
!
! Calling program - initialise_model
!
! Module calls - warning_reset_var
!
!************************************************************************
!
USE gridpars
USE iopars
USE morphpars
USE sedpars
USE sedswitches
USE switches
USE timepars
USE error_routines, ONLY: warning_reset_var
USE inout_routines, ONLY: open_file
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'reset_sed_params'
CALL log_timer_in()

!
!1. Reset switches if necessary
!------------------------------
!
!---flocculation module
IF (iopt_sed.EQ.2) THEN
   CALL warning_reset_var(iopt_dar,'iopt_dar',0)
   CALL warning_reset_var(iopt_sed_bbc,'iopt_sed_bbc',1)
   CALL warning_reset_var(iopt_sed_bedeq,'iopt_sed_bedeq',0)
   CALL warning_reset_var(iopt_sed_hiding,'iopt_sed_hiding',0)
   CALL warning_reset_var(iopt_sed_ws_hindset,'iopt_sed_ws_hindset',0)
   CALL warning_reset_var(iopt_sed_median,'iopt_sed_median',1)
   CALL warning_reset_var(iopt_sed_mode,'iopt_sed_mode',2)
   CALL warning_reset_var(iopt_sed_nodim,'iopt_sed_nodim',3)
   CALL warning_reset_var(iopt_sed_slope,'iopt_sed_slope',0)
   CALL warning_reset_var(iopt_sed_toteq,'iopt_sed_toteq',0)
   CALL warning_reset_var(iopt_sed_type,'iopt_sed_type',2)
   CALL warning_reset_var(iopt_sed_wave_diff,'iopt_sed_wave_diff',0)
   CALL warning_reset_var(iopt_sed_ws_floc,'iopt_sed_ws_floc',0)
   CALL warning_reset_var(iopt_tidal_accel,'iopt_tidal_accel',0)
ENDIF

!---sand
IF (iopt_sed_type.EQ.1) THEN
   IF (iopt_sed_bbc.GT.0.AND.iopt_sed_bbc.NE.2) THEN
      CALL warning_reset_var(iopt_sed_bbc,'iopt_sed_bbc',2)   
   ENDIF
   CALL warning_reset_var(iopt_sed_ws_floc,'iopt_sed_ws_floc',0)
   IF (iopt_sed_bedeq.EQ.4) THEN
      CALL warning_reset_var(iopt_sed_bstres_cr,'iopt_sed_bstres_cr',4)
      IF (nf.GT.1) THEN
         CALL warning_reset_var(iopt_sed_hiding,'iopt_sed_hiding',1)
      ENDIF
      CALL warning_reset_var(iopt_sed_ws,'iopt_sed_ws',6)
   ENDIF

!--- mud
ELSEIF (iopt_sed_type.EQ.2) THEN
   IF (iopt_sed_bbc.GT.0.AND.iopt_sed_bbc.NE.1) THEN
      CALL warning_reset_var(iopt_sed_bbc,'iopt_sed_bbc',1)   
   ENDIF
   IF (iopt_sed_beta.EQ.3) THEN
      CALL warning_reset_var(iopt_sed_beta,'iopt_sed_beta',1)
   ENDIF
   IF (iopt_sed_ws_hindset.GT.0) THEN
      CALL warning_reset_var(iopt_sed_ws_hindset,'iopt_sed_ws_hindset',2)
   ENDIF
   IF (iopt_sed_ws.GT.1.AND.iopt_sed_ws.NE.4 ) THEN
      CALL warning_reset_var(iopt_sed_ws,'iopt_sed_ws',3)
   ENDIF
   CALL warning_reset_var(iopt_sed_hiding,'iopt_sed_hiding',0)
   CALL warning_reset_var(iopt_sed_slope,'iopt_sed_slope',0)
ENDIF

!---advection scheme for settling
IF (iopt_adv_scal.GT.0.AND.iopt_sed_vadv.GT.0) THEN
   CALL warning_reset_var(iopt_sed_vadv,'iopt_sed_vadv',iopt_adv_scal)
ENDIF

!---2-D application
IF (iopt_grid_nodim.EQ.2) THEN
   CALL warning_reset_var(iopt_sed_wave_diff,'iopt_sed_wave_diff',0)
   CALL warning_reset_var(iopt_sed_nodim,'iopt_sed_nodim',2)
   CALL warning_reset_var(iopt_sed_vadv,'iopt_sed_vadv',0)
ENDIF
IF (iopt_sed_nodim.EQ.2) THEN
   CALL warning_reset_var(iopt_sed_ws_hindset,'iopt_sed_ws_hindset',0)
ENDIF

!---density
IF (iopt_dens.EQ.0.AND.iopt_sed_dens_grad.EQ.0) THEN
   CALL warning_reset_var(iopt_dens_grad,'iopt_dens_grad',0)
ENDIF

!---tidal acceleration
IF (iopt_tidal_accel.EQ.1) THEN
   CALL warning_reset_var(iopt_sed_nodim,'iopt_sed_nodim',2)
ENDIF

!---compatability with Wu et al. (2000)
IF (iopt_sed_bbc_eq.EQ.6.OR.iopt_sed_bedeq.EQ.4) THEN
   CALL warning_reset_var(iopt_sed_bstres_cr,'iopt_sed_bstres_cr',4)
   CALL warning_reset_var(iopt_sed_ws,'iopt_sed_ws',6)
   IF (nf.GT.1) THEN
      CALL warning_reset_var(iopt_sed_hiding,'iopt_sed_hiding',1)
   ENDIF
ENDIF

!---compatability with Van Rijn (2003/2007) 
IF (iopt_sed_bbc_eq.EQ.5.OR.iopt_sed_bedeq.GE.6.OR.iopt_sed_toteq.GT.5) &
     & THEN
   CALL warning_reset_var(iopt_sed_beta,'iopt_sed_beta',3)
   CALL warning_reset_var(iopt_sed_bbc_ref,'iopt_sed_bbc_ref',2)
   IF (nf.EQ.1) CALL warning_reset_var(iopt_sed_hiding,'iopt_sed_hiding',1)
ENDIF

!---when Wu total load is selected, the bed load formula should be Wu (2000)
IF (iopt_sed_toteq.EQ.5) THEN
   CALL warning_reset_var(iopt_sed_bedeq,'iopt_sed_bedeq',4)
   CALL warning_reset_var(iopt_sed_bstres_cr,'iopt_sed_bstres_cr',4)
   IF (nf.GT.1) THEN
      CALL warning_reset_var(iopt_sed_hiding,'iopt_sed_hiding',1)
   ENDIF
   CALL warning_reset_var(iopt_sed_ws,'iopt_sed_ws',6)
ENDIF

!---van Rijn (2003) total load
IF (iopt_sed_toteq.EQ.6) THEN
   CALL warning_reset_var(iopt_sed_bedeq,'iopt_sed_bedeq',6)
   CALL warning_reset_var(iopt_sed_bbc_ref,'iopt_sed_bbc_ref',2)
ENDIF

!---when van Rijn (2007) total load is selected, the bed load formula should be
!   van Rijn (2007)
IF (iopt_sed_toteq.EQ.7) THEN
   CALL warning_reset_var(iopt_sed_bedeq,'iopt_sed_bedeq',7)
ENDIF

!---no need for gaussian quadrature in current-wave sediment transport if
!   no waves
IF (iopt_waves.EQ.0) THEN
   CALL warning_reset_var(nrquad_wav,'nrquad_wav',1)
   CALL warning_reset_var(iopt_sed_wave_diff,'iopt_sed_wave_diff',0)
ENDIF

!---no flocculation with Camenen (Mud)
IF (iopt_sed_ws.EQ.3) THEN
    CALL warning_reset_var(iopt_sed_ws_floc,'iopt_sed_ws_floc',0)
ENDIF

!---bottom stress formulation
IF (iopt_bstres_form.NE.2) THEN
   CALL warning_reset_var(iopt_sed_bstres,'iopt_sed_bstres',0)
ENDIF
   
!
!2. Parameters
!--------------
!
!---fractions
IF (iopt_sed.EQ.2)  CALL warning_reset_var(nf,'nf',3)
nfp = MERGE (nf,1,iopt_sed.EQ.1)

!---number of bed layers
IF (iopt_morph.EQ.0) CALL warning_reset_var(nb,'nb',1)

!---time step
deltsed = icsed*delt2d

!---number of variables used at open boundaries and for nesting
maxsedvars = nf

CALL log_timer_out()


RETURN

END SUBROUTINE reset_sed_params

!=====================================================================

SUBROUTINE reset_morph_params
!************************************************************************
!
! *reset_morph_params* Reset parameters for the morphological model
!
! Author - Boudewijn Decrop and Alexander Breugem (IMDC)
!
! Version - @(COHERENS)reset_sediments.f90  V2.11.1
!
! Description -
!
! Calling program - initialise_model
!
! Module calls - warning_reset_var
!
!************************************************************************
!
USE iopars
USE morphpars
USE morphswitches
USE sedpars
USE sedswitches
USE switches
USE timepars
USE error_routines, ONLY: warning_reset_var
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'reset_morph_params'
CALL log_timer_in()

!
!1. Reset switches if necessary
!------------------------------
!
!---Ribberink exchange model only works for two fractions
IF (nf.NE.2) THEN
   CALL warning_reset_var(iopt_morph_vert_fluxes,&
                       & 'iopt_morph_vert_fluxes',1)
ENDIF

!---reset switches if morphology is off
IF (iopt_morph.EQ.0) THEN
   CALL warning_reset_var(iopt_morph_active_layer,'iopt_morph_active_layer',0)
   CALL warning_reset_var(iopt_morph_avalanching,'iopt_morph_avalanching',0)
   CALL warning_reset_var(iopt_morph_fixed_layer,'iopt_morph_fixed_layer',0)
   CALL warning_reset_var(iopt_morph_vert_fluxes,'iopt_morph_vert_fluxes',0)
   CALL warning_reset_var(iopt_morph_fixed_layer,'iopt_morph_fixed_layer',0)
ENDIF

!---fixed layer in case of multi layers
IF (nb.GT.1) THEN
   CALL warning_reset_var(iopt_morph_fixed_layer,'iopt_morph_fixed_layer',1)
ENDIF

!
!2. Parameters
!--------------
!
!---warning if vertical sorting and morphological factor are combined 
IF (nb.GT.1.AND.morph_factor.GT.1.AND.warnflag) THEN
   WRITE (iowarn,'(A)') 'WARNING: combining morphological factor with '//&
                      & 'vertical sorting.'
   WRITE (iowarn,'(A)') 'Check result files for instabilities or '//&
                        'unexpected results'
ENDIF

!---morphological acceleration factor
IF (iopt_tidal_accel.EQ.1) THEN
   CALL warning_reset_var(morph_factor,'morph_factor',1.0)
ENDIF

IF (icmorph.LT.icsed) icmorph = icsed
IF (icaval.LT.icmorph) icaval = icmorph

!---time step
deltmorph = icmorph*delt2d

!---counters
accelstep = 0

!---sediment error
sediment_err = 0.0

CALL log_timer_out()


RETURN

END SUBROUTINE reset_morph_params

!========================================================================

SUBROUTINE reset_dar_params
!************************************************************************
!
! *reset_dar_params* Reset parameters for the dredging/relocation model
!
! Author - Boudewijn Decrop and Alexander Breugem (IMDC)
!
! Version - @(COHERENS)reset_sediments.f90  V2.8
!
! Description -
!
! Calling program - initialise_model
!
! Module calls - warning_reset_var
!
!************************************************************************
!
USE darswitches
USE iopars
USE sedswitches
USE switches
USE error_routines, ONLY: warning_reset_var
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'reset_dar_params'
CALL log_timer_in()


!---temporal spreading
IF (iopt_dar_duration.EQ.1.OR.iopt_dar_duration.EQ.2) THEN
   CALL warning_reset_var(iopt_dar_effect,'iopt_dar_effect',1)
   IF (iopt_dar_time_dredge.EQ.3) THEN
      CALL warning_reset_var(iopt_dar_time_dredge,'iopt_dar_time_dredge',1)
   ENDIF
ENDIF
IF (iopt_dar_duration.EQ.3) THEN
   CALL warning_reset_var(iopt_dar_coupling_time,'iopt_dar_coupling_time',3)
ENDIF
IF (iopt_dar_coupling_time.EQ.3) THEN
   CALL warning_reset_var(iopt_dar_duration,'iopt_dar_duration',3)
ENDIF

!---spatial spreading
IF (iopt_dar_time_dredge.EQ.2.OR.iopt_dar_time_dredge.EQ.3) THEN
   CALL warning_reset_var(iopt_dar_duration,'iopt_dar_duration',3)
   CALL warning_reset_var(iopt_dar_coupling_time,'iopt_dar_coupling_time',3)
ENDIF
IF (iopt_dar_time_dredge.EQ.3) THEN
   CALL warning_reset_var(iopt_dar_effect,'iopt_dar_effect',2)
ENDIF

!---plume modeling 
IF (iopt_dar_effect.EQ.2) THEN
   CALL warning_reset_var(iopt_dar_time_dredge,'iopt_dar_time_dredge',3)
   CALL warning_reset_var(iopt_dar_duration,'iopt_dar_duration',3)
   CALL warning_reset_var(iopt_dar_coupling_time,'iopt_dar_coupling_time',3)
   IF (iopt_dar_relocation_criterion.EQ.1) THEN
      CALL warning_reset_var(iopt_dar_relocation_criterion,&
                          & 'iopt_dar_relocation_criterion',2)
   ENDIF
ENDIF


CALL log_timer_out()


RETURN

END SUBROUTINE reset_dar_params

!=====================================================================

SUBROUTINE reset_sedics
!************************************************************************
!
! *reset_sedics* Reset initial conditions for the sediment transport model
!
! Author - Boudewijn Decrop and Alexander Breugem (IMDC)
!
! Version - @(COHERENS)reset_sediments.f90  V2.9
!
! Description - 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE reset_sedics


END MODULE reset_sediments
