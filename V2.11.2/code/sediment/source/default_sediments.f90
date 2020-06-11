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

MODULE default_sediments
!************************************************************************
!
! *default_sediments* Default settings for sediment transport, morphological and
!                     dredging/relocation models
!
! Author - Boudewijn Decrop, Alexander Breugem (IMDC)
!
! Version - @(COHERENS)default_sediments.f90  V2.11.1
!
! $Date: 2018-05-23 10:09:17 +0200 (Wed, 23 May 2018) $
!
! $Revision: 1138 $
!
! Description - 
!
! Reference -
!
! Routines - default_dar_params, default_morphics, default_morph_params,
!            default_sedics, default_sed_params, default_sed_spec
!
!************************************************************************
!


CONTAINS

!====================================================================

SUBROUTINE default_sed_params
!************************************************************************
!
! *default_sed_params* Default settings of parameters for the sediment
!                      transport model
!
! Author - Boudewijn Decrop, Alexander Breugem (IMDC)
!
! Version - @(COHERENS)default_sediments.f90  V2.11.1
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE morphpars
USE physpars
USE sedpars
USE sedswitches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'default_sed_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!
!---type of transport
iopt_sed_bedeq = 1; iopt_sed_mode = 2; iopt_sed_nodim = 3
iopt_sed_toteq = 1; iopt_sed_type = 1

!---bottom boundary condition
iopt_sed_bbc = 1; iopt_sed_bbc_eq = 1; iopt_sed_bbc_ref = 1

!---sediment density
iopt_sed_dens_grad = 0

!---sediment diffusion
iopt_sed_beta = 1; iopt_sed_wave_diff = 0

!---bed stress
iopt_sed_bstres = 0; iopt_sed_bstres_cr = 1; iopt_sed_hiding = 0
iopt_sed_rough = 1

!---settling velocity
iopt_sed_ws_hindset = 0; iopt_sed_vadv = 3; iopt_sed_ws = 1
iopt_sed_ws_floc = 0; iopt_sed_ws_lim = 0

!---slope effects
iopt_sed_slope = 0

!---open boundary conditions
iopt_sed_obc_flux = 0

!---numerical
iopt_sed_filter = 0; iopt_sed_median = 1; iopt_sed_source_impl = 3

!
!2. Parameters
!-------------
!
!---number of fractions
nf = 1

!---number of bed layers
nb = 1

!---time counter
icsed = ic3d

!---Gauss-Lengendre integration
nrquad_sed = 7; nrquad_wav = 10; theta_sed_src = 0.501

!---diffusion coefficient
beta_sed_cst = 1.0

!---hiding factor
wu_exp = -0.6

!---bottom fluxes
cgel = 0.0; cmax = 0.65; coef_bed_grad = 1.3; height_c_cst = 0.01
n_RichZaki = 4.6; parth_coef = 2.0E-04; parth_exp = 1.0
zrough_grain = 0.05; zrough_sed_cst = 0.0; z0_coef = 30.0

!---fall velocity
alpha_VR = 2.19; a_leussen = 0.02; b_leussen = 0.0024; 
floc_VR_max = 10.0; floc_VR_min = 1.0 
wslimfac = 0.5

!---numerical
maxitbartnicki = 100

!---flocculation module
agg_alpha = 0.1; brk_es = 1.0E-04; brk_ffy = 1.0E-10; brk_frac = 1.0
brk_p = 1.0; brk_q = 0.7; floc_ncmax = 100.0; floc_ncmin = 2.0; nfrdim = 2.0

CALL log_timer_out()


RETURN

END SUBROUTINE default_sed_params

!====================================================================

SUBROUTINE default_morph_params
!************************************************************************
!
! *default_morph_params* Default settings of parameters for the morphological
!                        model
!
! Author - Boudewijn Decrop, Alexander Breugem (IMDC)
!
! Version - @(COHERENS)default_sediments.f90  V2.10
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE morphpars
USE morphswitches
USE sedpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'default_morph_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!
!---time integration
iopt_morph_time_int = 1
!---horizontal diffusion
iopt_morph_hdiff = 0
!---active layer
iopt_morph_active_layer = 0
!---vertical exchange in sea bed
iopt_morph_vert_fluxes = 0; 
!---fixed layers
iopt_morph_fixed_layer = 0
!---avalanching
iopt_morph_avalanching = 0
!---tidal averaging
iopt_morph_tidal_scheme = 1
!---mass correction
iopt_morph_corr = 0

!
!2. Parameters
!-------------
!
!---numerical parameters
k1_harris = 0.7; k2_harris = 6.0
alpha_jameson = 2.0; k2_jameson = 8.0; k4_jameson = 1.0/32.0
max_iterations_aval = 500

!---bed update
bed_porosity_cst = 0.4

!---vertical sorting parameters
average_bedform_length = 0.5
average_bedform_height = 0.05
dune_celerity = 0.001
similarity_range = 0.05
trough_probability = 0.1

!---bank erosion parameters
dyn_angle_dry = 25.0; dyn_angle_wet = 30.0
stat_angle_dry = 34.0; stat_angle_wet = 45.0

!---morphological acceleration         
accelstep = 0
morph_factor  = 1.0
morph_steps = 0
nstep_hydro = 0
number_tidal_steps = 20 

!---time counters
icmorph = icsed
icaval = icsed
icsedbal = 0

CALL log_timer_out()


RETURN

END SUBROUTINE default_morph_params

!========================================================================

SUBROUTINE default_dar_params
!************************************************************************
!
! *default_dar_params* Default settings of parameters for dredging and
!                      relocation
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)default_sediments.f90  V2.8
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE darpars
USE darswitches
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'default_dar_params'
CALL log_timer_in()

!
!1. Switches
!-----------
!
!--- type of modeling (morphology/plume)
iopt_dar_effect = 1; iopt_dar_scenario = 1

!--- coupling (time/space)
iopt_dar_coupling_space = 1; iopt_dar_coupling_time = 1
iopt_dar_relocation_criterion = 1; iopt_dar_user_volume = 1

!--- time spreading
iopt_dar_time_dredge = 1; iopt_dar_time_relocate = 1
iopt_dar_duration = 1

!--- dredging criterion
iopt_dar_dredging_criterion = 1

!
!2. Parameters
!-------------
!
!---number of sites/campaigns
nodredgingsites = 1
norelocationsites = 1
nocampaigns = 1
!--deepest point relocation
nocirclepoints = 16
!--track dredging
maxtrackpoints = 1000

CALL log_timer_out()


RETURN

END SUBROUTINE default_dar_params

!========================================================================

SUBROUTINE default_sedics
!************************************************************************
!
! *default_sedics* Default settings of initial conditions for the sediment
!                  transport model
!
! Author - Boudewijn Decrop, Alexander Breugem (IMDC)
!
! Version - @(COHERENS)default_sediments.f90  V2.8
!
! Description - 
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
USE sedarrays
USE sedpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: f, k


procname(pglev+1) = 'default_sedics'
CALL log_timer_in()

!---bed fractions
f_110: DO f=1,nf
k_110: DO k=1,nb
   WHERE (nodeatc(0:ncloc,0:nrloc).EQ.1)
      bed_fraction(:,:,k,f) = 1.0/REAL(nf)
   END WHERE
ENDDO k_110
ENDDO f_110

CALL log_timer_out()


RETURN

END SUBROUTINE default_sedics

!========================================================================

SUBROUTINE default_morphics
!************************************************************************
!
! *default_morphics* Default settings of initial conditions for the
!                    morphological model
!
! Author - Boudewijn Decrop, Alexander Breugem (IMDC)
!
! Version - @(COHERENS)default_sediments.f90  V2.8
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE morpharrays
USE morphpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'default_morphics'
CALL log_timer_in()

!---bed porosity
bed_porosity = bed_porosity_cst

CALL log_timer_out()


RETURN

END SUBROUTINE default_morphics

!========================================================================

SUBROUTINE default_sed_spec
!************************************************************************
!
! *default_sed_spec* Default settings of sediment particle attributes
!
! Author - Boudewijn Decrop, Alexander Breugem (IMDC)
!
! Version - @(COHERENS)default_sediments.f90  V2.9
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE sedarrays
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'default_sed_spec'
CALL log_timer_in()

!---particle attributes
dp = 0.0E-6; rhos = 2650.0
bstres_cr_cst = 0.0; ws_cst = 0.0

CALL log_timer_out()


RETURN

END SUBROUTINE default_sed_spec

END MODULE default_sediments
