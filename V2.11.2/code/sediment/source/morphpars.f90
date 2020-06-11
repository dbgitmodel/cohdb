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

MODULE morphpars
!************************************************************************
!
! *morphpars* Morphological parameters
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)morphpars.f90  V2.10
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description - 
!
!************************************************************************
!
!
USE syspars
 
IMPLICIT NONE

LOGICAL :: morphstep
INTEGER :: accelstep = 0, icaval = 1, icmorph = 1, icsedbal = 0, &
         & max_iterations_aval = 500, morph_counter, morph_steps = 0, &
         & nstep_hydro, number_tidal_steps = 20, period_end, period_start, &
         & tidal_counter
REAL :: alpha_jameson = 2.0, average_bedform_height = 0.5, &
      & average_bedform_length = 0.05, bed_porosity_cst = 0.4, deltmorph, &
      & dune_celerity = 0.001, dyn_angle_dry = 25.0, dyn_angle_wet = 30.0, &
      & k1_harris = 0.7, k2_harris = 6.0, k2_jameson = 8.0, &
      & k4_jameson = 1.0/32.0, morph_factor = 1.0, sediment_balance, &
      & sediment_err, sediment_obc, sediment_vol, similarity_range = 0.05, &
      & stat_angle_dry = 34.0, stat_angle_wet = 45.0, trough_probability = 0.1

SAVE

!
! Name                    Type    Purpose
!------------------------------------------------------------------------------
!*accelstep*              INTEGER Number of time steps for storing or
!                                 retrieving hydrodynamical data (used in case
!                                 of tidal acceleration)
!*alpha_jameson*          REAL    "alpha" factor in the Jameson et al. (1981)
!                                 parameterisation for the artificial diffusion
!                                 coefficient (4th order term) in the Exner
!                                 equation artifial diffusion (4th order)
!*average_bedform_height* REAL    Used in determination active layer thickness
!                                 in the presence of bed forms
!*average_bedform_length* REAL    Used in the calculation of the Riberrink
!                                 exchange flux
!*bed_porosity_cst*       REAL    User-defined spatially uniform value for the
!                                 bed porosity
!*deltmorph*              REAL    Time step for morphological update         [s]
!*dune_celerity*          REAL    The celerity of the average dune used in the
!                                 two-layer sediment sorting model
!*dyn_angle_dry*          REAL    Angle of repose after avalanching in dry cells
!*dyn_angle_wet*          REAL    Angle of repose after avalanching in wet cells
!*icmorph*                INTEGER Number of 3-D time steps per morphological
!                                 time step
!*icsedbal*               INTEGER Number 3-D time steps defining the
!                                 integration period for checking the global
!                                 sediment balance
!*k1_harris*              REAL    Factor in Harris & Wiberg formula for active
!                                 layer thickness
!*k2_harris*              REAL    Factor in Harris & Wiberg formula for active
!                                 layer thickness
!*k2_jameson*             REAL    "k2" factor in the Jameson et al. (1981)
!                                 parameterisation for the artificial diffusion
!                                 coefficient (4th order term) in the Exner
!                                 equation artifial diffusion (2nd order)
!*k4_jameson*             REAL    "k4" factor in the Jameson et al. (1981)
!                                 parameterisation for the artificial diffusion
!                                 coefficient (4th order term) in the Exner
!                                 equation artifial diffusion (4th order)
!*max_iterations_aval*    INTEGER Maximum number of iterations in the
!                                 calculation of avalanching
!*morphstep*              LOGICAL .TRUE. at morphological time step
!*morph_counter*          INTEGER Counter index of the current morphological
!                                 cycle
!*morph_factor*           REAL    Factor used for morphological acceleration
!                                 and resets the tidal_counter
!*morph_steps*            INTEGER Number of cycles in a morphological period
!                                 (used in case of tidal acceleration)
!*nr_periods*             INTEGER Number of morphological cycles (used in case
!                                 of tidal acceleration) 
!*nstep_hydro*            INTEGER Number of time steps in a hydrodynamic cycle
!                                 (used in case of tidal acceleration) 
!*number_tidal_steps*     INTEGER Number of steps at which hydrodynamic data
!                                 are stored in one characteristic tidal
!                                 averaging period
!*period*                 INTEGER Counter for morphological period (used in
!                                 case of tidal acceleration) 
!*period_end*             INTEGER End index of a characteristic period (used in
!                                 case of tidal acceleration) 
!*period_start*           INTEGER Start index of a characteristic period (used
!                                 in case of tidal acceleration) 
!*sediment_balance*       REAL    Domain integrated sediment volume error in
!                                 the sediment balance during a reference period
!                                                                        [m3/m3]
!*sediment_err*           REAL    Accumulated amount of non-allocated sediment
!                                                                           [m3]
!*sediment_obc*           REAL    Net amount of sediment volume
!                                 entering/leaving the domain at open
!                                 boundaries during a reference period   [m3/m3]
!*sediment_vol*           REAL    Net amount of sediment volume
!                                 entering/leaving the domain by
!                                 erosion/deposition during a reference period
!                                                                        [m3/m3]
!*similarity_range*       REAL    Expresses the range within which sediment
!                                 fractions are considered similar
!*stat_angle_dry*         REAL    Maximum angle of repose before avalanching
!                                 starts in dry cells
!*stat_angle_wet*         REAL    Maximum angle of repose before avalanching
!                                 starts in wet cells
!                                 characteristic period
!*tidal_counter*          INTEGER Tidal averaging counter that links the actual
!                                 time step to the apropriate step in the
!                                 characteristic period

!*trough_probability*     REAL    The probability of dune troughs reaching
!                                 below a certain threshold bed level
!
!************************************************************************
!


END MODULE morphpars
