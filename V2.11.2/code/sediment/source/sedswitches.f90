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

MODULE sedswitches
!************************************************************************
!
! *sediment* Sediment arrays
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)sediment.f90  V2.11.1
!
! $Date: 2018-05-23 10:09:17 +0200 (Wed, 23 May 2018) $
!
! $Revision: 1138 $
!
! Description - 
!
!************************************************************************
!
!
 
IMPLICIT NONE

!
!1. Switches
!-----------------------------------

!---type of transport
INTEGER :: iopt_sed_bedeq = 1, iopt_sed_mode = 2, iopt_sed_nodim = 3, &
         & iopt_sed_toteq = 1, iopt_sed_type = 1

!---bottom boundary condition
INTEGER :: iopt_sed_bbc = 1, iopt_sed_bbc_eq = 1, iopt_sed_bbc_ref = 1

!---sediment density
INTEGER :: iopt_sed_dens_grad = 0

!---sediment diffusion
INTEGER :: iopt_sed_beta = 1, iopt_sed_wave_diff = 0

!---bed stress
INTEGER :: iopt_sed_bstres = 0, iopt_sed_bstres_cr = 1, iopt_sed_hiding = 0, &
         & iopt_sed_rough = 0

!---settling velocity
INTEGER :: iopt_sed_ws_floc = 0, iopt_sed_vadv = 3, iopt_sed_ws_hindset = 0, &
         & iopt_sed_ws = 1, iopt_sed_ws_lim = 0

!---slope effects
INTEGER :: iopt_sed_slope = 0

!---open boundary conditions
INTEGER :: iopt_sed_obc_flux = 0

!---numerical
INTEGER ::  iopt_sed_filter = 0, iopt_sed_median = 1, iopt_sed_source_impl = 1

SAVE

! Name               Type    Purpose
!------------------------------------------------------------------------------
!*iopt_sed_bbc*       INTEGER Formulation for bed erosion
!                     = 0 => erosion disabled
!                     = 1 => "Standard" formulation for deposition,
!                            Partheniades-Krone for erosion
!                     = 2 => using equilibrium concentration/time in 2-D,
!                            reference concentration in 3-D
!*iopt_sed_bbc_eq*   INTEGER Selects model for equilibrium concentration
!                     = 1 => Rouse profile
!                     = 2 => q_sus/U: Engelund-Hansen (1967) standard
!                            formulation
!                     = 3 => q_sus/U: Engelund-Hansen (1967) using
!                            Chollet-Cunge (1979) formulation for the
!                            Shield parameter
!                     = 4 => q_sus/U: Ackers White (1973)
!                     = 5 => q_sus/U: Van Rijn (2003)
!                     = 6 => q_sus/U: Wu et al. (2000)
!*iopt_sed_bbc_ref*   INTEGER Selects model for reference concentation
!                     = 1 => Smith & McLean (1977)
!                     = 2 => Van Rijn (1984)
!*iopt_sed_bedeq*     INTEGER Selects formulation for bed load
!                     = 0 => disabled
!                     = 1 => Meyer-Peter Muller (1948)
!                     = 2 => Engelund & Fredsoe (1976)
!                     = 3 => Van Rijn (1984)
!                     = 4 => Wu et al (2000)
!                     = 5 => Souslby (1997)
!                     = 6 => Van Rijn (2003)
!                     = 7 => Van Rijn (2007)
!*iopt_sed_beta*      INTEGER Type of equation used for beta
!                     = 1 => beta = 1
!                     = 2 => user defined value
!                     = 3 => Van Rijn (1984)
!*iopt_sed_bstres*    INTEGER Selects type of formulation for the (skin) bottom
!                             shear stress in the sediment module
!                     = 0 => Taken from the physical model
!                     = 1 => Using a skin roughness model selected by
!                            iopt_sed_rough
!*iopt_sed_bstres_cr* INTEGER Selects critical bed shear stress
!                     = 1 => User defined uniform value for each fraction
!                     = 2 => Brownlie (1981)
!                     = 3 => Soulsby-Whitehouse (1997)
!                     = 4 => Wu et al. (2000) (NOT RECOMMENDED)
!*iopt_sed_dens_grad* INTEGER Selects density stratification model
!                     = 0 => No density stratification
!                     = 1 => Include density stratification by sediments
!                            (horizontal and vertical)
!*iopt_sed_filter*    INTEGER  Type of filter used to prevent negative
!                              concentrations
!                     = 0 => No filter
!                     = 1 => Bartnicki (1989) filter
!*iopt_sed_hiding*    INTEGER Selects hiding model used
!                     = 0 => No hiding factor
!                     = 1 => Wu et. al (2000)
!                     = 2 => Ashida & Michiue (1972)
!*iopt_sed_median*    INTEGER Selects type of calculation of d50
!                     = 1 => No interpolation
!                     = 2 => Linear interpolation (not recommended)
!*iopt_sed_mode*      INTEGER  Sediment transport mode
!                     = 1 => Bed load only (-> define iopt_sed_bedeq!)
!                     = 2 => Suspended load only (advection-diffusion)
!                     = 3 => Bed and suspended (via advection-diffusion) load
!                     = 4 => Total load formula (-> define iopt_sed_toteq!)
!*iopt_sed_nodim*     INTEGER Number of dimensions used for sediment transport
!                     = 2 => 2-D (depth averaged)
!                     = 3 => 3-D
!*iopt_sed_obc_flux*  INTEGER Disables/enables upwind correction for
!                             sediment loads at open boundaries (0/1)
!*iopt_sed_rough*     INTEGER Selects formulation for the skin bottom roughness
!                     = 1 => User-defined spatially uniform value 
!                     = 2 => User_defined spatially non-uniform value
!                     = 3 => As function of median grain-size
!*iopt_sed_slope*     INTEGER Selects bed slope effects model
!                     = 0 => No bed slope effects
!                     = 1 => Koch & Flokstra (1981)
!*iopt_sed_source_impl*INTEGER Selects discretization of the deposition 
!                              term in 2D
!                     = 0 => Explicit scheme
!                     = 1 => Semi-implict (using theta_sed_src)
!                     = 2 => Fully implicit
!                     = 3 => Adaptive (only implicit where the adaptation time
!                            scale is small)
!*iopt_sed_toteq*     INTEGER Selects total load equation
!                     = 0 => disabled
!                     = 1 => Engelund-Hansen (1967) standard formulation
!                     = 2 => Engelund-Hansen (1967) using Chollet-Cunge (1979)
!                            formulation for the Shield parameter
!                     = 3 => Ackers-White (1973)
!                     = 4 => Madsen-Grant (1976)
!                     = 5 => Wu et al (2000)
!                     = 6 => Van Rijn (2003)
!                     = 7 => Van Rijn (2007)
!*iopt_sed_type*      INTEGER Sediment type
!                     = 1 => sand
!                     = 2 => mud
!*iopt_sed_vadv*      INTEGER Selects settling and type of vertical advection
!                             scheme
!                     = 0 => settling disabled
!                     = 1 => upwind
!                     = 2 => central
!                     = 3 => TVD
!                 (if iopt_adv_scal>0 then iopt_sed_vadv=iopt_adv_scal or zero)
!*iopt_sed_wave_diff* INTEGER Selects turbulent diffusion coefficient due to
!                            waves
!                     = 0 => No diffusion coefficient
!                     = 1 => Diffusion according to Van Rijn (2007b)
!*iopt_sed_ws*        INTEGER Selects settling velocity equation
!                     = 1 => user defined uniform value for each fraction
!                     = 2 => Camenen (2007) sand
!                     = 3 => Camenen (2007) mud
!                     = 4 => Stokes (1851) small spherical particles
!                     = 5 => Soulsby (1997)
!                     = 6 => Zhang and Xie (1993)
!*iopt_sed_ws_floc*   INTEGER Selects flocculation effect on settling velocity
!                     = 0 => No flocculation
!                     = 1 => Van Leussen (1994)
!                     = 2 => Van Rijn (2007)
!*iopt_sed_ws_hindset*INTEGER Selects hindered settling equation
!                     = 0 => No hindered settling
!                     = 1 => Richardson & Zaki (1954
!                     = 2 => Winterwerp & van Kesteren (2004)
!                     = 3 => Integrated Van Leussen (1994) Van Rijn (2007)
!*iopt_sed_ws_lim*    INTEGER Disables/enables limitation of the settling
!                             velocity for shallow waters
!
!************************************************************************
!

END MODULE sedswitches
