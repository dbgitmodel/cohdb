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

MODULE sedpars
!************************************************************************
!
! *sedpars* Sediment parameters
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
USE syspars
 
IMPLICIT NONE

LOGICAL :: sedstep
INTEGER::  icsed = 1, maxitbartnicki = 100, nb = 1, nf = 1, nfp, &
         & nrquad_sed = 7, nrquad_wav = 10
REAL ::  alpha_VR = 2.19, a_leussen = 0.02, beta_sed_cst = 1.0, &
       & beta_sed_max = 1.5, beta_sed_min = 1.0, b_leussen = 0.0024, &
       & cgel = 0.0, cmax = 0.65, coef_bed_grad = 1.3, deltsed, &
       & floc_VR_max = 10.0, floc_VR_min = 1.0, height_c_cst = 0.01, &
       & n_RichZaki = 4.6, parth_coef = 2.0E-04, &
       & parth_exp = 1.0, theta_sed_src = 0.501, wslimfac = 0.5, &
       & wu_exp = -0.6, zrough_sed_cst = 0.0, zrough_grain = 0.05, &
       & z0_coef = 30.0

!---flocculation module
REAL :: agg_alpha = 0.1, brk_es = 1.0E-04, brk_ffy = 1.0E-10, brk_frac = 0.1, &
       & brk_p = 1.0, brk_q = 1.0, floc_ncmax = 100.0, floc_ncmin = 1.0001, &
       & nfrdim = 2.0

SAVE

! Name            Type     Purpose
!------------------------------------------------------------------------------
!*agg_alpha*      REAL     Collision efficiency factor                       [-]
!*alpha_VR*       REAL     Value of alpha in Van Rijn's 2007 flocculation 
!                          equation                         
!*a_leussen*      REAL     Coefficient a in Van Leussen (1994) flocculation
!                          equation                                          [s]
!*beta_sed_cst*   REAL     Constant ratio of sediment diffusivity and eddy
!                          viscosity                       
!*beta_sed_max*   REAL     Maximum value for the ratio of sediment
!                          diffusivity to eddy viscosity                     [-]
!*beta_sed_min*   REAL     Minimum value for the ratio of sediment
!                          diffusivity to eddy viscosity       
!*b_leussen*      REAL     Coeffcient b in Van Leussen (1994) flocculation
!                          equation                                        [s^2]
!*brk_es*         REAL     Efficiency factor for breaking              [s^0.5/m]
!*brk_ffy*        REAL     Yield strength of flocs                          [Pa]
!*brk_fac*        REAL     Fraction of mircoflocs by breakage                [-]
!*brk_p*          REAL     Empirical exponential factor of breakage kinetics [-]
!*brk_q*          REAL     Empirical exponential factor of breakage kinetics [-]
!*cgel*           REAL     Volumetric gelling concentration (hindered
!                          settlingof mud)                             [m^3/m^3]
!*cmax*           REAL     Volumetric maximum concentration (sand)     [m^3/m^3]
!*coef_bed_grad*  REAL     Coefficient in the bed gradient equation of Koch and
!                          Flokstra                                          [-]
!*deltsed*        REAL     Time step for sediment update                     [s]
!*floc_ncmax*     REAL     Maximum limit for floc_nc                         [-]
!*floc_ncmin*     REAL     Minimum limit for floc_nc                         [-]
!*floc_VR_max*    REAL     Maximum of Van Rijn (2007) flocculation factor    [-]
!*floc_VR_min*    REAL     Minimum of Van Rijn (2007) flocculation factor    [-]
!*height_c_cst*   REAL     Non-dimensional height of the near bed sediment flux
!                          (normalised with water depth)                     [-]
!*icsed*          INTEGER  Number of 3-D time steps per sediment time step
!*maxitbartnicki* INTEGER  Maximum number of iterations in bartnicki filter
!*nb*             INTEGER  Number of bed layers
!*nf*             INTEGER  Number of sediment fractions
!*nfp*            INTEGER  Number of primary particle fractions (=1 in case
!                          of flocculation, = nf otherwise)   
!*nfrdim*         REAL     Fractal dimension of macroflocs                   [-]
!*nrquad_sed*     INTEGER  Number of points used to calculate a depth average
!*nrquad_wav*     INTEGER  Number of timesteps used for numerical integration 
!                          of sediment transport over 1 wave period
!*n_RichZaki*     REAL     Exponent in Richardson and Zaki Equation (hindered
!                          settling                                          [-]
!*parth_coef*     REAL     Coefficient in the erosion equation of Partheniades
!                          (mud)                                      [kg/m^2/s]
!*parth_exp*      REAL     Exponent in the erosion equation of Partheniades
!                          (mud)                                             [-]
!*sedstep*        LOGICAL  .TRUE. at sediment transport time step
!*theta_sed_src*  REAL     Value of theta used in the semi-implicit
!                          discretization of the deposition term in 2D       [-]
!*wslimfac*       REAL     Factor used to set an upper limit for the fall
!                          velocity in shallow water                         [-]
!*wu_exp*         REAL     Exponent in the hiding factor exquation of Wu et al.
!                          (2000)                                            [-]
!*zrough_grain*   REAL     Proportionality factor in the definition of the
!                          grain size roughness
!*zrough_sed_cst* REAL     Uniform skin bottom roughnes                      [m]
!                          grain size roughness
!*z0_coef*        REAL     Factor with which z0 is multiplied to determine
!                          minimum depth for averaging bed boundary condition
!                                                                            [-]
!
!************************************************************************
!


END MODULE sedpars
