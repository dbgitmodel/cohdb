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

MODULE turbpars
!************************************************************************
!
! *turbpars* Turbulence parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)turbpars.f90  V2.1.1
!
! $Date: 2015-04-10 11:30:34 +0200 (Fri, 10 Apr 2015) $
!
! $Revision: 841 $
!
! Description - 
!
!************************************************************************
!
IMPLICIT NONE

!---stability functions
INTEGER :: ib22
REAL :: cfequil(5), cfstabtke(3), cfstab1(8), cfstab2(11),&
      & f0stabmom, f0stabscal, f0stabtke, eps0, skeps = 0.09, sq_my = 0.2,&
      & c_sk = 0.15

!---limiting conditions
REAL :: alphaM_max, alphaN_max, alphaN_min, alphaN_sl, tkelim = 1.0E-06,&
      & tkemin = 1.0E-014, dissipmin = 1.0E-012, zlmixmin = 1.7E-010

!---t.k.e.-equation
REAL :: sigma_k = 1.0, wfltke = 0.0

!---kl-equation
REAL :: e1_my = 1.8, e2_my = 1.33, e3_my = 1.0, sigma_kl

!---eps-equation
REAL :: c1_eps = 1.44, c2_eps = 1.92, c31_eps = 0.2, c32_eps = 1.0, sigma_eps

!---roughness lengths
REAL :: zrough_bot = 0.0, zrough_sur = 0.0

!---KPP scheme for internal waves
REAL :: riccrit_iw = 0.7, vdifmom_iw = 1.0E-04, vdifscal_iw = 5.0E-05,&
      & vdifshear_iw = 0.005

!---Pacanowski-Philander relations
REAL :: alpha_pp = 5.0, expmom_pp = 2.0, vbmom_pp = 1.0E-04,&
      & vbscal_pp = 1.0E-05, vmax_pp = 3.0, v0dif_pp = 0.01

!---Munk-Anderson relations
REAL :: alpha_ma = 10.0, beta_ma = 3.33, expmom_ma = 0.5, expscal_ma = 1.5,&
      & vmaxmom_ma = 3.0, vmaxscal_ma = 4.0, v0dif_ma = 0.06

!---algebraic flow-dependent relations
REAL :: delta1_ad = 0.0, delta2_ad = 0.0, r1_ad = 1.0, r2_ad = 1.0,&
      & k1_ad = 0.0025, k2_ad = 2.0E-05, lambda_ad = 0.0, cnu_ad = 2.0,&
      & omega1_ad = 1.0E-04

!---algebraic mixing length formulations
REAL :: alpha_Black = 0.2, beta_Xing = 2.0 

SAVE

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*ib22*         INTEGER Internal switch to simplify calculations
!*cfequil*      REAL    Coefficients equilibrium method
!*cfstabtke*    REAL    Coefficients in expression for t.k.e. stability
!                       function
!*cfstab1*      REAL    Coefficients in expressions for quasi-equilibrium
!                       stability functions
!*cfstab2*      REAL    Coefficients in expressions for non-equilibrium
!                       stability functions
!*f0stabmom*    REAL    Neutral stability coefficient for momentum
!*f0stabscal*   REAL    Neutral stability coefficient for scalars
!*f0stabtke*    REAL    Neutral stability coefficient for t.k.e.
!*eps0*         REAL    Parameter relating dissipation rate to mixing length
!*skeps*        REAL    Neutral stability coefficient for t.k.e. (k-eps model)
!*sq_my*        REAL    Mellor-Yamada coefficient S_q
!*c_sk*         REAL    Constant in Daly-Harlow formulation
!*alphaM_max*   REAL    Maximum value of alpha_M used in limiting conditions
!*alphaN_max*   REAL    Maximum value of alpha_N used in limiting conditions
!*alphaN_min*   REAL    Minimum value of alpha_N in unstable case
!*alphaN_sl*    REAL    Slope of limting curve (stable case)
!*tkelim*       REAL    Minimum value of t.k.e. in limiting conditions   [J/kg]
!*tkemin*       REAL    Numerical lower limit for turbulence energy      [J/kg]
!*dissipmin*    REAL    Numerical lower limit for dissipation rate       [W/kg]
!*zlmixmin*     REAL    Numerical lower limit for mixing length             [m]
!*sigma_k*      REAL    Ratio of diffusion coefficients for momentum and t.k.e.
!                       if iopt_turb_stab_tke=1
!*wfltke*       REAL    Surface wave factor in surface flux of turbulence
!                       energy
!*e1_my*        REAL    Mellor_Yamada coefficient E_1 
!*e2_my*        REAL    Mellor_Yamada coefficient E_2 
!*e3_my*        REAL    Mellor_Yamada coefficient E_3
!*sigma_kl*     REAL    Ratio of diffusion coefficients for t.k.e and kl
!*c1_eps*       REAL    Parameter c_1epsilon in dissipation equation
!*c2_eps*       REAL    Parameter c_2epsilon in dissipation equation
!*c31_eps*      REAL    Parameter c_3epsilon in dissipation equation
!                       (stable case)
!*c32_eps*      REAL    Parameter c_3epsilon in dissipation equation
!                       (unstable case)
!*sigma_eps*    REAL    Ratio of diffusion coefficients for t.k.e. and
!                       dissipation
!*zrough_bot*   REAL    Bottom roughness length                             [m]
!*zrough_sur*   REAL    Surface roughness length                            [m]
!*riccrit_iw*   REAL    Critical Richardson number in internal wave mixing
!                       scheme
!*vdifmom_iw*   REAL    Momentum diffusion coefficient due to internal wave
!                       breaking                                        [m^2/s]
!*vdifscal_iw*  REAL    Scalar diffusion coefficient due to internal wave
!                       breaking                                        [m^2/s]
!*vdifshear_iw* REAL    Maximum diffusion due to resolved shear below the
!                       mixed layer                                     [m^2/s]
!*alpha_pp*     REAL    Factor alpha_p in PP formulation
!*expmom_pp*    REAL    Exponential factor (n_p) in PP formulation
!*vbmom_pp*     REAL    Background momentum diffusion coefficient (nu_bp) in
!                       PP formulation                                  [m^2/s]
!*vbscal_pp*    REAL    Background scalar diffudion coefficient (lambda_bp) in 
!                       PP formulation                                  [m^2/s]
!*vmax_pp*      REAL    Limiting factor nu_max in PP formulation 
!*v0dif_pp*     REAL    Factor nu_0p in PP formulation                  [m^2/s]
!*alpha_ma*     REAL    Factor alpha_m im MA formulation
!*beta_ma*      REAL    Factor beta_m in MA formulation
!*expmom_ma*    REAL    Exponential factor n_1 in MA formulation
!*expscal_ma*   REAL    Exponential factor n_2 in MA formulation
!*vmaxmom_ma*   REAL    Factor nu_max in MA formulation
!*vmaxscal_ma*  REAL    Factor lambda_max in MA formulation
!*v0dif_ma*     REAL    Factor nu_0m in MA formulation                  [m^2/s]
!*delta1_ad*    REAL    Factor delta_1 in AD formulation
!*delta2_ad*    REAL    Factor delta_2 in AD formulation
!*r1_ad*        REAL    Factor r_1 in AD formulation
!*r2_ad*        REAL    Factor r_2 in AD formulation
!*k1_ad*        REAL    Factor K_1 in AD formulation
!*k2_ad*        REAL    Factor K_2 in AD formulation
!*lambda_ad*    REAL    Factor lambda_* in AD formulation                   [m]
!*cnu_ad*       REAL    Factor C_nu in AD formulation
!*omega1_ad*    REAL    Frequency omega_1 in AD formulation              [s^-1]
!*alpha_Black*  REAL    Constant in Blackadar's asymptotic length scale
!*beta_Xing*    REAL    Constant in Xing's exponential bottom factor
!
!  PP = Pacanowski and Philander scheme
!  MA = Munk and Anderson scheme
!  AD = Davies flow-dependent formulations
!
!************************************************************************
!

END MODULE turbpars
