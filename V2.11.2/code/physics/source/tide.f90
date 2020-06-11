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

MODULE tide
!************************************************************************
!
! *tide* Parameters and arrays for tidal forcing
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)tide.f90  V2.0
!
! $Date: 2015-04-10 11:30:34 +0200 (Fri, 10 Apr 2015) $
!
! $Revision: 841 $
!
! Description - 
!
!************************************************************************
!
USE syspars

IMPLICIT NONE


!---number of tidal constituents
INTEGER :: nconastro = 0, nconobc = 0

!---tidal indices
INTEGER, DIMENSION(MaxAstroTides) :: index_astro = 0
INTEGER, DIMENSION(MaxConstituents) :: index_obc = 0

!---nodal factors
REAL, ALLOCATABLE, DIMENSION(:) :: fnode_astro, fnode_obc

!---tidal phases
REAL, ALLOCATABLE, DIMENSION(:) :: phase_astro, phase_obc

!---tidal force
REAL, ALLOCATABLE, DIMENSION(:,:) :: fxastro, fyastro

!---key ids for tidal constituents
INTEGER, PARAMETER ::&
    & icon_MS0   =  1, icon_Sa    =  2, icon_Ssa   =  3, icon_058    =  4,&
    & icon_Msm   =  5, icon_Mm    =  6, icon_Msf   =  7, icon_Mf     =  8,&
    & icon_083   =  9, icon_Mt    = 10, icon_093   = 11, icon_095    = 12,&
    & icon_alpha1= 13, icon_2Q1   = 14, icon_sigma1= 15, icon_Q1     = 16,&
    & icon_rho1  = 17, icon_O1    = 18, icon_tau1  = 19, icon_beta1  = 20,&
    & icon_M1    = 21, icon_chi1  = 22, icon_pi1   = 23, icon_P1     = 24,&
    & icon_S1    = 25, icon_K1    = 26, icon_psi1  = 27, icon_phi1   = 28,&
    & icon_theta1= 29, icon_J1    = 30, icon_SO1   = 31, icon_OO1    = 32,&
    & icon_KQ1   = 33, icon_OQ2   = 34, icon_eps2  = 35, icon_2N2    = 36,&
    & icon_mu2   = 37, icon_238   = 38, icon_244   = 39, icon_N2     = 40,&
    & icon_246   = 41, icon_nu2   = 42, icon_248   = 43, icon_gamma2 = 44,&
    & icon_H1    = 45, icon_M2    = 46, icon_H2    = 47, icon_lambda2= 48,&
    & icon_L2    = 49, icon_T2    = 50, icon_S2    = 51, icon_R2     = 52,&
    & icon_K2    = 53, icon_eta2  = 54, icon_295   = 55, icon_M3     = 56,&
    & icon_2SM2  = 57, icon_2MK3  = 58, icon_SO3   = 59, icon_MK3    = 60,&
    & icon_SK3   = 61, icon_MN4   = 62, icon_M4    = 63, icon_MS4    = 64,&
    & icon_MK4   = 65, icon_S4    = 66, icon_2MN6  = 67, icon_M6     = 68,&
    & icon_MSN6  = 69, icon_2MS6  = 70, icon_2SM6  = 71, icon_S6     = 72,&
    & icon_M8    = 73, icon_2MSN8 = 74, icon_3MS8   = 75, icon_2MS8   = 76,&
    & icon_S8    = 77

!---tidal species index
INTEGER, DIMENSION(MaxConstituents) :: ispec_tides
DATA ispec_tides /12*0, 21*1, 22*2, 1*3, 1*2, 4*3, 5*4, 6*6, 5*8/

!---amplitudes of equilibrium tide
REAL, DIMENSION(MaxAstroTides) :: astro_ampl = &
 & (/0.198419, 0.003103, 0.019542, 0.001142, 0.004239, 0.022191, 0.003677, &
 &   0.042016, 0.001526, 0.008049, 0.001287, 0.001066, 0.000749, 0.002565, &
 &   0.003098, 0.019387, 0.003685, 0.101266, 0.001325, 0.000749, 0.007965, &
 &   0.001522, 0.002754, 0.047129, 0.001116, 0.142408, 0.001132, 0.002028, &
 &   0.001526, 0.007965, 0.001321, 0.004361, 0.000834, 0.000695, 0.001804, &
 &   0.006184, 0.007463, 0.000502, 0.000394, 0.046735, 0.000436, 0.008877, &
 &   0.000409, 0.000734, 0.000842, 0.244102, 0.000746, 0.001800, 0.006899, &
 &   0.006636, 0.113572, 0.000950, 0.030875, 0.001727, 0.000452, 0.003455/)    

!---amplitude corrections for earth tide
REAL, DIMENSION(MaxAstroTides) :: astro_earth
DATA astro_earth /15*0.693, 0.6946, 0.6948, 0.6950, 0.6956, 0.693, 0.6962, &
 & 0.6994, 0.7027, 0.7059, 0.7126, 0.7364, 0.5285, 0.6657, 0.6784, 0.6911, &
 & 0.693, 0.6925, 24*0.693/

!---tidal frequencies
REAL, DIMENSION(MaxConstituents) :: tidal_spectrum = &
&(/0.0000000E+00, 1.9909688E-07, 3.9821277E-07, 5.9730965E-07, 2.2859987E-06,&
 & 2.6392030E-06, 4.9252017E-06, 5.3234145E-06, 7.6094132E-06, 7.9626175E-06,&
 & 1.0248616E-05, 1.0601821E-05, 6.0033339E-05, 6.2319338E-05, 6.2672542E-05,&
 & 6.4958541E-05, 6.5311745E-05, 6.7597744E-05, 6.7995957E-05, 6.9883743E-05,&
 & 7.0281956E-05, 7.0635160E-05, 7.2323849E-05, 7.2522946E-05, 7.2722062E-05,&
 & 7.2921159E-05, 7.3120255E-05, 7.3319371E-05, 7.5207157E-05, 7.5560362E-05,&
 & 7.7846360E-05, 7.8244573E-05, 8.0883776E-05, 1.3260129E-04, 1.3295450E-04,&
 & 1.3524050E-04, 1.3559370E-04, 1.3579280E-04, 1.3768060E-04, 1.3787970E-04,&
 & 1.3807880E-04, 1.3823290E-04, 1.3843200E-04, 1.4016570E-04, 1.4031981E-04,&
 & 1.4051890E-04, 1.4071800E-04, 1.4280490E-04, 1.4315811E-04, 1.4524501E-04,&
 & 1.4544410E-04, 1.4564320E-04, 1.4584232E-04, 1.4848152E-04, 1.5116573E-04,&
 & 2.1077835E-04, 1.5036931E-04, 2.0811665E-04, 2.1304185E-04, 2.1344006E-04,&
 & 2.1836526E-04, 2.7839860E-04, 2.8103781E-04, 2.8596301E-04, 2.8636122E-04,&
 & 2.9088821E-04, 4.1891750E-04, 4.2155671E-04, 4.2384271E-04, 4.2648191E-04,&
 & 4.3140711E-04, 4.3633231E-04, 5.6207561E-04, 5.6436161E-04, 5.6700081E-04,&
 & 5.7192601E-04, 5.8177642E-04/)

!---names of tidal frequencies
CHARACTER (LEN=lenfreq), DIMENSION(MaxConstituents) :: tidal_freq_names = &
    & (/'MS0    ' ,'Sa     ' ,'Ssa    ', '058    ', 'Msm    ', 'Mm     ',&
    &   'Msf    ', 'Mf     ', '083    ', 'Mt     ', '093    ', '095    ',&
    &   'alpha1 ', '2Q1    ', 'sigma1 ', 'Q1     ', 'rho1   ', 'O1     ',&
    &   'tau1   ', 'beta1  ', 'M1     ', 'chi1   ', 'pi1    ', 'P1     ',&
    &   'S1     ', 'K1     ', 'psi1   ', 'phi1   ', 'theta1 ', 'J1     ',&
    &   'SO1    ', 'OO1    ', 'KQ1    ', 'OQ2    ', 'eps2   ', '2N2    ',&
    &   'mu2    ', '238    ', '244    ', 'N2     ', '246    ', 'nu2    ',&
    &   '248    ', 'gamma2 ', 'H1     ', 'M2     ', 'H2     ', 'lambda2',&
    &   'L2     ', 'T2     ', 'S2     ', 'R2     ', 'K2     ', 'eta2   ',&
    &   '295    ', 'M3     ', '2SM2   ', '2MK3   ', 'SO3    ', 'MK3    ',&
    &   'SK3    ', 'MN4    ', 'M4     ', 'MS4    ', 'MK4    ', 'S4     ',&
    &   '2MN6   ', 'M6     ', 'MSN6   ', '2MS6   ', '2SM6   ', 'S6     ',&
    &   'M8     ', '2MSN8  ', '3MS8   ', '2MS8   ', 'S8     '/)

SAVE

!
! Name             Type    Purpose
!------------------------------------------------------------------------------
!*nconastro*       INTEGER Number of tidal constituents for astronomical force
!*nconobc*         INTEGER Number of tidal constituents at open boundaries
!*index_astro*     INTEGER Indices of astronomical force constituents
!*index_obc*       INTEGER Indices of tidal constituents at open boundaries
!*fnode_astro*     REAL    Nodal factors in components of tidal force
!*fnode_obc*       REAL    Nodal factors of tidal constituents at open
!                          boundaries
!*phase_astro*     REAL    Astronomical phases of tidal force constituents
!                          with respect to the reference longitude        [rad]
!*phase_obc*       REAL    Astronomical phases of open boundary tidal
!                          constituents with respect to the reference
!                          longitude                                      [rad]
!*fxastro*         REAL    X-component of tidal force                   [m^2/s]
!*fyastro*         REAL    Y-component of tidal force                   [m^2/s]
!*ispec_tides*     INTEGER Tidal species
!*astro_ampl*      REAL    Amplitude of equilibrium tide                    [m]
!*astro_earth*     REAL    Elasticity factor for Earth tides (1+kn-hn)
!*tidal_spectrum*  REAL    Frequencies of tidal constituents            [rad/s]
!*tidal_freq_names*CHAR    Names of tidal frequencies
!
!************************************************************************
!

END MODULE tide
