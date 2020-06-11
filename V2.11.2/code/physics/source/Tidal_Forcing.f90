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
! *Tidal_Forcing* Astronomical parameters and tidal force
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Tidal_Forcing.f90  V2.11.2
!
! $Date: 2018-09-28 13:01:37 +0200 (Fri, 28 Sep 2018) $
!
! $Revision: 1188 $
!
! Description - 
!
! Reference -
!
! Subroutines - astro_params, astronomical_tides
!
!************************************************************************
!

!========================================================================

SUBROUTINE astro_params(index_tides,fnode_tides,phase_tides,notides,&
                      & dlonref_phase,cdatetime_tides)
!************************************************************************
!
! *astro_params* Update nodal factors and astronomical arguments
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Tidal_Forcing.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - astronomical_tides, harmonic_analysis_data,
!                   pressure_gradient_1d, update_2dobc_data
!
! External calls -
!
! Module calls - add_secs_to_date, complex_polar, convert_date, day_number
!
!************************************************************************
!
USE iopars
USE switches
USE syspars
USE tide
USE timepars
USE math_library, ONLY: complex_polar
USE time_routines, ONLY: add_secs_to_date, convert_date, day_number, &
                       & log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: cdatetime_tides
INTEGER, INTENT(IN) :: notides
REAL, INTENT(IN) :: dlonref_phase
INTEGER, INTENT(IN), DIMENSION(notides) :: index_tides
REAL, INTENT(OUT), DIMENSION(notides) :: fnode_tides, phase_tides

!
! Name             Type    Purpose
!------------------------------------------------------------------------------
!*index_tides*     INTEGER Frequency indices
!*fnode_tides*     REAL    Nodal correction factors
!*phase_tides*     REAL    Astronomical arguments                        [rad]
!*notides*         INTEGER Number of tidal constituents
!*dlonref_phase*   REAL    Reference longitude for astronomical phase
!                                                                    [degrees]
!*cdatetime_tides* CHAR    Date/time at local time zone
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: icon, iicon, ish, iyear, ndays, nyears
INTEGER :: ish_O1 = 1, ish_K1 = 2, ish_N2 = 3, ish_M2 = 4, ish_S2 = 5, &
         & ish_K2 = 6
REAL :: cosN, cosp_ps, cos2N, cos2p, cos2pN, cos2ps, cos2p2N, cos2p_N, &
      & daynum, dref, dsecs, fnode, phase, sinN, sinp_ps, sin2N, sin2p, &
      & sin2pN, sin2ps, sin2p2N, sin2p_N
REAL (KIND=kndrlong) :: hlon, Nlon, plon, pslon, slon, tjul
INTEGER, DIMENSION(7) :: idatetime_tides, IGDateTime
REAL, DIMENSION(6) :: fnode_sh, phase_sh, r1_sh, r2_sh
REAL, DIMENSION(notides) :: r1, r2


procname(pglev+1) = 'astro_params'
CALL log_timer_in()

!
!1. Initialise
!-------------
!

fnode_tides = 1.0; phase_tides = 0.0

!
!2. Astronomical parameters
!--------------------------
!
!---Greenwich time
idatetime_tides = convert_date(cdatetime_tides)
CALL add_secs_to_date(idatetime_tides,IGDateTime,-3600,time_zone)

!---time in Julian centuries since Greenwich noon on 1899/12/31
iyear = IGDateTime(1)
CALL day_number(IGDateTime,daynum)
IF (iyear.LT.1975) THEN
   nyears = IGDateTime(1) - 1900
   ndays = INT(daynum)+365*nyears+(nyears-1)/4
   tjul = (0.5+ndays)/36525.0_kndrlong
ELSE
   ndays = INT(daynum)+365*(iyear-1975)+(iyear-1973)/4
   tjul = (27393.5_kndrlong+ndays)/36525+&
        & (5.28E-04_kndrlong+3.56E-08_kndrlong*(ndays+1))/36525
ENDIF

!---astronomical longitudes
slon = 270.434358_kndrlong+MOD(481267*tjul,360.0_kndrlong)+&
     & 0.88314137_kndrlong*tjul-0.001133_kndrlong*tjul*tjul+&
     & 1.9E-06_kndrlong*(tjul**3)
hlon = 279.69668_kndrlong+MOD(36000*tjul,360.0_kndrlong)+&
     & 0.768925485_kndrlong*tjul+3.03E-04_kndrlong*tjul*tjul
plon = 334.329653_kndrlong+MOD(4069*tjul,360.0_kndrlong)+0.0340329575*tjul-&
     & 0.10325_kndrlong*tjul*tjul-1.2E-05_kndrlong*(tjul**3)
Nlon = -259.16_kndrlong+MOD(1934*tjul,360.0_kndrlong)+0.14_kndrlong*tjul-&
     & 0.0021_kndrlong*tjul*tjul
pslon = 281.22083_kndrlong+MOD(1.71902_kndrlong*tjul,360.0_kndrlong)+&
      & 0.00045_kndrlong*tjul*tjul+3.0E-06_kndrlong*(tjul**3)

!---convert to radians
slon = MOD(degtorad_d*slon,twopi_d)
hlon = MOD(degtorad_d*hlon,twopi_d)
plon = MOD(degtorad_d*plon,twopi_d)
Nlon = MOD(degtorad_d*Nlon,twopi_d)
pslon = MOD(degtorad_d*pslon,twopi_d)

!---auxiliary parameters
cosN = COS(Nlon); sinN = SIN(Nlon)
cos2N = COS(2.0*Nlon); sin2N = SIN(2.0*Nlon)
cos2p = COS(2.0*plon); sin2p = SIN(2.0*plon)
cos2ps = COS(2.0*pslon); sin2ps = SIN(2.0*pslon)
cos2pN = COS(2.0*plon+Nlon); sin2pN = SIN(2.0*plon+Nlon)
cos2p2N = COS(2.0*(plon+Nlon)); sin2p2N = SIN(2.0*(plon+Nlon))
cos2p_N = COS(2.0*plon-Nlon); sin2p_N = SIN(2.0*plon-Nlon)
cosp_ps = COS(plon-pslon); sinp_ps = SIN(plon-pslon)

!
!3. Astronomical phases without nodal correction at Greenwich midnigth
!---------------------------------------------------------------------
!

icon_300: DO icon=1,notides
   iicon = index_tides(icon)
   IF (iicon.GT.0) THEN
      SELECT CASE (iicon)
      CASE(icon_MS0); phase_tides(icon) = 0.0
      CASE(icon_Sa); phase_tides(icon) = hlon-pslon
      CASE(icon_Ssa); phase_tides(icon) = 2.0*hlon
      CASE(icon_058); phase_tides(icon) = 3.0*hlon-pslon
      CASE(icon_Msm); phase_tides(icon) = slon-2.0*hlon+plon
      CASE(icon_Mm); phase_tides(icon) =  slon-plon
      CASE(icon_Msf); phase_tides(icon) = 2.0*(slon-hlon)
      CASE(icon_Mf); phase_tides(icon) =  2.0*slon
      CASE(icon_083); phase_tides(icon) = 3.0*slon-2.0*hlon+plon
      CASE(icon_Mt); phase_tides(icon) = 3.0*slon-plon
      CASE(icon_093); phase_tides(icon) = 2.0*(2.0*slon-hlon)
      CASE(icon_095); phase_tides(icon) = 2.0*(2.0*slon-plon)
      CASE(icon_alpha1); phase_tides(icon) = -5.0*slon+3.0*hlon+plon-halfpi
      CASE(icon_2Q1); phase_tides(icon) = -4.0*slon+hlon+2.0*plon-halfpi
      CASE(icon_sigma1); phase_tides(icon) = -4.0*slon+3.0*hlon-halfpi
      CASE(icon_Q1); phase_tides(icon) = -3.0*slon+hlon+plon-halfpi
      CASE(icon_rho1); phase_tides(icon) = 3.0*(hlon-slon)-plon-halfpi
      CASE(icon_O1); phase_tides(icon) = -2.0*slon+hlon-halfpi
      CASE(icon_tau1); phase_tides(icon) = -2.0*slon+3.0*hlon+halfpi
      CASE(icon_beta1); phase_tides(icon) = -slon-hlon+plon+halfpi
      CASE(icon_M1); phase_tides(icon) = -slon+hlon+plon+halfpi
      CASE(icon_chi1); phase_tides(icon) = -slon+3.0*hlon-plon+halfpi
      CASE(icon_pi1); phase_tides(icon) = -2.0*hlon+pslon-halfpi
      CASE(icon_P1); phase_tides(icon) = -hlon-halfpi
      CASE(icon_S1); phase_tides(icon) = pslon+halfpi
      CASE(icon_K1); phase_tides(icon) = hlon+halfpi
      CASE(icon_psi1); phase_tides(icon) = 2.0*hlon-pslon+halfpi
      CASE(icon_phi1); phase_tides(icon) = 3.0*hlon+halfpi
      CASE(icon_theta1); phase_tides(icon) = slon-hlon+plon+halfpi
      CASE(icon_J1); phase_tides(icon) = slon+hlon+halfpi
      CASE(icon_SO1); phase_tides(icon) = 2.0*slon-hlon+halfpi
      CASE(icon_OO1); phase_tides(icon) = 2.0*slon+hlon+halfpi
      CASE(icon_KQ1); phase_tides(icon) = 3.0*slon+hlon-plon+halfpi
      CASE(icon_OQ2); phase_tides(icon) = -5.0*slon+2.0*hlon+3.0*plon
      CASE(icon_eps2); phase_tides(icon) = -5.0*slon+4.0*hlon+plon
      CASE(icon_2N2); phase_tides(icon) = -2.0*(2.0*slon-hlon-plon)
      CASE(icon_mu2); phase_tides(icon) = -4.0*(slon-hlon)
      CASE(icon_238); phase_tides(icon) = -4.0*slon+5.0*hlon-pslon
      CASE(icon_244); phase_tides(icon) = -3.0*slon+hlon+plon+pslon+pi
      CASE(icon_N2); phase_tides(icon) = -3.0*slon+2.0*hlon+plon
      CASE(icon_246); phase_tides(icon) = -3.0*(slon-hlon)+plon-pslon
      CASE(icon_nu2); phase_tides(icon) = -3.0*slon+4.0*hlon-plon
      CASE(icon_248); phase_tides(icon) = -3.0*slon+5.0*hlon-plon-pslon
      CASE(icon_gamma2); phase_tides(icon) =-2.0*(slon-plon)+pi
      CASE(icon_H1); phase_tides(icon) = -2.0*slon+hlon+pslon+pi
      CASE(icon_M2); phase_tides(icon) = -2.0*(slon-hlon)
      CASE(icon_H2); phase_tides(icon) = -2.0*slon+3.0*hlon-pslon
      CASE(icon_lambda2); phase_tides(icon) = -slon+plon+pi
      CASE(icon_L2); phase_tides(icon) = -slon+2.0*hlon-plon+pi
      CASE(icon_T2); phase_tides(icon) = -hlon+pslon
      CASE(icon_S2); phase_tides(icon) = 0.0
      CASE(icon_R2); phase_tides(icon) = hlon-pslon+pi
      CASE(icon_K2); phase_tides(icon) = 2.0*hlon
      CASE(icon_eta2); phase_tides(icon) = slon+2.0*hlon-plon
      CASE(icon_295); phase_tides(icon) = 2.0*(slon+hlon)
      CASE(icon_M3); phase_tides(icon) = -3.0*(slon-hlon)
      END SELECT
   ENDIF
ENDDO icon_300

!---phases for overtides
phase_sh(ish_O1) = -2.0*slon+hlon-halfpi
phase_sh(ish_K1) = hlon+halfpi
phase_sh(ish_N2) = -3.0*slon+2.0*hlon+plon
phase_sh(ish_M2) = -2.0*(slon-hlon)
phase_sh(ish_S2) = 0.0
phase_sh(ish_K2) = 2.0*hlon

!
!4. Nodal corrections
!--------------------
!
!4.1 Cosine factors
!------------------
!

icon_410: DO icon=1,notides
   iicon = index_tides(icon)
   IF (iicon.GT.0) THEN
      SELECT CASE (iicon)
      CASE(icon_MS0); r1(icon) = -0.0888*cosN
      CASE(icon_Sa); r1(icon) = -0.0529*cos2ps-0.0021*cosN
      CASE(icon_Ssa)
         r1(icon) = -0.0249*cosN+0.0103*cos2p-0.0055*cos2N+0.0039*cos2ps
      CASE(icon_058); r1(icon) = -0.0166*cosN
      CASE(icon_Msm); r1(icon) = -0.1369*cosN-0.0104*cos2pN-0.003*cos2p2N
      CASE(icon_Mm)
         r1(icon) = -0.1308*cosN-0.0534*cos2p-0.0219*cos2pN-0.006*cos2p2N
      CASE(icon_Msf); r1(icon) = 0.0068*cosN-0.0069*cos2p
      CASE(icon_Mf)
         r1(icon) = 0.4147*cosN+0.0432*cos2p+0.0387*cos2N-0.029*cos2p_N-&
                  & 0.0023*cos2pN
      CASE(icon_083)
         r1(icon) = 0.4132*cosN+0.3802*cos2p+0.0372*cos2pN+0.0372*cos2N-&
                  & 0.0248*cos2p_N
      CASE(icon_Mt)
         r1(icon) = 0.4146*cosN+0.04*cos2N+0.018*cos2p-0.0039*cos2p2N-&
                  & 0.0031*COS(2.0*plon-Nlon+pslon)-&
                  & 0.0031*COS(2.0*plon-Nlon-pslon)-&
                  & 0.0016*COS(2.0*plon+3.0*Nlon)
      CASE(icon_093); r1(icon) = 0.4118*cosN+0.0539*cos2p+0.0392*cos2N
      CASE(icon_095); r1(icon) = 0.4142*cosN+0.0414*cos2N
      CASE(icon_alpha1); r1(icon) = 0.9107*cosN
      CASE(icon_2Q1); r1(icon) = 0.1883*cosN-0.006*cos2p2N-0.0045*cos2N
      CASE(icon_sigma1); r1(icon) = 0.1883*cosN-0.0087*cos2p-0.005*cos2N
      CASE(icon_Q1)
         r1(icon) = 0.1887*cosN-0.0058*cos2N-0.0038*cos2pN-0.0028*cos2p
      CASE(icon_rho1)
         r1(icon) = 0.1887*cosN-0.0577*cos2p+0.0178*cos2pN-0.0052*cos2N
      CASE(icon_O1); r1(icon) = 0.1887*cosN-0.0065*cos2p-0.0058*cos2N
      CASE(icon_tau1); r1(icon) = -0.2497*cosN+0.0437*cos2p-0.0146*cos2N
      CASE(icon_beta1); r1(icon) = 0.2268*cosN
      CASE(icon_M1)
         r1(icon) = 0.3594*cos2p+0.1722*cosN+0.0664*cos2pN+0.0058*cos2p2N-&
                  & 0.0053*cos2N
      CASE(icon_chi1); r1(icon) = 0.1929*cosN
      CASE(icon_pi1); r1(icon) = -0.0084*cosN
      CASE(icon_P1); r1(icon) = -0.0112*cosN
      CASE(icon_S1); r1(icon) = 0.3529*cos2ps-0.0277*cosN
      CASE(icon_K1); r1(icon) = 0.1159*cosN-0.0029*cos2N
      CASE(icon_psi1); r1(icon) = 0.0171*cosN
      CASE(icon_phi1)
         r1(icon) = -0.0381*cosN+0.0343*cos2p-0.0191*cos2N+0.0133*cos2ps+&
                   & 0.095*cos2p_N
      CASE(icon_theta1); r1(icon) = 0.1671*cosN+0.0304*cos2pN
      CASE(icon_J1)
         r1(icon) = 0.1693*cosN-0.0155*cos2p-0.0097*cos2pN-0.0058*cos2p2N-&
                  & 0.0034*cos2N
      CASE(icon_SO1); r1(icon) = 0.1281*cosN
      CASE(icon_OO1)
         r1(icon) = 0.6404*cosN+0.1497*cos2p+0.1338*cos2N+0.0301*cos2p_N+&
                  & 0.0089*COS(3.0*Nlon)-0.0035*cos2pN
      CASE(icon_KQ1); r1(icon) = 0.6389*cosN+0.1343*cos2N+0.0602*cos2p
      CASE(icon_OQ2); r1(icon) = -0.0389*cosN
      CASE(icon_eps2); r1(icon) = -0.0364*cosN
      CASE(icon_2N2); r1(icon) = -0.0375*cosN-0.0063*cos2p2N
      CASE(icon_mu2); r1(icon) = -0.0373*cosN
      CASE(icon_238); r1(icon) = -0.0385*cosN-0.0308*cosp_ps
      CASE(icon_244); r1(icon) = -0.0294*cosN
      CASE(icon_N2); r1(icon) = -0.0373*cosN-0.0039*cos2p2N 
      CASE(icon_246)
         r1(icon) = -0.5752*cosp_ps-0.1947*COS(2.0*(plon-pslon))-0.0354*cosN
      CASE(icon_nu2); r1(icon) = -0.0374*cosN+0.0044*cos2p-0.0035*cos2pN
      CASE(icon_248); r1(icon) = -0.0377*cosN
      CASE(icon_gamma2); r1(icon) = 0.1474*cos2p2N-0.0368*cosN
      CASE(icon_H1); r1(icon) = -0.0413*cosp_ps-0.0229*cosN
      CASE(icon_M2); r1(icon) = -0.0373*cosN
      CASE(icon_H2); r1(icon) = -0.0207*cosN
      CASE(icon_lambda2); r1(icon) = -0.0451*cosN
      CASE(icon_L2)
         r1(icon) = -0.2503*cos2p-0.1103*cos2pN-0.037*cosN-0.0157*cos2p2N+&
                  & 0.0045*cos2p_N
      CASE(icon_T2); r1(icon) = 0.0
      CASE(icon_S2); r1(icon) = 0.0022*cosN
      CASE(icon_R2); r1(icon) = -0.252*cos2ps+0.0163*COS(Nlon+2.0*pslon)
      CASE(icon_K2); r1(icon) = 0.2582*cosN+0.0324*cos2N
      CASE(icon_eta2); r1(icon) = 0.4161*cosN+0.0492*cos2N-0.0067*cos2p
      CASE(icon_295); r1(icon) = 0.8633*cosN+0.2821*cos2N+0.0427*COS(3.0*Nlon)
      CASE(icon_M3); r1(icon) = 0.0
      END SELECT
   ENDIF
ENDDO icon_410

!---cosine factors for overtides
r1_sh(ish_O1) = 0.1887*cosN-0.0065*cos2p-0.0058*cos2N
r1_sh(ish_K1) = 0.1159*cosN-0.0029*cos2N
r1_sh(ish_N2) = -0.0373*cosN-0.0039*cos2p2N 
r1_sh(ish_M2) = -0.0373*cosN
r1_sh(ish_S2) = 0.0022*cosN
r1_sh(ish_K2) = 0.2582*cosN+0.0324*cos2N

!
!4.2 Sine factors
!----------------
!

icon_420: DO icon=1,notides
   iicon = index_tides(icon)
   IF (iicon.GT.0) THEN
      SELECT CASE (iicon)
      CASE(icon_MS0); r2(icon) = -0.0888*sinN
      CASE(icon_Sa); r2(icon) = -0.0529*sin2ps-0.0183*sinN
      CASE(icon_Ssa)
         r2(icon) = -0.0249*sinN-0.0103*sin2p-0.0055*sin2N-0.0039*sin2ps
      CASE(icon_058); r2(icon) = -0.0166*sinN
      CASE(icon_Msm); r2(icon) = 0.0059*sinN+0.0104*sin2pN+0.003*sin2p2N
      CASE(icon_Mm); r2(icon) = -0.0534*sin2p-0.0219*sin2pN-0.006*sin2p2N
      CASE(icon_Msf); r2(icon) = -0.1372*sinN-0.0069*sin2p
      CASE(icon_Mf)
         r2(icon) = 0.4147*sinN-0.0432*sin2p+0.0387*sin2N+0.0029*sin2p_N+&
                  & 0.0023*sin2pN
      CASE(icon_083)
         r2(icon) = 0.4132*sinN-0.3802*sin2p-0.0372*sin2pN+0.0372*sin2N+&
             & 0.0248*sin2p_N
      CASE(icon_Mt)
         r2(icon) = 0.4146*sinN+0.04*sin2N-0.018*sin2p-0.0039*sin2p2N+&
                  & 0.0031*SIN(2.0*plon-Nlon+pslon)+&
                  & 0.0031*SIN(2.0*plon-Nlon-pslon)-&
                  & 0.0016*SIN(2.0*plon+3.0*Nlon)
      CASE(icon_093); r2(icon) = 0.4118*sinN-0.0539*sin2p+0.0392*sin2N
      CASE(icon_095); r2(icon) = 0.4142*sinN+0.0414*sin2N
      CASE(icon_alpha1); r2(icon) = -0.1907*sinN
      CASE(icon_2Q1); r2(icon) = -0.1883*sinN+0.006*sin2p2N+0.0045*sin2N
      CASE(icon_sigma1); r2(icon) = -0.1883*sinN-0.0087*sin2p+0.005*sin2N
      CASE(icon_Q1)
         r2(icon) = -0.1887*sinN+0.0058*sin2N+0.0038*sin2pN-0.0028*sin2p
      CASE(icon_rho1)
         r2(icon) = -0.1887*sinN-0.0577*sin2p+0.0178*sin2pN+0.0052*sin2N
      CASE(icon_O1); r2(icon) = -0.1887*sinN-0.0065*sin2p+0.0058*sin2N
      CASE(icon_tau1); r2(icon) = -0.1895*sinN-0.0437*sin2p-0.0146*sin2N
      CASE(icon_beta1); r2(icon) = -0.2268*sinN
      CASE(icon_M1)
         r2(icon) = -0.3594*sin2p+0.2294*sinN-0.0664*sin2pN-0.0058*sin2p2N-&
                  & 0.0053*sin2N
      CASE(icon_chi1); r2(icon) = 0.2487*sinN
      CASE(icon_pi1); r2(icon) = 0.0084*sinN
      CASE(icon_P1); r2(icon) = 0.0112*sinN
      CASE(icon_S1); r2(icon) = -0.3529*sin2ps-0.0277*sinN
      CASE(icon_K1); r2(icon) = 0.1555*sinN-0.0029*sin2N
      CASE(icon_psi1); r2(icon) = 0.0171*sinN
      CASE(icon_phi1)
         r2(icon) = -0.0381*sinN-0.0343*sin2p-0.0191*sin2N-0.0133*sin2ps-&
                  & 0.095*sin2p_N
      CASE(icon_theta1); r2(icon) = 0.2279*sinN-0.0304*sin2pN
      CASE(icon_J1)
         r2(icon) = 0.2275*sinN-0.0155*sin2p-0.0097*sin2pN-0.0058*sin2p2N-&
                  & 0.0034*sin2N
      CASE(icon_SO1); r2(icon) = 0.1637*sinN
      CASE(icon_OO1)
         r2(icon) = 0.6404*sinN-0.1497*sin2p+0.1338*sin2N-0.0301*sin2p_N+&
                  & 0.0088*SIN(3.0*Nlon)+0.0035*sin2pN
      CASE(icon_KQ1); r2(icon) = 0.6389*sinN+0.1343*sin2N-0.0602*sin2p
      CASE(icon_OQ2); r2(icon) = 0.0389*sinN
      CASE(icon_eps2); r2(icon) = 0.0364*sinN
      CASE(icon_2N2); r2(icon) = 0.0375*sinN+0.0063*sin2p2N
      CASE(icon_mu2); r2(icon) = 0.0373*sinN
      CASE(icon_238); r2(icon) = 0.0385*sinN+0.0308*sinp_ps
      CASE(icon_244); r2(icon) = 0.0294*sinN
      CASE(icon_N2); r2(icon) = 0.0373*sinN+0.0039*sin2p2N 
      CASE(icon_246)
         r2(icon) = 0.5752*sinp_ps+0.1947*SIN(2.0*(plon-pslon))+0.0354*sinN
      CASE(icon_nu2); r2(icon) = 0.0374*sinN+0.0044*sin2p-0.0035*sin2pN
      CASE(icon_248); r2(icon) = 0.0377*sinN
      CASE(icon_gamma2); r2(icon) = -0.1474*sin2p2N-0.0368*sinN
      CASE(icon_H1); r2(icon) = -0.0413*sinp_ps+0.0229*sinN
      CASE(icon_M2); r2(icon) = 0.0373*sinN
      CASE(icon_H2); r2(icon) = 0.0207*sinN
      CASE(icon_lambda2); r2(icon) = 0.0451*sinN
      CASE(icon_L2)
         r2(icon) = -0.2503*sin2p-0.1103*sin2pN+0.037*sinN-0.0157*sin2p2N+&
                  & 0.0045*sin2p_N
      CASE(icon_T2); r2(icon) = 0.0
      CASE(icon_S2); r2(icon) = -0.0022*sinN
      CASE(icon_R2); r2(icon) = -0.252*sin2ps+0.0163*SIN(Nlon+2.0*pslon)
      CASE(icon_K2); r2(icon) = 0.3108*sinN+0.0324*sin2N
      CASE(icon_eta2); r2(icon) = 0.4563*sinN+0.0492*sin2N-0.0067*sin2p
      CASE(icon_295); r2(icon) = 0.8633*sinN+0.2821*sin2N+0.0427*SIN(3.0*Nlon)
      CASE(icon_M3); r2(icon) = 0.0
      END SELECT
   ENDIF
ENDDO icon_420

!---sine factors for overtides
r2_sh(ish_O1) = -0.1887*sinN-0.0065*sin2p+0.0058*sin2N
r2_sh(ish_K1) = 0.1555*sinN-0.0029*sin2N
r2_sh(ish_N2) = 0.0373*sinN+0.0039*sin2p2N 
r2_sh(ish_M2) = 0.0373*sinN
r2_sh(ish_S2) = -0.0022*sinN
r2_sh(ish_K2) = 0.3108*sinN+0.0324*sin2N

!
!5.3 Nodal corrections
!---------------------
!
!5.3.1 Astronomical constituents
!-------------------------------
!

icon_531: DO icon=1,notides
   iicon = index_tides(icon)
   IF (iicon.GT.0.AND.iicon.LE.MaxAstroTides) THEN
      CALL complex_polar(1.0+r1(icon),r2(icon),xamp=fnode,xpha=phase)
      IF (iopt_astro_pars.EQ.2) fnode_tides(icon) = fnode
      phase_tides(icon) = MOD(phase_tides(icon)+phase,twopi)
   ENDIF
ENDDO icon_531

!---nodal corrections for overtides
ish_532: DO ish=1,6
   CALL complex_polar(1.0+r1_sh(ish),r2_sh(ish),xamp=fnode_sh(ish),xpha=phase)
   phase_sh(ish) = MOD(phase_sh(ish)+phase,twopi)
ENDDO ish_532

!
!5.3.2 Overtides
!---------------
!
!---amplitudes
IF (iopt_astro_pars.EQ.2) THEN
   icon_5321: DO icon=1,notides
      iicon = index_tides(icon)
      IF (iicon.GT.0) THEN
         SELECT CASE (iicon)
         CASE(icon_2SM2); fnode_tides(icon) = fnode_sh(ish_M2)*fnode_sh(ish_S2)
         CASE(icon_2MK3)
            fnode_tides(icon) = fnode_sh(ish_K1)*fnode_sh(ish_M2)**2
         CASE(icon_SO3); fnode_tides(icon) = fnode_sh(ish_O1)*fnode_sh(ish_S2)
         CASE(icon_MK3); fnode_tides(icon) = fnode_sh(ish_M2)*fnode_sh(ish_K1)
         CASE(icon_SK3); fnode_tides(icon) = fnode_sh(ish_K1)*fnode_sh(ish_S2)
         CASE(icon_MN4); fnode_tides(icon) = fnode_sh(ish_M2)*fnode_sh(ish_N2)
         CASE(icon_M4); fnode_tides(icon) = fnode_sh(ish_M2)**2
         CASE(icon_MS4); fnode_tides(icon) = fnode_sh(ish_M2)*fnode_sh(ish_S2)
         CASE(icon_MK4); fnode_tides(icon) = fnode_sh(ish_M2)*fnode_sh(ish_K2)
         CASE(icon_S4); fnode_tides(icon) = fnode_sh(ish_S2)**2
         CASE(icon_2MN6)
            fnode_tides(icon) = fnode_sh(ish_N2)*fnode_sh(ish_M2)**2
         CASE(icon_M6); fnode_tides(icon) = fnode_sh(ish_M2)**3
         CASE(icon_MSN6)
            fnode_tides(icon) = fnode_sh(ish_M2)*fnode_sh(ish_N2)*&
                 & fnode_sh(ish_S2)
         CASE(icon_2MS6)
            fnode_tides(icon) = fnode_sh(ish_S2)*fnode_sh(ish_M2)**2
         CASE(icon_2SM6)
            fnode_tides(icon) = fnode_sh(ish_M2)*fnode_sh(ish_S2)**2
         CASE(icon_S6); fnode_tides(icon) = fnode_sh(ish_S2)**3
         CASE(icon_M8); fnode_tides(icon) = fnode_sh(ish_M2)**4
         CASE(icon_2MSN8)
            fnode_tides(icon) = fnode_sh(ish_N2)*fnode_sh(ish_S2)*&
                              & fnode_sh(ish_M2)**2
         CASE(icon_3MS8)
            fnode_tides(icon) = fnode_sh(ish_S2)*fnode_sh(ish_M2)**3
         CASE(icon_2MS8)
            fnode_tides(icon) = (fnode_sh(ish_S2)*fnode_sh(ish_M2))**2
         CASE(icon_S8); fnode_tides(icon) = fnode_sh(ish_S2)**4
         END SELECT
      ENDIF
   ENDDO icon_5321
ENDIF

!---phases
icon_5322: DO icon=1,notides
   iicon = index_tides(icon)
   IF (iicon.GT.MaxAstroTides) THEN
      SELECT CASE (iicon)
      CASE(icon_2SM2); phase_tides(icon) = phase_sh(ish_M2)
      CASE(icon_2MK3)
         phase_tides(icon) = 2.0*phase_sh(ish_M2) + phase_sh(ish_K1)
      CASE(icon_SO3); phase_tides(icon) = phase_sh(ish_O1)
      CASE(icon_MK3); phase_tides(icon) = phase_sh(ish_M2) + phase_sh(ish_K1)
      CASE(icon_SK3); phase_tides(icon) = phase_sh(ish_K1)
      CASE(icon_MN4); phase_tides(icon) = phase_sh(ish_M2) + phase_sh(ish_N2)
      CASE(icon_M4); phase_tides(icon) = 2.0*phase_sh(ish_M2)
      CASE(icon_MS4); phase_tides(icon) = phase_sh(ish_M2)
      CASE(icon_MK4); phase_tides(icon) = phase_sh(ish_M2) + phase_sh(ish_K2)
      CASE(icon_S4); phase_tides(icon) = 0.0
      CASE(icon_2MN6)
         phase_tides(icon) = 2.0*phase_sh(ish_M2) + phase_sh(ish_N2)
      CASE(icon_M6); phase_tides(icon) = 3.0*phase_sh(ish_M2)
      CASE(icon_MSN6); phase_tides(icon) = phase_sh(ish_M2) + phase_sh(ish_N2)
      CASE(icon_2MS6); phase_tides(icon) = 2.0*phase_sh(ish_M2)
      CASE(icon_2SM6); phase_tides(icon) = phase_sh(ish_M2)
      CASE(icon_S6); phase_tides(icon) = 0.0
      CASE(icon_M8); phase_tides(icon) = 4.0*phase_sh(ish_M2)
      CASE(icon_2MSN8)
         phase_tides(icon) = 2.0*phase_sh(ish_M2) + phase_sh(ish_N2)
      CASE(icon_3MS8); phase_tides(icon) = 3.0*phase_sh(ish_M2)
      CASE(icon_2MS8); phase_tides(icon) = 2.0*phase_sh(ish_M2)
      CASE(icon_S8); phase_tides(icon) = 0.0
      END SELECT
      phase_tides(icon) = MOD(phase_tides(icon),twopi)
   ENDIF
ENDDO icon_5322

!
!6. Phases at reference longitude
!--------------------------------
!

dsecs = 0.001*IGDateTime(7)+IGDateTime(6)+60.0*(IGDateTime(5)+&
                                              & 60.0*IGDateTime(4))
dref = degtorad*dlonref_phase
icon_610: DO icon=1,notides
   iicon = index_tides(icon)
   IF (iicon.GT.0) THEN
      phase_tides(icon) = phase_tides(icon) + dsecs*tidal_spectrum(iicon) + &
                        & ispec_tides(iicon)*dref
      phase_tides(icon) = MOD(phase_tides(icon),twopi)
   ENDIF
ENDDO icon_610

CALL log_timer_out()


RETURN

END SUBROUTINE astro_params

!========================================================================

SUBROUTINE astronomical_tides
!************************************************************************
!
! *astronomical_tides* Astronomical tidal force
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Tidal_Forcing.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - astro_params
!
! Module calls - error_alloc, exchange_mod, mult_index, UVarr_at_U, UVarr_at_V
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE modids
USE switches
USE syspars
USE tide
USE timepars
USE array_interp, ONLY: UVarr_at_U, UVarr_at_V
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: icon, iicon, npcc
INTEGER, DIMENSION(0:3) :: nospec
REAL, DIMENSION(nconastro) :: ampl_tides
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu, maskatv
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:) :: phase_astro_d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d1, array2d2, gxlonatu, &
                                         & gxlonatv, gylatatu, gylatatv
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:,:) :: phase


procname(pglev+1) = 'astronomical_tides'
CALL log_timer_in(npcc)

!
!1. Allocate/initialise arrays
!-----------------------------
!

ALLOCATE (maskatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (maskatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (phase(ncloc,nrloc),STAT=errstat)
CALL error_alloc('phase',2,(/ncloc,nrloc/),kndrlong)
ALLOCATE (array2d1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2d2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d2',2,(/ncloc,nrloc/),kndrtype)

maskatu = node2du(1:ncloc,1:nrloc).GT.1
maskatv = node2dv(1:ncloc,1:nrloc).GT.1

!
!2. Initialise on first call
!---------------------------
!

IF (nt.EQ.1) THEN

!
!2.1 Allocate arrays
!-------------------
!

   ALLOCATE (gxlonatu(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('gxlonatu',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (gxlonatv(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('gxlonatv',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (gylatatu(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('gylatatu',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (gylatatv(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('gylatatv',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE(phase_astro_d(nconastro),STAT=errstat)
   CALL error_alloc('phase_astro_d',1,(/nconastro/),kndrlong)
   IF (nconastro.GT.0) phase_astro_d = phase_astro

!
!2.2 Geographical coordinates at U- and V-nodes
!----------------------------------------------
!
!  ---longitude
   CALL UVarr_at_U(gxcoord(1:ncloc,1:nrloc+1),gxlonatu,0,0,(/1,1,nz/),&
                & (/ncloc,nrloc+1,nz/),1,0,.TRUE.)
   gxlonatu = degtorad*gxlonatu
   CALL UVarr_at_V(gxcoord(1:ncloc+1,1:nrloc),gxlonatv,0,0,(/1,1,nz/),&
                & (/ncloc+1,nrloc,nz/),1,0,.TRUE.)
   gxlonatv = degtorad*gxlonatv

!  ---latitude
   CALL UVarr_at_U(gycoord(1:ncloc,1:nrloc+1),gylatatu,0,0,(/1,1,nz/),&
                & (/ncloc,nrloc+1,nz/),1,0,.TRUE.)
   gylatatu = degtorad*gylatatu
   CALL UVarr_at_V(gycoord(1:ncloc+1,1:nrloc),gylatatv,0,0,(/1,1,nz/),&
                & (/ncloc+1,nrloc,nz/),1,0,.TRUE.)
   gylatatv = degtorad*gylatatv

ENDIF

!
!2. Amplitudes and phases
!------------------------
!

IF (nt.EQ.1.OR.mult_index(nt,icnodal)) THEN
   CALL astro_params(index_astro(1:nconastro),fnode_astro,phase_astro,&
                   & nconastro,0.0,CDateTime)
   phase_astro_d = phase_astro
ELSE
   icon_212: DO icon=1,nconastro
      iicon = index_astro(icon)
      phase_astro_d(icon) = MOD(phase_astro_d(icon)+&
                              & delt2d*tidal_spectrum(iicon),twopi_d)
   ENDDO icon_212
   phase_astro = phase_astro_d
ENDIF

nospec = 0
icon_211: DO icon=1,nconastro
   iicon = index_astro(icon)
   ampl_tides(icon) = astro_earth(iicon)*astro_ampl(iicon)*fnode_astro(icon)
   nospec(ispec_tides(iicon)) = nospec(ispec_tides(iicon)) + 1
ENDDO icon_211

!
!3. X-component of tidal force
!-----------------------------
!

fxastro = 0.0

!
!3.1 Long-period terms
!---------------------
!

IF (nospec(0).GT.0) THEN
   array2d2 = 0.0
   icon_310: DO icon=1,nconastro
      iicon = index_astro(icon)
      IF (ispec_tides(iicon).EQ.0) THEN
         phase = MOD(phase_astro_d(icon),twopi_d)
         array2d2 = array2d2 + ampl_tides(icon)*COS(phase)
      ENDIF
   ENDDO icon_310
   WHERE (maskatu)
      fxastro(1:ncloc,:) = fxastro(1:ncloc,:) - &
              & 1.5*gaccatu(1:ncloc,:)*SIN(2.0*gylatatu)*array2d2*&
              & (gylat(1:ncloc,1:nrloc)-gylat(0:ncloc-1,1:nrloc))&
              & /delxatu(1:ncloc,1:nrloc)
   END WHERE
ENDIF

!
!3.2 Diurnal terms
!-----------------
!

IF (nospec(1).GT.0) THEN
   array2d1 = 0.0; array2d2 = 0.0
   icon_320: DO icon=1,nconastro
      iicon = index_astro(icon)
      IF (ispec_tides(iicon).EQ.1) THEN
         phase = MOD(gxlonatu+phase_astro_d(icon),twopi_d)
         array2d1 = array2d1 + ampl_tides(icon)*SIN(phase)
         array2d2 = array2d2 + ampl_tides(icon)*COS(phase)
      ENDIF
   ENDDO icon_320
   WHERE (maskatu)
      fxastro(1:ncloc,:) = fxastro(1:ncloc,:) - gaccatu(1:ncloc,:)*&
              & (SIN(2.0*gylatatu)*array2d1*&
              & (gxlon(1:ncloc,1:nrloc)-gxlon(0:ncloc-1,1:nrloc))-&
              & 2.0*COS(2.0*gylatatu)*array2d2*&
              & (gylat(1:ncloc,1:nrloc)-gylat(0:ncloc-1,1:nrloc)))&
              & /delxatu(1:ncloc,1:nrloc)
   END WHERE
ENDIF

!
!3.3 Semi-diurnal terms
!----------------------
!

IF (nospec(2).GT.0) THEN
   array2d1 = 0.0; array2d2 = 0.0
   icon_330: DO icon=1,nconastro
      iicon = index_astro(icon)
      IF (ispec_tides(iicon).EQ.2) THEN
         phase = MOD(2.0*gxlonatu+phase_astro_d(icon),twopi_d)
         array2d1 = array2d1 + ampl_tides(icon)*SIN(phase)
         array2d2 = array2d2 + ampl_tides(icon)*COS(phase)
      ENDIF
   ENDDO icon_330
   WHERE (maskatu)
      fxastro(1:ncloc,:) = fxastro(1:ncloc,:) - gaccatu(1:ncloc,:)*&
              & ((1+COS(2.0*gylatatu))*array2d1*&
              & (gxlon(1:ncloc,1:nrloc)-gxlon(0:ncloc-1,1:nrloc))+&
              & SIN(2.0*gylatatu)*array2d2*&
              & (gylat(1:ncloc,1:nrloc)-gylat(0:ncloc-1,1:nrloc)))&
              & /delxatu(1:ncloc,1:nrloc)
   END WHERE
ENDIF

!
!3.4 Ter-diurnal terms
!---------------------
!

IF (nospec(3).GT.0) THEN
   array2d1 = 0.0; array2d2 = 0.0
   icon_340: DO icon=1,nconastro
      iicon = index_astro(icon)
      IF (ispec_tides(iicon).EQ.3) THEN
         phase = MOD(3.0*gxlonatu+phase_astro_d(icon),twopi_d)
         array2d1 = array2d1 + ampl_tides(icon)*SIN(phase)
         array2d2 = array2d2 + ampl_tides(icon)*COS(phase)
      ENDIF
   ENDDO icon_340
   WHERE (maskatu)
      fxastro(1:ncloc,:) = fxastro(1:ncloc,:) - 0.75*gaccatu(1:ncloc,:)*&
              & ((3.0*COS(gylatatu)+COS(3.0*gylatatu))*array2d1*&
              & (gxlon(1:ncloc,1:nrloc)-gxlon(0:ncloc-1,1:nrloc))+&
              & (SIN(gylatatu)+SIN(3.0*gylatatu))*array2d2*&
              & (gylat(1:ncloc,1:nrloc)-gylat(0:ncloc-1,1:nrloc)))&
              & /delxatu(1:ncloc,1:nrloc)
   END WHERE
ENDIF

!
!4. Y-component of tidal force
!-----------------------------
!

fyastro = 0.0

!
!4.1 Long-period terms
!---------------------
!

IF (nospec(0).GT.0) THEN
   array2d2 = 0.0
   icon_410: DO icon=1,nconastro
      iicon = index_astro(icon)
      IF (ispec_tides(iicon).EQ.0) THEN
         phase = MOD(phase_astro_d(icon),twopi_d)
         array2d2 = array2d2 + ampl_tides(icon)*COS(phase_astro_d(icon))
      ENDIF
   ENDDO icon_410
   WHERE (maskatv)
      fyastro(:,1:nrloc) = fyastro(:,1:nrloc) - &
              & 1.5*gaccatv(:,1:nrloc)*SIN(2.0*gylatatv)*array2d2*&
              & (gylat(1:ncloc,1:nrloc)-gylat(1:ncloc,0:nrloc-1))&
              & /delyatv(1:ncloc,1:nrloc)
   END WHERE
ENDIF

!
!4.2 Diurnal terms
!-----------------
!

IF (nospec(1).GT.0) THEN
   array2d1 = 0.0; array2d2 = 0.0
   icon_420: DO icon=1,nconastro
      iicon = index_astro(icon)
      IF (ispec_tides(iicon).EQ.1) THEN
         phase = MOD(gxlonatv+phase_astro_d(icon),twopi_d)
         array2d1 = array2d1 + ampl_tides(icon)*SIN(phase)
         array2d2 = array2d2 + ampl_tides(icon)*COS(phase)
      ENDIF
   ENDDO icon_420
   WHERE (maskatv)
      fyastro(:,1:nrloc) = fyastro(:,1:nrloc) - gaccatv(:,1:nrloc)*&
              & (SIN(2.0*gylatatv)*array2d1*&
              & (gxlon(1:ncloc,1:nrloc)-gxlon(1:ncloc,0:nrloc-1))-&
              & 2.0*COS(2.0*gylatatv)*array2d2*&
              & (gylat(1:ncloc,1:nrloc)-gylat(1:ncloc,0:nrloc-1)))&
              & /delyatv(1:ncloc,1:nrloc)
   END WHERE
ENDIF

!
!4.3 Semi-diurnal terms
!----------------------
!

IF (nospec(2).GT.0) THEN
   array2d1 = 0.0; array2d2 = 0.0
   icon_430: DO icon=1,nconastro
      iicon = index_astro(icon)
      IF (ispec_tides(iicon).EQ.2) THEN
         phase = MOD(2.0*gxlonatv+phase_astro_d(icon),twopi_d)
         array2d1 = array2d1 + ampl_tides(icon)*SIN(phase)
         array2d2 = array2d2 + ampl_tides(icon)*COS(phase)
      ENDIF
   ENDDO icon_430
   WHERE (maskatv)
      fyastro(:,1:nrloc) = fyastro(:,1:nrloc) - gaccatv(:,1:nrloc)*&
              & ((1+COS(2.0*gylatatv))*array2d1*&
              & (gxlon(1:ncloc,1:nrloc)-gxlon(1:ncloc,0:nrloc-1))+&
              & SIN(2.0*gylatatv)*array2d2*&
              & (gylat(1:ncloc,1:nrloc)-gylat(1:ncloc,0:nrloc-1)))&
              & /delyatv(1:ncloc,1:nrloc) 
   END WHERE
ENDIF

!
!4.4 Ter-diurnal terms
!---------------------
!

IF (nospec(3).GT.0) THEN
   array2d1 = 0.0; array2d2 = 0.0
   icon_440: DO icon=1,nconastro
      iicon = index_astro(icon)
      IF (ispec_tides(iicon).EQ.3) THEN
         phase = MOD(3.0*gxlonatv+phase_astro_d(icon),twopi_d)
         array2d1 = array2d1 + ampl_tides(icon)*SIN(phase)
         array2d2 = array2d2 + ampl_tides(icon)*COS(phase)
      ENDIF
   ENDDO icon_440
   WHERE (maskatv)
      fyastro(:,1:nrloc) = fyastro(:,1:nrloc) - 0.75*gaccatv(:,1:nrloc)*&
              & ((3.0*COS(gylatatv)+COS(3.0*gylatatv))*array2d1*&
              & (gxlon(1:ncloc,1:nrloc)-gxlon(1:ncloc,0:nrloc-1))+&
              & (SIN(gylatatv)+SIN(3.0*gylatatv))*array2d2*&
              & (gylat(1:ncloc,1:nrloc)-gylat(1:ncloc,0:nrloc-1)))&
              & /delyatv(1:ncloc,1:nrloc)
   END WHERE
ENDIF

!
!5. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL exchange_mod(fxastro,(/1,1/),(/0,1,0,0/),iarr_fxastro,corners=.FALSE.)
   CALL exchange_mod(fyastro,(/1,1/),(/0,0,0,1/),iarr_fyastro,corners=.FALSE.)
ENDIF

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatu,maskatv,phase,array2d1,array2d2)
IF (cold_start.OR.(iopt_tidal_accel.EQ.0.AND.nt.EQ.nstep).OR.&
 & (iopt_tidal_accel.EQ.1.AND.nt.EQ.hydro_end)) THEN
   DEALLOCATE (gxlonatu,gxlonatv,gylatatu,gylatatv,phase_astro_d)
ENDIF

CALL log_timer_out(npcc,itm_astro)


RETURN

END SUBROUTINE astronomical_tides
