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

MODULE physpars
!************************************************************************
!
! *physpars* Physical and numerical model parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)physpars.f90  V2.11.2
!
! $Date: 2018-07-23 10:21:02 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1167 $
!
! Description - 
!
!************************************************************************
!
USE iopars
USE syspars

IMPLICIT NONE


!---general
REAL, PARAMETER :: celtokelv = 273.15, kboltz = 1.38065E-23, &
                 & Rearth = 6371000.0, stefbolz = 5.67E-08
REAL :: rho_air = 1.2, specheat = 3987.5

!---reference and minimum values
REAL :: atmpres_ref = 101325.0, beta_sal_ref, beta_temp_ref, density_ref, &
      & dlat_ref = 0.0, dlon_ref = 0.0, dlon_ref_anal = 0.0, &
      & dlon_ref_obc = 0.0, gacc_mean, gacc_ref, sal_ref = 33.0, &
      & sst_ref = 12.0, temp_min = 0.0, temp_ref = 12.0

!---model grid
REAL :: b_SH = 0.1, dl_BB = 1.5, du_BB = 1.5, hcrit_SH = 200.0, &
      & sigstar_DJ = 0.0, sig0_DJ = 0.1, theta_SH = 8.0

!---diffusion coefficients
REAL :: hdifmom_cst = 0.0, hdifmom_fac, hdifscal_cst = 0.0, hdifscal_fac, &
      & kinvisc_cst = 1.0E-06, rho_crit_iso = 0.01, skewdiff_cst = 0.0, &
      & slopemax_iso = 0.01, smag_coef_mom = 0.1, smag_coef_scal = 0.1, &
      & vdifmom_cst = 1.0E-06, vdifscal_cst = 1.0E-06

!---water depths
REAL :: depmean_cst = 0.0, depmean_flag = 0.0, depmean_flag_eps

!---inundation schemes
INTEGER, PARAMETER :: nofldmasks = 11
REAL :: dcrit_fld = 0.1, dmin_fld = 0.02, dthd_fld = 0.1, ndwexp = 1.0
INTEGER, DIMENSION(nofldmasks) :: fld_mask

!---bottom/surface fluxes
REAL :: bdragcoef_cst = 0.0, bdraglin = 0.0, ccharno = 0.014, &
      & ces_scst = 0.00113, ces_ucst = 0.00113, chs_scst = 0.00113, &
      & chs_ucst = 0.00113, ckar = 0.4, zbtoz0lim = 2.0, zref_temp = 10.0, &
      & zref_wind = 10.0, zrough_cst = 0.0
REAL, DIMENSION(4) :: cdspars = (/0.0,0.00113,0.0,0.0/)

!---relaxation distance for momentum advection
REAL :: distrlx_obc = 0.0

!---open boundary conditions
REAL :: cgravratio = 0.03

!---optical parameters
REAL :: optattcoef1_cst = 10.0, optattcoef2_cst = 0.067, opt_frac = 0.54

!---implicit code
INTEGER :: itsimp = 0, maxitsimp = 1, maxmgiter, minmgiter, mgiteration, &
         & noitsimp, nomgiterations = 100, nomglevels = 1, nopostsweeps = 2, &
         & nopresweeps = 2, nosmoothsteps = 2
REAL :: dzetaresid, dzetaresid_conv = 1.0E-5, meanmgiter, mg_tol = 1.0E-07, &
      & residual_norm, ur_mg = 1.0, ur_smooth = 0.8

!---implicity factors
REAL :: theta_cor = 0.5, theta_sur = 0.5, theta_vadv = 0.501, theta_vdif = 1.0

!---numerics (other)
REAL, PARAMETER :: eps_adv = 1.0E-12

!---wave parameters
REAL :: wave_penetration_bed = 10.0, wave_penetration_surf = 2.0, &
      & wavethick_cst = 0.1

SAVE

!
! Name              Type    Purpose
!------------------------------------------------------------------------------
!*atmpres_ref*      REAL    Reference atmospheric pressure                  [Pa]
!*bdragcoef_cst*    REAL    Uniform bottom drag coefficient (iopt_bstres_drag=0)
!*bdraglin*         REAL    Linear bottom drag coefficient (iopt_bstres_form=1)
!                                                                          [m/s]
!*beta_sal_ref*     REAL    Reference value for salinity expansion coefficient
!                                                                       [PSU^-1]
!*beta_temp_ref*    REAL    Reference value for temperature expansion
!                           coefficient                                [1/deg C]
!*b_SH*             REAL    Parameter b in Song and Haidvogel (1994) vertical
!                           grid transformation
!*ccharno*          REAL    Charnock's constant
!*cdspars*          REAL    Parameters for constructing generic drag law as
!                           function of wind speed (U)/
!                             cds = cdspars(1) if U<U_crit
!                             cds = cdspars(2)+U*cdspars(3) if U>=Ucrit
!                             U_crit = cdspars(4)
!*celtokelv*      REAL    Factor to conver deg C to deg K                [deg C]
!*ces_scst*         REAL    Uniform exchange coefficient for latent heat in
!                           case of stable conditions (iopt_sflux_spars=1,4,7)
!*ces_ucst*         REAL    Uniform exchange coefficient for latent heat in
!                           case of unstable conditions (iopt_sflux_spars=1,4,7)
!*chs_scst*          REAL   Uniform exchange coefficient for sensible heat in
!                           case of stable conditions (iopt_sflux_spars=1,4,7)
!*chs_ucst*          REAL   Uniform exchange coefficient for sensible heat in
!                           case of unstable conditions (iopt_sflux_spars=1,4,7)
!*cgravratio*       REAL    Ratio of internal to external wave speed (at o.b.)
!*ckar*             REAL    Von Karman's constant
!*dcrit_fld*        REAL    Critical water depth for flooding/drying scheme
!*density_ref*      REAL    Reference density                           [kg/m^3]
!*distrlx_obc*      REAL    Maximum distance (from open boundaries) for
!                           application of relaxation of momentum advection  [m]
!*dlat_ref*         REAL    Reference latitude in case of Cartesian grid
!                           (decimal degrees, positive North)
!*dlon_ref*         REAL    Reference longitude in case of Cartesian grid
!                           (decimal degrees, positive East)
!*dlon_ref_anal*    REAL    If iopt_astro_pars>0, harmonically analysed phases
!                           are taken with respect to astronomical argument for
!                           this reference longitude at the central time
!                           (decimal degrees, positive East).
!*dlon_ref_obc*     REAL    If iopt_astro_pars>0, phases at open boundaries are
!                           assumed to be taken with respect to astronomical
!                           argument at this reference value (decimal degrees,
!                           positive East).
!*dl_BB*            REAL    Parameter dl in Burchard and Bolding (2002) vertical
!                           grid transformation
!*depmean_cst*      REAL    Mean water depth in case a uniform depth is taken
!*depmean_flag*   REAL    Data flag for marking land points in the bathymetry
!                                                                            [m]
!*dmin_fld*         REAL    Minimum water depth for flooding/drying scheme
!*dth_fld*          REAL    Threshold water depth for flooding/drying scheme
!*du_BB*            REAL    Parameter du in Burchard and Bolding (2002) vertical
!                           grid transformation
!*dzetaresid*       REAL    Maximum (over the whole domain) residual of surface
!                           elevation after the last iteration in the outer loop
!                           of the implicit scheme
!*dzetaresid_conv*REAL      Allowed  maximum (over the whole domain) residual of
!                           surface elevation after the last iteration in the
!                           outer loop of the implicit scheme
!*eps_adv*          REAL    Tolerance factor for calculation of flux ratios
!                           (TVD scheme)
!*fld_mask*         INTEGER Enables (0) or disables (1) a specific mask
!                           criterium
!*gacc_mean*        REAL    Domain-averaged acceleration of gravity      [m/s^2]
!*gacc_ref*         REAL    Uniform reference acceleration of gravity (if
!                           flagged, gravity acceleration is computed as
!                           function of latitude)                        [m/s^2]
!*hcrit_SH*         REAL    Critical water depth in Song and Haidvogel (1994)
!                           vertical grid transformation                     [m]
!*hdifmom_cst*      REAL    Uniform horizontal diffusion coefficient for
!                           momentum
!*hdifmom_fac*      REAL    Factor used for defining the horizontal diffusion
!                           coefficient for momentum 
!*hdifscal_cst*     REAL    Uniform horizontal diffusion coefficient for scalars
!                                                                        [m^2/s]
!*hdifscal_fac*     REAL    Factor used for defining the horizontal diffusion
!                           coefficient for scalars 
!*itsimp*           INTEGER Outer iteration number of the implicit scheme
!*kboltz*         REAL    Boltzmann's constant            [(m^2 kg)/(s^2 deg K)]
!                                                                   [J/kg/deg C]
!*kinvisc_cst*      REAL    Reference value for the kinematic viscosity  [m^2/s]
!*maxitsimp*        INTEGER Maximum allowed number of outer iterations for the
!                           implicit scheme
!*mgiteration*      INTEGER Iteration index in the multigrid iterative algorithm
!*mg_tol*           REAL    Relative tolerance used by multigrid for solving
!                           the  linear system
!*ndwexp*           REAL    Exponent in the power law used to define the
!                           alpha factor in the inundation scheme
!*nofldmasks*       INTEGER Number of available mask criteria
!*noitsimp*         INTEGER Number of outer iterations for the implicit scheme
!*optattcoef1_cst*  REAL    Attenuation coefficient for infrared radiation [1/m]
!*optattcoef2_cst*  REAL    Attenuation coefficient for short-wave radiation and
!                           pure sea water                                 [1/m]
!*opt_frac*         REAL    Infrared fraction of irradiance absorbed at sea
!                           surface
!*maxmgiter*        INTEGER Maximum number of multigrid iterations for the whole
!                           simulation
!*meanmgiter        REAL    Mean number of multigrid iterations over the whole
!                           simulation
!*minmgiter*        INTEGER Minimum number of multigrid iterations for the whole
!                           simulation
!*nomgiterations*   INTEGER Maximum number of iterations in the multigrid
!                           procedure
!*nomglevels*       INTEGER Number of multi-grid levels (including main grid)
!*nopostsweeps*     INTEGER Number of postrelaxation iterations
!*nopresweeps*      INTEGER Number of prerelaxation iterations
!*nosmoothsteps*    INTEGER Number of linear solver (smoother) iterations at
!                           the  coarsest level
!*Rearth*           REAL    Mean radius of the Earth                         [m]
!*residual_norm*    REAL    Residual norm used for checking the convergence of
!                           the multigrid algorithm
!*rho_air*          REAL    Air density                                 [kg/m^3]
!*rho_crit_iso*     REAL   Defines for the isoneutral diffusion scheme the
!                           depth of pycnoline as the deepest point where the
!                           density is higher that this critical threshold
!                                                                        [kg/m3]
!*sal_ref*          REAL    Reference salinity                             [PSU]
!*skewdiff_cst*     REAL    Coefficient used in the skew-diffusive formulation
!                           for the Gent and McWilliams (1990) eddy induced
!                           transport velocity                           [m^2/s]
!*slopemax_iso*     REAL    Maximum threshold for the density slopes in the
!                           isoneutral formulation
!*sigstar_DJ*       REAL    Parameter sigma* in the Davies and Jones (1991)
!                           vertical grid transformation
!*sig0_DJ*          REAL    Parameter sigma_0 in the Davies and Jones (1991)
!                           vertical grid transformation
!*smag_coef_mom*    REAL    Smagorinsky coefficient for momentum 
!*smag_coef_scal*   REAL    Smagorinsky coefficient for scalars
!*specheat*         REAL    Specific heat of seawater at constant pressure
!*sst_ref*          REAL    Reference sea surface temperature            [deg C]
!*stefboltz*        REAL    Stefan-Boltzmann constant [W/m^2/(deg K)^4
!*temp_min*         REAL    Minimum allowed value for temperature (if flagged,
!                           minimum temperature is set to freezing point
!                           temperature)                                 [deg C]
!*temp_ref*         REAL    Reference water temperature                  [deg C]
!*theta_cor*        REAL    Implicity factor for Coriolis terms
!*theta_sur*        REAL    Implicity factor for the surface slope
!*theta_SH*         REAL    Parameter theta in Song and Haidvogel (1994)
!                           vertical grid transformation
!*theta_vadv*       REAL    Implicity factor for vertical advection
!*theta_vdif*       REAL    Implicity factor for vertical diffusion
!*ur_mg*            EAL     Underrelaxation factor for multigrid updates
!*ur_smooth*        REAL    Underrelaxation factor for smoother
!*vdifmom_cst*      REAL    Uniform/background vertical diffusion coefficient
!                           for momentum                                 [m^2/s]
!*vdifscal_cst*     REAL    Uniform/background vertical diffusion coefficient
!                           for scalars                                  [m^2/s]
!*wave_penetration_bed* REAL Penetration length divided by wave boundary layer
!                            thickness for the 3D force profile due to bed
!                            dissipation
!*wave_penetration_surf*REAL Penetration length divided by ave lenght for the 3D
!                            force profile due to surface dissipation [-]
!*wavethick_cst*     REAL   Uniform thickness of the wave bottom boundary layer
!                           if iopt_waves_wfric = 0                          [m]
!*zbtoz0lim*         REAL   Minimum ratio for botttom cell height to bottom
!                           roughness
!*zref_temp*         REAL   Reference height for air temperature             [m]
!*zref_wind*         REAL   Reference height for surface wind                [m]
!*zrough_cst*        REAL   Uniform bottom roughness length (iopt_bstres_drag=2)

!
!************************************************************************
!

END MODULE physpars
