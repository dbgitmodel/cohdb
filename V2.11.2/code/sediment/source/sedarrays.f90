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

MODULE sedarrays
!************************************************************************
!
! *sedarrays* Sediment arrays
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)sediment.f90  V2.11.1
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
!************************************************************************
!
!
 
IMPLICIT NONE

!---particle attributes
REAL, ALLOCATABLE, DIMENSION(:) :: bstres_cr_cst, dp, massp, rhos, volp, ws_cst

!---sediment concentrations
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ceq, cref, ctot
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) ::  dar_sediment, cvol

!---flocculation module
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: floc_dens, floc_dia, floc_nc, floc_vol
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: cnump

!---sediment loads
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: qbedatu, qbedatu_int, qbedatv, &
     & qbedatv_int, qsusatc, qsusatu_int, qsusatv_int, qtotatu, qtotatu_int, &
     & qtotatv, qtotatv_int
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) ::  qsusatu, qsusatv

!---bottom fluxes
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: bottom_sed_dep, bottom_sed_ero, &
     & bottom_sed_flux, height_c, t_equil
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) ::  bed_fraction

!---bottom stress related arrays
REAL, ALLOCATABLE, DIMENSION(:,:) :: bdragcoefatc_sed, bdragcoefatu_sed, &
     & bdragcoefatv_sed, bstresatc_sed, bstresatc_mean_sed, bstresatc_wav_sed, &
     & bstresatu_sed, bstresatv_sed, d50_bed, fwave_sed, ubstresatu_sed, &
     & vbstresatv_sed, wavethickatc_sed, zaroughatc_sed, zroughatc_sed, &
     & zroughatu_sed, zroughatv_sed
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: bstres_cr

!---diffusion coefficients
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: beta_sed
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) ::  vdiffcoef_sed

!---bed slopes
REAL, ALLOCATABLE, DIMENSION(:,:) :: bed_slope_x, bed_slope_x_atu, &
     & bed_slope_x_atv, bed_slope_y, bed_slope_y_atu, bed_slope_y_atv

!---density
REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  densw
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) ::  beta_state_sed

!---fall velocity
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) ::  wfall

!---open boundary conditions
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: obccnumpatu, obccnumpatv, obcsedatu, &
                                       & obcsedatv

!
! Name              Type Purpose
!------------------------------------------------------------------------------
!
!*bdragcoefatc_sed*  REAL Bottom drag coefficient at C-nodes
!*bdragcoefatu_sed*  REAL Bottom drag coefficient at U-nodes
!*bdragcoefatv_sed*  REAL Bottom drag coefficient at V-nodes
!*bed_fraction*      REAL Amount of material in the bed for each sediment
!                         fraction                                           [-]
!*bed_slope_x*       REAL Slope of the bed in x direction at C-nodes       [m/m]
!*bed_slope_x_atu*   REAL Slope of the bed in x direction at U-nodes       [m/m]
!*bed_slope_x_atv*   REAL Slope of the bed in x direction at the V-nodes   [m/m]
!*bed_slope_y*       REAL Slope of the bed in y direction at the C-nodes   [m/m]
!*bed_slope_y_atu*   REAL Slope of the bed in y direction at the U-nodes   [m/m]
!*bed_slope_y_atv*   REAL Slope of the bed in y direction at the V-nodes   [m/m]
!*beta_sed *         REAL Ratio between eddy diffusivity and eddy viscosity for
!                         sediment                                           [-]
!*beta_state_sed*    REAL Expansion factor in equation of state (d rho/dC)
!                         (per sediment fraction)                            [-]
!*bottom_sed_dep*    REAL Amount of deposited sediment per unit area and time
!                                                                    [m^3/s/m^2]
!*bottom_sed_ero*    REAL Amount of eroded sediment per unit area and time
!                                                                    [m^3/s/m^2]
!*bottom_sed_flux*   REAL Amount of material suspended from the bed to the flow
!                        (per fraction)                              [m^3/s/m^2]
!*bstres_cr*         REAL Critical bed shear stress (per sediment fraction)
!                                                                     [m^2/s^2]]
!*bstres_cr_cst*     REAL Constant particle critical shear stress per fraction
!                                                                      [m^2/s^2]
!*bstresatc_sed*     REAL Maximum current-wave bottom stress at C-nodes
!                                                                      [m^2/s^2]
!*bstresatc_mean_sed*REAL Mean (current-wave) bottom stress at C-nodes [m^2/s^2]
!*bstresatc_waw_sed* REAL Bottom stress (wave only) at C-nodes         [m^2/s^2]
!*bstresatu_sed*     REAL (Maximum) bottom stress at U-nodes           [m^2/s^2]
!*bstresatv_sed*     REAL (Maximum) bottom stress at V-nodes           [m^2/s^2]
!*ceq*               REAL Equilibrium sediment concentration in 2D (per
!                         fraction)                                    [m^3/m^3]
!*cnump*             REAL Number concentration for the state variables in the
!                         flocculation model (microflocs, macroflocs,
!                         total number of microflocs bound in macroflocs)[1/m^3]
!*cref*              REAL Near bed reference (at height_c) sediment
!                         concentration (per fraction)                 [m^3/m^3]
!*ctot*              REAL Volumetric sediment concentration (sum over all
!                        fractions)                                    [m^3/m^3]
!*cvol*              REAL Volumetric sediment concentration            [m^3/m^3]
!                         Meaning depends on the value of iopt_sed
!                      = 1 => Volumetric sediment concentration (per sediment
!                             fraction)                                   
!                      = 2 => First fraction represents concentration of
!                             microflocs, second fraction concentration of
!                             macroflocs 
!*dar_sediment*      REAL Sediment sources due to dredging and relocation
!                         (per fraction)                             [m^3/m^3/s]
!*densw*             REAL Water density without sediment                [kg/m^3]
!*dp*                REAL Particle diameter (per sediment fraction)          [m]
!*d50_bed*           REAL Median grain size (by mass)  at the bed            [m]
!*floc_dens*         REAL Macrofloc density                             [kg/m^3]
!*floc_dia*          REAL Macrofloc diameter                                 [m]
!*floc_nc*           REAL Number of macroflocs bound in a nacrofloc          [-]
!*floc_tot*          REAL Volumetric concentration of microflocs bound in
!                         macroflocs                                   [m^3/m^3]
!*floc_vol*          REAL Macrofloc volume                             [m^3/m^3]
!*fwave_sed*         REAL Wave friction factor
!*height_c*          REAL The elevation at which the near bed boundary condition
!                         is applied (per fraction)                          [m]
!*massp*             REAL Particle mass (per sediment fraction)         [kg/m^3]
!*obccnumpatu*       REAL Storage array for sediment radiation o.b.c. at the
!                         U-nodes (flocculation model)
!*obccnumpatv*       REAL Storage array for sediment radiation o.b.c. at the
!                         V-nodes (flocculation model)
!*obcsedatu*         REAL Storage array for sediment radiation o.b.c. at the
!                         U-nodes
!*obcsedatv*         REAL Storage array for sediment radiation o.b.c. at the
!                         V-nodes
!*qbedatu*           REAL Bed load (per sediment fraction) in the X-direction
!                         at the U-nodes                                 [m^2/s]
!*qbedatu_int*       REAL Bed load (per sediment fraction) in the X-direction,
!                         averaged over a morphological time step, at the
!                         U-nodes                                        [m^2/s]
!*qbedatv*           REAL Bed load (per sediment fraction) in the Y-direction
!                         at the V-nodes                                 [m^2/s]
!*qbedatv_int*       REAL Bed load (per sediment fraction) in the X-direction,
!                         averaged over a morphological time step, at the
!                         U-nodes                                        [m^2/s]
!*qsusatc*           REAL Depth-integrated suspended load (per sediment
!                         fraction) at the C-nodes                       [m^2/s]
!*qsusatu*           REAL Suspended load (per sediment fraction) at U-nodes
!                                                                         [m/s]
!*qsusatu_int*       REAL Suspended load (per sediment fraction) in the
!                         X-direction, averaged over a morphological time and
!                         integrated over the water column, step at the U-nodes
!                                                                        [m^2/s]
!*qsusatv*           REAL Suspended load (per sediment fraction) at the
!                         V-nodes                                          [m/s]
!*qsusatv_int*       REAL Suspended load (per sediment fraction) in the
!                         Y-direction, averaged over a morphological time and
!                         integrated over the water column, step at the V-nodes
!                                                                        [m^2/s]
!*qtotatu*           REAL Total load (per sediment fraction) in the X-direction
!                         at the U-nodes                                 [m^2/s]
!*qtotatu_int*       REAL Total load (per sediment fraction) in the Y-direction,
!                         averaged over a morphological time step at the
!                         U-nodes                                        [m^2/s]
!*qtotatv*           REAL Total load (per sediment fraction) in the Y-direction
!                         at the V-nodes                                 [m^2/s]
!*qtotatv_int*       REAL Total load (per sediment fraction) in the Y-direction,
!                         averaged over a morphological time step at the
!                         V-nodes                                        [m^2/s]
!*rhos*              REAL Density of the sediment particles (per sediment
!                         fraction)                                     [kg/m^3]
!*t_equil*           REAL Time scale in which sediment concentrations in 2-D
!                         go to equilibrium (per fraction)                   [-]
!*ubstresatu_sed*    REAL X-component of the (maximum) bottom stress at the
!                         U-node                                       [m^2/s^2]
!*vbstresatv_sed*    REAL Y-component of the (maximum) bottom stress at the
!                         V-nodes                                      [m^2/s^2]
!*vdiffcoef_sed*     REAL Vertical eddy diffusivity for sediment (per fraction)
!                                                                        [m^2/s]
!*volp*              REAL Particle volume (per sediment fraction)          [m^3]
!*wavethickatc_sed*  REAL Thickness of the wave boundary layer at C-nodes    [m]
!*zaroughatc_sed*    REAL Apparent bottom roughness at C-nodes               [m]
!*zroughatc_sed*     REAL Bottom roughness at C-nodes                        [m]
!*zroughatu_sed*     REAL Bottom roughness at U-nodes                        [m]
!*zroughatv_sed*     REAL Bottom roughness at V-nodes                        [m]
!*wfall*             REAL Settling velocity (per sediment fraction)        [m/s]
!*ws_cst*            REAL Constant particle settling velocity (per sediment
!                         fraction                                         [m/s]
!
!******************************************************************************
!

SAVE

END MODULE sedarrays
