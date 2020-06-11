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

MODULE switches
!************************************************************************
!
! *switches* Model switches
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)switches.f90  V2.11.2
!
! $Date: 2018-07-23 10:21:02 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1167 $
!
! Description - 
!
!************************************************************************
!
IMPLICIT NONE

!---grid
INTEGER :: iopt_grid_htype = 1, iopt_grid_nodim = 3, iopt_grid_sph = 0, &
         & iopt_grid_vtype = 1, iopt_grid_vtype_transf = 0
!---interpolation
INTEGER :: iopt_arrint_depths = 1, iopt_arrint_hreg = 0, iopt_arrint_vreg = 0, &
         & iopt_arrint_3D = 0
!---hydrodynamics
INTEGER :: iopt_curr = 2, iopt_curr_wfall = 1, iopt_mode_2D, iopt_mode_3D
!---density
INTEGER :: iopt_dens = 0, iopt_dens_convect = 0, iopt_dens_grad = 1, &
         & iopt_sal = 0, iopt_temp = 0, iopt_temp_optic = 1, &
         & iopt_temp_sbc = 1
!---biology
INTEGER :: iopt_biolgy = 0
!---sediments
INTEGER ::  iopt_dar = 0, iopt_morph = 0, iopt_sed = 0, iopt_tidal_accel = 0
!---particle module
INTEGER :: iopt_part_model = 0, iopt_part_write = 0
!---bottom stress
INTEGER :: iopt_bstres_drag = 3, iopt_bstres_form = 2, iopt_bstres_nodim = 3, &
         & iopt_bstres_waves_bfric = 0
!---transport
INTEGER :: iopt_transp_full = 0
!---advection
INTEGER :: iopt_adv_scal = 3, iopt_adv_turb = 0, iopt_adv_tvd = 1, &
         & iopt_adv_2D = 1, iopt_adv_3D = 1
!---diffusion
INTEGER ::  iopt_hdif_coef = 0, iopt_hdif_lim = 0, iopt_hdif_scal = 0, &
          & iopt_hdif_turb = 0, iopt_hdif_2D = 0, iopt_hdif_3D = 0, &
          & iopt_kinvisc = 0, iopt_vdif_coef = 3, iopt_vdif_rot = 1
!---turbulence
INTEGER :: iopt_turb_alg = 1, iopt_turb_dis_bbc = 2, iopt_turb_dis_sbc = 2, &
         & iopt_turb_iwlim = 0, iopt_turb_kinvisc = 0, iopt_turb_lmix = 4, &
         & iopt_turb_ntrans = 1, iopt_turb_param = 2, iopt_turb_stab_form = 3, &
         & iopt_turb_stab_lev = 1, iopt_turb_stab_mod = 4, &
         & iopt_turb_stab_tke = 2, iopt_turb_tke_bbc = 2, &
         & iopt_turb_tke_sbc = 2
!---drying/wetting
INTEGER :: iopt_fld = 0, iopt_fld_alpha = 1
!---structures
INTEGER :: iopt_dischr = 0, iopt_dischr_land = 1, iopt_drycel = 0, &
         & iopt_thndam = 0, iopt_weibar = 0
!---explicit/implicit integration
INTEGER :: iopt_cor_impl = 1, iopt_hydro_impl = 0, iopt_scal_depos = 1, &
         & iopt_vadv_impl = 1, iopt_vdif_impl = 2
!---open boundary conditions
INTEGER :: iopt_obc_advflux = 1, iopt_obc_advrlx = 0, iopt_obc_bio = 0, &
         & iopt_obc_invbar = 0, iopt_obc_relax = 0, iopt_obc_sal = 0, &
         & iopt_obc_sed = 0, iopt_obc_temp = 0, iopt_obc_th = 0, &
         & iopt_obc_2D = 0, iopt_obc_2D_tang = 0, iopt_obc_3D = 0, &
         & iopt_obc_3D_tang = 0
!---astronomical tide
INTEGER :: iopt_astro_anal = 0, iopt_astro_pars = 0, iopt_astro_tide = 0
!---1-D applications
INTEGER:: iopt_sur_1D = 0
!---meteo surface forcing
INTEGER :: iopt_meteo = 0, iopt_meteo_data = 1, iopt_meteo_heat = 1, &
         & iopt_meteo_pres = 1, iopt_meteo_precip = 0, iopt_meteo_stres = 1
!---surface fluxes
INTEGER :: iopt_sflux_pars = 1, iopt_sflux_precip = 0, iopt_sflux_qlong = 1, &
         & iopt_sflux_qshort = 1, iopt_sflux_strat = 0
!---surface waves
INTEGER :: iopt_waves = 0, iopt_waves_couple = 0, iopt_waves_curr = 0, &
         & iopt_waves_dissip = 0, iopt_waves_extrapol = 1, iopt_waves_form = 1,&
         & iopt_waves_model = 0, iopt_waves_pres = 1
!---nesting
INTEGER :: iopt_nests = 0
!---parallel processing (MPI)
INTEGER :: iopt_MPI, iopt_MPI_abort = 0, iopt_MPI_comm_all = 2, &
         & iopt_MPI_comm_coll = 0, iopt_MPI_comm_exch = 2, &
         & iopt_MPI_comm_full = 0, iopt_MPI_comm_gath = 2, &
         & iopt_MPI_comm_scat = 2, iopt_MPI_partit = 1, iopt_MPI_sync = 0 
!---output
INTEGER :: iopt_out_anal = 0, iopt_out_avrgd = 0, iopt_out_tsers = 1
!---model coupling (MCT)
INTEGER :: iopt_MCT = 0
!---netCDF
INTEGER :: iopt_CDF, iopt_CDF_abort = 0, iopt_CDF_fill = 0, &
         & iopt_CDF_format = 1, iopt_CDF_shared = 0, iopt_CDF_sync = 0, &
         & iopt_CDF_tlim = 1
!---multigrid
INTEGER :: iopt_mg_smoother = 1, iopt_mg_cycle = 1, iopt_mg_prolong = 2
!---random generator
INTEGER :: iopt_rng_seed = 0
!---verification procedure
INTEGER :: iopt_verif

SAVE

!
! Name              Type    Purpose
!------------------------------------------------------------------------------
!*iopt_adv_scal*    INTEGER Selects advection scheme for scalar equations
!                    = 0 => No advection
!                    = 1 => Upwind
!                    = 2 => Lax-Wendroff or central
!                    = 3 => TVD
!*iopt_adv_tdv*     INTEGER Selects TVD limiting function
!                    = 1 => Superbee limiter
!                    = 2 => Monotone limiter
!*iopt_adv_turb*    INTEGER Selects advection scheme for turbulence equations
!                           (0/1/2/3)
!*iopt_adv_2D*      INTEGER Selects advection scheme for 2-D momentum
!                           equations (0/1/2/3)
!*iopt_adv_3D*      INTEGER Selects advection scheme for 3-D momentum
!                           equations (0/1/2/3)
!*iopt_arrint_depths*INTEGER Selects type of interpolation for water depths
!                    = 1 => Linear
!                    = 2 => Minimum depth
!*iopt_arrint_hreg* INTEGER Selects type of horizontal array interpolation
!                    = 0 => Uniform
!                    = 1 => With non-uniform weight factors
!*iopt_arrint_vreg* INTEGER Selects type of vertical array interpolation
!                    = 0 => Uniform
!                    = 1 => With non-uniform weight factors
!*iopt_arrint_3D*   INTEGER Selects dimension of mask or weight factor in some
!                           array interpolations
!                    = 0 => 2-D masks or weights
!                    = 1 => 3-D masks or weights
!*iopt_astro_anal*  INTEGER Disables/enables, if iopt_astro_pars > 0, the use
!                           of astronomical phases for harmonic analysis (0/1)
!*iopt_astro_pars*  INTEGER Type of evaluation of astronomical parameters
!                    = 0 => Astronomical argument set to zero, nodal factors
!                           set to 1
!                    = 1 => Evaluate astronomical phases at a given time and
!                           reference longitude, nodal factors are set to 1

!                    = 2 => Evaluate astronomical phases and nodal factors at a
!                           given time and reference longitude
!*iopt_astro_tide*  INTEGER Disable/enables inclusion of astronomical tidal
!                           force (0/1)
!*iopt_biolgy*      INTEGER Disables/enables activation of a biological module
!                           (0/1)
!*iopt_bstres_drag* INTEGER Type of formulation for bottom drag coefficient
!                    = 0 => not used
!                    = 1 => constant
!                    = 2 => from file
!                    = 3 => using constant roughness length
!                    = 4 => using roughness length from file
!*iopt_bstres_form* INTEGER Type of bottom stress formulation
!                    = 0 => zero bottom stress
!                    = 1 => linear law
!                    = 2 => quadratic law
!*iopt_bstres_nodim*INTEGER Type of current used in bottom stress formulation
!                    = 2 => depth-mean currents
!                    = 3 => 3-D bottom currents
!*iopt_bstres_waves_bfric*
!                   INTEGER Selects type of wave-current interaction model
!                           for the bottom stress
!                    = 0 => Disabled
!                    = 1 => Soulsby and Clarke (2005)
!                    = 2 => Malarkey and Davis (2012)
!                    = 3 => Grant and Madsen (1986)
!                    = 4 => Christoffersen and Johnsson (1985)
!*iopt_CDF*         INTEGER Disables/enables netCDF input/output (0/1)
!*iopt_CDF_abort*   INTEGER Disables/enables abortion of the program as soon
!                           as an error is detected in a netCDF routines (0/1)
!*iopt_CDF_fill*    INTEGER Enables/disables the use of fill values (0/1)
!*iopt_CDF_format*  INTEGER Selects file format of a netCDF file
!                    = 1 => classic format file
!                    = 2 => 64-bit offset format
!*iopt_CDF_shared*  INTEGER Disables/enables sharing of data in the same netCDF
!                           file between writing and reading processes at the
!                           expense of some larger CPU time (0/1).
!*iopt_CDF_sync*    INTEGER Disables/enables synchronisation of the disk copy of
!                           a netCDF user output file with in-memory buffers
!                           by calling NF90_SYNC before each read and after
!                           each write (0/1). Advantage is that data loss is
!                           minimised in case of abnormal termination.
!                           Disadvantage is a significantly larger CPU time.
!*iopt_CDF_tlim*    INTEGER Sets the type of time dimension in user output 
!                           netCDF files
!                    = 1 => fixed value. Data loss is prevented in case of
!                           abnormal termination.
!                    = 2 => unlimited time dimension. Advantage is that
!                           additional time records can be appended after
!                           completion of the run. Disadvantage is data loss is
!                           not prevented in case of abnormal termination
!                           (unless iopt_CDF_sync is set).
!*iopt_cor_impl*    INTEGER Type of time integration for Coriolis terms
!                    = 0 => fully explicit (theta_cor=0)
!                    = 1 => mixed implicit/explicit (theta_cor between 0 and 1)
!                    = 2 => fully implicit (theta_cor=1)
!*iopt_curr*        INTEGER Type of hydrodynamics
!                    = 0 => current and elevations are zero
!                    = 1 => currents are constant in time but may be
!                           non-uniform in space, elevations are zero
!                    = 2 => currents and elevations are updated in time
!*iopt_curr_wfall*  INTEGER Type of implementation of settling velocities for
!                           variables representing particulate matter
!                    = 1 => settling enabled without horizontal correction
!                    = 2 => settling enables with correction for horizontal
!                           advective currents
!*iopt_dar*         INTEGER Disables/enables activation of the dredging
!                           and relocation module (0/1)
!*iopt_dens*        INTEGER Selects equation of state
!                    = 0 => uniform density, zero expansion coefficients
!                    = 1 => linear equation of state
!                    = 2 => McDougall et al. (2006) without pressure effects
!                    = 3 => McDougall et al. (2006) including pressure effects
!*iopt_dens_convect*INTEGER Disables/enables convective adjustment (0/1)
!*iopt_dens_grad*   INTEGER Selects formulation for baroclinic gradient
!                    = 0 => No baroclinic gradient
!                    = 1 => Sigma second order
!                    = 2 => z-level method
!                    = 3 => Cubic H method
!*iopt_dischr*      INTEGER Disables/enables discharge module (0/1)
!*iopt_dischr_land* INTEGER Selects the procedure when discharges are located
!                           on a land cell
!                    = 1 => Land points are flagged (discharge disabled)
!                    = 2 => Discharge location is moved to the nearest wet
!                           point on the model grid 
!*iopt_drycel*      INTEGER Disables/enables dry cells (0/1)
!*iopt_fld*         INTEGER Selects drying/wetting algorithm (0/2)
!                    = 0 => drying/wetting disabled
!                    = 1 => without dynamic mask
!                    = 2 => with dynamic mask
!*iopt_fld_alpha*   INTEGER Type of formulation for the alpha factor in the
!                           drying/wetting algorithm
!                    = 0 => no reduction applied for advection (alpha=1)
!                    = 1 => power law (alpha=1 above critical depth; =0 at
!                           minimum depth)
!                    = 2 => no advection below the critical water depth
!*iopt_grid_htype*  INTEGER Type of horizontal model grid
!                    = 1 => uniform rectangular grid
!                    = 2 => non-uniform rectangular grid
!                    = 3 => arbitrary curvilinear grid
!*iopt_grid_node*   INTEGER Nodel type of global input grid
!                    = 1 => mean water depths at C-nodes,
!                           sigma-coordinates at W-nodes
!                    = 2 => mean water depths at UV-nodes,
!                           sigma-coordinates at UVW-nodes
!*iopt_grid_nodim*  INTEGER Dimension of model grid
!                    = 1 => 1-D vertical water column model grid
!                    = 2 => 2-D vertically intergrated grid
!                    = 3 => 3-D grid 
!*iopt_grid_sph*    INTEGER Type of coordinates
!                    = 0 => Cartesian
!                    = 1 => spherical
!*iopt_grid_vtype*   INTEGER Type of vertical grid
!                    = 1 => uniform sigma-grid
!                    = 2 => non-uniform sigma-grid
!                    = 3 => (non-uniform) s-grid
!*iopt_grid_vtype_transf*INTEGER Type of transformation for vertical grid
!                    = 0  => uniform (iopt_grid_vtype=1) or user-defined
!                    = 11 => log-transformation at the bottom following Davies
!                            and Jones (1991) if iopt_grid_vtype=2
!                    = 12 => log-transformation at the surface following Davies
!                            and Jones (1991) if iopt_grid_vtype=2
!                    = 13 => Burchard and Bolding (2002) if iopt_grid_vtype=2
!                    = 21 => Song and Haidvogel (1994) if iopt_grid_vtype=3
!*iopt_hdif_coef*   INTEGER Selects type of scheme for horizontal diffusion
!                           coefficients
!                    = 0 => horizontal diffusion disabled
!                    = 1 => constant coefficients
!                    = 2 => Smagorinsky coefficients
!*iopt_hdif_lim*   INTEGER  Selects type of limiting factor for scalar
!                           horizontal diffusion
!                    = 0 => disabled
!                    = 1 => Gerdes et aL. (1991)
!                    = 2 => Danabasoglu and McWilliams (1995)
!*iopt_hdif_scal*   INTEGER Selects type of horizontal diffusion scheme for
!                           scalars
!                    = 0 => disabled
!                    = 1 => Laplacian diffusion
!                    = 2 => diffusion along geopotential surfaces
!                    = 3 => isoneutral (along isopycnal) diffusion 
!*iopt_hdif_turb*   INTEGER Disables/enables horizontal diffusion scheme for
!                           turbulence equations (0/1)
!*iopt_hdif_2D*     INTEGER Disables/enables horizontal diffusion scheme for
!                           2-D momentum equations (0/1)
!*iopt_hdif_3D*     INTEGER Disables/enables horizontal diffusion scheme for
!                           3-D momentum equations (0/1)
!*iopt_hydro_impl*  INTEGER Type of time integration for the hydrodynamics
!                    = 0 => explicit scheme using mode splitting
!                    = 1 => implicit scheme using multigrid algorithm without 
!                           mode splitting
!*iopt_kinvisc*     INTEGER Selects type of kinematic viscosity
!                    = 0 => uniform value
!                    = 1 => from ITTC (1978) equation
!*iopt_MCT*         INTEGER Disables/enables model coupling using MCT 
!                           (0/1). This switch is only defined internally and
!                           cannot be rest by the user.
!*iopt_mg_cycle*    INTEGER Type of cycle for multigrid (1/2)
!                    = 1 => V-cycle
!                    = 2 => W-cycle
!*iopt_mg_prolong*  INTEGER Type of multigrid prolongation operator (1/2)
!                    = 1 => injection
!                    = 2 => bilinear interpolation
!*iopt_mg_smoother* INTEGER Type of smoother for multigrid (1/2)
!                    = 1 => Jacobi
!                    = 2 => Gauss-Seidel with Red-Black ordering
!*iopt_meteo*       INTEGER Disables/enables meteo input (0/1)
!*iopt_meteo_data*  INTEGER Type of meteo input
!                    = 1 => uwindatc, vwindatc, atmpres, airtemp, relhum,
!                           cloud_cover
!                    = 2 => usstresatc, vsstresatc, atmpres, qnonsol, qrad
!*iopt_meteo_heat*  INTEGER Disables/enables input of temperature or heat flux
!                           data (0/1)
!*iopt_meteo_precip*INTEGER Type of precipitation data input
!                    = 0 => no input
!                    = 1 => precipitation
!                    = 2 => evaporation minus precipitation
!*iopt_meteo_pres*  INTEGER Disables/enables input of atmospheric pressure (0/1)
!*iopt_meteo_stres* INTEGER Disables/enables input of wind or surface stress
!                           data (0/1)
!*iopt_mode_2D*     INTEGER Type of 2-D mode calculations
!                    = 0 => 2-D mode disabled
!                    = 1 => 2-D currents and elevations are initialised, but
!                           not updated in time
!                    = 2 => 2-D currents and elevations are initialised and
!                           updated in time
!*iopt_mode_3D*     INTEGER Type of 3-D mode (currents) calculations
!                    = 0 => 3-D mode disabled
!                    = 1 => 3-D currents are initialised, but not updated in
!                           time
!                    = 2 => 3-D currents are initialised and updated in time
!*iopt_morph*       INTEGER Disables/enables morphological calculation
!                           (0/1)
!*iopt_MPI*         INTEGER Disables/enables parallel processing with MPI
!                           (0/1). This switch is only defined internally and
!                           cannot be rest by the user.
!*iopt_MPI_abort*   INTEGER Disables/enables abortion of program when an error
!                           is detected in a MPI call (0/1)
!*iopt_MPI_comm_all*INTEGER Communication type for all to all operations
!                    = 1 => blocking, standard send
!                    = 2 => blocking, synchronous send
!                    = 3 => non-blocking, standard send
!                    = 4 => non-blocking, synchronous send
!*iopt_MPI_comm_coll*INTEGER Disables/enables use of collective calls where
!                            implemented (0/1)
!*iopt_MPI_comm_exch*INTEGER Communication type for exchange operations (1/5)
!                    = 1 => blocking, standard send
!                    = 2 => blocking, synchronous send
!                    = 3 => non-blocking, standard send
!                    = 4 => non-blocking, synchronous send
!                    = 5 => send-receive blocking calls
!*iopt_MPI_comm_full*INTEGER Disables/enables exchange of "full" (ensemble)
!                            arrays (0/1)
!*iopt_MPI_comm_gath*INTEGER Communication type for gather (all to one)
!                            operations (1/4)
!                    = 1 => blocking, standard send
!                    = 2 => blocking, synchronous send
!                    = 3 => non-blocking, standard send
!                    = 4 => non-blocking, synchronous send
!*iopt_MPI_comm_scat*INTEGER Communication type for scatter (one to all)
!                            operations (1/4)
!                    = 1 => blocking, standard send
!                    = 2 => blocking, synchronous send
!                    = 3 => non-blocking, standard send
!                    = 4 => non-blocking, synchronous send
!*iopt_MPI_partit*  INTEGER Type of domain decomposition
!                    = 1 => Regular partition according to the values of
!                           nprocsx, nprocsy
!                    = 2 => Externally defined
!*iopt_MPI_sync*    INTEGER Disables/enables synchronisation of blocking MPI
!                           calls (0/1)
!*iopt_nests*       INTEGER Disables/enables nested output (0/1)
!*iopt_obc_advflx*  INTEGER Type of open boundary condition for cross-stream
!                           advective fluxes in the momentum equations
!                    = 1 => Zero gradient condition
!                    = 2 => Upwind scheme
!*iopt_obc_advrlx*  INTEGER Disables/enables relaxation of advection for
!                           momentum in open boundary zones (0/1) 
!*iopt_obc_bio*     INTEGER General type of conditions at open boundaries for
!                           biological state variables
!                   = 0 => default conditions at all points
!                   = 1 => non-default conditions for all least one point
!*iopt_obc_invbar*  INTEGER Disables/enables inverse barometric effect at open
!                           boundaries (0/1)
!*iopt_obc_relax*   INTEGER Disables/enables relaxation at open boundaries
!                           (0/1)
!*iopt_obc_sal*     INTEGER General type of conditions at open boundaries for
!                           salinity
!                   = 0 => default conditions at all points
!                   = 1 => non-default conditions for all least one point
!*iopt_obc_sed*     INTEGER General type of conditions at open boundaries for
!                           sediments
!                   = 0 => default conditions at all points
!                   = 1 => non-default conditions for all least one point
!*iopt_obc_temp*    INTEGER General type of conditions at open boundaries for
!                           temperature
!                   = 0 => default conditions at all points
!                   = 1 => non-default conditions for all least one point
!*iopt_obc_th*      INTEGER Disables/enables the use of the Thatcher-Harleman
!                           open boundary condition (0/1)
!*iopt_obc_2D*      INTEGER General type of conditions at open boundaries for
!                           2-D mode
!                   = 0 => default conditions at all points
!                   = 1 => non-default conditions for all least one point
!*iopt_obc_2D_tang* INTEGER Disables/enables tangential open boundary
!                           conditions for the 2-D mode (0/1)
!*iopt_obc_3D*      INTEGER General type of conditions at open boundaries for
!                           the 3-D mode
!                   = 0 => default conditions at all points
!                   = 1 => non-default conditions for all least one point
!                          input for 2-D mode in initial conditions file
!*iopt_obc_3D_tang* INTEGER Disables/enables tangential open boundary
!                           conditions for the 3-D mode (0/1)
!*iopt_out_avrgd*   INTEGER Disables/enables time-averaged output (0/1)
!*iopt_out_anal*    INTEGER Disables/enables harmonic analysis and output (0/1)
!*iopt_out_tsers*   INTEGER Disables/enables time series output (0/1)
!*iopt_part_model*  INTEGER Type of mode for particle tracking
!                    = 0 => particle module is disabled
!                    = 1 => particle module runs online;
!                           COHERENS runs in serial
!                    = 2 => particle module runs online;
!                           COHERENS runs in parallel;
!                    = 3 => particle module runs offine in forward mode and
!                           reads data from previous COHERENS run
!                    = 4 => particle module runs offine in backward mode and
!                           reads data from previous COHERENS run
!*iopt_part_write*  INTEGEr Disables/enables writing of COHERENS data for the
!                           particle model (0/1)
!*iopt_rng_seed*    INTEGER Defines type of initial seed for all generators
!                    = 1 => The same default seed is used for all generators
!                           (the same sequence of random numbers is produced
!                           for all simulations
!                    = 2 => The initial seeds are generated depending on daytime
!                           (different runs produce different series of random
!                           numbers)   
!*iopt_sal*         INTEGER Type of salinity field
!                    = 0 => uniform in space and time
!                    = 1 => salinity is constant in time but may be non-uniform
!                           in space
!                    = 2 => salinity is updated in time
!*iopt_sed*         INTEGER Selects type of sediment module
!                    = 0 => Sediment model disabled
!                    = 1 => Sediment transport module without flocculation
!                    = 2 => Sediment transport module with flocculation
!*iopt_scal_depos   INTEGER Determines deposition flux in the transport module
!                    = 0 => No scalar deposition
!                    = 1 => First order discretization scalar deposition
!                    = 2 => Second order discretization scalar deposition
!*iopt_sflux_pars*   INTEGER Formulation for (neutral) surface drag and heat
!                           exchange coefficientS
!                    = 1 => generic
!                    = 2 => Large and Pond (1981)
!                    = 3 => Smith and Banke (1975)
!                    = 4 => Geernaert et al. (1986)
!                    = 5 => Kondo (1975)
!                    = 6 => Wu (1980)
!                    = 7 => Charnock's relation
!                    = 8 => North Sea formulation
!                    = 9 => Large ane Yeager (2004)
!*iopt_sflux_precip* INTEGER Selects how precipation/evaporation data are
!                            implemented if iopt_meteo_precip>0
!                    = 0 => disabled
!                    = 1 => as surface flux for salinity
!                    = 2 => as source/sink term in continuity equation
!                    = 3 => as surface flux for salinity and source/sink in
!                           continuity equation 
!*iopt_sflux_qlong*  INTEGER Selects formulation for surface longwave radiation
!                    = 1 => Brunt (1932), Gill (1982)
!                    = 2 => Clark et al. (1974)
!                    = 3 => Bignami et al. (1995)
!*iopt_sflux_qshort* INTEGER Selects formulation for surface shortwave
!                            radiative flux
!                    = 1 => Reed (1977), Rosati and Miyakoda (1988)
!                    = 2 => Zillman (1972)
!                    = 3 => Shine (1984)
!*iopt_sflux_strat* INTEGER Disables/enables stratification dependence for
!                           surface transfer coefficients
!*iopt_sur_1D*      INTEGER Selects surface forcing for 1-D applications
!                    = 0 => No surface forcing
!                    = 1 => Surface forcing for slopes and/or elevations
!*iopt_tidal_accel* INTEGER Disables/enables tidal acceleration for
!                           morphological calculations (0/1).
!*iopt_temp*        INTEGER Type of temperature field
!                    = 0 => uniform in space and time
!                    = 1 => temperature is constant in time but may be!
!                           non-uniform in space
!                    = 2 => temperature is updated in time
!*iopt_temp_optic*  INTEGER Disables/enables optical module (0/1)
!*iopt_temp_sbc*    INTEGER Type of surface boundary for temperature
!                    = 1 => Neumann
!                    = 2 => Dirichlet at one-half grid distance from the
!                           surface
!                    = 3 => Dirichlet at the surface
!*iopt_thndam*      INTEGER Disables/enables thin dams (0/1)
!*iopt_transp_full* INTEGER Type of update of a series of advection-diffusion
!                           equations
!                    = 0 => Each equation is updated separately
!                    = 1 => All equations are simultaneously updated
!*iopt_turb_alg*    INTEGER Type of algebraic formulation if iopt_vdif=1
!                    = 1 => Pacanowski-Philander
!                    = 2 => Munk-Anderson
!                    = 3 => Eddy coefficents proportional to the magnitude
!                           of the depth-averaged current times depth and
!                           a damping function for stratification
!                    = 4 => Eddy coefficients proportional to the squared
!                           magnitude of the depth-averaged current times
!                           a damping function for stratification
!                    = 5 => Eddy coefficients proportional to the magnitude
!                           of the depth-averaged current times the depth
!                           of the bottom boundary layer and a damping
!                           function for stratification
!                    = 6 => Eddy coefficients using parabolic eddy viscosity
!                           profile
!*iopt_turb_dis_bbc* INTEGER Type of bottom b.c. for dissipation rate
!                    = 1 => Neumann
!                    = 2 => Dirichlet
!*iopt_turb_dis_sbc*INTEGER Type of surface b.c. for dissipation rate
!                    = 1 => Neumann
!                    = 2 => Dirichlet
!*iopt_turb_iwlim*  INTEGER Selects type of internal wave mixing scheme
!                    = 0 => disabled
!                    = 1 => using limiting conditions
!                    = 2 => Large et al. scheme
!*iopt_turb_kinvisc*INTEGER Selects type of background mixing for momentum
!                    = 0 => User-defined background value
!                    = 1 => Kinematic viscosity coefficient
!*iopt_turb_lmix*   INTEGER Selects algebraic mixing length formulation
!                    = 1 => parabolic law
!                    = 2 => "modified" parabolic law
!                    = 3 => Xing's formulation
!                    = 4 => using Blackadar asymptotic form
!*iopt_turb_ntrans* INTEGER Selects number of transport equations to be solved
!                    = 0 => equilibrium method ("Mellor-Yamada level 2")
!                    = 1 => turbulence energy equation + mixing length
!                           selected by iopt_turb_lmix
!                    = 2 => k-eps or k-kl model
!*iopt_turb_param*  INTEGER Selects type of second turbulence variable for
!                           transport
!                    = 1 => mixing length 
!                    = 2 => dissipation rate
!*iopt_turb_stab_form*INTEGER Selects form of stability functions
!                    = 1 => constant
!                    = 2 => Munk-Andserson type
!                    = 3 => from Reynolds closure
!*iopt_turb_stab_lev* INTEGER Selects level for stability functions if
!                             iopt_turb_stab_form=3
!                    = 1 => quasi-equilibrium
!                    = 2 => non-equilibrium
!*iopt_turb_stab_mod* INTEGER Selects type of Reynolds closure
!                    = 1 => Mellor and Yamada
!                    = 2 => Kantha and Clayson
!                    = 3 => Burchard and Baumert
!                    = 4 => Hossain and Rodi
!                    = 5 => Canuto A
!                    = 6 => Canuto B
!*iopt_turb_stab_tke* INTEGER Selects type of stability function for t.k.e.
!                             diffusion
!                    = 1 => constant value
!                    = 2 => proportional to momentum stability function
!                    = 3 => Daly and Harlow formulation
!*iopt_turb_tke_bbc* INTEGER Type of bottom b.c. for turbulent energy
!                    = 1 => Neumann
!                    = 2 => Dirichlet
!*iopt_turb_tke_sbc*INTEGER Type of surface b.c. for turbulent energy
!                    = 1 => Neumann
!                    = 2 => Dirichlet
!*iopt_vadv_impl*   INTEGER Type of time integration for vertical advection
!                           (0/1/2)
!*iopt_vdif_coef*   INTEGER Selects type of vertical diffusion scheme
!                    = 0 => Vertical diffusion disabled
!                    = 1 => Uniform diffusion coefficients
!                    = 2 => Algebraic formulations
!                    = 3 => Second order turbulence closure
!*iopt_vdif_rot*    INTEGER Disables/enables inclusion of vertical diffusion
!                           due to rotated diffusion (0/1).
!                           Only used when iopt_hdif_scal >= 2.
!*iopt_vdif_impl*   INTEGER Type of time integration for vertical diffusion
!                           (0/1/2)
!*iopt_verif*       INTEGER Selects settings of test cases
!                    = 0 => original settings
!                    = 1 => modified settings for the verification procedure
!*iopt_waves*        INTEGER Disables/enables wave effects (0/1)
!*iopt_waves_couple*INTEGER Type of wave coupling
!                    = 0 => no coupling
!                    = 1 => one-way coupling (Coherens is receiving from the
!                           wave model)
!                    = 2 => one-way coupling (Coherens is sending to the wave
!                           model)
!                    = 3 => two-way coupling
!*iopt_waves_curr*  INTEGER Disables/enables wave/current interaction
!                           (Stokes drift) in the water column (0/1)
!*iopt_waves_dissip*INTEGER Disables/enables wave dissipation in the momentum
!                           equations (0/1)
!*iopt_waves_extrapol*INTEGER Disables/enables extrapolation if the model/data
!                           grid points are outside the data/model domain (0/1).
!*iopt_waves_form*  INTEGER Selects evaluation of wave parameters
!                     = 1 => using significant wave height and peak wave number
!                     = 2 => as random waves, obtained from the wave model by
!                            integration over the wave spectrum
!*iopt_waves_model* INTEGER Selects type of wave model
!                     = 0 => no coupling with external wave model
!                     = 1 => SWAN
!*iopt_waves_pres*  INTEGER Disables/enables wave pressure in the momentum
!                           equations (0/1)
!*iopt_weibar*      INTEGER Disables/enables weirs and barriers (0/1)
!
!************************************************************************
!

END MODULE switches
