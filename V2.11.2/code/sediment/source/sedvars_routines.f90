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

MODULE sedvars_routines
!************************************************************************
!
! *sedvars_routines* Utility routines for sediment model variables
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)sedvars_routines.f90  V2.11.1
!
! $Date: 2018-05-23 10:09:17 +0200 (Wed, 23 May 2018) $
!
! $Revision: 1138 $
!
! Description -
!
! Reference -
!
! Routines - inquire_sedvar, set_sedfiles_atts, set_sedvars_atts
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS

!============================================================================

SUBROUTINE inquire_sedvar(varid,f90_name,standard_name,comment,long_name,units,&
                        & node,vector_name,data_type,fill_value,nrank,&
                        & global_dims,local_dims,halo_dims,varatts)
!************************************************************************
!
! *inquire_sedvar* Obtain information about a specific sediment model variable
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)sedvars_routines.f90  V2.11.1
!
! Description -
!
! Calling program - inquire_var
!
! Module calls -
!
!************************************************************************
!
USE darpars
USE darswitches
USE datatypes
USE gridpars
USE iopars
USE morphpars
USE morphswitches
USE sedids
USE sedpars
USE switches
USE syspars

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lenname), INTENT(OUT) :: f90_name
CHARACTER (LEN=lendesc), INTENT(OUT) :: comment, long_name, standard_name, &
                                      & vector_name
CHARACTER (LEN=lenunit), INTENT(OUT) :: units
CHARACTER (LEN=*), INTENT(OUT) :: node
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(OUT) :: data_type, nrank
TYPE (VariableAtts), INTENT(OUT) :: varatts
INTEGER, INTENT(OUT), DIMENSION(4) :: halo_dims
INTEGER, INTENT(OUT), DIMENSION(5) :: global_dims, local_dims
REAL, INTENT(OUT), OPTIONAL :: fill_value

!
! Name           Type    Purpose
!------------------------------------------------------------------------------
!*varid*         INTEGER Variable id
!*f90_name*      CHAR    FORTRAN name
!*standard_name* CHAR    netCDF CF compliant standard name
!*comment*       CHAR    Miscellaneous information about the data or methods
!*long_name*     CHAR    Long descriptive name
!*units*         CHAR    Variable unit
!*node*          CHAR    Nodal type of model array
!*vector_name*   CHAR    Associated vector name
!*data_type*     INTEGER FORTRAN type of model variable
!*fill_value*    REAL    Fill value for missing data
!*nrank*         INTEGER Rank of model array
!*global_dims*   INTEGER Global shape of model array (excluding halos)
!*halo_dims*     INTEGER Halo size of model array (WESN directions)
!*local_dims*    INTEGER Local shape of model array (excluding halos)
!*varatts*       DERIVED Attributes of model array
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=lenname) :: fname
CHARACTER (LEN=lendesc) :: commt, lname, sname, vname
CHARACTER (LEN=lenunit) :: unit
CHARACTER (LEN=lennode) :: cnode
INTEGER :: dtype, nodim, nolocs
INTEGER,  DIMENSION(4) :: nhdims
INTEGER,  DIMENSION(5) :: ngdims, nldims
REAL :: fvalue

!
!1. Initialise
!-------------
!

unit = '_'; cnode = 'C'; vname = ''
dtype = float_type; nodim = 0
ngdims = 0; nhdims = 0; nldims = 0
fvalue = real_fill
!commt = 'standard name constructed following the netCDF CF-1.6 guidelines'
commt = ''
IF (iopt_dar.GT.0) THEN
   nolocs = MERGE(nocampaigns,nodredgingsites,iopt_dar_dredging_criterion.EQ.3)
ENDIF

!
!2. Sediments
!------------
!

SELECT CASE (varid)

!
!2.1 Sediment particle attributes
!--------------------------------
!

CASE (iarr_bstres_cr_cst)
   sname = 'particle_critical_stress_at_sea_floor'
   fname = 'bstres_cr_cst'
   lname = 'Bottom stress threshold for erosion'
   unit = 'm2 s-2'
   cnode  = ''
   nodim = 1
   ngdims(1) = nf
CASE (iarr_dp)
   sname = 'particle_diameter'
   fname = 'dp'
   lname = 'Sediment particle diameter'
   unit = 'm'
   cnode  = ''
   nodim = 1
   ngdims(1) = nf
CASE (iarr_rhos)
   sname = 'particle_density'
   fname = 'rhos'
   lname = 'Sediment particle density'
   unit = 'kg m-3'
   cnode  = ''
   nodim = 1
   ngdims(1) = nf
CASE (iarr_ws_cst)
   sname = 'particle_settling_velocity'
   fname = 'ws_cst'
   lname = 'Settling velocity'
   unit = 'm s-1'
   cnode  = ''
   nodim = 1
   ngdims(1) = nf

!
!2.2 Sediment concentrations
!---------------------------
!

CASE (iarr_ceq)
   sname = 'equilibrium_sediment_concentration_per_fraction'
   fname = 'ceq'
   lname = 'Equilibrium sediment concentration for 2-D transport'
   unit = 'kg m-3'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_cref)
   sname = 'reference_sediment_concentration_at_sea_floor_per_fraction'
   fname = 'cref'
   lname = 'Bottom reference concentration'
   unit = 'kg m-3'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
   nhdims = (/1,0,1,0/)
CASE (iarr_ctot)
   sname = 'total_sediment_mass_concentration'
   fname = 'ctot'
   lname = 'total sediment mass concentration'
   unit = 'kg m-3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,1,1,1/)
CASE (iarr_cvol)
   sname = 'sediment_mass_concentration_per_fraction'
   fname = 'cvol'
   lname = 'Sediment mass concentration per fraction'
   unit = 'kg m-3'
   nodim = 4
   nldims(1:4) = (/ncloc,nrloc,nz,nf/)
   ngdims(1:4) = (/nc,nr,nz,nf/)
   nhdims = nhalo
CASE (iarr_densw)
   sname = 'sea_water_density_without_sediment'
   fname = 'densw'
   lname = 'Sea water density without sediment'
   unit = 'kg m-3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
   
!
!2.3 Flocculation module
!-----------------------
!

CASE (iarr_cnump)
   sname = 'sediment_number_concentration_per_fraction'
   fname = 'cnump'
   lname = 'Sediment number concentration per fraction'
   unit = 'm-3'
   nodim = 4
   nldims(1:4) = (/ncloc,nrloc,nz,nf/)
   ngdims(1:4) = (/nc,nr,nz,nf/)
   nhdims = nhalo

CASE (iarr_floc_dens)
   sname = 'macrofloc_density'
   fname = 'floc_dens'
   lname = 'Macrofloc density'
   unit = 'kg m-3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1

CASE (iarr_floc_dia)
   sname = 'macrofloc_diameter'
   fname = 'floc_dia'
   lname = 'Macrofloc diameter'
   unit = 'micron'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)

CASE (iarr_floc_nc)
   sname = 'number_of_microfloc_in_macrofloc'
   fname = 'floc_nc'
   lname = 'Number of microflocs in a macrofloc'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1

CASE (iarr_floc_F)
   sname = 'macrofloc_mass_concentration'
   fname = 'floc_F'
   lname = 'Macrofloc mass concentration'
   unit = 'kg m-3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)

CASE (iarr_floc_P)
   sname = 'microfloc_mass_concentration'
   fname = 'floc_P'
   lname = 'Microfloc mass concentration'
   unit = 'kg m-3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)

CASE (iarr_floc_T)
   sname = 'number_concentration_microflocs_bound_in_macroflocs'
   fname = 'floc_T'
   lname = 'Number concentration of microflocs bound in macroflocs'
   unit = 'm-3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   
!
!2.4 Sediment loads
!------------------
!
!---bed load
CASE (iarr_qbedatc)
   sname = 'volumetric_bed_load_transport_per_fraction'
   fname = 'qbedatc'
   lname = 'Volumetric bed load transport'
   unit = 'kg m-1 s-1'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_qbedatu)
   sname = 'volumetric_bed_load_transport_per_fraction_at_u_nodes'
   fname = 'qbedatu'
   lname = 'X-component of volumetric bed load transport'
   vname = 'Volumetric bed load transport'
   unit = 'kg m-1 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
   nhdims = (/0,1,0,0/)
CASE (iarr_qbedatv)
   sname = 'volumetric_bed_load_transport_per_fraction_at_v_nodes'
   fname = 'qbedatv'
   lname = 'Y-component of volumetric bed load transport'
   vname  = 'Volumetric bed load transport'
   unit = 'kg m-1 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
   nhdims = (/0,0,0,1/)

!---suspended load
CASE (iarr_qsusatc)
   sname = 'volumetric_suspended_load_transport'
   fname = 'qsusatc'
   lname = 'Volumetric suspended load transport'
   unit = 'kg m-1 s-1'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_qsusatu)
   sname = 'volumetric_suspended_load_transport_at_u_nodes'
   fname = 'qsusatu'
   lname = 'Volumetric suspended load transport at U-nodes'
   unit = 'kg m-1 s-1'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
   nhdims = (/1,1,0,0/)
CASE (iarr_qsusatv)
   sname = 'volumetric_suspended_load_transport_at_v_nodes'
   fname = 'qsusatv'
   lname = 'Volumetric suspended load transport at V-nodes'
   unit = 'kg m-1 s-1'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
   nhdims = (/0,0,1,1/)

!---total load
CASE (iarr_qtotatc)
   sname = 'volumetric_total_load_transport_per_fraction'
   fname = 'qtotatc'
   lname = 'Volumetric total sediment transport'
   unit = 'kg m-1 s-1'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_qtotatu)
   sname = 'volumetric_total_load_transport_per_fraction_at_u_nodes'
   fname = 'qtotatu'
   lname = 'X-component of volumetric total sediment transport'
   vname ='Volumetric total sediment transport'
   unit = 'kg m-1 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
   nhdims = (/0,1,0,0/)
CASE (iarr_qtotatv)
   sname = 'volumetric_total_load_transport_per_fraction_at_v_nodes'
   fname = 'qtotatv'
   lname = 'Y-component of volumetric total sediment transport'
   vname = 'Volumetric total sediment transport'
   unit = 'kg m-1 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
   nhdims = (/0,0,0,1/)

!
!2.5 Bottom fluxes
!-----------------
!
   
CASE (iarr_bed_fraction)
   sname = 'sediment_fraction_at_sea_floor'
   fname = 'bed_fraction'
   lname =' Volume fraction of sediment in the bed'
   unit = '1'
   nodim = 4
   nldims(1:4) = (/ncloc,nrloc,nf,nb/)
   ngdims(1:4) = (/nc,nr,nf,nb/)
CASE (iarr_bottom_sed_dep)
   sname = 'sediment_volume_deposited_per_unit_area_and_time'
   fname = 'bottom_sed_dep'
   lname = 'Sediment volume deposited per unit area and time'
   unit = 'm s-1'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_bottom_sed_ero)
   sname = 'sediment_volume_eroded_per_unit_area_and_time'
   fname = 'bottom_sed_ero'
   lname = 'Sediment volume eroded per unit area and time'
   unit = 'm s-1'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_bottom_sed_flux)
   sname = 'fractional_sediment_flux_at_sea_floor'
   fname = 'bottom_sed_flux'
   lname = 'Net volumetric sediment flux per fraction at the bed'
   unit = 'm s-1'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_height_c)
   sname = 'reference_height_per_fraction'
   fname = 'heigth_c'
   lname = 'Reference height for bottom concentrations'
   unit = 'm'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_t_equil)
   sname = 'adaptation_time_scale_for_bottom_sediment_flux_per_fraction'
   fname = 't_equil'
   lname = 'Dimensionless adaptation timescale for bottom flux in 2-D transport'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)

!
!2.6 Bottom stress related arrays
!--------------------------------
!

CASE (iarr_bdragcoefatc_sed)
   sname = 'drag_coefficient_at_sea_floor_at_c_nodes'
   fname = 'bdragcoefatc_sed'
   lname = 'Bottom drag coefficient'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_bfricvel_sed)
   sname = 'friction_velocity_at_sea_floor'
   fname = 'bfricvel_sed'
   lname = 'Bottom friction velocity'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bfricvel_mean_sed)
   sname = 'mean_current-wave_wave_friction_velocity_at_sea_floor'
   fname = 'bfricvel_mean_sed'
   lname = 'Mean current-wave bottom friction velocity'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bfricvel_wav_sed)
   sname = 'wave_friction_velocity_at_sea_floor'
   fname = 'bfricvel_wav_sed'
   lname = 'Bottom wave friction'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bstres_cr)
   sname = 'critical_shear_stress_per_fraction'
   fname = 'bstres_cr'
   lname = 'Critical shear stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
   nhdims = (/1,0,1,0/)
CASE (iarr_bstresatc_sed)
   sname = 'stress_at_sea_floor'
   fname = 'bstresatc_sed'
   lname = 'Bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_bstresatc_mean_sed)
   sname = 'mean_current_wave_stress_at_sea_floor'
   fname = 'bstresatc_mean_sed'
   lname = 'Mean current-wave bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bstresatc_wav_sed)
   sname = 'wave_stress_at_sea_floor'
   fname = 'bstresatc_wave_sed'
   lname = 'Wave bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bstresatu_sed)
   sname = 'stress_at_sea_floor_at_U_nodes'
   fname = 'bstresatc_sed'
   lname = 'Bottom stress at U-nodes'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bstresatv_sed)
   sname = 'stress_at_sea_floor_at_V_nodes'
   fname = 'bstresatc_sed'
   lname = 'Bottom stress at V-nodes'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_d50_bed)
   sname = 'median_particle_diameter_at_sea_floor'
   fname = 'd50_bed'
   lname = 'Median particle diameter at the bed'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_fwave_sed)
   sname = 'wave_friction_factor'
   fname = 'fwave_sed'
   lname = 'Wave friction factor'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_ubstresatc_sed)
   sname = 'upward_x_stress_at_sea_floor_at_C_nodes'
   fname = 'ubstresatc_sed'
   lname = 'X-component of bottom stress at C-nodes'
   vname = 'Bottom stress'
   unit = 'Pa'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,0/)
CASE (iarr_ubstresatu_sed)
   sname = 'upward_x_stress_at_sea_floor_at_U_nodes'
   fname = 'ubstresatu_sed'
   lname = 'X-component of bottom stress at U-nodes'
   vname = 'Bottom stress'
   unit = 'Pa'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,0/)
CASE (iarr_vbstresatc_sed)
   sname = 'upward_y_stress_at_sea_floor_at_C_nodes'
   fname = 'vbstresatc_sed'
   lname = 'Y-component of bottom stress_at_C_nodes'
   vname = 'Bottom stres'
   unit = 'Pa'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims =(/0,0,0,1/)
CASE (iarr_vbstresatv_sed)
   sname = 'upward_y_stress_at_sea_floor_at_V_nodes'
   fname = 'vbstresatv_sed'
   lname = 'Y-component of bottom stress at V-nodes'
   vname = 'Bottom stres'
   unit = 'Pa'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims =(/0,0,0,1/)
CASE (iarr_wavethickatc_sed)
   sname = 'wave_boundary_layer_thickness'
   fname = 'wavethickatc_sed'
   lname = 'Wave boundary layer thickness '
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_zroughatc_sed)
   sname = 'sea_floor_roughness_length' 
   fname = 'zroughatc_sed'
   lname = 'Bottom roughness length'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_zroughatu_sed)
   sname = 'sea_floor_roughness_length_at_U_nodes' 
   fname = 'zroughatu_sed'
   lname = 'Bottom roughness length at U-nodes'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_zroughatv_sed)
   sname = 'sea_floor_roughness_length_at_V_nodes' 
   fname = 'zroughatv_sed'
   lname = 'Bottom roughness length at V-nodes'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
   
!
!2.7 Diffusion coefficients
!--------------------------
!

CASE (iarr_beta_sed)
   sname = 'ratio_momentum_to_sediment_diffusivity'
   fname = 'beta_sed'
   lname = 'Ratio between sediment and momentum eddy diffusivities'
   unit = '1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_vdiffcoef_sed)
   sname = 'sea_water_turbulent_sediment_diffusivity_per_fraction'
   fname = 'vdiffcoef_sed'
   lname = 'Vertical eddy diffusivity for sediments'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:4) = (/ncloc,nrloc,nz+1,nf/)
   ngdims(1:4) = (/nc,nr,nz+1,nf/)

!
!2.8 Bed slope arrays
!--------------------
!
   
CASE (iarr_bed_slope_x)
   sname = 'x_bed_slope'
   fname = 'bed_slope_x'
   lname = 'X-component of the bed slope at C-nodes'
   vname = 'Bed slope at C-nodes'
   nodim  = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bed_slope_x_atu)
   sname = 'x_bed_slope_at_u_nodes'
   fname = 'bed_slope_x_atu'
   lname = 'X-component of the bed slope at U-nodes'
   vname = 'Bed slope at U-nodes'
   cnode = 'U'
   nodim  = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bed_slope_x_atv)
   sname = 'x_bed_slope_at_v_nodes'
   fname = 'bed_slope_x_atv'
   lname = 'X-component of the bed slope at V-nodes'
   vname = 'Bed slope at V-nodes'
   cnode = 'V'
   nodim  = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bed_slope_y)
   sname = 'y_bed_slope'
   fname = 'bed_slope_y'
   lname = 'Y-component of the bed slope at C-nodes'
   vname = 'Bed slope at C-nodes'
   nodim  = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bed_slope_y_atu)
   sname = 'y_bed_slope_at_u_nodes'
   fname = 'bed_slope_y_atu'
   lname = 'Y-component of the bed slope at U-nodes'
   vname = 'Bed slope at C-nodes'
   cnode = 'U'
   nodim  = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bed_slope_y_atv)
   sname = 'y_bed_slope_at_v_nodes'
   fname = 'bed_slope_y_atv'
   lname = 'Y-component of the bed slope at V-nodes'
   vname = 'Bed slope at C-nodes'
   cnode = 'V'
   nodim  = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)

!
!2.9 Miscellaneous
!-----------------
!

CASE (iarr_beta_state_sed)
   sname = 'sediment_coefficient_for_sea_water_expansion_per_fraction'
   fname = 'beta_state_sed'
   lname = 'Sediment expansion coefficient'
   nodim = 3
   nldims(1:4) = (/ncloc,nrloc,nz,nf/)
   ngdims(1:4) = (/nc,nr,nz,nf/)
   nhdims(1:4) = (/1,1,1,1/)
CASE (iarr_dar_sediment)
   sname = 'dredging_relocation_sediment_sources'
   fname = 'dar_sediment'
   lname = 'Sediment sources due to dredging or relocation'
   unit = 'm3 s-1'
   nodim = 4
   nldims(1:4) = (/ncloc,nrloc,nz,nf/)
   ngdims(1:4) = (/nc,nr,nz,nf/)
CASE (iarr_wfall)
   sname = 'setlling_velocity_per_fraction'
   fname = 'wfall'
   lname = 'Settling velocity'
   unit = 'm s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:4) = (/ncloc,nrloc,nz+1,nf/)
   ngdims(1:4) = (/nc,nr,nz+1,nf/)
   nhdims = (/1,0,1,0/)

!
!2.10 Open boundary conditions
!-----------------------------
!
   
CASE (iarr_obccnumpatu)
   sname= 'storage_array_for_sediment_number_concentration_at_u_boundaries'
   fname = 'obccnumpatu'
   lname = 'Storage array for sediment number concentrations at U-open boundaries'
   unit = 'm-3'
   cnode = 'U'
   nodim = 3
   ngdims(1:4) = (/nobu,nz,3,nf/)
CASE (iarr_obccnumpatv)
   sname= 'storage_array_for_sediment_number_concentration_at_v_boundaries'
   fname = 'obccnumpatv'
   lname = 'Storage array for sediment number concentrations at V-open boundaries'
   unit = 'm-3'
   cnode = 'V'
   nodim = 3
   ngdims(1:4) = (/nobu,nz,3,nf/)
CASE (iarr_obcsedatu)
   sname= 'storage_array_for_sediment_per_fraction_at_u_boundaries'
   fname = 'obcsedatu'
   lname = 'Storage array for sediment concentrations at U-open boundaries'
   unit = 'm3 m-3'
   cnode = 'U'
   nodim = 3
   ngdims(1:4) = (/nobu,nz,3,nf/)
CASE (iarr_obcsedatv)
   sname= 'storage_array_for_sediment_per_fraction_at_v_boundaries'
   fname = 'obcsedatv'
   lname = 'Storage array for sediment concentrations at V-open boundaries'
   unit = 'm3 m-3'
   cnode = 'V'
   nodim = 3
   ngdims(1:4) = (/nobv,nz,3,nf/)

!
!3. Morphology
!-------------
!

CASE (iarr_active_layer_thickness)
   sname = 'active_layer_thickness'
   fname = 'active_layer_thickness'
   lname = 'Active layer thickness'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bed_layer_thickness)
   sname = 'bed_layer_thickness'
   fname = 'bed_layer_thickness'
   lname = 'Bed layer thickness'
   unit = 'm'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nb/)
   ngdims(1:3) = (/nc,nr,nb/)
CASE (iarr_bed_porosity)
   sname = 'bed_porosity'
   fname = 'bed_porosity'
   lname = 'Porosity of the bed layers'
   unit = '1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nb/)
   ngdims(1:3) = (/nc,nr,nb/)
CASE (iarr_bed_update)
   sname = 'bed_level_change'
   fname = 'bed_update'
   lname = 'Change in bed level elevation between two time steps'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bed_update_dep)
   sname = 'bed_level_change_per_fraction_due_to_deposition'
   fname = 'bed_update_dep'
   lname = 'Change in bed level elevation due to deposition per fraction'
   unit = 'm'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_bed_update_dep_old1)
   sname = 'bed_level_change_per_fraction_due_to_deposition_at_previous_'//&
         & 'time_level'
   fname = 'bed_update_dep_old1'
   lname = 'Change in bed level elevation due to deposition per fraction '//&
         & 'at previous time level'
   unit = 'm'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_bed_update_dep_old2)
   sname = 'bed_level_change_per_fraction_due_to_deposition_at_second_'//&
         & 'previous_time_level'
   fname = 'bed_update_dep_old2'
   lname = 'Change in bed level elevation due to deposition per fraction '//&
         & 'at second previous time level'
   unit = 'm'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_bed_update_ero)
   sname = 'bed_level_change_per_fraction_due_to_erosion'
   fname = 'bed_update_ero'
   lname = 'Change in bed level elevation due to erosion per fraction'
   unit = 'm'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_bed_update_ero_old1)
   sname = 'bed_level_change_per_fraction_due_to_erosion_at_previous_time_level'
   fname = 'bed_update_ero_old1'
   lname = 'Change in bed level elevation due to erosion per fraction '//&
         & 'at previous time level'
   unit = 'm'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_bed_update_ero_old2)
   sname = 'bed_level_change_per_fraction_due_to_erosion_at_second_'//&
         & 'previous_time_level'
   fname = 'bed_update_ero_old2'
   lname = 'Change in bed level elevation due to erpdion per fraction at '//&
         & 'second previous time level'
   unit = 'm'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_bed_update_int)
   sname = 'accumulated_bed_level_change'
   fname = 'bed_update_int'
   lname = 'Accumulated change in bed level elevation'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_dar_morph)
   sname = 'dredging_relocation_source_or_sinks_in_bed_layer'
   fname = 'dar_morph'
   lname = 'Amount of sediments removed or added in the sea bed by dredging '//&
         & 'or relocation'
   unit = 'm3 s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_depth_guess)
   sname = 'stored_water_depth'
   fname = 'depth_guess'
   lname = 'Stored total water depth for tidal acceleration' 
   unit = 'm'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,number_tidal_steps/)
   ngdims(1:3) = (/nc,nr,number_tidal_steps/)
   nhdims = 1
CASE (iarr_sed_avail)
   sname = 'sediment_layer_depth_per_fraction'
   fname = 'sed_avail'
   lname = 'Depth of sediment layer per fraction'
   unit = 'm'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nf/)
   ngdims(1:3) = (/nc,nr,nf/)
CASE (iarr_sed_avail_tot)
   sname = 'total_sediment_layer_depth'
   fname = 'sed_avail_tot'
   lname = 'Total depth of sediment layer'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_sediment_balance)
   sname = 'domain_integrated_sediment_balance_error'
   fname = 'sediment_balance'
   lname = 'Domain integrated sediment balance error'
   unit = 'm3 m-3'
   nodim = 0
CASE (iarr_sediment_obc)
   sname = 'net_amount_of_sediment_exchange_open_boundaries'
   fname = 'sediment_obc'
   lname = 'Net amount of volumetric sediment exchange at open boundaries'
   unit = 'm3 m-3'
   nodim = 0
CASE (iarr_sediment_vol)
   sname = 'net_amount_of_sediment_exchange_at_sea_bed'
   fname = 'sediment_vol'
   lname = 'Net amount of volumetric sediment exchange at sea bed'
   unit = 'm3 m-3'
   nodim = 0
CASE (iarr_tidalstep)
   sname = 'storage_time_steps_in_characteristic_period'
   fname = 'tidalstep'
   lname = 'Storage time steps in characteristic_period'
   dtype = int_type
   unit = '1'
   nodim = 1
   nldims(1) = number_tidal_steps
   ngdims(1) = number_tidal_steps
CASE (iarr_umvel_guess)
   sname = 'stored_depth_mean_x_sea_water_velocity'
   fname = 'umvel_guess'
   lname = 'Stored X-component of depth-mean current for tidal acceleration' 
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,number_tidal_steps/)
   ngdims(1:3) = (/nc,nr,number_tidal_steps/)
   nhdims  = (/nhalo,nhalo,1,1/)
CASE (iarr_vmvel_guess)
   sname = 'stored_depth_mean_y_sea_water_velocity'
   fname = 'vmvel_guess'
   lname = 'Stored Y-component of depth-mean current for tidal acceleration' 
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,number_tidal_steps/)
   ngdims(1:3) = (/nc,nr,number_tidal_steps/)
   nhdims = (/1,1,nhalo,nhalo/)
CASE (iarr_zeta_guess)
   sname = 'stored_sea_surface_height'
   fname = 'zeta_guess'
   lname = 'Stored surface elevation for tidal acceleration' 
   unit = 'm'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,number_tidal_steps/)
   ngdims(1:3) = (/nc,nr,number_tidal_steps/)
   nhdims = 1

!
!4. Dredging and relocation
!-------------------------- 
!
!---campaigns
CASE (iarr_camp_campdepth)
   sname = 'campaign_dredging_thickness'
   fname = 'camp_campdep'
   lname = 'Thickness of the dredged layer during a campaign'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_camp_campfractions)
   sname = 'campaign_fractions'
   fname = 'camp_campfractions'
   lname = 'Volume fractions of the dredged sediment during a campaign'
   unit = '1'
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,MaxFractions/)
   ngdims(1:2) = (/nolocs,MaxFractions/)
CASE (iarr_camp_camplevel)
   sname = 'campaign_dredging_level'
   fname = 'camp_camplevel'
   lname = 'Bed level up to which a site is dredged during a campaign'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_camp_CampaignStartDateTime)
   sname = 'campaign_starting_time'
   fname = 'camp_CampaignStartDateTime'
   lname = 'Campaign start time'
   dtype = char_type
   cnode = ''
   nodim = 2
   nldims(1:2) = (/lentime,nolocs/)
   ngdims(1:2) = (/lentime,nolocs/)
CASE (iarr_camp_campvolume)
   sname = 'campaign_volume'
   fname = 'camp_campvolume'
   lname = 'Volume of dredged sediment during a current campaign'
   unit = 'm3'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_camp_CoupleTime)
   sname = 'campaign_coupling_times'
   fname = 'camp_CoupleTime'
   lname = 'User-defined coupling times linked to a campaign'
   dtype = char_type
   cnode = ''
   nodim = 3
   nldims(1:3) = (/lentime,nolocs,MaxTimes/)
   ngdims(1:3) = (/lentime,nolocs,MaxTimes/)
CASE (iarr_camp_dredge)
   sname = 'campaign_dredging_flag'
   fname = 'camp_dredge'
   lname = 'Flag to decide for dredging in the current campaign'
   unit = 'log'
   dtype = log_type
   cnode = ''
   nodim = 1
CASE (iarr_camp_nr_dredgingsite)
   sname = 'dredging_site_id'
   fname = 'camp_nr_dredgingsite'
   lname = 'Number of the dredging site excavated in the current campaign'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_camp_nr_relocationsite)
   sname = 'relocation_sites_id'
   fname = 'camp_nr_relocationsite'
   lname = 'Number of the relocation site used in the current campaign'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_camp_nr_ship)
   sname = 'Campaign_ship_id'
   fname = 'camp_nr_ship_campaign'
   lname = 'Number of the ship that operates at the dredging site'
   dtype = int_type   
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_camp_nr_zones_coupled)
   sname = 'number_of_relocation_zones_coupled_to_campaign'
   fname = 'camp_nr_zonescoupled'
   lname = 'Number of relocation zones that are coupled to a campaign'
   cnode = ''
   dtype = int_type
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_camp_relocate)
   sname = 'campaign_relocation_flag'
   fname = 'camp_relocate'
   lname = 'Flag to decide for relocation in the current campaign'
   cnode = ''
   unit = 'log'
   dtype = log_type
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_camp_relocationsites)
   sname = 'number_of_relocation_sites_coupled_to_campaign'
   fname = 'camp_relocationsites'
   lname = 'Number of relocation sites visited during a campaign'
   dtype = int_type
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,MaxSites/)
   ngdims(1:2) = (/nolocs,MaxSites/)
CASE (iarr_camp_t_lag)
   sname = 'campaign_time_lag_dredging_relocation'
   fname = 'camp_t_lag'
   lname = 'Number of time steps in a campaign between start of dredging '//&
         & 'and start of relocation'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_camp_t_spread)
   sname = 'campaign_dredging_time_spread'
   fname = 'camp_t_spread'
   lname = 'Number of time steps to complete the dredging'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs

!---dredging sites
CASE (iarr_dredge_CoupleTime)
   sname = 'dredging_site_coupling_times'
   fname = 'dredge_CoupleTime'
   lname = 'Coupling times linked to a dredging site'
   dtype = char_type
   cnode = ''
   nodim = 3
   nldims(1:3) = (/lentime,nolocs,MaxTimes/)
   ngdims(1:3) = (/lentime,nolocs,MaxTimes/)
CASE (iarr_dredge_critical_volume)
   sname = 'dredging_criterion_critical_volume'
   fname = 'dredge_critical_volume'
   lname = 'Critical volume triggering a dredging campaign in case of '//&
         & 'criterion dredging'
   unit = 'm3'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_dredge_dredging_thickness)
   sname = 'maximum_dredging_depth_suctionhead'
   fname = 'dredge_dredging_thickness'
   lname = 'Maximum layer thickness a suctionhead dredger can remove from '//&
         & 'the bed during one pass'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_dredge_monitoring_frequency)
   sname = 'dredging_site_monitoring_frequency'
   fname = 'dredge_monitoring_frequency'
   lname = 'Number of time steps for monitoring of a dredging site'
   dtype = int_type
   cnode = ''
   nodim = 1
CASE (iarr_dredge_nr_coordpol)
   sname = 'number_of_points_on_polygon_dredging_site'
   fname = 'dredge_nr_coordpol'
   lname = 'Number of polygon points determining a dredging site'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_dredge_nr_coupled)
   sname = 'number_of_coupling_times'
   fname = 'dredge_nr_coupled'
   lname = 'Number of user-defined coupling times'
   dtype = int_type
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,MaxTimes/)
   ngdims(1:2) = (/nolocs,MaxTimes/)
CASE (iarr_dredge_nr_ship)
   sname = 'dredging_site_ship_id'
   fname = 'dredge_nr_ship_dredge'
   lname = 'Ship number operating at the dredging site'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_dredge_nr_zones_coupled)
   sname = 'number_of_relocation_zones_coupled_to_dredging_site'
   fname = 'dredge_nr_zones_coupled'
   lname = 'Number of relocation zones coupled to a dredging site'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_dredge_over_depth)
   sname = 'dredging_site_over_depth'
   fname = 'dredge_over_depth'
   lname = 'Over depth at the dredging site'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_dredge_relocationsites)
   sname = 'relocation_sites_coupled_to_dredging'
   fname = 'dredge_relocationsites'
   lname = 'Relocation sites coupled to the dredging site'
   dtype = int_type
   cnode = ''
   nodim  = 2
   nldims(1:2) = (/nolocs,MaxSites/)
   ngdims(1:2) = (/nolocs,MaxSites/)
CASE (iarr_dredge_target_depth)
   sname = 'dredging_site_target_depth'
   fname = 'dredge_target_depth'
   lname = 'Target depth at the dredging site'
   unit = 'm'
   cnode = ''
   nodim  = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_dredge_t_lag)
   sname = 'time_between_start_dredging_and_start_relocation'
   fname = 'dredge_t_lag'
   lname = 'Number of time steps between start of dredging and start of '//&
         & 'relocation'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_dredge_t_spread)
   sname = 'dredging_site_time_spread'
   fname = 'dredge_t_spread'
   lname = 'Number of time steps for completing a dredging'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_dredge_Xpol)
   sname = 'X_coordinates_polygon_dredging_site'
   fname = 'dredge_Xpol'
   lname = 'X-coordinates of the polygon determining a dredging site'
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,MaxCoordPol/)
   ngdims(1:2) = (/nolocs,MaxCoordPol/)
CASE (iarr_dredge_Ypol)
   sname = 'Y_coordinates_polygon_dredging_site'
   fname = 'dredge_Ypol'
   lname = 'Y-coordinates of the polygon determining a dredging site '
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,MaxCoordPol/)
   ngdims(1:2) = (/nolocs,MaxCoordPol/)

!---relocation sites
CASE (iarr_reloc_max_volume)
   sname = 'maximum_volume_relocation_site'
   fname = 'reloc_maxvol_rel'
   lname = 'Maximum volume that can be relocated at a relocation site'
   unit = 'm3'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_reloc_min_depth)
   sname = 'minimum_depth_relocation_site'
   fname = 'reloc_min_depth'
   lname = 'Minimum depth to be kept at a relocation site'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_reloc_nr_coordpol)
   sname = 'number_of_points_on_polygon_relocation_site'
   fname = 'reloc_nr_coordpol'
   lname = 'Number of polygon points determining a relocation site'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_reloc_Xpol)
   sname = 'X_coordinates_polygon_relocation_site'
   fname = 'reloc_Xpol'
   lname = 'X-coordinates of the polygon determining a relocation site '
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,MaxCoordPol/)
   ngdims(1:2) = (/nolocs,MaxCoordPol/)
CASE (iarr_reloc_Ypol)
   sname = 'Y_coordinates_polygon_relocation_site'
   fname = 'reloc_Ypol'
   lname = 'Y-coordinates of the polygon determining a relocation site '
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = ''
   nodim = 2
   nldims(1:2) =  (/nolocs,MaxCoordPol/)
   ngdims(1:2) =  (/nolocs,MaxCoordPol/)

!---ships
CASE (iarr_ship_D_dredge_max)
   sname = 'upper_heigth_above_sea_bed_for_plume_overflow_at_dredging_site'
   fname = 'ship_D_dredge_max'
   lname = 'Upper distance to the sea bed for overflow losses at the '//&
         & 'dredging site'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_D_dredge_min)
   sname = 'lower_heigth_above_sea_bed_for_plume_overflow_at_dredging_site'
   fname = 'ship_D_dredge_min'
   lname = 'Lower distance to the sea bed for overflow losses at the '//&
         & 'dredging site'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_D_relocate_max)
   sname = 'upper_heigth_above_sea_bed_for_plume_overflow_at_relocation_site'
   fname = 'ship_D_relocate_max'
   lname = 'Upper distance to the sea bed for overflow losses at the '//&
         & 'relocation site'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_D_relocate_min)
   sname = 'lower_heigth_above_sea_bed_for_plume_overflow_at_relocation_site'
   fname = 'ship_D_relocate_min'
   lname = 'Lower distance to the sea bed for overflow losses at the '//&
         & 'relocation site'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_H_dredge_max)
   sname = 'highest_depth_from_sea_surface_for_suctionhead_overflow_at_'//&
         & 'dredging_site'
   fname = 'ship_H_dredge_max'
   lname = 'Highest depth from the sea surface for suctionhead losses at '//&
         & 'the dredging site'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_H_dredge_min)
   sname = 'lowest_depth_from_sea_surface_for_suctionhead_overflow_at_'//&
         & 'dredging_site'
   fname = 'ship_H_dredge_min'
   lname = 'Lowest depth from the sea surface for suctionhead losses at '//&
         & 'the dredging site'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_H_relocate_max)
   sname = 'highest_depth_from_sea_surface_for_suctionhead_overflow_at_'//&
         & 'relocation_site'
   fname = 'ship_H_relocate_max'
   lname = 'Highest depth from the sea surface for suctionhead losses at '//&
         & 'the relocation site'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_H_relocate_min)
   sname = 'lowest_depth_from_sea_surface_for_suctionhead_overflow_at_'//&
         & 'relocation_site'
   fname = 'ship_H_relocate_min'
   lname = 'Lowest depth from the sea surface for suctionhead losses at '//&
         & 'the relocation site'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_p_overflow)
   sname = 'fractional_loss_dredged_sediment_by_ship_overflow'
   fname = 'ship_p_overflow'
   lname = 'Fractional loss of dredged sediment by overflow from the ship'
   unit = '1'
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,MaxFractions/)
   ngdims(1:2) = (/nolocs,MaxFractions/)
CASE (iarr_ship_p_passive)
   sname = 'fractional_loss_relocated_sediment_as_passive_plume'
   fname = 'ship_p_passive'
   lname = 'Fractional loss of relocated sediment as a passive plume'
   unit = '1'
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,MaxFractions/)
   ngdims(1:2) = (/nolocs,MaxFractions/)
CASE (iarr_ship_p_relocation)
   sname = 'fractional_loss_relocated_sediment_in_water_column'
   fname = 'ship_p_relocation'
   lname = 'Fractional loss of relocated sediment within the water column'
   unit = '1'
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,MaxFractions/)
   ngdims(1:2) = (/nolocs,MaxFractions/)
CASE (iarr_ship_p_suctionhead)
   sname = 'factional_loss_dredged_sediment_at_suction_head'
   fname = 'ship_p_suctionhead'
   lname = 'Fractional loss of dredged sediment at the suction head'
   unit = '1'
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,MaxFractions/)
   ngdims(1:2) = (/nolocs,MaxFractions/)
CASE (iarr_ship_R_circle)
   sname = 'radius_relocation_circle'
   fname = 'ship_R_circle'
   lname = 'Radius of the the circle in which the sediment is relocated'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_t_fill)
   sname = 'ship_filling_time'
   fname = 'ship_t_fill'
   lname = 'Number of time steps for filling a dredging ship'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_t_idle)
   sname = 'ship_idle_time'
   fname = 'ship_t_idle'
   lname = 'Number of time steps between ship arrival and start of dredging'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_t_relocate)
   sname = 'ship_relocation_time'
   fname = 'ship_t_relocate'
   lname = 'Number of time steps needed to empty the ship at the relocation '//&
         & 'site'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_t_sail_dredge)
   sname = 'ship_sailing_time_dredging_site'
   fname = 'ship_t_sail_dredge'
   lname = 'Number of time steps needed to sail to the dredging site'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_t_sail_relocate)
   sname = 'ship_sailing_time_to_relocation_site'
   fname = 'ship_t_sail_relocate'
   lname = 'Number of time steps needed to sail to the relocation site'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_t_wait)
   sname = 'ship_waiting_time'
   fname = 'ship_t_wait'
   lname = 'Number of time steps between the call for dredging and arrival '//&
         & 'of the dredging ship'
   dtype = int_type
   cnode = ''
   unit = 's'
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_U_ship)
   sname = 'ship_velocity'
   fname = 'ship_U_ship'
   lname = 'Sailing velocity of the ship'
   unit = 'm s-1'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_ship_V_ship)
   sname = 'ship_volume'
   fname = 'ship_V_ship'
   lname = 'Volume of the ship'
   unit = 'm3'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs

!----sub-zone dredging
CASE (iarr_subdred_nr_coordpol)
   sname = 'number_of_points_on_polygon_dredging_subzone'
   fname = 'subdred_nr_coordpol'
   lname = 'Number of polygon points determining a dredging sub-zone'
   dtype = int_type
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,MaxSubzones/)
   ngdims(1:2) = (/nolocs,MaxSubzones/)
CASE (iarr_subdred_nr_subzones)
   sname = 'number_of_subzones_inside_dredging_site'
   fname = 'subdred_nr_subzones'
   lname = 'Number of sub-zones inside a dredging site'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_subdred_Xpol)
   sname = 'X_coordinates_polygon_dredging_subzone'
   fname = 'subdred_Xpol'
   lname = 'X-coordinates of the polygon determining a dredging sub-zone'
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = ''
   nodim = 3
   nldims(1:3) = (/nolocs,MaxCoordPol,MaxSubzones/)
   ngdims(1:3) = (/nolocs,MaxCoordPol,MaxSubzones/)
CASE (iarr_subdred_Ypol)
   sname = 'Y_coordinates_polygon_dredging_subzone'
   fname = 'subdred_Ypol'
   lname = 'Y-coordinates of the polygon determining a dredging sub-zone'
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = ''
   nodim = 3
   nldims(1:3) = (/nolocs,MaxCoordPol,MaxSubzones/)
   ngdims(1:3) = (/nolocs,MaxCoordPol,MaxSubzones/)

!----sub-zone relocation
CASE (iarr_subrel_nr_coordpol)
   sname = 'number_of_points_on_polygon_relocation_subzone'
   fname = 'subrel_nr_coordpol'
   lname = 'Number of polygon points determining a relocation sub-zone'
   dtype = int_type
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,MaxSubzones/)
   ngdims(1:2) = (/nolocs,MaxSubzones/)
CASE (iarr_subrel_nr_subzones)
   sname = 'number_of_subzones_inside_relocation_site'
   fname = 'subrel_nr_subzones'
   lname = 'Number of sub-zones inside a relocation site'
   dtype = int_type
   nodim = 1
CASE (iarr_subrel_Xpol)
   sname = 'X_coordinates_polygon_relocation_subzone'
   fname = 'subrel_Xpol'
   lname = 'X-coordinates of the polygon determining a relocation sub-zone'
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = ''
   nodim = 3
   nldims(1:3) = (/nolocs,MaxCoordPol,MaxSubzones/)
   ngdims(1:3) = (/nolocs,MaxCoordPol,MaxSubzones/)
CASE (iarr_subrel_Ypol)
   sname = 'Y_coordinates_polygon_relocation_subzone'
   fname = 'subrel_Ypol'
   lname = 'Y-coordinates of the polygon determining a relocation sub-zone'
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = ''
   nodim = 3
   nldims(1:3) = (/nolocs,MaxCoordPol,MaxSubzones/)
   ngdims(1:3) = (/nolocs,MaxCoordPol,MaxSubzones/)

!---other arrays
CASE (iarr_dredging_depth)
   sname = 'depth_of_dredged_layer'
   fname = 'dredging_depth'
   lname = 'Depth of dredged sediment layer'
   unit = 'm'
   cnode = ''
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nolocs/)
   ngdims(1:3) = (/nc,nr,nolocs/)
CASE (iarr_D_dredge_max)
   sname = 'maximum_overloss_depth'
   fname = 'D_dredge_max'
   lname = 'Maximum depth for overloss during dredging'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_D_dredge_min)
   sname = 'minimum_overloss_depth'
   fname = 'D_dredge_min'
   lname = 'Minimum depth for overloss during dredging'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_D_relocate_max)
   sname = 'maximum_overflow_depth'
   fname = 'D_relocate_max'
   lname = 'Maximum depth for sediment overflow during relocation'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_D_relocate_min)
   sname = 'minimum_overflow_depth'
   fname = 'D_relocate_min'
   lname = 'Minimum depth for sediment overflow during relocation'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_fraction_ship)
   sname = 'ship_fraction_relocated_sediment'
   fname = 'fraction_ship'
   lname = 'Size fraction distribution of relocated sediment'
   unit =  '1'
   cnode = ''
   nodim = 2
   nldims(1:2) = (/noships,nf/)
   ngdims(1:2) = (/noships,nf/)
CASE (iarr_H_dredge_max)
   sname = 'maximum_height_above_sea_bed_for_suctionhead_losses'
   fname = 'H_dredge_max'
   lname = 'Maximum height above sea bed for suctionhead losses during dredging'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_H_dredge_min)
   sname = 'minimum_height_above_sea_bed_for_suctionhead_losses'
   fname = 'H_dredge_min'
   lname = 'Minimum height above sea bed for suctionhead losses during dredging'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_H_relocate_max)
   sname = 'maximum_height_above_sea_bed_for_passive_losses'
   fname = 'H_relocate_max'
   lname = 'Maximum height above sea bed for passive losses during relocation'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_H_relocate_min)
   sname = 'minimum_height_above_sea_bed_for_passive_losses'
   fname = 'H_relocate_min'
   lname = 'Minimum height above sea bed for passive losses during relocation'
   unit = 'm'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_itrack)
   sname = 'x_index_of_dredged_cells_during tracking'
   fname = 'itrack'
   lname = 'X-index of dredged cells during tracking'
   dtype = int_type
   cnode = ''
   nodim = 2
   nldims(1:2) = (/maxtrackpoints,nolocs/)
   ngdims(1:2) = (/maxtrackpoints,nolocs/)
CASE (iarr_jtrack)
   sname = 'y_index_of_dredged_cells_during tracking'
   fname = 'jtrack'
   lname = 'Y-index of dredged cells during tracking'
   dtype = int_type
   cnode = ''
   nodim = 2
   nldims(1:2) = (/maxtrackpoints,nolocs/)
   ngdims(1:2) = (/maxtrackpoints,nolocs/)
CASE (iarr_mask_dredging)
   sname = 'grid_mask_for_dredging_sites'
   fname = 'mask_dredging'
   lname = 'Grid mask for dredging sites'
   unit = 'log'
   dtype = log_type
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nodredgingsites/)
   ngdims(1:3) = (/nc,nr,nodredgingsites/)
CASE (iarr_mask_relocation)
   sname = 'grid_mask_for_relocation_sites'
   fname = 'mask_relocation'
   lname = 'Grid mask for relocation sites '
   unit = 'log'
   dtype = log_type
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,norelocationsites/)
   ngdims(1:3) = (/nc,nr,norelocationsites/)
CASE (iarr_mask_subzones_dredging)
   sname = 'grid_mask_for_dredging_sub_zones'
   fname = 'mask_subzones_dredging'
   lname = 'Grid mask for sub-zone dredging sites'
   unit = 'log'
   dtype = log_type
   nodim = 4
   nldims(1:4) = (/ncloc,nrloc,nodredgingsites,MaxSubzones/)
   ngdims(1:4) = (/nc,nr,nodredgingsites,MaxSubzones/)
CASE (iarr_mask_subzones_relocation)
   sname = 'grid_mask_for_relocation_sub_zones'
   fname = 'mask_subzones_relocation'
   lname = 'Grid mask for sub-zones relocation sites'
   unit = 'log'
   dtype = log_type
   nodim = 4
   nldims(1:4) = (/ncloc,nrloc,norelocationsites,MaxSubzones/)
   ngdims(1:4) = (/nc,nr,norelocationsites,MaxSubzones/)
CASE (iarr_nocouplingtimes)
   sname = 'number_of_dredging_relocation_coupling_times'
   fname = 'nocouplingtimes'
   lname = 'Number of coupling times for coupled dredging and relocation'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_notrackpoints)
   sname = 'number_of_ship_track_positions'
   fname = 'notrackpoints'
   lname = 'Number of ship track positions'
   dtype = int_type
   unit = ''
   cnode = ''
   nodim = 1
   nldims(1) = nodredgingsites
   ngdims(1) = nodredgingsites
CASE (iarr_p_overflow)
   sname = 'fractional_loss_of_overflow_sediment'
   fname = 'p_overflow'
   lname = 'Fractional loss of dredged sediment by overflow'
   unit = '1'
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,nf/)
   ngdims(1:2) = (/nolocs,nf/)
CASE (iarr_p_passive)
   sname = 'fractional_loss_passive_sediment_loss'
   fname = 'p_passive'
   lname = 'Fractional loss of relocated sediment giving passive plumes'
   unit = '1'
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,nf/)
   ngdims(1:2) = (/nolocs,nf/)
CASE (iarr_p_relocation)
   sname = 'fractional_loss_relocated_sediment'
   fname = 'p_relocation'
   lname = 'Fractional loss of relocated sediment by the ship'
   unit = '1'
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,nf/)
   ngdims(1:2) = (/nolocs,nf/)
CASE (iarr_p_suctionhead)
   sname = 'fractional_loss_suctionhead_sediment_loss'
   fname = 'p_suctionhead'
   lname = 'Suctionhead fractional losses of dredged sediment'
   unit = '1'
   cnode = ''
   nodim = 2
   nldims(1:2) = (/nolocs,nf/)
   ngdims(1:2) = (/nolocs,nf/)
CASE (iarr_track_points_counter)
   sname = 'number_of_track_positions'
   fname = 'track_points_counter'
   lname = 'Number of ship track positions'
   dtype = int_type
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_volume_ship)
   sname = 'ship_sediment_load'
   fname = 'volume_ship'
   lname = 'Amount of sediment loaded on the ship'
   unit = 'm3'
   cnode = ''
   nodim = 1
   nldims(1) = noships
   ngdims(1) = noships
CASE (iarr_V_dredge)
   sname = 'amount_of_dredged_sediment'
   fname = 'V_dredge'
   lname = 'Amount of dredged sediment'
   unit = 'm3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nolocs/)
   ngdims(1:3) = (/nc,nr,nolocs/)
CASE (iarr_V_dredge_max)
   sname = 'maximum_amount_of_tracked_dredged_sediment'
   fname = 'V_dredge_max'
   lname = 'Maximum amount of dredged sediment volume in case of tracked dredging'
   unit = 'm3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nolocs/)
   ngdims(1:3) = (/nc,nr,nolocs/)
CASE (iarr_V_dredge_rest)
   sname = 'amount_of_remaining_sediment_for_dredging'
   fname = 'V_dredge_rest'
   lname = 'Amount of sediment volume remaining for dredging'
   unit = 'm3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nolocs/)
   ngdims(1:3) = (/nc,nr,nolocs/)
CASE (iarr_V_dredge_total)
   sname = 'total_amount_of_dredged_sediment'
   fname = 'V_dredge_total'
   lname = 'Total amount of dredged sediment volume at the dredging site'
   unit = 'm3'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_V_lastcycle)
   sname = 'amount_of_dredged_sediment'
   fname = 'V_lastcycle'
   lname = 'Amount of sediment volume to be dredged during last cycle'
   unit = 'm3'
   cnode = ''
   nodim = 1
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_V_relocate_extra)
   sname = 'amount_of_non_relocated_sediment_at_relocation_site'
   fname = 'V_relocate_extra'
   lname = 'Amount of non-relocated sediment volume at the relocation site'
   unit = 'm3'
   cnode = ''
   nodim = 1
   nldims(1) = norelocationsites
   ngdims(1) = norelocationsites
CASE (iarr_V_relocate_total)
   sname = 'amount_of_relocated_sediment_at_relocation_site'
   fname = 'V_relocate_total'
   lname = 'Amount of relocated sediment volume at the relocation site'
   unit = 'm3'
   cnode = ''
   nodim = 1
   nldims(1) = norelocationsites
   ngdims(1) = norelocationsites
CASE (iarr_V_ship_rest)
   sname = 'amount_of_sediment_reamining_for_loading'
   fname = 'V_ship_rest'
   lname = 'Amount of sediment remaining volume for loading on the ship'
   unit = 'm3'
   cnode = ''
   nodim = 1 
   nldims(1) = nolocs
   ngdims(1) = nolocs
CASE (iarr_xtrack)
   sname = 'x_coordinate_of_ship_tracking_positions'
   fname = 'xtrack'
   lname = 'X-coordinates of tracking positions'
   unit = MERGE('m          ','degree_east',iopt_grid_sph.EQ.0)
   cnode = ''
   nodim = 2
   nldims(1:2) = (/maxtrackpoints,nolocs/)
   ngdims(1:2) = (/maxtrackpoints,nolocs/)
CASE (iarr_ytrack)
   sname = 'y_coordinate_of_ship_tracking_positions'
   fname = 'ytrack'
   lname = 'Y-coordinates of tracking positions'
   unit = MERGE('m           ','degree_north',iopt_grid_sph.EQ.0)
   cnode = ''
   nodim = 2
   nldims(1:2) = (/maxtrackpoints,nolocs/)
   ngdims(1:2) = (/maxtrackpoints,nolocs/)

END SELECT

!
!5. Store
!--------
!

f90_name = fname; standard_name = sname; comment = commt; long_name = lname
vector_name = vname; units = unit; node = cnode; data_type = dtype
fill_value = fvalue; nrank = nodim; global_dims = ngdims
local_dims = nldims; halo_dims = nhdims
varatts%standard_name = sname; varatts%comment = commt
varatts%f90_name = fname; varatts%long_name = lname
varatts%vector_name = vname; varatts%units = unit
varatts%data_type = dtype; varatts%nrank = nodim; varatts%global_dims = ngdims
varatts%local_dims = nldims; varatts%halo_dims = nhdims
varatts%fill_value = float_fill


RETURN

END SUBROUTINE inquire_sedvar

!============================================================================

SUBROUTINE set_sedfiles_atts(iddesc,ifil,iotype)
!************************************************************************
!
! *set_sedfiles_atts* Obtain global attributes of a sediment model file
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)sedvars_routines.f90  V2.11.1
!
! Description -
!
! Calling program - set_modfiles_atts
!
! Module calls -
!
!************************************************************************
!
USE darswitches
USE gridpars
USE iopars
USE morphswitches
USE sedpars
USE sedswitches
USE switches
USE syspars

IMPLICIT NONE
!
!* Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, iotype

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*   INTEGER File id
!*ifil*     INTEGER File number
!*iotype*   INTEGER I/O type of file
!             = 1 => input
!             = 2 => output
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=lendesc) :: title
INTEGER :: nocoords, novars


pglev = pglev + 1
procname(pglev) = 'set_sedfiles_atts'
IF (loglev1.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, TRIM(procname(pglev))
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, TRIM(procname(pglev))
   ENDIF
ENDIF

!
!1. Define attributes
!--------------------
!

SELECT CASE (iddesc)
!  ---initial conditions
   CASE (io_inicon,io_fincon)
      
      SELECT CASE (ifil)
!        --sediment water column model
         CASE (ics_sed)
            title = 'Sediment initial conditions'
            novars = 0
            IF (iopt_sed_mode.EQ.2.OR.iopt_sed_mode.EQ.3) THEN
               novars = novars + nf
            ENDIF
            IF (nfp.GT.1) novars = novars + nf
            IF (iopt_sed_rough.EQ.2) novars = novars + 1
            IF (iopt_obc_sed.EQ.1) THEN
               IF (nobu.GT.0) novars = novars + nf
               IF (nobv.GT.0) novars = novars + nf
            ENDIF
!        --morphology
         CASE (ics_morph)
            title = 'Morphological initial conditions'
            novars = 0
            IF (iopt_morph_fixed_layer.EQ.1) THEN
               novars = novars + 1
            ENDIF
            IF (iopt_morph_time_int.GT.1) THEN
               novars = novars + 2*nf
               IF (iopt_morph_time_int.EQ.3) novars = novars + 2*nf
            ENDIF
            novars = novars + 1
      END SELECT
         
      nocoords = 1

!  ---particle attributes
   CASE (io_sedspc)
      title = 'Particle attributes'
      novars = 4
      nocoords = 0

!  ---dredging/relocation attributes
   CASE (io_darspc)
      title = 'Dredging/relocation attributes'
      novars = MERGE (64,49,iopt_dar_dredging_criterion.EQ.3)
      nocoords = 0

END SELECT

!
!2. Store attributes
!-------------------
!

modfiles(iddesc,ifil,iotype)%title = title
modfiles(iddesc,ifil,iotype)%novars = novars
modfiles(iddesc,ifil,iotype)%nocoords = nocoords

IF (loglev2.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, logexit
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, logexit
   ENDIF
ENDIF
pglev = pglev - 1


RETURN

END SUBROUTINE set_sedfiles_atts

!========================================================================

SUBROUTINE set_sedvars_atts(iddesc,ifil,iotype,ivarid,nrank,nshape,data_type,&
                          & numvarsid,numvars)
!************************************************************************
!
! *set_sedvars_atts* Obtain variable attributes of a sediment model file
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)sedvars_routines.f90  V2.11.1
!
! Description -
!
! Calling progam - set_modvars_atts
!
! Module calls -
!
!************************************************************************
!
USE darpars
USE darswitches
USE datatypes
USE gridpars
USE iopars
USE modids
USE morphpars
USE morphswitches
USE sedids
USE sedpars
USE sedswitches
USE switches
USE syspars

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, iotype, numvars
INTEGER, INTENT(OUT), DIMENSION(numvars) :: data_type, ivarid, nrank, numvarsid
INTEGER, INTENT(OUT), DIMENSION(numvars,4) :: nshape 

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER File id
!*ifil*      INTEGER File number
!*iotype*    INTEGER I/O type of file
!             = 1 => input
!             = 2 => output
!*ivarid*    INTEGER Returned key variable ids
!*nrank*     INTEGER Returned variable ranks
!*nshape*    INTEGER Returned variable shapes
!*data_type* INTEGER Type of model variable
!*numvarsid* INTEGER Variable indeces in case the last dimension is a variable
!                    dimension
!*numvars*   INTEGER Number of variables (coordinate and data) in the forcing
!                    file
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: packing
INTEGER :: ivar, f, nolocs, n0
TYPE (FileParams) :: filepars


pglev = pglev + 1
procname(pglev) = 'set_sedvars_atts'
IF (loglev1.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, TRIM(procname(pglev))
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, TRIM(procname(pglev))
   ENDIF
ENDIF

!
!1. Initialise
!-------------
!

ivar = 2; ivarid = 0; nrank = -1; nshape = 1; numvarsid = -1
filepars = modfiles(iddesc,ifil,iotype)
data_type = MERGE(real_type,rlong_type,filepars%floattype.EQ.'S')
packing = filepars%packing

!
!2. Define attributes
!--------------------
!

SELECT CASE (iddesc)

!
!2.1 Initial conditions
!-----------------------
!

   CASE (io_inicon,io_fincon)

      SELECT CASE (ifil)
!        ---sediment water column model
         CASE (ics_sed)
            
!           --sediment concentrations
            IF (iopt_sed_mode.EQ.2.OR.iopt_sed_mode.EQ.3) THEN
               ivarid(ivar:ivar+nf-1) = MERGE(iarr_cvol,iarr_cnump,&
                                            & iopt_sed.EQ.1)
               nrank(ivar:ivar+nf-1) = MERGE(2,3,packing)
               nshape(ivar:ivar+nf-1,1) = MERGE(noseaatc,nc,packing)
               nshape(ivar:ivar+nf-1,2) = MERGE(nz,nr,packing)
               IF (.NOT.packing) nshape(ivar:ivar+nf-1,3) = nz
               numvarsid(ivar:ivar+nf-1) = (/(f,f=1,nf)/)
               ivar = ivar + nf
            ENDIF
!           --bed fractions
            IF (nfp.GT.1) THEN
               ivarid(ivar:ivar+nf-1) = iarr_bed_fraction
               nrank(ivar) = MERGE(2,3,packing)
               nshape(ivar:ivar+nf-1,1) = MERGE(noseaatc,nc,packing)
               nshape(ivar:ivar+nf-1,2) = MERGE(nb,nr,packing)
               IF (.NOT.packing) nshape(ivar:ivar+nf-1,3) = nb
               numvarsid(ivar:ivar+nf-1) = (/(f,f=1,nf)/)
               ivar = ivar + nf
            ENDIF
!           --skin roughness
            IF (iopt_sed_rough.EQ.2) THEN
               ivarid(ivar) = iarr_zroughatc_sed
               nrank(ivar) = MERGE(1,2,packing)
               nshape(ivar,1) = MERGE(noseaatc,nc,packing)
               IF (.NOT.packing) nshape(ivar,2) = nr
               ivar = ivar + 1
            ENDIF
               
!           --open boundary arrays
            IF (iopt_obc_sed.EQ.1) THEN
               IF (nobu.GT.0) THEN
                  ivarid(ivar:ivar+nf-1) = MERGE(iarr_obcsedatu,&
                                               & iarr_obccnumpatu,iopt_sed.EQ.1)
                  nrank(ivar:ivar+nf-1) = 3
                  nshape(ivar:ivar+nf-1,1) = nobu
                  nshape(ivar:ivar+nf-1,2) = nz
                  nshape(ivar:ivar+nf-1,3) = 3
                  ivar = ivar + nf
               ENDIF
               IF (nobv.GT.0) THEN 
                  ivarid(ivar:ivar+nf-1) = MERGE(iarr_obcsedatv,&
                                               & iarr_obccnumpatv,iopt_sed.EQ.1)
                  nrank(ivar:ivar+nf-1) = 3
                  nshape(ivar:ivar+nf-1,1) = nobv
                  nshape(ivar:ivar+nf-1,2) = nz
                  nshape(ivar:ivar+nf-1,3) = 3
                  ivar = ivar + nf
               ENDIF
            ENDIF

!        ---morphology
         CASE (ics_morph)
!           --bed layer thickness
            IF (iopt_morph_fixed_layer.EQ.1) THEN
               ivarid(ivar) = iarr_bed_layer_thickness
               nrank(ivar) = MERGE(2,3,packing)
               nshape(ivar,1) = MERGE(noseaatc,nc,packing)
               nshape(ivar,2) = MERGE(nb,nr,packing)
               IF (.NOT.packing) nshape(ivar,3) = nb
               ivar = ivar + 1
            ENDIF
!           --bed update
            IF (iopt_morph_time_int.GT.1) THEN
               ivarid(ivar:ivar+nf-1) = iarr_bed_update_dep_old1
               nrank(ivar:ivar+nf-1) = MERGE(1,2,packing)
               nshape(ivar:ivar+nf-1,1) = MERGE(noseaatc,nc,packing)
               IF (.NOT.packing) nshape(ivar:ivar+nf-1,2) = nr
               numvarsid(ivar:ivar+nf-1) = (/(f,f=1,nf)/)
               ivar = ivar + nf
               ivarid(ivar:ivar+nf-1) = iarr_bed_update_ero_old1
               nrank(ivar:ivar+nf-1) = MERGE(1,2,packing)
               nshape(ivar:ivar+nf-1,1) = MERGE(noseaatc,nc,packing)
               IF (.NOT.packing) nshape(ivar:ivar+nf-1,2) = nr
               numvarsid(ivar:ivar+nf-1) = (/(f,f=1,nf)/)
               ivar = ivar + nf
               IF (iopt_morph_time_int.EQ.3) THEN
                  ivarid(ivar:ivar+nf-1) = iarr_bed_update_dep_old2
                  nrank(ivar:ivar+nf-1) = MERGE(1,2,packing)
                  nshape(ivar:ivar+nf-1,1) = MERGE(noseaatc,nc,packing)
                  IF (.NOT.packing) nshape(ivar:ivar+nf-1,2) = nr
                  numvarsid(ivar:ivar+nf-1) = (/(f,f=1,nf)/)
                  ivar = ivar + nf
                  ivarid(ivar:ivar+nf-1) = iarr_bed_update_ero_old2
                  nrank(ivar:ivar+nf-1) = MERGE(1,2,packing)
                  nshape(ivar:ivar+nf-1,1) = MERGE(noseaatc,nc,packing)
                  IF (.NOT.packing) nshape(ivar:ivar+nf-1,2) = nr
                  numvarsid(ivar:ivar+nf-1) = (/(f,f=1,nf)/)
                  ivar = ivar + nf
               ENDIF
            ENDIF
            ivarid(ivar) = iarr_bed_update_int
            nrank(ivar) = MERGE(1,2,packing)
            nshape(ivar,1) = MERGE(noseaatc,nc,packing)
            IF (.NOT.packing) nshape(ivar,2) = nr
      END SELECT

!
!2.2 Particle characteristics
!----------------------------
!

   CASE (io_sedspc)
      ivarid = (/iarr_bstres_cr_cst,iarr_dp,iarr_rhos,iarr_ws_cst/)
      nrank = 1
      nshape(:,1) = nf

!
!2.3 Dredging/relocation attributes
!----------------------------------
!

   CASE (io_darspc)
      n0 = MERGE(15,0,iopt_dar_dredging_criterion.EQ.3)
      nolocs = MERGE(nocampaigns,nodredgingsites,&
                   & iopt_dar_dredging_criterion.EQ.3)

!     ---campaigns
      IF (n0.GT.0) THEN
         ivarid(1:n0) = (/&
           & iarr_camp_campdepth, iarr_camp_campfractions,&
           & iarr_camp_camplevel, iarr_camp_CampaignStartDateTime,&
           & iarr_camp_campvolume, iarr_camp_CoupleTime,&
           & iarr_camp_dredge, iarr_camp_nr_dredgingsite,&
           & iarr_camp_nr_relocationsite, iarr_camp_nr_ship,&
           & iarr_camp_nr_zones_coupled,iarr_camp_relocate,&
           & iarr_camp_relocationsites, iarr_camp_t_lag, iarr_camp_t_spread/)
         ivar_231: DO ivar=1,n0
            SELECT CASE(ivarid(ivar))
               CASE (iarr_camp_campfractions)
                  nrank(ivar) = 2
                  nshape(ivar,1:2) = (/nolocs,MaxFractions/)
               CASE (iarr_camp_CampaignStartDateTime)
                  nrank(ivar) = 2
                  nshape(ivar,1:2) = (/lentime,nolocs/)
               CASE (iarr_camp_CoupleTime)
                  nrank(ivar) = 3
                  nshape(ivar,1:3) = (/lentime,nolocs,MaxTimes/)
               CASE (iarr_camp_relocationsites)
                  nrank(ivar) = 2
                  nshape(ivar,1:2) = (/nolocs,MaxSites/)
               CASE DEFAULT
                  nrank(ivar) = 1
                  nshape(ivar,1) = nolocs
            END SELECT
            SELECT CASE(ivarid(ivar))
               CASE (iarr_camp_dredge,iarr_camp_relocate)
                  data_type(ivar) = log_type
               CASE (iarr_camp_nr_dredgingsite, iarr_camp_nr_relocationsite,&
                   & iarr_camp_nr_ship, iarr_camp_nr_zones_coupled, &
                   & iarr_camp_relocationsites, iarr_camp_t_lag,&
                   & iarr_camp_t_spread)
                  data_type(ivar) = int_type
               CASE (iarr_camp_CampaignStartDateTime,iarr_camp_CoupleTime)
                  data_type(ivar) = char_type
            END SELECT
         ENDDO ivar_231
      ENDIF

!     ---dredging sites
      ivarid(n0+1:n0+15) = (/&
          & iarr_dredge_CoupleTime,iarr_dredge_critical_volume,&    
          & iarr_dredge_dredging_thickness,iarr_dredge_monitoring_frequency,&
          & iarr_dredge_nr_coordpol,iarr_dredge_nr_coupled,&
          & iarr_dredge_nr_ship,iarr_dredge_nr_zones_coupled,&
          & iarr_dredge_over_depth,iarr_dredge_relocationsites,&
          & iarr_dredge_target_depth,iarr_dredge_t_lag,iarr_dredge_t_spread,&
          & iarr_dredge_Xpol,iarr_dredge_Ypol/)
      ivar_232: DO ivar=n0+1,n0+15
         SELECT CASE (ivarid(ivar))
            CASE (iarr_dredge_CoupleTime)
               nrank(ivar) = 3
               nshape(ivar,1:3) = (/lentime,nolocs,MaxTimes/)
            CASE (iarr_dredge_nr_coupled)
               nrank(ivar) = 2
               nshape(ivar,1:2) = (/nolocs,MaxTimes/)
            CASE (iarr_dredge_relocationsites)
               nrank(ivar)= 2
               nshape(ivar,1:2) = (/nolocs,MaxSites/)
            CASE (iarr_dredge_Xpol,iarr_dredge_Ypol)
               nrank(ivar)= 2
               nshape(ivar,1:2) = (/nolocs,MaxCoordPol/)
            CASE DEFAULT
               nrank(ivar) = 1
               nshape(ivar,1) = nolocs
         END SELECT
         SELECT CASE (ivarid(ivar))
            CASE (iarr_dredge_monitoring_frequency,&
                & iarr_dredge_nr_coordpol,iarr_dredge_nr_coupled,&
                & iarr_dredge_nr_ship,iarr_dredge_nr_zones_coupled,&
                & iarr_dredge_t_lag,iarr_dredge_t_spread)
               data_type(ivar) = int_type
            CASE (iarr_dredge_CoupleTime)
               data_type(ivar) = char_type
         END SELECT
      ENDDO ivar_232

!     ---relocation sites
      ivarid(n0+16:n0+20) = (/&
          & iarr_reloc_max_volume,iarr_reloc_min_depth, &
          & iarr_reloc_nr_coordpol,iarr_reloc_Xpol,iarr_reloc_Ypol/)
      ivar_233: DO ivar=n0+16,n0+20
         SELECT CASE (ivarid(ivar))
            CASE (iarr_reloc_Xpol,iarr_reloc_Ypol)
               nrank(ivar) = 2
               nshape(ivar,1:2) = (/nolocs,MaxCoordPol/)
            CASE DEFAULT
               nrank(ivar) = 1
               nshape(ivar,1) = nolocs
            END SELECT
            IF (ivarid(ivar).EQ.iarr_reloc_nr_coordpol) THEN
               data_type(ivar) = int_type
            ENDIF
      ENDDO ivar_233

!     ---ships
      ivarid(n0+21:n0+41) = (/&
          & iarr_ship_D_dredge_max,iarr_ship_D_dredge_min,&
          & iarr_ship_D_relocate_max,iarr_ship_D_relocate_min,&
          & iarr_ship_H_dredge_max,iarr_ship_H_dredge_min,&
          & iarr_ship_H_relocate_max,iarr_ship_H_relocate_min,&
          & iarr_ship_p_overflow,iarr_ship_p_passive,&
          & iarr_ship_p_relocation,iarr_ship_p_suctionhead,&
          & iarr_ship_R_circle,iarr_ship_t_fill,iarr_ship_t_idle,&
          & iarr_ship_t_relocate,iarr_ship_t_sail_dredge,&
          & iarr_ship_t_sail_relocate,iarr_ship_t_wait,iarr_ship_U_ship,&
          & iarr_ship_V_ship/)
      ivar_234: DO ivar=n0+21,n0+41
         SELECT CASE (ivarid(ivar))
            CASE (iarr_ship_p_overflow,iarr_ship_p_passive,&
                & iarr_ship_p_relocation,iarr_ship_p_suctionhead)
               nrank(ivar) = 2
               nshape(ivar,1:2) = (/nolocs,MaxFractions/)
            CASE DEFAULT
               nrank(ivar) = 1
               nshape(ivar,1) = nolocs
         END SELECT
         SELECT CASE (ivarid(ivar))
            CASE (iarr_ship_t_fill,iarr_ship_t_idle,iarr_ship_t_relocate,&
                & iarr_ship_t_sail_dredge,iarr_ship_t_sail_relocate,&
                & iarr_ship_t_wait)
         END SELECT
      ENDDO ivar_234

!     ---sub-zone dredging
      ivarid(n0+42:n0+45) = (/&
          & iarr_subdred_nr_coordpol,iarr_subdred_nr_subzones,&
          & iarr_subdred_Xpol,iarr_subdred_Ypol/)
      ivar_235: DO ivar=n0+42,n0+45
         SELECT CASE (ivarid(ivar))
            CASE (iarr_subdred_nr_coordpol)
               nrank(ivar) = 2
               nshape(ivar,1:2) = (/nolocs,MaxSubzones/)
            CASE (iarr_subdred_Xpol,iarr_subdred_Ypol)
               nrank(ivar) = 3
               nshape(ivar,1:3) = (/nolocs,MaxCoordPol,MaxSubzones/)
            CASE DEFAULT
               nrank(ivar) = 1
               nshape(ivar,1) = nolocs
         END SELECT
         SELECT CASE (ivarid(ivar))
            CASE (iarr_subdred_nr_coordpol,iarr_subdred_nr_subzones)
               data_type(ivar) = int_type
         END SELECT
      ENDDO ivar_235

!     ---sub-zone relocation
      ivarid(n0+46:n0+49)  =(/&
          & iarr_subrel_nr_coordpol,iarr_subrel_nr_subzones,&
          & iarr_subrel_Xpol,iarr_subrel_Ypol/)
      ivar_236: DO ivar=n0+46,n0+49
         SELECT CASE (ivarid(ivar))
            CASE (iarr_subrel_nr_coordpol)
               nrank(ivar) = 2
               nshape(ivar,1:2) = (/nolocs,MaxSubzones/)
            CASE (iarr_subrel_Xpol,iarr_subrel_Ypol)
               nrank(ivar) = 3
               nshape(ivar,1:3) = (/nolocs,MaxCoordPol,MaxSubzones/)
            CASE DEFAULT
               nrank(ivar) = 1
               nshape(ivar,1) = nolocs
         END SELECT
         SELECT CASE (ivarid(ivar))
            CASE (iarr_subrel_nr_coordpol,iarr_subrel_nr_subzones)
               data_type(ivar) = int_type
         END SELECT
      ENDDO ivar_236
END SELECT

IF (loglev2.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, logexit
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, logexit
   ENDIF
ENDIF
pglev = pglev - 1


RETURN

END SUBROUTINE set_sedvars_atts


END MODULE sedvars_routines
