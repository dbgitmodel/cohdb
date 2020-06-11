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

MODULE modvars_routines
!************************************************************************
!
! *modvars_routines* Utility routines for model variables and file attributes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modvars_routines.f90  V2.11.1
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description -
!
! Reference -
!
! Routines - inquire_var, set_modfiles_atts, set_modfiles_name, set_modvars_atts
!
! Functions - file_suffix
!
!************************************************************************
!
USE iopars
USE syspars


IMPLICIT NONE

CONTAINS

!========================================================================

FUNCTION file_suffix(form)
!************************************************************************
!
! *file_suffix* Returns the file suffix of a model file given its output form
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.7.1
!
! Description -
!
! Module calls - error_abort
!
!************************************************************************
!
!
!*Arguments
!
CHARACTER (LEN=lenform), INTENT(IN) :: form
CHARACTER (LEN=3) :: file_suffix


SELECT CASE (form)
   CASE ('A'); file_suffix = 'txt'
   CASE ('U'); file_suffix = 'bin'
   CASE ('N'); file_suffix = 'nc '
END SELECT


RETURN

END FUNCTION file_suffix

!========================================================================

SUBROUTINE inquire_var(varid,f90_name,standard_name,comment,long_name,units,&
                     & node,vector_name,data_type,fill_value,nrank,global_dims,&
                     & local_dims,halo_dims,numvar,varatts)
!************************************************************************
!
! *inquire_var* Obtain information about a specific model variable
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modvars_routines.f90  V2.11.1
!
! Description -
!
! Module calls - inquire_biovar, inquire_partvar, inquire_sedvar
!
!************************************************************************
!
USE datatypes
USE gridpars
USE modids
USE nestgrids
USE obconds
USE paralpars
USE physpars
USE relaxation
USE structures
USE switches
USE tide
USE wavepars
USE biovars_routines, ONLY: inquire_biovar
USE partvars_routines, ONLY: inquire_partvar
USE sedvars_routines, ONLY: inquire_sedvar

!
!*Arguments
!
CHARACTER (LEN=lenname), INTENT(OUT), OPTIONAL :: f90_name
CHARACTER (LEN=lendesc), INTENT(OUT), OPTIONAL :: comment, long_name, &
                                                & standard_name, vector_name
CHARACTER (LEN=lenunit), INTENT(OUT), OPTIONAL :: units
CHARACTER (LEN=lennode), INTENT(OUT), OPTIONAL :: node
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: numvar
INTEGER, INTENT(OUT), OPTIONAL :: data_type, nrank
TYPE (VariableAtts), INTENT(OUT), OPTIONAL :: varatts
INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(4) :: halo_dims
INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(5) :: global_dims, local_dims
REAL (KIND=kndrlong), INTENT(OUT), OPTIONAL :: fill_value

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
!*numvar*        INTEGER Variable number in case ivarid is a multi-variable
!                        key id
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
CHARACTER (LEN=12) :: cnum
INTEGER :: dtype, nodim
INTEGER,  DIMENSION(4) :: nhdims
INTEGER,  DIMENSION(5) :: ngdims, nldims
REAL :: fvalue
TYPE (VariableAtts) :: varattsx

!
!1. Non-physical parameters
!--------------------------
!

SELECT CASE (varid)
   CASE (MaxModArids+1:MaxModArids+MaxSedArids)
      CALL inquire_sedvar(varid,fname,sname,commt,lname,unit,cnode,vname,&
                        & dtype,fvalue,nodim,ngdims,nldims,nhdims,varattsx)
   CASE (MaxModArids+MaxSedArids+1:MaxModArids+MaxSedArids+MaxBioArids)
      CALL inquire_biovar(varid,fname,sname,commt,lname,unit,cnode,vname,&
           & dtype,fvalue,nodim,ngdims,nldims,nhdims,varattsx)
   CASE (MaxModArids+MaxSedArids+MaxBioArids+1:MaxTotArids)
      CALL inquire_partvar(varid,fname,sname,commt,lname,unit,cnode,vname,&
           & dtype,fvalue,nodim,ngdims,nldims,nhdims,varattsx)
END SELECT

IF (varid.GT.MaxModArids.AND.varid.LE.MaxTotArids) THEN
 
   IF (PRESENT(f90_name)) THEN
      IF (PRESENT(numvar)) THEN
         WRITE (cnum,'(I12)') numvar
         cnum = ADJUSTL(cnum)
         f90_name = TRIM(fname)//'_'//TRIM(cnum)
      ELSE
         f90_name = fname
      ENDIF
   ENDIF

   IF (PRESENT(standard_name)) standard_name = sname
   IF (PRESENT(comment)) comment = commt
   IF (PRESENT(long_name)) long_name = lname
   IF (PRESENT(units)) units = unit
   IF (PRESENT(node)) node = cnode
   IF (PRESENT(vector_name)) vector_name = vname
   IF (PRESENT(data_type)) data_type = dtype
   IF (PRESENT(fill_value)) fill_value = fvalue
   IF (PRESENT(nrank)) nrank = nodim
   IF (PRESENT(global_dims)) global_dims = ngdims
   IF (PRESENT(halo_dims)) halo_dims = nhdims
   IF (PRESENT(local_dims)) local_dims = nldims
   IF (PRESENT(varatts)) varatts = varattsx

   RETURN

ENDIF

!
!2. Initialise
!-------------
!

unit = '_'; cnode = 'C'; vname = ''
dtype = float_type; nodim = 0
ngdims = 0; nhdims = 0; nldims = 0
fvalue = float_fill
commt = 'standard name constructed following the netCDF CF-1.6 guidelines'

!
!3. Get model array info
!-----------------------
!

SELECT CASE (varid)

!
!3.1 Model grid
!--------------
!

CASE (iarr_alphatc_fld)
   sname = 'factor_for_flooding_drying_scheme'
   fname = 'alphatc_fld'
   lname = 'Drying/wetting factor at C-nodes'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_alphatu_fld)
   sname = 'factor_for_flooding_drying_scheme_at_u_nodes' 
   fname = 'alphatu_fld'
   lname = 'Drying/wetting factor at U-nodes'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_alphatv_fld)
   sname = 'factor_for_flooding_drying_scheme_at_v_nodes'
   fname = 'alphatv_fld'
   lname = 'Drying/wetting factor at V-nodes'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_coriolatu)
   sname = 'coriolis_frequency_at_u_nodes'
   fname = 'coriolatu'
   lname = 'Coriolis frequency at U-nodes'
   unit = 'radian s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,1,0/)
CASE (iarr_coriolatv)
   sname = 'coriolis_frequency_at_v_nodes'
   fname = 'coriolatv'
   lname = 'Coriolis frequency at V-nodes'
   unit = 'radian s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,0,1/)
CASE (iarr_delxatc)
   sname = 'magnitude_of_derivative_of_position_wrt_x_coordinate_index'
   commt = ''
   fname = 'delxatc'
   lname = 'X-spacing at C-nodes'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_delxatu)
   sname = 'magnitude_of_derivative_of_position_wrt_x_coordinate_index_'//&
         & 'at_u_nodes' 
   fname = 'delxatu'
   lname = 'X-spacing at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_delxatuv)
   sname = 'magnitude_of_derivative_of_position_wrt_x_coordinate_index_'//&
         & 'at_uv_nodes'
   fname = 'delxatuv'
   lname = 'X-spacing at UV-nodes'
   unit = 'm'
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_delxatv)
   sname = 'magnitude_of_derivative_of_position_wrt_x_coordinate_index_'//&
         & 'at_v_nodes' 
   fname = 'delxatv'
   lname = 'X-spacing at V-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_delyatc)
   sname = 'magnitude_of_derivative_of_position_wrt_y_coordinate_index'
   commt = ''
   fname = 'delyatc'
   lname = 'Y-spacing at C-nodes'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_delyatu)
   sname = 'magnitude_of_derivative_of_position_wrt_y_coordinate_index_'//&
         & 'at_u_nodes' 
   fname = 'delyatu'
   lname = 'Y-spacing at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_delyatuv)
   sname = 'magnitude_of_derivative_of_position_wrt_y_coordinate_index_'//&
         & 'at_uv_nodes' 
   fname = 'delyatuv'
   lname = 'Y-spacing at UV-nodes'
   unit = 'm'
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_delyatv)
   sname = 'magnitude_of_derivative_of_position_wrt_y_coordinate_index_'//&
         & 'at_v_nodes' 
   fname = 'delyatv'
   lname = 'Y-spacing at V-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   nhdims = nhalo
CASE (iarr_delzatc)
   sname = 'magnitude_of_derivative_of_position_wrt_model_level_number'
   commt = ''
   fname = 'delzatc'
   lname = 'Vertical grid spacing at C-nodes'
   unit = 'm'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
CASE (iarr_delzatu)
   sname = 'magnitude_of_derivative_of_position_wrt_model_level_number_'//&
         & 'at_u_nodes' 
   fname = 'delzatu'
   lname = 'Vertical grid spacing at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,1,0,0/)
CASE (iarr_delzatuv)
   sname = 'magnitude_of_derivative_of_position_wrt_model_level_number_'//&
         & 'at_uv_nodes' 
   fname = 'delzatuv'
   lname = 'Vertical grid spacing at UV-nodes'
   unit = 'm'
   cnode = 'UV'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,1,0,1/)
CASE (iarr_delzatuw)
   sname = 'magnitude_of_derivative_of_position_wrt_model_level_number_'//&
         & 'at_uw_nodes' 
   fname = 'delzatuw'
   lname = 'Vertical grid spacing at UW-nodes'
   unit = 'm'
   cnode = 'UW'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = (/0,1,0,0/)
CASE (iarr_delzatv)
   sname = 'magnitude_of_derivative_of_position_wrt_model_level_number_'//&
         & 'at_v_nodes' 
   fname = 'delzatv'
   lname = 'Vertical grid spacing at V-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,0,1,1/)
CASE (iarr_delzatvw)
   sname = 'magnitude_of_derivative_of_position_wrt_model_level_number_'//&
         & 'at_vw_nodes' 
   fname = 'delzatvw'
   lname = 'Vertical grid spacing at VW-nodes'
   unit = 'm'
   cnode = 'VW'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = (/0,0,0,1/)
CASE (iarr_delzatw)
   sname = 'magnitude_of_derivative_of_position_wrt_model_level_number_'//&
         & 'at_w_nodes' 
   fname = 'delzatw'
   lname = 'Vertical grid spacing at W-nodes'
   unit = 'm'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = 1
CASE (iarr_dryfac)
   sname = 'dry_area_fraction'
   fname = 'dryfac'
   lname = 'Dry area fraction'
   unit = '1'
   cnode = ''
   nodim = 0
CASE (iarr_gaccatc)
   sname = 'acceleration_at_c_nodes_due_to_gravity' 
   fname = 'gaccatc'
   lname = 'Acceleration of gravity at C-nodes'
   unit = 'm2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_gaccatu)
   sname = 'acceleration_at_u_nodes_due_to_gravity' 
   fname = 'gaccatu'
   lname = 'Acceleration of gravity at U-nodes'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,0/)
CASE (iarr_gaccatv)
   sname = 'acceleration_at_v_nodes_due_to_gravity' 
   fname = 'gaccatv'
   lname = 'Acceleration of gravity at V-nodes'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,0,1/)
CASE (iarr_gangleatc)
   sname = 'angle_of_model_grid_wrt_reference_grid'
   fname = 'gangleatc'
   lname = 'Grid angle at C-nodes'
   unit = 'radian'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_garea)
   sname = 'horizontal_cell_area'
   fname = 'garea'
   lname = 'Horizontal grid cell area'
   unit = 'm2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_gdelxglb)
   sname = 'rectangular_grid_spacings_in_x_direction'
   fname = 'gdelxglb'
   lname = 'Global grid spacings in X-direction for rectangular grid' 
   unit = MERGE('m          ','degree_east',iopt_grid_sph.EQ.0)
   nodim = 1
   ngdims(1) = nc
   nhdims(1:2) = 1
CASE (iarr_gdelyglb)
   sname = 'rectangular_grid_spacings_in_y_direction'
   fname = 'gdelyglb'
   lname = 'Global grid spacings in Y-direction for rectangular grid' 
   unit = MERGE('m           ','degree_north',iopt_grid_sph.EQ.0)
   nodim = 1
   ngdims(1) = nr
   nhdims(1:2) = 1
CASE (iarr_gscoordatc)
   sname = 'ocean_sigma_coordinate_at_c_nodes'
   fname = 'gscoordatc'
   lname = 'Sigma coordinates at C-nodes'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
CASE (iarr_gscoordatu)
   sname = 'ocean_sigma_coordinate_at_u_nodes' 
   fname = 'gscoordatu'
   lname = 'Sigma coordinates at U-nodes'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
CASE (iarr_gscoordatuvw)
   sname = 'local_ocean_sigma_coordinate_at_uvw_nodes' 
   fname = 'gscoordatuvw'
   lname = 'Local sigma coordinates at W-nodes'
   cnode = 'UVW'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = 1
CASE (iarr_gscoordatuw)
   sname = 'ocean_sigma_coordinate_at_uw_nodes' 
   fname = 'gscoordatuw'
   lname = 'Sigma coordinates at UW-nodes'
   cnode = 'UW'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = 1
CASE (iarr_gscoordatv)
   sname = 'ocean_sigma_coordinate_at_v_nodes' 
   fname = 'gscoordatv'
   lname = 'Sigma coordinates at V-nodes'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
CASE (iarr_gscoordatvw)
   sname = 'ocean_sigma_coordinate_at_vw_nodes' 
   fname = 'gscoordatvw'
   lname = 'Sigma coordinates at VW-nodes'
   cnode = 'VW'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = 1
CASE (iarr_gscoordatw)
   sname = 'ocean_sigma_coordinate_at_w_nodes' 
   fname = 'gscoordatw'
   lname = 'Sigma coordinates at W-nodes'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = 1
CASE (iarr_gscoordglb)
   sname = 'global_ocean_sigma_coordinate_at_w_nodes'
   fname = 'gscoordglb'
   lname = 'Global sigma coordinates at W-nodes'
   cnode = 'W'
   nodim = 3
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = 1
CASE (iarr_gsigcoordatc)
   sname = 'vertical_ocean_sigma_coordinate'
   fname = 'gsigcoordatc'
   lname = 'Sigma coordinates at centered nodes on uniform grid'
   cnode = 'C'
   nodim = 1
   ngdims(1) = nz
   CASE (iarr_gsigcoordatw)
   sname = 'vertical_ocean_sigma_coordinate'
   fname = 'gsigcoordatw'
   lname = 'Sigma coordinates at W-nodes on uniform grid'
   cnode = 'W'
   nodim = 1
   ngdims(1) = nz+1
CASE (iarr_gxcoord)
   sname = 'local_x_coordinate_at_uv_nodes' 
   fname = 'gxcoord'
   lname = 'Local X-coordinates at UV-nodes'
   unit = MERGE('m          ','degree_east',iopt_grid_sph.EQ.0)
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_gxcoordglb)
   sname = 'global_x_coordinate_at_uv_nodes' 
   fname = 'gxcoordglb'
   lname = 'Global X-coordinates at UV-nodes'
   unit = MERGE('m          ','degree_east',iopt_grid_sph.EQ.0)
   cnode = 'UV'
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_gxlon)
   sname = 'longitude'
   commt = ''
   fname = 'gxlon'
   lname = 'Longitude at C-nodes'
   unit = 'degree_east'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_gycoord)
   sname = 'local_y_coordinate_at_uv_nodes' 
   fname = 'gycoord'
   lname = 'Local Y-coordinates at UV-nodes'
   unit = MERGE('m           ','degree_north',iopt_grid_sph.EQ.0)
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_gycoordglb)
   sname = 'global_y_coordinate_at_uv_nodes' 
   fname = 'gycoordglb'
   lname = 'Global Y-coordinates at UV-nodes'
   unit = MERGE('m           ','degree_north',iopt_grid_sph.EQ.0)
   cnode = 'UV'
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_gylat)
   sname = 'latitude' 
   commt = ''
   fname = 'gylat'
   lname = 'Latitude at C-nodes'
   unit = 'degree_north'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_indexobu)
   sname = 'index_at_u_open_boundaries_for_local_to_global_mapping'
   fname = 'indexobu'
   lname = 'Index at U-open_boundaries for local to global mapping'
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = nobu
CASE (iarr_indexobuprocs)
   sname = 'index_at_u_open_boundaries_for_local_to_global_mapping_per_process'
   fname = 'indexobuprocs'
   lname = 'Index at U-open boundaries for local to global mapping per process'
   dtype = int_type
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,nprocs/)
CASE (iarr_indexobv)
   sname = 'index_at_v_open_boundaries_for_local_to_global_mapping'
   fname = 'indexobv'
   lname = 'Index at V-open_boundaries for local to global mapping'
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = nobv
CASE (iarr_indexobvprocs)
   sname = 'index_at_v_open_boundaries_for_local_to_global_mapping_per_process'
   fname = 'indexobvprocs'
   lname = 'Index at V-open boundaries for local to global mapping per process'
   dtype = int_type
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,nprocs/)
CASE (iarr_indexobx)
   sname = 'index_at_x_open_boundaries_for_local_to_global_mapping'
   fname = 'indexobx'
   lname = 'Index at X-open_boundaries for local to global mapping'
   dtype = int_type
   cnode = 'X'
   nodim = 1
   ngdims(1) = nobx
CASE (iarr_indexobxprocs)
   sname = 'index_at_x_open_boundaries_for_local_to_global_mapping_per_process'
   fname = 'indexobxprocs'
   lname = 'Index at V-open boundaries for local to global mapping per process'
   dtype = int_type
   cnode = 'X'
   nodim = 2
   ngdims(1:2) = (/nobx,nprocs/)
CASE (iarr_indexoby)
   sname = 'index_at_y_open_boundaries_for_local_to_global_mapping'
   fname = 'indexoby'
   lname = 'Index at Y-open_boundaries for local to global mapping'
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   ngdims(1) = noby
CASE (iarr_indexobyprocs)
   sname = 'index_at_y_open_boundaries_for_local_to_global_mapping_per_process'
   fname = 'indexobyprocs'
   lname = 'Index at Y-open boundaries for local to global mapping per process'
   dtype = int_type
   cnode = 'Y'
   nodim = 2
   ngdims(1:2) = (/noby,nprocs/)
CASE (iarr_iobu)
   sname = 'global_x_index_at_u_open_boundaries'
   fname = 'iobu'
   lname = 'Global X-index of U-open boundaries' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1:1) = nobu
CASE (iarr_iobuloc)
   sname = 'local_x_index_at_u_open_boundaries'
   fname = 'iobuloc'
   lname = 'Local X-index of U-open boundaries' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   nldims(1:1) = nobuloc_ext
CASE (iarr_iobv)
   sname = 'global_x_index_at_v_open_boundaries'
   fname = 'iobv'
   lname = 'Global X-index of V-open boundaries' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1:1) = nobv
CASE (iarr_iobvloc)
   sname = 'local_x_index_at_v_open_boundaries'
   fname = 'iobvloc'
   lname = 'Local X-index of V-open boundaries' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   nldims(1:1) = nobvloc_ext
CASE (iarr_iobx)
   sname = 'global_x_index_at_x_open_boundaries'
   fname = 'iobx'
   lname = 'Global X-index of X-open boundaries' 
   dtype = int_type
   cnode = 'X'
   nodim = 1
   ngdims(1:1) = nobx
CASE (iarr_iobxloc)
   sname = 'local_x_index_at_x_open_boundaries'
   fname = 'iobxloc'
   lname = 'Local X-index of X-open boundaries' 
   dtype = int_type
   cnode = 'X'
   nodim = 1
   nldims(1:1) = nobxloc_ext
CASE (iarr_ioby)
   sname = 'global_x_index_at_y_open_boundaries'
   fname = 'ioby'
   lname = 'Global X-index of Y-open boundaries' 
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   ngdims(1:1) = noby
CASE (iarr_iobyloc)
   sname = 'local_x_index_at_y_open_boundaries'
   fname = 'iobyloc'
   lname = 'Local X-index of Y-open boundaries' 
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   nldims(1:1) = nobyloc_ext
CASE (iarr_jobu)
   sname = 'global_y_index_at_u_open_boundaries'
   fname = 'jobu'
   lname = 'Global Y-index of U-open boundaries' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1:1) = nobu
CASE (iarr_jobuloc)
   sname = 'local_y_index_at_u_open_boundaries'
   fname = 'jobuloc'
   lname = 'Local Y-index of U-open boundaries' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   nldims(1:1) = nobuloc_ext
CASE (iarr_jobv)
   sname = 'global_y_index_at_v_open_boundaries'
   fname = 'jobv'
   lname = 'Global Y-index of V-open boundaries' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1:1) = nobv
CASE (iarr_jobvloc)
   sname = 'local_y_index_at_v_open_boundaries'
   fname = 'jobvloc'
   lname = 'Local Y-index of V-open boundaries' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   nldims(1:1) = nobvloc_ext
CASE (iarr_jobx)
   sname = 'global_y_index_at_x_open_boundaries'
   fname = 'jobx'
   lname = 'Global Y-index of X-open boundaries' 
   dtype = int_type
   cnode = 'X'
   nodim = 1
   ngdims(1:1) = nobx
CASE (iarr_jobxloc)
   sname = 'local_y_index_at_x_open_boundaries'
   fname = 'jobxloc'
   lname = 'Local Y-index of X-open boundaries' 
   dtype = int_type
   cnode = 'X'
   nodim = 1
   nldims(1:1) = nobxloc_ext
CASE (iarr_joby)
   sname = 'global_y_index_at_y_open_boundaries'
   fname = 'joby'
   lname = 'Global Y-index of Y-open boundaries' 
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   ngdims(1:1) = noby
CASE (iarr_jobyloc)
   sname = 'local_y_index_at_y_open_boundaries'
   fname = 'jobyloc'
   lname = 'Local Y-index of Y-open boundaries' 
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   nldims(1:1) = nobyloc_ext
CASE (iarr_maskatc_int)
   sname = 'Local dry_mask_at_c_nodes'
   fname = 'maskatc_int'
   lname = 'Local dry mask at C-nodes'
   unit = 'log'
   dtype = log_type
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 0   
CASE (iarr_mg_nc1procs)
   sname = 'global_x_index_of_local_sw_corner_per_process'
   fname = 'mgvars%nc1procs'
   lname = 'Global X-index of local SW corner per process' 
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nprocs,nomglevels/)
CASE (iarr_mg_nc2procs)
   sname = 'global_x_index_of_local_se_corner_per_process'
   fname = 'mgvars%nc2procs'
   lname = 'Global X-index of local SE corner per process' 
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nprocs,nomglevels/)
CASE (iarr_mg_nr1procs)
   sname = 'global_y_index_of_local_sw_corner_per_process'
   fname = 'mgvars%nr1procs'
   lname = 'Global Y-index of local SW corner per process' 
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nprocs,nomglevels/)
CASE (iarr_mg_nr2procs)
   sname = 'global_y_index_of_local_se_corner_per_process'
   fname = 'mgvars%nr2procs'
   lname = 'Global Y-index of local SE corner per process' 
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nprocs,nomglevels/)
CASE (iarr_ncprocs)
   sname = 'x_local_grid_dimension_per_process'
   fname = 'ncprocs'
   lname = 'Local grid dimension in X-direction per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs   
CASE (iarr_nobuprocs)
   sname = 'number_of_local_interior_u_open_boundaries_per_process' 
   fname = 'nobuprocs'
   lname = 'Number of local interior U-open boundaries'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nobvprocs)
   sname = 'number_of_local_interior_v_open_boundaries_per_process' 
   fname = 'nobvprocs'
   lname = 'Number of local interior V-open boundaries per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nobxprocs)
   sname = 'number_of_local_interior_x_open_boundaries_per_process'
   fname = 'nobxprocs'
   lname = 'Number of local interior X-open boundaries per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nobyprocs)
   sname = 'number_of_local_interior_y_open_boundaries_per_process' 
   fname = 'nobyprocs'
   lname = 'Number of local interior Y-open boundaries per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nodeatc)
   sname = 'grid_pointer_at_c_nodes'
   fname = 'nodeatc'
   lname = 'Grid pointer at C-nodes'
   dtype = int_type
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_nodeatu)
   sname = 'grid_pointer_at_u_nodes'
   fname = 'nodeatu'
   lname = 'Grid pointer at U-nodes'
   dtype = int_type
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_nodeatuv)
   sname = 'grid_pointer_at_uv_nodes'
   fname = 'nodeatuv'
   lname = 'Grid pointer at UV-nodes'
   dtype = int_type
   cnode = 'UV'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_nodeatuw)
   sname = 'grid_pointer_at_uw_nodes'
   fname = 'nodeatuw'
   lname = 'Grid pointer at UW-nodes'
   dtype = int_type
   cnode = 'UW'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo
CASE (iarr_nodeatv)
   sname = 'grid_pointer_at_v_nodes'
   fname = 'nodeatv'
   lname = 'Grid pointer at V-nodes'
   dtype = int_type
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_nodeatvw)
   sname = 'grid_pointer_at_vw_nodes'
   fname = 'nodeatvw'
   lname = 'Grid pointer at VW-nodes'
   dtype = int_type
   cnode = 'VW'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo
CASE (iarr_node2du)
   sname = '2d_grid_pointer_at_u_nodes'
   fname = 'node2du'
   lname = '2-D grid pointer at U-nodes'
   dtype = int_type
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_node2duv)
   sname = '2d_grid_pointer_at_uv_nodes'
   fname = 'node2duv'
   lname = '2-D grid pointer at UV-nodes'
   dtype = int_type
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_node2dv)
   sname = '2d_grid_pointer_at_v_nodes'
   fname = 'node2dv'
   lname = '2-D grid pointer at V-nodes'
   dtype = int_type
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_nosbuprocs)
   sname = 'number_of_local_interior_u_open_sea_boundaries_per_process' 
   fname = 'nosbuprocs'
   lname = 'Number of local interior U-open sea boundaries per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nosbvprocs)
   sname = 'number_of_local_interior_v_open_sea_boundaries_per_process' 
   fname = 'nosbvprocs'
   lname = 'Number of local interior V-open sea boundaries per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nrprocs)
   sname = 'y_local_grid_dimension_per_process'
   fname = 'nrprocs'
   lname = 'Local grid dimension in Y-direction'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nrvbuprocs)
   sname = 'number_of_local_interior_u_river_boundaries_per_process' 
   fname = 'nrvbuprocs'
   lname = 'Number of local interior U-river boundaries per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nrvbvprocs)
   sname = 'number_of_local_interior_v_river_boundaries_per_process' 
   fname = 'nrvbvprocs'
   lname = 'Number of local interior V-river boundaries per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_rlxobcatu)
   sname = 'relaxation_factor_for_momentum_advection_at_u_nodes'
   fname = 'rlxobcatu'
   lname = 'Relaxation factor for momentum advection at U-nodes'
   cnode ='U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_rlxobcatv)
   sname = 'relaxation_factor_for_momentum_advection_at_v_nodes'
   fname = 'rlxobcatv'
   lname = 'Relaxation factor for momentum advection at V-nodes'
   cnode ='V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_seapoint)
   sname = 'Local land_mask_at_c_nodes'
   fname = 'seapoint'
   lname = 'Local land mask'
   unit = 'log'
   dtype = log_type
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_seapointglb)
   sname = 'global_land_mask_at_c_nodes'
   fname = 'seapointglb'
   lname = 'Global land mask'
   unit = 'log'
   dtype = log_type
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_soutobv)
   sname = 'south_open_boundary_at_v_node'
   fname = 'soutobv'
   lname = 'South open boundary at V-node'
   unit = 'log'
   dtype = log_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = nobv
CASE (iarr_soutoby)
   sname = 'south_open_boundary_at_y_node'
   fname = 'soutoby'
   lname = 'South open boundary at Y-node'
   unit = 'log'
   dtype = log_type
   cnode = 'UV'
   nodim = 1
   ngdims(1) = noby
CASE (iarr_westobu)
   sname = 'west_open_boundary_at_u_node'
   fname = 'westobu'
   lname = 'West open boundary at U-node'
   unit = 'log'
   dtype = log_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = nobu
CASE (iarr_westobx)
   sname = 'west_open_boundary_at_x_node'
   fname = 'westobx'
   lname = 'West open boundary at X-node'
   unit = 'log'
   dtype = log_type
   cnode = 'UV'
   nodim = 1
   ngdims(1) = nobx

!
!3.2 Depths
!----------
!

CASE (iarr_depmeanatc)
   sname = 'sea_floor_depth_below_mean_sea_level'
   commt = ''
   fname = 'depmeanatc'
   lname = 'Mean water depth'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_depmeanatu)
   sname = 'sea_floor_depth_below_mean_sea_level_at_u_nodes'
   fname = 'depmeanatu'
   lname = 'Mean water depth at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,1,0,0/)
CASE (iarr_depmeanatuv)
   sname = 'sea_floor_depth_below_mean_sea_level_at_uv_nodes'
   fname = 'depmeanatuv'
   lname = 'Mean water depth at UV-nodes'
   unit = 'm'
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,1/)
CASE (iarr_depmeanatv)
   sname = 'sea_floor_depth_below_mean_sea_level_at_v_nodes'
   fname = 'depmeanatv'
   lname = 'Mean water depth at V-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,1,1/)
CASE (iarr_depmeanglb)
   sname = 'global_sea_floor_depth_below_mean_sea_level'
   fname = 'depmeanglb'
   lname = 'Global mean water depth'
   nhdims = 1
   unit = 'm'
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_deptotatc)
   sname = 'sea_floor_depth_below_sea_surface'
   commt = ''
   fname = 'deptotatc'
   lname = 'Total water depth'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_deptotatc_err)
   sname = 'error_in_sea_floor_depth_below_sea_surface_due_to_mass_'//&
         & 'conservation_violation'
   fname = 'deptotatc_err'
   lname = 'Total water depth error'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_deptotatc_old)
   sname = 'sea_floor_depth_below_sea_surface_at_old_time_level'
   fname = 'deptotatc_old'
   lname = 'Total water depth at old time level'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_deptotatu)
   sname = 'sea_floor_depth_below_sea_surface_at_u_nodes'
   fname = 'deptotatu'
   lname = 'Total water depth at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/nhalo,nhalo,1,1/)
CASE (iarr_deptotatu_old)
   sname = 'sea_floor_depth_below_sea_surface_at_u_nodes_at_old_time_level'
   fname = 'deptotatu_old'
   lname = 'Total water depth at U-nodes and old time level'
   unit = 'm'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,1,0,0/)
CASE (iarr_deptotatuv)
   sname = 'sea_floor_depth_below_sea_surface_at_uv_nodes'
   commt = ''
   fname = 'deptotatuv'
   lname = 'Total water depth at UV-nodes'
   unit = 'm'
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,1/)
CASE (iarr_deptotatv)
   sname = 'sea_floor_depth_below_sea_surface_at_v_nodes'
   commt = ''
   fname = 'deptotatv'
   lname = 'Total water depth at v-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,1,nhalo,nhalo/)
CASE (iarr_deptotatv_old)
   sname = 'sea_floor_depth_below_sea_surface_at_v_nodes_at_old_time_level'
   fname = 'deptotatv_old'
   lname = 'Total water depth at V-nodes and old time level'
   unit = 'm'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,1,1/)
CASE (iarr_dzeta)
   sname = 'sea_surface_elevation_difference_between_new_and_previous_'//&
         & 'iteration_level'
   fname = 'dzeta'
   lname = 'Surface elevation difference at new and previous iteration level'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_zeta)
   sname = 'sea_surface_height_above_geoid'
   commt = ''
   fname = 'zeta'
   lname = 'Surface elevation'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_zeta_old)
   sname = 'sea_surface_height_above_geoid_at_old_time_level'
   fname = 'zeta_old'
   lname = 'Surface elevation at old time level'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1

!
!3.3 Currents
!------------
!

CASE (iarr_hdvelmag)
   sname = 'integral_wrt_depth_of_horizontal_sea_water_velocity'
   fname = 'hdvelmag'
   lname = 'Magnitude of depth-integrated horizontal current'
   unit = 'm2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_hmvelmag)
   sname = 'depth_mean_magnitude_of_horizontal_sea_water_velocity'
   fname = 'hmvelmag'
   lname = 'Magnitude of depth-mean horizontal current'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_hmvelmag_hadv_cfl)
   sname = 'CFL_number_of_depth_mean_horizontal_current'
   fname = 'hmvelmag_hadv_cfl'
   lname = 'CFL number depth mean horizontal current'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_hvelmag)
   sname = 'magnitude_of_horizontal_sea_water_velocity'
   fname = 'hvelmag'
   lname = 'Magnitude of horizontal current'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_hvelmag_hadv_cfl)
   sname = 'CFL_number_horizontal_sea_water_velocity'
   fname = 'hvelmag_hadv_cfl'
   lname = 'CFL number 3-D horizontal current'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_p2dbcgradatu)
   sname = 'integral_wrt_depth_of_x_derivative_of_kinematic_baroclinic_'//&
         & 'ocean_pressure'
   fname = 'p2dbcgradatu'
   lname = 'X-component of depth-integrated kinematic baroclinic '//&
         & 'pressure gradient'
   vname = 'Depth-integrated kinematic baroclinic pressure gradient'
   unit = 'm2 s-2'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_p2dbcgradatv)
   sname = 'integral_wrt_depth_of_y_derivative_of_kinematic_baroclinic_'//&
         & 'ocean_pressure'
   fname = 'p2dbcgradatv'
   lname = 'Y-component of depth-integrated kinematic baroclinic '//&
         & 'pressure gradient'
   vname = 'Depth-integrated kinematic baroclinic pressure gradient'
   unit = 'm2 s-2'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_p3dbcgradatu)
   sname = 'x_derivative_of_kinematic_baroclinic_ocean_pressure'
   fname = 'p3dbcgradatu'
   lname = 'X-component of kinematic baroclinic pressure gradient'
   vname = 'Kinematic baroclonic pressure gradient'
   unit = 'm s-2'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_p3dbcgradatv)
   sname = 'y_derivative_of_kinematic_baroclinic_ocean_pressure'
   fname = 'p3dbcgradatv'
   lname = 'Y-component of kinematic baroclinic pressure gradient'
   vname = 'Kinematic baroclonic pressure gradient'
   unit = 'm s-2'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_udevint)
   sname = 'baroclinic_term_in_U_momentum_equation' 
   fname = 'udevint'
   lname = 'Baroclinic term in U-momentum equation'
   vname = 'Baroclinic term in 2-D momentum equations'
   unit = 'm2 s-2'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_udfvel)
   sname = 'integral_wrt_depth_of_x_sea_water_filtered_velocity'
   fname = 'udfvel'
   lname = 'X-component of filtered depth-integrated current'
   vname = 'Filtered depth-integrated current'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,0/)
CASE (iarr_udvel)
   sname = 'integral_wrt_depth_of_x_sea_water_velocity'
   fname = 'udvel'
   lname = 'X-component of depth-integrated current'
   vname = 'Depth-integrated current'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_udvel_old)
   sname = 'integral_wrt_depth_of_x_sea_water_velocity_at_old_time_level'
   fname = 'udvel_old'
   lname = 'X-component of depth-integrated current at old time level'
   vname = 'Depth-integrated current at old time level'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_ufvel)
   sname = 'x_sea_water_filtered_velocity'
   fname = 'ufvel'
   lname = 'X-component of filtered current'
   vname = 'Filtered current'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/nhalo-1,nhalo,0,0/)
CASE (iarr_umpred)
   sname = 'depth_mean_predicted_x_sea_water_velocity'
   fname = 'umpred'
   lname = 'X-component of depth-mean predicted current'
   vname = 'Depth-mean predicted current'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,1,0,0/)
CASE (iarr_umvel)
   sname = 'depth_mean_x_sea_water_velocity'
   fname = 'umvel'
   lname = 'X-component of depth-mean current'
   vname = 'Depth-mean current'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/nhalo,nhalo,1,1/)
CASE (iarr_umvel_hadv_cfl)
   sname = 'CFL_number_of_depth_mean_x_sea_water_velocity'
   fname = 'umvel_hadv_cfl'
   lname = 'CFL number X-component of depth-mean current'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/nhalo,nhalo,1,1/)
CASE (iarr_umvel_old)
   sname = 'depth_mean_x_sea_water_velocity_at_old_time_level'
   fname = 'umvel_old'
   lname = 'X-component of depth-mean current at old time level'
   vname = 'Depth-mean current'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/nhalo,nhalo,1,1/)
CASE (iarr_uvel)
   sname = 'x_sea_water_velocity'
   commt = ''
   fname = 'uvel'
   lname = 'X-component of current'
   vname = 'Current'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_uvel_hadv_cfl)
   sname = 'CFL_number_of_x_sea_water_velocity'
   fname = 'uvel_hadv_cfl'
   lname = 'CFL number X-component of 3-D current'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_uvel_old)
   sname = 'x_sea_water_velocity_at_old_time_level'
   fname = 'uvel_old'
   lname = 'X-component of current at old time level'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_vdevint)
   sname = 'baroclinic_term_in_V_momentum_equation' 
   fname = 'vdevint'
   lname = 'Baroclinic term in V-momentum equation'
   vname = 'Baroclinic term in 2-D momentum equations'
   unit = 'm2 s-2'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_vdfvel)
   sname = 'integral_wrt_depth_of_y_sea_water_filtered_velocity'
   fname = 'vdfvel'
   lname = 'Y-component of filtered depth-integrated current'
   vname = 'Filtered depth-integrated current'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,0,1/)
CASE (iarr_vdvel)
   sname = 'integral_wrt_depth_of_y_sea_water_velocity'
   fname = 'vdvel'
   lname = 'Y-component of depth-integrated current'
   vname = 'Depth-integrated current'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_vdvel_old)
   sname = 'integral_wrt_depth_of_y_sea_water_velocity_at_old_time_level'
   fname = 'vdvel_old'
   lname = 'Y-component of depth-integrated current at old time level'
   vname = 'Depth-integrated current at old time level'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_vel2d)
   sname = 'integral_wrt_depth_of_sea_water_velocity_at_nested_boundaries'  
   fname = 'vel2d'
   lname = 'Depth-integrated current at nested boundaries'
   unit = 'm2 s-1'
   cnode = ''
   nodim = 2
CASE (iarr_vel3d)
   sname = 'baroclinic_sea_water_velocity_at_nested_boundaries' 
   fname = 'vel3d'
   lname = 'Baroclinic current at nested boundaries'
   unit = 'm s-1'
   cnode = ''
   nodim = 3
CASE (iarr_vfvel)
   sname = 'y_sea_water_filtered_velocity'
   fname = 'vfvel'
   lname = 'Y-component of filtered current'
   vname = 'Filtered current'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,0,nhalo-1,nhalo/)
CASE (iarr_vmpred)
   sname = 'depth_mean_predicted_x_sea_water_velocity'
   fname = 'vmpred'
   lname = 'Y-component of depth-mean predicted current'
   vname = 'Depth-mean predicted current'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,1,1/)
CASE (iarr_vmvel)
   sname = 'depth_mean_y_sea_water_velocity'
   fname = 'vmvel'
   lname = 'Y-component of depth-mean current'
   vname = 'Depth-mean current'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,1,nhalo,nhalo/)
CASE (iarr_vmvel_hadv_cfl)
   sname = 'CFL_number_of_depth_mean_y_sea_water_velocity'
   fname = 'vmvel_hadv_cfl'
   lname = 'CFL number Y-component of depth-mean current'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/nhalo,nhalo,1,1/)   
CASE (iarr_vmvel_old)
   sname = 'depth_mean_y_sea_water_velocity_at_old_time_level' 
   fname = 'vmvel_old'
   lname = 'Y-component of depth-mean current at old time level'
   vname = 'Depth-mean current at old time level'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,1,nhalo,nhalo/)
CASE (iarr_vvel)
   sname = 'y_sea_water_velocity'
   commt = ''
   fname = 'vvel'
   lname = 'Y-component of current'
   vname = 'Current'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_vvel_hadv_cfl)
   sname = 'CFL_number_of_y_sea_water_velocity'
   fname = 'vvel_cfl'
   lname = 'CFL number Y-component of 3-D current'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo   
CASE (iarr_vvel_old)
   sname = 'y_sea_water_velocity_at_old_time_level'
   fname = 'vvel_old'
   lname = 'Y-component of current at old time level'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_wphys)
   sname = 'physical_upward_sea_water_velocity'
   fname = 'wphys'
   lname = 'Physical vertical velocity'
   vname = 'Current'
   unit = 'm s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_wvel)
   sname = 'transformed_upward_sea_water_velocity'
   fname = 'wvel'
   lname = 'Transformed vertical velocity'
   vname = 'Current'
   unit = 'm s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = (/1,0,1,0/)
CASE (iarr_wvel_vadv_cfl)
   sname = 'CFL_number_transformed_upward_sea_water_velocity'
   fname = 'wvel_vadv_cfl'
   lname = 'CFL number transformed vertical velocity'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
   
!
!3.4 Density
!-----------
!

CASE (iarr_beta_sal)
   sname = 'salinity_coefficient_for_sea_water_expansion'
   fname = 'beta_sal'
   lname = 'Salinity expansion coefficient'
   unit = 'PSU-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
CASE (iarr_beta_temp)
   sname = 'temperature_coefficient_for_sea_water_expansion'
   fname = 'beta_temp'
   lname = 'Temperature expansion coefficient'
   unit = 'degC-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
CASE (iarr_dens)
   sname = 'sea_water_density'
   commt = ''
   fname = 'dens'
   lname = 'Mass density'
   unit = 'kg m-3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
CASE (iarr_sal)
   sname = 'sea_water_salinity'
   commt = ''
   fname = 'sal'
   lname = 'Salinity'
   unit = 'PSU'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_temp)
   sname = 'sea_water_temperature'
   commt = ''
   fname = 'temp'
   lname = 'Temperature'
   unit = 'degC'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo

!
!3.5 Diffusion coefficients
!--------------------------
!

CASE (iarr_hdifcoef2datc)
   sname = 'depth_mean_ocean_momentum_diffusivity'
   fname = 'hdifcoef2datc'
   lname = 'Depth-mean diffusion coefficient at C-nodes'
   unit = 'm2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_hdifcoef2d_mom)
   sname = 'depth_mean_ocean_momentum_diffusivity'
   fname = 'hdifcoef2d_mom'
   lname = 'Depth-mean momentum diffusion coefficient at C-nodes'
   unit = 'm2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
   CASE (iarr_hdifcoef2d_scal)
   sname = 'depth_mean_ocean_scalar_diffusivity'
   fname = 'hdifcoef2d_scal'
   lname = 'Depth-mean scalar diffusion coefficient at C-nodes'
   unit = 'm2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_hdifcoef2datuv)
   sname = 'integral_wrt_depth_of_ocean_diffusivity_at_uv_nodes'
   fname = 'hdifcoef2datuv'
   lname = 'Depth-integrated diffusion coefficient at UV-nodes'
   unit = 'm3 s-1'
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,1/)
CASE (iarr_hdifcoef3datc)
   sname = 'ocean_horizontal_diffusivity'
   fname = 'hdifcoef3datc'
   lname = 'Horizontal 3-D diffusion coefficient at C-nodes'
   unit = 'm2 s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,0,1,0/)
   CASE (iarr_hdifcoef3d_mom)
   sname = 'ocean_horizontal_momentum_diffusivity'
   fname = 'hdifcoef3d_mom'
   lname = 'Horizontal 3-D momentum diffusion coefficient at C-nodes'
   unit = 'm2 s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,0,1,0/)
CASE (iarr_hdifcoef3d_scal)
   sname = 'ocean_horizontal_scalar_diffusivity'
   fname = 'hdifcoef3d_mom'
   lname = 'Horizontal 3-D scalar diffusion coefficient at C-nodes'
   unit = 'm2 s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,0,1,0/)
CASE (iarr_hdifcoef3datu)
   sname = 'ocean_horizontal_diffusivity_at_u_nodes'
   fname = 'hdifcoef3datu'
   lname = 'Horizontal 3-D diffusion coefficient at U-nodes'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,1,0,0/)
CASE (iarr_hdifcoef3datuv)
   sname = 'ocean_horizontal_diffusivity_at_uv_nodes' 
   fname = 'hdifcoef3datuv'
   lname = 'Horizontal 3-D diffusion coefficient at UV-nodes'
   unit = 'm2 s-1'
   cnode = 'UV'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,1,0,1/)
CASE (iarr_hdifcoef3datv)
   sname = 'ocean_horizontal_diffusivity_at_v_nodes'
   fname = 'hdifcoef3datv'
   lname = 'Horizontal 3-D diffusion coefficient at V-nodes'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,0,0,1/)
CASE (iarr_hdifcoef3datw)
   sname = 'ocean_horizontal_diffusivity_at_w_nodes'
   fname = 'hdifcoef3datv'
   lname = 'Horizontal 3-D diffusion coefficient at W-nodes'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 0
CASE (iarr_kinvisc)
   sname = 'sea_water_kinematic_viscosity'
   fname = 'kinvisc'
   lname = 'Kinematic viscosity'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = (/1,0,1,0/)
CASE (iarr_mom_vdif_cfl)
   sname = 'CFL_number_of_vertical_momentum_dfiffusion'
   fname = 'mom_vdif_cfl'
   lname = 'CFL number for vertical momentum diffusion'
   cnode = 'C'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_scal_vdif_cfl)
   sname = 'CFL_number_of_vertical_scalar_diffusion'
   fname = 'vvel_cfl'
   lname = 'CFL number for vertical scalar diffusion'
   cnode = 'C'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo   
CASE (iarr_vdifcoefmom)
   sname = 'sea_water_turbulent_viscosity'
   fname = 'vdifcoefmom'
   lname = 'Eddy viscosity'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = (/1,0,1,0/)
CASE (iarr_vdifcoefscal)
   sname = 'sea_water_turbulent_diffusivity'
   fname = 'vdifcoefscal'
   lname = 'Eddy diffusivity'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_vdifcoefscal_norot)
   sname = 'non_rotated_component_sea_water_turbulent_diffusivity'
   fname = 'vdifcoefscal_norot'
   lname = 'Non-rotated component of eddy diffusivity'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_vdifcoefscal_rot)
   sname = 'rotated_component_sea_water_turbulent_diffusivity'
   fname = 'vdifcoefscal_rot'
   lname = 'Rotated component of eddy diffusivity'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_vdifcoeftke)
   sname = 'vertical_diffusion_coefficient_for_turbulence_energy'
   fname = 'vdifcoeftke'
   lname = 'Vertical diffusion coefficient for turbulence energy'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_xslopeatu_geo)
   sname = 'isolevel_slopes_in_x_direction_with_respect_sigma_surfaces_at_u_'//&
         & 'nodes'
   fname = 'xslopeatu_geo'
   lname = 'Isolevel slopes in the X-direction at U-nodes'
   cnode  = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,1,0,0/)
CASE (iarr_xslopeatu_ziso)
   sname = 'isoneutral_slopes_in_x_direction_at_u_nodes'
   fname = 'xslopeatu_iso'
   lname = 'Isolevel slopes in the X-direction at U-nodes'
   cnode  = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,1,0,0/)
CASE (iarr_xslopeatw_geo)
   sname = 'isoneutral_slopes_in_x_direction_with_respect_sigma_surfaces_at_'//&
         & 'w_nodes'
   fname = 'xslopeatw_geo'
   lname = 'Isolevel slopes in the X-direction at W-nodes'
   cnode  = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 0
CASE (iarr_yslopeatv_geo)
   sname = 'isolevel_slopes_in_y_direction_with_respect_sigma_surfaces_at_'//&
         & 'v_nodes'
   fname = 'yslopeatv_geo'
   lname = 'Isolevel slopes in the Y-direction at V-nodes'
   cnode  = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,0,0,1/)
CASE (iarr_yslopeatv_ziso)
   sname = 'isoneutral_slopes_in_Y_direction_at_v_nodes'
   fname = 'yslopeatv_ziso'
   lname = 'Isoneutral slopes in the Y-direction at V-nodes'
   cnode  = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,0,1,1/)
CASE (iarr_yslopeatw_geo)
   sname = 'isoneutral_slopes_in_y_direction_with_respect_sigma_surfaces_at_'//&
         & 'w_nodes'
   fname = 'yslopeatw_geo'
   lname = 'Isolevel slopes in the Y-direction at W-nodes'
   cnode  = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_2D_hdif_cfl)
   sname = 'CFL_number_of_horizontal_dfiffusion'
   fname = '2D_hdif_cfl'
   lname = 'CFL number for horizontal diffusion'
   cnode = 'C'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)

!
!3.6 Turbulence
!--------------
!

CASE (iarr_buofreq2)
   sname = 'square_of_buoyancy_frequency_in_water'
   fname = 'buofreq2'
   lname = 'Squared buoyancy frequency'
   unit = 's-2'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_dissip)
   sname = 'kinetic_energy_dissipation_in_sea_water'
   fname = 'dissip'
   lname = 'Dissipation of turbulence kinetic energy'
   unit = 'W kg-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo
CASE (iarr_ricnum)
   sname = 'richardson_number_in_sea_water'
   commt = ''
   fname = 'ricnum'
   lname = 'Richardson number'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_shearfreq2)
   sname = 'square_of_shear_frequency_in_sea_water'
   fname = 'shearfreq2'
   lname = 'Squared shear frequency'
   unit = 's-2'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_tke)
   sname = 'turbulent_kinetic_energy_in_sea_water'
   fname = 'tke'
   lname = 'Turbulent kinetic energy'
   unit = 'J kg-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo
CASE (iarr_tke_old)
   sname = 'turbulent_kinetic_energy_in_sea_water_at_old_time_level'
   fname = 'tke_old'
   lname = 'Turbulent kinetic energy at old time level'
   unit = 'J kg-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo
CASE (iarr_tkezl)
   sname = 'product_of_turbulent_kinetic_energy_and_mixing_length'
   fname = 'tkezl'
   lname = 'Turbulent kinetic energy times mixing length'
   unit = 'J m kg-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo
CASE (iarr_zlmix)
   sname = 'mixing_length_in_sea_water'
   fname = 'zlmix'
   lname = 'Mixing length'
   unit = 'm'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo

!
!3.7 Tides
!---------
!

CASE (iarr_astro_earth)
   sname = 'elasticity_factor_for_earth_tides'
   fname = 'astro_earth'
   lname = 'Elasticity factor for Earth tides'
   cnode = ''
   nodim = 1
   ngdims(1) = MaxAstroTides
CASE (iarr_fnode_astro)
   sname = 'nodal_factor_of_tidal_force'
   fname = 'fnode_astro'
   lname = 'Nodal factor of tidal force'
   cnode = ''
   nodim = 1
   ngdims(1) = nconastro
CASE (iarr_fnode_obc)
   sname = 'nodal_factors_of_tidal_constituents_at_open_boundaries'
   fname = 'fnode_obc'
   lname = 'Nodal factors of tidal constituents at open boundaries'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_fxastro)
   sname = 'x_tidal_force'
   fname = 'fxastro'
   lname = 'X-component of tidal force'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_fyastro)
   sname = 'y_tidal_force'
   fname = 'fyastro'
   lname = 'Y-component of tidal force'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_index_astro)
   sname = 'key_ids_of_tidal_force_constituents'
   fname = 'index_astro'
   lname = 'Key ids of tidal force constituents'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nconastro
CASE (iarr_index_obc)
   sname = 'key_ids_of_tidal_force_constituents_at_open_boundaries'
   fname = 'index_obc'
   lname = 'Key ids of tidal force constituents at open boundaries'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_ispec_tides)
   sname = 'tidal_species' 
   fname = 'ispec_tides'
   lname = 'Tidal species'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = MaxConstituents
CASE (iarr_phase_astro)
   sname = 'astronomical_phase_of_tidal_force'
   fname = 'phase_astro'
   lname = 'Astronomical phase of tidal force'
   unit = 'radian'
   cnode = ''
   nodim = 1
   ngdims(1) = nconastro
CASE (iarr_phase_obc)
   sname = 'astronomical_phases_at_open_boundaries'
   fname = 'phase_obc'
   lname = 'Astronomical phases at open boundaries'
   unit = 'radian'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_tidal_spectrum)
   sname = 'tidal_frequency' 
   fname = 'tidal_spectrum'
   lname = 'Tidal frequencies'
   unit = 'radian s-1'
   cnode = ''
   nodim = 1
   ngdims(1) = MaxConstituents

!
!3.8 Meteo forcing
!-----------------
!

CASE (iarr_airtemp)
   sname = 'air_temperature'
   commt = ''
   fname = 'airtemp'
   lname = 'Air temperature'
   unit = 'degC'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_atmpres)
   sname = 'air_pressure_at_sea_level'
   commt = ''
   fname = 'atmpres'
   lname = 'Atmospheric pressure'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_cloud_cover)
   sname = 'cloud_area_fraction'
   commt = ''
   fname = 'cloud_cover'
   lname = 'Cloud cover'
   unit = '1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_evapminprec)
   sname = 'evaporation_minus_precipitation_flux'
   fname = 'evapminprec'
   lname = 'Evaporation minus precipitation rate'
   unit = 'kg m-2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_evaporation)
   sname = 'evaporation_flux'
   commt = ''
   fname = 'evaporation'
   lname = 'Evaporation rate'
   unit = 'kg m-2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_precipitation)
   sname = 'precipitation_flux' 
   fname = 'precipitation'
   lname = 'Precipitation rate'
   unit = 'kg m-2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_relhum)
   sname = 'relative_humidity'
   commt = ''
   fname = 'relhum'
   lname = 'Relative humidity'
   unit = '1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_sst)
   sname = 'sea_surface_temperature'
   commt = ''
   fname = 'sst'
   lname = 'Sea surface temperature'
   unit = 'degC'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_uwindatc)
   sname = 'x_wind' 
   commt = ''
   fname = 'uwindatc'
   lname = 'X-component of surface wind'
   vname = 'Surface wind'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_uwindatc_old)
   sname = 'x_wind_at_old_time_level'
   fname = 'uwindatc_old'
   lname = 'X-component of surface wind at old time level'
   unit = 'm s-1'
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_vappres_air)
   sname = 'saturated_vapour_pressure' 
   commt = ''
   fname = 'vappres'
   lname = 'Saturated vapour essure'
   unit = 'mbar'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_vwindatc)
   sname = 'y_wind' 
   commt = ''
   fname = 'vwindatc'
   lname = 'Y-component of surface wind'
   vname = 'Surface wind'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_vwindatc_old)
   sname = 'Y_wind_at_old_time_level'
   fname = 'vwindatc_old'
   lname = 'Y-component of surface wind at old time level'
   unit = 'm s-1'
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_windatc)
   sname = 'wind_speed'
   commt = ''
   fname = 'windatc'
   lname = 'Wind speed'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)

!
!3.9 Waves
!---------
!

CASE (iarr_gxcoordglbwav)
   sname = 'global_x_coordinates_wave_grid'
   fname = 'gxcoordglbwav'
   lname = 'Global X-coordinates on the wave model grid'
   unit = MERGE('m          ','degree_east',iopt_grid_sph.EQ.0)
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/ncwav,nrwav/)
CASE (iarr_gycoordglbwav)
   sname = 'global_y_coordinates_wave_grid'
   fname = 'gycoordglbwav'
   lname = 'Global Y-coordinates on the wave model grid'
   unit = MERGE('m           ','degree_north',iopt_grid_sph.EQ.0)
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/ncwav,nrwav/)
CASE (iarr_hbwdissipmag)
   sname = 'magnitude_of_bottom_wave_dissipation'
   fname = 'hbwdissipmag'
   lname = 'Magnitude of bottom wave dissipation'
   unit = 'm s-2'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_hmbwdissipmag)
   sname = 'magnitude_of_depth_mean_bottom_wave_dissipation'
   fname = 'hmbwdissipmag'
   lname = 'Magnitude of depth-mean bottom wave dissipation'
   unit = 'm s-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_hmswdissipmag)
   sname = 'magnitude_of_depth_mean_surface_wave_dissipation'
   fname = 'hmswdissipmag'
   lname = 'Magnitude of depth-mean surface wave dissipation'
   unit = 'm s-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_hswdissipmag)
   sname = 'magnitude_of_surface_wave_dissipation'
   fname = 'hswdissipmag'
   lname = 'Magnitude of surface wave dissipation'
   unit = 'm s-2'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_maskglbwav)
   sname = 'global_land_mask_wave_grid'
   fname = 'maskglbwav'
   lname = 'Global land mask on the wave model grid'
   unit = 'log'
   dtype = log_type
   nodim = 2
   ngdims(1:2) = (/ncwav,nrwav/)
CASE (iarr_ubwdissipatc)
   sname = 'x_bed_wave_dissipation_force_at_c_nodes'
   fname = 'ubwdissipatc'
   lname = 'X-component of bed wave dissipation force at C-nodes'
   unit = 'm s-2'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,0,0,0/)
CASE (iarr_umbwdissipatc)
   sname = 'integral_wrt_depth_of_x_bed_wave_dissipation_force_at_c_nodes'
   fname = 'umbwdissipatc'
   lname = 'X-component of depth-integrated bed wave dissipation force at '//&
           'C-nodes'
   unit = 'm s-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,0,0/)
CASE (iarr_umswdissipatc)
   sname = 'integral_wrt_depth_of_x_surface_wave_dissipation_force_at_c_nodes'
   fname = 'umswdissipatc'
   lname = 'X-component of depth-integrated surface wave dissipation force '//&
         & 'at C-nodes'
   unit = 'm s-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,0,0/)
CASE (iarr_uswdissipatc)
   sname = 'x_surface_wave_dissipation_force_at_c_nodes'
   fname = 'uswdissipatc'
   lname = 'X-component of surface wave dissipation force at C-nodes'
   unit = 'm s-2'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,0,0,0/)
CASE (iarr_vbwdissipatc)
   sname = 'y_bed_wave_dissipation_force_at_c_nodes'
   fname = 'vbwdissipatc'
   lname = 'Y-component of bed wave dissipation force at C-nodes'
   unit = 'm s-2'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,0,1,0/)
CASE (iarr_vmbwdissipatc)
   sname = 'integral_wrt_depth_of_y_bed_wave_dissipation_force_at_c_nodes'
   fname = 'vmbwdissipatc'
   lname = 'Y-component of depth-integrated bed wave dissipation force at '//&
         & 'C-nodes'
   unit = 'm2 s-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,1,0/)
CASE (iarr_vmswdissipatc)
   sname = 'integral_wrt_depth_of_y_surface_wave_dissipation_force_at_c_nodes'
   fname = 'vmswdissipatc'
   lname = 'Y-component of depth-integrated bed wave dissipation force at '//&
         & 'C-nodes'
   unit = 'm2 s-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,1,0/)
CASE (iarr_vswdissipatc)
   sname = 'y_surface_wave_dissipation_force_at_c_nodes'
   fname = 'vswdissipatc'
   lname = 'Y-component of surface wave dissipation force at C-nodes'
   unit = 'm s-2'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,0,1,0/)
CASE (iarr_wavedir)
   sname = 'surface_wave_direction_wrt_reference_grid'
   fname = 'wavedir'
   lname = 'Wave direction'
   unit = 'degrees'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_waveexcurs)
   sname = 'near_bed_orbital_wave_excursion'
   fname = 'waveexcurs'
   lname = 'Near bed orbital wave excursion'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_wavefreq)
   sname = 'wave_frequency'
   fname = 'wavefreq'
   lname = 'Peak wave frequency'
   unit = 'radian s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_waveheight)
   sname = 'wave_height'
   fname = 'waveheight'
   lname = 'Significant wave height'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_wavenum)
   sname = 'wave_number'
   fname = 'wavenum'
   lname = 'Wave number'
   unit = 'm -1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_waveperiod)
   sname = 'wave_period'
   fname = 'waveperiod'
   lname = 'Peak wave period'
   unit = 's'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_wavepres)
   sname = 'wave_induced_pressure'
   fname = 'wavepres'
   lname = 'Wave induced pressure'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_wavevel)
   sname = 'near_bed_wave_orbital_velocity'
   fname = 'wavevel'
   lname = 'Near bed wave orbital velocity'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
  
!
!3.10 Stokes velocities
!----------------------
!

CASE (iarr_hmstokesmag)
   sname = 'depth_mean_magnitude_of_stokes_velocity'
   fname = 'hmstokesmag'
   lname = 'Magnitude of depth-mean Stokes current'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_hmveltotmag)
   sname = 'depth_mean_magnitude_of_total_horizontal_sea_water_velocity'
   fname = 'hmveltotmag'
   lname = 'Magnitude of total depth-mean current'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_hstokesmag)
   sname = 'magnitude_of_horizontal_stokes_velocity'
   fname = 'hstokesmag'
   lname = 'Magnitude of Stokes current'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_hveltotmag)
   sname = 'magnitude_of_total_horizontal_velocity'
   fname = 'hveltotmag'
   lname = 'Magnitude of total current'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_stokessource2du)
   sname = 'stokes_source_terms_2d_u_momentum_equation'
   fname ='stokessource2du'
   lname = 'Stokes source term in the U-equation'
   unit = 'm2 s2'
   cnode  = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_stokessource2dv)
   sname = 'stokes_source_terms_2d_v_momentum_equation'
   fname ='stokessource2dv'
   lname = 'Stokes source term in the V-equation'
   unit = 'm2 s2'
   cnode  = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_udstokesatu)
   sname = 'integral_wrt_depth_of_x_sea_water_stokes_velocity_at_u_nodes'
   fname = 'udstokesatu'
   lname = 'X-component of depth-integrated Stokes current at U-nodes'
   vname = 'Depth-integrated Stokes current'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,0/)
CASE (iarr_umstokesatc)
   sname = 'depth_mean_x_sea_water_stokes_velocity_at_c_nodes'
   fname = 'umstokesatc'
   lname = 'X-component of depth-mean Stokes current at C-nodes'
   vname = 'Depth-mean Stokes current'
   unit = 'm s-1'
   cnode = 'C'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
   nhdims = (/nhalo,nhalo,1,0/)
CASE (iarr_umstokesatu)
   sname = 'depth_mean_x_sea_water_stokes_velocity_at_u_nodes'
   fname = 'umstokesatu'
   lname = 'X-component of depth-mean Stokes current at U-nodes'
   vname = 'Depth-mean Stokes current'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/nhalo,nhalo,1,0/)
CASE (iarr_umveltot)
   sname = 'depth_mean_x_sea_water_total_velocity'
   fname = 'umveltot'
   lname = 'X-component of depth-mean total current'
   vname = 'Depth-mean Stokes current'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ustokesatc)
   sname = 'x_sea_water_stokes_velocity_at_c_nodes'
   fname = 'ustokesatc'
   lname = 'X-component of Stokes current at C-nodes'
   vname = 'Stokes current'
   unit = 'm s-1'
   cnode = 'C'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims =  (/nhalo,nhalo,1,0/)
CASE (iarr_ustokesatu)
   sname = 'x_sea_water_stokes_velocity_at_u_nodes'
   fname = 'ustokesatu'
   lname = 'X-component of Stokes current at U-nodes'
   vname = 'Stokes current'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims =  (/nhalo,nhalo,1,0/)
CASE (iarr_uveltot)
   sname = 'x_sea_water_total_velocity'
   fname = 'uveltot'
   lname = 'X-component of total current'
   vname = 'Total current'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_vdstokesatv)
   sname = 'integral_wrt_depth_of_y_sea_water_stokes_velocity_at_v_nodes'
   fname = 'vdstokesatv'
   lname = 'Y-component of depth-integrated Stokes current at V-nodes'
   vname = 'Depth-integrated Stokes current'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,0,1/)
CASE (iarr_vmstokesatc)
   sname = 'depth_mean_y_sea_water_stokes_velocity'
   fname = 'vmstokesatc'
   lname = 'Y-component of depth-mean Stokes current'
   vname = 'Depth-mean Stokes current'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
   nhdims = (/1,0,nhalo,nhalo/)
CASE (iarr_vmstokesatv)
   sname = 'depth_mean_y_sea_water_stokes_velocity_at_v_nodes'
   fname = 'vmstokesatv'
   lname = 'Y-component of depth-mean Stokes current at V-nodes'
   vname = 'Depth-mean Stokes current'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,nhalo,nhalo/)
CASE (iarr_vmveltot)
   sname = 'depth_mean_y_sea_water_total_velocity'
   fname = 'vmveltot'
   lname = 'Y-component of depth-mean total current'
   vname = 'Depth-mean Stokes current'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_vstokesatc)
   sname = 'y_sea_water_stokes_velocity_at_c_nodes'
   fname = 'vstokesatc'
   lname = 'Y-component of Stokes current at C-nodes'
   vname = 'Stokes current'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims =  (/1,0,nhalo,nhalo/)
CASE (iarr_vstokesatv)
   sname = 'y_sea_water_stokes_velocity_at_v_nodes'
   fname = 'vstokesatv'
   lname = 'Y-component of Stokes current at V-nodes'
   vname = 'Stokes current'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims =  (/1,0,nhalo,nhalo/)
CASE (iarr_vveltot)
   sname = 'y_sea_water_total_velocity'
   fname = 'vveltot'
   lname = 'Y-component of total current'
   vname = 'Total current'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_wstokesatw)
   sname = 'transformed_upward_sea_water_stokes_velocity'
   fname = 'wstokesatw'
   lname = 'Transformed vertical Stokes velocity'
   vname = 'Stokes current'
   unit = 'm s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = (/1,0,1,0/)

!
!3.11 Optical arrays
!-------------------
!

CASE (iarr_optattcoef2)
   sname = 'volume_attenuation_coefficient_of_short_wave_radiative_flux_'//&
         & 'in_sea_water'
   fname = 'optattcoef2'
   lname = 'Inverse optical attenuation depth for short waves'
   unit = 'm-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_qrad)
   sname = 'surface_downwelling_spherical_irradiance_in_sea_water'
   fname = 'qrad'
   lname = 'Surface solar irradiance'
   unit = 'W m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_radiance)
   sname = 'downwelling_sperical_irradiance_in_sea_water'
   fname = 'radiance'
   lname = 'Solar irradiance'
   unit = 'W m-2'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)

!
!3.12 Bottom/surface fluxes
!--------------------------
!

CASE (iarr_bdragcoefatc)
   sname = 'drag_coefficient_at_sea_floor_at_c_nodes'
   fname = 'bdragcoefatc'
   lname = 'Bottom drag coefficient'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_bfricatu)
   sname = 'friction_velocity_term_at_sea_floor_at_u_nodes'
   fname = 'bfricatu'
   lname = 'Bottom friction velocity term at U-nodes'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bfricatv)
   sname = 'friction_velocity_term_at_sea_floor_at_v_nodes'
   fname = 'bfricatv'
   lname = 'Bottom friction velocity term at V-nodes'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bfricvel)
   sname = 'friction_velocity_at_sea_floor'
   fname = 'bfricvel'
   lname = 'Bottom friction velocity'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bfricvel_max)
   sname = 'maximum_current_wave_friction_velocity_at_sea_floor'
   fname = 'bfricvel_max'
   lname = 'Maximum current-wave bottom friction velocity'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bfricvel_wav)
   sname = 'wave_friction_velocity_at_sea_floor'
   fname = 'bfricvel_wav'
   lname = 'Bottom wave friction'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bstresatc)
   sname = 'stress_at_sea_floor'
   fname = 'bstresatc'
   lname = 'Bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bstresatc_max)
   sname = 'maximum_current_wave_stress_at_sea_floor'
   fname = 'bstresatc_max'
   lname = 'Maximum current-wave bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bstresatc_wav)
   sname = 'maximum_wave_stress_at_sea_floor'
   fname = 'bstresatc_wav'
   lname = 'Maximum wave bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bstresatu)
   sname = 'stress_at_sea_floor'
   fname = 'bstresatu'
   lname = 'Bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bstresatv)
   sname = 'stress_at_sea_floor'
   fname = 'bstresatv'
   lname = 'Bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_cds)
   sname = 'drag_coefficient_at_sea_surface'
   fname = 'cds'
   lname = 'Surface drag coefficient'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_ces)
   sname = 'exchange_coefficient_for_latent_heat_flux' 
   fname = 'ces'
   lname = 'Exchange coefficient for latent heat flux'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_chs)
   sname = 'exchange_coefficient_for_sensible_heat_flux'
   fname = 'chs'
   lname = 'Exchange coefficient for sensible heat flux'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_fwave)
   sname = 'wave_friction_factor'
   fname = 'fwave'
   lname = 'Wave friction factor'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_qlatent)
   sname = 'surface_upward_latent_heat_flux'
   commt = ''
   fname = 'qlatent'
   lname = 'Latent heat flux'
   unit = 'W m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_qlwave)
   sname = 'surface_net_upward_longwave_flux'
   fname = 'qlwave'
   lname = 'Long-wave heat flux'
   unit = 'W m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_qnonsol)
   sname = 'surface_upward_non_solar_heat_flux_in_sea_water'
   fname = 'qnonsol'
   lname = 'Non-solar heat flux'
   unit = 'W m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_qsensible)
   sname = 'surface_upward_sensible_heat_flux'
   fname = 'qsensible'
   lname = 'Sensible heat flux'
   unit = 'W m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_qtot)
   sname = 'surface_downward_heat_flux_in_sea_water'
   fname = 'qtot'
   lname = 'Total downward surface heat flux'
   unit = 'W m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_sfricatc)
   sname = 'surface_firction_velocity'
   fname = 'sfricatc'
   lname = 'Surface friction velocity'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ssalflux)
!  ? + units : PSU m s-1 instead of kg m-2 s -1!!!!! 
   sname = 'virtual_salt_flux_into_sea_water' 
   fname = 'ssalflux'
   lname = 'Salinity flux'
   unit = 'PSU m s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_sstresatc)
   sname = 'surface_downward_stress'
   fname = 'sstresatc'
   lname = 'Surface stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ubstresatc)
   sname = 'upward_x_stress_at_sea_floor'
   fname = 'ubstresatc'
   lname = 'X-component of bottom stress at C-nodes'
   vname = 'Bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ubstresatu)
   sname = 'upward_x_stress_at_sea_floor'
   fname = 'ubstresatu'
   lname = 'X-component of bottom stress at U-nodes'
   vname = 'Bottom stress'
   unit = 'Pa'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,0/)
CASE (iarr_usstresatc)
   sname = 'surface_downward_x_stress'
   fname = 'usstresatc'
   lname = 'X-component of surface stress'
   vname = 'Surface stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,0,0/)
CASE (iarr_usstresatu)
   sname = 'surface_downward_x_stress_at_u_nodes'
   fname = 'usstresatu'
   lname = 'X-component of surface stress at U-nodes'
   vname = 'Surface stress'
   unit = 'Pa'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_vbstresatc)
   sname = 'upward_y_stress_at_sea_floor'
   fname = 'vbstresatc'
   lname = 'Y-component of bottom stress at C-nodes'
   vname = 'Bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_vbstresatv)
   sname = 'upward_y_stress_at_sea_floor'
   fname = 'vbstresatv'
   lname = 'Y-component of bottom stress'
   vname = 'Bottom stres'
   unit = 'Pa'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims =(/0,0,0,1/)
CASE (iarr_vsstresatc)
   sname = 'surface_downward_x_stress'
   fname = 'vsstresatc'
   lname = 'X-component of surface stress'
   vname = 'Surface stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,0,0/)
CASE (iarr_vsstresatv)
   sname = 'surface_downward_y_stress_at_v_nodes'
   fname = 'vsstresatv'
   lname = 'Y-component of surface stress at V-nodes'
   vname = 'Surface stress'
   unit = 'Pa'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_wavethickatc)
   sname = 'wave_boundary_layer_thickness'
   fname = 'wavethickatc'
   lname = 'Wave boundary layer thickness '
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_zaroughatc)
   sname = 'sea_floor_apparent_roughness_length' 
   fname = 'zaroughatc'
   lname = 'Bottom apparent roughness length'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_zroughatc)
   sname = 'sea_floor_roughness_length' 
   fname = 'zroughatc'
   lname = 'Bottom roughness length'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)

!
!3.14 Open boundary forcing
!--------------------------
!

CASE (iarr_floutobu)
   sname = 'last_output_time_at_u_open_boundaries'
   fname = 'floutobu'
   lname = 'Last output time at U-open boundaries'
   unit = 's'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,nz/)
CASE (iarr_floutobv)
   sname = 'last_output_time_at_v_open_boundaries'
   fname = 'floutobv'
   lname = 'Last output time at V-open boundaries'
   unit = 's'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,nz/)
CASE (iarr_gxslope)
   sname = 'x_derivative_of_kinematic_ocean_pressure'
   fname = 'gxslope'
   lname = 'X-component of kinematic pressure gradient'
   unit = 'm s-2'
   cnode = ''
   nodim = 0
CASE (iarr_gxslope_amp)
   sname = 'amplitude_of_x_derivative_of_kinematic_ocean_pressure'
   fname = 'gxslope_amp'
   lname = 'Amplitude of X-component of kinematic pressure gradient'
   unit = 'm s-2'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_gxslope_pha)
   sname = 'phase_of_x_derivative_of_kinematic_ocean_pressure'
   fname = 'gxslope_pha'
   lname = 'Phase of X-component of kinematic pressure gradient'
   unit = 'radian'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_gyslope)
   sname = 'y_derivative_of_kinematic_ocean_pressure'
   fname = 'gyslope'
   lname = 'Y-component of kinematic pressure gradient'
   unit = 'm s-2'
   cnode = ''
   nodim = 0
CASE (iarr_gyslope_amp)
   sname = 'Amplitude_of_y_derivative_of_kinematic_ocean_pressure'
   fname = 'gyslope_amp'
   lname = 'Amplitude of Y-component of kinematic pressure gradient'
   unit = 'm s-2'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_gyslope_pha)
   sname = 'phase_of_y_derivative_of_kinematic_ocean_pressure'
   fname = 'gyslope_pha'
   lname = 'Phase of Y-component of kinematic pressure gradient'
   unit = 'radian'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_iloczobu)
   sname = 'position_of_elevation_points_with_respect_to_u_open_boundaries'
   fname = 'iloczobu'
   lname = 'Position of elevation points with respect to U-open boundaries'
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = nobu
CASE (iarr_iloczobv)
   sname = 'position_of_elevation_points_with_respect_to_v_open_boundaries'
   fname = 'iloczobv'
   lname = 'Position of elevation points with respect to V-open boundaries'
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = nobv
CASE (iarr_indexprof)
   sname = 'mapping_array_of_profile_numbers'
   fname = 'indexprof'
   lname = 'Map between of profile numbers'
   dtype = int_type
   cnode = ''
   nodim = 2
CASE (iarr_indexvar)
   sname = 'mapping_array_of_variable_indices'
   fname = 'indexvar'
   lname = 'Map of variable indices'
   dtype = int_type
   cnode = ''
   nodim = 2
CASE (iarr_index2dobuv)
   sname = 'mapping_array_of_data_points_to_normal_open_boundary_locations'
   fname = 'index2dobuv'
   lname = 'Map of data points to normal open boundary locations'  
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nobu+nobv,maxdatafiles(io_2uvobc,1)-1/)
CASE (iarr_index2dobxy)
   sname = 'mapping_array_of_data_points_to_tangential_open_boundary_locations'
   fname = 'index2dobxy'
   lname = 'Map of data points to tangential open boundary locations'  
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nobx+noby,maxdatafiles(io_2xyobc,1)-1/)
CASE (iarr_iobc2dtype)
   sname = 'type_2d_open_boundary_data_input'
   fname = 'iobc2dtype'
   lname = 'Type 2-D open boundary data input'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = maxdatafiles(io_2uvobc,1)-1
CASE (iarr_iprofobu)
   sname = 'profile_number_at_u_open_boundaries'
   fname = 'iprofobu'
   lname = 'Profile number at U-open boundaries'
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = nobu
CASE (iarr_iprofobv)
   sname = 'profile_number_at_v_open_boundaries'
   fname = 'iprofobv'
   lname = 'Profile number at V-open boundaries'
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = nobv
CASE (iarr_iprofobx)
   sname = 'profile_number_at_x_open_boundaries'
   fname = 'iprofobu'
   lname = 'Profile number at X-open boundaries'
   dtype = int_type
   cnode = 'X'
   nodim = 1
   ngdims(1) = nobx
CASE (iarr_iprofoby)
   sname = 'profile_number_at_y_open_boundaries'
   fname = 'iprofoby'
   lname = 'Profile number at Y-open boundaries'
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   ngdims(1) = noby
CASE (iarr_iqsecobu)
   sname = 'start_end_locations_at_u_node_discharge_boundaries'
   fname = 'iqsecobu'
   lname = 'Start/end locations at U-node discharge boundaries'
   dtype = int_type
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nqsecobu,2/)
CASE (iarr_isur1dtype)
   sname = 'number_of_1d_surface_data'
   fname = 'isur1dtype'
   lname = 'Number of 1-D surface data'
   dtype = int_type
   cnode = ''
   nodim = 0
CASE (iarr_itypobu)
   sname = 'type_of_open_boundary_condition_at_u_open_boundaries'
   fname = 'itypobu'
   lname = 'Type of open boundary condition at U-open boundaries'
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = nobu
CASE (iarr_itypobv)
   sname = 'type_of_open_boundary_condition_at_v_open_boundaries'
   fname = 'itypobv'
   lname = 'Type of open boundary condition at V-open boundaries'
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = nobv
CASE (iarr_itypobx)
   sname = 'type_of_open_boundary_condition_at_x_open_boundaries'
   fname = 'itypobx'
   lname = 'Type of open boundary condition at X-open boundaries'
   dtype = int_type
   cnode = 'X'
   nodim = 1
   ngdims(1) = nobx
CASE (iarr_itypoby)
   sname = 'type_of_open_boundary_condition_at_y_open_boundaries'
   fname = 'itypoby'
   lname = 'Type of open boundary condition at Y-open boundaries'
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   ngdims(1) = noby
CASE (iarr_ityp2dobu)
   sname = 'type_of_2d_open_boundary_condition_at_u_nodes'
   fname = 'ityp2dobu'
   lname = 'Type of 2-D open boundary condition at U-nodes'
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = nobu
CASE (iarr_ityp2dobv)
   sname = 'type_of_2d_open_boundary_condition_at_v_nodes'
   fname = 'ityp2dobv'
   lname = 'Type of 2-D open boundary condition at V-nodes'
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = nobv
CASE (iarr_ityp2dobx)
   sname = 'type_of_2d_open_boundary_condition_at_x_nodes'
   fname = 'ityp2dobx'
   lname = 'Type of 2-D open boundary condition at X-nodes'
   dtype = int_type
   cnode = 'X'
   nodim = 1
   ngdims(1) = nobx
CASE (iarr_ityp2doby)
   sname = 'type_of_2d_open_boundary_condition_at_y_nodes'
   fname = 'ityp2doby'
   lname = 'Type of 2-D open boundary condition at Y-nodes'
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   ngdims(1) = nobv
CASE (iarr_jqsecobv)
   sname = 'start_end_locations_at_v_node_discharge_boundaries'
   fname = 'jqsecobv'
   lname = 'Start/end locations at V-node discharge boundaries'
   dtype = int_type
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nqsecobv,2/)
CASE (iarr_noprofsd)
   sname = 'number_of_profiles_per_input_file'
   fname = 'noprofsd'
   lname = 'Number of profiles per input file'
   dtype = int_type
   cnode = ''
   nodim = 1
CASE (iarr_no2dobuv)
   sname = 'number_of_input_data_per_input_file_at_normal_open_boundaries'
   fname = 'no2dobuv'
   lname = 'Number of input data per input file'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = maxdatafiles(io_2uvobc,1)-1
CASE (iarr_no2dobxy)
   sname = 'number_of_input_data_per_input_file_at_tangential_open_boundaries'
   fname = 'no2dobxy'
   lname = 'Number of input data per input file'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = maxdatafiles(io_2xyobc,1)-1
CASE (iarr_obcbioatu)
   sname = 'storage_array_for_biology_at_u_open_boundaries'
   fname = 'obcsalatu'
   lname = 'Storage array for biology at U-open boundaries'
   unit = 'PSU'
   cnode = 'U'
   nodim = 4
CASE (iarr_obcbioatv)
   sname = 'storage_array_for_biology_at_v_open_boundaries'
   fname = 'obcsalatv'
   lname = 'Storage array for biology at V-open boundaries'
   unit = 'PSU'
   cnode = 'V'
   nodim = 4
CASE (iarr_obcsalatu)
   sname = 'storage_array_for_salinity_at_u_open_boundaries'
   fname = 'obcsalatu'
   lname = 'Storage array for salinity at U-open boundaries'
   unit = 'PSU'
   cnode = 'U'
   nodim = 4
   ngdims(1:4) = (/nobu,nz,3,1/)
CASE (iarr_obcsalatv)
   sname = 'storage_array_for_salinity_at_v_open_boundaries'
   fname = 'obcsalatv'
   lname = 'Storage array for salinity at V-open boundaries'
   unit = 'PSU'
   cnode = 'V'
   nodim = 4
   ngdims(1:4) = (/nobv,nz,3,1/)
CASE (iarr_obctmpatu)
   sname = 'storage_array_for_temperature_at_u_open_boundaries'
   fname = 'obctmpatu'
   lname = 'Storage array for temperature at U-open boundaries'
   unit = 'degC'
   cnode = 'U'
   nodim = 4
   ngdims(1:4) = (/nobu,nz,3,1/)
CASE (iarr_obctmpatv)
   sname = 'storage_array_for_temperature_at_v_open_boundaries'
   fname = 'obctmpatv'
   lname = 'Storage array for temperature at V-open boundaries'
   unit = 'degC'
   cnode = 'V'
   nodim = 4
   ngdims(1:4) = (/nobv,nz,3,1/)
CASE (iarr_obc2uvatu)
   sname = 'storage_array_for_2d_mode_at_u_open_boundaries'
   fname = 'obc2uvatu'
   lname = 'Storage array for 2-D mode at U-open boundaries'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,2/)
CASE (iarr_obc2uvatu_old)
   sname = 'storage_array_for_2d_mode_at_u_open_boundaries_at_old_time_level'
   fname = 'obc2uvatu_old'
   lname = 'Storage array for 2-D mode at U-open boundaries and old time level'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,2/)
CASE (iarr_obc2uvatv)
   sname = 'storage_array_for_2d_mode_at_v_open_boundaries'
   fname = 'obc2uvatv'
   lname = 'Storage array for 2-D mode at V-open boundaries'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,2/)
CASE (iarr_obc2uvatv_old)
   sname = 'storage_array_for_2d_mode_at_v_open_boundaries_at_old_time_level'
   fname = 'obc2uvatv_old'
   lname = 'Storage array for 2-D mode at V-open boundaries and old time level'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,2/)
CASE (iarr_obc3uvatu)
   sname = 'storage_array_for_3d_mode_at_u_open_boundaries'
   fname = 'obc3uvatu'
   lname = 'Storage array for 3-D mode at U-open boundaries'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 3
   ngdims(1:3) = (/nobu,nz,2/)
CASE (iarr_obc3uvatv)
   sname = 'storage_array_for_3d_mode_at_v_open_boundaries'
   fname = 'obc3uvatv'
   lname = 'Storage array for 3-D mode at V-open boundaries'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 3
   ngdims(1:3) = (/nobv,nz,2/)
CASE (iarr_profbio)
   sname = 'biological_concentrations_at_open_boundaries'
   fname = 'profbio'
   lname = 'Biology profiles at open boundaries'
   cnode = ''
   nodim = 3
   ngdims(3) = nz
CASE (iarr_profsal)
   sname = 'sea_water_salinity_at_open_boundaries'
   fname = 'profsal'
   lname = 'Salinity profiles at open boundaries'
   unit = 'PSU'
   cnode = ''
   nodim = 2
   ngdims(2) = nz
CASE (iarr_profsed)
   sname= 'sediment_concentrations_at_open_boundaries'
   fname = 'profsed'
   lname = 'Sediment concentration at open boundaries'
   cnode = ''
   nodim = 3
   ngdims(3) = nz
CASE (iarr_proftmp)
   sname = 'sea_water_temperature_at_open_boundaries'
   fname = 'proftmp'
   lname = 'Temperature profiles at open boundaries'
   unit = 'degC'
   cnode = ''
   nodim = 2
   ngdims(2) = nz
CASE (iarr_profvel)
   sname = 'baroclinic_sea_water_velocity_at_open_boundaries'
   fname = 'profvel'
   lname = 'Baroclinic current profiles at open boundaries'
   unit = 'degC'
   cnode = ''
   nodim = 2
   ngdims(2) = nz
CASE (iarr_return_time)
   sname = 'return_time_at_open_boundaries'
   fname = 'return_time'
   lname = 'Vertical profile of return times for TH-open boundary condition'
   unit = 's'
   cnode = ''
   nodim = 1
   ngdims(1) = nz
CASE (iarr_udatobu)
   sname = 'integral_wrt_depth_of_x_sea_water_velocity_at_u_open_boundaries'
   fname = 'udatobu'
   lname = 'X-component of depth-integrated current or discharge at '//&
         & 'U-normal boundaries'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,2/)
CASE (iarr_udatobu_amp)
   sname = 'amplitude_of_integral_wrt_depth_of_x_sea_water_velocity_at_'//&
         & 'u_open_boundaries'
   fname = 'udatobu_amp'
   lname = 'Amplitude of X-component of depth-integrated current or '//&
         & 'discharge at U-open boundaries'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,nconobc/)
CASE (iarr_udatobu_pha)
   sname = 'phase_of_integral_wrt_depth_of_x_sea_water_velocity_at_'//&
         & 'u_open_boundaries'
   fname = 'udatobu_pha'
   lname = 'Phase of X-component of depth-integrated current or discharge '//&
         & 'at U-open boundaries'
   unit = 'radian'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,nconobc/)
CASE (iarr_udatoby)
   sname = 'integral_wrt_depth_of_x_sea_water_velocity_wrt_depth_at_'//&
         & 'y_open_boundaries'
   fname = 'udatoby'
   lname = 'X-component of depth-integrated current at U-normal boundaries'
   unit = 'm2 s-1'
   cnode = 'Y'
   nodim = 2
   ngdims(1:2) = (/noby,2/)
CASE (iarr_udatoby_amp)
   sname = 'amplitude_of_integral_wrt_depth_of_x_sea_water_velocity_at_'//&
         & 'y_open_boundaries'
   fname = 'udatoby_amp'
   lname = 'Amplitude of X-component of depth-integrated current at Y-open '//&
         & 'boundaries'
   unit = 'm2 s-1'
   cnode = 'Y'
   nodim = 2
   ngdims(1:2) = (/noby,nconobc/)
CASE (iarr_udatoby_pha)
   sname = 'phase_of_integral_wrt_depth_of_x_sea_water_velocity_at_'//&
         & 'y_open_boundaries'
   fname = 'udatoby_pha'
   lname = 'Phase of X-component of depth-integrated current at Y-open '//&
         & 'boundaries'
   unit = 'radian'
   cnode = 'Y'
   nodim = 2
   ngdims(1:2) = (/noby,nconobc/)
CASE (iarr_vdatobv)
   sname = 'integral_wrt_depth_of_y_sea_water_velocity_at_v_open_boundaries' 
   fname = 'vdatobv'
   lname = 'Y-component of depth-integrated current or discharge at V-open '//&
         & 'boundaries'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,2/)
CASE (iarr_vdatobv_amp)
   sname = 'amplitude_of_integral_wrt_depth_of_y_sea_water_velocity_at_v_'//&
         & 'open_boundaries'
   fname = 'vdatobv_amp'
   lname = 'Amplitude of Y-component of depth-integrated current or '//&
         & 'discharge at V-open boundaries'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,nconobc/)
CASE (iarr_vdatobv_pha)
   sname = 'phase_of_integral_wrt_depth_of_y_sea_water_velocity_at_v_'//&
         & 'open_boundaries'
   fname = 'vdatobv_pha'
   lname = 'Phase of Y-component of depth-integrated current or discharge '//&
         & 'at V-open boundaries'
   unit = 'radian'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,nconobc/)
CASE (iarr_vdatobx)
   sname = 'integral_wrt_depth_of_y_sea_water_velocity_at_x_open_boundaries'
   fname = 'vdatobx'
   lname = 'Y-component of depth-integrated current at X-normal boundaries'
   unit = 'm2 s-1'
   cnode = 'X'
   nodim = 2
   ngdims(1:2) = (/nobx,2/)
CASE (iarr_vdatobx_amp)
   sname = 'amplitude_of_integral_wrt_depth_of_y_sea_water_velocity_at_x_'//&
         & 'open_boundaries'
   fname = 'vdatobx_amp'
   lname = 'Amplitude of Y-component of depth-integrated current at X-open '//&
         & 'boundaries'
   unit = 'm2 s-1'
   cnode = 'X'
   nodim = 2
   ngdims(1:2) = (/nobx,nconobc/)
CASE (iarr_vdatobx_pha)
   sname = 'phase_of_integral_wrt_depth_of_y_sea_water_velocity_at_x_'//&
         & 'open_boundaries'
   fname = 'vdatobx_pha'
   lname = 'Phase of Y-component of depth-integrated current at X-open '//&
         & 'boundaries'
   unit = 'radian'
   cnode = 'X'
   nodim = 2
   ngdims(1:2) = (/nobx,nconobc/)
CASE (iarr_vel2dobc)
   sname = 'transports_at_open_boundaries'
   fname = 'vel2dobc'
   lname = 'Transports at open boundaries'
   unit = 'm2 s-1'
   cnode = ''
   nodim = 1
CASE (iarr_zdatobu)
   sname = 'sea_surface_height_above_sea_level_at_u_open_boundaries'
   fname = 'zdatobu'
   lname = 'Surface elevation at U-open boundaries'
   unit = 'm'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,2/)
CASE (iarr_zdatobu_amp)
   sname = 'sea_surface_height_amplitude_at_u_open_boundaries'
   fname = 'zdatobu_amp'
   lname = 'Amplitude of surface elevation at U-open boundaries'
   unit = 'm'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) =(/nobu,nconobc/)
CASE (iarr_zdatobu_pha)
   sname = 'sea_surface_height_phase_at_u_open_boundaries'
   fname = 'zdatobu_pha'
   lname = 'Phase of surface elevation at U-open boundaries'
   unit = 'radian'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,nconobc/)
CASE (iarr_zdatobv)
   sname = 'sea_surface_height_at_v_open_boundaries' 
   fname = 'zdatobv'
   lname = 'Surface elevation at V-open boundaries'
   unit = 'm'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,2/)
CASE (iarr_zdatobv_amp)
   sname = 'sea_surface_height_amplitude_at_v_open_boundaries'
   fname = 'zdatobv_amp'
   lname = 'Amplitude of surface elevation at V-open boundaries'
   unit = 'm'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,nconobc/)
CASE (iarr_zdatobv_pha)
   sname = 'sea_surface_height_phase_at_v_open_boundaries'
   fname = 'zdatobv_pha'
   lname = 'Phase of surface elevation at V-open boundaries'
   unit = 'radian'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,nconobc/)
CASE (iarr_zetaobc)
   sname = 'sea_surface_height_above_sea_level_at_open_boundaries'
   fname = 'zetaobc'
   lname = 'Elevations at open boundaries'
   unit = 'm'
   cnode = ''
   nodim = 1
CASE (iarr_zetasur)
   fname = 'zetasur'
   lname = 'Elevations for surface forcing'
   unit = 'm'
   cnode = ''
   nodim = 0
CASE (iarr_zeta_amp)
   sname = 'sea_surface_height_amplitude'
   fname = 'zeta_amp'
   lname = 'Amplitude of surface elevation'
   unit = 'm'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_zeta_pha)
   sname = 'sea_surface_height_phase'
   fname = 'zeta_pha'
   lname = 'Phase of surface elevation'
   unit = 'radian'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc

!
!3.15 Structure module
!---------------------
!

CASE (iarr_idry)
   sname = 'global_x_index_at_dry_cells'
   fname = 'idry'
   lname = 'Global X-index of dry cells' 
   dtype = int_type
   nodim = 1
   ngdims(1) = numdry
CASE (iarr_indexwbaru)
   sname = 'index_at_u_node_weirs_barriers_for_local_to_global_mapping'
   fname ='indexwbaru'
   lname = 'Index at U-node weirs and barriers for local to global mapping'
   dtype = int_type
   cnode ='U'
   nodim = 1
   nldims(1) = numwbaruloc
CASE (iarr_indexwbaruprocs)
   sname = 'index_at_u_node_weirs_barriers_for_local_to_global_mapping_'//&
         & 'per_process'
   fname = 'indexwbaruprocs'
   lname = 'Index at U-node weirs and barriers for local to global '//&
         & 'mapping per process'
   dtype = int_type
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/numwbaru,nprocs/)
CASE (iarr_indexwbarv)
   sname = 'index_at_v_node_weirs__barriers_for_local_to_global_mapping'
   fname ='indexwbarv'
   lname = 'Index at V-node weirs and barriers for local to global mapping'
   dtype = int_type
   cnode ='V'
   nodim = 1
   nldims(1) = numwbarvloc
CASE (iarr_indexwbarvprocs)
   sname = 'index_at_v_node_weirs_barriers_for_local_to_global_mapping_'//&
         & 'per_process'
   fname = 'indexwbarvprocs'
   lname = 'Index at V-node weirs and barriers for local to global mapping '//&
         & 'per process'
   dtype = int_type
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/numwbarv,nprocs/)
CASE (iarr_ithinu)
   sname = 'global_x_index_at_thin_dams_at_u_nodes'
   fname = 'ithinu'
   lname = 'Global X-index of U-thin dams' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = numthinu
CASE (iarr_ithinuloc)
   sname = 'local_x_index_at_thin_dams_at_u_nodes'
   fname = 'ithinuloc'
   lname = 'Local X-index of U-thin dams' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   nldims(1) = numthinuloc
CASE (iarr_ithinv)
   sname = 'global_x_index_at_thin_dams_at_v_nodes'
   fname = 'ithinv'
   lname = 'Global X-index of V-thin dams' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = numthinv
CASE (iarr_ithinvloc)
   sname = 'local_x_index_at_thin_dams_at_v_nodes'
   fname = 'ithinvloc'
   lname = 'Local X-index of V-thin dams' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   nldims(1) = numthinvloc
CASE (iarr_iwbaru)
   sname = 'global_x_index_at_weirs_barriers_at_u_nodes'
   fname = 'ibwbaru'
   lname = 'Global X-index of U-node weirs/barriers' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_iwbaruloc)
   sname = 'local_x_index_at_weirs_barriers_at_u_nodes'
   fname = 'ibwbaruloc'
   lname = 'Local X-index of U-node weirs/barriers' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   nldims(1) = numwbaruloc
CASE (iarr_iwbarv)
   sname = 'global_x_index_at_weirs_barriers_at_v_nodes'
   fname = 'ibwbarv'
   lname = 'Global X-index of V-node weirs/barriers' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv
CASE (iarr_iwbarvloc)
   sname = 'local_x_index_at_weirs_barriers_at_v_nodes'
   fname = 'ibwbarvloc'
   lname = 'Local X-index of V-node weirs/barriers' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   nldims(1) = numwbarvloc
CASE (iarr_jdry)
   sname = 'global_y_index_at_dry_cells'
   fname = 'jdry'
   lname = 'Global Y-index of dry cells' 
   dtype = int_type
   nodim = 1
   ngdims(1) = numdry
CASE (iarr_jthinu)
   sname = 'global_y_index_at_thin_dams_at_u_nodes'
   fname = 'jthinu'
   lname = 'Global Y-index of U-thin dams' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = numthinu
CASE (iarr_jthinuloc)
   sname = 'local_y_index_at_thin_dams_at_u_nodes'
   fname = 'jthinuloc'
   lname = 'Local Y-index of U-thin dams' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   nldims(1) = numthinuloc
CASE (iarr_jthinv)
   sname = 'global_y_index_at_thin_dams_at_v_nodes'
   fname = 'jthinv'
   lname = 'Global Y-index of V-thin dams' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = numthinv
CASE (iarr_jthinvloc)
   sname = 'local_y_index_at_thin_dams_at_v_nodes'
   fname = 'jthinvloc'
   lname = 'Local Y-index of V-thin dams' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   nldims(1) = numthinvloc
CASE (iarr_jwbaru)
   sname = 'global_y_index_at_weirs_barriers_at_u_nodes'
   fname = 'jwbaru'
   lname = 'Global Y-index of U-weirs/barriers' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_jwbaruloc)
   sname = 'local_y_index_at_weirs_barriers_at_u_nodes'
   fname = 'jwbaruloc'
   lname = 'Global Y-index of U-node weirs/barriers' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   nldims(1) = numwbaruloc
CASE (iarr_jwbarv)
   sname = 'global_y_index_at_weirs_barriers_at_v_nodes'
   fname = 'jwbarv'
   lname = 'Global Y-index of V-node weirs/barriers' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv
CASE (iarr_jwbarvloc)
   sname = 'local_y_index_at_weirs_barriers_at_v_nodes'
   fname = 'jwbarvloc'
   lname = 'Local Y-index of V-node weirs/barriers' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   nldims(1) = numwbarvloc
CASE (iarr_nowbaruprocs)
   sname = 'number_of_local_u_weirs_barriers_per_process' 
   fname = 'nowbaruprocs'
   lname = 'Number of local U-node weirs/barriers'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nowbarvprocs)
   sname = 'number_of_local_v_weirs_barriers_per_process' 
   fname = 'nowbarvprocs'
   lname = 'Number of local V-node weirs/barriers'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_oricoefu)
   sname = 'discharge_coefficient_at_u_node_orifices'
   fname = 'oricoefu'
   lname = 'Discharge coefficient for orifices at U-nodes'
   unit = 'm1/2 s-1'
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_oricoefv)
   sname = 'discharge_coefficient_at_v_node_orifices'
   fname = 'oricoefv'
   lname = 'Discharge coefficient for orifices at V-nodes'
   unit = 'm1/2 s-1'
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv
CASE (iarr_oriheightu)
   sname = 'width_of_u_node_orifices'
   fname = 'oriheightu'
   lname = 'Orifice width at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_oriheightv)
   sname = 'width_of_u_node_orifices'
   fname = 'oriheightv'
   lname = 'Orifice width at V-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv
CASE (iarr_orisillu)
   sname = 'heigth_of_u_node_orifices_wtr_sea_bed'
   fname = 'orisillu'
   lname = 'Orifice height at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_orisillv)
   sname = 'heigth_of_v_node_orifices_wtr_sea_bed'
   fname = 'orisillv'
   lname = 'Orifice height at V-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv
CASE (iarr_wbarcoefu)
   sname = 'discharge_coefficient_at_u_node_weirs'
   fname = 'wbarcoefu'
   lname = 'Discharge coefficient for weirs/barriers at U-nodes'
   unit = 'm1/2 s-1'
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_wbarcoefv)
   sname = 'discharge_coefficient_at_v_node_weirs'
   fname = 'wbarcoefv'
   lname = 'Discharge coefficient for weirs/barriers at V-nodes'
   unit = 'm1/2 s-1'
   cnode = 'V'
   nodim = 1
   nldims(1) = numwbarv
CASE (iarr_wbarcrestu)
   sname = 'crest_heigth_wrt_sea_bed_at_u_node_weirs'
   fname = 'wbarcrestu'
   lname = 'Height of weir crest at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_wbarcrestv)
   sname = 'crest_heigth_wrt_sea_bed_at_v_node_weirs'
   fname = 'wbarcrestv'
   lname = 'Height of weir crest at V-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv
CASE (iarr_wbarmodlu)
   sname = 'modular_limit_at_u_node_weirs'
   fname = 'wbarmodlu'
   lname = 'Modular limit at U-node weirs'
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_wbarmodlv)
   sname = 'modular_limit_at_v_node_weirs'
   fname = 'wbarmodlv'
   lname = 'Modular limit at V-node weirs'
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv
CASE (iarr_wbarelossu)
   sname= 'energy_loss_at_u_node_weirs_barriers'
   fname = 'wbarelossu'
   lname = 'Energy loss sink term at U-node weirs/barriers'
   unit = 's-1'
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_wbarelossv)
   sname= 'energy_loss_at_v_node_weirs_barriers'
   fname = 'wbarelossv'
   lname = 'Energy loss sink term at V-node weirs/barriers'
   unit = 's-1'
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv

!
!3.16 Discharges
!---------------
!

CASE (iarr_disarea)
   sname = 'discharge_area'
   fname = 'disarea'
   lname = 'Discharge area'
   unit = 'm2'
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_disdir)
   sname = 'discharge_direction'
   fname = 'disdir'
   lname = 'Discharge direction'
   unit = 'radian'
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_disflag)
   sname = 'discharge_flag'
   fname = 'disflag'
   lname = 'Discharge flag'
   unit = 'log'
   dtype = log_type
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_dissal)
   sname = 'salinity_discharge'
   fname = 'disspeed'
   lname = 'Salinity discharge'
   unit = 'PSU'
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_disspeed)
   sname = 'discharge_speed'
   fname = 'disspeed'
   lname = 'Discharge speed'
   unit = 'm s-1'
   nodim = 1
   nldims(1) = numdisloc
CASE (iarr_distmp)
   sname = 'temperature_discharge'
   fname = 'distmp'
   lname = 'Temperature discharge'
   unit = 'degC'
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_disvol)
   sname = 'volume_discharge_rate'
   fname = 'disvol'
   lname = 'Volume discharge rate'
   unit = 'm3 s-1'
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_idis)
   sname = 'global_x_index_at_discharge_location'
   fname = 'idis'
   lname = 'Global X-index of discharge location'
   dtype = int_type
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_idisloc)
   sname = 'local_x_index_at_discharge_locations'
   fname = 'idisloc'
   lname = 'Local X-index of discharge locations'
   dtype = int_type
   nodim = 1
   nldims(1) = numdisloc_ext
CASE (iarr_indexdisloc)
   sname = 'index_at_discharge_locations_for_local_to_global_mapping'
   fname = 'indexdisloc'
   lname = 'Index at discharge locations for local to global mapping'
   dtype = int_type
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_indexdisprocs)
   sname = 'index_at_discharge_locations_for_local_to_global_mapping_'//&
         & 'per_process'
   fname = 'indexdisprocs'
   lname = 'Index at discharge locations for local to global mapping '//&
         & 'per process'
   dtype = int_type
   nodim = 2
   ngdims(1:2) = (/numdis,nprocs/)
CASE (iarr_jdis)
   sname = 'global_y_index_at_discharge_locations'
   fname = 'jdis'
   lname = 'Global Y-index of discharge locations'
   dtype = int_type
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_jdisloc)
   sname = 'local_y_index_at_discharge_locations'
   fname = 'jdisloc'
   lname = 'Local Y-index of discharge locations'
   dtype = int_type
   nodim = 1
   nldims(1) = numdisloc_ext
CASE (iarr_kdis)
   sname = 'global_z_index_at_discharge_locations'
   fname = 'kdis'
   lname = 'Vertical grid index of discharge locations'
   dtype = int_type
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_kdistype)
   sname = 'type_of_vertical_discharge_location'
   fname = 'kdistype'
   lname = 'Type of vertical discharge location'
   dtype = int_type
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_nodisprocs)
   sname = 'number_of_local_discharge_locations_per_process' 
   fname = 'nodisprocs'
   lname = 'Number of local discharge locations perv process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_xdiscoord)
   sname = 'x_coordinate_of_discharge_locations'
   fname = 'xdiscoord'
   lname = 'X-coordinates of discharge locations'
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_ydiscoord)
   sname = 'y_coordinate_of_discharge_locations'
   fname = 'ydiscoord'
   lname = 'Y-coordinates of discharge locations'
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_zdiscoord)
   sname = 'z_coordinate_of_discharge_locations'
   fname = 'zdiscoord'
   lname = 'Vertical coordinates of discharge locations'
   unit = 'm'
   nodim = 1
   ngdims(1) = numdis

!
!3.17 Parameters for parallel mode
!---------------------------------
!

CASE (iarr_idprocs)
   sname = 'process_ids'
   fname = 'idprocs'
   lname = 'Process ids'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs

!
!3.18 Energy equation, enstrophy, vorticity
!------------------------------------------
!

CASE (iarr_edens0d)
   sname = 'baroclinic_energy_in_sea_water'
   fname = 'edens0d'
   lname = 'Domain integrated baroclinic energy'
   unit = 'J'
   cnode = ''
   nodim = 0
CASE (iarr_edens2d)
   sname = 'integral_wrt_depth_of_baroclinic_energy__in_sea_water'
   fname = 'edens2d'
   lname = 'Vertically integrated baroclinic energy'
   unit = 'J m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_edens3d)
   sname = 'baroclinic_energy__in_sea_water'
   fname = 'edens3d'
   lname = 'Baroclinic energy'
   unit = 'J m-3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_edissip0d)
   sname = 'energy_dissipation_in_sea_water'
   fname = 'edissip0d'
   lname = 'Domain integrated energy dissipation'
   unit = 'W'
   cnode = ''
   nodim = 0
CASE (iarr_edissip2d)
   sname = 'integral_wrt_depth_of_energy_dissipation_in_sea_water'
   fname = 'edissip2d'
   lname = 'Vertically integrated energy dissipation'
   unit = 'W m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_edissip3d)
   sname = 'energy_dissipation_in_sea_water'
   fname = 'edissip3d'
   lname = 'Energy dissipation'
   unit = 'W m-3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_eflux2du)
   sname = 'integral_wrt_depth_of_x_energy_flux'
   fname = 'eflux2du'
   lname = 'X-component of depth-integrated energy flux'
   vname = 'Depth-integrated energy flux'
   unit = 'W m-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_eflux2dv)
   sname = 'integral_wrt_depth_of_y_energy_flux'
   fname = 'eflux2dv'
   lname = 'Y-component of depth-integrated energy flux'
   vname = 'Depth-integrated energy flux'
   unit = 'W m-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_eflux3du)
   sname = 'x_energy_flux' 
   fname = 'eflux3du'
   lname = 'X-component of energy flux'
   vname = 'Energy flux'
   unit = 'W m-2'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_eflux3dv)
   sname = 'y_energy_flux'
   fname = 'eflux3dv'
   lname = 'Y-component of energy flux'
   vname = 'Energy flux'
   unit = 'W m-2'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_eflux3dw)
   sname = 'upward_energy_flux'
   fname = 'eflux3dw'
   lname = 'Vertical component of energy flux'
   vname = 'Energy flux'
   unit = 'W m-2'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_ekin0d)
   sname = 'kinetic_energy_in_sea_water'
   fname = 'ekin0d'
   lname = 'Domain integrated kinetic energy'
   unit = 'J'
   cnode = ''
   nodim = 0
CASE (iarr_ekin2d)
   sname = 'integral_wrt_depth_of_kinetic_energy_in_sea_water'
   fname = 'ekin2d'
   lname = 'Vertically integrated kinetic energy'
   unit = 'J m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ekin3d)
   sname = 'kinetic_energy_in_sea_water'
   fname = 'ekin3d'
   lname = 'Kinetic energy'
   unit = 'J m-3'
   cnode = ''
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_enstr0d)
   sname = 'enstrophy_in_sea_water'
   fname = 'enstr'
   lname = 'Domain integrated enstrophy'
   unit = 'm3 s-2'
   cnode = ''
   nodim = 0
CASE (iarr_epot0d)
   sname = 'potential_energy_in_sea_water'
   fname = 'epot0d'
   lname = 'Domain integrated potential energy'
   unit = 'J'
   cnode = ''
   nodim = 0
CASE (iarr_epot2d)
   sname = 'integral_wrt_depth_of_potential_energy_in_sea_water'
   fname = 'epot2d'
   lname = 'Vertically integrated potential energy'
   unit = 'J m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_etot0d)
   sname = 'total_energy_in_sea_water'
   fname = 'etot0d'
   lname = 'Domain integrated total energy'
   unit = 'J'
   cnode = ''
   nodim = 0
CASE (iarr_etot2d)
   sname = 'integral_wrt_depth_of_total_energy_in_sea_water'
   fname = 'etot2d'
   lname = 'Vertically integrated total energy'
   unit = 'J m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_etot3d)
   sname = 'total_energy_density_in_sea_water'
   fname = 'etot3d'
   lname = 'Total energy'
   unit = 'J m-3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_vortic2d)
   sname = 'integral_wrt_depth_of_sea_water_relative_vorticity'
   fname = 'vortic2d'
   lname = 'Vertically integrated vorticity'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_vortic3d)
   sname = 'sea_water_relative_vorticity'
   fname = 'vortic3d'
   lname = 'Vorticity'
   unit = 's-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)

!
!3.19 Nesting
!------------
!

CASE (iarr_instbio)
   sname = 'biological_concentrations_for_nesting_per_set'
   fname = 'instbio'
   lname = 'Biological variables used for nesting per set'
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/maxbiovars,nonestsets/)
CASE (iarr_instsed)
   sname = 'sediment_fractions_for_nesting_per_set'
   fname = 'instsed'
   lname = 'Sediment fractions used for nesting per set'
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/maxsedvars,nonestsets/)
CASE (iarr_inst2dtype)
   sname = 'type_2d_open_boundary_data_output'
   fname = 'inst2dtype'
   lname = 'Type 2-D open boundary data output'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_lbhnstatc)
   sname = 'index_of_first_local_nest_c_node_point_in_full_array'
   fname = 'lbhnstatc'
   lname = 'Index of first local nest C-node point in full array'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_lbhnstatu)
   sname = 'index_of_first_local_nest_u_node_point_in_full_array'
   fname = 'lbhnstatu'
   lname = 'Index of first local nest U-node point in full array'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_lbhnstatv)
   sname = 'index_of_first_local_nest_v_node_point_in_full_array'
   fname = 'lbhnstatv'
   lname = 'Index of first local nest V-node point in full array'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_lbhnstatx)
   sname = 'index_of_first_local_nest_x_node_point_in_full_array'
   fname = 'lbhnstatx'
   lname = 'Index of first local nest X-node point in full array'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_lbhnstaty)
   sname = 'index_of_first_local_nest_y_node_point_in_full_array'
   fname = 'lbhnstaty'
   lname = 'Index of first local nest Y-node point in full array'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nestcoords)
   sname = 'type_of_grid_specification'
   fname = 'nestcoords'
   lname = 'Type of grid specification'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1:2) = nonestsets
CASE (iarr_nobionst)
   sname = 'number_nested_biological_concentrations_per_set'
   fname = 'nobionst'
   lname = 'Number of nested biological variables per set'
   dtype = int_type
   cnode =''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstatc)
   sname = 'local_number_of_nest_points_at_c_nodes_per_set'
   fname = 'nohnstatc'
   lname = 'Local number of nest points at C-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstatu)
   sname = 'local_number_of_nest_points_at_u_nodes_per_set'
   fname = 'nohnstatu'
   lname = 'Local number of nest points at U-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstatv)
   sname = 'local_number_of_nest_points_at_v_nodes_per_set'
   fname = 'nohnstatv'
   lname = 'Local number of nest points at V-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstatx)
   sname = 'local_number_of_nest_points_at_x_nodes_per_set'
   fname = 'nohnstatx'
   lname = 'Local number of nest points at X-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstaty)
   sname = 'local_number_of_nest_points_at_y_nodes_per_set'
   fname = 'nohnstaty'
   lname = 'Local number of nest points at Y-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstcprocs)
   sname = 'local_number_of_nest_points_at_c_nodes_per_process_and_set'
   fname = 'nohnstcprocs'
   lname = 'Local number of nest points at C-nodes per process and set'
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nprocs,nonestsets/)
CASE (iarr_nohnstglbc)
   sname = 'number_of_nest_points_at_c_nodes_per_set'
   fname = 'nohnstglbc'
   lname = 'Number of nest points at C-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstglbu)
   sname = 'number_of_nest_points_at_u_nodes_per_set'
   fname = 'nohnstglbu'
   lname = 'Number of nest points at U-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstglbv)
   sname = 'number_of_nest_points_at_v_nodes_per_set'
   fname = 'nohnstglbv'
   lname = 'Number of nest points at V-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstglbx)
   sname = 'number_of_nest_points_at_x_nodes_per_set'
   fname = 'nohnstglbx'
   lname = 'Number of nest points at X-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstglby)
   sname = 'number_of_nest_points_at_y_nodes_per_set'
   fname = 'nohnstglby'
   lname = 'Number of nest points at Y-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstuvprocs)
   sname = 'number_of_local_nest_points_at_u_and_v_nodes_per_process_and_set'
   fname = 'nohnstuvprocs'
   lname = 'Number of local nest points at U- and V-nodes per process and set'
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nprocs,nonestsets/)
CASE (iarr_nohnstxyprocs)
   sname = 'number_of_local_nest_points_at_x_and_y_nodes_per_process_and_set'
   fname = 'nohnstuvprocs'
   lname = 'Number of local nest points at X- and Y-nodes per process and set'
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nprocs,nonestsets/)
CASE (iarr_nosednst)
   sname = 'number_nested_sediment_concentrations_per_set'
   fname = 'nobionst'
   lname = 'Number of nested sediment fractions per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_novnst)
   sname = 'number_of_vertical_nest_points_per_set'
   fname = 'novnst'
   lname = 'Number of vertical nest points per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets

!
!3.20 Relaxation zones
!---------------------
!

CASE (iarr_idirrlx)
   sname = 'orientation_of_relaxation_zones'
   fname = 'idirrlx'
   lname = 'Orientation of relaxation zones'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = norlxzones
CASE (iarr_indexrlxatc)
   sname = 'open_boundary_index_where_an_external_profile_of_a_c_node_'//&
         & 'quantity_is_defined'
   fname = 'indexrlxatc'
   lname = 'Open boundary index where an external profile of a C-node '//&
         & 'quantity is defined'
   dtype = int_type
   cnode = ''
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,2/)
   ngdims(1:3) = (/nc,nr,2/)
CASE (iarr_indexrlxatuv)
   sname = 'open_boundary_index_where_an_external_profile_of_a_u_or_v_node_'//&
         & 'quantity_is_defined'
   fname = 'indexrlxatuv'
   lname = 'Open boundary index where an external profile of a U or V-node '//&
         & 'quantity is defined'
   dtype = int_type
   cnode = ''
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,2/)
   ngdims(1:3) = (/nc,nr,2/)
CASE (iarr_inodesrlx)
   sname = 'type_of_relaxation_zones'
   fname = 'inodesrlx'
   lname = 'Type of relaxation zones'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = 3
CASE (iarr_iposrlx)
   sname = 'x_index_of_sw_corner_of_relaxation_zone'
   fname = 'iposrlx'
   lname = 'X-index of SW corner of relaxation zone'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = norlxzones
CASE (iarr_iprofrlx)
   sname = '' 
   fname = 'iprofrlx'
   lname = ''
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = norlxzones
CASE (iarr_ityprlx)
   sname = 'type_of_relaxation_scheme'
   fname = 'ityprlx'
   lname = 'Type of relaxation scheme'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = norlxzones
CASE (iarr_jposrlx)
   sname = 'y_index_of_sw_corner_of_relaxation_zone'
   fname = 'jposrlx'
   lname = 'Y-index of SW corner of relaxation zone'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = norlxzones
CASE (iarr_ncrlx)
   sname = 'x_dimension_of_relaxation_zone'
   fname = 'ncrlx'
   lname = 'X-dimension of relaxation zone'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = norlxzones
CASE (iarr_nrrlx)
   sname = 'y_dimension_of_relaxation_zone'
   fname = 'nrrlx'
   lname = 'Y-dimension of relaxation zone'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = norlxzones
CASE (iarr_rlxwghtatc)
   sname = 'relaxation_factor_for_c_node_quantities'
   fname = 'rlxwghtatc'
   lname = 'Relaxation factor for C-node quantities'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,2/)
   ngdims(1:3) = (/nc,nr,2/)
CASE (iarr_rlxwghtatuv)
   sname = 'relaxation_factor_for_u_and_v_node_quantities'
   fname = 'rlxwghtatuv'
   lname = 'Relaxation factor for U- and V-node quantities'
   cnode = ''
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,2/)
   ngdims(1:3) = (/nc,nr,2/)

!
!3.21 Elliptic parameters
!------------------------
!

CASE (iarr_ellac2d)
   sname = 'anticyclonic_2d_current'
   fname = 'ellac2d'
   lname = 'Anticyclonic 2-D current'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ellac3d)
   sname = 'anticyclonic_3d_current'
   fname = 'ellac3d'
   lname = 'Anticyclonic 3-D current'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_ellcc2d)
   sname = 'cyclonic_2d_current'
   fname = 'ellcc2d'
   lname = 'Cyclonic 2-D current'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ellcc3d)
   sname = 'cyclonic_3d_current'
   fname = 'ellcc3d'
   lname = 'Cyclonic 3-D current'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_ellinc2d)
   sname = 'inclination_of_2d_tidal_ellipses'
   fname = 'ellinc2d'
   lname = 'Inclination of 2-D tidal ellipses '
   unit = MERGE('degrees','radian ',DegreesOut)
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ellinc3d)
   sname = 'inclination_of_3d_tidal_ellipses'
   fname = 'ellinc3d'
   lname = 'Inclination of 3-D tidal ellipses'
   unit = MERGE('degrees','radian ',DegreesOut)
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_ellip2d)
   sname = 'ellipticity_of_2d_tidal_ellipses'
   fname = 'ellip2d'
   lname = 'Ellipticity of 2-D tidal ellipses'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ellip3d)
   sname = 'ellipticity_of_3d_tidal_ellipses'
   fname = 'ellip3d'
   lname = 'Ellipticity of 3-D tidal ellipses'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_ellmaj2d)
   sname = 'major_axis_of_2d_tidal_ellipses'
   fname = 'ellmaj2d'
   lname = 'Major axis of 2-D tidal ellipses'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ellmaj3d)
   sname = 'major_axis_of_3d_tidal_ellipses'
   fname = 'ellmaj3d'
   lname = 'Major axis of 3-D tidal ellipses'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_ellmin2d)
   sname = 'minor_axis_of_2d_tidal_ellipses'
   fname = 'ellmin2d'
   lname = 'Minor axis of 2-D tidal ellipses'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ellmin3d)
   sname = 'minor_axis_of_3d_tidal_ellipses'
   fname = 'ellmin3d'
   lname = 'Minor axis of 3-D tidal ellipses'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_ellpha2d)
   sname = 'elliptic_phase_of_2d_tidal_ellipses'
   fname = 'ellpha2d'
   lname = 'Elliptic phase of 2-D tidal ellipses'
   unit = MERGE('degrees','radian ',DegreesOut)
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ellpha3d)
   sname = 'elliptic_phase_of_3d_tidal_ellipses'
   fname = 'ellpha3d'
   lname = 'Elliptic phase of 3-D tidal ellipses'
   unit = MERGE('degrees','radian ',DegreesOut)
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)

!
!3.22 Relative coordinates
!-------------------------
!
!---i-coordinate
CASE (iarr_icoordC)
   sname = 'i_interpolation_index_coordinate_from_external_to_c_node_grid'
   fname = 'icoordC'
   lname = 'I interpolation index coordinate from external to C-node grid'
   dtype = int_type
CASE (iarr_icoordCC)
   sname = 'i_interpolation_index_coordinate_from_external_c_node_to_model_'//&
         & 'c_node_grid' 
   fname = 'icoordCC'
   lname = 'I interpolation index coordinate from external C-node to model '//&
         & 'C-node grid'
   dtype = int_type
CASE (iarr_icoordCU)
   sname = 'i_interpolation_index_coordinate_from_external_c_node_to_model_'//&
         & 'u_node_grid'   
   fname = 'icoordCU'
   lname = 'I interpolation index coordinate from external C-node to model '//&
        & 'U-node grid'
   dtype = int_type
CASE (iarr_icoordCV)
   sname = 'i_interpolation_index_coordinate_from_external_c_node_to_model_'//&
         & 'v_node_grid'
   fname = 'icoordCV'
   lname = 'I interpolation index coordinate from external C-node to model '//&
         & 'V-node grid'
   dtype = int_type
CASE (iarr_icoordU)
   sname = 'i_interpolation_index_coordinate_from_external_to_u_node_grid'
   fname = 'icoordU'
   lname = 'I interpolation index coordinate from external to U-node grid'
   dtype = int_type
CASE (iarr_icoordUU)
   sname = 'i_interpolation_index_coordinate_from_external_u_node_to_model_'//&
         & 'u_node_grid'
   fname = 'icoordUU'
   lname = 'I interpolation index coordinate from external U-node to model '//&
         & 'U-node grid'
   dtype = int_type
CASE (iarr_icoordUY)
   sname = 'i_interpolation_index_coordinate_from_external_u_node_to_model_'//&
         & 'y_node_grid'
   fname = 'icoordUY'
   lname = 'I interpolation index coordinate from external U-node to model '//&
         & 'Y-node grid'
   dtype = int_type
CASE (iarr_icoordV)
   sname = 'i_interpolation_index_coordinate_from_external_to_v_node_grid'
   fname = 'icoordV'
   lname = 'I interpolation index coordinate from external to V-node grid'
   dtype = int_type
CASE (iarr_icoordVV)
   sname = 'i_interpolation_index_coordinate_from_external_v_node_to_model_'//&
         & 'v_node_grid'
   fname = 'icoordVV'
   lname = 'I interpolation index coordinate from external V-node to model '//&
         & 'V-node grid'
   dtype = int_type
CASE (iarr_icoordVX)
   sname = 'i_interpolation_index_coordinate_from_external_v_node_to_model_'//&
         & 'x_node_grid'
   fname = 'icoordVX'
   lname = 'I interpolation index coordinate from external V-node to model '//&
         & 'X-node grid'
   dtype = int_type

!---j-coordinate
CASE (iarr_jcoordC)
   sname = 'j_interpolation_index_coordinate_from_external_to_c_node_grid'
   fname = 'jcoordC'
   lname = 'J interpolation index coordinate from external to C-node grid'
   dtype = int_type
CASE (iarr_jcoordCC)
   sname = 'j_interpolation_index_coordinate_from_external_c_node_to_model_'//&
        & 'c_node_grid'
   fname = 'jcoordCC'
   lname = 'J interpolation index coordinate from external C-node to model '//&
         & 'C-node grid'
   dtype = int_type
CASE (iarr_jcoordCU)
   sname = 'j_interpolation_index_coordinate_from_external_c_node_to_model_'//&
         & 'u_node_grid'
   fname = 'jcoordCU'
   lname = 'J interpolation index coordinate from external C-node to model '//&
         & 'U-node grid'
   dtype = int_type
CASE (iarr_jcoordCV)
   sname = 'j_interpolation_index_coordinate_from_external_c_node_to_model_'//&
         & 'v_node_grid'
   fname = 'jcoordCV'
   lname = 'J interpolation index coordinate from external C-node to model '//&
         & 'V-node grid'
   dtype = int_type
CASE (iarr_jcoordU)
   sname = 'j_interpolation_index_coordinate_from_external_to_u_node_grid'
   fname = 'jcoordU'
   lname = 'J interpolation index coordinate from external to U-node grid'
   dtype = int_type
CASE (iarr_jcoordUU)
   sname = 'j_interpolation_index_coordinate_from_external_u_node_to_model_'//&
         & 'u_node_grid'
   fname = 'jcoordUU'
   lname = 'J interpolation index coordinate from external U-node to model '//&
         & 'U-node grid'
   dtype = int_type
CASE (iarr_jcoordUY)
   sname = 'j_interpolation_index_coordinate_from_external_u_node_to_model_'//&
         & 'y_node_grid'
   fname = 'jcoordUY'
   lname = 'J interpolation index coordinate from external U-node to model '//&
         & 'Y-node grid'
   dtype = int_type
CASE (iarr_jcoordV)
   sname = 'j_interpolation_index_coordinate_from_external_to_v_node_grid'
   fname = 'jcoordV'
   lname = 'J interpolation index coordinate from external to V-node grid'
   dtype = int_type
CASE (iarr_jcoordVV)
   sname = 'j_interpolation_index_coordinate_from_external_v_node_to_model_'//&
         & 'v_node_grid'
   fname = 'jcoordVV'
   lname = 'J interpolation index coordinate from external V-node to model '//&
         & 'V-node grid'
   dtype = int_type
CASE (iarr_jcoordVX)
   sname = 'j_interpolation_index_coordinate_from_external_v_node_to_model_'//&
         & 'x_node_grid'
   fname = 'jcoordVX'
   lname = 'J interpolation index coordinate from external V-node to model '//&
         & 'X-node grid'
   dtype = int_type

!---x-coordinate
CASE (iarr_weightsC)
   sname = 'weights_for_interpolation_from_external_to_c_node_grid'
   fname = 'weightsC'
   lname = 'Weights for interpolation from external to C-node grid'
CASE (iarr_weightsCC)
   sname = 'weights_for_interpolation_from_external_c_node_to_model_c_node_grid'
   fname = 'weightsCC'
   lname = 'Weights for interpolation from external C-node to model C-node grid'
CASE (iarr_weightsCU)
   sname = 'weights_for_interpolation_from_external_c_node_to_model_u_node_grid'
   fname = 'weightsCU'
   lname = 'Weights for interpolation from external C-node to model U-node grid'
CASE (iarr_weightsCV)
   sname = 'weights_for_interpolation_from_external_c_node_to_model_v_node_grid'
   fname = 'weightsCV'
   lname = 'Weights for interpolation from external C-node to model V-node grid'
CASE (iarr_weightsU)
   sname = 'weights_for_interpolation_from_external_to_u_node_grid'
   fname = 'weightsU'
   lname = 'Weights for interpolation from external to U-node grid'
CASE (iarr_weightsUU)
   sname = 'weights_for_interpolation_from_external_u_node_to_model_u_node_grid'
   lname = 'Weights for interpolation from external U-node to model U-node grid'
   fname = 'weightsUU'
CASE (iarr_weightsUY)
   sname = 'weights_for_interpolation_from_external_u_node_to_model_y_node_grid'
   fname = 'weightsUY'
   lname = 'Weights for interpolation from external U-node to model Y-node grid'
CASE (iarr_weightsV)
   sname = 'weights_for_interpolation_from_external_to_v_node_grid'
   fname = 'weightsV'
   lname = 'Weights for interpolation from external to V-node grid'
CASE (iarr_weightsVV)
   sname = 'weights_for_interpolation_from_external_v_node_to_model_v_node_grid'
   fname = 'weightsVV'
   lname = 'Weights for interpolation from external V-node to model V-node grid'
CASE (iarr_weightsVX)
   sname = 'weights_for_interpolation_from_external_v_node_to_model_x_node_grid'
   fname = 'weightsVX'
   lname = 'Weights for interpolation from external V-node to model X-node grid'

!
!3.23 Data coordinates
!---------------------
!

CASE (iarr_depout)
   sname = 'sea_floor_depth_below_sea_level'
   commt = ''
   fname = 'depout'
   lname = 'Mean water depth'
   unit = 'm'
   cnode = ''
CASE (iarr_levout)
   IF (iopt_grid_vtype.EQ.1) THEN
      fname = 'lev'
      sname = 'ocean_sigma_coordinate'
      lname = 'Sigma coordinate'
      nodim = 1
   ELSEIF (iopt_grid_vtype.EQ.2) THEN
      fname = 'lev'
      sname = 'ocean_s_coordinate'
      lname = 's-coordinate'
      nodim = 1
   ELSE
      fname = 'gsig'
      sname = 'ocean_s_coordinate_g'
      lname = 'General s-coordinate'
      nodim = 3
   ENDIF
CASE (iarr_time)
   sname = 'time'
   commt = ''
   fname = 'time'
   lname = 'Time'
   unit = 'date/time'
   dtype = char_type
   cnode = ''
CASE (iarr_xcoord)
   sname = MERGE('projection_x_coordinate','longitude              ',&
               & iopt_grid_sph.EQ.0)
   commt = ''
   fname =  'xcoord'
   lname =  MERGE('X-coordinate','longitude   ',iopt_grid_sph.EQ.0)
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = ''
CASE (iarr_xcoordatc)
   sname = MERGE('projection_x_coordinate_at_c_nodes',&
               & 'longitude                         ',iopt_grid_sph.EQ.0)
   commt = ''
   fname =  'xcoordatc'
   lname =  MERGE('X-coordinate at C-nodes','longitude at C-nodes   ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
CASE (iarr_xcoordatu)
   sname = MERGE('projection_x_coordinate_at_u_nodes',&
                & 'longitude_at_u_nodes              ',iopt_grid_sph.EQ.0)
   fname =  'xcoordatu'
   lname =  MERGE('X-coordinate at U-nodes','longitude at U-nodes   ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = 'U'
CASE (iarr_xcoordatv)
   sname = MERGE('projection_x_coordinate_at_v_nodes',&
               & 'longitude_at_v_nodes              ',iopt_grid_sph.EQ.0)
   fname =  'xcoordatv'
   lname =  MERGE('X-coordinate at V-nodes','longitude at V-nodes   ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = 'V'
CASE (iarr_xcoordatx)
   sname = MERGE('projection_x_coordinate_at_x_nodes',&
               & 'longitude_at_x_nodes              ',iopt_grid_sph.EQ.0)
   fname =  'xcoordatx'
   lname =  MERGE('X-coordinate at X-nodes','longitude at X-nodes   ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = 'X'
CASE (iarr_xcoordaty)
   sname = MERGE('projection_x_coordinate_at_y_nodes',&
               & 'longitude_at_y_nodes              ',iopt_grid_sph.EQ.0)
   fname =  'xcoordaty'
   lname =  MERGE('X-coordinate at Y-nodes','longitude at Y-nodes   ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = 'Y'
CASE (iarr_xout)
   sname =  MERGE('projection_x_coordinate','longitude              ',&
                & iopt_grid_sph.EQ.0)
   commt = ''
   fname =  MERGE('x  ','lon',iopt_grid_sph.EQ.0)
   lname =  MERGE('X-coordinate','longitude   ',iopt_grid_sph.EQ.0)
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = ''
CASE (iarr_ycoord)
   sname =  MERGE('projection_y_coordinate','latitude               ',&
                & iopt_grid_sph.EQ.0)
   commt = ''
   fname =  'ycoord'
   lname =  MERGE('Y-coordinate','latitude    ',iopt_grid_sph.EQ.0)
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = ''
CASE (iarr_ycoordatc)
   sname =  MERGE('projection_y_coordinate_at_c_nodes',&
                & 'latitude                          ',iopt_grid_sph.EQ.0)
   commt = ''
   fname =  'ycoordatc'
   lname =  MERGE('Y-coordinate at C-nodes','latitude at C-nodes    ',&
                 & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
CASE (iarr_ycoordatu)
   sname =  MERGE('projection_y_coordinate_at_u_nodes',&
                & 'latitude_at_u_nodes               ',iopt_grid_sph.EQ.0)
   fname =  'ycoordatu'
   lname =  MERGE('Y-coordinate at U-nodes','latitude at U-nodes    ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = 'U'
CASE (iarr_ycoordatv)
   sname =  MERGE('projection_y_coordinate_at_v_nodes',&
                & 'latitude_at_v_nodes               ',iopt_grid_sph.EQ.0)
   fname =  'ycoordatv'
   lname =  MERGE('Y-coordinate at V-nodes','latitude at V-nodes    ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = 'V'
CASE (iarr_ycoordatx)
   sname =  MERGE('projection_y_coordinate_at_x_nodes',&
                & 'latitude_at_x_nodes               ',iopt_grid_sph.EQ.0)
   fname =  'ycoordatx'
   lname =  MERGE('Y-coordinate at X-nodes','latitude at X-nodes    ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = 'X'
CASE (iarr_ycoordaty)
   sname =  MERGE('projection_y_coordinate_at_Y_nodes',&
                & 'latitude_at_y_nodes               ',iopt_grid_sph.EQ.0)
   fname =  'ycoordaty'
   lname =  MERGE('Y-coordinate at Y-nodes','latitude at Y-nodes    ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = 'Y'
CASE (iarr_yout)
   sname =  MERGE('projection_y_coordinate','latitude               ',&
                & iopt_grid_sph.EQ.0)
   commt = ''
   fname =  MERGE('y  ','lat',iopt_grid_sph.EQ.0)
   lname =  MERGE('Y-coordinate','latitude    ',iopt_grid_sph.EQ.0)
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = ''
CASE (iarr_zcmean)
   sname = 'ocean_z_coordinate_with_respect_to_mean_sea_level'
   fname = 'zout'
   lname = 'Z-coordinate'
   unit = 'm'
   cnode = ''
CASE (iarr_zcoord)
   sname = 'ocean_z_coordinate'
   fname = 'zout'
   lname = 'Z-coordinate'
   unit = 'm'
   cnode = ''
CASE (iarr_zcoordatc)
   sname = 'ocean_z_coordinate_at_z_nodes'
   fname = 'zcoordatc'
   lname = 'Z-coordinate at C-nodes'
   unit = 'm'
CASE (iarr_zcoordatu)
   sname = 'ocean_z_coordinate_at_u_nodes'
   fname = 'zcoordatu'
   lname = 'Z-coordinate at U-nodes'
   unit = 'm'
   cnode = 'U'
CASE (iarr_zcoordatv)
   sname = 'ocean_z_coordinate_at_v_nodes'
   fname = 'zcoordatv'
   lname = 'Z-coordinate at V-nodes'
   unit = 'm'
   cnode = 'V'
CASE (iarr_zcoordatx)
   sname = 'ocean_z_coordinate_at_x_nodes'
   fname = 'zcoordatx'
   lname = 'Z-coordinate at X-nodes'
   unit = 'm'
   cnode = 'X'
CASE (iarr_zcoordaty)
   sname = 'ocean_z_coordinate_at_y_nodes'
   fname = 'zcoordaty'
   lname = 'Z-coordinate at Y-nodes'
   unit = 'm'
   cnode = 'Y'
CASE (iarr_zetout)
   sname = 'sea_surface_height_above_sea_level'
   fname = 'zetout'
   lname = 'Surface elevation'
   unit = 'm'
   
!
!3.24 Model parameters
!---------------------
!

CASE (iarr_density_ref)
   sname = 'sea_water_reference_density '
   fname = 'density_ref'
   lname = 'Reference density'
   unit = 'kg m-3'
   nodim = 0
CASE (iarr_gacc_mean)
   sname = 'acceleration_due_to_gravity'
   fname = 'gacc_mean'
   lname = 'Mean acceleration of gravity'
   unit = 'm s-2'
   nodim = 0
END SELECT

!
!4. Store
!--------
!

IF (PRESENT(f90_name)) f90_name = fname
IF (PRESENT(standard_name)) standard_name = sname
IF (PRESENT(comment)) comment = commt
IF (PRESENT(long_name)) long_name = lname
IF (PRESENT(units)) units = unit
IF (PRESENT(node)) node = cnode
IF (PRESENT(vector_name)) vector_name = vname
IF (PRESENT(data_type)) data_type = dtype
IF (PRESENT(fill_value)) fill_value = fvalue
IF (PRESENT(nrank)) nrank = nodim
IF (PRESENT(global_dims)) global_dims = ngdims
IF (PRESENT(local_dims)) local_dims = nldims
IF (PRESENT(halo_dims)) halo_dims = nhdims
IF (PRESENT(varatts)) THEN
   varatts%standard_name = sname
   varatts%comment = commt
   varatts%f90_name = fname
   varatts%long_name = lname
   varatts%vector_name = vname
   varatts%units = unit
   varatts%data_type = dtype
   varatts%nrank = nodim
   varatts%global_dims = ngdims
   varatts%local_dims = nldims
   varatts%halo_dims = nhdims
   varatts%fill_value = float_fill
ENDIF


RETURN

END SUBROUTINE inquire_var

!========================================================================

SUBROUTINE set_modfiles_atts(iddesc,ifil,iotype)
!************************************************************************
!
! *set_modfiles_atts* Obtain the global attributes of a model file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modvars_routines.f90  V2.11
!
! Description -
!
! Module calls - inquire_var, set_biofiles_atts, set_modfiles_name,
!                set_partfiles_atts, set_sedfiles_atts
!
!************************************************************************
!
USE gridpars
USE nestgrids
USE obconds
USE relaxation
USE structures
USE switches
USE tide
USE biovars_routines, ONLY: set_biofiles_atts
USE sedvars_routines, ONLY: set_sedfiles_atts
USE partvars_routines, ONLY: set_partfiles_atts

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
CHARACTER (LEN=12) :: cfil
CHARACTER (LEN=lendesc) :: title
INTEGER :: nocoords, nofiles, novars


pglev = pglev + 1
procname(pglev) = 'set_modfiles_atts'
IF (loglev1.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, TRIM(procname(pglev))
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, TRIM(procname(pglev))
   ENDIF
ENDIF

!
!1. Initialise parameters
!------------------------
!

nocoords = 0
nofiles = maxdatafiles(iddesc,1)
WRITE (cfil,'(I12)') ifil; cfil = ADJUSTL(cfil)

!
!2. Define attributes
!--------------------
!

SELECT CASE (iddesc)

!
!2.1  Domain decomposition
!------------------------
!

CASE (io_mppmod)
   title = 'Domain decomposition'
   novars = 4

!
!2.2 Initial conditions
!----------------------
!

CASE (io_inicon,io_fincon)

   SELECT CASE (ifil)
      CASE (ics_phys)
         title = 'Physical initial conditions'
         novars = 0
         IF (iopt_mode_2D.GT.0) novars = novars + 3
         IF (iopt_mode_3D.GT.0) novars = novars + 2
         IF (iopt_mode_3D.GT.0.AND.iopt_grid_nodim.EQ.3) novars = novars + 1
         IF (iopt_temp.GT.0) novars = novars + 1
         IF (iopt_sal.GT.0) novars = novars + 1
         IF (iopt_vdif_coef.EQ.3) THEN
            novars = novars + 1
            IF (iopt_turb_ntrans.EQ.2.AND.iopt_turb_param.GT.0) THEN
               novars = novars + 1
            ENDIF
         ENDIF
         IF (iopt_bstres_form.EQ.2) THEN
            IF (iopt_bstres_drag.EQ.2.OR.iopt_bstres_drag.EQ.4) THEN
               novars = novars + 1
            ENDIF
         ENDIF
         IF (nconobc.GT.0) novars = novars + 2
         IF (iopt_astro_tide.EQ.1) novars = novars + 2
         IF (nobu.GT.0) THEN
            IF (iopt_obc_sal.EQ.1) novars = novars + 1
            IF (iopt_obc_temp.EQ.1) novars = novars + 1
            IF (iopt_obc_2D.EQ.1) novars = novars + 1
            IF (iopt_obc_3D.EQ.1) novars = novars + 1
         ENDIF
         IF (nobv.GT.0) THEN
            IF (iopt_obc_sal.EQ.1) novars = novars + 1
            IF (iopt_obc_temp.EQ.1) novars = novars + 1
            IF (iopt_obc_2D.EQ.1) novars = novars + 1
            IF (iopt_obc_3D.EQ.1) novars = novars + 1
         ENDIF
         nocoords = 1
     CASE (ics_sed,ics_morph)
        CALL set_sedfiles_atts(iddesc,ifil,iotype)
        GOTO 1000
     CASE (ics_bio)
        CALL set_biofiles_atts(iddesc,ifil,iotype)
        GOTO 1000
     CASE (ics_part)
        CALL set_partfiles_atts(iddesc,ifil,iotype)
        GOTO 1000
   END SELECT

!
!2.3 Model grid and data locations
!---------------------------------
!

CASE (io_modgrd)
   title = 'Model grid'
   novars = MERGE(1,3,iopt_grid_htype.EQ.1)
   IF (iopt_grid_vtype.GT.1) novars = novars + 1
   IF (nobu.GT.0) novars = novars + 2
   IF (nobv.GT.0) novars = novars + 2

CASE (io_metabs)
   title = 'Absolute coordinates of meteo grid'
   novars = 2

CASE (io_sstabs)
   title = 'Absolute coordinates of SST grid'
   novars = 2

CASE (io_wavabs)
   title = 'Absolute coordinates of wave grid'
   novars = 2

CASE (io_nstabs)
   title = 'Nested absolute grid locations: '//TRIM(cfil)
   novars = 0
   IF (nohnstglbc(ifil).GT.0) THEN
      novars = novars + 2
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF
   IF (nohnstglbu(ifil).GT.0) THEN
      novars = novars + 2
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF
   IF (nohnstglbv(ifil).GT.0) THEN
      novars = novars + 2
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF
   IF (nohnstglbx(ifil).GT.0) THEN
      novars = novars + 2
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF
   IF (nohnstglby(ifil).GT.0) THEN
      novars = novars + 2
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF

CASE (io_metrel)
   title = 'Relative coordinates of meteo grid'
   novars = 3

CASE (io_sstrel)
   title = 'Relative coordinates of SST grid'
   novars = 3

CASE (io_wavrel)
   title = 'Relative coordinates of wave grid'
   novars = 3

CASE (io_nstrel)
   title = 'Nested relative grid locations: '//TRIM(cfil)
   novars = 0
   IF (nohnstglbc(ifil).GT.0) THEN
      novars = novars + 3
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF
   IF (nohnstglbu(ifil).GT.0) THEN
      novars = novars + 6
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF
   IF (nohnstglbv(ifil).GT.0) THEN
      novars = novars + 6
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF
   IF (nohnstglbx(ifil).GT.0) THEN
      novars = novars + 3
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF
   IF (nohnstglby(ifil).GT.0) THEN
      novars = novars + 3
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF

!
!2.4 Open boundaries
!-------------------
!

CASE (io_1uvsur)
   IF (ifil.EQ.1) THEN
      title = 'Type of boundary conditions for 1-D mode'
      IF (isur1dtype.EQ.1.OR.isur1dtype.EQ.3) novars = 5
      IF (isur1dtype.EQ.1.OR.isur1dtype.EQ.2) novars = 3
   ELSEIF (ifil.EQ.2) THEN
      title = 'Boundary data for 1-D mode'
      SELECT CASE (isur1dtype)
         CASE (1); novars = 3
         CASE (2); novars = 1
         CASE (3); novars = 2
      END SELECT
      nocoords = 1
   ENDIF

CASE (io_2uvobc)
   IF (ifil.EQ.1) THEN
      title = 'Type of open boundary conditions for 2-D mode'
      novars = 0
      IF (nobu.GT.0) THEN
         novars = novars + MERGE(6,2,nconobc.GT.0)
      ENDIF
      IF (nobv.GT.0) THEN
         novars = novars + MERGE(6,2,nconobc.GT.0)
      ENDIF
      IF (nofiles.GT.1) novars = novars + 3
      IF (nqsecobu.GT.0) novars = novars + 1
      IF (nqsecobv.GT.0) novars = novars + 1
   ELSE
      title = 'Open boundary data for 2-D mode: '//TRIM(cfil)
      novars = MERGE(1,2,iobc2dtype(ifil).GT.1)
      nocoords = 1
   ENDIF

CASE (io_2xyobc)
   IF (ifil.EQ.1) THEN
      title = 'Type of tangential open boundary conditions for 2-D mode'
      novars = 0
      IF (nobx.GT.0) THEN
         novars = novars + MERGE(3,1,nconobc.GT.0)
      ENDIF
      IF (noby.GT.0) THEN
         novars = novars + MERGE(3,1,nconobc.GT.0)
      ENDIF
      IF (nofiles.GT.1) novars = novars + 2
   ELSE
      title = 'Tangential open boundary data for 2-D mode: '//TRIM(cfil)
      novars = 1
      nocoords = 1
   ENDIF

CASE (io_3uvobc)
   IF (ifil.EQ.1) THEN
      title = 'Type of open boundary conditions for 3-D current'
      novars = 0
      IF (nobu.GT.0) novars = novars + 2
      IF (nobv.GT.0) novars = novars + 2
      IF (nofiles.GT.1) novars = novars + 2*nofiles - 1
      IF (norlxzones.GT.0) novars = novars + 1
   ELSE
      title = 'Open boundary data for 3-D current: '//TRIM(cfil)
      novars = 1; nocoords = 1
   ENDIF

CASE (io_3xyobc)
   IF (ifil.EQ.1) THEN
      title = 'Type of tangential open boundary conditions for 3-D current'
      novars = 0
      IF (nobx.GT.0) novars = novars + 2
      IF (noby.GT.0) novars = novars + 2
      IF (nofiles.GT.1) novars = novars + 2*nofiles - 1
      IF (norlxzones.GT.0) novars = novars + 1
   ELSE
      title = 'Tangential open boundary data for 3-D current: '//TRIM(cfil)
      novars = 1; nocoords = 1
   ENDIF

CASE (io_salobc)
   IF (ifil.EQ.1) THEN
      title = 'Type of open boundary conditions for salinity'
      novars = 0
      IF (nobu.GT.0) novars = novars + 2
      IF (nobv.GT.0) novars = novars + 2
      IF (nofiles.GT.1) novars = novars + 2*nofiles - 1
      IF (norlxzones.GT.0) novars = novars + 1
   ELSE
      title = 'Open boundary data for salinity: '//TRIM(cfil)
      novars = 2; nocoords = 1
   ENDIF

CASE (io_tmpobc)
   IF (ifil.EQ.1) THEN
      title = 'Type of open boundary conditions for temperature'
      novars = 0
      IF (nobu.GT.0) novars = novars + 2
      IF (nobv.GT.0) novars = novars + 2
      IF (nofiles.GT.1) novars = novars + 2*nofiles - 1
      IF (norlxzones.GT.0) novars = novars + 1
   ELSE
      title = 'Open boundary data for temperature: '//TRIM(cfil)
      novars = 1; nocoords = 1
   ENDIF

CASE (io_sedobc)
   IF (ifil.EQ.1) THEN
      title = 'Type of open boundary conditions for sediments'
      novars = 0
      IF (nobu.GT.0) novars = novars + 2
      IF (nobv.GT.0) novars = novars + 2
      IF (nofiles.GT.1) novars = novars + 2*nofiles - 1
      IF (norlxzones.GT.0) novars = novars + 1
   ELSE
      title = 'Open boundary data for sediments: '//TRIM(cfil)
      novars = 1; nocoords = 1
   ENDIF

CASE (io_bioobc)
   IF (ifil.EQ.1) THEN
      title = 'Type of open boundary conditions for biology'
      novars = 0
      IF (nobu.GT.0) novars = novars + 2
      IF (nobv.GT.0) novars = novars + 2
      IF (nofiles.GT.1) novars = novars + 2*nofiles - 1
      IF (norlxzones.GT.0) novars = novars + 1
   ELSE
      title = 'Open boundary data for biology: '//TRIM(cfil)
      novars = 1; nocoords = 1
   ENDIF

CASE (io_rlxobc)
   title = 'Locations of relaxation zones' 
   novars = 7

!
!2.5 Nested output
!-----------------
!

CASE (io_nstspc)
   title = 'Number of nested grid locations'
   novars = 8

CASE (io_2uvnst)
   title = 'Open boundary data for 2-D mode: '//TRIM(cfil)
   novars = MERGE(2,1,inst2dtype(ifil).EQ.1)
   nocoords = 1

CASE (io_3uvnst)
   title = 'Open boundary data for 3-D current: '//TRIM(cfil)
   novars = 1; nocoords = 1

CASE (io_2xynst)
   title = 'Open boundary data for 2-D tangential mode: '//TRIM(cfil)
   novars = 1; nocoords = 1

CASE (io_3xynst)
   title = 'Open boundary data for 3-D tangential current: '//TRIM(cfil)
   novars = 1; nocoords = 1

CASE (io_salnst)
   title = 'Open boundary data for salinity: '//TRIM(cfil)
   novars = 1; nocoords = 1

CASE (io_tmpnst)
   title = 'Open boundary data for temperature: '//TRIM(cfil)
   novars = 1; nocoords = 1

CASE (io_sednst)
   title = 'Open boundary data for sediments: '//TRIM(cfil)
   novars = 1; nocoords = 1

CASE (io_bionst)
   title = 'Open boundary data for biology: '//TRIM(cfil)
   novars = 1; nocoords = 1

!
!2.6 Surface data
!----------------
!

CASE (io_metsur)
   IF (ifil.EQ.1) THEN
      title = 'Meteo input data'
      novars = 0
      IF (iopt_meteo_stres.EQ.1) novars = novars + 2
      IF (iopt_meteo_pres.EQ.1) novars = novars + 1
      IF (iopt_meteo_heat.EQ.1) THEN
         SELECT CASE (iopt_meteo_data)
            CASE (1); novars = novars + 3
            CASE (2); novars = novars + 2
         END SELECT
      ENDIF
      IF (iopt_meteo_precip.EQ.1) novars = novars + 1
   ELSEIF (ifil.EQ.2) THEN
      title = 'Meteo output data'
      novars = 0
   ENDIF
   nocoords = 1

CASE (io_sstsur)
   title = 'SST input data'
   novars = 1; nocoords = 1

CASE (io_wavsur)
   novars = 3
   IF (iopt_waves_form.EQ.2) THEN
      novars = novars + 2
      IF (iopt_waves_curr.EQ.1) THEN
         novars = novars + 2
         IF (iopt_waves_pres.EQ.1) novars = novars + 1
      ENDIF
   ENDIF
   IF (iopt_waves_dissip.EQ.1) novars = novars + 4
   IF (ifil.EQ.1) THEN
      title = 'Surface wave input data'
   ELSEIF (ifil.EQ.2) THEN
      title = 'Surface wave output data'
   ENDIF
   nocoords = 1

!
!2.7 Structures
!--------------
!

CASE (io_drycel)
   title = 'Dry cells'
   novars = 2

CASE (io_thndam)
   title = 'Thin dams'
   novars = 0
   IF (numthinu.GT.0) novars = novars + 2
   IF (numthinv.GT.0) novars = novars + 2

CASE (io_weibar)
   title = 'Weirs/barriers'
   novars = 0
   IF (numwbaru.GT.0) novars = novars + 8
   IF (numwbarv.GT.0) novars = novars + 8

!
!2.8 Discharges
!--------------
!

CASE (io_disspc)
   title = 'Discharge location type and flagging'
   novars = 1

CASE (io_disloc)
   title = 'Discharge locations'
   novars = 3; nocoords = 1

CASE (io_disvol)
   title = 'Volume discharge rates'
   novars = 1; nocoords = 1

CASE (io_discur)
   title = 'Discharge area and direction'
   novars = 2; nocoords = 1

CASE (io_dissal)
   title = 'Salinity discharge'
   novars = 1; nocoords = 1

CASE (io_distmp)
   title = 'Temperature discharge'
   novars = 1; nocoords = 1

!
!2.9 Sediment model files
!------------------------
!

CASE (io_sedspc,io_darspc)
   CALL set_sedfiles_atts(iddesc,ifil,iotype)
   GOTO 1000

!
!2.10 Biological model files
!---------------------------
!

CASE (io_bioabs,io_biorel,io_biospc,io_biosur)
   CALL set_biofiles_atts(iddesc,ifil,iotype)
   GOTO 1000

!
!2.11 Particle model files
!-------------------------
!   

CASE (io_parcld,io_pargrd,io_parphs,io_parspc)
   CALL set_partfiles_atts(iddesc,ifil,iotype)
   GOTO 1000
   
END SELECT

!
!3. Store attributes
!-------------------
!

CALL set_modfiles_name(iddesc,ifil,iotype)

modfiles(iddesc,ifil,iotype)%title = title
modfiles(iddesc,ifil,iotype)%novars = novars
modfiles(iddesc,ifil,iotype)%nocoords = nocoords

1000 CONTINUE
IF (loglev2.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, logexit
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, logexit
   ENDIF
ENDIF
pglev = pglev - 1


RETURN

END SUBROUTINE set_modfiles_atts

!========================================================================

SUBROUTINE set_modfiles_name(iddesc,ifil,iotype)
!************************************************************************
!
! *set_modfiles_name* Obtain the name of a model file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modvars_routines.f90  V2.11
!
! Description -
!
! Module calls - file_suffix
!
!************************************************************************
!
USE structures

!* Arguments
!
INTEGER, INTENT(IN) :: iddesc,ifil, iotype

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
CHARACTER (LEN=3) :: suffix
CHARACTER (LEN=12) :: cfil
CHARACTER (LEN=leniofile) :: namedesc


pglev = pglev + 1
procname(pglev) = 'set_modfiles_name'
IF (loglev1.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, TRIM(procname(pglev))
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, TRIM(procname(pglev))
   ENDIF
ENDIF

!
!1. Initialise parameters
!------------------------
!

WRITE (cfil,'(I12)') ifil; cfil = ADJUSTL(cfil)

!
!2. Define attributes
!--------------------
!

SELECT CASE (iddesc)

!
!2.1 Domain decomposition
!------------------------
!

CASE (io_mppmod)
   IF (ifil.EQ.1) THEN
      namedesc = 'mppmod'
   ELSE
      namedesc = 'mppmod'//TRIM(cfil)
   ENDIF

!
!2.2 Initial/final conditions
!----------------------------
!

CASE (io_inicon)
   SELECT CASE (ifil)
      CASE (ics_phys)
         namedesc = 'phsics'
      CASE (ics_sed)
         namedesc = 'sedics'
      CASE (ics_morph)
         namedesc = 'morics'
      CASE (ics_bio)
         namedesc = 'bioics'
      CASE (ics_part)
         namedesc = 'parics'
   END SELECT

CASE (io_fincon)
   SELECT CASE (ifil)
      CASE (ics_phys)
         namedesc = 'phsfin'
      CASE (ics_sed)
         namedesc = 'sedfin'
      CASE (ics_morph)
         namedesc = 'morfin'
      CASE (ics_bio)
         namedesc = 'biofin'
      CASE (ics_part)
         namedesc = 'parfin'
   END SELECT

!
!2.3 Model grid and data locations
!---------------------------------
!

CASE (io_modgrd)
   namedesc = 'modgrd'
CASE (io_metabs)
   namedesc = 'metabs'
CASE (io_sstabs)
   namedesc = 'sstabs'
CASE (io_wavabs)
   namedesc = 'wavabs'
CASE (io_nstabs)
   namedesc = 'nstabs'//TRIM(cfil)
CASE (io_metrel)
   namedesc = 'metrel'
CASE (io_sstrel)
   namedesc = 'sstrel'
CASE (io_wavrel)
   namedesc = 'wavrel'
CASE (io_nstrel)
   namedesc = 'nstrel'//TRIM(cfil)

!
!2.4 Open boundaries
!-------------------
!

CASE (io_1uvsur)
   namedesc = '1uvsur'//TRIM(cfil)
CASE (io_2uvobc)
   namedesc = '2uvobc'//TRIM(cfil)
CASE (io_2xyobc)
   namedesc = '2xyobc'//TRIM(cfil)
CASE (io_3uvobc)
   namedesc = '3uvobc'//TRIM(cfil)
CASE (io_3xyobc)
   namedesc = '3xyobc'//TRIM(cfil)
CASE (io_salobc)
   namedesc = 'salobc'//TRIM(cfil)
CASE (io_tmpobc)
   namedesc = 'tmpobc'//TRIM(cfil)
CASE (io_sedobc)
   namedesc = 'sedobc'//TRIM(cfil)
CASE (io_bioobc)
   namedesc = 'bioobc'//TRIM(cfil)
CASE (io_rlxobc)
   namedesc = 'rlxobc'

!
!2.5 Nested output
!-----------------
!

CASE (io_nstspc)
   namedesc = 'nstspc'
CASE (io_2uvnst)
   namedesc = '2uvnst'//TRIM(cfil)
CASE (io_3uvnst)
   namedesc = '3uvnst'//TRIM(cfil)
CASE (io_2xynst)
   namedesc = '2xynst'//TRIM(cfil)
CASE (io_3xynst)
   namedesc = '3xynst'//TRIM(cfil)
CASE (io_salnst)
   namedesc = 'salnst'//TRIM(cfil)
CASE (io_tmpnst)
   namedesc = 'tmpnst'//TRIM(cfil)
CASE (io_sednst)
   namedesc = 'sednst'//TRIM(cfil)
CASE (io_bionst)
   namedesc = 'bionst'//TRIM(cfil)

!
!2.6 Surface data
!----------------
!

CASE (io_metsur)
   namedesc = 'metsur'
CASE (io_sstsur)
   namedesc = 'sstsur'
CASE (io_wavsur)
   namedesc = 'wavsur'
   
!
!2.7 Structures
!--------------
!
CASE (io_drycel)
   namedesc = 'drycel'
CASE (io_thndam)
   namedesc = 'thndam'
CASE (io_weibar)
   namedesc = 'weibar'

!
!2.8 Discharges
!--------------
!

CASE (io_disspc)
   namedesc = 'disspc'
CASE (io_disloc)
   namedesc = 'disloc'//TRIM(cfil)
CASE (io_disvol)
   namedesc = 'disvol'//TRIM(cfil)
CASE (io_discur)
   namedesc = 'discur'//TRIM(cfil)
CASE (io_dissal)
   namedesc = 'dissal'//TRIM(cfil)
CASE (io_distmp)
   namedesc = 'distmp'//TRIM(cfil)

!
!2.9 Specifiers
!--------------
!

CASE (io_sedspc)
   namedesc = 'sedspc'
CASE (io_darspc)
   namedesc = 'darspc'
CASE (io_biospc)
   namedesc = 'biospc'
CASE (io_parspc)
   namedesc = 'parspc'
   
!
!2.10 Biological files
!---------------------
!

CASE (io_bioabs)
   namedesc = 'bioabs'
CASE (io_biorel)
   namedesc = 'biorel'
CASE (io_biosur)
   namedesc = 'biosur'

!
!2.11 Particle model
!-------------------
!

CASE (io_parcld)
   namedesc = 'parcld'//TRIM(cfil)
CASE (io_pargrd)
   namedesc = 'pargrd'
CASE (io_parphs)
   namedesc = 'parphs'
END SELECT

!
!3. Store attributes
!-------------------
!

IF (TRIM(modfiles(iddesc,ifil,iotype)%filename).EQ.'') THEN
   suffix = file_suffix(modfiles(iddesc,ifil,iotype)%form)
   SELECT CASE (modfiles(iddesc,ifil,iotype)%status)
      CASE ('N')
         modfiles(iddesc,ifil,iotype)%filename = 'NONE'
      CASE ('R','W')
         modfiles(iddesc,ifil,iotype)%filename = TRIM(intitle)//'.'//&
                        & TRIM(namedesc)//'.'//TRIM(suffix)
   END SELECT
ENDIF

IF (loglev2.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, logexit
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, logexit
   ENDIF
ENDIF
pglev = pglev - 1


RETURN

END SUBROUTINE set_modfiles_name

!========================================================================

SUBROUTINE set_modvars_atts(iddesc,ifil,iotype,varatts,numvars,noprofsd,&
                          & numprofs,novars)
!************************************************************************
!
! *set_modvars_atts* Obtain the variable attributes of a model file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modvars_routines.f90  V2.11
!
! Description -
!
! Module calls - inquire_var, lim_dims, set_biovars_atts, set_partvars_atts,
!                set_sedvars_atts, varatts_init
!
!************************************************************************
!
USE datatypes
USE gridpars
USE modids
USE nestgrids
USE obconds
USE paralpars
USE physpars
USE relaxation
USE structures
USE switches
USE tide
USE biovars_routines, ONLY: set_biovars_atts
USE partvars_routines, ONLY: set_partvars_atts
USE sedvars_routines, ONLY: set_sedvars_atts
USE datatypes_init, ONLY: varatts_init
USE utility_routines, ONLY: lim_dims

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, iotype, numvars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(2:) :: noprofsd
INTEGER, INTENT(IN), OPTIONAL :: novars, numprofs
TYPE (VariableAtts), INTENT(INOUT), DIMENSION(numvars) :: varatts

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER File id
!*ifil*      INTEGER File number
!*iotype*    INTEGER I/O type of file
!             = 1 => input
!             = 2 => output
!*varatts*   DERIVED Variable attributes
!*numvars*   INTEGER Number of variables in data file
!*noprofsd*  INTEGER Number of open boundary profiles in each data file
!                    (only for obc-files and ifil=1)
!*numprofs*  INTEGER Number of profiles in data file
!                    (only for obc-files and ifil>1)
!*novars*    INTEGER Number of "state" variables for which open boundary
!                    conditions are applied
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: packing
CHARACTER (LEN=12) :: cfil
INTEGER :: igrd, iifil, inode, ivar, nhtype, nocoords, nofiles, n1dat, n2dat
CHARACTER (LEN=lenname), DIMENSION(numvars) :: f90_name
CHARACTER (LEN=lendesc), DIMENSION(numvars) :: long_name, standard_name
CHARACTER (LEN=lenunit), DIMENSION(numvars) :: units
INTEGER, DIMENSION(numvars) :: data_type, ivarid, nrank, numvarsid
INTEGER, DIMENSION(numvars,5) :: nshape
TYPE (FileParams) :: filepars


pglev = pglev + 1
procname(pglev) = 'set_modvars_atts'

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

CALL varatts_init(varatts)

filepars = modfiles(iddesc,ifil,iotype)

standard_name = ''
f90_name = ''
long_name = ''
units = ''
data_type = MERGE(real_type,rlong_type,filepars%floattype.EQ.'S')
ivarid = 0
nrank = -1
nshape = 1
numvarsid = -1
nofiles = maxdatafiles(iddesc,1)
nocoords = filepars%nocoords
packing = filepars%packing

!
!2. Define attributes
!--------------------
!

SELECT CASE (iddesc)

!
!2.1 Domain decomposition
!------------------------
!

CASE (io_mppmod)
   ivarid(1:4) = (/iarr_mg_nc1procs,iarr_mg_nc2procs,iarr_mg_nr1procs,&
                 & iarr_mg_nr2procs/)
   nrank = 2
   nshape(1:4,1) = nprocs
   nshape(1:4,2) = nomglevels
   data_type(1:4) = int_type

!
!2.2 Initial conditions
!----------------------
!

CASE (io_inicon,io_fincon)

   SELECT CASE (ifil)

      CASE (ics_phys)
         
!        ---2-D mode
         ivar = 2
         IF (iopt_mode_2D.GT.0) THEN
            ivarid(ivar:ivar+2) = (/iarr_udvel,iarr_vdvel,iarr_zeta/)
            nrank(ivar:ivar+2) = MERGE(1,2,packing)
            nshape(ivar,1) = MERGE(noseaatu,nc,packing)
            nshape(ivar+1,1) = MERGE(noseaatv,nc,packing)
            nshape(ivar+2,1) = MERGE(noseaatc,nc,packing)
            IF (.NOT.packing) nshape(ivar:ivar+2,2) = nr 
            ivar = ivar + 3
         ENDIF

!        ---3-D currents
         IF (iopt_mode_3D.GT.0) THEN
            ivarid(ivar:ivar+1) = (/iarr_uvel,iarr_vvel/)
            nrank(ivar:ivar+1) = MERGE(2,3,packing)
            nshape(ivar,1) = MERGE(noseaatu,nc,packing)
            nshape(ivar+1,1) = MERGE(noseaatv,nc,packing)
            nshape(ivar:ivar+1,2) = MERGE(nz,nr,packing)
            IF (.NOT.packing) nshape(ivar:ivar+1,3) = nz    
            ivar = ivar + 2
            IF (iopt_grid_nodim.EQ.3) THEN
               ivarid(ivar) = iarr_wvel
               nrank(ivar) = MERGE(2,3,packing)
               nshape(ivar,1) = MERGE(noseaatc,nc,packing)
               nshape(ivar,2) = MERGE(nz+1,nr,packing)
               IF (.NOT.packing) nshape(ivar,3) = nz+1
               ivar = ivar + 1
            ENDIF
         ENDIF

!        ---density arrays
         IF (iopt_temp.GT.0) THEN
            ivarid(ivar) = iarr_temp
            nrank(ivar) = MERGE(2,3,packing)
            nshape(ivar,1) = MERGE(noseaatc,nc,packing)
            nshape(ivar,2) = MERGE(nz,nr,packing)
            IF (.NOT.packing) nshape(ivar,3) = nz
            ivar = ivar + 1
         ENDIF
         IF (iopt_sal.GT.0) THEN
            ivarid(ivar) = iarr_sal
            nrank(ivar) = MERGE(2,3,packing)
            nshape(ivar,1) = MERGE(noseaatc,nc,packing)
            nshape(ivar,2) = MERGE(nz,nr,packing)
            IF (.NOT.packing) nshape(ivar,3) = nz
            ivar = ivar + 1
         ENDIF
         
!        ---turbulence arrays
         IF (iopt_vdif_coef.EQ.3) THEN
            ivarid(ivar) = iarr_tke
            nrank(ivar) = MERGE(2,3,packing)
            nshape(ivar,1) = MERGE(noseaatc,nc,packing)
            nshape(ivar,2) = MERGE(nz+1,nr,packing)
            IF (.NOT.packing) nshape(ivar,3) = nz+1
            ivar = ivar + 1
            IF (iopt_turb_ntrans.EQ.2) THEN
               IF (iopt_turb_param.EQ.1) THEN
                  ivarid(ivar) = iarr_zlmix
               ELSEIF (iopt_turb_param.EQ.2) THEN
                  ivarid(ivar) = iarr_dissip
               ENDIF
               IF (iopt_turb_param.GT.0) THEN
                  nrank(ivar) = MERGE(2,3,packing)
                  nshape(ivar,1) = MERGE(noseaatc,nc,packing)
                  nshape(ivar,2) = MERGE(nz+1,nr,packing)
                  IF (.NOT.packing) nshape(ivar,3) = nz+1
                  ivar = ivar + 1
               ENDIF
            ENDIF
         ENDIF

!        ---bottom stress arrays
         IF (iopt_bstres_form.EQ.2) THEN
            IF (iopt_bstres_drag.EQ.2) THEN
               ivarid(ivar) = iarr_bdragcoefatc
               nrank(ivar) = MERGE(1,2,packing)
               nshape(ivar,1) = MERGE(noseaatc,nc,packing)
               IF (.NOT.packing) nshape(ivar,2) = nr
               ivar = ivar + 1
            ELSEIF (iopt_bstres_drag.EQ.4) THEN
               ivarid(ivar) = iarr_zroughatc
               nrank(ivar) = MERGE(1,2,packing)
               nshape(ivar,1) = MERGE(noseaatc,nc,packing)
               IF (.NOT.packing) nshape(ivar,2) = nr
               ivar = ivar + 1
            ENDIF

         ENDIF
         
!        ---tidal arrays
         IF (nconobc.GT.0) THEN
            ivarid(ivar:ivar+1) = (/iarr_fnode_obc,iarr_phase_obc/)
            nrank(ivar:ivar+1) = 1
            nshape(ivar:ivar+1,1) = nconobc
            ivar = ivar + 2
         ENDIF
         IF (iopt_astro_tide.EQ.1) THEN
            ivarid(ivar:ivar+1) = (/iarr_fnode_astro,iarr_phase_astro/)
            nrank(ivar:ivar+1) = 2
            nshape(ivar:ivar+1,1) = nconastro
            ivar = ivar + 2
         ENDIF

!        ---energy losses at weirs
         IF (iopt_weibar.EQ.1) THEN
            IF (numwbaru.GT.0) THEN
               ivarid(ivar) = iarr_wbarelossu
               nrank(ivar) = 1
               nshape(ivar,1) = numwbaru
               ivar = ivar + 1
            ENDIF
            IF (numwbarv.GT.0) THEN
               ivarid(ivar) = iarr_wbarelossv
               nrank(ivar) = 1
               nshape(ivar,1) = numwbarv
               ivar = ivar + 1
            ENDIF
         ENDIF
         
!        ---open boundary arrays
         IF (iopt_obc_sal.EQ.1) THEN
            IF (nobu.GT.0) THEN 
               ivarid(ivar) = iarr_obcsalatu
               nrank(ivar) = 4
               nshape(ivar,1:4) = (/nobu,nz,3,1/)
               ivar = ivar + 1
            ENDIF
            IF (nobv.GT.0) THEN 
               ivarid(ivar) = iarr_obcsalatv
               nrank(ivar) = 4
               nshape(ivar,1:4) = (/nobv,nz,3,1/)
               ivar = ivar + 1
            ENDIF
         ENDIF
         IF (iopt_obc_temp.EQ.1) THEN
            IF (nobu.GT.0) THEN 
               ivarid(ivar) = iarr_obctmpatu
               nrank(ivar) = 4
               nshape(ivar,1:4) = (/nobu,nz,3,1/)
               ivar = ivar + 1
            ENDIF
            IF (nobv.GT.0) THEN 
               ivarid(ivar) = iarr_obctmpatv
               nrank(ivar) = 4
               nshape(ivar,1:4) = (/nobv,nz,3,1/)
               ivar = ivar + 1
            ENDIF
         ENDIF
         IF (iopt_obc_2D.EQ.1) THEN
            IF (nobu.GT.0) THEN
               ivarid(ivar) = iarr_obc2uvatu
               nrank(ivar) = 2
               nshape(ivar,1:2) = (/nobu,2/)
               ivar = ivar + 1
            ENDIF
            IF (nobv.GT.0) THEN
               ivarid(ivar) = iarr_obc2uvatv
               nrank(ivar) = 2
               nshape(ivar,:2) = (/nobv,2/)
               ivar = ivar + 1
            ENDIF
         ENDIF
         IF (iopt_obc_3D.EQ.1) THEN
            IF (nobu.GT.0) THEN
               ivarid(ivar) = iarr_obc3uvatu
               nrank(ivar) = 3
               nshape(ivar,1:3) = (/nobu,nz,2/)
               ivar = ivar + 1
            ENDIF
            IF (nobv.GT.0) THEN
               ivarid(ivar) = iarr_obc3uvatv
               nrank(ivar) = 3
               nshape(ivar,1:3) = (/nobv,nz,2/)
            ENDIF

         ENDIF

      CASE (ics_sed,ics_morph)
         CALL set_sedvars_atts(iddesc,ifil,iotype,ivarid,nrank,nshape,&
                             & data_type,numvarsid,numvars)

      CASE (ics_bio)
         CALL set_biovars_atts(iddesc,ifil,iotype,ivarid,nrank,nshape,&
                             & data_type,numvars)

      CASE (ics_part)
         CALL set_partvars_atts(iddesc,ifil,iotype,ivarid,nrank,nshape,&
                              & data_type,numvars)

   END SELECT

!
!2.3 Model grid and data locations
!---------------------------------
!
!2.3.1 Model
!-----------
!

CASE (io_modgrd)
   nhtype = iopt_grid_htype
   ivar = 1
   IF (nhtype.EQ.2) THEN
      ivarid(1:2) = (/iarr_gdelxglb,iarr_gdelyglb/)
      nrank(1:2) = 1
      nshape(1,1) = nc
      nshape(2,1) = nr
      ivar= ivar + 2
   ELSEIF (nhtype.EQ.3) THEN
      ivarid(1:2) = (/iarr_gxcoordglb,iarr_gycoordglb/)
      nrank(1:2) = 2
      nshape(1,1:2) = (/nc,nr/)
      nshape(2,1:2) = (/nc,nr/)
      ivar= ivar + 2
   ENDIF
   IF (iopt_grid_vtype.EQ.2) THEN
      ivarid(ivar) = iarr_gsigcoordatw
      nrank(ivar) = 1
      nshape(ivar,1) = nz+1
      ivar = ivar + 1
   ELSEIF (iopt_grid_vtype.EQ.3) THEN
      ivarid(ivar) = iarr_gscoordglb
      nrank(ivar) = 3
      nshape(ivar,1:3) = (/nc-1,nr-1,nz+1/)
      ivar = ivar + 1
   ENDIF
   ivarid(ivar) = iarr_depmeanglb
   nrank(ivar) = 2
   nshape(ivar,1:2) = (/nc-1,nr-1/)
   ivar = ivar + 1
   IF (nobu.GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_iobu,iarr_jobu/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = nobu
      data_type(ivar:ivar+1) = int_type
      ivar = ivar + 2
   ENDIF
   IF (nobv.GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_iobv,iarr_jobv/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = nobv
      data_type(ivar:ivar+1) = int_type
   ENDIF

!
!2.3.2 Surface grids
!-------------------
!
!---using absolute coordinates
CASE (io_metabs,io_sstabs,io_wavabs,io_bioabs)
   IF (iddesc.EQ.io_metabs) igrd = igrd_meteo 
   IF (iddesc.EQ.io_sstabs) igrd = igrd_sst  
   IF (iddesc.EQ.io_wavabs) igrd = igrd_waves  
   IF (iddesc.EQ.io_bioabs) igrd = igrd_bio
   n1dat = surfacegrids(igrd,ifil)%n1dat
   n2dat = surfacegrids(igrd,ifil)%n2dat
   ivarid(1:2) = (/iarr_xcoord,iarr_ycoord/)
   nrank(1:2) = 2
   nshape(1:2,1) = n1dat
   nshape(1:2,2) = n2dat

!---using relative coordinates   
CASE (io_metrel,io_sstrel,io_wavrel,io_biorel)
   ivarid(1:3) = (/iarr_icoordC,iarr_jcoordC,iarr_weightsC/)
   nrank(1:2) = 2
   nshape(1:2,1) = nc
   nshape(1:2,2) = nr
   data_type(1:2) = int_type
   nrank(3) = 4
   nshape(3,1:4) = (/2,2,nc,nr/)

!
!2.3.3 Nested grids
!------------------
!
!2.3.3.1 using absolute coordinates
!----------------------------------
!

CASE (io_nstabs)
   
   ivar = 1

!  ---C-nodes
   IF (nohnstglbc(ifil).GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_xcoordatc,iarr_ycoordatc/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = nohnstglbc(ifil) 
      ivar = ivar + 2
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_zcoordatc
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/nohnstglbc(ifil),novnst(ifil)/)
         ivar = ivar + 1
      ENDIF
   ENDIF
   
!  ---U-nodes
   IF (nohnstglbu(ifil).GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_xcoordatu,iarr_ycoordatu/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = nohnstglbu(ifil)
      ivar = ivar + 2
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_zcoordatu
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/nohnstglbu(ifil),novnst(ifil)/)
         ivar = ivar + 1
      ENDIF
   ENDIF

!  ---V-nodes
   IF (nohnstglbv(ifil).GT.0) THEN
      inode_23312: DO inode=1,2
         ivarid(ivar:ivar+1) = (/iarr_xcoordatv,iarr_ycoordatv/)
         nrank(ivar:ivar+1) = 1
         nshape(ivar:ivar+1,1) = nohnstglbv(ifil)
         ivar = ivar + 2
      ENDDO inode_23312
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_zcoordatv
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/nohnstglbv(ifil),novnst(ifil)/)
         ivar = ivar + 1
      ENDIF
   ENDIF

!  ---X-nodes
   IF (nohnstglbx(ifil).GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_xcoordatx,iarr_ycoordatx/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = nohnstglbx(ifil)
      ivar = ivar + 2
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_zcoordatv
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/nohnstglbx(ifil),novnst(ifil)/)
         ivar = ivar + 1
      ENDIF
   ENDIF

!  ---Y-nodes
   IF (nohnstglby(ifil).GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_xcoordaty,iarr_ycoordaty/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = nohnstglby(ifil)
      ivar = ivar + 2
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_zcoordatu
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/nohnstglby(ifil),novnst(ifil)/)
         ivar = ivar + 1
      ENDIF
   ENDIF

!
!2.3.3.2 Using relative coordinates
!----------------------------------
!

CASE (io_nstrel)

   ivar = 1

!  ---C-nodes
   IF (nohnstglbc(ifil).GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_icoordCC,iarr_jcoordCC/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = nohnstglbc(ifil)
      data_type(ivar:ivar+1) = int_type
      ivar = ivar + 2
      ivarid(ivar) = iarr_weightsCC
      nrank(ivar) = 3
      nshape(ivar,1:3) = (/2,2,nohnstglbc(ifil)/)
      ivar = ivar + 1
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_zcoordatc
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/nohnstglbc(ifil),novnst(ifil)/)
         ivar = ivar + 1
      ENDIF
   ENDIF

!  ---U-nodes
   IF (nohnstglbu(ifil).GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_icoordUU,iarr_jcoordUU/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = nohnstglbu(ifil)
      data_type(ivar:ivar+1) = int_type
      ivar = ivar + 2
      ivarid(ivar) = iarr_weightsUU
      nrank(ivar) = 3
      nshape(ivar,1:3) = (/2,2,nohnstglbu(ifil)/)
      ivar = ivar + 1
      ivarid(ivar:ivar+1) = (/iarr_icoordCU,iarr_jcoordCU/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = nohnstglbu(ifil)
      ivar = ivar + 2
      ivarid(ivar) = iarr_weightsCU
      nrank(ivar) = 3
      nshape(ivar,1:3) = (/2,2,nohnstglbu(ifil)/)
      ivar = ivar + 1
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_zcoordatu
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/nohnstglbu(ifil),novnst(ifil)/)
         ivar = ivar + 1
      ENDIF
   ENDIF

!  ---V-nodes
   IF (nohnstglbv(ifil).GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_icoordVV,iarr_jcoordVV/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = nohnstglbv(ifil)
      data_type(ivar:ivar+1) = int_type
      ivar = ivar + 2
      ivarid(ivar) = iarr_weightsVV
      nrank(ivar) = 3
      nshape(ivar,1:3) = (/2,2,nohnstglbv(ifil)/)
      ivar = ivar + 1
      ivarid(ivar:ivar+1) = (/iarr_icoordCV,iarr_jcoordCV/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = nohnstglbv(ifil)
      ivar = ivar + 2
      ivarid(ivar) = iarr_weightsCV
      nrank(ivar) = 3
      nshape(ivar,1:3) = (/2,2,nohnstglbv(ifil)/)
      ivar = ivar + 1
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_zcoordatv
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/nohnstglbv(ifil),novnst(ifil)/)
         ivar = ivar + 1
      ENDIF
   ENDIF

!  ---X-nodes
   IF (nohnstglbx(ifil).GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_icoordVX,iarr_jcoordVX/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = nohnstglbx(ifil)
      data_type(ivar:ivar+1) = int_type
      ivar = ivar + 2
      ivarid(ivar) = iarr_weightsVX
      nrank(ivar) = 3
      nshape(ivar,1:3) = (/2,2,nohnstglbx(ifil)/)
      ivar = ivar + 1
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_zcoordatx
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/nohnstglbx(ifil),novnst(ifil)/)
         ivar = ivar + 1
      ENDIF
   ENDIF

!  ---Y-nodes
   IF (nohnstglby(ifil).GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_icoordUY,iarr_jcoordUY/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = nohnstglby(ifil)
      data_type(ivar:ivar+1) = int_type
      ivar = ivar + 2
      ivarid(ivar) = iarr_weightsUY
      nrank(ivar) = 3
      nshape(ivar,1:3) = (/2,2,nohnstglby(ifil)/)
      ivar = ivar + 1
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_zcoordaty
         nshape(ivar,1:2) = (/nohnstglby(ifil),novnst(ifil)/)
         nrank(ivar) = 2
      ENDIF
   ENDIF


!2.4 Open (surface) boundaries
!-----------------------------
!
!2.4.1 1-D case
!--------------
!

CASE (io_1uvsur)
   IF (ifil.EQ.1) THEN
      ivar = 1
      ivarid(ivar) = iarr_isur1dtype
      nrank(ivar) = 0
      ivar = ivar + 1
      IF (isur1dtype.EQ.1.OR.isur1dtype.EQ.3) THEN
         ivarid(ivar:ivar+3) = (/iarr_gxslope_amp,iarr_gxslope_pha,&
                               & iarr_gyslope_amp,iarr_gyslope_pha/)
         nrank(ivar:ivar+3) = 1
         nshape(ivar:ivar+3,1) = nconobc
         ivar = ivar + 4
      ENDIF
      IF (isur1dtype.EQ.1.OR.isur1dtype.EQ.2) THEN
         ivarid(ivar:ivar+1) = (/iarr_zeta_amp,iarr_zeta_pha/)
         nrank(ivar:ivar+1) = 1
         nshape(ivar:ivar+1,1) = nconobc
      ENDIF
   ELSEIF (ifil.EQ.2) THEN
      SELECT CASE (isur1dtype)
         CASE (1); ivarid(2:4) = (/iarr_gxslope,iarr_gyslope,iarr_zetasur/)
         CASE (2); ivarid(2) = iarr_zetasur
         CASE (3); ivarid(2:3) = (/iarr_gxslope,iarr_gyslope/)
      END SELECT
      nrank(2:numvars) = 1
      nshape(2:numvars,1) = 1
   ENDIF

!
!2.4.2 2-D case (normal)
!-----------------------
!

CASE (io_2uvobc)
   IF (ifil.EQ.1) THEN
      ivar = 1
      IF (nobu.GT.0) THEN
         ivarid(ivar:ivar+1) = (/iarr_ityp2dobu,iarr_iloczobu/)
         nrank(ivar:ivar+1) = 1
         nshape(ivar:ivar+1,1) = nobu
         data_type(ivar:ivar+1) = int_type
         ivar = ivar + 2
      ENDIF
      IF (nobv.GT.0) THEN
         ivarid(ivar:ivar+1) = (/iarr_ityp2dobv,iarr_iloczobv/)
         nrank(ivar:ivar+1) = 1
         nshape(ivar:ivar+1,1) = nobv
         data_type(ivar:ivar+1) = int_type
         ivar = ivar + 2
      ENDIF
      
!     ---number of data and index maps per input file
      IF (nofiles.GT.1) THEN
         ivarid(ivar:ivar+2) = (/iarr_no2dobuv,iarr_iobc2dtype,&
                               & iarr_index2dobuv/)
         nrank(ivar:ivar+2) = (/1,1,2/)
         nshape(ivar:ivar+2,1) = (/nofiles-1,nofiles-1,nobu+nobv/)
         nshape(ivar+2,2) = nofiles-1
         data_type(ivar:ivar+2) = int_type
         ivar = ivar + 3
      ENDIF

!     ---amplitudes and phases
      IF (nobu.GT.0.AND.nconobc.GT.0) THEN
         ivarid(ivar:ivar+3) = (/iarr_udatobu_amp,iarr_zdatobu_amp,&
                               & iarr_udatobu_pha,iarr_zdatobu_pha/)
         nrank(ivar:ivar+3) = 2
         nshape(ivar:ivar+3,1) = nobu
         nshape(ivar:ivar+3,2) = nconobc
         ivar = ivar + 4
      ENDIF
      IF (nobv.GT.0.AND.nconobc.GT.0) THEN
         ivarid(ivar:ivar+3) = (/iarr_vdatobv_amp,iarr_zdatobv_amp,&
                               & iarr_vdatobv_pha,iarr_zdatobv_pha/)
         nrank(ivar:ivar+3) = 2
         nshape(ivar:ivar+3,1) = nobv 
         nshape(ivar:ivar+3,2) = nconobc
         ivar = ivar + 4
      ENDIF

!     ---discharges along open boundary sections
      IF (nqsecobu.GT.0) THEN
         ivarid(ivar) = iarr_iqsecobu
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/nqsecobu,2/)
         data_type(ivar) = int_type
      ENDIF
      IF (nqsecobv.GT.0) THEN
         ivarid(ivar) = iarr_jqsecobv
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/nqsecobv,2/)
         data_type(ivar) = int_type
      ENDIF

   ELSE
!     ---data file
      SELECT CASE (iobc2dtype(ifil))
         CASE (1)
            ivarid(2:3) = (/iarr_vel2dobc,iarr_zetaobc/)
         CASE (2)
            ivarid(2) = iarr_zetaobc
         CASE (3)
            ivarid(2) = iarr_vel2dobc
      END SELECT
      nrank(2:numvars) = 1
      nshape(2:numvars,1) = no2dobuv(ifil)
   ENDIF

!
!2.4.3 2-D case (tangential)
!---------------------------
!

CASE (io_2xyobc)
   IF (ifil.EQ.1) THEN
      ivar = 1
      IF (nobx.GT.0) THEN
         ivarid(ivar) = iarr_ityp2dobx
         nrank(ivar) = 1
         nshape(ivar,1) = nobx
         data_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF
      IF (noby.GT.0) THEN
         ivarid(ivar) = iarr_ityp2doby
         nrank(ivar) = 1
         nshape(ivar,1) = noby
         data_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF
      
!     ---number of data and index maps per input file
      IF (nofiles.GT.1) THEN
         ivarid(ivar:ivar+1) = (/iarr_no2dobxy,iarr_index2dobxy/)
         nrank(ivar:ivar+1) = (/1,2/)
         nshape(ivar:ivar+1,1) = (/nofiles-1,nobx+noby/)
         nshape(ivar+1,2) = nofiles-1
         data_type(ivar:ivar+1) = int_type
         ivar = ivar + 2
      ENDIF

!     ---amplitudes and phases
      IF (nobx.GT.0.AND.nconobc.GT.0) THEN
         ivarid(ivar:ivar+1) = (/iarr_vdatobx_amp,iarr_vdatobx_pha/)
         nrank(ivar:ivar+1) = 2
         nshape(ivar:ivar+1,1) = nobx
         nshape(ivar:ivar+1,2) = nconobc
         ivar = ivar + 2
      ENDIF
      IF (noby.GT.0.AND.nconobc.GT.0) THEN
         ivarid(ivar:ivar+1) = (/iarr_udatoby_amp,iarr_udatoby_pha/)
         nrank(ivar:ivar+1) = 2
         nshape(ivar:ivar+1,1) = noby 
         nshape(ivar:ivar+1,2) = nconobc
      ENDIF

   ELSE
!     ---data file 
      ivarid(2) = iarr_vel2dobc
      nrank(2:numvars) = 1
      nshape(2:numvars,1) = no2dobxy(ifil)

   ENDIF

!
!2.4.4 3-D case (normal)
!-----------------------
!

CASE (io_3uvobc,io_salobc,io_tmpobc,io_sedobc,io_bioobc)
   IF (ifil.EQ.1) THEN
      ivar = 0
      IF (nobu.GT.0) THEN
         ivar = ivar + 1
         ivarid(ivar) = iarr_itypobu
         nrank(ivar) = 1
         nshape(ivar,1) = nobu
         data_type(ivar) = int_type
         ivar = ivar + 1
         ivarid(ivar) = iarr_iprofobu
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/nobu,novars/)
         data_type(ivar) = int_type
      ENDIF
      IF (nobv.GT.0) THEN
         ivar = ivar + 1
         ivarid(ivar) = iarr_itypobv
         nrank(ivar) = 1
         nshape(ivar,1) = nobv
         data_type(ivar) = int_type
         ivar = ivar + 1
         ivarid(ivar) = iarr_iprofobv
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/nobv,novars/)
         data_type(ivar) = int_type
      ENDIF
      IF (nofiles.GT.1) THEN
         ivar = ivar + 1
         ivarid(ivar) = iarr_noprofsd
         nrank(ivar) = 1
         nshape(ivar,1) = nofiles-1
         data_type(ivar) = int_type
      ENDIF
      iifil_244: DO iifil=2,nofiles
         ivar = ivar + 1
         WRITE (cfil,'(I12)') iifil; cfil = ADJUSTL(cfil)
         f90_name(ivar) = 'indexprof'//TRIM(cfil)
         long_name(ivar) = 'Profile number mapping vector'
         data_type(ivar) = int_type
         nrank(ivar) = 1
         nshape(ivar,1) = noprofsd(iifil)
         ivar = ivar + 1
         f90_name(ivar) = 'indexvar'//TRIM(cfil)
         long_name(ivar) = 'Variable number mapping array'
         data_type(ivar) = int_type
         nrank(ivar) = 1
         nshape(ivar,1) = noprofsd(iifil)
         ivar = ivar + 1
      ENDDO iifil_244
      IF (norlxzones.GT.0) THEN
         ivarid(ivar) = iarr_iprofrlx
         nshape(ivar,1) = norlxzones
         data_type(ivar) = int_type
      ENDIF
   ELSE
      nrank(2) = 2
      nshape(2,1:2) = (/numprofs,nz/)
      IF (iddesc.EQ.io_3uvobc) THEN
         ivarid(2) = iarr_profvel
      ELSEIF (iddesc.EQ.io_salobc) THEN
         ivarid(2) = iarr_profsal
      ELSEIF (iddesc.EQ.io_tmpobc) THEN
         ivarid(2) = iarr_proftmp
      ELSEIF (iddesc.EQ.io_sedobc) THEN
         ivarid(2) = iarr_profsed
      ELSEIF (iddesc.EQ.io_bioobc) THEN
         ivarid(2) = iarr_profbio
      ENDIF
   ENDIF

!
!2.4.5 3-D case (tangential)
!---------------------------
!

CASE (io_3xyobc)
   IF (ifil.EQ.1) THEN
      ivar = 0
      IF (nobx.GT.0) THEN
         ivar = ivar + 1
         ivarid(ivar) = iarr_itypobx
         nrank(ivar) = 1
         nshape(ivar,1) = nobx
         data_type(ivar) = int_type
         ivar = ivar + 1
         ivarid(ivar) = iarr_iprofobx
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/nobx,novars/)
         data_type(ivar) = int_type
      ENDIF
      IF (noby.GT.0) THEN
         ivar = ivar + 1
         ivarid(ivar) = iarr_itypoby
         nrank(ivar) = 1
         nshape(ivar,1) = noby
         data_type(ivar) = int_type
         ivar = ivar + 1
         ivarid(ivar) = iarr_iprofoby
         nrank(ivar) = 2
         nshape(ivar,1:2) = (/noby,novars/)
         data_type(ivar) = int_type
      ENDIF
      IF (nofiles.GT.1) THEN
         ivar = ivar + 1
         ivarid(ivar) = iarr_noprofsd
         nrank(ivar) = 1
         nshape(ivar,1) = nofiles-1
         data_type(ivar) = int_type
      ENDIF
      iifil_245: DO iifil=2,nofiles
         ivar = ivar + 1
         WRITE (cfil,'(I12)') iifil; cfil = ADJUSTL(cfil)
         f90_name(ivar) = 'indexprof'//TRIM(cfil)
         long_name(ivar) = 'Profile number mapping vector'
         nrank(ivar) = 1
         nshape(ivar,1) = noprofsd(iifil)
         data_type(ivar) = int_type
         ivar = ivar + 1
         f90_name(ivar) = 'indexvar'//TRIM(cfil)
         long_name(ivar) = 'Variable number mapping array'
         nrank(ivar) = 1
         nshape(ivar,1) = noprofsd(iifil)
         data_type(ivar) = int_type
         ivar = ivar + 1
      ENDDO iifil_245
      IF (norlxzones.GT.0) THEN
         ivarid(ivar) = iarr_iprofrlx
         nshape(ivar,1) = norlxzones
         data_type(ivar) = int_type
      ENDIF
   ELSE
      nrank(2) = 2
      nshape(2,1:2) = (/numprofs,nz/)
      ivarid(2) = iarr_profvel
   ENDIF

!
!2.4.6 Relaxation conditions
!---------------------------
!

CASE (io_rlxobc)
   ivarid(1:7) = (/iarr_inodesrlx,iarr_idirrlx,iarr_ityprlx,iarr_iposrlx,&
                 & iarr_jposrlx,iarr_ncrlx,iarr_nrrlx/)
   nrank = 1
   nshape(1,1) = 3
   nshape(1:7,1) = norlxzones
   data_type(1:7) = int_type   


!
!2.5 Nested output
!-----------------
!

CASE (io_nstspc)
   ivarid(1:8) = (/iarr_nestcoords,iarr_nohnstglbc,iarr_nohnstglbu,&
                 & iarr_nohnstglbv,iarr_novnst,iarr_inst2dtype,&
                 & iarr_nohnstglbx,iarr_nohnstglby/)
   nrank(1:8) = 1
   nshape(1:8,1) = nonestsets
   data_type(1:8) = int_type
   ivar = 8
   IF (iopt_sed.GT.0) THEN
      ivarid(ivar+1:ivar+2) = (/iarr_nosednst,iarr_instsed/)
      nrank(ivar+1) = 1; nrank(ivar+2) = 2
      nshape(1,ivar+1) = nonestsets
      nshape(1:2,ivar+2) = (/maxsedvars,nonestsets/)
      data_type(ivar+1:ivar+2) = int_type
   ENDIF
   IF (iopt_biolgy.GT.0) THEN
      ivarid(ivar+1:ivar+2) = (/iarr_nobionst,iarr_instbio/)
      nrank(ivar+1) = 1; nrank(ivar+2) = 2
      nshape(1,ivar+1) = nonestsets
      nshape(1:2,ivar+2) = (/maxbiovars,nonestsets/)
      data_type(ivar+1:ivar+2) = int_type
   ENDIF

CASE (io_2uvnst)
   SELECT CASE (inst2dtype(ifil))
      CASE (1)
         ivarid(2:3) = (/iarr_vel2d,iarr_zeta/)
      CASE (2)
         ivarid(2) = iarr_zeta
      CASE (3)
         ivarid(2) = iarr_vel2d
   END SELECT
   nrank(2:numvars) = 1
   nshape(2:numvars,1) = SUM(nohnstuvprocs(:,ifil))

CASE (io_2xynst)
   ivarid(2) = iarr_vel2d
   nrank(2:numvars) = 1
   nshape(2:numvars,1) = SUM(nohnstxyprocs(:,ifil))

CASE (io_3uvnst)
   ivarid(2) = iarr_profvel
   nrank(2) = 2
   nshape(2,1) = SUM(nohnstuvprocs(:,ifil))
   nshape(2,2) = novnst(ifil)

CASE (io_3xynst)
   ivarid(2) = iarr_profvel
   nrank(2) = 2
   nshape(2,1) = SUM(nohnstxyprocs(:,ifil))
   nshape(2,2) = novnst(ifil)

CASE (io_salnst)
   ivarid(2) = iarr_profsal
   nrank(2) = 2
   nshape(2,1) = SUM(nohnstcprocs(:,ifil))
   nshape(2,2) = novnst(ifil)

CASE (io_tmpnst)
   ivarid(2) = iarr_proftmp
   nrank(2) = 2
   nshape(2,1) = SUM(nohnstcprocs(:,ifil))
   nshape(2,2) = novnst(ifil)

CASE (io_sednst)
   ivarid(2) = iarr_profsed
   nrank(2) = 2
   nshape(2,1) = SUM(nohnstcprocs(:,ifil))
   nshape(2,2) = novnst(ifil)

CASE (io_bionst)
   ivarid(2) = iarr_profbio
   nrank(2) = 2
   nshape(2,1) = SUM(nohnstcprocs(:,ifil))
   nshape(2,2) = novnst(ifil)

!
!2.6 Surface data
!----------------
!
!---meteo
CASE (io_metsur)
   nrank(2:numvars) = 2
   nshape(2:numvars,1) = surfacegrids(igrd_meteo,1)%n1dat
   nshape(2:numvars,2) = surfacegrids(igrd_meteo,1)%n2dat
   ivar = 1
   IF (iopt_meteo_stres.EQ.1) THEN
      IF (iopt_meteo_data.EQ.1) THEN
         ivarid(ivar+1:ivar+2) = (/iarr_uwindatc,iarr_vwindatc/)
         ivar = ivar + 2
      ELSEIF (iopt_meteo_data.EQ.2) THEN
         ivarid(ivar+1:ivar+2) = (/iarr_usstresatc,iarr_vsstresatc/)
         ivar = ivar + 2
      ENDIF
   ENDIF
   IF (iopt_meteo_pres.EQ.1) THEN
      ivarid(ivar+1) = iarr_atmpres
      ivar = ivar + 1
   ENDIF
   IF (iopt_meteo_heat.EQ.1) THEN
      SELECT CASE (iopt_meteo_data)
         CASE (1)
            ivarid(ivar+1:ivar+3) = (/iarr_airtemp,iarr_relhum,iarr_cloud_cover/)
            ivar = ivar + 3
         CASE (2)
            ivarid(ivar+1:ivar+2) = (/iarr_qnonsol,iarr_qrad/)
            ivar = ivar + 2
      END SELECT
   ENDIF
   SELECT CASE (iopt_meteo_precip)
      CASE (1)
         ivarid(ivar+1) = iarr_evapminprec
      CASE (2)
         ivarid(ivar+1) = iarr_precipitation
   END SELECT

!---SST
CASE (io_sstsur)
   ivarid(2) = iarr_sst
   nrank(2) = 2
   nshape(2,1) = surfacegrids(igrd_sst,1)%n1dat
   nshape(2,2) = surfacegrids(igrd_sst,1)%n2dat
   
!---waves
CASE (io_wavsur)
   ivar = 1
   ivarid(ivar+1:ivar+3) = (/iarr_waveheight,iarr_waveperiod,iarr_wavedir/)
   ivar = ivar + 3
   IF (iopt_waves_form.EQ.2) THEN
      ivarid(ivar+1:ivar+2) = (/iarr_wavevel,iarr_waveexcurs/)
      ivar = ivar + 2
      IF (iopt_waves_curr.EQ.1) THEN
         ivarid(ivar+1:ivar+2) = (/iarr_umstokesatc,iarr_vmstokesatc/)
         ivar = ivar + 2
         IF (iopt_waves_pres.EQ.1) THEN
            ivarid(ivar+1) = iarr_wavepres
            ivar = ivar + 1
         ENDIF
      ENDIF
   ENDIF
   IF (iopt_waves_dissip.EQ.1) THEN
      ivarid(ivar+1:ivar+4) = (/iarr_umswdissipatc,iarr_vmswdissipatc,&
                              & iarr_umbwdissipatc,iarr_vmbwdissipatc/)
   ENDIF
   nrank(2:numvars) = 2
   IF (iopt_waves_couple.EQ.0) THEN
      nshape(2:numvars,1) = surfacegrids(igrd_waves,1)%n1dat
      nshape(2:numvars,2) = surfacegrids(igrd_waves,1)%n2dat
   ELSE
      nshape(2:numvars,1) = nc
      nshape(2:numvars,2) = nr
   ENDIF

!
!2.7 Structures
!--------------
!

CASE (io_drycel)
   ivarid(1:2) = (/iarr_idry,iarr_jdry/)
   nrank(1:2) = 1
   nshape(1:2,1) = numdry
   data_type(1:2) = int_type

CASE (io_thndam)
   ivar = 1
   IF (numthinu.GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_ithinu,iarr_jthinu/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = numthinu
      data_type(ivar:ivar+1) = int_type
      ivar = ivar + 2
   ENDIF
   IF (numthinv.GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_ithinv,iarr_jthinv/)
      nrank(ivar:ivar+1) = 1
      nshape(ivar:ivar+1,1) = numthinv
      data_type(ivar:ivar+1) = int_type
   ENDIF

CASE (io_weibar)
   ivar = 1
   IF (numwbaru.GT.0) THEN
      ivarid(ivar:ivar+7) = (/iarr_iwbaru,iarr_jwbaru,iarr_oricoefu,&
                            & iarr_oriheightu,iarr_orisillu,iarr_wbarcoefu,&
                            & iarr_wbarcrestu,iarr_wbarmodlu/)
      nrank(ivar:ivar+7) = 1
      nshape(ivar:ivar+7,1) = numwbaru
      data_type(ivar:ivar+1) = int_type
      ivar = ivar + 8
   ENDIF
   IF (numwbarv.GT.0) THEN
      ivarid(ivar:ivar+7)=(/iarr_iwbarv,iarr_jwbarv,iarr_oricoefv,&
                          & iarr_oriheightv,iarr_orisillv,iarr_wbarcoefv,&
                          & iarr_wbarcrestv,iarr_wbarmodlv/)
      nrank(ivar:ivar+7) = 1
      nshape(ivar:ivar+7,1) = numwbarv
      data_type(ivar:ivar+1) = int_type
   ENDIF

!
!2.8 Discharges
!--------------
!

CASE (io_disspc)
   ivarid(1) = iarr_kdistype
   nrank(1) = 1
   nshape(1,1) = numdis
   data_type(1) = int_type

CASE (io_disloc)
   ivarid(2:4) = (/iarr_xdiscoord,iarr_ydiscoord,iarr_zdiscoord/)
   nrank(2:4) = 1
   nshape(2:4,1) = numdis

CASE (io_disvol)
   ivarid(2) = iarr_disvol
   nrank(2) = 1
   nshape(2,1) = numdis

CASE (io_discur)
   ivarid(2:3) = (/iarr_disarea,iarr_disdir/)
   nrank(2:3) = 1
   nshape(2:3,1) = numdis

CASE (io_dissal)
   ivarid(2) = iarr_dissal
   nrank(2) = 1
   nshape(2,1) = numdis

CASE (io_distmp)
   ivarid(2) = iarr_distmp
   nrank(2) = 1
   nshape(2,1) = numdis

!
!2.9 Sediment model files
!------------------------
!

CASE (io_sedspc,io_darspc)
   CALL set_sedvars_atts(iddesc,ifil,iotype,ivarid,nrank,nshape,data_type,&
                       & numvarsid(1:numvars),numvars)

!
!2.10 Biological model files
!---------------------------
!

CASE (io_biospc,io_biosur)
   CALL set_biovars_atts(iddesc,ifil,iotype,ivarid,nrank,nshape,data_type,&
                       & numvars)

!
!2.11 Particle model files
!------------------------
!

CASE (io_parcld,io_pargrd,io_parphs,io_parspc)
   CALL set_partvars_atts(iddesc,ifil,iotype,ivarid,nrank,nshape,data_type,&
                        & numvars)
   
END SELECT

!
!2.12 Time coordinate
!--------------------
!

IF (nocoords.EQ.1) THEN
   ivarid(1) = iarr_time
   nrank(1) = 2
   nshape(1,1) = lentime
   data_type(1) = char_type
   IF (filepars%form.EQ.'N') THEN
      nshape(1,2) = unlimited_NF90
   ELSE
      nshape(1,2) = -1
   ENDIF
ENDIF

!
!3. Store attributes
!-------------------
!

ivar_310: DO ivar=1,numvars
   varatts(ivar)%ivarid = ivarid(ivar)
   varatts(ivar)%nrank = nrank(ivar)
   varatts(ivar)%global_dims = nshape(ivar,:)
   varatts(ivar)%data_type = data_type(ivar)
   varatts(ivar)%numvar = numvarsid(ivar)
   IF (ivarid(ivar).EQ.0) THEN
      varatts(ivar)%f90_name = f90_name(ivar)
      varatts(ivar)%standard_name = standard_name(ivar)
      varatts(ivar)%long_name = long_name(ivar)
      varatts(ivar)%units = units(ivar)
   ENDIF
ENDDO ivar_310

IF (loglev2.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, logexit
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, logexit
   ENDIF
ENDIF
pglev = pglev - 1


RETURN

END SUBROUTINE set_modvars_atts


END MODULE modvars_routines
