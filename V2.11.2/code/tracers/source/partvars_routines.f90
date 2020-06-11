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

MODULE partvars_routines
!************************************************************************
!
! *partvars_routines* Utility routines for sediment model variables
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)partvars_routines.f90  V2.11
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description -
!
! Reference -
!
! Routines - inquire_partvar, set_partfiles_atts, set_partvars_atts
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS

!============================================================================

SUBROUTINE inquire_partvar(varid,f90_name,standard_name,comment,long_name,&
                         & units,node,vector_name,data_type,fill_value,nrank,&
                         & global_dims,local_dims,halo_dims,varatts)
!************************************************************************
!
! *inquire_partvar* Obtain information about a specific sediment model variable
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)partvars_routines.f90  V2.11
!
! Description -
!
! Calling program - inquire_var
!
! Module calls -
!
!************************************************************************
!
USE datatypes
USE iopars
USE gridpars
USE partids
USE partpars
USE switches

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
CHARACTER (LEN=lenunit) :: tunits, unit
CHARACTER (LEN=lennode) :: cnode
INTEGER :: dtype, nodim
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
fvalue = float_fill
commt = ''

SELECT CASE (ptime_unit)
   CASE (1); tunits = 'sec'
   CASE (2); tunits = 'minute'
   CASE (3); tunits = 'hour'
   CASE (4); tunits = 'day'
END SELECT

!
!2. Get particle model array info
!--------------------------------
!

SELECT CASE (varid)

!
!2.1 Volume distributions
!------------------------
!
   
CASE (iarr_part_ageconc)
   sname = 'contaminant_age_concentration'
   fname = 'ageconc'
   lname = 'Contaminant age concentration'
   unit = TRIM(tunits)
   nodim = 4
   ngdims(1:4) = (/nc,nr,nz,nolabels/)
   
CASE (iarr_part_ageconc_tot)
   sname = 'contaminant_age_concentration'
   fname = 'ageconc_tot'
   lname = 'Contaminant age concentration'
   unit = TRIM(tunits)
   nodim = 3
   ngdims(1:3) = (/nc,nr,nz/)
      
CASE (iarr_part_conc)
   sname = 'concentration'
   fname = 'conc'
   lname = 'Contaminant concentration'
   unit = 'kg m-3'
   nodim = 4
   ngdims(1:4) = (/nc,nr,nz,nolabels/)
      
CASE (iarr_part_conc_tot)
   sname = 'concentration'
   fname = 'conc_tot'
   lname = 'Contaminant concentration'
   unit = 'kg m-3'
   nodim = 3
   ngdims(1:3) = (/nc,nr,nz/)

!
!2.2 Model grid arrays
!---------------------
!   

CASE (iarr_maskatc)
   sname = 'Local dry_mask_at_c_nodes'
   fname = 'maskatc'
   lname = 'Local dry mask at C-nodes'
   unit = 'log'
   dtype = log_type
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
   nhdims = 0

CASE (iarr_maskatu)
   sname = 'Local dry_mask_at_u_nodes'
   fname = 'maskatu'
   lname = 'Local dry mask at U-nodes'
   unit = 'log'
   dtype = log_type
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
   nhdims = 0

CASE (iarr_maskatv)
   sname = 'Local dry_mask_at_v_nodes'
   fname = 'maskatv'
   lname = 'Local dry mask at V-nodes'
   unit = 'log'
   dtype = log_type
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
   nhdims = 0
   
CASE (iarr_vdifcoefpart)
   sname = 'particle_turbulent_diffusivity'
   fname = 'vdifcoefpart'
   lname = 'Particle vertical diffusion coefficient'
   unit = 'm2 s-1'
   nodim = 3
   ngdims(1:3) = (/nc,nr,nz+1/)

!
!2.2 Particulate matter
!---------------------
!

CASE (iarr_diamconc)
   sname = 'particle_diameter'
   fname = 'diamconc'
   lname = 'Particle diameter'
   unit = 'm'
   cnode  = ''
   nodim = 1
   ngdims(1) = nolabels

CASE (iarr_massconc)
   sname = 'particle_mass'
   fname = 'massconc'
   lname = 'Particle mass'
   unit = 'kg'
   cnode  = ''
   nodim = 1
   ngdims(1) = nolabels

CASE (iarr_rhoconc)
   sname = 'particle_density'
   fname = 'rhoconc'
   lname = 'Particle mass density'
   unit = 'kg m-3'
   cnode  = ''
   nodim = 1
   ngdims(1) = nolabels
   
CASE (iarr_volumeconc)
   sname = 'particle_volume'
   fname = 'volumeconc'
   lname = 'Particle volume'
   unit = 'm-3'
   cnode  = ''
   nodim = 1
   ngdims(1) = nolabels

!
!2.3 Particulate attributes
!--------------------------
!

CASE (iarr_part_age)
   sname = 'particle_age'
   fname = 'part_age'
   lname = 'Particle age'
   unit = TRIM(tunits)
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart

CASE (iarr_part_center)
   sname = 'particle_at_cloud_center'
   fname = 'part_center'
   lname = 'Particle at cloud center'
   unit = 'log'
   dtype = log_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart

CASE (iarr_part_displx)
   sname = 'x_particle_displacement'
   fname = 'part_displx'
   lname = 'Particle displacement in X-direction'
   unit = 'm'
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart

CASE (iarr_part_disply)
   sname = 'y_particle_displacement'
   fname = 'part_disply'
   lname = 'Particle displacement in Y-direction'
   unit = 'm'
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart

CASE (iarr_part_displz)
   sname = 'z_particle_displacement'
   fname = 'part_displz'
   lname = 'Particle displacement in vertical direction'
   unit = 'm'
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart

CASE (iarr_part_drift_state)
   sname = 'particle_drift_state'
   fname = 'part_drift_state'
   lname = 'Particle drift state'
   dtype = int_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart

CASE (iarr_part_icoord)
   sname = 'i_interpolation_index_coordinate_to_particle location'
   fname = 'part_icoord'
   lname = 'I interpolation index coordinate to particle location'
   dtype = int_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart

CASE (iarr_part_jcoord)
   sname = 'j_interpolation_index_coordinate_to_particle location'
   fname = 'part_jcoord'
   lname = 'J interpolation index coordinate to particle location'
   dtype = int_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart

CASE (iarr_part_kcoord)
   sname = 'k_interpolation_index_coordinate_to_particle location'
   fname = 'part_kcoord'
   lname = 'Vertical interpolation index coordinate to particle location'
   dtype = int_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart

CASE (iarr_part_kdistype)
   sname = 'type_of_vertical_initial_location'
   fname = 'part_kdistype'
   lname = 'Type of vertical initial location'
   dtype = int_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart
   
CASE (iarr_part_label)
   sname = 'particle_label'
   fname = 'part_label'
   lname = 'Particle label'
   dtype = int_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart
   
CASE (iarr_part_ntstart)
   sname = 'particle_release_time_index'
   fname = 'part_ntstart'
   lname = 'Particle release time index'
   dtype = int_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart
   
CASE (iarr_part_state)
   sname = 'type_of_particle_drift'
   fname = 'part_state'
   lname = 'Type of particle drift'
   dtype = int_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart
   
CASE (iarr_part_tstart)
   sname = 'particle_release_time'
   fname = 'part_tstart'
   lname = 'Particle release time'
   unit = TRIM(tunits)
   dtype = rlong_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart
   
CASE (iarr_part_xcoord)
   sname = 'local_x_coordinate_at_particle_location'
   fname = 'part_xcoord'
   lname = 'Local X-coordinate at particle location'
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart
   
CASE (iarr_part_xpos)
   sname =  MERGE('projection_x_coordinate','longitude              ',&
                & iopt_grid_sph.EQ.0)
   fname = 'part_xpos'
   lname =  MERGE('X-coordinate','longitude   ',iopt_grid_sph.EQ.0)
   lname = TRIM(lname)//' at particle location'
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   dtype = rlong_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart
   
CASE (iarr_part_ycoord)
   sname = 'local_y_coordinate_at_particle_location'
   fname = 'part_ycoord'
   lname = 'Local Y-coordinate at particle location'
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart
   
CASE (iarr_part_ypos)
   sname =  MERGE('projection_y_coordinate','longitude              ',&
                & iopt_grid_sph.EQ.0)
   fname = 'part_ypos'
   lname =  MERGE('Y-coordinate','longitude   ',iopt_grid_sph.EQ.0)
   lname = TRIM(lname)//' at particle location'
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   dtype = rlong_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart
   
CASE (iarr_part_zpos)
   sname =  'ocean_z_coordinate_at_particle_location'
   fname = 'part_zpos'
   lname =  'Z-coordinate at particle location'
   unit =  'm'
   dtype = rlong_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = nopart

!
!2.4 Particulate cloud attributes
!--------------------------------
!

CASE (iarr_cloud_data)
   sname = 'cloud_center_locations'
   fname = 'cloud_data'
   lname = 'Locations of clod centers'
   dtype = rlong_type
   cnode = ''
   nodim = 1
   ngdims(1) = 3
   
CASE (iarr_cloud_kdistype)
   sname = 'cloud_type_of_vertical_initial_location'
   fname = 'cloud_kdistype'
   lname = 'Type of initial cloud vertical location'
   dtype = int_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = noclouds

CASE (iarr_cloud_label)
   sname = 'cloud_label'
   fname = 'cloud_label'
   lname = 'Cloud label'
   dtype = int_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = noclouds

CASE (iarr_cloud_length)
   sname = 'cloud_ellipse_length_along_y_axis'
   fname = 'cloud_length'
   lname = 'Length of cloud ellipse length in Y-direction'
   unit = 'm'
   cnode = 'P'
   nodim = 1
   ngdims(1) = noclouds

CASE (iarr_cloud_nopart)
   sname = 'cloud_number_of_particles_in_cloud_at_each_release'
   fname = 'cloud_nopart'
   lname = 'Cloud number of particles in cloud at each release'
   dtype = int_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = noclouds

CASE (iarr_cloud_noreleases)
   sname = 'number_of_releases_per_particle_cloud'
   fname = 'cloud_noreleases'
   lname = 'Number of releases per particle cloud'
   dtype = int_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = noclouds

CASE (iarr_cloud_orientation)
   sname = 'cloud_ellipse_orientation'
   fname = 'cloud_orientation'
   lname = 'Cloud ellipse orientation'
   unit = 'radian'
   cnode = 'P'
   nodim = 1
   ngdims(1) = noclouds

CASE (iarr_cloud_state)
   sname = 'cloud_particle_drift_type'
   fname = 'cloud_state'
   lname = 'Cloud particle drift state'
   dtype = int_type
   cnode = 'P'
   nodim = 1
   ngdims(1) = noclouds

CASE (iarr_cloud_thick)
   sname = 'cloud_ellipse_vertical_thickness'
   fname = 'cloud_thick'
   lname = 'Vertical thickness of cloud ellipse' 
   unit = 'm'
   cnode = 'P'
   nodim = 1
   ngdims(1) = noclouds

CASE (iarr_cloud_width)
   sname = 'cloud_ellipse_length_along_x_axis'
   fname = 'cloud_width'
   lname = 'Length of cloud ellipse in X-direction'
   unit = 'm'
   cnode = 'P'
   nodim = 1
   ngdims(1) = noclouds
   
END SELECT

!
!3. Store
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

END SUBROUTINE inquire_partvar

!============================================================================

SUBROUTINE set_partfiles_atts(iddesc,ifil,iotype)
!************************************************************************
!
! *set_sedfiles_atts* Obtain global attributes of a particle model file
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)partvars_routines.f90  V2.11
!
! Description -
!
! Calling program - set_modfiles_atts
!
! Module calls -
!
!************************************************************************
!
USE iopars  
USE partswitches
USE switches
USE syspars

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
procname(pglev) = 'set_partfiles_atts'
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
      title = 'Particle initial conditions'
      novars = 7
      IF (iopt_grid_nodim.EQ.3) novars = novars + 2
      nocoords = 1

!  ---specifier arrays
   CASE (io_parspc)
      title = 'Particle specifier arrays'
      novars = 0
      IF (iopt_part_conc.GT.1.OR.iopt_part_dens.GT.0) novars = novars + 2
      IF (iopt_part_cloud.EQ.1) THEN
         IF (iopt_grid_nodim.NE.2) novars = novars + 1
         novars = novars + 8
      ENDIF
      nocoords = 0

!  ---cloud data
   CASE (io_parcld)
      title = 'Particle cloud locations'
      novars = 1
      nocoords = 1
      
!  ---model grid
   CASE (io_pargrd)
      title = 'Particle grid data'
      novars = 12
      IF (iopt_grid_nodim.EQ.3) THEN
         novars = MERGE(novars+1,novars+4,iopt_grid_vtype.LT.3)
      ENDIF
      nocoords = 0

!  ---physical data for the particle model
   CASE (io_parphs)
      title = 'Particle physical data'
      novars = 5
      IF (iopt_grid_nodim.EQ.3) novars = novars + 1
      IF (iopt_part_vdif.EQ.2) novars = novars + 1
      IF (iopt_part_dens.EQ.2) novars = novars + 1
      IF (iopt_fld.EQ.2) novars = novars + 3
      nocoords = 1

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

END SUBROUTINE set_partfiles_atts

!========================================================================

SUBROUTINE set_partvars_atts(iddesc,ifil,iotype,ivarid,nrank,nshape,data_type,&
                           & numvars)
!************************************************************************
!
! *set_partvars_atts* Obtain variable attributes for a particle model file
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)partvars_routines.f90  V2.11
!
! Description -
!
! Calling progam - set_modvars_atts
!
! Module calls -
!
!************************************************************************
!
USE gridpars
USE iopars
USE modids
USE partids
USE partpars
USE partswitches
USE switches
USE syspars

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, iotype, numvars
INTEGER, INTENT(OUT), DIMENSION(numvars) :: data_type, ivarid, nrank
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
!*numvars*    INTEGER Number of variables (coordinate and data) in the forcing
!                    file
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar
TYPE (FileParams) :: filepars


pglev = pglev + 1
procname(pglev) = 'set_partvars_atts'
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

ivarid = 0; nrank = -1; nshape = 1
filepars = modfiles(iddesc,ifil,iotype)
data_type = MERGE(real_type,rlong_type,filepars%floattype.EQ.'S')

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
   ivar = 1
   ivarid(ivar+1:ivar+2) = (/iarr_part_state,iarr_part_drift_state/)
   data_type(ivar+1:ivar+2) = int_type
   ivar = ivar + 2
   IF (iopt_grid_nodim.NE.2) THEN
      ivarid(ivar+1) = iarr_part_kdistype
      data_type(ivar+1) = int_type
      ivar = ivar + 1
   ENDIF
   ivarid(ivar+1:ivar+3) = (/iarr_part_label,iarr_part_xpos,iarr_part_ypos/)
   data_type(ivar+1) = int_type
   data_type(ivar+2:ivar+3) = rlong_type
   ivar = ivar + 3
   IF (iopt_grid_nodim.NE.2) THEN
      ivarid(ivar+1) = iarr_part_zpos
      data_type(ivar+1) = rlong_type
      ivar = ivar + 1
   ENDIF
   ivarid(ivar+1:ivar+2) = (/iarr_part_tstart,iarr_part_age/)
   data_type(ivar+1:ivar+2) = rlong_type 
   ivar = ivar + 2
   nrank(2:ivar) = 1
   nshape(2:ivar,1) = nopart

!
!2.2 Specifier arrays
!--------------------
!

CASE (io_parspc)

   ivar = 0
   IF (iopt_part_conc.GT.1.OR.iopt_part_dens.GT.0) THEN
      ivarid(ivar+1:ivar+2) = (/iarr_diamconc,iarr_rhoconc/)
      nrank(ivar+1:ivar+2) = 1
      nshape(ivar+1:ivar+2,1) = nolabels
      ivar = ivar + 2
   ENDIF
   IF (iopt_part_cloud.EQ.1) THEN
      IF (iopt_grid_nodim.NE.2) THEN
         ivarid(ivar+1) = iarr_cloud_kdistype
         ivar = ivar + 1
      ENDIF
      ivarid(ivar+1:ivar+8) = (/iarr_cloud_label,iarr_cloud_length,&
                              & iarr_cloud_nopart,iarr_cloud_noreleases,&
                              & iarr_cloud_orientation,iarr_cloud_state,&
                              & iarr_cloud_thick,iarr_cloud_width/)
      nrank(ivar+1:ivar+8) = 1
      nshape(ivar+1:ivar+8,1) = noclouds
      ivar = ivar + 8
   ENDIF

!
!2.3 Particle grid
!-----------------
!

CASE (io_pargrd)

!  ---coordinate arrays   
   ivar = 0
   ivarid(ivar+1:ivar+2) = (/iarr_gxcoord,iarr_gycoord/)
   nrank(ivar+1:ivar+2) = 2
   nshape(ivar+1:ivar+2,1) = nc
   nshape(ivar+1:ivar+2,2) = nr
   ivar = ivar + 2
   IF (iopt_grid_nodim.EQ.3) THEN
      IF (iopt_grid_vtype.LT.3) THEN
         ivarid(ivar+1) = iarr_gsigcoordatw
         nrank(ivar+1) = 1
         nshape(ivar+1,1) = nz+1
         ivar = ivar + 1
      ELSE
         ivarid(ivar+1:ivar+4) = (/iarr_gscoordatc,iarr_gscoordatu,&
                                 & iarr_gscoordatv,iarr_gscoordatw/)
         nrank(ivar+1:ivar+4)= 3
         nshape(ivar+1:ivar+4,1) = nc
         nshape(ivar+1:ivar+4,2) = nr
         nshape(ivar+1:ivar+3,3) = nz
         nshape(ivar+4,3) = nz + 1
         ivar = ivar + 4
      ENDIF
   ENDIF

!  ---grid spacings
   ivarid(ivar+1:ivar+2) = (/iarr_delxatc,iarr_delyatc/)
   nrank(ivar+1:ivar+2) = 2
   nshape(ivar+1:ivar+2,1) = nc
   nshape(ivar+1:ivar+2,2) = nr
   ivar = ivar + 2
   
!  ---bathymetry
   ivarid(ivar+1:ivar+3) = (/iarr_depmeanatc,iarr_depmeanatu,iarr_depmeanatv/)
   nrank(ivar+1:ivar+3) = 2
   nshape(ivar+1:ivar+3,1) = nc
   nshape(ivar+1:ivar+3,2) = nr
   ivar = ivar + 3

!  ---land mask
   ivarid(ivar+1:ivar+3) = (/iarr_maskatc,iarr_maskatu,iarr_maskatv/)
   nrank(ivar+1:ivar+3) = 2
   nshape(ivar+1:ivar+3,1) = nc
   nshape(ivar+1:ivar+3,2) = nr
   data_type(ivar+1:ivar+3) = log_type
   ivar = ivar + 3

!  ---model parameters   
   ivarid(ivar+1:ivar+2) = (/iarr_gacc_mean,iarr_density_ref/)
   nrank(ivar+1:ivar+2) = 0
   
!
!2.4 Physical data for the particle model
!----------------------------------------
!
   
CASE (io_parphs)
   ivar = 1
   ivarid(ivar+1:ivar+5) = (/iarr_deptotatc,iarr_deptotatu,iarr_deptotatv,&
                           & iarr_uvel,iarr_vvel/)
   nrank(ivar+1:ivar+3) = 2
   nrank(ivar+4:ivar+5) = 3
   nshape(ivar+1:ivar+5,1) = nc
   nshape(ivar+1:ivar+5,2) = nr
   nshape(ivar+4:ivar+5,3) = nz
   ivar = ivar + 5
   IF (iopt_grid_nodim.EQ.3) THEN
      ivarid(ivar+1) = iarr_wvel
      nrank(ivar+1) = 3
      nshape(ivar+1,1:3) = (/nc,nr,nz+1/)
      ivar = ivar + 1
   ENDIF
   IF (iopt_part_vdif.EQ.2) THEN
      ivarid(ivar+1) = iarr_vdifcoefpart
      nrank(ivar+1) = 3
      nshape(ivar+1,1:3) = (/nc,nr,nz+1/)
      ivar = ivar + 1
   ENDIF
   IF (iopt_part_dens.EQ.2) THEN
      ivarid(ivar+1) = iarr_dens
      nrank(ivar+1) = 1
      nshape(ivar+1,1:3) = (/nc,nr,nz/)
   ENDIF
   IF (iopt_fld.EQ.2) THEN
      ivarid(ivar+1:ivar+3) = (/iarr_maskatc,iarr_maskatu,iarr_maskatv/)
      data_type(ivar+1:ivar+3) = log_type
      nrank(ivar+1:ivar+3) = 1
      nshape(ivar+1:ivar+3,1) = nc
      nshape(ivar+1:ivar+3,2) = nr
   ENDIF

!
!2.5 Particle cloud locations
!----------------------------
!   

CASE (io_parcld)
   ivarid(2) = iarr_cloud_data
   data_type(2) = rlong_type 
   nrank(2) = 1
   nshape(2,1) = 3

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

END SUBROUTINE set_partvars_atts

END MODULE partvars_routines
