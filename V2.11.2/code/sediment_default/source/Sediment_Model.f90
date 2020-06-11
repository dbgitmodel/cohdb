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
! *Sediment_Model* external sediment routines called within the main
!                  COHERENS code
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - routines are empty but need to be declared 
!
! Routines - allocate_dar_arrays, allocate_morph_arrays, allocate_sed_arrays,
!    assign_cif_vars_sed, baroclinic_gradient_sed_cubic,
!    baroclinic_gradient_sed_sigma, baroclinic_gradient_sed_z,
!    buoyancy_frequency_sed, dar_equation, dar_init, deallocate_dar_arrays,
!    deallocate_morph_arrays, deallocate_sed_arrays, equation_of_state_sed,
!    exchange_sedics, initialise_sediment_arrays, morphology_accel,
!    morphology_equation, read_dar_spec, read_morphics, read_sedics,
!    read_sed_spec, sediment_equation, write_cif_vars_sed, write_dar_spec,
!    write_morphics, write_sedics, write_sed_spec
!
!************************************************************************
!

!========================================================================

SUBROUTINE allocate_dar_arrays
!************************************************************************
!
! *allocate_dar_arrays* Allocate arrays for the dredging/relocation models
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE allocate_dar_arrays

!========================================================================

SUBROUTINE allocate_morph_arrays
!************************************************************************
!
! *allocate_morphology_arrays* Allocate arrays for the morphological model
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE allocate_morph_arrays

!========================================================================

SUBROUTINE allocate_sed_arrays
!************************************************************************
!
! *allocate_sed_arrays* Allocate arrays for the sediment transport model
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE allocate_sed_arrays

!========================================================================

SUBROUTINE assign_cif_vars_dar(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_sed* convert the string data from a sediment CIF input
!                       line to the appropriate numeric or non-numeric format 
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - read_cif_mod_params
!
!************************************************************************
!
USE syspars

IMPLICIT NONE
!
!*Arguments
!
CHARACTER (LEN=lenname), INTENT(IN), OPTIONAL :: cname
CHARACTER (LEN=lencifvar), INTENT(IN), DIMENSION(MaxCIFvars) :: cvals
INTEGER, INTENT(IN) :: numvars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cname*     CHAR    Name of the variable(s)
!*cvals*     CHAR    Data values on input data line
!*numvars*   INTEGER Number of data variables
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE assign_cif_vars_dar

!============================================================================

SUBROUTINE assign_cif_vars_morph(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_morph* convert the string data from a morphology CIF input
!                         line to the appropriate numeric or non-numeric format
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - read_cif_mod_params
!
!************************************************************************
!
USE syspars

IMPLICIT NONE
!
!*Arguments
!
CHARACTER (LEN=lenname), INTENT(IN), OPTIONAL :: cname
CHARACTER (LEN=lencifvar), INTENT(IN), DIMENSION(MaxCIFvars) :: cvals
INTEGER, INTENT(IN) :: numvars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cname*     CHAR    Name of the variable(s)
!*cvals*     CHAR    Data values on input data line
!*numvars*   INTEGER Number of data variables
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE assign_cif_vars_morph

!============================================================================

SUBROUTINE assign_cif_vars_sed(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_sed* convert the string data from a sediment CIF input
!                       line to the appropriate numeric or non-numeric format 
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - read_cif_mod_params
!
!************************************************************************
!
USE syspars

IMPLICIT NONE
!
!*Arguments
!
CHARACTER (LEN=lenname), INTENT(IN), OPTIONAL :: cname
CHARACTER (LEN=lencifvar), INTENT(IN), DIMENSION(MaxCIFvars) :: cvals
INTEGER, INTENT(IN) :: numvars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cname*     CHAR    Name of the variable(s)
!*cvals*     CHAR    Data values on input data line
!*numvars*   INTEGER Number of data variables
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE assign_cif_vars_sed

!========================================================================

SUBROUTINE baroclinic_gradient_sed_cubic(zcoord,dzx,dzy,dzz,cdir)
!************************************************************************
!
! *baroclinic_gradient_sed_cubic* Sediment contribution for the baroclinic
!                                 density gradient (cube-H method)
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
! Calling program - baroclinic_gradient
!
!************************************************************************
!
USE gridpars

IMPLICIT NONE
!
!*Arguments
!
CHARACTER (LEN=1) :: cdir
REAL, DIMENSION(0:ncloc+1,0:nrloc+1,nz), INTENT(IN) ::  dzx, dzy, dzz
REAL,DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz), INTENT(IN) :: zcoord

!
! Name      Type    Purpose
!-------------------------------------------------------------------------------
!*zcoord*   REAL    Z cordinate at W points                                  [m]
!*dzx*      REAL    Harmonic derivative of x coordinate
!*dzy*      REAL    Harmonic derivative of y coordinate
!*dzz*      REAL    Harmonic derivative of z coordinate
!*cdir*     CHAR    Direction of the derivative ('X','Y')
!
!*******************************************************************************
!


RETURN

END SUBROUTINE baroclinic_gradient_sed_cubic

!========================================================================

SUBROUTINE baroclinic_gradient_sed_sigma(zcoord,cdir)
!************************************************************************
!
! *baroclinic_grad_sed_sigma* Baroclinic density gradient due to sediments
!                             (sigma second order method)
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
! Calling program - baroclinic_gradient
!
!************************************************************************
!
USE gridpars

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=1) :: cdir
REAL,DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz), INTENT(IN) :: zcoord

!
! Name      Type    Purpose
!-------------------------------------------------------------------------------
!*zcoord*   REAL    Z-coordinates at W-nodes                                 [m]
!*cdir*     CHAR    Direction of the derivative ('X','Y')
!
!*******************************************************************************
!


RETURN

END SUBROUTINE baroclinic_gradient_sed_sigma

!========================================================================

SUBROUTINE baroclinic_gradient_sed_z(sigint,kint,cdir)
!************************************************************************
!
! *baroclinic_grad_sed_z* Baroclinic density gradient due to sediments
!                         (Z-level method)
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
! Calling program - baroclinic_gradient
!
!************************************************************************
!
USE gridpars

IMPLICIT NONE
!
!*Arguments
!
CHARACTER (LEN=1) :: cdir
INTEGER, DIMENSION(ncloc,nrloc,2:nz+1,2), INTENT(IN) :: kint
REAL, DIMENSION(ncloc,nrloc,2:nz+1,2), INTENT(IN) :: sigint

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*sigint*   REAL    Interpolated sigma coordinates
!*kint*     INTEGER Counter for the interpolated coordinates
!*cdir*     CHAR    Direction of the derivative ('X','Y')
!
!******************************************************************************
!


RETURN

END SUBROUTINE baroclinic_gradient_sed_z

!========================================================================

SUBROUTINE buoyancy_frequency_sed(bgrad)
!************************************************************************
!
! *buoyancy_frequency_sed* Sediment contribution for the squared buoyancy
!                          frequency
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
! Calling program - buoyancy_frequency
!
!************************************************************************
!
USE gridpars

IMPLICIT NONE

!
!*Arguments
!
REAL, DIMENSION(0:ncloc+1,0:nrloc+1,2:nz), INTENT(INOUT) :: bgrad

!
! Name      Type    Purpose
!-------------------------------------------------------------------------------
!*bgrad*       REAL    Non-averaged buoyancy gradient                    [1/s^2]
!
!*******************************************************************************
!


RETURN

END SUBROUTINE buoyancy_frequency_sed

!========================================================================

SUBROUTINE dar_equation
!************************************************************************
!
! *dar_equation* Main program unit for dredging and relocation
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - coherens_main, initialise_model
!
!************************************************************************
!



RETURN

END SUBROUTINE dar_equation

!========================================================================

SUBROUTINE dar_init
!************************************************************************
!
! *dar_init* Initialise dredging and relocation parameters and arrays
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - coherens_main, initialise_model
!
!************************************************************************
!



RETURN

END SUBROUTINE dar_init

!========================================================================

SUBROUTINE deallocate_dar_arrays
!************************************************************************
!
! *deallocate_dar_arrays* Deallocate arrays for dredging/relocation model
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Reference -
!
! Calling program - simulation_end
!
!************************************************************************
!


RETURN

END SUBROUTINE deallocate_dar_arrays

!========================================================================

SUBROUTINE deallocate_morph_arrays
!************************************************************************
!
! *deallocate_morph_arrays* Deallocate arrays for the morphological model
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Reference -
!
! Calling program - simulation_end
!
!************************************************************************
!


RETURN

END SUBROUTINE deallocate_morph_arrays

!========================================================================

SUBROUTINE deallocate_sed_arrays
!************************************************************************
!
! *deallocate_sed_arrays* Deallocate arrays for the sediment transport model
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
! Reference -
!
! Calling program - simulation_end
!
!************************************************************************
!


RETURN

END SUBROUTINE deallocate_sed_arrays

!======================================================================

SUBROUTINE equation_of_state_sed
!************************************************************************
!
! *equation_of_state_sed* Density and expansion coefficients of a
!                         water-sediment mixture
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
! Calling program - equation_of_state
!
!************************************************************************
!


RETURN

END SUBROUTINE equation_of_state_sed

!========================================================================

SUBROUTINE exchange_sedics
!************************************************************************
!
! *exchange_sedics* Exchange initial conditions for the sediments
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE exchange_sedics

!========================================================================

SUBROUTINE initialise_sediment_arrays
!************************************************************************
!
! *initialise_sediment_arrays* Initialise arrays for the sediment model
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE initialise_sediment_arrays

!========================================================================

SUBROUTINE morphology_accel
!************************************************************************
!
! *morphology_equation* Main program unit of the morphological model
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - coherens_main, initialise_model
!
!************************************************************************
!



RETURN

END SUBROUTINE morphology_accel

!========================================================================

SUBROUTINE morphology_equation
!************************************************************************
!
! *morphology_accel* Update sediment and morphology using the tidal
!                    acceleration method
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - coherens_main, initialise_model
!
!************************************************************************
!



RETURN

END SUBROUTINE morphology_equation

!========================================================================

SUBROUTINE read_dar_spec
!************************************************************************
!
! *read_dar_spec* Read dredging site and ship properties in standard format
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE read_dar_spec

!========================================================================

SUBROUTINE read_morphics
!************************************************************************
!
! *read_morphics* Read initial conditions for the mrophological model
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE read_morphics

!========================================================================

SUBROUTINE read_sedics
!************************************************************************
!
! *read_sedics* Read initial conditions for the sediment transport model
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE read_sedics

!========================================================================

SUBROUTINE read_sed_spec
!************************************************************************
!
! *read_sed_spec* Read particle properties in standard format
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE read_sed_spec

!========================================================================

SUBROUTINE sediment_equation
!************************************************************************
!
! *sediment_equation* Main unit of the sediment transport model
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
! Calling program - coherens_main, initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE sediment_equation

!========================================================================

SUBROUTINE write_cif_vars_dar
!************************************************************************
!
! *write_cif_vars_dar* write dredging and relocation model parameters to a CIF
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x

! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!

IMPLICIT NONE



RETURN

END SUBROUTINE write_cif_vars_dar

!========================================================================

SUBROUTINE write_cif_vars_sed
!************************************************************************
!
! *write_cif_vars_sed* write sediment model parameters to a CIF
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!

IMPLICIT NONE



RETURN

END SUBROUTINE write_cif_vars_sed

!========================================================================

SUBROUTINE write_cif_vars_morph
!************************************************************************
!
! *write_cif_vars_morph* write morphology model parameters to a CIF
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!

IMPLICIT NONE



RETURN

END SUBROUTINE write_cif_vars_morph

!========================================================================

SUBROUTINE write_dar_spec
!************************************************************************
!
! *write_dar_spec* Write dredging/relocation attributes in standard format
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE write_dar_spec

!========================================================================

SUBROUTINE write_morphics
!************************************************************************
!
! *write_morphics* Write initial conditions for the morphological model
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - coherens_main, initialise_model
!
!
!************************************************************************
!

RETURN

END SUBROUTINE write_morphics

!========================================================================

SUBROUTINE write_sedics
!************************************************************************
!
! *write_sedics* Write initial conditions for the sediment model
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
! Calling program - coherens_main, initialise_model
!
!
!************************************************************************
!

RETURN

END SUBROUTINE write_sedics

!========================================================================

SUBROUTINE write_sed_spec
!************************************************************************
!
! *write_sed_spec* Write particle properties in standard format
!
! Author -
!
! Version - @(COHERENS)Sediment_Model.f90  V2.5
!
! Description - empty default routine
!
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE write_sed_spec
