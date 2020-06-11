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
! *Particle_Model* External particle routines called within the main COHERENS
!                  code
!
! Author -
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! $Date: 2017-01-17 13:39:59 +0100 (Tue, 17 Jan 2017) $
!
! $Revision: 1000 $
!
! Description - routines are empty but need to be declared 
!
! Routines -  allocate_part_arrays, assign_cif_vars_part,
!             assign_cif_vars_tspart, combine_particle_grid,
!             combine_particle_phys_data, deallocate_part_arrays,
!             define_partics, define_particle_phys_data,
!             particle_concentrations, particle_drift, particle_model,
!             particle_trajects, particle_trajects_init, read_part_spec,
!             update_particle_phys_data, write_cif_vars_part,
!             write_cif_vars_tspart, write_partics, write_particle_grid,
!             write_particle_phys_data, write_part_spec
!
!************************************************************************
!

!========================================================================

SUBROUTINE allocate_part_arrays
!************************************************************************
!
! *allocate_part_arrays*  Allocate arrays for the particle module
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE allocate_part_arrays

!============================================================================

SUBROUTINE assign_cif_vars_part(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_part* convert the string data a particle model CIF
!                        input line to the appropriate numeric or non-numeric
!                        format 
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Reference -
!
! Calling program - read_cif_params
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

END SUBROUTINE assign_cif_vars_part

!============================================================================

SUBROUTINE assign_cif_vars_tspart(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_tspart* convert string data from an input line in the CIF
!                          block with parameters for particle trajectory output
!                          to the appropriate numeric or non-numeric format
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Reference -
!
! Calling program - read_cif_params
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

END SUBROUTINE assign_cif_vars_tspart

!========================================================================

SUBROUTINE combine_particle_grid
!************************************************************************
!
! *combine_particle_grid*  Send model grid data from COHERENS to the particle
!                          module
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE combine_particle_grid

!========================================================================

SUBROUTINE combine_particle_phys_data
!************************************************************************
!
! *combine_particle_phys_data* Send physical forcing data from COHERENS to the
!                              particle module
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - coherens_main, initialise_model, initialise_particle_model,
!                   particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE combine_particle_phys_data

!========================================================================

SUBROUTINE deallocate_part_arrays
!************************************************************************
!
! *deallocate_part_arrays*  Deallocate arrays for the particle module
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - simulation_end
!
!************************************************************************
!


RETURN

END SUBROUTINE deallocate_part_arrays

!========================================================================

SUBROUTINE define_partics
!************************************************************************
!
! *define_partics* Define initial conditions for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE define_partics

!========================================================================

SUBROUTINE define_particle_phys_data
!************************************************************************
!
! *define_particle_phys_data* Define physical arrays for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE define_particle_phys_data

!========================================================================

SUBROUTINE particle_concentrations
!************************************************************************
!
! *particle_concentration* convert particle densities into Eulerian
!                          concentrations
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - coherens_main, particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE particle_concentrations

!========================================================================

SUBROUTINE particle_drift
!************************************************************************
!
! *particle_drift* Particle drift model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - coherens_main, particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE particle_drift

!========================================================================

SUBROUTINE particle_model
!************************************************************************
!
! *particle_model* COHERENS particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.F90  V2.x
!
! $Date: 2017-01-17 13:39:59 +0100 (Tue, 17 Jan 2017) $
!
! $Revision: 1000 $
!
! Description - empty default routine
!
! Reference - routine is called if the particle model runs offline or online in
!             parallel mode
!
! Calling program - coherens_main
!
!************************************************************************
!

  
RETURN

END SUBROUTINE particle_model
  
!========================================================================

SUBROUTINE particle_trajects
!************************************************************************
!
! *particle_trajects_init* Initialise parameters for particle trajectory
!                          output
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - coherens_main, particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE particle_trajects
  
!========================================================================

SUBROUTINE particle_trajects_init
!************************************************************************
!
! *particle_trajects_init* Particle trajectory output
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - coherens_main, particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE particle_trajects_init

!========================================================================

SUBROUTINE read_part_spec
!************************************************************************
!
! *read_part_spec* Read specifiers arrays for particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE read_part_spec

!========================================================================

SUBROUTINE update_particle_phys_data
!************************************************************************
!
! *update_particle_phys_data* update physical data for the particle module
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model, initialise_particle_model, particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE update_particle_phys_data

!========================================================================

SUBROUTINE write_cif_vars_part
!************************************************************************
!
! *write_cif_vars_part* write particle model parameters to a CIF
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE write_cif_vars_part

!========================================================================

SUBROUTINE write_cif_vars_tspart
!************************************************************************
!
! *write_cif_vars_tspart* Write the CIF block with parameters for particle
!                         trajectory output
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE write_cif_vars_tspart

!========================================================================

SUBROUTINE write_partics
!************************************************************************
!
! *write_partics* Write initial conditions for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - coherens_main, initialise_model, initialise_particle_model,
!                   particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE write_partics

!========================================================================

SUBROUTINE write_particle_grid
!************************************************************************
!
! *write_particle_grid* Write model grid arrays for the particle module
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE write_particle_grid

!========================================================================

SUBROUTINE write_particle_phys_data(ciodatetime)
!************************************************************************
!
! *write_particle_phys_data* write physical data for the particle module
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model, coherens_main
!
!************************************************************************
!
USE syspars

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: ciodatetime

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*ciodatetime* CHAR     Date/time in data file
!
!------------------------------------------------------------------------------
!
RETURN

END SUBROUTINE write_particle_phys_data

!========================================================================

SUBROUTINE write_part_spec
!************************************************************************
!
! *write_part_spec* Write specifiers arrays for particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.x
!
! Description - empty default routine
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!


RETURN

END SUBROUTINE write_part_spec

