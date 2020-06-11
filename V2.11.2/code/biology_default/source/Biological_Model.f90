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
! *Biological_Model* external routines for the biological model called within
!                    the main COHERENS code
!
! Author -
!
! Version - @(COHERENS)Biology_Model.f90  V2.5.1
!
! $Date: 2015-01-16 16:39:38 +0100 (Fri, 16 Jan 2015) $
!
! $Revision: 796 $
!
! Description - routines are empty but need to be declared 
!
! Routines - allocate_bio_arrays, assign_cif_vars_bio, biological_model,
!    deallocate_bio_arrays, exchange_bioics, initialise_biology_arrays,
!    read_bioics, read_bio_spec, write_bioics, write_bio_spec,
!    write_cif_vars_bio
!
!************************************************************************
!

!========================================================================

SUBROUTINE allocate_bio_arrays
!************************************************************************
!
! *allocate_bio_arrays* Allocate biological arrays
!
! Author -
!
! Version - @(COHERENS)Biological_Model.f90  V2.5.1
!
! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE allocate_bio_arrays

!============================================================================

SUBROUTINE assign_cif_vars_bio(iddesc,cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_bio* convert the string data from a biological CIF input
!                       line to the appropriate numeric or non-numeric format 
!
! Author -
!
! Version - @(COHERENS)Biological_Model.f90  V2.5.1
!
! Description - empty default routine
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
INTEGER, INTENT(IN) :: iddesc, numvars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER CIF file key id
!*cname*     CHAR    Name of the variable(s)
!*cvals*     CHAR    Data values on input data line
!*numvars*   INTEGER Number of data variables
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE assign_cif_vars_bio

!========================================================================

SUBROUTINE biological_model
!************************************************************************
!
! *biological_model* Main unit of the biological model
!
! Author -
!
! Version - @(COHERENS)Biological_Model.f90  V2.5.1
!
! Description - empty default routine
!
! Calling program - coherens_main, initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE biological_model

!========================================================================

SUBROUTINE deallocate_bio_arrays
!************************************************************************
!
! *deallocate_bio_arrays* Deallocate biological arrays
!
! Author -
!
! Version - @(COHERENS)Biological_Model.f90  V2.5.1
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

END SUBROUTINE deallocate_bio_arrays

!========================================================================

SUBROUTINE exchange_bioics
!************************************************************************
!
! *exchange_bioics* Exchange initial conditions in the biological model
!
! Author -
!
! Version - @(COHERENS)Biological_Model.f90  V2.5.1
!
! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE exchange_bioics

!========================================================================

SUBROUTINE initialise_biology_arrays
!************************************************************************
!
! *initialise_biology_arrays* Initialise arrays for the biological model
!
! Author -
!
! Version - @(COHERENS)Biological_Model.f90  V2.5.1
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

END SUBROUTINE initialise_biology_arrays

!========================================================================

SUBROUTINE read_bioics
!************************************************************************
!
! *read_bioics* Read initial conditions for the biological model in standard 
!               format
!
! Author -
!
! Version - @(COHERENS)Biological_Model.f90  V2.5.1
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

END SUBROUTINE read_bioics

!========================================================================

SUBROUTINE read_bio_spec
!************************************************************************
!
! *read_bio_spec* Read specifgier arrays for the biological model in
!                 standard format
!
! Author -
!
! Version - @(COHERENS)Biological_Model.f90  V2.5.1
!
! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE read_bio_spec

!========================================================================

SUBROUTINE write_bioics
!************************************************************************
!
! *write_bioics* Write initial conditions for the biological model
!
! Author -
!
! Version - @(COHERENS)Biological_Model.f90  V2.5.1
!
! Description - empty default routine
!
! Calling program - coherens_main, initialise_model
!
!
!************************************************************************
!

RETURN

END SUBROUTINE write_bioics

!========================================================================

SUBROUTINE write_bio_spec
!************************************************************************
!
! *write_bio_spec* Write specifier arrays for the biological model in
!                  standard format
!
! Author -
!
! Version - @(COHERENS)Biological_Model.f90  V2.5.1
!
! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE write_bio_spec

!========================================================================

SUBROUTINE write_cif_vars_bio
!************************************************************************
!
! *write_cif_vars_bio* write biological model parameters to a CIF
!
! Author -
!
! Version - @(COHERENS)Biological_Model.f90  V2.5.1
!
! Description - empty default routine
!
! Calling program - initialise_model
!
!************************************************************************
!

IMPLICIT NONE



RETURN

END SUBROUTINE write_cif_vars_bio
