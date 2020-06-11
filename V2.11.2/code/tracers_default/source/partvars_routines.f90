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
! Version - @(COHERENS)partvars_routines.f90  V2.x
!
! $Date: 2017-12-22 12:37:21 +0100 (Fri, 22 Dec 2017) $
!
! $Revision: 1070 $
!
! Description - routines are empty but need to be declared
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
! Version - @(COHERENS)partvars_routines.f90  V2.x
!
! Description - empty default routine
!
! Calling program - inquire_var
!
! Module calls -
!
!************************************************************************
!  
USE datatypes
  
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
! Version - @(COHERENS)partvars_routines.f90  V2.x
!
! Description - empty default routine
!
! Calling program - set_modfiles_atts
!
! Module calls -
!
!************************************************************************
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
! Version - @(COHERENS)sedvars_routines.f90  V2.x
!
! Description - empty default routine
!
! Calling progam - set_modvars_atts
!
! Module calls -
!
!************************************************************************
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


RETURN

END SUBROUTINE set_partvars_atts

END MODULE partvars_routines
