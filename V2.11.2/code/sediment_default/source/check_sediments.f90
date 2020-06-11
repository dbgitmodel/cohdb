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

MODULE check_sediments
!************************************************************************
!
! *reset_check_sed* Check sediment model parameters and arrays for errors
!
! Author -
!
! Version - @(COHERENS)check_sediments.f90  V2.5
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - routines are empty but need to be declared
!
! Routines - check_dar_params, check_dar_spec, check_morphics,
!            check_morph_params, check_sedics, check_sed_params
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS


!========================================================================

SUBROUTINE check_sed_params
!************************************************************************
!
! *check_sed_params* Check parameters for the sediment transport model
!
! Author -
!
! Version - @(COHERENS)check_sediments.f90  V2.5
!
! Description - empty default routine 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE check_sed_params

!========================================================================

SUBROUTINE check_morph_params
!************************************************************************
!
! *check_morph_params* Check parameters for the morphological model
!
! Author -
!
! Version - @(COHERENS)check_sediments.f90  V2.x
!
! Description - empty default routine 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE check_morph_params

!========================================================================

SUBROUTINE check_dar_params
!************************************************************************
!
! *check_dar_params* Check parameters for dredging/relocation
!
! Author -
!
! Version - @(COHERENS)check_sediments.f90  V2.x
!
! Description - empty default routine 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE check_dar_params

!========================================================================

SUBROUTINE check_dar_spec
!************************************************************************
!
! *check_dar_spec* Check attribute arrays for dredging/relocation
!
! Author -
!
! Version - @(COHERENS)check_sediments.f90  V2.x
!
! Description - empty default routine 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE check_dar_spec

!========================================================================

SUBROUTINE check_sedics
!************************************************************************
!
! *check_sedics* Check initial conditions for the sediment transport model
!
! Author -
!
! Version - @(COHERENS)check_sediments.f90  V2.5
!
! Description - empty default routine 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE check_sedics

!========================================================================

SUBROUTINE check_morphics
!************************************************************************
!
! *check_morphics* Check initial conditions for the morphological model
!
! Author -
!
! Version - @(COHERENS)check_sediments.f90  V2.x
!
! Description - empty default routine 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE check_morphics


END MODULE check_sediments
