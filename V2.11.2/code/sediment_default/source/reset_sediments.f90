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

MODULE reset_sediments
!************************************************************************
!
! *reset_sediments* Reset sediment model parameters and arrays if needed
!
! Author -
!
! Version - @(COHERENS)reset_sediments.f90  V2.x
!
! $Date: 2015-01-16 16:39:38 +0100 (Fri, 16 Jan 2015) $
!
! $Revision: 796 $
!
! Description - routines are empty but need to be declared   
!
! Routines - reset_sed_params, reset_morph_params, reset_dar_params,
!            reset_sedics
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS


!========================================================================

SUBROUTINE reset_sed_params
!************************************************************************
!
! *reset_sed_params* Reset parameters for the sediment transport model
!
! Author -
!
! Version - @(COHERENS)reset_sediments.f90  V2.5
!
! Description - empty default routine 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE reset_sed_params

!========================================================================

SUBROUTINE reset_morph_params
!************************************************************************
!
! *reset_morph_params* Reset parameters for the morphological model
!
! Author -
!
! Version - @(COHERENS)reset_sediments.f90  V2.x
!
! Description - empty default routine 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE reset_morph_params

!========================================================================

SUBROUTINE reset_dar_params
!************************************************************************
!
! *reset_dar_params* Reset parameters for the dredging/relocation model
!
! Author -
!
! Version - @(COHERENS)reset_sediments.f90  V2.x
!
! Description - empty default routine 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE reset_dar_params

!========================================================================

SUBROUTINE reset_sedics
!************************************************************************
!
! *reset_sedics* Reset initial conditions for the sediment transport model
!
! Author -
!
! Version - @(COHERENS)reset_sediments.f90  V2.5
!
! Description - empty default routine 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE reset_sedics

END MODULE reset_sediments
