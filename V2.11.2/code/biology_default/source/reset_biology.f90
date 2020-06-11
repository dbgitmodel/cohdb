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

MODULE reset_biology
!************************************************************************
!
! *reset_biology* Reset biological model parameters and arrays if needed
!
! Author -
!
! Version - @(COHERENS)reset_biology.f90  V2.5.1
!
! $Date: 2015-01-16 16:39:38 +0100 (Fri, 16 Jan 2015) $
!
! $Revision: 796 $
!
! Description - routines are empty but need to be declared   
!
! Routines - reset_bioics, reset_bio_params
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS


!========================================================================

SUBROUTINE reset_bioics
!************************************************************************
!
! *reset_bioics* Reset initial conditions in the biological model
!
! Author -
!
! Version - @(COHERENS)reset_biology.f90  V2.5.1
!
! Description - empty default routine 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE reset_bioics

!========================================================================

SUBROUTINE reset_bio_params
!************************************************************************
!
! *reset_sed_params* Reset biological model parameters
!
! Author -
!
! Version - @(COHERENS)reset_biology.f90  V2.5.1
!
! Description - empty default routine 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE reset_bio_params

END MODULE reset_biology
