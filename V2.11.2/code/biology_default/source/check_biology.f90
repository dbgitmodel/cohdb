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

MODULE check_biology
!************************************************************************
!
! *check_biology* Check biological model parameters and arrays for errors
!
! Author -
!
! Version - @(COHERENS)check_biology.f90  V2.5.1
!
! $Date: 2015-01-16 16:39:38 +0100 (Fri, 16 Jan 2015) $
!
! $Revision: 796 $
!
! Description - routines are empty but need to be declared  
!
! Routines - check_bioics, check_bio_params
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS


!========================================================================

SUBROUTINE check_bioics
!************************************************************************
!
! *check_bioics* Check initial conditions
!
! Author -
!
! Version - @(COHERENS)check_biology.f90  V2.5.1
!
! Description - empty default routine 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE check_bioics

!========================================================================

SUBROUTINE check_bio_params
!************************************************************************
!
! *check_bio_params* Check biological parameters
!
! Author -
!
! Version - @(COHERENS)check_biology.f90  V2.5.1
!
! Description - empty default routine 
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE check_bio_params


END MODULE check_biology
