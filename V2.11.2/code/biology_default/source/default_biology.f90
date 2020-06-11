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

MODULE default_biology
!************************************************************************
!
! *default_biology* Default settings for the biological model
!
! Author -
!
! Version - @(COHERENS)default_biology.f90  V2.5.1
!
! $Date: 2015-01-16 16:39:38 +0100 (Fri, 16 Jan 2015) $
!
! $Revision: 796 $
!
! Description - routines are empty but need to be declared   
!
! Routines - default_bioics, default_bio_params
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS


!========================================================================

SUBROUTINE default_bioics
!************************************************************************
!
! *default_bioics* Default settings for biological initial conditions
!
! Author -
!
! Version - @(COHERENS)default_biology.f90  V2.5.1
!
! Description - empty default routine 
!
! Calling program - initialise model
!
!************************************************************************
!


RETURN

END SUBROUTINE default_bioics

!========================================================================

SUBROUTINE default_bio_params
!************************************************************************
!
! *default_bio_params* Default settings for biological model parameters
!
! Author -
!
! Version - @(COHERENS)default_biology.f90  V2.5.1
!
! $Date: 2015-01-16 16:39:38 +0100 (Fri, 16 Jan 2015) $
!
! $Revision: 796 $
!
! Description - empty default routine
!
! Calling program - initialise model
!
!************************************************************************
!


RETURN

END SUBROUTINE default_bio_params


END MODULE default_biology
