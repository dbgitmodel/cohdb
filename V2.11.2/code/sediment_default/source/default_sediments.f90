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

MODULE default_sediments
!************************************************************************
!
! *default_sediments* Default settings for the sediment model
!
! Author -
!
! Version - @(COHERENS)default_sediments.f90  V2.x
!
! $Date: 2015-01-16 16:39:38 +0100 (Fri, 16 Jan 2015) $
!
! $Revision: 796 $
!
! Description - routines are empty but need to be declared   
!
! Routines - default_dar_params, default_morphics, default_morph_params,
!            default_sedics, default_sed_params, default_sed_spec
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS


!========================================================================

SUBROUTINE default_sed_params
!************************************************************************
!
! *default_sed_params* Default settings of parameters for the sediment
!                      transport model
!
! Author -
!
! Version - @(COHERENS)default_sediments.f90  V2.5
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

END SUBROUTINE default_sed_params

!========================================================================

SUBROUTINE default_morph_params
!************************************************************************
!
! *default_morph_params* Default settings of parameters for the morphological
!                        model
!
! Author -
!
! Version - @(COHERENS)default_sediments.f90  V2.x
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

END SUBROUTINE default_morph_params

!========================================================================

SUBROUTINE default_dar_params
!************************************************************************
!
! *default_dar_params* Default settings of parameters for dredging and
!                      relocation
!
! Author -
!
! Version - @(COHERENS)default_sediments.f90  V2.x
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

END SUBROUTINE default_dar_params

!========================================================================

SUBROUTINE default_sedics
!************************************************************************
!
! *default_sedics* Default settings of initial conditions for the sediment
!                  transport model
!
! Author -
!
! Version - @(COHERENS)default_sediments.f90  V2.5
!
! Description - empty default routine 
!
! Calling program - initialise model
!
!************************************************************************
!


RETURN

END SUBROUTINE default_sedics

!========================================================================

SUBROUTINE default_morphics
!************************************************************************
!
! *default_morphics* Default settings of initial conditions for the sediment
!                    morphological model
!
! Author -
!
! Version - @(COHERENS)default_sediments.f90  V2.x
!
! Description - empty default routine 
!
! Calling program - initialise model
!
!************************************************************************
!


RETURN

END SUBROUTINE default_morphics

!========================================================================

SUBROUTINE default_sed_spec
!************************************************************************
!
! *default_sed_spec* Default settings of sediment particle attributes
!
! Author -
!
! Version - @(COHERENS)default_sediments.f90  V2.x
!
! Description - empty default routine 
!
! Calling program - initialise model
!
!************************************************************************
!


RETURN

END SUBROUTINE default_sed_spec


END MODULE default_sediments
