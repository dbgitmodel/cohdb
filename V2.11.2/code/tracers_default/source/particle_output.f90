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

MODULE particle_output
!************************************************************************
!
! *particle_output* Routines defining particle output data
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_output.f90  V2.x
!
! $Date: 2016-07-12 15:51:17 +0200 (Tue, 12 Jul 2016) $
!
! $Revision: 943 $
!
! Description - routines are empty but need to be declared
!
! Reference -
!
! Routines - define_part0d_vals, define_part2d_vals, define_part3d_vals
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS

!========================================================================

SUBROUTINE define_part0d_vals(outdat,ivarid,l,oopt)
!************************************************************************
!
! *define_part0d_vals* Define 0-D particle model output data
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_output.f90  V2.x
!
! Description - empty default routine
!
! Calling program - define_out0d_vals
!
! Module calls -
!  
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: ivarid, l, oopt
REAL, INTENT(OUT) :: outdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    0-D output data value
!*ivarid*  INTEGER Variable key id
!*l*       INTEGER Contaminant label
!*oopt*    INTEGER Type of operator
!
!------------------------------------------------------------------------------
!

RETURN

END SUBROUTINE define_part0d_vals

!========================================================================

SUBROUTINE define_part2d_vals(outdat,i,j,ivarid,l,oopt,kc,dep)
!************************************************************************
!
! *define_part2d_vals*  Define particle model output data at a 2-D grid location
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_output.f90  V2.x
!
! Description - empty default routine
!
! Calling program - define_out2d_vals
!
! Module calls -
!
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: i, ivarid, j, kc, l, oopt
REAL, INTENT(IN) :: dep
REAL, INTENT(OUT) :: outdat
!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    2-D output data value
!*i*       INTEGER (Local) X-index of the data location
!*j*       INTEGER (Local) Y-index of the data location
!*ivarid*  INTEGER Variable key id
!*f*       INTEGER Contaminant label
!*oopt*    INTEGER Type of operator
!*kc*      INTEGER Vertical level at which a 3-D variable is evaluated
!*dep*     REAL    Depth below the surface at which a 3-D variable is evaluated
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE define_part2d_vals

!========================================================================

SUBROUTINE define_part3d_vals(outdat,i,j,k,ivarid,l)
!************************************************************************
!
! *define_part3d_vals*  Define particle model output data at a 3-D grid location
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_output.f90  V2.x
!
! Description - empty default routine
!
! Calling program - define_out2d_vals, define_out3d_vals
!
! Module calls -
!
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: i, ivarid, j, k, l
REAL, INTENT(OUT) :: outdat
!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    3-D output data value
!*i*       INTEGER (Local) X-index of the data location
!*j*       INTEGER (Local) Y-index of the data location
!*k*       INTEGER Vertical index of the data location
!*ivarid*  INTEGER Variable key id
!*l*       INTEGER Contaminant label
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE define_part3d_vals


END MODULE particle_output
