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

MODULE sediment_output
!************************************************************************
!
! *sediment_output* Routines defining sediment output data
!
! Author -
!
! Version - @(COHERENS)sediment_output.f90  V2.x
!
! $Date: 2015-01-16 16:39:38 +0100 (Fri, 16 Jan 2015) $
!
! $Revision: 796 $
!
! Description - routines are empty but need to be declared 
!
! Routines - define_sed0d_vals, define_sed2d_vals, define_sed3d_vals
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS

!========================================================================

SUBROUTINE define_sed0d_vals(outdat,ivarid,f,oopt)
!************************************************************************
!
! *define_sed0d_vals* Define 0-D sediment output data
!
! Author -
!
! Version - @(COHERENS)sediment_output.f90  V2.x
!
! Calling program - define_out0d_vals
!
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: f, ivarid, oopt
REAL, INTENT(OUT) :: outdat
!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    0-D output data value
!*ivarid*  INTEGER Variable key id
!*f*       INTEGER Fraction number
!*oopt*    INTEGER Type of operator
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE define_sed0d_vals

!========================================================================

SUBROUTINE define_sed2d_vals(outdat,i,j,ivarid,f,oopt,kc,dep)
!************************************************************************
!
! *define_sed2d_vals*  Define sediment output data at a 2-D grid location
!
! Author -
!
! Version - @(COHERENS)sediment_output.f90  V2.x
!
! Description - empty default routine
!
! Calling program - define_out2d_vals
!
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: f, i, ivarid, j, kc, oopt
REAL, INTENT(IN) :: dep
REAL, INTENT(OUT) :: outdat
!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    2-D output data value
!*i*       INTEGER (Local) X-index of the data location
!*j*       INTEGER (Local) Y-index of the data location
!*ivarid*  INTEGER Variable key id
!*f*       INTEGER Fraction number
!*oopt*    INTEGER Type of operator
!*kc*      INTEGER Vertical level at which a 3-D variable is evaluated
!*dep*     REAL    Depth below the surface at which a 3-D variable is evaluated
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE define_sed2d_vals

!========================================================================

SUBROUTINE define_sed3d_vals(outdat,i,j,k,ivarid,f)
!************************************************************************
!
! *define_sed3d_vals*  Define sediment output data at a 3-D grid location
!
! Author -
!
! Version - @(COHERENS)sediment_output.f90  V2.x
!
! Description - empty default routine
!
! Calling program - define_out2d_vals, define_out3d_vals
!
!************************************************************************
!
USE syspars

!
!*Arguments
!
INTEGER, INTENT(IN) :: f, i, ivarid, j, k
REAL, INTENT(OUT) :: outdat
!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    3-D output data value
!*i*       INTEGER (Local) X-index of the data location
!*j*       INTEGER (Local) Y-index of the data location
!*k*       INTEGER Vertical index of the data location
!*ivarid*  INTEGER Variable key id
!*f*       INTEGER Fraction number
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE define_sed3d_vals


END MODULE sediment_output
