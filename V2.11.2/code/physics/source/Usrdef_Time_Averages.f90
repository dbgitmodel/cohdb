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
! *Usrdef_Time_Averages* Define time averaged output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.0
!
! $Date: 2015-04-10 11:30:34 +0200 (Fri, 10 Apr 2015) $
!
! $Revision: 841 $
!
! Description - empty default file
!
! Reference -
!
! Routines - usrdef_avr_params, usrdef_avr0d_vals, usrdef_avr2d_vals,
!            usrdef_avr3d_vals
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_avr_params
!************************************************************************
!
! *usrdef_avr_params* Specifiers for time averaged output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.0
!
! Description - empty default routine
!             - called by all processes
!
! Reference -
!
! Calling program - time_averages_init
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_avr_params

!========================================================================

SUBROUTINE usrdef_avr0d_vals(out0ddat,n0vars)
!************************************************************************
!
! *usrdef_avr0d_vals* 0-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.0
!
! Description - empty default routine
!             - called by all processes
!
! Reference -
!
! Calling program - time_averages_update
!
!************************************************************************
!
IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: n0vars
REAL, INTENT(OUT), DIMENSION(n0vars) :: out0ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out0ddat* REAL    output data
!*n0vars*   INTEGER number of 0-D variables
!
!------------------------------------------------------------------------------
!


IF (n0vars.GT.0) out0ddat = 0.0


RETURN

END SUBROUTINE usrdef_avr0d_vals

!========================================================================

SUBROUTINE usrdef_avr2d_vals(out2ddat,i,j,n2vars)
!************************************************************************
!
! *usrdef_avr2d_vals* 2-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.0
!
! Description - empty default routine
!             - called by all processes
!
! Reference -
!
! Calling program - time_averages_update
!
!************************************************************************
!

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, n2vars
REAL, INTENT(OUT), DIMENSION(n2vars) :: out2ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out2ddat* REAL    output data
!*i*        INTEGER X-index of output location
!*j*        INTEGER Y-index of output location
!*n2vars*   INTEGER number of 2-D variables
!
!------------------------------------------------------------------------------
!


IF (n2vars.GT.0) out2ddat = 0.0


RETURN

END SUBROUTINE usrdef_avr2d_vals

!========================================================================

SUBROUTINE usrdef_avr3d_vals(out3ddat,i,j,k,n3vars)
!************************************************************************
!
! *usrdef_avr3d_vals* 3-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.0
!
! Description - empty default routine
!             - called by all processes
!
! Reference -
!
! Calling program - time_averages_update
!
!************************************************************************
!

IMPLICIT NONE
!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, k, n3vars
REAL, INTENT(OUT), DIMENSION(n3vars) :: out3ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out3ddat* REAL    output data
!*i*        INTEGER X-index of output location
!*j*        INTEGER Y-index of output location
!*k*        INTEGER Vertical index of output location
!*n3vars*   INTEGER number of 3-D variables
!
!------------------------------------------------------------------------------
!


IF (n3vars.GT.0) out3ddat = 0.0


RETURN

END SUBROUTINE usrdef_avr3d_vals
