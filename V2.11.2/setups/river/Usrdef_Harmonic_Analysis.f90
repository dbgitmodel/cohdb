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
! *Usrdef_Harmonic_Analysis* Define harmonic analysis and output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.4
!
! $Date: 2016-10-10 16:36:15 +0200 (Mon, 10 Oct 2016) $
!
! $Revision: 973 $
!
! Description - test case river
!
! Reference -
!
! Routines - usrdef_anal_freqs, usrdef_anal_params, usrdef_anal0d_vals,
!            usrdef_anal2d_vals, usrdef_anal3d_vals
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_anal_freqs
!************************************************************************
!
! *usrdef_anal_freqs* Harmonic frequencies
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.0
!
! Description - test case river
!
! Reference -
!
! Calling program - harmonic_analysis_init
!
!************************************************************************
!
USE iopars
USE tide
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_anal_freqs'
CALL log_timer_in()

!
!1. Harmonic frequencies
!-----------------------
!

index_anal(1) = icon_S2

!
!3. Specifiers for harmonic analysis
!-----------------------------------
!
!---number of frequencies per set
nofreqsharm(1) = 1

!---frequency indices
ifreqsharm(1,1) = 1

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_anal_freqs

!========================================================================

SUBROUTINE usrdef_anal_params
!************************************************************************
!
! *usrdef_anal_params* Specifiers for harmonic output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.4
!
! Description - test case river
!
! Reference -
!
! Calling program - harmonic_analysis_init
!
!************************************************************************
!
USE gridpars
USE iopars
USE modids
USE switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_anal_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!

analvars%ivarid = (/iarr_uvel,iarr_vvel,iarr_wphys,iarr_sal/)
analvars%nrank = 3

!
!3. Variable indices
!-------------------
!

ivarsanal(1,1:4) = (/1,2,3,4/) 

!
!4. File parameters
!------------------
!

res3d(1)%defined = .TRUE.
res3d(1)%info = .TRUE.

!
!5. Output grid
!--------------
!

analgpars(1)%time_format = 4
analgpars(1)%tlims(1:2) = (/0,nstep/)
analgpars(1)%tlims(3) = 2880

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_anal_params

!========================================================================

SUBROUTINE usrdef_anal0d_vals(out0ddat,n0vars)
!************************************************************************
!
! *usrdef_anal0d_vals* 0-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.0
!
! Description - test case river
!
! Reference -
!
! Calling program - harmonic_analysis_update
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

      
RETURN

END SUBROUTINE usrdef_anal0d_vals

!========================================================================

SUBROUTINE usrdef_anal2d_vals(out2ddat,i,j,n2vars)
!************************************************************************
!
! *usrdef_anal2d_vals* 2-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.0
!
! Description - test case river
!
! Reference -
!
! Calling program - harmonic_analysis_update
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


RETURN

END SUBROUTINE usrdef_anal2d_vals

!========================================================================

SUBROUTINE usrdef_anal3d_vals(out3ddat,i,j,k,n3vars)
!************************************************************************
!
! *usrdef_anal3d_vals* 3-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.1.2
!
! Description - test case river
!
! Reference -
!
! Calling program - harmonic_analysis_update
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


RETURN

END SUBROUTINE usrdef_anal3d_vals
