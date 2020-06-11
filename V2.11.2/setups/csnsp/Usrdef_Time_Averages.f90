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
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.1.2
!
! $Date: 2016-10-10 16:36:15 +0200 (Mon, 10 Oct 2016) $
!
! $Revision: 973 $
!
! Description - test case csnsp
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
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.1.1
!
! Description - test case csnsp
!
! Reference -
!
! Calling program - time_averages_init
!
!************************************************************************
!
USE datatypes
USE gridpars
USE iopars
USE modids
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_avr_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!

avrvars%ivarid = (/iarr_sstresatc,iarr_qnonsol,iarr_cds,iarr_temp/)
avrvars%nrank = (/2,2,2,3/)

!
!2. Variable indices
!-------------------
!

ivarsavr(1,1:4) = (/1,2,3,4/)

!
!3. File parameters
!------------------
!

avr2d(1)%defined = .TRUE.
avr3d(1)%defined = .TRUE.

!
!4. Output grid
!--------------
!

avrgpars(1)%time_format = 4
avrgpars(1)%refdate = '1989/01/01;00:00:00:000'
avrgpars(1)%xlims = (/2,2,1/)
avrgpars(1)%ylims = (/2,2,1/)
avrgpars(1)%zlims = (/1,nz,1/)
avrgpars(1)%tlims = (/0,21600,216/)

CALL log_timer_out()


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
! Description - test case csnsp
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
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.1.2
!
! Description - test case csnsp
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
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.1.2
!
! Description - test case csnsp
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


RETURN

END SUBROUTINE usrdef_avr3d_vals
