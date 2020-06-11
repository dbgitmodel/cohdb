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
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.11.2
!
! $Date: 2016-10-10 16:36:15 +0200 (Mon, 10 Oct 2016) $
!
! $Revision: 973 $
!
! Description - test case papa
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
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.11.2
!
! Description - test case papa
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
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: l


procname(pglev+1) = 'usrdef_avr_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!
!---user-defined
avrvars(1)%f90_name = 'tmld'
avrvars(1)%standard_name = 'mixed_layer_depth'
avrvars(1)%long_name = 'Mixed layer depth'
avrvars(1)%units = 'm'

!---varids and ranks
avrvars%ivarid = (/0,iarr_qrad,iarr_qlwave,iarr_qlatent,iarr_qsensible,&
                 & iarr_qnonsol,iarr_qtot,iarr_temp/)
avrvars(1:novarstsr-1)%nrank = 2
avrvars(novarstsr)%nrank = 3

!
!2. Variable indices
!-------------------
!

ivarsavr(1,:) = (/(l,l=1,novarstsr)/)
ivarsavr(2,:) = (/(l,l=1,novarstsr)/)

!
!3. File parameters
!------------------
!

avr2d(1:2)%defined = .TRUE.
avr3d(1:2)%defined = .TRUE.

!
!4. Output grid
!--------------
!
!---first file (weekly)
avrgpars(1)%time_format = 5
avrgpars(1)%xlims = (/2,2,1/)
avrgpars(1)%ylims = (/2,2,1/)
avrgpars(1)%zlims = (/1,nz,1/)
avrgpars(1)%tlims = (/0,nstep-720,1008/)

!---second file (monthly)
avrgpars(2)%time_format = 6
avrgpars(2)%xlims = (/2,2,1/)
avrgpars(2)%ylims = (/2,2,1/)
avrgpars(2)%zlims = (/1,nz,1/)
avrgpars(2)%tlims = (/0,nstep,4383/)


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
! Description -
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
! Version - @(COHERENS)Usrdef_Time_Averages.f90  V2.11.2
!
! Description - test case papa
!
! Reference -
!
! Calling program - time_averages_update
!
!************************************************************************
!
USE density
USE depths  
USE grid
USE gridpars

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
!*Local variables
!
REAL, PARAMETER :: rho_crit = 0.01
INTEGER :: k, kk
REAL :: slo, sup


kk = 0
DO k=1,nz-1
   IF ((dens(i,j,k)-dens(i,j,nz)).GT.rho_crit) kk = k
ENDDO
IF (kk.EQ.0) THEN
   out2ddat(1) = deptotatc(i,j)
ELSE
   slo = 1.0 - gscoordatc(i,j,kk) 
   sup = 1.0 - gscoordatc(i,j,kk+1)
   out2ddat(1) = deptotatc(i,j)*((sup-slo)*(rho_crit+dens(i,j,nz))+&
               & slo*dens(i,j,kk+1)-sup*dens(i,j,kk))&
               & /(dens(i,j,kk+1)-dens(i,j,kk))
ENDIF


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
! Description -
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
