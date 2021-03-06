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
! *Usrdef_Time_Series* Define time series output
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.8
!
! $Date: 2015-10-06 17:42:25 +0200 (Tue, 06 Oct 2015) $
!
! $Revision: 886 $
!
! Last update - 25 May 2011  @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description - test case drythin
!
! Reference -
!
! Routines - usrdef_tsr_params, usrdef_tsr0d_vals, usrdef_tsr2d_vals,
!            usrdef_tsr3d_vals
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_tsr_params
!************************************************************************
!
! *usrdef_tsr_params* Specifiers for time series output
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.8
!
! Description - test case drythin
!
! Reference -
!
! Calling program - time_series_init
!
!************************************************************************
!
USE iopars
USE modids
USE switches
USE timepars 
USE modvars_routines, ONLY: inquire_var
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: icount

procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!
!---varids and ranks
tsrvars%ivarid = (/iarr_zeta,iarr_umvel,iarr_vmvel,iarr_vortic2d,&
                 & iarr_hmvelmag,iarr_uvel,iarr_vvel,iarr_wphys,iarr_hvelmag/)
tsrvars%nrank = (/2,2,2,2,2,3,3,3,3/)

!
!2. Variable indices
!-------------------
!

ivarstsr(1,1:9) = (/1,2,3,4,5,6,7,8,9/)

!
!3. File parameters
!------------------
!

tsr2d(1)%defined = .TRUE.
tsr3d(1)%defined = .TRUE.

!
!4. Output grid
!--------------
!
!---first file
icount = 1440
tsrgpars(1)%time_format = 3
tsrgpars(1)%tlims = (/0,nstep,icount/)

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_tsr_params

!========================================================================

SUBROUTINE usrdef_tsr0d_vals(out0ddat,n0vars)
!************************************************************************
!
! *usrdef_tsr0d_vals* 0-D output data
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description - test case drythin
!
! Reference -
!
! Calling program - time_series
!
! Module calls - combine_mod
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

END SUBROUTINE usrdef_tsr0d_vals

!========================================================================

SUBROUTINE usrdef_tsr2d_vals(out2ddat,i,j,n2vars)
!************************************************************************
!
! *usrdef_tsr2d_vals* 2-D output data
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description - test case drythin
!
! Reference -
!
! Calling program - time_series
!
! Module calls - Uvar_at_C
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

END SUBROUTINE usrdef_tsr2d_vals

!========================================================================

SUBROUTINE usrdef_tsr3d_vals(out3ddat,i,j,k,n3vars)
!************************************************************************
!
! *usrdef_tsr3d_vals* 3-D output data
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description - test case drythin
!
! Reference -
!
! Calling program - time_series
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
!*n3vars*   INTEGER number of 2-D variables
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_tsr3d_vals
