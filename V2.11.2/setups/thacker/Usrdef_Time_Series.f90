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
! Author - IMDC
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.5
!
! $Date: 2014-09-22 09:37:46 +0200 (Mon, 22 Sep 2014) $
!
! $Revision: 742 $
!
! Description - test case thacker
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
! Author - IMDC
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.5
!
! Description - test case thacker
!
! Reference -
!
! Calling program - time_series_init
!
!************************************************************************
!
USE gridpars
USE iopars
USE modids
USE sedids
USE switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i

procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

!
!1. Variable attributes
!----------------------
!
!---varids
tsrvars%ivarid = (/iarr_zeta,iarr_umvel,iarr_vmvel,iarr_deptotatc,&
                 & iarr_bstresatc,iarr_bottom_sed_flux,iarr_ctot,iarr_uvel,&
                 & iarr_vvel,iarr_wphys,iarr_ctot/)
!---rank
tsrvars(1:7)%nrank = 2
tsrvars(8:11)%nrank = 3

!---operator
tsrvars(7)%oopt = oopt_int

!---variable number
tsrvars(6)%numvar = 1

!
!2. Variable indices
!-------------------
!

ivarstsr(1,1:11) = (/(i,i=1,11)/)

!
!3. File parameters
!------------------
!


tsr2d(1)%defined = .TRUE.
tsr3d(1)%defined = iopt_grid_nodim.EQ.3

!
!4. Output grid
!--------------
!

tsrgpars(1)%time_format = 2
tsrgpars(1)%zlims = (/1,nz,1/)
tsrgpars(1)%tlims = (/0,nstep,60/)

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_tsr_params

!========================================================================

SUBROUTINE usrdef_tsr0d_vals(out0ddat,n0vars)
!************************************************************************
!
! *usrdef_tsr0d_vals* 0-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.0
!
! Description -
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


out0ddat = 0.0


RETURN

END SUBROUTINE usrdef_tsr0d_vals

!========================================================================

SUBROUTINE usrdef_tsr2d_vals(out2ddat,i,j,n2vars)
!************************************************************************
!
! *usrdef_tsr2d_vals* 2-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description -
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


out2ddat = 0.0


RETURN

END SUBROUTINE usrdef_tsr2d_vals

!========================================================================

SUBROUTINE usrdef_tsr3d_vals(out3ddat,i,j,k,n3vars)
!************************************************************************
!
! *usrdef_tsr3d_vals* 3-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description -
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


out3ddat = 0.0


RETURN

END SUBROUTINE usrdef_tsr3d_vals

