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
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.11
!
! $Date: 2016-11-25 12:04:00 +0100 (Fri, 25 Nov 2016) $
!
! $Revision: 994 $
!
! Description - test case partrhone
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
! *usrdef_out_params* Specifiers for time series output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.11
!
! Description - test case partrhone
!
! Reference -
!
! Calling program - time_series_init
!
!************************************************************************
!
USE iopars
USE modids  
USE partids
USE switches
USE timepars
use modvars_routines, ONLY: inquire_var
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
INTEGER :: icount, l

procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

flag = runtitle(10:10).EQ.'A'

!
!1. Output variables
!-------------------
!
!---varids and ranks
IF (flag) THEN
   tsrvars%ivarid = (/iarr_part_conc,iarr_part_conc,iarr_part_conc,&
                    & iarr_part_conc_tot,iarr_uvel,iarr_vvel,iarr_wphys,&
                    & iarr_sal/)
ELSE
   tsrvars%ivarid = (/iarr_part_conc,iarr_part_conc,iarr_part_conc,&
                    & iarr_part_conc_tot/)
ENDIF

tsrvars%nrank = 3

!---variable number
tsrvars(1:3)%numvar = (/1,2,3/)

!
!2. Variable indices
!-------------------
!

IF (flag) THEN
   ivarstsr(1,1:4) = (/1,2,3,4/)
   ivarstsr(2,1:4) = (/5,6,7,8/)
ELSE
   ivarstsr(1,1:novarstsr) = (/(l,l=1,novarstsr)/)
ENDIF

!
!3. File parameters
!------------------
!

tsr3d(1)%defined = .TRUE.
IF (flag) tsr3d(2)%defined = .TRUE.

!
!4. Output grid
!--------------
!
!---first file
tsrgpars(1)%time_format = 3
tsrgpars(1)%tlims(1) = MERGE(12960,720,iopt_hydro_impl.EQ.0)
tsrgpars(1)%tlims(2) = nstep
tsrgpars(1)%tlims(3) = MERGE(720,40,iopt_hydro_impl.EQ.0)

!---second file file
IF (flag) THEN
   tsrgpars(2)%time_format = 3
   icount = MERGE(4320,240,iopt_hydro_impl.EQ.0)
   tsrgpars(2)%tlims = (/icount,nstep,icount/)
ENDIF

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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description -
!
! Reference -
!
! Calling program - time_series
!
! Module calls - combine_mod, define_out0d_vals, max_vars, sum2_vars
!
!************************************************************************
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


RETURN

END SUBROUTINE usrdef_tsr3d_vals
