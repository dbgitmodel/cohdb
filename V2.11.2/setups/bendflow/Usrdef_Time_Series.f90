! This file is part of COHERENS. You can redistribute and/or modify it under
! the conditions of the COHERENS license. Details are found in the file
! COHERENS_License.

!************************************************************************
!
! *Usrdef_Time_Series* Define time series output
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.5
!
! $Date: 2015-10-06 17:42:25 +0200 (Tue, 06 Oct 2015) $
!
! $Revision: 886 $
!
! Description - test case bendflow
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
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.5
!
! Description - test case bendflow
!
! Reference -
!
! Calling program - time_series_init
!
!************************************************************************
!
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
IF (iopt_grid_nodim.EQ.2) THEN
   tsrvars%ivarid = (/iarr_zeta,iarr_umvel,iarr_vmvel,iarr_bstresatc,&
                    & iarr_qbedatu,iarr_qbedatv,iarr_bstres_cr,iarr_cref,&
                    & iarr_ceq,iarr_t_equil,iarr_ctot/)
ELSE
   tsrvars%ivarid = (/iarr_zeta,iarr_umvel,iarr_vmvel,iarr_bstresatc,&
                    & iarr_qbedatu,iarr_qbedatv,iarr_bstres_cr,iarr_cref,&
                    & iarr_uvel,iarr_vvel,iarr_wphys,iarr_ctot,iarr_wfall,&
                    & iarr_vdiffcoef_sed/)
ENDIF

!---rank
IF (iopt_grid_nodim.EQ.2) THEN
   tsrvars(1:10)%nrank = 2; tsrvars(11)%nrank = 3
ELSE
   tsrvars(1:8)%nrank = 2; tsrvars(9:14)%nrank = 3
ENDIF

!---operators
IF (iopt_grid_nodim.EQ.2) THEN
   tsrvars(11)%oopt = oopt_klev
   tsrvars(11)%klev = 1
ENDIF

!---variable number
IF (iopt_grid_nodim.EQ.2) THEN
   tsrvars(5:10)%numvar = 1
ELSE
   tsrvars(5:8)%numvar = 1
   tsrvars(13:14)%numvar = 1
ENDIF

!
!2. Variable indices
!-------------------
!

ivarstsr(1,1:novarstsr) = (/(i,i=1,novarstsr)/)

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

tsrgpars(1)%tlims = (/0,nstep,750/)

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
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.1.2
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
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.1.2
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
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.1.2
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

