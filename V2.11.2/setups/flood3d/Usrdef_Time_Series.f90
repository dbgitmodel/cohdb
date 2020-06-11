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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.7.1
!
! $Date: 2016-10-10 16:36:15 +0200 (Mon, 10 Oct 2016) $
!
! $Revision: 973 $
!
! Description - test case flood3d
!             - flooding/drying channel experiments
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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.7.1
!
! Description - test case flood3d
!             - flooding/drying channel experiments
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
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: icount

procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!
!---fortran name
tsrvars(1)%f90_name = 'depmeanerr'
tsrvars(2)%f90_name = 'depmaxerr'
tsrvars(9)%f90_name = 'mask'

!---standard name
tsrvars(1)%standard_name = 'mean_depth_error' 
tsrvars(2)%standard_name = 'maximum_depth_error' 
tsrvars(9)%standard_name = 'sea_area' 
tsrvars(1:2)%comment = &
            & 'standard name constructed following the netCDF CF-1.6 guidelines'
tsrvars(9)%comment = &
            & 'standard name constructed following the netCDF CF-1.6 guidelines'

!---long name
tsrvars(1)%long_name = 'Mean depth error' 
tsrvars(2)%long_name = 'Maximum depth error' 
tsrvars(9)%long_name = 'Sea area' 

!---variable attributes
tsrvars(1)%units = 'm'
tsrvars(2)%units = 'm'
tsrvars(9)%units = '_'

!---varids and ranks
tsrvars(1:9)%ivarid = (/0,0,iarr_dryfac,iarr_zeta,iarr_umvel,iarr_vmvel,&
                      & iarr_deptotatc,iarr_deptotatc_err,0/)
IF (iopt_grid_nodim.EQ.3) THEN
   tsrvars(10:12)%ivarid = (/iarr_uvel,iarr_vvel,iarr_wphys/)
ENDIF
tsrvars(1:3)%nrank = 0
tsrvars(4:9)%nrank = 2
IF (iopt_grid_nodim.EQ.3) tsrvars(10:12)%nrank = 3

!
!2. Variable indices
!-------------------
!

IF (iopt_verif.EQ.0) THEN
   ivarstsr(1,1:6) = (/4,5,6,7,8,9/)
   IF (iopt_grid_nodim.EQ.3) ivarstsr(1,7:9) = (/10,11,12/)
   ivarstsr(2,1:3) = (/1,2,3/)
ELSE
   ivarstsr(1,1:9) = (/1,2,3,4,5,6,7,8,9/)
   IF (iopt_grid_nodim.EQ.3) ivarstsr(1,10:12) = (/10,11,12/)
ENDIF

!
!3. File parameters
!------------------
!

IF (iopt_verif.EQ.0) THEN
   tsr0d(2)%defined = .TRUE.
ENDIF
tsr2d(1)%defined = .TRUE.
tsr3d(1)%defined = MERGE(.TRUE.,.FALSE.,iopt_grid_nodim.EQ.3)

!
!4. Output grid
!--------------
!
!---first file
tsrgpars(1)%time_format = 3
icount = 600
tsrgpars(1)%tlims = (/0,nstep,icount/)
tsrgpars(1)%time_grid = .TRUE.

!---second file
IF (iopt_verif.EQ.0) THEN
   icount = 20
   tsrgpars(2)%tlims = (/0,nstep,icount/)
   tsrgpars(2)%time_format = 3
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
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.3
!
! Description -  test case flood3d
!
! Reference -
!
! Calling program - time_series
!
! Module calls - max_vars, sum2_vars
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE paral_utilities, ONLY: max_vars, sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

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
!*Local variables
!
INTEGER, DIMENSION(4) :: nhdims = 0

procname(pglev+1) = 'usrdef_tsr0d_vals'
CALL log_timer_in()


CALL sum2_vars(deptotatc_err,out0ddat(1),nhdims,'C  ',0,commall=.TRUE.,&
             & mask=seapoint(1:ncloc,1:nrloc))
out0ddat(1) = out0ddat(1)/REAL(noseaatc)
CALL max_vars(deptotatc_err,out0ddat(2),0,commall=.TRUE.,&
            & mask=seapoint(1:ncloc,1:nrloc))

CALL log_timer_out()   


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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.3
!
! Description - test case flood3d
!
! Reference -
!
! Calling program - time_series
!
!************************************************************************
!
USE grid

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

!---mask
out2ddat(6) = MERGE(nodeatc(i,j)+1.0,0.0,seapoint(i,j))


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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.3
!
! Description - test case flood3d
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
