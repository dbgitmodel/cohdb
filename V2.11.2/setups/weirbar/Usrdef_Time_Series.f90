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
! Description - test case weirbar
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
! Description - test case weirbar
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
CHARACTER (LEN=12) :: cl
INTEGER :: l, numweirs, icount


procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!
!---fortran name, long name, units
numweirs = MERGE(1,7,LLT(runtitle(8:8),'E'))
IF (numweirs.EQ.1) THEN
   tsrvars(1)%f90_name = 'umvel_wb'
   tsrvars(2)%f90_name = 'deptotatu_wb'
   tsrvars(1)%standard_name = &
                 & 'integral_wrt_depth_of_horizontal_sea_water_velocity_at_weir'
   tsrvars(2)%standard_name = 'sea_floor_depth_below_sea_surface_at_weir'
   tsrvars(1:2)%comment = &
            & 'standard name constructed following the netCDF CF-1.6 guidelines'
   tsrvars(1)%long_name = 'Depth-mean current at weir '
   tsrvars(2)%long_name = 'Water depth at weir '
   tsrvars(1)%units = 'm s-1'
   tsrvars(2)%units = 'm'
ELSE
   l_110: DO l=1,numweirs
      WRITE(cl,'(I12)') l; cl = ADJUSTL(cl)
      tsrvars(2*l-1)%f90_name = 'umvel_wb'//TRIM(cl)
      tsrvars(2*l)%f90_name = 'deptotatu_wb'//TRIM(cl)
      tsrvars(2*l-1)%standard_name = &
      & 'integral_wrt_depth_of_horizontal_sea_water_velocity_at_weir_'//TRIM(cl)
      tsrvars(2*l)%standard_name = &
      & 'sea_floor_depth_below_sea_surface_at_weir_'//TRIM(cl)
      tsrvars(2*l-1:2*l)%comment = &
            & 'standard name constructed following the netCDF CF-1.6 guidelines'
      tsrvars(2*l-1)%long_name = 'Depth-mean current at weir '//TRIM(cl)
      tsrvars(2*l)%long_name = 'Water depth at weir '//TRIM(cl)
      tsrvars(2*l-1)%units = 'm s-1'
      tsrvars(2*l)%units = 'm'
   ENDDO l_110
ENDIF

!---varids and ranks
tsrvars(1:2*numweirs)%ivarid = 0
tsrvars(2*numweirs+1:novarstsr)%ivarid = (/iarr_zeta,iarr_umvel,iarr_vmvel,&
                                     & iarr_deptotatc,iarr_hmvelmag,iarr_uvel,&
                                     & iarr_vvel,iarr_wphys,iarr_hvelmag/)
tsrvars(1:2*numweirs)%nrank = 0
tsrvars(2*numweirs+1:novarstsr)%nrank = (/2,2,2,2,2,3,3,3,3/)

!
!2. Variable indices
!-------------------
!

ivarstsr(1,1:9) = (/(l,l=2*numweirs+1,novarstsr)/)
ivarstsr(2,1:2*numweirs) = (/(l,l=1,2*numweirs)/)

!
!3. File parameters
!------------------
!

tsr0d(2)%defined = .TRUE.
tsr2d(1)%defined = .TRUE.
tsr3d(1)%defined = MERGE(.TRUE.,.FALSE.,iopt_grid_nodim.EQ.3)

!
!4. Output grid
!--------------
!
!---first file
icount = 2880
tsrgpars(1)%time_format = 3
tsrgpars(1)%time_grid = .TRUE.
tsrgpars(1)%tlims = (/0,nstep,icount/)

!---second file
icount = 240
tsrgpars(2)%time_format = 3
tsrgpars(2)%tlims = (/0,nstep,icount/)

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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.6
!
! Description - test case weirbar
!
! Reference -
!
! Calling program - time_series
!
!************************************************************************
!
USE currents
USE depths
USE iopars
use timepars
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
INTEGER :: i, l


procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

IF (LLT(runtitle(8:8),'E')) THEN
   out0ddat(1) = umvel(101,1)
   out0ddat(2) = deptotatu(101,1)
ELSE
   l_110: DO l=1,7
      i = 71 + 10*(l-1)
      out0ddat(2*l-1)= 0.2*SUM(umvel(i,1:5))
      out0ddat(2*l)= 0.2*SUM(deptotatu(i,1:5))
   ENDDO l_110
ENDIF

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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description - test case weirbar
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
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description - test case weirbar
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
