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
!

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
! Description - test case discharges 
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
! Description - test case discharges 
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
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out
USE timepars 

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
!---fortran name
tsrvars(1)%f90_name = 'dvolume0d'
tsrvars(2)%f90_name = 'balance0d'

!---standard name
tsrvars(1)%standard_name = 'volume_increase'
tsrvars(2)%standard_name = 'volume_balance'
tsrvars(1:2)%comment = &
            & 'standard name constructed following the netCDF CF-1.6 guidelines'

!---long name
tsrvars(1)%long_name = 'Volume increase'
tsrvars(2)%long_name = 'Volume balance'

!---units
tsrvars(1)%units = 'percent'
tsrvars(2)%units = 'percent'

!---varids and ranks
tsrvars%ivarid = (/0,0,iarr_zeta,iarr_umvel,iarr_vmvel,iarr_hmvelmag,&
                 & iarr_uvel,iarr_vvel,iarr_wphys,iarr_hvelmag,iarr_sal/)
tsrvars%nrank = (/0,0,2,2,2,2,3,3,3,3,3/)

!
!2. Variable indices
!-------------------
!

ivarstsr(1,1:9) = (/3,4,5,6,7,8,9,10,11/)
ivarstsr(2,1:2) = (/1,2/)

!
!3. File parameters
!------------------
!

tsr0d(2)%defined = .TRUE.
tsr2d(1)%defined = .TRUE.
tsr3d(1)%defined = .TRUE.

!
!4. Output grid
!--------------
!
!---first file
icount = 6
tsrgpars(1)%time_format = 2
tsrgpars(1)%tlims = (/0,nstep,icount/)     
tsrgpars(1)%time_grid = .TRUE.

!---second file 
icount = 2
tsrgpars(2)%time_format = 2
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
! Description - test case discharges 
!
! Reference -
!
! Calling program - time_series
!
! Module calls - sum2_vars
!
!************************************************************************
!
USE depths
USE gridpars
USE iopars
USE modids
USE physpars
USE structures
USE timepars
USE paral_utilities, ONLY: sum2_vars
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
INTEGER, SAVE, DIMENSION(4) :: nhdims
REAL, SAVE :: delxdat, delydat, distot, volume_ini 
REAL :: dvolume


CALL log_timer_in()
procname(pglev+1) = 'usrdef_out0d_vals'

!
!1. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN
   delxdat = surfacegrids(igrd_model,1)%delxdat
   delydat = surfacegrids(igrd_model,1)%delydat
   nhdims = 1
   volume_ini = (nc-1)*(nr-1)*depmean_cst
   distot = 0.0
ENDIF

!
!2. Increase of total volume
!---------------------------
!

CALL sum2_vars(zeta,dvolume,nhdims,'C  ',iarr_zeta,commall=.TRUE.)

!
!3. Time integrated volume discharge
!-----------------------------------
!

IF (nt.GT.0.0) THEN
   distot = distot + 2.0*delt2d*SUM(disvol)/(delxdat*delydat)
ENDIF

!
!4. Volume balance
!-----------------
!

out0ddat(1) = 100.0*dvolume/volume_ini
out0ddat(2) = 100.0*(dvolume-distot)/volume_ini

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
! Author - ANTEA
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

