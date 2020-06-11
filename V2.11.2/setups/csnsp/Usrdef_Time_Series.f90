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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.7.1
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - test case csnsp
!
! Reference -
!
! Routines - usrdef_tsr_params, usrdef_tsr0d_vals, usrdef_tsr2d_vals,
!            usrdef_tsr3d_vals, mixed_layer
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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.7.1
!
! Description - test case csnsp
!
! Reference -
!
! Calling program - time_series_init
!
!************************************************************************
!
USE datatypes
USE gridpars
USE iopars
USE modids
USE switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!

tsrvars(1)%f90_name = 'tmld'
tsrvars(1)%standard_name = 'mixed_layer_depth'
tsrvars(1)%comment = &
            & 'standard name constructed following the netCDF CF-1.6 guidelines'
tsrvars(1)%long_name = 'Mixed layer depth'
tsrvars(1)%units = 'm'
IF(iopt_verif.EQ.0) THEN
   tsrvars%ivarid = (/0,iarr_temp/)
   tsrvars%nrank  = (/0,3/)
ELSE
   tsrvars%ivarid = (/0,iarr_zeta,iarr_temp,iarr_uvel,iarr_vvel,iarr_sal,&
                    & iarr_wphys/)
   tsrvars%nrank  = (/0,2,3,3,3,3,3/)
ENDIF

!
!2. Variable indices
!-------------------
!

IF (iopt_verif.EQ.0) THEN
    ivarstsr(1,1:2) = (/1,2/)
    ivarstsr(2,1) = 2
ELSE
    ivarstsr(1,:) = (/1,2,3,4,5,6,7/)
ENDIF

!
!3. File parameters
!------------------
!

tsr0d(1)%defined = .TRUE.

IF (iopt_verif.EQ.0) THEN
   tsr3d(1:2)%defined = .TRUE.
   tsr3d(2)%form = 'A'
ELSE
   tsr2d(1)%defined = .TRUE.
   tsr3d(1)%defined = .TRUE.
ENDIF

!
!4. Output grid
!--------------
!
!---first file
tsrgpars(1)%time_format = 4
tsrgpars(1)%refdate = '1989/01/01;00:00:00:000'
tsrgpars(1)%xlims = (/2,2,1/)
tsrgpars(1)%ylims = (/2,2,1/)
tsrgpars(1)%tlims(1:2) = (/0,nstep/)
tsrgpars(1)%tlims(3) = MERGE(9,135,iopt_verif.EQ.0)
    
!---second file
IF (iopt_verif.EQ.0) THEN
   tsrgpars(2)%time_format = 4
   tsrgpars(2)%refdate = '1989/01/01;00:00:00:000'
   tsrgpars(2)%xlims = (/2,2,1/)
   tsrgpars(2)%ylims = (/2,2,1/)
   tsrgpars(2)%zlims = (/1,nz,1/)
   tsrgpars(2)%tlims = (/15345,15345,1/)
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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.0
!
! Description - test case csnsp
!
! Reference -
!
! Calling program - time_series
!
! External calls - mixed_layer
!
!************************************************************************
!
USE iopars
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

procname(pglev+1) = 'usrdef_tsr0d_vals'
CALL log_timer_in()

!---mixed layer depth
CALL mixed_layer(out0ddat(1))

CALL log_timer_out()


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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.0
!
! Description - test case csnsp
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
! Description - test case csnsp
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

!========================================================================

SUBROUTINE mixed_layer(thdep)
!************************************************************************
!
! *mixed_layer* Mixed layer depth
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.0
!
! Description - test case csnsp
!
! Reference -
!
! Calling program - usrdef_tsr0d_vals, usrdef_output
!
!************************************************************************
!
USE density
USE depths
USE grid
USE gridpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!

REAL, INTENT(OUT) :: thdep

!
!*Local variables
!
LOGICAL :: flag
INTEGER :: k, kc


procname(pglev+1) = 'mixed_layer'
CALL log_timer_in()

flag = .FALSE.
k_110: DO k=2,nz
   IF ((.NOT.flag).AND.((temp(1,1,k)-temp(1,1,1)).GT.1.0)) THEN
      flag = .TRUE.
      kc = k
   ENDIF
   IF (.NOT.flag) THEN
      thdep  = 0.0
   ELSE
      thdep = depmeanatc(1,1)*(1.0-gscoordatc(1,1,kc-1))-delzatw(1,1,kc)*&
            & (temp(1,1,1)+1.0-temp(1,1,kc-1))/(temp(1,1,kc)-temp(1,1,kc-1))
   ENDIF
ENDDO k_110

CALL log_timer_out()


RETURN

END SUBROUTINE mixed_layer
