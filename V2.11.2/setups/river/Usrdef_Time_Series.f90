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
! $Date: 2016-10-10 16:36:15 +0200 (Mon, 10 Oct 2016) $
!
! $Revision: 973 $
!
! Description - test case river
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
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.7.1
!
! Description - test case river
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
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!
!---name and unit
tsrvars(1)%f90_name = 'hfront0d'
tsrvars(1)%long_name = 'Position of front'
tsrvars(1)%units = 'm'

!---standard name
tsrvars(1)%standard_name = 'distance_of_front_to_western_boundary'
tsrvars(1)%comment = &
            & 'standard name constructed following the netCDF CF-1.6 guidelines'

!---varids and ranks
tsrvars%ivarid = (/0,iarr_zeta,iarr_umvel,iarr_uvel,iarr_vvel,iarr_wphys,&
                 & iarr_sal/)
tsrvars%nrank = (/0,2,2,3,3,3,3/)

!
!2. Variable indices
!-------------------
!

IF (iopt_verif.EQ.0) THEN
    ivarstsr(1,1:6) = (/2,3,4,5,6,7/)
    ivarstsr(2,1) = 1
ELSE
    ivarstsr(1,:) = (/1,2,3,4,5,6,7/)
ENDIF

!
!3. File parameters
!------------------
!

IF (iopt_verif.EQ.0) THEN
    tsr2d(1)%defined = .TRUE.
    tsr3d(1)%defined = .TRUE.
    tsr0d(2)%defined = .TRUE.
ELSE
    tsr0d(1)%defined = .TRUE.
    tsr2d(1)%defined = .TRUE.
    tsr3d(1)%defined = .TRUE.
ENDIF

!
!4. Output grid
!--------------
!
!---first file
tsrgpars(1)%time_format = 3
tsrgpars(1)%tlims(1:2) = (/0,nstep/)
tsrgpars(1)%tlims(3) = 180
tsrgpars(1)%time_grid = .TRUE.
    
!---second file
IF (iopt_verif.EQ.0) THEN
   tsrgpars(2)%time_format = 3
   tsrgpars(2)%tlims(1:2) = (/0,nstep/)
   tsrgpars(2)%tlims(3) = 30
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
! Description - test case river
!
! Reference -
!
! Calling program - time_series
!
! Module calls - combine_mod
!
!************************************************************************
!
USE density
USE gridpars
USE iopars
USE modids
USE paralpars
USE paral_comms, ONLY: combine_mod
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
INTEGER :: i, ic
INTEGER, DIMENSION(4) :: nhdims
REAL, DIMENSION(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo,nz) :: salglb
REAL :: dhun, hfront, scrit, xc


CALL log_timer_in()
procname(pglev+1) = 'usrdef_out0d_vals'

nhdims = nhalo
CALL combine_mod(salglb,sal,(/1-nhalo,1-nhalo,1/),iarr_sal,0.0)

IF (master) THEN
   dhun = 1000.0
   scrit = 34.0
   ic = 0
   i_110: DO i=1,nc-1
      IF ((ic.EQ.0).AND.(salglb(i,1,nz).GT.scrit)) ic = i
   ENDDO i_110
   IF (ic.EQ.0) THEN
      hfront = (nc-1)*dhun
   ELSEIF (ic.EQ.1) THEN
      hfront = 0.0
   ELSE
      xc = dhun*(ic-0.5)
      hfront = xc+(scrit-salglb(ic-1,1,nz))*dhun/&
                & (salglb(ic,1,nz)-salglb(ic-1,1,nz))
   ENDIF
   out0ddat(1) = 0.001*hfront
ELSE
   out0ddat(1) = 0.0
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
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description - test case river
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
! Description - test case river
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
