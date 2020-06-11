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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.8
!
! $Date: 2018-07-23 16:34:19 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1170 $
!
! Description - test case plume
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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.8
!
! Description - test case plume
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
tsrvars(1)%f90_name = 'width'
tsrvars(2)%f90_name = 'front'

!---standard names
tsrvars(1)%standard_name = 'plume_width'
tsrvars(2)%standard_name = 'plume_length'
tsrvars(1:2)%comment = &
           &  'standard name constructed following the netCDF CF-1.6 guidelines'

!---long names
tsrvars(1)%long_name = 'Plume width'
tsrvars(2)%long_name = 'Plume length'

!---units
tsrvars(1)%units = 'km'
tsrvars(2)%units = 'km'

!---varids and ranks
IF (iopt_verif.EQ.0) THEN
    tsrvars%ivarid = (/0,0,iarr_umvel,iarr_vmvel,iarr_zeta,iarr_uvel,iarr_vvel,&
                     & iarr_wphys,iarr_sal/)
    tsrvars%nrank = (/0,0,2,2,2,3,3,3,3/)
ELSE
   tsrvars%ivarid = (/0,0,iarr_umvel,iarr_vmvel,iarr_zeta,&
                    & iarr_hdifcoef2d_scal,iarr_uvel,iarr_vvel,iarr_wphys,&
                    & iarr_sal/)
    tsrvars%nrank = (/0,0,2,2,2,2,3,3,3,3/)
ENDIF

!
!2. Variable indices
!-------------------
!

IF (iopt_verif.EQ.0) THEN
    ivarstsr(1,1:4) = (/6,7,8,9/)
    ivarstsr(2,1:4) = (/6,7,8,9/)
    ivarstsr(3,1:4) = (/6,7,8,9/)
    ivarstsr(4,1:5) = (/1,2,3,4,5/)
ELSE
    ivarstsr(1,:) = (/1,2,3,4,5,6,7,8,9,10/)
ENDIF

!
!3. File parameters
!------------------
!

IF (iopt_verif.EQ.0) THEN
    tsr3d(1:3)%defined = .TRUE.
    tsr0d(4)%defined = .TRUE.
    tsr2d(4)%defined = .TRUE.
ELSE
    tsr0d(1)%defined = .TRUE.
    tsr2d(1)%defined = .TRUE.
    tsr3d(1)%defined = .TRUE.
ENDIF

!
!4. Output grid
!--------------
!

IF (iopt_verif.EQ.0) THEN
   icount = MERGE(360,36,iopt_hydro_impl.EQ.0)
!  ---first file
   tsrgpars(1)%zlims = (/nz,nz,1/)
   tsrgpars(1)%tlims = (/0,nstep,icount/)
!  ---second file
   tsrgpars(2)%xlims = (/30,30,1/)
   tsrgpars(2)%tlims = (/0,nstep,icount/)
!  ---third file
   tsrgpars(3)%ylims = (/5,5,1/)
   tsrgpars(3)%tlims = (/0,nstep,icount/)
!  ---fourth file
   icount = MERGE(10,1,iopt_hydro_impl.EQ.0)
   tsrgpars(4)%xlims = (/30,30,1/)
   tsrgpars(4)%ylims = 1
   tsrgpars(4)%zlims = 1
   tsrgpars(4)%tlims = (/0,nstep,icount/)
ELSE
   icount = MERGE(360,36,iopt_hydro_impl.EQ.0)
   tsrgpars(1)%tlims = (/0,nstep,icount/)
ENDIF

tsrgpars(1:nosetstsr)%time_format = 3

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
! Description - test case plume
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
USE grid
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
INTEGER :: i, ib, ic, j, jc
REAL, DIMENSION(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo,nz) :: salglb
REAL :: dhun, gxc, gyc, hfront, scrit, width


CALL log_timer_in()
procname(pglev+1) = 'usrdef_out0d_vals'

!
!1. Evaluate parameters
!----------------------
!

CALL combine_mod(salglb,sal,(/1-nhalo,1-nhalo,1/),iarr_sal,0.0)

!
!1.1 Plume width
!---------------
!

IF (master) THEN
   dhun = 1000.0; scrit = 32.0
   ic = iobv(nobv)
   jc = 0
   j_110: DO j=1,nr-1
      IF ((jc.EQ.0).AND.(salglb(ic,j,nz).GT.scrit)) jc = j
   ENDDO j_110
   IF (jc.EQ.0) then
      width = dhun*(nr-1)
   ELSEIF (jc.EQ.1) THEN
      width = 0.0
   ELSE
      gyc = dhun*(jc-1.5)
      width = gyc+(scrit-salglb(ic,jc-1,nz))*dhun/&
                & (salglb(ic,jc,nz)-salglb(ic,jc-1,nz))
   ENDIF
   width = 0.001*width
ENDIF

!
!1.2 Plume length
!----------------
!

IF (master) THEN
   scrit = 32.0
   jc = 2
   ic = 0
   ib = iobv(nobv)
   i_120: DO i=ib,nc-1
      IF ((ic.EQ.0).AND.(salglb(i,jc,nz).GT.scrit)) ic = i
   ENDDO i_120
   IF (ic.EQ.0) THEN
      hfront = dhun*(nc-0.5-ib)
   ELSEIF (ic.EQ.ib) THEN
      hfront = 0.0
   ELSE
      gxc = dhun*(ic-ib-1)
      hfront = gxc+(scrit-salglb(ic-1,jc,nz))*dhun/&
                 & (salglb(ic,jc,nz)-salglb(ic-1,jc,nz))
   ENDIF
   hfront = 0.001*hfront
ENDIF

!
!2. Store parameters
!-------------------
!

IF (master) THEN
   out0ddat(1) = width
   out0ddat(2) = hfront
ELSE
   out0ddat(1:2) = 0.0
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
! Description - test case plume
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
! Description - test case plume
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
REAL, DIMENSION(n3vars) :: out3ddat

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
