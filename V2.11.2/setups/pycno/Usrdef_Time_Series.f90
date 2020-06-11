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
! Description - test case pycno
!
! Reference -
!
! Routines - usrdef_tsr_params, usrdef_tsr0d_vals, usrdef_tsr2d_vals
!             usrdef_tsr3d_vals, mixed_layer
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
! Description - test case pycno
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
USE syspars
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=lendesc) :: comment


procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!

IF (iopt_verif.EQ.0) THEN

!   ---fortran name
    tsrvars(1)%f90_name = 'tmld0d'
    tsrvars(2)%f90_name = 'usur'
    tsrvars(3)%f90_name = 'dr0'
    tsrvars(4)%f90_name = 'veddyd'
    tsrvars(6)%f90_name = 'drho'
    tsrvars(11)%f90_name = 'rich'

!   ---standard name
    tsrvars(1)%standard_name = 'mixed_layer_depth'
    tsrvars(2)%standard_name = 'surface_x_sea_water_velocity'
    tsrvars(3)%standard_name = 'surface_density_difference'
    tsrvars(4)%standard_name = 'sea_water_turbulent_diffusivity'
    tsrvars(6)%standard_name = 'density_difference'
    tsrvars(11)%standard_name = 'richardson_number_in_sea_water'

!   ---comment
    comment = 'standard name constructed following the netCDF CF-1.6 guidelines'
    tsrvars(1:4)%comment = comment
    tsrvars(6)%comment = comment
    tsrvars(11)%comment = comment

!   ---long name
    tsrvars(1)%long_name = 'Mixed layer depth'
    tsrvars(2)%long_name = 'Surface current'
    tsrvars(3)%long_name = 'Surface density difference'
    tsrvars(4)%long_name = 'Eddy diffusivity'
    tsrvars(6)%long_name = 'Density difference'
    tsrvars(11)%long_name = 'Richardson number'
    
!   ---unit
    tsrvars(1)%units = 'm'
    tsrvars(2)%units = 'm s-1'
    tsrvars(3)%units = 'kg m3'
    tsrvars(4)%units = 'm2 s-1'
    tsrvars(6)%units = 'kg m3'
    tsrvars(11)%units = ''

ELSE

!   ---fortran name
    tsrvars(1)%f90_name = 'tmld0d'
    tsrvars(2)%f90_name = 'dr0'
    tsrvars(3)%f90_name = 'veddyd'
    tsrvars(5)%f90_name = 'drho'
    tsrvars(6)%f90_name = 'rich'

!   ---standard name
    tsrvars(1)%standard_name = 'mixed_layer_depth'
    tsrvars(2)%standard_name = 'surface_density_difference'
    tsrvars(3)%standard_name = 'sea_water_turbulent_diffusivity'
    tsrvars(5)%standard_name = 'density_difference'
    tsrvars(6)%standard_name = 'richardson_number_in_sea_water'

!   ---comment
    comment = 'standard name constructed following the netCDF CF-1.6 guidelines'
    tsrvars(1:3)%comment = comment
    tsrvars(5:6)%comment = comment

!   ---long name
    tsrvars(1)%long_name = 'Mixed layer depth'
    tsrvars(2)%long_name = 'Surface density difference'
    tsrvars(3)%long_name = 'Eddy diffusivity'
    tsrvars(5)%long_name = 'Density difference'
    tsrvars(6)%long_name = 'Domain-averaged Richardson number'

!   ---unit
    tsrvars(1)%units = 'm'
    tsrvars(2)%units = 'kg m3'
    tsrvars(3)%units = 'm2 s-1'
    tsrvars(5)%units = 'kg m-3'
    tsrvars(6)%units = ''

ENDIF

!---vertical node
IF (iopt_verif.EQ.0) THEN
   tsrvars(7:11)%node = 'W'
ELSE
   tsrvars(6:7)%node = 'W'; tsrvars(12:15)%node = 'W'
ENDIF

!---varids and ranks
IF (iopt_verif.EQ.0) THEN
    tsrvars%ivarid = (/0,0,0,0,iarr_uvel,iarr_dens,iarr_vdifcoefscal,&
                     & iarr_tke,iarr_zlmix,iarr_dissip,iarr_ricnum/)
    tsrvars%nrank = (/0,0,0,0,3,3,3,3,3,3,3/)
ELSE
    tsrvars%ivarid = (/0,0,0,iarr_zeta,iarr_dens,iarr_ricnum,&
                     & iarr_vdifcoefmom,iarr_uvel,iarr_vvel,iarr_wphys,&
                     & iarr_sal,iarr_vdifcoefscal,iarr_tke,iarr_zlmix,&
                     & iarr_dissip/) 
    tsrvars%nrank = (/0,0,0,2,3,3,3,3,3,3,3,3,3,3,3/)
ENDIF

!
!2. Variable indices
!-------------------
!

IF (iopt_verif.EQ.0) THEN
    ivarstsr(1,1:4) = (/1,2,3,4/)
    ivarstsr(2,1:7) = (/5,6,7,8,9,10,11/)
ELSE
    ivarstsr(1,:) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
ENDIF

!
!3. File parameters
!------------------
!

tsr0d(1)%defined = .TRUE.

IF (iopt_verif.EQ.0) THEN
    tsr3d(2)%defined = .TRUE.
ELSE
    tsr3d(1)%defined = .TRUE.
ENDIF

IF (iopt_verif.EQ.1) tsr2d(1)%defined = .TRUE.

!
!4. Output grid
!--------------
!

IF (iopt_verif.EQ.0) THEN
!  ---first file
   tsrgpars(1)%time_format = 3
   tsrgpars(1)%tlims = (/0,nstep,1/)
!  ---second file
   tsrgpars(2)%time_format = 3
   tsrgpars(2)%xlims = 1
   tsrgpars(2)%ylims = 1
   tsrgpars(2)%tlims = (/nstep,nstep,1/)
ELSE
!  ---first file
   tsrgpars(1)%time_format = 1
   tsrgpars(1)%xlims = 1
   tsrgpars(1)%ylims = 1
   tsrgpars(1)%tlims = (/0,nstep,60/)
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
! Description - test case pycno
!
! Reference -
!
! Calling program - time_series
!
! External calls - mixed_layer
!
!************************************************************************
!
USE currents
USE density
USE diffusion
USE gridpars
USE iopars
USE physpars
USE timepars
USE switches
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
REAL :: riv, salin = 0, tmld


procname(pglev+1) = 'usrdef_tsr0d_vals'
CALL log_timer_in()

!---store initial surface salinity on first call
IF (nt.EQ.0) salin = sal(1,1,nz)

!---mixed layer depth
CALL mixed_layer(tmld,riv)

!---output fields
IF (iopt_verif.EQ.0) THEN
    out0ddat(1) = tmld
    out0ddat(2) = uvel(1,1,nz)
    out0ddat(3) = density_ref*beta_sal_ref*(sal(1,1,nz)-salin)
    out0ddat(4) = vdifcoefscal(1,1,180)
ELSE
    out0ddat(1) = tmld
    out0ddat(2) = density_ref*beta_sal_ref*(sal(1,1,nz)-salin)
    out0ddat(3) = vdifcoefscal(1,1,180)
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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.0
!
! Description - test case pycno
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
! Description - test case pycno
!
! Reference -
!
! Calling program - time_series
!
! Module calls - richardson number
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

SUBROUTINE mixed_layer(tmld,riv)
!************************************************************************
!
! *mixed_layer* Mixed layer parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.0
!
! Description - test case pycno
!
! Reference -
!
! Calling program - usrdef_tsr0d_vals, usrdef_output
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
REAL, INTENT(INOUT) :: riv, tmld

!
! Name    Type   Purpose
!------------------------------------------------------------------------------
!*riv*    REAL   Bulk Richardson number
!*tmld*   REAL   Mixed layer depth                                     [m]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
INTEGER :: k
REAL :: bf2 = 1.0E-04, udtur


procname(pglev+1) = 'mixed_layer'
CALL log_timer_in()

flag = .FALSE.
udtur = 0.0
k_201: DO k=nz,1,-1
   IF (.NOT.flag) THEN
      IF (uvel(1,1,k).GT.0.01*uvel(1,1,nz)) THEN
         udtur = udtur + uvel(1,1,k)*delzatc(1,1,k)
      ELSE
         flag = .TRUE.
         tmld = (1.0-gscoordatc(1,1,k))*depmeanatc(1,1)
         udtur = udtur + 0.5*uvel(1,1,k)*delzatc(1,1,k)
      ENDIF
   ENDIF
ENDDO k_201
IF (.NOT.flag) tmld = depmeanatc(1,1)
IF (udtur.GT.0.0001) THEN
   riv = 0.5*bf2*(tmld**2/udtur)**2
ELSE
   riv = 0.0
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE mixed_layer
