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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.11.2
!
! $Date: 2018-07-23 16:55:25 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1171 $
!
! Description - test case seamount
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
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.11.2
!
! Description - test case seamount
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
tsrvars(1)%f90_name = 'ekin'
tsrvars(2)%f90_name = 'ekin2d'
tsrvars(3)%f90_name = 'ekinbc'
tsrvars(4)%f90_name = 'ekinratio'
tsrvars(5)%f90_name = 'edens'
tsrvars(6)%f90_name = 'velmean'
tsrvars(7)%f90_name = 'velmean2d'
tsrvars(8)%f90_name = 'velmeandev'
tsrvars(9)%f90_name = 'velmax'
tsrvars(10)%f90_name = 'velmax2d'
tsrvars(11)%f90_name = 'velmaxdev'


!---standard name
tsrvars(1)%standard_name = 'kinetic_energy_in_sea_water'
tsrvars(2)%standard_name = 'barotropic_kinetic_energy_in_sea_water'
tsrvars(3)%standard_name = 'baroclinic_kinetic_energy_in_sea_water'
tsrvars(4)%standard_name = 'ratio_of_baroclonic_to_total_kinetic_energy'
tsrvars(5)%standard_name = 'available_potential_energy_in_sea_water'
tsrvars(6)%standard_name = 'mean_magnitude_of_horizontal_sea_water_velocity'
tsrvars(7)%standard_name = 'mean_magnitude_of_depth_mean_horizontal_sea_'//&
                         & 'water_velocity'
tsrvars(8)%standard_name = 'mean__magnitude_of_baroclinic_horizontal_sea_'//&
                           'water_velocity'
tsrvars(9)%standard_name = 'maximum_magnitude_of_horizontal_sea_water_velocity'
tsrvars(10)%standard_name = 'maximum_magnitude_of_depth_mean_horizontal_sea_'//&
                          & 'water_velocity'
tsrvars(11)%standard_name = 'maximum_magnitude_of_baroclinic_horizontal_sea_'//&
                            'water_velocity'

!---long name
tsrvars(1)%long_name = 'Kinetic energy'
tsrvars(2)%long_name = 'Barotropic kinetic energy'
tsrvars(3)%long_name = 'Baroclinic kinetic energy'
tsrvars(4)%long_name = 'Ratio of baroclinic to barotropic kinetic energy'
tsrvars(5)%long_name = 'Baroclinic potential energy'
tsrvars(6)%long_name = 'Mean current magnitude'
tsrvars(7)%long_name = 'Mean depth-mean current magnitude'
tsrvars(8)%long_name = 'Mean baroclinic current magnitude'
tsrvars(9)%long_name = 'Maximum current magnitude'
tsrvars(10)%long_name = 'Maximum depth-mean current magnitude'
tsrvars(11)%long_name = 'Maximum baroclinic current magnitude'

!---units
tsrvars(1)%units = 'PJ'
tsrvars(2)%units = 'PJ'
tsrvars(3)%units = 'PJ'
tsrvars(4)%units = '1'
tsrvars(5)%units = 'PJ'
tsrvars(6)%units = 'cm s-1'
tsrvars(7)%units = 'cm s-1'
tsrvars(8)%units = 'cm s-1'
tsrvars(9)%units = 'cm s-1'
tsrvars(10)%units = 'cm s-1'
tsrvars(11)%units = 'cm s-1'

!---varids and ranks
tsrvars%ivarid = (/0,0,0,0,0,0,0,0,0,0,0,iarr_umvel,iarr_vmvel,iarr_zeta,&
                 & iarr_uvel,iarr_vvel,iarr_wphys,iarr_temp/)
tsrvars%nrank = (/0,0,0,0,0,0,0,0,0,0,0,2,2,2,3,3,3,3/)

!
!2. Variable indices
!-------------------
!

ivarstsr(1,1:7) = (/12,13,14,15,16,17,18/)
ivarstsr(2,1:11) = (/1,2,3,4,5,6,7,8,9,10,11/)

!
!3. File parameters
!------------------
!

tsr2d(1)%defined = .TRUE.
tsr3d(1)%defined = .TRUE.
tsr0d(2)%defined = .TRUE.

!
!4. Output grid
!--------------
!

icount = NINT(43200/delt2d)
tsrgpars(1)%tlims = (/0,nstep,icount/)
icount = NINT(3600/delt2d)
tsrgpars(2)%tlims = (/0,nstep,icount/)
tsrgpars(1)%time_format = 4
tsrgpars(2)%time_format = 3
tsrgpars%vcoord = 1

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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.11.2
!
! Description - test case seamount
!
! Reference -
!
! Calling program - time_series
!
! Module calls - energy_0d, max_vars, sum2_vars, vector_mag_arr_atc
!
!************************************************************************
!
USE currents  
USE grid
USE gridpars
USE iopars
USE modids
USE diagnostic_routines, ONLY: energy_0d
USE math_library, ONLY: vector_mag_arr_atc
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
INTEGER :: k
INTEGER, DIMENSION(4) :: nhdims
REAL, DIMENSION(5) :: ecomps
REAL, DIMENSION(ncloc+1,nrloc,nz) :: array3du
REAL, DIMENSION(ncloc,nrloc+1,nz) :: array3dv
REAL, DIMENSION(ncloc,nrloc) :: curmag2d
REAL, DIMENSION(ncloc,nrloc,nz) :: curmag


procname(pglev+1) = 'usrdef_tsr0d_vals'
CALL log_timer_in()

nhdims = 0

!---energy components
CALL energy_0d(3,ecomps)
out0ddat(1) = 1.0E-09*ecomps(1)
out0ddat(5) = 1.0E-09*ecomps(4)
CALL energy_0d(2,ecomps)
out0ddat(2) = 1.0E-09*ecomps(1)
out0ddat(3) = out0ddat(1) - out0ddat(2)
IF (out0ddat(2).GT.0.0) THEN
   out0ddat(4) = out0ddat(3)/out0ddat(2)
ELSE
   out0ddat(4) = 0.0
ENDIF

!---3-D current mean and maximum
CALL vector_mag_arr_atc(uvel(1:ncloc+1,1:nrloc,:),vvel(1:ncloc,1:nrloc+1,:),&
                      & 2,1,nz,1,iarr_hvelmag,.TRUE.,vecmag=curmag)
CALL max_vars(curmag,out0ddat(9),iarr_hvelmag,mask=maskatc_int)
CALL sum2_vars(curmag,out0ddat(6),nhdims,'C  ',0)
out0ddat(6) = out0ddat(6)/(noseaatc*nz)

!---3-D current mean and maximum
CALL vector_mag_arr_atc(umvel(1:ncloc+1,1:nrloc),vmvel(1:ncloc,1:nrloc+1),&
                      & 2,1,1,1,iarr_hmvelmag,.TRUE.,vecmag=curmag2d)
CALL max_vars(curmag2d,out0ddat(10),iarr_hvelmag,mask=maskatc_int)
CALL sum2_vars(curmag2d,out0ddat(7),nhdims,'C  ',0)
out0ddat(7) = out0ddat(7)/noseaatc

!---3-D baroclinic current mean and maximum
k_110: DO k=1,nz
   array3du(:,:,k) = uvel(1:ncloc+1,1:nrloc,k)- umvel(1:ncloc+1,1:nrloc)
   array3dv(:,:,k) = vvel(1:ncloc,1:nrloc+1,k)- vmvel(1:ncloc,1:nrloc+1)
ENDDO k_110
CALL vector_mag_arr_atc(array3du,array3dv,2,1,nz,1,0,.TRUE.,vecmag=curmag)
CALL max_vars(curmag,out0ddat(11),0,mask=maskatc_int)
CALL sum2_vars(curmag,out0ddat(8),nhdims,'C  ',0)
out0ddat(8) = out0ddat(8)/(noseaatc*nz)

out0ddat(6:11) = 100.0*out0ddat(6:11)

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
