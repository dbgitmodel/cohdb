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
! Description - test case fredy
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
! Description - test case fredy
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
!---fortran name
tsrvars(1)%f90_name = 'ekin0d'
tsrvars(2)%f90_name = 'epot0d'
tsrvars(3)%f90_name = 'theta0d'
tsrvars(4)%f90_name = 'enstr0d'
tsrvars(5)%f90_name = 'a1pt0d'
tsrvars(6)%f90_name = 'salmin0d'
tsrvars(7)%f90_name = 'salmax0d'
tsrvars(8)%f90_name = 'smean'
tsrvars(14)%f90_name = 'vort'

!---standard name
tsrvars(1)%standard_name = 'kinetic_energy_in_sea_water'
tsrvars(2)%standard_name = 'available_potential_energy_in_sea_water'
tsrvars(3)%standard_name = &
                 & 'ratio_of_kinetic_to_available_potential_energy_in_sea_water'
tsrvars(4)%standard_name = 'enstrophy_in_sea_water'
tsrvars(5)%standard_name = 'area_bounded_by_surface_34.839_psu_contour'
tsrvars(6)%standard_name = 'salinity_minimum'
tsrvars(7)%standard_name = 'salinity_maximum'
tsrvars(14)%standard_name = 'sea_water_relative_vorticity'

!---comment
comment = 'standard name constructed following the netCDF CF-1.6 guidelines'
tsrvars(1:7)%comment = comment
tsrvars(14)%comment = comment

!---long name
tsrvars(1)%long_name = 'Kinetic energy'
tsrvars(2)%long_name = 'Available potential energy'
tsrvars(3)%long_name = 'Energy ratio'
tsrvars(4)%long_name = 'Enstrophy'
tsrvars(5)%long_name = 'A1%'
tsrvars(6)%long_name = 'Minimum salinity'
tsrvars(7)%long_name = 'Maximum salinity'
tsrvars(8)%long_name = 'Mean salinity deviation'
tsrvars(14)%long_name = 'Vorticity'

!---units
tsrvars(1)%units = 'GJ'
tsrvars(2)%units = 'GJ'
tsrvars(3)%units = '1'
tsrvars(4)%units = 'm3 s-2'
tsrvars(5)%units = '1e8 m2'
tsrvars(6)%units = 'PSU'
tsrvars(7)%units = 'PSU'
tsrvars(8)%units = 'PSU'
tsrvars(14)%units = 's-1'

!---varids and ranks
tsrvars%ivarid = (/0,0,0,0,0,0,0,0,iarr_zeta,iarr_uvel,iarr_vvel,iarr_wphys,&
                 & iarr_sal,0/)
tsrvars%nrank = (/0,0,0,0,0,0,0,0,2,3,3,3,3,3/)

!
!2. Variable indices
!-------------------
!

IF (iopt_verif.EQ.0) THEN
    ivarstsr(1,1:6) = (/9,10,11,12,13,14/)
    ivarstsr(2,1:8) = (/1,2,3,4,5,6,7,8/)
ELSE
    ivarstsr(1,1:14) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14/)
ENDIF

!
!3. File parameters
!------------------
!

IF (iopt_verif.EQ.0) THEN
    tsr0d(2)%defined = .TRUE.
ELSE
    tsr0d(1)%defined = .TRUE.
    tsr2d(1)%defined = .TRUE.
ENDIF
tsr3d(1)%defined = .TRUE.

!
!4. Output grid
!--------------
!
!---first file
tsrgpars(1)%time_format = 4
tsrgpars(1)%tlims(1:2) = (/0,nstep/)
tsrgpars(1)%tlims(3) = 2880

!---second file
IF (iopt_verif.EQ.0) THEN
   tsrgpars(2)%time_format = 3
   tsrgpars(2)%tlims(1:2) = (/0,nstep/)
   tsrgpars(2)%tlims(3) = 20
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
! Description - test case fredy
!
! Reference -
!
! Calling program - time_series
!
! Module calls - define_out0d_vals, sum_vars
!
!************************************************************************
!
USE density
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE model_output, ONLY: define_out0d_vals
USE paral_utilities, ONLY: sum_vars
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
INTEGER :: i1ptloc, i1pt
REAL :: scrit


procname(pglev+1) = 'usrdef_tsr0d_vals'
CALL log_timer_in()

!---energy components
CALL define_out0d_vals(out0ddat(1:2),2,ivarid=(/iarr_ekin0d,iarr_epot0d/))
out0ddat(1:2) = 0.001*out0ddat(1:2)
out0ddat(3) = out0ddat(1)/out0ddat(2)

!---enstrophy
CALL define_out0d_vals(out0ddat(4:4),1,ivarid=(/iarr_enstr0d/))
out0ddat(4) = 1.0E+06*out0ddat(4)

!---area bounded by surface 34.839 psu contour
scrit = 34.839
i1ptloc = COUNT(sal(1:ncloc,1:nrloc,nz).LT.scrit.AND.maskatc_int)
CALL sum_vars(i1ptloc,i1pt,0,commall=.TRUE.)
out0ddat(5) = 0.01*i1pt

!---minimum, maximum and mean salinity
CALL define_out0d_vals(out0ddat(6:8),3,ivarid=(/iarr_sal,iarr_sal,iarr_sal/),&
                     & oopt=(/oopt_min,oopt_max,oopt_mean/))
out0ddat(8) = out0ddat(8) - sal_ref

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
! Description - test case fredy
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
! Description - test case fredy
!
! Reference -
!
! Calling program - time_series
!
! Module calls - vorticity_3d
!
!************************************************************************
!
USE grid
USE modids
USE diagnostic_routines, ONLY: vorticity_3d
  
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
!*Local variables
!
INTEGER :: ii, jj


!---vorticity
out3ddat(5) = 0.0
ii_110: DO ii=i,i+1
jj_110: DO jj=j,j+1 
   out3ddat(5) = out3ddat(5) + 0.25*vorticity_3d(ii,jj,k)
ENDDO jj_110
ENDDO ii_110
out3ddat(5) = out3ddat(5)/coriolatu(1,1)


RETURN

END SUBROUTINE usrdef_tsr3d_vals
