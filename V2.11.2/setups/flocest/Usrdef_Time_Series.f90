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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.11.1
!
! $Date: 2014-12-05 17:20:26 +0100 (Fri, 05 Dec 2014) $
!
! $Revision: 780 $
!
! Description - test case flocest 
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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.11.1
!
! Description - test case flocest
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
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: l

procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!
!---fortran name
tsrvars(1)%f90_name = 'Pmass0d'
tsrvars(2)%f90_name = 'Fmass0d'
tsrvars(3)%f90_name = 'Tmass0d'

!---standard name
tsrvars(1)%standard_name = 'total_mass_flocculi_in_the_domain'
tsrvars(2)%standard_name = 'total_mass_macroflocs_in_the_domain'
tsrvars(3)%standard_name = 'total_sediment_mass_in_the_domain'

!---long name
tsrvars(1)%long_name = 'Total mass of flocculi in the domain'
tsrvars(2)%long_name = 'Total mass of macroflocs in the domain'
tsrvars(3)%long_name = 'Total sediment mass in the domain'

!---units
tsrvars(1)%units = 'kg'
tsrvars(2)%units = 'kg'
tsrvars(3)%units = 'kg'

!---varids and ranks
tsrvars%ivarid = (/0,0,0,iarr_zeta,iarr_umvel,iarr_vmvel,iarr_bstresatc,&
                 & iarr_uvel,iarr_vvel,iarr_wphys,iarr_sal,iarr_floc_P,&
                 & iarr_floc_F,iarr_floc_dia,iarr_floc_dens,iarr_floc_nc,&
                 & iarr_wfall/)
tsrvars(1:3)%nrank = 0
tsrvars(4:7)%nrank = 2
tsrvars(8:17)%nrank = 3
tsrvars(17)%numvar = 2

!
!2. Variable indices
!-------------------
!

ivarstsr(1,1:14) = (/(l,l=4,17)/)
ivarstsr(2,1:3) = (/1,2,3/)

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
!---first file
tsrgpars(1)%time_format = 3
tsrgpars(1)%tlims = (/0,nstep,1200/)

!---second file
tsrgpars(2)%time_format = 3
tsrgpars(2)%tlims = (/0,nstep,100/)

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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.11.1
!
! Description -
!
! Reference -
!
! Calling program - time_series
!
! Module calls - define_out0d_vals
!
!************************************************************************
!
USE iopars
USE sedids  
USE model_output, ONLY: define_out0d_vals  
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

CALL define_out0d_vals(out0ddat(1:2),2,ivarid=(/iarr_floc_P,iarr_floc_F/),&
                     & oopt=(/oopt_int,oopt_int/))
out0ddat(3) = out0ddat(1) + out0ddat(2)
out0ddat = 0.001*out0ddat

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

