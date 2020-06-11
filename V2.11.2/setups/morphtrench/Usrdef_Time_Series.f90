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
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.10
!
! $Date: 2016-10-10 16:36:15 +0200 (Mon, 10 Oct 2016) $
!
! $Revision: 973 $
!!
! Description - test case morphtrench
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
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.10
!
! Description - test case morphtrench
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
USE sedswitches
USE timepars
USE modvars_routines, ONLY: inquire_var
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
!---metadata
tsrvars(8)%f90_name = 'seafloor'
tsrvars(8)%long_name = 'bed topography'
tsrvars(8)%standard_name = 'bed_topography'
tsrvars(8)%units = 'm'
tsrvars(8)%comment = 'standard name constructed following the netCDF CF-1.6'//&
                   & ' guidelines'
IF (iopt_sed_mode.EQ.3) THEN
   CALL inquire_var(iarr_ctot,f90_name=tsrvars(10)%f90_name,&
               & standard_name=tsrvars(10)%standard_name,&
               & comment=tsrvars(10)%comment,long_name=tsrvars(10)%long_name,&
               & units=tsrvars(10)%units)
   tsrvars(10)%f90_name = TRIM(tsrvars(10)%f90_name)//'_1'
   tsrvars(10)%standard_name = TRIM(tsrvars(10)%standard_name)//'_1'
ENDIF

!---varids
IF (iopt_sed_mode.EQ.1) THEN
   tsrvars%ivarid = (/iarr_zeta,iarr_umvel,iarr_vmvel,iarr_bed_update_int,&
                    & iarr_bstresatc,iarr_qbedatu,iarr_qbedatv,0/)
ELSE
   tsrvars%ivarid = (/iarr_zeta,iarr_umvel,iarr_vmvel,iarr_bed_update_int,&
                    & iarr_bstresatc,iarr_qbedatu,iarr_qbedatv,&
                    & 0,iarr_bottom_sed_flux,0/)
ENDIF

!---rank
tsrvars%nrank = 2

!---variable number
tsrvars(6:7)%numvar = 1
IF (iopt_sed_mode.EQ.3) tsrvars(9)%numvar = 1

!
!2. Variable indices
!-------------------
!

ivarstsr(1,1:novarstsr) = (/(l,l=1,novarstsr)/)

!
!3. File parameters
!------------------
!
!---first file
tsr2d(1)%defined = .TRUE.

!
!4. Output grid
!--------------
!
!---first file
tsrgpars(1)%time_format = 3
tsrgpars(1)%tlims = (/0,nstep,21600/)

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_tsr_params

!========================================================================

SUBROUTINE usrdef_tsr0d_vals(out0ddat,n0vars)
!************************************************************************
!
! *usrdef_tsr0d_vals* 0-D output data
!
! Author -
!
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.1.2
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


RETURN

END SUBROUTINE usrdef_tsr0d_vals

!========================================================================

SUBROUTINE usrdef_tsr2d_vals(out2ddat,i,j,n2vars)
!************************************************************************
!
! *usrdef_tsr2d_vals* 2-D output data
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.10
!
! Description - test case morphtrench
!
! Reference -
!
! Calling program - time_series
!
!************************************************************************
!
USE depths
USE iopars
USE sedarrays
USE sedswitches
USE timepars

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
!*Local variables
!
REAL, SAVE :: depmax


IF (nt.EQ.0) depmax = 6.0

out2ddat(8) = depmax - depmeanatc(i,j)
IF (iopt_sed_mode.EQ.3) out2ddat(9) = ctot(i,j,1)


RETURN

END SUBROUTINE usrdef_tsr2d_vals

!========================================================================

SUBROUTINE usrdef_tsr3d_vals(out3ddat,i,j,k,n3vars)
!************************************************************************
!
! *usrdef_tsr3d_vals* 3-D output data
!
! Author -
!
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.1.2
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
!*Local variables
!


RETURN

END SUBROUTINE usrdef_tsr3d_vals
