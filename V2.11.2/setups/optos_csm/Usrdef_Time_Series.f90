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
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.8
!
! $Date: 2016-10-10 16:36:15 +0200 (Mon, 10 Oct 2016) $
!
! $Revision: 973 $
!
! Description - optos csm
!
! Reference -
!
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
! *usrdef_out_params* Specifiers for time series output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.8
!
! Description - optos csm
!
! Reference -
!
! Calling program - time_series_init
!
! Module calls - close_file, open_file 
!
!************************************************************************
!
USE gridpars
USE iopars
USE modids
USE timepars
USE switches
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iset, istat, iunit, l 

procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!
!---varids and ranks
tsrvars%ivarid = (/iarr_zeta,iarr_umvel,iarr_vmvel/)
tsrvars%nrank = 2

!
!2. Variable indices
!-------------------
!

ivarstsr(1,1:3) = (/1,2,3/)
IF (iopt_verif.EQ.0) ivarstsr(2,1:3) = (/1,2,3/)

!
!3. File parameters
!------------------
!

IF (iopt_verif.EQ.0) THEN
   tsr2d(1:nosetstsr)%defined = .TRUE.
ELSE
   tsr2d(1)%defined = .TRUE.
ENDIF

!
!4. Output grid
!--------------
!
!---first file
tsrgpars(1)%tlims(1:2) = (/0,nstep/)
tsrgpars(1)%tlims(3) = MERGE(540,36,iopt_hydro_impl.EQ.0)

!---second file
IF (iopt_verif.EQ.0) THEN
   tsrgpars(2)%gridded = .FALSE.
   tsrgpars(2)%nostats = nostatstsr
   tsrgpars(2)%tlims(1:2) = (/0,nstep/)
   tsrgpars(2)%tlims(3) = MERGE(30,2,iopt_hydro_impl.EQ.0)
ENDIF

tsrgpars(1:nosetstsr)%time_format = 4

!
!5. Station labels
!------------------
!

IF (iopt_verif.EQ.0) THEN
   lstatstsr(2,1:nostatstsr) = (/(l,l=1,nostatstsr)/)
ENDIF

!
!6. Station locations
!--------------------
!

IF (iopt_verif.EQ.0) THEN
   CALL open_file(iunit,'csm_stations','IN','A')
   istat_610: DO istat=1,nostatstsr
      READ (iunit,9001) tsrstatlocs(istat)%name, tsrstatlocs(istat)%ipos, &
                      & tsrstatlocs(istat)%jpos
   ENDDO istat_610
   CALL close_file(iunit,'A')
ENDIF

CALL log_timer_out()


RETURN

9001 FORMAT (A19,2I4)

END SUBROUTINE usrdef_tsr_params

!========================================================================

SUBROUTINE usrdef_tsr0d_vals(out0ddat,n0vars)
!************************************************************************
!
! *usrdef_tsr0d_vals* 0-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.0
!
! Description - optos csm
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
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.1.2
!
! Description - optos csm
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
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.0
!
! Description - optos csm
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
