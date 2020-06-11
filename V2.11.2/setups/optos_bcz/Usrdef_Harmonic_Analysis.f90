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
! *Usrdef_Harmonic_Analysis* Define harmonic analysis and output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.8
!
! $Date: 2016-10-10 16:36:15 +0200 (Mon, 10 Oct 2016) $
!
! $Revision: 973 $
!
! Description - test case optos bcz
!
! Reference -
!
! Routines - usrdef_anal_freqs, usrdef_anal_params, usrdef_anal0d_vals,
!            usrdef_anal2d_vals, usrdef_anal3d_vals
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_anal_freqs
!************************************************************************
!
! *usrdef_anal_freqs* Harmonic frequencies
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.3
!
! Description - test case optos bcz
!
! Reference -
!
! Calling program - harmonic_analysis_init
!
!************************************************************************
!
USE iopars
USE modids
USE switches
USE tide
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: icon

procname(pglev+1) = 'usrdef_anal_freqs'
CALL log_timer_in()

!
!1. Harmonic frequencies
!-----------------------
!

index_anal(1:nofreqsanal) = (/icon_M2, icon_S2, icon_N2, icon_K2, &
                            & icon_O1, icon_K1, icon_Q1, icon_P1/)

!
!3. Specifiers for harmonic analysis
!-----------------------------------
!
!---number of frequencies per set
nofreqsharm = 8

!---frequency indices
ifreqsharm(1,:) = (/(icon,icon=1,nofreqsanal)/)
ifreqsharm(2,:) = ifreqsharm(1,:)

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_anal_freqs

!========================================================================

SUBROUTINE usrdef_anal_params
!************************************************************************
!
! *usrdef_anal_params* Specifiers for harmonic output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.8
!
! Description - test case optos bcz
!
! Reference -
!
! Calling program - harmonic_analysis_init
!
! Module calls - close_file, inquire_var, open_file
!
!************************************************************************
!
USE gridpars
USE iopars
USE modids
USE switches
USE timepars
USE modvars_routines, ONLY: inquire_var
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: icount, istat, iunit, l


procname(pglev+1) = 'usrdef_anal_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!
!---fortran, standard and long name
CALL inquire_var(iarr_ekin0d,f90_name=analvars(1)%f90_name,&
               & standard_name=analvars(1)%standard_name,&
               & comment=analvars(1)%comment,long_name=analvars(1)%long_name)
CALL inquire_var(iarr_epot0d,f90_name=analvars(2)%f90_name,&
               & standard_name=analvars(2)%standard_name,&
               & comment=analvars(2)%comment,long_name=analvars(2)%long_name)
CALL inquire_var(iarr_etot0d,f90_name=analvars(3)%f90_name,&
               & standard_name=analvars(3)%standard_name,&
               & comment=analvars(3)%comment,long_name=analvars(3)%long_name)
CALL inquire_var(iarr_edissip0d,f90_name=analvars(4)%f90_name,&
               & standard_name=analvars(4)%standard_name,&
               & comment=analvars(4)%comment,long_name=analvars(4)%long_name)
   
!---units
analvars(1)%units = 'PJ'
analvars(2)%units = 'PJ'
analvars(3)%units = 'PJ'
analvars(4)%units = 'GW'

!---varids and ranks
analvars%ivarid = (/0,0,0,0,iarr_umvel,iarr_vmvel,iarr_zeta/)
analvars%nrank = (/0,0,0,0,2,2,2/)

!
!3. Variable indices
!-------------------
!
!---first file
ivarsanal(1,1:7) = (/1,2,3,4,5,6,7/)
ivarsell(1,1:4) = (/1,3,4,5/)

!---second file
ivarsanal(2,1:3) = (/5,6,7/)
ivarsell(2,1:4) = (/1,3,4,5/)

!---elliptic vectors
ivecell2d(1,:) = (/1,2/)
ivecell2d(2,:) = (/1,2/)

!
!4. File parameters
!------------------
!
!4.1 Residuals
!-------------
!
!---residuals
res0d(1)%defined = .TRUE.; res2d(1)%defined = .TRUE.

!---amplitudes, phases
amp2d(1:2,:)%defined = .TRUE.; pha2d(1:2,:)%defined = .TRUE.

!---elliptic parameters
ell2d(1:2,:)%defined = .TRUE.

!---file formats
res0d(1)%form = 'A'
amp2d(2,:)%form = 'A'; pha2d(2,:)%form = 'A'; ell2d(2,:)%form = 'A'

!
!5. Output grid
!--------------
!
!---full grid
analgpars(1)%time_format = 4
icount = MERGE(8640,576,iopt_hydro_impl.EQ.0)
analgpars(1)%tlims = (/icount,nstep,nstep-icount/)

!---stations
analgpars(2)%gridded = .FALSE.
analgpars(2)%nostats = nostatsanal
analgpars(2)%time_format = 4
analgpars(2)%tlims = (/icount,nstep,nstep-icount/)

!
!6. Stations
!-----------
!
!---labels
lstatsanal(2,:) = (/(l,l=1,nostatsanal)/)

!---locations
CALL open_file(iunit,'bcz_stations','IN','A')
istat_610: DO istat=1,nostatsanal
   READ (iunit,9001) analstatlocs(istat)%name, analstatlocs(istat)%ipos, &
                   & analstatlocs(istat)%jpos
ENDDO istat_610
CALL close_file(iunit,'A')

CALL log_timer_out()


RETURN

9001 FORMAT (A19,2I4)

END SUBROUTINE usrdef_anal_params

!========================================================================

SUBROUTINE usrdef_anal0d_vals(out0ddat,n0vars)
!************************************************************************
!
! *usrdef_anal0d_vals* 0-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.3
!
! Description - test case optos bcz
!
! Reference -
!
! Calling program - harmonic_analysis_update
!
! Module calls - define_out0d_vals
!
!************************************************************************
!
USE iopars
USE modids
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


procname(pglev+1) = 'usrdef_anal0d_vals'
CALL log_timer_in()

CALL define_out0d_vals(out0ddat,4,ivarid=(/iarr_ekin0d,iarr_epot0d,iarr_etot0d,&
                     & iarr_edissip0d/))

!---rescale to appropriate units
out0ddat(1:3) = 1.0E-09*out0ddat(1:3)
out0ddat(4) = 0.001*out0ddat(4)

CALL log_timer_out()   

      
RETURN

END SUBROUTINE usrdef_anal0d_vals

!========================================================================

SUBROUTINE usrdef_anal2d_vals(out2ddat,i,j,n2vars)
!************************************************************************
!
! *usrdef_anal2d_vals* 2-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.1.2
!
! Description - test case optos bcz
!
! Reference -
!
! Calling program - harmonic_analysis_update
!
! Module calls -
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

END SUBROUTINE usrdef_anal2d_vals

!========================================================================

SUBROUTINE usrdef_anal3d_vals(out3ddat,i,j,k,n3vars)
!************************************************************************
!
! *usrdef_anal3d_vals* 3-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.1.2
!
! Description - test case optos bcz
!
! Reference -
!
! Calling program - harmonic_analysis_update
!
! Module calls -
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
!*n3vars*   INTEGER number of 3-D variables
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_anal3d_vals
