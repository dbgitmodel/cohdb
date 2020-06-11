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
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.7
!
! $Date: 2016-10-10 16:36:15 +0200 (Mon, 10 Oct 2016) $
!
! $Revision: 973 $
!
! Description - test case obctang
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
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.7
!
! Description - test case obctang
!
! Reference -
!
! Calling program - harmonic_analysis_init
!
!************************************************************************
!
USE iopars
USE tide
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_anal_freqs'
CALL log_timer_in()

!
!1. Harmonic frequencies
!-----------------------
!

index_anal(1) = icon_S2

!
!3. Specifiers for harmonic analysis
!-----------------------------------
!
!---number of frequencies per set
nofreqsharm = 1

!---frequency indices
ifreqsharm(:,1) = 1

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
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.7
!
! Description - test case river
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
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: l


procname(pglev+1) = 'usrdef_anal_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!
!---varids
analvars%ivarid = (/iarr_ekin0d,iarr_epot0d,iarr_etot0d,iarr_edissip0d,&
                  & iarr_zeta,iarr_umvel,iarr_vmvel,iarr_uvel,iarr_vvel,&
                  & iarr_wphys/)
!---ranks
analvars%nrank = (/0,0,0,0,2,2,2,3,3,3/)

!
!3. Variable indices
!-------------------
!

IF (iopt_grid_nodim.EQ.2) THEN
   ivarsanal(1,1:7) = (/(l,l=1,7)/)
   ivarsanal(2,1:3) = (/(l,l=5,7)/)
   ivarsell(1,1:4) = (/1,3,4,5/)
   ivarsell(2,1:4) = (/1,3,4,5/)
ELSE
   ivarsanal(1,1:10) = (/(l,l=1,10)/)
   ivarsanal(2,1:6) = (/(l,l=5,10)/)
   ivarsell(1,1:8) = (/1,3,4,5,8,10,11,12/)
   ivarsell(2,1:8) = (/1,3,4,5,8,10,11,12/)
ENDIF
ivecell2d(1,:) = (/2,3/)
ivecell2d(2,:) = (/2,3/)
IF (iopt_grid_nodim.EQ.3) THEN
   ivecell3d(1,:) = (/1,2/)
   ivecell3d(2,:) = (/1,2/)
ENDIF

!
!4. File parameters
!------------------
!
!---defined
res0d(1)%defined = .TRUE.
amp2d(1:2,1)%defined = .TRUE.; amp3d(1:2,1)%defined = iopt_grid_nodim.EQ.3
pha2d(1:2,1)%defined = .TRUE.; pha3d(1:2,1)%defined = iopt_grid_nodim.EQ.3
ell2d(1:2,1)%defined = .TRUE.; ell3d(1:2,1)%defined = iopt_grid_nodim.EQ.3

!---formats
res0d(1)%form = 'A'
amp2d(2,1)%form = 'A'; pha2d(2,1)%form = 'A'
ell2d(2,1)%form = 'A'

!
!5. Output grid
!--------------
!

analgpars(2)%gridded = .FALSE.
analgpars(2)%nostats = nostatsanal
analgpars(1)%tlims = (/nstep/2,nstep,nstep/2/)
analgpars(2)%tlims = (/nstep/2,nstep,nstep/2/)
analgpars(1:2)%time_format = 3

!
!6. Stations
!-----------
!
!---labels
lstatsanal(2,:) = (/(l,l=1,nostatsanal)/)

!---locations
SELECT CASE (runtitle(8:8))
   CASE ('A','C','E','G')
      analstatlocs%ipos = (/55,130,70,140,100,25,175,10,160/)
      analstatlocs%jpos = (/70,60,65,70,50,10,20,60,80/)
   CASE ('B','D','F','H')
      analstatlocs%ipos = (/9,159,39,179,99/)
      analstatlocs%jpos = (/89,69,79,89,49/)
END SELECT

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_anal_params

!========================================================================

SUBROUTINE usrdef_anal0d_vals(out0ddat,n0vars)
!************************************************************************
!
! *usrdef_anal0d_vals* 0-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - harmonic_analysis_update
!
! Module calls - define_out0d_vals
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
! Description -
!
! Reference -
!
! Calling program - harmonic_analysis_update
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
! Description -
!
! Reference -
!
! Calling program - harmonic_analysis_update
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
