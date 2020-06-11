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
! Description - test case bohai
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
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.0
!
! Description - test case bohai
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


procname(pglev+1) = 'usrdef_anal_freqs'
CALL log_timer_in()

!
!1. Harmonic frequencies
!-----------------------
!

index_anal(1) = icon_M2
index_anal(2) = icon_S2

!
!3. Specifiers for harmonic analysis
!-----------------------------------
!
!---number of frequencies per set
nofreqsharm = 2

!---frequency indices
ifreqsharm(1,:) = (/1,2/)
ifreqsharm(2,:) = (/1,2/)

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
! Description - test case bohai
!
! Reference -
!
! Calling program - harmonic_analysis_init
!
! Module calls - inquire_var
!
!************************************************************************
!
USE gridpars
USE iopars
USE modids
USE switches
USE modvars_routines, ONLY: inquire_var
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
CALL inquire_var(iarr_eflux2du,f90_name=analvars(7)%f90_name,&
               & standard_name=analvars(7)%standard_name,&
               & comment=analvars(7)%comment,long_name=analvars(7)%long_name)
CALL inquire_var(iarr_eflux2dv,f90_name=analvars(8)%f90_name,&
               & standard_name=analvars(8)%standard_name,&
               & comment=analvars(8)%comment,long_name=analvars(8)%long_name)
CALL inquire_var(iarr_etot2d,f90_name=analvars(10)%f90_name,&
               & standard_name=analvars(10)%standard_name,&
               & comment=analvars(10)%comment,long_name=analvars(10)%long_name)
CALL inquire_var(iarr_edissip2d,f90_name=analvars(11)%f90_name,&
               & standard_name=analvars(11)%standard_name,&
               & comment=analvars(11)%comment,long_name=analvars(11)%long_name)
IF (iopt_grid_nodim.EQ.3) THEN
   CALL inquire_var(iarr_eflux3du,f90_name=analvars(15)%f90_name,&
               & standard_name=analvars(15)%standard_name,&
               & comment=analvars(15)%comment,long_name=analvars(15)%long_name)
   CALL inquire_var(iarr_eflux3dv,f90_name=analvars(16)%f90_name,&
               & standard_name=analvars(16)%standard_name,&
               & comment=analvars(16)%comment,long_name=analvars(16)%long_name)
   CALL inquire_var(iarr_eflux3dw,f90_name=analvars(17)%f90_name,&
               & standard_name=analvars(17)%standard_name,&
               & comment=analvars(17)%comment,long_name=analvars(17)%long_name)
   CALL inquire_var(iarr_etot3d,f90_name=analvars(18)%f90_name,&
               & standard_name=analvars(18)%standard_name,&
               & comment=analvars(18)%comment,long_name=analvars(18)%long_name)
   CALL inquire_var(iarr_edissip3d,f90_name=analvars(19)%f90_name,&
               & standard_name=analvars(19)%standard_name,&
               & comment=analvars(19)%comment,long_name=analvars(19)%long_name)
ENDIF 

!---vector name
CALL inquire_var(iarr_eflux2du,vector_name=analvars(7)%vector_name)
CALL inquire_var(iarr_eflux2dv,vector_name=analvars(8)%vector_name)
IF (iopt_grid_nodim.EQ.3) THEN
   CALL inquire_var(iarr_eflux3du,vector_name=analvars(15)%vector_name)
   CALL inquire_var(iarr_eflux3dv,vector_name=analvars(16)%vector_name)
   CALL inquire_var(iarr_eflux3dw,vector_name=analvars(17)%vector_name)
ENDIF
  
!---units
analvars(1)%units = 'PJ'
analvars(2)%units = 'PJ'
analvars(3)%units = 'PJ'
analvars(4)%units = 'GW'
analvars(7)%units = 'MW m-1'
analvars(8)%units = 'MW m-1'
analvars(10)%units = 'MJ m-2'
analvars(11)%units = 'W m-2'
IF (iopt_grid_nodim.EQ.3) THEN
   analvars(15)%units = 'MW m-2'
   analvars(16)%units = 'MW m-2'
   analvars(17)%units = 'MW m-2'
   analvars(18)%units = 'MJ m-3'
   analvars(19)%units = 'W m-3'
ENDIF

!---varids and ranks
IF (iopt_grid_nodim.EQ.2) THEN
   analvars%ivarid = (/0,0,0,0,iarr_umvel,iarr_vmvel,0,0,iarr_zeta,0,0/)
   analvars%nrank = (/0,0,0,0,2,2,2,2,2,2,2/)
ELSE
   analvars%ivarid = (/0,0,0,0,iarr_umvel,iarr_vmvel,0,0,iarr_zeta,0,0,&
                     & iarr_uvel,iarr_vvel,iarr_wphys,0,0,0,0,0/)
   analvars%nrank = (/0,0,0,0,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3/)
ENDIF

!
!3. Variable indices
!-------------------
!

IF (iopt_grid_nodim.EQ.2) THEN
   ivarsanal(1,1:11) = (/1,2,3,4,5,6,7,8,9,10,11/)
   ivarsell(1,1:4) = (/1,3,4,5/)
ELSE
   ivarsanal(1,1:19) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19/)
   ivarsell(1,1:8) = (/1,3,4,5,8,10,11,12/)
ENDIF
ivarsanal(2,1:3) = (/5,6,9/)
ivarsell(2,1:4) = (/1,3,4,5/)

ivecell2d(1,:) = (/1,2/)
ivecell2d(2,:) = (/1,2/)
ivecell3d(1,:) = (/1,2/)

!
!4. File parameters
!------------------
!
!4.1 Residuals
!-------------
!
!---residuals
res0d(1)%defined = .TRUE.; res2d(1)%defined = .TRUE.
res3d(1)%defined = iopt_grid_nodim.EQ.3
res0d(1)%form = 'A'

!---amplitudes, phases
amp2d(1:2,:)%defined = .TRUE.; pha2d(1:2,:)%defined = .TRUE.

!---elliptic parameters
ell2d(1:2,:)%defined = .TRUE.; ell3d(1,:)%defined = iopt_grid_nodim.EQ.3

!---file formats
res0d(1)%form = 'A'
amp2d(2,:)%form = 'A'; pha2d(2,:)%form = 'A'
ell2d(2,:)%form = 'A'

!
!5. Output grid
!--------------
!
!---full grid
analgpars(1)%nodim = MERGE(3,2,iopt_grid_nodim.EQ.3)
analgpars(1)%time_format = 4
IF (iopt_grid_nodim.EQ.2) THEN
   analgpars(1)%zlims = 1
ELSE
   analgpars(1)%zlims = (/1,nz,nz-1/)
ENDIF
IF (iopt_hydro_impl.EQ.0) THEN
   analgpars(1)%tlims =  (/5760,46080,40320/)
ELSE
   analgpars(1)%tlims =  (/576,4608,4032/)
ENDIF

!---stations
analgpars(2)%gridded = .FALSE.
analgpars(2)%nostats = nostatsanal
analgpars(2)%time_format = 3
IF (iopt_hydro_impl.EQ.0) THEN
   analgpars(2)%tlims =  (/5760,46080,40320/)
ELSE
   analgpars(2)%tlims =  (/576,4608,4032/)
ENDIF

!
!6. Stations
!-----------
!
!---labels
lstatsanal(2,:) = (/(l,l=1,nostatsanal)/)

!---locations
analstatlocs%ipos = (/10,88,91,49/)
analstatlocs%jpos = (/22,10,28,40/)


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
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.1.2
!
! Description - test case bohai
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
! Description - test case bohai
!
! Reference -
!
! Calling program - harmonic_analysis_update
!
! Module calls - define_out2d_vals
!
!************************************************************************
!
USE modids
USE model_output, ONLY: define_out2d_vals

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


CALL define_out2d_vals(out2ddat(3:4),i,j,2,&
                     & ivarid=(/iarr_eflux2du,iarr_eflux2dv/))
CALL define_out2d_vals(out2ddat(6:7),i,j,2,&
                     & ivarid=(/iarr_etot2d,iarr_edissip2d/))

!---rescale to appropriate units
out2ddat(3:4) = 1.0E-06*out2ddat(3:4)
out2ddat(6) = 1.0E-06*out2ddat(6)


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
! Description - test case bohai
!
! Reference -
!
! Calling program - harmonic_analysis_update
!
! Module calls - define_out3d_vals
!
!************************************************************************
!
USE modids
USE model_output, ONLY: define_out3d_vals
  
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


CALL define_out3d_vals(out3ddat(4:8),i,j,k,5,&
                     & ivarid=(/iarr_eflux3du,iarr_eflux3dv,iarr_eflux3dw,&
                     & iarr_etot3d,iarr_edissip3d/))

!---rescale to appropriate units
out3ddat(4:7) = 1.0E-06*out3ddat(4:7)


RETURN

END SUBROUTINE usrdef_anal3d_vals
