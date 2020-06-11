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
! Description - test case bohai
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
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.8
!
! Description - test case bohai
!
! Reference -
!
! Calling program - time_series_init
!
! Module calls - inquire_var
!
!************************************************************************
!
USE gridpars
USE iopars
USE modids
USE switches
USE timepars
USE modvars_routines, ONLY: inquire_var
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!
!---fortran, standard and long name
CALL inquire_var(iarr_ekin0d,f90_name=tsrvars(1)%f90_name,&
               & standard_name=tsrvars(1)%standard_name,&
               & comment=tsrvars(1)%comment,long_name=tsrvars(1)%long_name)
CALL inquire_var(iarr_epot0d,f90_name=tsrvars(2)%f90_name,&
               & standard_name=tsrvars(2)%standard_name,&
               & comment=tsrvars(2)%comment,long_name=tsrvars(2)%long_name)
CALL inquire_var(iarr_etot0d,f90_name=tsrvars(3)%f90_name,&
               & standard_name=tsrvars(3)%standard_name,&
               & comment=tsrvars(3)%comment,long_name=tsrvars(3)%long_name)
CALL inquire_var(iarr_edissip0d,f90_name=tsrvars(4)%f90_name,&
               & standard_name=tsrvars(4)%standard_name,&
               & comment=tsrvars(4)%comment,long_name=tsrvars(4)%long_name)
CALL inquire_var(iarr_eflux2du,f90_name=tsrvars(7)%f90_name,&
               & standard_name=tsrvars(7)%standard_name,&
               & comment=tsrvars(7)%comment,long_name=tsrvars(7)%long_name)
CALL inquire_var(iarr_eflux2dv,f90_name=tsrvars(8)%f90_name,&
               & standard_name=tsrvars(8)%standard_name,&
               & comment=tsrvars(8)%comment,long_name=tsrvars(8)%long_name)
CALL inquire_var(iarr_etot2d,f90_name=tsrvars(10)%f90_name,&
               & standard_name=tsrvars(10)%standard_name,&
               & comment=tsrvars(10)%comment,long_name=tsrvars(10)%long_name)
CALL inquire_var(iarr_edissip2d,f90_name=tsrvars(11)%f90_name,&
               & standard_name=tsrvars(11)%standard_name,&
               & comment=tsrvars(11)%comment,long_name=tsrvars(11)%long_name)
CALL inquire_var(iarr_eflux3du,f90_name=tsrvars(15)%f90_name,&
               & standard_name=tsrvars(15)%standard_name,&
               & comment=tsrvars(15)%comment,long_name=tsrvars(15)%long_name)
CALL inquire_var(iarr_eflux3dv,f90_name=tsrvars(16)%f90_name,&
               & standard_name=tsrvars(16)%standard_name,&
               & comment=tsrvars(16)%comment,long_name=tsrvars(16)%long_name)
CALL inquire_var(iarr_eflux3dw,f90_name=tsrvars(17)%f90_name,&
               & standard_name=tsrvars(17)%standard_name,&
               & comment=tsrvars(17)%comment,long_name=tsrvars(17)%long_name)
CALL inquire_var(iarr_etot3d,f90_name=tsrvars(18)%f90_name,&
               & standard_name=tsrvars(18)%standard_name,&
               & comment=tsrvars(18)%comment,long_name=tsrvars(18)%long_name)
CALL inquire_var(iarr_edissip3d,f90_name=tsrvars(19)%f90_name,&
               & standard_name=tsrvars(19)%standard_name,&
               & comment=tsrvars(19)%comment,long_name=tsrvars(19)%long_name)

!---vector name
CALL inquire_var(iarr_eflux2du,vector_name=tsrvars(7)%vector_name)
CALL inquire_var(iarr_eflux2dv,vector_name=tsrvars(8)%vector_name)
CALL inquire_var(iarr_eflux3du,vector_name=tsrvars(15)%vector_name)
CALL inquire_var(iarr_eflux3dv,vector_name=tsrvars(16)%vector_name)
CALL inquire_var(iarr_eflux3dw,vector_name=tsrvars(17)%vector_name)

!---units
tsrvars(1)%units = 'PJ'
tsrvars(2)%units = 'PJ'
tsrvars(3)%units = 'PJ'
tsrvars(4)%units = 'GW'
tsrvars(7)%units = 'MW m-1'
tsrvars(8)%units = 'MW m-1'
tsrvars(10)%units = 'MJ m2'
tsrvars(11)%units = 'W m-2'
tsrvars(15)%units = 'MW m-2'
tsrvars(16)%units = 'MW m-2'
tsrvars(17)%units = 'MW m-2'
tsrvars(18)%units = 'MJ m-3'
tsrvars(19)%units = 'W m-3'

!---varids and ranks
tsrvars%ivarid = (/0,0,0,0,iarr_umvel,iarr_vmvel,0,0,iarr_zeta,0,0,iarr_uvel,&
                 & iarr_vvel,iarr_wphys,0,0,0,0,0/)
tsrvars%nrank = (/0,0,0,0,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3/)

!
!2. Variable indices
!-------------------
!

IF (iopt_verif.EQ.0) THEN
    ivarstsr(1,1:15) = (/5,6,7,8,9,10,11,12,13,14,15,16,17,18,19/)
    ivarstsr(2,1:4) = (/1,2,3,4/)
ELSE
    ivarstsr(1,:) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19/)
ENDIF

!
!3. File parameters
!------------------
!

IF (iopt_verif.EQ.0) THEN
    tsr2d(1)%defined = .TRUE.
    tsr0d(2)%defined = .TRUE.
ELSE
   tsr0d(1)%defined = .TRUE.
   tsr2d(1)%defined = .TRUE.
ENDIF
tsr3d(1)%defined = MERGE(.TRUE.,.FALSE.,iopt_grid_nodim.EQ.3)

!
!4. Output grid
!--------------
!

tsrgpars(1:nosetstsr)%time_format = 3
IF (iopt_verif.EQ.0) THEN
!  ---first file
   IF (iopt_grid_nodim.EQ.2) THEN
      tsrgpars(1)%zlims = 1
   ELSE
      tsrgpars(1)%zlims = (/1,nz,nz-1/)
   ENDIF
   IF (iopt_hydro_impl.EQ.0) THEN
      tsrgpars(1)%tlims = (/43200,nstep,720/)
   ELSE
      tsrgpars(1)%tlims = (/4320,nstep,72/)
   ENDIF
!  ---second file
   tsrgpars(2)%nodim = 0
   IF (iopt_hydro_impl.EQ.0) THEN
      tsrgpars(2)%tlims = (/0,nstep,120/)
   ELSE
      tsrgpars(2)%tlims = (/0,nstep,12/)
   ENDIF
ELSE
   tsrgpars(1)%nodim = MERGE(3,2,iopt_grid_nodim.EQ.3)
   IF (iopt_hydro_impl.EQ.0) THEN
      tsrgpars(1)%tlims = (/43200,nstep,720/)
   ELSE
      tsrgpars(1)%tlims = (/4320,nstep,72/)
   ENDIF
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
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.1.2
!
! Description - test case bohai
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


procname(pglev+1) = 'usrdef_tsr0d_vals'
CALL log_timer_in()

CALL define_out0d_vals(out0ddat,4,ivarid=(/iarr_ekin0d,iarr_epot0d,iarr_etot0d,&
                     & iarr_edissip0d/))

!---rescale to appropriate units
out0ddat(1:3) = 1.0E-09*out0ddat(1:3)
out0ddat(4) = 0.001*out0ddat(4)

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
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.1.2
!
! Description - test case bohai
!
! Reference -
!
! Calling program - time_series
!
! Module calls - define_out_2d_vals
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


CALL define_out2d_vals(out2ddat(3:4),i,j,2,ivarid=(/iarr_eflux2du,&
                     & iarr_eflux2dv/))
CALL define_out2d_vals(out2ddat(6:7),i,j,2,&
                     & ivarid=(/iarr_etot2d,iarr_edissip2d/))

!---rescale to appropriate units
out2ddat(3:4) = 1.0E-06*out2ddat(3:4)
out2ddat(6) = 1.0E-06*out2ddat(6)


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
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.1.2
!
! Description - test case bohai
!
! Reference -
!
! Calling program - time_series
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
!*n3vars*   INTEGER number of 2-D variables
!
!------------------------------------------------------------------------------
!


CALL define_out3d_vals(out3ddat(4:8),i,j,k,8,ivarid=(/iarr_eflux3du,&
                     & iarr_eflux3dv,iarr_eflux3dw,iarr_etot3d,iarr_edissip3d/))

!---rescale to appropriate units
out3ddat(4:7) = 1.0E-06*out3ddat(4:7)


RETURN

END SUBROUTINE usrdef_tsr3d_vals
