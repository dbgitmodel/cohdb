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

SUBROUTINE usrdef_output
!************************************************************************
!
! *usrdef_output* User-formatted output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.11.2
!
! $Date: 2018-07-23 16:55:25 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1171 $
!
! Description - output parameters for test case seamount
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - 
!
! Module calls - close_file, energy_0d, max_vars, mult_index, open_file,
!                sum2_vars, vector_mag_arr_atc
!
!************************************************************************
!
USE currents  
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE syspars
USE timepars
USE diagnostic_routines, ONLY: energy_0d
USE math_library, ONLY: vector_mag_arr_atc
USE inout_routines, ONLY: close_file, open_file
USE paral_utilities, ONLY: max_vars, sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER, SAVE :: icout, iunit
INTEGER :: k
INTEGER, SAVE, DIMENSION(4) :: nhdims
REAL :: curmax, curmaxdev, curmax2d, curmean, curmeandev, curmean2d, edens, &
      & ekindev, ekintot, ekinrat, ekin2d, rday
REAL, DIMENSION(5) :: ecomps
REAL, DIMENSION(ncloc+1,nrloc,nz) :: array3du
REAL, DIMENSION(ncloc,nrloc+1,nz) :: array3dv
REAL, DIMENSION(ncloc,nrloc) :: curmag2d
REAL, DIMENSION(ncloc,nrloc,nz) :: curmag


!
!1. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN

!  ---open output file
   IF (master) THEN
      CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case seamount:'//&
           & ' simulation '//TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF

!  ---initialise parameters
   icout = NINT(21600/delt3d)
   nhdims = 0

ENDIF

IF (nt.GT.0.AND.(.NOT.mult_index(nt,icout))) RETURN

CALL log_timer_in()
procname(pglev+1) = 'usrdef_output'
 
!
!2. Evaluate parameters
!----------------------
!
!---day number
rday = nosecsrun/86400.0

!---energy components
CALL energy_0d(3,ecomps)
ekintot = 1.0E-09*ecomps(1)
edens = 1.0E-09*ecomps(4)
CALL energy_0d(2,ecomps)
ekin2d = 1.0E-09*ecomps(1)
ekindev = ekintot - ekin2d
IF (ekin2d.GT.0.0) THEN
   ekinrat = ekindev/ekin2d
ELSE
   ekinrat = 0.0
ENDIF

!---3-D current mean and maximum
CALL vector_mag_arr_atc(uvel(1:ncloc+1,1:nrloc,:),vvel(1:ncloc,1:nrloc+1,:),&
                      & 2,1,nz,1,iarr_hvelmag,.TRUE.,vecmag=curmag)
CALL max_vars(curmag,curmax,iarr_hvelmag,mask=maskatc_int)
CALL sum2_vars(curmag,curmean,nhdims,'C  ',0)
curmean = curmean/(noseaatc*nz)

!---2-D current mean and maximum
CALL vector_mag_arr_atc(umvel(1:ncloc+1,1:nrloc),vmvel(1:ncloc,1:nrloc+1),&
                      & 2,1,1,1,iarr_hmvelmag,.TRUE.,vecmag=curmag2d)
CALL max_vars(curmag2d,curmax2d,iarr_hvelmag,mask=maskatc_int)
CALL sum2_vars(curmag2d,curmean2d,nhdims,'C  ',0)
curmean2d = curmean2d/noseaatc

!---3-D baroclinic current mean and maximum
k_110: DO k=1,nz
   array3du(:,:,k) = uvel(1:ncloc+1,1:nrloc,k)- umvel(1:ncloc+1,1:nrloc)
   array3dv(:,:,k) = vvel(1:ncloc,1:nrloc+1,k)- vmvel(1:ncloc,1:nrloc+1)
ENDDO k_110
CALL vector_mag_arr_atc(array3du,array3dv,2,1,nz,1,0,.TRUE.,vecmag=curmag)
CALL max_vars(curmag,curmaxdev,0,mask=maskatc_int)
CALL sum2_vars(curmag,curmeandev,nhdims,'C  ',0)
curmeandev = curmeandev/(noseaatc*nz)

!
!3. Write output
!---------------
!

IF (master) THEN
   WRITE (iunit,9001) rday
   WRITE (iunit,9002) 'ekintot', ekintot
   WRITE (iunit,9002) 'ekin2d', ekin2d
   WRITE (iunit,9002) 'ekindev', ekindev
   WRITE (iunit,9002) 'ekinrat', ekinrat
   WRITE (iunit,9002) 'edens', edens
   WRITE (iunit,9002) 'curmean', 100.0*curmean
   WRITE (iunit,9002) 'curmean2d', 100.0*curmean2d
   WRITE (iunit,9002) 'curmeandev', 100.0*curmeandev
   WRITE (iunit,9002) 'curmax', 100.0*curmax
   WRITE (iunit,9002) 'curmax2d', 100.0*curmax2d
   WRITE (iunit,9002) 'curmaxdev', 100.0*curmaxdev
ENDIF

!
!4. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')

CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',F4.2,' days')
9002 FORMAT (T3,A,T13,': ',G12.4E3)

END SUBROUTINE usrdef_output
