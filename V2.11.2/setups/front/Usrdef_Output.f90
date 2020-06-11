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
! Version - @(COHERENS)Usrdef_Output.f90  V2.1.1
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - output parameters for test case front
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_mod, mult_index, open_file
!
!************************************************************************
!
USE density
USE gridpars
USE iopars
USE modids
USE paralpars
USE syspars
USE timepars
USE inout_routines, ONLY: close_file, open_file
USE paral_comms, ONLY: combine_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER :: i, iflag, ihour, il, ir, iunit = 0, k, kc1, kc2, kleft, kright
REAL :: cmax, cmin, dc, dcdz, depthl, depthr, dx, gradh, gradmax, gradv, &
      & hleng, zsur
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: depglb
REAL, DIMENSION(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo,nz) :: salglb


!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,1)%status = 'R'
   modfiles(io_modgrd,1,1)%form = modfiles(io_modgrd,1,2)%form
   modfiles(io_modgrd,1,1)%filename = TRIM(intitle)//'.modgrd.bin'
   modfiles(io_modgrd,1,2)%status = '0'
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0.AND.master) THEN
!  ---open output file
   IF (.NOT.cold_start) THEN
      CALL open_file(iunit,runtitle(1:6)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case front: simulation '&
                        & //TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF
!  ---global water depths
   ALLOCATE (depglb(nc-1))
   depthl = 50.0
   depthr = 30.0
   i_210: DO i=1,nc-1
      IF (i.GT.30) THEN
         depglb(i) = depthr
      ELSEIF (i.GT.25) THEN
         depglb(i) = depthl - (i-25)*(depthl-depthr)/5.0
      ELSE
         depglb(i) = depthl
      ENDIF
   ENDDO i_210
ENDIF

!
!3. Store parameters
!-------------------
!

IF (mult_index(nt,60)) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

   ihour = nosecsrun/3600

!  ---combine salinity at master
   CALL combine_mod(salglb,sal,(/1-nhalo,1-nhalo,1/),iarr_sal,0.0)
   IF (.NOT.master) GOTO 1000

!  ---horizontal parameters
   dx = 1.0
   gradmax = 0.0
   iflag = 0
   kright = NINT(nz+0.5-5.0*nz/depglb(nc-1))
   kleft = NINT(nz+0.5-5.0*nz/depglb(1))
   i_310: DO i=nc-1,2,-1
      kc1 = NINT(nz+0.5-5.0*nz/depglb(i))
      kc2 = NINT(nz+0.5-5.0*nz/depglb(i-1))
      dc = ABS(salglb(i,1,kc1)-salglb(i-1,1,kc2))
      gradmax = MAX(gradmax,dc)
      IF (iflag.EQ.0) THEN
         IF (ABS(salglb(i,1,kc1)-salglb(nc-1,1,kright)).GT.0.01) THEN
            iflag = 1
            ir = i
         ENDIF
      ELSEIF (iflag.EQ.1) THEN
         IF (ABS(salglb(i,1,kc1)-salglb(1,1,kleft)).LT.0.01) THEN
            iflag = 2
            il = i
         ENDIF
      ENDIF
   ENDDO i_310
   IF (iflag.EQ.1) il = 1
   gradh = gradmax/dx
   hleng = (ir-il)*dx

!  ---vertical parameters
   gradv = 0.0
   k_320: DO k=1,nz-1
      zsur = depglb(il)*(1.0-0.5*(2*k-1)/REAL(nz))
      IF (zsur.LT.30.0) THEN
         dcdz = nz*ABS(salglb(il,1,k+1)-salglb(il,1,k))/depglb(il)
         gradv = MAX(gradv,dcdz)
      ENDIF
   ENDDO k_320

!  ---max and min conc along vertical transect
   cmin = salglb(il,1,nz)
   cmax = salglb(il,1,nz)
   k_330: DO k=nz-1,1,-1
      zsur = depglb(il)*(1.0-0.5*(2*k-1)/REAL(nz))
      IF (zsur.LT.30.0) THEN
         IF (salglb(il,1,k).LT.cmin) cmin = salglb(il,1,k)
         IF (salglb(il,1,k).GT.cmax) cmax = salglb(il,1,k)
      ENDIF
   ENDDO k_330

!  ---write output
   IF (.NOT.cold_start) THEN
      WRITE (iunit,9001) ihour
      WRITE (iunit,9002) 'gradh', gradh
      WRITE (iunit,9002) 'gradv', gradv
      WRITE (iunit,9002) 'hleng', hleng
      WRITE (iunit,9002) 'cmin', cmin
      WRITE (iunit,9002) 'cmax', cmax
   ENDIF

1000 CALL log_timer_out

ENDIF

!
!4. Close output file
!--------------------
!

IF (master) THEN
   IF (nt.EQ.nstep) CALL close_file(iunit,'A')
   IF (nt.EQ.nstep.OR.cold_start) DEALLOCATE (depglb)
ENDIF

RETURN

9001 FORMAT ('Time: ',I2,' hours')
9002 FORMAT (T3,A,T12,': ',G12.4E3)

END SUBROUTINE usrdef_output
