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
! Version - @(COHERENS)Usrdef_Output.f90  V2.4
!
! $Date: 2017-09-08 09:49:20 +0200 (Fri, 08 Sep 2017) $
!
! $Revision: 1048 $
!
! Description - output parameters for test case river
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_mod, copy_vars, mult_index, num_proc,
!                open_file
!
!************************************************************************
!
USE currents
USE density
USE depths
USE gridpars
USE iopars
USE modids
USE obconds
USE paralpars
USE physpars
USE tide
USE timepars
USE grid_routines, ONLY: num_proc
USE inout_routines, ONLY: close_file, open_file
USE paral_comms, ONLY: combine_mod, copy_vars
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index


IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER :: i, ic, idproc, ih, ihloc, ihour, ip, k, kc, kp
INTEGER, SAVE :: icout, iunit, knt, mtot
REAL :: buamp, bupha, bures, denom, dep, dhun, fcos, freq, fsin, grad, &
      & gradmax, hfront, hgrad, pdep, phas, scrit, suamp, sumc, sumcc, sumss, &
      & supha, sures, vgrad, xc
REAL, SAVE :: ubsum, ubsumc, ubsums, ussum, ussumc, ussums
INTEGER, DIMENSION(4) :: nhdims
REAL, DIMENSION(0:nc+1,0:nr+1) :: deptotglb
REAL, DIMENSION(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo,nz) :: salglb


IF (runtitle(6:6).EQ.'0') RETURN

procname(pglev+1) = 'usrdef_output'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN
!  ---open output file
   IF (.NOT.cold_start.AND.master) THEN
      CALL open_file(iunit,TRIM(outtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case river: simulation '&
                          &//TRIM(outtitle)
      WRITE (iunit,*)
   ENDIF
!  ---initialise parameters for harmonic analysis
   mtot = 144
   knt = -mtot
   ubsum = 0.0; ubsumc = 0.0; ubsums = 0.0
   ussum = 0.0; ussumc = 0.0; ussums = 0.0
!  ---output frequency
   icout = 360
ENDIF

!
!2. Evaluate/write parameters
!----------------------------
!

ihour = nt/120
ic = 0; ip = 0; kc = nz/2; kp = 0

IF (nt.GT.0.AND.mult_index(nt,icout)) THEN

!
!2.1 Combine at master process
!-----------------------------
!

   nhdims = nhalo
   CALL combine_mod(salglb,sal,(/1-nhalo,1-nhalo,1/),iarr_sal,0.0)
   nhdims = 1
   CALL combine_mod(deptotglb,deptotatc,(/0,0/),iarr_deptotatc,0.0)

   IF (master) THEN

!
!2.2 Position of the front
!-------------------------
!

      dhun = 1000.0
      scrit = 34.0
      i_220: DO i=1,nc-1
         IF ((ic.EQ.0).AND.(salglb(i,1,nz).GT.scrit)) ic = i
      ENDDO i_220
      IF (ic.EQ.0) THEN
         hfront = dhun*(nc-1)
      ELSEIF (ic.EQ.1) THEN
         hfront = 0.0
      ELSE
         xc = dhun*(ic-1.5)
         hfront = xc+(scrit-salglb(ic-1,1,nz))*dhun&
                  & /(salglb(ic,1,nz)-salglb(ic-1,1,nz))
      ENDIF
      hfront = 0.001*hfront

!
!2.3 Halocline depth
!-------------------
!

      scrit = 34.5
      ip = MAX(hfront-5,1.0)
      kp = 0
      k_230: DO k=nz,1,-1
         IF ((kp.EQ.0).AND.(salglb(ip,1,k).GT.scrit)) kp = k
      ENDDO k_230
      IF (kp.EQ.0) THEN
         pdep = deptotglb(ip,1)
      ELSEIF (kp.EQ.nz) then
         pdep = 0.0
      ELSE
         dep = (1.0-(kp-0.5)/REAL(nz))*deptotglb(ip,1)
         pdep = dep-(deptotglb(ip,1)/REAL(nz))*(scrit-salglb(ip,1,kp))&
                 & /(salglb(ip,1,kp+1)-salglb(ip,1,kp))
      ENDIF

!
!2.4 Maximum horizontal gradient
!-------------------------------
!

      gradmax = 0.0
      i_240: DO i=2,nc-1
         grad = ABS(salglb(i,1,kc)-salglb(i-1,1,kc))
         IF (grad.GT.gradmax) gradmax = grad
      ENDDO i_240
      hgrad = gradmax*(1000.0/dhun)

!
!2.5 Maximum vertical gradient
!-----------------------------
!

      vgrad = 0.0
      k_250: DO k=nz,2,-1
         grad = REAL(nz)*ABS(salglb(ip,1,k)-salglb(ip,1,k-1))/deptotglb(ip,1)
         IF (grad.GT.vgrad) vgrad = grad
      ENDDO k_250

!
!2.6 Write output
!----------------
!

      WRITE (iunit,9001) ihour
      WRITE (iunit,9002) 'hfront', hfront
      WRITE (iunit,9002) 'pdep', pdep
      WRITE (iunit,9002) 'hgrad', hgrad
      WRITE (iunit,9002) 'vgrad', vgrad

   ENDIF

ENDIF

!
!3. Harmonic analysis
!--------------------
!
!3.1 Update sums
!---------------
!

IF (corrstep.AND.(ihour.GE.48).AND.(ihour.LE.72)) THEN

   freq = tidal_spectrum(index_obc(1))
   ih = 70
   idproc = num_proc(ih,1) - 1

   IF (idloc.EQ.idproc) THEN

      phas = freq*knt*delt3d
      ihloc = ih - nc1loc + 1

!     ---surface
      ussum  = ussum  + uvel(ihloc,1,nz)
      ussumc = ussumc + uvel(ihloc,1,nz)*cos(phas)
      ussums = ussums + uvel(ihloc,1,nz)*sin(phas)
!     ---bottom
      ubsum  = ubsum  + uvel(ihloc,1,kc)
      ubsumc = ubsumc + uvel(ihloc,1,kc)*cos(phas)
      ubsums = ubsums + uvel(ihloc,1,kc)*sin(phas)

      knt = knt + 1

   ENDIF

ENDIF

!
!3.2 Evaluate harmonic parameters
!--------------------------------
!

IF ((ihour.EQ.72).AND.(idloc.EQ.idproc)) THEN

   sumc = sin((mtot+0.5)*freq*delt3d)/sin(0.5*freq*delt3d)
   sumcc = 0.5*sin((2*mtot+1)*freq*delt3d)/sin(freq*delt3d)&
         & +  mtot + 0.5
   sumss = 2*mtot + 1 - sumcc
   denom = (2*mtot+1)*sumcc - sumc*sumc
!  ---surface
   sures = (sumcc*ussum-sumc*ussumc)/denom
   fcos = ((2*mtot+1)*ussumc-sumc*ussum)/denom
   fsin = ussums/sumss
   suamp = SQRT(fcos**2+fsin**2)
   supha = ATAN2(fsin,fcos)
   IF (supha.LT.0.0) supha = supha + twopi
   supha = supha*radtodeg
!  ---bottom
   bures = (sumcc*ubsum-sumc*ubsumc)/denom
   fcos = ((2*mtot+1)*ubsumc-sumc*ubsum)/denom
   fsin = ubsums/sumss
   buamp = SQRT(fcos**2+fsin**2)
   bupha = ATAN2(fsin,fcos)
   IF (bupha.LT.0.0) bupha = bupha + twopi
   bupha = bupha*radtodeg

ENDIF

!
!3.3 Write output parameters
!---------------------------
!

IF (ihour.EQ.72) THEN

   CALL copy_vars(sures,0,idroot=idproc)
   CALL copy_vars(suamp,0,idroot=idproc)
   CALL copy_vars(supha,0,idroot=idproc)
   CALL copy_vars(bures,0,idroot=idproc)
   CALL copy_vars(buamp,0,idroot=idproc)
   CALL copy_vars(bupha,0,idroot=idproc)

   IF (master) THEN
      WRITE (iunit,*)
      WRITE (iunit,'(A)') 'Harmonic parameters'
      WRITE (iunit,9002) 'sures', sures
      WRITE (iunit,9002) 'suamp', suamp
      WRITE (iunit,9002) 'supha', supha
      WRITE (iunit,9002) 'bures', bures
      WRITE (iunit,9002) 'buamp', buamp
      WRITE (iunit,9002) 'bupha', bupha
   ENDIF

ENDIF

!
!4. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')

CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',I2,' hours')
9002 FORMAT (T3,A,T12,': ',G12.4E3)

END SUBROUTINE usrdef_output
