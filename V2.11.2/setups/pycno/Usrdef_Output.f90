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
! $Date: 2018-07-23 10:21:02 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1167 $
!
! Description - output parameters for test case pycno
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - mixed_layer
!
! Module calls - close_file, least_squares_fit, mult_index, open_file
!
!************************************************************************
!
USE currents
USE density
USE diffusion
USE gridpars
USE iopars
USE physpars
USE syspars
USE timepars
USE turbulence
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: least_squares_fit, mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER, PARAMETER :: ndatmax = 200
INTEGER, SAVE :: iunit, nodat
REAL, SAVE :: bf0, rvmean, stotin, s0sur
REAL :: adr0, ah, arv, au2, bdr0, bh, brv, bu2, corrdr0, corrh, corru2, &
      & corrv, difmax, drosur, riv, rvmean2, sdev, stot, tkemax, &
      & tmld, vedmax
REAL, SAVE, DIMENSION(ndatmax) :: drodat, hdat, rivdat, udat, xdat


!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
!  ---initial conditions
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = modfiles(io_fincon,ics_phys,2)%form 
   modfiles(io_inicon,ics_phys,1)%filename = &
 & modfiles(io_fincon,ics_phys,2)%filename
   modfiles(io_fincon,ics_phys,2)%status = '0'
!  ---meteo forcing
   modfiles(io_metsur,1,1)%status = 'R'
   modfiles(io_metsur,1,1)%form = 'A'
   modfiles(io_metsur,1,1)%filename = modfiles(io_metsur,1,2)%filename
   modfiles(io_metsur,1,2)%status = '0'
ENDIF

IF (cold_start) RETURN

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN
!  ---open output file
   CALL open_file(iunit,runtitle(1:6)//'.'//suffix,'OUT','A')
   WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
   WRITE (iunit,'(A)') 'Output parameters for test case pycno: simulation '//&
                      & TRIM(runtitle)
   WRITE (iunit,*)
!  ---initialise parameters
   nodat = 0
   bf0 = 1.0E-02
   rvmean = 0.0
   s0sur = sal(1,1,nz)
   stotin = SUM(sal(1,1,:))/REAL(nz)
   drodat = 0.0; hdat = 0.0; rivdat = 0.0; udat = 0.0; xdat = 0.0
ENDIF

!
!3. Update parameters
!-------------------
!

IF (nt.GE.600.AND.mult_index(nt,60)) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

   nodat = nodat + 1
   CALL mixed_layer(tmld,riv)
   rvmean = rvmean + riv
   xdat(nodat) = LOG(nosecsrun*bf0)
   rivdat(nodat) = LOG(riv)
   hdat(nodat) = LOG(tmld)
   udat(nodat) = LOG(uvel(1,1,nz))
   drosur = density_ref*beta_sal_ref*(sal(1,1,nz)-s0sur)
   drodat(nodat) = LOG(drosur)

ENDIF

IF (nt.EQ.nstep) THEN
!  ---maximum values
   vedmax = MAXVAL(vdifcoefmom(1,1,:))
   difmax = MAXVAL(vdifcoefscal(1,1,:))
   tkemax = MAXVAL(tke(1,1,:))
!  ---mean salinity deviation
   stot = SUM(sal(1,1,:))/REAL(nz)
   sdev = stot/stotin - 1.0
ENDIF

!
!4. Output values at final step
!------------------------------
!

IF (nt.EQ.nstep) THEN

   WRITE (iunit,9001) nosecsrun/3600
   rvmean = rvmean/REAL(nodat)
   CALL least_squares_fit(xdat(1:nodat),rivdat(1:nodat),nodat,arv,brv,corrv)
   WRITE (iunit,9002) 'rvmean', rvmean
   WRITE (iunit,9002) 'brv', brv
   WRITE (iunit,9002) 'corrv', corrv

   CALL least_squares_fit(xdat(1:nodat),hdat(1:nodat),nodat,ah,bh,corrh)
   rvmean2 = 0.5*(EXP(ah))**4
   WRITE (iunit,9002) 'rvmean2', rvmean2
   WRITE (iunit,9002) 'bh', bh
   WRITE (iunit,9002) 'corrh', corrh

   CALL least_squares_fit(xdat(1:nodat),udat(1:nodat),nodat,au2,bu2,corru2)
   WRITE (iunit,9002) 'bu2', bu2
   WRITE (iunit,9002) 'corru2', corru2

   CALL least_squares_fit(xdat(1:nodat),drodat(1:nodat),nodat,adr0,bdr0,&
                       & corrdr0)
   WRITE (iunit,9002) 'bdr0', bdr0
   WRITE (iunit,9002) 'corrdr0', corrdr0

   WRITE (iunit,9002) 'vedmax', 100*vedmax
   WRITE (iunit,9002) 'difmax', 100*difmax
   WRITE (iunit,9002) 'tkemax', 10000*tkemax

   WRITE (iunit,9002) 'tmld', tmld
   WRITE (iunit,9002) 'usur', uvel(1,1,nz)
   WRITE (iunit,9002) 'drosur', drosur

   WRITE (iunit,9002) 'sdev', 1.0E+05*sdev

ENDIF

!
!5. Close output file
!--------------------
!

IF (nt.EQ.nstep) CALL close_file(iunit,'A')

IF (nt.GE.600.AND.mult_index(nt,60)) THEN
   CALL log_timer_out()
ENDIF


RETURN

9001 FORMAT ('Time: ',I2,' hours') 
9002 FORMAT (T3,A,T12,': ',G12.4E3)

END SUBROUTINE usrdef_output
