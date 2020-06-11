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
! *usrdef_output* User-formatted output (example routine)
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.5
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description - output parameters for test case sedhprof
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_mod, mult_index_, open_file
!
!************************************************************************
!
USE grid
USE gridpars
use fluxes
USE iopars
USE paralpars
USE sedarrays
USE sedids
use sedswitches
USE switches
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
INTEGER :: itime, j
INTEGER, SAVE :: iunit
REAL :: ce, ct, delh, diffe, diffrefe, diffreft, difft, e_scale, sedint, &
      & sedmax, sedmin, time_scale
REAL, DIMENSION(0:nc+1,0:nr+1,nz) :: cglb
REAL, DIMENSION(nr-1) :: cmean


!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,1)%status = 'R'
   modfiles(io_modgrd,1,1)%form = modfiles(io_modgrd,1,2)%form
   modfiles(io_modgrd,1,1)%filename = modfiles(io_modgrd,1,2)%filename
   modfiles(io_modgrd,1,2)%status = '0'
!  ---sediment open boundary conditions
   modfiles(io_sedobc,1:2,1)%status = 'R'
   modfiles(io_sedobc,1:2,1)%form = modfiles(io_sedobc,1:2,2)%form 
   modfiles(io_sedobc,1:2,1)%filename = modfiles(io_sedobc,1:2,2)%filename
   modfiles(io_sedobc,1:2,2)%status = '0'
!  ---sediment particle attributes
   modfiles(io_sedspc,1,1)%status = 'R'
   modfiles(io_sedspc,1,1)%form = modfiles(io_sedspc,1,2)%form
   modfiles(io_sedspc,1,1)%filename = modfiles(io_sedspc,1,2)%filename
   modfiles(io_sedspc,1,2)%status = '0'
ENDIF

IF (cold_start) RETURN

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN
!  ---open output file
   IF (master.AND.(.NOT.cold_start)) THEN
      CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case sedhprof: '&
                        & //'simulation '//TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF
ENDIF

!
!3. Store parameters
!-------------------
!

IF (mult_index(nt,30)) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

   itime = nosecsrun/60

!  ---parameters
   CALL combine_mod(cglb,ctot,(/0,0,1/),iarr_ctot,0.0,commall=.TRUE.)
   IF (iopt_grid_nodim.EQ.3) THEN
      cmean = SUM(cglb(1,1:nr-1,:),DIM=2)/REAL(nz)
   ELSE
      cmean = cglb(1,1:nr-1,1)
   ENDIF
   sedmax = MAXVAL(cmean)
   sedmin = MINVAL(cmean)
   delh = delyatc(1,1)
   sedint = delh*SUM(cmean)
   ce = cmean(1); ct = ce
   j_310: DO j=1,nr-1   
      diffrefe = 2.7183*(sedmax-sedmin) - cmean(j)
      diffe = 2.7183*(sedmax-sedmin) - ce
      diffreft = 0.99*sedmax - cmean(j)
      difft = 0.99*sedmax - ct
      e_scale = 0.0; time_scale = 0.0
      IF (ABS(diffreft).LT.ABS(difft)) THEN
         ct = cmean(j)
         time_scale = j*40 - 20
      ENDIF
      IF (ABS(diffrefe).LE.ABS(diffe)) THEN
         ce = cmean(j)
         e_scale = j*40 - 20
      ENDIF
   ENDDO j_310

!  ---write output
   IF (master) THEN
      WRITE (iunit,9001) itime
      IF (iopt_sed_type.EQ.1) THEN
         WRITE (iunit,9002) 'cheight', height_c(1,1,1)
         WRITE (iunit,9002) 'cref', cref(1,1,1)
         WRITE (iunit,9002) 'ceq', ceq(1,1,1)
      ENDIF
      WRITE (iunit,9002) 'sedmin', sedmin
      WRITE (iunit,9002) 'sedmax', sedmax
      WRITE (iunit,9002) 'sedint', sedint
      WRITE (iunit,9002) 't_scale', time_scale
      WRITE (iunit,9002) 'e_scale', e_scale
   ENDIF

   CALL log_timer_out()

ENDIF

!
!4. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')


RETURN

9001 FORMAT ('Time: ',I3,' minutes')
9002 FORMAT (T3,A,T12,': ',G12.4E3)

END SUBROUTINE usrdef_output
