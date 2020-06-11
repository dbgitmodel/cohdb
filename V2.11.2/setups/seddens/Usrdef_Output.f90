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
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.5
!
! $Date: 2018-07-23 10:21:02 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1167 $
!
! Description - test case seddens
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_mod, mult_index, open_file
!
!************************************************************************
!
USE currents
USE fluxes
USE gridpars
USE iopars
USE modids
USE paralpars
USE physpars
USE sedarrays
USE sedids
USE sedpars
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
LOGICAL :: flag
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER :: i, k
INTEGER, SAVE :: iunit
INTEGER, DIMENSION(1) :: maxpos
REAL :: area, delx, dely, delz, hcurmax, hgradmax, rtime, volume
REAL, DIMENSION(nc) :: array1d, hcur
REAL, DIMENSION(0:nc,0:nr) :: bstresglb
REAL, DIMENSION(0:nc+1,0:nr+1,nz) :: ctotglb
REAL, DIMENSION(nc,nr,nf) :: bfluxglb
REAL, DIMENSION(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo,nz) :: uvelglb

!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
!  ---sediment initial conditions
   modfiles(io_inicon,ics_sed,1)%status = 'R'
   modfiles(io_inicon,ics_sed,1)%form = modfiles(io_fincon,ics_sed,2)%form
   modfiles(io_inicon,ics_sed,1)%filename = &
 & modfiles(io_fincon,ics_sed,2)%filename
   modfiles(io_fincon,ics_sed,2)%status = '0'
!  ---sediment particle attributes
   modfiles(io_sedspc,1,1)%status = 'R'
   modfiles(io_sedspc,1,1)%form = modfiles(io_sedspc,1,2)%form
   modfiles(io_sedspc,1,1)%filename = modfiles(io_sedspc,1,2)%filename
   modfiles(io_sedspc,1,2)%status = '0'
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (cold_start) RETURN

IF (nt.EQ.0) THEN
!  ---open output file
   IF (master) THEN
      CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case seddens: '&
                      & //'simulation '//TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF
ENDIF

!
!3. Store parameters
!-------------------
!

IF (mult_index(nt,50)) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

   rtime = nt*delt2d

!  ---combine
   CALL combine_mod(ctotglb,ctot,(/0,0,1/),iarr_ctot,0.0,commall=.TRUE.)
   CALL combine_mod(bfluxglb,bottom_sed_flux,(/1,1,1/),iarr_bottom_sed_flux,&
                  & 0.0,commall=.TRUE.)
   CALL combine_mod(uvelglb,uvel,(/1-nhalo,1-nhalo,1/),iarr_uvel,0.0,&
                  & commall=.TRUE.)
   CALL combine_mod(bstresglb,bstresatc,(/0,0/),iarr_bstresatc,0.0,&
                  & commall=.TRUE.)

!  ---grid spacings
   delx = surfacegrids(igrd_model,1)%delxdat
   dely = surfacegrids(igrd_model,1)%delydat
   delz = depmean_cst/REAL(nz)
   area = (nc-1)*delx*dely
   volume = area*depmean_cst

!  ---current height
   hcur = 0.0
   i_310: DO i=1,nc-1
      flag = .TRUE.
      k_311: DO k=1,nz
         IF (flag.AND.ctotglb(i,1,k).LT.0.01) THEN
            hcur(i) = (k-0.5)*delz
            flag = .FALSE.
         ENDIF
      ENDDO k_311
      IF (flag) hcur(i) = (nz-0.5)*delz
   ENDDO i_310
   hcurmax = MAXVAL(hcur(1:nc-1))
   maxpos = MAXLOC(hcur(1:nc-1))
   i = maxpos(1)

!  ---maximum gradients
   array1d = delz*SUM(ctotglb(1:nc,1,:),DIM=2)/depmean_cst
   hgradmax = MAXVAL(array1d(2:nc)-array1d(1:nc-1))/delx

!  ---write output
   IF (master) THEN
      WRITE (iunit,9001) rtime
      WRITE (iunit,9002) 'bstresmax', density_ref*MAXVAL(bstresglb(1:nc,1:nr))
      WRITE (iunit,9002) 'ubot', uvelglb(i,1,1)
      WRITE (iunit,9002) 'umean', SUM(uvelglb(i,1,:))/depmean_cst
      WRITE (iunit,9002) 'hcurmax', hcurmax
      WRITE (iunit,9002) 'xmax', (i-0.5)*delx
      WRITE (iunit,9002) 'hgradmax', hgradmax
      WRITE (iunit,9002) 'sedmax', rhos(1)*MAXVAL(ctotglb(1:nc-1,1,:))
      WRITE (iunit,9002) 'sedmin', rhos(1)*MINVAL(ctotglb(1:nc-1,1,:))
      WRITE (iunit,9002) 'sedtot', rhos(1)*volume*SUM(ctotglb(1:nc-1,1,:))
      WRITE (iunit,9002) 'sedbot', rhos(1)*ctotglb(i,1,1)
      WRITE (iunit,9002) 'bflxtot', rhos(1)*area*SUM(bfluxglb(1:nc-1,1,1))
   ENDIF

   CALL log_timer_out()

ENDIF

!
!4. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')


RETURN

9001 FORMAT ('Time: ',F5.2,' seconds')
9002 FORMAT (T3,A,T12,': ',G12.4E3)

END SUBROUTINE usrdef_output
