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
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description - output parameters for test case seich
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_mod, maxloc_vars, minloc_vars,
!                mult_index, open_file, sum2_vars
!
!************************************************************************
!
USE currents
USE density
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE syspars
USE timepars
USE inout_routines, ONLY: close_file, open_file
USE paral_comms, ONLY: combine_mod
USE paral_utilities, ONLY: maxloc_vars, minloc_vars, sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER :: icount, j, jleft, jright, k, klo, kup
INTEGER, SAVE :: iunit
REAL :: dhun, hleft, hlen, hour, hright, sdev, stot, vmean, volume, &
      & v2lo, v2up, wmaxlo, wmaxup, wminlo, wminup
REAL, SAVE :: stotin
INTEGER, DIMENSION(3) :: lpos
INTEGER, DIMENSION(4) :: nhdims = 0
REAL, DIMENSION(ncloc,nrloc,nz) :: array3d
REAL, DIMENSION(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo,nz) :: vvelglb


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
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN
!  ---open output file
   IF (.NOT.cold_start.AND.master) THEN
      CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case seich: simulation '&
                         & //TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF
!  ---initial value of volume averaged salinity
   CALL sum2_vars(delzatc(1:ncloc,1:nrloc,:),volume,nhdims,'C  ',0)
   array3d = delzatc(1:ncloc,1:nrloc,:)*sal(1:ncloc,1:nrloc,:)
   CALL sum2_vars(array3d,stotin,nhdims,'C  ',0)
   IF (master) stotin = stotin/volume
ENDIF

!
!3. Store parameters
!-------------------
!

icount = 240
hour = nosecsrun/3600.0

IF ((hour.GT.0.0).AND.(mult_index(nt,icount)).AND.(hour.LE.6)) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

   CALL combine_mod(vvelglb,vvel,(/1-nhalo,1-nhalo,1/),iarr_vvel,0.0)

!
!3.1 upper layer
!---------------
!

   dhun = 500.0
   hlen = dhun*(nr-1)
   kup = 15

!  ---vertical current
   CALL minloc_vars(wphys(:,:,kup),wminup,lpos,iarr_wphys,mask=maskatc_int)
   jleft = lpos(2)
   CALL maxloc_vars(wphys(:,:,kup),wmaxup,lpos,iarr_wphys,mask=maskatc_int)
   jright = lpos(2)

!  ---location of the front
   IF (master) hleft = (0.5*hlen-(jleft-0.5)*dhun)/1000.0

!  ---horizontal current
   IF (master) THEN
      vmean = 0.0
      j_310: DO j=jleft,jright
         vmean = vmean + vvelglb(1,j,kup)
      ENDDO j_310
      v2up = vmean/(jright-jleft+1)
   ENDIF

!
!3.2 lower layer
!---------------
!

   klo = 5

!  ---vertical current
   CALL minloc_vars(wphys(:,:,klo),wminlo,lpos,iarr_wphys,mask=maskatc_int)
   jleft = lpos(2)
   CALL maxloc_vars(wphys(:,:,klo),wmaxlo,lpos,iarr_wphys,mask=maskatc_int)
   jright = lpos(2)

!  ---location of the front
   IF (master) hright = ((jright-0.5)*dhun-0.5*hlen)/1000.0

!  ---horizontal current
   IF (master) THEN
      vmean = 0.0
      j_320: DO j=jleft,jright
         vmean = vmean + vvelglb(1,j,klo)
      ENDDO j_320
      v2lo = vmean/(jright-jleft+1)
   ENDIF

!
!3.3 Averaged salinity
!---------------------
!

   CALL sum2_vars(delzatc(1:ncloc,1:nrloc,:),volume,nhdims,'C  ',0)
   array3d = delzatc(1:ncloc,1:nrloc,:)*sal(1:ncloc,1:nrloc,:)
   CALL sum2_vars(array3d,stot,nhdims,'C  ',0)
   IF (master) THEN
      sdev = (stot/volume-stotin)/stotin
   ENDIF

!
!3.4 Write output
!----------------
!

   IF (master) THEN
      WRITE (iunit,9001) hour
      WRITE (iunit,9002) 'hleft', hleft
      WRITE (iunit,9002) 'hright', hright
      WRITE (iunit,9002) 'v2up', v2up
      WRITE (iunit,9002) 'v2lo', v2lo
      WRITE (iunit,9002) 'wmaxup', wmaxup*1000.0
      WRITE (iunit,9002) 'wmaxlo', wmaxlo*1000.0
      WRITE (iunit,9002) 'wminup', wminup*1000.0
      WRITE (iunit,9002) 'wminlo', wminlo*1000.0
      WRITE (iunit,9002) 'sdev', sdev*1.0E+07
   ENDIF

   CALL log_timer_out()

ENDIF

!
!4. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')


RETURN

9001 FORMAT ('Time: ',F4.2,' hours')
9002 FORMAT (T3,A,T12,': ',G12.4E3)

END SUBROUTINE usrdef_output
