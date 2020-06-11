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
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.10
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description - test case morphtrench
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_mod, error_alloc, mult_index, open_file
!
!************************************************************************
!
USE depths
USE fluxes
USE gridpars
USE iopars
USE modids
USE morpharrays
USE morphpars
USE paralpars
USE sedarrays
USE sedids
USE sedpars
USE sedswitches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: close_file, open_file
USE paral_comms, ONLY: combine_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: sflag
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
CHARACTER (LEN=12) :: cl
INTEGER, PARAMETER :: nostats = 9
INTEGER :: i, j, l
INTEGER, SAVE :: icount, iunit
INTEGER, SAVE, DIMENSION(nostats) :: istat, jstat
REAL :: delh, rhour
REAL, SAVE :: spacetimefac1, spacetimefac2
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: bed_update_int_glb, bstresatc_glb, &
                                         & depmean_ini_glb, depmeanatc_glb, &
                                         & depero_src, qbed_src
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: bottom_sed_flux_glb, ctot_glb, &
                                           & qbedatu_glb, qbedatv_glb


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
!  ---open boundary conditions (2-D)
   modfiles(io_2uvobc,1:2,1)%status = 'R'
   modfiles(io_2uvobc,1:2,1)%form = modfiles(io_2uvobc,1:2,2)%form
   modfiles(io_2uvobc,1:2,1)%filename = &
 & modfiles(io_2uvobc,1:2,2)%filename
   modfiles(io_2uvobc,1:2,2)%status = '0'
!  ---sediment particle attributes
   modfiles(io_sedspc,1,1)%status = 'R'
   modfiles(io_sedspc,1,1)%form = 'N'
   modfiles(io_sedspc,1,1)%filename = modfiles(io_sedspc,1,2)%filename
   modfiles(io_sedspc,1,2)%status = '0'
!  ---initial conditions (morphology)
   IF (runtitle(9:9).NE.'0') THEN
      modfiles(io_inicon,ics_morph,1)%status = 'R'
      modfiles(io_inicon,ics_morph,1)%form = 'N'
      modfiles(io_inicon,ics_morph,1)%filename = 'morphics.nc'
      modfiles(io_fincon,ics_morph,2)%status = '0'
   ENDIF
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (cold_start.OR.runtitle(12:12).EQ.'0') RETURN

sflag = iopt_sed_mode.EQ.3

IF (nt.EQ.0) THEN

!  ---open output file
   IF (master) THEN
      CALL open_file(iunit,TRIM(outtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case morphtrench: '&
                      & //'simulation '//TRIM(outtitle)
      WRITE (iunit,*)
   ENDIF

!  ---station indices
   IF (LLT(runtitle(13:13),'E')) THEN
      istat = (/6,11,16,6,11,16,6,11,16/)
      jstat = (/21,21,21,25,25,25,30,30,30/)
   ELSE
      istat = (/2,2,2,3,3,3,4,4,4/)
      jstat = (/31,31,31,35,35,35,40,40,40/)
   ENDIF

!  ---allocate arrays
   ALLOCATE (bed_update_int_glb(nc,nr),STAT=errstat)
   CALL error_alloc('bed_update_int__glb',2,(/nc,nr/),real_type)
   ALLOCATE (bstresatc_glb(0:nc,0:nr),STAT=errstat)
   CALL error_alloc('bstresatc_glb',2,(/nc+1,nr+1/),real_type)
   ALLOCATE (depmean_ini_glb(nc,nr),STAT=errstat)
   CALL error_alloc('depmean_ini_glb',2,(/nc,nr/),real_type)
   ALLOCATE (depmeanatc_glb(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo),STAT=errstat)
   CALL error_alloc('depmeanatc_glb',2,(/nc+2*nhalo,nr+2*nhalo/),real_type)
   ALLOCATE (qbedatu_glb(0:nc+1,nr,nf),STAT=errstat)
   CALL error_alloc('qbedatu_glb',3,(/nc+2,nr,nf/),real_type)
   ALLOCATE (qbedatv_glb(nc,0:nr+1,nf),STAT=errstat)
   CALL error_alloc('qbedatv_glb',3,(/nc,nr+2,nf/),real_type)
   ALLOCATE (qbed_src(nc,nr),STAT=errstat)
   CALL error_alloc('qbed_src',2,(/nc,nr/),real_type)
   IF (sflag) THEN
      ALLOCATE (bottom_sed_flux_glb(nc,nr,nf),STAT=errstat)
      CALL error_alloc('bottom_sed_flux_glb',3,(/nc,nr,nf/),real_type)
      ALLOCATE (ctot_glb(0:nc+1,0:nr+1,nz),STAT=errstat)
      CALL error_alloc('ctot_glb',3,(/nc+1,nr+1,nz/),real_type)
      ALLOCATE (depero_src(nc,nr),STAT=errstat)
      CALL error_alloc('depero_src',2,(/nc,nr/),real_type)
   ENDIF

!  ---initialise
   IF (sflag) depero_src = 0.0
   qbed_src = 0.0
   depmean_ini_glb = depmeanglb(1:nc,1:nr)

!  ---initialise parameters
   delh = surfacegrids(igrd_model,1)%delxdat
   spacetimefac1 =  (morph_factor*deltmorph)/(1.0-bed_porosity_cst)
   spacetimefac2 = spacetimefac1/delh
   icount = 43200

ENDIF

!
!3. Write output parameters
!--------------------------
!

IF (morphstep) THEN
   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()
ENDIF

!
!3.1 Combine into global arrays
 !------------------------------
!

IF (mult_index(nt,icount)) THEN
   CALL combine_mod(bed_update_int_glb,bed_update_int,(/1,1/),&
                  & iarr_bed_update_int,0.0,commall=.TRUE.)
   CALL combine_mod(bstresatc_glb,bstresatc,(/0,0/),iarr_bstresatc,0.0,&
                  & commall=.TRUE.)
   CALL combine_mod(depmeanatc_glb,depmeanatc,(/1-nhalo,1-nhalo/),&
                  & iarr_depmeanatc,0.0,commall=.TRUE.)
ENDIF

IF (morphstep) THEN
   IF (sflag) THEN
      CALL combine_mod(bottom_sed_flux_glb,bottom_sed_flux,(/1,1,1/),&
                     & iarr_bottom_sed_flux,0.0,commall=.TRUE.)
      CALL combine_mod(ctot_glb,ctot,(/0,0,1/),iarr_ctot,0.0,commall=.TRUE.)
   ENDIF
   CALL combine_mod(qbedatu_glb,qbedatu,(/0,1,1/),iarr_qbedatu,0.0,&
                  & commall=.TRUE.)
   CALL combine_mod(qbedatv_glb,qbedatv,(/1,0,1/),iarr_qbedatv,0.0,&
                  & commall=.TRUE.)
ENDIF

!
!3.2 Update deposition/erosion bed updates
!-----------------------------------------
!

IF (morphstep) THEN
   IF (sflag) THEN
      depero_src(1:nc-1,1:nr-1) = depero_src(1:nc-1,1:nr-1) - &
                           & bottom_sed_flux_glb(1:nc-1,1:nr-1,1)*spacetimefac1
   ENDIF
   qbed_src(1:nc-1,1:nr-1) = qbed_src(1:nc-1,1:nr-1) + &
       & (qbedatu_glb(1:nc-1,1:nr-1,1)-qbedatu_glb(2:nc,1:nr-1,1)+&
       &  qbedatv_glb(1:nc-1,1:nr-1,1)-qbedatv_glb(1:nc-1,2:nr,1))*spacetimefac2
ENDIF

!
!3.3 Write data
!--------------
!

IF (mult_index(nt,icount).AND.master) THEN

   rhour = nosecsrun/3600.0
   WRITE (iunit,9001) rhour
   WRITE (iunit,*)

   l_320: DO l=1,nostats
      i = istat(l); j = jstat(l)
      WRITE (iunit,9003) 'Station', l
      WRITE (cl,'(I12)') l; cl = ADJUSTL(cl)
      WRITE (iunit,9002) 'depmeanatc_'//TRIM(cl), depmeanatc_glb(i,j)
      WRITE (iunit,9002) 'bstresatc_'//TRIM(cl), bstresatc_glb(i,j)
      WRITE (iunit,9002) 'qbedatu_'//TRIM(cl), rhos(1)*qbedatu_glb(i,j,1)
      WRITE (iunit,9002) 'qbedatv_'//TRIM(cl), rhos(1)*qbedatv_glb(i,j,1)
      WRITE (iunit,9002) 'bed_update_int_'//TRIM(cl), bed_update_int_glb(i,j)
      IF (sflag) THEN
         WRITE (iunit,9002) 'ctot_'//TRIM(cl), ctot_glb(i,j,1)
         WRITE (iunit,9002) 'depero_src_'//TRIM(cl), depero_src(i,j)
      ENDIF
      WRITE (iunit,9002) 'qbed_src_'//TRIM(cl), qbed_src(i,j)
      WRITE (iunit,9002) 'vert_sed_balance_'//TRIM(cl), &
         & depmeanatc_glb(i,j) - depmean_ini_glb(i,j) + bed_update_int_glb(i,j)
   ENDDO l_320

   IF (nt.LT.nstep) WRITE (iunit,*)

ENDIF

IF (morphstep) CALL log_timer_out()

!
!3. Finalise
!-----------
!

IF (nt.EQ.nstep) THEN

!  ---close output file
   IF (master) CALL close_file(iunit,'A')

!  ---deallocate
   DEALLOCATE (bed_update_int_glb,bstresatc_glb,depmean_ini_glb,depmeanatc_glb,&
             & qbedatu_glb,qbedatv_glb,qbed_src)
   IF (sflag) THEN
      DEALLOCATE (bottom_sed_flux_glb,ctot_glb,depero_src)
   ENDIF

ENDIF


RETURN

9001 FORMAT ('Time: ',F4.1,' hours')
9002 FORMAT (T3,A,T27,': ',G12.4E3)
9003 FORMAT (A,': ',I3)

END SUBROUTINE usrdef_output
