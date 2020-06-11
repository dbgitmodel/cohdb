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
! Author - Alexander Breugem Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.10
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description - test case morphctb
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
USE morphswitches
USE paralpars
USE physpars
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
INTEGER, PARAMETER :: nostats = 6
INTEGER :: i, l
INTEGER, SAVE :: icount, iunit
INTEGER, SAVE, DIMENSION(nostats) :: istat
REAL :: delh, rhour
REAL, SAVE :: spacetimefac1, spacetimefac2
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: bstresatc_glb, bed_update_int_glb, &
                                 & depmeanatc_glb, depmean_ini_glb, &
                                 & depero_src, qload_src, sed_avail_tot_glb
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: bed_update_dep_glb, &
                                 & bed_update_dep_int_glb, bed_update_ero_glb, &
                                 & bed_update_ero_int_glb, bottom_sed_flux_glb,&
                                 & ctot_glb, qloadatu_glb


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
   modfiles(io_2uvobc,1,1)%status = 'R'
   modfiles(io_2uvobc,1,1)%form = modfiles(io_2uvobc,1,2)%form
   modfiles(io_2uvobc,1,1)%filename = modfiles(io_2uvobc,1,2)%filename
   modfiles(io_2uvobc,1,2)%status = '0'   
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

IF (cold_start.OR.runtitle(9:9).EQ.'0') RETURN

sflag = iopt_sed_mode.EQ.3

IF (nt.EQ.0) THEN

!  ---open output file
   IF (master) THEN
      CALL open_file(iunit,TRIM(outtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case morphtcb: '&
                      & //'simulation '//TRIM(outtitle)
      WRITE (iunit,*)
   ENDIF

!  ---station indices
   istat = (/5,15,25,35,45,55/)

!  ---allocate arrays
   ALLOCATE (bed_update_dep_glb(nc,nr,nf),STAT=errstat)
   CALL error_alloc('bed_update_dep_glb',3,(/nc,nr,nf/),real_type)
   ALLOCATE (bed_update_dep_int_glb(nc,nr,nf),STAT=errstat)
   CALL error_alloc('bed_update_dep_int_glb',3,(/nc,nr,nf/),real_type)
   ALLOCATE (bed_update_ero_glb(nc,nr,nf),STAT=errstat)
   CALL error_alloc('bed_update_ero_glb',3,(/nc,nr,nf/),real_type)
   ALLOCATE (bed_update_ero_int_glb(nc,nr,nf),STAT=errstat)
   CALL error_alloc('bed_update_ero_int_glb',3,(/nc,nr,nf/),real_type)
   ALLOCATE (bed_update_int_glb(nc,nr),STAT=errstat)
   CALL error_alloc('bed_update_int_glb',2,(/nc,nr/),real_type)
   ALLOCATE (bstresatc_glb(0:nc,0:nr),STAT=errstat)
   CALL error_alloc('bstresatc_glb',2,(/nc+1,nr+1/),real_type)
   ALLOCATE (depmeanatc_glb(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo),STAT=errstat)
   CALL error_alloc('depmeanatc_glb',2,(/nc+2*nhalo,nr+2*nhalo/),real_type)
   ALLOCATE (depmean_ini_glb(nc,nr),STAT=errstat)
   CALL error_alloc('depmean_ini_glb',2,(/nc,nr/),real_type)
   ALLOCATE (qloadatu_glb(0:nc+1,nr,nf),STAT=errstat)
   CALL error_alloc('qloadatu_glb',3,(/nc+2,nr,nf/),real_type)
   ALLOCATE (qload_src(nc,nr),STAT=errstat)
   CALL error_alloc('qload_src',2,(/nc,nr/),real_type)
   IF (iopt_morph_fixed_layer.EQ.1) THEN
      ALLOCATE (sed_avail_tot_glb(nc,nr),STAT=errstat)
      CALL error_alloc('sed_avail_tot_glb',2,(/nc,nr/),real_type)
   ENDIF
   IF (sflag) THEN
      ALLOCATE (bottom_sed_flux_glb(nc,nr,nf),STAT=errstat)
      CALL error_alloc('bottom_sed_flux_glb',3,(/nc,nr,nf/),real_type)
      ALLOCATE (ctot_glb(0:nc+1,0:nr+1,nz),STAT=errstat)
      CALL error_alloc('ctot_glb',3,(/nc+1,nr+1,nz/),real_type)
      ALLOCATE (depero_src(nc,nr),STAT=errstat)
      CALL error_alloc('depero_src',2,(/nc,nr/),real_type)
   ENDIF

!  ---initialise
   bed_update_dep_int_glb = 0.0
   bed_update_ero_int_glb = 0.0
   depmean_ini_glb = depmeanglb(1:nc,1:nr)
   IF (sflag) depero_src = 0.0
   qload_src = 0.0

!  ---initialise parameters
   delh = surfacegrids(igrd_model,1)%delxdat
   spacetimefac1 =  (morph_factor*deltmorph)/(1.0-bed_porosity_cst)
   spacetimefac2 = spacetimefac1/delh
   icount = 720

ENDIF

!
!3. Write output parameters
!--------------------------
!

procname(pglev+1) = 'usrdef_output'
CALL log_timer_in()

!
!3.1 Combine into global arrays
!------------------------------
!

IF (nt.GT.0.AND.mult_index(nt,icount)) THEN

   CALL combine_mod(bstresatc_glb,bstresatc,(/0,0/),iarr_bstresatc,0.0,&
                  & commall=.TRUE.)
   CALL combine_mod(bed_update_int_glb,bed_update_int,(/1,1/),&
                  & iarr_bed_update_int,0.0,commall=.TRUE.)
   CALL combine_mod(depmeanatc_glb,depmeanatc,(/1-nhalo,1-nhalo/),&
                  & iarr_depmeanatc,depmean_flag,commall=.TRUE.)
   IF (iopt_morph_fixed_layer.EQ.1) THEN
      CALL combine_mod(sed_avail_tot_glb,sed_avail_tot,(/1,1/),&
                     & iarr_sed_avail_tot,0.0,commall=.TRUE.)
   ENDIF
   IF (sflag) THEN
      CALL combine_mod(ctot_glb,ctot,(/0,0,1/),iarr_ctot,0.0,commall=.TRUE.)
   ENDIF
ENDIF

IF (morphstep) THEN
   CALL combine_mod(bed_update_dep_glb,bed_update_dep,(/1,1,1/),&
                  & iarr_bed_update_dep,0.0,commall=.TRUE.)
   CALL combine_mod(bed_update_ero_glb,bed_update_ero,(/1,1,1/),&
                  & iarr_bed_update_ero,0.0,commall=.TRUE.)
   IF (sflag) THEN
      CALL combine_mod(bottom_sed_flux_glb,bottom_sed_flux,(/1,1,1/),&
                     & iarr_bottom_sed_flux,0.0,commall=.TRUE.)
      CALL combine_mod(qloadatu_glb,qbedatu,(/0,1,1/),iarr_qbedatu,0.0,&
                     & commall=.TRUE.)
   ELSE
      CALL combine_mod(qloadatu_glb,qtotatu,(/0,1,1/),iarr_qtotatu,0.0,&
                     & commall=.TRUE.)
   ENDIF
ENDIF

!
!3.2 Update deposition/erosion bed updates
!-----------------------------------------
!

IF (morphstep) THEN
   i_320: DO i=1,nc-1
      bed_update_dep_int_glb(i,1,:) = bed_update_dep_int_glb(i,1,:) + &
                                    & bed_update_dep_glb(i,1,:)
      bed_update_ero_int_glb(i,1,:) = bed_update_ero_int_glb(i,1,:) + &
                                    & bed_update_ero_glb(i,1,:)
      qload_src(i,1) = qload_src(i,1) + &
                    & (qloadatu_glb(i,1,1)-qloadatu_glb(i+1,1,1))*spacetimefac2
      IF (sflag) THEN
         depero_src(i,1) = depero_src(i,1) - &
                         & bottom_sed_flux_glb(i,1,1)*spacetimefac1
      ENDIF
   ENDDO i_320
ENDIF


!
!3.3 Write data
!--------------
!

IF (nt.GT.0.AND.mult_index(nt,icount).AND.master) THEN

   rhour = nosecsrun/3600.0
   WRITE (iunit,9001) rhour

   l_320: DO l=1,nostats
      i = istat(l)
      WRITE (iunit,9003) 'Station', l
      WRITE (cl,'(I12)') l; cl = ADJUSTL(cl)
      WRITE (iunit,9002) 'depmeanatc_'//TRIM(cl), depmeanatc_glb(i,1)
      WRITE (iunit,9002) 'bstresatc_'//TRIM(cl), bstresatc_glb(i,1)
      WRITE (iunit,9002) 'qloadatu_'//TRIM(cl), rhos(1)*qloadatu_glb(i,1,1)
      WRITE (iunit,9002) 'bed_update_int_'//TRIM(cl), bed_update_int_glb(i,1)
      WRITE (iunit,9002) 'bed_update_dep_int_'//TRIM(cl), &
                        & bed_update_dep_int_glb(i,1,1)
      WRITE (iunit,9002) 'bed_update_ero_int_'//TRIM(cl), &
                        & bed_update_ero_int_glb(i,1,1)
      IF (sflag) THEN
         WRITE (iunit,9002) 'ctot_'//TRIM(cl), ctot_glb(i,1,1)
         WRITE (iunit,9002) 'depero_src_'//TRIM(cl), depero_src(i,1)
      ENDIF
      WRITE (iunit,9002) 'qload_src_'//TRIM(cl), qload_src(i,1)
      IF (iopt_morph_fixed_layer.EQ.1) THEN
         WRITE (iunit,9002) 'sed_avail_tot_'//TRIM(cl), sed_avail_tot_glb(i,1)
      ENDIF
      WRITE (iunit,9002) 'vert_sed_balance_'//TRIM(cl), &
          & depmeanatc_glb(i,1) - depmean_ini_glb(i,1) + bed_update_int_glb(i,1)
   ENDDO l_320

   IF (nt.LT.nstep) WRITE (iunit,*)

ENDIF

CALL log_timer_out()

!
!3. Finalise
!-----------
!

IF (nt.EQ.nstep) THEN

!  ---close output file
   IF (master) CALL close_file(iunit,'A')

!  ---deallocate
   DEALLOCATE (bed_update_dep_glb,bed_update_dep_int_glb,bed_update_ero_glb,&
             & bed_update_ero_int_glb,bed_update_int_glb,&
             & bstresatc_glb,depmeanatc_glb,depmean_ini_glb,&
             & qloadatu_glb,qload_src)
   IF (iopt_morph_fixed_layer.EQ.1) DEALLOCATE (sed_avail_tot_glb)
   IF (sflag) DEALLOCATE (bottom_sed_flux_glb,ctot_glb,depero_src)

ENDIF


RETURN

9001 FORMAT ('Time: ',F4.1,' hours')
9002 FORMAT (T3,A,T27,': ',G12.4E3)
9003 FORMAT (A,': ',I3)

END SUBROUTINE usrdef_output
