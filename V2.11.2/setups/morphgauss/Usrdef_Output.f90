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
! Version - @(COHERENS)Usrdef_Output.f90  V2.10
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description - test case morphgauss
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
USE paralpars
USE morpharrays
USE morphpars
USE morphswitches
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
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
CHARACTER (LEN=12) :: cl
INTEGER, PARAMETER :: nostats = 5
INTEGER :: j, l
INTEGER, SAVE :: icount, iunit
INTEGER, SAVE, DIMENSION(nostats) :: jstat
REAL :: rhour
REAL, DIMENSION(0:nc,0:nr) :: bstresatc_glb
REAL, DIMENSION(nc,nr) :: bed_update_int_glb, sed_avail_tot_glb
REAL, DIMENSION(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo) :: depmeanatc_glb
REAL, DIMENSION(nc,nr,nf) :: bed_update_dep_glb, bed_update_ero_glb
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: depmean_ini_glb
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: bed_update_dep_int_glb, &
                                 & bed_update_ero_int_glb
REAL, DIMENSION(nc,0:nr+1,nf) :: qtotatv_glb


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
   modfiles(io_2uvobc,1:2,1)%filename = modfiles(io_2uvobc,1:2,2)%filename
   modfiles(io_2uvobc,1:2,2)%status = '0'
!  ---initial conditions (physics)
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = &
 & modfiles(io_fincon,ics_phys,2)%form
   modfiles(io_inicon,ics_phys,1)%filename = &
 & modfiles(io_fincon,ics_phys,2)%filename
   modfiles(io_fincon,ics_phys,2)%status = '0'
!  ---initial conditions (morphology)
   modfiles(io_inicon,ics_morph,1)%status = 'R'
   modfiles(io_inicon,ics_morph,1)%form = &
 & modfiles(io_fincon,ics_morph,2)%form
   modfiles(io_inicon,ics_morph,1)%filename = &
 & modfiles(io_fincon,ics_morph,2)%filename
   modfiles(io_fincon,ics_morph,2)%status = '0'
!  ---sediment particle attributes
   modfiles(io_sedspc,1,1)%status = 'R'
   modfiles(io_sedspc,1,1)%form = 'N'
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
      WRITE (iunit,'(A)') 'Output parameters for test case morphgauss: '&
                      & //'simulation '//TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF

!  ---station indices
   jstat = (/10,25,30,35,50/)

!  ---allocate arrays
   ALLOCATE (bed_update_dep_int_glb(nc,nr,nf),STAT=errstat)
   CALL error_alloc('bed_update_dep_int_glb',3,(/nc,nr,nf/),real_type)
   ALLOCATE (bed_update_ero_int_glb(nc,nr,nf),STAT=errstat)
   CALL error_alloc('bed_update_ero_int_glb',3,(/nc,nr,nf/),real_type)
   ALLOCATE (depmean_ini_glb(nc,nr),STAT=errstat)
   CALL error_alloc('depmean_ini_glb',2,(/nc,nr/),real_type)

!  ---initialise
   bed_update_dep_int_glb = 0.0
   bed_update_ero_int_glb = 0.0
   depmean_ini_glb = depmeanglb(1:nc,1:nr)

!  ---time counter
   icount = 2160

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

   CALL combine_mod(depmeanatc_glb,depmeanatc,(/1-nhalo,1-nhalo/),&
                  & iarr_depmeanatc,0.0,commall=.TRUE.)
   CALL combine_mod(bstresatc_glb,bstresatc,(/0,0/),iarr_bstresatc,0.0,&
                  & commall=.TRUE.)
   CALL combine_mod(qtotatv_glb,qtotatv,(/1,0,1/),iarr_qbedatv,0.0,&
                  & commall=.TRUE.)
   CALL combine_mod(bed_update_int_glb,bed_update_int,(/1,1/),&
                  & iarr_bed_update_int,0.0,commall=.TRUE.)
   CALL combine_mod(depmeanglb,depmeanatc,(/1-nhalo,1-nhalo/),&
                  & iarr_depmeanatc,depmean_flag,commall=.TRUE.)
   IF (iopt_morph_fixed_layer.EQ.1) THEN
      CALL combine_mod(sed_avail_tot_glb,sed_avail_tot,(/1,1/),&
                    & iarr_sed_avail_tot,0.0,commall=.TRUE.)
   ENDIF

ENDIF

IF (morphstep) THEN
   CALL combine_mod(bed_update_dep_glb,bed_update_dep,(/1,1,1/),&
                  & iarr_bed_update_dep,0.0,commall=.TRUE.)
   CALL combine_mod(bed_update_ero_glb,bed_update_ero,(/1,1,1/),&
                  & iarr_bed_update_ero,0.0,commall=.TRUE.)
ENDIF

!
!3.2 Update deposition/erosion bed updates
!-----------------------------------------
!

IF (morphstep) THEN
   j_320: DO j=1,nr-1
      bed_update_dep_int_glb(1,j,:) = bed_update_dep_int_glb(1,j,:) + &
                                    & bed_update_dep_glb(1,j,:)
      bed_update_ero_int_glb(1,j,:) = bed_update_ero_int_glb(1,j,:) + &
                                    & bed_update_ero_glb(1,j,:)
   ENDDO j_320
ENDIF


!
!3.2 Write data
!--------------
!

IF (nt.GT.0.AND.mult_index(nt,icount).AND.master) THEN

   rhour = nosecsrun/3600.0
   WRITE (iunit,9001) rhour
   WRITE (iunit,*)

   l_320: DO l=1,nostats
      j = jstat(l)
      WRITE (iunit,9003) 'Station', l
      WRITE (cl,'(I12)') l; cl = ADJUSTL(cl)
      WRITE (iunit,9002) 'depmeanatc_'//TRIM(cl), depmeanatc_glb(1,j)
      WRITE (iunit,9002) 'bstresatc_'//TRIM(cl), bstresatc_glb(1,j)
      WRITE (iunit,9002) 'qtotatv_'//TRIM(cl), qtotatv_glb(1,j,1)
      WRITE (iunit,9002) 'bed_update_int_'//TRIM(cl), bed_update_int_glb(1,j)
      WRITE (iunit,9002) 'bed_update_dep_int_'//TRIM(cl), &
                        & bed_update_dep_int_glb(1,j,1)
      WRITE (iunit,9002) 'bed_update_ero_int_'//TRIM(cl), &
                        & bed_update_ero_int_glb(1,j,1)
      IF (iopt_morph_fixed_layer.EQ.1) THEN
         WRITE (iunit,9002) 'sed_avail_tot_'//TRIM(cl), sed_avail_tot_glb(1,j)
      ENDIF
      WRITE (iunit,9002) 'vert_sed_balance_'//TRIM(cl), &
              & depmeanglb(1,j) - depmean_ini_glb(1,j) + bed_update_int_glb(1,j)
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
   DEALLOCATE (bed_update_dep_int_glb,bed_update_ero_int_glb,depmean_ini_glb)

ENDIF


RETURN

9001 FORMAT ('Time: ',F4.1,' hours')
9002 FORMAT (T3,A,T27,': ',G12.4E3)
9003 FORMAT (A,': ',I3)

END SUBROUTINE usrdef_output
