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
! Author - ANTEA/MUMM
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.4
!
! $Date: 2018-07-23 10:21:02 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1167 $
!
! Description - output parameters for test case flood2d
!             - flooding/drying channel experiments
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_obc, combine_stats_glb, error_alloc,
!                max_vars, min_vars, mult_index, open_file, sum2_vars
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: combine_stats_glb
USE paral_utilities, ONLY: max_vars, min_vars, sum2_vars
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: outflag
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask1, mask2, mask3
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER, PARAMETER :: novars = 11
INTEGER, SAVE :: icout, iunit
INTEGER :: i, iglb, ii, iiloc, j, nc1, nc2
INTEGER, SAVE, DIMENSION(4) :: nhdims
REAL, SAVE :: dhun, volerr
REAL :: ddepint, depmaxerr, depmeanerr, deptotint, deptotint1, &
      & deptotint2, deptotint3, dryarea, hour, transobu, uvelmax, uvelmin, &
      & wphysmax, wphysmin
REAL, DIMENSION(nobu) :: disobu
REAL, DIMENSION(novars) :: out0ddat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: dep_old
REAL, DIMENSION(ncloc,nrloc) :: ddeptot


procname(pglev+1) = 'usrdef_output'
CALL log_timer_in()

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
!  ---initial conditions
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = modfiles(io_fincon,ics_phys,2)%form
   modfiles(io_inicon,ics_phys,1)%filename = &
 & modfiles(io_fincon,ics_phys,2)%filename
   modfiles(io_fincon,ics_phys,2)%status = '0'
!  ---open boundary conditions (2-D)
   modfiles(io_2uvobc,1,1)%status = 'R'
   modfiles(io_2uvobc,1,1)%form = modfiles(io_2uvobc,1,2)%form
   modfiles(io_2uvobc,1,1)%filename = modfiles(io_2uvobc,1,2)%filename
   modfiles(io_2uvobc,1,2)%status = '0'
   GOTO 1001
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN

!
!2.1 Open files, write headers
!-----------------------------
!
!  ---.tst file
   IF (master) THEN
      CALL open_file(iunit,TRIM(outtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case flood2d:'//&
                         ' simulation '//TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF

!
!2.2 Allocate
!------------
!

   ALLOCATE (mask1(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('mask1',2,(/ncloc,nrloc/),log_type)
   ALLOCATE (mask2(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('mask2',2,(/ncloc,nrloc/),log_type)
   ALLOCATE (mask3(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('mask3',2,(/ncloc,nrloc/),log_type)
   ALLOCATE (dep_old(ncloc,nrloc),STAT=errstat)
   CALL error_alloc ('dep_old',2,(/ncloc,nrloc/),real_type)

!
!2.3 Initialise parameters and arrays
!------------------------------------
!

   SELECT CASE (runtitle(8:8))
      CASE ('A','B')
         icout = 600
         nc1 = 19; nc2 = 60
      CASE ('C','D')
         nc1 = 45; nc2 = 55
         icout= 600
   END SELECT
   nhdims = 0
   dhun = surfacegrids(igrd_model,1)%delydat
   volerr = 0.0
   dep_old = deptotatc(1:ncloc,1:nrloc)
   mask1(:,2) = .FALSE.; mask2(:,2) = .FALSE.; mask3(:,2) = .FALSE.
   i_240: DO i=1,ncloc
      iglb = i + nc1loc - 1
      mask1(i,1) = MERGE(.TRUE.,.FALSE.,iglb.LE.nc1) 
      mask2(i,1) = MERGE(.TRUE.,.FALSE.,iglb.GT.nc1.AND.iglb.LE.nc2)
      mask3(i,1) = MERGE(.TRUE.,.FALSE.,iglb.GT.nc2.AND.iglb.LT.nc)
   ENDDO i_240

ENDIF

IF (.NOT.corrstep) GOTO 1001

!
!3. Update parameters
!--------------------
!

outflag = mult_index(nt,icout).AND.nt.GT.0

!
!3.1 Volume error
!----------------
!
!---open boundary transport term
iiloc_310: DO iiloc=1,nobuloc
   i = iobuloc(iiloc); j = jobuloc(iiloc)
   ii = indexobu(iiloc)
   disobu(ii) = 0.5*(udvel_old(i,j)+udvel(i,j))
ENDDO iiloc_310
IF (iopt_MPI.EQ.1) THEN
   CALL combine_stats_glb(disobu,nobu,nobuprocs,indexobuprocs,0,commall=.TRUE.)
ENDIF
transobu = SUM(disobu)

!---volume error
IF (nt.GT.0) THEN
   CALL sum2_vars(deptotatc(1:ncloc,1:nrloc),deptotint,nhdims,'C  ',0,&
                & commall=.TRUE.,mask=seapoint(1:ncloc,1:nrloc))
   ddeptot = deptotatc(1:ncloc,1:nrloc) - dep_old
   CALL sum2_vars(ddeptot,ddepint,nhdims,'C  ',0,commall=.TRUE.,&
                & mask=seapoint(1:ncloc,1:nrloc))
   volerr = volerr + ddepint - delt3d*transobu/dhun
   dep_old = deptotatc(1:ncloc,1:nrloc)
ENDIF

!
!3.2 Depth error
!---------------
!

IF (outflag) THEN
   CALL sum2_vars(deptotatc_err,depmeanerr,nhdims,'C  ',0,commall=.TRUE.,&
                & mask=seapoint(1:ncloc,1:nrloc))
   depmeanerr = depmeanerr/REAL(noseaatc)
   CALL max_vars(deptotatc_err,depmaxerr,0,commall=.TRUE.,&
               & mask=seapoint(1:ncloc,1:nrloc)) 
ENDIF   

!
!3.3 Dry area
!------------
!

IF (outflag) dryarea = 1.0-nowetatc/REAL(noseaatc)

!
!3.4 Water volume in transects
!-----------------------------
!

IF (outflag) THEN
   CALL sum2_vars(deptotatc(1:ncloc,1:nrloc),deptotint1,nhdims,'C  ',0,&
                & commall=.TRUE.,mask=mask1)
   CALL sum2_vars(deptotatc(1:ncloc,1:nrloc),deptotint2,nhdims,'C  ',0,&
                & commall=.TRUE.,mask=mask2)
   CALL sum2_vars(deptotatc(1:ncloc,1:nrloc),deptotint3,nhdims,'C  ',0,&
                & commall=.TRUE.,mask=mask3)
ENDIF

!
!3.5 Current extrema
!-------------------
!

IF (outflag) THEN
   CALL max_vars(uvel(1:ncloc,1:nrloc,:),uvelmax,iarr_uvel,commall=.TRUE.,&
               & mask=node2du(1:ncloc,1:nrloc).GT.1)
   CALL min_vars(uvel(1:ncloc,1:nrloc,:),uvelmin,iarr_uvel,commall=.TRUE.,&
               & mask=node2du(1:ncloc,1:nrloc).GT.1)
   CALL max_vars(wphys,wphysmax,iarr_wphys,commall=.TRUE.,mask=maskatc_int)
   CALL min_vars(wphys,wphysmin,iarr_wphys,commall=.TRUE.,mask=maskatc_int)
ENDIF

!
!3.6 Store for output
!--------------------
!

IF (outflag) THEN
   out0ddat(1)  = 100.0*volerr/deptotint
!   out0ddat(2)  = 100.0*depmeanerr
!   out0ddat(3)  = 100.0*depmaxerr
   out0ddat(4)  = 100.0*dryarea
   out0ddat(5)  = 100.0*deptotint1/deptotint
   out0ddat(6)  = 100.0*deptotint2/deptotint
   out0ddat(7)  = 100.0*deptotint3/deptotint
   out0ddat(8)  = 100.0*uvelmax
   out0ddat(9)  = 100.0*uvelmin
   out0ddat(10) = 1000.0*wphysmax
   out0ddat(11) = 1000.0*wphysmin
ENDIF

!
!4. Write output
!---------------
!

IF (master.AND.outflag) THEN
   hour = nosecsrun/3600.0
   WRITE (iunit,9001) hour
   WRITE (iunit,9002) 'volerr',    out0ddat(1)
!   WRITE (iunit,9002) 'depmeanerr',out0ddat(2)
!   WRITE (iunit,9002) 'depmaxerr', out0ddat(3)
   WRITE (iunit,9002) 'dryarea',   out0ddat(4)
   WRITE (iunit,9002) 'volume1',   out0ddat(5)
   WRITE (iunit,9002) 'volume2',   out0ddat(6)
   WRITE (iunit,9002) 'volume3',   out0ddat(7)
   WRITE (iunit,9002) 'uvelmax',   out0ddat(8)
   WRITE (iunit,9002) 'uvelmin',   out0ddat(9)
   WRITE (iunit,9002) 'wphysmax',  out0ddat(10)
   WRITE (iunit,9002) 'wphysmin',  out0ddat(11)
ENDIF

!
!5. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')

!
!6. Deallocate
!-------------
!

IF (nt.EQ.nstep) DEALLOCATE (mask1,mask2,mask3,dep_old)

1001 CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',F4.1,' hours')
9002 FORMAT (T3,A,T13,': ',G12.4E3)

END SUBROUTINE usrdef_output
