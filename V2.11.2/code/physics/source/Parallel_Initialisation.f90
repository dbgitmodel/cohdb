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

!************************************************************************
!
! *Parallel_Initialisation* Initialise parallel setup
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Parallel_Initialisation.f90  V2.11
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - 
!
! Reference -
!
! Subroutines - define_model_comms, define_halo_comms, domain_coords,
!               domain_decomposition, read_partition, regular_partition,
!               set_domain_dims, set_procnums, write_partition
!
!************************************************************************
!

!========================================================================

SUBROUTINE define_halo_comms(mglevel)
!************************************************************************
!
! *define_halo_comms* Define parameters for halo send/receive communications
!
! Author - Patrick Luyten and Pieter Rauwoens
!
! Version - @(COHERENS)Parallel_Initialisation.f90  V2.8.1
!
! Description -
!
! Reference -
!
! Calling program - domain_decomposition
!
! Module  calls - mult_index
!
!************************************************************************
!
USE gridpars
USE iopars
USE multigrid
USE paralpars
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: mglevel

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*mglevel*   INTEGER multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j, nxloc, nyloc
TYPE(ExchComms), DIMENSION(MaxHaloComms) :: hcomms


procname(pglev+1) = 'define_halo_comms'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!

IF (mglevel.EQ.0) THEN
   nxloc = ncloc; nyloc = nrloc
ELSE
   nxloc = mgvars(mglevel)%ncloc; nyloc = mgvars(mglevel)%nrloc
ENDIF

i = icoordloc; j = jcoordloc

!
!2. West (S), East (R)
!---------------------
!

hcomms(1)%sfirst = .NOT.mult_index(i,2)
hcomms(1)%iddest = iddomain(i-1,j)
hcomms(1)%idsrce = iddomain(i+1,j)
hcomms(1)%tag = 1
hcomms(1)%i1snd = 1
hcomms(1)%i2snd = (/1,2/)
hcomms(1)%j1snd = 1
hcomms(1)%j2snd = nyloc
hcomms(1)%i1rcv = nxloc+1
hcomms(1)%i2rcv = (/nxloc+1,nxloc+2/)
hcomms(1)%j1rcv = 1
hcomms(1)%j2rcv = nyloc

!
!3. East (S), West (R)
!---------------------
!

hcomms(2)%sfirst = mult_index(i,2)
hcomms(2)%iddest = iddomain(i+1,j)
hcomms(2)%idsrce = iddomain(i-1,j)
hcomms(2)%tag = 2
hcomms(2)%i1snd = (/nxloc,nxloc-1/)
hcomms(2)%i2snd = nxloc
hcomms(2)%j1snd = 1
hcomms(2)%j2snd = nyloc
hcomms(2)%i1rcv = (/0,-1/)
hcomms(2)%i2rcv = 0
hcomms(2)%j1rcv = 1
hcomms(2)%j2rcv = nyloc

!
!4. South (S), North (R)
!-----------------------
!

hcomms(3)%sfirst = .NOT.mult_index(j,2)
hcomms(3)%iddest = iddomain(i,j-1)
hcomms(3)%idsrce = iddomain(i,j+1)
hcomms(3)%tag = 3
hcomms(3)%i1snd = 1
hcomms(3)%i2snd = nxloc
hcomms(3)%j1snd = 1
hcomms(3)%j2snd = (/1,2/)
hcomms(3)%i1rcv = 1
hcomms(3)%i2rcv = nxloc
hcomms(3)%j1rcv = nyloc+1
hcomms(3)%j2rcv = (/nyloc+1,nyloc+2/)

!
!5. North (S), South (R)
!-----------------------
!

hcomms(4)%sfirst = mult_index(j,2)
hcomms(4)%iddest = iddomain(i,j+1)
hcomms(4)%idsrce = iddomain(i,j-1)
hcomms(4)%tag = 4
hcomms(4)%i1snd = 1
hcomms(4)%i2snd = nxloc
hcomms(4)%j1snd = (/nyloc,nyloc-1/)
hcomms(4)%j2snd = nyloc
hcomms(4)%i1rcv = 1
hcomms(4)%i2rcv = nxloc
hcomms(4)%j1rcv = (/0,-1/)
hcomms(4)%j2rcv = 0

!
!6. Southwest (S), Northeast (R)
!-------------------------------
!

hcomms(5)%sfirst = .NOT.mult_index(i,2)
hcomms(5)%iddest = iddomain(i-1,j-1)
hcomms(5)%idsrce = iddomain(i+1,j+1)
hcomms(5)%tag = 5
hcomms(5)%i1snd = 1
hcomms(5)%i2snd = (/1,2/)
hcomms(5)%j1snd = 1
hcomms(5)%j2snd = (/1,2/)
hcomms(5)%i1rcv = nxloc+1
hcomms(5)%i2rcv = (/nxloc+1,nxloc+2/)
hcomms(5)%j1rcv = nyloc+1
hcomms(5)%j2rcv = (/nyloc+1,nyloc+2/)

!
!7. Northeast (S), Southwest (R)
!-------------------------------
!

hcomms(6)%sfirst = mult_index(i,2)
hcomms(6)%iddest = iddomain(i+1,j+1)
hcomms(6)%idsrce = iddomain(i-1,j-1)
hcomms(6)%tag = 6
hcomms(6)%i1snd = (/nxloc,nxloc-1/)
hcomms(6)%i2snd = nxloc
hcomms(6)%j1snd = (/nyloc,nyloc-1/)
hcomms(6)%j2snd = nyloc
hcomms(6)%i1rcv = (/0,-1/)
hcomms(6)%i2rcv = 0
hcomms(6)%j1rcv = (/0,-1/)
hcomms(6)%j2rcv = 0

!
!8. Northwest (S), Southeast (R)
!-------------------------------
!

hcomms(7)%sfirst = .NOT.mult_index(j,2)
hcomms(7)%iddest = iddomain(i-1,j+1)
hcomms(7)%idsrce = iddomain(i+1,j-1)
hcomms(7)%tag = 7
hcomms(7)%i1snd = 1
hcomms(7)%i2snd = (/1,2/)
hcomms(7)%j1snd = (/nyloc,nyloc-1/)
hcomms(7)%j2snd = nyloc
hcomms(7)%i1rcv = nxloc+1
hcomms(7)%i2rcv = (/nxloc+1,nxloc+2/)
hcomms(7)%j1rcv = (/0,-1/)
hcomms(7)%j2rcv = 0

!
!9. Southeast (S), Northwest (R)
!-------------------------------
!

hcomms(8)%sfirst = mult_index(j,2)
hcomms(8)%iddest = iddomain(i+1,j-1)
hcomms(8)%idsrce = iddomain(i-1,j+1)
hcomms(8)%tag = 8
hcomms(8)%i1snd = (/nxloc,nxloc-1/)
hcomms(8)%i2snd = nxloc
hcomms(8)%j1snd = 1
hcomms(8)%j2snd = (/1,2/)
hcomms(8)%i1rcv = (/0,-1/)
hcomms(8)%i2rcv = 0
hcomms(8)%j1rcv = nyloc+1
hcomms(8)%j2rcv = (/nyloc+1,nyloc+2/)

!
!10 .Store
!---------
!

IF (mglevel.EQ.0) THEN
   halocomms = hcomms
ELSE
   mgvars(mglevel)%halocomms = hcomms
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE define_halo_comms

!========================================================================

SUBROUTINE define_model_comms
!************************************************************************
!
! *define_model_comms* define communicators for each model component
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - simulation_start
!
! External calls - 
!
! Module calls - comms_comm_rank, comms_comm_split, error_abort
!
!************************************************************************
!
USE iopars
USE paralpars
USE switches
USE comms_MPI, ONLY : comms_comm_rank, comms_comm_split
USE error_routines, ONLY: error_alloc

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iproc, key


pglev = pglev + 1
procname(pglev) = 'define_model_comms'
WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//procname(pglev)


!---global model ids
modelidcoh = 1
modelidwav = MERGE(modelidcoh+1,0,iopt_waves_model.GT.0)
modelidpart = MERGE(modelidcoh+modelidwav+1,0,iopt_part_model.EQ.2)

!---model ids and keys for communicators
IF (idlocglb.LT.nprocscoh) THEN
   modelid = modelidcoh; nprocs = nprocscoh
   key = idlocglb
   IF (loglev1.GT.0) THEN
      WRITE(iolog,'(A)'), REPEAT(' ',pglev-1)//&
                   & 'This process is assigned to the COHERENS model component.'
   ENDIF
ELSEIF (iopt_waves_model.GT.0.AND.idlocglb.GE.nprocscoh.AND.&
      & idlocglb.LT.npworld-nprocspart) THEN
   modelid = modelidwav; nprocs = nprocswav
   key = idlocglb - nprocscoh
   IF (loglev1.GT.0) THEN
      WRITE(iolog,'(A)'), REPEAT(' ',pglev-1)//&
                   & 'This process is assigned to the wave model component.'
   ENDIF   
ELSEIF (iopt_part_model.EQ.2.AND.idlocglb.EQ.npworld-1) THEN
   modelid = modelidpart; nprocs = 1
   key = npworld - 1
   IF (loglev1.GT.0) THEN
      WRITE(iolog,'(A)'), REPEAT(' ',pglev-1)//&
                   & 'This process is assigned to the particle model component.'
   ENDIF
ENDIF

!---define communicator for each model component
IF (nprocswav.EQ.0.AND.nprocspart.EQ.0) THEN
   comm_model  = comm_world_MPI
ELSE
   CALL comms_comm_split(comm_world_MPI,modelid-1,key,comm_model)
ENDIF

!---local process id
CALL comms_comm_rank(comm_model,idloc)
iprocloc = idloc + 1

!---master process
idmastermod = 0; mastermod = idloc.EQ.idmastermod
idmastercoh = idmaster
IF (nprocspart.EQ.1) idmasterpart = nprocscoh
IF (nprocswav.GT.0) idmasterwav = nprocscoh + nprocspart

!---process ids   
ALLOCATE(idprocs(nprocs),STAT=errstat)
IF (mastermod) CALL error_alloc('idprocs',1,(/nprocs/),kndint)
idprocs = (/(iproc,iproc=0,nprocs-1)/)

pglev = pglev - 1


RETURN

END SUBROUTINE define_model_comms

!========================================================================

SUBROUTINE domain_coords
!************************************************************************
!
! *domain_coords*  Define domain coordinates and ranks
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Parallel_Initialisation.f90  V2.0
!
! Description - 
!
! Reference -
!
! Calling program - domain_decomposition
!
! External calls -
!
! Module calls - error_alloc
!
!************************************************************************
!
USE grid
USE iopars
USE paralpars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=12) :: cprocs
INTEGER :: i, iproc, j, ncmin, nrmin
LOGICAL, DIMENSION(nprocs) :: defflag


procname(pglev+1) = 'domain_coords'
CALL log_timer_in()

!
!1. Domain coordinates
!---------------------
!
!1.1 X-direction
!--------------
!

defflag = .FALSE.
i = 0
DO WHILE (.NOT.ALL(defflag))
   i = i+1
   ncmin = MINVAL(nc1procs,MASK=.NOT.defflag)
   iproc_110: DO iproc=1,nprocs
      IF (nc1procs(iproc).EQ.ncmin) THEN
         defflag(iproc) = .TRUE.
         icoordprocs(iproc) = i
         IF (idloc.EQ.idprocs(iproc)) icoordloc = i
      ENDIF
   ENDDO iproc_110
ENDDO
nprocsx = i

!
!1.2 Y-direction
!---------------
!

defflag = .FALSE.
j = 0
DO WHILE (.NOT.ALL(defflag))
   j = j+1
   nrmin = MINVAL(nr1procs,MASK=.NOT.defflag)
   iproc_120: DO iproc=1,nprocs
      IF (nr1procs(iproc).EQ.nrmin) THEN
         defflag(iproc) = .TRUE.
         jcoordprocs(iproc) = j
         IF (idloc.EQ.idprocs(iproc)) jcoordloc = j
      ENDIF
   ENDDO iproc_120
ENDDO
nprocsy = j

!
!2. Domain ranks
!---------------
!
!---allocate
ALLOCATE (iddomain(0:nprocsx+1,0:nprocsy+1),STAT=errstat)
CALL error_alloc('iddomain',2,(/nprocsx+2,nprocsy+2/),kndint)
iddomain = proc_null_MPI

!---define
i_210: DO iproc=1,nprocs
   i = icoordprocs(iproc); j = jcoordprocs(iproc)
   iddomain(i,j) = idprocs(iproc)
ENDDO i_210

!
!3. Write to log file
!--------------------
!

IF (loglev1.GT.0) THEN
   WRITE (cprocs,'(I12)') npworld; cprocs = ADJUSTL(cprocs)
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'Number of available processes: '//&
                     & TRIM(cprocs)
   WRITE (cprocs,'(I12)') nprocs; cprocs = ADJUSTL(cprocs)
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'Number of active processes: '//&
                     & TRIM(cprocs)
   WRITE (cprocs,'(I12)') nprocsx; cprocs = ADJUSTL(cprocs)
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'Number of processes in X-direction: '//&
                     & TRIM(cprocs)
   WRITE (cprocs,'(I12)') nprocsy; cprocs = ADJUSTL(cprocs)
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'Number of processes in Y-direction: '//&
                     & TRIM(cprocs)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE domain_coords

!========================================================================

SUBROUTINE domain_decomposition
!************************************************************************
!
! *domain_decomposition* Decompose domain for parallel processing
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Parallel_Initialisation.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - define_halo_comms, domain_coords, read_partition,
!                  regular_partition, usrdef_partition, write_partition
!
! Module calls - comms_send_int, error_abort, error_alloc,
!                error_alloc_struc_comp, error_limits_arr
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE modids
USE multigrid
USE paralpars
USE physpars
USE switches
USE syspars
USE comms_MPI, ONLY: comms_send_int
USE error_routines, ONLY: error_abort, error_alloc, error_alloc_struc_comp, &
                        & error_limits_arr
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iproc, lev


procname(pglev+1) = 'domain_decomposition'
CALL log_timer_in()

!
!.1 Allocate arrays
!------------------
!

ALLOCATE(nc1procs(nprocs),STAT=errstat)
CALL error_alloc('nc1procs',1,(/nprocs/),kndint)

ALLOCATE(nc2procs(nprocs),STAT=errstat)
CALL error_alloc('nc2procs',1,(/nprocs/),kndint)

ALLOCATE(ncprocs(nprocs),STAT=errstat)
CALL error_alloc('ncprocs',1,(/nprocs/),kndint)

ALLOCATE(nr1procs(nprocs),STAT=errstat)
CALL error_alloc('nr1procs',1,(/nprocs/),kndint)

ALLOCATE(nr2procs(nprocs),STAT=errstat)
CALL error_alloc('nr2procs',1,(/nprocs/),kndint)

ALLOCATE(nrprocs(nprocs),STAT=errstat)
CALL error_alloc('nrprocs',1,(/nprocs/),kndint)

ALLOCATE(icoordprocs(nprocs),STAT=errstat)
CALL error_alloc('icoordprocs',1,(/nprocs/),kndint)

ALLOCATE(jcoordprocs(nprocs),STAT=errstat)
CALL error_alloc('jcoordprocs',1,(/nprocs/),kndint)

ALLOCATE (nobuprocs(nprocs),STAT=errstat)
CALL error_alloc('nobuprocs',1,(/nprocs/),kndint)

ALLOCATE (nosbuprocs(nprocs),STAT=errstat)
CALL error_alloc('nosbuprocs',1,(/nprocs/),kndint)

ALLOCATE (nrvbuprocs(nprocs),STAT=errstat)
CALL error_alloc('nrvbuprocs',1,(/nprocs/),kndint)

ALLOCATE (nobvprocs(nprocs),STAT=errstat)
CALL error_alloc('nobvprocs',1,(/nprocs/),kndint)

ALLOCATE (nosbvprocs(nprocs),STAT=errstat)
CALL error_alloc('nosbvprocs',1,(/nprocs/),kndint)

ALLOCATE (nrvbvprocs(nprocs),STAT=errstat)
CALL error_alloc('nrvbvprocs',1,(/nprocs/),kndint)

ALLOCATE (nobxprocs(nprocs),STAT=errstat)
CALL error_alloc('nobxprocs',1,(/nprocs/),kndint)

ALLOCATE (nobyprocs(nprocs),STAT=errstat)
CALL error_alloc('nobyprocs',1,(/nprocs/),kndint)

ALLOCATE(comprocs(2*nprocs),STAT=errstat)
CALL error_alloc('comprocs',1,(/2*nprocs/),kndint)

!---multigrid
IF (iopt_hydro_impl.EQ.1) THEN

   lev_110: DO lev=0,nomglevels-1

      ALLOCATE(mgvars(lev)%nc1procs(nprocs),STAT=errstat)
      CALL error_alloc_struc_comp('mgvars','nc1procs',1,(/nprocs/),kndint,lev)

      ALLOCATE(mgvars(lev)%nc2procs(nprocs),STAT=errstat)
      CALL error_alloc_struc_comp('mgvars','nc2procs',1,(/nprocs/),kndint,lev)

      ALLOCATE(mgvars(lev)%nr1procs(nprocs),STAT=errstat)
      CALL error_alloc_struc_comp('mgvars','nr1procs',1,(/nprocs/),kndint,lev)

      ALLOCATE(mgvars(lev)%nr2procs(nprocs),STAT=errstat)
      CALL error_alloc_struc_comp('mgvars','nr2procs',1,(/nprocs/),kndint,lev)

   ENDDO lev_110
   
ENDIF

!
!2. Domain partition
!-------------------
!
!2.1 Serial case
!---------------
!

IF (iopt_MPI.EQ.0) THEN

   nc1procs = 1; nc2procs = nc
   nr1procs = 1; nr2procs = nr

!
!2.2 Regular partition
!---------------------
!

ELSEIF (iopt_MPI_partit.EQ.1) THEN

   CALL regular_partition

!
!2.3 Externally defined
!----------------------
!

ELSEIF (iopt_MPI_partit.EQ.2) THEN

!  ---read
   IF (modfiles(io_mppmod,1,1)%status.EQ.'R') THEN
      CALL read_partition
!  ---user-defined
   ELSEIF (modfiles(io_mppmod,1,1)%status.EQ.'N') THEN
      CALL usrdef_partition
   ENDIF

ENDIF

!
!3. Dimensions of local arrays
!-----------------------------
!

ncprocs = nc2procs - nc1procs + 1
nrprocs = nr2procs - nr1procs + 1

!
!4. Check values
!---------------
!

IF (iopt_MPI.EQ.1.AND.master) THEN
   iproc_41O: DO iproc=1,nprocs
      CALL error_limits_arr(nc1procs(iproc),'nc1procs',1,nc,1,(/iproc/))
      CALL error_limits_arr(nc2procs(iproc),'nc2procs',1,nc,1,(/iproc/))
      CALL error_limits_arr(nr1procs(iproc),'nr1procs',1,nr,1,(/iproc/))
      CALL error_limits_arr(nr2procs(iproc),'nr2procs',1,nr,1,(/iproc/))
      CALL error_limits_arr(ncprocs(iproc),'ncprocs',2,nc,1,(/iproc/))
      CALL error_limits_arr(nrprocs(iproc),'nrprocs',2,nr,1,(/iproc/))
   ENDDO iproc_41O
   CALL error_abort('domain_decomposition',ierrno_inival)
ENDIF

!
!5. Values on local grid
!-----------------------
!

iproc_510: DO iproc=1,nprocs
   IF (idloc.EQ.idprocs(iproc)) THEN
      nc1loc = nc1procs(iproc); nc2loc = nc2procs(iproc); ncloc = ncprocs(iproc)
      nr1loc = nr1procs(iproc); nr2loc = nr2procs(iproc); nrloc = nrprocs(iproc)
   ENDIF
ENDDO iproc_510

!
!6. Write data
!-------------
!

IF (iopt_MPI.EQ.1.AND.modfiles(io_mppmod,1,2)%status.EQ.'W') THEN
   CALL write_partition
ENDIF

!
!7. Domain coordinates and ranks
!-------------------------------
!

CALL domain_coords

!
!8. Define parameters for communications
!---------------------------------------
!
!---combine and collect
comprocs = (/0,(1,iproc=1,nprocs-1),-1,(2,iproc=1,nprocs-1)/)
comprocs = CSHIFT(comprocs,SHIFT=-idloc)

!---exchange
IF (iopt_MPI.EQ.1) CALL define_halo_comms(0)

!
!9. Send domain decomposition to particle module
!-----------------------------------------------
!

IF (master.AND.iopt_part_model.EQ.2) THEN
   CALL comms_send_int(nc1procs,nprocscoh,idmasterpart,1,comm_world_MPI,1,&
                     & iarr_mg_nc1procs)
   CALL comms_send_int(nc2procs,nprocscoh,idmasterpart,2,comm_world_MPI,1,&
                     & iarr_mg_nc2procs)
   CALL comms_send_int(nr1procs,nprocscoh,idmasterpart,3,comm_world_MPI,1,&
                     & iarr_mg_nr1procs)
   CALL comms_send_int(nr2procs,nprocscoh,idmasterpart,4,comm_world_MPI,1,&
                     & iarr_mg_nr2procs)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE domain_decomposition

!========================================================================

SUBROUTINE read_partition
!************************************************************************
!
! *read_partition* Read domain partition in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Parallel_Initialisation.f90  V2.1.0
!
! Description -
!
! Reference -
!
! Calling program - domain_decomposition
!
! External calls - close_filepars, error_alloc, error_alloc_struc,
!                  open_filepars, read_glbatts_mod, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE grid
USE iopars
USE multigrid
USE paralpars
USE physpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!* Local variables
!
INTEGER :: lev, ilev, numvars
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nx1procs, nx2procs, ny1procs, ny2procs
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_partition'
CALL log_timer_in()

filepars = modfiles(io_mppmod,1,1)

!
!1. File header
!--------------
!

CALL open_filepars(filepars)
CALL read_glbatts_mod(filepars)
numvars = filepars%novars
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL read_varatts_mod(filepars,varatts,numvars)

!
!2. Read partition
!-----------------
!
!---allocate
ALLOCATE (nx1procs(nprocs,nomglevels),STAT=errstat)
CALL error_alloc('nx1procs',2,(/nprocs,nomglevels/),kndint)
ALLOCATE (nx2procs(nprocs,nomglevels),STAT=errstat)
CALL error_alloc('nx2procs',2,(/nprocs,nomglevels/),kndint)
ALLOCATE (ny1procs(nprocs,nomglevels),STAT=errstat)
CALL error_alloc('ny1procs',2,(/nprocs,nomglevels/),kndint)
ALLOCATE (ny2procs(nprocs,nomglevels),STAT=errstat)
CALL error_alloc('ny2procs',2,(/nprocs,nomglevels/),kndint)

!---read
CALL read_vars(nx1procs,filepars,1,varatts=varatts)
CALL read_vars(nx2procs,filepars,2,varatts=varatts)
CALL read_vars(ny1procs,filepars,3,varatts=varatts)
CALL read_vars(ny2procs,filepars,4,varatts=varatts)

!---store
nc1procs = nx1procs(:,1); nc2procs = nx2procs(:,1)
nr1procs = ny1procs(:,1); nr2procs = ny2procs(:,1)
IF (iopt_hydro_impl.EQ.1) THEN
   lev_210: DO lev=0,nomglevels-1
      ilev = lev + 1
      mgvars(lev)%nc1procs = nx1procs(:,ilev)
      mgvars(lev)%nc2procs = nx2procs(:,ilev)
      mgvars(lev)%nr1procs = ny1procs(:,ilev)
      mgvars(lev)%nr2procs = ny2procs(:,ilev)
   ENDDO lev_210
ENDIF

!
!3. Finalise
!-----------
!
 
CALL close_filepars(filepars)
modfiles(io_mppmod,1,1) = filepars

!---deallocate
DEALLOCATE (varatts)
DEALLOCATE (nx1procs,nx2procs,ny1procs,ny2procs)

CALL log_timer_out()


RETURN

END SUBROUTINE read_partition

!========================================================================

SUBROUTINE regular_partition
!************************************************************************
!
! *regular_partition* Decomposes the computational domain using simple
!                     partitioning into nprocsx*nprocsy subdomains
!
! Author - Patrick Luyten and Pieter Rauwoens
!
! Version - @(COHERENS)Parallel_Initialisation.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - domain_decomposition
!
! External calls -
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE paralpars
USE physpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, iproc, j, group, ncmod, ncg, nrg, nrmod, nxdiv, nxmod, nydiv, &
         & nymod


procname(pglev+1) = 'regular_partition'
CALL log_timer_in()

group = 2**(nomglevels-1)
ncmod = MOD(nc,group); nrmod = MOD(nr,group)
IF (ncmod.EQ.0) THEN
   ncg = nc/group; ncmod = group
ELSE
   ncg = nc/group + 1
ENDIF
IF (nrmod.EQ.0) THEN
   nrg = nr/group; nrmod = group
ELSE
   nrg = nr/group + 1
ENDIF

nxdiv = ncg/nprocsx*group; nxmod = MOD(ncg,nprocsx)
nydiv = nrg/nprocsy*group; nymod = MOD(nrg,nprocsy)
j_110: DO j=1,nprocsy 
i_110: DO i=1,nprocsx
   iproc = (j-1)*nprocsx + i
   IF (i.LE.(nprocsx-nxmod)) THEN
      nc1procs(iproc) = 1 + (i-1)*nxdiv
      ncloc = nxdiv
   ELSE
      nc1procs(iproc) = 1 + (nprocsx-nxmod)*nxdiv + &
                          & (i-(nprocsx-nxmod+1))*(nxdiv+group)
      ncloc = nxdiv + group
   ENDIF
   IF (i.EQ.nprocsx) ncloc = ncloc - group + ncmod
   nc2procs(iproc) = nc1procs(iproc) + ncloc - 1 
   IF (j.LE.(nprocsy-nymod)) THEN
      nr1procs(iproc) = 1 + (j-1)*nydiv
      nrloc = nydiv
   ELSE
      nr1procs(iproc) = 1 + (nprocsy-nymod)*nydiv + &
                          & (j-(nprocsy-nymod+1))*(nydiv+group)
      nrloc = nydiv + group
   ENDIF
   IF( j.EQ.nprocsy) nrloc = nrloc - group + nrmod
   nr2procs(iproc) = nr1procs(iproc) + nrloc - 1
ENDDO i_110
ENDDO j_110

CALL log_timer_out()


RETURN

END SUBROUTINE regular_partition

!========================================================================

SUBROUTINE set_domain_dims(np,npx,npy)
!************************************************************************
!
! *set_domain_dims* Set domain dimensions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Parallel_Initialisation.f90  V2.0
!
! Description - 
!
! Reference -
!
! Calling program - set_procnums
!
! External calls -
!
! Module calls - comms_dims_create, error_abort
!
!************************************************************************
!
USE gridpars
USE iopars
USE comms_MPI, ONLY: comms_dims_create
USE error_routines, ONLY: error_abort
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
INTEGER, INTENT(IN) :: np
INTEGER, INTENT(INOUT) :: npx, npy

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*np*      INTEGER Total number of processes
!*npx*     INTEGER X-dimension of Cartesian grid
!*npy*     INTEGER Y-dimension of Cartesian grid
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(2) :: dims


procname(pglev+1) = 'set_domain_dims'
CALL log_timer_in()

dims = 0
CALL comms_dims_create(np,2,dims)

IF (nc.GE.nr) THEN
   npx = dims(1); npy = dims(2)
ELSE
   npx = dims(2); npy = dims(1)
ENDIF

CALL error_abort(procname(pglev),ierrno_MPI)

CALL log_timer_out()


RETURN

END SUBROUTINE set_domain_dims

!========================================================================

SUBROUTINE set_procnums
!************************************************************************
!
! *set_procnums* Number of processes used in the simulation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Parallel_Initialisation.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - set_domain_dims
!
! Module calls - error_abort, error_limits_var, error_mult, error_value_var,
!                mult_index, warning_reset_var
!
!************************************************************************
!
USE iopars
USE paralpars
USE switches
USE error_routines, ONLY: error_abort, error_limits_var, error_mult, &
                        & error_value_var, warning_reset_var
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY : mult_index

IMPLICIT NONE


procname(pglev+1) = 'set_procnums'
CALL log_timer_in()

!
!1. Serial case
!--------------
!

IF (.NOT.parallel_set.OR.iopt_grid_nodim.EQ.1) THEN

   IF (nprocsx.NE.1) CALL warning_reset_var(nprocsx,'nprocsx',1)
   IF (nprocsy.NE.1) CALL warning_reset_var(nprocsy,'nprocsy',1)

!
!2. Reset values for regular partition
!-------------------------------------
!

ELSEIF (iopt_MPI_partit.EQ.1) THEN

!
!2.1 Check values
!----------------
!

   CALL error_limits_var(nprocsx,'nprocsx',0,nprocscoh)
   CALL error_limits_var(nprocsy,'nprocsy',0,nprocscoh)

!
!2.2 Reset values (if necessary)
!-------------------------------
!

   IF ((nprocsx.GT.0).AND.(nprocsy.GT.0)) THEN
      nprocs = nprocsx*nprocsy
      CALL error_value_var(nprocs,'nprocs',nprocscoh)
   ENDIF
   IF (nprocs.EQ.1) THEN
      nprocsx = 1; nprocsy = 1
   ELSEIF ((nprocsx.EQ.0).AND.(nprocsy.EQ.0)) THEN
      CALL set_domain_dims(nprocs,nprocsx,nprocsy)
   ELSEIF ((nprocsx.GT.0).AND.(nprocsy.EQ.0)) THEN
      IF (mult_index(nprocs,nprocsx)) THEN
         nprocsy = nprocs/nprocsx
      ELSE
         CALL error_mult(nprocsx,'nprocsx',nprocs)
      ENDIF
   ELSEIF ((nprocsx.EQ.0).AND.(nprocsy.GT.0)) THEN
      IF (mult_index(nprocs,nprocsy)) THEN
         nprocsx = nprocs/nprocsy
      ELSE
         CALL error_mult(nprocsy,'nprocsy',nprocs)
      ENDIF
   ENDIF

ENDIF

CALL error_abort('set_procnums',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE set_procnums

!========================================================================

SUBROUTINE write_partition
!************************************************************************
!
! *write_partition* Write domain partition in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Parallel_Initialisation.f90  V2.10
!
! Description -
!
! Reference -
!
! Calling program - domain_decomposition
!
! External calls - close_filepars, error_alloc_struc, open_filepars,
!                  set_modfiles_atts, set_modvars_atts, write_atts_mod,
!                  write_vars
!
!************************************************************************
!
USE datatypes
USE grid
USE iopars
USE multigrid
USE paralpars
USE physpars
USE syspars
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!* Local variables
!
INTEGER :: ilev, lev, numvars
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nx1procs, nx2procs, ny1procs, ny2procs
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN 

procname(pglev+1) = 'write_partition'
CALL log_timer_in()

!
!1. Write file header
!--------------------
!
!---file attributes
CALL set_modfiles_atts(io_mppmod,1,2)
filepars = modfiles(io_mppmod,1,2)
numvars = filepars%novars

!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL set_modvars_atts(io_mppmod,1,2,varatts,numvars)

!---write
CALL write_atts_mod(filepars,varatts,numvars)

!
!2. Write partition
!------------------
!
!---allocate
ALLOCATE (nx1procs(nprocs,nomglevels),STAT=errstat)
CALL error_alloc('nx1procs',2,(/nprocs,nomglevels/),kndint)
ALLOCATE (nx2procs(nprocs,nomglevels),STAT=errstat)
CALL error_alloc('nx2procs',2,(/nprocs,nomglevels/),kndint)
ALLOCATE (ny1procs(nprocs,nomglevels),STAT=errstat)
CALL error_alloc('ny1procs',2,(/nprocs,nomglevels/),kndint)
ALLOCATE (ny2procs(nprocs,nomglevels),STAT=errstat)
CALL error_alloc('ny2procs',2,(/nprocs,nomglevels/),kndint)

!---store
nx1procs(:,1) = nc1procs; nx2procs(:,1) = nc2procs
ny1procs(:,1) = nr1procs; ny2procs(:,1) = nr2procs
IF (nomglevels.GT.1 )THEN
   lev_210: DO lev=1,nomglevels-1
      ilev = lev + 1
      nx1procs(:,ilev) = mgvars(lev)%nc1procs
      nx2procs(:,ilev) = mgvars(lev)%nc2procs
      ny1procs(:,ilev) = mgvars(lev)%nr1procs
      ny2procs(:,ilev) = mgvars(lev)%nr2procs
   ENDDO lev_210
ENDIF

!---write
CALL write_vars(nx1procs,filepars,1,varatts=varatts)
CALL write_vars(nx2procs,filepars,2,varatts=varatts)
CALL write_vars(ny1procs,filepars,3,varatts=varatts)
CALL write_vars(ny2procs,filepars,4,varatts=varatts)

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_mppmod,1,2) = filepars

!---deallocate
DEALLOCATE (varatts)
DEALLOCATE (nx1procs,nx2procs,ny1procs,ny2procs)

CALL log_timer_out()


RETURN

END SUBROUTINE write_partition
