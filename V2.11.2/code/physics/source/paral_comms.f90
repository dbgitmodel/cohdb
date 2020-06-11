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

MODULE paral_comms
!************************************************************************
!
! *paral_comms* Parallel communication library
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.11
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - collect, combine, copy, distribute and exchange operations
!
! Generic routines - collect_vars, combine_mod, combine_stats_glb,
!                    combine_stats_loc, combine_submod, copy_chars, copy_vars,
!                    distribute_mod, exchange_mod
!
! Routines - collect_log_expr, combine_mod_hgrid_2d, distribute_mod_hgrid_2d,
!            halo_exch_comms
!
!************************************************************************
!
USE gridpars
USE iopars
USE paralpars
USE switches
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

INTERFACE collect_vars
   MODULE PROCEDURE collect_vars_int_0d, collect_vars_int_1d, &
    & collect_vars_int_2d, collect_vars_int_3d, collect_vars_int_4d, &
    & collect_vars_log_0d, collect_vars_log_1d, collect_vars_log_2d, &
    & collect_vars_real_0d, collect_vars_real_1d, collect_vars_real_2d, &
    & collect_vars_real_3d, collect_vars_real_4d
END INTERFACE

INTERFACE combine_mod
   MODULE PROCEDURE combine_mod_int_2d, combine_mod_int_3d, &
  & combine_mod_int_4d, combine_mod_log_2d, combine_mod_real_2d, &
  & combine_mod_real_3d, combine_mod_real_4d
END INTERFACE

INTERFACE combine_stats_glb
   MODULE PROCEDURE combine_stats_glb_int_1d, combine_stats_glb_int_2d, &
    & combine_stats_glb_int_3d, combine_stats_glb_int_4d, &
    & combine_stats_glb_log_1d, combine_stats_glb_real_1d, &
    & combine_stats_glb_real_2d, combine_stats_glb_real_3d, &
    & combine_stats_glb_real_4d
END INTERFACE

INTERFACE combine_stats_loc
   MODULE PROCEDURE combine_stats_loc_real_1d, combine_stats_loc_real_2d, &
                  & combine_stats_loc_real_3d
END INTERFACE

INTERFACE combine_submod
   MODULE PROCEDURE combine_submod_real_2d, combine_submod_real_3d, &
                  & combine_submod_real_4d
END INTERFACE

INTERFACE copy_chars
   MODULE PROCEDURE copy_chars_0d, copy_chars_1d, copy_chars_2d, &
                  & copy_chars_3d, copy_chars_4d
END INTERFACE

INTERFACE copy_vars
   MODULE PROCEDURE copy_vars_int_0d, copy_vars_int_1d, copy_vars_int_2d, &
    & copy_vars_int_3d, copy_vars_int_4d, copy_vars_log_0d, copy_vars_log_1d, &
    & copy_vars_log_2d, copy_vars_log_3d, copy_vars_log_4d, copy_vars_real_0d, &
    & copy_vars_real_1d, copy_vars_real_2d, copy_vars_real_3d, copy_vars_real_4d
END INTERFACE

INTERFACE distribute_mod
   MODULE PROCEDURE distribute_mod_int_2d, distribute_mod_int_3d, &
    & distribute_mod_int_4d, distribute_mod_log_2d, distribute_mod_real_2d, &
    & distribute_mod_real_3d, distribute_mod_real_4d
END INTERFACE

INTERFACE exchange_mod
   MODULE PROCEDURE exchange_mod_int_2d, exchange_mod_int_3d, &
    & exchange_mod_int_4d, exchange_mod_log_2d, exchange_mod_real_2d, &
    & exchange_mod_real_3d, exchange_mod_real_4d
END INTERFACE

CONTAINS

!========================================================================

FUNCTION collect_log_expr(logvar,imode,comm,commtype)
!************************************************************************
!
! *collect_log_expr* Evaluates a logical expression on all processes
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)paral_comms.f90  V2.8
!
! Description - result depends on value imode
!                1 => returns .TRUE. if logvar is.TRUE. on at least one process 
!                2 => returns .TRUE. if logvar is.TRUE. on all processes 
!
! Module calls - collect_vars_log_0d
!
!************************************************************************
!
!* Arguments
!
LOGICAL, INTENT(IN) :: logvar
INTEGER, INTENT(IN) :: imode
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
LOGICAL :: collect_log_expr

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*logvar*    INTEGER Local integer data to be collected
!*imode*     INTEGER Type of result 
!*comm*      INTEGER Communicator for collect operation
!*commtype*  INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, itype
LOGICAL, DIMENSION(nprocs) :: logprocs


procname(pglev+1) = 'collect_log_expr'
CALL log_timer_in()

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   collect_log_expr = logvar
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
!

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!3. Collect logical expression on all processes
!----------------------------------------------
!

CALL collect_vars_log_0d(logvar,logprocs,nprocs,0,commop,itype)

!
!4. Return result
!----------------
!

IF (imode.EQ.1) THEN
   collect_log_expr = ANY(logprocs)
ELSEIF (imode.EQ.2) THEN
   collect_log_expr = ALL(logprocs)
ENDIF
   
1000 CALL log_timer_out()

RETURN

END FUNCTION collect_log_expr

!===========================================================================

SUBROUTINE collect_vars_int_0d(intloc,intprocs,noprocs,ivarid,comm,commtype)
!***************************************************************************
!
! *collect_vars_int_0d* Collect a global integer "process" array from local
!                       scalars
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_allgather_int, comms_barrier, comms_irecv_int,
!                comms_isend_int, comms_recv_int, comms_send_int,
!                comms_waitall, error_abort
!
!***************************************************************************
!
USE comms_MPI, ONLY: comms_allgather_int, comms_barrier, &
                   & comms_irecv_int, comms_isend_int, comms_recv_int, &
                   & comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: intloc, ivarid, noprocs
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
INTEGER, INTENT(INOUT), DIMENSION (:) :: intprocs

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*intloc*    INTEGER Local integer data to be collected
!*intprocs*  INTEGER Global process array
!*noprocs*   INTEGER Number of processes in communicator
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for collect operation
!*commtype*  INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iform, imode, iproc, itype, n, npcc, nreq
INTEGER, DIMENSION(2*noprocs-2) :: irequest
INTEGER, DIMENSION(1) :: int1d


procname(pglev+1) = 'collect_vars_int_0d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   intprocs(1) = intloc
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
!

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
int1d(1) = intloc

!
!4. Blocking
!-----------
!

IF (iform.EQ.1) THEN

   n_410: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         CALL comms_send_int(int1d,1,idprocs(iproc),idloc,commop,imode,ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         CALL comms_recv_int(intprocs(iproc:iproc),1,idprocs(iproc),&
                           & idprocs(iproc),commop,ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         intprocs(iproc) = intloc
      ENDIF
   ENDDO n_410

!
!5. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         nreq = nreq + 1
         CALL comms_isend_int(int1d,1,idprocs(iproc),idloc,commop,imode,&
                            & irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         nreq = nreq + 1
         CALL comms_irecv_int(intprocs(iproc:iproc),1,idprocs(iproc),&
                            & idprocs(iproc),commop,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         intprocs(iproc) = intloc
      ENDIF
   ENDDO n_510

!  ---wait for completion
   CALL comms_waitall(nreq,irequest)

!
!6. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_allgather_int(int1d,1,noprocs,intprocs,commop,ivarid)
ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_coll)


RETURN

END SUBROUTINE collect_vars_int_0d

!========================================================================

SUBROUTINE collect_vars_int_1d(intloc,intprocs,noprocs,ivarid,comm,commtype)
!************************************************************************
!
! *collect_vars_int_1d* Collect a global integer "process" array from local
!                       vector arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_allgather_int, comms_barrier, comms_irecv_int,
!                comms_isend_int, comms_recv_int, comms_send_int,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_allgather_int, comms_barrier, comms_irecv_int, &
                   & comms_isend_int, comms_recv_int, comms_send_int, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, noprocs
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
INTEGER, INTENT(IN), DIMENSION(:) :: intloc
INTEGER, INTENT(INOUT), DIMENSION (:,:) :: intprocs

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*intloc*    INTEGER Local integer data to be collected
!*intprocs*  INTEGER Global process array
!*noprocs*   INTEGER Number of processes in communicator
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for collect operation
!*commtype*  INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iform, imode, iproc, itype, n, ncount, npcc, nreq
INTEGER, DIMENSION(2*noprocs-2) :: irequest


ncount = SIZE(intloc)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'collect_vars_int_1d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   intprocs(:,1) = intloc
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
!

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Blocking
!-----------
!

IF (iform.EQ.1) THEN

   n_410: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         CALL comms_send_int(intloc,ncount,idprocs(iproc),idloc,commop,imode,&
                           & ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         CALL comms_recv_int(intprocs(:,iproc),ncount,idprocs(iproc),&
                           & idprocs(iproc),commop,ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         intprocs(:,iproc) = intloc
      ENDIF
   ENDDO n_410

!
!5. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         nreq = nreq + 1
         CALL comms_isend_int(intloc,ncount,idprocs(iproc),idloc,commop,&
                            & imode,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         nreq = nreq + 1
         CALL comms_irecv_int(intprocs(:,iproc),ncount,idprocs(iproc),&
                            & idprocs(iproc),commop,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         intprocs(:,iproc) = intloc
      ENDIF
   ENDDO n_510

!  ---wait for completion
   CALL comms_waitall(nreq,irequest)

!
!6. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_allgather_int(intloc,ncount,noprocs,intprocs,commop,ivarid)
ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_coll)


RETURN

END SUBROUTINE collect_vars_int_1d

!========================================================================

SUBROUTINE collect_vars_int_2d(intloc,intprocs,noprocs,ivarid,comm,commtype)
!************************************************************************
!
! *collect_vars_int_2d* Collect a global integer "process" array from local
!                       2-D arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_allgather_int, comms_barrier, comms_irecv_int,
!                comms_isend_int, comms_recv_int, comms_send_int,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_allgather_int, comms_barrier, comms_irecv_int, &
                   & comms_isend_int, comms_recv_int, comms_send_int, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, noprocs
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
INTEGER, INTENT(IN), DIMENSION(:,:) :: intloc
INTEGER, INTENT(INOUT), DIMENSION (:,:,:) :: intprocs

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*intloc*    INTEGER Local integer data to be collected
!*intprocs*  INTEGER Global process array
!*noprocs*   INTEGER Number of processes in communicator
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for collect operation
!*commtype*  INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iform, imode, iproc, itype, n, ncount, nreq, npcc
INTEGER, DIMENSION(2*noprocs-2) :: irequest


ncount = SIZE(intloc)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'collect_vars_int_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   intprocs(:,:,1) = intloc
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
!

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Blocking
!-----------
!

IF (iform.EQ.1) THEN

   n_410: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         CALL comms_send_int(intloc,ncount,idprocs(iproc),idloc,commop,imode,&
                           & ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         CALL comms_recv_int(intprocs(:,:,iproc),ncount,idprocs(iproc),&
                           & idprocs(iproc),commop,ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         intprocs(:,:,iproc) = intloc
      ENDIF
   ENDDO n_410

!
!5. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         nreq = nreq + 1
         CALL comms_isend_int(intloc,ncount,idprocs(iproc),idloc,commop,&
                            & imode,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         nreq = nreq + 1
         CALL comms_irecv_int(intprocs(:,:,iproc),ncount,idprocs(iproc),&
                            & idprocs(iproc),commop,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         intprocs(:,:,iproc) = intloc
      ENDIF
   ENDDO n_510

!  ---wait for completion
   CALL comms_waitall(nreq,irequest)

!
!6. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_allgather_int(intloc,ncount,noprocs,intprocs,commop,ivarid)
ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_coll)


RETURN

END SUBROUTINE collect_vars_int_2d

!========================================================================

SUBROUTINE collect_vars_int_3d(intloc,intprocs,noprocs,ivarid,comm,commtype)
!************************************************************************
!
! *collect_vars_int_3d* Collect a global integer "process" array from local
!                       3-D arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_allgather_int, comms_barrier, comms_irecv_int,
!                comms_isend_int, comms_recv_int, comms_send_int,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_allgather_int, comms_barrier, comms_irecv_int, &
                   & comms_isend_int, comms_recv_int, comms_send_int, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, noprocs
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
INTEGER, INTENT(IN), DIMENSION(:,:,:) :: intloc
INTEGER, INTENT(INOUT), DIMENSION (:,:,:,:) :: intprocs

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*intloc*    INTEGER Local integer data to be collected
!*intprocs*  INTEGER Global process array
!*noprocs*   INTEGER Number of processes in communicator
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for collect operation
!*commtype*  INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iform, imode, iproc, itype, n, ncount, nreq, npcc
INTEGER, DIMENSION(2*noprocs-2) :: irequest


ncount = SIZE(intloc)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'collect_vars_int_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   intprocs(:,:,:,1) = intloc
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
!

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Blocking
!-----------
!

IF (iform.EQ.1) THEN

   n_410: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         CALL comms_send_int(intloc,ncount,idprocs(iproc),idloc,commop,imode,&
                           & ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         CALL comms_recv_int(intprocs(:,:,:,iproc),ncount,idprocs(iproc),&
                           & idprocs(iproc),commop,ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         intprocs(:,:,:,iproc) = intloc
      ENDIF
   ENDDO n_410

!
!5. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         nreq = nreq + 1
         CALL comms_isend_int(intloc,ncount,idprocs(iproc),idloc,commop,&
                            & imode,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         nreq = nreq + 1
         CALL comms_irecv_int(intprocs(:,:,:,iproc),ncount,idprocs(iproc),&
                            & idprocs(iproc),commop,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         intprocs(:,:,:,iproc) = intloc
      ENDIF
   ENDDO n_510

!  ---wait for completion
   CALL comms_waitall(nreq,irequest)

!
!6. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_allgather_int(intloc,ncount,noprocs,intprocs,commop,ivarid)
ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_coll)


RETURN

END SUBROUTINE collect_vars_int_3d

!========================================================================

SUBROUTINE collect_vars_int_4d(intloc,intprocs,noprocs,ivarid,comm,commtype)
!************************************************************************
!
! *collect_vars_int_4d* Collect a global integer "process" array from local
!                       4-D arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_allgather_int, comms_barrier, comms_irecv_int,
!                comms_isend_int, comms_recv_int, comms_send_int,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_allgather_int, comms_barrier, comms_irecv_int, &
                   & comms_isend_int, comms_recv_int, comms_send_int, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, noprocs
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
INTEGER, INTENT(IN), DIMENSION(:,:,:,:) :: intloc
INTEGER, INTENT(INOUT), DIMENSION (:,:,:,:,:) :: intprocs

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*intloc*    INTEGER Local integer data to be collected
!*intprocs*  INTEGER Global process array
!*noprocs*   INTEGER Number of processes in communicator
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for collect operation
!*commtype*  INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iform, imode, iproc, itype, n, ncount, nreq, npcc
INTEGER, DIMENSION(2*noprocs-2) :: irequest


ncount = SIZE(intloc)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'collect_vars_int_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   intprocs(:,:,:,:,1) = intloc
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
!

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Blocking
!-----------
!

IF (iform.EQ.1) THEN

   n_410: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         CALL comms_send_int(intloc,ncount,idprocs(iproc),idloc,commop,imode,&
                           & ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         CALL comms_recv_int(intprocs(:,:,:,:,iproc),ncount,idprocs(iproc),&
                           & idprocs(iproc),commop,ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         intprocs(:,:,:,:,iproc) = intloc
      ENDIF
   ENDDO n_410

!
!5. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         nreq = nreq + 1
         CALL comms_isend_int(intloc,ncount,idprocs(iproc),idloc,commop,&
                            & imode,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         nreq = nreq + 1
         CALL comms_irecv_int(intprocs(:,:,:,:,iproc),ncount,idprocs(iproc),&
                            & idprocs(iproc),commop,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         intprocs(:,:,:,:,iproc) = intloc
      ENDIF
   ENDDO n_510

!  ---wait for completion
   CALL comms_waitall(nreq,irequest)

!
!6. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_allgather_int(intloc,ncount,noprocs,intprocs,commop,ivarid)
ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_coll)


RETURN

END SUBROUTINE collect_vars_int_4d

!========================================================================

SUBROUTINE collect_vars_log_0d(logloc,logprocs,noprocs,ivarid,comm,commtype)
!************************************************************************
!
! *collect_vars_log_0d* Collect a global logical "process" array from local
!                       scalars
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_allgather_log, comms_barrier, comms_irecv_log,
!                comms_isend_log, comms_recv_log, comms_send_log,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_allgather_log, comms_barrier, comms_irecv_log, &
                   & comms_isend_log, comms_recv_log, comms_send_log, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: logloc
INTEGER, INTENT(IN) :: ivarid, noprocs
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
LOGICAL, INTENT(INOUT), DIMENSION (:) :: logprocs

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*logloc*    LOGICAL Local logical data to be collected
!*logprocs*  LOGICAL Global process array
!*noprocs*   INTEGER Number of processes in communicator
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for collect operation
!*commtype*  INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iform, imode, iproc, itype, n, npcc, nreq
LOGICAL, DIMENSION(1) :: log1d
INTEGER, DIMENSION(2*noprocs-2) :: irequest


procname(pglev+1) = 'collect_vars_log_0d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   logprocs(1) = logloc
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
!

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
log1d(1) = logloc

!
!4. Blocking
!-----------
!

IF (iform.EQ.1) THEN

   n_410: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         CALL comms_send_log(log1d,1,idprocs(iproc),idloc,commop,imode,ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         CALL comms_recv_log(logprocs(iproc:iproc),1,idprocs(iproc),&
                           & idprocs(iproc),commop,ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         logprocs(iproc) = logloc
      ENDIF
   ENDDO n_410

!
!5. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         nreq = nreq + 1
         CALL comms_isend_log(log1d,1,idprocs(iproc),idloc,commop,imode,&
                            & irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         nreq = nreq + 1
         CALL comms_irecv_log(logprocs(iproc:iproc),1,idprocs(iproc),&
                            & idprocs(iproc),commop,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         logprocs(iproc) = logloc
      ENDIF
   ENDDO n_510

!  ---wait for completion
   CALL comms_waitall(nreq,irequest)

!
!6. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_allgather_log(log1d,1,noprocs,logprocs,commop,ivarid)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_coll)

!
!7. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF


RETURN

END SUBROUTINE collect_vars_log_0d

!========================================================================

SUBROUTINE collect_vars_log_1d(logloc,logprocs,noprocs,ivarid,comm,commtype)
!************************************************************************
!
! *collect_vars_log_1d* Collect a global logical "process" array from local
!                       vector arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_allgather_log, comms_barrier, comms_irecv_log,
!                comms_isend_log, comms_recv_log, comms_send_log,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_allgather_log, comms_barrier, comms_irecv_log, &
                   & comms_isend_log, comms_recv_log, comms_send_log, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, noprocs
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
LOGICAL, INTENT(IN), DIMENSION(:) :: logloc
LOGICAL, INTENT(INOUT), DIMENSION (:,:) :: logprocs

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*logloc*    LOGICAL Local logical data to be collected
!*logprocs*  LOGICAL Global process array
!*noprocs*   INTEGER Number of processes in communicator
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for collect operation
!*commtype*  INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iform, imode, iproc, itype, n, ncount, npcc, nreq
INTEGER, DIMENSION(2*noprocs-2) :: irequest


ncount = SIZE(logloc)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'collect_vars_log_1d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   logprocs(:,1) = logloc
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
!

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Blocking
!-----------
!

IF (iform.EQ.1) THEN

   n_410: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         CALL comms_send_log(logloc,ncount,idprocs(iproc),idloc,commop,imode,&
                           & ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         CALL comms_recv_log(logprocs(:,iproc),ncount,idprocs(iproc),&
                           & idprocs(iproc),commop,ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         logprocs(:,iproc) = logloc
      ENDIF
   ENDDO n_410

!
!5. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         nreq = nreq + 1
         CALL comms_isend_log(logloc,ncount,idprocs(iproc),idloc,commop,&
                            & imode,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         nreq = nreq + 1
         CALL comms_irecv_log(logprocs(:,iproc),ncount,idprocs(iproc),&
                            & idprocs(iproc),commop,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         logprocs(:,iproc) = logloc
      ENDIF
   ENDDO n_510

!  ---wait for completion
   CALL comms_waitall(nreq,irequest)

!
!6. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_allgather_log(logloc,ncount,noprocs,logprocs,commop,ivarid)
ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_coll)


RETURN

END SUBROUTINE collect_vars_log_1d

!========================================================================

SUBROUTINE collect_vars_log_2d(logloc,logprocs,noprocs,ivarid,comm,commtype)
!************************************************************************
!
! *collect_vars_log_2d* Collect a global logical "process" array from local
!                       2-D arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_allgather_log, comms_barrier, comms_irecv_log,
!                comms_isend_log, comms_recv_log, comms_send_log,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_allgather_log, comms_barrier, comms_irecv_log, &
                   & comms_isend_log, comms_recv_log, comms_send_log, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, noprocs
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
LOGICAL, INTENT(IN), DIMENSION(:,:) :: logloc
LOGICAL, INTENT(INOUT), DIMENSION (:,:,:) :: logprocs

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*logloc*    LOGICAL Local logical data to be collected
!*logprocs*  LOGICAL Global process array
!*noprocs*   INTEGER Number of processes in communicator
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for collect operation
!*commtype*  INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iform, imode, iproc, itype, n, ncount, npcc, nreq
INTEGER, DIMENSION(2*noprocs-2) :: irequest


ncount = SIZE(logloc)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'collect_vars_log_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   logprocs(:,:,1) = logloc
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
!

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Blocking
!-----------
!

IF (iform.EQ.1) THEN

   n_410: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         CALL comms_send_log(logloc,ncount,idprocs(iproc),idloc,commop,imode,&
                           & ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         CALL comms_recv_log(logprocs(:,:,iproc),ncount,idprocs(iproc),&
                           & idprocs(iproc),commop,ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         logprocs(:,:,iproc) = logloc
      ENDIF
   ENDDO n_410

!
!5. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         nreq = nreq + 1
         CALL comms_isend_log(logloc,ncount,idprocs(iproc),idloc,commop,&
                            & imode,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         nreq = nreq + 1
         CALL comms_irecv_log(logprocs(:,:,iproc),ncount,idprocs(iproc),&
                            & idprocs(iproc),commop,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         logprocs(:,:,iproc) = logloc
      ENDIF
   ENDDO n_510

!  ---wait for completion
   CALL comms_waitall(nreq,irequest)

!
!7. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_allgather_log(logloc,ncount,noprocs,logprocs,commop,ivarid)
ENDIF

!
!8. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_coll)


RETURN

END SUBROUTINE collect_vars_log_2d

!========================================================================

SUBROUTINE collect_vars_real_0d(realloc,realprocs,noprocs,ivarid,comm,commtype)
!************************************************************************
!
! *collect_vars_real_0d* Collect a global real "process" array from local
!                        scalars
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_allgather_real, comms_barrier, comms_irecv_real,
!                comms_isend_real, comms_recv_real, comms_send_real,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_allgather_real, comms_barrier, comms_irecv_real, &
                   & comms_isend_real, comms_recv_real, comms_send_real, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, noprocs
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
REAL, INTENT(IN) :: realloc
REAL, INTENT(INOUT), DIMENSION (:) :: realprocs

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*realloc*   REAL    Local real data to be collected
!*realprocs* REAL    Global process array
!*noprocs*   INTEGER Number of processes in communicator
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for collect operation
!*commtype*  INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iform, imode, iproc, itype, n, npcc, nreq
INTEGER, DIMENSION(2*noprocs-2) :: irequest
REAL, DIMENSION(1) :: real1d


procname(pglev+1) = 'collect_vars_real_0d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   realprocs(1) = realloc
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
!

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
real1d(1) = realloc

!
!4. Blocking
!-----------
!

IF (iform.EQ.1) THEN

   n_410: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         CALL comms_send_real(real1d,1,idprocs(iproc),idloc,commop,imode,&
                            & ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         CALL comms_recv_real(realprocs(iproc:iproc),1,idprocs(iproc),&
                            & idprocs(iproc),commop,ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         realprocs(iproc) = realloc
      ENDIF
   ENDDO n_410

!
!5. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         nreq = nreq + 1
         CALL comms_isend_real(real1d,1,idprocs(iproc),idloc,commop,imode,&
                            & irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         nreq = nreq + 1
         CALL comms_irecv_real(realprocs(iproc:iproc),1,idprocs(iproc),&
                             & idprocs(iproc),commop,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         realprocs(iproc) = realloc
      ENDIF
   ENDDO n_510

!  ---wait for completion
   CALL comms_waitall(nreq,irequest)

!
!6. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_allgather_real(real1d,1,noprocs,realprocs,commop,ivarid)
ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_coll)


RETURN

END SUBROUTINE collect_vars_real_0d

!========================================================================

SUBROUTINE collect_vars_real_1d(realloc,realprocs,noprocs,ivarid,comm,commtype)
!************************************************************************
!
! *collect_vars_real_1d* Collect a global real "process" array from local
!                        vector arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_allgather_real, comms_barrier, comms_irecv_real,
!                comms_isend_real, comms_recv_real, comms_send_real,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_allgather_real, comms_barrier, comms_irecv_real, &
                   & comms_isend_real, comms_recv_real, comms_send_real, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, noprocs
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
REAL, INTENT(IN), DIMENSION(:) :: realloc
REAL, INTENT(INOUT), DIMENSION (:,:) :: realprocs

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*realloc*   REAL    Local real array to be collected
!*realprocs* REAL    Global process array
!*noprocs*   INTEGER Number of processes in communicator
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for collect operation
!*commtype*  INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iform, imode, iproc, itype, n, ncount, npcc, nreq
INTEGER, DIMENSION(2*noprocs-2) :: irequest


ncount = SIZE(realloc)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'collect_vars_real_1d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   realprocs(:,1) = realloc
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
!

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Blocking
!-----------
!

IF (iform.EQ.1) THEN

   n_410: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         CALL comms_send_real(realloc,ncount,idprocs(iproc),idloc,commop,&
                            & imode,ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         CALL comms_recv_real(realprocs(:,iproc),ncount,idprocs(iproc),&
                            & idprocs(iproc),commop,ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         realprocs(:,iproc) = realloc
      ENDIF
   ENDDO n_410

!
!5. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         nreq = nreq + 1
         CALL comms_isend_real(realloc,ncount,idprocs(iproc),idloc,commop,&
                             & imode,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         nreq = nreq + 1
         CALL comms_irecv_real(realprocs(:,iproc),ncount,idprocs(iproc),&
                             & idprocs(iproc),commop,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         realprocs(:,iproc) = realloc
      ENDIF
   ENDDO n_510

!  ---wait for completion
   CALL comms_waitall(nreq,irequest)

!
!6. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_allgather_real(realloc,ncount,noprocs,realprocs,commop,ivarid)
ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_coll)


RETURN

END SUBROUTINE collect_vars_real_1d

!========================================================================

SUBROUTINE collect_vars_real_2d(realloc,realprocs,noprocs,ivarid,comm,commtype)
!************************************************************************
!
! *collect_vars_real_2d* Collect a global real "process" array from local 2-D
!                        arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_allgather_real, comms_barrier, comms_irecv_real,
!                comms_isend_real, comms_recv_real, comms_send_real,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_allgather_real, comms_barrier, comms_irecv_real, &
                   & comms_isend_real, comms_recv_real, comms_send_real, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, noprocs
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
REAL, INTENT(IN), DIMENSION(:,:) :: realloc
REAL, INTENT(INOUT), DIMENSION (:,:,:) :: realprocs

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*realloc*   REAL    Local real array to be collected
!*realprocs* REAL    Global process array
!*noprocs*   INTEGER Number of processes in communicator
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for collect operation
!*commtype*  INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iform, imode, iproc, itype, n, ncount, npcc, nreq
INTEGER, DIMENSION(2*noprocs-2) :: irequest


ncount = SIZE(realloc)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'collect_vars_real_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   realprocs(:,:,1) = realloc
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
!

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Blocking
!-----------
!

IF (iform.EQ.1) THEN

   n_410: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         CALL comms_send_real(realloc,ncount,idprocs(iproc),idloc,commop,&
                            & imode,ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         CALL comms_recv_real(realprocs(:,:,iproc),ncount,idprocs(iproc),&
                            & idprocs(iproc),commop,ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         realprocs(:,:,iproc) = realloc
      ENDIF
   ENDDO n_410

!
!5. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         nreq = nreq + 1
         CALL comms_isend_real(realloc,ncount,idprocs(iproc),idloc,commop,&
                             & imode,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         nreq = nreq + 1
         CALL comms_irecv_real(realprocs(:,:,iproc),ncount,idprocs(iproc),&
                             & idprocs(iproc),commop,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         realprocs(:,:,iproc) = realloc
      ENDIF
   ENDDO n_510

!  ---wait for completion
   CALL comms_waitall(nreq,irequest)

!
!6. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_allgather_real(realloc,ncount,noprocs,realprocs,commop,ivarid)
ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_coll)


RETURN

END SUBROUTINE collect_vars_real_2d

!========================================================================

SUBROUTINE collect_vars_real_3d(realloc,realprocs,noprocs,ivarid,comm,commtype)
!************************************************************************
!
! *collect_vars_real_3d* Collect a global real "process" array from local 3-D
!                        arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_allgather_real, comms_barrier, comms_irecv_real,
!                comms_isend_real, comms_recv_real, comms_send_real,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_allgather_real, comms_barrier, comms_irecv_real, &
                   & comms_isend_real, comms_recv_real, comms_send_real, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, noprocs
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
REAL, INTENT(IN), DIMENSION(:,:,:) :: realloc
REAL, INTENT(INOUT), DIMENSION (:,:,:,:) :: realprocs

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*realloc*   REAL    Local real array to be collected
!*realprocs* REAL    Global process array
!*noprocs*   INTEGER Number of processes in communicator
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for collect operation
!*commtype*  INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iform, imode, iproc, itype, n, ncount, npcc, nreq
INTEGER, DIMENSION(2*noprocs-2) :: irequest


ncount = SIZE(realloc)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'collect_vars_real_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   realprocs(:,:,:,1) = realloc
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
!

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Blocking
!-----------
!

IF (iform.EQ.1) THEN

   n_410: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         CALL comms_send_real(realloc,ncount,idprocs(iproc),idloc,commop,&
                            & imode,ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         CALL comms_recv_real(realprocs(:,:,:,iproc),ncount,idprocs(iproc),&
                            & idprocs(iproc),commop,ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         realprocs(:,:,:,iproc) = realloc
      ENDIF
   ENDDO n_410

!
!5. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         nreq = nreq + 1
         CALL comms_isend_real(realloc,ncount,idprocs(iproc),idloc,commop,&
                             & imode,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         nreq = nreq + 1
         CALL comms_irecv_real(realprocs(:,:,:,iproc),ncount,idprocs(iproc),&
                             & idprocs(iproc),commop,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         realprocs(:,:,:,iproc) = realloc
      ENDIF
   ENDDO n_510

!  ---wait for completion
   CALL comms_waitall(nreq,irequest)

!
!6. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_allgather_real(realloc,ncount,noprocs,realprocs,commop,ivarid)
ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_coll)


RETURN

END SUBROUTINE collect_vars_real_3d

!========================================================================

SUBROUTINE collect_vars_real_4d(realloc,realprocs,noprocs,ivarid,comm,commtype)
!************************************************************************
!
! *collect_vars_real_4d* Collect a global real "process" array from local 4-D
!                        arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_allgather_real, comms_barrier, comms_irecv_real,
!                comms_isend_real, comms_recv_real, comms_send_real,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_allgather_real, comms_barrier, comms_irecv_real, &
                   & comms_isend_real, comms_recv_real, comms_send_real, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, noprocs
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
REAL, INTENT(IN), DIMENSION(:,:,:,:) :: realloc
REAL, INTENT(INOUT), DIMENSION (:,:,:,:,:) :: realprocs

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*realloc*   REAL    Local real array to be collected
!*realprocs* REAL    Global process array
!*noprocs*   INTEGER Number of processes in communicator
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for collect operation
!*commtype*  INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iform, imode, iproc, itype, n, ncount, npcc, nreq
INTEGER, DIMENSION(2*noprocs-2) :: irequest


ncount = SIZE(realloc)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'collect_vars_real_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   realprocs(:,:,:,:,1) = realloc
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
!

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Blocking
!-----------
!

IF (iform.EQ.1) THEN

   n_410: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         CALL comms_send_real(realloc,ncount,idprocs(iproc),idloc,commop,&
                            & imode,ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         CALL comms_recv_real(realprocs(:,:,:,:,iproc),ncount,idprocs(iproc),&
                            & idprocs(iproc),commop,ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         realprocs(:,:,:,:,iproc) = realloc
      ENDIF
   ENDDO n_410

!
!5. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         nreq = nreq + 1
         CALL comms_isend_real(realloc,ncount,idprocs(iproc),idloc,commop,&
                             & imode,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.2) THEN
         nreq = nreq + 1
         CALL comms_irecv_real(realprocs(:,:,:,:,iproc),ncount,idprocs(iproc),&
                             & idprocs(iproc),commop,irequest(nreq),ivarid)
      ELSEIF (comprocs(n).EQ.0) THEN
         realprocs(:,:,:,:,iproc) = realloc
      ENDIF
   ENDDO n_510

!  ---wait for completion
   CALL comms_waitall(nreq,irequest)

!
!6. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_allgather_real(realloc,ncount,noprocs,realprocs,commop,ivarid)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_coll)

!
!7. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF


RETURN

END SUBROUTINE collect_vars_real_4d

!========================================================================

SUBROUTINE combine_mod_hgrid_2d(surfglb,surfloc,lbounds,ivarid,ifill,&
                              & rfill,idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_mod_hgrid_2d* Combine local derived type 2-D model arrays
!                        of type 'HRelativeCoords' into a global array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.9
!
! Description -
!
! Module calls - combine_mod_int_2d, combine_mod_real_2d, error_alloc
!
!************************************************************************
!
USE datatypes
USE error_routines, ONLY: error_alloc

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ifill, ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(IN) :: rfill
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds
TYPE (HRelativeCoords), INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):) :: &
                                                                      & surfglb
TYPE (HRelativeCoords), INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):) :: &
                                                                      & surfloc

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*surfglb*  DERIVED Global real grid array
!*surfloc*  DERIVED Distributed local array
!*lbounds*  INTEGER Lower bounds of arrays
!*ivarid*   INTEGER Variable id
!*ifill*    INTEGER Fill value for integer components
!*rfill*    REAL    Fill value for real components
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for distribute operation
!*commtype* INTEGER Type of communication (1/4)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, idsrce, itype, n1, n2
INTEGER, DIMENSION(2) :: uboundsglb, uboundsloc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: weightsglb, weightsloc


procname(pglev+1) = 'combine_mod_hgrid_2d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   surfglb = surfloc
   GOTO 1000
ENDIF

!
!2. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!3. Allocate and storage
!-----------------------
!

uboundsglb = UBOUND(surfglb)
ALLOCATE(weightsglb(lbounds(1):uboundsglb(1),lbounds(2):uboundsglb(2),2,2),&
       & STAT=errstat)
CALL error_alloc('weightsglb',4,(/uboundsglb(1)-lbounds(1)+1,&
                                & uboundsglb(2)-lbounds(2)+1,2,2/),kndrtype)

uboundsloc = UBOUND(surfloc)
ALLOCATE(weightsloc(lbounds(1):uboundsloc(1),lbounds(2):uboundsloc(2),2,2),&
       & STAT=errstat)
CALL error_alloc('weightsloc',4,(/uboundsloc(1)-lbounds(1)+1,&
                                & uboundsloc(2)-lbounds(2)+1,2,2/),kndrtype)
n2_310: DO n2=lbounds(2),uboundsloc(2) 
n1_310: DO n1=lbounds(1),uboundsloc(1) 
   weightsloc(n1,n2,:,:) = surfloc(n1,n2)%weights
ENDDO n1_310
ENDDO n2_310

!
!4. Distribute
!-------------
!

CALL combine_mod_int_2d(surfglb%icoord,surfloc%icoord,lbounds,ivarid,ifill,0,&
                      & idsrce,commop,itype,all)
CALL combine_mod_int_2d(surfglb%jcoord,surfloc%jcoord,lbounds,ivarid,ifill,0,&
                      & idsrce,commop,itype,all)
CALL combine_mod_real_4d(weightsglb,weightsloc,(/lbounds(1),lbounds(2),1,1/),&
                       & ivarid,rfill,0,idsrce,commop,itype,all)
n2_410: DO n2=lbounds(2),uboundsglb(2) 
n1_410: DO n1=lbounds(1),uboundsglb(1) 
   surfglb(n1,n2)%weights = weightsglb(n1,n2,:,:)
ENDDO n1_410
ENDDO n2_410

!
!5. Deallocate
!-------------
!

DEALLOCATE (weightsglb,weightsloc)

1000 CALL log_timer_out()


RETURN

END SUBROUTINE combine_mod_hgrid_2d

!========================================================================

SUBROUTINE combine_mod_int_2d(intglb,intloc,lbounds,ivarid,ifill,mglevel,&
                            & idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_mod_int_2d* Combine local integer 2-D model arrays into a global
!                      array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.11
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_waitall
!
!************************************************************************
!
USE grid
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ifill, ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds
INTEGER, INTENT(OUT), DIMENSION(lbounds(1):,lbounds(2):) :: intglb
INTEGER, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):) :: intloc

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intglb*   INTEGER Global integer grid array
!*intloc*   INTEGER Local integer grid array
!*lbounds*  INTEGER Lower bounds of arrays
!*ivarid*   INTEGER Variable id
!*ifill*    INTEGER Fill value for global array
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for combine operation
!*commtype* INTEGER Type of communication (1/4)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, idlocx, iform, imode, iproc, itype, i1, i2, j1, j2, &
         & lev, n, norcv, nosnd, npcc, nreq, nxloc, nx1loc, nx2loc, nyloc, &
         & ny1loc, ny2loc
INTEGER, DIMENSION(nprocscoh) :: nxprocs, nx1procs, nx2procs, nyprocs, &
                               & ny1procs, ny2procs
INTEGER, DIMENSION(npworld) :: idprocsx
INTEGER, DIMENSION(2*nprocs-2) :: irequest


procname(pglev+1) = 'combine_mod_int_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   intglb = intloc
   GOTO 1000
ENDIF

!
!2. Initialise global array
!--------------------------
!

intglb = ifill

!
!3. Optional arguments
!---------------------
!

IF (PRESENT(mglevel)) THEN
   lev = mglevel
ELSE
   lev = 0
ENDIF

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!4. Initialise parameters
!------------------------
!
!---(multi-)grid parameters
IF (lev.EQ.0) THEN
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
   nxprocs = ncprocs; nyprocs = nrprocs
   nx1procs = nc1procs; ny1procs = nr1procs
   nx2procs = nc2procs; ny2procs = nr2procs
ELSE
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
   nxprocs = mgvars(lev)%ncprocs; nyprocs = mgvars(lev)%nrprocs
   nx1procs = mgvars(lev)%nc1procs; ny1procs = mgvars(lev)%nr1procs
   nx2procs = mgvars(lev)%nc2procs; ny2procs = mgvars(lev)%nr2procs
ENDIF

!---process ids within communicator
IF (commop.EQ.comm_world_MPI) THEN
   idlocx = idlocglb
   idprocsx = idprocsglb
ELSE
   idlocx = idloc
   idprocsx(1:nprocs) = idprocs
ENDIF

!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
nosnd = nxloc*nyloc

!
!5. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!  ---send
   IF (idlocx.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_int(intloc(1:nxloc,1:nyloc),nosnd,iddest,idlocx,&
                           & commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_int(intloc(1:nxloc,1:nyloc),nosnd,iddest,idlocx,&
                            & commop,imode,irequest(nreq),ivarid)
      ENDIF
!  ---receive
   ELSE
      iproc_510: DO iproc=1,nprocscoh
         i1 = nx1procs(iproc); i2 = nx2procs(iproc)
         j1 = ny1procs(iproc); j2 = ny2procs(iproc)
         IF (idprocsx(iproc).EQ.iddest) THEN
            intglb(i1:i2,j1:j2) = intloc(1:nxloc,1:nyloc)
         ELSE
            norcv = nxprocs(iproc)*nyprocs(iproc)
            IF (iform.EQ.1) THEN
               CALL comms_recv_int(intglb(i1:i2,j1:j2),norcv,&
                                 & idprocsx(iproc),idprocsx(iproc),commop,&
                                 & ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_int(intglb(i1:i2,j1:j2),norcv,&
                                  & idprocsx(iproc),idprocsx(iproc),commop,&
                                  & irequest(nreq),ivarid)
            ENDIF
         ENDIF
      ENDDO iproc_510
   ENDIF

!
!6. Result on all processes
!--------------------------
!

ELSE

   n_610: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
!     ---send
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_int(intloc(1:nxloc,1:nyloc),nosnd,&
                              & idprocsx(iproc),idlocx,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_int(intloc(1:nxloc,1:nyloc),nosnd,&
                               & idprocsx(iproc),idlocx,commop,imode,&
                               & irequest(nreq),ivarid)
         ENDIF
!     ---receive
      ELSEIF (comprocs(n).EQ.2) THEN
         i1 = nx1procs(iproc); i2 = nc2procs(iproc)
         j1 = ny1procs(iproc); j2 = ny2procs(iproc)
         norcv = nxprocs(iproc)*nyprocs(iproc)
         IF (iform.EQ.1) THEN
            CALL comms_recv_int(intglb(i1:i2,j1:j2),norcv,idprocsx(iproc),&
                              & idprocsx(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_int(intglb(i1:i2,j1:j2),norcv,&
                               & idprocsx(iproc),idprocsx(iproc),commop,&
                               & irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         intglb(nc1loc:nc2loc,nr1loc:nr2loc) = intloc(1:ncloc,1:nrloc)
      ENDIF
   ENDDO n_610

ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_mod_int_2d

!========================================================================

SUBROUTINE combine_mod_int_3d(intglb,intloc,lbounds,ivarid,ifill,mglevel,&
                            & idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_mod_int_3d* Combine local integer 3-D model arrays into a global
!                      array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.11
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_waitall, error_abort
!
!************************************************************************
!
USE grid
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ifill, ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds
INTEGER, INTENT(OUT), DIMENSION(lbounds(1):,lbounds(2):,lbounds(3):) :: intglb
INTEGER, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):,lbounds(3):) :: intloc

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intglb*   INTEGER Global integer grid array
!*intloc*   INTEGER Local integer grid array
!*lbounds*  INTEGER Lower bounds of arrays
!*ivarid*   INTEGER Variable id
!*ifill*    INTEGER Fill value for global array
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for combine operation
!*commtype* INTEGER Type of communication (1/4)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, idlocx, iform, imode, iproc, itype, i1, i2, j1, j2, &
         & lev, n, norcv, nosnd, novals, npcc, nreq, nxloc, nx1loc, nx2loc, &
         & nyloc, ny1loc, ny2loc
INTEGER, DIMENSION(nprocscoh) :: nxprocs, nx1procs, nx2procs, nyprocs, &
                               & ny1procs, ny2procs
INTEGER, DIMENSION(npworld) :: idprocsx
INTEGER, DIMENSION(2*nprocs-2) :: irequest


procname(pglev+1) = 'combine_mod_int_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   intglb = intloc
   GOTO 1000
ENDIF

!
!2. Initialise global array
!--------------------------
!

intglb = ifill

!
!3. Optional arguments
!---------------------
!

IF (PRESENT(mglevel)) THEN
   lev = mglevel
ELSE
   lev = 0
ENDIF

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!4. Initialise parameters
!------------------------
!
!---(multi-)grid parameters
IF (lev.EQ.0) THEN
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
   nxprocs = ncprocs; nyprocs = nrprocs
   nx1procs = nc1procs; ny1procs = nr1procs
   nx2procs = nc2procs; ny2procs = nr2procs
ELSE
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
   nxprocs = mgvars(lev)%ncprocs; nyprocs = mgvars(lev)%nrprocs
   nx1procs = mgvars(lev)%nc1procs; ny1procs = mgvars(lev)%nr1procs
   nx2procs = mgvars(lev)%nc2procs; ny2procs = mgvars(lev)%nr2procs
ENDIF

!---process ids within communicator
IF (commop.EQ.comm_world_MPI) THEN
   idlocx = idlocglb
   idprocsx = idprocsglb
ELSE
   idlocx = idloc
   idprocsx(1:nprocs) = idprocs
ENDIF

!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
novals = SIZE(intloc,DIM=3)
nosnd = nxloc*nyloc*novals

!
!5. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!  ---send
   IF (idlocx.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_int(intloc(1:nxloc,1:nyloc,:),nosnd,iddest,idlocx,&
                           & commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_int(intloc(1:nxloc,1:nyloc,:),nosnd,iddest,idlocx,&
                            & commop,imode,irequest(nreq),ivarid)
      ENDIF
!  ---receive
   ELSE
      iproc_510: DO iproc=1,nprocscoh
         i1 = nx1procs(iproc); i2 = nx2procs(iproc)
         j1 = ny1procs(iproc); j2 = ny2procs(iproc)
         IF (idprocs(iproc).EQ.iddest) THEN
            intglb(i1:i2,j1:j2,:) = intloc(1:ncloc,1:nrloc,:)
         ELSE
            norcv = nxprocs(iproc)*nyprocs(iproc)
            IF (iform.EQ.1) THEN
               CALL comms_recv_int(intglb(i1:i2,j1:j2,:),norcv,&
                                & idprocsx(iproc),idprocsx(iproc),commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_int(intglb(i1:i2,j1:j2,:),norcv,&
                                  & idprocsx(iproc),idprocsx(iproc),commop,&
                                  & irequest(nreq),ivarid)
            ENDIF
         ENDIF
      ENDDO iproc_510
   ENDIF

!
!6. Result on all processes
!--------------------------
!

ELSE

   n_610: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
!     ---send
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_int(intloc(1:nxloc,1:nyloc,:),nosnd,&
                              & idprocsx(iproc),idlocx,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_int(intloc(1:nxloc,1:nyloc,:),nosnd,&
                               & idprocsx(iproc),idlocx,commop,imode,&
                               & irequest(nreq),ivarid)
         ENDIF
!     ---receive
      ELSEIF (comprocs(n).EQ.2) THEN
         i1 = nx1procs(iproc); i2 = nx2procs(iproc)
         j1 = ny1procs(iproc); j2 = ny2procs(iproc)
         norcv = nxprocs(iproc)*nyprocs(iproc)*novals
         IF (iform.EQ.1) THEN
            CALL comms_recv_int(intglb(i1:i2,j1:j2,:),norcv,idprocsx(iproc),&
                              & idprocsx(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_int(intglb(i1:i2,j1:j2,:),norcv,&
                               & idprocsx(iproc),idprocsx(iproc),commop,&
                               & irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         intglb(nx1loc:nx2loc,ny1loc:ny2loc,:) = intloc(1:nxloc,1:nyloc,:)
      ENDIF
   ENDDO n_610

ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_mod_int_3d

!========================================================================

SUBROUTINE combine_mod_int_4d(intglb,intloc,lbounds,ivarid,ifill,mglevel,&
                            & idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_mod_int_4d* Combine local integer 4-D model arrays into a global
!                      array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.11
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_waitall
!
!************************************************************************
!
USE grid
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ifill, ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
INTEGER, INTENT(IN), DIMENSION(4) :: lbounds
INTEGER, INTENT(OUT), DIMENSION(lbounds(1):,lbounds(2):,&
                              & lbounds(3):,lbounds(4):) :: intglb
INTEGER, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):,&
                             & lbounds(3):,lbounds(4):) :: intloc

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intglb*   INTEGER Global integer grid array
!*intloc*   INTEGER Local integer grid array
!*lbounds*  INTEGER Lower bounds of arrays
!*ivarid*   INTEGER Variable id
!*ifill*    INTEGER Fill value for global array
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for combine operation
!*commtype* INTEGER Type of communication (1/4)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, idlocx, iform, imode, iproc, itype, i1, i2, j1, j2, &
         & lev, n, norcv, nosnd, novals, npcc, nreq, nxloc, nx1loc, nx2loc, &
         & nyloc, ny1loc, ny2loc
INTEGER, DIMENSION(nprocscoh) :: nxprocs, nx1procs, nx2procs, nyprocs, &
                               & ny1procs, ny2procs
INTEGER, DIMENSION(npworld) :: idprocsx
INTEGER, DIMENSION(2*nprocs-2) :: irequest


procname(pglev+1) = 'combine_mod_int_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   intglb = intloc
   GOTO 1000
ENDIF

!
!2. Initialise global array
!--------------------------
!

intglb = ifill

!
!3. Optional arguments
!---------------------
!

IF (PRESENT(mglevel)) THEN
   lev = mglevel
ELSE
   lev = 0
ENDIF

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!4. Initialise parameters
!------------------------
!
!---(multi-)grid parameters
IF (lev.EQ.0) THEN
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
   nxprocs = ncprocs; nyprocs = nrprocs
   nx1procs = nc1procs; ny1procs = nr1procs
   nx2procs = nc2procs; ny2procs = nr2procs
ELSE
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
   nxprocs = mgvars(lev)%ncprocs; nyprocs = mgvars(lev)%nrprocs
   nx1procs = mgvars(lev)%nc1procs; ny1procs = mgvars(lev)%nr1procs
   nx2procs = mgvars(lev)%nc2procs; ny2procs = mgvars(lev)%nr2procs
ENDIF

!---process ids within communicator
IF (commop.EQ.comm_world_MPI) THEN
   idlocx = idlocglb
   idprocsx = idprocsglb
ELSE
   idlocx = idloc
   idprocsx(1:nprocs) = idprocs
ENDIF

!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
novals = SIZE(intloc,DIM=3)*SIZE(intloc,DIM=4)
nosnd = nxloc*nyloc*novals

!
!5. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!  ---send
   IF (idlocx.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_int(intloc(1:nxloc,1:nyloc,:,:),nosnd,iddest,idlocx,&
                           & commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_int(intloc(1:nxloc,1:nyloc,:,:),nosnd,iddest,idlocx,&
                            & commop,imode,irequest(nreq),ivarid)
      ENDIF
!  ---receive
   ELSE
      iproc_510: DO iproc=1,nprocscoh
         i1 = nx1procs(iproc); i2 = nx2procs(iproc)
         j1 = ny1procs(iproc); j2 = ny2procs(iproc)
         IF (idprocs(iproc).EQ.iddest) THEN
            intglb(i1:i2,j1:j2,:,:) = intloc(1:ncloc,1:nrloc,:,:)
         ELSE
            norcv = nxprocs(iproc)*nyprocs(iproc)*novals
            IF (iform.EQ.1) THEN
               CALL comms_recv_int(intglb(i1:i2,j1:j2,:,:),norcv,&
                                & idprocsx(iproc),idprocsx(iproc),commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_int(intglb(i1:i2,j1:j2,:,:),norcv,&
                                  & idprocsx(iproc),idprocsx(iproc),commop,&
                                  & irequest(nreq),ivarid)
            ENDIF
         ENDIF
      ENDDO iproc_510
   ENDIF

!
!6. Result on all processes
!--------------------------
!

ELSE

   n_610: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
!     ---send
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_int(intloc(1:nxloc,1:nyloc,:,:),nosnd,&
                              & idprocsx(iproc),idlocx,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_int(intloc(1:nxloc,1:nyloc,:,:),nosnd,&
                               & idprocsx(iproc),idlocx,commop,imode,&
                               & irequest(nreq),ivarid)
         ENDIF
!     ---receive
      ELSEIF (comprocs(n).EQ.2) THEN
         i1 = nx1procs(iproc); i2 = nx2procs(iproc)
         j1 = ny1procs(iproc); j2 = ny2procs(iproc)
         norcv = nxprocs(iproc)*nyprocs(iproc)*novals
         IF (iform.EQ.1) THEN
            CALL comms_recv_int(intglb(i1:i2,j1:j2,:,:),norcv,idprocsx(iproc),&
                              & idprocsx(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_int(intglb(i1:i2,j1:j2,:,:),norcv,&
                               & idprocsx(iproc),idprocsx(iproc),commop,&
                               & irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         intglb(nx1loc:nx2loc,ny1loc:ny2loc,:,:) = intloc(1:nxloc,1:nyloc,:,:)
      ENDIF
   ENDDO n_610

ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_mod_int_4d

!========================================================================

SUBROUTINE combine_mod_log_2d(logglb,logloc,lbounds,ivarid,lfill,mglevel,&
                            & idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_mod_log_2d* Combine local logical 2-D model arrays into a global
!                      array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.11
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_log, comms_isend_log,
!                comms_recv_log, comms_send_log, comms_waitall, error_abort
!
!************************************************************************
!
USE grid
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_log, comms_isend_log, &
                   & comms_recv_log, comms_send_log, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
LOGICAL, INTENT(IN) :: lfill
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds
LOGICAL, INTENT(OUT), DIMENSION(lbounds(1):,lbounds(2):) :: logglb
LOGICAL, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):) :: logloc

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*logglb*   LOGICAL Global logical grid array
!*logloc*   LOGICAL Local logical grid array
!*lbounds*  INTEGER Lower bounds of arrays
!*ivarid*   INTEGER Variable id
!*lfill*    LOGICAL Fill value for global array
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for combine operation
!*commtype* INTEGER Type of communication (1/4)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, idlocx, iform, imode, iproc, itype, i1, i2, j1, j2, &
         & lev, n, npcc, norcv, nosnd, nreq, nxloc, nx1loc, nx2loc, nyloc, &
         & ny1loc, ny2loc
INTEGER, DIMENSION(nprocscoh) :: nxprocs, nx1procs, nx2procs, nyprocs, &
                               & ny1procs, ny2procs
INTEGER, DIMENSION(npworld) :: idprocsx
INTEGER, DIMENSION(2*nprocs-2) :: irequest


procname(pglev+1) = 'combine_mod_log_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   logglb = logloc
   GOTO 1000
ENDIF

!
!2. Initialise global array
!--------------------------
!

logglb = lfill

!
!3. Optional arguments
!---------------------
!

IF (PRESENT(mglevel)) THEN
   lev = mglevel
ELSE
   lev = 0
ENDIF

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!4. Initialise parameters
!------------------------
!
!---(multi-)grid parameters
IF (lev.EQ.0) THEN
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
   nxprocs = ncprocs; nyprocs = nrprocs
   nx1procs = nc1procs; ny1procs = nr1procs
   nx2procs = nc2procs; ny2procs = nr2procs
ELSE
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
   nxprocs = mgvars(lev)%ncprocs; nyprocs = mgvars(lev)%nrprocs
   nx1procs = mgvars(lev)%nc1procs; ny1procs = mgvars(lev)%nr1procs
   nx2procs = mgvars(lev)%nc2procs; ny2procs = mgvars(lev)%nr2procs
ENDIF

!---process ids within communicator
IF (commop.EQ.comm_world_MPI) THEN
   idlocx = idlocglb
   idprocsx = idprocsglb
ELSE
   idlocx = idloc
   idprocsx(1:nprocs) = idprocs
ENDIF

!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
nosnd = nxloc*nyloc

!
!5. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!  ---send
   IF (idlocx.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_log(logloc(1:nxloc,1:nyloc),nosnd,iddest,idlocx,&
                           & commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_log(logloc(1:nxloc,1:nyloc),nosnd,iddest,idlocx,&
                            & commop,imode,irequest(nreq),ivarid)
      ENDIF
!  ---receive
   ELSE
      iproc_510: DO iproc=1,nprocscoh
         i1 = nx1procs(iproc); i2 = nx2procs(iproc)
         j1 = ny1procs(iproc); j2 = ny2procs(iproc)
         IF (idprocsx(iproc).EQ.iddest) THEN
            logglb(i1:i2,j1:j2) = logloc(1:nxloc,1:nyloc)
         ELSE
            norcv = nxprocs(iproc)*nyprocs(iproc)
            IF (iform.EQ.1) THEN
               CALL comms_recv_log(logglb(i1:i2,j1:j2),norcv,&
                                & idprocsx(iproc),idprocsx(iproc),commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_log(logglb(i1:i2,j1:j2),norcv,&
                                  & idprocsx(iproc),idprocsx(iproc),commop,&
                                  & irequest(nreq),ivarid)
            ENDIF
         ENDIF
      ENDDO iproc_510
   ENDIF

!
!6. Result on all processes
!--------------------------
!

ELSE

   n_610: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
!     ---send
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_log(logloc(1:nxloc,1:nyloc),nosnd,&
                              & idprocsx(iproc),idlocx,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_log(logloc(1:nxloc,1:nyloc),nosnd,&
                               & idprocsx(iproc),idlocx,commop,imode,&
                               & irequest(nreq),ivarid)
         ENDIF
!     ---receive
      ELSEIF (comprocs(n).EQ.2) THEN
         i1 = nx1procs(iproc); i2 = nx2procs(iproc)
         j1 = ny1procs(iproc); j2 = ny2procs(iproc)
         norcv = nxprocs(iproc)*nyprocs(iproc)
         IF (iform.EQ.1) THEN
            CALL comms_recv_log(logglb(i1:i2,j1:j2),norcv,idprocsx(iproc),&
                              & idprocsx(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_log(logglb(i1:i2,j1:j2),norcv,&
                               & idprocsx(iproc),idprocsx(iproc),commop,&
                               & irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         logglb(nx1loc:nx2loc,ny1loc:ny2loc) = logloc(1:nxloc,1:nyloc)
      ENDIF
   ENDDO n_610

ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_mod_log_2d

!========================================================================

SUBROUTINE combine_mod_real_2d(realglb,realloc,lbounds,ivarid,rfill,mglevel,&
                             & idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_mod_real_2d* Combine local real 2-D model arrays into a global array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.11
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE grid
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
REAL, INTENT(IN) :: rfill
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds
REAL, INTENT(OUT), DIMENSION(lbounds(1):,lbounds(2):) :: realglb
REAL, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):) :: realloc

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realglb*  REAL    Global real grid array
!*realloc*  REAL    Local real grid array
!*lbounds*  INTEGER Lower bounds of arrays
!*ivarid*   INTEGER Variable id
!*rfill*    REAL    Fill value for global array
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for combine operation
!*commtype* INTEGER Type of communication (1/4)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, idlocx, iform, imode, iproc, itype, i1, i2, j1, j2, &
         & lev, n, norcv, nosnd, npcc, nreq, nxloc, nx1loc, nx2loc, nyloc, &
         & ny1loc, ny2loc
INTEGER, DIMENSION(nprocscoh) :: nxprocs, nx1procs, nx2procs, nyprocs, &
                               & ny1procs, ny2procs
INTEGER, DIMENSION(npworld) :: idprocsx
INTEGER, DIMENSION(2*nprocs-2) :: irequest


procname(pglev+1) = 'combine_mod_real_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   realglb = realloc
   GOTO 1000
ENDIF

!
!2. Initialise global array
!--------------------------
!

realglb = rfill

!
!3. Optional arguments
!---------------------
!

IF (PRESENT(mglevel)) THEN
   lev = mglevel
ELSE
   lev = 0
ENDIF

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!4. Initialise parameters
!------------------------
!
!---(multi-)grid parameters
IF (lev.EQ.0) THEN
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
   nxprocs = ncprocs; nyprocs = nrprocs
   nx1procs = nc1procs; ny1procs = nr1procs
   nx2procs = nc2procs; ny2procs = nr2procs
ELSE
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
   nxprocs = mgvars(lev)%ncprocs; nyprocs = mgvars(lev)%nrprocs
   nx1procs = mgvars(lev)%nc1procs; ny1procs = mgvars(lev)%nr1procs
   nx2procs = mgvars(lev)%nc2procs; ny2procs = mgvars(lev)%nr2procs
ENDIF

!---process ids within communicator
IF (commop.EQ.comm_world_MPI) THEN
   idlocx = idlocglb
   idprocsx = idprocsglb
ELSE
   idlocx = idloc
   idprocsx(1:nprocs) = idprocs
ENDIF

!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
nosnd = nxloc*nyloc

!
!5. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!  ---send
   IF (idlocx.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_real(realloc(1:nxloc,1:nyloc),nosnd,iddest,&
                            & idlocx,commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_real(realloc(1:nxloc,1:nyloc),nosnd,iddest,&
                             & idlocx,commop,imode,irequest(nreq),ivarid)
      ENDIF

!  ---receive
   ELSE
      iproc_510: DO iproc=1,nprocscoh
         i1 = nx1procs(iproc); i2 = nx2procs(iproc)
         j1 = ny1procs(iproc); j2 = ny2procs(iproc)
         IF (idprocsx(iproc).EQ.iddest) THEN
            realglb(i1:i2,j1:j2) = realloc(1:ncloc,1:nrloc)
         ELSE
            norcv = nxprocs(iproc)*nyprocs(iproc)
            IF (iform.EQ.1) THEN
               CALL comms_recv_real(realglb(i1:i2,j1:j2),norcv,&
                                  & idprocsx(iproc),idprocsx(iproc),commop,&
                                  & ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_real(realglb(i1:i2,j1:j2),norcv,&
                                   & idprocsx(iproc),idprocsx(iproc),commop,&
                                   & irequest(nreq),ivarid)
            ENDIF
         ENDIF
      ENDDO iproc_510

   ENDIF

!
!6. Result on all processes
!--------------------------
!

ELSE

   n_610: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
!     ---send
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_real(realloc(1:nxloc,1:nyloc),nosnd,&
                               & idprocsx(iproc),idlocx,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_real(realloc(1:nxloc,1:nyloc),nosnd,&
                                & idprocsx(iproc),idlocx,commop,imode,&
                                & irequest(nreq),ivarid)
         ENDIF
!     ---receive
      ELSEIF (comprocs(n).EQ.2) THEN
         i1 = nx1procs(iproc); i2 = nx2procs(iproc)
         j1 = ny1procs(iproc); j2 = ny2procs(iproc)
         norcv = nxprocs(iproc)*nyprocs(iproc)
         IF (iform.EQ.1) THEN
            CALL comms_recv_real(realglb(i1:i2,j1:j2),norcv,&
                               & idprocsx(iproc),idprocsx(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_real(realglb(i1:i2,j1:j2),norcv,&
                                & idprocsx(iproc),idprocsx(iproc),commop,&
                                & irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         realglb(nx1loc:nx2loc,ny1loc:ny2loc) = realloc(1:nxloc,1:nyloc)
      ENDIF
   ENDDO n_610

ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_mod_real_2d

!========================================================================

SUBROUTINE combine_mod_real_3d(realglb,realloc,lbounds,ivarid,rfill,mglevel,&
                             & idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_mod_real_3d* Combine local real 3-D model arrays into a global array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.11
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE grid
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
REAL, INTENT(IN) :: rfill
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds
REAL, INTENT(OUT), DIMENSION(lbounds(1):,lbounds(2):,lbounds(3):) :: realglb
REAL, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):,lbounds(3):) :: realloc

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realglb*  REAL    Global real grid array
!*realloc*  REAL    Local real grid array
!*lbounds*  INTEGER Lower bounds of arrays
!*ivarid*   INTEGER Variable id
!*rfill*    REAL    Fill value for global array
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for combine operation
!*commtype* INTEGER Type of communication (1/4)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, idlocx, iform, imode, iproc, itype, i1, i2, j1, j2, &
         & lev, n, norcv, nosnd, novals, npcc, nreq, nxloc, nx1loc, nx2loc, &
         & nyloc, ny1loc, ny2loc
INTEGER, DIMENSION(nprocscoh) :: nxprocs, nx1procs, nx2procs, nyprocs, &
                               & ny1procs, ny2procs
INTEGER, DIMENSION(npworld) :: idprocsx
INTEGER, DIMENSION(2*nprocs-2) :: irequest


procname(pglev+1) = 'combine_mod_real_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   realglb = realloc
   GOTO 1000
ENDIF

!
!2. Initialise global array
!--------------------------
!

realglb = rfill

!
!3. Optional arguments
!---------------------
!

IF (PRESENT(mglevel)) THEN
   lev = mglevel
ELSE
   lev = 0
ENDIF

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!4. Initialise parameters and arrays
!-----------------------------------
!
!---(multi-)grid parameters
IF (lev.EQ.0) THEN
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
   nxprocs = ncprocs; nyprocs = nrprocs
   nx1procs = nc1procs; ny1procs = nr1procs
   nx2procs = nc2procs; ny2procs = nr2procs
ELSE
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
   nxprocs = mgvars(lev)%ncprocs; nyprocs = mgvars(lev)%nrprocs
   nx1procs = mgvars(lev)%nc1procs; ny1procs = mgvars(lev)%nr1procs
   nx2procs = mgvars(lev)%nc2procs; ny2procs = mgvars(lev)%nr2procs
ENDIF

!---process ids within communicator
IF (commop.EQ.comm_world_MPI) THEN
   idlocx = idlocglb
   idprocsx = idprocsglb
ELSE
   idlocx = idloc
   idprocsx(1:nprocs) = idprocs
ENDIF

!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
novals = SIZE(realloc,DIM=3)
nosnd = nxloc*nyloc*novals

!
!5. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!  ---send
   IF (idlocx.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_real(realloc(1:nxloc,1:nyloc,:),nosnd,iddest,&
                            & idlocx,commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_real(realloc(1:nxloc,1:nyloc,:),nosnd,iddest,&
                             & idlocx,commop,imode,irequest(nreq),ivarid)
      ENDIF
!  ---receive
   ELSE
      iproc_510: DO iproc=1,nprocscoh
         i1 = nx1procs(iproc); i2 = nx2procs(iproc)
         j1 = ny1procs(iproc); j2 = ny2procs(iproc)
         norcv = nxprocs(iproc)*nyprocs(iproc)*novals
         IF (idprocsx(iproc).EQ.iddest) THEN
            realglb(i1:i2,j1:j2,:) = realloc(1:nxloc,1:nyloc,:)
         ELSE
            IF (iform.EQ.1) THEN
               CALL comms_recv_real(realglb(i1:i2,j1:j2,:),norcv,&
                                  & idprocsx(iproc),idprocsx(iproc),commop,&
                                  & ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_real(realglb(i1:i2,j1:j2,:),norcv,&
                                   & idprocsx(iproc),idprocsx(iproc),commop,&
                                   & irequest(nreq),ivarid)
            ENDIF
         ENDIF
      ENDDO iproc_510
   ENDIF

!
!6. Result on all processes
!--------------------------
!

ELSE

   n_610: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
!     ---send
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_real(realloc(1:nxloc,1:nyloc,:),nosnd,&
                               & idprocsx(iproc),idlocx,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_real(realloc(1:nxloc,1:nyloc,:),nosnd,&
                                & idprocsx(iproc),idlocx,commop,imode,&
                                & irequest(nreq),ivarid)
         ENDIF
!     ---receive
      ELSEIF (comprocs(n).EQ.2) THEN
         i1 = nx1procs(iproc); i2 = nx2procs(iproc)
         j1 = ny1procs(iproc); j2 = ny2procs(iproc)
         norcv = nxprocs(iproc)*nyprocs(iproc)*novals
         IF (iform.EQ.1) THEN
            CALL comms_recv_real(realglb(i1:i2,j1:j2,:),norcv,&
                               & idprocsx(iproc),idprocsx(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_real(realglb(i1:i2,j1:j2,:),norcv,&
                                & idprocsx(iproc),idprocsx(iproc),commop,&
                                & irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         realglb(nx1loc:nx2loc,ny1loc:ny2loc,:) = realloc(1:nxloc,1:nyloc,:)
      ENDIF
   ENDDO n_610

ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_mod_real_3d

!========================================================================

SUBROUTINE combine_mod_real_4d(realglb,realloc,lbounds,ivarid,rfill,mglevel,&
                             & idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_mod_real_4d* Combine local real 4-D model arrays into a global array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.11
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE grid
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
REAL, INTENT(IN) :: rfill
INTEGER, INTENT(IN), DIMENSION(4) :: lbounds
REAL, INTENT(OUT), DIMENSION(lbounds(1):,lbounds(2):,&
                           & lbounds(3):,lbounds(4):) :: realglb
REAL, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):,&
                          & lbounds(3):,lbounds(4):) :: realloc

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realglb*  REAL    Global real grid array
!*realloc*  REAL    Local real grid array
!*lbounds*  INTEGER Lower bounds of arrays
!*ivarid*   INTEGER Variable id
!*rfill*    REAL    Fill value for global array
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for combine operation
!*commtype* INTEGER Type of communication (1/4)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, idlocx, iform, imode, iproc, itype, i1, i2, j1, j2, &
         & lev, n, norcv, nosnd, novals, npcc, nreq, nxloc, nx1loc, nx2loc, &
         & nyloc, ny1loc, ny2loc
INTEGER, DIMENSION(nprocscoh) :: nxprocs, nx1procs, nx2procs, nyprocs, &
                               & ny1procs, ny2procs
INTEGER, DIMENSION(npworld) :: idprocsx
INTEGER, DIMENSION(2*nprocs-2) :: irequest


procname(pglev+1) = 'combine_mod_real_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   realglb = realloc
   GOTO 1000
ENDIF

!
!2. Initialise global array
!--------------------------
!

realglb = rfill

!
!3. Optional arguments
!---------------------
!

IF (PRESENT(mglevel)) THEN
   lev = mglevel
ELSE
   lev = 0
ENDIF

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!4. Initialise parameters
!------------------------
!
!---(multi-)grid parameters
IF (lev.EQ.0) THEN
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
   nxprocs = ncprocs; nyprocs = nrprocs
   nx1procs = nc1procs; ny1procs = nr1procs
   nx2procs = nc2procs; ny2procs = nr2procs
ELSE
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
   nxprocs = mgvars(lev)%ncprocs; nyprocs = mgvars(lev)%nrprocs
   nx1procs = mgvars(lev)%nc1procs; ny1procs = mgvars(lev)%nr1procs
   nx2procs = mgvars(lev)%nc2procs; ny2procs = mgvars(lev)%nr2procs
ENDIF

!---process ids within communicator
IF (commop.EQ.comm_world_MPI) THEN
   idlocx = idlocglb
   idprocsx = idprocsglb
ELSE
   idlocx = idloc
   idprocsx(1:nprocs) = idprocs
ENDIF

!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
novals = SIZE(realloc,DIM=3)*SIZE(realloc,DIM=4)
nosnd = nxloc*nyloc*novals

!
!5. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!  ---send
   IF (idlocx.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_real(realloc(1:nxloc,1:nyloc,:,:),nosnd,iddest,&
                            & idlocx,commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_real(realloc(1:nxloc,1:nyloc,:,:),nosnd,iddest,&
                             & idlocx,commop,imode,irequest(nreq),ivarid)
      ENDIF
!  ---receive
   ELSE
      iproc_510: DO iproc=1,nprocscoh
         i1 = nx1procs(iproc); i2 = nx2procs(iproc)
         j1 = ny1procs(iproc); j2 = ny2procs(iproc)
         IF (idprocsx(iproc).EQ.iddest) THEN
            realglb(i1:i2,j1:j2,:,:) = realloc(1:nxloc,1:nyloc,:,:)
         ELSE
            norcv = nxprocs(iproc)*nyprocs(iproc)*novals
            IF (iform.EQ.1) THEN
               CALL comms_recv_real(realglb(i1:i2,j1:j2,:,:),norcv,&
                                  & idprocsx(iproc),idprocsx(iproc),commop,&
                                  & ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_real(realglb(i1:i2,j1:j2,:,:),norcv,&
                                   & idprocsx(iproc),idprocsx(iproc),commop,&
                                   & irequest(nreq),ivarid)
            ENDIF
         ENDIF
      ENDDO iproc_510
   ENDIF

!
!6. Result on all processes
!--------------------------
!

ELSE

   n_610: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
!     ---send
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_real(realloc(1:nxloc,1:nyloc,:,:),nosnd,&
                               & idprocsx(iproc),idlocx,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_real(realloc(1:nxloc,1:nyloc,:,:),nosnd,&
                                & idprocsx(iproc),idlocx,commop,imode,&
                                & irequest(nreq),ivarid)
         ENDIF
!     ---receive
      ELSEIF (comprocs(n).EQ.2) THEN
         i1 = nx1procs(iproc); i2 = nx2procs(iproc)
         j1 = ny1procs(iproc); j2 = ny2procs(iproc)
         norcv = nxprocs(iproc)*nyprocs(iproc)*novals
         IF (iform.EQ.1) THEN
            CALL comms_recv_real(realglb(i1:i2,j1:j2,:,:),norcv,&
                               & idprocsx(iproc),idprocsx(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_real(realglb(i1:i2,j1:j2,:,:),norcv,&
                                & idprocsx(iproc),idprocsx(iproc),commop,&
                                & irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         realglb(nx1loc:nx2loc,ny1loc:ny2loc,:,:) = &
      &  realloc(1:nxloc,1:nyloc,:,:)
      ENDIF
   ENDDO n_610

ENDIF

!
!7. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_comb)

RETURN

END SUBROUTINE combine_mod_real_4d

!========================================================================

SUBROUTINE combine_stats_glb_int_1d(intglb,maxstats,nostatprocs,lstatprocs,&
                                  & ivarid,idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_stats_glb_int_1d* Combine global integer vector data defined at
!                            stations on different domains
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.7
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_int,
!                comms_isend_int, comms_recv_int, comms_send_int,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid, maxstats
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
INTEGER, INTENT(INOUT), DIMENSION(:) :: intglb

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*intglb*      INTEGER Global integer array
!*maxstats*    INTEGER First dimension of array lstatprocs
!*nostatprocs* INTEGER Number of stations in each local array
!*lstatprocs*  INTEGER Indices of local stations in global array
!*ivarid*      INTEGER Variable id
!*idroot*      INTEGER id of root process
!*comm*        INTEGER Communicator for combine operation
!*commtype*    INTEGER Type of communication (1/4)
!*commall*     LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, l, lloc, n, norcv, &
         & nosnd, npcc, nreq
INTEGER, DIMENSION(2*nprocs-2) :: irequest
INTEGER, DIMENSION(SIZE(intglb)) :: idatrcv, idatsnd


IF (SIZE(intglb).EQ.0) RETURN  
procname(pglev+1) = 'combine_stats_glb_int_1d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Store data in local array
!----------------------------
!

IF (all.OR.idloc.NE.iddest) THEN
   iproc_310: DO iproc=1,nprocs
      IF (idloc.EQ.idprocs(iproc)) THEN
         nosnd = nostatprocs(iproc)
         lloc_311: DO lloc=1,nosnd
            l = lstatprocs(lloc,iproc)
            idatsnd(lloc) = intglb(l)
         ENDDO lloc_311
      ENDIF
   ENDDO iproc_310
ENDIF

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!4.1 Send
!--------
!

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_int(idatsnd(1:nosnd),nosnd,iddest,idloc,commop,&
                           & imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         IF (nosnd.GT.0) nreq = nreq + 1
         CALL comms_isend_int(idatsnd(1:nosnd),nosnd,iddest,idloc,commop,&
                            & imode,irequest(nreq),ivarid)
      ENDIF
   ENDIF

!
!4.2 Receive
!-----------
!

   IF (idloc.EQ.iddest) THEN
      iproc_420: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.iddest) THEN
!           ---receive
            norcv = nostatprocs(iproc)
            IF (iform.EQ.1) THEN 
               CALL comms_recv_int(idatrcv(1:norcv),norcv,idprocs(iproc),&
                                 & idprocs(iproc),commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               IF (norcv.GT.0) nreq = nreq + 1
               CALL comms_irecv_int(idatrcv(1:norcv),norcv,idprocs(iproc),&
                                  & idprocs(iproc),commop,irequest(nreq),&
                                  & ivarid)
            ENDIF
!           ---store
            lloc_421: DO lloc=1,norcv
               l = lstatprocs(lloc,iproc)
               intglb(l) = idatrcv(lloc)
            ENDDO lloc_421
         ENDIF
      ENDDO iproc_420
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE
   n_500: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1

!
!5.1 Send
!--------
!

      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_int(idatsnd(1:nosnd),nosnd,idprocs(iproc),idloc,&
                              & commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (nosnd.GT.0) nreq = nreq + 1
            CALL comms_isend_int(idatsnd(1:nosnd),nosnd,idprocs(iproc),idloc,&
                               & commop,imode,irequest(nreq),ivarid)
         ENDIF

!
!5.2 Receive
!-----------
!

      ELSEIF (comprocs(n).EQ.2) THEN
!        ---receive
         norcv = nostatprocs(iproc)
         IF (iform.EQ.1) THEN
            CALL comms_recv_int(idatrcv(1:norcv),norcv,idprocs(iproc),&
                              & idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (norcv.GT.0) nreq = nreq + 1
            CALL comms_irecv_int(idatrcv(1:norcv),norcv,idprocs(iproc),&
                               & idprocs(iproc),commop,irequest(nreq),ivarid)
         ENDIF
!           ---store
         lloc_521: DO lloc=1,norcv
            l = lstatprocs(lloc,iproc)
            intglb(l) = idatrcv(lloc)
         ENDDO lloc_521
      ENDIF

   ENDDO n_500

ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_stats_glb_int_1d

!========================================================================

SUBROUTINE combine_stats_glb_int_2d(intglb,maxstats,nostatprocs,lstatprocs,&
                                  & ivarid,idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_stats_glb_int_2d* Combine global integer 2-D data defined at
!                            stations on different domains
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.6
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_int,
!                comms_isend_int, comms_recv_int, comms_send_int,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid, maxstats
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
INTEGER, INTENT(INOUT), DIMENSION(:,:) :: intglb

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*intglb*      INTEGER Global integer array
!*maxstats*    INTEGER First dimension of array lstatprocs
!*nostatprocs* INTEGER Number of stations in each local array
!*lstatprocs*  INTEGER Indices of local stations in global array
!*ivarid*      INTEGER Variable id
!*idroot*      INTEGER id of root process
!*comm*        INTEGER Communicator for combine operation
!*commtype*    INTEGER Type of communication (1/4)
!*commall*     LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, l, lloc, n, norcv, &
         & nosnd, npcc, nreq, numdim
INTEGER, DIMENSION(2*nprocs-2) :: irequest
INTEGER, DIMENSION(SIZE(intglb,DIM=1),SIZE(intglb,DIM=2)) :: idatrcv, idatsnd


IF (SIZE(intglb).EQ.0) RETURN  
procname(pglev+1) = 'combine_stats_glb_int_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
numdim = SIZE(intglb,DIM=2)

!
!3. Store data in local array
!----------------------------
!

IF (all.OR.idloc.NE.iddest) THEN
   iproc_310: DO iproc=1,nprocs
      IF (idloc.EQ.idprocs(iproc)) THEN
         nosnd = nostatprocs(iproc)
         lloc_311: DO lloc=1,nosnd
            l = lstatprocs(lloc,iproc)
            idatsnd(lloc,:) = intglb(l,:)
         ENDDO lloc_311
      ENDIF
   ENDDO iproc_310
ENDIF

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!4.1 Send
!--------
!

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_int(idatsnd(1:nosnd,:),nosnd*numdim,iddest,idloc,&
                           & commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         IF (nosnd.GT.0) nreq = nreq + 1
         CALL comms_isend_int(idatsnd(1:nosnd,:),nosnd*numdim,iddest,idloc,&
                            & commop,imode,irequest(nreq),ivarid)
      ENDIF
   ENDIF

!
!4.2 Receive
!-----------
!

   IF (idloc.EQ.iddest) THEN
      iproc_420: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.iddest) THEN
!           ---receive
            norcv = nostatprocs(iproc)
            IF (iform.EQ.1) THEN 
               CALL comms_recv_int(idatrcv(1:norcv,:),norcv*numdim,&
                                 & idprocs(iproc),idprocs(iproc),commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               IF (norcv.GT.0) nreq = nreq + 1
               CALL comms_irecv_int(idatrcv(1:norcv,:),norcv*numdim,&
                                  & idprocs(iproc),idprocs(iproc),commop,&
                                  & irequest(nreq),ivarid)
            ENDIF
!           ---store
            lloc_421: DO lloc=1,norcv
               l = lstatprocs(lloc,iproc)
               intglb(l,:) = idatrcv(lloc,:)
            ENDDO lloc_421
         ENDIF
      ENDDO iproc_420
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE
   n_500: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1

!
!5.1 Send
!--------
!

      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_int(idatsnd(1:nosnd,:),nosnd*numdim,idprocs(iproc),&
                              & idloc,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (nosnd.GT.0) nreq = nreq + 1
            CALL comms_isend_int(idatsnd(1:nosnd,:),nosnd*numdim,&
                               & idprocs(iproc),idloc,commop,imode,&
                               & irequest(nreq),ivarid)
         ENDIF

!
!5.2 Receive
!-----------
!

      ELSEIF (comprocs(n).EQ.2) THEN
!        ---receive
         norcv = nostatprocs(iproc)
         IF (iform.EQ.1) THEN
            CALL comms_recv_int(idatrcv(1:norcv,:),norcv*numdim,idprocs(iproc),&
                              & idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (norcv.GT.0) nreq = nreq + 1
            CALL comms_irecv_int(idatrcv(1:norcv,:),norcv*numdim,&
                               & idprocs(iproc),idprocs(iproc),commop,&
                               & irequest(nreq),ivarid)
         ENDIF
!           ---store
         lloc_521: DO lloc=1,norcv
            l = lstatprocs(lloc,iproc)
            intglb(l,:) = idatrcv(lloc,:)
         ENDDO lloc_521
      ENDIF

   ENDDO n_500

ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_stats_glb_int_2d

!========================================================================

SUBROUTINE combine_stats_glb_int_3d(intglb,maxstats,nostatprocs,lstatprocs,&
                                  & ivarid,idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_stats_glb_int_3d* Combine global integer 3-D data defined at
!                            stations on different domains
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.6
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_int,
!                comms_isend_int, comms_recv_int, comms_send_int,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid, maxstats
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
INTEGER, INTENT(INOUT), DIMENSION(:,:,:) :: intglb

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*intglb*      INTEGER Global integer array
!*maxstats*    INTEGER First dimension of array lstatprocs
!*nostatprocs* INTEGER Number of stations in each local array
!*lstatprocs*  INTEGER Indices of local stations in global array
!*ivarid*      INTEGER Variable id
!*idroot*      INTEGER id of root process
!*comm*        INTEGER Communicator for combine operation
!*commtype*    INTEGER Type of communication (1/4)
!*commall*     LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, l, lloc, n, norcv, &
         & nosnd, npcc, nreq, numdim
INTEGER, DIMENSION(2*nprocs-2) :: irequest
INTEGER, DIMENSION(SIZE(intglb,DIM=1),SIZE(intglb,DIM=2),&
                 & SIZE(intglb,DIM=3)) :: idatrcv, idatsnd


IF (SIZE(intglb).EQ.0) RETURN  
procname(pglev+1) = 'combine_stats_glb_int_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
numdim = SIZE(intglb,DIM=2)*SIZE(intglb,DIM=3)

!
!3. Store data in local array
!----------------------------
!

IF (all.OR.idloc.NE.iddest) THEN
   iproc_310: DO iproc=1,nprocs
      IF (idloc.EQ.idprocs(iproc)) THEN
         nosnd = nostatprocs(iproc)
         lloc_311: DO lloc=1,nosnd
            l = lstatprocs(lloc,iproc)
            idatsnd(lloc,:,:) = intglb(l,:,:)
         ENDDO lloc_311
      ENDIF
   ENDDO iproc_310
ENDIF

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!4.1 Send
!--------
!

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_int(idatsnd(1:nosnd,:,:),nosnd*numdim,iddest,idloc,&
                           & commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         IF (nosnd.GT.0) nreq = nreq + 1
         CALL comms_isend_int(idatsnd(1:nosnd,:,:),nosnd*numdim,iddest,idloc,&
                            & commop,imode,irequest(nreq),ivarid)
      ENDIF
   ENDIF

!
!4.2 Receive
!-----------
!

   IF (idloc.EQ.iddest) THEN
      iproc_420: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.iddest) THEN
!           ---receive
            norcv = nostatprocs(iproc)
            IF (iform.EQ.1) THEN 
               CALL comms_recv_int(idatrcv(1:norcv,:,:),norcv*numdim,&
                                 & idprocs(iproc),idprocs(iproc),commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               IF (norcv.GT.0) nreq = nreq + 1
               CALL comms_irecv_int(idatrcv(1:norcv,:,:),norcv*numdim,&
                                  & idprocs(iproc),idprocs(iproc),commop,&
                                  & irequest(nreq),ivarid)
            ENDIF
!           ---store
            lloc_421: DO lloc=1,norcv
               l = lstatprocs(lloc,iproc)
               intglb(l,:,:) = idatrcv(lloc,:,:)
            ENDDO lloc_421
         ENDIF
      ENDDO iproc_420
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE
   n_500: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1

!
!5.1 Send
!--------
!

      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_int(idatsnd(1:nosnd,:,:),nosnd*numdim,&
                              & idprocs(iproc),idloc,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (nosnd.GT.0) nreq = nreq + 1
            CALL comms_isend_int(idatsnd(1:nosnd,:,:),nosnd*numdim,&
                               & idprocs(iproc),idloc,commop,imode,&
                               & irequest(nreq),ivarid)
         ENDIF

!
!5.2 Receive
!-----------
!

      ELSEIF (comprocs(n).EQ.2) THEN
!        ---receive
         norcv = nostatprocs(iproc)
         IF (iform.EQ.1) THEN
            CALL comms_recv_int(idatrcv(1:norcv,:,:),norcv*numdim,&
                              & idprocs(iproc),idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (norcv.GT.0) nreq = nreq + 1
            CALL comms_irecv_int(idatrcv(1:norcv,:,:),norcv*numdim,&
                               & idprocs(iproc),idprocs(iproc),commop,&
                               & irequest(nreq),ivarid)
         ENDIF
!           ---store
         lloc_521: DO lloc=1,norcv
            l = lstatprocs(lloc,iproc)
            intglb(l,:,:) = idatrcv(lloc,:,:)
         ENDDO lloc_521
      ENDIF

   ENDDO n_500

ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_stats_glb_int_3d

!========================================================================

SUBROUTINE combine_stats_glb_int_4d(intglb,maxstats,nostatprocs,lstatprocs,&
                                  & ivarid,idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_stats_glb_int_4d* Combine global integer 4-D data defined at
!                            stations on different domains
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.6
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_int,
!                comms_isend_int, comms_recv_int, comms_send_int,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid, maxstats
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
INTEGER, INTENT(INOUT), DIMENSION(:,:,:,:) :: intglb

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*intglb*      INTEGER Global integer array
!*maxstats*    INTEGER First dimension of array lstatprocs
!*nostatprocs* INTEGER Number of stations in each local array
!*lstatprocs*  INTEGER Indices of local stations in global array
!*ivarid*      INTEGER Variable id
!*idroot*      INTEGER id of root process
!*comm*        INTEGER Communicator for combine operation
!*commtype*    INTEGER Type of communication (1/4)
!*commall*     LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, l, lloc, n, norcv, &
         & nosnd, npcc, nreq, numdim
INTEGER, DIMENSION(2*nprocs-2) :: irequest
INTEGER, DIMENSION(SIZE(intglb,DIM=1),SIZE(intglb,DIM=2),SIZE(intglb,DIM=3),&
                 & SIZE(intglb,DIM=4)) :: idatrcv, idatsnd


IF (SIZE(intglb).EQ.0) RETURN  
procname(pglev+1) = 'combine_stats_glb_int_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
numdim = SIZE(intglb,DIM=2)*SIZE(intglb,DIM=3)*SIZE(intglb,DIM=4)

!
!3. Store data in local array
!----------------------------
!

IF (all.OR.idloc.NE.iddest) THEN
   iproc_310: DO iproc=1,nprocs
      IF (idloc.EQ.idprocs(iproc)) THEN
         nosnd = nostatprocs(iproc)
         lloc_311: DO lloc=1,nosnd
            l = lstatprocs(lloc,iproc)
            idatsnd(lloc,:,:,:) = intglb(l,:,:,:)
         ENDDO lloc_311
      ENDIF
   ENDDO iproc_310
ENDIF

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!4.1 Send
!--------
!

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_int(idatsnd(1:nosnd,:,:,:),nosnd*numdim,iddest,idloc,&
                           & commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         IF (nosnd.GT.0) nreq = nreq + 1
         CALL comms_isend_int(idatsnd(1:nosnd,:,:,:),nosnd*numdim,iddest,idloc,&
                            & commop,imode,irequest(nreq),ivarid)
      ENDIF
   ENDIF

!
!4.2 Receive
!-----------
!

   IF (idloc.EQ.iddest) THEN
      iproc_420: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.iddest) THEN
!           ---receive
            norcv = nostatprocs(iproc)
            IF (iform.EQ.1) THEN 
               CALL comms_recv_int(idatrcv(1:norcv,:,:,:),norcv*numdim,&
                                 & idprocs(iproc),idprocs(iproc),commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               IF (norcv.GT.0) nreq = nreq + 1
               CALL comms_irecv_int(idatrcv(1:norcv,:,:,:),norcv*numdim,&
                                  & idprocs(iproc),idprocs(iproc),commop,&
                                  & irequest(nreq),ivarid)
            ENDIF
!           ---store
            lloc_421: DO lloc=1,norcv
               l = lstatprocs(lloc,iproc)
               intglb(l,:,:,:) = idatrcv(lloc,:,:,:)
            ENDDO lloc_421
         ENDIF
      ENDDO iproc_420
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE
   n_500: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1

!
!5.1 Send
!--------
!

      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_int(idatsnd(1:nosnd,:,:,:),nosnd*numdim,&
                              & idprocs(iproc),idloc,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (nosnd.GT.0) nreq = nreq + 1
            CALL comms_isend_int(idatsnd(1:nosnd,:,:,:),nosnd*numdim,&
                               & idprocs(iproc),idloc,commop,imode,&
                               & irequest(nreq),ivarid)
         ENDIF

!
!5.2 Receive
!-----------
!

      ELSEIF (comprocs(n).EQ.2) THEN
!        ---receive
         norcv = nostatprocs(iproc)
         IF (iform.EQ.1) THEN
            CALL comms_recv_int(idatrcv(1:norcv,:,:,:),norcv*numdim,&
                              & idprocs(iproc),idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (norcv.GT.0) nreq = nreq + 1
            CALL comms_irecv_int(idatrcv(1:norcv,:,:,:),norcv*numdim,&
                               & idprocs(iproc),idprocs(iproc),commop,&
                               & irequest(nreq),ivarid)
         ENDIF
!           ---store
         lloc_521: DO lloc=1,norcv
            l = lstatprocs(lloc,iproc)
            intglb(l,:,:,:) = idatrcv(lloc,:,:,:)
         ENDDO lloc_521
      ENDIF

   ENDDO n_500

ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_stats_glb_int_4d

!========================================================================

SUBROUTINE combine_stats_glb_log_1d(logglb,maxstats,nostatprocs,lstatprocs,&
                                  & ivarid,idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_stats_glb_log_1d* Combine global logical vector data defined at
!                            stations on different domains
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.6
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_log,
!                comms_isend_log, comms_recv_log, comms_send_log,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_log, comms_isend_log, &
                   & comms_recv_log, comms_send_log, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid, maxstats
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
LOGICAL, INTENT(INOUT), DIMENSION(:) :: logglb

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*logglb*      LOGICAL Global logical array
!*maxstats*    INTEGER First dimension of array lstatprocs
!*nostatprocs* INTEGER Number of stations in each local array
!*lstatprocs*  INTEGER Indices of local stations in global array
!*ivarid*      INTEGER Variable id
!*idroot*      INTEGER id of root process
!*comm*        INTEGER Communicator for combine operation
!*commtype*    INTEGER Type of communication (1/4)
!*commall*     LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, l, lloc, n, norcv, &
         & nosnd, npcc, nreq
INTEGER, DIMENSION(2*nprocs-2) :: irequest
LOGICAL, DIMENSION(SIZE(logglb)) :: ldatrcv, ldatsnd


IF (SIZE(logglb).EQ.0) RETURN  
procname(pglev+1) = 'combine_stats_glb_log_1d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Store data in local array
!----------------------------
!

IF (all.OR.idloc.NE.iddest) THEN
   iproc_310: DO iproc=1,nprocs
      IF (idloc.EQ.idprocs(iproc)) THEN
         nosnd = nostatprocs(iproc)
         lloc_311: DO lloc=1,nosnd
            l = lstatprocs(lloc,iproc)
            ldatsnd(lloc) = logglb(l)
         ENDDO lloc_311
      ENDIF
   ENDDO iproc_310
ENDIF

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!4.1 Send
!--------
!

   IF (idloc.NE.iddest.AND.nosnd.GT.0) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_log(ldatsnd(1:nosnd),nosnd,iddest,idloc,commop,&
                           & imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         IF (nosnd.GT.0) nreq = nreq + 1
         CALL comms_isend_log(ldatsnd(1:nosnd),nosnd,iddest,idloc,commop,&
                            & imode,irequest(nreq),ivarid)
      ENDIF
   ENDIF

!
!4.2 Receive
!-----------
!

   IF (idloc.EQ.iddest) THEN
      iproc_420: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.iddest) THEN
!           ---receive
            norcv = nostatprocs(iproc)
            IF (iform.EQ.1) THEN 
               CALL comms_recv_log(ldatrcv(1:norcv),norcv,idprocs(iproc),&
                                 & idprocs(iproc),commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               IF (norcv.GT.0) nreq = nreq + 1
               CALL comms_irecv_log(ldatrcv(1:norcv),norcv,idprocs(iproc),&
                                  & idprocs(iproc),commop,irequest(nreq),ivarid)
            ENDIF
!           ---store
            lloc_421: DO lloc=1,norcv
               l = lstatprocs(lloc,iproc)
               logglb(l) = ldatrcv(lloc)
            ENDDO lloc_421
         ENDIF
      ENDDO iproc_420
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE
   n_500: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1

!
!5.1 Send
!--------
!

      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_log(ldatsnd(1:nosnd),nosnd,idprocs(iproc),idloc,&
                              & commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (nosnd.GT.0) nreq = nreq + 1
            CALL comms_isend_log(ldatsnd(1:nosnd),nosnd,idprocs(iproc),&
                               & idloc,commop,imode,irequest(nreq),ivarid)
         ENDIF

!
!5.2 Receive
!-----------
!

      ELSEIF (comprocs(n).EQ.2) THEN
!        ---receive
         norcv = nostatprocs(iproc)
         IF (iform.EQ.1) THEN
            CALL comms_recv_log(ldatrcv(1:norcv),norcv,idprocs(iproc),&
                              & idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (norcv.GT.0) nreq = nreq + 1
            CALL comms_irecv_log(ldatrcv(1:norcv),norcv,idprocs(iproc),&
                               & idprocs(iproc),commop,irequest(nreq),ivarid)
         ENDIF
!        ---store
         lloc_521: DO lloc=1,norcv
            l = lstatprocs(lloc,iproc)
            logglb(l) = ldatrcv(lloc)
         ENDDO lloc_521

      ENDIF

   ENDDO n_500

ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_stats_glb_log_1d

!========================================================================

SUBROUTINE combine_stats_glb_real_1d(realglb,maxstats,nostatprocs,lstatprocs,&
                                   & ivarid,idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_stats_glb_real_1d* Combine global real vector data defined at
!                             stations on different domains
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.6
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_real,
!                comms_isend_real, comms_recv_real, comms_send_real,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid, maxstats
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
REAL, INTENT(INOUT), DIMENSION(:) :: realglb

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realglb*     REAL    Global real array
!*maxstats*    INTEGER First dimension of array lstatprocs
!*nostatprocs* INTEGER Number of stations in each local array
!*lstatprocs*  INTEGER Indices of local stations in global array
!*ivarid*      INTEGER Variable id
!*idroot*      INTEGER id of root process
!*comm*        INTEGER Communicator for combine operation
!*commtype*    INTEGER Type of communication (1/4)
!*commall*     LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, l, lloc, n, norcv, &
         & nosnd, npcc, nreq
INTEGER, DIMENSION(2*nprocs-2) :: irequest
REAL, DIMENSION(SIZE(realglb)) :: rdatrcv, rdatsnd


IF (SIZE(realglb).EQ.0) RETURN  
procname(pglev+1) = 'combine_stats_glb_real_1d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Store data in local array
!----------------------------
!

IF (all.OR.idloc.NE.iddest) THEN
   iproc_310: DO iproc=1,nprocs
      IF (idloc.EQ.idprocs(iproc)) THEN
         nosnd = nostatprocs(iproc)
         lloc_311: DO lloc=1,nosnd
            l = lstatprocs(lloc,iproc)
            rdatsnd(lloc) = realglb(l)
         ENDDO lloc_311
      ENDIF
   ENDDO iproc_310
ENDIF

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!4.1 Send
!--------
!

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_real(rdatsnd(1:nosnd),nosnd,iddest,idloc,commop,&
                            & imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         IF (nosnd.GT.0) nreq = nreq + 1
         CALL comms_isend_real(rdatsnd(1:nosnd),nosnd,iddest,idloc,commop,&
                             & imode,irequest(nreq),ivarid)
      ENDIF
   ENDIF

!
!4.2 Receive
!-----------
!

   IF (idloc.EQ.iddest) THEN
      iproc_420: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.iddest) THEN
!           ---receive
            norcv = nostatprocs(iproc)
            IF (iform.EQ.1) THEN 
               CALL comms_recv_real(rdatrcv(1:norcv),norcv,idprocs(iproc),&
                                  & idprocs(iproc),commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               IF (norcv.GT.0) nreq = nreq + 1
               CALL comms_irecv_real(rdatrcv(1:norcv),norcv,idprocs(iproc),&
                                   & idprocs(iproc),commop,irequest(nreq),&
                                   & ivarid)
            ENDIF
!           ---store
            lloc_421: DO lloc=1,norcv
               l = lstatprocs(lloc,iproc)
               realglb(l) = rdatrcv(lloc)
            ENDDO lloc_421
         ENDIF
      ENDDO iproc_420
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE
   n_500: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1

!
!5.1 Send
!--------
!

      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_real(rdatsnd(1:nosnd),nosnd,idprocs(iproc),idloc,&
                               & commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (nosnd.GT.0) nreq = nreq + 1
            CALL comms_isend_real(rdatsnd(1:nosnd),nosnd,idprocs(iproc),idloc,&
                                & commop,imode,irequest(nreq),ivarid)
         ENDIF

!
!5.2 Receive
!-----------
!

      ELSEIF (comprocs(n).EQ.2) THEN
!        ---receive
         norcv = nostatprocs(iproc)
         IF (iform.EQ.1) THEN 
            CALL comms_recv_real(rdatrcv(1:norcv),norcv,idprocs(iproc),&
                               & idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (norcv.GT.0) nreq = nreq + 1
            CALL comms_irecv_real(rdatrcv(1:norcv),norcv,idprocs(iproc),&
                                & idprocs(iproc),commop,irequest(nreq),ivarid)
         ENDIF
!        ---store
         lloc_521: DO lloc=1,norcv
            l = lstatprocs(lloc,iproc)
            realglb(l) = rdatrcv(lloc)
         ENDDO lloc_521
      ENDIF

   ENDDO n_500

ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_stats_glb_real_1d

!========================================================================

SUBROUTINE combine_stats_glb_real_2d(realglb,maxstats,nostatprocs,lstatprocs,&
                                   & ivarid,idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_stats_glb_real_2d* Combine global real 2-D data defined at
!                             stations on different domains
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.6
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!             - array elements with the same first index are defined on the
!               same domain on input
!
! Module calls - comms_barrier, comms_irecv_real,
!                comms_isend_real, comms_recv_real, comms_send_real,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid, maxstats
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
REAL, INTENT(INOUT), DIMENSION(:,:) :: realglb

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realglb*     REAL    Global real array
!*maxstats*    INTEGER First dimension of array lstatprocs
!*nostatprocs* INTEGER Number of stations in each local array
!*lstatprocs*  INTEGER Indices of local stations in global array
!*ivarid*      INTEGER Variable id
!*idroot*      INTEGER id of root process
!*comm*        INTEGER Communicator for combine operation
!*commtype*    INTEGER Type of communication (1/4)
!*commall*     LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, l, lloc, n, norcv, &
         & nosnd, npcc, nreq, numdim
INTEGER, DIMENSION(2*nprocs-2) :: irequest
REAL, DIMENSION(SIZE(realglb,DIM=1),SIZE(realglb,DIM=2)) :: rdatrcv, rdatsnd


IF (SIZE(realglb).EQ.0) RETURN  
procname(pglev+1) = 'combine_stats_glb_real_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
numdim = SIZE(realglb,DIM=2)

!
!3. Store data in local array
!----------------------------
!

IF (all.OR.idloc.NE.iddest) THEN
   iproc_310: DO iproc=1,nprocs
      IF (idloc.EQ.idprocs(iproc)) THEN
         nosnd = nostatprocs(iproc)
         lloc_311: DO lloc=1,nosnd
            l = lstatprocs(lloc,iproc)
            rdatsnd(lloc,:) = realglb(l,:)
         ENDDO lloc_311
      ENDIF
   ENDDO iproc_310
ENDIF

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!4.1 Send
!--------
!

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_real(rdatsnd(1:nosnd,:),nosnd*numdim,iddest,idloc,&
                            & commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         IF (nosnd.GT.0) nreq = nreq + 1
         CALL comms_isend_real(rdatsnd(1:nosnd,:),nosnd*numdim,iddest,idloc,&
                             & commop,imode,irequest(nreq),ivarid)
      ENDIF
   ENDIF

!
!4.2 Receive
!-----------
!

   IF (idloc.EQ.iddest) THEN
      iproc_420: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.iddest) THEN
!           ---receive
            norcv = nostatprocs(iproc)
            IF (iform.EQ.1) THEN 
               CALL comms_recv_real(rdatrcv(1:norcv,:),norcv*numdim,&
                                  & idprocs(iproc),idprocs(iproc),commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               IF (norcv.GT.0) nreq = nreq + 1
               CALL comms_irecv_real(rdatrcv(1:norcv,:),norcv*numdim,&
                                   & idprocs(iproc),idprocs(iproc),commop,&
                                   & irequest(nreq),ivarid)
            ENDIF
!           ---store
            lloc_421: DO lloc=1,norcv
               l = lstatprocs(lloc,iproc)
               realglb(l,:) = rdatrcv(lloc,:)
            ENDDO lloc_421
         ENDIF
      ENDDO iproc_420
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE
   n_500: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1

!
!5.1 Send
!--------
!

      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_real(rdatsnd(1:nosnd,:),nosnd*numdim,&
                               & idprocs(iproc),idloc,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (nosnd.GT.0) nreq = nreq + 1
            CALL comms_isend_real(rdatsnd(1:nosnd,:),nosnd*numdim,&
                                & idprocs(iproc),idloc,commop,imode,&
                                & irequest(nreq),ivarid)
         ENDIF

!
!5.2 Receive
!-----------
!

      ELSEIF (comprocs(n).EQ.2) THEN
!        ---receive
         norcv = nostatprocs(iproc)
         IF (iform.EQ.1) THEN 
            CALL comms_recv_real(rdatrcv(1:norcv,:),norcv*numdim,&
                               & idprocs(iproc),idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (norcv.GT.0) nreq = nreq + 1
            CALL comms_irecv_real(rdatrcv(1:norcv,:),norcv*numdim,&
                                & idprocs(iproc),idprocs(iproc),commop,&
                                & irequest(nreq),ivarid)
         ENDIF
!           ---store
         lloc_521: DO lloc=1,norcv
            l = lstatprocs(lloc,iproc)
            realglb(l,:) = rdatrcv(lloc,:)
         ENDDO lloc_521
      ENDIF

   ENDDO n_500

ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_stats_glb_real_2d

!========================================================================

SUBROUTINE combine_stats_glb_real_3d(realglb,maxstats,nostatprocs,lstatprocs,&
                                   & ivarid,idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_stats_glb_real_3d* Combine global real 3-D data defined at
!                             stations on different domains
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.6
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!             - array elements with the same first index are defined on the
!               same domain on input
!
! Module calls - comms_barrier, comms_irecv_real,
!                comms_isend_real, comms_recv_real, comms_send_real,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid, maxstats
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
REAL, INTENT(INOUT), DIMENSION(:,:,:) :: realglb

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realglb*     REAL    Global real array
!*maxstats*    INTEGER First dimension of array lstatprocs
!*nostatprocs* INTEGER Number of stations in each local array
!*lstatprocs*  INTEGER Indices of local stations in global array
!*ivarid*      INTEGER Variable id
!*idroot*      INTEGER id of root process
!*comm*        INTEGER Communicator for combine operation
!*commtype*    INTEGER Type of communication (1/4)
!*commall*     LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, l, lloc, n, norcv, &
         & nosnd, npcc, nreq, numdim
INTEGER, DIMENSION(2*nprocs-2) :: irequest
REAL, DIMENSION(SIZE(realglb,DIM=1),SIZE(realglb,DIM=2),&
              & SIZE(realglb,DIM=3)) :: rdatrcv, rdatsnd


IF (SIZE(realglb).EQ.0) RETURN  
procname(pglev+1) = 'combine_stats_glb_real_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
numdim = SIZE(realglb,DIM=2)*SIZE(realglb,DIM=3)

!
!3. Store data in local array
!----------------------------
!

IF (all.OR.idloc.NE.iddest) THEN
   iproc_310: DO iproc=1,nprocs
      IF (idloc.EQ.idprocs(iproc)) THEN
         nosnd = nostatprocs(iproc)
         lloc_311: DO lloc=1,nosnd
            l = lstatprocs(lloc,iproc)
            rdatsnd(lloc,:,:) = realglb(l,:,:)
         ENDDO lloc_311
      ENDIF
   ENDDO iproc_310
ENDIF

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!4.1 Send
!--------
!

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_real(rdatsnd(1:nosnd,:,:),nosnd*numdim,iddest,idloc,&
                            & commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         IF (nosnd.GT.0) nreq = nreq + 1
         CALL comms_isend_real(rdatsnd(1:nosnd,:,:),nosnd*numdim,iddest,idloc,&
                             & commop,imode,irequest(nreq),ivarid)
      ENDIF
   ENDIF

!
!4.2 Receive
!-----------
!

   IF (idloc.EQ.iddest) THEN
      iproc_420: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.iddest) THEN
!           ---receive
            norcv = nostatprocs(iproc)
            IF (iform.EQ.1) THEN 
               CALL comms_recv_real(rdatrcv(1:norcv,:,:),norcv*numdim,&
                                  & idprocs(iproc),idprocs(iproc),commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               IF (norcv.GT.0) nreq = nreq + 1
               CALL comms_irecv_real(rdatrcv(1:norcv,:,:),norcv*numdim,&
                                   & idprocs(iproc),idprocs(iproc),commop,&
                                   & irequest(nreq),ivarid)
            ENDIF
!           ---store
            lloc_421: DO lloc=1,norcv
               l = lstatprocs(lloc,iproc)
               realglb(l,:,:) = rdatrcv(lloc,:,:)
            ENDDO lloc_421
         ENDIF
      ENDDO iproc_420
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE
   n_500: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1

!
!5.1 Send
!--------
!

      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_real(rdatsnd(1:nosnd,:,:),nosnd*numdim,&
                               & idprocs(iproc),idloc,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (nosnd.GT.0) nreq = nreq + 1
            CALL comms_isend_real(rdatsnd(1:nosnd,:,:),nosnd*numdim,&
                                & idprocs(iproc),idloc,commop,imode,&
                                & irequest(nreq),ivarid)
         ENDIF

!
!5.2 Receive
!-----------
!

      ELSEIF (comprocs(n).EQ.2) THEN
!        ---receive
         norcv = nostatprocs(iproc)
         IF (iform.EQ.1) THEN 
            CALL comms_recv_real(rdatrcv(1:norcv,:,:),norcv*numdim,&
                               & idprocs(iproc),idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (norcv.GT.0) nreq = nreq + 1
            CALL comms_irecv_real(rdatrcv(1:norcv,:,:),norcv*numdim,&
                                & idprocs(iproc),idprocs(iproc),commop,&
                                & irequest(nreq),ivarid)
         ENDIF
!           ---store
         lloc_521: DO lloc=1,norcv
            l = lstatprocs(lloc,iproc)
            realglb(l,:,:) = rdatrcv(lloc,:,:)
         ENDDO lloc_521
      ENDIF

   ENDDO n_500

ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_stats_glb_real_3d

!========================================================================

SUBROUTINE combine_stats_glb_real_4d(realglb,maxstats,nostatprocs,lstatprocs,&
                                   & ivarid,idroot,comm,commtype,commall)
!************************************************************************
!
! *combine_stats_glb_real_4d* Combine global real 4-D data defined at
!                             stations on different domains
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.6
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!             - array elements with the same first index are defined on the
!               same domain on input
!
! Module calls - comms_barrier, comms_irecv_real,
!                comms_isend_real, comms_recv_real, comms_send_real,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid, maxstats
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
REAL, INTENT(INOUT), DIMENSION(:,:,:,:) :: realglb

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realglb*     REAL    Global real array
!*maxstats*    INTEGER First dimension of array lstatprocs
!*nostatprocs* INTEGER Number of stations in each local array
!*lstatprocs*  INTEGER Indices of local stations in global array
!*ivarid*      INTEGER Variable id
!*idroot*      INTEGER id of root process
!*comm*        INTEGER Communicator for combine operation
!*commtype*    INTEGER Type of communication (1/4)
!*commall*     LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, l, lloc, n, norcv, &
         & nosnd, npcc, nreq, numdim
INTEGER, DIMENSION(2*nprocs-2) :: irequest
REAL, DIMENSION(SIZE(realglb,DIM=1),SIZE(realglb,DIM=2),&
              & SIZE(realglb,DIM=3),SIZE(realglb,DIM=4)) :: rdatrcv, rdatsnd


IF (SIZE(realglb).EQ.0) RETURN  
procname(pglev+1) = 'combine_stats_glb_real_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(iopt_MPI_comm_all,iopt_MPI_comm_gath,all)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
numdim = SIZE(realglb,DIM=2)*SIZE(realglb,DIM=3)*SIZE(realglb,DIM=4)

!
!3. Store data in local array
!----------------------------
!

IF (all.OR.idloc.NE.iddest) THEN
   iproc_310: DO iproc=1,nprocs
      IF (idloc.EQ.idprocs(iproc)) THEN
         nosnd = nostatprocs(iproc)
         lloc_311: DO lloc=1,nosnd
            l = lstatprocs(lloc,iproc)
            rdatsnd(lloc,:,:,:) = realglb(l,:,:,:)
         ENDDO lloc_311
      ENDIF
   ENDDO iproc_310
ENDIF

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!4.1 Send
!--------
!

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_real(rdatsnd(1:nosnd,:,:,:),nosnd*numdim,iddest,idloc,&
                            & commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         IF (nosnd.GT.0) nreq = nreq + 1
         CALL comms_isend_real(rdatsnd(1:nosnd,:,:,:),nosnd*numdim,iddest,&
                             & idloc,commop,imode,irequest(nreq),ivarid)
      ENDIF
   ENDIF

!
!4.2 Receive
!-----------
!

   IF (idloc.EQ.iddest) THEN
      iproc_420: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.iddest) THEN
!           ---receive
            norcv = nostatprocs(iproc)
            IF (iform.EQ.1) THEN 
               CALL comms_recv_real(rdatrcv(1:norcv,:,:,:),norcv*numdim,&
                                  & idprocs(iproc),idprocs(iproc),commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               IF (norcv.GT.0) nreq = nreq + 1
               CALL comms_irecv_real(rdatrcv(1:norcv,:,:,:),norcv*numdim,&
                                   & idprocs(iproc),idprocs(iproc),commop,&
                                   & irequest(nreq),ivarid)
            ENDIF
!           ---store
            lloc_421: DO lloc=1,norcv
               l = lstatprocs(lloc,iproc)
               realglb(l,:,:,:) = rdatrcv(lloc,:,:,:)
            ENDDO lloc_421
         ENDIF
      ENDDO iproc_420
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE
   n_500: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1

!
!5.1 Send
!--------
!

      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_real(rdatsnd(1:nosnd,:,:,:),nosnd*numdim,&
                               & idprocs(iproc),idloc,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (nosnd.GT.0) nreq = nreq + 1
            CALL comms_isend_real(rdatsnd(1:nosnd,:,:,:),nosnd*numdim,&
                                & idprocs(iproc),idloc,commop,imode,&
                                & irequest(nreq),ivarid)
         ENDIF

!
!5.2 Receive
!-----------
!

      ELSEIF (comprocs(n).EQ.2) THEN
!        ---receive
         norcv = nostatprocs(iproc)
         IF (iform.EQ.1) THEN 
            CALL comms_recv_real(rdatrcv(1:norcv,:,:,:),norcv*numdim,&
                               & idprocs(iproc),idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            IF (norcv.GT.0) nreq = nreq + 1
            CALL comms_irecv_real(rdatrcv(1:norcv,:,:,:),norcv*numdim,&
                                & idprocs(iproc),idprocs(iproc),commop,&
                                & irequest(nreq),ivarid)
         ENDIF
!           ---store
         lloc_521: DO lloc=1,norcv
            l = lstatprocs(lloc,iproc)
            realglb(l,:,:,:) = rdatrcv(lloc,:,:,:)
         ENDDO lloc_521
      ENDIF

   ENDDO n_500

ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_stats_glb_real_4d

!========================================================================

SUBROUTINE combine_stats_loc_real_1d(realglb,realloc,maxstats,nostatprocs,&
                                   & lstatprocs,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *combine_stats_loc_real_1d* Combine local vector real station data into a
!                             global array on the root process
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.6
!
! Description -
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall,
!                error_abort, error_alloc
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort, error_alloc

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, maxstats
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
REAL, INTENT(OUT), DIMENSION(:) :: realglb
REAL, INTENT(IN), DIMENSION(:) :: realloc

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*realglb*     REAL    Global station array
!*realloc*     REAL    Local station array with size given by number of local
!                      stations
!*maxstats*    INTEGER First dimension of array lstatprocs
!*nostatprocs* INTEGER Number of stations in each local array
!*lstatprocs*  INTEGER Indices of local stations in global array
!*ivarid*      INTEGER Variable id
!*idroot*      INTEGER id of root process
!*comm*        INTEGER Communicator for combine operation
!*commtype*    INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iddest, iform, imode, iproc, itype, l, lglb, norcv, nosnd, &
         & nostatsloc, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest
REAL, ALLOCATABLE, DIMENSION(:) :: realdat


procname(pglev+1) = 'combine_stats_loc_real_1d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   l_110: DO l=1,nostatprocs(1)
      lglb = lstatprocs(l,1)
      realglb(lglb) = realloc(l)
   ENDDO l_110
   GOTO 1000
ENDIF

!
!2. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_gath
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
nosnd = SIZE(realloc)

!
!4. Work space arrays
!--------------------
!

ALLOCATE (realdat(maxstats),STAT=errstat)
CALL error_alloc('realdat',1,(/maxstats/),kndrtype)
realdat = 0.0

!
!5. Combine
!----------
!
!5.1 Send
!--------
!

IF (idloc.NE.iddest) THEN

   IF (iform.EQ.1) THEN
      CALL comms_send_real(realloc,nosnd,iddest,idloc,commop,imode,ivarid)
   ELSEIF (iform.EQ.2) THEN
      IF (nosnd.GT.0) nreq = nreq + 1
      CALL comms_isend_real(realloc,nosnd,iddest,idloc,commop,imode,&
                     & irequest(nreq),ivarid)
   ENDIF

!
!5.2 Receive
!-----------
!

ELSE

   realglb = 0.0
   iproc_520: DO iproc=1,nprocs
      nostatsloc = nostatprocs(iproc)
      IF (idprocs(iproc).EQ.iddest) THEN
         l_521: DO l=1,nostatsloc
            lglb = lstatprocs(l,iproc)
            realglb(lglb) = realloc(l)
         ENDDO l_521
      ELSEIF (nostatsloc.GT.0) THEN
         norcv = nostatsloc
         IF (iform.EQ.1) THEN
            CALL comms_recv_real(realdat(1:nostatsloc),norcv,idprocs(iproc),&
                               & idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_real(realdat(1:nostatsloc),norcv,idprocs(iproc),&
                                & idprocs(iproc),commop,irequest(nreq),ivarid)
         ENDIF
         l_522: DO l=1,nostatsloc
            lglb = lstatprocs(l,iproc)
            realglb(lglb) = realdat(l)
         ENDDO l_522
      ENDIF
   ENDDO iproc_520
ENDIF

!
!5.3 Wait for completion
!-----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

!
!6. Deallocate
!-------------
!

DEALLOCATE (realdat)

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_stats_loc_real_1d

!========================================================================

SUBROUTINE combine_stats_loc_real_2d(realglb,realloc,maxstats,nostatprocs,&
                                   & lstatprocs,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *combine_stats_loc_real_2d* Combine local 2-D real station data into a global
!                             array on the root process
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.6
!
! Description -
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort,
!                error_alloc
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort, error_alloc

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, maxstats
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
REAL, INTENT(OUT), DIMENSION(:,:) :: realglb
REAL, INTENT(IN), DIMENSION(:,:) :: realloc

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*realglb*     REAL    Global station array
!*realloc*     REAL    Local station array
!                 first dimension: number of local stations
!                 second dimension: number of variables
!*maxstats*    INTEGER First dimension of array lstatprocs
!*nostatprocs* INTEGER Number of stations in each local array
!*lstatprocs*  INTEGER Indices of local stations in global array
!*ivarid*      INTEGER Variable id
!*idroot*      INTEGER id of root process
!*comm*        INTEGER Communicator for combine operation
!*commtype*    INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iddest, iform, imode, iproc, itype, l, lglb, norcv, nosnd, &
         & nostatsloc, novals, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest
REAL, ALLOCATABLE, DIMENSION(:,:) :: realdat


procname(pglev+1) = 'combine_stats_loc_real_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   l_110: DO l=1,nostatprocs(1)
      lglb = lstatprocs(l,1)
      realglb(lglb,:) = realloc(l,:)
   ENDDO l_110
   GOTO 1000
ENDIF

!
!2. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_gath
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
nosnd = SIZE(realloc)

!
!4. Work space arrays
!--------------------
!

novals = SIZE(realloc,DIM=2) 
ALLOCATE (realdat(maxstats,novals),STAT=errstat)
CALL error_alloc('realdat',2,(/maxstats,novals/),kndrtype)
realdat = 0.0

!
!5. Combine
!----------
!
!5.1 Send
!--------
!

IF (idloc.NE.iddest) THEN

   IF (iform.EQ.1) THEN
      CALL comms_send_real(realloc,nosnd,iddest,idloc,commop,imode,ivarid)
   ELSEIF (iform.EQ.2) THEN
      IF (nosnd.GT.0) nreq = nreq + 1
      CALL comms_isend_real(realloc,nosnd,iddest,idloc,commop,imode,&
                     & irequest(nreq),ivarid)
   ENDIF

!
!5.2 Receive
!-----------
!

ELSE

   realglb = 0.0
   iproc_520: DO iproc=1,nprocs
      nostatsloc = nostatprocs(iproc)
      IF (idprocs(iproc).EQ.iddest) THEN
         l_521: DO l=1,nostatsloc
            lglb = lstatprocs(l,iproc)
            realglb(lglb,:) = realloc(l,:)
         ENDDO l_521
      ELSEIF (nostatsloc.GT.0) THEN
         norcv = nostatsloc*novals
         IF (iform.EQ.1) THEN
            CALL comms_recv_real(realdat(1:nostatsloc,:),norcv,idprocs(iproc),&
                               & idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_real(realdat(1:nostatsloc,:),norcv,&
                                & idprocs(iproc),idprocs(iproc),commop,&
                                & irequest(nreq),ivarid)
         ENDIF
         l_522: DO l=1,nostatsloc
            lglb = lstatprocs(l,iproc)
            realglb(lglb,:) = realdat(l,:)
         ENDDO l_522
      ENDIF
   ENDDO iproc_520
ENDIF

!
!5.3 Wait for completion
!-----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

!
!6. Deallocate
!-------------
!

DEALLOCATE (realdat)

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_stats_loc_real_2d

!========================================================================

SUBROUTINE combine_stats_loc_real_3d(realglb,realloc,maxstats,nostatprocs,&
                                   & lstatprocs,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *combine_stats_loc_real_3d* Combine local 3-D real station data into a global
!                             array on the root process
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.6
!
! Description -
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort,
!                error_alloc
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort, error_alloc

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, maxstats
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
REAL, INTENT(OUT), DIMENSION(:,:,:) :: realglb
REAL, INTENT(IN), DIMENSION(:,:,:) :: realloc

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*realglb*     REAL    Global station array
!*realloc*     REAL    Local station array
!                 first dimension: number of local stations
!                 second dimension: vertical dimension
!                 third dimension: number of variables
!*maxstats*    INTEGER First dimension of array lstatprocs
!*nostatprocs* INTEGER Number of stations in each local array
!*lstatprocs*  INTEGER Indices of local stations in global array
!*ivarid*      INTEGER Variable id
!*idroot*      INTEGER id of root process
!*comm*        INTEGER Communicator for combine operation
!*commtype*    INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iddest, iform, imode, iproc, itype, l, lglb, norcv, nosnd, &
         & nostatsloc, novals, npcc, nreq, nzdim
INTEGER, DIMENSION(nprocs-1) :: irequest
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: realdat


procname(pglev+1) = 'combine_stats_loc_real_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   l_110: DO l=1,nostatprocs(1)
      lglb = lstatprocs(l,1)
      realglb(lglb,:,:) = realloc(l,:,:)
   ENDDO l_110
   GOTO 1000
ENDIF

!
!2. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_gath
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
nosnd = SIZE(realloc)

!
!4. Work space arrays
!--------------------
!

nzdim = SIZE(realloc,DIM=2)
novals = SIZE(realloc,DIM=3)
ALLOCATE (realdat(maxstats,nzdim,novals),STAT=errstat)
CALL error_alloc('realdat',3,(/maxstats,nzdim,novals/),kndrtype)
realdat = 0.0

!
!5. Combine
!----------
!
!5.1 Send
!--------
!

IF (idloc.NE.iddest) THEN

   IF (iform.EQ.1) THEN
      CALL comms_send_real(realloc,nosnd,iddest,idloc,commop,imode,ivarid)
   ELSEIF (iform.EQ.2) THEN
      IF (nosnd.GT.0) nreq = nreq + 1
      CALL comms_isend_real(realloc,nosnd,iddest,idloc,commop,imode,&
                     & irequest(nreq),ivarid)
   ENDIF

!
!5.2 Receive
!-----------
!

ELSE

   realglb = 0.0
   iproc_520: DO iproc=1,nprocs
      nostatsloc = nostatprocs(iproc)
      IF (idprocs(iproc).EQ.iddest) THEN
         l_521: DO l=1,nostatsloc
            lglb = lstatprocs(l,iproc)
            realglb(lglb,:,:) = realloc(l,:,:)
         ENDDO l_521
      ELSEIF (nostatsloc.GT.0) THEN
         norcv = nostatsloc*nzdim*novals
         IF (iform.EQ.1) THEN
            CALL comms_recv_real(realdat(1:nostatsloc,:,:),norcv,&
                               & idprocs(iproc),idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_real(realdat(1:nostatsloc,:,:),norcv,&
                                & idprocs(iproc),idprocs(iproc),commop,&
                                & irequest(nreq),ivarid)
         ENDIF
         l_522: DO l=1,nostatsloc
            lglb = lstatprocs(l,iproc)
            realglb(lglb,:,:) = realdat(l,:,:)
         ENDDO l_522
      ENDIF
   ENDDO iproc_520
ENDIF

!
!5.3 Wait for completion
!-----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

!
!6. Deallocate
!-------------
!

DEALLOCATE (realdat)

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_stats_loc_real_3d

!========================================================================

SUBROUTINE combine_submod_real_2d(realglb,realloc,limprocs,ivarid,rfill,&
                                & idroot,comm,commtype)
!************************************************************************
!
! *combine_submod_real_2d* Combine local real sub-model 2-D arrays into a
!                          global array on the root process
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), DIMENSION(2,2,nprocs) :: limprocs
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(IN) :: rfill
REAL, INTENT(OUT), DIMENSION(:,:) :: realglb
REAL, INTENT(IN), DIMENSION(:,:) :: realloc

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*realglb*   REAL    Global sub-grid array
!*realloc*   REAL    Local sub-grid array
!*limprocs*  INTEGER Start/end indices of local array section in global array
!*ivarid*    INTEGER Variable id
!*rfill*     REAL    Fill value for global array
!*idroot*    INTEGER id of root process
!*comm*      INTEGER Communicator for combine operation
!*commtype*  INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iddest, iform, imode, iproc, itype, i1, i2, j1, j2, norcv, &
         & nosnd, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


procname(pglev+1) = 'combine_submod_real_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   realglb = realloc
   GOTO 1000
ENDIF

!
!2. Initialise global array
!--------------------------
!

realglb = rfill

!
!3. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_gath
ENDIF

!
!4. Initialise parameters
!------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
nosnd = SIZE(realloc)

!
!5. Combine
!----------
!
!5.1 Send
!--------
!

IF (idloc.NE.iddest) THEN
   IF (iform.EQ.1) THEN 
      CALL comms_send_real(realloc,nosnd,iddest,idloc,commop,imode,ivarid)
   ELSEIF (iform.EQ.2) THEN
      IF (nosnd.GT.0) nreq = nreq + 1
      CALL comms_isend_real(realloc,nosnd,iddest,idloc,commop,imode,&
                          & irequest(nreq),ivarid)
   ENDIF

!
!5.2 Receive
!-----------
!

ELSE
   iproc_520: DO iproc=1,nprocs
      i1 = limprocs(1,1,iproc); i2 = limprocs(1,2,iproc)
      j1 = limprocs(2,1,iproc); j2 = limprocs(2,2,iproc)
      norcv = (i2-i1+1)*(j2-j1+1)
      IF (norcv.GT.0) THEN
         IF (idprocs(iproc).EQ.iddest) THEN
            realglb(i1:i2,j1:j2) = realloc
         ELSEIF (iform.EQ.1) THEN
            CALL comms_recv_real(realglb(i1:i2,j1:j2),norcv,idprocs(iproc),&
                               & idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_real(realglb(i1:i2,j1:j2),norcv,&
                                & idprocs(iproc),idprocs(iproc),commop,&
                                & irequest(nreq),ivarid)
         ENDIF
      ENDIF
   ENDDO iproc_520
ENDIF

!
!5.3 Wait for completion
!-----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_submod_real_2d

!========================================================================

SUBROUTINE combine_submod_real_3d(realglb,realloc,limprocs,ivarid,rfill,&
                                & idroot,comm,commtype)
!************************************************************************
!
! *combine_submod_real_3d* Combine local real sub-model 3-D arrays into a
!                          global array on the root process
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(IN) :: rfill
INTEGER, INTENT(IN), DIMENSION(2,2,nprocs) :: limprocs
REAL, INTENT(OUT), DIMENSION(:,:,:) :: realglb
REAL, INTENT(IN), DIMENSION(:,:,:) :: realloc

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*realglb*   REAL    Global sub-grid array
!*realloc*   REAL    Local sub-grid array
!*limprocs*  INTEGER Start/end indices of local array section in global array
!*ivarid*    INTEGER Variable id
!*rfill*     REAL    Fill value for global array
!*idroot*    INTEGER id of root process
!*comm*      INTEGER Communicator for combine operation
!*commtype*  INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iddest, iform, imode, iproc, itype, i1, i2, j1, j2, norcv, &
         & nosnd, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


procname(pglev+1) = 'combine_submod_real_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   realglb = realloc
   GOTO 1000
ENDIF

!
!2. Initialise global array
!--------------------------
!

realglb = rfill

!
!3. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_gath
ENDIF

!
!4. Initialise parameters
!------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
nosnd = SIZE(realloc)

!
!5. Combine
!----------
!
!5.1 Send
!--------
!

IF (idloc.NE.iddest) THEN
   IF (iform.EQ.1) THEN 
      CALL comms_send_real(realloc,nosnd,iddest,idloc,commop,imode,ivarid)
   ELSEIF (iform.EQ.2) THEN
      IF (nosnd.GT.0) nreq = nreq + 1
      CALL comms_isend_real(realloc,nosnd,iddest,idloc,commop,imode,&
                          & irequest(nreq),ivarid)
   ENDIF

!
!5.2 Receive
!-----------
!

ELSE
   iproc_520: DO iproc=1,nprocs
      i1 = limprocs(1,1,iproc); i2 = limprocs(1,2,iproc)
      j1 = limprocs(2,1,iproc); j2 = limprocs(2,2,iproc)
      norcv = (i2-i1+1)*(j2-j1+1)*SIZE(realloc,DIM=3)
      IF (norcv.GT.0) THEN
         IF (idprocs(iproc).EQ.iddest) THEN
            realglb(i1:i2,j1:j2,:) = realloc
         ELSEIF (iform.EQ.1) THEN
            CALL comms_recv_real(realglb(i1:i2,j1:j2,:),norcv,&
                               & idprocs(iproc),idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_real(realglb(i1:i2,j1:j2,:),norcv,&
                                & idprocs(iproc),idprocs(iproc),commop,&
                                & irequest(nreq),ivarid)
         ENDIF
      ENDIF
   ENDDO iproc_520
ENDIF

!
!5.3 Wait for completion
!-----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_submod_real_3d

!========================================================================

SUBROUTINE combine_submod_real_4d(realglb,realloc,limprocs,ivarid,rfill,&
                                & idroot,comm,commtype)
!************************************************************************
!
! *combine_submod_real_4d* Combine local real sub-model 4-D arrays into a
!                          global array on the root process
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(IN) :: rfill
INTEGER, INTENT(IN), DIMENSION(2,2,nprocs) :: limprocs
REAL, INTENT(OUT), DIMENSION(:,:,:,:) :: realglb
REAL, INTENT(IN), DIMENSION(:,:,:,:) :: realloc

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*realglb*   REAL    Global sub-grid array
!*realloc*   REAL    Local sub-grid array
!*limprocs*  INTEGER Start/end indices of local array section in global array
!*ivarid*    INTEGER Variable id
!*rfill*     REAL    Fill value for global array
!*idroot*    INTEGER id of root process
!*comm*      INTEGER Communicator for combine operation
!*commtype*  INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, iddest, iform, imode, iproc, itype, i1, i2, j1, j2, norcv, &
         & nosnd, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


procname(pglev+1) = 'combine_submod_real_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   realglb = realloc
   GOTO 1000
ENDIF

!
!2. Initialise global array
!--------------------------
!

realglb = rfill

!
!3. Optional arguments
!---------------------
!

IF (PRESENT(idroot)) THEN
   iddest = idroot
ELSE
   iddest = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_gath
ENDIF

!
!4. Initialise parameters and arrays
!-----------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
nosnd = SIZE(realloc)

!
!5. Combine
!----------
!
!5.1 Send
!--------
!

IF (idloc.NE.iddest) THEN
   IF (iform.EQ.1) THEN 
      CALL comms_send_real(realloc,nosnd,iddest,idloc,commop,imode,ivarid)
   ELSEIF (iform.EQ.2) THEN
      IF (nosnd.GT.0) nreq = nreq + 1
      CALL comms_isend_real(realloc,nosnd,iddest,idloc,commop,imode,&
                          & irequest(nreq),ivarid)
   ENDIF

!
!5.2 Receive
!-----------
!

ELSE
   iproc_520: DO iproc=1,nprocs
      i1 = limprocs(1,1,iproc); i2 = limprocs(1,2,iproc)
      j1 = limprocs(2,1,iproc); j2 = limprocs(2,2,iproc)
      norcv = (i2-i1+1)*(j2-j1+1)*SIZE(realloc,DIM=3)*SIZE(realloc,DIM=4)
      IF (norcv.GT.0) THEN
         IF (idprocs(iproc).EQ.iddest) THEN
            realglb(i1:i2,j1:j2,:,:) = realloc
         ELSEIF (iform.EQ.1) THEN
            CALL comms_recv_real(realglb(i1:i2,j1:j2,:,:),norcv,&
                               & idprocs(iproc),idprocs(iproc),commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_real(realglb(i1:i2,j1:j2,:,:),norcv,&
                                & idprocs(iproc),idprocs(iproc),commop,&
                                & irequest(nreq),ivarid)
         ENDIF
      ENDIF
   ENDDO iproc_520
ENDIF

!
!5.3 Wait for completion
!-----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_comb)


RETURN

END SUBROUTINE combine_submod_real_4d

!========================================================================

SUBROUTINE copy_chars_0d(chardat,lenstr,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_chars_0d* Copy a character string from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_char, comms_isend_char,
!                comms_recv_char, comms_send_char, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_char, comms_isend_char, &
                   & comms_recv_char, comms_send_char, comms_waitall 
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, lenstr
CHARACTER (LEN=lenstr), INTENT(INOUT) :: chardat
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*chardat*  CHAR    Character data to be copied
!*lenstr*   INTEGER Length of character string
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, l, npcc, nreq
CHARACTER (LEN=1), DIMENSION(lenstr) :: char1d
INTEGER, DIMENSION(nprocs-1) :: irequest


IF (lenstr.EQ.0) RETURN
procname(pglev+1) = 'copy_chars_0d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
IF (idloc.EQ.idsrce) THEN
   l_210: DO l=1,lenstr
      char1d(l) = chardat(l:l)
   ENDDO l_210
ENDIF

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_char(char1d,lenstr,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_char(char1d,lenstr,idsrce,idloc,commop,ivarid)
      l_320: DO l=1,lenstr
         chardat(l:l) = char1d(l)
      ENDDO l_320
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_char(char1d,lenstr,idprocs(iproc),idprocs(iproc),&
                                & commop,imode,irequest(nreq),ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_char(char1d,lenstr,idsrce,idloc,commop,ivarid)
      l_420: DO l=1,lenstr
         chardat(l:l) = char1d(l)
      ENDDO l_420
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_char(char1d,lenstr,idsrce,commop,ivarid)
   IF (idloc.NE.idsrce) THEN
      l_510: DO l=1,lenstr
         chardat(l:l) = char1d(l)
      ENDDO l_510
   ENDIF
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_chars_0d

!========================================================================

SUBROUTINE copy_chars_1d(chardat,lenstr,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_chars_1d* Copy a character vector array from root to local 
!                 processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_char, comms_isend_char,
!                comms_recv_char, comms_send_char, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_char, comms_isend_char, &
                   & comms_recv_char, comms_send_char, comms_waitall 
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, lenstr
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
CHARACTER (LEN=lenstr), INTENT(INOUT), DIMENSION(:) :: chardat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*chardat*  CHAR    Character data to be copied
!*lenstr*   INTEGER Length of character data
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = lenstr*SIZE(chardat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_chars_1d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_char(chardat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_char(chardat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_char(chardat,ncount,idprocs(iproc),&
                                & idprocs(iproc),commop,imode,irequest(nreq),&
                                & ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_char(chardat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_char(chardat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_chars_1d

!========================================================================

SUBROUTINE copy_chars_2d(chardat,lenstr,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_chars_2d* Copy a 2-D character array from root to local 
!                 processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_char, comms_isend_char,
!                comms_recv_char, comms_send_char, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_char, comms_isend_char, &
                   & comms_recv_char, comms_send_char, comms_waitall 
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, lenstr
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
CHARACTER (LEN=lenstr), INTENT(INOUT), DIMENSION(:,:) :: chardat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*chardat*  CHAR    Character data to be copied
!*lenstr*   INTEGER Length of character data
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = lenstr*SIZE(chardat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_chars_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_char(chardat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_char(chardat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_char(chardat,ncount,idprocs(iproc),&
                                & idprocs(iproc),commop,imode,irequest(nreq),&
                                & ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_char(chardat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_char(chardat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_chars_2d

!========================================================================

SUBROUTINE copy_chars_3d(chardat,lenstr,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_chars_3d* Copy a 3-D character array from root to local 
!                 processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_char, comms_isend_char,
!                comms_recv_char, comms_send_char, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_char, comms_isend_char, &
                   & comms_recv_char, comms_send_char, comms_waitall 
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, lenstr
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
CHARACTER (LEN=lenstr), INTENT(INOUT), DIMENSION(:,:,:) :: chardat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*chardat*  CHAR    Character data to be copied
!*lenstr*   INTEGER Length of character data
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = lenstr*SIZE(chardat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_chars_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_char(chardat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_char(chardat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_char(chardat,ncount,idprocs(iproc),&
                                & idprocs(iproc),commop,imode,irequest(nreq),&
                                & ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_char(chardat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_char(chardat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_chars_3d

!========================================================================

SUBROUTINE copy_chars_4d(chardat,lenstr,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_chars_4d* Copy a 4-D character array from root to local 
!                 processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_char, comms_isend_char,
!                comms_recv_char, comms_send_char, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_char, comms_isend_char, &
                   & comms_recv_char, comms_send_char, comms_waitall 
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid, lenstr
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
CHARACTER (LEN=lenstr), INTENT(INOUT), DIMENSION(:,:,:,:) :: chardat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*chardat*  CHAR    Character data to be copied
!*lenstr*   INTEGER Length of character data
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = lenstr*SIZE(chardat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_chars_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_char(chardat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_char(chardat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_char(chardat,ncount,idprocs(iproc),&
                                & idprocs(iproc),commop,imode,irequest(nreq),&
                                & ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_char(chardat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_char(chardat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_chars_4d

!========================================================================

SUBROUTINE copy_vars_int_0d(intdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_int_0d* Copy an integer scalar from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(INOUT) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Integer data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest
INTEGER, DIMENSION(1) :: int1d


procname(pglev+1) = 'copy_vars_int_0d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
IF (idloc.EQ.idsrce) int1d(1) = intdat

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_int(int1d,1,idprocs(iproc),idprocs(iproc),commop,&
                              & imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_int(int1d,1,idsrce,idloc,commop,ivarid)
      intdat = int1d(1)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_int(int1d,1,idprocs(iproc),idprocs(iproc),commop,&
                               & imode,irequest(nreq),ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_int(int1d,1,idsrce,idloc,commop,ivarid)
      intdat = int1d(1)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_int(int1d,1,idsrce,commop,ivarid)
   IF (idloc.NE.idsrce) intdat = int1d(1)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_int_0d

!========================================================================

SUBROUTINE copy_vars_int_1d(intdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_int_1d* Copy an integer vector array from root to local
!                    processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(INOUT), DIMENSION(:) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Integer data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = SIZE(intdat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_vars_int_1d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_int(intdat,ncount,idprocs(iproc),idprocs(iproc),&
                              & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_int(intdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_int(intdat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,irequest(nreq),ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_int(intdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_int(intdat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_int_1d

!========================================================================

SUBROUTINE copy_vars_int_2d(intdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_int_2d* Copy an integer 2-D array from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(INOUT), DIMENSION(:,:) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Integer data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = SIZE(intdat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_vars_int_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_int(intdat,ncount,idprocs(iproc),idprocs(iproc),&
                              & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_int(intdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_int(intdat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,irequest(nreq),ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_int(intdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_int(intdat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_int_2d

!========================================================================

SUBROUTINE copy_vars_int_3d(intdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_int_3d* Copy an integer 3-D array from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(INOUT), DIMENSION(:,:,:) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Integer data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = SIZE(intdat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_vars_int_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_int(intdat,ncount,idprocs(iproc),idprocs(iproc),&
                              & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_int(intdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_int(intdat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,irequest(nreq),ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_int(intdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_int(intdat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_int_3d

!========================================================================

SUBROUTINE copy_vars_int_4d(intdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_int_4d* Copy an integer 4-D array from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(INOUT), DIMENSION(:,:,:,:) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Integer data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = SIZE(intdat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_vars_int_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_int(intdat,ncount,idprocs(iproc),idprocs(iproc),&
                              & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_int(intdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_int(intdat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,irequest(nreq),ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_int(intdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_int(intdat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_int_4d

!========================================================================

SUBROUTINE copy_vars_log_0d(logdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_log_0d* Copy a logical scalar from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_log, comms_isend_log,
!                comms_recv_log, comms_send_log, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_log, comms_isend_log, &
                   & comms_recv_log, comms_send_log, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(INOUT) :: logdat
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*logdat*   LOGICAL Logical data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, npcc, nreq
LOGICAL, DIMENSION(1) :: log1d
INTEGER, DIMENSION(nprocs-1) :: irequest


procname(pglev+1) = 'copy_vars_log_0d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
IF (idloc.EQ.idsrce) log1d(1) = logdat

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_log(log1d,1,idprocs(iproc),idprocs(iproc),commop,&
                              & imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_log(log1d,1,idsrce,idloc,commop,ivarid)
      logdat = log1d(1)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_log(log1d,1,idprocs(iproc),idprocs(iproc),commop,&
                               & imode,irequest(nreq),ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_log(log1d,1,idsrce,idloc,commop,ivarid)
      logdat = log1d(1)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_log(log1d,1,idsrce,commop,ivarid)
   IF (idloc.NE.idsrce) log1d(1) = logdat
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_log_0d

!========================================================================

SUBROUTINE copy_vars_log_1d(logdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_log_1d* Copy a logical vector array from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_log, comms_isend_log,
!                comms_recv_log, comms_send_log, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_log, comms_isend_log, &
                   & comms_recv_log, comms_send_log, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
LOGICAL, INTENT(INOUT), DIMENSION(:) :: logdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*logdat*   LOGICAL Logical data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = SIZE(logdat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_vars_log_1d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_log(logdat,ncount,idprocs(iproc),idprocs(iproc),&
                              & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_log(logdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_log(logdat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,irequest(nreq),ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
       CALL comms_recv_log(logdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_log(logdat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_log_1d

!========================================================================

SUBROUTINE copy_vars_log_2d(logdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_log_2d* Copy a logical 2-D array from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_log, comms_isend_log,
!                comms_recv_log, comms_send_log, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_log, comms_isend_log, &
                   & comms_recv_log, comms_send_log, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
LOGICAL, INTENT(INOUT), DIMENSION(:,:) :: logdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*logdat*   LOGICAL Logical data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = SIZE(logdat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_vars_log_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_log(logdat,ncount,idprocs(iproc),idprocs(iproc),&
                              & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_log(logdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_log(logdat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,irequest(nreq),ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
       CALL comms_recv_log(logdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_log(logdat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_log_2d

!========================================================================

SUBROUTINE copy_vars_log_3d(logdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_log_3d* Copy a logical 3-D array from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_log, comms_isend_log,
!                comms_recv_log, comms_send_log, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_log, comms_isend_log, &
                   & comms_recv_log, comms_send_log, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
LOGICAL, INTENT(INOUT), DIMENSION(:,:,:) :: logdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*logdat*   LOGICAL Logical data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = SIZE(logdat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_vars_log_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_log(logdat,ncount,idprocs(iproc),idprocs(iproc),&
                              & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_log(logdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_log(logdat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,irequest(nreq),ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
       CALL comms_recv_log(logdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_log(logdat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_log_3d

!========================================================================

SUBROUTINE copy_vars_log_4d(logdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_log_4d* Copy a logical 4-D array from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_log, comms_isend_log,
!                comms_recv_log, comms_send_log, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_log, comms_isend_log, &
                   & comms_recv_log, comms_send_log, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
LOGICAL, INTENT(INOUT), DIMENSION(:,:,:,:) :: logdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*logdat*   LOGICAL Logical data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = SIZE(logdat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_vars_log_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_log(logdat,ncount,idprocs(iproc),idprocs(iproc),&
                              & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_log(logdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_log(logdat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,irequest(nreq),ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
       CALL comms_recv_log(logdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_log(logdat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_log_4d

!========================================================================

SUBROUTINE copy_vars_real_0d(realdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_real_0d* Copy a real scalar from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall 
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(INOUT) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Real data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest
REAL, DIMENSION(1) :: real1d


procname(pglev+1) = 'copy_vars_real_0d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0
IF (idloc.EQ.idsrce) real1d(1) = realdat

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_real(real1d,1,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_real(real1d,1,idsrce,idloc,commop,ivarid)
      realdat = real1d(1)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_real(real1d,1,idprocs(iproc),idprocs(iproc),&
                                & commop,imode,irequest(nreq),ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_real(real1d,1,idsrce,idloc,commop,ivarid)
      realdat = real1d(1)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_real(real1d,1,idsrce,commop,ivarid)
   IF (idloc.NE.idsrce) realdat = real1d(1)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_real_0d

!========================================================================

SUBROUTINE copy_vars_real_1d(realdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_real_1d* Copy a real vector array from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall 
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(INOUT), DIMENSION(:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Real data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = SIZE(realdat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_vars_real_1d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_real(realdat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_real(realdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_real(realdat,ncount,idprocs(iproc),&
                                & idprocs(iproc),commop,imode,irequest(nreq),&
                                & ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_real(realdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_real(realdat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_real_1d

!========================================================================

SUBROUTINE copy_vars_real_2d(realdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_real_2d* Copy a real 2-D array from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall 
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(INOUT), DIMENSION(:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Real data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = SIZE(realdat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_vars_real_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_real(realdat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_real(realdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_real(realdat,ncount,idprocs(iproc),&
                                & idprocs(iproc),commop,imode,irequest(nreq),&
                                & ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_real(realdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_real(realdat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_real_2d

!========================================================================

SUBROUTINE copy_vars_real_3d(realdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_real_3d* Copy a real 3-D array from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall 
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(INOUT), DIMENSION(:,:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Real data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = SIZE(realdat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_vars_real_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_real(realdat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_real(realdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_real(realdat,ncount,idprocs(iproc),&
                                & idprocs(iproc),commop,imode,irequest(nreq),&
                                & ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_real(realdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_real(realdat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_real_3d

!========================================================================

SUBROUTINE copy_vars_real_4d(realdat,ivarid,idroot,comm,commtype)
!************************************************************************
!
! *copy_vars_real_4d* Copy a real 4-D array from root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.0
!
! Description -
!
! Module calls - comms_barrier, comms_bcast_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_bcast_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall 
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(INOUT), DIMENSION(:,:,:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Real data to be copied
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for copy operation
!*commtype* INTEGER Type of communication (1/5)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, idsrce, iform, imode, iproc, itype, ncount, npcc, nreq
INTEGER, DIMENSION(nprocs-1) :: irequest


ncount = SIZE(realdat)
IF (ncount.EQ.0) RETURN
procname(pglev+1) = 'copy_vars_real_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!--------------------
!

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = MERGE(5,iopt_MPI_comm_scat,iopt_MPI_comm_coll.EQ.1)
ENDIF

!
!2. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!3. Blocking
!-----------
!

IF (iform.EQ.1) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_310: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            CALL comms_send_real(realdat,ncount,idprocs(iproc),idprocs(iproc),&
                               & commop,imode,ivarid)
         ENDIF
      ENDDO iproc_310

!  ---receive
   ELSE
      CALL comms_recv_real(realdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!
!4. Non-blocking
!---------------
!

ELSEIF (iform.EQ.2) THEN

!  ---send
   IF (idloc.EQ.idsrce) THEN
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).NE.idsrce) THEN
            nreq = nreq + 1
            CALL comms_isend_real(realdat,ncount,idprocs(iproc),&
                                & idprocs(iproc),commop,imode,irequest(nreq),&
                                & ivarid)
         ENDIF
      ENDDO iproc_410

!  ---receive
   ELSE
      CALL comms_recv_real(realdat,ncount,idsrce,idloc,commop,ivarid)
   ENDIF

!  ---wait for completion
   IF (nreq.GT.0) CALL comms_waitall(nreq,irequest(1:nreq))

!
!5. Collective
!-------------
!

ELSEIF (iform.EQ.3) THEN
   CALL comms_bcast_real(realdat,ncount,idsrce,commop,ivarid)
ENDIF

!
!6. Wait for completion
!----------------------
!

IF (iopt_MPI_sync.EQ.1) CALL comms_barrier(commop)

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_copy)


RETURN

END SUBROUTINE copy_vars_real_4d

!========================================================================

SUBROUTINE distribute_mod_int_2d(intglb,intloc,lbounds,nhdist,ivarid,ifill,&
                               & mglevel,shared,idroot,comm,commtype)
!************************************************************************
!
! *distribute_mod_int_2d* Distribute a global integer 2-D model array from
!                         root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
! Module calls - comms_barrier, comms_isend_int, comms_recv_int,
!                comms_send_int, comms_waitall, error_abort
!
!************************************************************************
!
USE grid
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_isend_int, comms_recv_int, &
                   & comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: shared
INTEGER, INTENT(IN) :: ifill, ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds
INTEGER, INTENT(IN), DIMENSION(4) :: nhdist
INTEGER, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):) :: intglb
INTEGER, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):) :: intloc

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intglb*   INTEGER Global integer grid array
!*intloc*   INTEGER Distributed local array
!*lbounds*  INTEGER Lower bounds of arrays
!*nhdist*   INTEGER Halo sizes for distribution (WESN directions)
!*ivarid*   INTEGER Variable id
!*ifill*    INTEGER Fill value for local array
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*shared*   LOGICAL Global array is known to all processes if .TRUE., to root
!                   process only if .FALSE.
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for distribute operation
!*commtype* INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: sharedx
INTEGER :: commop, idim, idlocx, idsrce, iform, imode, iproc, itype, i1, i2, &
         & jdim, j1, j2, lev, nhe, nhn, nhs, nhw, norcv, nosnd, npcc, nreq, &
         & nx, nxloc, nx1loc, nx2loc, ny, nyloc, ny1loc, ny2loc
INTEGER, DIMENSION(nprocscoh) :: nxprocs, nx1procs, nx2procs, nyprocs, &
                               & ny1procs, ny2procs
INTEGER, DIMENSION(npworld) :: idprocsx
INTEGER, DIMENSION(nprocs-1) :: irequest


procname(pglev+1) = 'distribute_mod_int_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_scat
ENDIF

IF (PRESENT(shared)) THEN
   sharedx = shared
ELSE
   sharedx = .TRUE.
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!---halo sizes
nhw = nhdist(1); nhe = nhdist(2); nhs = nhdist(3); nhn = nhdist(4)

!---fill values
intloc = ifill

!---(multi-)grid parameters
IF (lev.EQ.0) THEN
   nx = nc; ny = nr
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
   nxprocs = ncprocs; nyprocs = nrprocs
   nx1procs = nc1procs; ny1procs = nr1procs
   nx2procs = nc2procs; ny2procs = nr2procs
ELSE
   nx = mgvars(lev)%nc; ny = mgvars(lev)%nr
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
   nxprocs = mgvars(lev)%ncprocs; nyprocs = mgvars(lev)%nrprocs
   nx1procs = mgvars(lev)%nc1procs; ny1procs = mgvars(lev)%nr1procs
   nx2procs = mgvars(lev)%nc2procs; ny2procs = mgvars(lev)%nr2procs
ENDIF

!---process ids within communicator
IF (commop.EQ.comm_world_MPI) THEN
   idlocx = idlocglb
   idprocsx = idprocsglb
ELSE
   idlocx = idloc
   idprocsx(1:nprocs) = idprocs
ENDIF

!
!3. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   intloc(1-nhdist(1):nx+nhdist(2),1-nhdist(3):nx+nhdist(4)) = &
 & intglb(1-nhdist(1):nx+nhdist(2),1-nhdist(3):ny+nhdist(4))
   GOTO 1000
ENDIF

!
!4. Shared global array
!----------------------
!

IF (sharedx) THEN
   intloc(1-nhdist(1):nxloc+nhdist(2),1-nhdist(3):nyloc+nhdist(4)) = &
 & intglb(nx1loc-nhdist(1):nx2loc+nhdist(2),ny1loc-nhdist(3):ny2loc+nhdist(4))
   GOTO 1000
ENDIF

!
!5. Distribute
!-------------
!
!5.1 Send
!--------
!

IF (idlocx.EQ.idsrce) THEN
   iproc_510: DO iproc=1,nprocscoh
      i1 = nx1procs(iproc); j1 = ny1procs(iproc)
      i2 = nx2procs(iproc); j2 = ny2procs(iproc)
      idim = nxprocs(iproc); jdim = nyprocs(iproc)
      IF (idprocsx(iproc).EQ.idsrce) THEN
         intloc(1-nhw:idim+nhe,1-nhs:jdim+nhn) = &
      &  intglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn)
      ELSE
         nosnd = (idim+nhw+nhe)*(jdim+nhs+nhn)
         IF (iform.EQ.1) THEN
            CALL comms_send_int(intglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn),&
                              & nosnd,idprocsx(iproc),idprocsx(iproc),commop,&
                              & imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_int(intglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn),&
                               & nosnd,idprocsx(iproc),idprocsx(iproc),commop,&
                               & imode,irequest(nreq),ivarid)
         ENDIF
      ENDIF
   ENDDO iproc_510

!
!5.2 Receive
!-----------
!

ELSE
   norcv = (nxloc+nhw+nhe)*(nyloc+nhs+nhn)
   CALL comms_recv_int(intloc(1-nhw:ncloc+nhe,1-nhs:nrloc+nhn),norcv,&
                     & idsrce,idlocx,commop,ivarid)
ENDIF

!
!5.3 Wait for completion
!-----------------------
!

IF (iform.EQ.2.AND.nreq.GT.0) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_dist)


RETURN

END SUBROUTINE distribute_mod_int_2d

!========================================================================

SUBROUTINE distribute_mod_int_3d(intglb,intloc,lbounds,nhdist,ivarid,ifill,&
                               & mglevel,shared,idroot,comm,commtype)
!************************************************************************
!
! *distribute_mod_int_3d* Distribute a global integer 3-D model array from
!                         root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
! Module calls - comms_barrier, comms_isend_int, comms_recv_int,
!                comms_send_int, comms_waitall, error_abort
!
!************************************************************************
!
USE grid
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_isend_int, comms_recv_int, &
                   & comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: shared
INTEGER, INTENT(IN) :: ifill, ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds
INTEGER, INTENT(IN), DIMENSION(4) :: nhdist
INTEGER, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):,lbounds(3):) :: intglb
INTEGER, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):,&
                                & lbounds(3):) :: intloc
!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intglb*   INTEGER Global integer grid array
!*intloc*   INTEGER Distributed local array
!*lbounds*  INTEGER Lower bounds of arrays
!*nhdist*   INTEGER Halo sizes for distribution (WESN directions)
!*ivarid*   INTEGER Variable id
!*ifill*    INTEGER Fill value for local array
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*shared*   LOGICAL Global array is known to all processes if .TRUE., to root
!                   process only if .FALSE.
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for distribute operation
!*commtype* INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: sharedx
INTEGER :: commop, idim, idlocx, idsrce, iform, imode, iproc, itype, i1, i2, &
         & jdim, j1, j2, lev, nhe, nhn, nhs, nhw, norcv, nosnd, novals, npcc, &
         & nreq, nx, nxloc, nx1loc, nx2loc, ny, nyloc, ny1loc, ny2loc
INTEGER, DIMENSION(nprocscoh) :: nxprocs, nx1procs, nx2procs, nyprocs, &
                               & ny1procs,ny2procs
INTEGER, DIMENSION(npworld) :: idprocsx
INTEGER, DIMENSION(nprocs-1) :: irequest


procname(pglev+1) = 'distribute_mod_int_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_scat
ENDIF

IF (PRESENT(shared)) THEN
   sharedx = shared
ELSE
   sharedx = .TRUE.
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---size of last dimension
novals = SIZE(intloc,DIM=3)

!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!---halo sizes
nhw = nhdist(1); nhe = nhdist(2); nhs = nhdist(3); nhn = nhdist(4)

!---fill values
intloc = ifill

!---(multi-)grid parameters
IF (lev.EQ.0) THEN
   nx = nc; ny = nr
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
   nxprocs = ncprocs; nyprocs = nrprocs
   nx1procs = nc1procs; ny1procs = nr1procs
   nx2procs = nc2procs; ny2procs = nr2procs
ELSE
   nx = mgvars(lev)%nc; ny = mgvars(lev)%nr
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
   nxprocs = mgvars(lev)%ncprocs; nyprocs = mgvars(lev)%nrprocs
   nx1procs = mgvars(lev)%nc1procs; ny1procs = mgvars(lev)%nr1procs
   nx2procs = mgvars(lev)%nc2procs; ny2procs = mgvars(lev)%nr2procs
ENDIF

!---process ids within communicator
IF (commop.EQ.comm_world_MPI) THEN
   idlocx = idlocglb
   idprocsx = idprocsglb
ELSE
   idlocx = idloc
   idprocsx(1:nprocs) = idprocs
ENDIF

!
!3. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   intloc(1-nhdist(1):nx+nhdist(2),1-nhdist(3):ny+nhdist(4),:) = &
 & intglb(1-nhdist(1):nx+nhdist(2),1-nhdist(3):ny+nhdist(4),:)
   GOTO 1000
ENDIF

!
!4. Shared global array
!----------------------
!

IF (sharedx) THEN
   intloc(1-nhdist(1):nxloc+nhdist(2),1-nhdist(3):nyloc+nhdist(4),:) = &
 & intglb(nx1loc-nhdist(1):nx2loc+nhdist(2),ny1loc-nhdist(3):ny2loc+nhdist(4),:)
   GOTO 1000
ENDIF

!
!5. Distribute
!-------------
!
!5.1 Send
!--------
!

IF (idlocx.EQ.idsrce) THEN
   iproc_510: DO iproc=1,nprocscoh
      i1 = nx1procs(iproc); j1 = ny1procs(iproc)
      i2 = nx2procs(iproc); j2 = ny2procs(iproc)
      idim = nxprocs(iproc); jdim = nyprocs(iproc)
      IF (idprocsx(iproc).EQ.idsrce) THEN
         intloc(1-nhw:idim+nhe,1-nhs:jdim+nhn,:) = &
      &  intglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn,:)
      ELSE
         nosnd = (idim+nhw+nhe)*(jdim+nhs+nhn)*novals
         IF (iform.EQ.1) THEN
            CALL comms_send_int(intglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn,:),&
                              & nosnd,idprocsx(iproc),idprocsx(iproc),commop,&
                              & imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_int(intglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn,:),&
                               & nosnd,idprocsx(iproc),idprocsx(iproc),commop,&
                               & imode,irequest(nreq),ivarid)
         ENDIF
      ENDIF
   ENDDO iproc_510

!
!5.2 Receive
!-----------
!

ELSE
   norcv = (nxloc+nhw+nhe)*(nyloc+nhs+nhn)*novals
   CALL comms_recv_int(intloc(1-nhw:nxloc+nhe,1-nhs:nyloc+nhn,:),norcv,&
                     & idsrce,idlocx,commop,ivarid)
ENDIF

!
!5.3 Wait for completion
!-----------------------
!

IF (iform.EQ.2.AND.nreq.GT.0) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_dist)


RETURN

END SUBROUTINE distribute_mod_int_3d

!========================================================================

SUBROUTINE distribute_mod_int_4d(intglb,intloc,lbounds,nhdist,ivarid,ifill,&
                               & mglevel,shared,idroot,comm,commtype)
!************************************************************************
!
! *distribute_mod_int_4d* Distribute a global integer 4-D model array from
!                         root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
! Module calls - comms_barrier, comms_isend_int, comms_recv_int,
!                comms_send_int, comms_waitall, error_abort
!
!************************************************************************
!
USE grid
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_isend_int, comms_recv_int, &
                   & comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: shared
INTEGER, INTENT(IN) :: ifill, ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
INTEGER, INTENT(IN), DIMENSION(4) :: lbounds, nhdist
INTEGER, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):,&
                             & lbounds(3):,lbounds(4):) :: intglb
INTEGER, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):,&
                                & lbounds(3):,lbounds(4):) :: intloc
!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intglb*   INTEGER Global integer grid array
!*intloc*   INTEGER Distributed local array
!*lbounds*  INTEGER Lower bounds of arrays
!*nhdist*   INTEGER Halo sizes for distribution (WESN directions)
!*ivarid*   INTEGER Variable id
!*ifill*    INTEGER Fill value for local array
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*shared*   LOGICAL Global array is known to all processes if .TRUE., to root
!                   process only if .FALSE.
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for distribute operation
!*commtype* INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: sharedx
INTEGER :: commop, idim, idlocx, idsrce, iform, imode, iproc, itype, i1, i2, &
         & jdim, j1, j2, lev, nhe, nhn, nhs, nhw, norcv, nosnd, novals, npcc, &
         & nreq, nx, nxloc, nx1loc, nx2loc, ny, nyloc, ny1loc, ny2loc
INTEGER, DIMENSION(nprocscoh) :: nxprocs, nx1procs, nx2procs, nyprocs, &
                               & ny1procs, ny2procs
INTEGER, DIMENSION(npworld) :: idprocsx
INTEGER, DIMENSION(nprocs-1) :: irequest


procname(pglev+1) = 'distribute_mod_int_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_scat
ENDIF

IF (PRESENT(shared)) THEN
   sharedx = shared
ELSE
   sharedx = .TRUE.
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---size of last dimension
novals = SIZE(intloc,DIM=3)*SIZE(intloc,DIM=4)

!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!---halo sizes
nhw = nhdist(1); nhe = nhdist(2); nhs = nhdist(3); nhn = nhdist(4)

!---fill values
intloc = ifill

!---(multi-)grid parameters
IF (lev.EQ.0) THEN
   nx = nc; ny = nr
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
   nxprocs = ncprocs; nyprocs = nrprocs
   nx1procs = nc1procs; ny1procs = nr1procs
   nx2procs = nc2procs; ny2procs = nr2procs
ELSE
   nx = mgvars(lev)%nc; ny = mgvars(lev)%nr
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
   nxprocs = mgvars(lev)%ncprocs; nyprocs = mgvars(lev)%nrprocs
   nx1procs = mgvars(lev)%nc1procs; ny1procs = mgvars(lev)%nr1procs
   nx2procs = mgvars(lev)%nc2procs; ny2procs = mgvars(lev)%nr2procs
ENDIF

!---process ids within communicator
IF (commop.EQ.comm_world_MPI) THEN
   idlocx = idlocglb
   idprocsx = idprocsglb
ELSE
   idlocx = idloc
   idprocsx(1:nprocs) = idprocs
ENDIF

!
!3. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   intloc(1-nhdist(1):nx+nhdist(2),1-nhdist(3):ny+nhdist(4),:,:) = &
 & intglb(1-nhdist(1):nx+nhdist(2),1-nhdist(3):ny+nhdist(4),:,:)
   GOTO 1000
ENDIF

!
!4. Shared global array
!----------------------
!

IF (sharedx) THEN
   intloc(1-nhdist(1):nxloc+nhdist(2),1-nhdist(3):nyloc+nhdist(4),:,:) = &
 & intglb(nx1loc-nhdist(1):nx2loc+nhdist(2),&
        & ny1loc-nhdist(3):ny2loc+nhdist(4),:,:)
   GOTO 1000
ENDIF

!
!5. Distribute
!-------------
!
!5.1 Send
!--------
!

IF (idlocx.EQ.idsrce) THEN
   iproc_510: DO iproc=1,nprocscoh
      i1 = nx1procs(iproc); j1 = ny1procs(iproc)
      i2 = nx2procs(iproc); j2 = ny2procs(iproc)
      idim = nxprocs(iproc); jdim = nyprocs(iproc)
      IF (idprocsx(iproc).EQ.idsrce) THEN
         intloc(1-nhw:idim+nhe,1-nhs:jdim+nhn,:,:) = &
      &  intglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn,:,:)
      ELSE
         nosnd = (idim+nhw+nhe)*(jdim+nhs+nhn)*novals
         IF (iform.EQ.1) THEN
            CALL comms_send_int(intglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn,:,:),&
                              & nosnd,idprocsx(iproc),idprocsx(iproc),commop,&
                              & imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_int(intglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn,:,:),&
                               & nosnd,idprocsx(iproc),idprocsx(iproc),commop,&
                               & imode,irequest(nreq),ivarid)
         ENDIF
      ENDIF
   ENDDO iproc_510

!
!5.2 Receive
!-----------
!

ELSE
   norcv = (nxloc+nhw+nhe)*(nyloc+nhs+nhn)*novals
   CALL comms_recv_int(intloc(1-nhw:nxloc+nhe,1-nhs:nyloc+nhn,:,:),norcv,&
                     & idsrce,idlocx,commop,ivarid)
ENDIF

!
!5.3 Wait for completion
!-----------------------
!

IF (iform.EQ.2.AND.nreq.GT.0) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_dist)


RETURN

END SUBROUTINE distribute_mod_int_4d

!========================================================================

SUBROUTINE distribute_mod_log_2d(logglb,logloc,lbounds,nhdist,ivarid,lfill,&
                               & mglevel,shared,idroot,comm,commtype)
!************************************************************************
!
! *distribute_mod_log_2d* Distribute a global 2-D logical model array from
!                         root to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
! Module calls - comms_barrier, comms_isend_log, comms_recv_log,
!                comms_send_log, comms_waitall, error_abort
!
!************************************************************************
!
USE grid
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_isend_log, comms_recv_log, &
                   & comms_send_log, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: lfill
LOGICAL, INTENT(IN), OPTIONAL :: shared
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds
INTEGER, INTENT(IN), DIMENSION(4) :: nhdist
LOGICAL, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):) :: logglb
LOGICAL, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):) :: logloc
!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*logglb*   LOGICAL Global logical grid array
!*logloc*   LOGICAL Distributed local array
!*lbounds*  INTEGER Lower bounds of arrays
!*nhdist*   INTEGER Halo sizes for distribution (WESN directions)
!*ivarid*   INTEGER Variable id
!*lfill*    LOGICAL Fill value for local array
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*shared*   LOGICAL Global array is known to all processes if .TRUE., to root
!                   process only if .FALSE.
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for distribute operation
!*commtype* INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: sharedx
INTEGER :: commop, idim, idlocx, idsrce, iform, imode, iproc, itype, i1, i2, &
         & jdim, j1, j2, lev, nhe, nhn, nhs, nhw, norcv, nosnd, npcc, nreq, &
         & nx, nxloc, nx1loc, nx2loc, ny, nyloc, ny1loc, ny2loc
INTEGER, DIMENSION(nprocscoh) :: nxprocs, nx1procs, nx2procs, nyprocs, &
                               & ny1procs, ny2procs
INTEGER, DIMENSION(npworld) :: idprocsx
INTEGER, DIMENSION(nprocs-1) :: irequest


procname(pglev+1) = 'distribute_mod_log_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_scat
ENDIF

IF (PRESENT(shared)) THEN
   sharedx = shared
ELSE
   sharedx = .TRUE.
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!---halo sizes
nhw = nhdist(1); nhe = nhdist(2); nhs = nhdist(3); nhn = nhdist(4)

!---fill values
logloc = lfill

!---(multi-)grid parameters
IF (lev.EQ.0) THEN
   nx = nc; ny = nr
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
   nxprocs = ncprocs; nyprocs = nrprocs
   nx1procs = nc1procs; ny1procs = nr1procs
   nx2procs = nc2procs; ny2procs = nr2procs
ELSE
   nx = mgvars(lev)%nc; ny = mgvars(lev)%nr
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
   nxprocs = mgvars(lev)%ncprocs; nyprocs = mgvars(lev)%nrprocs
   nx1procs = mgvars(lev)%nc1procs; ny1procs = mgvars(lev)%nr1procs
   nx2procs = mgvars(lev)%nc2procs; ny2procs = mgvars(lev)%nr2procs
ENDIF

!---process ids within communicator
IF (commop.EQ.comm_world_MPI) THEN
   idlocx = idlocglb
   idprocsx = idprocsglb
ELSE
   idlocx = idloc
   idprocsx(1:nprocs) = idprocs
ENDIF

!
!3. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   logloc(1-nhdist(1):nx+nhdist(2),1-nhdist(3):nx+nhdist(4)) = &
 & logglb(1-nhdist(1):nx+nhdist(2),1-nhdist(3):ny+nhdist(4))
   GOTO 1000
ENDIF

!
!4. Shared global array
!----------------------
!

IF (sharedx) THEN
   logloc(1-nhdist(1):nxloc+nhdist(2),1-nhdist(3):nyloc+nhdist(4)) = &
 & logglb(nx1loc-nhdist(1):nx2loc+nhdist(2),ny1loc-nhdist(3):ny2loc+nhdist(4))
   GOTO 1000
ENDIF

!
!5. Distribute
!-------------
!
!5.1 Send
!--------
!

IF (idlocx.EQ.idsrce) THEN
   iproc_510: DO iproc=1,nprocscoh
      i1 = nx1procs(iproc); j1 = ny1procs(iproc)
      i2 = nx2procs(iproc); j2 = ny2procs(iproc)
      idim = nxprocs(iproc); jdim = nyprocs(iproc)
      IF (idprocsx(iproc).EQ.idsrce) THEN
         logloc(1-nhw:idim+nhe,1-nhs:jdim+nhn) = &
      &  logglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn)
      ELSE
         nosnd = (idim+nhw+nhe)*(jdim+nhs+nhn)
         IF (iform.EQ.1) THEN
            CALL comms_send_log(logglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn),&
                              & nosnd,idprocsx(iproc),idprocsx(iproc),commop,&
                              & imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_log(logglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn),&
                               & nosnd,idprocsx(iproc),idprocsx(iproc),commop,&
                               & imode,irequest(nreq),ivarid)
         ENDIF
      ENDIF
   ENDDO iproc_510

!
!5.2 Receive
!-----------
!

ELSE
   norcv = (nxloc+nhw+nhe)*(nyloc+nhs+nhn)
   CALL comms_recv_log(logloc(1-nhw:nxloc+nhe,1-nhs:nyloc+nhn),norcv,&
                     & idsrce,idlocx,commop,ivarid)
ENDIF

!
!5.3 Wait for completion
!-----------------------
!

IF (iform.EQ.2.AND.nreq.GT.0) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_dist)


RETURN

END SUBROUTINE distribute_mod_log_2d

!========================================================================

SUBROUTINE distribute_mod_real_2d(realglb,realloc,lbounds,nhdist,ivarid,rfill,&
                                & mglevel,shared,idroot,comm,commtype)
!************************************************************************
!
! *distribute_mod_real_2d* Distribute a global real 2-D model array from root
!                          to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
! Module calls - comms_barrier, comms_isend_real, comms_recv_real,
!                comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE grid
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_isend_real, comms_recv_real, &
                   & comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: shared
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
REAL, INTENT(IN) :: rfill
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds
INTEGER, INTENT(IN), DIMENSION(4) :: nhdist
REAL, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):) :: realglb
REAL, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):) :: realloc

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realglb*  REAL    Global real grid array
!*realloc*  REAL    Distributed local array
!*lbounds*  INTEGER Lower bounds of arrays
!*nhdist*   INTEGER Halo sizes for distribution (WESN directions)
!*ivarid*   INTEGER Variable id
!*rfill*    REAL    Fill value for local array
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*shared*   LOGICAL Global array is known to all processes if .TRUE., to root
!                   process only if .FALSE.
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for distribute operation
!*commtype* INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: sharedx
INTEGER :: commop, idim, idlocx, idsrce, iform, imode, iproc, itype, i1, i2, &
         & jdim, j1, j2, lev, nhe, nhn, nhs, nhw, norcv, nosnd, npcc, nreq, &
         & nx, nxloc, nx1loc, nx2loc, ny, nyloc, ny1loc, ny2loc
INTEGER, DIMENSION(nprocscoh) :: nxprocs, nx1procs, nx2procs, nyprocs, &
                               & ny1procs, ny2procs
INTEGER, DIMENSION(npworld) :: idprocsx
INTEGER, DIMENSION(nprocs-1) :: irequest


procname(pglev+1) = 'distribute_mod_real_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmaster
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_scat
ENDIF

IF (PRESENT(shared)) THEN
   sharedx = shared
ELSE
   sharedx = .TRUE.
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!---halo sizes
nhw = nhdist(1); nhe = nhdist(2); nhs = nhdist(3); nhn = nhdist(4)

!---fill values
realloc = rfill

!---(multi-)grid parameters
IF (lev.EQ.0) THEN
   nx = nc; ny = nr
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
   nxprocs = ncprocs; nyprocs = nrprocs
   nx1procs = nc1procs; ny1procs = nr1procs
   nx2procs = nc2procs; ny2procs = nr2procs
ELSE
   nx = mgvars(lev)%nc; ny = mgvars(lev)%nr
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
   nxprocs = mgvars(lev)%ncprocs; nyprocs = mgvars(lev)%nrprocs
   nx1procs = mgvars(lev)%nc1procs; ny1procs = mgvars(lev)%nr1procs
   nx2procs = mgvars(lev)%nc2procs; ny2procs = mgvars(lev)%nr2procs
ENDIF

!---process ids within communicator
IF (commop.EQ.comm_world_MPI) THEN
   idlocx = idlocglb
   idprocsx = idprocsglb
ELSE
   idlocx = idloc
   idprocsx(1:nprocs) = idprocs
ENDIF

!
!3. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   realloc(1-nhdist(1):nx+nhdist(2),1-nhdist(3):ny+nhdist(4)) = &
 & realglb(1-nhdist(1):nx+nhdist(2),1-nhdist(3):ny+nhdist(4))
   GOTO 1000
ENDIF

!
!4. Shared global array
!----------------------
!

IF (sharedx) THEN
   realloc(1-nhdist(1):nxloc+nhdist(2),1-nhdist(3):nyloc+nhdist(4)) = &
 & realglb(nx1loc-nhdist(1):nx2loc+nhdist(2),ny1loc-nhdist(3):ny2loc+nhdist(4))
   GOTO 1000
ENDIF

!
!5. Distribute
!-------------
!
!5.1 Send
!--------
!

IF (idlocx.EQ.idsrce) THEN
   iproc_510: DO iproc=1,nprocscoh
      i1 = nx1procs(iproc); j1 = ny1procs(iproc)
      i2 = nx2procs(iproc); j2 = ny2procs(iproc)
      idim = nxprocs(iproc); jdim = nyprocs(iproc)
      IF (idprocsx(iproc).EQ.idsrce) THEN
         realloc(1-nhw:idim+nhe,1-nhs:jdim+nhn) = &
      &  realglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn)
      ELSE
         nosnd = (idim+nhw+nhe)*(jdim+nhs+nhn)
         IF (iform.EQ.1) THEN
            CALL comms_send_real(realglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn),&
                               & nosnd,idprocsx(iproc),idprocsx(iproc),commop,&
                               & imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_real(realglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn),&
                                & nosnd,idprocsx(iproc),idprocsx(iproc),commop,&
                                & imode,irequest(nreq),ivarid)
         ENDIF
      ENDIF
   ENDDO iproc_510

!
!5.2 Receive
!-----------
!

ELSE
   norcv = (nxloc+nhw+nhe)*(nyloc+nhs+nhn)
   CALL comms_recv_real(realloc(1-nhw:nxloc+nhe,1-nhs:nyloc+nhn),&
                      & norcv,idsrce,idlocx,commop,ivarid)
ENDIF

!
!5.3 Wait for completion
!-----------------------
!

IF (iform.EQ.2.AND.nreq.GT.0) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_dist)


RETURN

END SUBROUTINE distribute_mod_real_2d

!========================================================================

SUBROUTINE distribute_mod_real_3d(realglb,realloc,lbounds,nhdist,ivarid,rfill,&
                                & mglevel,shared,idroot,comm,commtype)
!************************************************************************
!
! *distribute_mod_real_3d* Distribute a global real 3-D model array from root
!                          to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
! Module calls - comms_barrier, comms_isend_real, comms_recv_real,
!                comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE grid
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_isend_real, comms_recv_real, &
                   & comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: shared
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
REAL, INTENT(IN) :: rfill
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds
INTEGER, INTENT(IN), DIMENSION(4) :: nhdist
REAL, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):,lbounds(3):) :: realglb
REAL, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):,lbounds(3):) :: realloc
!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realglb*  REAL    Global real grid array
!*realloc*  REAL    Distributed local array
!*lbounds*  INTEGER Lower bounds of arrays
!*nhdist*   INTEGER Halo sizes for distribution (WESN directions)
!*ivarid*   INTEGER Variable id
!*rfill*    REAL    Fill value for local array
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*shared*   LOGICAL Global array is known to all processes if .TRUE., to root
!                   process only if .FALSE.
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for distribute operation
!*commtype* INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: sharedx
INTEGER :: commop, idim, idlocx, idsrce, iform, imode, iproc, itype, i1, i2, &
         & jdim, j1, j2, lev, nhe, nhn, nhs, nhw, norcv, nosnd, novals, npcc, &
         & nreq, nx, nxloc, nx1loc, nx2loc, ny, nyloc, ny1loc, ny2loc
INTEGER, DIMENSION(nprocscoh) :: nxprocs, nx1procs, nx2procs, nyprocs, &
                               & ny1procs, ny2procs
INTEGER, DIMENSION(npworld) :: idprocsx
INTEGER, DIMENSION(nprocs-1) :: irequest


procname(pglev+1) = 'distribute_mod_real_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_scat
ENDIF

IF (PRESENT(shared)) THEN
   sharedx = shared
ELSE
   sharedx = .TRUE.
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---size of last dimension
novals = SIZE(realloc,DIM=3)

!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!---halo sizes
nhw = nhdist(1); nhe = nhdist(2); nhs = nhdist(3); nhn = nhdist(4)

!---fill values
realloc = rfill

!---(multi-)grid parameters
IF (lev.EQ.0) THEN
   nx = nc; ny = nr
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
   nxprocs = ncprocs; nyprocs = nrprocs
   nx1procs = nc1procs; ny1procs = nr1procs
   nx2procs = nc2procs; ny2procs = nr2procs
ELSE
   nx = mgvars(lev)%nc; ny = mgvars(lev)%nr
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
   nxprocs = mgvars(lev)%ncprocs; nyprocs = mgvars(lev)%nrprocs
   nx1procs = mgvars(lev)%nc1procs; ny1procs = mgvars(lev)%nr1procs
   nx2procs = mgvars(lev)%nc2procs; ny2procs = mgvars(lev)%nr2procs
ENDIF

!---process ids within communicator
IF (commop.EQ.comm_world_MPI) THEN
   idlocx = idlocglb
   idprocsx = idprocsglb
ELSE
   idlocx = idloc
   idprocsx(1:nprocs) = idprocs
ENDIF

!
!3. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   realloc(1-nhdist(1):nx+nhdist(2),1-nhdist(3):ny+nhdist(4),:) = &
 & realglb(1-nhdist(1):nx+nhdist(2),1-nhdist(3):ny+nhdist(4),:)
   GOTO 1000
ENDIF

!
!4. Shared global array
!----------------------
!

IF (sharedx) THEN
   realloc(1-nhdist(1):nxloc+nhdist(2),1-nhdist(3):nyloc+nhdist(4),:) = &
 & realglb(nx1loc-nhdist(1):nx2loc+nhdist(2),&
         & ny1loc-nhdist(3):ny2loc+nhdist(4),:)
   GOTO 1000
ENDIF

!
!5. Distribute
!-------------
!
!5.1 Send
!--------
!

IF (idlocx.EQ.idsrce) THEN
   iproc_510: DO iproc=1,nprocscoh
      i1 = nx1procs(iproc); j1 = ny1procs(iproc)
      i2 = nx2procs(iproc); j2 = ny2procs(iproc)
      idim = nxprocs(iproc); jdim = nyprocs(iproc)
      IF (idprocsx(iproc).EQ.idsrce) THEN
         realloc(1-nhw:idim+nhe,1-nhs:jdim+nhn,:) = &
      &  realglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn,:)
      ELSE
         nosnd = (idim+nhw+nhe)*(jdim+nhs+nhn)*novals
         IF (iform.EQ.1) THEN
            CALL comms_send_real(realglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn,:),&
                               & nosnd,idprocsx(iproc),idprocsx(iproc),commop,&
                               & imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_real(realglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn,:),&
                                & nosnd,idprocsx(iproc),idprocsx(iproc),commop,&
                                & imode,irequest(nreq),ivarid)
         ENDIF
      ENDIF
   ENDDO iproc_510

!
!5.2 Receive
!-----------
!

ELSE
   norcv = (nxloc+nhw+nhe)*(nyloc+nhs+nhn)*novals
   CALL comms_recv_real(realloc(1-nhw:nxloc+nhe,1-nhs:nyloc+nhn,:),&
                      & norcv,idsrce,idlocx,commop,ivarid)
ENDIF

!
!5.3 Wait for completion
!-----------------------
!

IF (iform.EQ.2.AND.nreq.GT.0) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)

1000 CALL log_timer_out(npcc,itm_com_dist)


RETURN

END SUBROUTINE distribute_mod_real_3d

!========================================================================

SUBROUTINE distribute_mod_real_4d(realglb,realloc,lbounds,nhdist,ivarid,rfill,&
                                & mglevel,shared,idroot,comm,commtype)
!************************************************************************
!
! *distribute_mod_real_4d* Distribute a global real 4-D model array from root
!                          to local processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
! Module calls - comms_barrier, comms_isend_real, comms_recv_real,
!                comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE grid
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_isend_real, comms_recv_real, &
                   & comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: shared
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
REAL, INTENT(IN) :: rfill
INTEGER, INTENT(IN), DIMENSION(4) :: lbounds, nhdist
REAL, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):,&
                          & lbounds(3):,lbounds(4):) :: realglb
REAL, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):,&
                             & lbounds(3):,lbounds(4):) :: realloc
!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realglb*  REAL    Global real grid array
!*realloc*  REAL    Distributed local array
!*lbounds*  INTEGER Lower bounds of arrays
!*nhdist*   INTEGER Halo sizes for distribution (WESN directions)
!*ivarid*   INTEGER Variable id
!*rfill*    REAL    Fill value for local array
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*shared*   LOGICAL Global array is known to all processes if .TRUE., to root
!                   process only if .FALSE.
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for distribute operation
!*commtype* INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: sharedx
INTEGER :: commop, idim, idlocx, idsrce, iform, imode, iproc, itype, i1, i2, &
         & jdim, j1, j2, lev, nhe, nhn, nhs, nhw, norcv, nosnd, novals, npcc, &
         & nreq, nx, nxloc, nx1loc, nx2loc, ny, nyloc, ny1loc, ny2loc
INTEGER, DIMENSION(nprocscoh) :: nxprocs, nx1procs, nx2procs, nyprocs, &
                               & ny1procs, ny2procs
INTEGER, DIMENSION(npworld) :: idprocsx
INTEGER, DIMENSION(nprocs-1) :: irequest


procname(pglev+1) = 'distribute_mod_real_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_scat
ENDIF

IF (PRESENT(shared)) THEN
   sharedx = shared
ELSE
   sharedx = .TRUE.
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---size of last dimension
novals = SIZE(realloc,DIM=3)*SIZE(realloc,DIM=4)

!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!---halo sizes
nhw = nhdist(1); nhe = nhdist(2); nhs = nhdist(3); nhn = nhdist(4)

!---fill values
realloc = rfill

!---(multi-)grid parameters
IF (lev.EQ.0) THEN
   nx = nc; ny = nr
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
   nxprocs = ncprocs; nyprocs = nrprocs
   nx1procs = nc1procs; ny1procs = nr1procs
   nx2procs = nc2procs; ny2procs = nr2procs
ELSE
   nx = mgvars(lev)%nc; ny = mgvars(lev)%nr
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
   nxprocs = mgvars(lev)%ncprocs; nyprocs = mgvars(lev)%nrprocs
   nx1procs = mgvars(lev)%nc1procs; ny1procs = mgvars(lev)%nr1procs
   nx2procs = mgvars(lev)%nc2procs; ny2procs = mgvars(lev)%nr2procs
ENDIF

!---process ids within communicator
IF (commop.EQ.comm_world_MPI) THEN
   idlocx = idlocglb
   idprocsx = idprocsglb
ELSE
   idlocx = idloc
   idprocsx(1:nprocs) = idprocs
ENDIF

!
!3. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   realloc(1-nhdist(1):nx+nhdist(2),1-nhdist(3):ny+nhdist(4),:,:) = &
 & realglb(1-nhdist(1):nx+nhdist(2),1-nhdist(3):ny+nhdist(4),:,:)
   GOTO 1000
ENDIF

!
!4. Shared global array
!----------------------
!

IF (sharedx) THEN
   realloc(1-nhdist(1):nxloc+nhdist(2),1-nhdist(3):nyloc+nhdist(4),:,:) = &
 & realglb(nx1loc-nhdist(1):nx2loc+nhdist(2),&
         & ny1loc-nhdist(3):ny2loc+nhdist(4),:,:)
   GOTO 1000
ENDIF

!
!5. Distribute
!-------------
!
!5.1 Send
!--------
!

IF (idlocx.EQ.idsrce) THEN
   iproc_510: DO iproc=1,nprocscoh
      i1 = nx1procs(iproc); j1 = ny1procs(iproc)
      i2 = nx2procs(iproc); j2 = ny2procs(iproc)
      idim = nxprocs(iproc); jdim = nyprocs(iproc)
      IF (idprocsx(iproc).EQ.idsrce) THEN
         realloc(1-nhw:idim+nhe,1-nhs:jdim+nhn,:,:) = &
      &  realglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn,:,:)
      ELSE
         nosnd = (idim+nhw+nhe)*(jdim+nhs+nhn)*novals
         IF (iform.EQ.1) THEN
            CALL comms_send_real(realglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn,:,:),&
                               & nosnd,idprocsx(iproc),idprocsx(iproc),commop,&
                               & imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_real(realglb(i1-nhw:i2+nhe,j1-nhs:j2+nhn,:,:),&
                                & nosnd,idprocsx(iproc),idprocsx(iproc),commop,&
                                & imode,irequest(nreq),ivarid)
         ENDIF
      ENDIF
   ENDDO iproc_510

!
!5.2 Receive
!-----------
!

ELSE
   norcv = (nxloc+nhw+nhe)*(nyloc+nhs+nhn)*novals
   CALL comms_recv_real(realloc(1-nhw:nxloc+nhe,1-nhs:nyloc+nhn,:,:),&
                      & norcv,idsrce,idlocx,commop,ivarid)
ENDIF

!
!5.3 Wait for completion
!-----------------------
!

IF (iform.EQ.2.AND.nreq.GT.0) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
1000 CALL log_timer_out(npcc,itm_com_dist)


RETURN

END SUBROUTINE distribute_mod_real_4d

!========================================================================

SUBROUTINE distribute_mod_hgrid_2d(surfglb,surfloc,lbounds,nhdist,ivarid,ifill,&
                                 & rfill,mglevel,shared,idroot,comm,commtype)
!************************************************************************
!
! *distribute_mod_hgrid_2d* Distribute a global derived type 2-D model array
!                           of type 'HRelativeCoords' from root to local
!                           processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.9
!
! Description -
!
! Module calls - distribute_mod_int_2d, distribute_mod_real_4d, error_alloc
!
!************************************************************************
!
USE datatypes
USE multigrid
USE error_routines, ONLY: error_alloc
!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: shared
INTEGER, INTENT(IN) :: ifill, ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
REAL, INTENT(IN) :: rfill
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds
INTEGER, INTENT(IN), DIMENSION(4) :: nhdist
TYPE (HRelativeCoords), INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):) :: &
                                                                      & surfglb
TYPE (HRelativeCoords), INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):) :: &
                                                                      & surfloc

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*surfglb*  DERIVED Global real grid array
!*surfloc*  DERIVED Distributed local array
!*lbounds*  INTEGER Lower bounds of arrays
!*nhdist*   INTEGER Halo sizes for distribution (WESN directions)
!*ivarid*   INTEGER Variable id
!*ifill*    INTEGER Fill value for integer components
!*rfill*    REAL    Fill value for real components
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*shared*   LOGICAL Global array is known to all processes if .TRUE., to root
!                   process only if .FALSE.
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for distribute operation
!*commtype* INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: sharedx
INTEGER :: commop, idsrce, itype, lev, nx, nxloc, nx1loc, nx2loc, ny, nyloc, &
         & ny1loc, ny2loc, n1, n2
INTEGER, DIMENSION(2) :: uboundsglb, uboundsloc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: weightsglb, weightsloc


procname(pglev+1) = 'distribute_mod_hgrid_2d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (PRESENT(idroot)) THEN
   idsrce = idroot
ELSE
   idsrce = idmastermod
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_scat
ENDIF

IF (PRESENT(shared)) THEN
   sharedx = shared
ELSE
   sharedx = .TRUE.
ENDIF

!
!2. Initialise parameters multi-grid parameters
!----------------------------------------------
!

IF (lev.EQ.0) THEN
   nx = nc; ny = nr
   nxloc = ncloc; nyloc = nrloc
   nx1loc = nc1loc; nx2loc = nc2loc
   ny1loc = nr1loc; ny2loc = nr2loc
ELSE
   nx = mgvars(lev)%nc; ny = mgvars(lev)%nr
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
   nx1loc = mgvars(lev)%nc1loc; nx2loc = mgvars(lev)%nc2loc
   ny1loc = mgvars(lev)%nr1loc; ny2loc = mgvars(lev)%nr2loc
ENDIF

!
!3. Serial processing
!--------------------
!

IF (iopt_MPI.EQ.0) THEN
   surfloc(1-nhdist(1):nx+nhdist(2),1-nhdist(3):ny+nhdist(4)) = &
 & surfglb(1-nhdist(1):nx+nhdist(2),1-nhdist(3):ny+nhdist(4))
   GOTO 1000
ENDIF

!
!4. Shared global array
!----------------------
!

IF (sharedx) THEN
   surfloc(1-nhdist(1):nxloc+nhdist(2),1-nhdist(3):nyloc+nhdist(4)) = &
 & surfglb(nx1loc-nhdist(1):nx2loc+nhdist(2),ny1loc-nhdist(3):ny2loc+nhdist(4))
   GOTO 1000
ENDIF

!
!5. Allocate and storage
!-----------------------
!

uboundsglb = UBOUND(surfglb)
ALLOCATE(weightsglb(lbounds(1):uboundsglb(1),lbounds(2):uboundsglb(2),2,2),&
       & STAT=errstat)
CALL error_alloc('weightsglb',4,(/uboundsglb(1)-lbounds(1)+1,&
                                & uboundsglb(2)-lbounds(2)+1,2,2/),kndrtype)

uboundsloc = UBOUND(surfloc)
ALLOCATE(weightsloc(lbounds(1):uboundsloc(1),lbounds(2):uboundsloc(2),2,2),&
       & STAT=errstat)
CALL error_alloc('weightsloc',4,(/uboundsglb(1)-lbounds(1)+1,&
                                & uboundsglb(2)-lbounds(2)+1,2,2/),kndrtype)
n2_510: DO n2=lbounds(2),uboundsglb(2) 
n1_510: DO n1=lbounds(1),uboundsglb(1) 
   weightsglb(n1,n2,:,:) = surfglb(n1,n2)%weights
ENDDO n1_510
ENDDO n2_510

!
!6. Distribute
!-------------
!

CALL distribute_mod_int_2d(surfglb%icoord,surfloc%icoord,lbounds,nhdist,&
                         & ivarid,ifill,lev,shared,idsrce,commop,itype)
CALL distribute_mod_int_2d(surfglb%jcoord,surfloc%jcoord,lbounds,nhdist,&
                         & ivarid,ifill,lev,shared,idsrce,commop,itype)
CALL distribute_mod_real_4d(weightsglb,weightsloc,&
                         & (/lbounds(1),lbounds(2),2,2/),nhdist,ivarid,rfill,&
                         & lev,shared,idsrce,commop,itype)
n2_610: DO n2=lbounds(2),uboundsloc(2) 
n1_610: DO n1=lbounds(1),uboundsloc(1) 
   surfloc(n1,n2)%weights = weightsloc(n1,n2,:,:)
ENDDO n1_610
ENDDO n2_610

!
!7. Deallocate
!-------------
!

DEALLOCATE (weightsglb,weightsloc)


1000 CALL log_timer_out()


RETURN

END SUBROUTINE distribute_mod_hgrid_2d

!========================================================================

SUBROUTINE exchange_mod_int_2d(intdat,lbounds,nhexch,ivarid,comm,corners,&
                             & commtype,mglevel)
!************************************************************************
!
! *exchange_mod_int_2d* Exchange the halos of 2-D integer model arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
! Module calls - comms_barrier, comms_irecv_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_sendrecv_int,
!                comms_waitall, error_abort, halo_exch_comms
!
!************************************************************************
!
USE datatypes
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_sendrecv_int, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: corners
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, mglevel
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds
INTEGER, INTENT(IN), DIMENSION(4) :: nhexch
INTEGER, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):) :: intdat

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*intdat*    INTEGER Integer array variable for halo exchange
!*lbounds*   INTEGER Lower bounds of array
!*nhexch*    INTEGER Halo sizes in exchange communications (WESN directions)
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for halo exchange
!*corners*   LOGICAL .TRUE. if corner communications are enabled
!*commtype*  INTEGER Type of communication (1/5)
!*mglevel*   INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: commcorners, sfirst
INTEGER :: commop, iddest, idsrce, iform, imode, itype, i1rcv, i1snd, i2rcv, &
         & i2snd, j1rcv, j1snd, j2rcv, j2snd, lev, n, norcv, nosnd, npcc, &
         & nreq, tag
INTEGER, DIMENSION(2*MaxHaloComms) :: irequest
TYPE(ExchComms), DIMENSION(MaxHaloComms) :: arrcomms


IF (ALL(nhexch.EQ.0)) RETURN
procname(pglev+1) = 'exchange_mod_int_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(corners)) THEN
   commcorners = corners
ELSE
   commcorners = .TRUE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_exch
ENDIF

!
!2. Set parameters
!-----------------
!
!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!---ranks of source and destination processes
IF (lev.EQ.0) THEN
   CALL halo_exch_comms(arrcomms,nhexch,commcorners)
ELSE
   CALL halo_exch_comms(arrcomms,nhexch,commcorners,mglevel)
ENDIF

!
!3. Exchange data
!----------------
!

n_300: DO n=1,MaxHaloComms
   sfirst = arrcomms(n)%sfirst; tag = arrcomms(n)%tag
   iddest = arrcomms(n)%iddest; idsrce = arrcomms(n)%idsrce
   i1snd = arrcomms(n)%i1snd(1); i2snd = arrcomms(n)%i2snd(1)
   j1snd = arrcomms(n)%j1snd(1); j2snd = arrcomms(n)%j2snd(1)
   i1rcv = arrcomms(n)%i1rcv(1); i2rcv = arrcomms(n)%i2rcv(1)
   j1rcv = arrcomms(n)%j1rcv(1); j2rcv = arrcomms(n)%j2rcv(1)
   nosnd = (i2snd-i1snd+1)*(j2snd-j1snd+1)
   norcv = (i2rcv-i1rcv+1)*(j2rcv-j1rcv+1)

!
!3.1 Blocking
!------------
!

   IF (iform.EQ.1) THEN
      IF (sfirst) THEN
         CALL comms_send_int(intdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,tag,&
                           & commop,imode,ivarid)
         CALL comms_recv_int(intdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,idsrce,tag,&
                           & commop,ivarid)
      ELSE
         CALL comms_recv_int(intdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,idsrce,tag,&
                           & commop,ivarid)
         CALL comms_send_int(intdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,tag,&
                           & commop,imode,ivarid)
      ENDIF

!
!3.2 Non-blocking
!----------------
!

   ELSEIF (iform.EQ.2) THEN
      IF (sfirst) THEN
         nreq = nreq + 1
         CALL comms_isend_int(intdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,&
                            & tag,commop,imode,irequest(nreq),ivarid)
         nreq = nreq + 1
         CALL comms_irecv_int(intdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,idsrce,&
                            & tag,commop,irequest(nreq),ivarid)
      ELSE
         nreq = nreq + 1
         CALL comms_irecv_int(intdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,idsrce,&
                            & tag,commop,irequest(nreq),ivarid)
         nreq = nreq + 1
         CALL comms_isend_int(intdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,&
                            & tag,commop,imode,irequest(nreq),ivarid)
      ENDIF

!
!3.3 Send-receive
!----------------
!

   ELSEIF (iform.EQ.3) THEN
      CALL comms_sendrecv_int(intdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,&
                            & tag,intdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,&
                            & idsrce,tag,commop,ivarid)
   ENDIF

ENDDO n_300

!
!4. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_exch)


RETURN

END SUBROUTINE exchange_mod_int_2d

!========================================================================

SUBROUTINE exchange_mod_int_3d(intdat,lbounds,nhexch,ivarid,comm,corners,&
                             & commtype,mglevel)
!************************************************************************
!
! *exchange_mod_int_3d* Exchange the halos of 3-D integer model arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
! Module calls - comms_barrier, comms_irecv_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_sendrecv_int,
!                comms_waitall, error_abort, halo_exch_comms
!
!************************************************************************
!
USE datatypes
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_sendrecv_int, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: corners
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, mglevel
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds
INTEGER, INTENT(IN), DIMENSION(4) :: nhexch
INTEGER, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):,&
                                & lbounds(3):) :: intdat
!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*intdat*    INTEGER Integer array variable for halo exchange
!*lbounds*   INTEGER Lower bounds of array
!*nhexch*    INTEGER Halo sizes in exchange communications (WESN directions)
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for halo exchange
!*corners*   LOGICAL .TRUE. if corner communications are enabled
!*commtype*  INTEGER Type of communication (1/5)
!*mglevel*   INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: commcorners, sfirst
INTEGER :: commop, iddest, idsrce, iform, imode, itype, i1rcv, i1snd, i2rcv, &
         & i2snd, j1rcv, j1snd, j2rcv, j2snd, lev, n, norcv, nosnd, novals, &
         & npcc, nreq, tag
INTEGER, DIMENSION(2*MaxHaloComms) :: irequest
TYPE(ExchComms), DIMENSION(MaxHaloComms) :: arrcomms


IF (ALL(nhexch.EQ.0)) RETURN
procname(pglev+1) = 'exchange_mod_int_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(corners)) THEN
   commcorners = corners
ELSE
   commcorners = .TRUE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_exch
ENDIF

!
!2. Set parameters
!-----------------
!
!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!---ranks of source and destination processes
IF (lev.EQ.0) THEN
   CALL halo_exch_comms(arrcomms,nhexch,commcorners)
ELSE
   CALL halo_exch_comms(arrcomms,nhexch,commcorners,mglevel)
ENDIF

!---size of last dimension
novals = SIZE(intdat,DIM=3)

!
!3. Exchange data
!----------------
!

n_300: DO n=1,MaxHaloComms
   sfirst = arrcomms(n)%sfirst; tag = arrcomms(n)%tag
   iddest = arrcomms(n)%iddest; idsrce = arrcomms(n)%idsrce
   i1snd = arrcomms(n)%i1snd(1); i2snd = arrcomms(n)%i2snd(1)
   j1snd = arrcomms(n)%j1snd(1); j2snd = arrcomms(n)%j2snd(1)
   i1rcv = arrcomms(n)%i1rcv(1); i2rcv = arrcomms(n)%i2rcv(1)
   j1rcv = arrcomms(n)%j1rcv(1); j2rcv = arrcomms(n)%j2rcv(1)
   nosnd = (i2snd-i1snd+1)*(j2snd-j1snd+1)*novals
   norcv = (i2rcv-i1rcv+1)*(j2rcv-j1rcv+1)*novals

!
!3.1 Blocking
!------------
!

   IF (iform.EQ.1) THEN
      IF (sfirst) THEN
         CALL comms_send_int(intdat(i1snd:i2snd,j1snd:j2snd,:),nosnd,iddest,&
                           & tag,commop,imode,ivarid)
         CALL comms_recv_int(intdat(i1rcv:i2rcv,j1rcv:j2rcv,:),norcv,idsrce,&
                           & tag,commop,ivarid)
      ELSE
         CALL comms_recv_int(intdat(i1rcv:i2rcv,j1rcv:j2rcv,:),norcv,idsrce,&
                           & tag,commop,ivarid)
         CALL comms_send_int(intdat(i1snd:i2snd,j1snd:j2snd,:),nosnd,iddest,&
                           & tag,commop,imode,ivarid)
      ENDIF

!
!3.2 Non-blocking
!----------------
!

   ELSEIF (iform.EQ.2) THEN
      IF (sfirst) THEN
         nreq = nreq + 1
         CALL comms_isend_int(intdat(i1snd:i2snd,j1snd:j2snd,:),nosnd,iddest,&
                            & tag,commop,imode,irequest(nreq),ivarid)
         nreq = nreq + 1
         CALL comms_irecv_int(intdat(i1rcv:i2rcv,j1rcv:j2rcv,:),norcv,idsrce,&
                            & tag,commop,irequest(nreq),ivarid)
      ELSE
         nreq = nreq + 1
         CALL comms_irecv_int(intdat(i1rcv:i2rcv,j1rcv:j2rcv,:),norcv,idsrce,&
                            & tag,commop,irequest(nreq),ivarid)
         nreq = nreq + 1
         CALL comms_isend_int(intdat(i1snd:i2snd,j1snd:j2snd,:),nosnd,iddest,&
                            & tag,commop,imode,irequest(nreq),ivarid)
      ENDIF

!
!3.3 Send-receive
!----------------
!

   ELSEIF (iform.EQ.3) THEN
      CALL comms_sendrecv_int(intdat(i1snd:i2snd,j1snd:j2snd,:),nosnd,iddest,&
                            & tag,intdat(i1rcv:i2rcv,j1rcv:j2rcv,:),norcv,&
                            & idsrce,tag,commop,ivarid)
   ENDIF

ENDDO n_300

!
!4. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_exch)


RETURN

END SUBROUTINE exchange_mod_int_3d

!========================================================================

SUBROUTINE exchange_mod_int_4d(intdat,lbounds,nhexch,ivarid,comm,corners,&
                             & commtype,mglevel)
!************************************************************************
!
! *exchange_mod_int_4d* Exchange the halos of 4-D integer model arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
! Module calls - comms_barrier, comms_irecv_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_sendrecv_int,
!                comms_waitall, error_abort, halo_exch_comms
!
!************************************************************************
!
USE datatypes
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_sendrecv_int, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: corners
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, mglevel
INTEGER, INTENT(IN), DIMENSION(4) :: lbounds,nhexch
INTEGER, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):,&
                                & lbounds(3):,lbounds(4):) :: intdat
!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*intdat*    INTEGER Integer array variable for halo exchange
!*lbounds*   INTEGER Lower bounds of array
!*nhexch*    INTEGER Halo sizes in exchange communications (WESN directions)
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for halo exchange
!*corners*   LOGICAL .TRUE. if corner communications are enabled
!*commtype*  INTEGER Type of communication (1/5)
!*mglevel*    INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: commcorners, sfirst
INTEGER :: commop, iddest, idsrce, iform, imode, itype, i1rcv, i1snd, i2rcv, &
         & i2snd, j1rcv, j1snd, j2rcv, j2snd, lev, n, norcv, nosnd, novals, &
         & npcc, nreq, tag
INTEGER, DIMENSION(2*MaxHaloComms) :: irequest
TYPE(ExchComms), DIMENSION(MaxHaloComms) :: arrcomms


IF (ALL(nhexch.EQ.0)) RETURN
procname(pglev+1) = 'exchange_mod_int_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(corners)) THEN
   commcorners = corners
ELSE
   commcorners = .TRUE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_exch
ENDIF

!
!2. Set parameters
!-----------------
!
!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!---ranks of source and destination processes
IF (lev.EQ.0) THEN
   CALL halo_exch_comms(arrcomms,nhexch,commcorners)
ELSE
   CALL halo_exch_comms(arrcomms,nhexch,commcorners,mglevel)
ENDIF

!---size of last dimension
novals = SIZE(intdat,DIM=3)*SIZE(intdat,DIM=4)

!
!3. Exchange data
!----------------
!

n_300: DO n=1,MaxHaloComms
   sfirst = arrcomms(n)%sfirst; tag = arrcomms(n)%tag
   iddest = arrcomms(n)%iddest; idsrce = arrcomms(n)%idsrce
   i1snd = arrcomms(n)%i1snd(1); i2snd = arrcomms(n)%i2snd(1)
   j1snd = arrcomms(n)%j1snd(1); j2snd = arrcomms(n)%j2snd(1)
   i1rcv = arrcomms(n)%i1rcv(1); i2rcv = arrcomms(n)%i2rcv(1)
   j1rcv = arrcomms(n)%j1rcv(1); j2rcv = arrcomms(n)%j2rcv(1)
   nosnd = (i2snd-i1snd+1)*(j2snd-j1snd+1)*novals
   norcv = (i2rcv-i1rcv+1)*(j2rcv-j1rcv+1)*novals

!
!3.1 Blocking
!------------
!

   IF (iform.EQ.1) THEN
      IF (sfirst) THEN
         CALL comms_send_int(intdat(i1snd:i2snd,j1snd:j2snd,:,:),nosnd,iddest,&
                           & tag,commop,imode,ivarid)
         CALL comms_recv_int(intdat(i1rcv:i2rcv,j1rcv:j2rcv,:,:),norcv,idsrce,&
                           & tag,commop,ivarid)
      ELSE
         CALL comms_recv_int(intdat(i1rcv:i2rcv,j1rcv:j2rcv,:,:),norcv,idsrce,&
                           & tag,commop,ivarid)
         CALL comms_send_int(intdat(i1snd:i2snd,j1snd:j2snd,:,:),nosnd,iddest,&
                           & tag,commop,imode,ivarid)
      ENDIF

!
!3.2 Non-blocking
!----------------
!

   ELSEIF (iform.EQ.2) THEN
      IF (sfirst) THEN
         nreq = nreq + 1
         CALL comms_isend_int(intdat(i1snd:i2snd,j1snd:j2snd,:,:),nosnd,&
                            & iddest,tag,commop,imode,irequest(nreq),ivarid)
         nreq = nreq + 1
         CALL comms_irecv_int(intdat(i1rcv:i2rcv,j1rcv:j2rcv,:,:),norcv,&
                            & idsrce,tag,commop,irequest(nreq),ivarid)
      ELSE
         nreq = nreq + 1
         CALL comms_irecv_int(intdat(i1rcv:i2rcv,j1rcv:j2rcv,:,:),norcv,&
                            & idsrce,tag,commop,irequest(nreq),ivarid)
         nreq = nreq + 1
         CALL comms_isend_int(intdat(i1snd:i2snd,j1snd:j2snd,:,:),nosnd,&
                            & iddest,tag,commop,imode,irequest(nreq),ivarid)
      ENDIF

!
!3.3 Send-receive
!----------------
!

   ELSEIF (iform.EQ.3) THEN
      CALL comms_sendrecv_int(intdat(i1snd:i2snd,j1snd:j2snd,:,:),nosnd,&
                            & iddest,tag,intdat(i1rcv:i2rcv,j1rcv:j2rcv,:,:),&
                            & norcv,idsrce,tag,commop,ivarid)
   ENDIF

ENDDO n_300

!
!4. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_exch)


RETURN

END SUBROUTINE exchange_mod_int_4d

!========================================================================

SUBROUTINE exchange_mod_log_2d(logdat,lbounds,nhexch,ivarid,comm,corners,&
                             & commtype,mglevel)
!************************************************************************
!
! *exchange_mod_log_2d* Exchange the halos of 2-D logical model arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
! Module calls - comms_barrier, comms_irecv_log, comms_isend_log,
!                comms_recv_log, comms_send_log, comms_sendrecv_log,
!                comms_waitall, error_abort, halo_exch_comms
!
!************************************************************************
!
USE datatypes
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_log, comms_isend_log, &
                   & comms_recv_log, comms_send_log, comms_sendrecv_log, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: corners
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, mglevel
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds
INTEGER, INTENT(IN), DIMENSION(4) :: nhexch
LOGICAL, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):) :: logdat

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*logdat*    LOGICAL Logical array variable for halo exchange
!*lbounds*   INTEGER Lower bounds of array
!*nhexch*    INTEGER Halo sizes in exchange communications (WESN directions)
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for halo exchange
!*corners*   LOGICAL .TRUE. if corner communications are enabled
!*commtype*  INTEGER Type of communication (1/5)
!*mglevel*   INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: commcorners, sfirst
INTEGER :: commop, iddest, idsrce, iform, imode, itype, i1rcv, i1snd, i2rcv, &
         & i2snd, j1rcv, j1snd, j2rcv, j2snd, lev, n, norcv, nosnd, npcc, &
         & nreq, tag
INTEGER, DIMENSION(2*MaxHaloComms) :: irequest
TYPE(ExchComms), DIMENSION(MaxHaloComms) :: arrcomms


IF (ALL(nhexch.EQ.0)) RETURN
procname(pglev+1) = 'exchange_mod_log_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(corners)) THEN
   commcorners = corners
ELSE
   commcorners = .TRUE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_exch
ENDIF

!
!2. Set parameters
!-----------------
!
!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!---ranks of source and destination processes
IF (lev.EQ.0) THEN
   CALL halo_exch_comms(arrcomms,nhexch,commcorners)
ELSE
   CALL halo_exch_comms(arrcomms,nhexch,commcorners,mglevel)
ENDIF

!
!3. Exchange data
!----------------
!

n_300: DO n=1,MaxHaloComms
   sfirst = arrcomms(n)%sfirst; tag = arrcomms(n)%tag
   iddest = arrcomms(n)%iddest; idsrce = arrcomms(n)%idsrce
   i1snd = arrcomms(n)%i1snd(1); i2snd = arrcomms(n)%i2snd(1)
   j1snd = arrcomms(n)%j1snd(1); j2snd = arrcomms(n)%j2snd(1)
   i1rcv = arrcomms(n)%i1rcv(1); i2rcv = arrcomms(n)%i2rcv(1)
   j1rcv = arrcomms(n)%j1rcv(1); j2rcv = arrcomms(n)%j2rcv(1)
   nosnd = (i2snd-i1snd+1)*(j2snd-j1snd+1)
   norcv = (i2rcv-i1rcv+1)*(j2rcv-j1rcv+1)

!
!3.1 Blocking
!------------
!

   IF (iform.EQ.1) THEN
      IF (sfirst) THEN
         CALL comms_send_log(logdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,tag,&
                           & commop,imode,ivarid)
         CALL comms_recv_log(logdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,idsrce,tag,&
                           & commop,ivarid)
      ELSE
         CALL comms_recv_log(logdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,idsrce,tag,&
                           & commop,ivarid)
         CALL comms_send_log(logdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,tag,&
                           & commop,imode,ivarid)
      ENDIF

!
!3.2 Non-blocking
!----------------
!

   ELSEIF (iform.EQ.2) THEN
      IF (sfirst) THEN
         nreq = nreq + 1
         CALL comms_isend_log(logdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,&
                            & tag,commop,imode,irequest(nreq),ivarid)
         nreq = nreq + 1
         CALL comms_irecv_log(logdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,idsrce,&
                            & tag,commop,irequest(nreq),ivarid)
      ELSE
         nreq = nreq + 1
         CALL comms_irecv_log(logdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,idsrce,&
                            & tag,commop,irequest(nreq),ivarid)
         nreq = nreq + 1
         CALL comms_isend_log(logdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,&
                            & tag,commop,imode,irequest(nreq),ivarid)
      ENDIF

!
!3.3 Send-receive
!----------------
!

   ELSEIF (iform.EQ.3) THEN
      CALL comms_sendrecv_log(logdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,&
                            & tag,logdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,&
                            & idsrce,tag,commop,ivarid)
   ENDIF

ENDDO n_300

!
!4. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_exch)


RETURN

END SUBROUTINE exchange_mod_log_2d

!========================================================================

SUBROUTINE exchange_mod_real_2d(realdat,lbounds,nhexch,ivarid,comm,corners,&
                              & commtype,mglevel)
!************************************************************************
!
! *exchange_mod_real_2d* Exchange the halos of 2-D real model arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_sendrecv_real,
!                comms_waitall, error_abort, halo_exch_comms
!
!************************************************************************
!
USE datatypes
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_sendrecv_real, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: corners
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, mglevel
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds
INTEGER, INTENT(IN), DIMENSION(4) :: nhexch
REAL, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):) :: realdat

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*realdat*   REAL    Real array variable for halo exchange
!*lbounds*   INTEGER Lower bounds of array
!*nhexch*    INTEGER Halo sizes in exchange communications (WESN directions)
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for halo exchange
!*corners*   LOGICAL .TRUE. if corner communications are enabled
!*commtype*  INTEGER Type of communication (1/5)
!*mglevel*   INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: commcorners, sfirst
INTEGER :: commop, iddest, idsrce, iform, imode, itype, i1rcv, i1snd, i2rcv, &
         & i2snd, j1rcv, j1snd, j2rcv, j2snd, lev, n, norcv, nosnd, npcc, &
         & nreq, tag
INTEGER, DIMENSION(2*MaxHaloComms) :: irequest
TYPE(ExchComms), DIMENSION(MaxHaloComms) :: arrcomms


IF (ALL(nhexch.EQ.0)) RETURN
procname(pglev+1) = 'exchange_mod_real_2d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(corners)) THEN
   commcorners = corners
ELSE
   commcorners = .TRUE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_exch
ENDIF

!
!2. Set parameters
!-----------------
!
!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!---ranks of source and destination processes
IF (lev.EQ.0) THEN
   CALL halo_exch_comms(arrcomms,nhexch,commcorners)
ELSE
   CALL halo_exch_comms(arrcomms,nhexch,commcorners,mglevel)
ENDIF

!
!3. Exchange data
!----------------
!

n_300: DO n=1,MaxHaloComms
   sfirst = arrcomms(n)%sfirst; tag = arrcomms(n)%tag
   iddest = arrcomms(n)%iddest; idsrce = arrcomms(n)%idsrce
   i1snd = arrcomms(n)%i1snd(1); i2snd = arrcomms(n)%i2snd(1)
   j1snd = arrcomms(n)%j1snd(1); j2snd = arrcomms(n)%j2snd(1)
   i1rcv = arrcomms(n)%i1rcv(1); i2rcv = arrcomms(n)%i2rcv(1)
   j1rcv = arrcomms(n)%j1rcv(1); j2rcv = arrcomms(n)%j2rcv(1)
   nosnd = (i2snd-i1snd+1)*(j2snd-j1snd+1)
   norcv = (i2rcv-i1rcv+1)*(j2rcv-j1rcv+1)

!
!3.1 Blocking
!------------
!

   IF (iform.EQ.1) THEN
      IF (sfirst) THEN
         CALL comms_send_real(realdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,&
                            & tag,commop,imode,ivarid)
         CALL comms_recv_real(realdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,idsrce,&
                            & tag,commop,ivarid)
      ELSE
         CALL comms_recv_real(realdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,idsrce,&
                            & tag,commop,ivarid)
         CALL comms_send_real(realdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,&
                            & tag,commop,imode,ivarid)
      ENDIF

!
!3.2 Non-blocking
!----------------
!

   ELSEIF (iform.EQ.2) THEN
      IF (sfirst) THEN
         nreq = nreq + 1
         CALL comms_isend_real(realdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,&
                             & tag,commop,imode,irequest(nreq),ivarid)
         nreq = nreq + 1
         CALL comms_irecv_real(realdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,idsrce,&
                             & tag,commop,irequest(nreq),ivarid)
      ELSE
         nreq = nreq + 1
         CALL comms_irecv_real(realdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,idsrce,&
                             & tag,commop,irequest(nreq),ivarid)
         nreq = nreq + 1
         CALL comms_isend_real(realdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,&
                             & tag,commop,imode,irequest(nreq),ivarid)
      ENDIF

!
!3.3 Send-receive
!----------------
!

   ELSEIF (iform.EQ.3) THEN
      CALL comms_sendrecv_real(realdat(i1snd:i2snd,j1snd:j2snd),nosnd,iddest,&
                             & tag,realdat(i1rcv:i2rcv,j1rcv:j2rcv),norcv,&
                             & idsrce,tag,commop,ivarid)
   ENDIF

ENDDO n_300

!
!4. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_exch)


RETURN

END SUBROUTINE exchange_mod_real_2d

!========================================================================

SUBROUTINE exchange_mod_real_3d(realdat,lbounds,nhexch,ivarid,comm,corners,&
                              & commtype,mglevel)
!************************************************************************
!
! *exchange_mod_real_3d* Exchange the halos of 3-D real model arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_sendrecv_real,
!                comms_waitall, error_abort, halo_exch_comms
!
!************************************************************************
!
USE datatypes
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_sendrecv_real, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: corners
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, mglevel
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds
INTEGER, INTENT(IN), DIMENSION(4) :: nhexch
REAL, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):,lbounds(3):) :: realdat

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*realdat*   REAL    Real array variable for halo exchange
!*lbounds*   INTEGER Lower bounds of array
!*nhexch*    INTEGER Halo sizes in exchange communications (WESN directions)
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for halo exchange
!*corners*   LOGICAL .TRUE. if corner communications are enabled
!*commtype*  INTEGER Type of communication (1/5)
!*mglevel*   INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: commcorners, sfirst
INTEGER :: commop, iddest, idsrce, iform, imode, itype, i1rcv, i1snd, i2rcv, &
         & i2snd, j1rcv, j1snd, j2rcv, j2snd, lev, n, norcv, nosnd, novals, &
         & npcc, nreq, tag
INTEGER, DIMENSION(2*MaxHaloComms) :: irequest
TYPE(ExchComms), DIMENSION(MaxHaloComms) :: arrcomms


IF (ALL(nhexch.EQ.0)) RETURN
procname(pglev+1) = 'exchange_mod_real_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(corners)) THEN
   commcorners = corners
ELSE
   commcorners = .TRUE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_exch
ENDIF

!
!2. Set parameters
!-----------------
!
!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!---ranks of source and destination processes
IF (lev.EQ.0) THEN
   CALL halo_exch_comms(arrcomms,nhexch,commcorners)
ELSE
   CALL halo_exch_comms(arrcomms,nhexch,commcorners,mglevel)
ENDIF

!---size of last dimension
novals = SIZE(realdat,DIM=3)

!
!3. Exchange data
!----------------
!

n_300: DO n=1,MaxHaloComms
   sfirst = arrcomms(n)%sfirst; tag = arrcomms(n)%tag
   iddest = arrcomms(n)%iddest; idsrce = arrcomms(n)%idsrce
   i1snd = arrcomms(n)%i1snd(1); i2snd = arrcomms(n)%i2snd(1)
   j1snd = arrcomms(n)%j1snd(1); j2snd = arrcomms(n)%j2snd(1)
   i1rcv = arrcomms(n)%i1rcv(1); i2rcv = arrcomms(n)%i2rcv(1)
   j1rcv = arrcomms(n)%j1rcv(1); j2rcv = arrcomms(n)%j2rcv(1)
   nosnd = (i2snd-i1snd+1)*(j2snd-j1snd+1)*novals
   norcv = (i2rcv-i1rcv+1)*(j2rcv-j1rcv+1)*novals

!
!3.1 Blocking
!------------
!

   IF (iform.EQ.1) THEN
      IF (sfirst) THEN
         CALL comms_send_real(realdat(i1snd:i2snd,j1snd:j2snd,:),nosnd,iddest,&
                            & tag,commop,imode,ivarid)
         CALL comms_recv_real(realdat(i1rcv:i2rcv,j1rcv:j2rcv,:),norcv,idsrce,&
                            & tag,commop,ivarid)
      ELSE
         CALL comms_recv_real(realdat(i1rcv:i2rcv,j1rcv:j2rcv,:),norcv,idsrce,&
                            & tag,commop,ivarid)
         CALL comms_send_real(realdat(i1snd:i2snd,j1snd:j2snd,:),nosnd,iddest,&
                            & tag,commop,imode,ivarid)
      ENDIF

!
!3.2 Non-blocking
!----------------
!

   ELSEIF (iform.EQ.2) THEN
      IF (sfirst) THEN
         nreq = nreq + 1
         CALL comms_isend_real(realdat(i1snd:i2snd,j1snd:j2snd,:),nosnd,&
                             & iddest,tag,commop,imode,irequest(nreq),ivarid)
         nreq = nreq + 1
         CALL comms_irecv_real(realdat(i1rcv:i2rcv,j1rcv:j2rcv,:),norcv,&
                             & idsrce,tag,commop,irequest(nreq),ivarid)
      ELSE
         nreq = nreq + 1
         CALL comms_irecv_real(realdat(i1rcv:i2rcv,j1rcv:j2rcv,:),norcv,&
                             & idsrce,tag,commop,irequest(nreq),ivarid)
         nreq = nreq + 1
         CALL comms_isend_real(realdat(i1snd:i2snd,j1snd:j2snd,:),nosnd,&
                             & iddest,tag,commop,imode,irequest(nreq),ivarid)
      ENDIF

!
!3.3 Send-receive
!----------------
!

   ELSEIF (iform.EQ.3) THEN
      CALL comms_sendrecv_real(realdat(i1snd:i2snd,j1snd:j2snd,:),nosnd,&
                             & iddest,tag,realdat(i1rcv:i2rcv,j1rcv:j2rcv,:),&
                             & norcv,idsrce,tag,commop,ivarid)
   ENDIF

ENDDO n_300

!
!4. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_exch)


RETURN

END SUBROUTINE exchange_mod_real_3d

!========================================================================

SUBROUTINE exchange_mod_real_4d(realdat,lbounds,nhexch,ivarid,comm,corners,&
                              & commtype,mglevel)
!************************************************************************
!
! *exchange_mod_real_4d* Exchange the halos of 4-D real model arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_sendrecv_real,
!                comms_waitall, error_abort, halo_exch_comms
!
!************************************************************************
!
USE datatypes
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_sendrecv_real, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: corners
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, mglevel
INTEGER, INTENT(IN), DIMENSION(4) :: lbounds, nhexch
REAL, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):,&
                             & lbounds(3):,lbounds(4):) :: realdat
!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*realdat*   REAL    Real array variable for halo exchange
!*lbounds*   INTEGER Lower bounds of array
!*nhexch*    INTEGER Halo sizes in exchange communications (WESN directions)
!*ivarid*    INTEGER Variable id
!*comm*      INTEGER Communicator for halo exchange
!*corners*   LOGICAL .TRUE. if corner communications are enabled
!*commtype*  INTEGER Type of communication (1/5)
!*mglevel*   INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: commcorners, sfirst
INTEGER :: commop, iddest, idsrce, iform, imode, itype, i1rcv, i1snd, i2rcv, &
         & i2snd, j1rcv, j1snd, j2rcv, j2snd, lev, n, norcv, nosnd, novals, &
         & npcc, nreq, tag
INTEGER, DIMENSION(2*MaxHaloComms) :: irequest
TYPE(ExchComms), DIMENSION(MaxHaloComms) :: arrcomms


IF (ALL(nhexch.EQ.0)) RETURN
procname(pglev+1) = 'exchange_mod_real_4d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(mglevel)) THEN
   lev = 0
ELSE
   lev = mglevel
ENDIF

IF (ALL(nhexch.EQ.0)) RETURN
IF (PRESENT(comm)) THEN
   commop = comm
ELSE
   commop = comm_model
ENDIF

IF (PRESENT(corners)) THEN
   commcorners = corners
ELSE
   commcorners = .TRUE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   itype = iopt_MPI_comm_exch
ENDIF

!
!2. Set parameters
!-----------------
!
!---communication
imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!---ranks of source and destination processes
IF (lev.EQ.0) THEN
   CALL halo_exch_comms(arrcomms,nhexch,commcorners)
ELSE
   CALL halo_exch_comms(arrcomms,nhexch,commcorners,mglevel)
ENDIF

!---size of last dimension
novals = SIZE(realdat,DIM=3)*SIZE(realdat,DIM=4)

!
!3. Exchange data
!----------------
!

n_300: DO n=1,MaxHaloComms
   sfirst = arrcomms(n)%sfirst; tag = arrcomms(n)%tag
   iddest = arrcomms(n)%iddest; idsrce = arrcomms(n)%idsrce
   i1snd = arrcomms(n)%i1snd(1); i2snd = arrcomms(n)%i2snd(1)
   j1snd = arrcomms(n)%j1snd(1); j2snd = arrcomms(n)%j2snd(1)
   i1rcv = arrcomms(n)%i1rcv(1); i2rcv = arrcomms(n)%i2rcv(1)
   j1rcv = arrcomms(n)%j1rcv(1); j2rcv = arrcomms(n)%j2rcv(1)
   nosnd = (i2snd-i1snd+1)*(j2snd-j1snd+1)*novals
   norcv = (i2rcv-i1rcv+1)*(j2rcv-j1rcv+1)*novals

!
!3.1 Blocking
!------------
!

   IF (iform.EQ.1) THEN
      IF (sfirst) THEN
         CALL comms_send_real(realdat(i1snd:i2snd,j1snd:j2snd,:,:),nosnd,&
                            & iddest,tag,commop,imode,ivarid)
         CALL comms_recv_real(realdat(i1rcv:i2rcv,j1rcv:j2rcv,:,:),norcv,&
                            & idsrce,tag,commop,ivarid)
      ELSE
         CALL comms_recv_real(realdat(i1rcv:i2rcv,j1rcv:j2rcv,:,:),norcv,&
                            & idsrce,tag,commop,ivarid)
         CALL comms_send_real(realdat(i1snd:i2snd,j1snd:j2snd,:,:),nosnd,&
                            & iddest,tag,commop,imode,ivarid)
      ENDIF

!
!3.2 Non-blocking
!----------------
!

   ELSEIF (iform.EQ.2) THEN
      IF (sfirst) THEN
         nreq = nreq + 1
         CALL comms_isend_real(realdat(i1snd:i2snd,j1snd:j2snd,:,:),nosnd,&
                             & iddest,tag,commop,imode,irequest(nreq),ivarid)
         nreq = nreq + 1
         CALL comms_irecv_real(realdat(i1rcv:i2rcv,j1rcv:j2rcv,:,:),norcv,&
                             & idsrce,tag,commop,irequest(nreq),ivarid)
      ELSE
         nreq = nreq + 1
         CALL comms_irecv_real(realdat(i1rcv:i2rcv,j1rcv:j2rcv,:,:),norcv,&
                             & idsrce,tag,commop,irequest(nreq),ivarid)
         nreq = nreq + 1
         CALL comms_isend_real(realdat(i1snd:i2snd,j1snd:j2snd,:,:),nosnd,&
                             & iddest,tag,commop,imode,irequest(nreq),ivarid)
      ENDIF

!
!3.3 Send-receive
!----------------
!

   ELSEIF (iform.EQ.3) THEN
      CALL comms_sendrecv_real(realdat(i1snd:i2snd,j1snd:j2snd,:,:),nosnd,&
                             & iddest,tag,&
                             & realdat(i1rcv:i2rcv,j1rcv:j2rcv,:,:),norcv,&
                             & idsrce,tag,commop,ivarid)
   ENDIF

ENDDO n_300

!
!4. Wait for completion
!----------------------
!

IF (iform.EQ.2) THEN
   CALL comms_waitall(nreq,irequest(1:nreq))
ELSEIF (iopt_MPI_sync.EQ.1) THEN
   CALL comms_barrier(commop)
ENDIF

CALL error_abort(procname(pglev),ierrno_comms)
CALL log_timer_out(npcc,itm_com_exch)


RETURN

END SUBROUTINE exchange_mod_real_4d

!========================================================================

SUBROUTINE halo_exch_comms(arrcomms,nhexch,corners,mglevel)
!************************************************************************
!
! *halo_exch_comms* Parameters for halo exchanges
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_comms.f90  V2.8
!
! Description -
!
!************************************************************************
!
USE datatypes
USE multigrid

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: corners
INTEGER, INTENT(IN), DIMENSION(4) :: nhexch
TYPE(ExchComms), INTENT(OUT), DIMENSION(MaxHaloComms) :: arrcomms
INTEGER, INTENT(IN), OPTIONAL :: mglevel
!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*arrcoms*    DERIVED Properties of halo communications
!*nhexch*     INTEGER Halo sizes in exchange communications
!*corners*    LOGICAL .TRUE. if corner communications are enabled
!*mglevel*    INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: lev, nheast, nhnort, nhsout, nhwest
TYPE(ExchComms), DIMENSION(MaxHaloComms) :: hcomms

procname(pglev+1) = 'halo_exch_comms'
CALL log_timer_in()

!
!1. Optional argument
!---------------------
!

IF (PRESENT(mglevel)) THEN
   lev = mglevel
ELSE
   lev = 0
ENDIF

!
!2. Initialise
!-------------
!

IF (lev.EQ.0) THEN
   hcomms = halocomms
ELSE
   hcomms = mgvars(lev)%halocomms
ENDIF

arrcomms(:)%sfirst = hcomms(:)%sfirst 
arrcomms(:)%tag = hcomms(:)%tag
arrcomms(:)%iddest = proc_null_MPI; arrcomms(:)%idsrce = proc_null_MPI
arrcomms(:)%i1snd(1) = 1; arrcomms(:)%i2snd(1) = 1
arrcomms(:)%j1snd(1) = 1; arrcomms(:)%j2snd(1) = 1
arrcomms(:)%i1rcv(1) = 1; arrcomms(:)%i2rcv(1) = 1
arrcomms(:)%j1rcv(1) = 1; arrcomms(:)%j2rcv(1) = 1

nhwest = nhexch(1); nheast = nhexch(2)
nhsout = nhexch(3); nhnort = nhexch(4)

!
!2. Reset where needed
!---------------------
!
!---send West, receive East
IF (nheast.GT.0) THEN
   arrcomms(1)%iddest = hcomms(1)%iddest
   arrcomms(1)%idsrce = hcomms(1)%idsrce
   arrcomms(1)%i1snd(1) = hcomms(1)%i1snd(nheast)
   arrcomms(1)%i2snd(1) = hcomms(1)%i2snd(nheast)
   arrcomms(1)%j1snd(1) = hcomms(1)%j1snd(1)
   arrcomms(1)%j2snd(1) = hcomms(1)%j2snd(1)
   arrcomms(1)%i1rcv(1) = hcomms(1)%i1rcv(nheast)
   arrcomms(1)%i2rcv(1) = hcomms(1)%i2rcv(nheast)
   arrcomms(1)%j1rcv(1) = hcomms(1)%j1rcv(1)
   arrcomms(1)%j2rcv(1) = hcomms(1)%j2rcv(1)
ENDIF

!---send East, receive West
IF (nhwest.GT.0) THEN
   arrcomms(2)%iddest = hcomms(2)%iddest
   arrcomms(2)%idsrce = hcomms(2)%idsrce
   arrcomms(2)%i1snd(1) = hcomms(2)%i1snd(nhwest)
   arrcomms(2)%i2snd(1) = hcomms(2)%i2snd(nhwest)
   arrcomms(2)%j1snd(1) = hcomms(2)%j1snd(1)
   arrcomms(2)%j2snd(1) = hcomms(2)%j2snd(1)
   arrcomms(2)%i1rcv(1) = hcomms(2)%i1rcv(nhwest)
   arrcomms(2)%i2rcv(1) = hcomms(2)%i2rcv(nhwest)
   arrcomms(2)%j1rcv(1) = hcomms(2)%j1rcv(1)
   arrcomms(2)%j2rcv(1) = hcomms(2)%j2rcv(1)
ENDIF

!---send South, receive North
IF (nhnort.GT.0) THEN
   arrcomms(3)%iddest = hcomms(3)%iddest
   arrcomms(3)%idsrce = hcomms(3)%idsrce
   arrcomms(3)%i1snd(1) = hcomms(3)%i1snd(1)
   arrcomms(3)%i2snd(1) = hcomms(3)%i2snd(1)
   arrcomms(3)%j1snd(1) = hcomms(3)%j1snd(nhnort)
   arrcomms(3)%j2snd(1) = hcomms(3)%j2snd(nhnort)
   arrcomms(3)%i1rcv(1) = hcomms(3)%i1rcv(1)
   arrcomms(3)%i2rcv(1) = hcomms(3)%i2rcv(1)
   arrcomms(3)%j1rcv(1) = hcomms(3)%j1rcv(nhnort)
   arrcomms(3)%j2rcv(1) = hcomms(3)%j2rcv(nhnort)
ENDIF

!---send North, receive South
IF (nhsout.GT.0) THEN
   arrcomms(4)%iddest = hcomms(4)%iddest
   arrcomms(4)%idsrce = hcomms(4)%idsrce
   arrcomms(4)%i1snd(1) = hcomms(4)%i1snd(1)
   arrcomms(4)%i2snd(1) = hcomms(4)%i2snd(1)
   arrcomms(4)%j1snd(1) = hcomms(4)%j1snd(nhsout)
   arrcomms(4)%j2snd(1) = hcomms(4)%j2snd(nhsout)
   arrcomms(4)%i1rcv(1) = hcomms(4)%i1rcv(1)
   arrcomms(4)%i2rcv(1) = hcomms(4)%i2rcv(1)
   arrcomms(4)%j1rcv(1) = hcomms(4)%j1rcv(nhsout)
   arrcomms(4)%j2rcv(1) = hcomms(4)%j2rcv(nhsout)
ENDIF

!---send Southwest, receive Northeast
IF (nheast*nhnort.GT.0.AND.corners) THEN
   arrcomms(5)%iddest = hcomms(5)%iddest
   arrcomms(5)%idsrce = hcomms(5)%idsrce
   arrcomms(5)%i1snd(1) = hcomms(5)%i1snd(nheast)
   arrcomms(5)%i2snd(1) = hcomms(5)%i2snd(nheast)
   arrcomms(5)%j1snd(1) = hcomms(5)%j1snd(nhnort)
   arrcomms(5)%j2snd(1) = hcomms(5)%j2snd(nhnort)
   arrcomms(5)%i1rcv(1) = hcomms(5)%i1rcv(nheast)
   arrcomms(5)%i2rcv(1) = hcomms(5)%i2rcv(nheast)
   arrcomms(5)%j1rcv(1) = hcomms(5)%j1rcv(nhnort)
   arrcomms(5)%j2rcv(1) = hcomms(5)%j2rcv(nhnort)
ENDIF

!---send Northeast, receive Southwest
IF (nhwest*nhsout.GT.0.AND.corners) THEN
   arrcomms(6)%iddest = hcomms(6)%iddest
   arrcomms(6)%idsrce = hcomms(6)%idsrce
   arrcomms(6)%i1snd(1) = hcomms(6)%i1snd(nhwest)
   arrcomms(6)%i2snd(1) = hcomms(6)%i2snd(nhwest)
   arrcomms(6)%j1snd(1) = hcomms(6)%j1snd(nhsout)
   arrcomms(6)%j2snd(1) = hcomms(6)%j2snd(nhsout)
   arrcomms(6)%i1rcv(1) = hcomms(6)%i1rcv(nhwest)
   arrcomms(6)%i2rcv(1) = hcomms(6)%i2rcv(nhwest)
   arrcomms(6)%j1rcv(1) = hcomms(6)%j1rcv(nhsout)
   arrcomms(6)%j2rcv(1) = hcomms(6)%j2rcv(nhsout)
ENDIF

!---send Northwest, receive Southeast
IF (nheast*nhsout.GT.0.AND.corners) THEN
   arrcomms(7)%iddest = hcomms(7)%iddest
   arrcomms(7)%idsrce = hcomms(7)%idsrce
   arrcomms(7)%i1snd(1) = hcomms(7)%i1snd(nheast)
   arrcomms(7)%i2snd(1) = hcomms(7)%i2snd(nheast)
   arrcomms(7)%j1snd(1) = hcomms(7)%j1snd(nhsout)
   arrcomms(7)%j2snd(1) = hcomms(7)%j2snd(nhsout)
   arrcomms(7)%i1rcv(1) = hcomms(7)%i1rcv(nheast)
   arrcomms(7)%i2rcv(1) = hcomms(7)%i2rcv(nheast)
   arrcomms(7)%j1rcv(1) = hcomms(7)%j1rcv(nhsout)
   arrcomms(7)%j2rcv(1) = hcomms(7)%j2rcv(nhsout)
ENDIF

!---send Southeast, receive Northwest
IF (nhwest*nhnort.GT.0.AND.corners) THEN
   arrcomms(8)%iddest = hcomms(8)%iddest
   arrcomms(8)%idsrce = hcomms(8)%idsrce
   arrcomms(8)%i1snd(1) = hcomms(8)%i1snd(nhwest)
   arrcomms(8)%i2snd(1) = hcomms(8)%i2snd(nhwest)
   arrcomms(8)%j1snd(1) = hcomms(8)%j1snd(nhnort)
   arrcomms(8)%j2snd(1) = hcomms(8)%j2snd(nhnort)
   arrcomms(8)%i1rcv(1) = hcomms(8)%i1rcv(nhwest)
   arrcomms(8)%i2rcv(1) = hcomms(8)%i2rcv(nhwest)
   arrcomms(8)%j1rcv(1) = hcomms(8)%j1rcv(nhnort)
   arrcomms(8)%j2rcv(1) = hcomms(8)%j2rcv(nhnort)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE halo_exch_comms


END MODULE paral_comms
