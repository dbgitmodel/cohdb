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

MODULE paral_utilities
!************************************************************************
!
! *paral_utilities* Library of parallel utility routines
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.8
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description -
!
! Generic routines - max_vars, maxloc_vars, min_vars, minloc_vars, sum_vars,
!                    sum2_vars
!
!************************************************************************
!
USE iopars
USE paralpars
USE switches
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc

INTERFACE max_vars
   MODULE PROCEDURE max_vars_int_0d, max_vars_int_2d, max_vars_int_3d, &
                  & max_vars_real_0d, max_vars_real_2d, max_vars_real_3d
END INTERFACE

INTERFACE maxloc_vars
   MODULE PROCEDURE maxloc_vars_int_2d, maxloc_vars_int_3d, &
                  & maxloc_vars_real_2d, maxloc_vars_real_3d
END INTERFACE

INTERFACE min_vars
   MODULE PROCEDURE min_vars_int_0d, min_vars_int_2d, min_vars_int_3d, &
                  & min_vars_real_0d, min_vars_real_2d, min_vars_real_3d
END INTERFACE

INTERFACE minloc_vars
   MODULE PROCEDURE minloc_vars_int_2d, minloc_vars_int_3d, &
                  & minloc_vars_real_2d, minloc_vars_real_3d
END INTERFACE

INTERFACE sum_vars
   MODULE PROCEDURE sum_vars_int_0d, sum_vars_int_2d, sum_vars_int_3d, &
                  & sum_vars_real_0d, sum_vars_real_2d, sum_vars_real_3d
END INTERFACE

INTERFACE sum2_vars
   MODULE PROCEDURE sum2_vars_real_2d, sum2_vars_real_3d, sum2_vars_real_4d
END INTERFACE

CONTAINS

!========================================================================

SUBROUTINE max_vars_int_0d(intdat,intmax,ivarid,idroot,comm,commtype,commall)
!************************************************************************
!
! *max_vars_int_0d* Compute the global maximum of an integer scalar
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.7.1
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_int, comms_isend_int,
!                comms_recv_int, comms_reduce_int, comms_send_int,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_reduce_int, comms_send_int, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: intdat, ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(OUT) :: intmax

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Local value of integer input scalar
!*intmax*   INTEGER Global maximum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, n, npcc, nreq
INTEGER, DIMENSION(2*nprocs-2) :: irequest
INTEGER, DIMENSION(1) :: int1d


procname(pglev+1) = 'max_vars_int_0d'
CALL log_timer_in(npcc,ivarid)

!
!1. Serial case
!----------------
!

IF (iopt_MPI.EQ.0) THEN
   intmax = intdat
   GOTO 1000
END IF

!
!2. Optional arguments
!--------------------
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
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!4.1 Blocking/non-blocking
!-------------------------
!

   IF (iform.EQ.1.OR.iform.EQ.2) THEN

      IF (idloc.NE.iddest) THEN
         int1d(1) = intdat
         IF (iform.EQ.1) THEN
            CALL comms_send_int(int1d,1,iddest,idloc,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_int(int1d,1,iddest,idloc,commop,imode,&
                               & irequest(nreq),ivarid)
         ENDIF
      ELSE
         intmax = int_fill
         iproc_410: DO iproc=1,nprocs
            IF (idprocs(iproc).EQ.iddest) THEN
               intmax = MAX(intdat,intmax)
            ELSE
               IF (iform.EQ.1) THEN
                  CALL comms_recv_int(int1d,1,idprocs(iproc),idprocs(iproc),&
                                    & commop,ivarid)
               ELSEIF (iform.EQ.2) THEN
                  nreq = nreq + 1
                  CALL comms_irecv_int(int1d,1,idprocs(iproc),idprocs(iproc),&
                                     & commop,irequest(nreq),ivarid)
               ENDIF
               intmax = MAX(int1d(1),intmax)
            ENDIF
         ENDDO iproc_410
      ENDIF

!
!4.2 Reduce
!----------
!
      
   ELSEIF (iform.EQ.3) THEN
      CALL comms_reduce_int(intdat,intmax,1,commop,iddest,ivarid)
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE

!
!5.1 Blocking/non-blocking
!-------------------------
!

   IF (iform.EQ.1.OR.iform.EQ.2) THEN

      intmax = int_fill
      n_510: DO n=1,2*nprocs
         iproc = MOD(n-1,nprocs)+1
         IF (comprocs(n).EQ.1) THEN
            int1d(1) = intdat
            IF (iform.EQ.1) THEN
               CALL comms_send_int(int1d,1,idprocs(iproc),idloc,commop,imode,&
                                 & ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_isend_int(int1d,1,idprocs(iproc),idloc,commop,imode,&
                                  & irequest(nreq),ivarid)
            ENDIF
         ELSEIF (comprocs(n).EQ.2) THEN
            IF (iform.EQ.1) THEN
               CALL comms_recv_int(int1d,1,idprocs(iproc),idprocs(iproc),&
                                 & commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_int(int1d,1,idprocs(iproc),idprocs(iproc),&
                                  & commop,irequest(nreq),ivarid)
            ENDIF
            intmax = MAX(int1d(1),intmax)
         ELSEIF (comprocs(n).EQ.0) THEN
            intmax = MAX(intdat,intmax)
         ENDIF
      ENDDO n_510

!
!5.2 Reduce
!----------
!

   ELSEIF (iform.EQ.3) THEN
      CALL comms_reduce_int(intdat,intmax,1,commop,iddest,ivarid)
   ENDIF

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
1000 CALL log_timer_out(npcc,itm_com_util)


RETURN

END SUBROUTINE max_vars_int_0d

!========================================================================

SUBROUTINE max_vars_int_2d(intdat,intmax,ivarid,idroot,comm,commtype,commall,&
                         & mask)
!************************************************************************
!
! *max_vars_int_2d* Compute the global maximum of an integer 2-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.7.1
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - max_vars_int_0d
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(OUT) :: intmax
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
INTEGER, INTENT(IN), DIMENSION(:,:) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Local value of integer input array
!*intmax*   INTEGER Global maximum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, imaxval, itype


procname(pglev+1) = 'max_vars_int_2d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Local maximum
!----------------
!

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      imaxval = MAXVAL(intdat,MASK=mask)
   ELSE
      imaxval = -ABS(int_fill)
   ENDIF
ELSE
   imaxval = MAXVAL(intdat)
ENDIF

!
!2. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   intmax = imaxval
   GOTO 1000
ENDIF

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

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!4. Domain maximum
!-----------------
!

CALL max_vars_int_0d(imaxval,intmax,ivarid,iddest,commop,itype,all)

1000 CALL log_timer_out()


RETURN

END SUBROUTINE max_vars_int_2d

!========================================================================

SUBROUTINE max_vars_int_3d(intdat,intmax,ivarid,idroot,comm,commtype,commall,&
                         & mask)
!************************************************************************
!
! *max_vars_int_3d* Compute the global maximum of an integer 3-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.7.1
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - max_vars_int_0d
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(OUT) :: intmax
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
INTEGER, INTENT(IN), DIMENSION(:,:,:) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Local value of integer input array
!*intmax*   INTEGER Global maximum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, imaxval, itype, k
INTEGER, DIMENSION(SIZE(intdat,DIM=3)) :: izmax


procname(pglev+1) = 'max_vars_int_3d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Local maximum
!----------------
!

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      k_110: DO k=1,SIZE(intdat,DIM=3)
         izmax(k) = MAXVAL(intdat(:,:,k),MASK=mask)
      ENDDO k_110
      imaxval = MAXVAL(izmax)
   ELSE
      imaxval = -ABS(int_fill)
   ENDIF
ELSE
   imaxval = MAXVAL(intdat)
ENDIF

!
!2. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   intmax = imaxval
   GOTO 1000
ENDIF

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

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!4. Domain maximum
!-----------------
!

CALL max_vars_int_0d(imaxval,intmax,ivarid,iddest,commop,itype,all)

1000 CALL log_timer_out()


RETURN

END SUBROUTINE max_vars_int_3d

!========================================================================

SUBROUTINE max_vars_real_0d(realdat,realmax,ivarid,idroot,comm,commtype,&
                          & commall)
!************************************************************************
!
! *max_vars_real_0d* Compute the global maximum of a real scalar
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.7.1
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_reduce_real, comms_send_real,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_reduce_real, comms_send_real, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(IN) :: realdat
REAL, INTENT(OUT) :: realmax

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local value of real input scalar
!*realmax*  REAL    Global maximum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, n, npcc, nreq
INTEGER, DIMENSION(2*nprocs-2) :: irequest
REAL, DIMENSION(1) :: real1d


procname(pglev+1) = 'max_vars_real_0d'
CALL log_timer_in(npcc,ivarid)

!
!1. Serial case
!----------------
!

IF (iopt_MPI.EQ.0) THEN
   realmax = realdat
   GOTO 1000
END IF

!
!2. Optional arguments
!--------------------
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
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!4.1 Blocking/non-blocking
!-------------------------
!

   IF (iform.EQ.1.OR.iform.EQ.2) THEN

      IF (idloc.NE.iddest) THEN
         real1d(1) = realdat
         IF (iform.EQ.1) THEN
            CALL comms_send_real(real1d,1,iddest,idloc,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_real(real1d,1,iddest,idloc,commop,imode,&
                                & irequest(nreq),ivarid)
         ENDIF
      ELSE
         realmax = -ABS(float_fill)
         iproc_410: DO iproc=1,nprocs
            IF (idprocs(iproc).EQ.iddest) THEN
               realmax = MAX(realdat,realmax)
            ELSE
               IF (iform.EQ.1) THEN
                  CALL comms_recv_real(real1d,1,idprocs(iproc),idprocs(iproc),&
                                     & commop,ivarid)
               ELSEIF (iform.EQ.2) THEN
                  nreq = nreq + 1
                  CALL comms_irecv_real(real1d,1,idprocs(iproc),&
                                      & idprocs(iproc),commop,irequest(nreq),&
                                      & ivarid)
               ENDIF
               realmax = MAX(real1d(1),realmax)
            ENDIF
         ENDDO iproc_410
      ENDIF

!
!4.2 Reduce
!----------
!
      
   ELSEIF (iform.EQ.3) THEN
      CALL comms_reduce_real(realdat,realmax,1,commop,iddest,ivarid)
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE

!
!5.1 Blocking/non-blocking
!-------------------------
!

   IF (iform.EQ.1.OR.iform.EQ.2) THEN

      realmax = -ABS(float_fill)
      n_510: DO n=1,2*nprocs
         iproc = MOD(n-1,nprocs)+1
         IF (comprocs(n).EQ.1) THEN
            real1d(1) = realdat
            IF (iform.EQ.1) THEN
               CALL comms_send_real(real1d,1,idprocs(iproc),idloc,commop,&
                                  & imode,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_isend_real(real1d,1,idprocs(iproc),idloc,commop,&
                                   & imode,irequest(nreq),ivarid)
            ENDIF
         ELSEIF (comprocs(n).EQ.2) THEN
            IF (iform.EQ.1) THEN
               CALL comms_recv_real(real1d,1,idprocs(iproc),idprocs(iproc),&
                                  & commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_real(real1d,1,idprocs(iproc),idprocs(iproc),&
                                   & commop,irequest(nreq),ivarid)
            ENDIF
            realmax = MAX(real1d(1),realmax)
         ELSEIF (comprocs(n).EQ.0) THEN
            realmax = MAX(realdat,realmax)
         ENDIF
      ENDDO n_510

!
!5.2 Reduce
!----------
!

   ELSEIF (iform.EQ.3) THEN
      CALL comms_reduce_real(realdat,realmax,1,commop,iddest,ivarid)
   ENDIF

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
1000 CALL log_timer_out(npcc,itm_com_util)


RETURN

END SUBROUTINE max_vars_real_0d

!========================================================================

SUBROUTINE max_vars_real_2d(realdat,realmax,ivarid,idroot,comm,commtype,&
                          & commall,mask)
!************************************************************************
!
! *max_vars_real_2d* Compute the global maximum of a real 2-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.7.1
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - max_vars_real_0d
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(OUT) :: realmax
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
REAL, INTENT(IN), DIMENSION(:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local value of real input array
!*realmax*  REAL    Global maximum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for reduce or send/receive operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, itype
REAL :: rmaxval


procname(pglev+1) = 'max_vars_real_2d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Local maximum
!----------------
!

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      rmaxval = MAXVAL(realdat,MASK=mask)
   ELSE
      rmaxval = -ABS(float_fill)
   ENDIF
ELSE
   rmaxval = MAXVAL(realdat)
ENDIF

!
!2. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   realmax = rmaxval
   GOTO 1000
ENDIF

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

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!4. Domain maximum
!-----------------
!

CALL max_vars_real_0d(rmaxval,realmax,ivarid,iddest,commop,itype,all)

1000 CALL log_timer_out()


RETURN

END SUBROUTINE max_vars_real_2d

!========================================================================

SUBROUTINE max_vars_real_3d(realdat,realmax,ivarid,idroot,comm,commtype,&
                          & commall,mask)
!************************************************************************
!
! *max_vars_real_3d* Compute the global maximum of a real 3-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.7.1
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - max_vars_real_0d
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(OUT) :: realmax
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
REAL, INTENT(IN), DIMENSION(:,:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local value of real input array
!*realmax*  REAL    Global maximum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for reduce or send/receive operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, itype, k
REAL :: rmaxval
REAL, DIMENSION(SIZE(realdat,DIM=3)) :: rzmax


procname(pglev+1) = 'max_vars_real_3d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Local maximum
!----------------
!

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      k_110: DO k=1,SIZE(realdat,DIM=3)
         rzmax(k) = MAXVAL(realdat(:,:,k),MASK=mask)
      ENDDO k_110
      rmaxval = MAXVAL(rzmax)
   ELSE
      rmaxval = -ABS(float_fill)
   ENDIF
ELSE
   rmaxval = MAXVAL(realdat)
ENDIF

!
!2. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   realmax = rmaxval
   GOTO 1000
ENDIF

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

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!4. Domain maximum
!-----------------
!

CALL max_vars_real_0d(rmaxval,realmax,ivarid,iddest,commop,itype,all)

1000 CALL log_timer_out()


RETURN

END SUBROUTINE max_vars_real_3d

!========================================================================

SUBROUTINE maxloc_vars_int_2d(intdat,intmax,maxpos,ivarid,mglevel,idroot,comm,&
                            & commtype,commall,mask)
!************************************************************************
!
! *maxloc_vars_int_2d* Compute the value and location of the global maximum of
!                       an integer 2-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.8
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_waitall, error_abort
!
!************************************************************************
!
USE gridpars
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(OUT) :: intmax
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
INTEGER, INTENT(OUT), DIMENSION(2) :: maxpos
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
INTEGER, INTENT(IN), DIMENSION(:,:) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Local value of integer input array
!*intmax*   INTEGER Global maximum on output
!*maxpos*   INTEGER Global grid indices of maximum
!*ivarid*   INTEGER Variable id
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/8)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, n, npcc, nreq, nx1loc, &
         & ny1loc
INTEGER, DIMENSION(3) :: imaxtmp, imaxvals
INTEGER, DIMENSION(2*nprocs-2) :: irequest


procname(pglev+1) = 'maxloc_vars_int_2d'
CALL log_timer_in(npcc,ivarid)

!
!1. Value and location of local maximum
!--------------------------------------
!

IF (PRESENT(mglevel)) THEN
   nx1loc = mgvars(mglevel)%nc1loc; ny1loc = mgvars(mglevel)%nr1loc
ELSE
   nx1loc = nc1loc; ny1loc = nr1loc
ENDIF

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      imaxvals(3) = MAXVAL(intdat,MASK=mask)
      imaxvals(1:2) = MAXLOC(intdat,MASK=mask)
      imaxvals(1) = imaxvals(1) + nx1loc - 1
      imaxvals(2) = imaxvals(2) + ny1loc - 1
   ELSE
      imaxvals(3) = -ABS(int_fill)
      imaxvals(1:2) = 0
   ENDIF
ELSE
   imaxvals(3) = MAXVAL(intdat)
   imaxvals(1:2) = MAXLOC(intdat)
   imaxvals(1) = imaxvals(1) + nx1loc - 1
   imaxvals(2) = imaxvals(2) + ny1loc - 1
ENDIF

!---serial case
IF (iopt_MPI.EQ.0) THEN
   intmax = imaxvals(3)
   maxpos = imaxvals(1:2)
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
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
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_int(imaxvals,3,iddest,idloc,commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_int(imaxvals,3,iddest,idloc,commop,imode,&
                            & irequest(nreq),ivarid)
      ENDIF
   ELSE
      intmax = -ABS(int_fill)
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).EQ.iddest) THEN
            IF (imaxvals(3).GT.intmax) THEN
               intmax = imaxvals(3)
               maxpos = imaxvals(1:2)
            ENDIF
         ELSE
            IF (iform.EQ.1) THEN
               CALL comms_recv_int(imaxtmp,3,idprocs(iproc),idprocs(iproc),&
                                 & commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_int(imaxtmp,3,idprocs(iproc),idprocs(iproc),&
                                  & commop,irequest(nreq),ivarid)
            ENDIF
            IF (imaxtmp(3).GT.intmax) THEN
               intmax = imaxtmp(3)
               maxpos = imaxtmp(1:2)
            ENDIF
         ENDIF
      ENDDO iproc_410
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE

   intmax = -ABS(int_fill)
   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_int(imaxvals,3,idprocs(iproc),idloc,commop,imode,&
                              & ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_int(imaxvals,3,idprocs(iproc),idloc,commop,imode,&
                               & irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.2) THEN
         IF (iform.EQ.1) THEN
            CALL comms_recv_int(imaxtmp,3,idprocs(iproc),idprocs(iproc),&
                              & commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_int(imaxtmp,3,idprocs(iproc),idprocs(iproc),&
                               & commop,irequest(nreq),ivarid)
         ENDIF
         IF (imaxtmp(3).GT.intmax) THEN
            intmax = imaxtmp(3)
            maxpos = imaxtmp(1:2)
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         IF (imaxvals(3).GT.intmax) THEN
            intmax = imaxvals(3)
            maxpos = imaxvals(1:2)
         ENDIF
      ENDIF
   ENDDO n_510

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

1000 CALL log_timer_out(npcc,itm_com_util)


RETURN

END SUBROUTINE maxloc_vars_int_2d

!========================================================================

SUBROUTINE maxloc_vars_int_3d(intdat,intmax,maxpos,ivarid,mglevel,idroot,comm,&
                            & commtype,commall,mask)
!************************************************************************
!
! *maxloc_vars_int_3d* Compute the value and location of the global maximum of
!                      an integer 3-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.8
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_waitall, error_abort
!
!************************************************************************
!
USE gridpars
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(OUT) :: intmax
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
INTEGER, INTENT(OUT), DIMENSION(3) :: maxpos
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
INTEGER, INTENT(IN), DIMENSION(:,:,:) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Local value of integer input array
!*intmax*   INTEGER Global maximum on output
!*maxpos*   INTEGER Global grid indices of maximum
!*ivarid*   INTEGER Variable id
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/8)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, k, n, npcc, nreq, &
          & nx1loc, ny1loc
INTEGER, DIMENSION(1) :: kmax
INTEGER, DIMENSION(4) :: imaxtmp, imaxvals
INTEGER, DIMENSION(2*nprocs-2) :: irequest
INTEGER, DIMENSION(SIZE(intdat,DIM=3)) :: izmax
INTEGER, DIMENSION(2,SIZE(intdat,DIM=3)) :: maxlocs


procname(pglev+1) = 'maxloc_vars_int_3d'
CALL log_timer_in(npcc,ivarid)

!
!1. Value and location of local maximum
!--------------------------------------
!

IF (PRESENT(mglevel)) THEN
   nx1loc = mgvars(mglevel)%nc1loc; ny1loc = mgvars(mglevel)%nr1loc
ELSE
   nx1loc = nc1loc; ny1loc = nr1loc
ENDIF

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      k_110: DO k=1,SIZE(intdat,DIM=3)
         izmax(k) = MAXVAL(intdat(:,:,k),MASK=mask)
         maxlocs(:,k) = MAXLOC(intdat(:,:,k),MASK=mask)
      ENDDO k_110
      imaxvals(4)  = MAXVAL(izmax)
      kmax = MAXLOC(izmax)
      imaxvals(1) = maxlocs(1,kmax(1)) + nx1loc - 1
      imaxvals(2) = maxlocs(2,kmax(1)) + ny1loc - 1
      imaxvals(3) = kmax(1)
   ELSE
      imaxvals(4) = -ABS(int_fill)
      imaxvals(1:3) = 0
   ENDIF
ELSE
   imaxvals(4) = MAXVAL(intdat)
   imaxvals(1:3) = MAXLOC(intdat)
   imaxvals(1) = imaxvals(1) + nx1loc - 1
   imaxvals(2) = imaxvals(2) + ny1loc - 1
ENDIF

!---serial case
IF (iopt_MPI.EQ.0) THEN
   intmax = imaxvals(4)
   maxpos = imaxvals(1:3)
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
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
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_int(imaxvals,4,iddest,idloc,commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_int(imaxvals,4,iddest,idloc,commop,imode,&
                            & irequest(nreq),ivarid)
      ENDIF
   ELSE
      intmax = -ABS(int_fill)
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).EQ.iddest) THEN
            IF (imaxvals(4).GT.intmax) THEN
               intmax = imaxvals(4)
               maxpos = imaxvals(1:3)
            ENDIF
         ELSE
            IF (iform.EQ.1) THEN
               CALL comms_recv_int(imaxtmp,4,idprocs(iproc),idprocs(iproc),&
                                 & commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_int(imaxtmp,4,idprocs(iproc),idprocs(iproc),&
                                  & commop,irequest(nreq),ivarid)
            ENDIF
            IF (imaxtmp(4).GT.intmax) THEN
               intmax = imaxtmp(4)
               maxpos = imaxtmp(1:3)
            ENDIF
         ENDIF
      ENDDO iproc_410
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE

   intmax = -ABS(int_fill)
   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_int(imaxvals,4,idprocs(iproc),idloc,commop,imode,&
                              & ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_int(imaxvals,4,idprocs(iproc),idloc,commop,imode,&
                               & irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.2) THEN
         IF (iform.EQ.1) THEN
            CALL comms_recv_int(imaxtmp,4,idprocs(iproc),idprocs(iproc),&
                              & commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_int(imaxtmp,4,idprocs(iproc),idprocs(iproc),&
                               & commop,irequest(nreq),ivarid)
         ENDIF
         IF (imaxtmp(4).GT.intmax) THEN
            intmax = imaxtmp(4)
            maxpos = imaxtmp(1:3)
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         IF (imaxvals(4).GT.intmax) THEN
            intmax = imaxvals(4)
            maxpos = imaxvals(1:3)
         ENDIF
      ENDIF
   ENDDO n_510

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
1000 CALL log_timer_out(npcc,itm_com_util)


RETURN

END SUBROUTINE maxloc_vars_int_3d

!========================================================================

SUBROUTINE maxloc_vars_real_2d(realdat,realmax,maxpos,ivarid,mglevel,idroot,&
                             & comm,commtype,commall,mask)
!************************************************************************
!
! *maxloc_vars_real_2d* Compute the value and location of the global maximum
!                       of a 2-D real array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.8
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE gridpars
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
REAL, INTENT(OUT) :: realmax
INTEGER, INTENT(OUT), DIMENSION(2) :: maxpos
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
REAL, INTENT(IN), DIMENSION(:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local value of real input array
!*realmax*  REAL    Global maximum on output
!*maxpos*   INTEGER Global grid indices of maximum
!*ivarid*   INTEGER Variable id
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/8)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, n, npcc, nreq, nx1loc, &
         & ny1loc 
INTEGER, DIMENSION(2*nprocs-2) :: irequest
REAL, DIMENSION(3) :: rmaxtmp, rmaxvals


procname(pglev+1) = 'maxloc_vars_real_2d'
CALL log_timer_in(npcc,ivarid)

!
!1. Value and location of local maximum
!--------------------------------------
!

IF (PRESENT(mglevel)) THEN
   nx1loc = mgvars(mglevel)%nc1loc; ny1loc = mgvars(mglevel)%nr1loc
ELSE
   nx1loc = nc1loc; ny1loc = nr1loc
ENDIF

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      rmaxvals(3)  = MAXVAL(realdat,MASK=mask)
      rmaxvals(1:2) = MAXLOC(realdat,MASK=mask)
      rmaxvals(1) = rmaxvals(1) + nx1loc - 1
      rmaxvals(2) = rmaxvals(2) + ny1loc - 1
   ELSE
      rmaxvals(3) = -ABS(float_fill)
      rmaxvals(1:2) = 0.0
   ENDIF
ELSE
   rmaxvals(3) = MAXVAL(realdat)
   rmaxvals(1:2) = MAXLOC(realdat)
   rmaxvals(1) = rmaxvals(1) + nx1loc - 1
   rmaxvals(2) = rmaxvals(2) + ny1loc - 1
ENDIF

!---serial case
IF (iopt_MPI.EQ.0) THEN
   realmax = rmaxvals(3)
   maxpos = NINT(rmaxvals(1:2))
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
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
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_real(rmaxvals,3,iddest,idloc,commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_real(rmaxvals,3,iddest,idloc,commop,imode,&
                             & irequest(nreq),ivarid)
      ENDIF
   ELSE
      realmax = -ABS(float_fill)
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).EQ.iddest) THEN
            IF (rmaxvals(3).GT.realmax) THEN
               realmax = rmaxvals(3)
               maxpos = NINT(rmaxvals(1:2))
            ENDIF
         ELSE
            IF (iform.EQ.1) THEN
               CALL comms_recv_real(rmaxtmp,3,idprocs(iproc),idprocs(iproc),&
                                  & commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_real(rmaxtmp,3,idprocs(iproc),idprocs(iproc),&
                                   & commop,irequest(nreq),ivarid)
            ENDIF
            IF (rmaxtmp(3).GT.realmax) THEN
               realmax = rmaxtmp(3)
               maxpos = NINT(rmaxtmp(1:2))
            ENDIF
         ENDIF
      ENDDO iproc_410

   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE

   realmax = -ABS(float_fill)
   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_real(rmaxvals,3,idprocs(iproc),idloc,commop,imode,&
                               & ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_real(rmaxvals,3,idprocs(iproc),idloc,commop,&
                                & imode,irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.2) THEN
         IF (iform.EQ.1) THEN
            CALL comms_recv_real(rmaxtmp,3,idprocs(iproc),idprocs(iproc),&
                               & commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_real(rmaxtmp,3,idprocs(iproc),idprocs(iproc),&
                                & commop,irequest(nreq),ivarid)
         ENDIF
         IF (rmaxtmp(3).GT.realmax) THEN
            realmax = rmaxtmp(3)
            maxpos = NINT(rmaxtmp(1:2))
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         IF (rmaxvals(3).GT.realmax) THEN
            realmax = rmaxvals(3)
            maxpos = NINT(rmaxvals(1:2))
         ENDIF
      ENDIF
   ENDDO n_510

ENDIF

!
!6. Flag result if undefined
!---------------------------
!

IF (ABS(realmax+float_fill).LE.float_fill_eps) realmax = float_fill

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
1000 CALL log_timer_out(npcc,itm_com_util)


RETURN

END SUBROUTINE maxloc_vars_real_2d

!========================================================================

SUBROUTINE maxloc_vars_real_3d(realdat,realmax,maxpos,ivarid,mglevel,idroot,&
                             & comm,commtype,commall,mask)
!************************************************************************
!
! *maxloc_vars_real_3d* Compute the value and location of the global maximum
!                       of a real 3-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.8
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE gridpars
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
REAL, INTENT(OUT) :: realmax
INTEGER, INTENT(OUT), DIMENSION(3) :: maxpos
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
REAL, INTENT(IN), DIMENSION(:,:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local value of real input array
!*realmax*  REAL    Global maximum on output
!*maxpos*   INTEGER Global grid indices of maximum
!*ivarid*   INTEGER Variable id

!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/8)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, k, n, npcc, nreq, &
         & nx1loc, ny1loc
INTEGER, DIMENSION(1) :: kmax
INTEGER, DIMENSION(2,SIZE(realdat,DIM=3)) :: maxlocs
INTEGER, DIMENSION(2*nprocs-2) :: irequest
REAL, DIMENSION(SIZE(realdat,DIM=3)) :: rzmax
REAL, DIMENSION(4) :: rmaxtmp, rmaxvals


procname(pglev+1) = 'maxloc_vars_real_3d'
CALL log_timer_in(npcc,ivarid)

!
!1. Value and location of local maximum
!--------------------------------------
!

IF (PRESENT(mglevel)) THEN
   nx1loc = mgvars(mglevel)%nc1loc; ny1loc = mgvars(mglevel)%nr1loc
ELSE
   nx1loc = nc1loc; ny1loc = nr1loc
ENDIF

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      k_110: DO k=1,SIZE(realdat,DIM=3)
         rzmax(k) = MAXVAL(realdat(:,:,k),MASK=mask)
         maxlocs(:,k) = MAXLOC(realdat(:,:,k),MASK=mask)
      ENDDO k_110
      rmaxvals(4)  = MAXVAL(rzmax)
      kmax = MAXLOC(rzmax)
      rmaxvals(1) = maxlocs(1,kmax(1)) + nx1loc - 1
      rmaxvals(2) = maxlocs(2,kmax(1)) + ny1loc - 1
      rmaxvals(3) = kmax(1)
   ELSE
      rmaxvals(4) = -ABS(float_fill)
      rmaxvals(1:3) = 0.0
   ENDIF
ELSE
   rmaxvals(4) = MAXVAL(realdat)
   rmaxvals(1:3) = MAXLOC(realdat)
   rmaxvals(1) = rmaxvals(1) + nx1loc - 1
   rmaxvals(2) = rmaxvals(2) + ny1loc - 1
ENDIF

!---serial case
IF (iopt_MPI.EQ.0) THEN
   realmax = rmaxvals(4)
   maxpos = NINT(rmaxvals(1:3))
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
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
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_real(rmaxvals,4,iddest,idloc,commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_real(rmaxvals,4,iddest,idloc,commop,imode,&
                             & irequest(nreq),ivarid)
      ENDIF
   ELSE
      realmax = -ABS(float_fill)
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).EQ.iddest) THEN
            IF (rmaxvals(4).GT.realmax) THEN
               realmax = rmaxvals(4)
               maxpos = NINT(rmaxvals(1:3))
            ENDIF
         ELSE
            IF (iform.EQ.1) THEN
               CALL comms_recv_real(rmaxtmp,4,idprocs(iproc),idprocs(iproc),&
                                  & commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_real(rmaxtmp,4,idprocs(iproc),idprocs(iproc),&
                                   & commop,irequest(nreq),ivarid)
            ENDIF
            IF (rmaxtmp(4).GT.realmax) THEN
               realmax = rmaxtmp(4)
               maxpos = NINT(rmaxtmp(1:3))
            ENDIF
         ENDIF
      ENDDO iproc_410
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE

   realmax = -ABS(float_fill)
   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_real(rmaxvals,4,idprocs(iproc),idloc,commop,imode,&
                               & ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_real(rmaxvals,4,idprocs(iproc),idloc,commop,&
                                & imode,irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.2) THEN
         IF (iform.EQ.1) THEN
            CALL comms_recv_real(rmaxtmp,4,idprocs(iproc),idprocs(iproc),&
                               & commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_real(rmaxtmp,4,idprocs(iproc),idprocs(iproc),&
                                & commop,irequest(nreq),ivarid)
         ENDIF
         IF (rmaxtmp(4).GT.realmax) THEN
            realmax = rmaxtmp(4)
            maxpos = NINT(rmaxtmp(1:3))
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         IF (rmaxvals(4).GT.realmax) THEN
            realmax = rmaxvals(4)
            maxpos = NINT(rmaxvals(1:3))
         ENDIF
      ENDIF
   ENDDO n_510

ENDIF

!
!6. Flag result if undefined
!---------------------------
!

IF (ABS(realmax+float_fill).LE.float_fill_eps) realmax = float_fill

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
1000 CALL log_timer_out(npcc,itm_com_util)


RETURN

END SUBROUTINE maxloc_vars_real_3d

!========================================================================

SUBROUTINE min_vars_int_0d(intdat,intmin,idroot,ivarid,comm,commtype,commall)
!************************************************************************
!
! *min_vars_int_0d* Compute the global minimum of an integer scalar
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.7.1
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_int, comms_isend_int,
!                comms_recv_int, comms_reduce_int, comms_send_int,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_reduce_int, comms_send_int, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: intdat, ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
INTEGER, INTENT(OUT) :: intmin

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Local value of integer input variable
!*intmin*   INTEGER Global minimum on output
!*idroot*   INTEGER id of root process
!*ivarid*   INTEGER Variable id
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, n, npcc, nreq
INTEGER, DIMENSION(2*nprocs-2) :: irequest
INTEGER, DIMENSION(1) :: int1d


procname(pglev+1) = 'min_vars_int_0d'
CALL log_timer_in(npcc,ivarid)

!
!1. Serial case
!----------------
!

IF (iopt_MPI.EQ.0) THEN
   intmin = intdat
   GOTO 1000
END IF

!
!2. Optional arguments
!--------------------
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
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!4.1 Blocking/non-blocking
!-------------------------
!

   IF (iform.EQ.1.OR.iform.EQ.2) THEN

      IF (idloc.NE.iddest) THEN
         int1d(1) = intdat
         IF (iform.EQ.1) THEN
            CALL comms_send_int(int1d,1,iddest,idloc,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_int(int1d,1,iddest,idloc,commop,imode,&
                               & irequest(nreq),ivarid)
         ENDIF
      ELSE
         intmin = ABS(int_fill)
         iproc_410: DO iproc=1,nprocs
            IF (idprocs(iproc).EQ.iddest) THEN
               intmin = MIN(intdat,intmin)
            ELSE
               IF (iform.EQ.1) THEN
                  CALL comms_recv_int(int1d,1,idprocs(iproc),idprocs(iproc),&
                                    & commop,ivarid)
               ELSEIF (iform.EQ.2) THEN
                  nreq = nreq + 1
                  CALL comms_irecv_int(int1d,1,idprocs(iproc),idprocs(iproc),&
                                     & commop,irequest(nreq),ivarid)
               ENDIF
               intmin = MIN(int1d(1),intmin)
            ENDIF
         ENDDO iproc_410
      ENDIF

!
!4.2 Reduce
!----------
!
      
   ELSEIF (iform.EQ.3) THEN
      CALL comms_reduce_int(intdat,intmin,2,commop,iddest,ivarid)
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE

!
!5.1 Blocking/non-blocking
!-------------------------
!

   IF (iform.EQ.1.OR.iform.EQ.2) THEN

      intmin = ABS(int_fill)
      n_510: DO n=1,2*nprocs
         iproc = MOD(n-1,nprocs)+1
         IF (comprocs(n).EQ.1) THEN
            int1d(1) = intdat
            IF (iform.EQ.1) THEN
               CALL comms_send_int(int1d,1,idprocs(iproc),idloc,commop,imode,&
                                 & ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_isend_int(int1d,1,idprocs(iproc),idloc,commop,imode,&
                                  & irequest(nreq),ivarid)
            ENDIF
         ELSEIF (comprocs(n).EQ.2) THEN
            IF (iform.EQ.1) THEN
               CALL comms_recv_int(int1d,1,idprocs(iproc),idprocs(iproc),&
                                 & commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_int(int1d,1,idprocs(iproc),idprocs(iproc),&
                                  & commop,irequest(nreq),ivarid)
            ENDIF
            intmin = MIN(int1d(1),intmin)
         ELSEIF (comprocs(n).EQ.0) THEN
            intmin = MIN(intdat,intmin)
         ENDIF
      ENDDO n_510

!
!5.2 Reduce
!----------
!

   ELSEIF (iform.EQ.3) THEN
      CALL comms_reduce_int(intdat,intmin,2,commop,iddest,ivarid)
   ENDIF

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
1000 CALL log_timer_out(npcc,itm_com_util)


RETURN

END SUBROUTINE min_vars_int_0d

!========================================================================

SUBROUTINE min_vars_int_2d(intdat,intmin,ivarid,idroot,comm,commtype,commall,&
                         & mask)
!************************************************************************
!
! *min_vars_int_2d* Compute the global minimum of an integer 2-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.7.1
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - min_vars_int_0d
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(OUT) :: intmin
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
INTEGER, INTENT(IN), DIMENSION(:,:) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Local value of integer input array
!*intmin*   INTEGER Global minimum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iminval, itype


procname(pglev+1) = 'min_vars_int_2d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Local minimum
!----------------
!

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      iminval = MINVAL(intdat,MASK=mask)
   ELSE
      iminval = ABS(int_fill)
   ENDIF
ELSE
   iminval = MINVAL(intdat)
ENDIF

!
!2. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   intmin = iminval
   GOTO 1000
ENDIF

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

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!4. Domain maximum
!-----------------
!

CALL min_vars_int_0d(iminval,intmin,ivarid,iddest,commop,itype,all)

1000 CALL log_timer_out()


RETURN

END SUBROUTINE min_vars_int_2d

!========================================================================

SUBROUTINE min_vars_int_3d(intdat,intmin,ivarid,idroot,comm,commtype,commall,&
                         & mask)
!************************************************************************
!
! *min_vars_int_3d* Compute the global minimum of an integer 3-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.7.1
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - min_vars_int_0d
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(OUT) :: intmin
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
INTEGER, INTENT(IN), DIMENSION(:,:,:) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Local value of integer input array
!*intmin*   INTEGER Global minimum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iminval, itype, k
INTEGER, DIMENSION(SIZE(intdat,DIM=3)) :: izmin


procname(pglev+1) = 'min_vars_int_3d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Local minimum
!----------------
!

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      k_110: DO k=1,SIZE(intdat,DIM=3)
         izmin(k) = MINVAL(intdat(:,:,k),MASK=mask)
      ENDDO k_110
      iminval = MINVAL(izmin)
   ELSE
      iminval = ABS(int_fill)
   ENDIF
ELSE
   iminval = MINVAL(intdat)
ENDIF

!
!2. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   intmin = iminval
   GOTO 1000
ENDIF

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

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!4. Domain maximum
!-----------------
!

CALL min_vars_int_0d(iminval,intmin,ivarid,iddest,commop,itype,all)

1000 CALL log_timer_out()


RETURN

END SUBROUTINE min_vars_int_3d

!========================================================================

SUBROUTINE min_vars_real_0d(realdat,realmin,ivarid,idroot,comm,commtype,&
                          & commall)
!************************************************************************
!
! *min_vars_real_0d* Compute the global minimum of a real scalar
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.7.1
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_reduce_real, comms_send_real,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_reduce_real, comms_send_real, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(IN) :: realdat
REAL, INTENT(OUT) :: realmin

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local value of real input variable
!*realmin*  REAL    Global minimum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, n, npcc, nreq
INTEGER, DIMENSION(2*nprocs-2) :: irequest
REAL, DIMENSION(1) :: real1d


procname(pglev+1) = 'min_vars_real_0d'
CALL log_timer_in(npcc,ivarid)

!
!1. Serial case
!----------------
!

IF (iopt_MPI.EQ.0) THEN
   realmin = realdat
   GOTO 1000
END IF

!
!2. Optional arguments
!--------------------
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
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!4.1 Blocking/non-blocking
!-------------------------
!

   IF (iform.EQ.1.OR.iform.EQ.2) THEN
      IF (idloc.NE.iddest) THEN
         real1d(1) = realdat
         IF (iform.EQ.1) THEN
            CALL comms_send_real(real1d,1,iddest,idloc,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_real(real1d,1,iddest,idloc,commop,imode,&
                                & irequest(nreq),ivarid)
         ENDIF
      ELSE
         realmin = ABS(float_fill)
         iproc_410: DO iproc=1,nprocs
            IF (idprocs(iproc).EQ.iddest) THEN
               realmin = MIN(realdat,realmin)
            ELSE
               IF (iform.EQ.1) THEN
                  CALL comms_recv_real(real1d,1,idprocs(iproc),idprocs(iproc),&
                                     & commop,ivarid)
               ELSEIF (iform.EQ.2) THEN
                  nreq = nreq + 1
                  CALL comms_irecv_real(real1d,1,idprocs(iproc),&
                                      & idprocs(iproc),commop,irequest(nreq),&
                                      & ivarid)
               ENDIF
               realmin = MIN(real1d(1),realmin)
            ENDIF
         ENDDO iproc_410
      ENDIF

!
!4.2 Reduce
!----------
!
      
   ELSEIF (iform.EQ.3) THEN
      CALL comms_reduce_real(realdat,realmin,2,commop,iddest,ivarid)
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE

!
!5.1 Blocking/non-blocking
!-------------------------
!

   IF (iform.EQ.1.OR.iform.EQ.2) THEN

      realmin = ABS(float_fill)
      n_510: DO n=1,2*nprocs
         iproc = MOD(n-1,nprocs)+1
         IF (comprocs(n).EQ.1) THEN
            real1d(1) = realdat
            IF (iform.EQ.1) THEN
               CALL comms_send_real(real1d,1,idprocs(iproc),idloc,commop,&
                                  & imode,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_isend_real(real1d,1,idprocs(iproc),idloc,commop,&
                                   & imode,irequest(nreq),ivarid)
            ENDIF
         ELSEIF (comprocs(n).EQ.2) THEN
            IF (iform.EQ.1) THEN
               CALL comms_recv_real(real1d,1,idprocs(iproc),idprocs(iproc),&
                                  & commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_real(real1d,1,idprocs(iproc),idprocs(iproc),&
                                   & commop,irequest(nreq),ivarid)
            ENDIF
            realmin = MIN(real1d(1),realmin)
         ELSEIF (comprocs(n).EQ.0) THEN
            realmin = MIN(realdat,realmin)
         ENDIF
      ENDDO n_510

!
!5.2 Reduce
!----------
!

   ELSEIF (iform.EQ.3) THEN
      CALL comms_reduce_real(realdat,realmin,2,commop,iddest,ivarid)
   ENDIF

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
1000 CALL log_timer_out(npcc,itm_com_util)


RETURN

END SUBROUTINE min_vars_real_0d

!========================================================================

SUBROUTINE min_vars_real_2d(realdat,realmin,ivarid,idroot,comm,commtype,&
                          & commall,mask)
!************************************************************************
!
! *min_vars_real_2d* Compute the global minimum of a real 2-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.7.1
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - min_vars_real_0d
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(OUT) :: realmin
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
REAL, INTENT(IN), DIMENSION(:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local value of real input array
!*realmin*  REAL    Global minimum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, itype
REAL :: rminval


procname(pglev+1) = 'min_vars_real_2d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Local minimum
!----------------
!

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      rminval = MINVAL(realdat,MASK=mask)
   ELSE
      rminval = ABS(float_fill)
   ENDIF
ELSE
   rminval = MINVAL(realdat)
ENDIF

!
!2. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   realmin = rminval
   GOTO 1000
ENDIF

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

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!4. Domain maximum
!-----------------
!

CALL min_vars_real_0d(rminval,realmin,ivarid,iddest,commop,itype,all)

1000 CALL log_timer_out()


RETURN

END SUBROUTINE min_vars_real_2d

!========================================================================

SUBROUTINE min_vars_real_3d(realdat,realmin,ivarid,idroot,comm,commtype,&
                          & commall,mask)
!************************************************************************
!
! *min_vars_real_3d* Compute the global minimum of a real array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.7.1
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - min_vars_real_0d
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(OUT) :: realmin
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
REAL, INTENT(IN), DIMENSION(:,:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local value of real input array
!*realmin*  REAL    Global minimum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, itype, k
REAL :: rminval
REAL, DIMENSION(SIZE(realdat,DIM=3)) :: rzmin


procname(pglev+1) = 'min_vars_real_3d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Local minimum
!----------------
!

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      k_110: DO k=1,SIZE(realdat,DIM=3)
         rzmin(k) = MINVAL(realdat(:,:,k),MASK=mask)
      ENDDO k_110
      rminval = MINVAL(rzmin)
   ELSE
      rminval = ABS(float_fill)
   ENDIF
ELSE
   rminval = MINVAL(realdat)
ENDIF

!
!2. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   realmin = rminval
   GOTO 1000
ENDIF

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

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!4. Domain maximum
!-----------------
!

CALL min_vars_real_0d(rminval,realmin,ivarid,iddest,commop,itype,all)

1000 CALL log_timer_out()


RETURN

END SUBROUTINE min_vars_real_3d

!========================================================================

SUBROUTINE minloc_vars_int_2d(intdat,intmin,minpos,ivarid,mglevel,idroot,comm,&
                            & commtype,commall,mask)
!************************************************************************
!
! *minloc_vars_int_2d* Compute the value and location of the global minimum of
!                      an integer 2-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.8
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_waitall, error_abort
!
!************************************************************************
!
USE gridpars
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(OUT) :: intmin
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
INTEGER, INTENT(OUT), DIMENSION(2) :: minpos
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
INTEGER, INTENT(IN), DIMENSION(:,:) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Local value of integer input array
!*intmin*   INTEGER Global minimum on output
!*minpos*   INTEGER Global grid indices of minimum
!*ivarid*   INTEGER Variable id
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, n, npcc, nreq, nx1loc, &
        &  ny1loc
INTEGER, DIMENSION(3) :: imintmp, iminvals
INTEGER, DIMENSION(2*nprocs-2) :: irequest


procname(pglev+1) = 'minloc_vars_int_2d'
CALL log_timer_in(npcc,ivarid)

!
!1. Value and location of local minimum
!--------------------------------------
!

IF (PRESENT(mglevel)) THEN
   nx1loc = mgvars(mglevel)%nc1loc; ny1loc = mgvars(mglevel)%nr1loc
ELSE
   nx1loc = nc1loc; ny1loc = nr1loc
ENDIF

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      iminvals(3) = MINVAL(intdat,MASK=mask)
      iminvals(1:2) = MINLOC(intdat,MASK=mask)
      iminvals(1) = iminvals(1) + nx1loc - 1
      iminvals(2) = iminvals(2) + ny1loc - 1
   ELSE
      iminvals(3) = ABS(int_fill)
      iminvals(1:2) = 0
   ENDIF
ELSE
   iminvals(3) = MINVAL(intdat)
   iminvals(1:2) = MINLOC(intdat)
   iminvals(1) = iminvals(1) + nx1loc - 1
   iminvals(2) = iminvals(2) + ny1loc - 1
ENDIF

!---serial case
IF (iopt_MPI.EQ.0) THEN
   intmin = iminvals(3)
   minpos = iminvals(1:2)
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
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
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_int(iminvals,3,iddest,idloc,commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_int(iminvals,3,iddest,idloc,commop,imode,&
                            & irequest(nreq),ivarid)
      ENDIF
   ELSE
      intmin = ABS(int_fill)
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).EQ.iddest) THEN
            IF (iminvals(3).LT.intmin) THEN
               intmin = iminvals(3)
               minpos = iminvals(1:2)
            ENDIF
         ELSE
            IF (iform.EQ.1) THEN
               CALL comms_recv_int(imintmp,3,idprocs(iproc),idprocs(iproc),&
                                 & commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_int(imintmp,3,idprocs(iproc),idprocs(iproc),&
                                  & commop,irequest(nreq),ivarid)
            ENDIF
            IF (imintmp(3).LT.intmin) THEN
               intmin = imintmp(3)
               minpos = imintmp(1:2)
            ENDIF
         ENDIF
      ENDDO iproc_410
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE

   intmin = ABS(int_fill)
   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_int(iminvals,3,idprocs(iproc),idloc,commop,imode,&
                              & ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_int(iminvals,3,idprocs(iproc),idloc,commop,imode,&
                               & irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.2) THEN
         IF (iform.EQ.1) THEN
            CALL comms_recv_int(imintmp,3,idprocs(iproc),idprocs(iproc),&
                              & commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_int(imintmp,3,idprocs(iproc),idprocs(iproc),&
                               & commop,irequest(nreq),ivarid)
         ENDIF
         IF (imintmp(3).LT.intmin) THEN
            intmin = imintmp(3)
            minpos = imintmp(1:2)
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         IF (iminvals(3).LT.intmin) THEN
            intmin = iminvals(3)
            minpos = iminvals(1:2)
         ENDIF
      ENDIF
   ENDDO n_510

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
1000 CALL log_timer_out(npcc,itm_com_util)


RETURN

END SUBROUTINE minloc_vars_int_2d

!========================================================================

SUBROUTINE minloc_vars_int_3d(intdat,intmin,minpos,ivarid,mglevel,idroot,comm,&
                            & commtype,commall,mask)
!************************************************************************
!
! *minloc_vars_int_3d* Compute the value and location of the global minimum of
!                      an integer 3-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.8
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_int, comms_isend_int,
!                comms_recv_int, comms_send_int, comms_waitall, error_abort
!
!************************************************************************
!
USE gridpars
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_send_int, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(OUT) :: intmin
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
INTEGER, INTENT(OUT), DIMENSION(3) :: minpos
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
INTEGER, INTENT(IN), DIMENSION(:,:,:) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Local value of integer input array
!*intmin*   INTEGER Global minimum on output
!*minpos*   INTEGER Global grid indices of minimum
!*ivarid*   INTEGER Variable id
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, k, n, npcc, nreq, &
         & nx1loc, ny1loc
INTEGER, DIMENSION(1) :: kmin
INTEGER, DIMENSION(4) :: imintmp, iminvals
INTEGER, DIMENSION(SIZE(intdat,DIM=3)) :: izmin
INTEGER, DIMENSION(2*nprocs-2) :: irequest
INTEGER, DIMENSION(2,SIZE(intdat,DIM=3)) :: minlocs


procname(pglev+1) = 'minloc_vars_int_3d'
CALL log_timer_in(npcc,ivarid)

!
!1. Value and location of local minimum
!--------------------------------------
!

IF (PRESENT(mglevel)) THEN
   nx1loc = mgvars(mglevel)%nc1loc; ny1loc = mgvars(mglevel)%nr1loc
ELSE
   nx1loc = nc1loc; ny1loc = nr1loc
ENDIF

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      k_110: DO k=1,SIZE(intdat,DIM=3)
         izmin(k) = MINVAL(intdat(:,:,k),MASK=mask)
         minlocs(:,k) = MINLOC(intdat(:,:,k),MASK=mask)
      ENDDO k_110
      iminvals(4)  = MINVAL(izmin)
      kmin = MINLOC(izmin)
      iminvals(1) = minlocs(1,kmin(1)) + nx1loc - 1
      iminvals(2) = minlocs(2,kmin(1)) + ny1loc - 1
      iminvals(3) = kmin(1)
   ELSE
      iminvals(4) = ABS(int_fill)
      iminvals(1:3) = 0
   ENDIF
ELSE
   iminvals(4) = MINVAL(intdat)
   iminvals(1:3) = MINLOC(intdat)
   iminvals(1) = iminvals(1) + nx1loc - 1
   iminvals(2) = iminvals(2) + ny1loc - 1
ENDIF

!---serial case
IF (iopt_MPI.EQ.0) THEN
   intmin = iminvals(4)
   minpos = iminvals(1:3)
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
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
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_int(iminvals,4,iddest,idloc,commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_int(iminvals,4,iddest,idloc,commop,imode,&
                            & irequest(nreq),ivarid)
      ENDIF
   ELSE
      intmin = ABS(int_fill)
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).EQ.iddest) THEN
            IF (iminvals(4).LT.intmin) THEN
               intmin = iminvals(4)
               minpos = iminvals(1:3)
            ENDIF
         ELSE
            IF (iform.EQ.1) THEN
               CALL comms_recv_int(imintmp,4,idprocs(iproc),idprocs(iproc),&
                                 & commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_int(imintmp,4,idprocs(iproc),idprocs(iproc),&
                                  & commop,irequest(nreq),ivarid)
            ENDIF
            IF (imintmp(4).LT.intmin) THEN
               intmin = imintmp(4)
               minpos = imintmp(1:3)
            ENDIF
         ENDIF
      ENDDO iproc_410
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE

   intmin = ABS(int_fill)
   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_int(iminvals,4,idprocs(iproc),idloc,commop,imode,&
                              & ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_int(iminvals,4,idprocs(iproc),idloc,commop,imode,&
                               & irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.2) THEN
         IF (iform.EQ.1) THEN
            CALL comms_recv_int(imintmp,4,idprocs(iproc),idprocs(iproc),&
                              & commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_int(imintmp,4,idprocs(iproc),idprocs(iproc),&
                               & commop,irequest(nreq),ivarid)
         ENDIF
         IF (imintmp(4).LT.intmin) THEN
            intmin = imintmp(4)
            minpos = imintmp(1:3)
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         IF (iminvals(4).LT.intmin) THEN
            intmin = iminvals(4)
            minpos = iminvals(1:3)
         ENDIF
      ENDIF
   ENDDO n_510

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
1000 CALL log_timer_out(npcc,itm_com_util)


RETURN

END SUBROUTINE minloc_vars_int_3d

!========================================================================

SUBROUTINE minloc_vars_real_2d(realdat,realmin,minpos,ivarid,mglevel,idroot,&
                             & comm,commtype,commall,mask)
!************************************************************************
!
! *minloc_vars_real_2d* Compute the value and location of the global minimum
!                       of a 2-D real array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.8
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall, error_abort
!
!************************************************************************
!
USE gridpars
USE multigrid
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_send_real, comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
INTEGER, INTENT(IN) :: ivarid
REAL, INTENT(OUT) :: realmin
INTEGER, INTENT(OUT), DIMENSION(2) :: minpos
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
REAL, INTENT(IN), DIMENSION(:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local value of real input array
!*realmin*  REAL    Global minimum on output
!*minpos*   INTEGER Global grid indices of minimum
!*ivarid*   INTEGER Variable id
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/8)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, n, npcc, nreq, nx1loc, &
         & ny1loc
INTEGER, DIMENSION(2*nprocs-2) :: irequest
REAL, DIMENSION(3) :: rmintmp, rminvals


procname(pglev+1) = 'minloc_vars_real_2d'
CALL log_timer_in(npcc,ivarid)

!
!1. Value and location of local minimum
!--------------------------------------
!

IF (PRESENT(mglevel)) THEN
   nx1loc = mgvars(mglevel)%nc1loc; ny1loc = mgvars(mglevel)%nr1loc
ELSE
   nx1loc = nc1loc; ny1loc = nr1loc
ENDIF

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      rminvals(3) = MINVAL(realdat,MASK=mask)
      rminvals(1:2) = MINLOC(realdat,MASK=mask)
      rminvals(1) = rminvals(1) + nx1loc - 1
      rminvals(2) = rminvals(2) + ny1loc - 1
   ELSE
      rminvals(3) = ABS(float_fill)
      rminvals(1:2) = 0.0
   ENDIF
ELSE
   rminvals(3) = MINVAL(realdat)
   rminvals(1:2) = MINLOC(realdat)
   rminvals(1) = rminvals(1) + nx1loc - 1
   rminvals(2) = rminvals(2) + ny1loc - 1
ENDIF

!---serial case
IF (iopt_MPI.EQ.0) THEN
   realmin = rminvals(3)
   minpos = NINT(rminvals(1:2))
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
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
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_real(rminvals,3,iddest,idloc,commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_real(rminvals,3,iddest,idloc,commop,imode,&
                             & irequest(nreq),ivarid)
      ENDIF
   ELSE
      realmin = ABS(float_fill)
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).EQ.iddest) THEN
            IF (rminvals(3).LT.realmin) THEN
               realmin = rminvals(3)
               minpos = NINT(rminvals(1:2))
            ENDIF
         ELSE
            IF (iform.EQ.1) THEN
               CALL comms_recv_real(rmintmp,3,idprocs(iproc),idprocs(iproc),&
                                  & commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_real(rmintmp,3,idprocs(iproc),idprocs(iproc),&
                                   & commop,irequest(nreq),ivarid)
            ENDIF
            IF (rmintmp(3).LT.realmin) THEN
               realmin = rmintmp(3)
               minpos = NINT(rmintmp(1:2))
            ENDIF
         ENDIF
      ENDDO iproc_410
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE

   realmin = ABS(float_fill)
   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_real(rminvals,3,idprocs(iproc),idloc,commop,imode,&
                               & ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_real(rminvals,3,idprocs(iproc),idloc,commop,&
                                & imode,irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.2) THEN
         IF (iform.EQ.1) THEN
            CALL comms_recv_real(rmintmp,3,idprocs(iproc),idprocs(iproc),&
                              & commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_real(rmintmp,3,idprocs(iproc),idprocs(iproc),&
                                & commop,irequest(nreq),ivarid)
         ENDIF
         IF (rmintmp(3).LT.realmin) THEN
            realmin = rmintmp(3)
            minpos = NINT(rmintmp(1:2))
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         IF (rminvals(3).LT.realmin) THEN
            realmin = rminvals(3)
            minpos = NINT(rminvals(1:2))
         ENDIF
      ENDIF
   ENDDO n_510

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
1000 CALL log_timer_out(npcc,itm_com_util)


RETURN

END SUBROUTINE minloc_vars_real_2d

!========================================================================

SUBROUTINE minloc_vars_real_3d(realdat,realmin,minpos,ivarid,mglevel,idroot,&
                             & comm,commtype,commall,mask)
!************************************************************************
!
! *minloc_vars_real_3d* Compute the value and location of the global minimum
!                       of a real 3-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.8
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_send_real, comms_waitall
!
!************************************************************************
!
USE gridpars
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
REAL, INTENT(OUT) :: realmin
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
INTEGER, INTENT(OUT), DIMENSION(3) :: minpos
REAL, INTENT(IN), DIMENSION(:,:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local value of real input array
!*realmin*  REAL    Global minimum on output
!*minpos*   INTEGER Global grid indices of minimum
!*ivarid*   INTEGER Variable id
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/8)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, k, n, npcc, nreq, &
         & nx1loc, ny1loc
INTEGER, DIMENSION(1) :: kmin
INTEGER, DIMENSION(2*nprocs-2) :: irequest
INTEGER, DIMENSION(2,SIZE(realdat,DIM=3)) :: minlocs
REAL, DIMENSION(4) :: rmintmp, rminvals
REAL, DIMENSION(SIZE(realdat,DIM=3)) :: rzmin


procname(pglev+1) = 'minloc_vars_real_3d'
CALL log_timer_in(npcc,ivarid)

!
!1. Value and location of local minimum
!--------------------------------------
!

IF (PRESENT(mglevel)) THEN
   nx1loc = mgvars(mglevel)%nc1loc; ny1loc = mgvars(mglevel)%nr1loc
ELSE
   nx1loc = nc1loc; ny1loc = nr1loc
ENDIF

IF (PRESENT(mask)) THEN
   IF (ANY(mask)) THEN
      k_110: DO k=1,SIZE(realdat,DIM=3)
         rzmin(k) = MINVAL(realdat(:,:,k),MASK=mask)
         minlocs(:,k) = MINLOC(realdat(:,:,k),MASK=mask)
      ENDDO k_110
      rminvals(4)  = MINVAL(rzmin)
      kmin = MINLOC(rzmin)
      rminvals(1) = minlocs(1,kmin(1)) + nx1loc - 1
      rminvals(2) = minlocs(2,kmin(1)) + ny1loc - 1
      rminvals(3) = kmin(1)
   ELSE
      rminvals(4) = ABS(float_fill)
      rminvals(1:3) = 0.0
   ENDIF
ELSE
   rminvals(4) = MINVAL(realdat)
   rminvals(1:3) = MINLOC(realdat)
   rminvals(1) = rminvals(1) + nx1loc - 1
   rminvals(2) = rminvals(2) + ny1loc - 1
ENDIF

!---serial case
IF (iopt_MPI.EQ.0) THEN
   realmin = rminvals(4)
   minpos = NINT(rminvals(1:3))
   GOTO 1000
ENDIF

!
!2. Optional arguments
!--------------------
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
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

   IF (idloc.NE.iddest) THEN
      IF (iform.EQ.1) THEN
         CALL comms_send_real(rminvals,4,iddest,idloc,commop,imode,ivarid)
      ELSEIF (iform.EQ.2) THEN
         nreq = nreq + 1
         CALL comms_isend_real(rminvals,4,iddest,idloc,commop,imode,&
                             & irequest(nreq),ivarid)
      ENDIF
   ELSE
      realmin = ABS(float_fill)
      iproc_410: DO iproc=1,nprocs
         IF (idprocs(iproc).EQ.iddest) THEN
            IF (rminvals(4).LT.realmin) THEN
               realmin = rminvals(4)
               minpos = NINT(rminvals(1:3))
            ENDIF
         ELSE
            IF (iform.EQ.1) THEN
               CALL comms_recv_real(rmintmp,4,idprocs(iproc),idprocs(iproc),&
                                  & commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_real(rmintmp,4,idprocs(iproc),idprocs(iproc),&
                                   & commop,irequest(nreq),ivarid)
            ENDIF
            IF (rmintmp(4).LT.realmin) THEN
               realmin = rmintmp(4)
               minpos = NINT(rmintmp(1:3))
            ENDIF
         ENDIF
      ENDDO iproc_410
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE

   realmin = ABS(float_fill)
   n_510: DO n=1,2*nprocs
      iproc = MOD(n-1,nprocs)+1
      IF (comprocs(n).EQ.1) THEN
         IF (iform.EQ.1) THEN
            CALL comms_send_real(rminvals,4,idprocs(iproc),idloc,commop,imode,&
                               & ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_real(rminvals,4,idprocs(iproc),idloc,commop,&
                                & imode,irequest(nreq),ivarid)
         ENDIF
      ELSEIF (comprocs(n).EQ.2) THEN
         IF (iform.EQ.1) THEN
            CALL comms_recv_real(rmintmp,4,idprocs(iproc),idprocs(iproc),&
                               & commop,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_irecv_real(rmintmp,4,idprocs(iproc),idprocs(iproc),&
                                & commop,irequest(nreq),ivarid)
         ENDIF
         IF (rmintmp(4).LT.realmin) THEN
            realmin = rmintmp(4)
            minpos = NINT(rmintmp(1:3))
         ENDIF
      ELSEIF (comprocs(n).EQ.0) THEN
         IF (rminvals(4).LT.realmin) THEN
            realmin = rminvals(4)
            minpos = NINT(rminvals(1:3))
         ENDIF
      ENDIF
   ENDDO n_510

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
1000 CALL log_timer_out(npcc,itm_com_util)


RETURN

END SUBROUTINE minloc_vars_real_3d

!========================================================================

SUBROUTINE sum_vars_int_0d(intdat,intsum,ivarid,idroot,comm,commtype,commall)
!************************************************************************
!
! *sum_vars_int_0d* Compute the global sum of an integer scalar
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.0
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_int, comms_isend_int,
!                comms_recv_int, comms_reduce_int, comms_send_int,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_int, comms_isend_int, &
                   & comms_recv_int, comms_reduce_int, comms_send_int, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: intdat, ivarid
INTEGER, INTENT(OUT) :: intsum
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Local value of integer input variable
!*intsum*   INTEGER Global sum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, n, npcc, nreq
INTEGER, DIMENSION(2*nprocs-2) :: irequest
INTEGER, DIMENSION(1) :: int1d


procname(pglev+1) = 'sum_vars_int_0d'
CALL log_timer_in(npcc,ivarid)

!
!1. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   intsum = intdat
   GOTO 1000
END IF

!
!2. Optional arguments
!--------------------
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
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!4.1 Blocking/non-blocking
!-------------------------
!

   IF (iform.EQ.1.OR.iform.EQ.2) THEN

      IF (idloc.NE.iddest) THEN
         int1d(1) = intdat
         IF (iform.EQ.1) THEN
            CALL comms_send_int(int1d,1,iddest,idloc,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_int(int1d,1,iddest,idloc,commop,imode,&
                               & irequest(nreq),ivarid)
         ENDIF
      ELSE
         intsum = 0
         iproc_410: DO iproc=1,nprocs
            IF (idprocs(iproc).EQ.iddest) THEN
               intsum = intsum + intdat
            ELSE
               IF (iform.EQ.1) THEN
                  CALL comms_recv_int(int1d,1,idprocs(iproc),idprocs(iproc),&
                                    & commop,ivarid)
               ELSEIF (iform.EQ.2) THEN
                  nreq = nreq + 1
                  CALL comms_irecv_int(int1d,1,idprocs(iproc),idprocs(iproc),&
                                     & commop,irequest(nreq),ivarid)
               ENDIF
               intsum = intsum + int1d(1)
            ENDIF
         ENDDO iproc_410
      ENDIF

!
!4.2 Reduce
!----------
!
      
   ELSEIF (iform.EQ.3) THEN
      CALL comms_reduce_int(intdat,intsum,3,commop,iddest,ivarid)
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE

!
!5.1 Blocking/non-blocking
!-------------------------
!

   IF (iform.EQ.1.OR.iform.EQ.2) THEN
      intsum = 0
      n_510: DO n=1,2*nprocs
         iproc = MOD(n-1,nprocs)+1
         IF (comprocs(n).EQ.1) THEN
            int1d(1) = intdat
            IF (iform.EQ.1) THEN
               CALL comms_send_int(int1d,1,idprocs(iproc),idloc,commop,imode,&
                                 & ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_isend_int(int1d,1,idprocs(iproc),idloc,commop,imode,&
                                  & irequest(nreq),ivarid)
            ENDIF
         ELSEIF (comprocs(n).EQ.2) THEN
            IF (iform.EQ.1) THEN
               CALL comms_recv_int(int1d,1,idprocs(iproc),idprocs(iproc),&
                                 & commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_int(int1d,1,idprocs(iproc),idprocs(iproc),&
                                  & commop,irequest(nreq),ivarid)
            ENDIF
            intsum = intsum + int1d(1)
         ELSEIF (comprocs(n).EQ.0) THEN
            intsum = intsum + intdat
         ENDIF
      ENDDO n_510

!
!5.2 Reduce
!----------
!
      
   ELSEIF (iform.EQ.3) THEN
      CALL comms_reduce_int(intdat,intsum,3,commop,iddest,ivarid)
   ENDIF

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
1000 CALL log_timer_out(npcc,itm_com_util)


RETURN

END SUBROUTINE sum_vars_int_0d

!========================================================================

SUBROUTINE sum_vars_int_2d(intdat,intsum,ivarid,idroot,comm,commtype,commall,&
                         & mask)
!************************************************************************
!
! *sum_vars_int_2d* Compute the global sum of an integer 2-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.0
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - sum_vars_int_0d
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(OUT) :: intsum
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
INTEGER, INTENT(IN), DIMENSION(:,:) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Local value of integer input array
!*intsum*   INTEGER Global sum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, itype, ival


procname(pglev+1) = 'sum_vars_int_2d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Local sum
!------------
!

IF (PRESENT(mask)) THEN
   ival = MERGE(SUM(intdat,MASK=mask),0,ANY(mask))
ELSE
   ival = SUM(intdat)
ENDIF

!
!2. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   intsum = ival
   GOTO 1000
ENDIF

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

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!4. Domain maximum
!-----------------
!

CALL sum_vars_int_0d(ival,intsum,ivarid,iddest,commop,itype,all)

1000 CALL log_timer_out()


RETURN

END SUBROUTINE sum_vars_int_2d

!========================================================================

SUBROUTINE sum_vars_int_3d(intdat,intsum,ivarid,idroot,comm,commtype,commall,&
                         & mask)
!************************************************************************
!
! *sum_vars_int_3d* Compute the global sum of an integer 3-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.0
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - sum_vars_int_0d
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(OUT) :: intsum
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
INTEGER, INTENT(IN), DIMENSION(:,:,:) :: intdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER Local value of integer input array
!*intsum*   INTEGER Global sum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, itype, ival, k
INTEGER, DIMENSION(SIZE(intdat,DIM=3)) :: izsum


procname(pglev+1) = 'sum_vars_int_3d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Local sum
!------------
!

k_110: DO k=1,SIZE(intdat,DIM=3)
   IF (PRESENT(mask)) THEN
      izsum(k) = MERGE(SUM(intdat(:,:,k),MASK=mask),0,ANY(mask))
   ELSE
      izsum(k) = SUM(intdat(:,:,k))
   ENDIF
ENDDO k_110
ival = SUM(izsum)

!
!2. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   intsum = ival
   GOTO 1000
ENDIF

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

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!4. Domain maximum
!-----------------
!

CALL sum_vars_int_0d(ival,intsum,ivarid,iddest,commop,itype,all)

1000 CALL log_timer_out()


RETURN

END SUBROUTINE sum_vars_int_3d

!========================================================================

SUBROUTINE sum_vars_real_0d(realdat,realsum,ivarid,idroot,comm,commtype,&
                          & commall)
!************************************************************************
!
! *sum_vars_real_0d* Compute the global sum of a real scalar
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.0
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - comms_barrier, comms_irecv_real, comms_isend_real,
!                comms_recv_real, comms_reduce_real, comms_send_real,
!                comms_waitall, error_abort
!
!************************************************************************
!
USE comms_MPI, ONLY: comms_barrier, comms_irecv_real, comms_isend_real, &
                   & comms_recv_real, comms_reduce_real, comms_send_real, &
                   & comms_waitall
USE error_routines, ONLY: error_abort

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(IN) :: realdat
REAL, INTENT(OUT) :: realsum

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local value of real input variable
!*realsum*  REAL    Global sum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, iform, imode, iproc, itype, n, npcc, nreq
INTEGER, DIMENSION(2*nprocs-2) :: irequest
REAL, DIMENSION(1) :: real1d


procname(pglev+1) = 'sum_vars_real_0d'
CALL log_timer_in(npcc,ivarid)

!
!1. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   realsum = realdat
   GOTO 1000
END IF

!
!2. Optional arguments
!--------------------
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
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!3. Parameters for communication
!-------------------------------
!

imode = MOD(itype-1,2) + 1
iform = (itype-1)/2 + 1
nreq = 0

!
!4. Result at root process only
!------------------------------
!

IF (.NOT.all) THEN

!
!5.1 Blocking/non-blocking
!-------------------------
!

   IF (iform.EQ.1.OR.iform.EQ.2) THEN

      IF (idloc.NE.iddest) THEN
         IF (iform.EQ.1) THEN
            real1d(1) = realdat
            CALL comms_send_real(real1d,1,iddest,idloc,commop,imode,ivarid)
         ELSEIF (iform.EQ.2) THEN
            nreq = nreq + 1
            CALL comms_isend_real(real1d,1,iddest,idloc,commop,imode,&
                                & irequest(nreq),ivarid)
         ENDIF
      ELSE
         realsum = 0.0
         iproc_410: DO iproc=1,nprocs
            IF (idprocs(iproc).EQ.iddest) THEN
               realsum = realsum + realdat
            ELSE
               IF (iform.EQ.1) THEN
                  CALL comms_recv_real(real1d,1,idprocs(iproc),idprocs(iproc),&
                                     & commop,ivarid)
               ELSEIF (iform.EQ.2) THEN
                  nreq = nreq + 1
                  CALL comms_irecv_real(real1d,1,idprocs(iproc),&
                                      & idprocs(iproc),commop,irequest(nreq),&
                                      & ivarid)
               ENDIF
               realsum = realsum + real1d(1)
            ENDIF
         ENDDO iproc_410
      ENDIF

!
!4.2 Reduce
!----------
!
      
   ELSEIF (iform.EQ.3) THEN
      CALL comms_reduce_real(realdat,realsum,3,commop,iddest,ivarid)
   ENDIF

!
!5. Result on all processes
!--------------------------
!

ELSE

!
!5.1 Blocking/non-blocking
!-------------------------
!

   IF (iform.EQ.1.OR.iform.EQ.2) THEN
      realsum = 0.0
      n_510: DO n=1,2*nprocs
         iproc = MOD(n-1,nprocs)+1
         IF (comprocs(n).EQ.1) THEN
            real1d(1) = realdat
            IF (iform.EQ.1) THEN
               CALL comms_send_real(real1d,1,idprocs(iproc),idloc,commop,&
                                  & imode,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_isend_real(real1d,1,idprocs(iproc),idloc,commop,&
                                   & imode,irequest(nreq),ivarid)
            ENDIF
         ELSEIF (comprocs(n).EQ.2) THEN
            IF (iform.EQ.1) THEN
               CALL comms_recv_real(real1d,1,idprocs(iproc),idprocs(iproc),&
                                  & commop,ivarid)
            ELSEIF (iform.EQ.2) THEN
               nreq = nreq + 1
               CALL comms_irecv_real(real1d,1,idprocs(iproc),idprocs(iproc),&
                                   & commop,irequest(nreq),ivarid)
            ENDIF
            realsum = realsum + real1d(1)
         ELSEIF (comprocs(n).EQ.0) THEN
            realsum = realsum + realdat
         ENDIF
      ENDDO n_510

!
!5.2 Reduce
!----------
!
      
   ELSEIF (iform.EQ.3) THEN
      CALL comms_reduce_real(realdat,realsum,3,commop,iddest,ivarid)
   ENDIF

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
1000 CALL log_timer_out(npcc,itm_com_util)


RETURN

END SUBROUTINE sum_vars_real_0d

!========================================================================

SUBROUTINE sum_vars_real_2d(realdat,realsum,ivarid,idroot,comm,commtype,&
                          & commall,mask)
!************************************************************************
!
! *sum_vars_real_2d* Compute the global sum of a real 2-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.0
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - sum_vars_real_0d
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(OUT) :: realsum
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
REAL, INTENT(IN), DIMENSION(:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local value of real input array
!*realsum*  REAL    Global sum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all 
INTEGER :: commop, iddest, itype
REAL :: rval


procname(pglev+1) = 'sum_vars_real_2d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Local sum
!------------
!

IF (PRESENT(mask)) THEN
   rval = MERGE(SUM(realdat,MASK=mask),0.0,ANY(mask))
ELSE
   rval = SUM(realdat)
ENDIF

!
!2. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   realsum = rval
   GOTO 1000
ENDIF

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

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!4. Domain sum
!-------------
!

CALL sum_vars_real_0d(rval,realsum,ivarid,iddest,commop,itype,all)

1000 CALL log_timer_out()


RETURN

END SUBROUTINE sum_vars_real_2d

!========================================================================

SUBROUTINE sum_vars_real_3d(realdat,realsum,ivarid,idroot,comm,commtype,&
                          & commall,mask)
!************************************************************************
!
! *sum_vars_real_3d* Compute the global sum of a real 3-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.0
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - sum_vars_real_0d
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot
REAL, INTENT(OUT) :: realsum
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
REAL, INTENT(IN), DIMENSION(:,:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local value of real input array
!*realsum*  REAL    Global sum on output
!*ivarid*   INTEGER Variable id
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for send/receive or reduce operations
!*commtype* INTEGER Type of communication (1/9)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all 
INTEGER :: commop, iddest, itype, k
REAL :: rval
REAL, DIMENSION(SIZE(realdat,DIM=3)) :: rzsum


procname(pglev+1) = 'sum_vars_real_3d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Local sum
!------------
!

k_110: DO k=1,SIZE(realdat,DIM=3)
   IF (PRESENT(mask)) THEN
      rzsum(k) = MERGE(SUM(realdat(:,:,k),MASK=mask),0.0,ANY(mask))
   ELSE
      rzsum(k) = SUM(realdat(:,:,k))
   ENDIF
ENDDO k_110
rval = SUM(rzsum)

!
!2. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   realsum = rval
   GOTO 1000
ENDIF

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

IF (PRESENT(commall)) THEN
   all = commall
ELSE
   all = .FALSE.
ENDIF

IF (PRESENT(commtype)) THEN
   itype = commtype
ELSE
   IF (all) THEN
      itype = MERGE(5,iopt_MPI_comm_all,iopt_MPI_comm_coll.EQ.1)
   ELSE
      itype = MERGE(5,iopt_MPI_comm_gath,iopt_MPI_comm_coll.EQ.1)
   ENDIF
ENDIF

!
!4. Domain sum
!-------------
!

CALL sum_vars_real_0d(rval,realsum,ivarid,iddest,commop,itype,all)

1000 CALL log_timer_out()


RETURN

END SUBROUTINE sum_vars_real_3d

!========================================================================

SUBROUTINE sum2_vars_real_2d(realdat,realsum,nhdims,cnode,ivarid,mglevel,&
                           & idroot,comm,commtype,commall,mask)
!************************************************************************
!
! *sum2_vars_real_2d* Process independent computation of the global sum of a
!                     real 2-D model array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.8
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - combine_mod, error_alloc
!
!************************************************************************
!
USE grid
USE gridpars
USE multigrid
USE syspars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: combine_mod

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
REAL, INTENT(OUT) :: realsum
INTEGER, INTENT(IN), DIMENSION(4) :: nhdims
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mask
REAL, INTENT(IN), DIMENSION(:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local model array
!*realsum*  REAL    Global sum on output
!*nhdims*   INTEGER Halo sizes (WESN directions)   
!*cnode*    CHAR    Grid node of model array ('0', 'C', 'U', 'V', 'E')
!*ivarid*   INTEGER Variable id
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for combine operations
!*commtype* INTEGER Type of communication (1/8)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points if PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, itype, lev, nx, nxloc, ny, nyloc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskloc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: realloc, realglb


procname(pglev+1) = 'sum2_vars_real_2d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
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
!2. Initialise parameters
!------------------------
!

IF (lev.EQ.0) THEN
   nx = nc; ny = nr
   nxloc = ncloc; nyloc = nrloc
ELSE
   nx = mgvars(lev)%nc; ny = mgvars(lev)%nr
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
ENDIF

!
!3. Local mask array
!-------------------
!
!---allocate
ALLOCATE (maskloc(nxloc,nyloc),STAT=errstat)
CALL error_alloc('maskloc',2,(/nxloc,nyloc/),kndlog)

!---define
IF (PRESENT(mask)) THEN
   maskloc = mask
ELSE
   IF (lev.EQ.0) THEN
      SELECT CASE (TRIM(cnode))
         CASE ('0')
            maskloc = .TRUE.
         CASE ('C','W')
            maskloc = nodeatc(1:nxloc,1:nyloc).GT.0
         CASE ('U')
            maskloc = node2du(1:nxloc,1:nyloc).GT.1
         CASE ('V')
            maskloc = node2dv(1:nxloc,1:nyloc).GT.1
         CASE ('UV')
            maskloc = node2duv(1:nxloc,1:nyloc).GT.0
      END SELECT
   ELSE
      SELECT CASE (TRIM(cnode))
      CASE ('0')
         maskloc = .TRUE.
      CASE ('C','W')
         maskloc = mgvars(lev)%nodeatc(1:nxloc,1:nyloc).GT.0
      CASE ('U')
         maskloc = mgvars(lev)%node2du(1:nxloc,1:nyloc).GT.1
      CASE ('V')
         maskloc = mgvars(lev)%node2dv(1:nxloc,1:nyloc).GT.1
      END SELECT
   ENDIF
ENDIF

!
!4. Serial case
!--------------
!

IF (iopt_MPI.EQ.0) THEN
   realsum = MERGE(SUM(realdat(1:nxloc,1:nyloc),MASK=maskloc),0.0,ANY(maskloc))
   GOTO 1000
ENDIF

!
!5. Local storage array
!----------------------
!
!---allocate
ALLOCATE (realloc(1-nhdims(1):nxloc+nhdims(2),1-nhdims(3):nyloc+nhdims(4)),&
        & STAT=errstat)
CALL error_alloc('realloc',2,(/nxloc+nhdims(1)+nhdims(2),&
                             & nyloc+nhdims(3)+nhdims(4)/),kndrtype)

!---define
realloc = realdat
WHERE (.NOT.maskloc)
   realloc(1:nxloc,1:nyloc) = float_fill
END WHERE

!
!6. Global array
!---------------
!
!---allocate
ALLOCATE (realglb(1-nhdims(1):nx+nhdims(2),1-nhdims(3):ny+nhdims(4)),&
        & STAT=errstat)
CALL error_alloc('realglb',2,(/nx+nhdims(1)+nhdims(2),&
                             & ny+nhdims(3)+nhdims(4)/),kndrtype)

!---combine
CALL combine_mod(realglb,realloc,(/1-nhdims(1),1-nhdims(3)/),ivarid,0.0,&
               & lev,iddest,commop,itype,all)

!
!7. Global sum
!-------------
!

IF (all.OR.idloc.EQ.iddest) THEN
   IF (ANY(ABS(realglb(1:nc,1:nr)-float_fill).GT.float_fill_eps)) THEN
      realsum = SUM(realglb(1:nc,1:nr),&
                  & MASK=ABS(realglb(1:nc,1:nr)-float_fill).GT.float_fill_eps)
   ELSE
      realsum = 0.0
   ENDIF
ENDIF

!
!8. Deallocate arrays
!--------------------
!

DEALLOCATE (realloc,realglb)
1000 DEALLOCATE (maskloc)

CALL log_timer_out()


RETURN

END SUBROUTINE sum2_vars_real_2d

!========================================================================

SUBROUTINE sum2_vars_real_3d(realdat,realsum,nhdims,cnode,ivarid,mglevel,&
                           & idroot,comm,commtype,commall,mask)
!************************************************************************
!
! *sum2_vars_real_3d* Process independent computation of the global sum of a
!                     real 3-D model array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.8
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - combine_mod, error_alloc
!
!************************************************************************
!
USE grid
USE gridpars
USE multigrid
USE syspars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: combine_mod

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
REAL, INTENT(OUT) :: realsum
INTEGER, INTENT(IN), DIMENSION(4) :: nhdims
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) :: mask
REAL, INTENT(IN), DIMENSION(:,:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local model array
!*realsum*  REAL    Global sum on output
!*nhdims*   INTEGER Halo sizes (WESN directions)   
!*cnode*    CHAR    Grid node of model array ('0', 'C', 'U', 'V', 'E')
!*ivarid*   INTEGER Variable id
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for combine operations
!*commtype* INTEGER Type of communication (1/8)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points if PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, itype, k, lev, nx, nxloc, ny, nyloc, n3dim
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskloc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realloc, realglb


procname(pglev+1) = 'sum2_vars_real_3d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
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
!2. Initialise parameters
!------------------------
!

IF (lev.EQ.0) THEN
   nx = nc; ny = nr
   nxloc = ncloc; nyloc = nrloc
ELSE
   nx = mgvars(lev)%nc; ny = mgvars(lev)%nr
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
ENDIF

!
!3. Local mask array
!-------------------
!
!---allocate
ALLOCATE (maskloc(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskloc',3,(/ncloc,nrloc,nz/),kndlog)

!---define
IF (PRESENT(mask)) THEN
   maskloc = mask
ELSE
   IF (lev.EQ.0) THEN
      SELECT CASE (TRIM(cnode))
         CASE ('0')
            maskloc = .TRUE.
         CASE ('C')
            maskloc(:,:,1) = nodeatc(1:nxloc,1:nyloc).GT.0
            k_310: DO k=2,nz
               maskloc(:,:,k) = maskloc(:,:,1)
            ENDDO k_310
         CASE ('U')
            maskloc = nodeatu(1:nxloc,1:nyloc,:).GT.1
         CASE ('V')
            maskloc = nodeatv(1:nxloc,1:nyloc,:).GT.1
         CASE ('UV')
            maskloc = nodeatuv(1:nxloc,1:nyloc,:).GT.0
      END SELECT
   ELSE
      SELECT CASE (TRIM(cnode))
         CASE ('0')
            maskloc = .TRUE.
         CASE ('C')
            maskloc(:,:,1) = mgvars(lev)%nodeatc(1:nxloc,1:nyloc).GT.0
         CASE ('U')
            maskloc(:,:,1) = mgvars(lev)%node2du(1:nxloc,1:nyloc).GT.1
         CASE ('V')
            maskloc(:,:,1) = mgvars(lev)%node2dv(1:nxloc,1:nyloc).GT.1
      END SELECT
      k_320: DO k=2,nz
         maskloc(:,:,k) = maskloc(:,:,1)
      ENDDO k_320
   ENDIF
ENDIF

!
!4. Serial case
!--------------
!

n3dim = SIZE(realdat,DIM=3)
IF (iopt_MPI.EQ.0) THEN
   realsum = MERGE(SUM(realdat(1:nxloc,1:nyloc,:),MASK=maskloc),0.0,&
                 & ANY(maskloc))
   GOTO 1000
ENDIF

!
!5. Local storage array
!----------------------
!
!---allocate
ALLOCATE (realloc(1-nhdims(1):nxloc+nhdims(2),1-nhdims(3):nyloc+nhdims(4),&
            & n3dim),STAT=errstat)
CALL error_alloc('realloc',3,(/nxloc+nhdims(1)+nhdims(2),&
                             & nyloc+nhdims(3)+nhdims(4),n3dim/),kndrtype)

!---define
realloc = realdat
WHERE (.NOT.maskloc)
   realloc(1:nxloc,1:nyloc,:) = float_fill
END WHERE

!
!6. Global array
!---------------
!
!---allocate
ALLOCATE (realglb(1-nhdims(1):nx+nhdims(2),1-nhdims(3):ny+nhdims(4),n3dim),&
        & STAT=errstat)
CALL error_alloc('realglb',3,(/nx+nhdims(1)+nhdims(2),&
                             & ny+nhdims(3)+nhdims(4),n3dim/),kndrtype)

!---combine
CALL combine_mod(realglb,realloc,(/1-nhdims(1),1-nhdims(3),1/),ivarid,0.0,&
               & lev,iddest,commop,itype,all)

!
!7. Global sum
!-------------
!

IF (all.OR.idloc.EQ.iddest) THEN
   IF (ANY(ABS(realglb(1:nc,1:nr,:)-float_fill).GT.float_fill_eps)) THEN
      realsum = SUM(realglb(1:nc,1:nr,:),&
                  & MASK=ABS(realglb(1:nc,1:nr,:)-float_fill).GT.float_fill_eps)
   ELSE
      realsum = 0.0
   ENDIF
ENDIF

!
!8. Deallocate arrays
!--------------------
!

DEALLOCATE (realloc,realglb)
1000 DEALLOCATE (maskloc)

CALL log_timer_out()


RETURN

END SUBROUTINE sum2_vars_real_3d

!========================================================================

SUBROUTINE sum2_vars_real_4d(realdat,realsum,nhdims,cnode,ivarid,mglevel,&
                           & idroot,comm,commtype,commall,mask)
!************************************************************************
!
! *sum2_vars_real_4d* Process independent computation of the global sum of a
!                     real 4-D model array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paral_utilities.f90  V2.8
!
! Description - result is returned on all processes if commall is PRESENT
!               and .TRUE. or to the process with rank idroot otherwise
!
! Module calls - combine_mod, error_alloc
!
!************************************************************************
!
USE grid
USE gridpars
USE multigrid
USE syspars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: combine_mod

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: commall
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: ivarid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, idroot, mglevel
REAL, INTENT(OUT) :: realsum
INTEGER, INTENT(IN), DIMENSION(4) :: nhdims
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) :: mask
REAL, INTENT(IN), DIMENSION(:,:,:,:) :: realdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL    Local model array
!*realsum*  REAL    Global sum on output
!*nhdims*   INTEGER Halo sizes (WESN directions)   
!*cnode*    CHAR    Grid node of model array ('0', 'C', 'U', 'V', 'E')
!*ivarid*   INTEGER Variable id
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!*idroot*   INTEGER id of root process
!*comm*     INTEGER Communicator for combine operations
!*commtype* INTEGER Type of communication (1/8)
!*commall*  LOGICAL Result available on all processes if .TRUE. and PRESENT
!*mask*     LOGICAL Local mask array to exclude dry points if PRESENT
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: all
INTEGER :: commop, iddest, itype, k, l, lev, nx, nxloc, ny, nyloc, n3dim, n4dim
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskloc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: realloc, realglb


procname(pglev+1) = 'sum2_vars_real_4d'
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
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
!2. Initialise parameters
!------------------------
!

IF (lev.EQ.0) THEN
   nx = nc; ny = nr
   nxloc = ncloc; nyloc = nrloc
ELSE
   nx = mgvars(lev)%nc; ny = mgvars(lev)%nr
   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc
ENDIF

!
!3. Local mask array
!-------------------
!
!---allocate
ALLOCATE (maskloc(nxloc,nyloc,nz),STAT=errstat)
CALL error_alloc('maskloc',3,(/nxloc,nyloc,nz/),kndlog)

!---define
IF (PRESENT(mask)) THEN
   maskloc = mask
ELSE
   IF (lev.EQ.0) THEN
      SELECT CASE (TRIM(cnode))
         CASE ('0')
            maskloc = .TRUE.
         CASE ('C')
            maskloc(:,:,1) = nodeatc(1:ncloc,1:nrloc).GT.0
            k_310: DO k=2,nz
               maskloc(:,:,k) = maskloc(:,:,1)
            ENDDO k_310
         CASE ('U')
            maskloc = nodeatu(1:ncloc,1:nrloc,:).GT.1
         CASE ('V')
            maskloc = nodeatv(1:ncloc,1:nrloc,:).GT.1
         CASE ('UV')
            maskloc = nodeatuv(1:ncloc,1:nrloc,:).GT.0
      END SELECT
   ELSE
      SELECT CASE (TRIM(cnode))
         CASE ('0')
            maskloc = .TRUE.
         CASE ('C')
            maskloc(:,:,1) = mgvars(lev)%nodeatc(1:nxloc,1:nyloc).GT.0
         CASE ('U')
            maskloc(:,:,1) = mgvars(lev)%node2du(1:nxloc,1:nyloc).GT.1
         CASE ('V')
            maskloc(:,:,1) = mgvars(lev)%node2dv(1:nxloc,1:nyloc).GT.1
      END SELECT
      k_320: DO k=2,nz
         maskloc(:,:,k) = maskloc(:,:,1)
      ENDDO k_320
   ENDIF
ENDIF

!
!4. Serial case
!--------------
!

n3dim = SIZE(realdat,DIM=3); n4dim = SIZE(realdat,DIM=4)
IF (iopt_MPI.EQ.0) THEN
   realsum = 0.0
   IF (ANY(maskloc)) THEN
      l_210: DO l=1,n4dim
         realsum = realsum + SUM(realdat(:,:,:,l),MASK=maskloc)
      ENDDO l_210
   ENDIF
   GOTO 1000
ENDIF

!
!5. Local storage array
!----------------------
!
!---allocate
ALLOCATE (realloc(1-nhdims(1):nxloc+nhdims(2),1-nhdims(3):nyloc+nhdims(4),&
                & n3dim,n4dim),STAT=errstat)
CALL error_alloc('realloc',4,(/nxloc+nhdims(1)+nhdims(2),&
                             & nyloc+nhdims(3)+nhdims(4),n3dim,n4dim/),&
                             & kndrtype)

!---define
realloc = realdat
l_410: DO l=1,n4dim
   WHERE (.NOT.maskloc)
      realloc(1:nxloc,1:nyloc,:,l) = float_fill
   END WHERE
ENDDO l_410

!
!6. Global array
!---------------
!
!---allocate
ALLOCATE (realglb(1-nhdims(1):nx+nhdims(2),1-nhdims(3):ny+nhdims(4),&
                & n3dim,n4dim),STAT=errstat)
CALL error_alloc('realglb',4,(/nx+nhdims(1)+nhdims(2),&
                             & ny+nhdims(3)+nhdims(4),n3dim,n4dim/),kndrtype)

!---combine
CALL combine_mod(realglb,realloc,(/1-nhdims(1),1-nhdims(3),n3dim,n4dim/),&
               & ivarid,0.0,lev,iddest,commop,itype,all)

!
!7. Global sum
!-------------
!

IF (all.OR.idloc.EQ.iddest) THEN
   realsum = 0.0
   l_610: DO l=1,n4dim
      IF (ANY(ABS(realglb(1:nc,1:nr,:,l)-float_fill).GT.float_fill_eps)) THEN
         realsum = realsum + SUM(realglb(1:nc,1:nr,:,l),&
                & MASK=ABS(realglb(1:nc,1:nr,:,l)-float_fill).GT.float_fill_eps)
      ENDIF
   ENDDO l_610
ENDIF

!
!8. Deallocate arrays
!--------------------
!

DEALLOCATE (realloc,realglb)
1000 DEALLOCATE (maskloc)

CALL log_timer_out()


RETURN

END SUBROUTINE sum2_vars_real_4d


END MODULE paral_utilities
