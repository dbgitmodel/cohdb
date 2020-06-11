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

MODULE comms_MPI
!************************************************************************
!
! *comms_MPI* MPI routine library
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.10.1
!
! $Date: 2016-10-10 16:36:15 +0200 (Mon, 10 Oct 2016) $
!
! $Revision: 973 $
!
! Description -
!
! Routines - comms_allgather_int, comms_allgather_log,
!   comms_allgather_real, comms_barrier, comms_bcast_char, comms_bcast_int,
!   comms_bcast_log, comms_bcast_real, comms_comm_free, comms_comm_rank,
!   comms_comm_size, comms_comm_split, comms_dims_create, comms_errhandler,
!   comms_finalize, comms_initialise, comms_irecv_char, comms_irecv_int,
!   comms_irecv_log, comms_irecv_real, comms_isend_char, comms_isend_int,
!   comms_isend_log, comms_isend_real, comms_recv_char, comms_recv_int,
!   comms_recv_log, comms_recv_real, comms_reduce_int, comms_reduce_2int,
!   comms_reduce_log, comms_reduce_real, comms_reduce_2real, comms_send_char,
!   comms_send_int, comms_send_log, comms_send_real, conns_sendrecv_char,
!   comms_sendrecv_int, comms_sendrecv_log, comms_sendrecv_real,
!   comms_waitall, error_MPI, error_MPI_comm
!
! Module calls - inquire_var
!
!************************************************************************
!
USE iopars
USE paralpars
USE syspars
USE modvars_routines, ONLY: inquire_var
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

#ifdef MPI
INCLUDE 'mpif.h'
#define MPI_max MPI_MAX
#define MPI_min MPI_MIN
#define MPI_sum MPI_SUM
#define MPI_prod MPI_PROD
#define MPI_land MPI_LAND
#define MPI_band MPI_BAND
#define MPI_lor MPI_LOR
#define MPI_bor MPI_BOR
#define MPI_lxor MPI_LXOR
#define MPI_maxloc MPI_MAXLOC
#define MPI_minloc MPI_MINLOC
!
!*Local variables
!
CHARACTER (LEN=lenname) :: f90_name
INTEGER, DIMENSION(MPI_status_size) :: status
INTEGER, DIMENSION(12) :: ops_reduce = &
       & (/MPI_max,MPI_min,MPI_sum,MPI_prod,MPI_land,MPI_band,&
         & MPI_lor,MPI_bor,MPI_lxor,MPI_bxor,MPI_maxloc,MPI_minloc/)
#endif /*MPI*/
LOGICAL :: info = .FALSE.
INTEGER :: npcc


CONTAINS

!========================================================================

SUBROUTINE comms_allgather_int(sendbuf,count,nprocs,recvbuf,comm,ivarid)
!************************************************************************
!
! *comms_allgather_int* Gather integer data on all processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_allgather
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, nprocs
INTEGER, INTENT(IN), DIMENSION(count) :: sendbuf
INTEGER, INTENT(OUT), DIMENSION(count,nprocs) :: recvbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*sendbuf* INTEGER Starting address of send buffer
!*count*   INTEGER Number of data
!*nprocs*  INTEGER Number of processes
!*recvbuf* INTEGER Starting address of receive buffer
!*comm*    INTEGER Communicator
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_allgather_int'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_allgather(sendbuf,count,integer_MPI,recvbuf,count,integer_MPI,comm,&
                 & errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_allgather')
#else
recvbuf = 0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_allgather_int

!========================================================================

SUBROUTINE comms_allgather_log(sendbuf,count,nprocs,recvbuf,comm,ivarid)
!************************************************************************
!
! *comms_allgather_log* Gather logical data on all processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_allgather
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, nprocs
LOGICAL, INTENT(IN), DIMENSION(count) :: sendbuf
LOGICAL, INTENT(OUT), DIMENSION(count,nprocs) :: recvbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*sendbuf* LOGICAL Starting address of send buffer
!*count*   INTEGER Number of data
!*nprocs*  INTEGER Number of processes
!*recvbuf* LOGICAL Starting address of receive buffer
!*comm*    INTEGER Communicator
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_allgather_log'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_allgather(sendbuf,count,logical_MPI,recvbuf,count,logical_MPI,comm,&
                 & errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_allgather')
#else
recvbuf = .FALSE.
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_allgather_log

!========================================================================

SUBROUTINE comms_allgather_real(sendbuf,count,nprocs,recvbuf,comm,ivarid)
!************************************************************************
!
! *comms_allgather_real* Gather real data on all processes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_allgather
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, nprocs
REAL, INTENT(IN), DIMENSION(count) :: sendbuf
REAL, INTENT(OUT), DIMENSION(count,nprocs) :: recvbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*sendbuf* REAL    Starting address of send buffer
!*count*   INTEGER Number of data
!*nprocs*  INTEGER Number of processes
!*recvbuf* REAL    Starting address of receive buffer
!*comm*    INTEGER Communicator
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_allgather_real'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_allgather(sendbuf,count,float_MPI,recvbuf,count,float_MPI,comm,&
                 & errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_allgather')
#else
recvbuf = 0.0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_allgather_real

!========================================================================

SUBROUTINE comms_barrier(comm)
!************************************************************************
!
! *comms_barrier* Synchronise processes in a communicator
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_barrier
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*comm*    INTEGER Communicator
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'comms_barrier'
CALL log_timer_in(npcc,info=info)

#ifdef MPI
CALL MPI_barrier(comm,errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_barrier',abort=.TRUE.)
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_barrier

!========================================================================

SUBROUTINE comms_bcast_char(cbuf,count,root,comm,ivarid)
!************************************************************************
!
! *comms_bcast_char* Broadcast character data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_bcast
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, root
CHARACTER, INTENT(INOUT), DIMENSION(count) :: cbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*cbuf*    CHAR    Starting address of buffer
!*count*   INTEGER Number of data
!*root*    INTEGER Rank of broadcast root
!*comm*    INTEGER Communicator
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_bcast_char'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_bcast(cbuf,count,character_MPI,root,comm,errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_bcast')
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_bcast_char

!========================================================================

SUBROUTINE comms_bcast_int(ibuf,count,root,comm,ivarid)
!************************************************************************
!
! *comms_bcast_int* Broadcast integer data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_bcast
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, root
INTEGER, INTENT(INOUT), DIMENSION(count) :: ibuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ibuf*    INTEGER Starting address of buffer
!*count*   INTEGER Number of data
!*root*    INTEGER Rank of broadcast root
!*comm*    INTEGER Communicator
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_bcast_int'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_bcast(ibuf,count,integer_MPI,root,comm,errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_bcast')
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_bcast_int

!========================================================================

SUBROUTINE comms_bcast_log(lbuf,count,root,comm,ivarid)
!************************************************************************
!
! *comms_bcast_log* Broadcast logical data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_bcast
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, root
LOGICAL, INTENT(INOUT), DIMENSION(count) :: lbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*lbuf*    LOGICAL Starting address of buffer
!*count*   INTEGER Number of data
!*root*    INTEGER Rank of broadcast root
!*comm*    INTEGER Communicator
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_bcast_log'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_bcast(lbuf,count,logical_MPI,root,comm,errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_bcast')
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_bcast_log

!========================================================================

SUBROUTINE comms_bcast_real(rbuf,count,root,comm,ivarid)
!************************************************************************
!
! *comms_bcast_real* Broadcast real data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_bcast
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, root
REAL, INTENT(INOUT), DIMENSION(count) :: rbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*rbuf*    REAL    Starting address of buffer
!*count*   INTEGER Number of data
!*root*    INTEGER Rank of broadcast root
!*comm*    INTEGER Communicator
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_bcast_real'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_bcast(rbuf,count,float_MPI,root,comm,errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_bcast')
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_bcast_real

!========================================================================

SUBROUTINE comms_comm_free(comm)
!************************************************************************
!
! *comms_comm_free* Deallocate communicator
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_comm_free
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(INOUT) :: comm

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*comm*    INTEGER Communicator to be deallocated
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'comms_comm_free'
CALL log_timer_in(npcc,info=info)

#ifdef MPI
CALL MPI_comm_free(comm,errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_comm_free')
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_comm_free

!========================================================================

SUBROUTINE comms_comm_rank(comm,rank)
!************************************************************************
!
! *comms_comm_rank* Return rank of calling process
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_comm_rank
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm
INTEGER, INTENT(OUT) :: rank

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*comm*    INTEGER Communicator
!*rank*    INTEGER Rank of calling process
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'comms_comm_rank'
CALL log_timer_in(npcc,info=info)

#ifdef MPI
CALL MPI_comm_rank(comm,rank,errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_comm_rank')
#else
rank = 0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)

RETURN

END SUBROUTINE comms_comm_rank

!========================================================================

SUBROUTINE comms_comm_size(comm,size)
!************************************************************************
!
! *comms_comm_size* Number of processes in the group of comm
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_comm_size
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm
INTEGER, INTENT(OUT) :: size

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*comm*    INTEGER Communicator
!*size*    INTEGER Number of processes in communicator
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'comms_comm_size'
CALL log_timer_in(npcc,info=info)

#ifdef MPI
CALL MPI_comm_size(comm,size,errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_comm_size')
#else
size = 0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_comm_size

!========================================================================

SUBROUTINE comms_comm_split(comm,color,key,newcomm)
!************************************************************************
!
! *comms_comm_split* Create a new communicator by partinioning of comm
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_comm_split
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: color, comm, key
INTEGER, INTENT(OUT) :: newcomm

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*comm*    INTEGER Communicator
!*color*   INTEGER Control of subset assignment
!*key*     INTEGER Control of rank assignment
!*newcomm* INTEGER New communicator
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'comms_comm_split'
CALL log_timer_in(npcc,info=info)

#ifdef MPI
CALL MPI_comm_split(comm,color,key,newcomm,errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_comm_split')
#else
newcomm = 0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_comm_split

!========================================================================

SUBROUTINE comms_dims_create(nnodes,ndims,dims)
!************************************************************************
!
! *comms_dims_create* Returns in 'dims' a Cartesian grid of dimension 'ndims'
!                     and a total of 'nnodes' nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_dims_create
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: ndims, nnodes
INTEGER, INTENT(INOUT), DIMENSION(ndims) :: dims

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*nnodes*  INTEGER Number of nodes in a grid
!*ndims*   INTEGER Number of grid dimensions
!*dims*    INTEGER Number of nodes in each dimension
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'comms_dimms_create'
CALL log_timer_in(npcc,info=info)

#ifdef MPI
CALL MPI_dims_create(nnodes,ndims,dims,errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_dims_create')
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_dims_create

!========================================================================

SUBROUTINE comms_errhandler(comm)
!************************************************************************
!
! *comms_errhandler* Set error handler to MPI_errors_return
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description - reset error handler if necessary
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_errhandler_free, MPI_errhandler_get, MPI_errhandler_set
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*comm*    INTEGER Communicator
!
!------------------------------------------------------------------------------
!
!*Local variables
!
#ifdef MPI
INTEGER :: errhandler


CALL MPI_errhandler_get(comm,errhandler,errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_errhandler_get')
IF (errhandler.NE.MPI_errors_return) THEN
   CALL MPI_errhandler_free(errhandler,errstat)
   IF (errstat.NE.MPI_success) CALL error_MPI('MPI_errhandler_free')
   CALL MPI_errhandler_set(comm,MPI_errors_return,errstat)
   IF (errstat.NE.MPI_success) CALL error_MPI('MPI_errhandler_set')
ENDIF

#endif /*MPI*/


RETURN

END SUBROUTINE comms_errhandler

!========================================================================

SUBROUTINE comms_finalize
!************************************************************************
!
! *comms_finalize* Finalise MPI
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description - verify first whether all pending communications are completed
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_finalize
!
!************************************************************************
!


#ifdef MPI
CALL MPI_finalize(errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_finalize')
#endif /*MPI*/


RETURN

END SUBROUTINE comms_finalize

!========================================================================

SUBROUTINE comms_initialise
!************************************************************************
!
! *comms_initialise* Initialise MPI
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_init, MPI_initialized
!
!************************************************************************
!
!*Local variables
!
#ifdef MPI
LOGICAL :: flag


CALL MPI_initialized(flag,errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_initialized')
IF (.NOT.flag) THEN
   CALL MPI_init(errstat)
   IF (errstat.NE.MPI_success) CALL error_MPI('MPI_init')
ENDIF

#endif /*MPI*/


RETURN

END SUBROUTINE comms_initialise

!========================================================================

SUBROUTINE comms_irecv_char(cbuf,count,source,tag,comm,request,ivarid)
!************************************************************************
!
! *comms_irecv_char* Non-blocking receive of character data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_irecv
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, source, tag
INTEGER, INTENT(OUT) :: request
CHARACTER, INTENT(OUT), DIMENSION(count) :: cbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*cbuf*    CHAR    Starting address of receive buffer
!*count*   INTEGER Number of data
!*source*  INTEGER Rank of source process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*request* INTEGER Communication request
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN

procname(pglev+1) = 'comms_irecv_char'
CALL log_timer_in(npcc,info=info)

#ifdef MPI
CALL MPI_irecv(cbuf,count,character_MPI,source,tag,comm,request,errstat)
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to receive character data: '//&
                            & TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to receive character data'
      ENDIF
   ENDIF
   CALL error_MPI_comm('MPI_irecv',source,idloc)
ENDIF
#else
request = 0
cbuf = ''
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_irecv_char

!========================================================================

SUBROUTINE comms_irecv_int(ibuf,count,source,tag,comm,request,ivarid)
!************************************************************************
!
! *comms_irecv_int* Non-blocking receive of integer data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_irecv
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, source, tag
INTEGER, INTENT(OUT) :: request
INTEGER, INTENT(OUT), DIMENSION(count) :: ibuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ibuf*    INTEGER Starting address of receive buffer
!*count*   INTEGER Number of data
!*source*  INTEGER Rank of source process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*request* INTEGER Communication request
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_irecv_int'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_irecv(ibuf,count,integer_MPI,source,tag,comm,request,errstat)
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to receive integer data: '//TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to receive integer data'
      ENDIF
   ENDIF
   CALL error_MPI_comm('MPI_irecv',source,idloc)
ENDIF
#else
request = 0
ibuf = 0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_irecv_int

!========================================================================

SUBROUTINE comms_irecv_log(lbuf,count,source,tag,comm,request,ivarid)
!************************************************************************
!
! *comms_irecv_log* Non-blocking receive of logical data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_irecv
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, source, tag
INTEGER, INTENT(OUT) :: request
LOGICAL, INTENT(OUT), DIMENSION(count) :: lbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*lbuf*    LOGICAL Starting address of receive buffer
!*count*   INTEGER Number of data
!*source*  INTEGER Rank of source process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*request* INTEGER Communication request
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_irecv_log'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_irecv(lbuf,count,logical_MPI,source,tag,comm,request,errstat)
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to receive logical data: '//TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to receive logical data'
      ENDIF
   ENDIF
   CALL error_MPI_comm('MPI_irecv',source,idloc)
ENDIF
#else
request = 0
lbuf = .FALSE.
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_irecv_log

!========================================================================

SUBROUTINE comms_irecv_real(rbuf,count,source,tag,comm,request,ivarid)
!************************************************************************
!
! *comms_irecv_real* Non-blocking receive of real data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_irecv
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, source, tag
INTEGER, INTENT(OUT) :: request
REAL, INTENT(OUT), DIMENSION(count) :: rbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*rbuf*    REAL    Starting address of receive buffer
!*count*   INTEGER Number of data
!*source*  INTEGER Rank of source process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*request* INTEGER Communication request
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_irecv_real'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_irecv(rbuf,count,float_MPI,source,tag,comm,request,errstat)
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to receive real data: '//TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to receive real data'
      ENDIF
   ENDIF
   CALL error_MPI_comm('MPI_irecv',source,idloc)
ENDIF
#else
request = 0
rbuf = 0.0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_irecv_real

!========================================================================

SUBROUTINE comms_isend_char(cbuf,count,dest,tag,comm,mode,request,ivarid)
!************************************************************************
!
! *comms_isend_char* Non-blocking send of character data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_isend, MPI_issend
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, dest, ivarid, mode, tag
INTEGER, INTENT(OUT) :: request
CHARACTER, INTENT(IN), DIMENSION(count) :: cbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*cbuf*    CHAR    Starting address of send buffer
!*count*   INTEGER Number of data
!*dest*    INTEGER Rank of destination process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*request* INTEGER Communication request
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
#ifdef MPI
CHARACTER (len=lenname) :: mpiname
#endif /*MPI*/

IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_isend_char'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
SELECT CASE (mode)
CASE (1)
   CALL MPI_isend(cbuf,count,character_MPI,dest,tag,comm,request,errstat)
   mpiname = 'MPI_isend'
CASE (2)   
   CALL MPI_issend(cbuf,count,character_MPI,dest,tag,comm,request,errstat)
   mpiname = 'MPI_issend'
END SELECT
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to send character data: '//TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to send character data'
      ENDIF
   ENDIF
   CALL error_MPI_comm(mpiname,idloc,dest)
ENDIF
#else
request = 0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_isend_char

!========================================================================

SUBROUTINE comms_isend_int(ibuf,count,dest,tag,comm,mode,request,ivarid)
!************************************************************************
!
! *comms_isend_int* Non-blocking send of integer data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_isend, MPI_issend
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, dest, ivarid, mode, tag
INTEGER, INTENT(OUT) :: request
INTEGER, INTENT(IN), DIMENSION(count) :: ibuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ibuf*    INTEGER Starting address of send buffer
!*count*   INTEGER Number of data
!*dest*    INTEGER Rank of destination process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*request* INTEGER Communication request
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
#ifdef MPI
CHARACTER (len=lenname) :: mpiname
#endif /*MPI*/


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_isend_int'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
SELECT CASE (mode)
CASE (1)
   CALL MPI_isend(ibuf,count,integer_MPI,dest,tag,comm,request,errstat)
   mpiname = 'MPI_isend'
CASE (2)   
   CALL MPI_issend(ibuf,count,integer_MPI,dest,tag,comm,request,errstat)
   mpiname = 'MPI_issend'
END SELECT
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to send integer data: '//TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to send integer data'
      ENDIF
   ENDIF
   CALL error_MPI_comm(mpiname,idloc,dest)
ENDIF
#else
request = 0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_isend_int

!========================================================================

SUBROUTINE comms_isend_log(lbuf,count,dest,tag,comm,mode,request,ivarid)
!************************************************************************
!
! *comms_isend_log* Non-blocking send of logical data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_isend, MPI_issend
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, dest, ivarid, mode, tag
INTEGER, INTENT(OUT) :: request
LOGICAL, INTENT(IN), DIMENSION(count) :: lbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*lbuf*    LOGICAL Starting address of send buffer
!*count*   INTEGER Number of data
!*dest*    INTEGER Rank of destination process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*request* INTEGER Communication request
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
#ifdef MPI
CHARACTER (len=lenname) :: mpiname
#endif /*MPI*/


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_isend_log'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
SELECT CASE (mode)
CASE (1)
   CALL MPI_isend(lbuf,count,logical_MPI,dest,tag,comm,request,errstat)
   mpiname = 'MPI_isend'
CASE (2)   
   CALL MPI_issend(lbuf,count,logical_MPI,dest,tag,comm,request,errstat)
   mpiname = 'MPI_issend'
END SELECT
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to send logical data: '//TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to send logical data'
      ENDIF
   ENDIF
   CALL error_MPI_comm(mpiname,idloc,dest)
ENDIF
#else
request = 0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_isend_log

!========================================================================

SUBROUTINE comms_isend_real(rbuf,count,dest,tag,comm,mode,request,ivarid)
!************************************************************************
!
! *comms_isend_real* Non-blocking send of real data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_isend, MPI_issend
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, dest, ivarid, mode, tag
INTEGER, INTENT(OUT) :: request
REAL, INTENT(IN), DIMENSION(count) :: rbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*rbuf*    REAL    Starting address of send buffer
!*count*   INTEGER Number of data
!*dest*    INTEGER Rank of destination process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*request* INTEGER Communication request
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
#ifdef MPI
CHARACTER (len=lenname) :: mpiname
#endif /*MPI*/


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_isend_real'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
SELECT CASE (mode)
CASE (1)
   CALL MPI_isend(rbuf,count,float_MPI,dest,tag,comm,request,errstat)
   mpiname = 'MPI_isend'
CASE (2)   
   CALL MPI_issend(rbuf,count,float_MPI,dest,tag,comm,request,errstat)
   mpiname = 'MPI_issend'
END SELECT
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to send real data: '//TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to send real data'
      ENDIF
   ENDIF
   CALL error_MPI_comm(mpiname,idloc,dest)
ENDIF
#else
request = 0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_isend_real

!========================================================================

SUBROUTINE comms_recv_char(cbuf,count,source,tag,comm,ivarid)
!************************************************************************
!
! *comms_recv_char* Blocking receive of character data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_recv
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, source, tag
CHARACTER, INTENT(OUT), DIMENSION(count) :: cbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*cbuf*    CHAR    Starting address of receive buffer
!*count*   INTEGER Number of data
!*source*  INTEGER Rank of source process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_recv_char'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_recv(cbuf,count,character_MPI,source,tag,comm,status,errstat)
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to receive character data: '//&
                           & TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to receive character data'
      ENDIF
   ENDIF
   CALL error_MPI_comm('MPI_recv',source,idloc)
ENDIF
#else
cbuf = ''
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_recv_char

!========================================================================

SUBROUTINE comms_recv_int(ibuf,count,source,tag,comm,ivarid)
!************************************************************************
!
! *comms_recv_int* Blocking receive of integer data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_recv
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, source, tag
INTEGER, INTENT(OUT), DIMENSION(count) :: ibuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ibuf*    INTEGER Starting address of receive buffer
!*count*   INTEGER Number of data
!*source*  INTEGER Rank of source process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_recv_int'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_recv(ibuf,count,integer_MPI,source,tag,comm,status,errstat)
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to receive integer data: '//TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to receive integer data'
      ENDIF
   ENDIF
   CALL error_MPI_comm('MPI_recv',source,idloc)
ENDIF
#else
ibuf = 0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_recv_int

!========================================================================

SUBROUTINE comms_recv_log(lbuf,count,source,tag,comm,ivarid)
!************************************************************************
!
! *comms_recv_log* Blocking receive of logical data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_recv
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, source, tag
LOGICAL, INTENT(OUT), DIMENSION(count) :: lbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*lbuf*    LOGICAL Starting address of receive buffer
!*count*   INTEGER Number of data
!*source*  INTEGER Rank of source process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_recv_log'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_recv(lbuf,count,logical_MPI,source,tag,comm,status,errstat)
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to receive logical data: '//TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to receive logical data'
      ENDIF
   ENDIF
   CALL error_MPI_comm('MPI_recv',source,idloc)
ENDIF
#else
lbuf = .FALSE.
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_recv_log

!========================================================================

SUBROUTINE comms_recv_real(rbuf,count,source,tag,comm,ivarid)
!************************************************************************
!
! *comms_recv_real* Blocking receive of real data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_recv
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, ivarid, source, tag
REAL, INTENT(OUT), DIMENSION(count) :: rbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*rbuf*    REAL    Starting address of receive buffer
!*count*   INTEGER Number of data
!*source*  INTEGER Rank of source process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_recv_real'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_recv(rbuf,count,float_MPI,source,tag,comm,status,errstat)
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to receive real data: '//TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to receive real data'
      ENDIF
   ENDIF
   CALL error_MPI_comm('MPI_recv',source,idloc)
ENDIF
#else
rbuf = 0.0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_recv_real

!========================================================================

SUBROUTINE comms_reduce_int(sendbuf,recvbuf,iop,comm,root,ivarid)
!************************************************************************
!
! *comms_reduce_int* Global reduce operation on integer data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description - result is returned on root or on all processes
!
! Reference -
!
! Module calls - error_abort, error_MPI, error_vals_var_int
!
! MPI calls - MPI_allreduce, MPI_reduce
!
!************************************************************************
!
USE error_routines, ONLY: error_abort, error_vals_var_int

!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, ivarid, iop
INTEGER, INTENT(IN), OPTIONAL :: root
INTEGER, INTENT(IN) :: sendbuf
INTEGER, INTENT(OUT) :: recvbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*sendbuf* INTEGER Send buffer
!*recvbuf* INTEGER Receive buffer
!*iop*     INTEGER Reduce operation
!*comm*    INTEGER Communicator
!*root*    INTEGER Rank of root process
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'comms_reduce_int'
CALL log_timer_in(npcc,ivarid,info=info)

SELECT CASE (iop)
CASE (1:4,6,8,10)
#ifdef MPI
IF (PRESENT(root)) THEN
   CALL MPI_reduce(sendbuf,recvbuf,1,integer_MPI,ops_reduce(iop),root,comm,&
                 & errstat)
   IF (errstat.NE.MPI_success) CALL error_MPI('MPI_reduce')
ELSE
   CALL MPI_allreduce(sendbuf,recvbuf,1,integer_MPI,ops_reduce(iop),comm,&
                    & errstat)
   IF (errstat.NE.MPI_success) CALL error_MPI('MPI_allreduce')
ENDIF
#else
recvbuf = 0
#endif /*MPI*/
CASE DEFAULT
   CALL error_vals_var_int(iop,'iop',7,(/1,2,3,4,6,8,10/))
   CALL error_abort('comms_reduce_int',ierrno_arg)
END SELECT

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_reduce_int

!========================================================================

SUBROUTINE comms_reduce_2int(sendbuf,recvbuf,iop,comm,root,ivarid)
!************************************************************************
!
! *comms_reduce_2int* Global reduce operation on pairs of integer data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Reference -
!
! Module calls - error_abort, error_MPI, error_vals_var_int
!
! MPI calls - MPI_allreduce, MPI_reduce
!
!************************************************************************
!
USE error_routines, ONLY: error_abort, error_vals_var_int

!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, ivarid, iop
INTEGER, INTENT(IN), OPTIONAL :: root
INTEGER, INTENT(IN), DIMENSION(2) :: sendbuf
INTEGER, INTENT(OUT), DIMENSION(2) :: recvbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*sendbuf* INTEGER Send buffer
!*recvbuf* INTEGER Receive buffer
!*iop*     INTEGER Reduce operation
!*comm*    INTEGER Communicator
!*root*    INTEGER Rank of root process
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'comms_reduce_2int'
CALL log_timer_in(npcc,ivarid,info=info)

SELECT CASE (iop)
CASE (11,12)
#ifdef MPI
IF (PRESENT(root)) THEN
   CALL MPI_reduce(sendbuf,recvbuf,1,MPI_2integer,ops_reduce(iop),root,comm,&
                 & errstat)
   IF (errstat.NE.MPI_success) CALL error_MPI('MPI_reduce')
ELSE
   CALL MPI_allreduce(sendbuf,recvbuf,1,MPI_2integer,ops_reduce(iop),comm,&
                    & errstat)
   IF (errstat.NE.MPI_success) CALL error_MPI('MPI_allreduce')
ENDIF
#else
recvbuf = 0
#endif /*MPI*/
CASE DEFAULT
   CALL error_vals_var_int(iop,'iop',2,(/11,12/))
   CALL error_abort('comms_reduce_2int',ierrno_arg)
END SELECT

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_reduce_2int

!========================================================================

SUBROUTINE comms_reduce_log(sendbuf,recvbuf,iop,comm,root,ivarid)
!************************************************************************
!
! *comms_reduce_log* Global reduce operation on logical data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description - result is returned on root or on all processes
!
! Reference -
!
! Module calls - error_abort, error_MPI, error_vals_var_int
!
! MPI calls - MPI_allreduce, MPI_reduce
!
!************************************************************************
!
USE error_routines, ONLY: error_abort, error_vals_var_int

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: sendbuf
LOGICAL, INTENT(OUT) :: recvbuf
INTEGER, INTENT(IN) :: comm, ivarid, iop
INTEGER, INTENT(IN), OPTIONAL :: root

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*sendbuf* LOGICAL Send buffer
!*recvbuf* LOGICAL Receive buffer
!*iop*     INTEGER Reduce operation
!*comm*    INTEGER Communicator
!*root*    INTEGER Rank of root process
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'comms_reduce_log'
CALL log_timer_in(npcc,ivarid,info=info)

SELECT CASE (iop)
CASE (5,7,9)
#ifdef MPI
IF (PRESENT(root)) THEN
   CALL MPI_reduce(sendbuf,recvbuf,1,logical_MPI,ops_reduce(iop),root,comm,&
                 & errstat)
   IF (errstat.NE.MPI_success) CALL error_MPI('MPI_reduce')
ELSE
   CALL MPI_allreduce(sendbuf,recvbuf,1,logical_MPI,ops_reduce(iop),comm,&
                    & errstat)
   IF (errstat.NE.MPI_success) CALL error_MPI('MPI_allreduce')
ENDIF
#else
recvbuf = .FALSE.
#endif /*MPI*/
CASE DEFAULT
   CALL error_vals_var_int(iop,'iop',3,(/5,7,9/))
   CALL error_abort('comms_reduce_log',ierrno_arg)
END SELECT

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_reduce_log

!========================================================================

SUBROUTINE comms_reduce_real(sendbuf,recvbuf,iop,comm,root,ivarid)
!************************************************************************
!
! *comms_reduce_real* Global reduce operation on real data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description - result is returned on root or on all processes
!
! Reference -
!
! Module calls - error_abort, error_MPI, error_vals_var_int
!
! MPI calls - MPI_allreduce, MPI_reduce
!
!************************************************************************
!
USE error_routines, ONLY: error_abort, error_vals_var_int

!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, ivarid,iop
INTEGER, INTENT(IN), OPTIONAL :: root
REAL, INTENT(IN) :: sendbuf
REAL, INTENT(OUT) :: recvbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*sendbuf* REAL    Send buffer
!*recvbuf* REAL    Receive buffer
!*iop*     INTEGER Reduce operation
!*comm*    INTEGER Communicator
!*root*    INTEGER Rank of root process
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'comms_reduce_real'
CALL log_timer_in(npcc,ivarid,info=info)

SELECT CASE (iop)
CASE (1:4)
#ifdef MPI
IF (PRESENT(root)) THEN
   CALL MPI_reduce(sendbuf,recvbuf,1,float_MPI,ops_reduce(iop),root,comm,&
                 & errstat)
   IF (errstat.NE.MPI_success) CALL error_MPI('MPI_reduce')
ELSE
   CALL MPI_allreduce(sendbuf,recvbuf,1,float_MPI,ops_reduce(iop),comm,errstat)
   IF (errstat.NE.MPI_success) CALL error_MPI('MPI_allreduce')
ENDIF
#else
recvbuf = 0.0
#endif /*MPI*/
CASE DEFAULT
   CALL error_vals_var_int(iop,'iop',4,(/1,2,3,4/))
   CALL error_abort('comms_reduce_real',ierrno_arg)
END SELECT

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_reduce_real

!========================================================================

SUBROUTINE comms_reduce_2real(sendbuf,recvbuf,iop,comm,root,ivarid)
!************************************************************************
!
! *comms_reduce_2real* Global reduce operation on pairs of real data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description - result is returned on root or on all processes
!
! Reference -
!
! Module calls - error_abort, error_MPI, error_vals_var_int
!
! MPI calls - MPI_allreduce, MPI_reduce
!
!************************************************************************
!
USE error_routines, ONLY: error_abort, error_vals_var_int

!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, iop, ivarid
INTEGER, INTENT(IN), OPTIONAL :: root
REAL, INTENT(IN), DIMENSION(2) :: sendbuf
REAL, INTENT(OUT), DIMENSION(2) :: recvbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*sendbuf* REAL    Send buffer
!*recvbuf* REAL    Receive buffer
!*iop*     INTEGER Reduce operation
!*comm*    INTEGER Communicator
!*root*    INTEGER Rank of root process
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'comms_reduce_2real'
CALL log_timer_in(npcc,ivarid,info=info)

SELECT CASE (iop)
CASE (11,12)
#ifdef MPI
IF (PRESENT(root)) THEN
   CALL MPI_reduce(sendbuf,recvbuf,1,float2_MPI,ops_reduce(iop),root,comm,&
                 & errstat)
   IF (errstat.NE.MPI_success) CALL error_MPI('MPI_reduce')
ELSE
   CALL MPI_allreduce(sendbuf,recvbuf,1,float2_MPI,ops_reduce(iop),comm,errstat)
   IF (errstat.NE.MPI_success) CALL error_MPI('MPI_allreduce')
ENDIF
#else
recvbuf = 0.0
#endif /*MPI*/
CASE DEFAULT
   CALL error_vals_var_int(iop,'iop',2,(/11,12/))
   CALL error_abort('comms_reduce_2real',ierrno_arg)
END SELECT

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_reduce_2real

!========================================================================

SUBROUTINE comms_send_char(cbuf,count,dest,tag,comm,mode,ivarid)
!************************************************************************
!
! *comms_send_char* Blocking send of character data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_send, MPI_ssend
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, dest, ivarid, mode, tag
CHARACTER, INTENT(IN), DIMENSION(count) :: cbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*cbuf*    CHAR    Starting address of send buffer
!*count*   INTEGER Number of data
!*dest*    INTEGER Rank of destination process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*mode*    INTEGER Type of send operation
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
#ifdef MPI
CHARACTER (len=lenname) :: mpiname
#endif /*MPI*/


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_send_char'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
SELECT CASE (mode)
CASE (1)
   CALL MPI_send(cbuf,count,character_MPI,dest,tag,comm,errstat)
   mpiname = 'MPI_send'
CASE (2)   
   CALL MPI_ssend(cbuf,count,character_MPI,dest,tag,comm,errstat)
   mpiname = 'MPI_ssend'
END SELECT
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to send character data: '//TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to send character data'
      ENDIF
   ENDIF
   CALL error_MPI_comm(mpiname,idloc,dest)
ENDIF
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_send_char

!========================================================================

SUBROUTINE comms_send_int(ibuf,count,dest,tag,comm,mode,ivarid)
!************************************************************************
!
! *comms_send_int* Blocking send of integer data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_send, MPI_ssend
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, dest, ivarid, mode, tag
INTEGER, INTENT(IN), DIMENSION(count) :: ibuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ibuf*    INTEGER Starting address of send buffer
!*count*   INTEGER Number of data
!*dest*    INTEGER Rank of destination process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*mode*    INTEGER Type of send operation
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
#ifdef MPI
CHARACTER (len=lenname) :: mpiname
#endif /*MPI*/

IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_send_int'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
SELECT CASE (mode)
CASE (1)
   CALL MPI_send(ibuf,count,integer_MPI,dest,tag,comm,errstat)
   mpiname = 'MPI_send'
CASE (2)   
   CALL MPI_ssend(ibuf,count,integer_MPI,dest,tag,comm,errstat)
   mpiname = 'MPI_ssend'
END SELECT
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to send integer data: '//TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to send integer data'
      ENDIF
   ENDIF
   CALL error_MPI_comm(mpiname,idloc,dest)
ENDIF
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_send_int

!========================================================================

SUBROUTINE comms_send_log(lbuf,count,dest,tag,comm,mode,ivarid)
!************************************************************************
!
! *comms_send_log* Blocking send of logical data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_send, MPI_ssend
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, dest, ivarid, mode, tag
LOGICAL, INTENT(IN), DIMENSION(count) :: lbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*lbuf*    LOGICAL Starting address of send buffer
!*count*   INTEGER Number of data
!*dest*    INTEGER Rank of destination process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*mode*    INTEGER Type of send operation
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
#ifdef MPI
CHARACTER (len=lenname) :: mpiname
#endif /*MPI*/


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_send_log'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
SELECT CASE (mode)
CASE (1)
   CALL MPI_send(lbuf,count,logical_MPI,dest,tag,comm,errstat)
   mpiname = 'MPI_send'
CASE (2)   
   CALL MPI_ssend(lbuf,count,logical_MPI,dest,tag,comm,errstat)
   mpiname = 'MPI_ssend'
END SELECT
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to send logical data: '//TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to send logical data'
      ENDIF
   ENDIF
   CALL error_MPI_comm(mpiname,idloc,dest)
ENDIF
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_send_log

!========================================================================

SUBROUTINE comms_send_real(rbuf,count,dest,tag,comm,mode,ivarid)
!************************************************************************
!
! *comms_send_real* Blocking send of real data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_send, MPI_ssend
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, count, dest, ivarid, mode, tag
REAL, INTENT(IN), DIMENSION(count) :: rbuf

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*rbuf*    REAL    Starting address of send buffer
!*count*   INTEGER Number of data
!*dest*    INTEGER Rank of destination process
!*tag*     INTEGER Message tag
!*comm*    INTEGER Communicator
!*mode*    INTEGER Type of send operation
!*ivarid*  INTEGER Variable id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
#ifdef MPI
CHARACTER (len=lenname) :: mpiname
#endif /*MPI*/


IF (count.EQ.0) RETURN
procname(pglev+1) = 'comms_send_real'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
SELECT CASE (mode)
CASE (1)
   CALL MPI_send(rbuf,count,float_MPI,dest,tag,comm,errstat)
   mpiname = 'MPI_send'
CASE (2)   
   CALL MPI_ssend(rbuf,count,float_MPI,dest,tag,comm,errstat)
   mpiname = 'MPI_ssend'
END SELECT
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to send real data: '//TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to send real data'
      ENDIF
   ENDIF
   CALL error_MPI_comm(mpiname,idloc,dest)
ENDIF
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_send_real

!========================================================================

SUBROUTINE comms_sendrecv_char(sendbuf,sendcount,dest,sendtag,recvbuf,&
                             & recvcount,source,recvtag,comm,ivarid)
!************************************************************************
!
! *comms_sendrecv_char*  Send-receive of character data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_sendrecv
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, dest, ivarid, recvcount, recvtag, sendcount, &
                     & sendtag, source
CHARACTER, INTENT(IN), DIMENSION(sendcount) :: sendbuf
CHARACTER, INTENT(OUT), DIMENSION(recvcount) :: recvbuf

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*sendbuf*   CHAR    Starting address of send buffer
!*sendcount* INTEGER Number of data in send buffer
!*dest*      INTEGER Rank of destination process
!*sendtag*   INTEGER Message tag of send operation
!*recvbuf*   CHAR    Starting address of receive buffer
!*recvcount* INTEGER Number of data in receive buffer
!*source*    INTEGER Rank of source process
!*recvtag*   INTEGER Message tag of receive operation
!*comm*      INTEGER Communicator
!*ivarid*    INTEGER Variable id
!
!------------------------------------------------------------------------------
!

IF (sendcount.EQ.0.AND.recvcount.EQ.0) RETURN
procname(pglev+1) = 'comms_sendrecv_char'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_sendrecv(sendbuf,sendcount,character_MPI,dest,sendtag,recvbuf,&
                & recvcount,character_MPI,source,recvtag,comm,status,errstat)
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to send or receive character data: '//&
                           & TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to send or receive character data'
      ENDIF
   ENDIF
   CALL error_MPI_comm('MPI_sendrecv',source,dest)
ENDIF
#else
recvbuf = ''
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_sendrecv_char

!========================================================================

SUBROUTINE comms_sendrecv_int(sendbuf,sendcount,dest,sendtag,recvbuf,&
                            & recvcount,source,recvtag,comm,ivarid)
!************************************************************************
!
! *comms_sendrecv_int*  Send-receive of integer data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_sendrecv
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, dest, ivarid, recvcount, recvtag, sendcount, &
                     & sendtag, source
INTEGER, INTENT(IN), DIMENSION(sendcount) :: sendbuf
INTEGER, INTENT(OUT), DIMENSION(recvcount) :: recvbuf

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*sendbuf*   INTEGER Starting address of send buffer
!*sendcount* INTEGER Number of data in send buffer
!*dest*      INTEGER Rank of destination process
!*sendtag*   INTEGER Message tag of send operation
!*recvbuf*   INTEGER Starting address of receive buffer
!*recvcount* INTEGER Number of data in receive buffer
!*source*    INTEGER Rank of source process
!*recvtag*   INTEGER Message tag of receive operation
!*comm*      INTEGER Communicator
!*ivarid*    INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (sendcount.EQ.0.AND.recvcount.EQ.0) RETURN
procname(pglev+1) = 'comms_sendrecv_int'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_sendrecv(sendbuf,sendcount,integer_MPI,dest,sendtag,recvbuf,&
                & recvcount,integer_MPI,source,recvtag,comm,status,errstat)
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to send or receive integer data: '//&
                           & TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to send or receive integer data'
      ENDIF
   ENDIF
   CALL error_MPI_comm('MPI_sendrecv',source,dest)
ENDIF
#else
recvbuf = 0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_sendrecv_int

!========================================================================

SUBROUTINE comms_sendrecv_log(sendbuf,sendcount,dest,sendtag,recvbuf,&
                            & recvcount,source,recvtag,comm,ivarid)
!************************************************************************
!
! *comms_sendrecv_log*  Send-receive of logical data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_sendrecv
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, dest, ivarid, recvcount, recvtag, sendcount, &
                     & sendtag, source
LOGICAL, INTENT(IN), DIMENSION(sendcount) :: sendbuf
LOGICAL, INTENT(OUT), DIMENSION(recvcount) :: recvbuf

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*sendbuf*   LOGICAL Starting address of send buffer
!*sendcount* INTEGER Number of data in send buffer
!*dest*      INTEGER Rank of destination process
!*sendtag*   INTEGER Message tag of send operation
!*recvbuf*   LOGICAL Starting address of receive buffer
!*recvcount* INTEGER Number of data in receive buffer
!*source*    INTEGER Rank of source process
!*recvtag*   INTEGER Message tag of receive operation
!*comm*      INTEGER Communicator
!*ivarid*    INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (sendcount.EQ.0.AND.recvcount.EQ.0) RETURN
procname(pglev+1) = 'comms_sendrecv_log'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_sendrecv(sendbuf,sendcount,logical_MPI,dest,sendtag,recvbuf,&
                & recvcount,logical_MPI,source,recvtag,comm,status,errstat)
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to send or receive logical data: '//&
                           & TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to send or receive logical data'
      ENDIF
   ENDIF
   CALL error_MPI_comm('MPI_sendrecv',source,dest)
ENDIF
#else
recvbuf = .FALSE.
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_sendrecv_log

!========================================================================

SUBROUTINE comms_sendrecv_real(sendbuf,sendcount,dest,sendtag,recvbuf,&
                             & recvcount,source,recvtag,comm,ivarid)
!************************************************************************
!
! *comms_sendrecv_real*  Send-receive of real data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI_comm
!
! MPI calls - MPI_sendrecv
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: comm, dest, ivarid, recvcount, recvtag, sendcount, &
                     & sendtag, source
REAL, INTENT(IN), DIMENSION(sendcount) :: sendbuf
REAL, INTENT(OUT), DIMENSION(recvcount) :: recvbuf

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*sendbuf*   REAL    Starting address of send buffer
!*sendcount* INTEGER Number of data in send buffer
!*dest*      INTEGER Rank of destination process
!*sendtag*   INTEGER Message tag of send operation
!*recvbuf*   REAL    Starting address of receive buffer
!*recvcount* INTEGER Number of data in receive buffer
!*source*    INTEGER Rank of source process
!*recvtag*   INTEGER Message tag of receive operation
!*comm*      INTEGER Communicator
!*ivarid*    INTEGER Variable id
!
!------------------------------------------------------------------------------
!


IF (sendcount.EQ.0.AND.recvcount.EQ.0) RETURN
procname(pglev+1) = 'comms_sendrecv_real'
CALL log_timer_in(npcc,ivarid,info=info)

#ifdef MPI
CALL MPI_sendrecv(sendbuf,sendcount,float_MPI,dest,sendtag,recvbuf,&
                & recvcount,float_MPI,source,recvtag,comm,status,errstat)
IF (errstat.NE.MPI_success) THEN
   IF (errchk) THEN
      IF (ivarid.GT.0) THEN
         CALL inquire_var(ivarid,f90_name=f90_name)
         WRITE (ioerr,'(A)') 'Unable to send or receive real data: '//&
                           & TRIM(f90_name)
      ELSE
         WRITE (ioerr,'(A)') 'Unable to send or receive real data'
      ENDIF
   ENDIF
   CALL error_MPI_comm('MPI_sendrecv',source,dest)
ENDIF
#else
recvbuf = 0.0
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_sendrecv_real

!========================================================================

SUBROUTINE comms_waitall(count,requests)
!************************************************************************
!
! *comms_waitall* Complete communications associated with array of requests
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Reference -
!
! Module calls - error_MPI
!
! MPI calls - MPI_waitall
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: count
INTEGER, INTENT(INOUT), DIMENSION(count) :: requests

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*count*    INTEGER Number of requests
!*requests* INTEGER Array of requests
!
!------------------------------------------------------------------------------
!
!*Local variables
!
#ifdef MPI
INTEGER, DIMENSION(MPI_status_size,count) :: statuses


procname(pglev+1) = 'comms_waitall'
CALL log_timer_in(npcc,info=info)

CALL MPI_waitall(count,requests,statuses,errstat)
IF (errstat.NE.MPI_success) CALL error_MPI('MPI_waitall')
#endif /*MPI*/

CALL log_timer_out(npcc,itm_MPI,info)


RETURN

END SUBROUTINE comms_waitall

!========================================================================

SUBROUTINE error_MPI(mpiname,abort)
!************************************************************************
!
! *error_MPI* Report error in MPI routine
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Module calls - error_abort
!
! MPI calls - MPI_error_string
!
!************************************************************************
!
USE switches
USE error_routines, ONLY: error_abort

!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: abort
CHARACTER (LEN=*), INTENT(IN) :: mpiname

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*mpiname*   CHAR    Name of MPI routine
!*abort*     LOGICAL Abort program if present and .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: abortx
#ifdef MPI
   CHARACTER (LEN=MPI_max_error_string) :: errmsg
   INTEGER :: ierr, lenerr
#endif /*MPI*/


IF (PRESENT(abort)) THEN
   abortx = abort
ELSEIF (iopt_MPI_abort.EQ.0) THEN
   abortx = .FALSE.
ELSE
   abortx = .TRUE.
ENDIF

#ifdef MPI
IF (errstat.NE.MPI_success) THEN
#else
IF (errstat.NE.0) THEN
#endif /*MPI*/
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      IF (.NOT.abortx) THEN
         WRITE (ioerr,'(A)') 'Error occurred in MPI routine '//TRIM(mpiname)
      ENDIF
#ifdef MPI
      CALL MPI_error_string(errstat,errmsg,lenerr,ierr)
      IF (ierr.NE.MPI_success) CALL error_abort('MPI_error_string',ierrno_MPI)
      WRITE (ioerr,'(A)') errmsg(1:lenerr)
#endif /*MPI*/
   ENDIF
   IF (abortx) CALL error_abort(mpiname,ierrno_MPI)
ENDIF


RETURN

END SUBROUTINE error_MPI

!========================================================================

SUBROUTINE error_MPI_comm(mpiname,idsrce,iddest)
!************************************************************************
!
! *error_MPI_comm* Report error in MPI send or receive communication
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)comms_MPI.F90  V2.0
!
! Description -
!
! Module calls - error_MPI
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: mpiname
INTEGER, INTENT(IN) :: iddest, idsrce

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*mpiname*   CHAR    Name of MPI routine
!*idsrce*    INTEGER Rank of source process
!*iddest*    INTEGER Rank of destination process
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cdest, csrce
     

IF (errchk) THEN
   WRITE (csrce,'(I12)') idsrce; csrce = ADJUSTL(csrce)
   WRITE (cdest,'(I12)') iddest; cdest = ADJUSTL(cdest)
   WRITE (ioerr,'(A)') 'Source process: '//csrce
   WRITE (ioerr,'(A)') 'Destination process: '//cdest
ENDIF
CALL error_MPI(mpiname)


RETURN

END SUBROUTINE error_MPI_comm


END MODULE comms_MPI
