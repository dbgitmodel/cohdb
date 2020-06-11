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
! *Relaxation_Zones* Relaxation scheme at open boundaries
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Relaxation_Zones.f90  V2.11
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - 
!
! Reference -
!
! Subroutines - define_rlxobc_spec, read_rlxobc_spec, relaxation_at_C,
!               relaxation_at_U, relaxation_at_V, relaxation_weights,
!               write_rlxobc_spec
!
!************************************************************************
!
!========================================================================

SUBROUTINE define_rlxobc_spec
!************************************************************************
!
! *define_rlxobc_spec* Define specifier arrays for relaxation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Relaxation_Zones.f90  V2.1.0
!
! Description - 
!
! Reference -
!
! Calling program - initialise_setup
!
! External calls - read_rlxobc_spec, usrdef_rlxobc_spec, write_rlxobc_spec
!
! Module calls - error_alloc
!
!************************************************************************
!
USE gridpars
USE iopars
USE paralpars
USE relaxation
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'define_rlxobc_spec'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (idirrlx(norlxzones),STAT=errstat)
CALL error_alloc('idirrlx',1,(/norlxzones/),kndint)
ALLOCATE (iposrlx(norlxzones),STAT=errstat)
CALL error_alloc('iposrlx',1,(/norlxzones/),kndint)
ALLOCATE (ityprlx(norlxzones),STAT=errstat)
CALL error_alloc('ityprlx',1,(/norlxzones/),kndint)
ALLOCATE (jposrlx(norlxzones),STAT=errstat)
CALL error_alloc('jposrlx',1,(/norlxzones/),kndint)
ALLOCATE (ncrlx(norlxzones),STAT=errstat)
CALL error_alloc('ncrlx',1,(/norlxzones/),kndint)
ALLOCATE (nrrlx(norlxzones),STAT=errstat)
CALL error_alloc('nrrlx',1,(/norlxzones/),kndint)

!
!2. Obtain specifier arrays
!--------------------------
!
!---read
IF (modfiles(io_rlxobc,1,1)%status.EQ.'R') THEN
   CALL read_rlxobc_spec
!---user-defined
ELSEIF (modfiles(io_rlxobc,1,1)%status.EQ.'N') THEN
   CALL usrdef_rlxobc_spec
ENDIF
   
!---write
IF (modfiles(io_rlxobc,1,2)%defined.AND.master) CALL write_rlxobc_spec

!
!3. Allocate weight factor arrays
!--------------------------------
!

IF (inodesrlx(1).EQ.1) THEN
   ALLOCATE (indexrlxatc(ncloc,nrloc,2),STAT=errstat)
   CALL error_alloc('indexrlxatc',3,(/ncloc,nrloc,2/),kndint)
   indexrlxatc = 0
   ALLOCATE (rlxwghtatc(ncloc,nrloc,2),STAT=errstat)
   CALL error_alloc('rlxwghtatc',3,(/ncloc,nrloc,2/),kndrtype)
   rlxwghtatc = 0.0
ENDIF
IF (inodesrlx(2).EQ.1) THEN
   ALLOCATE (indexrlxatuv(ncloc,nrloc,2),STAT=errstat)
   CALL error_alloc('indexrlxatuv',3,(/ncloc,nrloc,2/),kndint)
   indexrlxatuv = 0
   ALLOCATE (rlxwghtatuv(ncloc,nrloc,2),STAT=errstat)
   CALL error_alloc('rlxwghtatuv',3,(/ncloc,nrloc,2/),kndrtype)
   rlxwghtatuv = 0.0
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE define_rlxobc_spec

!========================================================================

SUBROUTINE read_rlxobc_spec
!************************************************************************
!
! *read_rlxobc_spec* Read specifier arrays for relaxation in standard
!                    format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Relaxation_Zones.f90  V2.1.0
!
! Description -
!
! Reference -
!
! Calling program - define_rlxobc_spec
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE relaxation
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: numvars
TYPE (FileParams):: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_rlxobc_spec'
CALL log_timer_in()

filepars = modfiles(io_rlxobc,1,1)

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
!
!2. Read specifier arrays
!------------------------
!

CALL read_vars(inodesrlx,filepars,1,varatts)
CALL read_vars(idirrlx,filepars,2,varatts)
CALL read_vars(ityprlx,filepars,3,varatts)
CALL read_vars(iposrlx,filepars,4,varatts)
CALL read_vars(jposrlx,filepars,5,varatts)
CALL read_vars(ncrlx,filepars,6,varatts)
CALL read_vars(nrrlx,filepars,7,varatts)

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_rlxobc,1,1) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE read_rlxobc_spec

!========================================================================

SUBROUTINE relaxation_at_C(psi,psiobu,psiobv,novars,iprofrlx)
!************************************************************************
!
! *relaxation_at_C* Apply relaxation scheme for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Relaxation_Zones.f90  V2.7.1
!
! Description -
!
! Reference -
!
! Calling program - transport_at_C_3d, transport_at_C_4d1, transport_at_C_4d2
!
! External calls -
!
! Module calls - combine_stats_glb, error_alloc
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE relaxation
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: combine_stats_glb
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: novars
INTEGER, INTENT(IN), DIMENSION(norlxzones) :: iprofrlx
REAL, INTENT(INOUT), DIMENSION(1-nhalo:ncloc+nhalo,&
                             & 1-nhalo:nrloc+nhalo,nz,novars) :: psi
REAL, INTENT(INOUT), DIMENSION(nobu,nz,novars) :: psiobu
REAL, INTENT(INOUT), DIMENSION(nobv,nz,novars) :: psiobv

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*novars*   INTEGER Number of variables
!*iprofrlx* INTEGER Array to select open boundary zones
!*psi*      REAL    Model array                                           [psi]
!*psiobu*   REAL    Vertical profiles at U-open boundaries
!                   (zero gradient condition if flagged)                  [psi]
!*psiobv*   REAL    Vertical profiles at V-open boundaries
!                   (zero gradient condition if flagged)                  [psi]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ii, ivar, j, jj, k, l, n
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: psix, psiy
REAL :: rx, ry


IF (ALL(iprofrlx.EQ.0)) RETURN

procname(pglev+1) = 'relaxation_at_C'
CALL log_timer_in()

!
!1. Allocate
!-----------
!

ALLOCATE (psix(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('psix',3,(/ncloc,nrloc,nz/),kndrtype)
psix = float_fill
ALLOCATE (psiy(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('psiy',3,(/ncloc,nrloc,nz/),kndrtype)
psiy = float_fill

!
!2. Combine open boundary profiles
!---------------------------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL combine_stats_glb(psiobu,nobu,nobuprocs,indexobuprocs,0,commall=.TRUE.)
   CALL combine_stats_glb(psiobv,nobv,nobvprocs,indexobvprocs,0,commall=.TRUE.)
ENDIF

!
!3. Apply scheme
!---------------
!
!3.1 X-direction
!---------------
!

ivar_310: DO ivar=1,novars
j_310: DO j=1,nrloc
i_310: DO i=1,ncloc

!  ---X-direction
   l = indexrlxatc(i,j,1)
   IF (l.NE.0) THEN
      n = MOD(l-1,norlxzones) + 1; ii = (l-1)/norlxzones + 1
      IF (iprofrlx(n).GT.0) THEN
         rx = rlxwghtatc(i,j,1)
         WHERE (ABS(psiobu(ii,:,ivar)-float_fill).GT.float_fill_eps)
            psix(i,j,:) = rx*psiobu(ii,:,ivar)+(1.0-rx)*psi(i,j,:,ivar)
         END WHERE
      ENDIF
   ENDIF

!  ---Y-direction
   l = indexrlxatc(i,j,2)
   IF (l.NE.0) THEN
      n = MOD(l-1,norlxzones) + 1; jj = (l-1)/norlxzones + 1
      IF (iprofrlx(n).GT.0) THEN
         ry = rlxwghtatc(i,j,2)
         WHERE (ABS(psiobv(jj,:,ivar)-float_fill).GT.float_fill_eps)
            psiy(i,j,:) = ry*psiobv(jj,:,ivar)+(1.0-ry)*psi(i,j,:,ivar)
         END WHERE
      ENDIF
   ENDIF

!  ---XY-directions
   k_311: DO k=1,nz
      IF ((ABS(psix(i,j,k)-float_fill).GT.float_fill_eps).AND.&
        & (ABS(psiy(i,j,k)-float_fill).GT.float_fill_eps)) THEN
         psi(i,j,k,ivar) = (rx*psix(i,j,k)+ry*psiy(i,j,k))/(rx+ry)
      ELSEIF (ABS(psix(i,j,k)-float_fill).GT.float_fill_eps) THEN
         psi(i,j,k,ivar) = psix(i,j,k)
      ELSEIF (ABS(psiy(i,j,k)-float_fill).GT.float_fill_eps) THEN
         psi(i,j,k,ivar) = psiy(i,j,k)
      ENDIF
   ENDDO k_311

ENDDO i_310
ENDDO j_310
ENDDO ivar_310

!
!4. Deallocate
!-------------
!

DEALLOCATE (psix,psiy)

CALL log_timer_out()


RETURN

END SUBROUTINE relaxation_at_C

!========================================================================

SUBROUTINE relaxation_at_U(psi,psiobu,nzdim,iprofrlx)
!************************************************************************
!
! *relaxation_at_U* Apply relaxation scheme for a quantity at U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Relaxation_Zones.f90  V2.7.1
!
! Description -
!
! Reference -
!
! Calling program - current_corr
!
! External calls -
!
! Module calls - combine_stats_glb
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE relaxation
USE switches
USE paral_comms, ONLY: combine_stats_glb
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: nzdim
INTEGER, INTENT(IN), DIMENSION(norlxzones) :: iprofrlx
REAL, INTENT(INOUT), DIMENSION(1-nhalo:ncloc+nhalo,&
                             & 1-nhalo:nrloc+nhalo,nzdim) :: psi
REAL, INTENT(INOUT), DIMENSION(nobu,nzdim) :: psiobu

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*nzdim*    INTEGER Vertical dimension of input array
!*iprofrlx* INTEGER Array to select open boundary zones
!*psi*      REAL    Model array                                           [psi]
!*psiobu*   REAL    Vertical profiles at U-open boundaries
!                   (zero gradient condition if flagged)                  [psi]
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ii, j, l , n
REAL :: rx


IF (ALL(iprofrlx.EQ.0)) RETURN

procname(pglev+1) = 'relaxation_at_U'
CALL log_timer_in()

!
!1. Combine open boundary profiles
!---------------------------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL combine_stats_glb(psiobu,nobu,nobuprocs,indexobuprocs,0,commall=.TRUE.)
ENDIF

!
!2. Apply scheme
!---------------
!

j_210: DO j=1,nrloc
i_210: DO i=1,ncloc
   l = indexrlxatuv(i,j,1)
   IF (l.NE.0) THEN
      n = MOD(l-1,norlxzones) + 1; ii = (l-1)/norlxzones + 1
      IF (iprofrlx(n).GT.0) THEN
         rx = rlxwghtatuv(i,j,1)
         WHERE (ABS(psiobu(ii,:)-float_fill).GT.float_fill_eps)
            psi(i,j,:) = rx*psiobu(ii,:)+(1.0-rx)*psi(i,j,:)
         END WHERE
      ENDIF
   ENDIF
ENDDO i_210
ENDDO j_210

CALL log_timer_out()


RETURN

END SUBROUTINE relaxation_at_U

!========================================================================

SUBROUTINE relaxation_at_V(psi,psiobv,nzdim,iprofrlx)
!************************************************************************
!
! *relaxation_at_VY* Apply relaxation scheme for a quantity at V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Relaxation_Zones.f90  V2.7.1
!
! Description -
!
! Reference -
!
! Calling program - current_corr
!
! External calls -
!
! Module calls - combine_stats_glb, error_alloc
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE relaxation
USE switches
USE paral_comms, ONLY: combine_stats_glb
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: nzdim
INTEGER, INTENT(IN), DIMENSION(norlxzones) :: iprofrlx
REAL, INTENT(INOUT), DIMENSION(1-nhalo:ncloc+nhalo,&
                             & 1-nhalo:nrloc+nhalo,nzdim) :: psi
REAL, INTENT(INOUT), DIMENSION(nobv,nzdim) :: psiobv

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*nzdim*    INTEGER Vertical dimension of input array
!*iprofrlx* INTEGER Array to select open boundary zones
!*psi*      REAL    Model array                                           [psi]
!*psiobv*   REAL    Vertical profiles at V-open boundaries
!                   (zero gradient condition if flagged)                  [psi]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j, jj, l, n
REAL :: ry


IF (ALL(iprofrlx.EQ.0)) RETURN

procname(pglev+1) = 'relaxation_at_V'
CALL log_timer_in()

!
!1. Combine open boundary profiles
!---------------------------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL combine_stats_glb(psiobv,nobv,nobvprocs,indexobvprocs,0,commall=.TRUE.)
ENDIF

!
!2. Apply scheme
!---------------
!

j_210: DO j=1,nrloc
i_210: DO i=1,ncloc
   l = indexrlxatuv(i,j,2)
   IF (l.NE.0) THEN
      n = MOD(l-1,norlxzones) + 1; jj = (l-1)/norlxzones + 1
      IF (iprofrlx(n).GT.0) THEN
         ry = rlxwghtatuv(i,j,2)
         WHERE (ABS(psiobv(jj,:)-float_fill).GT.float_fill_eps)
            psi(i,j,:) = ry*psiobv(jj,:)+(1.0-ry)*psi(i,j,:)
         END WHERE
      ENDIF
   ENDIF
ENDDO i_210
ENDDO j_210

CALL log_timer_out()


RETURN

END SUBROUTINE relaxation_at_V

!========================================================================

SUBROUTINE relaxation_weights
!************************************************************************
!
! *relaxation_weights* Weight factors for relaxation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Relaxation_Zones.f90  V2.2
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - find_obc_index_glb, local_proc, relax_factor
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE relaxation
USE grid_routines, ONLY: find_obc_index_glb, local_proc
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: relax_factor

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, idir, iglb, ii, irlx, ityp, i1glb, j, jglb, jj, jrlx, j1glb, &
         & n, nx, ny
REAL :: width
LOGICAL, DIMENSION(ncloc,nrloc) :: maskloc


procname(pglev+1) = 'relaxation_weights'
CALL log_timer_in()

!
!1. 'C'-type variables
!---------------------
!

IF (inodesrlx(1).GT.0) THEN

   maskloc = nodeatc(1:ncloc,1:nrloc).GT.0

   n_100: DO n=1,norlxzones
      ityp = ityprlx(n); idir = idirrlx(n)
      i1glb = iposrlx(n); j1glb = jposrlx(n)
      nx = ncrlx(n); ny = nrrlx(n)

!     ---West
      IF (idir.EQ.1) THEN
         jrlx_110: DO jrlx=1,ny
            jglb = jrlx + j1glb - 1; j = jglb - nr1loc + 1
            width = nx + 1
            ii = find_obc_index_glb(i1glb,jglb,'U  ')
            irlx_111: DO irlx=1,nx
               iglb = irlx + i1glb - 1; i = iglb - nc1loc + 1
               IF (local_proc(i,j)) THEN
                  IF (maskloc(i,j)) THEN
                     indexrlxatc(i,j,1) = n + (ii-1)*norlxzones
                     rlxwghtatc(i,j,1) = relax_factor(irlx,ityp,width)
                  ELSE
                     EXIT irlx_111
                  ENDIF
               ENDIF
            ENDDO irlx_111
         ENDDO jrlx_110

!     ---East
      ELSEIF (idir.EQ.2) THEN
         jrlx_120: DO jrlx=1,ny
            jglb = jrlx + j1glb - 1; j = jglb - nr1loc + 1
            width = nx + 1
            ii = find_obc_index_glb(i1glb,jglb,'U  ')
            irlx_121: DO irlx=1,nx
               iglb = i1glb - irlx; i = iglb - nc1loc + 1
               IF (local_proc(i,j)) THEN
                  IF (maskloc(i,j)) THEN
                     indexrlxatc(i,j,1) = n + (ii-1)*norlxzones
                     rlxwghtatc(i,j,1) = relax_factor(irlx,ityp,width)
                  ELSE
                     EXIT irlx_121
                  ENDIF
               ENDIF
            ENDDO irlx_121
         ENDDO jrlx_120

!     ---South
      ELSEIF (idir.EQ.3) THEN
         irlx_130: DO irlx=1,nx
            iglb = irlx + i1glb - 1; i = iglb - nc1loc + 1
            width = ny + 1
            jj = find_obc_index_glb(iglb,j1glb,'V  ')
            jrlx_131: DO jrlx=1,ny
               jglb = jrlx + j1glb - 1; j = jglb - nr1loc + 1
               IF (local_proc(i,j)) THEN
                  IF (maskloc(i,j)) THEN
                     indexrlxatc(i,j,2) = n + (jj-1)*norlxzones
                     rlxwghtatc(i,j,2) = relax_factor(jrlx,ityp,width)
                  ELSE
                     EXIT jrlx_131
                  ENDIF
               ENDIF
            ENDDO jrlx_131
         ENDDO irlx_130

!     ---North
      ELSEIF (idir.EQ.4) THEN
         irlx_140: DO irlx=1,nx
            iglb = irlx + i1glb - 1; i = iglb - nc1loc + 1
            width = ny + 1
            jj = find_obc_index_glb(iglb,j1glb,'V  ')
            jrlx_141: DO jrlx=1,ny
               jglb = j1glb - jrlx; j = jglb - nr1loc + 1
               IF (local_proc(i,j)) THEN
                  IF (maskloc(i,j)) THEN
                     indexrlxatc(i,j,2) = n + (jj-1)*norlxzones
                     rlxwghtatc(i,j,2) = relax_factor(jrlx,ityp,width)
                  ELSE
                     EXIT jrlx_141
                  ENDIF
               ENDIF
            ENDDO jrlx_141
         ENDDO irlx_140
      ENDIF

   ENDDO n_100

ENDIF

!
!2. 'U'-type variables
!---------------------
!

IF (inodesrlx(2).GT.0) THEN

   maskloc = node2du(1:ncloc,1:nrloc).EQ.2

   n_200: DO n=1,norlxzones
      ityp = ityprlx(n); idir = idirrlx(n)
      i1glb = iposrlx(n); j1glb = jposrlx(n)
      nx = ncrlx(n); ny = nrrlx(n)

!     ---West
      IF (idir.EQ.1) THEN
         jrlx_210: DO jrlx=1,ny
            jglb = jrlx + j1glb - 1; j = jglb - nr1loc + 1
            width = nx + 1
            ii = find_obc_index_glb(i1glb,jglb,'U  ')
            irlx_211: DO irlx=1,nx
               iglb = irlx + i1glb; i = iglb - nc1loc + 1
               IF (local_proc(i,j)) THEN
                  IF (maskloc(i,j)) THEN
                     indexrlxatuv(i,j,1) = n + (ii-1)*norlxzones
                     rlxwghtatuv(i,j,1) = relax_factor(irlx,ityp,width)
                  ELSE
                     EXIT irlx_211
                  ENDIF
               ENDIF
            ENDDO irlx_211
         ENDDO jrlx_210

!     ---East
      ELSEIF (idir.EQ.2) THEN
         jrlx_220: DO jrlx=1,ny
            jglb = jrlx + j1glb - 1; j = jglb - nr1loc + 1
            width = nx + 1
            ii = find_obc_index_glb(i1glb,jglb,'U  ')
            irlx_221: DO irlx=1,nx
               iglb = i1glb - irlx; i = iglb - nc1loc + 1
               IF (local_proc(i,j)) THEN
                  IF (maskloc(i,j)) THEN
                     indexrlxatuv(i,j,1) = n + (ii-1)*norlxzones
                     rlxwghtatuv(i,j,1) = relax_factor(irlx,ityp,width)
                  ELSE
                     EXIT irlx_221
                  ENDIF
               ENDIF
            ENDDO irlx_221
         ENDDO jrlx_220
      ENDIF

   ENDDO n_200
   
ENDIF

!
!3. 'V'-type variables
!---------------------
!

IF (inodesrlx(2).GT.0) THEN

   maskloc = node2dv(1:ncloc,1:nrloc).EQ.2

   n_300: DO n=1,norlxzones
      ityp = ityprlx(n); idir = idirrlx(n)
      i1glb = iposrlx(n); j1glb = jposrlx(n)
      nx = ncrlx(n); ny = nrrlx(n)

!     ---South
      IF (idir.EQ.3) THEN
         irlx_310: DO irlx=1,nx
            iglb = irlx + i1glb - 1; i = iglb - nc1loc + 1
            width = ny + 1
            jj = find_obc_index_glb(iglb,j1glb,'V  ')
            jrlx_311: DO jrlx=1,ny
               jglb = jrlx + j1glb; j = jglb - nr1loc + 1
               IF (local_proc(i,j)) THEN
                  IF (maskloc(i,j)) THEN
                     indexrlxatuv(i,j,2) = n + (jj-1)*norlxzones
                     rlxwghtatuv(i,j,2) = relax_factor(jrlx,ityp,width)
                  ELSE
                     EXIT jrlx_311
                  ENDIF
               ENDIF
            ENDDO jrlx_311
         ENDDO irlx_310

!     ---North
      ELSEIF (idir.EQ.4) THEN
         irlx_320: DO irlx=1,nx
            iglb = irlx + i1glb - 1; i = iglb - nc1loc + 1
            width = ny + 1
            jj = find_obc_index_glb(iglb,j1glb,'V  ')
            jrlx_321: DO jrlx=1,ny
               jglb = j1glb - jrlx; j = jglb - nr1loc + 1
               IF (local_proc(i,j)) THEN
                  IF (maskloc(i,j)) THEN
                     indexrlxatuv(i,j,2) = n + (jj-1)*norlxzones
                     rlxwghtatuv(i,j,2) = relax_factor(jrlx,ityp,width)
                  ELSE
                     EXIT jrlx_321
                  ENDIF
               ENDIF
            ENDDO jrlx_321
         ENDDO irlx_320
      ENDIF

   ENDDO n_300

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE relaxation_weights

!========================================================================

SUBROUTINE write_rlxobc_spec
!************************************************************************
!
! *write_rlxobc_spec* Write specifier arrays for relaxation in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Relaxation_Zones.f90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - define_rlxobc_spec
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, write_atts_mod,
!                write_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE paralpars
USE relaxation
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: numvars
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_rlxobc_spec'
CALL log_timer_in()

filepars = modfiles(io_rlxobc,1,2)

!
!1. Write file header
!--------------------
!
!---file attributes
CALL set_modfiles_atts(io_rlxobc,1,2)
filepars = modfiles(io_rlxobc,1,2)
numvars = filepars%novars

!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL set_modvars_atts(io_rlxobc,1,2,varatts,numvars)

!---write
CALL write_atts_mod(filepars,varatts,numvars)

!
!2. Write specifier arrays
!-------------------------
!

CALL write_vars(inodesrlx,filepars,1,varatts)
CALL write_vars(idirrlx,filepars,2,varatts)
CALL write_vars(ityprlx,filepars,3,varatts)
CALL write_vars(iposrlx,filepars,4,varatts)
CALL write_vars(jposrlx,filepars,5,varatts)
CALL write_vars(ncrlx,filepars,6,varatts)
CALL write_vars(nrrlx,filepars,7,varatts)

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_rlxobc,1,2) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE write_rlxobc_spec
