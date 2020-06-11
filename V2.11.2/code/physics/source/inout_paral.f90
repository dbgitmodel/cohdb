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

MODULE inout_paral
!************************************************************************
!
! *inout_paral* Routines for parallel I/O operations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.11
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description -
!
! Generic routines - combine_write_mod, combine_write_stats_glb,
!                    combine_write_stats_loc, combine_write_submod,
!                    read_distribute_mod
!
!************************************************************************
!
USE iopars
USE paralpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc

INTERFACE combine_write_mod
   MODULE PROCEDURE combine_write_mod_real_2d, combine_write_mod_real_3d, &
                  & combine_write_mod_real_4d, combine_write_mod_log_2d
END INTERFACE

INTERFACE combine_write_stats_glb
   MODULE PROCEDURE combine_write_stats_glb_real_1d, &
                  & combine_write_stats_glb_real_2d, &
                  & combine_write_stats_glb_real_3d, &
                  & combine_write_stats_glb_real_4d
END INTERFACE

INTERFACE combine_write_stats_loc
   MODULE PROCEDURE combine_write_stats_loc_real_1d, &
                  & combine_write_stats_loc_real_2d, &
                  & combine_write_stats_loc_real_3d
END INTERFACE

INTERFACE combine_write_submod
   MODULE PROCEDURE combine_write_submod_real_2d, &
                  & combine_write_submod_real_3d, combine_write_submod_real_4d
END INTERFACE

INTERFACE read_distribute_mod
   MODULE PROCEDURE read_distribute_mod_real_2d, read_distribute_mod_real_3d, &
                  & read_distribute_mod_real_4d, read_distribute_mod_log_2d
END INTERFACE


CONTAINS

!========================================================================

SUBROUTINE combine_write_mod_real_2d(realloc,filepars,varid,lbounds,ubounds,&
                                   & cnode,varatts,comm,commtype)
!************************************************************************
!
! *combine_write_mod_real_2d* Write a global real 2-D model array combined on
!                             the root process from local arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.10.1
!
! Description -
!
! Module calls - combine_mod, error_alloc, write_vars
!
!************************************************************************
!
USE datatypes
USE gridpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: global_mask
USE inout_routines, ONLY: write_vars
USE paral_comms, ONLY: combine_mod

!
!*  Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
TYPE(FileParams), INTENT(INOUT) :: filepars
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):) :: realloc
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*realloc*   REAL     Local data array
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*lbounds*   INTEGER  Lower bounds of local array
!*ubounds*   INTEGER  Upper bounds of local array
!*cnode*     CHAR     Nodal type of the data on the model grid
!                     ('C','W','U','V')
!*varatts*   DERIVED  Variable attributes
!*comm*      INTEGER  Communicator for combine operation
!*commtype*  INTEGER  Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: packing
INTEGER :: commop, itype, ivarid, nowet
INTEGER, DIMENSION(2) :: lbnds
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskglb
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: realdat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: realglb


procname(pglev+1) = 'combine_write_mod_real_2d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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
!2. Initialise parameters and arrays
!-----------------------------------
!
!---packing
packing = filepars%packing.AND.master

!---number of wet points
IF (packing) THEN
   SELECT CASE (TRIM(cnode))
      CASE ('C','W')
         nowet = noseaatc
      CASE ('U')
         nowet = noseaatu
      CASE ('V')
         nowet = noseaatv
   END SELECT
ENDIF      

!---allocate global data array
ALLOCATE (realglb(nc,nr),STAT=errstat)
CALL error_alloc('realglb',2,(/nc,nr/),kndrtype)

!---data array for packing
IF (packing) THEN
   ALLOCATE (maskglb(nc,nr),STAT=errstat)
   CALL error_alloc('maskglb',2,(/nc,nr/),kndlog)
   ALLOCATE (realdat(nowet),STAT=errstat)
   CALL error_alloc('realdat',1,(/nowet/),kndrtype)
ENDIF

!---global mask
IF (packing) CALL global_mask(maskglb,cnode,1)

!
!3. Combine data
!---------------
!

lbnds = 1
CALL combine_mod(realglb,realloc(1:ncloc,1:nrloc),lbnds,ivarid,0.0,0,&
               & idmaster,commop,itype)

!
!4. Pack data using land mask
!----------------------------
!

IF (packing) realdat = PACK(realglb,MASK=maskglb)

!
!5. Write data
!-------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realdat,filepars,varid,varatts=varatts)
   ELSE
      CALL write_vars(realdat,filepars,varid)
   ENDIF
ELSEIF (master) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realglb,filepars,varid,varatts=varatts)
   ELSE
      CALL write_vars(realglb,filepars,varid)
   ENDIF
ENDIF

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (realglb)
IF (packing) DEALLOCATE (maskglb,realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE combine_write_mod_real_2d

!========================================================================

SUBROUTINE combine_write_mod_real_3d(realloc,filepars,varid,lbounds,ubounds,&
                                   & cnode,varatts,vecids,comm,commtype)
!************************************************************************
!
! *combine_write_mod_real_3d* Write a global real 3-D model array combined on
!                             the root process from local arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.10.1
!
! Description -
!
! Module calls - combine_mod, error_alloc, write_vars
!
!************************************************************************
!
USE datatypes
USE gridpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: global_mask
USE inout_routines, ONLY: write_vars
USE paral_comms, ONLY: combine_mod

!
!*  Arguments
!
CHARACTER (LEN=lennode) :: cnode
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
TYPE(FileParams), INTENT(INOUT) :: filepars
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):,lbounds(3):) :: realloc
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realloc*  REAL     Local data array
!*filepars* DERIVED  File attributes
!*nzdim*    INTEGER
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*lbounds*  INTEGER  Lower bounds of local array
!*ubounds*  INTEGER  Upper bounds of local array
!*cnode*    CHAR     Nodal type of the data on the model grid
!                    ('C','W','U','V')
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!*comm*     INTEGER  Communicator for combine operation
!*commtype* INTEGER  Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: packing
INTEGER :: commop, itype, ivar, ivarid, k, nowet
INTEGER, DIMENSION(3) :: lbnds
INTEGER, DIMENSION(SIZE(realloc,DIM=3)) :: idvars
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskglb
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: realdat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realglb


procname(pglev+1) = 'combine_write_mod_real_3d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSE
   idvars = (/(ivar,ivar=1,SIZE(realloc,DIM=3))/)
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---packing
packing = filepars%packing.AND.master

!---number of wet points
IF (packing) THEN
   SELECT CASE (TRIM(cnode))
      CASE ('C','W')
         nowet = noseaatc
      CASE ('U')
         nowet = noseaatu
      CASE ('V')
         nowet = noseaatv
   END SELECT
ENDIF

!---allocate global data array
ALLOCATE (realglb(nc,nr,lbounds(3):ubounds(3)),STAT=errstat)
CALL error_alloc('realglb',3,(/nc,nr,ubounds(3)-lbounds(3)+1/),kndrtype)

!---data array for packing
IF (packing) THEN
   ALLOCATE (maskglb(nc,nr),STAT=errstat)
   CALL error_alloc('maskglb',2,(/nc,nr/),kndlog)
   ALLOCATE (realdat(nowet,lbounds(3):ubounds(3)),STAT=errstat)
   CALL error_alloc('realdat',2,(/nowet,ubounds(3)-lbounds(3)+1/),kndrtype)
ENDIF

!---global mask
IF (packing) CALL global_mask(maskglb,cnode,1)

!
!3. Combine data
!---------------
!

lbnds = 1
CALL combine_mod(realglb,realloc(1:ncloc,1:nrloc,:),lbnds,ivarid,0.0,0,&
               & idmaster,commop,itype)

!
!4. Pack data using land mask
!----------------------------
!

IF (packing) THEN
   k_510: DO k=lbounds(3),ubounds(3)
      realdat(:,k) = PACK(realglb(:,:,k),MASK=maskglb)
   ENDDO k_510
ENDIF

!
!5. Write data
!-------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realdat,filepars,varid,varatts=varatts,vecids=idvars)
   ELSE
      CALL write_vars(realdat,filepars,varid,vecids=idvars)
   ENDIF
ELSEIF (master) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realglb,filepars,varid,varatts=varatts,vecids=idvars)
   ELSE
      CALL write_vars(realglb,filepars,varid,vecids=idvars)
   ENDIF
ENDIF

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (realglb)
IF (packing) DEALLOCATE (maskglb,realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE combine_write_mod_real_3d

!========================================================================

SUBROUTINE combine_write_mod_real_4d(realloc,filepars,varid,lbounds,ubounds,&
                                   & cnode,varatts,vecids,comm,commtype)
!************************************************************************
!
! *combine_write_mod_real_4d* Write a global real 4-D model array combined on
!                             the root process from local arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.10.1
!
! Description -
!
! Module calls - combine_mod, error_alloc, write_vars
!
!************************************************************************
!
USE datatypes
USE gridpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: global_mask
USE inout_routines, ONLY: write_vars
USE paral_comms, ONLY: combine_mod

!
!*  Arguments
!
CHARACTER (LEN=lennode) :: cnode
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
TYPE(FileParams), INTENT(INOUT) :: filepars
INTEGER, INTENT(IN), DIMENSION(4) :: lbounds, ubounds
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):,lbounds(3):,&
                          & lbounds(4):) :: realloc
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realloc*  REAL     Local data array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*lbounds*  INTEGER  Lower bounds of local array
!*ubounds*  INTEGER  Upper bounds of local array
!*cnode*    CHAR     Nodal type of the data on the model grid
!                    ('C','W','U','V')
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!*comm*     INTEGER  Communicator for combine operation
!*commtype* INTEGER  Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: packing
INTEGER :: commop, itype, ivar, ivarid, k, l, nowet
INTEGER, DIMENSION(4) :: lbnds
INTEGER, DIMENSION(SIZE(realloc,DIM=4)) :: idvars
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskglb
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realdat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: realglb


procname(pglev+1) = 'combine_write_mod_real_4d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSE
   idvars = (/(ivar,ivar=1,SIZE(realloc,DIM=4))/)
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---packing
packing = filepars%packing.AND.master

!---number of wet points
IF (packing) THEN
   SELECT CASE (TRIM(cnode))
      CASE ('C','W')
         nowet = noseaatc
      CASE ('U')
         nowet = noseaatu
      CASE ('V')
         nowet = noseaatv
   END SELECT
ENDIF

!---allocate global data array
ALLOCATE (realglb(nc,nr,lbounds(3):ubounds(3),lbounds(4):ubounds(4)),&
        & STAT=errstat)
CALL error_alloc('realglb',4,&
              & (/nc,nr,ubounds(3)-lbounds(3)+1,ubounds(4)-lbounds(4)+1/),&
              & kndrtype)

!---data array for packing
IF (packing) THEN
   ALLOCATE (maskglb(nc,nr),STAT=errstat)
   CALL error_alloc('maskglb',2,(/nc,nr/),kndlog)
   ALLOCATE (realdat(nowet,lbounds(3):ubounds(3),lbounds(4):ubounds(4)),&
           & STAT=errstat)
   CALL error_alloc('realdat',3,&
                 & (/nowet,ubounds(3)-lbounds(3)+1,ubounds(4)-lbounds(4)+1/),&
                 & kndrtype)
ENDIF

!---global mask
IF (packing) CALL global_mask(maskglb,cnode,1)

!
!3. Combine data
!---------------
!

lbnds = 1
CALL combine_mod(realglb,realloc(1:ncloc,1:nrloc,:,:),lbnds,ivarid,0.0,0,&
               & idmaster,commop,itype)

!
!4. Pack data using land mask
!----------------------------
!

IF (packing) THEN
   l_510: DO l=lbounds(4),ubounds(4)
   k_510: DO k=lbounds(3),ubounds(3)
      realdat(:,k,l) = PACK(realglb(:,:,k,l),MASK=maskglb)
   ENDDO k_510
   ENDDO l_510
ENDIF

!
!5. Write data
!-------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realdat,filepars,varid,varatts=varatts,vecids=idvars)
   ELSE
      CALL write_vars(realdat,filepars,varid,vecids=idvars)
   ENDIF
ELSEIF (master) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realglb,filepars,varid,varatts=varatts,vecids=idvars)
   ELSE
      CALL write_vars(realglb,filepars,varid,vecids=idvars)
   ENDIF
ENDIF

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (realglb)
IF (packing) DEALLOCATE (maskglb,realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE combine_write_mod_real_4d

!========================================================================

SUBROUTINE combine_write_mod_log_2d(logloc,filepars,varid,lbounds,ubounds,&
                                  & cnode,varatts,comm,commtype)
!************************************************************************
!
! *combine_write_mod_log_2d* Write a global logical 2-D model array combined on
!                            the root process from local arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.11
!
! Description -
!
! Module calls - combine_mod, error_alloc, write_vars
!
!************************************************************************
!
USE datatypes
USE gridpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: global_mask
USE inout_routines, ONLY: write_vars
USE paral_comms, ONLY: combine_mod

!
!*  Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
TYPE(FileParams), INTENT(INOUT) :: filepars
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds, ubounds
LOGICAL, INTENT(IN), DIMENSION(lbounds(1):,lbounds(2):) :: logloc
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*logloc*    LOGICAL  Local data array
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*lbounds*   INTEGER  Lower bounds of local array
!*ubounds*   INTEGER  Upper bounds of local array
!*cnode*     CHAR     Nodal type of the data on the model grid
!                     ('C','W','U','V')
!*varatts*   DERIVED  Variable attributes
!*comm*      INTEGER  Communicator for combine operation
!*commtype*  INTEGER  Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: packing
INTEGER :: commop, itype, ivarid, nowet
INTEGER, DIMENSION(2) :: lbnds
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskglb
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:) :: logdat
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: logglb


procname(pglev+1) = 'combine_write_mod_log_2d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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
!2. Initialise parameters and arrays
!-----------------------------------
!
!---packing
packing = filepars%packing.AND.master

!---number of wet points
IF (packing) THEN
   SELECT CASE (TRIM(cnode))
      CASE ('C','W')
         nowet = noseaatc
      CASE ('U')
         nowet = noseaatu
      CASE ('V')
         nowet = noseaatv
   END SELECT
ENDIF      

!---allocate global data array
ALLOCATE (logglb(nc,nr),STAT=errstat)
CALL error_alloc('realglb',2,(/nc,nr/),kndlog)

!---data array for packing
IF (packing) THEN
   ALLOCATE (maskglb(nc,nr),STAT=errstat)
   CALL error_alloc('maskglb',2,(/nc,nr/),kndlog)
   ALLOCATE (logdat(nowet),STAT=errstat)
   CALL error_alloc('realdat',1,(/nowet/),kndlog)
ENDIF

!---global mask
IF (packing) CALL global_mask(maskglb,cnode,1)

!
!3. Combine data
!---------------
!

lbnds = 1
CALL combine_mod(logglb,logloc(1:ncloc,1:nrloc),lbnds,ivarid,.FALSE.,0,&
               & idmaster,commop,itype)

!
!4. Pack data using land mask
!----------------------------
!

IF (packing) logdat = PACK(logglb,MASK=maskglb)

!
!5. Write data
!-------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(logdat,filepars,varid,varatts=varatts)
   ELSE
      CALL write_vars(logdat,filepars,varid)
   ENDIF
ELSEIF (master) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(logglb,filepars,varid,varatts=varatts)
   ELSE
      CALL write_vars(logglb,filepars,varid)
   ENDIF
ENDIF

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (logglb)
IF (packing) DEALLOCATE (maskglb,logdat)

CALL log_timer_out()


RETURN

END SUBROUTINE combine_write_mod_log_2d

!========================================================================

SUBROUTINE combine_write_stats_glb_real_1d(realglb,filepars,varid,maxstats,&
                                         & nostatprocs,lstatprocs,varatts,&
                                         & comm,commtype)
!************************************************************************
!
! *combine_write_stats_glb_real_1d* Combine and write a global real vector
!                                   defined at stations on different domains
!                          
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.6
!
! Description -
!
! Module calls - combine_stats_glb, write_vars
!
!************************************************************************
!
USE datatypes
USE switches
USE inout_routines, ONLY: write_vars
USE paral_comms, ONLY: combine_stats_glb

!
!*  Arguments
!
INTEGER, INTENT(IN) :: maxstats, varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
TYPE(FileParams), INTENT(INOUT) :: filepars
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
REAL, INTENT(INOUT), DIMENSION(:) :: realglb
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type      Purpose
!------------------------------------------------------------------------------
!*realglb*    REAL    Global real array
!*filepars*   DERIVED File attributes
!*varid*      INTEGER If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*maxstats*   INTEGER First dimension of array lstatprocs
!*nostatprocs*INTEGER Number of stations in each local array
!*lstatprocs* INTEGER Indices of local stations in global array
!*varatts*    DERIVED Variable attributes
!*comm*       INTEGER Communicator for combine operation
!*commtype*   INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, itype, ivarid


IF (SIZE(realglb).EQ.0) RETURN
procname(pglev+1) = 'combine_write_stats_glb_real_1d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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
!2. Combine data
!---------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL combine_stats_glb(realglb,maxstats,nostatprocs,lstatprocs,ivarid,&
                        & idmaster,commop,itype)
ENDIF

!
!3. Write data
!-------------
!

IF (master) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realglb,filepars,varid,varatts=varatts)
   ELSE
      CALL write_vars(realglb,filepars,varid)
   ENDIF
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE combine_write_stats_glb_real_1d

!========================================================================

SUBROUTINE combine_write_stats_glb_real_2d(realglb,filepars,varid,maxstats,&
                                         & nostatprocs,lstatprocs,varatts,&
                                         & vecids,comm,commtype)
!************************************************************************
!
! *combine_write_stats_glb_real_2d* Combine and write a global real 2-D array
!                                   defined at stations on different domains
!                          
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.6
!
! Description -
!
! Module calls - combine_stats_glb, write_vars
!
!************************************************************************
!
USE datatypes
USE switches
USE inout_routines, ONLY: write_vars
USE paral_comms, ONLY: combine_stats_glb

!
!*  Arguments
!
INTEGER, INTENT(IN) :: maxstats, varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
TYPE(FileParams), INTENT(INOUT) :: filepars
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL, INTENT(INOUT), DIMENSION(:,:) :: realglb
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type      Purpose
!------------------------------------------------------------------------------
!*realglb*    REAL    Global real array
!*filepars*   DERIVED File attributes
!*varid*      INTEGER If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*maxstats*   INTEGER First dimension of array lstatprocs
!*nostatprocs*INTEGER Number of stations in each local array
!*lstatprocs* INTEGER Indices of local stations in global array
!*varatts*    DERIVED Variable attributes
!*vecids*     INTEGER Arrays of variable IDs if varid = 0
!*comm*       INTEGER Communicator for combine operation
!*commtype*   INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, itype, ivar, ivarid
INTEGER, DIMENSION(SIZE(realglb,DIM=2)) :: idvars


IF (SIZE(realglb).EQ.0) RETURN
procname(pglev+1) = 'combine_write_stats_glb_real_2d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSE
   idvars = (/(ivar,ivar=1,SIZE(realglb,DIM=2))/)
ENDIF

!
!2. Combine data
!---------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL combine_stats_glb(realglb,maxstats,nostatprocs,lstatprocs,ivarid,&
                        & idmaster,commop,itype)
ENDIF

!
!3. Write data
!-------------
!

IF (master) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realglb,filepars,varid,varatts=varatts)
   ELSE
      CALL write_vars(realglb,filepars,varid,vecids=idvars)
   ENDIF
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE combine_write_stats_glb_real_2d

!========================================================================

SUBROUTINE combine_write_stats_glb_real_3d(realglb,filepars,varid,maxstats,&
                                         & nostatprocs,lstatprocs,varatts,&
                                         & vecids,comm,commtype)
!************************************************************************
!
! *combine_write_stats_glb_real_3d* Combine and write a global real 3-D array
!                                   defined at stations on different domains
!                          
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.6
!
! Description -
!
! Module calls - combine_stats_glb, write_vars
!
!************************************************************************
!
USE datatypes
USE switches
USE inout_routines, ONLY: write_vars
USE paral_comms, ONLY: combine_stats_glb

!
!*  Arguments
!
INTEGER, INTENT(IN) :: maxstats, varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
TYPE(FileParams), INTENT(INOUT) :: filepars
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL, INTENT(INOUT), DIMENSION(:,:,:) :: realglb
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type      Purpose
!------------------------------------------------------------------------------
!*realglb*    REAL    Global real array
!*filepars*   DERIVED File attributes
!*varid*      INTEGER If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*maxstats*   INTEGER First dimension of array lstatprocs
!*nostatprocs*INTEGER Number of stations in each local array
!*lstatprocs* INTEGER Indices of local stations in global array
!*varatts*    DERIVED Variable attributes
!*vecids*     INTEGER Arrays of variable IDs if varid = 0
!*comm*       INTEGER Communicator for combine operation
!*commtype*   INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, itype, ivar, ivarid
INTEGER, DIMENSION(SIZE(realglb,DIM=3)) :: idvars


IF (SIZE(realglb).EQ.0) RETURN
procname(pglev+1) = 'combine_write_stats_glb_real_3d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSE
   idvars = (/(ivar,ivar=1,SIZE(realglb,DIM=3))/)
ENDIF

!
!2. Combine data
!---------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL combine_stats_glb(realglb,maxstats,nostatprocs,lstatprocs,ivarid,&
                        & idmaster,commop,itype)
ENDIF

!
!3. Write data
!-------------
!

IF (master) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realglb,filepars,varid,varatts=varatts)
   ELSE
      CALL write_vars(realglb,filepars,varid,vecids=idvars)
   ENDIF
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE combine_write_stats_glb_real_3d

!========================================================================

SUBROUTINE combine_write_stats_glb_real_4d(realglb,filepars,varid,maxstats,&
                                         & nostatprocs,lstatprocs,varatts,&
                                         & vecids,comm,commtype)
!************************************************************************
!
! *combine_write_stats_glb_real_4d* Combine and write a global real 4-D array
!                                   defined at stations on different domains
!                          
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.6
!
! Description -
!
! Module calls - combine_stats_glb, write_vars
!
!************************************************************************
!
USE datatypes
USE switches
USE inout_routines, ONLY: write_vars
USE paral_comms, ONLY: combine_stats_glb

!
!*  Arguments
!
INTEGER, INTENT(IN) :: maxstats, varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
TYPE(FileParams), INTENT(INOUT) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
REAL, INTENT(INOUT), DIMENSION(:,:,:,:) :: realglb
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type      Purpose
!------------------------------------------------------------------------------
!*realglb*    REAL    Global real array
!*filepars*   DERIVED File attributes
!*varid*      INTEGER If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*maxstats*   INTEGER First dimension of array lstatprocs
!*nostatprocs*INTEGER Number of stations in each local array
!*lstatprocs* INTEGER Indices of local stations in global array
!*varatts*    DERIVED Variable attributes
!*vecids*     INTEGER Arrays of variable IDs if varid = 0
!*comm*       INTEGER Communicator for combine operation
!*commtype*   INTEGER Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, itype, ivar, ivarid
INTEGER, DIMENSION(SIZE(realglb,DIM=4)) :: idvars


IF (SIZE(realglb).EQ.0) RETURN
procname(pglev+1) = 'combine_write_stats_glb_real_4d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSE
   idvars = (/(ivar,ivar=1,SIZE(realglb,DIM=4))/)
ENDIF

!
!2. Combine data
!---------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL combine_stats_glb(realglb,maxstats,nostatprocs,lstatprocs,ivarid,&
                        & idmaster,commop,itype)
ENDIF

!
!3. Write data
!-------------
!

IF (master) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realglb,filepars,varid,varatts=varatts)
   ELSE
      CALL write_vars(realglb,filepars,varid,vecids=idvars)
   ENDIF
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE combine_write_stats_glb_real_4d

!========================================================================

SUBROUTINE combine_write_stats_loc_real_1d(realloc,filepars,varid,maxstats,&
                                         & nostatsglb,nostatprocs,lstatprocs,&
                                         & varatts,comm,commtype)
!************************************************************************
!
! *combine_write_stats_loc_real_1d* Write global station data combined on the
!                                   root process from local vector data arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.6
!
! Description -
!
! Module calls - combine_stats_loc, error_alloc, write_vars
!
!************************************************************************
!
USE datatypes
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: write_vars
USE paral_comms, ONLY: combine_stats_loc

!
!*  Arguments
!
INTEGER, INTENT(IN) :: maxstats, nostatsglb, varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
TYPE(FileParams), INTENT(INOUT) :: filepars
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
REAL, INTENT(IN), DIMENSION(:) :: realloc
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*realloc*     REAL     Local station array with size given by number of local
!                       stations
!*filepars*    DERIVED  File attributes
!*varid*       INTEGER  If >0, variable id of input array.
!                       If =0, last dimension is a variable dimension
!*maxstats*    INTEGER  First dimension of array lstatprocs
!*nostatsglb*  INTEGER  Number of stations in global array
!*nostatprocs* INTEGER  Number of stations in each local array
!*lstatprocs*  INTEGER  Indices of local stations in global array
!*varatts*     DERIVED  Variable attributes
!*comm*        INTEGER  Communicator for combine operation
!*commtype*    INTEGER  Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: commop, itype, ivarid
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: realglb


procname(pglev+1) = 'combine_write_stats_loc_real_1d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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
!2. Allocate array
!-----------------
!

ALLOCATE (realglb(nostatsglb),STAT=errstat)
CALL error_alloc('realglb',1,(/nostatsglb/),kndrtype)
realglb = 0.0

!
!3. Combine data
!---------------
!

CALL combine_stats_loc(realglb,realloc,maxstats,nostatprocs,lstatprocs,ivarid,&
                     & idmaster,commop,itype)

!
!4. Write data
!-------------
!

IF (master) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realglb,filepars,varid,varatts=varatts)
   ELSE
      CALL write_vars(realglb,filepars,varid)
   ENDIF
ENDIF

!
!5. Deallocate array
!-------------------
!

DEALLOCATE (realglb)

CALL log_timer_out()


RETURN

END SUBROUTINE combine_write_stats_loc_real_1d

!========================================================================

SUBROUTINE combine_write_stats_loc_real_2d(realloc,filepars,varid,maxstats,&
                                         & nostatsglb,nostatprocs,lstatprocs,&
                                         & reduced,varatts,vecids,comm,commtype)
!************************************************************************
!
! *combine_write_stats_loc_real_2d* Write global station data combined on the
!                                   root process from local 2-D data arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.6
!
! Description -
!
! Module calls - combine_stats_loc, error_alloc, write_vars
!
!************************************************************************
!
USE datatypes
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: write_vars
USE paral_comms, ONLY: combine_stats_loc

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: reduced
INTEGER, INTENT(IN) :: maxstats, nostatsglb, varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
TYPE(FileParams), INTENT(INOUT) :: filepars
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL, INTENT(IN), DIMENSION(:,:) :: realloc
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*realloc*     REAL     Local array on input
!                 first dimension: number of local stations
!                 second dimension: number of variables
!*filepars*    DERIVED  File attributes
!*varid*       INTEGER  If >0, variable id of input array.
!                       If =0, last dimension is a variable dimension
!*maxstats*    INTEGER  First dimension of array lstatprocs
!*nostatsglb*  INTEGER  Number of stations in global array
!*nostatprocs* INTEGER  Number of stations in each local array
!*lstatprocs*  INTEGER  Indices of local stations in global array
!*reduced*     LOGICAL  If .TRUE., the last array dimension is taken within
!                       the first one so that the rank of the output array is
!                       reduced by 1 
!*varatts*     DERIVED  Variable attributes
!*vecids*      INTEGER  Arrays of variable IDs if varid = 0
!*comm*        INTEGER  Communicator for combine operation
!*commtype*    INTEGER  Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: reducedx
INTEGER :: commop, itype, ivarid, l, novals
INTEGER, DIMENSION(2) :: ind
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: realglb2
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: realglb


procname(pglev+1) = 'combine_write_stats_loc_real_2d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (PRESENT(reduced)) THEN
   reducedx = reduced
ELSE
   reducedx = .TRUE.
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
!2. Allocate array
!-----------------
!

novals = SIZE(realloc,DIM=2)
ALLOCATE (realglb(nostatsglb,novals),STAT=errstat)
CALL error_alloc('realglb',2,(/nostatsglb,novals/),kndrtype)
realglb = 0.0
IF (reduced) THEN
   ALLOCATE (realglb2(nostatsglb*novals),STAT=errstat)
   CALL error_alloc('realglb',2,(/nostatsglb*novals/),kndrtype)
   realglb2 = 0.0
ENDIF

!
!3. Combine data
!---------------
!

CALL combine_stats_loc(realglb,realloc,maxstats,nostatprocs,lstatprocs,ivarid,&
                     & idmaster,commop,itype)

!
!4. Reshape if needed
!--------------------
!

IF (reduced) THEN
   ind(1) = 1; ind(2) = nostatsglb
   l_410: DO l=1,novals
      realglb2(ind(1):ind(2)) = realglb(:,l)
      ind = ind + nostatsglb
   ENDDO l_410
ENDIF

!
!5. Write data
!-------------
!

IF (master) THEN
   IF (reduced) THEN
      IF (PRESENT(varatts)) THEN
         CALL write_vars(realglb2,filepars,varid,varatts=varatts)
      ELSE
         CALL write_vars(realglb2,filepars,varid)
      ENDIF
   ELSE
      IF (PRESENT(varatts).AND.PRESENT(vecids)) THEN
         CALL write_vars(realglb,filepars,varid,varatts=varatts,vecids=vecids)
      ELSEIF (PRESENT(varatts).AND.(.NOT.PRESENT(vecids))) THEN
         CALL write_vars(realglb,filepars,varid,varatts=varatts)
      ELSEIF (.NOT.PRESENT(varatts).AND.PRESENT(vecids)) THEN
         CALL write_vars(realglb,filepars,varid,vecids=vecids)
      ELSEIF (.NOT.(PRESENT(varatts).OR.PRESENT(vecids))) THEN
         CALL write_vars(realglb,filepars,varid)
      ENDIF
   ENDIF
ENDIF

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (realglb)
IF (reduced) DEALLOCATE (realglb2)

CALL log_timer_out()


RETURN

END SUBROUTINE combine_write_stats_loc_real_2d

!========================================================================

SUBROUTINE combine_write_stats_loc_real_3d(realloc,filepars,varid,maxstats,&
                                         & nostatsglb,nostatprocs,lstatprocs,&
                                         & reduced,varatts,vecids,comm,commtype)
!************************************************************************
!
! *combine_write_stats_loc_real_3d* Write global station data combined on the
!                                   root process from local 3-D data arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.6
!
! Description -
!
! Module calls - combine_stats_loc, error_alloc, write_vars
!
!************************************************************************
!
USE datatypes
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: write_vars
USE paral_comms, ONLY: combine_stats_loc

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: reduced
INTEGER, INTENT(IN) :: maxstats, nostatsglb, varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
TYPE(FileParams), INTENT(INOUT) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nostatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: lstatprocs
REAL, INTENT(IN), DIMENSION(:,:,:) :: realloc
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*realloc*     REAL     Local array on input
!                 first dimension: number of local stations
!                 second dimension: vertical dimension
!                 third dimension: number of variables
!*filepars*    DERIVED  File attributes
!*varid*       INTEGER  If >0, variable id of input array.
!                       If =0, last dimension is a variable dimension
!*maxstats*    INTEGER  First dimension of array lstatprocs
!*nostatsglb*  INTEGER  Number of stations in global array
!*nostatprocs* INTEGER  Number of stations in each local array
!*lstatprocs*  INTEGER  Indices of local stations in global array
!*reduced*     LOGICAL  If .TRUE., the last array dimension is taken within
!                       the first one so that the rank of the output array is
!                       reduced by 1
!*varatts*     DERIVED  Variable attributes
!*vecids*      INTEGER  Arrays of variable IDs if varid = 0
!*comm*        INTEGER  Communicator for combine operation
!*commtype*    INTEGER  Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: reducedx
INTEGER :: commop, itype, ivarid, l, novals, nzdim
INTEGER, DIMENSION(2) :: ind
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: realglb2
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realglb


procname(pglev+1) = 'combine_write_stats_loc_real_3d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

IF (PRESENT(reduced)) THEN
   reducedx = reduced
ELSE
   reducedx = .TRUE.
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
!2. Allocate array
!-----------------
!

nzdim = SIZE(realloc,DIM=2)
novals = SIZE(realloc,DIM=3)
ALLOCATE (realglb(nostatsglb,nzdim,novals),STAT=errstat)
CALL error_alloc('realglb',3,(/nostatsglb,nzdim,novals/),kndrtype)
realglb = 0.0
IF (reduced) THEN
   ALLOCATE (realglb2(nostatsglb*novals,nzdim),STAT=errstat)
   CALL error_alloc('realglb',2,(/nostatsglb*novals,nzdim/),kndrtype)
   realglb2 = 0.0
ENDIF

!
!3. Combine data
!---------------
!

CALL combine_stats_loc(realglb,realloc,maxstats,nostatprocs,lstatprocs,ivarid,&
                     & idmaster,commop,itype)

!
!4. Reshape if needed
!--------------------
!

IF (reduced) THEN
   ind(1) = 1; ind(2) = nostatsglb
   l_410: DO l=1,novals
      realglb2(ind(1):ind(2),:) = realglb(:,:,l)
      ind = ind + nostatsglb
   ENDDO l_410
ENDIF

!
!5. Write data
!-------------
!

IF (master) THEN
   IF (reduced) THEN
      IF (PRESENT(varatts)) THEN
         CALL write_vars(realglb2,filepars,varid,varatts=varatts)
      ELSE
         CALL write_vars(realglb2,filepars,varid)
      ENDIF
   ELSE
      IF (PRESENT(varatts).AND.PRESENT(vecids)) THEN
         CALL write_vars(realglb,filepars,varid,varatts=varatts,vecids=vecids)
      ELSEIF (PRESENT(varatts).AND.(.NOT.PRESENT(vecids))) THEN
         CALL write_vars(realglb,filepars,varid,varatts=varatts)
      ELSEIF (.NOT.PRESENT(varatts).AND.PRESENT(vecids)) THEN
         CALL write_vars(realglb,filepars,varid,vecids=vecids)
      ELSEIF (.NOT.(PRESENT(varatts).OR.PRESENT(vecids))) THEN
         CALL write_vars(realglb,filepars,varid)
      ENDIF
   ENDIF
ENDIF

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (realglb)
IF (reduced) DEALLOCATE (realglb2)

CALL log_timer_out()


RETURN

END SUBROUTINE combine_write_stats_loc_real_3d

!========================================================================

SUBROUTINE combine_write_submod_real_2d(realloc,filepars,varid,ndimsglb,&
                                      & limprocs,nowet,varatts,comm,commtype)
!************************************************************************
!
! *combine_write_submod_real_2d* Write a global sub-model real 2-D array
!                                combined on the root process from local arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.10.1
!
! Description -
!
! Module calls - combine_submod, error_alloc, write_vars
!
!************************************************************************
!
USE datatypes
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: write_vars
USE paral_comms, ONLY: combine_submod

!
!*  Arguments
!
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, nowet
TYPE(FileParams), INTENT(INOUT) :: filepars
INTEGER, INTENT(IN), DIMENSION(2) :: ndimsglb
INTEGER, INTENT(IN), DIMENSION(2,2,nprocs) :: limprocs
REAL, INTENT(IN), DIMENSION(:,:) :: realloc
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type      Purpose
!------------------------------------------------------------------------------
!*realloc*   REAL     Local array on input
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*ndimsglb*  INTEGER  Shape of global array
!*limprocs*  INTEGER  Start/end indices of local array section in global array
!*nowet*     INTEGER  Number of no-wet horizontal data points in case the data
!                     are witten in "packed" (vector) format
!*varatts*   DERIVED  Variable attributes
!*comm*      INTEGER  Communicator for combine operation
!*commtype*  INTEGER  Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: fill, packing
INTEGER :: commop, itype, ivarid, nowetx
REAL :: fill_value, fill_value_eps
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: realdat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: realglb


IF (PRODUCT(ndimsglb).EQ.0) RETURN
procname(pglev+1) = 'combine_write_submod_real_2d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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

IF (PRESENT(nowet)) THEN
   nowetx = nowet
ELSE
   nowetx = 0
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---fill parameters
fill = filepars%fill
IF (fill) THEN
   fill_value = double_fill
   fill_value_eps = 0.00001*ABS(fill_value)
ELSE
   fill_value = 0.0
ENDIF

!---packing
packing = filepars%packing.AND.nowetx.GT.0.AND.master

!---allocate global data array
ALLOCATE (realglb(ndimsglb(1),ndimsglb(2)),STAT=errstat)
CALL error_alloc('realglb',2,ndimsglb,kndrtype)

!---data array for packing
IF (packing) THEN
   ALLOCATE (realdat(nowetx),STAT=errstat)
   CALL error_alloc('realdat',1,(/nowetx/),kndrtype)
ENDIF

!
!3. Combine data
!---------------
!

CALL combine_submod(realglb,realloc,limprocs,ivarid,fill_value,idmaster,commop,&
                  & itype)

!
!4. Pack data using land mask
!----------------------------
!

IF (packing) THEN
   realdat = PACK(realglb,MASK=ABS(realglb-fill_value).GT.fill_value_eps)
ENDIF

!
!5. Write data
!-------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realdat,filepars,varid,varatts=varatts)
   ELSE
      CALL write_vars(realdat,filepars,varid)
   ENDIF
ELSEIF (master) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realglb,filepars,varid,varatts=varatts)
   ELSE
      CALL write_vars(realglb,filepars,varid)
   ENDIF
ENDIF

!
!6. Deallocate array
!-------------------
!

DEALLOCATE (realglb)
IF (packing) DEALLOCATE (realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE combine_write_submod_real_2d

!========================================================================

SUBROUTINE combine_write_submod_real_3d(realloc,filepars,varid,ndimsglb,&
                                      & limprocs,nowet,varatts,vecids,comm,&
                                      & commtype)
!************************************************************************
!
! *combine_write_submod_real_3d* Write a global sub-model real 3-D array
!                                combined on the root process from local arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.10.1
!
! Description -
!
! Module calls - combine_submod, error_alloc, write_vars
!
!************************************************************************
!
USE datatypes
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: write_vars
USE paral_comms, ONLY: combine_submod

!
!*  Arguments
!
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, nowet
TYPE(FileParams), INTENT(INOUT) :: filepars
INTEGER, INTENT(IN), DIMENSION(3) :: ndimsglb
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
INTEGER, INTENT(IN), DIMENSION(2,2,nprocs) :: limprocs
REAL, INTENT(IN), DIMENSION(:,:,:) :: realloc
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type      Purpose
!------------------------------------------------------------------------------
!*realloc*   REAL     Local array on input
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*ndimsglb*  INTEGER  Shape of global array
!*limprocs*  INTEGER  Start/end indices of local array section in global array
!*nowet*     INTEGER  Number of no-wet horizontal data points in case the data
!                     are witten in "packed" (vector) format
!*varatts*   DERIVED  Variable attributes
!*vecids*    INTEGER  Arrays of variable IDs if varid = 0

!*comm*      INTEGER  Communicator for combine operation
!*commtype*  INTEGER  Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: fill, packing
INTEGER :: commop, itype, ivar, ivarid, k, nowetx
REAL :: fill_value, fill_value_eps
INTEGER, DIMENSION(ndimsglb(3)) :: idvars 
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: realdat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realglb


IF (PRODUCT(ndimsglb).EQ.0) RETURN
procname(pglev+1) = 'combine_write_submod_real_3d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSE
   idvars = (/(ivar,ivar=1,ndimsglb(3))/)
ENDIF

IF (PRESENT(nowet)) THEN
   nowetx = nowet
ELSE
   nowetx = 0
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---fill parameters
fill = filepars%fill
IF (fill) THEN
   fill_value = double_fill
   fill_value_eps = 0.00001*ABS(fill_value)
ELSE
   fill_value = 0.0
ENDIF

!---packing
packing = filepars%packing.AND.nowetx.GT.0.AND.master

!---allocate global data array
ALLOCATE (realglb(ndimsglb(1),ndimsglb(2),ndimsglb(3)),STAT=errstat)
CALL error_alloc('realglb',3,ndimsglb,kndrtype)

!---data array for packing
IF (packing) THEN
   ALLOCATE (realdat(nowetx,ndimsglb(3)),STAT=errstat)
   CALL error_alloc('realdat',2,(/nowetx,ndimsglb(3)/),kndrtype)
ENDIF

!
!3. Combine data
!---------------
!

CALL combine_submod(realglb,realloc,limprocs,ivarid,fill_value,idmaster,commop,&
                  & itype)

!
!4. Pack data using land mask
!----------------------------
!

IF (packing) THEN
   k_410: DO k=1,ndimsglb(3)
      realdat(:,k) = PACK(realglb(:,:,k),&
                   & MASK=ABS(realglb(:,:,k)-fill_value).GT.fill_value_eps)
   ENDDO k_410
ENDIF

!
!5. Write data
!-------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realdat,filepars,varid,varatts=varatts,vecids=idvars)
   ELSE
      CALL write_vars(realdat,filepars,varid,vecids=idvars)
   ENDIF
ELSEIF (master) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realglb,filepars,varid,varatts=varatts,vecids=idvars)
   ELSE
      CALL write_vars(realglb,filepars,varid,vecids=idvars)
   ENDIF
ENDIF

!
!6. Deallocate array
!-------------------
!

DEALLOCATE (realglb)
IF (packing) DEALLOCATE (realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE combine_write_submod_real_3d

!========================================================================

SUBROUTINE combine_write_submod_real_4d(realloc,filepars,varid,ndimsglb,&
                                      & limprocs,nowet,varatts,vecids,comm,&
                                      & commtype)
!************************************************************************
!
! *combine_write_submod_real_4d* Write a global sub-model real 4-D array
!                                combined on the root process from local arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.10.1
!
! Description -
!
! Module calls - combine_submod, error_alloc, write_vars
!
!************************************************************************
!
USE datatypes
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: write_vars
USE paral_comms, ONLY: combine_submod

!
!*  Arguments
!
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype, nowet
TYPE(FileParams), INTENT(INOUT) :: filepars
INTEGER, INTENT(IN), DIMENSION(4) :: ndimsglb
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
INTEGER, INTENT(IN), DIMENSION(2,2,nprocs) :: limprocs
REAL, INTENT(IN), DIMENSION(:,:,:,:) :: realloc
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type      Purpose
!------------------------------------------------------------------------------
!*realloc*   REAL     Local array on input
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*ndimsglb*  INTEGER  Shape of global array
!*limprocs*  INTEGER  Start/end indices of local array section in global array
!*nowet*     INTEGER  Number of no-wet horizontal data points in case the data
!                     are written in "packed" (vector) format
!*varatts*   DERIVED  Variable attributes
!*vecids*    INTEGER  Arrays of variable IDs if varid = 0
!*comm*      INTEGER  Communicator for combine operation
!*commtype*  INTEGER  Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: fill, packing
INTEGER :: commop, itype, ivar, ivarid, k, l, nowetx
REAL :: fill_value, fill_value_eps
INTEGER, DIMENSION(ndimsglb(4)) :: idvars
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realdat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: realglb


IF (PRODUCT(ndimsglb).EQ.0) RETURN
procname(pglev+1) = 'combine_write_submod_real_4d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSE
   idvars = (/(ivar,ivar=1,ndimsglb(4))/)
ENDIF

IF (PRESENT(nowet)) THEN
   nowetx = nowet
ELSE
   nowetx = 0
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---fill parameters
fill = filepars%fill
IF (fill) THEN
   fill_value = double_fill
   fill_value_eps = 0.00001*ABS(fill_value)
ELSE
   fill_value = 0.0
ENDIF

!---packing
packing = filepars%packing.AND.nowetx.GT.0.AND.master

!---allocate global data array
ALLOCATE (realglb(ndimsglb(1),ndimsglb(2),ndimsglb(3),ndimsglb(4)),&
        & STAT=errstat)
CALL error_alloc('realglb',4,ndimsglb,kndrtype)

!---data array for packing
IF (packing) THEN
   ALLOCATE (realdat(nowetx,ndimsglb(3),ndimsglb(4)),STAT=errstat)
   CALL error_alloc('realdat',3,(/nowetx,ndimsglb(3),ndimsglb(4)/),kndrtype)
ENDIF

!
!3. Combine data
!---------------
!

CALL combine_submod(realglb,realloc,limprocs,ivarid,fill_value,idmaster,commop,&
                  & itype)

!
!4. Pack data using land mask
!----------------------------
!

IF (packing) THEN
   l_410: DO l=1,ndimsglb(4)
   k_410: DO k=1,ndimsglb(3)
      realdat(:,k,l) = PACK(realglb(:,:,k,l),&
                     & MASK=ABS(realglb(:,:,k,l)-fill_value).GT.fill_value_eps)
   ENDDO k_410
   ENDDO l_410
ENDIF

!
!5. Write data
!-------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realdat,filepars,varid,varatts=varatts,vecids=idvars)
   ELSE
      CALL write_vars(realdat,filepars,varid,vecids=idvars)
   ENDIF
ELSEIF (master) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realglb,filepars,varid,varatts=varatts,vecids=idvars)
   ELSE
      CALL write_vars(realglb,filepars,varid,vecids=idvars)
   ENDIF
ENDIF

!
!6. Deallocate array
!-------------------
!

DEALLOCATE (realglb)
IF (packing) DEALLOCATE (realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE combine_write_submod_real_4d

!========================================================================

SUBROUTINE read_distribute_mod_real_2d(realloc,filepars,varid,lbounds,ubounds,&
                                     & nhdims,cnode,fdist,varatts,comm,commtype)
!************************************************************************
!
! *read_distribute_mod_real_2d* Read/distribute a global 2-D real model array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.10.1
!
! Description -
!
! Module calls - distribute_mod, error_alloc, global_mask, read_vars
!
!************************************************************************
!
USE datatypes
USE grid  
USE gridpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: global_mask
USE inout_routines, ONLY: read_vars
USE paral_comms, ONLY: distribute_mod

!
!*  Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
LOGICAL, INTENT(IN) :: fdist
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
INTEGER, INTENT(IN), DIMENSION(4) :: nhdims
TYPE(FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds, ubounds
REAL, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):) :: realloc
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type      Purpose
!------------------------------------------------------------------------------
!*realloc*   REAL     Distributed local array
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*lbounds*   INTEGER  Lower bounds of local/global array
!*ubounds*   INTEGER  Lower bounds of global array
!*nhdims*    INTEGER  Halo sizes of data array
!*cnode*     CHAR     Nodal type of the data on the model grid
!                     ('C','W','U','V')
!*fdist*     LOGICAL  Distribute data if .TRUE.
!*varatts*   DERIVED  Variable attributes
!*comm*      INTEGER  Communicator for distribute operation
!*commtype*  INTEGER  Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: packing
INTEGER :: commop, itype, ivarid, nowet
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskglb
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: realdat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: realglb


procname(pglev+1) = 'read_distribute_mod_real_2d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---packing
packing = filepars%packing

!---number of wet points
IF (packing) THEN
   SELECT CASE (TRIM(cnode))
      CASE ('C','W')
         nowet = noseaatc
      CASE ('U')
         nowet = noseaatu
      CASE ('V')
         nowet = noseaatv
   END SELECT
ENDIF

!---allocate global data array
ALLOCATE (realglb(lbounds(1):ubounds(1),lbounds(2):ubounds(2)),STAT=errstat)
CALL error_alloc('realglb',2,ubounds-lbounds+1,kndrtype)
realglb = 0.0

!---allocate global arrays for packing
IF (packing) THEN 
   ALLOCATE (maskglb(nc,nr),STAT=errstat)
   CALL error_alloc('maskglb',2,(/nc,nr/),kndlog)
   ALLOCATE (realdat(nowet),STAT=errstat)
   CALL error_alloc('realdat',1,(/nowet/),kndrtype)
ENDIF
   
!
!3. Global mask
!--------------
!
 
IF (packing) CALL global_mask(maskglb,cnode,1)

!
!4. Read data
!------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL read_vars(realdat,filepars,varid,varatts=varatts)
   ELSE
      CALL read_vars(realdat,filepars,varid)
   ENDIF
ELSE
   IF (PRESENT(varatts)) THEN
      CALL read_vars(realglb(1:nc,1:nr),filepars,varid,varatts=varatts)
   ELSE
      CALL read_vars(realglb(1:nc,1:nr),filepars,varid)
   ENDIF
ENDIF

!
!5. Unpack data using land mask
!------------------------------
!

IF (packing) realglb(1:nc,1:nr) = UNPACK(realdat,MASK=maskglb,FIELD=0.0)

!
!6. Distribute data
!------------------
!

IF (fdist) THEN
   CALL distribute_mod(realglb,realloc,lbounds,nhdims,ivarid,0.0,0,&
                    & .TRUE.,idmaster,commop,itype)
ENDIF

!
!7. Zero values at land cells/nodes
!----------------------------------
!

SELECT CASE (TRIM(cnode))
   CASE ('C','W')
      WHERE (.NOT.maskatc_int)
         realloc(1:ncloc,1:nrloc) = 0.0
      END WHERE
   CASE ('U')
      WHERE (node2du(1:ncloc,1:nrloc).LE.1)
         realloc(1:ncloc,1:nrloc) = 0.0
      END WHERE
   CASE ('V')
      WHERE (node2dv(1:ncloc,1:nrloc).LE.1)
         realloc(1:ncloc,1:nrloc) = 0.0
      END WHERE
END SELECT

!
!8. Deallocate array
!-------------------
!

DEALLOCATE (realglb)
IF (packing) DEALLOCATE (maskglb,realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE read_distribute_mod_real_2d

!========================================================================

SUBROUTINE read_distribute_mod_real_3d(realloc,filepars,varid,lbounds,ubounds,&
                                     & nhdims,cnode,fdist,varatts,vecids,comm,&
                                     & commtype)
!************************************************************************
!
! *read_distribute_mod_real_3d* Read/distribute a global 3-D real model array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.10.1
!
! Description -
!
! Module calls - distribute_mod, error_alloc, global_mask, read_vars
!
!************************************************************************
!
USE datatypes
USE grid
USE gridpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: global_mask
USE inout_routines, ONLY: read_vars
USE paral_comms, ONLY: distribute_mod

!
!*  Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
LOGICAL, INTENT(IN) :: fdist
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
TYPE(FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
INTEGER, INTENT(IN), DIMENSION(4) :: nhdims
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):,lbounds(3):) :: realloc
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type      Purpose
!------------------------------------------------------------------------------
!*realloc*   REAL     Distributed local array
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*lbounds*   INTEGER  Lower bounds of local/global array
!*ubounds*   INTEGER  Lower bounds of global array
!*nhdims*    INTEGER  Halo sizes of data array
!*cnode*     CHAR     Nodal type of the data on the model grid
!                     ('C','W','U','V')
!*fdist*     LOGICAL  Distribute data if .TRUE.
!*varatts*   DERIVED  Variable attributes
!*vecids*    INTEGER  Arrays of variable IDs if varid = 0
!*comm*      INTEGER  Communicator for distribute operation
!*commtype*  INTEGER  Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: packing
INTEGER :: commop, itype, ivar, ivarid, k, nowet
INTEGER, DIMENSION(SIZE(realloc,DIM=3)) :: idvars
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskglb
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: realdat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realglb


procname(pglev+1) = 'read_distribute_mod_real_3d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSE
   idvars = (/(ivar,ivar=1,SIZE(realloc,DIM=3))/)
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---packing
packing = filepars%packing

!---number of wet points
IF (packing) THEN
   SELECT CASE (TRIM(cnode))
      CASE ('C','W')
         nowet = noseaatc
      CASE ('U')
         nowet = noseaatu
      CASE ('V')
         nowet = noseaatv
   END SELECT
ENDIF

!---allocate global data array
ALLOCATE (realglb(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                & lbounds(3):ubounds(3)),STAT=errstat)
CALL error_alloc('realglb',3,ubounds-lbounds+1,kndrtype)
realglb = 0.0

!---allocate global arrays for packing
IF (packing) THEN
   ALLOCATE (maskglb(nc,nr),STAT=errstat)
   CALL error_alloc('maskglb',2,(/nc,nr/),kndlog)
   ALLOCATE (realdat(nowet,lbounds(3):ubounds(3)),STAT=errstat)
   CALL error_alloc('realdat',2,(/nowet,ubounds(3)-lbounds(3)+1/),kndrtype)
ENDIF

!
!3. Global mask
!--------------
!
 
IF (packing) CALL global_mask(maskglb,cnode,1)

!
!4. Read data
!------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL read_vars(realdat,filepars,varid,varatts=varatts,vecids=idvars)
   ELSE
      CALL read_vars(realdat,filepars,varid,vecids=idvars)
   ENDIF
ELSE
   IF (PRESENT(varatts)) THEN
      CALL read_vars(realglb(1:nc,1:nr,:),filepars,varid,varatts=varatts,&
                   & vecids=idvars)
   ELSE
      CALL read_vars(realglb(1:nc,1:nr,:),filepars,varid,vecids=idvars)
   ENDIF
ENDIF

!
!5. Unpack data using land mask
!------------------------------
!

IF (packing) THEN
   k_510: DO k=lbounds(3),ubounds(3)
      realglb(1:nc,1:nr,k) = UNPACK(realdat(:,k),MASK=maskglb,FIELD=0.0)
   ENDDO k_510
ENDIF

!
!6. Distribute data
!------------------
!

IF (fdist) THEN
   CALL distribute_mod(realglb,realloc,lbounds,nhdims,ivarid,0.0,0,.TRUE.,&
                     & idmaster,commop,itype)
ENDIF

!
!7. Zero values at land cells/nodes
!----------------------------------
!

k_710: DO k=lbounds(3),ubounds(3)
   SELECT CASE (TRIM(cnode))
      CASE ('C','W')
         WHERE (.NOT.maskatc_int)
            realloc(1:ncloc,1:nrloc,k) = 0.0
         END WHERE
      CASE ('U')
         WHERE (node2du(1:ncloc,1:nrloc).LE.1)
            realloc(1:ncloc,1:nrloc,k) = 0.0
         END WHERE
      CASE ('V')
         WHERE (node2dv(1:ncloc,1:nrloc).LE.1)
            realloc(1:ncloc,1:nrloc,k) = 0.0
         END WHERE
    END SELECT
ENDDO k_710
 
!
!8. Deallocate array
!-------------------
!

DEALLOCATE (realglb)
IF (packing) DEALLOCATE (maskglb,realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE read_distribute_mod_real_3d

!========================================================================

SUBROUTINE read_distribute_mod_real_4d(realloc,filepars,varid,lbounds,ubounds,&
                                     & nhdims,cnode,fdist,varatts,vecids,comm,&
                                     & commtype)
!************************************************************************
!
! *read_distribute_mod_real_4d* Read/distribute a global 4-D real model array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.10.1
!
! Description -
!
! Module calls - distribute_mod, error_alloc, global_mask, read_vars
!
!************************************************************************
!
USE datatypes
USE grid
USE gridpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: global_mask
USE inout_routines, ONLY: read_vars
USE paral_comms, ONLY: distribute_mod

!
!*  Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
LOGICAL, INTENT(IN) :: fdist
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
TYPE(FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), DIMENSION(4) :: lbounds, nhdims, ubounds
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):,&
                             & lbounds(3):,lbounds(4):) :: realloc
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type      Purpose
!------------------------------------------------------------------------------
!*realloc*   REAL     Distributed local array
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*lbounds*   INTEGER  Lower bounds of local/global array
!*ubounds*   INTEGER  Lower bounds of global array
!*nhdims*    INTEGER  Halo sizes of data array
!*cnode*     CHAR     Nodal type of the data on the model grid
!                     ('C','W','U','V')
!*fdist*     LOGICAL  Distribute data if .TRUE.
!*varatts*   DERIVED  Variable attributes
!*vecids*    INTEGER  Arrays of variable IDs if varid = 0
!*comm*      INTEGER  Communicator for distribute operation
!*commtype*  INTEGER  Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: packing
INTEGER :: commop, itype, ivar, ivarid, k, l, nowet
INTEGER, DIMENSION(SIZE(realloc,DIM=4)) :: idvars
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskglb
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realdat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: realglb


procname(pglev+1) = 'read_distribute_mod_real_4d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSE
   idvars = (/(ivar,ivar=1,SIZE(realloc,DIM=4))/)
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---packing
packing = filepars%packing

!---number of wet points
IF (packing) THEN
   SELECT CASE (TRIM(cnode))
      CASE ('C','W')
         nowet = noseaatc
      CASE ('U')
         nowet = noseaatu
      CASE ('V')
         nowet = noseaatv
   END SELECT
ENDIF

!---allocate global data array
ALLOCATE (realglb(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                & lbounds(3):ubounds(3),lbounds(4):ubounds(4)),STAT=errstat)
CALL error_alloc('realglb',4,ubounds-lbounds+1,kndrtype)
realglb = 0.0

!---allocate global arrays for packing
IF (packing) THEN
   ALLOCATE (maskglb(nc,nr),STAT=errstat)
   CALL error_alloc('maskglb',2,(/nc,nr/),kndlog)
   ALLOCATE (realdat(nowet,lbounds(3):ubounds(3),lbounds(4):ubounds(4)),&
           & STAT=errstat)
   CALL error_alloc('realdat',4,(/nowet,ubounds(3:4)-lbounds(3:4)+1/),kndrtype)
ENDIF

!
!3. Initialise global array
!--------------------------
!

IF (packing) CALL global_mask(maskglb,cnode,1)

!
!4. Read data
!------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL read_vars(realdat,filepars,varid,varatts=varatts,vecids=idvars)
   ELSE
      CALL read_vars(realdat,filepars,varid,vecids=idvars)
   ENDIF
ELSE
   IF (PRESENT(varatts)) THEN
      CALL read_vars(realglb(1:nc,1:nr,:,:),filepars,varid,varatts=varatts,&
                   & vecids=idvars)
   ELSE
      CALL read_vars(realglb(1:nc,1:nr,:,:),filepars,varid,vecids=idvars)
   ENDIF
ENDIF

!
!5. Unpack data using land mask
!------------------------------
!

IF (packing) THEN
   l_510: DO l=lbounds(4),ubounds(4)
   k_510: DO k=lbounds(3),ubounds(3)
      realglb(1:nc,1:nr,k,l) = UNPACK(realdat(:,k,l),MASK=maskglb,FIELD=0.0)
   ENDDO k_510
   ENDDO l_510
ENDIF

!
!6. Distribute data
!------------------
!

IF (fdist) THEN
   CALL distribute_mod(realglb,realloc,lbounds,nhdims,ivarid,0.0,0,.TRUE.,&
                     & idmaster,commop,itype)
ENDIF

!
!7. Zero values at land cells/nodes
!----------------------------------
!

l_710: DO l=lbounds(4),ubounds(4)
k_710: DO k=lbounds(3),ubounds(3)
   SELECT CASE (TRIM(cnode))
      CASE ('C','W')
         WHERE (.NOT.maskatc_int)
            realloc(1:ncloc,1:nrloc,k,l) = 0.0
         END WHERE
      CASE ('U')
         WHERE (node2du(1:ncloc,1:nrloc).LE.1)
            realloc(1:ncloc,1:nrloc,k,l) = 0.0
         END WHERE
      CASE ('V')
         WHERE (node2dv(1:ncloc,1:nrloc).LE.1)
            realloc(1:ncloc,1:nrloc,k,l) = 0.0
         END WHERE
    END SELECT
ENDDO k_710
ENDDO l_710

!
!8. Deallocate array
!-------------------
!

DEALLOCATE (realglb)
IF (packing) DEALLOCATE (maskglb,realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE read_distribute_mod_real_4d

!========================================================================

SUBROUTINE read_distribute_mod_log_2d(logloc,filepars,varid,lbounds,ubounds,&
                                    & nhdims,cnode,fdist,varatts,comm,commtype)
!************************************************************************
!
! *read_distribute_mod_log_2d* Read/distribute a global 2-D logical model array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.11
!
! Description -
!
! Module calls - distribute_mod, error_alloc, global_mask, read_vars
!
!************************************************************************
!
USE datatypes
USE grid  
USE gridpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: global_mask
USE inout_routines, ONLY: read_vars
USE paral_comms, ONLY: distribute_mod

!
!*  Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
LOGICAL, INTENT(IN) :: fdist
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: comm, commtype
INTEGER, INTENT(IN), DIMENSION(4) :: nhdims
TYPE(FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds, ubounds
LOGICAL, INTENT(INOUT), DIMENSION(lbounds(1):,lbounds(2):) :: logloc
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*logloc*    LOGICAL  Distributed local array
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*lbounds*   INTEGER  Lower bounds of local/global array
!*ubounds*   INTEGER  Lower bounds of global array
!*nhdims*    INTEGER  Halo sizes of data array
!*cnode*     CHAR     Nodal type of the data on the model grid
!                     ('C','W','U','V')
!*fdist*     LOGICAL  Distribute data if .TRUE.
!*varatts*   DERIVED  Variable attributes
!*comm*      INTEGER  Communicator for distribute operation
!*commtype*  INTEGER  Type of communication (1/4)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: packing
INTEGER :: commop, itype, ivarid, nowet
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskglb
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:) :: logdat
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: logglb


procname(pglev+1) = 'read_distribute_mod_log_2d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

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

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---packing
packing = filepars%packing

!---number of wet points
IF (packing) THEN
   SELECT CASE (TRIM(cnode))
      CASE ('C','W')
         nowet = noseaatc
      CASE ('U')
         nowet = noseaatu
      CASE ('V')
         nowet = noseaatv
   END SELECT
ENDIF

!---allocate global data array
ALLOCATE (logglb(lbounds(1):ubounds(1),lbounds(2):ubounds(2)),STAT=errstat)
CALL error_alloc('logglb',2,ubounds-lbounds+1,kndlog)
logglb = .FALSE.

!---allocate global arrays for packing
IF (packing) THEN 
   ALLOCATE (maskglb(nc,nr),STAT=errstat)
   CALL error_alloc('maskglb',2,(/nc,nr/),kndlog)
   ALLOCATE (logdat(nowet),STAT=errstat)
   CALL error_alloc('logdat',1,(/nowet/),kndlog)
ENDIF
   
!
!3. Global mask
!--------------
!
 
IF (packing) CALL global_mask(maskglb,cnode,1)

!
!4. Read data
!------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL read_vars(logdat,filepars,varid,varatts=varatts)
   ELSE
      CALL read_vars(logdat,filepars,varid)
   ENDIF
ELSE
   IF (PRESENT(varatts)) THEN
      CALL read_vars(logglb(1:nc,1:nr),filepars,varid,varatts=varatts)
   ELSE
      CALL read_vars(logglb(1:nc,1:nr),filepars,varid)
   ENDIF
ENDIF

!
!5. Unpack data using land mask
!------------------------------
!

IF (packing) logglb(1:nc,1:nr) = UNPACK(logdat,MASK=maskglb,FIELD=.FALSE.)

!
!6. Distribute data
!------------------
!

IF (fdist) THEN
   CALL distribute_mod(logglb,logloc,lbounds,nhdims,ivarid,.FALSE.,0,&
                    & .TRUE.,idmaster,commop,itype)
ENDIF

!
!7. Zero values at land cells/nodes
!----------------------------------
!

SELECT CASE (TRIM(cnode))
   CASE ('C','W')
      WHERE (.NOT.maskatc_int)
         logloc(1:ncloc,1:nrloc) = .FALSE.
      END WHERE
   CASE ('U')
      WHERE (node2du(1:ncloc,1:nrloc).LE.1)
         logloc(1:ncloc,1:nrloc) = .FALSE.
      END WHERE
   CASE ('V')
      WHERE (node2dv(1:ncloc,1:nrloc).LE.1)
         logloc(1:ncloc,1:nrloc) = .FALSE.
      END WHERE
END SELECT

!
!8. Deallocate array
!-------------------
!

DEALLOCATE (logglb)
IF (packing) DEALLOCATE (maskglb,logdat)

CALL log_timer_out()


RETURN

END SUBROUTINE read_distribute_mod_log_2d


END MODULE inout_paral
