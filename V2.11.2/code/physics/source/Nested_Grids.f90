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
! *Nested_Grids* Write open boundary conditions for nested grids 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.11
!
! $Date: 2018-06-05 16:03:05 +0200 (Tue, 05 Jun 2018) $
!
! $Revision: 1142 $
!
! Description - 
!
! Reference -
!
! Routines - define_nstgrd_locs, define_nstgrd_spec, nest_locations,
!            read_nstgrd_abs, read_nstgrd_rel, read_nstgrd_spec,
!            update_nest_data_2d, update_nest_data_3d, update_nest_data_prof,
!            write_nest_data_2d, write_nest_data_3d, write_nest_data_prof,
!            write_nstgrd_abs, write_nstgrd_rel, write_nstgrd_spec
!
!************************************************************************
!

!========================================================================

SUBROUTINE define_nstgrd_locs(iset,nhdat,nzdat,nonodes,hnests,zcoord,cnode)
!************************************************************************
!
! *define_nstgrd_locs* Geographical positions of nested open boundary locations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.8
!
! Description - 
!
! Reference -
!
! Calling program - nest_locations
!
! External calls - read_nstgrd_abs, read_nstgrd_rel, usrdef_nstgrd_abs,
!                  usrdef_nstgrd_rel, write_nstgrd_abs, write_nstgrd_rel
!
! Module calls - data_to_model_hcoords, error_alloc
!
!************************************************************************
!
USE datatypes
USE iopars
USE nestgrids
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_interp, ONLY: data_to_model_hcoords
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: iset, nhdat, nonodes, nzdat
REAL, INTENT(OUT), DIMENSION(nhdat,nzdat) :: zcoord
TYPE (HRelativeCoords), INTENT(OUT), DIMENSION(nhdat,nonodes) :: hnests

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iset*      INTEGER Number of output set
!*nhdat*     INTEGER Number of horizontal data locations
!*nzdat*     INTEGER Number of vertical data points
!*nonodes*   INTEGER Number of nodes for interpolation
!*hnests*    DERIVED Relative coordinates of nest locations
!*zcoord*    REAL    Z-coordinates of nest locations
!*cnode*     CHAR    Nodal type of nest locations
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, PARAMETER :: intflag = 1
LOGICAL :: flag
REAL, ALLOCATABLE, DIMENSION(:) :: xcoord, ycoord


!
!1. Initialise
!-------------
!

SELECT CASE (TRIM(cnode))
   CASE ('C'); flag = nohnstglbc(iset).GT.0
   CASE ('U'); flag = nohnstglbu(iset).GT.0
   CASE ('V'); flag = nohnstglbv(iset).GT.0
   CASE ('X'); flag = nohnstglbx(iset).GT.0
   CASE ('Y'); flag = nohnstglby(iset).GT.0
END SELECT
IF (.NOT.flag) RETURN

procname(pglev+1) = 'define_nstgrd_locs: '//TRIM(cnode)
CALL log_timer_in()

!
!2. Absolute coordinates
!-----------------------
!

IF (nestcoords(iset).EQ.1) THEN

!
!2.1 Allocate
!------------
!

   ALLOCATE (xcoord(nhdat),STAT=errstat)
   CALL error_alloc('xcoord',1,(/nhdat/),kndrtype)
   ALLOCATE (ycoord(nhdat),STAT=errstat)
   CALL error_alloc('ycoord',1,(/nhdat/),kndrtype)

!
!2.2 Define
!----------
!

!  ---read
   IF (modfiles(io_nstabs,iset,1)%status.EQ.'R') THEN
      CALL read_nstgrd_abs(iset,nhdat,nzdat,xcoord,ycoord,zcoord)
!  ---user-defined
   ELSEIF (modfiles(io_nstabs,iset,1)%status.EQ.'N') THEN
      CALL usrdef_nstgrd_abs(iset,nhdat,nzdat,xcoord,ycoord,zcoord,cnode)
   ENDIF
!  ---write locations
   IF (modfiles(io_nstabs,iset,2)%status.EQ.'W') THEN
      CALL write_nstgrd_abs(iset,nhdat,nzdat,xcoord,ycoord,zcoord)
   ENDIF

!
!2.3 Convert to relative coordinates
!-----------------------------------
!

   CALL data_to_model_hcoords(hnests(:,1),cnode,nhdat,1,xcoord,ycoord)
   IF (TRIM(cnode).EQ.'U'.OR.TRIM(cnode).EQ.'V') THEN
      CALL data_to_model_hcoords(hnests(:,2),'C  ',nhdat,1,xcoord,ycoord)
   ENDIF

!   
!2.4 Write relative coordinates
!------------------------------
!
   
   IF (modfiles(io_nstrel,iset,2)%status.EQ.'W') THEN
      CALL write_nstgrd_rel(iset,nhdat,nzdat,nonodes,hnests,zcoord)
   ENDIF

!
!2.5 Deallocate
!--------------
!

   DEALLOCATE (xcoord,ycoord)

!
!3. Relative coordinates
!-----------------------
!

ELSEIF (nestcoords(iset).EQ.2) THEN

!  ---read
   IF (modfiles(io_nstrel,iset,1)%status.EQ.'R') THEN
      CALL read_nstgrd_rel(iset,nhdat,nzdat,nonodes,hnests,zcoord)
!  ---user-defined
   ELSEIF (modfiles(io_nstrel,iset,1)%status.EQ.'N') THEN
      CALL usrdef_nstgrd_rel(iset,nhdat,nzdat,nonodes,hnests,zcoord,cnode)
   ENDIF
!  ---write
   IF (modfiles(io_nstrel,iset,2)%status.EQ.'W') THEN
      CALL write_nstgrd_rel(iset,nhdat,nzdat,nonodes,hnests,zcoord)
   ENDIF

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE define_nstgrd_locs

!========================================================================

SUBROUTINE define_nstgrd_spec
!************************************************************************
!
! *define_nstgrd_spec* Number of output locations and specifier arrays for
!                      nesting
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.11.1
!
! Description - 
!
! Reference -
!
! Calling program - nest_locations
!
! External calls - read_nstgrd_spec, usrdef_nstgrd_spec, write_nstgrd_spec
!
! Module calls - error_alloc, error_limits_arr, error_vals_arr_struc_arr,
!                warning_reset_arr, warning_reset_arr_struc
!
!************************************************************************
!
USE gridpars
USE iopars
USE nestgrids
USE paralpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc, error_limits_arr, error_vals_arr_struc_char,&
                        & warning_reset_arr, warning_reset_arr_struc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: ireset, iset


procname(pglev+1) = 'define_nstgrd_spec'
CALL log_timer_in()

!
!1. Allocate/initialise specifier arrays
!---------------------------------------
!
!---type of grids
ALLOCATE (nestcoords(nonestsets),STAT=errstat)
CALL error_alloc('nestcoords',1,(/nonestsets/),kndint)

!---global number of horizontal nest points
ALLOCATE (nohnstglbc(nonestsets),STAT=errstat)
CALL error_alloc('nohnstglbc',1,(/nonestsets/),kndint)
nohnstglbc = 0
ALLOCATE (nohnstglbu(nonestsets),STAT=errstat)
CALL error_alloc('nohnstglbu',1,(/nonestsets/),kndint)
nohnstglbu = 0
ALLOCATE (nohnstglbv(nonestsets),STAT=errstat)
CALL error_alloc('nohnstglbv',1,(/nonestsets/),kndint)
nohnstglbv = 0
ALLOCATE (nohnstglbx(nonestsets),STAT=errstat)
CALL error_alloc('nohnstglbx',1,(/nonestsets/),kndint)
nohnstglbx = 0
ALLOCATE (nohnstglby(nonestsets),STAT=errstat)
CALL error_alloc('nohnstglby',1,(/nonestsets/),kndint)
nohnstglby = 0

!---number of vertical nest points per area
ALLOCATE (novnst(nonestsets),STAT=errstat)
CALL error_alloc('novnst',1,(/nonestsets/),kndint)
novnst = 0

!---local number of horizontal nest points per area
ALLOCATE (nohnstatc(nonestsets),STAT=errstat)
CALL error_alloc('nohnstatc',1,(/nonestsets/),kndint)
nohnstatc = 0
ALLOCATE (nohnstatu(nonestsets),STAT=errstat)
CALL error_alloc('nohnstatu',1,(/nonestsets/),kndint)
nohnstatu = 0
ALLOCATE (nohnstatv(nonestsets),STAT=errstat)
CALL error_alloc('nohnstatv',1,(/nonestsets/),kndint)
nohnstatv = 0
ALLOCATE (nohnstatx(nonestsets),STAT=errstat)
CALL error_alloc('nohnstatx',1,(/nonestsets/),kndint)
nohnstatx = 0
ALLOCATE (nohnstaty(nonestsets),STAT=errstat)
CALL error_alloc('nohnstaty',1,(/nonestsets/),kndint)
nohnstaty = 0

!---buffer index
ALLOCATE (lbhnstatc(nonestsets),STAT=errstat)
CALL error_alloc('lbhnstatc',1,(/nonestsets/),kndint)
lbhnstatc = 1
ALLOCATE (lbhnstatu(nonestsets),STAT=errstat)
CALL error_alloc('lbhnstatu',1,(/nonestsets/),kndint)
lbhnstatu = 1
ALLOCATE (lbhnstatv(nonestsets),STAT=errstat)
CALL error_alloc('lbhnstatv',1,(/nonestsets/),kndint)
lbhnstatv = 1
ALLOCATE (lbhnstatx(nonestsets),STAT=errstat)
CALL error_alloc('lbhnstatx',1,(/nonestsets/),kndint)
lbhnstatx = 1
ALLOCATE (lbhnstaty(nonestsets),STAT=errstat)
CALL error_alloc('lbhnstaty',1,(/nonestsets/),kndint)
lbhnstaty = 1

!---number of points per process and set
ALLOCATE (nohnstcprocs(nprocs,nonestsets),STAT=errstat)
CALL error_alloc('nohnstcprocs',2,(/nprocs,nonestsets/),kndint)
nohnstcprocs = 0
ALLOCATE (nohnstuvprocs(nprocs,nonestsets),STAT=errstat)
CALL error_alloc('nohnstuvprocs',2,(/nprocs,nonestsets/),kndint)
nohnstuvprocs = 0
ALLOCATE (nohnstxyprocs(nprocs,nonestsets),STAT=errstat)
CALL error_alloc('nohnstxyprocs',2,(/nprocs,nonestsets/),kndint)
nohnstxyprocs = 0

!---type of 2-D data
ALLOCATE (inst2dtype(nonestsets),STAT=errstat)
CALL error_alloc('inst2dtype',1,(/nonestsets/),kndint) 
inst2dtype = 0

!---sediment fractions used in nesting
IF (iopt_sed.GT.0) THEN
   ALLOCATE (nosednst(nonestsets),STAT=errstat)
   CALL error_alloc('nosednst',1,(/nonestsets/),kndint)
   ALLOCATE (instsed(maxsedvars,nonestsets),STAT=errstat)
   CALL error_alloc('instsed',2,(/maxsedvars,nonestsets/),kndint)
ENDIF

!---biological variables used in nesting
IF (iopt_biolgy.GT.0) THEN
   ALLOCATE (nobionst(nonestsets),STAT=errstat)
   CALL error_alloc('nobionst',1,(/nonestsets/),kndint)
   ALLOCATE (instbio(maxbiovars,nonestsets),STAT=errstat)
   CALL error_alloc('instbio',2,(/maxbiovars,nonestsets/),kndint)
ENDIF

!
!2. Obtain specifier arrays
!--------------------------
!
!---read
IF (modfiles(io_nstspc,1,1)%status.EQ.'R') THEN
   CALL read_nstgrd_spec
!---user-defined
ELSEIF (modfiles(io_nstspc,1,1)%status.EQ.'N') THEN
   CALL usrdef_nstgrd_spec
ENDIF
!---write
IF (modfiles(io_nstspc,1,2)%status.EQ.'W') CALL write_nstgrd_spec

!
!3. Check specifier arrays
!-------------------------
!

IF (master) THEN
   iset_310: DO iset=1,nonestsets
      CALL error_limits_arr(nestcoords(iset),'nestcoords',1,2,1,(/iset/))
      CALL error_limits_arr(inst2dtype(iset),'inst2dtype',1,3,1,(/iset/))
   ENDDO iset_310
ENDIF

!
!4. Check coordinate files
!-------------------------
!

iset_410: DO iset=1,nonestsets
   IF (nestcoords(iset).EQ.1) THEN
      CALL warning_reset_arr_struc(modfiles(io_nstrel,iset,1)%status,&
                                & 'modfiles','%status','0',3,&
                                & (/io_nstrel,iset,1/))
      CALL warning_reset_arr_struc(modfiles(io_nstrel,iset,2)%status,&
                                & 'modfiles','%status','0',3,&
                                & (/io_nstrel,iset,2/))
      CALL error_vals_arr_struc_char(modfiles(io_nstabs,iset,1)%status,&
                                  & 'modfiles','%status','"R" "N"',3,&
                                  & (/io_nstabs,iset,1/))
   ELSEIF (nestcoords(iset).EQ.2) THEN
      CALL warning_reset_arr_struc(modfiles(io_nstabs,iset,1)%status,&
                                & 'modfiles','%status','0',3,&
                                & (/io_nstabs,iset,1/))
      CALL warning_reset_arr_struc(modfiles(io_nstabs,iset,2)%status,&
                                & 'modfiles','%status','0',3,&
                                & (/io_nstabs,iset,2/))
      CALL error_vals_arr_struc_char(modfiles(io_nstrel,iset,1)%status,&
                                  & 'modfiles','%status','"R" "N"',3,&
                                  & (/io_nstrel,iset,1/))
   ENDIF
ENDDO iset_410

!
!5. Reset parameters if needed
!-----------------------------
!

iset_510: DO iset=1,nonestsets
   IF (nohnstglbc(iset).GT.0) THEN
      ireset = nohnstglbu(iset) + nohnstglbv(iset)
     CALL warning_reset_arr(nohnstglbc(iset),'nohnstglbc',ireset,1,(/iset/))
   ENDIF
ENDDO iset_510

CALL log_timer_out()


RETURN

END SUBROUTINE define_nstgrd_spec

!========================================================================

SUBROUTINE nest_locations
!************************************************************************
!
! *nest_locations* Define locations of nested open boundaries
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.7
!
! Description - initialise_model
!
! Reference -
!
! Calling program - 
!
! External calls - define_nstgrd_spec, define_nstgrd_locs
!
! Module calls - data_to_model_vcoords_3d,
!                error_alloc, error_alloc_struc, hrelativecoords_init,
!                local_proc, vrelativecoords_init
!
!************************************************************************
!
USE datatypes
USE gridpars
USE iopars
USE nestgrids
USE paralpars
USE syspars
USE datatypes_init, ONLY: hrelativecoords_init, vrelativecoords_init 
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE grid_interp, ONLY: data_to_model_vcoords_3d
USE grid_routines, ONLY: local_proc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, iproc, iset, j, lbuf, lglb, lloc, lsub, l1, l2, nhdat, nodat, &
         & npcc, numhnst, numhnstglb, numvnst, nzdat
INTEGER, ALLOCATABLE, DIMENSION(:) :: lbnstglb
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nohnstuprocs, nohnstvprocs
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nohnstxprocs, nohnstyprocs
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: indexnstu, indexnstv
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: indexnstx, indexnsty
REAL, ALLOCATABLE, DIMENSION(:,:) :: znstglb, znstloc
TYPE (HRelativeCoords), ALLOCATABLE, DIMENSION(:,:) :: hnstglb


procname(pglev+1) = 'nest_locations'
CALL log_timer_in(npcc)

!
!1. Allocate work space array
!----------------------------
!

ALLOCATE (lbnstglb(nonestsets),STAT=errstat)
CALL error_alloc('lbnstglb',1,(/nonestsets/),kndint)

!
!2. Number of output locations
!-----------------------------
!

CALL define_nstgrd_spec

!
!3. C-node points
!----------------
!
!3.1 Global buffer index
!-----------------------
!

lbuf = 1
iset_310: DO iset=1,nonestsets
   lbnstglb(iset) = lbuf
   lbuf = lbuf + nohnstglbc(iset)
ENDDO iset_310

!
!3.2 Allocate work space arrays
!------------------------------
!

numhnstglb = SUM(nohnstglbc); numvnst = MAXVAL(novnst)
ALLOCATE (hnstglb(numhnstglb,1),STAT=errstat)
CALL error_alloc_struc('hnstglb',2,(/numhnstglb,1/),'HRelativeCoords')
CALL hrelativecoords_init(hnstglb,.TRUE.)
ALLOCATE (znstglb(numhnstglb,numvnst),STAT=errstat)
CALL error_alloc('znstglb',2,(/numhnstglb,numvnst/),kndrtype)

!
!3.3 Define coordinate structure
!-------------------------------
!

iset_330: DO iset=1,nonestsets
   nhdat = nohnstglbc(iset)
   IF (nhdat.GT.0) THEN
      l1 = lbnstglb(iset); l2 = l1 + nhdat - 1
      nzdat = novnst(iset)
      CALL define_nstgrd_locs(iset,nhdat,nzdat,1,hnstglb(l1:l2,:),&
                            & znstglb(l1:l2,1:nzdat),'C  ')
   ENDIF
ENDDO iset_330

!
!3.4 Local number of data points and buffer index
!------------------------------------------------
!

iproc_340: DO iproc=1,nprocs
   lloc = 0
   lglb_341: DO lglb=1,numhnstglb
      i = hnstglb(lglb,1)%icoord; j = hnstglb(lglb,1)%jcoord
      IF (.NOT.local_proc(i,j,iproc=iproc).OR.(i.EQ.int_fill)) CYCLE lglb_341
      lloc = lloc + 1
      iset_3411: DO iset=1,nonestsets
         l1 = lbnstglb(iset); l2 = l1 + nohnstglbc(iset)
         IF (lglb.GE.l1.AND.lglb.LT.l2) THEN
            nohnstcprocs(iproc,iset) = nohnstcprocs(iproc,iset) + 1
            IF (idloc.EQ.idprocs(iproc)) THEN
               nohnstatc(iset) = nohnstatc(iset) + 1
               IF (nohnstatc(iset).EQ.1) lbhnstatc(iset) = lloc 
            ENDIF
         ENDIF
      ENDDO iset_3411
   ENDDO lglb_341
ENDDO iproc_340
   
!
!3.5 Index maps
!--------------
!
!3.5.1 Allocate
!--------------
!

nodat = MAXVAL(nohnstcprocs)
ALLOCATE (indexnstc(nodat,nprocs,nonestsets),STAT=errstat)
CALL error_alloc('indexnstc',3,(/nodat,nprocs,nonestsets/),kndint)
IF (nodat.GT.0) indexnstc = 0

!
!3.5.2 Define
!-----------
!

iproc_352: DO iproc=1,nprocs
iset_352: DO iset=1,nonestsets
   l1 = lbnstglb(iset); l2 = l1 + nohnstglbc(iset) - 1
   lsub = 0
   lglb_3521: DO lglb=l1,l2
      i = hnstglb(lglb,1)%icoord; j = hnstglb(lglb,1)%jcoord
      IF (local_proc(i,j,iproc=iproc).AND.(i.NE.int_fill)) THEN
         lsub = lsub + 1
         indexnstc(lsub,iproc,iset) = lglb - l1 + 1
      ENDIF
   ENDDO lglb_3521
ENDDO iset_352
ENDDO iproc_352

!
!3.6 Local coordinate structures
!-------------------------------
!
!3.6.1 Allocate/initialise
!-------------------------
!
!---allocate
numhnst = SUM(nohnstatc)
ALLOCATE (hnstctoc(numhnst),STAT=errstat)
CALL error_alloc_struc('hnstctoc',1,(/numhnst/),'HRelativeCoords')
ALLOCATE (vnstctoc(2,2,numhnst,numvnst),STAT=errstat)
CALL error_alloc_struc('vnstctoc',4,(/2,2,numhnst,numvnst/),'VRelativeCoords')
nodat = MAXVAL(nohnstatc)
ALLOCATE (znstloc(nodat,numvnst),STAT=errstat)
CALL error_alloc('znstloc',2,(/nodat,numvnst/),kndrtype)

!---initialise
CALL hrelativecoords_init(hnstctoc,.TRUE.)
CALL vrelativecoords_init(vnstctoc,.TRUE.)

!
!3.6.2 Define
!------------
!

iproc_362: DO iproc=1,nprocs

   IF (idloc.EQ.idprocs(iproc)) THEN

      iset_3621: DO iset=1,nonestsets
         nhdat = nohnstatc(iset); nzdat = novnst(iset)
         l1 = lbhnstatc(iset); l2 = l1 + nhdat - 1
         lsub = 0
         lloc_36211: DO lloc=l1,l2
            lsub = lsub + 1
            lglb = indexnstc(lsub,iproc,iset) + lbnstglb(iset) - 1
            IF (hnstglb(lglb,1)%icoord.NE.int_fill) THEN
               hnstctoc(lloc)%icoord = hnstglb(lglb,1)%icoord - nc1loc + 1
               hnstctoc(lloc)%jcoord = hnstglb(lglb,1)%jcoord - nr1loc + 1
               hnstctoc(lloc)%weights = hnstglb(lglb,1)%weights
            ENDIF
            IF (nzdat.GT.0) znstloc(lsub,1:nzdat) = znstglb(lglb,1:nzdat)
         ENDDO lloc_36211

         IF (nhdat*nzdat.GT.0) THEN
            CALL data_to_model_vcoords_3d(znstloc(1:nhdat,1:nzdat),&
                                & hnstctoc(l1:l2),vnstctoc(:,:,l1:l2,1:nzdat),&
                                & nhdat,1,nzdat,'C  ')
         ENDIF

      ENDDO iset_3621

   ENDIF

ENDDO iproc_362

!
!3.7 Deallocate
!--------------
!

DEALLOCATE (hnstglb,znstglb,znstloc)

!
!4. U-node points
!----------------
!
!4.1 Global buffer index
!-----------------------
!

lbuf = 1
iset_410: DO iset=1,nonestsets
   lbnstglb(iset) = lbuf
   lbuf = lbuf + nohnstglbu(iset)
ENDDO iset_410

!
!4.2 Allocate work space arrays
!------------------------------
!

numhnstglb = SUM(nohnstglbu)
ALLOCATE (hnstglb(numhnstglb,2),STAT=errstat)
CALL error_alloc_struc('hnstglb',2,(/numhnstglb,2/),'HRelativeCoords')
CALL hrelativecoords_init(hnstglb,.TRUE.)
ALLOCATE (znstglb(numhnstglb,numvnst),STAT=errstat)
CALL error_alloc('znstglb',2,(/numhnstglb,numvnst/),kndrtype)
ALLOCATE (nohnstuprocs(nprocs,nonestsets),STAT=errstat)
CALL error_alloc('nohnstuprocs',2,(/nprocs,nonestsets/),kndint)
nohnstuprocs = 0

!
!4.3 Define coordinate structure
!-------------------------------
!

iset_430: DO iset=1,nonestsets
   nhdat = nohnstglbu(iset)
   IF (nhdat.GT.0) THEN
      l1 = lbnstglb(iset); l2 = l1 + nhdat - 1
      nzdat = novnst(iset)
      CALL define_nstgrd_locs(iset,nhdat,nzdat,2,hnstglb(l1:l2,:),&
                            & znstglb(l1:l2,1:nzdat),'U  ')
   ENDIF
ENDDO iset_430

!
!4.4 Local number of data points and buffer index
!------------------------------------------------
!

iproc_440: DO iproc=1,nprocs
   lloc = 0
   lglb_441: DO lglb=1,numhnstglb
      i = hnstglb(lglb,1)%icoord; j = hnstglb(lglb,1)%jcoord
      IF (.NOT.local_proc(i,j,iproc=iproc).OR.(i.EQ.int_fill)) CYCLE lglb_441
      lloc = lloc + 1
      iset_4411: DO iset=1,nonestsets
         l1 = lbnstglb(iset); l2 = l1 + nohnstglbu(iset)
         IF (lglb.GE.l1.AND.lglb.LT.l2) THEN
            nohnstuprocs(iproc,iset) = nohnstuprocs(iproc,iset) + 1
            IF (idloc.EQ.idprocs(iproc)) THEN
               nohnstatu(iset) = nohnstatu(iset) + 1
               IF (nohnstatu(iset).EQ.1) lbhnstatu(iset) = lloc 
            ENDIF
         ENDIF
      ENDDO iset_4411
   ENDDO lglb_441
ENDDO iproc_440
   
!
!4.5 Index map (U-points)
!------------------------
!
!4.5.1 Allocate
!--------------
!

nodat = MAXVAL(nohnstuprocs)
ALLOCATE (indexnstu(nodat,nprocs,nonestsets),STAT=errstat)
CALL error_alloc('indexnstu',3,(/nodat,nprocs,nonestsets/),kndint)
IF (nodat.GT.0) indexnstu = 0

!
!4.5.2 Define
!-----------
!

iproc_452: DO iproc=1,nprocs
iset_452: DO iset=1,nonestsets
   l1 = lbnstglb(iset); l2 = l1 + nohnstglbu(iset) - 1
   lsub = 0
   lglb_4521: DO lglb=l1,l2
      i = hnstglb(lglb,1)%icoord; j = hnstglb(lglb,1)%jcoord
      IF (local_proc(i,j,iproc=iproc).AND.(i.NE.int_fill)) THEN
         lsub = lsub + 1
         indexnstu(lsub,iproc,iset) = lglb - l1 + 1
      ENDIF
   ENDDO lglb_4521
ENDDO iset_452
ENDDO iproc_452

!
!4.6 Local coordinate structures
!-------------------------------
!
!4.6.1 Allocate/initialise
!-------------------------
!
!---allocate
numhnst = SUM(nohnstatu)
ALLOCATE (hnstutou(numhnst),STAT=errstat)
CALL error_alloc_struc('hnstutou',1,(/numhnst/),'HRelativeCoords')
ALLOCATE (hnstctou(numhnst),STAT=errstat)
CALL error_alloc_struc('hnstctou',1,(/numhnst/),'HRelativeCoords')
ALLOCATE (vnstutou(2,2,numhnst,numvnst),STAT=errstat)
CALL error_alloc_struc('vnstutou',4,(/2,2,numhnst,numvnst/),'VRelativeCoords')
nodat = MAXVAL(nohnstatu)
ALLOCATE (znstloc(nodat,numvnst),STAT=errstat)
CALL error_alloc('znstloc',2,(/nodat,numvnst/),kndrtype)

!---initialise
CALL hrelativecoords_init(hnstutou,.TRUE.)
CALL hrelativecoords_init(hnstctou,.TRUE.)
CALL vrelativecoords_init(vnstutou,.TRUE.)

!
!4.6.2 Define
!------------
!

iproc_462: DO iproc=1,nprocs

   IF (idloc.EQ.idprocs(iproc)) THEN

      iset_4621: DO iset=1,nonestsets
         nhdat = nohnstatu(iset); nzdat = novnst(iset)
         l1 = lbhnstatu(iset); l2 = l1 + nhdat - 1
         lsub = 0
         lloc_46211: DO lloc=l1,l2
            lsub = lsub + 1
            lglb = indexnstu(lsub,iproc,iset) + lbnstglb(iset) - 1
            IF (hnstglb(lglb,1)%icoord.NE.int_fill) THEN
               hnstutou(lloc)%icoord = hnstglb(lglb,1)%icoord - nc1loc + 1
               hnstutou(lloc)%jcoord = hnstglb(lglb,1)%jcoord - nr1loc + 1
               hnstutou(lloc)%weights = hnstglb(lglb,1)%weights
               hnstctou(lloc)%icoord = hnstglb(lglb,2)%icoord - nc1loc + 1
               hnstctou(lloc)%jcoord = hnstglb(lglb,2)%jcoord - nr1loc + 1
               hnstctou(lloc)%weights = hnstglb(lglb,2)%weights
            ENDIF
            IF (nzdat.GT.0) znstloc(lsub,1:nzdat) = znstglb(lglb,1:nzdat)
         ENDDO lloc_46211

         IF (nhdat*nzdat.GT.0) THEN
            CALL data_to_model_vcoords_3d(znstloc(1:nhdat,1:nzdat),&
                                & hnstutou(l1:l2),vnstutou(:,:,l1:l2,1:nzdat),&
                                & nhdat,1,nzdat,'U  ')
         ENDIF

      ENDDO iset_4621

   ENDIF

ENDDO iproc_462

!
!4.7 Deallocate
!--------------
!

DEALLOCATE (hnstglb,znstglb,znstloc)

!
!5. V-node points
!----------------
!
!5.1 Global buffer index
!-----------------------
!

lbuf = 1
iset_510: DO iset=1,nonestsets
   lbnstglb(iset) = lbuf
   lbuf = lbuf + nohnstglbv(iset)
ENDDO iset_510

!
!5.2 Allocate work space arrays
!------------------------------
!

numhnstglb = SUM(nohnstglbv)
ALLOCATE (hnstglb(numhnstglb,2),STAT=errstat)
CALL error_alloc_struc('hnstglb',2,(/numhnstglb,2/),'HRelativeCoords')
CALL hrelativecoords_init(hnstglb,.TRUE.)
ALLOCATE (znstglb(numhnstglb,numvnst),STAT=errstat)
CALL error_alloc('znstglb',2,(/numhnstglb,numvnst/),kndrtype)
ALLOCATE (nohnstvprocs(nprocs,nonestsets),STAT=errstat)
CALL error_alloc('nohnstvprocs',2,(/nprocs,nonestsets/),kndint)
nohnstvprocs = 0

!
!5.3 Define coordinate structure
!-------------------------------
!

iset_530: DO iset=1,nonestsets
   nhdat = nohnstglbv(iset)
   IF (nhdat.GT.0) THEN
      l1 = lbnstglb(iset); l2 = l1 + nhdat - 1
      nzdat = novnst(iset)
      CALL define_nstgrd_locs(iset,nhdat,nzdat,2,hnstglb(l1:l2,:),&
                            & znstglb(l1:l2,1:nzdat),'V  ')
   ENDIF
ENDDO iset_530

!
!5.4 Local number of data points and buffer index
!------------------------------------------------
!

iproc_540: DO iproc=1,nprocs
   lloc = 0
   lglb_541: DO lglb=1,numhnstglb
      i = hnstglb(lglb,1)%icoord; j = hnstglb(lglb,1)%jcoord
      IF (.NOT.local_proc(i,j,iproc=iproc).OR.(i.EQ.int_fill)) CYCLE lglb_541
      lloc = lloc + 1
      iset_5411: DO iset=1,nonestsets
         l1 = lbnstglb(iset); l2 = l1 + nohnstglbv(iset)
         IF (lglb.GE.l1.AND.lglb.LT.l2) THEN
            nohnstvprocs(iproc,iset) = nohnstvprocs(iproc,iset) + 1
            IF (idloc.EQ.idprocs(iproc)) THEN
               nohnstatv(iset) = nohnstatv(iset) + 1
               IF (nohnstatv(iset).EQ.1) lbhnstatv(iset) = lloc 
            ENDIF
         ENDIF
      ENDDO iset_5411
   ENDDO lglb_541
ENDDO iproc_540
   
!
!5.5 Index maps (V-points)
!-------------------------
!
!5.5.1 Allocate
!--------------
!

nodat = MAXVAL(nohnstvprocs)
ALLOCATE (indexnstv(nodat,nprocs,nonestsets),STAT=errstat)
CALL error_alloc('indexnstv',3,(/nodat,nprocs,nonestsets/),kndint)
IF (nodat.GT.0) indexnstv = 0

!
!5.5.2 Define
!------------
!

iproc_552: DO iproc=1,nprocs
iset_552: DO iset=1,nonestsets
   l1 = lbnstglb(iset); l2 = l1 + nohnstglbv(iset) - 1
   lsub = 0
   lglb_5521: DO lglb=l1,l2
      i = hnstglb(lglb,1)%icoord; j = hnstglb(lglb,1)%jcoord
      IF (local_proc(i,j,iproc=iproc).AND.(i.NE.int_fill)) THEN
         lsub = lsub + 1
         indexnstv(lsub,iproc,iset) = lglb - l1 + 1
      ENDIF
   ENDDO lglb_5521
ENDDO iset_552
ENDDO iproc_552

!
!5.6 Local coordinate structures
!-------------------------------
!
!5.6.1 Allocate/initialise
!-------------------------
!
!---allocate
numhnst = SUM(nohnstatv)
ALLOCATE (hnstvtov(numhnst),STAT=errstat)
CALL error_alloc_struc('hnstvtov',1,(/numhnst/),'HRelativeCoords')
ALLOCATE (hnstctov(numhnst),STAT=errstat)
CALL error_alloc_struc('hnstctov',1,(/numhnst/),'HRelativeCoords')
ALLOCATE (vnstvtov(2,2,numhnst,numvnst),STAT=errstat)
CALL error_alloc_struc('vnstvtov',4,(/2,2,numhnst,numvnst/),'VRelativeCoords')
nodat = MAXVAL(nohnstatv)
ALLOCATE (znstloc(nodat,numvnst),STAT=errstat)
CALL error_alloc('znstloc',2,(/nodat,numvnst/),kndrtype)

!---initialise
CALL hrelativecoords_init(hnstvtov,.TRUE.)
CALL hrelativecoords_init(hnstctov,.TRUE.)
CALL vrelativecoords_init(vnstvtov,.TRUE.)

!
!5.6.2 Define
!------------
!

iproc_562: DO iproc=1,nprocs

   IF (idloc.EQ.idprocs(iproc)) THEN

      iset_5621: DO iset=1,nonestsets
         nhdat = nohnstatv(iset); nzdat = novnst(iset)
         l1 = lbhnstatv(iset); l2 = l1 + nhdat - 1
         lsub = 0
         lloc_56211: DO lloc=l1,l2
            lsub = lsub + 1
            lglb = indexnstv(lsub,iproc,iset) + lbnstglb(iset) - 1
            IF (hnstglb(lglb,1)%icoord.NE.int_fill) THEN
               hnstvtov(lloc)%icoord = hnstglb(lglb,1)%icoord - nc1loc + 1
               hnstvtov(lloc)%jcoord = hnstglb(lglb,1)%jcoord - nr1loc + 1
               hnstvtov(lloc)%weights = hnstglb(lglb,1)%weights
               hnstctov(lloc)%icoord = hnstglb(lglb,2)%icoord - nc1loc + 1
               hnstctov(lloc)%jcoord = hnstglb(lglb,2)%jcoord - nr1loc + 1
               hnstctov(lloc)%weights = hnstglb(lglb,2)%weights
            ENDIF
            IF (nzdat.GT.0) znstloc(lsub,1:nzdat) = znstglb(lglb,1:nzdat)
         ENDDO lloc_56211

         IF (nhdat*nzdat.GT.0) THEN
            CALL data_to_model_vcoords_3d(znstloc(1:nhdat,1:nzdat),&
                               & hnstvtov(l1:l2),vnstvtov(:,:,l1:l2,1:nzdat),&
                               & nhdat,1,nzdat,'V  ')
         ENDIF

      ENDDO iset_5621

   ENDIF

ENDDO iproc_562

!
!5.7 Index map (U,V-points)
!--------------------------
!

nohnstuvprocs = nohnstuprocs + nohnstvprocs
nodat = MAXVAL(nohnstuvprocs)
ALLOCATE (indexnstuv(nodat,nprocs,nonestsets),STAT=errstat)
CALL error_alloc('indexnstuv',3,(/nodat,nprocs,nonestsets/),kndint)
IF (nodat.GT.0) indexnstuv = 0

iset_570: DO iset=1,nonestsets
iproc_570: DO iproc=1,nprocs
   lloc_571: DO lloc=1,nohnstuprocs(iproc,iset)
      indexnstuv(lloc,iproc,iset) = indexnstu(lloc,iproc,iset)
   ENDDO lloc_571
   l1 = nohnstuprocs(iproc,iset)
   lloc_572: DO lloc=1,nohnstvprocs(iproc,iset)
      indexnstuv(l1+lloc,iproc,iset) = nohnstglbu(iset) + &
                                     & indexnstv(lloc,iproc,iset)
   ENDDO lloc_572
ENDDO iproc_570
ENDDO iset_570

!
!5.8 Deallocate
!--------------
!

DEALLOCATE (hnstglb,indexnstu,indexnstv,nohnstuprocs,nohnstvprocs,znstglb,&
          & znstloc)

!
!6. X-node points
!----------------
!
!6.1 Global buffer index
!-----------------------
!

lbuf = 1
iset_1410: DO iset=1,nonestsets
   lbnstglb(iset) = lbuf
   lbuf = lbuf + nohnstglbx(iset)
ENDDO iset_1410

!
!6.2 Allocate work space arrays
!------------------------------
!

numhnstglb = SUM(nohnstglbx)
ALLOCATE (hnstglb(numhnstglb,1),STAT=errstat)
CALL error_alloc_struc('hnstglb',2,(/numhnstglb,1/),'HRelativeCoords')
CALL hrelativecoords_init(hnstglb,.TRUE.)
ALLOCATE (znstglb(numhnstglb,numvnst),STAT=errstat)
CALL error_alloc('znstglb',2,(/numhnstglb,numvnst/),kndrtype)
ALLOCATE (nohnstxprocs(nprocs,nonestsets),STAT=errstat)
CALL error_alloc('nohnstxprocs',2,(/nprocs,nonestsets/),kndint)
nohnstxprocs = 0

!
!6.3 Define coordinate structure
!-------------------------------
!

iset_1430: DO iset=1,nonestsets
   nhdat = nohnstglbx(iset)
   IF (nhdat.GT.0) THEN
      l1 = lbnstglb(iset); l2 = l1 + nhdat - 1
      nzdat = novnst(iset)
      CALL define_nstgrd_locs(iset,nhdat,nzdat,1,hnstglb(l1:l2,:),&
                            & znstglb(l1:l2,1:nzdat),'X  ')
   ENDIF
ENDDO iset_1430

!
!6.4 Local number of data points and buffer index
!------------------------------------------------
!

iproc_1440: DO iproc=1,nprocs
   lloc = 0
   lglb_1441: DO lglb=1,numhnstglb
      i = hnstglb(lglb,1)%icoord; j = hnstglb(lglb,1)%jcoord
      IF (.NOT.local_proc(i,j,iproc=iproc).OR.(i.EQ.int_fill)) CYCLE lglb_1441
      lloc = lloc + 1
      iset_14411: DO iset=1,nonestsets
         l1 = lbnstglb(iset); l2 = l1 + nohnstglbx(iset)
         IF (lglb.GE.l1.AND.lglb.LT.l2) THEN
            nohnstxprocs(iproc,iset) = nohnstxprocs(iproc,iset) + 1
            IF (idloc.EQ.idprocs(iproc)) THEN
               nohnstatx(iset) = nohnstatx(iset) + 1
               IF (nohnstatx(iset).EQ.1) lbhnstatx(iset) = lloc 
            ENDIF
         ENDIF
      ENDDO iset_14411
   ENDDO lglb_1441
ENDDO iproc_1440
 
!
!6.5 Index map (X-points)
!------------------------
!
!6.5.1 Allocate
!--------------
!

nodat = MAXVAL(nohnstxprocs)
ALLOCATE (indexnstx(nodat,nprocs,nonestsets),STAT=errstat)
CALL error_alloc('indexnstx',3,(/nodat,nprocs,nonestsets/),kndint)
IF (nodat.GT.0) indexnstx = 0

!
!6.5.2 Define
!-----------
!

iproc_1452: DO iproc=1,nprocs
iset_1452: DO iset=1,nonestsets
   l1 = lbnstglb(iset); l2 = l1 + nohnstglbx(iset) - 1
   lsub = 0
   lglb_14521: DO lglb=l1,l2
      i = hnstglb(lglb,1)%icoord; j = hnstglb(lglb,1)%jcoord
      IF (local_proc(i,j,iproc=iproc).AND.(i.NE.int_fill)) THEN
         lsub = lsub + 1
         indexnstx(lsub,iproc,iset) = lglb - l1 + 1
      ENDIF
   ENDDO lglb_14521
ENDDO iset_1452
ENDDO iproc_1452

!
!6.6 Local coordinate structures
!-------------------------------
!
!6.6.1 Allocate/initialise
!-------------------------
!
!---allocate
numhnst = SUM(nohnstatx)
ALLOCATE (hnstvtox(numhnst),STAT=errstat)
CALL error_alloc_struc('hnstvtox',1,(/numhnst/),'HRelativeCoords')
ALLOCATE (vnstvtox(2,2,numhnst,numvnst),STAT=errstat)
CALL error_alloc_struc('vnstvtox',4,(/2,2,numhnst,numvnst/),'VRelativeCoords')
nodat = MAXVAL(nohnstatx)
ALLOCATE (znstloc(nodat,numvnst),STAT=errstat)
CALL error_alloc('znstloc',2,(/nodat,numvnst/),kndrtype)

!---initialise
CALL hrelativecoords_init(hnstvtox,.TRUE.)
CALL vrelativecoords_init(vnstvtox,.TRUE.)

!
!6.6.2 Define
!------------
!

iproc_1462: DO iproc=1,nprocs

   IF (idloc.EQ.idprocs(iproc)) THEN

      iset_14621: DO iset=1,nonestsets
         nhdat = nohnstatx(iset); nzdat = novnst(iset)
         l1 = lbhnstatx(iset); l2 = l1 + nhdat - 1
         lsub = 0
         lloc_146211: DO lloc=l1,l2
            lsub = lsub + 1
            lglb = indexnstx(lsub,iproc,iset) + lbnstglb(iset) - 1
            IF (hnstglb(lglb,1)%icoord.NE.int_fill) THEN
               hnstvtox(lloc)%icoord = hnstglb(lglb,1)%icoord - nc1loc + 1
               hnstvtox(lloc)%jcoord = hnstglb(lglb,1)%jcoord - nr1loc + 1
               hnstvtox(lloc)%weights = hnstglb(lglb,1)%weights
            ENDIF
            IF (nzdat.GT.0) znstloc(lsub,1:nzdat) = znstglb(lglb,1:nzdat)
         ENDDO lloc_146211

         IF (nhdat*nzdat.GT.0) THEN
            CALL data_to_model_vcoords_3d(znstloc(1:nhdat,1:nzdat),&
                                & hnstvtox(l1:l2),vnstvtox(:,:,l1:l2,1:nzdat),&
                                & nhdat,1,nzdat,'X  ')
         ENDIF

      ENDDO iset_14621

   ENDIF

ENDDO iproc_1462

!
!6.7 Deallocate
!--------------
!

DEALLOCATE (hnstglb,znstglb,znstloc)

!
!7. Y-node points
!----------------
!
!7.1 Global buffer index
!-----------------------
!

lbuf = 1
iset_1510: DO iset=1,nonestsets
   lbnstglb(iset) = lbuf
   lbuf = lbuf + nohnstglby(iset)
ENDDO iset_1510

!
!7.2 Allocate work space arrays
!------------------------------
!

numhnstglb = SUM(nohnstglby)
ALLOCATE (hnstglb(numhnstglb,1),STAT=errstat)
CALL error_alloc_struc('hnstglb',2,(/numhnstglb,1/),'HRelativeCoords')
CALL hrelativecoords_init(hnstglb,.TRUE.)
ALLOCATE (znstglb(numhnstglb,numvnst),STAT=errstat)
CALL error_alloc('znstglb',2,(/numhnstglb,numvnst/),kndrtype)
ALLOCATE (nohnstyprocs(nprocs,nonestsets),STAT=errstat)
CALL error_alloc('nohnstyprocs',2,(/nprocs,nonestsets/),kndint)
nohnstyprocs = 0

!
!7.3 Define coordinate structure
!-------------------------------
!

iset_1530: DO iset=1,nonestsets
   nhdat = nohnstglby(iset)
   IF (nhdat.GT.0) THEN
      l1 = lbnstglb(iset); l2 = l1 + nhdat - 1
      nzdat = novnst(iset)
      CALL define_nstgrd_locs(iset,nhdat,nzdat,1,hnstglb(l1:l2,:),&
                            & znstglb(l1:l2,1:nzdat),'Y  ')
   ENDIF
ENDDO iset_1530

!
!7.4 Local number of data points and buffer index
!------------------------------------------------
!

iproc_1540: DO iproc=1,nprocs
   lloc = 0
   lglb_1541: DO lglb=1,numhnstglb
      i = hnstglb(lglb,1)%icoord; j = hnstglb(lglb,1)%jcoord
      IF (.NOT.local_proc(i,j,iproc=iproc).OR.(i.EQ.int_fill)) CYCLE lglb_1541
      lloc = lloc + 1
      iset_15411: DO iset=1,nonestsets
         l1 = lbnstglb(iset); l2 = l1 + nohnstglby(iset)
         IF (lglb.GE.l1.AND.lglb.LT.l2) THEN
            nohnstyprocs(iproc,iset) = nohnstyprocs(iproc,iset) + 1
            IF (idloc.EQ.idprocs(iproc)) THEN
               nohnstaty(iset) = nohnstaty(iset) + 1
               IF (nohnstaty(iset).EQ.1) lbhnstaty(iset) = lloc 
            ENDIF
         ENDIF
      ENDDO iset_15411
   ENDDO lglb_1541
ENDDO iproc_1540
   
!
!7.5 Index maps (Y-points)
!-------------------------
!
!7.5.1 Allocate
!--------------
!

nodat = MAXVAL(nohnstyprocs)
ALLOCATE (indexnsty(nodat,nprocs,nonestsets),STAT=errstat)
CALL error_alloc('indexnsty',3,(/nodat,nprocs,nonestsets/),kndint)
IF (nodat.GT.0) indexnsty = 0

!
!7.5.2 Define
!------------
!

iproc_1552: DO iproc=1,nprocs
iset_1552: DO iset=1,nonestsets
   l1 = lbnstglb(iset); l2 = l1 + nohnstglby(iset) - 1
   lsub = 0
   lglb_15521: DO lglb=l1,l2
      i = hnstglb(lglb,1)%icoord; j = hnstglb(lglb,1)%jcoord
      IF (local_proc(i,j,iproc=iproc).AND.(i.NE.int_fill)) THEN
         lsub = lsub + 1
         indexnsty(lsub,iproc,iset) = lglb - l1 + 1
      ENDIF
   ENDDO lglb_15521
ENDDO iset_1552
ENDDO iproc_1552

!
!7.6 Local coordinate structures
!-------------------------------
!
!7.6.1 Allocate/initialise
!-------------------------
!
!---allocate
numhnst = SUM(nohnstaty)
ALLOCATE (hnstutoy(numhnst),STAT=errstat)
CALL error_alloc_struc('hnstutoy',1,(/numhnst/),'HRelativeCoords')
ALLOCATE (vnstutoy(2,2,numhnst,numvnst),STAT=errstat)
CALL error_alloc_struc('vnstutoy',4,(/2,2,numhnst,numvnst/),'VRelativeCoords')
nodat = MAXVAL(nohnstaty)
ALLOCATE (znstloc(nodat,numvnst),STAT=errstat)
CALL error_alloc('znstloc',2,(/nodat,numvnst/),kndrtype)

!---initialise
CALL hrelativecoords_init(hnstutoy,.TRUE.)
CALL vrelativecoords_init(vnstutoy,.TRUE.)

!
!7.6.2 Define
!------------
!

iproc_1562: DO iproc=1,nprocs

   IF (idloc.EQ.idprocs(iproc)) THEN

      iset_15621: DO iset=1,nonestsets
         nhdat = nohnstaty(iset); nzdat = novnst(iset)
         l1 = lbhnstaty(iset); l2 = l1 + nhdat - 1
         lsub = 0
         lloc_156211: DO lloc=l1,l2
            lsub = lsub + 1
            lglb = indexnsty(lsub,iproc,iset) + lbnstglb(iset) - 1
            IF (hnstglb(lglb,1)%icoord.NE.int_fill) THEN
               hnstutoy(lloc)%icoord = hnstglb(lglb,1)%icoord - nc1loc + 1
               hnstutoy(lloc)%jcoord = hnstglb(lglb,1)%jcoord - nr1loc + 1
               hnstutoy(lloc)%weights = hnstglb(lglb,1)%weights
            ENDIF
            IF (nzdat.GT.0) znstloc(lsub,1:nzdat) = znstglb(lglb,1:nzdat)
         ENDDO lloc_156211

         IF (nhdat*nzdat.GT.0) THEN
            CALL data_to_model_vcoords_3d(znstloc(1:nhdat,1:nzdat),&
                               & hnstutoy(l1:l2),vnstutoy(:,:,l1:l2,1:nzdat),&
                               & nhdat,1,nzdat,'Y  ')
         ENDIF

      ENDDO iset_15621

   ENDIF

ENDDO iproc_1562

!
!7.7 Index map (X,Y-points)
!--------------------------
!

nohnstxyprocs = nohnstxprocs + nohnstyprocs
nodat = MAXVAL(nohnstxyprocs)
ALLOCATE (indexnstxy(nodat,nprocs,nonestsets),STAT=errstat)
CALL error_alloc('indexnstxy',3,(/nodat,nprocs,nonestsets/),kndint)
IF (nodat.GT.0) indexnstxy = 0

iset_1570: DO iset=1,nonestsets
iproc_1570: DO iproc=1,nprocs
   lloc_1571: DO lloc=1,nohnstxprocs(iproc,iset)
      indexnstxy(lloc,iproc,iset) = indexnstx(lloc,iproc,iset)
   ENDDO lloc_1571
   l1 = nohnstxprocs(iproc,iset)
   lloc_1572: DO lloc=1,nohnstyprocs(iproc,iset)
      indexnstxy(l1+lloc,iproc,iset) = nohnstglbx(iset) + &
                                     & indexnsty(lloc,iproc,iset)
   ENDDO lloc_1572
ENDDO iproc_1570
ENDDO iset_1570

!
!7.8 Deallocate
!--------------
!

DEALLOCATE (hnstglb,indexnstx,indexnsty,nohnstxprocs,nohnstyprocs,znstglb,&
          & znstloc)


!
!8. Deallocate work space array
!------------------------------
!

DEALLOCATE (lbnstglb)

CALL log_timer_out(npcc,itm_nest)


RETURN

END SUBROUTINE nest_locations

!========================================================================

SUBROUTINE read_nstgrd_abs(iset,nhdat,nzdat,xcoord,ycoord,zcoord)
!************************************************************************
!
! *read_nstgrd_abs* Read absolute geographical positions of nested open
!                   boundaries in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - define_nstgrd_locs
!
! External calls -
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iset, nhdat, nzdat
REAL, INTENT(OUT), DIMENSION(nhdat) :: xcoord, ycoord
REAL, INTENT(OUT), DIMENSION(nhdat,nzdat) :: zcoord

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iset*      INTEGER Output file number
!*nhdat*     INTEGER Number of horizontal data locations
!*nzdat*     INTEGER Number of vertical data points
!*xcoord*    REAL    X-coordinates of nest locations
!*ycoord*    REAL    Y-coordinates of nest locations
!*zcoord*    REAL    Z-coordinates of nest locations
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, SAVE :: numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_nstgrd_abs'
CALL log_timer_in()

filepars = modfiles(io_nstabs,iset,1)

!
!1. File header
!--------------
!

IF (filepars%iostat.EQ.0) THEN 
   CALL open_filepars(filepars)
   CALL read_glbatts_mod(filepars)
   numvars = filepars%novars
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
   CALL read_varatts_mod(filepars,varatts,numvars)
   varid = 0
ENDIF

!
!2. Read locations
!-----------------
!

varid = varid + 1
CALL read_vars(xcoord,filepars,varid,varatts=varatts)
varid = varid + 1
CALL read_vars(ycoord,filepars,varid,varatts=varatts)
IF (nzdat.GT.0) THEN
   varid = varid + 1
   CALL read_vars(zcoord,filepars,varid,varatts=varatts)
ENDIF

!
!3. Finalise
!-----------
!

IF (varid.EQ.numvars) THEN
   CALL close_filepars(filepars)
   DEALLOCATE(varatts)
ENDIF
modfiles(io_nstabs,iset,1) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE read_nstgrd_abs

!========================================================================

SUBROUTINE read_nstgrd_rel(iset,nhdat,nzdat,nonodes,hnests,zcoord)
!************************************************************************
!
! *read_nstgrd_rel* Read relative geographical positions of nested open
!                   boundaries in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - define_nstgrd_locs
!
! External calls -
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE nestgrids
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iset, nhdat, nonodes, nzdat
REAL, INTENT(OUT), DIMENSION(nhdat,nzdat) :: zcoord
TYPE (HRelativeCoords), INTENT(OUT), DIMENSION(nhdat,nonodes) :: hnests

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iset*      INTEGER Output file number
!*nhdat*     INTEGER Number of horizontal data locations
!*nzdat*     INTEGER Number of vertical data points
!*nonodes*   INTEGER Number of nodes for interpolation
!*hnests*    DERIVED Relative coordinates of nest locations
!*zcoord*    REAL    Z-coordinates of nest locations
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: inode, l
INTEGER, SAVE :: numvars, varid
TYPE (FileParams) :: filepars
REAL, DIMENSION(2,2,nhdat) :: weightsdat
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_nstgrd_rel'
CALL log_timer_in()

filepars = modfiles(io_nstrel,iset,1)

!
!1. File header
!--------------
!

IF (filepars%iostat.EQ.0) THEN

!  ---open file
   CALL open_filepars(filepars)

!  ---read global attributes
   CALL read_glbatts_mod(filepars)

!  ---read variable attributes
   numvars = filepars%novars
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
   CALL read_varatts_mod(filepars,varatts,numvars)

!  ---initialise
   varid = 0

ENDIF

!
!2. Read locations
!-----------------
!

inode_210: DO inode=1,nonodes
   
   varid = varid + 1
   CALL read_vars(hnests(:,inode)%icoord,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(hnests(:,inode)%jcoord,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(weightsdat,filepars,varid,varatts=varatts)
   l_211: DO l=1,nhdat
      hnests(l,inode)%weights = weightsdat(:,:,l)
   ENDDO l_211
   IF (inode.EQ.nonodes.AND.novnst(iset).GT.0) THEN
      varid = varid + 1
      CALL read_vars(zcoord,filepars,varid,varatts=varatts)
   ENDIF

ENDDO inode_210
   
!
!3. Finalise
!-----------
!

IF (varid.EQ.numvars) THEN
   CALL close_filepars(filepars)
   DEALLOCATE(varatts)
ENDIF

modfiles(io_nstrel,iset,1) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE read_nstgrd_rel

!========================================================================

SUBROUTINE read_nstgrd_spec
!************************************************************************
!
! *read_nstgrd_spec* Read number of output locations and specifier arrays
!                    for nesting in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.5
!
! Description - 
!
! Reference -
!
! Calling program - define_nstgrd_spec
!
! External calls -
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE nestgrids
USE switches
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_nstgrd_spec'
CALL log_timer_in()

filepars = modfiles(io_nstspc,1,1)

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
!2. Read data
!------------
!
!---type of grids
CALL read_vars(nestcoords,filepars,1,varatts=varatts)

!---number of horizontal nest points
CALL read_vars(nohnstglbc,filepars,2,varatts=varatts)
CALL read_vars(nohnstglbu,filepars,3,varatts=varatts)
CALL read_vars(nohnstglbv,filepars,4,varatts=varatts)

!---number of vertical nest points
CALL read_vars(novnst,filepars,5,varatts=varatts)

!---type of 2-D output
CALL read_vars(inst2dtype,filepars,6,varatts=varatts)

!---tangential conditions
CALL read_vars(nohnstglbx,filepars,7,varatts=varatts)
CALL read_vars(nohnstglby,filepars,8,varatts=varatts)

!---sediment fractions used in nesting
varid = 8
IF (iopt_sed.GT.0) THEN
   varid = varid + 1
   CALL read_vars(nosednst,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(instsed,filepars,varid,varatts=varatts)
ENDIF

!---biological variables used in nesting
IF (iopt_biolgy.GT.0) THEN
   varid = varid + 1
   CALL read_vars(nobionst,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(instbio,filepars,varid,varatts=varatts)
ENDIF

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_nstspc,1,1) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE read_nstgrd_spec

!========================================================================

SUBROUTINE update_nest_data_2d(iddesc)
!************************************************************************
!
! *update_nest_data_2d* Interpolate and write 2-D data at nested open
!                       boundaries
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.7
!
! Description - 
!
! Reference -
!
! Calling program - current_2d
!
! External calls - write_nest_data_2d
!
! Module calls - error_alloc, intpol2d_model_to_data_2d, loop_index
!
!************************************************************************
!
USE currents
USE depths
USE gridpars
USE iopars
USE nestgrids
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE grid_interp, ONLY: intpol2d_model_to_data_2d
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: loop_index

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc

!
! Name      Type      Purpose
!-------------------------------------------------------------------------
!*iddesc*   INTEGER   Data file id
!
!-------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: iset, l1, l2, maxstatsuv, maxstatsxy, nhdatlocuv, nhdatlocxy, &
         & nodatu, nodatv, nodatx, nodaty, nostatsglbuv, nostatsglbxy, npcc, &
         & varid
INTEGER, DIMENSION(4) :: nhdims
REAL, ALLOCATABLE, DIMENSION(:) :: outvals2d


IF (.NOT.ANY(modfiles(iddesc,1:nonestsets,2)%defined)) RETURN

procname(pglev+1) = 'update_nest_data_2d'
CALL log_timer_in(npcc)

iset_1000: DO iset=1,nonestsets

   IF (modfiles(iddesc,iset,2)%defined.AND.&
     & loop_index(modfiles(iddesc,iset,2)%tlims,nt)) THEN

      SELECT CASE (iddesc)

!
!1. Normal components
!--------------------
!

         CASE (io_2uvnst)

!
!1.1. Allocate
!-------------
!

            maxstatsuv = MAXVAL(nohnstuvprocs)
            nostatsglbuv = SUM(nohnstuvprocs(:,iset))
            nodatu = nohnstatu(iset); nodatv = nohnstatv(iset)
            nhdatlocuv = nodatu + nodatv
            ALLOCATE (outvals2d(nhdatlocuv),STAT=errstat)
            CALL error_alloc('outvals2d',1,(/nhdatlocuv/),kndrtype)

!
!1.2. Transports
!---------------
!

            IF (inst2dtype(iset).NE.2) THEN
               nhdims = nhalo
               l1 = lbhnstatu(iset); l2 = l1 + nodatu - 1
               IF (nodatu.GT.0) THEN
                  CALL intpol2d_model_to_data_2d(udvel,outvals2d(1:nodatu),&
                                        & hnstutou(l1:l2),nhdims,nodatu,1)
               ENDIF
               l1 = lbhnstatv(iset); l2 = l1 + nodatv - 1
               IF (nodatv.GT.0) THEN
                  CALL intpol2d_model_to_data_2d(vdvel,&
                                        & outvals2d(nodatu+1:nhdatlocuv),&
                                        & hnstvtov(l1:l2),nhdims,nodatv,1)
               ENDIF
               CALL write_nest_data_2d(iset,1,outvals2d,maxstatsuv,nhdatlocuv,&
                                     & nostatsglbuv,nohnstuvprocs(:,iset),&
                                     & indexnstuv(:,:,iset),iddesc)
            ENDIF

!
!1.3. Surface elevation
!----------------------
!

            IF (inst2dtype(iset).NE.3) THEN
               varid = MERGE(1,2,inst2dtype(iset).EQ.2)
               nhdims = 1
               l1 = lbhnstatu(iset); l2 = l1 + nodatu - 1
               IF (nodatu.GT.0) THEN
                  CALL intpol2d_model_to_data_2d(zeta,outvals2d(1:nodatu),&
                                        & hnstctou(l1:l2),nhdims,nodatu,1)
               ENDIF
               l1 = lbhnstatv(iset); l2 = l1 + nodatv - 1
               IF (nodatv.GT.0) THEN
                  CALL intpol2d_model_to_data_2d(zeta,&
                                        & outvals2d(nodatu+1:nhdatlocuv),&
                                        & hnstctov(l1:l2),nhdims,nodatv,1)
               ENDIF
               CALL write_nest_data_2d(iset,varid,outvals2d,maxstatsuv,&
                                     & nhdatlocuv,nostatsglbuv,&
                                     & nohnstuvprocs(:,iset),&
                                     & indexnstuv(:,:,iset),iddesc)
            ENDIF

!
!1.4. Deallocate
!---------------
!

            DEALLOCATE (outvals2d)

!
!2. Tangential components
!------------------------
!

         CASE (io_2xynst)

!
!2.1. Allocate
!-------------
!

            maxstatsxy = MAXVAL(nohnstxyprocs)
            nostatsglbxy = SUM(nohnstxyprocs(:,iset))
            nodatx = nohnstatx(iset); nodaty = nohnstaty(iset)
            nhdatlocxy = nodatx + nodaty
            ALLOCATE (outvals2d(nhdatlocxy),STAT=errstat)
            CALL error_alloc('outvals2d',1,(/nhdatlocxy/),kndrtype)

!
!2.2. Transports
!---------------
!

            nhdims = nhalo
            l1 = lbhnstatx(iset); l2 = l1 + nodatx - 1
            IF (nodatx.GT.0) THEN
               CALL intpol2d_model_to_data_2d(vdvel,outvals2d(1:nodatx),&
                                        & hnstvtox(l1:l2),nhdims,nodatx,1)
            ENDIF
            l1 = lbhnstaty(iset); l2 = l1 + nodaty - 1
            IF (nodaty.GT.0) THEN
               CALL intpol2d_model_to_data_2d(udvel,&
                                  & outvals2d(nodatx+1:nhdatlocxy),&
                                  & hnstutoy(l1:l2),nhdims,nodaty,1)
            ENDIF
            CALL write_nest_data_2d(iset,1,outvals2d,maxstatsxy,nhdatlocxy,&
                                  & nostatsglbxy,nohnstxyprocs(:,iset),&
                                  & indexnstxy(:,:,iset),iddesc)

!
!2.3. Deallocate
!---------------
!

            DEALLOCATE (outvals2d)

      END SELECT

   ENDIF

ENDDO iset_1000

CALL log_timer_out(npcc,itm_nest)


RETURN

END SUBROUTINE update_nest_data_2d

!========================================================================

SUBROUTINE update_nest_data_3d(iddesc)
!************************************************************************
!
! *update_nest_data_3d* Interpolate and write 3-D current data at nested open
!                       boundaries
!
! Version - @(COHERENS)fNested_Grids.f90  V2.7
!
! Author - Patrick Luyten
!
! Description - 
!
! Reference -
!
! Calling program - current_corr
!
! External calls - write_nest_data_3d
!
! Module calls - error_alloc, intpol3d_model_to_data_3d, loop_index
!
!************************************************************************
!
USE currents
USE gridpars
USE iopars
USE nestgrids
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE grid_interp, ONLY: intpol3d_model_to_data_3d
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: loop_index

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc

!
! Name      Type      Purpose
!-------------------------------------------------------------------------
!*iddesc*   INTEGER   Data file id
!
!-------------------------------------------------------------------------
!
!
!*Local variables
!
INTEGER :: iset, k, l1, l2, maxstats, nhdatloc, nodatu, nodatv, nodatx, &
         & nodaty, nostatsglb, npcc, nzdat
INTEGER, DIMENSION(4) :: nhdims
REAL, ALLOCATABLE, DIMENSION(:,:) :: outvals3d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dc


IF (.NOT.ANY(modfiles(iddesc,1:nonestsets,2)%defined)) RETURN

procname(pglev+1) = 'update_nest_data_3d'
CALL log_timer_in(npcc)

!
!1. Allocate array
!-----------------
!

IF (ANY(modfiles(io_3uvnst,1:nonestsets,2)%defined)) THEN
   ALLOCATE (array3dc(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
   CALL error_alloc('array3dc',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
ENDIF

!
!2. Interpolate and write
!------------------------
!

iset_200: DO iset=1,nonestsets

   IF (modfiles(iddesc,iset,2)%defined.AND.&
     & loop_index(modfiles(iddesc,iset,2)%tlims,nt)) THEN

      SELECT CASE (iddesc)

!
!2.1 Normal components
!---------------------
!

         CASE (io_3uvnst)

!
!2.1.1 Allocate
!--------------
!

            nhdims = 1
            maxstats = MAXVAL(nohnstuvprocs); nzdat = novnst(iset)
            nostatsglb = SUM(nohnstuvprocs(:,iset))
            nodatu = nohnstatu(iset); nodatv = nohnstatv(iset)
            nhdatloc = nodatu + nodatv
            ALLOCATE (outvals3d(nhdatloc,nzdat),STAT=errstat)
            CALL error_alloc('outvals3d',2,(/nhdatloc,nzdat/),kndrtype)

!
!2.1.2 Interpolate at U-nodes
!----------------------------
!

            k_220: DO k=1,nz
               array3dc(:,:,k) = uvel(0:ncloc+1,0:nrloc+1,k) - &
                               & umvel(0:ncloc+1,0:nrloc+1)
            ENDDO k_220
            l1 = lbhnstatu(iset); l2 = l1 + nodatu - 1
            IF (nodatu.GT.0) THEN
               CALL intpol3d_model_to_data_3d(array3dc,outvals3d(1:nodatu,:),&
                              & hnstutou(l1:l2),vnstutou(:,:,l1:l2,1:nzdat),&
                              & nhdims,nz,nodatu,nzdat,1)
            ENDIF

!
!2.1.3 Interpolate at V-nodes
!----------------------------
!

            k_230: DO k=1,nz
               array3dc(:,:,k) = vvel(0:ncloc+1,0:nrloc+1,k) - &
                               & vmvel(0:ncloc+1,0:nrloc+1)
            ENDDO k_230
            l1 = lbhnstatv(iset); l2 = l1 + nodatv - 1
            IF (nodatv.GT.0) THEN
               CALL intpol3d_model_to_data_3d(array3dc,&
                              & outvals3d(nodatu+1:nhdatloc,:),hnstvtov(l1:l2),&
                              & vnstvtov(:,:,l1:l2,1:nzdat),nhdims,nz,nodatv,&
                              & nzdat,1)
            ENDIF

!
!2.1.4 Combine and write
!-----------------------
!

            CALL write_nest_data_3d(iset,outvals3d,maxstats,nzdat,nhdatloc,&
                                  & nostatsglb,nohnstuvprocs(:,iset),&
                                  & indexnstuv(:,:,iset),iddesc)

!           ---deallocate
            DEALLOCATE (outvals3d)


!
!2.2 Tangential components
!-------------------------
!

         CASE (io_3xynst)

!
!2.2.1 Interpolate and write
!---------------------------
!
!2.2.1 Allocate
!--------------
!

            nhdims = nhalo
            maxstats = MAXVAL(nohnstxyprocs); nzdat = novnst(iset)
            nostatsglb = SUM(nohnstxyprocs(:,iset))
            nodatx = nohnstatx(iset); nodaty = nohnstaty(iset)
            nhdatloc = nodatx + nodaty
            ALLOCATE (outvals3d(nhdatloc,nzdat),STAT=errstat)
            CALL error_alloc('outvals3d',2,(/nhdatloc,nzdat/),kndrtype)

!
!2.2.2 Interpolate at X-nodes
!----------------------------
!

            l1 = lbhnstatx(iset); l2 = l1 + nodatx - 1
            IF (nodatx.GT.0) THEN
               CALL intpol3d_model_to_data_3d(vvel,outvals3d(1:nodatx,:),&
                              & hnstvtox(l1:l2),vnstvtox(:,:,l1:l2,1:nzdat),&
                              & nhdims,nz,nodatx,nzdat,1)
            ENDIF

!
!2.2.3 Interpolate at Y-nodes
!----------------------------
!

            l1 = lbhnstaty(iset); l2 = l1 + nodaty - 1
            IF (nodaty.GT.0) THEN
               CALL intpol3d_model_to_data_3d(uvel,&
                              & outvals3d(nodatx+1:nhdatloc,:),hnstutoy(l1:l2),&
                              & vnstutoy(:,:,l1:l2,1:nzdat),nhdims,nz,nodaty,&
                              & nzdat,1)
            ENDIF


!
!2.2.4 Combine and write
!-----------------------
!

            CALL write_nest_data_3d(iset,outvals3d,maxstats,nzdat,nhdatloc,&
                                  & nostatsglb,nohnstxyprocs(:,iset),&
                                  & indexnstxy(:,:,iset),iddesc)

!           ---deallocate
            DEALLOCATE (outvals3d)

      END SELECT

   ENDIF

ENDDO iset_200

!
!3. Deallocate array
!-------------------
!

IF (ALLOCATED(array3dc)) DEALLOCATE (array3dc)

CALL log_timer_out(npcc,itm_nest)


RETURN

END SUBROUTINE update_nest_data_3d

!========================================================================

SUBROUTINE update_nest_data_prof(psic,novars,novarsnst,instvars,iddesc)
!************************************************************************
!
! *update_nest_data_prof* Interpolate and write scalar profile data at nested
!                         open boundaries
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.1.2
!
! Description - 
!
! Reference -
!
! Calling program - salinity_equation, temperature_equation
!
! External calls - write_nest_data_prof
!
! Module calls - error_alloc, intpol3d_model_to_data_3d, loop_index
!
!************************************************************************
!
USE gridpars
USE iopars
USE nestgrids
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE grid_interp, ONLY: intpol3d_model_to_data_3d
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: loop_index

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: iddesc, novars
INTEGER, INTENT(IN), DIMENSION(nonestsets) :: novarsnst 
INTEGER, INTENT(IN), DIMENSION(novars,nonestsets) :: instvars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,&
                          & 1-nhalo:nrloc+nhalo,nz,novars) :: psic

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*psic*     REAL    C-node quantity for nesting
!*novars*   INTEGER Total number of 3-D ("state") variables
!*instvars* INTEGER Indices of the variables to be nested
!*iddesc*   INTEGER File id of nest file

!
!------------------------------------------------------------------------------
!*Local variables
!
INTEGER :: iivar, iset, ivar, l1, l2, maxstats, nhdatloc, nostatsglb, npcc, &
         & nzdat
INTEGER, DIMENSION(4) :: nhdims
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: outvals


IF (.NOT.ANY(modfiles(iddesc,1:nonestsets,2)%defined)) RETURN

procname(pglev+1) = 'update_nest_data_prof'
CALL log_timer_in(npcc)

iset_1000: DO iset=1,nonestsets

   IF (modfiles(iddesc,iset,2)%defined.AND.&
     & loop_index(modfiles(iddesc,iset,2)%tlims,nt)) THEN

!
!1. Allocate
!-----------
!


      nhdims = nhalo
      maxstats = MAXVAL(nohnstcprocs); nzdat = novnst(iset)
      nostatsglb = SUM(nohnstcprocs(:,iset))
      nhdatloc = nohnstatc(iset)
      ALLOCATE (outvals(nhdatloc,nzdat,novarsnst(iset)),STAT=errstat)
      CALL error_alloc('outvals',3,(/nhdatloc,nzdat,novarsnst(iset)/),kndrtype)
!
!2. Interpolate
!--------------
!

      l1 = lbhnstatc(iset); l2 = l1 + nhdatloc - 1
      IF (nhdatloc.GT.0) THEN
         ivar_210: DO ivar=1,novarsnst(iset)
            iivar = instvars(ivar,iset)
            CALL intpol3d_model_to_data_3d(psic(:,:,:,iivar),&
                                         & outvals(:,:,ivar),hnstctoc(l1:l2),&
                                         & vnstctoc(:,:,l1:l2,1:nzdat),nhdims,&
                                         & nz,nhdatloc,nzdat,1)
         ENDDO ivar_210
      ENDIF

!
!3. Combine and write
!--------------------
!

      CALL write_nest_data_prof(iddesc,iset,outvals,maxstats,nzdat,&
                              & novarsnst(iset),nhdatloc,nostatsglb,&
                              & nohnstcprocs(:,iset),indexnstc(:,:,iset))

!
!4. Deallocate
!-------------
!

      DEALLOCATE (outvals)

   ENDIF

ENDDO iset_1000

CALL log_timer_out(npcc,itm_nest)


RETURN

END SUBROUTINE update_nest_data_prof

!========================================================================

SUBROUTINE write_nest_data_2d(iset,varid,outvals,maxstats,nhdatloc,nostatsglb,&
                            & nhdatprocs,indexprocs,iddesc)
!************************************************************************
!
! *write_nest_data_2d* Write 2-D interpolated nest data in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - update_nest_data_2d
!
! External calls -
!
! Module calls - combine_write_stats_loc, error_alloc_struc, lim_dims,
!                open_filepars, output_flag, set_modfiles_atts,
!                set_modvars_atts, write_atts_mod, write_time
!
!************************************************************************
!
USE datatypes
USE iopars
USE paralpars
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc_struc
USE inout_paral, ONLY: combine_write_stats_loc
USE inout_routines, ONLY: open_filepars, output_flag, write_atts_mod, write_time
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines
USE utility_routines, ONLY: lim_dims

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, iset, maxstats, nhdatloc, nostatsglb, varid
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nhdatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: indexprocs
REAL, INTENT(IN), DIMENSION(nhdatloc) :: outvals

!
! Name        Type     Purpose
!------------------------------------------------------------------------------
!*iset*       INTEGER  Set index
!*varid*      INTEGER  Variable number in the data file
!*outvals*    REAL     Output data
!*maxstats*   INTEGER  First dimension of array indexprocs
!*nhdatloc*   INTEGER  Local number of horizontal output locations
!*nostatsglb* INTEGER  Number of stations in global array
!*nhdatprocs* INTEGER  Number of output locations per process
!*indexprocs* INTEGER  Indices of local output array within global output array
!*iddesc*     INTEGER  Data file id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
INTEGER :: numvars, tlimit
TYPE (FileParams) :: filepars
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'write_nest_data_2d'
CALL log_timer_in()

!
!1. Disable writing if needed
!----------------------------
!

filepars = modfiles(iddesc,iset,2)
flag = output_flag(CDatetime,filepars%tskips)
IF (.NOT.flag) GOTO 1000

!
!2. Write file header on first call
!----------------------------------
!

IF (master.AND.filepars%iostat.EQ.0) THEN

!  ---file attributes
   CALL set_modfiles_atts(iddesc,iset,2)
   filepars = modfiles(iddesc,iset,2)
   numvars = filepars%novars + 1
   tlimit = lim_dims(filepars%tlims)

!  ---open file
   CALL open_filepars(filepars)

!  ---allocate
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')

!  ---variable attributes
   CALL set_modvars_atts(iddesc,iset,2,varatts,numvars)

!  ---write attributes
   CALL write_atts_mod(filepars,varatts,numvars,tdim=tlimit)

!  ---deallocate
   DEALLOCATE (varatts)   

ENDIF

!
!3. Write data
!-------------
!
!---date/time
IF (varid.EQ.1) CALL write_time(CDateTime,filepars)

!---nest data
CALL combine_write_stats_loc(outvals,filepars,varid+1,maxstats,nostatsglb,&
                           & nhdatprocs,indexprocs)

modfiles(iddesc,iset,2) = filepars

1000 CALL log_timer_out()


RETURN

END SUBROUTINE write_nest_data_2d

!========================================================================

SUBROUTINE write_nest_data_3d(iset,outvals,maxstats,nzdat,nhdatloc,&
                            & nostatsglb,nhdatprocs,indexprocs,iddesc)
!************************************************************************
!
! *write_nest_data_prof* Write 3-D interpolated nest data in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - update_nest_data_3d
!
! External calls -
!
! Module calls - combine_write_stats_loc, error_alloc_struc, lim_dims,
!                open_filepars, output_flag, set_modfiles_atts,
!                set_modvars_atts, write_atts_mod, write_time
!
!************************************************************************
!
USE datatypes
USE iopars
USE paralpars
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc_struc
USE inout_paral, ONLY: combine_write_stats_loc
USE inout_routines, ONLY: open_filepars, output_flag, write_atts_mod, write_time
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines
USE utility_routines, ONLY: lim_dims

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, iset, maxstats, nhdatloc, nostatsglb, nzdat
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nhdatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: indexprocs
REAL, INTENT(IN), DIMENSION(nhdatloc,nzdat) :: outvals

!
! Name        Type     Purpose
!------------------------------------------------------------------------------
!*iset*       INTEGER  Set index
!*outvals*    REAL     Output data
!*maxstats*   INTEGER  First dimension of array indexprocs
!*nzdat*      INTEGER  Number of data points in the vertical
!*nhdatloc*   INTEGER  Local number of horizontal output locations
!*nostatsglb* INTEGER  Number of stations in global array
!*nhdatprocs* INTEGER  Number of output locations per process
!*indexprocs* INTEGER  Indices of local output array within global output array
!*iddesc*     INTEGER  Data file id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
INTEGER :: numvars, tlimit
TYPE (FileParams) :: filepars
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'write_nest_data_3d'
CALL log_timer_in()

!
!1. Disable writing if needed
!----------------------------
!

filepars = modfiles(iddesc,iset,2)
flag = output_flag(CDateTime,filepars%tskips)
IF (.NOT.flag) GOTO 1000

!
!2. Write file header on first call
!----------------------------------
!

IF (master.AND.filepars%iostat.EQ.0) THEN

!  ---file attributes
   CALL set_modfiles_atts(iddesc,iset,2)
   filepars = modfiles(iddesc,iset,2)
   numvars = filepars%novars + 1
   tlimit = lim_dims(filepars%tlims)
   
!  ---open file
   CALL open_filepars(filepars)

!  ---allocate
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')

!  ---variable attributes
   CALL set_modvars_atts(iddesc,iset,2,varatts,numvars)

!  ---write attributes
   CALL write_atts_mod(filepars,varatts,numvars,tdim=tlimit)

!  ---deallocate
   DEALLOCATE (varatts)

ENDIF

!
!3. Write data
!-------------
!
!---date/time
CALL write_time(CDateTime,filepars)

!---nest data
CALL combine_write_stats_loc(outvals,filepars,2,maxstats,nostatsglb,nhdatprocs,&
                           & indexprocs,reduced=.FALSE.)

modfiles(iddesc,iset,2) = filepars

1000 CALL log_timer_out()


RETURN

END SUBROUTINE write_nest_data_3d

!========================================================================

SUBROUTINE write_nest_data_prof(iddesc,iset,outvals,maxstats,nzdat,nonst,&
                              & nhdatloc,nostatsglb,nhdatprocs,indexprocs)
!************************************************************************
!
! *write_nest_data_prof* Write 3-D interpolated nest data in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - update_nest_data_prof
!
! External calls -
!
! Module calls - combine_write_stats_loc, error_alloc_struc, lim_dims,
!                open_filepars, output_flag, set_modfiles_atts,
!                set_modvars_atts, write_atts_mod, write_time
!
!************************************************************************
!
USE datatypes
USE iopars
USE paralpars
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc_struc
USE inout_paral, ONLY: combine_write_stats_loc
USE inout_routines, ONLY: open_filepars, output_flag, write_atts_mod, write_time
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines
USE utility_routines, ONLY: lim_dims

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, iset, maxstats, nhdatloc, nonst, nostatsglb, &
                     & nzdat
INTEGER, INTENT(IN), DIMENSION(nprocs) :: nhdatprocs
INTEGER, INTENT(IN), DIMENSION(maxstats,nprocs) :: indexprocs
REAL, INTENT(IN), DIMENSION(nhdatloc,nzdat,nonst) :: outvals

!
! Name        Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*     INTEGER  Output data file id
!*iset*       INTEGER  Set index
!*outvals*    REAL     Output data
!*maxstats*   INTEGER  First dimension of array indexprocs
!*nzdat*      INTEGER  Number of data points in the vertical
!*nonst*      INTEGER  Number of "state" variables to be nested
!*nhdatloc*   INTEGER  Local number of horizontal output locations
!*nostatsglb* INTEGER  Number of stations in global array
!*nhdatprocs* INTEGER  Number of output locations per process
!*indexprocs* INTEGER  Indices of local output array within global output array
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
INTEGER :: numvars, tlimit
TYPE (FileParams) :: filepars
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'write_nest_data_prof'
CALL log_timer_in()

!
!1. Disable writing if needed
!----------------------------
!

filepars = modfiles(iddesc,iset,2)
flag = output_flag(CDateTime,filepars%tskips)
IF (.NOT.flag) GOTO 1000

!
!2. Write file header on first call
!----------------------------------
!

IF (master.AND.filepars%iostat.EQ.0) THEN

!  ---file attributes
   CALL set_modfiles_atts(iddesc,iset,2)
   filepars = modfiles(iddesc,iset,2)
   numvars = filepars%novars + 1
   tlimit = lim_dims(filepars%tlims)
   
!  ---open file
   CALL open_filepars(filepars)

!  ---allocate
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')

!  ---variable attributes
   CALL set_modvars_atts(iddesc,iset,2,varatts,numvars)

!  ---write attributes
   CALL write_atts_mod(filepars,varatts,numvars,tdim=tlimit)

!  ---deallocate
   DEALLOCATE (varatts)

ENDIF

!
!3. Write data
!-------------
!
!---date/time
CALL write_time(CDateTime,filepars)

!---nest data
CALL combine_write_stats_loc(outvals,filepars,2,maxstats,nostatsglb,nhdatprocs,&
                           & indexprocs,reduced=.TRUE.)

modfiles(iddesc,iset,2) = filepars

1000 CALL log_timer_out()


RETURN

END SUBROUTINE write_nest_data_prof

!========================================================================

SUBROUTINE write_nstgrd_abs(iset,nhdat,nzdat,xcoord,ycoord,zcoord)
!************************************************************************
!
! *write_nstgrd_abs* Write absolute geographical positions of nested open
!                    boundaries in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - define_nstgrd_locs
!
! External calls -
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
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iset, nhdat, nzdat
REAL, INTENT(IN), DIMENSION(nhdat) :: xcoord, ycoord
REAL, INTENT(IN), DIMENSION(nhdat,nzdat) :: zcoord

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iset*      INTEGER Output file number
!*nhdat*     INTEGER Number of horizontal data locations
!*nzdat*     INTEGER Number of vertical data points
!*xcoord*    REAL    X-coordinates of nest locations
!*ycoord*    REAL    Y-coordinates of nest locations
!*zcoord*    REAL    Z-coordinates of nest locations
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, SAVE :: numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_nstgrd_abs'
CALL log_timer_in()

filepars = modfiles(io_nstabs,iset,2)

!
!1. Write file header on first call
!----------------------------------
!

IF (filepars%iostat.EQ.0) THEN

!  ---file attributes
   CALL set_modfiles_atts(io_nstabs,iset,2)
   filepars = modfiles(io_nstabs,iset,2)
   
!  ---open file
   CALL open_filepars(filepars)

!  ---variable attributes
   numvars = filepars%novars
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
   CALL set_modvars_atts(io_nstabs,iset,2,varatts,numvars)

!  ---write
   CALL write_atts_mod(filepars,varatts,numvars)

!  ---initialise
   varid = 0

ENDIF

!
!2. Write data
!-------------
!

varid = varid + 1
CALL write_vars(xcoord,filepars,varid,varatts=varatts)
varid = varid + 1
CALL write_vars(ycoord,filepars,varid,varatts=varatts)
IF (nzdat.GT.0) THEN
   varid = varid + 1
   CALL write_vars(zcoord,filepars,varid,varatts=varatts)
ENDIF

!
!3. Finalise
!-----------
!

IF (varid.EQ.numvars) THEN
   CALL close_filepars(filepars)
   DEALLOCATE(varatts)
ENDIF

modfiles(io_nstabs,iset,2) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE write_nstgrd_abs

!========================================================================

SUBROUTINE write_nstgrd_rel(iset,nhdat,nzdat,nonodes,hnests,zcoord)
!************************************************************************
!
! *write_nstgrd_rel* Write relative geographical positions of nested open
!                    boundaries in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.9
!
! Description - 
!
! Reference -
!
! Calling program - define_nstgrd_locs
!
! External calls -
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, write_atts_mod, write_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE nestgrids
USE paralpars
USE syspars
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iset, nhdat, nonodes, nzdat
REAL, INTENT(IN), DIMENSION(nhdat,nzdat) :: zcoord
TYPE (HRelativeCoords), INTENT(OUT), DIMENSION(nhdat,nonodes) :: hnests

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iset*      INTEGER Output file number
!*nhdat*     INTEGER Number of horizontal data locations
!*nzdat*     INTEGER Number of vertical data points
!*nonodes*   INTEGER Number of nodes for interpolation
!*hnests*    DERIVED Relative coordinates of nest locations
!*zcoord*    REAL    Z-coordinates of nest locations
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, SAVE :: numvars, varid
INTEGER :: inode, l
TYPE (FileParams) :: filepars
REAL, DIMENSION(2,2,nhdat) :: weightsdat
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_nstgrd_rel'
CALL log_timer_in()

filepars = modfiles(io_nstrel,iset,2)

!
!1. Write file header on first call
!----------------------------------
!

IF (filepars%iostat.EQ.0) THEN

!  ---file attributes
   CALL set_modfiles_atts(io_nstrel,iset,2)
   filepars = modfiles(io_nstrel,iset,2)

!  ---open file
   CALL open_filepars(filepars)

!  ---variable attributes
   numvars = filepars%novars
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
   CALL set_modvars_atts(io_nstrel,iset,2,varatts,numvars)

!  ---write
   CALL write_atts_mod(filepars,varatts,numvars)

!  ---initialise
   varid = 0

ENDIF

!
!2. Write data
!-------------
!

inode_220: DO inode=1,nonodes
   
   varid = varid + 1
   CALL write_vars(hnests(:,inode)%icoord,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(hnests(:,inode)%jcoord,filepars,varid,varatts=varatts)
   l_221: DO l=1,nhdat
      weightsdat(:,:,l) = hnests(l,inode)%weights
   ENDDO l_221
   varid = varid + 1
   CALL write_vars(weightsdat,filepars,varid,varatts=varatts)
   IF (inode.EQ.nonodes.AND.novnst(iset).GT.0) THEN
      varid = varid + 1
      CALL write_vars(zcoord,filepars,varid,varatts=varatts)
   ENDIF

ENDDO inode_220

!
!3. Finalise
!-----------
!

IF (varid.EQ.numvars) THEN
   CALL close_filepars(filepars)
   DEALLOCATE(varatts)
ENDIF

modfiles(io_nstrel,iset,2) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE write_nstgrd_rel

!========================================================================

SUBROUTINE write_nstgrd_spec
!************************************************************************
!
! *write_nstgrd_spec* Write number of output locations and specifier
!                     arrays for nesting in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Nested_Grids.f90  V2.5
!
! Description - 
!
! Reference -
!
! Calling program - define_nstgrd_spec
!
! External calls -
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, write_atts_mod,
!                write_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE nestgrids
USE paralpars
USE switches
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_nstgrd_spec'
CALL log_timer_in()

filepars = modfiles(io_nstspc,1,2)

!
!1. Write file header
!--------------------
!
!---file attributes
CALL set_modfiles_atts(io_nstspc,1,2)
filepars = modfiles(io_nstspc,1,2)
numvars = filepars%novars

!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL set_modvars_atts(io_nstspc,1,2,varatts,numvars)

!---write
CALL write_atts_mod(filepars,varatts,numvars)

!
!2. Write data
!-------------
!
!---type of grids
CALL write_vars(nestcoords,filepars,1,varatts=varatts)

!---number of horizontal nest points
CALL write_vars(nohnstglbc,filepars,2,varatts=varatts)
CALL write_vars(nohnstglbu,filepars,3,varatts=varatts)
CALL write_vars(nohnstglbv,filepars,4,varatts=varatts)

!---number of vertical nest points
CALL write_vars(novnst,filepars,5,varatts=varatts)

!---type of 2-D output
CALL write_vars(inst2dtype,filepars,6,varatts=varatts)

!---tangential conditions
CALL write_vars(nohnstglbx,filepars,7,varatts=varatts)
CALL write_vars(nohnstglby,filepars,8,varatts=varatts)

!---sediment fractions used in nesting
varid = 8
IF (iopt_sed.GT.0) THEN
   varid = varid + 1
   CALL write_vars(nosednst,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(instsed,filepars,varid,varatts=varatts)
ENDIF

!---biological variables used in nesting
IF (iopt_biolgy.GT.0) THEN
   varid = varid + 1
   CALL write_vars(nobionst,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(instbio,filepars,varid,varatts=varatts)
ENDIF

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_nstspc,1,2) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE write_nstgrd_spec
