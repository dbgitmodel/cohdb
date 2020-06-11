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
! *MultiGrid_Initialisation* Series of routines for initialisation of the
!                            multi-grid scheme
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Initialisation.f90  V2.8.1
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - 
!
! Reference -
!
! Routines - allocate_coarse_grid__arrays, boundary_coarsening,
!            deallocate_mg_arrays, initialise_multigrid,
!            mg_domain_decomposition, mg_grid_arrays, mg_pointer_weights,
!            open_boundary_mg_arrays
!
!************************************************************************
!

!========================================================================

SUBROUTINE allocate_mg_arrays
!************************************************************************
!
! *allocate_mg_arrays* Allocate arrays on the multigrids
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Initialisation.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_multigrid
!
! External calls -
!
! Module calls - error_alloc_struc_comp
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE multigrid
USE physpars
USE syspars
USE error_routines, ONLY: error_alloc_struc_comp
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: lev, nxloc, nyloc


procname(pglev+1) = 'allocate_mg_arrays'
CALL log_timer_in()

lev_2000: DO lev=0,nomglevels-1

   nxloc = mgvars(lev)%ncloc; nyloc = mgvars(lev)%nrloc 

!
!1. Allocate
!-----------
!
!1.1 Grid arrays
!---------------
!

   ALLOCATE (mgvars(lev)%delxatc(0:nxloc,nyloc),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','delxatc',2,(/nxloc+1,nyloc/),&
                             & kndrtype,lev)

   ALLOCATE (mgvars(lev)%delxatu(0:nxloc+1,nyloc),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','delxatu',2,(/nxloc+2,nyloc/),&
                             & kndrtype,lev)

   ALLOCATE (mgvars(lev)%delxatv(nxloc,nyloc+1),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','delxatv',2,(/nxloc,nyloc+1/),&
                             & kndrtype,lev)

   ALLOCATE (mgvars(lev)%delyatc(nxloc,0:nyloc),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','delyatc',2,(/nxloc,nyloc+1/),&
                             & kndrtype,lev)

   ALLOCATE (mgvars(lev)%delyatu(nxloc+1,nyloc),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','delyatu',2,(/nxloc+1,nyloc/),&
                             & kndrtype,lev)

   ALLOCATE (mgvars(lev)%delyatv(nxloc,0:nyloc+1),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','delyatv',2,(/nxloc,nyloc+2/),&
                             & kndrtype,lev)

   ALLOCATE (mgvars(lev)%gaccatc(0:nxloc,0:nyloc),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','gaccatc',2,(/nxloc+1,nyloc+1/),&
                             & kndrtype,lev)

   ALLOCATE (mgvars(lev)%gaccatu(0:nxloc+1,nyloc),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','gaccatu',2,(/nxloc+2,nyloc/),&
                             & kndrtype,lev)

   ALLOCATE (mgvars(lev)%gaccatv(nxloc,0:nyloc+1),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','gaccatv',2,(/nxloc,nyloc+2/),&
                             & kndrtype,lev)

   ALLOCATE (mgvars(lev)%maskatc_int(nxloc,nyloc),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','maskatc_int',2,(/nxloc,nyloc/),&
                             & kndlog,lev)

   ALLOCATE (mgvars(lev)%nodeatc(0:nxloc+1,0:nyloc+1),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','nodeatc',2,(/nxloc+2,nyloc+2/),&
                             & kndrtype,lev)

   ALLOCATE (mgvars(lev)%node2du(nxloc+1,nyloc+1),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','node2du',2,(/nxloc+1,nyloc+1/),&
                             & kndrtype,lev)

   ALLOCATE (mgvars(lev)%node2dv(nxloc+1,nyloc+1),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','node2dv',2,(/nxloc+1,nyloc+1/),&
                             & kndrtype,lev)


!
!1.2 Water depths
!----------------
!

   ALLOCATE (mgvars(lev)%deptotatc(0:nxloc+1,0:nyloc+1),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','deptotatc',2,(/nxloc+2,nyloc+2/),&
                             & kndrtype,lev)

   ALLOCATE (mgvars(lev)%deptotatu(0:nxloc+1,nyloc+1),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','deptotatu',2,(/nxloc+2,nyloc+1/),&
                             & kndrtype,lev)

   ALLOCATE (mgvars(lev)%deptotatv(nxloc+1,0:nyloc+1),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','deptotatv',2,(/nxloc+1,nyloc+2/),&
                             & kndrtype,lev)

!
!1.3 Work space arrays
!---------------------
!

   ALLOCATE (mgvars(lev)%blacknode(nxloc,nyloc),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','blacknode',2,(/nxloc,nyloc/),&
                             & kndlog,lev)
   mgvars(lev)%blacknode = .FALSE.

   ALLOCATE (mgvars(lev)%correction(0:nxloc+1,0:nyloc+1),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','correction',2,(/nxloc+2,nyloc+2/),&
                             & kndrtype,lev)
   mgvars(lev)%correction = 0.0

   ALLOCATE (mgvars(lev)%precon(nxloc,nyloc),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','precon',2,(/nxloc,nyloc/),&
                             & kndrtype,lev)
   mgvars(lev)%precon = 0.0

   ALLOCATE (mgvars(lev)%rednode(nxloc,nyloc),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','rednode',2,(/nxloc,nyloc/),&
                             & kndlog,lev)
   mgvars(lev)%rednode = .FALSE.

   ALLOCATE (mgvars(lev)%residual(nxloc+1,nyloc+1),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','residual',2,(/nxloc+1,nyloc+1/),&
                             & kndrtype,lev)
   mgvars(lev)%residual = 0.0

   ALLOCATE (mgvars(lev)%rhs(nxloc,nyloc),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','rhs',2,(/nxloc,nyloc/),&
                             & kndrtype,lev)
   mgvars(lev)%rhs = 0.0

   ALLOCATE (mgvars(lev)%solvec(0:nxloc+1,0:nyloc+1),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','solvec',2,(/nxloc+2,nyloc+2/),&
                             & kndrtype,lev)
   mgvars(lev)%solvec = 0.0

   ALLOCATE (mgvars(lev)%subdiags(nxloc,nyloc,-2:2),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','subdiags',3,(/nxloc,nyloc,5/),&
                             & kndrtype,lev)
   mgvars(lev)%subdiags = 0.0

   ALLOCATE (mgvars(lev)%weights(0:nxloc+1,0:nyloc+1),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','weights',2,(/nxloc+2,nyloc+2/),&
                             & kndrtype,lev)
   mgvars(lev)%weights = 0.0

!
!2. Initialise
!-------------
!

   IF (lev.EQ.0) THEN

      mgvars(lev)%delxatc = delxatc(0:ncloc,1:nrloc)
      mgvars(lev)%delxatu = delxatu(0:ncloc+1,1:nrloc)
      mgvars(lev)%delxatv = delxatv(1:ncloc,1:nrloc+1)
      mgvars(lev)%delyatc = delyatc(1:ncloc,0:nrloc)
      mgvars(lev)%delyatu = delyatu(1:ncloc+1,1:nrloc)
      mgvars(lev)%delyatv = delyatv(1:ncloc,0:nrloc+1)
      mgvars(lev)%gaccatc = gaccatc(0:ncloc,0:nrloc)
      mgvars(lev)%gaccatu = gaccatu
      mgvars(lev)%gaccatv = gaccatv
      mgvars(lev)%maskatc_int = maskatc_int
      mgvars(lev)%nodeatc = nodeatc(0:ncloc+1,0:nrloc+1)
      mgvars(lev)%node2du = node2du(1:ncloc+1,1:nrloc+1)
      mgvars(lev)%node2dv = node2dv(1:ncloc+1,1:nrloc+1)
      mgvars(lev)%deptotatc = deptotatc(0:ncloc+1,0:nrloc+1)
      mgvars(lev)%deptotatu = deptotatu(0:ncloc+1,1:nrloc+1)
      mgvars(lev)%deptotatv = deptotatv(1:ncloc+1,0:nrloc+1)

   ELSE

      mgvars(lev)%delxatc = 0.0; mgvars(lev)%delxatu = 0.0
      mgvars(lev)%delxatv = 0.0; mgvars(lev)%delyatc = 0.0
      mgvars(lev)%delyatu = 0.0; mgvars(lev)%delyatv = 0.0
      mgvars(lev)%gaccatc = 0.0; mgvars(lev)%gaccatu = 0.0
      mgvars(lev)%gaccatv = 0.0; mgvars(lev)%maskatc_int = .FALSE.
      mgvars(lev)%nodeatc = 0; mgvars(lev)%node2du = 0
      mgvars(lev)%node2dv = 0; mgvars(lev)%deptotatc = 0.0
      mgvars(lev)%deptotatu = 0.0; mgvars(lev)%deptotatv = 0.0

   ENDIF

ENDDO lev_2000

CALL log_timer_out()


RETURN

END SUBROUTINE allocate_mg_arrays

!========================================================================

SUBROUTINE boundary_coarsening
!************************************************************************
!
! *boundary_coarsening* Parameters and arrays for coarse grid open boundaries
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Initialisation.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_multigrid
!
! External calls - error_alloc, error_alloc_struc_comp
!
! Module calls -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE multigrid
USE obconds
USE paralpars
USE physpars
USE syspars
USE error_routines, ONLY: error_alloc, error_alloc_struc_comp
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag, sout, west
INTEGER :: i, iiF, iloc, ityp, i_C, i_F, j, jjF, j_C, j_F, levC, levF, ncC, &
         & ncF, nobuC, nobuF, nobvC, nobvF, nrC, nrF
REAL :: eps
LOGICAL, DIMENSION(nobu+nobv) :: obWS
INTEGER, DIMENSION(nobu+nobv) :: iob, iloczob, ityp2dob, job
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: boundaryfaceU, boundaryfaceV, &
                                            & seapointC, seapointF


procname(pglev+1) = 'boundary_coarsening'
CALL log_timer_in()

!
!1. Fine grid
!------------
!
!1.1 Sea mask
!------------
!

IF (nomglevels.GT.1) THEN
   eps = 0.00001*ABS(depmean_flag)
   ALLOCATE(seapointF(0:nc+1,0:nr+1),STAT=errstat)
   CALL error_alloc('seapointF',2,(/nc+2,nr+2/),kndlog)
   seapointF = seapointglb
ENDIF

! 
!1.2 Allocate/define open boundary locations and arrays
!------------------------------------------------------
!
!---U-nodes
mgvars(0)%nobu = nobu; mgvars(0)%nobv = nobv;
ALLOCATE (mgvars(0)%iobu(nobu),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','iobu',1,(/nobu/),kndint,0)
mgvars(0)%iobu = iobu
ALLOCATE (mgvars(0)%jobu(nobu),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','jobu',1,(/nobu/),kndint,0)
mgvars(0)%jobu = jobu
ALLOCATE (mgvars(0)%ityp2dobu(nobu),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','ityp2dobu',1,(/nobu/),kndint,0)
mgvars(0)%ityp2dobu = ityp2dobu
ALLOCATE (mgvars(0)%iloczobu(nobu),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','iloczobu',1,(/nobu/),kndint,0)
mgvars(0)%iloczobu = iloczobu
ALLOCATE (mgvars(0)%westobu(nobu),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','westobu',1,(/nobu/),kndlog,0)
mgvars(0)%westobu = westobu
ALLOCATE (mgvars(0)%indexobu(nobu),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','indexobu',1,(/nobu/),kndint,0)
mgvars(0)%indexobu = indexobu
ALLOCATE (mgvars(0)%indexobuprocs(nobu,nprocs),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','indexobuprocs',2,(/nobu,nprocs/),&
                             & kndint,0)
mgvars(0)%indexobuprocs = indexobuprocs

!---V-nodes
ALLOCATE (mgvars(0)%iobv(nobv),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','iobv',1,(/nobv/),kndint,0)
mgvars(0)%iobv = iobv
ALLOCATE (mgvars(0)%jobv(nobv),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','jobv',1,(/nobv/),kndint,0)
mgvars(0)%jobv = jobv
ALLOCATE (mgvars(0)%ityp2dobv(nobv),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','ityp2dobv',1,(/nobv/),kndint,0)
mgvars(0)%ityp2dobv = ityp2dobv
ALLOCATE (mgvars(0)%iloczobv(nobv),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','iloczobv',1,(/nobv/),kndint,0)
mgvars(0)%iloczobv = iloczobv
ALLOCATE (mgvars(0)%soutobv(nobv),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','soutobv',1,(/nobv/),kndlog,0)
mgvars(0)%soutobv = soutobv
ALLOCATE (mgvars(0)%indexobv(nobv),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','indexobv',1,(/nobv/),kndint,0)
mgvars(0)%indexobv = indexobv 
ALLOCATE (mgvars(0)%indexobvprocs(nobv,nprocs),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','indexobvprocs',2,(/nobv,nprocs/),&
                          & kndint,0)
mgvars(0)%indexobvprocs = indexobvprocs

!
!2. Coarse grids
!---------------
!

levC_2000: DO levC=1,nomglevels-1

!
!2.1 Sea mask on coarse grid
!---------------------------
!

   levF = levC - 1
   ncF = mgvars(levF)%nc; ncC = mgvars(levC)%nc
   nrF = mgvars(levF)%nr; nrC = mgvars(levC)%nr

   ALLOCATE(seapointC(0:ncC+1,0:nrC+1),STAT=errstat)
   CALL error_alloc('seapointC',2,(/ncC+2,nrC+2/),kndlog)
   seapointC = .FALSE.
   seapointC(1:ncC,1:nrC) = seapointF(1:ncF:2,1:nrF:2).OR. &
                          & seapointF(1:ncF:2,2:nrF+1:2).OR. &
                          & seapointF(2:ncF+1:2,1:nrF:2).OR. &
                          & seapointF(2:ncF+1:2,2:nrF+1:2)

!
!2.2 Marked locations of open and coastal velocity boundaries
!------------------------------------------------------------
!
!  ---U-nodes
   ALLOCATE (boundaryfaceU(ncC,nrC-1),STAT=errstat)
   boundaryfaceU = seapointC(0:ncC-1,1:nrC-1).NEQV.seapointC(1:ncC,1:nrC-1)

!  ---V-nodes
   ALLOCATE (boundaryfaceV(ncC-1,nrC),STAT=errstat)
   boundaryfaceV = seapointC(1:ncC-1,0:nrC-1).NEQV.seapointC(1:ncC-1,1:nrC)

!
!2.3 Check marked locations for open boundaries
!----------------------------------------------
!
!2.3.1 U-nodes
!-------------
!

   nobuF = mgvars(levF)%nobu
   nobuC = 0

   IF (nobuF.GT.0) THEN
      ityp2dob = 0; iloczob  = 0

!     ---loop over all coarse U-faces
      i_C_231: DO i_C=1,ncC
      j_C_231: DO j_C=1,nrC-1
         IF (boundaryfaceU(i_C,j_C)) THEN
            i_F = 2*i_C-1; j_F = 2*j_C-1
            flag = .FALSE.; ityp = 0; iloc = 0

!           ---loop over fine U-faces
            iiF_2311: DO iiF=1,nobuF
               i = mgvars(levF)%iobu(iiF); j = mgvars(levF)%jobu(iiF)
               west = mgvars(levF)%westobu(iiF)
               IF (j.EQ.j_F.OR.j.EQ.(j_F+1)) THEN
                  IF (i.EQ.i_F.OR.(i.EQ.(i_F+1).AND.west).OR.&
                   & (i.EQ.(i_F-1).AND.(.NOT.west))) THEN
                     flag = .TRUE.
                     ityp = MAX(ityp,mgvars(levF)%ityp2dobu(iiF))
                     iloc = MAX(iloc,mgvars(levF)%iloczobu(iiF))
                  ENDIF
               ENDIF
            ENDDO iiF_2311

            IF (flag) THEN
               nobuC = nobuC + 1
               iob(nobuC) = i_C; job(nobuC) = j_C
               ityp2dob(nobuC) = ityp; iloczob(nobuC) = iloc
               obWS(nobuC) = .NOT.seapointC(i_C-1,j_C).OR.i_C.EQ.1
            ENDIF
         ENDIF
         
      ENDDO j_C_231
      ENDDO i_C_231

   ENDIF

!  ---allocate open boundary arrays
   mgvars(levC)%nobu = nobuC
   ALLOCATE (mgvars(levC)%iobu(nobuC),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','iobu',1,(/nobuC/),kndint,levC)
   ALLOCATE (mgvars(levC)%jobu(nobuC),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','jobu',1,(/nobuC/),kndint,levC)
   ALLOCATE (mgvars(levC)%ityp2dobu(nobuC),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','ityp2dobu',1,(/nobuC/),kndint,levC)
   ALLOCATE (mgvars(levC)%iloczobu(nobuC),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','iloczobu',1,(/nobuC/),kndint,levC)
   ALLOCATE (mgvars(levC)%westobu(nobuC),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','westobu',1,(/nobuC/),kndlog,levC)
   ALLOCATE (mgvars(levC)%indexobu(nobuC),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','indexobu',1,(/nobuC/),kndint,levC)
   ALLOCATE (mgvars(levC)%indexobuprocs(nobuC,nprocs),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','indexobuprocs',2,(/nobuC,nprocs/),&
                             & kndint,levC)

!  ---define
   IF (nobuC.GT.0) THEN
      mgvars(levC)%iobu = iob(1:nobuC)
      mgvars(levC)%jobu = job(1:nobuC)
      mgvars(levC)%ityp2dobu = ityp2dob(1:nobuC)
      mgvars(levC)%iloczobu = iloczob(1:nobuC)
      mgvars(levC)%westobu = obWS(1:nobuC)
   ENDIF

!
!2.3.2 V-nodes
!-------------
!

   nobvF = mgvars(levF)%nobv
   nobvC = 0
   IF (nobvF.GT.0) THEN
      ityp2dob = 0; iloczob  = 0

!     ---loop over all coarse V-faces
      i_C_232: DO i_C=1,ncC-1
      j_C_232: DO j_C=1,nrC
         IF (boundaryfaceV(i_C,j_C)) THEN
            i_F = 2*i_C-1; j_F = 2*j_C-1
            flag = .FALSE.; ityp = 0; iloc = 0

!           ---loop over fine V-faces
            jjF_2321: DO jjF=1,nobvF
               i = mgvars(levF)%iobv(jjF); j = mgvars(levF)%jobv(jjF)
               sout = mgvars(levF)%soutobv(jjF)
               IF (i.EQ.i_F.OR.i.EQ.(i_F+1)) THEN
                  IF (j.EQ.j_F.OR.(j.EQ.(j_F+1).AND.sout).OR.&
                   & (j.EQ.(j_F-1).AND.(.NOT.sout))) THEN
                     flag = .TRUE.
                     ityp = MAX(ityp,mgvars(levF)%ityp2dobv(jjF))
                     iloc = MAX(iloc,mgvars(levF)%iloczobv(jjF))
                  ENDIF
               ENDIF
            ENDDO jjF_2321
            
            IF (flag) THEN
               nobvC = nobvC + 1
               iob(nobvC) = i_C; job(nobvC) = j_C
               ityp2dob(nobvC) = ityp; iloczob(nobvC) = iloc
               obWS(nobvC) = .NOT.seapointC(i_C,j_C-1).OR.j_C.EQ.1 
            ENDIF
         ENDIF
         
      ENDDO j_C_232
      ENDDO i_C_232

   ENDIF

!  ---allocate and define ob arrays
   mgvars(levC)%nobv = nobvC
   ALLOCATE (mgvars(levC)%iobv(nobvC),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','iobv',1,(/nobvC/),kndint,levC)
   ALLOCATE (mgvars(levC)%jobv(nobvC),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','jobv',1,(/nobvC/),kndint,levC)
   ALLOCATE (mgvars(levC)%ityp2dobv(nobvC),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','ityp2dobv',1,(/nobvC/),kndint,levC)
   ALLOCATE (mgvars(levC)%iloczobv(nobvC),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','iloczobv',1,(/nobvC/),kndint,levC)
   ALLOCATE (mgvars(levC)%soutobv(nobvC),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','soutobv',1,(/nobvC/),kndint,levC)
   ALLOCATE (mgvars(levC)%indexobv(nobvC),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','indexobv',1,(/nobvC/),kndint,levC)
   ALLOCATE (mgvars(levC)%indexobvprocs(nobvC,nprocs),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','indexobvprocs',2,(/nobvC,nprocs/),&
                             & kndint,levC)
      
!  ---define
   IF (nobvC.GT.0) THEN
      mgvars(levC)%iobv = iob(1:nobvC)
      mgvars(levC)%jobv = job(1:nobvC)
      mgvars(levC)%ityp2dobv = ityp2dob(1:nobvC)
      mgvars(levC)%iloczobv = iloczob(1:nobvC)
      mgvars(levC)%soutobv = obWS(1:nobvC)
   ENDIF

!
!2.4 Reset sea masks
!-------------------
!

   DEALLOCATE (seapointF)
   IF (levC.LT.nomglevels-1) THEN
      ALLOCATE(seapointF(0:ncC+1,0:nrC+1),STAT=errstat)
      CALL error_alloc('seapointF',2,(/ncC+2,nrC+2/),kndlog)
      seapointF = seapointC
   ENDIF
   DEALLOCATE (boundaryfaceU,boundaryfaceV,seapointC)

ENDDO levC_2000


CALL log_timer_out()

RETURN

END SUBROUTINE boundary_coarsening

!========================================================================

SUBROUTINE deallocate_mg_arrays
!************************************************************************
!
! *deallocate_mg_arrays* Deallocate/nullify multigrid arrays
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Initialisation.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_multigrid
!
! Module calls -
!
!************************************************************************
!
USE iopars
USE multigrid
USE physpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: lev


procname(pglev+1) = 'deallocate_mg_arrays'
CALL log_timer_in()

lev_100: DO lev=0,nomglevels-1

!  ---grid arrays
   DEALLOCATE (mgvars(lev)%westobu,mgvars(lev)%soutobv,mgvars(lev)%iobu,&
             & mgvars(lev)%iobv,mgvars(lev)%jobu,mgvars(lev)%jobv,&
             & mgvars(lev)%iobuloc,mgvars(lev)%iobvloc,mgvars(lev)%jobuloc,&
             & mgvars(lev)%jobvloc,mgvars(lev)%indexobu,mgvars(lev)%indexobv)
   DEALLOCATE (mgvars(lev)%delxatc,mgvars(lev)%delxatu,mgvars(lev)%delxatv,&
             & mgvars(lev)%delyatc,mgvars(lev)%delyatu,mgvars(lev)%delyatv,&
             & mgvars(lev)%gaccatc,mgvars(lev)%gaccatu,mgvars(lev)%gaccatv,&
             & mgvars(lev)%maskatc_int,mgvars(lev)%nodeatc,mgvars(lev)%node2du,&
             & mgvars(lev)%node2dv)    

!  ---water depths
   DEALLOCATE (mgvars(lev)%deptotatc,mgvars(lev)%deptotatu,&
             & mgvars(lev)%deptotatv)

!  ---open boundary conditions
   DEALLOCATE (mgvars(lev)%iloczobu,mgvars(lev)%iloczobv,mgvars(lev)%ityp2dobu,&
             & mgvars(lev)%ityp2dobv)

!  ---domain decomposition
   DEALLOCATE (mgvars(lev)%ncprocs,mgvars(lev)%nrprocs,mgvars(lev)%nc1procs,&
             & mgvars(lev)%nr1procs,mgvars(lev)%nc2procs,mgvars(lev)%nr2procs,&
             & mgvars(lev)%nobuprocs,mgvars(lev)%nobvprocs,&
             & mgvars(lev)%indexobuprocs,mgvars(lev)%indexobvprocs,&
             & mgvars(lev)%halocomms)

!  ---work space arrays
   DEALLOCATE (mgvars(lev)%blacknode,mgvars(lev)%correction,mgvars(lev)%precon,&
             & mgvars(lev)%rednode,mgvars(lev)%residual,mgvars(lev)%rhs,&
             & mgvars(lev)%solvec,mgvars(lev)%subdiags,mgvars(lev)%weights)

ENDDO lev_100

CALL log_timer_out()


RETURN

END SUBROUTINE deallocate_mg_arrays

!========================================================================

SUBROUTINE initialise_multigrid
!************************************************************************
!
! *initialise_multigrid* Allocate and initialise parameters annd arrays for
!                        the multi-grid procedure
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Initialisation.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - allocate_coarse_grid__arrays, assign_fine_grid_arrays,
!                  boundary_coarsening, mg_domain_decomposition,
!                  mg_grid_arrays, mg_pointers_weights, mg_water_depths,
!                  open_boundary_mg_arrays
!
! Module calls -
!
!************************************************************************
!
USE gridpars
USE iopars
USE multigrid
USE physpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: lev, npcc


procname(pglev+1) = 'initialise_multigrid'
CALL log_timer_in(npcc)

!
!1. Global grid dimensions
!-------------------------
!

mgvars(0)%nc = nc; mgvars(0)%nr = nr
lev_210: DO lev=1,nomglevels-1
   mgvars(lev)%nc = (mgvars(lev-1)%nc+1)/2
   mgvars(lev)%nr = (mgvars(lev-1)%nr+1)/2
ENDDO lev_210

!
!2. Open boundary parameters and arrays
!--------------------------------------
!

CALL boundary_coarsening

!
!3. Domain decomposition
!-----------------------
!

CALL mg_domain_decomposition

!
!4. Allocate arrays
!------------------
!

CALL allocate_mg_arrays

!
!5. Grid arrays
!--------------
!

IF (nomglevels.GT.1) CALL mg_grid_arrays

!
!6. Open boundary locations
!--------------------------
!

CALL open_boundary_mg_arrays

!
!7. Pointers and weight factors
!------------------------------
!

CALL mg_pointers_weights 

!
!8. Water depths
!---------------
!

CALL mg_water_depths

CALL log_timer_out(npcc,itm_mg)


RETURN

END SUBROUTINE initialise_multigrid

!========================================================================

SUBROUTINE mg_domain_decomposition
!************************************************************************
!
! *mg_domain_decomposition* Decompose multigrid sub-domains for parallel
!                           processing 
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Initialisation.f90  V2.10
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - define_halo_comms, mg_regular_partition
!
! Module calls - error_abort, error_alloc_struc_comp, error_limits_arr
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE multigrid
USE paralpars
USE physpars
USE switches
USE syspars
USE error_routines, ONLY: error_abort, error_alloc_struc_comp, error_limits_arr
USE time_routines, ONLY: log_timer_in, log_timer_out


IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iproc, lev


procname(pglev+1) = 'mg_domain_decomposition'
CALL log_timer_in()


lev_1000: DO lev=0,nomglevels-1

!
!.1 Allocate arrays
!------------------
!

   ALLOCATE(mgvars(lev)%ncprocs(nprocs),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','ncprocs',1,(/nprocs/),kndint,lev)
   mgvars(lev)%ncprocs = MERGE(ncprocs,0,lev.EQ.0)

   ALLOCATE(mgvars(lev)%nrprocs(nprocs),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','nrprocs',1,(/nprocs/),kndint,lev)
   mgvars(lev)%nrprocs = MERGE(nrprocs,0,lev.EQ.0)

   ALLOCATE (mgvars(lev)%nobuprocs(nprocs),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','nobuprocs',1,(/nprocs/),kndint,lev)
   mgvars(lev)%nobuprocs = MERGE(nobuprocs,0,lev.EQ.0)

   ALLOCATE (mgvars(lev)%nobvprocs(nprocs),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','nobvprocs',1,(/nprocs/),kndint,lev)
   mgvars(lev)%nobvprocs = MERGE(nobvprocs,0,lev.EQ.0)

   ALLOCATE (mgvars(lev)%halocomms(MaxHaloComms),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','MaxHaloComms',1,(/MaxHaloComms*20/),&
                             & kndint,lev)
   IF (lev.EQ.0) mgvars(0)%halocomms = halocomms

!
!2. Domain partition
!-------------------
!

   IF (lev.EQ.0) THEN

      mgvars(0)%nc1procs = nc1procs; mgvars(0)%nc2procs= nc2procs
      mgvars(0)%nr1procs = nr1procs; mgvars(0)%nr2procs= nr2procs
      
   ELSE

!
!2.1 Serial case
!---------------
!

      IF (iopt_MPI.EQ.0) THEN
         mgvars(lev)%nc1procs = 1; mgvars(lev)%nc2procs = mgvars(lev)%nc
         mgvars(lev)%nr1procs = 1; mgvars(lev)%nr2procs = mgvars(lev)%nr

!
!2.2 Regular partition
!---------------------
!

      ELSEIF (iopt_MPI_partit.EQ.1) THEN
         CALL mg_regular_partition(lev)
      ENDIF

!
!3. Dimensions of local arrays
!-----------------------------
!

      mgvars(lev)%ncprocs = mgvars(lev)%nc2procs - mgvars(lev)%nc1procs + 1
      mgvars(lev)%nrprocs = mgvars(lev)%nr2procs - mgvars(lev)%nr1procs + 1

!
!4. Check values
!---------------
!

      IF (iopt_MPI.EQ.1.AND.master) THEN
         iproc_411: DO iproc=1,nprocs
            CALL error_limits_arr(mgvars(lev)%nc1procs(iproc),&
                            & 'nc1procs',1,mgvars(lev)%nc,2,(/lev,iproc/))
            CALL error_limits_arr(mgvars(lev)%nc2procs(iproc),&
                            & 'nc2procs',1,mgvars(lev)%nc,2,(/lev,iproc/))
            CALL error_limits_arr(mgvars(lev)%nr1procs(iproc),&
                            & 'nr1procs',1,mgvars(lev)%nr,2,(/lev,iproc/))
            CALL error_limits_arr(mgvars(lev)%nr2procs(iproc),&
                            & 'nr2procs',1,mgvars(lev)%nr,2,(/lev,iproc/))
            CALL error_limits_arr(mgvars(lev)%ncprocs(iproc),&
                            & 'ncprocs',2,mgvars(lev)%nc,2,(/lev,iproc/))
            CALL error_limits_arr(mgvars(lev)%nrprocs(iproc),&
                            & 'nrprocs',2,mgvars(lev)%nr,2,(/lev,iproc/))
         ENDDO iproc_411
         CALL error_abort('mg_domain_decomposition',ierrno_inival)
      ENDIF
   ENDIF

!
!5. Values on local grid
!-----------------------
!

   IF (lev.EQ.0) THEN
      mgvars(0)%nc1loc = nc1loc; mgvars(0)%nr1loc = nr1loc
      mgvars(0)%nc2loc = nc2loc; mgvars(0)%nr2loc = nr2loc
      mgvars(0)%ncloc = ncloc; mgvars(0)%nrloc = nrloc
   ELSE
      iproc_510: DO iproc=1,nprocs
         IF (idloc.EQ.idprocs(iproc)) THEN
            mgvars(lev)%nc1loc = mgvars(lev)%nc1procs(iproc)
            mgvars(lev)%nc2loc = mgvars(lev)%nc2procs(iproc)
            mgvars(lev)%ncloc  = mgvars(lev)%ncprocs(iproc)
            mgvars(lev)%nr1loc = mgvars(lev)%nr1procs(iproc)
            mgvars(lev)%nr2loc = mgvars(lev)%nr2procs(iproc)
            mgvars(lev)%nrloc  = mgvars(lev)%nrprocs(iproc)
         ENDIF
      ENDDO iproc_510
   ENDIF

!
!6. Define parameters for communication
!--------------------------------------
!

   IF (lev.GT.0.AND.iopt_MPI.EQ.1) THEN
      CALL define_halo_comms(lev)
   ENDIF

ENDDO lev_1000

CALL log_timer_out()


RETURN

END SUBROUTINE mg_domain_decomposition

!========================================================================

SUBROUTINE mg_grid_arrays
!************************************************************************
!
! *mg_grid_arrays* Grid spacings and other arrays related to the coarse grids
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Initialisation.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_multigrid
!
! Module calls - distribute_mod, error_alloc, exchange_mod, grid_spacings_curv,
!                grid_spacings_rectang
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE modids
USE multigrid
USE physpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: construct_rectgrid, grid_spacings_curv, &
                       & grid_spacings_rectang
USE paral_comms, ONLY: distribute_mod, exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: levC, levF, ncC, ncF, nclocC, nc1locC, nc2locC, nrC, nrF, nrlocC, &
         & nr1locC, nr2locC
INTEGER, DIMENSION(2) :: lbounds
INTEGER, DIMENSION(4) :: nhdist, nhexch
REAL :: xstart, ystart
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: gdelxglbC, gdelxglbF, gdelyglbC, &
                                       & gdelyglbF
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: delglb, gxcoordglbC, gxcoordglbF, &
                      & gycoordglbC, gycoordglbF, xcoord, ycoord, ycoordloc


procname(pglev+1) = 'mg_grid_arrays'
CALL log_timer_in()

!
!1. Grid arrays on the fine grid
!-------------------------------
!

ALLOCATE (gxcoordglbF(0:nc+1,0:nr+1),STAT=errstat)
CALL error_alloc('gxcoordglbF',2,(/nc+2,nr+2/),kndrtype)

ALLOCATE (gycoordglbF(0:nc+1,0:nr+1),STAT=errstat)
CALL error_alloc('gycoordglbF',2,(/nc+2,nr+2/),kndrtype)

ALLOCATE(gdelxglbF(0:nc+1),STAT=errstat)
CALL error_alloc('gdelxglbF',1,(/nc+2/),kndrtype)

ALLOCATE(gdelyglbF(0:nr+1),STAT=errstat)
CALL error_alloc('gdelyglbF',1,(/nr+2/),kndrtype)

IF (iopt_grid_htype.LE.2.AND.surfacegrids(igrd_model,1)%rotated) THEN
   xstart = surfacegrids(igrd_model,1)%x0rot - gdelxglb(0)
   ystart = surfacegrids(igrd_model,1)%y0rot - gdelyglb(0)
   CALL construct_rectgrid(xstart,ystart,gdelxglb,gdelyglb,&
                         & gxcoordglbF,gycoordglbF,nc+2,nr+2)
ELSE
   gxcoordglbF = gxcoordglb; gycoordglbF = gycoordglb
ENDIF

gdelxglbF = gdelxglb; gdelyglbF = gdelyglb
 
!
!2. Coarse grids
!---------------
!

levC_2000: DO levC=1,nomglevels-1

   levF = levC - 1
   ncF = mgvars(levF)%nc; ncC = mgvars(levC)%nc
   nrF = mgvars(levF)%nr; nrC = mgvars(levC)%nr
   nclocC = mgvars(levC)%ncloc; nrlocC = mgvars(levC)%nrloc
   nc1locC = mgvars(levC)%nc1loc; nr1locC = mgvars(levC)%nr1loc
   nc2locC = mgvars(levC)%nc2loc; nr2locC = mgvars(levC)%nr2loc

!
!2.1 Allocate arrays
!-------------------
!
!  ---grid coordinates
   ALLOCATE (gxcoordglbC(0:ncC+1,0:nrC+1),STAT=errstat)
   CALL error_alloc('gxcoordglbC',2,(/ncC+2,nrC+2/),kndrtype)
   ALLOCATE (gycoordglbC(0:ncC+1,0:nrC+1),STAT=errstat)
   CALL error_alloc('gycoordglbC',2,(/ncC+2,nrC+2/),kndrtype)

!  ---grid spacings
   ALLOCATE(gdelxglbC(0:ncC+1),STAT=errstat)
   CALL error_alloc('gdelxglbC',1,(/ncC+2/),kndrtype)
   ALLOCATE(gdelyglbC(0:nrC+1),STAT=errstat)
   CALL error_alloc('gdelyglbC',1,(/nrC+2/),kndrtype)

!  ---work space arrays
   ALLOCATE (delglb(0:ncC+1,0:nrC+1),STAT=errstat)
   CALL error_alloc('delglb',2,(/ncC+2,nrC+2/),kndrtype)
   delglb = 0.0
   ALLOCATE (xcoord(0:ncC+1,0:nrC+1),STAT=errstat)
   CALL error_alloc('xcoord',2,(/ncC+2,nrC+2/),kndrtype)
   ALLOCATE (ycoord(0:ncC+1,0:nrC+1),STAT=errstat)
   CALL error_alloc('ycoord',2,(/ncC+2,nrC+2/),kndrtype)
   ALLOCATE (ycoordloc(nclocC,nrlocC),STAT=errstat)
   CALL error_alloc('ycoordloc',2,(/nclocC,nrlocC/),kndrtype)

!
!2.2 Horizontal grid coordinates
!-------------------------------
!
!  ---interior points
   gxcoordglbC(1:ncC,1:nrC) =  gxcoordglbF(1:ncF:2,1:nrF:2)
   gycoordglbC(1:ncC,1:nrC) =  gycoordglbF(1:ncF:2,1:nrF:2)

!  ---extend X-coordinates outside domain
   gxcoordglbC(0,1:nrC) = gxcoordglbC(1,1:nrC) + &
                       & (gxcoordglbC(1,1:nrC)-gxcoordglbC(2,1:nrC))
   gxcoordglbC(ncC+1,1:nrC) = gxcoordglbC(ncC,1:nrC) + &
                           & (gxcoordglbC(ncC,1:nrC)-gxcoordglbC(ncC-1,1:nrC))
   gxcoordglbC(1:ncC,0) = gxcoordglbC(1:ncC,1) + &
                       & (gxcoordglbC(1:ncC,1)-gxcoordglbC(1:ncC,2))
   gxcoordglbC(1:ncC,nrC+1) = gxcoordglbC(1:ncC,nrC) + &
                           & (gxcoordglbC(1:ncC,nrC)-gxcoordglbC(1:ncC,nrC-1))
   gxcoordglbC(0,0) = gxcoordglbC(1,0) + (gxcoordglbC(1,0)-gxcoordglbC(2,0))
   gxcoordglbC(0,nrC+1) = gxcoordglbC(1,nrC+1) + &
                      & (gxcoordglbC(1,nrC+1)-gxcoordglbC(2,nrC+1))
   gxcoordglbC(ncC+1,0) = gxcoordglbC(ncC+1,1) + &
                       & (gxcoordglbC(ncC+1,1)-gxcoordglbC(ncC+1,2))
   gxcoordglbC(ncC+1,nrC+1) = gxcoordglbC(ncC,nrC+1) + &
                           & (gxcoordglbC(ncC,nrC+1)-gxcoordglbC(ncC-1,nrC+1))

!  ---extend Y-coordinates outside domain
   gycoordglbC(0,1:nrC) = gycoordglbC(1,1:nrC) + &
                       & (gycoordglbC(1,1:nrC)-gycoordglbC(2,1:nrC))
   gycoordglbC(ncC+1,1:nrC) = gycoordglbC(ncC,1:nrC) + &
                           & (gycoordglbC(ncC,1:nrC)-gycoordglbC(ncC-1,1:nrC))
   gycoordglbC(1:ncC,0) = gycoordglbC(1:ncC,1) + &
                       & (gycoordglbC(1:ncC,1)-gycoordglbC(1:ncC,2))
   gycoordglbC(1:ncC,nrC+1) = gycoordglbC(1:ncC,nrC) + &
                           & (gycoordglbC(1:ncC,nrC)-gycoordglbC(1:ncC,nrC-1))
   gycoordglbC(0,0) = gycoordglbC(1,0) + (gycoordglbC(1,0)-gycoordglbC(2,0))
   gycoordglbC(0,nrC+1) = gycoordglbC(1,nrC+1) + &
                      & (gycoordglbC(1,nrC+1)-gycoordglbC(2,nrC+1))
   gycoordglbC(ncC+1,0) = gycoordglbC(ncC+1,1) + &
                       & (gycoordglbC(ncC+1,1)-gycoordglbC(ncC+1,2))
   gycoordglbC(ncC+1,nrC+1) = gycoordglbC(ncC,nrC+1) + &
                           & (gycoordglbC(ncC,nrC+1)-gycoordglbC(ncC-1,nrC+1))

!  ---grid spacings
   gdelxglbC(1:ncC) = gdelxglbF(1:ncF:2) + gdelxglbF(2:ncF+1:2)
   gdelxglbC(0) = gdelxglbC(1); gdelxglbC(ncC+1) = gdelxglbC(ncC)
   gdelyglbC(1:nrC) = gdelyglbF(1:nrF:2) + gdelyglbF(2:nrF+1:2)
   gdelyglbC(0) = gdelyglbC(1); gdelyglbC(nrC+1) = gdelyglbC(nrC)

!
!2.3 Horizontal grid spacings
!----------------------------
!
!2.3.1 X-direction
!-----------------
!
!  ---V-nodes
   IF (iopt_grid_htype.LE.2) THEN
      CALL grid_spacings_rectang(delglb(1:ncC,1:nrC+1),gdelxglbC(1:ncC),&
                               & gycoordglbC(1,1:nrC+1),2,3,'X')
   ELSE
      CALL grid_spacings_curv(delglb(1:ncC,1:nrC+1),gxcoordglbC(1:ncC,1:nrC+1),&
                            & gxcoordglbC(2:ncC+1,1:nrC+1),&
                            & gycoordglbC(1:ncC,1:nrC+1),&
                            & gycoordglbC(2:ncC+1,1:nrC+1),2,3)
   ENDIF
   lbounds = 1; nhdist = (/0,0,0,1/)
   CALL distribute_mod(delglb(1:ncC,1:nrC+1),mgvars(levC)%delxatv,lbounds,&
                     & nhdist,iarr_delxatv,0.0,mglevel=levC)

!  ---C-nodes
   IF (iopt_grid_htype.LE.2) THEN
      ycoord(1,1:nrC) = 0.5*(gycoordglbC(1,1:nrC)+gycoordglbC(1,2:nrC+1))
      CALL grid_spacings_rectang(delglb(0:ncC,1:nrC),gdelxglbC(0:ncC),&
                               & ycoord(1,1:nrC),2,3,'X')
   ELSE
      xcoord(0:ncC+1,1:nrC) = 0.5*(gxcoordglbC(0:ncC+1,1:nrC)+&
                                 & gxcoordglbC(0:ncC+1,2:nrC+1))
      ycoord(0:ncC+1,1:nrC) = 0.5*(gycoordglbC(0:ncC+1,1:nrC)+&
                                 & gycoordglbC(0:ncC+1,2:nrC+1))
      CALL grid_spacings_curv(delglb(0:ncC,1:nrC),xcoord(0:ncC,1:nrC),&
                            & xcoord(1:ncC+1,1:nrC),ycoord(0:ncC,1:nrC),&
                            & ycoord(1:ncC+1,1:nrC),2,3)
   ENDIF
   lbounds = (/0,1/); nhdist = (/1,0,0,0/)
   CALL distribute_mod(delglb(0:ncC,1:nrC),mgvars(levC)%delxatc,lbounds,nhdist,&
                     & iarr_delxatc,0.0,mglevel=levC)

!  ---U-nodes
   IF (iopt_grid_htype.LE.2) THEN
      ycoord(1,1:nrC) = 0.5*(gycoordglbC(1,1:nrC)+gycoordglbC(1,2:nrC+1))
      CALL grid_spacings_rectang(delglb(1:ncC,1:nrC),gdelxglbC(1:ncC),&
                               & ycoord(1,1:nrC),2,3,'X')
   ELSE
      xcoord(0:ncC,1:nrC) = 0.25*(gxcoordglbC(0:ncC,1:nrC)+&
                                & gxcoordglbC(1:ncC+1,1:nrC)+&
                                & gxcoordglbC(0:ncC,2:nrC+1)+&
                                & gxcoordglbC(1:ncC+1,2:nrC+1))
      ycoord(0:ncC,1:nrC) = 0.25*(gycoordglbC(0:ncC,1:nrC)+&
                                & gycoordglbC(1:ncC+1,1:nrC)+&
                                & gycoordglbC(0:ncC,2:nrC+1)+&
                                & gycoordglbC(1:ncC+1,2:nrC+1))
      CALL grid_spacings_curv(delglb(1:ncC,1:nrC),xcoord(0:ncC-1,1:nrC),&
                            & xcoord(1:ncC,1:nrC),ycoord(0:ncC-1,1:nrC),&
                            & ycoord(1:ncC,1:nrC),2,3)
   ENDIF
   delglb(0,1:nrC) = delglb(1,1:nrC); delglb(ncC+1,1:nrC) = delglb(ncC,1:nrC)
   lbounds = (/0,1/);  nhdist = (/1,1,0,0/)
   CALL distribute_mod(delglb(0:ncC+1,1:nrC),mgvars(levC)%delxatu,lbounds,&
                     & nhdist,iarr_delxatu,0.0,mglevel=levC)

!
!2.3.2 Y-direction
!------------------
!
!  ---U-nodes
   IF (iopt_grid_htype.LE.2) THEN
      ycoord(1,1:nrC) = 0.5*(gycoordglbC(1,1:nrC)+gycoordglbC(1,2:nrC+1))
      CALL grid_spacings_rectang(delglb(1:ncC+1,1:nrC),gdelyglbC(1:nrC),&
                               & ycoord(1,1:nrC),2,3,'Y')
   ELSE
      CALL grid_spacings_curv(delglb(1:ncC+1,1:nrC),gxcoordglbC(1:ncC+1,1:nrC),&
                            & gxcoordglbC(1:ncC+1,2:nrC+1),&
                            & gycoordglbC(1:ncC+1,1:nrC),&
                            & gycoordglbC(1:ncC+1,2:nrC+1),2,3) 
   ENDIF
   lbounds = 1; nhdist = (/0,1,0,0/)
   CALL distribute_mod(delglb(1:nCC+1,1:nrC),mgvars(levC)%delyatu,lbounds,&
                     & nhdist,iarr_delyatu,0.0,mglevel=levC)

!  ---C-nodes
   IF (iopt_grid_htype.LE.2) THEN
      ycoord(1,0:nrC) = 0.5*(gycoordglbC(1,0:nrC)+gycoordglbC(1,1:nrC+1))
      CALL grid_spacings_rectang(delglb(1:ncC,0:nrC),gdelyglbC(0:nrC),&
                               & ycoord(1,0:nrC),2,3,'Y')
   ELSE
      xcoord(1:ncC,0:nrC+1) = 0.5*(gxcoordglbC(1:ncC,0:nrC+1)+&
                               & gxcoordglbC(2:ncC+1,0:nrC+1))
      ycoord(1:ncC,0:nrC+1) = 0.5*(gycoordglbC(1:ncC,0:nrC+1)+&
                               & gycoordglbC(2:ncC+1,0:nrC+1))
      CALL grid_spacings_curv(delglb(1:ncC,0:nrC),xcoord(1:ncC,0:nrC),&
                            & xcoord(1:ncC,1:nrC+1),ycoord(1:ncC,0:nrC),&
                            & ycoord(1:ncC,1:nrC+1),2,3)
   ENDIF
   lbounds = (/1,0/); nhdist = (/0,0,1,0/)
   CALL distribute_mod(delglb(1:ncC,0:nrC),mgvars(levC)%delyatc,lbounds,nhdist,&
                     & iarr_delyatc,0.0,mglevel=levC)

!  ---V-nodes
   IF (iopt_grid_htype.LE.2) THEN
      CALL grid_spacings_rectang(delglb(1:ncC,1:nrC),gdelyglbC(1:nrC),&
                               & gycoordglbC(1,1:nrC),2,3,'Y')
   ELSE
      xcoord(1:ncC,0:nrC) = 0.25*(gxcoordglbC(1:ncC,0:nrC)+&
                                & gxcoordglbC(2:ncC+1,0:nrC)+&
                                & gxcoordglbC(1:ncC,1:nrC+1)+&
                                & gxcoordglbC(2:ncC+1,1:nrC+1))
      ycoord(1:ncC,0:nrC) = 0.25*(gycoordglbC(1:ncC,0:nrC)+&
                                & gycoordglbC(2:ncC+1,0:nrC)+&
                                & gycoordglbC(1:ncC,1:nrC+1)+&
                                & gycoordglbC(2:ncC+1,1:nrC+1))
      CALL grid_spacings_curv(delglb(1:ncC,1:nrC),xcoord(1:ncC,0:nrC-1),&
                            & xcoord(1:ncC,1:nrC),ycoord(1:ncC,0:nrC-1),&
                            & ycoord(1:ncC,1:nrC),2,3)
   ENDIF
   delglb(1:ncC,0) = delglb(1:ncC,1); delglb(1:ncC,nrC+1) = delglb(1:ncC,nrC)
   lbounds  = (/1,0/); nhdist = (/0,0,1,1/)
   CALL distribute_mod(delglb(1:ncC,0:nrC+1),mgvars(levC)%delyatv,lbounds,&
                     & nhdist,iarr_delyatv,0.0,mglevel=levC)

!
!2.4 Acceleration of gravity
!---------------------------
!
!  ---interior domain
   IF (ABS(gacc_ref-float_fill).GT.float_fill_eps.OR.iopt_grid_sph.EQ.0) THEN
      mgvars(levC)%gaccatc = gacc_mean
      mgvars(levC)%gaccatu = gacc_mean
      mgvars(levC)%gaccatv = gacc_mean
   ELSE
      ycoordloc = 0.25*degtorad*&
                & (gycoordglbC(nc1locC:nc2locC,nr1locC:nr2locC)+&
                &  gycoordglbC(nc1locC+1:nc2locC+1,nr1locC:nr2locC)+&
                &  gycoordglbC(nc1locC:nc2locC,nr1locC+1:nr2locC+1)+&
                &  gycoordglbC(nc1locC+1:nc2locC+1,nr1locC+1:nr2locC+1)) 
      mgvars(levC)%gaccatc(1:nclocC,1:nrlocC) = 9.78032 + &
                & 0.005172*(SIN(ycoordloc))**2 - 0.00006*(SIN(2.0*ycoordloc))**2
      ycoordloc = 0.5*(gycoordglbC(nc1locC:nc2locC,nr1locC:nr2locC)+&
                     & gycoordglbC(nc1locC:nc2locC,nr1locC+1:nr2locC+1))
      mgvars(levC)%gaccatu(1:nclocC,1:nrlocC) = 9.78032 + &
                & 0.005172*(SIN(ycoordloc))**2 - 0.00006*(SIN(2.0*ycoordloc))**2
      ycoordloc = 0.5*(gycoordglbC(nc1locC:nc2locC,nr1locC:nr2locC)+&
                     & gycoordglbC(nc1locC+1:nc2locC+1,nr1locC:nr2locC))
      mgvars(levC)%gaccatv(1:nclocC,1:nrlocC) = 9.78032 + &
                & 0.005172*(SIN(ycoordloc))**2 - 0.00006*(SIN(2.0*ycoordloc))**2
   ENDIF

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1.AND.iopt_grid_sph.EQ.1.AND.&
     & ABS(gacc_ref-float_fill).LE.float_fill_eps) THEN
      lbounds  = 0; nhexch = (/1,0,1,0/)
      CALL exchange_mod(mgvars(levC)%gaccatc,lbounds,nhexch,iarr_gaccatc,&
                      & mglevel=levC)
      lbounds  = (/0,1/); nhexch = (/1,1,0,0/)
      CALL exchange_mod(mgvars(levC)%gaccatu,lbounds,nhexch,iarr_gaccatu,&
                      & mglevel=levC)
      lbounds  = (/1,0/); nhexch = (/0,0,1,1/)
      CALL exchange_mod(mgvars(levC)%gaccatv,lbounds,nhexch,iarr_gaccatv,&
                      & mglevel=levC)
   ENDIF

!
!2.5 Deallocate/allocate arrays
!------------------------------
!

   DEALLOCATE(gxcoordglbF,gycoordglbF,gdelxglbF,gdelyglbF)
   IF (levC.LT.nomglevels-1) THEN
      ALLOCATE (gxcoordglbF(0:ncC+1,0:nrC+1),STAT=errstat)
      CALL error_alloc('gxcoordglbF',2,(/ncC+2,nrC+2/),kndrtype)
      gxcoordglbF = gxcoordglbC
      ALLOCATE (gycoordglbF(0:ncC+1,0:nrC+1),STAT=errstat)
      CALL error_alloc('gycoordglbF',2,(/ncC+2,nrC+2/),kndrtype)
      gycoordglbF = gycoordglbC
      ALLOCATE(gdelxglbF(0:ncC+1),STAT=errstat)
      CALL error_alloc('gdelxglbF',1,(/ncC+2/),kndrtype)
      gdelxglbF = gdelxglbC
      ALLOCATE(gdelyglbF(0:nrC+1),STAT=errstat)
      CALL error_alloc('gdelyglbF',1,(/nrC+2/),kndrtype)
      gdelyglbF = gdelyglbC
   ENDIF
   DEALLOCATE (delglb,gxcoordglbC,gycoordglbC,gdelxglbC,gdelyglbC,xcoord,&
             & ycoord,ycoordloc)

ENDDO levC_2000

CALL log_timer_out()


RETURN

END SUBROUTINE mg_grid_arrays

!========================================================================

SUBROUTINE mg_pointers_weights
!************************************************************************
!
! *mg_pointers_weights* Pointer arrays and weight factors in the multi-grid
!                       procedure
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Initialisation.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_multigrid
!
! Module calls - exchange_mod
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE multigrid
USE physpars
USE structures
USE switches
USE timepars
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flowdir
INTEGER :: i, ic, idir, iglb, ii, iiloc, j, jc, jdir, jglb, jj, jjloc, levC, &
         & levF, nclocC, nclocF, nrlocC, nrlocF
LOGICAL, DIMENSION(3,3) :: land
INTEGER, DIMENSION(2) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
REAL :: sigcrest

procname(pglev+1) = 'mg_pointers_weights'
CALL log_timer_in()

!
!1. Weight factors on fine grid
!------------------------------
!
!1.1 Interior points
!-------------------
!

mgvars(0)%weights = MERGE(1.0,0.0,nodeatc(0:ncloc+1,0:nrloc+1).GT.0)

!
!1.2 Corections for thin dams
!----------------------------
!

IF (iopt_thndam.EQ.1) THEN

!  ---U-nodes
   iiloc_110: DO iiloc=1,numthinuloc
      i = ithinuloc(iiloc); j = jthinuloc(iiloc)
      mgvars(0)%weights(i-1:i,j) = 0.0
   ENDDO iiloc_110

!  ---V-nodes
   jjloc_120: DO jjloc=1,numthinvloc
      i = ithinvloc(jjloc); j = jthinvloc(jjloc)
      mgvars(0)%weights(i,j-1:j) = 0.0
   ENDDO jjloc_120

ENDIF

!
!1.3 Corrections for weirs/barriers
!---------------------------------
!

IF (iopt_weibar.EQ.1) THEN

!  ---U-nodes
   iiloc_131: DO iiloc=1,numwbaruloc
      i = iwbaruloc(iiloc); j = jwbaruloc(iiloc)
      ii = indexwbaru(iiloc)
      idir_1311: DO idir=1,2
         flowdir = MERGE(udvel(i-1,j).GE.0.0,udvel(i+1,j).LT.0.0,idir.EQ.1)
         IF (flowdir) THEN
            ic = MERGE(i-1,i,idir.EQ.1)
            sigcrest = wbarcrestu(ii)/deptotatc(ic,j)
            IF (sigcrest.GE.1.0) mgvars(0)%weights(i-1:i,j) = 0.0
         ENDIF 
      ENDDO idir_1311
   ENDDO iiloc_131

!  ---V-nodes
   jjloc_132: DO jjloc=1,numwbarvloc
      i = iwbarvloc(jjloc); j = jwbarvloc(jjloc)
      jj = indexwbarv(jjloc)
      jdir_1321: DO jdir=1,2
         flowdir = MERGE(vdvel(i,j-1).GT.0.0,vdvel(i,j+1).LT.0.0,idir.EQ.1)
         IF (flowdir) THEN
            jc = MERGE(j-1,j,jdir.EQ.1)
            sigcrest = wbarcrestv(jj)/deptotatc(i,jc)
            IF (sigcrest.GE.1.0) mgvars(0)%weights(i,j-1:j) = 0.0
         ENDIF 
      ENDDO jdir_1321
   ENDDO jjloc_132

ENDIF

!
!1.4 Black/red nodes
!-------------------
!

IF (nt.EQ.0) THEN
   j_140: DO j=1,nrloc
   i_140: DO i=1,ncloc
      iglb = nc1loc+i-1; jglb = nr1loc+j-1
      IF (MOD(iglb+jglb,2).EQ.0) THEN
         mgvars(0)%blacknode(i,j) = .TRUE.
         mgvars(0)%rednode(i,j) = .FALSE.
      ELSE
         mgvars(0)%blacknode(i,j) = .FALSE.
         mgvars(0)%rednode(i,j) = .TRUE.
      ENDIF
   ENDDO i_140
   ENDDO j_140
ENDIF

!
!1.5 Exchange halo sections
!--------------------------
!

IF (iopt_MPI.EQ.1.AND.(iopt_thndam.EQ.1.OR.iopt_weibar.EQ.1)) THEN
   lbounds = 0; nhexch = 1
   CALL exchange_mod(mgvars(0)%weights,lbounds,nhexch,0,mglevel=0)
ENDIF

!
!2. Coarse grids
!---------------
!

IF (iopt_fld.EQ.2.OR.nt.EQ.0) THEN

   levC_2000: DO levC=1,nomglevels-1

      levF  = levC - 1
      nclocF = mgvars(levF)%ncloc; nrlocF = mgvars(levF)%nrloc
      nclocC = mgvars(levC)%ncloc; nrlocC = mgvars(levC)%nrloc

!
!2.1 Weight factors
!------------------
!
!     ---interior points
      mgvars(levC)%weights(1:nclocC,1:nrlocC) = 0.25*&
          & (mgvars(levF)%weights(1:nclocF:2,1:nrlocF:2)+&
          &  mgvars(levF)%weights(2:nclocF+1:2,1:nrlocF:2)+&
          &  mgvars(levF)%weights(1:nclocF:2,2:nrlocF+1:2)+&
          &  mgvars(levF)%weights(2:nclocF+1:2,2:nrlocF+1:2))
      mgvars(levC)%nodeatc(1:nclocC,1:nrlocC) = &
          & MERGE(0,1,mgvars(levC)%weights(1:nclocC,1:nrlocC).EQ.0)
 
!     ---exchange halo sections
      IF (iopt_MPI.EQ.1) THEN
         lbounds = 0; nhexch = 1
         CALL exchange_mod(mgvars(levC)%weights,lbounds,nhexch,0,mglevel=levC)
      ENDIF
   
!     ---small island correction
      j_211: DO j=1,nrlocC
      i_211: DO i=1,nclocC
         IF (mgvars(levC)%weights(i,j).GT.0.AND.&
           & mgvars(levC)%weights(i,j).LE.0.5) THEN
            land = mgvars(levC)%weights(i-1:i+1,j-1:j+1).EQ.0
            IF (.NOT.(ANY(land))) mgvars(levC)%nodeatc(i,j) = 0
         ENDIF
      ENDDO i_211
      ENDDO j_211
      WHERE (mgvars(levC)%nodeatc(1:nclocC,1:nrlocC).EQ.0)
         mgvars(levC)%weights(1:nclocC,1:nrlocC) = 0.0
      END WHERE
 
!     ---mask array
      mgvars(levC)%maskatc_int = mgvars(levC)%nodeatc(1:nclocC,1:nrlocC).GT.0

!     ---exchange halo sections
      IF (iopt_MPI.EQ.1) THEN
         lbounds = 0; nhexch = 1
         CALL exchange_mod(mgvars(levC)%weights,lbounds,nhexch,0,mglevel=levC)
         CALL exchange_mod(mgvars(levC)%nodeatc,lbounds,nhexch,iarr_nodeatc,&
                         & mglevel=levC)
      ENDIF

!     ---black/red nodes
      IF (nt.EQ.0) THEN
         j_112: DO j=1,nrlocC
         i_112: DO i=1,nclocC
            iglb = mgvars(levC)%nc1loc+i-1; jglb = mgvars(levC)%nr1loc+j-1
            IF (MOD(iglb+jglb,2).EQ.0) THEN
               mgvars(levC)%blacknode(i,j) = .TRUE.
               mgvars(levC)%rednode(i,j) = .FALSE.
            ELSE
               mgvars(levC)%blacknode(i,j) = .FALSE.
               mgvars(levC)%rednode(i,j) = .TRUE.
            ENDIF
         ENDDO i_112
         ENDDO j_112
      ENDIF

!
!2.2 Pointers at U-nodes
!-----------------------
!
!     ---interior points
      j_221: DO j=1,nrlocC
      i_221: DO i=1,nclocC
         mgvars(levC)%node2du(i,j) = COUNT(mgvars(levC)%nodeatc(i-1:i,j).GT.0)
      ENDDO i_221
      ENDDO j_221

!     ---open boundaries
      iiloc_222: DO iiloc=1,mgvars(levC)%nobuloc
         i = mgvars(levC)%iobuloc(iiloc); j = mgvars(levC)%jobuloc(iiloc)
         mgvars(levC)%node2du(i,j) = 3
      ENDDO iiloc_222
 
!     ---exchange halo sections
      IF (iopt_MPI.EQ.1) THEN
         lbounds = 1; nhexch = (/0,1,0,1/)
         CALL exchange_mod(mgvars(levC)%node2du,lbounds,nhexch,0,&
                         & corners=.FALSE.,mglevel=levC)
      ENDIF

!
!2.3 Pointers at V-nodes
!-----------------------
!
!     ---interior points
      j_231: DO j=1,nrlocC
      i_231: DO i=1,nclocC
         mgvars(levC)%node2dv(i,j) = COUNT(mgvars(levC)%nodeatc(i,j-1:j).GT.0)
      ENDDO i_231
      ENDDO j_231

!     ---open boundaries
      jjloc_232: DO jjloc=1,mgvars(levC)%nobvloc
         i = mgvars(levC)%iobvloc(jjloc); j = mgvars(levC)%jobvloc(jjloc)
         mgvars(levC)%node2dv(i,j) = 3
      ENDDO jjloc_232
 
!     ---exchange halo sections
      IF (iopt_MPI.EQ.1) THEN
         lbounds = 1; nhexch = (/0,1,0,1/)
         CALL exchange_mod(mgvars(levC)%node2dv,lbounds,nhexch,0,&
                         & corners=.FALSE.,mglevel=levC)
      ENDIF

   ENDDO levC_2000

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE mg_pointers_weights

!========================================================================

SUBROUTINE mg_regular_partition(mglevel)
!************************************************************************
!
! *mg_regular_partition* Define the domain decomposition for a coarse grid
!                        level in the multigrid procedure partitioning into
!                        nprocsx*nprocsy subdomains for
!                        sub-grids in the multigrid procedure
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)Parallel_Initialisation.f90  V2.10
!
! Description - method is based on simple partition
!             - (eventual )idle domains (as defined on the main grid level)
!               are removed 
!
! Reference -
!
! Calling program - mg_domain_decomposition
!
! External calls -
!
!************************************************************************
!
USE iopars
USE multigrid
USE paralpars
USE physpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
INTEGER, INTENT(IN) :: mglevel

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: group, i, iproc, j, lev, ncmod, ncg, nrg, nrmod, nxdiv, nxloc, &
         & nxmod, nydiv, nyloc, nymod
INTEGER, DIMENSION(nprocsx*nprocsy) :: nx1procs, nx2procs, ny1procs, ny2procs


procname(pglev+1) = 'mg_regular_partition'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!

lev = mglevel

group = 2**(nomglevels-lev-1)
ncmod = MOD(mgvars(lev)%nc,group); nrmod = MOD(mgvars(lev)%nr,group)
IF (ncmod.EQ.0) THEN
   ncg = mgvars(lev)%nc/group; ncmod = group
ELSE
   ncg = mgvars(lev)%nc/group + 1
ENDIF
IF (nrmod.EQ.0) THEN
   nrg = mgvars(lev)%nr/group; nrmod = group
ELSE
   nrg = mgvars(lev)%nr/group + 1
ENDIF

nxdiv = ncg/nprocsx*group; nxmod = MOD(ncg,nprocsx)
nydiv = nrg/nprocsy*group; nymod = MOD(nrg,nprocsy)

!
!2. Simple domain decomposition
!------------------------------
!

j_210: DO j=1,nprocsy
i_210: DO i=1,nprocsx
   iproc = (j-1)*nprocsx + i
   IF (i.LE.(nprocsx-nxmod)) THEN
      nx1procs(iproc) = 1 + (i-1)*nxdiv
      nxloc = nxdiv
   ELSE
      nx1procs(iproc) = 1 + (nprocsx-nxmod)*nxdiv + &
                      & (i-(nprocsx-nxmod+1))*(nxdiv+group)
      nxloc = nxdiv + group
   ENDIF
   IF (i.EQ.nprocsx) nxloc = nxloc - group + ncmod
   nx2procs(iproc) = nx1procs(iproc) + nxloc - 1
   IF (j.LE.(nprocsy-nymod)) THEN
      ny1procs(iproc) = 1 + (j-1)*nydiv
      nyloc = nydiv
   ELSE
      ny1procs(iproc) = 1 + (nprocsy-nymod)*nydiv + &
                      & (j-(nprocsy-nymod+1))*(nydiv+group)
      nyloc = nydiv + group
   ENDIF
   IF (j.EQ.nprocsy) nyloc = nyloc - group + nrmod
   ny2procs(iproc) = ny1procs(iproc) + nyloc - 1
ENDDO i_210
ENDDO j_210

!
!3. Store
!--------
!

mgvars(lev)%nc1procs = nx1procs; mgvars(lev)%nc2procs = nx2procs
mgvars(lev)%nr1procs = ny1procs; mgvars(lev)%nr2procs = ny2procs

CALL log_timer_out()


RETURN

END SUBROUTINE mg_regular_partition

!========================================================================

SUBROUTINE mg_water_depths
!************************************************************************
!
! *mg_water_depths* Total water depths on the multi level grids
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Initialisation.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_multigrid, multi_grid_scheme
!
! Module calls - error_alloc, exchange_mod
!
!************************************************************************
!
USE depths
USE gridpars
USE iopars
USE modids
USE multigrid
USE physpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: levC, levF, ncC, ncF, nrC, nrF
INTEGER, DIMENSION(2) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: weightsu, weightsv, wsum 


procname(pglev+1) = 'mg_water_depths'
CALL log_timer_in()

!
!1. Fine grid
!------------
!

mgvars(0)%deptotatc = deptotatc(0:ncloc+1,0:nrloc+1)
mgvars(0)%deptotatu = deptotatu(0:ncloc+1,1:nrloc+1)
mgvars(0)%deptotatv = deptotatv(1:ncloc+1,0:nrloc+1)

!
!2. Coarse grids
!---------------
!

levC_200: DO levC=1,nomglevels-1

!
!2.1 Initialise and allocate
!---------------------------
!

   levF = levC - 1
   ncC = mgvars(levC)%ncloc; nrC = mgvars(levC)%nrloc
   ncF = mgvars(levF)%ncloc; nrF = mgvars(levF)%nrloc

   ALLOCATE (wsum(ncC,nrC),STAT=errstat)
   CALL error_alloc('wsum',2,(/ncC,nrC/),kndrtype)

!
!2.2 C-nodes
!-----------
!
!  ---interior points
   wsum = mgvars(levF)%weights(1:ncF:2,1:nrF:2) + &
        & mgvars(levF)%weights(1:ncF:2,2:nrF+1:2) + &
        & mgvars(levF)%weights(2:ncF+1:2,1:nrF:2) + &
        & mgvars(levF)%weights(2:ncF+1:2,2:nrF+1:2)
   WHERE (mgvars(levC)%maskatc_int)
      mgvars(levC)%deptotatc(1:ncC,1:nrC) = (&
        & mgvars(levF)%weights(1:ncF:2,1:nrF:2)*&
        & mgvars(levF)%deptotatc(1:ncF:2,1:nrF:2) +&
        & mgvars(levF)%weights(1:ncF:2,2:nrF+1:2)*&
        & mgvars(levF)%deptotatc(1:ncF:2,2:nrF+1:2) +&
        & mgvars(levF)%weights(2:ncF+1:2,1:nrF:2)*&
        & mgvars(levF)%deptotatc(2:ncF+1:2,1:nrF:2) +&
        & mgvars(levF)%weights(2:ncF+1:2,2:nrF+1:2)*&
        & mgvars(levF)%deptotatc(2:ncF+1:2,2:nrF+1:2))/wsum
   ELSEWHERE
      mgvars(levC)%deptotatc(1:ncC,1:nrC) = 0.0
   END WHERE

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = 0; nhexch = 1
      CALL exchange_mod(mgvars(levC)%deptotatc,lbounds,nhexch,iarr_deptotatc,&
                      & mglevel=levC)
   ENDIF

!
!2.3 U-nodes
!-----------
!
!  ---allocate
   ALLOCATE (weightsu(ncF,nrF+1),STAT=errstat)
   CALL error_alloc('weightsu',2,(/ncF,nrF+1/),kndrtype)

!  ---interior points
   weightsu = MERGE(1.0,0.0,mgvars(levF)%node2du(1:ncF,:).GT.1)
   wsum = weightsu(1:ncF:2,1:nrF:2) + weightsu(1:ncF:2,2:nrF+1:2)
   WHERE (wsum.GT.0.0)
      mgvars(levC)%deptotatu(1:ncC,1:nrC) = (&
     & weightsu(1:ncF:2,1:nrF:2)*mgvars(levF)%deptotatu(1:ncF:2,1:nrF:2)+&
     & weightsu(1:ncF:2,2:nrF+1:2)*mgvars(levF)%deptotatu(1:ncF:2,2:nrF+1:2))&
     & /wsum
   ELSEWHERE
      mgvars(levC)%deptotatu(1:ncC,1:nrC) = 0.0
   END WHERE

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = (/0,1/); nhexch = (/1,1,0,1/)
      CALL exchange_mod(mgvars(levC)%deptotatu,lbounds,nhexch,iarr_deptotatu,&
                      & mglevel=levC)
   ENDIF

!  ---deallocate
   DEALLOCATE (weightsu)

!
!2.4 V-nodes
!-----------
!
!  ---allocate
   ALLOCATE (weightsv(ncF+1,nrF),STAT=errstat)
   CALL error_alloc('weightsv',2,(/ncF+1,nrF/),kndrtype)

!  ---interior points
   weightsv = MERGE(1.0,0.0,mgvars(levF)%node2dv(:,1:nrF).GT.1)
   wsum = weightsv(1:ncF:2,1:nrF:2) + weightsv(2:ncF+1:2,1:nrF:2)
   WHERE (wsum.GT.0.0)
      mgvars(levC)%deptotatv(1:ncC,1:nrC) = (&
     & weightsv(1:ncF:2,1:nrF:2)*mgvars(levF)%deptotatv(1:ncF:2,1:nrF:2)+&
     & weightsv(2:ncF+1:2,1:nrF:2)*mgvars(levF)%deptotatv(2:ncF+1:2,1:nrF:2))&
     & /wsum
   ELSEWHERE
      mgvars(levC)%deptotatv(1:ncC,1:nrC) = 0.0
   END WHERE

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = (/1,0/); nhexch = (/0,1,1,1/)
      CALL exchange_mod(mgvars(levC)%deptotatv,lbounds,nhexch,iarr_deptotatv,&
                      & mglevel=levC)
   ENDIF

!  ---deallocate
   DEALLOCATE (weightsv)

!
!2.5 Deallocate
!--------------
!

   DEALLOCATE (wsum)

ENDDO levC_200

CALL log_timer_out()


RETURN

END SUBROUTINE mg_water_depths

!========================================================================

SUBROUTINE open_boundary_mg_arrays
!************************************************************************
!
! *open_boundary_mg_arrays* Locations of open boundaries on the coarse grids 
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Initialisation.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_multigrid
!
! Module calls - combine_stats_glb, error_alloc, error_alloc_struc_comp,
!                local_proc
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
USE syspars
USE error_routines, ONLY: error_alloc, error_alloc_struc_comp
USE grid_routines, ONLY: local_proc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, ii, iiloc, iproc, j, jj, jjloc, lev, nobuC, nobulocC_ext, &
         & nobvC, nobvlocC_ext
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: indu, indv, jndu, jndv


procname(pglev+1) = 'open_boundary_mg_arrays'
CALL log_timer_in()

!
!1. Fine grid
!------------
!
!---U-nodes
mgvars(0)%nobuloc = nobuloc
mgvars(0)%nobuloc_ext = nobuloc_ext
ALLOCATE (mgvars(0)%iobuloc(nobuloc_ext),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','iobuloc',1,(/nobuloc_ext/),kndint,0)
mgvars(0)%iobuloc = iobuloc
ALLOCATE (mgvars(0)%jobuloc(nobuloc_ext),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','jobuloc',1,(/nobuloc_ext/),kndint,0)
mgvars(0)%jobuloc = jobuloc

!---V-nodes
mgvars(0)%nobvloc = nobvloc
mgvars(0)%nobvloc_ext = nobvloc_ext
ALLOCATE (mgvars(0)%iobvloc(nobvloc_ext),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','iobvloc',1,(/nobvloc_ext/),kndint,0)
mgvars(0)%iobvloc = iobvloc
ALLOCATE (mgvars(0)%jobvloc(nobvloc_ext),STAT=errstat)
CALL error_alloc_struc_comp('mgvars','jobvloc',1,(/nobvloc_ext/),kndint,0)
mgvars(0)%jobvloc = jobvloc

!
!2. Coarse grids
!---------------
!

lev_200: DO lev=1,nomglevels-1

   nobuC = mgvars(lev)%nobu; nobvC = mgvars(lev)%nobv

!
!2.1 Allocate work space arrays
!------------------------------
!

   ALLOCATE (indu(nobuC),STAT=errstat)
   CALL error_alloc('indu',1,(/nobuC/),kndint)
   ALLOCATE (jndu(nobuC),STAT=errstat)
   CALL error_alloc('jndu',1,(/nobuC/),kndint)
   ALLOCATE (indv(nobvC),STAT=errstat)
   CALL error_alloc('indv',1,(/nobvC/),kndint)
   ALLOCATE (jndv(nobvC),STAT=errstat)
   CALL error_alloc('jndv',1,(/nobvC/),kndint)

!
!2.2 Open boundary locations (local)
!----------------------------------
!

   iproc_220: DO iproc=1,nprocs

!
!2.2.1 U-nodes
!-------------
!
!     ---internal locations
      iiloc = 0; mgvars(lev)%nobuprocs(iproc) = 0
      ii_2211: DO ii=1,nobuC
         i = mgvars(lev)%iobu(ii); j = mgvars(lev)%jobu(ii)
         IF (local_proc(i,j,iproc=iproc,mglevel=lev)) THEN
            iiloc = iiloc + 1
            indu(iiloc) = i-mgvars(lev)%nc1procs(iproc)+1
            jndu(iiloc) = j-mgvars(lev)%nr1procs(iproc)+1
            mgvars(lev)%indexobuprocs(iiloc,iproc) = ii
         ENDIF
      ENDDO ii_2211
      mgvars(lev)%nobuprocs(iproc) = iiloc
      IF (idloc.EQ.idprocs(iproc)) mgvars(lev)%nobuloc = iiloc
      
!     ---add location in eastern halo      
      ii_2212: DO ii=1,nobuC
         i = mgvars(lev)%iobu(ii); j = mgvars(lev)%jobu(ii)
         IF (i.EQ.(mgvars(lev)%nc2procs(iproc)+1).AND.&
           & j.GE.mgvars(lev)%nr1procs(iproc).AND.&
           & j.LE.mgvars(lev)%nr2procs(iproc)) THEN
            iiloc = iiloc + 1
            indu(iiloc) = mgvars(lev)%ncprocs(iproc)+1
            jndu(iiloc) = j-mgvars(lev)%nr1procs(iproc)+1
            mgvars(lev)%indexobuprocs(iiloc,iproc) = ii
         ENDIF
      ENDDO ii_2212

!     ---allocate and store
      IF (idloc.EQ.idprocs(iproc)) THEN
         nobulocC_ext = iiloc
         mgvars(lev)%nobuloc_ext = nobulocC_ext
         ALLOCATE (mgvars(lev)%iobuloc(nobulocC_ext),STAT=errstat)
         CALL error_alloc_struc_comp('mgvars','iobuloc',1,(/nobulocC_ext/),&
                                   & kndint,lev)
         ALLOCATE (mgvars(lev)%jobuloc(nobulocC_ext),STAT=errstat)
         CALL error_alloc_struc_comp('mgvars','jobuloc',1,(/nobulocC_ext/),&
                                   & kndint,lev)
         IF (nobulocC_ext.GT.0) THEN
            mgvars(lev)%iobuloc = indu(1:nobulocC_ext)
            mgvars(lev)%jobuloc = jndu(1:nobulocC_ext)
            mgvars(lev)%indexobu(1:nobulocC_ext) = &
          & mgvars(lev)%indexobuprocs(1:nobulocC_ext,iproc)
         ENDIF
      ENDIF

!
!2.2.2 V-nodes
!-------------
!
!     ---internal locations
      jjloc = 0; mgvars(lev)%nobvprocs(iproc) = 0
      jj_2221: DO jj=1,nobvC
         i = mgvars(lev)%iobv(jj); j = mgvars(lev)%jobv(jj)
         IF (local_proc(i,j,iproc=iproc,mglevel=lev)) THEN
            jjloc = jjloc + 1
            indv(jjloc) = i-mgvars(lev)%nc1procs(iproc)+1
            jndv(jjloc) = j-mgvars(lev)%nr1procs(iproc)+1
            mgvars(lev)%indexobvprocs(jjloc,iproc) = jj
         ENDIF
      ENDDO jj_2221
      mgvars(lev)%nobvprocs(iproc) = jjloc
      IF (idloc.EQ.idprocs(iproc)) mgvars(lev)%nobvloc = jjloc

!     ---add locations in northern halo
      jj_2222: DO jj=1,nobvC
         i = mgvars(lev)%iobv(jj); j = mgvars(lev)%jobv(jj)
         IF (i.GE.mgvars(lev)%nc1procs(iproc).AND.&
           & i.LE.mgvars(lev)%nc2procs(iproc).AND.&
           & j.EQ.(mgvars(lev)%nr2procs(iproc)+1)) THEN
            jjloc = jjloc + 1
            indv(jjloc) = i-mgvars(lev)%nc1procs(iproc)+1
            jndv(jjloc) = mgvars(lev)%nrprocs(iproc)+1
            mgvars(lev)%indexobvprocs(jjloc,iproc) = jj
         ENDIF
      ENDDO jj_2222
      
!     ---allocate and store
      IF (idloc.EQ.idprocs(iproc)) THEN
         nobvlocC_ext = jjloc
         mgvars(lev)%nobvloc_ext = nobvlocC_ext
         ALLOCATE (mgvars(lev)%iobvloc(nobvlocC_ext),STAT=errstat)
         CALL error_alloc_struc_comp('mgvars','iobvloc',1,(/nobvlocC_ext/),&
                                   & kndint,lev)
         ALLOCATE (mgvars(lev)%jobvloc(nobvlocC_ext),STAT=errstat)
         CALL error_alloc_struc_comp('mgvars','jobvloc',1,(/nobvlocC_ext/),&
                                    & kndint,lev)
         IF (nobvlocC_ext.GT.0) THEN
            mgvars(lev)%iobvloc = indv(1:nobvlocC_ext)
            mgvars(lev)%jobvloc = jndv(1:nobvlocC_ext)
            mgvars(lev)%indexobv(1:nobvlocC_ext) = &
          & mgvars(lev)%indexobvprocs(1:nobvlocC_ext,iproc)
         ENDIF
      ENDIF

   ENDDO iproc_220

!
!2.3 Deallocate work arrays
!-------------------------
!

   DEALLOCATE (indu,jndu,indv,jndv)

ENDDO lev_200

CALL log_timer_out()


RETURN

END SUBROUTINE open_boundary_mg_arrays
