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

MODULE multigrid
!************************************************************************
!
! *multigrid* Derived type definitions and arrays for the multigrid method
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)multigrid.f90  V2.8
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - 
!
!************************************************************************
!
USE datatypes
USE syspars

IMPLICIT NONE


TYPE :: MultiGridVars

!  ---grid parameters and arrays
   INTEGER :: nc, nobu, nobv, nr
   INTEGER :: ncloc, nc1loc, nc2loc, nrloc, nobuloc, nobuloc_ext, nobvloc, &
            & nobvloc_ext, nr1loc, nr2loc
   LOGICAL, POINTER, DIMENSION(:) :: soutobv, westobu
   INTEGER, POINTER, DIMENSION(:) :: iobu, iobv, jobu, jobv
   INTEGER, POINTER, DIMENSION(:) :: iobuloc, iobvloc, jobuloc, jobvloc
   INTEGER, POINTER, DIMENSION(:) :: indexobu, indexobv
   LOGICAL, POINTER, DIMENSION(:,:) :: maskatc_int
   INTEGER, POINTER, DIMENSION(:,:) :: nodeatc, node2du, node2dv
   REAL, POINTER, DIMENSION(:,:) :: gaccatc, gaccatu, gaccatv
   REAL, POINTER, DIMENSION(:,:) :: delxatc, delxatu, delxatv, &
                                  & delyatc, delyatu, delyatv

!  ---water depths
   REAL, POINTER, DIMENSION(:,:) :: deptotatc, deptotatu, deptotatv
 
!  ---open boundary conditions
   INTEGER, POINTER, DIMENSION(:) :: iloczobu, iloczobv, ityp2dobu, ityp2dobv

!  ---domain decomposition
   INTEGER, POINTER, DIMENSION(:) :: ncprocs, nc1procs, nc2procs, &
                                   & nrprocs, nr1procs, nr2procs
   INTEGER, POINTER, DIMENSION(:) :: nobuprocs, nobvprocs
   INTEGER, POINTER, DIMENSION(:,:) :: indexobuprocs, indexobvprocs
   TYPE(ExchComms), POINTER, DIMENSION(:) :: halocomms

!  ---work space
   LOGICAL :: linearization_done
   LOGICAL, POINTER, DIMENSION(:,:) :: blacknode, rednode
   REAL, POINTER, DIMENSION(:,:) :: correction, precon, residual, rhs, &
                                  & solvec, weights
   REAL, POINTER, DIMENSION(:,:,:) :: subdiags

END type MultiGridVars

TYPE (MultiGridVars), DIMENSION(0:MaxMGLevels-1) :: mgvars

SAVE

!
! Name          Type    Purpose
!------------------------------------------------------------------------------

!*blacknode*    LOGICAL .TRUE. at "black" nodes for the Gauss-Seidel solver
!*correction*   REAL    Entries of the correction vector in the multigrid
!                       algorithm
!*deptotatc*    REAL    Total water depth at C-nodes                        [m]
!*deptotatu*    REAL    Total water depth at U-nodes                        [m]
!*deptotatv*    REAL    Total water depth at V-nodes                        [m]
!*delxatc*      REAL    Grid spacing in X-direction at C-nodes              [m]
!*delxatu*      REAL    Grid spacing in X-direction at U-nodes              [m]
!*delxatv*      REAL    Grid spacing in X-direction at V-nodes              [m]
!*delyatc*      REAL    Grid spacing in Y-direction at C-nodes              [m]
!*delyatu*      REAL    Grid spacing in Y-direction at U-nodes              [m]
!*delyatv*      REAL    Grid spacing in Y-direction at V-nodes              [m]
!*gaccatc*      REAL    Acceleration of gravity at C-nodes              [m^2/s]
!*gaccatu*      REAL    Acceleration of gravity at U-nodes              [m^2/s]
!*gaccatv*      REAL    Acceleration of gravity at V-nodes              [m^2/s]
!*halocomms*    DERIVED Properties of halo (send/receive) communications
!*iloczobu*     INTEGER Position of elevation points used for o.b.c. at U-nodes
!                       (0/2)
!                 =  0 => not specified
!                 =  1 => at U-open boundary
!                 =  2 => at external C-node
!*iloczobv*     INTEGER Position of elevation points used for o.b.c. at V-nodes
!                      (0/2)
!*indexobu*     INTEGER Global indices of the local open boundary points
!                       at U-nodes
!*indexobuprocs*INTEGER Global indices per process of the local open boundary
!                       points at U-nodes
!*indexobv*     INTEGER Global indices of the local open boundary points
!                       at V-nodes
!*indexobvprocs*INTEGER Global indices per process of the local open boundary
!                       points at V-nodes
!*iobu*         INTEGER X-indices of the global open boundary points at U-nodes
!*iobuloc*      INTEGER X-indices of the local open boundary points at U-nodes
!                       including points within the first column of the
!                       eastern halo
!*iobv*         INTEGER X-indices of the global open boundary points at V-nodes
!*iobvloc*      INTEGER X-indices of the local open boundary points at V-nodes
!                       including points within the first column of the
!                       northern halo
!*ityp2dobu*    INTEGER Type of 2-D o.b.c. at U-nodes (0/17)
!*ityp2dobv*    INTEGER Type of 2-D o.b.c. at V-nodes (0/17)
!*jobu*         INTEGER Y-indices of the global open boundary points at U-nodes
!*jobuloc*      INTEGER Y-indices of the local open boundary points at U-nodes
!                       including points within the first column of the
!                       eastern halo
!*jobv*         INTEGER Y-indices of the global open boundary points at V-nodes
!*jobvloc*      INTEGER Y-indices of the local open boundary points at V-nodes
!                       including points within the first column of the
!                       northern halo
!*linearization_done* LOGICAL .TRUE. if linearization is already performed
!                     in the present iteration step at the present grid level
!*nc*           INTEGER Global number of grid cells in the X-direction
!*ncloc*        INTEGER Local number of grid cells in the X-direction
!*ncprocs*      INTEGER Array with number of grid points per process in
!                       X-direction
!*nc1loc*       INTEGER Global X-index of first column on the local domain grid
!*nc1procs*     INTEGER Array with values of nc1loc per process
!*nc2loc*       INTEGER Global X-index of last column on the local domain grid
!*nc2procs*     INTEGER Array with values of nc2loc per process
!*nobu*         INTEGER Number of open boundary points at U-nodes
!*nobuloc*      INTEGER Local number of interior open boundary points at U-nodes
!*nobuloc_ext*  INTEGER Local number of interior open boundary points at U-nodes
!                       including points within the first column of the eastern
!                       halo
!*nobuprocs*    INTEGER Array with values of nobuloc per process
!*nobv*         INTEGER Number of open boundary points at V-nodes
!*nobvloc*      INTEGER Local number of interior open boundary points at V-nodes
!*nobvloc_ext*  INTEGER Local number of interior open boundary points at V-nodes
!                       including points within the first row of the northern
!                       halo
!*nobvprocs*    INTEGER Array with values of nobvloc per process
!*nodeatc*      INTEGER Pointers at C-nodes
!                  = 0 => dry
!                  = 1 => wet
!*nodeatu*      INTEGER Pointers at U-nodes
!                  = 0 => dry (land) cell face
!                  = 1 => coastal boundary
!                  = 2 => wet velocity node
!                  = 3 => open sea boundary
!                  = 4 => river open boundary
!*nodeatv*      INTEGER Pointers at V-nodes
!                  = 0 => dry (land) cell face
!                  = 1 => coastal boundary
!                  = 2 => wet velocity node
!                  = 3 => open sea boundary
!                  = 4 => river open boundary
!*nr*           INTEGER Global number of grid cells in the Y-direction
!*nrloc*        INTEGER Local number of grid cells in the Y-direction
!*nrprocs*      INTEGER Array with number of grid points per process in
!                       Y-direction
!*nr1procs*     INTEGER Array with values of ncrloc per process
!*nr1loc*       INTEGER Global Y-index of first row on the local domain grid
!*nr2procs*     INTEGER Array with values of ncrloc per process
!*nr2loc*       INTEGER Global Y-index of last row on the local domain grid
!*precon*       REAL    Diagonal entries of the preconditioner matrix
!*rednode*      LOGICAL .TRUE. at "red" nodes for the Gauss-Seidel solver
!*residual*     REAL    Entries of the residual vector
!*rhs*          REAL    Right hand side of the implicit matrix
!*solvec*       REAL    Solution vector
!*soutobv*      LOGICAL .TRUE. (.FALSE.) at South (North) V-open (global)
!                       boundaries
!*subdiags*     REAL    Sub-diagonal entries of the 5-diagonal implicit matrix
!                (:,:,-2) => sub-diogonal with western entries
!                (:,:,-1) => sub-diogonal with southern entries
!                (:,:,0)  => sub-diogonal with main diagonal entries
!                (:,:,1)  => sub-diogonal with eastern entries
!                (:,:,2)  => sub-diogonal with northern entries
!*westobu*      LOGICAL .TRUE. (.FALSE.) at West (East) V-open (global)
!                       boundaries
!*weights*      REAL   Weight factors for bilinear interpolation in the
!                      multigrid algorithm
!
!******************************************************************************
!

END MODULE multigrid
