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

MODULE grid
!************************************************************************
!
! *grid* Model grid arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid.f90  V2.11
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
!************************************************************************
!
 
IMPLICIT NONE

LOGICAL, ALLOCATABLE, DIMENSION(:) ::  soutobv, soutoby, westobu, westobx
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: maskatc_int, seapoint, seapointglb
INTEGER, ALLOCATABLE, DIMENSION(:) :: iobu, iobv, iobx, ioby,&
                                    & jobu, jobv, jobx, joby
INTEGER, ALLOCATABLE, DIMENSION(:) :: iobuloc, iobvloc, iobxloc, iobyloc, &
                                    & jobuloc, jobvloc, jobxloc, jobyloc
INTEGER, ALLOCATABLE, DIMENSION(:) :: indexobu, indexobv, indexobx, indexoby
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nodeatc, node2du, node2duv, node2dv
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: nodeatu, nodeatuv, nodeatuw, &
                                        & nodeatv, nodeatvw
REAL, ALLOCATABLE, DIMENSION(:) :: gdelxglb, gdelyglb
REAL, ALLOCATABLE, DIMENSION(:,:) :: alphatc_fld, alphatu_fld, alphatv_fld, &
                                   & coriolatu, coriolatv
REAL, ALLOCATABLE, DIMENSION(:,:) :: gaccatc, gaccatu, gaccatv, gangleatc, &
                                   & garea, gxlon, gylat, rmaskatc
REAL, ALLOCATABLE, DIMENSION(:,:) :: gxcoord, gxcoordglb, gxcoordglbatc, &
                                   & gycoord, gycoordglb, gycoordglbatc
REAL, ALLOCATABLE, DIMENSION(:,:) :: rlxobcatu, rlxobcatv
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: gscoordglb, gscoordatc, gscoordatu, &
                                     & gscoordatuvw, gscoordatv, gscoordatw, &
                                     & gscoordatuw, gscoordatvw
REAL, ALLOCATABLE, DIMENSION(:) :: gsigcoordatc, gsigcoordatw
REAL, ALLOCATABLE, DIMENSION(:,:) :: delxatc, delxatu, delxatv, delxatuv
REAL, ALLOCATABLE, DIMENSION(:,:) :: delyatc, delyatu, delyatv, delyatuv
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: delzatc, delzatu, delzatv, delzatuv, &
                                     & delzatuw, delzatvw, delzatw

!---parallel processing
INTEGER, ALLOCATABLE, DIMENSION(:) :: ncprocs, nc1procs, nc2procs
INTEGER, ALLOCATABLE, DIMENSION(:) :: nrprocs, nr1procs, nr2procs
INTEGER, ALLOCATABLE, DIMENSION(:) :: nobuprocs, nosbuprocs, nrvbuprocs
INTEGER, ALLOCATABLE, DIMENSION(:) :: nobvprocs, nosbvprocs, nrvbvprocs
INTEGER, ALLOCATABLE, DIMENSION(:) :: nobxprocs, nobyprocs
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: indexobuprocs, indexobvprocs,&
                                      & indexobxprocs, indexobyprocs

SAVE

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*alphatc_fld*  REAL    Factor multiplying advecive terms terms in turbulence
!                       transport equations for flooding/drying scheme
!*alphatu_fld*  REAL    Factor multiplying terms in U-momentum equations for
!                       flooding/drying scheme
!*alphatv_fld*  REAL    Factor multiplying terms in V-momentum equations for
!                       flooding/drying scheme
!*coriolatu*    REAL    Coriolis frequency at U-nodes                    [rad/s]
!*coriolatv*    REAL    Coriolis frequency at V-nodes                    [rad/s]
!*delxatc*      REAL    Grid spacing in X-direction at C-nodes               [m]
!*delxatu*      REAL    Grid spacing in X-direction at U-nodes               [m]
!*delxatuv*     REAL    Grid spacing in X-direction at UV-nodes              [m]
!*delxatv*      REAL    Grid spacing in X-direction at V-nodes               [m]
!*delyatc*      REAL    Grid spacing in Y-direction at C-nodes               [m]
!*delyatu*      REAL    Grid spacing in Y-direction at U-nodes               [m]
!*delyatuv*     REAL    Grid spacing in Y-direction at UV-nodes              [-]
!*delyatv*      REAL    Grid spacing in Y-direction at V-nodes               [m]
!*delzatc*      REAL    Vertical grid spacing at C-nodes                     [m]
!*delzatu*      REAL    Vertical grid spacing at U-nodes                     [m]
!*delzatuv*     REAL    Vertical grid spacing at UV-nodes                    [m]
!*delzatuw*     REAL    Vertical grid spacing at UW-nodes                    [m]
!*delzatv*      REAL    Vertical grid spacing at V-nodes                     [m]
!*delzatvw*     REAL    Vertical grid spacing at VW-nodes                    [m]
!*delzatw*      REAL    Vertical grid spacing at W-nodes                     [m]
!*gaccatc*      REAL    Acceleration of gravity at C-nodes               [m^2/s]
!*gaccatu*      REAL    Acceleration of gravity at U-nodes               [m^2/s]
!*gaccatv*      REAL    Acceleration of gravity at V-nodes               [m^2/s]
!*gangleatc*    REAL    Angle of X-axis to the reference X-axis            [rad]
!*garea*        REAL    Area of the model grid                             [m^2]
!*gdelxglb*     REAL    Grid spacings in the X-direction of a non-uniform
!                       rectangular grid (global array)
!                                                       [m or degrees longitude]
!*gdelyglb*     REAL    Grid spacings in the Y-direction of a non-uniform
!                       rectangular grid (global array)  [m or degrees latitude]
!*gscoordglb*   REAL    Sigma coordinates at C-nodes
!*gscoordatu*   REAL    Sigma coordinates at U-nodes
!*gscoordatuvw* REAL    Sigma coordinates (local) at UVW-nodes
!*gscoordatuw*  REAL    Sigma coordinates at UW-nodes
!*gscoordatv*   REAL    Sigma coordinates at V-nodes
!*gscoordatvw*  REAL    Sigma coordinates at VW-nodes
!*gscoordatw*   REAL    Sigma coordinates at W-nodes
!*gscoordglb*   REAL    Sigma coordinates (global) at W- or UVW-nodes
!*gsigcoordatc* REAL    Sigma cordinates at C-nodes in case of horizontal
!                       uniform grid
!*gsigcoordatw* REAL    Sigma cordinates at W-nodes in case of horizontal
!                       uniform grid
!*gxcoord*      REAL    X-coordinates (local) of cell corners
!                                                       [m or degrees longitude]
!*gxcoordglb*   REAL    X-coordinates (global) of cell corners
!                                                       [m or degrees longitude]
!*gxcoordglbatc*REAL    X-coordinates (global) at C-nodes
!                                                       [m or degrees longitude]
!*gxlon*        REAL    Longitude at C-nodes                               [rad]
!*gycoord*      REAL    Y-coordinates (local) of cell corners
!                                                        [m or degrees latitude]
!*gycoordglb*   REAL    Y-coordinates (global)of cell corners
!                                                        [m or degrees latitude]
!*gycoordglbatc*REAL    Y-coordinates (global) at C-nodes
!                                                        [m or degrees latitude]
!*gylat*        REAL    Latitude at C-nodes                                [rad]
!*indexobu*     INTEGER Global indices of the local open boundary points
!                       at U-nodes
!*indexobuprocs*INTEGER Global indices per process of the local open boundary
!                       points at U-nodes
!*indexobv*     INTEGER Global indices of the local open boundary points
!                       at V-nodes
!*indexobvprocs*INTEGER Global indices per process of the local open boundary
!                       points at V-nodes
!*indexobx*     INTEGER Global indices of the local open boundary points
!                       at X-nodes
!*indexobxprocs*INTEGER Global indices per process of the local open boundary
!                       points at X-nodes
!*indexoby*     INTEGER Global indices of the local open boundary points
!                       at Y-nodes
!*indexobyprocs*INTEGER Global indices per process of the local open boundary
!                       points at Y-nodes
!*iobu*         INTEGER X-indices of the global open boundary points at U-nodes
!*iobuloc*      INTEGER X-indices of the local open boundary points at U-nodes
!                       including points within the first column of the
!                       eastern halo
!*iobv*         INTEGER X-indices of the global open boundary points at V-nodes
!*iobvloc*      INTEGER X-indices of the local open boundary points at V-nodes
!                       including points within the first column of the
!                       northern halo
!*iobx*         INTEGER X-indices of the global open boundary points at X-nodes
!*iobxloc*      INTEGER X-indices of the local open boundary points at X-nodes
!                       including points within the first column of the
!                       eastern halo
!*ioby*         INTEGER X-indices of the global open boundary points at Y-nodes
!*iobyloc*      INTEGER X-indices of the local open boundary points at Y-nodes
!                       including points within the first column of the
!                       northern halo
!*jobu*         INTEGER Y-indices of the global open boundary points at U-nodes
!*jobuloc*      INTEGER Y-indices of the local open boundary points at U-nodes
!                       including points within the first column of the
!                       eastern halo
!*jobv*         INTEGER Y-indices of the global open boundary points at V-nodes
!*jobvloc*      INTEGER Y-indices of the local open boundary points at V-nodes
!                       including points within the first column of the
!                       northern halo
!*jobx*         INTEGER Y-indices of the global open boundary points at X-nodes
!*jobxloc*      INTEGER Y-indices of the local open boundary points at X-nodes
!                       including points within the first column of the
!                       eastern halo
!*joby*         INTEGER Y-indices of the global open boundary points at Y-nodes
!*jobyloc*      INTEGER Y-indices of the local open boundary points at Y-nodes
!                       including points within the first column of the
!                       northern halo
!*maskatc_int*  LOGICAL .TRUE. (.FALSE.) at wet (dry) interior points
!*ncprocs*      INTEGER Array with number of grid points per process in
!                       X-direction
!*nc1procs*     INTEGER Array with values of nc1loc per process
!*nc2procs*     INTEGER Array with values of nc2loc per process
!*nobuprocs*    INTEGER Array with values of nobuloc per process
!*nobvprocs*    INTEGER Array with values of nobvloc per process
!*nobxprocs*    INTEGER Array with values of nobxloc per process
!*nobyprocs*    INTEGER Array with values of nobyloc per process
!*nodeatc*      INTEGER Pointers at C-nodes
!                  = 0 => dry
!                  = 1 => wet
!*nodeatu*      INTEGER Pointers at U-nodes
!                  = 0 => dry (land) cell face
!                  = 1 => coastal boundary
!                  = 2 => wet velocity node
!                  = 3 => open sea boundary
!                  = 4 => river open boundary
!*nodeatuv*     INTEGER Pointers at corner nodes
!                  = 0 => at least two surrrounding U-nodes or at least two
!                         surrrounding V-nodes are dry
!                  = 1 => at most one surrounding U-node and at most one
!                         surrounding V-node is dry and none of the four
!                         surrounding velocity nodes are open boundaries
!                  = 2 => X-node open boundary (at least one of the surrounding
!                         U-nodes is an open boundary), but not Y-node open
!                         boundary 
!                  = 3 => Y-node open boundary (at least one of the surrounding
!                         V-nodes is an open boundary), but not X-node open
!                         boundary
!                  = 4 => X- and Y-node boundary
!*nodeatuw*     INTEGER Pointers at UW-nodes
!                  = 0 => dry (land) cell face
!                  = 1 => coastal boundary
!                  = 2 => wet velocity node
!                  = 3 => open sea boundary
!                  = 4 => river open boundary
!*nodeatv*      INTEGER Pointers at V-nodes
!                  = 0 => dry (land) cell face
!                  = 1 => coastal boundary
!                  = 2 => wet velocity node
!                  = 3 => bottom boundary
!                  = 4 => open sea boundary
!                  = 5 => river open boundary
!*nodeatvw*     INTEGER Pointers at VW-nodes
!                  = 0 => dry (land) cell face
!                  = 1 => coastal boundary
!                  = 2 => wet velocity node
!                  = 3 => bottom boundary
!                  = 4 => open sea boundary
!                  = 5 => river open boundary
!*nodeisblack*  LOGICAL .TRUE. at black nodes (for Gauss-Seidel RB solver)
!*nodeisred*    LOGICAL .TRUE. at red nodes (for Gauss-Seidel RB solver)
!*node2du*      INTEGER Pointers at U-nodes (2-D case)
!                  = 0 => dry (land) cell face
!                  = 1 => all points in the vertical are coastal boundaries
!                  = 2 => at leat one point in the vertical is a wet velocity
!                         node
!                  = 3 => open sea boundary
!                  = 4 => river open boundary
!*node2duv*     INTEGER Pointers at corner nodes (2-D case)
!                  = 0 => at least two surrrounding U-nodes or at least two
!                         surrrounding V-nodes are dry (for all vertical nodes)
!                  = 1 => at most one surrounding U-node and at most one
!                         surrounding V-node is dry and none of the four
!                         surrounding velocity nodes are open boundaries
!                         (for at least one vertical node)
!                  = 2 => X-node open boundary (at least one of the surrounding
!                         U-nodes is an open boundary), but not Y-node open
!                         boundary 
!                  = 3 => Y-node open boundary (at least one of the surrounding
!                         V-nodes is an open boundary), but not X-node open
!                         boundary
!                  = 4 => X- and Y-node boundary (for at least one vertical
!                         node)
!*node2dv*      INTEGER Pointers at V-nodes (2-D case)
!                  = 0 => dry (land) cell face
!                  = 1 => all points in the vertical are coastal boundaries
!                  = 2 => at leat one point in the vertical is a wet velocity
!                         node
!                  = 3 => open sea boundary
!                  = 4 => river open boundary
!*nosbuprocs*   INTEGER Array with values of nosbuloc per process
!*nosbvprocs*   INTEGER Array with values of nosbvloc per process
!*nrprocs*      INTEGER Array with number of grid points per process in
!                       Y-direction
!*nrvbuprocs*   INTEGER Array with values of nrvbuloc per process
!*nrvbvprocs*   INTEGER Array with values of nrvbvloc per process
!*nr1procs*     INTEGER Array with values of ncrloc per process
!*nr2procs*     INTEGER Array with values of ncrloc per process
!*rmaskatc*     REAL    Equals 1.0 at wet and 0.0 at dry cells
!*rlxobcatu*    REAL    Relaxation factor for momentum advection at U-nodes
!*rlxobcatv*    REAL    Relaxation factor for momentum advection at V-nodes
!*seapoint*     LOGICAL Local mask array set to .FALSE. for C-points on
!                       permanent land
!*seapointglb*  LOGICAL Global mask array set to .FALSE. for C-points on
!                       permanent land
!*soutobv*      LOGICL  .TRUE. (.FALSE.) at South (North) V-open (global)
!                       boundaries
!*soutoby*      LOGICAL .TRUE. (.FALSE.) at South (North) corner-open (global)
!                       boundaries
!*westobu*      LOGICAL .TRUE. (.FALSE.) at West (East) V-open (global)
!                       boundaries
!*westobx*      LOGICAL .TRUE. (.FALSE.) at West (East) corner-open (global)
!                       boundaries
!
!************************************************************************
!

END MODULE grid
