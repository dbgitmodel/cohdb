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

MODULE gridpars
!************************************************************************
!
! *gridpars* Model grid parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)gridpars.f90  V2.10.1
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

!---grid dimensions
INTEGER :: nc = 0, nr = 0
INTEGER :: nz = 0
INTEGER :: ncloc, nc1loc, nc2loc, nrloc, nr1loc, nr2loc

!---number of open boundary points
INTEGER :: nobu, nobv
INTEGER :: nobx, noby
INTEGER :: nosbu = 0, nosbv = 0, nrvbu = 0, nrvbv = 0
INTEGER :: nobuloc, nobvloc
INTEGER :: nobxloc, nobyloc
INTEGER :: nobuloc_ext, nobvloc_ext, nobxloc_ext, nobyloc_ext
INTEGER :: nosbuloc, nosbvloc, nrvbuloc, nrvbvloc

!---shortcuts for regular grid spacings
LOGICAL :: dXregX, dXregY, dYregX, dYregY, dXYreg, dZregZ, rotate_gvecs

!---halo dimensions
INTEGER, PARAMETER :: nhalo = 2
INTEGER :: nhdens, nhfvel, nhscal, nhturb, nh2vel, nh3vel

!---number of sea/active grid pounts
INTEGER :: noseaatc, noseaatcloc, noseaatu, noseaatuloc, noseaatv, noseaatvloc,&
         & nowetatc, nowetatcloc, nowetatu, nowetatuloc, nowetatv, nowetatvloc

!---number of sediment/biological variables used at open boundaries and for
!   nesting
INTEGER :: maxbiovars = 0, maxsedvars = 0

!---mean, maximum and minimum grid spacing
REAL :: delhatc_mean, delxatc_max, delxatc_min, delyatc_max, delyatc_min

!---grid area
REAL :: gareatot

SAVE

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*dehatc_mean* REAL    Mean grid resolution                                 [m]
!*delxatc_max  REAL    Maximum grid spacing in the X-direction              [m]
!*delxatc_min  REAL    Minimum grid spacing in the X-direction              [m]
!*delyatc_max  REAL    Maximum grid spacing in the Y-direction              [m]
!*delyatc_min  REAL    Minimum grid spacing in the YX-direction             [m]
!*dXregX*      LOGICAL .TRUE. if grid spacing in X-direction is uniform in X
!*dXregY*      LOGICAL .TRUE. if grid spacing in X-direction is uniform in Y
!*dXYreg*      LOGICAL .TRUE. if all grid spacings are uniform in X and Y
!*dYregX*      LOGICAL .TRUE. if grid spacing in Y-direction is uniform in X
!*dYregY*      LOGICAL .TRUE. if grid spacing in Y-direction is uniform in Y
!*dZregX*      LOGICAL .TRUE. for uniform vertical sigma-spacings
!*gareatot*    REAL    Total surface area of the model grid
!*maxbiovars*  INTEGER Maximum number of biological variables used at open
!                      boundaries and for nesting
!*maxsedvars*  INTEGER Maximum number of sediment fractions used at open
!                      boundaries and for nesting
!*nc*          INTEGER Number of grid cells in the X-direction
!*ncloc*       INTEGER Local number of grid cells in the X-direction
!*nc1loc*      INTEGER Global X-index of first column on the local domain grid
!*nc2loc*      INTEGER Global X-index of last column on the local domain grid
!*nhalo*       INTEGER Number of halo rows and columns for arrray dimensioning
!*nhdens*      INTEGER Halo size of density arrays                        (1/2)
!*nhfvel*      INTEGER Halo size of horizontal advective currents         (1/2)
!*nhscal*      INTEGER Halo size of scalar (non-density) arrays           (1/2)
!*nhturb*      INTEGER Halo size of turbulence arrays                   (0/1/2)
!*nh2vel*      INTEGER Halo size of 2-D velocity vectors                  (1/2)
!*nh3vel*      INTEGER Halo size of 3-D velocity vectors                  (1/2)
!*nobu*        INTEGER Number of open boundary points at U-nodes
!*nobuloc*     INTEGER Local number of interior open boundary points at U-nodes
!*nobuloc_ext* INTEGER Local number of interior open boundary points at U-nodes
!                      including points within the first column of the eastern
!                      halo
!*nobv*        INTEGER Number of open boundary points at V-nodes
!*nobvloc*     INTEGER Local number of interior open boundary points at V-nodes
!*nobvloc_ext* INTEGER Local number of interior open boundary points at V-nodes
!                      including points within the first row of the northern
!                      halo
!*nobx*        INTEGER Number of open boundary points at X-nodes
!*nobxloc*     INTEGER Local number of interior open boundary points at X-nodes
!*nobxloc_ext* INTEGER Local number of interior open boundary points at X-nodes
!                      including points within the first column of the eastern
!                      halo
!*noby*        INTEGER Number of open boundary points at Y-nodes
!*nobyloc*     INTEGER Local number of interior open boundary points at y-nodes
!*nobyloc_ext* INTEGER Local number of interior open boundary points at Y-nodes
!                      including points within the first column of the northern
!                      halo
!*nosbu*       INTEGER Global number of open sea U-open boundary points
!*nosbuloc*    INTEGER Local number of open sea U-open boundary points inside
!                      domain
!*nosbv*       INTEGER Global number of open sea V-open boundary points
!*nosbvloc*    INTEGER Local number of open sea V-open boundary points inside
!                      domain
!*noseaatc*    INTEGER Number of sea (wet or dry) points on global domain
!*noseaatcloc* INTEGER Number of sea (wet or dry) points on local domain
!*noseaatu*    INTEGER Number of sea (wet or dry) U-points on global domain
!*noseaatuloc* INTEGER Number of sea (wet or dry) U-points on local domain
!*noseaatv*    INTEGER Number of sea (wet or dry) V-points on global domain
!*noseaatvloc* INTEGER Number of sea (wet or dry) V-points on local domain
!*nowetatc*    INTEGER Number of active C-nodes on global domain
!*nowetcloc*   INTEGER Number of active C-nodes on local domain
!*nowetatu*    INTEGER Number of active U-nodes on global domain
!*nowetatuloc* INTEGER Number of active U-nodes on local domain
!*nowetatv*    INTEGER Number of active V-nodes on global domain
!*nowetatvloc* INTEGER Number of active V-nodes on local domain
!*nr*          INTEGER Number of grid cells in the Y-direction
!*nrloc*        INTEGER Local number of grid cells in the Y-direction
!*nrvbu*       INTEGER Global number of river U-open boundary points
!*nrvbuloc*    INTEGER Local number of river U-open boundary points inside
!                      domain
!*nrvbv*       INTEGER Global number of river V-open boundary points
!*nrvbvloc*    INTEGER Local number of river V-open boundary points inside
!                      domain
!*nr1loc*      INTEGER Global Y-index of first row on the local domain grid
!*nr2loc*      INTEGER Global Y-index of last row on the local domain grid
!*nz*          INTEGER Number of grid cells in the vertical direction
!*rotate_gvecs*LOGICAL .TRUE. if vectors in the model grid are rotated with
!                      respect to the reference coordinate system
!
!************************************************************************
!

END MODULE gridpars
