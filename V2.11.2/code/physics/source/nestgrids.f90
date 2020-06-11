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

MODULE nestgrids
!************************************************************************
!
! *nestgrids* Parameters and arrays for nesting
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)nestgrids.f90  V2.7
!
! $Date: 2015-04-10 11:30:34 +0200 (Fri, 10 Apr 2015) $
!
! $Revision: 841 $
!
! Description - 
!
!************************************************************************
!
USE datatypes

IMPLICIT NONE


INTEGER :: nonestsets = 0
INTEGER, ALLOCATABLE, DIMENSION(:) :: nestcoords
INTEGER, ALLOCATABLE, DIMENSION(:) :: nohnstglbc, nohnstglbu, nohnstglbv, &
                                    & nohnstglbx, nohnstglby
INTEGER, ALLOCATABLE, DIMENSION(:) :: nohnstatc, nohnstatu, nohnstatv, &
                                    & nohnstatx, nohnstaty
INTEGER, ALLOCATABLE, DIMENSION(:) :: novnst
INTEGER, ALLOCATABLE, DIMENSION(:) :: lbhnstatc, lbhnstatu, lbhnstatv, &
                                    & lbhnstatx, lbhnstaty
INTEGER, ALLOCATABLE, DIMENSION(:) :: inst2dtype
INTEGER, ALLOCATABLE, DIMENSION(:) :: nobionst, nosednst
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: instbio, instsed
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nohnstcprocs, nohnstuvprocs, &
                                      & nohnstxyprocs
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: indexnstc, indexnstuv, indexnstxy

TYPE (HRelativeCoords), ALLOCATABLE, DIMENSION(:) :: &
     & hnstctoc, hnstctou, hnstctov, hnstutou, hnstvtov, hnstvtox, hnstutoy
TYPE (VRelativeCoords), ALLOCATABLE, DIMENSION(:,:,:,:) :: &
     & vnstctoc, vnstutou, vnstvtov, vnstvtox, vnstutoy

SAVE

!
! Name           Type    Purpose
!------------------------------------------------------------------------------
!*hnstctoc*      DERIVED Relative horizontal coordinates of C-node nest
!                        locations with respect to local model C-nodes
!*hnstctou*      DERIVED Relative horizontal coordinates of U-node nest
!                        locations for normal open boundary conditions with
!                        respect to local model C-nodes
!*hnstctov*      DERIVED Relative horizontal coordinates of V-node nest
!                        locations for normal open boundary conditions with
!                        respect to local model C-nodes
!*hnstutou*      DERIVED Relative horizontal coordinates of U-node nest
!                        locations for normal open boundary conditions with
!                        respect to local model U-nodes
!*hnstvtov*      DERIVED Relative horizontal coordinates of V-node nest
!                        locations for normal open boundary conditions with
!                        respect to local model V-nodes
!*hnstvtox*      DERIVED Relative horizontal coordinates of V-node nest
!                        locations for tangential open boundary conditions with
!                        respect to local model X-nodes
!*hnstutoy*      DERIVED Relative horizontal coordinates of U-node nest
!                        locations for tangential open boundary conditions with
!                        respect to local model Y-nodes
!*indexnstc*     INTEGER Local to global index map (C-node)
!*indexnstuv*    INTEGER Local to global index map (U- and V-nodes)
!*indexnstxy*    INTEGER Local to global index map (X- and Y-nodes)
!*instbio*       INTEGER Variable indices of the biological state variables
!                        used for nesting per set
!*instsed*       INTEGER Variable indices of the sediment fractions used for
!                        nesting per set
!*inst2dtype*    INTEGER Type of data in 2-D nest files for (U,V)-nodes
!                      = 1 => transports and elevations
!                      = 2 => elevations
!                      = 3 => transports
!*lbhnstatc*     INTEGER Buffer index of a C-node nest within a full array
!                        containing all nest points on the local domain
!*lbhnstatu*     INTEGER Buffer index of a U-node nest within a full array
!                        containing all nest points on the local domain
!*lbhnstatv*     INTEGER Buffer index of a V-node nest within a full array
!                        containing all nest points on the local domain
!*lbhnstatx*     INTEGER Buffer index of a X-node nest within a full array
!                        containing all nest points on the local domain
!*lbhnstaty*     INTEGER Buffer index of a Y-node nest within a full array
!                        containing all nest points on the local domain
!*nestcoords*    INTEGER Type of coordinates for grid locations (1/2)
!*nobionst*      INTEGER Number of nested biological state variables per nested
!                        sub-grid
!*nohnstatc*     INTEGER Local number of horizontal C-node nest points per set
!*nohnstatu*     INTEGER Local number of horizontal U-node nest points per set
!*nohnstatv*     INTEGER Local number of horizontal V-node nest points per set
!*nohnstatx*     INTEGER Local number of horizontal X-node nest points per set
!*nohnstaty*     INTEGER Local number of horizontal Y-node nest points per set
!*nohnstcprocs*  INTEGER Number of C-node points per process and set
!*nohnstglbc*    INTEGER Global number of horizontal C-node nest points per set
!*nohnstglbu*    INTEGER Global number of horizontal U-node nest points per set
!*nohnstglbv*    INTEGER Global number of horizontal V-node nest points per set
!*nohnstglbx*    INTEGER Global number of horizontal X-node nest points per set
!*nohnstglby*    INTEGER Global number of horizontal Y-node nest points per set
!*nohnstuvprocs* INTEGER Number of U- and V-node points per process and set
!*nohnstxyprocs* INTEGER Number of X- and Y-node points per process and set
!*nonestsets*    INTEGER Total number of nested areas
!*nosednst*      INTEGER Number of sediment fractions per nested sub-grid
!*novnst*        INTEGER Number of vertical nest points per set
!*vnstctoc*      DERIVED Relative vertical coordinates of C-node locations
!                        with respect to local model C-nodes
!*vnstutou*      DERIVED Relative vertical coordinates of U-node locations for
!                        normal open boundary conditions with respect to local
!                        model U-nodes
!*vnstutoy*      DERIVED Relative vertical coordinates of U-node locations for
!                        tangential open boundary conditions with respect to
!                        local model Y-nodes
!*vnstvtov*      DERIVED Relative vertical coordinates of V-node locations for
!                        normal open boundary conditions with respect to local
!                        model V-nodes
!*vnstvtox*      DERIVED Relative vertical coordinates of V-node locations for
!                        tangential open boundary conditions with respect to
!                        local model X-nodes
!
!------------------------------------------------------------------------------
!

END MODULE nestgrids
