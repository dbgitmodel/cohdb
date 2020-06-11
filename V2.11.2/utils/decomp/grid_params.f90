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

MODULE grid_params
!************************************************************************
!
! *grid_params* parameters for grid generators
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_params.f90  V2.9
!
! $Date: 2015-10-06 17:42:25 +0200 (Tue, 06 Oct 2015) $
!
! $Revision: 886 $
!
! Description - 
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

CHARACTER (LEN=leniofile) :: decomp_file, info_file, post_file
INTEGER :: iopt_partit_reset = 1, iopt_partit_post = 0, npwork = 0
INTEGER, ALLOCATABLE, DIMENSION(:) :: nowetcprocs

SAVE

!
! Name                Type    Purpose
!------------------------------------------------------------------------------
!*decomp_file*       CHAR    Domain grid file
!*info_file*         CHAR    Output file with info about partitioning
!*post_file*         CHAR    Output file for postprocessing
!*iopt_partit_reset* INTEGER Disables/enables removal of empty domains
!*iopt_partit_post*  INTEGER Disables/enables writing of mumap files
!*npwork*            INTEGER Initial (trial) number of processes
!*icoord*            INTEGER X-index coordinates of domain grid
!*jcoord*            INTEGER Y-index coordinates of domain grid
!*nowetcprocs*       INTEGER Array with number of local active C-nodes per
!                            process
!
!************************************************************************
!

END MODULE grid_params
