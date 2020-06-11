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

MODULE paralpars
!************************************************************************
!
! *paralpars* Parameters for parallel setup
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)paralpars.f90  V2.11
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


!---compiled/not compiled with MPI
LOGICAL :: parallel_set = .FALSE.

!---communicator
INTEGER :: comm_model

!---model ids
INTEGER :: modelid, modelidcoh, modelidpart, modelidwav

!---number of coupled models
INTEGER :: nocoupledmodels

!---process numbers
INTEGER ::  npworld, nprocs, nprocscoh, nprocspart, nprocswav, nprocsx, nprocsy

!---process ids
LOGICAL :: master, mastermod
INTEGER, PARAMETER :: idmaster = 0
INTEGER :: icoordloc, idloc, idlocglb, idmastercoh, idmastermod, idmasterpart,&
         & idmasterwav, iprocglb, iprocloc, jcoordloc

!---tags
INTEGER, PARAMETER :: itagcohtowav = 1001, itagwavtocoh = 1002

!---array definitions
INTEGER, ALLOCATABLE, DIMENSION(:) :: comprocs, idprocs, idprocsglb
INTEGER, ALLOCATABLE, DIMENSION(:) :: icoordprocs, jcoordprocs
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iddomain
TYPE(ExchComms), DIMENSION(MaxHaloComms) :: halocomms

!---MPI definitions
INTEGER :: any_source_MPI = 0, any_tag_MPI = 0, bsend_overhead_MPI = 0, &
         & character_MPI = 0, comm_null_MPI = 0, comm_world_MPI = 0, &
         & complex_MPI = 0, double_precision_MPI = 0, double_precision2_MPI, &
         & float_MPI, float2_MPI, integer_MPI = 0, logical_MPI = 0, &
         & proc_null_MPI = 0, real_MPI = 0, &
         & real2_MPI = 0, undefined_MPI = 0

SAVE

!
! Name           Type    Purpose
!------------------------------------------------------------------------------
!*comprocs*       INTEGER Rank order for collective communications
!*comm_model*     INTEGER MPI communicator for each model componet
!*halocomms*      DERIVED Properties of halo (send/receive) communications
!*icoordloc*      INTEGER X-domain coordinate of local process
!*icoordprocs*    INTEGER Array with values of icoordloc per process
!*iddomain*       INTEGER Process ranks as defined on Cartesian topology
!*idloc*          INTEGER Rank of local process (within the local communicator)
!*idlocglb*       INTEGER Rank of local process (within the global communicator)
!*idmaster*       INTEGER Rank of master process in global communicator
!*idmastercoh*    INTEGER Global rank of the COHERENS master within the global
!                         communicator
!*idmastermod*    INTEGER Local rank of master process within the local
!                         communicator
!*idmasterpart*   INTEGER Global rank of particle model process
!*idmasterwav*    INTEGER Global rank of the wave master within the global
!                         communicator
!*idprocs*        INTEGER Array of process ranks within the local communicator
!*idprocsglb*     INTEGER Array of process ranks within the global communicator
!*iprocglb*       INTEGER Global process index number
!*iprocloc*       INTEGER Local process index number
!*itagcohtowav*   INTEGER Tag used for sending data by COHERENS to the wave
!                         model
!*itagwavtocoh*   INTEGER Tag used for sending data by the wave model to
!                         the hydrodynamic model
!*jcoordloc*      INTEGER Y-domain coordinate of local process
!*jcoordprocs*    INTEGER Array with values of jcoordloc per process
!*master*         LOGICAL .TRUE. for the master (rank = 0) process within the
!                         global ("world") communicator
!*mastermod*      LOGICAL .TRUE. for the master process within the local model
!                         communicator
!*modelid*        INTEGER Model id within local communicator
!*modelidcoh*     INTEGER Global id of the hydrodynamic model
!*modelidpart*    INTEGER Global id of the particle model
!*modelidwave*    INTEGER Global id of the wave model
!*nocoupledmodels*INTEGER Number of coupled models
!*nprocs*         INTEGER Number of active (non-idle) processes used by the
!                         hydrodynamic model
!*nprocspart*     INTEGER Number of processors used by the particle model
!*nprocswav*      INTEGER Number of processors used by the wave model
!*nprocsx*        INTEGER Number of processes in X-direction
!*nprocsy*        INTEGER Number of processes in Y-direction
!*npworld*        INTEGER Total number of available processes in MPI_comm_world
!*parallel_set*   LOGICAL .TRUE. for enabling parallel mode
!
!************************************************************************
!

END MODULE paralpars
