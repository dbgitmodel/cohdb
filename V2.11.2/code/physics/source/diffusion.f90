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

MODULE diffusion
!************************************************************************
!
! *diffusion* Diffusion coefficients
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)diffusion.f90  V2.11.2
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

REAL, ALLOCATABLE, DIMENSION(:,:) :: hdifcoef2datc,  hdifcoef2datuv
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: hdifcoef3datc,  hdifcoef3datu,&
                                     & hdifcoef3datv, hdifcoef3datuv, &
                                     & hdifcoef3datw
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: kinvisc, vdifcoefmom, vdifcoefscal, &
                                     & vdifcoefscal_rot, vdifcoeftke

!---slopes for isolevel or isoneutral diffusion
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: xslopeatu_geo, xslopeatw_geo, &
                                     & yslopeatv_geo, yslopeatw_geo 
REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: xslopeatu_siso, yslopeatv_siso
REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: xslopeatu_ziso, yslopeatv_ziso

SAVE

!
! Name              Type Purpose
!----------------------------------------------------------------------------
!*hdifcoef2datc*    REAL Horizontal 2-D diffusion coefficient at C-nodes [m^2/s]
!*hdifcoef2datuv*   REAL Horizontal 2-D diffusion coefficient at UV-nodes
!                                                                        [m^2/s]
!*hdifcoef3datc*    REAL Horizontal 3-D diffusion coefficient at C-nodes [m^2/s]
!*hdifcoef3datu*    REAL Horizontal 3-D diffusion coefficient at U-nodes [m^2/s]
!*hdifcoef3datuv*   REAL Horizontal 3-D diffusion coefficient at UV-nodes
!                                                                        [m^2/s]
!*hdifcoef3datv*    REAL Horizontal 3-D diffusion coefficient at V-nodes [m^2/s]
!*hdifcoef3datw*    REAL Horizontal 3-D diffusion coefficient at W-nodes [m^2/s]
!*kinvisc*          REAL Kinematic viscosity  at W-nodes                 [m^2/s]
!*vdifcoefmom*      REAL Vertical diffusion coefficient for momentum at W-nodes
!                                                                        [m^2/s]
!*vdifcoefscal*     REAL Vertical diffusion coefficient for scalars at W-nodes
!                                                                        [m^2/s]
!*vdifcoefscal_rot* REAL Isoneutral ("rotated" contribution to the vertical
!                        diffusion coefficient for scalars at W-nodes    [m^2/s]
!*vdifcoeftke*      REAL Vertical diffusion coefficient for t.k.e. at W-nodes
!                                                                        [m^2/s]
!*xslopeatu_geo*    REAL X-xomponent of the isolevel slope at U-nodes
!*xslopeatu_siso*   REAL X-xomponent of the isoneutral slope with respect
!                        to s-coordinate surfaces at U-nodes
!*xslopeatu_ziso*   REAL X-xomponent of the isoneutral slope with respect
!                        to z-coordinate surfaces at U-nodes
!*xslopeatw_geo*    REAL X-xomponent of the isolevel slope at W-nodes
!*yslopeatv_geo*    REAL Y-xomponent of the isolevel slope at V-nodes
!*yslopeatv_siso*   REAL Y-xomponent of the isoneutral slope with respect
!                        to s-coordinate surfaces at V-nodes
!*yslopeatu_ziso*   REAL X-xomponent of the isoneutral slope with respect to
!                        z-coordinate surfaces at V-nodes
!*yslopeatw_geo*    REAL Y-xomponent of the isolevel slope at W-nodes
!
!************************************************************************
!

END MODULE diffusion
