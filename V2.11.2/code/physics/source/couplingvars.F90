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

MODULE couplingvars
!************************************************************************
!
! *couplingvars* MCT types for coupling COHERENS and wave model
!
! Author - Alexander Breugem
!
! Version - @(COHERENS)couplingvars.f90  V2.9
!
! $Date: 2018-02-19 10:06:10 +0100 (Mon, 19 Feb 2018) $
!
! $Revision: 1092 $
!
! Description - 
!
!************************************************************************
!
!---MCT datatypes
#ifdef MCT
   USE m_AttrVect, ONLY: AttrVect
   USE m_GlobalSegMap, ONLY: GlobalSegMap
   USE m_Router, ONLY: Router
   USE m_SparseMatrix, ONLY : SparseMatrix
   USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus
#endif /*MCT*/

IMPLICIT NONE

!---MCT type definitions
#ifdef MCT
   TYPE(AttrVect) :: Av_CohRecvCoh, Av_CohRecvWav, Av_CohSendWav, &
                   & Av_WavRecvCoh, Av_WavRecvWav, Av_WavSendCoh
   TYPE(GlobalSegMap) :: Gsm_CohRecvCoh, Gsm_CohRecvWav, Gsm_CohSendWav, &
                       & Gsm_WavRecvCoh, Gsm_WavRecvWav, Gsm_WavSendCoh
   TYPE(Router) :: Rout_CohtoWav, Rout_WavtoCoh
   TYPE(SparseMatrix) :: Sm_CohinttoWav, Sm_WavinttoCoh
   TYPE(SparseMatrixPlus) :: Smp_CohinttoWav, Smp_WavinttoCoh
#endif /*MCT*/

SAVE

!
! Name           Type     Purpose
!------------------------------------------------------------------------------
!*Av_CohRecvCoh*   DERIVED Attribute vector for the wave data interpolated on
!                          the COHERENS grid
!*Av_CohRecvWav*   DERIVED Attribute vector for the data received by COHERENS
!                          (on the wave grid) from the wave model
!*Av_CohSendWav*   DERIVED Attribute vector for the COHERENS data send to the
!                          wave model by COHERENS
!*Av_WavRecvCoh*   DERIVED Attribute vector for the data received by the wave
!                          model (on the COHERENS grid) from COHERENS
!*Av_WavRecvWav*   DERIVED Attribute vector for the COHERENS data interpolated!
!                          on the wave grid
!*Av_WavSendCoh*   DERIVED Attribute vector for the wave data (on the wave
!                          grid) send to COHERENS by the wave model
!
!*Gsm_CohRecvCoh*  DERIVED Global segmentation map for interpolating wave data
!                          on the COHERENS grid
!*Gsm_CohRecvWav*  DERIVED Global segmentation map for receiving wave data by
!                          COHERENS on the wave model grid
!*Gsm_CohSendWav*  DERIVED Global segmentation map for sending COHERENS data to
!                          the wave model on the COHERENS grid
!*Gsm_WavRecvCoh*  DERIVED Global segmentation map for receiving COHERENS data
!                          by the wave model on the COHERENS grid 
!*Gsm_WavRecvWav*  DERIVED Global segmentation map for interpolating COHERENS
!                          data on the wave model grid
!*Gsm_WavSendCoh*  DERIVED Global segmentation map for sending wave data to
!                          COHERENS on the wave model grid
!
!*Rout_CohtoWav*   DERIVED Router for sending data from COHERENS to the wave
!                          model
!*Rout_WavtoCoh*   DERIVED Router for sending data from the wave model to
!                          COHERENS
!
!Sm_CohinttoWav*   DERIVED Sparse matrix for interpolation of COHERENS data by
!                          the wave model on the wave model grid
!Sm_WavinttoCoh*   DERIVED Sparse matrix for interpolation of wave data by
!                          COHERENS on the COHEREND grid
!
!*Smp_CohinttoWav* DERIVED SparseMatrixPlus definitions for interpolation of
!                          COHERENS data by the wave model on the wave model
!                          grid
!*Smp_Wavinttocoh* DERIVED SparseMatrixPlus definitions for interpolation of
!                          wave data by COHERENS on the COHERENS grid
!
!************************************************************************
!

END MODULE couplingvars

