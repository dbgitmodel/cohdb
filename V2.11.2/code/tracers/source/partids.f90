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

MODULE partids
!************************************************************************
!
! *partids* Key ids for the particle model
!
! Author - Valerie Duliere
!
! Version - @(COHERENS)partids.f90  V2.11
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description -
!
! Reference -
!
!************************************************************************
!

USE syspars

IMPLICIT NONE


!---start key id number
INTEGER, PARAMETER, PRIVATE :: n0 = MaxModArids + MaxSedArids + MaxBioArids

!
!1. Volume distributions
!-----------------------
!

INTEGER, PARAMETER :: &
     & iarr_part_ageconc = n0+1, iarr_part_ageconc_tot = n0+2, &
     & iarr_part_conc = n0+3, iarr_part_conc_tot = n0+4 

!
!2. Model grid arrays
!--------------------
!

INTEGER, PARAMETER :: &
& iarr_maskatc = n0+5, iarr_maskatu = n0+6, iarr_maskatv = n0+7, &
& iarr_vdifcoefpart = n0+8

!
!3. Particulate matter
!---------------------
!

INTEGER, PARAMETER :: &
& iarr_diamconc = n0+9, iarr_massconc = n0+10, iarr_rhoconc = n0+11, &
& iarr_volumeconc = n0+12

!
!4. Particle attributes
!----------------------
!

INTEGER, PARAMETER :: &
& iarr_part_age = n0+13, iarr_part_center = n0+14, iarr_part_displx = n0+15, &
& iarr_part_disply = n0+16, iarr_part_displz = n0+17, &
& iarr_part_drift_state = n0+18, iarr_part_icoord = n0+19, &
& iarr_part_jcoord = n0+20, iarr_part_kcoord = n0+21, &
& iarr_part_kdistype = n0+22, iarr_part_label = n0+23, &
& iarr_part_ntstart = n0+24, iarr_part_state = n0+25, &
& iarr_part_tstart = n0+26, iarr_part_xcoord = n0+27, iarr_part_xpos = n0+28, &
& iarr_part_ycoord = n0+29, iarr_part_ypos = n0+30, iarr_part_zpos = n0+31

!
!5. Cloud attributes
!-------------------
!

INTEGER, PARAMETER :: &
& iarr_cloud_data = n0+32, iarr_cloud_kdistype = n0+33, &
& iarr_cloud_label = n0+34, iarr_cloud_length = n0+35, &
& iarr_cloud_nopart = n0+36, iarr_cloud_noreleases = n0+37, &
& iarr_cloud_orientation = n0+38, iarr_cloud_state = n0+39, &
& iarr_cloud_thick = n0+40, iarr_cloud_width = n0+41


END MODULE partids
