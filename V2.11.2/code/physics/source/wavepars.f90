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


MODULE wavepars
!************************************************************************
!
! *wavepars* Surface wave parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)wavepars.f90  V2.9
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description - 
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

INTEGER, PARAMETER :: MaxCohtoWavvars = 5, MaxWavtoCohvars = 13
CHARACTER (LEN=lenname), DIMENSION(MaxCohtoWavvars) :: cohtowav_names
CHARACTER (LEN=lenname), DIMENSION(MaxWavtoCohvars) :: wavtocoh_names
CHARACTER (LEN=lentime) :: CEndDateTimeWav, CStartDateTimeWav
INTEGER, DIMENSION(7) ::  IEndDateTimeWav, IStartDateTimeWav
INTEGER :: ncwav, ncwavloc, ncwavloc_ext, nc1wavloc, nc2wavloc, nrwav, &
         & nrwavloc, nr1wavloc, nr2wavloc, nrwavloc_ext, haloE, haloN, haloS, &
         & haloW
REAL :: deltwav


SAVE

!
! Name             Type    Purpose
!------------------------------------------------------------------------------
!*CEndDateTimeWav   CHAR    End date/time for wave input
!*CStartDateTimeWav CHAR    Start date/time for wave input
!*cohtowav_names*   CHAR    Names of the variables send by COHERENS to the wave
!                           model
!*deltwav*          REAL    Time step for wave input/coupling                [s]
!*IEndDateTimeWav   INTEGER End date/time for wave input/coupling
!*IStartDateTimeWav INTEGER Start date/time for wave input/coupling
!*MaxCohtoWavvars*  INTEGER Maximum number of wave variables used in the
!                           coupling
!*MaxWavtoCohvars*  INTEGER Maximum number of COHERENS variables used in the
!                           coupling
!*ncwav*            INTEGER Global number of wave model grid cells in the
!                           X-direction
!*ncwavloc*         INTEGER Local number of wave model grid cells in the
!                           X-direction
!*ncwav_ext*        INTEGER Global number (including halos) of wave model grid
!                           cells in the X-direction
!*nc1wavloc*        INTEGER Global index of the first column (excluding halos)
!                           on the local domain wave grid
!*nc2wavloc*        INTEGER Global index of the last column (excluding halos)
!                           on the local domain wave grid
!*nrwav*            INTEGER Global number of wave model grid cells in the
!                           Y-direction
!*nrwavloc*         INTEGER Local number of wave model grid cells in the
!                           Y-direction
!*nrwav_ext*        INTEGER Global number (including halos) of wave model grid
!                           cells in the Y-direction
!*nr1wavloc*        INTEGER Global index of the first row (excluding halos) on
!                           the local domain wave grid
!*nr2wavloc*        INTEGER Global index of the last row (excluding halos) on
!                           the local domain wave grid
!*wavtocoh_names*   CHAR    Names of the variables send by the wave model to
!                           COHERENS
!*whaloE*           UNTEGER Halo size on the eastern boundary of the local
!                           domain grid
!*whaloN*           UNTEGER Halo size on the northern boundary of the local
!                           domain grid
!*whaloS*           UNTEGER Halo size on the southern boundary of the local
!                           domain grid
!*whaloW*           UNTEGER Halo size on the western boundary of the local
!                           domain grid
!
!************************************************************************
!

END MODULE wavepars
