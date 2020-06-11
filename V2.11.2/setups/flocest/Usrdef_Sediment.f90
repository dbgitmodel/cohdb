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

!************************************************************************
!
! *Usrdef_Sediment* User-defined model setup
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.11.1
!
! $Date: 2016-10-26 13:58:08 +0000 (Wed, 26 Oct 2016) $
!
! $Revision: 983 $
!
! Description - test case flocest
!
! Reference -
!
! Routines - usrdef_sed_params, usrdef_morph_params, usrdef_dar_params,
!            usrdef_sedics, usrdef_morphics, usrdef_sed_spec, usrdef_dar_spec
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_sed_params
!************************************************************************
!
! *usrdef_sed_params* Define parameters sediment transport
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.11.1
!
! Description - test case flocest
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE sedpars
USE sedswitches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_sed_params'
CALL log_timer_in()

!
!1. Sediment switches
!--------------------
!
!---sediment density stratification
iopt_sed_dens_grad = MERGE(0,1,runtitle(8:8).EQ.'E')

!
!2 Sediment parameters
!---------------------
!
!---erosion coefficient
parth_coef = 0.001

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_sed_params

!=========================================================

SUBROUTINE usrdef_morph_params
!************************************************************************
!
! *usrdef_morph_params* Define parameters for the morphological model
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_morph_params

!========================================================================

SUBROUTINE usrdef_dar_params
!************************************************************************
!
! *usrdef_dar_params* Define parameters for dredging and relocation
!
! Author - Alexander Breugem and  Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_dar_params

!========================================================================

SUBROUTINE usrdef_sedics
!************************************************************************
!
! *usrdef_sedics* Define arrays for initial conditions sediment transport
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.11.1
!
! Description - test case flocest
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE gridpars
USE iopars
USE sedarrays
USE sedpars
USE syspars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, k
REAL :: cnump_ini, flocdia_ini, flocnc_ini
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: cnumpglb 


IF (runtitle(8:8).NE.'B') RETURN

procname(pglev+1) = 'usrdef_sedics'
CALL log_timer_in()
   
ALLOCATE (cnumpglb(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo,nz,3),STAT=errstat)
CALL error_alloc('cnumpglb',4,(/ncloc+2*nhalo,nr+2*nhalo,nz,3/),kndrtype)
cnumpglb = 0.0

massp(1) = rhos(1)*pi*dp(1)**3/6.0
cnump_ini = 0.3/massp(1)
flocdia_ini = 100.0E-06 
flocnc_ini = (flocdia_ini/dp(1))**nfrdim

i_110: DO i=1,nc-1
   IF (i.LE.50) THEN
      cnumpglb(i,1:nr-1,:,1) = cnump_ini
      cnumpglb(i,1:nr-1,:,2:3) = 0.0
   ELSE
      cnumpglb(i,1:nr-1,:,1) = 0.0
      cnumpglb(i,1:nr-1,:,3) = cnump_ini
      cnumpglb(i,1:nr-1,:,2) = cnump_ini/flocnc_ini 
   ENDIF
ENDDO i_110

cnump = cnumpglb(nc1loc-nhalo:nc2loc+nhalo,nr1loc-nhalo:nr2loc+nhalo,:,:)

DEALLOCATE (cnumpglb)

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_sedics

!=======================================================================

SUBROUTINE usrdef_morphics
!************************************************************************
!
! *usrdef_morphics* Define initial conditions for the morphological model
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_morphics

!========================================================================

SUBROUTINE usrdef_sed_spec
!************************************************************************
!
! *usrdef_sed_spec* Define characteristics of sediment fractions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.11.1
!
! Description - test case flocest
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE physpars
USE sedarrays
USE sedpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_sed_spec'
CALL log_timer_in()

!---diameter [m]
dp(1) = MERGE(30.0,15.0,runtitle(8:8).EQ.'C')
dp(1) = dp(1)*1.0E-06

!---density [m^3/s]
rhos(1) = 1600.0
      
!---(uniform) critical shear stress [m2/s2]
bstres_cr_cst = MERGE(0.002,0.001,runtitle(8:8).EQ.'D')

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_sed_spec

!========================================================================

SUBROUTINE usrdef_dar_spec
!************************************************************************
!
! *usrdef_dar_spec* Define site and ship properties for dredging and relocation
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_dar_spec
