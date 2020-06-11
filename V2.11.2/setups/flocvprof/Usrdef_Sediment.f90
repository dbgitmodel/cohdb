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
! Author - Boudewijn Decrop
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.11.1
!
! $Date: 2016-10-26 13:58:08 +0000 (Wed, 26 Oct 2016) $
!
! $Revision: 983 $
!
! Description - test case flocvprof
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
! Author - Boudewijn Decrop
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.11.1
!
! Description - test case flocvprof
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE sedpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=1) :: ctest


procname(pglev+1) = 'usrdef_sed_params'
CALL log_timer_in()

!
!2 Sediment parameters
!---------------------
!
!---erosion coefficient
parth_coef = 0.001

!---fractal dimension
ctest = runtitle(10:10)
SELECT CASE (ctest)
   CASE ('F'); nfrdim = 1.9
   CASE ('G'); nfrdim = 2.1
   CASE DEFAULT; nfrdim = 2.0
END SELECT

!---collision efficiency factor
SELECT CASE (ctest)
   CASE ('B'); agg_alpha = 0.05
   CASE ('C'); agg_alpha = 0.15
   CASE DEFAULT; agg_alpha = 0.1
END SELECT

!---efficiency for breaking
SELECT CASE (ctest)
   CASE ('D'); brk_es = 2.5E-04
   CASE ('E'); brk_es = 0.7E-04
   CASE DEFAULT; brk_es = 1.0E-04
END SELECT

!---fraction of microflocs by breaking 
SELECT CASE (ctest)
   CASE ('J'); brk_frac = 0.0
   CASE ('K'); brk_frac = 0.3
   CASE DEFAULT; brk_frac = 0.1
END SELECT

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
! Author - Boudewijn Decrop
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.11.1
!
! Description - test case flocvprof
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE sedarrays
USE sedpars
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: k
REAL :: cnump_ini, flocdia_ini, flocnc_ini


procname(pglev+1) = 'usrdef_sedics'
CALL log_timer_in()

massp(1) = rhos(1)*pi*dp(1)**3/6.0
cnump_ini = 0.32/massp(1)
flocdia_ini = 100.0E-06 
flocnc_ini = (flocdia_ini/dp(1))**nfrdim

k_110: DO k=1,nz
   WHERE (maskatc_int)
      cnump(1:ncloc,1:nrloc,k,1) = 0.9*cnump_ini
      cnump(1:ncloc,1:nrloc,k,3) = 0.1*cnump_ini
      cnump(1:ncloc,1:nrloc,k,2) = cnump(1:ncloc,1:nrloc,k,3)/flocnc_ini
   END WHERE
ENDDO k_110

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
! Author - Boudewijn Decrop
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.11.1
!
! Description - test case flocvprof
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

!
!*Local variables
!
CHARACTER (LEN=1) :: ctest


procname(pglev+1) = 'usrdef_sed_spec'
CALL log_timer_in()

ctest = runtitle(10:10)

!---diameter [m]
dp(1) = 15.0E-6

!---density [m^3/s]
SELECT CASE (ctest)
   CASE ('H'); rhos(1) = 1600.0
   CASE ('I'); rhos(1) = 1400.0   
   CASE DEFAULT; rhos(1) = 1800.0
END SELECT
      
!---(uniform) setling velocity [m/s]
bstres_cr_cst = 1.0/density_ref

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
