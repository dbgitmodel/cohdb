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
! *Usrdef_Harmonic_Analysis* Define harmonic analysis and output
!                            (example routines)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.1.2
!
! $Date: 2013-05-28 15:26:47 +0200 (Tue, 28 May 2013) $
!
! $Revision: 574 $
!
! Description - example file
!
! Reference -
!
! Routines - usrdef_anal_freqs, usrdef_anal_params, usrdef_anal0d_vals,
!            usrdef_anal2d_vals, usrdef_anal3d_vals
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_anal_freqs
!************************************************************************
!
! *usrdef_anal_freqs* Harmonic frequencies (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.1.2
!
! Description - if the index ID in index_anal is positive, names and values of
!               frequencies do need to be defined unles the user wants to use
!               non-default values
! Reference -
!
! Calling program - harmonic_analysis_init
!
! External calls -
!
! Module calls -
!
!************************************************************************
!
USE iopars
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: ifreq, iset


procname(pglev+1) = 'usrdef_anal_freqs'
CALL log_timer_in()

!
!1. Harmonic frequencies
!-----------------------
!

ifreq_110: DO ifreq=1,nofreqsanal
!  ---harmonic tidal index ID
   index_anal(ifreq) = ?
!  ---frequency name
   harm_freq_names(ifreq) = ?
!  ---frequency value [rad/s]
   harm_freq(ifreq) = ?
ENDDO ifreq_110

!
!2. Reference times
!------------------
!

iset_210: DO iset=1,nosetsanal
   cdate_time_ref(iset) = cdatetime_undef
ENDDO iset_210

!
!3. Specifiers for harmonic analysis
!-----------------------------------
!

iset_310: DO iset=1,nosetsanal

!  ---number of frequencies per set
   nofreqsharm(iset) = ?

!  ---frequency indices
   ifreq_311: DO ifreq=1,nofreqsanal
      ifreqsharm(iset,ifreq) = ?
   ENDDO ifreq_311

ENDDO iset_310

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_anal_freqs

!========================================================================

SUBROUTINE usrdef_anal_params
!************************************************************************
!
! *usrdef_anal_params* Specifiers for harmonic output (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.1.2
!
! Description - if %ivarid is defined (non-zero value), the %f90_name,
!               %long_name, %units and %vector_name attributes do not need to
!               be defined, unless the user wants to use non-default values
!             - variables need to be defined in the order of increasing ranks:
!               0, 2, 3
!
! Reference -
!
! Calling program - harmonic_analysis_init
!
! External calls - 
!
! Module calls -
!
!************************************************************************
!
USE iopars
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: ifreq, iset, istat, ivar, nostats


procname(pglev+1) = 'usrdef_anal_params'
CALL log_timer_in()

!
!1. Variable attributes
!----------------------
!

ivar_110: DO ivar=1,novarsanal
!  ---variable id
   analvars(ivar)%ivarid  = 0
!  ---variable rank (0/2/3)
   analvars(ivar)%nrank = ?
!  ---variable number
   analvars(ivar)%numvar = -1
!  ---output operator (oopt_null,oopt_max,oopt_min,oopt_mean,oopt_klev,oopt_dep)
   analvars(ivar)%oopt = oopt_null
!  ---vertical level for 2-D output of 3-D variables
   analvars(ivar)%klev = 0
!  ---vertical depth for 2-D output of 3-D variables
   analvars(ivar)%dep = 0.0
   IF (tsrvars(ivar)%ivarid.GT.0) THEN
!     ---fortran name
      analvars(ivar)%f90_name = ?
!     ---long description name
      analvars(ivar)%long_name = ?
!     ---variable unit
      analvars(ivar)%units = ?
!     ---vector name (vector quantities only)
      analvars(ivar)%vector_name = ?
   ENDIF
ENDDO ivar_110

!
!2. Variable indices
!-------------------
!

iset_210: DO iset=1,nosetsanal
   ivarsanal(iset,:) = ?
ENDDO iset_210

!
!3. Elliptic variables
!---------------------
!

iset_310: DO iset=1,nosetsanal
!  ---variable numbers
   ivarsell(iset,:) = ?
!  ---components of elliptic vector (2-D)
   ivecell2d(iset,1:2) = ?
!  ---components of elliptic vector (3-D)
   ivecell3d(iset,1:2) = ?
ENDDO iset_310

!
!4. File attributes
!-------------------
!
! 4.1 Residuals and external grid file
!-------------------------------------
!

iset_410: DO iset=1,nosetsanal

!  ---defined
   res0d(iset)%defined = .FALSE.; res2d(iset)%defined = .FALSE.
   res3d(iset)%defined = .FALSE.; analgrd(iset)%defined = .FALSE.

!  ---file format ('A','U','N')
   res0d(iset)%form = 'U'; res2d(iset)%form = 'U'
   res3d(iset)%form = 'U'; analgrd(iset)%form = 'U'

!  ---file name
   res0d(iset)%filename = ''; res2d(iset)%filename = ''
   res3d(iset)%filename = ''; analgrd(iset)%filename = ''

!  ---info file
   res0d(iset)%info = .TRUE.; res2d(iset)%info = .TRUE.
   res3d(iset)%info = .TRUE.; analgrd(iset)%info = .TRUE.

ENDDO iset_410

!
!4.2 Amplitudes, phases and elliptic parameters
!----------------------------------------------
!

iset_420: DO iset=1,nosetsanal
ifreq_420: DO ifreq=1,nofreqsharm(iset)

!  ---defined
   amp0d(iset,ifreq)%defined = .FALSE.; amp2d(iset,ifreq)%defined = .FALSE.
   amp3d(iset,ifreq)%defined = .FALSE.; pha0d(iset,ifreq)%defined = .FALSE.
   pha2d(iset,ifreq)%defined = .FALSE.; pha3d(iset,ifreq)%defined = .FALSE.
   ell2d(iset,ifreq)%defined = .FALSE.; ell3d(iset,ifreq)%defined = .FALSE.

!  ---file format ('A','U','N')
   amp0d(iset,ifreq)%form = 'U'; amp2d(iset,ifreq)%form = 'U'
   amp3d(iset,ifreq)%form = 'U'; pha0d(iset,ifreq)%form = 'U'
   pha2d(iset,ifreq)%form = 'U'; pha3d(iset,ifreq)%form = 'U'
   ell2d(iset,ifreq)%form = 'U'; ell3d(iset,ifreq)%form = 'U'

!  ---file name
   amp0d(iset,ifreq)%filename = ''; amp2d(iset,ifreq)%filename = ''
   amp3d(iset,ifreq)%filename = ''; amp3d(iset,ifreq)%filename = ''
   pha2d(iset,ifreq)%filename = ''; pha3d(iset,ifreq)%filename = ''
   ell2d(iset,ifreq)%filename = ''; ell3d(iset,ifreq)%filename = ''

!  ---info file
   amp0d(iset,ifreq)%info = .TRUE.; amp2d(iset,ifreq)%info = .TRUE.
   amp3d(iset,ifreq)%info = .TRUE.; pha0d(iset,ifreq)%info = .TRUE.
   pha2d(iset,ifreq)%info = .TRUE.; pha3d(iset,ifreq)%info = .TRUE.
   ell2d(iset,ifreq)%info = .TRUE.; ell3d(iset,ifreq)%info = .TRUE.

ENDDO ifreq_420
ENDDO iset_420

!
!5. Output grid (space/time) attributes
!--------------------------------------
!

iset_510: DO iset=1,nosetsanal
!  ---regular or irregular (stations) grid
   analgpars(iset)%gridded = .TRUE.
!  ---separate grid file
   analgpars(iset)%grid_file = .FALSE.
!  ---apply land mask (regular grid)
   analgpars(iset)%land_mask = .FALSE.
!  ---time-dependent spatial grid
   analgpars(iset)%time_grid = .FALSE.
!  ---time resolution (start/end/increment)
   analgpars(iset)%tlims(1:3) = ?
!  ---date/time of first output (optional)
   analgpars(iset)%startdate = cdatetime_undef
!  ---date/time of last output (optional)
   analgpars(iset)%enddate = cdatetime_undef
!  ---reference date in case of a numerical time coordinate
   analgpars(iset)%refdate = cdatetime_undef
!  ---grid dimension (0/2/3)
   analgpars(iset)%nodim = ?
!  ---number of stations (irregular grid)
   analgpars(iset)%nostats = ?
!  ---time format (0/1/2/3/4/5/6/7)
   analgpars(iset)%time_format = 0
!  ---resolution in X-direction (regular grid)
   analgpars(iset)%xlims(1:3) = ?
!  ---resolution in Y-direction (regular grid)
   analgpars(iset)%ylims(1:3) = ?
!  ---resolution in  vertical direction 
   analgpars(iset)%zlims(1:3) = ?

ENDDO iset_510

!
!6. Station attributes
!---------------------
!

istat_610: DO istat=1,nostatsanal
!  ---X-index
   analstatlocs(istat)%ipos = ?
!  ---Y-index
   analstatlocs(istat)%jpos = ?
!  ---station name
   analstatlocs(istat)%name = ?
ENDDO istat_610

!
!7. Station labels
!-----------------
!

iset_710: DO iset=1,nosetsanal
   nostats = analgpars(iset)%nostats
   lstatsanal(iset,1:nostats) = ?
ENDDO iset_710

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_anal_params

!========================================================================

SUBROUTINE usrdef_anal0d_vals(out0ddat,n0vars)
!************************************************************************
!
! *usrdef_anal0d_vals* 0-D output data (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.1.2
!
! Description - define output only for variables whose key id is set to zero
!
! Reference -
!
! Calling program - harmonic_analysis_update
!
! External calls - 
!
! Module calls-
!
!************************************************************************
!
USE iopars

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: n0vars
REAL, INTENT(OUT), DIMENSION(n0vars) :: out0ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out0ddat* REAL    output data
!*n0vars*   INTEGER number of 0-D variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar


ivar_110: DO ivar=1,n0vars
   out0ddat(ivar) = ?
ENDDO ivar_110

      
RETURN

END SUBROUTINE usrdef_anal0d_vals

!========================================================================

SUBROUTINE usrdef_anal2d_vals(out2ddat,i,j,n2vars)
!************************************************************************
!
! *usrdef_anal2d_vals* 2-D output data (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.1.2
!
! Description - define output only for variables whose key id is set to zero
!
! Reference -
!
! Calling program - harmonic_analysis_update
!
! External calls - 
!
! Module calls -
!
!************************************************************************
!

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, n2vars
REAL, INTENT(OUT), DIMENSION(n2vars) :: out2ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out2ddat* REAL    output data
!*i*        INTEGER X-index of output location
!*j*        INTEGER Y-index of output location
!*n2vars*   INTEGER number of 2-D variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar


ivar_110: DO ivar=1,n2vars
   out2ddat(ivar) = ?
ENDDO ivar_110


RETURN

END SUBROUTINE usrdef_anal2d_vals

!========================================================================

SUBROUTINE usrdef_anal3d_vals(out3ddat,i,j,k,n3vars)
!************************************************************************
!
! *usrdef_anal3d_vals* 3-D output data (example routine)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Harmonic_Analysis.f90  V2.1.2
!
! Description - define output only for variables whose key id is set to zero
!
! Reference -
!
! Calling program - harmonic_analysis_update
!
! External calls - 
!
! Module calls-
!
!************************************************************************
!
  
IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, k, n3vars
REAL, INTENT(OUT), DIMENSION(n3vars) :: out3ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out3ddat* REAL    output data
!*i*        INTEGER X-index of output location
!*j*        INTEGER Y-index of output location
!*k*        INTEGER Vertical index of output location
!*n3vars*   INTEGER number of 3-D variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar


ivar_110: DO ivar=1,n3vars
   out3ddat(ivar) = ?
ENDDO ivar_110


RETURN

END SUBROUTINE usrdef_anal3d_vals
