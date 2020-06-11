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
! *Surface_Waves* Define surface wave data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Waves.f90  V2.11.2
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description - 
!
! Reference -
!
! Subroutines - update_wave_params, wave_input, wave_sources_momentum
!
!************************************************************************
!

!========================================================================

SUBROUTINE update_wave_params
!************************************************************************
!
! *update_wave_params* update wave parameters not obtained from input
!
! Author - Alexander Breugem and Patrick Luyten (IMDC)
!
! Version - @(COHERENS)Surface_Waves.f90  V2.11.2
!
! Description - wave number obtained from the dispersion relation using the
!               approximate method of Hunt (1979)
!
! Reference -
!
! Calling program -
!
! External calls - wave_boundary_layer
!
! Module calls - Carr_at_U, Carr_at_V, exchange_mod
!
!************************************************************************
!
USE depths
USE fluxes  
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
use timepars
USE wavevars
USE array_interp, ONLY: Carr_at_U, Carr_at_V
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: info = .TRUE.
INTEGER :: i, j, k, npcc
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
REAL, PARAMETER :: epsmax = 50.0
REAL :: alpha, xlim
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc1, array2dc2, wnumdep
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: zprof


procname(pglev+1) = 'update_wave_params'
CALL log_timer_in(npcc)

!
!1. Initialise
!-------------
!
!---allocate
ALLOCATE (array2dc1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dc2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (zprof(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('zprof',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (wnumdep(ncloc,nrloc),STAT=errstat)
CALL error_alloc('wnumdep',2,(/ncloc,nrloc/),kndrtype)

zprof = 0.0

!
!2. Wave number
!--------------
!
!---solve dispersion relation
WHERE (maskatc_int)
   array2dc1 = deptotatc(1:ncloc,1:nrloc)*wavefreq(1:ncloc,1:nrloc)**2/&
             & gaccatc(1:ncloc,1:nrloc)
ELSEWHERE
   array2dc1 = 0.0
END WHERE
j_210: DO j=1,nrloc
i_210: DO i=1,ncloc
   xlim = SQRT(array2dc1(i,j))
   IF (.NOT.maskatc_int(i,j)) THEN
      wavenum(i,j) = 0.0
   ELSEIF (xlim.LE.1.0E-06) THEN
      wavenum(i,j) = wavefreq(i,j)/SQRT(gaccatc(i,j)*deptotatc(i,j))
   ELSEIF (xlim.LT.2.5) THEN
      alpha = wavefreq(i,j)**2*deptotatc(i,j)/gaccatc(i,j)
      wavenum(i,j) = wavefreq(i,j)*SQRT((alpha+1/(1.0+0.666*alpha+&
                              & 0.445*alpha**2-0.105*alpha**3+0.272*alpha**4))&
                              & /(gaccatc(i,j)*deptotatc(i,j)))
   ELSE
      wavenum(i,j) = wavefreq(i,j)**2/gaccatc(i,j)
   ENDIF
ENDDO i_210
ENDDO j_210

!---work space array
WHERE (maskatc_int)
   wnumdep = MIN(wavenum(1:ncloc,1:nrloc)*deptotatc(1:ncloc,1:nrloc),epsmax)
ELSEWHERE
   wnumdep = 0.0
END WHERE

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   lbounds(1:2) = 0; nhexch = (/1,0,1,0/)
   CALL exchange_mod(wavenum,lbounds(1:2),nhexch,iarr_wavenum)
ENDIF

!
!3. Wave orbital velocities and excursion amplitude
!--------------------------------------------------
!

IF (iopt_waves_form.EQ.1) THEN

!  ---define
   WHERE (wnumdep.GT.0.0)
      wavevel(1:ncloc,1:nrloc) = 0.5*waveheight(1:ncloc,1:nrloc)*&
           & wavefreq(1:ncloc,1:nrloc)/SINH(wnumdep)
      waveexcurs(1:ncloc,1:nrloc) = 0.5*waveheight(1:ncloc,1:nrloc)/&
                                     & SINH(wnumdep)
   ELSEWHERE
      wavevel(1:ncloc,1:nrloc) = 0.0
      waveexcurs(1:ncloc,1:nrloc) = 0.0
   END WHERE

!---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds(1:2) = 0; nhexch = (/1,0,1,0/)
      CALL exchange_mod(wavevel,lbounds(1:2),nhexch,iarr_wavevel)
      CALL exchange_mod(waveexcurs,lbounds(1:2),nhexch,iarr_waveexcurs)
   ENDIF

ENDIF

!
!4. Stokes velocities
!--------------------
!

IF (iopt_waves_curr.EQ.1) THEN

!
!4.1 Depth-mean current
!----------------------
!

   IF (iopt_waves_form.EQ.1) THEN
      WHERE (wnumdep.GT.0.0)
         array2dc1 = 0.125*wavefreq(1:ncloc,1:nrloc)*&
                   & waveheight(1:ncloc,1:nrloc)**2/&
                   & TANH(wnumdep)/deptotatc(1:ncloc,1:nrloc)
         umstokesatc(1:ncloc,1:nrloc) = COS(wavedir(1:ncloc,1:nrloc))*array2dc1
         vmstokesatc(1:ncloc,1:nrloc) = SIN(wavedir(1:ncloc,1:nrloc))*array2dc1
      ELSEWHERE
         umstokesatc(1:ncloc,1:nrloc) = 0.0
         vmstokesatc(1:ncloc,1:nrloc) = 0.0
      END WHERE
   ENDIF

!
!4.2 3-D current
!---------------
!

   IF (iopt_grid_nodim.EQ.3) THEN

      WHERE (wnumdep.GT.0.0)
         array2dc1 = 2.0*wnumdep
         array2dc2 = array2dc1/SINH(array2dc1)
      END WHERE
      k_421: DO k=1,nz
         WHERE (wnumdep.GT.0.0)
            zprof(:,:,k) = array2dc2*&
                    & COSH(MIN(array2dc1*gscoordatc(1:ncloc,1:nrloc,k),epsmax))
         END WHERE
      ENDDO k_421
      WHERE (wnumdep.GT.0.0)
         array2dc1 = SUM(zprof*delzatc(1:ncloc,1:nrloc,:),DIM=3)         
         array2dc1 = array2dc1/deptotatc(1:ncloc,1:nrloc)
      END WHERE
      k_422: DO k=1,nz
         WHERE (wnumdep.GT.0.0)
            zprof(:,:,k) = zprof(:,:,k)/array2dc1
         END WHERE
      ENDDO k_422
      k_423: DO k=1,nz
         WHERE (maskatc_int)
            ustokesatc(1:ncloc,1:nrloc,k) = zprof(:,:,k)*&
                                          & umstokesatc(1:ncloc,1:nrloc)
            vstokesatc(1:ncloc,1:nrloc,k) = zprof(:,:,k)*&
                                          & vmstokesatc(1:ncloc,1:nrloc)
         END WHERE
      ENDDO k_423
      
   ENDIF

ENDIF

!
!5. Dissipation forces
!---------------------
!

IF (iopt_waves_dissip.EQ.1.AND.iopt_grid_nodim.EQ.3) THEN

!
!5.1 Surface
!-----------
!

   WHERE (wnumdep.GT.0.0)
      array2dc1 = wave_penetration_surf*wnumdep
      array2dc2 = array2dc1/SINH(wnumdep)
   END WHERE
   k_511: DO k=1,nz
      WHERE (wnumdep.GT.0.0)
         zprof(:,:,k) = array2dc2*&
                    & COSH(MIN(array2dc1*gscoordatc(1:ncloc,1:nrloc,k),epsmax))
      END WHERE
   ENDDO k_511
   WHERE (wnumdep.GT.0.0)
      array2dc1 = SUM(zprof*delzatc(1:ncloc,1:nrloc,:),DIM=3)
      array2dc1 = array2dc1/deptotatc(1:ncloc,1:nrloc)
   END WHERE
   k_512: DO k=1,nz
      WHERE (wnumdep.GT.0.0)
         zprof(:,:,k) = zprof(:,:,k)/array2dc1
      END WHERE
   ENDDO k_512
   k_513: DO k=1,nz
      WHERE (maskatc_int)
         uswdissipatc(1:ncloc,1:nrloc,k) = zprof(:,:,k)*&
              & (umswdissipatc(1:ncloc,1:nrloc)/deptotatc(1:ncloc,1:nrloc))
         vswdissipatc(1:ncloc,1:nrloc,k) = zprof(:,:,k)*&
              & (vmswdissipatc(1:ncloc,1:nrloc)/deptotatc(1:ncloc,1:nrloc))
      END WHERE
   ENDDO k_513

!
!5.2 Bottom
!----------
!

   WHERE (wavethickatc.GT.0.0)
      array2dc1 = deptotatc(1:ncloc,1:nrloc)/(wave_penetration_bed*wavethickatc)
      array2dc2 = array2dc1/SINH(MIN(array2dc1,epsmax))
   END WHERE
   k_521: DO k=1,nz
      WHERE (wavethickatc.GT.0.0)
         zprof(:,:,k) = array2dc2*&
              & COSH(MIN(array2dc1*(1.0-gscoordatc(1:ncloc,1:nrloc,k)),epsmax))
      END WHERE
   ENDDO k_521
   WHERE (wavethickatc.GT.0.0)
      array2dc1 = SUM(zprof*delzatc(1:ncloc,1:nrloc,:),DIM=3)
      array2dc1 = array2dc1/deptotatc(1:ncloc,1:nrloc)
   END WHERE
   k_522: DO k=1,nz
      WHERE (wavethickatc.GT.0.0)
         zprof(:,:,k) = zprof(:,:,k)/array2dc1
      END WHERE
   ENDDO k_522
   k_523: DO k=1,nz
      WHERE (maskatc_int)
         ubwdissipatc(1:ncloc,1:nrloc,k) = zprof(:,:,k)*&
              & (umbwdissipatc(1:ncloc,1:nrloc)/deptotatc(1:ncloc,1:nrloc))
         vbwdissipatc(1:ncloc,1:nrloc,k) = zprof(:,:,k)*&
              & (vmbwdissipatc(1:ncloc,1:nrloc)/deptotatc(1:ncloc,1:nrloc))
      END WHERE
   ENDDO k_523

ENDIF

!
!6. Wave induced pressure
!------------------------
!

IF (iopt_waves_pres.EQ.1.AND.iopt_waves_form.EQ.1) THEN

   WHERE (wnumdep.GT.0.0)
      wavepres(1:ncloc,1:nrloc) = 0.125*gaccatc(1:ncloc,1:nrloc)*&
                 & wavenum(1:ncloc,1:nrloc)*waveheight(1:ncloc,1:nrloc)**2/&
                 & SINH(2.0*wnumdep)
   ELSEWHERE
      wavepres(1:ncloc,1:nrloc) = 0.0
   END WHERE

ENDIF

!
!7. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1.AND.iopt_waves_curr.EQ.1) THEN

!
!7.1 2-D arrays
!--------------
!

   IF (iopt_waves_form.EQ.1) THEN
!     ---Stokes drift
      lbounds(1:2) = (/1-nhalo,0/); nhexch = (/nhalo,nhalo,1,0/)
      CALL exchange_mod(umstokesatc,lbounds(1:2),nhexch,iarr_umstokesatc)
      lbounds(1:2) = (/0,1-nhalo/); nhexch = (/1,0,nhalo,nhalo/)
      CALL exchange_mod(vmstokesatc,lbounds(1:2),nhexch,iarr_vmstokesatc)
!     ---wave pressure
      IF (iopt_waves_pres.EQ.1) THEN
         lbounds(1:2) = 0; nhexch = (/1,0,1,0/)
         CALL exchange_mod(wavepres,lbounds(1:2),nhexch,iarr_wavepres)
      ENDIF
   ENDIF

!
!7.2 3-D arrays
!--------------
!

   IF (iopt_grid_nodim.EQ.3) THEN
!     ---Stokes drift
      lbounds = (/1-nhalo,0,1/); nhexch = (/nhalo,nhalo,1,0/)
      CALL exchange_mod(ustokesatc,lbounds,nhexch,iarr_ustokesatc)
      lbounds = (/0,1-nhalo,1/); nhexch = (/1,0,nhalo,nhalo/)
      CALL exchange_mod(vstokesatc,lbounds,nhexch,iarr_vstokesatc)
!     ---dissipation forces
      IF (iopt_waves_dissip.EQ.1) THEN
         lbounds = (/0,1,1/); nhexch = (/1,0,0,0/)
         CALL exchange_mod(uswdissipatc,lbounds,nhexch,iarr_uswdissipatc)
         CALL exchange_mod(ubwdissipatc,lbounds,nhexch,iarr_ubwdissipatc)
         lbounds = (/1,0,1/); nhexch = (/0,0,1,0/)
         CALL exchange_mod(vswdissipatc,lbounds,nhexch,iarr_vswdissipatc)
         CALL exchange_mod(vbwdissipatc,lbounds,nhexch,iarr_vbwdissipatc)
      ENDIF
   ENDIF

ENDIF

!
!8. Interpolate at velocity nodes
!--------------------------------
!

IF (iopt_waves_curr.EQ.1) THEN

!  ---2-D Stokes drift
   IF (iopt_waves_form.EQ.1) THEN
      CALL Carr_at_U(umstokesatc(0:ncloc,1:nrloc),umstokesatu(1:ncloc,1:nrloc),&
                   & 1,1,(/0,1,nz/),(/ncloc,nrloc,nz/),1,iarr_umstokesatc,info)
      CALL Carr_at_V(vmstokesatc(1:ncloc,0:nrloc),vmstokesatv(1:ncloc,1:nrloc),&
                   & 1,1,(/1,0,nz/),(/ncloc,nrloc,nz/),1,iarr_vmstokesatc,info)
   ENDIF

!  ---3-D Stokes drift
   IF (iopt_grid_nodim.EQ.3) THEN
      CALL Carr_at_U(ustokesatc(0:ncloc,1:nrloc,:),&
                   & ustokesatu(1:ncloc,1:nrloc,:),1,1,(/0,1,1/),&
                   & (/ncloc,nrloc,nz/),1,iarr_ustokesatc,info)
      CALL Carr_at_V(vstokesatc(1:ncloc,0:nrloc,:),&
                   & vstokesatv(1:ncloc,1:nrloc,:),1,1,(/1,0,1/),&
                   & (/ncloc,nrloc,nz/),1,iarr_vstokesatc,info)
   ENDIF

ENDIF

!
!9. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1.AND.iopt_waves_curr.EQ.1) THEN

!  ---2-D Stokes drift
   IF (iopt_waves_form.EQ.1) THEN
      lbounds(1:2) = (/1-nhalo,0/); nhexch = (/nhalo,nhalo,1,0/)
      CALL exchange_mod(umstokesatu,lbounds(1:2),nhexch,iarr_umstokesatu)
      lbounds(1:2) = (/0,1-nhalo/); nhexch = (/1,0,nhalo,nhalo/)
      CALL exchange_mod(vmstokesatv,lbounds(1:2),nhexch,iarr_vmstokesatv)
   ENDIF

!  ---3-D Stokes drift
   IF (iopt_grid_nodim.EQ.3) THEN
      lbounds = (/1-nhalo,0,1/); nhexch = (/nhalo,nhalo,1,0/)
      CALL exchange_mod(ustokesatu,lbounds,nhexch,iarr_ustokesatu)
      lbounds = (/0,1-nhalo,1/); nhexch = (/1,0,nhalo,nhalo/)
      CALL exchange_mod(vstokesatv,lbounds,nhexch,iarr_vstokesatv)
   ENDIF

ENDIF

!
!10. Stokes transports
!---------------------
!

IF (iopt_waves_curr.EQ.1) THEN

   WHERE (node2du(1:ncloc+1,1:nrloc).GT.0)
      udstokesatu = deptotatu(1:ncloc+1,1:nrloc)*umstokesatu(1:ncloc+1,1:nrloc)
   END WHERE
   WHERE (node2dv(1:ncloc,1:nrloc+1).GT.0)
      vdstokesatv = deptotatv(1:ncloc,1:nrloc+1)*vmstokesatv(1:ncloc,1:nrloc+1)
   END WHERE

ENDIF

!
!11. Deallocate
!--------------
!

DEALLOCATE (array2dc1,array2dc2,wnumdep,zprof)

CALL log_timer_out(npcc,itm_waves)


RETURN

END SUBROUTINE update_wave_params

!==============================================================================

SUBROUTINE wave_input
!************************************************************************
!
! *wave_input* update wave input (if needed)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Waves.f90  V2.11.1
!
! Description - in case of coupling with a wave model, a two-way exchange is
!               performed  
!
! Reference -
!
! Calling program - coherens_main, initialise_model
!
! External calls - coh2wav_recv, define_surface_input_grid, update_surface_data,
!                  update_wave_params, wave_number, write_surface_data
!
! Module calls - Carr_at_U, Carr_at_V, combine_mod, error_alloc,
!                error_alloc_struc, exchange_mod, hrelativecoords_init,
!                set_modfiles_atts
!
!************************************************************************
!
USE datatypes
USE grid
USE gridpars
USE iopars
USE modids
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Carr_at_U, Carr_at_V
USE datatypes_init, ONLY: hrelativecoords_init
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE modvars_routines, ONLY: set_modfiles_atts
USE paral_comms, ONLY: combine_mod, exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: info = .TRUE.
INTEGER :: npcc
INTEGER, SAVE :: ivar, nhtype, novars, n1dat, n1grd, n2dat, n2grd
INTEGER, DIMENSION(2) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
INTEGER (KIND=kndilong), SAVE, DIMENSION(2) :: nosecsdat
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:) :: angle
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: wavevals, wavevalsglb
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: wavedat
TYPE (HRelativeCoords), SAVE, ALLOCATABLE, DIMENSION(:,:) :: wavegrid

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*nhtype*     INTEGER Type of data grid (0/1/2)
!*n1dat*      INTEGER X-dimension of data grid
!*n2dat*      INTEGER Y-dimension of data grid
!*novars*     INTEGER Number of input parameters
!*nosecsdat*  LONGINT No. of seconds between start of simulation and old/new
!                     data time
!*wavedat*    REAL    Wave data at old and new data time
!*wavevals*   REAL    Wave data interpolated in time and on model grid
!*wavegrid*   DERIVED Model grid relative coordinates with respect to wave grid
!
!------------------------------------------------------------------------------
!
!1. Return if no update required
!-------------------------------
!

IF (.NOT.wavestep) RETURN

procname(pglev+1) = 'wave_input'
CALL log_timer_in(npcc)

nhtype = surfacegrids(igrd_waves,1)%nhtype

!
!2. Initialise on first call
!---------------------------
!

IF (nt.EQ.0) THEN

!
!2.1 Without coupling
!--------------------
! 
   
   IF (iopt_waves_couple.EQ.0) THEN
   
!
!2.1.1 Surface grid
!------------------
!
!     ---grid parameters
      n1dat = surfacegrids(igrd_waves,1)%n1dat
      n2dat = surfacegrids(igrd_waves,1)%n2dat
      n1grd  = MERGE(ncloc,0,nhtype.GT.0.AND.nhtype.LT.4)
      n2grd  = MERGE(nrloc,0,nhtype.GT.0.AND.nhtype.LT.4)

!     ---allocate coordinate structure
      ALLOCATE (wavegrid(n1grd,n2grd),STAT=errstat)
      CALL error_alloc_struc('wavegrid',2,(/n1grd,n2grd/),'HRelativeCoords')
      CALL hrelativecoords_init(wavegrid,.FALSE.)

!     ---define surface grid
      IF (nhtype.GT.0.AND.nhtype.LT.4) THEN
         CALL define_surface_input_grid(igrd_waves,1,wavegrid)
      ENDIF

!
!2.1.2 Number of wave variables
!------------------------------
!

      CALL set_modfiles_atts(io_wavsur,1,1)
      novars = modfiles(io_wavsur,1,1)%novars

!
!2.1.3 Allocate memory for data arrays
!-------------------------------------
!

      ALLOCATE (wavedat(n1dat,n2dat,novars,2),STAT=errstat)
      CALL error_alloc('wavedat',4,(/n1dat,n2dat,novars,2/),kndrtype)
      wavedat = 0.0
      ALLOCATE (wavevals(ncloc,nrloc,novars),STAT=errstat)
      CALL error_alloc('wavevals',3,(/ncloc,nrloc,novars/),kndrtype)
      wavevals = 0.0
      ALLOCATE (angle(novars),STAT=errstat)
      CALL error_alloc('angle',1,(/novars/),kndlog)
      angle = .FALSE.; angle(3) = .TRUE.

!      
!2.2 With coupling
!-----------------
!      

   ELSEIF (modfiles(io_wavsur,1,2)%status.EQ.'W') THEN

!
!2.2.1 Data parameters
!---------------------
!
!     ---surface grid
      surfacegrids(igrd_waves,2)%n1dat = nc
      surfacegrids(igrd_waves,2)%n2dat = nr
      surfacegrids(igrd_waves,2)%nhtype = iopt_grid_htype

!     ---output times
      modfiles(io_wavsur,1,2)%tlims = modfiles(io_wavsur,1,1)%tlims

!     ---number of wave variables
      CALL set_modfiles_atts(io_wavsur,1,2)
      novars = modfiles(io_wavsur,1,2)%novars

!
!2.2.3 Allocate memory for data arrays
!-------------------------------------
!

      ALLOCATE (wavevals(ncloc,nrloc,novars),STAT=errstat)
      CALL error_alloc('wavevals',3,(/ncloc,nrloc,novars/),kndrtype)
      wavevals = 0.0
      ALLOCATE (wavevalsglb(nc,nr,novars),STAT=errstat)
      CALL error_alloc('wavevalsglb',3,(/nc,nr,novars/),kndrtype)
      wavevalsglb = 0.0

   ENDIF
   
ENDIF

!
!3. Update input parameters
!--------------------------
!

IF (iopt_waves_couple.EQ.0) THEN

!
!3.1 Coupling with the wave model is disabled
!--------------------------------------------
!

   CALL update_surface_data(io_wavsur,1,wavedat,wavevals,wavegrid,n1dat,&
                          & n2dat,novars,n1grd,n2grd,nosecsdat,angle)

!
!3.2 Wave variables
!-------------------
!
   ivar = 1
   waveheight(1:ncloc,1:nrloc) = wavevals(:,:,ivar)
   ivar = ivar + 1
   waveperiod(1:ncloc,1:nrloc) = wavevals(:,:,ivar)
   ivar = ivar + 1
   wavedir(1:ncloc,1:nrloc) = wavevals(:,:,ivar)
   IF (iopt_waves_form.EQ.2) THEN
      ivar = ivar + 1
      wavevel(1:ncloc,1:nrloc) = wavevals(:,:,ivar)
      ivar = ivar + 1
      waveexcurs(1:ncloc,1:nrloc) = wavevals(:,:,ivar)
   ENDIF
   IF (iopt_waves_curr.EQ.1.AND.iopt_waves_form.EQ.2) THEN
      ivar = ivar + 1
      umstokesatc(1:ncloc,1:nrloc) = wavevals(:,:,ivar)
      ivar = ivar + 1
      vmstokesatc(1:ncloc,1:nrloc) = wavevals(:,:,ivar)
      IF (iopt_waves_pres.EQ.1) THEN
         ivar = ivar + 1
         wavepres(1:ncloc,1:nrloc) = wavevals(:,:,ivar)
      ENDIF
   ENDIF
   IF (iopt_waves_dissip.EQ.1) THEN
      ivar = ivar + 1
      umswdissipatc(1:ncloc,1:nrloc) = wavevals(:,:,ivar)
      ivar = ivar + 1
      vmswdissipatc(1:ncloc,1:nrloc) = wavevals(:,:,ivar)
      ivar = ivar + 1
      umbwdissipatc(1:ncloc,1:nrloc) = wavevals(:,:,ivar)
      ivar = ivar + 1
      vmbwdissipatc(1:ncloc,1:nrloc) = wavevals(:,:,ivar)
   ENDIF

!
!3.2 Wave coupling enabled
!-------------------------
!

ELSE

!
!3.1.1 Sending data to wave model
!--------------------------------
!   

   CALL coh2wav_send
   IF (iopt_waves_couple.EQ.2) GOTO 1000
   
!
!3.1.2  Receiving data from wave model
!-------------------------------------
!   

   CALL coh2wav_recv

!
!3.1.3 Write data (if requested)
!-------------------------------   
!

   IF (modfiles(io_wavsur,1,2)%status.EQ.'W') THEN

!     ---store data locally
      ivar = 1
      wavevals(:,:,ivar) = waveheight(1:ncloc,1:nrloc)
      ivar = ivar + 1
      wavevals(:,:,ivar) = waveperiod(1:ncloc,1:nrloc)
      ivar = ivar + 1
      wavevals(:,:,ivar) = wavedir(1:ncloc,1:nrloc)
      IF (iopt_waves_form.EQ.2) THEN
         ivar = ivar + 1
         wavevals(:,:,ivar) = wavevel(1:ncloc,1:nrloc)
         ivar = ivar + 1
         wavevals(:,:,ivar) = waveexcurs(1:ncloc,1:nrloc)
      ENDIF
      IF (iopt_waves_curr.EQ.1.AND.iopt_waves_form.EQ.2) THEN
         ivar = ivar + 1
         wavevals(:,:,ivar) = umstokesatc(1:ncloc,1:nrloc)
         ivar = ivar + 1
         wavevals(:,:,ivar) = vmstokesatc(1:ncloc,1:nrloc)
         IF (iopt_waves_pres.EQ.1) THEN
            ivar = ivar + 1
            wavevals(:,:,ivar) = wavepres(1:ncloc,1:nrloc)
         ENDIF
      ENDIF
      IF (iopt_waves_dissip.EQ.1) THEN
         ivar = ivar + 1
         wavevals(:,:,ivar) = umswdissipatc(1:ncloc,1:nrloc)
         ivar = ivar + 1
         wavevals(:,:,ivar) = vmswdissipatc(1:ncloc,1:nrloc)
         ivar = ivar + 1
         wavevals(:,:,ivar) = umbwdissipatc(1:ncloc,1:nrloc)
         ivar = ivar + 1
         wavevals(:,:,ivar) = vmbwdissipatc(1:ncloc,1:nrloc)
      ENDIF

!     ---combine on master
      CALL combine_mod(wavevalsglb,wavevals,(/1,1,1/),0,0.0)

!     ---write data
      CALL write_surface_data(io_wavsur,1,2,CDateTime,wavevalsglb,nc,nr,novars)

   ENDIF

ENDIF

!
!4. Related wave parameters
!--------------------------
!
!---frequency
WHERE (waveperiod(1:ncloc,1:nrloc).GT.0.0)
   wavefreq(1:ncloc,1:nrloc) = twopi/waveperiod(1:ncloc,1:nrloc)
ELSEWHERE
   wavefreq(1:ncloc,1:nrloc) = 0.0
END WHERE

!
!5. Exchange halos
!-----------------
!

IF (iopt_MPI.EQ.1) THEN

!  ---wave parameters
   lbounds = 0; nhexch = (/1,0,1,0/)
   CALL exchange_mod(waveheight,lbounds,nhexch,iarr_waveheight)
   CALL exchange_mod(waveperiod,lbounds,nhexch,iarr_waveperiod)
   CALL exchange_mod(wavefreq,lbounds,nhexch,iarr_wavefreq)
   CALL exchange_mod(wavedir,lbounds,nhexch,iarr_wavedir)
   IF (iopt_waves_form.EQ.2) THEN
      CALL exchange_mod(wavevel,lbounds,nhexch,iarr_wavevel)
      CALL exchange_mod(waveexcurs,lbounds,nhexch,iarr_waveexcurs)
   ENDIF

!  ---Stokes drift and wave induced pressure
   IF (iopt_waves_curr.EQ.1.AND.iopt_waves_form.EQ.2) THEN
      lbounds = (/1-nhalo,0/); nhexch = (/nhalo,nhalo,1,0/)
      CALL exchange_mod(umstokesatc,lbounds,nhexch,iarr_umstokesatc)
      lbounds = (/0,1-nhalo/); nhexch = (/1,0,nhalo,nhalo/)
      CALL exchange_mod(vmstokesatc,lbounds,nhexch,iarr_vmstokesatc)
      !     ---wave pressure
      IF (iopt_waves_pres.EQ.1) THEN
         lbounds = 0; nhexch = (/1,0,1,0/)
         CALL exchange_mod(wavepres,lbounds,nhexch,iarr_wavepres)
      ENDIF
   ENDIF

!  ---dissipation forces
   IF (iopt_waves_dissip.EQ.1) THEN
      lbounds = (/0,1/); nhexch = (/1,0,0,0/) 
      CALL exchange_mod(umswdissipatc,lbounds,nhexch,iarr_umswdissipatc)
      CALL exchange_mod(umbwdissipatc,lbounds,nhexch,iarr_umbwdissipatc)
      lbounds = (/1,0/); nhexch = (/0,0,1,0/) 
      CALL exchange_mod(vmswdissipatc,lbounds,nhexch,iarr_vmswdissipatc)
      CALL exchange_mod(vmbwdissipatc,lbounds,nhexch,iarr_vmbwdissipatc)
   ENDIF
   
ENDIF

!
!6. Interpolate at velocity nodes
!--------------------------------
!

IF (iopt_waves_curr.EQ.1.AND.iopt_waves_form.EQ.2) THEN
   CALL Carr_at_U(umstokesatc(0:ncloc,1:nrloc),umstokesatu(1:ncloc,1:nrloc),&
                & 1,1,(/0,1,nz/),(/ncloc,nrloc,nz/),1,iarr_umstokesatc,info)
   CALL Carr_at_V(vmstokesatc(1:ncloc,0:nrloc),vmstokesatv(1:ncloc,1:nrloc),&
                & 1,1,(/1,0,nz/),(/ncloc,nrloc,nz/),1,iarr_vmstokesatc,info)
ENDIF

!
!7. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1.AND.iopt_waves_curr.EQ.1.AND.iopt_waves_form.EQ.2) THEN
   lbounds = (/1-nhalo,0/); nhexch = (/nhalo,nhalo,1,0/)
   CALL exchange_mod(umstokesatu,lbounds,nhexch,iarr_umstokesatu)
   lbounds = (/0,1-nhalo/); nhexch = (/1,0,nhalo,nhalo/)
   CALL exchange_mod(vmstokesatv,lbounds,nhexch,iarr_vmstokesatv)
ENDIF

!
!8. Update wave parameters
!-------------------------
!


IF (.NOT.cold_start.OR.iopt_waves_couple.EQ.0) CALL update_wave_params

!
!9. Deallocate arrays
!--------------------
!

IF (cold_start.OR.nt+ABS(modfiles(io_wavsur,1,1)%tlims(3)).GT.nstep) THEN
   IF (iopt_waves_couple.EQ.0) THEN
      DEALLOCATE (angle,wavegrid,wavedat,wavevals)
   ELSEIF (modfiles(io_wavsur,1,2)%status.EQ.'W') THEN
      DEALLOCATE (wavevals,wavevalsglb)
   ENDIF
ENDIF

1000 CALL log_timer_out(npcc,itm_waves)


RETURN

END SUBROUTINE wave_input

!========================================================================

SUBROUTINE wave_sources_momentum(sourceatu,sourceatv,nzdim)
!************************************************************************
!
! *wave_sources_momentum* calculate the wave source terms in the momentum
!                         equations
!
! Author - Alexander Breugem and Patrick Luyten (IMDC)
!
! Version - @(COHERENS)Surface_Waves.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program -
!
! External calls -
!
! Module calls - Carr_at_U, Carr_at_V, Uarr_at_C, Varr_at_C, error_alloc
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Carr_at_U, Carr_at_V, Uarr_at_C, Varr_at_C 
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

INTEGER, INTENT(IN) :: nzdim
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nzdim) :: sourceatu, sourceatv

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*nzdim*      INTEGER Vertical dimension (1 for 2-D, nz for 3-D case)
!*sourceatu*  REAL    Source term in the U-momentum equation    [m/s] or [m^2/s]
!*sourceatv*  REAL    Source term in the V-momentum equation    [m/s] or [m^2/s]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: info = .TRUE.
INTEGER :: k, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu, maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: umbwdissipatu, umstokesatv, &
        & umswdissipatu, umvelatc, vmbwdissipatv, vmstokesatu, vmswdissipatv, &
        & vmvelatc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3d, ubwdissipatu, &
        & ustokesatv, uswdissipatu, uvelatc, vbwdissipatv, vstokesatu, &
        & vswdissipatv, vvelatc


procname(pglev+1) = 'wave_sources_momentum'
CALL log_timer_in(npcc)

!
!1. Initialise
!--------------
!
!1.1 Allocate
!------------
!

ALLOCATE (maskatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (maskatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc/),kndlog)

!
!1.2 Mask arrays
!---------------
!

maskatu = node2du(1:ncloc,1:nrloc).EQ.2
maskatv = node2dv(1:ncloc,1:nrloc).EQ.2

!
!2. 3-D momentum equations
!-------------------------
!

IF (nzdim.GT.1) THEN

!
!2.1 Allocate
!------------
!

 
   ALLOCATE (array3d(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('array3d',3,(/ncloc,nrloc,nz/),kndrtype)
   ALLOCATE (uvelatc(0:ncloc,0:nrloc,nz),STAT=errstat)
   CALL error_alloc('uvelatc',3,(/ncloc+1,nrloc+1,nz/),kndrtype)
   ALLOCATE (vvelatc(0:ncloc,0:nrloc,nz),STAT=errstat)
   CALL error_alloc('vvelatc',3,(/ncloc+1,nrloc+1,nz/),kndrtype)
   ALLOCATE (vstokesatu(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('vstokesatu',3,(/ncloc,nrloc,nz/),kndrtype)
   ALLOCATE (ustokesatv(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('ustokesatv',3,(/ncloc,nrloc,nz/),kndrtype)
   IF (iopt_waves_dissip.EQ.1) THEN
      ALLOCATE (uswdissipatu(ncloc,nrloc,nz),STAT=errstat)
      CALL error_alloc('uswdissipatu',3,(/ncloc,nrloc,nz/),kndrtype)
      ALLOCATE (ubwdissipatu(ncloc,nrloc,nz),STAT=errstat)
      CALL error_alloc('ubwdissipatu',3,(/ncloc,nrloc,nz/),kndrtype)
      ALLOCATE (vbwdissipatv(ncloc,nrloc,nz),STAT=errstat)
      CALL error_alloc('vbwdissipatv',3,(/ncloc,nrloc,nz/),kndrtype)
      ALLOCATE (vswdissipatv(ncloc,nrloc,nz),STAT=errstat)
      CALL error_alloc('vswdissipatv',3,(/ncloc,nrloc,nz/),kndrtype)
   ENDIF

!
!2.2 Interpolate
!---------------
!
!  ---currents
   CALL Uarr_at_C(uvel(0:ncloc+1,0:nrloc,:),uvelatc,1,1,(/0,0,1/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_uvel,info)
   CALL Varr_at_C(vvel(0:ncloc,0:nrloc+1,:),vvelatc,1,1,(/0,0,1/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vvel,info)
!  ---Stokes drift
   CALL Carr_at_U(vstokesatc(:,1:nrloc,:),vstokesatu,1,1,(/0,1,1/),&
               & (/ncloc,nrloc,nz/),1,iarr_vstokesatc,info)
   CALL Carr_at_V(ustokesatc(1:ncloc,:,:),ustokesatv,1,1,(/1,0,1/),&
               & (/ncloc,nrloc,nz/),1,iarr_ustokesatc,info)
!  ---dissipation forces
   IF (iopt_waves_dissip.EQ.1) THEN
      CALL Carr_at_U(uswdissipatc,uswdissipatu,1,1,(/0,1,1/),&
                  & (/ncloc,nrloc,nz/),1,iarr_uswdissipatc,info)
      CALL Carr_at_U(ubwdissipatc,ubwdissipatu,1,1,(/0,1,1/),&
                  & (/ncloc,nrloc,nz/),1,iarr_ubwdissipatc,info)
      CALL Carr_at_V(vswdissipatc,vswdissipatv,1,1,(/1,0,1/),&
                  & (/ncloc,nrloc,nz/),1,iarr_vswdissipatc,info)
      CALL Carr_at_V(vbwdissipatc,vbwdissipatv,1,1,(/1,0,1/),&
                  & (/ncloc,nrloc,nz/),1,iarr_vbwdissipatc,info)
   ENDIF

!
!2.3 U-nodes
!-----------
!
 
   k_231: DO k=1,nz

      WHERE (maskatu)
!       ---Coriolis
        sourceatu(:,:,k) = sourceatu(:,:,k) + &
                         & coriolatu(1:ncloc,1:nrloc)*vstokesatu(:,:,k)
!        ---Stokes terms
         array3d(:,:,k) = (ustokesatu(1:ncloc,1:nrloc,k)*&
                 & (uvelatc(1:ncloc,1:nrloc,k)-uvelatc(0:ncloc-1,1:nrloc,k))+&
                 & vstokesatu(:,:,k)*&
                 & (vvelatc(1:ncloc,1:nrloc,k)-vvelatc(0:ncloc-1,1:nrloc,k)))/&
                 & delxatu(1:ncloc,1:nrloc)
         sourceatu(:,:,k) = sourceatu(:,:,k) + array3d(:,:,k)
      END WHERE

!     ---pressure force
      IF (iopt_waves_pres.EQ.1) THEN
         WHERE (maskatu)
            sourceatu(:,:,k) = sourceatu(:,:,k) + &
                 & (wavepres(0:ncloc-1,1:nrloc)-wavepres(1:ncloc,1:nrloc))/&
                 & delxatu(1:ncloc,1:nrloc)
         END WHERE
      ENDIF
      
!     ---dissipation forces
      IF (iopt_waves_dissip.EQ.1) THEN
         WHERE (maskatu)
            sourceatu(:,:,k) = sourceatu(:,:,k) +  uswdissipatu(:,:,k) + &
                             & ubwdissipatu(:,:,k)
            
         END WHERE
      ENDIF

   ENDDO k_231

!  ---integrated Stokes term
   IF (iopt_hydro_impl.EQ.0) THEN
      stokessource2du = 0.0
      k_232: DO k=1,nz
         WHERE (maskatu)
            stokessource2du = stokessource2du + &
                            & delzatu(1:ncloc,:,k)*array3d(:,:,k)
         END WHERE
      ENDDO k_232
   ENDIF

!
!2.4 V-nodes
!-----------
!

   k_241: DO k=1,nz

      WHERE (maskatv)
!       ---Coriolis
        sourceatv(:,:,k) = sourceatv(:,:,k) - &
                        & coriolatv(1:ncloc,1:nrloc)*ustokesatv(:,:,k)
!        ---Stokes terms
         array3d(:,:,k) = (ustokesatv(:,:,k)*&
                 & (uvelatc(1:ncloc,1:nrloc,k)-uvelatc(1:ncloc,0:nrloc-1,k))+&
                 & vstokesatv(1:ncloc,1:nrloc,k)*&
                 & (vvelatc(1:ncloc,1:nrloc,k)-vvelatc(1:ncloc,0:nrloc-1,k)))/&
                 & delyatv(1:ncloc,1:nrloc)
         sourceatv(:,:,k) = sourceatv(:,:,k) + array3d(:,:,k)
      END WHERE

!     ---pressure force
      IF (iopt_waves_pres.EQ.1) THEN
         WHERE (maskatv)
            sourceatv(:,:,k) = sourceatv(:,:,k) + &
                 & (wavepres(1:ncloc,0:nrloc-1)-wavepres(1:ncloc,1:nrloc))/&
                 & delyatv(1:ncloc,1:nrloc)
         END WHERE
      ENDIF
      
!     ---dissipation forces
      IF (iopt_waves_dissip.EQ.1) THEN
         WHERE (maskatv)
            sourceatv(:,:,k) = sourceatv(:,:,k) +  vswdissipatv(:,:,k) + &
                             & vbwdissipatv(:,:,k)
         END WHERE
      ENDIF

   ENDDO k_241

!  ---integrated Stokes term
   IF (iopt_hydro_impl.EQ.0) THEN
      stokessource2dv = 0.0
      k_242: DO k=1,nz
         WHERE (maskatv)
            stokessource2dv = stokessource2dv + &
                            & delzatv(:,1:nrloc,k)*array3d(:,:,k)
         END WHERE
      ENDDO k_242
   ENDIF

!
!2.5 Deallocate
!--------------
!
   
   DEALLOCATE (array3d,uvelatc,vvelatc,vstokesatu,ustokesatv)
   IF (iopt_waves_dissip.EQ.1) THEN
      DEALLOCATE(uswdissipatu,ubwdissipatu,vswdissipatv,vbwdissipatv)
   ENDIF
   
!
!3. 2-D momentum equations
!-------------------------
!

ELSEIF (nzdim.EQ.1) THEN
 
!
!3.1 Allocate
!------------
!

   ALLOCATE (vmstokesatu(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('vmstokesatu',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (umstokesatv(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('umstokesatv',2,(/ncloc,nrloc/),kndrtype)
   IF (iopt_waves_dissip.EQ.1) THEN
      ALLOCATE (umswdissipatu(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('umswdissipatu',2,(/ncloc,nrloc/),kndrtype)
      ALLOCATE (umbwdissipatu(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('umbwdissipatu',2,(/ncloc,nrloc/),kndrtype)
      ALLOCATE (vmswdissipatv(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('vmswdissipatv',2,(/ncloc,nrloc/),kndrtype)
      ALLOCATE (vmbwdissipatv(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('vmbwdissipatv',2,(/ncloc,nrloc/),kndrtype)
   ENDIF
   IF (iopt_grid_nodim.EQ.2) THEN
      ALLOCATE (umvelatc(0:ncloc,0:nrloc),STAT=errstat)
      CALL error_alloc('umvelatc',2,(/ncloc+1,nrloc+1/),kndrtype)
      ALLOCATE (vmvelatc(0:ncloc,0:nrloc),STAT=errstat)
      CALL error_alloc('vmvelatc',2,(/ncloc+1,nrloc+1/),kndrtype)
   ENDIF

!
!3.2 Interpolate at velocity nodes
!---------------------------------
!
!  ---currents
   IF (iopt_grid_nodim.EQ.2) THEN
      CALL Uarr_at_C(umvel(0:ncloc+1,0:nrloc),umvelatc,1,1,(/0,0,nz/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_umvel,info)
      CALL Varr_at_C(vmvel(0:ncloc,0:nrloc+1),vmvelatc,1,1,(/0,0,nz/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vmvel,info)
   ENDIF
 
!  ---Stokes drift
   CALL Carr_at_U(vmstokesatc(:,1:nrloc),vmstokesatu,1,1,(/0,1,nz/),&
               & (/ncloc,nrloc,nz/),1,iarr_vmstokesatc,info)
   CALL Carr_at_V(umstokesatc(1:ncloc,:),umstokesatv,1,1,(/1,0,nz/),&
               & (/ncloc,nrloc,nz/),1,iarr_umstokesatc,info)
 
!  ---dissipation forces
   IF (iopt_waves_dissip.EQ.1) THEN
      CALL Carr_at_U(umswdissipatc,umswdissipatu,1,1,(/0,1,nz/),&
               & (/ncloc,nrloc,nz/),1,iarr_umswdissipatc,info)
      CALL Carr_at_U(umbwdissipatc,umbwdissipatu,1,1,(/0,1,nz/),&
               & (/ncloc,nrloc,nz/),1,iarr_umbwdissipatc,info)
      CALL Carr_at_V(vmswdissipatc,vmswdissipatv,1,1,(/1,0,nz/),&
               & (/ncloc,nrloc,nz/),1,iarr_vmswdissipatc,info)
      CALL Carr_at_V(vmbwdissipatc,vmbwdissipatv,1,1,(/1,0,nz/),&
               & (/ncloc,nrloc,nz/),1,iarr_vmbwdissipatc,info)
   ENDIF

!
!3.3 U-nodes
!-----------
!
!  ---Coriolis
   WHERE (maskatu)
      sourceatu(:,:,1) = sourceatu(:,:,1) + coriolatu(1:ncloc,1:nrloc)*&
                                       & deptotatu(1:ncloc,1:nrloc)*vmstokesatu
   END WHERE
   
!  ---Stokes terms
   IF (iopt_grid_nodim.EQ.3.AND.iopt_hydro_impl.EQ.0) THEN
      WHERE (maskatu)
         sourceatu(:,:,1) = sourceatu(:,:,1) + stokessource2du
      END WHERE
   ELSE
      WHERE (maskatu)
         sourceatu(:,:,1) = sourceatu (:,:,1) + &
                 & deptotatu(1:ncloc,1:nrloc)*(umstokesatu(1:ncloc,1:nrloc)*&
                 & (umvelatc(1:ncloc,1:nrloc)-umvelatc(0:ncloc-1,1:nrloc))+&
                 & vmstokesatu*&
                 & (vmvelatc(1:ncloc,1:nrloc)-vmvelatc(0:ncloc-1,1:nrloc)))/&
                 & delxatu(1:ncloc,1:nrloc)
      END WHERE
   ENDIF

!  ---pressure force
   IF (iopt_waves_pres.EQ.1) THEN
      WHERE (maskatu)
         sourceatu(:,:,1) = sourceatu(:,:,1) + deptotatu(1:ncloc,1:nrloc)*&
                    & (wavepres(0:ncloc-1,1:nrloc)-wavepres(1:ncloc,1:nrloc))/&
                    & delxatu(1:ncloc,1:nrloc)
      END WHERE
   ENDIF
   
!  ---dissipation forces
   IF (iopt_waves_dissip.EQ.1) THEN
      WHERE (maskatu)
         sourceatu(:,:,1) = sourceatu(:,:,1) +  umswdissipatu + umbwdissipatu
      END WHERE
   ENDIF

!
!3.4 V-nodes
!-----------
!
!  ---Coriolis
   WHERE (maskatv)
      sourceatv(:,:,1) = sourceatv(:,:,1) - coriolatv(1:ncloc,1:nrloc)*&
                                       & deptotatv(1:ncloc,1:nrloc)*umstokesatv
   END WHERE

!  ---Stokes terms
   IF (iopt_grid_nodim.EQ.3.AND.iopt_hydro_impl.EQ.0) THEN
      WHERE (maskatv)
         sourceatv(:,:,1) = sourceatv(:,:,1) + stokessource2dv
      END WHERE
   ELSE
      WHERE (maskatv)
         sourceatv(:,:,1) = sourceatv(:,:,1) + &
                & deptotatv(1:ncloc,1:nrloc)*(umstokesatv*&
                & (umvelatc(1:ncloc,1:nrloc)-umvelatc(1:ncloc,0:nrloc-1))+&
                & vmstokesatv(1:ncloc,1:nrloc)*&
                & (vmvelatc(1:ncloc,1:nrloc)-vmvelatc(1:ncloc,0:nrloc-1)))/&
                & delyatv(1:ncloc,1:nrloc)
      END WHERE
   ENDIF

!  ---pressure force
   IF (iopt_waves_pres.EQ.1) THEN
      WHERE (maskatv)
         sourceatv(:,:,1) = sourceatv(:,:,1) + deptotatv(1:ncloc,1:nrloc)*&
                    & (wavepres(1:ncloc,0:nrloc-1)-wavepres(1:ncloc,1:nrloc))/&
                    & delyatv(1:ncloc,1:nrloc)
      END WHERE
   ENDIF
   
!  ---dissipation forces
   IF (iopt_waves_dissip.EQ.1) THEN
      WHERE (maskatv)
         sourceatv(:,:,1) = sourceatv(:,:,1) +  vmswdissipatv + vmbwdissipatv
      END WHERE
   ENDIF

!
!3.5 Deallocate
!--------------
!

   IF (iopt_grid_nodim.EQ.2) DEALLOCATE (umvelatc,vmvelatc)
   DEALLOCATE (vmstokesatu,umstokesatv)
   IF (iopt_waves_dissip.EQ.1) THEN
      DEALLOCATE (umswdissipatu,umbwdissipatu,vmswdissipatv,vmbwdissipatv)
   ENDIF

ENDIF

!
!4. Deallocate
!-------------
!

DEALLOCATE (maskatu,maskatv)

CALL log_timer_out(npcc,itm_waves)


RETURN

END SUBROUTINE wave_sources_momentum
