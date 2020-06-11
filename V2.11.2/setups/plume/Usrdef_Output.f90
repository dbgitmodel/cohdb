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

SUBROUTINE usrdef_output
!************************************************************************
!
! *usrdef_output* User-formatted output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.8
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description - output parameters for test case plume
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_mod, copy_vars, mult_index, num_proc,
!                open_file
!
!************************************************************************
!
USE currents
USE density
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE obconds
USE paralpars
USE switches
USE syspars
USE tide
USE timepars
USE grid_routines, ONLY: num_proc
USE inout_routines, ONLY: close_file, open_file
USE paral_comms, ONLY: combine_mod, copy_vars
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER :: i, ib, ic, icloc1, icloc2, ic1, ic2, idproc1, idproc2, ihour,&
         & j, jc, jcloc1, jcloc2, jc1, jc2, k, kc
INTEGER, SAVE :: icout, iunit, knt, mtot
REAL :: bellip, belmaj, belmin, bures, bvres, denom, dhun, duamp, dupha,&
      & dures, fcos, freq, fsin, gxc, gyc, hbulge, hfront, hwidth, pdep, &
      & phase, rmin, rmini, rminr, rplus, rplusi, rplusr, scrit, sellip, &
      & selmaj, selmin, sumc, sumcc, sumss, sures, svres, umean, zetamp, &
      & zetpha, zetres
REAL :: rflag = -999.9
REAL, DIMENSION(nc-1) :: width
REAL, SAVE, DIMENSION(0:2,6) :: sum
REAL, DIMENSION(0:nc+1,0:nr+1) :: deptotglb
REAL, DIMENSION(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo,nz) :: salglb

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*first_call* LOGICAL Flag to determine if this is the first call
!*suffix*     CHAR    Suffix of output file name
!*ihour*      INTEGER Number of hours since start of simulation
!*iunit*      INTEGER Unit number for output file
!*knt*        INTEGER Time counter for harmonic analysis
!*sum*        REAL    Local arrays of sums for harmonic analysis
!*width*      REAL    Plume width                                          [km]
!*bellip*     REAL    Ellipticity of tidal ellipse at P and 10m depth
!*belmaj*     REAL    Semi-major axis of tidal ellipse at P and 10m depth [m/s]
!*belmin*     REAL    Semi-minor axis of tidal ellipse at P and 10m depth [m/s]
!*bures*      REAL    Residual u-current at P and 10m depth               [m/s]
!*bvres*      REAL    Residual v-current at P and 10m depth               [m/s]
!*duamp*      REAL    Amplitude of depth averaged current at Q            [m/s]
!*dupha*      REAL    Phase of depth averaged current at Q                [deg]
!*dures*      REAL    Residual value of depth averaged current at Q       [m/s]
!*hbulge*     REAL    Width of the plume bulge                             [km]
!*hfront*     REAL    Plume length                                         [km]
!*hwidth*     REAL    Plume width along BB transect                        [km]
!*pdep*       REAL    Plume depth at P                                      [m]
!*sellip*     REAL    Ellipticity of surface tidal ellipse at P
!*selmaj*     REAL    Semi-major axis of surface tidal ellipse at P       [m/s]
!*selmin*     REAL    Semi-minor axis of surface tidal ellipse at P       [m/s]
!*sures*      REAL    Residual surface u-current at P                     [m/s]
!*svres*      REAL    Residual surface v-current at P                     [m/s]
!*zetamp*     REAL    Amplitude of sea surface elevation at Q               [m]
!*zetpha*     REAL    Phase of sea surface elevation at Q                 [deg]
!*zetres*     REAL    Residual sea surface elevation at Q                   [m]
!
!------------------------------------------------------------------------------
!
!2. Initialise parameters
!------------------------
!

IF (iopt_sal.EQ.1) RETURN
procname(pglev+1) = 'usrdef_output'
CALL log_timer_in()

IF (nt.EQ.0) THEN
!  ---open output file
   IF (.NOT.cold_start.AND.master) THEN
      CALL open_file(iunit,TRIM(outtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case plume: simulation '&
                          &//TRIM(outtitle)
      WRITE (iunit,*)
   ENDIF
!  ---initialise parameters for harmonic analysis
   mtot = 72
   knt = -mtot
   sum  = 0.0
!  ---output frequency
   icout = MERGE(360,36,iopt_hydro_impl.EQ.0)
ENDIF

!
!3. Evaluate/write parameters
!----------------------------
!

ihour = MERGE(nt/120,nt/12,iopt_hydro_impl.EQ.0)

IF (nt.GT.0.AND.mult_index(nt,icout)) THEN

!
!3.1 Combine at master process
!-----------------------------
!

   CALL combine_mod(salglb,sal,(/1-nhalo,1-nhalo,1/),iarr_sal,0.0)
   CALL combine_mod(deptotglb,deptotatc,(/0,0/),iarr_deptotatc,0.0)

   IF (master) THEN

!
!3.2 Plume widths
!----------------
!

      dhun = 1000.0; scrit = 32.0
      i_321: DO i=1,nc-1
         jc = 0
         j_3211: DO j=1,nr-1
            IF ((jc.EQ.0).AND.(salglb(i,j,nz).GT.scrit)) jc = j
         ENDDO j_3211
         IF (jc.EQ.0) THEN
            width(i) = dhun*(nr-1)
         ELSEIF (jc.EQ.1) THEN
            width(i) = 0.0
         ELSE
            gyc = dhun*(jc-1.5)
            width(i) = gyc+(scrit-salglb(i,jc-1,nz))*dhun/&
                         & (salglb(i,jc,nz)-salglb(i,jc-1,nz))
         ENDIF
         width(i) = 0.001*width(i)
      ENDDO i_321

!     ---bulge width
      hbulge = MAXVAL(width)

!     ---width of coastal plume
      ic = 50
      hwidth = width(ic)

!
!3.3 Halocline depth
!-------------------
!

      scrit = 30.0
      ic = 50; jc = 5; kc = 0
      k_330: DO k=nz,1,-1
         IF ((kc.EQ.0).AND.(salglb(ic,jc,k).GT.scrit)) kc = k
      ENDDO k_330
      IF (kc.EQ.0) THEN
         pdep = deptotglb(ic,jc)
      ELSEIF  (kc.EQ.nz) THEN
         pdep = 0.0
      ELSE
         pdep = (deptotglb(ic,jc)/REAL(nz))*(nz-kc+0.5-&
                  & (scrit-salglb(ic,jc,kc))/&
                  & (salglb(ic,jc,kc+1)-salglb(ic,jc,kc)))
      ENDIF

!
!3.4 Plume length
!----------------
!

      scrit = 32.0
      jc = 2; ic = 0; ib = iobv(nobv)
      i_340: DO i=ib,nc-1
         IF ((ic.EQ.0).AND.(salglb(i,jc,nz).GT.scrit)) ic = i
      ENDDO i_340
      IF (ic.EQ.0) THEN
         hfront = dhun*(nc-0.5-ib)
      ELSEIF (ic.EQ.ib) THEN
         hfront = 0.0
      ELSE
         gxc = dhun*(ic-ib-1)
         hfront = gxc+(scrit-salglb(ic-1,jc,nz))*dhun/&
                    & (salglb(ic,jc,nz)-salglb(ic-1,jc,nz))
      ENDIF
      hfront = 0.001*hfront

!
!3.5 Write output parameters
!---------------------------
!

      WRITE (iunit,9001) ihour
      WRITE (iunit,9002) 'hbulge', hbulge
      WRITE (iunit,9002) 'hwidth', hwidth
      WRITE (iunit,9002) 'pdep', pdep
      WRITE (iunit,9002) 'hfront', hfront

   ENDIF

ENDIF

!
!4. Harmonic analysis
!--------------------
!
!4.1 Update sums
!---------------
!

IF (corrstep.AND.(ihour.GE.60).AND.(ihour.LE.72)) THEN

   freq = tidal_spectrum(index_obc(1))
   ic1 = 50; jc1 = 5; kc = 10
   IF (nconobc.GT.0) THEN
      phase = freq*knt*delt3d
   ELSE
      phase = 0.0
   ENDIF
   idproc1 = num_proc(ic1,jc1) - 1

   IF (idloc.EQ.idproc1) THEN

      icloc1 = ic1 - nc1loc + 1; jcloc1 = jc1 - nr1loc + 1

!     ---surface
      sum(0,1) = sum(0,1) + uvel(icloc1,jcloc1,nz)
      sum(1,1) = sum(1,1) + uvel(icloc1,jcloc1,nz)*COS(phase)
      sum(2,1) = sum(2,1) + uvel(icloc1,jcloc1,nz)*SIN(phase)
      sum(0,2) = sum(0,2) + vvel(icloc1,jcloc1,nz)
      sum(1,2) = sum(1,2) + vvel(icloc1,jcloc1,nz)*COS(phase)
      sum(2,2) = sum(2,2) + vvel(icloc1,jcloc1,nz)*SIN(phase)

!     ---bottom
      sum(0,3) = sum(0,3) + uvel(icloc1,jcloc1,kc)
      sum(1,3) = sum(1,3) + uvel(icloc1,jcloc1,kc)*COS(phase)
      sum(2,3) = sum(2,3) + uvel(icloc1,jcloc1,kc)*SIN(phase)
      sum(0,4) = sum(0,4) + vvel(icloc1,jcloc1,kc)
      sum(1,4) = sum(1,4) + vvel(icloc1,jcloc1,kc)*COS(phase)
      sum(2,4) = sum(2,4) + vvel(icloc1,jcloc1,kc)*SIN(phase)
      
   ENDIF

!  ---depth averaged
   ic2 = 50; jc2 = 30
   idproc2 = num_proc(ic2,jc2) - 1
   IF (idloc.EQ.idproc2) THEN
      icloc2 = ic2 - nc1loc + 1; jcloc2 = jc2 - nr1loc + 1
      umean = udvel(icloc2,jcloc2)/deptotatc(icloc2,jcloc2)
      sum(0,5) = sum(0,5) + umean
      sum(1,5) = sum(1,5) + umean*COS(phase)
      sum(2,5) = sum(2,5) + umean*SIN(phase)
   ENDIF

!  ---surface elevation
   IF (idloc.EQ.idproc2) THEN
      sum(0,6) = sum(0,6) + zeta(icloc2,jcloc2)
      sum(1,6) = sum(1,6) + zeta(icloc2,jcloc2)*COS(phase)
      sum(2,6) = sum(2,6) + zeta(icloc2,jcloc2)*SIN(phase)
   ENDIF
   
   knt = knt + 1

ENDIF

!
!4.2 Evaluate harmonic parameters
!--------------------------------
!
 
IF (ihour.EQ.72) THEN

!  ---with tidal forcing

   IF (nconobc.GT.0) THEN

!     ---harmonic sums
      sumc = SIN((mtot+0.5)*freq*delt3d)/SIN(0.5*freq*delt3d)
      sumcc = 0.5*SIN((2*mtot+1)*freq*delt3d)/&
             & SIN(freq*delt3d) + mtot + 0.5
      sumss = 2*mtot + 1 - sumcc
      denom = (2*mtot+1)*sumcc - sumc*sumc

!     ---surface
      IF (idloc.EQ.idproc1) THEN
         sures = (sumcc*sum(0,1)-sumc*sum(1,1))/denom
         svres = (sumcc*sum(0,2)-sumc*sum(1,2))/denom
         rplusr = 0.5*(sum(1,1)+sum(2,2))
         rplusi = 0.5*(sum(1,2)-sum(2,1))
         rminr  = 0.5*(sum(1,1)-sum(2,2))
         rmini  =-0.5*(sum(1,2)+sum(2,1))
         rplus = SQRT(rplusr*rplusr + rplusi*rplusi)
         rmin =  SQRT(rminr*rminr   + rmini*rmini)
         selmaj = rplus + rmin
         selmin = rplus - rmin
         IF (selmaj.GT.1.0E-20) THEN
            sellip = selmin/selmaj
         ELSE
            sellip = 0.0
         ENDIF
      ENDIF

!     ---bottom
      IF (idloc.EQ.idproc1) THEN
         bures = (sumcc*sum(0,3)-sumc*sum(1,3))/denom
         bvres = (sumcc*sum(0,4)-sumc*sum(1,4))/denom
         rplusr = 0.5*(sum(1,3)+sum(2,4))
         rplusi = 0.5*(sum(1,4)-sum(2,3))
         rminr  = 0.5*(sum(1,3)-sum(2,4))
         rmini  =-0.5*(sum(1,4)+sum(2,3))
         rplus = SQRT(rplusr*rplusr + rplusi*rplusi)
         rmin =  SQRT(rminr*rminr   + rmini*rmini)
         belmaj = rplus + rmin
         belmin = rplus - rmin
         IF (belmaj.GT.1.0E-20) THEN
            bellip = belmin/belmaj
         ELSE
            bellip = 0.0
         ENDIF
      ENDIF

!     ---depth averaged
      IF (idloc.EQ.idproc2) THEN
         dures = (sumcc*sum(0,5)-sumc*sum(1,5))/denom
         fcos = ((2*mtot+1)*sum(1,5)-sumc*sum(0,5))/denom
         fsin = sum(2,5)/sumss
         duamp = SQRT(fcos*fcos+fsin*fsin)
         dupha = ATAN2(fsin,fcos)
         IF (dupha.LT.0.0) dupha = dupha + twopi
         dupha = radtodeg*dupha
      ENDIF

!     ---surface elevation
      IF (idloc.EQ.idproc2) THEN
         zetres = (sumcc*sum(0,6)-sumc*sum(1,6))/denom
         fcos = ((2*mtot+1)*sum(1,6)-sumc*sum(0,6))/denom
         fsin = sum(2,6)/sumss
         zetamp = SQRT(fcos*fcos+fsin*fsin)
         zetpha = ATAN2(fsin,fcos)
         IF (zetpha.LT.0.0) zetpha = zetpha + twopi
         zetpha = radtodeg*zetpha
      ENDIF

!  ---without tidal forcing

   ELSE
!     ---surface
      IF (idloc.EQ.idproc1) THEN
         sures = sum(0,1)/REAL(mtot)
         svres = sum(0,2)/REAL(mtot)
         selmaj = rflag
         sellip = rflag
      ENDIF
!     ---bottom
      IF (idloc.EQ.idproc1) THEN
         bures = sum(0,3)/REAL(mtot)
         bvres = sum(0,4)/REAL(mtot)
         belmaj = rflag
         bellip = rflag
      ENDIF
!     ---depth averaged
      IF (idloc.EQ.idproc2) THEN
         dures = sum(0,5)/REAL(mtot)
         duamp = rflag
         dupha = rflag
      ENDIF
!     ---surface elevation
      IF (idloc.EQ.idproc2) THEN
         zetres = sum(0,6)/REAL(mtot)
         zetamp = rflag
         zetpha = rflag
      ENDIF

   ENDIF

ENDIF

!
!4.3 Write output parameters
!---------------------------
!

IF (ihour.EQ.72) THEN

!
!4.3.1 Copy parameters to root process
!-------------------------------------
!

   IF (iopt_MPI.EQ.1) THEN
!     ---surface
      CALL copy_vars(sures,0,idroot=idproc1)
      CALL copy_vars(svres,0,idroot=idproc1)
      CALL copy_vars(selmaj,0,idroot=idproc1)
      CALL copy_vars(sellip,0,idroot=idproc1)
!     ---bottom
      CALL copy_vars(bures,0,idroot=idproc1)
      CALL copy_vars(bvres,0,idroot=idproc1)
      CALL copy_vars(belmaj,0,idroot=idproc1)
      CALL copy_vars(bellip,0,idroot=idproc1)
!     ---depth averaged
      CALL copy_vars(dures,0,idroot=idproc2)
      CALL copy_vars(duamp,0,idroot=idproc2)
      CALL copy_vars(dupha,0,idroot=idproc2)
!     ---surface elevation
      CALL copy_vars(zetres,0,idroot=idproc2)
      CALL copy_vars(zetamp,0,idroot=idproc2)
      CALL copy_vars(zetpha,0,idroot=idproc2)
   ENDIF

!
!4.3.2 Write parameters
!----------------------
!

   IF (master) THEN
      WRITE (iunit,*)
      WRITE (iunit,'(A)') 'Harmonic parameters'
!     ---surface
      WRITE (iunit,9002) 'sures', sures
      WRITE (iunit,9002) 'svres', svres
      WRITE (iunit,9002) 'selmaj', selmaj
      WRITE (iunit,9002) 'sellip', sellip
!     ---bottom
      WRITE (iunit,9002) 'bures', bures
      WRITE (iunit,9002) 'bvres', bvres
      WRITE (iunit,9002) 'belmaj', belmaj
      WRITE (iunit,9002) 'bellip', bellip
!     ---depth averaged
      WRITE (iunit,9002) 'dures', dures
      WRITE (iunit,9002) 'duamp', duamp
      WRITE (iunit,9002) 'dupha', dupha
!     ---surface elevation
      WRITE (iunit,9002) 'zetres', zetres
      WRITE (iunit,9002) 'zetamp', zetamp
      WRITE (iunit,9002) 'zetpha', zetpha
   ENDIF

ENDIF

!
!5. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')

CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',I2,' hours')
9002 FORMAT (T3,A,T12,': ',G12.4E3)

END SUBROUTINE usrdef_output
