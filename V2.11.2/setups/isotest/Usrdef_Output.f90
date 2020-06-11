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
! Version - @(COHERENS)Usrdef_Output.f90  V2.11.2
!
! $Date: 2018-07-23 16:55:25 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1171 $
!
! Description - output parameters for test case isotest
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
USE density
USE depths
USE diffusion
USE grid
USE gridpars
USE iopars
USE paralpars
USE physpars
USE switches
USE timepars
USE inout_routines, ONLY: close_file, open_file
USE paral_utilities, ONLY: max_vars, min_vars, sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out  
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL zflag
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER, SAVE :: icout, iunit
INTEGER :: i, iday, iip, ip, ivar, k, kkp, kp, l
INTEGER, DIMENSION(4) :: nhdims
REAL :: psisum
REAL, DIMENSION(3) :: fldiffprojmax, fldiffprojmean, fldiffprojmin, psivar, &
                    & xdiffmax, xdiffmean, xdiffmin, zdiffmax, zdiffmean, &
                    & zdiffmin
REAL, DIMENSION(2,nz) :: xdiff
REAL, DIMENSION(nz+1) :: zdiff, zgrad
REAL, DIMENSION(ncloc+1,nz) :: xgrad
REAL, DIMENSION(ncloc,nrloc) :: array2d
REAL, DIMENSION(ncloc,nrloc,nz) :: array3d
REAL, DIMENSION(0:ncloc+1,nz,3) :: xdifflux
REAL, DIMENSION(0:ncloc+1,nz,3) :: psic
REAL, DIMENSION(ncloc,nz,3) :: fldiffproj
REAL, DIMENSION(ncloc,nz+1,3) :: zdifflux
REAL, DIMENSION(ncloc+1,nz,0:1,0:1) :: slopefac


!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,1)%status = 'R'
   modfiles(io_modgrd,1,1)%form = modfiles(io_modgrd,1,2)%form 
   modfiles(io_modgrd,1,1)%filename = modfiles(io_modgrd,1,2)%filename 
   modfiles(io_modgrd,1,2)%status = '0'
!  ---initial conditions
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = modfiles(io_inicon,ics_phys,2)%form
   modfiles(io_inicon,ics_phys,1)%filename = &
 & modfiles(io_fincon,ics_phys,2)%filename
   modfiles(io_fincon,ics_phys,2)%status = '0'
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN
!  ---open output file
   IF (master) THEN
      CALL open_file(iunit,TRIM(outtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case isotest: '//&
                        & 'simulation '//TRIM(outtitle)
      WRITE (iunit,*)
   ENDIF
!  ---output frequency
   icout = 864000/NINT(delt2d)
ENDIF

IF (.NOT.mult_index(nt,icout)) RETURN
zflag = runtitle(8:8).GT.'C'

procname(pglev+1) = 'usrdef_output'
CALL log_timer_in()

!
!3. Evaluate parameters
!----------------------
!

psic(:,:,1) = temp(0:ncloc+1,1,:)
psic(:,:,2) = sal(0:ncloc+1,1,:)
psic(:,:,3) = dens(0:ncloc+1,1,:) - 1000.0
xdifflux = 0.0; zdifflux = 0.0

!
!3.1 Diffusive fluxes of T and S
!------------------------------
!

ivar_310: DO ivar=1,2

!      
!3.1.1 Laplacian diffusion
!-------------------------
!
      
   IF (iopt_hdif_scal.EQ.1) THEN
      i_311: DO i=1,ncloc+1
         IF (node2du(i,1).EQ.2) THEN
            k_3111: DO k=1,nz
               xdifflux(i,k,ivar) = -hdifscal_fac*hdifcoef3datu(i,1,k)*&
                                 & (psic(i,k,ivar)-psic(i-1,k,ivar))&
                                 & /delxatu(i,1)
            ENDDO k_3111
         ENDIF
      ENDDO i_311

!         
!3.1.2 Isolevel diffusion
!------------------------
!         
         
   ELSEIF (iopt_hdif_scal.EQ.2) THEN
      i_3121: DO i=1,ncloc+1
         
!        ---X-component               
         IF (node2du(i,1).EQ.2) THEN
            k_31211: DO k=2,nz
               zdiff(k) = 0.5*(psic(i-1,k,ivar)+psic(i,k,ivar)-&
                             & psic(i-1,k-1,ivar)-psic(i,k-1,ivar))&
                             & /delzatuw(i,1,k)
            ENDDO k_31211
            zdiff(1) = zdiff(2)
            zdiff(nz+1) = zdiff(nz)
            k_31212: DO k=2,nz            
               xdifflux(i,k,ivar) = -hdifscal_fac*hdifcoef3datu(i,1,k)*&
                                  & (psic(i,k,ivar)-psic(i-1,k,ivar))&
                                  & /delxatu(i,1)
               xdifflux(i,k,ivar) = xdifflux(i,k,ivar) + xslopeatu_geo(i,1,k)*&
                                  & 0.5*(zdiff(k)+zdiff(k+1))
            ENDDO k_31212
         ENDIF

      ENDDO i_3121
         
!     ---Z-component
      i_3122: DO i=1,ncloc
         IF (maskatc_int(i,1).AND.zflag) THEN
            xdiff(:,1) = 0.0
            l_31221: DO l=1,2
            k_31221: DO k=2,nz
               xdiff(l,k) = 0.5*(psic(i+l-1,k-1,ivar)+psic(i+l-1,k,ivar)-&
                               & psic(i+l-2,k-1,ivar)-psic(i+l-2,k,ivar))&
                               & /delxatu(i+l-1,1)
            ENDDO k_31221
            ENDDO l_31221
            k_31222: DO k=2,nz
               zdifflux(i,k,ivar) = -vdifcoefscal_rot(i,1,k)*&
                                  & (psic(i,k,ivar)-psic(i,k-1,ivar))&
                                  & /delzatw(i,1,k)
               zdifflux(i,k,ivar) = zdifflux(i,k,ivar) + &
                                  & hdifscal_fac*hdifcoef3datw(i,1,k)*&
                                  & xslopeatw_geo(i,1,k)*&
                                  &  0.5*(xdiff(1,k)+xdiff(2,k))
             ENDDO k_31222
         ENDIF

      ENDDO i_3122
         
!         
!3.1.3 Isoneutral diffusion
!--------------------------
!
         
   ELSEIF (iopt_hdif_scal.EQ.3) THEN

      i_3131: DO i=1,ncloc+1
         
!        ---X-component                  
         IF (node2du(i,1).EQ.2) THEN
            kp_3131: DO kp=0,1
            ip_3131: DO ip=0,1
               iip = i+ip-1
               k_31311: DO k=1,nz
                  kkp = k+1-kp
                  IF (kkp.GT.1.AND.kkp.LE.nz) THEN
                     slopefac(i,k,ip,kp) = 0.25*hdifscal_fac*&
                                         & (hdifcoef3datu(i,1,k)*&
                                         & xslopeatu_siso(iip,1,k,1-ip,kp))&
                                         & /delzatw(iip,1,kkp)
                  ENDIF
               ENDDO k_31311
               k_31312: DO k=1,nz
                  kkp = k+1-kp
                  IF (kkp.GT.1.AND.kkp.LE.nz) THEN
                     zdiff(k) = slopefac(i,k,ip,kp)*&
                             & (psic(iip,kkp,ivar)-psic(iip,kkp-1,ivar))
                     xdifflux(i,k,ivar) = xdifflux(i,k,ivar) + zdiff(k)
                  ENDIF
               ENDDO k_31312
               IF (kkp.EQ.nz+1) THEN
                  xdifflux(i,nz,ivar) = xdifflux(i,nz,ivar) + zdiff(nz-1)
               ELSEIF (kkp.EQ.1) THEN
                  xdifflux(i,1,ivar) = xdifflux(i,1,ivar) + zdiff(2)
               ENDIF
            ENDDO ip_3131
            ENDDO kp_3131
            k_3132: DO k=1,nz
               xdifflux(i,k,ivar) = xdifflux(i,k,ivar) - &
                                  & hdifscal_fac*hdifcoef3datu(i,1,k)*&
                                  & (psic(i,k,ivar)-psic(i-1,k,ivar))&
                                  & /delxatu(i,1)
            ENDDO k_3132
         ENDIF

      ENDDO i_3131
         
!     ---Z-component
      i_3132: DO i=1,ncloc
         IF (maskatc_int(i,1).AND.zflag) THEN
            kp_3133: DO kp=0,1
            ip_3133: DO ip=0,1
               iip = i+1-ip
               k_31331: DO k=2,nz
                  kkp = k+kp-1
                  slopefac(i,k-1,ip,kp) = -hdifscal_fac*hdifcoef3datu(iip,1,k)*&
                                        & xslopeatu_siso(i,1,kkp,1-ip,kp)
               ENDDO k_31331
               k_31332: DO k=2,nz
                  zdifflux(i,k,ivar) = zdifflux(i,k,ivar) + &
                                     & slopefac(i,k-1,ip,kp)*&
                                    & (psic(iip,kkp,ivar)-psic(iip-1,kkp,ivar))&
                                    & /delxatu(iip,1)           
               ENDDO k_31332
            ENDDO ip_3133
            ENDDO kp_3133
            k_3134: DO k=2,nz
               zdifflux(i,k,ivar) = zdifflux(i,k,ivar) - &
                                  & vdifcoefscal_rot(i,1,k)*&
                                  & (psic(i,k,ivar)-psic(i,k-1,ivar))&
                                  & /delzatw(i,1,k)
            ENDDO k_3134
         ENDIF
         
      ENDDO i_3132
      
   ENDIF
   
ENDDO ivar_310
   
!
!3.2 Diffusive fluxes of density
!-------------------------------
!

i_321: DO i=1,ncloc+1
   IF (node2du(i,1).EQ.2) THEN
      k_3211: DO k=1,nz
         xdifflux(i,k,3) = -beta_temp(i,1,k)*xdifflux(i,k,1)+&
                         &  beta_sal(i,1,k)*xdifflux(i,k,2)
      ENDDO k_3211
   ENDIF
ENDDO i_321

i_322: DO i=1,ncloc
   IF (maskatc_int(i,1).AND.zflag) THEN
      k_3221: DO k=2,nz
         zdifflux(i,k,3) = -beta_temp(i,1,k)*zdifflux(i,k,1)+&
                          & beta_sal(i,1,k)*zdifflux(i,k,2)
      ENDDO k_3221
   ENDIF
ENDDO i_322
   
!
!3.3 Projection of diffusive fluxes of T and S onto normal of level surface
!--------------------------------------------------------------------------
!      

ivar_330: DO ivar=1,3

   i_331: DO i=1,ncloc+1
      IF (node2du(i,1).EQ.2) THEN
         k_33311: DO k=1,nz
            xgrad(i,k) = (psic(i,k,ivar)-psic(i-1,k,ivar))/delxatu(i,1)
         ENDDO k_33311
      ELSE
         xgrad(i,:) = 0.0
      ENDIF
   ENDDO i_331
   i_332: DO i=1,ncloc
      IF (maskatc_int(i,1)) THEN
         k_3321: DO k=1,nz
            fldiffproj(i,k,ivar) = 0.5*(xdifflux(i,k,ivar)*xgrad(i,k)+&
                                      & xdifflux(i+1,k,ivar)*xgrad(i+1,k))
         ENDDO k_3321
      ENDIF
      IF (maskatc_int(i,1).AND.zflag) THEN
         k_3322: DO k=2,nz
            zgrad(k) = (psic(i,k,ivar)-psic(i,k-1,ivar))/delzatw(i,1,k)
         ENDDO k_3322
         zgrad(1) = 0.0; zgrad(nz+1) = 0.0
         k_3323: DO k=1,nz
            fldiffproj(i,k,ivar) = fldiffproj(i,k,ivar) + 0.5*&
                                 & (zdifflux(i,k,ivar)*zgrad(k)+&
                                  & zdifflux(i,k+1,ivar)*zgrad(k+1))
         ENDDO k_3323
      ENDIF
   ENDDO i_332
            
ENDDO ivar_330
   
!
!3.4 Variance            
!------------
!

nhdims = 0
ivar_340: DO ivar=1,3
!  ---global mean
   array2d = 0.0
   i_341: DO i=1,ncloc
      IF (maskatc_int(i,1)) THEN
         array2d(i,1) = SUM(psic(i,:,ivar)*delzatc(i,1,:))/depmeanatc(i,1)
      ENDIF
   ENDDO i_341
   CALL sum2_vars(array2d,psisum,nhdims,'C  ',0,commall=.TRUE.)
   psisum = psisum/(nc-1)
!  ---variance
   array2d = 0.0
   i_342: DO i=1,ncloc
      IF (maskatc_int(i,1)) THEN
         array2d(i,1) = SUM(delzatc(i,1,:)*(psic(i,:,ivar)-psisum)**2)&
                      & /depmeanatc(i,1)
      ENDIF
   ENDDO i_342
   CALL sum2_vars(array2d,psivar(ivar),nhdims,'C  ',0,commall=.TRUE.)
   IF (master) psivar(ivar) = psivar(ivar)/(nc-1)
ENDDO ivar_340

!
!3.5 Combine on master process
!-----------------------------
!      

ivar_350: DO ivar=1,3
         
!  ---global mean
   array2d = 0.0
   i_351: DO i=1,ncloc
      IF (node2du(i,1).EQ.2) THEN
         array2d(i,1) = SUM(xdifflux(i,:,ivar)*delzatu(i,1,:))/depmeanatu(i,1)
      ENDIF
   ENDDO i_351
   CALL sum2_vars(array2d,xdiffmean(ivar),nhdims,'U  ',0)
   IF (master) xdiffmean(ivar) = xdiffmean(ivar)/(nc-2)
   IF (zflag) THEN
      array2d = 0.0
      i_352: DO i=1,ncloc
         IF (maskatc_int(i,1)) THEN
            array2d(i,1) = SUM(zdifflux(i,2:nz,ivar)*delzatw(i,1,2:nz))&
                            & /depmeanatc(i,1)
         ENDIF
      ENDDO i_352
      CALL sum2_vars(array2d,zdiffmean(ivar),nhdims,'W  ',0)
      IF (master) zdiffmean(ivar) = zdiffmean(ivar)/(nc-1)
   ENDIF
   array2d = 0.0
   i_353: DO i=1,ncloc
      IF (maskatc_int(i,1)) THEN
         array2d(i,1) = SUM(fldiffproj(i,:,ivar)*delzatc(i,1,:))/depmeanatc(i,1)
      ENDIF
   ENDDO i_353
   CALL sum2_vars(array2d,fldiffprojmean(ivar),nhdims,'C  ',0)
   IF (master) fldiffprojmean(ivar) = fldiffprojmean(ivar)/(nc-1)
   
!  ---global maximum
   array3d = 0.0
   array3d(:,1,:) = xdifflux(1:ncloc,:,ivar)
   CALL max_vars(array3d,xdiffmax(ivar),0,&
               & mask=node2du(1:ncloc,1:nrloc).EQ.2)
   IF (zflag) THEN
      array3d = 0.0
      array3d(:,1,:) = zdifflux(:,1:nz,ivar)
      CALL max_vars(array3d,zdiffmax(ivar),0,mask=maskatc_int)
   ENDIF
   array3d = 0.0
   array3d(:,1,:) = fldiffproj(:,:,ivar)
   CALL max_vars(array3d,fldiffprojmax(ivar),0,mask=maskatc_int)
   
!  ---global minimum
   array3d = 0.0
   array3d(:,1,:) = xdifflux(1:ncloc,:,ivar)
   CALL min_vars(array3d,xdiffmin(ivar),0,&
               & mask=node2du(1:ncloc,1:nrloc).EQ.2)
   IF (zflag) THEN
      array3d = 0.0
      array3d(:,1,:) = zdifflux(:,1:nz,ivar)
      CALL min_vars(array3d,zdiffmin(ivar),0,mask=maskatc_int)
   ENDIF
   array3d = 0.0
   array3d(:,1,:) = fldiffproj(:,:,ivar)
   CALL min_vars(array3d,fldiffprojmin(ivar),0,mask=maskatc_int)
      
ENDDO ivar_350
   
!
!4. Write parameters
!-------------------
!

IF (master) THEN
   iday = nosecsrun/86400

   WRITE (iunit,9001) iday
   WRITE (iunit,9002) 'xdiffmean_temp', xdiffmean(1)
   WRITE (iunit,9002) 'xdiffmean_sal', xdiffmean(2)
   WRITE (iunit,9002) 'xdiffmean_dens', xdiffmean(3)
   IF (zflag) THEN
      WRITE (iunit,9002) 'zdiffmean_temp', zdiffmean(1)
      WRITE (iunit,9002) 'zdiffmean_sal', zdiffmean(2)
      WRITE (iunit,9002) 'zdiffmean_dens', zdiffmean(3)
   ENDIF
   WRITE (iunit,9002) 'diffprojmean_temp', fldiffprojmean(1)
   WRITE (iunit,9002) 'diffprojmean_sal', fldiffprojmean(2)
   WRITE (iunit,9002) 'diffprojmean_dens', fldiffprojmean(3)
   WRITE (iunit,9002) 'xdiffmax_temp', xdiffmax(1)
   WRITE (iunit,9002) 'xdiffmax_sal', xdiffmax(2)
   WRITE (iunit,9002) 'xdiffmax_dens', xdiffmax(3)
   IF (zflag) THEN
      WRITE (iunit,9002) 'zdiffmax_temp', zdiffmax(1)
      WRITE (iunit,9002) 'zdiffmax_sal', zdiffmax(2)
      WRITE (iunit,9002) 'zdiffmax_dens', zdiffmax(3)
   ENDIF
   WRITE (iunit,9002) 'diffprojmax_temp', fldiffprojmax(1)
   WRITE (iunit,9002) 'diffprojmax_sal', fldiffprojmax(2)
   WRITE (iunit,9002) 'diffprojmax_dens', fldiffprojmax(3)
   WRITE (iunit,9002) 'xdiffmin_temp', xdiffmin(1)
   WRITE (iunit,9002) 'xdiffmin_sal', xdiffmin(2)
   WRITE (iunit,9002) 'xdiffmin_dens', xdiffmin(3)
   IF (zflag) THEN
      WRITE (iunit,9002) 'zdiffmin_temp', zdiffmin(1)
      WRITE (iunit,9002) 'zdiffmin_sal', zdiffmin(2)
      WRITE (iunit,9002) 'zdiffmin_dens', zdiffmin(3)
   ENDIF
   WRITE (iunit,9002) 'diffprojmin_temp', fldiffprojmin(1)
   WRITE (iunit,9002) 'diffprojmin_sal', fldiffprojmin(2)
   WRITE (iunit,9002) 'diffprojmin_dens', fldiffprojmin(3)
   WRITE (iunit,9002) 'tempvar', psivar(1)
   WRITE (iunit,9002) 'salvar', psivar(2)
   WRITE (iunit,9002) 'densvar', psivar(3)

ENDIF

!
!5. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')

CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',I3,' days')
9002 FORMAT (T3,A,T20,': ',G12.4E3)

END SUBROUTINE usrdef_output
