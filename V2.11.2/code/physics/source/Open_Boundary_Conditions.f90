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
! *Open_boundary_conditions* Apply open boundary conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Conditions.f90  V2.9
!
! $Date: 2017-10-10 13:57:00 +0200 (Tue, 10 Oct 2017) $
!
! $Revision: 1055 $
!
! Description -
!
! Routines - open_boundary_conds_2d, open_boundary_conds_impl,
!            open_boundary_conds_3d, open_boundary_conds_prof,
!            open_boundary_outflow_time
!
!************************************************************************
!

!========================================================================

SUBROUTINE open_boundary_conds_2d
!************************************************************************
!
! *open_boundary_conds_2d* Apply open boundary conditions for 2-D mode
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_conditions.f90  V2.9
!
! Description - 
!
! Reference -
!
! Calling program - correct_free_surf, current_2d
!
! External calls -
!
! Module calls - cfl_orlan, Uvar_at_C, Vvar_at_C
!
!************************************************************************
!
USE currents
USE depths
USE density
USE fluxes
USE grid
USE gridpars
USE iopars
USE obconds
USE physpars
USE switches
USE tide
USE timepars
USE wavevars
USE array_interp, ONLY: Uvar_at_C, Vvar_at_C
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: cfl_orlan

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, ic, ii, iiloc, is, isign, iu1, iu2, j, jc, jj, jjloc, js, &
         & jsign, jv1, jv2, npcc
REAL :: a, ac, b, c, cfl, chin, dxu, dyv, grx, gry, obrhs, udx, uvelatc, vdy, &
      & vvelatc, zetobu, zetobudat, zetobv, zetobvdat
REAL, DIMENSION(2) :: udvals


procname(pglev+1) = 'open_boundary_conds_2d'
CALL log_timer_in(npcc)

!
!1. U-nodes
!----------
!

iiloc_110: DO iiloc=1,nobuloc_ext
   i = iobuloc(iiloc); j = jobuloc(iiloc)
   ii = indexobu(iiloc)
   IF (node2du(i,j).LE.1.OR.(westobu(ii).AND.i.EQ.ncloc+1)) CYCLE iiloc_110
   ic = MERGE(i,i-1,westobu(ii))
   isign = MERGE(1,-1,westobu(ii)); iu1 = MERGE(i+1,i-1,westobu(ii))

!  ---define parameters
   SELECT CASE (ityp2dobu(ii))
      CASE (3,8:13)
         IF (iloczobu(ii).EQ.1) THEN
            is = 2
         ELSEIF (iloczobu(ii).EQ.2) THEN
            is = 1
         ENDIF
   END SELECT
   SELECT CASE (ityp2dobu(ii))
      CASE (3,5,8:13)
         c = SQRT(gaccatc(ic,j)*deptotatc(ic,j))
         a = delt2d*c/delxatc(ic,j)
         ac = delt2d*gaccatc(ic,j)*deptotatc(ic,j)/delxatc(ic,j)
   END SELECT
   SELECT CASE (ityp2dobu(ii))
      CASE (3,9:13)
         grx = delxatc(ic,j)/delxatu(i,j)
   END SELECT
   SELECT CASE (ityp2dobu(ii))
   CASE (3,9,12)
      IF (iopt_hydro_impl.EQ.0) THEN
         zetobu = zeta(ic,j); zetobudat = zdatobu(ii,2)
      ELSE
         zetobu = (1.0-theta_sur)*zeta_old(ic,j)+theta_sur*zeta(ic,j)
         zetobudat = (1.0-theta_sur)*zdatobu_old(ii)+theta_sur*zdatobu(ii,2)
      ENDIF
   END SELECT
   SELECT CASE (ityp2dobu(ii))
      CASE (11:13)
         b = 1.0/(1.0+1.5*a*grx)
   END SELECT
   SELECT CASE (ityp2dobu(ii))
      CASE (1,3,9:13,17)
         vvelatc = theta_cor*Vvar_at_C(vdvel(ic,j:j+1),ic,j,nz,1,1)+&
                   (1.0-theta_cor)*Vvar_at_C(vdvel_old(ic,j:j+1),ic,j,nz,1,1)
         obrhs = delt2d*(coriolatu(i,j)*vvelatc+deptotatu(i,j)*fxastro(i,j)+&
                       & usstresatu(i,j)-ubstresatu(i,j)) 
         IF (iopt_waves_curr.EQ.1) THEN
            obrhs = obrhs + delt2d*(&
                    & coriolatu(i,j)*deptotatc(ic,j)*vmstokesatc(ic,j))
         ENDIF
         IF (iopt_waves_dissip.EQ.1) THEN
            obrhs = obrhs + delt2d*(umswdissipatc(ic,j)+umbwdissipatc(ic,j))
         ENDIF
   END SELECT
   SELECT CASE (ityp2dobu(ii))
      CASE (11:13)
         udx = Uvar_at_C(udvel_old(ic:ic+1,j),ic,j,nz,1,1)*&
             & (delyatu(iu1,j)-delyatu(i,j))/delyatc(ic,j)
   END SELECT
   SELECT CASE (ityp2dobu(ii))
      CASE (10:13)
         dyv = (delxatv(ic,j+1)*vdvel_old(ic,j+1)-&
              & delxatv(ic,j)*vdvel_old(ic,j))/delyatc(ic,j)
   END SELECT

!  ---local solution for momentum
   SELECT CASE (ityp2dobu(ii))
      CASE (9,12)
         udatobu(ii,2) = obc2uvatu_old(ii,2) - isign*is*ac*grx*&
                       & (zetobu-zetobudat) + obrhs
         obc2uvatu(ii,2) = udatobu(ii,2)
      CASE (10)
         udatobu(ii,2) = obc2uvatu_old(ii,2) + obrhs
         obc2uvatu(ii,2) = udatobu(ii,2)
   END SELECT

!  ---apply open boundary condition
   SELECT CASE (ityp2dobu(ii))
      CASE (1); udvel(i,j) = udvel_old(i,j) + obrhs
      CASE (2); udvel(i,j) = udvel(iu1,j)
      CASE (3)
         udvel(i,j) = udvel_old(i,j) - &
                    & isign*is*ac*grx*(zetobu-zetobudat) + obrhs
      CASE (4); udvel(i,j) = udatobu(ii,2)
      CASE (5); udvel(i,j) = (udvel_old(i,j)+a*udvel(iu1,j))/(1.0+a)
      CASE (6)
         iu2 = MERGE(i+2,i-2,westobu(ii))
         cfl = cfl_orlan(udvel_old(iu1,j),obc2uvatu_old(ii,1),&
                       & obc2uvatu_old(ii,2))
         obc2uvatu(ii,1) = udvel_old(iu1,j)
         obc2uvatu(ii,2) = udvel_old(iu2,j)
         udvel(i,j) = (1.0-cfl)*udvel_old(i,j) + cfl*udvel_old(iu1,j)
      CASE (7)
         iu2 = MERGE(i+2,i-2,westobu(ii))
         udvel(i,j) = MERGE(udvel_old(i,j),udvel_old(iu1,j),&
                          & udvel_old(iu1,j).LT.obc2uvatu_old(ii,2))
         obc2uvatu(ii,2) = udvel_old(iu2,j)
      CASE (8:9)
         udvel(i,j) = udatobu(ii,2)-0.5*isign*is*c*(zeta(ic,j)-zdatobu(ii,2))
      CASE (10)
         obc2uvatu(ii,1) = obc2uvatu_old(ii,1) - delt2d*dyv/delxatc(ic,j)
         udvel(i,j) = udatobu(ii,2) - isign*c*(zeta(ic,j)-obc2uvatu(ii,1))
      CASE (11:12)
         chin = udatobu(ii,2)+0.5*isign*is*c*(zdatobu(ii,2)+(2-is)*zeta(ic,j))
         obc2uvatu(ii,1) = b*(obc2uvatu_old(ii,1)+a*grx*(udvel(iu1,j)+&
                          & 0.5*chin-2.0*isign*c*zeta(ic,j))+&
                          & a*(isign*dyv+udx)+obrhs)
         udvel(i,j) = 0.5*(obc2uvatu(ii,1)+chin)
      CASE (13)
         obc2uvatu(ii,2) = obc2uvatu_old(ii,2) - a*(isign*dyv+udx)+obrhs
         obc2uvatu(ii,1) = b*(obc2uvatu_old(ii,1)+a*grx*(udvel(iu1,j)+&
                          & 0.5*obc2uvatu(ii,2)-2.0*isign*c*zeta(ic,j))+&
                          & a*(isign*dyv+udx)+obrhs)
         obc2uvatu(ii,2) = obc2uvatu_old(ii,2) + obrhs
         udvel(i,j) = 0.5*(obc2uvatu(ii,1)+obc2uvatu(ii,2))
      CASE (14)                                      
         udvel(i,j) = udatobu(ii,2)*deptotatu(i,j) 
      CASE (15:16)                                      
         udvel(i,j) = udatobu(ii,2)/delyatu(i,j)    
      CASE (17)                                      
         udvel(i,j) = udvel_old(i,j) - &
                    & delt2d*gaccatu(i,j)*deptotatu(i,j)*zdatobu(ii,2)+obrhs 
   END SELECT

ENDDO iiloc_110

!
!2. V-nodes
!----------
!

jjloc_210: DO jjloc=1,nobvloc_ext
   i = iobvloc(jjloc); j = jobvloc(jjloc)
   jj = indexobv(jjloc)
   IF (node2dv(i,j).LE.1.OR.(soutobv(jj).AND.j.EQ.nrloc+1)) CYCLE jjloc_210
   jc = MERGE(j,j-1,soutobv(jj))
   jsign = MERGE(1,-1,soutobv(jj)); jv1 = MERGE(j+1,j-1,soutobv(jj))

!  ---define parameters
   SELECT CASE (ityp2dobv(jj))
      CASE (3,8:13)
         IF (iloczobv(jj).EQ.1) THEN
            js = 2
         ELSEIF (iloczobv(jj).EQ.2) THEN
            js = 1
         ENDIF
   END SELECT
   SELECT CASE (ityp2dobv(jj))
      CASE (3,5,8:13)
         c = SQRT(gaccatc(i,jc)*deptotatc(i,jc))
         a = delt2d*c/delyatc(i,jc)
         ac = delt2d*gaccatc(i,jc)*deptotatc(i,jc)/delyatc(i,jc)
   END SELECT
   SELECT CASE (ityp2dobv(jj))
      CASE (3,9:13)
         gry = delyatc(i,jc)/delyatv(i,j)
   END SELECT
   SELECT CASE (ityp2dobv(jj))
      CASE (3,9,12)
         IF (iopt_hydro_impl.EQ.0) THEN
            zetobv = zeta(i,jc); zetobvdat = zdatobv(jj,2)
         ELSE
            zetobv = (1.0-theta_sur)*zeta_old(i,jc)+theta_sur*zeta(i,jc)
            zetobvdat = (1.0-theta_sur)*zdatobv_old(jj)+theta_sur*zdatobv(jj,2)
         ENDIF
   END SELECT
   SELECT CASE (ityp2dobv(jj))
      CASE (11:13)
         b = 1.0/(1.0+1.5*a*gry)
   END SELECT
   SELECT CASE (ityp2dobv(jj))
      CASE (1,3,9:13,17)
         udvals = MERGE(udvel_old(i:i+1,jc),udvel(i:i+1,jc),&
                      & node2du(i,jc).GT.2)
         uvelatc = 0.5*theta_cor*(udvals(1)+udvals(2))+&
                   (1.0-theta_cor)*Uvar_at_C(udvel_old(i:i+1,jc),i,jc,nz,1,1)
         obrhs = delt2d*(-coriolatv(i,j)*uvelatc+deptotatv(i,j)*fyastro(i,j)+&
                       & vsstresatv(i,j)-vbstresatv(i,j))
         IF (iopt_waves_curr.EQ.1) THEN
            obrhs = obrhs + delt2d*(&
                    & -coriolatv(i,j)*deptotatc(i,jc)*umstokesatc(i,jc))
         ENDIF
         IF (iopt_waves_dissip.EQ.1) THEN
            obrhs = obrhs + delt2d*(vmswdissipatc(i,jc)+vmbwdissipatc(i,jc))
         ENDIF
         
   END SELECT
   SELECT CASE (ityp2dobv(jj))
      CASE (11:13)
         vdy = Vvar_at_C(vdvel_old(i,jc:jc+1),i,jc,nz,1,1)*&
             & (delxatv(i,jv1)-delxatv(i,j))/delxatc(i,jc)
   END SELECT
   SELECT CASE (ityp2dobv(jj))
      CASE (10:13)
         dxu = (delyatu(i+1,jc)*udvel_old(i+1,jc)-&
              & delyatu(i,jc)*udvel_old(i,jc))/delxatc(i,jc)
   END SELECT

!  ---local solution for momentum
   SELECT CASE (ityp2dobv(jj))
      CASE (9,12)
         vdatobv(jj,2) = obc2uvatv_old(jj,2)-jsign*js*a*c*gry*&
                       & (zetobv-zetobvdat)+obrhs
         obc2uvatv(jj,2) = vdatobv(jj,2)
      CASE (10)
         vdatobv(jj,2) = obc2uvatv_old(jj,2)+obrhs
         obc2uvatv(jj,2) = vdatobv(jj,2)
   END SELECT

!  ---apply open boundary condition
   SELECT CASE (ityp2dobv(jj))
      CASE (1); vdvel(i,j) = vdvel_old(i,j) + obrhs
      CASE (2); vdvel(i,j) = vdvel(i,jv1)
      CASE (3)
         vdvel(i,j) = vdvel_old(i,j) -&
                    & jsign*js*a*c*gry*(zetobv-zetobvdat) + obrhs
      CASE (4); vdvel(i,j) = vdatobv(jj,2)
      CASE (5); vdvel(i,j) = (vdvel_old(i,j)+a*vdvel(i,jv1))/(1.0+a)
      CASE (6)
         jv2 = MERGE(j+2,j-2,soutobv(jj))
         cfl = cfl_orlan(vdvel_old(i,jv1),obc2uvatv_old(jj,1),&
                       & obc2uvatv_old(jj,2))
         obc2uvatv(jj,1) = vdvel_old(i,jv1)
         obc2uvatv(jj,2) = vdvel_old(i,jv2)
         vdvel(i,j) = (1.0-cfl)*vdvel_old(i,j) + cfl*vdvel_old(i,jv1)
      CASE (7)
         jv2 = MERGE(j+2,j-2,soutobv(jj))
         vdvel(i,j) = MERGE(vdvel_old(i,j),vdvel_old(i,jv1),&
                          & vdvel_old(i,jv1).LT.obc2uvatv_old(jj,2))
         obc2uvatv(jj,2) = vdvel_old(i,jv2)
      CASE (8:9)
         vdvel(i,j) = vdatobv(jj,2)-0.5*jsign*js*c*(zeta(i,jc)-zdatobv(jj,2))
      CASE (10)
         obc2uvatv(jj,1) = obc2uvatv_old(jj,1) - delt2d*dxu/delyatc(i,jc)
         vdvel(i,j) = vdatobv(jj,2) - jsign*c*(zeta(i,jc)-obc2uvatv(jj,1))
      CASE (11:12)
         chin = vdatobv(jj,2)+0.5*jsign*js*c*(zdatobv(jj,2)+(2-js)*zeta(i,jc))
         obc2uvatv(jj,1) = b*(obc2uvatv_old(jj,1)+a*gry*(vdvel(i,jv1)+&
                          & 0.5*chin-2.0*jsign*c*zeta(i,jc))+&
                          & a*(jsign*dxu+vdy)+obrhs)
         vdvel(i,j) = 0.5*(obc2uvatv(jj,1)+chin)
      CASE (13)
         obc2uvatv(jj,2) = obc2uvatv_old(jj,2) - a*(jsign*dxu+vdy)+obrhs
         obc2uvatv(jj,1) = b*(obc2uvatv_old(jj,1)+a*gry*(vdvel(i,jv1)+&
                          & 0.5*obc2uvatv(jj,2)-2.0*jsign*c*zeta(i,jc))+&
                          & a*(jsign*dxu+vdy)+obrhs)
         obc2uvatv(jj,2) = obc2uvatv_old(jj,2) + obrhs
         vdvel(i,j) = 0.5*(obc2uvatv(jj,1)+obc2uvatv(jj,2))
      CASE (14)                                     
         vdvel(i,j) = vdatobv(jj,2)*deptotatv(i,j)
      CASE (15:16)                                     
         vdvel(i,j) = vdatobv(jj,2)/delxatv(i,j)
      CASE (17)                                     
         vdvel(i,j) = vdvel_old(i,j) - &
                    & delt2d*gaccatv(i,j)*deptotatv(i,j)*zdatobv(jj,2)+obrhs 

   END SELECT

ENDDO jjloc_210

CALL log_timer_out(npcc,itm_bconds)


RETURN

END SUBROUTINE open_boundary_conds_2d

!========================================================================

SUBROUTINE open_boundary_conds_2d_impl(mglevel)
!************************************************************************
!
! *open_boundary_conds_2d_impl* Store the implicit terms in the solution matrix
!                               and vector arising from the 2-D open boundary
!                               conditions in case an implicit scheme is used
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)Open_Boundary_Conditions.f90  V2.9
!
! Description - 
!
! Reference -
!
! Calling program - correct_free_surf
!
! External calls -
!
! Module calls - cfl_orlan, mat_add, vec_add, Uvar_at_C, Vvar_at_C
!
!************************************************************************
!
USE currents
USE depths
USE fluxes
USE grid
USE gridpars
USE iopars
USE multigrid
USE obconds
USE physpars
USE switches
USE syspars
USE tide
USE timepars
USE wavevars
USE array_interp, ONLY: Uvar_at_C, Vvar_at_C
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: cfl_orlan

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: mglevel

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ic, ii, iiloc, is, isign, iu1, iu2, j, jc, jj, jjloc, js, jsign, &
         & jv1, jv2, lev, npcc
REAL :: a, b, c, cfl, chin, dxu, dyv, grx, gry, Macoeff, Mvcoeff, obrhs, udx, &
       & uvelatc, vdy, vvelatc, zetobudat, zetobvdat
REAL (KIND=kndrlong) :: Ma, Ma1, Mv


IF (maxdatafiles(io_2uvobc,1).EQ.0) RETURN

procname(pglev+1) = 'open_boundary_conds_2d_impl'
CALL log_timer_in(npcc)

lev = mglevel

!
!1. Fine grid level
!------------------
!

IF (lev.EQ.0) THEN

!
!1.1 U-nodes
!-----------
!

   iiloc_110: DO iiloc=1,nobuloc_ext
      i = iobuloc(iiloc); j = jobuloc(iiloc)
      ii = indexobu(iiloc)
      IF (node2du(i,j).LE.1.OR.&
       & (westobu(ii).AND.i.EQ.ncloc+1).OR.(.NOT.westobu(ii).AND.i.EQ.1)) THEN
         CYCLE iiloc_110
      ENDIF
      ic = MERGE(i,i-1,westobu(ii))
      isign = MERGE(1,-1,westobu(ii)); iu1 = MERGE(i+1,i-1,westobu(ii))
!     ---define parameters
      SELECT CASE (ityp2dobu(ii))
         CASE (3,8:13)
            IF (iloczobu(ii).EQ.1) THEN
               is = 2
            ELSEIF (iloczobu(ii).EQ.2) THEN
               is = 1
            ENDIF
      END SELECT
      SELECT CASE (ityp2dobu(ii))
         CASE (3,5,8:13)
            c = SQRT(gaccatc(ic,j)*deptotatc(ic,j))
            a = delt2d*c/delxatc(ic,j)
      END SELECT
      SELECT CASE (ityp2dobu(ii))
         CASE (3,9:13)
            grx = delxatc(ic,j)/delxatu(i,j)
      END SELECT
      SELECT CASE (ityp2dobu(ii))
         CASE (3,9,12)
            zetobudat = (1.0-theta_sur)*zdatobu_old(ii) + &
                      & theta_sur*zdatobu(ii,2)
      END SELECT
      SELECT CASE (ityp2dobu(ii))
         CASE (11:13)
            b = 1.0/(1.0+1.5*a*grx)
      END SELECT
      SELECT CASE (ityp2dobu(ii))
         CASE (1,3,9:13,17)
            vvelatc = Vvar_at_C(vdvel_old(ic,j:j+1),ic,j,nz,1,1)
            obrhs = delt2d*(coriolatu(i,j)*vvelatc+deptotatu(i,j)*fxastro(i,j)+&
                          & usstresatu(i,j)-ubstresatu(i,j)) 
            IF (iopt_waves_curr.EQ.1) THEN
               obrhs = obrhs + delt2d*(&
                     & coriolatu(i,j)*deptotatc(ic,j)*vmstokesatc(ic,j)+&
                     & umswdissipatc(ic,j)+umbwdissipatc(ic,j))
            ENDIF
      END SELECT
      SELECT CASE (ityp2dobu(ii))
         CASE (11:13)
            udx = Uvar_at_C(udvel_old(ic:ic+1,j),ic,j,nz,1,1)*&
               & (delyatu(iu1,j)-delyatu(i,j))/delyatc(ic,j)
      END SELECT
      SELECT CASE (ityp2dobu(ii))
         CASE (10:13)
            dyv = (delxatv(ic,j+1)*vdvel_old(ic,j+1)-&
                 & delxatv(ic,j)*vdvel_old(ic,j))/delyatc(ic,j)
      END SELECT

!     ---local solution for momentum
      SELECT CASE (ityp2dobu(ii))
         CASE (9,12)
            udatobu(ii,2) = obc2uvatu_old(ii,2)-isign*is*a*c*grx*&
                         & (zeta(ic,j)-zetobudat)+obrhs
            obc2uvatu(ii,2) = udatobu(ii,2)
         CASE (10)
            udatobu(ii,2) = obc2uvatu_old(ii,2)+obrhs
            obc2uvatu(ii,2) = udatobu(ii,2)
      END SELECT

!     ---calculate entries for matrix system
      Mvcoeff = delyatu(i,j)*isign
      Macoeff = Mvcoeff*gaccatu(iu1,j)*delt3d/delxatu(iu1,j)*deptotatu(iu1,j)*&
              & isign
      Mv = 0.0; Ma = 0.0; Ma1= 0.0
      SELECT CASE (ityp2dobu(ii))
         CASE (0)
            Mv = Mvcoeff*udvel_old(i,j)
         CASE (1) 
            Mv = Mvcoeff*(udvel_old(i,j) + obrhs)
         CASE (2)
            Mv = Mvcoeff*deptotatu(iu1,j)*umpred(iu1,j)
            Ma1= Macoeff 
         CASE (3)
            Mv = Mvcoeff*(udvel_old(i,j)-&
               & isign*is*a*c*grx*(zeta(ic,j)-zetobudat)+ obrhs)
            Ma = Mvcoeff*isign*is*a*c*grx
         CASE (4) 
            Mv = Mvcoeff*udatobu(ii,2)
         CASE (5) 
            Mv = Mvcoeff*(udvel_old(i,j)+&
                        & a*deptotatu(iu1,j)*umpred(iu1,j))/(1.0+a)
            Ma1= Macoeff*a/(1.0+a)
         CASE (6)
            iu2 = MERGE(i+2,i-2,westobu(ii))
            cfl = cfl_orlan(udvel_old(iu1,j),obc2uvatu_old(ii,1),&
                          & obc2uvatu_old(ii,2))
            obc2uvatu(ii,1) = udvel_old(iu1,j)
            obc2uvatu(ii,2) = udvel_old(iu2,j)
            Mv = Mvcoeff*((1.0-cfl)*udvel_old(i,j)+cfl*udvel_old(iu1,j))
         CASE (7)
            iu2 = MERGE(i+2,i-2,westobu(ii))
            obc2uvatu(ii,2) = udvel_old(iu2,j)
            Mv = Mvcoeff*MERGE(udvel_old(i,j),udvel_old(iu1,j),&
                             & udvel_old(iu1,j).LT.obc2uvatu_old(ii,2))
         CASE (8)
            Mv = Mvcoeff*(udatobu(ii,2)-0.5*isign*is*c*(zeta(ic,j)-&
                        & zdatobu(ii,2)))
            Ma = Mvcoeff*0.5*isign*is*c
         CASE (9)
            Mv = Mvcoeff*(udatobu(ii,2)-0.5*isign*is*c*(zeta(ic,j)-&
                        & zdatobu(ii,2)))
            Ma = Mvcoeff*isign*is*c*(0.5+a*grx)
         CASE (10)
            obc2uvatu(ii,1) = obc2uvatu_old(ii,1) - delt2d*dyv/delxatc(ic,j)
            Mv = Mvcoeff*(udatobu(ii,2)-isign*c*(zeta(ic,j)-obc2uvatu(ii,1)))
            Ma = Mvcoeff*isign*c*(1.0 + is*a*grx)
         CASE (11)
            chin = udatobu(ii,2)+&
                 & 0.5*isign*is*c*(zdatobu(ii,2)+(2-is)*zeta(ic,j))
            obc2uvatu(ii,1) = b*(obc2uvatu_old(ii,1)+&
                            & a*grx*(deptotatu(iu1,j)*umpred(iu1,j)+&
                            & 0.5*chin-2.0*isign*c*zeta(ic,j))+&
                            & a*(isign*dyv+udx)+obrhs)
            Mv = Mvcoeff*0.5*(obc2uvatu(ii,1)+chin)
            Ma = -Mvcoeff*0.25*isign*is*c*(2-is)
            Ma = Ma + 0.5*b*a*grx*Ma 
            Ma = Ma + Mvcoeff*b*a*grx*isign*c
            Ma1= Macoeff*0.5*b*a*grx
         CASE (12)
            chin = udatobu(ii,2)+0.5*isign*is*c*&
                               & (zdatobu(ii,2)+(2-is)*zeta(ic,j))
            obc2uvatu(ii,1) = b*(obc2uvatu_old(ii,1)+&
                               & a*grx*(deptotatu(iu1,j)*umpred(iu1,j)+&
                               & 0.5*chin-2.0*isign*c*zeta(ic,j))+&
                               & a*(isign*dyv+udx)+obrhs)
            Mv = Mvcoeff*0.5*(obc2uvatu(ii,1)+chin)
            Ma = Mvcoeff*0.5*isign*is*c*(0.5*(is-2)+a*grx )
            Ma = Ma + 0.5*b*a*grx*Ma 
            Ma = Ma + Mvcoeff*b*a*grx*isign*c
            Ma1= Macoeff*0.5*b*a*grx
         CASE (13)
            obc2uvatu(ii,2) = obc2uvatu_old(ii,2) - a*(isign*dyv+udx)+obrhs
            obc2uvatu(ii,1) = b*(obc2uvatu_old(ii,1)+&
                               & a*grx*(deptotatu(iu1,j)*umpred(iu1,j)+&
                               & 0.5*obc2uvatu(ii,2)-2.0*isign*c*zeta(ic,j))+&
                               & a*(isign*dyv+udx)+obrhs)
            obc2uvatu(ii,2) = obc2uvatu_old(ii,2) + obrhs
            Mv = Mvcoeff*0.5*(obc2uvatu(ii,1)+obc2uvatu(ii,2))
            Ma = Mvcoeff*b*a*grx*isign*c
            Ma1= Macoeff*0.5*b*a*grx
         CASE (14)                                      
            Mv = Mvcoeff*udatobu(ii,2)*deptotatu(i,j) 
         CASE (15:16)                                      
            Mv = Mvcoeff*udatobu(ii,2)/delyatu(i,j)    
         CASE (17)                                      
            Mv = Mvcoeff*(udvel_old(i,j)-&
                       & delt2d*gaccatu(i,j)*deptotatu(i,j)*zdatobu(ii,2)+obrhs)
      END SELECT

!     ---insert in matrix
      Ma = Ma*theta_sur; Ma1 = Ma1*theta_sur
      mgvars(0)%rhs(ic,j) = mgvars(0)%rhs(ic,j) + Mv
      mgvars(0)%subdiags(ic,j,0) = mgvars(0)%subdiags(ic,j,0) + Ma - Ma1
      IF (westobu(ii)) THEN
         mgvars(0)%subdiags(ic,j,1) = mgvars(0)%subdiags(ic,j,1) + Ma1
      ELSE
         mgvars(0)%subdiags(ic,j,-2) = mgvars(0)%subdiags(ic,j,-2) + Ma1
      ENDIF
   ENDDO iiloc_110

!
!1.2 V-nodes
!-----------
!

   jjloc_120: DO jjloc=1,nobvloc_ext
      i = iobvloc(jjloc); j = jobvloc(jjloc)
      jj = indexobv(jjloc)
      IF (node2dv(i,j).LE.1.OR.&
       & (soutobv(jj).AND.j.EQ.nrloc+1).OR.(.NOT.soutobv(jj).AND.j.EQ.1)) THEN
         CYCLE jjloc_120
      ENDIF
      jc = MERGE(j,j-1,soutobv(jj))
      jsign = MERGE(1,-1,soutobv(jj)); jv1 = MERGE(j+1,j-1,soutobv(jj))
!     ---define parameters
      SELECT CASE (ityp2dobv(jj))
         CASE (3,8:13)
            IF (iloczobv(jj).EQ.1) THEN
               js = 2
            ELSEIF (iloczobv(jj).EQ.2) THEN
               js = 1
            ENDIF
      END SELECT
      SELECT CASE (ityp2dobv(jj))
         CASE (3,5,8:13)
            c = SQRT(gaccatc(i,jc)*deptotatc(i,jc))
            a = delt2d*c/delyatc(i,jc)
      END SELECT
      SELECT CASE (ityp2dobv(jj))
         CASE (3,9:13)
            gry = delyatc(i,jc)/delyatv(i,j)
      END SELECT
      SELECT CASE (ityp2dobv(jj))
         CASE (3,9,12)
            zetobvdat = (1.0-theta_sur)*zdatobv_old(jj) + &
                      & theta_sur*zdatobv(jj,2)
      END SELECT
      SELECT CASE (ityp2dobv(jj))
         CASE (11:13)
            b = 1.0/(1.0+1.5*a*gry)
      END SELECT
      SELECT CASE (ityp2dobv(jj))
         CASE (1,3,9:13,17)
            uvelatc = Uvar_at_C(udvel_old(i:i+1,jc),i,jc,nz,1,1)
            obrhs = delt2d*(-coriolatv(i,j)*uvelatc+deptotatv(i,j)*fyastro(i,j)&
                          & +vsstresatv(i,j)-vbstresatv(i,j))
            IF (iopt_waves_curr.EQ.1) THEN
               obrhs = obrhs + delt2d*(&
                     & -coriolatv(i,j)*deptotatc(i,jc)*umstokesatc(i,jc)+&
                     & vmswdissipatc(i,jc)+vmbwdissipatc(i,jc))
            ENDIF
      END SELECT
      SELECT CASE (ityp2dobv(jj))
         CASE (11:13)
            vdy = Vvar_at_C(vdvel_old(i,jc:jc+1),i,jc,nz,1,1)*&
             & (delxatv(i,jv1)-delxatv(i,j))/delxatc(i,jc)
      END SELECT
      SELECT CASE (ityp2dobv(jj))
         CASE (10:13)
            dxu = (delyatu(i+1,jc)*udvel_old(i+1,jc)-&
                 & delyatu(i,jc)*udvel_old(i,jc))/delxatc(i,jc)
      END SELECT

!     ---local solution for momentum
      SELECT CASE (ityp2dobv(jj))
         CASE (9,12)
            vdatobv(jj,2) = obc2uvatv_old(jj,2)-jsign*js*a*c*gry*&
                         & (zeta(i,jc)-zetobvdat)+obrhs
            obc2uvatv(jj,2) = vdatobv(jj,2)
         CASE (10)
            vdatobv(jj,2) = obc2uvatv_old(jj,2)+obrhs
            obc2uvatv(jj,2) = vdatobv(jj,2)
      END SELECT

!     ---calculate entries for matrix system
      Mvcoeff = delxatv(i,j)*jsign
      Macoeff = Mvcoeff*gaccatv(i,jv1)*delt3d/delyatv(i,jv1)*deptotatv(i,jv1)*&
              & jsign
      Mv = 0.0; Ma = 0.0; Ma1= 0.0
      SELECT CASE (ityp2dobv(jj))
         CASE (0)
            Mv = Mvcoeff*vdvel_old(i,j)
         CASE (1)
            Mv = Mvcoeff*(vdvel_old(i,j)+obrhs)
         CASE (2)
            Mv = Mvcoeff*deptotatv(i,jv1)*vmpred(i,jv1)
            Ma1 = Macoeff
         CASE (3)
            Mv = Mvcoeff*(vdvel_old(i,j)-&
               & jsign*js*a*c*gry*(zeta(i,jc)-zetobvdat)+obrhs)
            Ma = Mvcoeff*jsign*js*a*c*gry
         CASE (4)
            Mv = Mvcoeff*vdatobv(jj,2)
         CASE (5)
            Mv = Mvcoeff*(vdvel_old(i,j)+&
                        & a*deptotatv(i,jv1)*vmpred(i,jv1))/(1.0+a)
            Ma1= Macoeff*a/(1.0+a)
         CASE (6)
            jv2 = MERGE(j+2,j-2,soutobv(jj))
            cfl = cfl_orlan(vdvel_old(i,jv1),obc2uvatv_old(jj,1),&
                          & obc2uvatv_old(jj,2))
            obc2uvatv(jj,1) = vdvel_old(i,jv1)
            obc2uvatv(jj,2) = vdvel_old(i,jv2)
            Mv = Mvcoeff*((1.0-cfl)*vdvel_old(i,j)+cfl*vdvel_old(i,jv1))
         CASE (7)
            jv2 = MERGE(j+2,j-2,soutobv(jj))
            obc2uvatv(jj,2) = vdvel_old(i,jv2)
            Mv = Mvcoeff*MERGE(vdvel_old(i,j),vdvel_old(i,jv1),&
                             & vdvel_old(i,jv1).LT.obc2uvatv_old(jj,2))
         CASE (8)
            Mv = Mvcoeff*(vdatobv(jj,2)-&
                        & 0.5*jsign*js*c*(zeta(i,jc)-zdatobv(jj,2)))
            Ma = Mvcoeff*0.5*jsign*js*c
         CASE (9)
            Mv = Mvcoeff*(vdatobv(jj,2)-&
                        & 0.5*jsign*js*c*(zeta(i,jc)-zdatobv(jj,2)))
            Ma = Mvcoeff*jsign*js*c*(0.5+a*gry)
         CASE (10)
            obc2uvatv(jj,1) = obc2uvatv_old(jj,1) - delt2d*dxu/delyatc(i,jc)
            Mv = Mvcoeff*(vdatobv(jj,2)-jsign*c*(zeta(i,jc)-obc2uvatv(jj,1)))
            Ma = Mvcoeff*jsign*c*(1.0+js*a*gry) 
         CASE (11)
            chin = vdatobv(jj,2)+&
                 & 0.5*jsign*js*c*(zdatobv(jj,2)+(2-js)*zeta(i,jc))
            obc2uvatv(jj,1) = b*(obc2uvatv_old(jj,1)+&
                                & a*gry*(deptotatv(i,jv1)*vmpred(i,jv1)+&
                                & 0.5*chin-2.0*jsign*c*zeta(i,jc))+&
                                & a*(jsign*dxu+vdy)+obrhs)
            Mv = Mvcoeff*0.5*(obc2uvatv(jj,1)+chin)
            Ma = -Mvcoeff*0.25*jsign*js*c*(2-js) !Ri
            Ma = Ma + 0.5*b*a*gry*Ma
            Ma = Ma + Mvcoeff*b*a*gry*jsign*c
            Ma1= Macoeff*0.5*b*a*gry
         CASE (12)
            chin = vdatobv(jj,2)+0.5*jsign*js*c*(zdatobv(jj,2)+&
                 & (2-js)*zeta(i,jc))
            obc2uvatv(jj,1) = b*(obc2uvatv_old(jj,1)+&
                               & a*gry*(deptotatv(i,jv1)*vmpred(i,jv1)+&
                               & 0.5*chin-2.0*jsign*c*zeta(i,jc))+&
                               & a*(jsign*dxu+vdy)+obrhs)
            Mv = Mvcoeff*0.5*(obc2uvatv(jj,1)+chin)
            Ma = Mvcoeff*0.5*jsign*js*c*( 0.5*(js-2) + a*gry )
            Ma = Ma + 0.5*b*a*gry*Ma 
            Ma = Ma + Mvcoeff*b*a*gry*jsign*c
            Ma1= Macoeff*0.5*b*a*gry
         CASE (13)
            obc2uvatv(jj,2) = obc2uvatv_old(jj,2) - a*(jsign*dxu+vdy) + obrhs
            obc2uvatv(jj,1) = b*(obc2uvatv_old(jj,1)+&
                               & a*gry*(deptotatv(i,jv1)*vmpred(i,jv1)+&
                               & 0.5*obc2uvatv(jj,2)-2.0*jsign*c*zeta(i,jc))+&
                               & a*(jsign*dxu+vdy)+obrhs)
            obc2uvatv(jj,2) = obc2uvatv_old(jj,2) + obrhs
            Mv = Mvcoeff*0.5*(obc2uvatv(jj,1)+obc2uvatv(jj,2))
            Ma = Mvcoeff*b*a*gry*jsign*c
            Ma1= Macoeff*0.5*b*a*gry
         CASE (14)                                     
            Mv = Mvcoeff*vdatobv(jj,2)*deptotatv(i,j)
         CASE (15:16)                                     
            Mv = Mvcoeff*vdatobv(jj,2)/delxatv(i,j)
         CASE (17)                                     
            Mv = Mvcoeff*(vdvel_old(i,j)-&
                       & delt2d*gaccatv(i,j)*deptotatv(i,j)*zdatobv(jj,2)+obrhs)
      END SELECT

!     ---insert in matrix
      Ma = Ma*theta_sur; Ma1 = Ma1*theta_sur
      mgvars(0)%rhs(i,jc) = mgvars(0)%rhs(i,jc) + Mv
      mgvars(0)%subdiags(i,jc,0) = mgvars(0)%subdiags(i,jc,0) + Ma - Ma1
      IF (soutobv(jj)) THEN
         mgvars(0)%subdiags(i,jc,2) = mgvars(0)%subdiags(i,jc,2) + Ma1
      ELSE
         mgvars(0)%subdiags(i,jc,-1) = mgvars(0)%subdiags(i,jc,-1) + Ma1
      ENDIF
   ENDDO jjloc_120

!
!2. Coarse grid levels
!---------------------
!

ELSEIF (lev.GT.0) THEN

!
!2.1 U-nodes
!-----------
!

   iiloc_210: DO iiloc=1,mgvars(lev)%nobuloc_ext
      i = mgvars(lev)%iobuloc(iiloc); j = mgvars(lev)%jobuloc(iiloc)
      ii = mgvars(lev)%indexobu(iiloc)
      IF ((mgvars(lev)%westobu(ii).AND.i.EQ.mgvars(lev)%ncloc+1).OR.&
        & (.NOT.mgvars(lev)%westobu(ii).AND.i.EQ.1)) THEN
         CYCLE iiloc_210
      ENDIF
      ic = MERGE(i,i-1,mgvars(lev)%westobu(ii))
      isign = MERGE(1,-1,mgvars(lev)%westobu(ii))
      iu1 = MERGE(i+1,i-1,mgvars(lev)%westobu(ii))
!     ---define parameters
      SELECT CASE (mgvars(lev)%ityp2dobu(ii))
         CASE (3,8:13)
            IF (mgvars(lev)%iloczobu(ii).EQ.1) THEN
               is = 2
            ELSEIF (mgvars(lev)%iloczobu(ii).EQ.2) THEN
               is = 1
            ENDIF
         END SELECT
      SELECT CASE (mgvars(lev)%ityp2dobu(ii))
         CASE (3,5,8:13)
            c = SQRT(mgvars(lev)%gaccatc(ic,j)*mgvars(lev)%deptotatc(ic,j))
            a = delt2d*c/mgvars(lev)%delxatc(ic,j)
      END SELECT
      SELECT CASE (mgvars(lev)%ityp2dobu(ii))
         CASE (3,9:13)
            grx = mgvars(lev)%delxatc(ic,j)/mgvars(lev)%delxatu(i,j)
      END SELECT
      SELECT CASE (mgvars(lev)%ityp2dobu(ii))
         CASE (11:13)
            b = 1.0/(1.0+1.5*a*grx)
      END SELECT

!     ---calculate entries for matrix system
      Mvcoeff = mgvars(lev)%delyatu(i,j)*isign
      Macoeff = Mvcoeff*mgvars(lev)%gaccatu(iu1,j)*&
          & delt3d/mgvars(lev)%delxatu(iu1,j)*mgvars(lev)%deptotatu(iu1,j)*isign
      Mv = 0.0; Ma = 0.0; Ma1= 0.0
      SELECT CASE (mgvars(lev)%ityp2dobu(ii))
         CASE (2)
            Ma1= Macoeff 
         CASE (3)
            Ma = Mvcoeff*isign*is*a*c*grx
         CASE (5) 
            Ma1= Macoeff*a/(1.0+a)
         CASE (8)
            Ma = Mvcoeff*0.5*isign*is*c
         CASE (9)
            Ma = Mvcoeff*isign*is*c*(0.5+a*grx)
         CASE (10)
            Ma = Mvcoeff*isign*c*(1.0+is*a*grx)
         CASE (11)
            Ma = -Mvcoeff*0.25*isign*is*c*(2-is) !Ri
            Ma = Ma + 0.5*b*a*grx*Ma 
            Ma = Ma + Mvcoeff*b*a*grx*isign*c
            Ma1= Macoeff*0.5*b*a*grx
         CASE (12)
            Ma = Mvcoeff*0.5*isign*is*c*(0.5*(is-2)+a*grx )
            Ma = Ma + 0.5*b*a*grx*Ma 
            Ma = Ma + Mvcoeff*b*a*grx*isign*c
            Ma1= Macoeff*0.5*b*a*grx
         CASE (13)
            Ma = Mvcoeff*b*a*grx*isign*c
            Ma1= Macoeff*0.5*b*a*grx
      END SELECT

!     ---insert in matrix
      Ma = Ma*theta_sur; Ma1 = Ma1*theta_sur
      mgvars(lev)%subdiags(ic,j,0) = mgvars(lev)%subdiags(ic,j,0) + Ma - Ma1
      IF (mgvars(lev)%westobu(ii)) THEN
         mgvars(lev)%subdiags(ic,j,1) = mgvars(lev)%subdiags(ic,j,1) + Ma1
      ELSE
         mgvars(lev)%subdiags(ic,j,-2) = mgvars(lev)%subdiags(ic,j,-2) + Ma1
      ENDIF

   ENDDO iiloc_210

!
!2.2 V-nodes
!-----------
!

   jjloc_220: DO jjloc=1,mgvars(lev)%nobvloc_ext
      i = mgvars(lev)%iobvloc(jjloc); j = mgvars(lev)%jobvloc(jjloc)
      jj = mgvars(lev)%indexobv(jjloc)
      IF ((mgvars(lev)%soutobv(jj).AND.j.EQ.mgvars(lev)%nrloc+1).OR.&
        & (.NOT.mgvars(lev)%soutobv(jj).AND.j.EQ.1)) THEN
         CYCLE jjloc_220
      ENDIF
      jc = MERGE(j,j-1,mgvars(lev)%soutobv(jj))
      jsign = MERGE(1,-1,mgvars(lev)%soutobv(jj))
      jv1 = MERGE(j+1,j-1,mgvars(lev)%soutobv(jj))
!     ---define parameters
      SELECT CASE (mgvars(lev)%ityp2dobv(jj))
         CASE (3,8:13)
            IF (mgvars(lev)%iloczobv(jj).EQ.1) THEN
               js = 2
            ELSEIF (mgvars(lev)%iloczobv(jj).EQ.2) THEN
               js = 1
            ENDIF
      END SELECT
      SELECT CASE (mgvars(lev)%ityp2dobv(jj))
         CASE (3,5,8:13)
            c = SQRT(mgvars(lev)%gaccatc(i,jc)*mgvars(lev)%deptotatc(i,jc))
            a = delt2d*c/mgvars(lev)%delyatc(i,jc)
      END SELECT
      SELECT CASE (mgvars(lev)%ityp2dobv(jj))
         CASE (3,9:13)
            gry = mgvars(lev)%delyatc(i,jc)/mgvars(lev)%delyatv(i,j)
      END SELECT
      SELECT CASE (mgvars(lev)%ityp2dobv(jj))
         CASE (11:13)
            b = 1.0/(1.0+1.5*a*gry)
      END SELECT

!     ---calculate entries for matrix system
      Mvcoeff = mgvars(lev)%delxatv(i,j)*jsign
      Macoeff = Mvcoeff*mgvars(lev)%gaccatv(i,jv1)*&
             & delt3d/mgvars(lev)%delyatv(i,jv1)*&
             & mgvars(lev)%deptotatv(i,jv1)*jsign
      Mv = 0.0; Ma = 0.0; Ma1= 0.0
      SELECT CASE (mgvars(lev)%ityp2dobv(jj))
         CASE (2)
            Ma1 = Macoeff
         CASE (3)
            Ma = Mvcoeff*jsign*js*a*c*gry
         CASE (5)
            Ma1= Macoeff*a/(1.0+a)
         CASE (8)
            Ma = Mvcoeff*0.5*jsign*js*c
         CASE (9)
            Ma = Mvcoeff*jsign*js*c*(0.5+a*gry)
         CASE (10)
            Ma = Mvcoeff*jsign*c*(1.0+js*a*gry) 
         CASE (11)
            Ma = -Mvcoeff*0.25*jsign*js*c*(2-js) !Ri
            Ma = Ma + 0.5*b*a*gry*Ma
            Ma = Ma + Mvcoeff*b*a*gry*jsign*c
            Ma1= Macoeff*0.5*b*a*gry
         CASE (12)
            Ma = Mvcoeff*0.5*jsign*js*c*(0.5*(js-2)+a*gry )
            Ma = Ma + 0.5*b*a*gry*Ma 
            Ma = Ma + Mvcoeff*b*a*gry*jsign*c
            Ma1= Macoeff*0.5*b*a*gry
         CASE (13)
            Ma = Mvcoeff*b*a*gry*jsign*c
            Ma1= Macoeff*0.5*b*a*gry
      END SELECT

!     ---insert in matrix
      Ma = Ma*theta_sur; Ma1 = Ma1*theta_sur
      mgvars(lev)%subdiags(i,jc,0) = mgvars(lev)%subdiags(i,jc,0) + Ma - Ma1
      IF (mgvars(lev)%soutobv(jj)) THEN
         mgvars(lev)%subdiags(i,jc,2) = mgvars(lev)%subdiags(i,jc,2) + Ma1
      ELSE
         mgvars(lev)%subdiags(i,jc,-1) = mgvars(lev)%subdiags(i,jc,-1) + Ma1
      ENDIF

   ENDDO jjloc_220

ENDIF

CALL log_timer_out(npcc,itm_bconds)


RETURN

END SUBROUTINE open_boundary_conds_2d_impl

!========================================================================

SUBROUTINE open_boundary_conds_3d(psiobux,psiobvy,obcdata,iddesc,itypobux,&
                                & itypobvy,iprofobux,iprofobvy,noprofs,nobux,&
                                & nobvy)
!************************************************************************
!
! *open_boundary_conds_3d* Obtain external profiles at open boundaries
!                          for 2-D or 3-D currents
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Conditions.f90  V2.2
!
! Description - 
!
! Reference -
!
! Calling program - current_corr
!
! External calls -
!
! Module calls - cfl_orlan, tridiag_vert_1d, Uarr_at_C, Varr_at_C
!
!************************************************************************
!
USE currents
USE depths
USE diffusion
USE fluxes
USE grid
USE gridpars
USE iopars
USE modids
USE obconds
USE physpars
USE switches
USE timepars
USE array_interp, ONLY: Uarr_at_C, Varr_at_C
USE nla_library, ONLY: tridiag_vert_1d
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: cfl_orlan

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, nobux, nobvy, noprofs
INTEGER, INTENT(IN), DIMENSION(nobux) :: iprofobux, itypobux
INTEGER, INTENT(IN), DIMENSION(nobvy) :: iprofobvy, itypobvy
REAL, INTENT(INOUT), DIMENSION(nobux,nz) :: psiobux
REAL, INTENT(INOUT), DIMENSION(nobvy,nz) :: psiobvy
REAL, INTENT(INOUT), DIMENSION(0:noprofs,nz) :: obcdata

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*psiobux*   REAL    Updated profile at external grid point (U- or X-node
!                    boundaries)
!*psiobvy*   REAL    Updated profile at external grid point (V- or Y-node
!                    boundaries)
!*obcdata*   REAL    Profiles from data file
!*iddesc*    INTEGER Data file id (io_3uvobc or io_3xyobc)
!*itypobux*  INTEGER Type of U- or X-open boundary condition
!*itypobvy*  INTEGER Type of V- or Y-open boundary condition
!*iprofobux* INTEGER Profile number at U- or X-open boundaries
!*iprofobvy* INTEGER Profile number at V- or Y-open boundaries
!*noprofs*   INTEGER Number of data profiles
!*nobux*     INTEGER Number of nodes at U- or X-open boundaries
!*nobvy*     INTEGER Number of nodes at V- or Y-open boundaries
!
!************************************************************************
!
!*Local variables
!
INTEGER :: i, ic1, ic2, ii, iiloc, iprof, isign, iu1, iu2, j, jc1, jc2, jj, &
         & jjloc, jprof, jsign, jv1, jv2, k, npcc
REAL :: cfl0, dxy, dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, fac, theta_vdif1, &
      & udev, uint0, uint1, uint2, vdev, vint0, vint1, vint2, xexp, ximp
REAL, DIMENSION(nz) :: uvelatc, vvelatc
REAL, DIMENSION(2,nz) :: uvals
REAL, DIMENSION(nz,4) :: tridcf
REAL, DIMENSION(2:nz) :: array1d1
REAL, DIMENSION(nz+1) :: array1d2, difflux


procname(pglev+1) = 'open_boundary_conds_3d'
CALL log_timer_in(npcc)



SELECT CASE (iddesc)

!
!1. Conditions at (U,V)-open boundaries
!--------------------------------------
!

   CASE (io_3uvobc)

      fac = cgravratio*delt3d

!
!
!1.1. U-nodes
!------------
!

      iiloc_110: DO iiloc=1,nobuloc

         i = iobuloc(iiloc); j = jobuloc(iiloc)
         ii = indexobu(iiloc)

!
!1.1.1 External profile
!----------------------
!

         IF (iprofobux(ii).GT.0) THEN

            iprof = iprofobux(ii)
            uint0 = SUM(obcdata(iprof,:)*delzatu(i,j,:))/deptotatu(i,j)
            uvel(i,j,:) = obcdata(iprof,:) - uint0

!
!1.1.2 First order zero gradient
!-------------------------------
!

         ELSEIF (itypobux(ii).EQ.0) THEN

            iu1 = MERGE(i+1,i-1,westobu(ii))
            dy1 = MERGE(1.0,delyatu(iu1,j)/delyatu(i,j),dYregX)
            WHERE (nodeatu(iu1,j,:).GT.1)
               uvel(i,j,:) = dy1*(delzatu(iu1,j,:)/delzatu(i,j,:))*uvel(iu1,j,:)
            ELSEWHERE
               uvel(i,j,:) = 0.0
            END WHERE

!
!1.1.3 Second order zero gradient
!--------------------------------
!

         ELSEIF (itypobux(ii).EQ.1) THEN
     
            ic1 = MERGE(i,i-1,westobu(ii))
            ic2 = MERGE(i+1,i-2,westobu(ii))
            iu1 = MERGE(i+1,i-1,westobu(ii))
            iu2 = MERGE(i+2,i-2,westobu(ii))

            k_113: DO k=1,nz
               IF (nodeatu(iu1,j,k).LE.1) THEN
                  uvel(i,j,k) = 0.0
               ELSEIF (nodeatu(iu2,j,k).NE.2) THEN
                  dy1 = MERGE(1.0,delyatu(iu1,j)/delyatu(i,j),dYregX)
                  dz1 = delzatu(iu1,j,k)/delzatu(i,j,k)
                  uvel(i,j,k) = dy1*dz1*uvel(iu1,j,k)
               ELSE
                  dx1 = MERGE(1.0,delxatc(ic1,j)/delxatc(ic2,j),dXregX)
                  dy1 = MERGE(1.0,delyatu(iu1,j)/delyatu(i,j),dYregX)
                  dy2 = MERGE(1.0,delyatc(ic1,j)/delyatu(ic2,j),dYregX)
                  dxy = dx1*dy2
                  dy3 = MERGE(1.0,delyatu(iu2,j)/delyatu(i,j),dYregX)
                  dz1 = delzatu(iu1,j,k)/delzatu(i,j,k)
                  dz2 = delzatu(iu2,j,k)/delzatu(i,j,k)
                  uvel(i,j,k) = dy1*dz1*(1.0+dxy)*uvel(iu1,j,k)-&
                              & dy3*dz2*dxy*uvel(iu2,j,k)
               ENDIF
            ENDDO k_113

!
!1.1.4 Local solution
!--------------------
!
           
         ELSEIF (itypobux(ii).EQ.2) THEN

            ic1 = MERGE(i,i-1,westobu(ii))
            iu1 = MERGE(i+1,i-1,westobu(ii))
            tridcf = 0.0

!
!1.1.4.1 Time derivative
!-----------------------
!

            tridcf(:,2) = 1.0
            tridcf(:,4) = (deptotatu_old(i,j)/deptotatu(i,j))*uvel(i,j,:)

!
!1.1.4.2 Surface and bottom stress
!---------------------------------
!

            tridcf(:,4) = tridcf(:,4) + &
                       & delt3d*(ubstresatu(i,j)-usstresatu(i,j))/deptotatu(i,j)

!
!1.1.4.3 Coriolis
!----------------
!

            vvelatc = 0.5*(delzatv(ic1,j,:)*vvel(ic1,j,:)/deptotatv(ic1,j)+&
                    & delzatv(ic1,j+1,:)*vvel(ic1,j+1,:)/deptotatv(ic1,j+1))
            tridcf(:,4) = tridcf(:,4) + delt3d*coriolatu(i,j)*deptotatu(i,j)*&
                     & (theta_cor*vvelatc+(1.0-theta_cor)*obc3uvatu(ii,:,2))/&
                     & delzatu(i,j,:)
            obc3uvatu(ii,:,2) = vvelatc

!
!1.1.4.4 Baroclinic pressure
!---------------------------
!

            WHERE (nodeatu(iu1,j,:).GT.1)
               tridcf(:,4) = tridcf(:,4) + delt3d*(delzatu(iu1,j,:)&
                  & /delzatu(i,j,:))*&
                  & (p3dbcgradatu(iu1,j,:)-p2dbcgradatu(iu1,j)/deptotatu(iu1,j))
            END WHERE

!
!1.1.4.5 Vertical diffusion
!--------------------------
!

            IF (iopt_vdif_coef.GT.0) THEN

!              ---time factors
               theta_vdif1 = 1.0-theta_vdif
               xexp = delt3d*theta_vdif1
               ximp = delt3d*theta_vdif

!              ---initialise fluxes
               difflux = 0.0
 
!              ---work space array
               array1d1 = vdifcoefmom(ic1,j,2:nz)/delzatuw(i,j,2:nz)

!              ---explicit fluxes
               IF (iopt_vdif_impl.NE.2) THEN
                  difflux(2:nz) = xexp*array1d1*&
                               & (uvel(i,j,2:nz)-uvel(i,j,1:nz-1))

!              ---diffusive term (explicit)
                  tridcf(:,4) = tridcf(:,4)+&
                              & (difflux(2:nz+1)-difflux(1:nz))/delzatu(i,j,:)
               ENDIF

!              ---diffusive term (implicit)
               IF (iopt_vdif_impl.NE.0) THEN
                  array1d2(2:nz) = ximp*array1d1/delzatu(i,j,2:nz)
                  tridcf(2:nz,1) = tridcf(2:nz,1) - array1d2(2:nz)
                  tridcf(2:nz,2) = tridcf(2:nz,2) + array1d2(2:nz)
                  array1d2(1:nz-1) = ximp*array1d1/delzatu(i,j,1:nz-1)
                  tridcf(1:nz-1,2) = tridcf(1:nz-1,2) + array1d2(1:nz-1)
                  tridcf(1:nz-1,3) = tridcf(1:nz-1,3) - array1d2(1:nz-1)
               ENDIF

!              ---surface boundary condition
               tridcf(nz,4) = tridcf(nz,4) + &
                            & delt3d*usstresatu(i,j)/delzatu(i,j,nz)

!              ---bottom boundary condition
               tridcf(1,4) = tridcf(1,4) - delt3d*ubstresatu(i,j)/delzatu(i,j,1)

            ENDIF

!
!1.1.4.6 Solve tridiagonal system
!--------------------------------
!

            CALL tridiag_vert_1d(tridcf,uvel(i,j,:),nz,1,.FALSE.)
            uvel(i,j,:) = uvel(i,j,:) - &
                        & SUM(delzatu(i,j,:)*uvel(i,j,:))/deptotatu(i,j)
            obc3uvatu(ii,:,1) = uvel(i,j,:)

!
!1.1.5 Radiation condition
!-------------------------
!

         ELSEIF (itypobux(ii).EQ.3) THEN

            ic1 = MERGE(i,i-1,westobu(ii))
            iu1 = MERGE(i+1,i-1,westobu(ii))
            isign = MERGE(1,-1,westobu(ii))
            uint0 = SUM(uvel_old(i,j,:)*delzatu(i,j,:))/deptotatu(i,j)
            uint1 = SUM(uvel_old(iu1,j,:)*delzatu(iu1,j,:))/deptotatu(iu1,j)
            cfl0 = isign*fac*SQRT(gaccatc(ic1,j)*deptotatc(ic1,j))&
                 & /delxatc(ic1,j)
            uvel(i,j,:) = (1.0-cfl0)*(uvel_old(i,j,:)-uint0)+&
                         & cfl0*(uvel_old(iu1,j,:)-uint1)

!
!1.1.6 Orlanski condition
!------------------------
!

         ELSEIF (itypobux(ii).EQ.4) THEN

            iu1 = MERGE(i+1,i-1,westobu(ii))
            iu2 = MERGE(i+2,i-2,westobu(ii))
            uint0 = SUM(uvel_old(i,j,:)*delzatu(i,j,:))/deptotatu(i,j)
            uint1 = SUM(uvel_old(iu1,j,:)*delzatu(iu1,j,:))/deptotatu(iu1,j)
            uint2 = SUM(uvel_old(iu2,j,:)*delzatu(iu2,j,:))/deptotatu(iu2,j)
            k_116: DO k=1,nz
               udev = uvel_old(iu1,j,k) - uint1
               cfl0 = cfl_orlan(udev,obc3uvatu(ii,k,1),obc3uvatu(ii,k,2))
               obc3uvatu(ii,k,1) = udev
               obc3uvatu(ii,k,2)  = uvel_old(iu2,j,k) - uint2
               uvel(i,j,k) = (1.0-cfl0)*(uvel_old(i,j,k)-uint0)+cfl0*udev
            ENDDO k_116

!
!1.1.7 Discharge condition
!-------------------------
!

         ELSEIF (itypobux(ii).EQ.5) THEN

            iprof = iprofobux(ii)
            uvel(i,j,:) = obcdata(iprof,:)/(delyatu(i,j)*delzatu(i,j,:)) - &
                        & udvel(i,j)/deptotatu(i,j)

         ENDIF
   
         psiobux(ii,:) = uvel(i,j,:)

      ENDDO iiloc_110

!
!1.2. V-nodes
!------------
!

      jjloc_120: DO jjloc=1,nobvloc

         i = iobvloc(jjloc); j = jobvloc(jjloc)
         jj = indexobv(jjloc)

!
!1.2.1 External profile
!----------------------
!

         IF (iprofobvy(jj).GT.0) THEN

            jprof = iprofobvy(jj)
            vint0 = SUM(obcdata(jprof,:)*delzatv(i,j,:))/deptotatv(i,j)
            vvel(i,j,:) = obcdata(jprof,:) - vint0

!
!1.2.2 First order zero gradient
!-------------------------------
!

         ELSEIF (itypobvy(jj).EQ.0) THEN

            jv1 = MERGE(j+1,j-1,soutobv(jj))
            dx1 = MERGE(1.0,delxatv(i,jv1)/delxatv(i,j),dXregY)
            WHERE (nodeatv(i,jv1,:).GT.1)
               vvel(i,j,:) = dx1*(delzatv(i,jv1,:)/delzatv(i,j,:))*vvel(i,jv1,:)
            ELSEWHERE
               vvel(i,j,:) = 0.0
            END WHERE

!
!1.2.3 Second order zero gradient
!--------------------------------
!

         ELSEIF (itypobvy(jj).EQ.1) THEN

            jc1 = MERGE(j,j-1,soutobv(jj))
            jc2 = MERGE(j+1,j-2,soutobv(jj))
            jv1 = MERGE(j+1,j-1,soutobv(jj))
            jv2 = MERGE(j+2,j-2,soutobv(jj))

            k_123: DO k=1,nz
               IF (nodeatv(i,jv1,k).LE.1) THEN
                  vvel(i,j,k) = 0.0
               ELSEIF (nodeatv(i,jv2,k).NE.2) THEN
                  dx1 = MERGE(1.0,delxatv(i,jv1)/delxatv(i,j),dXregY)
                  dz1 = delzatv(i,jv1,k)/delzatv(i,j,k)
                  vvel(i,j,k) = dx1*dz1*vvel(i,jv1,k)
               ELSE
                  dx1 = MERGE(1.0,delxatv(i,jv1)/delxatv(i,j),dXregY)
                  dy1 = MERGE(1.0,delyatc(i,jc1)/delyatc(i,jc2),dYregY)
                  dx2 = MERGE(1.0,delxatc(i,jc1)/delxatv(i,jc2),dXregY)
                  dxy = dy1*dx2
                  dx3 = MERGE(1.0,delxatv(i,jv2)/delxatv(i,j),dXregY)
                  dz1 = delzatv(i,jv1,k)/delzatv(i,j,k)
                  dz2 = delzatv(i,jv2,k)/delzatv(i,j,k)
                  vvel(i,j,k) = dx1*dz1*(1.0+dxy)*vvel(i,jv1,k)-&
                           & dx3*dz2*dxy*vvel(i,jv2,k)
               ENDIF
            ENDDO k_123

!
!1.2.4 Local solution
!--------------------
!
           
         ELSEIF (itypobvy(jj).EQ.2) THEN

            jc1 = MERGE(j,j-1,soutobv(jj))
            jv1 = MERGE(j+1,j-1,soutobv(jj))
            tridcf = 0.0

!
!1.2.4.1 Time derivative
!-----------------------
!

            tridcf(:,2) = 1.0
            tridcf(:,4) = (deptotatv_old(i,j)/deptotatv(i,j))*vvel(i,j,:)

!
!1.2.4.2 Surface and bottom stress
!---------------------------------
!

            tridcf(:,4) = tridcf(:,4) + &
                    & delt3d*(vbstresatv(i,j)-vsstresatv(i,j))/deptotatv(i,j)

!
!1.2.4.3 Coriolis
!----------------
!

            uvals(1,:) = (delzatu(i,jc1,:)/deptotatu(i,jc1))*&
                        & MERGE(uvel_old(i,jc1,:),uvel(i,jc1,:),&
                        & nodeatu(i,jc1,:).GT.2)
            uvals(2,:) = (delzatu(i+1,jc1,:)/deptotatu(i+1,jc1))*&
                        & MERGE(uvel_old(i+1,jc1,:),uvel(i+1,jc1,:),&
                        & nodeatu(i+1,jc1,:).GT.2)
            uvelatc = 0.5*(uvals(1,:)+uvals(2,:))
            tridcf(:,4) = tridcf(:,4) - delt3d*coriolatv(i,j)*deptotatv(i,j)*&
                     & (theta_cor*uvelatc+(1.0-theta_cor)*obc3uvatv(jj,:,2))/&
                     & delzatv(i,j,:)
            obc3uvatv(jj,:,2) = uvelatc

!
!1.2.4.4 Baroclinic pressure
!---------------------------
!

            WHERE (nodeatv(i,jv1,:).GT.1)
               tridcf(:,4) = tridcf(:,4) + delt3d*(delzatv(i,jv1,:)&
                  & /delzatv(i,j,:))*&
                  & (p3dbcgradatv(i,jv1,:)-p2dbcgradatv(i,jv1)/deptotatv(i,jv1))
            END WHERE

!
!1.2.4.5 Vertical diffusion
!--------------------------
!

            IF (iopt_vdif_coef.GT.0) THEN

!              ---time factors
               theta_vdif1 = 1.0-theta_vdif
               xexp = delt3d*theta_vdif1
               ximp = delt3d*theta_vdif

!              ---initialise fluxes
               difflux = 0.0
 
!              ---work space array
               array1d1 = vdifcoefmom(i,jc1,2:nz)/delzatvw(i,j,2:nz)

!              ---explicit fluxes
               IF (iopt_vdif_impl.NE.2) THEN
                  difflux(2:nz) = xexp*array1d1*&
                                & (vvel(i,j,2:nz)-vvel(i,j,1:nz-1))

!              ---diffusive term (explicit)
                  tridcf(:,4) = tridcf(:,4)+&
                              & (difflux(2:nz+1)-difflux(1:nz))/delzatv(i,j,:)
               ENDIF

!              ---diffusive term (implicit)
               IF (iopt_vdif_impl.NE.0) THEN
                  array1d2(2:nz) = ximp*array1d1/delzatv(i,j,2:nz)
                  tridcf(2:nz,1) = tridcf(2:nz,1) - array1d2(2:nz)
                  tridcf(2:nz,2) = tridcf(2:nz,2) + array1d2(2:nz)
                  array1d2(1:nz-1) = ximp*array1d1/delzatv(i,j,1:nz-1)
                  tridcf(1:nz-1,2) = tridcf(1:nz-1,2) + array1d2(1:nz-1)
                  tridcf(1:nz-1,3) = tridcf(1:nz-1,3) - array1d2(1:nz-1)
               ENDIF

!              ---surface boundary condition
               tridcf(nz,4) = tridcf(nz,4) + &
                            & delt3d*vsstresatv(i,j)/delzatv(i,j,nz)

!              ---bottom boundary condition
               tridcf(1,4) = tridcf(1,4) - delt3d*vbstresatv(i,j)/delzatv(i,j,1)

            ENDIF

!
!1.2.4.6 Solve tridiagonal system
!--------------------------------
!

            CALL tridiag_vert_1d(tridcf,vvel(i,j,:),nz,1,.FALSE.)
            vvel(i,j,:) = vvel(i,j,:) - &
                        & SUM(delzatv(i,j,:)*vvel(i,j,:))/deptotatv(i,j)
            obc3uvatv(jj,:,1) = vvel(i,j,:)

!
!1.2.5 Radiation condition
!-------------------------
!

         ELSEIF (itypobvy(jj).EQ.3) THEN

            jc1 = MERGE(j,j-1,soutobv(jj))
            jv1 = MERGE(j+1,j-1,soutobv(jj))
            jsign = MERGE(1,-1,soutobv(jj))
            vint0 = SUM(vvel_old(i,j,:)*delzatv(i,j,:))/deptotatv(i,j)
            vint1 = SUM(vvel_old(i,jv1,:)*delzatv(i,jv1,:))/deptotatv(i,jv1)
            cfl0 = jsign*fac*SQRT(gaccatc(i,jc1)*deptotatc(i,jc1))&
                 & /delyatc(i,jc1)
            vvel(i,j,:) = (1.0-cfl0)*(vvel_old(i,j,:)-vint0)+&
                        & cfl0*(vvel_old(i,jv1,:)-vint1)
!
!1.2.6 Orlanski condition
!------------------------
!

         ELSEIF (itypobvy(jj).EQ.4) THEN

            jv1 = MERGE(j+1,j-1,soutobv(jj))
            jv2 = MERGE(j+2,j-2,soutobv(jj))
            vint0 = SUM(vvel_old(i,j,:)*delzatv(i,j,:))/deptotatv(i,j)
            vint1 = SUM(vvel_old(i,jv1,:)*delzatv(i,jv1,:))/deptotatv(i,jv1)
            vint2 = SUM(vvel_old(i,jv2,:)*delzatv(i,jv2,:))/deptotatv(i,jv2)
            k_126: DO k=1,nz
               vdev = vvel_old(i,jv1,k) - vint1
               cfl0 = cfl_orlan(vdev,obc3uvatv(jj,k,1),obc3uvatv(jj,k,2))
               obc3uvatv(jj,k,1) = vdev
               obc3uvatv(jj,k,2)  = vvel_old(i,jv2,k) - vint2
               vvel(i,j,k) = (1.0-cfl0)*(vvel_old(i,j,k)-vint0)+cfl0*vdev
            ENDDO k_126

!
!1.2.7 Discharge condition
!--------------------------
!

         ELSEIF (itypobvy(jj).EQ.5) THEN

            jprof = iprofobvy(jj)
            vvel(i,j,:) = obcdata(jprof,:)/(delxatv(i,j)*delzatv(i,j,:)) - &
                        & vdvel(i,j)/deptotatv(i,j)

         ENDIF

         psiobvy(jj,:) = vvel(i,j,:)

      ENDDO jjloc_120


!
!2. Tangential conditions
!------------------------
!

   CASE (io_3xyobc)

!
!2.1. X-nodes
!------------
!

      iiloc_210: DO iiloc=1,nobxloc

         i = iobxloc(iiloc); j = jobxloc(iiloc)
         ii = indexobx(iiloc)
         IF (itypobux(ii).EQ.2.AND.iprofobux(ii).GT.0) THEN
            iprof = iprofobux(ii)
            vint0 = SUM(obcdata(iprof,:)*delzatuv(i,j,:))/deptotatuv(i,j)
            psiobux(ii,:) = obcdata(iprof,:) + vdatobx(ii,2)/deptotatuv(i,j) - &
                          & vint0
         ENDIF

      ENDDO iiloc_210

!
!2.2. Y-nodes
!------------
!

      jjloc_220: DO jjloc=1,nobyloc

         i = iobyloc(jjloc); j = jobyloc(jjloc)
         jj = indexoby(jjloc)
         IF (itypobvy(jj).EQ.2.AND.iprofobvy(jj).GT.0) THEN
            jprof = iprofobvy(jj)
            uint0 = SUM(obcdata(jprof,:)*delzatuv(i,j,:))/deptotatuv(i,j)
            psiobvy(jj,:) = obcdata(jprof,:) + udatoby(jj,2)/deptotatuv(i,j) - &
                          & uint0
         ENDIF

      ENDDO jjloc_220

END SELECT

CALL log_timer_out(npcc,itm_bconds)


RETURN

END SUBROUTINE open_boundary_conds_3d

!========================================================================

SUBROUTINE open_boundary_conds_prof(psi,psiobu,psiobv,obcdata,obcpsiatu,&
                                  & obcpsiatv,itypobu,itypobv,iprofobu,&
                                  & iprofobv,noprofs,novars)
!************************************************************************
!
! *open_boundary_conds_prof* Obtain external profiles at open boundaries
!                            for 3-D scalars
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Conditions.f90  V2.7.1
!
! Description - 
!
! Reference -
!
! Calling program - salinity_equation, sediment_advdiff, temperature_equation 
!
! External calls -
!
! Module calls - cfl_orlan
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE obconds
USE physpars
USE timepars
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: cfl_orlan

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: noprofs, novars
INTEGER, INTENT(INOUT), DIMENSION(nobu) :: itypobu
INTEGER, INTENT(INOUT), DIMENSION(nobv) :: itypobv
INTEGER, INTENT(INOUT), DIMENSION(nobu,novars) :: iprofobu
INTEGER, INTENT(INOUT), DIMENSION(nobv,novars) :: iprofobv
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,&
                          & 1-nhalo:nrloc+nhalo,nz,novars) :: psi
REAL, INTENT(INOUT), DIMENSION(nobu,nz,novars) :: psiobu
REAL, INTENT(INOUT), DIMENSION(nobv,nz,novars) :: psiobv
REAL, INTENT(IN), DIMENSION(0:noprofs,nz,novars) :: obcdata
REAL, INTENT(INOUT), DIMENSION(nobu,nz,0:2,novars) :: obcpsiatu
REAL, INTENT(INOUT), DIMENSION(nobv,nz,0:2,novars) :: obcpsiatv

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*psi*       REAL    Full model array
!*psiobu*    REAL    Updated profile at external grid point (U-node boundaries)
!*psiobv*    REAL    Updated profile at external grid point (V-node boundaries)
!*obcdata*   REAL    Profiles from data file
!*obcpsiatu* REAL    Storage array for radiation o.b.c at U-nodes
!*obcpsiatv* REAL    Storage array for radiation o.b.c at V-node
!*itypobu*   INTEGER Type of o.b.c. at U-nodes
!*itypobv*   INTEGER Type of o.b.c. at V-nodes
!*iprofobu*  INTEGER Profile number at U-nodes
!*iprofobv*  INTEGER Profile number at V-nodes
!*noprofs*   INTEGER Number of data profiles
!*novars*    INTEGER Number of variables
!
!************************************************************************
!
!*Local variables
!
INTEGER :: i, ic1, ic2, ii, iiloc, iprof, ivar, j, jc1, jc2, jj, jjloc, jprof, &
         & k, npcc
REAL :: cfl, fac


procname(pglev+1) = 'open_boundary_conds_prof'
CALL log_timer_in(npcc)

fac = cgravratio*delt3d

!
!1. U-nodes
!----------
!

ivar_100: DO ivar=1,novars
iiloc_100: DO iiloc=1,nobuloc
   i = iobuloc(iiloc); j = jobuloc(iiloc)
   ii = indexobu(iiloc)
   iprof = iprofobu(ii,ivar)
   IF (iprof.GT.0) THEN
      WHERE (ABS(obcdata(iprof,:,ivar)-float_fill).GT.float_fill_eps)
         psiobu(ii,:,ivar) = obcdata(iprof,:,ivar)
      END WHERE
   ENDIF

   SELECT CASE (itypobu(ii))

!     ---zero gradient condition
      CASE(0)
         WHERE (ABS(obcdata(iprof,:,ivar)-float_fill).LE.float_fill_eps)
            psiobu(ii,:,ivar) = float_fill
         END WHERE

!     ---radiation condition using internal wave speed
      CASE (1)
         ic1 = MERGE(i,i-1,westobu(ii))
         cfl = fac*SQRT(gaccatu(i,j)*deptotatu(i,j))/delxatu(i,j)
         WHERE (ABS(obcdata(iprof,:,ivar)-float_fill).LE.float_fill_eps)
            psiobu(ii,:,ivar) = (1.0-cfl)*psiobu(ii,:,ivar)+&
                               & cfl*psi(ic1,j,:,ivar)
         END WHERE

!     ---Orlanski radiation condition
      CASE (2)
         ic1 = MERGE(i,i-1,westobu(ii))
         ic2 = MERGE(i+1,i-2,westobu(ii))
         k_120: DO k=1,nz
            IF (ABS(obcdata(iprof,k,ivar)-float_fill).LE.float_fill_eps) THEN
               cfl = cfl_orlan(psi(ic1,j,k,ivar),obcpsiatu(ii,k,1,ivar),&
                             & obcpsiatu(ii,k,ivar,2))
               obcpsiatu(ii,k,1,ivar) = psi(ic1,j,k,ivar)
               obcpsiatu(ii,k,2,ivar) = psi(ic2,j,k,ivar)
               psiobu(ii,k,ivar) = (1.0-cfl)*obcpsiatu(ii,k,0,ivar)+&
                                 & cfl*psi(ic1,j,k,ivar)
               obcpsiatu(ii,k,0,ivar) = psiobu(ii,k,ivar)
            ENDIF
         ENDDO k_120

!     ---Thatcher-Harleman condition
      CASE (3)
         ic1 = MERGE(i,i-1,westobu(ii))
         k_130: DO k=1,nz
            IF (floutobu(ii,k).GT.return_time(k)) THEN
               psiobu(ii,k,ivar) = obcdata(iprof,k,ivar)
            ELSEIF (floutobu(ii,k).GT.0.0) THEN
               psiobu(ii,k,ivar) = obcpsiatu(ii,k,0,ivar) + &
                         & 0.5*(obcdata(iprof,k,ivar)-obcpsiatu(ii,k,0,ivar))*&
                         & (COS(pi*(1.0-floutobu(ii,k)/return_time(k)))+1.0)
            ELSE
               obcpsiatu(ii,k,0,ivar) = psi(ic1,j,k,ivar)
            ENDIF
         ENDDO k_130
 
   END SELECT

ENDDO iiloc_100
ENDDO ivar_100

!
!2. V-nodes
!----------
!

ivar_200: DO ivar=1,novars
jjloc_200: DO jjloc=1,nobvloc
   i = iobvloc(jjloc); j = jobvloc(jjloc)
   jj = indexobv(jjloc)
   jprof = iprofobv(jj,ivar)
   IF (jprof.GT.0) THEN
      WHERE (ABS(obcdata(jprof,:,ivar)-float_fill).GT.float_fill_eps)
         psiobv(jj,:,ivar) = obcdata(jprof,:,ivar)
      END WHERE
   ENDIF

   SELECT CASE (itypobv(jj))

!     ---zero gradient condition
      CASE(0)
         WHERE (ABS(obcdata(jprof,:,ivar)-float_fill).LE.float_fill_eps)
            psiobv(jj,:,ivar) = float_fill
         END WHERE

!     ---radiation condition using internal wave speed
      CASE (1)
         jc1 = MERGE(j,j-1,soutobv(jj))
         cfl = fac*SQRT(gaccatv(i,j)*deptotatv(i,j))/delxatv(i,j)
         WHERE (ABS(obcdata(jprof,:,ivar)-float_fill).LE.float_fill_eps)
            psiobv(jj,:,ivar) = (1.0-cfl)*psiobv(jj,:,ivar)+&
                              & cfl*psi(i,jc1,:,ivar)
         END WHERE

!     ---Orlanski radiation condition
      CASE (2)
         jc1 = MERGE(j,j-1,soutobv(jj))
         jc2 = MERGE(j+1,j-2,soutobv(jj))
         k_220: DO k=1,nz
            IF (ABS(obcdata(jprof,k,ivar)-float_fill).LE.float_fill_eps) THEN
               cfl = cfl_orlan(psi(i,jc1,k,ivar),obcpsiatv(jj,k,1,ivar),&
                             & obcpsiatv(jj,k,2,ivar))
               obcpsiatv(jj,k,1,ivar) = psi(i,jc1,k,ivar)
               obcpsiatv(jj,k,2,ivar) = psi(i,jc2,k,ivar)
               psiobv(jj,k,ivar) = (1.0-cfl)*obcpsiatv(jj,k,0,ivar)+&
                                 & cfl*psi(i,jc1,k,ivar)
               obcpsiatv(jj,k,0,ivar) = psiobv(jj,k,ivar)
            ENDIF
         ENDDO k_220

!     ---Thatcher-Harleman condition
      CASE (3)
         jc1 = MERGE(J,J-1,soutobv(jj))
         k_230: DO k=1,nz
            IF (floutobv(jj,k).GT.return_time(k)) THEN
               psiobv(jj,k,ivar) = obcdata(iprof,k,ivar)
            ELSEIF (floutobv(jj,k).GT.0.0) THEN
               psiobv(jj,k,ivar) = obcpsiatv(jj,k,0,ivar) + &
                         & 0.5*(obcdata(iprof,k,ivar)-obcpsiatu(ii,k,0,ivar))*&
                         & (COS(pi*(1.0-floutobv(jj,k)/return_time(k)))+1.0)
            ELSE
               obcpsiatv(jj,k,0,ivar) = psi(i,jc1,k,ivar)
            ENDIF
         ENDDO k_230

   END SELECT

ENDDO jjloc_200
ENDDO ivar_200

CALL log_timer_out(npcc,itm_bconds)


RETURN

END SUBROUTINE open_boundary_conds_prof

!========================================================================

SUBROUTINE open_boundary_outflow_time
!************************************************************************
!
! *open_boundary_outflow_time* Update the time of last outflow at open
!                              boundaries
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Conditions.f90  V2.7
!
! Description - routine is called in connection with the Thatcher-Harleman
!               open boundary condition for scalars
!
! Reference -
!
! Calling program - hydrodynamic_equations
!
! External calls -
!
! Module calls - cfl_orlan
!
!************************************************************************
!
USE currents
USE grid
USE gridpars
USE iopars
USE obconds
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, ii, iiloc, j, jj, jjloc, npcc


procname(pglev+1) = 'open_boundary_conds_prof'
CALL log_timer_in(npcc)

!
!1.U-open boundaries
!-------------------
!

ii_110: DO iiloc=1,nobuloc
   i = iobuloc(iiloc); j = jobuloc(iiloc)
   ii = indexobu(iiloc)
   IF (westobu(ii)) THEN
      WHERE (uvel(i,j,:).LT.0.0)
         floutobu(ii,:)  = floutobu(ii,:) + delt3d
      ELSEWHERE
         floutobu(ii,:) = 0.0
      END WHERE
   ELSE
      WHERE (uvel(i,j,:).GT.0.0)
         floutobu(ii,:)  = floutobu(ii,:) + delt3d
      ELSEWHERE
         floutobu(ii,:) = 0.0
      END WHERE
   ENDIF
ENDDO ii_110

!
!2.V-open boundaries
!-------------------
!

jj_210: DO jjloc=1,nobvloc
   i = iobvloc(jjloc); j = jobvloc(jjloc)
   jj = indexobv(jjloc)
   IF (soutobv(jj)) THEN
      WHERE (vvel(i,j,:).LT.0.0)
         floutobv(jj,:)  = floutobv(jj,:) + delt3d
      ELSEWHERE
         floutobv(jj,:) = 0.0
      END WHERE
   ELSE
      WHERE (vvel(i,j,:).GT.0.0)
         floutobv(jj,:)  = floutobv(jj,:) + delt3d
      ELSEWHERE
         floutobv(jj,:) = 0.0
      END WHERE
   ENDIF
ENDDO jj_210

CALL log_timer_out(npcc,itm_bconds)


RETURN

END SUBROUTINE open_boundary_outflow_time
