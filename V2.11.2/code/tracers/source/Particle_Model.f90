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
! *Particle_Model* Routines for the  particle drift model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.11
!
! $Date: 2018-01-19 09:43:56 +0100 (Fri, 19 Jan 2018) $
!
! $Revision: 1074 $
!
! Description -
!
! Reference -
!
! Routines - horizontal_part_drift, particle_concentrations, particle_drift,
!            particle_model, vertical_particle_drift
!
!************************************************************************
!

  
!========================================================================

SUBROUTINE horizontal_particle_drift
!************************************************************************
!
! *horizontal_particle_drift* horizontal particle drift
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.11
!
! $Date: 2018-01-19 09:43:56 +0100 (Fri, 19 Jan 2018) $
!
! $Revision: 1074 $
!
! Description -
!
! Reference -
!
! Calling program - particle_drift
!
! External calls -
!
! Module calls - Cint_at_p, complex_polar, error_abort, hrel_coords_curv_rel,
!                hrel_coords_rect_nonunif_rel, hrel_coords_rect_unif_rel,
!                rng_standard_uniform, Uint_at_p, Vint_at_p
!
!************************************************************************
!
USE currents
USE grid  
USE gridpars  
USE iopars
USE meteo
USE partpars
USE partswitches
USE partvars
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_abort
USE math_library, ONLY: complex_polar
USE grid_interp, ONLY: hrel_coords_curv_rel, hrel_coords_rect_nonunif_rel, &
                     & hrel_coords_rect_unif_rel
USE particle_routines, ONLY: Cint_at_p, Uint_at_p, Vint_at_p
USE rng_library, ONLY: rng_standard_uniform
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, istat, j, n, p, state
REAL :: delxdat, delydat, displx, disply, hconv, hfacx, hfacy, lway, &
      & theta, uvel_adv, uvelp, uvelp_new, uvelp_old, uwindp, uwindp_new, &
      & uwindp_old, vvel_adv, vvelp, vvelp_new, vvelp_old, vwindp, vwindp_new, &
      & vwindp_old, winddir, windmag, x, xran, y
REAL (KIND=kndrlong) :: xpos, ypos, zpos
REAL, DIMENSION(4) :: rkcoef, rktime, uvel_rk, vvel_rk

  
procname(pglev+1) = 'horizontal_particle_drift'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!
!---horizontal diffusion
IF (iopt_part_hdif.GT.0) THEN
   hfacx = SQRT(4.0*xdifpart_cst/delt3d); hfacy = SQRT(4.0*ydifpart_cst/delt3d)
ENDIF

!---runge-kutta integration factors
IF (iopt_part_hadv.EQ.2) THEN
   rkcoef = (/1.0,2.0,2.0,1.0/)
   rkcoef = rkcoef/6.0
   rktime = (/0.0,0.5,0.5,1.0/)
ENDIF

!---conversion factor for spherical coordinates
hconv = degtorad_d*Rearth

!---lee-way cross-wind direction change factor
IF (iopt_part_leeway.EQ.0) THEN
   lway = 1.0
ELSE
   CALL rng_standard_uniform(xran,prangenlw)
   lway = SIGN(1.0,xran-0.04*(delt3d/3600.0))
ENDIF

!---displacements
part%displx = 0.0; part%disply = 0.0

!---work space
delxdat = surfacegrids(igrd_model,1)%delxdat
delydat = surfacegrids(igrd_model,1)%delydat

!
!2. Horizontal diffusive current
!-------------------------------
!

IF (iopt_part_hdif.GT.0) THEN

   p_210: DO p=1,nopart
      IF (part(p)%drift_state.EQ.2) THEN
         part(p)%displx = deltpart*hfacx*ranmagpart(p)*COS(randirpart(p))
         part(p)%disply = deltpart*hfacy*ranmagpart(p)*SIN(randirpart(p))
      ENDIF
   ENDDO p_210

ENDIF
   
!
!3. Horizontal advective current
!--------------------------------
!
!3.1 First order Euler time integration
!--------------------------------------
!   

IF (iopt_part_hadv.EQ.1) THEN
   p_310: DO p=1,nopart
      IF (part(p)%drift_state.EQ.2) THEN

!      
!3.1.1 Initialise parameters
!---------------------------
!         

         i = part(p)%icoord; j = part(p)%jcoord
         x = part(p)%xcoord; y = part(p)%ycoord
         zpos = part(p)%zpos
         state = part(p)%state

!         
!3.1.2 Current driven drift
!----------------------------         
!
!        ---interpolate currents at particle locations
         IF (state.EQ.1) THEN
            uvelp = Uint_at_p(uvel(i:i+1,j-1:j+1,nz),i,j,x,y)
            vvelp = Vint_at_p(vvel(i-1:i+1,j:j+1,nz),i,j,x,y)
         ELSE
            uvelp = Uint_at_p(uvel(i:i+1,j-1:j+1,:),i,j,x,y,zpos)
            vvelp = Vint_at_p(vvel(i-1:i+1,j:j+1,:),i,j,x,y,zpos)
         ENDIF

!        ---drift current
         uvel_adv = cdrift_fac*uvelp
         vvel_adv = cdrift_fac*vvelp
         
!
!3.1.3 Wind driven drift
!-------------------------
!

         IF (iopt_part_wind.EQ.1.AND.state.EQ.1) THEN

!           ---interpolate wind currents at particle location
            uwindp = Cint_at_p(uwindatc(i-1:i+1,j-1:j+1),i,j,x,y)
            vwindp = Cint_at_p(vwindatc(i-1:i+1,j-1:j+1),i,j,x,y)
            
!           ---drift angle            
            CALL complex_polar(uwindp,vwindp,xamp=windmag,xpha=winddir)
            theta = lway*degtorad*(MAX(0.0,wdrift_angle1-wdrift_angle2*&
                                 & SQRT(windmag)))
            theta = MOD(theta,twopi)
            
!           ---add wind-driven drift
            uvel_adv = uvel_adv + wdrift_slope*windmag*COS(winddir-theta)
            vvel_adv = vvel_adv + wdrift_slope*windmag*SIN(winddir-theta)

         ENDIF

         part(p)%displx = part(p)%displx + deltpart*uvel_adv
         part(p)%disply = part(p)%disply + deltpart*vvel_adv

      ENDIF
         
   ENDDO p_310
            
!
!3.2 Fourth order Runge-Kutta
!----------------------------
!   

ELSEIF (iopt_part_hadv.EQ.2) THEN
            
   p_320: DO p=1,nopart
      IF (part(p)%drift_state.EQ.2) THEN

!      
!3.2.1 Initialise parameters
!---------------------------
!         

         i = part(p)%icoord; j = part(p)%jcoord
         x = part(p)%xcoord; y = part(p)%ycoord
         xpos = part(p)%xpos; ypos = part(p)%ypos; zpos = part(p)%zpos
         displx = 0.0; disply = 0.0
         state = part(p)%state

!
!3.2.2 Runge-Kutta algorithm
!----------------------------
!         

         n_322: DO n=1,4

!
!3.2.2.1 Update locations            
!--------------------------
!

            IF (n.GT.1) THEN
               SELECT CASE (iopt_grid_htype)
                  CASE (1)
                     CALL hrel_coords_rect_unif_rel(xpos,ypos,displx,disply,&
                                     & delxdat,delydat,(/1,1/),(/nc,nr/),istat,&
                                     & i,j,xcoord=x,ycoord=y)
                  CASE (2)
                     CALL hrel_coords_rect_nonunif_rel(xpos,ypos,displx,&
                                     & disply,gxcoord(1:nc,1),gycoord(1,1:nr),&
                                     & nc,nr,(/1,1/),(/nc,nr/),istat,i,j,&
                                     & xcoord=x,ycoord=y)
                  CASE (3)
                     CALL hrel_coords_curv_rel(xpos,ypos,displx,disply,&
                                     & gxcoord(1:nc,1:nr),gycoord(1:nc,1:nr),&
                                     & nc,nr,(/1,1/),(/nc,nr/),nsize_sb,istat,&
                                     & 2,i,j,xcoord=x,ycoord=y)
               END SELECT
               IF (istat.EQ.0.AND.(.NOT.maskatc(i,j))) THEN
                  part(p)%drift_state = 3
                  CYCLE p_320
               ELSEIF (istat.EQ.1) THEN
                  part(p)%drift_state = 4
                  CYCLE p_320
               ENDIF
            ENDIF
            
!            
!3.2.2.2 Current driven drift
!------------------------------
!
!           ---interpolate currents at particle locations and old/new time step
            IF (state.EQ.1) THEN
               uvelp_old = Uint_at_p(uvel_old(i:i+1,j-1:j+1,nz),i,j,x,y)
               vvelp_old = Vint_at_p(vvel_old(i-1:i+1,j:j+1,nz),i,j,x,y)
               uvelp_new = Uint_at_p(uvel(i:i+1,j-1:j+1,nz),i,j,x,y)
               vvelp_new = Vint_at_p(vvel(i-1:i+1,j:j+1,nz),i,j,x,y)
            ELSE
               uvelp_old = Uint_at_p(uvel_old(i:i+1,j-1:j+1,:),i,j,x,y,zpos)
               vvelp_old = Vint_at_p(vvel_old(i-1:i+1,j:j+1,:),i,j,x,y,zpos)
               uvelp_new = Uint_at_p(uvel(i:i+1,j-1:j+1,:),i,j,x,y,zpos)
               vvelp_new = Vint_at_p(vvel(i-1:i+1,j:j+1,:),i,j,x,y,zpos)
            ENDIF

!           ---drift current
            uvel_rk(n) = cdrift_fac*(uvelp_old+rktime(n)*(uvelp_new-uvelp_old))
            vvel_rk(n) = cdrift_fac*(vvelp_old+rktime(n)*(vvelp_new-vvelp_old))

!
!3.2.2.3 Wind driven drift
!-------------------------
!

            IF (iopt_part_wind.EQ.1.AND.state.EQ.1) THEN

!              ---interpolate wind at particle locations and old time step
               uwindp_old = Cint_at_p(uwindatc_old(i-1:i+1,j-1:j+1),i,j,x,y)
               vwindp_old = Cint_at_p(vwindatc_old(i-1:i+1,j-1:j+1),i,j,x,y)

!              ---interpolate wind at particle locations and new time step
               uwindp_new = Cint_at_p(uwindatc(i-1:i+1,j-1:j+1),i,j,x,y)
               vwindp_new = Cint_at_p(vwindatc(i-1:i+1,j-1:j+1),i,j,x,y)

!              ---at current Runge-Kutta step
               uwindp = uwindp_old+rktime(n)*(uwindp_new-uwindp_old)
               vwindp = vwindp_old+rktime(n)*(vwindp_new-vwindp_old)

!              ---drift angle
               CALL complex_polar(uwindp,vwindp,xamp=windmag,xpha=winddir)
               theta = lway*degtorad*&
                     & (MAX(0.0,wdrift_angle1-wdrift_angle2*SQRT(windmag)))

!              ---add wind-driven drift            
               uvel_rk(n) = uvel_rk(n) + wdrift_slope*windmag*COS(winddir-theta)
               vvel_rk(n) = vvel_rk(n) + wdrift_slope*windmag*SIN(winddir-theta)

            ENDIF

!
!3.2.2.4 Displacements
!-----------------------
!

            IF (n.LT.4) THEN
               displx = deltpart*rktime(n+1)*uvel_rk(n)
               disply = deltpart*rktime(n+1)*vvel_rk(n)
               IF (iopt_grid_sph.EQ.1) THEN
                  displx = displx/(hconv*COS(degtorad_d*ypos))
                  disply = disply/hconv
               ENDIF
            ENDIF
            
         ENDDO n_322

!
!3.2.3 Drift current at new time step
!----------------------------------      
!

         uvel_adv = SUM(rkcoef*uvel_rk)
         vvel_adv = SUM(rkcoef*vvel_rk)

!
!3.2.4 Drift displacements
!-----------------------
!

         part(p)%displx = part(p)%displx + deltpart*uvel_adv
         part(p)%disply = part(p)%disply + deltpart*vvel_adv

      ENDIF
      
   ENDDO p_320

ENDIF

CALL error_abort('horizontal_particle_drift',ierrno_runval)

!
!4. Correct displacements for spherical grid
!-------------------------------------------
!

IF (iopt_grid_sph.EQ.1) THEN
   p_410: DO p=1,nopart
      IF (part(p)%drift_state.EQ.2) THEN
         part(p)%displx = part(p)%displx/(hconv*COS(degtorad_d*part(p)%ypos))
         part(p)%disply = part(p)%disply/hconv
      ENDIF
   ENDDO p_410
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE horizontal_particle_drift

!========================================================================

SUBROUTINE particle_concentrations
!************************************************************************
!
! *particle_concentration* convert particle densities into Eulerian
!                          concentrations 
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.11
!
! $Date: 2018-01-19 09:43:56 +0100 (Fri, 19 Jan 2018) $
!
! $Revision: 1074 $
!
! Description -
!
! Reference -
!
! Calling program - coherens_main, particle_model
!
! External calls -
!
! Module calls - distribute_mod, error_alloc
!
!************************************************************************
!
USE gridpars  
USE iopars
USE paralpars
USE partids
USE partpars
USE partswitches
USE partvars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: distribute_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j, k, l, npcc, p
INTEGER, DIMENSION(4) :: lbounds, nhdist
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: realglb


IF (iopt_part_conc.EQ.0) RETURN

procname(pglev+1) = 'particle_concentrations'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

IF (iopt_part_model.EQ.2) THEN
   ALLOCATE (realglb(nc,nr,nz,nolabels),STAT=errstat)
   CALL error_alloc('realglb',4,(/nc,nr,nz,nolabels/),kndrtype)
   realglb = 0.0
ENDIF

!
!2. Update Eulerian concentrations
!---------------------------------
!

IF (iopt_part_model.NE.2.OR.&
 & (iopt_part_model.EQ.2.AND.modelid.EQ.modelidpart)) THEN

!
!2.1 Number densities
!--------------------
!

   conc = 0.0; conc_tot = 0.0
   p_210: DO p=1,nopart
      IF (part(p)%drift_state.EQ.2.OR.part(p)%drift_state.EQ.3) THEN
         i = part(p)%icoord; j = part(p)%jcoord
         k = part(p)%kcoord; l = part(p)%label
         conc(i,j,k,l) = conc(i,j,k,l) + 1.0
         conc_tot(i,j,k) = conc_tot(i,j,k) + 1.0
      ENDIF
   ENDDO p_210
   
!
!2.2 Age concentrations
!---------------------   
!

   ageconc = 0.0; ageconc_tot = 0.0
   p_221: DO p=1,nopart
      IF (part(p)%drift_state.EQ.2.OR.part(p)%drift_state.EQ.3) THEN
         i = part(p)%icoord; j = part(p)%jcoord
         k = part(p)%kcoord; l = part(p)%label
         IF (conc(i,j,k,l).GT.0.0) THEN
            ageconc(i,j,k,l) = ageconc(i,j,k,l) + part(p)%age/conc(i,j,k,l)
         ENDIF
         IF (conc_tot(i,j,k).GT.0.0) THEN
            ageconc_tot(i,j,k) = ageconc_tot(i,j,k) + &
                               & part(p)%age/conc_tot(i,j,k)
         ENDIF
      ENDIF
   ENDDO p_221

!
!2.3 Convert into volume concentrations
!--------------------------------------
!

   IF (iopt_part_conc.EQ.2) THEN

      l_231: DO l=1,nolabels
         conc(:,:,:,l) = conc(:,:,:,l)*volumeconc(l)
      ENDDO l_231
      conc_tot = 0.0
      l_232: DO l=1,nolabels
         conc_tot = conc_tot + conc(:,:,:,l)
      END DO l_232

!
!2.4 Convert into mass concentrations
!------------------------------------
!

   ELSEIF (iopt_part_conc.EQ.3) THEN

      l_241: DO l=1,nolabels
         conc(:,:,:,l) = conc(:,:,:,l)*massconc(l)
      ENDDO l_241
      conc_tot = 0.0
      l_242: DO l=1,nolabels
         conc_tot = conc_tot + conc(:,:,:,l)
      ENDDO l_242
   ENDIF

ENDIF

!
!3. Distribute data to COHERENS process domains
!-----------------------------------------------
!

IF (iopt_part_model.EQ.2) THEN
   lbounds = 1; nhdist = 0
   IF (modelid.EQ.modelidpart) realglb = conc
   CALL distribute_mod(realglb,conc,lbounds,nhdist,iarr_part_conc,0.0,&
                     & shared=.FALSE.,idroot=idmasterpart,comm=comm_world_MPI)
   IF (modelid.EQ.modelidpart) realglb = ageconc
   CALL distribute_mod(realglb,ageconc,lbounds,nhdist,iarr_part_ageconc,0.0,&
                     & shared=.FALSE.,idroot=idmasterpart,comm=comm_world_MPI)
   IF (modelid.EQ.modelidpart) realglb(:,:,:,1) = conc_tot
   CALL distribute_mod(realglb(:,:,:,1),conc_tot,lbounds(1:3),nhdist,&
                     & iarr_part_conc_tot,0.0,shared=.FALSE.,&
                     & idroot=idmasterpart,comm=comm_world_MPI)
   IF (modelid.EQ.modelidpart) realglb(:,:,:,1) = ageconc_tot
   CALL distribute_mod(realglb(:,:,:,1),ageconc_tot,lbounds(1:3),nhdist,&
                     & iarr_part_ageconc_tot,0.0,shared=.FALSE.,&
                     & idroot=idmasterpart,comm=comm_world_MPI)
ENDIF

!
!4. Deallocate
!-------------
!

IF (iopt_part_model.EQ.2) DEALLOCATE (realglb)

IF (modelid.EQ.modelidpart) THEN
   CALL log_timer_out()
ELSE
   CALL log_timer_out(npcc,itm_part)
ENDIF


RETURN

END SUBROUTINE particle_concentrations

!========================================================================

SUBROUTINE particle_drift
!************************************************************************
!
! *particle_drift* Particle drift model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.11
!
! $Date: 2018-01-19 09:43:56 +0100 (Fri, 19 Jan 2018) $
!
! $Revision: 1074 $
!
! Description -
!
! Reference -
!
! Calling program - coherens_main, particle_model
!
! External calls - horizontal_particle_drift, particle_concentrations,
!                  particle_releases, vertical_particle_drift
!
! Module calls - error_abort, hrel_coords_curv_rel,
!                hrel_coords_rect_nonunif_rel, hrel_coords_rect_unif_rel,
!                rng_close, rng_normal_arr, rng_opened, vert_coord_index_part
!
!************************************************************************
!
USE grid
USE gridpars  
USE iopars
USE paralpars
USE partpars
USE partswitches
USE partvars
USE physpars
USE switches
USE timepars
USE syspars
USE error_routines, ONLY: error_abort
USE grid_interp, ONLY: hrel_coords_curv_rel, hrel_coords_rect_nonunif_rel, &
                     & hrel_coords_rect_unif_rel
USE particle_routines, ONLY: vert_coord_index_part
USE rng_library, ONLY: rng_close, rng_normal_arr, rng_opened
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, istat, j, npcc, p
REAL :: delxdat, delydat


procname(pglev+1) = 'particle_drift'
CALL log_timer_in(npcc)

!
!1. Initialise
!-------------
!
!---parameters
delxdat = surfacegrids(igrd_model,1)%delxdat
delydat = surfacegrids(igrd_model,1)%delydat

!---new releases
CALL particle_releases

!---activate, new releases, apply land mask
p_110: DO p=1,nopart
   IF (part(p)%drift_state.EQ.2.OR.part(p)%drift_state.EQ.3) THEN
      i = part(p)%icoord; j = part(p)%jcoord
      IF (maskatc(i,j)) THEN
         IF (part(p)%drift_state.EQ.3) part(p)%drift_state = 2
      ELSEIF (part(p)%drift_state.EQ.2) THEN
         part(p)%drift_state = 3
         part(p)%zpos = 0.0
      ENDIF
   ENDIF
ENDDO p_110

!
!2. Random numbers for diffusive drift
!-------------------------------------
!

IF (iopt_part_hdif.GT.0.OR.iopt_part_vdif.GT.0) THEN   
   CALL rng_normal_arr(ranmagpart,nopart,prangenwalk,0.0,1.0,xlo=-1.0,xhi=1.0)
ENDIF
IF (iopt_part_hdif.GT.0) THEN
   CALL rng_normal_arr(randirpart,nopart,prangenwalk,0.0,pi,xlo=-pi,xhi=pi)
ENDIF   

!
!3. Horizontal transport
!-----------------------
!

IF (iopt_part_hdif.GT.0.OR.iopt_part_hadv.GT.0) THEN
   CALL horizontal_particle_drift
ENDIF

!
!4. Update new particle locations
!--------------------------------
!

p_410: DO p=1,nopart
   
!  ---horizontal location cell indices
   IF (part(p)%drift_state.EQ.2) THEN
!     --age      
      part(p)%age = part(p)%age + deltpart/time_convert(ptime_unit)
      SELECT CASE (iopt_grid_htype)
         CASE (1)
            CALL hrel_coords_rect_unif_rel(part(p)%xpos,part(p)%ypos,&
                                  & part(p)%displx,part(p)%disply,delxdat,&
                                  & delydat,(/1,1/),(/nc,nr/),istat,&
                                  & part(p)%icoord,part(p)%jcoord,&
                                  & part(p)%xcoord,part(p)%ycoord,update=.TRUE.)
         CASE (2)
            CALL hrel_coords_rect_nonunif_rel(part(p)%xpos,part(p)%ypos,&
                                  & part(p)%displx,part(p)%disply,&
                                  & gxcoord(1:nc,1),gycoord(1,1:nr),nc,nr,&
                                  & (/1,1/),(/nc,nr/),istat,&
                                  & part(p)%icoord,part(p)%jcoord,&
                                  & xcoord=part(p)%xcoord,&
                                  & ycoord=part(p)%ycoord,update=.TRUE.)
         CASE (3)
            CALL hrel_coords_curv_rel(part(p)%xpos,part(p)%ypos,part(p)%displx,&
                                  & part(p)%disply,gxcoord(1:nc,1:nr),&
                                  & gycoord(1:nc,1:nr),nc,nr,&
                                  & (/1,1/),(/nc,nr/),nsize_sb,istat,2,&
                                  & part(p)%icoord,part(p)%jcoord,&
                                  & xcoord=part(p)%xcoord,&
                                  & ycoord=part(p)%ycoord,update=.TRUE.)
         END SELECT
         
      IF (istat.EQ.1) THEN
         part(p)%drift_state = 4
      ELSEIF (.NOT.maskatc(i,j)) THEN
         part(p)%drift_state = 3
      ENDIF
   ENDIF

ENDDO p_410

CALL error_abort('particle_drift',ierrno_runval)
   
!
!5. Vertical transport
!---------------------
!

IF (iopt_part_vdif.GT.0.OR.iopt_part_vadv.GT.0) THEN
   CAll vertical_particle_drift
ENDIF

!
!6. Activate new releases for next time step
!-------------------------------------------
!

p_610: DO p=1,nopart
   IF (part(p)%ntstart.EQ.nt.AND.part(p)%drift_state.EQ.1) THEN
      part(p)%drift_state = 2
   ENDIF
ENDDO p_610

!
!7. Vertical location/cell index
!-------------------------------
!

p_710: DO p=1,nopart
   IF (part(p)%state.EQ.2) THEN
      IF (part(p)%drift_state.EQ.2) THEN
         CALL vert_coord_index_part(part(p)%icoord,part(p)%jcoord,&
                                  & part(p)%kcoord,part(p)%xcoord,&
                                  & part(p)%ycoord,part(p)%zpos)
      ELSEIF (part(p)%drift_state.GT.2) THEN
         part(p)%zpos = 0.0
      ENDIF
   ENDIF
ENDDO p_710

!
!8. Close random generators
!--------------------------
!

IF (nt.EQ.nstep) THEN
   IF (rng_opened(prangenwalk)) CALL rng_close(prangenwalk)
   IF (rng_opened(prangenlw)) CALL rng_close(prangenlw)
ENDIF

IF (modelid.EQ.modelidpart) THEN
   CALL log_timer_out()
ELSE
   CALL log_timer_out(npcc,itm_part)
ENDIF


RETURN

END SUBROUTINE particle_drift

!========================================================================

SUBROUTINE particle_model
!************************************************************************
!
! *particle_model* COHERENS particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.11
!
! $Date: 2018-01-19 09:43:56 +0100 (Fri, 19 Jan 2018) $
!
! $Revision: 1074 $
!
! Description -
!
! Reference - routine is called if the particle model runs offline or online in
!             parallel mode
!
! Calling program - coherens_main
!
! External calls - combine_particle_phys_data, initialise_particle_model,
!                  meteo_input, particle_concentrations, particle_drift,
!                  particle_trajects, time_series, update_particle_phys_data,
!                  usrdef_part_output, write_partics
!
! Module calls - monitor_files, update_time
!
!************************************************************************
!
USE currents
USE depths  
USE iopars
USE meteo
USE partpars
USE partswitches
USE switches
USE timepars
USE inout_routines, ONLY: monitor_files
USE time_routines, ONLY: log_timer_in, log_timer_out, update_time

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc


procname(pglev+1) = 'particle_model'
CALL log_timer_in(npcc)

!
!1. Initialise particle model
!----------------------------
!

CALL initialise_particle_model  

!---terminate in case of cold start
IF (cold_start) GOTO 1000

!
!2. Time loop
!------------
!

nt_210: DO nt=1,nstep
   
!
!4.1 Restart log file
!--------------------
!

   CALL monitor_files

!
!4.2 Date/time
!-------------
!

   CALL update_time
   IF (.NOT.corrstep) CYCLE nt_210
   
!
!4.3 Store old values
!--------------------
!

   uvel_old = uvel; vvel_old = vvel
   deptotatc_old = deptotatc
   
!
!4.4 Update meteo forcing
!------------------------
!   
   
   IF (iopt_part_wind.EQ.1) CALL meteo_input
      
!
!4.5 Update physical data
!------------------------
!

   IF (iopt_part_model.EQ.2) THEN
      CALL combine_particle_phys_data
   ELSE
      CALL update_particle_phys_data
   ENDIF

!
!4.6 Lagrangian drift
!---------------------
!

   CALL particle_drift

!
!4.7 Eulerian distributions
!--------------------------
!

   CALL particle_concentrations

!
!4.8 Particle output
!-------------------
!   

   IF (iopt_part_model.GT.2.AND.iopt_out_tsers.GT.0) CALL time_series
   CALL particle_trajects
   CALL usrdef_part_output

!
!4.9 Restart file
!----------------
!   

   CALL write_partics
   
ENDDO nt_210
   
1000 CALL log_timer_out(npcc,itm_part)


RETURN

END SUBROUTINE particle_model

!========================================================================

SUBROUTINE vertical_particle_drift
!************************************************************************
!
! *vertical_part_drift* vertical particle drift
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.11
!
! $Date: 2018-01-19 09:43:56 +0100 (Fri, 19 Jan 2018) $
!
! $Revision: 1074 $
!
! Description -
!
! Reference -
!
! Calling program - particle_drift
!
! External calls -
!
! Module calls - Cint_at_p, Wint_at_p
!
!************************************************************************
!
USE currents
USE density
USE depths
use grid
use gridpars
USE iopars
USE partpars
USE partswitches
USE partpars
USE partvars
USE physpars
USE syspars
USE timepars
USE particle_routines, ONLY: Cint_at_p, Wint_at_p
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j, l, p
REAL :: bfac1, bfac2, bfac3, buo, dcrit, densp, depmeanp, deptotp, deptotp_old,&
      & diam, displ_adv, displ_dif, gsig, gsigold, vfac, x, y, vdifcoefp, &
      & wvel_adv, wvel_buo, wvel_dif
REAL (KIND=kndrlong) :: zpos


procname(pglev+1) = 'vertical_particle_drift'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!

IF (iopt_part_vdif.EQ.1) THEN
   vfac = SQRT(2.0*zdifpart_cst/delt3d)
ELSEIF (iopt_part_vdif.EQ.2) THEN
   vfac = SQRT(2.0/delt3d)
ENDIF

IF (iopt_part_dens.GT.0) THEN
   bfac1 = gacc_ref/(18.0*kinvisc_cst)
   bfac2 = SQRT(8.0*gacc_ref/3.0)
   bfac3 = 9.52*(kinvisc_cst*kinvisc_cst/gacc_ref)**(1.0/3.0)
ENDIF
   
!
!2. Vertical transport
!---------------------
!   
   
   p_200: DO p=1,nopart

      IF (part(p)%state.EQ.2.AND.part(p)%drift_state.EQ.2) THEN

         i = part(p)%icoord; j = part(p)%jcoord
         x = part(p)%xcoord; y = part(p)%ycoord
         zpos = part(p)%zpos
         wvel_adv = 0.0; wvel_buo= 0.0
         buo = 0.0

!
!2.1 Vertical diffusive current   
!------------------------------
!
      
         IF (iopt_part_vdif.EQ.1) THEN
            wvel_dif = ranmagpart(p)*vfac
         ELSEIF (iopt_part_vdif.EQ.2) THEN
            vdifcoefp = Wint_at_p(vdifcoefpart(i-1:i+1,j-1:j+1,:),i,j,x,y,zpos)
            wvel_dif = vfac*ranmagpart(p)*SQRT(vdifcoefp)
         ELSE
            wvel_dif = 0.0
         ENDIF

!
!2.2 Vertical advective current
!------------------------------
!

         IF (iopt_part_vadv.EQ.1) THEN
            wvel_adv = Wint_at_p(wvel(i-1:i+1,j-1:j+1,:),i,j,x,y,zpos)
         ENDIF

!
!2.3 Vertical buoyancy current
!-----------------------------
!

         IF (iopt_part_dens.GT.0) THEN
            l = part(p)%label
            IF (iopt_part_dens.EQ.1) THEN
               buo = 1.0-rhoconc(l)/density_ref
            ELSE
               densp = Cint_at_p(dens(i-1:i+1,j-1:j+1,:),i,j,x,y,zpos)
               buo = 1.0-rhoconc(l)/densp            
            ENDIF
            diam = diamconc(l)
            dcrit = bfac3/buo**(1.0/3.0)
            IF (buo.LT.0.0.OR.diam.LE.dcrit) THEN
               wvel_buo = bfac1*diam**2*buo
            ELSE
               wvel_buo = bfac2*SQRT(diam*buo)
            ENDIF
         ENDIF

!
!2.4 Update vertical positions
!-----------------------------
!      
!        ---initialise parameters
         depmeanp = Cint_at_p(depmeanatc(i-1:i+1,j-1:j+1),i,j,x,y)
         deptotp = Cint_at_p(deptotatc(i-1:i+1,j-1:j+1),i,j,x,y)
         deptotp_old = Cint_at_p(deptotatc_old(i-1:i+1,j-1:j+1),i,j,x,y)
         gsigold = (zpos+depmeanp)/deptotp_old
         displ_adv = deltpart*(wvel_adv+wvel_buo)/deptotp
         displ_adv = MIN(1.0-gsigold,MAX(displ_adv,-gsigold))
         displ_dif = deltpart*wvel_dif/deptotp
!        ---add displacements         
         gsig = gsigold + displ_adv + displ_dif
!        ---diffusive reflection at the surface
         IF (displ_dif.GT.0.0.AND.gsig.GT.1.0) THEN
            gsig = 2.0 - gsig
         ENDIF
!        ---diffusive reflection at the bottom
         IF (displ_dif.LT.0.0.AND.gsig.LT.0.0) THEN
            gsig = -gsig
         ENDIF
         gsig = MIN(1.0,MAX(gsig,0.0))
         part(p)%zpos = gsig*deptotp - depmeanp
!        ---update vertical location
         IF (gsig.EQ.0.0) THEN
            part(p)%drift_state = 5
            part(p)%zpos = -depmeanp
         ENDIF
         part(p)%displz = part(p)%zpos - zpos
      ENDIF
      
ENDDO p_200

CALL log_timer_out()


RETURN

END SUBROUTINE vertical_particle_drift
