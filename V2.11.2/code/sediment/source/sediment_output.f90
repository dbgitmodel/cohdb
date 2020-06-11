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

MODULE sediment_output
!************************************************************************
!
! *sediment_output* Routines defining sediment output data
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)sediment_output.f90  V2.11.1
!
! $Date: 2018-05-23 10:09:17 +0200 (Wed, 23 May 2018) $
!
! $Revision: 1138 $
!
! Description -
!
! Reference -
!
! Routines - define_sed0d_vals, define_sed2d_vals, define_sed3d_vals
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS

!========================================================================

SUBROUTINE define_sed0d_vals(outdat,ivarid,f,oopt)
!************************************************************************
!
! *define_sed0d_vals* Define 0-D sediment output data
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)sediment_output.f90  V2.11.1
!
! Description -
!
! Calling program - define_out0d_vals
!
! Module calls - error_alloc, inquire_var, max_vars, min_vars, max_vars,
!                vector_mag_arr, Uarr_at_C, Varr_at_C, Warr_at_C   
!  
!************************************************************************
!
USE currents  
USE dararrays
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphpars
USE physpars
USE sedarrays
USE sedids
USE sedpars
USE switches
USE syspars
USE array_interp, ONLY: Uarr_at_C, Varr_at_C, Warr_at_C
USE error_routines, ONLY: error_alloc
USE math_library, ONLY: vector_mag_arr_atc
USE modvars_routines, ONLY: inquire_var
USE paral_utilities, ONLY: max_vars, min_vars, sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
INTEGER, INTENT(IN) :: f, ivarid, oopt
REAL, INTENT(OUT) :: outdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    0-D output data value
!*ivarid*  INTEGER Variable key id
!*f*       INTEGER Fraction number
!*oopt*    INTEGER Type of operator
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ff, k, nodim, nzdim
INTEGER, DIMENSION(4) :: nhdims = 0
INTEGER, DIMENSION(5) :: ngdims
REAL :: atot, vtot
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: area, out2d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: out3d, volume


procname(pglev+1) = 'define_sed0d_vals'
CALL log_timer_in()

!
!1. Without operator
!-------------------
!

IF (oopt.EQ.oopt_null) THEN

   SELECT CASE (ivarid)
      CASE (iarr_sediment_balance); outdat = sediment_balance
      CASE (iarr_sediment_obc); outdat = sediment_obc
      CASE (iarr_sediment_vol); outdat = sediment_vol
   END SELECT

!
!2. With operator
!----------------
!

ELSE

   CALL inquire_var(ivarid,nrank=nodim,global_dims=ngdims)

!
!2.1 Allocate arrays
!-------------------
!

   IF (nodim.EQ.2) THEN
      ALLOCATE (out2d(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('out2d',2,(/ncloc,nrloc/),kndrtype)
      out2d = 0.0
   ELSEIF (nodim.EQ.3) THEN
      nzdim = ngdims(3)
      ALLOCATE (out3d(ncloc,nrloc,nzdim),STAT=errstat)
      CALL error_alloc('out3d',3,(/ncloc,nrloc,nzdim/),kndrtype)
      out3d = 0.0
   ENDIF

!
!2.2 Area and volume of the global domain without land
!-----------------------------------------------------
!

   IF (oopt.EQ.oopt_mean.OR.oopt.EQ.oopt_int) THEN
      ALLOCATE (area(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('area',2,(/ncloc,nrloc/),kndrtype)
      area = delxatc(1:ncloc,1:nrloc)*delyatc(1:ncloc,1:nrloc)
      IF (oopt.EQ.oopt_mean) THEN 
         CALL sum2_vars(area,atot,nhdims,'C  ',0,commall=.TRUE.)
      ENDIF
      IF (nodim.EQ.3) THEN
         ALLOCATE (volume(ncloc,nrloc,nzdim),STAT=errstat)
         CALL error_alloc('volume',3,(/ncloc,nrloc,nzdim/),kndrtype)
!         IF (nzdim.EQ.nb) THEN
!            k_221: DO k=1,nb
!               volume(:,:,k) = area*bed_layer_thickness(:,:,k)
!            ENDDO k_221
!         ELSE
            k_310: DO k=1,nzdim
               volume(:,:,k) = area*delzatc(1:ncloc,1:nrloc,k)
            ENDDO k_310
!         ENDIF
         IF (oopt.EQ.oopt_mean) THEN
            CALL sum2_vars(volume,vtot,nhdims,'C  ',0,commall=.TRUE.)
         ENDIF
      ENDIF
   ENDIF

!
!2.3 Store data
!---------------
!

   SELECT CASE (ivarid)

!     ---2-D variable
      CASE (iarr_active_layer_thickness); out2d = active_layer_thickness
      CASE (iarr_bdragcoefatc_sed); out2d = bdragcoefatc_sed
      CASE (iarr_bed_update); out2d = bed_update
      CASE (iarr_bed_update_dep); out2d = bed_update_dep(:,:,f)
      CASE (iarr_bed_update_ero); out2d = bed_update_ero(:,:,f)
      CASE (iarr_bed_update_int); out2d = bed_update_int
      CASE (iarr_beta_sed); out2d = beta_sed(:,:,f)
      CASE (iarr_bfricvel_sed); out2d = SQRT(bstresatc_sed(1:ncloc,1:nrloc))
      CASE (iarr_bfricvel_mean_sed); out2d = SQRT(bstresatc_mean_sed)
      CASE (iarr_bfricvel_wav_sed)
         out2d = SQRT(bstresatc_wav_sed(1:ncloc,1:nrloc))
      CASE (iarr_bottom_sed_dep); out2d = bottom_sed_dep(:,:,f)
      CASE (iarr_bottom_sed_ero); out2d = bottom_sed_ero(:,:,f)
      CASE (iarr_bottom_sed_flux); out2d = bottom_sed_flux(:,:,f)
      CASE (iarr_bstres_cr); out2d = density_ref*bstres_cr(1:ncloc,1:nrloc,f)
      CASE (iarr_bstresatc_sed)
         out2d = density_ref*bstresatc_sed(1:ncloc,1:nrloc)
      CASE (iarr_bstresatc_mean_sed)
         out2d = density_ref*bstresatc_mean_sed
      CASE (iarr_bstresatc_wav_sed)
         out2d = density_ref*bstresatc_wav_sed(1:ncloc,1:nrloc)
      CASE (iarr_ceq); out2d = rhos(f)*ceq(:,:,f)
      CASE (iarr_cref); out2d = rhos(f)*cref(1:ncloc,1:nrloc,f)
      CASE (iarr_dar_morph); out2d = dar_morph(:,:,f)
      CASE (iarr_d50_bed); out2d = d50_bed(1:ncloc,1:nrloc)
      CASE (iarr_fwave_sed); out2d = fwave_sed(1:ncloc,1:nrloc)
      CASE (iarr_height_c); out2d = height_c(:,:,f)
      CASE (iarr_qbedatc)
         CALL vector_mag_arr_atc(qbedatu(1:ncloc+1,:,f),qbedatv(:,1:nrloc+1,f),&
                               & 1,1,1,1,ivarid,.FALSE.,vecmag=out2d)
         out2d = rhos(f)*out2d
      CASE (iarr_qbedatu)
         CALL Uarr_at_C(qbedatu(1:ncloc+1,:,f),out2d,1,1,(/1,1,1/),&
                     & (/ncloc+1,nrloc,1/),1,ivarid,.FALSE.)
         out2d = rhos(f)*out2d
      CASE (iarr_qbedatv)
         CALL Varr_at_C(qbedatv(:,1:nrloc+1,f),out2d,1,1,(/1,1,1/),&
                     & (/ncloc,nrloc+1,1/),1,ivarid,.FALSE.)
         out2d = rhos(f)*out2d
      CASE (iarr_qtotatc)
         CALL vector_mag_arr_atc(qtotatu(1:ncloc+1,:,f),qtotatv(:,1:nrloc+1,f),&
                               & 1,1,1,1,ivarid,.FALSE.,vecmag=out2d)
         out2d = rhos(f)*out2d
      CASE (iarr_qtotatu)
         CALL Uarr_at_C(qtotatu(1:ncloc+1,:,f),out2d,1,1,(/1,1,1/),&
                     & (/ncloc+1,nrloc,1/),1,ivarid,.FALSE.)
         out2d = rhos(f)*out2d
      CASE (iarr_qtotatv)
         CALL Varr_at_C(qtotatv(:,1:nrloc+1,f),out2d,1,1,(/1,1,1/),&
                     & (/ncloc,nrloc+1,1/),1,ivarid,.FALSE.)
         out2d = rhos(f)*out2d
      CASE (iarr_sed_avail); out2d = sed_avail(:,:,f)
      CASE (iarr_sed_avail_tot); out2d = sed_avail_tot
      CASE (iarr_t_equil); out2d = t_equil(:,:,f)
      CASE (iarr_ubstresatc_sed)
         CALL Uarr_at_C(uvel(1:ncloc+1,1:nrloc,1),out2d,1,1,(/1,1,1/),&
                     & (/ncloc+1,nrloc,1/),1,ivarid,.FALSE.)
         out2d = density_ref*SQRT(bstresatc_sed(1:ncloc,1:nrloc))*out2d
      CASE (iarr_vbstresatc_sed)
         CALL Varr_at_C(vbstresatv_sed,out2d,1,1,(/1,1,1/),&
                     & (/ncloc,nrloc+1,1/),1,ivarid,.FALSE.)
         out2d = density_ref*out2d
      CASE (iarr_wavethickatc_sed); out2d = wavethickatc_sed
      CASE (iarr_zroughatc_sed); out2d = zroughatc_sed(1:ncloc,1:nrloc)
         
!     ---3-D variable
      CASE (iarr_bed_fraction); out3d = bed_fraction(1:ncloc,1:nrloc,:,f)
      CASE (iarr_bed_layer_thickness); out3d = bed_layer_thickness
      CASE (iarr_bed_porosity); out3d = bed_porosity(:,:,:)
      CASE (iarr_beta_state_sed); out3d = beta_state_sed(1:ncloc,1:ncloc,:,f)
      CASE (iarr_cnump); out3d = cnump(1:ncloc,1:nrloc,:,f)
      CASE (iarr_ctot)
         out3d = 0.0
         IF (iopt_sed.EQ.1) THEN
            k_231: DO k=1,nz
            ff_231: DO ff=1,nf
               WHERE (maskatc_int)
                  out3d(:,:,k) = out3d(:,:,k) + &
                               & rhos(ff)*cvol(1:ncloc,1:nrloc,k,ff)
               END WHERE
            ENDDO ff_231
            ENDDO k_231
         ELSE
            k_232: DO k=1,nz
               WHERE (maskatc_int)
                  out3d(:,:,k) = rhos(1)*cvol(1:ncloc,1:nrloc,k,1) + &
                        & floc_dens(1:ncloc,1:nrloc,k)*cvol(1:ncloc,1:nrloc,k,2)
               END WHERE
            ENDDO k_232
         ENDIF
      CASE (iarr_cvol); out3d = rhos(f)*cvol(1:ncloc,1:ncloc,:,f)
      CASE (iarr_dar_sediment); out3d = dar_sediment(:,:,:,f)
      CASE (iarr_densw); out3d = densw(1:ncloc,1:nrloc,:)
      CASE (iarr_floc_dens); out3d = floc_dens(1:ncloc,1:nrloc,:)
      CASE (iarr_floc_dia); out3d = 1.0E+06*floc_dia
      CASE (iarr_floc_F)
         out3d = floc_dens(1:ncloc,1:nrloc,:)*cvol(1:ncloc,1:nrloc,:,2)
      CASE (iarr_floc_nc); out3d = floc_nc(1:ncloc,1:nrloc,:)
      CASE (iarr_floc_P)
         out3d = rhos(1)*cvol(1:ncloc,1:nrloc,:,1)
      CASE (iarr_floc_T)
         out3d = rhos(1)*cvol(1:ncloc,1:nrloc,:,3)
      CASE (iarr_floc_vol); out3d = floc_vol
      CASE (iarr_qsusatc)
         CALL vector_mag_arr_atc(qsusatu(:,:,:,f),qsusatv(:,:,:,f),1,1,nz,1,&
                               & ivarid,.FALSE.,vecmag=out3d)
         out3d = rhos(f)*out3d
      CASE (iarr_qsusatu)
         CALL Uarr_at_C(qsusatu(1:ncloc+1,:,:,f),out3d,1,1,(/1,1,1/),&
                     & (/ncloc+1,nrloc,nz/),1,ivarid,.FALSE.)
         out3d = rhos(f)*out3d
      CASE (iarr_qsusatv)
         CALL Varr_at_C(qsusatv(:,1:nrloc+1,:,f),out3d,1,1,(/1,1,1/),&
                     & (/ncloc,nrloc+1,nz/),1,ivarid,.FALSE.)
         out3d = rhos(f)*out3d
      CASE (iarr_vdiffcoef_sed)
         CALL Warr_at_C(vdiffcoef_sed(:,:,:,f),out3d,(/1,1,1/),&
                     & (/ncloc,nrloc,nz+1/),1,ivarid,.FALSE.)
      CASE (iarr_wfall)
         CALL Warr_at_C(wfall(1:ncloc,1:nrloc,:,f),out3d,(/1,1,1/),&
                     & (/ncloc,nrloc,nz+1/),1,ivarid,.FALSE.)

   END SELECT

!
!2.4 Apply operator
!------------------
!

   IF (nodim.EQ.2) THEN

      SELECT CASE (oopt)

         CASE (oopt_min)
            CALL min_vars(out2d,outdat,ivarid,commall=.TRUE.,mask=maskatc_int)
         CASE (oopt_max)
            CALL max_vars(out2d,outdat,ivarid,commall=.TRUE.,mask=maskatc_int)
         CASE (oopt_mean)
            CALL sum2_vars(area*out2d,outdat,nhdims,'C  ',ivarid,commall=.TRUE.)
            outdat = outdat/atot
         CASE (oopt_int)
            CALL sum2_vars(area*out2d,outdat,nhdims,'C  ',ivarid,commall=.TRUE.)

      END SELECT

   ELSEIF (nodim.EQ.3) THEN

      SELECT CASE (oopt)

         CASE (oopt_min)
            CALL min_vars(out3d,outdat,ivarid,commall=.TRUE.,mask=maskatc_int)
         CASE (oopt_max)
            CALL max_vars(out3d,outdat,ivarid,commall=.TRUE.,mask=maskatc_int)
         CASE (oopt_mean)
            CALL sum2_vars(volume*out3d,outdat,nhdims,'C  ',ivarid,&
                         & commall=.TRUE.)
            outdat = outdat/vtot

         CASE (oopt_int)
            CALL sum2_vars(volume*out3d,outdat,nhdims,'C  ',ivarid,&
                         & commall=.TRUE.)

      END SELECT

   ENDIF

!
!2.5 Deallocate arrays
!---------------------
!

   IF (ALLOCATED(out2d)) DEALLOCATE (out2d)
   IF (ALLOCATED(out3d)) DEALLOCATE (out3d)
   IF (ALLOCATED(area)) DEALLOCATE (area)
   IF (ALLOCATED(volume)) DEALLOCATE (volume)

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE define_sed0d_vals

!========================================================================

SUBROUTINE define_sed2d_vals(outdat,i,j,ivarid,f,oopt,kc,dep)
!************************************************************************
!
! *define_sed2d_vals*  Define sediment output data at a 2-D grid location
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)sediment_output.f90  V2.11.1
!
! Description -
!
! Calling program - define_out2d_vals
!
! Module calls - intpol1d_model_to_dep, inquire_var, Uarr_at_C, Uvar_at_C,
!                Varr_at_C, vector_mag_var_atc, Vvar_at_C, Warr_at_C, Wvar_at_C
!
!************************************************************************
!
USE currents  
USE dararrays
USE depths
USE grid
USE gridpars
USE morpharrays
USE physpars
USE sedarrays
USE sedids
USE sedpars
USE switches
USE array_interp, ONLY: Uarr_at_C, Uvar_at_C, Varr_at_C, Vvar_at_C, Warr_at_C, &
                      & Wvar_at_C
USE grid_interp, ONLY: intpol1d_model_to_dep
USE math_library, ONLY: vector_mag_var_atc
USE modvars_routines, ONLY: inquire_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: f, i, ivarid, j, kc, oopt
REAL, INTENT(IN) :: dep
REAL, INTENT(OUT) :: outdat
!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    2-D output data value
!*i*       INTEGER (Local) X-index of the data location
!*j*       INTEGER (Local) Y-index of the data location
!*ivarid*  INTEGER Variable key id
!*f*       INTEGER Fraction number
!*oopt*    INTEGER Type of operator
!*kc*      INTEGER Vertical level at which a 3-D variable is evaluated
!*dep*     REAL    Depth below the surface at which a 3-D variable is evaluated
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: bedflag
INTEGER :: ff, k, nzdim
INTEGER, DIMENSION(5) :: ngdims
REAL, DIMENSION(MAX(nz,nb)) :: out1d
REAL, DIMENSION(1,1,MAX(nz,nb)) :: out3d


!
!1. Without operator
!-------------------
!

IF (oopt.EQ.oopt_null) THEN

   SELECT CASE (ivarid)

      CASE (iarr_active_layer_thickness); outdat = active_layer_thickness(i,j)
      CASE (iarr_bdragcoefatc_sed); outdat = bdragcoefatc_sed(i,j)
      CASE (iarr_bed_update); outdat = bed_update(i,j)
      CASE (iarr_bed_update_dep); outdat = bed_update_dep(i,j,f)
      CASE (iarr_bed_update_ero); outdat = bed_update_ero(i,j,f)
      CASE (iarr_bed_update_int); outdat = bed_update_int(i,j)
      CASE (iarr_beta_sed); outdat = beta_sed(i,j,f)
      CASE (iarr_bfricvel_sed); outdat = SQRT(bstresatc_sed(i,j))
      CASE (iarr_bfricvel_mean_sed); outdat = SQRT(bstresatc_mean_sed(i,j))
      CASE (iarr_bfricvel_wav_sed); outdat = SQRT(bstresatc_wav_sed(i,j)) 
      CASE (iarr_bottom_sed_dep); outdat = bottom_sed_dep(i,j,f)
      CASE (iarr_bottom_sed_ero); outdat = bottom_sed_ero(i,j,f)
      CASE (iarr_bottom_sed_flux); outdat = bottom_sed_flux(i,j,f)
      CASE (iarr_bstres_cr); outdat = density_ref*bstres_cr(i,j,f)
      CASE (iarr_bstresatc_sed); outdat = density_ref*bstresatc_sed(i,j)
      CASE (iarr_bstresatc_mean_sed)
         outdat = density_ref*bstresatc_mean_sed(i,j)
      CASE (iarr_bstresatc_wav_sed); outdat = density_ref*bstresatc_wav_sed(i,j)
      CASE (iarr_ceq); outdat = rhos(f)*ceq(i,j,f)
      CASE (iarr_cref); outdat = rhos(f)*cref(i,j,f)
      CASE (iarr_dar_morph); outdat = dar_morph(i,j,f)
      CASE (iarr_d50_bed); outdat = d50_bed(i,j)
      CASE (iarr_fwave_sed); outdat = fwave_sed(i,j)
      CASE (iarr_height_c); outdat = height_c(i,j,f)
      CASE (iarr_qbedatc)
         CALL vector_mag_var_atc(qbedatu(i:i+1,j,f),qbedatv(i,j:j+1,f),i,j,&
                               & 1,1,1,vecmag=outdat)
         outdat = rhos(f)*outdat
      CASE (iarr_qbedatu)
         outdat = rhos(f)*Uvar_at_C(qbedatu(i:i+1,j,f),i,j,nz,1,1)
      CASE (iarr_qbedatv)
         outdat = rhos(f)*Vvar_at_C(qbedatv(i,j:j+1,f),i,j,nz,1,1)
      CASE (iarr_qtotatc)
         CALL vector_mag_var_atc(qtotatu(i:i+1,j,f),qtotatv(i,j:j+1,f),i,j,&
                               & 1,1,1,vecmag=outdat)
         outdat = rhos(f)*outdat
      CASE (iarr_qtotatu)
         outdat = rhos(f)*Uvar_at_C(qtotatu(i:i+1,j,f),i,j,1,1,1)
      CASE (iarr_qtotatv)
         outdat = rhos(f)*Vvar_at_C(qtotatv(i,j:j+1,f),i,j,1,1,1)
      CASE (iarr_sed_avail); outdat = sed_avail(i,j,f)
      CASE (iarr_sed_avail_tot); outdat = sed_avail_tot(i,j)
      CASE (iarr_t_equil); outdat = t_equil(i,j,f)
      CASE (iarr_ubstresatc_sed)
         outdat = density_ref*Uvar_at_C(ubstresatu_sed(i:i+1,j),i,j,nz,1,1)
      CASE (iarr_vbstresatc_sed)
         outdat = density_ref*Vvar_at_C(vbstresatv_sed(i,j:j+1),i,j,nz,1,1)
      CASE (iarr_wavethickatc_sed); outdat = wavethickatc_sed(i,j)
      CASE (iarr_zroughatc_sed); outdat = zroughatc_sed(i,j)
   END SELECT

!
!2. 3-D variable at a given level
!--------------------------------
!

ELSEIF (oopt.EQ.oopt_klev) THEN

   SELECT CASE (ivarid)
      CASE (iarr_bed_fraction); outdat = bed_fraction(i,j,kc,f)
      CASE (iarr_bed_layer_thickness); outdat = bed_layer_thickness(i,j,kc)
      CASE (iarr_bed_porosity); outdat = bed_porosity(i,j,kc)
      CASE (iarr_beta_state_sed); outdat = beta_state_sed(i,j,kc,f)
      CASE (iarr_cnump); outdat = cnump(i,j,kc,f)
      CASE (iarr_ctot)
         IF (iopt_sed.EQ.1) THEN
            outdat = SUM(rhos*cvol(i,j,kc,:))
         ELSE
            outdat = rhos(1)*cvol(i,j,kc,1) + floc_dens(i,j,kc)*cvol(i,j,kc,2)
         ENDIF
      CASE (iarr_cvol); outdat = rhos(f)*cvol(i,j,kc,f)
      CASE (iarr_dar_sediment); outdat = dar_sediment(i,j,kc,f)
      CASE (iarr_densw); outdat = densw(i,j,kc)
      CASE (iarr_floc_dens); outdat = floc_dens(i,j,kc)
      CASE (iarr_floc_dia); outdat = 1.0E+06*floc_dia(i,j,kc)
      CASE (iarr_floc_F); outdat = floc_dens(i,j,kc)*cvol(i,j,kc,2)
      CASE (iarr_floc_nc); outdat = floc_nc(i,j,kc)
      CASE (iarr_floc_P); outdat = rhos(1)*cvol(i,j,kc,1)
      CASE (iarr_floc_T); outdat = rhos(1)*cvol(i,j,kc,3)
      CASE (iarr_floc_vol); outdat = floc_vol(i,j,kc)
      CASE (iarr_qsusatc)
         CALL vector_mag_var_atc(qsusatu(i:i+1,j,kc,f),qsusatv(i,j:j+1,kc,f),&
                               & i,j,kc,1,1,vecmag=outdat)
         outdat = rhos(f)*outdat
      CASE (iarr_qsusatu)
         outdat = rhos(f)*Uvar_at_C(qsusatu(i:i+1,j,kc,f),i,j,kc,1,1)
      CASE (iarr_qsusatv)
         outdat = rhos(f)*Vvar_at_C(qsusatv(i,j:j+1,kc,f),i,j,kc,1,1)
      CASE (iarr_vdiffcoef_sed)
         outdat = Wvar_at_C(vdiffcoef_sed(i,j,kc:kc+1,f),i,j)
      CASE (iarr_wfall); outdat = Wvar_at_C(wfall(i,j,kc:kc+1,f),i,j)

   END SELECT

!
!3. Other operators
!------------------
!

ELSE

!
!3.1 Store data
!--------------
!

   bedflag = .FALSE.; nzdim = nz

   SELECT CASE (ivarid)

      CASE (iarr_bed_fraction)
         nzdim = nb; bedflag = .TRUE.
         out1d(1:nb) = bed_fraction(i,j,:,f)
      CASE (iarr_bed_layer_thickness)
         nzdim = nb; bedflag = .TRUE.
         out1d(1:nb) = bed_layer_thickness(i,j,:)
      CASE (iarr_bed_porosity)
         nzdim = nb; bedflag = .TRUE.
         out1d(1:nb) = bed_porosity(i,j,:)
      CASE (iarr_beta_state_sed); out1d(1:nz) = beta_state_sed(i,j,:,f)
      CASE (iarr_cnump); out1d(1:nz) = cnump(i,j,:,f)
      CASE (iarr_ctot)
         out1d = 0.0
         IF (iopt_sed.EQ.1) THEN
            k_311: DO k=1,nz
            ff_311: DO ff=1,nf
               out1d(k) = out1d(k) + rhos(ff)*cvol(i,j,k,ff)
            ENDDO ff_311
            ENDDO k_311
         ELSE
            k_312: do k=1,nz
               out1d(k) = out1d(k) + rhos(1)*cvol(i,j,k,1) + &
                        & floc_dens(i,j,k)*cvol(i,j,k,2)
            ENDDO k_312
         ENDIF
      CASE (iarr_cvol)
         out1d(1:nz) = rhos(f)*cvol(i,j,:,f)
      CASE (iarr_dar_sediment); out1d(1:nz) = dar_sediment(i,j,:,f)
      CASE (iarr_densw); out1d(1:nz) = densw(i,j,:)
      CASE (iarr_floc_dens); out1d(1:nz) = floc_dens(i,j,:)
      CASE (iarr_floc_dia); out1d(1:nz) = 1.0E+06*floc_dia(i,j,:)
      CASE (iarr_floc_F); out1d(1:nz) = floc_dens(i,j,:)*cvol(i,j,:,2)
      CASE (iarr_floc_nc); out1d(1:nz) = floc_nc(i,j,:)
      CASE (iarr_floc_P); out1d(1:nz) = rhos(1)*cvol(i,j,:,1)
      CASE (iarr_floc_T); out1d(1:nz) = rhos(1)*cvol(i,j,:,3)
      CASE (iarr_floc_vol); out1d(1:nz) = floc_vol(i,j,:)
      CASE (iarr_qsusatc)
         k_313: DO k=1,nz
            CALL vector_mag_var_atc(qsusatu(i:i+1,j,k,f),qsusatv(i,j:j+1,k,f),&
                                   & i,j,k,1,1,vecmag=out1d(k))
         ENDDO k_313
         out1d(1:nz) = rhos(f)*out1d(1:nz)
      CASE (iarr_qsusatu)
         CALL Uarr_at_C(qsusatu(i:i+1,j:j,:,f),out3d(:,:,1:nz),1,1,(/i,j,1/),&
                  & (/i+1,j,nz/),1,ivarid,.FALSE.)
         out1d(1:nz) = rhos(f)*out3d(1,1,:)
      CASE (iarr_qsusatv)
         CALL Uarr_at_C(qsusatv(i:i,j:j+1,:,f),out3d(:,:,1:nz),1,1,(/i,j,1/),&
                     & (/i,j+1,nz/),1,ivarid,.FALSE.)
         out1d(1:nz) = rhos(f)*out3d(1,1,:)
      CASE (iarr_vdiffcoef_sed)
         CALL Warr_at_C(vdiffcoef_sed(i:i,j:j,:,f),out3d(:,:,1:nz),(/i,j,1/),&
                     & (/i,j,nz+1/),1,ivarid,.FALSE.)
         out1d(1:nz) = out3d(1,1,1:nz)
      CASE (iarr_wfall)
         CALL Warr_at_C(wfall(i:i,j:j,:,f),out3d(:,:,1:nz),(/i,j,1/),&
                     & (/i,j,nz+1/),1,ivarid,.FALSE.)
         out1d(1:nz) = out3d(1,1,1:nz)
   END SELECT

!
!3.2 Apply operator
!------------------
!

   SELECT CASE (oopt)

      CASE (oopt_min); outdat = MINVAL(out1d)
      CASE (oopt_max); outdat = MAXVAL(out1d)
      CASE (oopt_mean)
         CALL inquire_var(ivarid,global_dims=ngdims)
         nzdim = ngdims(3)
         IF (nzdim.EQ.nb) THEN
            outdat = SUM(bed_layer_thickness(i,j,:)*out1d(1:nb))&
                      & /sed_avail_tot(i,j)
         ELSE
            outdat = SUM(delzatc(i,j,:)*out1d(1:nz))/deptotatc(i,j)
         ENDIF
      CASE (oopt_int)
        CALL inquire_var(ivarid,global_dims=ngdims)
        IF (bedflag) THEN
           outdat = SUM(bed_layer_thickness(i,j,:)*out1d(1:nb))
        ELSE
           outdat = SUM(delzatc(i,j,:)*out1d(1:nz))
        ENDIF
      CASE (oopt_dep)
         CALL inquire_var(ivarid,global_dims=ngdims)
         IF (.NOT.bedflag) THEN
            CALL intpol1d_model_to_dep(out1d,outdat,i,j,dep,'C  ')
         ENDIF

   END SELECT

ENDIF


RETURN

END SUBROUTINE define_sed2d_vals

!========================================================================

SUBROUTINE define_sed3d_vals(outdat,i,j,k,ivarid,f)
!************************************************************************
!
! *define_sed3d_vals*  Define sediment output data at a 3-D grid location
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)sediment_output.f90  V2.11.1
!
! Description -
!
! Calling program - define_out2d_vals, define_out3d_vals
!
! Module calls - Uvar_at_C, vector_mag_var_atc, Vvar_at_C, Wvar_at_C
!
! Reference -
!
!************************************************************************
!
USE morpharrays
USE sedarrays
USE sedids
USE switches
USE syspars
USE array_interp, ONLY: Uvar_at_C, Vvar_at_C, Wvar_at_C
USE math_library, ONLY: vector_mag_var_atc

!
!*Arguments
!
INTEGER, INTENT(IN) :: f, i, ivarid, j, k
REAL, INTENT(OUT) :: outdat
!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    3-D output data value
!*i*       INTEGER (Local) X-index of the data location
!*j*       INTEGER (Local) Y-index of the data location
!*k*       INTEGER Vertical index of the data location
!*ivarid*  INTEGER Variable key id
!*f*       INTEGER Fraction number
!
!------------------------------------------------------------------------------
!


SELECT CASE (ivarid)

   CASE (iarr_bed_fraction); outdat = bed_fraction(i,j,k,f)
   CASE (iarr_bed_layer_thickness); outdat = bed_layer_thickness(i,j,k)
   CASE (iarr_bed_porosity); outdat = bed_porosity(i,j,k)
   CASE (iarr_beta_state_sed); outdat = beta_state_sed(i,j,k,f)
   CASE (iarr_cnump); outdat = cnump(i,j,k,f)
   CASE (iarr_ctot)
      IF (iopt_sed.EQ.1.OR.f.NE.2) THEN
         outdat = SUM(rhos*cvol(i,j,k,:))
      ELSE
         outdat = rhos(1)*cvol(i,j,k,1)+floc_dens(i,j,k)*cvol(i,j,k,2)
      ENDIF
   CASE (iarr_cvol); outdat = rhos(f)*cvol(i,j,k,f)
   CASE (iarr_dar_sediment); outdat = dar_sediment(i,j,k,f)
   CASE (iarr_densw); outdat = densw(i,j,k)
   CASE (iarr_floc_dens); outdat = floc_dens(i,j,k)
   CASE (iarr_floc_dia); outdat = 1.0E+06*floc_dia(i,j,k)
   CASE (iarr_floc_F); outdat = floc_dens(i,j,k)*cvol(i,j,k,2)
   CASE (iarr_floc_nc); outdat = floc_nc(i,j,k)
   CASE (iarr_floc_P); outdat = rhos(1)*cvol(i,j,k,1)
   CASE (iarr_floc_T); outdat = rhos(1)*cvol(i,j,k,3)
   CASE (iarr_floc_vol); outdat = floc_vol(i,j,k)
   CASE (iarr_qsusatc)
      CALL vector_mag_var_atc(qsusatu(i:i+1,j,k,f),qsusatv(i,j:j+1,k,f),&
                            & i,j,k,1,1,vecmag=outdat)
      outdat = rhos(f)*outdat
   CASE (iarr_qsusatu)
      outdat = rhos(f)*Uvar_at_C(qsusatu(i:i+1,j,k,f),i,j,k,1,1)
   CASE (iarr_qsusatv)
      outdat = rhos(f)*Vvar_at_C(qsusatv(i,j:j+1,k,f),i,j,k,1,1)
   CASE (iarr_vdiffcoef_sed); outdat = Wvar_at_C(vdiffcoef_sed(i,j,k:k+1,f),i,j)
   CASE (iarr_wfall); outdat = Wvar_at_C(wfall(i,j,k:k+1,f),i,j)

END SELECT


RETURN

END SUBROUTINE define_sed3d_vals

END MODULE sediment_output
