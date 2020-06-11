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

MODULE model_output
!************************************************************************
!
! *model_output* Routines defining user output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)model_output.f90  V2.11.2
!
! $Date: 2018-09-28 13:01:37 +0200 (Fri, 28 Sep 2018) $
!
! $Revision: 1188 $
!
! Description -
!
! Routines - define_out0d_vals, define_out2d_vals, define_out_3d_vals 
!
!
!************************************************************************
!


CONTAINS

!========================================================================

SUBROUTINE define_out0d_vals(outdat,novars,outvars,ivarid,numvar,oopt)
!************************************************************************
!
! *define_out0d_vals* Return 0-D output data values using the variable key ids
!                     and output operator ids
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)model_output.f90  V2.11.2
! Description - either outvars or ivarid must be present but not both 
!             - numvar may only be present if ivarid is present
!             - oopt may only be present if ivarid is present
!             - the following operations can be performed depending on the
!               value of oopt
!                oopt_max  = domain (area) maximum of a 3-D (2-D) variable
!                oopt_min  = domain (area) minimum of a 3-D (2-D) variable
!                oopt_mean = domain (area) mean of a 3-D (2-D) variable
!                oopt_int = domain (area) integrated value of a 3-D (2-D)
!                           variable
!
! Module calls - complex_polar, define_bio0d_vals, define_part0d_vals,
!                define_sed0d_vals, energy_force_0d, energy_0d, enstrophy_0d,
!                error_alloc, inquire_var, min_vars, max_vars, sum_vars,
!                Uarr_at_C, Varr_at_C, vector_mag_arr_atc, Warr_at_C
!
!************************************************************************
!
USE currents
USE datatypes
USE density
USE depths
USE diffusion
USE fluxes
USE grid
USE gridpars
USE iopars
USE meteo
USE modids
USE optics
USE physpars
USE switches
USE syspars
USE turbulence
USE wavevars
USE array_interp, ONLY: Uarr_at_C, Varr_at_C, Warr_at_C
USE biology_output, ONLY: define_bio0d_vals
USE diagnostic_routines, ONLY: energy_force_0d, energy_0d, enstrophy_0d
USE error_routines, ONLY: error_alloc
USE math_library, ONLY: complex_polar, vector_mag_arr_atc
USE modvars_routines, ONLY: inquire_var
USE paral_utilities, ONLY: min_vars, max_vars, sum2_vars
USE particle_output, ONLY: define_part0d_vals
USE sediment_output, ONLY: define_sed0d_vals
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: novars
REAL, INTENT(OUT), DIMENSION(novars) :: outdat
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(novars) :: ivarid, numvar, oopt
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(novars) :: outvars

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    Returned output values
!*novars*  INTEGER Number of variables
!*outvars* DERIVED Variable's attributes (if ivarid is not present)
!*ivarid*  INTEGER Variable's key ids (if outvars is not present)
!*numvar*  INTEGER Variable number (sediment module only)
!*oopt*    INTEGER Variable's operator ids (if outvars is not present)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: eflag0d
INTEGER :: ivar, ivaridx, k, nodim, numvarx, ooptx
INTEGER, DIMENSION(4) :: nhdims = 0
REAL :: atot, vtot
REAL, DIMENSION(5) :: ecomps0d
REAL, DIMENSION(6) :: wcomps0d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d, out2d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3d, out3d, volume


procname(pglev+1) = 'define_out0d_vals'
CALL log_timer_in()

eflag0d = .FALSE.

ivar_1000: DO ivar=1,novars

!  ---optional arguments
   IF (PRESENT(outvars)) THEN
      ivaridx = outvars(ivar)%ivarid; numvarx = outvars(ivar)%numvar
      ooptx = outvars(ivar)%oopt
   ELSE
      ivaridx = ivarid(ivar)
      IF (PRESENT(oopt)) THEN
         ooptx = oopt(ivar)
      ELSE
         ooptx = oopt_null
      ENDIF
      IF (PRESENT(numvar)) THEN
         numvarx = numvar(ivar)
      ELSEIF (PRESENT(outvars)) THEN
         numvarx = outvars(ivar)%numvar
      ELSE
         numvarx = 0
      ENDIF
   ENDIF

   IF (ivaridx.EQ.0) CYCLE ivar_1000 

!
!1. Non-physical variables
!-------------------------
!

   SELECT CASE (ivaridx)
!     ---sediment variables
      CASE (MaxModArids+1:MaxModArids+MaxSedArids)
         CALL define_sed0d_vals(outdat(ivar),ivaridx,numvarx,ooptx)
         CYCLE ivar_1000
!     ---biological variables
      CASE (MaxModArids+MaxSedArids+1:MaxModArids+MaxSedArids+MaxBioArids)
         CALL define_bio0d_vals(outdat(ivar),ivaridx,ooptx)
         CYCLE ivar_1000
!     ---particle output
      CASE (MaxModArids+MaxSedArids+MaxBioArids+1:MaxTotArids)
         CALL define_part0d_vals(outdat(ivar),ivaridx,numvarx,ooptx)
         CYCLE ivar_1000
   END SELECT

!
!2. Physical variables
!---------------------
!
!2.1 Without operator
!--------------------
!

   IF (ooptx.EQ.oopt_null) THEN

      IF (.NOT.eflag0d) THEN
         SELECT CASE (ivaridx)
            CASE (iarr_edens0d,iarr_ekin0d,iarr_epot0d,iarr_etot0d)
               CALL energy_0d(iopt_grid_nodim,ecomps0d)
               eflag0d = .TRUE.
         END SELECT
      ENDIF

      SELECT CASE (ivaridx)
         CASE (iarr_dryfac); outdat(ivar) =  1.0 - nowetatc/REAL(noseaatc)
         CASE (iarr_edens0d); outdat(ivar) = ecomps0d(4)
         CASE (iarr_edissip0d)
            CALL energy_force_0d(iopt_grid_nodim,wcomps0d)
            outdat(ivar) = -wcomps0d(6)
         CASE (iarr_ekin0d); outdat(ivar) = ecomps0d(1)
         CASE (iarr_enstr0d); outdat(ivar) = enstrophy_0d(iopt_grid_nodim)
         CASE (iarr_epot0d); outdat(ivar) = SUM(ecomps0d(2:4))
         CASE (iarr_etot0d); outdat(ivar) = ecomps0d(5)
      END SELECT
      
      CYCLE ivar_1000

   ENDIF

!
!2.2 Allocate arrays
!-------------------
!

   CALL inquire_var(ivaridx,nrank=nodim)

   IF (nodim.EQ.2) THEN
      ALLOCATE (out2d(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('out2d',2,(/ncloc,nrloc/),kndrtype)
      out2d = 0.0
      ALLOCATE (array2d(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('array2d',2,(/ncloc,nrloc/),kndrtype)
      array2d = 0.0
   ELSEIF (nodim.EQ.3) THEN
      ALLOCATE (out3d(ncloc,nrloc,nz),STAT=errstat)
      CALL error_alloc('out3d',3,(/ncloc,nrloc,nz/),kndrtype)
      out3d = 0.0
      ALLOCATE (array3d(ncloc,nrloc,nz),STAT=errstat)
      CALL error_alloc('array3d',3,(/ncloc,nrloc,nz/),kndrtype)
      array3d = 0.0
   ENDIF

!
!2.3 Area and volume of the global domain without land
!-----------------------------------------------------
!

   IF (ooptx.EQ.oopt_mean.OR.ooptx.EQ.oopt_int) THEN
      IF (ooptx.EQ.oopt_mean) THEN 
         CALL sum2_vars(garea,atot,nhdims,'C  ',0,commall=.TRUE.)
      ENDIF
      IF (nodim.EQ.3) THEN
         ALLOCATE (volume(ncloc,nrloc,nz),STAT=errstat)
         CALL error_alloc('volume',3,(/ncloc,nrloc,nz/),kndrtype)
         k_310: DO k=1,nz
            volume(:,:,k) = garea*delzatc(1:ncloc,1:nrloc,k)
         ENDDO k_310
         IF (ooptx.EQ.oopt_mean) THEN
            CALL sum2_vars(volume,vtot,nhdims,'C  ',0,commall=.TRUE.)
         ENDIF
      ENDIF
   ENDIF

!
!2.4 Store data
!---------------
!

   SELECT CASE (ivaridx)

!     ---2-D variable
      CASE (iarr_alphatc_fld); out2d = alphatc_fld
      CASE (iarr_airtemp); out2d = airtemp
      CASE (iarr_atmpres); out2d = atmpres(1:ncloc,1:nrloc)
      CASE (iarr_bdragcoefatc); out2d = bdragcoefatc(1:ncloc,1:nrloc)
      CASE (iarr_bfricvel); out2d = SQRT(bstresatc(1:ncloc,1:nrloc))
      CASE (iarr_bfricvel_max); out2d = SQRT(bstresatc_max)
      CASE (iarr_bfricvel_wav); out2d = SQRT(bstresatc_wav)
      CASE (iarr_bstresatc); out2d = density_ref*bstresatc(1:ncloc,1:nrloc)
      CASE (iarr_bstresatc_max); out2d = density_ref*bstresatc_max
      CASE (iarr_bstresatc_wav); out2d = density_ref*bstresatc_wav
      CASE (iarr_cds); out2d = cds(1:ncloc,1:nrloc)
      CASE (iarr_ces); out2d = ces
      CASE (iarr_chs); out2d = chs
      CASE (iarr_cloud_cover); out2d = cloud_cover
      CASE (iarr_depmeanatc)
         WHERE (ABS(depmeanatc(1:ncloc,1:nrloc)-depmean_flag).GT.&
                  & depmean_flag_eps)
            out2d = depmeanatc(1:ncloc,1:nrloc)
         ELSEWHERE
            out2d = double_fill
         END WHERE
      CASE (iarr_deptotatc); out2d = deptotatc(1:ncloc,1:nrloc)
      CASE (iarr_deptotatc_err); out2d = deptotatc_err
      CASE (iarr_fwave); out2d = fwave
      CASE (iarr_evapminprec); out2d = evapminprec
      CASE (iarr_hdifcoef2d_mom);
         out2d = hdifmom_fac*hdifcoef2datc(1:ncloc,1:nrloc)&
               & /deptotatc(1:ncloc,1:nrloc)
      CASE (iarr_hdifcoef2d_scal);
         out2d = hdifscal_fac*hdifcoef2datc(1:ncloc,1:nrloc)&
               & /deptotatc(1:ncloc,1:nrloc)
      CASE (iarr_hdvelmag)
         CALL vector_mag_arr_atc(udvel(1:ncloc+1,1:nrloc),&
                               & vdvel(1:ncloc,1:nrloc+1),&
                               & 1,1,1,1,ivaridx,.FALSE.,vecmag=out2d)
      CASE (iarr_hmbwdissipmag)
         CALL complex_polar(umbwdissipatc(1:ncloc,:),vmbwdissipatc(:,1:nrloc),&
                          & xamp=array2d,maskvals=maskatc_int)
         WHERE (maskatc_int)
            out2d = array2d/deptotatc(1:ncloc,1:nrloc)
         END WHERE
      CASE (iarr_hmstokesmag)
         CALL complex_polar(umstokesatc(1:ncloc,1:nrloc),&
                          & vmstokesatc(1:ncloc,1:nrloc),xamp=out2d,&
                          & maskvals=maskatc_int)
      CASE (iarr_hmswdissipmag)
         CALL complex_polar(umswdissipatc(1:ncloc,:),vmswdissipatc(:,1:nrloc),&
                          & xamp=array2d,maskvals=maskatc_int)
         WHERE (maskatc_int)
            out2d = array2d/deptotatc(1:ncloc,1:nrloc)
         END WHERE
      CASE (iarr_hmvelmag)
         CALL vector_mag_arr_atc(umvel(1:ncloc+1,1:nrloc),&
                               & vmvel(1:ncloc,1:nrloc+1),&
                               & 1,1,1,1,ivaridx,.FALSE.,vecmag=out2d)
      CASE (iarr_hmveltotmag)
         CALL vector_mag_arr_atc(umvel(1:ncloc+1,1:nrloc),&
                               & vmvel(1:ncloc,1:nrloc+1),&
                               & 1,1,1,1,ivaridx,.FALSE.,vecmag=array2d)
         CALL complex_polar(umstokesatc(1:ncloc,1:nrloc),&
                          & vmstokesatc(1:ncloc,1:nrloc),xamp=out2d,&
                          & maskvals=maskatc_int)
         out2d = array2d + out2d
      CASE (iarr_precipitation); out2d = precipitation
      CASE (iarr_qlatent); out2d = qlatent
      CASE (iarr_qlwave); out2d = qlwave
      CASE (iarr_qnonsol); out2d = qnonsol
      CASE (iarr_qrad); out2d = qrad
      CASE (iarr_qsensible); out2d = qsensible
      CASE (iarr_qtot); out2d = qrad - qnonsol
      CASE (iarr_relhum); out2d = relhum
      CASE (iarr_ssalflux); out2d = ssalflux
      CASE (iarr_sfricatc); out2d = SQRT(sstresatc)
      CASE (iarr_sstresatc); out2d = density_ref*sstresatc
      CASE (iarr_ubstresatc)
         CALL Uarr_at_C(ubstresatu,out2d,1,1,(/1,1,1/),(/ncloc+1,nrloc,1/),1,&
                      & ivaridx,.FALSE.)
         out2d  = density_ref*out2d
      CASE (iarr_udstokesatu)
         CALL Uarr_at_C(udstokesatu(1:ncloc+1,:),out2d,1,1,(/1,1,nz/),&
                     & (/ncloc+1,nrloc,nz/),1,ivaridx,.FALSE.)
      CASE (iarr_udvel)
         CALL Uarr_at_C(udvel(1:ncloc+1,1:nrloc),out2d,1,1,(/1,1,nz/),&
                     & (/ncloc+1,nrloc,nz/),1,ivaridx,.FALSE.)
      CASE (iarr_umbwdissipatc)
         WHERE (maskatc_int)
         out2d = umbwdissipatc(1:ncloc,:)/deptotatc(1:ncloc,1:nrloc)
      END WHERE
      CASE (iarr_umstokesatc)
         out2d = umstokesatc(1:ncloc,1:nrloc)
      CASE (iarr_umswdissipatc)
         WHERE (maskatc_int)
            out2d = umswdissipatc(1:ncloc,:)
         END WHERE
      CASE (iarr_umvel)
         CALL Uarr_at_C(umvel(1:ncloc+1,1:nrloc),out2d,1,1,(/1,1,nz/),&
                     & (/ncloc+1,nrloc,nz/),1,ivaridx,.FALSE.)
      CASE (iarr_umveltot)
         CALL Uarr_at_C(umvel(1:ncloc+1,1:nrloc),out2d,1,1,(/1,1,nz/),&
                     & (/ncloc+1,nrloc,nz/),1,ivaridx,.FALSE.)
         out2d = out2d + umstokesatc(1:ncloc,1:nrloc)
      CASE (iarr_usstresatc)
         CALL Uarr_at_C(usstresatu,out2d,1,1,(/1,1,1/),(/ncloc+1,nrloc,1/),1,&
                      & ivaridx,.FALSE.)
         out2d  = density_ref*out2d
      CASE (iarr_uwindatc); out2d = uwindatc(1:ncloc,1:nrloc)
      CASE (iarr_vbstresatc)
         CALL Varr_at_C(vbstresatv,out2d,1,1,(/1,1,1/),(/ncloc,nrloc+1,1/),1,&
                      & ivaridx,.FALSE.)
         out2d = density_ref*out2d
      CASE (iarr_vdstokesatv)
         CALL Varr_at_C(vdstokesatv(:,1:nrloc+1),out2d,1,1,(/1,1,nz/),&
                     & (/ncloc,nrloc+1,nz/),1,ivaridx,.FALSE.)
      CASE (iarr_vdvel)
         CALL Varr_at_C(vdvel(1:ncloc,1:nrloc+1),out2d,1,1,(/1,1,nz/),&
                     & (/ncloc,nrloc+1,nz/),1,ivaridx,.FALSE.)
      CASE (iarr_vmbwdissipatc)
         WHERE (maskatc_int)
            out2d = vmbwdissipatc(:,1:nrloc)
         END WHERE
      CASE (iarr_vmstokesatc)
         out2d = vmstokesatc(1:ncloc,1:nrloc)
      CASE (iarr_vmswdissipatc)
         WHERE (maskatc_int)
            out2d = vmswdissipatc(:,1:nrloc)
         END WHERE
      CASE (iarr_vmvel)
         CALL Varr_at_C(vmvel(1:ncloc,1:nrloc+1),out2d,1,1,(/1,1,nz/),&
                     & (/ncloc,nrloc+1,nz/),1,ivaridx,.FALSE.)
      CASE (iarr_vmveltot)
         CALL Varr_at_C(vmvel(1:ncloc,1:nrloc+1),out2d,1,1,(/1,1,nz/),&
                     & (/ncloc,nrloc+1,nz/),1,ivaridx,.FALSE.)
         out2d = out2d + vmstokesatc(1:ncloc,1:nrloc)
      CASE (iarr_vsstresatc)
         CALL Varr_at_C(vsstresatv,out2d,1,1,(/1,1,1/),(/ncloc,nrloc+1,1/),1,&
                      & ivaridx,.FALSE.)
         out2d  = density_ref*out2d
      CASE (iarr_vwindatc); out2d = vwindatc(1:ncloc,1:nrloc)
      CASE (iarr_wavedir); out2d = radtodeg*wavedir(1:ncloc,1:nrloc)
      CASE (iarr_waveexcurs); out2d = waveexcurs(1:ncloc,1:nrloc)
      CASE (iarr_wavefreq); out2d = wavefreq(1:ncloc,1:nrloc)
      CASE (iarr_waveheight); out2d = waveheight(1:ncloc,1:nrloc)
      CASE (iarr_wavenum); out2d = wavenum(1:ncloc,1:nrloc)
      CASE (iarr_waveperiod); out2d = waveperiod(1:ncloc,1:nrloc)
      CASE (iarr_wavepres); out2d = wavepres(1:ncloc,1:nrloc)
      CASE (iarr_wavethickatc); out2d = wavethickatc
      CASE (iarr_wavevel); out2d = wavevel(1:ncloc,1:nrloc)
      CASE (iarr_windatc)
         CALL complex_polar(uwindatc(1:ncloc,1:nrloc),&
                          & vwindatc(1:ncloc,1:nrloc),xamp=out2d,&
                          & maskvals=maskatc_int)
      CASE (iarr_zeta); out2d = zeta(1:ncloc,1:nrloc)
      CASE (iarr_zroughatc); out2d = zroughatc(1:ncloc,1:nrloc)
      CASE (iarr_zaroughatc); out2d = zaroughatc(1:ncloc,1:nrloc)

!     ---3-D variable
      CASE (iarr_beta_sal); out3d = beta_sal(1:ncloc,1:nrloc,:)
      CASE (iarr_beta_temp); out3d = beta_temp(1:ncloc,1:nrloc,:)
      CASE (iarr_buofreq2)
         out3d(:,:,1) = buofreq2(:,:,2)
         out3d(:,:,2:nz-1) = 0.5*(buofreq2(:,:,2:nz-1)+buofreq2(:,:,3:nz))
         out3d(:,:,nz) = buofreq2(:,:,nz)
      CASE (iarr_dens); out3d = dens(1:ncloc,1:nrloc,:) - 1000.0
      CASE (iarr_dissip)
         out3d(:,:,1) = dissip(1:ncloc,1:nrloc,2)
         out3d(:,:,2:nz-1) = 0.5*(dissip(1:ncloc,1:nrloc,2:nz-1)+&
                                & dissip(1:ncloc,1:nrloc,3:nz))
         out3d(:,:,nz) = dissip(1:ncloc,1:nrloc,nz)
      CASE (iarr_hbwdissipmag)
         CALL complex_polar(ubwdissipatc(1:ncloc,:,:),&
                          & vbwdissipatc(:,1:nrloc,:),xamp=out3d,&
                          & maskvals=maskatc_int)
      CASE (iarr_hdifcoef3d_mom)
         out3d = hdifmom_fac*hdifcoef3datc(1:ncloc,1:nrloc,:)
      CASE (iarr_hdifcoef3d_scal)
         out3d = hdifscal_fac*hdifcoef3datc(1:ncloc,1:nrloc,:)
      CASE (iarr_hstokesmag)
         CALL complex_polar(ustokesatc(1:ncloc,1:nrloc,:),&
                          & vstokesatc(1:ncloc,1:nrloc,:),xamp=out3d,&
                          & maskvals=maskatc_int)
      CASE (iarr_hswdissipmag)
         CALL complex_polar(uswdissipatc(1:ncloc,:,:),&
                          & vswdissipatc(:,1:nrloc,:),xamp=out3d,&
                          & maskvals=maskatc_int)
      CASE (iarr_hvelmag)
         CALL vector_mag_arr_atc(uvel(1:ncloc+1,1:nrloc,:),&
                               & vvel(1:ncloc,1:nrloc+1,:),&
                               & 1,1,nz,1,ivaridx,.FALSE.,vecmag=out3d)
      CASE (iarr_hveltotmag)
         CALL vector_mag_arr_atc(uvel(1:ncloc+1,1:nrloc,:),&
                               & vvel(1:ncloc,1:nrloc+1,:),&
                               & 1,1,1,1,ivaridx,.FALSE.,vecmag=array3d)
         CALL complex_polar(ustokesatc(1:ncloc,1:nrloc,:),&
                          & vstokesatc(1:ncloc,1:nrloc,:),xamp=out3d,&
                          & maskvals=maskatc_int)
         out3d = array3d + out3d
      CASE (iarr_kinvisc)
         out3d = kinvisc(1:ncloc,1:nrloc,:)
      CASE (iarr_radiance)
         out3d = 0.5*(radiance(:,:,1:nz)+radiance(:,:,2:nz+1))
      CASE (iarr_sal); out3d = sal(1:ncloc,1:nrloc,:)
      CASE (iarr_shearfreq2)
         out3d(:,:,1) = shearfreq2(:,:,2)
         out3d(:,:,2:nz-1) = 0.5*(shearfreq2(:,:,2:nz-1)+shearfreq2(:,:,3:nz))
         out3d(:,:,nz) = shearfreq2(:,:,nz)
      CASE (iarr_temp); out3d = temp(1:ncloc,1:nrloc,:)
      CASE (iarr_tke)
         IF (iopt_turb_tke_bbc.EQ.1) THEN
            out3d(:,:,1) = tke(1:ncloc,1:nrloc,2)
         ELSE
            out3d(:,:,1) = 0.5*(tke(1:ncloc,1:nrloc,1)+tke(1:ncloc,1:nrloc,2))
         ENDIF
         out3d(:,:,2:nz-1) = 0.5*(tke(1:ncloc,1:nrloc,2:nz-1)+&
                                & tke(1:ncloc,1:nrloc,3:nz))
         IF (iopt_turb_tke_sbc.EQ.1) THEN
            out3d(:,:,nz) = tke(1:ncloc,1:nrloc,nz)
         ELSE
            out3d(:,:,nz) = 0.5*(tke(1:ncloc,1:nrloc,nz)+&
                               & tke(1:ncloc,1:nrloc,nz+1))
         ENDIF
      CASE (iarr_ubwdissipatc)
         out3d = ubwdissipatc(1:ncloc,:,:)
      CASE (iarr_ustokesatc)
         out3d = ustokesatc(1:ncloc,1:nrloc,:)
      CASE (iarr_uswdissipatc)
         out3d = uswdissipatc(1:ncloc,:,:)
      CASE (iarr_uvel)
         CALL Uarr_at_C(uvel(1:ncloc+1,1:nrloc,:),out3d,1,1,(/1,1,1/),&
                     & (/ncloc+1,nrloc,nz/),1,ivaridx,.FALSE.)
      CASE (iarr_uveltot)
         CALL Uarr_at_C(uvel(1:ncloc+1,1:nrloc,k),out3d,1,1,(/1,1,1/),&
                     & (/ncloc+1,nrloc,nz/),1,ivaridx,.FALSE.)
         out3d = out3d + ustokesatc(1:ncloc,1:nrloc,:)
      CASE (iarr_vbwdissipatc)
         out3d = vbwdissipatc(:,1:nrloc,:)
      CASE (iarr_vdifcoefmom)
         out3d = 0.5*(vdifcoefmom(1:ncloc,1:nrloc,1:nz)+&
                    & vdifcoefmom(1:ncloc,1:nrloc,2:nz+1))
      CASE (iarr_vdifcoefscal)
         out3d = 0.5*(vdifcoefscal(:,:,1:nz)+vdifcoefscal(:,:,2:nz+1))
      CASE (iarr_vdifcoefscal_rot)
         out3d = 0.5*(vdifcoefscal_rot(:,:,1:nz)+vdifcoefscal_rot(:,:,2:nz+1))
      CASE (iarr_vdifcoefscal_norot)
         out3d = 0.5*(vdifcoefscal(:,:,1:nz)-vdifcoefscal_rot(:,:,1:nz)+&
                    & vdifcoefscal(:,:,2:nz+1)-vdifcoefscal_rot(:,:,2:nz+1))
      CASE (iarr_vdifcoeftke)
         out3d = 0.5*(vdifcoeftke(:,:,1:nz)+vdifcoeftke(:,:,2:nz+1))
      CASE (iarr_vstokesatc)
         out3d = vstokesatc(1:ncloc,1:nrloc,:)
      CASE (iarr_vswdissipatc)
         out3d = vswdissipatc(:,1:nrloc,:)
      CASE (iarr_vvel)
         CALL Varr_at_C(vvel(1:ncloc,1:nrloc+1,:),out3d,1,1,(/1,1,1/),&
                     & (/ncloc,nrloc+1,nz/),1,ivaridx,.FALSE.)
      CASE (iarr_vveltot)
         CALL Varr_at_C(vvel(1:ncloc,1:nrloc+1,:),out3d,1,1,(/1,1,1/),&
                     & (/ncloc,nrloc+1,nz/),1,ivaridx,.FALSE.)
         out3d = out3d + vstokesatc(1:ncloc,1:nrloc,:)
      CASE (iarr_wphys); out3d = wphys
      CASE (iarr_wstokesatw)
         CALL Warr_at_C(wstokesatw(1:ncloc,1:nrloc,:),out3d,(/1,1,1/),&
                     & (/ncloc,nrloc,nz+1/),1,ivaridx,.FALSE.)
      CASE (iarr_wvel)
         CALL Warr_at_C(wvel(1:ncloc,1:nrloc,:),out3d,(/1,1,1/),&
                     & (/ncloc,nrloc,nz+1/),1,ivaridx,.FALSE.)
      CASE (iarr_zlmix)
         out3d = 0.5*(zlmix(1:ncloc,1:nrloc,1:nz)+zlmix(1:ncloc,1:nrloc,2:nz+1))

   END SELECT

!
!2.5 Apply operator
!------------------
!

   IF (nodim.EQ.2) THEN

      SELECT CASE (ooptx)

         CASE (oopt_min)
            CALL min_vars(out2d,outdat(ivar),ivaridx,commall=.TRUE.,&
                        & mask=maskatc_int)
         CASE (oopt_max)
            CALL max_vars(out2d,outdat(ivar),ivaridx,commall=.TRUE.,&
                        & mask=maskatc_int)
         CASE (oopt_mean)
            CALL sum2_vars(garea*out2d,outdat(ivar),nhdims,'C  ',ivaridx,&
                         & commall=.TRUE.)
            outdat(ivar) = outdat(ivar)/atot

         CASE (oopt_int)
            CALL sum2_vars(garea*out2d,outdat(ivar),nhdims,'C  ',ivaridx,&
                         & commall=.TRUE.)

      END SELECT

   ELSEIF (nodim.EQ.3) THEN

      SELECT CASE (ooptx)

         CASE (oopt_min)
            CALL min_vars(out3d,outdat(ivar),ivaridx,commall=.TRUE.,&
                        & mask=maskatc_int)
         CASE (oopt_max)
            CALL max_vars(out3d,outdat(ivar),ivaridx,commall=.TRUE.,&
                        & mask=maskatc_int)
         CASE (oopt_mean)
            CALL sum2_vars(volume*out3d,outdat(ivar),nhdims,'C  ',ivaridx,&
                         & commall=.TRUE.)
            outdat(ivar) = outdat(ivar)/vtot

         CASE (oopt_int)
            CALL sum2_vars(volume*out3d,outdat(ivar),nhdims,'C  ',ivaridx,&
                         & commall=.TRUE.)

      END SELECT

   ENDIF

!
!2.6 Deallocate arrays
!---------------------
!

   IF (ALLOCATED(out2d)) DEALLOCATE (array2d,out2d)
   IF (ALLOCATED(out3d)) DEALLOCATE (array3d,out3d)
   IF (ALLOCATED(volume)) DEALLOCATE (volume)

ENDDO ivar_1000

CALL log_timer_out()


RETURN

END SUBROUTINE define_out0d_vals

!========================================================================

SUBROUTINE define_out2d_vals(outdat,i,j,novars,outvars,ivarid,numvar,oopt,klev,&
                           & dep,node)
!************************************************************************
!
! *define_out2d_vals* Return 2-D output data values using the variable key ids
!                     and output operator ids
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)model_output.f90  V2.11.2
!
! Description - either outvars or ivarid need to be present but not both
!             - numvar may only be present if ivarid is present
!             - oopt may only be present if ivarid is present
!             - klev must be present if ivarid is present and
!               oopt = oopt_klev for at least one variable
!             - dep must be present if ivarid is present and
!               oopt = oopt_kdep for at least one variable
!             - node is the vertical grid location for the data: 'C' (default)
!               or 'W'
!             - the following operations can be performed depending on the
!               value of oopt
!                oopt_max  = vertical maximum of a 3-D variable
!                oopt_min  = vertical minimum of a 3-D variable
!                oopt_mean = vertical mean of a 3-D variable
!                oopt_int = vertically integrated value of a 3-D variable
!                oopt_klev = 3-D variable at a given vertical grid level
!                oopt_dep  = 3-D variable at a depth below the free surface
!
! Module calls - complex_polar, define_bio2d_vals, define_bio3d_vals,
!                define_part_2d_vals, define_part_3d_vals, define_sed_2dvals,
!                define_sed3d_vals, energy_flux_2d, energy_flux_3d,
!                energy_force_2d, energy_force_3d, energy_2d, energy_3d,
!                intpol1d_model_to_dep, inquire_var, richardson_number_0d,
!                Uarr_at_C, Uvar_at_C, vector_mag_arr_atc, vector_mag_var_atc,
!                vorticity_2d, Varr_at_C, Vvar_at_C, Wvar_at_C
!
!************************************************************************
!
USE currents
USE datatypes
USE density
USE depths
USE diffusion
USE fluxes
USE grid
USE gridpars
USE iopars
USE meteo
USE modids
USE optics
USE physpars
USE switches
USE syspars
USE timepars
USE turbulence
USE wavevars
USE array_interp, ONLY: Uarr_at_C, Uvar_at_C, Varr_at_C, Vvar_at_C, Wvar_at_C
USE biology_output, ONLY: define_bio2d_vals, define_bio3d_vals
USE diagnostic_routines, ONLY: energy_flux_2d, energy_flux_3d, &
                             & energy_force_2d, energy_force_3d, energy_2d, &
                             & energy_3d, vorticity_2d, vorticity_3d
USE grid_interp, ONLY: intpol1d_model_to_dep
USE math_library, ONLY: complex_polar, vector_mag_arr_atc, vector_mag_var_atc
USE modvars_routines, ONLY: inquire_var
USE particle_output, ONLY: define_part2d_vals, define_part3d_vals
USE sediment_output, ONLY: define_sed2d_vals, define_sed3d_vals
USE turbulence_routines, ONLY: richardson_number_0d

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, novars
REAL, INTENT(OUT), DIMENSION(novars) :: outdat
CHARACTER (LEN=lennode), OPTIONAL, DIMENSION(novars) :: node
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(novars) :: ivarid, klev, numvar, oopt
REAL, INTENT(IN), OPTIONAL, DIMENSION(novars) :: dep
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(novars) :: outvars

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    Returned output values
!*i*       INTEGER (Local) X-index of the data location
!*j*       INTEGER (Local) Y-index of the data location
!*novars*  INTEGER Number of variables
!*outvars* DERIVED Variable's attributes (if ivarid is not present)
!*ivarid*  INTEGER Variable's key ids (if outvars is not present)
!*numvar*  INTEGER Variable number (sediment module only)
!*oopt*    INTEGER Variable's operator ids (if outvars is not present)
!*klev*    INTEGER Variable's vertical level attributes (if outvars is not
!                  present)
!*dep*     REAL    Variable's depth attribute (if outvars is not present)    [m]
!*node*    CHAR    Vertical grid node of the output data ('C', 'W')
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: eflag2d, eflag3d, fflag2d, fflag3d
CHARACTER (LEN=lennode) :: cnode, nodex
INTEGER :: ivar, ivaridx, k, kc, numvarx, nzmax, nzmin, ooptx
REAL :: depx, varx
REAL, DIMENSION(4) :: ecomps2d
REAL, DIMENSION(3) :: ecomps3d
REAL, DIMENSION(6) :: wcomps
REAL, DIMENSION(4,2) :: fcomps2d
REAL, DIMENSION(4,3) :: fcomps3d
REAL, DIMENSION(2) :: out1d
REAL, DIMENSION(nz) :: out1dc
REAL, DIMENSION(nz+1) :: out1dw
REAL, DIMENSION(1,1,nz) :: out3dc


!---initialise
IF (.NOT.maskatc_int(i,j)) THEN
   outdat = 0.0
   RETURN
ENDIF

eflag2d = .FALSE.; eflag3d = .FALSE.
fflag2d = .FALSE.; fflag3d = .FALSE.

ivar_200: DO ivar=1,novars

!  ---optional arguments
   IF (PRESENT(outvars)) THEN
      ivaridx = outvars(ivar)%ivarid; numvarx = outvars(ivar)%numvar
      ooptx = outvars(ivar)%oopt; nodex = TRIM(outvars(ivar)%node)
      IF (ooptx.EQ.oopt_klev) kc = outvars(ivar)%klev
      IF (ooptx.EQ.oopt_dep) depx = outvars(ivar)%dep
   ELSE
      ivaridx = ivarid(ivar)
      IF (PRESENT(oopt)) THEN
         ooptx = oopt(ivar)
      ELSE
         ooptx = oopt_null
      ENDIF
      IF (ooptx.EQ.oopt_klev) kc = klev(ivar)
      IF (ooptx.EQ.oopt_dep) depx = dep(ivar)
      IF (PRESENT(node)) THEN
         nodex = node(ivar)
      ELSE
         nodex = 'C'
      ENDIF
      IF (PRESENT(numvar)) THEN
         numvarx = numvar(ivar)
      ELSEIF (PRESENT(outvars)) THEN
         numvarx = outvars(ivar)%numvar
      ELSE
         numvarx = 0
      ENDIF
   ENDIF

   IF (ivaridx.EQ.0) CYCLE ivar_200

!
!1. Non-physical variables
!-------------------------
!

   SELECT CASE (ivaridx)
!     ---sediment variables
      CASE (MaxModArids+1:MaxModArids+MaxSedArids)
         CALL define_sed2d_vals(outdat(ivar),i,j,ivaridx,numvarx,ooptx,kc,depx)
         CYCLE ivar_200
!     ---biological variables
      CASE (MaxModArids+MaxSedArids+1:MaxModArids+MaxSedArids+MaxBioArids)
         CALL define_bio2d_vals(outdat(ivar),i,j,ivaridx,ooptx,kc,depx)
         CYCLE ivar_200
!     ---particle model variables
      CASE (MaxModArids+MaxSedArids+MaxBioArids+1:MaxtotArids)
         CALL define_part2d_vals(outdat(ivar),i,j,ivaridx,numvarx,ooptx,kc,depx)
         CYCLE ivar_200
   END SELECT

!
!2. Physical variables
!---------------------
!
!2.1 Energy calls
!----------------
!

   SELECT CASE (ivaridx)
      CASE (iarr_edens2d,iarr_ekin2d,iarr_epot2d,iarr_etot2d)
         IF (.NOT.eflag2d) THEN
            CALL energy_2d(i,j,iopt_grid_nodim,ecomps2d)
            eflag2d = .TRUE.
         ENDIF
      CASE (iarr_eflux2du,iarr_eflux2dv)
         IF (.NOT.fflag2d) THEN
            CALL energy_flux_2d(i,j,iopt_grid_nodim,fcomps2d)
            fflag2d = .TRUE.
         ENDIF
   END SELECT

!
!2.2 Without operator
!--------------------
!

   IF (ooptx.EQ.oopt_null) THEN

      SELECT CASE (ivaridx)

         CASE (iarr_alphatc_fld); outdat(ivar) = alphatc_fld(i,j)
         CASE (iarr_airtemp); outdat(ivar) = airtemp(i,j)
         CASE (iarr_atmpres); outdat(ivar) = atmpres(i,j)
         CASE (iarr_bdragcoefatc); outdat(ivar) = bdragcoefatc(i,j)
         CASE (iarr_bfricvel); outdat(ivar) = SQRT(bstresatc(i,j))
         CASE (iarr_bfricvel_max); outdat(ivar) = SQRT(bstresatc_max(i,j))
         CASE (iarr_bfricvel_wav); outdat(ivar) = SQRT(bstresatc_wav(i,j))
         CASE (iarr_bstresatc); outdat(ivar) = density_ref*bstresatc(i,j)
         CASE (iarr_bstresatc_max)
            outdat(ivar) = density_ref*bstresatc_max(i,j)
         CASE (iarr_bstresatc_wav)
            outdat(ivar) = density_ref*bstresatc_wav(i,j)
         CASE (iarr_cds); outdat(ivar) = cds(i,j)
         CASE (iarr_ces); outdat(ivar) = ces(i,j)
         CASE (iarr_chs); outdat(ivar) = chs(i,j)
         CASE (iarr_cloud_cover); outdat(ivar) = cloud_cover(i,j)
         CASE (iarr_depmeanatc)
            outdat(ivar) = depmeanatc(i,j)
         CASE (iarr_deptotatc); outdat(ivar) = deptotatc(i,j)
         CASE (iarr_deptotatc_err); outdat(ivar) = deptotatc_err(i,j)
         CASE (iarr_edens2d)
            outdat(ivar) = ecomps2d(3)
         CASE (iarr_edissip2d)
            CALL energy_force_2d(i,j,iopt_grid_nodim,wcomps)
            outdat(ivar) = -wcomps(6)
         CASE (iarr_eflux2du)
            outdat(ivar) = fcomps2d(4,1)
         CASE (iarr_eflux2dv)
            outdat(ivar) = fcomps2d(4,2)
         CASE (iarr_ekin2d)
            outdat(ivar) = ecomps2d(1)
         CASE (iarr_epot2d)
            outdat(ivar) = ecomps2d(2) + ecomps2d(3)
         CASE (iarr_etot2d)
            outdat(ivar) = ecomps2d(4)
         CASE (iarr_fwave); outdat(ivar)= fwave(i,j)
         CASE (iarr_evapminprec); outdat(ivar) = evapminprec(i,j)
         CASE (iarr_hmbwdissipmag)
            CALL complex_polar(umbwdissipatc(i,j),vmbwdissipatc(i,j),&
                             & xamp=outdat(ivar))
         CASE (iarr_hdifcoef2d_mom)
            outdat(ivar) = hdifmom_fac*hdifcoef2datc(i,j)/deptotatc(i,j)
         CASE (iarr_hdifcoef2d_scal)
            outdat(ivar) = hdifscal_fac*hdifcoef2datc(i,j)/deptotatc(i,j)
         CASE (iarr_hdvelmag)
            CALL vector_mag_var_atc(udvel(i:i+1,j),vdvel(i,j:j+1),i,j,1,1,1,&
                                  & vecmag=outdat(ivar))
         CASE (iarr_hmstokesmag)
         CALL complex_polar(umstokesatc(i,j),vmstokesatc(i,j),&
                          & xamp=outdat(ivar))
         CASE (iarr_hmswdissipmag)
            CALL complex_polar(umswdissipatc(i,j),vmswdissipatc(i,j),&
                             & xamp=outdat(ivar))
         CASE (iarr_hmvelmag)
            CALL vector_mag_var_atc(umvel(i:i+1,j),vmvel(i,j:j+1),i,j,1,1,1,&
                 & vecmag=outdat(ivar))
         CASE (iarr_hmvelmag_hadv_cfl)
            CALL vector_mag_var_atc(umvel(i:i+1,j),vmvel(i,j:j+1),i,j,1,1,1,&
                                  & vecmag=outdat(ivar))
            outdat(ivar) = ABS(outdat(ivar))*&
                            & (delt2d/MIN(delxatc(i,j),delyatc(i,j)))
         CASE (iarr_hmveltotmag)
            CALL vector_mag_var_atc(umvel(i:i+1,j),vmvel(i,j:j+1),i,j,1,1,1,&
                                  & vecmag=outdat(ivar))
            CALL complex_polar(umstokesatc(i,j),vmstokesatc(i,j),xamp=varx)
            outdat(ivar) = outdat(ivar) + varx
         CASE (iarr_precipitation); outdat(ivar) = precipitation(i,j)
         CASE (iarr_p2dbcgradatu)
            outdat(ivar) = Uvar_at_C(p2dbcgradatu(i:i+1,j),i,j,nz,2,1)
         CASE (iarr_p2dbcgradatv)
            outdat(ivar) = Vvar_at_C(p2dbcgradatv(i,j:j+1),i,j,nz,2,1)
         CASE (iarr_qlatent); outdat(ivar) = qlatent(i,j)
         CASE (iarr_qlwave); outdat(ivar) = qlwave(i,j)
         CASE (iarr_qnonsol); outdat(ivar) = qnonsol(i,j)
         CASE (iarr_qrad); outdat(ivar) = qrad(i,j)
         CASE (iarr_qsensible); outdat(ivar) = qsensible(i,j)
         CASE (iarr_qtot); outdat(ivar) = qrad(i,j) - qnonsol(i,j)
         CASE (iarr_relhum); outdat(ivar) = relhum(i,j)
         CASE (iarr_sfricatc); outdat(ivar) = SQRT(sstresatc(i,j))
         CASE (iarr_ssalflux); outdat(ivar) = ssalflux(i,j)
         CASE (iarr_sstresatc); outdat(ivar) = density_ref*sstresatc(i,j)
         CASE (iarr_ubstresatc)
            outdat(ivar) = density_ref*Uvar_at_C(ubstresatu(i:i+1,j),i,j,nz,1,1)
         CASE (iarr_udstokesatu)
            outdat(ivar) = Uvar_at_C(udstokesatu(i:i+1,j),i,j,nz,1,1)
         CASE (iarr_udvel); outdat(ivar) = Uvar_at_C(udvel(i:i+1,j),i,j,nz,1,1)
         CASE (iarr_umbwdissipatc)
            outdat(ivar) = umbwdissipatc(i,j)
         CASE (iarr_umstokesatc); outdat(ivar) = umstokesatc(i,j)
         CASE (iarr_umswdissipatc)
            outdat(ivar) = umswdissipatc(i,j)
         CASE (iarr_umvel); outdat(ivar) = Uvar_at_C(umvel(i:i+1,j),i,j,nz,1,1)
         CASE (iarr_umvel_hadv_cfl)
            outdat(ivar) = Uvar_at_C(umvel(i:i+1,j),i,j,nz,1,1)
            outdat(ivar) = ABS(outdat(ivar))*(delt2d/delxatc(i,j))
         CASE (iarr_umveltot)
            outdat(ivar) = Uvar_at_C(umvel(i:i+1,j),i,j,nz,1,1)
            outdat(ivar) = outdat(ivar) + umstokesatc(i,j)
         CASE (iarr_usstresatc)
            outdat(ivar) = density_ref*Uvar_at_C(usstresatu(i:i+1,j),i,j,nz,1,1)
         CASE (iarr_uwindatc); outdat(ivar) = uwindatc(i,j)
         CASE (iarr_vbstresatc)
            outdat(ivar) = density_ref*Vvar_at_C(vbstresatv(i,j:j+1),i,j,nz,1,1)
         CASE (iarr_vdstokesatv)
            outdat(ivar) = Vvar_at_C(vdstokesatv(i,j:j+1),i,j,nz,1,1)
         CASE (iarr_vdvel); outdat(ivar) = Vvar_at_C(vdvel(i,j:j+1),i,j,nz,1,1)
         CASE (iarr_vmbwdissipatc)
            outdat(ivar) = vmbwdissipatc(i,j)
         CASE (iarr_vmstokesatc); outdat(ivar) = vmstokesatc(i,j)
         CASE (iarr_vmswdissipatc)
            outdat(ivar) = vmswdissipatc(i,j)
         CASE (iarr_vmvel); outdat(ivar) = Vvar_at_C(vmvel(i,j:j+1),i,j,nz,1,1)
         CASE (iarr_vmvel_hadv_cfl)
            outdat(ivar) = Vvar_at_C(vmvel(i,j:j+1),i,j,nz,1,1)
            outdat(ivar) = ABS(outdat(ivar))*(delt2d/delyatc(i,j))
         CASE (iarr_vmveltot)
            outdat(ivar) = Vvar_at_C(vmvel(i,j:j+1),i,j,nz,1,1)
            outdat(ivar) = outdat(ivar) + vmstokesatc(i,j)
         CASE (iarr_vortic2d); outdat(ivar) = vorticity_2d(i,j,iopt_grid_nodim)
         CASE (iarr_vsstresatc)
            outdat(ivar) = density_ref*Vvar_at_C(vsstresatv(i,j:j+1),i,j,nz,1,1)
         CASE (iarr_vwindatc); outdat(ivar) = vwindatc(i,j)
         CASE (iarr_wavedir); outdat(ivar) = radtodeg*wavedir(i,j)
         CASE (iarr_waveexcurs); outdat(ivar) = waveexcurs(i,j)
         CASE (iarr_wavefreq); outdat(ivar) = wavefreq(i,j)
         CASE (iarr_waveheight); outdat(ivar) = waveheight(i,j)
         CASE (iarr_wavenum); outdat(ivar) = wavenum(i,j)
         CASE (iarr_waveperiod); outdat(ivar) = waveperiod(i,j)
         CASE (iarr_wavepres); outdat(ivar) = wavepres(i,j)
         CASE (iarr_wavethickatc); outdat(ivar) = wavethickatc(i,j)
         CASE (iarr_wavevel); outdat(ivar) = wavevel(i,j)
         CASE (iarr_windatc)
            CALL complex_polar(uwindatc(i,j),vwindatc(i,j),&
                             & xamp=outdat(ivar))
         CASE (iarr_zeta); outdat(ivar) = zeta(i,j)
         CASE (iarr_zaroughatc); outdat(ivar) = zaroughatc(i,j)
         CASE (iarr_zroughatc); outdat(ivar) = zroughatc(i,j)
         CASE (iarr_2d_hdif_cfl)
            outdat(ivar) = MAX(hdifmom_fac,hdifscal_fac)*hdifcoef2datc(i,j)*&
                         & (delt2d/MIN(delxatc(i,j),delyatc(i,j))**2)
      END SELECT

!
!2.3 3-D variable at a given level
!---------------------------------
!

   ELSEIF (ooptx.EQ.oopt_klev) THEN

!     ---energy calls
      SELECT CASE (ivaridx)
         CASE (iarr_edens3d,iarr_ekin3d,iarr_etot3d)
            IF (.NOT.eflag3d) THEN
               CALL energy_3d(i,j,kc,ecomps3d)
               eflag3d = .TRUE.
            ENDIF
         CASE (iarr_eflux3du,iarr_eflux3dv,iarr_eflux3dw)
            IF (.NOT.fflag3d) THEN
               CALL energy_flux_3d(i,j,kc,fcomps3d)
               fflag3d = .TRUE.
            ENDIF
      END SELECT

!
!2.3.1 C-node level
!------------------
!

      IF (TRIM(nodex).EQ.'C') THEN 

         SELECT CASE(ivaridx)
  
            CASE (iarr_beta_sal); outdat(ivar) = beta_sal(i,j,kc)
            CASE (iarr_beta_temp); outdat(ivar) = beta_temp(i,j,kc)
            CASE (iarr_buofreq2)
               IF (kc.EQ.1) THEN
                  outdat(ivar) = buofreq2(i,j,2)
               ELSEIF (kc.EQ.nz) THEN
                  outdat(ivar) = buofreq2(i,j,nz)
               ELSE
                  outdat(ivar) = Wvar_at_C(buofreq2(i,j,kc:kc+1),i,j)
               ENDIF
            CASE (iarr_dens); outdat(ivar) = dens(i,j,kc) - 1000.0
            CASE (iarr_dissip)
               IF (kc.EQ.1) THEN
                  outdat(ivar) = dissip(i,j,2)
               ELSEIF (kc.EQ.nz) THEN
                  outdat(ivar) = dissip(i,j,nz)
               ELSE
                  outdat(ivar) = Wvar_at_C(dissip(i,j,kc:kc+1),i,j)
               ENDIF
            CASE (iarr_edens3d)
               outdat(ivar) = ecomps3d(2)
            CASE (iarr_edissip3d)
               CALL energy_force_3d(i,j,kc,wcomps(1:4))
               outdat(ivar) = -wcomps(4)
            CASE (iarr_eflux3du)
               outdat(ivar) = fcomps3d(4,1)
            CASE (iarr_eflux3dv)
               outdat(ivar) = fcomps3d(4,2)
            CASE (iarr_eflux3dw)
               outdat(ivar) = fcomps3d(4,3)
            CASE (iarr_ekin3d)
               outdat(ivar) = ecomps3d(1)
            CASE (iarr_etot3d)
               outdat(ivar) = ecomps3d(3)
            CASE (iarr_hbwdissipmag)
               CALL complex_polar(ubwdissipatc(i,j,kc),vbwdissipatc(i,j,kc),&
                                & xamp=outdat(ivar))
            CASE (iarr_hdifcoef3d_mom)
               outdat(ivar) = hdifmom_fac*hdifcoef3datc(i,j,kc)
            CASE (iarr_hdifcoef3d_scal)
               outdat(ivar) = hdifscal_fac*hdifcoef3datc(i,j,kc)
            CASE (iarr_hstokesmag)
               CALL complex_polar(ustokesatc(i,j,kc),vstokesatc(i,j,kc),&
                                & xamp=outdat(ivar))
            CASE (iarr_hswdissipmag)
               CALL complex_polar(uswdissipatc(i,j,kc),vswdissipatc(i,j,kc),&
                                & xamp=outdat(ivar))
            CASE (iarr_hvelmag)
               CALL vector_mag_var_atc(uvel(i:i+1,j,kc),vvel(i,j:j+1,kc),&
                    & i,j,kc,1,1,vecmag=outdat(ivar))
            CASE (iarr_hvelmag_hadv_cfl)
               CALL vector_mag_var_atc(uvel(i:i+1,j,kc),vvel(i,j:j+1,kc),&
                                     & i,j,kc,1,1,vecmag=outdat(ivar))
               outdat(ivar) = ABS(outdat(ivar))*&
                               & (delt3d/MIN(delxatc(i,j),delyatc(i,j)))
            CASE (iarr_hveltotmag)
               CALL vector_mag_var_atc(uvel(i:i+1,j,kc),vvel(i,j:j+1,kc),&
                                     & i,j,kc,1,1,vecmag=outdat(ivar))
               CALL complex_polar(ustokesatc(i,j,kc),vstokesatc(i,j,kc),&
                                & xamp=varx)
               outdat(ivar) = outdat(ivar) + varx
            CASE (iarr_kinvisc)
               outdat(ivar) = Wvar_at_C(kinvisc(i,j,kc:kc+1),i,j)
            CASE (iarr_mom_vdif_cfl)
               IF (kc.EQ.1) THEN
                  outdat(ivar) = vdifcoefmom(i,j,1)*(delt3d/(delzatc(i,j,1)**2))
               ELSEIF (kc.EQ.nz) THEN
                  outdat(ivar) = vdifcoefmom(i,j,nz+1)*&
                               & (delt3d/(delzatc(i,j,nz)**2))
               ELSE
                  outdat(ivar) = Wvar_at_C(vdifcoefmom(i,j,kc:kc+1),i,j)*&
                                        & (delt3d/delzatc(i,j,kc)**2)
               ENDIF
            CASE (iarr_p3dbcgradatu)
               outdat(ivar) = Uvar_at_C(p3dbcgradatu(i:i+1,j,kc),i,j,kc,2,1)
            CASE (iarr_p3dbcgradatv)
               outdat(ivar) = Vvar_at_C(p3dbcgradatv(i,j:j+1,kc),i,j,kc,2,1)
            CASE (iarr_radiance)
               outdat(ivar) = Wvar_at_C(radiance(i,j,kc:kc+1),i,j)
            CASE (iarr_ricnum)
               IF (kc.EQ.1) THEN
                  outdat(ivar) = richardson_number_0d(i,j,2)
               ELSEIF (kc.EQ.nz) THEN
                  outdat(ivar) = richardson_number_0d(i,j,nz)
               ELSE
                  out1d(1) =  richardson_number_0d(i,j,k)
                  out1d(2) =  richardson_number_0d(i,j,k+1)
                  outdat(ivar) = Wvar_at_C(out1d,i,j)
               ENDIF
            CASE (iarr_sal); outdat(ivar) = sal(i,j,kc)
            CASE (iarr_scal_vdif_cfl)
               IF (kc.EQ.1) THEN
                  outdat(ivar) = vdifcoefscal(i,j,1)*&
                               & (delt3d/(delzatc(i,j,1)**2))
               ELSEIF (kc.EQ.nz) THEN
                  outdat(ivar) = vdifcoefscal(i,j,nz)*&
                              & (delt3d/(delzatc(i,j,nz)**2))
               ELSE
                  outdat(ivar) = Wvar_at_C(vdifcoefscal(i,j,kc:kc+1),i,j)*&
                                        & (delt3d/delzatc(i,j,kc)**2)
               ENDIF
            CASE (iarr_shearfreq2)
               IF (kc.EQ.1) THEN
                  outdat(ivar) = shearfreq2(i,j,2)
               ELSEIF (kc.EQ.nz) THEN
                  outdat(ivar) = shearfreq2(i,j,nz)
               ELSE
                  outdat(ivar) = Wvar_at_C(shearfreq2(i,j,kc:kc+1),i,j)
               ENDIF
            CASE (iarr_temp); outdat(ivar) = temp(i,j,kc)
            CASE (iarr_tke)
               IF (iopt_turb_tke_bbc.EQ.1.AND.kc.EQ.1) THEN
                  outdat(ivar) = tke(i,j,2)
               ELSEIF (iopt_turb_tke_sbc.EQ.1.AND.kc.EQ.nz) THEN
                  outdat(ivar) = tke(i,j,nz)
               ELSE
                  outdat(ivar) = Wvar_at_C(tke(i,j,kc:kc+1),i,j)
               ENDIF
            CASE (iarr_ubwdissipatc); outdat(ivar) = ubwdissipatc(i,j,kc)
            CASE (iarr_ustokesatc); outdat(ivar) = ustokesatc(i,j,kc)
            CASE (iarr_uswdissipatc); outdat(ivar) = uswdissipatc(i,j,kc)
            CASE (iarr_uvel)
               outdat(ivar) = Uvar_at_C(uvel(i:i+1,j,kc),i,j,kc,1,1)
            CASE (iarr_uvel_hadv_cfl)
               outdat(ivar) = Uvar_at_C(uvel(i:i+1,j,kc),i,j,kc,1,1)
               outdat(ivar) = ABS(outdat(ivar))*(delt3d/delxatc(i,j))
           CASE (iarr_uveltot)
               outdat(ivar) = Uvar_at_C(uvel(i:i+1,j,kc),i,j,kc,1,1)
               outdat(ivar) = outdat(ivar) + ustokesatc(i,j,kc)
            CASE (iarr_vbwdissipatc); outdat(ivar) = vbwdissipatc(i,j,kc)
            CASE (iarr_vdifcoefmom)
               outdat(ivar) = Wvar_at_C(vdifcoefmom(i,j,kc:kc+1),i,j)
            CASE (iarr_vdifcoefscal)
               outdat(ivar) = Wvar_at_C(vdifcoefscal(i,j,kc:kc+1),i,j)
            CASE (iarr_vdifcoefscal_norot)
               outdat(ivar) = Wvar_at_C(vdifcoefscal(i,j,kc:kc+1),i,j) - &
                            & Wvar_at_C(vdifcoefscal_rot(i,j,kc:kc+1),i,j)
            CASE (iarr_vdifcoefscal_rot)
               outdat(ivar) = Wvar_at_C(vdifcoefscal_rot(i,j,kc:kc+1),i,j)
            CASE (iarr_vdifcoeftke)
               outdat(ivar) = Wvar_at_C(vdifcoeftke(i,j,kc:kc+1),i,j)
            CASE (iarr_vortic3d); outdat(ivar) = vorticity_3d(i,j,kc)
            CASE (iarr_vstokesatc); outdat(ivar) = vstokesatc(i,j,kc)
            CASE (iarr_vswdissipatc); outdat(ivar) = vswdissipatc(i,j,kc)
            CASE (iarr_vvel)
               outdat(ivar) = Vvar_at_C(vvel(i,j:j+1,kc),i,j,kc,1,1)
            CASE (iarr_vvel_hadv_cfl)
               outdat(ivar) = Vvar_at_C(vvel(i,j:j+1,kc),i,j,kc,1,1)
               outdat(ivar) = ABS(outdat(ivar))*(delt3d/delyatc(i,j))
            CASE (iarr_vveltot)
               outdat(ivar) = Vvar_at_C(vvel(i,j:j+1,kc),i,j,kc,1,1)
               outdat(ivar) = outdat(ivar) + vstokesatc(i,j,kc)
            CASE (iarr_wphys); outdat(ivar) = wphys(i,j,kc)
            CASE (iarr_wstokesatw)
               outdat(ivar) = Wvar_at_C(wstokesatw(i,j,kc:kc+1),i,j)
            CASE (iarr_wvel)
               outdat(ivar) = Wvar_at_C(wvel(i,j,kc:kc+1),i,j)
            CASE (iarr_wvel_vadv_cfl)
               outdat(ivar) = Wvar_at_C(wvel(i,j,kc:kc+1),i,j)
               outdat(ivar) = ABS(outdat(ivar))*(delt3d/delzatc(i,j,kc))
            CASE (iarr_xslopeatu_geo)
               outdat(ivar) = Uvar_at_C(xslopeatu_geo(i:i+1,j,kc),i,j,kc,2,1)
            CASE (iarr_xslopeatu_ziso)
               out1d = 0.25*(xslopeatu_ziso(i-1:i,j,kc,1,0)+&
                           & xslopeatu_ziso(i:i+1,j,kc,0,0)+&
                           & xslopeatu_ziso(i-1:i,j,kc,1,1)+&
                           & xslopeatu_ziso(i:i+1,j,kc,0,1))
               outdat(ivar) = 0.5*(out1d(1)+out1d(2))
            CASE (iarr_xslopeatw_geo)
               IF (kc.EQ.1) THEN
                  outdat(ivar) = xslopeatw_geo(i,j,2)
               ELSEIF (k.EQ.nz) THEN
                  outdat(ivar) = xslopeatw_geo(i,j,nz)
               ELSE
                  outdat(ivar) = Wvar_at_C(xslopeatw_geo(i,j,kc:kc+1),i,j)
               ENDIF
            CASE (iarr_yslopeatv_geo)
               outdat(ivar) = Vvar_at_C(xslopeatu_geo(i,j:j+1,kc),i,j,kc,2,1)
            CASE (iarr_yslopeatv_ziso)
               out1d = 0.25*(yslopeatv_ziso(i,j-1:j,kc,1,0)+&
                           & yslopeatv_ziso(i,j:j+1,kc,0,0)+&
                           & yslopeatv_ziso(i,j-1:j,kc,1,1)+&
                           & yslopeatv_ziso(i,j:j+1,kc,0,1))
               outdat(ivar) = 0.5*(out1d(1)+out1d(2))
            CASE (iarr_yslopeatw_geo)
               IF (kc.EQ.1) THEN
                  outdat(ivar) = yslopeatw_geo(i,j,2)
               ELSEIF (k.EQ.nz) THEN
                  outdat(ivar) = yslopeatw_geo(i,j,nz)
               ELSE
                  outdat(ivar) = Wvar_at_C(yslopeatw_geo(i,j,kc:kc+1),i,j)
               ENDIF
            CASE (iarr_zlmix)
               outdat(ivar) = Wvar_at_C(zlmix(i,j,kc:kc+1),i,j)
!           ---sediment variables
            CASE (MaxModArids+1:MaxModArids+MaxSedArids)
               CALL define_sed3d_vals(outdat(ivar),i,j,kc,ivaridx,numvarx)
!           ---biological variables
            CASE (MaxModArids+MaxSedArids+1:MaxModArids+MaxSedArids+MaxBioArids)
               CALL define_bio3d_vals(outdat(ivar),i,j,kc,ivaridx)
!           ---particle model variables
            CASE (MaxModArids+MaxSedArids+MaxBioArids+1:MaxtotArids)
               CALL define_part3d_vals(outdat(ivar),i,j,kc,ivaridx,numvarx)
         END SELECT

!
!2.3.2 W-node level
!------------------
!

      ELSEIF (TRIM(nodex).EQ.'W') THEN 

         SELECT CASE(ivaridx)
  
            CASE (iarr_buofreq2)
               IF (kc.EQ.1) THEN
                  outdat(ivar) = buofreq2(i,j,2)
               ELSEIF (kc.EQ.nz+1) THEN
                  outdat(ivar) = buofreq2(i,j,nz)
               ELSE
                  outdat(ivar) = buofreq2(i,j,kc)
               ENDIF
            CASE (iarr_dissip)
               IF (kc.EQ.1) THEN
                  outdat(ivar) = dissip(i,j,2)
               ELSEIF (kc.EQ.nz+1) THEN
                  outdat(ivar) = dissip(i,j,nz)
               ELSE
                  outdat(ivar) = dissip(i,j,kc)
               ENDIF
            CASE (iarr_kinvisc)
               outdat(ivar) = kinvisc(i,j,kc)
            CASE (iarr_radiance)
               outdat(ivar) = radiance(i,j,kc)
            CASE (iarr_mom_vdif_cfl)
               IF (kc.EQ.1) THEN
                  outdat(ivar) = vdifcoefmom(i,j,1)*(delt3d/(delzatc(i,j,1)**2))
               ELSEIF (kc.EQ.nz) THEN
                  outdat(ivar) = vdifcoefmom(i,j,nz+1)*&
                               & (delt3d/(delzatc(i,j,nz)**2))
               ELSE
                  outdat(ivar) = vdifcoefmom(i,j,kc)*(delt3d/delzatw(i,j,kc)**2)
               ENDIF
            CASE (iarr_ricnum)
               IF (kc.EQ.1) THEN
                  outdat(ivar) = richardson_number_0d(i,j,2)
               ELSEIF (kc.EQ.nz+1) THEN
                  outdat(ivar) = richardson_number_0d(i,j,nz)
               ELSE
                  outdat(ivar) =  richardson_number_0d(i,j,kc)
               ENDIF
            CASE (iarr_scal_vdif_cfl)
               IF (kc.EQ.1) THEN
                  outdat(ivar) = vdifcoefscal(i,j,1)*&
                               & (delt3d/(delzatc(i,j,1)**2))
               ELSEIF (kc.EQ.nz) THEN
                  outdat(ivar) = vdifcoefscal(i,j,nz+1)*&
                               & (delt3d/(delzatc(i,j,nz)**2))
               ELSE
                  outdat(ivar) = vdifcoefscal(i,j,kc)*&
                               & (delt3d/delzatw(i,j,kc)**2)
               ENDIF
            CASE (iarr_shearfreq2)
               IF (kc.EQ.1) THEN
                  outdat(ivar) = shearfreq2(i,j,2)
               ELSEIF (kc.EQ.nz+1) THEN
                  outdat(ivar) = shearfreq2(i,j,nz)
               ELSE
                  outdat(ivar) = shearfreq2(i,j,kc)
               ENDIF
            CASE (iarr_tke)
               IF (iopt_turb_tke_bbc.EQ.1.AND.kc.EQ.1) THEN
                  outdat(ivar) = tke(i,j,2)
               ELSEIF (iopt_turb_tke_sbc.EQ.1.AND.kc.EQ.nz) THEN
                  outdat(ivar) = tke(i,j,nz)
               ELSE
                  outdat(ivar) = tke(i,j,kc)
               ENDIF
            CASE (iarr_vdifcoefmom)
               outdat(ivar) = vdifcoefmom(i,j,kc)
            CASE (iarr_vdifcoefscal)
               outdat(ivar) = vdifcoefscal(i,j,kc)
            CASE (iarr_vdifcoeftke)
               outdat(ivar) = vdifcoeftke(i,j,kc)
            CASE (iarr_wvel)
               outdat(ivar) = wvel(i,j,kc)
            CASE (iarr_wvel_vadv_cfl)
               IF (kc.EQ.1.OR.kc.EQ.nz+1) THEN
                  outdat(ivar) = 0.0
               ELSE
                  outdat(ivar) = ABS(wvel(i,j,kc))*(delt3d/delzatw(i,j,kc))
               ENDIF
            CASE (iarr_zlmix)
               outdat(ivar) = zlmix(i,j,kc)
         END SELECT

      ENDIF

!
!2.4 Other operators
!-------------------
!

   ELSE

      CALL inquire_var(ivaridx,node=cnode)

!
!2.4.1 C-node variable
!---------------------
!

      IF (TRIM(cnode).EQ.'C') THEN

!
!2.4.1.1 Store data
!------------------   
!

         SELECT CASE (ivaridx)
  
            CASE (iarr_beta_sal); out1dc = beta_sal(i,j,:)
            CASE (iarr_beta_temp); out1dc = beta_temp(i,j,:)
            CASE (iarr_dens); out1dc = dens(i,j,:) - 1000.0
            CASE (iarr_hbwdissipmag)
               CALL complex_polar(ubwdissipatc(i,j,:),vbwdissipatc(i,j,:),&
                                & xamp=out1dc)
            CASE (iarr_hdifcoef3d_mom)
               out1dc = hdifmom_fac*hdifcoef3datc(i,j,:)
            CASE (iarr_hdifcoef3d_scal)
               out1dc = hdifscal_fac*hdifcoef3datc(i,j,:)
            CASE (iarr_hstokesmag)
               CALL complex_polar(ustokesatc(i,j,:),vstokesatc(i,j,:),&
                                & xamp=out1dc)
            CASE (iarr_hswdissipmag)
               CALL complex_polar(uswdissipatc(i,j,:),vswdissipatc(i,j,:),&
                                & xamp=out1dc)
            CASE (iarr_hvelmag)
               k_2411: DO k=1,nz
                  CALL vector_mag_var_atc(uvel(i:i+1,j,k),vvel(i,j:j+1,k),&
                                        & i,j,k,1,1,vecmag=out1dc(k))
               ENDDO k_2411
            CASE (iarr_hvelmag_hadv_cfl)
               k_2412: DO k=1,nz
                  CALL vector_mag_var_atc(uvel(i:i+1,j,k),vvel(i,j:j+1,k),&
                       & i,j,k,1,1,vecmag=out1dc(k))
                  out1dc(k) = ABS(out1dc(k))*&
                               & (delt3d/MIN(delxatc(i,j),delyatc(i,j)))
               ENDDO k_2412

            CASE (iarr_hveltotmag)
               k_2413: DO k=1,nz
                  CALL vector_mag_var_atc(uvel(i:i+1,j,k),vvel(i,j:j+1,k),&
                                        & i,j,k,1,1,vecmag=out1dc(k))
                  CALL complex_polar(ustokesatc(i,j,k),vstokesatc(i,j,k),&
                                   & xamp=varx)
                  out1dc(k) = out1dc(k) + varx 
               ENDDO k_2413
            CASE (iarr_sal); out1dc = sal(i,j,:)
            CASE (iarr_temp); out1dc = temp(i,j,:)
            CASE (iarr_ubwdissipatc); out1dc = ubwdissipatc(i,j,:)
            CASE (iarr_ustokesatc); out1dc = ustokesatc(i,j,:)
            CASE (iarr_uswdissipatc); out1dc = uswdissipatc(i,j,:)
            CASE (iarr_uvel)
               CALL Uarr_at_C(uvel(i:i+1,j:j,:),out3dc,1,1,(/i,j,1/),&
                            & (/i+1,j,nz/),1,ivaridx,.FALSE.)
               out1dc = out3dc(1,1,:)
            CASE (iarr_uvel_hadv_cfl)
               CALL Uarr_at_C(uvel(i:i+1,j:j,:),out3dc,1,1,(/i,j,1/),&
                            & (/i+1,j,nz/),1,ivaridx,.FALSE.)
               out1dc = ABS(out3dc(1,1,:))*(delt3d/delxatc(i,j))
            CASE (iarr_uveltot)
               CALL Uarr_at_C(uvel(i:i+1,j:j,:),out3dc,1,1,(/i,j,1/),&
                           & (/i+1,j,nz/),1,ivaridx,.FALSE.)
               out1dc = out3dc(1,1,:) + ustokesatc(i,j,:)
            CASE (iarr_vbwdissipatc); out1dc = vbwdissipatc(i,j,:)
            CASE (iarr_vortic3d)
               k_2414: DO k=1,nz
                  out1dc(k) = vorticity_3d(i,j,k)
               END DO k_2414
            CASE (iarr_vstokesatc); out1dc = vstokesatc(i,j,:)
            CASE (iarr_vswdissipatc); out1dc = vswdissipatc(i,j,:)
            CASE (iarr_vvel)
               CALL Varr_at_C(vvel(i:i,j:j+1,:),out3dc,1,1,(/i,j,1/),&
                           & (/i,j+1,nz/),1,ivaridx,.FALSE.)
               out1dc = out3dc(1,1,:)
            CASE (iarr_vvel_hadv_cfl)
               CALL Varr_at_C(vvel(i:i,j:j+1,:),out3dc,1,1,(/i,j,1/),&
                           & (/i,j+1,nz/),1,ivaridx,.FALSE.)
               out1dc = ABS(out3dc(1,1,:))*(delt3d/delyatc(i,j))
            CASE (iarr_vveltot)
               CALL Uarr_at_C(vvel(i:i,j:j+1,:),out3dc,1,1,(/i,j,1/),&
                            & (/i+1,j,nz/),1,ivaridx,.FALSE.)
               out1dc = out3dc(1,1,:) + vstokesatc(i,j,:)
            CASE (iarr_wphys); out1dc = wphys(i,j,:)
         END SELECT

!
!2.4.1.2 Apply operator
!----------------------
!

         SELECT CASE (ooptx)
            CASE (oopt_min); outdat(ivar) = MINVAL(out1dc)
            CASE (oopt_max); outdat(ivar) = MAXVAL(out1dc)
            CASE (oopt_mean)
               outdat(ivar) = SUM(delzatc(i,j,:)*out1dc)/deptotatc(i,j)
            CASE (oopt_int)
               outdat(ivar) = SUM(delzatc(i,j,:)*out1dc)
            CASE (oopt_dep)
               CALL intpol1d_model_to_dep(out1dc,outdat(ivar),i,j,depx,'C  ')
         END SELECT

!
!2.4.2 W-node variable
!---------------------
!

      ELSEIF (TRIM(cnode).EQ.'W') THEN

!
!2.4.2.1 Store data
!------------------   
!

         SELECT CASE(ivaridx)
            CASE (iarr_buofreq2)
               nzmin = 2; nzmax= nz
               out1dw(2:nz) = buofreq2(i,j,2:nz)
            CASE (iarr_dissip)
               nzmin = 2; nzmax= nz
               out1dw(2:nz) = dissip(i,j,2:nz)
            CASE (iarr_kinvisc)
               nzmin = 1; nzmax = nz+1
               out1dw = kinvisc(i,j,:)
            CASE (iarr_mom_vdif_cfl)
               nzmin = 1; nzmax = nz+1
               out1dw(1) = vdifcoefmom(i,j,1)*(delt3d/(delzatc(i,j,1)**2))
               out1dw(nz+1) = vdifcoefmom(i,j,nz+1)*&
                            & (delt3d/(delzatc(i,j,nz)**2))
               out1dw(2:nz) = vdifcoefmom(i,j,2:nz)*&
                            & (delt3d/delzatw(i,j,2:nz)**2)
            CASE (iarr_radiance)
               nzmin = 1; nzmax = nz+1
               out1dw = radiance(i,j,:)
            CASE (iarr_ricnum)
               nzmin = 2; nzmax= nz
               k_2421: DO k=2,nz
                  out1dw(k) = richardson_number_0d(i,j,k)
               END DO k_2421
            CASE (iarr_scal_vdif_cfl)
               nzmin = 1; nzmax = nz+1
               out1dw(1) = vdifcoefscal(i,j,1)*(delt3d/(delzatc(i,j,1)**2))
               out1dw(nz+1) = vdifcoefscal(i,j,nz+1)*&
                            & (delt3d/(delzatc(i,j,nz)**2))
               out1dw(2:nz) = vdifcoefscal(i,j,2:nz)*&
                            & (delt3d/delzatw(i,j,2:nz)**2)
            CASE (iarr_shearfreq2)
               nzmin = 2; nzmax= nz
               out1dw(2:nz) = shearfreq2(i,j,2:nz)
            CASE (iarr_tke)
               nzmin = MERGE(2,1,iopt_turb_tke_bbc.EQ.1)
               nzmax = MERGE(nz,nz+1,iopt_turb_tke_sbc.EQ.1)
               out1dw(nzmin:nzmax) = tke(i,j,nzmin:nzmax)
            CASE (iarr_vdifcoefmom)
               nzmin = 1; nzmax = nz+1
               out1dw = vdifcoefmom(i,j,:)
            CASE (iarr_vdifcoefscal)
               nzmin = 1; nzmax = nz+1
               out1dw = vdifcoefscal(i,j,:)
            CASE (iarr_vdifcoeftke)
               nzmin = 1; nzmax = nz+1
               out1dw = vdifcoeftke(i,j,:)
            CASE (iarr_wvel)
               nzmin = 1; nzmax = nz+1
               out1dw = wvel(i,j,:)
            CASE (iarr_wvel_vadv_cfl)
               nzmin = 1; nzmax = nz+1
               out1dw(1) = 0.0; out1dw(nz+1) = 0.0
               out1dw(2:nz) = ABS(wvel(i,j,2:nz))*(delt3d/delzatw(i,j,2:nz))
            CASE (iarr_zlmix)
               nzmin = 1; nzmax = nz+1
               out1dw = zlmix(i,j,:)

      END SELECT

!
!2.4.2.2 Apply operator
!----------------------
!

         SELECT CASE (ooptx)
            CASE (oopt_min); outdat(ivar) = MINVAL(out1dw(nzmin:nzmax))
            CASE (oopt_max); outdat(ivar) = MAXVAL(out1dw(nzmin:nzmax))
            CASE (oopt_mean)
               outdat(ivar) = 0.5*delzatw(i,j,nzmin)*out1dw(nzmin) + &
                            & 0.5*delzatw(i,j,nzmax)*out1dw(nzmax)
               IF (nzmax-nzmin.GE.2) THEN
                  outdat(ivar) = outdat(ivar) + &
                               & SUM(delzatw(i,j,nzmin+1:nzmax-1)*&
                               & out1dw(nzmin+1:nzmax-1))/deptotatc(i,j)
               ENDIF
            CASE (oopt_int)
               outdat(ivar) = 0.5*(delzatw(i,j,nzmin)*out1dw(nzmin) + &
                            & delzatw(i,j,nzmax)*out1dw(nzmax))/deptotatc(i,j)
               IF (nzmax-nzmin.GE.2) THEN
                  outdat(ivar) = outdat(ivar) + &
                               & SUM(delzatw(i,j,nzmin+1:nzmax-1)*&
                               & out1dw(nzmin+1:nzmax-1))
               ENDIF
            CASE (oopt_dep)
               CALL intpol1d_model_to_dep(out1dw,outdat(ivar),i,j,depx,'W  ',&
                                        & nzmin=nzmin,nzmax=nzmax)
         END SELECT

      ENDIF

   ENDIF

ENDDO ivar_200


RETURN

END SUBROUTINE define_out2d_vals

!========================================================================

SUBROUTINE define_out3d_vals(outdat,i,j,k,novars,outvars,ivarid,numvar,node)
!************************************************************************
!
! *define_out3d_vals* Return 3-D output data values using the variable key ids
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)model_output.f90  V2.11.2
!
! Description - either outvars or ivarid need to be present but not both
!             - numvar may only be present if ivarid is present
!
! Module calls - complex_polar, define_bio3d_vals, define_sed3d_vals, 
!                define_part3d_vals, energy_flux_3d, energy_force_3d, energy_3d,
!                richardson_number_0d, Uvar_at_C, vector_mag_var_atc,
!                vorticity_3d, Vvar_at_C, Wvar_at_C
!
!************************************************************************
!
USE currents
USE datatypes
USE density
USE diffusion
USE grid
USE gridpars
USE modids
USE optics
USE physpars
USE switches
USE syspars
USE timepars
USE turbulence
USE wavevars
USE array_interp, ONLY: Uvar_at_C, Vvar_at_C, Wvar_at_C
USE biology_output, ONLY: define_bio3d_vals
USE diagnostic_routines, ONLY: energy_flux_3d, energy_force_3d, energy_3d, &
                             & vorticity_3d
USE math_library, ONLY: complex_polar, vector_mag_var_atc
USE particle_output, ONLY: define_part3d_vals
USE sediment_output, ONLY: define_sed3d_vals
USE turbulence_routines, ONLY: richardson_number_0d

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, k, novars
REAL, INTENT(OUT), DIMENSION(novars) :: outdat
CHARACTER (LEN=lennode), OPTIONAL, DIMENSION(novars) :: node
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(novars) :: ivarid, numvar
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(novars) :: outvars

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    Returned output values
!*i*       INTEGER (Local) X-index of the data location
!*j*       INTEGER (Local) Y-index of the data location
!*k*       INTEGER Vertical index of the data location
!*novars*  INTEGER Number of variables
!*outvars* DERIVED Variable's attributes (if varid is not present)
!*ivarid*  INTEGER Variable's key ids (if outvars is not present)
!*numvar*  INTEGER Variable number (sediment module only)
!*node*    CHAR    Vertical grid node of the output data:'C' (default) or 'W'
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: eflag3d, fflag3d
CHARACTER (LEN=lennode):: nodex
INTEGER :: ivar, ivaridx, numvarx
REAL :: varx, wsum
REAL, DIMENSION(2,2) :: weights
REAL, DIMENSION(2) :: out1d
REAL, DIMENSION(3) :: ecomps3d
REAL, DIMENSION(4) :: wcomps3d
REAL, DIMENSION(4,3) :: fcomps3d


IF (.NOT.maskatc_int(i,j)) THEN
   outdat = 0.0
   RETURN
ENDIF

eflag3d = .FALSE.; fflag3d = .FALSE.

ivar_1000: DO ivar=1,novars

   IF (PRESENT(outvars)) THEN
      ivaridx = outvars(ivar)%ivarid
      numvarx = outvars(ivar)%numvar
      nodex = outvars(ivar)%node
   ELSE
      ivaridx = ivarid(ivar)
      IF (PRESENT(node)) THEN
         nodex = node(ivar)
      ELSE
         nodex = 'C'
      ENDIF
      IF (PRESENT(numvar)) THEN
         numvarx = numvar(ivar)
      ELSEIF (PRESENT(outvars)) THEN
         numvarx = outvars(ivar)%numvar
      ELSE
         numvarx = 0
      ENDIF
   ENDIF

   IF (ivaridx.EQ.0) CYCLE ivar_1000 

!  ---energy calls
   SELECT CASE (ivaridx)
      CASE (iarr_edens3d,iarr_ekin3d,iarr_etot3d)
         IF (.NOT.eflag3d) THEN
            CALL energy_3d(i,j,k,ecomps3d)
            eflag3d = .TRUE.
         ENDIF
      CASE (iarr_eflux3du,iarr_eflux3dv,iarr_eflux3dw)
         IF (.NOT.fflag3d) THEN
            CALL energy_flux_3d(i,j,k,fcomps3d)
            fflag3d = .TRUE.
         ENDIF
   END SELECT

!  ---evaluate data (C-nodes)
   IF (TRIM(nodex).EQ.'C') THEN
      SELECT CASE (ivaridx)
         CASE (iarr_beta_sal); outdat(ivar) = beta_sal(i,j,k)
         CASE (iarr_beta_temp); outdat(ivar) = beta_temp(i,j,k)
         CASE (iarr_buofreq2)
            IF (k.EQ.1) THEN
               outdat(ivar) = buofreq2(i,j,2)
            ELSEIF (k.EQ.nz) THEN
               outdat(ivar) = buofreq2(i,j,nz)
            ELSE
               outdat(ivar) = Wvar_at_C(buofreq2(i,j,k:k+1),i,j)
            ENDIF
         CASE (iarr_dens); outdat(ivar) = dens(i,j,k) - 1000.0
         CASE (iarr_dissip)
            IF (k.EQ.1) THEN
               outdat(ivar) = dissip(i,j,2)
            ELSEIF (k.EQ.nz) THEN
               outdat(ivar) = dissip(i,j,nz)
            ELSE
               outdat(ivar) = Wvar_at_C(dissip(i,j,k:k+1),i,j)
            ENDIF
         CASE (iarr_edens3d)
            outdat(ivar) = ecomps3d(2)
         CASE (iarr_edissip3d)
            CALL energy_force_3d(i,j,k,wcomps3d)
            outdat(ivar) = -wcomps3d(4)
         CASE (iarr_eflux3du)
            outdat(ivar) = fcomps3d(4,1)
         CASE (iarr_eflux3dv)
            outdat(ivar) = fcomps3d(4,2)
         CASE (iarr_eflux3dw)
            outdat(ivar) = fcomps3d(4,3)
         CASE (iarr_ekin3d)
            outdat(ivar) = ecomps3d(1)
         CASE (iarr_etot3d)
            outdat(ivar) = ecomps3d(3)
         CASE (iarr_hdifcoef3d_mom)
            outdat(ivar)= hdifmom_fac*hdifcoef3datc(i,j,k)
         CASE (iarr_hdifcoef3d_scal)
            outdat(ivar)= hdifscal_fac*hdifcoef3datc(i,j,k)
         CASE (iarr_hbwdissipmag)
            CALL complex_polar(ubwdissipatc(i,j,k),vbwdissipatc(i,j,k),&
                             & xamp=outdat(ivar))
         CASE (iarr_hstokesmag)
            CALL complex_polar(ustokesatc(i,j,k),vstokesatc(i,j,k),&
                             & xamp=outdat(ivar))
         CASE (iarr_hswdissipmag)
            CALL complex_polar(uswdissipatc(i,j,k),vswdissipatc(i,j,k),&
                             & xamp=outdat(ivar))
         CASE (iarr_hvelmag)
            CALL vector_mag_var_atc(uvel(i:i+1,j,k),vvel(i,j:j+1,k),&
                 & i,j,k,1,1,vecmag=outdat(ivar))
         CASE (iarr_hvelmag_hadv_cfl)
            CALL vector_mag_var_atc(uvel(i:i+1,j,k),vvel(i,j:j+1,k),&
                                   & i,j,k,1,1,vecmag=outdat(ivar))
            outdat(ivar) = ABS(outdat(ivar))*&
                            & (delt3d/MIN(delxatc(i,j),delyatc(i,j)))
         CASE (iarr_hveltotmag)
            CALL vector_mag_var_atc(uvel(i:i+1,j,k),vvel(i,j:j+1,k),&
                                  & i,j,k,1,1,vecmag=outdat(ivar))
            CALL complex_polar(ustokesatc(i,j,k),vstokesatc(i,j,k),&
                             & xamp=varx)
            outdat(ivar) = outdat(ivar) + varx
         CASE (iarr_kinvisc)
            outdat(ivar) = kinvisc(i,j,k)
         CASE (iarr_mom_vdif_cfl)
            outdat(ivar) = Wvar_at_C(vdifcoefmom(i,j,k:k+1),i,j)*&
                 & (delt3d/delzatc(i,j,k)**2)
         CASE (iarr_p3dbcgradatu)
            outdat(ivar) = Uvar_at_C(p3dbcgradatu(i:i+1,j,k),i,j,k,2,1)
         CASE (iarr_p3dbcgradatv)
            outdat(ivar) = Vvar_at_C(p3dbcgradatv(i,j:j+1,k),i,j,k,2,1)
         CASE (iarr_radiance);
            outdat(ivar) = Wvar_at_C(radiance(i,j,k+1),i,j)
         CASE (iarr_ricnum)
            IF (k.EQ.1) THEN
               outdat(ivar) = richardson_number_0d(i,j,2)
            ELSEIF (k.EQ.nz) THEN
               outdat(ivar) = richardson_number_0d(i,j,nz)
            ELSE
               out1d(1) =  richardson_number_0d(i,j,k)
               out1d(2) =  richardson_number_0d(i,j,k+1)
               outdat(ivar) = Wvar_at_C(out1d,i,j)
            ENDIF
         CASE (iarr_sal); outdat(ivar) = sal(i,j,k)
         CASE (iarr_scal_vdif_cfl)
            outdat(ivar) = Wvar_at_C(vdifcoefscal(i,j,k:k+1),i,j)*&
                         & (delt3d/delzatc(i,j,k)**2)
         CASE (iarr_shearfreq2)
            IF (k.EQ.1) THEN
               outdat(ivar) = shearfreq2(i,j,2)
            ELSEIF (k.EQ.nz) THEN
               outdat(ivar) = shearfreq2(i,j,nz)
            ELSE
               outdat(ivar) = Wvar_at_C(shearfreq2(i,j,k:k+1),i,j)
            ENDIF
         CASE (iarr_temp); outdat(ivar) = temp(i,j,k)
         CASE (iarr_tke)
            IF (iopt_turb_tke_bbc.EQ.1.AND.k.EQ.1) THEN
               outdat(ivar) = tke(i,j,2)
            ELSEIF (iopt_turb_tke_sbc.EQ.1.AND.k.EQ.nz) THEN
               outdat(ivar) = tke(i,j,nz)
            ELSE
               outdat(ivar) = Wvar_at_C(tke(i,j,k:k+1),i,j)
            ENDIF
         CASE (iarr_ubwdissipatc); outdat(ivar) = ubwdissipatc(i,j,k)
         CASE (iarr_ustokesatc); outdat(ivar) = ustokesatc(i,j,k)
         CASE (iarr_uswdissipatc); outdat(ivar) = uswdissipatc(i,j,k)
         CASE (iarr_uvel); outdat(ivar) = Uvar_at_C(uvel(i:i+1,j,k),i,j,k,1,1)
         CASE (iarr_uvel_hadv_cfl)
            outdat(ivar) = Uvar_at_C(uvel(i:i+1,j,k),i,j,k,1,1)
            outdat(ivar) = ABS(outdat(ivar))*(delt3d/delxatc(i,j))
         CASE (iarr_uveltot)
            outdat(ivar) = Uvar_at_C(uvel(i:i+1,j,k),i,j,k,1,1)
            outdat(ivar) = outdat(ivar) + ustokesatc(i,j,k)
         CASE (iarr_vbwdissipatc); outdat(ivar) = vbwdissipatc(i,j,k)
         CASE (iarr_vdifcoefmom)
            outdat(ivar) = Wvar_at_C(vdifcoefmom(i,j,k:k+1),i,j)
         CASE (iarr_vdifcoefscal)
            outdat(ivar) = Wvar_at_C(vdifcoefscal(i,j,k:k+1),i,j)
         CASE (iarr_vdifcoefscal_norot)
            outdat(ivar) = Wvar_at_C(vdifcoefscal(i,j,k:k+1),i,j) - &
                        &  Wvar_at_C(vdifcoefscal_rot(i,j,k:k+1),i,j)
         CASE (iarr_vdifcoefscal_rot)
            outdat(ivar) = Wvar_at_C(vdifcoefscal_rot(i,j,k:k+1),i,j)
         CASE (iarr_vdifcoeftke)
            outdat(ivar) = Wvar_at_C(vdifcoeftke(i,j,k:k+1),i,j)
         CASE (iarr_vortic3d); outdat(ivar) = vorticity_3d(i,j,k)
         CASE (iarr_vstokesatc); outdat(ivar) = vstokesatc(i,j,k)
         CASE (iarr_vswdissipatc); outdat(ivar) = vswdissipatc(i,j,k)
         CASE (iarr_vvel); outdat(ivar) = Vvar_at_C(vvel(i,j:j+1,k),i,j,k,1,1)
         CASE (iarr_vvel_hadv_cfl)
            outdat(ivar) = Vvar_at_C(vvel(i,j:j+1,k),i,j,k,1,1)
            outdat(ivar) = ABS(outdat(ivar))*(delt3d/delyatc(i,j))
         CASE (iarr_vveltot)
            outdat(ivar) = Vvar_at_C(vvel(i,j:j+1,k),i,j,k,1,1)
            outdat(ivar) = outdat(ivar) + vstokesatc(i,j,k)
         CASE (iarr_wphys); outdat(ivar) = wphys(i,j,k)
         CASE (iarr_wstokesatw)
            outdat(ivar) = Wvar_at_C(wstokesatw(i,j,k:k+1),i,j)
         CASE (iarr_wvel)
            outdat(ivar) = Wvar_at_C(wvel(i,j,k:k+1),i,j)
         CASE (iarr_wvel_vadv_cfl)
            outdat(ivar) = Wvar_at_C(wvel(i,j,k:k+1),i,j)
            outdat(ivar) = ABS(outdat(ivar))*(delt3d/delzatc(i,j,k))
         CASE (iarr_xslopeatu_geo)
            outdat(ivar) = Uvar_at_C(xslopeatu_geo(i:i+1,j,k),i,j,k,1,1)
         CASE (iarr_xslopeatu_ziso)
            weights = MERGE(1.0,0.0,xslopeatu_ziso(i,j,k,1:0:-1,:).NE.0.0)
            wsum = SUM(weights)
            IF (wsum.GT.0.0) THEN
               outdat(ivar) = SUM(xslopeatu_ziso(i,j,k,:,:)*weights)/wsum
            ENDIF
         CASE (iarr_xslopeatw_geo)
            IF (k.EQ.1) THEN
               outdat(ivar)= xslopeatw_geo(i,j,2)
            ELSEIF (k.EQ.nz) THEN
               outdat(ivar)= xslopeatw_geo(i,j,nz)
            ELSE
               outdat(ivar) = Wvar_at_C(xslopeatu_geo(i,j,k:k+1),i,j)
            ENDIF
         CASE (iarr_yslopeatv_geo)
            outdat(ivar) = Vvar_at_C(yslopeatv_geo(i,j:j+1,k),i,j,k,1,1)
         CASE (iarr_yslopeatv_ziso)
            weights = MERGE(1.0,0.0,yslopeatv_ziso(i,j,k,1:0:-1,:).NE.0.0)
            wsum = SUM(weights)
            IF (wsum.GT.0.0) THEN
               outdat(ivar) = SUM(yslopeatv_ziso(i,j,k,:,:)*weights)/wsum
            ENDIF
         CASE (iarr_yslopeatw_geo)
            IF (k.EQ.1) THEN
               outdat(ivar)= yslopeatw_geo(i,j,2)
            ELSEIF (k.EQ.nz) THEN
               outdat(ivar)= yslopeatw_geo(i,j,nz)
            ELSE
               outdat(ivar) = Wvar_at_C(yslopeatw_geo(i,j,k:k+1),i,j)
            ENDIF
         CASE (iarr_zlmix)
            outdat(ivar) = Wvar_at_C(zlmix(i,j,k:k+1),i,j)

!        ---sediment variables
         CASE (MaxModArids+1:MaxModArids+MaxSedArids)
            CALL define_sed3d_vals(outdat(ivar),i,j,k,ivaridx,numvarx)
!        ---biological variables
         CASE (MaxModArids+MaxSedArids+1:MaxModArids+MaxSedArids+MaxBioArids)
            CALL define_bio3d_vals(outdat(ivar),i,j,k,ivaridx)
!        ---particle model variables
         CASE (MaxModArids+MaxSedArids+MaxBioArids+1:MaxTotArids)
            CALL define_part3d_vals(outdat(ivar),i,j,k,ivaridx,numvarx)

      END SELECT

!  ---evaluate data (W-nodes)
   ELSEIF (TRIM(nodex).EQ.'W') THEN
      SELECT CASE (ivaridx)
         CASE (iarr_buofreq2)
            IF (k.EQ.1) THEN
               outdat(ivar) = buofreq2(i,j,2)
            ELSEIF (k.EQ.nz+1) THEN
               outdat(ivar) = buofreq2(i,j,nz)
            ELSE
               outdat(ivar) = buofreq2(i,j,k)
            ENDIF
         CASE (iarr_dissip)
            IF (k.EQ.1) THEN
               outdat(ivar) = dissip(i,j,2)
            ELSEIF (k.EQ.nz+1) THEN
               outdat(ivar) = dissip(i,j,nz)
            ELSE
               outdat(ivar) = dissip(i,j,k)
            ENDIF
         CASE (iarr_kinvisc)
            outdat(ivar) = kinvisc(i,j,k)
         CASE (iarr_mom_vdif_cfl)
            IF (k.EQ.1) THEN
               outdat(ivar) = vdifcoefmom(i,j,1)*(delt3d/delzatc(i,j,2)**2) 
            ELSEIF (k.EQ.nz+1) THEN
               outdat(ivar) = vdifcoefmom(i,j,nz+1)*(delt3d/delzatc(i,j,nz)**2) 
            ELSE
               outdat(ivar) = vdifcoefmom(i,j,k)*(delt3d/delzatw(i,j,k)**2)
            ENDIF
         CASE (iarr_radiance); outdat(ivar) = radiance(i,j,k)
         CASE (iarr_ricnum)
            IF (k.EQ.1) THEN
               outdat(ivar) = richardson_number_0d(i,j,2)
            ELSEIF (k.EQ.nz+1) THEN
               outdat(ivar) = richardson_number_0d(i,j,nz)
            ELSE
               outdat(ivar) =  richardson_number_0d(i,j,k)
            ENDIF
         CASE (iarr_scal_vdif_cfl)
            IF (k.EQ.1) THEN
               outdat(ivar) = vdifcoefscal(i,j,1)*(delt3d/delzatc(i,j,2)**2) 
            ELSEIF (k.EQ.nz+1) THEN
               outdat(ivar) = vdifcoefscal(i,j,nz+1)*&
                            & (delt3d/delzatc(i,j,nz)**2) 
            ELSE
               outdat(ivar) = vdifcoefscal(i,j,k)*(delt3d/delzatw(i,j,k)**2)
            ENDIF
         CASE (iarr_shearfreq2)
            IF (k.EQ.1) THEN
               outdat(ivar) = shearfreq2(i,j,2)
            ELSEIF (k.EQ.nz+1) THEN
               outdat(ivar) = shearfreq2(i,j,nz)
            ELSE
               outdat(ivar) = shearfreq2(i,j,k)
            ENDIF
         CASE (iarr_tke)
            IF (iopt_turb_tke_bbc.EQ.1.AND.k.EQ.1) THEN
               outdat(ivar) = tke(i,j,2)
            ELSEIF (iopt_turb_tke_sbc.EQ.1.AND.k.EQ.nz+1) THEN
               outdat(ivar) = tke(i,j,nz)
            ELSE
               outdat(ivar) = tke(i,j,k)
            ENDIF
         CASE (iarr_vdifcoefmom)
            outdat(ivar) = vdifcoefmom(i,j,k)
         CASE (iarr_vdifcoefscal)
            outdat(ivar) = vdifcoefscal(i,j,k)
         CASE (iarr_vdifcoeftke); outdat(ivar) = vdifcoeftke(i,j,k)
         CASE (iarr_wvel); outdat(ivar) = wvel(i,j,k)
         CASE (iarr_wvel_vadv_cfl)
            IF (k.EQ.1.OR.k.EQ.nz+1) THEN
               outdat(ivar) = 0.0
            ELSE
               outdat(ivar) = ABS(wvel(i,j,k))*(delt3d/delzatw(i,j,k))
            ENDIF
         CASE (iarr_zlmix); outdat(ivar) = zlmix(i,j,k)

      END SELECT
   ENDIF

ENDDO ivar_1000


RETURN

END SUBROUTINE define_out3d_vals


END MODULE model_output
