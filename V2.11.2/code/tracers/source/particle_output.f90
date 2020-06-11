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

MODULE particle_output
!************************************************************************
!
! *particle_output* Routines defining particle output data
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_output.f90  V2.11
!
! $Date: 2017-12-22 12:37:21 +0100 (Fri, 22 Dec 2017) $
!
! $Revision: 1070 $
!
! Description -
!
! Reference -
!
! Routines - define_part0d_vals, define_part2d_vals, define_part3d_vals
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS

!========================================================================

SUBROUTINE define_part0d_vals(outdat,ivarid,l,oopt)
!************************************************************************
!
! *define_part0d_vals* Define 0-D particle model output data
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_output.f90  V2.11
!
! Description -
!
! Calling program - define_out0d_vals
!
! Module calls - error_alloc, inquire_var, max_vars, min_vars, sum2_vars,
!                Warr_at_C
!  
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE partids
USE partvars
USE syspars
USE array_interp, ONLY: Warr_at_C
USE error_routines, ONLY: error_alloc
USE modvars_routines, ONLY: inquire_var
USE paral_utilities, ONLY: max_vars, min_vars, sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
INTEGER, INTENT(IN) :: ivarid, l, oopt
REAL, INTENT(OUT) :: outdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    0-D output data value
!*ivarid*  INTEGER Variable key id
!*l*       INTEGER Contaminant label
!*oopt*    INTEGER Type of operator
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, nodim, nzdim
INTEGER, DIMENSION(4) :: nhdims = 0
INTEGER, DIMENSION(5) :: ngdims
REAL :: atot, vtot
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: area, out2d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: out3d, volume


procname(pglev+1) = 'define_part0d_vals'
CALL log_timer_in()


!
!1. Allocate arrays
!------------------
!

CALL inquire_var(ivarid,nrank=nodim,global_dims=ngdims)

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
!2. Area and volume of the global domain without land
!----------------------------------------------------
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
      k_210: DO k=1,nzdim
         volume(:,:,k) = area*delzatc(1:ncloc,1:nrloc,k)
      ENDDO k_210
      IF (oopt.EQ.oopt_mean) THEN
         CALL sum2_vars(volume,vtot,nhdims,'C  ',0,commall=.TRUE.)
      ENDIF
   ENDIF
ENDIF

!
!3. Store data
!-------------
!

SELECT CASE (ivarid)
   CASE (iarr_part_ageconc); out3d = ageconc(:,:,:,l)
   CASE (iarr_part_ageconc_tot); out3d = ageconc_tot
   CASE (iarr_part_conc); out3d = conc(:,:,:,l)
   CASE (iarr_part_conc_tot); out3d = conc_tot
   CASE (iarr_vdifcoefpart)
      CALL Warr_at_C(vdifcoefpart,out3d,(/1,1,1/),(/nc,nr,nz+1/),1,ivarid,&
                  & .FALSE.)
END SELECT

!
!4. Apply operator
!------------------
!

IF (nodim.EQ.2) THEN

   SELECT CASE (oopt)

      CASE (oopt_min)
         CALL min_vars(out2d,outdat,ivarid,commall=.TRUE.,&
                     & mask=maskatc(1:nc,1:nr))
      CASE (oopt_max)
         CALL max_vars(out2d,outdat,ivarid,commall=.TRUE.,&
                     & mask=maskatc(1:nc,1:nr))
      CASE (oopt_mean)
         CALL sum2_vars(area*out2d,outdat,nhdims,'C  ',ivarid,commall=.TRUE.)
         outdat = outdat/atot
      CASE (oopt_int)
         CALL sum2_vars(area*out2d,outdat,nhdims,'C  ',ivarid,commall=.TRUE.)
         
   END SELECT

ELSEIF (nodim.EQ.3) THEN

   SELECT CASE (oopt)

      CASE (oopt_min)
         CALL min_vars(out3d,outdat,ivarid,commall=.TRUE.,&
                     & mask=maskatc(1:nc,1:nr))
      CASE (oopt_max)
         CALL max_vars(out3d,outdat,ivarid,commall=.TRUE.,&
                     & mask=maskatc(1:nc,1:nr))
      CASE (oopt_mean)
         CALL sum2_vars(volume*out3d,outdat,nhdims,'C  ',ivarid,commall=.TRUE.)
         outdat = outdat/vtot
      CASE (oopt_int)
         CALL sum2_vars(volume*out3d,outdat,nhdims,'C  ',ivarid,commall=.TRUE.)

   END SELECT

ENDIF

!
!5. Deallocate arrays
!--------------------
!

IF (ALLOCATED(out2d)) DEALLOCATE (out2d)
IF (ALLOCATED(out3d)) DEALLOCATE (out3d)
IF (ALLOCATED(area)) DEALLOCATE (area)
IF (ALLOCATED(volume)) DEALLOCATE (volume)

CALL log_timer_out()


RETURN

END SUBROUTINE define_part0d_vals

!========================================================================

SUBROUTINE define_part2d_vals(outdat,i,j,ivarid,l,oopt,kc,dep)
!************************************************************************
!
! *define_part2d_vals*  Define particle model output data at a 2-D grid location
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_output.f90  V2.11
!
! Description -
!
! Calling program - define_out2d_vals
!
! Module calls - intpol1d_model_to_dep, Warr_at_C, Wvar_at_C
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE partids
USE partvars
USE array_interp, ONLY: Warr_at_C, Wvar_at_C
USE grid_interp, ONLY: intpol1d_model_to_dep

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, ivarid, j, kc, l, oopt
REAL, INTENT(IN) :: dep
REAL, INTENT(OUT) :: outdat
!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    2-D output data value
!*i*       INTEGER (Local) X-index of the data location
!*j*       INTEGER (Local) Y-index of the data location
!*ivarid*  INTEGER Variable key id
!*f*       INTEGER Contaminant label
!*oopt*    INTEGER Type of operator
!*kc*      INTEGER Vertical level at which a 3-D variable is evaluated
!*dep*     REAL    Depth below the surface at which a 3-D variable is evaluated
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(nz) :: out1d
REAL, DIMENSION(1,1,nz) :: out3d


!
!1. 3-D variable at a given level
!--------------------------------
!

IF (oopt.EQ.oopt_klev) THEN

   SELECT CASE (ivarid)
      CASE (iarr_part_ageconc); outdat = ageconc(i,j,kc,l)
      CASE (iarr_part_ageconc_tot); outdat = ageconc_tot(i,j,kc)
      CASE (iarr_part_conc); outdat = conc(i,j,kc,l)
      CASE (iarr_part_conc_tot); outdat = conc_tot(i,j,kc)
      CASE (iarr_vdifcoefpart)
         outdat = Wvar_at_C(vdifcoefpart(i,j,kc:kc+1),i,j)
   END SELECT

!
!2. Other operators
!------------------
!

ELSE

!
!2.1 Store data
!--------------
!

   SELECT CASE (ivarid)

      CASE (iarr_part_ageconc)
         out1d(1:nz) = ageconc(i,j,:,l)
      CASE (iarr_part_ageconc_tot)
         out1d(1:nz) = ageconc_tot(i,j,:)
      CASE (iarr_part_conc)
         out1d = conc(i,j,:,l)
      CASE (iarr_part_conc_tot)
         out1d = conc_tot(i,j,:)
      CASE (iarr_vdifcoefpart)
         CALL Warr_at_C(vdifcoefpart(i:i,j:j,:),out3d,(/i,j,1/),&
                     & (/i,j,nz+1/),1,ivarid,.FALSE.)
         out1d = out3d(1,1,:)
   END SELECT

!
!2.2 Apply operator
!------------------
!

   SELECT CASE (oopt)
      CASE (oopt_min); outdat = MINVAL(out1d)
      CASE (oopt_max); outdat = MAXVAL(out1d)
      CASE (oopt_mean)
         outdat = SUM(delzatc(i,j,:)*out1d)/deptotatc(i,j)
      CASE (oopt_int)
        outdat = SUM(delzatc(i,j,:)*out1d)
      CASE (oopt_dep)
         CALL intpol1d_model_to_dep(out1d,outdat,i,j,dep,'C  ')
   END SELECT

ENDIF


RETURN

END SUBROUTINE define_part2d_vals

!========================================================================

SUBROUTINE define_part3d_vals(outdat,i,j,k,ivarid,l)
!************************************************************************
!
! *define_part3d_vals*  Define particle model output data at a 3-D grid location
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)particle_output.f90  V2.11
!
! Description -
!
! Calling program - define_out2d_vals, define_out3d_vals
!
! Module calls - Wvar_at_C
!
!************************************************************************
!
USE partids
USE partvars  
USE array_interp, ONLY: Wvar_at_C

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, ivarid, j, k, l
REAL, INTENT(OUT) :: outdat
!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    3-D output data value
!*i*       INTEGER (Local) X-index of the data location
!*j*       INTEGER (Local) Y-index of the data location
!*k*       INTEGER Vertical index of the data location
!*ivarid*  INTEGER Variable key id
!*l*       INTEGER Contaminant label
!
!------------------------------------------------------------------------------
!


SELECT CASE (ivarid)

   CASE (iarr_part_ageconc); outdat = ageconc(i,j,k,l)
   CASE (iarr_part_ageconc_tot); outdat = ageconc_tot(i,j,k)
   CASE (iarr_part_conc); outdat = conc(i,j,k,l)
   CASE (iarr_part_conc_tot); outdat = conc_tot(i,j,k)
   CASE (iarr_vdifcoefpart); outdat = Wvar_at_C(vdifcoefpart(i,j,k:k+1),i,j)

END SELECT


RETURN

END SUBROUTINE define_part3d_vals


END MODULE particle_output
