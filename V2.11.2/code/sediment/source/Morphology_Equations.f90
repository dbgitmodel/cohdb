!
! Copyright 2013 RBINS-MUMM
!
! Licensed under the EUPL, Version 1.1 or - as soon they will be approved by
! the European Commission - subsequent versions of the EUPL (the "Licence");
! You may not use this work except in compliance with the Licence.
! You may obtain a copy fof the Licence at:
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
! *Morphology_Equations* morphological model
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.11.2
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description - This module contains the equations for calculating the 
!               morphological changes associated with different modes of
!               sediment transport, including parametrisations of processes
!               such as vertical sediment sorting and avalanching
!
! Routines - adjust_layer_top, adjust_one_layer, adjust_multi_layers,
!            avalanching, bed_art_diff, bed_elevation_equation, bed_time_int,
!            calculate_active_layer_thickness, exchange_fluxes_ribberink,
!            mass_balance, mass_correct_update, morphology_accel,
!            morphology_equation, morphological_step,
!            update_hydrodynamic_estimation
!
!************************************************************************
!

!========================================================================

SUBROUTINE adjust_layer_top
!************************************************************************
!
! *adjust_layer_top* adjust the thickness of the bed layers to the bed update 
!                    or the calculated active layer thickness
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.10
!
! Description -
!
! Reference -
!
! Calling program - adjust_multi_layers
!
! External calls - calculate_active_layer_thickness, median_particle_diameter
!
! Module calls -
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphswitches
USE sedarrays
USE sedpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: f, i, j, k, layer_number
REAL, DIMENSION(nb) :: bed_layer_depth, bed_layer_thick_old
REAL, DIMENSION(nb,nf) :: bed_frac_old


procname(pglev+1) = 'adjust_layer_top'
CALL log_timer_in()

!
!1. Median particle diameter
!---------------------------
!

CALL median_particle_diameter(d50_bed,50.0)

!
!2. Active layer thickness
!-------------------------
!

CALL calculate_active_layer_thickness

!
!3. Adjust top layer to active layer_thickness
!---------------------------------------------
!

j_300: DO j=1,nrloc
i_300: DO i=1,ncloc

   IF (maskatc_int(i,j).AND.bed_layer_thickness(i,j,1).GT.0.0.AND.&
     & active_layer_thickness(i,j).GT.0.0) THEN

!
!3.1. Initialise
!---------------
!
!     ---depth of the bed layers with respect to the bed surface
      bed_layer_depth(1) = bed_layer_thickness(i,j,1)
      k_311: DO k=2,nb
         bed_layer_depth(k) = bed_layer_depth(k-1) + bed_layer_thickness(i,j,k)
      ENDDO k_311

!     ---old values
      bed_layer_thick_old = bed_layer_thickness(i,j,:)
      bed_frac_old = bed_fraction(i,j,:,:)

!
!3.2 Number of layers making up the new top layer
!------------------------------------------------
!
!     ---number of layers
      IF (bed_layer_depth(nb).LE.active_layer_thickness(i,j)) THEN
         layer_number = nb
      ELSE
         k = 1; layer_number = 0
         DO WHILE (bed_layer_depth(k).LE.active_layer_thickness(i,j))
            layer_number = layer_number + 1
            k = k + 1
         ENDDO
      ENDIF

!
!3.2 Adjust bed layer thickness
!------------------------------
!
!     ---first layer
      bed_layer_thickness(i,j,1) = active_layer_thickness(i,j)

!     ---second layer
      IF (layer_number.LE.nb-1) THEN
          bed_layer_thickness(i,j,2) = bed_layer_depth(1+layer_number) - &
                                     & active_layer_thickness(i,j)
      ENDIF

!     ---other layers
      IF (layer_number.EQ.nb) THEN
         k_321: DO k=2,nb
            bed_layer_thickness(i,j,k) = 0.0
         ENDDO k_321

      ELSEIF (layer_number.GT.1) THEN
!        --assigning correct index to shifted layers
         k_322: DO k=3,nb-layer_number+1
            bed_layer_thickness(i,j,k) = bed_layer_thick_old(k-1+layer_number)
         ENDDO k_322
!        --creating dummy layers
         k_323: DO k=nb-layer_number+2,nb
            bed_layer_thickness(i,j,k) = 0.0
         ENDDO k_323

      ELSEIF (layer_number.EQ.0) THEN
!        --assigning correct index to shifted layers
         k_324: DO k=3,nb-1
            bed_layer_thickness(i,j,k) = bed_layer_thick_old(k-1)
         ENDDO k_324
!        ---merging bottom layers
         bed_layer_thickness(i,j,nb)= bed_layer_thick_old(nb-1) + &
                                    & bed_layer_thick_old(nb)

      ENDIF
      
!
!3.3 Adjust bed layer fraction
!-----------------------------
!

      f_330: DO f=1,nf

!        ---first layer
         IF (layer_number.EQ.nb) THEN
            bed_fraction(i,j,1,f) = sed_avail(i,j,f)/sed_avail_tot(i,j)
         ELSEIF (layer_number.GT.0) THEN
            bed_fraction(i,j,1,f) = &
               & (SUM(bed_layer_thick_old(1:layer_number)*&
               &      bed_frac_old(1:layer_number,f))+ &
               & (active_layer_thickness(i,j)-bed_layer_depth(layer_number))*&
               & bed_frac_old(layer_number+1,f))/bed_layer_thickness(i,j,1)
         ENDIF

!        ---other layers
         IF (layer_number.EQ.nb) THEN
            k_331: DO k=2,nb
               bed_fraction(i,j,k,:) = 0.0
            ENDDO k_331

         ELSEIF (layer_number.GT.1) THEN
!           --assigning correct index to shifted layers
            k_332: DO k=2,nb-layer_number+1
               bed_fraction(i,j,k,:) = bed_frac_old(k-1+layer_number,:)
            ENDDO k_332
!           --creating dummy layers
            k_333: DO k=nb-layer_number+2,nb
               bed_fraction(i,j,k,:) = 0.0
            ENDDO k_333

         ELSEIF (layer_number.EQ.0) THEN
!           --assigning correct index to shifted layers
            k_334: DO k=2,nb-1
               bed_fraction(i,j,k,:) = bed_frac_old(k-1,:)
            ENDDO k_334
!           ---merging bottom layers
            IF (bed_layer_thick_old(nb).GT.0.0.AND.&
              & bed_layer_thick_old(nb-1).GT.0.0) THEN
               bed_fraction(i,j,nb,:) = &
                           & (bed_frac_old(nb,:)*bed_layer_thick_old(nb)+&
                           &  bed_frac_old(nb-1,:)*bed_layer_thick_old(nb-1))/&
                           &  bed_layer_thickness(i,j,nb)
            ELSEIF (bed_layer_thick_old(nb).EQ.0.0.AND.&
                    bed_layer_thick_old(nb-1).GT.0.0) THEN
               bed_fraction(i,j,nb,:) = bed_frac_old(nb-1,:)
            ELSE
               bed_fraction(i,j,nb,:) = 0.0
            ENDIF

         ENDIF
         
      ENDDO f_330

   ENDIF

ENDDO i_300
ENDDO j_300

!
!4. Median particle diameter
!---------------------------
!

CALL median_particle_diameter(d50_bed,50.0)

CALL log_timer_out()


RETURN

END SUBROUTINE adjust_layer_top

!========================================================================

SUBROUTINE adjust_one_layer
!************************************************************************
!
! *adjust_bed_layer* adjust bed layer thickness and fractions in case of
!                    a one-layer bed model
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - morphology_equation
!
! External calls - mass_correct_update, median_particle_diameter
!
! Module calls - collect_log_expr, error_alloc
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphswitches
USE sedarrays
USE sedpars
USE syspars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: collect_log_expr
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
INTEGER :: f, i, j
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: sedamount_tot 
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: bed_update_err, sedamount


procname(pglev+1) = 'adjust_one_layer'
CALL log_timer_in()


!
!1. Initialise arrays
!--------------------
!
!---allocate
ALLOCATE (mask(ncloc,nrloc),STAT=errstat)
CALL error_alloc('mask',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (bed_update_err(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('bed_update_err',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (sedamount(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('bed_update_err',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (sedamount_tot(ncloc,nrloc),STAT=errstat)
CALL error_alloc('sedamount_tot',2,(/ncloc,nrloc/),kndrtype)

bed_update_err = 0.0; sedamount = 0.0; sedamount_tot = 0.0

!
!2. check for mass deficit 
!-------------------------
!
!2.1 Sediment availability per fraction in upper layer
!-----------------------------------------------------
!
!---using active layer
IF (iopt_morph_active_layer.GT.0) THEN
   j_211: DO j=1,nrloc
   i_211: DO i=1,ncloc
      IF (maskatc_int(i,j)) THEN
         f_2111: DO f=1,nf
            sedamount(i,j,f) = MERGE(sed_avail(i,j,f),&
                          & bed_fraction(i,j,1,f)*active_layer_thickness(i,j),&
                          & active_layer_thickness(i,j).EQ.0.0) 
         ENDDO f_2111
      ENDIF
   ENDDO i_211
   ENDDO j_211

!---using whole bed layer
ELSE
   f_212: DO f=1,nf
      WHERE (maskatc_int)
         sedamount(:,:,f) = sed_avail(:,:,f)
      END WHERE
   ENDDO f_212
ENDIF

!---total amount
WHERE (maskatc_int)
   sedamount_tot = SUM(sedamount,DIM=3)
END WHERE

!
!2.2 Check where bed update cuts below the fixed layer
!-----------------------------------------------------
!

f_220: DO f=1,nf
   WHERE (maskatc_int)
      bed_update_err(:,:,f) = -(bed_update_dep(:,:,f)+bed_update_ero(:,:,f)+&
                              & sedamount(:,:,f))
   ELSEWHERE
      bed_update_err(:,:,f) = 0.0
   END WHERE
ENDDO f_220

!
!3. Apply bed update correction scheme (if needed)
!-------------------------------------------------
!

flag = collect_log_expr(ANY(bed_update_err.GT.0.0),1)
IF (flag) CALL mass_correct_update(bed_update_err,sedamount)

!
!4. Bed layer update
!-------------------
!
!4.1 Store old values
!--------------------
!

bed_layer_thickness_old = bed_layer_thickness
bed_fraction_old = bed_fraction(1:ncloc,1:nrloc,:,:)

!
!4.2 Bed update
!--------------
!

WHERE (maskatc_int)
   bed_update = MAX(SUM(bed_update_dep+bed_update_ero,DIM=3),-sedamount_tot)
ELSEWHERE
   bed_update = 0.0
END WHERE

!
!4.3 Case of empty sediment layer
!--------------------------------
!

mask = bed_update+sed_avail_tot.LE.0.0.AND.maskatc_int
WHERE (mask)
   bed_layer_thickness(:,:,1) = 0.0
   sed_avail_tot = 0.0
END WHERE

f_430: DO f=1,nf
   WHERE (mask)
      bed_fraction(1:ncloc,1:nrloc,1,f) = 0.0
      sed_avail(:,:,f) = 0.0
   END WHERE
END DO f_430
WHERE (mask)
   sed_avail_tot = 0.0 
END WHERE

!
!4.4 Bed layer thickness
!-----------------------
!

mask = bed_update+sed_avail_tot.GT.0.0.AND.maskatc_int
WHERE (mask)
   bed_layer_thickness(:,:,1) = bed_layer_thickness_old(:,:,1) + bed_update
   sed_avail_tot = bed_layer_thickness(:,:,1)
END WHERE

!
!4.5 Bed layer fractions and sediment amounts
!--------------------------------------------
!

f_450: DO f=1,nf
   WHERE (mask)
      bed_fraction(1:ncloc,1:nrloc,1,f) = &
          & (bed_fraction_old(:,:,1,f)*bed_layer_thickness_old(:,:,1)+&
          &  bed_update_dep(:,:,f)+bed_update_ero(:,:,f))/&
          &  bed_layer_thickness(:,:,1)
      sed_avail(:,:,f) = bed_fraction(1:ncloc,1:nrloc,1,f)*&
                       & bed_layer_thickness(:,:,1)
   END WHERE
ENDDO f_450
WHERE (mask)
   sed_avail_tot = bed_layer_thickness(:,:,1)
END WHERE

!
!5. Median particle diameter
!---------------------------
!

CALL median_particle_diameter(d50_bed,50.0)

!
!6. Active layer thickness
!-------------------------
!

IF (iopt_morph_active_layer.EQ.1) CALL calculate_active_layer_thickness

!
!7. Deallocate arrays
!--------------------
!

DEALLOCATE (bed_update_err,mask,sedamount,sedamount_tot)

CALL log_timer_out()


RETURN

END SUBROUTINE adjust_one_layer

!========================================================================

SUBROUTINE adjust_multi_layers
!************************************************************************
!
! *adjust_multi_layers* adjust bed layer thickness and bed fractions for
!                       deposition and erosion separately and apply sorting in
!                       case of a multi-layer bed model
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.11.2
!
! Description - procedure is invoked to conservation of fraction masses
!
! Reference -
!
! Calling program - morphology_equation
!
! External calls - adjust_layer_top, exchange_fluxes_rubberink,
!                  mass_corrector_update
!
! Moduile calls - collect_log_expr, error_alloc
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphpars
USE morphswitches
USE sedarrays
USE sedpars
USE syspars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: collect_log_expr 
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag, similarity
INTEGER :: i, j, k, f, layers_removed
INTEGER, DIMENSION(nf) :: layers_removed_frac
REAL, PARAMETER :: thickmin = 1.0E-06
REAL :: bed_update_dep_all, bed_update_ero_all, thickness
REAL, DIMENSION(nb) :: bed_layer_depth, bed_layer_thick_old 
REAL, DIMENSION(nf) :: bed_fraction_dep
REAL, DIMENSION(nb,nf) :: bed_frac_old
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: sedamount_tot
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: bed_update_err, sedamount


procname(pglev+1) = 'adjust_multi_layers'
CALL log_timer_in()

!
!1. Initialise arrays
!--------------------
!
!---allocate
ALLOCATE (bed_update_err(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('bed_update_err',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (sedamount(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('sedamount',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (sedamount_tot(ncloc,nrloc),STAT=errstat)
CALL error_alloc('sedamount_tot',2,(/ncloc,nrloc/),kndrtype)

bed_update_err = 0.0; sedamount = 0.0; sedamount_tot = 0.0

!
!2. Adjust deposition rates in case net erosion exceeds the active layer
!-----------------------------------------------------------------------
!
!2.1 Amount of sediment within active or fixed layer
!---------------------------------------------------
!
!2.1.1 Active layer
!------------------
!

IF (iopt_morph_active_layer.GT.0) THEN

   j_210: DO j=1,nrloc
   i_210: DO i=1,ncloc

      IF (maskatc_int(i,j)) THEN

         IF (active_layer_thickness(i,j).GT.0.0.AND.&
           & sed_avail_tot(i,j).GT.active_layer_thickness(i,j)) THEN

            bed_layer_depth(1) = bed_layer_thickness(i,j,1)
            sedamount(i,j,:) = 0.0
            k = 1
            DO WHILE (bed_layer_depth(k).LT.active_layer_thickness(i,j))
               sedamount(i,j,:) = sedamount(i,j,:) + &
                              & bed_fraction(i,j,k,:)*bed_layer_thickness(i,j,k)
               k = k + 1
               bed_layer_depth(k) = bed_layer_depth(k-1) + &
                                  & bed_layer_thickness(i,j,k)
            ENDDO
            IF (k.EQ.1) THEN
               sedamount(i,j,:) = bed_fraction(i,j,1,:)*&
                                & active_layer_thickness(i,j)
            ELSE
               sedamount(i,j,:) = sedamount(i,j,:) + bed_fraction(i,j,k,:)*&
                            & (active_layer_thickness(i,j)-bed_layer_depth(k-1))
            ENDIF
            sedamount_tot(i,j) = SUM(sedamount(i,j,:))

         ELSE
            sedamount(i,j,:) = sed_avail(i,j,:)
            sedamount_tot(i,j) = sed_avail_tot(i,j)
         ENDIF

      ENDIF

   ENDDO i_210
   ENDDO j_210

!
!2.1.2 Fixed layer
!-----------------
!

ELSE

   f_212: DO f=1,nf
      WHERE (maskatc_int)
         sedamount(:,:,f) = sed_avail(:,:,f)
      END WHERE
   ENDDO f_212
   WHERE (maskatc_int)
      sedamount_tot = sed_avail_tot
   END WHERE

ENDIF

!
!2.2 Check where bed update cuts below the fixed layer
!-----------------------------------------------------
!

f_220: DO f=1,nf
   WHERE (maskatc_int)
      bed_update_err(:,:,f) = -(bed_update_dep(:,:,f)+bed_update_ero(:,:,f)+&
                              & sedamount(:,:,f))
   ELSEWHERE
      bed_update_err(:,:,f) = 0.0
   END WHERE
ENDDO f_220

!
!2.3 Apply bed update correction scheme (if needed)
!-------------------------------------------------
!

flag = collect_log_expr(ANY(bed_update_err.GT.0.0),1)
IF (flag) CALL mass_correct_update(bed_update_err,sedamount)

!
!3. Bed layer update
!-------------------
!

WHERE (maskatc_int)
   bed_update = MAX(SUM(bed_update_dep+bed_update_ero,DIM=3),-sedamount_tot)
ELSEWHERE
   bed_update = 0.0
END WHERE

!
!4. Adjust layers
!----------------
!

j_400: DO j=1,nrloc
i_400: DO i=1,ncloc

   IF (maskatc_int(i,j)) THEN

!
!4.1 Adjust for deposition
!-------------------------
!
!4.1.1 Store old values
!----------------------
!

      bed_update_dep_all = SUM(bed_update_dep(i,j,:))

      IF (bed_update_dep_all.GT.0.0) THEN

         bed_layer_thick_old = bed_layer_thickness(i,j,:)
         bed_frac_old = bed_fraction(i,j,:,:)
         bed_fraction_dep = bed_update_dep(i,j,:)/bed_update_dep_all

!
!4.1.2 Check similarity
!----------------------
!

         IF (bed_layer_thick_old(1).GT.0.0) THEN
            similarity = .TRUE.
            f_412: DO f=1,nf
               IF (bed_fraction_dep(f).GT.0.0.AND.bed_fraction_dep(f).LT.&
                & (bed_fraction(i,j,1,f)-similarity_range).OR.&
                & bed_fraction_dep(f).GT.&
                & (bed_fraction(i,j,1,f)+similarity_range)) THEN
                  similarity = .FALSE.
               ENDIF
            ENDDO f_412
         ELSE
            similarity = .FALSE.
         ENDIF

!
!4.1.3 First layer
!-----------------
!

         IF (similarity) THEN
!           ---deposition in upper layer
            bed_layer_thickness(i,j,1) = bed_layer_thick_old(1) + &
                                       & bed_update_dep_all
            bed_fraction(i,j,1,:) = (bed_frac_old(1,:)*bed_layer_thick_old(1)+&
                                   & bed_update_dep(i,j,:))/&
                                  & (bed_layer_thick_old(1)+bed_update_dep_all)
         ELSE
!           ---deposition in new layer
            bed_layer_thickness(i,j,1) = bed_update_dep_all
            bed_fraction(i,j,1,:) = bed_fraction_dep
         ENDIF

!
!4.1.4 Shifting layers 2 to nb-1
!-------------------------------
!

         IF (.NOT.similarity) THEN
            k_414: DO k=2,nb-1
               bed_layer_thickness(i,j,k) = bed_layer_thick_old(k-1)
               bed_fraction(i,j,k,:) = bed_frac_old(k-1,:)
            ENDDO k_414

!
!4.1.5 Merging of bottom layers
!------------------------------
!

            bed_layer_thickness(i,j,nb) = bed_layer_thick_old(nb) + &
                                        & bed_layer_thick_old(nb-1)
            IF (bed_layer_thick_old(nb).GT.0.0.AND.&
              & bed_layer_thick_old(nb-1).GT.0.0) THEN
               bed_fraction(i,j,nb,:) = &
                           & (bed_frac_old(nb,:)*bed_layer_thick_old(nb)+&
                           &  bed_frac_old(nb-1,:)*bed_layer_thick_old(nb-1))/&
                           &  bed_layer_thickness(i,j,nb)
            ELSEIF (bed_layer_thick_old(nb).EQ.0.0.AND.&
                    bed_layer_thick_old(nb-1).GT.0.0) THEN
               bed_fraction(i,j,nb,:) = bed_frac_old(nb-1,:)
            ELSE
               bed_fraction(i,j,nb,:) = 0.0
            ENDIF
         ENDIF

      ENDIF

!
!4.2 Adjust for erosion
!----------------------
!

      bed_update_ero_all = SUM(bed_update_ero(i,j,:))

      IF (bed_update_ero_all.LT.0.0) THEN

         bed_layer_thick_old = bed_layer_thickness(i,j,:)
         bed_frac_old = bed_fraction(i,j,:,:)

!
!4.2.1 Number of removed layers per fraction
!-------------------------------------------
!

         f_421: DO f=1,nf
            layers_removed_frac(f) = 0
            IF (bed_update_ero(i,j,f).LT.0.0) THEN
               thickness = bed_frac_old(1,f)*bed_layer_thick_old(1) + &
                         & bed_update_ero(i,j,f)
               k = 2
               DO WHILE (k.LE.nb.AND.thickness.LE.thickmin)
                  thickness = thickness + &
                            & bed_frac_old(k,f)*bed_layer_thick_old(k)
                  layers_removed_frac(f) = layers_removed_frac(f) + 1
                  k = k + 1
               ENDDO
            ENDIF
         ENDDO f_421
         layers_removed = MAXVAL(layers_removed_frac)

!
!4.2.2 Eroded layer is smaller than top layer thickness
!------------------------------------------------------
!

         IF (layers_removed.EQ.0) THEN

            bed_layer_thickness(i,j,1) = bed_layer_thick_old(1) + &
                                       & bed_update_ero_all
            IF (bed_layer_thickness(i,j,1).GT.0.0) THEN
               bed_fraction(i,j,1,:) = &
                   & (bed_frac_old(1,:)*bed_layer_thick_old(1)+&
                   &  bed_update_ero(i,j,:))/bed_layer_thickness(i,j,1)
            ELSE
               bed_fraction(i,j,1,:) = 0.0
            ENDIF

!
!4.2.3 Eroded layer is larger than top layer thickness
!-----------------------------------------------------
!

         ELSEIF (layers_removed.LE.nb-1) THEN

!           ---layer_thicknesses
            bed_layer_thickness(i,j,1) = &
             & SUM(bed_layer_thick_old(1:layers_removed+1)) + bed_update_ero_all
            k_4231: DO k=2,nb-layers_removed
               bed_layer_thickness(i,j,k) = &
             & bed_layer_thick_old(k+layers_removed)
            ENDDO k_4231
            k_4232: DO k=nb-layers_removed+1,nb
               bed_layer_thickness(i,j,k) = 0.0
            ENDDO k_4232

!           ---bed fractions
            f_4233: DO f=1,nf
               bed_fraction(i,j,1,f) = &
                & (SUM(bed_layer_thick_old(1:layers_removed+1)*&
                & bed_frac_old(1:layers_removed+1,f))+bed_update_ero(i,j,f))/&
                & bed_layer_thickness(i,j,1)
               k_42331: DO k=2,nb-layers_removed
                  bed_fraction(i,j,k,f) = bed_frac_old(k+layers_removed,f)
               ENDDO k_42331
               k_42332: DO k=nb-layers_removed+1,nb
                  bed_fraction(i,j,k,f) = 0.0
               ENDDO k_42332
            ENDDO f_4233

!
!4.2.4 At least one fraction is removed from the bed layer
!---------------------------------------------------------
!

         ELSE

!           ---layer thicknesses
            bed_layer_thickness(i,j,1) = MAX(0.0,&
                                  & SUM(bed_layer_thick_old)+bed_update_ero_all)
            bed_layer_thickness(i,j,2:nb) = 0.0

!           ---bed fractions
            IF (bed_layer_thickness(i,j,1).GT.0.0) THEN
               f_424: DO f=1,nf
                  bed_fraction(i,j,1,f) = &
                           & (SUM(bed_layer_thick_old*bed_frac_old(:,f))+&
                                & bed_update_ero(i,j,f))/&
                                & bed_layer_thickness(i,j,1)
                  bed_fraction(i,j,2:nb,f) = 0.0
               ENDDO f_424
            ELSE
               bed_fraction(i,j,:,:) = 0.0
            ENDIF

         ENDIF

      ENDIF

   ENDIF

ENDDO i_400
ENDDO j_400

!
!5. Sediment availability
!------------------------
!
!---per fraction
f_510: DO f=1,nf
   WHERE (maskatc_int)
      sed_avail(:,:,f) = SUM(bed_fraction(1:ncloc,1:nrloc,:,f)*&
                           & bed_layer_thickness,DIM=3)
   END WHERE
ENDDO f_510

!---total
WHERE (maskatc_int)
   sed_avail_tot = SUM(bed_layer_thickness,DIM=3)
END WHERE

!
!6. Ribberink exchange fluxes in case of migrating dunes
!-------------------------------------------------------
!

IF (iopt_morph_vert_fluxes.EQ.1) CALL exchange_fluxes_ribberink

!
!7. Adjust top layer to active layer
!-----------------------------------
!

CALL adjust_layer_top

!
!8. Deallocate
!-------------
!

DEALLOCATE (bed_update_err,sedamount,sedamount_tot)

CALL log_timer_out()


RETURN

END SUBROUTINE adjust_multi_layers

!========================================================

SUBROUTINE avalanching
!************************************************************************
!
! *avalanching* smooth out steep bed slopes, enabling coastal or river bank 
!               erosion
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.10
!
! Description -
!
! Reference -
!
! Calling program - morphology_equation
!
! External calls - depth_at_nodes, median_particle_diameter, water_depths
!
! Module calls - combine_mod, distribute_mod, error_alloc
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE morpharrays
USE morphpars
USE morphswitches
USE physpars
USE sedarrays
USE sedids
USE sedpars
USE syspars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: combine_mod, distribute_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag, similarity, sortxdir
CHARACTER (LEN=12) :: cnum
INTEGER :: f, i, idep, iend, iero, iside, istart, iteration, i1, i2, j, jdep, &
         & jend, jero, jside, jstart, j1, j2, k, layers_removed
REAL :: bed_layer_thickness_rest, deltadep, deposition, dyn_slope_dry, &
      & dyn_slope_wet, erosion, real_large, stat_slope_dry, stat_slope_wet, &
      & thickness_rest
INTEGER, DIMENSION(2) :: location
REAL, DIMENSION(2) :: area, update
REAL, DIMENSION(nb) :: bed_layer_thickness_dep, bed_layer_thickness_ero
REAL, DIMENSION(nf) :: bed_fraction_aval
REAL, DIMENSION(nb,nf) :: bed_fraction_dep, bed_frac_height, bed_fraction_ero
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask_slope_atc, mask_slope_atu, &
                                            & mask_slope_atv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: bed_slope_x_atu_glb, &
 & bed_slope_y_atv_glb, bed_update_glb, delxatcglb, delxatuglb, delyatcglb, &
 & delyatvglb, dyn_slope_atu, dyn_slope_atu_glb, dyn_slope_atv, &
 & dyn_slope_atv_glb, sed_avail_tot_glb, stat_slope_atu, stat_slope_atu_glb, &
 & stat_slope_atv, stat_slope_atv_glb
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: bed_layer_thickness_glb
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: bed_fraction_glb


procname(pglev+1) = 'avalanching'
CALL log_timer_in()

!
!1. Initialise parameters and arrays
!-----------------------------------
!
!---parameters
real_large = ABS(real_fill)
stat_slope_dry = TAN(degtorad*stat_angle_dry)
dyn_slope_dry = TAN(degtorad*dyn_angle_dry)
stat_slope_wet = TAN(degtorad*stat_angle_wet)
dyn_slope_wet = TAN(degtorad*dyn_angle_wet)


ALLOCATE (bed_fraction_glb(0:nc,0:nr,nb,nf),STAT=errstat)
CALL error_alloc('bed_fraction_glb',4,(/nc+1,nr+1,nb,nf/),kndrtype)
bed_fraction_glb = 0.0

ALLOCATE (bed_layer_thickness_glb(nc,nr,nb),STAT=errstat)
CALL error_alloc('bed_layer_thickness_glb',3,(/nc,nr,nb/),kndrtype)
bed_layer_thickness_glb = 0.0
   
ALLOCATE (bed_slope_x_atu_glb(nc,nr),STAT=errstat)
CALL error_alloc('bed_slope_x_atu_glb',2,(/nc,nr/),kndrtype)
bed_slope_x_atu_glb = 0.0
   
ALLOCATE (bed_slope_y_atv_glb(nc,nr),STAT=errstat)
CALL error_alloc('bed_slope_y_atv_glb',2,(/nc,nr/),kndrtype)
bed_slope_y_atv_glb = 0.0
      
ALLOCATE (bed_update_glb(nc,nr),STAT=errstat)
CALL error_alloc('bed_update',2,(/nc,nr/),kndrtype)
bed_update_glb = 0.0

ALLOCATE (delxatcglb(nc,nr),STAT=errstat)
CALL error_alloc('delxatcglb',2,(/nc,nr/),kndrtype)
delxatcglb = 0.0

ALLOCATE (delxatuglb(nc,nr),STAT=errstat)
CALL error_alloc('delxatuglb',2,(/nc,nr/),kndrtype)
delxatuglb = 0.0

ALLOCATE (delyatcglb(nc,nr),STAT=errstat)
CALL error_alloc('delyatcglb',2,(/nc,nr/),kndrtype)
delyatcglb = 0.0

ALLOCATE (delyatvglb(nc,nr),STAT=errstat)
CALL error_alloc('delyatvglb',2,(/nc,nr/),kndrtype)
delyatvglb = 0.0

ALLOCATE (dyn_slope_atu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('dyn_slope_atu',2,(/ncloc,nrloc/),kndrtype)
dyn_slope_atu = 0.0

ALLOCATE (dyn_slope_atu_glb(nc,nr),STAT=errstat)
CALL error_alloc('dyn_slope_atu_glb',2,(/nc,nr/),kndrtype)
dyn_slope_atu_glb = 0.0

ALLOCATE (dyn_slope_atv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('dyn_slope_atv',2,(/ncloc,nrloc/),kndrtype)
dyn_slope_atv = 0.0

ALLOCATE (dyn_slope_atv_glb(nc,nr),STAT=errstat)
CALL error_alloc('dyn_slope_atv_glb',2,(/nc,nr/),kndrtype)
dyn_slope_atv_glb = 0.0

ALLOCATE (mask_slope_atc(nc-1,nr-1),STAT=errstat)
CALL error_alloc('mask_slope_atc',2,(/nc-1,nr-1/),kndlog)
mask_slope_atc = .FALSE.

ALLOCATE (mask_slope_atu(nc,nr),STAT=errstat)
CALL error_alloc('mask_slope_atu',2,(/nc,nr/),kndlog)
mask_slope_atu = .FALSE.

ALLOCATE (mask_slope_atv(nc,nr),STAT=errstat)
CALL error_alloc('mask_slope_atv',2,(/nc,nr/),kndlog)
mask_slope_atv = .FALSE.

ALLOCATE (stat_slope_atu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('stat_slope_atu',2,(/ncloc,nrloc/),kndrtype)
stat_slope_atu = 0.0

ALLOCATE (stat_slope_atu_glb(nc,nr),STAT=errstat)
CALL error_alloc('stat_slope_atu',2,(/nc,nr/),kndrtype)
stat_slope_atu_glb = 0.0

ALLOCATE (stat_slope_atv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('stat_slope_atv',2,(/ncloc,nrloc/),kndrtype)
stat_slope_atv = 0.0

ALLOCATE (stat_slope_atv_glb(nc,nr),STAT=errstat)
CALL error_alloc('stat_slope_atv_glb',2,(/nc,nr/),kndrtype)
stat_slope_atv_glb = 0.0

IF (iopt_morph_fixed_layer.EQ.1) THEN   
   ALLOCATE (sed_avail_tot_glb(nc,nr),STAT=errstat)
   CALL error_alloc('sed_avail_tot_glb',2,(/nc,nr/),kndrtype)
   sed_avail_tot_glb = 0.0
ENDIF

!
!2. Combine local to global arrays
!---------------------------------
!

CALL combine_mod(bed_fraction_glb,bed_fraction,(/0,0,1,1/),&
     & iarr_bed_fraction,0.0,commall=.TRUE.)
CALL combine_mod(bed_layer_thickness_glb,bed_layer_thickness,&
     & (/1,1,1/),iarr_bed_layer_thickness,0.0,commall=.TRUE.)
CALL combine_mod(bed_slope_x_atu_glb,bed_slope_x_atu,(/1,1/),&
     & iarr_bed_slope_x_atu,0.0,commall=.TRUE.)
CALL combine_mod(bed_slope_y_atv_glb,bed_slope_y_atv,(/1,1/),&
     & iarr_bed_slope_y_atv,0.0,commall=.TRUE.)
CALL combine_mod(bed_update_glb,bed_update,(/1,1/),&
     & iarr_bed_update,0.0,commall=.TRUE.)
CALL combine_mod(delxatcglb,delxatc(1:ncloc,1:nrloc),(/1,1/),&
     & iarr_delxatc,0.0,commall=.TRUE.)
CALL combine_mod(delxatuglb,delxatu(1:ncloc,1:nrloc),(/1,1/),&
     & iarr_delxatu,0.0,commall=.TRUE.)
CALL combine_mod(delyatcglb,delyatc(1:ncloc,1:nrloc),(/1,1/),&
     & iarr_delyatc,0.0,commall=.TRUE.)
CALL combine_mod(delyatvglb,delyatv(1:ncloc,1:nrloc),(/1,1/),&
     & iarr_delyatv,0.0,commall=.TRUE.)
CALL combine_mod(depmeanglb,depmeanatc,(/1-nhalo,1-nhalo/),iarr_depmeanatc,&
               & depmean_flag,commall=.TRUE.)
IF (iopt_morph_fixed_layer.EQ.1) THEN
   CALL combine_mod(sed_avail_tot_glb,sed_avail_tot,(/1,1/),iarr_sed_avail_tot,&
                  & 0.0,commall=.TRUE.)
ENDIF

!
!3. Critical slope angles
!------------------------
!


j_310: DO j=1,nrloc
i_310: DO i=1,ncloc
   IF (.NOT.ALL(seapoint(i-1:i,j)).OR.&
    & (nodeatu(i,j,1).EQ.1.AND.ALL(seapoint(i-1:i,j)))) THEN
      stat_slope_atu(i,j) = real_large
      dyn_slope_atu(i,j) = real_large
   ELSEIF (.NOT.ANY(nodeatc(i-1:i,j).GT.0)) THEN
      stat_slope_atu(i,j) = stat_slope_dry
      dyn_slope_atu(i,j) = dyn_slope_dry
   ELSE
      stat_slope_atu(i,j) = stat_slope_wet
      dyn_slope_atu(i,j) = dyn_slope_wet
   ENDIF
   IF (.NOT.ALL(seapoint(i,j-1:j)).OR.&
    & (nodeatv(i,j,1).EQ.1.AND.ALL(seapoint(i,j-1:j)))) THEN
      stat_slope_atv(i,j) = real_large
      dyn_slope_atv(i,j) = real_large
   ELSEIF (.NOT.ANY(nodeatc(i,j-1:j).GT.0)) THEN
      stat_slope_atv(i,j) = stat_slope_dry
      dyn_slope_atv(i,j) = dyn_slope_dry
   ELSE
      stat_slope_atv(i,j) = stat_slope_wet
      dyn_slope_atv(i,j) = dyn_slope_wet
   ENDIF
ENDDO i_310
ENDDO j_310

!--combine to global
CALL combine_mod(dyn_slope_atu_glb,dyn_slope_atu,(/1,1/),0,real_fill,&
               & commall=.TRUE.)
CALL combine_mod(dyn_slope_atv_glb,dyn_slope_atv,(/1,1/),0,real_fill,&
               & commall=.TRUE.)
CALL combine_mod(stat_slope_atu_glb,stat_slope_atu,(/1,1/),0,real_fill,&
               & commall=.TRUE.)
CALL combine_mod(stat_slope_atv_glb,stat_slope_atv,(/1,1/),0,real_fill,&
               & commall=.TRUE.)

!
!4. Mask arrays
!--------------
!
!4.1 Create masks at velocity nodes with supercritical bed slopes
!----------------------------------------------------------------
!

mask_slope_atu = ABS(bed_slope_x_atu_glb).GT.stat_slope_atu_glb
mask_slope_atv = ABS(bed_slope_y_atv_glb).GT.stat_slope_atv_glb

!
!4.2 Exclude erosion in case of an empty bed layer
!-------------------------------------------------
!
   
IF (iopt_morph_fixed_layer.EQ.1) THEN
   j_420: DO j=1,nr-1
   i_420: DO i=1,nc-1
      IF (sed_avail_tot_glb(i,j).EQ.0.0) THEN
         IF (mask_slope_atu(i,j).AND.bed_slope_x_atu_glb(i,j).GT.0.0) THEN
            mask_slope_atu(i,j) = .FALSE.
         ENDIF
         IF (mask_slope_atu(i+1,j).AND.bed_slope_x_atu_glb(i+1,j).LT.0.0) THEN
            mask_slope_atu(i+1,j) = .FALSE.
         ENDIF
         IF (mask_slope_atv(i,j).AND.bed_slope_y_atv_glb(i,j).GT.0.0) THEN
            mask_slope_atv(i,j) = .FALSE.
         ENDIF
         IF (mask_slope_atv(i,j+1).AND.bed_slope_y_atv_glb(i,j+1).LT.0.0) THEN
            mask_slope_atv(i,j+1) = .FALSE.
         ENDIF
      ENDIF
   ENDDO i_420
   ENDDO j_420
ENDIF

!
!4.3 Mask at cell centers
!------------------------
!

mask_slope_atc = mask_slope_atu(1:nc-1,1:nr-1).OR.mask_slope_atu(2:nc,1:nr-1)& 
           & .OR.mask_slope_atv(1:nc-1,1:nr-1).OR.mask_slope_atv(1:nc-1,2:nr)

!
!5. Apply avalanching at deepest grid points by iteration
!--------------------------------------------------------
!
!---max. number of iterations
iteration = 0

DO WHILE (ANY(mask_slope_atc).AND.iteration.LT.max_iterations_aval)
   iteration = iteration + 1
   
!
!5.1 Location of deepest point with critical slope
!-------------------------------------------------
!
   
   location = MAXLOC(depmeanglb(1:nc-1,1:nr-1),&
                   & MASK=mask_slope_atc(1:nc-1,1:nr-1))
   i = location(1); j = location(2)
   
!
!5.2 Locate surrounding velocity node with steepest slope
!--------------------------------------------------------
!
   
   IF (MAXVAL(ABS(bed_slope_x_atu_glb(i:i+1,j))).GE.&
     & MAXVAL(ABS(bed_slope_y_atv_glb(i,j:j+1)))) THEN
      flag = ABS(bed_slope_x_atu_glb(i+1,j)).GT.&
           & ABS(bed_slope_x_atu_glb(i,j))
      sortxdir = .TRUE.
      istart  = MERGE(i,i-1,flag)
      iend = MERGE(i+1,i,flag)
      iside = MERGE(i+1,i-1,flag)
      jstart = j; jend = j; jside = j
      
   ELSE
      flag = ABS(bed_slope_y_atv_glb(i,j+1)).GT.&
           & ABS(bed_slope_y_atv_glb(i,j))
      sortxdir = .FALSE.
      jstart  = MERGE(j,j-1,flag)
      jend = MERGE(j+1,j,flag)
      jside = MERGE(j+1,j-1,flag)
      istart = i; iend = i; iside = i
   ENDIF
   
!
!5.3  Location and amount of deposition and erosion
!--------------------------------------------------
!
   
   area(1) = delxatcglb(istart,jstart)*delyatcglb(istart,jstart)
   area(2) = delxatcglb(iend,jend)*delyatcglb(iend,jend)

!  ---X-direction
   IF (sortxdir) THEN
      deltadep = delxatuglb(iend,j)*dyn_slope_atu_glb(iend,j)
      update(1) = area(2)*(ABS(depmeanglb(i,j)-&
           & depmeanglb(iside,jside))-deltadep)/(area(1)+area(2))*&
           & SIGN(1.0,depmeanglb(istart,jstart)-depmeanglb(iend,jend))
      update(2) = -update(1)*area(1)/area(2)
      IF (update(1).LT.0.0) THEN
         iero = istart; idep = iend
         erosion = update(1); deposition = update(2)
      ELSE
         iero = iend; idep = istart
         erosion = update(2); deposition = update(1)
      ENDIF
      jero = j; jdep = j

!  ---Y-direction
   ELSE
      deltadep = delyatvglb(i,jend)*dyn_slope_atv_glb(i,jend)
      update(1) = area(2)*(ABS(depmeanglb(i,j)-&
           & depmeanglb(iside,jside))-deltadep)/(area(1)+area(2))*&
           & SIGN(1.0,depmeanglb(istart,jstart)-depmeanglb(iend,jend))
      update(2) = -update(1)*area(1)/area(2)
      IF (update(1).LT.0.0) THEN
         jero = jstart; jdep = jend
         erosion = update(1); deposition = update(2)
      ELSE
         jero = jend; jdep = jstart
         erosion = update(2); deposition = update(1)
      ENDIF
      iero = i; idep = i

   ENDIF

!  ---erosion limited by available sediment
   IF (iopt_morph_fixed_layer.EQ.1) THEN
      erosion = MAX(bed_update_glb(iero,jero)+erosion,&
                 & -sed_avail_tot_glb(iero,jero)) - bed_update_glb(iero,jero)
   ENDIF

!
!5.4 Bed update
!--------------
!

   bed_update_glb(iero,jero) = bed_update_glb(iero,jero) + erosion
   bed_update_glb(idep,jdep) = bed_update_glb(idep,jdep) + deposition

!
!5.5 Update layer depth and fractions
!------------------------------------
!
!5.5.1 One-layer model
!---------------------
!

   IF (nb.EQ.1) THEN

!     ---erosion
      bed_layer_thickness_glb(iero,jero,1) = &
                         & bed_layer_thickness_glb(iero,jero,1) + erosion

      IF (iopt_morph_fixed_layer.EQ.1) THEN
         sed_avail_tot_glb(iero,jero) = sed_avail_tot_glb(iero,jero) + erosion
      ENDIF

!     ---deposition
      bed_layer_thickness_dep = bed_layer_thickness_glb(idep,jdep,:)
      bed_layer_thickness_glb(idep,jdep,1) = bed_layer_thickness_dep(1) + &
                                           & deposition
      bed_fraction_glb(idep,jdep,1,:) = &
       & (bed_layer_thickness_glb(i,j,1)*bed_fraction_glb(idep,jdep,1,:)+&
        & deposition*bed_fraction_glb(iero,jero,1,:))/&
       & (bed_layer_thickness_glb(idep,jdep,1)+deposition)
      IF (iopt_morph_fixed_layer.EQ.1) THEN
         sed_avail_tot_glb(idep,jdep) = sed_avail_tot_glb(idep,jdep) + &
                                      & deposition
      ENDIF
 
!
!5.5.2 Multi-layer model
!------------------------  
!

   ELSE
   
!
!5.5.2.1 Initialise
!------------------
!
      
      bed_layer_thickness_dep = bed_layer_thickness_glb(idep,jdep,:)
      bed_fraction_dep = bed_fraction_glb(idep,jdep,:,:)
      bed_layer_thickness_ero = bed_layer_thickness_glb(iero,jero,:)
      bed_fraction_ero = bed_fraction_glb(iero,jero,:,:)
      
!
!5.5.2.1 Erosion
!---------------
!
!     ---number of layers removed by erosion
      layers_removed = 0
      thickness_rest = 0.0
      bed_layer_thickness_rest = bed_layer_thickness_ero(1) + erosion
      thickness_rest = ABS(bed_layer_thickness_rest)
      
      IF (bed_layer_thickness_rest.LT.0.0) THEN
         layers_removed = 1
         k = 2
         DO WHILE (k.LE.nb.AND.&
          & thickness_rest.GE.bed_layer_thickness_ero(k))
            thickness_rest = thickness_rest - bed_layer_thickness_ero(k)
            layers_removed = layers_removed + 1
            k = k + 1
         END DO
      ENDIF
 
 !    ---bed layer thickness and fraction
      IF (layers_removed.EQ.0) THEN
         bed_layer_thickness_glb(iero,jero,1) = thickness_rest
      ELSE
         IF (layers_removed.LT.nb) THEN
            bed_layer_thickness_glb(iero,jero,1) = &
                 & bed_layer_thickness_ero(1+layers_removed) - thickness_rest
         ENDIF
         k_55211: DO k=2,nb-layers_removed
            bed_layer_thickness_glb(iero,jero,k) = &
                 & bed_layer_thickness_ero(k+layers_removed)
            bed_fraction_glb(iero,jero,k,:) = &
                 & bed_fraction_ero(k+layers_removed,:)
         ENDDO k_55211
         k_55212: DO k=nb-layers_removed+1,nb
            bed_layer_thickness_glb(iero,jero,k) = 0.0
            bed_fraction_glb(iero,jero,k,:) = 0
         ENDDO k_55212
      ENDIF
 
!     ---sediment availability
      IF (iopt_morph_fixed_layer.EQ.1) THEN
         sed_avail_tot_glb(iero,jero) = sed_avail_tot_glb(iero,jero) + erosion
      ENDIF
      
!
!5.5.2.2 Deposition
!-----------------
!
!     ---fraction distribution of the eroded material
      IF (layers_removed.GT.0) THEN
         k_55221: DO k=1,layers_removed
            bed_frac_height(k,:) = bed_layer_thickness_ero(k)*&
                                 & bed_fraction_ero(k,:)
         ENDDO k_55221
      ENDIF
      f_55222: DO f=1,nf
         IF (layers_removed.EQ.0) THEN
            bed_fraction_aval(f) = bed_fraction_ero(1,f)
         ELSEIF (layers_removed.LT.nb) THEN
            bed_fraction_aval(f) = (SUM(bed_frac_height(1:layers_removed,f))+&
                 & thickness_rest*bed_fraction_ero(layers_removed+1,f))/&
               & (SUM(bed_layer_thickness_ero(1:layers_removed))+thickness_rest)
         ELSE
            bed_fraction_aval(f) = (SUM(bed_frac_height(:,f)))/&
                                 & (SUM(bed_layer_thickness_ero))
         ENDIF
      ENDDO f_55222

!     ---check similarity
      similarity = .TRUE.
      f_55223: DO f=1,nf
         IF (bed_fraction_aval(f).LT.&
              & (bed_fraction_dep(1,f)-similarity_range) .OR.&
              & bed_fraction_aval(f).GT.&
              & (bed_fraction_dep(1,f)+similarity_range)) THEN
            similarity = .FALSE.
         ENDIF
      ENDDO f_55223

!     ---deposition in top layer
      IF (similarity) THEN
!        --bed layer thickness
         bed_layer_thickness_glb(idep,jdep,1) = bed_layer_thickness_dep(1) + &
                                              & deposition
!        --bed fraction
         bed_fraction_glb(idep,jdep,1,:) = (bed_fraction_dep(1,:)*&
              & bed_layer_thickness_dep(1)+bed_fraction_aval*deposition)/&
              & (bed_layer_thickness_dep(1)+deposition)

!     ---deposition in new layer
      ELSE
!        --top layer
         bed_layer_thickness_glb(idep,jdep,1) = deposition
         bed_fraction_glb(idep,jdep,1,:) = bed_fraction_aval
         
!        --shifting layers 2 to nb-1
         k_55224: DO k=2,nb-1
            bed_layer_thickness_glb(idep,jdep,k) = bed_layer_thickness_dep(k-1)
            bed_fraction_glb(idep,jdep,k,:) = bed_fraction_dep(k-1,:)
         ENDDO k_55224
         
!        ----merging of bottom layers
         IF (ALL(bed_layer_thickness_dep(nb-1:nb).EQ.0.0) .OR.&
           & ALL(bed_layer_thickness_dep(nb-1:nb).GT.0.0)) THEN
            bed_layer_thickness_glb(idep,jdep,nb) = &
                & bed_layer_thickness_dep(nb) + bed_layer_thickness_dep(nb-1)
            bed_fraction_glb(idep,jdep,nb,:) = &
                & (bed_fraction_dep(nb,:)*bed_layer_thickness_dep(nb)+&
                &  bed_fraction_dep(nb-1,:)*bed_layer_thickness_dep(nb-1))/ &
                &  (bed_layer_thickness_dep(nb) + bed_layer_thickness_dep(nb-1))
         ELSEIF (bed_layer_thickness_dep(nb).EQ.0.0.AND.&
               & bed_layer_thickness_dep(nb-1).GT.0.0) THEN
            bed_layer_thickness_glb(idep,jdep,nb)= bed_layer_thickness_dep(nb-1)
            bed_fraction_glb(idep,jdep,nb,:) = bed_fraction_dep(nb-1,:)
         ENDIF
         
      ENDIF

!     ---sediment availability
      IF (iopt_morph_fixed_layer.EQ.1) THEN      
         sed_avail_tot_glb(idep,jdep) = sed_avail_tot_glb(idep,jdep) + &
                                      & deposition
      ENDIF
      
   ENDIF
   
!
!5.6 Update bathymery
!--------------------
!

   depmeanglb(iero,jero) = depmeanglb(iero,jero) - erosion
   depmeanglb(idep,jdep) = depmeanglb(idep,jdep) - deposition
   
!
!5.7 Adjust bed slopes and mask arrays
!-------------------------------------
!

   IF (sortxdir) THEN

      bed_slope_x_atu_glb(iend,j) = &
         & (depmeanglb(istart,j)-depmeanglb(iend,j))/delxatuglb(iend,j)
      mask_slope_atu(iend,j) = ABS(bed_slope_x_atu_glb(iend,j)).GT.&
                                 & stat_slope_atu_glb(iend,j)
      IF (seapointglb(istart-1,j)) THEN
         bed_slope_x_atu_glb(istart,j) = &
         & (depmeanglb(istart-1,j)-depmeanglb(istart,j))/delxatuglb(istart,j)
         mask_slope_atu(istart,j) = ABS(bed_slope_x_atu_glb(istart,j)).GT.&
                                      & stat_slope_atu_glb(istart,j)
      ENDIF
      IF (seapointglb(iend+1,j)) THEN
         bed_slope_x_atu_glb(iend+1,j) = &
         & (depmeanglb(iend,j)-depmeanglb(iend+1,j))/delxatuglb(iend+1,j)
         mask_slope_atu(iend+1,j) = ABS(bed_slope_x_atu_glb(iend+1,j)).GT.&
                                      & stat_slope_atu_glb(iend+1,j)
      ENDIF
      
   ELSE

      bed_slope_y_atv_glb(i,jend) = &
         & (depmeanglb(i,jstart)-depmeanglb(i,jend))/delyatvglb(i,jend)
      mask_slope_atv(i,jend) = ABS(bed_slope_y_atv_glb(i,jend)).GT.&
                                 & stat_slope_atv_glb(i,jend)
      IF (seapointglb(i,jstart-1)) THEN
         bed_slope_y_atv_glb(i,jstart) = &
         & (depmeanglb(i,jstart-1)-depmeanglb(i,jstart))/delyatvglb(i,jstart)
         mask_slope_atv(i,jstart) = ABS(bed_slope_y_atv_glb(i,jstart)).GT.&
                                      & stat_slope_atv_glb(i,jstart)
      ENDIF
      IF (seapointglb(i,jend+1)) THEN
         bed_slope_y_atv_glb(i,jend+1) = &
         & (depmeanglb(i,jend)-depmeanglb(i,jend+1))/delyatvglb(i,jend+1)
         mask_slope_atv(i,jend+1) = ABS(bed_slope_y_atv_glb(i,jend+1)).GT.&
                                      & stat_slope_atv_glb(i,jend+1)
      ENDIF

   ENDIF

!  ---exclude erosion in case of an empty bed layer
   IF (iopt_morph_fixed_layer.EQ.1) THEN
      IF (sed_avail_tot_glb(iero,jero).EQ.0.0) THEN
         IF (mask_slope_atu(iero,jero).AND.&
           & bed_slope_x_atu_glb(iero,jero).GT.0.0) THEN
            mask_slope_atu(iero,jero) = .FALSE.
         ENDIF
         IF (mask_slope_atu(iero+1,jero).AND.&
           & bed_slope_x_atu_glb(iero+1,jero).LT.0.0) THEN
            mask_slope_atu(iero+1,jero) = .FALSE.
         ENDIF
         IF (mask_slope_atv(iero,jero).AND.&
           & bed_slope_y_atv_glb(iero,jero).GT.0.0) THEN
            mask_slope_atv(iero,jero) = .FALSE.
         ENDIF
         IF (mask_slope_atv(iero,jero+1).AND.&
           & bed_slope_y_atv_glb(iero,jero+1).LT.0.0) THEN
            mask_slope_atv(iero,jero+1) = .FALSE.
         ENDIF
      ENDIF
   ENDIF

!  ---mask at cell centres
   IF (sortxdir) THEN
      i1 = MAX(istart-1,1); i2 = MIN(iend+1,nc-1)
      mask_slope_atc(i1:i2,j) = mask_slope_atu(i1:i2,j).OR.&
                              & mask_slope_atu(i1+1:i2+1,j).OR.&
                              & mask_slope_atv(i1:i2,j).OR.&
                              & mask_slope_atu(i1:i2,j+1)
   ELSE
      j1 = MAX(jstart-1,1); j2 = MIN(jend+1,nr-1)
      mask_slope_atc(i,j1:j2) = mask_slope_atu(i,j1:j2).OR.&
                              & mask_slope_atu(i,j1+1:j2+1).OR.&
                              & mask_slope_atv(i,j1:j2).OR.&
                              & mask_slope_atu(i,j1+1:j2+1)
   ENDIF

ENDDO

!---write to log file
IF (sedflag) THEN
   IF (iteration.GE.max_iterations_aval) THEN
      WRITE (cnum,'(I12)') COUNT(mask_slope_atc); cnum  = ADJUSTL(cnum)
      WRITE(iosed,*) 'WARNING: Maximum iterations reached in avalanching'
      WRITE(iosed,'(A)') 'Number of remaining grid cells: '//TRIM(cnum)
   ELSEIF (iteration.GT.0) THEN
      WRITE (cnum,'(I12)') iteration; cnum  = ADJUSTL(cnum)
      WRITE (iosed,'(A)') 'Number of iterations performed for &
           & avalanching: '//TRIM(cnum)
   ENDIF
ENDIF

!
!6. Distribute data
!--------------------
!

CALL distribute_mod(bed_slope_x_atu_glb,bed_slope_x_atu,(/1,1/),(/0,0,0,0/),&
     & iarr_bed_slope_x_atu,0.0)
CALL distribute_mod(bed_slope_y_atv_glb,bed_slope_y_atv,(/1,1/),(/0,0,0,0/),&
     & iarr_bed_slope_y_atv,0.0)
CALL distribute_mod(bed_layer_thickness_glb,bed_layer_thickness,&
     & (/1,1,1/),(/0,0,0,0/),iarr_bed_layer_thickness,0.0)
CALL distribute_mod(bed_fraction_glb,bed_fraction,(/0,0,1,1/),&
     & (/1,0,1,0/),iarr_bed_fraction,0.0)
CALL distribute_mod(bed_update_glb,bed_update,(/1,1/),(/0,0,0,0/),&
     & iarr_bed_update,0.0)
CALL distribute_mod(depmeanglb,depmeanatc,(/1-nhalo,1-nhalo/),&
     & (/nhalo,nhalo,nhalo,nhalo/),iarr_depmeanatc,depmean_flag)

!
!7. Update sediment availability
!-------------------------------
!

IF (iopt_morph_fixed_layer.EQ.1) THEN
   !  ---sediment availability per fraction
   f_810: DO f=1,nf
      WHERE (maskatc_int)
         sed_avail(:,:,f) = SUM(bed_fraction(1:ncloc,1:nrloc,:,f)*&
              & bed_layer_thickness,DIM=3)
      END WHERE
   ENDDO f_810

!  ---total
   WHERE (maskatc_int)
      sed_avail_tot = SUM(bed_layer_thickness,DIM=3)
   END WHERE
ENDIF

!
!8. Median particle diameter
!---------------------------
!

CALL median_particle_diameter(d50_bed,50.0)

!
!9. Deallocate arrays
!--------------------
!

DEALLOCATE (bed_fraction_glb,bed_layer_thickness_glb,bed_slope_x_atu_glb,&
          & bed_slope_y_atv_glb,bed_update_glb,delxatcglb,delxatuglb,&
          & delyatcglb,delyatvglb,dyn_slope_atu,dyn_slope_atu_glb,&
          & dyn_slope_atv,dyn_slope_atv_glb,mask_slope_atc,mask_slope_atu,&
          & mask_slope_atv,stat_slope_atu,stat_slope_atu_glb,stat_slope_atv,&
          & stat_slope_atv_glb)
IF (iopt_morph_fixed_layer.EQ.1) DEALLOCATE (sed_avail_tot_glb)

CALL log_timer_out()


RETURN

END SUBROUTINE avalanching

!=============================================================================

SUBROUTINE bed_art_diff(depflux,eroflux)
!*****************************************************************************
!
! *bed_art_diff* add the artificial diffusive flux divergence terms in
!                the Exner equation
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - bed_elevation_equation
!
! External calls - 
!
! Module calls - error_alloc, exchange_mod
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE sedarrays
USE sedpars
USE sedswitches
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
REAL, DIMENSION(ncloc,nrloc,nf), INTENT(INOUT):: depflux, eroflux

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*depflux*    REAL    Source terms in the Exner equation(s)
!*eroflux*    REAL    sink terms in the Exner equation(s)
!
!------------------------------------------------------------------------------
!
!*Local variables

INTEGER :: f, i, ii, iiloc, iu, j, jj, jjloc, jv 
REAL :: eps4, flux
REAL, PARAMETER :: alpha_jameson = 2.0, eps = 1.0E-06, k2_jameson = 8.0, &
                 & k4_jameson = 1./32.
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu, maskatv 
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d, art_diffatu, &
                               & art_diffatv, cbedatu, cbedatv, eps2, &
                               & nuatc, nuatu, nuatv, qsedatu, qsedatv


procname(pglev+1) = 'bed_art_diff'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatu(0:ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc+2,nrloc/),kndlog)
ALLOCATE (maskatv(ncloc,0:nrloc+1),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc+2/),kndlog)
ALLOCATE (array2d(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE(art_diffatu(0:ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('art_diffatu',2,(/ncloc+2,nrloc/),kndrtype)
ALLOCATE(art_diffatv(ncloc,0:nrloc+1),STAT=errstat)
CALL error_alloc('art_diffatv',2,(/ncloc,nrloc+2/),kndrtype)
ALLOCATE(cbedatu(0:ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('cbedatu',2,(/ncloc+2,nrloc/),kndrtype)
ALLOCATE(cbedatv(ncloc,0:nrloc+1),STAT=errstat)
CALL error_alloc('cbedatv',2,(/ncloc,nrloc+2/),kndrtype)
ALLOCATE(eps2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('eps2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE(nuatc(0:ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('nuatu',2,(/ncloc+1,nrloc+1/),kndrtype)
ALLOCATE(nuatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('nuatu',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE(nuatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('nuatv',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (qsedatu(0:ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('qsedatu',2,(/ncloc+2,nrloc/),kndrtype)
ALLOCATE (qsedatv(ncloc,0:nrloc+1),STAT=errstat)
CALL error_alloc('qsedatv',2,(/ncloc,nrloc+2/),kndrtype)

!
!2. Initialise
!-------------
!

art_diffatu = 0.0; art_diffatv = 0.0
cbedatu = 0.0; cbedatv = 0.0

!
!3. Sediment loads
!-----------------
!

SELECT CASE (iopt_sed_mode)

   CASE (1)
!     ---bed load only
      qsedatu = SUM(qbedatu_int,DIM=3)
      qsedatv = SUM(qbedatv_int,DIM=3)
   CASE (2)
!     ---suspended load only
      qsedatu = qsedatu + SUM(qsusatu_int,DIM=3)
      qsedatv = qsedatv + SUM(qsusatv_int,DIM=3)
   CASE (3)
!     ---bed and suspended load
      qsedatu = SUM(qbedatu_int+qsusatu_int,DIM=3)
      qsedatv = SUM(qbedatv_int+qsusatv_int,DIM=3)
    CASE (4)
!     ---total load
      qsedatu = SUM(qtotatu_int,DIM=3)
      qsedatv = SUM(qtotatv_int,DIM=3)

END SELECT

!
!4. Bed propagation speed
!------------------------
!
!4.1 U-nodes
!-----------
!
!4.1.1 Mask array
!----------------
!
!---interior points
maskatu = nodeatu(0:ncloc+1,1:nrloc,1).EQ.2
WHERE (maskatu(1:ncloc,:))
   array2d = depmeanatc(1:ncloc,1:nrloc) - depmeanatc(0:ncloc-1,1:nrloc)
ELSEWHERE
   array2d = 0.0
END WHERE
maskatu(1:ncloc,:) = maskatu(1:ncloc,:).AND.ABS(array2d).GT.eps

!---exchange halos
IF (iopt_MPI.EQ.1) CALL exchange_mod(maskatu,(/0,1/),(/1,1,0,0/),0)

!
!4.1.2 Bed propagation speed
!---------------------------
!
!---interior points
WHERE (maskatu(1:ncloc,:))
   cbedatu(1:ncloc,:) = delxatu(1:ncloc,1:nrloc)/&
                      & (delxatc(0:ncloc-1,1:nrloc)+delxatc(1:ncloc,1:nrloc))* &
                      & (qsedatu(2:ncloc+1,:)-qsedatu(0:ncloc-1,:))/array2d
END WHERE

!---exchange halo sections
IF (nobu.GT.0.AND.iopt_MPI.EQ.1) THEN
   CALL exchange_mod(cbedatu,(/0,1/),(/1,1,0,0/),0)
ENDIF

!---open boundaries
iiloc_412: DO iiloc=1,nobuloc
   i = iobuloc(iiloc); j = jobuloc(iiloc)
   ii = indexobu(iiloc); iu = MERGE(i+1,i-1,westobu(ii))
   cbedatu(i,j) = cbedatu(iu,j)
ENDDO iiloc_412

!
!4.2 V-nodes
!-----------
!
!4.2.1 Mask array
!----------------
!
!---interior points
maskatv = nodeatv(1:ncloc,0:nrloc+1,1).EQ.2
WHERE (maskatv(:,1:nrloc))
   array2d = depmeanatc(1:ncloc,1:nrloc) - depmeanatc(1:ncloc,0:nrloc-1)
ELSEWHERE
   array2d = 0.0
END WHERE
maskatv(:,1:nrloc) = maskatv(:,1:nrloc).AND.ABS(array2d).GT.eps

!---exchange halos
IF (iopt_MPI.EQ.1) CALL exchange_mod(maskatv,(/1,0/),(/0,0,1,1/),0)

!
!4.2.2 Bed propagation speed
!---------------------------
!
!---interior points
WHERE (maskatv(:,1:nrloc))
   cbedatv(:,1:nrloc)= delyatv(1:ncloc,1:nrloc)/&
                     & (delyatc(1:ncloc,0:nrloc-1)+delyatc(1:ncloc,1:nrloc))* &
                     & (qsedatv(:,2:nrloc+1)-qsedatv(:,0:nrloc-1))/array2d
END WHERE

!---exchange halo sections
IF (nobv.GT.0.AND.iopt_MPI.EQ.1) THEN
   CALL exchange_mod(cbedatv,(/1,0/),(/0,0,1,1/),0)
ENDIF

!---open boundaries
jjloc_422: DO jjloc=1,nobvloc
   i = iobvloc(jjloc); j = jobvloc(jjloc)
   jj = indexobv(jjloc); jv = MERGE(j+1,j-1,soutobv(jj))
   cbedatv(i,j) = cbedatv(i,jv)
ENDDO jjloc_422

!
!5. Diffusive fluxes
!-------------------
!
!5.1 U-nodes
!-----------
!
!5.1.1 nu-factor at centres
!--------------------------
!

j_511: DO j=1,nrloc
i_511: DO i=0,ncloc
   IF (nodeatc(i,j).GT.0) THEN
      IF (maskatu(i,j).AND.maskatu(i+1,j)) THEN
         nuatc(i,j) = ABS((depmeanatc(i-1,j)-2.0*depmeanatc(i,j)+&
                         & depmeanatc(i+1,j))/&
                        & (depmeanatc(i-1,j)+2.0*depmeanatc(i,j)+&
                        &  depmeanatc(i+1,j)))
      ELSEIF (maskatu(i,j)) THEN
         nuatc(i,j) = ABS((depmeanatc(i-1,j)-depmeanatc(i,j))/&
                        & (depmeanatc(i-1,j)+3.0*depmeanatc(i,j)))
      ELSEIF (maskatu(i+1,j)) THEN
         nuatc(i,j) = ABS((depmeanatc(i+1,j)-depmeanatc(i,j))/&
                        & (depmeanatc(i+1,j)+3.0*depmeanatc(i,j)))
      ELSE
         nuatc(i,j) = 0.0
      ENDIF
   ENDIF
ENDDO i_511
ENDDO j_511

!
!5.1.2 nu-factor at U-nodes
!------------------------
!

WHERE (maskatu(1:ncloc,:))
   nuatu = MAX(nuatc(0:ncloc-1,1:nrloc),nuatc(1:ncloc,1:nrloc))
END WHERE

!
!5.1.3 Diffusive flux (second order)
!-----------------------------------
!

WHERE (maskatu(1:ncloc,:))
   eps2 = MIN(0.5,k2_jameson*nuatu)
   art_diffatu(1:ncloc,:) = eps2*cbedatu(1:ncloc,:)*&
                   & (depmeanatc(0:ncloc-1,1:nrloc)-depmeanatc(1:ncloc,1:nrloc))
END WHERE

!
!5.1.4 Diffusive flux (fourth order) 
!-----------------------------------
!

j_514: DO j=1,nrloc
i_514: DO i=1,ncloc
   IF (maskatu(i,j).AND.ALL(nodeatc(i-2:i+1,j).GT.0)) THEN
      eps4 = MAX(0.0,k4_jameson-alpha_jameson*nuatu(i,j))
      art_diffatu(i,j) = art_diffatu(i,j) + eps4*cbedatu(i,j)*&
                   & (depmeanatc(i+1,j)-3.0*depmeanatc(i,j)+&
                    & 3.0*depmeanatc(i-1,j)-depmeanatc(i-2,j))
   ENDIF
ENDDO i_514
ENDDO j_514

!
!5.1.5 Exchange halos
!--------------------
!

IF (iopt_MPI.EQ.1) CALL exchange_mod(art_diffatu,(/0,1/),(/1,1,0,0/),0)

!
!5.2 V-nodes
!-----------
!
!5.2.1 nu-factor at centres
!--------------------------
!

j_521: DO j=0,nrloc
i_521: DO i=1,ncloc
   IF (nodeatc(i,j).GT.0) THEN
      IF (maskatv(i,j).AND.maskatv(i,j+1)) THEN
         nuatc(i,j) = ABS((depmeanatc(i,j-1)-2.0*depmeanatc(i,j)+&
                         & depmeanatc(i,j+1))/&
                        & (depmeanatc(i,j-1)+2.0*depmeanatc(i,j)+&
                         & depmeanatc(i,j+1)))
      ELSEIF (maskatv(i,j)) THEN
         nuatc(i,j) = ABS((depmeanatc(i,j-1)-depmeanatc(i,j))/&
                        & (depmeanatc(i,j-1)+3*depmeanatc(i,j)))
      ELSEIF (maskatv(i,j+1)) THEN
         nuatc(i,j) = ABS((depmeanatc(i,j+1)-depmeanatc(i,j))/&
                        & (depmeanatc(i,j+1)+3*depmeanatc(i,j)))
      ELSE
         nuatc(i,j) = 0.0
      ENDIF
   ENDIF
ENDDO i_521
ENDDO j_521

!
!5.2.2 nu-factor at V-nodes
!----- --------------------
!

WHERE (maskatv(:,1:nrloc))
   nuatv = MAX(nuatc(1:ncloc,0:nrloc-1),nuatc(1:ncloc,1:nrloc))
END WHERE

!
!5.2.3 Diffusive flux (second order)
!-----------------------------------
!

WHERE (maskatv(:,1:nrloc))
   eps2 = MIN(0.5,k2_jameson*nuatv)
   art_diffatv(:,1:nrloc) = eps2*cbedatv(:,1:nrloc)*&
                   & (depmeanatc(1:ncloc,0:nrloc-1)-depmeanatc(1:ncloc,1:nrloc))
END WHERE

!
!5.2.4 Diffusive flux (fourth order) 
!-----------------------------------
!

j_524: DO j=1,nrloc
i_524: DO i=1,ncloc
   IF (maskatv(i,j).AND.ALL(nodeatc(i,j-2:j+1).GT.0)) THEN
      eps4 = MAX(0.0,k4_jameson-alpha_jameson*nuatv(i,j))
      art_diffatv(i,j) = art_diffatv(i,j) + eps4*cbedatv(i,j)*&
                      & (depmeanatc(i,j+1)-3.0*depmeanatc(i,j)+&
                       & 3.0*depmeanatc(i,j-1)-depmeanatc(i,j-2))
   ENDIF
ENDDO i_524
ENDDO j_524

!
!5.2.5 Exchange halos
!--------------------
!

IF (iopt_MPI.EQ.1) CALL exchange_mod(art_diffatv,(/1,0/),(/0,0,1,1/),0)

!
!6. Add diffusive fluxes to source/sink terms
!--------------------------------------------
!

f_600: DO f=1,nf

!
!6.1 X-direction
!---------------
!

   j_610: DO j=1,nrloc
   i_610: DO i=1,ncloc+1
      flux = delyatu(i,j)*art_diffatu(i,j)
      IF (flux.GT.0) THEN
         IF (i.LT.ncloc+1) THEN
            depflux(i,j,f) = depflux(i,j,f) + bed_fraction(i,j,1,f)*flux
         ENDIF
         IF (i.GT.1) THEN
            eroflux(i-1,j,f) = eroflux(i-1,j,f) - bed_fraction(i-1,j,1,f)*flux
         ENDIF
      ELSEIF (flux.LT.0.0) THEN
         IF (i.GT.1) THEN
            depflux(i-1,j,f) = depflux(i-1,j,f) - bed_fraction(i-1,j,1,f)*flux
         ENDIF
         IF (i.LT.ncloc+1) THEN
            eroflux(i,j,f) = eroflux(i,j,f) + bed_fraction(i,j,1,f)*flux
         ENDIF
      ENDIF
   ENDDO i_610
   ENDDO j_610

!
!6.2 Y-direction
!---------------
!

   j_620: DO j=1,nrloc+1
   i_620: DO i=1,ncloc
      flux = delxatv(i,j)*art_diffatv(i,j)
      IF (flux.GT.0.0) THEN
         IF (j.LT.nrloc+1) THEN
            depflux(i,j,f) = depflux(i,j,f) + bed_fraction(i,j,1,f)*flux
         ENDIF
         IF (j.GT.1) THEN
            eroflux(i,j-1,f) = eroflux(i,j-1,f) - bed_fraction(i,j-1,1,f)*flux
         ENDIF
      ELSEIF (flux.LT.0.0) THEN
         IF (j.GT.1) THEN
            depflux(i,j-1,f) = depflux(i,j-1,f) - bed_fraction(i,j-1,1,f)*flux
         ENDIF
         IF (j.LT.nrloc+1) THEN
            eroflux(i,j,f) = eroflux(i,j,f) + bed_fraction(i,j,1,f)*flux
         ENDIF
      ENDIF
   ENDDO i_620
   ENDDO j_620

ENDDO f_600

!
!7. Deallocate arrays
!--------------------
!

DEALLOCATE (array2d,art_diffatu,art_diffatv,cbedatu,cbedatv,eps2,nuatc,&
          & maskatu,maskatv,nuatu,nuatv,qsedatu,qsedatv)

CALL log_timer_out()


RETURN

END SUBROUTINE bed_art_diff

!========================================================================

SUBROUTINE bed_elevation_equation
!************************************************************************
!
! *bed_elevation_equation* updates bed level by the Exner equation
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - sediment_sorting_predictor
!
! External calls - bed_art_diff, bed_time_int
!
! Module calls - error_alloc, log_timer_in, log_timer_out
!
!************************************************************************
!
USE dararrays
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphpars
USE morphswitches
USE sedarrays
USE sedpars
USE switches
USE sedswitches
USE syspars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j, f
REAL :: dtmorph, flux
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: depdiff, spacetimefac
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: source, depflux, eroflux


procname(pglev+1) = 'bed_elevation_equation'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (source(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('source',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (depdiff(ncloc,nrloc),STAT=errstat)
CALL error_alloc('depdiff',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (depflux(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('depflux',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (eroflux(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('eroflux',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (spacetimefac(ncloc,nrloc),STAT=errstat)
CALL error_alloc('spacetimefac',2,(/ncloc,nrloc/),kndrtype)

!
!2. Initialise parameters and arrays
!-----------------------------------
!

dtmorph = morph_factor*deltmorph

source = 0.0; depdiff = 0.0; depflux = 0.0; eroflux = 0.0
bed_update_dep = 0.0; bed_update_ero = 0.0; bed_update = 0.0

WHERE (maskatc_int)
   spacetimefac = dtmorph/(garea*(1.0-bed_porosity(:,:,1)))
END WHERE

!
!3. Update bed elevation equation
!--------------------------------
!
!3.1 Flux divergence (source/sink) per fraction
!----------------------------------------------
!

SELECT CASE (iopt_sed_mode)
   
CASE (1,3)
   
!
!3.1.1 Bed load
!--------------
!
   
   f_311: DO f=1,nf
      
!     ---X-direction
      j_3111: DO j=1,nrloc
      i_3111: DO i=1,ncloc+1
         flux = delyatu(i,j)*qbedatu_int(i,j,f)
         IF (flux.GT.0.0) THEN
            IF (i.LT.ncloc+1) depflux(i,j,f) = depflux(i,j,f) + flux
            IF (i.GT.1) eroflux(i-1,j,f) = eroflux(i-1,j,f) - flux
         ELSEIF (flux.LT.0.0) THEN
            IF (i.GT.1) depflux(i-1,j,f) = depflux(i-1,j,f) - flux
            IF (i.LT.ncloc+1) eroflux(i,j,f) = eroflux(i,j,f) + flux
         ENDIF
      ENDDO i_3111
      ENDDO j_3111
      
!     ---Y-direction
      j_3112: DO j=1,nrloc+1
      i_3112: DO i=1,ncloc
         flux = delxatv(i,j)*qbedatv_int(i,j,f)
         IF (flux.GT.0.0) THEN
            IF (j.LT.nrloc+1) depflux(i,j,f) = depflux(i,j,f) + flux
            IF (j.GT.1) eroflux(i,j-1,f) = eroflux(i,j-1,f) - flux
         ELSEIF (flux.LT.0.0) THEN
            IF (j.GT.1) depflux(i,j-1,f) = depflux(i,j-1,f) - flux
            IF (j.LT.nrloc+1) eroflux(i,j,f) = eroflux(i,j,f) + flux
         ENDIF
      ENDDO i_3112
      ENDDO j_3112
      
   ENDDO f_311

CASE (4)
   
!
!3.1.2 Total load
!----------------
!
   
   f_312: DO f=1,nf
      
!     ---X-direction
      j_3121: DO j=1,nrloc
      i_3121: DO i=1,ncloc+1
         flux = delyatu(i,j)*qtotatu_int(i,j,f)
         IF (flux.GT.0.0) THEN
            IF (i.LT.ncloc+1) depflux(i,j,f) = depflux(i,j,f) + flux 
            IF (i.GT.1) eroflux(i-1,j,f) = eroflux(i-1,j,f) - flux
         ELSEIF (flux.LT.0.0) THEN
            IF (i.GT.1) depflux(i-1,j,f) = depflux(i-1,j,f) - flux
            IF (i.LT.ncloc+1) eroflux(i,j,f) = eroflux(i,j,f) + flux
         ENDIF
      ENDDO i_3121
      ENDDO j_3121
      
!     ---Y-direction
      j_3122: DO j=1,nrloc+1
      i_3122: DO i=1,ncloc
         flux = delxatv(i,j)*qtotatv_int(i,j,f)
         IF (flux.GT.0.0) THEN
            IF (j.LT.nrloc+1) depflux(i,j,f) = depflux(i,j,f) + flux
            IF (j.GT.1) eroflux(i,j-1,f) = eroflux(i,j-1,f) - flux
         ELSEIF (flux.LT.0.0) THEN
            IF (j.GT.1) depflux(i,j-1,f) = depflux(i,j-1,f) - flux
            IF (j.LT.nrloc+1) eroflux(i,j,f) = eroflux(i,j,f) + flux
         ENDIF
      ENDDO i_3122
      ENDDO j_3122
      
   ENDDO f_312
   
END SELECT

f_313: DO f=1,nf
   WHERE (maskatc_int)
      depflux(:,:,f) = spacetimefac*depflux(:,:,f)
      eroflux(:,:,f) = spacetimefac*eroflux(:,:,f)
   ELSEWHERE
      depflux(:,:,f) = 0.0
      eroflux(:,:,f) = 0.0
   END WHERE
ENDDO f_313

!
!3.2 Add artificial diffusion term
!---------------------------------
!

IF (iopt_morph_hdiff.EQ.1) CALL bed_art_diff(depflux,eroflux)

!
!3.3 Other sources/sinks
!-----------------------
!

f_330: DO f=1,nf

!  ----depôsition/erosion
   IF (iopt_sed_mode.EQ.2.OR.iopt_sed_mode.EQ.3) THEN
      WHERE (maskatc_int)
         source(:,:,f) = -bottom_sed_flux(:,:,f)*dtmorph/&
                        & (1.0-bed_porosity(:,:,1))
      END WHERE
   ENDIF

!  ---dredging/relocation
   IF (iopt_dar.EQ.1) THEN 
      WHERE (maskatc_int)
         source(:,:,f) = source(:,:,f) + dar_morph(:,:,f)/garea
      END WHERE
   ENDIF

   !  ---add source/sinks
   IF (iopt_sed_mode.EQ.2.OR.iopt_sed_mode.EQ.3.OR.iopt_dar.EQ.1) THEN
      WHERE (source(:,:,f).GT.0.0)
         depflux(:,:,f) = depflux(:,:,f) + source(:,:,f)
      END WHERE
      WHERE (source(:,:,f).LT.0.0)
         eroflux(:,:,f) = eroflux(:,:,f) + source(:,:,f)
      END WHERE
   ENDIF

ENDDO f_330

!
!3.4 Time integration
!--------------------
!

CALL bed_time_int(bed_update_dep,depflux,bed_update_dep_old1,&
                & bed_update_dep_old2)
CALL bed_time_int(bed_update_ero,eroflux,bed_update_ero_old1,&
                & bed_update_ero_old2)

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (depdiff,depflux,eroflux,source,spacetimefac)

CALL log_timer_out()


RETURN

END SUBROUTINE bed_elevation_equation

!========================================================================

SUBROUTINE bed_time_int(bedupdate,flux,bedupdate_old1,bedupdate_old2)
!************************************************************************
!
! *bed_time_int* integrate Exner equation in time using Adams-Bashforth scheme
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.10
!
! Description - type (order) of time-integration scheme is selected by
!               iopt_morph_time_int
!             - routine is separately called for update of source
!               ("deposition") and sink ("erosion") term 
!
! Reference -
!
! Calling program - bed_elevation_equation, bed_elevation_equation
!
! External calls - 
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE morphswitches
USE sedpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
REAL, DIMENSION(ncloc,nrloc,nf), INTENT(IN) :: flux
REAL, DIMENSION(ncloc,nrloc,nf), INTENT(INOUT) :: bedupdate, bedupdate_old1, &
                                                & bedupdate_old2
!
! Name            Type    Purpose
!------------------------------------------------------------------------------
!*bedupdate*      REAL    bed update at the new time step
!*flux*           REAL    rhs terms at t^n on input, t^(n+1) on output
!*bedupdate_old1* REAL    rhs terms at t^(n-1) on input, t^n on output 
!*bedupdate_old2* REAL    rhs terms at t^(n-2) on input, t^(n-1) on output 
!
!*******************************************************************************
!
!*Local variables
!
INTEGER :: f
REAL, PARAMETER :: ab2fac1 = 1.5, ab2fac2 = -0.5
REAL, PARAMETER :: ab3fac1 = 23./12., ab3fac2 = -4./3., ab3fac3 = 5./12.


procname(pglev+1) = 'bed_time_int'
CALL log_timer_in()

f_110: DO f=1,nf

   SELECT CASE (iopt_morph_time_int)
!     ---first order (Euler)
      CASE (1)
         WHERE (maskatc_int)
            bedupdate(:,:,f) = flux(:,:,f)
         END WHERE
!     ---second order Adams-Bashforth
      CASE (2)
         WHERE (maskatc_int)
            bedupdate(:,:,f) =  ab2fac1*flux(:,:,f) + &
                              & ab2fac2*bedupdate_old1(:,:,f)
            bedupdate_old1(:,:,f) = flux(:,:,f)
         END WHERE
!     ---third order Adams-Bashforth
      CASE (3)
         WHERE (maskatc_int)
            bedupdate(:,:,f) = ab3fac1*flux(:,:,f) + &
                             & ab3fac2*bedupdate_old1(:,:,f) + &
                             & ab3fac3*bedupdate_old2(:,:,f)
            bedupdate_old1(:,:,f) = flux(:,:,f)
            bedupdate_old2(:,:,f) = bedupdate_old1(:,:,f)
         END WHERE
   END SELECT

ENDDO f_110

CALL log_timer_out()


RETURN

END SUBROUTINE bed_time_int

!==============================================================================

SUBROUTINE calculate_active_layer_thickness
!******************************************************************************
!
! *calculate_active_layer_thickness* active layer thickness
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.9
!
! Description -
!
! Reference -
!
! Calling program - sediment_sorting_corrector, sediment_sorting_predictor
!
! External calls - 
!
! Module calls - log_timer_in, log_timer_out
!
!************************************************************************
!
USE fluxes
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphpars
USE morphswitches
USE sedarrays
USE sedpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

procname(pglev+1) = 'calculate_active_layer_thickness'
CALL log_timer_in()

!
!1. Active layer thickness
!-------------------------
!

SELECT CASE (iopt_morph_active_layer)

!  ---Armanini (1995)
   CASE (1)
      WHERE (maskatc_int)
         active_layer_thickness = 4.5*d50_bed(1:ncloc,1:nrloc)
      END WHERE
!  ---Ribberink (1987)
   CASE (2)
      WHERE (maskatc_int)
         active_layer_thickness = 0.5*average_bedform_height
      END WHERE
!  ---Harris and Wiberg (1997)
   CASE (3)
      WHERE (maskatc_int)
         active_layer_thickness = MAX(k1_harris*&
                       & (bstresatc_sed(1:ncloc,1:nrloc)-&
                       & SUM(bstres_cr(1:ncloc,1:nrloc,:),3)/REAL(nf)),0.0) + & 
                       & k2_harris*d50_bed(1:ncloc,1:nrloc)
      END WHERE

END SELECT

!
!2. Check whether active layer does not cut into a fixed layer
!-------------------------------------------------------------
!

IF (iopt_morph_fixed_layer.EQ.1) THEN
   WHERE (maskatc_int)
      active_layer_thickness = MIN(active_layer_thickness,sed_avail_tot)
   END WHERE
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE calculate_active_layer_thickness

!=======================================================

SUBROUTINE exchange_fluxes_ribberink
!************************************************************************
!
! *exchange_fluxes_ribberink* exchange of sediment between 
!                             top layers due to migrating dunes
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - sediment_sorting
!
! External calls - log_timer_in, log_timer_out
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphpars
USE sedarrays
USE sedpars
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j
REAL :: exchange_flux, xtmp
REAL, DIMENSION(nb,nf) :: bed_frac_old


procname(pglev+1) = 'exchange_fluxes_ribberink'
CALL log_timer_in()

IF (nf.EQ.2) THEN

   xtmp =  deltmorph*1.22*dune_celerity*trough_probability/&
         & average_bedform_length
   j_210: DO j=1,nrloc
   i_210: DO i=1,ncloc

      IF (maskatc_int(i,j)) THEN
         bed_frac_old = bed_fraction(i,j,:,:)
         exchange_flux = xtmp*active_layer_thickness(i,j)

!        ---change bed fraction distribution in layer 1
         bed_fraction(i,j,1,1) = &
               & (bed_frac_old(1,1)*active_layer_thickness(i,j)-&
               & (0.7*bed_frac_old(1,1)-bed_frac_old(2,1))*exchange_flux)/&
               & active_layer_thickness(i,j)
         bed_fraction(i,j,1,2) = &
               & (bed_frac_old(1,2)*active_layer_thickness(i,j)-&
               & ((1-0.7*bed_frac_old(1,1))*bed_frac_old(1,2)-&
               & bed_frac_old(2,2))*exchange_flux)/active_layer_thickness(i,j)

!       ---change bed fraction distribution in layer 2
         bed_fraction(i,j,2,1) = &
               & (bed_frac_old(2,1)*bed_layer_thickness(i,j,2)+&
              & (0.7*bed_frac_old(1,1)-bed_frac_old(2,1))*exchange_flux)/&
               & bed_layer_thickness(i,j,2)
         bed_fraction(i,j,2,2) = &
               & (bed_frac_old(2,2)*bed_layer_thickness(i,j,2)+&
               & ((1-0.7*bed_frac_old(1,1))*bed_frac_old(1,2)-&
                & bed_frac_old(2,2))*exchange_flux)/bed_layer_thickness(i,j,2)
      ENDIF

   ENDDO i_210
   ENDDO j_210

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE exchange_fluxes_ribberink

!=======================================================

SUBROUTINE mass_balance 
!************************************************************************
!
! *mass_balance* update sediment mass balance over the whole domain
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.10
!
! Description -
!
! Reference -
!
! Calling program - morphology_equation
!
! External calls - 
!
! Module calls - mult_index, sum_vars
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphpars
USE sedarrays
USE sedpars
USE sedswitches
USE switches
USE timepars
USE paral_utilities, ONLY : sum_vars
USE time_routines, ONLY: log_timer_in, log_timer_out, write_time_log
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL, SAVE :: bload, sload, tload
INTEGER :: i, ii, iiloc, isign, j, jj, jjloc, jsign, k
REAL, SAVE :: bedvol_old, ctot_old, dsedobc, dt
REAL :: bedvol, ctotal, dsedvol


procname(pglev+1) = 'mass_balance'
CALL log_timer_in()

!
!1. Initalise at initial time
!----------------------------
!

IF (nt.EQ.0) THEN

!
!1.1 Initialise parameters
!-------------------------
!

   bload = iopt_sed_mode.EQ.1.OR.iopt_sed_mode.EQ.3
   tload = iopt_sed_mode.EQ.4
   sload = iopt_sed_mode.EQ.2.OR.iopt_sed_mode.EQ.3
   dt = deltmorph*morph_factor
   dsedobc = 0.0

!
!1.2 Amount of sediment volume
!----------------------------
!
!   ---sea bed
   bedvol_old = SUM((1.0-bed_porosity(:,:,1))*bed_update_int*garea,&
                   & mask=maskatc_int)

!  ---water column
!   IF (sload) THEN
!      ctot_old = 0.0
!      k_130: DO k=1,nz
!         ctot_old =  ctot_old + SUM(delzatc(1:ncloc,1:nrloc,k)*garea*&
!                                  & ctot(1:ncloc,1:nrloc,k),maskatc_int)
!      ENDDO k_130
!   ENDIF

ENDIF

!
!2. Transport through open boundaries
!------------------------------------
!

IF (nt.GT.0) THEN

!
!2.1 U-open boundaries
!---------------------
!

   iiloc_210: DO iiloc=1,nobuloc
      i = iobuloc(iiloc); j = jobuloc(iiloc)
      ii = indexobu(iiloc); isign = MERGE(1,-1,westobu(ii))

!     ---bed load
      IF (bload) THEN
         dsedobc = dsedobc + dt*isign*delyatu(i,j)*SUM(qbedatu_int(i,j,:))
      ENDIF
         
!     ---total load
      IF (tload) THEN
         dsedobc = dsedobc + dt*isign*delyatu(i,j)*SUM(qtotatu_int(i,j,:))
      ENDIF

!     ---suspended load
!      IF (sload) THEN
!         dsedobc = dsedobc + deltmorph*isign*delyatu(i,j)*&
!                 & SUM(qsusatu_int(i,j,:))
!      ENDIF

   ENDDO iiloc_210

!
!2.2 V-open boundaries
!---------------------
!

   jjloc_220: DO jjloc=1,nobvloc
      i = iobvloc(jjloc); j = jobvloc(jjloc)
      jj = indexobv(jjloc); jsign = MERGE(1,-1,soutobv(jj))

!     ---bed load
      IF (bload) THEN
         dsedobc = dsedobc + dt*jsign*delxatv(i,j)*SUM(qbedatv_int(i,j,:))
      ENDIF
         
!     ---total load
      IF (tload) THEN
         dsedobc = dsedobc + dt*jsign*delxatv(i,j)*SUM(qtotatv_int(i,j,:))
      ENDIF

!     ---suspended load
      IF (sload) THEN
         dsedobc = dsedobc + deltmorph*jsign*delxatv(i,j)*&
                 & SUM(qsusatv_int(i,j,:))
      ENDIF

   ENDDO jjloc_220

ENDIF

!
!3. Mass balance
!---------------
!

IF (nt.GT.0.AND.mult_index(nt,icsedbal)) THEN

!  ---bed sediment
   bedvol = SUM((1.0-bed_porosity(:,:,1))*bed_update_int*garea,mask=maskatc_int)
   dsedvol = bedvol - bedvol_old

!  ---water column sediment
!   IF (sload) THEN
!      ctotal = 0.0
!      k_310: DO k=1,nz
!         ctotal =  ctotal + SUM(delzatc(1:ncloc,1:nrloc,k)*garea*&
!                              & ctot(1:ncloc,1:nrloc,k),maskatc_int)
!      ENDDO k_310
!      dsedvol = dsedvol + (ctotal-ctot_old)
!   ENDIF

!  ---sum over all process domains
   IF (iopt_MPI.EQ.1) THEN
      CALL sum_vars(dsedvol,sediment_vol,0,commall=.TRUE.)
      CALL sum_vars(dsedobc,sediment_obc,0,commall=.TRUE.)
   ELSE
      sediment_vol  = dsedvol
      sediment_obc  = dsedobc
   ENDIF

!  ---sediment balance
   sediment_balance = sediment_vol - sediment_obc

!  ---write to log file 
   IF (sedflag) THEN
      CALL write_time_log(iosed,sedlog_date)
      WRITE (iosed,9001) 'Net change of sediment amount within the domain: ', &
                        & sediment_vol, ' m^3'
      WRITE (iosed,9001) 'Net amount of sediments transported across '//&
                       & 'open boundaries: ', sediment_obc, ' m^3'
      WRITE (iosed,9001) 'Sediment balance error: ', sediment_balance, ' m^3'
   ENDIF

!  ---reset
   dsedobc = 0.0
   bedvol_old = bedvol
!   IF (sload) ctot_old = ctotal

ENDIF

CALL log_timer_out()


RETURN

9001 FORMAT (A,G15.7,A)

END SUBROUTINE mass_balance

!========================================================================

SUBROUTINE mass_correct_update(bed_update_err,sedamount)
!************************************************************************
!
! *mass_correct_update* spread net excess sediment over neighbouring cells
!                       to ensure mass conservation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.10
!
! Description -
!
! Reference -
!
! Calling program - adjust_multi_layer, adjust_one_layer
!
! Module calls - combine_mod, distribute_mod, error_abort, error_alloc
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphpars
USE morphswitches
USE sedids
USE sedpars
USE syspars
USE error_routines, ONLY: error_abort, error_alloc
USE paral_comms, ONLY: combine_mod, distribute_mod
USE paral_utilities, ONLY: sum_vars
USE time_routines, ONLY: log_timer_in, log_timer_out, write_time_log

IMPLICIT NONE

!
!* Arguments
!
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nf) :: bed_update_err, sedamount

! Name              Type  Purpose
!------------------------------------------------------------------------------
!*bed_update_err*  REAL Excess amount of sediment to be spread over
!                       neighbouring cells                            [m^3/m^3]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: search
CHARACTER (LEN=12) :: citer
INTEGER :: f, ff, i, iteration, i1, i2, j, j1, j2
REAL :: sederr, weightssum
INTEGER, DIMENSION(2) :: location
INTEGER, DIMENSION(4) :: nhdist
REAL, DIMENSION(2,2) :: weights
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskglb
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: bed_update_dep_glb, &
                                           & bed_update_err_glb, &
                                           & bed_update_ero_glb, sedamount_glb


procname(pglev+1) = 'mass_correct_update'
CALL log_timer_in()

!
!1. No correction is applied
!---------------------------
!

IF (iopt_morph_corr.EQ.0) THEN
   WHERE (bed_update_err.GT.0.0)
      bed_update_ero = -bed_update_dep - sedamount
   END WHERE
   CALL sum_vars(bed_update_err,sederr,0,mask=maskatc_int)
   sediment_err = sediment_err + sederr
   IF (sedflag) THEN
      CALL write_time_log(iosed,sedlog_date)
      WRITE (iosed,9001) 'Amount of non-allocated sediment: ', sederr, ' m^3'
      WRITE (iosed,9001) 'Total amount of non-allocated sediment: ', &
                       & sediment_err, ' m^3'
   ENDIF
   GOTO 1000
ENDIF

!
!2. Allocate arrays
!------------------
!

ALLOCATE (maskglb(0:nc+1,0:nr+1),STAT=errstat)
CALL error_alloc('maskglb',2,(/nc+2,nr+2/),kndlog)
maskglb = .FALSE.
ALLOCATE (bed_update_dep_glb(nc,nr,nf),STAT=errstat)
CALL error_alloc('bed_update_dep_glb',3,(/nc,nr,nf/),kndrtype)
bed_update_dep_glb = 0.0
ALLOCATE (bed_update_err_glb(nc,nr,nf),STAT=errstat)
CALL error_alloc('bed_update_err_glb',3,(/nc,nr,nf/),kndrtype)
bed_update_err_glb = 0.0
ALLOCATE (bed_update_ero_glb(nc,nr,nf),STAT=errstat)
CALL error_alloc('bed_update_ero_glb',3,(/nc,nr,nf/),kndrtype)
bed_update_ero_glb = 0.0
ALLOCATE (sedamount_glb(nc,nr,nf),STAT=errstat)
CALL error_alloc('sedamount_glb',3,(/nc,nr,nf/),kndrtype)
sedamount_glb = 0.0

!
!3. Combine
!----------
!

CALL combine_mod(maskglb(1:nc,1:nr),maskatc_int,(/1,1/),0,.FALSE.,&
               & commall=.TRUE.)
CALL combine_mod(bed_update_dep_glb,bed_update_dep,(/1,1,1/),0,0.0,&
               & commall=.TRUE.)
CALL combine_mod(bed_update_err_glb,bed_update_err,(/1,1,1/),0,0.0,&
               & commall=.TRUE.)
CALL combine_mod(bed_update_ero_glb,bed_update_ero,(/1,1,1/),0,0.0,&
               & commall=.TRUE.)
CALL combine_mod(sedamount_glb,sedamount,(/1,1,1/),0,0.0,commall=.TRUE.)

!
!4. Extract excess erosion over neighbouring cells
!-------------------------------------------------
!

iteration = 0

f_400: DO f=1,nf

   DO WHILE (ANY(bed_update_err_glb(:,:,f).GT.0.0))
      iteration = iteration + 1
      location = MAXLOC(bed_update_err_glb(:,:,f),&
                      & MASK=bed_update_err_glb(:,:,f).GT.0.0)
      i = location(1); j = location(2)

!
!4.1 Search for neighbouring cells
!---------------------------------
!
!     ---West
      i1 = i - 1; search = .TRUE.
      DO WHILE (search)
         IF (.NOT.maskglb(i1,j)) THEN
            i1 = 0; search = .FALSE.
         ELSEIF (bed_update_err_glb(i1,j,f).GE.0.0) THEN
            i1 = i1 - 1
         ELSE
            search = .FALSE.
         ENDIF
      ENDDO

!     ---East
      i2 = i + 1; search = .TRUE.
      DO WHILE (search)
         IF (.NOT.maskglb(i2,j)) THEN
            i2 = 0; search = .FALSE.
         ELSEIF (bed_update_err_glb(i2,j,f).GE.0.0) THEN
            i2 = i2 + 1
         ELSE
            search = .FALSE.
         ENDIF
      ENDDO

!     ---South
      j1 = j - 1; search = .TRUE.
      DO WHILE (search)
         IF (.NOT.maskglb(i,j1)) THEN
            j1 = 0; search = .FALSE.
         ELSEIF (bed_update_err_glb(i,j1,f).GE.0.0) THEN
            j1 = j1 - 1
         ELSE
            search = .FALSE.
         ENDIF
      ENDDO

!     ---North
      j2 = j + 1; search = .TRUE.
      DO WHILE (search)
         IF (.NOT.maskglb(i,j2)) THEN
            j2 = 0; search = .FALSE.
         ELSEIF (bed_update_err_glb(i,j2,f).GE.0.0) THEN
            j2 = j2 + 1
         ELSE
            search = .FALSE.
         ENDIF
      ENDDO

!
!4.2  No neighbouring cells found
!--------------------------------
!

      IF (i1+i2+j1+j2.EQ.0) THEN

!        ---sediment error
         IF (iopt_morph_corr.EQ.1) THEN
            WHERE (bed_update_err_glb(:,:,f).GT.0.0)
               bed_update_ero_glb(:,:,f) = &
                & - bed_update_dep_glb(:,:,f) - sedamount_glb(:,:,f)
            END WHERE
         ENDIF
         sederr = 0.0
         ff_421: DO ff=1,nf
            sederr = sederr + SUM(bed_update_err_glb(:,:,ff),&
                                & MASK=maskglb(1:nc,1:nr))
         ENDDO ff_421
         sediment_err = sediment_err + sederr

!        ---search is terminated, no abort
         IF (iopt_morph_corr.EQ.1) THEN
            IF (sedflag) THEN
                CALL write_time_log(iosed,sedlog_date)
               WRITE (iosed,'(A)') 'All grid cells have a negative mass deficit'
               WRITE (iosed,'(A)') 'Sediment mass conservation cannot be '//&
                                 & 'maintained'
               WRITE (iosed,9001) 'Amount of non-allocated sediment: ', &
                                 & sederr, ' m^3'
               WRITE (iosed,9001) 'Total amount of non-allocated sediment: ', &
                                & sediment_err, ' m^3'
            ENDIF
            CYCLE f_400

!        ---search is terminated, propgram aborts         
         ELSEIF (iopt_morph_corr.EQ.2) THEN
            nerrs = 1
            IF (errchk) THEN
               WRITE (ioerr,'(A)') 'All grid cells have a negative mass deficit'
               WRITE (iosed,'(A)') 'Sediment mass conservation cannot be '//&
                                 & 'maintained'
               WRITE (iosed,9001) 'Amount of non-allocated sediment: ', &
                                & sederr, ' m^3'
               WRITE (iosed,9001) 'Total amount of non-allocated sediment: ', &
                                & sediment_err, ' m^3'
            ENDIF
            CALL error_abort('mass_correct_update',ierrno_runval)
         ENDIF

      ENDIF

!
!4.3 Extract at neighbouring cells
!---------------------------------
!

      weights = 0.0
      IF (i1.GT.0) weights(1,1) = bed_update_err_glb(i1,j,f)
      IF (i2.GT.0) weights(2,1) = bed_update_err_glb(i2,j,f)
      IF (j1.GT.0) weights(1,2) = bed_update_err_glb(i,j1,f)
      IF (j2.GT.0) weights(2,1) = bed_update_err_glb(i,j2,f)
      weightssum = SUM(weights)
      weights = bed_update_err_glb(i,j,f)*weights/weightssum
      IF (i1.GT.0) THEN
         bed_update_ero_glb(i1,j,f) = bed_update_ero_glb(i1,j,f) - weights(1,1)
         bed_update_err_glb(i1,j,f) = -(bed_update_dep_glb(i1,j,f)+&
                                      & bed_update_ero_glb(i1,j,f)+&
                                      & sedamount_glb(i1,j,f))
      ENDIF
      IF (i2.GT.0) THEN
         bed_update_ero_glb(i2,j,f) = bed_update_ero_glb(i2,j,f) - weights(2,1)
         bed_update_err_glb(i2,j,f) = -(bed_update_dep_glb(i2,j,f)+&
                                      & bed_update_ero_glb(i2,j,f)+&
                                      & sedamount_glb(i2,j,f))
      ENDIF
      IF (j1.GT.0) THEN
         bed_update_ero_glb(i,j1,f) = bed_update_ero_glb(i,j1,f) - weights(1,2)
         bed_update_err_glb(i,j1,f) = -(bed_update_dep_glb(i,j1,f)+&
                                      & bed_update_ero_glb(i,j1,f)+&
                                      & sedamount_glb(i,j1,f))
      ENDIF
      IF (j2.GT.0) THEN
         bed_update_ero_glb(i,j2,f) = bed_update_ero_glb(i,j2,f) - weights(2,2)
         bed_update_err_glb(i,j2,f) = -(bed_update_dep_glb(i,j2,f)+&
                                      & bed_update_ero_glb(i,j2,f)+&
                                      & sedamount_glb(i,j2,f))
      ENDIF
      bed_update_ero_glb(i,j,f) = -(bed_update_dep_glb(i,j,f)+&
                                  & sedamount_glb(i,j,f))
      bed_update_err_glb(i,j,f) = 0.0

   ENDDO

ENDDO f_400

!---write to log file
IF (iteration.GT.0.AND.sedflag) THEN
   WRITE (citer,'(I12)') iteration; citer = ADJUSTL(citer)
   CALL write_time_log(iosed,sedlog_date)
   WRITE (iosed,'(A)') 'Number of iterations to correct for mass deficit: '&
                     & //TRIM(citer)
ENDIF

!
!5. Redistribute on local processes
!-----------------------------------
!

nhdist = 0
CALL distribute_mod(bed_update_ero_glb,bed_update_ero,(/1,1,1/),nhdist,&
                  & iarr_bed_update_ero,0.0)

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (maskglb,bed_update_dep_glb,bed_update_err_glb,bed_update_ero_glb,&
          & sedamount_glb)

1000 CALL log_timer_out()


RETURN

9001 FORMAT (A,G15.7,A)

END SUBROUTINE mass_correct_update

!============================================================================

SUBROUTINE morphology_accel
!****************************************************************************
!
! *morphology_accel* update sediment and morphology using the tidal
!                    acceleration method
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - astronomical_tides, define_phsics, hydrodynamic_equations,
!                  mask_function, morphological_step, store_depths_old,
!                  time_series, usrdef_output, weirs_mask, write_morphics,
!                  write_phsics, write_sedics
!
! Module calls - monitor_files, mult_index, update_time
!
!************************************************************************
!
USE iopars
USE morphpars
USE morpharrays
USE switches
USE timepars
USE inout_routines, ONLY: monitor_files
USE time_routines, ONLY: log_timer_in, log_timer_out,update_time
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=12) :: cval
INTEGER, SAVE :: morph_end, morph_start
INTEGER :: n, nosteps, npcc, ntmod, nperiod, noperiods


procname(pglev+1) = 'morphology_accel'
CALL log_timer_in(npcc)

!
!1. Initialise
!-------------
!
!---number of characteristic periods
IF (mult_index(nstep,nstep_hydro*morph_steps)) THEN
  noperiods = nstep/(nstep_hydro*morph_steps)
ELSE
  noperiods = nstep/(nstep_hydro*morph_steps) + 1
ENDIF

!---time index at the end of the last hydrodynamic period 
hydro_end = (noperiods-1)*nstep_hydro*morph_steps+nstep_hydro

!---time steps for storing hydrodynamical data
nosteps = nstep_hydro/number_tidal_steps
ntmod = MOD(nstep_hydro,number_tidal_steps)
n_110: DO n=1,number_tidal_steps
   IF (n.LE.(number_tidal_steps-ntmod)) THEN
      tidalstep(n) = nosteps
   ELSE
      tidalstep(n) = nosteps + 1
   ENDIF
ENDDO n_110

!---updqate time for sediments and morphology
accelstep = 0

!
!2.Time loop
!-----------
!

nperiod_200: DO nperiod=1,noperiods

!   ---initialise
   tidal_counter = 1
   accelstep = accelstep + tidalstep(tidal_counter)

   IF (loglev1.GT.0) THEN
      WRITE (cval,'(I12)') nperiod; cval = ADJUSTL(cval)
      WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//&
                        & 'Entering tidal averaging cycle: '//TRIM(cval)
   ENDIF

   period_start = (nperiod-1)*(nstep_hydro*morph_steps) + 1
   period_end = (nperiod-1)*(nstep_hydro*morph_steps) + nstep_hydro

!
!2.1 Hydrodynamical steps
!------------------------
!

   nt_210: DO nt=period_start,period_end
   
!
!2.1.1 Restart log file
!----------------------
!
         
      CALL monitor_files
         
!
!2.1.2 Date/time
!-------------
!
         
      CALL update_time
         
!
!2.1.3 Define new initial conditions
!-----------------------------------
!

      IF (physinit) CALL define_phsics

!
!2.1.4 Update hydrodynamics
!--------------------------
!
!     ---store old water depths
      CALL store_depths_old
!     ---dynamic mask
      IF (iopt_fld.EQ.2) CALL mask_function
      IF (iopt_weibar.EQ.1) CALL weirs_mask
!     ---astronomical tide
      IF (iopt_astro_tide.EQ.1) CALL astronomical_tides
!     ---currents and elevations
      CALL hydrodynamic_equations
         
!
!2.1.5 Update sediment transport and morphology
!----------------------------------------------
!

      IF (nt.EQ.accelstep) CALL morphological_step
      
!
!2.1.6 Model output
!----------------
!
!     ---user-defined
      IF (iopt_out_tsers.EQ.1) CALL time_series
      CALL usrdef_output

!     ---restart files
      CALL write_phsics
      IF (iopt_sed.GT.0) CALL write_sedics
      IF (iopt_morph.GT.0) CALL write_morphics
         
   ENDDO nt_210

!
!2.2 Morphological steps
!-----------------------
!
!  ---initialise at first morphological time step
   morph_start = period_end
   morph_end = nperiod*nstep_hydro*morph_steps
   morph_counter = 0

   nt_220: DO nt=morph_start+1,morph_end

!
!2.2.1 Reinitialise at new morphological cycle
!---------------------------------------------
!

      IF (mult_index(nt-1,nstep_hydro)) THEN
         tidal_counter = 1
         accelstep = accelstep + tidalstep(tidal_counter)
         morph_counter = morph_counter + 1
         IF (loglev1.GT.0) THEN
            WRITE (cval,'(I12)') morph_counter; cval = ADJUSTL(cval)
            WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'morph_counter: '//&
                              & TRIM(cval)
         ENDIF
      ENDIF

!         
!2.2.2 Restart log file
!----------------------
!
 
      CALL monitor_files

!
!2.2.3 Date/time
!---------------         
!    

      CALL update_time
        
!
!2.2.4 Update sediment transport and morphology
!----------------------------------------------
!

      IF (nt.EQ.accelstep) CALL morphological_step

!
!2.2.5 Model output
!------------------
!
!     ---user-defined
      IF (iopt_out_tsers.EQ.1) CALL time_series
      CALL usrdef_output

!     ---restart files
      CALL write_phsics
      IF (iopt_sed.GT.0) CALL write_sedics
      IF (iopt_morph.GT.0) CALL write_morphics

   ENDDO nt_220

ENDDO nperiod_200

CALL log_timer_out(npcc,itm_morph)


RETURN

END SUBROUTINE morphology_accel

!==============================================================================

SUBROUTINE morphology_equation
!******************************************************************************
!
! *morphology_equation* Main program unit of the morphological model
!
! Author - Alexander Breugem and Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.10
!
! Description -
!
! Reference -
!
! Calling program - coherens_main, initialise_model
!
! External calls - adjust_layer_top, adjust_multi_layers, adjust_one_layer,
!                  avalanching, bed_elevation_equation, depth_at_nodes,
!                  calculate_active_layer_thicknes, mass_balance, water_depths
!
! Module calls - exchange_mod, mult_index
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE morpharrays
USE morphpars
USE morphswitches
USE sedarrays
USE sedids
USE sedpars
USE switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out
USE paral_comms, ONLY: exchange_mod
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: f, npcc
INTEGER, DIMENSION(4) :: nhexch


morphstep = (nt/icmorph)*icmorph.EQ.nt
IF (.NOT.morphstep) RETURN

procname(pglev+1) = 'morphology_equation'
CALL log_timer_in(npcc)

!
!1. Initialise parameters for mass balance
!----------------------------------------
!

IF (nt.EQ.0) THEN

!
!1.1 Initial mass balance
!------------------------
!

   IF (icsedbal.GT.1) CALL mass_balance

!
!1.2 Sediment availablity
!------------------------
!

   IF (iopt_morph_fixed_layer.EQ.1) THEN

!     ---per fraction
      f_120: DO f=1,nf
         WHERE (maskatc_int)
            sed_avail(:,:,f) = SUM(bed_layer_thickness*&
                                 & bed_fraction(1:ncloc,1:nrloc,:,f),DIM=3)
         END WHERE
      ENDDO f_120

!     ---total
      WHERE (maskatc_int)
         sed_avail_tot = SUM(bed_layer_thickness,DIM=3)
      END WHERE

   ENDIF

!
!1.3 Adjust top layer to active layer
!------------------------------------
!

   IF (nb.GT.1.AND.iopt_morph_active_layer.GT.0) CALL adjust_layer_top

   GOTO 1000

ENDIF

!
!2. Bed update
!-------------
!

CALL bed_elevation_equation

!
!3. Adjust bed parameters
!------------------------
!
!3.1 One layer model
!-------------------
!

IF (nb.EQ.1) THEN

   IF (iopt_morph_fixed_layer.EQ.1) THEN
      CALL adjust_one_layer
   ELSE
      WHERE (maskatc_int)
         bed_update = SUM(bed_update_dep+bed_update_ero,DIM=3)
      ELSEWHERE
         bed_update = 0.0
      END WHERE
   ENDIF

!
!3.2 Multi-layer model
!---------------------
!

ELSE

   CALL adjust_multi_layers

ENDIF

!
!4. Adjust bathymetry
!--------------------
!
!---bathymetry
WHERE (maskatc_int)
   depmeanatc(1:ncloc,1:nrloc) = depmeanatc(1:ncloc,1:nrloc) - bed_update
END WHERE

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = nhalo
   CALL exchange_mod(depmeanatc,(/1-nhalo,1-nhalo/),nhexch,iarr_depmeanatc)
ENDIF

!
!5. Bed fractions
!----------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = (/1,0,1,0/)
   CALL exchange_mod(bed_fraction,(/0,0,1,1/),nhexch,iarr_bed_fraction)
ENDIF

!
!6. Avalanching
!--------------
!

IF (iopt_morph_avalanching.GT.0.AND.mult_index(nt,icaval)) THEN
   CALL avalanching
ENDIF

!
!7. Water depths
!---------------
!
!---bathymetry at other nodes
CALL depth_at_nodes

!---water depths
CALL water_depths

!
!8. Median particle diameter
!---------------------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = (/1,0,1,0/)
   CALL exchange_mod(d50_bed,(/0,0/),nhexch,iarr_d50_bed)
ENDIF

!
!9. Accumulated bed update
!-------------------------
!

WHERE (maskatc_int)
   bed_update_int = bed_update_int + bed_update
END WHERE

!
!10. Mass balance
!----------------
!

IF (icsedbal.GT.1) CALL mass_balance

1000 CALL log_timer_out(npcc,itm_morph)


RETURN

END SUBROUTINE morphology_equation

!========================================================================

SUBROUTINE morphological_step 
!************************************************************************
!
! *morphological_step2* update sediment transport and morphodynamics at the
!                       appropriate time steps when tidal acceleration is
!                       activated
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.10
!
! Description -
!
! Reference -
!
! Calling program - morphology_accel
!
! External calls - morphology_equation, sediment_equation,
!                  update_hydrodynamic_estimation
!
! Module calls - log_timer_in, log_timer_out, mult_index
!
!************************************************************************
!
USE currents
USE depths
USE gridpars
USE iopars
USE morpharrays
USE morphpars
USE switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE


procname(pglev+1) = 'morphological_step'
CALL log_timer_in()

!
!1. Hydrodynamical step
!----------------------
!

IF (nt.GT.period_start.AND.nt.LE.period_end) THEN

!  ---save hydrodynamical data
   umvel_guess(:,:,tidal_counter) = umvel(1:ncloc+1,1:nrloc)
   vmvel_guess(:,:,tidal_counter) = vmvel(1:ncloc,1:nrloc+1)
   zeta_guess(:,:,tidal_counter) = zeta
   depth_guess(:,:,tidal_counter) = deptotatc(1:ncloc,1:nrloc)

!  ---update sediments
   IF (iopt_sed.GT.0) CALL sediment_equation

!   ---update morphology
   IF (iopt_morph.GT.0) CALL morphology_equation

!  ---update counters for hydrodynamics
   IF (nt.LT.period_end) THEN
      accelstep = accelstep + tidalstep(tidal_counter)
      tidal_counter = tidal_counter + 1
   ENDIF

!
!2. Morphological steps
!----------------------
!

ELSE

!  ---update hydrodynamics  
   CALL update_hydrodynamic_estimation

!  ---update sediments
   IF (iopt_sed.GT.0) CALL sediment_equation

!  ---update morphology
   IF (iopt_morph.GT.0) CALL morphology_equation

!  ---update couters for hydrodynamics
   IF (.NOT.mult_index(nt,nstep_hydro)) THEN
      accelstep = accelstep + tidalstep(tidal_counter)
      tidal_counter = tidal_counter + 1
   ENDIF

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE morphological_step

!========================================================================

SUBROUTINE update_hydrodynamic_estimation
!************************************************************************
!
! *update_hydrodynamic_estimation* upload and correct hydrodynamical data
!                                  during a morphological cycle in case
!                                  morphological acceleration is activated
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Morphology_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - morphological_step
!
! External calls - bottom_stress, water_depths
!
! Module calls - error_alloc, exchange_mod, Carr_at_U, Carr_at_V, Uarr_at_C,
!                Varr_at_C
!
!************************************************************************
!
USE currents
USE depths
USE fluxes
USE grid
USE gridpars
USE iopars
USE modids
USE morpharrays
USE morphpars
USE morphswitches
USE switches
USE syspars
USE array_interp, ONLY:  Carr_at_U, Carr_at_V, Uarr_at_C, Varr_at_C
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: n
INTEGER, DIMENSION(4) :: nhexch
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: umvelatc, velcorr, vmvelatc


procname(pglev+1) = 'update_hydrodynamic_estimation'
CALL log_timer_in()

n = tidal_counter

!
!1. Initialise arrays
!--------------------
!
!---allocate
ALLOCATE (velcorr(ncloc,nrloc),STAT=errstat)
CALL error_alloc('velcorr',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (umvelatc(0:ncloc,nrloc),STAT=errstat)
CALL error_alloc('umvelatc',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (vmvelatc(ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('vmvelatc',2,(/ncloc,nrloc+1/),kndrtype)

!---initialise
velcorr = 0.0; umvelatc = 0.0; vmvelatc = 0.0

!
!2. Without velocity correction
!------------------------------
!
!---elevations
zeta = zeta_guess(:,:,n)

!---currents
IF (iopt_morph_tidal_scheme.EQ.0) THEN
   umvel(1:ncloc+1,1:nrloc) = umvel_guess(:,:,n)
   vmvel(1:ncloc,1:nrloc+1) = vmvel_guess(:,:,n)
ENDIF

!
!3. Apply velocity correction
!----------------------------
!

IF (iopt_morph_tidal_scheme.GT.0) THEN

!
!3.1 Velocities at C-nodes
!--------------------------
!

   CALL Uarr_at_C(umvel_guess(:,:,n),umvelatc(1:ncloc,:),1,1,&
               & (/1,1,1/),(/ncloc+1,nrloc,1/),1,iarr_umvel,.FALSE.)
   CALL Varr_at_C(vmvel_guess(:,:,n),vmvelatc(:,1:nrloc),1,1,&
               & (/1,1,1/),(/ncloc,nrloc+1,1/),1,iarr_vmvel,.FALSE.)

!
!3.2 Correction factor
!---------------------
!

   SELECT CASE (iopt_morph_tidal_scheme)

      CASE (1)
         WHERE (maskatc_int)
            velcorr = depth_guess(:,:,n)/deptotatc(1:ncloc,1:nrloc)
         END WHERE

      CASE (2)   
         WHERE (maskatc_int)
            velcorr = SQRT(depth_guess(:,:,n)/deptotatc(1:ncloc,1:nrloc))*&
              & (1.0+LOG(depth_guess(:,:,n)/zroughatc(1:ncloc,1:nrloc)))/&
              & (1.0+LOG(deptotatc(1:ncloc,1:nrloc)/zroughatc(1:ncloc,1:nrloc)))
         END WHERE

      CASE (3)
         WHERE (maskatc_int)
            velcorr = SQRT(depth_guess(:,:,n)/deptotatc(1:ncloc,1:nrloc))*&
              & (1.0+LOG(depth_guess(:,:,n)/zroughatc(1:ncloc,1:nrloc)))/&
              & (1.0+LOG(deptotatc(1:ncloc,1:nrloc)/zroughatc(1:ncloc,1:nrloc)))
            velcorr = SQRT(velcorr)
         END WHERE

   END SELECT

!
!3.3 Depth-mean current at C-nodes
!---------------------------------
!

   WHERE (maskatc_int)
      umvelatc(1:ncloc,:) = velcorr*umvelatc(1:ncloc,:)
      vmvelatc(:,1:nrloc) = velcorr*vmvelatc(:,1:nrloc)
   END WHERE

!
!3.4 Exchange halo sections
!--------------------------
!

   IF (iopt_MPI.EQ.1) THEN
      nhexch = (/1,0,0,0/)
      CALL exchange_mod(umvelatc,(/0,1/),nhexch,0)
      nhexch = (/0,0,1,0/)
      CALL exchange_mod(vmvelatc,(/1,0/),nhexch,0)
   ENDIF

!
!3.5. Interpolate at velocity nodes
!---------------------------------
!

   CALL Carr_at_U(umvelatc,umvel(1:ncloc,1:nrloc),1,1,(/0,1,nz/),&
               & (/ncloc,nrloc,nz/),1,iarr_umvel,.FALSE.)
   CALL Carr_at_V(vmvelatc,vmvel(1:ncloc,1:nrloc),1,1,(/1,0,nz/),&
               & (/ncloc,nrloc,nz/),1,iarr_umvel,.FALSE.)

ENDIF

!
!4. Exchange halo sections
!--------------------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = (/nhalo,nhalo,1,1/)
   CALL exchange_mod(umvel,(/1-nhalo,0/),nhexch,iarr_umvel)
   nhexch = (/1,1,nhalo,nhalo/)
   CALL exchange_mod(vmvel,(/0,1-nhalo/),nhexch,iarr_vmvel)
ENDIF

!
!4. Update physical variables
!----------------------------
!
!---water depths
CALL water_depths

!---bottom stress
IF (iopt_bstres_form.GT.0) CALL bottom_stress

!
!5. Deallocate arrays
!--------------------
!

DEALLOCATE (velcorr,umvelatc,vmvelatc)

CALL log_timer_out()


RETURN

END SUBROUTINE update_hydrodynamic_estimation
