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

!*****************************************************************
!
! *Inundation_Schemes* Routines ror applying inundation schemes
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Inundation_Schemes.f90 V2.11
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description -
!
! Reference -
!
! Routines - drying_factor, mask_function, minimum_depth
!            
!
!*****************************************************************
!

!=================================================================

SUBROUTINE mask_function
!*****************************************************************
!
! *mask_function* Apply dynamic mask
!
! Author - ANTEA
!
! Version - @(COHERENS)Inundation_Schemes.f90 V2.11
!
! Description - apply mask criterion to set wet cells to dry, or dry cells
!               to wet
!             - reset pointer arrays
!             - set currents to zero at blocked velocity nodes
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - thin_dams, update_pointer_arrays
!
! Module calls - error_alloc, exchange_mod
!
!*****************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
LOGICAL :: drycell
INTEGER :: ii, i, j, jj, k, l
LOGICAL, DIMENSION(nofldmasks) :: flmask
INTEGER, DIMENSION(2) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
REAL :: sum1, sum2
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: weight


procname(pglev+1) = 'mask_function'
CALL log_timer_in()

!
!1. Initialise open boundary conditions on first call
!----------------------------------------------------
!

IF (nt.EQ.1) THEN
  ALLOCATE (weight(0:ncloc+1,0:nrloc+1),STAT=errstat)
  CALL error_alloc('weight',2,(/ncloc+2,nrloc+2/),kndrtype)
  weight = MERGE(1.0,0.0,seapoint)
ENDIF

!
!2. (Re)set cell pointer at C-nodes
!----------------------------------
!

j_200: DO j=1,nrloc
i_200: DO i=1,ncloc

!
!2.1 Evaluate mask criteria
!--------------------------
!

   l_211: DO l=1,nofldmasks
      IF (fld_mask(l).NE.0.AND.weight(i,j).EQ.1.0) THEN
         SELECT CASE (l)
            CASE (1)
               flmask(1) = MAX(deptotatc(i-1,j),deptotatc(i+1,j),&
                             & deptotatc(i,j-1),deptotatc(i,j+1),&
                             & deptotatc(i,j)).LT.dthd_fld
            CASE (2)
               flmask(2) = MIN(deptotatc(i-1,j),deptotatc(i+1,j),&
                             & deptotatc(i,j-1), deptotatc(i,j+1),&
                             & deptotatc(i,j)).LT.dthd_fld
            CASE(3)
               sum1 = weight(i-1,j)*deptotatc(i-1,j) + &
                    & weight(i+1,j)*deptotatc(i+1,j) + &
                    & weight(i,j-1)*deptotatc(i,j-1) + &
                    & weight(i,j+1)*deptotatc(i,j+1) + &
                    & weight(i,j)*deptotatc(i,j)
               sum2 = weight(i-1,j) + weight(i+1,j) + weight(i,j-1) + &
                    & weight(i,j+1) + weight(i,j)
               flmask(3) = sum1/sum2.LT.dthd_fld
            CASE (4)
               flmask(4) = MAX(deptotatc(i-1,j),deptotatc(i+1,j),&
                             & deptotatc(i,j-1), deptotatc(i,j+1)).LT.dthd_fld
            CASE (5)
               flmask(5) = MIN(deptotatc(i-1,j),deptotatc(i+1,j),&
                             & deptotatc(i,j-1), deptotatc(i,j+1)).LT.dthd_fld
            CASE (6)
               sum1 = weight(i-1,j)*deptotatc(i-1,j) + &
                    & weight(i+1,j)*deptotatc(i+1,j) + &
                    & weight(i,j-1)*deptotatc(i,j-1) + &
                    & weight(i,j+1)*deptotatc(i,j+1)
               sum2 = weight(i-1,j) + weight(i+1,j) + &
                    & weight(i,j-1) + weight(i,j+1)
               IF (sum2.GT.0.0) flmask(6) = sum1/sum2.LT.dthd_fld
            CASE (7)
               flmask(7) = MAX(nodeatc(i-1,j),nodeatc(i+1,j),nodeatc(i,j-1),&
                             & nodeatc(i,j+1)).EQ.0
            CASE (8)
               flmask(8) = MIN(nodeatc(i-1,j),nodeatc(i+1,j),nodeatc(i,j-1),&
                             & nodeatc(i,j+1)).EQ.0
            CASE (9)
               flmask(9) = MAX(nodeatc(i-1,j),nodeatc(i+1,j),nodeatc(i,j-1), &
                             & nodeatc(i,j+1)).EQ.0.AND.&
                           & deptotatc(i,j).LT.dthd_fld
            CASE (10)
               flmask(10) = MIN(nodeatc(i-1,j),nodeatc(i+1,j),&
                               & nodeatc(i,j-1),nodeatc(i,j+1)).EQ.0.AND.&
                           & deptotatc(i,j).LT.dthd_fld
            CASE (11)
               IF (.NOT.(seapoint(i,j-1).OR.seapoint(i,j+1))) THEN
                  flmask(11) = MIN((depmeanatc(i-1,j)-deptotatc(i-1,j)),&
                       & (depmeanatc(i+1,j)-deptotatc(i+1,j))).GT.&
                       & depmeanatc(i,j)
               ELSEIF (.NOT.(seapoint(i-1,j).OR.seapoint(i+1,j))) THEN
                  flmask(11) = MIN((depmeanatc(i,j-1)-deptotatc(i,j-1)),&
                       & (depmeanatc(i,j+1)-deptotatc(i,j+1))).GT.&
                       & depmeanatc(i,j)
               ENDIF

         END SELECT

      ELSE
         flmask(l) = .FALSE.
      ENDIF
      
   ENDDO l_211

!
!2.2 Apply criteria at C-nodes
!-----------------------------
!

   drycell=.FALSE.
   l_220: DO l=1,nofldmasks
      drycell = drycell.OR.flmask(l)
   ENDDO l_220
   nodeatc(i,j) = MERGE(0,1,weight(i,j).EQ.0.0.OR.drycell)
   maskatc_int(i,j) = nodeatc(i,j).EQ.1
   rmaskatc(i,j) = nodeatc(i,j)

ENDDO i_200
ENDDO j_200

!
!2.3 Exchange halo sections
!--------------------------
!

IF (iopt_MPI.EQ.1) THEN
   lbounds(1:2) = 1-nhalo; nhexch = nhalo
   CALL exchange_mod(nodeatc,lbounds,nhexch,iarr_nodeatc)
ENDIF

!
!3. (Re)set pointer arrays
!-------------------------
!
!3.1 U-nodes
!-----------
!
!---interior faces
j_311: DO j=1,nrloc
i_311: DO i=1,ncloc
   WHERE (nodeatu(i,j,:).LE.2)
      nodeatu(i,j,:) = COUNT(nodeatc(i-1:i,j).GT.0)
   END WHERE
ENDDO i_311
ENDDO j_311

!---open sea boundaries
ii_312: DO ii=1,nobuloc_ext
   i = iobuloc(ii); j = jobuloc(ii)
   IF (ALL(nodeatc(i-1:i,j).EQ.0)) THEN
      nodeatu(i,j,:) = 1
   ELSE
      nodeatu(i,j,:) = MERGE(3,4,ii.LE.nosbu)
   ENDIF
ENDDO ii_312
   
!---thin dams
IF (iopt_thndam.EQ.1) CALL thin_dams('U')

!
!
!3.2 V-nodes
!-----------
!
!---interior faces
j_321: DO j=1,nrloc
i_321: DO i=1,ncloc
   WHERE (nodeatv(i,j,:).LE.2)
      nodeatv(i,j,:) = COUNT(nodeatc(i,j-1:j).GT.0)
   END WHERE
ENDDO i_321
ENDDO j_321

!---open sea boundaries
jj_322: DO jj=1,nobvloc_ext
   i = iobvloc(jj); j = jobvloc(jj)
   IF (ALL(nodeatc(i,j-1:j).EQ.0)) THEN
      nodeatv(i,j,:) = 1
   ELSE
      nodeatv(i,j,:) = MERGE(3,4,jj.LE.nosbv)
   ENDIF
ENDDO jj_322

!---thin dams
IF (iopt_thndam.EQ.1) CALL thin_dams('V')

!
!3.3 Other nodes/reset currents
!------------------------------
!

CALL update_pointer_arrays

!
!4. Block vertical currents at C-nodes
!-------------------------------------
!

IF (iopt_grid_nodim.EQ.3) THEN
   k_410: DO k=1,nz
      WHERE (nodeatc(0:ncloc,0:nrloc).EQ.0)
         wvel(:,:,k) = 0.0
      END WHERE
      WHERE (maskatc_int)
         wphys(:,:,k) = 0.0
      END WHERE
   ENDDO k_410
ENDIF
   
!
!5. Deallocate arrays
!--------------------
!

IF (nt+ic3d.GT.nstep) DEALLOCATE (weight)

CALL log_timer_out()


RETURN

END SUBROUTINE mask_function

!========================================================================

SUBROUTINE minimum_depths
!*****************************************************************
!
! *minimum_depths* Impose minimum depth to avoid negative depths
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Inundation_Schemes.f90 V2.11
!
! Description -
!
! Reference -
!
! Calling program - water_depths
!
! External calls - thin_dams, update_pointer_arrays   
!
! Module calls - exchange_mod
!  
!*****************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE physpars
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'minimum_depths'
CALL log_timer_in()

!
!1. Update error in water depth by imposing minimum depth
!---------------------------------------------------------
!

IF (nt.GT.0) THEN
   WHERE (maskatc_int.AND.deptotatc_prev.GT.dmin_fld.AND.&
        & deptotatc(1:ncloc,1:nrloc).LT.dmin_fld)
      deptotatc_err = deptotatc_err - deptotatc(1:ncloc,1:nrloc) + dmin_fld
   END WHERE
ENDIF

!
!2. Minimùum depths
!------------------
!
!2.1 C-nodes
!-----------
!
   
WHERE (seapoint)
   deptotatc = MAX(deptotatc,dmin_fld)
END WHERE

!
!2.2 U-nodes
!-----------
!

WHERE (node2du(1:ncloc,1:nrloc).LE.2.AND.&
    & (deptotatc(0:ncloc-1,1:nrloc).LT.dcrit_fld.OR.&
     & deptotatc(1:ncloc,1:nrloc).LT.dcrit_fld))
   deptotatu(1:ncloc,1:nrloc) = MIN(deptotatc(0:ncloc-1,1:nrloc),&
                                  & deptotatc(1:ncloc,1:nrloc))
END WHERE

WHERE (node2du(1:ncloc,1:nrloc).GT.2)
   deptotatu(1:ncloc,1:nrloc) = MAX(deptotatu(1:ncloc,1:nrloc),dmin_fld)
END WHERE

!
!2.3 V-nodes
!----------
!

WHERE (node2dv(1:ncloc,1:nrloc).LE.2.AND.&
    & (deptotatc(1:ncloc,0:nrloc-1).LT.dcrit_fld.OR.&
    &  deptotatc(1:ncloc,1:nrloc).LT.dcrit_fld))
   deptotatv(1:ncloc,1:nrloc) = MIN(deptotatc(1:ncloc,0:nrloc-1),&
                                  & deptotatc(1:ncloc,1:nrloc))
END WHERE

WHERE (node2dv(1:ncloc,1:nrloc).GT.2)
   deptotatv(1:ncloc,1:nrloc) = MAX(deptotatv(1:ncloc,1:nrloc),dmin_fld)
END WHERE

!
!2.4 UV-nodes
!-----------
!

WHERE (node2duv(1:ncloc+1,1:nrloc+1).GT.0)
   deptotatuv = MAX(deptotatuv,dmin_fld)
END WHERE

CALL log_timer_out()


RETURN

END SUBROUTINE minimum_depths

!========================================================================

SUBROUTINE drying_factor
!*****************************************************************
!
! *drying_factor* Update the drying ("alpha") factor
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Inundation_Schemes.f90 V2.11
!
! Description -
!
! Reference -
!
! Calling program - hydrodynamic_equations, initialise_model
!
!*****************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE physpars
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
REAL :: depdif


procname(pglev+1) = 'drying_factor'
CALL log_timer_in()


depdif = dcrit_fld - dmin_fld

!
!1. Power law
!-------------   
!

IF (iopt_fld_alpha.EQ.1) THEN
   
!  ---C-nodes
   WHERE (maskatc_int)
      alphatc_fld = MAX(MIN(1.0,&
                  & ((deptotatc(1:ncloc,1:nrloc)-dmin_fld)/depdif)**ndwexp),0.0)
   END WHERE

!   ---U-nodes
   WHERE (node2du(1:ncloc,1:nrloc).GT.1)
      alphatu_fld = MAX(MIN(1.0,&
                  & ((deptotatu(1:ncloc,1:nrloc)-dmin_fld)/depdif)**ndwexp),0.0)
   END WHERE

!   ---V-nodes
   WHERE (node2dv(1:ncloc,1:nrloc).GT.1)
      alphatv_fld = MAX(MIN(1.0,&
                  & ((deptotatv(1:ncloc,1:nrloc)-dmin_fld)/depdif)**ndwexp),0.0)
   END WHERE

!
!2. Advection disabled below critical depth
!------------------------------------------
!   

ELSEIF (iopt_fld_alpha.EQ.2) THEN

!  ---C-nodes
   WHERE (maskatc_int)
      alphatc_fld = MERGE(0.0,1.0,deptotatc(1:ncloc,1:nrloc).LE.dcrit_fld)
   END WHERE
!  ---U-nodes
   WHERE (node2du(1:ncloc,1:nrloc).GT.1)
      alphatu_fld = MERGE(0.0,1.0,deptotatu(1:ncloc,1:nrloc).LE.dcrit_fld)
   END WHERE
!  ---V-nodes
   WHERE (node2dv(1:ncloc,1:nrloc).GT.1)
      alphatv_fld = MERGE(0.0,1.0,deptotatv(1:ncloc,1:nrloc).LE.dcrit_fld)
   END WHERE

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE drying_factor
