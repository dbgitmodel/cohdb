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

MODULE array_interp
!************************************************************************
!
! *array_interp* Array and scalar interpolation routines on the model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - optimisation for 2-D masks implemented by ANTEA
!
! Subroutines - Carr_at_U, Carr_at_UV, Carr_at_UW, Carr_at_V, Carr_at_VW,
!               Carr_at_W, Uarr_at_C, Uarr_at_UV, Uarr_at_UW, Uarr_at_V,
!               Uarr_at_W, UVarr_at_C, UVarr_at_U, UVarr_at_V, UWarr_at_U,
!               UWarr_at_W, Varr_at_C, Varr_at_U, Varr_at_UV, Varr_at_VW,
!               Varr_at_W, VWarr_at_V, VWarr_at_W, Warr_at_C, Warr_at_U,
!               Warr_at_UW, Warr_at_V, Warr_at_VW
!
! Functions - Cvar_at_U, Cvar_at_UV, Cvar_at_UW, Cvar_at_V, Cvar_at_VW,
!             Cvar_at_W, Uvar_at_C, Uvar_at_UV, Uvar_at_UW, Uvar_at_V,
!             Uvar_at_W, UVvar_at_C, UVvar_at_U, UVvar_at_V, UWvar_at_U,
!             UWvar_at_W, Vvar_at_C, Vvar_at_U, Vvar_at_UV, Vvar_at_VW,
!             Vvar_at_W, VWvar_at_V, VWvar_at_W, Wvar_at_C, Wvar_at_U,
!             Wvar_at_UW, Wvar_at_V, Wvar_at_VW
!
! Module calls - error_alloc
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: fdest, hreg, vreg
INTEGER :: l, l1, l2, l3, npcc, u1, u2, u3
REAL :: flagval, wsum
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: fdest2
LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: fdest3
REAL, ALLOCATABLE, DIMENSION(:,:) :: sflag2, weights2, wsum2
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: sflag3, weights3, wsum3

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*hreg*      LOGICAL .TRUE. for area averaging in the horizontal
!*vreg*      LOGICAL .TRUE. for area averaging in the vertical
!*l1*        INTEGER Lower bound of input array in first dimension
!*l2*        INTEGER Lower bound of input array in second dimension
!*l3*        INTEGER Lower bound of input array in third dimension
!*u1*        INTEGER Upper bound of input array in first dimension
!*u2*        INTEGER Upper bound of input array in second dimension
!*u3*        INTEGER Upper bound of input array in third dimension
!*flagval*   REAL    Replacement value at dry (or invalid) points
!*fdest2*    REAL    Data flag for elimination of land points or cell faces at
!                    destination center node
!*fdest3*    REAL    Data flag for elimination of land points or cell faces at
!                    destination velocity node
!*sflag2*    REAL    Data flag for elimination of land points or cell faces at
!                    source center node
!*sflag3*    REAL    Data flag for elimination of land points or cell faces at
!                    source velocity node
!                    Values are set to 1 or 0 depending on value of intnode.
!*weights2*  REAL    Non-uniform weight factors in nominator (2-D case)
!*weights3*  REAL    Non-uniform weight factors in nominator (3-D case)
!*wsum2*     REAL    Sum of non-uniform weight factors in denominator
!                    (2-D case)
!*wsum3*     REAL    Sum of non-uniform weight factors in denominator
!                    (3-D case)
!
!------------------------------------------------------------------------------
!

CONTAINS

!========================================================================

SUBROUTINE Carr_at_U(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                   & info,mask,outflag,hregular,angle)
!************************************************************************
!
! *Carr_at_U* Interpolate an array from C- to U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.11.1
!
! Description -
!
!************************************************************************

USE grid
USE gridpars
!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: angle, hregular
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(lbounds(1):ubounds(1),&
                                       & lbounds(2):ubounds(2)) :: mask
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1)+1:ubounds(1),lbounds(2):ubounds(2),&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at C-nodes
!*xout*      REAL    Output array at U-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => wet points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array
!*ubounds*   INTEGER Upper bounds of input array
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*mask*      LOGICAL Mask array at source (C-) nodes
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!*angle*     LOGICAL Flag to select whether the input data have an angle
!                    dimension (in radians) 
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: fangle
INTEGER :: k


procname(pglev+1) = 'Carr_at_U'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXregX
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

IF (.NOT.PRESENT(angle)) THEN
   fangle = .FALSE.
ELSE
   fangle = angle
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (fdest2(l1+1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('fdest2',2,(/u1-l1,u2-l2+1/),kndlog)
ELSE
   ALLOCATE (fdest3(l1+1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('fdest3',3,(/u1-l1,u2-l2+1,u3-l3+1/),kndlog)
ENDIF

!---source
SELECT CASE (intsrce)
   CASE (0); sflag2 = 1.0
   CASE (1)
      IF (PRESENT(mask)) THEN
         sflag2 = MERGE(1.0,0.0,mask.AND.nodeatc(l1:u1,l2:u2).GT.0)
      ELSE
         sflag2 = MERGE(1.0,0.0,nodeatc(l1:u1,l2:u2).GT.0)
      ENDIF
END SELECT

!---destination
IF (iopt_arrint_3D.EQ.0) THEN
   SELECT CASE (intdest)
      CASE (0); fdest2 = .TRUE.
      CASE (1); fdest2 = node2du(l1+1:u1,l2:u2).GT.0
      CASE (2); fdest2 = node2du(l1+1:u1,l2:u2).EQ.2
      CASE (3); fdest2 = node2du(l1+1:u1,l2:u2).GT.1
      CASE (4); fdest2 = node2du(l1+1:u1,l2:u2).EQ.1.OR.&
                       & node2du(l1+1:u1,l2:u2).EQ.2
   END SELECT
ELSE
   SELECT CASE (intdest)
      CASE (0); fdest3 = .TRUE.
      CASE (1); fdest3 = nodeatu(l1+1:u1,l2:u2,l3:u3).GT.0
      CASE (2); fdest3 = nodeatu(l1+1:u1,l2:u2,l3:u3).EQ.2
      CASE (3); fdest3 = nodeatu(l1+1:u1,l2:u2,l3:u3).GT.1
      CASE (4); fdest3 = nodeatu(l1+1:u1,l2:u2,l3:u3).EQ.1.OR.&
                       & nodeatu(l1+1:u1,l2:u2,l3:u3).EQ.2
   END SELECT
ENDIF

!
!4. Weight factors
!-----------------
!
!---allocate
IF (.NOT.hreg) THEN
   ALLOCATE (weights2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('weights2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ENDIF
ALLOCATE (wsum2(l1+1:u1,l2:u2),STAT=errstat)
CALL error_alloc('wsum2',2,(/u1-l1,u2-l2+1/),kndrtype)

!---define
IF (hreg) THEN
   wsum2 = sflag2(l1+1:u1,:) + sflag2(l1:u1-1,:)
ELSE
   weights2 = delxatc(l1:u1,l2:u2)
   wsum2 = sflag2(l1+1:u1,:)*weights2(l1:u1-1,:) + &
         & sflag2(l1:u1-1,:)*weights2(l1+1:u1,:)
ENDIF

IF (iopt_arrint_3D.EQ.0) THEN
   WHERE (wsum2.LE.0.0) fdest2 = .FALSE.
ELSE
   k_410: DO k=l3,u3
      WHERE (wsum2.LE.0.0) fdest3(:,:,k) = .FALSE.
   END DO k_410
ENDIF

!
!5. Interpolate
!--------------
!
!5.1 Using 2-D mask
!------------------
!

IF (iopt_arrint_3D.EQ.0) THEN

!
!5.1.1 Without weights
!---------------------
!

   IF (hreg) THEN
      l_511: DO l=1,nosize
      k_511: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = (sflag2(l1:u1-1,:)*xin(l1:u1-1,:,k,l)+&
                             sflag2(l1+1:u1,:)*xin(l1+1:u1,:,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_511
      ENDDO l_511
      
!
!5.1.2 With weights
!------------------
!

   ELSE
      l_512: DO l=1,nosize
      k_512: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = &
               &(sflag2(l1:u1-1,:)*weights2(l1+1:u1,:)*xin(l1:u1-1,:,k,l)+&
               & sflag2(l1+1:u1,:)*weights2(l1:u1-1,:)*xin(l1+1:u1,:,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_512
      ENDDO l_512
    
   ENDIF

!
!5.2 Using 3-D mask
!------------------
!

ELSE

!
!5.2.1 Without weights
!---------------------
!

   IF (hreg) THEN
      l_521: DO l=1,nosize
      k_521: DO k=l3,u3
         WHERE (fdest3(:,:,k))
            xout(:,:,k,l) = (sflag2(l1:u1-1,:)*xin(l1:u1-1,:,k,l)+&
                             sflag2(l1+1:u1,:)*xin(l1+1:u1,:,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_521
      ENDDO l_521

!
!5.2.2 With weights
!------------------
!

   ELSE

      l_522: DO l=1,nosize
      k_522: DO k=l3,u3
         WHERE (fdest3(:,:,k))
            xout(:,:,k,l) = &
               &(sflag2(l1:u1-1,:)*weights2(l1+1:u1,:)*xin(l1:u1-1,:,k,l)+&
               & sflag2(l1+1:u1,:)*weights2(l1:u1-1,:)*xin(l1+1:u1,:,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_522
      ENDDO l_522

   ENDIF

ENDIF

!
!5.3 Correction for variables with angle dimensions
!--------------------------------------------------
!

IF (fangle) THEN

!
!5.3.1 Using 2-D mask
!--------------------
!

   IF (iopt_arrint_3D.EQ.0) THEN
   
      l_531: DO l=1,nosize
      k_531: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = xout(:,:,k,l) + &
                          & pi*(1.0-&
                          & SIGN(1.0,pi-&
                               & ABS(xin(l1+1:u1,:,k,l)-xin(l1:u1-1,:,k,l))))
         END WHERE
      ENDDO k_531
      ENDDO l_531

!
!5.3.2 Using 3-D mask
!--------------------
!

   ELSE
   
      l_532: DO l=1,nosize
      k_532: DO k=l3,u3
         WHERE (fdest3(:,:,k))
            xout(:,:,k,l) = xout(:,:,k,l) + &
                          & pi*(1.0-&
                          & SIGN(1.0,pi-&
                               & ABS(xin(l1+1:u1,:,k,l)-xin(l1:u1-1,:,k,l))))
         END WHERE
      ENDDO k_532
      ENDDO l_532
      
   ENDIF

ENDIF

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag2,wsum2)
IF (iopt_arrint_3D.EQ.0) THEN
   DEALLOCATE (fdest2)
ELSE
   DEALLOCATE (fdest3)
ENDIF
IF (.NOT.hreg) DEALLOCATE (weights2)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Carr_at_U

!========================================================================

SUBROUTINE Carr_at_UV(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                    & info,outflag,hregular)
!************************************************************************
!
! *Carr_at_UV* Interpolate an array from C- to UV-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.6
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1)+1:ubounds(1),lbounds(2)+1:ubounds(2),&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at C-nodes
!*xout*      REAL    Output array at UV-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => wet points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => wet points only
!*lbounds*   INTEGER Lower bounds of input array
!*ubounds*   INTEGER Upper bounds of input array
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Carr_at_UV'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXYreg
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (fdest2(l1+1:u1,l2+1:u2),STAT=errstat)
   CALL error_alloc('fdest2',2,(/u1-l1,u2-l2/),kndlog)
ELSE
   ALLOCATE (fdest3(l1+1:u1,l2+1:u2,l3:u3),STAT=errstat)
   CALL error_alloc('fdest3',3,(/u1-l1,u2-l2,u3-l3+1/),kndlog)
ENDIF

!---source
SELECT CASE (intsrce)
   CASE (0); sflag2 = 1.0
   CASE (1); sflag2 = MERGE(1.0,0.0,nodeatc(l1:u1,l2:u2).GT.0)
END SELECT

!---destination
IF (iopt_arrint_3D.EQ.0) THEN
   SELECT CASE (intdest)
      CASE (0); fdest2 = .TRUE.
      CASE (1); fdest2 = node2duv(l1+1:u1,l2+1:u2).GT.0
    END SELECT
ELSE
   SELECT CASE (intdest)
      CASE (0); fdest3 = .TRUE.
      CASE (1); fdest3 = nodeatuv(l1+1:u1,l2+1:u2,l3:u3).GT.0
   END SELECT
ENDIF

!
!4. Weight factors
!-----------------
!
!---allocate
IF (.NOT.hreg) THEN
   ALLOCATE (weights2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('weights2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ENDIF
ALLOCATE (wsum2(l1+1:u1,l2+1:u2),STAT=errstat)
CALL error_alloc('wsum2',2,(/u1-l1,u2-l2/),kndrtype)

!---define
IF (hreg) THEN
   wsum2 = sflag2(l1+1:u1,l2+1:u2) + sflag2(l1:u1-1,l2+1:u2) + &
         & sflag2(l1+1:u1,l2:u2-1) + sflag2(l1:u1-1,l2:u2-1)
ELSE
   IF (dXregX.AND.dXregY) THEN 
      weights2 = delyatc(l1:u1,l2:u2)
   ELSEIF (dYregX.AND.dYregY) THEN
      weights2 = delxatc(l1:u1,l2:u2)
   ELSE
      weights2 = delxatc(l1:u1,l2:u2)*delyatc(l1:u1,l2:u2)
   ENDIF
   wsum2 = sflag2(l1+1:u1,l2+1:u2)*weights2(l1:u1-1,l2:u2-1) + &
         & sflag2(l1:u1-1,l2+1:u2)*weights2(l1+1:u1,l2:u2-1) + &
         & sflag2(l1+1:u1,l2:u2-1)*weights2(l1:u1-1,l2+1:u2) + &
         & sflag2(l1:u1-1,l2:u2-1)*weights2(l1+1:u1,l2+1:u2)
ENDIF

IF (iopt_arrint_3D.EQ.0) THEN
   WHERE (wsum2.LE.0.0) fdest2(:,:) = .FALSE.
ELSE
   k_410: DO k=l3,u3
      WHERE (wsum2.LE.0.0) fdest3(:,:,k) = .FALSE.
   END DO k_410
ENDIF

!
!5. Interpolate
!--------------
!
!5.1 Using 2-D mask
!------------------
!

IF (iopt_arrint_3D.EQ.0) THEN

!
!5.1.1 Without weights
!---------------------
!

   IF (hreg) THEN
      l_511: DO l=1,nosize
      k_511: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = (sflag2(l1:u1-1,l2:u2-1)*xin(l1:u1-1,l2:u2-1,k,l)+&
                           & sflag2(l1+1:u1,l2:u2-1)*xin(l1+1:u1,l2:u2-1,k,l)+&
                           & sflag2(l1:u1-1,l2+1:u2)*xin(l1:u1-1,l2+1:u2,k,l)+&
                           & sflag2(l1+1:u1,l2+1:u2)*xin(l1+1:u1,l2+1:u2,k,l))&
                           & /wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval            
         END WHERE
      ENDDO k_511
      ENDDO l_511

!
!5.1.2 With weights
!------------------
!

   ELSE
      l_512: DO l=1,nosize
      k_512: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = (sflag2(l1:u1-1,l2:u2-1)*&
                    & weights2(l1+1:u1,l2+1:u2)*xin(l1:u1-1,l2:u2-1,k,l)+&
                           & sflag2(l1+1:u1,l2:u2-1)*&
                    & weights2(l1:u1-1,l2+1:u2)*xin(l1+1:u1,l2:u2-1,k,l)+&
                           & sflag2(l1:u1-1,l2+1:u2)*&
                    & weights2(l1+1:u1,l2:u2-1)*xin(l1:u1-1,l2+1:u2,k,l)+&
                           & sflag2(l1+1:u1,l2+1:u2)*&
                    & weights2(l1:u1-1,l2:u2-1)*xin(l1+1:u1,l2+1:u2,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval            
         END WHERE
      ENDDO k_512
      ENDDO l_512
   ENDIF

!
!5.2 Using 3-D mask
!------------------
!

ELSE

!
!5.2.1 Without weights
!---------------------
!

   IF (hreg) THEN
      l_521: DO l=1,nosize
      k_521: DO k=l3,u3
         WHERE (fdest3(:,:,k))
            xout(:,:,k,l) = (sflag2(l1:u1-1,l2:u2-1)*xin(l1:u1-1,l2:u2-1,k,l)+&
                           & sflag2(l1+1:u1,l2:u2-1)*xin(l1+1:u1,l2:u2-1,k,l)+&
                           & sflag2(l1:u1-1,l2+1:u2)*xin(l1:u1-1,l2+1:u2,k,l)+&
                           & sflag2(l1+1:u1,l2+1:u2)*xin(l1+1:u1,l2+1:u2,k,l))&
                           & /wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval            
         END WHERE
      ENDDO k_521
      ENDDO l_521

!
!5.2.2 With weights
!------------------
!

   ELSE
      l_522: DO l=1,nosize
      k_522: DO k=l3,u3
         WHERE (fdest3(:,:,k))
            xout(:,:,k,l) = (sflag2(l1:u1-1,l2:u2-1)*&
                    & weights2(l1+1:u1,l2+1:u2)*xin(l1:u1-1,l2:u2-1,k,l)+&
                           & sflag2(l1+1:u1,l2:u2-1)*&
                    & weights2(l1:u1-1,l2+1:u2)*xin(l1+1:u1,l2:u2-1,k,l)+&
                           & sflag2(l1:u1-1,l2+1:u2)*&
                    & weights2(l1+1:u1,l2:u2-1)*xin(l1:u1-1,l2+1:u2,k,l)+&
                           & sflag2(l1+1:u1,l2+1:u2)*&
                    & weights2(l1:u1-1,l2:u2-1)*xin(l1+1:u1,l2+1:u2,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval            
         END WHERE
      ENDDO k_522
      ENDDO l_522
   ENDIF

ENDIF

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag2,wsum2)
IF (iopt_arrint_3D.EQ.0) THEN
   DEALLOCATE (fdest2)
ELSE
   DEALLOCATE (fdest3)
ENDIF
IF (.NOT.hreg) DEALLOCATE (weights2)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Carr_at_UV

!========================================================================

SUBROUTINE Carr_at_UW(xin,xout,intdest,lbounds,ubounds,nosize,ivarid,info,&
                    & outflag,hregular,vregular)
!************************************************************************
!
! *Carr_at_UW* Interpolate an array from C- to UW-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: hregular, vregular
INTEGER, INTENT(IN) :: intdest, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1)+1:ubounds(1),lbounds(2):ubounds(2),&
                           & lbounds(3)+1:ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at C-nodes
!*xout*      REAL    Output array at UW-nodes
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!*vregular*  LOGICAL Flag to select uniform or area averaging in the vertical
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: hvreg
INTEGER :: k


procname(pglev+1) = 'Carr_at_UW'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXregX
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(vregular)) THEN
   vreg = iopt_arrint_vreg.EQ.0
ELSE
   vreg = vregular
ENDIF
hvreg = hreg.AND.vreg

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ALLOCATE (fdest3(l1+1:u1,l2:u2,l3+1:u3),STAT=errstat)
CALL error_alloc('fdest3',3,(/u1-l1,u2-l2+1,u3-l3/),kndlog)

!---source
sflag2 = MERGE(1.0,0.0,nodeatc(l1:u1,l2:u2).GT.0)

!---destination
SELECT CASE (intdest)
   CASE (1); fdest3 = nodeatuw(l1+1:u1,l2:u2,l3+1:u3).GT.0
   CASE (2); fdest3 = nodeatuw(l1+1:u1,l2:u2,l3+1:u3).EQ.2
   CASE (3); fdest3 = nodeatuw(l1+1:u1,l2:u2,l3+1:u3).GT.1
   CASE (4); fdest3 = nodeatuw(l1+1:u1,l2:u2,l3+1:u3).EQ.1.OR.&
                    & nodeatuw(l1+1:u1,l2:u2,l3+1:u3).EQ.2
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
IF (hvreg) THEN
   ALLOCATE (wsum2(l1+1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('wsum2',2,(/u1-l1,u2-l2+1/),kndrtype)
ELSE
   ALLOCATE (weights3(l1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('weights3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
   ALLOCATE (wsum3(l1+1:u1,l2:u2,l3+1:u3),STAT=errstat)
   CALL error_alloc('wsum3',3,(/u1-l1,u2-l2+1,u3-l3/),kndrtype)
ENDIF

!---define
IF (hvreg) THEN
   wsum2 = 2.0*(sflag2(l1+1:u1,:)+sflag2(l1:u1-1,:))
   k_411: DO k=l3+1,u3
      WHERE (wsum2.LE.0.0)
         fdest3(:,:,k) = .FALSE.
      END WHERE
   ENDDO k_411
ELSE
   k_412: DO k=l3,u3
      weights3(:,:,k) = 1.0
      IF (.NOT.hreg) THEN
         weights3(:,:,k) = weights3(:,:,k)*delxatc(l1:u1,l2:u2)
      ENDIF
      IF (.NOT.vreg) THEN
         weights3(:,:,k) = weights3(:,:,k)*delzatc(l1:u1,l2:u2,k)
      ENDIF
      WHERE (nodeatu(l1+1:u1,l2:u2,k).GE.2.AND.nodeatc(l1:u1-1,l2:u2).GT.0)
         weights3(l1+1:u1,:,k) = weights3(l1:u1-1,:,k)
      END WHERE
      WHERE (nodeatu(l1+1:u1,l2:u2,k).GE.2.AND.nodeatc(l1+1:u1,l2:u2).GT.0)
         weights3(l1:u1-1,:,k) = weights3(l1+1:u1,:,k)
      END WHERE
   ENDDO k_412
   k_413: DO k=l3+1,u3
      wsum3(:,:,k) =  sflag2(l1+1:u1,:)*&
                    &(weights3(l1:u1-1,:,k-1)+weights3(l1:u1-1,:,k))+ &
                    & sflag2(l1:u1-1,:)*&
                    &(weights3(l1+1:u1,:,k-1)+weights3(l1+1:u1,:,k))
   ENDDO k_413
   WHERE (wsum3.LE.0.0) fdest3 = .FALSE.
ENDIF

!
!5. Interpolate
!--------------
!
!5.1 Without weights
!-------------------
!

IF (hvreg) THEN

   l_511: DO l=1,nosize
   k_511: DO k=l3+1,u3
      WHERE (fdest3(:,:,k))
         xout(:,:,k,l) =(sflag2(l1:u1-1,:)*&
                       &(xin(l1:u1-1,:,k-1,l)+xin(l1:u1-1,:,k,l))+&
                       & sflag2(l1+1:u1,:)*&
                       &(xin(l1+1:u1,:,k-1,l)+xin(l1+1:u1,:,k,l)))/wsum2
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_511
   ENDDO l_511

!
!5.2 With weights
!----------------
!

ELSE

   l_520: DO l=1,nosize
   k_520: DO k=l3+1,u3
      WHERE (fdest3(:,:,k))
         xout(:,:,k,l) = &
               &(sflag2(l1:u1-1,:)*weights3(l1+1:u1,:,k)*xin(l1:u1-1,:,k-1,l)+&
               & sflag2(l1+1:u1,:)*weights3(l1:u1-1,:,k)*xin(l1+1:u1,:,k-1,l)+&
               & sflag2(l1:u1-1,:)*weights3(l1+1:u1,:,k-1)*xin(l1:u1-1,:,k,l)+&
               & sflag2(l1+1:u1,:)*weights3(l1:u1-1,:,k-1)*xin(l1+1:u1,:,k,l))&
               & /wsum3(:,:,k)
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_520
   ENDDO l_520

ENDIF

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag2,fdest3)
IF (hvreg) THEN
   DEALLOCATE (wsum2)
ELSE
   DEALLOCATE (weights3,wsum3)
ENDIF

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Carr_at_UW

!========================================================================

SUBROUTINE Carr_at_V(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                   & info,mask,outflag,hregular,angle)
!************************************************************************
!
! *Carr_at_V* Interpolate an array from C- to V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.11.1
!
! Description -
!
!************************************************************************
!
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: angle, hregular
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(lbounds(1):ubounds(1),&
                                       & lbounds(2):ubounds(2)) :: mask
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2)+1:ubounds(2),&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at C-nodes
!*xout*      REAL    Output array at V-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => wet points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array
!*ubounds*   INTEGER Upper bounds of input array
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*mask*      LOGICAL Mask array at source (C-) nodes
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!*angle*     LOGICAL Flag to select whether the input data have an angle
!                    dimension (in radians)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: fangle
INTEGER :: k


procname(pglev+1) = 'Carr_at_V'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dYregY
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

IF (.NOT.PRESENT(angle)) THEN
   fangle = .FALSE.
ELSE
   fangle = angle
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (fdest2(l1:u1,l2+1:u2),STAT=errstat)
   CALL error_alloc('fdest2',2,(/u1-l1+1,u2-l2/),kndlog)
ELSE
   ALLOCATE (fdest3(l1:u1,l2+1:u2,l3:u3),STAT=errstat)
   CALL error_alloc('fdest3',3,(/u1-l1+1,u2-l2,u3-l3+1/),kndlog)
ENDIF

!---source
SELECT CASE (intsrce)
   CASE (0); sflag2 = 1.0
   CASE (1)
      IF (PRESENT(mask)) THEN
         sflag2 = MERGE(1.0,0.0,mask.AND.nodeatc(l1:u1,l2:u2).GT.0)
      ELSE
         sflag2 = MERGE(1.0,0.0,nodeatc(l1:u1,l2:u2).GT.0)
      ENDIF
END SELECT

!---destination
IF (iopt_arrint_3D.EQ.0) THEN
   SELECT CASE (intdest)
      CASE (0); fdest2 = .TRUE.
      CASE (1); fdest2 = node2dv(l1:u1,l2+1:u2).GT.0
      CASE (2); fdest2 = node2dv(l1:u1,l2+1:u2).EQ.2
      CASE (3); fdest2 = node2dv(l1:u1,l2+1:u2).GT.1
      CASE (4); fdest2 = node2dv(l1:u1,l2+1:u2).EQ.1.OR.&
                       & node2dv(l1:u1,l2+1:u2).EQ.2
   END SELECT
ELSE
   SELECT CASE (intdest)
      CASE (0); fdest3 = .TRUE.
      CASE (1); fdest3 = nodeatv(l1:u1,l2+1:u2,l3:u3).GT.0
      CASE (2); fdest3 = nodeatv(l1:u1,l2+1:u2,l3:u3).EQ.2
      CASE (3); fdest3 = nodeatv(l1:u1,l2+1:u2,l3:u3).GT.1
      CASE (4); fdest3 = nodeatv(l1:u1,l2+1:u2,l3:u3).EQ.1.OR.&
                       & nodeatv(l1:u1,l2+1:u2,l3:u3).EQ.2
   END SELECT
ENDIF

!
!4. Weight factors
!-----------------
!
!---allocate
IF (.NOT.hreg) THEN
   ALLOCATE (weights2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('weights2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ENDIF
ALLOCATE (wsum2(l1:u1,l2+1:u2),STAT=errstat)
CALL error_alloc('wsum2',2,(/u1-l1+1,u2-l2/),kndrtype)

!---define
IF (hreg) THEN
   wsum2 = sflag2(:,l2+1:u2) + sflag2(:,l2:u2-1)
ELSE
   weights2 = delyatc(l1:u1,l2:u2)
   wsum2 = sflag2(:,l2+1:u2)*weights2(:,l2:u2-1) + &
         & sflag2(:,l2:u2-1)*weights2(:,l2+1:u2)
ENDIF

IF (iopt_arrint_3D.EQ.0) THEN
   WHERE (wsum2.LE.0.0) fdest2(:,:) = .FALSE.
ELSE
   k_410: DO k=l3,u3
      WHERE (wsum2.LE.0.0) fdest3(:,:,k) = .FALSE.
   END DO k_410
ENDIF

!
!5. Interpolate
!--------------
!
!5.1 Using 2-D mask
!------------------
!

IF (iopt_arrint_3D.EQ.0) THEN

!
!5.1.1 Without weights
!---------------------
!

   IF (hreg) THEN
      l_511: DO l=1,nosize
      k_511: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = (sflag2(:,l2:u2-1)*xin(:,l2:u2-1,k,l)+&
                             sflag2(:,l2+1:u2)*xin(:,l2+1:u2,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_511
      ENDDO l_511

!
!5.1.2 With weights
!------------------
!

   ELSE
      l_512: DO l=1,nosize
      k_512: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = &
               &(sflag2(:,l2:u2-1)*weights2(:,l2+1:u2)*xin(:,l2:u2-1,k,l)+&
               & sflag2(:,l2+1:u2)*weights2(:,l2:u2-1)*xin(:,l2+1:u2,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_512
      ENDDO l_512
   ENDIF

!
!5.2 Using 3-D mask
!------------------
!

ELSE

!
!5.2.1 Without weights
!---------------------
!

   IF (hreg) THEN
      l_521: DO l=1,nosize
      k_521: DO k=l3,u3
         WHERE (fdest3(:,:,k))
            xout(:,:,k,l) = (sflag2(:,l2:u2-1)*xin(:,l2:u2-1,k,l)+&
                             sflag2(:,l2+1:u2)*xin(:,l2+1:u2,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_521
      ENDDO l_521

!
!5.2.2 With weights
!------------------
!

   ELSE
      l_522: DO l=1,nosize
      k_522: DO k=l3,u3
         WHERE (fdest3(:,:,k))
            xout(:,:,k,l) = &
               &(sflag2(:,l2:u2-1)*weights2(:,l2+1:u2)*xin(:,l2:u2-1,k,l)+&
               & sflag2(:,l2+1:u2)*weights2(:,l2:u2-1)*xin(:,l2+1:u2,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_522
      ENDDO l_522
   ENDIF
   
ENDIF

!
!5.3 Correction for variables with angle dimensions
!--------------------------------------------------
!

IF (fangle) THEN

!
!5.3.1 Using 2-D mask
!--------------------
!

   IF (iopt_arrint_3D.EQ.0) THEN
   
      l_531: DO l=1,nosize
      k_531: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = xout(:,:,k,l) + &
                          & pi*(1.0-&
                          & SIGN(1.0,pi-&
                               & ABS(xin(:,l2+1:u2,k,l)-xin(:,l2:u2-1,k,l))))
         END WHERE
      ENDDO k_531
      ENDDO l_531

!
!5.3.2 Using 3-D mask
!--------------------
!

   ELSE
   
      l_532: DO l=1,nosize
      k_532: DO k=l3,u3
         WHERE (fdest3(:,:,k))
            xout(:,:,k,l) = xout(:,:,k,l) + &
                          & pi*(1.0-&
                          & SIGN(1.0,pi-&
                               & ABS(xin(:,l2+1:u2,k,l)-xin(:,l2:u2-1,k,l))))
         END WHERE
      ENDDO k_532
      ENDDO l_532
      
   ENDIF

ENDIF

!
!6. Deallocate
!-------------
!

DEALLOCATE (sflag2,wsum2)
IF (iopt_arrint_3D.EQ.0) THEN
   DEALLOCATE (fdest2)
ELSE
   DEALLOCATE (fdest3)
ENDIF
IF (.NOT.hreg) DEALLOCATE (weights2)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Carr_at_V

!========================================================================

SUBROUTINE Carr_at_VW(xin,xout,intdest,lbounds,ubounds,nosize,ivarid,info,&
                    & outflag,hregular,vregular)
!************************************************************************
!
! *Carr_at_VW* Interpolate an array from C- to VW-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: hregular, vregular
INTEGER, INTENT(IN) :: intdest, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2)+1:ubounds(2),&
                           & lbounds(3)+1:ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at C-nodes
!*xout*      REAL    Output array at VW-nodes
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!*vregular*  LOGICAL Flag to select uniform or area averaging in the vertical
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: hvreg
INTEGER :: k


procname(pglev+1) = 'Carr_at_VW'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXregX
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(vregular)) THEN
   vreg = iopt_arrint_vreg.EQ.0
ELSE
   vreg = vregular
ENDIF
hvreg = hreg.AND.vreg

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ALLOCATE (fdest3(l1:u1,l2+1:u2,l3+1:u3),STAT=errstat)
CALL error_alloc('fdest3',3,(/u1-l1+1,u2-l2,u3-l3/),kndlog)

!---source
sflag2 = MERGE(1.0,0.0,nodeatc(l1:u1,l2:u2).GT.0)

!---destination
SELECT CASE (intdest)
   CASE (1); fdest3 = nodeatvw(l1:u1,l2+1:u2,l3+1:u3).GT.0
   CASE (2); fdest3 = nodeatvw(l1:u1,l2+1:u2,l3+1:u3).EQ.2
   CASE (3); fdest3 = nodeatvw(l1:u1,l2+1:u2,l3+1:u3).GT.1
   CASE (4); fdest3 = nodeatvw(l1:u1,l2+1:u2,l3+1:u3).EQ.1.OR.&
                    & nodeatvw(l1:u1,l2+1:u2,l3+1:u3).EQ.2
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
IF (hvreg) THEN
   ALLOCATE (wsum2(l1:u1,l2+1:u2),STAT=errstat)
   CALL error_alloc('wsum2',2,(/u1-l1+1,u2-l2/),kndrtype)
ELSE
   ALLOCATE (weights3(l1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('weights3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
   ALLOCATE (wsum3(l1:u1,l2+1:u2,l3+1:u3),STAT=errstat)
   CALL error_alloc('wsum3',3,(/u1-l1+1,u2-l2,u3-l3/),kndrtype)
ENDIF

!---define
IF (hvreg) THEN
   wsum2 = 2.0*(sflag2(:,l2+1:u2)+sflag2(:,l2:u2-1))
   k_411: DO k=l3+1,u3
      WHERE (wsum2.LE.0.0) fdest3(:,:,k) = .FALSE.
   ENDDO k_411
ELSE
   k_412: DO k=l3,u3
      weights3(:,:,k) = 1.0
      IF (.NOT.hreg) THEN
         weights3(:,:,k) = weights3(:,:,k)*delyatc(l1:u1,l2:u2)
      ENDIF
      IF (.NOT.vreg) THEN
         weights3(:,:,k) = weights3(:,:,k)*delzatc(l1:u1,l2:u2,k)
      ENDIF
      WHERE (nodeatv(l1:u1,l2+1:u2,k).GE.2.AND.nodeatc(l1:u1,l2:u2-1).GT.0)
         weights3(:,l2+1:u2,k) = weights3(:,l2:u2-1,k)
      END WHERE
      WHERE (nodeatv(l1:u1,l2+1:u2,k).GE.2.AND.nodeatc(l1:u1,l2+1:u2).GT.0)
         weights3(:,l2:u2-1,k) = weights3(:,l2+1:u2,k)
      END WHERE
   ENDDO k_412
   k_413: DO k=l3+1,u3
      wsum3(:,:,k) = sflag2(:,l2+1:u2)*&
                   &(weights3(:,l2:u2-1,k-1)+weights3(:,l2:u2-1,k)) + &
                   & sflag2(:,l2:u2-1)*&
                   &(weights3(:,l2+1:u2,k-1)+weights3(:,l2+1:u2,k))
   ENDDO k_413
   WHERE (wsum3.LE.0.0)
      fdest3 = .FALSE.
   END WHERE
ENDIF

!
!5. Interpolate
!--------------
!
!5.1 Without weights
!-------------------
!

IF (hvreg) THEN

   l_510: DO l=1,nosize
   k_510: DO k=l3+1,u3
      WHERE (fdest3(:,:,k))
         xout(:,:,k,l) =(sflag2(:,l2:u2-1)*&
                       &(xin(:,l2:u2-1,k-1,l)+xin(:,l2:u2-1,k,l))+&
                       & sflag2(:,l2+1:u2)*&
                       &(xin(:,l2+1:u2,k-1,l)+xin(:,l2+1:u2,k,l)))/wsum2
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_510
   ENDDO l_510

!
!5.2 With weights
!----------------
!

ELSE

   l_520: DO l=1,nosize
   k_520: DO k=l3+1,u3
      WHERE (fdest3(:,:,k))
         xout(:,:,k,l) = &
               &(sflag2(:,l2:u2-1)*weights3(:,l2+1:u2,k)*xin(:,l2:u2-1,k-1,l)+&
               & sflag2(:,l2+1:u2)*weights3(:,l2:u2-1,k)*xin(:,l2+1:u2,k-1,l)+&
               & sflag2(:,l2:u2-1)*weights3(:,l2+1:u2,k-1)*xin(:,l2:u2-1,k,l)+&
               & sflag2(:,l2+1:u2)*weights3(:,l2:u2-1,k-1)*xin(:,l2+1:u2,k,l))&
               & /wsum3(:,:,k)
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_520
   ENDDO l_520

ENDIF

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag2,fdest3)
IF (hvreg) THEN
   DEALLOCATE (wsum2)
ELSE
   DEALLOCATE (weights3,wsum3)
ENDIF

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Carr_at_VW

!========================================================================

SUBROUTINE Carr_at_W(xin,xout,lbounds,ubounds,nosize,ivarid,info,outflag,&
                   & vregular)
!************************************************************************
!
! *Carr_at_W* Interpolate an array from C- to W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.0
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: vregular
INTEGER, INTENT(IN) :: ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                           & lbounds(3)+1:ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at C-nodes
!*xout*      REAL    Output array at W-nodes
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*vregular*  LOGICAL Flag to select uniform or area averaging in the vertical
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Carr_at_W'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(vregular)) THEN
   vreg = iopt_arrint_vreg.EQ.0.OR.dZregZ
ELSE
   vreg = vregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Weight factors
!-----------------
!

IF (.NOT.vreg) THEN
!  ---allocate
   ALLOCATE (weights3(l1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('weights3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
   ALLOCATE (wsum3(l1:u1,l2:u2,l3+1:u3),STAT=errstat)
   CALL error_alloc('wsum3',3,(/u1-l1+1,u2-l2+1,u3-l3/),kndrtype)
!  ---define
   weights3 = delzatc(l1:u1,l2:u2,l3:u3)
   wsum3 = weights3(:,:,l3:u3-1) + weights3(:,:,l3+1:u3)
ENDIF

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (vreg) THEN

   l_410: DO l=1,nosize
   k_410: DO k=l3+1,u3
      WHERE (nodeatc(l1:u1,l2:u2).GT.0)
         xout(:,:,k,l) = 0.5*(xin(:,:,k-1,l)+xin(:,:,k,l))
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_410
   ENDDO l_410

!
!4.2 With weights
!----------------
!
      
ELSE

   l_420: DO l=1,nosize
   k_420: DO k=l3+1,u3
      WHERE (nodeatc(l1:u1,l2:u2).GT.0)
         xout(:,:,k,l) = (weights3(:,:,k)*xin(:,:,k-1,l)+&
                        & weights3(:,:,k-1)*xin(:,:,k,l))/wsum3(:,:,k)
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_420
   ENDDO l_420

ENDIF

!
!5. Deallocate
!-------------
!        

IF (.NOT.vreg) DEALLOCATE (weights3,wsum3)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Carr_at_W

!========================================================================

SUBROUTINE Uarr_at_C(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                   & info,outflag)
!************************************************************************
!
! *Uarr_at_C* Interpolate an array from U- to C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.6
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!

LOGICAL, INTENT(IN) :: info
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1)-1,lbounds(2):ubounds(2),&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*xout*      REAL    Output array at C-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => wet points only
!*lbounds*   INTEGER Lower bounds of input array
!*ubounds*   INTEGER Upper bounds of input array
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Uarr_at_C'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ELSE
   ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ENDIF
ALLOCATE (fdest2(l1:u1-1,l2:u2),STAT=errstat)
CALL error_alloc('fdest2',2,(/u1-l1,u2-l2+1/),kndlog)

!---source
IF (iopt_arrint_3D.EQ.0) THEN
   SELECT CASE (intsrce)
      CASE (0); sflag2 = 1.0
      CASE (1); sflag2 = MERGE(1.0,0.0,node2du(l1:u1,l2:u2).GT.0)
      CASE (2); sflag2 = MERGE(1.0,0.0,node2du(l1:u1,l2:u2).EQ.2)
      CASE (3); sflag2 = MERGE(1.0,0.0,node2du(l1:u1,l2:u2).GT.1)
      CASE (4); sflag2 = MERGE(1.0,0.0,node2du(l1:u1,l2:u2).EQ.1.OR.&
                                     & node2du(l1:u1,l2:u2).EQ.2)
   END SELECT
ELSE
   SELECT CASE (intsrce)
      CASE (0); sflag3 = 1.0
      CASE (1); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).GT.0)
      CASE (2); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).EQ.2)
      CASE (3); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).GT.1)
      CASE (4); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                  & nodeatu(l1:u1,l2:u2,l3:u3).EQ.2)
   END SELECT
ENDIF

!---destination
SELECT CASE (intdest)
   CASE (0); fdest2 = .TRUE.
   CASE (1); fdest2 = nodeatc(l1:u1-1,l2:u2).GT.0
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (wsum2(l1:u1-1,l2:u2),STAT=errstat)
   CALL error_alloc('wsum2',2,(/u1-l1,u2-l2+1/),kndrtype)
ELSE
   ALLOCATE (wsum3(l1:u1-1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('wsum3',3,(/u1-l1,u2-l2+1,u3-l3+1/),kndrtype)
ENDIF

!---define
IF (iopt_arrint_3D.EQ.0) THEN
   wsum2 = sflag2(l1:u1-1,:)+ sflag2(l1+1:u1,:)
ELSE
   wsum3 = sflag3(l1:u1-1,:,:) + sflag3(l1+1:u1,:,:)
ENDIF

!
!5. Interpolate
!--------------
!

IF (iopt_arrint_3D.EQ.0) THEN
   l_501: DO l=1,nosize
   k_501: DO k=l3,u3

      WHERE (fdest2.AND.wsum2.GT.0.0)
         xout(:,:,k,l) = (sflag2(l1:u1-1,:)*xin(l1:u1-1,:,k,l)+&
                        & sflag2(l1+1:u1,:)*xin(l1+1:u1,:,k,l))/wsum2
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_501
   ENDDO l_501

ELSE

   l_502: DO l=1,nosize
   k_502: DO k=l3,u3
      WHERE (fdest2.AND.wsum3(:,:,k).GT.0.0)
         xout(:,:,k,l) = (sflag3(l1:u1-1,:,k)*xin(l1:u1-1,:,k,l)+&
                        & sflag3(l1+1:u1,:,k)*xin(l1+1:u1,:,k,l))/wsum3(:,:,k)
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_502
   ENDDO l_502

ENDIF

!
!6. Deallocate
!-------------
!

DEALLOCATE (fdest2)        
IF (iopt_arrint_3D.EQ.0) THEN
   DEALLOCATE (sflag2,wsum2)
ELSE
   DEALLOCATE (sflag3,wsum3)
ENDIF

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Uarr_at_C

!========================================================================

SUBROUTINE Uarr_at_UV(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                    & info,outflag,hregular)
!************************************************************************
!
! *Uarr_at_UV* Interpolate an array from U- to UV-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.6
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2)+1:ubounds(2),&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*xout*      REAL    Output array at UV-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => interior points and open boundaries only
!                    = 2 => interior points and UV-node open boundaries only
!*lbounds*   INTEGER Lower bounds of input array
!*ubounds*   INTEGER Upper bounds of input array
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Uarr_at_UV'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dYregY
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
   ALLOCATE (fdest2(l1:u1,l2+1:u2),STAT=errstat)
   CALL error_alloc('fdest2',2,(/u1-l1+1,u2-l2/),kndlog)
ELSE
   ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('sflag3',2,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
   ALLOCATE (fdest3(l1:u1,l2+1:u2,l3:u3),STAT=errstat)
   CALL error_alloc('fdest3',3,(/u1-l1+1,u2-l2,u3-l3+1/),kndlog)
ENDIF

!---source
IF (iopt_arrint_3D.EQ.0) THEN
   SELECT CASE (intsrce)
      CASE (0); sflag2 = 1.0
      CASE (1); sflag2 = MERGE(1.0,0.0,node2du(l1:u1,l2:u2).GT.0)
      CASE (2); sflag2 = MERGE(1.0,0.0,node2du(l1:u1,l2:u2).EQ.2)
      CASE (3); sflag2 = MERGE(1.0,0.0,node2du(l1:u1,l2:u2).GT.1)
      CASE (4); sflag2 = MERGE(1.0,0.0,node2du(l1:u1,l2:u2).EQ.1.OR.&
                                     & node2du(l1:u1,l2:u2).EQ.2)
   END SELECT
ELSE
   SELECT CASE (intsrce)
      CASE (0); sflag3 = 1.0
      CASE (1); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).GT.0)
      CASE (2); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).EQ.2)
      CASE (3); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).GT.1)
      CASE (4); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                     & nodeatu(l1:u1,l2:u2,l3:u3).EQ.2)
   END SELECT
ENDIF

!---destination
IF (iopt_arrint_3D.EQ.0) THEN
   SELECT CASE (intdest)
      CASE (0); fdest2 = .TRUE.
      CASE (1); fdest2 = node2duv(l1:u1,l2+1:u2).GT.0
      CASE (2); fdest2 = node2duv(l1:u1,l2+1:u2).EQ.1.OR.&
                       & node2duv(l1:u1,l2+1:u2).EQ.3
   END SELECT
ELSE
   SELECT CASE (intdest)
      CASE (0); fdest3 = .TRUE.
      CASE (1); fdest3 = nodeatuv(l1:u1,l2+1:u2,l3:u3).GT.0
      CASE (2); fdest3 = nodeatuv(l1:u1,l2+1:u2,l3:u3).EQ.1.OR.&
                       & nodeatuv(l1:u1,l2+1:u2,l3:u3).EQ.3
   END SELECT
ENDIF

!
!4. Weight factors
!-----------------
!
!---allocate
IF (.NOT.hreg) THEN
   ALLOCATE (weights2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('weights2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ENDIF
IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (wsum2(l1:u1,l2+1:u2),STAT=errstat)
   CALL error_alloc('wsum2',2,(/u1-l1+1,u2-l2/),kndrtype)
ELSE
   ALLOCATE (wsum3(l1:u1,l2+1:u2,l3:u3),STAT=errstat)
   CALL error_alloc('wsum3',3,(/u1-l1+1,u2-l2,u3-l3+1/),kndrtype)
ENDIF

!---define
IF (iopt_arrint_3D.EQ.0) THEN
   IF (hreg) THEN
      wsum2 = sflag2(:,l2+1:u2) + sflag2(:,l2:u2-1)
   ELSE
      weights2 = delyatu(l1:u1,l2:u2)
      wsum2 = sflag2(:,l2+1:u2)*weights2(:,l2:u2-1) + &
            & sflag2(:,l2:u2-1)*weights2(:,l2+1:u2)
   ENDIF
   fdest2 = fdest2.AND.wsum2.GT.0.0
ELSE
   IF (hreg) THEN
      wsum3 = sflag3(:,l2+1:u2,:) + sflag3(:,l2:u2-1,:)
   ELSE
      weights2 = delyatu(l1:u1,l2:u2)
      k_410: DO k=l3,u3
         wsum3(:,:,k) = sflag3(:,l2+1:u2,k)*weights2(:,l2:u2-1) + &
                       & sflag3(:,l2:u2-1,k)*weights2(:,l2+1:u2)
      ENDDO k_410
   ENDIF
   WHERE (wsum3.LE.0.0) fdest3 = .FALSE.
ENDIF

!
!5. Interpolate
!--------------
!
!5.1 Using 2-D mask
!------------------
!

IF (iopt_arrint_3D.EQ.0) THEN

!
!5.1.1 Without weights
!---------------------
!

   IF (hreg) THEN
      l_511: DO l=1,nosize
      k_511: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = (sflag2(:,l2:u2-1)*xin(:,l2:u2-1,k,l)+&
                           & sflag2(:,l2+1:u2)*xin(:,l2+1:u2,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_511
      ENDDO l_511

!
!5.1.2 With weights
!------------------
!

   ELSE
      l_512: DO l=1,nosize
      k_512: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = &
               &(sflag2(:,l2:u2-1)*weights2(:,l2+1:u2)*xin(:,l2:u2-1,k,l)+&
               & sflag2(:,l2+1:u2)*weights2(:,l2:u2-1)*xin(:,l2+1:u2,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_512
      ENDDO l_512
   ENDIF

!
!5.2 Using 3-D mask
!------------------
!

ELSE

!
!5.2.1 Without weights
!---------------------
!

   IF (hreg) THEN
      l_521: DO l=1,nosize
         WHERE (fdest3)
            xout(:,:,:,l) = (sflag3(:,l2:u2-1,:)*xin(:,l2:u2-1,:,l)+&
                           & sflag3(:,l2+1:u2,:)*xin(:,l2+1:u2,:,l))/wsum3
         ELSEWHERE
            xout(:,:,:,l) = flagval
         END WHERE
      ENDDO l_521

!
!5.2.2 With weights
!------------------
!

   ELSE
      l_522: DO l=1,nosize
      k_522: DO k=l3,u3
         WHERE (fdest3(:,:,k))
            xout(:,:,k,l) = &
               &(sflag3(:,l2:u2-1,k)*weights2(:,l2+1:u2)*xin(:,l2:u2-1,k,l)+&
               & sflag3(:,l2+1:u2,k)*weights2(:,l2:u2-1)*xin(:,l2+1:u2,k,l))&
               & /wsum3(:,:,k)
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_522
      ENDDO l_522
   ENDIF

ENDIF

!
!6. Deallocate
!-------------
!

IF (iopt_arrint_3D.EQ.0) THEN
   DEALLOCATE (sflag2,fdest2,wsum2)
ELSE
   DEALLOCATE (sflag3,fdest3,wsum3)
ENDIF
IF (.NOT.hreg) DEALLOCATE (weights2)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Uarr_at_UV

!========================================================================

SUBROUTINE Uarr_at_UW(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                    & info,outflag,vregular)
!************************************************************************
!
! *Uarr_at_UW* Interpolate an array from U- to UW-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: vregular
INTEGER, INTENT(IN) :: intdest, intsrce,ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                           & lbounds(3)+1:ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*xout*      REAL    Output array at UW-nodes
!*intnode*   INTEGER Selects vqlid points at source and destination node
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*vregular*  LOGICAL Flag to select uniform or area averaging in the vertical
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Uarr_at_UW'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(vregular)) THEN
   vreg = iopt_arrint_vreg.EQ.0.OR.dZregZ
ELSE
   vreg = vregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ALLOCATE (fdest3(l1:u1,l2:u2,l3+1:u3),STAT=errstat)
CALL error_alloc('fdest3',3,(/u1-l1+1,u2-l2+1,u3-l3/),kndlog)

!---source
SELECT CASE (intsrce)
   CASE (1); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).GT.0)
   CASE (2); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).EQ.2)
   CASE (3); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).GT.1)
   CASE (4); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                  & nodeatu(l1:u1,l2:u2,l3:u3).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (1); fdest3 = nodeatuw(l1:u1,l2:u2,l3+1:u3).GT.0
   CASE (2); fdest3 = nodeatuw(l1:u1,l2:u2,l3+1:u3).EQ.2
   CASE (3); fdest3 = nodeatuw(l1:u1,l2:u2,l3+1:u3).GT.1
   CASE (4); fdest3 = nodeatuw(l1:u1,l2:u2,l3+1:u3).EQ.1.OR.&
                    & nodeatuw(l1:u1,l2:u2,l3+1:u3).EQ.2
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
IF (.NOT.vreg) THEN
   ALLOCATE (weights3(l1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('weights3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ENDIF
ALLOCATE (wsum3(l1:u1,l2:u2,l3+1:u3),STAT=errstat)
CALL error_alloc('wsum3',3,(/u1-l1+1,u2-l2+1,u3-l3/),kndrtype)

!---define
IF (vreg) THEN
   wsum3 = sflag3(:,:,l3:u3-1) + sflag3(:,:,l3+1:u3)
ELSE
   weights3 = delzatu(l1:u1,l2:u2,l3:u3)
   wsum3 = sflag3(:,:,l3+1:u3)*weights3(:,:,l3:u3-1) + &
         & sflag3(:,:,l3:u3-1)*weights3(:,:,l3+1:u3)
ENDIF

WHERE (wsum3.LE.0.0) fdest3 = .FALSE.

!
!5. Interpolate
!--------------
!
!5.1 Without weights
!-------------------
!

IF (vreg) THEN
   
   l_510: DO l=1,nosize
      WHERE (fdest3)
         xout(:,:,:,l) = sflag3(:,:,l3:u3-1)*xin(:,:,l3:u3-1,l) + &
                       & sflag3(:,:,l3+1:u3)*xin(:,:,l3+1:u3,l)
      ELSEWHERE
         xout(:,:,:,l) = flagval
      END WHERE
   ENDDO l_510

!
!5.2 With weights
!----------------
!

ELSE

   l_520: DO l=1,nosize
      WHERE (fdest3)
         xout(:,:,:,l) = &
               & (sflag3(:,:,l3:u3-1)*weights3(:,:,l3+1:u3)*xin(:,:,l3:u3-1,l)+&
               &  sflag3(:,:,l3+1:u3)*weights3(:,:,l3:u3-1)*xin(:,:,l3+1:u3,l))&
               & /wsum3
      ELSEWHERE
         xout(:,:,:,l) = flagval
      END WHERE
   ENDDO l_520

ENDIF

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag3,fdest3,wsum3)
IF (.NOT.vreg) DEALLOCATE (weights3)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Uarr_at_UW

!========================================================================

SUBROUTINE Uarr_at_V(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                   & info,outflag,hregular)
!************************************************************************
!
! *Uarr_at_V* Interpolate an array from U- to V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.6
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!


LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: intsrce, intdest, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1)-1,lbounds(2)+1:ubounds(2),&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*xout*      REAL    Output array at V-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array
!*ubounds*   INTEGER Upper bounds of input array
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Uarr_at_V'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXYreg
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
   ALLOCATE (fdest2(l1:u1-1,l2+1:u2),STAT=errstat)
   CALL error_alloc('fdest2',2,(/u1-l1,u2-l2/),kndlog)
ELSE
   ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('sflag3',2,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
   ALLOCATE (fdest3(l1:u1-1,l2+1:u2,l3:u3),STAT=errstat)
   CALL error_alloc('fdest3',3,(/u1-l1,u2-l2,u3-l3+1/),kndlog)
ENDIF

!---source
IF (iopt_arrint_3D.EQ.0) THEN
   SELECT CASE (intsrce)
      CASE (0); sflag2 = 1.0
      CASE (1); sflag2 = MERGE(1.0,0.0,node2du(l1:u1,l2:u2).GT.0)
      CASE (2); sflag2 = MERGE(1.0,0.0,node2du(l1:u1,l2:u2).EQ.2)
      CASE (3); sflag2 = MERGE(1.0,0.0,node2du(l1:u1,l2:u2).GT.1)
      CASE (4); sflag2 = MERGE(1.0,0.0,node2du(l1:u1,l2:u2).EQ.1.OR.&
                                     & node2du(l1:u1,l2:u2).EQ.2)
   END SELECT
ELSE
   SELECT CASE (intsrce)
      CASE (0); sflag3 = 1.0
      CASE (1); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).GT.0)
      CASE (2); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).EQ.2)
      CASE (3); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).GT.1)
      CASE (4); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                  & nodeatu(l1:u1,l2:u2,l3:u3).EQ.2)
   END SELECT
ENDIF

!---destination
IF (iopt_arrint_3D.EQ.0) THEN
   SELECT CASE (intdest)
      CASE (0); fdest2 = .TRUE.
      CASE (1); fdest2 = node2dv(l1:u1-1,l2+1:u2).GT.0
      CASE (2); fdest2 = node2dv(l1:u1-1,l2+1:u2).EQ.2
      CASE (3); fdest2 = node2dv(l1:u1-1,l2+1:u2).GT.1
      CASE (4); fdest2 = node2dv(l1:u1-1,l2+1:u2).EQ.1.OR.&
                       & node2dv(l1:u1-1,l2+1:u2).EQ.2
   END SELECT
ELSE
   SELECT CASE (intdest)
      CASE (0); fdest3 = .TRUE.
      CASE (1); fdest3 = nodeatv(l1:u1-1,l2+1:u2,l3:u3).GT.0
      CASE (2); fdest3 = nodeatv(l1:u1-1,l2+1:u2,l3:u3).EQ.2
      CASE (3); fdest3 = nodeatv(l1:u1-1,l2+1:u2,l3:u3).GT.1
      CASE (4); fdest3 = nodeatv(l1:u1-1,l2+1:u2,l3:u3).EQ.1.OR.&
                       & nodeatv(l1:u1-1,l2+1:u2,l3:u3).EQ.2
   END SELECT
ENDIF

!
!4. Weight factors
!-----------------
!
!---allocate
IF (.NOT.hreg) THEN
   ALLOCATE (weights2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('weights2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ENDIF
IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (wsum2(l1:u1-1,l2+1:u2),STAT=errstat)
   CALL error_alloc('wsum2',2,(/u1-l1,u2-l2/),kndrtype)
ELSE
   ALLOCATE (wsum3(l1:u1-1,l2+1:u2,l3:u3),STAT=errstat)
   CALL error_alloc('wsum3',3,(/u1-l1,u2-l2,u3-l3+1/),kndrtype)
ENDIF

!---define
IF (iopt_arrint_3D.EQ.0) THEN
   IF (hreg) THEN
      wsum2 = sflag2(l1+1:u1,l2+1:u2) + sflag2(l1:u1-1,l2+1:u2) + &
            & sflag2(l1+1:u1,l2:u2-1) + sflag2(l1:u1-1,l2:u2-1)
   ELSE
      IF (dXregX.AND.dXregY) THEN 
         weights2 = delyatv(l1:u1,l2:u2)
      ELSEIF (dYregX.AND.dYregY) THEN
         weights2 = delxatu(l1:u1,l2:u2)
      ELSE
         weights2 = delxatu(l1:u1,l2:u2)*delyatv(l1:u1,l2:u2)
      ENDIF
      wsum2 = sflag2(l1+1:u1,l2+1:u2)*weights2(l1:u1-1,l2:u2-1) + &
            & sflag2(l1:u1-1,l2+1:u2)*weights2(l1+1:u1,l2:u2-1) + &
            & sflag2(l1+1:u1,l2:u2-1)*weights2(l1:u1-1,l2+1:u2) + &
            & sflag2(l1:u1-1,l2:u2-1)*weights2(l1+1:u1,l2+1:u2)
   ENDIF
   WHERE (wsum2.LE.0.0) fdest2 = .FALSE.
ELSE
   IF (hreg) THEN
      wsum3 = sflag3(l1+1:u1,l2+1:u2,:) + sflag3(l1:u1-1,l2+1:u2,:) + &
            & sflag3(l1+1:u1,l2:u2-1,:) + sflag3(l1:u1-1,l2:u2-1,:)
   ELSE
      IF (dXregX.AND.dXregY) THEN 
         weights2 = delyatv(l1:u1,l2:u2)
      ELSEIF (dYregX.AND.dYregY) THEN
         weights2 = delxatu(l1:u1,l2:u2)
      ELSE
         weights2 = delxatu(l1:u1,l2:u2)*delyatv(l1:u1,l2:u2)
      ENDIF
      k_410: DO k=l3,u3
         wsum3(:,:,k) = sflag3(l1+1:u1,l2+1:u2,k)*weights2(l1:u1-1,l2:u2-1) + &
                      & sflag3(l1:u1-1,l2+1:u2,k)*weights2(l1+1:u1,l2:u2-1) + &
                      & sflag3(l1+1:u1,l2:u2-1,k)*weights2(l1:u1-1,l2+1:u2) + &
                      & sflag3(l1:u1-1,l2:u2-1,k)*weights2(l1+1:u1,l2+1:u2)
      ENDDO k_410
   ENDIF
   WHERE (wsum3.LE.0.0) fdest3 = .FALSE.
ENDIF

!
!5. Interpolate
!--------------
!
!5.1 Using 2-D mask
!------------------
!

IF (iopt_arrint_3D.EQ.0) THEN

!
!5.1.1 Without weights
!---------------------
!

   IF (hreg) THEN
      l_511: DO l=1,nosize
      k_511: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = &
                       &(sflag2(l1:u1-1,l2:u2-1)*xin(l1:u1-1,l2:u2-1,k,l)+&
                       & sflag2(l1+1:u1,l2:u2-1)*xin(l1+1:u1,l2:u2-1,k,l)+&
                       & sflag2(l1:u1-1,l2+1:u2)*xin(l1:u1-1,l2+1:u2,k,l)+&
                       & sflag2(l1+1:u1,l2+1:u2)*xin(l1+1:u1,l2+1:u2,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_511
      ENDDO l_511

!
!5.1.2 With weights
!------------------
!

   ELSE
      l_512: DO l=1,nosize
      k_512: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = &
                   &(sflag2(l1:u1-1,l2:u2-1)*weights2(l1+1:u1,l2+1:u2)*&
                   & xin(l1:u1-1,l2:u2-1,k,l)+&
                   & sflag2(l1+1:u1,l2:u2-1)*weights2(l1:u1-1,l2+1:u2)*&
                   & xin(l1+1:u1,l2:u2-1,k,l)+&
                   & sflag2(l1:u1-1,l2+1:u2)*weights2(l1+1:u1,l2:u2-1)*&
                   & xin(l1:u1-1,l2+1:u2,k,l)+&
                   & sflag2(l1+1:u1,l2+1:u2)*weights2(l1:u1-1,l2:u2-1)*&
                   & xin(l1+1:u1,l2+1:u2,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_512
      ENDDO l_512
   ENDIF

!
!5.2 Using 3-D mask
!------------------
!

ELSE

!
!5.2.1 Without weights
!---------------------
!

   IF (hreg) THEN
      l_521: DO l=1,nosize
         WHERE (fdest3)
            xout(:,:,:,l) = &
                     &(sflag3(l1:u1-1,l2:u2-1,:)*xin(l1:u1-1,l2:u2-1,:,l)+&
                     & sflag3(l1+1:u1,l2:u2-1,:)*xin(l1+1:u1,l2:u2-1,:,l)+&
                     & sflag3(l1:u1-1,l2+1:u2,:)*xin(l1:u1-1,l2+1:u2,:,l)+&
                     & sflag3(l1+1:u1,l2+1:u2,:)*xin(l1+1:u1,l2+1:u2,:,l))/wsum3
         ELSEWHERE
            xout(:,:,:,l) = flagval
         END WHERE
      ENDDO l_521

!
!5.2.2 With weights
!------------------
!

   ELSE
      l_522: DO l=1,nosize
      k_522: DO k=l3,u3
         WHERE (fdest3(:,:,k))
            xout(:,:,k,l) = &
                &(sflag3(l1:u1-1,l2:u2-1,k)*weights2(l1+1:u1,l2+1:u2)*&
                & xin(l1:u1-1,l2:u2-1,k,l)+&
                & sflag3(l1+1:u1,l2:u2-1,k)*weights2(l1:u1-1,l2+1:u2)*&
                & xin(l1+1:u1,l2:u2-1,k,l)+&
                & sflag3(l1:u1-1,l2+1:u2,k)*weights2(l1+1:u1,l2:u2-1)*&
                & xin(l1:u1-1,l2+1:u2,k,l)+&
                & sflag3(l1+1:u1,l2+1:u2,k)*weights2(l1:u1-1,l2:u2-1)*&
                & xin(l1+1:u1,l2+1:u2,k,l))/wsum3(:,:,k)
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_522
      ENDDO l_522
   ENDIF

ENDIF

!
!6. Deallocate
!-------------
!        

IF (iopt_arrint_3D.EQ.0) THEN
   DEALLOCATE (sflag2,fdest2,wsum2)
ELSE
   DEALLOCATE (sflag3,fdest3,wsum3)
ENDIF
IF (.NOT.hreg) DEALLOCATE (weights2)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Uarr_at_V

!========================================================================

SUBROUTINE Uarr_at_W(xin,xout,intsrce,lbounds,ubounds,nosize,ivarid,info,&
                   & outflag,vregular)
!************************************************************************
!
! *Uarr_at_W* Interpolate an array from U- to W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: vregular
INTEGER, INTENT(IN) :: intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1)-1,lbounds(2):ubounds(2),&
                           & lbounds(3)+1:ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*xout*      REAL    Output array at W-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*vregular*  LOGICAL Flag to select uniform or area averaging in the vertical
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Uarr_at_W'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(vregular)) THEN
   vreg = iopt_arrint_vreg.EQ.0.OR.dZregZ
ELSE
   vreg = vregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ALLOCATE (fdest2(l1:u1-1,l2:u2),STAT=errstat)
CALL error_alloc('fdest2',2,(/u1-l1,u2-l2+1/),kndlog)

!---source
SELECT CASE (intsrce)
   CASE (1); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).GT.0)
   CASE (2); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).EQ.2)
   CASE (3); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).GT.1)
   CASE (4); sflag3 = MERGE(1.0,0.0,nodeatu(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                  & nodeatu(l1:u1,l2:u2,l3:u3).EQ.2)
END SELECT

!---destination
fdest2 = nodeatc(l1:u1-1,l2:u2).GT.0

!
!4. Weight factors
!-----------------
!
!---allocate
IF (.NOT.vreg) THEN
   ALLOCATE (weights3(l1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('weights3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ENDIF
ALLOCATE (wsum3(l1:u1-1,l2:u2,l3+1:u3),STAT=errstat)
CALL error_alloc('wsum3',3,(/u1-l1,u2-l2+1,u3-l3/),kndrtype)

!---define
IF (vreg) THEN
   wsum3 = sflag3(l1:u1-1,:,l3:u3-1)+sflag3(l1+1:u1,:,l3:u3-1)+&
         & sflag3(l1:u1-1,:,l3+1:u3)+sflag3(l1+1:u1,:,l3+1:u3)
   WHERE (wsum3.LE.0.0) fdest3 = .FALSE.
ELSE
   k_410: DO k=l3,u3
      weights3(:,:,k) = delzatu(l1:u1,l2:u2,k)
      WHERE (nodeatu(l1+1:u1,l2:u2,k).GE.2.AND.nodeatc(l1:u1-1,l2:u2).GT.0)
         weights3(l1+1:u1,:,k) = delzatc(l1:u1-1,l2:u2,k)
      END WHERE
      WHERE (nodeatu(l1:u1-1,l2:u2,k).GE.2.AND.nodeatc(l1:u1-1,l2:u2).GT.0)
         weights3(l1:u1-1,:,k) = delzatc(l1:u1-1,l2:u2,k)
      END WHERE
   ENDDO k_410
   wsum3(:,:,l3+1:u3) = &
                     & sflag3(l1:u1-1,:,l3:u3-1)*weights3(l1+1:u1,:,l3+1:u3) + &
                     & sflag3(l1+1:u1,:,l3:u3-1)*weights3(l1:u1-1,:,l3+1:u3) + &
                     & sflag3(l1:u1-1,:,l3+1:u3)*weights3(l1+1:u1,:,l3:u3-1) + &
                     & sflag3(l1+1:u1,:,l3+1:u3)*weights3(l1:u1-1,:,l3:u3-1)
ENDIF

!
!5. Interpolate
!--------------
!
!5.1 Without weights
!-------------------
!

IF (vreg) THEN

   l_510: DO l=1,nosize
   k_510: DO k=l3+1,u3
      WHERE (fdest2.AND.wsum3(:,:,k).GT.0.0)
         xout(:,:,k,l) = (sflag3(l1:u1-1,:,k-1)*xin(l1:u1-1,:,k-1,l)+&
                        & sflag3(l1:u1-1,:,k)*xin(l1:u1-1,:,k,l)+&
                        & sflag3(l1+1:u1,:,k-1)*xin(l1+1:u1,:,k-1,l)+&
                        & sflag3(l1+1:u1,:,k)*xin(l1+1:u1,:,k,l))&
                        & /wsum3(:,:,k)
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_510
   ENDDO l_510

!
!5.2 With weights
!----------------
!

ELSE

   l_520: DO l=1,nosize
   k_520: DO k=l3+1,u3
      WHERE (fdest2.AND.wsum3(:,:,k).GT.0.0)
         xout(:,:,k,l) = &
            &(sflag3(l1:u1-1,:,k-1)*weights3(l1+1:u1,:,k)*xin(l1:u1-1,:,k-1,l)+&
            & sflag3(l1+1:u1,:,k-1)*weights3(l1:u1-1,:,k)*xin(l1+1:u1,:,k-1,l)+&
            & sflag3(l1:u1-1,:,k)*weights3(l1+1:u1,:,k-1)*xin(l1:u1-1,:,k,l)+&
            & sflag3(l1+1:u1,:,k)*weights3(l1:u1-1,:,k-1)*xin(l1+1:u1,:,k,l))&
            & /wsum3(:,:,k)
       ELSEWHERE
          xout(:,:,k,l) = flagval
       END WHERE
    ENDDO k_520
    ENDDO l_520

ENDIF

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag3,fdest2,wsum3)
IF (.NOT.vreg) DEALLOCATE (weights3)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Uarr_at_W

!========================================================================

SUBROUTINE UVarr_at_C(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                    & info,outflag)
!************************************************************************
!
! *UVarr_at_C* Interpolate an array from UV- to C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*Arguments
!
LOGICAL, INTENT(IN) :: info
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1)-1,lbounds(2):ubounds(2)-1,&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at UV-nodes
!*xout*      REAL    Output array at C-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => wet points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => wet points only
!*lbounds*   INTEGER Lower bounds of input array 
!*ubounds*   INTEGER Upper bounds of input array
!*nosize*    INTEGER Fourth dimension of input array
!*intnode*   INTEGER Dry/wet weight factor for interpolation
!                    = 0 => all points
!                    = 1 => wet corner points only
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'UVarr_at_C'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ALLOCATE (fdest2(l1:u1-1,l2:u2-1),STAT=errstat)
CALL error_alloc('fdest2',2,(/u1-l1,u2-l2/),kndlog)

!---source
SELECT CASE (intsrce)
   CASE (0); sflag3 = 1.0
   CASE (1); sflag3 = MERGE(1.0,0.0,nodeatuv(l1:u1,l2:u2,l3:u3).GT.0)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest2 = .TRUE.
   CASE (1); fdest2 = nodeatc(l1:u1-1,l2:u2-1).GT.0
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
ALLOCATE (wsum3(l1:u1-1,l2:u2-1,l3:u3),STAT=errstat)
CALL error_alloc('wsum3',3,(/u1-l1,u2-l2,u3-l3+1/),kndrtype)

!---define
wsum3 = sflag3(l1:u1-1,l2:u2-1,:) + sflag3(l1+1:u1,l2:u2-1,:) + &
      & sflag3(l1:u1-1,l2+1:u2,:) + sflag3(l1+1:u1,l2+1:u2,:) 

!
!5. Interpolate
!--------------
!

l_510: DO l=1,nosize
k_510: DO k=l3,u3
   WHERE (fdest2.AND.wsum3(:,:,k).GT.0.0)
      xout(:,:,k,l) = &
             &(sflag3(l1:u1-1,l2:u2-1,k)*xin(l1:u1-1,l2:u2-1,k,l)+&
             & sflag3(l1+1:u1,l2:u2-1,k)*xin(l1+1:u1,l2:u2-1,k,l)+&
             & sflag3(l1:u1-1,l2+1:u2,k)*xin(l1:u1-1,l2+1:u2,k,l)+&
             & sflag3(l1+1:u1,l2+1:u2,k)*xin(l1+1:u1,l2+1:u2,k,l))/wsum3(:,:,k)
   ELSEWHERE
      xout(:,:,k,l) = flagval
   END WHERE
ENDDO k_510
ENDDO l_510

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag3,fdest2,wsum3)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE UVarr_at_C

!========================================================================

SUBROUTINE UVarr_at_U(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                    & info,outflag)
!************************************************************************
!
! *UVarr_at_U* Interpolate an array from UV- to U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*Arguments
!
LOGICAL, INTENT(IN) :: info
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2)-1,&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at UV-nodes
!*xout*      REAL    Output array at U-nodes
!*intsrce*   INTEGER Selects valid points at source and destination node
!                    = 0 => all points
!                    = 1 => interior points and open boundaries only
!                    = 2 => interior points and UV-node open boundaries only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array
!*ubounds*   INTEGER Upper bounds of input array
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'UVarr_at_U'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ALLOCATE (fdest3(l1:u1,l2:u2-1,l3:u3),STAT=errstat)
CALL error_alloc('fdest3',3,(/u1-l1+1,u2-l2,u3-l3+1/),kndlog)

!---source
SELECT CASE (intsrce)
   CASE (0); sflag3 = 1.0
   CASE (1); sflag3 = MERGE(1.0,0.0,nodeatuv(l1:u1,l2:u2,l3:u3).GT.0)
   CASE (2); sflag3 = MERGE(1.0,0.0,nodeatuv(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                  & nodeatuv(l1:u1,l2:u2,l3:u3).EQ.3)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest3 = .TRUE.
   CASE (1); fdest3 = nodeatu(l1:u1,l2:u2-1,l3:u3).GT.0
   CASE (2); fdest3 = nodeatu(l1:u1,l2:u2-1,l3:u3).EQ.2
   CASE (3); fdest3 = nodeatu(l1:u1,l2:u2-1,l3:u3).GT.1
   CASE (4); fdest3 = nodeatu(l1:u1,l2:u2-1,l3:u3).EQ.1.OR.&
                    & nodeatu(l1:u1,l2:u2-1,l3:u3).EQ.2
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
ALLOCATE (wsum3(l1:u1,l2:u2-1,l3:u3),STAT=errstat)
CALL error_alloc('wsum3',3,(/u1-l1+1,u2-l2,u3-l3+1/),kndrtype)

!---define
wsum3 = sflag3(:,l2:u2-1,:) + sflag3(:,l2+1:u2,:)

WHERE (wsum3.LE.0.0) fdest3 = .FALSE.

!
!5. Interpolate
!--------------
!

l_510: DO l=1,nosize
   WHERE (fdest3)
      xout(:,:,:,l) = (sflag3(:,l2:u2-1,:)*xin(:,l2:u2-1,:,l)+&
                     & sflag3(:,l2+1:u2,:)*xin(:,l2+1:u2,:,l))/wsum3
   ELSEWHERE
      xout(:,:,:,l) = flagval
   END WHERE
ENDDO l_510

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag3,fdest3,wsum3)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE UVarr_at_U

!========================================================================

SUBROUTINE UVarr_at_V(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                    & info,outflag)
!************************************************************************
!
! *UVarr_at_V* Interpolate an array from UV- to V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*Arguments
!
LOGICAL, INTENT(IN) :: info
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1)-1,lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at UV-nodes
!*xout*      REAL    Output array at V-nodes
!*intsrce*   INTEGER Selects valid points at source and destination node
!                    = 0 => all points
!                    = 1 => interior points and open boundaries only
!                    = 2 => interior points and UV-node open boundaries only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array
!*ubounds*   INTEGER Upper bounds of input array
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'UVarr_at_V'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ALLOCATE (fdest3(l1:u1-1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('fdest3',3,(/u1-l1,u2-l2+1,u3-l3+1/),kndlog)

!---source
SELECT CASE (intsrce)
   CASE (0); sflag3 = 1.0
   CASE (1); sflag3 = MERGE(1.0,0.0,nodeatuv(l1:u1,l2:u2,l3:u3).GT.0)
   CASE (2); sflag3 = MERGE(1.0,0.0,nodeatuv(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                  & nodeatuv(l1:u1,l2:u2,l3:u3).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest3 = .TRUE.
   CASE (1); fdest3 = nodeatv(l1:u1-1,l2:u2,l3:u3).GT.0
   CASE (2); fdest3 = nodeatv(l1:u1-1,l2:u2,l3:u3).EQ.2
   CASE (3); fdest3 = nodeatv(l1:u1-1,l2:u2,l3:u3).GT.1
   CASE (4); fdest3 = nodeatv(l1:u1-1,l2:u2,l3:u3).EQ.1.OR.&
                    & nodeatv(l1:u1-1,l2:u2,l3:u3).EQ.2
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
ALLOCATE (wsum3(l1:u1-1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('wsum3',3,(/u1-l1,u2-l2+1,u3-l3+1/),kndrtype)

!---define
wsum3 = sflag3(l1:u1-1,:,:) + sflag3(l1+1:u1,:,:)

WHERE (wsum3.LE.0.0) fdest3 = .FALSE.

!
!5. Interpolate
!--------------
!

l_510: DO l=1,nosize
   WHERE (fdest3)
      xout(:,:,:,l) = (sflag3(l1:u1-1,:,:)*xin(l1:u1-1,:,:,l)+&
                     & sflag3(l1+1:u1,:,:)*xin(l1+1:u1,:,:,l))/wsum3
   ELSEWHERE
      xout(:,:,:,l) = flagval
   END WHERE
ENDDO l_510

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag3,fdest3,wsum3)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE UVarr_at_V

!========================================================================

SUBROUTINE UWarr_at_U(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                    & info,outflag)
!************************************************************************
!
! *UWarr_at_U* Interpolate an array from UW- to U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                           & lbounds(3):ubounds(3)-1,nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at UW-nodes
!*xout*      REAL    Output array at U-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'UWarr_at_U'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ALLOCATE (fdest3(l1:u1,l2:u2,l3:u3-1),STAT=errstat)
CALL error_alloc('fdest3',3,(/u1-l1+1,u2-l2+1,u3-l3/),kndlog)

!---source
SELECT CASE (intsrce)
   CASE (1); sflag3 = MERGE(1.0,0.0,nodeatuw(l1:u1,l2:u2,l3:u3).GT.0)
   CASE (2); sflag3 = MERGE(1.0,0.0,nodeatuw(l1:u1,l2:u2,l3:u3).EQ.2)
   CASE (3); sflag3 = MERGE(1.0,0.0,nodeatuw(l1:u1,l2:u2,l3:u3).GT.1)
   CASE (4); sflag3 = MERGE(1.0,0.0,nodeatuw(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                  & nodeatuw(l1:u1,l2:u2,l3:u3).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE(1); fdest3 = nodeatu(l1:u1,l2:u2,l3:u3-1).GT.0
   CASE(2); fdest3 = nodeatu(l1:u1,l2:u2,l3:u3-1).EQ.2
   CASE(3); fdest3 = nodeatu(l1:u1,l2:u2,l3:u3-1).GT.1
   CASE(4); fdest3 = nodeatu(l1:u1,l2:u2,l3:u3-1).EQ.1.OR.&
                   & nodeatu(l1:u1,l2:u2,l3:u3-1).EQ.2
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
ALLOCATE (wsum3(l1:u1,l2:u2,l3:u3-1),STAT=errstat)
CALL error_alloc('wsum3',3,(/u1-l1+1,u2-l2+1,u3-l3/),kndrtype)

!---define
wsum3 = sflag3(:,:,l3:u3-1) + sflag3(:,:,l3+1:u3)

WHERE (wsum3.LE.0.0) fdest3 = .FALSE.

!
!5. Interpolate
!--------------
!

l_510: DO l=1,nosize
   WHERE (fdest3)
      xout(:,:,:,l) = (sflag3(:,:,l3:u3-1)*xin(:,:,l3:u3-1,l)+&
                     & sflag3(:,:,l3+1:u3)*xin(:,:,l3+1:u3,l))/wsum3
   ELSEWHERE
      xout(:,:,:,l) = flagval
   END WHERE
ENDDO l_510

!
!6. Deallocate
!-------------
!

DEALLOCATE (sflag3,fdest3,wsum3)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE UWarr_at_U

!========================================================================

SUBROUTINE UWarr_at_W(xin,xout,intsrce,lbounds,ubounds,nosize,ivarid,info,&
                    & outflag)
!************************************************************************
!
! *UWarr_at_W* Interpolate an array from UW- to W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
INTEGER, INTENT(IN) :: intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1)-1,lbounds(2):ubounds(2),&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at UW-nodes
!*xout*      REAL    Output array at W-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'UWarr_at_W'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ALLOCATE (fdest2(l1:u1-1,l2:u2),STAT=errstat)
CALL error_alloc('fdest2',2,(/u1-l1,u2-l2+1/),kndlog)

!---source
SELECT CASE (intsrce)
   CASE (1); sflag3 = MERGE(1.0,0.0,nodeatuw(l1:u1,l2:u2,l3:u3).GT.0)
   CASE (2); sflag3 = MERGE(1.0,0.0,nodeatuw(l1:u1,l2:u2,l3:u3).EQ.2)
   CASE (3); sflag3 = MERGE(1.0,0.0,nodeatuw(l1:u1,l2:u2,l3:u3).GT.1)
   CASE (4); sflag3 = MERGE(1.0,0.0,nodeatuw(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                  & nodeatuw(l1:u1,l2:u2,l3:u3).EQ.2)
END SELECT

!---destination
fdest2 = nodeatc(l1:u1-1,l2:u2).GT.0

!
!4. Weight factors
!-----------------
!
!---allocate
ALLOCATE (wsum3(l1:u1-1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('wsum3',3,(/u1-l1,u2-l2+1,u3-l3+1/),kndrtype)

!---define
wsum3 = sflag3(l1:u1-1,:,:) + sflag3(l1+1:u1,:,:)

!
!5. Interpolate
!--------------
!

l_510: DO l=1,nosize
k_510: DO k=l3,u3
   WHERE (fdest2.AND.wsum3(:,:,k).GT.0.0)
      xout(:,:,k,l) = (sflag3(l1:u1-1,:,k)*xin(l1:u1-1,:,k,l)+&
                     & sflag3(l1+1:u1,:,k)*xin(l1+1:u1,:,k,l))/wsum3(:,:,k)
   ELSEWHERE
      xout(:,:,k,l) = flagval
   END WHERE
ENDDO k_510
ENDDO l_510

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag3,fdest2,wsum3)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE UWarr_at_W

!========================================================================

SUBROUTINE Varr_at_C(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                   & info,outflag)
!************************************************************************
!
! *Varr_at_C* Interpolate an array from V- to C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.6
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!


LOGICAL, INTENT(IN) :: info
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2)-1,&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at V-nodes
!*xout*      REAL    Output array at C-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => wet points only
!*lbounds*   INTEGER Lower bounds of input array
!*ubounds*   INTEGER Upper bounds of input array
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Varr_at_C'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ELSE
   ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ENDIF
ALLOCATE (fdest2(l1:u1,l2:u2-1),STAT=errstat)
CALL error_alloc('fdest2',2,(/u1-l1+1,u2-l2/),kndlog)

!---source
IF (iopt_arrint_3D.EQ.0) THEN
   SELECT CASE (intsrce)
      CASE (0); sflag2 = 1.0
      CASE (1); sflag2 = MERGE(1.0,0.0,node2dv(l1:u1,l2:u2).GT.0)
      CASE (2); sflag2 = MERGE(1.0,0.0,node2dv(l1:u1,l2:u2).EQ.2)
      CASE (3); sflag2 = MERGE(1.0,0.0,node2dv(l1:u1,l2:u2).GT.1)
      CASE (4); sflag2 = MERGE(1.0,0.0,node2dv(l1:u1,l2:u2).EQ.1.OR.&
                                     & node2dv(l1:u1,l2:u2).EQ.2)
   END SELECT
ELSE
   SELECT CASE (intsrce)
      CASE (0); sflag3 = 1.0
      CASE (1); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).GT.0)
      CASE (2); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).EQ.2)
      CASE (3); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).GT.1)
      CASE (4); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                  & nodeatv(l1:u1,l2:u2,l3:u3).EQ.2)
   END SELECT
ENDIF

!---destination
SELECT CASE (intdest)
   CASE (0); fdest2 = .TRUE.
   CASE (1); fdest2 = nodeatc(l1:u1,l2:u2-1).GT.0
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (wsum2(l1:u1,l2:u2-1),STAT=errstat)
   CALL error_alloc('wsum2',2,(/u1-l1+1,u2-l2/),kndrtype)
ELSE
   ALLOCATE (wsum3(l1:u1,l2:u2-1,l3:u3),STAT=errstat)
   CALL error_alloc('wsum3',3,(/u1-l1+1,u2-l2,u3-l3+1/),kndrtype)
ENDIF

!---define
IF (iopt_arrint_3D.EQ.0) THEN
   wsum2 = sflag2(:,l2:u2-1) + sflag2(:,l2+1:u2)
ELSE
   wsum3 = sflag3(:,l2:u2-1,:) + sflag3(:,l2+1:u2,:)
ENDIF

!
!5. Interpolate
!--------------
!

IF (iopt_arrint_3D.EQ.0) THEN

   l_501: DO l=1,nosize
   k_501: DO k=l3,u3
      WHERE (fdest2.AND.wsum2.GT.0.0)
         xout(:,:,k,l) = (sflag2(:,l2:u2-1)*xin(:,l2:u2-1,k,l)+&
                        & sflag2(:,l2+1:u2)*xin(:,l2+1:u2,k,l))/wsum2
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_501
   ENDDO l_501

ELSE

   l_502: DO l=1,nosize
   k_502: DO k=l3,u3
      WHERE (fdest2.AND.wsum3(:,:,k).GT.0.0)
         xout(:,:,k,l) = (sflag3(:,l2:u2-1,k)*xin(:,l2:u2-1,k,l)+&
                        & sflag3(:,l2+1:u2,k)*xin(:,l2+1:u2,k,l))/wsum3(:,:,k)
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_502
   ENDDO l_502

ENDIF

!
!6. Deallocate
!-------------
!

DEALLOCATE (fdest2)        
IF (iopt_arrint_3D.EQ.0) THEN
   DEALLOCATE (sflag2,wsum2)
ELSE
   DEALLOCATE (sflag3,wsum3)
ENDIF

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Varr_at_C

!========================================================================

SUBROUTINE Varr_at_U(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                   & info,outflag,hregular)
!************************************************************************
!
! *Varr_at_U* Interpolate an array from V- to U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.6
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1)+1:ubounds(1),lbounds(2):ubounds(2)-1,&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at V-nodes
!*xout*      REAL    Output array at U-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array
!*ubounds*   INTEGER Upper bounds of input array
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Varr_at_U'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXYreg
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
   ALLOCATE (fdest2(l1+1:u1,l2:u2-1),STAT=errstat)
   CALL error_alloc('fdest2',2,(/u1-l1,u2-l2/),kndlog)
ELSE
   ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
   ALLOCATE (fdest3(l1+1:u1,l2:u2-1,l3:u3),STAT=errstat)
   CALL error_alloc('fdest3',3,(/u1-l1,u2-l2,u3-l3+1/),kndlog)
ENDIF

!---source
IF (iopt_arrint_3D.EQ.0) THEN
   SELECT CASE (intsrce)
      CASE (0); sflag2 = 1.0
      CASE (1); sflag2 = MERGE(1.0,0.0,node2dv(l1:u1,l2:u2).GT.0)
      CASE (2); sflag2 = MERGE(1.0,0.0,node2dv(l1:u1,l2:u2).EQ.2)
      CASE (3); sflag2 = MERGE(1.0,0.0,node2dv(l1:u1,l2:u2).GT.1)
      CASE (4); sflag2 = MERGE(1.0,0.0,node2dv(l1:u1,l2:u2).EQ.1.OR.&
                                     & node2dv(l1:u1,l2:u2).EQ.2)
   END SELECT
ELSE
   SELECT CASE (intsrce)
      CASE (0); sflag3 = 1.0
      CASE (1); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).GT.0)
      CASE (2); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).EQ.2)
      CASE (3); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).GT.1)
      CASE (4); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                   & nodeatv(l1:u1,l2:u2,l3:u3).EQ.2)
   END SELECT
ENDIF

!---destination
IF (iopt_arrint_3D.EQ.0) THEN
   SELECT CASE (intdest)
      CASE (0); fdest2 = .TRUE.
      CASE (1); fdest2 = node2du(l1+1:u1,l2:u2-1).GT.0
      CASE (2); fdest2 = node2du(l1+1:u1,l2:u2-1).EQ.2
      CASE (3); fdest2 = node2du(l1+1:u1,l2:u2-1).GT.1
      CASE (4); fdest2 = node2du(l1+1:u1,l2:u2-1).EQ.1.OR.&
                       & node2du(l1+1:u1,l2:u2-1).EQ.2
   END SELECT
ELSE
   SELECT CASE (intdest)
      CASE (0); fdest3 = .TRUE.
      CASE (1); fdest3 = nodeatu(l1+1:u1,l2:u2-1,l3:u3).GT.0
      CASE (2); fdest3 = nodeatu(l1+1:u1,l2:u2-1,l3:u3).EQ.2
      CASE (3); fdest3 = nodeatu(l1+1:u1,l2:u2-1,l3:u3).GT.1
      CASE (4); fdest3 = nodeatu(l1+1:u1,l2:u2-1,l3:u3).EQ.1.OR.&
                    & nodeatu(l1+1:u1,l2:u2-1,l3:u3).EQ.2
   END SELECT
ENDIF

!
!4. Weight factors
!-----------------
!
!---allocate
IF (.NOT.hreg) THEN
   ALLOCATE (weights2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('weights2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ENDIF 

IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (wsum2(l1+1:u1,l2:u2-1),STAT=errstat)
   CALL error_alloc('wsum2',2,(/u1-l1,u2-l2/),kndrtype)
ELSE
   ALLOCATE (wsum3(l1+1:u1,l2:u2-1,l3:u3),STAT=errstat)
   CALL error_alloc('wsum3',3,(/u1-l1,u2-l2,u3-l3+1/),kndrtype)
ENDIF

!---define
IF (iopt_arrint_3D.EQ.0) THEN
   IF (hreg) THEN
      wsum2 = sflag2(l1+1:u1,l2+1:u2) + sflag2(l1:u1-1,l2+1:u2) + &
            & sflag2(l1+1:u1,l2:u2-1) + sflag2(l1:u1-1,l2:u2-1)
   ELSE
      IF (dYregY.AND.dYregX) THEN 
         weights2 = delxatu(l1:u1,l2:u2)
      ELSEIF (dXregY.AND.dXregX) THEN
         weights2 = delyatv(l1:u1,l2:u2)
      ELSE
         weights2 = delxatu(l1:u1,l2:u2)*delyatv(l1:u1,l2:u2)
      ENDIF
      wsum2 = sflag2(l1+1:u1,l2+1:u2)*weights2(l1:u1-1,l2:u2-1) + &
            & sflag2(l1:u1-1,l2+1:u2)*weights2(l1+1:u1,l2:u2-1) + &
            & sflag2(l1+1:u1,l2:u2-1)*weights2(l1:u1-1,l2+1:u2) + &
            & sflag2(l1:u1-1,l2:u2-1)*weights2(l1+1:u1,l2+1:u2)
   ENDIF
   WHERE (wsum2.LE.0.0) fdest2 = .FALSE.
ELSE
   IF (hreg) THEN
      wsum3 = sflag3(l1+1:u1,l2+1:u2,:) + sflag3(l1:u1-1,l2+1:u2,:) + &
            & sflag3(l1+1:u1,l2:u2-1,:) + sflag3(l1:u1-1,l2:u2-1,:)
   ELSE
      IF (dYregY.AND.dYregX) THEN 
         weights2 = delxatu(l1:u1,l2:u2)
      ELSEIF (dXregY.AND.dXregX) THEN
         weights2 = delyatv(l1:u1,l2:u2)
      ELSE
         weights2 = delxatu(l1:u1,l2:u2)*delyatv(l1:u1,l2:u2)
      ENDIF
      k_410: DO k=l3,u3
         wsum3(:,:,k) = sflag3(l1+1:u1,l2+1:u2,k)*weights2(l1:u1-1,l2:u2-1) + &
                      & sflag3(l1:u1-1,l2+1:u2,k)*weights2(l1+1:u1,l2:u2-1) + &
                      & sflag3(l1+1:u1,l2:u2-1,k)*weights2(l1:u1-1,l2+1:u2) + &
                      & sflag3(l1:u1-1,l2:u2-1,k)*weights2(l1+1:u1,l2+1:u2)
      ENDDO k_410
   ENDIF
   WHERE (wsum3.LE.0.0) fdest3 = .FALSE.
ENDIF

!
!5. Interpolate
!--------------
!
!5.1 Using 2-D mask
!------------------
!

IF (iopt_arrint_3D.EQ.0) THEN

!
!5.1.1 Without weights
!-------------------
!

   IF (hreg) THEN
      l_511: DO l=1,nosize
      k_511: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = &
                       &(sflag2(l1:u1-1,l2:u2-1)*xin(l1:u1-1,l2:u2-1,k,l)+&
                       & sflag2(l1+1:u1,l2:u2-1)*xin(l1+1:u1,l2:u2-1,k,l)+&
                       & sflag2(l1:u1-1,l2+1:u2)*xin(l1:u1-1,l2+1:u2,k,l)+&
                       & sflag2(l1+1:u1,l2+1:u2)*xin(l1+1:u1,l2+1:u2,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_511
      ENDDO l_511

!
!5.1.2 With weights
!----------------
!

   ELSE
      l_512: DO l=1,nosize
      k_512: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = &
                &(sflag2(l1:u1-1,l2:u2-1)*weights2(l1+1:u1,l2+1:u2)*&
                & xin(l1:u1-1,l2:u2-1,k,l)+&
                & sflag2(l1+1:u1,l2:u2-1)*weights2(l1:u1-1,l2+1:u2)*&
                & xin(l1+1:u1,l2:u2-1,k,l)+&
                & sflag2(l1:u1-1,l2+1:u2)*weights2(l1+1:u1,l2:u2-1)*&
                & xin(l1:u1-1,l2+1:u2,k,l)+&
                & sflag2(l1+1:u1,l2+1:u2)*weights2(l1:u1-1,l2:u2-1)*&
                & xin(l1+1:u1,l2+1:u2,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_512
      ENDDO l_512
   ENDIF

!
!5.2 Using 3-D mask
!------------------
!

ELSE

!
!5.2.1 Without weights
!---------------------
!

   IF (hreg) THEN
      l_521: DO l=1,nosize
         WHERE (fdest3)
            xout(:,:,:,l) = &
                     &(sflag3(l1:u1-1,l2:u2-1,:)*xin(l1:u1-1,l2:u2-1,:,l)+&
                     & sflag3(l1+1:u1,l2:u2-1,:)*xin(l1+1:u1,l2:u2-1,:,l)+&
                     & sflag3(l1:u1-1,l2+1:u2,:)*xin(l1:u1-1,l2+1:u2,:,l)+&
                     & sflag3(l1+1:u1,l2+1:u2,:)*xin(l1+1:u1,l2+1:u2,:,l))/wsum3
         ELSEWHERE
            xout(:,:,:,l) = flagval
         END WHERE
      ENDDO l_521
      
!
!5.2.2 With weights
!------------------
!

   ELSE
      l_522: DO l=1,nosize
      k_522: DO k=l3,u3
         WHERE (fdest3(:,:,k))
            xout(:,:,k,l) = &
                &(sflag3(l1:u1-1,l2:u2-1,k)*weights2(l1+1:u1,l2+1:u2)*&
                & xin(l1:u1-1,l2:u2-1,k,l)+&
                & sflag3(l1+1:u1,l2:u2-1,k)*weights2(l1:u1-1,l2+1:u2)*&
                & xin(l1+1:u1,l2:u2-1,k,l)+&
                & sflag3(l1:u1-1,l2+1:u2,k)*weights2(l1+1:u1,l2:u2-1)*&
                & xin(l1:u1-1,l2+1:u2,k,l)+&
                & sflag3(l1+1:u1,l2+1:u2,k)*weights2(l1:u1-1,l2:u2-1)*&
                & xin(l1+1:u1,l2+1:u2,k,l))/wsum3(:,:,k)
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_522
      ENDDO l_522

   ENDIF

ENDIF

!
!6. Deallocate
!-------------
!        
IF (iopt_arrint_3D.EQ.0) THEN
   DEALLOCATE (sflag2,fdest2,wsum2)
ELSE
   DEALLOCATE (sflag3,fdest3,wsum3)
ENDIF
IF (.NOT.hreg) DEALLOCATE (weights2)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Varr_at_U

!========================================================================

SUBROUTINE Varr_at_UV(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                    & info,outflag,hregular)
!************************************************************************
!
! *Varr_at_UV* Interpolate an array from V- to UV-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.6
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1)+1:ubounds(1),lbounds(2):ubounds(2),&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at V-nodes
!*xout*      REAL    Output array at UV-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => interior points and open boundaries only
!                    = 2 => interior points and UV-node open boundaries only
!*lbounds*   INTEGER Lower bounds of input array
!*ubounds*   INTEGER Upper bounds of input array
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Varr_at_UV'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXregX
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
   ALLOCATE (fdest2(l1+1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('fdest2',2,(/u1-l1,u2-l2+1/),kndlog)
ELSE
   ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
   ALLOCATE (fdest3(l1+1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('fdest3',3,(/u1-l1,u2-l2+1,u3-l3+1/),kndlog)
ENDIF

!---source
IF (iopt_arrint_3D.EQ.0) THEN

   SELECT CASE (intsrce)
      CASE (0); sflag2 = 1.0
      CASE (1); sflag2 = MERGE(1.0,0.0,node2dv(l1:u1,l2:u2).GT.0)
      CASE (2); sflag2 = MERGE(1.0,0.0,node2dv(l1:u1,l2:u2).EQ.2)
      CASE (3); sflag2 = MERGE(1.0,0.0,node2dv(l1:u1,l2:u2).GT.1)
      CASE (4); sflag2 = MERGE(1.0,0.0,node2dv(l1:u1,l2:u2).EQ.1.OR.&
                                     & node2dv(l1:u1,l2:u2).EQ.2)
   END SELECT
ELSE
   SELECT CASE (intsrce)
      CASE (0); sflag3 = 1.0
      CASE (1); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).GT.0)
      CASE (2); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).EQ.2)
      CASE (3); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).GT.1)
      CASE (4); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                     & nodeatv(l1:u1,l2:u2,l3:u3).EQ.2)
   END SELECT
ENDIF

!---destination
IF (iopt_arrint_3D.EQ.0) THEN
   SELECT CASE (intdest)
      CASE (0); fdest2 = .TRUE.
      CASE (1); fdest2 = node2duv(l1+1:u1,l2:u2).GT.0
      CASE (2); fdest2 = node2duv(l1+1:u1,l2:u2).EQ.1.OR.&
                       & node2duv(l1+1:u1,l2:u2).EQ.2
   END SELECT
ELSE
   SELECT CASE (intdest)
      CASE (0); fdest3 = .TRUE.
      CASE (1); fdest3 = nodeatuv(l1+1:u1,l2:u2,l3:u3).GT.0
      CASE (2); fdest3 = nodeatuv(l1+1:u1,l2:u2,l3:u3).EQ.1.OR.&
                       & nodeatuv(l1+1:u1,l2:u2,l3:u3).EQ.2
   END SELECT
ENDIF

!
!4. Weight factors
!-----------------
!

IF (.NOT.hreg) THEN
   ALLOCATE (weights2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('weights2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ENDIF

IF (iopt_arrint_3D.EQ.0) THEN
   ALLOCATE (wsum2(l1+1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('wsum2',2,(/u1-l1,u2-l2+1/),kndrtype)
ELSE
   ALLOCATE (wsum3(l1+1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('wsum3',3,(/u1-l1,u2-l2+1,u3-l3+1/),kndrtype)
ENDIF

!---define
IF (iopt_arrint_3D.EQ.0) THEN
   IF (hreg) THEN
      wsum2 = sflag2(l1+1:u1,:) + sflag2(l1:u1-1,:)
   ELSE
      weights2 = delxatv(l1:u1,l2:u2)
      wsum2 = sflag2(l1+1:u1,:)*weights2(l1:u1-1,:) + &
            & sflag2(l1:u1-1,:)*weights2(l1+1:u1,:)
   ENDIF
   WHERE (wsum2.LE.0.0) fdest2 = .FALSE.
ELSE
   IF (hreg) THEN
      wsum3 = sflag3(l1+1:u1,:,:) + sflag3(l1:u1-1,:,:)
   ELSE
      weights2 = delxatv(l1:u1,l2:u2)
      k_410: DO k=l3,u3
         wsum3(:,:,k) = sflag3(l1+1:u1,:,k)*weights2(l1:u1-1,:) + &
                      & sflag3(l1:u1-1,:,k)*weights2(l1+1:u1,:)
      ENDDO k_410
   ENDIF
   WHERE (wsum3.LE.0.0) fdest3 = .FALSE.
ENDIF

!
!5. Interpolate
!--------------
!
!5.1 Using 2-D mask
!------------------
!

IF (iopt_arrint_3D.EQ.0) THEN

!
!5.1.1 Without weights
!-------------------
!

   IF (hreg) THEN
      l_511: DO l=1,nosize
      k_511: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = (sflag2(l1:u1-1,:)*xin(l1:u1-1,:,k,l)+&
                           & sflag2(l1+1:u1,:)*xin(l1+1:u1,:,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_511
      ENDDO l_511

!
!5.1.2 With weights
!------------------
!

   ELSE

      l_512: DO l=1,nosize
      k_521: DO k=l3,u3
         WHERE (fdest2)
            xout(:,:,k,l) = &
              &(sflag2(l1:u1-1,:)*weights2(l1+1:u1,:)*xin(l1:u1-1,:,k,l)+&
              & sflag2(l1+1:u1,:)*weights2(l1:u1-1,:)*xin(l1+1:u1,:,k,l))/wsum2
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_521
      ENDDO l_512
   ENDIF

!
!5.2 Using 3-D mask
!------------------
!

ELSE

!
!5.2.1 Without weights
!---------------------
!

   IF (hreg) THEN
      l_521: DO l=1,nosize
         WHERE (fdest3)
            xout(:,:,:,l) = (sflag3(l1:u1-1,:,:)*xin(l1:u1-1,:,:,l)+&
                           & sflag3(l1+1:u1,:,:)*xin(l1+1:u1,:,:,l))/wsum3
         ELSEWHERE
            xout(:,:,:,l) = flagval
         END WHERE
      ENDDO l_521

!
!5.2.2 With weights
!------------------
!

   ELSE

      l_522: DO l=1,nosize
      k_522: DO k=l3,u3
         WHERE (fdest3(:,:,k))
            xout(:,:,k,l) = &
               &(sflag3(l1:u1-1,:,k)*weights2(l1+1:u1,:)*xin(l1:u1-1,:,k,l)+&
               & sflag3(l1+1:u1,:,k)*weights2(l1:u1-1,:)*xin(l1+1:u1,:,k,l))&
               & /wsum3(:,:,k)
         ELSEWHERE
            xout(:,:,k,l) = flagval
         END WHERE
      ENDDO k_522
      ENDDO l_522
   ENDIF

ENDIF

!
!6. Deallocate
!-------------
!        
IF (iopt_arrint_3D.EQ.0) THEN
   DEALLOCATE (sflag2,fdest2,wsum2)
ELSE
   DEALLOCATE (sflag3,fdest3,wsum3)
ENDIF
IF (.NOT.hreg) DEALLOCATE (weights2)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Varr_at_UV

!========================================================================

SUBROUTINE Varr_at_VW(xin,xout,intsrce,intdest,lbounds,ubounds,nosize, &
                    & ivarid,info,outflag,vregular)
!************************************************************************
!
! *Varr_at_VW* Interpolate an array from V- to VW-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: vregular
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                           & lbounds(3)+1:ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at V-nodes
!*xout*      REAL    Output array at VW-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*vregular*  LOGICAL Flag to select uniform or area averaging in the vertical
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Varr_at_VW'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(vregular)) THEN
   vreg = iopt_arrint_vreg.EQ.0.OR.dZregZ
ELSE
   vreg = vregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ALLOCATE (fdest3(l1:u1,l2:u2,l3+1:u3),STAT=errstat)
CALL error_alloc('fdest3',3,(/u1-l1+1,u2-l2+1,u3-l3/),kndlog)

!---source
SELECT CASE (intsrce)
   CASE (1); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).GT.0)
   CASE (2); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).EQ.2)
   CASE (3); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).GT.1)
   CASE (4); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                  & nodeatv(l1:u1,l2:u2,l3:u3).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (1); fdest3 = nodeatvw(l1:u1,l2:u2,l3+1:u3).GT.0
   CASE (2); fdest3 = nodeatvw(l1:u1,l2:u2,l3+1:u3).EQ.2
   CASE (3); fdest3 = nodeatvw(l1:u1,l2:u2,l3+1:u3).GT.1
   CASE (4); fdest3 = nodeatvw(l1:u1,l2:u2,l3+1:u3).EQ.1.OR.&
                    & nodeatvw(l1:u1,l2:u2,l3+1:u3).EQ.2
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
IF (.NOT.vreg) THEN
   ALLOCATE (weights3(l1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('weights3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ENDIF
ALLOCATE (wsum3(l1:u1,l2:u2,l3+1:u3),STAT=errstat)
CALL error_alloc('wsum3',3,(/u1-l1+1,u2-l2+1,u3-l3/),kndrtype)

!---define
IF (vreg) THEN
   wsum3 = sflag3(:,:,l3:u3-1) + sflag3(:,:,l3+1:u3)
ELSE
   weights3 = delzatv(l1:u1,l2:u2,l3:u3)
   wsum3 = sflag3(:,:,l3+1:u3)*weights3(:,:,l3:u3-1) + &
         & sflag3(:,:,l3:u3-1)*weights3(:,:,l3+1:u3)
ENDIF

WHERE (wsum3.LE.0.0) fdest3 = .FALSE.

!
!5. Interpolate
!--------------
!
!5.1 Without weights
!-------------------
!

IF (vreg) THEN

   l_510: DO l=1,nosize
      WHERE (fdest3)
         xout(:,:,:,l) = sflag3(:,:,l3:u3-1)*xin(:,:,l3:u3-1,l) + &
                       & sflag3(:,:,l3+1:u3)*xin(:,:,l3+1:u3,l)
      ELSEWHERE
         xout(:,:,:,l) = flagval
      END WHERE
   ENDDO l_510

!
!5.2 With weights
!----------------
!

ELSE

   l_520: DO l=1,nosize
      WHERE (fdest3)
         xout(:,:,:,l) = &
               & (sflag3(:,:,l3:u3-1)*weights3(:,:,l3+1:u3)*xin(:,:,l3:u3-1,l)+&
               &  sflag3(:,:,l3+1:u3)*weights3(:,:,l3:u3-1)*xin(:,:,l3+1:u3,l))&
               & /wsum3
      ELSEWHERE
         xout(:,:,:,l) = flagval
      END WHERE
   ENDDO l_520

ENDIF

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag3,fdest3,wsum3)
IF (.NOT.vreg) DEALLOCATE (weights3)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Varr_at_VW

!========================================================================

SUBROUTINE Varr_at_W(xin,xout,intsrce,lbounds,ubounds,nosize,ivarid,info,&
                   & outflag,vregular)
!************************************************************************
!
! *Varr_at_W* Interpolate an array from V- to W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: vregular
INTEGER, INTENT(IN) :: intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2)-1,&
                           & lbounds(3)+1:ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at V-nodes
!*xout*      REAL    Output array at W-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*vregular*  LOGICAL Flag to select uniform or area averaging in the vertical
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Varr_at_W'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(vregular)) THEN
   vreg = iopt_arrint_vreg.EQ.0
ELSE
   vreg = vregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ALLOCATE (fdest2(l1:u1,l2:u2-1),STAT=errstat)
CALL error_alloc('fdest3',2,(/u1-l1+1,u2-l2/),kndlog)

!---source
SELECT CASE (intsrce)
   CASE (1); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).GT.0)
   CASE (2); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).EQ.2)
   CASE (3); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).GT.1)
   CASE (4); sflag3 = MERGE(1.0,0.0,nodeatv(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                  & nodeatv(l1:u1,l2:u2,l3:u3).EQ.2)
END SELECT

!---destination
fdest2 = nodeatc(l1:u1,l2:u2-1).GT.0

!
!4. Weight factors
!-----------------
!
!---allocate
IF (.NOT.vreg) THEN
   ALLOCATE (weights3(l1:u1,l2:u2,l3:u3),STAT=errstat)
   CALL error_alloc('weights3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ENDIF
ALLOCATE (wsum3(l1:u1,l2:u2-1,l3+1:u3),STAT=errstat)
CALL error_alloc('wsum3',3,(/u1-l1+1,u2-l2,u3-l3/),kndrtype)

!---define
IF (vreg) THEN
   wsum3 = sflag3(:,l2:u2-1,l3:u3-1)+sflag3(:,l2+1:u2,l3:u3-1)+&
         & sflag3(:,l2:u2-1,l3+1:u3)+sflag3(:,l2+1:u2,l3+1:u3)
   WHERE (wsum3.LE.0.0) fdest3 = .FALSE.
ELSE
   k_410: DO k=l3,u3
      weights3(:,:,k) = delzatv(l1:u1,l2:u2,k)
      WHERE (nodeatv(l1:u1,l2+1:u2,k).GE.2.AND.nodeatc(l1:u1,l2:u2-1).GT.0)
         weights3(:,l2+1:u2,k) = delzatc(l1:u1,l2:u2-1,k)
      END WHERE
      WHERE (nodeatv(l1:u1,l2:u2-1,k).GE.2.AND.nodeatc(l1:u1,l2:u2-1).GT.0)
         weights3(:,l2:u2-1,k) = delzatc(l1:u1,l2:u2-1,k)
      END WHERE
   ENDDO k_410
   wsum3(:,:,l3+1:u3) = &
                     & sflag3(:,l2:u2-1,l3:u3-1)*weights3(:,l2+1:u2,l3+1:u3) + &
                     & sflag3(:,l2+1:u2,l3:u3-1)*weights3(:,l2:u2-1,l3+1:u3) + &
                     & sflag3(:,l2:u2-1,l3+1:u3)*weights3(:,l2+1:u2,l3:u3-1) + &
                     & sflag3(:,l2+1:u2,l3+1:u3)*weights3(:,l2:u2-1,l3:u3-1)
ENDIF

!
!5. Interpolate
!--------------
!
!5.1 Without weights
!-------------------
!

IF (vreg) THEN

   l_510: DO l=1,nosize
   k_510: DO k=l3+1,u3
      WHERE (fdest2.AND.wsum3(:,:,k).GT.0.0)
         xout(:,:,k,l) = (sflag3(:,l2:u2-1,k-1)*xin(:,l2:u2-1,k-1,l)+&
                        & sflag3(:,l2:u2-1,k)*xin(:,l2:u2-1,k,l)+&
                        & sflag3(:,l2+1:u2,k-1)*xin(:,l2+1:u2,k-1,l)+&
                        & sflag3(:,l2+1:u2,k)*xin(:,l2+1:u2,k,l))&
                        & /wsum3(:,:,k)
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_510
   ENDDO l_510

!
!5.2 With weights
!----------------
!

ELSE

   l_520: DO l=1,nosize
   k_520: DO k=l3+1,u3
      WHERE (fdest2.AND.wsum3(:,:,k).GT.0.0)
         xout(:,:,k,l) = &
           & (sflag3(:,l2:u2-1,k-1)*weights3(:,l2+1:u2,k)*xin(:,l2:u2-1,k-1,l)+&
           &  sflag3(:,l2+1:u2,k-1)*weights3(:,l2:u2-1,k)*xin(:,l2+1:u2,k-1,l)+&
           &  sflag3(:,l2:u2-1,k)*weights3(:,l2+1:u2,k-1)*xin(:,l2:u2-1,k,l)+&
           &  sflag3(:,l2+1:u2,k)*weights3(:,l2:u2-1,k-1)*xin(:,l2+1:u2,k,l))&
           & /wsum3(:,:,k)
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_520
   ENDDO l_520

ENDIF

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag3,fdest2,wsum3)
IF (.NOT.vreg) DEALLOCATE (weights3)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Varr_at_W

!========================================================================

SUBROUTINE VWarr_at_V(xin,xout,intsrce,intdest,lbounds,ubounds,nosize,ivarid,&
                    & info,outflag)
!************************************************************************
!
! *VWarr_at_V* Interpolate an array from VW- to V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                           & lbounds(3):ubounds(3)-1,nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at VW-nodes
!*xout*      REAL    Output array at V-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'VWarr_at_V'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ALLOCATE (fdest3(l1:u1,l2:u2,l3:u3-1),STAT=errstat)
CALL error_alloc('fdest3',3,(/u1-l1+1,u2-l2+1,u3-l3/),kndlog)

!---source
SELECT CASE (intsrce)
   CASE (1); sflag3 = MERGE(1.0,0.0,nodeatvw(l1:u1,l2:u2,l3:u3).GT.0)
   CASE (2); sflag3 = MERGE(1.0,0.0,nodeatvw(l1:u1,l2:u2,l3:u3).EQ.2)
   CASE (3); sflag3 = MERGE(1.0,0.0,nodeatvw(l1:u1,l2:u2,l3:u3).GT.1)
   CASE (4); sflag3 = MERGE(1.0,0.0,nodeatvw(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                  & nodeatvw(l1:u1,l2:u2,l3:u3).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (1); fdest3 = nodeatv(l1:u1,l2:u2,l3:u3-1).GT.0
   CASE (2); fdest3 = nodeatv(l1:u1,l2:u2,l3:u3-1).EQ.2
   CASE (3); fdest3 = nodeatv(l1:u1,l2:u2,l3:u3-1).GT.1
   CASE (4); fdest3 = nodeatv(l1:u1,l2:u2,l3:u3-1).EQ.1.OR.&
                   & nodeatv(l1:u1,l2:u2,l3:u3-1).EQ.2
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
ALLOCATE (wsum3(l1:u1,l2:u2,l3:u3-1),STAT=errstat)
CALL error_alloc('wsum3',3,(/u1-l1+1,u2-l2+1,u3-l3/),kndrtype)

!---define
wsum3 = sflag3(:,:,l3:u3-1) + sflag3(:,:,l3+1:u3)

WHERE (wsum3.LE.0.0) fdest3 = .FALSE.

!
!5. Interpolate
!--------------
!

l_410: DO l=1,nosize
   WHERE (fdest3)
      xout(:,:,:,l) = (sflag3(:,:,l3:u3-1)*xin(:,:,l3:u3-1,l)+&
                     & sflag3(:,:,l3+1:u3)*xin(:,:,l3+1:u3,l))/wsum3
   ELSEWHERE
      xout(:,:,:,l) = flagval
   END WHERE
ENDDO l_410

!
!6. Deallocate
!-------------
!

DEALLOCATE (sflag3,fdest3,wsum3)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE VWarr_at_V

!========================================================================

SUBROUTINE VWarr_at_W(xin,xout,intsrce,lbounds,ubounds,nosize,ivarid,info,&
                    & outflag)
!************************************************************************
!
! *VWarr_at_W* Interpolate an array from VW- to W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.11.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
INTEGER, INTENT(IN) :: intsrce, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2)-1,&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at VW-nodes
!*xout*      REAL    Output array at W-nodes
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'VWarr_at_W'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag3(l1:u1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('sflag3',3,(/u1-l1+1,u2-l2+1,u3-l3+1/),kndrtype)
ALLOCATE (fdest2(l1:u1,l2:u2-1),STAT=errstat)
CALL error_alloc('fdest2',2,(/u1-l1+1,u2-l2/),kndlog)

!---source
SELECT CASE (intsrce)
   CASE (1); sflag3 = MERGE(1.0,0.0,nodeatvw(l1:u1,l2:u2,l3:u3).GT.0)
   CASE (2); sflag3 = MERGE(1.0,0.0,nodeatvw(l1:u1,l2:u2,l3:u3).EQ.2)
   CASE (3); sflag3 = MERGE(1.0,0.0,nodeatvw(l1:u1,l2:u2,l3:u3).GT.1)
   CASE (4); sflag3 = MERGE(1.0,0.0,nodeatvw(l1:u1,l2:u2,l3:u3).EQ.1.OR.&
                                  & nodeatvw(l1:u1,l2:u2,l3:u3).EQ.2)
END SELECT

!---destination
fdest2 = nodeatc(l1:u1,l2:u2-1).GT.0

!
!4. Weight factors
!-----------------
!
!---allocate
ALLOCATE (wsum3(l1:u1,l2:u2-1,l3:u3),STAT=errstat)
CALL error_alloc('wsum3',3,(/u1-l1+1,u2-l2,u3-l3+1/),kndrtype)

!---define
wsum3 = sflag3(:,l2:u2-1,:) + sflag3(:,l2+1:u2,:)

!
!5. Interpolate
!--------------
!

l_510: DO l=1,nosize
k_510: DO k=l3,u3
   WHERE (fdest2.AND.wsum3(:,:,k).GT.0.0)
      xout(:,:,k,l) = (sflag3(:,l2:u2-1,k)*xin(:,l2:u2-1,k,l)+&
                     & sflag3(:,l2+1:u2,k)*xin(:,l2+1:u2,k,l))/wsum3(:,:,k)
   ELSEWHERE
      xout(:,:,k,l) = flagval
   END WHERE
ENDDO k_510
ENDDO l_510

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag3,fdest2,wsum3)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE VWarr_at_W

!========================================================================

SUBROUTINE Warr_at_C(xin,xout,lbounds,ubounds,nosize,ivarid,info,outflag)
!************************************************************************
!
! *Warr_at_C* Interpolate an array from W- to C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.0
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
INTEGER, INTENT(IN) :: ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                           & lbounds(3):ubounds(3)-1,nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at W-nodes
!*xout*      REAL    Output array at C-nodes
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Warr_at_C'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Interpolate
!--------------
!

l_310: DO l=1,nosize
k_310: DO k=l3,u3-1
   WHERE (nodeatc(l1:u1,l2:u2).GT.0)
      xout(:,:,k,l) = 0.5*(xin(:,:,k,l)+xin(:,:,k+1,l))
   ELSEWHERE
      xout(:,:,k,l) = flagval
   END WHERE
ENDDO k_310
ENDDO l_310

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Warr_at_C

!========================================================================

SUBROUTINE Warr_at_U(xin,xout,intdest,lbounds,ubounds,nosize,ivarid,info,&
                   & outflag,hregular)
!************************************************************************
!
! *Warr_at_U* Interpolate an array from W- to U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: intdest, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1)+1:ubounds(1),lbounds(2):ubounds(2),&
                           & lbounds(3):ubounds(3)-1,nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at W-nodes
!*xout*      REAL    Output array at U-nodes
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Warr_at_U'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXregX
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ALLOCATE (fdest3(l1+1:u1,l2:u2,l3:u3-1),STAT=errstat)
CALL error_alloc('fdest3',3,(/u1-l1,u2-l2+1,u3-l3/),kndlog)

!---source
sflag2 = MERGE(1.0,0.0,nodeatc(l1:u1,l2:u2).GT.0)

!---destination
SELECT CASE (intdest)
   CASE (1); fdest3 = nodeatu(l1+1:u1,l2:u2,l3:u3-1).GT.0
   CASE (2); fdest3 = nodeatu(l1+1:u1,l2:u2,l3:u3-1).EQ.2
   CASE (3); fdest3 = nodeatu(l1+1:u1,l2:u2,l3:u3-1).GT.1
   CASE (4); fdest3 = nodeatu(l1+1:u1,l2:u2,l3:u3-1).EQ.1.OR.&
                    & nodeatu(l1+1:u1,l2:u2,l3:u3-1).EQ.2
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
IF (.NOT.hreg) THEN
   ALLOCATE (weights2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('weights2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ENDIF
ALLOCATE (wsum2(l1+1:u1,l2:u2),STAT=errstat)
CALL error_alloc('wsum2',2,(/u1-l1,u2-l2+1/),kndrtype)

!---define
IF (hreg) THEN
   wsum2 = 2.0*(sflag2(l1+1:u1,:)+sflag2(l1:u1-1,:))
ELSE
   weights2 = delxatc(l1:u1,l2:u2)
   wsum2 = 2.0*(sflag2(l1+1:u1,:)*weights2(l1:u1-1,:) + &
              & sflag2(l1:u1-1,:)*weights2(l1+1:u1,:))
ENDIF

!
!5. Interpolate
!--------------
!
!5.1 Without weights
!-------------------
!

IF (hreg) THEN

   l_510: DO l=1,nosize
   k_510: DO k=l3,u3-1
      WHERE (fdest3(:,:,k).AND.wsum2.GT.0.0)
         xout(:,:,k,l) = &
              &(sflag2(l1:u1-1,:)*(xin(l1:u1-1,:,k,l)+xin(l1:u1-1,:,k+1,l))+&
              & sflag2(l1+1:u1,:)*(xin(l1+1:u1,:,k,l)+xin(l1+1:u1,:,k+1,l)))&
              & /wsum2
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_510
   ENDDO l_510

!
!5.2 With weights
!----------------
!

ELSE

   l_520: DO l=1,nosize
   k_520: DO k=l3,u3-1
      WHERE (fdest3(:,:,k).AND.wsum2.GT.0.0)
         xout(:,:,k,l) = (sflag2(l1:u1-1,:)*weights2(l1+1:u1,:)*&
                        &(xin(l1:u1-1,:,k,l)+xin(l1:u1-1,:,k+1,l))+&
                        &(sflag2(l1+1:u1,:)*weights2(l1:u1-1,:)*&
                        &(xin(l1+1:u1,:,k,l)+xin(l1+1:u1,:,k+1,l))))/wsum2
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_520
   ENDDO l_520

ENDIF

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag2,fdest3,wsum2)
IF (.NOT.hreg) DEALLOCATE (weights2)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Warr_at_U

!========================================================================

SUBROUTINE Warr_at_UW(xin,xout,intdest,lbounds,ubounds,nosize,ivarid,info,&
                    & outflag,hregular)
!************************************************************************
!
! *Warr_at_UW* Interpolate an array from W- to UW-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: intdest, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1)+1:ubounds(1),lbounds(2):ubounds(2),&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at W-nodes
!*xout*      REAL    Output array at UW-nodes
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Warr_at_UW'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXregX
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ALLOCATE (fdest3(l1+1:u1,l2:u2,l3:u3),STAT=errstat)
CALL error_alloc('fdest3',3,(/u1-l1,u2-l2+1,u3-l3+1/),kndlog)

!---source
sflag2 = MERGE(1.0,0.0,nodeatc(l1:u1,l2:u2).GT.0)

!---destination
SELECT CASE (intdest)
   CASE (1); fdest3 = nodeatuw(l1+1:u1,l2:u2,l3:u3).GT.0
   CASE (2); fdest3 = nodeatuw(l1+1:u1,l2:u2,l3:u3).EQ.2
   CASE (3); fdest3 = nodeatuw(l1+1:u1,l2:u2,l3:u3).GT.1
   CASE (4); fdest3 = nodeatuw(l1+1:u1,l2:u2,l3:u3).EQ.1.OR.&
                    & nodeatuw(l1+1:u1,l2:u2,l3:u3).EQ.2
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
IF (.NOT.hreg) THEN
   ALLOCATE (weights2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('weights2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ENDIF
ALLOCATE (wsum2(l1+1:u1,l2:u2),STAT=errstat)
CALL error_alloc('wsum2',2,(/u1-l1,u2-l2+1/),kndrtype)

!---define
IF (hreg) THEN
   wsum2 = sflag2(l1+1:u1,:) + sflag2(l1:u1-1,:)
ELSE
   weights2 = delxatc(l1:u1,l2:u2)
   wsum2 = sflag2(l1+1:u1,:)*weights2(l1:u1-1,:) + &
         & sflag2(l1:u1-1,:)*weights2(l1+1:u1,:)
ENDIF

k_410: DO k=l3,u3
   WHERE (wsum2.LE.0.0) fdest3(:,:,k) = .FALSE.
ENDDO k_410

!
!5. Interpolate
!--------------
!
!5.1 Without weights
!-------------------
!

IF (hreg) THEN

   l_510: DO l=1,nosize
   k_510: DO k=l3,u3
      WHERE (fdest3(:,:,k))
         xout(:,:,k,l) = (sflag2(l1:u1-1,:)*xin(l1:u1-1,:,k,l)+&
                          sflag2(l1+1:u1,:)*xin(l1+1:u1,:,k,l))/wsum2
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_510
   ENDDO l_510

!
!5.2 With weights
!----------------
!

ELSE

   l_520: DO l=1,nosize
   k_520: DO k=l3,u3
      WHERE (fdest3(:,:,k))
         xout(:,:,k,l) = &
            &(sflag2(l1:u1-1,:)*weights2(l1+1:u1,:)*xin(l1:u1-1,:,k,l)+&
            & sflag2(l1+1:u1,:)*weights2(l1:u1-1,:)*xin(l1+1:u1,:,k,l))/wsum2
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_520
   ENDDO l_520

ENDIF

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag2,fdest3,wsum2)
IF (.NOT.hreg) DEALLOCATE (weights2)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Warr_at_UW

!========================================================================

SUBROUTINE Warr_at_V(xin,xout,intdest,lbounds,ubounds,nosize,ivarid,info,&
                   & outflag,hregular)
!************************************************************************
!
! *Warr_at_V* Interpolate an array from W- to V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: intdest, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2)+1:ubounds(2),&
                           & lbounds(3):ubounds(3)-1,nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at W-nodes
!*xout*      REAL    Output array at V-nodes
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Warr_at_V'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dYregY
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ALLOCATE (fdest3(l1:u1,l2+1:u2,l3:u3-1),STAT=errstat)
CALL error_alloc('fdest3',3,(/u1-l1+1,u2-l2,u3-l3/),kndlog)

!---source
sflag2 = MERGE(1.0,0.0,nodeatc(l1:u1,l2:u2).GT.0)

!---destination
SELECT CASE (intdest)
   CASE (1); fdest3 = nodeatv(l1:u1,l2+1:u2,l3:u3-1).GT.0
   CASE (2); fdest3 = nodeatv(l1:u1,l2+1:u2,l3:u3-1).EQ.2
   CASE (3); fdest3 = nodeatv(l1:u1,l2+1:u2,l3:u3-1).GT.1
   CASE (4); fdest3 = nodeatv(l1:u1,l2+1:u2,l3:u3-1).EQ.1.OR.&
                    & nodeatv(l1:u1,l2+1:u2,l3:u3-1).EQ.2
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
IF (.NOT.hreg) THEN
   ALLOCATE (weights2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('weights2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ENDIF
ALLOCATE (wsum2(l1:u1,l2+1:u2),STAT=errstat)
CALL error_alloc('wsum2',2,(/u1-l1+1,u2-l2/),kndrtype)

!---define
IF (hreg) THEN
   wsum2 = 2.0*(sflag2(:,l2+1:u2)+sflag2(:,l2:u2-1))
ELSE
   weights2 = delyatc(l1:u1,l2:u2)
   wsum2 = 2.0*(sflag2(:,l2+1:u2)*weights2(:,l2:u2-1)+&
              & sflag2(:,l2:u2-1)*weights2(:,l2+1:u2))
ENDIF

!
!5. Interpolate
!--------------
!
!5.1 Without weights
!-------------------
!

IF (hreg) THEN

   l_510: DO l=1,nosize
   k_510: DO k=l3,u3-1
      WHERE (fdest3(:,:,k).AND.wsum2.GT.0.0)
         xout(:,:,k,l) = &
           &(sflag2(:,l2:u2-1)*(xin(:,l2:u2-1,k,l)+xin(:,l2:u2-1,k+1,l))+&
           & sflag2(:,l2+1:u2)*(xin(:,l2+1:u2,k,l)+xin(:,l2+1:u2,k+1,l)))/wsum2
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_510
   ENDDO l_510

!
!5.2 With weights
!----------------
!

ELSE

   l_520: DO l=1,nosize
   k_520: DO k=l3,u3-1
      WHERE (fdest3(:,:,k).AND.wsum2.GT.0.0)
         xout(:,:,k,l) = (sflag2(:,l2:u2-1)*weights2(:,l2+1:u2)*&
                        &(xin(:,l2:u2-1,k,l)+xin(:,l2:u2-1,k+1,l))+&
                        &(sflag2(:,l2+1:u2)*weights2(:,l2:u2-1)*&
                        &(xin(:,l2+1:u2,k,l)+xin(:,l2+1:u2,k+1,l))))/wsum2
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_520
   ENDDO l_520

ENDIF

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag2,fdest3,wsum2)
IF (.NOT.hreg) DEALLOCATE (weights2)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Warr_at_V

!========================================================================

SUBROUTINE Warr_at_VW(xin,xout,intdest,lbounds,ubounds,nosize,ivarid,info,&
                    & outflag,hregular)
!************************************************************************
!
! *Warr_at_VW* Interpolate an array from W- to VW-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: intdest, ivarid, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(IN), DIMENSION(lbounds(1):ubounds(1),lbounds(2):ubounds(2),&
                          & lbounds(3):ubounds(3),nosize) :: xin
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1),lbounds(2)+1:ubounds(2),&
                           & lbounds(3):ubounds(3),nosize) :: xout

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at W-nodes
!*xout*      REAL    Output array at VW-nodes
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*lbounds*   INTEGER Lower bounds of input array (first three dimensions) 
!*ubounds*   INTEGER Upper bounds of input array (first three dimensions)
!*nosize*    INTEGER Fourth dimension of input array
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'Warr_at_VW'
CALL log_timer_in(npcc,ivarid,info=info)

!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dYregY
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!3. Data flags
!-------------
!
!---allocate
ALLOCATE (sflag2(l1:u1,l2:u2),STAT=errstat)
CALL error_alloc('sflag2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ALLOCATE (fdest3(l1:u1,l2+1:u2,l3:u3),STAT=errstat)
CALL error_alloc('fdest3',3,(/u1-l1+1,u2-l2,u3-l3+1/),kndlog)

!---source
sflag2 = MERGE(1.0,0.0,nodeatc(l1:u1,l2:u2).GT.0)

!---destination
SELECT CASE (intdest)
   CASE (1); fdest3 = nodeatvw(l1:u1,l2+1:u2,l3:u3).GT.0
   CASE (2); fdest3 = nodeatvw(l1:u1,l2+1:u2,l3:u3).EQ.2
   CASE (3); fdest3 = nodeatvw(l1:u1,l2+1:u2,l3:u3).GT.1
   CASE (4); fdest3 = nodeatvw(l1:u1,l2+1:u2,l3:u3).EQ.1.OR.&
                    & nodeatvw(l1:u1,l2+1:u2,l3:u3).EQ.2
END SELECT

!
!4. Weight factors
!-----------------
!
!---allocate
IF (.NOT.hreg) THEN
   ALLOCATE (weights2(l1:u1,l2:u2),STAT=errstat)
   CALL error_alloc('weights2',2,(/u1-l1+1,u2-l2+1/),kndrtype)
ENDIF
ALLOCATE (wsum2(l1:u1,l2+1:u2),STAT=errstat)
CALL error_alloc('wsum2',2,(/u1-l1+1,u2-l2/),kndrtype)

!---define
IF (hreg) THEN
   wsum2 = sflag2(:,l2+1:u2) + sflag2(:,l2:u2-1)
ELSE
   weights2 = delyatc(l1:u1,l2:u2)
   wsum2 = sflag2(:,l2+1:u2)*weights2(:,l2:u2-1) + &
         & sflag2(:,l2:u2-1)*weights2(:,l2+1:u2)
ENDIF

k_410: DO k=l3,u3
   WHERE (wsum2.LE.0.0) fdest3(:,:,k) = .FALSE.
ENDDO k_410

!
!5. Interpolate
!--------------
!
!5.1 Without weights
!-------------------
!

IF (hreg) THEN

   l_510: DO l=1,nosize
   k_510: DO k=l3,u3
      WHERE (fdest3(:,:,k))
         xout(:,:,k,l) = (sflag2(:,l2:u2-1)*xin(:,l2:u2-1,k,l)+&
                          sflag2(:,l2+1:u2)*xin(:,l2+1:u2,k,l))/wsum2
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_510
   ENDDO l_510

!
!5.2 With weights
!----------------
!

ELSE

   l_520: DO l=1,nosize
   k_520: DO k=l3,u3
      WHERE (fdest3(:,:,k))
         xout(:,:,k,l) = &
            &(sflag2(:,l2:u2-1)*weights2(:,l2+1:u2)*xin(:,l2:u2-1,k,l)+&
            & sflag2(:,l2+1:u2)*weights2(:,l2:u2-1)*xin(:,l2+1:u2,k,l))/wsum2
      ELSEWHERE
         xout(:,:,k,l) = flagval
      END WHERE
   ENDDO k_520
   ENDDO l_520

ENDIF

!
!6. Deallocate
!-------------
!        

DEALLOCATE (sflag2,fdest3,wsum2)
IF (.NOT.hreg) DEALLOCATE (weights2)

CALL log_timer_out(npcc,itm_arrint,info)


RETURN

END SUBROUTINE Warr_at_VW

!========================================================================

FUNCTION Cvar_at_U(xin,i,j,k,intsrce,intdest,outflag,hregular)
!************************************************************************
!
! *Cvar_at_U* Interpolate a C-node array at U-node grid point (i+1,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: Cvar_at_U

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER Lower X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => wet points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the
!                     horizontal
!*Cvar_at_U* REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag, weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXregX
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE (intsrce)
   CASE (0); sflag = 1.0
   CASE (1); sflag = MERGE(1.0,0.0,nodeatc(i:i+1,j).GT.0)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest = .TRUE.
   CASE (1); fdest = nodeatu(i+1,j,k).GT.0
   CASE (2); fdest = nodeatu(i+1,j,k).EQ.2
   CASE (3); fdest = nodeatu(i+1,j,k).GT.1
   CASE (4); fdest = nodeatu(i+1,j,k).EQ.1.OR.nodeatu(i+1,j,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

IF (hreg) THEN
   wsum = SUM(sflag)
ELSE
   weights = delxatc(i:i+1,j)
   wsum = sflag(2)*weights(1)+sflag(1)*weights(2)
ENDIF

IF (wsum.LE.0.0) fdest = .FALSE.

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (hreg) THEN
   IF (fdest) THEN
      Cvar_at_U = (sflag(1)*xin(1)+sflag(2)*xin(2))/wsum
   ELSE
      Cvar_at_U = flagval
   ENDIF

!
!4.2 With weights
!----------------
!

ELSE
   IF (fdest) THEN
      Cvar_at_U = (sflag(1)*weights(2)*xin(1)+sflag(2)*weights(1)*xin(2))/wsum
   ELSE
      Cvar_at_U = flagval
   ENDIF
ENDIF


RETURN

END FUNCTION Cvar_at_U

!========================================================================

FUNCTION Cvar_at_UV(xin,i,j,k,intsrce,intdest,outflag,hregular)
!************************************************************************
!
! *Cvar_at_UV* Interpolate a C-node array at UV-node grid point (i+1,j+1)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2,2) :: xin
REAL :: Cvar_at_UV

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER Lower X-index of interpolating array
!*j*         INTEGER Lower Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => wet points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => wet points only
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the
!                    horizontal
!*Cvar_at_UV*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2,2) :: sflag, weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXYreg
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE (intsrce)
   CASE (0); sflag = 1.0
   CASE (1); sflag = MERGE(1.0,0.0,nodeatc(i:i+1,j:j+1).GT.0)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest = .TRUE.
   CASE (1); fdest = nodeatuv(i+1,j+1,k).GT.0
END SELECT

!
!3. Weight factors
!-----------------
!

IF (hreg) THEN
   wsum = SUM(sflag)
ELSE
   IF (dXregX.AND.dXregY) THEN
      weights = delyatc(i:i+1,j:j+1)
   ELSEIF (dYregX.AND.dYregY) THEN
      weights = delxatc(i:i+1,j:j+1)
   ELSE
      weights = delxatc(i:i+1,j:j+1)*delyatc(i:i+1,j:j+1)
   ENDIF
   wsum = sflag(2,2)*weights(1,1)+sflag(1,2)*weights(2,1)+&
        & sflag(2,1)*weights(1,2)+sflag(1,1)*weights(2,2) 
ENDIF

IF (wsum.LE.0.0) fdest = .FALSE.

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (hreg) THEN
   IF (fdest) THEN
      Cvar_at_UV = (sflag(1,1)*xin(1,1)+sflag(2,1)*xin(2,1)+&
                  & sflag(2,1)*xin(2,1)+sflag(2,2)*xin(2,2))/wsum
   ELSE
      Cvar_at_UV = flagval
   ENDIF

!
!4.2 With weights
!----------------
!

ELSE
   IF (fdest) THEN
      Cvar_at_UV = &
          &(sflag(1,1)*weights(2,2)*xin(1,1)+sflag(2,1)*weights(1,2)*xin(2,1)+&
          & sflag(1,2)*weights(2,1)*xin(1,2)+sflag(2,2)*weights(1,1)*xin(2,2))&
          & /wsum
   ELSE
      Cvar_at_UV = flagval
   ENDIF
ENDIF


RETURN

END FUNCTION Cvar_at_UV

!========================================================================

FUNCTION Cvar_at_UW(xin,i,j,k,intdest,outflag,hregular,vregular)
!************************************************************************
!
! *Cvar_at_UW* Interpolate a C-node array at UW-node grid point (i+1,j,k+1)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular, vregular
INTEGER, INTENT(IN) :: i, intdest, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2,2) :: xin
REAL :: Cvar_at_UW

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER Lower X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Lower Z-index of interpolating array
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the
!                     horizontal
!*vregular*  LOGICAL Flag to select uniform or area averaging in the vertical
!*Cvar_at_UW*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: hvreg
INTEGER :: kk
REAL, DIMENSION(2) :: sflag
REAL, DIMENSION(2,2) :: weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXregX
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(vregular)) THEN
   vreg = iopt_arrint_vreg.EQ.0
ELSE
   vreg = vregular
ENDIF
hvreg = hreg.AND.vreg

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
sflag = MERGE(1.0,0.0,nodeatc(i:i+1,j).GT.0)

!---destination
SELECT CASE (intdest)
   CASE (1); fdest = nodeatuw(i+1,j,k).GT.0
   CASE (2); fdest = nodeatuw(i+1,j,k).EQ.2
   CASE (3); fdest = nodeatuw(i+1,j,k).GT.1
   CASE (4); fdest = nodeatuw(i+1,j,k).EQ.1.OR.nodeatuw(i+1,j,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

IF (hvreg) THEN
   wsum = 2.0*SUM(sflag)
ELSE
   k_310: DO kk=1,2
      weights(:,kk) = 1.0
      IF (.NOT.hreg) THEN
         weights(:,kk) = weights(:,kk)*delxatc(i:i+1,j)
      ENDIF
      IF (.NOT.vreg) THEN
         weights(:,kk) = weights(:,kk)*delzatc(i:i+1,j,k+kk-1)
      ENDIF
   ENDDO k_310
   IF (node2du(i+1,j).GE.2) THEN
      IF (nodeatc(i,j).GT.0) THEN
         weights(2,:) = weights(1,:)
      ELSE
         weights(1,:) = weights(2,:)
      ENDIF
   ENDIF
   wsum =  sflag(2)*(weights(1,1)+weights(1,2))+ &
         & sflag(1)*(weights(2,1)+weights(2,2))
ENDIF

IF (wsum.LE.0.0) fdest = .FALSE.

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (hvreg) THEN
      IF (fdest) THEN
         Cvar_at_UW = sflag(1)*(xin(1,1)+xin(2,1))+sflag(2)*(xin(1,2)+xin(2,2))
      ELSE
         Cvar_at_UW = flagval
      ENDIF

!
!4.2 With weights
!----------------
!

ELSE
   IF (fdest) THEN
      Cvar_at_UW = &
          &(sflag(1)*weights(2,2)*xin(1,1)+sflag(2)*weights(1,2)*xin(2,1) +&
          & sflag(1)*weights(2,1)*xin(1,2)+sflag(2)*weights(1,1)*xin(2,2))/wsum
   ELSE
      Cvar_at_UW = flagval
   ENDIF
ENDIF


RETURN

END FUNCTION Cvar_at_UW

!========================================================================

FUNCTION Cvar_at_V(xin,i,j,k,intsrce,intdest,outflag,hregular)
!************************************************************************
!
! *Cvar_at_V* Interpolate a C-node array at V-node grid point (i,j+1)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: Cvar_at_V

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Lower Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                     = 1 => wet points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the
!                     horizontal
!*Cvar_at_V* REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag, weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dYregY
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!

!---source
SELECT CASE (intsrce)
   CASE (0); sflag = 1.0
   CASE (1); sflag = MERGE(1.0,0.0,nodeatc(i,j:j+1).GT.0)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest = .TRUE.
   CASE (1); fdest = nodeatv(i,j+1,k).GT.0
   CASE (2); fdest = nodeatv(i,j+1,k).EQ.2
   CASE (3); fdest = nodeatv(i,j+1,k).GT.1
   CASE (4); fdest = nodeatv(i,j+1,k).EQ.1.OR.nodeatv(i,j+1,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

IF (hreg) THEN
   wsum = SUM(sflag)
ELSE
   weights = delyatc(i,j:j+1)
   wsum = sflag(2)*weights(1)+sflag(1)*weights(2)
ENDIF

IF (wsum.LE.0.0) fdest = .FALSE.

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (hreg) THEN
   IF (fdest) THEN
      Cvar_at_V = (sflag(1)*xin(1)+sflag(2)*xin(2))/wsum
   ELSE
      Cvar_at_V = flagval
   ENDIF

!
!4.2 With weights
!----------------
!

ELSE
   IF (fdest) THEN
      Cvar_at_V = (sflag(1)*weights(2)*xin(1)+sflag(2)*weights(1)*xin(2))/wsum
   ELSE
      Cvar_at_V = flagval
   ENDIF
ENDIF


RETURN

END FUNCTION Cvar_at_V

!========================================================================

FUNCTION Cvar_at_VW(xin,i,j,k,intdest,outflag,hregular,vregular)
!************************************************************************
!
! *Cvar_at_VW* Interpolate a C-node array at VW-node grid point (i,j+1,k+1)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular, vregular
INTEGER, INTENT(IN) :: i, intdest, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2,2) :: xin
REAL :: Cvar_at_VW

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Lower Y-index of interpolating array
!*k*         INTEGER Lower Z-index of interpolating array
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the
!                    horizontal
!*vregular*  LOGICAL Flag to select uniform or area averaging in the vertical
!*Cvar_at_VW*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: hvreg
INTEGER :: kk
REAL, DIMENSION(2) :: sflag
REAL, DIMENSION(2,2) :: weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXregX
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(vregular)) THEN
   vreg = iopt_arrint_vreg.EQ.0
ELSE
   vreg = vregular
ENDIF
hvreg = hreg.AND.vreg

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flag
!------------
!
!---source
sflag = MERGE(1.0,0.0,nodeatc(i,j:j+1).GT.0)

!---destination
SELECT CASE (intdest)
   CASE (1); fdest = nodeatvw(i,j+1,k).GT.0
   CASE (2); fdest = nodeatvw(i,j+1,k).EQ.2
   CASE (3); fdest = nodeatvw(i,j+1,k).GT.1
   CASE (4); fdest = nodeatvw(i,j+1,k).EQ.1.OR.nodeatvw(i,j+1,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

IF (hvreg) THEN
   wsum = 2.0*SUM(sflag)
ELSE
   k_310: DO kk=1,2
      weights(:,kk) = 1.0
      IF (.NOT.hreg) THEN
         weights(:,kk) = weights(:,kk)*delyatc(i,j:j+1)
      ENDIF
      IF (.NOT.vreg) THEN
         weights(:,kk) = weights(:,kk)*delzatc(i,j:j+1,k+kk-1)
      ENDIF
   ENDDO k_310
   IF (node2dv(i,j+1).GE.2) THEN
      IF (nodeatc(i,j).GT.0) THEN
         weights(2,:) = weights(1,:)
      ELSE
         weights(1,:) = weights(2,:)
      ENDIF
   ENDIF
   wsum =  sflag(2)*(weights(1,1)+weights(1,2))+ &
         & sflag(1)*(weights(2,1)+weights(2,2))
ENDIF

IF (wsum.LE.0.0) fdest = .FALSE. 

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (hvreg) THEN
   IF (fdest) THEN
      Cvar_at_VW = sflag(1)*(xin(1,1)+xin(1,2))+sflag(2)*(xin(2,1)+xin(2,2))
   ELSE
      Cvar_at_VW = flagval
   ENDIF

!
!4.2 With weights
!----------------
!

ELSE
   IF (fdest) THEN
      Cvar_at_VW = &
          &(sflag(1)*weights(2,2)*xin(1,1)+sflag(2)*weights(1,2)*xin(2,1) +&
          & sflag(1)*weights(2,1)*xin(1,2)+sflag(2)*weights(1,1)*xin(2,2))/wsum
   ELSE
      Cvar_at_VW = flagval
   ENDIF
ENDIF


RETURN

END FUNCTION Cvar_at_VW

!========================================================================

FUNCTION Cvar_at_W(xin,i,j,k,outflag,vregular)
!************************************************************************
!
! *Cvar_at_W* Interpolate a C-node array at W-node grid point (i,j,k+1)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: vregular
INTEGER, INTENT(IN) :: i, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: Cvar_at_W

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Lower Z-index of interpolating array
!*outflag*   REAL    Output flag for dry points
!*vregular*  LOGICAL Flag to select uniform or area averaging in the vertical
!*Cvar_at_W* REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(vregular)) THEN
   vreg = iopt_arrint_vreg.EQ.0.OR.dZregZ
ELSE
   vreg = vregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Weight factors
!-----------------
!

IF (.NOT.vreg) THEN
   weights = delzatc(i,j,k:k+1)
   wsum = SUM(weights)
ENDIF

!
!3. Interpolate
!--------------
!
!3.1 Without weights
!-------------------
!

IF (vreg) THEN
   Cvar_at_W = MERGE(0.5*SUM(xin),flagval,nodeatc(i,j).GT.0)

!
!4.2 With weights
!----------------
!

ELSE
   IF (nodeatc(i,j).GT.0) THEN
      Cvar_at_W = (weights(2)*xin(1)+weights(1)*xin(2))/wsum
   ELSE
      Cvar_at_W = flagval
   ENDIF
ENDIF


RETURN

END FUNCTION Cvar_at_W

!========================================================================

FUNCTION Uvar_at_C(xin,i,j,k,intsrce,intdest,outflag)
!************************************************************************
!
! *Uvar_at_C* Interpolate a U-node array at C-node grid point (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: Uvar_at_C

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => wet points only
!*outflag*   REAL    Output flag for dry points
!*Uvar_at_C* REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE (intsrce)
   CASE (0); sflag = 1.0
   CASE (1); sflag = MERGE(1.0,0.0,nodeatu(i:i+1,j,k).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatu(i:i+1,j,k).EQ.2)
   CASE (3); sflag = MERGE(1.0,0.0,nodeatu(i:i+1,j,k).GT.1)
   CASE (4); sflag = MERGE(1.0,0.0,nodeatu(i:i+1,j,k).EQ.1.OR.&
                                 & nodeatu(i:i+1,j,k).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest = .TRUE.
   CASE (1); fdest = nodeatc(i,j).GT.0
END SELECT

!
!3. Weight factors
!-----------------
!

wsum = SUM(sflag)

!
!4. Interpolate
!--------------
!

IF (fdest.AND.wsum.GT.0.0) THEN
   Uvar_at_C = SUM(sflag*xin)/wsum
ELSE
   Uvar_at_C = flagval
ENDIF


RETURN

END FUNCTION Uvar_at_C

!========================================================================

FUNCTION Uvar_at_UV(xin,i,j,k,intsrce,intdest,outflag,hregular)
!************************************************************************
!
! *Uvar_at_UV* Interpolate a U-node array at UV-node grid point (i,j+1)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: Uvar_at_UV

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Lower Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => interior points and open boundaries only
!                    = 2 => interior points and UV-node open boundaries only
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the
!                    horizontal
!*Uvar_at_UV*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL :: wsum
REAL, DIMENSION(2) :: sflag, weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dYregY
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE (intsrce)
   CASE (0); sflag = 1.0
   CASE (1); sflag = MERGE(1.0,0.0,nodeatu(i,j:j+1,k).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatu(i,j:j+1,k).EQ.2)
   CASE (3); sflag = MERGE(1.0,0.0,nodeatu(i,j:j+1,k).GT.1)
   CASE (4); sflag = MERGE(1.0,0.0,nodeatu(i,j:j+1,k).EQ.1.OR.&
                                 & nodeatu(i,j:j+1,k).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest = .TRUE.
   CASE (1); fdest = nodeatuv(i,j+1,k).GT.0
   CASE (2); fdest = nodeatuv(i,j+1,k).EQ.1.OR.nodeatuv(i,j+1,k).EQ.3
END SELECT

!
!3. Weight factors
!-----------------
!

IF (hreg) THEN
   wsum = SUM(sflag)
ELSE
   weights = delyatu(i,j:j+1)
   wsum = sflag(2)*weights(1) + sflag(1)*weights(2)
ENDIF

IF (wsum.LE.0.0) fdest = .FALSE.

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (hreg) THEN

   IF (fdest) THEN
      Uvar_at_UV = SUM(sflag*xin)/wsum
   ELSE
      Uvar_at_UV = flagval
   ENDIF

!
!4.2 With weights
!----------------
!

ELSE

   IF (fdest) THEN
      Uvar_at_UV = (sflag(1)*weights(2)*xin(1)+&
                  & sflag(2)*weights(1)*xin(2))/wsum
   ELSE
      Uvar_at_UV = flagval
   ENDIF

ENDIF


RETURN

END FUNCTION Uvar_at_UV

!========================================================================

FUNCTION Uvar_at_UW(xin,i,j,k,intsrce,intdest,outflag,vregular)
!************************************************************************
!
! *Uvar_at_UW* Interpolate a U-node array at UW-node grid point (i,j,k+1)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: vregular
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: Uvar_at_UW

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Lower Z-index of interpolating array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*vregular*  LOGICAL Flag to select uniform or area averaging in the vertical
!*Uvar_at_UW*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag, weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(vregular)) THEN
   vreg = iopt_arrint_vreg.EQ.0.OR.dZregZ
ELSE
   vreg = vregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE (intsrce)
   CASE (1); sflag = MERGE(1.0,0.0,nodeatu(i,j,k:k+1).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatu(i,j,k:k+1).EQ.2)
   CASE (3); sflag = MERGE(1.0,0.0,nodeatu(i,j,k:k+1).GT.1)
   CASE (4); sflag = MERGE(1.0,0.0,nodeatu(i,j,k:k+1).EQ.1.OR.&
                                 & nodeatu(i,j,k:k+1).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (1); fdest = nodeatuw(i,j,k).GT.0
   CASE (2); fdest = nodeatuw(i,j,k).EQ.2
   CASE (3); fdest = nodeatuw(i,j,k).GT.1
   CASE (4); fdest = nodeatuw(i,j,k).EQ.1.OR.nodeatuw(i,j,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

IF (vreg) THEN
   wsum = SUM(sflag)
ELSE
   weights = delzatu(i,j,k:k+1)
   wsum = sflag(2)*weights(1) + sflag(1)*weights(2)
ENDIF

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (vreg) THEN

   Uvar_at_UW = MERGE(SUM(sflag*xin)/wsum,flagval,fdest)

!
!4.2 With weights
!----------------
!

ELSE

   IF (fdest) THEN
      Uvar_at_UW = (sflag(1)*weights(2)*xin(1)+sflag(2)*weights(1)*xin(2))/wsum
   ELSE
      Uvar_at_UW = flagval
   ENDIF

ENDIF


RETURN

END FUNCTION Uvar_at_UW

!========================================================================

FUNCTION Uvar_at_V(xin,i,j,k,intsrce,intdest,outflag,hregular)
!************************************************************************
!
! *Uvar_at_V* Interpolate a U-node array at V-node grid point (i,j+1)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2,2) :: xin
REAL :: Uvar_at_V

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Lower Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the
!                    horizontal
!*Uvar_at_V* REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2,2) :: sflag, weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXYreg
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE (intsrce)
   CASE (0); sflag = 1.0
   CASE (1); sflag = MERGE(1.0,0.0,nodeatu(i:i+1,j:j+1,k).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatu(i:i+1,j:j+1,k).EQ.2)
   CASE (3); sflag = MERGE(1.0,0.0,nodeatu(i:i+1,j:j+1,k).GT.1)
   CASE (4); sflag = MERGE(1.0,0.0,nodeatu(i:i+1,j:j+1,k).EQ.1.OR.&
                                 & nodeatu(i:i+1,j:j+1,k).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest = .TRUE.
   CASE (1); fdest = nodeatv(i,j+1,k).GT.0
   CASE (2); fdest = nodeatv(i,j+1,k).EQ.2
   CASE (3); fdest = nodeatv(i,j+1,k).GT.1
   CASE (4); fdest = nodeatv(i,j+1,k).EQ.1.OR.nodeatv(i,j+1,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

IF (hreg) THEN
   wsum = SUM(sflag)
ELSE
   IF (dXregX.AND.dXregY) THEN 
      weights = delyatv(i:i+1,j:j+1)
   ELSEIF (dYregX.AND.dYregY) THEN
      weights = delxatu(i:i+1,j:j+1)
   ELSE
      weights = delxatu(i:i+1,j:j+1)*delyatv(i:i+1,j:j+1)
   ENDIF
   wsum = sflag(2,2)*weights(1,1)+sflag(1,2)*weights(2,1)+&
        & sflag(2,1)*weights(1,2)+sflag(1,1)*weights(2,2)
ENDIF
fdest = fdest.AND.wsum.GT.0.0

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (hreg) THEN

   IF (fdest) THEN
      Uvar_at_V = SUM(sflag*xin)/wsum
   ELSE
      Uvar_at_V = flagval
   ENDIF

!
!4.2 With weights
!----------------
!

ELSE

   IF (fdest) THEN
      Uvar_at_V = (sflag(1,1)*weights(2,2)*xin(1,1) + &
                 & sflag(2,1)*weights(1,2)*xin(2,1) + &
                 & sflag(1,2)*weights(2,1)*xin(1,2) + &
                 & sflag(2,2)*weights(1,1)*xin(2,2))/wsum
   ELSE
      Uvar_at_V = flagval
   ENDIF

ENDIF


RETURN

END FUNCTION Uvar_at_V

!========================================================================

FUNCTION Uvar_at_W(xin,i,j,k,intsrce,outflag,vregular)
!************************************************************************
!
! *Uvar_at_W* Interpolate a U-node array at W-node grid point (i,j,k+1)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: vregular
INTEGER, INTENT(IN) :: i, intsrce, j ,k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2,2) :: xin
REAL :: Uvar_at_W

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER Lower X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Lower Z-index of interpolating array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*vregular*  LOGICAL Flag to select uniform or area averaging in the vertical
!*Uvar_at_W* REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: kk, k1
REAL, DIMENSION(2,2) :: sflag
REAL, DIMENSION(2,2) :: weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(vregular)) THEN
   vreg = iopt_arrint_vreg.EQ.0
ELSE
   vreg = vregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE (intsrce)
   CASE (1); sflag = MERGE(1.0,0.0,nodeatu(i:i+1,j,k:k+1).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatu(i:i+1,j,k:k+1).EQ.2)
   CASE (3); sflag = MERGE(1.0,0.0,nodeatu(i:i+1,j,k:k+1).GT.1)
   CASE (4); sflag = MERGE(1.0,0.0,nodeatu(i:i+1,j,k:k+1).EQ.1.OR.&
                                 & nodeatu(i:i+1,j,k:k+1).EQ.2)
END SELECT

!---destination
fdest = nodeatc(i,j).GT.0

!
!3. Weight factors
!-----------------
!

IF (vreg) THEN
   wsum = SUM(sflag)
ELSE
   kk_310: DO kk=1,2
      k1 = k+kk-1
      weights(:,kk) = delzatu(i:i+1,j,k1)
      IF (nodeatu(i+1,j,k1).GE.2.AND.nodeatc(i,j).GT.0) THEN
         weights(2,kk) = delzatc(i,j,k1)
      ENDIF
      IF (nodeatu(i,j,k1).GE.2.AND.nodeatc(i,j).GT.0) THEN
         weights(1,kk) = delzatc(i,j,k1)
      ENDIF
   ENDDO kk_310
   wsum = sflag(1,1)*weights(2,2)+sflag(2,1)*weights(1,2)+&
        & sflag(1,2)*weights(2,1)+sflag(2,2)*weights(1,1)
ENDIF

IF (wsum.LE.0.0) fdest = .FALSE. 

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (vreg) THEN

   IF (fdest) THEN
      Uvar_at_W = (sflag(1,1)*xin(1,1)+sflag(1,2)*xin(1,2)+&
                &  sflag(2,1)*xin(2,1)+sflag(2,2)*xin(2,2))/wsum
   ELSE
      Uvar_at_W = flagval
   ENDIF

!
!4.2 With weights
!----------------
!

ELSE

   IF (fdest) THEN
      Uvar_at_W = (sflag(1,1)*weights(2,2)*xin(1,1)+&
                 & sflag(2,1)*weights(1,2)*xin(2,1)+&
                 & sflag(1,2)*weights(2,1)*xin(1,2)+&
                 & sflag(2,2)*weights(1,1)*xin(2,2))/wsum
   ELSE
      Uvar_at_W = flagval
   ENDIF

ENDIF


RETURN

END FUNCTION Uvar_at_W

!========================================================================

FUNCTION UVvar_at_C(xin,i,j,k,intsrce,intdest,outflag)
!************************************************************************
!
! *UVvar_at_C* Interpolate a UV-node array at C-node grid point (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2,2) :: xin
REAL :: UVvar_at_C

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER Lower X-index of interpolating array
!*j*         INTEGER Lower Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => wet points only
!*intdest*   INTEGER Selects valid points at destination node
!                     = 0 => all points
!                     = 1 => wet points only
!*outflag*   REAL    Output flag for dry points
!*UVvar_at_C*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2,2) :: sflag


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE (intsrce)
   CASE (0); sflag = 1.0
   CASE (1); sflag = MERGE(1.0,0.0,nodeatuv(i:i+1,j:j+1,k).GT.0)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest = .TRUE.
   CASE (1); fdest = nodeatc(i,j).GT.0
END SELECT

!
!3. Weight factors
!-----------------
!

wsum = SUM(sflag)
IF (wsum.LE.0.0) fdest = .FALSE.

!
!4. Interpolate
!--------------
!

IF (fdest) THEN
   UVvar_at_C = SUM(sflag*xin)/wsum
ELSE
   UVvar_at_C = flagval
ENDIF


RETURN

END FUNCTION UVvar_at_C

!========================================================================

FUNCTION UVvar_at_U(xin,i,j,k,intsrce,intdest,outflag)
!************************************************************************
!
! *UVvar_at_U* Interpolate a UV-node array at U-node grid point (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: UVvar_at_U

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Lower Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => interior points and open boundaries only
!                    = 2 => interior points and UV-node open boundaries only
!*intdest*   INTEGER Selects valid points at source and destination node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*UVvar_at_U*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE (intsrce)
   CASE (0); sflag = 1.0
   CASE (1); sflag = MERGE(1.0,0.0,nodeatuv(i,j:j+1,k).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatuv(i,j:j+1,k).EQ.1.OR.&
                                 & nodeatuv(i,j:j+1,k).EQ.3)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest = .TRUE.
   CASE (1); fdest = nodeatu(i,j,k).GT.0
   CASE (2); fdest = nodeatu(i,j,k).EQ.2
   CASE (3); fdest = nodeatu(i,j,k).GT.1
   CASE (4); fdest = nodeatu(i,j,k).EQ.1.OR.nodeatu(i,j,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

wsum = SUM(sflag)

!
!4. Interpolate
!--------------
!

IF (fdest.AND.wsum.GT.0.0) THEN
   UVvar_at_U = SUM(sflag*xin)/wsum
ELSE
   UVvar_at_U = flagval
ENDIF


RETURN

END FUNCTION UVvar_at_U

!========================================================================

FUNCTION UVvar_at_V(xin,i,j,k,intsrce,intdest,outflag)
!************************************************************************
!
! *UVvar_at_V* Interpolate a UV-node array at V-node grid point (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: UVvar_at_V

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => interior points and open boundaries only
!                    = 2 => interior points and UV-node open boundaries only
!*intdest*   INTEGER Selects valid points at source and destination node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*UVvar_at_V*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE (intsrce)
   CASE (0); sflag = 1.0
   CASE (1); sflag = MERGE(1.0,0.0,nodeatuv(i:i+1,j,k).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatuv(i:i+1,j,k).EQ.1.OR.&
                                 & nodeatuv(i:i+1,j,k).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest = .TRUE.
   CASE (1); fdest = nodeatv(i,j,k).GT.0
   CASE (2); fdest = nodeatv(i,j,k).EQ.2
   CASE (3); fdest = nodeatv(i,j,k).GT.1
   CASE (4); fdest = nodeatv(i,j,k).EQ.1.OR.nodeatv(i,j,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

wsum = SUM(sflag)

!
!4. Interpolate
!--------------
!

IF (fdest.AND.wsum.GT.0.0) THEN
   UVvar_at_V = SUM(sflag*xin)/wsum
ELSE
   UVvar_at_V = flagval
ENDIF


RETURN

END FUNCTION UVvar_at_V

!========================================================================

FUNCTION UWvar_at_U(xin,i,j,k,intsrce,intdest,outflag)
!************************************************************************
!
! *UWvar_at_U* Interpolate a UW-node array at U-node grid point (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: UWvar_at_U

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*UWvar_at_U*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag


!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE( intsrce)
   CASE (1); sflag = MERGE(1.0,0.0,nodeatuw(i,j,k:k+1).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatuw(i,j,k:k+1).EQ.2)
   CASE (3); sflag = MERGE(1.0,0.0,nodeatuw(i,j,k:k+1).GT.1)
   CASE (4); sflag = MERGE(1.0,0.0,nodeatuw(i,j,k:k+1).EQ.1.OR.&
                                 & nodeatuw(i,j,k:k+1).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (1); fdest = nodeatu(i,j,k).GT.0
   CASE (2); fdest = nodeatu(i,j,k).EQ.2
   CASE (3); fdest = nodeatu(i,j,k).GT.1
   CASE (4); fdest = nodeatu(i,j,k).EQ.1.OR.nodeatu(i,j,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

wsum = SUM(sflag)

!
!4. Interpolate
!--------------
!

IF (fdest.AND.wsum.GT.0) THEN
   UWvar_at_U = SUM(sflag*xin)/wsum
ELSE
   UWvar_at_U = flagval
ENDIF


RETURN

END FUNCTION UWvar_at_U

!========================================================================

FUNCTION UWvar_at_W(xin,i,j,k,intsrce,outflag)
!************************************************************************
!
! *UWvar_at_W* Interpolate a UW-node array at W-node grid point (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: UWvar_at_W

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER Lower X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*UWvar_at_W*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flag
!------------
!
!---source
SELECT CASE (intsrce)
   CASE (1); sflag = MERGE(1.0,0.0,nodeatuw(i:i+1,j,k).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatuw(i:i+1,j,k).EQ.2)
   CASE (3); sflag = MERGE(1.0,0.0,nodeatuw(i:i+1,j,k).GT.1)
   CASE (4); sflag = MERGE(1.0,0.0,nodeatuw(i:i+1,j,k).EQ.1.OR.&
                                 & nodeatuw(i:i+1,j,k).EQ.2)
END SELECT

!---destination
fdest = nodeatc(i,j).GT.0

!
!3. Weight factors
!-----------------
!

wsum = SUM(sflag)

!
!4. Interpolate
!--------------
!

IF (fdest.AND.wsum.GT.0.0) THEN
   UWvar_at_W = SUM(sflag*xin)/wsum
ELSE
   UWvar_at_W = flagval
ENDIF


RETURN

END FUNCTION UWvar_at_W

!========================================================================

FUNCTION Vvar_at_C(xin,i,j,k,intsrce,intdest,outflag)
!************************************************************************
!
! *Vvar_at_C* Interpolate a V-node array at C-node grid point (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: Vvar_at_C

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                     = 0 => all points
!                     = 1 => wet points only
!*outflag*   REAL    Output flag for dry points
!*Vvar_at_C* REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flag and weight factors
!-------------------------------
!
!---source
SELECT CASE (intsrce)
   CASE (0); sflag = 1.0
   CASE (1); sflag = MERGE(1.0,0.0,nodeatv(i,j:j+1,k).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatv(i,j:j+1,k).EQ.2)
   CASE (3); sflag = MERGE(1.0,0.0,nodeatv(i,j:j+1,k).GT.1)
   CASE (4); sflag = MERGE(1.0,0.0,nodeatv(i,j:j+1,k).EQ.1.OR.&
                                 & nodeatv(i,j:j+1,k).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest = .TRUE.
   CASE (1); fdest = nodeatc(i,j).GT.0
END SELECT

!
!3. Weight factors
!-----------------
!

wsum = SUM(sflag)

!
!4. Interpolate
!--------------
!

IF (nodeatc(i,j).GT.0.AND.wsum.GT.0.0) THEN
   Vvar_at_C = SUM(sflag*xin)/wsum
ELSE
   Vvar_at_C = flagval
ENDIF


RETURN

END FUNCTION Vvar_at_C

!========================================================================

FUNCTION Vvar_at_U(xin,i,j,k,intsrce,intdest,outflag,hregular)
!************************************************************************
!
! *Vvar_at_U* Interpolate a V-node array at U-node grid point (i+1,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2,2) :: xin
REAL :: Vvar_at_U

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER Lower X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the
!                    horizontal
!*Vvar_at_U* REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2,2) :: sflag, weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXYreg
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!

!---source
SELECT CASE (intsrce)
   CASE (0); sflag = 1.0
   CASE (1); sflag = MERGE(1.0,0.0,nodeatv(i:i+1,j:j+1,k).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatv(i:i+1,j:j+1,k).EQ.2)
   CASE (3); sflag = MERGE(1.0,0.0,nodeatv(i:i+1,j:j+1,k).GT.1)
   CASE (4); sflag = MERGE(1.0,0.0,nodeatv(i:i+1,j:j+1,k).EQ.1.OR.&
                                 & nodeatv(i:i+1,j:j+1,k).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest = .TRUE.
   CASE (1); fdest = nodeatu(i+1,j,k).GT.0
   CASE (2); fdest = nodeatu(i+1,j,k).EQ.2
   CASE (3); fdest = nodeatu(i+1,j,k).GT.1
   CASE (4); fdest = nodeatu(i+1,j,k).EQ.1.OR.nodeatu(i+1,j,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

IF (hreg) THEN
   wsum = SUM(sflag)
ELSE
   IF (dXregX.AND.dXregY) THEN 
      weights = delyatv(i:i+1,j:j+1)
   ELSEIF (dYregX.AND.dYregY) THEN
      weights = delxatu(i:i+1,j:j+1)
   ELSE
      weights = delxatu(i:i+1,j:j+1)*delyatv(i:i+1,j:j+1)
   ENDIF
   wsum = sflag(2,2)*weights(1,1)+sflag(1,2)*weights(2,1)+&
        & sflag(2,1)*weights(1,2)+sflag(1,1)*weights(2,2)
ENDIF
fdest = fdest.AND.wsum.GT.0.0

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (hreg) THEN

   IF (fdest) THEN
      Vvar_at_U = SUM(sflag*xin)/wsum
   ELSE
      Vvar_at_U = flagval
   ENDIF

!
!4.2 With weights
!----------------
!

ELSE

   IF (fdest) THEN
      Vvar_at_U = (sflag(1,1)*weights(2,2)*xin(1,1) + &
                 & sflag(2,1)*weights(1,2)*xin(2,1) + &
                 & sflag(1,2)*weights(2,1)*xin(1,2) + &
                 & sflag(2,2)*weights(1,1)*xin(2,2))/wsum
   ELSE
      Vvar_at_U = flagval
   ENDIF

ENDIF


RETURN

END FUNCTION Vvar_at_U

!========================================================================

FUNCTION Vvar_at_UV(xin,i,j,k,intsrce,intdest,outflag,hregular)
!************************************************************************
!
! *Vvar_at_UV* Interpolate a V-node array at UV-node grid point (i+1,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: i, intsrce, intdest, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: Vvar_at_UV

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER Lower X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 0 => all points
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 0 => all points
!                    = 1 => interior points and open boundaries only
!                    = 2 => interior points and UV-node open boundaries only
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the
!                     horizontal
!*Vvar_at_UV*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag, weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXregX
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE (intsrce)
   CASE (0); sflag = 1.0
   CASE (1); sflag = MERGE(1.0,0.0,nodeatv(i:i+1,j,k).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatv(i:i+1,j,k).EQ.2)
   CASE (3); sflag = MERGE(1.0,0.0,nodeatv(i:i+1,j,k).GT.1)
   CASE (4); sflag = MERGE(1.0,0.0,nodeatv(i:i+1,j,k).EQ.1.OR.&
                                 & nodeatv(i:i+1,j,k).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (0); fdest = .TRUE.
   CASE (1); fdest = nodeatuv(i+1,j,k).GT.0
   CASE (2); fdest = nodeatuv(i+1,j,k).EQ.1.OR.nodeatuv(i+1,j,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

IF (hreg) THEN
   wsum = SUM(sflag)
ELSE
   weights = delxatv(i:i+1,j)
   wsum = sflag(2)*weights(1) + sflag(1)*weights(2)
ENDIF

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (hreg) THEN

   IF (fdest.AND.wsum.GT.0.0) THEN
      Vvar_at_UV = SUM(sflag*xin)
   ELSE
      Vvar_at_UV = flagval
   ENDIF

!
!4.2 With weights
!----------------
!

ELSE

   IF (fdest.AND.wsum.GT.0.0) THEN
      Vvar_at_UV = (sflag(1)*weights(2)*xin(1)+&
                  & sflag(2)*weights(1)*xin(2))/wsum
   ELSE
      Vvar_at_UV = flagval
   ENDIF

ENDIF

RETURN

END FUNCTION Vvar_at_UV

!========================================================================

FUNCTION Vvar_at_VW(xin,i,j,k,intsrce,intdest,outflag,vregular)
!************************************************************************
!
! *Vvar_at_VW* Interpolate a V-node array at VW-node grid point (i,j,k+1)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: vregular
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: Vvar_at_VW

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Lower Z-index of interpolating array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*vregular*  LOGICAL Flag to select uniform or area averaging in the vertical
!*Vvar_at_VW*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag, weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(vregular)) THEN
   vreg = iopt_arrint_vreg.EQ.0.OR.dZregZ
ELSE
   vreg = vregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE (intsrce)
   CASE (1); sflag = MERGE(1.0,0.0,nodeatv(i,j,k:k+1).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatv(i,j,k:k+1).EQ.2)
   CASE (3); sflag = MERGE(1.0,0.0,nodeatv(i,j,k:k+1).GT.1)
   CASE (4); sflag = MERGE(1.0,0.0,nodeatv(i,j,k:k+1).EQ.1.OR.&
                                 & nodeatv(i,j,k:k+1).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (1); fdest = nodeatvw(i,j,k).GT.0
   CASE (2); fdest = nodeatvw(i,j,k).EQ.2
   CASE (3); fdest = nodeatvw(i,j,k).GT.1
   CASE (4); fdest = nodeatvw(i,j,k).EQ.1.OR.nodeatvw(i,j,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

IF (vreg) THEN
   wsum = SUM(sflag)
ELSE
   weights = delzatv(i,j,k:k+1)
   wsum = sflag(2)*weights(1) + sflag(1)*weights(2)
ENDIF

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (vreg) THEN

   Vvar_at_VW = MERGE(SUM(sflag*xin)/wsum,flagval,fdest)

!
!4.2 With weights
!----------------
!

ELSE

   IF (fdest) THEN
      Vvar_at_VW = (sflag(1)*weights(2)*xin(1)+sflag(2)*weights(1)*xin(2))/wsum
   ELSE
      Vvar_at_VW = flagval
   ENDIF

ENDIF


RETURN

END FUNCTION Vvar_at_VW

!========================================================================

FUNCTION Vvar_at_W(xin,i,j,k,intsrce,outflag,vregular)
!************************************************************************
!
! *Vvar_at_W* Interpolate a V-node array at W-node grid point (i,j,k+1)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: vregular
INTEGER, INTENT(IN) :: i, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2,2) :: xin
REAL :: Vvar_at_W

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at V-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Lower Z-index of interpolating array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*vregular*  LOGICAL Flag to select uniform or area averaging in the vertical
!*Vvar_at_W* REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: kk, k1
REAL, DIMENSION(2,2) :: sflag
REAL, DIMENSION(2,2) :: weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(vregular)) THEN
   vreg = iopt_arrint_vreg.EQ.0
ELSE
   vreg = vregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE (intsrce)
   CASE (1); sflag = MERGE(1.0,0.0,nodeatv(i,j:j+1,k:k+1).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatv(i,j:j+1,k:k+1).EQ.2)
   CASE (3); sflag = MERGE(1.0,0.0,nodeatv(i,j:j+1,k:k+1).GT.1)
   CASE (4); sflag = MERGE(1.0,0.0,nodeatv(i,j:j+1,k:k+1).EQ.1.OR.&
                                 & nodeatv(i,j:j+1,k:k+1).EQ.2)
END SELECT

!---destination
fdest = nodeatc(i,j).GT.0

!
!3. Weight factors
!-----------------
!

IF (vreg) THEN
   wsum = SUM(sflag)
ELSE
   kk_310: DO kk=1,2
      k1 = k+kk-1
      weights(:,kk) = delzatv(i,j:j+1,k1)
      IF (nodeatv(i,j+1,k1).GE.2.AND.nodeatc(i,j).GT.0) THEN
         weights(2,kk) = delzatc(i,j,k1)
      ENDIF
      IF (nodeatv(i,j,k1).GE.2.AND.nodeatc(i,j).GT.0) THEN
         weights(1,kk) = delzatc(i,j,k1)
      ENDIF
   ENDDO kk_310
   wsum = sflag(1,1)*weights(2,2)+sflag(2,1)*weights(1,2)+&
        & sflag(1,2)*weights(2,1)+sflag(2,2)*weights(1,1)
ENDIF

IF (wsum.LE.0.0) fdest = .FALSE. 

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (vreg) THEN

   IF (fdest) THEN
      Vvar_at_W = (sflag(1,1)*xin(1,1)+sflag(1,2)*xin(1,2)+&
                &  sflag(2,1)*xin(2,1)+sflag(2,2)*xin(2,2))/wsum
   ELSE
      Vvar_at_W = flagval
   ENDIF

!
!4.2 With weights
!----------------
!

ELSE

   IF (fdest) THEN
      Vvar_at_W = (sflag(1,1)*weights(2,2)*xin(1,1)+&
                 & sflag(2,1)*weights(1,2)*xin(2,1)+&
                 & sflag(1,2)*weights(2,1)*xin(1,2)+&
                 & sflag(2,2)*weights(1,1)*xin(2,2))/wsum
   ELSE
      Vvar_at_W = flagval
   ENDIF

ENDIF


RETURN

END FUNCTION Vvar_at_W

!========================================================================

FUNCTION VWvar_at_V(xin,i,j,k,intsrce,intdest,outflag)
!************************************************************************
!
! *VWvar_at_V* Interpolate a VW-node array at V-node grid point (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: VWvar_at_V

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*VWvar_at_V*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag


!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
SELECT CASE( intsrce)
   CASE (1); sflag = MERGE(1.0,0.0,nodeatvw(i,j,k:k+1).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatvw(i,j,k:k+1).EQ.2)
   CASE (3); sflag = MERGE(1.0,0.0,nodeatvw(i,j,k:k+1).GT.1)
   CASE (4); sflag = MERGE(1.0,0.0,nodeatvw(i,j,k:k+1).EQ.1.OR.&
                                 & nodeatvw(i,j,k:k+1).EQ.2)
END SELECT

!---destination
SELECT CASE (intdest)
   CASE (1); fdest = nodeatv(i,j,k).GT.0
   CASE (2); fdest = nodeatv(i,j,k).EQ.2
   CASE (3); fdest = nodeatv(i,j,k).GT.1
   CASE (4); fdest = nodeatv(i,j,k).EQ.1.OR.nodeatv(i,j,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

wsum = SUM(sflag)

!
!4. Interpolate
!--------------
!

IF (fdest.AND.wsum.GT.0.0) THEN
   VWvar_at_V = SUM(sflag*xin)/wsum
ELSE
   VWvar_at_V = flagval
ENDIF


RETURN

END FUNCTION VWvar_at_V

!========================================================================

FUNCTION VWvar_at_W(xin,i,j,k,intsrce,outflag)
!************************************************************************
!
! *VWvar_at_W* Interpolate a VW-node array at W-node grid point (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: VWvar_at_W

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Lower Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intsrce*   INTEGER Selects valid points at source node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*VWvar_at_W*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flag
!------------
!
!---source
SELECT CASE (intsrce)
   CASE (1); sflag = MERGE(1.0,0.0,nodeatvw(i,j:j+1,k).GT.0)
   CASE (2); sflag = MERGE(1.0,0.0,nodeatvw(i,j:j+1,k).EQ.2)
   CASE (3); sflag = MERGE(1.0,0.0,nodeatvw(i,j:j+1,k).GT.1)
   CASE (4); sflag = MERGE(1.0,0.0,nodeatvw(i,j:j+1,k).EQ.1.OR.&
                                 & nodeatvw(i,j:j+1,k).EQ.2)
END SELECT

!---destination
fdest = nodeatc(i,j).GT.0

!
!3. Weight factors
!-----------------
!

wsum = SUM(sflag)

!
!4. Interpolate
!--------------
!

IF (fdest.AND.wsum.GT.0.0) THEN
   VWvar_at_W = SUM(sflag*xin)/wsum
ELSE
   VWvar_at_W = flagval
ENDIF


RETURN

END FUNCTION VWvar_at_W

!========================================================================

FUNCTION Wvar_at_C(xin,i,j,outflag)
!************************************************************************
!
! *Wvar_at_C* Interpolate a W-node array at C-node grid point (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.0
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
INTEGER, INTENT(IN) :: i, j
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: Wvar_at_C

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*outflag*   REAL    Output flag for dry points
!*regular*   LOGICAL Flag to select uniform or area averaging
!*Wvar_at_C* REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Interpolate
!--------------
!

IF (nodeatc(i,j).GT.0) THEN
   Wvar_at_C = 0.5*SUM(xin)
ELSE
   Wvar_at_C = flagval
ENDIF


RETURN

END FUNCTION Wvar_at_C

!========================================================================

FUNCTION Wvar_at_U(xin,i,j,k,intdest,outflag,hregular)
!************************************************************************
!
! *Wvar_at_U* Interpolate a W-node array at U-node grid point (i+1,j,k)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: i, intdest, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2,2) :: xin
REAL :: Wvar_at_U

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER Lower X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the
!                    horizontal
!*Wvar_at_U* REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag, weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXregX
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!

!---source
sflag = MERGE(1.0,0.0,nodeatc(i:i+1,j).GT.0)

!---destination
SELECT CASE (intdest)
   CASE (1); fdest = nodeatu(i+1,j,k).GT.0
   CASE (2); fdest = nodeatu(i+1,j,k).EQ.2
   CASE (3); fdest = nodeatu(i+1,j,k).GT.1
   CASE (4); fdest = nodeatu(i+1,j,k).EQ.1.OR.nodeatu(i+1,j,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

IF (hreg) THEN
   wsum = 2.0*SUM(sflag)
ELSE
   weights = delxatc(i:i+1,j)
   wsum = 2.0*(sflag(2)*weights(1)+sflag(1)*weights(2))
ENDIF

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (hreg) THEN

   IF (fdest.AND.wsum.GT.0.0) THEN
      Wvar_at_U = (sflag(1)*(xin(1,1)+xin(1,2))+&
                 & sflag(2)*(xin(2,1)+xin(2,2)))/wsum
   ELSE
      Wvar_at_U = flagval
   ENDIF


!
!4.2 With weights
!----------------
!

ELSE

   IF (fdest.AND.wsum.GT.0.0) THEN
      Wvar_at_U = (sflag(1)*weights(2)*(xin(1,1)+xin(1,2))+&
                 & sflag(2)*weights(1)*(xin(2,1)+xin(2,2)))/wsum
   ELSE
      Wvar_at_U = flagval
   ENDIF

ENDIF


RETURN

END FUNCTION Wvar_at_U

!========================================================================

FUNCTION Wvar_at_UW(xin,i,j,k,intdest,outflag,hregular)
!************************************************************************
!
! *Wvar_at_UW* Interpolate a W-node array at UW-node grid point (i+1,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: i, intdest, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: Wvar_at_UW

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER Lower X-index of interpolating array
!*j*         INTEGER Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the
!                    horizontal
!*Wvar_at_UW*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag, weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dXregX
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
sflag = MERGE(1.0,0.0,nodeatc(i:i+1,j).GT.0)

!---destination
SELECT CASE (intdest)
   CASE (1); fdest = nodeatuw(i+1,j,k).GT.0
   CASE (2); fdest = nodeatuw(i+1,j,k).EQ.2
   CASE (3); fdest = nodeatuw(i+1,j,k).GT.1
   CASE (4); fdest = nodeatuw(i+1,j,k).EQ.1.OR.nodeatuw(i+1,j,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

IF (hreg) THEN
   wsum  = SUM(sflag)
ELSE
   weights = delxatc(i:i+1,j)
   wsum = sflag(2)*weights(1)+sflag(1)*weights(2)
ENDIF

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (hreg) THEN
   IF (fdest.AND.wsum.GT.0.0) THEN
      Wvar_at_UW = SUM(sflag*xin)/wsum
   ELSE
      Wvar_at_UW = flagval
   ENDIF

!
!4.2 With weights
!----------------
!

ELSE

   IF (fdest.AND.wsum.GT.0.0) THEN
      Wvar_at_UW = (sflag(1)*weights(2)*xin(1)+&
                  & sflag(2)*weights(1)*xin(2))/wsum
   ELSE
      Wvar_at_UW = flagval
   ENDIF

ENDIF


RETURN

END FUNCTION Wvar_at_UW

!========================================================================

FUNCTION Wvar_at_V(xin,i,j,k,intdest,outflag,hregular)
!************************************************************************
!
! *Wvar_at_V* Interpolate a W-node array at V-node grid point (i,j+1)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: i, intdest, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2,2) :: xin
REAL :: Wvar_at_V

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Lower Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the
!                    horizontal
!*Wvar_at_V* REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag, weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dYregY
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
sflag = MERGE(1.0,0.0,nodeatc(i,j:j+1).GT.0)

!---destination
SELECT CASE (intdest)
   CASE (1); fdest = nodeatv(i,j+1,k).GT.0
   CASE (2); fdest = nodeatv(i,j+1,k).EQ.2
   CASE (3); fdest = nodeatv(i,j+1,k).GT.1
   CASE (4); fdest = nodeatv(i,j+1,k).EQ.1.OR.nodeatv(i,j+1,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

IF (hreg) THEN
   wsum = 2.0*SUM(sflag)
ELSE
   weights = delyatc(i,j:j+1)
   wsum = 2.0*(sflag(2)*weights(1)+sflag(1)*weights(2))
ENDIF

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (hreg) THEN

   IF (fdest.AND.wsum.GT.0.0) THEN
      Wvar_at_V = (sflag(1)*(xin(1,1)+xin(1,2))+&
                 & sflag(2)*(xin(2,1)+xin(2,2)))/wsum
   ELSE
      Wvar_at_V = flagval
   ENDIF

!
!4.2 With weights
!----------------
!

ELSE

   IF (fdest.AND.wsum.GT.0.0) THEN
      Wvar_at_V = (sflag(1)*weights(2)*(xin(1,1)+xin(1,2))+&
                 & sflag(2)*weights(1)*(xin(2,1)+xin(2,2)))/wsum
   ELSE
      Wvar_at_V = flagval
   ENDIF

ENDIF


RETURN

END FUNCTION Wvar_at_V

!========================================================================

FUNCTION Wvar_at_VW(xin,i,j,k,intdest,outflag,hregular)
!************************************************************************
!
! *Wvar_at_VW* Interpolate a W-node array at VW-node grid point (i,j+1)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)array_interp.f90  V2.2
!
! Description -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars

!
!*  Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: i, intdest, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xin
REAL :: Wvar_at_VW

!
! Name       Type    Purpose
!-----------------------------------------------------------------------------
!*xin*       REAL    Input array at U-nodes
!*i*         INTEGER X-index of interpolating array
!*j*         INTEGER Lower Y-index of interpolating array
!*k*         INTEGER Vertical index of interpolating array or nz for a 2-D 
!                    array
!*intdest*   INTEGER Selects valid points at destination node
!                    = 1 => coastal boundaries, interior points and open
!                           boundaries only
!                    = 2 => interior points only
!                    = 3 => interior points and open boundaries only
!                    = 4 => coastal boundaries and interior points only
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the
!                    horizontal
!*Wvar_at_VW*REAL    Interpolated value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, DIMENSION(2) :: sflag, weights


!
!1. Optional arguments
!---------------------
!

IF (.NOT.PRESENT(hregular)) THEN
   hreg = iopt_arrint_hreg.EQ.0.OR.dYregY
ELSE
   hreg = hregular
ENDIF

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Data flags
!-------------
!
!---source
sflag = MERGE(1.0,0.0,nodeatc(i,j:j+1).GT.0)

!---destination
SELECT CASE (intdest)
   CASE (1); fdest = nodeatvw(i,j+1,k).GT.0
   CASE (2); fdest = nodeatvw(i,j+1,k).EQ.2
   CASE (3); fdest = nodeatvw(i,j+1,k).GT.1
   CASE (4); fdest = nodeatvw(i,j+1,k).EQ.1.OR.nodeatvw(i,j+1,k).EQ.2
END SELECT

!
!3. Weight factors
!-----------------
!

IF (hreg) THEN
   wsum  = SUM(sflag)
ELSE
   weights = delyatc(i,j:j+1)
   wsum = sflag(2)*weights(1)+sflag(1)*weights(2)
ENDIF

!
!4. Interpolate
!--------------
!
!4.1 Without weights
!-------------------
!

IF (hreg) THEN
   IF (fdest.AND.wsum.GT.0) THEN
      Wvar_at_VW = SUM(sflag*xin)/wsum
   ELSE
      Wvar_at_VW = flagval
   ENDIF

!
!4.2 With weights
!----------------
!

ELSE

   IF (fdest.AND.wsum.GT.0.0) THEN
      Wvar_at_VW = (sflag(1)*weights(2)*xin(1)+&
                  & sflag(2)*weights(1)*xin(2))/wsum
   ELSE
      Wvar_at_VW = flagval
   ENDIF

ENDIF


RETURN

END FUNCTION Wvar_at_VW


END MODULE array_interp
