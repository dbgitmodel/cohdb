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
! *Corrector_Terms* Ensemble of routines for evaluation of corrector terms
!                   in transport equations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Corrector_Terms.f90  V2.11.2
!
! $Date: 2016-10-10 16:36:15 +0200 (Mon, 10 Oct 2016) $
!
! $Revision: 973 $
!
! Description - 
!
! Reference -
!
! Routines - Xcorr_at_C, Ycorr_at_C, Zcorr_at_C,
!            Xcorr_at_W, Ycorr_at_W, Zcorr_at_W  
!
!************************************************************************
!

!========================================================================

SUBROUTINE Xcorr_at_C(ucorratc,klo,kup)
!************************************************************************
!
! *Xcorr_at_C* Corrector term in X-direction for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Corrector_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_C_3d, transport_at_C_4d1, transport_at_C_4d2
!
! Module calls - error_alloc
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,klo:kup) :: ucorratc

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ucorratc*  REAL    Corrector term                                        [-]
!*klo*       INTEGER Lower vertical array bound
!*kup*       INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, npcc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, array2du
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dc


procname(pglev+1) = 'Xcorr_at_C'
CALL log_timer_in(npcc)

!
!1. Allocate
!-----------
!

ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2du(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('array2du',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (array3dc(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('array3dc',3,(/ncloc,nrloc,nz/),kndrtype)

!
!2. Work space arrays
!--------------------
!

IF (dYregX) THEN
   WHERE (maskatc_int)
      array2dc = delxatc(1:ncloc,1:nrloc)/delt3d
   END WHERE
ELSE
   WHERE (maskatc_int)
      array2dc = garea/delt3d
   END WHERE
ENDIF

k_210: DO k=klo,kup
   WHERE (maskatc_int)
      array3dc(:,:,k) = array2dc*delzatc(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_210

!
!3. Corrector term
!-----------------
!

IF (dYregX) THEN
   k_310: DO k=klo,kup
      WHERE (nodeatu(1:ncloc+1,1:nrloc,k).GT.0)
         array2du = delzatu(1:ncloc+1,:,k)*ufvel(1:ncloc+1,:,k)
      END WHERE
      WHERE (maskatc_int)
         ucorratc(:,:,k) = (array2du(2:ncloc+1,:)-array2du(1:ncloc,:))&
                              & /array3dc(:,:,k)
      END WHERE
   ENDDO k_310

ELSE
   k_320: DO k=klo,kup
      WHERE (nodeatu(1:ncloc+1,1:nrloc,k).GT.0)
         array2du = delyatu(1:ncloc+1,1:nrloc)*delzatu(1:ncloc+1,:,k)*&
                  & ufvel(1:ncloc+1,:,k)
      END WHERE
      WHERE (maskatc_int)
         ucorratc(:,:,k) = (array2du(2:ncloc+1,:)-array2du(1:ncloc,:))&
                         & /array3dc(:,:,k)
      END WHERE
   ENDDO k_320

ENDIF

!
!4. Deallocate
!-------------
!

DEALLOCATE (array2dc,array2du,array3dc)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Xcorr_at_C

!========================================================================

SUBROUTINE Ycorr_at_C(vcorratc,klo,kup)
!************************************************************************
!
! *Ycorr_at_C* Corrector term in Y-direction for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Corrector_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_C_3d, transport_at_C_4d1, transport_at_C_4d2
!
! Module calls - error_alloc
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,klo:kup) :: vcorratc

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*vcorratc*  REAL    Corrector term                                        [-]
!*klo*       INTEGER Lower vertical array bound
!*kup*       INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, npcc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, array2dv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dc


procname(pglev+1) = 'Ycorr_at_C'
CALL log_timer_in(npcc)

!
!1. Allocate
!-----------
!

ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('array2dv',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (array3dc(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('array3dc',3,(/ncloc,nrloc,nz/),kndrtype)

!
!2. Work space arrays
!--------------------
!

IF (dXregY) THEN
   WHERE (maskatc_int)
      array2dc = delyatc(1:ncloc,1:nrloc)/delt3d
   END WHERE
ELSE
   WHERE (maskatc_int)
      array2dc = garea/delt3d
   END WHERE
ENDIF

k_210: DO k=klo,kup
   WHERE (maskatc_int)
      array3dc(:,:,k) = array2dc*delzatc(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_210

!
!3. Corrector term
!-----------------
!

IF (dXregY) THEN
   k_310: DO k=klo,kup
      WHERE (nodeatv(1:ncloc,1:nrloc+1,k).GT.0)
         array2dv = delzatv(:,1:nrloc+1,k)*vfvel(:,1:nrloc+1,k)
      END WHERE
      WHERE (maskatc_int)
         vcorratc(:,:,k) = (array2dv(:,2:nrloc+1)-array2dv(:,1:nrloc))&
                              & /array3dc(:,:,k)
      END WHERE
   ENDDO k_310

ELSE
   k_320: DO k=klo,kup
      WHERE (nodeatv(1:ncloc,1:nrloc+1,k).GT.0)
         array2dv = delxatv(1:ncloc,1:nrloc+1)*delzatv(:,1:nrloc+1,k)*&
                  & vfvel(:,1:nrloc+1,k)
      END WHERE
      WHERE (maskatc_int)
         vcorratc(:,:,k) = (array2dv(:,2:nrloc+1)-array2dv(:,1:nrloc))&
                         & /array3dc(:,:,k)
     END WHERE
   ENDDO k_320

ENDIF
   
!
!4. Deallocate
!-------------
!

DEALLOCATE (array2dc,array2dv,array3dc)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Ycorr_at_C

!========================================================================

SUBROUTINE Zcorr_at_C(wcorratc,klo,kup)
!************************************************************************
!
! *Zcorr_at_C* Corrector term in Z-direction for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Corrector_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_C_3d, transport_at_C_4d1, transport_at_C_4d2
!
! Module calls - error_alloc
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,klo:kup) :: wcorratc

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*wcorratc*  REAL    Corrector term                                        [-]
!*klo*       INTEGER Lower vertical array bound
!*kup*       INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, npcc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dc


procname(pglev+1) = 'Zcorr_at_C'
CALL log_timer_in(npcc)

!
!1. Allocate
!-----------
!

ALLOCATE (array3dc(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('array3dc',3,(/ncloc,nrloc,nz/),kndrtype)

!
!2. Work space array
!-------------------
!

k_210: DO k=klo,kup
   WHERE (maskatc_int)
      array3dc(:,:,k) = delzatc(1:ncloc,1:nrloc,k)/delt3d
   END WHERE
ENDDO k_210

!
!3. Corrector term
!-----------------
!

k_310: DO k=klo,kup
   WHERE (maskatc_int)
      wcorratc(:,:,k) = (wfvel(1:ncloc,1:nrloc,k+1)-wfvel(1:ncloc,1:nrloc,k))&
                      & /array3dc(:,:,k)
   END WHERE
ENDDO k_310

!
!4. Deallocate
!-------------
!

DEALLOCATE (array3dc)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Zcorr_at_C

!========================================================================

SUBROUTINE Xcorr_at_W(ucorratw,uvelatuw,klo,kup)
!************************************************************************
!
! *Xcorr_at_W* Corrector term in X-direction for a quantity at W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Corrector_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_W
!
! Module calls - error_alloc
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,klo:kup) :: ucorratw
REAL, INTENT(IN), DIMENSION(2-nhalo:ncloc+nhalo,nrloc,2:nz) :: uvelatuw

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ucorratw*  REAL    Corrector term                                        [-]
!*uvelatuw*  REAL    U-velocity at UW-nodes                              [m/s]
!*klo*       INTEGER Lower vertical array bound
!*kup*       INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, npcc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, array2du


procname(pglev+1) = 'Xcorr_at_W'
CALL log_timer_in(npcc)

!
!1. Allocate
!-----------
!

ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2du(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('array2du',2,(/ncloc+1,nrloc/),kndrtype)

!
!2. Initialise
!-------------
!

array2du = 0.0

!
!3. Regular grid
!---------------
!

IF (dYregX) THEN
   WHERE (maskatc_int)
      array2dc = delxatc(1:ncloc,1:nrloc)/delt3d
   END WHERE
   k_310: DO k=klo,kup
      WHERE (nodeatuw(1:ncloc+1,1:nrloc,k).GT.1)
         array2du = delzatuw(:,:,k)*uvelatuw(1:ncloc+1,:,k)
      END WHERE
      WHERE (maskatc_int)
         ucorratw(:,:,k) = (array2du(2:ncloc+1,:)-array2du(1:ncloc,:))&
                        & /(array2dc*delzatw(1:ncloc,1:nrloc,k))
      END WHERE
   ENDDO k_310

!
!4. Irregular grid
!-----------------
!

ELSE
   WHERE (maskatc_int)
      array2dc = garea/delt3d
   END WHERE
   k_410: DO k=klo,kup
      WHERE (nodeatuw(1:ncloc+1,1:nrloc,k).GT.1)
         array2du = delyatu(1:ncloc+1,1:nrloc)*delzatuw(:,:,k)*&
                  & uvelatuw(1:ncloc+1,:,k)
      END WHERE
      WHERE (maskatc_int)
         ucorratw(:,:,k) = (array2du(2:ncloc+1,:)-array2du(1:ncloc,:))&
                        & /(array2dc*delzatw(1:ncloc,1:nrloc,k))
      END WHERE
   ENDDO k_410

ENDIF

!
!5. Deallocate
!-------------
!

DEALLOCATE (array2dc,array2du)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Xcorr_at_W

!========================================================================

SUBROUTINE Ycorr_at_W(vcorratw,vvelatvw,klo,kup)
!************************************************************************
!
! *Ycorr_at_W* Corrector term in Y-direction for a quantity at W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Corrector_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_W
!
! Module calls - error_alloc
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,klo:kup) :: vcorratw
REAL, INTENT(IN), DIMENSION(ncloc,2-nhalo:nrloc+nhalo,2:nz) :: vvelatvw

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*vcorratw*  REAL    Corrector term                                        [-]
!*vvelatvw*  REAL    V-velocity at VW-nodes                              [m/s]
!*klo*       INTEGER Lower vertical array bound
!*kup*       INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, npcc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, array2dv


procname(pglev+1) = 'Ycorr_at_W'
CALL log_timer_in(npcc)

!
!1. Allocate
!-----------
!

ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('array2dv',2,(/ncloc,nrloc+1/),kndrtype)

!
!2. Initialise
!-------------
!

array2dv = 0.0

!
!3. Regular grid
!---------------
!

IF (dXregY) THEN
   WHERE (maskatc_int)
      array2dc = delyatc(1:ncloc,1:nrloc)/delt3d
   END WHERE
   k_310: DO k=klo,kup
      WHERE (nodeatvw(1:ncloc,1:nrloc+1,k).GT.1)
         array2dv = delzatvw(:,:,k)*vvelatvw(:,1:nrloc+1,k)
      END WHERE
      WHERE (maskatc_int)
         vcorratw(:,:,k) = (array2dv(:,2:nrloc+1)-array2dv(:,1:nrloc))&
                        & /(array2dc*delzatw(1:ncloc,1:nrloc,k))
      END WHERE
   ENDDO k_310

!
!4. Irregular grid
!-----------------
!

ELSE
   WHERE (maskatc_int)
      array2dc = garea/delt3d
   END WHERE
   k_410: DO k=klo,kup
      WHERE (nodeatvw(1:ncloc,1:nrloc+1,k).GT.1)
         array2dv = delxatv(1:ncloc,1:nrloc+1)*delzatvw(:,:,k)*&
                  & vvelatvw(:,1:nrloc+1,k)
      END WHERE
      WHERE (maskatc_int)
         vcorratw(:,:,k) = (array2dv(:,2:nrloc+1)-array2dv(:,1:nrloc))&
                        & /(array2dc*delzatw(1:ncloc,1:nrloc,k))
      END WHERE
   ENDDO k_410

ENDIF

!
!5. Deallocate
!-------------
!

DEALLOCATE (array2dc,array2dv)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Ycorr_at_W

!========================================================================

SUBROUTINE Zcorr_at_W(wcorratw,wvelatc,klo,kup)
!************************************************************************
!
! *Zcorr_at_W* Corrector term in Z-direction for a quantity at W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Corrector_Terms.f90  V2.10.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_W
!
! Module calls -
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,klo:kup) :: wcorratw
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz) :: wvelatc

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*wcorratw*  REAL    Corrector term                                        [-]
!*wvelatc*   REAL    Vertical velocity at W-nodes                        [m/s]
!*klo*       INTEGER Lower vertical array bound
!*kup*       INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, npcc


procname(pglev+1) = 'Zcorr_at_W'
CALL log_timer_in(npcc)

k_110: DO k=klo,kup
   WHERE (maskatc_int)
      wcorratw(:,:,k) = (wvelatc(:,:,k+1)-wvelatc(:,:,k))&
                      & /delzatw(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_110

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Zcorr_at_W
