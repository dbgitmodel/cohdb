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

MODULE turbulence_routines
!************************************************************************
!
! *turbulence_routines* Utility routines for turbulence closures
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)turbulence_routines.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Routines - mom_strat_MA, richardson_number, richardson_number_0d,
!            scal_strat_MA
!
!************************************************************************
!
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

LOGICAL :: info = .FALSE.

CONTAINS

!=======================================================================

FUNCTION mom_strat_MA(ricnum,mask)
!************************************************************************
!
! *mom_strat_MA* Stratification factor for momentum in Munk-Anderson relations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)turbulence_routines.f90  V2.11.2
!
! Description - 
!
! Module calls - error_alloc
!
!************************************************************************
!
USE gridpars
USE syspars
USE turbpars
USE error_routines, ONLY: error_alloc

!
!*Arguments
!

LOGICAL, INTENT(IN), DIMENSION(ncloc,nrloc) :: mask
REAL, INTENT(IN), DIMENSION(ncloc,nrloc) :: ricnum
REAL, DIMENSION(ncloc,nrloc) :: mom_strat_MA

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*ricnum*       REAL    Richardson number
!*mask*         LOGICAL .TRUE. for wet cells
!*mom_strat_MA* REAL    Stratification factor for momentum
!
!************************************************************************
!
!*Local variables
!
INTEGER :: i, j
REAL :: epsmin
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: denom


procname(pglev+1) = 'mom_strat_MA'
CALL log_timer_in(info=info)

!---allocate array
ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)

!---stratification factor
epsmin = vmaxmom_ma**(-1.0/expmom_ma)
j_110: DO j=1,nrloc
i_110: DO i=1,ncloc
   IF (mask(i,j)) THEN
      denom(i,j) = MAX(1.0+alpha_ma*ricnum(i,j),epsmin)
      mom_strat_MA(i,j) = 1.0/(denom(i,j)**expmom_ma)
   ENDIF
ENDDO i_110
ENDDO j_110

!---deallocate array
DEALLOCATE (denom)

CALL log_timer_out(info=info)


RETURN

END FUNCTION mom_strat_MA

!========================================================================

FUNCTION richardson_number(k,mask)
!************************************************************************
!
! *richardson_number* Richardson number
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)turbulence_routines.f90  V2.1.2
!
! Description - 
!
! Module calls -
!
!************************************************************************
!
USE gridpars
USE turbulence

!
!*Arguments
!
INTEGER, INTENT(IN) :: k
LOGICAL, INTENT(IN), DIMENSION(:,:) :: mask
REAL, DIMENSION(LBOUND(mask,1):UBOUND(mask,1),&
              & LBOUND(mask,2):UBOUND(mask,2)) :: richardson_number

!
! Name               Type    Purpose
!------------------------------------------------------------------------------
!*k*                 INTEGER Vertical grid index
!*mask*              LOGICAL .TRUE. for wet cells
!*richardson_number* REAL    Richardson number
!
!************************************************************************
!
!*Local variables
!
INTEGER :: i, j, l1, l2, u1, u2
REAL, PARAMETER :: epstol = 1.0E-20, ricmax = 1.0E+03


procname(pglev+1) = 'richardson_number'
CALL log_timer_in(info=info)

l1 = LBOUND(mask,1)
l2 = LBOUND(mask,2)
u1 = UBOUND(mask,1)
u2 = UBOUND(mask,2)

IF ((k.EQ.1).OR.(k.EQ.(nz+1))) THEN
   WHERE (mask) richardson_number = 0.0
ELSE
   j_110: DO j=l2,u2
   i_110: DO i=l1,u1
      IF (mask(i,j)) THEN
         richardson_number(i,j) = MERGE(ricmax*SIGN(1.0,buofreq2(i,j,k)),&
              & buofreq2(i,j,k)/shearfreq2(i,j,k),&
              & (ABS(buofreq2(i,j,k)).GT.shearfreq2(i,j,k)*ricmax).OR.&
              & (shearfreq2(i,j,k).LT.epstol))
      ENDIF
   ENDDO i_110
   ENDDO j_110
ENDIF

CALL log_timer_out(info=info)


RETURN

END FUNCTION richardson_number

!========================================================================

FUNCTION richardson_number_0d(i,j,k)
!************************************************************************
!
! *richardson_number_0d* Richardson number at a specific grid location
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)turbulence_routines.f90  V2.1.2
!
! Description - 
!
! Module calls -
!
!************************************************************************
!
USE gridpars
USE turbulence

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, k
REAL :: richardson_number_0d

!
! Name               Type    Purpose
!------------------------------------------------------------------------------
!*i*                 INTEGER X-index on the model grid
!*j*                 INTEGER Y-index on the model grid
!*k*                 INTEGER Vertical grid index
!*richardson_number* REAL    Richardson number
!
!************************************************************************
!
!*Local variables
!
REAL, PARAMETER :: epstol = 1.0E-20, ricmax = 1.0E+03


IF ((k.EQ.1).OR.(k.EQ.(nz+1))) THEN
   richardson_number_0d = 0.0
ELSE
   richardson_number_0d = MERGE(ricmax*SIGN(1.0,buofreq2(i,j,k)),&
        & buofreq2(i,j,k)/shearfreq2(i,j,k),&
        & (ABS(buofreq2(i,j,k)).GT.shearfreq2(i,j,k)*ricmax).OR.&
        & (shearfreq2(i,j,k).LT.epstol))
ENDIF


RETURN

END FUNCTION richardson_number_0d

!========================================================================

FUNCTION scal_strat_MA(ricnum,mask)
!************************************************************************
!
! *scal_strat_MA* Stratification factor for scalars in Munk-Anderson relations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)turbulence_routines.f90  V2.11.2
!
! Description - 
!
! Module calls - error_alloc
!
!************************************************************************
!
USE gridpars
USE syspars
USE turbpars
USE error_routines, ONLY: error_alloc

!
!*Arguments
!
LOGICAL, INTENT(IN), DIMENSION(ncloc,nrloc) :: mask
REAL, INTENT(IN), DIMENSION(ncloc,nrloc) :: ricnum
REAL, DIMENSION(ncloc,nrloc) :: scal_strat_MA

!
! Name           Type    Purpose
!------------------------------------------------------------------------------
!*ricnum*        REAL    Richardson number
!*mask*          LOGICAL .TRUE. for wet cells
!*scal_strat_MA* REAL    Stratification factor for scalars
!
!************************************************************************
!
!*Local variables
!
INTEGER :: i, j
REAL :: epsmin
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: denom


procname(pglev+1) = 'scal_strat_MA'
CALL log_timer_in(info=info)

!---allocate array
ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)

!---stratification factor
epsmin = vmaxscal_ma**(-1.0/expscal_ma)
j_110: DO j=1,nrloc 
i_110: DO i=1,ncloc
   IF (mask(i,j)) THEN
      denom(i,j) = MAX(1.0+beta_ma*ricnum(i,j),epsmin)
      scal_strat_MA(i,j) = 1.0/(denom(i,j)**expscal_ma)
   ENDIF
ENDDO i_110
ENDDO j_110

!---deallocate array
DEALLOCATE (denom)

CALL log_timer_out(info=info)


RETURN

END FUNCTION scal_strat_MA


END MODULE turbulence_routines
