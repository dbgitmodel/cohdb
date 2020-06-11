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
! *MultiGrid_Schemes* Series of routines for the multigrid solution procedure
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Schemes.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
! Reference -
!
! Routines - check_convergence, construct_mg_matrix, mg_cycles, mg_smoothers,
!            multi_grid_scheme, prolongate, restrict, update_residual,
!            update_solution
!
!************************************************************************
!

!========================================================================

SUBROUTINE check_convergence(iteration,maxiterations,mglevel)
!************************************************************************
!
! *check_convergence* check convergence of the iteration process
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Schemes.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - multi_grid_scheme
!
! External calls - update_residual
!
! Module calls - error_alloc, sum_vars
!
!************************************************************************
!
USE iopars
USE multigrid
USE physpars
USE syspars
USE error_routines, ONLY: error_alloc
USE paral_utilities, ONLY: sum_vars
!USE paral_utilities, ONLY: sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iteration, maxiterations, mglevel

!
! Name           Type    Purpose
!-----------------------------------------------------------------------------
!*iteration*     INTEGER Iteration index
!*maxiterations* INTEGER Number of iterations
!*mglevel*       INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: citer, clev
INTEGER :: lev, nx, ny
INTEGER, DIMENSION(2) :: maxpos
REAL :: residual_log, reference_sum, residual_sum
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d


lev = mglevel
IF ((lev.EQ.0.AND.(iteration.LT.maxiterations.AND.monlog.LT.2)).OR.&
  & (lev.GT.0.AND.monlog.LT.3)) RETURN

procname(pglev+1) = 'check_convergence'
CALL log_timer_in()

!
!1. Allocate
!-----------
!

nx = mgvars(lev)%ncloc; ny = mgvars(lev)%nrloc
ALLOCATE (array2d(nx,ny),STAT=errstat)
CALL error_alloc('array2d',2,(/nx,ny/),kndrtype)

!
!2. Calculate the residual
!-------------------------
!

CALL update_residual(0,lev)

!
!3. Residual norm
!----------------
!

array2d = ABS(mgvars(lev)%residual(1:nx,1:ny))
CALL sum_vars(array2d,residual_sum,0,commall=.TRUE.,&
            & mask=mgvars(lev)%maskatc_int)
!CALL sum2_vars(array2d,residual_sum,(/0,0,0,0/),'C  ',0,mglevel=lev,&
!            & commall=.TRUE.,mask=mgvars(lev)%maskatc_int)
array2d = ABS(mgvars(lev)%subdiags(:,:,0)*mgvars(lev)%solvec(1:nx,1:ny))
maxpos = MAXLOC(mgvars(lev)%residual(1:nx,1:ny),mask=mgvars(lev)%maskatc_int)
CALL sum_vars(array2d,reference_sum,0,commall=.TRUE.,&
            & mask=mgvars(lev)%maskatc_int)
!CALL sum2_vars(array2d,reference_sum,(/0,0,0,0/),'C  ',0,mglevel=lev,&
!             & commall=.TRUE.,mask=mgvars(lev)%maskatc_int)
IF (residual_sum.EQ.0.0) THEN
   residual_norm = 0.0
ELSE
   residual_norm = residual_sum/reference_sum
ENDIF
residual_log = LOG10(ABS(residual_norm))

!
!4. Write monitoring info
!------------------------
!

IF (monflag) THEN
   IF ((lev.EQ.0.AND.monlog.GT.1).OR.(lev.GT.0.AND.monlog.EQ.3)) THEN
      WRITE (clev,'(I12)') lev; clev = ADJUSTL(clev)
      WRITE (citer,'(I12)') iteration; citer = ADJUSTL(citer)
      WRITE (iomon,9001) 'L'//TRIM(clev)//': ', iteration, residual_log, maxpos
   ENDIF
ENDIF

!
!5. Deallocate
!-------------
!

DEALLOCATE (array2d)

CALL log_timer_out()


RETURN

9001 FORMAT (A,I4,1X,G15.7,2I4)

END SUBROUTINE check_convergence

!========================================================================

SUBROUTINE construct_mg_matrix(itsweep,mglevel)
!************************************************************************
!
! *construct_mg_matrix* update solution from previous iteration
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Schemes.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program -
!
! External calls - open_boundary_conds_impl, open_boundary_conds_impl_level
!
! Module calls - exchange_mod
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE iopars
USE meteo
USE multigrid
USE physpars
USE structures
USE switches
USE syspars
USE timepars
USE wavevars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: itsweep, mglevel

!
! Name      Type    Purpose
!-----------------------------------------------------------------------------
!*itsweep*  INTEGER Defines which matrix elements in the matrix multiplication
!            = 0 => all elements (Jacobi method)
!            = 1 => "red" points only (Gauus-Seidel) method
!            = 2 => "black" points only (Gauus-Seidel) method
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j, lev, n, nx, nloc, ny
REAL :: fac
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu, maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, array2du, array2dv


IF (itsweep.EQ.2.OR.mgvars(mglevel)%linearization_done) RETURN

procname(pglev+1)= 'construct_mg_matrix'
CALL log_timer_in()

lev = mglevel
nx = mgvars(lev)%ncloc; ny = mgvars(lev)%nrloc

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatu(nx+1,ny),STAT=errstat)
CALL error_alloc('maskatu',2,(/nx+1,ny/),kndlog)

ALLOCATE (maskatv(nx,ny+1),STAT=errstat)

CALL error_alloc('maskatv',2,(/nx,ny+1/),kndlog)

ALLOCATE (array2dc(nx,ny),STAT=errstat)
CALL error_alloc('array2dc',2,(/nx,ny/),kndrlong)

ALLOCATE (array2du(nx+1,ny),STAT=errstat)
CALL error_alloc('array2du',2,(/nx+1,ny/),kndrtype)

ALLOCATE (array2dv(nx,ny+1),STAT=errstat)
CALL error_alloc('array2dv',2,(/nx,ny+1/),kndrtype)

!
!2. Initialise
!-------------
!
!---mask arrays
maskatu = mgvars(lev)%node2du(:,1:ny).EQ.2
maskatv = mgvars(lev)%node2dv(1:nx,:).EQ.2

!---work space arrays
fac = theta_sur*delt3d
array2dc = mgvars(lev)%delxatc(1:nx,:)*mgvars(lev)%delyatc(:,1:ny)/delt3d
WHERE (maskatu)
   array2du = fac*mgvars(lev)%gaccatu(1:nx+1,:)*&
            & mgvars(lev)%deptotatu(1:nx+1,1:ny)*&
            & mgvars(lev)%delyatu/mgvars(lev)%delxatu(1:nx+1,:)
ELSEWHERE
   array2du = 0.0
END WHERE
WHERE (maskatv)
   array2dv = fac*mgvars(lev)%gaccatv(:,1:ny+1)* &
            & mgvars(lev)%deptotatv(1:nx,1:ny+1)*&
            & mgvars(lev)%delxatv/mgvars(lev)%delyatv(:,1:ny+1)
ELSEWHERE
   array2dv = 0.0
END WHERE

!
!3. Build the coefficient matrix
!-------------------------------
!
!3.1 Time derivative
!-------------------
!

WHERE (mgvars(lev)%maskatc_int)
   mgvars(lev)%subdiags(:,:,0) = array2dc
END WHERE

IF (lev.EQ.0) THEN
   WHERE (maskatc_int)
      mgvars(0)%rhs = array2dc*(zeta_old(1:nx,1:ny)-zeta(1:nx,1:ny))
   END WHERE
ENDIF

!
!3.2 Fluxes at U-nodes
!---------------------
!
!3.2.1 East
!----------
!

WHERE (maskatu(2:nx+1,:))
   mgvars(lev)%subdiags(:,:,1) = -array2du(2:nx+1,:)
   mgvars(lev)%subdiags(:,:,0) = mgvars(lev)%subdiags(:,:,0) + &
                                     & array2du(2:nx+1,:)
END WHERE

IF (lev.EQ.0) THEN
   WHERE (maskatu(2:nx+1,:))
      mgvars(0)%rhs = mgvars(0)%rhs - &
               & delyatu(2:nx+1,1:ny)*deptotatu(2:nx+1,1:ny)*umpred(2:nx+1,1:ny)
   END WHERE
   IF (iopt_waves_curr.EQ.1) THEN
      WHERE (maskatu(2:nx+1,:))
         mgvars(0)%rhs = mgvars(0)%rhs - &
               & delyatu(2:nx+1,1:ny)*deptotatu(2:nx+1,1:ny)*&
               & umstokesatu(2:nx+1,1:ny)
      END WHERE
   ENDIF
ENDIF

!
!3.2.2 West
!----------
!

WHERE (maskatu(1:nx,:))
   mgvars(lev)%subdiags(:,:,-2) = -array2du(1:nx,:)
   mgvars(lev)%subdiags(:,:,0) = mgvars(lev)%subdiags(:,:,0) + &
                                     & array2du(1:nx,:)
END WHERE

IF (lev.EQ.0) THEN
   WHERE (maskatu(1:nx,:))
      mgvars(0)%rhs = mgvars(0)%rhs + &
                & delyatu(1:nx,1:ny)*deptotatu(1:nx,1:ny)*umpred(1:nx,1:ny)
   END WHERE
   IF (iopt_waves_curr.EQ.1) THEN
      WHERE (maskatu(1:nx,:))
         mgvars(0)%rhs = mgvars(0)%rhs + &
              & delyatu(1:nx,1:ny)*deptotatu(1:nx,1:ny)*umstokesatu(1:nx,1:ny)
      END WHERE
   ENDIF
ENDIF

!
!3.3 Fluxes at V-nodes
!---------------------
!
!3.3.1 North
!-----------
!

WHERE (maskatv(:,2:ny+1))
   mgvars(lev)%subdiags(:,:,2) = -array2dv(:,2:ny+1)
   mgvars(lev)%subdiags(:,:,0) = mgvars(lev)%subdiags(:,:,0) + &
                                     & array2dv(:,2:ny+1) 
END WHERE

IF (lev.EQ.0) THEN
   WHERE (maskatv(:,2:ny+1))
      mgvars(0)%rhs = mgvars(0)%rhs - &
               & delxatv(1:nx,2:ny+1)*deptotatv(1:nx,2:ny+1)*vmpred(1:nx,2:ny+1)
   END WHERE
   IF (iopt_waves_curr.EQ.1) THEN
      WHERE (maskatv(:,2:ny+1))
         mgvars(0)%rhs = mgvars(0)%rhs - &
               & delxatv(1:nx,2:ny+1)*deptotatv(1:nx,2:ny+1)*&
               & vmstokesatv(1:nx,2:ny+1)
      END WHERE
   ENDIF
ENDIF

!
!3.3.2 South
!-----------
!

WHERE (maskatv(:,1:ny))
   mgvars(lev)%subdiags(:,:,-1) = -array2dv(:,1:ny)
   mgvars(lev)%subdiags(:,:,0) = mgvars(lev)%subdiags(:,:,0) + &
                                    &  array2dv(:,1:ny)
END WHERE

IF (lev.EQ.0) THEN
   WHERE (maskatv(:,1:ny))
      mgvars(0)%rhs = mgvars(0)%rhs + &
                & delxatv(1:nx,1:ny)*deptotatv(1:nx,1:ny)*vmpred(1:nx,1:ny)
   END WHERE
   IF (iopt_waves_curr.EQ.1) THEN
      WHERE (maskatv(:,1:ny))
         mgvars(0)%rhs = mgvars(0)%rhs + &
           & delxatv(1:nx,1:ny)*deptotatv(1:nx,1:ny)*vmstokesatv(1:nx,1:ny)
      END WHERE
   ENDIF
ENDIF

!
!3.4 Add discharges
!------------------
!

IF (lev.EQ.0.AND.iopt_dischr.EQ.1.AND.lev.EQ.0) THEN
   nloc_340: DO nloc=1,numdisloc
      n = indexdisloc(nloc)
      IF (disflag(n)) THEN
         i = idisloc(nloc); j = jdisloc(nloc)
         mgvars(0)%rhs(i,j) = mgvars(0)%rhs(i,j) + disspeed(n)*&
                            & mgvars(0)%delxatc(i,j)*mgvars(lev)%delyatc(i,j)
      ENDIF
   ENDDO nloc_340
ENDIF

!
!3.5 Add precipitation minus evaporation
!---------------------------------------
!

IF (lev.EQ.0.AND.iopt_sflux_precip.GT.1) THEN
   IF (iopt_meteo_precip.EQ.1) THEN
      mgvars(0)%rhs = mgvars(0)%rhs - rmaskatc*delt2d*&
                    & (evaporation-precipitation)/density_ref
   ELSEIF (iopt_meteo_precip.EQ.2) THEN
      mgvars(0)%rhs = mgvars(0)%rhs - rmaskatc*delt2d*evapminprec/density_ref
   ENDIF
ENDIF

!
!3.6 Boundary Conditions
!-----------------------
!

CALL open_boundary_conds_2d_impl(lev)

!
!4. Preconditioner
!-----------------
!

WHERE (mgvars(lev)%maskatc_int)
   mgvars(lev)%precon = 1.0/mgvars(lev)%subdiags(:,:,0)
END WHERE

!
!5. Set the flag for speedup in case of linear equations
!-------------------------------------------------------
!

mgvars(lev)%linearization_done = .TRUE.

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE(maskatu, maskatv,array2dc, array2du, array2dv)

CALL log_timer_out()


RETURN

END SUBROUTINE construct_mg_matrix

!========================================================================

RECURSIVE SUBROUTINE mg_cycles(numiter,mglevel)
!************************************************************************
!
! *mg_cycles* apply method using V- or W-cycling method
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Schemes.f90  V2.8
!
! Description - routine is called recursively
!
! Calling program - multi_grid_scheme
!
! External calls - construct_mg_matrix, mg_cycles, mg_smoother, prolongate,
!                  restrict, update_residual, update_solution
!
!************************************************************************
!
USE iopars
USE physpars
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: mglevel, numiter

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*numiter*  INTEGER Type of cycle
!            = 1 => On first call or V-cycle
!            = 2 => Subsequent (recursive) calls for W-cycle
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!

INTEGER :: iter, levC, levF, nocycle


procname(pglev+1) = 'mg_cycles'
CALL log_timer_in()

levF = mglevel; levC = levF + 1
nocycle = MERGE(1,2,iopt_mg_cycle.EQ.1)    

iter_110: DO iter=1,numiter 

!  ---update solution on "fine" grid
   CALL mg_smoothers(nopresweeps,levF)
   CALL construct_mg_matrix(0,levF)
!  ---residual on "fine" grid
   CALL update_residual(0,levF)
!  ---restrict from "fine" to "coarse" grid
   CALL restrict(levF,levC)
!  ---multi-grid cycling
   IF (levC.EQ.nomglevels-1) THEN
!     --at the coarsest level -> call for a direct solver
      CALL mg_smoothers(nosmoothsteps,levC)
   ELSE
!     --at intermediate level -> call for mg_cycles recursively
      pglev = pglev - 1
      CALL mg_cycles(nocycle,levC)
      pglev = pglev + 1
   ENDIF
!  ---prolongate from "coarse" to "fine" grid
   CALL prolongate(levC,levF)
!  ---update solution on "fine" grid
   CALL mg_smoothers(nopostsweeps,levF)

ENDDO iter_110

CALL log_timer_out()


RETURN

END SUBROUTINE mg_cycles

!========================================================================

SUBROUTINE mg_smoothers(nosweeps,mglevel)
!************************************************************************
!
! *mg_smoothers* update solution on the specified grid level and a given
!                number of sweeps
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Schemes.f90  V2.8
!
! Description -
!
! Calling program - mg_cycles, multigrid_scheme
!
! External calls - update_solution
!
!************************************************************************
!
USE iopars
USE physpars
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: mglevel, nosweeps

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*nosweeps* INTEGER Number of sweeps
!*mglevel*  INTEGER Multi-grid level for updating the solution (0 for the main
!                   grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: isweep

procname(pglev+1) = 'mg_smoothers'
CALL log_timer_in()

SELECT CASE(iopt_mg_smoother)

!
!1. Jacobi method
!----------------
!
   CASE (1)
      isweep_110: DO isweep=1,nosweeps
         CALL update_solution(0,isweep,nosweeps,mglevel,ur_smooth)
      ENDDO isweep_110

!
!2. Gauss-Seidel method
!----------------------
!

   CASE (2)
      isweep_210: DO isweep=1,nosweeps
!        ---red cells only
         CALL update_solution(1,isweep,nosweeps,mglevel,ur_smooth)
!        ---black cells only
         CALL update_solution(2,isweep,nosweeps,mglevel,ur_smooth)
      ENDDO isweep_210

END SELECT

CALL log_timer_out()


RETURN

END SUBROUTINE mg_smoothers

!========================================================================

SUBROUTINE multi_grid_scheme
!************************************************************************
!
! *multi_grid_scheme* multigrid solution procedure
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Schemes.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program -
!
! External calls - check_convergence, exchange_mod, mg_cycles,
!                  mg_pointer_weights, mg_smoother, mg_water_depths
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE multigrid
USE physpars
USE switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=12) :: citer
INTEGER :: lev, npcc


procname(pglev+1) = 'multi_grid_scheme'
CALL log_timer_in(npcc)

!
!1. Initialise
!-------------
!

mgvars%linearization_done = .FALSE.

lev_110: DO lev=0,nomglevels-1
   mgvars(lev)%correction = 0.0
   mgvars(lev)%precon = 0.0
   mgvars(lev)%residual = 0.0
   mgvars(lev)%rhs = 0.0
   mgvars(lev)%solvec = 0.0
   mgvars(lev)%subdiags = 0.0
ENDDO lev_110

!
!2. Update pointers and weight factors
!-------------------------------------
!

IF (iopt_fld.EQ.2.OR.iopt_thndam.EQ.1.OR.iopt_weibar.EQ.1) THEN
   CALL mg_pointers_weights
ENDIF

!
!3. Update water depths
!----------------------
!

CALL mg_water_depths
   
!
!4. Perform number of pre-defined iterations
!-------------------------------------------
!
!

mgiteration_410: DO mgiteration=1,nomgiterations
   IF (monlog.GT.1) THEN
      WRITE (citer,'(I12)') mgiteration; citer = ADJUSTL(citer)
      WRITE (iomon,'(A)') 'I: '//TRIM(citer) 
   ENDIF
   IF (nomglevels.EQ.1) THEN
!     ---one(fine)-level only
      CALL mg_smoothers(nosmoothsteps,0)
   ELSE
!     ---multi levels
      CALL mg_cycles(1,0)
   ENDIF
!  ---check for convergence 
   IF (residual_norm.LT.mg_tol) EXIT mgiteration_410
ENDDO mgiteration_410

mgiteration = MIN(mgiteration,nomgiterations)

!
!5. Return solution
!------------------
!
!---interior grid points
WHERE (nodeatc(0:ncloc+1,0:nrloc+1).GT.0)
   dzeta = mgvars(0)%solvec
END WHERE

CALL log_timer_out(npcc,itm_mg)


RETURN

END SUBROUTINE multi_grid_scheme

!========================================================================

SUBROUTINE prolongate(mglevelC,mglevelF)
!************************************************************************
!
! *prolongate* apply prolongation operator from coarse grid to next finer grid
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Schemes.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - 
!
! Module calls - error_alloc, exchange_mod, local_proc
!
!************************************************************************
!
USE gridpars
USE iopars
USE multigrid
USE structures
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: local_proc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
INTEGER, INTENT(IN) :: mglevelC, mglevelF

!
! Name        Type     Purpose
!-----------------------------------------------------------------------------
!*mglevelC*   INTEGER  Current coarse grid level
!*mglevelF*   INTEGER  Next fine grid level
!
!------------------------------------------------------------------------------
!*Local variables
!
INTEGER :: i, iglb, ii, i_C, i_F, j, jglb, jj, j_C, j_F, levC, levF, ncC, ncF, &
         & nrC, nrF
INTEGER, DIMENSION(4) :: nhexch
REAL :: denom
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, wcprod, weights


procname(pglev+1) = 'prolongate'
CALL log_timer_in()

!
!1. Exchange halo sections
!-------------------------
!

levC = mglevelC; levF = mglevelF
IF (iopt_MPI.EQ.1) THEN
   nhexch = 1
   CALL exchange_mod(mgvars(levC)%correction,(/0,0/),nhexch,0,corners=.TRUE.,&
                   & mglevel=levC)
ENDIF

!
!2. Initialise/allocate parameters and arrays
!--------------------------------------------
!

ncC = mgvars(levC)%ncloc; nrC = mgvars(levC)%nrloc
ncF = mgvars(levF)%ncloc; nrF = mgvars(levF)%nrloc

ALLOCATE (array2dc(0:ncF+1,0:nrF+1),STAT=errstat)
array2dc = 0.0
CALL error_alloc('wsum',2,(/ncF+2,nrF+2/),kndrtype)

IF (iopt_mg_prolong.EQ.2) THEN
   ALLOCATE (weights(0:ncC+1,0:nrC+1),STAT=errstat)
   CALL error_alloc('weights',2,(/ncC+2,nrC+2/),kndrtype)
   weights = mgvars(levC)%weights
   ALLOCATE (wcprod(0:ncC+1,0:nrC+1),STAT=errstat)
   CALL error_alloc('wcprod',2,(/ncC+2,nrC+2/),kndrtype)
   wcprod = mgvars(levC)%weights*mgvars(levC)%correction
ENDIF

!
!3. Interpolate correction
!-------------------------
!

SELECT CASE (iopt_mg_prolong)

!
!3.1 Pure injection
!------------------
!

   CASE (1)

      WHERE (mgvars(levC)%maskatc_int) 
         array2dc(1:ncF:2,1:nrF:2) = mgvars(levC)%correction(1:ncC,1:nrC) 
         array2dc(1:ncF:2,2:nrF+1:2) = mgvars(levC)%correction(1:ncC,1:nrC) 
         array2dc(2:ncF+1:2,1:nrF:2) = mgvars(levC)%correction(1:ncC,1:nrC) 
         array2dc(2:ncF+1:2,2:nrF+1:2) = mgvars(levC)%correction(1:ncC,1:nrC) 
      END WHERE

!
!3.2 Bi-linear interpolation
!---------------------------
!

   CASE (2) 

!
!3.2.1 Wet points
!----------------
!

      WHERE (mgvars(levC)%maskatc_int) 

         array2dc(1:ncF:2,1:nrF:2) = ( &
           &  9.0*wcprod(1:ncC,1:nrC)+3.0*wcprod(0:ncC-1,1:nrC) +& 
           &  3.0*wcprod(1:ncC,0:nrC-1)+wcprod(0:ncC-1,0:nrC-1)) &
           &/(9.0*weights(1:ncC,1:nrC)+3.0*weights(0:ncC-1,1:nrC) +&
           &  3.0*weights(1:nCC,0:nrC-1)+weights(0:ncC-1,0:nrC-1))
    
         array2dc(1:ncF:2,2:nrF+1:2) = ( &
           &  9.0*wcprod(1:ncC,1:nrC)+3.0*wcprod(0:ncC-1,1:nrC)+& 
           &  3.0*wcprod(1:ncC,2:nrC+1)+wcprod(0:ncC-1,2:nrC+1)) &
           &/(9.0*weights(1:ncC,1:nrC)+3.0*weights(0:ncC-1,1:nrC) +&
           &  3.0*weights(1:ncC,2:nrC+1)+weights(0:ncC-1,2:nrC+1))

         array2dc(2:ncF+1:2,1:nrF:2) = ( &
           &  9.0*wcprod(1:ncC,1:nrC)+3.0*wcprod(2:ncC+1,1:nrC) +& 
           &  3.0*wcprod(1:ncC,0:nrC-1)+wcprod(2:ncC+1,0:nrC-1)) &
           &/(9.0*weights(1:ncC,1:nrC)+3.0*weights(2:ncC+1,1:nrC) +&
           &  3.0*weights(1:ncC,0:nrC-1)+weights(2:ncC+1,0:nrC-1))

         array2dc(2:ncF+1:2,2:nrF+1:2) = ( &
           &  9.0*wcprod(1:ncC,1:nrC)+3.0*wcprod(2:ncC+1,1:nrC) +& 
           &  3.0*wcprod(1:ncC,2:nrC+1)+wcprod(2:ncC+1,2:nrC+1)) &
           &/(9.0*weights(1:ncC,1:nrC)+3.0*weights(2:ncC+1,1:nrC) +&
           &      weights(1:ncC,2:nrC+1)+weights(2:ncC+1,2:nrC+1))    

        END WHERE

!
!3.2.2 Case of small islands interpolation also at wet fine nodes
!----------------------------------------------------------------
!

        j_C_322: DO j_C=1,nrC
        i_C_322: DO i_C=1,ncC
           IF (.NOT.mgvars(levC)%maskatc_int(i_C,j_C)) THEN
              i_F = 2*i_C-1; j_F = 2*j_C-1
              IF (mgvars(levF)%nodeatc(i_F,j_F).GT.0) THEN
                 denom = 9.0*weights(i_C,j_C) + 3.0*weights(i_C-1,j_C) + &
                       & 3.0*weights(i_C,j_C-1) + weights(i_C-1,j_C-1)
                 IF (denom.GT.0.0) THEN
                    array2dc(i_F,j_F) = &
                         & (9.0*wcprod(i_C,j_C)+3.0*wcprod(i_C-1,j_C)+&
                         &  3.0*wcprod(i_C,j_C-1)+wcprod(i_C-1,j_C-1))/denom
                 ENDIF
              ENDIF
              i_F = 2*i_C-1; j_F = 2*j_C
              IF (mgvars(levF)%nodeatc(i_F,j_F).GT.0) THEN
                 denom = 9.0*weights(i_C,j_C) + 3.0*weights(i_C-1,j_C) + &
                       & 3.0*weights(i_C,j_C+1) + weights(i_C-1,j_C+1)
                 IF (denom.GT.0.0) THEN
                    array2dc(i_F,j_F) = &
                         & (9.0*wcprod(i_C,j_C)+3.0*wcprod(i_C-1,j_C)+&
                         &  3.0*wcprod(i_C,j_C+1)+wcprod(i_C-1,j_C+1))/denom
                 ENDIF
              ENDIF
              i_F = 2*i_C; j_F = 2*j_C-1
              IF (mgvars(levF)%nodeatc(i_F,j_F).GT.0) THEN
                 denom = 9.0*weights(i_C,j_C) + 3.0*weights(i_C+1,j_C) + &
                       & 3.0*weights(i_C,j_C-1) + weights(i_C+1,j_C-1)
                 IF (denom.GT.0.0) THEN
                    array2dc(i_F,j_F) = &
                         & (9.0*wcprod(i_C,j_C)+3.0*wcprod(i_C+1,j_C)+&
                         &  3.0*wcprod(i_C,j_C-1)+wcprod(i_C+1,j_C-1))/denom
                 ENDIF
              ENDIF
              i_F = 2*i_C; j_F = 2*j_C
              IF (mgvars(levF)%nodeatc(i_F,j_F).GT.0) THEN
                 denom = 9.0*weights(i_C,j_C) + 3.0*weights(i_C+1,j_C) + &
                       & 3.0*weights(i_C,j_C+1) + weights(i_C+1,j_C+1)
                 IF (denom.GT.0.0) THEN
                    array2dc(i_F,j_F) = &
                         & (9.0*wcprod(i_C,j_C)+3.0*wcprod(i_C+1,j_C)+&
                         &  3.0*wcprod(i_C,j_C+1)+wcprod(i_C+1,j_C+1))/denom
                 ENDIF
              ENDIF
           ENDIF
        ENDDO i_C_322
        ENDDO j_C_322

END SELECT

!
!4. Thin dams correction
!-----------------------
!

IF (levF.EQ.0.AND.iopt_thndam.EQ.1) THEN

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      nhexch = (/1,0,1,0/)
      CALL exchange_mod(array2dc,(/0,0/),nhexch,0,corners=.FALSE.)
   ENDIF

!  ---U-nodes
   ii_410: DO ii=1,numthinu
      iglb = ithinu(ii); jglb = jthinu(ii)
      i = iglb - nc1loc + 1; j = jglb - nr1loc + 1
      IF ((iglb.GT.nc1loc.AND.iglb.LE.(nc2loc+1)).AND.&
        & (jglb.GE.nr1loc.AND.jglb.LE.nr2loc)) THEN
         array2dc(i-1,j) = array2dc(i-2,j)
      ENDIF
      IF (local_proc(i,j)) array2dc(i,j) = array2dc(i+1,j)
   ENDDO ii_410

!  ---V-nodes
   jj_420: DO jj=1,numthinv
      iglb = ithinv(jj); jglb = jthinv(jj)
      i = iglb - nc1loc + 1; j = jglb - nr1loc + 1
      IF ((iglb.GE.nc1loc.AND.iglb.LE.nc2loc).AND.&
        & (jglb.GT.nr1loc.AND.jglb.LE.(nr2loc+1))) THEN
         array2dc(i,j-1) = array2dc(i,j-2)
      ENDIF
      IF (local_proc(i,j)) array2dc(i,j) = array2dc(i,j+1)
   ENDDO jj_420

ENDIF

!
!5. Add interpolated correction
!------------------------------
!

WHERE (mgvars(levF)%maskatc_int) 
   mgvars(levF)%correction(1:ncF,1:nrF) = mgvars(levF)%correction(1:ncF,1:nrF) &
                                      & + array2dc(1:ncF,1:nrF)
   mgvars(levF)%solvec(1:ncF,1:nrF) = mgvars(levF)%solvec(1:ncF,1:nrF) &
                                      & + array2dc(1:ncF,1:nrF)
END WHERE

!
!6. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = 1
   CALL exchange_mod(mgvars(levF)%solvec,(/0,0/),nhexch,0,corners=.TRUE.,&
                   & mglevel=levF)
ENDIF

!
!7. Deallocate arrays
!--------------------
!

DEALLOCATE(array2dc)
IF (iopt_mg_prolong.EQ.2) DEALLOCATE (weights,wcprod)

CALL log_timer_out()


RETURN

END SUBROUTINE prolongate

!========================================================================

SUBROUTINE restrict(mglevelF,mglevelC)
!************************************************************************
!
! *restrict* apply restriction operator to next coarser grid level
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Schemes.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program -
!
! External calls - construct_mg_matrix, update_residual

! Module calls - error_alloc, exchange_mod 
!
!************************************************************************
!
USE iopars
USE multigrid
USE switches
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: mglevelC, mglevelF

!
! Name      Type     Purpose
!-----------------------------------------------------------------------------
!*mglevelC* INTEGER  Level of next coarse grid
!*mglevelF* INTEGER  Level of current fine grid
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: levC, levF, ncC, ncF, nrC, nrF
INTEGER, DIMENSION(4) :: nhexch
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: wsum


procname(pglev+1) = 'restrict'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!

levC = mglevelC; levF = mglevelF
ncC = mgvars(levC)%ncloc; nrC = mgvars(levC)%nrloc
ncF = mgvars(levF)%ncloc; nrF = mgvars(levF)%nrloc

ALLOCATE (wsum(ncC,nrC),STAT=errstat)
CALL error_alloc('wsum',2,(/ncC,nrC/),kndrtype)

!
!1. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = (/0,1,0,1/)
   CALL exchange_mod(mgvars(levF)%residual,(/1,1/),nhexch,0,corners=.TRUE.,&
                   & mglevel=levF)
ENDIF

!
!2. Apply restrictor operator
!----------------------------
!

WHERE (mgvars(levC)%maskatc_int)

   wsum = mgvars(levF)%weights(1:ncF:2,1:nrF:2) + &
        & mgvars(levF)%weights(1:ncF:2,2:nrF+1:2) + &
        & mgvars(levF)%weights(2:ncF+1:2,1:nrF:2) + &
        & mgvars(levF)%weights(2:ncF+1:2,2:nrF+1:2)

   mgvars(levC)%residual(1:ncC,1:nrC) = &
     & mgvars(levF)%weights(1:ncF:2,1:nrF:2)*&
     & mgvars(levF)%residual(1:ncF:2,1:nrF:2) +&
     & mgvars(levF)%weights(1:ncF:2,2:nrF+1:2)*&
     & mgvars(levF)%residual(1:ncF:2,2:nrF+1:2) +&
     & mgvars(levF)%weights(2:ncF+1:2,1:nrF:2)*&
     & mgvars(levF)%residual(2:ncF+1:2,1:nrF:2) +&
     & mgvars(levF)%weights(2:ncF+1:2,2:nrF+1:2)*&
     & mgvars(levF)%residual(2:ncF+1:2,2:nrF+1:2)

   mgvars(levC)%solvec(1:ncC,1:nrC) = (&
     & mgvars(levF)%weights(1:ncF:2,1:nrF:2)*&
     & mgvars(levF)%solvec(1:ncF:2,1:nrF:2) +&
     & mgvars(levF)%weights(1:ncF:2,2:nrF+1:2)*&
     & mgvars(levF)%solvec(1:ncF:2,2:nrF+1:2) +&
     & mgvars(levF)%weights(2:ncF+1:2,1:nrF:2)*&
     & mgvars(levF)%solvec(2:ncF+1:2,1:nrF:2) +&
     & mgvars(levF)%weights(2:ncF+1:2,2:nrF+1:2)*&
     & mgvars(levF)%solvec(2:ncF+1:2,2:nrF+1:2))/wsum
   
END WHERE

!
!3. Deallocate
!-------------
!

DEALLOCATE (wsum)

!
!4. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = 1
   CALL exchange_mod(mgvars(levC)%solvec,(/0,0/),nhexch,0,corners=.FALSE.,&
                   & mglevel=levC)
ENDIF

!
!5. (Re)set correction
!---------------------
!

mgvars(levC)%correction = 0.0

!
!6. Set the right hand side of the coarse grid correction equation
!-----------------------------------------------------------------
!
!---update matrix entries
CALL construct_mg_matrix(0,levC)
!---store restricted defect in rhs
WHERE (mgvars(levC)%maskatc_int)
   mgvars(levC)%rhs = -mgvars(levC)%residual(1:ncC,1:nrC)
END WHERE
CALL update_residual(0,levC)
WHERE (mgvars(levC)%maskatc_int)
   mgvars(levC)%rhs = -mgvars(levC)%residual(1:ncC,1:nrC)
END WHERE

CALL log_timer_out()


RETURN

END SUBROUTINE restrict

!========================================================================

SUBROUTINE update_residual(itsweep,mglevel)
!************************************************************************
!
! *update_residual* update residual
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Schemes.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program -
!
! Module calls - error_alloc
!
!************************************************************************
!
USE iopars
USE multigrid
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: itsweep, mglevel

!
! Name      Type    Purpose
!-----------------------------------------------------------------------------
!*itsweep*  INTEGER Defines which matrix elements in the matrix multiplication
!            = 0 => all elements (Jacobi method)
!            = 1 => "red" points only (Gauus-Seidel) method
!            = 2 => "black" points only (Gauus-Seidel) method
!*mglevel*  INTEGER Multi-grid level (0 for the main grid)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: lev, nx, ny
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc


procname(pglev+1) = 'update_residual'
CALL log_timer_in()

lev = mglevel
nx = mgvars(lev)%ncloc; ny = mgvars(lev)%nrloc

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatc(nx,ny),STAT=errstat)
CALL error_alloc('maskatc',2,(/nx,ny/),kndlog)

!
!2. Mask arrays
!--------------
!

SELECT CASE (itsweep)
   CASE (0)
      maskatc = mgvars(lev)%maskatc_int
   CASE (1)
      maskatc = mgvars(lev)%maskatc_int.AND.mgvars(lev)%rednode
   CASE (2)
      maskatc = mgvars(lev)%maskatc_int.AND.mgvars(lev)%blacknode
END SELECT

!
!3. Calculate residual
!---------------------
!

WHERE (maskatc)
   mgvars(lev)%residual(1:nx,1:ny) = mgvars(lev)%rhs - &
    & mgvars(lev)%subdiags(:,:,0)*mgvars(lev)%solvec(1:nx,1:ny) - &
    & mgvars(lev)%subdiags(:,:,1)*mgvars(lev)%solvec(2:nx+1,1:ny) - &
    & mgvars(lev)%subdiags(:,:,-2)*mgvars(lev)%solvec(0:nx-1,1:ny) - &
    & mgvars(lev)%subdiags(:,:,2)*mgvars(lev)%solvec(1:nx,2:ny+1) - &
    & mgvars(lev)%subdiags(:,:,-1)*mgvars(lev)%solvec(1:nx,0:ny-1)   
END WHERE

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE(maskatc)

CALL log_timer_out()


RETURN

END SUBROUTINE update_residual

!========================================================================

SUBROUTINE update_solution(itsweep,iteration,maxiterations,mglevel,damping)
!************************************************************************
!
! *update_solution* update solution from previous iteration
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)MultiGrid_Schemes.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program -
!
! External calls - construct_mg_matrix, update_residual
!
! Module calls - error_alloc, exchange_mod
!
!************************************************************************
!
USE iopars
USE multigrid
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iteration, itsweep, maxiterations, mglevel
REAL, INTENT(IN) :: damping

!
! Name           Type    Purpose
!-----------------------------------------------------------------------------
!*itsweep*       INTEGER Defines which matrix elements in the matrix
!                        multiplication are activated
!                 = 0 => all elements (Jacobi method)
!                 = 1 => "red" points only (Gauss-Seidel) method
!                 = 2 => "black" points only (Gauss-Seidel) method
!*iteration*     INTEGER iteration index
!*maxiterations* INTEGER Number of (smoothing) iterations
!*mglevel*       INTEGER Multi-grid level (0 for the main grid)
!*damping*       REAL    Underrelaxation factor
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: lev, nx, ny
INTEGER, DIMENSION(4) :: nhexch
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc


procname(pglev+1) = 'update_solution'
CALL log_timer_in()

lev = mglevel
nx = mgvars(lev)%ncloc; ny = mgvars(lev)%nrloc

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array2dc(nx,ny),STAT=errstat)
CALL error_alloc('array2dc',2,(/nx,ny/),kndrlong)
array2dc = 0.0
ALLOCATE (maskatc(nx,ny),STAT=errstat)
CALL error_alloc('maskatc',2,(/nx,ny/),kndlog)

!
!2. Mask arrays
!--------------
!

SELECT CASE (itsweep)
   CASE(0)
      maskatc = mgvars(lev)%maskatc_int
   CASE(1)
      maskatc = mgvars(lev)%maskatc_int.AND.mgvars(lev)%rednode
   CASE(2)
      maskatc = mgvars(lev)%maskatc_int.AND.mgvars(lev)%blacknode
END SELECT
     
!
!3. Build the matrix components
!------------------------------
!
!---matrix elements 
CALL construct_mg_matrix(itsweep,lev)
!---residual
CALL update_residual(itsweep,lev)

!
!4. Solution and correction vector at next iteration
!---------------------------------------------------
!

WHERE (maskatc)
   array2dc = damping*mgvars(lev)%precon*&
            & mgvars(lev)%residual(1:nx,1:ny)
   mgvars(lev)%correction(1:nx,1:ny) = &
 & mgvars(lev)%correction(1:nx,1:ny) + array2dc
   mgvars(lev)%solvec(1:nx,1:ny) = &
 & mgvars(lev)%solvec(1:nx,1:ny) + array2dc
END WHERE

!
!5. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = 1
   IF (itsweep.NE.1.AND.iteration.EQ.maxiterations) THEN
      CALL exchange_mod(mgvars(lev)%solvec,(/0,0/),nhexch,0,corners=.TRUE.,&
                      & mglevel=lev)
   ELSE
      CALL exchange_mod(mgvars(lev)%solvec,(/0,0/),nhexch,0,corners=.FALSE.,&
                      & mglevel=lev)
   ENDIF
ENDIF

!
!6. Write monitoring info
!------------------------
!

IF (itsweep.NE.1) THEN
   CALL check_convergence(iteration,maxiterations,lev)
ENDIF

!
!7. Deallocate arrays
!--------------------
!

DEALLOCATE (array2dc,maskatc)

CALL log_timer_out()


RETURN

END SUBROUTINE update_solution
