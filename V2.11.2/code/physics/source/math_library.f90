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

MODULE math_library
!************************************************************************
!
! *math_library* Library of diverse mathematical routines
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description -
!
! Reference - Press W.H., Flannery B.P., Teukolsky S.A. and Vetterling W.T.,
!             1992. Numerical Recipes. The art of scientific computing.
!             Cambridge University Press, Cambridge, 702 pp.
!
! Generic routines - complex_polar
!
! Routines - brent_root, complex_sqrt_arr, complex_sqrt_var, gauss_quad,
!            polygon_check, polygon_mask, poly_all_roots, poly_div, poly_root,
!            quadratic_root_arr, quadratic_root_var, vector_mag_arr_atc,
!            vector_mag_arr_atu, vector_mag_arr_atv, vector_mag_var_atc,
!            vector_mag_var_atu, vector_mag_var_atv
!
!************************************************************************
!
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

INTERFACE complex_polar
   MODULE PROCEDURE complex_polar_0d, complex_polar_1d, complex_polar_2d, &
                  & complex_polar_3d, complex_polar_4d
END INTERFACE

CONTAINS

!========================================================================

FUNCTION brent_root(func,varname,ipars,rpars,noint,noreal,x1,x2,tol,fval,ierr)
!************************************************************************
!
! *brent_root* Find the root of the function "func" using Brent's method
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.0
!
! Description - inital location of root is assumed to be within the interval
!               [x1,x2]
!
! Reference - Numerical Recipes (routine "zbrent")
!
! Module calls - error_proc
!
!************************************************************************
!
USE error_routines, ONLY: error_proc

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(OUT) :: ierr
INTEGER, INTENT(IN) :: noint, noreal
REAL, INTENT(IN)  :: tol
REAL, INTENT(IN) :: x1, x2
REAL, INTENT(OUT) :: fval
INTEGER, INTENT(IN), DIMENSION(noint) :: ipars
REAL, INTENT(IN), DIMENSION(noreal) :: rpars

INTERFACE
   FUNCTION func(x,ipars,rpars,noint,noreal)
     INTEGER, INTENT(IN) :: noint, noreal
     INTEGER, INTENT(IN), DIMENSION(noint) :: ipars
     REAL, INTENT(IN), DIMENSION(noreal) :: rpars
     REAL,INTENT(IN) :: x
     REAL :: func
   END FUNCTION func
END INTERFACE

REAL :: brent_root

!
! Name       Type  Purpose
!------------------------------------------------------------------------------
!*func*      PROC  Input function
!*varname*   CHAR  FORTRAN name of function
!*ipars*     INTEGER Integer parameters in function call
!*rpars*     REAL    Real parameters in function call
!*noint*     INTEGER Number of integer parameters in function call
!*noreal*    INTEGER Number of real parameters in function call
!*x1*        REAL    Left value of bounding interval
!*x2*        REAL    Right value of bounding interval
!*tol*       REAL    Accuracy of root
!*fval*      REAL    Value of "func" at root (for consistency checking by
!                    calling program)
!*ierr*      INTEGER Output error code
!                  = 0 => no error; root has been found
!                  = 1 => x1 is equal to or larger than x2 
!                  = 2 => unable to to adjust interval [x1,x2] after maxcount
!                         iterations
!                  = 3 => unable to locate root after maxiter iterations
!*brent_root* REAL Root of function 'func'
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: citer
INTEGER :: maxcount = 20, maxit = 100
INTEGER :: icount, iter
REAL :: factor = 10.0
REAL :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm


procname(pglev+1) = 'brent_root'
CALL log_timer_in(logname=TRIM(procname(pglev+1))//': '//TRIM(varname))

ierr = 0

!
!1. Check and adjust bounding interval (if necessary)
!----------------------------------------------------
!

IF (x1.GE.x2) THEN
   ierr = 1
   GOTO 1001
ENDIF

a = x1; b = x2
fa = func(a,ipars,rpars,noint,noreal)
fb = func(b,ipars,rpars,noint,noreal)
icount = 0
DO WHILE ((fa*fb.GT.0).AND.(icount.LT.maxcount))
   icount = icount + 1
   IF (ABS(fa).LT.ABS(fb)) THEN
      a = a + factor*(a-b)
      fa = func(a,ipars,rpars,noint,noreal)
   ELSE
      b = b + factor*(b-a)
      fb = func(b,ipars,rpars,noint,noreal)
   ENDIF
ENDDO

IF (icount.EQ.maxcount) THEN
   ierr = 2
   GOTO 1001
ENDIF

c = b; fc = fb

!
!2. Locate zero by iteration
!---------------------------
!

iter_200: DO iter=1,maxit
   IF ((fb.GT.0.0.AND.fc.GT.0.0).OR.(fb.LT.0.0.AND.fc.LT.0.0)) THEN
      c = a; fc = fa
      d = b-a; e = d
   ENDIF
   IF (ABS(fc).LT.ABS(fb)) THEN
      a = b; b = c; c = a
      fa = fb; fb = fc; fc = fa
   ENDIF
   tol1 = 2.0*EPSILON(x1)*ABS(b) + 0.5*tol
   xm = 0.5*(c-b)
   IF ((ABS(xm).LE.tol1).OR.(fb.EQ.0.0)) THEN
      brent_root = b; fval = fb
      GOTO 1002
   ENDIF
   IF ((ABS(e).GE.tol1).AND.(ABS(fa).GT.ABS(fb))) THEN
      s = fb/fa
      IF (a.EQ.c) THEN
         p = 2.0*xm*s
         q = 1.0-s
      ELSE
         q = fa/fc; r = fb/fc
         p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
         q = (q-1.0)*(r-1.0)*(s-1.0)
      ENDIF
      IF (p.GT.0.0) q = -q
      p = ABS(p)
      IF (2.0*p.LT.MIN(3.0*xm*q-ABS(tol1*q),ABS(e*q))) THEN
         e = d; d = p/q
      ELSE
         d = xm; e = d
      ENDIF
   ELSE
      d = xm; e = d
   ENDIF
   a = b; fa = fb
   b = b + MERGE(d,SIGN(tol1,xm),ABS(d).GT.tol1)
   fb = func(b,ipars,rpars,noint,noreal)
ENDDO iter_200
      
ierr = 3
GOTO 1001

1002 CALL log_timer_out()


RETURN

1001 CONTINUE
IF (errchk) THEN
   WRITE (citer,'(I12)') iter; citer = ADJUSTL(citer)
   WRITE (ioerr,'(A)') 'Unable to find root of function '//TRIM(varname)//&
                     & ' after '//TRIM(citer)//' iterations'
ENDIF
CALL error_proc(ierr,abort=.TRUE.)

RETURN

END FUNCTION brent_root

!========================================================================

SUBROUTINE complex_polar_0d(xreal,ximag,xamp,xpha)
!************************************************************************
!
! *complex_polar_0d* Returns the amplitude and phase of a complex number
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.11.2
!
! Description - Phase is in radians and only returned if xpha is present
!
!************************************************************************
!
USE syspars

!
!*  Arguments
!
REAL, INTENT(IN) :: xreal, ximag
REAL, INTENT(OUT), OPTIONAL :: xamp, xpha

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*xreal*    REAL    Real part of complex number
!*ximag*    REAL    Imaginary part of complex number
!*xamp*     REAL    Amplitude of complex number
!*xpha*     REAL    Phase of complex number                              [rad]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL, PARAMETER :: epslim = 1.0E-20
REAL :: absi, absr


absr = ABS(xreal); absi = ABS(ximag)

!---amplitude
IF (PRESENT(xamp)) THEN
   IF (absr.LT.epslim.OR.absi.LT.epslim) THEN
      xamp = MAX(absr,absi)
   ELSEIF (absr.GT.absi) THEN
      xamp = absr*SQRT(1.0+(ximag/xreal)**2)
   ELSE
      xamp = absi*SQRT(1.0+(xreal/ximag)**2)
   ENDIF
ENDIF

!---phase
IF (PRESENT(xpha)) THEN
   IF (absr.GT.epslim.OR.absi.GT.epslim) THEN
      xpha = ATAN2(ximag,xreal)
      IF (xpha.LT.0.0) xpha = xpha + twopi
   ELSE
      xpha = MERGE(0.0,pi,absr.GE.absi)
   ENDIF
ENDIF


RETURN

END SUBROUTINE complex_polar_0d

!========================================================================

SUBROUTINE complex_polar_1d(xreal,ximag,xamp,xpha)
!************************************************************************
!
! *complex_polar_1d* Returns the amplitude and phase of a complex vector
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.11.2
!
! Description - Phases is in radians and only returned if xpha is present
!
!************************************************************************
!
USE syspars

!
!*  Arguments
!
REAL, INTENT(IN), DIMENSION(:) :: xreal, ximag
REAL, INTENT(OUT), OPTIONAL, DIMENSION(:) :: xamp, xpha

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*xreal*    REAL    Real part of complex array
!*ximag*    REAL    Imaginary part of complex array
!*xamp*     REAL    Amplitude of complex array
!*xpha*     REAL    Phase of complex array                               [rad]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: l
INTEGER, DIMENSION(1) :: ndims
REAL, PARAMETER :: epslim = 1.0E-20
REAL :: absi, absr


ndims = SIZE(xreal)

l_100: DO l=1,ndims(1)
   absr = ABS(xreal(l)); absi = ABS(ximag(l))
   
!
!1. Amplitudes
!-------------
!

   IF (PRESENT(xamp)) THEN
      IF (ABS(xreal(l)).LT.epslim.OR.ABS(ximag(l)).LT.epslim) THEN
         xamp(l) = MAX(ABS(xreal(l)),ABS(ximag(l)))
      ELSEIF (absr.GT.absi) THEN
         xamp(l) = absr*SQRT(1.0+(ximag(l)/xreal(l))**2)
      ELSE
         xamp(l) = absi*SQRT(1.0+(xreal(l)/ximag(l))**2)
      ENDIF
   ENDIF

!
!2. Phases
!---------
!

   IF (PRESENT(xpha)) THEN
      IF (absr.GT.epslim.OR.absi.GT.epslim) THEN
         xpha(l) = ATAN2(ximag(l),xreal(l))
         IF (xpha(l).LT.0.0) xpha(l) = xpha(l) + twopi
      ELSE
         xpha(l) = MERGE(0.0,pi,absr.GE.absi)
      ENDIF
   ENDIF

ENDDO l_100


RETURN

END SUBROUTINE complex_polar_1d

!========================================================================

SUBROUTINE complex_polar_2d(xreal,ximag,xamp,xpha,maskvals,outflag)
!************************************************************************
!
! *complex_polar_2d* Returns the amplitude and phase of a complex 2-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.11.2
!
! Description - Phases is in radians and only returned if xpha is present
!
!************************************************************************
!
USE syspars
USE error_routines, ONLY: error_alloc

!
!*  Arguments
!
REAL, INTENT(IN), OPTIONAL :: outflag
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: maskvals
REAL, INTENT(IN), DIMENSION(:,:) :: xreal, ximag
REAL, INTENT(OUT), OPTIONAL, DIMENSION(:,:) :: xamp, xpha

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*xreal*    REAL    Real part of complex array
!*ximag*    REAL    Imaginary part of complex array
!*xamp*     REAL    Amplitude of complex array
!*xpha*     REAL    Phase of complex array                               [rad]
!*maskvals* LOGICAL Horizontal array of masked (dry) points
!*outflag*  REAL    Flag for masked points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask
INTEGER :: i, j
REAL, PARAMETER :: epslim = 1.0E-20
REAL :: absi, absr, xflag
INTEGER, DIMENSION(2) :: ndims


!
!1. Initialise
!-------------
!
!---allocate
ndims = SHAPE(xreal)
ALLOCATE(mask(ndims(1),ndims(2)),STAT=errstat)
CALL error_alloc('mask',2,(/ndims(1),ndims(2)/),kndlog)

!---mask array
IF (PRESENT(maskvals)) THEN
   mask = maskvals
ELSE
   mask = .TRUE.
ENDIF

!
!2. Amplitudes and phases
!------------------------
!

IF (PRESENT(outflag)) THEN
   xflag = outflag
ELSE
   xflag = float_fill
ENDIF

j_210: DO j=1,ndims(2)
i_210: DO i=1,ndims(1)
   IF (mask(i,j)) THEN
      absr = ABS(xreal(i,j)); absi = ABS(ximag(i,j))
      IF (PRESENT(xamp)) THEN
         IF (absr.LT.epslim.OR.absi.LT.epslim) THEN
            xamp(i,j) = MAX(absr,absi)
         ELSEIF (absr.GT.absi) THEN
            xamp(i,j) = absr*SQRT(1.0+(ximag(i,j)/xreal(i,j))**2)
         ELSE
            xamp(i,j) = absi*SQRT(1.0+(xreal(i,j)/ximag(i,j))**2)
         ENDIF
      ENDIF
      IF (PRESENT(xpha)) THEN
         IF (absr.GT.epslim.OR.absi.GT.epslim) THEN
            xpha(i,j) = ATAN2(ximag(i,j),xreal(i,j))
            IF (xpha(i,j).LT.0.0) xpha(i,j) = xpha(i,j) + twopi
         ELSE
            xpha(i,j) = MERGE(0.0,pi,absr.GE.absi)
         ENDIF
      ENDIF
   ELSE
      IF (PRESENT(xamp)) xamp(i,j) = xflag  
      IF (PRESENT(xpha)) xpha(i,j) = xflag
   ENDIF
ENDDO i_210
ENDDO j_210

!
!3. Deallocate
!-------------
!

DEALLOCATE (mask)


RETURN

END SUBROUTINE complex_polar_2d

!========================================================================

SUBROUTINE complex_polar_3d(xreal,ximag,xamp,xpha,maskvals,outflag)
!************************************************************************
!
! *complex_polar_3d* Returns the amplitude and phase of a complex 3-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.11.2
!
! Description - Phases is in radians and only returned if xpha is present
!
!************************************************************************
!
USE syspars
USE error_routines, ONLY: error_alloc

!
!*  Arguments
!
REAL, INTENT(IN), OPTIONAL :: outflag
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: maskvals
REAL, INTENT(IN), DIMENSION(:,:,:) :: xreal, ximag
REAL, INTENT(OUT), OPTIONAL, DIMENSION(:,:,:) :: xamp, xpha

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*xreal*    REAL    Real part of complex array
!*ximag*    REAL    Imaginary part of complex array
!*xamp*     REAL    Amplitude of complex array
!*xpha*     REAL    Phase of complex array                               [rad]
!*maskvals* LOGICAL Horizontal array of masked (dry) points
!*outflag*  REAL    Flag for masked points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask
INTEGER :: i, j, k
REAL, PARAMETER :: epslim = 1.0E-20
REAL :: absi, absr, xflag
INTEGER, DIMENSION(3) :: ndims


!
!1. Initialise
!-------------
!
!---allocate
ndims = SHAPE(xreal)
ALLOCATE(mask(ndims(1),ndims(2)),STAT=errstat)
CALL error_alloc('mask',2,(/ndims(1),ndims(2)/),kndlog)

!---mask array
IF (PRESENT(maskvals)) THEN
   mask = maskvals
ELSE
   mask = .TRUE.
ENDIF

!
!2. Amplitudes and phases
!------------------------
!

IF (PRESENT(outflag)) THEN
   xflag = outflag
ELSE
   xflag = float_fill
ENDIF

k_210: DO k=1,ndims(3)
j_210: DO j=1,ndims(2)
i_210: DO i=1,ndims(1)
   IF (mask(i,j)) THEN
      absr = ABS(xreal(i,j,k)); absi = ABS(ximag(i,j,k))
      IF (PRESENT(xamp)) THEN
         IF (absr.LT.epslim.OR.absi.LT.epslim) THEN
            xamp(i,j,k) = MAX(absr,absi)
         ELSEIF (absr.GT.absi) THEN
            xamp(i,j,k) = absr*SQRT(1.0+(ximag(i,j,k)/xreal(i,j,k))**2)
         ELSE
            xamp(i,j,k) = absi*SQRT(1.0+(xreal(i,j,k)/ximag(i,j,k))**2)
         ENDIF
      ENDIF
      IF (PRESENT(xpha)) THEN
         IF (absr.GT.epslim.OR.absi.GT.epslim) THEN
            xpha(i,j,k) = ATAN2(ximag(i,j,k),xreal(i,j,k))
            IF (xpha(i,j,k).LT.0.0) xpha(i,j,k) = xpha(i,j,k) + twopi
         ELSE
            xpha(i,j,k) = MERGE(0.0,pi,absr.GE.absi)
         ENDIF
      ENDIF
   ELSE
      IF (PRESENT(xamp)) xamp(i,j,k) = xflag
      IF (PRESENT(xpha)) xpha(i,j,k) = xflag
   ENDIF
ENDDO i_210
ENDDO j_210
ENDDO k_210

!
!3. Deallocate
!-------------
!

DEALLOCATE (mask)


RETURN

END SUBROUTINE complex_polar_3d

!========================================================================

SUBROUTINE complex_polar_4d(xreal,ximag,xamp,xpha,maskvals,outflag)
!************************************************************************
!
! *complex_polar_4d* Returns the amplitude and phase of a complex 4-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.11.2
!
! Description - Phases is in radians and only returned if xpha is present
!
!************************************************************************
!
USE syspars
USE error_routines, ONLY: error_alloc

!
!*  Arguments
!
REAL, INTENT(IN), OPTIONAL :: outflag
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: maskvals
REAL, INTENT(IN), DIMENSION(:,:,:,:) :: xreal, ximag
REAL, INTENT(OUT), OPTIONAL, DIMENSION(:,:,:,:) :: xamp, xpha

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*xreal*    REAL    Real part of complex array
!*ximag*    REAL    Imaginary part of complex array
!*xamp*     REAL    Amplitude of complex array
!*xpha*     REAL    Phase of complex array                               [rad]
!*maskvals* LOGICAL Horizontal array of masked (dry) points
!*outflag*  REAL    Flag for masked points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask
INTEGER :: i, j, k, l
REAL, PARAMETER :: epslim = 1.0E-20
REAL :: absi, absr, xflag
INTEGER, DIMENSION(4) :: ndims


!
!1. Initialise
!-------------
!
!---allocate
ndims = SHAPE(xreal)
ALLOCATE(mask(ndims(1),ndims(2)),STAT=errstat)
CALL error_alloc('mask',2,(/ndims(1),ndims(2)/),kndlog)

!---mask array
IF (PRESENT(maskvals)) THEN
   mask = maskvals
ELSE
   mask = .TRUE.
ENDIF

!
!2. Ampltitude and phase
!-----------------------
!

IF (PRESENT(outflag)) THEN
   xflag = outflag
ELSE
   xflag = float_fill
ENDIF

l_210: DO l=1,ndims(4)
k_210: DO k=1,ndims(3)
j_210: DO j=1,ndims(2)
i_210: DO i=1,ndims(1)
   IF (mask(i,j)) THEN
      absr = ABS(xreal(i,j,k,l)); absi = ABS(ximag(i,j,k,l))
      IF (PRESENT(xamp)) THEN
         IF (absr.LT.epslim.OR.absi.LT.epslim) THEN
            xamp(i,j,k,l) = MAX(absr,absi)
         ELSEIF (absr.GT.absi) THEN
            xamp(i,j,k,l) = absr*SQRT(1.0+(ximag(i,j,k,l)/xreal(i,j,k,l))**2)
         ELSE
            xamp(i,j,k,l) = absi*SQRT(1.0+(xreal(i,j,k,l)/ximag(i,j,k,l))**2)
         ENDIF
      ENDIF
      IF (PRESENT(xpha)) THEN
         IF (absr.GT.epslim.OR.absi.GT.epslim) THEN
            xpha(i,j,k,l) = ATAN2(ximag(i,j,k,l),xreal(i,j,k,l))
            IF (xpha(i,j,k,l).LT.0.0) xpha(i,j,k,l) = xpha(i,j,k,l) + twopi
         ELSE
            xpha(i,j,k,l) =  MERGE(0.0,pi,absr.GE.absi)
         ENDIF
      ENDIF
   ELSE
      IF (PRESENT(xamp)) xamp(i,j,k,l) = xflag         
      IF (PRESENT(xpha)) xpha(i,j,k,l) = xflag
   ENDIF
ENDDO i_210
ENDDO j_210
ENDDO k_210
ENDDO l_210

!
!3. Deallocate
!-------------
!

DEALLOCATE (mask)


RETURN

END SUBROUTINE complex_polar_4d

!========================================================================

FUNCTION complex_sqrt_arr(z,ndims,mask,maskvals,outflag)
!************************************************************************
!
! *complex_sqrt_arr* Square root of a complex array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.7.1
!
! Description -
!
! Reference - Numerical Recipes
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: mask
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(4) :: ndims
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(ndims(1),ndims(2)) :: maskvals
COMPLEX, INTENT(IN), DIMENSION(ndims(1),ndims(2),ndims(3),ndims(4)) :: z
COMPLEX, DIMENSION(ndims(1),ndims(2),ndims(3),ndims(4)) :: complex_sqrt_arr

!
! Name              Type    Purpose
!------------------------------------------------------------------------------
!*z*                COMPLEX Complex input array
!*ndims*            INTEGER Shape of input array
!*mask*             LOGICAL Enables/disables flagging of dry array points
!*maskvals*         LOGICAL Horizontal array marking dy points
!*outflag*          REAL    Flag for dry array points
!*complex_sqrt_arr* COMPLEX Square root
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j, k, l
REAL :: a, b, w, xflag
COMPLEX :: s


procname(pglev+1) = 'complex_sqrt_arr'
CALL log_timer_in()

!
!1. No flagging
!--------------
!

IF (.NOT.mask) THEN

   l_110: DO l=1,ndims(4)
   k_110: DO k=1,ndims(3)
   j_110: DO j=1,ndims(2)
   i_110: DO i=1,ndims(1)
      a = REAL(z(i,j,k,l))
      b = AIMAG(z(i,j,k,l))
      IF ((a.EQ.0.0).AND.(b.EQ.0.0)) THEN
         complex_sqrt_arr(i,j,k,l) = 0.0
         CYCLE i_110
      ENDIF
      IF (ABS(a).GE.ABS(b)) THEN
         w = SQRT(ABS(a))*SQRT(0.5*(1.0+SQRT(1.0+(b/a)**2)))
      ELSE
         w = SQRT(ABS(b))*SQRT(0.5*(ABS(a/b)+SQRT(1.0+(a/b)**2)))
      ENDIF
      IF (a.GE.0.0) THEN
         s = CMPLX(w,0.5*b/w)
      ELSEIF (b.GE.0.0) THEN
         s = CMPLX(0.5*ABS(b)/w,w)
      ELSE
         s = CMPLX(0.5*ABS(b)/w,-w)
      ENDIF
      complex_sqrt_arr(i,j,k,l) = s
   ENDDO i_110
   ENDDO j_110
   ENDDO k_110
   ENDDO l_110

!
!2. With flagging
!----------------
!

ELSE

   IF (PRESENT(outflag)) THEN
      xflag = outflag
   ELSE
      xflag = float_fill
   ENDIF

   l_210: DO l=1,ndims(4)
   k_210: DO k=1,ndims(3)
   j_210: DO j=1,ndims(2)
   i_210: DO i=1,ndims(1)
      IF (maskvals(i,j)) THEN
         a = REAL(z(i,j,k,l))
         b = AIMAG(z(i,j,k,l))
         IF ((a.EQ.0.0).AND.(b.EQ.0.0)) THEN
            complex_sqrt_arr(i,j,k,l) = 0.0
            CYCLE i_210
         ENDIF
         IF (ABS(a).GE.ABS(b)) THEN
            w = SQRT(ABS(a))*SQRT(0.5*(1.0+SQRT(1.0+(b/a)**2)))
         ELSE
            w = SQRT(ABS(b))*SQRT(0.5*(ABS(a/b)+SQRT(1.0+(a/b)**2)))
         ENDIF
         IF (a.GE.0.0) THEN
            s = CMPLX(w,0.5*b/w)
         ELSEIF (b.GE.0.0) THEN
            s = CMPLX(0.5*ABS(b)/w,w)
         ELSE
            s = CMPLX(0.5*ABS(b)/w,-w)
         ENDIF
         complex_sqrt_arr(i,j,k,l) = s
      ENDIF
   ENDDO i_210
   ENDDO j_210
   ENDDO k_210
   ENDDO l_210

ENDIF

CALL log_timer_out()


RETURN

END FUNCTION complex_sqrt_arr

!========================================================================

FUNCTION complex_sqrt_var(z)
!************************************************************************
!
! *complex_sqrt_var* Square root of a complex number
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.0
!
! Description -
!
! Reference - Numerical Recipes
!
!************************************************************************
!
!*Arguments
!
COMPLEX, INTENT(IN) :: z
COMPLEX :: complex_sqrt_var

!
! Name              Type    Purpose
!------------------------------------------------------------------------------
!*z*                COMPLEX Complex input scalar
!*complex_sqrt_var* COMPLEX Square root
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL :: a, b, w
COMPLEX :: s


a = REAL(z)
b = AIMAG(z)

IF ((a.EQ.0.0).AND.(b.EQ.0.0)) THEN
   complex_sqrt_var = 0.0
   RETURN
ENDIF
IF (ABS(a).GE.ABS(b)) THEN
   w = SQRT(ABS(a))*SQRT(0.5*(1.0+SQRT(1.0+(b/a)**2)))
ELSE
   w = SQRT(ABS(b))*SQRT(0.5*(ABS(a/b)+SQRT(1.0+(a/b)**2)))
ENDIF

IF (a.GE.0.0) THEN
   s = CMPLX(w,0.5*b/w)
ELSEIF (b.GE.0.0) THEN
   s = CMPLX(0.5*ABS(b)/w,w)
ELSE
   s = CMPLX(0.5*ABS(b)/w,-w)
ENDIF

complex_sqrt_var = s


RETURN

END FUNCTION complex_sqrt_var

!=======================================================================

SUBROUTINE gauss_quad(nrpoints,weights,locations)
!***********************************************************************
!
! *qauss_quad* locations and weights for Gauss-Legendre quadrature
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)math_library.f90  V2.5
!
! Description - all polynomial coefficients are given as (/a_n, a_n-1, .., a_0/)
!
! Reference -
!
! Calling program - 
!
! External calls - 
!
! Module calls - poly_all_roots
!
!************************************************************************
!
USE iopars
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
INTEGER, INTENT(IN) :: nrpoints
REAL (KIND=kndrlong), INTENT(OUT), DIMENSION(nrpoints) :: locations, weights

!
! Name         Type       Purpose
!------------------------------------------------------------------------------
! *nrpoints*    INTEGER    The number of nodes used in the gaussian quadrature
! *weights*     REAL       The weight at each node
! *locations*   REAL       The locations of each node
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL (KIND=kndrlong), DIMENSION(nrpoints) :: legendre_diff_coef
REAL (KIND=kndrlong), DIMENSION(nrpoints+1) :: legendre_coef, v0, v1, v2 
COMPLEX (KIND=kndclong), DIMENSION(nrpoints) :: complex_locs
INTEGER :: ierr, p
REAL (KIND=kndrlong) :: n


procname(pglev+1) = 'gauss_quad'
CALL log_timer_in()

!
!1. Calculation of polynomial coefficients
!-----------------------------------------
!
!1.1 Legendre coefficients using Bonnet's recursion equation
!-----------------------------------------------------------
!

v0(1:nrpoints) = 0.0D0; v0(nrpoints+1)  = 1.0D0 
v1(1:nrpoints+1) = 0.0D0; v1(nrpoints)  = 1.0D0 

p_110: DO p = 2,nrpoints
    n = p-1.0D0;
    v2(1:nrpoints) = v1(2:nrpoints+1) 
    v2(nrpoints+1) = 0.0D0
    legendre_coef =  (2.0D0*n+1.0D0)/(n+1.0D0)*v2 
    legendre_coef = legendre_coef - n/(n+1.0D0)*v0
    v0 = v1
    v1 = legendre_coef
ENDDO p_110

!
!1.2 Coefficients of the Legendre polynomial's derivatives
!---------------------------------------------------------
!

legendre_diff_coef(1) = 0.0D0;
p_120: DO p = 2,nrpoints+1
    legendre_diff_coef(p-1) = (nrpoints-p+2.0D0)*legendre_coef(p-1)
 ENDDO p_120

!
!3. Quadrature locations
!-----------------------
!

CALL poly_all_roots(legendre_coef(nrpoints+1:1:-1),complex_locs,nrpoints,ierr)
locations = REAL(complex_locs,KIND=kndrlong)

IF (ierr.NE.0) THEN
   WRITE(iowarn,'(A)') 'Error when calculating roots of legendre equation.'
   WRITE(iowarn,'(A)') 'Integrals may not be calculated correctly.'
ENDIF

!
!3. Quadrature weight factors
!----------------------------
!

v0(1:nrpoints+1) = 0.0D0
v0(1:nrpoints) = v0(1:nrpoints) + legendre_diff_coef(nrpoints)
p_310: DO p = 2,nrpoints
   v0(1:nrpoints) = v0(1:nrpoints) + &
               & legendre_diff_coef(nrpoints+1-p)*(locations(1:nrpoints)**(p-1))
ENDDO p_310

weights = 2.0D0/((1.0D0-locations**2)*v0(1:nrpoints)**2)

CALL log_timer_out()


RETURN 

END SUBROUTINE gauss_quad

!========================================================================

SUBROUTINE polygon_check(xpol,ypol,npol,ierr)
!************************************************************************
!
! *polygon_test* check whether a polygon has double points or intersects itself
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)math_library.f90  V2.9
!
! Description -
!
! Reference - Numerical Recipes  par 21.4(routine "ispolysimple")
!
! Module calls - error_proc
!
!************************************************************************
!
USE paralpars
USE error_routines, ONLY: write_error_message

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: npol
REAL, DIMENSION(npol), INTENT(IN) :: xpol,ypol
INTEGER,  INTENT(OUT) :: ierr

!
! Name   Type    Purpose
!------------------------------------------------------------------------------
!*xpol*  REAL    X-Coordinates of the polygon to check
!*ypol*  REAL    Y-Coordinates of the polygon to check
!*npol*  INTEGER Number of points in the polygon
!*ierr*  INTEGER Error output code
!                = 0 => no error
!                = 1 => polygon contains double points
!                = 2 => polygon intersects itself
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: m, n
REAL, PARAMETER :: eps = 1.0E-15
REAL :: apx, apy, aqx, aqy, bpx, bpy, bqx, bqy, dist, t1, t2, t3, t4


procname(pglev+1) = 'polygon_check'
CALL log_timer_in()

ierr = 0

!
!1. Check whether one of the points occurs twice
!-----------------------------------------------
!
!---loop over each point in the polygon
n_110: DO n=1,npol-1
m_110: DO m =n+1,npol-1
   CALL complex_polar(xpol(n)-xpol(m),ypol(n)-ypol(m),dist)
   IF (dist.LE.eps) THEN
      ierr = 1
      RETURN
   ENDIF
ENDDO m_110
ENDDO n_110

!---check whether the first point of the polygon coincide with the last one
CALL complex_polar(xpol(npol)-xpol(1),ypol(npol)-ypol(1),dist)
IF (dist.GT.eps) THEN
   IF (master) THEN
      CALL write_error_message('First and last point of polygon are not equal')
   ENDIF
   ierr = 1
ENDIF

!
!2. Check all pairs of intersections
!-----------------------------------
!

n_210: DO n=1,npol-1
!  ---start and end point of the first section
   apx = xpol(n); apy = ypol(n)
   aqx = xpol(n+1); aqy = ypol(n+1)
   m_211: DO m= n+1,npol
!     --start and end point of the second section
      bpx = xpol(m); bpy = ypol(m)
      bqx = xpol(m+1); bqy = ypol(m+1)
!     --check if bp and bq are on different sides of ap-aq
      t1 = (aqx-apx)*(bpy-apy) - (bpx-apx)*(aqy-apy)
      t2 = (aqx-apx)*(bpy-apy) - (bpx-apx)*(aqy-apy)
!     --check if ap and aq are on different sides of bp-bq
      t3 = (bqx-bpx)*(apy-bpy) - (apx-bpx)*(bqy-bpy)
      t4 = (bqx-bpx)*(apy-bpy) - (apx-bpx)*(bqy-bpy)
!     ---found an intersection
      IF (t1*t2.GT.0.0.AND.t3*t4.GT.0.0) THEN
         ierr = 2
         RETURN
      ENDIF
   ENDDO m_211
ENDDO n_210

CALL log_timer_out()


RETURN

END SUBROUTINE polygon_check

!===============================================================================

SUBROUTINE polygon_mask(xpoly,ypoly,npoly,xcoord,ycoord,nxgrd,nygrd,mask,info)
!*******************************************************************************
!
! *polygon_mask* check which points on the 2-D input grid are inside a given
!                polygon
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)math_library.f90  V2.8
!
! Description - returns result in array 'mask'
!
! Reference - Numerical Recipes  par 21.4 (routine "polywind")
!
! Module calls - error_proc
!
!************************************************************************
!
USE iopars
USE syspars
USE error_routines, ONLY: error_proc

IMPLICIT NONE
!
!*Arguments
!
LOGICAL, INTENT(IN) :: info
INTEGER, INTENT(IN) :: npoly, nxgrd, nygrd
REAL, DIMENSION(npoly), INTENT(IN) :: xpoly, ypoly
REAL, DIMENSION(nxgrd,nygrd), INTENT(IN) :: xcoord ,ycoord
LOGICAL, DIMENSION(nxgrd,nygrd), INTENT(OUT) :: mask

!
! Name    Type    Purpose
!------------------------------------------------------------------------------
!*xpoly*   REAL   X-Coordinates the points defining the the polygon.
!                 The last point must be equal to the first one
!*ypoly*   REAL   Y-Coordinates of the points defining the polygon.
!                 The last point must be equal to the first one
!*npoly*  INTEGER Number of edge points defining the polygon
!*xcoord* REAL    X-Coordinates of the input grid
!*ycoord* REAL    Y-Coordinates of the input grid
!*nxgrd*  INTEGER Number of cells in the X-direction of the input grid
!*nygrd*  INTEGER Number of cells in the Y-direction of the input grid 
!*mask*   REAL    Returned grid mask (.TRUE. at grid points inside the polygon,
!                 .FALSE. otherwise
!*info*   LOGICAL Disables/enables writing of log info
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j, winding
REAL :: x, y
REAL, DIMENSION(npoly-1) :: signed_area


procname(pglev+1) = 'polygon_mask'
CALL log_timer_in(info=info)

mask = .FALSE.

j_110: DO i=1,nxgrd
i_110: DO j=1,nygrd
   x = xcoord(i,j); y = ycoord(i,j)
!  ---area of the triangle with the point
   signed_area = (xpoly(1:npoly-1)-x)*(ypoly(2:npoly)-y) - &
               & (ypoly(1:npoly-1)-y)*(xpoly(2:npoly)-x)
!  ---winding angle
   winding = COUNT((signed_area.GT.0).AND.(ypoly(2:npoly).GT.y).AND.&
                &  (ypoly(1:npoly-1).LE.y)) - &
           & COUNT((signed_area.LT.0).AND.(ypoly(2:npoly).LE.y).AND.&
           &       (ypoly(1:npoly-1).GT.y))
!  ---check whether the point is inside the polygon
   mask(i,j) =  winding.EQ.1.OR.winding.EQ.-1
ENDDO i_110
ENDDO j_110

CALL log_timer_out(info=info)


RETURN

END SUBROUTINE polygon_mask

!========================================================================

SUBROUTINE poly_all_roots(a,x,m,ierr)
!************************************************************************
!
! *poly_all_roots* Find the all the roots of a polynomial by deflation
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)math_library.f90  V2.5
!
! Description - coefficients "a" of the polynomial are given as a(1:m+1) where
!               m is the order of the polynomial 
!
! Reference - Numerical Recipes
!
! Module calls - error_proc, poly_div, poly_root
!
!************************************************************************
!
USE syspars
USE error_routines, ONLY: error_proc

!
!*Arguments
!
INTEGER, INTENT(IN) :: m
INTEGER, INTENT(OUT) :: ierr
COMPLEX (KIND=kndclong), INTENT(OUT), DIMENSION(m)  :: x
REAL (KIND=kndrlong), INTENT(IN), DIMENSION(m+1) :: a

!
! Name    Type    Purpose
!------------------------------------------------------------------------------
!*a*      REAL    Coefficients of polynomial (a(1) for constant term
!                 a(m+1) for highest power x**m)
!*x*      COMPLEX Roots of the equation
!*m*      INTEGER Degree of polynomial
!*ierr*   INTEGER Output error code (0 for successfull exit)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, npcc
REAL (KIND=kndrlong), DIMENSION(m+1)  :: aa, dummy
REAL (KIND=kndrlong), DIMENSION(2)  :: v2
REAL (KIND=kndrlong), DIMENSION(3)  :: v3
COMPLEX (KIND=kndclong)  :: x0
REAL (KIND=kndrlong), PARAMETER :: small = 1.0D-15


procname(pglev+1) = 'poly_all_roots'
CALL log_timer_in(npcc)

!
!1. Calling the root solver and deflating the polynome
!-----------------------------------------------------
!

aa = a
x0 = 0.0D0
k = 1
k_110: DO WHILE (k.LE.m)
!  ---solving
   CALL poly_root(aa,x0,m-k+1,ierr)
!  ---deflating
   IF (ABS(IMAG(x0)).LT.small) THEN
!     --real Root
      v2(1) = -REAL(x0,KIND=kndrlong)
      v2(2) = 1
      x0 = CMPLX(REAL(x0,KIND=kndrlong),0.0D0,KIND=kndrlong)
      x(k) = x0
      CALL poly_div(aa,m+1,v2,2,dummy)
      k = k+1
   ELSE
!     --complex root
      v3(1) = REAL(x0,KIND=kndrlong)**2 + IMAG(x0)**2
      v3(2) = -2.0D0*REAL(x0,KIND=kndrlong)
      v3(3) = 1.0D0
      x(k) = x0
      x(k+1) = CONJG(x0)
      CALL poly_div(aa,m+1,v3,3,dummy)
      k = k+2
   ENDIF
ENDDO k_110

!
!2. Polishing the polynome using the found roots
!-----------------------------------------------
!

k_210: DO k= 1,m
   CALL poly_root(a,x(k),m,ierr)
ENDDO k_210

CALL log_timer_out(npcc,itm_libs)


RETURN

END SUBROUTINE poly_all_roots

!========================================================================

SUBROUTINE poly_div(u,nu,v,nv,r)
!************************************************************************
!
! *poly_div* Divide two polynomials
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)math_library.f90  V2.5
!
! Description -
!
! Reference - Numerical Recipes (routine "poldiv")
!
! Module calls - error_proc
!
!************************************************************************
!
USE syspars
USE error_routines, ONLY: error_proc

!
!*Arguments
!
INTEGER, INTENT(IN) :: nu, nv
REAL (KIND=kndrlong), INTENT(OUT), DIMENSION(nu) :: r
REAL (KIND=kndrlong), INTENT(INOUT), DIMENSION(nu) :: u
REAL (KIND=kndrlong), INTENT(IN), DIMENSION(nv) :: v

!
! Name    Type    Purpose
!------------------------------------------------------------------------------
!*u*      REAL    Coefficients of nominator polynomial (u(1) for constant term
!                 u(nu) for highest power x**nu)
!*v*      REAL    Coefficients of denominator polynomial
!*x*      COMPLEX Initial guess on input; root on output
!*nu*     INTEGER Degree+1 of polynomial u 
!*nv*     INTEGER Degree+1 of polynomial v
!*r*      REAL    Coefficients of polynomial with remainder
!
!------------------------------------------------------------------------------
!
!*Local variables
!

INTEGER ::  j, k

procname(pglev+1) = 'poly_div'
CALL log_timer_in()

!
!1. Initialisation
!-----------------
!

r = u; u = 0.0D0

!
!2 Division
!----------
!

k_210: DO k = nu-nv,0,-1
   u(k+1) = r(nv+k)/v(nv)
   j_211: DO j = nv+k-1,k+1,-1
      r(j) = r(j) - u(k+1)*v(j-k)
   ENDDO j_211
ENDDO k_210

!
!3. Finalisation
!---------------
!

j_310: DO j = nv,nu
   r(j) = 0.0D0
ENDDO j_310

CALL log_timer_out()


RETURN

END SUBROUTINE poly_div

!========================================================================

SUBROUTINE poly_root(a,x,m,ierr)
!************************************************************************
!
! *poly_root* Find the root of polynomial with initial guess 'x' using
!             Laguerre's method
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.0
!
! Description - coefficients "a" of the polynomial are given as a(1:m+1) where
!               m is the order of the polynomial 
!
! Reference - Numerical Recipes (routine "laguer")
!
! Module calls - error_proc
!
!************************************************************************
!
USE syspars
USE error_routines, ONLY: error_proc

!
!*Arguments
!
INTEGER, INTENT(IN) :: m
INTEGER, INTENT(OUT) :: ierr
COMPLEX (KIND=kndclong), INTENT(INOUT) :: x
REAL (KIND=kndrlong), INTENT(IN), DIMENSION(m+1) :: a

!
! Name    Type    Purpose
!------------------------------------------------------------------------------
!*a*      REAL    Coefficients of polynomial (a(1) for constant term
!                 a(m+1) for highest power x**m)
!*x*      COMPLEX Initial guess on input; root on output
!*m*      INTEGER Degree of polynomial
!*ierr*   INTEGER Output error code (0 for successfull exit)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: citer
INTEGER, PARAMETER :: mr = 8, mt = 10, maxit = mt*mr
INTEGER :: iter, j, npcc
REAL (KIND=kndrlong) :: abm, abp, abx, err
COMPLEX (KIND=kndclong) :: b, cm, d, dx, f, g, gm, gp, g2, h, sq, x1
REAL (KIND=kndrlong), DIMENSION(mr) :: frac = (/0.5D0,0.25D0,0.75D0,0.13D0,&
                                            & 0.38D0,0.62D0,0.88D0,1.0D0/)

!
! Name    Type    Purpose
!------------------------------------------------------------------------------
!*maxit*  INTEGER Maximum number of allowed iterations
!*mt*     INTEGER Number of iterations to break any limit cycle
!*frac*   REAL    Fractional values used after each mt iterations
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'poly_root'
CALL log_timer_in(npcc)

ierr = 0

iter_100: DO iter=1, maxit
   d = CMPLX(0.0,KIND=kndrlong); f = CMPLX(0.0,KIND=kndrlong)
   b = a(m+1)
   err = ABS(b)
   abx = ABS(x)
   j_110: DO j=m,1,-1
      f = x*f + d
      d = x*d + b
      b = x*b + a(j)
      err = ABS(b) + abx*err
   ENDDO j_110
   err = EPSILON(1.0)*err
   IF (ABS(b).LE.err) GOTO 1000
   g = d/b
   g2 = g*g
   h = g2 - 2.0D0*f/b
   cm = CMPLX(m,KIND=kndrlong)
   sq = SQRT((cm-1.0D0)*(cm*h-g2))
   gp = g + sq
   gm = g - sq
   abp = ABS(gp)
   abm = ABS(gm)
   IF (abp.LT.abm) gp = gm
   IF (MAX(abp,abm).GT.0.0D0) THEN
      dx = m/gp
   ELSE
      dx = EXP(CMPLX(LOG(1.0D0+abx),REAL(iter),KIND=kndrlong))
   ENDIF
   x1 = x - dx
   IF (x.EQ.x1) GOTO 1000
   IF (MOD(iter,mt).NE.0) THEN
      x = x1
   ELSE
      x = x - dx*frac(iter/mt)
   ENDIF
ENDDO iter_100

ierr = 1
IF (errchk) THEN
   WRITE (citer,'(I12)') iter; citer = ADJUSTL(citer)
   WRITE (ioerr,'(A)') 'Unable to find root of polynomial after '&
                     & //TRIM(citer)//' iterations'
ENDIF
CALL error_proc(ierr,abort=.TRUE.)

1000 CALL log_timer_out(npcc,itm_libs)


RETURN

END SUBROUTINE poly_root

!========================================================================

SUBROUTINE quadratic_root_arr(coef,root1,root2,ndims,outflag,mask,maskvals)
!************************************************************************
!
! *quadratic_root_arr* Solve the quadratic equation with array coefficients
!                      coef(*,3)*x^2 + coef(*,2)*x + coef(*,1) = 0
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.0
!
! Description - 
!
! Reference - 
!
! Module calls - complex_sqrt_var
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: mask
REAL, INTENT(IN) :: outflag
INTEGER, INTENT(IN), DIMENSION(4) :: ndims
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(ndims(1),ndims(2)) :: maskvals
REAL, INTENT(IN), DIMENSION(ndims(1),ndims(2),ndims(3),ndims(4),3) :: coef
COMPLEX, INTENT(OUT), DIMENSION(ndims(1),ndims(2),ndims(3),ndims(4)) :: &
                              & root1, root2

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*coef*     REAL     Array coefficients of quadratic equation
!*root1*    COMPLEX  First root
!*root2*    COMPLEX Second root
!*ndims*    INTEGER  Shape of coeeficient arrays
!*outflag*  REAL     Flag for undefined output
!*mask*     LOGICAL  Enables/disables flagging of dry array points
!*maskvals* LOGICAL  Horizontal array marking dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j, k, l, npcc
REAL :: a, b, c, s
COMPLEX :: discrim, q


procname(pglev+1) = 'quadratic_root_arr'
CALL log_timer_in(npcc)

!
!1. No flagging
!--------------
!

IF (.NOT.mask) THEN

   l_110: DO l=1,ndims(4)
   k_110: DO k=1,ndims(3)
   j_110: DO j=1,ndims(2)
   i_110: DO i=1,ndims(1)

      a = coef(i,j,k,l,3)
      IF (a.EQ.0.0) THEN
         root1(i,j,k,l) = outflag
         root2(i,j,k,l) = outflag
         CYCLE i_110
      ENDIF
      b = coef(i,j,k,l,2)
      c = coef(i,j,k,l,1)
      discrim = CMPLX(b*b-4.0*c*a,0.0)
      discrim = complex_sqrt_var(discrim)
      s = SIGN(1.0,b)
      q = -0.5*(b+s*discrim)
      root1(i,j,k,l) = q/a
      root2(i,j,k,l) = c/q
      
   ENDDO i_110
   ENDDO j_110
   ENDDO k_110
   ENDDO l_110

!
!2. With flagging
!----------------
!

ELSE

   l_120: DO l=1,ndims(4)
   k_120: DO k=1,ndims(3)
   j_120: DO j=1,ndims(2)
   i_120: DO i=1,ndims(1)
      IF (maskvals(i,j)) THEN
         a = coef(i,j,k,l,3)
         IF (a.EQ.0.0) THEN
            root1(i,j,k,l) = outflag
            root2(i,j,k,l) = outflag
            CYCLE i_120
         ENDIF
         b = coef(i,j,k,l,2)
         c = coef(i,j,k,l,1)
         discrim = CMPLX(b*b-4.0*c*a,0.0)
         discrim = complex_sqrt_var(discrim)
         s = SIGN(1.0,b)
         q = -0.5*(b+s*discrim)
         root1(i,j,k,l) = q/a
         root2(i,j,k,l) = c/q
      ENDIF
      
   ENDDO i_120
   ENDDO j_120
   ENDDO k_120
   ENDDO l_120

ENDIF

CALL log_timer_out(npcc,itm_libs)


RETURN

END SUBROUTINE quadratic_root_arr

!========================================================================

SUBROUTINE quadratic_root_var(coef,root1,root2,outflag)
!************************************************************************
!
! *quadratic_root_var* Solve the quadratic equation with scalar coefficients
!                      coef(3)*x^2 + coef(2)*x + coef(1) = 0
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.0
!
! Description - 
!
! Reference - 
!
! Module calls - complex_sqrt_var
!
!************************************************************************
!
!*Arguments
!
REAL, INTENT(IN), DIMENSION(3) :: coef
REAL, INTENT(IN) :: outflag 
COMPLEX, INTENT(OUT) :: root1, root2

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*coef*    REAL     Coeficients of quadratic equation
!*root1*   COMPLEX  First root
!*root2*   COMPLEX  Second root
!*outflag* REAL     Flag for undefined output
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: npcc
REAL :: s
COMPLEX :: discrim, q


CALL log_timer_in(npcc)

IF (coef(3).EQ.0.0) THEN
   root1 = outflag
   root2 = outflag
ELSE
   discrim = CMPLX(coef(2)*coef(2)-4.0*coef(1)*coef(3),0.0)
   discrim = complex_sqrt_var(discrim)
   s = SIGN(1.0,coef(2))
   q = -0.5*(coef(2)+s*discrim)
   root1 = q/coef(3)
   root2 = coef(1)/q
ENDIF

CALL log_timer_out(npcc,itm_libs)


RETURN

END SUBROUTINE quadratic_root_var

!========================================================================

SUBROUTINE vector_mag_arr_atc(xcomp,ycomp,intsrce,intdest,nzdim,&
                            & nosize,ivarid,info,vecmag,vecpha,outflag)
!************************************************************************
!
! *vector_mag_arr_atc* Calculate the (horizontal) magnitude of a vector at the
!                      C-nodes of the model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.11.2
!
! Description -
!
! Reference - 
!
! Module calls - complex_polar, Uarr_at_C, Varr_at_C
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE array_interp, ONLY: Uarr_at_C, Varr_at_C
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
LOGICAL, INTENT(IN) :: info
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize, nzdim
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(ncloc+1,nrloc,nzdim,nosize) :: xcomp
REAL, INTENT(IN), DIMENSION(ncloc,nrloc+1,nzdim,nosize) :: ycomp
REAL, INTENT(OUT), OPTIONAL, DIMENSION(ncloc,nrloc,nzdim,nosize) :: vecmag, &
                                                                  & vecpha

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xcomp*     REAL    X-component of input vector at the U-nodes
!*ycomp*     REAL    Y-component of input vector at the V-nodes
!*vecmag*    REAL    Returned vector magnitude at the C-nodes
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
!*nzdim*     INTEGER Third dimension of input vector
!*nosize*    INTEGER Fourth dimension of input vector
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL :: flagval
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskvals
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: array4dc1, array4dc2


procname(pglev+1) = 'vector_mag_arr_atc'
CALL log_timer_in()

!
!1. Optional argument
!--------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Allocate arrays
!------------------
!

ALLOCATE (maskvals(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskvals',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (array4dc1(ncloc,nrloc,nzdim,nosize),STAT=errstat)
CALL error_alloc('array4dc1',4,(/ncloc,nrloc,nzdim,nosize/),kndrtype)
ALLOCATE (array4dc2(ncloc,nrloc,nzdim,nosize),STAT=errstat)
CALL error_alloc('array4dc2',4,(/ncloc,nrloc,nzdim,nosize/),kndrtype)

!
!3. Interpolate components at C-nodes
!-------------------------------------
!

CALL Uarr_at_C(xcomp,array4dc1,intsrce,intdest,(/1,1,1/),&
            & (/ncloc+1,nrloc,nzdim/),nosize,ivarid,info,flagval)
CALL Varr_at_C(ycomp,array4dc2,intsrce,intdest,(/1,1,1/),&
            & (/ncloc,nrloc+1,nzdim/),nosize,ivarid,info,flagval)

!
!4. Vector magnitude
!-------------------
!

maskvals = MERGE(maskatc_int,.TRUE.,intdest.EQ.1)
IF (.NOT.PRESENT(vecpha)) THEN
   CALL complex_polar(array4dc1,array4dc2,xamp=vecmag,maskvals=maskvals,&
                    & outflag=flagval)
ELSEIF(.NOT.PRESENT(vecmag)) THEN
   CALL complex_polar(array4dc1,array4dc2,xpha=vecpha,maskvals=maskvals,&
                    & outflag=flagval)
ELSE
   CALL complex_polar(array4dc1,array4dc2,xamp=vecmag,xpha=vecpha,&
                    & maskvals=maskvals,outflag=flagval)
ENDIF

!
!5. Deallocate
!-------------
!

DEALLOCATE (maskvals,array4dc1,array4dc2)

CALL log_timer_out()


RETURN

END SUBROUTINE vector_mag_arr_atc

!========================================================================

SUBROUTINE vector_mag_arr_atu(xcomp,ycomp,intsrce,intdest,nzdim,nosize,ivarid,&
                            & info,vecmag,vecpha,outflag,hregular)
!************************************************************************
!
! *vector_mag_arr_atu* Calculate the (horizontal) magnitude of a vector at the
!                      U-nodes of the model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.11.2
!
! Description -
!
! Reference - 
!
! Module calls - complex_polar, Varr_at_U
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE switches
USE array_interp, ONLY: Varr_at_U
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize, nzdim
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nzdim,nosize) :: xcomp
REAL, INTENT(IN), DIMENSION(0:ncloc,nrloc+1,nzdim,nosize) :: ycomp
REAL, INTENT(OUT), OPTIONAL, DIMENSION(ncloc,nrloc,nzdim,nosize) :: vecmag, &
                                                                  & vecpha

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xcomp*     REAL    X-component of input vector at the U-nodes
!*ycomp*     REAL    Y-component of input vector at the V-nodes
!*vecmag*    REAL    Returned vector magnitude at the U-nodes
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
!*nzdim*     INTEGER Third dimension of input vector
!*nosize*    INTEGER Fourth dimension of input vector
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: hreg
REAL :: flagval
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskvals
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: array4du


procname(pglev+1) = 'vector_mag_arr_atu'
CALL log_timer_in()

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
!2. Allocate arrays
!------------------
!

ALLOCATE (maskvals(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskvals',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (array4du(ncloc,nrloc,nzdim,nosize),STAT=errstat)
CALL error_alloc('array4du',4,(/ncloc,nrloc,nzdim,nosize/),kndrtype)

!
!3. Interpolate components at U-nodes
!-------------------------------------
!

CALL Varr_at_U(ycomp,array4du,intsrce,intdest,(/0,1,1/),&
            & (/ncloc,nrloc+1,nzdim/),nosize,ivarid,info,flagval,hreg)

!
!4. Mask array
!-------------
!

SELECT CASE (intdest)
   CASE (0); maskvals = .TRUE.
   CASE (1); maskvals = node2du(1:ncloc,1:nrloc).GT.0
   CASE (2); maskvals = node2du(1:ncloc,1:nrloc).EQ.2
   CASE (3); maskvals = node2du(1:ncloc,1:nrloc).GT.1
   CASE (4); maskvals = node2du(1:ncloc,1:nrloc).EQ.1.OR.&
                      & node2du(1:ncloc,1:nrloc).EQ.2
END SELECT

!
!5. Vector magnitude
!-------------------
!

IF (.NOT.PRESENT(vecpha)) THEN
   CALL complex_polar(xcomp,array4du,xamp=vecmag,maskvals=maskvals,&
                    & outflag=flagval)
ELSEIF(.NOT.PRESENT(vecmag)) THEN
   CALL complex_polar(xcomp,array4du,xpha=vecpha,maskvals=maskvals,&
                    & outflag=flagval)
ELSE
   CALL complex_polar(xcomp,array4du,xamp=vecmag,xpha=vecpha,maskvals=maskvals,&
                    & outflag=flagval)
ENDIF

!
!6. Deallocate
!-------------
!

DEALLOCATE (maskvals,array4du)

CALL log_timer_out()


RETURN

END SUBROUTINE vector_mag_arr_atu

!========================================================================

SUBROUTINE vector_mag_arr_atv(xcomp,ycomp,intsrce,intdest,nzdim,nosize,ivarid,&
                            & info,vecmag,vecpha,outflag,hregular)
!************************************************************************
!
! *vector_mag_arr_atv* Calculate the (horizontal) magnitude of a vector at the
!                      V-nodes of the model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.11.2
!
! Description -
!
! Reference - 
!
! Module calls - complex_polar, Varr_at_U
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE switches
USE array_interp, ONLY: Uarr_at_V
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
LOGICAL, INTENT(IN) :: info
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: intdest, intsrce, ivarid, nosize, nzdim
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(ncloc+1,0:nrloc,nzdim,nosize) :: xcomp
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nzdim,nosize) :: ycomp
REAL, INTENT(OUT), OPTIONAL, DIMENSION(ncloc,nrloc,nzdim,nosize) :: vecmag, &
                                                                  & vecpha

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*xcomp*     REAL    X-component of input vector at the U-nodes
!*ycomp*     REAL    Y-component of input vector at the V-nodes
!*vecmag*    REAL    Returned vector magnitude at the V-nodes
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
!*nzdim*     INTEGER Third dimension of input vector
!*nosize*    INTEGER Fourth dimension of input vector
!*ivarid*    INTEGER Variable id
!*info*      LOGICAL Disables/enables writing of log info
!*outflag*   REAL    Output flag for dry points
!*hregular*  LOGICAL Flag to select uniform or area averaging in the horizontal
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: hreg
REAL :: flagval
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskvals
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: array4dv


procname(pglev+1) = 'vector_mag_arr_atv'
CALL log_timer_in()

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
!2. Allocate arrays
!------------------
!

ALLOCATE (maskvals(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskvals',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (array4dv(ncloc,nrloc,nzdim,nosize),STAT=errstat)
CALL error_alloc('array4dv',4,(/ncloc,nrloc,nzdim,nosize/),kndrtype)

!
!3. Interpolate components at V-nodes
!-------------------------------------
!

CALL Uarr_at_V(xcomp,array4dv,intsrce,intdest,(/1,0,1/),&
            & (/ncloc+1,nrloc,nzdim/),nosize,ivarid,info,flagval)

!
!4. Mask array
!-------------
!

SELECT CASE (intdest)
   CASE (0); maskvals = .TRUE.
   CASE (1); maskvals = node2dv(1:ncloc,1:nrloc).GT.0
   CASE (2); maskvals = node2dv(1:ncloc,1:nrloc).EQ.2
   CASE (3); maskvals = node2dv(1:ncloc,1:nrloc).GT.1
   CASE (4); maskvals = node2dv(1:ncloc,1:nrloc).EQ.1.OR.&
                      & node2dv(1:ncloc,1:nrloc).EQ.2
END SELECT

!
!5. Vector magnitude
!-------------------
!

IF (.NOT.PRESENT(vecpha)) THEN
   CALL complex_polar(array4dv,ycomp,xamp=vecmag,maskvals=maskvals,&
                    & outflag=flagval)
ELSEIF (.NOT.PRESENT(vecmag)) THEN
   CALL complex_polar(array4dv,ycomp,xpha=vecpha,maskvals=maskvals,&
                    & outflag=flagval)
ELSE
   CALL complex_polar(array4dv,ycomp,xamp=vecmag,xpha=vecpha,&
                    & maskvals=maskvals,outflag=flagval)
ENDIF

!
!6. Deallocate
!-------------
!

DEALLOCATE (maskvals,array4dv)

CALL log_timer_out()


RETURN

END SUBROUTINE vector_mag_arr_atv

!========================================================================

SUBROUTINE vector_mag_var_atc(xcomp,ycomp,i,j,k,intsrce,intdest,vecmag,vecpha,&
                            & outflag)
!************************************************************************
!
! *vector_mag_var_atc* Calculate the (horizontal) magnitude and/or phase
!                      of a vector at a C-node grid location
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.5
!
! Description -
!
! Reference - 
!
! Module calls - complex_polar, Uvar_at_C, Vvar_at_C
!
!************************************************************************
!
USE array_interp, ONLY: Uvar_at_C, Vvar_at_C 

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2) :: xcomp
REAL, INTENT(IN), DIMENSION(2) :: ycomp
REAL, INTENT(OUT), OPTIONAL :: vecmag, vecpha

!
! Name         Type    Purpose
!-----------------------------------------------------------------------------
!*xcomp*       REAL    X-component of input vector at the U-node
!*ycomp*       REAL    Y-component of input vector at the V-node
!*i*           INTEGER X-index of input vector on the model grid
!*j*           INTEGER Y-index of input vector on the model grid
!*k*           INTEGER Vertical index of input vector on the model grid
!*intsrce*     INTEGER Selects valid points at source node
!                 = 0 => all points
!                 = 1 => coastal boundaries, interior points and open
!                        boundaries only
!                 = 2 => interior points only
!                 = 3 => interior points and open boundaries only
!                 = 4 => coastal boundaries and interior points only
!*intdest*     INTEGER Selects valid points at destination node
!                 = 0 => all points
!                 = 1 => wet points only
!*outflag*     REAL    Output flag for dry points
!*vecmag*      REAL    Returned vector magnitude at C-node
!*vecpha*      REAL    Returned vector phase at C-node                    [rad]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL :: flagval, xcompatc, ycompatc


!
!1. Optional argument
!--------------------
!

IF (.NOT.PRESENT(outflag)) THEN
   flagval = 0.0
ELSE
   flagval = outflag
ENDIF

!
!2. Interpolate components at C-nodes
!-------------------------------------
!

xcompatc = Uvar_at_C(xcomp,i,j,k,intsrce,intdest,flagval)
ycompatc = Vvar_at_C(ycomp,i,j,k,intsrce,intdest,flagval)

!
!3. Vector magnitude
!-------------------
!

IF (.NOT.PRESENT(vecpha)) THEN
   CALL complex_polar(xcompatc,ycompatc,xamp=vecmag)
ELSEIF(.NOT.PRESENT(vecmag)) THEN
   CALL complex_polar(xcompatc,ycompatc,xpha=vecpha)
ELSE
   CALL complex_polar(xcompatc,ycompatc,xamp=vecmag,xpha=vecpha)
ENDIF


RETURN

END SUBROUTINE vector_mag_var_atc

!========================================================================

SUBROUTINE vector_mag_var_atu(xcomp,ycomp,i,j,k,intsrce,intdest,vecmag,vecpha,&
                            & outflag,hregular)
!************************************************************************
!
! *vector_mag_var_atu* Calculate the (horizontal) magnitude and/or phase
!                      of a vector at a U-node grid location
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.5
!
! Description -
!
! Reference - 
!
! Module calls - complex_polar, Uvar_at_V
!
!************************************************************************
!
USE gridpars
USE switches
USE array_interp, ONLY: Vvar_at_U 

!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN) :: xcomp
REAL, INTENT(IN), DIMENSION(2,2) :: ycomp
REAL, INTENT(OUT), OPTIONAL :: vecmag, vecpha

!
! Name         Type    Purpose
!-----------------------------------------------------------------------------
!*xcomp*       REAL    X-component of input vector at the U-node
!*ycomp*       REAL    Y-component of input vector at the V-node
!*i*           INTEGER X-index of input vector on the model grid
!*j*           INTEGER Y-index of input vector on the model grid
!*k*           INTEGER Vertical index of input vector on the model grid
!*intsrce*     INTEGER Selects valid points at source node
!                 = 0 => all points
!                 = 1 => coastal boundaries, interior points and open
!                        boundaries only
!                 = 2 => interior points only
!                 = 3 => interior points and open boundaries only
!                 = 4 => coastal boundaries and interior points only
!*intdest*     INTEGER Selects valid points at destination node
!                  = 0 => all points
!                  = 1 => coastal boundaries, interior points and open
!                         boundaries only
!                  = 2 => interior points only
!                  = 3 => interior points and open boundaries only
!                  = 4 => coastal boundaries and interior points only
!*vecmag*      REAL    Returned vector magnitude at the U-node
!*vecpha*      REAL    Returned vector phase at the U-node                [rad]
!*outflag*     REAL    Output flag for dry points
!*hregular*    LOGICAL Flag to select uniform or area averaging in the
!                      horizontal
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: hreg
REAL :: flagval, ycompatu


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
!2. Interpolate components at C-nodes
!-------------------------------------
!

ycompatu = Vvar_at_U(ycomp,i,j,k,intsrce,intdest,flagval,hreg)

!
!3. Vector magnitude
!-------------------
!

IF (.NOT.PRESENT(vecpha)) THEN
   CALL complex_polar(xcomp,ycompatu,xamp=vecmag)
ELSEIF(.NOT.PRESENT(vecmag)) THEN
   CALL complex_polar(xcomp,ycompatu,xpha=vecpha)
ELSE
   CALL complex_polar(xcomp,ycompatu,xamp=vecmag,xpha=vecpha)
ENDIF


RETURN

END SUBROUTINE vector_mag_var_atu

!========================================================================

SUBROUTINE vector_mag_var_atv(xcomp,ycomp,i,j,k,intsrce,intdest,vecmag,vecpha,&
                            & outflag,hregular)
!************************************************************************
!
! *vector_mag_var_atv* Calculate the (horizontal) magnitude and/or phase
!                      of a vector at a V-node grid location
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)math_library.f90  V2.5
!
! Description -
!
! Reference - 
!
! Module calls - complex_polar, Vvar_at_U
!
!************************************************************************
!
USE gridpars
USE switches
USE array_interp, ONLY: Uvar_at_V 

!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: hregular
INTEGER, INTENT(IN) :: i, intdest, intsrce, j, k
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(2,2) :: xcomp
REAL, INTENT(IN) :: ycomp
REAL, INTENT(OUT), OPTIONAL :: vecmag, vecpha

!
! Name         Type    Purpose
!-----------------------------------------------------------------------------
!*xcomp*       REAL    X-component of input vector at the U-node
!*ycomp*       REAL    Y-component of input vector at the V-node
!*i*           INTEGER X-index of input vector on the model grid
!*j*           INTEGER Y-index of input vector on the model grid
!*k*           INTEGER Vertical index of input vector on the model grid
!*intsrce*     INTEGER Selects valid points at source node
!                 = 0 => all points
!                 = 1 => coastal boundaries, interior points and open
!                        boundaries only
!                 = 2 => interior points only
!                 = 3 => interior points and open boundaries only
!                 = 4 => coastal boundaries and interior points only
!*intdest*     INTEGER Selects valid points at destination node
!                  = 0 => all points
!                  = 1 => coastal boundaries, interior points and open
!                         boundaries only
!                  = 2 => interior points only
!                  = 3 => interior points and open boundaries only
!                  = 4 => coastal boundaries and interior points only
!*vecmag*      REAL    Returned vector magnitude at the C-node
!*vecpha*      REAL    Returned vector phase at the C-node                [rad]
!*outflag*     REAL    Output flag for dry points
!*hregular*    LOGICAL Flag to select uniform or area averaging in the
!                      horizontal
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: hreg
REAL :: flagval, xcompatv


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
!2. Interpolate components at C-nodes
!-------------------------------------
!

xcompatv = Uvar_at_V(xcomp,i,j,k,intsrce,intdest,flagval,hreg)

!
!3. Vector magnitude
!-------------------
!

IF (.NOT.PRESENT(vecpha)) THEN
   CALL complex_polar(xcompatv,ycomp,xamp=vecmag)
ELSEIF(.NOT.PRESENT(vecmag)) THEN
   CALL complex_polar(xcompatv,ycomp,xpha=vecpha)
ELSE
   CALL complex_polar(xcompatv,ycomp,xamp=vecmag,xpha=vecpha)
ENDIF


RETURN

END SUBROUTINE vector_mag_var_atv


END MODULE math_library
