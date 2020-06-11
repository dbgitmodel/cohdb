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

MODULE nla_library
!************************************************************************
!
! *nla_library* Linear algebra routine library
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)nla_library.f90  V2.11.2
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
! Routines - cholesky_decomp, cholesky_solve, LU_decomp, LU_solve, svd_decomp,
!            symm_eigen, tridiag_reduce, tridiag_vert, tridiag_vert_1d
!
!************************************************************************
!
USE iopars
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc


CONTAINS

!========================================================================

SUBROUTINE cholesky_decomp(a,n,ndim,varname,ierr,info)
!************************************************************************
!
! *cholesky_decomp* Cholesky decomposition of a symmetric positive
!                   definite matrix "a"
!
! Author -
!
! Version - @(COHERENS)nla_library.f90  V2.0
!
! Description - Cholesky matrix is returned in the lower triangular
!               and diagonal part of "a"
!
! Reference - Numerical Recipes (routine "choldc")
!
! Module calls - error_proc
!
!************************************************************************
!
USE error_routines, ONLY: error_proc

!
!*Arguments
!
LOGICAL, INTENT(IN) :: info
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: n, ndim
INTEGER, INTENT(OUT) :: ierr
REAL, INTENT(INOUT), DIMENSION(ndim,ndim) :: a

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*a*       REAL    On input, symmetric positive definite matrix
!                  On output, Cholesky matrix
!*n*       INTEGER Logical dimension of square matrix 'a'
!*ndim*    INTEGER Physical dimension of square matrix 'a'
!*varname* CHAR    FORTRAN name of array
!*ierr*    INTEGER Error output code
!                  = 0 => exited normally
!                  = 1 => singular matrix
!*info*    LOGICAL Disables/enables writing of log info
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i
REAL :: summ


procname(pglev+1) = 'cholesky_decomp'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(varname),&
                & info=info)

!
!1. Case n = 1
!-------------
!

IF (n.EQ.1) THEN
   IF (a(1,1).EQ.0.0) ierr = 1
   GOTO 1000
ENDIF

!
!2. Construct Cholesky matrix
!----------------------------
!

i_210: DO i=1,n
   summ = a(i,i) - DOT_PRODUCT(a(i,1:i-1),a(i,1:i-1))
   IF (summ.LE.0.0) THEN
      ierr = 1
      GOTO 1001
   ENDIF
   a(i,i) = SQRT(summ)
   a(i+1:n,i) = (a(i,i+1:n)-MATMUL(a(i+1:n,1:i-1),a(i,1:i-1)))/a(i,i)
ENDDO i_210

1000 CALL log_timer_out(npcc,itm_libs,info)


RETURN

1001 CONTINUE
IF (errchk) THEN
   WRITE (ioerr,'(A)') 'Unable to compute Cholesky decomposition of array: '&
                     & //TRIM(varname)
ENDIF
CALL error_proc(ierr,abort=.TRUE.)

RETURN

END SUBROUTINE cholesky_decomp

!========================================================================

SUBROUTINE cholesky_solve(a,b,n,ndim,varname,info)
!************************************************************************
!
! *cholesky_solve* Solve a linear system of equations using Cholesky
!                  decomposition
!
! Author -
!
! Version - @(COHERENS)nla_library.f90  V2.0
!
! Description - input matrix is the Cholesky matrix obtained from a previous
!               call to cholesky_decomp
!
! Reference - Numerical Recipes (routine "cholsl")
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: info
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: n, ndim
REAL, INTENT(INOUT), DIMENSION(ndim) :: b
REAL, INTENT(IN), DIMENSION(ndim,ndim) :: a

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*a*       REAL    Cholesky matrix returned by cholesky_decomp
!*b*       REAL    R.h.s. on input, solution vector on output
!*n*       INTEGER Number of equations (logical dimension of square matrix 'a')
!*ndim*    INTEGER Physical dimension of square matrix 'a'
!*varname* CHAR    FORTRAN name of array
!*info*    LOGICAL Disables/enables writing of log info
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i


procname(pglev+1) = 'cholesky_solve'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(varname),&
                & info=info)

!
!1. Case n = 1
!-------------
!

IF (n.EQ.1) THEN
   b(1) = b(1)/a(1,1)
   GOTO 1000
ENDIF

!
!2. Solve
!--------
!

i_210: DO i=1,n
   b(i) = (b(i)-DOT_PRODUCT(a(i,1:i-1),b(1:i-1)))/a(i,i)
ENDDO i_210
i_220: DO i=n,1,-1
   b(i) = (b(i)-DOT_PRODUCT(a(i+1:n,i),b(i+1:n)))/a(i,i)
ENDDO i_220

1000 CALL log_timer_out(npcc,itm_libs,info)


RETURN

END SUBROUTINE cholesky_solve

!========================================================================

SUBROUTINE LU_decomp(a,indx,n,ndim,varname,ierr,info)
!************************************************************************
!
! *LU_decomp* Decompose the matrix 'a' into a lower (L) and
!             upper (U) triangular part
!
! Author -
!
! Version - @(COHERENS)nla_library.f90  V2.0
!
! Description - replaces the matrix a by its LU decomposition using partial
!               pivoting
!             - on output L is stored into the lower triangular part of "a"
!               (excluding diagonal), U is stored into the upper triangular
!               part (including diagonal) 
!             - indx is an output vector which records the row permutation
!               effected by partial pivoting
!             - the routine can be used in combination with LU_solve to solve
!               system(s) of linear equations with matrix "a"
!
! Reference - Numerical Recipes (routine "ludcmp")
!
! Module calls - error_proc, swap_data
!
!************************************************************************
!
USE error_routines, ONLY: error_proc
USE utility_routines, ONLY : swap_data

!
!*Arguments
!
LOGICAL, INTENT(IN) :: info
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: n, ndim
INTEGER, INTENT(OUT) :: ierr
INTEGER, INTENT(OUT), DIMENSION(ndim) :: indx
REAL, INTENT(INOUT), DIMENSION(ndim,ndim) :: a

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*a*       REAL    On input, matrix to be decomposed.
!                  On output, L is substituted below the diagonal,
!                  U at and above the diagonal
!*indx*    INTEGER Output vector recording row permutation
!*n*       INTEGER Logical dimension of square matrix 'a'
!*ndim*    INTEGER Physical dimension of square matrix 'a'
!*varname* CHAR    FORTRAN name of array
!*ierr*    INTEGER Error output code
!                  = 0 => exited normally
!                  = 1 => row of zeros in 'a'
!                  = 2 => singular matrix
!*info*    LOGICAL Disables/enables writing of log info
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: imax, j, numcp
INTEGER, DIMENSION(1) :: maxinds
REAL, DIMENSION(n) :: vv


procname(pglev+1) = 'LU_decomp'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(varname),&
                & info=info)

ierr = 0

!
!1. Case n = 1
!-------------
!

IF (n.EQ.1) THEN
   indx = 0
   IF (a(1,1).EQ.0.0) ierr = 1
   GOTO 1000
ENDIF

!
!2. Scaling factor
!-----------------
!

vv = MAXVAL(ABS(a(1:n,1:n)),DIM=2)
IF (ANY(vv.EQ.0.0)) THEN
   ierr = 2
   GOTO 1001
ENDIF
vv = 1.0/vv

!
!3. Calculate elements of L- and U-matrix
!----------------------------------------
!

j_310: DO j=1,n
   maxinds = MAXLOC(vv(j:n)*ABS(a(j:n,j)))
   imax = (j-1)+ maxinds(1)
   IF (j.NE.imax) THEN
      CALL swap_data(a(imax,1:n),a(j,1:n))
      vv(imax) = vv(j)
   ENDIF
   indx(j) = imax
   IF (a(j,j).EQ.0) THEN
      ierr = 1
      GOTO 1001
   ENDIF
   a(j+1:n,j) = a(j+1:n,j)/a(j,j)
   numcp = n-j
   a(j+1:n,j+1:n) = a(j+1:n,j+1:n) - &
                  & SPREAD(a(j+1:n,j),DIM=2,NCOPIES=numcp)*&
                  & SPREAD(a(j,j+1:n),DIM=1,NCOPIES=numcp)
ENDDO j_310

1000 CALL log_timer_out(npcc,itm_libs,info)


RETURN

1001 CONTINUE
IF (errchk) THEN
   WRITE (ioerr,'(A)') 'Unable to compute LU decomposition of array: '&
                     & //TRIM(varname)
ENDIF
CALL error_proc(ierr,abort=.TRUE.)

RETURN

END SUBROUTINE LU_decomp

!========================================================================

SUBROUTINE LU_solve(a,indx,b,n,ndim,varname,info)
!************************************************************************
!
! *LU_solve* Solve a linear system of equations using LU decomposition
!
! Author -
!
! Version - @(COHERENS)nla_library.f90  V2.0
!
! Description - solve linear system of n equations of the form ax = b
!               where a is written as the product of a lower and upper
!               triangular matrix by a previous call to LU_decomp
!
! Reference - Numerical Recipes (routine "lubksb")
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: info
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: n, ndim
INTEGER, INTENT(IN), DIMENSION(ndim) :: indx
REAL, INTENT(INOUT), DIMENSION(ndim) :: b
REAL, INTENT(IN), DIMENSION(ndim,ndim) :: a

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*a*       REAL    Coeffient matrix in LU-form returned by LU_decomp
!*indx*    INTEGER Permutation vector returned by LU_decomp
!*b*       REAL    R.h.s. on input, solution vector on output
!*n*       INTEGER Number of equations (logical dimension of square matrix 'a')
!*ndim*    INTEGER Physical dimension of square matrix 'a'
!*varname* CHAR    FORTRAN name of array
!*info*    LOGICAL Disables/enables writing of log info
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ii, ll
REAL :: summ


procname(pglev+1) = 'LU_solve'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(varname),&
                & info=info)

!
!1. Case n = 1
!-------------
!

IF (n.EQ.1) THEN
   b(1) = b(1)/a(1,1)
   GOTO 1000
ENDIF

!
!2. Forward substitution
!-----------------------
!

ii = 0
i_210: DO i=1,n
   ll = indx(i)
   summ = b(ll)
   b(ll) = b(i)
   IF (ii.NE.0) THEN
      summ = summ - DOT_PRODUCT(a(i,ii:i-1),b(ii:i-1))
   ELSEIF (summ.NE.0.0) THEN
      ii = i
   ENDIF
   b(i) = summ
ENDDO i_210
  
!
!3. Backward substitution
!------------------------
!

i_310: DO i=n,1,-1
   b(i) = (b(i)-DOT_PRODUCT(a(i,i+1:n),b(i+1:n)))/a(i,i)
ENDDO i_310

1000 CALL log_timer_out(npcc,itm_libs,info)
  
  
RETURN
  
END SUBROUTINE LU_solve

!=============================================================================

SUBROUTINE svd_decomp(a,w,v,m,n,varname,ierr,info)
!*****************************************************************************
!
! *svd_decomp* Singular value decomposition of the matrix 'a'
!
! Author -
!
! Version - @(COHERENS)nla_library.f90  V2.0
!
! Description - computes the "reduced" singular value decomposition of the mxn
!               matrix a, given by a = uwv^t
!             - u(m,n) and v(n,n) are unitary matrices
!             - the matrix u replaces a on output
!             - the diagonal matrix of singular values w is output as the
!               n-dimensional vector w
!             - the n*n matrix v (not the transpose v^t) is output as v.
!
! Reference - Numerical Recipes (routine "svdcmp")
!
! Module calls - complex_polar, error_proc, outer_product, qsort_index
!
!*****************************************************************************
!
USE error_routines, ONLY: error_proc
USE math_library, ONLY: complex_polar
USE utility_routines, ONLY: outer_product, qsort_index

!
!* Arguments
!
LOGICAL, INTENT(IN) :: info
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: m,n
INTEGER, INTENT(OUT) :: ierr
REAL, INTENT(OUT), DIMENSION(n) :: w
REAL, INTENT(INOUT), DIMENSION(m,n) :: a   
REAL, INTENT(OUT), DIMENSION(n,n) :: v

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*a*       REAL    Matrix to be decomposed on input, unitary matrix 'u'
!                  on output
!*w*       REAL    Vector of singular values (with n-m zeros at the end if n>m)
!*v*       REAL    Unitary matrix 'v' on output
!*m*       INTEGER First dimension of matrix 'a'
!*n*       INTEGER Second dimension of matrix 'a'
!*varname* CHAR    FORTRAN name of array
!*ierr*    INTEGER Error code on exit (0 for successfull)
!*info*    LOGICAL Disables/enables writing of log info
!
!------------------------------------------------------------------------------
! 
!* Local variables
!
INTEGER, PARAMETER :: itermax = 30
INTEGER :: i, iter, j, k, l, nm
REAL :: anorm, c, f, g, h, s, scale, x, y, z
INTEGER, DIMENSION(n) :: indx
REAL, DIMENSION(m) :: tempm
REAL, DIMENSION(n) :: rv, tempn


procname(pglev+1) = 'svd_decomp'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(varname),&
                & info=info)

ierr = 0

!
!1. Householder reduction to bidiagonal form
!-------------------------------------------
!

g = 0.0
scale = 0.0

i_110: DO i=1,n

   l = i + 1
   rv(i) = scale*g
   g = 0.0
   scale = 0.0

   IF (i.LE.m) THEN
      scale = SUM(ABS(a(i:m,i)))
      IF (scale.NE.0.0) THEN
         a(i:m,i) = a(i:m,i)/scale
         s = DOT_PRODUCT(a(i:m,i),a(i:m,i))
         f = a(i,i) 
         g = -SIGN(SQRT(s),f)
         h = f*g - s 
         a(i,i) = f - g
         tempn(l:n) = MATMUL(a(i:m,i),a(i:m,l:n))/h
         a(i:m,l:n) = a(i:m,l:n) + outer_product(a(i:m,i),tempn(l:n))
         a(i:m,i) = scale*a(i:m,i) 
      ENDIF
   ENDIF

   w(i) = scale*g
   g = 0.0
   scale = 0.0
   IF ((i.LE.m).AND.(i.NE.n)) THEN
      scale = SUM(ABS(a(i,l:n)))
      IF (scale.NE.0.0) THEN
         a(i,l:n) = a(i,l:n)/scale
         s = DOT_PRODUCT(a(i,l:n),a(i,l:n))
         f = a(i,l)
         g = -SIGN(SQRT(s),f)
         h = f*g - s
         a(i,l) = f - g
         rv(l:n) = a(i,l:n)/h
         tempm(l:m) = MATMUL(a(l:m,l:n),a(i,l:n))
         a(l:m,l:n) = a(l:m,l:n) + outer_product(tempm(l:m),rv(l:n))
         a(i,l:n) = scale*a(i,l:n) 
      ENDIF
   ENDIF

ENDDO i_110

anorm = MAXVAL(ABS(w)+ABS(rv))

!
!2. Accumulation of right-hand side transformations
!--------------------------------------------------
!

i_210: DO i=n,1,-1
   IF (i.LT.n) THEN
      IF (g.NE.0.0) THEN
!        ---double division to avoid possible underflow
         v(l:n,i) = (a(i,l:n)/a(i,l))/g  
         tempn(l:n) = MATMUL(a(i,l:n),v(l:n,l:n))
         v(l:n,l:n) = v(l:n,l:n) + outer_product(v(l:n,i),tempn(l:n))
      ENDIF
      v(i,l:n) = 0.0
      v(l:n,i) = 0.0
   ENDIF
   v(i,i) = 1.0
   g = rv(i)
   l = i
ENDDO i_210

!
!3. Accumulation of left-hand side transformations
!-------------------------------------------------
!

i_310: DO i=MIN(m,n),1,-1
   l = i+1
   g = w(i)
   a(i,l:n) = 0.0
   IF (g.NE.0.0) THEN
      g = 1.0/g
      tempn(l:n) = (MATMUL(a(l:m,i),a(l:m,l:n))/ a(i,i))*g
      a(i:m,l:n) = a(i:m,l:n) + outer_product(a(i:m,i),tempn(l:n))
      a(i:m,i) = a(i:m,i)*g
   ELSE
      a(i:m,i) = 0.0
   ENDIF
   a(i,i) = a(i,i) + 1.0
ENDDO i_310

!
!4. Diagonalization of the bidiagonal form
!   Loop over singular values, and over allowed iterations
!---------------------------------------------------------
!

k_400: DO k=n,1,-1
   iter_410: DO iter=1,itermax

!     ---test for splitting
      l_411: DO l=k,1,-1
         nm = l-1
         IF((ABS(rv(l))+anorm).EQ.anorm) EXIT l_411
!        ---Note that rv(1) is always zero, so can never fall through bottom
!           of loop
         IF ((ABS(w(nm))+anorm).EQ.anorm) THEN
!           ---cancellation of rv(l) if l > 1
            c = 0.0
            s = 1.0
            i_4111: DO i=l,k
               f = s*rv(i)
               rv(i) = c*rv(i)
               IF ((ABS(f)+anorm).EQ.anorm) EXIT i_4111
               g = w(i)
               CALL complex_polar(f,g,xamp=h)
               w(i) = h
               h = 1.0/h
               c = g*h
               s = -(f*h)
               tempm(1:m) = a(1:m, nm)
               a(1:m,nm) = a(1:m,nm)*c + a(1:m,i)*s
               a(1:m,i) = -tempm(1:m)*s + a(1:m,i)*c 
            ENDDO i_4111
            EXIT l_411
         ENDIF
      ENDDO l_411

      z = w(k)
!     ---convergence
      IF (l.EQ.k) THEN
!        --singular value is made non negative
         IF (z.LT.0.0) THEN
            w(k) = -z
            v(1:n,k) = -v(1:n,k)
         ENDIF
         EXIT iter_410
      ENDIF
      IF (iter.EQ.itermax) THEN
         ierr = 1
         GOTO 1001
      ENDIF

!     ---shift from bottom 2-by-2 minor
      x = w(l)
      nm = k-1
      y = w(nm)
      g = rv(nm)
      h = rv(k)
      f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)      
      CALL complex_polar(f,1.0,xamp=g)
      f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x

!     ---next QR transformation
      c = 1.0
      s = 1.0

      j_412: DO j=l,nm
         i = j+1
         g = rv(i)
         y = w(i)
         h = s*g
         g = c*g
         CALL complex_polar(f,h,xamp=z)
         rv(j) = z
         c = f/z
         s = h/z
         f = x*c + g*s
         g = g*c - x*s
         h = y*s
         y = y*c
         tempn(1:n) = v(1:n,j)
         v(1:n,j) = v(1:n,j)*c + v(1:n,i)*s
         v(1:n,i) = v(1:n,i)*c - tempn(1:n)*s
         CALL complex_polar(f,h,xamp=z)
         w(j) = z
!        ---rotation can be arbitrary if z = 0
         IF (z.NE.0.0) THEN
            z = 1.0/z
            c = f*z
            s = h*z
         ENDIF
         f = c*g + s*y
         x = c*y - s*g
         tempm(1:m) = a(1:m,j)
         a(1:m,j) = a(1:m,j)*c + a(1:m,i)*s
         a(1:m,i) = a(1:m,i)*c - tempm(1:m)*s
      ENDDO j_412

      rv(l) = 0.0
      rv(k) = f
      w(k) = x

   ENDDO iter_410
ENDDO k_400

!
!5. Sort singular values in descending order
!-------------------------------------------
!

CALL qsort_index(w,indx,-1)
tempn = w
w = tempn(indx)
l_510: DO l=1,m
   tempn = a(l,:)
   a(l,:) = tempn(indx)
ENDDO l_510
l_520: DO l=1,n
   tempn = v(l,:)
   v(l,:) = tempn(indx)
ENDDO l_520

CALL log_timer_out(npcc,itm_libs,info)


RETURN

1001 CONTINUE
IF (errchk) THEN
   WRITE (ioerr,'(A)') 'Unable to compute singular value decomposition of&
                      & array: '//TRIM(varname)
ENDIF
CALL error_proc(ierr,abort=.TRUE.)

RETURN

END SUBROUTINE svd_decomp

!========================================================================

SUBROUTINE symm_eigen(a,d,n,vectors,varname,ierr,info)
!************************************************************************
!
! *symm_eigen* Eigenvalues of a real symmetric matrix
!
! Author -
!
! Version - @(COHERENS)nla_library.f90  V2.0
!
! Description - input matrix is first reduced to tridiagonal form by calling
!               tridiag_reduce
!
! Calling program -
!
! Reference - Numerical recipes (routine "tqli")
!
! Module calls - complex_polar, error_proc, qsort_index, tridiag_reduce
!
!************************************************************************
!
USE error_routines, ONLY: error_proc
USE math_library, ONLY: complex_polar
USE utility_routines, ONLY: qsort_index

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info, vectors
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(OUT) :: ierr
REAL, INTENT(OUT), DIMENSION(n) :: d
REAL, INTENT(INOUT), DIMENSION(n,n) :: a 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*a*       REAL    On input, the symmetric matrix 'a'.
!                  On output the matrix of eigenvectors if vectors  = .TRUE.,
!                  or contains no useful information if vectors  = .FALSE.
!*d*       REAL    Vector of eigenvalues
!*n*       INTEGER Dimension of matrix 'a'
!*vectors* LOGICAL If .TRUE. eigenvalues and eigenvectors have to found.
!                  If .FALSE. only eigenvalues are needed.
!*varname* CHAR    FORTRAN name of array
!*ierr*    INTEGER Error code on exit (0 for successfull)
!*info*    LOGICAL Disables/enables writing of log info
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, PARAMETER :: itermax = 30
INTEGER :: i, iter, l, m
REAL :: b, c, dd, f, g, p, r, s
INTEGER, DIMENSION(n) :: indx
REAL, DIMENSION(n) :: e, ff


procname(pglev+1) = 'symm_eigen'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(varname),&
                & info=info)

!
!1. Reduce matrix to tridiagonal form
!------------------------------------
!

CALL tridiag_reduce(a,d,e,n,vectors,varname,info)
e = EOSHIFT(e,1)

!
!2. Eigenvalues and eigenvectors by iteration
!--------------------------------------------
!

l_200: DO l=1,n
   iter = 0
   iter_210: DO
      m_211: DO m=l,n-1
         dd = ABS(d(m)) + ABS(d(m+1))
         IF ((ABS(e(m))+dd).EQ.dd) EXIT m_211
      ENDDO m_211
      IF (m.EQ.l) EXIT iter_210
      IF (iter.EQ.itermax) THEN
         ierr = 1
         GOTO 1001
      ENDIF
      iter = iter + 1
      g = (d(l+1)-d(l))/(2.0*e(l))
      CALL complex_polar(g,1.0,xamp=r)
      g = d(m) - d(l) + e(l)/(g+SIGN(r,g))
      s = 1.0
      c = 1.0
      p = 0.0
      i_212: DO i=m-1,l,-1
         f = s*e(i)
         b = c*e(i)
         CALL complex_polar(f,g,xamp=r)
         e(i+1) = r
         IF (r.EQ.0.0) THEN
            d(i+1) = d(i+1) - p
            e(m) = 0.0
            CYCLE iter_210
         ENDIF
         s = f/r
         c = g/r
         g = d(i+1) - p
         r = (d(i)-g)*s+2.0*c*b
         p = s*r
         d(i+1) = g + p
         g = c*r - b
         IF (vectors) THEN
            ff = a(:,i+1)
            a(:,i+1) = s*a(:,i) + c*ff
            a(:,i) = c*a(:,i) - s*ff
         ENDIF
      ENDDO i_212
      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0
   ENDDO iter_210
ENDDO l_200

!
!3. Sort eigenvalues/eigenvectors in descending order
!----------------------------------------------------
!

CALL qsort_index(d,indx,-1)
ff = d
d = ff(indx)

IF (vectors) THEN
   l_310: DO l=1,n
      ff = a(l,:)
      a(l,:) = ff(indx)
   ENDDO l_310
ENDIF

CALL log_timer_out(npcc,itm_libs,info)


RETURN

1001 CONTINUE
IF (errchk) THEN
   WRITE (ioerr,'(A)') 'Unable to compute eigenvalue decomposition of&
                      & array: '//TRIM(varname)
ENDIF
CALL error_proc(ierr,abort=.TRUE.)

RETURN

END SUBROUTINE symm_eigen

!========================================================================

SUBROUTINE tridiag_reduce(a,d,e,n,vectors,varname,info)
!************************************************************************
!
! *tridiag_reduce* Householder reduction of a real symmetric matrix to
!                  tridiagonal form
!
! Author -
!
! Version - @(COHERENS)nla_library.f90  V2.0
!
! Description -
!
! Calling program - tridiag_eigen
!
! Reference - Numerical recipes (routine "tred2")
!
! Module calls - outer_product
!
!************************************************************************
!
USE utility_routines, ONLY: outer_product

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: info
CHARACTER (LEN=*), INTENT(IN) :: varname
LOGICAL, INTENT(IN) :: vectors 
INTEGER, INTENT(IN) :: n
REAL, INTENT(OUT), DIMENSION(n) :: d, e
REAL, INTENT(INOUT), DIMENSION(n,n) :: a 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*a*       REAL    On input, the symmetric matrix 'a'.
!                  On output the orthogonal matrix effecting the transformation
!*d*       REAL    Diagonal elements of the output tridiagonal matrix
!*e*       REAL    Off-diogonal elements of the output tridiagonal matrix
!                  (e(1) = 0)
!*n*       INTEGER Dimension of matrix 'a'
!*vectors* LOGICAL If .TRUE. eigenvalues and eigenvectors have to found.
!                  If .FALSE. only eigenvalues are needed.
!*varname* CHAR    FORTRAN name of matrix
!*info*    LOGICAL Disables/enables writing of log info
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j, l
REAL :: f, g, h, hh, scale
REAL, DIMENSION(n) :: gg


procname(pglev+1) = 'tridiag_reduce'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(varname),&
                & info=info)

!
!1. Reduce matrix
!----------------
!

i_110: DO i=n,2,-1
   l = i - 1
   h = 0.0
   IF (l.GT.1) THEN
      scale = SUM(ABS(a(i,1:l)))
      IF (scale.EQ.0.0) THEN
         e(i) = a(i,l)
      ELSE
         a(i,1:l) = a(i,1:l)/scale
         h = SUM(a(i,1:l)**2)
         f = a(i,l)
         g = -SIGN(SQRT(h),f)
         e(i) = scale*g
         h = h - f*g
         a(i,l) = f - g
         IF (vectors) a(1:l,i) = a(i,1:l)/h
         j_111: DO j=1,l
            e(j) = (DOT_PRODUCT(a(j,1:j),a(i,1:j)) + &
                  & DOT_PRODUCT(a(j+1:l,j),a(i,j+1:l)))/h
         ENDDO j_111
         f = DOT_PRODUCT(e(1:l),a(i,1:l))
         hh = f/(h+h)
         e(1:l) = e(1:l) - hh*a(i,1:l)
         j_112: DO j=1,l
            a(j,1:j) = a(j,1:j) - a(i,j)*e(1:j) - e(j)*a(i,1:j)
         ENDDO j_112
      ENDIF
   ELSE
      e(i) = a(i,l)
   ENDIF
   d(i) = h
ENDDO i_110

!
!2. Orthogonal transformation matrix
!-----------------------------------
!

IF (vectors) d(1) = 0.0
e(1) = 0.0
i_210: DO i=1,n
   IF (vectors) THEN
      l = i - 1
      IF (d(i).NE.0.0) THEN
         gg(1:l) = MATMUL(a(i,1:l),a(1:l,1:l))
         a(1:l,1:l) = a(1:l,1:l) - outer_product(a(1:l,i),gg(1:l))
      ENDIF
      d(i) = a(i,i)
      a(i,i) = 1.0
      a(i,1:l) = 0.0
      a(1:l,i) = 0.0
   ELSE
      d(i) = a(i,i)
   ENDIF
ENDDO i_210

CALL log_timer_out(npcc,itm_libs,info)


RETURN

END SUBROUTINE tridiag_reduce

!========================================================================

SUBROUTINE tridiag_vert(tridcf,psi,nhdims,nzdim,novars,cnode)
!************************************************************************
!
! *tridiag_vert* Solve a vertical tridiagonal system of linear equations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)nla_library.f90  V2.11.2
!
! Description -
!
! Reference - Numerical Recipes (routine "tridag")
!
! Module calls - error_alloc
!
!************************************************************************
!
USE grid
USE gridpars
USE syspars
USE error_routines, ONLY: error_alloc

!
!*Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: novars, nzdim
INTEGER, INTENT(IN), DIMENSION(4) :: nhdims
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nzdim,4,novars) :: tridcf
REAL, INTENT(INOUT), DIMENSION(1-nhdims(1):ncloc+nhdims(2),&
                             & 1-nhdims(3):nrloc+nhdims(4),nzdim,novars) :: psi

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*tridcf*   REAL     Coefficients of tridiagonal matrix system
!*psi*      REAL     Solution of tridiagonal system
!*nhdims*   INTEGER  Halo sizes of solution array (WESN directions)
!*nzdim*    INTEGER  Vertical array dimension
!*novars*   INTEGER  Number of variables
!*cnode*    CHAR     Type of node ('C', 'U', 'V')
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar, k
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskvals
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: beta
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: gamma


procname(pglev+1) = 'tridiag_vert'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskvals(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskvals',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (beta(ncloc,nrloc),STAT=errstat)
CALL error_alloc('beta',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (gamma(ncloc,nrloc,nzdim),STAT=errstat)
CALL error_alloc('gamma',3,(/ncloc,nrloc,nzdim/),kndrtype)

!
!2. Mask array
!-------------
!

SELECT CASE (TRIM(cnode))
   CASE ('C')
      maskvals = nodeatc(1:ncloc,1:nrloc).GT.0
   CASE ('U')
      maskvals = node2du(1:ncloc,1:nrloc).LE.2
   CASE ('V')
      maskvals = node2dv(1:ncloc,1:nrloc).LE.2
END SELECT

!
!3. Solve tridiagonal system
!---------------------------
!

ivar_310: DO ivar=1,novars
   WHERE (maskvals)
      beta = tridcf(:,:,1,2,ivar)
      psi(1:ncloc,1:nrloc,1,ivar) = tridcf(:,:,1,4,ivar)/beta
   END WHERE
   k_311: DO k=2,nzdim
      WHERE (maskvals)
         gamma(:,:,k) = tridcf(:,:,k-1,3,ivar)/beta
         beta = tridcf(:,:,k,2,ivar) - tridcf(:,:,k,1,ivar)*gamma(:,:,k)
         psi(1:ncloc,1:nrloc,k,ivar) = (tridcf(:,:,k,4,ivar)-&
                                      & tridcf(:,:,k,1,ivar)*&
                                      & psi(1:ncloc,1:nrloc,k-1,ivar))/beta
      END WHERE
   ENDDO k_311
   k_312: DO k=nzdim-1,1,-1
      WHERE (maskvals)
         psi(1:ncloc,1:nrloc,k,ivar) = psi(1:ncloc,1:nrloc,k,ivar) -&
                                     & gamma(:,:,k+1)*&
                                     & psi(1:ncloc,1:nrloc,k+1,ivar)
      END WHERE
   ENDDO k_312
ENDDO ivar_310

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (beta,gamma,maskvals)

CALL log_timer_out(npcc,itm_libs)


RETURN

END SUBROUTINE tridiag_vert

!========================================================================

SUBROUTINE tridiag_vert_1d(tridcf,psi,nzdim,novars,info)
!************************************************************************
!
! *tridiag_vert_1d* Solve a vertical tridiagonal system of linear equations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)nla_library.f90  V2.1.2
!
! Description - Equations are solved at one horizontal location
!
! Reference - Numerical Recipes (routine "tridag")
!
! Module calls - error_alloc
!
!************************************************************************
!
USE grid
USE gridpars

!
!*Arguments
!
LOGICAL, INTENT(IN) :: info
INTEGER, INTENT(IN) :: novars, nzdim
REAL, INTENT(IN), DIMENSION(nzdim,4,novars) :: tridcf
REAL, INTENT(INOUT), DIMENSION(nzdim,novars) :: psi

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*tridcf*   REAL     Coefficients of tridiagonal matrix system
!*psi*      REAL     Solution of tridiagonal system
!*nzdim*    INTEGER  Vertical array dimension
!*novars*   INTEGER  Number of variables
!*info*     INTEGER  Writes log info if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar, k
REAL :: beta
REAL, DIMENSION(nzdim) :: gamma


procname(pglev+1) = 'tridiag_vert_1d'
CALL log_timer_in(npcc,info=info)


ivar_110: DO ivar=1,novars
   beta = tridcf(1,2,ivar)
   psi(1,ivar) = tridcf(1,4,ivar)/beta
   k_111: DO k=2,nzdim
      gamma(k) = tridcf(k-1,3,ivar)/beta
      beta = tridcf(k,2,ivar) - tridcf(k,1,ivar)*gamma(k)
      psi(k,ivar) = (tridcf(k,4,ivar)-tridcf(k,1,ivar)*psi(k-1,ivar))/beta
   ENDDO k_111
   k_112: DO k=nzdim-1,1,-1
      psi(k,ivar) = psi(k,ivar) - gamma(k+1)*psi(k+1,ivar)
   ENDDO k_112
ENDDO ivar_110

CALL log_timer_out(npcc,itm_libs,info=info)


RETURN

END SUBROUTINE tridiag_vert_1d


END MODULE nla_library
