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

MODULE fft_library
!************************************************************************
!
! *fft_library* Library for computing fast Fourier transforms
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)fft_library.f90  V2.0
!
! $Date: 2016-04-04 09:25:08 +0200 (Mon, 04 Apr 2016) $
!
! $Revision: 931 $
!
! Description -
!
! Reference - Press W.H., Flannery B.P., Teukolsky S.A. and Vetterling W.T.,
!             1992. Numerical Recipes. The art of scientific computing.
!             Cambridge University Press, Cambridge, 702 pp.
!
! Routines - fft_fourrow2, fft_fourrow3, fft_four1, fft_four2, fft_four3,
!            fft_realft1, fft_realft2, fft_realft3
!
!************************************************************************
!
USE iopars
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


CONTAINS

!===========================================================================

SUBROUTINE fft_fourrow2(data,ndims,iform)
!**************************************************************************
!
! *fft_fourrow2* 2-D Fourier transform
!
! Author -
!
! Version - @(COHERENS)fft_library.f90  V2.0
!
! Description 
!  - Replaces each row (constant first index) of "data(1:M,1:N)" by its
!    discrete Fourier transform (transform on second index), if iform = 1;
!    or replaces each row of data by N times its inverse discrete Fourier
!    transform, if iform = -1.
!  - N must be an integer power of 2.
!
! Reference - Numerical recipes (routine "fourrow")
!
! Module calls - swap_data
!
!*****************************************************************************
!
USE utility_routines, ONLY: swap_data

!
!*Arguments
!
INTEGER, INTENT(IN) :: iform
INTEGER, INTENT(IN), DIMENSION(2) :: ndims
COMPLEX, INTENT(INOUT), DIMENSION(ndims(1),ndims(2)) :: data

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*data*     COMPLEX Transformed data on output
!*ndims*    INTEGER Shape of array data
!*iform*    INTEGER Type of transform
!                   =  1 : forward
!                   = -1 : inverse
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, istep, j, m, mmax, n, n2
REAL (KIND=kndrlong) :: theta
COMPLEX :: ws
COMPLEX (KIND=kndrlong) :: w, wp
COMPLEX, DIMENSION(ndims(1)) :: tmp


procname(pglev+1) = 'fft_fourrow2'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!

n = ndims(2)
n2 = n/2
j = n2

!
!2. Bit-reversal
!---------------
!

i_210: DO i=1,n-2
   IF (j.GT.i) CALL swap_data(data(:,j+1),data(:,i+1))
   m = n2
   DO WHILE ((m.GE.2).AND.(j.GE.m))
      j = j - m; m = m/2
   ENDDO
   j = j + m
ENDDO i_210

!
!3. Danielson-Lanczos algorithm
!------------------------------
!

mmax = 1
DO WHILE (n.GT.mmax)
   istep = 2*mmax
   theta = pi_d/(iform*mmax)
   wp = CMPLX(-2.0_kndrlong*SIN(0.5_kndrlong*theta)**2,SIN(theta),KIND=kndrlong)
   w = CMPLX(1.0_kndrlong,KIND=kndrlong)
   m_310: DO m=1,mmax
      ws = w
      i_311: DO i=m,n,istep
         j = i + mmax
         tmp = ws*data(:,j)
         data(:,j) = data(:,i) - tmp
         data(:,i) = data(:,i) + tmp
      ENDDO i_311
      w = w*wp + w
   ENDDO m_310
   mmax = istep
ENDDO

CALL log_timer_out()


RETURN

END SUBROUTINE fft_fourrow2

!===========================================================================

SUBROUTINE fft_fourrow3(data,ndims,iform)
!**************************************************************************
!
! *fft_fourrow3* 3-D Fourier transform
!
! Author -
!
! Version - @(COHERENS)fft_library.f90  V2.0
!
! Description 
!  - Replaces each third-index section (constant first and second indices) of
!    "data(1:L,1:M,1:N)" by its discrete Fourier transform (transform on third
!    index), if iform = 1; or replaces each third-index section of data by N
!    times its inverse discrete Fourier transform, if iform = -1.
!  - N must be an integer power of 2.
!
! Reference - Numerical recipes (routine "fourrow_3d")
!
! Module calls - swap_data
!
!*****************************************************************************
!
USE utility_routines, ONLY: swap_data

!
!*Arguments
!
INTEGER, INTENT(IN) :: iform
INTEGER, INTENT(IN), DIMENSION(3) :: ndims
COMPLEX, INTENT(INOUT), DIMENSION(ndims(1),ndims(2),ndims(3)) :: data

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*data*     COMPLEX Transformed data on output
!*ndims*    INTEGER Shape of array data
!*iform*    INTEGER Type of transform
!                   =  1 : forward
!                   = -1 : inverse
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, istep, j, m, mmax, n, n2
REAL (KIND=kndrlong) :: theta
COMPLEX :: ws
COMPLEX (KIND=kndrlong) :: w, wp
COMPLEX, DIMENSION(ndims(1),ndims(2)) :: tmp


procname(pglev+1) = 'fft_fourrow3'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!

n = ndims(3)
n2 = n/2
j = n2

!
!2. Bit-reversal
!---------------
!

i_210: DO i=1,n-2
   IF (j.GT.i) CALL swap_data(data(:,:,j+1),data(:,:,i+1))
   m = n2
   DO WHILE ((m.GE.2).AND.(j.GE.m))
      j = j - m; m = m/2
   ENDDO
   j = j + m
ENDDO i_210

!
!3. Danielson-Lanczos algorithm
!------------------------------
!

mmax = 1
DO WHILE (n.GT.mmax)
   istep = 2*mmax
   theta = pi_d/(iform*mmax)
   wp = CMPLX(-2.0_kndrlong*SIN(0.5_kndrlong*theta)**2,SIN(theta),KIND=kndrlong)
   w = CMPLX(1.0_kndrlong,KIND=kndrlong)
   m_310: DO m=1,mmax
      ws = w
      i_311: DO i=m,n,istep
         j = i + mmax
         tmp = ws*data(:,:,j)
         data(:,:,j) = data(:,:,i) - tmp
         data(:,:,i) = data(:,:,i) + tmp
      ENDDO i_311
      w = w*wp + w
   ENDDO m_310
   mmax = istep
ENDDO

CALL log_timer_out()


RETURN

END SUBROUTINE fft_fourrow3

!===========================================================================

SUBROUTINE fft_four1(data,ndim,iform)
!**************************************************************************
!
! *fft_four1* 1-D Fourier transform
!
! Author -
!
! Version - @(COHERENS)fft_library.f90  V2.0
!
! Description 
!  - Replaces the 1-D complex array "data" by its discrete Fourier transform,
!    if iform =1; or replaces "data" by its inverse 1-D discrete Fourier
!    transform times the size of "data", if iform = -1.
!  - The size of "data" must be an integer power of 2.
!
! Reference - Numerical recipes (routine "four1")
!
! Module calls - error_alloc, fft_fourrow2
!
!*****************************************************************************
!
USE error_routines, ONLY: error_alloc

!
!*Arguments
!
INTEGER, INTENT(IN) :: iform
INTEGER, INTENT(IN) :: ndim
COMPLEX, INTENT(INOUT), DIMENSION(ndim) :: data

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*data*     COMPLEX Transformed data on output
!*ndim*     INTEGER Size of array data
!*iform*    INTEGER Type of transform
!                   =  1 : forward
!                   = -1 : inverse
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: j, k, m1, m2, n
REAL (KIND=kndrlong), ALLOCATABLE, DIMENSION(:) :: theta
COMPLEX (KIND=kndrlong), ALLOCATABLE, DIMENSION(:) :: w, wp
COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: dat, tmp


procname(pglev+1) = 'fft_four1'
CALL log_timer_in()

!---find dimensions as close to square as possible
n = ndim
m1 = 2**CEILING(0.5*LOG(REAL(n))/LOG(2.0))
m2 = n/m1

!---reshape input array
ALLOCATE (dat(m1,m2),STAT=errstat)
CALL error_alloc('dat',2,(/m1,m2/),kndcmplx)
ALLOCATE (theta(m1),STAT=errstat)
CALL error_alloc('theta',1,(/m1/),kndrlong)
ALLOCATE (w(m1),STAT=errstat)
CALL error_alloc('w',1,(/m1/),kndcmplx)
ALLOCATE (wp(m1),STAT=errstat)
CALL error_alloc('wp',1,(/m1/),kndcmplx)
ALLOCATE (tmp(m2,m1),STAT=errstat)
CALL error_alloc('tmp',2,(/m2,m1/),kndcmplx)
dat = RESHAPE(data,SHAPE(dat))

!---transform on second index
CALL fft_fourrow2(dat,SHAPE(dat),iform)

!---initialise trigonometric recurrence
theta = (/((k-1)*iform,k=1,m1)/)
theta = theta*twopi_d/REAL(n)
wp = CMPLX(-2.0_kndrlong*SIN(0.5_kndrlong*theta)**2,SIN(theta),KIND=kndrlong)
w = CMPLX(1.0_kndrlong,KIND=kndrlong)

!---multiply by the extra phase factor
j_110: DO j=2,m2
   w = w*wp + w
   dat(:,j) = dat(:,j)*w
ENDDO j_110

!---transpose and transform on (original) first index
tmp = TRANSPOSE(dat)
CALL fft_fourrow2(tmp,SHAPE(tmp),iform)

!---reshape the result back to one dimension
data = RESHAPE(tmp,SHAPE(data))
DEALLOCATE (dat,w,wp,theta,tmp)

CALL log_timer_out()


RETURN

END SUBROUTINE fft_four1

!===========================================================================

SUBROUTINE fft_four2(data,ndims,iform)
!**************************************************************************
!
! *fft_four2* 2-D Fourier transform
!
! Author -
!
! Version - @(COHERENS)fft_library.f90  V2.0
!
! Description 
!  - Replaces the 2-D complex array "data" by its discrete Fourier transform,
!    if iform =1; or replaces "data" by its inverse 2-D discrete Fourier
!    transform times the product of its two sizes, if iform = -1.
!  - Both sizes of "data" must be integer powers of 2.
!
! Reference - Numerical recipes (routine "four2")
!
! Module calls - error_alloc, fft_fourrow2
!
!*****************************************************************************
!
USE error_routines, ONLY: error_alloc

!
!*Arguments
!
INTEGER, INTENT(IN) :: iform
INTEGER, INTENT(IN), DIMENSION(2) :: ndims
COMPLEX, INTENT(INOUT), DIMENSION(ndims(1),ndims(2)) :: data

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*data*     COMPLEX Transformed data on output
!*ndims*    INTEGER Shape of array data
!*iform*    INTEGER Type of transform
!                   =  1 : forward
!                   = -1 : inverse
!------------------------------------------------------------------------------
!
!*Local variables
!
COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: data1


procname(pglev+1) = 'fft_four2'
CALL log_timer_in()

ALLOCATE (data1(ndims(2),ndims(1)),STAT=errstat)
CALL error_alloc('data1',2,(/ndims(2),ndims(1)/),kndcmplx)
!---transform in second dimension
CALL fft_fourrow2(data,ndims,iform)
!---transpose
data1 = TRANSPOSE(data)
!---transform in (original) first dimension
CALL fft_fourrow2(data1,SHAPE(data1),iform)
!---transpose into data
data = TRANSPOSE(data1)
DEALLOCATE (data1)

CALL log_timer_out()


RETURN

END SUBROUTINE fft_four2

!===========================================================================

SUBROUTINE fft_four3(data,ndims,iform)
!**************************************************************************
!
! *fft_four3* 3-D Fourier transform
!
! Author -
!
! Version - @(COHERENS)fft_library.f90  V2.0
!
! Description 
!  - Replaces the 3-D complex array "data" by its discrete Fourier transform,
!    if iform =1; or replaces "data" by its inverse 3-D discrete Fourier
!    transform times the product of its three sizes, if iform = -1.
!  - All sizes of "data" must be integer powers of 2.
!
! Reference - Numerical recipes (four3)
!
! Module calls - error_alloc, fft_fourrow3
!
!*****************************************************************************
!
USE error_routines, ONLY: error_alloc

!
!*Arguments
!
INTEGER, INTENT(IN) :: iform
INTEGER, INTENT(IN), DIMENSION(3) :: ndims 
COMPLEX, INTENT(INOUT), DIMENSION(ndims(1),ndims(2),ndims(3)) :: data

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*data*     COMPLEX Transformed data on output
!*ndims*    INTEGER Shape of array data
!*iform*    INTEGER Type of transform
!                   =  1 : forward
!                   = -1 : inverse
!
!------------------------------------------------------------------------------
!
!*Local variables
!
COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: dat2, dat3


procname(pglev+1) = 'fft_four3'
CALL log_timer_in()

!---transform in third dimension
CALL fft_fourrow3(data,ndims,iform)
!---transpose
ALLOCATE (dat2(ndims(2),ndims(3),ndims(1)),STAT=errstat)
CALL error_alloc('dat2',3,(/ndims(2),ndims(3),ndims(1)/),kndcmplx)
dat2 = RESHAPE(data,SHAPE=SHAPE(dat2),order=(/3,1,2/))
!---transform in (original) first dimension
CALL fft_fourrow3(dat2,SHAPE(dat2),iform)
!---transpose
ALLOCATE (dat3(ndims(3),ndims(1),ndims(2)),STAT=errstat)
CALL error_alloc('dat3',3,(/ndims(3),ndims(1),ndims(2)/),kndcmplx)
dat3 = RESHAPE(dat2,SHAPE=SHAPE(dat3),order=(/3,1,2/))
DEALLOCATE (dat2)
!---transform in (original) second dimension
CALL fft_fourrow3(dat3,SHAPE(dat3),iform)
!---transpose back to output order
data = RESHAPE(dat3,SHAPE=SHAPE(data),order=(/3,1,2/))
DEALLOCATE (dat3)

CALL log_timer_out()


RETURN

END SUBROUTINE fft_four3

!===========================================================================

SUBROUTINE fft_realft1(data,spec,nodat,iform)
!**************************************************************************
!
! *fft_realft1* 1-D Fourier transform
!
! Author -
!
! Version - @(COHERENS)fft_library.f90  V2.0
!
! Description 
!  - If iform=1, returns the complex fast Fourier transform
!    of a 1-D real input "data" array by the complex array "spec".
!  - If iform = -1, the inverse transform (times nodat/2) of a complex array 
!    "spec" is performed and returned as a real array "data"
!  - Array sizes must be integer powers of 2.
!
! Reference - Numerical recipes (routine "realft")
!
! Module calls - error_abort, fft_four1
!
!*****************************************************************************
!
USE error_routines, ONLY: error_abort

!
!*Arguments
!
INTEGER, INTENT(IN) :: iform, nodat
REAL, INTENT(INOUT), DIMENSION(nodat) :: data
COMPLEX, INTENT(INOUT), DIMENSION(nodat/2) :: spec

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*data*     REAL    Input data if iform = 1, output data (inverse transform)
!                   if iform = -1
!*spec*     COMPLEX Output data (forward transform) if iform = 1, input data
!                   if iform = 1
!*nodat*    INTEGER Number of data points
!*iform*    INTEGER Type of transform
!                   =  1 : forward
!                   = -1 : inverse
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cval
INTEGER :: k, nh, npcc, nq
REAL :: c1=0.5, c2, theta
COMPLEX :: z
COMPLEX, DIMENSION(nodat/4) :: w
COMPLEX, DIMENSION(nodat/4-1) :: h1, h2


procname(pglev+1) = 'fft_realft1'
CALL log_timer_in(npcc)

!
!1. Check array bounds
!---------------------
!

IF (IAND(nodat,nodat-1).NE.0) THEN
   nerrs = 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') nodat; cval = ADJUSTL(cval)
      WRITE (ioerr,'(A)') 'Wrong value for first data dimension: '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be an integer power of 2'
   ENDIF
   CALL error_abort('fft_realft1',ierrno_runval)
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!

nh = nodat/2; nq = nodat/4
c2 = -iform*0.5

theta = iform*twopi/REAL(nodat)
z = CMPLX(COS(theta),SIN(theta))
w(1) = CMPLX(1.0,0.0)
k_220: DO k=2,nq
   w(k) = w(k-1)*z
ENDDO k_220

!
!3. Forward transform
!--------------------
!

IF (iform.EQ.1) THEN
   spec = CMPLX(data(1:nodat-1:2),data(2:nodat:2))
   CALL fft_four1(spec,nodat/2,1)
ENDIF

!
!4. Complex data array
!---------------------
!

h1 = c1*(spec(2:nq)+CONJG(spec(nh:nq+2:-1)))
h2 = c2*(spec(2:nq)-CONJG(spec(nh:nq+2:-1)))

spec(2:nq) = h1 + w(2:nq)*h2
spec(nh:nq+2:-1) = CONJG(h1-w(2:nq)*h2)
z = spec(1)

!
!5. Apply inverse transform
!--------------------------
!

IF (iform.EQ.1) THEN
   spec(1) = CMPLX(REAL(z)+AIMAG(z),REAL(z)-AIMAG(z))
ELSE
   spec(1) = CMPLX(c1*(REAL(z)+AIMAG(z)),c1*(REAL(z)-AIMAG(z)))
   CALL fft_four1(spec,nodat/2,-1)
ENDIF

!
!6. Transform to real format
!---------------------------
!

IF (iform.EQ.-1) THEN 
   data(1:nodat-1:2) = REAL(spec)
   data(2:nodat:2) = AIMAG(spec)
ENDIF

CALL log_timer_out(npcc,itm_libs)


RETURN

END SUBROUTINE fft_realft1

!===========================================================================

SUBROUTINE fft_realft2(data,spec,speq,ndims,iform)
!**************************************************************************
!
! *fft_realft2* 2-D Fourier transform
!
! Author -
!
! Version - @(COHERENS)fft_library.f90  V2.0
!
! Description 
!  - If iform=1, returns the complex fast Fourier transform
!    of a 2-D real input "data" array as two complex arrays : "spec" with the
!    zero and positive frequency values of the first frequency component and
!    "speq" with the Nyquist critical frequency values of the first frequency
!    component. The second frequency components are stored for zero,
!    positive and negative frequencies, in standard wrap-around order.
!  - If form = -1, the inverse transform (times ndims(1)*ndims(2)/2) is
!    performed, with output "data" deriving from input "spec" and "speq".
!  - Array sizes must be integer powers of 2.
!
! Reference - Numerical recipes (routine "rlft2")
!
! Module calls - error_abort, fft_four2
!
!*****************************************************************************
!
USE error_routines, ONLY: error_abort

!
!*Arguments
!
INTEGER, INTENT(IN) :: iform
INTEGER, INTENT(IN), DIMENSION(2) :: ndims
REAL, INTENT(INOUT), DIMENSION(ndims(1),ndims(2)) :: data
COMPLEX, INTENT(INOUT), DIMENSION(ndims(2)) :: speq
COMPLEX, INTENT(INOUT), DIMENSION(ndims(1)/2,ndims(2)) :: spec

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*data*     REAL    Input data if iform = 1, output data (inverse transform)
!                   if iform = -1
!*spec*     COMPLEX Output data (forward transform) if iform = 1, input data
!                   if iform = -1
!*speq*     COMPLEX Nyquist critical frequencies
!*ndims*    INTEGER Shape of real "data" array and complex "spec" array
!                   (with first dimension divided by two)
!*iform*    INTEGER Type of transform
!                   =  1 : forward
!                   = -1 : inverse
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cval
INTEGER :: i1, j1, npcc, n2
REAL (KIND=kndrlong) :: theta
COMPLEX :: c1 = (0.5,0.0), c2, w
COMPLEX (KIND=kndrlong) :: wp, ww
COMPLEX, DIMENSION(ndims(2)) :: h1, h2


procname(pglev+1) = 'fft_realft2'
CALL log_timer_in(npcc)

!
!1. Check array bounds
!---------------------
!

IF (IAND(ndims(1),ndims(1)-1).NE.0) THEN
   nerrs = 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ndims(1); cval = ADJUSTL(cval)
      WRITE (ioerr,'(A)') 'Wrong value for first data dimension: '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be an integer power of 2'
   ENDIF
ENDIF
IF (IAND(ndims(2),ndims(2)-1).NE.0) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ndims(2); cval = ADJUSTL(cval)
      WRITE (ioerr,'(A)') 'Wrong value for second data dimension: '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be an integer power of 2'
   ENDIF
ENDIF
CALL error_abort('fft_realft2',ierrno_runval)

!
!2. Initialise parameters
!------------------------
!

n2 = ndims(2)
c2 = CMPLX(0.0,-0.5*iform)
theta = twopi_d/(iform*ndims(1))
wp = CMPLX(-2.0_kndrlong*SIN(0.5_kndrlong*theta)**2,SIN(theta),KIND=kndrlong)

!
!3. Case of forward transform
!----------------------------
!

IF (iform.EQ.1) THEN
   spec = CMPLX(data(1:ndims(1)-1:2,:),data(2:ndims(1):2,:))
   CALL fft_four2(spec,SHAPE(spec),iform)
   speq = spec(1,:)
ENDIF

!
!4. Complex data array
!---------------------
!

h1(1) = c1*(spec(1,1)+CONJG(speq(1)))
h1(2:n2) = c1*(spec(1,2:n2)+CONJG(speq(n2:2:-1)))
h2(1) = c2*(spec(1,1)-CONJG(speq(1)))
h2(2:n2) = c2*(spec(1,2:n2)-CONJG(speq(n2:2:-1)))
spec(1,:) = h1 + h2
speq(1) = CONJG(h1(1)-h2(1))
speq(n2:2:-1) = CONJG(h1(2:n2)-h2(2:n2))
ww = CMPLX(1.0_kndrlong)
i1_410: DO i1=2,ndims(1)/4+1
   j1 = ndims(1)/2-i1+2
   ww = ww*wp + ww
   w = ww
   h1(1) = c1*(spec(i1,1)+CONJG(spec(j1,1)))
   h1(2:n2) = c1*(spec(i1,2:n2)+CONJG(spec(j1,n2:2:-1)))
   h2(1) = c2*(spec(i1,1)-CONJG(spec(j1,1)))
   h2(2:n2) = c2*(spec(i1,2:n2)-CONJG(spec(j1,n2:2:-1)))
   spec(i1,:) = h1+w*h2
   spec(j1,1) = CONJG(h1(1)-w*h2(1))
   spec(j1,n2:2:-1) = CONJG(h1(2:n2)-w*h2(2:n2))
ENDDO i1_410

!
!5. Apply inverse transform and convert to real format
!-----------------------------------------------------
!

IF (iform.EQ.-1) THEN
   CALL fft_four2(spec,SHAPE(spec),iform)
   data(1:ndims(1)-1:2,:) = REAL(spec)
   data(2:ndims(1):2,:) = AIMAG(spec)
ENDIF

CALL log_timer_out(npcc,itm_libs)


RETURN

END SUBROUTINE fft_realft2

!===========================================================================

SUBROUTINE fft_realft3(data,spec,speq,ndims,iform)
!**************************************************************************
!
! *fft_realft3* 3-D Fourier transform
!
! Author -
!
! Version - @(COHERENS)fft_library.f90  V2.0
!
! Description 
!  - If iform=1, returns the complex fast Fourier transform
!    of a 3-D real input "data" array as two complex arrays : "spec" with the
!    zero and positive frequency values of the first frequency component and
!    "speq" with the Nyquist critical frequency values of the first frequency
!    component. The second and third frequency components are stored for zero,
!    positive and negative frequencies, in standard wrap-around order.
!  - If form = -1, the inverse transform (times ndims(1)*ndims(2)*ndims(3)/2)
!    is performed, with output "data" deriving from input "spec" and "speq".
!  - Array sizes must be integer powers of 2.
!
! Reference - Numerical recipes (routine "rlft3")
!
! Module calls - error_abort, fft_four3
!
!*****************************************************************************
!
USE error_routines, ONLY: error_abort

!
!*Arguments
!
INTEGER, INTENT(IN) :: iform
INTEGER, INTENT(IN), DIMENSION(3) :: ndims
REAL, INTENT(INOUT), DIMENSION(ndims(1),ndims(2),ndims(3)) :: data
COMPLEX, INTENT(INOUT), DIMENSION(ndims(2),ndims(3)) :: speq
COMPLEX, INTENT(INOUT), DIMENSION(ndims(1)/2,ndims(2),ndims(3)) :: spec

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*data*     REAL    Input data if iform = 1, output data (inverse transform)
!                   if iform = -1
!*spec*     COMPLEX Output data (forward transform) if iform = 1, input data
!                   if iform = -1
!*speq*     COMPLEX Nyquist critical frequencies
!*ndims*    INTEGER Shape of real "data" array and complex "spec" array
!                   (with first dimension divided by two)
!*iform*    INTEGER Type of transform
!                   =  1 : forward
!                   = -1 : inverse
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cval
INTEGER :: i1, i3, j1, j3, npcc, n2
REAL (KIND=kndrlong) :: theta
COMPLEX :: c1 = (0.5,0.0), c2, w
COMPLEX (KIND=kndrlong) :: wp, ww
COMPLEX, DIMENSION(ndims(2)) :: h1, h2


procname(pglev+1) = 'fft_realft3'
CALL log_timer_in(npcc)

!
!1. Check array bounds
!---------------------
!

IF (IAND(ndims(1),ndims(1)-1).NE.0) THEN
   nerrs = 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ndims(1); cval = ADJUSTL(cval)
      WRITE (ioerr,'(A)') 'Wrong value for first data dimension: '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be an integer power of 2'
   ENDIF
ENDIF
IF (IAND(ndims(2),ndims(2)-1).NE.0) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ndims(2); cval = ADJUSTL(cval)
      WRITE (ioerr,'(A)') 'Wrong value for second data dimension: '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be an integer power of 2'
   ENDIF
ENDIF
IF (IAND(ndims(3),ndims(3)-1).NE.0) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ndims(3); cval = ADJUSTL(cval)
      WRITE (ioerr,'(A)') 'Wrong value for third data dimension: '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be an integer power of 2'
   ENDIF
ENDIF
CALL error_abort('fft_realft3',ierrno_runval)

!
!2. Initialise parameters
!------------------------
!

n2 = ndims(2)
c2 = CMPLX(0.0,-0.5*iform)
theta = twopi_d/(iform*ndims(1))
wp = CMPLX(-2.0_kndrlong*SIN(0.5_kndrlong*theta)**2,SIN(theta),KIND=kndrlong)

!
!3. Case of forward transform
!----------------------------
!

IF (iform.EQ.1) THEN
   spec = CMPLX(data(1:ndims(1)-1:2,:,:),data(2:ndims(1):2,:,:))
   CALL fft_four3(spec,SHAPE(spec),iform)
   speq = spec(1,:,:)
ENDIF

!
!4. Complex data array
!---------------------
!

i3_410: DO i3=1,ndims(3)
   j3 = MERGE(1,ndims(3)-i3+2,i3.EQ.1)
   h1(1) = c1*(spec(1,1,i3)+CONJG(speq(1,j3)))
   h1(2:n2) = c1*(spec(1,2:n2,i3)+CONJG(speq(n2:2:-1,j3)))
   h2(1) = c2*(spec(1,1,i3)-CONJG(speq(1,j3)))
   h2(2:n2) = c2*(spec(1,2:n2,i3)-CONJG(speq(n2:2:-1,j3)))
   spec(1,:,i3) = h1 + h2
   speq(1,j3) = CONJG(h1(1)-h2(1))
   speq(n2:2:-1,j3) = CONJG(h1(2:n2)-h2(2:n2))
   ww = CMPLX(1.0_kndrlong,KIND=kndrlong)
   i1_410: DO i1=2,ndims(1)/4+1
      j1 = ndims(1)/2-i1+2
      ww = ww*wp + ww
      w = ww
      h1(1) = c1*(spec(i1,1,i3)+CONJG(spec(j1,1,j3)))
      h1(2:n2) = c1*(spec(i1,2:n2,i3)+CONJG(spec(j1,n2:2:-1,j3)))
      h2(1) = c2*(spec(i1,1,i3)-CONJG(spec(j1,1,j3)))
      h2(2:n2) = c2*(spec(i1,2:n2,i3)-CONJG(spec(j1,n2:2:-1,j3)))
      spec(i1,:,i3) = h1 + w*h2
      spec(j1,1,j3) = CONJG(h1(1)-w*h2(1))
      spec(j1,n2:2:-1,j3) = CONJG(h1(2:n2)-w*h2(2:n2))
   ENDDO i1_410
ENDDO i3_410

!
!5. Apply inverse transform and convert to real format
!-----------------------------------------------------
!

IF (iform.EQ.-1) THEN
   CALL fft_four3(spec,SHAPE(spec),iform)
   data(1:ndims(1)-1:2,:,:) = REAL(spec)
   data(2:ndims(1):2,:,:) = AIMAG(spec)
ENDIF

CALL log_timer_out(npcc,itm_libs)


RETURN

END SUBROUTINE fft_realft3


END MODULE fft_library
