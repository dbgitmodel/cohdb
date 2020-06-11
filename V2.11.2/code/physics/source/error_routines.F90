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

MODULE error_routines
!************************************************************************
!
! *error_routines* Error reporting
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.11
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
! Generic routines - error_arg_var, error_lbound_arr, error_lbound_arr_struc,
!    error_lbound_var, error_limits_arr, error_limits_arr_struc,
!    error_limits_var, error_ubound_arr, error_ubound_arr_struc,
!    error_ubound_var, error_value_arr, error_value_arr_struc, error_value_var,
!    warning_reset_arr, warning_reset_arr_struc, warning_reset_var
!
! Routines - check_mult_counter, check_space_limits,
!    check_space_limits_arr_struc, check_time_limits,
!    check_time_limits_arr_struc, error_abort, error_alloc, error_alloc_struc,
!    error_alloc_struc_comp, error_array_index, error_diff_vals_arrlist,
!    error_diff_vals_varlist, error_dim_arr, error_file, error_mult,
!    error_proc, error_shape, error_vals_arr_char, error_vals_arr_int,
!    error_vals_arr_struc_char, error_vals_arr_struc_int, error_vals_var_char,
!    error_vals_var_int, warning_ubound_var_real, warning_ubound_var_double,
!    write_error_message
!
!************************************************************************
!
USE iopars
USE syspars

IMPLICIT NONE

INTERFACE error_arg_var
   MODULE PROCEDURE error_arg_var_char, error_arg_var_int, &
                  & error_arg_var_log
END INTERFACE

INTERFACE error_lbound_arr
   MODULE PROCEDURE error_lbound_arr_int, error_lbound_arr_real, &
                  & error_lbound_arr_double
END INTERFACE

INTERFACE error_lbound_arr_struc
   MODULE PROCEDURE error_lbound_arr_struc_int, error_lbound_arr_struc_real, &
                  & error_lbound_arr_struc_double
END INTERFACE

INTERFACE error_lbound_var
   MODULE PROCEDURE error_lbound_var_int, error_lbound_var_real, &
                  & error_lbound_var_double
END INTERFACE

INTERFACE error_limits_arr
   MODULE PROCEDURE error_limits_arr_int, error_limits_arr_real, &
                  & error_limits_arr_double
END INTERFACE

INTERFACE error_limits_arr_struc
   MODULE PROCEDURE error_limits_arr_struc_int, error_limits_arr_struc_real, &
                  & error_limits_arr_struc_double
END INTERFACE

INTERFACE error_limits_var
   MODULE PROCEDURE error_limits_var_int, error_limits_var_real, &
                  & error_limits_var_double
END INTERFACE

INTERFACE error_ubound_arr
   MODULE PROCEDURE error_ubound_arr_int, error_ubound_arr_real, &
                  & error_ubound_arr_double
END INTERFACE

INTERFACE error_ubound_arr_struc
   MODULE PROCEDURE error_ubound_arr_struc_int, error_ubound_arr_struc_real, &
                  & error_ubound_arr_struc_double
END INTERFACE

INTERFACE error_ubound_var
   MODULE PROCEDURE error_ubound_var_int, error_ubound_var_real, &
                  & error_ubound_var_double
END INTERFACE

INTERFACE error_value_arr
   MODULE PROCEDURE error_value_arr_int, error_value_arr_log
END INTERFACE

INTERFACE error_value_arr_struc
   MODULE PROCEDURE error_value_arr_struc_int, error_value_arr_struc_log
END INTERFACE

INTERFACE error_value_var
   MODULE PROCEDURE error_value_var_char, error_value_var_int, &
                  & error_value_var_log, error_value_var_real, &
                  & error_value_var_double
END INTERFACE

INTERFACE warning_reset_arr
   MODULE PROCEDURE warning_reset_arr_char, warning_reset_arr_int, &
                  & warning_reset_arr_log
END INTERFACE

INTERFACE warning_reset_arr_struc
   MODULE PROCEDURE warning_reset_arr_struc_char, warning_reset_arr_struc_int, &
                  & warning_reset_arr_struc_log
END INTERFACE

INTERFACE warning_reset_var
   MODULE PROCEDURE warning_reset_var_char, warning_reset_var_int, &
                  & warning_reset_var_log, warning_reset_var_real, &
                  & warning_reset_var_double
END INTERFACE

CONTAINS

!========================================================================

SUBROUTINE check_mult_counter(ival,varname,ifix,matchzero)
!************************************************************************
!
! *check_mult_counter* Check whether counter ival is an integer multiple of ifix
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.8
!
! Description - displays error message if needed
!
! Module calls - mult_index
!
!************************************************************************
!
USE utility_routines, ONLY: mult_index

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
LOGICAL, INTENT(IN), OPTIONAL :: matchzero
INTEGER, INTENT(IN) :: ifix, ival

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Start/end/step spatial indices
!*varname*   CHAR    Counter name
!*ifix*      INTEGER Counter value of which ival must be a multiple
!*matchzero* LOGICAL If PRESENT and .TRUE., result is .TRUE. if 'iy' equals 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=12) :: cfix, cval


IF (PRESENT(matchzero)) THEN
   flag = matchzero
ELSE
   flag = .FALSE.
ENDIF
   
IF (.NOT.mult_index(ival,ifix,matchzero=flag)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; WRITE (cfix,'(I12)') ifix
      cval = ADJUSTL(cval); cfix = ADJUSTL(cfix)
      WRITE (ioerr,'(A)') 'Invalid value for counter '//TRIM(varname)//': '&
                        & //TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be an integer multiple of : '//TRIM(cfix)
   ENDIF
ENDIF


RETURN

END SUBROUTINE check_mult_counter

!========================================================================

SUBROUTINE check_space_limits(slims,arrname,limmax)
!************************************************************************
!
! *check_space_limits* Check Start/End/Step values of the space-limits array 
!                      'slims'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays warning message if needed
!
! Module calls - error_limits_arr, mult_index
!
!************************************************************************
!
USE utility_routines, ONLY: mult_index

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: limmax
INTEGER, INTENT(IN), DIMENSION(3) :: slims

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*slims*    INTEGER Start/end/step spatial indices
!*arrname*  CHAR    Array name
!*limmax*   INTEGER Maximum allowed end value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(3) :: clims


!---end value
CALL error_limits_arr(slims(2),arrname,1,limmax,1,indx=(/2/))
!---start value
CALL error_limits_arr(slims(1),arrname,1,slims(2),1,indx=(/1/))
!---step value
IF (.NOT.mult_index(slims(2)-slims(1),slims(3))) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      n_110: DO n=1,3
         WRITE (clims(n),'(I12)') slims(n); clims(n) = ADJUSTL(clims(n))
      ENDDO n_110
      WRITE (ioerr,8001) 'Invalid values for space-limit array '//&
                       & TRIM(arrname), (TRIM(clims(n)),n=1,3)
      WRITE (ioerr,'(A)') '  Interval must contain an integer number of&
                           & space steps'
   ENDIF
ENDIF


RETURN

8001 FORMAT(A,':',3(' ',A))

END SUBROUTINE check_space_limits

!========================================================================

SUBROUTINE check_space_limits_arr_struc(slims,arrname,compname,limmax,ndims,&
                                      & indx)
!************************************************************************
!
! *check_space_limits_arr_struc* Check Start/End/Step values of space-limits
!                                component 'slims' of a structure array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays warning message if needed
!
! Module calls - error_limits_arr_struc, mult_index
!
!************************************************************************
!
USE utility_routines, ONLY: mult_index

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
INTEGER, INTENT(IN) :: limmax, ndims
INTEGER, INTENT(IN), DIMENSION(3) :: slims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*slims*    INTEGER Start/end/step spatial indices
!*arrname*  CHAR    Name of derived type array
!*compname* CHAR    Name of derive type component
!*limmax*   INTEGER Maximum allowed end value
!*ndims*    INTEGER Array rank
!*indx*     INTEGER Array vector index
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(3) :: clims
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


!---end value
CALL error_limits_arr_struc(slims(2),arrname,compname,1,limmax,1,indx=(/2/))

!---start value
CALL error_limits_arr_struc(slims(1),arrname,compname,1,slims(2),1,indx=(/1/))

!---step value
IF (.NOT.mult_index(slims(2)-slims(1),slims(3))) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      n_110: DO n=1,3
         WRITE (clims(n),'(I12)') slims(n); clims(n) = ADJUSTL(clims(n))
      ENDDO n_110
      n_120: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_120
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,8001) 'Invalid value for space-limit component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array '//TRIM(arrname),(TRIM(clims(n)),n=1,3)
      CASE (2)
         WRITE (ioerr,8001) 'Invalid value for space-limit component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array '//TRIM(arrname),&
                 & (TRIM(clims(n)),n=1,3)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for space-limit component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '&
                 & //TRIM(arrname), (TRIM(clims(n)),n=1,3)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for space-limit component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//TRIM(cindx(4))//&
                 & ') of array '//TRIM(arrname), (TRIM(clims(n)),n=1,3)
      END SELECT
      WRITE (ioerr,'(A)') 'Interval must contain an integer number of&
                           & space steps'
   ENDIF
ENDIF


RETURN

8001 FORMAT(A,':',3(' ',A))

END SUBROUTINE check_space_limits_arr_struc

!========================================================================

SUBROUTINE check_time_limits(tlims,arrname,minstep,maxstep)
!************************************************************************
!
! *check_time_limits* Check Start/End/Step values of time-limits array 'tlims'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays warning message if needed
!
! Module calls - error_lbound_arr, error_limits_arr, error_ubound_arr,
!                mult_index, warning_reset_arr
!
!************************************************************************
!
USE utility_routines, ONLY: mult_index

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: maxstep, minstep
INTEGER, INTENT(INOUT), DIMENSION(3) :: tlims

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*tlims*     INTEGER Start/end/step time indices
!*arrname*   CHAR    Array name
!*minstep*   INTEGER Minimum start value
!*maxstep*   INTEGER Maximum end value
!*end_reset* LOGICAL Reset end value to its maximum if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n
CHARACTER (LEN=22), DIMENSION(3) :: clims


!---end value
CALL error_lbound_arr(tlims(2),arrname,minstep,.TRUE.,1,indx=(/2/))

!---start value
CALL error_limits_arr(tlims(1),arrname,minstep,tlims(2),1,indx=(/1/))

!---step value
IF (.NOT.mult_index(tlims(2)-tlims(1),tlims(3))) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      n_110: DO n=1,3
         WRITE (clims(n),'(I12)') tlims(n); clims(n) = ADJUSTL(clims(n))
      ENDDO n_110
      WRITE (ioerr,8001) 'Invalid values for time-limit array '//TRIM(arrname),&
                        & (TRIM(clims(n)),n=1,3)
      WRITE (ioerr,'(A)') '  Interval must contain an integer number of time&
                          & steps'
   ENDIF
ENDIF
!---reset end value if necessary
IF (maxstep.GT.0.AND.tlims(2).GT.maxstep) THEN
   CALL error_ubound_arr(tlims(2),arrname,maxstep,.TRUE.,1,(/2/))
ENDIF


RETURN

8001 FORMAT(A,':',3(' ',A))

END SUBROUTINE check_time_limits

!========================================================================

SUBROUTINE check_time_limits_arr_struc(tlims,arrname,compname,minstep,maxstep,&
                                     & ndims,indx,flag3d)
!************************************************************************
!
! *check_time_limits_arr_struc* Check Start/End/Step values of time-limits
!                               component 'tlims' of a structure array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.7
!
! Description - displays warning message if needed
!
! Module calls - error_lbound_arr_struc, error_limits_arr_struc,
!                error_ubound_arr_struc, mult_index, warning_reset_arr_struc
!
!************************************************************************
!
USE timepars
USE utility_routines, ONLY: mult_index

!
!*Arguments
!
LOGICAL, INTENT(IN) , OPTIONAL :: flag3d
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
INTEGER, INTENT(IN) :: maxstep, minstep, ndims
INTEGER, INTENT(INOUT), DIMENSION(3) :: tlims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*tlims*     INTEGER Start/end/step time indices
!*arrname*   CHAR    Name of derived type array
!*compname*  CHAR    Name of derived type component
!*minstep*   INTEGER Minimum start value
!*maxstep*   INTEGER Maximum end value
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
!*flag3d*    LOGICAL Tests whether a time interval contains an integer number of
!                    3-D time steps
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag1, flag2, flag3dx
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(3) :: clims
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


!---optional argument
IF (PRESENT(flag3d)) THEN
   flag3dx = flag3d
ELSE
   flag3dx = .FALSE.
ENDIF

!---end value
CALL error_lbound_arr_struc(tlims(2),arrname,'%tlims(2)',minstep,.TRUE.,1,&
                          & indx)

!---start value
CALL error_limits_arr_struc(tlims(1),arrname,'%tlims(1)',minstep,tlims(2),1,&
                          & indx)

!---step value
flag1 = .FALSE.; flag2 = .FALSE.
IF (tlims(3).NE.0) THEN
   IF (.NOT.mult_index(tlims(2)-tlims(1),tlims(3))) THEN
      flag1 = .TRUE.
      nerrs = nerrs + 1
   ENDIF
   IF (.NOT.mult_index(tlims(3),ic3d).AND.flag3dx) THEN
      flag2 = .TRUE.
      nerrs = nerrs + 1
   ENDIF
   IF (errchk.AND.(flag1.OR.flag2).AND.nerrs.LE.maxerrors) THEN
      n_110: DO n=1,3
         WRITE (clims(n),'(I12)') tlims(n); clims(n) = ADJUSTL(clims(n))
      ENDDO n_110
      n_120: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_120
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,8001) 'Invalid value for time-limit component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array '//TRIM(arrname), (TRIM(clims(n)),n=1,3)
      CASE (2)
         WRITE (ioerr,8001) 'Invalid value for time-limit component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array '//TRIM(arrname),&
                 & (TRIM(clims(n)),n=1,3)
      CASE (3)
         WRITE (ioerr,8001) 'Invalid value for time-limit component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '&
                 & //TRIM(arrname), (TRIM(clims(n)),n=1,3)
      CASE DEFAULT
         WRITE (ioerr,8001) 'Invalid value for time-limit component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//TRIM(cindx(4))//&
                 & ') of array '//TRIM(arrname), (TRIM(clims(n)),n=1,3)
      END SELECT
      IF (flag1) THEN
         WRITE (ioerr,'(A)') 'Interval must contain an integer number of '//&
                           & 'time steps'
      ENDIF
      IF (flag2) THEN
         WRITE (ioerr,'(A)') 'Interval must contain an integer number of '//&
                           & '3-D time steps'
      ENDIF
   ENDIF
ENDIF

!---reset end value if necessary
IF (maxstep.GT.0.AND.tlims(2).GT.maxstep) THEN
   CALL error_ubound_arr_struc(tlims(2),arrname,compname,maxstep,.TRUE.,1,(/2/))
ENDIF


RETURN

8001 FORMAT(A,':',3(' ',A))

END SUBROUTINE check_time_limits_arr_struc

!========================================================================

SUBROUTINE error_abort(prname,icode)
!************************************************************************
!
! *error_abort* Display error message and abort program if necessary
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description -
!
! MPI calls - MPI_abort
!
!************************************************************************
!
USE paralpars
USE switches
#ifdef MPI
   INCLUDE 'mpif.h'
#endif /*MPI*/

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: prname
INTEGER, INTENT(IN) :: icode

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*prname*    CHAR    Name of routine where error occurred
!*icode*     INTEGER Error code number
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: ccode, cnerrs
#ifdef MPI
   INTEGER :: ierr
#endif /*MPI*/


IF (nerrs.GT.0) THEN
   IF (loglev1.GT.0) THEN
      WRITE (ccode,'(I12)') icode; ccode = ADJUSTL(ccode)
      WRITE (cnerrs,'(I12)') nerrs; cnerrs = ADJUSTL(cnerrs)
      WRITE (iolog,9001) TRIM(cnerrs), prname
      WRITE (iolog,'(A)') 'PROGRAM TERMINATED ABNORMALLY'
   ENDIF
   IF (errchk) THEN
      WRITE (ccode,'(I12)') icode; ccode = ADJUSTL(ccode)
      WRITE (cnerrs,'(I12)') nerrs; cnerrs = ADJUSTL(cnerrs)
      WRITE (ioerr,9001) TRIM(cnerrs), prname
      WRITE (ioerr,8001) TRIM(ccode), TRIM(error_code(icode))
      WRITE (ioerr,'(A)') 'PROGRAM TERMINATED ABNORMALLY'
   ENDIF
   IF (iopt_MPI.EQ.0) THEN
      STOP
   ELSE
#ifdef MPI
      CALL MPI_abort(comm_world_MPI,icode,ierr)
#endif /*MPI*/
   ENDIF
ENDIF


RETURN

9001 FORMAT ('A total of ',A,' errors occurred in ',A)
8001 FORMAT ('Error type ',A,' : ',A)

END SUBROUTINE error_abort

!========================================================================

SUBROUTINE error_alloc(arrname,ndims,nshape,kind_type,lenstr,abort)
!************************************************************************
!
! *error_alloc* Displays an error message when a dynamic array cannot be
!               allocated
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.10.1
!
! Description - abort program if abort is not present or present and .TRUE.
!
! Module calls - error_abort
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: abort
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: kind_type, ndims
INTEGER, INTENT(IN), OPTIONAL :: lenstr
INTEGER, INTENT(IN), DIMENSION(ndims) :: nshape

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*arrname*   CHAR    Array name
!*ndims*     INTEGER Array rank
!*nshape*    INTEGER Array shape
!*kind_type* INTEGER KIND type of the data
!*lenstr*    INTEGER String length in case of character array
!*abort*     LOGICAL Abort program if PRESENT and .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cerrstat
CHARACTER (LEN=22) :: csize
INTEGER :: n
INTEGER (KIND=kndilong) :: nsize
CHARACTER (LEN=12), DIMENSION(ndims) :: cshape


IF (errstat.NE.0) THEN

   nsize = PRODUCT(nshape)
   IF (kind_type.EQ.kndchar) THEN
      nsize = kndchar*lenstr*nsize
   ELSE
      nsize = kind_type*nsize
   ENDIF

   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (ioerr,'(2A)') 'Unable to allocate array: ', arrname
      n_110: DO n=1,ndims
         WRITE (cshape(n),'(I12)') nshape(n)
         cshape(n) = ADJUSTL(cshape(n))
      ENDDO n_110
      WRITE (ioerr,'(8A)') 'Shape:', (' '//TRIM(cshape(n)),n=1,ndims)
      WRITE (csize,'(I12)') nsize
      csize = ADJUSTL(csize)
      WRITE (ioerr,'(3A)') 'Size: ', TRIM(csize), ' bytes'
      WRITE (cerrstat,'(I12)') errstat
      cerrstat = ADJUSTL(cerrstat)
      WRITE (ioerr,'(A)') 'Exit code '//TRIM(cerrstat)
   ENDIF
   
   IF (.NOT.PRESENT(abort)) THEN
      CALL error_abort(procname(pglev),ierrno_alloc)
   ELSEIF (abort) THEN
      CALL error_abort(procname(pglev),ierrno_alloc)
   ENDIF

ENDIF


RETURN

END SUBROUTINE error_alloc

!========================================================================

SUBROUTINE error_alloc_struc(arrname,ndims,nshape,structype,abort)
!************************************************************************
!
! *error_alloc_struc* Displays an error message when a dynamic array of derived
!                     type cannot be allocated
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.11
!
! Description -
!
! Module calls - error_abort, error_arg_var
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: abort
CHARACTER (LEN=*), INTENT(IN) :: arrname, structype
INTEGER, INTENT(IN) :: ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: nshape

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*arrname*   CHAR    Array name
!*ndims*     INTEGER Array rank
!*nshape*    INTEGER Array shape
!*structype* CHAR    Type of derived structure
!*lenstr*    INTEGER String length in case of character array
!*abort*     LOGICAL Abort program if PRESENT and .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!

CHARACTER (LEN=12) :: cerrstat
CHARACTER (LEN=22) :: csize
INTEGER :: n
INTEGER (KIND=kndilong) :: nsize
CHARACTER (LEN=12), DIMENSION(ndims) :: cshape


IF (errstat.NE.0) THEN

   SELECT CASE (TRIM(structype))
      CASE ('ExchComms')
         nsize = kndlog+kndint*19
      CASE ('FileParams')
         nsize = kndlog*6+kndchar*(2+lenform+lencomm+leniofile*2+lendesc)+&
               & kndint*13+kndrtype*3
      CASE ('GlbAtts')
         nsize = kndchar*lendesc*2
      CASE ('GridParams')
         nsize = kndlog*3+kndint*3+kndrtype*8
      CASE ('HRelativeCoords')
         nsize = kndint*2+kndrtype*4
      CASE ('OutGridParams')
         nsize = kndlog*6+kndchar*(1+lentime*3)+kndint*20+kndrtype
      CASE ('PartAtts')
         nsize = kndlog + kndint*8 + kndrtype*5 + kndrlong*5
      CASE ('OutPartParams')
         nsize = kndlog + kndchar*lentime*3 + kndint*5 + kndrtype*4
      CASE ('PartCloud')
         nsize = kndint*4 + kndrtype*4 
      CASE ('StationLocs')
         nsize = kndint*2+kndchar*lenname
      CASE ('VariableAtts')
         nsize = kndlog+kndchar*(1+lenname+lendesc*5+lenunit+lennode)+&
               & kndint*22+kndrtype+kndrlong
      CASE ('VRelativeCoords')
         nsize = kndint+kndrtype*2
      CASE DEFAULT
         procname(pglev+1) = 'error_alloc_struc'
         CALL error_arg_var(structype,'structype')
   END SELECT
   nsize = nsize*PRODUCT(nshape)

   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (ioerr,'(2A)') 'Unable to allocate derived type array: ', arrname
      n_110: DO n=1,ndims
         WRITE (cshape(n),'(I12)') nshape(n)
         cshape(n) = ADJUSTL(cshape(n))
      ENDDO n_110
      WRITE (ioerr,'(8A)') 'Shape:', (' '//TRIM(cshape(n)),n=1,ndims)
      WRITE (csize,'(I12)') nsize
      csize = ADJUSTL(csize)
      WRITE (ioerr,'(3A)') 'Size: ', TRIM(csize), ' bytes'
      WRITE (cerrstat,'(I12)') errstat
      cerrstat = ADJUSTL(cerrstat)
      WRITE (ioerr,'(A)') 'Exit code '//TRIM(cerrstat)
   ENDIF
   
   IF (.NOT.PRESENT(abort)) THEN
      CALL error_abort(procname(pglev),ierrno_alloc)
   ELSEIF (abort) THEN
      CALL error_abort(procname(pglev),ierrno_alloc)
   ENDIF

ENDIF


RETURN

END SUBROUTINE error_alloc_struc

!========================================================================

SUBROUTINE error_alloc_struc_comp(strucname,arrname,ndims,nshape,kind_type,&
                                & indx,lenstr,abort)
!************************************************************************
!
! *error_alloc_struc_comp* Displays an error message when an array component of
!                          a derived type scalar or vector cannot be allocated
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.10.1
!
! Description -
!
! Module calls - error_abort, error_arg_var
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: abort
CHARACTER (LEN=*), INTENT(IN) :: arrname, strucname
INTEGER, INTENT(IN) :: kind_type, ndims
INTEGER, INTENT(IN), OPTIONAL :: indx, lenstr
INTEGER, INTENT(IN), DIMENSION(ndims) :: nshape

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*strucname* CHAR    Name of the derived type scalar or vector array
!*arrname*   CHAR    Name of the derive type component
!*ndims*     INTEGER Rank of the component array
!*nshape*    INTEGER Shape of the component array
!*kind_type* INTEGER KIND type of the data
!*indx*      INTEGER If PRESENT, index of the derived type vector
!*abort*     LOGICAL Abort program if PRESENT and .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cerrstat, cindx
CHARACTER (LEN=22) :: csize
INTEGER :: n
INTEGER (KIND=kndilong) :: nsize
CHARACTER (LEN=12), DIMENSION(ndims) :: cshape

IF (errstat.NE.0) THEN

   nsize = PRODUCT(nshape)
   IF (kind_type.EQ.kndchar) THEN
      nsize = kndchar*lenstr*nsize
   ELSE
      nsize = kind_type*nsize
   ENDIF

   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      IF (.NOT.PRESENT(indx)) THEN
         WRITE (ioerr,'(2A)') 'Unable to allocate array: ', &
                            & strucname//'%'//arrname
      ELSE
         WRITE (cindx,'(I12)') indx; cindx = ADJUSTL(cindx)
         WRITE (ioerr,'(2A)') 'Unable to allocate array: ', &
                             & strucname//'('//TRIM(cindx)//')'//'%'//arrname
      ENDIF
      n_110: DO n=1,ndims
         WRITE (cshape(n),'(I12)') nshape(n)
         cshape(n) = ADJUSTL(cshape(n))
      ENDDO n_110
      WRITE (ioerr,'(8A)') 'Shape:', (' '//TRIM(cshape(n)),n=1,ndims)
      WRITE (csize,'(I12)') nsize
      csize = ADJUSTL(csize)
      WRITE (ioerr,'(3A)') 'Size: ', TRIM(csize), ' bytes'
      WRITE (cerrstat,'(I12)') errstat
      cerrstat = ADJUSTL(cerrstat)
      WRITE (ioerr,'(A)') 'Exit code '//TRIM(cerrstat)
   ENDIF
   
   IF (.NOT.PRESENT(abort)) THEN
      CALL error_abort(procname(pglev),ierrno_alloc)
   ELSEIF (abort) THEN
      CALL error_abort(procname(pglev),ierrno_alloc)
   ENDIF

ENDIF


RETURN

END SUBROUTINE error_alloc_struc_comp

!========================================================================

SUBROUTINE error_array_index(ival,arrname,minval,maxval,ndim)
!************************************************************************
!
! *error_array_index* Checks whether the value the index ival belonging to
!                     dimension 'ndim' of an array is between the bounds
!                     'minval' and 'maxval'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.1.0
!
! Description - displays error message if needed
!
! Module calls - error_abort
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: ival, minval, maxval, ndim

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ival*    INTEGER Array index
!*arrname* CHAR    Array name
!*minval*  INTEGER Lower array bound
!*maxval*  INTEGER Upper array bound
!*ndim*    INTEGER Dimensio of array index
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cdim, cmaxval, cminval, cval


IF (ival.LT.minval.OR.ival.GT.maxval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      SELECT CASE (ndim)
         CASE(1); cdim = 'first'
         CASE(2); cdim = 'second'
         CASE(3); cdim = 'third'
         CASE(4); cdim = 'fourth'
      END SELECT
      WRITE (cval,'(I12)') ival; cval = ADJUSTL(cval)
      WRITE (cminval,'(I12)') minval; cminval = ADJUSTL(cminval)
      WRITE (cmaxval,'(I12)') maxval; cmaxval = ADJUSTL(cmaxval)
      WRITE (ioerr,'(A)') 'Invalid value for the '//TRIM(cdim)//' array index '&
                        & //TRIM(cval)//' of array '//TRIM(arrname)
      WRITE (ioerr,8001) 'Must be between: ', TRIM(cminval), TRIM(cmaxval)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_inival)
ENDIF


RETURN

8001 FORMAT (A,A,' and ',A)

END SUBROUTINE error_array_index

!========================================================================

SUBROUTINE error_arg_var_char(cval,argname,cfix)
!************************************************************************
!
! *error_arg_var_char* Checks whether the value 'cval' of a string argument
!                      in a  routine call equals the specified string 'cfix'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!             - if cfix is not present, cval is taken as an invalid string
!               argument
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: argname, cval
CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: cfix

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cval*      CHAR    Input string vaue
!*argname*   CHAR    Name of argument
!*cfix*      CHAR    Comparison string value
!
!------------------------------------------------------------------------------
!

IF (PRESENT(cfix)) THEN
   IF (TRIM(cval).NE.TRIM(cfix)) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (ioerr,'(A)') 'Invalid return/input value for string argument '&
                           & //TRIM(argname)//' in call of routine '&
                           & //TRIM(procname(pglev+1))//': '//TRIM(cval)
         WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cfix)
      ENDIF
   ENDIF
ELSE
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (ioerr,'(A)') 'Invalid return/input value for string argument '&
                        & //TRIM(argname)//' in call of routine '&
                        & //TRIM(procname(pglev+1))//': '//TRIM(cval)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_arg_var_char

!========================================================================

SUBROUTINE error_arg_var_int(ival,argname,ifix)
!************************************************************************
!
! *error_arg_var_int* Checks whether the value of an integer argument 'ival' in
!                     a routine call equals 'ifix'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: argname
INTEGER, INTENT(IN) :: ival, ifix

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Input string value
!*argname*   CHAR    Name of argument
!*ifix*      INTEGER Comparison integer value
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cfix, cval


IF (ival.NE.ifix) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; WRITE(cfix,'(I12)') ifix
      cval = ADJUSTL(cval); cfix = ADJUSTL(cfix)
      WRITE (ioerr,'(A)') 'Invalid return/input value for integer argument '&
                        & //TRIM(argname)//' in call of routine '&
                        & //TRIM(procname(pglev+1))//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cfix)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_arg_var_int

!========================================================================

SUBROUTINE error_arg_var_log(lval,argname,lfix)
!************************************************************************
!
! *error_arg_var_log* Checks whether the value of a logical argument 'lval' in
!                     a routine call equals 'lfix'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: lval, lfix
CHARACTER (LEN=*), INTENT(IN) :: argname

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*lval*      LOGICAL Input logical variable
!*argname*   CHAR    Name of argument
!*lfix*      LOGICAL Comparison logical variable
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=7) :: cfix, cval


IF (.NOT.(lval.EQV.lfix)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      cval = MERGE('.TRUE. ','.FALSE.',lval)
      cfix = MERGE('.TRUE. ','.FALSE.',lfix)
      WRITE (ioerr,'(A)') 'Invalid return/input value for logical argument '&
                        & //TRIM(argname)//' in call of routine '&
                        & //TRIM(procname(pglev+1))//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cfix)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_arg_var_log

!========================================================================

SUBROUTINE error_diff_vals_arrlist(ilist,arrname,ndims,idim,indx,nozero)
!************************************************************************
!
! *error_diff_vals_arrlist* Check whether all elements of the vector section
!                           'ilist' within an array are all different
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - if 'nozero' is present and .TRUE., zero values are excluded
!
! Module calls - diff_vals, digit_number_int
!
!************************************************************************
!
USE utility_routines, ONLY: diff_vals, digit_number_int
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: nozero
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: idim, ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx
INTEGER, INTENT(IN), DIMENSION(:) :: ilist

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ilist*     INTEGER Array section of specified dimension
!*arrname*   CHAR    Array name
!*ndims*     INTEGER Array rank
!*idim*      INTEGER Array dimension to be checked
!*indx*      INTEGER Array index vector (specified dimension is discarded)
!*nozero*    LOGICAL Zeros are discarded if present and .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: nozerox
CHARACTER (LEN=12) :: cdec, cdim
CHARACTER (LEN=29) :: cfmt
INTEGER :: i, ii, n, ndec, ndim, nmax
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx
INTEGER, DIMENSION(SIZE(ilist)) :: ilistx


!---optional argument
IF (PRESENT(nozero)) THEN
   nozerox = nozero
ELSE
   nozerox = .FALSE.
ENDIF

!---subtract sub-list without zeros (if needed)
ndim = SIZE(ilist)
IF (nozerox) THEN
   nmax = COUNT(ilist.NE.0)
   ii = 0
   i_110: DO i=1,ndim
      IF (ilist(i).NE.0) THEN
         ii = ii + 1
         ilistx(ii) = ilist(i)
      ENDIF
   ENDDO i_110
ELSE
   nmax = ndim
   ilistx = ilist
ENDIF

!---check
IF (.NOT.diff_vals(ilistx(1:nmax))) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      n_120: DO n=1,ndims
         IF (n.EQ.idim) THEN
            cindx(n) = ':'
         ELSE
            WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
         ENDIF
      ENDDO n_120
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid values for section '//TRIM(cindx(1))//&
                 & ' of integer array '//TRIM(arrname)//':'
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid values for section ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of integer array '//&
                 & TRIM(arrname)//':'
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid values for section ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of integer array '//TRIM(arrname)//':'
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid values for section ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//','//&
                 & TRIM(cindx(4))//') of integer array '//TRIM(arrname)//':'
      END SELECT
      ndec = 1
      n_130: DO n=1,ndim
         ndec = MAX(ndec,digit_number_int(ilist(n)))
      ENDDO n_130
      ndec = ndec + 1
      WRITE (cdec,'(I12)') ndec; cdec = ADJUSTL(cdec)
      WRITE (cdim,'(I12)') ndim; cdim = ADJUSTL(cdim)
      cfmt = '('//TRIM(cdim)//'I'//TRIM(cdec)//')'
      WRITE (ioerr,cfmt) ilist
      IF (nozerox) THEN
         WRITE (ioerr,'(A)') 'All values must be different or zero'
      ELSE
         WRITE (ioerr,'(A)') 'All values must be different'
      ENDIF
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_diff_vals_arrlist

!========================================================================

SUBROUTINE error_diff_vals_varlist(ilist,listname,nozero)
!************************************************************************
!
! *error_diff_vals_varlist* Check whether the elements of the vector 'ilist'
!                           are all different
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - if nozero is present and .TRUE., zero values are excluded
!
! Module calls - diff_vals, digit_number_int
!
!************************************************************************
!
USE utility_routines, ONLY: diff_vals, digit_number_int
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: nozero
CHARACTER (LEN=*), INTENT(IN) :: listname
INTEGER, INTENT(IN), DIMENSION(:) :: ilist

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ilist*     INTEGER Array vector
!*listname*  CHAR    Vector name
!*nozero*    LOGICAL Zeros are discarded if present and .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: nozerox
CHARACTER (LEN=12) :: cdec, cdim
CHARACTER (LEN=29) :: cfmt
INTEGER :: i, ii, n, ndec, ndim, nmax
INTEGER, DIMENSION(SIZE(ilist)) :: ilistx


!---optional argument
IF (PRESENT(nozero)) THEN
   nozerox = nozero
ELSE
   nozerox = .FALSE.
ENDIF

!---subtract sub-list without zeros (if needed)
ndim = SIZE(ilist)
IF (nozerox) THEN
   nmax = COUNT(ilist.NE.0)
   ii = 0
   i_110: DO i=1,ndim
      IF (ilist(i).NE.0) THEN
         ii = ii + 1
         ilistx(ii) = ilist(i)
      ENDIF
   ENDDO i_110
ELSE
   nmax = ndim
   ilistx = ilist
ENDIF

!---check
IF (.NOT.diff_vals(ilistx(1:nmax))) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (ioerr,'(A)') 'Invalid values for list '//TRIM(listname)//':'
      ndec = 1
      n_120: DO n=1,ndim
         ndec = MAX(ndec,digit_number_int(ilist(n)))
      ENDDO n_120
      ndec = ndec + 1
      WRITE (cdec,'(I12)') ndec; cdec = ADJUSTL(cdec)
      WRITE (cdim,'(I12)') ndim; cdim = ADJUSTL(cdim)
      cfmt = '('//TRIM(cdim)//'I'//TRIM(cdec)//')'
      WRITE (ioerr,cfmt) ilist
      IF (nozerox) THEN
         WRITE (ioerr,'(A)') 'All values must be different or zero'
      ELSE
         WRITE (ioerr,'(A)') 'All values must be different'
      ENDIF
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_diff_vals_varlist

!========================================================================

SUBROUTINE error_dim_arr(nsize,arrname,nfix)
!************************************************************************
!
! *error_dim_arr* Checks whether an array has the size given by 'nfix' 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description -
!
! Module calls -
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: nfix, nsize

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*nsize*     INTEGER Array size
!*arrname*   CHAR    Array name
!*nfix*      INTEGER Value to which nrank must be equal
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cfix, csize


IF (nsize.NE.nfix) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (csize,'(I12)') nsize; WRITE(cfix,'(I12)') nfix
      csize = ADJUSTL(csize); cfix = ADJUSTL(cfix)
      WRITE (ioerr,'(A)') 'Invalid size of variable '//TRIM(arrname)//': '&
                        & //TRIM(csize)
      WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cfix)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_dim_arr

!========================================================================

SUBROUTINE error_file(ierr,iounit,filepars,abort)
!************************************************************************
!
! *error_file* Report I/O file error 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description -
!
! Module calls - error_abort
!
!********************************************************************7****
!
USE datatypes

!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: abort
INTEGER, INTENT(IN) :: ierr
INTEGER, INTENT(IN), OPTIONAL :: iounit
TYPE(FileParams), INTENT(IN), OPTIONAL :: filepars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ierr*      INTEGER File error code
!*iounit*    INTEGER File unit number (present if filepars is not present)
!*filepars*  DERIVED File attributes  (present if iounit is not present)
!*abort*     LOGICAL Abort program if .TRUE.
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cunit
INTEGER :: iunit


IF (nerrs.EQ.0) nerrs = 1
IF (errchk.AND.nerrs.LE.maxerrors) THEN
   IF (PRESENT(filepars)) THEN
      WRITE (ioerr,'(A)') 'File: ', TRIM(filepars%filename)
      iunit = filepars%iunit
   ELSE
      iunit = iounit
   ENDIF
   WRITE (cunit,'(I12)') iunit; cunit = ADJUSTL(cunit)
   SELECT CASE (ierr)
   CASE (ierrno_fopen)
      WRITE (ioerr,'(A)') 'No file connected to unit '//TRIM(cunit)
   CASE (ierrno_read)
      WRITE (ioerr,'(A)') 'Read error occurred on file unit '//TRIM(cunit)
   CASE (ierrno_write)
      WRITE (ioerr,'(A)') 'Write error occurred on file unit '//TRIM(cunit)
   CASE (ierrno_fend)
      WRITE (ioerr,'(A)') 'End of file condition occurred on unit '//&
           & TRIM(cunit)
   END SELECT
ENDIF
IF (.NOT.PRESENT(abort)) THEN
   CALL error_abort(procname(pglev),ierr)
ELSEIF (abort) THEN
   CALL error_abort(procname(pglev),ierr)
ENDIF


RETURN

END SUBROUTINE error_file

!========================================================================

SUBROUTINE error_lbound_arr_int(ival,arrname,minval,matchmin,ndims,indx)
!************************************************************************
!
! *error_lbound_arr_int* Checks whether element 'ival' of an integer array
!                        is greater than or equal to 'minval' if 'matchmin' is
!                        .TRUE. or greater than 'minval' if 'matchmin' is
!                        .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmin
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: ival, minval, ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Input array element
!*arrname*   CHAR    Array name
!*minval*    INTEGER Imposed lower bound
!*matchmin*  LOGICAL Minimum is allowed if .TRUE.
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array index vector
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cminval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF ((matchmin.AND.ival.LT.minval).OR.&
  & (.NOT.matchmin.AND.ival.LE.minval)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; WRITE (cminval,'(I12)') minval
      cval = ADJUSTL(cval); cminval = ADJUSTL(cminval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for element '//TRIM(cindx(1))//&
                 & ' of integer array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of integer array '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of integer array '//TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//','//&
                 & TRIM(cindx(4))//') of integer array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      END SELECT
      IF (matchmin) THEN
         WRITE (ioerr,'(A)') 'Must be greater than or equal to: '&
                           & //TRIM(cminval)
      ELSE
         WRITE (ioerr,'(A)') 'Must be strictly greater than: '//TRIM(cminval)
      ENDIF
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_lbound_arr_int

!========================================================================

SUBROUTINE error_lbound_arr_real(rval,arrname,minval,matchmin,ndims,indx)
!************************************************************************
!
! *error_lbound_arr_real* Checks whether element 'rval' of a real array is
!                         greater than or equal to 'minval' if 'matchmin' is
!                         .TRUE. or greater than 'minval' if 'matchmin' is
!                         .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmin
CHARACTER (LEN=*), INTENT(IN) :: arrname
REAL (KIND=kndreal), INTENT(IN) :: rval, minval
INTEGER, INTENT(IN) :: ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      REAL    Input array element
!*arrname*   CHAR    Array name
!*minval*    REAL    Imposed lower bound
!*matchmin*  LOGICAL Minimum is allowed if .TRUE.
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array index vector
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cminval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF ((matchmin.AND.rval.LT.minval).OR.&
  & (.NOT.matchmin.AND.rval.LE.minval)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; WRITE (cminval,'(G15.7)') minval
      cval = ADJUSTL(cval); cminval = ADJUSTL(cminval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for element '&
                 & //TRIM(cindx(1))//' of real array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of real array '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of real array '//TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//','&
                 & //TRIM(cindx(4))//') of real array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      END SELECT
      IF (matchmin) THEN
         WRITE (ioerr,'(A)') 'Must be greater than or equal to: '&
                           & //TRIM(cminval)
      ELSE
         WRITE (ioerr,'(A)') 'Must be strictly greater than: '//TRIM(cminval)
      ENDIF
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_lbound_arr_real

!========================================================================

SUBROUTINE error_lbound_arr_double(rval,arrname,minval,matchmin,ndims,indx)
!************************************************************************
!
! *error_lbound_arr_double* Checks whether element 'rval' of a double
!                           precision array is greater than or equal to
!                           'minval' if 'matchmin' is .TRUE. or greater than
!                           'minval' if 'matchmin' is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.11
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmin
CHARACTER (LEN=*), INTENT(IN) :: arrname
REAL (KIND=kndrlong), INTENT(IN) :: rval, minval
INTEGER, INTENT(IN) :: ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      DOUBLE  Input array element
!*arrname*   CHAR    Array name
!*minval*    DOUBLE  Imposed lower bound
!*matchmin*  LOGICAL Minimum is allowed if .TRUE.
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array index vector
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cminval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF ((matchmin.AND.rval.LT.minval).OR.&
  & (.NOT.matchmin.AND.rval.LE.minval)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; WRITE (cminval,'(G15.7)') minval
      cval = ADJUSTL(cval); cminval = ADJUSTL(cminval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for element '&
                 & //TRIM(cindx(1))//' of real array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of real array '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of real array '//TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//','&
                 & //TRIM(cindx(4))//') of real array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      END SELECT
      IF (matchmin) THEN
         WRITE (ioerr,'(A)') 'Must be greater than or equal to: '&
                           & //TRIM(cminval)
      ELSE
         WRITE (ioerr,'(A)') 'Must be strictly greater than: '//TRIM(cminval)
      ENDIF
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_lbound_arr_double

!========================================================================

SUBROUTINE error_lbound_arr_struc_int(ival,arrname,compname,minval,matchmin,&
                                     & ndims,indx)
!************************************************************************
!
! *error_lbound_arr_struc_int* Checks whether integer component 'ival' of a
!        structure array is greater than or equal to 'minval' if 'matchmin' is
!        .TRUE. or greater than 'minval' if 'matchmin' is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmin
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
INTEGER, INTENT(IN) :: ival, minval, ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Component of derived type array element
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*minval*    INTEGER Imposed lower bound
!*matchmin*  LOGICAL Minimum is allowed if .TRUE.
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array index vector
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cminval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF ((matchmin.AND.ival.LT.minval).OR.&
  & (.NOT.matchmin.AND.ival.LE.minval)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; WRITE (cminval,'(I12)') minval
      cval = ADJUSTL(cval); cminval = ADJUSTL(cminval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '//&
                 & TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//&
                 & TRIM(cindx(4))//') of array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      END SELECT
      IF (matchmin) THEN
         WRITE (ioerr,'(A)') 'Must be greater than or equal to: '&
                           & //TRIM(cminval)
      ELSE
         WRITE (ioerr,'(A)') 'Must be strictly greater than: '//TRIM(cminval)
      ENDIF
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_lbound_arr_struc_int

!========================================================================

SUBROUTINE error_lbound_arr_struc_real(rval,arrname,compname,minval,matchmin,&
                                     & ndims,indx)
!************************************************************************
!
! *error_lbound_arr_struc_real* Checks whether real component 'rval' of a
!        structure array is greater than or equal to 'minval' if 'matchmin' is
!        .TRUE. or greater than 'minval' if 'matchmin' is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmin
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
INTEGER, INTENT(IN) :: ndims
REAL (KIND=kndreal), INTENT(IN) :: rval, minval
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      REAL    Component of derived type array element
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*minval*    REAL    Imposed lower bound
!*matchmin*  LOGICAL Minimum is allowed if .TRUE.
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array index vector
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cminval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF ((matchmin.AND.rval.LT.minval).OR.&
  & (.NOT.matchmin.AND.rval.LE.minval)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; WRITE (cminval,'(G15.7)') minval
      cval = ADJUSTL(cval); cminval = ADJUSTL(cminval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '//&
                 & TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//TRIM(cindx(4))//&
                 & ') of array '//TRIM(arrname)//': '//TRIM(cval)
      END SELECT
      IF (matchmin) THEN
         WRITE (ioerr,'(A)') 'Must be greater than or equal to: '&
                           & //TRIM(cminval)
      ELSE
         WRITE (ioerr,'(A)') 'Must be strictly greater than: '//TRIM(cminval)
      ENDIF
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_lbound_arr_struc_real

!========================================================================

SUBROUTINE error_lbound_arr_struc_double(rval,arrname,compname,minval,matchmin,&
                                       & ndims,indx)
!************************************************************************
!
! *error_lbound_arr_struc_double* Checks whether double precision component
!        'rval' of a structure array is greater than or equal to 'minval' if
!        'matchmin' is .TRUE. or greater than 'minval' if 'matchmin' is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.11
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmin
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
INTEGER, INTENT(IN) :: ndims
REAL (KIND=kndrlong), INTENT(IN) :: rval, minval
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      DOUBLE  Component of derived type array element
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*minval*    REAL    Imposed lower bound
!*matchmin*  LOGICAL Minimum is allowed if .TRUE.
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array index vector
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cminval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF ((matchmin.AND.rval.LT.minval).OR.&
  & (.NOT.matchmin.AND.rval.LE.minval)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; WRITE (cminval,'(G15.7)') minval
      cval = ADJUSTL(cval); cminval = ADJUSTL(cminval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '//&
                 & TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//TRIM(cindx(4))//&
                 & ') of array '//TRIM(arrname)//': '//TRIM(cval)
      END SELECT
      IF (matchmin) THEN
         WRITE (ioerr,'(A)') 'Must be greater than or equal to: '&
                           & //TRIM(cminval)
      ELSE
         WRITE (ioerr,'(A)') 'Must be strictly greater than: '//TRIM(cminval)
      ENDIF
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_lbound_arr_struc_double

!========================================================================

SUBROUTINE error_lbound_var_int(ival,varname,minval,matchmin)
!************************************************************************
!
! *error_lbound_var_int* Checks whether integer parameter 'ival' is greater
!                        than or equal to 'minval' if 'matchmin' is .TRUE. or
!                        greater than 'minval' if 'matchmin' is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmin
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: ival, minval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Integer variable
!*varname*   CHAR    Variable name
!*minval*    INTEGER Imposed lower bound
!*matchmin*  LOGICAL Minimum is allowed if .TRUE.
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cval, cminval


IF (matchmin) THEN
   IF (ival.LT.minval) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cval,'(I12)') ival; WRITE (cminval,'(I12)') minval
         cval = ADJUSTL(cval); cminval = ADJUSTL(cminval)
         WRITE (ioerr,'(A)') 'Invalid value for integer parameter '&
                           & //TRIM(varname)//': '//TRIM(cval)
         WRITE (ioerr,'(A)') 'Must be greater than or equal to: '&
                           & //TRIM(cminval)
      ENDIF
   ENDIF
ELSEIF (ival.LE.minval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; WRITE (cminval,'(I12)') minval
      WRITE (ioerr,'(A)') 'Invalid value for integer parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be strictly greater than: '//TRIM(cminval)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_lbound_var_int

!========================================================================

SUBROUTINE error_lbound_var_real(rval,varname,minval,matchmin)
!************************************************************************
!
! *error_lbound_var_real* Checks whether real parameter 'rval' is greater than
!                         or equal to 'minval' if 'matchmin' is .TRUE. or
!                         greater than 'minval' if 'matchmin' is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmin
CHARACTER (LEN=*), INTENT(IN) :: varname
REAL (KIND=kndreal), INTENT(IN) :: rval, minval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      REAL    Real variable
!*varname*   CHAR    Variable name
!*minval*    REAL    Imposed lower bound
!*matchmin*  LOGICAL Minimum is allowed if .TRUE.
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cval, cminval
    

IF (matchmin) THEN
   IF (rval.LT.minval) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cval,'(G15.7)') rval; WRITE (cminval,'(G15.7)') minval
         cval = ADJUSTL(cval); cminval = ADJUSTL(cminval)
         WRITE (ioerr,'(A)') 'Invalid value for real parameter '&
                           & //TRIM(varname)//': '//TRIM(cval)
         WRITE (ioerr,'(A)') 'Must be greater than or equal to: '&
                           & //TRIM(cminval)
      ENDIF
   ENDIF
ELSEIF (rval.LE.minval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; WRITE (cminval,'(G15.7)') minval
      WRITE (ioerr,'(A)') 'Invalid value for real parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be strictly greater than: '//TRIM(cminval)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_lbound_var_real

!========================================================================

SUBROUTINE error_lbound_var_double(rval,varname,minval,matchmin)
!************************************************************************
!
! *error_lbound_var_double* Checks whether double precision parameter 'rval'
!                           is greater than or equal to 'minval' if 'matchmin'
!                           is .TRUE. or greater than 'minval' if 'matchmin'
!                           is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.11
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmin
CHARACTER (LEN=*), INTENT(IN) :: varname
REAL (KIND=kndrlong), INTENT(IN) :: rval, minval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      DOUBLE  Real variable
!*varname*   CHAR    Variable name
!*minval*    DOUBLE  Imposed lower bound
!*matchmin*  LOGICAL Minimum is allowed if .TRUE.
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cval, cminval
    

IF (matchmin) THEN
   IF (rval.LT.minval) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cval,'(G15.7)') rval; WRITE (cminval,'(G15.7)') minval
         cval = ADJUSTL(cval); cminval = ADJUSTL(cminval)
         WRITE (ioerr,'(A)') 'Invalid value for real parameter '&
                           & //TRIM(varname)//': '//TRIM(cval)
         WRITE (ioerr,'(A)') 'Must be greater than or equal to: '&
                           & //TRIM(cminval)
      ENDIF
   ENDIF
ELSEIF (rval.LE.minval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; WRITE (cminval,'(G15.7)') minval
      WRITE (ioerr,'(A)') 'Invalid value for real parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be strictly greater than: '//TRIM(cminval)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_lbound_var_double

!========================================================================

SUBROUTINE error_limits_arr_int(ival,arrname,minval,maxval,ndims,indx)
!************************************************************************
!
! *error_limits_arr_int* Checks whether element 'ival' of an integer array
!                        is between bounds 'minval' and 'maxval'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: ival, minval, maxval, ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Integer array element
!*arrname*   CHAR    Array name
!*minval*    INTEGER Imposed lower bound
!*maxval*    INTEGER Imposed upper bound
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cminval, cmaxval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (ival.LT.minval.OR.ival.GT.maxval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; cval = ADJUSTL(cval)
      WRITE (cminval,'(I12)') minval; cminval = ADJUSTL(cminval)
      WRITE (cmaxval,'(I12)') maxval; cmaxval = ADJUSTL(cmaxval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for element '//TRIM(cindx(1))//&
                 & ' of integer array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of integer array '//arrname//': '//&
                 & TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of integer array '//TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//','&
                 & //TRIM(cindx(4))//') of integer array '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      END SELECT
      WRITE (ioerr,8001) 'Must be between: ', TRIM(cminval), TRIM(cmaxval)
   ENDIF
ENDIF


RETURN

8001 FORMAT (A,A,' and ',A)

END SUBROUTINE error_limits_arr_int

!========================================================================

SUBROUTINE error_limits_arr_real(rval,arrname,minval,maxval,ndims,indx)
!************************************************************************
!
! *error_limits_arr_real* Checks whether element 'rval' of a real array is
!                         between bounds 'minval' and 'maxval'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: ndims
REAL (KIND=kndreal), INTENT(IN) :: rval, minval, maxval
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      REAL    Real array element
!*arrname*   CHAR    Array name
!*minval*    REAL    Imposed lower bound
!*maxval*    REAL    Imposed upper bound
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cminval, cmaxval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (rval.LT.minval.OR.rval.GT.maxval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; cval = ADJUSTL(cval)
      WRITE (cminval,'(G15.7)') minval; cminval = ADJUSTL(cminval)
      WRITE (cmaxval,'(G15.7)') maxval; cmaxval = ADJUSTL(cmaxval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for element '//TRIM(cindx(1))//&
                 & ' of real array '//TRIM(arrname)//': '//cval
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of real array '//TRIM(arrname)//&
                 & ': '//cval
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of real array '//TRIM(arrname)//': '//cval
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//','&
                 & //TRIM(cindx(4))//') of real array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      END SELECT
      WRITE (ioerr,8001) 'Must be between: ', TRIM(cminval), TRIM(cmaxval)
   ENDIF
ENDIF


RETURN

8001 FORMAT (A,A,' and ',A)

END SUBROUTINE error_limits_arr_real

!========================================================================

SUBROUTINE error_limits_arr_double(rval,arrname,minval,maxval,ndims,indx)
!************************************************************************
!
! *error_limits_arr_double* Checks whether element 'rval' of a double precision
!                           array is between bounds 'minval' and 'maxval'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.11
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: ndims
REAL (KIND=kndrlong), INTENT(IN) :: rval, minval, maxval
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      DOUBLE  Real array element
!*arrname*   CHAR    Array name
!*minval*    DOUBLE  Imposed lower bound
!*maxval*    DOUBLE  Imposed upper bound
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cminval, cmaxval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (rval.LT.minval.OR.rval.GT.maxval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; cval = ADJUSTL(cval)
      WRITE (cminval,'(G15.7)') minval; cminval = ADJUSTL(cminval)
      WRITE (cmaxval,'(G15.7)') maxval; cmaxval = ADJUSTL(cmaxval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for element '//TRIM(cindx(1))//&
                 & ' of real array '//TRIM(arrname)//': '//cval
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of real array '//TRIM(arrname)//&
                 & ': '//cval
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of real array '//TRIM(arrname)//': '//cval
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//','&
                 & //TRIM(cindx(4))//') of real array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      END SELECT
      WRITE (ioerr,8001) 'Must be between: ', TRIM(cminval), TRIM(cmaxval)
   ENDIF
ENDIF


RETURN

8001 FORMAT (A,A,' and ',A)

END SUBROUTINE error_limits_arr_double

!========================================================================

SUBROUTINE error_limits_arr_struc_int(ival,arrname,compname,minval,maxval,&
                                    & ndims,indx)
!************************************************************************
!
! *error_limits_arr_struc_int* Checks whether integer component 'ival' of a
!                              structure array is between bounds 'minval' and
!                              'maxval'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
INTEGER, INTENT(IN) :: ival, minval, maxval, ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Value of derived type array component
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*minval*    INTEGER Imposed lower bound
!*maxval*    INTEGER Imposed upper bound
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cminval, cmaxval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (ival.LT.minval.OR.ival.GT.maxval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; cval = ADJUSTL(cval)
      WRITE (cminval,'(I12)') minval; cminval = ADJUSTL(cminval)
      WRITE (cmaxval,'(I12)') maxval; cmaxval = ADJUSTL(cmaxval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)&
                 & //' and element '//TRIM(cindx(1))//' of array '//&
                 & TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '//&
                 & TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//TRIM(cindx(4))//&
                 & ') of array '//TRIM(arrname)//': '//TRIM(cval)
      END SELECT
      WRITE (ioerr,8001) 'Must be between: ', TRIM(cminval), TRIM(cmaxval)
   ENDIF
ENDIF


RETURN

8001 FORMAT (A,A,' and ',A)

END SUBROUTINE error_limits_arr_struc_int

!========================================================================

SUBROUTINE error_limits_arr_struc_real(rval,arrname,compname,minval,maxval,&
                                     & ndims,indx)
!************************************************************************
!
! *error_limits_arr_struc_real* Checks whether real component 'rval' of a
!                               derived type array is between bounds 'minval'
!                               and 'maxval'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
INTEGER, INTENT(IN) :: ndims
REAL (KIND=kndreal), INTENT(IN) :: rval, minval, maxval
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      REAL    Value of derived type array component
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*minval*    REAL    Imposed lower bound
!*maxval*    REAL    Imposed upper bound
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cminval, cmaxval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (rval.LT.minval.OR.rval.GT.maxval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; cval = ADJUSTL(cval)
      WRITE (cminval,'(G15.7)') minval; cminval = ADJUSTL(cminval)
      WRITE (cmaxval,'(G15.7)') maxval; cmaxval = ADJUSTL(cmaxval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array '//TRIM(arrname)//': '//cval
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//' of array '//TRIM(arrname)//': '//cval
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '//&
                 & TRIM(arrname)//': '//cval
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//TRIM(cindx(4))//&
                 & ') of array '//TRIM(arrname)//': '//TRIM(cval)
      END SELECT
      WRITE (ioerr,8001) 'Must be between: ', TRIM(cminval), TRIM(cmaxval)
   ENDIF
ENDIF


RETURN

8001 FORMAT (A,A,' and ',A)

END SUBROUTINE error_limits_arr_struc_real

!========================================================================

SUBROUTINE error_limits_arr_struc_double(rval,arrname,compname,minval,maxval,&
                                       & ndims,indx)
!************************************************************************
!
! *error_limits_arr_struc_double* Checks whether double precision component
!                                 'rval' of a derived type array is between
!                                 bounds 'minval' and 'maxval'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.11
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
INTEGER, INTENT(IN) :: ndims
REAL (KIND=kndrlong), INTENT(IN) :: rval, minval, maxval
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      DOUBLE  Value of derived type array component
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*minval*    DOUBLE  Imposed lower bound
!*maxval*    DOUBLE  Imposed upper bound
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cminval, cmaxval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (rval.LT.minval.OR.rval.GT.maxval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; cval = ADJUSTL(cval)
      WRITE (cminval,'(G15.7)') minval; cminval = ADJUSTL(cminval)
      WRITE (cmaxval,'(G15.7)') maxval; cmaxval = ADJUSTL(cmaxval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array '//TRIM(arrname)//': '//cval
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//' of array '//TRIM(arrname)//': '//cval
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '//&
                 & TRIM(arrname)//': '//cval
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//TRIM(cindx(4))//&
                 & ') of array '//TRIM(arrname)//': '//TRIM(cval)
      END SELECT
      WRITE (ioerr,8001) 'Must be between: ', TRIM(cminval), TRIM(cmaxval)
   ENDIF
ENDIF


RETURN

8001 FORMAT (A,A,' and ',A)

END SUBROUTINE error_limits_arr_struc_double

!========================================================================

SUBROUTINE error_limits_var_int(ival,varname,minval,maxval)
!************************************************************************
!
! *error_limits_var_int* Checks whether integer parameter 'ival' is between
!                        bounds 'minval' and 'maxval'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: ival, minval, maxval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Integer variable
!*varname*   CHAR    Variable name
!*minval*    INTEGER Imposed lower bound
!*maxval*    INTEGER Imposed upper bound
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cminval, cmaxval, cval


IF (ival.LT.minval.OR.ival.GT.maxval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; cval = ADJUSTL(cval)
      WRITE (cminval,'(I12)') minval; cminval = ADJUSTL(cminval)
      WRITE (cmaxval,'(I12)') maxval; cmaxval = ADJUSTL(cmaxval)
      WRITE (ioerr,'(A)') 'Invalid value for integer parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,8001) 'Must be between: ', TRIM(cminval), TRIM(cmaxval)
   ENDIF
ENDIF


RETURN

8001 FORMAT (A,A,' and ',A)

END SUBROUTINE error_limits_var_int

!========================================================================

SUBROUTINE error_limits_var_real(rval,varname,minval,maxval)
!************************************************************************
!
! *error_limits_var_real* Checks whether real parameter 'rval' is between
!                         bounds 'minval' and 'maxval'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
REAL (KIND=kndreal), INTENT(IN) :: rval, minval, maxval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      REAL    Real variable
!*varname*   CHAR    Variable name
!*minval*    REAL    Imposed lower bound
!*maxval*    REAL    Imposed upper bound
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cminval, cmaxval, cval


IF (rval.LT.minval.OR.rval.GT.maxval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; cval = ADJUSTL(cval)
      WRITE (cminval,'(G15.7)') minval; cminval = ADJUSTL(cminval)
      WRITE (cmaxval,'(G15.7)') maxval; cmaxval = ADJUSTL(cmaxval)
      WRITE (ioerr,'(A)') 'Invalid value for real parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,8001) 'Must be between: ', TRIM(cminval), TRIM(cmaxval)
   ENDIF
ENDIF


RETURN

8001 FORMAT (A,A,' and ',A)

END SUBROUTINE error_limits_var_real

!========================================================================

SUBROUTINE error_limits_var_double(rval,varname,minval,maxval)
!************************************************************************
!
! *error_limits_var_double* Checks whether double precision parameter 'rval'
!                           is between bounds 'minval' and 'maxval'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.11
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
REAL (KIND=kndrlong), INTENT(IN) :: rval, minval, maxval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      DOUBLE  Real variable
!*varname*   CHAR    Variable name
!*minval*    DOUBLE  Imposed lower bound
!*maxval*    DOUBLE  Imposed upper bound
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cminval, cmaxval, cval


IF (rval.LT.minval.OR.rval.GT.maxval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; cval = ADJUSTL(cval)
      WRITE (cminval,'(G15.7)') minval; cminval = ADJUSTL(cminval)
      WRITE (cmaxval,'(G15.7)') maxval; cmaxval = ADJUSTL(cmaxval)
      WRITE (ioerr,'(A)') 'Invalid value for real parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,8001) 'Must be between: ', TRIM(cminval), TRIM(cmaxval)
   ENDIF
ENDIF


RETURN

8001 FORMAT (A,A,' and ',A)

END SUBROUTINE error_limits_var_double

!========================================================================

SUBROUTINE error_mult(ival,varname,multval)
!************************************************************************
!
! *error_mult* Checks whether integer parameter 'ival' is a divisor of
!              'multval' 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
! Module calls - mult_index
!
!************************************************************************
!
USE utility_routines, ONLY: mult_index

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: ival, multval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Integer variable
!*varname*   CHAR    Variable name
!*multval*   INTEGER Divisor value
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cval, cmultval


IF (.NOT.mult_index(multval,ival)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; WRITE (cmultval,'(I12)') multval
      cval = ADJUSTL(cval); cmultval = ADJUSTL(cmultval)
      WRITE (ioerr,'(A)') 'Invalid value for integer parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be a divisor of: '//TRIM(cmultval)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_mult

!========================================================================

SUBROUTINE error_proc(ierr,abort)
!************************************************************************
!
! *error_proc* Display error number which occurred in a routine and abort
!              program if necessary
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - abort program if abort is not present or present and .TRUE.
!
! Module calls - error_abort
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: abort
INTEGER, INTENT(IN) :: ierr

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ierr*      INTEGER Error number
!*abort*     LOGICAL Abort program if .TRUE.
! 
!------------------------------------------------------------------------------
!
!*Local variables
! 
CHARACTER (LEN=12) :: cerr


IF (ierr.NE.0) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cerr,'(I12)') ierr; cerr = ADJUSTL(cerr)
      WRITE (ioerr,'(A)') 'Error '//TRIM(cerr)//' occurred in '//&
                        & TRIM(procname(pglev))
   ENDIF
   IF (PRESENT(abort)) THEN
      IF (abort) CALL error_abort(procname(pglev),ierrno_runval)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_proc

!========================================================================

SUBROUTINE error_shape(arrshape,arrname,fixshape,ndim)
!************************************************************************
!
! *error_shape* Checks whether 'arrshape' (shape of an array with name
!               'arrname') equals the shape given by 'fixshape'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: ndim
INTEGER, INTENT(IN) , DIMENSION(ndim) :: arrshape, fixshape

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*arrshape*  INTEGER Shape of input array
!*arrname*   CHAR    Array name
!*fixshape*  INTEGER Imposed shape
!*ndim*      INTEGER Array rank
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndim) :: cdim 
CHARACTER (LEN=200) :: cline


IF (ANY(arrshape-fixshape.NE.0)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      n_110: DO n=1,ndim
         WRITE (cdim(n),'(I12)') arrshape(n)
         IF (n.EQ.1) THEN
            cline = TRIM(ADJUSTL(cdim(n)))
         ELSE
            cline = TRIM(cline)//' '//TRIM(ADJUSTL(cdim(n)))
         ENDIF
      ENDDO n_110
      WRITE (ioerr,'(A)') 'Wrong shape for array '//TRIM(arrname)//': '//&
                        & TRIM(cline)
      n_120: DO n=1,ndim
         WRITE (cdim(n),'(I12)') fixshape(n)
         cdim(n) = ADJUSTL(cdim(n))
         IF (n.EQ.1) THEN
            cline = TRIM(ADJUSTL(cdim(n)))
         ELSE
            cline = TRIM(cline)//' '//TRIM(ADJUSTL(cdim(n)))
         ENDIF
      ENDDO n_120
      WRITE (ioerr,'(9A)') 'Should be: '//TRIM(cline)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_shape

!========================================================================

SUBROUTINE error_ubound_arr_int(ival,arrname,maxval,matchmax,ndims,indx)
!************************************************************************
!
! *error_ubound_arr_int* Checks whether element 'ival' of an integer array
!                        is lower than or equal to 'maxval' if 'matchmax' is
!                        .TRUE. or lower than 'maxval if 'matchmax' is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmax
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: ival, maxval, ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Integer array element
!*arrname*   CHAR    Array name
!*maxval*    INTEGER Imposed upper bound
!*matchmax*  LOGICAL Maximum is allowed if .TRUE.
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cmaxval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF ((matchmax.AND.ival.GT.maxval).OR.&
  & (.NOT.matchmax.AND.ival.GE.maxval)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; WRITE (cmaxval,'(I12)') maxval
      cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for element '//TRIM(cindx(1))//&
                 & ' of integer array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of integer array '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//')'//&
                 & ' of integer array '//TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 &','//TRIM(cindx(2))//','//TRIM(cindx(3))//','&
                 & //TRIM(cindx(4))//') of integer array '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      END SELECT
      IF (matchmax) THEN
         WRITE (ioerr,'(A)') 'Must be lower than or equal to: '&
                           & //TRIM(cmaxval)
      ELSE
         WRITE (ioerr,'(A)') 'Must be strictly lower than: '//TRIM(cmaxval)
      ENDIF
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_ubound_arr_int

!========================================================================

SUBROUTINE error_ubound_arr_real(rval,arrname,maxval,matchmax,ndims,indx)
!************************************************************************
!
! *error_ubound_arr_real* Checks whether element 'rval' of a real array is
!                         lower than or equal to 'maxval' if 'matchmax' is
!                         .TRUE. or lower than 'maxval' if 'matchmax' is
!                         .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmax
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: ndims
REAL (KIND=kndreal), INTENT(IN) :: rval, maxval
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      REAL    Real array element
!*arrname*   CHAR    Array name
!*maxval*    REAL    Imposed upper bound
!*matchmax*  LOGICAL Maximum is allowed if .TRUE.
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cmaxval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF ((matchmax.AND.rval.GT.maxval).OR.&
  & (.NOT.matchmax.AND.rval.GE.maxval)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; WRITE (cmaxval,'(G15.7)') maxval
      cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for element '//TRIM(cindx(1))//&
                 & ' of real array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of real array '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of real array '//TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for element ('&
                 & //TRIM(cindx(1))//','//TRIM(cindx(2))//','//&
                 & TRIM(cindx(3))//','//TRIM(cindx(4))//') of real array '//&
                 & TRIM(arrname)//': '//TRIM(cval)
      END SELECT
      IF (matchmax) THEN
         WRITE (ioerr,'(A)') 'Must be lower than or equal to: '&
                           & //TRIM(cmaxval)
      ELSE
         WRITE (ioerr,'(A)') 'Must be strictly lower than: '//TRIM(cmaxval)
      ENDIF
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_ubound_arr_real

!========================================================================

SUBROUTINE error_ubound_arr_double(rval,arrname,maxval,matchmax,ndims,indx)
!************************************************************************
!
! *error_ubound_arr_double* Checks whether element 'rval' of a double precision
!                           array is lower than or equal to 'maxval' if
!                           'matchmax' is .TRUE. or lower than 'maxval'
!                           if 'matchmax' is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.11
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmax
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: ndims
REAL (KIND=kndrlong), INTENT(IN) :: rval, maxval
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      DOUBLE  Real array element
!*arrname*   CHAR    Array name
!*maxval*    DOUBLE  Imposed upper bound
!*matchmax*  LOGICAL Maximum is allowed if .TRUE.
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cmaxval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF ((matchmax.AND.rval.GT.maxval).OR.&
  & (.NOT.matchmax.AND.rval.GE.maxval)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; WRITE (cmaxval,'(G15.7)') maxval
      cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for element '//TRIM(cindx(1))//&
                 & ' of real array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of real array '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of real array '//TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for element ('&
                 & //TRIM(cindx(1))//','//TRIM(cindx(2))//','//&
                 & TRIM(cindx(3))//','//TRIM(cindx(4))//') of real array '//&
                 & TRIM(arrname)//': '//TRIM(cval)
      END SELECT
      IF (matchmax) THEN
         WRITE (ioerr,'(A)') 'Must be lower than or equal to: '&
                           & //TRIM(cmaxval)
      ELSE
         WRITE (ioerr,'(A)') 'Must be strictly lower than: '//TRIM(cmaxval)
      ENDIF
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_ubound_arr_double

!========================================================================

SUBROUTINE error_ubound_arr_struc_int(ival,arrname,compname,maxval,matchmax,&
                                    & ndims,indx)
!************************************************************************
!
! *error_ubound_arr_struc_int* Checks whether integer
!      component 'ival' of a structure array is lower than or equal to 'maxval'
!      if 'matchmax' is .TRUE. or lower than 'maxval' if 'matchmax' is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmax
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
INTEGER, INTENT(IN) :: ival, maxval, ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Component of derived type array
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*maxval*    INTEGER Imposed upper bound
!*matchmax*  LOGICAL Maximum is allowed if .TRUE.
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cmaxval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF ((matchmax.AND.ival.GT.maxval).OR.&
  & (.NOT.matchmax.AND.ival.GE.maxval)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; WRITE (cmaxval,'(I12)') maxval
      cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '//&
                 & TRIM(arrname)//': '//&
                 & TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//&
                 & TRIM(cindx(4))//') of array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      END SELECT
      IF (matchmax) THEN
         WRITE (ioerr,'(A)') 'Must be lower than or equal to: '&
                           & //TRIM(cmaxval)
      ELSE
         WRITE (ioerr,'(A)') 'Must be strictly lower than: '//TRIM(cmaxval)
      ENDIF
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_ubound_arr_struc_int

!========================================================================

SUBROUTINE error_ubound_arr_struc_real(rval,arrname,compname,maxval,matchmax,&
                                     & ndims,indx)
!************************************************************************
!
! *error_ubound_arr_struc_real* Checks whether component 'rval'
!       of a structure array is lower than or equal to 'maxval' if 'matchmax'
!       is .TRUE. or lower than 'maxval' if 'matchmax' is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmax
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
REAL (KIND=kndreal), INTENT(IN) :: rval, maxval
INTEGER, INTENT(IN) :: ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      REAL    Component of derived type array
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*maxval*    REAL    Imposed upper bound
!*matchmax*  LOGICAL Maximum is allowed if .TRUE.
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cmaxval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF ((matchmax.AND.rval.GT.maxval).OR.&
  & (.NOT.matchmax.AND.rval.GE.maxval)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; WRITE (cmaxval,'(G15.7)') maxval
      cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '//&
                 & TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//TRIM(cindx(4))//&
                 & ') of array '//TRIM(arrname)//': '//TRIM(cval)
      END SELECT
      IF (matchmax) THEN
         WRITE (ioerr,'(A)') 'Must be lower than or equal to: '&
                           & //TRIM(cmaxval)
      ELSE
         WRITE (ioerr,'(A)') 'Must be strictly lower than: '//TRIM(cmaxval)
      ENDIF
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_ubound_arr_struc_real

!========================================================================

SUBROUTINE error_ubound_arr_struc_double(rval,arrname,compname,maxval,matchmax,&
                                       & ndims,indx)
!************************************************************************
!
! *error_ubound_arr_struc_double* Checks whether double precision component
!       'rval' of a structure array is lower than or equal to 'maxval' if
!       'matchmax' is .TRUE. or lower than 'maxval' if 'matchmax' is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.11
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmax
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
REAL (KIND=kndrlong), INTENT(IN) :: rval, maxval
INTEGER, INTENT(IN) :: ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      DOUBLE  Component of derived type array
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*maxval*    DOUBLE  Imposed upper bound
!*matchmax*  LOGICAL Maximum is allowed if .TRUE.
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cmaxval, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF ((matchmax.AND.rval.GT.maxval).OR.&
  & (.NOT.matchmax.AND.rval.GE.maxval)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; WRITE (cmaxval,'(G15.7)') maxval
      cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '//&
                 & TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for real component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//TRIM(cindx(4))//&
                 & ') of array '//TRIM(arrname)//': '//TRIM(cval)
      END SELECT
      IF (matchmax) THEN
         WRITE (ioerr,'(A)') 'Must be lower than or equal to: '&
                           & //TRIM(cmaxval)
      ELSE
         WRITE (ioerr,'(A)') 'Must be strictly lower than: '//TRIM(cmaxval)
      ENDIF
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_ubound_arr_struc_double

!========================================================================

SUBROUTINE error_ubound_var_int(ival,varname,maxval,matchmax)
!************************************************************************
!
! *error_ubound_var_int* Checks whether integer parameter ival is lower than
!                        or equal to 'maxval' if 'matchmax' is .TRUE. or lower
!                        than 'maxval' if 'matchmax' is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmax
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: ival, maxval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Integer variable
!*varname*   CHAR    Variable name
!*maxval*    INTEGER Imposed upper bound
!*matchmax*  LOGICAL Maximum is allowed if .TRUE.
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cval, cmaxval


IF (matchmax) THEN
   IF (ival.GT.maxval) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cval,'(I12)') ival; WRITE (cmaxval,'(I12)') maxval
         cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
         WRITE (ioerr,'(A)') 'Invalid value for integer parameter '&
                           & //TRIM(varname)//': '//TRIM(cval)
         WRITE (ioerr,'(A)') 'Must be lower than or equal to: '&
                           & //TRIM(cmaxval)
      ENDIF
   ENDIF
ELSEIF (ival.GE.maxval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; WRITE (cmaxval,'(I12)') maxval
      cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
      WRITE (ioerr,'(A)') 'Invalid value for integer parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be strictly lower than: '//TRIM(cmaxval)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_ubound_var_int

!========================================================================

SUBROUTINE error_ubound_var_real(rval,varname,maxval,matchmax)
!************************************************************************
!
! *error_ubound_var_real* Checks whether real parameter 'rval' is lower than
!                         or equal to 'maxval' if 'matchmax' is .TRUE. or lower
!                         than 'maxval' if 'matchmax' is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmax
CHARACTER (LEN=*), INTENT(IN) :: varname
REAL (KIND=kndreal), INTENT(IN) :: rval, maxval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      REAL    Real variable
!*varname*   CHAR    Variable name
!*maxval*    REAL    Imposed upper bound
!*matchmax*  LOGICAL Maximum is allowed if .TRUE.
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cval, cmaxval


IF (matchmax) THEN
   IF (rval.GT.maxval) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cval,'(G15.7)') rval; WRITE (cmaxval,'(G15.7)') maxval
         cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
         WRITE (ioerr,'(A)') 'Invalid value for real parameter '&
                           & //TRIM(varname)//': '//TRIM(cval)
         WRITE (ioerr,'(A)') 'Must be lower than or equal to: '&
                           & //TRIM(cmaxval)
      ENDIF
   ENDIF
ELSEIF (rval.GE.maxval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; WRITE (cmaxval,'(G15.7)') maxval
      cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
      WRITE (ioerr,'(A)') 'Invalid value for real parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be strictly greater than: '//TRIM(cmaxval)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_ubound_var_real

!========================================================================

SUBROUTINE error_ubound_var_double(rval,varname,maxval,matchmax)
!************************************************************************
!
! *error_ubound_var_double* Checks whether double preicsion parameter 'rval' is
!                           lower than or equal to 'maxval' if 'matchmax' is
!                           .TRUE. or lower than 'maxval' if 'matchmax' is
!                           .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.11
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmax
CHARACTER (LEN=*), INTENT(IN) :: varname
REAL (KIND=kndrlong), INTENT(IN) :: rval, maxval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      DOUBLE  Real variable
!*varname*   CHAR    Variable name
!*maxval*    DOUBLE  Imposed upper bound
!*matchmax*  LOGICAL Maximum is allowed if .TRUE.
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cval, cmaxval


IF (matchmax) THEN
   IF (rval.GT.maxval) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cval,'(G15.7)') rval; WRITE (cmaxval,'(G15.7)') maxval
         cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
         WRITE (ioerr,'(A)') 'Invalid value for real parameter '&
                           & //TRIM(varname)//': '//TRIM(cval)
         WRITE (ioerr,'(A)') 'Must be lower than or equal to: '&
                           & //TRIM(cmaxval)
      ENDIF
   ENDIF
ELSEIF (rval.GE.maxval) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') rval; WRITE (cmaxval,'(G15.7)') maxval
      cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
      WRITE (ioerr,'(A)') 'Invalid value for real parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be strictly greater than: '//TRIM(cmaxval)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_ubound_var_double

!========================================================================

SUBROUTINE error_vals_arr_char(cval,arrname,string,ndims,indx)
!************************************************************************
!
! *error_vals_arr_char* Checks whether element 'cval' of a character array
!                       has a allowed value given by 'string'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.4
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname, cval, string
INTEGER, INTENT(IN) :: ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cval*      CHAR    Element of character array
!*arrname*   CHAR    Array name
!*string*    CHAR    String to which the array element must be equal
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (SCAN(string,TRIM(cval)).EQ.0) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for element '//TRIM(cindx(1))//&
                 & ' of character array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of character array '//&
                 & TRIM(arrname)//': '//TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of character array '//TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//','&
                 & //TRIM(cindx(4))//') of character array '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      END SELECT
      WRITE (ioerr,'(A)') 'Allowed values are: '//TRIM(string)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_vals_arr_char

!========================================================================

SUBROUTINE error_vals_arr_int(ival,arrname,nlist,ndims,ifix,indx)
!************************************************************************
!
! *error_vals_arr_int* Checks whether element 'ival' of an integer array
!                      equals one of the values in the list vector 'ifix'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
! Module calls - digit_number_int, index_position
!
!************************************************************************
!
USE utility_routines, ONLY: digit_number_int, index_position

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: ival, ndims, nlist
INTEGER, INTENT(IN), DIMENSION(nlist) :: ifix
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Element of integer array
!*arrname*   CHAR    Array name
!*nlist*     INTEGER Size of ifix
!*ndims*     INTEGER Array rank
!*ifix*      INTEGER List of integers which should include ival
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cdec, clist, cval
CHARACTER (LEN=29) :: cfmt
INTEGER :: ipos, n, ndec
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (nlist.EQ.0) RETURN

ipos = index_position(ival,ifix)

IF (ipos.EQ.0) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; cval = ADJUSTL(cval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for element '//TRIM(cindx(1))//&
                           & ' of integer array '//TRIM(arrname)//': '//&
                           & TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))&
                 & //','//TRIM(cindx(2))//') of integer array '//&
                 & TRIM(arrname)//': '//TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of integer array '//TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//','&
                 & //TRIM(cindx(4))//') of integer array '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      END SELECT
      ndec = 1
      n_120: DO n=1,nlist
         ndec = MAX(ndec,digit_number_int(ifix(n)))
      ENDDO n_120
      WRITE (cdec,'(I12)') ndec+1; cdec = ADJUSTL(cdec)
      WRITE (clist,'(I12)') nlist; clist = ADJUSTL(clist)
      cfmt = '(A,'//TRIM(clist)//'I'//TRIM(cdec)//')'
      WRITE (ioerr,cfmt) 'Allowed values are: ', ifix  
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_vals_arr_int

!========================================================================

SUBROUTINE error_vals_arr_struc_char(cval,arrname,compname,string,ndims,indx)
!************************************************************************
!
! *error_vals_arr_struc_char* Checks whether character component 'cval' of a
!                             structure array has a allowed value given by
!                             'string'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname, cval, string
INTEGER, INTENT(IN) :: ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cval*      CHAR    Component of derived type array
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*string*    CHAR    String to which the character component must be equal
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx

IF (SCAN(string,'"'//TRIM(cval)//'"').EQ.0) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for character component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for character component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for character component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '//&
                 & TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for character component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//&
                 & TRIM(cindx(4))//') of array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      END SELECT
      WRITE (ioerr,'(A)') 'Allowed values are: '//TRIM(string)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_vals_arr_struc_char

!========================================================================

SUBROUTINE error_vals_arr_struc_int(ival,arrname,compname,nlist,ndims,ifix,indx)
!************************************************************************
!
! *error_vals_arr_struc_int* Checks whether integer component 'ival' of a
!                            structure array equals one of the values in the
!                            list 'ifix'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
! Module calls - digit_number_int, index_position
!
!************************************************************************
!
USE utility_routines, ONLY: digit_number_int, index_position

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
INTEGER, INTENT(IN) :: ival, ndims, nlist
INTEGER, INTENT(IN), DIMENSION(nlist) :: ifix
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Component of derived type array
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*nlist*     INTEGER Size of ifix
!*ndims*     INTEGER Size of indx
!*ifix*      INTEGER Integer value to which the array component must be equal
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cdec, clist, cval
CHARACTER (LEN=29) :: cfmt
INTEGER :: ipos, n, ndec
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (nlist.EQ.0) RETURN

ipos = index_position(ival,ifix)

IF (ipos.EQ.0) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; cval = ADJUSTL(cval)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '//&
                 & TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//TRIM(cindx(4))//&
                 & ') of array '//TRIM(arrname)//': '//TRIM(cval)
      END SELECT
      ndec = 1
      n_120: DO n=1,nlist
         ndec = MAX(ndec,digit_number_int(ifix(n)))
      ENDDO n_120
      WRITE (cdec,'(I12)') ndec+1; cdec = ADJUSTL(cdec)
      WRITE (clist,'(I12)') nlist; clist = ADJUSTL(clist)
      cfmt = '(A,'//TRIM(clist)//'I'//TRIM(cdec)//')'
      WRITE (ioerr,cfmt) 'Allowed values are: ', ifix  
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_vals_arr_struc_int

!========================================================================

SUBROUTINE error_vals_var_char(cval,varname,string)
!************************************************************************
!
! *error_vals_var_char* Checks whether character parameter 'cval' has an
!                       allowed value as given by 'string'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - allowed values in string must be between ""
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: cval, string, varname

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cval*      CHAR    Character variable
!*varname*   CHAR    Variable name
!*string*    CHAR    Character string to which the variable must be equal
! 
!------------------------------------------------------------------------------
!

IF (SCAN(string,'"'//TRIM(cval)//'"').EQ.0) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (ioerr,'(A)') 'Invalid value for character parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Allowed values are: '//TRIM(string)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_vals_var_char

!========================================================================

SUBROUTINE error_vals_var_int(ival,varname,nlist,ifix)
!************************************************************************
!
! *error_vals_var_int* Checks whether integer parameter 'ival' equals one
!                      of the values in the vector list 'ifix'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
! Module calls - digit_number_int, index_position
!
!************************************************************************
!
USE utility_routines, ONLY: digit_number_int, index_position

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: ival, nlist
INTEGER, INTENT(IN), DIMENSION(nlist) :: ifix

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Integer variable
!*varname*   CHAR    Variable name
!*nlist*     INTEGER Size of the vector list ifix
!*ifix*      INTEGER Integer value to which the variable must be equal
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cdec, clist, cval
CHARACTER (LEN=29) :: cfmt
INTEGER :: ipos, n, ndec


IF (nlist.EQ.0) RETURN

ipos = index_position(ival,ifix)

IF (ipos.EQ.0) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; cval = ADJUSTL(cval)
      WRITE (ioerr,'(A)') 'Invalid value for integer parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      ndec = 1
      n_110: DO n=1,nlist
         ndec = MAX(ndec,digit_number_int(ifix(n)))
      ENDDO n_110
      WRITE (cdec,'(I12)') ndec+1; cdec = ADJUSTL(cdec)
      WRITE (clist,'(I12)') nlist; clist = ADJUSTL(clist)
      cfmt = '(A,'//TRIM(clist)//'I'//TRIM(cdec)//')'
      WRITE (ioerr,cfmt) 'Allowed values are: ', ifix  
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_vals_var_int

!========================================================================

SUBROUTINE error_value_arr_int(ival,arrname,ifix,ndims,indx)
!************************************************************************
!
! *error_value_arr_int* Checks whether element 'ival' of an integer array
!                       equals 'ifix'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: ival, ifix, ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx(ndims)

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Integer array element
!*arrname*   CHAR    Array name
!*ifix*      INTEGER Integer value to which the array element must be equal
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Arry vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cfix, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (ival.NE.ifix) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; WRITE(cfix,'(I12)') ifix
      cval = ADJUSTL(cval); cfix = ADJUSTL(cfix)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for element '//TRIM(cindx(1))//&
                 & ' of integer array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of integer array '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of integer array '//TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//','&
                 & //TRIM(cindx(4))//') of integer array '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      END SELECT
      WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cfix)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_value_arr_int

!========================================================================

SUBROUTINE error_value_arr_log(lval,arrname,lfix,ndims,indx)
!************************************************************************
!
! *error_value_arr_log* Checks whether element 'lval' of a logical array
!                       equals 'lfix'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
LOGICAL, INTENT(IN) :: lval, lfix
INTEGER, INTENT(IN) :: ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*lval*      LOGICAL Integer array element
!*arrname*   CHAR    Array name
!*lfix*      LOGICAL Logical value to which the array element must be equal
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Arry vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cfix, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (.NOT.(lval.EQV.lfix)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      cval = MERGE('.TRUE. ','.FALSE.',lval)
      cfix = MERGE('.TRUE. ','.FALSE.',lfix)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for element '//TRIM(cindx(1))//&
                 & ' of logical array '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of logical array '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of logical array '//TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//','&
                 & //TRIM(cindx(4))//') of logical array '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      END SELECT
      WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cfix)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_value_arr_log

!========================================================================

SUBROUTINE error_value_arr_struc_int(ival,arrname,compname,ifix,ndims,indx)
!************************************************************************
!
! *error_value_arr_struc_int* Checks whether integer component 'ival' of a
!                             structure array equals 'ifix'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
INTEGER, INTENT(IN) :: ival, ifix, ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Component of derived type array
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*ifix*      INTEGER Integer value to which the array component must be equal
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Arry vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cfix, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (ival.NE.ifix) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; WRITE(cfix,'(I12)') ifix
      cval = ADJUSTL(cval); cfix = ADJUSTL(cfix)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array structure '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array '//TRIM(arrname)//': '//&
                 & TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '//&
                 & TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//TRIM(cindx(4))//&
                 & ') of array '//TRIM(arrname)//': '//TRIM(cval)
      END SELECT
      WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cfix)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_value_arr_struc_int

!========================================================================

SUBROUTINE error_value_arr_struc_log(lval,arrname,compname,lfix,ndims,indx)
!************************************************************************
!
! *error_value_arr_struc_log* Checks whether logical component 'lval' of a
!                             structure array equals 'lfix'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
LOGICAL, INTENT(IN) :: lval, lfix
INTEGER, INTENT(IN) ::  ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*lval*      LOGICAL Component of derived type array
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*lfix*      LOGICAL Logical value to which the array component must be equal
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Arry vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cfix, cval
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (.NOT.lval.EQV.lfix) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      cval = MERGE('.TRUE. ','.FALSE.',lval)
      cfix = MERGE('.TRUE. ','.FALSE.',lfix)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (ioerr,'(A)') 'Invalid value for logical component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array structure '//TRIM(arrname)//': '//TRIM(cval)
      CASE (2)
         WRITE (ioerr,'(A)') 'Invalid value for logical component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array structure '//TRIM(arrname)//&
                 & ': '//TRIM(cval)
      CASE (3)
         WRITE (ioerr,'(A)') 'Invalid value for logical component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of array structure '//TRIM(arrname)//': '//TRIM(cval)
      CASE DEFAULT
         WRITE (ioerr,'(A)') 'Invalid value for logical component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//TRIM(cindx(4))//&
                 & ') of array structure '//TRIM(arrname)//': '//TRIM(cval)
      END SELECT
      WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cfix)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_value_arr_struc_log

!========================================================================

SUBROUTINE error_value_var_char(cval,varname,cfix)
!************************************************************************
!
! *error_value_var_char* Checks whether character string 'cval' equals 'cfix'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: cval, varname, cfix 

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cval*      CHAR    Character variable
!*varname*   CHAR    Variable name
!*cfix*      CHAR    String to which the character variable must be equal
! 
!------------------------------------------------------------------------------
!

IF (TRIM(cval).NE.TRIM(cfix)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (ioerr,'(A)') 'Invalid value for character parameter '& 
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cfix)
   ENDIF
ENDIF

  
RETURN

END SUBROUTINE error_value_var_char

!========================================================================

SUBROUTINE error_value_var_int(ival,varname,ifix)
!************************************************************************
!
! *error_value_var_int* Checks whether integer parameter 'ival' equals 'ifix'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: ifix, ival

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Integer variable
!*varname*   CHAR    Variable name
!*ifix*      INTEGER Integer value to which the variable must be equal
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cfix, cval


IF (ival.NE.ifix) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') ival; WRITE(cfix,'(I12)') ifix
      cval = ADJUSTL(cval); cfix = ADJUSTL(cfix)
      WRITE (ioerr,'(A)') 'Invalid value for integer parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cfix)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_value_var_int

!========================================================================

SUBROUTINE error_value_var_log(lval,varname,lfix)
!************************************************************************
!
! *error_value_var_log* Checks whether logical parameter 'lval' equals 'lfix'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
LOGICAL, INTENT(IN) :: lval, lfix

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*lval*      LOGICAL Logical variable
!*varname*   CHAR    Variable name
!*lfix*      LOGICAL Logical value to which the variable must be equal
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=7) :: cfix, cval


IF (.NOT.(lval.EQV.lfix)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      cval = MERGE('.TRUE. ','.FALSE.',lval)
      cfix = MERGE('.TRUE. ','.FALSE.',lfix)
      WRITE (ioerr,'(A)') 'Invalid value for logical parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cfix)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_value_var_log

!========================================================================

SUBROUTINE error_value_var_real(rval,varname,rfix)
!************************************************************************
!
! *error_value_var_real* Checks whether real parameter 'rval' equals 'rfix'
!
! Author - Alexander Breugem
!
! Version - @(COHERENS)error_routines.F90  V2.9
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
REAL (KIND=kndreal), INTENT(IN) :: rval, rfix

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      REAL    Integer variable
!*varname*   CHAR    Variable name
!*rfix*      INTEGER Integer value to which the variable must be equal
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cfix, cval
REAL:: eps


eps = 2.0*SPACING(rval)
IF (ABS(rval-rfix).GT.eps) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(E12.4)') rval; WRITE(cfix,'(E12.4)') rfix
      cval = ADJUSTL(cval); cfix = ADJUSTL(cfix)
      WRITE (ioerr,'(A)') 'Invalid value for real parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cfix)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_value_var_real

!========================================================================

SUBROUTINE error_value_var_double(rval,varname,rfix)
!************************************************************************
!
! *error_value_var_double* Checks whether double precision parameter 'rval'
!                          equals 'rfix'
!
! Author - Alexander Breugem
!
! Version - @(COHERENS)error_routines.F90  V2.11
!
! Description - displays error message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
REAL (KIND=kndrlong), INTENT(IN) :: rval, rfix

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      DOUBLE  Integer variable
!*varname*   CHAR    Variable name
!*rfix*      INTEGER Integer value to which the variable must be equal
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cfix, cval
REAL:: eps


eps = 2.0*SPACING(rval)
IF (ABS(rval-rfix).GT.eps) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(E12.4)') rval; WRITE(cfix,'(E12.4)') rfix
      cval = ADJUSTL(cval); cfix = ADJUSTL(cfix)
      WRITE (ioerr,'(A)') 'Invalid value for real parameter '&
                        & //TRIM(varname)//': '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cfix)
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_value_var_double

!========================================================================

SUBROUTINE warning_reset_arr_char(cval,arrname,cset,ndims,indx)
!************************************************************************
!
! *warning_reset_arr_char* Resets the element 'cval' of a character array to
!                         'cset' if needed
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays warning message if needed
!             - cval and cset are assumed to have the same length
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname, cset
CHARACTER (LEN=*), intent(INOUT) :: cval
INTEGER, INTENT(IN) :: ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cval*      CHAR    Array element
!*arrname*   CHAR    Array name
!*cset*      INTEGER String to which the element is set
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (TRIM(cval).NE.TRIM(cset)) THEN
   IF (warnflag) THEN
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (iowarn,'(A)') 'WARNING: value of element '//TRIM(cindx(1))//&
                 & ' of character array '//TRIM(arrname)//' is set from '//&
                 & TRIM(cval)//' to '//TRIM(cset)
      CASE (2)
         WRITE (iowarn,'(A)') 'WARNING: value of element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of character array '//&
                 & TRIM(arrname)//' is set from '//TRIM(cval)//' to '//&
                 & TRIM(cset)
      CASE (3)
         WRITE (iowarn,'(A)') 'WARNING:  value of element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of character array '//TRIM(arrname)//' is set from '//&
                 & TRIM(cval)//' to '//TRIM(cset)
      CASE DEFAULT
         WRITE (iowarn,'(A)') 'WARNING: value of element ('//&
                 & TRIM(cindx(1))//','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ','//TRIM(cindx(4))//') of character array '//&
                 & TRIM(arrname)//' is set from '//TRIM(cval)//' to '//&
                 & TRIM(cset)
      END SELECT
   ENDIF
   cval = cset
ENDIF


RETURN

END SUBROUTINE warning_reset_arr_char

!========================================================================

SUBROUTINE warning_reset_arr_int(ival,arrname,iset,ndims,indx)
!************************************************************************
!
! *warning_reset_arr_int* Resets the element 'ival' of an integer array to
!                         'iset' if needed
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays warning message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN) :: iset, ndims
INTEGER, INTENT(INOUT) :: ival
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Array element
!*arrname*   CHAR    Array name
!*iset*      INTEGER Integer value to which the element is set
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cval, cset
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (ival.NE.iset) THEN
   IF (warnflag) THEN
      WRITE (cval,'(I12)') ival; WRITE(cset,'(I12)') iset
      cval = ADJUSTL(cval); cset = ADJUSTL(cset)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (iowarn,'(A)') 'WARNING: value of element '//TRIM(cindx(1))//&
                 & ' of integer array '//TRIM(arrname)//' is set from '//&
                 & TRIM(cval)//' to '//TRIM(cset)
      CASE (2)
         WRITE (iowarn,'(A)') 'WARNING: value of element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of integer array '//TRIM(arrname)//&
                 & ' is set from '//TRIM(cval)//' to '//TRIM(cset)
      CASE (3)
         WRITE (iowarn,'(A)') 'WARNING:  value of element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of integer array '//TRIM(arrname)//' is set from '//&
                 & TRIM(cval)//' to '//TRIM(cset)
      CASE DEFAULT
         WRITE (iowarn,'(A)') 'WARNING: value of element ('//&
                 & TRIM(cindx(1))//','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ','//TRIM(cindx(4))//') of integer array '//TRIM(arrname)//&
                 & ' is set from '//TRIM(cval)//' to '//TRIM(cset)
      END SELECT
   ENDIF
   ival = iset
ENDIF


RETURN

END SUBROUTINE warning_reset_arr_int

!========================================================================

SUBROUTINE warning_reset_arr_log(lval,arrname,lset,ndims,indx)
!************************************************************************
!
! *warning_reset_arr_log* Resets the element 'lval' of a logical array to
!                         'lset' if needed
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays warning message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: lset
CHARACTER (LEN=*), INTENT(IN) :: arrname
LOGICAL, INTENT(INOUT) :: lval
INTEGER, INTENT(IN) :: ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*lval*      LOGICAL Array element
!*arrname*   CHAR    Array name
!*lset*      LOGICAL Logical value to which the element is set
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cval, cset
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (.NOT.(lval.EQV.lset)) THEN
   IF (warnflag) THEN
      cval = MERGE ('.TRUE. ','.FALSE.',lval)
      cset = MERGE ('.TRUE. ','.FALSE.',lset)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (iowarn,'(A)') 'WARNING: value of element '//TRIM(cindx(1))//&
                 & ' of logical array '//TRIM(arrname)//' is set from '//&
                 & TRIM(cval)//' to '//TRIM(cset)
      CASE (2)
         WRITE (iowarn,'(A)') 'WARNING: value of element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//') of logical array '//TRIM(arrname)//&
                 & ' is set from '//TRIM(cval)//' to '//TRIM(cset)
      CASE (3)
         WRITE (iowarn,'(A)') 'WARNING:  value of element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//&
                 & ') of logical array '//TRIM(arrname)//' is set from '//&
                 & TRIM(cval)//' to '//TRIM(cset)
      CASE DEFAULT
         WRITE (iowarn,'(A)') 'WARNING: value of element ('//TRIM(cindx(1))//&
                 & ','//TRIM(cindx(2))//','//TRIM(cindx(3))//','//&
                 & TRIM(cindx(4))//') of logical array '//TRIM(arrname)//&
                 & ' is set from '//TRIM(cval)//' to '//TRIM(cset)
      END SELECT
   ENDIF
   lval = lset
ENDIF


RETURN

END SUBROUTINE warning_reset_arr_log

!========================================================================

SUBROUTINE warning_reset_arr_struc_char(cval,arrname,compname,cset,ndims,indx)
!************************************************************************
!
! *warning_reset_arr_struc_int* Resets a character component 'cval' of a
!                               structure array to 'cset' if needed
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays warning message if needed
!             - cval and cset are assumed to have the same length
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname, cset
CHARACTER (LEN=*), INTENT(INOUT) :: cval
INTEGER, INTENT(IN) :: ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cval*      CHAR    Component of derived type array
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*cset*      CHAR    String to which the component is set
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (TRIM(cval).NE.TRIM(cset)) THEN
   IF (warnflag) THEN
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (iowarn,'(A)') 'WARNING: value of character component '&
                 & //TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array '//TRIM(arrname)//' is set from '//TRIM(cval)//&
                 & ' to '//TRIM(cset)
      CASE (2)
         WRITE (iowarn,'(A)') 'WARNING: value of character component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','&
                 & //TRIM(cindx(2))//') of array '//TRIM(arrname)//&
                 & ' is set from '//TRIM(cval)//' to '//TRIM(cset)
      CASE (3)
         WRITE (iowarn,'(A)') 'WARNING: value of character component '&
                 & //TRIM(compname)//' and element ('//TRIM(cindx(1))//','&
                 & //TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '&
                 & //TRIM(arrname)//' is set from '//TRIM(cval)//' to '//&
                 & TRIM(cset)
      CASE DEFAULT
         WRITE (iowarn,'(A)') 'WARNING: value of character component '&
                 & //TRIM(compname)//' and element ('//TRIM(cindx(1))//','&
                 & //TRIM(cindx(2))//','//TRIM(cindx(3))//','//TRIM(cindx(4))&
                 & //') of '//TRIM(arrname)//' is set from '//TRIM(cval)//&
                 & ' to '//TRIM(cset)
      END SELECT
   ENDIF
   cval = cset
ENDIF


RETURN

END SUBROUTINE warning_reset_arr_struc_char

!========================================================================

SUBROUTINE warning_reset_arr_struc_int(ival,arrname,compname,iset,ndims,indx)
!************************************************************************
!
! *warning_reset_arr_struc_int* Resets an integer component 'ival' of a
!                               structure array to 'iset' if needed
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays warning message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
INTEGER, INTENT(IN) :: iset, ndims
INTEGER, INTENT(INOUT) :: ival
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Component of derived type array
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*iset*      INTEGER Integer value to which the component is set
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cval, cset
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (ival.NE.iset) THEN
   IF (warnflag) THEN
      WRITE (cval,'(I12)') ival; WRITE(cset,'(I12)') iset
      cval = ADJUSTL(cval); cset = ADJUSTL(cset)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (iowarn,'(A)') 'WARNING: value of integer component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' of array '//TRIM(arrname)//' is set from '//TRIM(cval)//&
                 & ' to '//TRIM(cset)
      CASE (2)
         WRITE (iowarn,'(A)') 'WARNING: value of integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array '//TRIM(arrname)//&
                 & ' is set from '//TRIM(cval)//' to '//TRIM(cset)
      CASE (3)
         WRITE (iowarn,'(A)') 'WARNING: value of integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '//&
                 & TRIM(arrname)//' is set from '//TRIM(cval)//' to '//&
                 & TRIM(cset)
      CASE DEFAULT
         WRITE (iowarn,'(A)') 'WARNING: value of integer component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//&
                 & TRIM(cindx(4))//') of '//TRIM(arrname)//' is set from '//&
                 & TRIM(cval)//' to '//TRIM(cset)
      END SELECT
   ENDIF
   ival = iset
ENDIF


RETURN

END SUBROUTINE warning_reset_arr_struc_int

!========================================================================

SUBROUTINE warning_reset_arr_struc_log(lval,arrname,compname,lset,ndims,indx)
!************************************************************************
!
! *warning_reset_arr_struc_log* Resets a logical component 'lval' of a
!                               structure array to 'lset' if needed
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays warning message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: lset
LOGICAL, INTENT(INOUT) :: lval
CHARACTER (LEN=*), INTENT(IN) :: arrname, compname
INTEGER, INTENT(IN) :: ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: indx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*lval*      LOGICAL Component of derived type array
!*arrname*   CHAR    Array name
!*compname*  CHAR    Component name
!*lset*      LOGICAL Logical value to which the component is set
!*ndims*     INTEGER Array rank
!*indx*      INTEGER Array vector index
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cval, cset
INTEGER :: n
CHARACTER (LEN=12), DIMENSION(ndims) :: cindx


IF (.NOT.(lval.EQV.lset)) THEN
   IF (warnflag) THEN
      cval = MERGE ('.TRUE. ','.FALSE.',lval)
      cset = MERGE ('.TRUE. ','.FALSE.',lset)
      n_110: DO n=1,ndims
         WRITE (cindx(n),'(I12)') indx(n); cindx(n) = ADJUSTL(cindx(n))
      ENDDO n_110
      SELECT CASE (ndims)
      CASE (1)
         WRITE (iowarn,'(A)') 'WARNING: value of logical component '//&
                 & TRIM(compname)//' and element '//TRIM(cindx(1))//&
                 & ' array '//TRIM(arrname)//' is set from '//TRIM(cval)//&
                 & ' to '//TRIM(cset)
      CASE (2)
         WRITE (iowarn,'(A)') 'WARNING: value of logical component '//&
                 & TRIM(compname)// ' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//') of array '//TRIM(arrname)//&
                 & ' is set from '//TRIM(cval)//' to '//TRIM(cset)
      CASE (3)
         WRITE (iowarn,'(A)') 'WARNING: value of logical component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//') of array '//&
                 & TRIM(arrname)//' is set from '//TRIM(cval)//' to '//&
                 & TRIM(cset)
      CASE DEFAULT
         WRITE (iowarn,'(A)') 'WARNING: value of logical component '//&
                 & TRIM(compname)//' and element ('//TRIM(cindx(1))//','//&
                 & TRIM(cindx(2))//','//TRIM(cindx(3))//','//TRIM(cindx(4))//&
                 & ') of array '//TRIM(arrname)//' is set from '//TRIM(cval)//&
                 & ' to '//TRIM(cset)
      END SELECT
   ENDIF
   lval = lset
ENDIF


RETURN

END SUBROUTINE warning_reset_arr_struc_log

!========================================================================

SUBROUTINE warning_reset_var_char(cval,varname,cset)
!************************************************************************
!
! *warning_reset_var_int* Resets character parameter 'cval' to 'cset' if needed
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays warning message if needed
!             - cval and cset are assumed to have the same length
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: cset, varname
CHARACTER (LEN=*), INTENT(INOUT) :: cval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cval*      INTEGER Character variable
!*varname*   CHAR    Variable name
!*cset*      INTEGER Character string to which the variable is set
! 
!------------------------------------------------------------------------------
!


IF (cval.NE.cset) THEN
   IF (warnflag) THEN
      WRITE (iowarn,'(6A)') 'WARNING: value of character parameter ',&
                           & TRIM(varname), ' is set from ', TRIM(cval),&
                           &' to ', TRIM(cset)
   ENDIF
   cval = cset
ENDIF


RETURN

END SUBROUTINE warning_reset_var_char

!========================================================================

SUBROUTINE warning_reset_var_int(ival,varname,iset)
!************************************************************************
!
! *warning_reset_var_int* Resets integer parameter 'ival' to 'iset' if needed
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays warning message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: iset
INTEGER, INTENT(INOUT) :: ival

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ival*      INTEGER Integer variable
!*varname*   CHAR    Variable name
!*iset*      INTEGER Integer value to which the variable is set
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cval, cset


IF (ival.NE.iset) THEN
   IF (warnflag) THEN
      WRITE (cval,'(I12)') ival; WRITE(cset,'(I12)') iset
      cval = ADJUSTL(cval); cset = ADJUSTL(cset)
      WRITE (iowarn,'(6A)') 'WARNING: value of integer parameter ',&
                           & TRIM(varname), ' is set from ', TRIM(cval),&
                           &' to ', TRIM(cset)
   ENDIF
   ival = iset
ENDIF


RETURN

END SUBROUTINE warning_reset_var_int

!========================================================================

SUBROUTINE warning_reset_var_log(lval,varname,lset)
!************************************************************************
!
! *warning_reset_var_log* Resets logical parameter 'lval' to 'lset' if needed
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays warning message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: lset
LOGICAL, INTENT(INOUT) :: lval
CHARACTER (LEN=*), INTENT(IN) :: varname

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*lval*      LOGICAL Logical variable
!*varname*   CHAR    Variable name
!*lset*      LOGICAL Logical value to which the variable is set
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=7) :: cval, cset


IF (.NOT.(lval.EQV.lset)) THEN
   IF (warnflag) THEN
      cval = MERGE ('.TRUE. ','.FALSE.',lval)
      cset = MERGE ('.TRUE. ','.FALSE.',lset)
      WRITE (iowarn,'(6A)') 'WARNING: value of logical parameter ',&
                           & TRIM(varname),' is set from ', TRIM(cval),&
                           &' to ', TRIM(cset)
   ENDIF
   lval = lset
ENDIF


RETURN

END SUBROUTINE warning_reset_var_log

!========================================================================

SUBROUTINE warning_reset_var_real(rval,varname,rset)
!************************************************************************
!
! *warning_reset_var_real* Resets real parameter 'rval' to 'rset'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays warning message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
REAL, INTENT(IN) :: rset
REAL (KIND=kndreal), INTENT(INOUT) :: rval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      REAL    Real variable
!*varname*   CHAR    Variable name
!*rset*      REAL    Real value to which the variable is set
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cval, cset


IF (rval.NE.rset) THEN
   IF (warnflag) THEN
      WRITE (cval,'(G15.7)') rval; WRITE(cset,'(G15.7)') rset
      WRITE (iowarn,'(6A)') 'WARNING: value of real parameter ',&
                           & TRIM(varname), ' is set from ', TRIM(cval),&
                           &' to ', TRIM(cset)
   ENDIF
   rval = rset
ENDIF


RETURN

END SUBROUTINE warning_reset_var_real

!========================================================================

SUBROUTINE warning_reset_var_double(rval,varname,rset)
!************************************************************************
!
! *warning_reset_var_double* Resets double precision parameter 'rval' to 'rset'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.11
!
! Description - displays warning message if needed
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
REAL, INTENT(IN) :: rset
REAL (KIND=kndrlong), INTENT(INOUT) :: rval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      DOUBLE  Real variable
!*varname*   CHAR    Variable name
!*rset*      REAL    Real value to which the variable is set
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cval, cset


IF (rval.NE.rset) THEN
   IF (warnflag) THEN
      WRITE (cval,'(G15.7)') rval; WRITE(cset,'(G15.7)') rset
      WRITE (iowarn,'(6A)') 'WARNING: value of real parameter ',&
                           & TRIM(varname), ' is set from ', TRIM(cval),&
                           &' to ', TRIM(cset)
   ENDIF
   rval = rset
ENDIF


RETURN

END SUBROUTINE warning_reset_var_double

!========================================================================

SUBROUTINE warning_ubound_var_real(rval,varname,maxval,matchmax)
!************************************************************************
!
! *warning_ubound_var_real* Checks whether real parameter 'rval' is lower than
!                           or equal to 'maxval' if 'matchmax' is .TRUE. or
!                           lower than 'maxval' if 'matchmax' is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.0
!
! Description - displays warning message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmax
CHARACTER (LEN=*), INTENT(IN) :: varname
REAL, INTENT(IN) :: rval, maxval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      REAL    Real variable
!*varname*   CHAR    Variable name
!*maxval*    REAL    Imposed upper bound
!*matchmax*  LOGICAL Maximum is allowed if .TRUE.
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cval, cmaxval


IF (matchmax) THEN
   IF (rval.GT.maxval) THEN
      IF (warnflag) THEN
         WRITE (cval,'(G15.7)') rval; WRITE (cmaxval,'(G15.7)') maxval
         cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
         WRITE (iowarn,'(A)') 'WARNING: suspect value for real parameter '&
                            & //TRIM(varname)//': '//TRIM(cval)
         WRITE (iowarn,'(A)') 'Should be lower than or equal to: '&
                            & //TRIM(cmaxval)
      ENDIF
   ENDIF
ELSEIF (rval.GE.maxval) THEN
   IF (warnflag) THEN
      WRITE (cval,'(G15.7)') rval; WRITE (cmaxval,'(G15.7)') maxval
      cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
      WRITE (iowarn,'(A)') 'WARNING: suspect value for real parameter '&
                         & //TRIM(varname)//': '//TRIM(cval)
      WRITE (iowarn,'(A)') 'Should be strictly greater than: '//TRIM(cmaxval)
   ENDIF
ENDIF


RETURN

END SUBROUTINE warning_ubound_var_real

!========================================================================

SUBROUTINE warning_ubound_var_double(rval,varname,maxval,matchmax)
!************************************************************************
!
! *warning_ubound_var_double* Checks whether double precision parameter 'rval'
!                             is lower than or equal to 'maxval' if 'matchmax'
!                             is .TRUE. or lower than 'maxval' if 'matchmax'
!                             is .FALSE.
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.11
!
! Description - displays warning message if needed
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmax
CHARACTER (LEN=*), INTENT(IN) :: varname
REAL (KIND=kndrlong), INTENT(IN) :: rval, maxval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*rval*      DOUBLE  Real variable
!*varname*   CHAR    Variable name
!*maxval*    REAL    Imposed upper bound
!*matchmax*  LOGICAL Maximum is allowed if .TRUE.
! 
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=15) :: cval, cmaxval


IF (matchmax) THEN
   IF (rval.GT.maxval) THEN
      IF (warnflag) THEN
         WRITE (cval,'(G15.7)') rval; WRITE (cmaxval,'(G15.7)') maxval
         cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
         WRITE (iowarn,'(A)') 'WARNING: suspect value for real parameter '&
                            & //TRIM(varname)//': '//TRIM(cval)
         WRITE (iowarn,'(A)') 'Should be lower than or equal to: '&
                            & //TRIM(cmaxval)
      ENDIF
   ENDIF
ELSEIF (rval.GE.maxval) THEN
   IF (warnflag) THEN
      WRITE (cval,'(G15.7)') rval; WRITE (cmaxval,'(G15.7)') maxval
      cval = ADJUSTL(cval); cmaxval = ADJUSTL(cmaxval)
      WRITE (iowarn,'(A)') 'WARNING: suspect value for real parameter '&
                         & //TRIM(varname)//': '//TRIM(cval)
      WRITE (iowarn,'(A)') 'Should be strictly greater than: '//TRIM(cmaxval)
   ENDIF
ENDIF


RETURN

END SUBROUTINE warning_ubound_var_double

!========================================================================

SUBROUTINE write_error_message(string)
!************************************************************************
!
! *write_error_message* write an error message
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.8
!
! Description -
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: string

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*meaasge*   CHAR    Error message
! 
!------------------------------------------------------------------------------
!


nerrs = nerrs + 1
IF (errchk.AND.nerrs.LE.maxerrors) THEN
   WRITE (ioerr,'(A)') TRIM(string)
ENDIF

END SUBROUTINE write_error_message


END MODULE error_routines
