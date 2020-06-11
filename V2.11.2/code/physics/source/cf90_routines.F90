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

MODULE cf90_routines
!************************************************************************
!
! *cf90_routines* netcdf routine library
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.11
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description -
!
! Generic routines - cf90_get_att_chars, cf90_get_att, cf90_get_var_chars,
!                    cf90_get_var, cf90_put_att_chars, cf90_put_att,
!                    cf90_put_var_chars, cf90_put_var
!
! Routines - cf90_abort, cf90_close, cf90_copy_att, cf90_create,
!            cf90_def_dim, cf90_def_var, cf90_del_att, cf90_enddef, cf90_error,
!            cf90_inq_attname, cf90_inq_dimid, cf90_inq_libvers,
!            cf90_inq_varid, cf90_inquire, cf90_inquire_attribute,
!            cf90_inquire_dimension, cf90_inquire_variable, cf90_open,
!            cf90_redef, cf90_rename_att, cf90_rename_dim, cf90_rename_var,
!            cf90_set_fill, cf90_sync
!
!************************************************************************
!

USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out
USE syspars
#ifdef CDF
   USE netcdf
#endif /*CDF*/
   
IMPLICIT NONE

!
!*Local variables
!

INTEGER :: npcc

INTERFACE cf90_get_att_chars
   MODULE PROCEDURE cf90_get_att_chars_0d, cf90_get_att_chars_1d
END INTERFACE

INTERFACE cf90_get_att
   MODULE PROCEDURE cf90_get_att_int_0d, cf90_get_att_int_1d, &
                  & cf90_get_att_real_0d, cf90_get_att_real_1d, &
                  & cf90_get_att_double_0d, cf90_get_att_double_1d
END INTERFACE cf90_get_att

INTERFACE cf90_get_var_chars
   MODULE PROCEDURE cf90_get_var_chars_0d, cf90_get_var_chars_1d, &
                  & cf90_get_var_chars_2d 
END INTERFACE cf90_get_var_chars

INTERFACE cf90_get_var
   MODULE PROCEDURE cf90_get_var_int_0d, cf90_get_var_int_1d, &
                  & cf90_get_var_int_2d, cf90_get_var_int_3d, &
                  & cf90_get_var_int_4d, cf90_get_var_log_0d, &
                  & cf90_get_var_log_1d, cf90_get_var_log_2d, &
                  & cf90_get_var_real_0d, cf90_get_var_real_1d, &
                  & cf90_get_var_real_2d, cf90_get_var_real_3d, &
                  & cf90_get_var_real_4d, cf90_get_var_double_0d, &
                  & cf90_get_var_double_1d, cf90_get_var_double_2d, &
                  & cf90_get_var_double_3d, cf90_get_var_double_4d
END INTERFACE

INTERFACE cf90_put_att_chars
   MODULE PROCEDURE cf90_put_att_chars_0d, cf90_put_att_chars_1d
END INTERFACE

INTERFACE cf90_put_att
   MODULE PROCEDURE cf90_put_att_int_0d, cf90_put_att_int_1d, &
                   & cf90_put_att_real_0d, cf90_put_att_real_1d, &
                   & cf90_put_att_double_0d, cf90_put_att_double_1d
END INTERFACE

INTERFACE cf90_put_var_chars
   MODULE PROCEDURE cf90_put_var_chars_0d, cf90_put_var_chars_1d, &
                  & cf90_put_var_chars_2d  
END INTERFACE cf90_put_var_chars

INTERFACE cf90_put_var 
   MODULE PROCEDURE cf90_put_var_int_0d, cf90_put_var_int_1d, &
                  & cf90_put_var_int_2d, cf90_put_var_int_3d, &
                  & cf90_put_var_int_4d, cf90_put_var_log_0d, &
                  & cf90_put_var_log_1d, cf90_put_var_log_2d, &
                  & cf90_put_var_real_0d, cf90_put_var_real_1d, &
                  & cf90_put_var_real_2d, cf90_put_var_real_3d, &
                  & cf90_put_var_real_4d, cf90_put_var_double_0d, &
                  & cf90_put_var_double_1d, cf90_put_var_double_2d, &
                  & cf90_put_var_double_3d, cf90_put_var_double_4d
   
END INTERFACE

CONTAINS

!========================================================================

SUBROUTINE cf90_abort(ncid)
!************************************************************************
!
! *cf90_abort* Abort netCDF
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description - 
!
! Reference -
!
! netCDF calls - NF90_abort
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: ncid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_abort'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_abort(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_abort')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_abort

!========================================================================

SUBROUTINE cf90_close(ncid)
!************************************************************************
!
! *cf90_close* Close an open netCDF dataset
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! Reference -
!
! netCDF calls - NF90_close
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: ncid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_close'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_close(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_close')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_close

!========================================================================

SUBROUTINE cf90_copy_att(ncid_in,varid_in,name,ncid_out,varid_out)
!************************************************************************
!
! *cf90_copy_att* Copy an attribute from one netCDF file to another
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! Reference -
!
! netCDF calls - NF90_copy_att
!
!************************************************************************
!
!*  Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid_in, ncid_out, varid_in, varid_out

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ncid_in*   INTEGER NetCDF ID of input file
!*varid_in*  INTEGER ID of associated variable in input file
!*name*      CHAR    Attribute name
!*ncid_out*  INTEGER NetCDF ID of output file
!*varid_out* INTEGER ID of associated variable in output file
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_copy_att'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_copy_att(ncid_in,varid_in,name,ncid_out,varid_out)
#endif /*CDF*/
IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_copy_att')


CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_copy_att

!========================================================================

SUBROUTINE cf90_create(path,cmode,ncid,initialsize,chunksize)
!************************************************************************
!
! *cf90_create* Open a new netCDF dataset
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! Reference -
!
! netCDF calls - NF90_create
!
!************************************************************************
!
!*  Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: path
INTEGER, INTENT(IN) :: cmode
INTEGER, INTENT(OUT) :: ncid
INTEGER, INTENT(IN), OPTIONAL :: initialsize
INTEGER, INTENT(INOUT), OPTIONAL :: chunksize

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*path*        CHAR    Name of new netCDF file
!*cmode*       INTEGER Creation mode
!*ncid*        INTEGER NetCDF ID
!*initialsize* INTEGER Initial file size in bytes
!*chunksize*   INTEGER Chunksize 
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: chunk, initsize


procname(pglev+1) = 'cf90_create'
CALL log_timer_in(npcc)

IF (PRESENT(initialsize)) THEN
   initsize = initialsize
ELSE
   initsize = 0
ENDIF
IF (PRESENT(chunksize)) THEN
   chunk = chunksize
ELSE
   chunk = sizehint_default_NF90
ENDIF

#ifdef CDF
errstat = NF90_create(path,cmode,ncid,initialsize=initsize,chunksize=chunk)
#else
ncid = -1
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_create')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_create

!========================================================================

SUBROUTINE cf90_def_dim(ncid,name,len,dimid)
!************************************************************************
!
! *cf90_def_dim* Adds a new dimension to an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! Reference -
!
! netCDF calls - NF90_def_dim
!
!************************************************************************
!
!*  Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: len, ncid
INTEGER, INTENT(OUT) :: dimid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*name*    CHAR    Dimension name
!*len*     INTEGER Dimension length
!*dimid*   INTEGER Returned dimension ID
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_def_dim'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(name))

#ifdef CDF
errstat = NF90_def_dim(ncid,name,len,dimid)
#else
dimid = 0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_def_dim')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_def_dim

!========================================================================

SUBROUTINE cf90_def_var(ncid,name,xtype,dimids,varid)
!************************************************************************
!
! *cf90_def_var* Adds a new variable to an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! Reference -
!
! netCDF calls - NF90_def_var
!
!************************************************************************
!
!*  Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, xtype
INTEGER, INTENT(OUT) :: varid
INTEGER, INTENT(IN), DIMENSION(:) :: dimids

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*name*    CHAR    Variable name
!*xtype*   INTEGER NetCDF data type of variable
!*dimids*  INTEGER Dimension IDs of variable
!*varid*   INTEGER Returned variable ID
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_def_var'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(name))

#ifdef CDF
errstat = NF90_def_var(ncid,name,xtype,dimids,varid)
#else
varid = 0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_def_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_def_var

!========================================================================

SUBROUTINE cf90_del_att(ncid,varid,name)
!************************************************************************
!
! *cf90_del_atts* Delete an attribute from a netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! Reference -
!
! netCDF calls - NF90_del_att
!
!************************************************************************
!
!*  Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, varid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_del_att'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_del_att(ncid,varid,name)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_del_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_del_att

!========================================================================

SUBROUTINE cf90_enddef(ncid,h_minfree,v_align,v_minfree,r_align)
!************************************************************************
!
! *cf90_enddef* Take an open netCDF dataset out of define mode
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! Reference -
!
! netCDF calls - NF90_enddef
!
!************************************************************************
!
!*  Arguments
!
INTEGER, INTENT(IN) :: ncid
INTEGER, INTENT(IN), OPTIONAL :: h_minfree, v_align, v_minfree, r_align

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*ncid*      INTEGER NetCDF ID
!*h_minfree* INTEGER Size in bytes of the pad at the end of the header
!*v_align*   INTEGER Alignment of the start of the data section for fixed size
!                    variables
!*v_minfree* INTEGER Size in bytes of the pad at the end of the data section
!                    for fixed size variables
!*r_align*   INTEGER Alignment of the start of the data section in case of an
!                    unlimited dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: hmin, valign, vmin, ralign


procname(pglev+1) = 'cf90_enddef'
CALL log_timer_in(npcc)

IF (PRESENT(h_minfree)) THEN
   hmin = h_minfree
ELSE
   hmin = 0
ENDIF
IF (PRESENT(v_minfree)) THEN
   vmin = v_minfree
ELSE
   vmin = 0
ENDIF
IF (PRESENT(v_align)) THEN
   valign = v_align
ELSE
   valign = 4
ENDIF
IF (PRESENT(r_align)) THEN
   ralign = r_align
ELSE
   ralign = 4
ENDIF

#ifdef CDF
errstat = NF90_enddef(ncid,h_minfree=hmin,v_align=valign,v_minfree=vmin,&
                    & r_align=ralign)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_enddef')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_enddef

!========================================================================

SUBROUTINE cf90_error(cf90name)
!************************************************************************
!
! *cf90_error* Report an error in a netCDF routine
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! Module calls - error_abort
!
! netCDF calls - NF90_strerror
!
!************************************************************************
!
USE switches
USE error_routines, ONLY: error_abort

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: cf90name

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cf90name*  CHAR    Name of the netCDF routine
! 
!------------------------------------------------------------------------------
!


#ifdef CDF
IF (errstat.NE.noerr_NF90) THEN
#else
IF (errstat.NE.0) THEN
#endif /*CDF*/
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      IF (iopt_CDF_abort.EQ.0) THEN
         WRITE (ioerr,'(A)') 'Error occurred in netCDF routine '//&
                            & TRIM(cf90name)
#ifdef CDF
         WRITE (ioerr,'(A)') TRIM(NF90_strerror(errstat))
#endif /*CDF*/
      ENDIF
   ENDIF
   IF (iopt_CDF_abort.EQ.1) CALL error_abort(cf90name,ierrno_CDF)
ENDIF


RETURN

END SUBROUTINE cf90_error

!========================================================================

SUBROUTINE cf90_get_att_chars_0d(ncid,varid,name,lenstr,chardat)
!************************************************************************
!
! *cf90_get_att_chars_0d* Get the string value of a scalar attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description - value is a fixed-length string
!
! netCDF calls - NF90_get_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, lenstr, varid
CHARACTER (LEN=lenstr), INTENT(OUT) :: chardat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*lenstr*  INTEGER Length of character string
!*chardat* CHAR    Returned attribute character value
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_get_att_chars_0d'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_get_att(ncid,varid,name,chardat)
#else
chardat = ''
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_att_chars_0d

!========================================================================

SUBROUTINE cf90_get_att_chars_1d(ncid,varid,name,lenstr,chardat)
!************************************************************************
!
! *cf90_get_att_chars_1d* Get the string value of a vector attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description - value is a vector of fixed-length strings
!
! netCDF calls - NF90_get_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, lenstr, varid
CHARACTER (LEN=lenstr), INTENT(OUT), DIMENSION(:) :: chardat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*lenstr*  INTEGER Lenght of character string vector
!*chardat* CHAR    Returned attribute character values
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=lenstr*SIZE(chardat)) :: chartmp
INTEGER :: isize, l1, nosize


procname(pglev+1) = 'cf90_get_att_chars_1d'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_get_att(ncid,varid,name,chartmp)
#else
chardat = ''
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_att')

nosize = SIZE(chardat)
isize_110: DO isize=1,nosize
   l1 = (isize-1)*lenstr
   chardat(isize) = chartmp(l1+1:l1+lenstr)
ENDDO isize_110

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_att_chars_1d

!========================================================================

SUBROUTINE cf90_get_att_int_0d(ncid,varid,name,intdat)
!************************************************************************
!
! *cf90_get_att_real_0d* Get the value of an integer scalar attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.11
!
! Description - 
!
! netCDF calls - NF90_get_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, varid
INTEGER, INTENT(OUT) :: intdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*intdat*  INTEGER Returned attribute integer value
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_get_att_int_0d'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_get_att(ncid,varid,name,intdat)
#else
intdat = 0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_att_int_0d

!========================================================================

SUBROUTINE cf90_get_att_int_1d(ncid,varid,name,intdat)
!************************************************************************
!
! *cf90_get_att_int_1d* Get the values of an integer vector attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.11
!
! Description - 
!
! netCDF calls - NF90_get_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, varid
INTEGER, INTENT(OUT), DIMENSION(:) :: intdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*intdat*  INTEGER Returned attribute integer values
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_get_att_int_1d'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_get_att(ncid,varid,name,intdat)
#else
intdat = 0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_att_int_1d

!========================================================================

SUBROUTINE cf90_get_att_real_0d(ncid,varid,name,realdat)
!************************************************************************
!
! *cf90_get_att_real_0d* Get the value of a single precision real scalar
!                        attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description - 
!
! netCDF calls - NF90_get_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, varid
REAL (KIND=kndreal), INTENT(OUT) :: realdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*realdat* REAL    Returned attribute real value
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_get_att_real_0d'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_get_att(ncid,varid,name,realdat)
#else
realdat = 0.0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_att_real_0d

!========================================================================

SUBROUTINE cf90_get_att_real_1d(ncid,varid,name,realdat)
!************************************************************************
!
! *cf90_get_att_real_1d* Get the values of single precision real vector
!                        attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description - 
!
! netCDF calls - NF90_get_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, varid
REAL (KIND=kndreal), INTENT(OUT), DIMENSION(:) :: realdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*realdat* REAL    Returned attribute real values
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_get_att_real_1d'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_get_att(ncid,varid,name,realdat)
#else
realdat = 0.0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_att_real_1d

!========================================================================

SUBROUTINE cf90_get_att_double_0d(ncid,varid,name,realdat)
!************************************************************************
!
! *cf90_get_att_double_0d* Get the value of a double precision real scalar
!                          attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description - 
!
! netCDF calls - NF90_get_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, varid
REAL (KIND=kndrlong), INTENT(OUT) :: realdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*realdat* DOUBLE  Returned attribute real value
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_get_att_double_0d'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_get_att(ncid,varid,name,realdat)
#else
realdat = 0.0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_att_double_0d

!========================================================================

SUBROUTINE cf90_get_att_double_1d(ncid,varid,name,realdat)
!************************************************************************
!
! *cf90_get_att_double_1d* Get the values of double precision real vector
!                          attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description - 
!
! netCDF calls - NF90_get_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, varid
REAL (KIND=kndrlong), INTENT(OUT), DIMENSION(:) :: realdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*realdat* DOUBLE  Returned attribute real values
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_get_att_double_1d'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_get_att(ncid,varid,name,realdat)
#else
realdat = 0.0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_att_double_1d

!========================================================================

SUBROUTINE cf90_get_var_chars_0d(ncid,varid,chardat,timerec)
!************************************************************************
!
! *cf90_get_var_chars_0d* Read a string variable from an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
CHARACTER (LEN=*), INTENT(OUT) :: chardat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*chardat* CHAR    Returned string
!*timerec* INTEGER Time record number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: l
INTEGER, DIMENSION(2) :: nstart


procname(pglev+1) = 'cf90_get_var_chars_0d'
CALL log_timer_in(npcc)

l = LEN(chardat)
nstart = (/1,timerec/)

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN 
   errstat = NF90_get_var(ncid,varid,chardat,start=nstart(1:1),count=(/l/))
ELSE
   errstat = NF90_get_var(ncid,varid,chardat,start=nstart,count=(/l,1/))
ENDIF
#else
chardat = ''
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_chars_0d

!========================================================================

SUBROUTINE cf90_get_var_chars_1d(ncid,varid,chardat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_chars_1d* Read a vector of string variables from an open
!                         netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(1) :: count, start, stride
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(:) :: chardat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*chardat* CHAR    Returned string
!*timerec* INTEGER Time record number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(3) :: ncount, nstart, nstride


procname(pglev+1) = 'cf90_get_var_chars_1d'
CALL log_timer_in(npcc)


!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/1,start,timerec/)
ELSE
   nstart = (/1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/LEN(chardat(1)),count,1/)
ELSE
   ncount = (/LEN(chardat(1)),SIZE(chardat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/1,stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN 
   errstat = NF90_get_var(ncid,varid,chardat,start=nstart(1:2),&
                        & count=ncount(1:2),stride=nstride(1:2))
ELSE
   errstat = NF90_get_var(ncid,varid,chardat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
chardat = ''
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_chars_1d

!========================================================================

SUBROUTINE cf90_get_var_chars_2d(ncid,varid,chardat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_chars_2d* Read a 2-D array of string variables from an open
!                         netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.8
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(:,:) :: chardat
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(2) :: count, start, stride

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*chardat* CHAR    Returned character 2-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for input buffer
!*count*   INTEGER Number of input data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(4) :: ncount, nstart, nstride


procname(pglev+1) = 'cf90_get_var_chars_2d'
CALL log_timer_in(npcc)


!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/1,start,timerec/)
ELSE
   nstart = (/1,1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/LEN(chardat(1,1)),count,1/)
ELSE
   ncount = (/LEN(chardat(1,1)),SHAPE(chardat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/1,stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN 
   errstat = NF90_get_var(ncid,varid,chardat,start=nstart(1:3),&
                        & count=ncount(1:3),stride=nstride(1:3))
ELSE
   errstat = NF90_get_var(ncid,varid,chardat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
chardat = ''
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_chars_2d

!========================================================================

SUBROUTINE cf90_get_var_int_0d(ncid,varid,intdat,timerec)
!************************************************************************
!
! *cf90_get_var_int_0d* Read an integer scalar variable from an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(OUT) :: intdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*intdat*  INTEGER Returned integer scalar
!*timerec* INTEGER Time record number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(1) :: nstart


procname(pglev+1) = 'cf90_get_var_int_0d'
CALL log_timer_in(npcc)

nstart = timerec

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,intdat)
ELSE
   errstat = NF90_get_var(ncid,varid,intdat,start=nstart)
ENDIF
#else
intdat = 0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_int_0d

!========================================================================

SUBROUTINE cf90_get_var_int_1d(ncid,varid,intdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_int_1d* Read a vector of integer data from an open netCDF file
!
! Author - Patrick Luyten
!
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(1) :: count, start, stride
INTEGER, INTENT(OUT), DIMENSION(:) :: intdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*intdat*  INTEGER Returned integer vector array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start index for input buffer
!*count*   INTEGER Number of input data
!*stride*  INTEGER Sampling interval
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(2) :: ncount, nstart, nstride


procname(pglev+1) = 'cf90_get_var_int_1d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(intdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,intdat,start=nstart(1:1),&
                        & count=ncount(1:1),stride=nstride(1:1))
ELSE
   errstat = NF90_get_var(ncid,varid,intdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
intdat = 0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_int_1d

!========================================================================

SUBROUTINE cf90_get_var_int_2d(ncid,varid,intdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_int_2d* Read a 2-D array of integer data from an open netCDF
!                       file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(2) :: count, start, stride
INTEGER, INTENT(OUT), DIMENSION(:,:) :: intdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*intdat*  INTEGER Returned integer 2-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for input buffer
!*count*   INTEGER Number of input data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(3) :: ncount, nstart, nstride


procname(pglev+1) = 'cf90_get_var_int_2d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(intdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,intdat,start=nstart(1:2),&
                        & count=ncount(1:2),stride=nstride(1:2))
ELSE
   errstat = NF90_get_var(ncid,varid,intdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
intdat = 0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_int_2d

!========================================================================

SUBROUTINE cf90_get_var_int_3d(ncid,varid,intdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_int_3d* Read a 3-D array of integer data from an open netCDF
!                       file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(3) :: count, start, stride
INTEGER, INTENT(OUT), DIMENSION(:,:,:) :: intdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*intdat*  INTEGER Returned integer 3-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for input buffer
!*count*   INTEGER Number of input data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(4) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_get_var_int_3d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(intdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,intdat,start=nstart(1:3),&
                        & count=ncount(1:3),stride=nstride(1:3))
ELSE
   errstat = NF90_get_var(ncid,varid,intdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
intdat = 0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_int_3d

!========================================================================

SUBROUTINE cf90_get_var_int_4d(ncid,varid,intdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_int_4d* Read a 4-D array of integer data from an open netCDF
!                       file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(4) :: count, start, stride
INTEGER, INTENT(OUT), DIMENSION(:,:,:,:) :: intdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*intdat*  INTEGER Returned integer 4-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for input buffer
!*count*   INTEGER Number of input data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(5) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_get_var_int_4d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(intdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,intdat,start=nstart(1:4),&
                        & count=ncount(1:4),stride=nstride(1:4))
ELSE
   errstat = NF90_get_var(ncid,varid,intdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
intdat = 0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_int_4d

!========================================================================

SUBROUTINE cf90_get_var_log_0d(ncid,varid,logdat,timerec)
!************************************************************************
!
! *cf90_get_var_log_0d* Read a logical scalar variable from an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
LOGICAL, INTENT(OUT) :: logdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*logdat*  LOGICAL Returned logical scalar
!*timerec* INTEGER Time record number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: intdat
INTEGER, DIMENSION(1) :: nstart


procname(pglev+1) = 'cf90_get_var_int_0d'
CALL log_timer_in(npcc)

nstart = timerec

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,intdat)
ELSE
   errstat = NF90_get_var(ncid,varid,intdat,start=nstart)
ENDIF
#else
intdat = 0
#endif /*CDF*/

logdat = intdat.EQ.1

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_log_0d

!========================================================================

SUBROUTINE cf90_get_var_log_1d(ncid,varid,logdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_log_1d* Read a vector of logical data from an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.11
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(1) :: count, start, stride
LOGICAL, INTENT(OUT), DIMENSION(:) :: logdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*logdat*  LOGICAL Returned logical vector array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start index for input buffer
!*count*   INTEGER Number of input data
!*stride*  INTEGER Sampling interval
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(2) :: ncount, nstart, nstride
INTEGER, DIMENSION(SIZE(logdat)) :: intdat


procname(pglev+1) = 'cf90_get_var_log_1d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(logdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,intdat,start=nstart(1:1),&
                        & count=ncount(1:1),stride=nstride(1:1))
ELSE
   errstat = NF90_get_var(ncid,varid,intdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
intdat = 0
#endif /*CDF*/

logdat = intdat.EQ.1

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_log_1d

!========================================================================

SUBROUTINE cf90_get_var_log_2d(ncid,varid,logdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_log_2d* Read a 2-D array of logical data from an open netCDF
!                       file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.11
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(2) :: count, start, stride
LOGICAL, INTENT(OUT), DIMENSION(:,:) :: logdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*logdat*  LOGICAL Returned logical 2-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for input buffer
!*count*   INTEGER Number of input data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(3) :: ncount, nstart, nstride
INTEGER, DIMENSION(LBOUND(logdat,DIM=1):UBOUND(logdat,DIM=1),&
                 & LBOUND(logdat,DIM=2):UBOUND(logdat,DIM=2)) :: intdat

procname(pglev+1) = 'cf90_get_var_log_2d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(logdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,intdat,start=nstart(1:2),&
                        & count=ncount(1:2),stride=nstride(1:2))
ELSE
   errstat = NF90_get_var(ncid,varid,intdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
intdat = 0
#endif /*CDF*/

WHERE (intdat.EQ.1)
   logdat = .TRUE.
ELSEWHERE
   logdat = .FALSE.
END WHERE

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_log_2d

!========================================================================

SUBROUTINE cf90_get_var_real_0d(ncid,varid,realdat,timerec)
!************************************************************************
!
! *cf90_get_var_real_0d* Read a single precision real scalar variable from an
!                        open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
REAL (KIND=kndreal), INTENT(OUT) :: realdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Returned real scalar
!*timerec* INTEGER Time record number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(1) :: nstart


procname(pglev+1) = 'cf90_get_var_real_0d'
CALL log_timer_in(npcc)

nstart = timerec

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,realdat)
ELSE
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart)
ENDIF
#else
realdat = 0.0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_real_0d

!========================================================================

SUBROUTINE cf90_get_var_real_1d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_real_1d* Read a vector of single precision real data from an
!                        open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(1) :: count, start, stride
REAL (KIND=kndreal), INTENT(OUT), DIMENSION(:) :: realdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Returned real vector array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start index for input buffer
!*count*   INTEGER Number of input data
!*stride*  INTEGER Sampling interval
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(2) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_get_var_real_1d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart(1:1),&
                        & count=ncount(1:1),stride=nstride(1:1))
ELSE
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
realdat = 0.0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_real_1d

!========================================================================

SUBROUTINE cf90_get_var_real_2d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_real_2d* Read a 2-D array of single precision real data from an
!                        open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(2) :: count, start, stride
REAL (KIND=kndreal), INTENT(OUT), DIMENSION(:,:) :: realdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Returned real 2-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for input buffer
!*count*   INTEGER Number of input data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(3) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_get_var_real_2d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart(1:2),&
                        & count=ncount(1:2),stride=nstride(1:2))
ELSE
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
realdat = 0.0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_real_2d

!========================================================================

SUBROUTINE cf90_get_var_real_3d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_real_3d* Read a 3-D array of single precision real data from an
!                        open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!
  
USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(3) :: count, start, stride
REAL (KIND=kndreal), INTENT(OUT), DIMENSION(:,:,:) :: realdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Returned real 3-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for input buffer
!*count*   INTEGER Number of input data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(4) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_get_var_real_3d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart(1:3),&
                        & count=ncount(1:3),stride=nstride(1:3))
ELSE
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
realdat = 0.0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_real_3d

!========================================================================

SUBROUTINE cf90_get_var_real_4d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_real_4d* Read a 4-D array of single precision real data from an
!                        open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!********************* ***************************************************
!
  
USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(4) :: count, start, stride
REAL (KIND=kndreal), INTENT(OUT), DIMENSION(:,:,:,:) :: realdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Vairalbe ID
!*realdat* REAL    Returned real 4-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for input buffer
!*count*   INTEGER Number of input data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(5) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_get_var_real_4d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart(1:4),&
                        & count=ncount(1:4),stride=nstride(1:4))
ELSE
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
realdat = 0.0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_real_4d

!========================================================================

SUBROUTINE cf90_get_var_double_0d(ncid,varid,realdat,timerec)
!************************************************************************
!
! *cf90_get_var_real_0d* Read a double precision real scalar variable from an
!                        open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
REAL (KIND=kndrlong), INTENT(OUT) :: realdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Returned real scalar
!*timerec* INTEGER Time record number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(1) :: nstart


procname(pglev+1) = 'cf90_get_var_double_0d'
CALL log_timer_in(npcc)

nstart = timerec

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,realdat)
ELSE
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart)
ENDIF
#else
realdat = 0.0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_double_0d

!========================================================================

SUBROUTINE cf90_get_var_double_1d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_double_1d* Read a vector of double precision real data from an
!                          open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(1) :: count, start, stride
REAL (KIND=kndrlong), INTENT(OUT), DIMENSION(:) :: realdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Returned real vector array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start index for input buffer
!*count*   INTEGER Number of input data
!*stride*  INTEGER Sampling interval
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(2) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_get_var_double_1d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart(1:1),&
                        & count=ncount(1:1),stride=nstride(1:1))
ELSE
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
realdat = 0.0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_double_1d

!========================================================================

SUBROUTINE cf90_get_var_double_2d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_double_2d* Read a 2-D array of double precision real data from
!                          an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(2) :: count, start, stride
REAL (KIND=kndrlong), INTENT(OUT), DIMENSION(:,:) :: realdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Returned real 2-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for input buffer
!*count*   INTEGER Number of input data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(3) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_get_var_double_2d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart(1:2),&
                        & count=ncount(1:2),stride=nstride(1:2))
ELSE
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
realdat = 0.0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_double_2d

!========================================================================

SUBROUTINE cf90_get_var_double_3d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_double_3d* Read a 3-D array of double precision real data from
!                          an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(3) :: count, start, stride
REAL (KIND=kndrlong), INTENT(OUT), DIMENSION(:,:,:) :: realdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Returned real 3-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for input buffer
!*count*   INTEGER Number of input data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(4) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_get_var_double_3d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart(1:3),&
                        & count=ncount(1:3),stride=nstride(1:3))
ELSE
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
realdat = 0.0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_double_3d

!========================================================================

SUBROUTINE cf90_get_var_double_4d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_get_var_double_4d* Read a 4-D array of double precision real data from
!                          an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_get_var
!
!********************* ***************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(4) :: count, start, stride
REAL (KIND=kndrlong), INTENT(OUT), DIMENSION(:,:,:,:) :: realdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Vairalbe ID
!*realdat* REAL    Returned real 4-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for input buffer
!*count*   INTEGER Number of input data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(5) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_get_var_double_4d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
IF (timerec.EQ.0) THEN
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart(1:4),&
                        & count=ncount(1:4),stride=nstride(1:4))
ELSE
   errstat = NF90_get_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
#else
realdat = 0.0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_get_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_get_var_double_4d

!========================================================================

SUBROUTINE cf90_inq_attname(ncid,varid,attnum,name,ierr)
!************************************************************************
!
! *cf90_inq_attname* Returns the name of an attribute given its variable ID and
!                    number
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - NF90_inq_attname
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT) :: name
INTEGER, INTENT(IN) :: attnum, ncid, varid
INTEGER, INTENT(OUT), OPTIONAL :: ierr

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*attnum*  INTEGER Attribute number
!*name*    CHAR    Attribute name
!*ierr*    INTEGER If PRESENT, the program does not abort and the error code
!                  number is passed to the calling program
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_inq_attname'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_inq_attname(ncid,varid,attnum,name)
#else
name = ''
errstat = -1
#endif /*CDF*/

IF (PRESENT(ierr)) THEN
   ierr = errstat
ELSEIF (errstat.NE.noerr_NF90) THEN
   CALL cf90_error('NF90_inq_attname')
ENDIF

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_inq_attname

!========================================================================

SUBROUTINE cf90_inq_dimid(ncid,name,dimid)
!************************************************************************
!
! *cf90_inq_dimid* Returns the dimension ID given the dimension name
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! netCDF calls - NF90_inq_dimid
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid
INTEGER, INTENT(OUT) :: dimid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*name*    CHAR    Dimension name
!*dimid*   INTEGER Returned dimension ID
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_inq_dimid'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_inq_dimid(ncid,name,dimid)
#else
dimid = 0
#endif /*CDF*/
IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_inq_dimid')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_inq_dimid

!========================================================================

FUNCTION cf90_inq_libvers()
!************************************************************************
!
! *cf90_inq_libvers* Returns the version of netCDF library
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! netCDF calls - NF90_inq_libvers
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=80) :: cf90_inq_libvers

!
! Name              Type   Purpose
!------------------------------------------------------------------------------
!*cf90_inq_libvers* CHAR   NetCDF version
!
!------------------------------------------------------------------------------
!


#ifdef CDF
cf90_inq_libvers = NF90_inq_libvers()
#else
cf90_inq_libvers = ''
#endif /*CDF*/


RETURN

END FUNCTION cf90_inq_libvers

!========================================================================

SUBROUTINE cf90_inq_varid(ncid,name,varid)
!************************************************************************
!
! *cf90_inq_varid* Returns the ID of a netCDF variable given its name
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! netCDF calls - NF90_inq_varid
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid
INTEGER, INTENT(OUT) :: varid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*name*    CHAR    Variable name
!*varid*   INTEGER Returned variable ID
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_inq_varid'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_inq_varid(ncid,name,varid)
#else
varid = 0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_inq_varid')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_inq_varid

!========================================================================

SUBROUTINE cf90_inquire(ncid,ndimensions,nvariables,nattributes,&
                      & unlimiteddimid,abort)
!************************************************************************
!
! *cf90_inquire* Returns information about an open netCDF file given its
!                netcdf ID
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.11
!
! Description -
!
! netCDF calls - NF90_inquire
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: abort
INTEGER, INTENT(IN) :: ncid
INTEGER, INTENT(OUT), OPTIONAL :: ndimensions, nvariables, nattributes, &
                                & unlimiteddimid

!
! Name            Type     Purpose
!------------------------------------------------------------------------------
!*ncid*           INTEGER NetCDF ID
!*ndimensions*    INTEGER Returned number of dimensions
!*nvariables*     INTEGER Returned number of variables

!*nattributes*    INTEGER Returned number of global attributes
!*unlimiteddimid* INTEGER Returned ID of the unlimited dimension (if any)
!*abort*          LOGICAL If PRESENT and .TRUE. or not PRESENT, an error message is
!                         written, abortion of the program then depends on the
!                         setting of iopt_CDF_abort.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: abortx


procname(pglev+1) = 'cf90_inquire'
CALL log_timer_in(npcc)

#ifdef CDF
IF (PRESENT(abort)) THEN
   abortx = abort
ELSE
   abortx = .TRUE.
ENDIF
errstat = NF90_inquire(ncid,ndimensions,nvariables,nattributes,unlimiteddimid)
#else
IF (PRESENT(ndimensions)) ndimensions = 0
IF (PRESENT(nvariables)) nvariables = 0
IF (PRESENT(nattributes)) nattributes = 0
IF (PRESENT(unlimiteddimid)) unlimiteddimid = 0
IF (abortx.AND.errstat.NE.noerr_NF90) CALL cf90_error('NF90_inquire')
#endif /*CDF*/

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_inquire

!========================================================================

SUBROUTINE cf90_inquire_attribute(ncid,varid,name,xtype,len,attnum,ierr)
!************************************************************************
!
! *cf90_inquire_attribute* Returns information about a netCDF attribute
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - NF90_inquire_attribute
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, varid
INTEGER, INTENT(OUT), OPTIONAL :: attnum, len, xtype
INTEGER, INTENT(OUT), OPTIONAL :: ierr

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of associated variable
!*name*    CHAR    Attribute name
!*xtype*   INTEGER Returned netCDF data type of attribute
!*len*     INTEGER Returned number of values stored in attribute
!*attnum*  INTEGER Returned attribute number
!*ierr*    INTEGER If PRESENT, the program does not abort and the error code
!                  number is passed to the calling program
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_inquire_attribute'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_inquire_attribute(ncid,varid,name,xtype,len,attnum)
#else
IF (PRESENT(attnum)) attnum = 0
IF (PRESENT(len)) len = 0
IF (PRESENT(xtype)) xtype = 0
#endif /*CDF*/

IF (PRESENT(ierr)) THEN
   ierr = errstat
ELSEIF (errstat.NE.noerr_NF90) THEN
   CALL cf90_error('NF90_inquire_attribute')
ENDIF

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_inquire_attribute

!========================================================================

SUBROUTINE cf90_inquire_dimension(ncid,dimid,name,len)
!************************************************************************
!
! *cf90_inquire_dimension* Returns information about a netCDF dimension
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! netCDF calls - NF90_inquire_dimension
!
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: dimid, ncid
CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: name
INTEGER, INTENT(OUT), OPTIONAL :: len

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*dimid*   INTEGER Dimension ID
!*name*    CHAR    Returned dimension name
!*len*     INTEGER Returned dimension length
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_inquire_dimension'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_inquire_dimension(ncid,dimid,name,len)
#else
IF (PRESENT(name)) name = ''
IF (PRESENT(len)) len = 0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_inquire_dimension')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_inquire_dimension

!========================================================================

SUBROUTINE cf90_inquire_variable(ncid,varid,name,xtype,ndims,dimids,natts)
!************************************************************************
!
! *cf90_inquire_variable* Returns information about a netCDF variable
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! netCDF calls - NF90_inquire_variable
!
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, varid
CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: name
INTEGER, INTENT(OUT), OPTIONAL :: natts, ndims, xtype
INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(:) :: dimids

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*name*    CHAR    Returned variable name
!*xtype*   INTEGER Returned netCDF datatype of variable
!*ndims*   INTEGER Returned variable rank
!*dimids*  INTEGER Returned IDs of the variable's dimensions
!*natts*   INTEGER Returned number of associated attributes
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_inquire_variable'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_inquire_variable(ncid,varid,name,xtype,ndims,dimids,natts)
#else
IF (PRESENT(name)) name = ''
IF (PRESENT(natts)) natts = 0
IF (PRESENT(ndims)) ndims = 0
IF (PRESENT(xtype)) xtype = 0
IF (PRESENT(dimids)) dimids = 0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_inquire_variable')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_inquire_variable

!========================================================================

SUBROUTINE cf90_open(path,mode,ncid,chunksize)
!************************************************************************
!
! *cf90_open* Open an existing netCDF dataset
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! Reference -
!
! netCDF calls - NF90_open
!
!************************************************************************
!
!*  Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: path
INTEGER, INTENT(IN) :: mode
INTEGER, INTENT(OUT) :: ncid
INTEGER, INTENT(INOUT), OPTIONAL :: chunksize

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*path*      CHAR    NetCDF file name
!*mode*      INTEGER Access mode
!*ncid*      INTEGER NetCDF ID
!*chunksize* INTEGER Chunksize
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: chunk


procname(pglev+1) = 'cf90_open'
CALL log_timer_in(npcc)

IF (PRESENT(chunksize)) THEN
   chunk = chunksize
ELSE
   chunk = sizehint_default_NF90
ENDIF

#ifdef CDF
errstat = NF90_open(path,mode,ncid,chunksize=chunk)
#else
ncid = -1
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_open')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_open

!========================================================================

SUBROUTINE cf90_put_att_chars_0d(ncid,varid,name,lenstr,chardat)
!************************************************************************
!
! *cf90_put_att_chars_0d* Add or change the string value of a scalar attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description - value is a fixed-length string
!
! netCDF calls - NF90_put_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, lenstr, varid
CHARACTER (LEN=lenstr), INTENT(IN) :: chardat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*lenstr*  INTEGER Length of character string
!*chardat* CHAR    Attribute character value
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_put_att_chars_0d'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(name))

#ifdef CDF
errstat = NF90_put_att(ncid,varid,name,chardat)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_att_chars_0d

!========================================================================

SUBROUTINE cf90_put_att_chars_1d(ncid,varid,name,lenstr,chardat)
!************************************************************************
!
! *cf90_put_att_chars_1d* Add or change the string value of a vector attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description - value is a vector of fixed-length strings
!
! netCDF calls - NF90_put_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, lenstr, varid
CHARACTER (LEN=lenstr), INTENT(IN), DIMENSION(:) :: chardat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*lenstr*  INTEGER Length of character string vector
!*chardat* CHAR    Attribute character values
!
!------------------------------------------------------------------------------
!
!
!*Local variables
!
CHARACTER (LEN=lenstr*SIZE(chardat)) :: chartmp
INTEGER :: isize, l1, nosize


procname(pglev+1) = 'cf90_put_att_chars_1d'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(name))

nosize = SIZE(chardat)
isize_110: DO isize=1,nosize
   l1 = (isize-1)*lenstr
   chartmp(l1+1:l1+lenstr) =  chardat(isize)
ENDDO isize_110

#ifdef CDF
errstat = NF90_put_att(ncid,varid,name,chartmp)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_att_chars_1d

!========================================================================

SUBROUTINE cf90_put_att_int_0d(ncid,varid,name,intdat)
!************************************************************************
!
! *cf90_put_att_int_0d* Add or change the integer value of an attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description - 
!
! netCDF calls - NF90_put_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: intdat, ncid, varid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*intdat*  INTEGER Attribute integer value
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_put_att_int_0d'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(name))

#ifdef CDF
errstat = NF90_put_att(ncid,varid,name,intdat)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_att_int_0d

!========================================================================

SUBROUTINE cf90_put_att_int_1d(ncid,varid,name,intdat)
!************************************************************************
!
! *cf90_put_att_int_1d* Add or change the integer values of a vector attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description - 
!
! netCDF calls - NF90_put_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, varid
INTEGER, INTENT(IN), DIMENSION(:) :: intdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*intdat*  INTEGER Aattribute integer values
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_put_att_int_1d'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(name))

#ifdef CDF
errstat = NF90_put_att(ncid,varid,name,intdat)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_att_int_1d

!========================================================================

SUBROUTINE cf90_put_att_real_0d(ncid,varid,name,realdat)
!************************************************************************
!
! *cf90_put_att_real_0d* Add or change the value of a single precision real
!                        scalar attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description - 
!
! netCDF calls - NF90_put_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, varid
REAL (KIND=kndreal), INTENT(IN) :: realdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*realdat* REAL    Attribute real value
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_put_att_real_0d'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(name))

#ifdef CDF
errstat = NF90_put_att(ncid,varid,name,realdat)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_att_real_0d

!========================================================================

SUBROUTINE cf90_put_att_real_1d(ncid,varid,name,realdat)
!************************************************************************
!
! *cf90_put_att_real_1d* Add or change the value of a single precisino real
!                        vector attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description - 
!
! netCDF calls - NF90_put_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, varid
REAL (KIND=kndreal), INTENT(IN), DIMENSION(:) :: realdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*realdat* REAL    Attribute real values
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_put_att_real_1d'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(name))

#ifdef CDF
errstat = NF90_put_att(ncid,varid,name,realdat)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_att_real_1d

!========================================================================

SUBROUTINE cf90_put_att_double_0d(ncid,varid,name,realdat)
!************************************************************************
!
! *cf90_put_att_double_0d* Add or change the value of a double precision
!                          scalar attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description - 
!
! netCDF calls - NF90_put_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, varid
REAL (KIND=kndrlong), INTENT(IN) :: realdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*realdat* DOUBLE  Attribute real value
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_put_att_double_0d'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(name))

#ifdef CDF
errstat = NF90_put_att(ncid,varid,name,realdat)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_att_double_0d

!========================================================================

SUBROUTINE cf90_put_att_double_1d(ncid,varid,name,realdat)
!************************************************************************
!
! *cf90_put_att_double_1d* Add or change the value of a double precision real
!                          vector attribute 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description - 
!
! netCDF calls - NF90_put_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: ncid, varid
REAL (KIND=kndrlong), INTENT(IN), DIMENSION(:) :: realdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*name*    CHAR    Attribute name
!*realdat* DOUBLE  Attribute real values
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_put_att_double_1d'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//TRIM(name))

#ifdef CDF
errstat = NF90_put_att(ncid,varid,name,realdat)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_att_double_1d

!========================================================================

SUBROUTINE cf90_put_var_chars_0d(ncid,varid,chardat,timerec)
!************************************************************************
!
! *cf90_put_var_chars_0d* Write a string variable to an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: chardat
INTEGER, INTENT(IN) :: ncid, timerec, varid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*chardat* CHAR    Output string
!*timerec* INTEGER Time record number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: l
INTEGER, DIMENSION(2) :: nstart


procname(pglev+1) = 'cf90_put_var_chars_0d'
CALL log_timer_in(npcc)

l = LEN(chardat)
nstart = (/1,timerec/)

#ifdef CDF
IF (timerec.EQ.0) THEN 
   errstat = NF90_put_var(ncid,varid,chardat,start=nstart(1:1),count=(/l/))
ELSE
   errstat = NF90_put_var(ncid,varid,chardat,start=nstart,count=(/l,1/))
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_chars_0d

!========================================================================

SUBROUTINE cf90_put_var_chars_1d(ncid,varid,chardat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_chars_1d* Write a vector of string variables to an open netCDF
!                         file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN), DIMENSION(:) :: chardat
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(1) :: count, start, stride
INTEGER, INTENT(IN) :: ncid, timerec, varid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*chardat* CHAR    Output string
!*timerec* INTEGER Time record number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(3) :: ncount, nstart, nstride


procname(pglev+1) = 'cf90_put_var_chars_1d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/1,start,timerec/)
ELSE
   nstart = (/1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/LEN(chardat(1)),count,1/)
ELSE
   ncount = (/LEN(chardat(1)),SIZE(chardat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/1,stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (timerec.EQ.0) THEN 
   errstat = NF90_put_var(ncid,varid,chardat,start=nstart(1:2),&
                        & count=ncount(1:2),stride=nstride(1:2))
ELSE
   errstat = NF90_put_var(ncid,varid,chardat,start=nstart,&
                        & count=ncount,stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_chars_1d

!========================================================================

SUBROUTINE cf90_put_var_chars_2d(ncid,varid,chardat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_chars_2d* Write a 2-D array of string variables to an open
!                         netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.8
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN), DIMENSION(:,:) :: chardat
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(2) :: count, start, stride
INTEGER, INTENT(IN) :: ncid, timerec, varid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*chardat* CHAR    Output string
!*timerec* INTEGER Time record number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(4) :: ncount, nstart, nstride


procname(pglev+1) = 'cf90_put_var_chars_2d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/1,start,timerec/)
ELSE
   nstart = (/1,1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/LEN(chardat(1,1)),count,1/)
ELSE
   ncount = (/LEN(chardat(1,1)),SHAPE(chardat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/1,stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (timerec.EQ.0) THEN 
   errstat = NF90_put_var(ncid,varid,chardat,start=nstart(1:3),&
                        & count=ncount(1:2),stride=nstride(1:3))
ELSE
   errstat = NF90_put_var(ncid,varid,chardat,start=nstart,&
                        & count=ncount,stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_chars_2d

!========================================================================

SUBROUTINE cf90_put_var_int_0d(ncid,varid,intdat,timerec)
!************************************************************************
!
! *cf90_put_var_int_0d* Write an integer scalar variable to an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: intdat, ncid, timerec, varid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*intdat*  INTEGER Output integer scalar
!*timerec* INTEGER Time record number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(1) :: nstart


procname(pglev+1) = 'cf90_put_var_int_0d'
CALL log_timer_in(npcc)

nstart = timerec

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,intdat)
ELSE
   errstat = NF90_put_var(ncid,varid,intdat,start=nstart)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_int_0d

!========================================================================

SUBROUTINE cf90_put_var_int_1d(ncid,varid,intdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_int_1d* Write a vector of integer data to an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(1) :: count, start, stride
INTEGER, INTENT(IN), DIMENSION(:) :: intdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*intdat*  INTEGER Output integer vector array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start index for output buffer
!*count*   INTEGER Number of output data
!*stride*  INTEGER Sampling interval
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(2) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_put_var_int_1d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(intdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,intdat,start=nstart(1:1),&
                        & count=ncount(1:1),stride=nstride(1:1))
ELSE
   errstat = NF90_put_var(ncid,varid,intdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_int_1d

!========================================================================

SUBROUTINE cf90_put_var_int_2d(ncid,varid,intdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_int_2d* Write a 2-D array of integer data to an open netCDF
!                       file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(2) :: count, start, stride
INTEGER, INTENT(IN), DIMENSION(:,:) :: intdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*intdat*  INTEGER Output integer 2-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for output buffer
!*count*   INTEGER Number of output data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(3) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_put_var_int_2d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(intdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,intdat,start=nstart(1:2),&
                        & count=ncount(1:2),stride=nstride(1:2))
ELSE
   errstat = NF90_put_var(ncid,varid,intdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_int_2d

!========================================================================

SUBROUTINE cf90_put_var_int_3d(ncid,varid,intdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_int_3d* Write a 3-D array of integer data to an open netCDF
!                       file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(3) :: count, start, stride
INTEGER, INTENT(IN), DIMENSION(:,:,:) :: intdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*intdat*  INTEGER Output integer 3-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for output buffer
!*count*   INTEGER Number of output data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(4) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_put_var_int_3d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(intdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,intdat,start=nstart(1:3),&
                        & count=ncount(1:3),stride=nstride(1:3))
ELSE
   errstat = NF90_put_var(ncid,varid,intdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_int_3d

!========================================================================

SUBROUTINE cf90_put_var_int_4d(ncid,varid,intdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_int_4d* Write a 4-D array of integer data to an open netCDF
!                       file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.7.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(4) :: count, start, stride
INTEGER, INTENT(IN), DIMENSION(:,:,:,:) :: intdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*intdat*  INTEGER Output integer 4-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for output buffer
!*count*   INTEGER Number of output data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(5) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_put_var_int_4d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(intdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,intdat,start=nstart(1:4),&
                        & count=ncount(1:4),stride=nstride(1:4))
ELSE
   errstat = NF90_put_var(ncid,varid,intdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_int_4d

!========================================================================

SUBROUTINE cf90_put_var_log_0d(ncid,varid,logdat,timerec)
!************************************************************************
!
! *cf90_put_var_log_0d* Write an logical scalar variable to an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.11
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
LOGICAL, INTENT(IN) :: logdat
INTEGER, INTENT(IN) :: ncid, timerec, varid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*logdat*  LOGICAL Output logical scalar
!*timerec* INTEGER Time record number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: intdat
INTEGER, DIMENSION(1) :: nstart


procname(pglev+1) = 'cf90_put_var_log_0d'
CALL log_timer_in(npcc)

nstart = timerec

#ifdef CDF
intdat = MERGE (1,0,logdat)
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,intdat)
ELSE
   errstat = NF90_put_var(ncid,varid,intdat,start=nstart)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_log_0d

!========================================================================

SUBROUTINE cf90_put_var_log_1d(ncid,varid,logdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_log_1d* Write a vector of logical data to an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.11
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(1) :: count, start, stride
LOGICAL, INTENT(IN), DIMENSION(:) :: logdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*logdat*  LOGICAL Output logical vector array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start index for output buffer
!*count*   INTEGER Number of output data
!*stride*  INTEGER Sampling interval
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(2) ::ncount, nstart, nstride
INTEGER, DIMENSION(SIZE(logdat)) :: intdat

procname(pglev+1) = 'cf90_put_var_log_1d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(intdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
intdat = MERGE (1,0,logdat)
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,intdat,start=nstart(1:1),&
                        & count=ncount(1:1),stride=nstride(1:1))
ELSE
   errstat = NF90_put_var(ncid,varid,intdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_log_1d

!========================================================================

SUBROUTINE cf90_put_var_log_2d(ncid,varid,logdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_log_2d* Write a 2-D array of logical data to an open netCDF
!                       file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.11
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(2) :: count, start, stride
LOGICAL, INTENT(IN), DIMENSION(:,:) :: logdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*logdat*  LOGICAL Output logical 2-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for output buffer
!*count*   INTEGER Number of output data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(3) ::ncount, nstart, nstride
INTEGER, DIMENSION(LBOUND(logdat,DIM=1):UBOUND(logdat,DIM=1),&
                 & LBOUND(logdat,DIM=2):UBOUND(logdat,DIM=2)) :: intdat

procname(pglev+1) = 'cf90_put_var_log_2d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(logdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
intdat = MERGE (1,0,logdat)
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,intdat,start=nstart(1:2),&
                        & count=ncount(1:2),stride=nstride(1:2))
ELSE
   errstat = NF90_put_var(ncid,varid,intdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_log_2d

!========================================================================

SUBROUTINE cf90_put_var_real_0d(ncid,varid,realdat,timerec)
!************************************************************************
!
! *cf90_put_var_real_0d* Write a single precision real scalar variable to an
!                        open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!
USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
REAL (KIND=kndreal), INTENT(IN) :: realdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Output real scalar
!*timerec* INTEGER Time record number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(1) :: nstart


procname(pglev+1) = 'cf90_put_var_real_0d'
CALL log_timer_in(npcc)

nstart = timerec

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,realdat)
ELSE
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_real_0d

!========================================================================

SUBROUTINE cf90_put_var_real_1d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_real_1d* Write a vector of single precision real data to an
!                        open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(1) :: count, start, stride
REAL (KIND=kndreal), INTENT(IN), DIMENSION(:) :: realdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Output real vector array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start index for output buffer
!*count*   INTEGER Number of output data
!*stride*  INTEGER Sampling interval
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(2) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_put_var_real_1d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart(1:1),&
                        & count=ncount(1:1),stride=nstride(1:1))
ELSE
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_real_1d

!========================================================================

SUBROUTINE cf90_put_var_real_2d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_real_2d* Write a 2-D array of single precision real data to an
!                        open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!
USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(2) :: count, start, stride
REAL (KIND=kndreal), INTENT(IN), DIMENSION(:,:) :: realdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Output real 2-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for output buffer
!*count*   INTEGER Number of output data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(3) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_put_var_real_2d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart(1:2),&
                        & count=ncount(1:2),stride=nstride(1:2))
ELSE
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_real_2d

!========================================================================

SUBROUTINE cf90_put_var_real_3d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_real_3d* Write a 3-D array of single precision real data to an
!                        open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!
USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(3) :: count, start, stride
REAL (KIND=kndreal), INTENT(IN), DIMENSION(:,:,:) :: realdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Output real 3-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for output buffer
!*count*   INTEGER Number of output data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(4) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_put_var_real_3d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart(1:3),&
                        & count=ncount(1:3),stride=nstride(1:3))
ELSE
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_real_3d

!========================================================================

SUBROUTINE cf90_put_var_real_4d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_real_4d* Write a 4-D array of single precision real data to an
!                        open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches
!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(4) :: count, start, stride
REAL (KIND=kndreal), INTENT(IN), DIMENSION(:,:,:,:) :: realdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Vairalbe ID
!*realdat* REAL    Output real 4-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for output buffer
!*count*   INTEGER Number of output data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(5) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_put_var_real_4d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart(1:4),&
                        & count=ncount(1:4),stride=nstride(1:4))
ELSE
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_real_4d

!========================================================================

SUBROUTINE cf90_put_var_double_0d(ncid,varid,realdat,timerec)
!************************************************************************
!
! *cf90_put_var_double_0d* Write a double precision real scalar variable to an
!                          open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
REAL (KIND=kndrlong), INTENT(IN) :: realdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Output real scalar
!*timerec* INTEGER Time record number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(1) :: nstart


procname(pglev+1) = 'cf90_put_var_double_0d'
CALL log_timer_in(npcc)

nstart = timerec

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,realdat)
ELSE
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_double_0d

!========================================================================

SUBROUTINE cf90_put_var_double_1d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_double_1d* Write a vector of double precision real data to an
!                          open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(1) :: count, start, stride
REAL (KIND=kndrlong), INTENT(IN), DIMENSION(:) :: realdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Output real vector array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start index for output buffer
!*count*   INTEGER Number of output data
!*stride*  INTEGER Sampling interval
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(2) :: ncount, nstart, nstride


procname(pglev+1) = 'cf90_put_var_double_1d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart(1:1),&
                        & count=ncount(1:1),stride=nstride(1:1))
ELSE
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_double_1d

!========================================================================

SUBROUTINE cf90_put_var_double_2d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_double_2d* Write a 2-D array of double precision real data to
!                          an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(2) :: count, start, stride
REAL (KIND=kndrlong), INTENT(IN), DIMENSION(:,:) :: realdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Output real 2-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for output buffer
!*count*   INTEGER Number of output data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(3) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_put_var_double_2d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart(1:2),&
                        & count=ncount(1:2),stride=nstride(1:2))
ELSE
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_double_2d

!========================================================================

SUBROUTINE cf90_put_var_double_3d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_double_3d* Write a 3-D array of double precision real data to
!                          an open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(3) :: count, start, stride
REAL (KIND=kndrlong), INTENT(IN), DIMENSION(:,:,:) :: realdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*realdat* REAL    Output real 3-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for output buffer
!*count*   INTEGER Number of output data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(4) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_put_var_double_3d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart(1:3),&
                        & count=ncount(1:3),stride=nstride(1:3))
ELSE
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_double_3d

!========================================================================

SUBROUTINE cf90_put_var_double_4d(ncid,varid,realdat,timerec,start,count,stride)
!************************************************************************
!
! *cf90_put_var_real_4d* Write a 4-D array of double precision real data to an
!                        open netCDF file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.10.1
!
! Description -
!
! netCDF calls - cf90_sync, NF90_put_var
!
!************************************************************************
!

USE switches

!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, timerec, varid
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(4) :: count, start, stride
REAL (KIND=kndrlong), INTENT(IN), DIMENSION(:,:,:,:) :: realdat 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Vairalbe ID
!*realdat* REAL    Output real 4-D array
!*timerec* INTEGER Time record number
!*start*   INTEGER Start indices for output buffer
!*count*   INTEGER Number of output data along each dimension
!*stride*  INTEGER Sampling interval along each dimension
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(5) ::ncount, nstart, nstride


procname(pglev+1) = 'cf90_put_var_double_4d'
CALL log_timer_in(npcc)

!---optional arguments for netcdf call
IF (PRESENT(start)) THEN
   nstart = (/start,timerec/)
ELSE
   nstart = (/1,1,1,1,timerec/)
ENDIF
IF (PRESENT(count)) THEN
   ncount = (/count,1/)
ELSE
   ncount = (/SHAPE(realdat),1/)
ENDIF
IF (PRESENT(stride)) THEN
   nstride = (/stride,1/)
ELSE
   nstride = 1
ENDIF

#ifdef CDF
IF (timerec.EQ.0) THEN
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart(1:4),&
                        & count=ncount(1:4),stride=nstride(1:4))
ELSE
   errstat = NF90_put_var(ncid,varid,realdat,start=nstart,count=ncount,&
                        & stride=nstride)
ENDIF
IF (iopt_CDF_sync.EQ.1) CALL cf90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_put_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_put_var_double_4d

!========================================================================

SUBROUTINE cf90_redef(ncid)
!************************************************************************
!
! *cf90_redef* Puts an open netCDF file in define mode
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! netCDF calls - NF90_redef
!
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_redef'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_redef(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_redef')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_redef

!========================================================================

SUBROUTINE cf90_rename_att(ncid,varid,curname,newname)
!************************************************************************
!
! *cf90_rename_att* Rename an existing attribute
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! netCDF calls - NF90_rename_att
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: curname, newname
INTEGER, INTENT(IN) :: ncid, varid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER ID of the associated variable
!*curname* CHAR    Current attribute's name
!*newname* CHAR    New attribute's name
!
!------------------------------------------------------------------------------
!

procname(pglev+1) = 'cf90_rename_att'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_rename_att(ncid,varid,curname,newname)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_rename_att')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_rename_att

!========================================================================

SUBROUTINE cf90_rename_dim(ncid,dimid,name)
!************************************************************************
!
! *cf90_rename_dim* Rename an existing dimension
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! netCDF calls - NF90_rename_dim
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: name
INTEGER, INTENT(IN) :: dimid, ncid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*dimid*   INTEGER Dimension ID
!*name*    CHAR    New dimension name
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_rename_dim'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_rename_dim(ncid,dimid,name)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_rename_dim')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_rename_dim

!========================================================================

SUBROUTINE cf90_rename_var(ncid,varid,newname)
!************************************************************************
!
! *cf90_rename_var* Change the name of a netCDF variable
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! netCDF calls - NF90_rename_var
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: newname
INTEGER, INTENT(IN) :: ncid, varid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!*varid*   INTEGER Variable ID
!*newname* CHAR    New variable name
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_rename_var'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_rename_var(ncid,varid,newname)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_rename_var')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_rename_var

!========================================================================

SUBROUTINE cf90_set_fill(ncid,fillmode,old_mode)
!************************************************************************
!
! *cf90_set_fill* Set fill mode for a netCDF dataset open for writing
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! netCDF calls - NF90_set_fill
!
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid, fillmode
INTEGER, INTENT(OUT) :: old_mode

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*ncid*     INTEGER NetCDF ID
!*fillmode* INTEGER New fill mode for the netCDF file
!*old_mode* INTEGER Old fill mode of the netCDF file
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_set_fill'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_set_fill(ncid,fillmode,old_mode)
#else
old_mode = 0
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_set_fill')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_set_fill

!========================================================================

SUBROUTINE cf90_sync(ncid)
!************************************************************************
!
! *cf90_sync* Synchronize disk copy of a netCDF dataset
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)cf90_routines.F90  V2.0
!
! Description -
!
! netCDF calls - NF90_sync
!
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: ncid

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ncid*    INTEGER NetCDF ID
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'cf90_sync'
CALL log_timer_in(npcc)

#ifdef CDF
errstat = NF90_sync(ncid)
#endif /*CDF*/

IF (errstat.NE.noerr_NF90) CALL cf90_error('NF90_sync')

CALL log_timer_out(npcc,itm_CDF)


RETURN

END SUBROUTINE cf90_sync


END MODULE cf90_routines
