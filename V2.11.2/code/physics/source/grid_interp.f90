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

MODULE grid_interp
!************************************************************************
!
! *grid_interp* Interpolations between model and data grids 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! $Date: 2018-06-14 15:25:04 +0200 (Thu, 14 Jun 2018) $
!
! $Revision: 1155 $
!
! Routines - data_to_data_vcoords_2d, data_to_model_hcoords,
!            data_to_model_vcoords_1d, data_to_model_vcoords_2d,
!            data_to_model_vcoords_3d, hrel_coords_curv_abs,
!            hrel_coords_curv_rel, hrel_coords_rect_nonunif_abs,
!            hrel_coords_rect_nonunif_rel, hrel_coords_rect_unif_abs,
!            hrel_coords_rect_unif_rel, intpol_data_to_model_2d,
!            intpol1d_model_to_dep, intpol2d_model_to_data_2d,
!            intpol3d_data_to_data_2d, intpol3d_model_to_data_1d,
!            intpol3d_model_to_data_2d, intpol3d_model_to_data_3d,
!            model_to_data_hcoords
!
!************************************************************************
!
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
! Interfaces
!

INTERFACE hrel_coords_curv_abs
   MODULE PROCEDURE hrel_coords_curv_abs_s, hrel_coords_curv_abs_d
END INTERFACE

INTERFACE hrel_coords_curv_rel
   MODULE PROCEDURE hrel_coords_curv_rel_s, hrel_coords_curv_rel_d
END INTERFACE
 
INTERFACE hrel_coords_rect_nonunif_abs
   MODULE PROCEDURE hrel_coords_rect_nonunif_abs_s, &
                  & hrel_coords_rect_nonunif_abs_d
END INTERFACE

INTERFACE hrel_coords_rect_nonunif_rel
   MODULE PROCEDURE hrel_coords_rect_nonunif_rel_s, &
                  & hrel_coords_rect_nonunif_rel_d
END INTERFACE
 
INTERFACE hrel_coords_rect_unif_abs
   MODULE PROCEDURE hrel_coords_rect_unif_abs_s, hrel_coords_rect_unif_abs_d
END INTERFACE

INTERFACE hrel_coords_rect_unif_rel
   MODULE PROCEDURE hrel_coords_rect_unif_rel_s, hrel_coords_rect_unif_rel_d
END INTERFACE


CONTAINS

  
!========================================================================

SUBROUTINE data_to_data_vcoords_2d(zin,zout,vcoords,nhdat,nzdatin,nzdatout,&
                                 & nzeff)
!************************************************************************
!
! *data_to_data_vcoords_2d* Relative vertical coordinates of an input data
!                           grid with respect to an output data grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.9
!
! Description - the two data grids share the same horizontal locations
!             - vertical dimensions may be different
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
INTEGER, INTENT(IN) :: nhdat, nzdatin, nzdatout, nzeff
REAL, INTENT(IN), DIMENSION(nhdat,nzdatin) :: zin
REAL, INTENT(IN), DIMENSION(nhdat,nzdatout) :: zout
TYPE (VRelativeCoords), INTENT(OUT), DIMENSION(nhdat,nzdatout) :: vcoords

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*zin*        REAL    Z-coordinates (surface distance) of input data points 
!*zout*       REAL    Z-coordinates (surface distance) of output data points 
!*vcoords*    DERIVED Vertical coordinate array structure for interpolation
!*nhdat*      INTEGER Number of horizontal data locations
!*nzdatin*    INTEGER Vertical dimension of input data
!*nzdatout*   INTEGER Vertical dimension of output grid
!*nzeff*      INTEGER Effective vertical dimension of output grid
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: kc, kdat, ldat
REAL :: dsig, z


IF ((SIZE(zin).EQ.0).OR.(nzeff.EQ.0)) RETURN

procname(pglev+1) = 'data_to_data_vcoords_2d'
CALL log_timer_in()

ldat_110: DO ldat=1,nhdat
kdat_110: DO kdat=1,nzeff
   z = zout(ldat,kdat)
   IF (ABS(z-float_fill).GT.float_fill_eps) THEN
      IF (z.LE.zin(ldat,1)) THEN
         vcoords(ldat,kdat)%kcoord = 1
         vcoords(ldat,kdat)%weights = (/1.0,0.0/)
      ELSEIF (z.GE.zin(ldat,nzdatin)) THEN
         vcoords(ldat,kdat)%kcoord = nzdatin-1
         vcoords(ldat,kdat)%weights = (/0.0,1.0/)
      ELSE
         kc = 1
         DO WHILE (z.GE.zin(ldat,kc+1))
            kc = kc + 1
         END DO
         vcoords(ldat,kdat)%kcoord = kc
         dsig = (z-zin(ldat,kc))/(zin(ldat,kc+1)-zin(ldat,kc))
         vcoords(ldat,kdat)%weights = (/1-dsig,dsig/)
      ENDIF
   ELSE
      vcoords(ldat,kdat)%kcoord = int_fill
      vcoords(ldat,kdat)%weights = 0.0
   ENDIF
ENDDO kdat_110
ENDDO ldat_110

CALL log_timer_out()


RETURN

END SUBROUTINE data_to_data_vcoords_2d

!========================================================================

SUBROUTINE data_to_model_hcoords(surfcoords,cnode,nxdat,nydat,xdat,ydat,&
                               & maskdat,extrapol,land)
!************************************************************************
!
! *data_to_model_hcoords* Relative coordinates and weight factors used to
!                         interpolate external (gridded or non-gridded)
!                         data on the model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description - coordinates are calculated with respect to the model grid node
!               selected by cnode
!             - model and data grid are rotated in case of a rotated model grid
!             - either error checking or extrapolation is performed for data
!               points outside the model grid or surrounded by land cells
!               on the model grid
!
! Module calls - distance_minloc, error_abort, error_alloc, global_mask,
!                grid_rotation_transf, hrel_coords_curv_abs,
!                hrel_coords_rect_nounif_abs, hrel_coords_rect_unif_abs,
!                UVarr_at_C, UVarr_at_U, UVarr_at_V
!
!************************************************************************
!
USE datatypes
USE grid
USE gridpars
USE modids
USE paralpars
USE switches
USE syspars
USE array_interp, ONLY: UVarr_at_C, UVarr_at_U, UVarr_at_V
USE error_routines, ONLY: error_abort, error_alloc
USE grid_routines, ONLY: distance_minloc, global_mask, &
                       & grid_rotation_transf

!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: extrapol, land
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: nxdat, nydat
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(nxdat,nydat) :: maskdat
REAL, INTENT(IN), DIMENSION(nxdat,nydat) :: xdat, ydat
TYPE (HRelativeCoords), INTENT(INOUT), DIMENSION(nxdat,nydat) :: surfcoords

!
! Name         Type   Purpose
!-------------------------------------------------------------------------------
!*surfcoords* DERIVED Relative coordinates and weight factors for interpolation
!                     on the model grid
!*cnode*      CHAR    Nodal type of the model grid ('C',' U', 'V', X', 'Y','UV')
!*nxdat*      INTEGER X-dimension of the external grid 
!*nydat*      INTEGER Y-dimension of the external grid 
!*xdat*       REAL    X-coordinates of the locations on the external grid
!                                                       [m or degrees longitude]
!*ydat*       REAL    Y-coordinates of the locations on the external grid
!                                                        [m or degrees latitude]
!*maskdat*    LOGICAL Mask for land points on the external grid
!*extrapol*   LOGICAL Allows extrapolation to the nearest inside grid point
!                     for data points outside the domain if PRESENT AND .TRUE.
!                     Default is .FALSE.
!*land*       LOGICAL Allows extrapolation for data points inside a land
!                     cell on the model grid  if PRESENT and .TRUE.
!                     Default is .FALSE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: extrap, flagdat, landx
CHARACTER (LEN=12) :: cidat, cjdat
CHARACTER (LEN=15) :: cxdat, cydat
INTEGER :: icoord, idat, ierr, istat, jcoord, jdat, nhtype
REAL :: delxdat, delydat, xpos, ypos, wsum
INTEGER, DIMENSION(2) :: lbounds, locmin, ubounds 
REAL, DIMENSION(2,2) :: weights
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskglb
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: nstatus
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: xcoord, xdatrot, ycoord, ydatrot  


procname(pglev+1) = 'data_to_model_hcoords: '//TRIM(cnode)
CALL log_timer_in()

!
!1. Initialise
!-------------
!
!1.1 Optional arguments
!----------------------
!

flagdat = PRESENT(maskdat)

IF (PRESENT(extrapol)) THEN
   extrap = extrapol
ELSE
   extrap = .FALSE.
ENDIF

IF (PRESENT(land)) THEN
   landx = land
ELSE
   landx = .FALSE.
ENDIF

!
!1.2 Initialise parameters
!-------------------------
!

nhtype = iopt_grid_htype
IF (nhtype.EQ.1) THEN
   IF (iopt_grid_sph.EQ.0) THEN
      delxdat = surfacegrids(igrd_model,1)%delxdat
      delydat = surfacegrids(igrd_model,1)%delydat
   ELSE
      delxdat = degtorad*surfacegrids(igrd_model,1)%delxdat
      delydat = degtorad*surfacegrids(igrd_model,1)%delydat
   ENDIF
ENDIF

lbounds = 1
SELECT CASE(TRIM(cnode))
   CASE ('C'); ubounds = (/nc-1,nr-1/)
   CASE ('U'); ubounds = (/nc,nr-1/)
   CASE ('V'); ubounds = (/nc-1,nr/)
   CASE ('UV'); ubounds = (/nc,nr/)
END SELECT
   
!
!1.3 Allocate arrays
!-------------------
!

ALLOCATE (maskglb(nc,nr),STAT=errstat)
CALL error_alloc('maskglb',2,(/nc,nr/),kndlog)
ALLOCATE (xcoord(nc,nr),STAT=errstat)
CALL error_alloc('xcoord',2,(/nc,nr/),kndrtype)
ALLOCATE (ycoord(nc,nr),STAT=errstat)
CALL error_alloc('ycoord',2,(/nc,nr/),kndrtype)
ALLOCATE (xdatrot(nxdat,nydat),STAT=errstat)
CALL error_alloc('xdatrot',2,(/nxdat,nydat/),kndrtype)
ALLOCATE (ydatrot(nxdat,nydat),STAT=errstat)
CALL error_alloc('ydatrot',2,(/nxdat,nydat/),kndrtype)
ALLOCATE (nstatus(nxdat,nydat),STAT=errstat)
CALL error_alloc('nstatus',2,(/nxdat,nydat/),kndint)
nstatus = 0

!
!1.4 Global mask
!---------------
!

IF (TRIM(cnode).EQ.'C') THEN
   maskglb = seapointglb(1:nc,1:nr)
ELSE
   CALL global_mask(maskglb,cnode,2)
ENDIF

!
!1.5 Coordinate arrays
!---------------------
!
!1.5.1 At specified node
!-----------------------
!

SELECT CASE(TRIM(cnode))
CASE ('C')
   CALL UVarr_at_C(gxcoordglb(1:nc+1,1:nr+1),xcoord,0,0,(/1,1,1/),&
                & (/nc+1,nr+1,1/),1,iarr_gxcoordglb,.TRUE.)
   CALL UVarr_at_C(gycoordglb(1:nc+1,1:nr+1),ycoord,0,0,(/1,1,1/),&
                & (/nc+1,nr+1,1/),1,iarr_gycoordglb,.TRUE.)
CASE ('U','Y')
   CALL UVarr_at_U(gxcoordglb(1:nc,1:nr+1),xcoord,0,0,(/1,1,nz/),&
                & (/nc,nr+1,nz/),1,iarr_gxcoordglb,.TRUE.)
   CALL UVarr_at_U(gycoordglb(1:nc,1:nr+1),ycoord,0,0,(/1,1,nz/),&
                & (/nc,nr+1,nz/),1,iarr_gycoordglb,.TRUE.)
CASE ('V','X')
   CALL UVarr_at_V(gxcoordglb(1:nc+1,1:nr),xcoord,0,0,(/1,1,nz/),&
                & (/nc+1,nr,nz/),1,iarr_gxcoordglb,.TRUE.)
   CALL UVarr_at_V(gycoordglb(1:nc+1,1:nr),ycoord,0,0,(/1,1,nz/),&
                & (/nc+1,nr,nz/),1,iarr_gycoordglb,.TRUE.)
CASE ('UV')
   xcoord = gxcoordglb(1:nc,1:nr)
   ycoord = gycoordglb(1:nc,1:nr)
END SELECT
IF (iopt_grid_sph.EQ.1) THEN
   xcoord = degtorad*xcoord; ycoord = degtorad*ycoord
ENDIF

!
!1.5.2 Rotate (if needed)
!------------------------
!

xdatrot = xdat; ydatrot = ydat
IF (surfacegrids(igrd_model,1)%rotated) THEN
   CALL grid_rotation_transf(xcoord,ycoord,1,1,igrd_model)
   CALL grid_rotation_transf(xdatrot,ydatrot,1,1,igrd_model)
ENDIF

!
!2. Relative coordinates and weights
!-----------------------------------
!

jdat_210: DO jdat=1,nydat
idat_210: DO idat=1,nxdat

!  ---land mask   
   IF (flagdat) THEN
      IF (.NOT.maskdat(idat,jdat)) THEN
         nstatus(idat,jdat) = -1
         CYCLE idat_210
      ENDIF
   ENDIF

!  ---relative coordinates and weights
   IF (iopt_grid_sph.EQ.0) THEN
      xpos = xdatrot(idat,jdat); ypos = ydat(idat,jdat)
   ELSE
      xpos = degtorad*xdatrot(idat,jdat); ypos = degtorad*ydatrot(idat,jdat)
   ENDIF
   SELECT CASE (nhtype)
      CASE (1)
         CALL hrel_coords_rect_unif_abs(xpos,ypos,xcoord(1,1),ycoord(1,1),&
              & delxdat,delydat,lbounds,ubounds,istat,icoord,jcoord,&
              & weights=weights,extrapol=extrap)
      CASE (2)
         CALL hrel_coords_rect_nonunif_abs(xpos,ypos,xcoord(:,1),ycoord(1,:),&
              & nc,nr,lbounds,ubounds,istat,icoord,jcoord,weights=weights,&
              & extrapol=extrap)
      CASE (3)
         CALL hrel_coords_curv_abs(xpos,ypos,xcoord,ycoord,nc,nr,lbounds,&
              & ubounds,istat,1,icoord,jcoord,weights=weights,extrapol=extrap)
   END SELECT

!  ---status of location (outside domain or land), store data
   IF (istat.EQ.1.AND.(.NOT.extrap)) THEN
      nstatus(idat,jdat) = 1
   ELSE
      surfcoords(idat,jdat)%icoord = icoord
      surfcoords(idat,jdat)%jcoord = jcoord
      weights = MERGE(weights,0.0,maskglb(icoord:icoord+1,jcoord:jcoord+1))
      wsum = SUM(weights)
      IF (wsum.EQ.0.0) THEN
         nstatus(idat,jdat) = 2
      ELSE
         surfcoords(idat,jdat)%weights = weights/wsum
      ENDIF
   ENDIF
ENDDO idat_210
ENDDO jdat_210

!
!3. Perform extrapolation in case the location is located on land
!----------------------------------------------------------------
!

jdat_310: DO jdat=1,nydat
idat_310: DO idat=1,nxdat
   IF (nstatus(idat,jdat).EQ.2.AND.landx) THEN
      IF (iopt_grid_sph.EQ.0) THEN
         xpos = xdatrot(idat,jdat); ypos = ydat(idat,jdat)
      ELSE
         xpos = degtorad*xdatrot(idat,jdat); ypos = degtorad*ydatrot(idat,jdat)
      ENDIF
      CALL distance_minloc(xpos,ypos,xcoord,ycoord,nc,nr,1,3,locmin=locmin,&
                         & mask=maskglb)
      surfcoords(idat,jdat)%icoord = locmin(1)
      surfcoords(idat,jdat)%jcoord = locmin(2)
      surfcoords(idat,jdat)%weights(1:2,1) = (/1.0,0.0/)
      surfcoords(idat,jdat)%weights(1:2,2) = 0.0
      nstatus(idat,jdat) = 0
   ENDIF
ENDDO idat_310
ENDDO jdat_310

!
!4. Check locations
!------------------
!
   
nerrs = COUNT(nstatus.GT.0)

IF (master.AND.nerrs.GT.0.AND.errchk) THEN

   ierr = 0
   idat_410: DO idat=1,nxdat
   jdat_410: DO jdat=1,nydat
      IF (nstatus(idat,jdat).GT.0.AND.ierr.LE.maxerrors) THEN
         WRITE (cidat,'(I12)') idat; cidat = ADJUSTL(cidat)
         WRITE (cjdat,'(I12)') jdat; cjdat = ADJUSTL(cjdat)
         WRITE (cxdat,'(G15.7)') xdat(idat,jdat); cxdat = ADJUSTL(cxdat)
         WRITE (cydat,'(G15.7)') ydat(idat,jdat); cydat = ADJUSTL(cydat)
         IF (nstatus(idat,jdat).EQ.1) THEN
            IF (nydat.GT.1) THEN
               WRITE (ioerr,'(A)') 'Data point ('//TRIM(cidat)//','//&
                       & TRIM(cjdat)//') at ('//TRIM(cxdat)//','//TRIM(cydat)//&
                       & ') is outside the data grid'
            ELSE
               WRITE (ioerr,'(A)') 'Data point '//TRIM(cidat)//' at ('&
                       & //TRIM(cxdat)//','//TRIM(cydat)//&
                       & ') is outside the data grid'
            ENDIF
         ELSEIF (nstatus(idat,jdat).EQ.2) THEN
            IF (nydat.GT.1) THEN
               WRITE (ioerr,'(A)') 'Data point ('//TRIM(cidat)//','//&
                       & TRIM(cjdat)//') at ('//TRIM(cxdat)//','//TRIM(cydat)//&
                       & ') is located on land'
            ELSE
               WRITE (ioerr,'(A)') 'Data point '//TRIM(cidat)//' at ('&
                       & //TRIM(cxdat)//','//TRIM(cydat)//') is located on land'
            ENDIF
         ENDIF
         ierr = ierr + 1
      ENDIF
   ENDDO jdat_410
   ENDDO idat_410

ENDIF

CALL error_abort('data_to_model_hcoords',ierrno_inival)

!
!5. Deallocate arrays
!--------------------
!

DEALLOCATE (maskglb,nstatus,xcoord,xdatrot,ycoord,ydatrot)

CALL log_timer_out()


RETURN

END SUBROUTINE data_to_model_hcoords

!==============================================================================

SUBROUTINE data_to_model_vcoords_1d(zcoord,vcoords,i,j,nzdat,cnode,checkdata)
!******************************************************************************
!
! *data_to_model_vcoords_1d* Relative vertical coordinates of a 1-D vertical
!                            data grid with respect to the local model grid
!
! Version - @(COHERENS)grid_interp.f90  V2.9
!
! Author - Patrick Luyten
!
! Description - coordinates are calculated at grid point (i,j)
!
! Module calls - error_abort
!
!******************************************************************************
!
USE datatypes
USE depths
USE grid
USE gridpars
USE syspars
USE error_routines, ONLY: error_abort

!
!*Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
LOGICAL, INTENT(IN), OPTIONAL :: checkdata
INTEGER, INTENT(IN) :: i, j, nzdat
REAL, INTENT(IN), DIMENSION(nzdat) :: zcoord
TYPE (VRelativeCoords), INTENT(OUT), DIMENSION(nzdat) :: vcoords

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*zcoord*      REAL    Z-coordinates (surface distance) of data points
!*vcoords*     DERIVED Vertical coordinate array structure for interpolation
!*i*           INTEGER X-index coordinate of model/data point
!*j*           INTEGER Y-index coordinate of model/data point
!*nzdat*       INTEGER Vertical dimension of data grid
!*cnode*       CHAR    Nodal type of model grid ('C','U','V',X','Y')
!*checkdata*   LOGICAL Performs error checking if .TRUE. (default value)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cind
CHARACTER (LEN=15) :: cpos
LOGICAL :: check
INTEGER :: k, kc
REAL :: dsig, sigk
REAL, DIMENSION(nz) :: gsig


IF (SIZE(zcoord).EQ.0) RETURN

procname(pglev+1) = 'data_to_model_vcoords_1d: '//TRIM(cnode)
CALL log_timer_in()

!
!1. Optional argument
!--------------------
!

IF (PRESENT(checkdata)) THEN
   check = checkdata
ELSE
   check = .TRUE.
ENDIF

!
!2. C-nodes
!----------
!

SELECT CASE (TRIM(cnode))
CASE('C')

   IF (maskatc_int(i,j)) THEN
      gsig = gscoordatc(i,j,:)
      k_210: DO k=1,nzdat
         IF (ABS(zcoord(k)-float_fill).GT.float_fill_eps) THEN
            sigk = 1.0-zcoord(k)/depmeanatc(i,j)
            IF (sigk.LE.gsig(1)) THEN
               vcoords(k)%kcoord = 1
               vcoords(k)%weights = (/1.0,0.0/)
            ELSEIF (sigk.GE.gsig(nz)) THEN
               vcoords(k)%kcoord = nz-1
               vcoords(k)%weights = (/0.0,1.0/)
            ELSE
               kc = 1
               DO WHILE (sigk.GE.gsig(kc+1))
                  kc = kc + 1
               ENDDO
               vcoords(k)%kcoord = kc
               dsig = (sigk-gsig(kc))/(gsig(kc+1)-gsig(kc))
               vcoords(k)%weights = (/1-dsig,dsig/)
            ENDIF
         ENDIF
      ENDDO k_210
   ENDIF

!
!3. U-nodes
!----------
!

CASE('U','Y')

   IF (node2du(i,j).GT.1) THEN
      gsig = gscoordatu(i,j,:)
      k_310: DO k=1,nzdat
         IF (ABS(zcoord(k)-float_fill).GT.float_fill_eps) THEN
            sigk = 1.0-zcoord(k)/depmeanatu(i,j)
            IF (sigk.LE.gsig(1)) THEN
               vcoords(k)%kcoord = 1
               vcoords(k)%weights = (/1.0,0.0/)
            ELSEIF (sigk.GE.gsig(nz)) THEN
               vcoords(k)%kcoord = nz-1
               vcoords(k)%weights = (/0.0,1.0/)
            ELSE
               kc = 1
               DO WHILE (sigk.GE.gsig(kc+1))
                  kc = kc + 1
               ENDDO
               vcoords(k)%kcoord = kc
               dsig = (sigk-gsig(kc))/(gsig(kc+1)-gsig(kc))
               vcoords(k)%weights = (/1-dsig,dsig/)
            ENDIF
         ENDIF
      ENDDO k_310
   ENDIF

!
!4. V-nodes
!----------
!

CASE('V','X')

   IF (node2dv(i,j).GT.1) THEN
      gsig = gscoordatv(i,j,:)
      k_410: DO k=1,nzdat
         IF (ABS(zcoord(k)-float_fill).GT.float_fill_eps) THEN
            sigk = 1.0-zcoord(k)/depmeanatv(i,j)
            IF (sigk.LE.gsig(1)) THEN
               vcoords(k)%kcoord = 1
               vcoords(k)%weights = (/1.0,0.0/)
            ELSEIF (sigk.GE.gsig(nz)) THEN
               vcoords(k)%kcoord = nz-1
               vcoords(k)%weights = (/0.0,1.0/)
            ELSE
               kc = 1
               DO WHILE (sigk.GE.gsig(kc+1))
                  kc = kc + 1
               ENDDO
               vcoords(k)%kcoord = kc
               dsig = (sigk-gsig(kc))/(gsig(kc+1)-gsig(kc))
               vcoords(k)%weights = (/1-dsig,dsig/)
            ENDIF
         ENDIF
      ENDDO k_410
   ENDIF

END SELECT

!
!5. Error checking
!-----------------
!

IF (check) THEN
   nerrs = MERGE(1,0,ANY(vcoords%kcoord.EQ.int_fill))
   IF (nerrs.GT.0.AND.errchk) THEN
      WRITE (ioerr,'(A)') 'Unable to obtain vertical relative coordinates'//&
                       & ' at data points:'
      k_510: DO k=1,nzdat
         IF (vcoords(k)%kcoord.EQ.int_fill) THEN
            WRITE (cind,'(I12)') k; cind = ADJUSTL(cind)
            WRITE (cpos,'(G15.7)') zcoord(k); cpos = ADJUSTL(cpos)
            WRITE (ioerr,'(A)') TRIM(cind)//', '//TRIM(cpos)
         ENDIF
      ENDDO k_510
   ENDIF
   CALL error_abort('data_to_model_vcoords_1d',ierrno_inival)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE data_to_model_vcoords_1d

!========================================================================

SUBROUTINE data_to_model_vcoords_2d(zcoord,vcoords,nzdat,cnode,checkdata)
!************************************************************************
!
! *data_to_model_vcoords_2d* Relative vertical coordinates of a 3-D data grid
!                            with respect to the local model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.9
!
! Description - model and data are defined along the same horizontal grid
!             - coordinates are calculated at all horizontal grid points
!             - dry cells are taken as invalid data points
!
! Module calls - error_abort
!
!************************************************************************
!
USE datatypes
USE depths
USE grid
USE gridpars
USE syspars
USE error_routines, ONLY: error_abort

!
!*Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
LOGICAL, INTENT(IN), OPTIONAL :: checkdata
INTEGER, INTENT(IN) :: nzdat
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nzdat) :: zcoord
TYPE (VRelativeCoords), INTENT(OUT), DIMENSION(ncloc,nrloc,nzdat) :: vcoords

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*zcoord*      REAL    Z-coordinates (surface distance) of data points 
!*vcoords*     DERIVED Vertical coordinate array structure for interpolation
!*nzdat*       INTEGER Vertical dimension of data array
!*cnode*       CHAR    Nodal type of model grid ('C','U','V','X','Y')
!*checkdata*   LOGICAL Performs error checking if .TRUE. (default value)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12), DIMENSION(3) :: cind
CHARACTER (LEN=15) :: cpos
LOGICAL :: check
INTEGER :: i, j, k, kc
REAL :: dsig, sigk
REAL, DIMENSION(nz) :: gsig


IF (nzdat.EQ.0) RETURN

procname(pglev+1) = 'data_to_model_vcoords_2d: '//TRIM(cnode)
CALL log_timer_in()

!
!1. Optional argument
!--------------------
!

IF (PRESENT(checkdata)) THEN
   check = checkdata
ELSE
   check = .TRUE.
ENDIF

!
!2. C-nodes
!----------
!

SELECT CASE(TRIM(cnode))

CASE ('C')

   j_210: DO j=1,nrloc
   i_210: DO i=1,ncloc
      IF (maskatc_int(i,j)) THEN
         gsig = gscoordatc(i,j,:)
         k_211: DO k=1,nzdat
            IF (ABS(zcoord(i,j,k)-float_fill).GT.float_fill_eps) THEN
               sigk = 1.0-zcoord(i,j,k)/depmeanatc(i,j)
               IF (sigk.LE.gsig(1)) THEN
                  vcoords(i,j,k)%kcoord = 1
                  vcoords(i,j,k)%weights = (/1.0,0.0/)
               ELSEIF (sigk.GE.gsig(nz)) THEN
                  vcoords(i,j,k)%kcoord = nz-1
                  vcoords(i,j,k)%weights = (/0.0,1.0/)
               ELSE
                  kc = 1
                  DO WHILE (sigk.GE.gsig(kc+1))
                     kc = kc + 1
                  ENDDO
                  vcoords(i,j,k)%kcoord = kc
                  dsig = (sigk-gsig(kc))/(gsig(kc+1)-gsig(kc))
                  vcoords(i,j,k)%weights = (/1.0-dsig,dsig/)
               ENDIF
            ENDIF
         ENDDO k_211
      ENDIF
   ENDDO i_210
   ENDDO j_210

!
!3. U-nodes
!----------
!

CASE('U','Y')

   j_310: DO j=1,nrloc
   i_310: DO i=1,ncloc
      IF (node2du(i,j).GT.1) THEN
         gsig = gscoordatu(i,j,:)
         k_311: DO k=1,nzdat
            IF (ABS(zcoord(i,j,k)-float_fill).GT.float_fill_eps) THEN
               sigk = 1.0-zcoord(i,j,k)/depmeanatu(i,j)
               IF (sigk.LE.gsig(1)) THEN
                  vcoords(i,j,k)%kcoord = 1
                  vcoords(i,j,k)%weights = (/1.0,0.0/)
               ELSEIF (sigk.GE.gsig(nz)) THEN
                  vcoords(i,j,k)%kcoord = nz-1
                  vcoords(i,j,k)%weights = (/0.0,1.0/)
               ELSE
                  kc = 1
                  DO WHILE (sigk.GE.gsig(kc+1))
                     kc = kc + 1
                  ENDDO
                  vcoords(i,j,k)%kcoord = kc
                  dsig = (sigk-gsig(kc))/(gsig(kc+1)-gsig(kc))
                  vcoords(i,j,k)%weights = (/1.0-dsig,dsig/)
               ENDIF
            ENDIF
         ENDDO k_311
      ENDIF
   ENDDO i_310
   ENDDO j_310

!
!4. V-nodes
!----------
!

CASE('V','X')

   j_410: DO j=1,nrloc
   i_410: DO i=1,ncloc
      IF (node2dv(i,j).GT.1) THEN
         gsig = gscoordatv(i,j,:)
         k_411: DO k=1,nzdat
            IF (ABS(zcoord(i,j,k)-float_fill).GT.float_fill_eps) THEN
               sigk = 1.0-zcoord(i,j,k)/depmeanatv(i,j)
               IF (sigk.LE.gsig(1)) THEN
                  vcoords(i,j,k)%kcoord = 1
                  vcoords(i,j,k)%weights = (/1.0,0.0/) 
               ELSEIF (sigk.GE.gsig(nz)) THEN
                  vcoords(i,j,k)%kcoord = nz-1
                  vcoords(i,j,k)%weights = (/0.0,1.0/)
               ELSE
                  kc = 1
                  DO WHILE (sigk.GE.gsig(kc+1))
                     kc = kc + 1
                  ENDDO
                  vcoords(i,j,k)%kcoord = kc
                  dsig = (sigk-gsig(kc))/(gsig(kc+1)-gsig(kc))
                  vcoords(i,j,k)%weights = (/1.0-dsig,dsig/)
               ENDIF
            ENDIF
         ENDDO k_411
      ENDIF
   ENDDO i_410
   ENDDO j_410

END SELECT

!
!5. Error checking
!-----------------
!

IF (check) THEN
   nerrs = MERGE(1,0,ANY(vcoords%kcoord.EQ.int_fill))
   IF (nerrs.GT.0.AND.errchk) THEN
      WRITE (ioerr,'(A)') 'Unable to obtain vertical relative coordinates'//&
                       & ' at data points:'
      i_510: DO i=1,ncloc 
      j_510: DO j=1,nrloc 
      k_510: DO k=1,nzdat
         IF (vcoords(i,j,k)%kcoord.EQ.int_fill) THEN
            WRITE (cind(1),'(I12)') i; cind(1) = ADJUSTL(cind(1))
            WRITE (cind(2),'(I12)') j; cind(2) = ADJUSTL(cind(2))
            WRITE (cind(3),'(I12)') k; cind(3) = ADJUSTL(cind(3))
            WRITE (cpos,'(G15.7)') zcoord(i,j,k); cpos = ADJUSTL(cpos)
            WRITE (ioerr,'(A)') TRIM(cind(1))//', '//TRIM(cind(2))//', '//&
                              & TRIM(cind(3))//', '//TRIM(cpos)
         ENDIF
      ENDDO k_510
      ENDDO j_510
      ENDDO i_510
   ENDIF
   CALL error_abort('data_to_model_vcoords_2d',ierrno_inival)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE data_to_model_vcoords_2d

!========================================================================

SUBROUTINE data_to_model_vcoords_3d(zcoord,hcoords,vcoords,n1dat,n2dat,n3dat,&
                                  & cnode,checkdata)
!************************************************************************
!
! *data_to_model_vcoords_3d* Relative vertical coordinates of a 3-D data grid
!                            with respect to the local model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.9
!
! Description - model and data may have different horizontal grids
!
! Module calls - error_abort
!  
!************************************************************************
!
USE datatypes
USE depths
USE grid
USE gridpars
USE syspars
USE error_routines, ONLY: error_abort

!
!*Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
LOGICAL, INTENT(IN), OPTIONAL :: checkdata
INTEGER, INTENT(IN) :: n1dat, n2dat, n3dat
REAL, INTENT(IN), DIMENSION(n1dat,n2dat,n3dat) :: zcoord
TYPE (HRelativeCoords), INTENT(IN), DIMENSION(n1dat*n2dat) :: hcoords
TYPE (VRelativeCoords), INTENT(OUT), DIMENSION(2,2,n1dat*n2dat,n3dat) :: &
                                             & vcoords

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*zcoord*      REAL    Z-coordinates (surface distance) of data points 
!*hcoords*     DERIVED Horizontal coordinate array structure for interpolation
!*vcoords*     DERIVED Vertical coordinate array structure for interpolation
!*n1dat*       INTEGER First (X-) dimension of vertical coordinate array
!*n2dat*       INTEGER Second (Y-) dimension of vertical coordinate array
!*n3dat*       INTEGER Third (vertical) dimension of vertical coordinate array
!*cnode*       CHAR    Nodal type of model grid ('C','U','V','X','Y')
!*checkdata*   LOGICAL Performs error checking if .TRUE. (default value)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12), DIMENSION(4) :: cind
CHARACTER (LEN=15) :: cpos
LOGICAL :: check
INTEGER :: i, idat, ii, isub, j, jdat, jj, jsub, kdat, kc, l
REAL :: dsig, sigk
REAL, DIMENSION(nz) :: gsig


IF (SIZE(zcoord).EQ.0) RETURN

procname(pglev+1) = 'data_to_model_vcoords_3d: '//TRIM(cnode)
CALL log_timer_in()

l = 0

!
!1. Optional argument
!--------------------
!

IF (PRESENT(checkdata)) THEN
   check = checkdata
ELSE
   check = .TRUE.
ENDIF

!
!2. C-nodes
!----------
!

SELECT CASE (TRIM(cnode))

CASE ('C')

   j_210: DO jdat=1,n2dat
   i_210: DO idat=1,n1dat
      l = l + 1
      i = hcoords(l)%icoord; j = hcoords(l)%jcoord
      ii_211: DO ii=i,i+1
      jj_211: DO jj=j,j+1
         isub = ii-i+1; jsub = jj-j+1
         IF (nodeatc(ii,jj).EQ.0) THEN
            vcoords(isub,jsub,l,:)%kcoord = 0
            vcoords(isub,jsub,l,:)%weights(1) = 0.0
            vcoords(isub,jsub,l,:)%weights(2) = 0.0
            CYCLE jj_211
         ENDIF
         gsig = gscoordatc(ii,jj,:)
         kdat_2111: DO kdat=1,n3dat
            IF (ABS(zcoord(idat,jdat,kdat)-float_fill).LE.float_fill_eps) THEN
               CYCLE kdat_2111
            ENDIF
            sigk = 1.0-zcoord(idat,jdat,kdat)/depmeanatc(ii,jj)
            IF (sigk.LE.gsig(1)) THEN
               vcoords(isub,jsub,l,kdat)%kcoord = 1
               vcoords(isub,jsub,l,kdat)%weights = (/1.0,0.0/)
            ELSEIF (sigk.GE.gsig(nz)) THEN
               vcoords(isub,jsub,l,kdat)%kcoord = nz-1
               vcoords(isub,jsub,l,kdat)%weights = (/0.0,1.0/)
            ELSE
               kc = 1
               DO WHILE (sigk.GE.gsig(kc+1))
                  kc = kc + 1
               ENDDO
               vcoords(isub,jsub,l,kdat)%kcoord = kc
               dsig = (sigk-gsig(kc))/(gsig(kc+1)-gsig(kc))
               vcoords(isub,jsub,l,kdat)%weights = (/1-dsig,dsig/)
            ENDIF
         ENDDO kdat_2111
      ENDDO jj_211
      ENDDO ii_211
   ENDDO i_210
   ENDDO j_210

!
!3. U-nodes
!----------
!

CASE('U','Y')

   j_310: DO jdat=1,n2dat
   i_310: DO idat=1,n1dat
      l = l + 1
      i = hcoords(l)%icoord; j = hcoords(l)%jcoord
      ii_311: DO ii=i,i+1
      jj_311: DO jj=j,j+1
         isub = ii-i+1; jsub = jj-j+1
         IF (node2du(ii,jj).LE.1) THEN
            vcoords(isub,jsub,l,:)%kcoord = 0
            vcoords(isub,jsub,l,:)%weights(1) = 0.0
            vcoords(isub,jsub,l,:)%weights(2) = 0.0
            CYCLE jj_311
         ENDIF
         gsig = gscoordatu(ii,jj,:)
         kdat_3111: DO kdat=1,n3dat
            IF (ABS(zcoord(idat,jdat,kdat)-float_fill).LE.float_fill_eps) THEN
               CYCLE kdat_3111
            ENDIF
            sigk = 1.0-zcoord(idat,jdat,kdat)/depmeanatu(ii,jj)
            IF (sigk.LE.gsig(1)) THEN
               vcoords(isub,jsub,l,kdat)%kcoord = 1
               vcoords(isub,jsub,l,kdat)%weights = (/1.0,0.0/)
            ELSEIF (sigk.GE.gsig(nz)) THEN
               vcoords(isub,jsub,l,kdat)%kcoord = nz-1
               vcoords(isub,jsub,l,kdat)%weights = (/0.0,1.0/)
            ELSE
               kc = 1
               DO WHILE (sigk.GE.gsig(kc+1))
                  kc = kc + 1
               ENDDO
               vcoords(isub,jsub,l,kdat)%kcoord = kc
               dsig = (sigk-gsig(kc))/(gsig(kc+1)-gsig(kc))
               vcoords(isub,jsub,l,kdat)%weights = (/1-dsig,dsig/)
            ENDIF
         ENDDO kdat_3111
      ENDDO jj_311
      ENDDO ii_311
   ENDDO i_310
   ENDDO j_310

!
!4. V-nodes
!----------
!

CASE('V','X')

   j_410: DO jdat=1,n2dat
   i_410: DO idat=1,n1dat
      l = l + 1
      i = hcoords(l)%icoord; j = hcoords(l)%jcoord
      ii_411: DO ii=i,i+1
      jj_411: DO jj=j,j+1
         isub = ii-i+1; jsub = jj-j+1
         IF (node2dv(ii,jj).LE.1) THEN
            vcoords(isub,jsub,l,:)%kcoord = 0
            vcoords(isub,jsub,l,:)%weights(1) = 0.0
            vcoords(isub,jsub,l,:)%weights(2) = 0.0
            CYCLE jj_411
         ENDIF
         gsig = gscoordatv(ii,jj,:)
         kdat_4111: DO kdat=1,n3dat
            IF (ABS(zcoord(idat,jdat,kdat)-float_fill).LE.float_fill_eps) THEN
               CYCLE kdat_4111
            ENDIF
            sigk = 1.0-zcoord(idat,jdat,kdat)/depmeanatv(ii,jj)
            IF (sigk.LE.gsig(1)) THEN
               vcoords(isub,jsub,l,kdat)%kcoord = 1
               vcoords(isub,jsub,l,kdat)%weights = (/1.0,0.0/)
            ELSEIF (sigk.GE.gsig(nz)) THEN
               vcoords(isub,jsub,l,kdat)%kcoord = nz-1
               vcoords(isub,jsub,l,kdat)%weights = (/0.0,1.0/)
            ELSE
               kc = 1
               DO WHILE (sigk.GE.gsig(kc+1))
                  kc = kc + 1
               ENDDO
               vcoords(isub,jsub,l,kdat)%kcoord = kc
               dsig = (sigk-gsig(kc))/(gsig(kc+1)-gsig(kc))
               vcoords(isub,jsub,l,kdat)%weights = (/1-dsig,dsig/)
            ENDIF
         ENDDO kdat_4111
      ENDDO jj_411
      ENDDO ii_411
   ENDDO i_410
   ENDDO j_410

END SELECT

!
!5. Error checking
!-----------------
!

IF (check) THEN

   nerrs = MERGE(1,0,ANY(vcoords%kcoord.EQ.int_fill))
   IF (nerrs.GT.0.AND.errchk) THEN

      WRITE (ioerr,'(A)') 'Unable to obtain vertical relative coordinates'//&
                       & ' at data points:'
      i_510: DO i=1,2
      j_510: DO j=1,2
      l_510: DO l=1,n1dat*n2dat
      k_510: DO kdat=1,n3dat
         IF (vcoords(i,j,l,kdat)%kcoord.EQ.int_fill) THEN
            WRITE (cind(1),'(I12)') i; cind(1) = ADJUSTL(cind(1))
            WRITE (cind(2),'(I12)') j; cind(2) = ADJUSTL(cind(2))
            WRITE (cind(3),'(I12)') l; cind(3) = ADJUSTL(cind(3))
            WRITE (cind(4),'(I12)') kdat; cind(4) = ADJUSTL(cind(4))
            WRITE (cpos,'(G15.7)') zcoord(i,j,kdat); cpos = ADJUSTL(cpos)
            WRITE (ioerr,'(A)') TRIM(cind(1))//', '//TRIM(cind(2))//', '//&
                              & TRIM(cind(3))//', '//TRIM(cind(4))//', '//&
                              & TRIM(cpos)
         ENDIF
      ENDDO k_510
      ENDDO l_510
      ENDDO j_510
      ENDDO i_510
   ENDIF
   CALL error_abort('data_to_model_vcoords_3d',ierrno_inival)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE data_to_model_vcoords_3d

!========================================================================

SUBROUTINE hrel_coords_curv_abs_s(xpos,ypos,xdat,ydat,nxdat,nydat,lbounds,&
                                & ubounds,istat,inunit,icoord,jcoord,xcoord,&
                                & ycoord,weights,extrapol)
!************************************************************************
!
! *hrel_coords_curv_abs_s* calculates the horizontal relative coordinates and
!                          weight factors for interpolating data from a
!                          curvilinear grid to a given location
!                    
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description - used in case of a curvilinear grid
!             - icoord, jcoord, xcoord, ycoord and weights are undefined
!               if istat=1 and extrapol is .FALSE.
!             - xcoord and ycoord must be both  or not present
!             - xpos and ypos given in single precision
!  
! Module calls - check_cell_location, distance_minloc
!
!************************************************************************
!
USE syspars  
USE grid_routines, ONLY: check_cell_location, distance_minloc    
  
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: extrapol
INTEGER, INTENT(IN) :: inunit, nxdat, nydat
INTEGER, INTENT(OUT) :: icoord, istat, jcoord
INTEGER, DIMENSION(2), INTENT(IN) :: lbounds, ubounds  
REAL (KIND=kndreal), INTENT(IN) :: xpos, ypos
REAL, INTENT(IN), DIMENSION(nxdat,nydat) :: xdat, ydat
REAL, INTENT(OUT), OPTIONAL :: xcoord, ycoord
REAL, INTENT(OUT), OPTIONAL, DIMENSION(2,2) :: weights

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*xpos*     REAL     X-coordinate of the input location [m or degrees longitude]
!*ypos*     REAL     Y-coordinate of the input location  [m or degrees latitude]
!*xdat*     REAL     X-coordinates of the input data grid
!                                                       [m or degrees longitude]
!*ydat*     REAL     Y-coordinate of the input data grid [m or degrees latitude]
!*nxdat*    INTEGER  X-dimension of the input data grid
!*nydat*    INTEGER  Y-dimension of the input data grid
!*lbounds*  INTEGER  Lower bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*ubounds*  INTEGER  Upper bounds of the bounding box delineating the inside
!*istat*    INTEGER  0/1 if location is inside/outside the input grid
!*inunit*   INTEGER  unit of input locations (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!*icoord*   INTEGER  Returned X-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*jcoord*   INTEGER  Returned Y-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*xcoord*   REAL     Returned normalised X-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*ycoord*   REAL     Returned normalised Y-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*weights*  REAL     Returned weight factors for interpolation
!*extrapol* LOGICAL  Allows extrapolation to the nearest inside grid point
!                    for data points outside the domain if PRESENT AND .TRUE.
!                    Default is .FALSE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: extrap, flag, frel
CHARACTER (LEN=15) :: cxpos, cypos
INTEGER, PARAMETER :: maxit = 100
INTEGER :: i, ii, imin, iter, j, jj, jmin
INTEGER, DIMENSION(2) :: locmin
REAL, PARAMETER :: epsconv = 1.0E-06
REAL :: denom, dx, dxh1, dxh2, dxh3, dy, dyh1, dyh2, dyh3, x, y
REAL, DIMENSION(2) :: a, b, c
REAL, DIMENSION(4) :: xvals, yvals


!
!1. Optional arguments
!---------------------
!

IF (PRESENT(extrapol)) THEN
   extrap = extrapol
ELSE
   extrap = .FALSE.
ENDIF

frel = PRESENT(xcoord).AND.PRESENT(ycoord)

!
!2. Search for containing grid cell
!----------------------------------
!
!---location of nearest grid point
CALL distance_minloc(xpos,ypos,xdat,ydat,nxdat,nydat,inunit,3,locmin=locmin)
imin = locmin(1); jmin = locmin(2) 

!---search loop
istat = 1
jj_210: DO jj=jmin-1,jmin
ii_210: DO ii=imin-1,imin
   IF (ii.GE.lbounds(1).AND.ii.LT.ubounds(1).AND.&
     & jj.GE.lbounds(2).AND.jj.LT.ubounds(2)) THEN
      xvals = (/xdat(ii,jj),xdat(ii+1,jj),xdat(ii+1,jj+1),xdat(ii,jj+1)/)
      yvals = (/ydat(ii,jj),ydat(ii+1,jj),ydat(ii+1,jj+1),ydat(ii,jj+1)/)
      CALL check_cell_location(xpos,ypos,xvals,yvals,flag)
      IF (flag) THEN
         istat = 0; i = ii; j = jj
         EXIT jj_210
      ENDIF
   ENDIF
ENDDO ii_210
ENDDO jj_210

!
!3. Extrapolate if needed
!------------------------
!

IF (istat.EQ.1.AND.extrap) THEN
   i = imin; j = jmin
   x = 0.0; y = 0.0
ENDIF

!
!4. Relative coordinates by iteration
!------------------------------------
!

IF ((frel.OR.PRESENT(weights)).AND.istat.EQ.0) THEN

   dxh1 = xdat(i+1,j) - xdat(i,j)
   dxh2 = xdat(i,j+1) - xdat(i,j)
   dxh3 = xdat(i+1,j+1) - xdat(i+1,j) - dxh2

   dyh1 = ydat(i+1,j) - ydat(i,j)
   dyh2 = ydat(i,j+1) - ydat(i,j)
   dyh3 = ydat(i+1,j+1) - ydat(i+1,j) - dyh2

   x = 0.5; y = 0.5
   
   iter_410: DO iter=1,maxit
      c(1) = xpos - xdat(i,j) - dxh1*x - dxh2*y - dxh3*x*y
      c(2) = ypos - ydat(i,j) - dyh1*x - dyh2*y - dyh3*x*y 
      a(1) = dxh1 + dxh3*y; a(2) = dyh1 + dyh3*y
      b(1) = dxh2 + dxh3*x; b(2) = dyh2 + dyh3*x
      denom = a(1)*b(2) - a(2)*b(1)
      dx = (c(1)*b(2)-c(2)*b(1))/denom
      dy = (a(1)*c(2)-a(2)*c(1))/denom
      IF (ABS(dx).LT.epsconv.AND.ABS(dy).LT.epsconv) EXIT iter_410
      x = x + dx; y = y + dy
   ENDDO iter_410

   IF (iter.GT.maxit) GOTO 1001

ENDIF

!
!5. Return grid location, relative coordinates and/or weights
!------------------------------------------------------------


IF (istat.EQ.0.OR.extrap) THEN
   icoord = i; jcoord = j
   IF (frel) THEN
      xcoord = x; ycoord = y
   ENDIF
   IF (PRESENT(weights)) THEN
      weights(1,1) = (1.0-x)*(1.0-y)
      weights(2,1) = x*(1.0-y)
      weights(1,2) = (1.0-x)*y
      weights(2,2) = x*y
   ENDIF
ENDIF


RETURN

1001 CONTINUE
nerrs = nerrs + 1
IF (errchk) THEN
   WRITE (cxpos,'(G15.7)') xpos; cxpos = ADJUSTL(cxpos)
   WRITE (cypos,'(G15.7)') ypos; cypos = ADJUSTL(cypos)
   WRITE (ioerr,'(A)') 'Iteration algorithm fails at location ('&
                     & //TRIM(cxpos)//','//TRIM(cypos)//')'
ENDIF

END SUBROUTINE hrel_coords_curv_abs_s

!========================================================================

SUBROUTINE hrel_coords_curv_abs_d(xpos,ypos,xdat,ydat,nxdat,nydat,lbounds,&
                                & ubounds,istat,inunit,icoord,jcoord,xcoord,&
                                & ycoord,weights,extrapol)
!************************************************************************
!
! *hrel_coords_curv_abs_d* calculates the horizontal relative coordinates and
!                          weight factors for interpolating data from a
!                          curvilinear grid to a given location
!                    
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description - used in case of a curvilinear grid
!             - icoord, jcoord, xcoord, ycoord and weights are undefined
!               if istat=1 and extrapol is .FALSE.
!             - xcoord and ycoord must be both  or not present
!             - xpos and ypos given in double precision
!  
! Module calls - check_cell_location, distance_minloc
!
!************************************************************************
!
USE syspars  
USE grid_routines, ONLY: check_cell_location, distance_minloc    
  
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: extrapol
INTEGER, INTENT(IN) :: inunit, nxdat, nydat
INTEGER, INTENT(OUT) :: icoord, istat, jcoord
INTEGER, DIMENSION(2), INTENT(IN) :: lbounds, ubounds  
REAL (KIND=kndrlong), INTENT(IN) :: xpos, ypos
REAL, INTENT(IN), DIMENSION(nxdat,nydat) :: xdat, ydat
REAL, INTENT(OUT), OPTIONAL :: xcoord, ycoord
REAL, INTENT(OUT), OPTIONAL, DIMENSION(2,2) :: weights

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*xpos*     REAL     X-coordinate of the input location [m or degrees longitude]
!*ypos*     REAL     Y-coordinate of the input location  [m or degrees latitude]
!*xdat*     REAL     X-coordinates of the input data grid
!                                                       [m or degrees longitude]
!*ydat*     REAL     Y-coordinate of the input data grid [m or degrees latitude]
!*nxdat*    INTEGER  X-dimension of the input data grid
!*nydat*    INTEGER  Y-dimension of the input data grid
!*lbounds*  INTEGER  Lower bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*ubounds*  INTEGER  Upper bounds of the bounding box delineating the inside
!*istat*    INTEGER  0/1 if location is inside/outside the input grid
!*inunit*   INTEGER  unit of input locations (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!*icoord*   INTEGER  Returned X-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*jcoord*   INTEGER  Returned Y-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*xcoord*   REAL     Returned normalised X-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*ycoord*   REAL     Returned normalised Y-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*weights*  REAL     Returned weight factors for interpolation
!*extrapol* LOGICAL  Allows extrapolation to the nearest inside grid point
!                    for data points outside the domain if PRESENT AND .TRUE.
!                    Default is .FALSE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: extrap, flag, frel
CHARACTER (LEN=15) :: cxpos, cypos
INTEGER, PARAMETER :: maxit = 100
INTEGER :: i, ii, imin, iter, j, jj, jmin
INTEGER, DIMENSION(2) :: locmin
REAL, PARAMETER :: epsconv = 1.0E-06
REAL :: denom, dx, dxh1, dxh2, dxh3, dy, dyh1, dyh2, dyh3, x, y
REAL, DIMENSION(2) :: a, b, c
REAL, DIMENSION(4) :: xvals, yvals


!
!1. Optional arguments
!---------------------
!

IF (PRESENT(extrapol)) THEN
   extrap = extrapol
ELSE
   extrap = .FALSE.
ENDIF

frel = PRESENT(xcoord).AND.PRESENT(ycoord)

!
!2. Search for containing grid cell
!----------------------------------
!
!---location of nearest grid point
CALL distance_minloc(xpos,ypos,xdat,ydat,nxdat,nydat,inunit,3,locmin=locmin)
imin = locmin(1); jmin = locmin(2) 

!---search loop
istat = 1
jj_210: DO jj=jmin-1,jmin
ii_210: DO ii=imin-1,imin   
   IF (ii.GE.lbounds(1).AND.ii.LT.ubounds(1).AND.&
     & jj.GE.lbounds(2).AND.jj.LT.ubounds(2)) THEN
      xvals = (/xdat(ii,jj),xdat(ii+1,jj),xdat(ii+1,jj+1),xdat(ii,jj+1)/)
      yvals = (/ydat(ii,jj),ydat(ii+1,jj),ydat(ii+1,jj+1),ydat(ii,jj+1)/)
      CALL check_cell_location(xpos,ypos,xvals,yvals,flag)
      IF (flag) THEN
         istat = 0; i = ii; j = jj
         EXIT jj_210
      ENDIF
   ENDIF
ENDDO ii_210
ENDDO jj_210

!
!3. Extrapolate if needed
!------------------------
!

IF (istat.EQ.1.AND.extrap) THEN
   i = imin; j = jmin
   x = 0.0; y = 0.0
ENDIF

!
!4. Relative coordinates by iteration
!------------------------------------
!

IF ((frel.OR.PRESENT(weights)).AND.istat.EQ.0) THEN

   dxh1 = xdat(i+1,j) - xdat(i,j)
   dxh2 = xdat(i,j+1) - xdat(i,j)
   dxh3 = xdat(i+1,j+1) - xdat(i+1,j) - dxh2

   dyh1 = ydat(i+1,j) - ydat(i,j)
   dyh2 = ydat(i,j+1) - ydat(i,j)
   dyh3 = ydat(i+1,j+1) - ydat(i+1,j) - dyh2

   x = 0.5; y = 0.5
   
   iter_410: DO iter=1,maxit
      c(1) = xpos - xdat(i,j) - dxh1*x - dxh2*y - dxh3*x*y
      c(2) = ypos - ydat(i,j) - dyh1*x - dyh2*y - dyh3*x*y 
      a(1) = dxh1 + dxh3*y; a(2) = dyh1 + dyh3*y
      b(1) = dxh2 + dxh3*x; b(2) = dyh2 + dyh3*x
      denom = a(1)*b(2) - a(2)*b(1)
      dx = (c(1)*b(2)-c(2)*b(1))/denom
      dy = (a(1)*c(2)-a(2)*c(1))/denom
      IF (ABS(dx).LT.epsconv.AND.ABS(dy).LT.epsconv) EXIT iter_410
      x = x + dx; y = y + dy
   ENDDO iter_410

   IF (iter.GT.maxit) GOTO 1001

ENDIF

!
!5. Return grid location, relative coordinates and/or weights
!------------------------------------------------------------


IF (istat.EQ.0.OR.extrap) THEN
   icoord = i; jcoord = j
   IF (frel) THEN
      xcoord = x; ycoord = y
   ENDIF
   IF (PRESENT(weights)) THEN
      weights(1,1) = (1.0-x)*(1.0-y)
      weights(2,1) = x*(1.0-y)
      weights(1,2) = (1.0-x)*y
      weights(2,2) = x*y
   ENDIF
ENDIF


RETURN

1001 CONTINUE
nerrs = nerrs + 1
IF (errchk) THEN
   WRITE (cxpos,'(G15.7)') xpos; cxpos = ADJUSTL(cxpos)
   WRITE (cypos,'(G15.7)') ypos; cypos = ADJUSTL(cypos)
   WRITE (ioerr,'(A)') 'Iteration algorithm fails at location ('&
                     & //TRIM(cxpos)//','//TRIM(cypos)//')'
ENDIF

END SUBROUTINE hrel_coords_curv_abs_d

!========================================================================

SUBROUTINE hrel_coords_curv_rel_s(xpos,ypos,displx,disply,xdat,ydat,nxdat,&
                                & nydat,lbounds,ubounds,lsize,istat,inunit,&
                                & icoord,jcoord,xcoord,ycoord,weights,&
                                & extrapol,update)
!************************************************************************
!
! *hrel_coords_curv_rel_s* updates the horizontal relative coordinates and
!                          weight factors after applying a given displacement
!                          to the old location
!                    
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description - used in case of a curvilinear grid
!             - icoord, jcoord, xcoord, ycoord and weights are undefined
!               if istat=1 and extrapol is .FALSE.
!             - xcoord and ycoord must be both  or not present
!             - xpos and ypos given in single precision  
!
! Module calls - check_cell_location, hrel_coords_curv_abs
!  
!************************************************************************
!
USE syspars  
USE grid_routines, ONLY: check_cell_location  
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: extrapol, update
INTEGER, INTENT(IN) :: inunit, lsize, nxdat, nydat
INTEGER, INTENT(INOUT) :: icoord, jcoord
INTEGER, INTENT(OUT) :: istat
INTEGER, DIMENSION(2), INTENT(IN) :: lbounds, ubounds  
REAL, INTENT(IN) :: displx, disply
REAL (KIND=kndreal), INTENT(INOUT) :: xpos, ypos
REAL, INTENT(IN), DIMENSION(nxdat,nydat) :: xdat, ydat
REAL, INTENT(OUT), OPTIONAL :: xcoord, ycoord
REAL, INTENT(OUT), OPTIONAL, DIMENSION(2,2) :: weights

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*xpos*     REAL     X-coordinate of the input location [m or degrees longitude]
!*ypos*     REAL     Y-coordinate of the input location  [m or degrees latitude]
!*displx*   REAL     Applied displacement in the X-direction
!                                                       [m or degrees longitude]
!*disply*   REAL     Applied displacement in the Y-direction
!                                                        [m or degrees latitude]
!*xdat*     REAL     X-coordinates of the input data grid
!                                                       [m or degrees longitude]
!*ydat*     REAL     Y-coordinate of the input data grid [m or degrees latitude]
!*nxdat*    INTEGER  X-dimension of the input data grid
!*nydat*    INTEGER  Y-dimension of the input data grid
!*lbounds*  INTEGER  Lower bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*ubounds*  INTEGER  Upper bounds of the bounding box delineating the inside
!*lsize*    INTEGER  Distance in grid cells to the current cell of the
!                    surrounding cell where the search is first performed 
!*istat*    INTEGER  0/1 if location is inside/outside the input grid
!*inunit*   INTEGER  unit of input locations (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!*icoord*   INTEGER  Returned X-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*jcoord*   INTEGER  Returned Y-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*xcoord*   REAL     Returned normalised X-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*ycoord*   REAL     Returned normalised Y-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*weights*  REAL     Returned weight factors for interpolation
!*extrapol* LOGICAL  Allows extrapolation to the nearest inside grid point
!                    for data points outside the domain if PRESENT AND .TRUE.
!                    Default is .FALSE.
!*update*   LOGICAL  If PRESENT and .TRUE., the displaced locations are returned
!                    in xpos and ypos
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: extrap, flag, frel, updatex
CHARACTER (LEN=15) :: cxpos, cypos
INTEGER, PARAMETER :: maxit = 100
INTEGER :: i, ii, iter, j, jj, l
REAL, PARAMETER :: epsconv = 1.0E-06
REAL :: denom, dx, dxh1, dxh2, dxh3, dy, dyh1, dyh2, &
      & dyh3, x, y
REAL (KIND=kndreal) :: xp, yp
REAL, DIMENSION(2) :: a, b, c
REAL, DIMENSION(4) :: xvals, yvals


!
!1. Optional arguments
!---------------------
!

IF (PRESENT(extrapol)) THEN
   extrap = extrapol
ELSE
   extrap = .FALSE.
ENDIF

IF (PRESENT(update)) THEN
   updatex = update
ELSE
   updatex = .FALSE.
ENDIF

frel = PRESENT(xcoord).AND.PRESENT(ycoord)

!
!2. Apply search criterion to the current and (if needed) to surroundig cells
!----------------------------------------------------------------------------
!

i = icoord; j = jcoord
xp = xpos + displx; yp = ypos + disply

!
!2.1 Current cell
!----------------
!

xvals = (/xdat(i,j),xdat(i+1,j),xdat(i+1,j+1),xdat(i,j+1)/)
yvals = (/ydat(i,j),ydat(i+1,j),ydat(i+1,j+1),ydat(i,j+1)/)
flag = .TRUE.
IF (i.GE.lbounds(1).AND.i.LT.ubounds(1).AND.&
  & j.GE.lbounds(2).AND.j.LT.ubounds(2)) THEN
   xvals = (/xdat(i,j),xdat(i+1,j),xdat(i+1,j+1),xdat(i,j+1)/)
   yvals = (/ydat(i,j),ydat(i+1,j),ydat(i+1,j+1),ydat(i,j+1)/)
   CALL check_cell_location(xp,yp,xvals,yvals,flag)
ENDIF

!
!2.2 Neighbouring cells
!----------------------
!

l_220: DO l=1,lsize

!  ---cells to the South/North
   IF (.NOT.flag) THEN
      jj_221: DO jj=j-l,j+l,2*l
      ii_221: DO ii=i-l,i+l
         IF (ii.GE.lbounds(1).AND.ii.LT.ubounds(1).AND.&
           & jj.GE.lbounds(2).AND.jj.LT.ubounds(2)) THEN
            xvals = (/xdat(ii,jj),xdat(ii+1,jj),xdat(ii+1,jj+1),xdat(ii,jj+1)/)
            yvals = (/ydat(ii,jj),ydat(ii+1,jj),ydat(ii+1,jj+1),&
                    & ydat(ii,jj+1)/)  
            CALL check_cell_location(xp,yp,xvals,yvals,flag)
            IF (flag) THEN
               i = ii; j = jj
               EXIT l_220
            ENDIF
         ENDIF
      ENDDO ii_221
      ENDDO jj_221
   ENDIF

!  ---cells to the West/East
   IF (.NOT.flag) THEN
      ii_222: DO ii=i-l,i+l,2*l
      jj_222: DO jj=j-l+1,j+l-1
         IF (ii.GE.lbounds(1).AND.ii.LT.ubounds(1).AND.&
           & jj.GE.lbounds(2).AND.jj.LT.ubounds(2)) THEN
            xvals = (/xdat(ii,jj),xdat(ii+1,jj),xdat(ii+1,jj+1),xdat(ii,jj+1)/)
            yvals = (/ydat(ii,jj),ydat(ii+1,jj),ydat(ii+1,jj+1),&
                    & ydat(ii,jj+1)/)  
            CALL check_cell_location(xp,yp,xvals,yvals,flag)
            IF (flag) THEN
               i = ii; j = jj
               EXIT l_220
            ENDIF
         ENDIF
      ENDDO jj_222
      ENDDO ii_222
   ENDIF

ENDDO l_220

!
!3. Extend search to the whole domain if needed
!----------------------------------------------
!

IF (flag) THEN
   istat = 0
ELSE
   CALL hrel_coords_curv_abs(xp,yp,xdat,ydat,nxdat,nydat,lbounds,ubounds,istat,&
                           & inunit,i,j,extrapol=extrap)
   IF (istat.EQ.1) THEN
      x = 0.0; y = 0.0
   ENDIF
ENDIF

!
!4. Relative coordinates by iteration
!------------------------------------
!

IF ((frel.OR.PRESENT(weights)).AND.istat.EQ.0) THEN

   dxh1 = xdat(i+1,j) - xdat(i,j)
   dxh2 = xdat(i,j+1) - xdat(i,j)
   dxh3 = xdat(i+1,j+1) - xdat(i+1,j) - dxh2

   dyh1 = ydat(i+1,j) - ydat(i,j)
   dyh2 = ydat(i,j+1) - ydat(i,j)
   dyh3 = ydat(i+1,j+1) - ydat(i+1,j) - dyh2

   x = 0.5; y = 0.5
   
   iter_410: DO iter=1,maxit
      c(1) = xp - xdat(i,j) - dxh1*x - dxh2*y - dxh3*x*y
      c(2) = yp - ydat(i,j) - dyh1*x - dyh2*y - dyh3*x*y 
      a(1) = dxh1 + dxh3*y; a(2) = dyh1 + dyh3*y
      b(1) = dxh2 + dxh3*x; b(2) = dyh2 + dyh3*x
      denom = a(1)*b(2) - a(2)*b(1)
      dx = (c(1)*b(2)-c(2)*b(1))/denom
      dy = (a(1)*c(2)-a(2)*c(1))/denom
      IF (ABS(dx).LT.epsconv.AND.ABS(dy).LT.epsconv) EXIT iter_410
      x = x + dx; y = y + dy
   ENDDO iter_410

   IF (iter.GT.maxit) GOTO 1001

ENDIF

!
!5. Return grid location, relative coordinates and/or weights
!------------------------------------------------------------


IF (istat.EQ.0.OR.extrap) THEN
   icoord = i; jcoord = j
   IF (frel) THEN
      xcoord = x; ycoord = y
   ENDIF
   IF (PRESENT(weights)) THEN
      weights(1,1) = (1.0-x)*(1.0-y)
      weights(2,1) = x*(1.0-y)
      weights(1,2) = (1.0-x)*y
      weights(2,2) = x*y
   ENDIF
   IF (updatex) THEN
      xpos = xp; ypos = yp
   ENDIF
ENDIF


RETURN

1001 CONTINUE
nerrs = nerrs + 1
IF (errchk) THEN
   WRITE (cxpos,'(G15.7)') xpos; cxpos = ADJUSTL(cxpos)
   WRITE (cypos,'(G15.7)') ypos; cypos = ADJUSTL(cypos)
   WRITE (ioerr,'(A)') 'Iteration algorithm fails at location ('&
                     & //TRIM(cxpos)//','//TRIM(cypos)//')'
ENDIF

END SUBROUTINE hrel_coords_curv_rel_s

!========================================================================

SUBROUTINE hrel_coords_curv_rel_d(xpos,ypos,displx,disply,xdat,ydat,nxdat,&
                                & nydat,lbounds,ubounds,lsize,istat,inunit,&
                                & icoord,jcoord,xcoord,ycoord,weights,extrapol,&
                                & update)
!************************************************************************
!
! *hrel_coords_curv_rel_d* updates the horizontal relative coordinates and
!                          weight factors after applying a given displacement
!                          to the old location
!                    
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description - used in case of a curvilinear grid
!             - icoord, jcoord, xcoord, ycoord and weights are undefined
!               if istat=1 and extrapol is .FALSE.
!             - xcoord and ycoord must be both  or not present
!             - xpos and ypos given in double precision  
!
! Module calls - check_cell_location, hrel_coords_curv_abs
!  
!************************************************************************
!
USE syspars  
USE grid_routines, ONLY: check_cell_location  
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: extrapol, update
INTEGER, INTENT(IN) :: inunit, lsize, nxdat, nydat
INTEGER, INTENT(INOUT) :: icoord, jcoord
INTEGER, INTENT(OUT) :: istat
INTEGER, DIMENSION(2), INTENT(IN) :: lbounds, ubounds  
REAL, INTENT(IN) :: displx, disply
REAL (KIND=kndrlong), INTENT(INOUT) :: xpos, ypos
REAL, INTENT(IN), DIMENSION(nxdat,nydat) :: xdat, ydat
REAL, INTENT(OUT), OPTIONAL :: xcoord, ycoord
REAL, INTENT(OUT), OPTIONAL, DIMENSION(2,2) :: weights

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*xpos*     REAL     X-coordinate of the input location [m or degrees longitude]
!*ypos*     REAL     Y-coordinate of the input location  [m or degrees latitude]
!*displx*   REAL     Applied displacement in the X-direction
!                                                       [m or degrees longitude]
!*disply*   REAL     Applied displacement in the Y-direction
!                                                        [m or degrees latitude]
!*xdat*     REAL     X-coordinates of the input data grid
!                                                       [m or degrees longitude]
!*ydat*     REAL     Y-coordinate of the input data grid [m or degrees latitude]
!*nxdat*    INTEGER  X-dimension of the input data grid
!*nydat*    INTEGER  Y-dimension of the input data grid
!*lbounds*  INTEGER  Lower bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*ubounds*  INTEGER  Upper bounds of the bounding box delineating the inside
!*lsize*    INTEGER  Distance in grid cells to the current cell of the
!                    surrounding cell where the search is first performed 
!*istat*    INTEGER  0/1 if location is inside/outside the input grid
!*inunit*   INTEGER  unit of input locations (spherical grid only)
!                   = 1 => radians
!                   = 2 => degrees
!*icoord*   INTEGER  Returned X-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*jcoord*   INTEGER  Returned Y-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*xcoord*   REAL     Returned normalised X-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*ycoord*   REAL     Returned normalised Y-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*weights*  REAL     Returned weight factors for interpolation
!*extrapol* LOGICAL  Allows extrapolation to the nearest inside grid point
!                    for data points outside the domain if PRESENT AND .TRUE.
!                    Default is .FALSE.
!*update*   LOGICAL  If PRESENT and .TRUE., the displaced locations are returned
!                    in xpos and ypos
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: extrap, flag, frel, updatex
CHARACTER (LEN=15) :: cxpos, cypos
INTEGER, PARAMETER :: maxit = 100
INTEGER :: i, ii, iter, j, jj, l
REAL, PARAMETER :: epsconv = 1.0E-06
REAL :: denom, dx, dxh1, dxh2, dxh3, dy, dyh1, dyh2, dyh3, x, y
REAL (KIND=kndrlong) :: xp, yp
REAL, DIMENSION(2) :: a, b, c
REAL, DIMENSION(4) :: xvals, yvals


!
!1. Optional arguments
!---------------------
!

IF (PRESENT(extrapol)) THEN
   extrap = extrapol
ELSE
   extrap = .FALSE.
ENDIF

IF (PRESENT(update)) THEN
   updatex = update
ELSE
   updatex = .FALSE.
ENDIF

frel = PRESENT(xcoord).AND.PRESENT(ycoord)

!
!2. Apply search criterion to the current and (if needed) to surroundig cells
!----------------------------------------------------------------------------
!

i = icoord; j = jcoord
xp = xpos + displx; yp = ypos + disply

!
!2.1 Current cell
!----------------
!

xvals = (/xdat(i,j),xdat(i+1,j),xdat(i+1,j+1),xdat(i,j+1)/)
yvals = (/ydat(i,j),ydat(i+1,j),ydat(i+1,j+1),ydat(i,j+1)/)
flag = .TRUE.
IF (i.GE.lbounds(1).AND.i.LT.ubounds(1).AND.&
  & j.GE.lbounds(2).AND.j.LT.ubounds(2)) THEN
   xvals = (/xdat(i,j),xdat(i+1,j),xdat(i+1,j+1),xdat(i,j+1)/)
   yvals = (/ydat(i,j),ydat(i+1,j),ydat(i+1,j+1),ydat(i,j+1)/)
   CALL check_cell_location(xp,yp,xvals,yvals,flag)
ENDIF

!
!2.2 Neighbouring cells
!----------------------
!

l_220: DO l=1,lsize

!  ---cells to the South/North
   IF (.NOT.flag) THEN
      jj_221: DO jj=j-l,j+l,2*l
      ii_221: DO ii=i-l,i+l
         IF (ii.GE.lbounds(1).AND.ii.LT.ubounds(1).AND.&
           & jj.GE.lbounds(2).AND.jj.LT.ubounds(2)) THEN
            xvals = (/xdat(ii,jj),xdat(ii+1,jj),xdat(ii+1,jj+1),xdat(ii,jj+1)/)
            yvals = (/ydat(ii,jj),ydat(ii+1,jj),ydat(ii+1,jj+1),&
                    & ydat(ii,jj+1)/)  
            CALL check_cell_location(xp,yp,xvals,yvals,flag)
            IF (flag) THEN
               i = ii; j = jj
               EXIT l_220
            ENDIF
         ENDIF
      ENDDO ii_221
      ENDDO jj_221
   ENDIF

!  ---cells to the West/East
   IF (.NOT.flag) THEN
      ii_222: DO ii=i-l,i+l,2*l
      jj_222: DO jj=j-l+1,j+l-1
         IF (ii.GE.lbounds(1).AND.ii.LT.ubounds(1).AND.&
           & jj.GE.lbounds(2).AND.jj.LT.ubounds(2)) THEN
            xvals = (/xdat(ii,jj),xdat(ii+1,jj),xdat(ii+1,jj+1),xdat(ii,jj+1)/)
            yvals = (/ydat(ii,jj),ydat(ii+1,jj),ydat(ii+1,jj+1),&
                    & ydat(ii,jj+1)/)  
            CALL check_cell_location(xp,yp,xvals,yvals,flag)
            IF (flag) THEN
               i = ii; j = jj
               EXIT l_220
            ENDIF
         ENDIF
      ENDDO jj_222
      ENDDO ii_222
   ENDIF

ENDDO l_220

!
!3. Extend search to the whole domain if needed
!----------------------------------------------
!

IF (flag) THEN
   istat = 0
ELSE
   CALL hrel_coords_curv_abs(xp,yp,xdat,ydat,nxdat,nydat,lbounds,ubounds,istat,&
                           & inunit,i,j,extrapol=extrap)
   IF (istat.EQ.1) THEN
      x = 0.0; y = 0.0
   ENDIF
ENDIF

!
!4. Relative coordinates by iteration
!------------------------------------
!

IF ((frel.OR.PRESENT(weights)).AND.istat.EQ.0) THEN

   dxh1 = xdat(i+1,j) - xdat(i,j)
   dxh2 = xdat(i,j+1) - xdat(i,j)
   dxh3 = xdat(i+1,j+1) - xdat(i+1,j) - dxh2

   dyh1 = ydat(i+1,j) - ydat(i,j)
   dyh2 = ydat(i,j+1) - ydat(i,j)
   dyh3 = ydat(i+1,j+1) - ydat(i+1,j) - dyh2

   x = 0.5; y = 0.5
   
   iter_410: DO iter=1,maxit
      c(1) = xp - xdat(i,j) - dxh1*x - dxh2*y - dxh3*x*y
      c(2) = yp - ydat(i,j) - dyh1*x - dyh2*y - dyh3*x*y 
      a(1) = dxh1 + dxh3*y; a(2) = dyh1 + dyh3*y
      b(1) = dxh2 + dxh3*x; b(2) = dyh2 + dyh3*x
      denom = a(1)*b(2) - a(2)*b(1)
      dx = (c(1)*b(2)-c(2)*b(1))/denom
      dy = (a(1)*c(2)-a(2)*c(1))/denom
      IF (ABS(dx).LT.epsconv.AND.ABS(dy).LT.epsconv) EXIT iter_410
      x = x + dx; y = y + dy
   ENDDO iter_410

   IF (iter.GT.maxit) GOTO 1001

ENDIF

!
!5. Return grid location, relative coordinates and/or weights
!------------------------------------------------------------


IF (istat.EQ.0.OR.extrap) THEN
   icoord = i; jcoord = j
   IF (frel) THEN
      xcoord = x; ycoord = y
   ENDIF
   IF (PRESENT(weights)) THEN
      weights(1,1) = (1.0-x)*(1.0-y)
      weights(2,1) = x*(1.0-y)
      weights(1,2) = (1.0-x)*y
      weights(2,2) = x*y
   ENDIF
   IF (updatex) THEN
      xpos = xp; ypos = yp
   ENDIF
ENDIF


RETURN

1001 CONTINUE
nerrs = nerrs + 1
IF (errchk) THEN
   WRITE (cxpos,'(G15.7)') xpos; cxpos = ADJUSTL(cxpos)
   WRITE (cypos,'(G15.7)') ypos; cypos = ADJUSTL(cypos)
   WRITE (ioerr,'(A)') 'Iteration algorithm fails at location ('&
                     & //TRIM(cxpos)//','//TRIM(cypos)//')'
ENDIF

END SUBROUTINE hrel_coords_curv_rel_d

!========================================================================

SUBROUTINE hrel_coords_rect_nonunif_abs_s(xpos,ypos,xdat,ydat,nxdat,nydat,&
                                        & lbounds,ubounds,istat,icoord,jcoord,&
                                        & xcoord,ycoord,weights,extrapol)
!************************************************************************
!
! *hrel_coords_nonunif_abs_s* calculates the horizontal relative coordinates and
!                             weight factors for interpolating data from a
!                             non-uniform rectangular grid to a given location
!                    
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description - used in case of a non-uniform rectangular grid
!             - icoord, jcoord, xcoord, ycoord and weights are undefined
!               if istat=1 and extrapol is .FALSE.
!             - xcoord and ycoord must be both  or not present  
!             - xpos and ypos given in single precision  
!
!************************************************************************
!
USE syspars
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: extrapol
INTEGER, INTENT(IN) :: nxdat, nydat
INTEGER, INTENT(OUT) :: icoord, istat, jcoord
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds, ubounds  
REAL (KIND=kndreal), INTENT(IN) :: xpos, ypos
REAL, INTENT(IN), DIMENSION(nxdat) :: xdat
REAL, INTENT(IN), DIMENSION(nydat) :: ydat
REAL, INTENT(OUT), OPTIONAL :: xcoord, ycoord
REAL, INTENT(OUT), OPTIONAL, DIMENSION(2,2) :: weights

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*xpos*     REAL     X-coordinate of the input location [m or degrees longitude]
!*ypos*     REAL     Y-coordinate of the input location  [m or degrees latitude]
!*xdat*     REAL     X-coordinates of the input data grid
!                                                       [m or degrees longitude]
!*ydat*     REAL     Y-coordinate of the input data grid [m or degrees latitude]
!*nxdat*    INTEGER  X-dimension of the input data grid
!*nydat*    INTEGER  Y-dimension of the input data grid
!*lbounds*  INTEGER  Lower bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*ubounds*  INTEGER  Upper bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*istat*    INTEGER  0/1 if location is inside/outside the input grid
!*icoord*   INTEGER  Returned X-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*jcoord*   INTEGER  Returned Y-index coordinate of the location with respect
!                    to the  lower left corner of the containing cell in the
!                    input grid
!*xcoord*   REAL     Returned normalised X-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*ycoord*     REAL   Returned normalised Y-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*weights*  REAL     Returned weight factors for interpolation
!*extrapol* LOGICAL  Allows extrapolation to the nearest inside grid point
!                    for data points outside the domain if PRESENT AND .TRUE.
!                    Default is .FALSE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: extrap, frel, search
INTEGER :: i, j
REAL :: x, y


!
!1. Optional arguments
!---------------------
!

IF (PRESENT(extrapol)) THEN
   extrap = extrapol
ELSE
   extrap = .FALSE.
ENDIF

frel = PRESENT(xcoord).AND.PRESENT(ycoord)

!
!2. Search for containing grid cell
!----------------------------------
!
!---X-direction
i = lbounds(1)-1; x = 0.0
IF (xpos.GE.xdat(lbounds(1))) THEN
   i = lbounds(1); search = .TRUE.
ELSE
   i = lbounds(1)-1; search = .FALSE.
ENDIF
DO WHILE (search.AND.i.LT.ubounds(1))
   IF (xpos.GE.xdat(i).AND.xpos.LT.xdat(i+1)) THEN
      search = .FALSE.
      x = (xpos-xdat(i))/(xdat(i+1)-xdat(i))
   ELSE
      i = i + 1
   ENDIF
ENDDO
IF (.NOT.search) x = MIN(x,1.0)

!---Y-direction
j = lbounds(2)-1; y = 0.0
IF (ypos.GE.ydat(lbounds(2))) THEN
   j = lbounds(2); search = .TRUE.
ELSE
   j = lbounds(2)-1; search = .FALSE.
ENDIF
DO WHILE (search.AND.j.LT.ubounds(2))
   IF (ypos.GE.ydat(j).AND.ypos.LT.ydat(j+1)) THEN
      search = .FALSE.
      y = (ypos-ydat(j))/(ydat(j+1)-ydat(j))
   ELSE
      j = j + 1
   ENDIF
ENDDO
IF (.NOT.search) y = MIN(y,1.0)

!
!3. Check out of domain condition
!--------------------------------
!

IF (i.LT.lbounds(1).OR.i.GT.ubounds(1).OR.&
  & j.LT.lbounds(2).OR.j.GT.ubounds(2)) THEN
   istat = 1
ELSE
   istat = 0
ENDIF

!
!4. Extrapolate if needed
!------------------------
!

IF (istat.EQ.1.AND.extrap) THEN
   IF (i.LT.lbounds(1)) THEN
      i = lbounds(1); x = 0.0
   ELSEIF (i.GE.ubounds(1)) THEN
      i = ubounds(1)-1; x = 1.0
   ENDIF
   IF (j.LT.lbounds(2)) THEN
      j = lbounds(2); y = 0.0
   ELSEIF (j.GE.ubounds(2)) THEN
      j = ubounds(2)-1; y = 1.0
   ENDIF
ENDIF

!
!5. Relative coordinates and/or weights
!--------------------------------------
!

IF (istat.EQ.0.OR.extrap) THEN
   icoord = i; jcoord = j
   IF (frel) THEN
      xcoord = x; ycoord = y
   ENDIF
   IF (PRESENT(weights)) THEN
      weights(1,1) = (1.0-x)*(1.0-y)
      weights(2,1) = x*(1.0-y)
      weights(1,2) = (1.0-x)*y
      weights(2,2) = x*y
   ENDIF
ENDIF


RETURN

END SUBROUTINE hrel_coords_rect_nonunif_abs_s

!========================================================================

SUBROUTINE hrel_coords_rect_nonunif_abs_d(xpos,ypos,xdat,ydat,nxdat,nydat,&
                                        & lbounds,ubounds,istat,icoord,jcoord,&
                                        & xcoord,ycoord,weights,extrapol)
!************************************************************************
!
! *hrel_coords_nonunif_abs_d* calculates the horizontal relative coordinates and
!                             weight factors for interpolating data from a
!                             non-uniform rectangular grid to a given location
!                    
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description - used in case of a non-uniform rectangular grid
!             - icoord, jcoord, xcoord, ycoord and weights are undefined
!               if istat=1 and extrapol is .FALSE.
!             - xcoord and ycoord must be both  or not present  
!             - xpos and ypos given in double precision  
!
!************************************************************************
!
USE syspars
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: extrapol
INTEGER, INTENT(IN) :: nxdat, nydat
INTEGER, INTENT(OUT) :: icoord, istat, jcoord
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds, ubounds  
REAL (KIND=kndrlong), INTENT(IN) :: xpos, ypos
REAL, INTENT(IN), DIMENSION(nxdat) :: xdat
REAL, INTENT(IN), DIMENSION(nydat) :: ydat
REAL, INTENT(OUT), OPTIONAL :: xcoord, ycoord
REAL, INTENT(OUT), OPTIONAL, DIMENSION(2,2) :: weights

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*xpos*     REAL     X-coordinate of the input location [m or degrees longitude]
!*ypos*     REAL     Y-coordinate of the input location  [m or degrees latitude]
!*xdat*     REAL     X-coordinates of the input data grid
!                                                       [m or degrees longitude]
!*ydat*     REAL     Y-coordinate of the input data grid [m or degrees latitude]
!*nxdat*    INTEGER  X-dimension of the input data grid
!*nydat*    INTEGER  Y-dimension of the input data grid
!*lbounds*  INTEGER  Lower bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*ubounds*  INTEGER  Upper bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*istat*    INTEGER  0/1 if location is inside/outside the input grid
!*icoord*   INTEGER  Returned X-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*jcoord*   INTEGER  Returned Y-index coordinate of the location with respect
!                    to the  lower left corner of the containing cell in the
!                    input grid
!*xcoord*   REAL     Returned normalised X-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*ycoord*   REAL     Returned normalised Y-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*weights*  REAL     Returned weight factors for interpolation
!*extrapol* LOGICAL  Allows extrapolation to the nearest inside grid point
!                    for data points outside the domain if PRESENT AND .TRUE.
!                    Default is .FALSE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: extrap, frel, search
INTEGER :: i, j
REAL :: x, y


!
!1. Optional arguments
!---------------------
!

IF (PRESENT(extrapol)) THEN
   extrap = extrapol
ELSE
   extrap = .FALSE.
ENDIF

frel = PRESENT(xcoord).AND.PRESENT(ycoord)

!
!2. Search for containing grid cell
!----------------------------------
!
!---X-direction
i = lbounds(1)-1; x = 0.0
IF (xpos.GE.xdat(lbounds(1))) THEN
   i = lbounds(1); search = .TRUE.
ELSE
   i = lbounds(1)-1; search = .FALSE.
ENDIF
DO WHILE (search.AND.i.LT.ubounds(1))
   IF (xpos.GE.xdat(i).AND.xpos.LT.xdat(i+1)) THEN
      search = .FALSE.
      x = (xpos-xdat(i))/(xdat(i+1)-xdat(i))
   ELSE
      i = i + 1
   ENDIF
ENDDO
IF (.NOT.search) x = MIN(x,1.0)

!---Y-direction
j = lbounds(2)-1; y = 0.0
IF (ypos.GE.ydat(lbounds(2))) THEN
   j = lbounds(2); search = .TRUE.
ELSE
   j = lbounds(2)-1; search = .FALSE.
ENDIF
DO WHILE (search.AND.j.LT.ubounds(2))
   IF (ypos.GE.ydat(j).AND.ypos.LT.ydat(j+1)) THEN
      search = .FALSE.
      y = (ypos-ydat(j))/(ydat(j+1)-ydat(j))
   ELSE
      j = j + 1
   ENDIF
ENDDO
IF (.NOT.search) y = MIN(y,1.0)

!
!3. Check out of domain condition
!--------------------------------
!

IF (i.LT.lbounds(1).OR.i.GT.ubounds(1).OR.&
  & j.LT.lbounds(2).OR.j.GT.ubounds(2)) THEN
   istat = 1
ELSE
   istat = 0
ENDIF

!
!4. Extrapolate if needed
!------------------------
!

IF (istat.EQ.1.AND.extrap) THEN
   IF (i.LT.lbounds(1)) THEN
      i = lbounds(1); x = 0.0
   ELSEIF (i.GE.ubounds(1)) THEN
      i = ubounds(1)-1; x = 1.0
   ENDIF
   IF (j.LT.lbounds(2)) THEN
      j = lbounds(2); y = 0.0
   ELSEIF (j.GE.ubounds(2)) THEN
      j = ubounds(2)-1; y = 1.0
   ENDIF
ENDIF

!
!5. Relative coordinates and/or weights
!--------------------------------------
!

IF (istat.EQ.0.OR.extrap) THEN
   icoord = i; jcoord = j
   IF (frel) THEN
      xcoord = x; ycoord = y
   ENDIF
   IF (PRESENT(weights)) THEN
      weights(1,1) = (1.0-x)*(1.0-y)
      weights(2,1) = x*(1.0-y)
      weights(1,2) = (1.0-x)*y
      weights(2,2) = x*y
   ENDIF
ENDIF


RETURN

END SUBROUTINE hrel_coords_rect_nonunif_abs_d

!========================================================================

SUBROUTINE hrel_coords_rect_nonunif_rel_s(xpos,ypos,displx,disply,xdat,ydat,&
                                        & nxdat,nydat,lbounds,ubounds,istat,&
                                        & icoord,jcoord,xcoord,ycoord,weights,&
                                        & extrapol,update)
!************************************************************************
!
! *hrel_coords_nonunif_rel_s* updates the horizontal relative coordinates and
!                             weight factors after applying a given displacement
!                             to the old location
!                     
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description - used in case of a non-uniform rectangular grid
!             - icoord, jcoord, xcoord, ycoord and weights are undefined
!               if istat=1 and extrapol is .FALSE.
!             - xcoord and ycoord must be both  or not present  
!             - xpos and ypos given in single precision  
!
!************************************************************************
!
USE syspars
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: extrapol, update
INTEGER, INTENT(IN) :: nxdat, nydat
INTEGER, INTENT(INOUT) :: icoord, jcoord
INTEGER, INTENT(OUT) :: istat
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds, ubounds  
REAL, INTENT(IN) :: displx, disply
REAL (KIND=kndreal), INTENT(INOUT) :: xpos, ypos
REAL, INTENT(IN), DIMENSION(nxdat) :: xdat
REAL, INTENT(IN), DIMENSION(nydat) :: ydat
REAL, INTENT(OUT), OPTIONAL :: xcoord, ycoord
REAL, INTENT(OUT), OPTIONAL, DIMENSION(2,2) :: weights

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*xpos*     REAL     X-coordinate of the input location [m or degrees longitude]
!*ypos*     REAL     Y-coordinate of the input location  [m or degrees latitude]
!*displx*   REAL     Applied displacement in the X-direction
!                                                       [m or degrees longitude]
!*disply*   REAL     Applied displacement in the Y-direction
!                                                        [m or degrees latitude]
!*xdat*     REAL     X-coordinates of the input data grid
!                                                       [m or degrees longitude]
!*ydat*     REAL     Y-coordinate of the input data grid [m or degrees latitude]
!*nxdat*    INTEGER  X-dimension of the input data grid
!*nydat*    INTEGER  Y-dimension of the input data grid
!*lbounds*  INTEGER  Lower bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*ubounds*  INTEGER  Upper bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*istat*    INTEGER  0/1 if location is inside/outside the input grid
!*icoord*   INTEGER  Returned X-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*jcoord*   INTEGER  Returned Y-index coordinate of the location with respect
!                    to the  lower left corner of the containing cell in the
!                    input grid
!*xcoord*   REAL     Returned normalised X-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*ycoord*   REAL     Returned normalised Y-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*weights*  REAL     Returned weight factors for interpolation
!*extrapol* LOGICAL  Allows extrapolation to the nearest inside grid point
!                    for data points outside the domain if PRESENT AND .TRUE.
!                    Default is .FALSE.
!*update*   LOGICAL  If PRESENT and .TRUE., the displaced locations are returned
!                    in xpos and ypos
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: extrap, frel, search, updatex
INTEGER :: i, j
REAL :: x, y
REAL (KIND=kndreal) :: xp, yp


!
!1. Optional arguments
!---------------------
!

IF (PRESENT(extrapol)) THEN
   extrap = extrapol
ELSE
   extrap = .FALSE.
ENDIF

IF (PRESENT(update)) THEN
   updatex = update
ELSE
   updatex = .FALSE.
ENDIF

frel = PRESENT(xcoord).AND.PRESENT(ycoord)

!
!2. Search for containing grid cell
!----------------------------------
!
!2.1 X-direction
!---------------
!
!---eastward displacement
xp = xpos + displx
IF (displx.GE.0.0) THEN
   IF (xp.GT.xdat(ubounds(1))) THEN
      search = .FALSE.; i = ubounds(1)
   ELSE
      search = .TRUE.; i = icoord
   ENDIF
   DO WHILE (search.AND.i.LT.ubounds(1))
      IF (xp.GE.xdat(i).AND.xp.LT.xdat(i+1)) THEN
         search = .FALSE.
         x = (xp-xdat(i))/(xdat(i+1)-xdat(i))
      ELSE
         i = i + 1
      ENDIF
   ENDDO
   IF (.NOT.search) x = MIN(x,1.0)

!---westward displacement
ELSE
   IF (xp.LT.xdat(lbounds(1))) THEN
      search = .FALSE.; i = lbounds(1) - 1
   ELSE
      search = .TRUE.; i = icoord
   ENDIF
   DO WHILE (search.AND.i.GE.lbounds(1))
      IF (xp.GE.xdat(i).AND.xp.LT.xdat(i+1)) THEN
         search = .FALSE.
         x = (xp-xdat(i))/(xdat(i+1)-xdat(i))
      ELSE
         i = i - 1
      ENDIF
   ENDDO
   IF (.NOT.search) x = MIN(x,1.0)
ENDIF

!
!2.2 Y-direction
!---------------
!
!---northtward displacement
yp = ypos + disply
IF (disply.GE.0.0) THEN
   IF (yp.GT.ydat(ubounds(2))) THEN
      search = .FALSE.; j = ubounds(2)
   ELSE
      search = .TRUE.; j = jcoord
   ENDIF
   DO WHILE (search.AND.j.LT.ubounds(2))
      IF (yp.GE.ydat(j).AND.yp.LT.ydat(j+1)) THEN
         search = .FALSE.
         y = (yp-ydat(j))/(ydat(j+1)-ydat(j))
      ELSE
         j = j + 1
      ENDIF
   ENDDO
   IF (.NOT.search) y = MIN(y,1.0)

!---southward displacement
ELSE
   IF (yp.LT.ydat(lbounds(2))) THEN
      search = .FALSE.; j = lbounds(2) - 1
   ELSE
      search = .TRUE.; j = jcoord
   ENDIF
   DO WHILE (search.AND.j.GE.lbounds(2))
      IF (yp.GE.ydat(j).AND.yp.LT.ydat(j+1)) THEN
         search = .FALSE.
         y = (yp-ydat(j))/(ydat(j+1)-ydat(j))
      ELSE
         j = j - 1
      ENDIF
   ENDDO
   IF (.NOT.search) y = MIN(y,1.0)
ENDIF

!
!3. Check out of domain condition
!--------------------------------
!

IF (i.LT.lbounds(1).OR.i.GT.ubounds(1).OR.&
  & j.LT.lbounds(2).OR.j.GT.ubounds(2)) THEN
   istat = 1
ELSE
   istat = 0
ENDIF

!
!4. Extrapolate if needed
!------------------------
!

IF (istat.EQ.1.AND.extrap) THEN
   IF (i.LT.lbounds(1)) THEN
      i = lbounds(1); x = 0.0
   ELSEIF (i.GE.ubounds(1)) THEN
      i = ubounds(1)-1; x = 1.0
   ENDIF
   IF (j.LT.lbounds(2)) THEN
      j = lbounds(2); y = 0.0
   ELSEIF (j.GE.ubounds(2)) THEN
      j = ubounds(2)-1; y = 1.0
   ENDIF
ENDIF

!
!5. Relative coordinates and/or weights
!--------------------------------------
!

IF (istat.EQ.0.OR.extrap) THEN
   icoord = i; jcoord = j
   IF (frel) THEN
      xcoord = x; ycoord = y
   ENDIF
   IF (PRESENT(weights)) THEN
      weights(1,1) = (1.0-x)*(1.0-y)
      weights(2,1) = x*(1.0-y)
      weights(1,2) = (1.0-x)*y
      weights(2,2) = x*y
   ENDIF
   IF (updatex) THEN
      xpos = xp; ypos = yp
   ENDIF
ENDIF


RETURN

END SUBROUTINE hrel_coords_rect_nonunif_rel_s

!========================================================================

SUBROUTINE hrel_coords_rect_nonunif_rel_d(xpos,ypos,displx,disply,xdat,ydat,&
                                        & nxdat,nydat,lbounds,ubounds,istat,&
                                        & icoord,jcoord,xcoord,ycoord,weights,&
                                        & extrapol,update)
!************************************************************************
!
! *hrel_coords_nonunif_rel_d* updates the horizontal relative coordinates and
!                             weight factors after applying a given displacement
!                             to the old location
!                     
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description - used in case of a non-uniform rectangular grid
!             - icoord, jcoord, xcoord, ycoord and weights are undefined
!               if istat=1 and extrapol is .FALSE.
!             - xcoord and ycoord must be both  or not present  
!             - xpos and ypos given in double precision  
!
!************************************************************************
!
USE syspars
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: extrapol, update
INTEGER, INTENT(IN) :: nxdat, nydat
INTEGER, INTENT(INOUT) :: icoord, jcoord
INTEGER, INTENT(OUT) :: istat
INTEGER, INTENT(IN), DIMENSION(2) :: lbounds, ubounds  
REAL, INTENT(IN) :: displx, disply
REAL (KIND=kndrlong), INTENT(INOUT) :: xpos, ypos
REAL, INTENT(IN), DIMENSION(nxdat) :: xdat
REAL, INTENT(IN), DIMENSION(nydat) :: ydat
REAL, INTENT(OUT), OPTIONAL :: xcoord, ycoord
REAL, INTENT(OUT), OPTIONAL, DIMENSION(2,2) :: weights

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*xpos*     REAL     X-coordinate of the input location [m or degrees longitude]
!*ypos*     REAL     Y-coordinate of the input location  [m or degrees latitude]
!*displx*   REAL     Applied displacement in the X-direction
!                                                       [m or degrees longitude]
!*disply*   REAL     Applied displacement in the Y-direction
!                                                        [m or degrees latitude]
!*xdat*     REAL     X-coordinates of the input data grid
!                                                       [m or degrees longitude]
!*ydat*     REAL     Y-coordinate of the input data grid [m or degrees latitude]
!*nxdat*    INTEGER  X-dimension of the input data grid
!*nydat*    INTEGER  Y-dimension of the input data grid
!*lbounds*  INTEGER  Lower bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*ubounds*  INTEGER  Upper bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*istat*    INTEGER  0/1 if location is inside/outside the input grid
!*icoord*   INTEGER  Returned X-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*jcoord*   INTEGER  Returned Y-index coordinate of the location with respect
!                    to the  lower left corner of the containing cell in the
!                    input grid
!*xcoord*   REAL     Returned normalised X-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*ycoord*     REAL     Returned normalised Y-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*weights*  REAL     Returned weight factors for interpolation
!*extrapol* LOGICAL  Allows extrapolation to the nearest inside grid point
!                    for data points outside the domain if PRESENT AND .TRUE.
!                    Default is .FALSE.
!*update*   LOGICAL  If PRESENT and .TRUE., the displaced locations are returned
!                    in xpos and ypos
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: extrap, frel, search, updatex
INTEGER :: i, j
REAL :: x, y
REAL (KIND=kndrlong) :: xp, yp


!
!1. Optional arguments
!---------------------
!

IF (PRESENT(extrapol)) THEN
   extrap = extrapol
ELSE
   extrap = .FALSE.
ENDIF

IF (PRESENT(update)) THEN
   updatex = update
ELSE
   updatex = .FALSE.
ENDIF

frel = PRESENT(xcoord).AND.PRESENT(ycoord)

!
!2. Search for containing grid cell
!----------------------------------
!
!2.1 X-direction
!---------------
!
!---eastward displacement
xp = xpos + displx
IF (displx.GE.0.0) THEN
   IF (xp.GT.xdat(ubounds(1))) THEN
      search = .FALSE.; i = ubounds(1)
   ELSE
      search = .TRUE.; i = icoord
   ENDIF
   DO WHILE (search.AND.i.LT.ubounds(1))
      IF (xp.GE.xdat(i).AND.xp.LT.xdat(i+1)) THEN
         search = .FALSE.
         x = (xp-xdat(i))/(xdat(i+1)-xdat(i))
      ELSE
         i = i + 1
      ENDIF
   ENDDO
   IF (.NOT.search) x = MIN(x,1.0)

!---westward displacement
ELSE
   IF (xp.LT.xdat(lbounds(1))) THEN
      search = .FALSE.; i = lbounds(1) - 1
   ELSE
      search = .TRUE.; i = icoord
   ENDIF
   DO WHILE (search.AND.i.GE.lbounds(1))
      IF (xp.GE.xdat(i).AND.xp.LT.xdat(i+1)) THEN
         search = .FALSE.
         x = (xp-xdat(i))/(xdat(i+1)-xdat(i))
      ELSE
         i = i - 1
      ENDIF
   ENDDO
   IF (.NOT.search) x = MIN(x,1.0)
ENDIF

!
!2.2 Y-direction
!---------------
!
!---northtward displacement
yp = ypos + disply
IF (disply.GE.0.0) THEN
   IF (yp.GT.ydat(ubounds(2))) THEN
      search = .FALSE.; j = ubounds(2)
   ELSE
      search = .TRUE.; j = jcoord
   ENDIF
   DO WHILE (search.AND.j.LT.ubounds(2))
      IF (yp.GE.ydat(j).AND.yp.LT.ydat(j+1)) THEN
         search = .FALSE.
         y = (yp-ydat(j))/(ydat(j+1)-ydat(j))
      ELSE
         j = j + 1
      ENDIF
   ENDDO
   IF (.NOT.search) y = MIN(y,1.0)

!---southward displacement
ELSE
   IF (yp.LT.ydat(lbounds(2))) THEN
      search = .FALSE.; j = lbounds(2) - 1
   ELSE
      search = .TRUE.; j = jcoord
   ENDIF
   DO WHILE (search.AND.j.GE.lbounds(2))
      IF (yp.GE.ydat(j).AND.yp.LT.ydat(j+1)) THEN
         search = .FALSE.
         y = (yp-ydat(j))/(ydat(j+1)-ydat(j))
      ELSE
         j = j - 1
      ENDIF
   ENDDO
   IF (.NOT.search) y = MIN(y,1.0)
ENDIF

!
!3. Check out of domain condition
!--------------------------------
!

IF (i.LT.lbounds(1).OR.i.GT.ubounds(1).OR.&
  & j.LT.lbounds(2).OR.j.GT.ubounds(2)) THEN
   istat = 1
ELSE
   istat = 0
ENDIF

!
!4. Extrapolate if needed
!------------------------
!

IF (istat.EQ.1.AND.extrap) THEN
   IF (i.LT.lbounds(1)) THEN
      i = lbounds(1); x = 0.0
   ELSEIF (i.GE.ubounds(1)) THEN
      i = ubounds(1)-1; x = 1.0
   ENDIF
   IF (j.LT.lbounds(2)) THEN
      j = lbounds(2); y = 0.0
   ELSEIF (j.GE.ubounds(2)) THEN
      j = ubounds(2)-1; y = 1.0
   ENDIF
ENDIF

!
!5. Relative coordinates and/or weights
!--------------------------------------
!

IF (istat.EQ.0.OR.extrap) THEN
   icoord = i; jcoord = j
   IF (frel) THEN
      xcoord = x; ycoord = y
   ENDIF
   IF (PRESENT(weights)) THEN
      weights(1,1) = (1.0-x)*(1.0-y)
      weights(2,1) = x*(1.0-y)
      weights(1,2) = (1.0-x)*y
      weights(2,2) = x*y
   ENDIF
   IF (updatex) THEN
      xpos = xp; ypos = yp
   ENDIF
ENDIF


RETURN

END SUBROUTINE hrel_coords_rect_nonunif_rel_d

!========================================================================

SUBROUTINE hrel_coords_rect_unif_abs_s(xpos,ypos,x0dat,y0dat,delxdat,delydat,&
                                     & lbounds,ubounds,istat,icoord,jcoord,&
                                     & xcoord,ycoord,weights,extrapol)
!************************************************************************
!
! *hrel_coords_rect_unif_abs_s* calculates the horizontal relative coordinates
!                               and weight factors for interpolating data from a
!                               uniform rectangular grid to a given location
!                    
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description - used in case of a uniform rectangular grid
!             - icoord, jcoord, xcoord, ycoord and weights are undefined
!               if istat=1 and extrapol is .FALSE.
!             - xcoord and ycoord must be both  or not presen
!             - xpos and ypos given in single precision
!
!************************************************************************
!  
USE syspars
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: extrapol
INTEGER, INTENT(OUT) :: icoord, istat, jcoord
INTEGER, DIMENSION(2), INTENT(IN) :: lbounds, ubounds  
REAL, INTENT(IN) :: delxdat, delydat, x0dat, y0dat
REAL (KIND=kndreal), INTENT(IN) :: xpos, ypos
REAL, INTENT(OUT), OPTIONAL :: xcoord, ycoord
REAL, INTENT(OUT), OPTIONAL, DIMENSION(2,2) :: weights

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*xpos*     REAL     X-coordinate of the input location [m or degrees longitude]
!*ypos*     REAL     Y-coordinate of the input location  [m or degrees latitude]
!*x0dat*    REAL     X-coordinate of lower left corner of the input data grid
!                                                       [m or degrees longitude]
!*y0dat*    REAL     Y-coordinate of lower left corner of the input data grid
!                                                        [m or degrees latitude]
!*delxdat*  REAL     Uniform grid spacing in the X-direction
!                                                       [m or degrees longitude]
!*delydat*  REAL     Uniform grid spacing in the Y-direction
!                                                        [m or degrees latitude]
!*lbounds*  INTEGER  Lower bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*ubounds*  INTEGER  Upper bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*istat*    INTEGER  0/1 if location is inside/outside the input grid
!*icoord*   INTEGER  Returned X-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*jcoord*   INTEGER  Returned Y-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*xcoord*   REAL     Returned normalised X-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*ycoord*   REAL     Returned normalised Y-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*weights*  REAL     Returned weight factors for interpolation
!*extrapol* LOGICAL  Allows extrapolation to the nearest inside grid point
!                    for data points outside the domain if PRESENT AND .TRUE.
!                    Default is .FALSE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: frel, extrap
INTEGER :: i, j
REAL :: x, y


!
!1. Optional arguments
!---------------------
!

IF (PRESENT(extrapol)) THEN
   extrap = extrapol
ELSE
   extrap = .FALSE.
ENDIF

frel = PRESENT(xcoord).AND.PRESENT(ycoord)

!
!2. Containing grid cell
!-----------------------
!

x = (xpos-x0dat)/delxdat + 1.0
y = (ypos-y0dat)/delydat + 1.0
i = INT(x); j = INT(y)
x = MIN(x-i,1.0); y = MIN(y-j,1.0)
IF (i.EQ.ubounds(1)) THEN
   i = ubounds(1) - 1; x = 1.0
ENDIF
IF (j.EQ.ubounds(2)) THEN
   j = ubounds(2) - 1; y = 1.0
ENDIF

!
!3. Check out of domain condition
!--------------------------------
!

IF (i.LT.lbounds(1).OR.i.GT.ubounds(1).OR.&
  & j.LT.lbounds(2).OR.j.GT.ubounds(2)) THEN
   istat = 1
ELSE
   istat = 0
ENDIF

!
!4. Extrapolate if needed
!------------------------
!

IF (istat.EQ.1.AND.extrap) THEN
   IF (i.LT.lbounds(1)) THEN
      i = lbounds(1); x = 0.0
   ELSEIF (i.GT.ubounds(1)) THEN
      i = ubounds(1)-1; x = 1.0
   ENDIF
   IF (j.LT.lbounds(2)) THEN
      j = lbounds(2); y = 0.0
   ELSEIF (j.GT.ubounds(2)) THEN
      j = ubounds(2)-1; y = 1.0
   ENDIF
ENDIF

!
!5. Relative coordinates and/or weights
!--------------------------------------
!

IF (istat.EQ.0.OR.extrap) THEN
   icoord = i; jcoord = j
   IF (frel) THEN
      xcoord = x; ycoord = y
   ENDIF
   IF (PRESENT(weights)) THEN
      weights(1,1) = (1.0-x)*(1.0-y)
      weights(2,1) = x*(1.0-y)
      weights(1,2) = (1.0-x)*y
      weights(2,2) = x*y
   ENDIF
ENDIF

   
RETURN

END SUBROUTINE hrel_coords_rect_unif_abs_s

!========================================================================

SUBROUTINE hrel_coords_rect_unif_abs_d(xpos,ypos,x0dat,y0dat,delxdat,delydat,&
                                     & lbounds,ubounds,istat,icoord,jcoord,&
                                     & xcoord,ycoord,weights,extrapol)
!************************************************************************
!
! *hrel_coords_rect_unif_abs_d* calculates the horizontal relative coordinates
!                               and weight factors for interpolating data from a
!                               uniform rectangular grid to a given location
!                    
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description - used in case of a uniform rectangular grid
!             - icoord, jcoord, xcoord, ycoord and weights are undefined
!               if istat=1 and extrapol is .FALSE.
!             - xcoord and ycoord must be both  or not presen
!             - xpos and ypos given in double precision
!
!************************************************************************
!
USE syspars  
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: extrapol
INTEGER, INTENT(OUT) :: icoord, istat, jcoord
INTEGER, DIMENSION(2), INTENT(IN) :: lbounds, ubounds  
REAL, INTENT(IN) :: delxdat, delydat, x0dat, y0dat
REAL (KIND=kndrlong), INTENT(IN) :: xpos, ypos
REAL, INTENT(OUT), OPTIONAL :: xcoord, ycoord
REAL, INTENT(OUT), OPTIONAL, DIMENSION(2,2) :: weights

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*xpos*     REAL     X-coordinate of the input location [m or degrees longitude]
!*ypos*     REAL     Y-coordinate of the input location  [m or degrees latitude]
!*x0dat*    REAL     X-coordinate of lower left corner of the input data grid
!                                                       [m or degrees longitude]
!*y0dat*    REAL     Y-coordinate of lower left corner of the input data grid
!                                                        [m or degrees latitude]
!*delxdat*  REAL     Uniform grid spacing in the X-direction
!                                                       [m or degrees longitude]
!*delydat*  REAL     Uniform grid spacing in the Y-direction
!                                                        [m or degrees latitude]
!*lbounds*  INTEGER  Lower bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*ubounds*  INTEGER  Upper bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*istat*    INTEGER  0/1 if location is inside/outside the input grid
!*icoord*   INTEGER  Returned X-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*jcoord*   INTEGER  Returned Y-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*xcoord*     REAL     Returned normalised X-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*ycoord*     REAL     Returned normalised Y-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*weights*  REAL     Returned weight factors for interpolation
!*extrapol* LOGICAL  Allows extrapolation to the nearest inside grid point
!                    for data points outside the domain if PRESENT AND .TRUE.
!                    Default is .FALSE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: frel, extrap
INTEGER :: i, j
REAL :: x, y


!
!1. Optional arguments
!---------------------
!

IF (PRESENT(extrapol)) THEN
   extrap = extrapol
ELSE
   extrap = .FALSE.
ENDIF

frel = PRESENT(xcoord).AND.PRESENT(ycoord)

!
!2. Containing grid cell
!-----------------------
!

x = (xpos-x0dat)/delxdat + 1.0
y = (ypos-y0dat)/delydat + 1.0
i = INT(x); j = INT(y)
x = MIN(x-i,1.0); y = MIN(y-j,1.0)
IF (i.EQ.ubounds(1)) THEN
   i = ubounds(1) - 1; x = 1.0
ENDIF
IF (j.EQ.ubounds(2)) THEN
   j = ubounds(2) - 1; y = 1.0
ENDIF

!
!3. Check out of domain condition
!--------------------------------
!

IF (i.LT.lbounds(1).OR.i.GT.ubounds(1).OR.&
  & j.LT.lbounds(1).OR.j.GT.ubounds(2)) THEN
   istat = 1
ELSE
   istat = 0
ENDIF

!
!4. Extrapolate if needed
!------------------------
!

IF (istat.EQ.1.AND.extrap) THEN
   IF (i.LT.lbounds(1)) THEN
      i = lbounds(1); x = 0.0
   ELSEIF (i.GT.ubounds(1)) THEN
      i = ubounds(1)-1; x = 1.0
   ENDIF
   IF (j.LT.lbounds(2)) THEN
      j = lbounds(2); y = 0.0
   ELSEIF (j.GT.ubounds(2)) THEN
      j = ubounds(2)-1; y = 1.0
   ENDIF
ENDIF

!
!5. Relative coordinates and/or weights
!--------------------------------------
!

IF (istat.EQ.0.OR.extrap) THEN
   icoord = i; jcoord = j
   IF (frel) THEN
      xcoord = x; ycoord = y
   ENDIF
   IF (PRESENT(weights)) THEN
      weights(1,1) = (1.0-x)*(1.0-y)
      weights(2,1) = x*(1.0-y)
      weights(1,2) = (1.0-x)*y
      weights(2,2) = x*y
   ENDIF
ENDIF

   
RETURN

END SUBROUTINE hrel_coords_rect_unif_abs_d

!========================================================================

SUBROUTINE hrel_coords_rect_unif_rel_s(xpos,ypos,displx,disply,delxdat,delydat,&
                                     & lbounds,ubounds,istat,icoord,jcoord,&
                                     & xcoord,ycoord,weights,extrapol,update)
!************************************************************************
!
! *hrel_coords_rect_unif_rel_s* updates the horizontal relative coordinates and
!                               weight factors after applying a given
!                               displacement to the old location
!                    
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description - used in case of a uniform rectangular grid
!             - icoord, jcoord, xcoord, ycoord and weights are undefined
!               if istat=1 and extrapol is .FALSE.
!             - xcoord and ycoord must be both  or not present  
!             - xpos and ypos given in single precision
!  
!************************************************************************
!
USE syspars  
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: extrapol, update
INTEGER, INTENT(INOUT) :: icoord, jcoord
INTEGER, INTENT(OUT) :: istat
INTEGER, DIMENSION(2), INTENT(IN) :: lbounds, ubounds  
REAL, INTENT(IN) :: delxdat, delydat, displx, disply
REAL (KIND=kndreal), INTENT(INOUT) :: xpos, ypos
REAL, INTENT(INOUT) :: xcoord, ycoord
REAL, INTENT(OUT), OPTIONAL, DIMENSION(2,2) :: weights

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*xpos*     REAL     X-coordinate of the input location [m or degrees longitude]
!*ypos*     REAL     Y-coordinate of the input location  [m or degrees latitude]
!*displx*   REAL     Applied displacement in the X-direction
!                                                       [m or degrees longitude]
!*disply*   REAL     Applied displacement in the Y-direction
!*delxdat*  REAL     Uniform grid spacing in the X-direction
!                                                       [m or degrees longitude]
!*delydat*  REAL     Uniform grid spacing in the Y-direction
!                                                        [m or degrees latitude]
!*lbounds*  INTEGER  Lower bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*ubounds*  INTEGER  Upper bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*istat*    INTEGER  0/1 if location is inside/outside the input grid
!*icoord*   INTEGER  Returned X-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*jcoord*   INTEGER  Returned Y-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*xcoord*   REAL     Returned normalised X-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*ycoord*   REAL     Returned normalised Y-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*weights*  REAL     Returned weight factors for interpolation
!*extrapol* LOGICAL  Allows extrapolation to the nearest inside grid point
!                    for data points outside the domain if PRESENT AND .TRUE.
!                    Default is .FALSE.
!*update*   LOGICAL  If PRESENT and .TRUE., the displaced locations are returned
!                    in xpos and ypos
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: extrap, updatex
INTEGER :: i, j
REAL (KIND=kndreal) :: x, y


!
!1. Optional argument
!--------------------
!

IF (PRESENT(extrapol)) THEN
   extrap = extrapol
ELSE
   extrap = .FALSE.
ENDIF

IF (PRESENT(update)) THEN
   updatex = update
ELSE
   updatex = .FALSE.
ENDIF

!
!2. Containing grid cell
!-----------------------
!

x = xcoord + displx/delxdat
y = ycoord + disply/delydat
i = MERGE(icoord+INT(x),icoord+INT(x)-1,x.GE.0.0)
j = MERGE(jcoord+INT(y),jcoord-INT(y)-1,y.GE.0.0)
x = MERGE(x-INT(x),1+x-INT(x),x.GE.0.0)
y = MERGE(y-INT(y),1+y-INT(y),y.GE.0.0)
x = MIN(x,1.0); y = MIN(y,1.0)

!
!3. Check out of domain condition
!--------------------------------
!

IF (i.LT.lbounds(1).OR.i.GE.ubounds(1).OR.&
  & j.LT.lbounds(1).OR.j.GE.ubounds(2)) THEN
   istat = 1
ELSE
   istat = 0
ENDIF

!
!4. Extrapolate if needed
!------------------------
!

IF (istat.EQ.1.AND.extrap) THEN
   IF (i.LT.lbounds(1)) THEN
      i = lbounds(1); x = 0.0
   ELSEIF (icoord.GT.ubounds(1)) THEN
      i = ubounds(1)-1; x = 1.0
   ENDIF
   IF (j.LT.lbounds(2)) THEN
      j = lbounds(2); y = 0.0
   ELSEIF (jcoord.GT.ubounds(2)) THEN
      j = ubounds(2)-1; y = 1.0
   ENDIF
ENDIF

!
!5. Relative coordinates and/or weights
!--------------------------------------
!

IF (istat.EQ.0.OR.extrap) THEN
   icoord = i; jcoord = j
   xcoord = x; ycoord = y
   IF (PRESENT(weights)) THEN
      weights(1,1) = (1.0-x)*(1.0-y)
      weights(2,1) = x*(1.0-y)
      weights(1,2) = (1.0-x)*y
      weights(2,2) = x*y
   ENDIF
   IF (updatex) THEN
      xpos = xpos + displx; ypos = ypos + disply
   ENDIF
ENDIF

   
RETURN

END SUBROUTINE hrel_coords_rect_unif_rel_s

!========================================================================

SUBROUTINE hrel_coords_rect_unif_rel_d(xpos,ypos,displx,disply,delxdat,delydat,&
                                     & lbounds,ubounds,istat,icoord,jcoord,&
                                     & xcoord,ycoord,weights,extrapol,update)
!************************************************************************
!
! *hrel_coords_rect_unif_rel_d* updates the horizontal relative coordinates and
!                               weight factors after applying a given
!                               displacement to the old location
!                    
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description - used in case of a uniform rectangular grid
!             - icoord, jcoord, xcoord, ycoord and weights are undefined
!               if istat=1 and extrapol is .FALSE.
!             - xcoord and ycoord must be both  or not present  
!             - xpos and ypos given in double precision
!  
!************************************************************************
!
USE syspars  
!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: extrapol, update
INTEGER, INTENT(INOUT) :: icoord, jcoord
INTEGER, INTENT(OUT) :: istat
INTEGER, DIMENSION(2), INTENT(IN) :: lbounds, ubounds  
REAL, INTENT(IN) :: delxdat, delydat, displx, disply
REAL (KIND=kndrlong), INTENT(INOUT) :: xpos, ypos
REAL, INTENT(INOUT) :: xcoord, ycoord
REAL, INTENT(OUT), OPTIONAL, DIMENSION(2,2) :: weights

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*xpos*     REAL     X-coordinate of the input location [m or degrees longitude]
!*ypos*     REAL     Y-coordinate of the input location  [m or degrees latitude]
!*displx*   REAL     Applied displacement in the X-direction
!                                                       [m or degrees longitude]
!*disply*   REAL     Applied displacement in the Y-direction
!*delxdat*  REAL     Uniform grid spacing in the X-direction
!                                                       [m or degrees longitude]
!*delydat*  REAL     Uniform grid spacing in the Y-direction
!                                                        [m or degrees latitude]
!*lbounds*  INTEGER  Lower bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*ubounds*  INTEGER  Upper bounds of the bounding box delineating the inside
!                    domain of the data grid (W/E/S/N)
!*istat*    INTEGER  0/1 if location is inside/outside the input grid
!*icoord*   INTEGER  Returned X-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*jcoord*   INTEGER  Returned Y-index coordinate of the location with respect
!                    to the lower left corner of the containing cell in the
!                    input grid
!*xcoord*   REAL     Returned normalised X-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*ycoord*   REAL     Returned normalised Y-coordinate (between 0 and 1) with
!                    respect to the origin of the containing grid cell at the
!                    input location
!*weights*  REAL     Returned weight factors for interpolation
!*extrapol* LOGICAL  Allows extrapolation to the nearest inside grid point
!                    for data points outside the domain if PRESENT AND .TRUE.
!                    Default is .FALSE.
!*update*   LOGICAL  If PRESENT and .TRUE., the displaced locations are returned
!                    in xpos and ypos
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: extrap, updatex
INTEGER :: i, j
REAL (KIND=kndrlong) :: x, y


!
!1. Optional argument
!--------------------
!

IF (PRESENT(extrapol)) THEN
   extrap = extrapol
ELSE
   extrap = .FALSE.
ENDIF

IF (PRESENT(update)) THEN
   updatex = update
ELSE
   updatex = .FALSE.
ENDIF

!
!2. Containing grid cell
!-----------------------
!

x = xcoord + displx/delxdat
y = ycoord + disply/delydat
i = MERGE(icoord+INT(x),icoord+INT(x)-1,x.GE.0.0)
j = MERGE(jcoord+INT(y),jcoord-INT(y)-1,y.GE.0.0)
x = MERGE(x-INT(x),1+x-INT(x),x.GE.0.0)
y = MERGE(y-INT(y),1+y-INT(y),y.GE.0.0)
x = MIN(x,1.0); y = MIN(y,1.0)

!
!3. Check out of domain condition
!--------------------------------
!

IF (i.LT.lbounds(1).OR.i.GE.ubounds(1).OR.&
  & j.LT.lbounds(1).OR.j.GE.ubounds(2)) THEN
   istat = 1
ELSE
   istat = 0
ENDIF

!
!4. Extrapolate if needed
!------------------------
!

IF (istat.EQ.1.AND.extrap) THEN
   IF (i.LT.lbounds(1)) THEN
      i = lbounds(1); x = 0.0
   ELSEIF (icoord.GT.ubounds(1)) THEN
      i = ubounds(1)-1; x = 1.0
   ENDIF
   IF (j.LT.lbounds(2)) THEN
      j = lbounds(2); y = 0.0
   ELSEIF (jcoord.GT.ubounds(2)) THEN
      j = ubounds(2)-1; y = 1.0
   ENDIF
ENDIF

!
!5. Relative coordinates and/or weights
!--------------------------------------
!

IF (istat.EQ.0.OR.extrap) THEN
   icoord = i; jcoord = j
   xcoord = x; ycoord = y
   IF (PRESENT(weights)) THEN
      weights(1,1) = (1.0-x)*(1.0-y)
      weights(2,1) = x*(1.0-y)
      weights(1,2) = (1.0-x)*y
      weights(2,2) = x*y
   ENDIF
   IF (updatex) THEN
      xpos = xpos + displx; ypos = ypos + disply
   ENDIF
ENDIF

   
RETURN

END SUBROUTINE hrel_coords_rect_unif_rel_d

!========================================================================

SUBROUTINE intpol_data_to_model_2d(invals,outvals,ncin,nrin,nosize,surfgrid,&
                                 & datflag,land,inflag,outflag)
!************************************************************************
!
! *intpol_data_to_model_2d* Interpolate a 2-D horizontal data array on the
!                           model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description - data are interpolated at model C-nodes
!
! Module calls - error_alloc
!
!************************************************************************
!
USE datatypes
USE grid
USE gridpars
USE syspars
USE error_routines, ONLY: error_alloc

!
!*Arguments
!
LOGICAL, INTENT(IN) :: land
INTEGER, INTENT(IN) :: datflag, ncin, nosize, nrin
REAL, INTENT(IN), OPTIONAL :: inflag, outflag
REAL, INTENT(IN), DIMENSION(ncin,nrin,nosize) :: invals
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,nosize) :: outvals
TYPE (HRelativeCoords), INTENT(IN), DIMENSION(ncloc,nrloc) :: surfgrid

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*invals*   REAL     Non-interpolated input data array
!*outvals*  REAL     Interpolated output data array
!*ncin*     INTEGER  X-dimension of input data grid
!*nrin*     INTEGER  Y-dimension of input data grid
!*nosize*   INTEGER  Number of output data parameters
!*surfgrid* DERIVED  Relative coordinates of model grid with respect to
!                    surface data grid
!*datflag*  INTEGER  Enables flagging in case of invalid data points
!                    = 0 : no flagging
!                    = 1 : interpolation performed if all interpolating data
!                          values are not flagged
!                    = 2 : interpolation performed if at least one
!                          interpolating data value is not flagged
!*land*     LOGICAL  Includes land points if .TRUE.
!*inflag*   REAL     Flag for invalid input data (must be present if datflag=1)
!*outflag*  REAL     Replacement flag at invalid data locations
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ii, isize, j, jj
REAL :: flagin, flagin_eps, flagout, wsum
REAL, DIMENSION(2,2) :: weights
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskvals


IF (SIZE(invals).EQ.0) RETURN

procname(pglev+1) = 'intpol_data_to_model_2d'
CALL log_timer_in()

!
!1. Allocate array
!-----------------
!

ALLOCATE (maskvals(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskvals',2,(/ncloc,nrloc/),kndlog)

!
!2. Output flags
!---------------
!

IF (PRESENT(outflag)) THEN
   flagout = outflag
ELSE
   flagout = 0.0
ENDIF
IF (PRESENT(inflag)) THEN
   flagin = inflag
ELSE
   flagin = 0.0
ENDIF
flagin_eps = 0.00001*ABS(flagin)

!
!3. Include/exclude land points
!------------------------------
!

IF (land) THEN
   maskvals = .TRUE.
ELSE
   maskvals = seapoint(1:ncloc,1:nrloc)
ENDIF

!
!4. Interpolate
!--------------
!
!4.1 No data flags
!-----------------
!

IF (datflag.EQ.0) THEN

   j_410: DO j=1,nrloc
   i_410: DO i=1,ncloc
      IF (maskvals(i,j)) THEN
         ii = surfgrid(i,j)%icoord; jj = surfgrid(i,j)%jcoord
         weights = surfgrid(i,j)%weights
         isize_411: DO isize=1,nosize
            outvals(i,j,isize) = SUM(weights*invals(ii:ii+1,jj:jj+1,isize))
         ENDDO isize_411
      ELSE
         outvals(i,j,:) = flagout
      ENDIF
   ENDDO i_410
   ENDDO j_410

!
!4.2 Only when all data are not flagged
!--------------------------------------
!

ELSEIF (datflag.EQ.1) THEN

   j_420: DO j=1,nrloc
   i_420: DO i=1,ncloc
      IF (maskvals(i,j)) THEN
         ii = surfgrid(i,j)%icoord; jj = surfgrid(i,j)%jcoord
         weights = surfgrid(i,j)%weights
         isize_421: DO isize=1,nosize
            IF (ALL(ABS((invals(ii:ii+1,jj:jj+1,isize)-flagin)).GT.&
                       & flagin_eps)) THEN
               outvals(i,j,isize) = SUM(weights*invals(ii:ii+1,jj:jj+1,isize))
            ELSE
               outvals(i,j,isize) = flagout
            ENDIF
         ENDDO isize_421
      ELSE
         outvals(i,j,:) = flagout
      ENDIF
   ENDDO i_420
   ENDDO j_420

!
!4.3 Only when at least one data value is not flagged
!----------------------------------------------------
!

ELSEIF (datflag.EQ.2) THEN

   j_430: DO j=1,nrloc
   i_430: DO i=1,ncloc
      IF (maskvals(i,j)) THEN
         ii = surfgrid(i,j)%icoord; jj = surfgrid(i,j)%jcoord
         weights = surfgrid(i,j)%weights
         isize_431: DO isize=1,nosize
            weights = MERGE(weights,0.0,&
                         & ABS((invals(ii:ii+1,jj:jj+1,isize)-flagin)).GT.&
                         & flagin_eps)
            wsum = SUM(weights)
            IF (wsum.GT.0) THEN
               outvals(i,j,isize) = SUM(weights*invals(ii:ii+1,jj:jj+1,isize))
            ELSE
               outvals(i,j,isize) = flagout
            ENDIF
         ENDDO isize_431
      ELSE
         outvals(i,j,:) = flagout
      ENDIF
   ENDDO i_430
   ENDDO j_430

ENDIF

!
!5. Deallocate arrays
!--------------------
!

DEALLOCATE (maskvals)

CALL log_timer_out()


RETURN

END SUBROUTINE intpol_data_to_model_2d

!========================================================================

SUBROUTINE intpol1d_model_to_dep(invals,outval,i,j,dep,cnode,nzmin,nzmax,&
                               & outflag)
!************************************************************************
!
! *intpol1d_model_to_dep* Interpolate a (1-D) vertical model array at a
!                         given depth
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.2
!
! Description -
!
! Module calls - error_alloc
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE syspars

!
!*Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: i, j
INTEGER, INTENT(IN), OPTIONAL :: nzmax, nzmin
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN) :: dep
REAL, INTENT(IN), DIMENSION(:) :: invals
REAL, INTENT(OUT) :: outval

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*invals*   REAL    Vertical model input array (of size nz or nz+1)
!*outval*   REAL    Returned interpolated value
!*i*        INTEGER (Local) X-index of the data location
!*j*        INTEGER (Local) Y-index of the data location
!*dep*      REAL    Interpolation depth (positive distance from the free
!                   surface)                                                 [m]
!*cnode*    CHAR    Vertical grid node where invals is defined ('C' or 'W')
!*nzmin*    INTEGER Lowest physical vertical level in case invals is defined
!                   at W-nodes (1 otherwise)
!*nzmax*    INTEGER Highest physical vertical level in case invals is defined
!                   at W-nodes (nz otherwise)
!*outflag*  REAL    Replacement flag if dep is larger than the water depth
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k
REAL :: sigk
REAL, DIMENSION(nz+1) :: gsig


!
!1. Interpolating depth below sea bed
!------------------------------------
!

IF (dep.GT.deptotatc(i,j)) THEN
   IF (PRESENT(outflag)) THEN
      outval = outflag
   ELSE
      outval = 0.0
   ENDIF
ENDIF

!
!2. C-node 
!---------
!

IF (TRIM(cnode).EQ.'C') THEN
   sigk = 1.0 - dep/deptotatc(i,j)
   gsig(1:nz) = gscoordatc(i,j,:)
   IF (sigk.LE.gsig(1)) THEN
      outval = invals(1)
   ELSEIF (sigk.GE.gsig(nz)) THEN
      outval = invals(nz)
   ELSE
      k = 1
      DO WHILE(sigk.GE.gsig(k+1))
         k = k + 1
      ENDDO
      outval = ((gsig(k+1)-sigk)*invals(k)+(sigk-gsig(k))*invals(k+1))/&
              & (gsig(k+1)-gsig(k))
   ENDIF
   
!
!3. W-node
!---------
!

ELSEIF (TRIM(cnode).EQ.'C') THEN
   sigk = 1.0-dep/deptotatc(i,j)
   gsig(nzmin:nzmax) = gscoordatw(i,j,nzmin:nzmax)
   IF (sigk.LE.gsig(nzmin)) THEN
      outval = invals(nzmin)
   ELSEIF (sigk.GE.gsig(nzmax)) THEN
      outval = invals(nzmax)
   ELSE
      k = 1
      DO WHILE(sigk.GE.gsig(k+1))
         k = k + 1
      ENDDO
      outval = ((gsig(k+1)-sigk)*invals(k)+(sigk-gsig(k))*invals(k+1))/&
              & (gsig(k+1)-gsig(k))
   ENDIF

ENDIF

END SUBROUTINE intpol1d_model_to_dep

!========================================================================

SUBROUTINE intpol2d_model_to_data_2d(invals,outvals,hcoords,nhdims,nhdat,&
                                   & nosize,datvals,outflag)
!************************************************************************
!
! *intpol2d_model_to_data_2d* Interpolate a 2-D model array at a number of 
!                             horizontal data locations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.9
!
! Description -
!
! Module calls - error_alloc
!
!************************************************************************
!
USE datatypes
USE grid
USE gridpars
USE syspars
USE error_routines, ONLY: error_alloc

!
!*Arguments
!
INTEGER, INTENT(IN) :: nhdat, nosize
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(4) :: nhdims
REAL, INTENT(INOUT), OPTIONAL, DIMENSION(nhdat) :: datvals
REAL, INTENT(OUT), DIMENSION(nhdat,nosize) :: outvals
REAL, INTENT(IN), DIMENSION(1-nhdims(1):ncloc+nhdims(2),&
                          & 1-nhdims(3):nrloc+nhdims(4),nosize) :: invals
TYPE (HRelativeCoords), INTENT(IN), DIMENSION(nhdat) :: hcoords

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*invals*   REAL     Input model array
!*outvals*  REAL     Model array interpolated at data points
!*hcoords*  DERIVED  Relative horizontal grid coordinate structure of data
!                    points with respect to model grid
!*nhdims*   INTEGER  Halo sizes of model array (WESN directions)
!*nhdat*    INTEGER  Number of horizontal data points
!*nosize*   INTEGER  Number of variables for interpolation
!*datvals*  REAL     Optional array of data values for flagging
!*outflag*  REAL     Replacement flag if interpolation cannot be performed
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, isize, j, ldat
REAL :: flagval
REAL, DIMENSION(2,2) :: weights
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:) :: maskvals


IF (SIZE(outvals).EQ.0) RETURN

procname(pglev+1) = 'intpol2d_model_to_data_2d'
CALL log_timer_in()

!
!1. Allocate array
!-----------------
!

ALLOCATE (maskvals(nhdat),STAT=errstat)
CALL error_alloc('maskvals',1,(/nhdat/),kndlog)

!
!2. Output flag
!--------------
!

IF (PRESENT(outflag)) THEN
   flagval = outflag
ELSE
   flagval = float_fill
ENDIF

!
!3. Data flags
!-------------
!

IF (PRESENT(datvals)) THEN
   maskvals = ABS(datvals-float_fill).GT.float_fill_eps
   WHERE (.NOT.maskvals) datvals = flagval
ELSE
   maskvals = .TRUE.
ENDIF

!
!4.Interpolate
!-------------
!

ldat_410: DO ldat=1,nhdat
   IF (maskvals(ldat)) THEN
      i = hcoords(ldat)%icoord; j = hcoords(ldat)%jcoord
      weights = hcoords(ldat)%weights
      isize_411: DO isize=1,nosize
         outvals(ldat,isize) = SUM(weights*invals(i:i+1,j:j+1,isize))
      ENDDO isize_411
   ELSE
      outvals(ldat,:) = flagval
   ENDIF
ENDDO ldat_410

!
!5. Deallocate array
!-------------------
!

DEALLOCATE (maskvals)

CALL log_timer_out()


RETURN

END SUBROUTINE intpol2d_model_to_data_2d

!========================================================================

SUBROUTINE intpol3d_data_to_data_2d(invals,outvals,vcoords,nhdat,nzdatin,&
                                  & nzdatout,nzeff,nosize,outflag)
!************************************************************************
!
! *intpol3d_data_to_data_2d* Interpolate a 3-D input data array on a 3-D output
!                            data array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.9
!
! Description - data grids share the same horizontal locations
!             - vertical dimensions may be different
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
INTEGER, INTENT(IN) :: nhdat, nosize, nzdatin, nzdatout, nzeff
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(IN), DIMENSION(nhdat,nzdatin,nosize) :: invals
REAL, INTENT(OUT), DIMENSION(nhdat,nzdatout,nosize) :: outvals
TYPE (VRelativeCoords), INTENT(IN), DIMENSION(nhdat,nzdatout) :: vcoords

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*invals*   REAL     Input data array
!*outvals*  REAL     Input data array interpolated at output data points
!*vcoords*  DERIVED  Relative vertical grid coordinate structure of input data
!                    data grid with respect to output data grid
!*nhdat*    INTEGER  Number of horizontal data locations
!*nzdatin*  INTEGER  Vertical dimension of input data
!*nzdatout* INTEGER  Vertical dimension of output grid
!*nzeff*    INTEGER  Effective vertical dimension of output grid
!*nosize*   INTEGER  Number of variables for interpolation
!*outflag*  REAL     Replacement flag if interpolation cannot be performed
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, kdat, ldat
REAL :: flagval
REAL, DIMENSION(2) :: weights


IF ((SIZE(invals).EQ.0).OR.(nzeff.EQ.0)) RETURN

procname(pglev+1) = 'intpol3d_data_to_data_2d'
CALL log_timer_in()

!
!1. Output flag
!--------------
!

IF (PRESENT(outflag)) THEN
   flagval = outflag
ELSE
   flagval = float_fill
ENDIF

!
!2. Interpolate
!--------------
!

ldat_210: DO ldat=1,nhdat
kdat_210: DO kdat=1,nzeff
   k = vcoords(ldat,kdat)%kcoord
   weights = vcoords(ldat,kdat)%weights
   outvals(ldat,kdat,:) = (weights(1)*invals(ldat,k,:)+&
                         & weights(2)*invals(ldat,k+1,:))
ENDDO kdat_210
ENDDO ldat_210

CALL log_timer_out()


RETURN

END SUBROUTINE intpol3d_data_to_data_2d

!========================================================================

SUBROUTINE intpol3d_model_to_data_1d(invals,outvals,vcoords,nzdim,nzdat,&
                                   & nosize,datvals,outflag)
!************************************************************************
!
! *intpol3d_model_to_data_1d* Interpolate a vertical model array at a data
!                             vertical profile array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.9
!
! Description - model and data are defined along the same horizontal grid
!             - interpolation is performed at grid point (i,j)
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
INTEGER, INTENT(IN) :: nosize, nzdat, nzdim
REAL, INTENT(IN), OPTIONAL :: outflag
REAL, INTENT(INOUT), OPTIONAL, DIMENSION(nzdat) :: datvals
REAL, INTENT(IN), DIMENSION(nzdim,nosize) :: invals
REAL, INTENT(OUT), DIMENSION(nzdat,nosize) :: outvals
TYPE (VRelativeCoords), INTENT(IN), DIMENSION(nzdat) :: vcoords

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*invals*   REAL     Model array 
!*outvals*  REAL     Model array interpolated at vertical locations of data
!*vcoords*  DERIVED  Relative vertical grid coordinate structure of data
!                    points with respect to model grid
!*nzdim*    INTEGER  Vertical dimension of model array
!*nzdat*    INTEGER  Vertical dimension of data array
!*nosize*   INTEGER  Number of variables for interpolation
!*datvals*  REAL     Optional array of data values for flagging
!*outflag*  REAL     Optional replacement flag for flagged output
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, kdat
REAL :: flagval
LOGICAL, DIMENSION(nzdat) :: maskvals
REAL, DIMENSION(2) :: weights


IF (nzdat.EQ.0) RETURN

procname(pglev+1) = 'intpol3d_model_to_data_1d'
CALL log_timer_in()

!
!1. Output flag
!--------------
!

IF (PRESENT(outflag)) THEN
   flagval = outflag
ELSE
   flagval = float_fill
ENDIF

!
!2. Mask for invalid data
!------------------------
!

IF (PRESENT(datvals)) THEN
   maskvals = ABS(datvals-float_fill).GT.float_fill_eps
   WHERE (.NOT.maskvals) datvals = flagval
ELSE
   maskvals = .TRUE.
ENDIF

!
!3.Interpolate
!-------------
!

kdat_310: DO kdat=1,nzdat
   IF (maskvals(kdat)) THEN
      k = vcoords(kdat)%kcoord
      weights = vcoords(kdat)%weights
      outvals(kdat,:) = weights(1)*invals(k,:) + weights(2)*invals(k+1,:)
   ELSE
      outvals(kdat,:) = flagval
   ENDIF
ENDDO kdat_310

CALL log_timer_out()


RETURN

END SUBROUTINE intpol3d_model_to_data_1d

!========================================================================

SUBROUTINE intpol3d_model_to_data_2d(invals,outvals,vcoords,nhdims,nzdim,&
                                   & nzdat,nosize,datvals,outflag)
!************************************************************************
!
! *intpol3d_model_to_data_2d* Interpolate a 3-D model array at a number of
!                             3-D data locations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.9
!
! Description - model and data are defined on the same horizontal grid
!
! Module calls - error_alloc
!
!************************************************************************
!
USE datatypes
USE gridpars
USE syspars
USE error_routines, ONLY: error_alloc

!
!*Arguments
!
INTEGER, INTENT(IN) :: nosize, nzdat, nzdim
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(4) :: nhdims
REAL, INTENT(IN), DIMENSION(1-nhdims(1):ncloc+nhdims(2),&
                          & 1-nhdims(3):nrloc+nhdims(4),nzdim,nosize) :: invals
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,nzdat,nosize) :: outvals
REAL, INTENT(INOUT), OPTIONAL, DIMENSION(ncloc,nrloc,nzdat) :: datvals
TYPE (VRelativeCoords), INTENT(IN), DIMENSION(ncloc,nrloc,nzdat) :: vcoords

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*invals*   REAL     Input model array
!*outvals*  REAL     Model array interpolated at data points
!*vcoords*  DERIVED  Relative vertical coordinats of data points with respect
!                    to model grid
!*nhdims*   INTEGER  Halo sizes of model array (WESN directions)
!*nzdim*    INTEGER  Vertical dimension of model array
!*nzdat*    INTEGER  Number of data points in the vertical
!*nosize*   INTEGER  Number of variables for interpolation
!*datvals*  REAL     Optional array of data values for flagging
!*outflag*  REAL     Replacement flag if interpolation cannot be performed
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j, k, kdat
REAL :: flagval
REAL, DIMENSION(2) :: weights
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskvals


IF (nzdat.EQ.0) RETURN

procname(pglev+1) = 'intpol3d_model_to_data_2d'
CALL log_timer_in()

!
!1. Allocate array
!-----------------
!

ALLOCATE (maskvals(ncloc,nrloc,nzdat),STAT=errstat)
CALL error_alloc('maskvals',3,(/ncloc,nrloc,nzdat/),kndlog)

!
!2. Output flag
!--------------
!

IF (PRESENT(outflag)) THEN
   flagval = outflag
ELSE
   flagval = float_fill
ENDIF

!
!3. Mask array for invalid data points
!-------------------------------------
!

IF (PRESENT(datvals)) THEN
   maskvals = ABS(datvals-float_fill).GT.float_fill_eps
   WHERE (.NOT.maskvals) datvals = flagval
ELSE
   maskvals = .TRUE.
ENDIF

!
!4.Interpolate
!-------------
!

j_410: DO j=1,nrloc
i_410: DO i=1,ncloc
kdat_410: DO kdat=1,nzdat
   IF (maskvals(i,j,kdat)) THEN
      k = vcoords(i,j,kdat)%kcoord
      weights = vcoords(i,j,kdat)%weights
      outvals(i,j,kdat,:) = weights(1)*invals(i,j,k,:) + &
                          & weights(2)*invals(i,j,k+1,:)
   ELSE
      outvals(i,j,kdat,:) = flagval
   ENDIF
ENDDO kdat_410
ENDDO i_410
ENDDO j_410

!
!5. Deallocate
!-------------
!

DEALLOCATE (maskvals)

CALL log_timer_out()


RETURN

END SUBROUTINE intpol3d_model_to_data_2d

!========================================================================

SUBROUTINE intpol3d_model_to_data_3d(invals,outvals,hcoords,vcoords,nhdims,&
                                   & nzdim,nhdat,nzdat,nosize,datvals,outflag)
!************************************************************************
!
! *intpol3d_model_to_data_3d* Interpolate a 3-D model array at a number of
!                             3-D data locations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.9
!
! Description - model and data may be definecon different grids
!
! Module calls - error_alloc
!
!************************************************************************
!
USE datatypes
USE grid
USE gridpars
USE syspars
USE error_routines, ONLY: error_alloc

!
!*Arguments
!
INTEGER, INTENT(IN) :: nhdat, nosize, nzdat, nzdim
REAL, INTENT(IN), OPTIONAL :: outflag
INTEGER, INTENT(IN), DIMENSION(4) :: nhdims
REAL, INTENT(INOUT), OPTIONAL, DIMENSION(nhdat,nzdat) :: datvals
REAL, INTENT(OUT), DIMENSION(nhdat,nzdat,nosize) :: outvals
REAL, INTENT(IN), DIMENSION(1-nhdims(1):ncloc+nhdims(2),&
                          & 1-nhdims(3):nrloc+nhdims(4),nzdim,nosize) :: invals
TYPE (HRelativeCoords), INTENT(IN), DIMENSION(nhdat) :: hcoords
TYPE (VRelativeCoords), INTENT(IN), DIMENSION(2,2,nhdat,nzdat) :: vcoords

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*invals*   REAL     Input model array
!*outvals*  REAL     Model array interpolated at data points
!*hcoords*  DERIVED  Relative horizontal grid coordinate structure of data
!                    points with respect to model grid
!*vcoords*  DERIVED  Relative vertical grid coordinate structure of data
!                    points with respect to model grid
!*nhdims*   INTEGER  Halo sizes of model array (WESN directions)
!*nzdim*    INTEGER  Vertical dimension of model array
!*nhdat*    INTEGER  Number of horizontal data points
!*nzdat*    INTEGER  Number of data points in the vertical
!*nosize*   INTEGER  Number of variables for interpolation
!*datvals*  REAL     Optional array of data values for flagging
!*outflag*  REAL     Optional replacement flag for flagged output
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ii, isize, isub, j, jj, jsub, k, kdat, ldat
REAL :: flagval
REAL, DIMENSION(2) :: weightsz
REAL, DIMENSION(2,2) :: weights
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskvals
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: zval


IF (SIZE(outvals).EQ.0) RETURN

procname(pglev+1) = 'intpol3d_model_to_data_3d'
CALL log_timer_in()

!
!1. Allocate array
!-----------------
!

ALLOCATE (maskvals(nhdat,nzdat),STAT=errstat)
CALL error_alloc('maskvals',2,(/nhdat,nzdat/),kndlog)
ALLOCATE (zval(2,2,nosize),STAT=errstat)
CALL error_alloc('zval',3,(/2,2,nosize/),kndrtype)

!
!2. Output flag
!--------------
!

IF (PRESENT(outflag)) THEN
   flagval = outflag
ELSE
   flagval = float_fill
ENDIF

!
!3. Mask array for invalid data points
!-------------------------------------
!

IF (PRESENT(datvals)) THEN
   ldat_310: DO ldat=1,nhdat
   kdat_310: DO kdat=1,nzdat
      IF (ABS(datvals(ldat,kdat)-float_fill).GT.float_fill_eps) THEN
         maskvals(ldat,kdat) = .TRUE.
      ELSE
         maskvals(ldat,kdat) = .FALSE.
         datvals(ldat,kdat) = flagval
      ENDIF
   ENDDO kdat_310
   ENDDO ldat_310
ELSE
   maskvals = .TRUE.
ENDIF

!
!4.Interpolate
!-------------
!

outvals = flagval
ldat_410: DO ldat=1,nhdat
   i = hcoords(ldat)%icoord; j = hcoords(ldat)%jcoord
   weights = hcoords(ldat)%weights
   kdat_411: DO kdat=1,nzdat
      IF (maskvals(ldat,kdat)) THEN
         ii_4111: DO ii=i,i+1
         jj_4111: DO jj=j,j+1
            isub = ii-i+1; jsub = jj-j+1
            k = vcoords(isub,jsub,ldat,kdat)%kcoord
            weightsz = vcoords(isub,jsub,ldat,kdat)%weights
            IF (weights(isub,jsub).GT.0.0) THEN
               zval(isub,jsub,:) =  weightsz(1)*invals(ii,jj,k,:) + &
                                  & weightsz(2)*invals(ii,jj,k+1,:)
            ELSE
               zval(isub,jsub,:) = 0.0
            ENDIF
         ENDDO jj_4111
         ENDDO ii_4111
         isize_4112: DO isize=1,nosize
            outvals(ldat,kdat,isize) = SUM(weights*zval(:,:,isize))
         ENDDO isize_4112
      ELSE
         outvals(ldat,kdat,:) = flagval
      ENDIF
   ENDDO kdat_411
ENDDO ldat_410

!
!5. Deallocate arrays
!--------------------
!

DEALLOCATE (maskvals,zval)

CALL log_timer_out()


RETURN

END SUBROUTINE intpol3d_model_to_data_3d

!========================================================================

SUBROUTINE model_to_data_hcoords(surfgrid,surfcoords,cnode,nxdat,nydat,&
                               & xdat,ydat,maskdat)
!************************************************************************
!
! *model_to_data_hcoords* Relative coordinates and weight factors used to
!                         interpolate local model data on an external grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)grid_interp.f90  V2.11
!
! Description -
!
! Module calls - collect_log_expr, distance_minloc, error_abort, error_alloc,
!                hrel_coords_curv_abs, hrel_coords_nounif_abs,
!                hrel_coords_curv_abs, UVarr_at_C, UVarr_at_U, UVarr_at_V
!
!************************************************************************
!
USE datatypes
USE grid
USE gridpars
USE modids
USE switches
USE syspars
USE array_interp, ONLY: UVarr_at_C, UVarr_at_U, UVarr_at_V
USE error_routines, ONLY: error_abort, error_alloc
USE grid_routines, ONLY: distance_minloc

!
!*Arguments
!
CHARACTER (LEN=lennode), INTENT(IN) :: cnode
INTEGER, INTENT(IN) :: nxdat, nydat
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(nxdat,nydat) :: maskdat
REAL, INTENT(IN), DIMENSION(nxdat,nydat) :: xdat, ydat
TYPE (GridParams), INTENT(IN) :: surfgrid
TYPE (HRelativeCoords), INTENT(INOUT), DIMENSION(ncloc,nrloc) :: surfcoords

!
! Name        Type     Purpose
!------------------------------------------------------------------------------
!*surfgrid*   DERIVED Attributes of the external grid
!*surfcoords* DERIVED Relative coordinates and weight factors for interpolation
!                     of the local model grid on the external grid
!*cnode*      CHAR    Nodal type of the model grid ('C','U','V', 'X','Y')
!*nxdat*      INTEGER X-dimension of the external grid
!*nydat*      INTEGER Y-dimension of the external grid
!*xdat*       REAL    X-coordinates of the external grid
!*ydat*       REAL    Y-coordinates of the external grid
!*maskdat*    LOGICAL Mask for land points on the external grid
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: extrapol, flagdat, land
CHARACTER (LEN=12) :: ci, cj
CHARACTER (LEN=15) :: cx, cy
INTEGER :: i, icoord, ierr, istat, j, jcoord, nhtype
INTEGER, DIMENSION(2) :: lbounds, ubounds
REAL :: delxdat, delydat, wsum
INTEGER, DIMENSION(2) :: locmin
REAL, DIMENSION(2,2) :: weights
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskloc
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: nstatus
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: xcoord, ycoord


procname(pglev+1) = 'model_to_data_hcoords'
CALL log_timer_in()

!
!1. Initialise
!-------------
!
!1.1 Optional arguments
!----------------------
!

flagdat = PRESENT(maskdat)

!
!1.2 Initialise parameters
!-------------------------
!

nhtype = surfgrid%nhtype
delxdat = surfgrid%delxdat
delydat = surfgrid%delydat
lbounds = 1; ubounds = (/nxdat,nydat/)
extrapol = surfgrid%extrapol
land = surfgrid%land

!
!1.3 Allocate arrays
!-------------------
!

ALLOCATE (maskloc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskloc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (xcoord(ncloc,nrloc),STAT=errstat)
CALL error_alloc('xcoord',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (ycoord(ncloc,nrloc),STAT=errstat)
CALL error_alloc('ycoord',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (nstatus(ncloc,nrloc),STAT=errstat)
CALL error_alloc('nstatus',2,(/ncloc,nrloc/),kndint)

nstatus = 0
   
!
!1.3 Local mask and coordinate arrays
!------------------------------------
!

SELECT CASE(TRIM(cnode))
   CASE ('C')
      maskloc = maskatc_int
      CALL UVarr_at_C(gxcoord(1:ncloc+1,1:nrloc+1),xcoord,0,0,(/1,1,nz/),&
                  & (/ncloc+1,nrloc+1,nz/),1,iarr_gxcoord,.TRUE.)
      CALL UVarr_at_C(gycoord(1:ncloc+1,1:nrloc+1),ycoord,0,0,(/1,1,nz/),&
                   & (/ncloc+1,nrloc+1,nz/),1,iarr_gycoord,.TRUE.)
   CASE ('U','Y')
      maskloc = node2du(1:ncloc,1:nrloc).GT.1
      CALL UVarr_at_U(gxcoord(1:ncloc,1:nrloc+1),xcoord,0,0,(/1,1,nz/),&
                   & (/ncloc,nrloc+1,nz/),1,iarr_gxcoord,.TRUE.)
      CALL UVarr_at_U(gycoord(1:ncloc,1:nrloc+1),ycoord,0,0,(/1,1,nz/),&
                   & (/ncloc,nrloc+1,nz/),1,iarr_gycoord,.TRUE.)
   CASE ('V','X')
      maskloc = node2dv(1:ncloc,1:nrloc).GT.1
      CALL UVarr_at_V(gxcoord(1:ncloc+1,1:nrloc),xcoord,0,0,(/1,1,nz/),&
                   & (/ncloc+1,nrloc,nz/),1,iarr_gxcoord,.TRUE.)
      CALL UVarr_at_V(gycoord(1:ncloc+1,1:nrloc),ycoord,0,0,(/1,1,nz/),&
                   & (/ncloc+1,nrloc,nz/),1,iarr_gycoord,.TRUE.)
END SELECT

!
!2. Relative coordinates and weights
!-----------------------------------
!

j_210: DO j=1,nrloc
i_210: DO i=1,ncloc

!  ---land mask
   IF (.NOT.maskloc(i,j)) THEN
      nstatus(i,j) = -1
      CYCLE i_210
   ENDIF

!  ---relative coordinates and weights
   SELECT CASE (nhtype)
      CASE (1)
         CALL hrel_coords_rect_unif_abs(xcoord(i,j),ycoord(i,j),xdat(1,1),&
              & ydat(1,1),delxdat,delydat,lbounds,ubounds,istat,icoord,jcoord,&
              & weights=weights,extrapol=extrapol)
      CASE (2)
         CALL hrel_coords_rect_nonunif_abs(xcoord(i,j),ycoord(i,j),xdat(:,1),&
              & ydat(1,:),nxdat,nydat,lbounds,ubounds,istat,icoord,jcoord,&
              & weights=weights,extrapol=extrapol)
      CASE (3)
         CALL hrel_coords_curv_abs(xcoord(i,j),ycoord(i,j),xdat,ydat,nxdat,&
              & nydat,lbounds,ubounds,istat,1,icoord,jcoord,weights=weights,&
              & extrapol=extrapol)
   END SELECT
   
!  ---status of location (outside domain or land), store data
   IF (istat.EQ.1.AND.(.NOT.extrapol)) THEN
      nstatus(i,j) = 1
   ELSE
      surfcoords(i,j)%icoord = icoord
      surfcoords(i,j)%jcoord = jcoord
      IF (flagdat) THEN
         weights = MERGE(weights,0.0,maskdat(icoord:icoord+1,jcoord:jcoord+1))
      ENDIF
      wsum = SUM(weights)
      IF (wsum.EQ.0.0) THEN
         nstatus(i,j) = 2
      ELSE
         weights = weights/wsum
      ENDIF
      surfcoords(i,j)%weights = weights
   ENDIF
   
ENDDO i_210
ENDDO j_210

!
!3. Perform extrapolation in case of model grid point is located inside
!   a land cell on the data grid
!----------------------------------------------------------------------
!

IF (land) THEN
   j_310: DO j=1,nrloc
   i_310: DO i=1,ncloc
      IF (nstatus(i,j).EQ.2) THEN
         IF (flagdat) THEN
            CALL distance_minloc(xcoord(i,j),ycoord(i,j),xdat,ydat,nxdat,nydat,&
                            & 2,3,locmin=locmin,mask=maskdat)
         ELSE
            CALL distance_minloc(xcoord(i,j),ycoord(i,j),xdat,ydat,nxdat,nydat,&
                            & 2,3,locmin=locmin)
         ENDIF
         surfcoords(i,j)%icoord = locmin(1)
         surfcoords(i,j)%jcoord = locmin(2)
         surfcoords(i,j)%weights(1:2,1) = (/1.0,0.0/)
         surfcoords(i,j)%weights(1:2,2) = 0.0
         nstatus(i,j) = 0
      ENDIF
   ENDDO i_310
   ENDDO j_310
ENDIF

!
!4. Check locations
!------------------
!
   
nerrs = COUNT(nstatus.GT.0)

IF (nerrs.GT.0.AND.errchk) THEN

   ierr = 0
   i_410: DO i=1,ncloc
   j_410: DO j=1,nrloc
      IF (nstatus(i,j).GT.0.AND.ierr.LE.maxerrors) THEN
         WRITE (ci,'(I12)') i+nc1loc-1; ci = ADJUSTL(ci)
         WRITE (cj,'(I12)') j+nr1loc-1; cj = ADJUSTL(cj)
         WRITE (cx,'(G15.7)') xcoord(i,j); cx = ADJUSTL(cx)
         WRITE (cy,'(G15.7)') ycoord(i,j); cy = ADJUSTL(cy)
         IF (nstatus(i,j).EQ.1) THEN
            WRITE (ioerr,'(A)') 'Model grid point ('//TRIM(ci)//','//&
                       & TRIM(cj)//') at ('//TRIM(cx)//','//TRIM(cy)//&
                       & ') is outside the data grid'
         ELSEIF (nstatus(i,j).EQ.2) THEN
            WRITE (ioerr,'(A)') 'Model grid point ('//TRIM(ci)//','//&
                       & TRIM(cj)//') at ('//TRIM(cx)//','//TRIM(cy)//&
                       & ') is located on land'
         ENDIF
         ierr = ierr + 1
      ENDIF
   ENDDO j_410
   ENDDO i_410
ENDIF

CALL error_abort('model_to_data_hcoords',ierrno_inival)

!
!5. Deallocate arrays
!--------------------
!

DEALLOCATE (maskloc,nstatus,xcoord,ycoord)

CALL log_timer_out()


RETURN

END SUBROUTINE model_to_data_hcoords


END MODULE grid_interp
