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
! *Grid_Arrays* Ensemble of routines to define model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Grid_Arrays.f90  V2.11.2
!
! $Date: 2018-04-17 12:38:28 +0200 (Tue, 17 Apr 2018) $
!
! $Revision: 1121 $
!
! Description - 
!
! Reference -
!
! Subroutines - depth_at_nodes, define_global_grid, define_local_grid,
!               grid_arrays, open_boundary_arrays, pointer_arrays, read_grid,
!               store_depths_old, update_pointer_arrays, water_depths,
!               write_grid 
!
!************************************************************************
!

!========================================================================

SUBROUTINE depth_at_nodes
!************************************************************************
!
! *depth_at_nodes* Interpolate mean water depths at nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Grid_Arrays.f90  V2.5.1
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! Externals -
!
! Module calls - exchange_mod, Carr_at_U, Carr_at_V, Carr_at_UV
!
!************************************************************************
!
USE depths
USE gridpars
USE iopars
USE modids
USE switches
USE array_interp, ONLY: Carr_at_U, Carr_at_V, Carr_at_UV
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER, DIMENSION(2) :: lbounds
INTEGER, DIMENSION(4) :: nhexch


procname(pglev+1) = 'depth_at_nodes'
CALL log_timer_in()

!
!1. Interpolate
!-------------
!
!---U-nodes
CALL Carr_at_U(depmeanatc(0:ncloc,1:nrloc),depmeanatu(1:ncloc,1:nrloc),1,1,&
            & (/0,1,nz/),(/ncloc,nrloc,nz/),1,iarr_depmeanatu,.TRUE.)

!---V-nodes
CALL Carr_at_V(depmeanatc(1:ncloc,0:nrloc),depmeanatv(1:ncloc,1:nrloc),1,1,&
            & (/1,0,nz/),(/ncloc,nrloc,nz/),1,iarr_depmeanatv,.TRUE.)

!---UV-nodes
CALL Carr_at_UV(depmeanatc(0:ncloc+1,0:nrloc+1),&
            & depmeanatuv(1:ncloc+1,1:nrloc+1),1,1,(/0,0,nz/),&
            & (/ncloc+1,nrloc+1,nz/),1,iarr_depmeanatuv,.TRUE.)

!
!2. Exchange halo sections
!------------------------
!

IF (iopt_MPI.EQ.1) THEN
   lbounds = 0; nhexch = 1
   CALL exchange_mod(depmeanatu,lbounds,nhexch,iarr_depmeanatu)
   CALL exchange_mod(depmeanatv,lbounds,nhexch,iarr_depmeanatv)
ENDIF

CALL log_timer_out


RETURN

END SUBROUTINE depth_at_nodes

!========================================================================

SUBROUTINE define_global_grid
!************************************************************************
!
! *define_global_grid* Define global grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Grid_Arrays.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! Externals - dry_cells, read_dry_cells, usrdef_dry_cells, write_dry_cells
!
! Module calls - construct_rectgrid, error_abort, error_alloc,
!                grid_rotation_params, UVarr_at_C
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE array_interp, ONLY: UVarr_at_C
USE error_routines, ONLY: error_abort, error_alloc
USE grid_routines, ONLY: construct_rectgrid, grid_rotation_params
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
CHARACTER (LEN=12) :: clev, citer
INTEGER :: icount, i, j, k
REAL :: alpha, c1, c2, c3, df, f, signew, sigold, sigreg, xstart, ystart
REAL, PARAMETER :: epsconv = 1.0E-06, maxits = 10


procname(pglev+1) = 'define_global_grid'
CALL log_timer_in()

!
!1. Horizontal grid
!------------------
!
!1.1 Rectangular grid
!--------------------
!

IF (iopt_grid_htype.LE.2) THEN

!
!1.1.1 Parameters for grid rotation
!----------------------------------
!

   IF (surfacegrids(igrd_model,1)%rotated.AND.iopt_grid_sph.EQ.1) THEN
      CALL grid_rotation_params(igrd_model)
   ELSE
      surfacegrids(igrd_model,1)%x0rot = surfacegrids(igrd_model,1)%x0dat 
      surfacegrids(igrd_model,1)%y0rot = surfacegrids(igrd_model,1)%y0dat
      surfacegrids(igrd_model,1)%longpole = 0.0
   ENDIF

!
!1.1.2 Grid spacings
!-------------------
!

   IF (iopt_grid_htype.EQ.1) THEN
      gdelxglb = surfacegrids(igrd_model,1)%delxdat
      gdelyglb = surfacegrids(igrd_model,1)%delydat
   ELSE
      gdelxglb(0) = gdelxglb(1); gdelxglb(nc+1) = gdelxglb(nc)
      gdelyglb(0) = gdelyglb(1); gdelyglb(nr+1) = gdelyglb(nr)
   ENDIF

!
!1.1.3 Construct grid
!--------------------
!

   xstart = surfacegrids(igrd_model,1)%x0rot - gdelxglb(0)
   ystart = surfacegrids(igrd_model,1)%y0rot - gdelyglb(0)
   CALL construct_rectgrid(xstart,ystart,gdelxglb,gdelyglb,&
                         & gxcoordglb,gycoordglb,nc+2,nr+2)
ENDIF

!
!1.2 Extend curvilinear grid outside the global domain
!-----------------------------------------------------
!

IF (iopt_grid_htype.EQ.3) THEN

!
!1.2.1 X-coordinates
!-------------------
!

   gxcoordglb(0,1:nr) = gxcoordglb(1,1:nr) + &
                     & (gxcoordglb(1,1:nr)-gxcoordglb(2,1:nr))
   gxcoordglb(nc+1,1:nr) = gxcoordglb(nc,1:nr) + &
                        & (gxcoordglb(nc,1:nr)-gxcoordglb(nc-1,1:nr))
   gxcoordglb(1:nc,0) = gxcoordglb(1:nc,1) + &
                     & (gxcoordglb(1:nc,1)-gxcoordglb(1:nc,2))
   gxcoordglb(1:nc,nr+1) = gxcoordglb(1:nc,nr) + &
                        & (gxcoordglb(1:nc,nr)-gxcoordglb(1:nc,nr-1))
   gxcoordglb(0,0) = gxcoordglb(1,0) + (gxcoordglb(1,0)-gxcoordglb(2,0))
   gxcoordglb(0,nr+1) = gxcoordglb(1,nr+1) + &
                     & (gxcoordglb(1,nr+1)-gxcoordglb(2,nr+1))
   gxcoordglb(nc+1,0) = gxcoordglb(nc+1,1) + &
                     & (gxcoordglb(nc+1,1)-gxcoordglb(nc+1,2))
   gxcoordglb(nc+1,nr+1) = gxcoordglb(nc,nr+1) + &
                        & (gxcoordglb(nc,nr+1)-gxcoordglb(nc-1,nr+1))
      
!
!1.2.2 Y-coordinates
!-------------------
!

   gycoordglb(0,1:nr) = gycoordglb(1,1:nr) + &
                     & (gycoordglb(1,1:nr)-gycoordglb(2,1:nr))
   gycoordglb(nc+1,1:nr) = gycoordglb(nc,1:nr) + &
                        & (gycoordglb(nc,1:nr)-gycoordglb(nc-1,1:nr))
   gycoordglb(1:nc,0) = gycoordglb(1:nc,1) + &
                     & (gycoordglb(1:nc,1)-gycoordglb(1:nc,2))
   gycoordglb(1:nc,nr+1) = gycoordglb(1:nc,nr) + &
                        & (gycoordglb(1:nc,nr)-gycoordglb(1:nc,nr-1))
   gycoordglb(0,0) = gycoordglb(1,0) + (gycoordglb(1,0)-gycoordglb(2,0))
   gycoordglb(0,nr+1) = gycoordglb(1,nr+1) + &
                     & (gycoordglb(1,nr+1)-gycoordglb(2,nr+1))
   gycoordglb(nc+1,0) = gycoordglb(nc+1,1) + &
                     & (gycoordglb(nc+1,1)-gycoordglb(nc+1,2))
   gycoordglb(nc+1,nr+1) = gycoordglb(nc,nr+1) + &
                        & (gycoordglb(nc,nr+1)-gycoordglb(nc-1,nr+1))

ENDIF

!
!1.3 Dry cells
!-------------
!

IF (iopt_drycel.EQ.1) THEN
   IF (modfiles(io_drycel,1,1)%status.EQ.'R') THEN
      CALL read_dry_cells
   ELSEIF (modfiles(io_drycel,1,1)%status.EQ.'N') THEN
      CALL usrdef_dry_cells
   ENDIF
   IF (modfiles(io_drycel,1,2)%status.EQ.'W') CALL write_dry_cells
   CALL dry_cells
ENDIF

!
!2. Water depths
!---------------
!

depmeanglb(1-nhalo:0,:) = depmean_flag
depmeanglb(:,1-nhalo:0) = depmean_flag
depmeanglb(nc:nc+nhalo,:) = depmean_flag
depmeanglb(:,nr:nr+nhalo) = depmean_flag

!---global mask
seapointglb = MERGE(.TRUE.,.FALSE.,&
              & ABS(depmeanglb(0:nc+1,0:nr+1)-depmean_flag).GT.depmean_flag_eps)
      
!
!3. Vertical grid
!----------------
!
!3.1 Horizontally uniform
!------------------------
!

IF (iopt_grid_vtype.LE.2) THEN
   gsigcoordatw(1) = 0.0
   gsigcoordatw(nz+1) = 1.0
ENDIF

!
!3.1.1 Vertically uniform
!------------------------
!

IF (iopt_grid_vtype.EQ.1) THEN

   k_311: DO k=2,nz
      gsigcoordatw(k) = REAL(k-1)/REAL(nz)
   ENDDO k_311

!
!3.1.2 Vertically non-uniform
!----------------------------
!

ELSEIF (iopt_grid_vtype.EQ.2) THEN

!  ---log-transformation (bottom)
   IF (iopt_grid_vtype_transf.EQ.11) THEN
      alpha = LOG(1.0+1.0/sig0_DJ)+1.0/sigstar_DJ
      k_3121: DO k=2,nz
         sigreg = (k-1)/REAL(nz)
         icount = 0
         signew = gsigcoordatw(k-1)
         sigold = 10.0
         DO WHILE (ABS(signew-sigold).GT.epsconv)
            icount = icount + 1
            IF (icount.GT.maxits) GOTO 1000
            sigold = signew
            f = alpha*sigreg - LOG(1.0+sigold/sig0_DJ) - sigold/sigstar_DJ
            df = -1.0/(sigold+sig0_DJ)-1.0/sigstar_DJ
            signew = sigold-f/df
         ENDDO
         gsigcoordatw(k) = signew
      ENDDO k_3121

!  ---log-transformation (surface)
   ELSEIF (iopt_grid_vtype_transf.EQ.12) THEN
      alpha = LOG(1.0+1.0/sig0_DJ)+1.0/sigstar_DJ
      k_3122: DO k=nz,2,-1
         sigreg = (k-1)/REAL(nz)
         icount = 0
         signew = gsigcoordatw(k+1)
         sigold = 10.0
         DO WHILE (ABS(signew-sigold).GT.epsconv)
            icount = icount + 1
            IF (icount.GT.maxits) GOTO 1000
            sigold = signew
            f = alpha*(sigreg-1.0) + LOG((1.0-sigold)/sig0_DJ+1.0) + &
             & (1.0-sigold)/sigstar_DJ
            df = -1.0/(1.0-sigold+sig0_DJ)-1.0/sigstar_DJ
            signew = sigold-f/df
         ENDDO
         gsigcoordatw(k) = signew
      ENDDO k_3122

!  ---Burchard and Bolding (2002)
   ELSEIF (iopt_grid_vtype_transf.EQ.13) THEN
      c1 = dl_BB + du_BB; c2 = TANH(dl_BB); c3 = TANH(dl_BB)+TANH(du_BB) 
      k_3123: DO k=2,nz
         sigreg = (k-1)/REAL(nz)
         gsigcoordatw(k) = TANH((c1*sigreg-dl_BB)+c2)/c3
      ENDDO k_3123
   ENDIF

ENDIF

!
!3.1.3 Insert in global array
!----------------------------
!

IF (iopt_grid_vtype.LE.2) THEN
   k_313: DO k=1,nz+1
      WHERE (seapointglb)
         gscoordglb(:,:,k) = gsigcoordatw(k)
      ELSEWHERE
         gscoordglb(:,:,k) = 0.0
      END WHERE
   ENDDO k_313
ENDIF

!
!3.1.4 Sigma coordinates at centered nodes
!-----------------------------------------
!

IF (iopt_grid_vtype.LE.2) THEN
   gsigcoordatc = 0.5*(gsigcoordatw(1:nz)+gsigcoordatw(2:nz+1))
ENDIF

!
!3.2 Horizontally non-uniform
!----------------------------
!

IF (iopt_grid_vtype.EQ.3) THEN

!
!3.2.1 Initialise
!----------------
!

   k_321: DO k=1,nz+1
      WHERE (.NOT.seapointglb)
         gscoordglb(:,:,k) = 0.0
      END WHERE
   ENDDO k_321
   WHERE (seapointglb)
      gscoordglb(:,:,1) = 0.0
      gscoordglb(:,:,nz+1) = 1.0
   END WHERE

!
!3.2.2 Song and Haidvogel (1994)    
!-------------------------------
!

   IF (iopt_grid_vtype_transf.EQ.21) THEN

      gscoordglb(:,:,1) = 0.0
      k_322: DO k=2,nz
         sigreg = (k-1)/REAL(nz)
         c1 = (1.0-b_SH)*SINH(theta_SH*(sigreg-1.0))/SINH(theta_SH)
         c2 = 0.5*b_SH*(TANH(theta_SH*(sigreg-0.5))/TANH(0.5*theta_SH)-1.0)
         c3 = c1 + c2 + 1.0 - sigreg
         j_3221: DO j=1,nr-1
         i_3221: DO i=1,nc-1
            IF (seapointglb(i,j)) then
               f = MERGE(((depmeanglb(i,j)-hcrit_SH)*c3)/depmeanglb(i,j),0.0,&
                         & depmeanglb(i,j).GT.hcrit_SH)
               gscoordglb(i,j,k) = sigreg + f
            ELSE
               gscoordglb(i,j,k) = 0.0
            ENDIF
         ENDDO i_3221
         ENDDO j_3221
      ENDDO k_322

   ENDIF

ENDIF

!
!4. Open Boundary locations (1-D)
!--------------------------------
!

IF (iopt_grid_nodim.EQ.1) THEN
   iobu(1) = 1; jobu(1) = 1; iobu(2) = 1; jobu(2) = 2
   iobu(3) = 3; jobu(3) = 1; iobu(4) = 3; jobu(4) = 2
   iobv(1) = 1; jobv(1) = 1; iobv(2) = 2; jobv(2) = 1
   iobv(3) = 1; jobv(3) = 3; iobv(4) = 2; jobv(4) = 3
ENDIF

CALL log_timer_out()


RETURN

1000 CONTINUE
nerrs = 1
IF (errchk) THEN
   WRITE (clev,'(I12)') k; clev = ADJUSTL(clev)
   WRITE (citer,'(I12)') icount-1; citer = ADJUSTL(citer)
   WRITE (ioerr,'(A)') 'Algorithm for calculating value of gsigcoordatw at '// &
                     & ' level '//TRIM(clev)//' did not converge after '//&
                     & TRIM(citer)//' iterations'
ENDIF
CALL error_abort('define_global_grid',ierrno_inival)

END SUBROUTINE define_global_grid

!========================================================================

SUBROUTINE define_local_grid
!************************************************************************
!
! *define_local_grid* Define local grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Grid_Arrays.f90  V2.5.1
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! Externals -
!
! Module calls -
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'define_local_grid'
CALL log_timer_in()

gxcoord = gxcoordglb(nc1loc-1:nc2loc+1,nr1loc-1:nr2loc+1)
gycoord = gycoordglb(nc1loc-1:nc2loc+1,nr1loc-1:nr2loc+1)

depmeanatc = depmeanglb(nc1loc-nhalo:nc2loc+nhalo,nr1loc-nhalo:nr2loc+nhalo)
gscoordatw = gscoordglb(nc1loc-1:nc2loc+1,nr1loc-1:nr2loc+1,:)

CALL log_timer_out()


RETURN

END SUBROUTINE define_local_grid

!========================================================================

SUBROUTINE grid_arrays
!************************************************************************
!
! *grid_arrays* Define arrays related to model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Grid_Arrays.f90  V2.11.2
!
! Description - horizontal grid spacings, sigma-coordinates,
!               geographical coordinates, Coriolis frequency and acceleration
!               of gravity
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - complex_polar, distance_minloc, distribute_mod, error_alloc,
!                error_lbound_var, exchange_mod, grid_rotation_transf,
!                grid_spacings_curv, grid_spacings_rectang, max_vars, min_vars,
!                sum2_vars, Carr_at_UV, UVarr_at_C, UVarr_at_U, UVarr_at_V,
!                UWarr_at_U, VWarr_at_V, Warr_at_C, Warr_at_UW, Warr_at_VW
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE array_interp, ONLY: Carr_at_UV, UVarr_at_C, UVarr_at_U, UVarr_at_V, &
                      & UWarr_at_U, VWarr_at_V, Warr_at_C, Warr_at_UW, &
                      & Warr_at_VW
USE error_routines, ONLY: error_alloc, error_lbound_var
USE grid_routines, ONLY: distance_minloc, grid_rotation_transf, &
                       & grid_spacings_curv, grid_spacings_rectang
USE math_library, ONLY: complex_polar
USE paral_comms, ONLY: distribute_mod, exchange_mod
USE paral_utilities, ONLY: max_vars, min_vars, sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: i, ii, j, jj, k, l, nhtype
REAL :: c1, c2, c3, delxatc_mean, delyatc_mean, distmin, f, hgrid, rlatref, &
      & rotfreq, sigreg, x, y, y0dat
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhdims, nhdist, nhexch
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: xcoordobc, ycoordobc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, array2du, array2dv,&
                                         & delglb, delx, dely, xcoord, ycoord


procname(pglev+1) = 'grid_arrays'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2du(ncloc+1,0:nrloc),STAT=errstat)
CALL error_alloc('array2du',2,(/ncloc+1,nrloc+1/),kndrtype)
ALLOCATE (array2dv(0:ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('array2dv',2,(/ncloc+1,nrloc+1/),kndrtype)
ALLOCATE (delglb(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo),STAT=errstat)
CALL error_alloc('delglb',2,(/nc+2*nhalo,nr+2*nhalo/),kndrtype)
ALLOCATE (delx(ncloc,nrloc),STAT=errstat)
CALL error_alloc('delx',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (dely(ncloc,nrloc),STAT=errstat)
CALL error_alloc('dely',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (xcoord(0:nc+1,0:nr+1),STAT=errstat)
CALL error_alloc('xcoord',2,(/nc+2,nr+2/),kndrtype)
ALLOCATE (ycoord(0:nc+1,0:nr+1),STAT=errstat)
CALL error_alloc('ycoord',2,(/nc+2,nr+2/),kndrtype)

!
!2. Shortcuts
!------------
!

nhtype = iopt_grid_htype
y0dat = surfacegrids(igrd_model,1)%y0dat
dXregX = nhtype.EQ.1
dXregY = nhtype.LT.3.AND.iopt_grid_sph.EQ.0
dYregX = nhtype.LT.3
dYregY = nhtype.EQ.1
dXYreg = nhtype.EQ.1.AND.iopt_grid_sph.EQ.0
dZregZ = iopt_grid_vtype.EQ.1
rotate_gvecs = surfacegrids(igrd_model,1)%rotated.OR.iopt_grid_htype.EQ.3

!
!3. Horizontal grid spacings
!---------------------------
!

lbounds(1:2) = 1-nhalo; nhdist = nhalo; nhdims = 0

delglb = 0.0

!
!3.1 X-direction
!---------------
!
!---V-nodes
IF (iopt_grid_htype.LE.2) THEN
   CALL grid_spacings_rectang(delglb(0:nc,0:nr+1),gdelxglb(0:nc),&
                            & gycoordglb(1,:),2,3,'X')
ELSE
   CALL grid_spacings_curv(delglb(0:nc,0:nr+1),gxcoordglb(0:nc,:),&
                         & gxcoordglb(1:nc+1,:),gycoordglb(0:nc,:),&
                         & gycoordglb(1:nc+1,:),2,3)
ENDIF
delglb(nc+1,:) = delglb(nc,:)
CALL distribute_mod(delglb,delxatv,lbounds(1:2),nhdist,iarr_delxatv,0.0,&
                  & mglevel=0)

!---C-nodes
IF (iopt_grid_htype.LE.2) THEN
   ycoord(1,0:nr) = 0.5*(gycoordglb(1,0:nr)+gycoordglb(1,1:nr+1))
   CALL grid_spacings_rectang(delglb(0:nc,0:nr),gdelxglb(0:nc),&
                            & ycoord(1,0:nr),2,3,'X')
ELSE
   xcoord(:,0:nr) = 0.5*(gxcoordglb(:,0:nr)+gxcoordglb(:,1:nr+1))
   ycoord(:,0:nr) = 0.5*(gycoordglb(:,0:nr)+gycoordglb(:,1:nr+1))
   CALL grid_spacings_curv(delglb(0:nc,0:nr),xcoord(0:nc,0:nr),&
                         & xcoord(1:nc+1,0:nr),ycoord(0:nc,0:nr),&
                         & ycoord(1:nc+1,0:nr),2,3)
ENDIF
delglb(nc+1,0:nr) = delglb(nc,0:nr)
delglb(1:nc+1,nr+1) = delglb(1:nc+1,nr)
delglb(0,nr+1) = delglb(0,nr)
CALL distribute_mod(delglb,delxatc,lbounds(1:2),nhdist,iarr_delxatc,0.0,&
                  & mglevel=0)

!---U-nodes
IF (iopt_grid_htype.LE.2) THEN
   ycoord(1,0:nr) = 0.5*(gycoordglb(1,0:nr)+gycoordglb(1,1:nr+1))
   CALL grid_spacings_rectang(delglb(1:nc,0:nr),gdelxglb(1:nc),&
                            & ycoord(1,0:nr),2,3,'X')
ELSE
   xcoord(0:nc,0:nr) = 0.25*(gxcoordglb(0:nc,0:nr)+gxcoordglb(1:nc+1,0:nr)+&
                           & gxcoordglb(0:nc,1:nr+1)+gxcoordglb(1:nc+1,1:nr+1))
   ycoord(0:nc,0:nr) = 0.25*(gycoordglb(0:nc,0:nr)+gycoordglb(1:nc+1,0:nr)+&
                           & gycoordglb(0:nc,1:nr+1)+gycoordglb(1:nc+1,1:nr+1))
   CALL grid_spacings_curv(delglb(1:nc,0:nr),xcoord(0:nc-1,0:nr),&
                         & xcoord(1:nc,0:nr),ycoord(0:nc-1,0:nr),&
                         & ycoord(1:nc,0:nr),2,3)
ENDIF
delglb(0,0:nr) = delglb(1,0:nr)
delglb(nc+1,0:nr) = delglb(nc,0:nr)
delglb(0:nc+1,nr+1) = delglb(0:nc+1,nr)
CALL distribute_mod(delglb,delxatu,lbounds(1:2),nhdist,iarr_delxatu,0.0,&
                  & mglevel=0)

!---UV-nodes
IF (iopt_grid_htype.LE.2) THEN
   CALL grid_spacings_rectang(delglb(1:nc,0:nr+1),gdelxglb(1:nc),&
                            & gycoordglb(1,:),2,3,'X')
ELSE
   xcoord(0:nc,:) = 0.5*(gxcoordglb(0:nc,:)+gxcoordglb(1:nc+1,:))
   ycoord(0:nc,:) = 0.5*(gycoordglb(0:nc,:)+gycoordglb(1:nc+1,:))
   CALL grid_spacings_curv(delglb(1:nc,0:nr+1),xcoord(0:nc-1,:),&
                         & xcoord(1:nc,:),ycoord(0:nc-1,:),ycoord(1:nc,:),2,3)
ENDIF
delglb(0,:) = delglb(1,:)
delglb(nc+1,:) = delglb(nc,:)
CALL distribute_mod(delglb,delxatuv,lbounds(1:2),nhdist,iarr_delxatuv,0.0,&
                  & mglevel=0)

!
!3.2 Y-direction
!---------------
!
!---U-nodes
IF (iopt_grid_htype.LE.2) THEN
   ycoord(1,0:nr) = 0.5*(gycoordglb(1,0:nr)+gycoordglb(1,1:nr+1))
   CALL grid_spacings_rectang(delglb(0:nc+1,0:nr),gdelyglb(0:nr),&
                            & ycoord(1,0:nr),2,3,'Y')
ELSE
   CALL grid_spacings_curv(delglb(0:nc+1,0:nr),gxcoordglb(:,0:nr),&
                         & gxcoordglb(:,1:nr+1),gycoordglb(:,0:nr),&
                         & gycoordglb(:,1:nr+1),2,3)
ENDIF
delglb(:,nr+1) = delglb(:,nr)
CALL distribute_mod(delglb,delyatu,lbounds(1:2),nhdist,iarr_delyatu,0.0,&
                  & mglevel=0)

!---C-nodes
IF (iopt_grid_htype.LE.2) THEN
   ycoord(1,0:nr) = 0.5*(gycoordglb(1,0:nr)+gycoordglb(1,1:nr+1))
   CALL grid_spacings_rectang(delglb(0:nc,0:nr),gdelyglb(0:nr),ycoord(1,0:nr),&
                            & 2,3,'Y')
ELSE
   xcoord(0:nc,:) = 0.5*(gxcoordglb(0:nc,:)+gxcoordglb(1:nc+1,:))
   ycoord(0:nc,:) = 0.5*(gycoordglb(0:nc,:)+gycoordglb(1:nc+1,:))
   CALL grid_spacings_curv(delglb(0:nc,0:nr),xcoord(0:nc,0:nr),&
                         & xcoord(0:nc,1:nr+1),ycoord(0:nc,0:nr),&
                         & ycoord(0:nc,1:nr+1),2,3)
ENDIF
delglb(0:nc,nr+1) = delglb(0:nc,nr)
delglb(nc+1,1:nr+1) = delglb(nc,1:nr+1)
delglb(nc+1,0) = delglb(nc,0)
CALL distribute_mod(delglb,delyatc,lbounds(1:2),nhdist,iarr_delyatc,0.0,&
                  & mglevel=0)

!---V-nodes
IF (iopt_grid_htype.LE.2) THEN
   CALL grid_spacings_rectang(delglb(0:nc,1:nr),gdelyglb(1:nr),&
                            & gycoordglb(1,1:nr),2,3,'Y')
ELSE
   xcoord(0:nc,0:nr) = 0.25*(gxcoordglb(0:nc,0:nr)+gxcoordglb(1:nc+1,0:nr)+&
                           & gxcoordglb(0:nc,1:nr+1)+gxcoordglb(1:nc+1,1:nr+1))
   ycoord(0:nc,0:nr) = 0.25*(gycoordglb(0:nc,0:nr)+gycoordglb(1:nc+1,0:nr)+&
                           & gycoordglb(0:nc,1:nr+1)+gycoordglb(1:nc+1,1:nr+1))
   CALL grid_spacings_curv(delglb(0:nc,1:nr),xcoord(0:nc,0:nr-1),&
                         & xcoord(0:nc,1:nr),ycoord(0:nc,0:nr-1),&
                         & ycoord(0:nc,1:nr),2,3)
ENDIF
delglb(0:nc,0) = delglb(0:nc,1)
delglb(0:nc,nr+1) = delglb(0:nc,nr)
delglb(nc+1,0:nr+1) = delglb(nc,0:nr+1)
CALL distribute_mod(delglb,delyatv,lbounds(1:2),nhdist,iarr_delyatv,0.0,&
                  & mglevel=0)

!---UV-nodes
IF (iopt_grid_htype.LE.2) THEN
   CALL grid_spacings_rectang(delglb(0:nc+1,1:nr),gdelyglb(1:nr),&
                            & gycoordglb(1,1:nr),2,3,'Y')
ELSE
   xcoord(:,0:nr) = 0.5*(gxcoordglb(:,0:nr)+gxcoordglb(:,1:nr+1))
   ycoord(:,0:nr) = 0.5*(gycoordglb(:,0:nr)+gycoordglb(:,1:nr+1))
   CALL grid_spacings_curv(delglb(0:nc+1,1:nr),xcoord(:,0:nr-1),xcoord(:,1:nr),&
                         & ycoord(:,0:nr-1),ycoord(:,1:nr),2,3)
ENDIF
delglb(:,0) = delglb(:,1)
delglb(:,nr+1) = delglb(:,nr)
CALL distribute_mod(delglb,delyatuv,lbounds(1:2),nhdist,iarr_delyatuv,0.0,&
                  & mglevel=0)

!
!4. Rotation angle
!------------------
!
!4.1. Rectangular grid
!----------------------
!

IF (nhtype.LE.2) THEN
   IF (surfacegrids(igrd_model,1)%rotated) THEN
      gangleatc = degtorad*surfacegrids(igrd_model,1)%gridangle
   ELSE
      gangleatc = 0.0
   ENDIF

!
!4.2 Curvilinear grid
!--------------------
!

ELSEIF (nhtype.EQ.3) THEN
   delx = 0.25*(gxcoord(2:ncloc+1,1:nrloc)+gxcoord(2:ncloc+1,2:nrloc+1)-&
              & gxcoord(1:ncloc,1:nrloc)-gxcoord(1:ncloc,2:nrloc+1))
   dely = 0.25*(gycoord(2:ncloc+1,1:nrloc)+gycoord(2:ncloc+1,2:nrloc+1)-&
              & gycoord(1:ncloc,1:nrloc)-gycoord(1:ncloc,2:nrloc+1))
   IF (iopt_grid_sph.EQ.1) THEN
      array2dc = y0dat+0.125*degtorad*(3.0*gycoord(2:ncloc+1,1:nrloc)+&
               & 3.0*gycoord(2:ncloc+1,2:nrloc+1)+gycoord(1:ncloc,1:nrloc)+&
               & gycoord(1:ncloc,2:nrloc+1))
      delx = delx*COS(array2dc)
   ENDIF
   CALL complex_polar(delx,dely,xpha=gangleatc(1:ncloc,1:nrloc))
ENDIF

!
!4.3 Halo exchange
!--------------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL exchange_mod(gangleatc,(/0,0/),(/1,1,1,1/),0)
ENDIF

!
!5. Coordinate arrays in reference (geographical) coordinates
!------------------------------------------------------------
!
!5.1 Rotated grid
!----------------
!

IF (nhtype.LE.2.AND.surfacegrids(igrd_model,1)%rotated) THEN

!  ---apply transformation
   CALL grid_rotation_transf(gxcoordglb,gycoordglb,-1,2,igrd_model)

!  ---local arrays
   gxcoord = gxcoordglb(nc1loc-1:nc2loc+1,nr1loc-1:nr2loc+1)
   gycoord = gycoordglb(nc1loc-1:nc2loc+1,nr1loc-1:nr2loc+1)

!
!5.2 Curvilinear grid
!--------------------
!

ELSEIF (nhtype.EQ.3) THEN

   gxcoord = gxcoordglb(nc1loc-1:nc2loc+1,nr1loc-1:nr2loc+1)
   gycoord = gycoordglb(nc1loc-1:nc2loc+1,nr1loc-1:nr2loc+1)

ENDIF 

!
!5.3 At grid centers
!-------------------
!

CALL UVarr_at_C(gxcoordglb(1:nc+1,1:nr+1),gxcoordglbatc,0,0,(/1,1,nz/),&
             & (/nc+1,nr+1,nz/),1,iarr_gxcoordglb,.TRUE.)
CALL UVarr_at_C(gycoordglb(1:nc+1,1:nr+1),gycoordglbatc,0,0,(/1,1,nz/),&
             & (/nc+1,nr+1,nz/),1,iarr_gycoordglb,.TRUE.)

!
!6. Sigma coordinates
!--------------------
!
!6.1 Horizontally uniform grid
!-----------------------------
!

IF (iopt_grid_vtype.LT.3) THEN

   k_610: DO k=1,nz+1
!     ---UW-nodes
      gscoordatuw(:,:,k) = gsigcoordatw(k)
!     ---VW-nodes
      gscoordatvw(:,:,k) = gsigcoordatw(k)
!     ---UVW-nodes
      gscoordatuvw(:,:,k) = gsigcoordatw(k)

   ENDDO k_610

!
!6.2 Non-uniform grid
!--------------------
!

ELSE

!  ---UW-nodes
   CALL Warr_at_UW(gscoordatw,gscoordatuw(1:ncloc+1,:,:),3,(/0,0,1/),&
                & (/ncloc+1,nrloc+1,nz+1/),1,iarr_gscoordatuw,.TRUE.)
!  ---VW-nodes
   CALL Warr_at_VW(gscoordatw,gscoordatvw(:,1:nrloc+1,:),3,(/0,0,1/),&
                & (/ncloc+1,nrloc+1,nz+1/),1,iarr_gscoordatvw,.TRUE.)
!  ---UVW-nodes
   CALL Carr_at_UV(gscoordatw,gscoordatuvw(1:ncloc+1,1:nrloc+1,:),1,1,&
                & (/0,0,1/),(/ncloc+1,nrloc+1,nz+1/),1,iarr_gscoordatuvw,.TRUE.)

!
!6.3 Exchange halo sections
!----------------------------
!

   IF (iopt_MPI.EQ.1) THEN
      lbounds = (/0,0,1/); nhexch = (/1,0,0,0/)
      CALL exchange_mod(gscoordatuw,lbounds,nhexch,iarr_gscoordatuw)
      nhexch = (/0,0,1,0/)
      CALL exchange_mod(gscoordatvw,lbounds,nhexch,iarr_gscoordatvw)
      nhexch = (/1,0,1,0/)
      CALL exchange_mod(gscoordatuvw,lbounds,nhexch,iarr_gscoordatuvw)
   ENDIF

ENDIF

!
!6.4 Sigma at non-staggered nodes
!--------------------------------
!
!6.4.1 Horizontally uniform grid
!-------------------------------
!

IF (iopt_grid_vtype.LT.3) THEN

   k_641: DO k=1,nz
!     ---C-nodes
      gscoordatc(:,:,k) = gsigcoordatc(k)
!     ---U-nodes
      gscoordatu(:,:,k) = gsigcoordatc(k)
!     ---V-nodes
      gscoordatv(:,:,k) = gsigcoordatc(k)
   ENDDO k_641

ELSE

!  ---C-nodes
   IF (iopt_grid_vtype_transf.EQ.21) THEN
      k_642: DO k=1,nz
         sigreg = (k-0.5)/REAL(nz)
         c1 = (1.0-b_SH)*SINH(theta_SH*(sigreg-1.0))/SINH(theta_SH)
         c2 = 0.5*b_SH*(TANH(theta_SH*(sigreg-0.5))/TANH(0.5*theta_SH)-1.0)
         c3 = c1 + c2 + 1.0 - sigreg
         j_6421: DO j=0,nrloc+1
         i_6421: DO i=0,ncloc+1
            IF (seapoint(i,j)) THEN
               f = MERGE(((depmeanatc(i,j)-hcrit_SH)*c3)/depmeanatc(i,j),0.0,&
                        & depmeanatc(i,j).GT.hcrit_SH)
               gscoordatc(i,j,k) = sigreg + f
            ELSE
               gscoordatc(i,j,k) = 0.0
            ENDIF
         ENDDO i_6421
         ENDDO j_6421
      ENDDO k_642
   ELSE
      CALL Warr_at_C(gscoordatw,gscoordatc,(/0,0,1/),(/ncloc+1,nrloc+1,nz+1/),&
                   & 1,iarr_gscoordatc,.TRUE.)
   ENDIF

!  ---U-nodes
   CALL UWarr_at_U(gscoordatuw,gscoordatu,3,3,(/0,0,1/),&
                & (/ncloc+1,nrloc+1,nz+1/),1,iarr_gscoordatuw,.TRUE.)

!  ---V-nodes
   CALL VWarr_at_V(gscoordatvw,gscoordatv,3,3,(/0,0,1/),&
                & (/ncloc+1,nrloc+1,nz+1/),1,iarr_gscoordatvw,.TRUE.)

ENDIF

!
!7. Grid resolution
!------------------
!
!---mean grid resolution
CALL sum2_vars(delxatc(1:ncloc,1:nrloc),delxatc_mean,nhdims,'0  ',&
             & iarr_delxatc,commall=.TRUE.)
delxatc_mean = delxatc_mean/REAL(nc*nr)
CALL sum2_vars(delyatc(1:ncloc,1:nrloc),delyatc_mean,nhdims,'0  ',iarr_delyatc,&
             & commall=.TRUE.)
delyatc_mean = delyatc_mean/REAL(nc*nr)
delhatc_mean = 0.5*(delxatc_mean+delyatc_mean)

!---maximum grid resolution
CALL max_vars(delxatc(1:ncloc,1:nrloc),delxatc_max,iarr_delxatc,&
            & mask=seapoint(1:ncloc,1:nrloc),commall=.TRUE.)
CALL max_vars(delyatc(1:ncloc,1:nrloc),delyatc_max,iarr_delyatc,&
            & mask=seapoint(1:ncloc,1:nrloc),commall=.TRUE.)

!---minimum grid resolution
CALL min_vars(delxatc(1:ncloc,1:nrloc),delxatc_min,iarr_delxatc,&
            & mask=seapoint(1:ncloc,1:nrloc),commall=.TRUE.)
CALL min_vars(delyatc(1:ncloc,1:nrloc),delyatc_min,iarr_delyatc,&
            & mask=seapoint(1:ncloc,1:nrloc),commall=.TRUE.)

!
!8. Grid area
!------------
!

garea = delxatc(1:ncloc,1:nrloc)*delyatc(1:ncloc,1:nrloc)
CALL sum2_vars(garea,gareatot,nhdims,'0  ',iarr_garea,commall=.TRUE.)

!
!9. Longitude and latitude
!-------------------------
!

IF (iopt_grid_sph.EQ.0) THEN
   gxlon = degtorad*dlon_ref
   gylat = degtorad*dlat_ref
ELSE
   CALL UVarr_at_C(gxcoord,gxlon,0,0,(/0,0,nz/),(/ncloc+1,nrloc+1,nz/),1,&
                 & iarr_gxcoord,.TRUE.)
   CALL UVarr_at_C(gycoord,gylat,0,0,(/0,0,nz/),(/ncloc+1,nrloc+1,nz/),1,&
                 & iarr_gycoord,.TRUE.)
   gxlon = gxlon*degtorad
   gylat = gylat*degtorad
ENDIF

!
!10. Coriolis frequency
!----------------------
!

rotfreq = 4.0*pi/86164.0
IF (iopt_grid_sph.EQ.0) THEN
   rlatref = degtorad*dlat_ref
   coriolatu = rotfreq*SIN(rlatref)
   coriolatv = rotfreq*SIN(rlatref)
ELSEIF (iopt_grid_sph.EQ.1) THEN
   CALL UVarr_at_U(gycoord(1:ncloc+1,0:nrloc+1),array2du,0,0,&
                & (/1,0,nz/),(/ncloc+1,nrloc+1,nz/),1,iarr_gycoord,.TRUE.)
   CALL UVarr_at_V(gycoord(0:ncloc+1,1:nrloc+1),array2dv,0,0,&
                & (/0,1,nz/),(/ncloc+1,nrloc+1,nz/),1,iarr_gycoord,.TRUE.)
   array2du = array2du*degtorad
   array2dv = array2dv*degtorad
   coriolatu = rotfreq*SIN(array2du)
   coriolatv = rotfreq*SIN(array2dv)
ENDIF 

!
!11. Acceleration of gravity
!---------------------------
!
!11.1 Interior domain
!--------------------
!

IF (ABS(gacc_ref-float_fill).GT.float_fill_eps) THEN
   gaccatc = gacc_ref; gaccatu = gacc_ref; gaccatv = gacc_ref
   gacc_mean = gacc_ref
ELSE
   IF (iopt_grid_sph.EQ.0) THEN
      gacc_mean = 9.78032 + 0.005172*(SIN(rlatref))**2 - &
                & 0.00006*(SIN(2.0*rlatref))**2
      gaccatc = gacc_mean; gaccatu = gacc_mean; gaccatv = gacc_mean
   ELSEIF (iopt_grid_sph.EQ.1) THEN
      CALL UVarr_at_C(gycoord(1:ncloc+1,1:nrloc+1),array2dc,0,0,(/1,1,nz/),&
                   & (/ncloc+1,nrloc+1,nz/),1,iarr_gycoord,.TRUE.)
      array2dc = degtorad*array2dc
      gaccatc(1:ncloc,1:nrloc) = 9.78032 + &
              & 0.005172*(SIN(array2dc))**2 - 0.00006*(SIN(2.0*array2dc))**2
      gaccatu(1:ncloc+1,:) = 9.78032 + 0.005172*(SIN(array2du(:,1:nrloc)))**2 -&
              & 0.00006*(SIN(2.0*array2du(:,1:nrloc)))**2
      gaccatv(:,1:nrloc+1) = 9.78032 + 0.005172*(SIN(array2dv(1:ncloc,:)))**2 -&
              & 0.00006*(SIN(2.0*array2dv(1:ncloc,:)))**2
      CALL sum2_vars(gaccatc(1:ncloc,1:nrloc),gacc_mean,nhdims,'0  ',&
                   & iarr_gaccatc,commall=.TRUE.)
      gacc_mean = gacc_mean/REAL(nc*nr)
   ENDIF
ENDIF

!
!11.2 Exchange halo sections
!---------------------------
!

IF ((iopt_MPI.EQ.1).AND.(ABS(gacc_ref-float_fill).LE.float_fill_eps).AND.&
  & (iopt_grid_sph.EQ.1)) THEN
   lbounds(1:2) = 0; nhexch = 1
   CALL exchange_mod(gaccatc,lbounds(1:2),nhexch,iarr_gaccatc)
   lbounds(1:2) = (/0,1/); nhexch = (/1,1,0,0/)
   CALL exchange_mod(gaccatu,lbounds(1:2),nhexch,iarr_gaccatu)
   lbounds(1:2) = (/1,0/); nhexch = (/0,0,1,1/)
   CALL exchange_mod(gaccatv,lbounds(1:2),nhexch,iarr_gaccatv)
ENDIF

!
!12. Relaxation factor for momentum advection
!--------------------------------------------
!

IF (iopt_obc_advrlx.EQ.1) THEN

!  ---check relaxation distance
   hgrid = MAX(delxatc_mean,delyatc_mean)
   CALL error_lbound_var(distrlx_obc,'distrlx_obc',hgrid,.TRUE.)

!  ---allocate
   ALLOCATE(xcoordobc(nosbu+nosbv,1),STAT=errstat)
   CALL error_alloc('xcoordobc',2,(/nosbu+nosbv,1/),kndrtype)
   ALLOCATE(ycoordobc(nosbu+nosbv,1),STAT=errstat)
   CALL error_alloc('ycoordobc',2,(/nosbu+nosbv,1/),kndrtype)

!  ---positions of open boundaries
   l = 0
   ii_1210: DO ii=1,nosbu
      i = iobu(ii); j = jobu(ii); l = l + 1
      xcoordobc(l,1) = 0.5*(gxcoordglb(i,j)+gxcoordglb(i,j+1))
      ycoordobc(l,1) = 0.5*(gycoordglb(i,j)+gycoordglb(i,j+1))
   ENDDO ii_1210
   jj_1220: DO jj=1,nosbv
      i = iobv(jj); j = jobv(jj); l = l + 1
      xcoordobc(l,1) = 0.5*(gxcoordglb(i,j)+gxcoordglb(i+1,j))
      ycoordobc(l,1) = 0.5*(gycoordglb(i,j)+gycoordglb(i+1,j))
   ENDDO jj_1220

!  ---relaxation factor at U-nodes
   j_1230: DO j=1,nrloc
   i_1230: DO i=1,ncloc
      IF (node2du(i,j).EQ.2.OR.node2du(i,j).EQ.4) THEN
         x = 0.5*(gxcoord(i,j)+gxcoord(i,j+1))
         y = 0.5*(gycoord(i,j)+gycoord(i,j+1))
         CALL distance_minloc(x,y,xcoordobc,ycoordobc,nosbu+nosbv,1,2,3,&
                            & distmin=distmin)
         rlxobcatu(i,j) = MIN(distmin/distrlx_obc,1.0)
      ENDIF
   ENDDO i_1230
   ENDDO j_1230

!  ---relaxation factor at V-nodes
   j_1240: DO j=1,nrloc
   i_1240: DO i=1,ncloc
      IF (node2dv(i,j).EQ.2.OR.node2dv(i,j).EQ.4) THEN
         x = 0.5*(gxcoord(i,j)+gxcoord(i+1,j))
         y = 0.5*(gycoord(i,j)+gycoord(i+1,j))
         CALL distance_minloc(x,y,xcoordobc,ycoordobc,nosbu+nosbv,1,2,3,&
                            & distmin=distmin)
         rlxobcatv(i,j) = MIN(distmin/distrlx_obc,1.0)
      ENDIF
   ENDDO i_1240
   ENDDO j_1240

ENDIF

!
!13. Deallocate arrays
!---------------------
!

DEALLOCATE (array2dc,array2du,array2dv,delglb,delx,dely,xcoord,ycoord)
IF (iopt_obc_advrlx.EQ.1) DEALLOCATE(xcoordobc,ycoordobc)

CALL log_timer_out()


RETURN

END SUBROUTINE grid_arrays

!========================================================================

SUBROUTINE open_boundary_arrays
!************************************************************************
!
! *open_boundary_arrays* Allocate and evaluate arrays at open boundaries
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Grid_Arrays.f90  V2.6
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - combi_mod, combine_stats_glb, error_alloc, local_proc
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: local_proc
USE paral_comms, ONLY: combine_mod, combine_stats_glb
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: i, ii, iiloc, iob, iproc, j, jj, jjloc, job
INTEGER, ALLOCATABLE, DIMENSION(:) :: indu, indv, jndu, jndv
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: node2dglb
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: flag

procname(pglev+1) = 'open_boundary_arrays'
CALL log_timer_in()

!
!1. Allocate work space arrays
!-----------------------------
!

ALLOCATE (node2dglb(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo),STAT=errstat)
CALL error_alloc('node2dglb',2,(/nc+2*nhalo,nr+2*nhalo/),kndint)
ALLOCATE (indu(2*nobu),STAT=errstat)
CALL error_alloc('indu',1,(/2*nobu/),kndint)
ALLOCATE (jndu(2*nobu),STAT=errstat)
CALL error_alloc('jndu',1,(/2*nobu/),kndint)
ALLOCATE (indv(2*nobv),STAT=errstat)
CALL error_alloc('indv',1,(/2*nobv/),kndint)
ALLOCATE (jndv(2*nobv),STAT=errstat)
CALL error_alloc('jndv',1,(/2*nobv/),kndint)
ALLOCATE (flag(nc,nr),STAT=errstat)
CALL error_alloc('flag',2,(/nc,nr/),kndlog)

!
!2. Global index arrays at X- and Y-open boundaries
!--------------------------------------------------
!
!2.1 Global number and locations
!-------------------------------
!
!---number of X-nodes
CALL combine_mod(node2dglb,node2dv,(/1-nhalo,1-nhalo/),iarr_node2dv,0,&
               & mglevel=0,commall=.TRUE.)
iob = 0; flag = .FALSE.
ii_210: DO ii=1,nobu
   i = iobu(ii); j = jobu(ii)
!  ---corner node below
   IF (.NOT.flag(i,j).AND.ANY(node2dglb(i-1:i,j).EQ.2)) THEN
      iob = iob + 1
      indu(iob) = i; jndu(iob) = j
      flag(i,j) = .TRUE.
   ENDIF
!  ---corner node above
   IF (.NOT.flag(i,j+1).AND.ANY(node2dglb(i-1:i,j+1).EQ.2)) THEN
      iob = iob + 1
      indu(iob) = i; jndu(iob) = j+1
      flag(i,j+1) = .TRUE.
   ENDIF
ENDDO ii_210
nobx = iob

CALL combine_mod(node2dglb,node2du,(/1-nhalo,1-nhalo/),iarr_node2du,0,&
               & mglevel=0,commall=.TRUE.)

job = 0; flag = .FALSE.
jj_210: DO jj=1,nobv
   i = iobv(jj); j = jobv(jj)
!  ---corner node to the left
   IF (.NOT.flag(i,j).AND.ANY(node2dglb(i,j-1:j).EQ.2)) THEN
      job = job + 1
      indv(job) = i; jndv(job) = j
      flag(i,j) = .TRUE.
   ENDIF
!  ---corner node to the right
   IF (.NOT.flag(i+1,j).AND.ANY(node2dglb(i+1,j-1:j).EQ.2)) THEN
      job = job + 1
      indv(job) = i+1; jndv(job) = j
      flag(i+1,j) = .TRUE.
   ENDIF
ENDDO jj_210
noby = job

!
!2.2 Allocate arrays
!-------------------
!

ALLOCATE (iobx(nobx),STAT=errstat)
CALL error_alloc('iobx',1,(/nobx/),kndint)

ALLOCATE (jobx(nobx),STAT=errstat)
CALL error_alloc('jobx',1,(/nobx/),kndint)

ALLOCATE (ioby(noby),STAT=errstat)
CALL error_alloc('ioby',1,(/noby/),kndint)

ALLOCATE (joby(noby),STAT=errstat)
CALL error_alloc('joby',1,(/noby/),kndint)

ALLOCATE (westobu(nobu),STAT=errstat)
CALL error_alloc('westobu',1,(/nobu/),kndlog)
IF (nobu.GT.0) westobu = .FALSE.

ALLOCATE (soutobv(nobv),STAT=errstat)
CALL error_alloc('soutobv',1,(/nobv/),kndlog)
IF (nobv.GT.0) soutobv = .FALSE.

ALLOCATE (westobx(nobx),STAT=errstat)
CALL error_alloc('westobx',1,(/nobx/),kndlog)

ALLOCATE (soutoby(noby),STAT=errstat)
CALL error_alloc('soutoby',1,(/noby/),kndlog)

ALLOCATE (indexobuprocs(nobu,nprocs),STAT=errstat)
CALL error_alloc('indexobuprocs',2,(/nobu,nprocs/),kndint)

ALLOCATE (indexobvprocs(nobv,nprocs),STAT=errstat)
CALL error_alloc('indexobvprocs',2,(/nobv,nprocs/),kndint)

ALLOCATE (indexobxprocs(nobx,nprocs),STAT=errstat)
CALL error_alloc('indexobxprocs',2,(/nobx,nprocs/),kndint) 

ALLOCATE (indexobyprocs(noby,nprocs),STAT=errstat)
CALL error_alloc('indexoby',2,(/noby,nprocs/),kndint) 

ALLOCATE (indexobu(nobu),STAT=errstat)
CALL error_alloc('indexobu',1,(/nobu/),kndint)
indexobu = 0

ALLOCATE (indexobv(nobv),STAT=errstat)
CALL error_alloc('indexobv',1,(/nobv/),kndint)
indexobv = 0

ALLOCATE (indexobx(nobx),STAT=errstat)
CALL error_alloc('indexobx',1,(/nobx/),kndint) 
indexobx = 0

ALLOCATE (indexoby(noby),STAT=errstat)
CALL error_alloc('indexoby',1,(/noby/),kndint) 
indexoby = 0

!
!2.3 Store global index arrays
!-----------------------------
!

IF (nobx.GT.0) THEN
   iobx = indu(1:nobx); jobx = jndu(1:nobx)
ENDIF
IF (noby.GT.0) THEN
   ioby = indv(1:noby); joby = jndv(1:noby)
ENDIF

!
!3. Open boundary locations (local)
!----------------------------------
!

iproc_300: DO iproc=1,nprocs

!
!3.1 U-nodes
!-----------
!
!  ---internal locations
   iiloc = 0; nosbuprocs(iproc) = 0
   ii_311: DO ii=1,nobu
      i = iobu(ii); j = jobu(ii)
      IF (local_proc(i,j,iproc=iproc)) THEN
         iiloc = iiloc + 1
         indu(iiloc) = i-nc1procs(iproc)+1
         jndu(iiloc) = j-nr1procs(iproc)+1
         indexobuprocs(iiloc,iproc) = ii
         IF (ii.LE.nosbu) nosbuprocs(iproc) = nosbuprocs(iproc) + 1
      ENDIF
   ENDDO ii_311
   nobuprocs(iproc) = iiloc
   nrvbuprocs(iproc) = nobuprocs(iproc) - nosbuprocs(iproc)
   IF (idloc.EQ.idprocs(iproc)) THEN
      nobuloc = iiloc
      nosbuloc = nosbuprocs(iproc)
      nrvbuloc = nrvbuprocs(iproc)
   ENDIF

!  ---add locations in eastern halo
   ii_312: DO ii=1,nobu
      i = iobu(ii); j = jobu(ii)
      IF (i.EQ.(nc2procs(iproc)+1).AND.&
        & j.GE.nr1procs(iproc).AND.j.LE.nr2procs(iproc)) THEN
         iiloc = iiloc + 1
         indu(iiloc) = ncprocs(iproc)+1
         jndu(iiloc) = j-nr1procs(iproc)+1
         indexobuprocs(iiloc,iproc) = ii
      ENDIF
   ENDDO ii_312

!  ---allocate and store   
   IF (idloc.EQ.idprocs(iproc)) THEN
      nobuloc_ext = iiloc
      ALLOCATE(iobuloc(nobuloc_ext),STAT=errstat)
      CALL error_alloc('iobuloc',1,(/nobuloc_ext/),kndint)
      ALLOCATE(jobuloc(nobuloc_ext),STAT=errstat)
      CALL error_alloc('jobuloc',1,(/nobuloc_ext/),kndint)
      IF (nobuloc_ext.GT.0) THEN
         iobuloc = indu(1:nobuloc_ext); jobuloc = jndu(1:nobuloc_ext)
         indexobu(1:nobuloc_ext) = indexobuprocs(1:nobuloc_ext,iproc)
      ENDIF
   ENDIF

!
!3.2 V-nodes
!-----------
!
!  ---internal locations
   jjloc = 0; nosbvprocs(iproc) = 0
   jj_321: DO jj=1,nobv
      i = iobv(jj); j = jobv(jj)
      IF (local_proc(i,j,iproc=iproc)) THEN
         jjloc = jjloc + 1
         indv(jjloc) = i-nc1procs(iproc)+1
         jndv(jjloc) = j-nr1procs(iproc)+1
         indexobvprocs(jjloc,iproc) = jj
         IF (jj.LE.nosbv) nosbvprocs(iproc) = nosbvprocs(iproc) + 1
      ENDIF
   ENDDO jj_321
   nobvprocs(iproc) = jjloc
   nrvbvprocs(iproc) = nobvprocs(iproc) - nosbvprocs(iproc)
   IF (idloc.EQ.idprocs(iproc)) THEN
      nobvloc = jjloc
      nosbvloc = nosbvprocs(iproc)
      nrvbvloc = nrvbvprocs(iproc)
   ENDIF

!  ---add locations in northern halo
   jj_322: DO jj=1,nobv
      i = iobv(jj); j = jobv(jj)
      IF (i.GE.nc1procs(iproc).AND.i.LE.nc2procs(iproc).AND.&
        & j.EQ.(nr2procs(iproc)+1)) THEN
         jjloc = jjloc + 1
         indv(jjloc) = i-nc1procs(iproc)+1
         jndv(jjloc) = nrprocs(iproc)+1
         indexobvprocs(jjloc,iproc) = jj
      ENDIF
   ENDDO jj_322

!  ---allocate and store
   IF (idloc.EQ.idprocs(iproc)) THEN
      nobvloc_ext = jjloc
      ALLOCATE(iobvloc(nobvloc_ext),STAT=errstat)
      CALL error_alloc('iobvloc',1,(/nobvloc_ext/),kndint)
      ALLOCATE(jobvloc(nobvloc_ext),STAT=errstat)
      CALL error_alloc('jobvloc',1,(/nobvloc_ext/),kndint)
      IF (nobvloc_ext.GT.0) THEN
         iobvloc = indv(1:nobvloc_ext); jobvloc = jndv(1:nobvloc_ext)
         indexobv(1:nobvloc_ext) = indexobvprocs(1:nobvloc_ext,iproc)
      ENDIF
   ENDIF

!
!3.3 X-nodes
!-----------
!
!  ---internal locations
   iiloc = 0; nobxprocs(iproc) = 0
   ii_331: DO ii=1,nobx
      i = iobx(ii); j = jobx(ii)
      IF (local_proc(i,j,iproc=iproc)) THEN
         iiloc = iiloc + 1
         indu(iiloc) = i-nc1procs(iproc)+1
         jndu(iiloc) = j-nr1procs(iproc)+1
         indexobxprocs(iiloc,iproc) = ii
      ENDIF
   ENDDO ii_331
   nobxprocs(iproc) = iiloc
   IF (idloc.EQ.idprocs(iproc)) nobxloc = iiloc

!  ---add locations in eastern halo
   ii_332: DO ii=1,nobx
      i = iobx(ii); j = jobx(ii)
      IF (i.EQ.(nc2procs(iproc)+1).AND.&
        & j.GE.nr1procs(iproc).AND.j.LE.nr2procs(iproc)) THEN
         iiloc = iiloc + 1
         indu(iiloc) = ncprocs(iproc)+1
         jndu(iiloc) = j-nr1procs(iproc)+1
         indexobxprocs(iiloc,iproc) = ii
      ENDIF
   ENDDO ii_332

!  ---allocate and store   
   IF (idloc.EQ.idprocs(iproc)) THEN
      nobxloc_ext = iiloc
      ALLOCATE(iobxloc(nobxloc_ext),STAT=errstat)
      CALL error_alloc('iobxloc',1,(/nobxloc_ext/),kndint)
      ALLOCATE(jobxloc(nobxloc_ext),STAT=errstat)
      CALL error_alloc('jobxloc',1,(/nobxloc_ext/),kndint)
      IF (nobxloc_ext.GT.0) THEN
         iobxloc = indu(1:nobxloc_ext); jobxloc = jndu(1:nobxloc_ext)
         indexobx(1:nobxloc_ext) = indexobxprocs(1:nobxloc_ext,iproc)
      ENDIF
   ENDIF

!
!3.4 Y-nodes
!-----------
!
!  ---internal locations
   jjloc = 0; nobyprocs(iproc) = 0
   jj_341: DO jj=1,noby
      i = ioby(jj); j = joby(jj)
      IF (local_proc(i,j,iproc=iproc)) THEN
         jjloc = jjloc + 1
         indv(jjloc) = i-nc1procs(iproc)+1
         jndv(jjloc) = j-nr1procs(iproc)+1
         indexobyprocs(jjloc,iproc) = jj
      ENDIF
   ENDDO jj_341
   nobyprocs(iproc) = jjloc
   IF (idloc.EQ.idprocs(iproc)) nobyloc = jjloc

!  ---add locations in northern halo
   jj_342: DO jj=1,noby
      i = ioby(jj); j = joby(jj)
      IF (i.GE.nc1procs(iproc).AND.i.LE.nc2procs(iproc).AND.&
        & j.EQ.(nr2procs(iproc)+1)) THEN
         jjloc = jjloc + 1
         indv(jjloc) = i-nc1procs(iproc)+1
         jndv(jjloc) = nrprocs(iproc)+1
         indexobyprocs(jjloc,iproc) = jj
      ENDIF
   ENDDO jj_342

!  ---allocate and store   
   IF (idloc.EQ.idprocs(iproc)) THEN
      nobyloc_ext = jjloc
      ALLOCATE(iobyloc(nobyloc_ext),STAT=errstat)
      CALL error_alloc('iobyloc',1,(/nobyloc_ext/),kndint)
      ALLOCATE(jobyloc(nobyloc_ext),STAT=errstat)
      CALL error_alloc('jobyloc',1,(/nobyloc_ext/),kndint)
      IF (nobyloc_ext.GT.0) THEN
         iobyloc = indv(1:nobyloc_ext); jobyloc = jndv(1:nobyloc_ext)
         indexoby(1:nobyloc_ext) = indexobyprocs(1:nobyloc_ext,iproc)
      ENDIF
   ENDIF

ENDDO iproc_300

!
!4. Orientation of open boundary cell faces
!------------------------------------------
!
!4.1 Local
!---------
!
!---U-nodes
ii_411: DO iiloc=1,nobuloc
   i = iobuloc(iiloc); j = jobuloc(iiloc)
   ii = indexobu(iiloc)
   IF ((i+nc1loc-1).EQ.1) THEN
      westobu(ii) = .TRUE.
   ELSEIF (nodeatc(i-1,j).EQ.0) THEN
      westobu(ii) = .TRUE.
   ELSE
      westobu(ii) = .FALSE.
   ENDIF
ENDDO ii_411

!--- V-nodes
jjloc_412: DO jjloc=1,nobvloc
   i = iobvloc(jjloc); j = jobvloc(jjloc)
   jj = indexobv(jjloc)
   IF ((j+nr1loc-1).EQ.1) THEN
      soutobv(jj) = .TRUE.
   ELSEIF (nodeatc(i,j-1).EQ.0) THEN
      soutobv(jj) = .TRUE.
   ELSE
      soutobv(jj) = .FALSE.
   ENDIF
ENDDO jjloc_412

!---X-nodes
iiloc_413: DO iiloc=1,nobxloc
   i = iobxloc(iiloc); j = jobxloc(iiloc)
   ii = indexobx(iiloc)
   IF ((i+nc1loc-1).EQ.1) THEN
      westobx(ii) = .TRUE.
   ELSEIF (node2dv(i-1,j).NE.2) THEN
      westobx(ii) = .TRUE.
   ELSE
      westobx(ii) = .FALSE.
   ENDIF
ENDDO iiloc_413

!---Y-nodes
jjloc_414: DO jjloc=1,nobyloc
   i = iobyloc(jjloc); j = jobyloc(jjloc)
   jj = indexoby(jjloc)
   IF ((j+nr1loc-1).EQ.1) THEN
      soutoby(jj) = .TRUE.
   ELSEIF (node2du(i,j-1).NE.2) THEN
      soutoby(jj) = .TRUE.
   ELSE
      soutoby(jj) = .FALSE.
   ENDIF
ENDDO jjloc_414

!
!4.2 Global arrays
!-----------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL combine_stats_glb(westobu,nobu,nobuprocs,indexobuprocs,iarr_westobu,&
                        & commall=.TRUE.)
   CALL combine_stats_glb(soutobv,nobv,nobvprocs,indexobvprocs,iarr_soutobv,&
                        & commall=.TRUE.)
   CALL combine_stats_glb(westobx,nobx,nobxprocs,indexobxprocs,iarr_westobx,&
                        & commall=.TRUE.)
   CALL combine_stats_glb(soutoby,noby,nobyprocs,indexobyprocs,iarr_soutoby,&
                        & commall=.TRUE.)
ENDIF

!
!5. Deallocate work space arrays
!-------------------------------
!

DEALLOCATE (indu,indv,jndu,jndv,node2dglb,flag)

CALL log_timer_out()


RETURN

END SUBROUTINE open_boundary_arrays

!========================================================================

SUBROUTINE pointer_arrays
!************************************************************************
!
! *pointer_arrays* Cell pointers
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Grid_Arrays.f90  V2.6
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - thin_dams, update_pointer_arrays
!
! Module calls - exchange_mod, local_proc, sum_vars
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE multigrid
USE physpars
USE switches
USE grid_routines, ONLY: local_proc
USE paral_comms, ONLY: exchange_mod
USE paral_utilities, ONLY: sum_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: i, ii, j, jj
INTEGER, DIMENSION(2) :: lbounds
INTEGER, DIMENSION(4) :: nhexch


procname(pglev+1) = 'pointer_arrays'
CALL log_timer_in()

!
!1. Cell pointers
!----------------
!
!1.1 Centers
!-----------
!
!---interior points
seapoint = seapointglb(nc1loc-1:nc2loc+1,nr1loc-1:nr2loc+1)
nodeatc(1:ncloc,1:nrloc) = MERGE(1,0,seapoint(1:ncloc,1:nrloc))
rmaskatc = REAL(nodeatc(1:ncloc,1:nrloc))

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   lbounds = 1-nhalo; nhexch = nhalo
   CALL exchange_mod(nodeatc,lbounds,nhexch,iarr_nodeatc)
ENDIF

!---mask arrays
maskatc_int = nodeatc(1:ncloc,1:nrloc).GT.0

!
!1.2 U-nodes
!-----------
!
!---interior faces
j_121: DO j=1,nrloc
i_121: DO i=1,ncloc
   nodeatu(i,j,:) = COUNT(nodeatc(i-1:i,j).GT.0)
ENDDO i_121
ENDDO j_121

!---open boundaries
ii_122: DO ii=1,nobu
   i = iobu(ii) - nc1loc + 1; j = jobu(ii) - nr1loc + 1
   IF (local_proc(i,j)) THEN
      nodeatu(i,j,:) = MERGE(3,4,ii.LE.nosbu)
   ENDIF
ENDDO ii_122

!---U thin-dams
IF (iopt_thndam.EQ.1) CALL thin_dams('U')

!
!1.3 V-nodes
!-----------
!
!---interior faces
j_131: DO j=1,nrloc
i_131: DO i=1,ncloc
   nodeatv(i,j,:) = COUNT(nodeatc(i,j-1:j).GT.0)
ENDDO i_131
ENDDO j_131

!---open boundaries
jj_132: DO jj=1,nobv
   i = iobv(jj) - nc1loc + 1; j = jobv(jj) - nr1loc + 1
   IF (local_proc(i,j)) THEN
      nodeatv(i,j,:) = MERGE(3,4,jj.LE.nosbv)
   ENDIF
ENDDO jj_132

!---V thin-dams
IF (iopt_thndam.EQ.1) CALL thin_dams('V')

!
!1.4 Other nodes/reset currents
!------------------------------
!

CALL update_pointer_arrays

!
!
!2. Number of active points
!--------------------------
!
!---C-nodes
noseaatcloc = COUNT(seapoint(1:ncloc,1:nrloc))
noseaatc = COUNT(seapointglb(1:nc,1:nr))

!---U-nodes
noseaatuloc = COUNT(node2du(1:ncloc,1:nrloc).GT.0)
CALL sum_vars(noseaatuloc,noseaatu,0,commall=.TRUE.)

!---U-nodes
noseaatvloc = COUNT(node2dv(1:ncloc,1:nrloc).GT.0)
CALL sum_vars(noseaatvloc,noseaatv,0,commall=.TRUE.)

CALL log_timer_out()


RETURN

END SUBROUTINE pointer_arrays

!========================================================================

SUBROUTINE read_grid
!************************************************************************
!
! *read_grid* Read model grid arrays in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Grid_Arrays.f90  V2.5.1
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE depths
USE grid
USE gridpars
USE iopars
USE switches
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_grid'
CALL log_timer_in()

filepars = modfiles(io_modgrd,1,1)

!
!1. File header
!--------------
!

CALL open_filepars(filepars)
CALL read_glbatts_mod(filepars)
numvars = filepars%novars
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL read_varatts_mod(filepars,varatts,numvars)

!
!2. Coordinate arrays
!--------------------
!
!---horizontal
IF (iopt_grid_htype.EQ.2) THEN
   CALL read_vars(gdelxglb(1:nc),filepars,1,varatts=varatts)
   CALL read_vars(gdelyglb(1:nr),filepars,2,varatts=varatts)
   varid = 3
ELSEIF (iopt_grid_htype.EQ.3) THEN
   CALL read_vars(gxcoordglb(1:nc,1:nr),filepars,1,varatts=varatts)
   CALL read_vars(gycoordglb(1:nc,1:nr),filepars,2,varatts=varatts)
   varid = 3
ELSE
   varid = 1
ENDIF

!---vertical
IF (iopt_grid_vtype.EQ.2) THEN
   CALL read_vars(gsigcoordatw,filepars,varid,varatts=varatts)
   varid = varid + 1
ELSEIF (iopt_grid_vtype.EQ.3) THEN
   CALL read_vars(gscoordglb(1:nc-1,1:nr-1,:),filepars,varid,varatts=varatts)
   varid = varid + 1
ENDIF

!
!3. Water depths
!---------------
!

CALL read_vars(depmeanglb(1:nc-1,1:nr-1),filepars,varid,varatts=varatts)
varid = varid + 1

!
!4. Open boundary locations
!--------------------------
!
!---U-nodes
IF (nobu.GT.0) THEN
   CALL read_vars(iobu,filepars,varid,varatts=varatts)
   CALL read_vars(jobu,filepars,varid+1,varatts=varatts)
   varid = varid + 2
ENDIF

!---V-nodes
IF (nobv.GT.0) THEN
   CALL read_vars(iobv,filepars,varid,varatts=varatts)
   CALL read_vars(jobv,filepars,varid+1,varatts=varatts)
ENDIF

!
!5. Finalise
!-----------
!
 
CALL close_filepars(filepars)
modfiles(io_modgrd,1,1) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE read_grid

!========================================================================

SUBROUTINE store_depths_old
!************************************************************************
!
! *store_depths_old* Store total water depths at old time
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Grid_Arrays.f90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - coherens_main, initialise_physical_arrays
!
! External calls -
! 
! Module calls - exchange_mod, Carr_at_U, Carr_at_UV, Carr_at_V
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'store_depths_old'
CALL log_timer_in()

!---C-nodes
WHERE (maskatc_int)
   deptotatc_old(1:ncloc,1:nrloc) = deptotatc(1:ncloc,1:nrloc)
END WHERE

!---U-nodes
WHERE (node2du(0:ncloc+1,1:nrloc).GT.0)
   deptotatu_old = deptotatu(0:ncloc+1,1:nrloc)
END WHERE

!---V-nodes
WHERE (node2dv(1:ncloc,0:nrloc+1).GT.0)
   deptotatv_old = deptotatv(1:ncloc,0:nrloc+1)
END WHERE

CALL log_timer_out()


RETURN

END SUBROUTINE store_depths_old

!========================================================================

SUBROUTINE update_pointer_arrays
!************************************************************************
!
! *update_pointer_arrays* Update pointer arrays at velocity and corner
!                         nodes for weirs or dry cells
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Grid_Arrays.f90  V2.6
!
! Description - set currents to zero at blocking interfaces
!
! Reference -
!
! Calling program - mask_function, pointer_arrays, weirs_mask
!
! External calls -
! 
! Module calls - exchange_mod
!
!************************************************************************
!
USE currents
USE grid
USE gridpars
USE iopars
USE modids
USE switches
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j, k
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhexch


procname(pglev+1) = 'update_pointer_arrays'
CALL log_timer_in()

!
!1. (Re)set pointer arrays
!-------------------------
!
!1.1 U-nodes
!-----------
!
!---exchange halos
IF (iopt_MPI.EQ.1) THEN
   lbounds = (/1-nhalo,1-nhalo,1/); nhexch = nhalo
   CALL exchange_mod(nodeatu,lbounds,nhexch,iarr_nodeatu)
ENDIF

!---2-D case
WHERE (ANY(nodeatu.EQ.2,DIM=3))
   node2du = 2
ELSEWHERE
   node2du = nodeatu(:,:,nz)
END WHERE

!
!1.2 UW-nodes
!------------
!
!---interior faces
j_120: DO j=1-nhalo,nrloc+nhalo
i_120: DO i=1-nhalo,ncloc+nhalo
   IF (ALL(nodeatu(i,j,1:nz-1).EQ.nodeatu(i,j,nz))) THEN
      nodeatuw(i,j,2:nz) = nodeatu(i,j,2:nz)
   ELSE
      k_121: DO k=2,nz
         nodeatuw(i,j,k) = MERGE(2,1,ANY(nodeatu(i,j,k-1:k).EQ.2))
      ENDDO k_121
   ENDIF
ENDDO i_120
ENDDO j_120

!--bottom, surface
nodeatuw(:,:,1) = node2du
nodeatuw(:,:,nz+1) = node2du

!
!1.3 V-nodes
!-----------
!
!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   CALL exchange_mod(nodeatv,lbounds,nhexch,iarr_nodeatv)
ENDIF

!---2-D case
WHERE (ANY(nodeatv.EQ.2,DIM=3))
   node2dv = 2
ELSEWHERE
   node2dv = nodeatv(:,:,nz)
END WHERE

!
!1.4 VW-nodes
!------------
!
!---interior faces
j_140: DO j=1-nhalo,nrloc+nhalo
i_140: DO i=1-nhalo,ncloc+nhalo
   IF (ALL(nodeatv(i,j,1:nz-1).EQ.nodeatv(i,j,nz))) THEN
      nodeatvw(i,j,2:nz) = nodeatv(i,j,2:nz)
   ELSE
      k_141: DO k=2,nz
         nodeatvw(i,j,k) = MERGE(2,1,ANY(nodeatv(i,j,k-1:k).EQ.2))
      ENDDO k_141
   ENDIF
ENDDO i_140
ENDDO j_140

!--bottom, surface
nodeatvw(:,:,1) = node2dv
nodeatvw(:,:,nz+1) = node2dv

!
!1.5 UV-nodes
!------------
!

k_151: DO k=1,nz
j_151: DO j=1,nrloc
i_151: DO i=1,ncloc
   IF (ALL(nodeatu(i,j-1:j,k).LE.1).OR.ALL(nodeatv(i-1:i,j,k).LE.1)) THEN
      nodeatuv(i,j,k) = 0
   ELSEIF (ANY(nodeatu(i,j-1:j,k).GT.2).AND.ANY(nodeatv(i-1:i,j,k).GT.2)) THEN
      IF (ANY(nodeatu(i,j-1:j,k).LE.1).AND.ANY(nodeatv(i-1:i,j,k).LE.1)) THEN
         nodeatuv(i,j,k) = 0
      ELSEIF (ANY(nodeatu(i,j-1:j,k).EQ.2).AND.&
            & ANY(nodeatv(i-1:i,j,k).EQ.2)) THEN
         nodeatuv(i,j,k) = 4
      ENDIF
   ELSEIF (ANY(nodeatu(i,j-1:j,k).GT.2)) THEN
      nodeatuv(i,j,k) = 2
   ELSEIF (ANY(nodeatv(i-1:i,j,k).GT.2)) THEN
      nodeatuv(i,j,k) = 3
   ELSE
      nodeatuv(i,j,k) = 1
   ENDIF
ENDDO i_151
ENDDO j_151
ENDDO k_151

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   CALL exchange_mod(nodeatuv,lbounds,nhexch,iarr_nodeatuv)
ENDIF

!---2-D case
j_152: DO j=1-nhalo,nrloc+nhalo
i_152: DO i=1-nhalo,ncloc+nhalo
   IF (ANY(nodeatuv(i,j,:).EQ.1)) THEN
      node2duv(i,j) = 1
   ELSEIF (ANY(nodeatuv(i,j,:).EQ.4)) THEN
      node2duv(i,j) = 4
   ELSE
      node2duv(i,j) = nodeatuv(i,j,nz)
   ENDIF
ENDDO i_152
ENDDO j_152

!
!2. Blocking currents at blocked cells and interfaces
!----------------------------------------------------
!
!2.1 U-nodes
!-----------
!

WHERE (node2du.EQ.1)
   udvel = 0.0
END WHERE
WHERE (node2du(:,0:nrloc+1).EQ.1)
   umvel = 0.0
END WHERE
WHERE (node2du(0:ncloc+1,1:nrloc).EQ.1)
   umpred = 0.0
END WHERE
k_210: DO k=1,nz
   WHERE (nodeatu(:,:,k).EQ.1)
      uvel(:,:,k) = 0.0
   END WHERE
   WHERE (nodeatu(2-nhalo:ncloc+nhalo,1:nrloc,k).EQ.1)
      ufvel(:,:,k) = 0.0
   END WHERE
ENDDO k_210

!
!2.2 V-nodes
!-----------
!

WHERE (node2dv.EQ.1)
   vdvel = 0.0
END WHERE
WHERE (node2dv(0:ncloc+1,:).EQ.1)
   vmvel = 0.0
END WHERE
WHERE (node2dv(1:ncloc,0:nrloc+1).EQ.1)
   vmpred = 0.0
END WHERE
k_220: DO k=1,nz
   WHERE (nodeatv(:,:,k).EQ.1)
      vvel(:,:,k) = 0.0
   END WHERE
   WHERE (nodeatv(1:ncloc,2-nhalo:nrloc+nhalo,k).EQ.1)
      vfvel(:,:,k) = 0.0
   END WHERE
ENDDO k_220

CALL log_timer_out()


RETURN

END SUBROUTINE update_pointer_arrays

!========================================================================

SUBROUTINE water_depths
!************************************************************************
!
! *water_depths* Update total water depths
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Grid_Arrays.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - initialise_model, pressure_gradient_1d, surface_elevation
!
! External calls - minimum_depths, weirs_depth
! 
! Module calls - Carr_at_U, Carr_at_UV, Carr_at_V, exchange_mod
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE switches
USE timepars
USE array_interp, ONLY: Carr_at_U, Carr_at_UV, Carr_at_V
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, ic, ii, iiloc, j, jc, jj, jjloc, k
INTEGER, DIMENSION(4) :: nhexch


procname(pglev+1) = 'water_depths'
CALL log_timer_in()

!
!1. Water depths at cell centers
!-------------------------------
!
!---store value at previous time
IF (iopt_fld.GT.0) THEN
   WHERE (maskatc_int)
      deptotatc_prev = deptotatc(1:ncloc,1:nrloc)
   END WHERE
ENDIF
   
!---centers
WHERE (nodeatc(0:ncloc+1,0:nrloc+1).GT.0)
   deptotatc = depmeanatc(0:ncloc+1,0:nrloc+1) + zeta
END WHERE

!
!2. Water depths at velocity nodes
!---------------------------------
!
!2.1 Using linear interpolation
!------------------------------
!

IF (iopt_arrint_depths.EQ.1) THEN

!  ---U-nodes
   CALL Carr_at_U(zeta(0:ncloc,1:nrloc),deptotatu(1:ncloc,1:nrloc),1,3,&
               & (/0,1,nz/),(/ncloc,nrloc,nz/),1,iarr_zeta,.TRUE.)
   WHERE (node2du(1:ncloc,1:nrloc).GT.0)
      deptotatu(1:ncloc,1:nrloc) = depmeanatu(1:ncloc,1:nrloc) + &
                                 & deptotatu(1:ncloc,1:nrloc)
   END WHERE
   
!  ---V-nodes
   CALL Carr_at_V(zeta(1:ncloc,0:nrloc),deptotatv(1:ncloc,1:nrloc),1,3,&
               & (/1,0,nz/),(/ncloc,nrloc,nz/),1,iarr_zeta,.TRUE.)
   WHERE (node2dv(1:ncloc,1:nrloc).GT.0)
      deptotatv(1:ncloc,1:nrloc) = depmeanatv(1:ncloc,1:nrloc) + &
                                 & deptotatv(1:ncloc,1:nrloc)
   END WHERE

!
!2.2 Using minimum depth of adjacent cells
!-----------------------------------------
!

ELSEIF (iopt_arrint_depths.EQ.2) THEN

!  ---interior U-nodes
   WHERE (node2du(1:ncloc,1:nrloc).EQ.2)
      deptotatu(1:ncloc,1:nrloc) = MIN(depmeanatc(0:ncloc-1,1:nrloc), &
                                     & depmeanatc(1:ncloc,1:nrloc))
   END WHERE

!  ---U-coastal boundaries
   WHERE (node2du(1:ncloc,1:nrloc).EQ.1)
      deptotatu(1:ncloc,1:nrloc) = MERGE (deptotatc(1:ncloc,1:nrloc),&
                                    & deptotatc(0:ncloc-1,1:nrloc),maskatc_int)
   END WHERE

!  ---U-open boundaries
   iiloc_221: DO iiloc=1,nobuloc
      i = iobuloc(iiloc); j = jobuloc(ii)
      ii =indexobu(iiloc); ic = MERGE(i,i-1,westobu(ii))
      deptotatu(i,j) = deptotatc(ic,j)
   ENDDO iiloc_221
   
!  ---interior V-nodes
   WHERE (node2dv(1:ncloc,1:nrloc).EQ.2)
      deptotatv(1:ncloc,1:nrloc) = MIN(depmeanatc(1:ncloc,0:nrloc-1), &
                                     & depmeanatc(1:ncloc,1:nrloc))
   END WHERE

!  ---V-coastal boundaries
   WHERE (node2dv(1:ncloc,1:nrloc).EQ.1)
      deptotatv(1:ncloc,1:nrloc) = MERGE (deptotatc(1:ncloc,1:nrloc),&
                                     & deptotatc(1:ncloc,0:nrloc-1),maskatc_int)
   END WHERE

!  ---V-open boundaries
   jjloc_222: DO jjloc=1,nobvloc
      i = iobvloc(jjloc); j = jobvloc(jj)
      jj = indexobv(jjloc); jc = MERGE(j,j-1,soutobv(jj))
      deptotatv(i,j) = deptotatc(i,jc)
   ENDDO jjloc_222

ENDIF

!---UV-nodes
CALL Carr_at_UV(zeta,deptotatuv,1,1,(/0,0,nz/),(/ncloc+1,nrloc+1,nz/),1,&
              & iarr_zeta,.TRUE.)
WHERE (node2duv(1:ncloc+1,1:nrloc+1).GT.0)
   deptotatuv = depmeanatuv(1:ncloc+1,1:nrloc+1) + deptotatuv
END WHERE

!
!3. Water depths at weirs/barriers
!------- -------------------------
!

IF (iopt_weibar.EQ.1) CALL weirs_depth

!
!4. Avoid negative depths in case of drying
!------------------------------------------
!

IF (iopt_fld.GT.0) CALL minimum_depths

!
!5. Halo exchanges
!-----------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = (/nhalo,nhalo,1,1/)
   CALL exchange_mod(deptotatu,(/1-nhalo,0/),nhexch,iarr_deptotatu)
   nhexch = (/1,1,nhalo,nhalo/)
   CALL exchange_mod(deptotatv,(/0,1-nhalo/),nhexch,iarr_deptotatv)
ENDIF

!
!6. Vertical grid spacings
!-------------------------
!

IF (nt.EQ.0.OR.corrstep) THEN

!
!6.1 Verically non-staggered nodes
!---------------------------------
!

   k_610: DO k=1,nz

!     ---C-nodes
      WHERE (seapoint)
         delzatc(:,:,k) = deptotatc*(gscoordatw(:,:,k+1)-gscoordatw(:,:,k))
      END WHERE

!     ---U-nodes
      WHERE (nodeatu(0:ncloc+1,1:nrloc,k).GT.0)  
         delzatu(:,:,k) = deptotatu(0:ncloc+1,1:nrloc)*&
                        & (gscoordatuw(:,1:nrloc,k+1)-gscoordatuw(:,1:nrloc,k))
      END WHERE

!     ---V-nodes
      WHERE (nodeatv(1:ncloc,0:nrloc+1,k).GT.0)  
         delzatv(:,:,k) = deptotatv(1:ncloc,0:nrloc+1)*&
                        & (gscoordatvw(1:ncloc,:,k+1)-gscoordatvw(1:ncloc,:,k))
      END WHERE

!     ---UV-nodes
      WHERE (nodeatuv(1:ncloc+1,1:nrloc+1,k).GT.0)  
         delzatuv(:,:,k) = deptotatuv*&
              & (gscoordatuvw(1:ncloc+1,1:nrloc+1,k+1)-&
              &  gscoordatuvw(1:ncloc+1,1:nrloc+1,k))
      END WHERE

   ENDDO k_610

!
!6.2 Vertically staggered nodes
!------------------------------
!

   k_620: DO k=2,nz

!     ---W-nodes
      WHERE (seapoint)
         delzatw(:,:,k) = deptotatc*(gscoordatc(:,:,k)-gscoordatc(:,:,k-1))
      END WHERE

!     ---UW-nodes
      WHERE (nodeatuw(1:ncloc+1,1:nrloc,k).GT.0)  
         delzatuw(:,:,k) = deptotatu(1:ncloc+1,1:nrloc)*&
           & (gscoordatu(1:ncloc+1,1:nrloc,k)-gscoordatu(1:ncloc+1,1:nrloc,k-1))
      END WHERE

!     ---VW-nodes
      WHERE (nodeatvw(1:ncloc,1:nrloc+1,k).GT.0)  
         delzatvw(:,:,k) = deptotatv(1:ncloc,1:nrloc+1)*&
           & (gscoordatv(1:ncloc,1:nrloc+1,k)-gscoordatv(1:ncloc,1:nrloc+1,k-1))
      END WHERE

   ENDDO k_620

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE water_depths

!========================================================================

SUBROUTINE write_grid
!************************************************************************
!
! *write_grid* Write model grid arrays in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Grid_Arrays.f90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, write_atts_mod,
!                write_vars
!
!************************************************************************
!
USE datatypes
USE depths
USE grid
USE gridpars
USE iopars
USE paralpars
USE switches
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN 

procname(pglev+1) = 'write_grid'
CALL log_timer_in()


!
!1. Write file header
!--------------------
!
!---file attributes
CALL set_modfiles_atts(io_modgrd,1,2)
filepars = modfiles(io_modgrd,1,2)
numvars = filepars%novars

!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL set_modvars_atts(io_modgrd,1,2,varatts,numvars)

!---write
CALL write_atts_mod(filepars,varatts,numvars)

!
!2. Coordinate arrays
!--------------------
!
!---horizontal
IF (iopt_grid_htype.EQ.2) THEN
   CALL write_vars(gdelxglb(1:nc),filepars,1,varatts=varatts)
   CALL write_vars(gdelyglb(1:nr),filepars,2,varatts=varatts)
   varid = 3
ELSEIF (iopt_grid_htype.EQ.3) THEN
   CALL write_vars(gxcoordglb(1:nc,1:nr),filepars,1,varatts=varatts)
   CALL write_vars(gycoordglb(1:nc,1:nr),filepars,2,varatts=varatts)
   varid = 3
ELSE
   varid = 1
ENDIF

!---vertical
IF (iopt_grid_vtype.EQ.2) THEN
   CALL write_vars(gsigcoordatw,filepars,varid,varatts=varatts)
   varid = varid + 1
ELSEIF (iopt_grid_vtype.EQ.3) THEN
   CALL write_vars(gscoordglb(1:nc-1,1:nr-1,:),filepars,varid,varatts=varatts)
   varid = varid + 1
ENDIF

!
!3. Water depths
!---------------
!

CALL write_vars(depmeanglb(1:nc-1,1:nr-1),filepars,varid,varatts=varatts)
varid = varid + 1

!
!4. Open boundary locations
!--------------------------
!
!---U-nodes
IF (nobu.GT.0) THEN
   CALL write_vars(iobu,filepars,varid,varatts=varatts)
   CALL write_vars(jobu,filepars,varid+1,varatts=varatts)
   varid = varid + 2
ENDIF

!---V-nodes
IF (nobv.GT.0) THEN
   CALL write_vars(iobv,filepars,varid,varatts=varatts)
   CALL write_vars(jobv,filepars,varid+1,varatts=varatts)
ENDIF

!
!5. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_modgrd,1,2) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE write_grid
