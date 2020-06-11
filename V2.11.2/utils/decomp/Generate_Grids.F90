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
! *Generate_Grids* Utility program for constructing a domain decomposition
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Generate_Grids.f90  V2.10
!
! $Date: 2014-03-13 15:42:28 +0100 (Thu, 13 Mar 2014) $
!
! $Revision: 682 $
!
! Description - 
!
! Reference -
!
! Routines - decomp_mumap, domain_decomp, grid_pointers, reset_partition,
!            set_partition, sqrt_factoring, write_partition
!
!************************************************************************
!

!========================================================================

PROGRAM generate_grids
!************************************************************************
!
! *generate_grids* Test program for generating domain decomposition
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Generate_Grids.F90  V2.10
!  
! Description - 
!
! Reference -
!
! External calls - decomp_mumap, domain_decomp, grid_pointers,
!                  usrdef_grid_params, usrdef_model_grid
!
! Module calls - close_file, construct_regular_grid, error_alloc, open_file
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE grid_params
USE iopars
USE multigrid
USE paralpars
USE physpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: construct_rectgrid
USE inout_routines, ONLY: close_file, open_file

#ifdef CDF
   USE netcdf
#endif /*CDF*/

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j, lev
REAL :: delxdat, delydat, x0dat, y0dat


pglev = 1
procname(1) = 'generate_grids'

!
!1. Initialise parameters
!------------------------
!
!---no parallel mode
parallel_set = .FALSE.; iopt_MPI = 0; master  = .TRUE.; mastermod = .TRUE.

!---netCDF
#ifdef CDF
   iopt_CDF = 1
#else
   iopt_CDF = 0
#endif /*CDF*/

!---log file
exitlog = .TRUE.
loglev1 = 7; loglev2 = MERGE(loglev1,0,exitlog)
runlog_file = 'gridlog'
IF (loglev1.GT.0) THEN
   CALL open_file(iolog,runlog_file,'OUT','A')
ENDIF

!---error file
errchk = .TRUE.; maxerrors = MaxErrMesgs
errlog_file = 'griderr'
CALL open_file(ioerr,errlog_file,'OUT','A')

!---default file names
info_file = ''
post_file = ''
decomp_file = ''

!---default file parameters
modfiles(io_mppmod,1,2)%form = 'A'
modfiles(io_mppmod,1,2)%status = 'W'
modfiles(io_mppmod,1,2)%filename = ' '

!
!2. Define generator setup 
!-------------------------
!
!---user-defined parameters
CALL usrdef_grid_params

!---reset output file names (if needed)
intitle = runtitle
IF (LEN_TRIM(info_file).EQ.0) THEN
   info_file  = TRIM(runtitle)//'.partitI'
ENDIF
IF (LEN_TRIM(decomp_file).EQ.0) THEN
   decomp_file  = TRIM(runtitle)//'.decomp'
ENDIF
IF (LEN_TRIM(post_file).EQ.0) THEN
   post_file  = TRIM(runtitle)
ENDIF

!
!3. Model grid
!-------------
!
!3.1 Allocate
!------------
!

ALLOCATE (gxcoord(nc+1,nr+1),STAT=errstat)
CALL error_alloc('gxcoord',2,(/nc+1,nr+1/),real_type)
gxcoord = 0.0
ALLOCATE (gycoord(nc+1,nr+1),STAT=errstat)
CALL error_alloc('gycoord',2,(/nc+1,nr+1/),real_type)
gycoord = 0.0
ALLOCATE (depmeanatc(nc,nr),STAT=errstat)
CALL error_alloc('depmeanatc',2,(/nc,nr/),real_type)
depmeanatc = 0.0
ALLOCATE (iobu(nobu),STAT=errstat)
CALL error_alloc('iobu',1,(/nobu/),int_type)
iobu = 0
ALLOCATE (jobu(nobu),STAT=errstat)
CALL error_alloc('jobu',1,(/nobu/),int_type)
jobu = 0
ALLOCATE (iobv(nobv),STAT=errstat)
CALL error_alloc('iobv',1,(/nobv/),int_type)
iobv = 0
ALLOCATE (jobv(nobv),STAT=errstat)
CALL error_alloc('jobv',1,(/nobv/),int_type)
jobv = 0
ALLOCATE (nodeatc(0:nc,0:nr),STAT=errstat)
CALL error_alloc('nodeatc',2,(/nc+1,nr+1/),int_type)
nodeatc = 0
ALLOCATE (nodeatu(nc+1,nr,1),STAT=errstat)
CALL error_alloc('nodeatu',3,(/nc+1,nr,1/),int_type)
nodeatu = 0
ALLOCATE (nodeatv(nc,nr+1,1),STAT=errstat)
CALL error_alloc('nodeatv',3,(/nc,nr+1,1/),int_type)
nodeatv = 0

!
!3.2 Define grid
!---------------
!
!---user-defined
CALL usrdef_model_grid

!---uniform rectangular grid
IF (iopt_grid_htype.EQ.1) THEN
   delxdat = surfacegrids(igrd_model,1)%delxdat
   delydat = surfacegrids(igrd_model,1)%delydat
   x0dat = surfacegrids(igrd_model,1)%x0dat
   y0dat = surfacegrids(igrd_model,1)%y0dat
   CALL construct_rectgrid(x0dat,y0dat,delxdat,delydat,gxcoord,gycoord,&
                         & nc+1,nr+1)
!---non-uniform rectangular grid
ELSEIF (iopt_grid_htype.EQ.2) THEN
   j_321: DO j=2,nr
      gxcoord(1:nc,j) = gxcoord(1:nc,1)
   ENDDO j_321
   i_322: DO i=2,nc
      gycoord(i,1:nr) = gycoord(1,1:nr)
   ENDDO i_322
ENDIF

!---extend grid at eastern and northern edges
IF (iopt_grid_htype.GT.1) THEN
   gxcoord(nc+1,1:nr) = 2.0*gxcoord(nc,1:nr) - gxcoord(nc-1,1:nr)
   gxcoord(1:nc,nr+1) = 2.0*gxcoord(1:nc,nr) - gxcoord(1:nc,nr-1)
   gxcoord(nc+1,nr+1) = 2.0*gxcoord(nc,nr+1) - gxcoord(nc-1,nr+1)
   gycoord(nc+1,1:nr) = 2.0*gycoord(nc,1:nr) - gycoord(nc-1,1:nr)
   gycoord(1:nc,nr+1) = 2.0*gycoord(1:nc,nr) - gycoord(1:nc,nr-1)
   gycoord(nc+1,nr+1) = 2.0*gycoord(nc+1,nr) - gycoord(nc+1,nr-1)
ENDIF

!---water depths at edges
depmeanatc(nc,:) = 0.0; depmeanatc(:,nr) = 0.0

!---multigrid dimensions
mgvars(0)%nc = nc; mgvars(0)%nr = nr
lev_322: DO lev=1,nomglevels-1
   mgvars(lev)%nc = (mgvars(lev-1)%nc+1)/2
   mgvars(lev)%nr = (mgvars(lev-1)%nr+1)/2
ENDDO lev_322

!
!3.3 Pointer arrays
!------------------
!

CALL grid_pointers

!
!4. Domain decomposition
!-----------------------
!

CALL domain_decomp

IF (iopt_partit_post.EQ.1) CALL decomp_mumap

!
!5. Deallocate
!-------------
!

DEALLOCATE (gxcoord,gycoord)
DEALLOCATE (depmeanatc)
DEALLOCATE (iobu,iobv,jobu,jobv)
DEALLOCATE (nodeatc,nodeatu,nodeatv)

lev_510: DO lev=0,nomglevels-1
   DEALLOCATE (mgvars(lev)%nc1procs,mgvars(lev)%nc2procs,&
             & mgvars(lev)%nr1procs,mgvars(lev)%nr2procs)
ENDDO lev_510

!
!6. Close files
!--------------
!

CALL close_file(ioerr,'A',fildel=.TRUE.)
IF (loglev1.GT.0) CALL close_file(iolog,'A')


STOP 'Domain generator terminated'


END PROGRAM generate_grids

!========================================================================

SUBROUTINE domain_decomp
!************************************************************************
!
! *domain_decomp* Define domain partition
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Generate_Grids.F90  V2.10
!  
! Description -
!
! Reference -
!
! Calling program - generate_grids
!
! External calls - reset_partition, set_partition, write_partition
!
! Module calls - close_file, convert_loc_to_char, open_file
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE grid_params
USE iopars
USE multigrid
USE paralpars
USE physpars
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: convert_loc_to_char
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=1) :: outform
CHARACTER (LEN=21) :: charloc1, charloc2
INTEGER :: iproc, iunit, i1, i2, j1, j2, maxsize, maxwet, minsize, minwet, &
         & ntotsize, ntotwet
REAL :: meansize, meanwet, nptglbwet, npttotsize, npttotwet, nsizerange, &
      & nwetrange
INTEGER, ALLOCATABLE, DIMENSION(:) :: nsize


procname(pglev+1) = 'domain_decomp'
CALL log_timer_in()

!
!1. Define decomposition
!-----------------------
!
!---regular partition
CALL set_partition

!---reset if necessary
IF (iopt_partit_reset.EQ.1) CALL reset_partition

!
!2. Define domain info
!---------------------
!
!---allocate
ALLOCATE (nowetcprocs(nprocs),STAT=errstat)
CALL error_alloc('nowetcprocs',1,(/nprocs/),real_type)
ALLOCATE (nsize(nprocs),STAT=errstat)
CALL error_alloc('nsize',1,(/nprocs/),real_type)

!---define parameters and arrays
iproc_210: DO iproc=1,nprocs
   i1 = mgvars(0)%nc1procs(iproc); i2 = mgvars(0)%nc2procs(iproc)
   j1 = mgvars(0)%nr1procs(iproc); j2 = mgvars(0)%nr2procs(iproc)
   nsize(iproc) = (i2-i1+1)*(j2-j1+1)
   nowetcprocs(iproc) = COUNT(depmeanatc(i1:i2,j1:j2).GT.0.0)
ENDDO iproc_210

ntotsize = SUM(nsize(1:nprocs))
maxsize = MAXVAL(nsize(1:nprocs))
minsize = MINVAL(nsize(1:nprocs))
meansize = ntotsize/nprocs
nsizerange = 100.0*(maxsize-minsize)/REAL(meansize)
npttotsize = 100.0*(ntotsize/REAL(nc*nr))

ntotwet = SUM(nowetcprocs(1:nprocs))
maxwet = MAXVAL(nowetcprocs(1:nprocs))
minwet = MINVAL(nowetcprocs(1:nprocs))
meanwet = ntotwet/nprocs
nwetrange = 100.0*(maxwet-minwet)/REAL(meanwet)
npttotwet = 100.0*(ntotwet/REAL(ntotsize))
nptglbwet = 100.0*(COUNT(depmeanatc(1:nc,1:nr).GT.0.0)/REAL(nc*nr))

!
!3. Write info
!-------------
!

CALL open_file(iunit,info_file,'OUT','A')
WRITE (iunit,9001) 'Title: ', TRIM(runtitle)
WRITE (iunit,9002) 'npwork', npwork
WRITE (iunit,9002) 'nprocs', nprocs
WRITE (iunit,9002) 'nprocsx', nprocsx
WRITE (iunit,9002) 'nprocsy', nprocsy
WRITE (iunit,9002) 'nomglevels', nomglevels
WRITE (iunit,9002) 'Total size', ntotsize
WRITE (iunit,9002) 'Max size', maxsize
WRITE (iunit,9002) 'Min size', minsize
WRITE (iunit,9003) 'Max size ratio', maxsize/REAL(meansize)
WRITE (iunit,9003) 'Min size ratio', minsize/REAL(meansize)
WRITE (iunit,9003) 'Size range (%)', nsizerange
WRITE (iunit,9003) '% of total domain', npttotsize
WRITE (iunit,9002) 'Max no. of wet points', maxwet
WRITE (iunit,9002) 'Min no. of wet points',minwet
WRITE (iunit,9003) 'Max wet size ratio', maxwet/REAL(meanwet)
WRITE (iunit,9003) 'Min wet size ratio', minwet/REAL(meanwet)
WRITE (iunit,9003) 'Wet size range (%)', nwetrange
WRITE (iunit,9003) '% wet in global domain', nptglbwet
WRITE (iunit,9003) '% wet in new domain', npttotwet

WRITE (iunit,*)
WRITE (iunit,'(A)') 'PID  size   %total %wet    nc1  nc2  nr1  nr2      S&
  &          W        N         E&
  &      '
WRITE (iunit,'(A)') '--------------------------------------------------------&
  &--------------------------------'
iproc_310: DO iproc=1,nprocs
   i1 = mgvars(0)%nc1procs(iproc); i2 = mgvars(0)%nc2procs(iproc)
   j1 = mgvars(0)%nr1procs(iproc); j2 = mgvars(0)%nr2procs(iproc)
   charloc1 = convert_loc_to_char(gxcoord(i1,j1),gycoord(i1,j1))
   charloc2 = convert_loc_to_char(gxcoord(i2,j2),gycoord(i2,j2))
   WRITE (iunit,9004) iproc, nsize(iproc), &
                    & 100.0*(nsize(iproc)/REAL(ntotsize)), &
                    & 100.0*(nowetcprocs(iproc)/REAL(nsize(iproc))), &
                    & i1, i2, j1, j2, charloc1(1:9), charloc1(12:21), &
                    & charloc2(1:9), charloc2(12:21)
ENDDO iproc_310

CALL close_file(iunit,'A')

!
!4. Write 'mppmod'-file
!----------------------
!

CALL write_partition

!
!5. Deallocate
!-------------
!

DEALLOCATE (nowetcprocs,nsize)

CALL log_timer_out()


RETURN

9001 FORMAT (A,': ',T26,A)
9002 FORMAT (A,': ',T26,I10)
9003 FORMAT (A,': ',T26,F6.2)
9004 FORMAT (I4,1X,I6,2(1X,F6.2),4(1X,I4),1X,4(1X,A))

END SUBROUTINE domain_decomp

!========================================================================

SUBROUTINE grid_pointers
!************************************************************************
!
! *grid_pointers* Define grid pointer arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Generate_Grids.F90  V2.9
!  
! Description -
!
! Reference -
!
! Calling program - generate_grids
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, ii, j, jj


procname(pglev+1) = 'grid_pointers'
CALL log_timer_in()

!
!1. Centres
!----------
!

nodeatc(1:nc,1:nr) = MERGE(1,0,depmeanatc(1:nc,1:nr).GT.0)

!
!2. U-nodes
!----------
!
!---interior faces
j_210: DO j=1,nr
i_210: DO i=1,nc
   nodeatu(i,j,1) = COUNT(nodeatc(i-1:i,j).GT.0)
ENDDO i_210
ENDDO j_210

!---open boundaries
ii_220: DO ii=1,nobu
   i = iobu(ii); j = jobu(ii)
   nodeatu(i,j,1) = 3
ENDDO ii_220

!
!3. V-nodes
!----------
!
!---interior faces
j_310: DO j=1,nr
i_310: DO i=1,nc
   nodeatv(i,j,1) = COUNT(nodeatc(i,j-1:j).GT.0)
ENDDO i_310
ENDDO j_310

!---open boundaries
jj_320: DO jj=1,nobv
   i = iobv(jj); j = jobv(jj)
   nodeatv(i,j,1) = 3
ENDDO jj_320

CALL log_timer_out()


RETURN

END SUBROUTINE grid_pointers

!========================================================================

SUBROUTINE set_partition
!************************************************************************
!
! *set_partition* Decomposes the computational domain using simple
!                 partitioning into nprocsx*nprocsy subdomains
!
! Author - Patrick Luyten and Pieter Rauwoens
!
! Version - @(COHERENS)Generate_Grids.F90  V2.10
!
! Description - account is made of sub-grids levels in case of a multigrid
!               procedure
!
! Reference -
!
! Calling program - domain_decomp
!
! External calls - sqrt_factoring
!
!************************************************************************
!
USE gridpars
USE grid_params
USE iopars
USE multigrid
USE paralpars
USE physpars
USE syspars
USE error_routines, ONLY: error_alloc_struc_comp
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, iproc, j, group, lev, ncmod, ncg, nrg, nrmod, nxdiv, nxmod, &
         & nydiv, nymod


procname(pglev+1) = 'set_partition'
CALL log_timer_in()

!
!1. Domain sizes
!----------------
!

IF (nprocsx.EQ.0.AND.nprocsy.EQ.0) THEN
   CALL sqrt_factoring(npwork,nprocsx,nprocsy,nc,nr)
ELSEIF (nprocsx.GT.0.AND.nprocsy.EQ.0) THEN
   nprocsy = npwork/nprocsx
ELSEIF (nprocsx.EQ.0.AND.nprocsy.GT.0) THEN
   nprocsx = npwork/nprocsy
ELSEIF (nprocsx.GT.0.AND.nprocsy.GT.0) THEN
   npwork = nprocsx*nprocsy
ENDIF

!
!2. Domain decomposition
!-----------------------
!

lev_200: DO lev=0,nomglevels-1

!
!2.1 Allocate
!------------
!

   ALLOCATE(mgvars(lev)%nc1procs(npwork),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','nc1procs',1,(/npwork/),int_type,lev)
   mgvars(lev)%nc1procs = 0

   ALLOCATE(mgvars(lev)%nc2procs(npwork),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','nc2procs',1,(/npwork/),int_type,lev)
   mgvars(lev)%nc2procs = 0

   ALLOCATE(mgvars(lev)%nr1procs(npwork),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','nr1procs',1,(/npwork/),int_type,lev)
   mgvars(lev)%nr1procs = 0

   ALLOCATE(mgvars(lev)%nr2procs(npwork),STAT=errstat)
   CALL error_alloc_struc_comp('mgvars','nr2procs',1,(/npwork/),int_type,lev)
   mgvars(lev)%nr2procs = 0

!
!2.2 Initialise parameters
!-------------------------
!

   group = 2**(nomglevels-lev-1)
   ncmod = MOD(mgvars(lev)%nc,group); nrmod = MOD(mgvars(lev)%nr,group)
   IF (ncmod.EQ.0) THEN
      ncg = mgvars(lev)%nc/group; ncmod = group
   ELSE
      ncg = mgvars(lev)%nc/group + 1
   ENDIF
   IF (nrmod.EQ.0) THEN
      nrg = mgvars(lev)%nr/group; nrmod = group
   ELSE
      nrg = mgvars(lev)%nr/group + 1
   ENDIF
   
   nxdiv = ncg/nprocsx*group; nxmod = MOD(ncg,nprocsx)
   nydiv = nrg/nprocsy*group; nymod = MOD(nrg,nprocsy)

!
!2.3 Simple decomposition
!------------------------
!

   j_230: DO j=1,nprocsy
   i_230: DO i=1,nprocsx
      iproc = (j-1)*nprocsx + i
      IF (i.LE.(nprocsx-nxmod)) THEN
         mgvars(lev)%nc1procs(iproc) = 1 + (i-1)*nxdiv
         ncloc = nxdiv
      ELSE
         mgvars(lev)%nc1procs(iproc) = 1 + (nprocsx-nxmod)*nxdiv + &
                                     & (i-(nprocsx-nxmod+1))*(nxdiv+group)
         ncloc = nxdiv + group
      ENDIF
      IF (i.EQ.nprocsx) ncloc = ncloc - group + ncmod
      mgvars(lev)%nc2procs(iproc) = mgvars(lev)%nc1procs(iproc) + ncloc - 1
      IF (j.LE.(nprocsy-nymod)) THEN
         mgvars(lev)%nr1procs(iproc) = 1 + (j-1)*nydiv
         nrloc = nydiv
      ELSE
         mgvars(lev)%nr1procs(iproc) = 1 + (nprocsy-nymod)*nydiv + &
                                     & (j-(nprocsy-nymod+1))*(nydiv+group)
         nrloc = nydiv + group
      ENDIF
      IF (j.EQ.nprocsy) nrloc = nrloc - group + nrmod
      mgvars(lev)%nr2procs(iproc) = mgvars(lev)%nr1procs(iproc) + nrloc - 1
   ENDDO i_230
   ENDDO j_230

ENDDO lev_200

CALL log_timer_out()


RETURN

END SUBROUTINE set_partition

!========================================================================

SUBROUTINE reset_partition
!************************************************************************
!
! *reset_partition* Reset domain partition
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Generate_Grids.F90  V2.10
!
! Description -
!
! Reference -
!
! Calling program - domain_decomp
!
! External calls -
!
! Module calls -
!
!************************************************************************
!
USE grid
USE grid_params
USE iopars
USE multigrid
USE paralpars
USE physpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iiproc, iproc, i1, i2, j1, j2, lev
LOGICAL, DIMENSION(npwork) :: procdel
INTEGER, DIMENSION(npwork) :: nc1tmp, nc2tmp, nr1tmp, nr2tmp 


procname(pglev+1) = 'reset_partition'
CALL log_timer_in()

!
!1. Mark dry domains for deletion
!--------------------------------
!

iproc_100: DO iproc=1,npwork
   i1 = mgvars(0)%nc1procs(iproc); i2 = mgvars(0)%nc2procs(iproc)
   j1 = mgvars(0)%nr1procs(iproc); j2 = mgvars(0)%nr2procs(iproc)
   IF (ANY(nodeatc(i1:i2,j1:j2).EQ.1)) THEN
      procdel(iproc)  = .FALSE.
   ELSEIF (ANY(nodeatu(i1,j1:j2,1).GT.0).OR.&
         & ANY(nodeatu(i2+1,j1:j2,1).GT.0)) THEN
      procdel(iproc)  = .FALSE.
   ELSEIF (ANY(nodeatv(i1:i2,j1,1).GT.0).OR.&
         & ANY(nodeatv(i1:i2,j2+1,1).GT.0)) THEN
      procdel(iproc)  = .FALSE.
   ELSE
      procdel(iproc) = .TRUE.
   ENDIF
ENDDO iproc_100

!
!2. Remove empty domains
!----------------------
!

nprocs = COUNT(.NOT.procdel)

lev_200: DO lev=0,nomglevels-1
   nc1tmp = mgvars(lev)%nc1procs; nc2tmp = mgvars(lev)%nc2procs
   nr1tmp = mgvars(lev)%nr1procs; nr2tmp = mgvars(lev)%nr2procs
   iiproc = 0
   iproc_200: DO iproc=1,npwork
      IF (.NOT.procdel(iproc)) THEN
         iiproc = iiproc + 1
         mgvars(lev)%nc1procs(iiproc) = nc1tmp(iproc)
         mgvars(lev)%nc2procs(iiproc) = nc2tmp(iproc)
         mgvars(lev)%nr1procs(iiproc) = nr1tmp(iproc)
         mgvars(lev)%nr2procs(iiproc) = nr2tmp(iproc)
      ENDIF
   ENDDO iproc_200
ENDDO lev_200

CALL log_timer_out()


RETURN

END SUBROUTINE reset_partition

!========================================================================

SUBROUTINE sqrt_factoring(nps,npx,npy,nx,ny)
!************************************************************************
!
! *sqrt_factoring* Split nps into factors npx and npy
!
! Version - @(COHERENS)Generate_Grids.F90  V2.9
!
! Author - Patrick Luyten
!
! Reference -
!
! Calling program - set_partition
!
!************************************************************************
!
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: nps, nx, ny
INTEGER, INTENT(INOUT) :: npx, npy

!
!*Local variables
!
INTEGER :: i


procname(pglev+1) = 'sqrt_factoring'
CALL log_timer_in()

i = SQRT(REAL(nps))
DO WHILE ((nps/i)*i.NE.nps)
   i = i - 1
ENDDO
IF (nx.LE.ny) THEN
   npx = i; npy = nps/npx
ELSE     
   npy = i; npx = nps/npy
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE sqrt_factoring

!========================================================================

SUBROUTINE decomp_mumap
!************************************************************************
!
! *decomp_mumap* Write domain decomposition to mumap file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Generate_Grids.F90  V2.10
!  
! Description -
!
! Reference -
!
! Calling program - generate_grids
!
! Module calls - close_fil,e open_file
!
!************************************************************************
!
USE iopars
USE grid
USE grid_params
USE multigrid
USE paralpars
USE switches
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iproc, iunit1, iunit2, i1, i2, j1, j2
REAL :: xc, yc


procname(pglev+1) = 'decomp_mumap'
CALL log_timer_in()

!---write MUMAP class 'c' file
CALL open_file(iunit1,TRIM(post_file)//'.cdat','OUT','A')
CALL open_file(iunit2,TRIM(post_file)//'.adat','OUT','A')

iproc_110: DO iproc=1,nprocs
   i1 = mgvars(0)%nc1procs(iproc); i2 = mgvars(0)%nc2procs(iproc) + 1
   j1 = mgvars(0)%nr1procs(iproc); j2 = mgvars(0)%nr2procs(iproc) + 1
   IF (iopt_grid_sph.EQ.0) THEN
      WRITE (iunit1,9001) gxcoord(i1,j1), gycoord(i1,j1), 0
      WRITE (iunit1,9001) gxcoord(i1,j1), gycoord(i2,j2), 1
      WRITE (iunit1,9001) gxcoord(i2,j2), gycoord(i2,j2), 1
      WRITE (iunit1,9001) gxcoord(i2,j2), gycoord(i1,j1), 1
      WRITE (iunit1,9001) gxcoord(i1,j1), gycoord(i1,j1), 1
   ELSE
      WRITE (iunit1,9001) gycoord(i1,j1), gxcoord(i1,j1), 0
      WRITE (iunit1,9001) gycoord(i1,j1), gxcoord(i2,j2), 1
      WRITE (iunit1,9001) gycoord(i2,j2), gxcoord(i2,j2), 1
      WRITE (iunit1,9001) gycoord(i2,j2), gxcoord(i1,j1), 1
      WRITE (iunit1,9001) gycoord(i1,j1), gxcoord(i1,j1), 1
   ENDIF
   xc = 0.5*(gxcoord(i1,j1)+gxcoord(i2,j1))
   yc = 0.5*(gycoord(i1,j1)+gycoord(i1,j2))
   IF (iopt_grid_sph.EQ.0) THEN
      WRITE (iunit2,9002) xc, yc, 2, iproc
   ELSE
      WRITE (iunit2,9002) yc, xc, 2, iproc
   ENDIF
ENDDO iproc_110

CALL close_file(iunit1,'A')
CALL close_file(iunit2,'A')

CALL log_timer_out()


RETURN

9001 FORMAT ((2G15.7,1X),I1)
9002 FORMAT ((2G15.7,1X),I2,1X,I3)

END SUBROUTINE decomp_mumap

!========================================================================

SUBROUTINE write_partition
!************************************************************************
!
! *write_partition* Write domain partition
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Generate_Grids.F90  V2.10
!  
! Description -
!
! Reference -
!
! Calling program - domain_decomposition
!
! Module calls - close_filepars, error_alloc, error_alloc_struc,
!                open_filepars, set_modfiles_atts, set_modvars_atts,
!                write_atts_mod, write_vars
!
!************************************************************************
!
USE datatypes
USE grid
USE iopars
USE multigrid
USE paralpars
USE physpars
USe switches
USE syspars
USE cf90_routines, ONLY: cf90_inq_libvers
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: convert_date, log_timer_in, log_timer_out

!
!* Local variables
!
CHARACTER (LEN=lentime) :: cdatetimex
INTEGER :: iglb, ilev, lev, numvars
INTEGER, DIMENSION(8) :: intdate
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nx1procs, nx2procs, ny1procs, ny2procs
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'write_partition'
CALL log_timer_in()

!
!1. CF global attributes
!-----------------------
!

iglb_710: DO iglb=1,numglbatts
   SELECT CASE (iglb)
      CASE (1)
         glbatts(iglb)%name = 'Conventions'
         glbatts(iglb)%value = CF_version
      CASE (2)
         glbatts(iglb)%name = 'title'
      CASE (3)
         glbatts(iglb)%name = 'history'
         CALL DATE_AND_TIME(VALUES=intdate)
         intdate(4)= intdate(5); intdate(5) = intdate(6)
         intdate(6)= intdate(7); intdate(7) = intdate(8)
         cdatetimex = convert_date(intdate(1:7),seps='- ::')
         IF (TRIM(history_CF).EQ.'') THEN
            history_CF = 'created '//cdatetimex(1:19)
         ENDIF
         glbatts(iglb)%value = ''
      CASE (4)
         glbatts(iglb)%name = 'institution'
         glbatts(iglb)%value = &
                         & 'RBINS – Operational Directorate Natural Environment'
      CASE (5)
         glbatts(iglb)%name = 'source'
         glbatts(iglb)%value = 'Coherens version '//TRIM(model_version)
      CASE (6)
         glbatts(iglb)%name = 'comment'
         glbatts(iglb)%value = ''
      CASE (7)
         glbatts(iglb)%name = 'references'
         glbatts(iglb)%value = ''
      CASE (8)
         IF (iopt_CDF.EQ.1) THEN
            glbatts(iglb)%name = 'netcdf'
            glbatts(iglb)%value = cf90_inq_libvers()
         ENDIF
   END SELECT
ENDDO iglb_710

!
!2. Open data file
!-----------------
!
!---file attributes
CALL set_modfiles_atts(io_mppmod,1,2)
filepars = modfiles(io_mppmod,1,2)
filepars%floattype = 'S'

numvars = filepars%novars

!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL set_modvars_atts(io_mppmod,1,2,varatts,numvars)

!---write
CALL write_atts_mod(filepars,varatts,numvars)

!
!3. Write partition
!------------------
!
!---allocate
ALLOCATE (nx1procs(nprocs,nomglevels),STAT=errstat)
CALL error_alloc('nx1procs',2,(/nprocs,nomglevels/),int_type)
ALLOCATE (nx2procs(nprocs,nomglevels),STAT=errstat)
CALL error_alloc('nx2procs',2,(/nprocs,nomglevels/),int_type)
ALLOCATE (ny1procs(nprocs,nomglevels),STAT=errstat)
CALL error_alloc('ny1procs',2,(/nprocs,nomglevels/),int_type)
ALLOCATE (ny2procs(nprocs,nomglevels),STAT=errstat)
CALL error_alloc('ny2procs',2,(/nprocs,nomglevels/),int_type)

!---store
lev_310: DO lev=0,nomglevels-1
   ilev = lev + 1
   nx1procs(:,ilev) = mgvars(lev)%nc1procs(1:nprocs)
   nx2procs(:,ilev) = mgvars(lev)%nc2procs(1:nprocs)
   ny1procs(:,ilev) = mgvars(lev)%nr1procs(1:nprocs)
   ny2procs(:,ilev) = mgvars(lev)%nr2procs(1:nprocs)
ENDDO lev_310

!---write
CALL write_vars(nx1procs,filepars,1,varatts)
CALL write_vars(nx2procs,filepars,2,varatts)
CALL write_vars(ny1procs,filepars,3,varatts)
CALL write_vars(ny2procs,filepars,4,varatts)

!
!4. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_mppmod,1,2) = filepars

!---deallocate
DEALLOCATE (varatts)
DEALLOCATE (nx1procs,nx2procs,ny1procs,ny2procs)

CALL log_timer_out()


RETURN

END SUBROUTINE write_partition
