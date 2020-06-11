!************************************************************************
!
! *Post_Program* Coherens posprocessing program
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Post_Program.F90  V2.11.2
!
! $Date: 2015-04-23 17:11:33 +0200 (Thu, 23 Apr 2015) $
!
! $Revision: 845 $
!
! Description - 
!
! Reference -
!
! Routines - coords_to_index, HP_plot, HT_plot, PD_plot, PH_plot, post_end,
!            post_read, post_start, PV_plot, read_line, TH_plot, time_limits,
!            time_locate, TS_plot, TZ_plot, vert_retrieve_1d, vert_retrieve_3d,
!            VT_plot, write_mumap_coast, write_mumap_cont, write_mumap_par,
!            write_mumap_vec, write_muplot_par, ZP_plot
!
!************************************************************************
!

!========================================================================

PROGRAM post_program
!************************************************************************
!
! *post_program* Main postprocessing program
!
! Author - Patrick Luyten
!
! Description - prepares files for plotting using idl
!
! Reference -
!
! External calls - HP_plot, HT_plot, PD_plot, PH_plot, post_end, post_start,
!                  PV_plot, read_line, TH_plot, TS_plot, TZ_plot, VT_plot,
!                  ZP_plot
!
! Module calls - close_file, open_file
!
!************************************************************************
!
USE iopars
USE paralpars
USE plotpars
USE switches
USE syspars
USE inout_routines, ONLY: close_file, open_file

#ifdef CDF
   USE netcdf
#endif /*CDF*/

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iser, iunit, l, l1, l2, n, nolevels, notracks
REAL :: test
LOGICAL, DIMENSION(Maxpars) :: lval
CHARACTER (LEN=MaxString), DIMENSION(Maxpars) :: cval
INTEGER, DIMENSION(Maxpars) :: ival
REAL, DIMENSION(Maxpars) :: rval


pglev = 1
procname(1) = 'post_program'

!
!1. Initialise parameters
!------------------------
!
!---no parallel mode
parallel_set = .FALSE.; iopt_MPI = 0; iopt_CDF_abort = 1
master  = .TRUE.

!---kind parameter for real data
kndrtype = KIND(test)

!---data type for real data
float_type = MERGE(real_type,rlong_type,kndrtype.EQ.kndreal)

!---universal parameters
IF (kndrtype.EQ.kndreal) THEN
   pi = pi_s
   halfpi = halfpi_s
   twopi = twopi_s
   enap = enap_s
   degtorad = degtorad_s
   radtodeg = radtodeg_s
ELSE
   pi = pi_d
   halfpi = halfpi_d
   twopi = twopi_d
   enap = enap_d
   degtorad = degtorad_d
   radtodeg = radtodeg_d
ENDIF

!---netcdf data types
#ifdef CDF
   char_NF90 = NF90_char; clobber_NF90 = NF90_clobber
   double_NF90 = NF90_double; fill_NF90 = NF90_fill
   fill_double_NF90 = NF90_fill_double
   fill_int_NF90 = NF90_fill_int; fill_real_NF90 = NF90_fill_real 
   global_NF90 = NF90_global; int_NF90 = NF90_int
   noerr_NF90 = NF90_noerr; nofill_NF90 = NF90_nofill
   nowrite_NF90 = NF90_nowrite; offset_64bit_NF90 = NF90_64bit_offset
   real_NF90 = NF90_real; share_NF90 = NF90_share
   sizehint_default_NF90 = NF90_sizehint_default
   unlimited_NF90 = NF90_unlimited; write_NF90 = NF90_write
#endif /*CDF*/

!---data flags
#ifdef CDF
   int_fill = fill_int_NF90
   real_fill = fill_real_NF90;
   double_fill = fill_double_NF90; 
#else
   int_fill = fill_int_def
   real_fill = fill_real_def
   double_fill = fill_double_def
#endif /*CDF*/
real_fill_eps = 0.00001*ABS(real_fill)
double_fill_eps = 0.00001*ABS(double_fill)
IF (float_type.EQ.real_type) THEN
   float_fill = real_fill
   float_fill_eps = real_fill_eps
ELSE
   float_fill = double_fill
   float_fill_eps = double_fill_eps
ENDIF

!---log file
exitlog = .TRUE.
loglev1 = 7; loglev2 = MERGE(loglev1,0,exitlog)
runlog_file = 'postlog'
IF (loglev1.GT.0) THEN
   CALL open_file(iolog,runlog_file,'OUT','A')
ENDIF

!---error file
errchk = .TRUE.; maxerrors = MaxErrMesgs
errlog_file = 'posterr'
CALL open_file(ioerr,errlog_file,'OUT','A')

!---files.vis file
CALL open_file(io_filevis,'files.vis','OUT','A')

!
!2. Process post parameter files
!-------------------------------
!

OPEN (5,FILE='defposts',STATUS='OLD')
9999 READ (5,'(A)',END=9998) postfile
IF (LEN_TRIM(postfile).EQ.0.OR.postfile(1:1).EQ.'!') GOTO 9999
CALL open_file(iunit,postfile,'IN','A')
plinenum = 0
CALL read_line(iunit,1,(/1/),cval,lval,ival,rval,'runtitle')
runtitle = TRIM(cval(1))


CALL open_file(io_parvis,TRIM(runtitle)//'.vis','OUT','A')
WRITE (io_filevis,'(A)') TRIM(runtitle)//'.vis'

!
!3. Read plot parameter file
!---------------------------
!
!3.1 Coast line
!--------------
!

CALL read_line(iunit,1,(/1/),cval,lval,ival,rval,'coast_name')
coast_name = cval(1)
IF (coast_name(1:1).EQ.'%') THEN
   coast_def = .TRUE.
ELSE
   coast_def = .FALSE.
ENDIF

!
!3.2 Data file
!-------------
!
!---read name/form of data file
999 CALL read_line(iunit,2,(/1,1/),cval,lval,ival,rval,'datafile form')
outfile = TRIM(cval(1)); outform = TRIM(cval(2))
filepars%filename = outfile; filepars%form = outform; filepars%status = 'R'

!---generic file name for plotting
l2 = LEN_TRIM(outfile) - MERGE(3,4,outform.EQ.'N')
l1 = SCAN(outfile(1:l2),'/',back=.TRUE.)
l1 = l1 + 1
plotfile = outfile(l1:l2)
l = LEN_TRIM(runtitle)
plotfile(1:l) = runtitle(1:l)

!---initialise post-processing of datafile
IF (outform.EQ.'N') CALL post_start

!
!3.3 Plot type
!-------------
!

998 CALL read_line(iunit,1,(/1/),cval,lval,ival,rval,'plottype')
plottype = TRIM(cval(1))

!
!4. Read and process data
!------------------------
!

SELECT CASE (TRIM(plottype))

!
!4.1 Horizontal transects
!------------------------
!

   CASE ('HT')

!     ---number of time series
      CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'noseries')
      noseries = ival(1); noplots = 0

      iser_410: DO iser=1,noseries

!        ---(re)initialise post-processing of datafile
         IF (outform.NE.'N') CALL post_start

!        ---time location
         CALL read_line(iunit,4,(/1,1,3,2/),cval,lval,ival,rval,&
                     & 'PStartDateTime PEndDateTime norecskip ztime')
         PStartDateTime = TRIM(cval(1))//',000'
         PEndDateTime = TRIM(cval(2))//',000'
         norecskip = ival(1)
         ztime = MERGE(.FALSE.,lval(1),.NOT.time_grid)
         ptime_format = 0
!        ---number of lists
         CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'nolists')
         nolists = ival(1)
!        ---allocate
         ALLOCATE (postatts(nolists))

         n_411: DO n=1,nolists

!           ---animation
            CALL read_line(iunit,1,(/2/),cval,lval,ival,rval,'anim')
            postatts(n)%anim = lval(1)
!           ---horizontal location
            CALL read_line(iunit,6,(/4,4,4,4,3,4/),cval,lval,ival,rval,&
                        & 'gxstart gystart gxend gyend horz_unit hlen_unit')
            postatts(n)%gxstart = rval(1); postatts(n)%gystart = rval(2)
            postatts(n)%gxend = rval(3); postatts(n)%gyend = rval(4)
            postatts(n)%horz_unit = ival(1); postatts(n)%hlen_unit = rval(5)
            IF (iopt_grid_sph.EQ.1) postatts(n)%hlen_unit = 1.0
!           ---scalar and vector index
            CALL read_line(iunit,3,(/3,3,3/),cval,lval,ival,rval,&
                        & 'numvar numvec1 numvec2')
            postatts(n)%numvar = ival(1)
            postatts(n)%numvec1 = ival(2); postatts(n)%numvec2 = ival(3)
!           ---parameters for contouring
            IF (postatts(n)%numvar.GT.0) THEN
               CALL read_line(iunit,2,(/3,3/),cval,lval,ival,rval,&
                           & 'novarskip1 novarskip2')
               postatts(n)%novarskip1 = ival(1)
               postatts(n)%novarskip2 = ival(2)
               CALL read_line(iunit,3,(/2,3,3/),cval,lval,ival,rval,&
                           & 'contreg nolevels icontstyle')
               postatts(n)%contreg = lval(1); postatts(n)%nolevels = ival(1)
               nolevels = ival(1)
               postatts(n)%icontstyle = ival(2)
               IF (postatts(n)%contreg) THEN
                  CALL read_line(iunit,2,(/4,4/),cval,lval,ival,rval,&
                              & 'contmin contmax')
                  postatts(n)%contmin = rval(1); postatts(n)%contmax = rval(2)
               ELSE
                  CALL read_line(iunit,nolevels,(/(4,l=1,nolevels)/),&
                               & cval, lval,ival,rval,'contlevels')
                  postatts(n)%contlevels(1:nolevels) = rval(1:nolevels)
               ENDIF
            ENDIF
!           ---parameters for vector plot
            IF (postatts(n)%numvec1.GT.0) THEN
               CALL read_line(iunit,3,(/3,3,4/),cval,lval,ival,rval,&
                           & 'novecskip1 novecskip2 hrefvec')
               postatts(n)%novecskip1 = ival(1)
               postatts(n)%novecskip2 = ival(2)
               postatts(n)%hrefvec = rval(1)
            ENDIF
!           ---vertical plot specifications
            IF (nodim.EQ.3) THEN
               CALL read_line(iunit,3,(/3,3,4/),cval,lval,ival,rval,&
                           & 'iplotform kpos depplot')
               postatts(n)%iplotform = ival(1); postatts(n)%kpos = ival(2)
               postatts(n)%depplot = rval(1)
            ENDIF

         ENDDO n_411

!        ---prepare plot
         CALL HT_plot

!        ---deallocate
         DEALLOCATE (postatts)

      ENDDO iser_410

      GOTO 998

!
!4.2 Vertical transects
!----------------------
!

   CASE ('VT')

!     ---number of time series
      CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'noseries')
      noseries = ival(1); noplots = 0

      iser_420: DO iser=1,noseries

!        ---(re)initialise post-processing of datafile
         IF (outform.NE.'N') CALL post_start

!        ---time location
         CALL read_line(iunit,4,(/1,1,3,2/),cval,lval,ival,rval,&
                     & 'PStartDateTime PEndDateTime norecskip ztime')
         PStartDateTime = TRIM(cval(1))//',000'
         PEndDateTime = TRIM(cval(2))//',000'
         norecskip = ival(1)
         ztime = MERGE(.FALSE.,lval(1),.NOT.time_grid)
         ptime_format = 0
!        ---number of lists
         CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'nolists')
         nolists = ival(1)
!        ---allocate
         ALLOCATE (postatts(nolists))

         n_421: DO n=1,nolists

!           ---animation
            CALL read_line(iunit,1,(/2/),cval,lval,ival,rval,'anim')
            postatts(n)%anim = lval(1)
!           ---horizontal location
            CALL read_line(iunit,6,(/4,4,4,4,3,4/),cval,lval,ival,rval,&
                        & 'gxstart gystart gxend gyend horz_unit hlen_unit')
            postatts(n)%gxstart = rval(1); postatts(n)%gystart = rval(2)
            postatts(n)%gxend = rval(3); postatts(n)%gyend = rval(4)
            postatts(n)%horz_unit = ival(1); postatts(n)%hlen_unit = rval(5)
            IF (iopt_grid_sph.EQ.1) postatts(n)%hlen_unit = 1.0
!           ---scalar and vector index
            CALL read_line(iunit,4,(/3,3,3,3/),cval,lval,ival,rval,&
                        & 'numvar numvec1 numvec2 numvec3')
            postatts(n)%numvar = ival(1)
            postatts(n)%numvec1 = ival(2)
            postatts(n)%numvec2 = ival(3)
            postatts(n)%numvec3 = ival(4)
!           ---parameters for contouring
            IF (postatts(n)%numvar.GT.0) THEN
               CALL read_line(iunit,2,(/3,3/),cval,lval,ival,rval,&
                           & 'novarskip1 novarskip2')
               postatts(n)%novarskip1 = ival(1)
               postatts(n)%novarskip2 = ival(2)
               CALL read_line(iunit,3,(/2,3,3/),cval,lval,ival,rval,&
                           & 'contreg nolevels icontstyle')
               postatts(n)%contreg = lval(1); postatts(n)%nolevels = ival(1)
               nolevels = ival(1); postatts(n)%icontstyle = ival(2)
               IF (postatts(n)%contreg) THEN
                  CALL read_line(iunit,2,(/4,4/),cval,lval,ival,rval,&
                              & 'contmin contmax')
                  postatts(n)%contmin = rval(1); postatts(n)%contmax = rval(2)
               ELSE
                  CALL read_line(iunit,nolevels,(/(4,l=1,nolevels)/),&
                               & cval,lval,ival,rval,'contlevels')
                  postatts(n)%contlevels(1:nolevels) = rval(1:nolevels)
               ENDIF
            ENDIF
!           ---parameters for vector plot
            IF (postatts(n)%numvec1.GT.0) THEN
               CALL read_line(iunit,4,(/3,3,4,4/),cval,lval,ival,rval,&
                           & 'novecskip1 novecskip2 hrefvec vrefvec')
               postatts(n)%novecskip1 = ival(1)
               postatts(n)%novecskip2 = ival(2)
               postatts(n)%hrefvec = rval(1); postatts(n)%vrefvec = rval(2)
            ENDIF
!           ---vertical location
            CALL read_line(iunit,3,(/4,4,3/),cval,lval,ival,rval,&
                        & 'depmin depmax vert_unit')
            postatts(n)%depmin = rval(1); postatts(n)%depmax = rval(2)
            postatts(n)%vert_unit = ival(1)
            
         ENDDO n_421

!        ---prepare plot
         CALL VT_plot

!        ---deallocate
         DEALLOCATE (postatts)

      ENDDO iser_420

      GOTO 998

!
!4.3 Time series
!---------------
!

   CASE ('TS')

!     ---number of time series
      CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'noseries')
      noseries = ival(1); noplots = 0

      iser_430: DO iser=1,noseries

!        ---(re)initialise post-processing of datafile
         IF (outform.NE.'N') CALL post_start

!        ---time location
         CALL read_line(iunit,6,(/1,1,1,3,3,2/),cval,lval,ival,rval,&
     & 'PStartDateTime PEndDateTime PRefDateTime norecskip ptime_format ztime')
         PStartDateTime = TRIM(cval(1))//',000'
         PEndDateTime = TRIM(cval(2))//',000'
         PRefDateTime = TRIM(cval(3))//',000'
         norecskip = ival(1); ptime_format = ival(2)
         ztime = MERGE(.FALSE.,lval(1),.NOT.time_grid)
!        ---number of lists
         CALL read_line(iunit,2,(/3,3/),cval,lval,ival,rval,'nolists nofigs')
         nolists = ival(1); nofigs = ival(2)
!        ---allocate
         ALLOCATE (postatts(nolists))

         n_431: DO n=1,nolists
!           ---horizontal location
            IF (nodim.GT.0) THEN
               IF (gridded) THEN
                  CALL read_line(iunit,3,(/4,4,3/),cval,lval,ival,rval,&
                              & 'gxpos gypos horz_unit')
                  postatts(n)%gxpos = rval(1); postatts(n)%gypos = rval(2)
                  postatts(n)%horz_unit = ival(1)
               ELSE
                  CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'numstat')
                  postatts(n)%numstat = ival(1)
               ENDIF
            ENDIF

!           ---variable index
            CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'numvar')
            postatts(n)%numvar = ival(1)
!           ---vertical plot specifications
            IF (nodim.EQ.3) THEN
               CALL read_line(iunit,3,(/3,3,4/),cval,lval,ival,rval,&
                           & 'iplotform kpos depplot')
               postatts(n)%iplotform = ival(1)
               postatts(n)%kpos = ival(2)
               postatts(n)%depplot = rval(1)
            ENDIF
!           ---plot style
            CALL read_line(iunit,2,(/3,3/),cval,lval,ival,rval,&
                        & 'numfig linepsyms')
            postatts(n)%numfig = ival(1); postatts(n)%linepsyms = ival(2)

         ENDDO n_431

!        ---prepare plot
         CALL TS_plot

!        ---deallocate
         DEALLOCATE (postatts)

      ENDDO iser_430

      GOTO 998

!
!4.4 Depth-time contours
!-----------------------
!

   CASE ('TZ')

!     ---number of time series
      CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'noseries')
      noseries = ival(1); noplots = 0

      iser_440: DO iser=1,noseries

!        ---(re)initialise post-processing of datafile
         IF (outform.NE.'N') CALL post_start

!        ---time location
         CALL read_line(iunit,6,(/1,1,1,3,3,2/),cval,lval,ival,rval,&
     & 'PStartDateTime PEndDateTime PRefDateTime norecskip ptime_format ztime')
         PStartDateTime = TRIM(cval(1))//',000'
         PEndDateTime = TRIM(cval(2))//',000'
         PRefDateTime = TRIM(cval(3))//',000'
         norecskip = ival(1); ptime_format = ival(2)
         ztime = MERGE(.FALSE.,lval(1),.NOT.time_grid)
!        ---number of lists
         CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'nolists')
         nolists = ival(1)
!        ---allocate
         ALLOCATE (postatts(nolists))

         n_441: DO n=1,nolists
!           ---horizontal location
            IF (gridded) THEN
               CALL read_line(iunit,3,(/4,4,3/),cval,lval,ival,rval,&
                           & 'gxpos gypos horz_unit')
               postatts(n)%gxpos = rval(1); postatts(n)%gypos = rval(2)
               postatts(n)%horz_unit = ival(1)
            ELSE
               CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'numstat')
               postatts(n)%numstat = ival(1)
            ENDIF
!           ---variable index
            CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'numvar')
            postatts(n)%numvar = ival(1)
!           ---parameters for contouring
            CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'novarskip1')
            postatts(n)%novarskip1 = ival(1)
            CALL read_line(iunit,3,(/2,3,3/),cval,lval,ival,rval,&
                        & 'contreg nolevels icontstyle')
            postatts(n)%contreg = lval(1); postatts(n)%nolevels = ival(1)
            nolevels = ival(1); postatts(n)%icontstyle = ival(2)
            IF (postatts(n)%contreg) THEN
               CALL read_line(iunit,2,(/4,4/),cval,lval,ival,rval,&
                           & 'contmin contmax')
               postatts(n)%contmin = rval(1); postatts(n)%contmax = rval(2)
            ELSE
               CALL read_line(iunit,nolevels,(/(4,l=1,nolevels)/),&
                            & cval,lval,ival,rval,'contlevels')
               postatts(n)%contlevels(1:nolevels) = rval(1:nolevels)
            ENDIF
!           ---vertical location
            CALL read_line(iunit,3,(/4,4,3/),cval,lval,ival,rval,&
                        & 'depmin depmax vert_unit')
            postatts(n)%depmin = rval(1); postatts(n)%depmax = rval(2)
            postatts(n)%vert_unit = ival(1)

         ENDDO n_441

!        ---prepare plot
         CALL TZ_plot

!        ---deallocate
         DEALLOCATE (postatts)

      ENDDO iser_440
      
      GOTO 998

!
!4.5 Transect-time contours
!--------------------------
!

   CASE ('TH')

!     ---number of time series
      CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'noseries')
      noseries = ival(1); noplots = 0

      iser_450: DO iser=1,noseries

!        ---(re)initialise post-processing of datafile
         IF (outform.NE.'N') CALL post_start

!        ---time location
         CALL read_line(iunit,6,(/1,1,1,3,3,2/),cval,lval,ival,rval,&
     & 'PStartDateTime PEndDateTime PRefDateTime norecskip ptime_format ztime')
         PStartDateTime = TRIM(cval(1))//',000'
         PEndDateTime = TRIM(cval(2))//',000'
         PRefDateTime = TRIM(cval(3))//',000'
         norecskip = ival(1); ptime_format = ival(2)
         ztime = MERGE(.FALSE.,lval(1),.NOT.time_grid)

!        ---number of lists
         CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'nolists')
         nolists = ival(1)
!        ---allocate
         ALLOCATE (postatts(nolists))

         n_451: DO n=1,nolists
!           ---horizontal location
            CALL read_line(iunit,6,(/4,4,4,4,3,4/),cval,lval,ival,rval,&
                        & 'gxstart gystart gxend gyend horz_unit hlen_unit')
            postatts(n)%gxstart = rval(1); postatts(n)%gystart = rval(2)
            postatts(n)%gxend = rval(3); postatts(n)%gyend = rval(4)
            postatts(n)%horz_unit = ival(1); postatts(n)%hlen_unit = rval(5)
            IF (iopt_grid_sph.EQ.1) postatts(n)%hlen_unit = 1.0
!           ---variable index
            CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'numvar')
            postatts(n)%numvar = ival(1)
!           ---parameters for contouring
            CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'novarskip1')
            postatts(n)%novarskip1 = ival(1)
            CALL read_line(iunit,3,(/2,3,3/),cval,lval,ival,rval,&
                        & 'contreg nolevels icontstyle')
            postatts(n)%contreg = lval(1); postatts(n)%nolevels = ival(1)
            nolevels = ival(1); postatts(n)%icontstyle = ival(2)
            IF (postatts(n)%contreg) THEN
               CALL read_line(iunit,2,(/4,4/),cval,lval,ival,rval,&
                           & 'contmin contmax')
               postatts(n)%contmin = rval(1); postatts(n)%contmax = rval(2)
            ELSE
               CALL read_line(iunit,nolevels,(/(4,l=1,nolevels)/),&
                            & cval,lval,ival,rval,'contlevels')
               postatts(n)%contlevels(1:nolevels) = rval(1:nolevels)
            ENDIF
!           ---vertical plot specifications
            IF (nodim.EQ.3) THEN
               CALL read_line(iunit,3,(/3,3,4/),cval,lval,ival,rval,&
                           & 'iplotform kpos depplot')
               postatts(n)%iplotform = ival(1); postatts(n)%kpos = ival(2)
               postatts(n)%depplot = rval(1)
            ENDIF
         ENDDO n_451

!        ---prepare plot
         CALL TH_plot

!        ---deallocate
         DEALLOCATE (postatts)

      ENDDO iser_450

      GOTO 998

!
!4.6 Vertical profiles
!---------------------
!

   CASE ('ZP')

!    ---(re)initialise post-processing of datafile
      IF (outform.NE.'N') CALL post_start

!    ---time grid
      CALL read_line(iunit,1,(/2/),cval,lval,ival,rval,'ztime')
      ztime = MERGE(.FALSE.,lval(1),.NOT.time_grid)
      ptime_format = 0

!     ---number of lists
      CALL read_line(iunit,2,(/3,3/),cval,lval,ival,rval,'nolists nofigs')
      nolists = ival(1); nofigs = ival(2); noplots = 0

!     ---allocate
      ALLOCATE (postatts(nolists))

      n_460: DO n=1,nolists

!        ---time location
         CALL read_line(iunit,1,(/1/),cval,lval,ival,rval,'PDateTime')
         postatts(n)%PDateTime = TRIM(cval(1))//',000'
!        ---horizontal location
         IF (gridded) THEN
            CALL read_line(iunit,3,(/4,4,3/),cval,lval,ival,rval,&
                        & 'gxpos gypos horz_unit')
            postatts(n)%gxpos = rval(1); postatts(n)%gypos = rval(2)
            postatts(n)%horz_unit = ival(1)
         ELSE
            CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'numstat')
            postatts(n)%numstat = ival(1)
         ENDIF
!        ---vertical location
         CALL read_line(iunit,3,(/4,4,3/),cval,lval,ival,rval,&
                     & 'depmin depmax vert_unit')
         postatts(n)%depmin = rval(1); postatts(n)%depmax = rval(2)
         postatts(n)%vert_unit = ival(1)

!        ---variable index
         CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'numvar')
         postatts(n)%numvar = ival(1)
!        ---plot style
         CALL read_line(iunit,2,(/3,3/),cval,lval,ival,rval,'numfig linepsyms')
         postatts(n)%numfig = ival(1); postatts(n)%linepsyms = ival(2)

      ENDDO n_460

!     ---prepare plot
      CALL ZP_plot

!     ---deallocate
      DEALLOCATE (postatts)

      GOTO 998

!
!4.7 Horizontal profiles
!-----------------------
!

   CASE ('HP')

!    ---(re)initialise post-processing of datafile
      IF (outform.NE.'N') CALL post_start

!    ---time grid
      CALL read_line(iunit,1,(/2/),cval,lval,ival,rval,'ztime')
      ztime = MERGE(.FALSE.,lval(1),.NOT.time_grid)
      ptime_format = 0

!     ---number of lists
      CALL read_line(iunit,2,(/3,3/),cval,lval,ival,rval,'nolists nofigs')
      nolists = ival(1); nofigs = ival(2); noplots = 0

!     ---allocate
      ALLOCATE (postatts(nolists))

      n_470: DO n=1,nolists

!        ---time location
         CALL read_line(iunit,1,(/1/),cval,lval,ival,rval,'PDateTime')
         postatts(n)%PDateTime = TRIM(cval(1))//',000'
!        ---horizontal location
         CALL read_line(iunit,6,(/4,4,4,4,3,4/),cval,lval,ival,rval,&
                     & 'gxstart gystart gxend gyend horz_unit hlen_unit')
         postatts(n)%gxstart = rval(1); postatts(n)%gystart = rval(2)
         postatts(n)%gxend = rval(3); postatts(n)%gyend = rval(4)
         postatts(n)%horz_unit = ival(1); postatts(n)%hlen_unit = rval(5)
         IF (iopt_grid_sph.EQ.1) postatts(n)%hlen_unit = 1.0
!        ---variable index
         CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'numvar')
         postatts(n)%numvar = ival(1)
!        ---vertical plot specifications
         IF (nodim.EQ.3) THEN
            CALL read_line(iunit,3,(/3,3,4/),cval,lval,ival,rval,&
                        & 'iplotform kpos depplot')
            postatts(n)%iplotform = ival(1); postatts(n)%kpos = ival(2)
            postatts(n)%depplot = rval(1)
         ENDIF
!        ---plot style
         CALL read_line(iunit,2,(/3,3/),cval,lval,ival,rval,&
                     & 'numfig linepsyms')
         postatts(n)%numfig = ival(1); postatts(n)%linepsyms = ival(2)

      ENDDO n_470

!     ---prepare plot
      CALL HP_plot

!     ---deallocate
      DEALLOCATE (postatts)

      GOTO 998

!
!4.8 Horizontal particle trajectories
!------------------------------------     
!

   CASE ('PH')

!     ---number of time series
      CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'noseries')
      noseries = ival(1); noplots = 0

      iser_480: DO iser=1,noseries

!        ---(re)initialise post-processing of datafile
         IF (outform.NE.'N') CALL post_start

!        ---time location
         CALL read_line(iunit,4,(/1,1,3,3/),cval,lval,ival,rval,&
                     & 'PStartDateTime PEndDateTime norecskip ptime_format')
         PStartDateTime = TRIM(cval(1))//',000'
         PEndDateTime = TRIM(cval(2))//',000'
         norecskip = ival(1); ptime_format = ival(2)
!        ---number of lists
         CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'nolists')
         nolists = ival(1)
!        ---allocate
         ALLOCATE (postatts(nolists))
         
         n_481: DO n=1,nolists

!           ---horizontal location
            CALL read_line(iunit,6,(/4,4,4,4,3,4/),cval,lval,ival,rval,&
                        & 'gxstart gystart gxend gyend horz_unit hlen_unit')
            postatts(n)%gxstart = rval(1); postatts(n)%gystart = rval(2)
            postatts(n)%gxend = rval(3); postatts(n)%gyend = rval(4)
            postatts(n)%horz_unit = ival(1); postatts(n)%hlen_unit = rval(5)
            IF (iopt_grid_sph.EQ.1) postatts(n)%hlen_unit = 1.0
!           ---number of trajectories
            CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'notracks')
            notracks = ival(1)
            IF (notracks.EQ.0) THEN
               notracks = nopartout
            ELSE
               notracks = ival(1)
            ENDIF
            postatts(n)%notracks = notracks
!           ---particle indices
            IF (ival(1).EQ.0) THEN
               postatts(n)%ptrack(1:notracks) = (/(l,l=1,nopartout)/)
            ELSE
               CALL read_line(iunit,notracks,(/(3,l=1,notracks)/),cval,lval,&
                            & ival,rval,'ptrack')
               postatts(n)%ptrack(1:notracks) = ival(1:notracks)
            ENDIF
            CALL read_line(iunit,notracks,(/(3,l=1,notracks)/),cval,lval,&
                         & ival,rval,'linestyles')
            postatts(n)%linestyles(1:notracks) = ival(1:notracks)
         ENDDO n_481

!        ---prepare plot         
         CALL PH_plot

!        ---deallocate
         DEALLOCATE (postatts)

      ENDDO iser_480

      GOTO 998

!
!4.9 Particle trajectories along vertical transects
!--------------------------------------------------
!      

   CASE ('PV')

!     ---number of time series
      CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'noseries')
      noseries = ival(1); noplots = 0

      iser_490: DO iser=1,noseries

!        ---(re)initialise post-processing of datafile
         IF (outform.NE.'N') CALL post_start

!        ---time location
         CALL read_line(iunit,4,(/1,1,3,3/),cval,lval,ival,rval,&
                     & 'PStartDateTime PEndDateTime norecskip ptime_format')
         PStartDateTime = TRIM(cval(1))//',000'
         PEndDateTime = TRIM(cval(2))//',000'
         norecskip = ival(1); ptime_format = ival(2)
!        ---number of lists
         CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'nolists')
         nolists = ival(1)
!        ---allocate
         ALLOCATE (postatts(nolists))

         n_491: DO n=1,nolists
!           ---horizontal location
            CALL read_line(iunit,6,(/4,4,4,4,3,4/),cval,lval,ival,rval,&
                        & 'gxstart gystart gxend gyend horz_unit hlen_unit')
            postatts(n)%gxstart = rval(1); postatts(n)%gystart = rval(2)
            postatts(n)%gxend = rval(3); postatts(n)%gyend = rval(4)
            postatts(n)%horz_unit = ival(1); postatts(n)%hlen_unit = rval(5)
            IF (iopt_grid_sph.EQ.1) postatts(n)%hlen_unit = 1.0
!           ---number of trajectories
            CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'notracks')
            notracks = ival(1)
            IF (notracks.EQ.0) THEN
               notracks = nopartout
            ELSE
               notracks = ival(1)
            ENDIF
            postatts(n)%notracks = notracks
!           ---particle indices
            IF (ival(1).EQ.0) THEN
               postatts(n)%ptrack(1:notracks) = (/(l,l=1,nopartout)/)
            ELSE
               CALL read_line(iunit,notracks,(/(3,l=1,notracks)/),cval,lval,&
                            & ival,rval,'ptrack')
               postatts(n)%ptrack(1:notracks) = ival(1:notracks)
            ENDIF
            CALL read_line(iunit,notracks,(/(3,l=1,notracks)/),cval,lval,&
                         & ival,rval,'linestyles')
            postatts(n)%linestyles(1:notracks) = ival(1:notracks)

         ENDDO n_491

!        ---prepare plot         
         CALL PV_plot

!        ---deallocate
         DEALLOCATE (postatts)

      ENDDO iser_490

      GOTO 998

!
!4.10 Particle propersties as function of travelled distance
!-----------------------------------------------------------
!      

   CASE ('PD')

!     ---number of time series
      CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'noseries')
      noseries = ival(1); noplots = 0

      iser_4100: DO iser=1,noseries

!        ---(re)initialise post-processing of datafile
         IF (outform.NE.'N') CALL post_start

!        ---time location
         CALL read_line(iunit,4,(/1,1,3,3/),cval,lval,ival,rval,&
                     & 'PStartDateTime PEndDateTime norecskip ptime_format')
         PStartDateTime = TRIM(cval(1))//',000'
         PEndDateTime = TRIM(cval(2))//',000'
         norecskip = ival(1); ptime_format = ival(2)
!        ---number of lists
         CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'nolists')
         nolists = ival(1)
!        ---allocate
         ALLOCATE (postatts(nolists))

         n_4101: DO n=1,nolists
!           ---length unit
            CALL read_line(iunit,1,(/4/),cval,lval,ival,rval,'hlen_unit')
            postatts(n)%hlen_unit = rval(1)
            IF (iopt_grid_sph.EQ.1) postatts(n)%hlen_unit = 1.0
!           ---variable index
            CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'numvar')
            postatts(n)%numvar = ival(1)
!           ---number of trajectories
            CALL read_line(iunit,1,(/3/),cval,lval,ival,rval,'notracks')
            notracks = ival(1)
            IF (notracks.EQ.0) THEN
               notracks = nopartout
            ELSE
               notracks = ival(1)
            ENDIF
            postatts(n)%notracks = notracks
!           ---particle indices
            IF (ival(1).EQ.0) THEN
               postatts(n)%ptrack(1:notracks) = (/(l,l=1,nopartout)/)
            ELSE
               CALL read_line(iunit,notracks,(/(3,l=1,notracks)/),cval,lval,&
                            & ival,rval,'ptrack')
               postatts(n)%ptrack(1:notracks) = ival(1:notracks)
            ENDIF
            CALL read_line(iunit,notracks,(/(3,l=1,notracks)/),cval,lval,&
                         & ival,rval,'linestyles')
            postatts(n)%linestyles(1:notracks) = ival(1:notracks)
            
         ENDDO n_4101

!        ---prepare plot         
         CALL PD_plot

!        ---deallocate
         DEALLOCATE (postatts)

      ENDDO iser_4100

      GOTO 998

!
!4.11 New file
!-------------
!

   CASE ('N')

      CALL post_end
      GOTO 999

!
!4.12 Quit
!---------
!

   CASE ('Q')

      CALL post_end
      GOTO 997

END SELECT

!
!5. Close files
!--------------
!

997 CALL close_file(iunit,'A')
CALL close_file(io_parvis,'A')

!
!6. New parameter file
!---------------------
!

GOTO 9999

9998 CONTINUE
WRITE (io_filevis,'(A)') 'end'
CALL close_file(io_filevis,'A')
CALL close_file(ioerr,'A',fildel=.TRUE.)
IF (loglev1.GT.0) CALL close_file(iolog,'A')


STOP 'Postprocessor terminated'

END PROGRAM post_program

!========================================================================

SUBROUTINE coords_to_index(x,y,iloc,jloc)
!************************************************************************
!
! *coords_to_index* Convert geographical into index coordinates
!
! Author - Patrick Luyten
!
! Description - 
!
! Reference -
!
! Calling program - HP_plot, HT_plot, TH_plot, TS_plot, TZ_plot, VT_plot,
!                   write_mumap_ocats, ZP_plot
!
! Module calls - distance_minloc, error_abort, error_limits_var
!
!************************************************************************
!
USE grid
USE iopars
USE plotpars
USE error_routines, ONLY: error_abort, error_limits_var
USE grid_routines, ONLY: distance_minloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(OUT) :: iloc, jloc
REAL, INTENT(INOUT) :: x, y

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*x*         REAL    Geographical coordinate in X-direction
!*y*         REAL    Geographical coordinate in Y-direction
!*iloc*      INTEGER Index coordinate in X-direction
!*jloc*      INTEGER Index coordinate in Y-direction
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, DIMENSION(2) :: indloc


procname(pglev+1) = 'coords_to_index'
CALL log_timer_in()

!
!1. Index locations
!------------------
!

CALL distance_minloc(x,y,gxcoord,gycoord,ncout,nrout,2,3,locmin=indloc)
iloc = indloc(1); jloc = indloc(2)
x = gxcoord(iloc,jloc); y = gycoord(iloc,jloc)

!
!2. Check
!--------
!

CALL error_limits_var(iloc,'iloc',1,ncout)
CALL error_limits_var(jloc,'jloc',1,nrout)
CALL error_abort(procname(pglev),ierrno_runval)

CALL log_timer_out()


RETURN

END SUBROUTINE coords_to_index

!========================================================================

SUBROUTINE HP_plot
!************************************************************************
!
! *HP_plot* Horizontal profiles along transects
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - post_program
!
! External calls - coords_to_index, post_read, time_locate,
!                  vert_retrieve_1d, write_muplot_par
!
! Module calls - close_file, convert_loc_to_char, error_alloc, open_file
!
!************************************************************************
!
USE grid
USE iopars
USE physpars
USE plotpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: convert_loc_to_char
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=12) :: ccurv, cfig
CHARACTER (LEN=leniofile) :: paramfile
CHARACTER (LEN=lentime) :: ctime
CHARACTER (LEN=lendesc) :: cloc
CHARACTER (LEN=12), DIMENSION(4) :: cpos
CHARACTER (LEN=21), DIMENSION(2) :: geoloc
INTEGER :: horz_unit, i, icurv, iend, ifig, ilist, irec, istart, iunit, ivar, &
         & j, jend, jstart, n, nopoints, numcurves
REAL :: dattrans, deli, delj, dels, fac, gxend, gxstart, gyend, gystart
REAL (KIND=kndrlong) :: rtime

CHARACTER (LEN=leniofile), DIMENSION(nolists) :: datfiles
INTEGER, DIMENSION(2) :: irange, istyle
INTEGER, DIMENSION(0:nolists) :: norec
INTEGER, DIMENSION(nolists) :: iodat, itrtype
INTEGER, DIMENSION(nofigs) :: nocurves
INTEGER, DIMENSION(nofigs,nolists) :: indcurv
REAL, DIMENSION(2) :: rangemax, rangemin

INTEGER, ALLOCATABLE, DIMENSION(:) :: indlist, itrans, jtrans
REAL, ALLOCATABLE, DIMENSION(:) :: dtrans, xtrans, ytrans


procname(pglev+1) = 'HP_plot'
CALL log_timer_in()

!
!1. Number of curves and curve index array
!-----------------------------------------
!

ifig_110: DO ifig=1,nofigs
   numcurves = 0
   ilist_111: DO ilist=1,nolists
      IF (postatts(ilist)%numfig.EQ.ifig) THEN
         numcurves = numcurves + 1
         indcurv(ifig,numcurves) = ilist
      ENDIF
   ENDDO ilist_111
   nocurves(ifig) = numcurves
ENDDO ifig_110

!
!2. Time location
!----------------
!

norec(0) = 0
ilist_210: DO ilist=1,nolists
   CALL time_locate(postatts(ilist)%PDateTime,norec(ilist))
ENDDO ilist_210

!
!3. Create/open data files
!-------------------------
!

ifig_310: DO ifig=1,nofigs
icurv_310: DO icurv=1,nocurves(ifig)
   ilist = indcurv(ifig,icurv)
   WRITE (cfig,'(I12)') noplots+ifig; cfig = ADJUSTL(cfig)
   WRITE (ccurv,'(I12)') icurv; ccurv = ADJUSTL(ccurv)
   datfiles(ilist) = TRIM(plotfile)//'.'//TRIM(ccurv)//'HPdat'//TRIM(cfig)
   CALL open_file(iodat(ilist),datfiles(ilist),'OUT','A')
ENDDO icurv_310
ENDDO ifig_310

!
!4. Write data files
!-------------------
!

ilist_410: DO ilist=1,nolists

!
!4.1 Read data
!-------------
!
!  ---read data
   irec_411: DO irec=norec(ilist-1)+1,norec(ilist)
      IF (outform.NE.'N'.OR.irec.EQ.norec(ilist)) THEN
         CALL post_read(irec,ctime,rtime,.TRUE.)
      ENDIF
   ENDDO irec_411

!
!4.2 Define transect
!-------------------
!
!  ---start/end locations
   IF (postatts(ilist)%horz_unit.EQ.1) THEN
      istart = NINT(postatts(ilist)%gxstart)
      iend = NINT(postatts(ilist)%gxend)
      jstart = NINT(postatts(ilist)%gystart)
      jend = NINT(postatts(ilist)%gyend)
      postatts(ilist)%gxstart = gxcoord(istart,jstart)
      postatts(ilist)%gxend = gxcoord(iend,jend)
      postatts(ilist)%gystart = gycoord(istart,jstart)
      postatts(ilist)%gyend = gycoord(iend,jend)
      postatts(ilist)%istart = istart
      postatts(ilist)%iend = iend
      postatts(ilist)%jstart = jstart
      postatts(ilist)%jend = jend
   ELSE
      CALL coords_to_index(postatts(ilist)%gxstart,postatts(ilist)%gystart,&
                         & istart,jstart)
      CALL coords_to_index(postatts(ilist)%gxend,postatts(ilist)%gyend,&
                         & iend,jend)
      postatts(ilist)%istart = istart
      postatts(ilist)%iend = iend
      postatts(ilist)%jstart = jstart
      postatts(ilist)%jend = jend
   ENDIF

!  ---number of locations along transect
   nopoints = MAX(ABS(iend-istart),ABS(jend-jstart)) + 1

!  ---allocate arrays
   ALLOCATE (dtrans(nopoints),STAT=errstat)
   CALL error_alloc('dtrans',1,(/nopoints/),kndrtype)
   ALLOCATE (itrans(nopoints),STAT=errstat)
   CALL error_alloc('itrans',1,(/nopoints/),kndint)
   ALLOCATE (jtrans(nopoints),STAT=errstat)
   CALL error_alloc('jtrans',1,(/nopoints/),kndint)
   ALLOCATE (xtrans(nopoints),STAT=errstat)
   CALL error_alloc('xtrans',1,(/nopoints/),kndrtype)
   ALLOCATE (ytrans(nopoints),STAT=errstat)
   CALL error_alloc('ytrans',1,(/nopoints/),kndrtype)

!  ---type of transect
   IF (jstart.EQ.jend) THEN
      itrtype(ilist) = 1
   ELSEIF (istart.EQ.iend) THEN
      itrtype(ilist) = 2
   ELSE
      itrtype(ilist) = 3
   ENDIF

!  ---coordinate arrays along transect
   SELECT CASE (itrtype(ilist))
   CASE (1)
      itrans = (/(i,i=istart,iend)/)
      jtrans = jstart
      dtrans = (/(gxcoord(i,jstart),i=istart,iend)/)
   CASE (2)
      itrans = istart
      jtrans = (/(j,j=jstart,jend)/)
      dtrans = (/(gycoord(istart,j),j=jstart,jend)/)
   CASE (3)
      deli = (iend-istart)/REAL(nopoints-1)
      delj = (jend-jstart)/REAL(nopoints-1)
      itrans = (/(NINT(istart+(n-1)*deli),n=1,nopoints)/)
      jtrans = (/(NINT(jstart+(n-1)*delj),n=1,nopoints)/)
      xtrans = (/(gxcoord(itrans(n),jtrans(n)),n=1,nopoints)/)
      ytrans = (/(gycoord(itrans(n),jtrans(n)),n=1,nopoints)/)
      dtrans(1) = 0.0
      n_420: DO n=2,nopoints
         IF (iopt_grid_sph.EQ.0) THEN
            dels = SQRT((xtrans(n)-xtrans(n-1))**2+(ytrans(n)-ytrans(n-1))**2)
         ELSEIF (iopt_grid_sph.EQ.1) THEN
            fac = COS(0.5*degtorad*(ytrans(n-1)+ytrans(n)))
            dels = Rearth*degtorad*SQRT((fac*(xtrans(n)-xtrans(n-1)))**2+&
                                      & (ytrans(n)-ytrans(n-1))**2)
         ENDIF
         dtrans(n) = dtrans(n-1) + dels
      ENDDO n_420
   END SELECT
   IF ((iopt_grid_sph.EQ.0).OR.(itrtype(ilist).EQ.3)) THEN
      dtrans= postatts(ilist)%hlen_unit*dtrans
   ENDIF
   
!
!4.3 Store and write data
!------------------------
!

   ivar = postatts(ilist)%numvar
   n_430: DO n=1,nopoints
      i = itrans(n); j = jtrans(n)
      IF (nodim.EQ.2) THEN
         dattrans = outvals3d(i,j,ivar)
      ELSEIF (nodim.EQ.3) THEN
         CALL vert_retrieve_1d(dattrans,outvals4d(i,j,:,ivar),i,j,ilist)
      ENDIF
      IF (ABS(dattrans-fill_value).GT.fill_value_eps) THEN
         WRITE (iodat(ilist),9001) dtrans(n), dattrans
      ENDIF
   ENDDO n_430
   CALL close_file(iodat(ilist),'A')

!
!4.4 Deallocate
!--------------
!

   DEALLOCATE (dtrans,itrans,jtrans,xtrans,ytrans)

ENDDO ilist_410

!
!5. Write parameter file
!-----------------------
!

ifig_510: DO ifig=1,nofigs

!  ---select curves
   numcurves = nocurves(ifig)
   ALLOCATE (indlist(numcurves),STAT=errstat)
   CALL error_alloc('indlist',1,(/numcurves/),kndint)
   indlist = indcurv(ifig,1:numcurves)
   ilist = indlist(1)

!  ---initialise parameters
   horz_unit = postatts(ilist)%horz_unit
   gxstart = postatts(ilist)%gxstart; gxend = postatts(ilist)%gxend
   gystart = postatts(ilist)%gystart; gyend = postatts(ilist)%gyend
   istart = postatts(ilist)%istart; gxend = postatts(ilist)%iend
   jstart = postatts(ilist)%jstart; gyend = postatts(ilist)%jend

!  ---open data file
   WRITE (cfig,'(I12)') noplots + ifig; cfig = ADJUSTL(cfig)
   paramfile = TRIM(plotfile)//'.HPpar'//TRIM(cfig)
   CALL open_file(iunit,paramfile,'OUT','A')

!  ---axis titles
   SELECT CASE (itrtype(ilist))
   CASE (1)
      plottitles(1) = TRIM(coordvars(2)%f90_name)
   CASE (2)
      plottitles(1) = TRIM(coordvars(3)%f90_name)
   CASE (3)
      plottitles(1) = 'Distance along transect'
   END SELECT
   ivar = postatts(ilist)%numvar
   plottitles(2) = TRIM(outvars(ivar)%long_name)//' ('//&
                 & TRIM(outvars(ivar)%units)//')'

!  ---main title
   plottitles(3) = TRIM(outfile)

!  ---subtitle
   SELECT CASE (itrtype(ilist))
   CASE (1)
      IF (iopt_grid_sph.EQ.0.OR.horz_unit.EQ.1) THEN
         WRITE (cpos(1),'(I12)') jstart; cpos(1) = ADJUSTL(cpos(1))
         cloc = 'Transect at j = '//TRIM(cpos(1))
      ELSEIF (iopt_grid_sph.EQ.1) THEN
         geoloc(1) = convert_loc_to_char(gxstart,gystart)
         cloc = 'Transect at '//geoloc(1)(1:9)
      ENDIF
   CASE (2)
      IF (iopt_grid_sph.EQ.0.OR.horz_unit.EQ.1) THEN
         WRITE (cpos(1),'(I12)') istart; cpos(1) = ADJUSTL(cpos(1))
         cloc = 'Transect at i = '//TRIM(cpos(1))
      ELSEIF (iopt_grid_sph.EQ.1) THEN
         geoloc(1) = convert_loc_to_char(gxstart,gystart)
         cloc = 'Transect at '//geoloc(1)(12:21)
      ENDIF
   CASE (3)
      IF (iopt_grid_sph.EQ.0.OR.horz_unit.EQ.1) THEN
         WRITE (cpos(1),'(I12)') istart; cpos(1) = ADJUSTL(cpos(1))
         WRITE (cpos(2),'(I12)') jstart; cpos(2) = ADJUSTL(cpos(2))
         WRITE (cpos(3),'(I12)') iend; cpos(3) = ADJUSTL(cpos(3))
         WRITE (cpos(4),'(I12)') jend; cpos(4) = ADJUSTL(cpos(4))
         cloc = 'From ('//TRIM(cpos(1))//','//TRIM(cpos(2))//' to ('//&
              & TRIM(cpos(3))//','//TRIM(cpos(4))//')'
      ELSEIF (iopt_grid_sph.EQ.1) THEN
         geoloc(1) = convert_loc_to_char(gxstart,gystart)
         geoloc(2) = convert_loc_to_char(gxend,gyend)
         cloc = 'From '//TRIM(geoloc(1))//' to '//TRIM(geoloc(2))
      ENDIF
   END SELECT
   plottitles(4) = TRIM(cloc)//'; '//postatts(ilist)%PDateTime(1:19)

!  ---axis styles and ranges
   istyle(1) = 0; irange(1) = 0
!  IF (nocurves(ifig).EQ.1) THEN
      istyle(2) = 0; irange(2) = 0
!  ELSE
!     istyle(2) = 2; irange(2) = 1
!  ENDIF
   rangemin = 0.0; rangemax = 0.0

!  ---write parameter file
   CALL write_muplot_par(iunit,numcurves,datfiles(indlist),&
                       & istyle,irange,rangemin,rangemax,&
                       & postatts(indlist)%linepsyms)

!  ---close parameter file
   CALL close_file(iunit,'A')

!  ---deallocate
   DEALLOCATE (indlist)

!  ---write to '.vis'-file
   WRITE (io_parvis,9002) 'muplot', TRIM(paramfile), 'y', &
                        & TRIM(outvars(ivar)%f90_name)


ENDDO ifig_510

noplots = noplots + nofigs

CALL log_timer_out()


RETURN

9001 FORMAT (2(G15.7,1X))
9002 FORMAT (4(A,1X))

END SUBROUTINE HP_plot

!========================================================================

SUBROUTINE HT_plot
!************************************************************************
!
! *HT_plot* Horizontal transects
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - post_program
!
! External calls - coords_to_index, post_read, time_limits, vert_retrieve_3d,
!                  write_mumap_coast, write_mumap_cont, write_mumap_par,
!                  write_mumap_vec
!
! Module calls - close_file, error_alloc, lim_dims, loop_index, open_file
!
!************************************************************************
!
USE depths
USE grid
USE iopars
USE plotpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims, loop_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: anim, flag, header
CHARACTER (LEN=1) :: vecclass
CHARACTER (LEN=12) :: cfig, cfil, cnt, cplot
CHARACTER (LEN=lendesc) :: cform
CHARACTER (LEN=lentime) :: ctime
CHARACTER (LEN=leniofile) :: linefile
INTEGER :: iend, ifil, ilist, ireg, istart, iunit, ivar, ivec, ivec1, &
         & ivec2, jend, jstart, l, n, ncvar, ncvec, nofiles, nolevels, &
         & norec, notimes, novarskip1, novarskip2, novecskip1, novecskip2, &
         & nrvar, nrvec, ntout
REAL :: aspect, contmax, contmin, gxend, gxstart, gyend, gystart, hlen_unit, &
      & valmax, valmin, xborder, yborder
REAL (KIND=kndrlong) :: rtime

INTEGER, DIMENSION(3) :: tlims, xlims, ylims
INTEGER, DIMENSION(3,nolists) :: xlimsvar, xlimsvec, ylimsvar, ylimsvec
REAL, DIMENSION(MaxContLevs) :: contlevels

CHARACTER (LEN=lentime), ALLOCATABLE, DIMENSION(:) :: outtime
CHARACTER (LEN=leniofile), ALLOCATABLE, DIMENSION(:) :: paramfiles
CHARACTER (LEN=leniofile), ALLOCATABLE, DIMENSION(:,:) :: contfiles, vecfiles
INTEGER, ALLOCATABLE, DIMENSION(:) :: iopar
REAL, ALLOCATABLE, DIMENSION(:,:) :: contdat, depvar, depvec, xvar, xvec, &
                                   & yvar, yvec
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: vecdat, zvar, zvec


procname(pglev+1) = 'HT_plot'
CALL log_timer_in()

!
!1.Initialise parameters and arrays
!----------------------------
!
!---output time limiters
CALL time_limits(tlims)
!---number of plot time steps
notimes = lim_dims(tlims)
!---allocate plot time series array
ALLOCATE (outtime(notimes),STAT=errstat) 
CALL error_alloc('outtime',1,(/notimes/),kndchar,lenstr=lentime)
outtime = ''
ALLOCATE (contfiles(nolists,notimes),STAT=errstat)
CALL error_alloc('contfiles',2,(/nolists,notimes/),kndchar,lenstr=leniofile)
contfiles = ''
ALLOCATE (vecfiles(nolists,notimes),STAT=errstat)
CALL error_alloc('vecfiles',2,(/nolists,notimes/),kndchar,lenstr=leniofile)
vecfiles = ''
ALLOCATE (paramfiles(notimes),STAT=errstat)
CALL error_alloc('paramfiles',1,(/notimes/),kndchar,lenstr=leniofile)
paramfiles = ''
ALLOCATE (iopar(notimes),STAT=errstat)
CALL error_alloc('iopar',1,(/notimes/),kndint)

!
!2. Horizontal locations
!-----------------------
!

ilist_200: DO ilist=1,nolists

   IF (postatts(ilist)%horz_unit.EQ.1) THEN
      istart = NINT(postatts(ilist)%gxstart)
      iend = NINT(postatts(ilist)%gxend)
      jstart = NINT(postatts(ilist)%gystart)
      jend = NINT(postatts(ilist)%gyend)
      postatts(ilist)%gxstart = gxcoord(istart,jstart)
      postatts(ilist)%gxend = gxcoord(iend,jend)
      postatts(ilist)%gystart = gycoord(istart,jstart)
      postatts(ilist)%gyend = gycoord(iend,jend)
      postatts(ilist)%istart = istart
      postatts(ilist)%iend = iend
      postatts(ilist)%jstart = jstart
      postatts(ilist)%jend = jend
   ELSE
      CALL coords_to_index(postatts(ilist)%gxstart,postatts(ilist)%gystart,&
                         & istart,jstart)
      CALL coords_to_index(postatts(ilist)%gxend,postatts(ilist)%gyend,&
                         & iend,jend)
      postatts(ilist)%istart = istart
      postatts(ilist)%iend = iend
      postatts(ilist)%jstart = jstart
      postatts(ilist)%jend = jend
   ENDIF

   IF (postatts(ilist)%numvar.NE.0) THEN
      novarskip1 = postatts(ilist)%novarskip1
      novarskip2 = postatts(ilist)%novarskip2
      xlimsvar(:,ilist) = (/istart,iend,novarskip1/)
      ylimsvar(:,ilist) = (/jstart,jend,novarskip2/)
   ENDIF
   IF (postatts(ilist)%numvec1.NE.0) THEN
      novecskip1 = postatts(ilist)%novecskip1
      novecskip2 = postatts(ilist)%novecskip2
      xlimsvec(:,ilist) = (/istart,iend,novecskip1/)
      ylimsvec(:,ilist) = (/jstart,jend,novecskip2/)
   ENDIF
ENDDO ilist_200

!
!3. Time loop
!------------
!

ntout = 0
norec_300: DO norec=tlims(1),tlims(2)

!
!3.1 Read data
!-------------
!

   flag = loop_index(tlims,norec)
   IF (flag.OR.outform.NE.'N') THEN
      CALL post_read(norec,ctime,rtime,.TRUE.)
   ENDIF
   IF (.NOT.flag) CYCLE norec_300
   ntout = ntout + 1
   outtime(ntout) = ctime
   WRITE (cnt,'(I12)') ntout; cnt = ADJUSTL(cnt)

!
!3.2 Process data
!----------------
!

   ilist_320: DO ilist=1,nolists

!
!3.2.1 Initialise parameters
!---------------------------
!

      anim = postatts(ilist)%anim
      gxstart = postatts(ilist)%gxstart
      gxend = postatts(ilist)%gxend
      gystart = postatts(ilist)%gystart
      gyend = postatts(ilist)%gyend
      ivar = postatts(ilist)%numvar
      ivec1 = postatts(ilist)%numvec1
      ivec2 = postatts(ilist)%numvec2
      hlen_unit = postatts(ilist)%hlen_unit
      nolevels = postatts(ilist)%nolevels
      contmin = postatts(ilist)%contmin
      contmax = postatts(ilist)%contmax
      vecclass = MERGE('h','g',iopt_grid_htype.EQ.1)
      WRITE (cfig,'(I12)') noplots+ilist; cfig = ADJUSTL(cfig)

!
!3.2.2 Allocate arrays
!---------------------
!
!     ---scalar data
      IF (ivar.NE.0.AND.(norec.EQ.tlims(1).OR.nolists.GT.1)) THEN
         ncvar = lim_dims(xlimsvar(:,ilist))
         nrvar = lim_dims(ylimsvar(:,ilist))
         ALLOCATE (xvar(ncvar,nrvar),STAT=errstat)
         CALL error_alloc('xvar',2,(/ncvar,nrvar/),kndrtype)
         ALLOCATE (yvar(ncvar,nrvar),STAT=errstat)
         CALL error_alloc('yvar',2,(/ncvar,nrvar/),kndrtype)
         ALLOCATE (zvar(ncvar,nrvar,nzout),STAT=errstat)
         CALL error_alloc('zvar',3,(/ncvar,nrvar,nzout/),kndrtype)
         ALLOCATE (depvar(ncvar,nrvar),STAT=errstat)
         CALL error_alloc('depvar',2,(/ncvar,nrvar/),kndrtype)
         ALLOCATE (contdat(ncvar,nrvar),STAT=errstat)
         CALL error_alloc('contdat',2,(/ncvar,nrvar/),kndrtype)
      ENDIF

!     ---vector data
      IF (ivec1.NE.0.AND.(norec.EQ.tlims(1).OR.nolists.GT.1)) THEN
         ncvec = lim_dims(xlimsvec(:,ilist))
         nrvec = lim_dims(ylimsvec(:,ilist))
         ALLOCATE (xvec(ncvec,nrvec),STAT=errstat)
         CALL error_alloc('xvec',2,(/ncvec,nrvec/),kndrtype)
         ALLOCATE (yvec(ncvec,nrvec),STAT=errstat)
         CALL error_alloc('yvec',2,(/ncvec,nrvec/),kndrtype)
         ALLOCATE (zvec(ncvec,nrvec,nzout),STAT=errstat)
         CALL error_alloc('zvec',3,(/ncvec,nrvec,nzout/),kndrtype)
         ALLOCATE (depvec(ncvec,nrvec),STAT=errstat)
         CALL error_alloc('depvec',2,(/ncvec,nrvec/),kndrtype)
         ALLOCATE (vecdat(ncvec,nrvec,2),STAT=errstat)
         CALL error_alloc('vecdat',3,(/ncvec,nrvec,2/),kndrtype)
      ENDIF

!
!3.2.3 Grid for plotting
!-----------------------
!
!     ---scalar grid
      IF (ivar.NE.0) THEN
         xlims = xlimsvar(:,ilist); ylims = ylimsvar(:,ilist)
         xvar = hlen_unit*gxcoord(xlims(1):xlims(2):xlims(3),&
                                & ylims(1):ylims(2):ylims(3))
         yvar = hlen_unit*gycoord(xlims(1):xlims(2):xlims(3),&
                                & ylims(1):ylims(2):ylims(3))
         zvar = gzcoord(xlims(1):xlims(2):xlims(3),&
                      & ylims(1):ylims(2):ylims(3),:)
         depvar = depout(xlims(1):xlims(2):xlims(3),&
                       & ylims(1):ylims(2):ylims(3))
      ENDIF

!     ---vector grid
      IF (ivec1.NE.0) THEN
         xlims = xlimsvec(:,ilist); ylims = ylimsvec(:,ilist)
         xvec = hlen_unit*gxcoord(xlims(1):xlims(2):xlims(3),&
                                & ylims(1):ylims(2):ylims(3))
         yvec = hlen_unit*gycoord(xlims(1):xlims(2):xlims(3),&
                                & ylims(1):ylims(2):ylims(3))
         zvec = gzcoord(xlims(1):xlims(2):xlims(3),&
                      & ylims(1):ylims(2):ylims(3),:)
         depvec = depout(xlims(1):xlims(2):xlims(3),&
                           & ylims(1):ylims(2):ylims(3))
      ENDIF
      
!
!3.2.4 Data for contouring
!-------------------------
!

      IF (ivar.NE.0) THEN

!        ---evaluate data
         xlims = xlimsvar(:,ilist); ylims = ylimsvar(:,ilist)
         IF (nodim.EQ.2) THEN
            contdat = outvals3d(xlims(1):xlims(2):xlims(3),&
                              & ylims(1):ylims(2):ylims(3),ivar)
         ELSEIF (nodim.EQ.3) THEN
            CALL vert_retrieve_3d(contdat,outvals4d(:,:,:,ivar),&
                               & (/ncvar,nrvar/),xlims,ylims,ilist)
         ENDIF

!        ---contour levels
         IF (postatts(ilist)%contreg) THEN
            IF (contmin.EQ.0.0.AND.contmax.EQ.0.0) THEN
               valmin = MINVAL(contdat,&
                             & MASK=ABS(contdat-fill_value).GT.fill_value_eps)
               valmax = MAXVAL(contdat,&
                             & MASK=ABS(contdat-fill_value).GT.fill_value_eps)
            ELSE
               valmin = contmin; valmax = contmax
            ENDIF
            contlevels(1:nolevels) = &
                 &(/(valmin+(n-1)*(valmax-valmin)/(nolevels-1),n=1,nolevels)/)
         ELSE
            contlevels(1:nolevels) = postatts(ilist)%contlevels(1:nolevels)
         ENDIF

!        ---create contour file
         contfiles(ilist,ntout) = TRIM(plotfile)//'.'//TRIM(cnt)//'HTidat'&
                              & //TRIM(cfig)
         CALL open_file(iunit,contfiles(ilist,ntout),'OUT','A')

!        ---write data
         ireg = MERGE(0,1,iopt_grid_htype.EQ.1)
         CALL write_mumap_cont(iunit,ncvar,nrvar,xvar,yvar,contdat,&
              & nolevels,contlevels(1:nolevels),&
              & postatts(ilist)%icontstyle,iopt_grid_sph,ireg)
         CALL close_file(iunit,'A')

      ELSE
         contfiles(ilist,ntout) = ''
      ENDIF

!
!3.2.5 Data for vector plots
!---------------------------
!

      IF (ivec1.NE.0) THEN

!        ---evaluate data
         xlims = xlimsvec(:,ilist); ylims = ylimsvec(:,ilist)
         l_325: DO l=1,2
            ivec = MERGE(ivec1,ivec2,l.EQ.1)
            IF (nodim.EQ.2) THEN
               vecdat(:,:,l) = outvals3d(xlims(1):xlims(2):xlims(3),&
                                       & ylims(1):ylims(2):ylims(3),&
                                       & ivec)
            ELSEIF (nodim.EQ.3) THEN
               CALL vert_retrieve_3d(vecdat(:,:,l),outvals4d(:,:,:,ivec),&
                                  & (/ncvec,nrvec/),xlims,ylims,ilist)
            ENDIF

         ENDDO l_325

!        ---create vector file
         vecfiles(ilist,ntout) = TRIM(plotfile)//'.'//TRIM(cnt)//'HT'//&
                               & vecclass//'dat'//TRIM(cfig)
         CALL open_file(iunit,vecfiles(ilist,ntout),'OUT','A')

!        ---write data
         CALL write_mumap_vec(iunit,vecclass,ncvec,nrvec,xvec,yvec,&
                            & vecdat,postatts(ilist)%hrefvec,&
                            & postatts(ilist)%hrefvec,iopt_grid_sph)
         CALL close_file(iunit,'A')

      ELSE
         vecfiles(ilist,ntout) = ''
      ENDIF

!
!3.2.6 Coast line data
!---------------------
!

      IF (norec.EQ.tlims(2)) THEN
         
!        ---create coast line file
         IF (coast_def) THEN
            linefile = TRIM(plotfile)//'.HTcdat'//TRIM(cfig)
            CALL open_file(iunit,linefile,'OUT','A')
         ELSE
            linefile = TRIM(coast_name)
         ENDIF

!        ---write data
         IF (coast_def) THEN
            CALL write_mumap_coast(iunit,ilist,linefile)
         ENDIF

      ENDIF

!
!3.2.7 Deallocate
!----------------
!

      IF (norec.EQ.tlims(2).OR.nolists.GT.1) THEN
         IF (ivar.NE.0) DEALLOCATE(xvar,yvar,zvar,depvar,contdat)
         IF (ivec1.NE.0) DEALLOCATE(xvec,yvec,zvec,depvec,vecdat)
      ENDIF

!
!3.2.8 Write parameter file(s)
!----------------------------
!

      IF (norec.EQ.tlims(2)) THEN

!
!3.2.8.1 Create/open file(s)
!---------------------------
!

         nofiles = MERGE(1,notimes,anim)
         IF (anim) THEN
            paramfiles(1) = TRIM(plotfile)//'.HTpar'//TRIM(cfig)
            CALL open_file(iopar(1),paramfiles(1),'OUT','A')
         ELSE
             ifil_3281: DO ifil=1,nofiles
               WRITE (cfil,'(I12)') ifil; cfil = ADJUSTL(cfil)
                paramfiles(ifil) = TRIM(plotfile)//'.'//TRIM(cfil)//'HTpar'//&
                                 & TRIM(cfig)
                CALL open_file(iopar(ifil),paramfiles(ifil),'OUT','A')
             ENDDO ifil_3281
         ENDIF

!
!3.2.8.2 Write parameter file(s)
!-------------------------------
!
!        ---plot borders
         xborder = border*(gxend-gxstart)
         yborder = border*(gyend-gystart)
         plotcorners(1) = gxstart - xborder
         plotcorners(2) = gystart - yborder
         plotcorners(3) = gxend + xborder
         plotcorners(4) = gyend + yborder
         plotcorners = hlen_unit*plotcorners
         aspect = 1.0

!        ---plot titles
         plottitles(1) = TRIM(outfile)
         IF (ivar.NE.0) THEN
            plottitles(2) = TRIM(outvars(ivar)%long_name)//' ('//&
                          & TRIM(outvars(ivar)%units)//')'
         ELSE
            plottitles(2) = TRIM(outvars(ivec1)%vector_name)//' ('//&
                          & TRIM(outvars(ivec1)%units)//')'
         ENDIF

!        ---write parameter file(s)
         n_32821: DO n=1,notimes
            ifil = MERGE(1,n,anim)
            header = (.NOT.anim.OR.n.EQ.1)
            IF (nodim.EQ.2) THEN
               plottitles(3) = outtime(n)(1:19)
            ELSEIF (nodim.EQ.3) THEN
               SELECT CASE (postatts(ilist)%iplotform)
               CASE(1)
                  WRITE (cplot,'(I12)') postatts(ilist)%kpos
                  cplot = ADJUSTL(cplot)
                  cform = 'k = '//TRIM(cplot)
               CASE(2)
                  WRITE (cplot,'(F8.2)') postatts(ilist)%depplot
                  cplot = ADJUSTL(cplot)
                  cform = 'depth = '//TRIM(cplot)
               CASE(3); cform = 'depth-averaged'
               CASE(4); cform = 'depth-integrated'
               CASE(5); cform = 'max. value'
               END SELECT
               plottitles(3) = outtime(n)(1:19)//'; '//TRIM(cform)
            ENDIF
            CALL write_mumap_par(iopar(ifil),header,iopt_grid_sph,aspect,3,&
                              & (/'i',vecclass,'c'/),&
                              & (/contfiles(ilist,n),vecfiles(ilist,n),&
                              & linefile/))
            IF (anim.AND.n.LT.notimes) WRITE (iopar(ifil),'(A)') 'z'
         ENDDO n_32821

!        ---close parameter file(s)
         ifil_32822: DO ifil=1,nofiles
            CALL close_file(iopar(ifil),'A')
         ENDDO ifil_32822

!
!3.2.8.3 Write to '.vis'-file
!----------------------------
!

         ifil_3283: DO ifil=1,nofiles
            IF (ivar.NE.0) THEN
               WRITE (io_parvis,9001) 'mumap2', TRIM(paramfiles(ifil)), 'y',&
                                              & TRIM(outvars(ivar)%f90_name)
            ELSEIF (ivec1.NE.0) THEN
               WRITE (io_parvis,9001) 'mumap2', TRIM(paramfiles(ifil)), 'y',&
                                             & TRIM(outvars(ivec1)%vector_name)
            ENDIF
         ENDDO ifil_3283

      ENDIF

   ENDDO ilist_320

ENDDO norec_300

!
!4. Update plot number
!---------------------
!

noplots = noplots + nolists

!
!5. Deallocate time series arrays
!--------------------------------
!

DEALLOCATE (outtime,contfiles,vecfiles,paramfiles,iopar)

CALL log_timer_out()


RETURN

9001 FORMAT(4(A,1X))

END SUBROUTINE HT_plot

!========================================================================

SUBROUTINE PD_plot
!************************************************************************
!
! *PD_plot* Particle attributes versus travelled distance
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - post_program
!
! External calls - post_read, time_limits, write_muplot_par
!
! Module calls - close_file, distance_pts, error_alloc, loop_index, open_file
!
!************************************************************************
!
USE iopars
USE plotpars
USE syspars  
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: distance_pts
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: loop_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=12) :: cfig, cn
CHARACTER (LEN=lentime) :: ctime
CHARACTER (LEN=leniofile) :: paramfile
INTEGER :: ilist, iunit, n, norec, notracks, ntrackmax, numvar, p
REAL :: dist, drange, hlen_unit, pdat, x, y
REAL (KIND=kndrlong) :: rtime
CHARACTER (LEN=leniofile), ALLOCATABLE, DIMENSION(:,:) :: trackfiles
INTEGER, DIMENSION(2) :: istyle, irange
REAL, DIMENSION(2) :: rangemax, rangemin
INTEGER, DIMENSION(3) :: tlims
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: first
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iodat
REAL, DIMENSION(nolists) :: ymax, ymin
REAL, ALLOCATABLE, DIMENSION(:,:) :: dtrack, xtrack_old, ytrack_old

procname(pglev+1) = 'PD_plot'
CALL log_timer_in()



!
!1.Initialise parameters
!-----------------------
!

CALL time_limits(tlims)

!
!2. Allocate arrays
!------------------
!

ntrackmax = MAXVAL(postatts%notracks)
ALLOCATE (trackfiles(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('trackfiles',2,(/nolists,ntrackmax/),kndchar,lenstr=leniofile)
trackfiles = ''
ALLOCATE (iodat(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('iodat',2,(/nolists,ntrackmax/),kndint)
ALLOCATE (xtrack_old(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('xtrack_old',2,(/nolists,ntrackmax/),kndrtype)
ALLOCATE (ytrack_old(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('ytrack_old',2,(/nolists,ntrackmax/),kndrtype)
ALLOCATE (dtrack(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('dtrack',2,(/nolists,ntrackmax/),kndrtype)
dtrack = double_fill
ALLOCATE (first(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('first',2,(/nolists,ntrackmax/),kndlog)
first = .TRUE.

!
!3. Time loop
!------------
!

norec_300: DO norec=tlims(1),tlims(2)

!
!3.1 Read data
!-------------
!

   flag = loop_index(tlims,norec)
   IF (flag.OR.outform.NE.'N') THEN
      CALL post_read(norec,ctime,rtime,.TRUE.)
   ENDIF
   IF (.NOT.flag) CYCLE norec_300
   
!
!3.2 Process data
!----------------
!

   ilist_320: DO ilist=1,nolists

!
!3.2.1 Initialise parameters
!---------------------------
!

      notracks = postatts(ilist)%notracks
      numvar = postatts(ilist)%numvar
      hlen_unit = postatts(ilist)%hlen_unit
      notracks = postatts(ilist)%notracks
      
!
!3.2.2 Create trajectory files
!-----------------------------
!      

      IF (norec.EQ.tlims(1)) THEN
         WRITE (cfig,'(I12)') noplots+ilist; cfig = ADJUSTL(cfig)
         n_322: DO n=1,notracks
            WRITE (cn,'(I12)') n; cn = ADJUSTL(cn)
            trackfiles(ilist,n) = TRIM(plotfile)//'.'//TRIM(cn)//'PDcdat'//&
                                & TRIM(cfig)
            CALL open_file(iodat(ilist,n),trackfiles(ilist,n),'OUT','A')
         ENDDO n_322
      ENDIF
            
!      
!3.2.3 Write trajectory data
!---------------------------
!      

      p_323: DO n=1,notracks
         p = postatts(ilist)%ptrack(n)
         x = outvals2d(p,1)
         y = outvals2d(p,2)
         pdat = outvals2d(p,numvar)
         flag = ABS(pdat-fill_value).GT.fill_value_eps
         IF (flag) THEN
            IF (first(ilist,n)) THEN
               dtrack(ilist,n) = 0.0
               WRITE (iodat(ilist,n),9001) 0.0, pdat
               xtrack_old(ilist,n) = x; ytrack_old(ilist,n) = y
               ymin(ilist) = pdat; ymax(ilist) = pdat
               first(ilist,n) = .FALSE.
            ELSE
               dist = distance_pts(xtrack_old(ilist,n),x,ytrack_old(ilist,n),&
                                 & y,2,3)
               dtrack(ilist,n) = dtrack(ilist,n) + hlen_unit*dist
               xtrack_old(ilist,n) = x; ytrack_old(ilist,n) = y
               WRITE (iodat(ilist,n),9001) dtrack(ilist,n), pdat
               ymin(ilist) = MIN(pdat,ymin(ilist))
               ymax(ilist) = MAX(pdat,ymax(ilist))
            ENDIF
         ENDIF
      ENDDO p_323

!
!3.2.4 Close trajectory files
!----------------------------
!
      
      IF (norec.EQ.tlims(2)) THEN
         n_324: DO n=1,notracks
            CALL close_file(iodat(ilist,n),'A')
         ENDDO n_324
      ENDIF
            
!
!3.2.5 Write parameter files
!---------------------------
!      

      IF (norec.EQ.tlims(2)) THEN

!        ---create parameter file
         WRITE (cfig,'(I12)') noplots+ilist; cfig = ADJUSTL(cfig)
         paramfile = TRIM(plotfile)//'.PDpar'//TRIM(cfig)
         CALL open_file(iunit,paramfile,'OUT','A')

!        ---axis titles
         plottitles(1) = 'distance'
         plottitles(2) = TRIM(outvars(numvar)%long_name)

!        ---main title
         plottitles(3) = TRIM(outfile)

!        ---subtitle
         plottitles(4) = ''

!        ---axis styles and ranges
         istyle = 1; irange = 2
         rangemin(1) = 0.0
         rangemax(1) = MAXVAL(dtrack(ilist,:),&
                     & MASK=dtrack(ilist,:).LT.double_fill)
         drange = ymin(ilist) - ymax(ilist)
         rangemin(2) = ymin(ilist) - 0.05*drange
         rangemax(2) = ymax(ilist) + 0.05*drange

!        ---write parameter file
         CALL write_muplot_par(iunit,notracks,trackfiles(ilist,1:notracks),&
                             & istyle,irange,rangemin,rangemax,&
                             & postatts(ilist)%linestyles(1:notracks))

!        ---close parameter file
         CALL close_file(iunit,'A')
                     
!        --- write to '.vis'-file
         WRITE (io_parvis,9002) 'muplot', TRIM(paramfile), 'y', &
                              & TRIM(plottitles(2))

      ENDIF

   ENDDO ilist_320

ENDDO norec_300

!
!4. Update plot number
!---------------------
!

noplots = noplots + nolists

!
!5. Deallocate arrays
!--------------------
!

DEALLOCATE (trackfiles,iodat,xtrack_old,ytrack_old,dtrack,first)

CALL log_timer_out()


RETURN

9001 FORMAT (2(G15.7,1X))
9002 FORMAT (4(A,1X))

END SUBROUTINE PD_plot

!========================================================================

SUBROUTINE PH_plot
!************************************************************************
!
! *PH_plot* Horizontal particle trajectories
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - post_program
!
! External calls - post_read, write_mumap_coast, write_mumap_par
!
! Module calls - close_file, error_alloc, lim_dims, loop_index, open_file
!
!************************************************************************
!
USE grid  
USE iopars
USE plotpars
USE switches
USE syspars  
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims, loop_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=12) :: cfig, cp
CHARACTER (LEN=lentime) :: ctime
CHARACTER (LEN=leniofile) :: linefile, paramfile
INTEGER :: iend, ilist, istart, iunit, jend, jstart, n, norec, notracks, &
         & ntrackmax, p
REAL :: aspect, gxend, gxstart, gyend, gystart, hlen_unit, x, xborder, xdat, &
      & y, yborder, ydat
REAL (KIND=kndrlong) :: rtime
INTEGER, DIMENSION(3) :: tlims
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: first
CHARACTER (LEN=leniofile), ALLOCATABLE, DIMENSION(:,:) :: symbfiles, trackfiles

INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iodat, iosymb


procname(pglev+1) = 'PH_plot'
CALL log_timer_in()

!
!1.Initialise parameters
!-----------------------
!

CALL time_limits(tlims)

!
!2. Allocate arrays
!------------------
!

ntrackmax = MAXVAL(postatts%notracks)
ALLOCATE (trackfiles(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('trackfiles',2,(/nolists,ntrackmax/),kndchar,lenstr=leniofile)
trackfiles = ''
ALLOCATE (symbfiles(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('symbfiles',2,(/nolists,ntrackmax/),kndchar,lenstr=leniofile)
symbfiles = ''
ALLOCATE (iodat(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('iodat',2,(/nolists,ntrackmax/),kndint)
ALLOCATE (iosymb(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('iosymb',2,(/nolists,ntrackmax/),kndint)
ALLOCATE (first(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('first',2,(/nolists,ntrackmax/),kndlog)
first = .TRUE.

!
!3. Time loop
!------------
!

norec_300: DO norec=tlims(1),tlims(2)

!
!3.1 Read data
!-------------
!

   flag = loop_index(tlims,norec)
   IF (flag.OR.outform.NE.'N') THEN
      CALL post_read(norec,ctime,rtime,.TRUE.)
   ENDIF
   IF (.NOT.flag) CYCLE norec_300
   
!
!3.2 Process data
!----------------
!

   ilist_320: DO ilist=1,nolists

!
!3.2.1 Horizontal locations
!--------------------------
!      

      IF (norec.EQ.tlims(1)) THEN
         IF (postatts(ilist)%horz_unit.EQ.1) THEN
            istart = NINT(postatts(ilist)%gxstart)
            iend = NINT(postatts(ilist)%gxend)
            jstart = NINT(postatts(ilist)%gystart)
            jend = NINT(postatts(ilist)%gyend)
            postatts(ilist)%gxstart = gxcoord(istart,jstart)
            postatts(ilist)%gxend = gxcoord(iend,jend)
            postatts(ilist)%gystart = gycoord(istart,jstart)
            postatts(ilist)%gyend = gycoord(iend,jend)
         ELSE
            CALL coords_to_index(postatts(ilist)%gxstart,&
                               & postatts(ilist)%gystart,istart,jstart)
            CALL coords_to_index(postatts(ilist)%gxend,postatts(ilist)%gyend,&
                               & iend,jend)
         ENDIF
         postatts(ilist)%istart = istart; postatts(ilist)%iend = iend
         postatts(ilist)%jstart = jstart; postatts(ilist)%jend = jend
      ENDIF
      
!
!3.2.1 Create trajectory files
!-----------------------------
!      

      notracks = postatts(ilist)%notracks
      IF (norec.EQ.tlims(1)) THEN
         WRITE (cfig,'(I12)') noplots+ilist; cfig = ADJUSTL(cfig)
         n_321: DO n=1,notracks
            p = postatts(ilist)%ptrack(n)
            WRITE (cp,'(I12)') p; cp = ADJUSTL(cp)
            trackfiles(ilist,n) = TRIM(plotfile)//'.'//TRIM(cp)//'PHcdat'//&
                                & TRIM(cfig)
            CALL open_file(iodat(ilist,n),trackfiles(ilist,n),'OUT','A')
            symbfiles(ilist,n) = TRIM(plotfile)//'.'//TRIM(cp)//'PHbdat'//&
                               & TRIM(cfig)
            CALL open_file(iosymb(ilist,n),symbfiles(ilist,n),'OUT','A')
            WRITE (iosymb(ilist,n),'(A)') '#'
            WRITE (iosymb(ilist,n),'(///)')
            WRITE (iosymb(ilist,n),'(A)') '#'
         ENDDO n_321

      ENDIF
      
!      
!3.2.2 Write trajectory data
!---------------------------
!      

      hlen_unit = postatts(ilist)%hlen_unit
      p_322: DO n=1,notracks
         p = postatts(ilist)%ptrack(n)
         x = outvals2d(p,1); y = outvals2d(p,2)
         flag = ABS(x-fill_value).GT.fill_value_eps
         IF (flag) THEN
            IF (iopt_grid_sph.EQ.0) THEN
               xdat = x*hlen_unit; ydat = y*hlen_unit
            ELSE
               xdat = y; ydat = x
            ENDIF
            IF (first(ilist,n)) THEN
               WRITE (iodat(ilist,n),9001) xdat, ydat, 0, 1.0, &
                    & postatts(ilist)%linestyles(n)
!              ---start symbol               
               WRITE (iosymb(ilist,n),9001) xdat, ydat, 14
               first(ilist,n) = .FALSE.
            ELSE
               WRITE (iodat(ilist,n),9001) xdat, ydat, 1
               IF (norec.EQ.tlims(2)) THEN
                  WRITE (iosymb(ilist,n),9001) xdat, ydat, 15
               ENDIF
            ENDIF
         ENDIF
      ENDDO p_322

!
!3.2.3 Close trajectory files
!----------------------------
!
      
      IF (norec.EQ.tlims(2)) THEN
         n_323: DO n=1,notracks
            CALL close_file(iodat(ilist,n),'A')
            CALL close_file(iosymb(ilist,n),'A')
         ENDDO n_323
      ENDIF

!
!3.2.4 Coast line data
!---------------------
!

      IF (norec.EQ.tlims(2)) THEN
         WRITE (cfig,'(I12)') noplots+ilist; cfig = ADJUSTL(cfig)
!        ---create coast line file
         IF (coast_def) THEN
            linefile = TRIM(plotfile)//'.PHcdat'//TRIM(cfig)
            CALL open_file(iunit,linefile,'OUT','A')
         ELSE
            linefile = TRIM(coast_name)
         ENDIF
!        ---write data
         IF (coast_def) THEN
            CALL write_mumap_coast(iunit,ilist,linefile)
         ENDIF

      ENDIF
      
!
!3.2.5 Write parameter files
!---------------------------
!      

      IF (norec.EQ.tlims(2)) THEN
         
!        ---plot borders
         gxstart = postatts(ilist)%gxstart
         gxend = postatts(ilist)%gxend
         gystart = postatts(ilist)%gystart
         gyend = postatts(ilist)%gyend
         xborder = border*(gxend-gxstart)
         yborder = border*(gyend-gystart)
         plotcorners(1) = gxstart - xborder
         plotcorners(2) = gystart - yborder
         plotcorners(3) = gxend + xborder
         plotcorners(4) = gyend + yborder
         plotcorners = hlen_unit*plotcorners
         aspect = 1.0

!        ---plot titles
         plottitles(1) = TRIM(outfile)
         plottitles(2) = 'Particle trajectory'
         plottitles(3) = ''

!        ---parameter file         
         n_325: DO n=1,notracks
!           --create
            p = postatts(ilist)%ptrack(n)
            WRITE (cp,'(I12)') p; cp = ADJUSTL(cp)
            paramfile = TRIM(plotfile)//'.'//TRIM(cp)//'PHpar'//TRIM(cfig)
            CALL open_file(iunit,paramfile,'OUT','A')
!           ---plot title
            plottitles(2) = 'Trajectory particle '//TRIM(cp)
!           ---write         
            CALL write_mumap_par(iunit,.TRUE.,iopt_grid_sph,aspect,3,&
                               & (/'c','b','c'/),(/trackfiles(ilist,n),&
                               & symbfiles(ilist,n),linefile/))
!           --close parameter files
            CALL close_file(iunit,'A')
!           -- write to '.vis'-file
            WRITE (io_parvis,9002) 'mumap2', TRIM(paramfile), 'y', &
                                 & TRIM(plottitles(2))
         ENDDO n_325

      ENDIF

   ENDDO ilist_320

ENDDO norec_300

!
!4. Update plot number
!---------------------
!

noplots = noplots + nolists

!
!5. Deallocate arrays
!--------------------
!

DEALLOCATE (trackfiles,symbfiles,iodat,iosymb,first)

CALL log_timer_out()


RETURN

9001 FORMAT (2(G15.7,1X),I2,1X,G15.7,1X,I2)
9002 FORMAT (4(A,1X))

END SUBROUTINE PH_plot

!========================================================================

SUBROUTINE post_end
!************************************************************************
!
! *post_end* Finalise postprocessing of a data file
!
! Author - Patrick Luyten
!
! Description - read global, coordinate and variable attributes
!             - allocate output grid arrays
!             - read coordinate arrays 
!
! Reference -
!
! Calling program - post_program, post_start
!
! Module calls - close_filepars
!
!************************************************************************
!
USE depths
USE grid
USE iopars
USE plotpars
USE switches
USE inout_routines, ONLY: close_filepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'post_end'
CALL log_timer_in()

!
!1. Close data files
!-------------------
!

CALL close_filepars(filepars)

!
!2. Deallocate
!-------------
!

DEALLOCATE (coordvars,outvars)
IF (nodim.GT.0) THEN
   DEALLOCATE (depout,gxcoord,gycoord,gzcoord,gzcoord_ini,maskatc_int)
   IF (outgpars%vcoord.EQ.2) THEN
      IF (iopt_grid_vtype.EQ.2) THEN
         DEALLOCATE (gsigcoordatc)
      ELSEIF (iopt_grid_vtype.EQ.3) THEN
         DEALLOCATE (gsout)
      ENDIF
   ENDIF
   IF (outgpars%vcoord.EQ.2.AND.time_grid) DEALLOCATE(zeta)
   IF (packing) DEALLOCATE (lmaskout)
   IF (.NOT.gridded) DEALLOCATE (station_names)
ENDIF
IF (.NOT.traject) THEN
   SELECT CASE (nodim)
      CASE (0); DEALLOCATE (outvals1d)
      CASE (2)
         IF (gridded) THEN
            DEALLOCATE (outvals3d)
         ELSE
            DEALLOCATE (outvals2d)
         ENDIF
      CASE (3)
         IF (gridded) THEN
            DEALLOCATE (outvals4d)
         ELSE
            DEALLOCATE (outvals3d)
         ENDIF
   END SELECT
ELSE
   DEALLOCATE (outvals2d)
ENDIF
   
DEALLOCATE (vecids)

CALL log_timer_out()


RETURN

END SUBROUTINE post_end

!========================================================================

SUBROUTINE post_read(timerec,ctime,rtime,info)
!************************************************************************
!
! *post_read* Read output data for postprocessing
!
! Author - Patrick Luyten
!
! Description -
!
! Reference -
!
! Calling program - HP_plot, HT_plot, TH_plot, TS_plot, TZ_plot, VT_plot,
!                   ZP_plot
!
! Module calls - add_secs_to_date, date_to_year, diff_dates, read_submod,
!                read_time, read_vars, year_to_date
!
!************************************************************************
!
USE depths
USE grid
USE iopars
USE physpars
USE plotpars
USE switches
USE syspars
USE timepars
USE inout_routines, ONLY: read_submod, read_time, read_vars
USE time_routines, ONLY: add_secs_to_date, date_to_year, diff_dates, &
                       & log_timer_in, log_timer_out, year_to_date

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: timerec
LOGICAL, INTENT(IN) :: info
CHARACTER (LEN=lentime), INTENT(OUT) :: ctime
REAL (KIND=kndrlong), INTENT(OUT) :: rtime

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*timerec*   INTEGER Time record number
!*ctime*     CHAR    Date/time in character format
!*rtime*     REAL    Time in real format
!*info*      LOGICAL Enables/disables info for log file
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=lentime) :: ctimeout
INTEGER :: ivar, k, lloglev1, lloglev2, nosecs
REAL :: dtime, tconv
REAL (KIND=kndrlong) :: rtimeout


IF (.NOT.info) THEN
   lloglev1 = loglev1; lloglev2 = loglev2
   loglev1 = 0; loglev2 = 0
ELSE
   procname(pglev+1) = 'post_read'
   CALL log_timer_in()
ENDIF

!
!1. New data time
!----------------
!
!1.1 Read time
!-------------
!

IF (gtime_format.EQ.0) THEN
   CALL read_time(ctimeout,filepars,timerec)
   IF (info.AND.loglev1.GE.pglev) THEN
      WRITE (iolog,'(A)') 'time: '//ctimeout
   ENDIF
ELSE
   CALL read_time(rtimeout,filepars,timerec)
   IF (info.AND.loglev1.GE.pglev) THEN
      WRITE (iolog,'(A,G12.3)') 'time:', rtimeout
   ENDIF
ENDIF

!
!1.2 Convert to appropriate unit
!-------------------------------
!

ctime = ''; rtime = 0.0

!---input date in string format
IF (gtime_format.EQ.0) THEN
   IF (ptime_format.EQ.0) THEN
      ctime = ctimeout
   ELSEIF (ptime_format.EQ.8) THEN
      CALL date_to_year(ctimeout,rtime)
   ELSE
      CALL diff_dates(PRefDateTime,ctimeout,ptime_format,rtime=rtime)
   ENDIF

!---input date in real format
ELSEIF (gtime_format.LT.8) THEN
   IF (ptime_format.EQ.0.OR.ptime_format.EQ.8) THEN
      nosecs = time_convert(gtime_format)*INT(rtimeout)
      CALL add_secs_to_date(outgpars%refdate,ctimeout,nosecs,1.0)
      dtime = time_convert(gtime_format)*(rtimeout-INT(rtimeout))
      CALL add_secs_to_date(ctimeout,ctime,1,dtime)
      IF (ptime_format.EQ.8) CALL date_to_year(ctime,rtime)
   ELSE
      tconv = time_convert(gtime_format)/REAL(time_convert(ptime_format))
      rtime = tconv*rtimeout - offset
   ENDIF

!---input date as year date
ELSE
   IF (ptime_format.EQ.0) THEN
      CALL year_to_date(rtimeout,ctime)
   ELSEIF (ptime_format.LT.8) THEN
      CALL year_to_date(rtimeout,ctimeout)
      CALL diff_dates(PRefDateTime,ctimeout,ptime_format,rtime=rtime)
   ELSE
      rtime = rtimeout
   ENDIF
ENDIF

!
!2. Time grid
!------------
!

IF (time_grid) THEN

   filepars%timerec = timerec
   
!
!2.1 Sigma coordinate grid
!-------------------------
!   

   IF (outgpars%vcoord.EQ.2) THEN
!     ---elevation data
      CALL read_submod(zeta,filepars,filepars%zvarid,varatts=coordvars,&
                  & maskvals=maskatc_int)
!     ---update Z-coordinates
      IF (iopt_grid_vtype.EQ.2) THEN
         k_211: DO k=1,nzout
            WHERE (maskatc_int)
               gzcoord(:,:,k) = gsigcoordatc(k)*(depout+zeta) + zeta
            END WHERE
         ENDDO k_211
      ELSEIF (iopt_grid_vtype.EQ.3) THEN
         k_212: DO k=1,nzout
            WHERE (maskatc_int)
               gzcoord(:,:,k) = gsout(:,:,k)*(depout+zeta) + zeta
            END WHERE
         ENDDO k_212
      ENDIF

!
!2.2 Z-coordinate grid
!---------------------
!

   ELSEIF (outgpars%vcoord.EQ.1) THEN
      CALL read_submod(gzcoord,filepars,filepars%zvarid,varatts=coordvars,&
           & maskvals=maskatc_int)
      IF (packing) THEN
         k_220: DO k=1,nzout
            WHERE (.NOT.maskatc_int)
               gzcoord(:,:,k) = 0.0
            END WHERE
         ENDDO k_220
      ENDIF
   ENDIF

ELSEIF (nodim.EQ.3) THEN

   gzcoord = gzcoord_ini

ENDIF

!
!3. Read data
!------------
!

filepars%timerec = timerec

IF (.NOT.traject) THEN
   SELECT CASE (nodim)
      CASE (0); CALL read_vars(outvals1d,filepars,0,vecids=vecids)
      CASE (2)
         IF (gridded) THEN
            CALL read_submod(outvals3d,filepars,0,vecids=vecids,&
                           & maskvals=maskatc_int)
         ELSE
            CALL read_vars(outvals2d,filepars,0,vecids=vecids)
         ENDIF
      CASE (3)
         IF (gridded) THEN
            CALL read_submod(outvals4d,filepars,0,vecids=vecids,&
                           & maskvals=maskatc_int)
         ELSE
            CALL read_vars(outvals3d,filepars,0,vecids=vecids)
         ENDIF
   END SELECT
ELSE
   CALL read_vars(outvals2d,filepars,0,vecids=vecids)
   ivar_310: DO ivar=1,novars
      IF (TRIM(outvars(ivar)%f90_name).EQ.'part_age') THEN
         outvals2d(:,ivar) = tconv*outvals2d(:,ivar)
      ENDIF
   END DO ivar_310
ENDIF

!
!4. Mask dry points
!------------------
!

IF (.NOT.(traject.OR.packing).AND.gridded) THEN
   ivar_410: DO ivar=1,novars 
      IF (nodim.EQ.2) THEN
         WHERE (.NOT.maskatc_int)
            outvals3d(:,:,ivar) = fill_value
         END WHERE
      ELSEIF (nodim.EQ.3) THEN
         k_411: DO k=1,nzout
            WHERE (.NOT.maskatc_int)
               outvals4d(:,:,k,ivar) = fill_value
            END WHERE
         ENDDO k_411
      ENDIF
   ENDDO ivar_410
ENDIF

IF (.NOT.info) THEN
   loglev1 = lloglev1; loglev2 = lloglev2
ELSE
   CALL log_timer_out()
ENDIF


RETURN

END SUBROUTINE post_read

!========================================================================

SUBROUTINE post_start
!************************************************************************
!
! *post_start* Initialise postprocessing of a data file
!
! Author - Patrick Luyten
!
! Description - read global, coordinate and variable attributes
!             - allocate output grid arrays
!             - read coordinate arrays
!             - allocate array for storing output data
!
! Reference -
!
! Calling program - post_program
!
! External calls - post_end
!
! Module calls - add_secs_to_date, cf90_inq_varid, cf90_get_var, convert_date,
!                error_alloc, error_alloc_struc, inquire_modfil_atts,
!                inquire_modfil_dims, inquire_modfil_vars, open_filepars,
!                read_metadata_line, read_submod, read_time, read_vars,
!                varatts_init
!
!************************************************************************
!
USE depths
USE grid
USE iopars
USE plotpars
USE physpars
USE switches
USE syspars
USE timepars
USE cf90_routines, ONLY: cf90_inq_varid, cf90_get_var
USE cif_routines, ONLY: conv_from_chars
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE inout_routines, ONLY: inquire_modfil_atts, inquire_modfil_dims, &
                        & inquire_modfil_vars, open_filepars, &
                        & read_metadata_line, read_submod, read_time, read_vars
USE time_routines, ONLY: add_secs_to_date, convert_date, log_timer_in, &
                       & log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lenname) :: cname
CHARACTER (LEN=lendesc) :: ftype, standard_name
CHARACTER (LEN=lenunit) :: tunit, tunitx
CHARACTER (LEN=lentime) :: refdate
INTEGER :: i, icoord, idim, isep, istat, iunit, ivar, j, k, l, linenum, n, &
         & nodims, nowetout, nrank, numvars, varid
REAL (KIND=kndrlong) :: rtime1, rtime2
CHARACTER (LEN=lenname), SAVE, ALLOCATABLE, DIMENSION(:) :: dimnames
CHARACTER (LEN=lendesc), DIMENSION(15) :: cvals
INTEGER, DIMENSION(3) :: ndims
INTEGER, DIMENSION(7) :: irefdate
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: outdims


procname(pglev+1) = 'post_start'
CALL log_timer_in()

!
!1. Open data file
!-----------------
!
!---close data file if needed
IF (outform.NE.'N') THEN
   INQUIRE (FILE=filepars%filename,OPENED=flag)
   IF (flag) CALL post_end
ENDIF

!---open data file
CALL open_filepars(filepars)

!
!2. Number of variables, dimensions
!----------------------------------
!

iunit = filepars%iunit
outform = filepars%form
CALL inquire_modfil_dims(iunit,outform,ndimensions=nodims,nvariables=novars,&
                       & ncoordinates=nocoords,floattype=filepars%floattype)
filepars%nocoords = nocoords
filepars%novars = novars

!
!3. Allocate
!-----------
!

ALLOCATE (outdims(nodims),STAT=errstat)
CALL error_alloc('outdims',1,(/nodims/),kndint)
outdims = 0

ALLOCATE (dimnames(nodims),STAT=errstat)
CALL error_alloc('dimnames',1,(/nodims/),kndchar)
dimnames = ''

ALLOCATE (coordvars(nocoords),STAT=errstat)
CALL error_alloc_struc('coordvars',1,(/nocoords/),'VariableAtts')
CALL varatts_init(coordvars)

ALLOCATE (outvars(novars),STAT=errstat)
CALL error_alloc_struc('outvars',1,(/novars/),'VariableAtts')
CALL varatts_init(outvars)

!
!4. Array dimensions
!-------------------
!
!---dimension values and names
CALL inquire_modfil_dims(iunit,outform,dimensions=outdims,dimnames=dimnames,&
                       & vcoord=vcoord,ftype=ftype)

!---output dimension and type of of output grid 
ncout = 0; nrout = 0; nowetout = 0; nostats = 0; nzout = 0
gridded = .TRUE.; packing = .FALSE.; traject = .FALSE. 
idim_410: DO idim=1,nodims
   SELECT CASE (TRIM(dimnames(idim)))
      CASE ('time')
         nstepout = outdims(idim)
         nodim = 0
      CASE ('x','lon')
         ncout = outdims(idim)
         iopt_grid_htype = 2
         nodim = 2
      CASE ('y','lat')
         nrout = outdims(idim)
      CASE ('nc')
         ncout = outdims(idim)
         iopt_grid_htype = 3
         nodim = 2
      CASE ('nr')
         nrout = outdims(idim)
      CASE ('seapoint')
         nowetout = outdims(idim)
         packing = .TRUE.
      CASE ('station')
         nostats = outdims(idim)
      CASE ('lev','gsig','z')
         nzout = outdims(idim)
         nodim = 3
      CASE ('trajectory')
         nopartout = outdims(idim)
         nodim = 2
   END SELECT
ENDDO idim_410

IF (TRIM(ftype).NE.'') THEN
   SELECT CASE (TRIM(ftype))
      CASE ('timeSeries'); gridded = .TRUE.
      CASE ('timeSeriesProfile'); gridded = .TRUE.
      CASE ('trajectory'); traject = .TRUE.
   END SELECT
ENDIF

!---store
outgpars%ncout = ncout; outgpars%nrout = nrout; outgpars%nzout = nzout
outgpars%nowetout = nowetout; outgpars%nostats = nostats
outgpars%nstepout = nstepout; outgpars%nodim = nodim
outgpars%gridded = gridded; outgpars%packing = packing
outgpars%vcoord = vcoord
filepars%packing = packing

!
!5. Coordinate and variable attributes
!-------------------------------------
!
!---coordinates
icoord_511: DO icoord=1,nocoords
   CALL inquire_modfil_atts(iunit,outform,icoord,&
                          & nrank=coordvars(icoord)%nrank,&
                          & data_type=coordvars(icoord)%data_type,&
                          & global_dims=coordvars(icoord)%global_dims,&
                          & f90_name=coordvars(icoord)%f90_name,&
                          & standard_name=coordvars(icoord)%standard_name,&
                          & long_name=coordvars(icoord)%long_name,&
                          & units=coordvars(icoord)%units,&
                          & fill=coordvars(icoord)%fill,&
                          & fill_value=coordvars(icoord)%fill_value)
ENDDO icoord_511

icoord_512: DO icoord=2,nocoords
   nrank = coordvars(icoord)%nrank
   SELECT CASE (coordvars(icoord)%f90_name)
      CASE('eta')
         IF (nodim.EQ.3) THEN
            coordvars(icoord)%nrank = coordvars(icoord)%nrank - 1
         ENDIF
      CASE('x')
         iopt_grid_sph = 0
      CASE ('lon')
         iopt_grid_sph = 1
   END SELECT
ENDDO icoord_512

!---variables
ivar_521: DO ivar=1,novars
   CALL inquire_modfil_atts(iunit,outform,nocoords+ivar,&
                          & nrank=outvars(ivar)%nrank,&
                          & data_type=outvars(ivar)%data_type,&
                          & global_dims=outvars(ivar)%global_dims,&
                          & f90_name=outvars(ivar)%f90_name,&
                          & standard_name=outvars(ivar)%standard_name,&
                          & long_name=outvars(ivar)%long_name,&
                          & units=outvars(ivar)%units,&
                          & vector_name=outvars(ivar)%vector_name,&
                          & fill=outvars(ivar)%fill,&
                          & fill_value=outvars(ivar)%fill_value)
   outvars(ivar)%nrank = outvars(ivar)%nrank - 1
   IF (outvars(ivar)%fill) fill_value = outvars(ivar)%fill_value
ENDDO ivar_521

!---fill value
fill_value_eps = 0.00001*ABS(fill_value)
filepars%fill = ANY(outvars%fill)

!
!6. Time parameters
!------------------
!
!---check for time grid
time_grid = .FALSE.
IF (nodim.EQ.3) THEN
   icoord_610: DO icoord=1,nocoords
      IF (TRIM(coordvars(icoord)%f90_name).EQ.'eta') time_grid = .TRUE.
   ENDDO icoord_610
ENDIF

!---time unit
icoord_620: DO icoord=1,nocoords
   SELECT CASE (TRIM(coordvars(icoord)%f90_name))
      CASE ('time')
         tunit = coordvars(icoord)%units
         filepars%timeid = icoord
      CASE ('eta','zout')
         filepars%zvarid = icoord
   END SELECT
ENDDO icoord_620
tunitx = tunit
isep = INDEX(tunitx,' ')
tunitx = tunit(1:isep-1)
SELECT CASE (TRIM(tunitx))
   CASE ('seconds'); gtime_format = 1
   CASE ('minutes'); gtime_format = 2
   CASE ('hours'); gtime_format = 3
   CASE ('days'); gtime_format = 4
   CASE ('weeks'); gtime_format = 5
   CASE ('months'); gtime_format = 6
   CASE ('years'); gtime_format = 7
END SELECT

!---time step
filepars%timerec = 0; filepars%maxrecs = nstepout
IF (outform.NE.'N') THEN
   REWIND (iunit)
   istat = 1
   DO WHILE (istat.EQ.1)
      CALL read_metadata_line(iunit,outform,cvals,numvars,cname,linenum,&
                            & istat)
      IF (TRIM(cname).EQ.'time_step') THEN
         CALL conv_from_chars(cvals(1),deltout,1)
      ENDIF
   ENDDO
ELSE
   CALL read_time(rtime1,filepars)
   CALL read_time(rtime2,filepars)
   deltout = rtime2 - rtime1
   filepars%timerec = 0
ENDIF

!---reference date
tunitx = tunit
l = LEN_TRIM(tunitx)
refdate(1:19) = tunitx(l-18:l)
refdate(20:23) = ' 000'
irefdate = convert_date(refdate)
outgpars%refdate = convert_date(irefdate)

!---start date
IF (outform.NE.'N') THEN
   REWIND (iunit)
   istat = 1
   DO WHILE (istat.EQ.1)
      CALL read_metadata_line(iunit,outform,cvals,numvars,cname,linenum,&
                            & istat)
      IF (TRIM(cname).EQ.'StartDate') outgpars%startdate = TRIM(cvals(1))
   ENDDO
ELSE
   CALL add_secs_to_date(outgpars%refdate,outgpars%startdate,&
                       & time_convert(gtime_format),REAL(rtime1))
ENDIF

!
!7. Vertical grid
!----------------
!
!---type of vertical grid
iopt_grid_vtype = 0
icoord_710: DO icoord=1,nocoords
   SELECT CASE (TRIM(coordvars(icoord)%f90_name))
      CASE ('lev','gsig')
         standard_name = coordvars(icoord)%standard_name
         IF (TRIM(standard_name).EQ.'ocean_sigma_coordinate') THEN
            iopt_grid_vtype = 2
         ELSEIF (TRIM(standard_name).EQ.'ocean_s_coordinate') THEN
            iopt_grid_vtype = 3
         ENDIF
   END SELECT
ENDDO icoord_710

!---parameters for SH formulation
IF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
   IF (outform.NE.'N') THEN
      REWIND (iunit)
      istat = 1
      DO WHILE (istat.EQ.1)
         CALL read_metadata_line(iunit,outform,cvals,numvars,cname,linenum,&
                               & istat)
         IF (TRIM(cname).EQ.'SH_parameters') THEN
            CALL conv_from_chars(cvals(1),theta_SH,1)
            CALL conv_from_chars(cvals(2),b_SH,2)
            CALL conv_from_chars(cvals(3),hcrit_SH,3)
         ENDIF
      ENDDO
   ELSE
      CALL cf90_inq_varid(iunit,'theta_SH',varid)
      CALL cf90_get_var(iunit,varid,theta_SH,0)
      CALL cf90_inq_varid(iunit,'b_SH',varid)
      CALL cf90_get_var(iunit,varid,b_SH,0)
      CALL cf90_inq_varid(iunit,'hcrit_SH',varid)
      CALL cf90_get_var(iunit,varid,hcrit_SH,0)
   ENDIF
ENDIF

!
!8. Station names
!----------------
!

IF (.NOT.gridded.AND.nodim.GT.0) THEN
   ALLOCATE (station_names(nostats),STAT=errstat)
   CALL error_alloc('station_names',1,(/nostats/),kndchar,lenstr=lenname)
   CALL inquire_modfil_vars(iunit,outform,station_names=station_names)
ENDIF

!
!9. Allocate coordinate arrays
!-----------------------------
!

IF (nodim.GT.0) THEN
   IF (gridded) THEN
      ndims = (/ncout,nrout,nzout/)
   ELSE
      ndims = (/nostats,1,nzout/)
   ENDIF
   ALLOCATE (gxcoord(ndims(1),ndims(2)),STAT=errstat)
   CALL error_alloc('gxcoord',2,ndims(1:2),kndrtype)
   gxcoord = 0.0
   ALLOCATE (gycoord(ndims(1),ndims(2)),STAT=errstat)
   CALL error_alloc('gycoord',2,ndims(1:2),kndrtype)
   gycoord = 0.0
   IF (vcoord.EQ.2) THEN
      IF (iopt_grid_vtype.EQ.2) THEN
         ALLOCATE (gsigcoordatc(ndims(3)),STAT=errstat)
         CALL error_alloc('gsigcoordatc',1,ndims(3),kndrtype)
         IF (nzout.GT.0) gsigcoordatc = 0.0
      ELSEIF (iopt_grid_vtype.EQ.3) THEN
         ALLOCATE (gsout(ndims(1),ndims(2),ndims(3)),STAT=errstat)
         CALL error_alloc('gsigcoordatc',3,ndims,kndrtype)
         IF (nzout.GT.0) gsout = 0.0
      ENDIF
   ENDIF
   ALLOCATE (depout(ndims(1),ndims(2)),STAT=errstat)
   CALL error_alloc('depout',2,ndims(1:2),kndrtype)
   depout = 0.0
   IF (vcoord.EQ.2.AND.time_grid) THEN
      ALLOCATE (zeta(ndims(1),ndims(2)),STAT=errstat)
      CALL error_alloc('zeta',2,ndims(1:2),kndrtype)
      zeta = 0.0
   ENDIF
   ALLOCATE (gzcoord_ini(ndims(1),ndims(2),ndims(3)),STAT=errstat)
   CALL error_alloc('gzcoord_ini',3,ndims(1:3),kndrtype)
   IF (nzout.GT.0) gzcoord_ini = 0.0
   ALLOCATE (gzcoord(ndims(1),ndims(2),ndims(3)),STAT=errstat)
   CALL error_alloc('gzcoord',3,ndims(1:3),kndrtype)
   IF (nzout.GT.0) gzcoord = 0.0
   ALLOCATE (maskatc_int(ndims(1),ndims(2)),STAT=errstat)
   CALL error_alloc('maskatc_int',2,ndims(1:2),kndlog)
   maskatc_int = .FALSE.
   IF (packing) THEN
      ALLOCATE (lmaskout(nowetout),STAT=errstat)
      CALL error_alloc('lmaskout',1,(/nowetout/),kndint)
   ENDIF
ENDIF

!
!10. Read coordinate arrays
!--------------------------
!

varid = 1
IF (nodim.GT.0) THEN
   IF (gridded) THEN
!     ---horizontal coordinates (regular grid)
      IF (iopt_grid_htype.LE.2) THEN
         varid = varid + 1
         CALL read_vars(gxcoord(:,1),filepars,varid,coordvars)
          j_1010: DO j=2,nrout
            gxcoord(:,j) = gxcoord(:,1)
         ENDDO j_1010
         varid = varid + 1
         CALL read_vars(gycoord(1,:),filepars,varid,coordvars)
          i_1020: DO i=2,ncout
            gycoord(i,:) = gycoord(1,:)
         ENDDO i_1020
!     ---horizontal coordinates (curvilinear grid)
      ELSE
         varid = varid + 1
         CALL read_vars(gxcoord,filepars,varid,coordvars)
         varid = varid + 1
         CALL read_vars(gycoord,filepars,varid,coordvars)
      ENDIF
!     ---land mask array
      IF (packing) THEN
         varid = varid + 1
         CALL read_vars(lmaskout,filepars,varid,coordvars)
         l_1030: DO l=1,nowetout
            n = lmaskout(l)
            j = (n-1)/ncout + 1
            i = n - (j-1)*ncout
            maskatc_int(i,j) = .TRUE.
         ENDDO l_1030
         varid = varid + 1
         CALL read_submod(depout,filepars,varid,varatts=coordvars,&
                        & maskvals=maskatc_int)
      ENDIF
!     ---bathymetry
      IF (.NOT.packing) THEN
         varid = varid + 1
         CALL read_vars(depout,filepars,varid,varatts=coordvars)
          WHERE (ABS(depout-fill_value).GT.fill_value_eps)
            maskatc_int = .TRUE.
         ELSEWHERE
            depout = 0.0
         END WHERE
      ENDIF
      
   ELSE
!     ---horizontal coordinates and bathymetry (irregular grid)
      varid = varid + 1
      CALL read_vars(gxcoord(:,1),filepars,varid,coordvars)
      varid = varid + 1
      CALL read_vars(gycoord(:,1),filepars,varid,coordvars)
      varid = varid + 1
      CALL read_vars(depout(:,1),filepars,varid,coordvars)
   ENDIF

!  ---vertical coordinate
   IF (nzout.GT.0) THEN
      varid = varid + 1
!     --sigma-coordinates
      IF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.2) THEN
         CALL read_vars(gsigcoordatc,filepars,varid,coordvars)
      ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
         CALL read_submod(gsout,filepars,varid,varatts=coordvars,&
                        & maskvals=maskatc_int)
         IF (packing) THEN
            k_1040: DO k=1,nzout
               WHERE (.NOT.maskatc_int)
                  gsout(:,:,k) = 0.0
               END WHERE
            ENDDO k_1040
         ENDIF
      ELSEIF (vcoord.EQ.1) THEN
         CALL read_submod(gzcoord_ini,filepars,varid,varatts=coordvars,&
              & maskvals=maskatc_int)
         IF (packing) THEN
            k_1050: DO k=1,nzout
               WHERE (.NOT.maskatc_int)
                  gzcoord_ini(:,:,k) = 0.0
               END WHERE
            ENDDO k_1050
         ENDIF
      ENDIF
   ENDIF
   
ENDIF

!
!11. Z-coordinates
!-----------------
!

IF (nzout.GT.0.AND.vcoord.EQ.2) THEN
   IF (iopt_grid_vtype.EQ.2) THEN
      k_1110: DO k=1,nzout
         WHERE (maskatc_int)
            gzcoord_ini(:,:,k) = gsigcoordatc(k)*depout
         END WHERE
      ENDDO k_1110
   ELSEIF (iopt_grid_vtype.EQ.3) THEN
      k_1120: DO k=1,nzout
         gzcoord_ini(:,:,k) = gsout(:,:,k)*depout
      ENDDO k_1120
   ENDIF
ENDIF

!
!12. Allocate data output array
!------------------------------
!

IF (.NOT.traject) THEN
   
   SELECT CASE (nodim)
      CASE (0)
         ALLOCATE (outvals1d(novars),STAT=errstat)
         CALL error_alloc('outvals1d',1,(/novars/),kndrtype)
      CASE (2)
         IF (gridded) THEN
            ALLOCATE (outvals3d(ncout,nrout,novars),STAT=errstat)
            CALL error_alloc('outvals3d',3,(/ncout,nrout,novars/),kndrtype)
         ELSE
            ALLOCATE (outvals2d(nostats,novars),STAT=errstat)
            CALL error_alloc('outvals2d',2,(/nostats,novars/),kndrtype)
         ENDIF
      CASE (3)
         IF (gridded) THEN
            ALLOCATE (outvals4d(ncout,nrout,nzout,novars),STAT=errstat)
            CALL error_alloc('outvals4d',4,(/ncout,nrout,nzout,novars/),&
                           & kndrtype)
         ELSE
            ALLOCATE (outvals3d(nostats,nzout,novars),STAT=errstat)
            CALL error_alloc('outvals3d',3,(/nostats,nzout,novars/),kndrtype)
         ENDIF
   END SELECT

ELSE

   ALLOCATE (outvals2d(nopartout,novars),STAT=errstat)
   CALL error_alloc('outvals2d',2,(/nopartout,novars/),kndrtype)

ENDIF

!
!13. Variables ids
!-----------------
!

ALLOCATE (vecids(novars),STAT=errstat)
CALL error_alloc('vecids',1,(/novars/),kndint)
vecids = (/(nocoords+ivar,ivar=1,novars)/)

!
!14. Deallocate
!--------------
!

DEALLOCATE (dimnames,outdims)

CALL log_timer_out()


RETURN

END SUBROUTINE post_start

!========================================================================

SUBROUTINE PV_plot
!************************************************************************
!
! *PV_plot* Particle trajectories along vertical transects
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - post_program
!
! External calls - coords_to_index, post_read, time_limits, write_muplot_par
!
! Module calls - close_file, convert_loc_to_char, distance_pts, error_alloc,
!                loop_index, open_file
!
!************************************************************************
!
USE grid
USE iopars
USE plotpars
USE switches
USE syspars
USE error_routines, ONLY : error_alloc
USE grid_routines, ONLY: convert_loc_to_char, distance_pts
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: loop_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=12) :: cfig, cp
CHARACTER (LEN=lendesc) :: cloc
CHARACTER (LEN=leniofile) :: linefile, paramfile
CHARACTER (LEN=lentime) :: ctime
CHARACTER (LEN=12), DIMENSION(4) :: cpos
CHARACTER (LEN=21), DIMENSION(2) :: geoloc
INTEGER :: horz_unit, i, iend, ilist, istart, iunit, j, jend, jstart, &
         & maxpoints, n, norec, notracks, ntrackmax, numpoints, p
REAL :: aspect, deli, delj, dels, dist, gxend, gxstart, gyend, gystart, &
      & hlen_unit, phi, s, x, xborder, y, yborder, z, z1, z2
REAL (KIND=kndrlong) :: rtime
CHARACTER (LEN=leniofile), ALLOCATABLE, DIMENSION(:,:) :: symbfiles, trackfiles
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: first
INTEGER, DIMENSION(nolists) :: itrtype, nopoints
INTEGER, DIMENSION(3) :: tlims
REAL, DIMENSION(nolists) :: dmax, dmin, trangle, zmax, zmin
INTEGER, ALLOCATABLE, DIMENSION(:) :: itrans, jtrans
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iodat, iosymb
REAL, ALLOCATABLE, DIMENSION(:) :: deptrans, xtrans, ytrans
REAL, ALLOCATABLE, DIMENSION(:,:) :: dtrans 


procname(pglev+1) = 'PV_plot'
CALL log_timer_in()

!
!1.Time limiters
!---------------
!

CALL time_limits(tlims)

!
!2. Initialise
!-------------
!

ilist_200: DO ilist=1,nolists

!
!2.1  Horizontal grid locations
!------------------------------
!
   
   IF (postatts(ilist)%horz_unit.EQ.1) THEN
      istart = NINT(postatts(ilist)%gxstart)
      iend = NINT(postatts(ilist)%gxend)
      jstart = NINT(postatts(ilist)%gystart)
      jend = NINT(postatts(ilist)%gyend)
      postatts(ilist)%gxstart = gxcoord(istart,jstart)
      postatts(ilist)%gxend = gxcoord(iend,jend)
      postatts(ilist)%gystart = gycoord(istart,jstart)
      postatts(ilist)%gyend = gycoord(iend,jend)
      postatts(ilist)%istart = istart
      postatts(ilist)%iend = iend
      postatts(ilist)%jstart = jstart
      postatts(ilist)%jend = jend
   ELSE
      CALL coords_to_index(postatts(ilist)%gxstart,postatts(ilist)%gystart,&
                         & istart,jstart)
      CALL coords_to_index(postatts(ilist)%gxend,postatts(ilist)%gyend,&
                         & iend,jend)
      postatts(ilist)%istart = istart
      postatts(ilist)%iend = iend
      postatts(ilist)%jstart = jstart
      postatts(ilist)%jend = jend
   ENDIF

!
!2.2 Number of locations along transect
!--------------------------------------
!

   nopoints(ilist) = MAX(ABS(iend-istart),ABS(jend-jstart)) + 1

ENDDO ilist_200

!
!3. Allocate arrays
!------------------
!

maxpoints = MAXVAL(nopoints)
ALLOCATE (dtrans(nolists,maxpoints),STAT=errstat)
CALL error_alloc('dtrans',2,(/nolists,maxpoints/),kndrtype)
ntrackmax = MAXVAL(postatts%notracks)
ALLOCATE (trackfiles(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('trackfiles',2,(/nolists,ntrackmax/),kndchar,lenstr=leniofile)
trackfiles = ''
ALLOCATE (symbfiles(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('symbfiles',2,(/nolists,ntrackmax/),kndchar,lenstr=leniofile)
symbfiles = ''
ALLOCATE (iodat(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('iodat',2,(/nolists,ntrackmax/),kndint)
ALLOCATE (iosymb(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('iosymb',2,(/nolists,ntrackmax/),kndint)
ALLOCATE (first(nolists,ntrackmax),STAT=errstat)
CALL error_alloc('first',2,(/nolists,ntrackmax/),kndlog)
first = .TRUE.

!
!4. Define transect data
!-----------------------
!

ilist_400: DO ilist=1,nolists

!
!4.1 Type of transect
!--------------------
!   

   istart = postatts(ilist)%istart
   iend = postatts(ilist)%iend
   jstart = postatts(ilist)%jstart
   jend = postatts(ilist)%jend
   IF (jstart.EQ.jend) THEN
      itrtype(ilist) = 1
   ELSEIF (istart.EQ.iend) THEN
      itrtype(ilist) = 2
   ELSE
      itrtype(ilist) = 3
   ENDIF

!
!4.2 Allocate arrays
!-------------------
!

   numpoints = nopoints(ilist)
   ALLOCATE (itrans(numpoints),STAT=errstat)
   CALL error_alloc('itrans',1,(/numpoints/),kndint)
   ALLOCATE (jtrans(numpoints),STAT=errstat)
   CALL error_alloc('jtrans',1,(/numpoints/),kndint)
   ALLOCATE (xtrans(numpoints),STAT=errstat)
   CALL error_alloc('xtrans',1,(/numpoints/),kndrtype)
   ALLOCATE (ytrans(numpoints),STAT=errstat)
   CALL error_alloc('ytrans',1,(/numpoints/),kndrtype)
   ALLOCATE (deptrans(numpoints),STAT=errstat)
   CALL error_alloc('deptrans',1,(/numpoints/),kndrtype)
   
!
!4.3 Coordinate arrays along transect
!------------------------------------
!

   gxstart = postatts(ilist)%gxstart
   gxend = postatts(ilist)%gxend
   gystart = postatts(ilist)%gystart
   gyend = postatts(ilist)%gyend
   hlen_unit = postatts(ilist)%hlen_unit

   SELECT CASE (itrtype(ilist))
      CASE (1)
         itrans = (/(i,i=istart,iend)/)
         jtrans = jstart
         dtrans(ilist,1:numpoints)  = (/(gxcoord(i,jstart),i=istart,iend)/)
         trangle(ilist) = 0.0
         dmin(ilist) = gxstart; dmax(ilist) = gxend
      CASE (2)
         itrans = istart
         jtrans = (/(j,j=jstart,jend)/)
         dtrans(ilist,1:numpoints)  = (/(gycoord(istart,j),j=jstart,jend)/)
         trangle(ilist) = halfpi
         dmin(ilist) = gystart; dmax(ilist) = gyend 
      CASE (3)
         deli = (iend-istart)/REAL(numpoints-1)
         delj = (jend-jstart)/REAL(numpoints-1)
         itrans = (/(NINT(istart+(n-1)*deli),n=1,numpoints)/)
         jtrans = (/(NINT(jstart+(n-1)*delj),n=1,numpoints)/)
         xtrans = (/(gxcoord(itrans(n),jtrans(n)),n=1,numpoints)/)
         ytrans = (/(gycoord(itrans(n),jtrans(n)),n=1,numpoints)/)
         dtrans(ilist,1) = 0.0
         n_431: DO n=2,numpoints
            dels = distance_pts(xtrans(n-1),xtrans(n),&
                              & ytrans(n-1),ytrans(n),2,3)
            dtrans(ilist,n) = dtrans(ilist,n-1) + dels
         ENDDO n_431
         trangle(ilist) = ATAN2(gyend-gystart,gxend-gxstart)
         dmin(ilist) = MINVAL(dtrans(ilist,1:numpoints))
         dmax(ilist) = MAXVAL(dtrans(ilist,1:numpoints))
   END SELECT
   deptrans = (/(depout(itrans(n),jtrans(n)),n=1,numpoints)/)
   IF (iopt_grid_sph.EQ.0) THEN
      dtrans(ilist,1:numpoints) = hlen_unit*dtrans(ilist,1:numpoints)
      dmin(ilist) = hlen_unit*dmin(ilist)
      dmax(ilist) = hlen_unit*dmax(ilist)
   ENDIF
   zmax(ilist) = MAXVAL(deptrans)
   zmin(ilist)= 0.0

!   
!4.4 Write bottom line
!---------------------
!   
!  ---create file
   WRITE (cfig,'(I12)') noplots+ilist; cfig = ADJUSTL(cfig)
   linefile = TRIM(plotfile)//'.'//'PVcdat'//TRIM(cfig)
   CALL open_file(iunit,linefile,'OUT','A')

!  ---write data
   n_440: DO n=1,numpoints
      IF (n.EQ.1) THEN
         WRITE (iunit,9001) dtrans(ilist,1), -deptrans(n), 0
      ELSE
         WRITE (iunit,9001) dtrans(ilist,n), -deptrans(n), 1
      ENDIF
   ENDDO n_440

!  ---close file
   CALL close_file(iunit,'(A)')

!
!4.5 Deallocate
!---------------
!

   DEALLOCATE (itrans,jtrans,xtrans,ytrans,deptrans)
   
ENDDO ilist_400

!
!5. Time loop
!------------
!

norec_500: DO norec=tlims(1),tlims(2)

!
!5.1 Read data
!-------------
!

   flag = loop_index(tlims,norec)
   IF (flag.OR.outform.NE.'N') THEN
      CALL post_read(norec,ctime,rtime,.TRUE.)
   ENDIF
   IF (.NOT.flag) CYCLE norec_500

!
!5.2 Process data
!----------------
!

   ilist_520: DO ilist=1,nolists

!
!5.2.1 Initialise parameters
!---------------------------
!      

      gxstart = postatts(ilist)%gxstart
      gxend = postatts(ilist)%gxend
      gystart = postatts(ilist)%gystart
      gyend = postatts(ilist)%gyend
      istart = postatts(ilist)%istart
      jstart = postatts(ilist)%jstart
      gyend = postatts(ilist)%jend
      horz_unit = postatts(ilist)%horz_unit
      hlen_unit = postatts(ilist)%hlen_unit
      notracks = postatts(ilist)%notracks

!      
!5.2.2 First time 
!----------------
!

      IF (norec.EQ.tlims(1)) THEN
!        ---initialise
         first(ilist,:) = .TRUE.
!        ---create trajectory files
         WRITE (cfig,'(I12)') noplots+ilist; cfig = ADJUSTL(cfig)
         n_522: DO n=1,notracks
            p = postatts(ilist)%ptrack(n)
            WRITE (cp,'(I12)') p; cp = ADJUSTL(cp)
            trackfiles(ilist,n) = TRIM(plotfile)//'.'//TRIM(cp)//'PVcdat'//&
                                & TRIM(cfig)
            CALL open_file(iodat(ilist,n),trackfiles(ilist,n),'OUT','A')
            symbfiles(ilist,n) = TRIM(plotfile)//'.'//TRIM(cp)//'PVbdat'//&
                               & TRIM(cfig)
            CALL open_file(iosymb(ilist,n),symbfiles(ilist,n),'OUT','A')
            WRITE (iosymb(ilist,n),'(A)') '#'
            WRITE (iosymb(ilist,n),'(///)')
            WRITE (iosymb(ilist,n),'(A)') '#'
         ENDDO n_522
      ENDIF
      
!      
!5.2.3 Write trajectory data
!---------------------------
!      

      p_523: DO n=1,notracks
         p = postatts(ilist)%ptrack(n)
         x = outvals2d(p,1); y = outvals2d(p,2); z = outvals2d(p,3)
         flag = ABS(x-fill_value).GT.fill_value_eps 
         IF (flag) THEN
            SELECT CASE (itrtype(ilist))
               CASE (1); s = x
               CASE (2); s = y
               CASE (3)
                  dist = distance_pts(gxstart,gystart,x,y,2,3)
                  phi = ATAN2(y-gystart,x-gxstart)
                  s = dist*COS(trangle(ilist)-phi)
            END SELECT
            IF (iopt_grid_sph.EQ.0) s = hlen_unit*s
            IF (s.GE.dmin(ilist).AND.s.LE.dmax(ilist)) THEN
               zmin(ilist) = MAX(zmin(ilist),z)
               IF (first(ilist,n)) THEN
                  WRITE (iodat(ilist,n),9001) s, z, 0, 1.0, &
                       & postatts(ilist)%linestyles(n)
                  WRITE (iosymb(ilist,n),9001) s, z, 14
                  first(ilist,n) = .FALSE.
               ELSE
                  WRITE (iodat(ilist,n),9001) s, z, 1
                  IF (norec.EQ.tlims(2)) THEN
                     WRITE (iosymb(ilist,n),9001) s, z, 15
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO p_523

!
!5.2.4 Close trajectort files
!----------------------------
!      

           
      IF (norec.EQ.tlims(2)) THEN
         n_524: DO n=1,notracks
            CALL close_file(iodat(ilist,n),'A')
            CALL close_file(iosymb(ilist,n),'A')
         ENDDO n_524
      ENDIF
      
!
!5.2.5 Write parameter file
!----------------------------
!

      IF (norec.EQ.tlims(2)) THEN

!        ---plot borders         
         z1 = zmin(ilist); z2 = zmax(ilist)
         xborder = border*(dtrans(ilist,numpoints)-dtrans(ilist,1))
         yborder = border*(z2+z1)
         plotcorners(1) = dtrans(ilist,1) - xborder
         plotcorners(2) = -z2 - yborder
         plotcorners(3) = dtrans(ilist,numpoints) + xborder
         plotcorners(4) = yborder + z1
         aspect = (dtrans(ilist,numpoints)-dtrans(ilist,1))/(z2+z1)

!        ---main title
         plottitles(1) = TRIM(outfile)

!        ---first subtitle
         plottitles(2) = TRIM(coordvars(3)%long_name)

!        ---subtitles
         SELECT CASE (itrtype(ilist))
            CASE (1)
               IF (iopt_grid_sph.EQ.0.OR.horz_unit.EQ.1) THEN
                  WRITE (cpos(1),'(I12)') jstart; cpos(1) = ADJUSTL(cpos(1))
                  cloc = 'Transect at j = '//TRIM(cpos(1))
               ELSEIF (iopt_grid_sph.EQ.1) THEN
                  geoloc(1) = convert_loc_to_char(gxstart,gystart)
                  cloc = 'Transect at '//geoloc(1)(1:9)
               ENDIF
            CASE (2)
               IF (iopt_grid_sph.EQ.0.OR.horz_unit.EQ.1) THEN
                  WRITE (cpos(1),'(I12)') istart; cpos(1) = ADJUSTL(cpos(1))
                  cloc = 'Transect at i = '//TRIM(cpos(1))
               ELSEIF (iopt_grid_sph.EQ.1) THEN
                  geoloc(1) = convert_loc_to_char(gxstart,gystart)
                  cloc = 'Transect at '//geoloc(1)(12:21)
               ENDIF
            CASE (3)
               IF (iopt_grid_sph.EQ.0.OR.horz_unit.EQ.1) THEN
                  WRITE (cpos(1),'(I12)') istart; cpos(1) = ADJUSTL(cpos(1))
                  WRITE (cpos(2),'(I12)') jstart; cpos(2) = ADJUSTL(cpos(2))
                  WRITE (cpos(3),'(I12)') iend; cpos(3) = ADJUSTL(cpos(3))
                  WRITE (cpos(4),'(I12)') jend; cpos(4) = ADJUSTL(cpos(4))
                  cloc = 'From ('//TRIM(cpos(1))//','//TRIM(cpos(2))//') to ('&
                      & //TRIM(cpos(3))//','//TRIM(cpos(4))//')'
               ELSEIF (iopt_grid_sph.EQ.1) THEN
                  geoloc(1) = convert_loc_to_char(gxstart,gystart)
                  geoloc(2) = convert_loc_to_char(gxend,gyend)
                  cloc = 'From '//TRIM(geoloc(1))//' to '//TRIM(geoloc(2))
               ENDIF
         END SELECT
         plottitles(3) =  TRIM(cloc)

!       ---parameter file
         n_525: DO n=1,notracks
!           --create
            p = postatts(ilist)%ptrack(n)
            WRITE (cp,'(I12)') p; cp = ADJUSTL(cp)
            WRITE (cfig,'(I12)') noplots + ilist; cfig = ADJUSTL(cfig)
            paramfile = TRIM(plotfile)//'.PVpar'//TRIM(cfig)
            CALL open_file(iunit,paramfile,'OUT','A')
!           --plot title
            plottitles(2) = 'Trajectory particle '//TRIM(cp)               
!           --write
            CALL write_mumap_par(iunit,.TRUE.,iopt_grid_sph,aspect,3,&
                            & (/'c','b','c'/),(/trackfiles(ilist,n),&
                            & symbfiles(ilist,n),linefile/))
!           --close parameter file
            CALL close_file(iunit,'A')
!           -- write to '.vis'-file
            WRITE (io_parvis,9002) 'mumap2', TRIM(paramfile), 'y', &
                                 & TRIM(plottitles(2))
         ENDDO n_525

      ENDIF

   ENDDO ilist_520
   
ENDDO norec_500


!
!6. Update plot number
!---------------------
!

noplots = noplots + nolists

!
!7. Deallocate arrays
!--------------------
!

DEALLOCATE (dtrans,symbfiles,trackfiles,iodat,iosymb,first)

CALL log_timer_out()


RETURN

9001 FORMAT (2(G15.7,1X),I2,1X,G15.7,1X,I2)
9002 FORMAT (4(A,1X))

END SUBROUTINE PV_plot

!========================================================================

SUBROUTINE read_line(iunit,nvars,itype,cval,lval,ival,rval,names)
!************************************************************************
!
! *read_line* read data values from input line
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - post_program
! 
! Module calls - error_abort
!
!************************************************************************
!
USE iopars
USE plotpars
USE switches
USE error_routines, ONLY: nerrs, error_abort

IMPLICIT NONE

!
!* Arguments
!
INTEGER, INTENT(IN) :: iunit, nvars
INTEGER, INTENT(IN), DIMENSION(nvars) :: itype
CHARACTER (LEN=120), INTENT(OUT), DIMENSION(Maxpars) :: cval
LOGICAL, INTENT(OUT), DIMENSION(MaxPars) :: lval
INTEGER, INTENT(OUT), DIMENSION(MaxPars) :: ival
REAL, INTENT(OUT), DIMENSION(MaxPars) :: rval
CHARACTER (LEN=*), INTENT(IN) :: names

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*iunit*   INTEGER File unit
!*nvars*   INTEGER Number of parameters on input line
!*itype*   INTEGER Data type of input parameters
!*cval*    CHAR    Character input parameters
!*lval*    LOGICAL Logical input parameters
!*ival*    INTEGER Integer input parameters
!*rval*    REAL    Real input parameters
!*names*   CHAR    String with names of parameters
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cnum
CHARACTER (LEN=600) :: cline
LOGICAL :: lblank
INTEGER :: indc, indi, indl, indr, ivar, l, lpos
CHARACTER (LEN=120), DIMENSION(nvars) :: cdat


pglev = pglev + 1

!
!1. Read until first data line
!-----------------------------
!

cline = '#'
DO WHILE (cline(1:1).EQ.'#')
   plinenum = plinenum + 1
   WRITE (cnum,'(I12)') plinenum; cnum = ADJUSTL(cnum)
   READ (iunit,'(A)',END=1001) cline
END DO

!
!2. Echo data names on log file
!------------------------------
!

IF (loglev1.GT.0) THEN
    WRITE (iolog,'(A)') TRIM(names)
    WRITE (iolog,'(A)') TRIM(cline)
ENDIF

!
!3. Extract data as strings
!--------------------------
!

lblank = .TRUE.; ivar = 0
l_310: DO l=1,LEN_TRIM(cline)
   IF (cline(l:l).EQ.' ') THEN
      lblank = .TRUE.
   ELSEIF (lblank) THEN
      lblank = .FALSE.
      ivar = ivar + 1; cdat(ivar) = ' '; lpos = 1
      cdat(ivar)(lpos:lpos) = cline(l:l)
   ELSE
      lpos = lpos + 1
      cdat(ivar)(lpos:lpos) = cline(l:l)
   ENDIF
ENDDO l_310
IF (ivar.NE.nvars) GOTO 1002

!
!4. Convert strings to data values
!---------------------------------
!

indc = 0; indl = 0; indi = 0; indr = 0; cval = ' '
ivar_410: DO ivar=1,nvars
   SELECT CASE(itype(ivar))
      CASE(1)
         indc = indc + 1
         cval(indc) = cdat(ivar)
      CASE(2)
         indl = indl + 1
         READ (cdat(ivar),'(L1)',ERR=1003) lval(indl)
      CASE(3)
         indi = indi + 1
         READ (cdat(ivar),*,ERR=1003) ival(indi)
      CASE(4)
         indr = indr + 1
         READ (cdat(ivar),*,ERR=1003) rval(indr)
   END SELECT
ENDDO ivar_410

pglev = pglev - 1


RETURN

1001 nerrs = 1
WRITE (0,'(A)') 'End of file condition occurred on line '//TRIM(cnum)//&
              & ' of file '//TRIM(postfile)
CALL error_abort('read_line',ierrno_fend)
1002 nerrs = 1
WRITE (0,'(A,I2)') 'Invalid number of parameters on line '//TRIM(cnum)//&
              & ' of file '//TRIM(postfile)//': ', ivar
WRITE (0,'(A,I2)') '  Should be: ', nvars
CALL error_abort('read_line',ierrno_inival)
1003 nerrs = 1
WRITE (0,'(A)') 'Invalid value for parameter on line '//TRIM(cnum)//&
              & ' of file '//TRIM(postfile)
WRITE (0,'(A)') 'Reads: '//TRIM(cdat(ivar))
CALL error_abort('read_line',ierrno_inival)


RETURN

END SUBROUTINE read_line

!========================================================================

SUBROUTINE TH_plot
!************************************************************************
!
! *TH_plot* Time-distance contour plots along transect
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - post_program
!
! External calls - coords_to_index, post_read, time_limits, vert_retrieve_1d,
!                  write_mumap_cont, write_mumap_par
!
! Module calls - close_file, convert_date, convert_loc_to_char, diff_dates,
!                error_alloc, lim_dims, loop_index, open_file
!
!************************************************************************
!
USE grid
USE iopars
USE physpars
USE plotpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: convert_loc_to_char
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: convert_date, diff_dates, log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims, loop_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=12) :: cfig
CHARACTER (LEN=lendesc) :: cloc
CHARACTER (LEN=leniofile) :: paramfile
CHARACTER (LEN=lentime) :: ctime
INTEGER :: horz_unit, i, iend, ilist, istart, iunit, ivar, j, jend, jstart, &
         & maxpoints, n, &
         & nhdat, nn, nolevels, norec, notimes, novarskiph, ntout, numpoints
REAL :: aspect, contmax, contmin, deli, delj, dels, fac, gxend, gxstart, &
      & gyend, gystart, valmax, valmin, xborder, yborder
REAL (KIND=kndrlong) :: rtime

CHARACTER (LEN=12), DIMENSION(4) :: cpos
CHARACTER (LEN=21), DIMENSION(2) :: geoloc
CHARACTER (LEN=leniofile), DIMENSION(nolists) :: contfiles
INTEGER, DIMENSION(3) :: tlims
INTEGER, DIMENSION(nolists) :: itrtype, nopoints
REAL, DIMENSION(MaxContLevs) :: contlevels

INTEGER, ALLOCATABLE, DIMENSION(:,:) :: itrans, jtrans
REAL, ALLOCATABLE, DIMENSION(:) :: timedat
REAL, ALLOCATABLE, DIMENSION(:,:) :: dtrans, htrans, trangle, ttrans, xtrans, &
                                   & ytrans
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: conttrans


procname(pglev+1) = 'TH_plot'
CALL log_timer_in()

!
!1. Initialise time parameters and arrays
!----------------------------------------
!
!---output time limiters
CALL time_limits(tlims)

!---offset
CALL diff_dates(outgpars%refdate,PRefDateTime,ptime_format,rtime=offset)

!---number of time steps
notimes = lim_dims(tlims)

!
!2. Define transects
!-------------------
!
!2.1 Horizontal locations
!------------------------
!

ilist_210: DO ilist=1,nolists
   
   IF (postatts(ilist)%horz_unit.EQ.1) THEN
      istart = NINT(postatts(ilist)%gxstart)
      iend = NINT(postatts(ilist)%gxend)
      jstart = NINT(postatts(ilist)%gystart)
      jend = NINT(postatts(ilist)%gyend)
      postatts(ilist)%gxstart = gxcoord(istart,jstart)
      postatts(ilist)%gxend = gxcoord(iend,jend)
      postatts(ilist)%gystart = gycoord(istart,jstart)
      postatts(ilist)%gyend = gycoord(iend,jend)
      postatts(ilist)%istart = istart
      postatts(ilist)%iend = iend
      postatts(ilist)%jstart = jstart
      postatts(ilist)%jend = jend
   ELSE
      CALL coords_to_index(postatts(ilist)%gxstart,postatts(ilist)%gystart,&
                         & istart,jstart)
      CALL coords_to_index(postatts(ilist)%gxend,postatts(ilist)%gyend,&
                         & iend,jend)
      postatts(ilist)%istart = istart
      postatts(ilist)%iend = iend
      postatts(ilist)%jstart = jstart
      postatts(ilist)%jend = jend
   ENDIF
   
!  ---number of data ponts along transect
   nopoints(ilist) = MAX(ABS(iend-istart),ABS(jend-jstart)) + 1
   
ENDDO ilist_210

!
!2.2 Alllocate arrays
!--------------------
!

maxpoints = MAXVAL(nopoints)
ALLOCATE (timedat(notimes),STAT=errstat)
CALL error_alloc('timedat',1,(/notimes/),kndrtype)
ALLOCATE (conttrans(notimes,maxpoints,nolists),STAT=errstat)
CALL error_alloc('conttrans',3,(/notimes,maxpoints,nolists/),kndrtype)
conttrans = 0.0
ALLOCATE (itrans(maxpoints,nolists),STAT=errstat)
CALL error_alloc('itrans',2,(/maxpoints,nolists/),kndint)
ALLOCATE (jtrans(maxpoints,nolists),STAT=errstat)
CALL error_alloc('jtrans',2,(/maxpoints,nolists/),kndint)
ALLOCATE (htrans(maxpoints,nolists),STAT=errstat)
CALL error_alloc('htrans',2,(/maxpoints,nolists/),kndrtype)
ALLOCATE (trangle(maxpoints,nolists),STAT=errstat)
CALL error_alloc('trangle',2,(/maxpoints,nolists/),kndrtype)
ALLOCATE (xtrans(maxpoints,nolists),STAT=errstat)
CALL error_alloc('xtrans',2,(/maxpoints,nolists/),kndrtype)
ALLOCATE (ytrans(maxpoints,nolists),STAT=errstat)
CALL error_alloc('ytrans',2,(/maxpoints,nolists/),kndrtype)

!
!2.3 Define transect parameters
!------------------------------
!

ilist_230: DO ilist=1,nolists

!  ---initialise parameters
   istart = postatts(ilist)%istart; iend = postatts(ilist)%iend
   jstart = postatts(ilist)%jstart; jend = postatts(ilist)%jend
   numpoints = nopoints(ilist)

!  ---type of transect
   IF (jstart.EQ.jend) THEN
      itrtype(ilist) = 1
   ELSEIF (istart.EQ.iend) THEN
      itrtype(ilist) = 2
   ELSE
      itrtype(ilist) = 3
   ENDIF

!  ---coordinate arrays along transect
   SELECT CASE (itrtype(ilist))
   CASE (1)
      itrans(1:numpoints,ilist) = (/(i,i=istart,iend)/)
      jtrans(1:numpoints,ilist) = jstart
      htrans(1:numpoints,ilist) = (/(gxcoord(i,jstart),i=istart,iend)/)
      trangle(1:numpoints,ilist) = 0.0
   CASE (2)
      itrans(1:numpoints,ilist) = istart
      jtrans(1:numpoints,ilist) = (/(j,j=jstart,jend)/)
      htrans(1:numpoints,ilist) = (/(gycoord(istart,j),j=jstart,jend)/)
      trangle(1:numpoints,ilist) = halfpi
   CASE (3)
      deli = (iend-istart)/REAL(numpoints-1)
      delj = (jend-jstart)/REAL(numpoints-1)
      itrans(1:numpoints,ilist) = (/(NINT(istart+(n-1)*deli),n=1,numpoints)/)
      jtrans(1:numpoints,ilist) = (/(NINT(jstart+(n-1)*delj),n=1,numpoints)/)
      xtrans(1:numpoints,ilist) = (/(gxcoord(itrans(n,ilist),jtrans(n,ilist)),&
                                   & n=1,numpoints)/)
      ytrans(1:numpoints,ilist) = (/(gycoord(itrans(n,ilist),jtrans(n,ilist)),&
                                   & n=1,numpoints)/)
      htrans(1,ilist) = 0.0
      n_231: DO n=2,numpoints
         IF (iopt_grid_sph.EQ.0) THEN
            dels = SQRT((xtrans(n,ilist)-xtrans(n-1,ilist))**2+&
                      & (ytrans(n,ilist)-ytrans(n-1,ilist))**2)
         ELSE
            fac = COS(0.5*degtorad*(ytrans(n-1,ilist)+ytrans(n,ilist)))
            dels = Rearth*degtorad*&
                 & SQRT((fac*(xtrans(n,ilist)-xtrans(n-1,ilist)))**2+&
                           & (ytrans(n,ilist)-ytrans(n-1,ilist))**2)
         ENDIF
         htrans(n,ilist) = htrans(n-1,ilist) + dels
         trangle(n-1,ilist) = ATAN2(ytrans(n,ilist)-ytrans(n-1,ilist),&
                                  & xtrans(n,ilist)-xtrans(n-1,ilist))
      ENDDO n_231
      trangle(numpoints,ilist) = trangle(numpoints-1,ilist)
   END SELECT
   IF (iopt_grid_sph.EQ.0.OR.itrtype(ilist).EQ.3) THEN
      htrans(1:numpoints,ilist) = postatts(ilist)%hlen_unit*&
                                & htrans(1:numpoints,ilist)
   ENDIF

ENDDO ilist_230

!
!3. Skip data reads
!------------------
!

IF (outform.NE.'N') THEN
   norec_310: DO norec=1,tlims(1)-1
      CALL post_read(norec,ctime,rtime,fullinfo)
   ENDDO norec_310
ENDIF

!
!4. Time loop
!------------
!

ntout = 0
norec_400: DO norec=tlims(1),tlims(2)

!
!4.1 Read data
!-------------
!

   flag = loop_index(tlims,norec)
   IF (flag.OR.outform.NE.'N') THEN
      CALL post_read(norec,ctime,rtime,fullinfo)
   ENDIF
   IF (.NOT.flag) CYCLE norec_400
   ntout = ntout + 1
!  ---store time
   timedat(ntout) = rtime

!
!4.2 Process data
!----------------
!

   ilist_420: DO ilist=1,nolists

!
!4.2.1 Initialise parameters
!---------------------------
!

      ivar = postatts(ilist)%numvar
      novarskiph = postatts(ilist)%novarskip1
      numpoints = nopoints(ilist)

!
!4.2.2 Evaluate and store data
!-----------------------------
!

      nn = 0
      n_422: DO n=1,numpoints,novarskiph
         nn = nn + 1
         i = itrans(n,ilist); j = jtrans(n,ilist)
         IF (nodim.EQ.2) THEN
            conttrans(ntout,nn,ilist) = outvals3d(i,j,ivar)
         ELSEIF (nodim.EQ.3) THEN
            CALL vert_retrieve_1d(conttrans(ntout,nn,ilist),&
                                & outvals4d(i,j,:,ivar),i,j,ilist)
         ENDIF
      ENDDO n_422

   ENDDO ilist_420

ENDDO norec_400

!
!5. Write contour data
!---------------------
!

ilist_500: DO ilist=1,nolists

!
!5.1 Initialise parameters
!-------------------------
!

   contmin = postatts(ilist)%contmin
   contmax = postatts(ilist)%contmax
   nolevels = postatts(ilist)%nolevels
   novarskiph = postatts(ilist)%novarskip1
   numpoints = nopoints(ilist)
   nhdat = (numpoints-1)/novarskiph + 1

!            
!5.2 Contour levels
!------------------
!

   IF (postatts(ilist)%contreg) THEN
      IF (contmin.EQ.0.0.AND.contmax.EQ.0.0) THEN
         valmin = MINVAL(conttrans(:,1:nhdat,ilist),&
            & MASK=ABS(conttrans(:,1:nhdat,ilist)-fill_value).GT.fill_value_eps)
         valmax = MAXVAL(conttrans(:,1:nhdat,ilist),&
            & MASK=ABS(conttrans(:,1:nhdat,ilist)-fill_value).GT.fill_value_eps)
      ELSE
         valmin = contmin; valmax = contmax
      ENDIF
      contlevels(1:nolevels) = (/(valmin+(n-1)*(valmax-valmin)/(nolevels-1),&
                              & n=1,nolevels)/)
   ELSE
      contlevels(1:nolevels) = postatts(ilist)%contlevels(1:nolevels)
   ENDIF

!
!5.3 Dummy arrays
!----------------
!
!  ---allocate
   ALLOCATE (ttrans(notimes,nhdat),STAT=errstat)
   CALL error_alloc('ttrans',2,(/notimes,nhdat/),kndrtype)
   ALLOCATE (dtrans(notimes,nhdat),STAT=errstat)
   CALL error_alloc('dtrans',2,(/notimes,nhdat/),kndrtype)

!  ---evaluate
   ntout_530: DO ntout=1,notimes
      ttrans(ntout,:) = timedat(ntout)
      IF (ntout.EQ.1) THEN
         dtrans(1,:) = htrans(1:numpoints:novarskiph,ilist)
      ELSE
         dtrans(ntout,:) = dtrans(1,:)
      ENDIF
   ENDDO ntout_530

!
!5.4 Create data file
!--------------------
!

   WRITE (cfig,'(I12)') noplots + ilist; cfig = ADJUSTL(cfig)
   contfiles(ilist) = TRIM(plotfile)//'.THidat'//TRIM(cfig)
   CALL open_file(iunit,contfiles(ilist),'OUT','A')

!
!5.5 Write data
!--------------
!

   CALL write_mumap_cont(iunit,notimes,nhdat,ttrans,dtrans,&
                       & conttrans(:,1:nhdat,ilist),nolevels,&
                       & contlevels(1:nolevels),postatts(ilist)%icontstyle,0,2)
   CALL close_file(iunit,'A')

!
!5.6 Deallocate
!--------------
!

   DEALLOCATE (ttrans,dtrans)

ENDDO ilist_500

!
!6. Write parameter file
!------------------------
!

ilist_600: DO ilist=1,nolists

!
!6.1 Initialise parameters
!-------------------------
!

   ivar = postatts(ilist)%numvar
   horz_unit = postatts(ilist)%horz_unit
   novarskiph = postatts(ilist)%novarskip1
   numpoints = nopoints(ilist)
   nhdat = (numpoints-1)/novarskiph + 1
   istart = postatts(ilist)%istart; iend  = postatts(ilist)%iend
   jstart = postatts(ilist)%jstart; jend  = postatts(ilist)%jend
   gxstart = postatts(ilist)%gxstart; gxend = postatts(ilist)%gxend
   gystart = postatts(ilist)%gystart; gyend = postatts(ilist)%gyend
   WRITE (cfig,'(I12)') noplots + ilist; cfig = ADJUSTL(cfig)   

!
!6.2 Create and open parameter file
!----------------------------------
!

   paramfile = TRIM(plotfile)//'.THpar'//TRIM(cfig)
   CALL open_file(iunit,paramfile,'OUT','A')

!
!6.3 Plot borders
!----------------
!

   xborder = border*(timedat(notimes)-timedat(1))
   yborder = border*(htrans(nhdat,ilist)-htrans(1,ilist))
   plotcorners(1) = timedat(1) - xborder
   plotcorners(2) = htrans(1,ilist) - yborder
   plotcorners(3) = timedat(notimes) + xborder
   plotcorners(4) = htrans(nhdat,ilist) + yborder
   aspect = (timedat(notimes)-timedat(1))/&
           & ABS(htrans(nhdat,ilist)-htrans(1,ilist))

!
!6.4 Titles
!----------
!
!  ---main title
   plottitles(1) = TRIM(outfile)

!  ---first subtitle
   plottitles(2) = TRIM(outvars(ivar)%long_name)//' ('//&
                 & TRIM(outvars(ivar)%units)//')'

!  ---second subtitle
   SELECT CASE (itrtype(ilist))
   CASE (1)
      IF (iopt_grid_sph.EQ.0.OR.horz_unit.EQ.1) THEN
         WRITE (cpos(1),'(I12)') jstart; cpos(1) = ADJUSTL(cpos(1))
         cloc = 'Transect at j = '//TRIM(cpos(1))
      ELSEIF (iopt_grid_sph.EQ.1) THEN
         geoloc(1) = convert_loc_to_char(gxstart,gystart)
         cloc = 'Transect at '//geoloc(1)(1:9)
      ENDIF
   CASE (2)
      IF (iopt_grid_sph.EQ.0.OR.horz_unit.EQ.1) THEN
         WRITE (cpos(1),'(I12)') istart; cpos(1) = ADJUSTL(cpos(1))
         cloc = 'Transect at i = '//TRIM(cpos(1))
      ELSEIF (iopt_grid_sph.EQ.1) THEN
         geoloc(1) = convert_loc_to_char(gxstart,gystart)
         cloc = 'Transect at '//geoloc(1)(12:21)
      ENDIF
   CASE (3)
      IF (iopt_grid_sph.EQ.0.OR.horz_unit.EQ.1) THEN
         WRITE (cpos(1),'(I12)') istart; cpos(1) = ADJUSTL(cpos(1))
         WRITE (cpos(2),'(I12)') jstart; cpos(2) = ADJUSTL(cpos(2))
         WRITE (cpos(3),'(I12)') iend; cpos(3) = ADJUSTL(cpos(3))
         WRITE (cpos(4),'(I12)') jend; cpos(4) = ADJUSTL(cpos(4))
         cloc = 'From ('//TRIM(cpos(1))//','//TRIM(cpos(2))//' to ('//&
              & TRIM(cpos(3))//','//TRIM(cpos(4))//')'
      ELSEIF (iopt_grid_sph.EQ.1) THEN
         geoloc(1) = convert_loc_to_char(gxstart,gystart)
         geoloc(2) = convert_loc_to_char(gxend,gyend)
         cloc = 'From '//TRIM(geoloc(1))//' to '//TRIM(geoloc(2))
      ENDIF
   END SELECT
   plottitles(3) = TRIM(cloc)

!  ---write parameter file
   CALL write_mumap_par(iunit,.TRUE.,0,aspect,1,(/'i'/),(/contfiles(ilist)/))
   CALL close_file(iunit,'A')

!
!6.5 Write to '.vis'-file
!------------------------
!

   WRITE (io_parvis,9001) 'mumap2', TRIM(paramfile), 'y', &
                        & TRIM(outvars(ivar)%f90_name)

ENDDO ilist_600

!
!7. Deallocate arrays
!----------------------
!

DEALLOCATE (timedat,conttrans,itrans,jtrans,htrans,trangle,xtrans,ytrans)

!
!8. Update figure number
!-----------------------
!

noplots = noplots + nolists

CALL log_timer_out()


RETURN

9001 FORMAT (4(A,1X))

END SUBROUTINE TH_plot

!========================================================================

SUBROUTINE time_limits(tlims)
!************************************************************************
!
! *time_limits* start/end end increment record number for time series plots 
!
! Author - Patrick Luyten
!
! Description -
!
! Reference -
!
! Calling program - HT_plot, TH_plot, TS_plot, TZ_plot, VT_plot
!
! External calls - time_locate
!
! Module calls - check_time_limits, error_abort
!
!************************************************************************
!
USE iopars
USE plotpars
USE timepars
USE error_routines, ONLY: check_time_limits, error_abort
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(OUT), DIMENSION(3) :: tlims

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*tlims*     INTEGER Start/end/increment time record numbers
!
!------------------------------------------------------------------------------
!

procname(pglev+1) = 'time_limits'
CALL log_timer_in()

!
!1. Evaluate limiting times
!--------------------------
!

CALL time_locate(PStartDateTime,tlims(1))
CALL time_locate(PEndDateTime,tlims(2))
tlims(3) = norecskip

!
!2. Check
!--------
!

CALL check_time_limits(tlims,'tlims',0,nstepout)
CALL error_abort('time_limits',ierrno_runval)

CALL log_timer_out()


RETURN

END SUBROUTINE time_limits

!========================================================================

SUBROUTINE time_locate(ctime,norec)
!************************************************************************
!
! *time_locate* Return the record number within a time series output data set
!               which is closest to the input date/time 'ctime'
!
! Author - Patrick Luyten
!
! Description - 
!
! Reference -
!
! Calling program - HP_plot, time_limits, ZP_plot
!
! Module calls - add_secs_to_date, diff_dates
!
!************************************************************************
!
USE iopars
USE plotpars
USE timepars
USE syspars
USE time_routines, ONLY: add_secs_to_date, diff_dates, log_timer_in, &
                       & log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ctime
INTEGER, INTENT(OUT) :: norec

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*ctime*   CHAR    Input date/time
!*norec*   INTEGER Time record number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL (KIND=kndrlong) :: rtime


procname(pglev+1) = 'time_locate'
CALL log_timer_in()

!
!1. One time output only
!-----------------------
!

IF (nstepout.EQ.1) THEN
   norec = 1
   GOTO 1000
ENDIF

!
!2. Regular time output
!----------------------
!

CALL diff_dates(outgpars%startdate,ctime,gtime_format,rtime=rtime)
norec = INT(ABS(rtime/deltout)) + 1
IF ((norec.LT.nstepout).AND.(ABS(rtime).GT.ABS(deltout)*(norec-0.5))) THEN
   norec = norec + 1
   CALL add_secs_to_date(outgpars%startdate,ctime,gtime_format,&
                      & (norec-1)*deltout)
ENDIF

1000 CALL log_timer_out()


RETURN

END SUBROUTINE time_locate

!========================================================================

SUBROUTINE TS_plot
!************************************************************************
!
! *TS_plot* Time series plots
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - post_program
!
! External calls - coords_to_index, plot_time_limitss, post_read,
!                  vert_retrieve_1d, write_muplot_par
!
! Module calls - close_file, convert_date, convert_loc_to_char, diff_dates,
!                error_alloc, loop_index, open_file
!
!************************************************************************
!
USE grid
USE iopars
USE plotpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: convert_loc_to_char
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: convert_date, diff_dates, log_timer_in, &
                       & log_timer_out
USE utility_routines, ONLY: loop_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=12) :: ccurv, cfig, cipos, cjpos, cplot
CHARACTER (LEN=lendesc) :: cdate, cform, cloc
CHARACTER (LEN=lentime) :: ctime
CHARACTER (LEN=leniofile) :: paramfile
INTEGER :: icurv, ifig, ilist, ipos, iunit, ivar, jpos, norec, &
         & numcurves, numstat
REAL :: outdat
REAL, DIMENSION(2) :: rangemax, rangemin
REAL (KIND=kndrlong) :: rtime

CHARACTER (LEN=leniofile), DIMENSION(nolists) :: datfiles
INTEGER, DIMENSION(2) :: irange, istyle
INTEGER, DIMENSION(3) :: tlims 
INTEGER, DIMENSION(nolists) :: iodat
INTEGER, DIMENSION(nofigs) :: nocurves
INTEGER, DIMENSION(nofigs,nolists) :: indcurv
INTEGER, ALLOCATABLE, DIMENSION(:) :: indlist


procname(pglev+1) = 'TS_plot'
CALL log_timer_in()

!
!1. Number of curves and curve index array
!-----------------------------------------
!

ifig_110: DO ifig=1,nofigs
   numcurves = 0
   ilist_111: DO ilist=1,nolists
      IF (postatts(ilist)%numfig.EQ.ifig) THEN
         numcurves = numcurves + 1
         indcurv(ifig,numcurves) = ilist
      ENDIF
   ENDDO ilist_111
   nocurves(ifig) = numcurves
ENDDO ifig_110

!
!2. Horizontal location
!----------------------
!

IF (nodim.GT.0) THEN
   ilist_210: DO ilist=1, nolists
      IF (gridded) THEN
         IF (postatts(ilist)%horz_unit.EQ.1) THEN
            ipos = NINT(postatts(ilist)%gxpos)
            jpos = NINT(postatts(ilist)%gypos)
            postatts(ilist)%gxpos = gxcoord(ipos,jpos)
            postatts(ilist)%gypos = gycoord(ipos,jpos)
            postatts(ilist)%ipos = ipos
            postatts(ilist)%jpos = jpos
         ELSE
            CALL coords_to_index(postatts(ilist)%gxpos,postatts(ilist)%gypos,&
                               & postatts(ilist)%ipos,postatts(ilist)%jpos)
         ENDIF
      ELSE
         numstat = postatts(ilist)%numstat
         postatts(ilist)%gxpos = gxcoord(numstat,1)
         postatts(ilist)%gypos = gycoord(numstat,1)
         postatts(ilist)%ipos = numstat
         postatts(ilist)%jpos = 1
      ENDIF
   ENDDO ilist_210
ELSE
   postatts%ipos = 0; postatts%jpos = 0
ENDIF

!
!3. Time parameters
!------------------
!
!---output time limiters
CALL time_limits(tlims)

!---offset
CALL diff_dates(outgpars%startdate,PStartDateTime,ptime_format,rtime=offset)

!
!4. Create/open data files
!-------------------------
!

ifig_410: DO ifig=1,nofigs
icurv_410: DO icurv=1,nocurves(ifig)
   ilist = indcurv(ifig,icurv)
   WRITE (cfig,'(I12)') noplots+ifig; cfig = ADJUSTL(cfig)
   WRITE (ccurv,'(I12)') icurv; ccurv = ADJUSTL(ccurv)
   datfiles(ilist) = TRIM(plotfile)//'.'//TRIM(ccurv)//'TSdat'//TRIM(cfig)
   CALL open_file(iodat(ilist),datfiles(ilist),'OUT','A')
ENDDO icurv_410
ENDDO ifig_410

!
!5. Skip data reads
!------------------
!

IF (outform.NE.'N') THEN
   norec_510: DO norec=1,tlims(1)-1
      CALL post_read(norec,ctime,rtime,fullinfo)
   ENDDO norec_510
ENDIF

!
!6. Read data, write data files
!------------------------------
!

norec_610: DO norec=tlims(1),tlims(2)

!  ---read data
   flag = loop_index(tlims,norec)
   IF (flag.OR.outform.NE.'N') THEN
      CALL post_read(norec,ctime,rtime,fullinfo)
   ENDIF
   IF (.NOT.flag) CYCLE norec_610
 
!  ---store data for plotting
   ilist_611: DO ilist=1,nolists
      ipos = postatts(ilist)%ipos; jpos = postatts(ilist)%jpos
      ivar = postatts(ilist)%numvar
      SELECT CASE(nodim)
      CASE(0); outdat = outvals1d(ivar)
      CASE(2)
         IF (gridded) THEN
            outdat = outvals3d(ipos,jpos,ivar)
         ELSE
            outdat = outvals2d(ipos,ivar)
         ENDIF
      CASE(3)
         IF (gridded) THEN
            CALL vert_retrieve_1d(outdat,outvals4d(ipos,jpos,:,ivar),ipos,jpos,&
                                & ilist)
         ELSE
            CALL vert_retrieve_1d(outdat,outvals3d(ipos,:,ivar),ipos,1,ilist)
         ENDIF
      END SELECT
      IF (ABS(outdat-fill_value).GT.fill_value_eps) THEN
         WRITE (iodat(ilist),9001) rtime, outdat
      ENDIF
   ENDDO ilist_611

ENDDO norec_610

!
!7. Close data files
!-------------------
!

ilist_710: DO ilist=1,nolists
   CALL close_file(iodat(ilist),'A')
ENDDO ilist_710

!
!8. Write parameter file
!-----------------------
!

ifig_810: DO ifig=1,nofigs

!  ---select curves
   numcurves = nocurves(ifig)
   ALLOCATE (indlist(numcurves),STAT=errstat)
   CALL error_alloc('indlist',1,(/numcurves/),kndint)
   indlist = indcurv(ifig,1:numcurves)
   ilist = indlist(1)

!  ---open data file
   WRITE (cfig,'(I12)') noplots + ifig; cfig = ADJUSTL(cfig)
   paramfile = TRIM(plotfile)//'.TSpar'//TRIM(cfig)
   CALL open_file(iunit,paramfile,'OUT','A')

!  ---axis titles
   plottitles(1) = TRIM(time_unit(ptime_format))
   ivar = postatts(ilist)%numvar
   plottitles(2) = TRIM(outvars(ivar)%long_name)//' ('//&
                 & TRIM(outvars(ivar)%units)//')'

!  ---main title
   plottitles(3) = TRIM(outfile)

!  ---subtitle
   IF (nodim.GT.0) THEN
      ipos = postatts(ilist)%ipos; jpos = postatts(ilist)%jpos
      IF (gridded) THEN
         IF (iopt_grid_sph.EQ.0.OR.postatts(ilist)%horz_unit.EQ.1) THEN
            WRITE (cipos,'(I12)') ipos; cipos = ADJUSTL(cipos)
            WRITE (cjpos,'(I12)') jpos; cjpos = ADJUSTL(cjpos)
            cloc = '('//TRIM(cipos)//','//TRIM(cjpos)//')'
         ELSEIF (iopt_grid_sph.EQ.1) THEN
            cloc(1:21) = convert_loc_to_char(postatts(ilist)%gxpos,&
                                           & postatts(ilist)%gypos)
         ENDIF
      ELSE
         numstat = postatts(ilist)%numstat
         cloc = TRIM(station_names(numstat))
      ENDIF
   ENDIF
   cdate = PStartDateTime(1:19)//' to '//PEndDateTime(1:19)
   IF (nodim.EQ.3) THEN
      SELECT CASE (postatts(ilist)%iplotform)
      CASE(1)
         WRITE (cplot,'(I12)') postatts(ilist)%kpos; cplot = ADJUSTL(cplot)
         cform = 'k = '//TRIM(cplot)
      CASE(2)
         WRITE (cplot,'(F6.2)') postatts(ilist)%depplot
         cplot = ADJUSTL(cplot)
         cform = 'z = '//TRIM(cplot)
      CASE(3); cform = 'depth-averaged'
      CASE(4); cform = 'depth-integrated'
      CASE(5); cform = 'max. value'
      END SELECT
   ENDIF
   SELECT CASE (nodim)
      CASE (0); plottitles(4) = TRIM(cdate)
      CASE (2); plottitles(4) = TRIM(cloc)//'; '//TRIM(cdate)
      CASE (3); plottitles(4) = TRIM(cloc)//'; '//TRIM(cdate)//'; '//&
                              & TRIM(cform)
   END SELECT

!  ---axis styles and ranges
   istyle(1) = 0; irange(1) = 0
   IF (nocurves(ifig).EQ.1) THEN
      istyle(2) = 0; irange(2) = 0
   ELSE
      istyle(2) = 2; irange(2) = 1
   ENDIF
   rangemin = 0.0; rangemax = 0.0

!  ---write parameter file
   CALL write_muplot_par(iunit,numcurves,datfiles(indlist),&
                       & istyle,irange,rangemin,rangemax,&
                       & postatts(indlist)%linepsyms)

!  ---close parameter file
   CALL close_file(iunit,'A')

!  ---deallocate
   DEALLOCATE (indlist)

!  ---write to '.vis'-file
   WRITE (io_parvis,9002) 'muplot', TRIM(paramfile), 'y',&
                         & TRIM(outvars(ivar)%f90_name)

ENDDO ifig_810

noplots = noplots + nofigs

CALL log_timer_out()


RETURN

9001 FORMAT (2(G15.7,1X))
9002 FORMAT (4(A,1X))

END SUBROUTINE TS_plot

!========================================================================

SUBROUTINE TZ_plot
!************************************************************************
!
! *TZ_plot* Time depth contour plots
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - post_program
!
! External calls - coords_to_index, post_read, time_limits, write_mumap_cont,
!                  write_mumap_par
!
! Module calls - close_file, convert_date, convert_loc_to_char, diff_dates,
!                error_alloc, lim_dims, loop_index, open_file
!
!************************************************************************
!
USE depths
USE grid
USE iopars
USE plotpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: convert_loc_to_char
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: convert_date, diff_dates, log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims, loop_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=7) :: cdep
CHARACTER (LEN=12) :: cfig, cipos, cjpos
CHARACTER (LEN=lentime) :: ctime
CHARACTER (LEN=leniofile) :: paramfile
INTEGER :: ilist, ipos, iunit, ivar, jpos, k, kk, k1, k2, n, nolevels, &
         & norec, notimes, nskipz, ntout, numstat, nzdim
REAL :: aspect, contmax, contmin, valmax, valmin, xborder, yborder, &
      & z1, z2
REAL (KIND=kndrlong) :: rtime

CHARACTER (LEN=leniofile), DIMENSION(nolists) :: contfiles, linefiles
INTEGER, DIMENSION(3) :: tlims
INTEGER, DIMENSION(nolists) :: iocont, kmax, kmin
REAL, DIMENSION(nolists) :: depdat, zmax, zmin
REAL, DIMENSION(MaxContLevs) :: contlevels

REAL, ALLOCATABLE, DIMENSION(:,:) :: timedat, zetdat
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: contdat, zdat


procname(pglev+1) = 'TZ_plot'
CALL log_timer_in()

!
!1. Time parameters
!------------------
!
!---output time limiters
CALL time_limits(tlims)

!---offset
CALL diff_dates(outgpars%startdate,PStartDateTime,ptime_format,rtime=offset)

!---number of time steps
notimes = lim_dims(tlims)

!
!2. Allocate arrays
!------------------
!

ALLOCATE (timedat(notimes,nzout),STAT=errstat)
CALL error_alloc('timedat',2,(/notimes,nzout/),kndrtype)
ALLOCATE (contdat(notimes,nzout,nolists),STAT=errstat)
CALL error_alloc('contdat',3,(/notimes,nzout,nolists/),kndrtype)
ALLOCATE (zdat(0:notimes,nzout,nolists),STAT=errstat)
CALL error_alloc('zdat',3,(/notimes+1,nzout,nolists/),kndrtype)
ALLOCATE (zetdat(notimes,nolists),STAT=errstat)
CALL error_alloc('zetdat',2,(/notimes,nolists/),kndrtype)

!
!3. Spatial location and water depth
!-----------------------------------
!

ilist_310: DO ilist=1, nolists

!  ---horizontal location
   IF (gridded) THEN
      IF (postatts(ilist)%horz_unit.EQ.1) THEN
         ipos = NINT(postatts(ilist)%gxpos)
         jpos = NINT(postatts(ilist)%gypos)
         postatts(ilist)%gxpos = gxcoord(ipos,jpos)
         postatts(ilist)%gypos = gycoord(ipos,jpos)
         postatts(ilist)%ipos = ipos
         postatts(ilist)%jpos = jpos
      ELSE
         CALL coords_to_index(postatts(ilist)%gxpos,postatts(ilist)%gypos,&
                            & postatts(ilist)%ipos,postatts(ilist)%jpos)
      ENDIF
   ELSE
      numstat = postatts(ilist)%numstat
      postatts(ilist)%gxpos = gxcoord(numstat,1)
      postatts(ilist)%gypos = gycoord(numstat,1)
      postatts(ilist)%ipos = numstat
      postatts(ilist)%jpos = 1
   ENDIF

!  ---water depth
   ipos = postatts(ilist)%ipos; jpos = postatts(ilist)%jpos
   depdat(ilist) = depout(ipos,jpos)

!  ---vertical location
   IF (postatts(ilist)%vert_unit.EQ.1) THEN
      kmin(ilist) = postatts(ilist)%depmin
      kmax(ilist) = postatts(ilist)%depmax
   ELSE
      kmin(ilist) = 1; kmax(ilist) = nzout
   ENDIF   

ENDDO ilist_310

!
!4. Vertical grid
!----------------
!

ilist_410: DO ilist=1,nolists
   k1 = kmin(ilist); k2 = kmax(ilist)
   nskipz = postatts(ilist)%novarskip1
   nzdim = (k2-k1)/nskipz + 1
   ipos = postatts(ilist)%ipos; jpos = postatts(ilist)%jpos
   k_411: DO k=1,nzdim
      kk = k1 + (k-1)*nskipz
      zdat(0,k,ilist) = gzcoord_ini(ipos,jpos,kk)
   ENDDO k_411
   IF (.NOT.ztime) THEN
      zetdat = 0.0
      ntout_411: DO ntout=1,notimes
         zdat(ntout,1:nzdim,ilist) = zdat(0,1:nzdim,ilist)
      ENDDO ntout_411
   ENDIF
   IF (postatts(ilist)%vert_unit.EQ.1) THEN
      zmax(ilist) = MERGE(depdat(ilist),ABS(zdat(0,1,ilist)),k1.EQ.1)
      zmin(ilist) = MERGE(0.0,ABS(zdat(0,k2,ilist)),k2.EQ.nzout)
   ELSEIF (postatts(ilist)%vert_unit.EQ.2) THEN
      zmax(ilist) = MERGE(depdat(ilist),postatts(ilist)%depmax,&
                        & postatts(ilist)%depmax.EQ.0.0)
      zmin(ilist) = postatts(ilist)%depmin
   ENDIF
ENDDO ilist_410

!
!5. Create/open data files
!-------------------------
!

ilist_510: DO ilist=1,nolists
   WRITE (cfig,'(I12)') noplots + ilist; cfig = ADJUSTL(cfig)
   contfiles(ilist) = TRIM(plotfile)//'.TZidat'//TRIM(cfig)
   CALL open_file(iocont(ilist),contfiles(ilist),'OUT','A')
ENDDO ilist_510

!
!6. Skip data reads
!------------------
!

IF (outform.NE.'N') THEN
   norec_610: DO norec=1,tlims(1)-1
      CALL post_read(norec,ctime,rtime,fullinfo)
   ENDDO norec_610
ENDIF

!
!7. Read and store data for contouring
!-------------------------------------
!

ntout = 0
norec_710: DO norec=tlims(1),tlims(2)

!  ---read data
   flag = loop_index(tlims,norec)
   IF (flag.OR.outform.NE.'N') THEN
      CALL post_read(norec,ctime,rtime,fullinfo)
   ENDIF
   IF (.NOT.flag) CYCLE norec_710
   ntout = ntout + 1

!  ---store time
   timedat(ntout,:) = rtime

!  ---store data
   ilist_711: DO ilist=1,nolists
      ivar = postatts(ilist)%numvar
      k1 = kmin(ilist); k2 = kmax(ilist)
      nskipz = postatts(ilist)%novarskip1; nzdim = (k2-k1)/nskipz + 1
      ipos = postatts(ilist)%ipos; jpos = postatts(ilist)%jpos
      IF (gridded) THEN
         contdat(ntout,1:nzdim,ilist) = outvals4d(ipos,jpos,k1:k2:nskipz,ivar)
      ELSE
         contdat(ntout,1:nzdim,ilist) = outvals3d(ipos,k1:k2:nskipz,ivar)
      ENDIF
      IF (ztime) THEN
         zetdat(ntout,ilist) = zeta(ipos,jpos)
         zdat(ntout,1:nzdim,ilist) = gzcoord(ipos,jpos,k1:k2:nskipz)
      ENDIF
      IF (postatts(ilist)%vert_unit.EQ.2) THEN
         WHERE ((ABS(zdat(0,1:nzdim,ilist)).LT.zmin(ilist)).OR.&
              & (ABS(zdat(0,1:nzdim,ilist)).GT.zmax(ilist)))
            contdat(ntout,1:nzdim,ilist) = fill_value
         END WHERE
      ENDIF
   ENDDO ilist_711

ENDDO norec_710

!
!8. Write contour files
!----------------------
!

ilist_810: DO ilist=1,nolists

!  ---initialise
   nskipz = postatts(ilist)%novarskip1
   k1 = kmin(ilist); k2 = kmax(ilist)
   nzdim = (k2-k1)/nskipz + 1
   contmin = postatts(ilist)%contmin
   contmax = postatts(ilist)%contmax
   nolevels = postatts(ilist)%nolevels

!  ---contour levels
   IF (postatts(ilist)%contreg) THEN
      IF (contmin.EQ.0.0.AND.contmax.EQ.0.0) THEN
         valmin = MINVAL(contdat(:,1:nzdim,ilist),&
              & MASK=ABS(contdat(:,1:nzdim,ilist)-fill_value).GT.fill_value_eps)
         valmax = MAXVAL(contdat(:,1:nzdim,ilist),&
              & MASK=ABS(contdat(:,1:nzdim,ilist)-fill_value).GT.fill_value_eps)
      ELSE
         valmin = contmin; valmax = contmax
      ENDIF
      contlevels(1:nolevels) = &
              &(/(valmin+(n-1)*(valmax-valmin)/(nolevels-1),&
              & n=1,nolevels)/)
   ELSE
      contlevels(1:nolevels) = postatts(ilist)%contlevels(1:nolevels)
   ENDIF

!  ---write data
   CALL write_mumap_cont(iocont(ilist),notimes,nzdim,timedat(:,1:nzdim),&
                       & zdat(1:notimes,1:nzdim,ilist),&
                       & contdat(:,1:nzdim,ilist),nolevels,&
                       & contlevels(1:nolevels),postatts(ilist)%icontstyle,0,2)
   CALL close_file(iocont(ilist),'A')

ENDDO ilist_810

!
!9. Draw bottom/surface lines
!----------------------------
!

ilist_910: DO ilist=1,nolists

   IF (postatts(ilist)%depmax.GT.0.0.AND.ztime) THEN
      linefiles(ilist) = ' '
   ELSE
      WRITE (cfig,'(I12)') noplots+ilist; cfig = ADJUSTL(cfig)
      linefiles(ilist) = TRIM(plotfile)//'.TZcdat'//TRIM(cfig)
      CALL open_file(iunit,linefiles(ilist),'OUT','A')
!     ---bottom line
      IF (postatts(ilist)%depmax.EQ.0.0) THEN
         WRITE (iunit,9001) timedat(1,1), -zmax(ilist), 0
         WRITE (iunit,9001) timedat(notimes,1), -zmax(ilist), 1

      ENDIF
!     ---surface line
      IF (ztime) THEN
         WRITE (iunit,9001) timedat(1,1), zetdat(1,ilist), 0
         ntout_911: DO ntout=2,notimes
            WRITE (iunit,9001) timedat(ntout,1), zetdat(ntout,ilist), 1
         ENDDO ntout_911
      ENDIF
      CALL close_file(iunit,'A')
   ENDIF
         
ENDDO ilist_910

!
!10. Write parameter files
!-------------------------
!

ilist_1010: DO ilist=1,nolists

!
!10.1 Create/open file
!---------------------
!

   WRITE (cfig,'(I12)') noplots+ilist; cfig = ADJUSTL(cfig)
   paramfile = TRIM(plotfile)//'.TZpar'//TRIM(cfig)
   CALL open_file(iunit,paramfile,'OUT','A')

!
!10.2 Plot borders
!-----------------
!

   k1 = kmin(ilist); k2 = kmax(ilist)
   nskipz = postatts(ilist)%novarskip1; nzdim = (k2-k1)/nskipz + 1
   IF (postatts(ilist)%depmin.EQ.0.0) THEN
      z1 = -MAXVAL(zdat(1:notimes,1:nzdim,ilist))
   ELSE
      z1 = zmin(ilist)
   ENDIF
   z2 = zmax(ilist)
   xborder = border*(timedat(notimes,1)-timedat(1,1))
   yborder = border*(z2-z1)
   plotcorners(1) = timedat(1,1)- xborder
   plotcorners(2) = -z2 - yborder
   plotcorners(3) = timedat(notimes,1) + xborder
   plotcorners(4) = yborder - z1
   aspect = (timedat(notimes,1)-timedat(1,1))/(z2-z1)

!
!10.3 Plot titles
!----------------
!
!  ---main title
   plottitles(1) = TRIM(outfile)

!  ---first subtitle
   ivar = postatts(ilist)%numvar
   plottitles(2) = TRIM(outvars(ivar)%long_name)//' ('//&
                 & TRIM(outvars(ivar)%units)//')'

!  ---second subtitle
   ipos = postatts(ilist)%ipos; jpos = postatts(ilist)%jpos
   IF (gridded) THEN
      IF (iopt_grid_sph.EQ.0.OR.postatts(ilist)%horz_unit.EQ.1) THEN
         WRITE (cipos,'(I12)') ipos; cipos = ADJUSTL(cipos)
         WRITE (cjpos,'(I12)') jpos; cjpos = ADJUSTL(cjpos)
         WRITE (cdep,'(F7.2)') depout(ipos,jpos)
         cdep = ADJUSTL(cdep)
         plottitles(3) = '('//TRIM(cipos)//','//TRIM(cjpos)//')'//&
                       & '; depth = '//TRIM(cdep)
      ELSEIF (iopt_grid_sph.EQ.1) THEN
         plottitles(3)(1:21) = convert_loc_to_char(postatts(ilist)%gxpos,&
                                                 & postatts(ilist)%gypos)
         plottitles(3) = TRIM(plottitles(3))//'; depth = '//TRIM(cdep)
      ENDIF
   ELSE
      WRITE (cipos,'(I12)') ipos; cipos = ADJUSTL(cipos)
      WRITE (cdep,'(F7.2)') depout(ipos,1); cdep = ADJUSTL(cdep)
      plottitles(3) = 'Station '//TRIM(cipos)//'; depth = '//TRIM(cdep)
   ENDIF

!
!10.4 Write parameter file
!-------------------------
!

   CALL write_mumap_par(iunit,.TRUE.,0,aspect,2,(/'i','c'/),&
                     & (/contfiles(ilist),linefiles(ilist)/))

!
!10.5 Close parameter file
!-------------------------
!

   CALL close_file(iunit,'A')

!
!10.6 Write to '.vis'-file
!-------------------------
!

   WRITE (io_parvis,9002) 'mumap2', TRIM(paramfile), 'y',&
                         & TRIM(outvars(ivar)%f90_name)

ENDDO ilist_1010

!
!11. Update figure number
!------------------------
!

noplots = noplots + nolists

!
!12. Deallocate arrays
!-- ------------------
!

DEALLOCATE (timedat,contdat,zdat,zetdat)

CALL log_timer_out()


RETURN

9001 FORMAT (2(G15.7,1X),I1)
9002 FORMAT (4(A,1X))

END SUBROUTINE TZ_plot

!========================================================================

SUBROUTINE vert_retrieve_1d(outdat,profdat,i,j,ilist)
!************************************************************************
!
! *vert_retrieve_1d* Select output for plotting from a vertical profile
!
! Author - Patrick Luyten
!
! Description - plot values are selected by specifiers iform, kplot, dplot
!
! Calling program - HP_plot, TH_plot, TS_plot
!
!************************************************************************
!
USE depths
USE grid
USE iopars
USE plotpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
INTEGER, INTENT(IN) :: i, ilist, j
REAL, INTENT(OUT) :: outdat
REAL, INTENT(IN), DIMENSION(nzout) :: profdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    Retrieved output value
!*profdat* REAL    Vertical data profile
!*i*       INTEGER X-index of data point
!*j*       INTEGER Y-index of data point
!*ilist*   INTEGER Plot list number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: iplotform, k, kpos, k1, k2
REAL :: depplot, deptot, ratio
INTEGER, DIMENSION(1) :: klo, kup
REAL, DIMENSION(nzout) :: zdat


procname(pglev+1) = 'vert_retrieve_1d'
CALL log_timer_in()

!
!1. Return at dry location
!-------------------------
!

IF (.NOT.maskatc_int(i,j)) THEN
   outdat = fill_value
   RETURN
ENDIF

!
!2. Initialise
!-------------
!

iplotform = postatts(ilist)%iplotform
kpos = postatts(ilist)%kpos
depplot = postatts(ilist)%depplot
IF (iplotform.GT.1) THEN
   deptot = depout(i,j) + zeta(i,j)
   zdat = MERGE(zeta(i,j)-gzcoord(i,j,:),-gzcoord_ini(i,j,:),ztime)
ENDIF

!
!2. Evaluate according to value of iform
!---------------------------------------
!

SELECT CASE (iplotform)

!
!2.1 At level kplot
!------------------
!

CASE (1)
   outdat = profdat(kpos)

!
!2.2 At depth dplot
!------------------
!

CASE (2)
   IF (deptot.LT.depplot) THEN
      outdat = fill_value
   ELSEIF (zdat(1).LT.depplot) THEN
      outdat = profdat(1)
   ELSEIF (zdat(nzout).GT.depplot) THEN
      outdat = profdat(nzout)
   ELSE
      klo = MINLOC(zdat-depplot,MASK=zdat.GE.depplot)
      kup = MINLOC(depplot-zdat,MASK=zdat.LT.depplot)
      k1 = klo(1); k2 = kup(1)
      ratio = (depplot-zdat(k2))/(zdat(k1)-zdat(k2))
      outdat = ratio*profdat(k1)+(1.0-ratio)*profdat(k2)
   ENDIF

!
!2.3 Depth-averaged
!------------------
!

CASE (3)
   outdat = 0.5*(zdat(nzout)+zdat(nzout-1))*profdat(nzout) +&
               &(deptot-0.5*(zdat(1)+zdat(2)))*profdat(1)
   k_230: DO k=2,nzout-1
      outdat = outdat + 0.5*(zdat(k-1)-zdat(k+1))*profdat(k)
   ENDDO k_230
   outdat = outdat/deptot

!
!2.4 Depth-integrated
!--------------------
!

CASE (4)
   outdat = 0.5*(zdat(nzout)+zdat(nzout-1))*profdat(nzout) +&
               &(deptot-0.5*(zdat(1)+zdat(2)))*profdat(1)
   k_240: DO k=2,nzout-1
      outdat = outdat + 0.5*(zdat(k-1)-zdat(k+1))*profdat(k)
   ENDDO k_240

!
!2.5 Maximum value in the vertical
!---------------------------------
!

CASE (5)
   outdat = MAXVAL(profdat)

END SELECT

CALL log_timer_out()


RETURN

END SUBROUTINE vert_retrieve_1d

!========================================================================

SUBROUTINE vert_retrieve_3d(outdat,profdat,ndims,xlims,ylims,ilist)
!************************************************************************
!
! *vert_retrieve_3d* Select output for plotting from a 2-D array of vertical
!                    profiles
!
! Author - Patrick Luyten
!
! Description - plot values are selected by specifiers iform, kplot, dplot
!
! Calling program - HT_plot
!
!************************************************************************
!
USE depths
USE grid
USE iopars
USE plotpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
INTEGER, INTENT(IN) :: ilist
INTEGER, INTENT(IN), DIMENSION(2) :: ndims
INTEGER, INTENT(IN), DIMENSION(3) :: xlims, ylims
REAL, INTENT(IN), DIMENSION(ncout,nrout,nzout) :: profdat
REAL, INTENT(OUT), DIMENSION(ndims(1),ndims(2)) :: outdat

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*outdat*  REAL    Retrieved output values
!*profdat* REAL    Vertical data profiles
!*ndims*   INTEGER Shape of retrieved output array
!*xlims*   INTEGER Start/end/increment of selected output data in X-direction
!*ylims*   INTEGER Start/end/increment of selected output data in Y-direction
!*ilist*   INTEGER Plot list number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ii, iplotform, j, jj, k, kpos, k1, k2
REAL :: depplot, deptot, ratio
INTEGER, DIMENSION(1) :: klo, kup
REAL, DIMENSION(nzout) :: zdat


procname(pglev+1) = 'vert_retrieve_3d'
CALL log_timer_in()

!
!1. Initialise
!-------------
!

iplotform = postatts(ilist)%iplotform
kpos = postatts(ilist)%kpos
depplot = postatts(ilist)%depplot

!
!2. Evaluate
!-----------
!

i_200:DO i=xlims(1),xlims(2),xlims(3)
j_200:DO j=ylims(1),ylims(2),ylims(3)
   ii = (i-xlims(1))/xlims(3) + 1
   jj = (j-ylims(1))/ylims(3) + 1
   IF (maskatc_int(i,j)) THEN

!
!2.1 Initialise
!--------------
!

      IF (iplotform.GT.1) THEN
         IF (ztime) THEN
            deptot = depout(i,j) + zeta(i,j)
            zdat = MERGE(zeta(i,j)-gzcoord(i,j,:),-gzcoord_ini(i,j,:),ztime)
         ELSE
            deptot = depout(i,j)
            zdat = -gzcoord_ini(i,j,:)
         ENDIF
      ENDIF

!
!2.2 Evaluate according to value of iform
!----------------------------------------
!

      SELECT CASE (iplotform)

!
!2.2.1 At level kplot
!------------------
!

      CASE (1)
         outdat(ii,jj) = profdat(i,j,kpos)

!
!2.2.2 At depth dplot
!--------------------
!

      CASE (2)
         IF (deptot.LT.depplot) THEN
            outdat(ii,jj) = fill_value
         ELSEIF (zdat(1).LT.depplot) THEN
            outdat(ii,jj) = profdat(i,j,1)
         ELSEIF (zdat(nzout).GT.depplot) THEN
            outdat(ii,jj) = profdat(i,j,nzout)
         ELSE
            klo = MINLOC(zdat-depplot,MASK=zdat.GE.depplot)
            kup = MINLOC(depplot-zdat,MASK=zdat.LT.depplot)
            k1 = klo(1); k2 = kup(1)
            ratio = (depplot-zdat(k2))/(zdat(k1)-zdat(k2))
            outdat(ii,jj) = ratio*profdat(i,j,k1)+(1.0-ratio)*profdat(i,j,k2)
         ENDIF

!
!2.2.3 Depth-averaged
!--------------------
!

      CASE (3)
         outdat(ii,jj) = 0.5*(zdat(nzout)+zdat(nzout-1))*profdat(i,j,nzout) +&
                            &(deptot-0.5*(zdat(1)+zdat(2)))*profdat(i,j,1)
         k_223: DO k=2,nzout-1
            outdat(ii,jj) = outdat(ii,jj) + &
                          & 0.5*(zdat(k-1)-zdat(k+1))*profdat(i,j,k)
         ENDDO k_223
         outdat = outdat/deptot

!
!2.2.4 Depth-integrated
!----------------------
!

      CASE (4)
         outdat(ii,jj) = 0.5*(zdat(nzout)+zdat(nzout-1))*profdat(i,j,nzout) +&
                        &(deptot-0.5*(zdat(1)+zdat(2)))*profdat(i,j,1)
         k_224: DO k=2,nzout-1
            outdat(ii,jj) = outdat(ii,jj) + &
                          & 0.5*(zdat(k-1)-zdat(k+1))*profdat(i,j,k)
         ENDDO k_224

!
!2.2.5 Maximum value in the vertical
!-----------------------------------
!

      CASE (5)
         outdat(ii,jj) = MAXVAL(profdat(i,j,:))

      END SELECT

   ELSE
   
      outdat(ii,jj) = fill_value

   ENDIF

ENDDO j_200
ENDDO i_200


CALL log_timer_out()


RETURN

END SUBROUTINE vert_retrieve_3d

!========================================================================

SUBROUTINE VT_plot
!************************************************************************
!
! *VT_plot* Vertical transects
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - post_program
!
! External calls - coords_to_index, post_read, time_limits,
!                  write_mumap_cont, write_mumap_par, write_mumap_vec
!
! Module calls - close_file, convert_loc_to_char, error_alloc, lim_dims,
!                loop_index, open_file
!
!************************************************************************
!
USE depths
USE grid
USE iopars
USE physpars
USE plotpars
USE switches
USE syspars
USE error_routines, ONLY : error_alloc
USE grid_routines, ONLY: convert_loc_to_char
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims, loop_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: anim, flag, header
CHARACTER (LEN=12) :: cfig, cfil, cnt
CHARACTER (LEN=lendesc) :: cloc
CHARACTER (LEN=lentime) :: ctime
INTEGER :: horz_unit, i, iend, ifil, ilist, istart, itrtype, ivar, ivec1, &
         & ivec2, ivec3, iunit, j, jend, jstart, k, kmax, kmin, l, n, nofiles, &
         & nolevels, norec, notimes, novarskiph, novarskipz, novecskiph, &
         & novecskipz, ntout, numpoints, n1dim, n2dim, vert_unit
REAL :: aspect, contmax, contmin, deli, delj, depmax, depmin, dels, fac, &
      & gxend, gxstart, gyend, gystart, hlen_unit, trcos, trsin, &
      & valmax, valmin, xborder, yborder, zmax, zmin, z1, z2
REAL (KIND=kndrlong) :: rtime

CHARACTER (LEN=12), DIMENSION(4) :: cpos
CHARACTER (LEN=21), DIMENSION(2) :: geoloc
INTEGER, DIMENSION(3) :: slims, tlims, zlims
INTEGER, DIMENSION(nolists) :: nopoints
REAL, DIMENSION(MaxContLevs) :: contlevels
REAL, DIMENSION(nzout) :: xcomp, ycomp

CHARACTER (LEN=lentime), ALLOCATABLE, DIMENSION(:) :: outtime
CHARACTER (LEN=leniofile), ALLOCATABLE, DIMENSION(:) :: paramfiles
CHARACTER (LEN=leniofile), ALLOCATABLE, DIMENSION(:,:) :: contfiles, &
                                                        & linefiles, vecfiles
INTEGER, ALLOCATABLE, DIMENSION(:) :: iopar, itrans, jtrans
REAL, ALLOCATABLE, DIMENSION(:) :: deptrans, trangle, xtrans, ytrans, zettrans
REAL, ALLOCATABLE, DIMENSION(:,:) :: conttrans, dtrans, ztrans, ztransout
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: vectrans


procname(pglev+1) = 'VT_plot'
CALL log_timer_in()

!
!1.Time parameters and arrays
!----------------------------
!
!---output time limiters
CALL time_limits(tlims)

!---number of time steps
notimes = lim_dims(tlims)

!---allocate time series arrays
ALLOCATE (outtime(notimes),STAT=errstat) 
CALL error_alloc('outtime',1,(/notimes/),kndchar,lenstr=lentime)
ALLOCATE (contfiles(nolists,notimes),STAT=errstat)
CALL error_alloc('contfiles',2,(/nolists,notimes/),kndchar,lenstr=leniofile)
ALLOCATE (vecfiles(nolists,notimes),STAT=errstat)
CALL error_alloc('vecfiles',2,(/nolists,notimes/),kndchar,lenstr=leniofile)
ALLOCATE (linefiles(nolists,notimes),STAT=errstat)
CALL error_alloc('linefiles',2,(/nolists,notimes/),kndchar,lenstr=leniofile)
ALLOCATE (paramfiles(notimes),STAT=errstat)
CALL error_alloc('paramfiles',1,(/notimes/),kndchar,lenstr=leniofile)
ALLOCATE (iopar(notimes),STAT=errstat)
CALL error_alloc('iopar',1,(/notimes/),kndint)

!
!2. Initialise
!-------------
!

ilist_200: DO ilist=1,nolists

!
!2.1  Horizontal grid locations
!------------------------------
!
   
   IF (postatts(ilist)%horz_unit.EQ.1) THEN
      istart = NINT(postatts(ilist)%gxstart)
      iend = NINT(postatts(ilist)%gxend)
      jstart = NINT(postatts(ilist)%gystart)
      jend = NINT(postatts(ilist)%gyend)
      postatts(ilist)%gxstart = gxcoord(istart,jstart)
      postatts(ilist)%gxend = gxcoord(iend,jend)
      postatts(ilist)%gystart = gycoord(istart,jstart)
      postatts(ilist)%gyend = gycoord(iend,jend)
      postatts(ilist)%istart = istart
      postatts(ilist)%iend = iend
      postatts(ilist)%jstart = jstart
      postatts(ilist)%jend = jend
   ELSE
      CALL coords_to_index(postatts(ilist)%gxstart,postatts(ilist)%gystart,&
                         & istart,jstart)
      CALL coords_to_index(postatts(ilist)%gxend,postatts(ilist)%gyend,&
                         & iend,jend)
      postatts(ilist)%istart = istart
      postatts(ilist)%iend = iend
      postatts(ilist)%jstart = jstart
      postatts(ilist)%jend = jend
   ENDIF

!
!2.2 Number of locations along transect
!--------------------------------------
!

   nopoints(ilist) = MAX(ABS(iend-istart),ABS(jend-jstart)) + 1

ENDDO ilist_200

!
!3. Skip data reads
!------------------
!

IF (outform.NE.'N') THEN
   norec_310: DO norec=1,tlims(1)-1
      CALL post_read(norec,ctime,rtime,.TRUE.)
   ENDDO norec_310
ENDIF

!
!4. Time loop
!------------
!

ntout = 0
norec_400: DO norec=tlims(1),tlims(2)

!
!4.1 Read data
!-------------
!

   flag = loop_index(tlims,norec)
   IF (flag.OR.outform.NE.'N') THEN
      CALL post_read(norec,ctime,rtime,.TRUE.)
   ENDIF
   IF (.NOT.flag) CYCLE norec_400
   ntout = ntout + 1; outtime(ntout) = ctime
   WRITE (cnt,'(I12)') ntout; cnt = ADJUSTL(cnt)

!
!4.2 Process data
!----------------
!

   ilist_420: DO ilist=1,nolists

!
!4.2.1 Initialise parameters
!---------------------------
!

      anim = postatts(ilist)%anim
      numpoints = nopoints(ilist)
      ivar = postatts(ilist)%numvar
      ivec1 = postatts(ilist)%numvec1
      ivec2 = postatts(ilist)%numvec2
      ivec3 = postatts(ilist)%numvec3
      hlen_unit = postatts(ilist)%hlen_unit
      horz_unit = postatts(ilist)%horz_unit
      vert_unit = postatts(ilist)%vert_unit
      istart = postatts(ilist)%istart; iend = postatts(ilist)%iend
      jstart = postatts(ilist)%jstart; jend = postatts(ilist)%jend
      novarskiph = postatts(ilist)%novarskip1
      novarskipz = postatts(ilist)%novarskip2
      novecskiph = postatts(ilist)%novecskip1
      novecskipz = postatts(ilist)%novecskip2
      nolevels = postatts(ilist)%nolevels
      contmin = postatts(ilist)%contmin
      contmax = postatts(ilist)%contmax
      depmin = postatts(ilist)%depmin
      depmax = postatts(ilist)%depmax
      WRITE (cfig,'(I12)') noplots+ilist; cfig = ADJUSTL(cfig)

!
!4.2.2 Allocate arrays
!---------------------
!
!

      IF (norec.EQ.tlims(1).OR.nolists.GT.1) THEN
         ALLOCATE (dtrans(numpoints,nzout),STAT=errstat)
         CALL error_alloc('dtrans',2,(/numpoints,nzout/),kndrtype)
         ALLOCATE (itrans(numpoints),STAT=errstat)
         CALL error_alloc('itrans',1,(/numpoints/),kndint)
         ALLOCATE (jtrans(numpoints),STAT=errstat)
         CALL error_alloc('jtrans',2,(/numpoints/),kndint)
         ALLOCATE (xtrans(numpoints),STAT=errstat)
         CALL error_alloc('xtrans',1,(/numpoints/),kndrtype)
         ALLOCATE (ytrans(numpoints),STAT=errstat)
         CALL error_alloc('ytrans',2,(/numpoints/),kndrtype)
         ALLOCATE (ztrans(numpoints,nzout),STAT=errstat)
         CALL error_alloc('ztrans',2,(/numpoints,nzout/),kndrtype)
         ALLOCATE (ztransout(numpoints,nzout),STAT=errstat)
         CALL error_alloc('ztransout',2,(/numpoints,nzout/),kndrtype)
         ALLOCATE (trangle(numpoints),STAT=errstat)
         CALL error_alloc('trangle',2,(/numpoints/),kndrtype)
         ALLOCATE (deptrans(numpoints),STAT=errstat)
         CALL error_alloc('deptrans',1,(/numpoints/),kndrtype)
         ALLOCATE (zettrans(numpoints),STAT=errstat)
         CALL error_alloc('zettrans',1,(/numpoints/),kndrtype)
         IF (ivar.NE.0) THEN
            ALLOCATE (conttrans(numpoints,nzout),STAT=errstat)
            CALL error_alloc('conttrans',2,(/numpoints,nzout/),kndrtype)
         ENDIF
         IF (ivec1.NE.0) THEN
            ALLOCATE (vectrans(numpoints,nzout,3),STAT=errstat)
            CALL error_alloc('vectrans',3,(/numpoints,nzout,3/),kndrtype)
         ENDIF
      ENDIF

!
!4.2.3 Define transect
!---------------------
!

!     ---type of transect
      IF (jstart.EQ.jend) THEN
         itrtype = 1
      ELSEIF (istart.EQ.iend) THEN
         itrtype = 2
      ELSE
         itrtype = 3
      ENDIF

!     ---coordinate arrays along transect
      SELECT CASE (itrtype)
      CASE (1)
         itrans = (/(i,i=istart,iend)/)
         jtrans = jstart
         dtrans(:,1) = (/(gxcoord(i,jstart),i=istart,iend)/)
         trangle = 0.0
      CASE (2)
         itrans = istart
         jtrans = (/(j,j=jstart,jend)/)
         dtrans(:,1) = (/(gycoord(istart,j),j=jstart,jend)/)
         trangle = halfpi
      CASE (3)
         deli = (iend-istart)/REAL(numpoints-1)
         delj = (jend-jstart)/REAL(numpoints-1)
         itrans = (/(NINT(istart+(n-1)*deli),n=1,numpoints)/)
         jtrans = (/(NINT(jstart+(n-1)*delj),n=1,numpoints)/)
         xtrans = (/(gxcoord(itrans(n),jtrans(n)),n=1,numpoints)/)
         ytrans = (/(gycoord(itrans(n),jtrans(n)),n=1,numpoints)/)
         dtrans(1,1) = 0.0
         n_4231: DO n=2,numpoints
            IF (iopt_grid_sph.EQ.0) THEN
               dels = SQRT((xtrans(n)-xtrans(n-1))**2+&
                    & (ytrans(n)-ytrans(n-1))**2)
            ELSE
               fac = COS(0.5*degtorad*(ytrans(n-1)+ytrans(n)))
               dels = Rearth*degtorad*SQRT((fac*(xtrans(n)-xtrans(n-1)))**2+&
                    & (ytrans(n)-ytrans(n-1))**2)
            ENDIF
            dtrans(n,1) = dtrans(n-1,1) + dels
            trangle(n-1) = ATAN2(ytrans(n)-ytrans(n-1),xtrans(n)-xtrans(n-1))
         ENDDO n_4231
         trangle(numpoints) = trangle(numpoints-1)
      END SELECT
      k_4232: DO k=2,nzout
         dtrans(:,k) = dtrans(:,1)
      ENDDO k_4232
      IF (iopt_grid_sph.EQ.0.OR.itrtype.EQ.3) dtrans= hlen_unit*dtrans

!
!4.2.4 Vertical grid and depth arrays
!------------------------------------
!

      n_424: DO n=1,numpoints
         i = itrans(n); j = jtrans(n)
         deptrans(n) = depout(i,j)
         ztrans(n,:) = gzcoord_ini(i,j,:)
         zettrans(n) = 0.0
      ENDDO n_424
      IF (vert_unit.EQ.1) THEN
         kmin = depmin; kmax = depmax
         zmax = MERGE(MAXVAL(deptrans),MAXVAL(ABS(ztrans(:,kmin))),kmin.EQ.1)
         zmin = MERGE(0.0,MINVAL(ABS(ztrans(:,kmax))),kmax.EQ.nzout)
      ELSE
         kmin = 1; kmax = nzout
         zmax = MERGE(MAXVAL(deptrans),depmax,depmax.EQ.0.0)
         zmin = -depmin
      ENDIF

!
!4.2.5 Z-coordinates and surface elevations
!------------------------------------------
!

      IF (ztime) THEN
         n_425: DO n=1,numpoints
            i = itrans(n); j = jtrans(n)
            ztransout(n,:) = gzcoord(i,j,:)
            zettrans(n) = zeta(i,j)
         ENDDO n_425
      ELSE
         ztransout = ztrans
      ENDIF

!
!4.2.6 Vector data
!-----------------
!

      IF (ivec1.NE.0) THEN

!        ---evaluate vector components along transect
         n_4261: DO n=1,numpoints
            i = itrans(n); j = jtrans(n)
            trcos = COS(trangle(n)); trsin = SIN(trangle(n))
            xcomp = outvals4d(i,j,:,ivec1)
            ycomp = outvals4d(i,j,:,ivec2)
            IF (trsin.EQ.0.0) THEN
               vectrans(n,:,1) = xcomp
               vectrans(n,:,2) = ycomp
            ELSEIF (trcos.EQ.0.0) THEN
               vectrans(n,:,1) = ycomp
               vectrans(n,:,2) = -xcomp
            ELSE
               WHERE ((ABS(xcomp-fill_value).GT.fill_value_eps).AND.&
                    & (ABS(ycomp-fill_value).GT.fill_value_eps))
                  vectrans(n,:,1) = trcos*xcomp + trsin*ycomp
                  vectrans(n,:,2) = trcos*ycomp - trsin*xcomp
               ELSEWHERE
                  vectrans(n,:,1) = fill_value; vectrans(n,:,2) = fill_value
               END WHERE
            ENDIF
            vectrans(n,:,3) = outvals4d(i,j,:,ivec3)
         ENDDO n_4261

!        ---eliminate data outside depth range
         IF (vert_unit.EQ.2) THEN
            l_4262: DO l=1,2
               WHERE (ABS(ztrans).LT.zmin.OR.ABS(ztrans).GT.zmax)
                  vectrans(:,:,l) = fill_value
               END WHERE
            ENDDO l_4262
         ENDIF

!        ---create vector file
         vecfiles(ilist,ntout) = TRIM(plotfile)//'.'//TRIM(cnt)//'VTgdat'&
                             & //TRIM(cfig)
         CALL open_file(iunit,vecfiles(ilist,ntout),'OUT','A')

!        ---write data
         slims = (/1,numpoints,novecskiph/)
         zlims = (/kmin,kmax,novecskipz/)
         n1dim = lim_dims(slims); n2dim = lim_dims(zlims)
         CALL write_mumap_vec(iunit,'g',n1dim,n2dim,&
           & dtrans(slims(1):slims(2):slims(3),zlims(1):zlims(2):zlims(3)),&
           & ztransout(slims(1):slims(2):slims(3),zlims(1):zlims(2):zlims(3)),&
           & vectrans(slims(1):slims(2):slims(3),zlims(1):zlims(2):zlims(3),&
           & 1:3:2),postatts(ilist)%hrefvec, postatts(ilist)%vrefvec,0)
         CALL close_file(iunit,'A')

      ELSE
         vecfiles(ilist,ntout) = ''
      ENDIF

!
!4.2.7 Contour data
!------------------
!

      IF (ivar.NE.0) THEN

!        ---scalar data
         IF (ivar.GT.0) THEN
            n_4271: DO n=1,numpoints
               i = itrans(n); j = jtrans(n)
               conttrans(n,:) = outvals4d(i,j,:,ivar)
            ENDDO n_4271
         ELSEIF (ivar.EQ.-1.AND.vcoord.EQ.2) THEN
            n_4272: DO n=1,numpoints
               IF (iopt_grid_vtype.EQ.2) THEN
                  conttrans(n,:) = gsigcoordatc
               ELSE
                  i = itrans(n); j = jtrans(n)
                  conttrans(n,:) = gsout(i,j,:)
               ENDIF
            ENDDO n_4272
         ENDIF
         IF (vert_unit.EQ.2) THEN
            WHERE (ABS(ztrans).LT.zmin.OR.ABS(ztrans).GT.zmax)
               conttrans = fill_value
            END WHERE
         ENDIF

!        ---contour levels
         slims = (/1,numpoints,novarskiph/)
         zlims = (/kmin,kmax,novarskipz/)
         IF (postatts(ilist)%contreg) THEN
            IF (contmin.EQ.0.0.AND.contmax.EQ.0.0) THEN
               valmin = conttrans(1,1); valmax = conttrans(1,1)
               n_4273: DO n=slims(1),slims(2),slims(3)
               k_4273: DO k=zlims(1),zlims(2),zlims(3)
                  IF (ABS(conttrans(n,k)-fill_value).GT.fill_value_eps) THEN
                     valmin = MIN(valmin,conttrans(n,k))
                     valmax = MAX(valmax,conttrans(n,k))
                  ENDIF
               ENDDO k_4273
               ENDDO n_4273
            ELSE
               valmin = contmin; valmax = contmax
            ENDIF
            contlevels(1:nolevels) = (/(valmin+(n-1)*(valmax-valmin)/&
                                     & (nolevels-1),n=1,nolevels)/)
         ELSE
            contlevels(1:nolevels) = postatts(ilist)%contlevels(1:nolevels)
         ENDIF

!        ---create contour file
         contfiles(ilist,ntout) = TRIM(plotfile)//'.'//TRIM(cnt)//'VTidat'&
                              & //TRIM(cfig)
         CALL open_file(iunit,contfiles(ilist,ntout),'OUT','A')

!        ---write data
         n1dim = lim_dims(slims); n2dim = lim_dims(zlims)
         CALL write_mumap_cont(iunit,n1dim,n2dim,&
           & dtrans(slims(1):slims(2):slims(3),zlims(1):zlims(2):zlims(3)),&
           & ztransout(slims(1):slims(2):slims(3),zlims(1):zlims(2):zlims(3)),&
           & conttrans(slims(1):slims(2):slims(3),zlims(1):zlims(2):zlims(3)),&
           & nolevels,contlevels(1:nolevels),postatts(ilist)%icontstyle,0,2)
         CALL close_file(iunit,'A')

      ELSE
         contfiles(ilist,ntout) = ''
      ENDIF

!
!4.2.8 Draw bottom/surface lines
!-------------------------------
!
!     ---create file
      linefiles(ilist,ntout) = TRIM(plotfile)//'.'//TRIM(cnt)//'VTcdat'//&
                             & TRIM(cfig)
      CALL open_file(iunit,linefiles(ilist,ntout),'OUT','A')

!     ---write data
      WRITE (iunit,9001) dtrans(1,1), -deptrans(1), 0
      n_4281: DO n=2,numpoints
         WRITE (iunit,9001) dtrans(n,1), -deptrans(n), 1
      ENDDO n_4281
      IF (ztime) THEN
         WRITE (iunit,9001) dtrans(1,1), zettrans(1), 0
         n_4282: DO n=2,numpoints
            WRITE (iunit,9001) dtrans(n,1), zettrans(n), 1
         ENDDO n_4282
      ENDIF
      CALL close_file(iunit,'A')

!
!4.2.9 Write parameter file(s)
!----------------------------
!

      IF (norec.EQ.tlims(2)) THEN

!
!4.2.9.1 Initalise parameters
!----------------------------
!

         istart = postatts(ilist)%istart; iend = postatts(ilist)%iend
         jstart = postatts(ilist)%jstart; jend = postatts(ilist)%jend
         gxstart = postatts(ilist)%gxstart; gxend = postatts(ilist)%gxend
         gystart = postatts(ilist)%gystart; gyend = postatts(ilist)%gyend

!
!4.2.9.2 Create/open file(s)
!-------------------------
!

         nofiles = MERGE(1,notimes,anim)
         IF (anim) THEN
            paramfiles(1) = TRIM(plotfile)//'.VTpar'//TRIM(cfig)
            CALL open_file(iopar(1),paramfiles(1),'OUT','A')
         ELSE
            ifil_4292: DO ifil=1,nofiles
               WRITE (cfil,'(I12)') ifil; cfil = ADJUSTL(cfil)
               paramfiles(ifil) = TRIM(plotfile)//'.'//TRIM(cfil)//'VTpar'//&
                                & TRIM(cfig)
               CALL open_file(iopar(ifil),paramfiles(ifil),'OUT','A')
            ENDDO ifil_4292
         ENDIF

!
!4.2.9.3 Write parameter file(s)
!---------------------------
!
!        ---plot borders
         IF (depmin.EQ.0.0) THEN
            z1 = MAXVAL(zettrans)
         ELSE
            z1 = zmin
         ENDIF
         z2 = zmax
         xborder = border*(dtrans(numpoints,1)-dtrans(1,1))
         yborder = border*(z2+z1)
         plotcorners(1) = dtrans(1,1) - xborder
         plotcorners(2) = -z2 - yborder
         plotcorners(3) = dtrans(numpoints,1) + xborder
         plotcorners(4) = yborder + z1
         aspect = (dtrans(numpoints,1)-dtrans(1,1))/(z2+z1)

!        ---main titles
         plottitles(1) = TRIM(outfile)

!        ---first subtitle
         IF (ivar.NE.0) THEN
            plottitles(2) = TRIM(outvars(ivar)%long_name)//' ('//&
                          & TRIM(outvars(ivar)%units)//')'
         ELSE
            plottitles(2) = TRIM(outvars(ivec1)%vector_name)//' ('//&
                          & TRIM(outvars(ivec1)%units)//')'
         ENDIF

!        ---second subtitle 
         SELECT CASE (itrtype)
         CASE (1)
            IF (iopt_grid_sph.EQ.0.OR.horz_unit.EQ.1) THEN
               WRITE (cpos(1),'(I12)') jstart
               cpos(1) = ADJUSTL(cpos(1))
               cloc = 'Transect at j = '//TRIM(cpos(1))
            ELSEIF (iopt_grid_sph.EQ.1) THEN
               geoloc(1) = convert_loc_to_char(gxstart,gystart)
               cloc = 'Transect at '//geoloc(1)(1:9)
            ENDIF
         CASE (2)
            IF (iopt_grid_sph.EQ.0.OR.horz_unit.EQ.1) THEN
               WRITE (cpos(1),'(I12)') istart
               cpos(1) = ADJUSTL(cpos(1))
               cloc = 'Transect at i = '//TRIM(cpos(1))
            ELSEIF (iopt_grid_sph.EQ.1) THEN
               geoloc(1) = convert_loc_to_char(gxstart,gystart)
               cloc = 'Transect at '//geoloc(1)(12:21)
            ENDIF
         CASE (3)
            IF (iopt_grid_sph.EQ.0.OR.horz_unit.EQ.1) THEN
               WRITE (cpos(1),'(I12)') istart
               cpos(1) = ADJUSTL(cpos(1))
               WRITE (cpos(2),'(I12)') jstart
               cpos(2) = ADJUSTL(cpos(2))
               WRITE (cpos(3),'(I12)') iend; cpos(3) = ADJUSTL(cpos(3))
               WRITE (cpos(4),'(I12)') jend; cpos(4) = ADJUSTL(cpos(4))
               cloc = 'From ('//TRIM(cpos(1))//','//TRIM(cpos(2))//' to ('//&
              & TRIM(cpos(3))//','//TRIM(cpos(4))//')'
            ELSEIF (iopt_grid_sph.EQ.1) THEN
               geoloc(1) = convert_loc_to_char(gxstart,gystart)
               geoloc(2) = convert_loc_to_char(gxend,gyend)
               cloc = 'From '//TRIM(geoloc(1))//' to '//TRIM(geoloc(2))
            ENDIF
         END SELECT

!        ---write parameter file(s)
         n_42931: DO n=1,notimes
            ifil = MERGE(1,n,anim)
            header = .NOT.anim.OR.n.EQ.1
            plottitles(3) = TRIM(cloc)//'; '//outtime(n)(1:19)
            IF (ivec1.NE.0) THEN
               CALL write_mumap_par(iopar(ifil),header,0,aspect,3,&
                                 & (/'i','g','c'/),&
                                 & (/contfiles(ilist,n),vecfiles(ilist,n),&
                                 & linefiles(ilist,n)/))
            ELSE
               CALL write_mumap_par(iopar(ifil),header,0,aspect,2,&
                                 & (/'i','c'/),(/contfiles(ilist,n),&
                                 & linefiles(ilist,n)/))
            ENDIF
            IF (anim.AND.n.LT.notimes) WRITE (iopar(ifil),'(A)') 'z'
         ENDDO n_42931

!        ---close parameter file(s)
         ifil_42932: DO ifil=1,nofiles
            CALL close_file(iopar(ifil),'A')
         ENDDO ifil_42932

!
!4.2.9.4 Write to '.vis'-file
!----------------------------
!

         ifil_4294: DO ifil=1,nofiles
            IF (ivar.NE.0) THEN
               WRITE (io_parvis,9002) 'mumap2', TRIM(paramfiles(ifil)), 'y',&
                                              & TRIM(outvars(ivar)%f90_name)
            ELSEIF (ivec1.NE.0) THEN
               WRITE (io_parvis,9002) 'mumap2', TRIM(paramfiles(ifil)), 'y',&
                                             & TRIM(outvars(ivec1)%vector_name)
            ENDIF
         ENDDO ifil_4294

      ENDIF

!
!4.2.10 Deallocate arrays
!------------------------
!

      IF (norec.EQ.tlims(2).OR.nolists.GT.1) THEN
         DEALLOCATE (dtrans,xtrans,ytrans,ztrans,ztransout,trangle,deptrans,&
                   & zettrans)
         DEALLOCATE (itrans,jtrans)
         IF (ivar.NE.0) DEALLOCATE (conttrans)
         IF (ivec1.NE.0) DEALLOCATE (vectrans)
      ENDIF


   ENDDO ilist_420

ENDDO norec_400

!
!5. Update plot number
!---------------------
!

noplots = noplots + nolists

!
!6. Deallocate time series arrays
!--------------------------------
!

DEALLOCATE (outtime,contfiles,vecfiles,linefiles,paramfiles,iopar)

CALL log_timer_out()


RETURN

9001 FORMAT (2(G15.7,1X),I1)
9002 FORMAT (4(A,1X))

END SUBROUTINE VT_plot

!========================================================================

SUBROUTINE write_mumap_coast(iunit,ilist,linefile)
!************************************************************************
!
! *write_mumap_coast* Write mumap2 'class c' coast-line file
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - HT_plot
!
! External calls - coords_to_index
!
! Module calls - close_file, error_alloc, open_file, Carr_at_UV
!
!************************************************************************
!
USE grid
USE iopars
USE modids
USE plotpars
USE switches
USE syspars
USE array_interp, ONLY: Carr_at_UV
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
CHARACTER (LEN=leniofile), INTENT(INOUT) :: linefile
INTEGER, INTENT(INOUT) :: ilist, iunit

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iunit*     INTEGER File unit
!*ilist*     INTEGER Plot list number
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=120) :: cline
INTEGER :: i, iend, im, iocoast, istart, j, jend, jstart, nolines
REAL :: gxend, gxstart, gyend, gystart, x, y
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: xcoord, ycoord


procname(pglev+1) = 'write_mumap_coast'
CALL log_timer_in()

!
!1. Locate plot borders on output grid
!--------------------------------  ---
!

istart = postatts(ilist)%istart; iend = postatts(ilist)%iend
jstart = postatts(ilist)%jstart; jend = postatts(ilist)%jend
gxstart = postatts(ilist)%gxstart; gxend = postatts(ilist)%gxend
gystart = postatts(ilist)%gystart; gyend = postatts(ilist)%gyend

!
!2. Default case
!---------------
!

IF (coast_def) THEN

!
!2.1 Coordinate arrays at cell corners
!------------------------------------
!

!  ---allocate
   ALLOCATE (maskatc(0:ncout+1,0:nrout+1),STAT=errstat)
   CALL error_alloc('maskatc',2,(/ncout+2,nrout+2/),kndlog)
   ALLOCATE (xcoord(ncout+1,nrout+1),STAT=errstat)
   CALL error_alloc('xcoord',2,(/ncout+1,nrout+1/),kndrtype)
   ALLOCATE (ycoord(ncout+1,nrout+1),STAT=errstat)
   CALL error_alloc('ycoord',2,(/ncout+1,nrout+1/),kndrtype)

!  ---mask array
   maskatc = .FALSE.
   maskatc(1:ncout,1:nrout) = maskatc_int

!  ---interpolate at interior points
   CALL Carr_at_UV(gxcoord,xcoord(2:ncout,2:nrout),0,0,(/1,1,nzout/),&
                & (/ncout,nrout,nzout/),1,iarr_gxcoord,.TRUE.,hregular=.TRUE.)
   CALL Carr_at_UV(gycoord,ycoord(2:ncout,2:nrout),0,0,(/1,1,nzout/),&
                & (/ncout,nrout,nzout/),1,iarr_gycoord,.TRUE.,hregular=.TRUE.)

!  ---west boundary
   j_211: DO j=2,nrout
      xcoord(1,j) = 0.25*(3.0*gxcoord(1,j)+3.0*gxcoord(1,j-1)-&
                        & gxcoord(2,j)-gxcoord(2,j-1))
      ycoord(1,j) = 0.25*(3.0*gycoord(1,j)+3.0*gycoord(1,j-1)-&
                        & gycoord(2,j)-gycoord(2,j-1))
   ENDDO j_211

!  ---east boundary
   j_212: DO j=2,nrout
      xcoord(ncout+1,j) = 0.25*(3.0*gxcoord(ncout,j)+3.0*gxcoord(ncout,j-1)-&
                              & gxcoord(ncout-1,j)-gxcoord(ncout-1,j-1))
      ycoord(ncout+1,j) = 0.25*(3.0*gycoord(ncout,j)+3.0*gycoord(ncout,j-1)-&
                              & gycoord(ncout-1,j)-gycoord(ncout-1,j-1))
   ENDDO j_212

!  ---south boundary   
   i_213: DO i=2,ncout
      xcoord(i,1) = 0.25*(3.0*gxcoord(i,1)+3.0*gxcoord(i-1,1)-&
                        & gxcoord(i,2)-gxcoord(i-1,2))
      ycoord(i,1) = 0.25*(3.0*gycoord(i,1)+3.0*gycoord(i-1,1)-&
                        & gycoord(i,2)-gycoord(i-1,2))
   ENDDO i_213

!  ---north boundary
   i_214: DO i=2,ncout
      xcoord(i,nrout+1) = 0.25*(3.0*gxcoord(i,nrout)+3.0*gxcoord(i-1,nrout)-&
                              & gxcoord(i,nrout-1)-gxcoord(i-1,nrout-1))
      ycoord(i,nrout+1) = 0.25*(3.0*gycoord(i,nrout)+3.0*gycoord(i-1,nrout)-&
                              & gycoord(i,nrout-1)-gycoord(i-1,nrout-1))
   ENDDO i_214

!  ---southwest corner
   xcoord(1,1) = 0.25*(9.0*gxcoord(1,1)-3.0*gxcoord(1,2)-&
                     & 3.0*gxcoord(2,1)+gxcoord(2,2))
   ycoord(1,1) = 0.25*(9.0*gycoord(1,1)-3.0*gycoord(1,2)-&
                     & 3.0*gycoord(2,1)+gycoord(2,2))

!  ---southeast corner
   xcoord(ncout+1,1) = 0.25*(9.0*gxcoord(ncout,1)-3.0*gxcoord(ncout,2)-&
                           & 3.0*gxcoord(ncout-1,1)+gxcoord(ncout-1,2))
   ycoord(ncout+1,1) = 0.25*(9.0*gycoord(ncout,1)-3.0*gycoord(ncout,2)-&
                           & 3.0*gycoord(ncout-1,1)+gycoord(ncout-1,2))

!  ---Northwest corner
   xcoord(1,nrout+1) = 0.25*(9.0*gxcoord(1,nrout)-3.0*gxcoord(2,nrout)-&
                           & 3.0*gxcoord(1,nrout-1)+gxcoord(2,nrout-1))
   ycoord(1,nrout+1) = 0.25*(9.0*gycoord(1,nrout)-3.0*gycoord(2,nrout)-&
                           & 3.0*gycoord(1,nrout-1)+gycoord(2,nrout-1))

!  ---Northeast corner
   xcoord(ncout+1,nrout+1) = 0.25*(9.0*gxcoord(ncout,nrout)-&
                   & 3.0*gxcoord(ncout-1,nrout)-3.0*gxcoord(ncout,nrout-1)+&
                   & gxcoord(ncout-1,nrout-1))
   ycoord(ncout+1,nrout+1) = 0.25*(9.0*gycoord(ncout,nrout)-&
                   & 3.0*gycoord(ncout-1,nrout)-3.0*gycoord(ncout,nrout-1)+&
                   & gycoord(ncout-1,nrout-1))

!  ---apply length scale unit
   IF (iopt_grid_sph.EQ.0) THEN
      xcoord = postatts(ilist)%hlen_unit*xcoord
      ycoord = postatts(ilist)%hlen_unit*ycoord
   ENDIF

!
!2.2 Cartesian grid
!------------------
!

   nolines = 0

   IF (iopt_grid_sph.EQ.0) THEN

!    ---draw U-faces
     i_221: DO i=istart,iend+1
     j_221: DO j=jstart,jend
        IF ((i.EQ.1.AND.maskatc(i,j)).OR.&
          & (i.EQ.ncout+1.AND.maskatc(i-1,j))) THEN
           WRITE (iunit,9001) xcoord(i,j), ycoord(i,j), 0
           WRITE (iunit,9001) xcoord(i,j), ycoord(i,j+1), 1
           nolines = nolines + 1
        ELSEIF ((i.GT.1.AND.i.LT.ncout+1).AND.&
              & ((.NOT.maskatc(i-1,j).AND.maskatc(i,j)).OR.&
              & (maskatc(i-1,j).AND.(.NOT.maskatc(i,j))))) THEN
           WRITE (iunit,9001) xcoord(i,j), ycoord(i,j), 0
           WRITE (iunit,9001) xcoord(i,j), ycoord(i,j+1), 1
           nolines = nolines + 1
        ENDIF
     ENDDO j_221
     ENDDO i_221

!    ---draw V-faces
     i_222: DO i=istart,iend
     j_222: DO j=jstart,jend+1
        IF ((j.EQ.1.AND.maskatc(i,j)).OR.&
          & (j.EQ.nrout+1.AND.maskatc(i,j-1))) THEN
           WRITE (iunit,9001) xcoord(i,j), ycoord(i,j), 0
           WRITE (iunit,9001) xcoord(i+1,j), ycoord(i,j), 1
           nolines = nolines + 1
        ELSEIF ((j.GT.1.AND.j.LT.nrout+1).AND.&
              & ((.NOT.maskatc(i,j-1).AND.maskatc(i,j)).OR.&
              & (maskatc(i,j-1).AND.(.NOT.maskatc(i,j))))) THEN
           WRITE (iunit,9001) xcoord(i,j), ycoord(i,j), 0
           WRITE (iunit,9001) xcoord(i+1,j), ycoord(i,j), 1
           nolines = nolines + 1
        ENDIF
     ENDDO j_222
     ENDDO i_222

!
!2.3 Spherical grid
!------------------
!

  ELSEIF (iopt_grid_sph.EQ.1) THEN

!    ---draw U-faces
     i_231: DO i=istart,iend+1
     j_231: DO j=jstart,jend
        IF ((i.EQ.1.AND.maskatc(i,j)).OR.&
          & (i.EQ.ncout+1.AND.maskatc(i-1,j))) THEN
           WRITE (iunit,9001) ycoord(i,j), xcoord(i,j), 0
           WRITE (iunit,9001) ycoord(i,j+1), xcoord(i,j), 1
           nolines = nolines + 1
        ELSEIF ((i.GT.1.AND.i.LT.ncout+1).AND.&
              & ((.NOT.maskatc(i-1,j).AND.maskatc(i,j)).OR.&
              & (maskatc(i-1,j).AND.(.NOT.maskatc(i,j))))) THEN
           WRITE (iunit,9001) ycoord(i,j), xcoord(i,j), 0
           WRITE (iunit,9001) ycoord(i,j+1), xcoord(i,j), 1
           nolines = nolines + 1
        ENDIF
     ENDDO j_231
     ENDDO i_231

!    ---draw V-faces
     i_232: DO i=istart,iend
     j_232: DO j=jstart,jend+1
        IF ((j.EQ.1.AND.maskatc(i,j)).OR.&
          & (j.EQ.nrout+1.AND.maskatc(i,j-1))) THEN
           WRITE (iunit,9001) ycoord(i,j), xcoord(i,j), 0
           WRITE (iunit,9001) ycoord(i,j), xcoord(i+1,j), 1
           nolines = nolines + 1
        ELSEIF ((j.GT.1.AND.j.LT.nrout+1).AND.&
              & ((.NOT.maskatc(i,j-1).AND.maskatc(i,j)).OR.&
              & (maskatc(i,j-1).AND.(.NOT.maskatc(i,j))))) THEN
           WRITE (iunit,9001) ycoord(i,j), xcoord(i,j), 0
           WRITE (iunit,9001) ycoord(i,j), xcoord(i+1,j), 1
           nolines = nolines + 1
        ENDIF
     ENDDO j_232
     ENDDO i_232

  ENDIF

!
!2.4 Deallocate
!--------------
!

  DEALLOCATE (maskatc,xcoord,ycoord)

!
!2.5 Close file
!--------------
!

  IF (nolines.GT.0) THEN
     CALL close_file(iunit,'A')
  ELSE
     linefile = ''
     CALL close_file(iunit,'A',fildel=.TRUE.)
  ENDIF

!
!3. Non-default case
!-------------------
!

ELSE

   CALL open_file(iocoast,TRIM(coast_name),'IN','A')
99 READ (iocoast,'(A)',END=98) cline
   IF (cline(1:1).NE.'#') THEN
      READ (cline,*) x, y, im
      IF (x.LT.gxstart.OR.x.GT.gxend.OR.y.LT.gystart.OR.y.GT.gyend) im = 0
      WRITE (iunit,9001) x, y, im
   ELSE
      WRITE (iunit,'(A)') TRIM(cline)
   ENDIF
   GOTO 99
98 CALL close_file(iocoast,'A')

ENDIF

CALL log_timer_out()


RETURN

9001 FORMAT (2(G15.7,1X),I1)

END SUBROUTINE write_mumap_coast

!========================================================================

SUBROUTINE write_mumap_cont(iunit,nxdim,nydim,xcoord,ycoord,vcon,ncurv,cval,&
                          & istyle,icoord,ireg)
!************************************************************************
!
! *write_mumap_cont* Write mumap2 class 'i' file for contouring
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - HT_plot, TH_plot, TZ_plot, VT_plot
!
!************************************************************************
!
USE iopars
USE plotpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
INTEGER, INTENT(IN) :: icoord, ireg, istyle, iunit, ncurv, nxdim, nydim
REAL, INTENT(IN), DIMENSION(nxdim,nydim) :: vcon, xcoord, ycoord
REAL, INTENT(IN), DIMENSION(ncurv) :: cval

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iunit*     INTEGER File unit
!*nxdim*     INTEGER X-dimension of data file
!*nydim*     INTEGER Y-dimension of data file
!*xcoord*    REAL    X-coordinates of plot data
!*ycoord*    REAL    Y-coordinates of plot data
!*vcon*      REAL    Contour data
!*ncurv*     INTEGER Number of levels for contouring
!*cval*      CHAR    Contour values
!*istyle*    INTEGER Contour style
!                = 1 => isolines annotated and drawn, no filling
!                = 2 => isolines not annotated, not drawn, with filling
!                = 3 => isolines not annotated, drawn, with filling
!*icoord*    INTEGER Coordinate type
!                = 0 => Cartesian
!                = 1 => spherical
!*ireg*      INTEGER 'ireg' parameter for class "i" file
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, icolfi, ihilo, ilab, inter, j, lidraw
REAL :: conv, dX, dY, vref, X1o, X2o


procname(pglev+1) = 'write_mumap_cont'
CALL log_timer_in()

!
!1.Evaluate parameters
!---------------------
!
!---corner coordinates and grid spacings
IF (icoord.EQ.0) THEN
   X1o = xcoord(1,1); X2o = ycoord(1,1)
ELSEIF (icoord.GT.0) THEN
   X1o = ycoord(1,1); X2o = xcoord(1,1)
ENDIF
dX = xcoord(2,1) - xcoord(1,1); dY = ycoord(1,2) - ycoord(1,1)
!---style
SELECT CASE(istyle)
CASE(1); ilab = 1; lidraw = 1; icolfi = 0
CASE(2); ilab = 0; lidraw = 0; icolfi = 1
CASE(3); ilab = 0; lidraw = 1; icolfi = 1
END SELECT
ihilo = 0; inter = 0
!---reference value
vref = 0.0; conv = 1.0

!
!2. Write file
!-------------
!
!---contour parameters
WRITE (iunit,9001) icoord, nxdim, nydim, X1o, X2o, dX, dY, ireg
WRITE (iunit,'(A)') '#'
WRITE (iunit,9002) 1, fill_value
WRITE (iunit,'(6I3)') ncurv, ihilo, ilab, lidraw, icolfi, inter
WRITE (iunit,9003) vref, conv, ' '
!---contour levels
WRITE (iunit,9004) cval
!---contour data
SELECT CASE (ireg)
CASE(0)
   i_210: DO j=1,nydim
      WRITE (iunit,9004) vcon(:,j)
   ENDDO i_210
CASE(1)
   i_221: DO j=1,nydim
      WRITE (iunit,9004) vcon(:,j)
   ENDDO i_221
   IF (icoord.EQ.0) THEN
      i_222: DO j=1,nydim
         WRITE (iunit,9004) xcoord(:,j)
      ENDDO i_222
      i_223: DO j=1,nydim
         WRITE (iunit,9004) ycoord(:,j)
      ENDDO i_223
   ELSE
      i_224: DO j=1,nydim
         WRITE (iunit,9004) ycoord(:,j)
      ENDDO i_224
      i_225: DO j=1,nydim
         WRITE (iunit,9004) xcoord(:,j)
      ENDDO i_225
   ENDIF
CASE(2)
   IF (icoord.EQ.0) THEN
      j_231: DO j=1,nydim
      i_231: DO i=1,nxdim
         WRITE (iunit,9004) xcoord(i,j), ycoord(i,j), vcon(i,j)
      ENDDO i_231
      ENDDO j_231
   ELSE
      j_232: DO j=1,nydim
      i_232: DO i=1,nxdim
         WRITE (iunit,9004) ycoord(i,j), xcoord(i,j), vcon(i,j)
      ENDDO i_232
      ENDDO j_232
   ENDIF
END SELECT

CALL log_timer_out()


RETURN

9001 FORMAT (I1,2(1X,I4),4(1X,G15.7),1X,I1)
9002 FORMAT (I1,1X,G15.7)
9003 FORMAT (2(G15.7,1X),A)
9004 FORMAT (50(G15.7,1X))

END SUBROUTINE write_mumap_cont

!========================================================================

SUBROUTINE write_mumap_par(iunit,header,icoord,aspect,nofiles,fileclass,&
                         & filenames)
!************************************************************************
!
! *write_mumap_par* Write mumap2 parameter file
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - HT_plot, TH_plot, TZ_plot, VT_plot
!
!************************************************************************
!
USE iopars
USE plotpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
LOGICAL, INTENT(IN) :: header
INTEGER, INTENT(IN) :: icoord, iunit, nofiles
CHARACTER (LEN=1), INTENT(IN), DIMENSION(nofiles) :: fileclass
CHARACTER (LEN=*), INTENT(IN), DIMENSION(nofiles) :: filenames
REAL, INTENT(IN) :: aspect

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iunit*     INTEGER File unit
!*header*    LOGICAL Enables/disables header lines in parameter file
!*icoord*    INTEGER Coordinate type
!                = 0 => Cartesian
!                = 1 => spherical
!*aspect*    REAL    Aspect ratio
!*nofiles*   INTEGER Number of data files
!*fileclass* CHAR    Data file classes
!*filenames* CHAR    Data file names
!
!------------------------------------------------------------------------------
!
!
!*Local variables
!
INTEGER :: icolbg, icside, iorien, iscale, isfram, isquad, itform, ivquad, n
REAL :: ticsx, ticsy


procname(pglev+1) = 'write_mumap_par'
CALL log_timer_in()

!---set parameters by default
iorien = 0; isfram = 2; icolbg = 9; iscale = 0; isquad = 0; ivquad = 5
ticsx = 0.0; ticsy = 0.0; itform = 1; icside = 5

!---write parameter file
IF (header) THEN
   WRITE (iunit,9001) icoord, aspect, iorien
   WRITE (iunit,9002) isfram, icolbg, iscale, isquad, ivquad, ticsx, ticsy,&
                    & itform, icside
   WRITE (iunit,9003) plotcorners
ENDIF
n_110: DO n=1,3
   IF (TRIM(plottitles(n)).EQ.'') THEN
      WRITE (iunit,*)
   ELSE
      WRITE (iunit,'(A)') TRIM(plottitles(n))
   ENDIF
ENDDO n_110
n_120: DO n=1,nofiles
   IF (TRIM(filenames(n)).NE.'') THEN
      WRITE (iunit,'(3A)') fileclass(n), ' ', TRIM(filenames(n))
   ENDIF
ENDDO n_120

CALL log_timer_out()


RETURN

9001 FORMAT (I1,1X,G15.7,1X,I1)
9002 FORMAT (5I2,2(1X,G15.7),2I2)
9003 FORMAT (4(G15.7,1X))

END SUBROUTINE write_mumap_par

!========================================================================

SUBROUTINE write_mumap_vec(iunit,vecclass,nxdim,nydim,xcoord,ycoord,vecdat,&
                         & xvref,yvref,icoord)
!************************************************************************
!
! *write_mumap_vec* Write mumap2 class 'g' or class 'h' file for vector plots
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - HT_plot, VT_plot
!
!************************************************************************
!
USE iopars
USE plotpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
INTEGER, INTENT(IN) :: icoord, iunit, nxdim, nydim
CHARACTER (LEN=1), INTENT(IN) :: vecclass
REAL, INTENT(IN), DIMENSION(nxdim,nydim) :: xcoord, ycoord
REAL, INTENT(IN), DIMENSION(nxdim,nydim,2) :: vecdat
REAL, INTENT(IN) :: xvref, yvref

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iunit*     INTEGER File unit
!*vecclass*  CHAR    Vector class of data file
!*nxdim*     INTEGER X-dimension of data file
!*nydim*     INTEGER Y-dimension of data file
!*xcoord*    REAL    X-coordinates of plot data
!*ycoord*    REAL    Y-coordinates of plot data
!*vecdat*    REAL    Vector data
!*xvref*     REAL    Reference value for X-component of vector
!*yvref*     REAL    Reference value for Y-component of vector
!*icoord*    INTEGER Coordinate type
!                = 0 => Cartesian
!                = 1 => spherical
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ifthr, j, n
CHARACTER (LEN=15) :: legvec
REAL :: dX, dY, X1o, X2o


procname(pglev+1) = 'write_mumap_vec'
CALL log_timer_in()

!
!1. Set parameters
!-----------------
!

IF (vecclass.EQ.'h') THEN
   IF (icoord.EQ.0) THEN
      X1o = xcoord(1,1); X2o = ycoord(1,1)
   ELSEIF (icoord.GT.0) THEN
      X1o = ycoord(1,1); X2o = xcoord(1,1)
   ENDIF
   dX = xcoord(2,1) - xcoord(1,1); dY = ycoord(1,2) - ycoord(1,1)
ENDIF
ifthr = 0

!
!2. Legends
!----------
!

WRITE (legvec,'(G15.7)') xvref; legvec = ADJUSTL(legvec)
WRITE (iunit,'(A)') TRIM(legvec)
WRITE (legvec,'(G15.7)') yvref; legvec = ADJUSTL(legvec)
WRITE (iunit,'(A)') TRIM(legvec)
WRITE (iunit,'(A)') '#'

!
!3. Plot parameters
!------------------
!

IF (vecclass.EQ.'g') THEN
   WRITE (iunit,9001) xvref, yvref, ifthr
ELSEIF (vecclass.EQ.'h') THEN
   WRITE (iunit,9002) icoord, nxdim, nydim, X1o, X2o, dX, dY,&
                    & xvref, yvref, fill_value, ifthr
ENDIF

!
!4. Plot data
!------------
!
!4.1 Class-'g'
!-------------
!

IF (vecclass.EQ.'g') THEN
   WRITE (iunit,'(A)') '#'

!  ---Cartesian grid
   IF (icoord.EQ.0) THEN
      j_411: DO j=1,nydim
      i_411: DO i=1,nxdim
         IF ((ABS(vecdat(i,j,1)-fill_value).GT.fill_value_eps).AND.&
           & (ABS(vecdat(i,j,2)-fill_value).GT.fill_value_eps)) THEN
            WRITE (iunit,9003) xcoord(i,j), ycoord(i,j),&
                             & vecdat(i,j,1), vecdat(i,j,2)
         ENDIF
      ENDDO i_411
      ENDDO j_411

!   ---spherical grid
   ELSEIF (icoord.GT.0) THEN
      i_412: DO i=1,nxdim
      j_412: DO j=1,nydim
         IF ((ABS(vecdat(i,j,1)-fill_value).GT.fill_value_eps).AND.&
           & (ABS(vecdat(i,j,2)-fill_value).GT.fill_value_eps)) THEN
            WRITE (iunit,9003) ycoord(i,j), xcoord(i,j),&
                             & vecdat(i,j,1), vecdat(i,j,2)
         ENDIF
      ENDDO j_412
      ENDDO i_412
   ENDIF

!
!4.2 Class-'h'
!-------------
!

ELSEIF (vecclass.EQ.'h') THEN
   n_420: DO n=1,2
   j_420: DO j=1,nydim
      WRITE (iunit,9003) vecdat(:,j,n)
   ENDDO j_420
   ENDDO n_420
ENDIF

CALL log_timer_out()


RETURN

9001 FORMAT (2(G15.7,1X),I1)
9002 FORMAT (I1,1X,2I4,7(1X,G15.7),I2)
9003 FORMAT (50(G15.7,1X))

END SUBROUTINE write_mumap_vec
         
!========================================================================

SUBROUTINE write_muplot_par(iunit,ncurves,datfils,istyle,irange,rangemin,&
                          & rangemax,lstyles)
!************************************************************************
!
! *write_muplot_par* Write parameter file used by IDL routine muplot
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - HP_plot, TS_plot, ZP_plot
!
!************************************************************************
!
USE iopars
USE plotpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
INTEGER, INTENT(IN) :: iunit, ncurves
INTEGER, INTENT(IN), DIMENSION(2) :: irange, istyle
INTEGER, INTENT(IN), DIMENSION(ncurves) :: lstyles
CHARACTER (LEN=*), INTENT(IN), DIMENSION(ncurves) :: datfils
REAL, INTENT(IN), DIMENSION(2) :: rangemax, rangemin

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iunit*     INTEGER File unit
!*ncurves*   INTEGER Number of curves
!*datfils*   CHAR    Data file names
!*istyle*    INTEGER Style for axis plots
!*irange*    INTEGER Range for axis plots
!*rangemin*  REAL    Minimum range
!*rangemax*  REAL    Maximum range
!*lstyles*   INTEGER Line styles
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n


procname(pglev+1) = 'write_muplot_par'
CALL log_timer_in()

WRITE (iunit,9001) ncurves, fill_value
n_110: DO n=1,ncurves
   WRITE (iunit,'(A)') TRIM(datfils(n))
ENDDO n_110
n_120: DO n=1,4
   IF (TRIM(plottitles(n)).EQ.'') THEN
      WRITE (iunit,*)
   ELSE
      WRITE (iunit,'(A)') TRIM(plottitles(n))
   ENDIF
ENDDO n_120
n_130: DO n=1,2
   WRITE (iunit,9002) istyle(n), irange(n), rangemin(n), rangemax(n)
ENDDO n_130
WRITE (iunit,'(100I3)') (lstyles(n),n=1,ncurves)

CALL log_timer_out()


RETURN

9001 FORMAT (I3,1X,G15.7)
9002 FORMAT (2I2,2(1X,G15.7))

END SUBROUTINE write_muplot_par

!========================================================================

SUBROUTINE ZP_plot
!************************************************************************
!
! *ZP_plot* Vertical profile plots
!
! Author - Patrick Luyten
!
! Description -
!
! Calling program - post_program
!
! External calls - coords_to_index, post_read, time_locate,
!                  write_muplot_par
!
! Module calls - close_file, convert_loc_to_char, error_alloc, open_file
!
!************************************************************************
!
USE depths
USE grid
USE iopars
USE plotpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: convert_loc_to_char
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=12) :: ccurv, cfig, cipos, cjpos
CHARACTER (LEN=lendesc) :: cloc
CHARACTER (LEN=leniofile) :: paramfile
CHARACTER (LEN=lentime) :: ctime
INTEGER :: icurv, ifig, ilist, ipos, irec, iunit, ivar, jpos, k, kmax, &
         & kmin, numcurves, numstat, vert_unit
REAL :: depmin, depmax, zetdat, zmax
REAL (KIND=kndrlong) :: rtime

CHARACTER (LEN=leniofile), DIMENSION(nolists) :: datfiles
INTEGER, DIMENSION(2) :: irange, istyle
INTEGER, DIMENSION(0:nolists) :: norec
INTEGER, DIMENSION(nolists) :: iodat
INTEGER, DIMENSION(nofigs) :: nocurves
INTEGER, DIMENSION(nofigs,nolists) :: indcurv
REAL, DIMENSION(2) :: rangemin, rangemax
REAL, DIMENSION(nzout) :: profdat, zdat
INTEGER, ALLOCATABLE, DIMENSION(:) :: indlist


procname(pglev+1) = 'ZP_plot'
CALL log_timer_in()

!
!1. Number of curves and curve index array
!-----------------------------------------
!

ifig_110: DO ifig=1,nofigs
   numcurves = 0
   ilist_111: DO ilist=1,nolists
      IF (postatts(ilist)%numfig.EQ.ifig) THEN
         numcurves = numcurves + 1
         indcurv(ifig,numcurves) = ilist
      ENDIF
   ENDDO ilist_111
   nocurves(ifig) = numcurves
ENDDO ifig_110

!
!2. Horizontal location
!----------------------
!

ilist_210: DO ilist=1,nolists

   IF (gridded) THEN
      IF (postatts(ilist)%horz_unit.EQ.1) THEN
         ipos = NINT(postatts(ilist)%gxpos)
         jpos = NINT(postatts(ilist)%gypos)
         postatts(ilist)%gxpos = gxcoord(ipos,jpos)
         postatts(ilist)%gypos = gycoord(ipos,jpos)
         postatts(ilist)%ipos = ipos
         postatts(ilist)%jpos = jpos
      ELSE
         CALL coords_to_index(postatts(ilist)%gxpos,postatts(ilist)%gypos,&
                            & postatts(ilist)%ipos,postatts(ilist)%jpos)
      ENDIF
   ELSE
      numstat = postatts(ilist)%numstat
      postatts(ilist)%gxpos = gxcoord(numstat,1)
      postatts(ilist)%gypos = gycoord(numstat,1)
      postatts(ilist)%ipos = numstat
      postatts(ilist)%jpos = 1
   ENDIF

ENDDO ilist_210

!
!3. Time location
!----------------
!

norec(0) = 0
ilist_310: DO ilist=1,nolists
   CALL time_locate(postatts(ilist)%PDateTime,norec(ilist))
ENDDO ilist_310

!
!4. Create/open files
!--------------------
!

ifig_410: DO ifig=1,nofigs
icurv_410: DO icurv=1,nocurves(ifig)
   ilist = indcurv(ifig,icurv)
   WRITE (cfig,'(I12)') noplots+ifig; cfig = ADJUSTL(cfig)
   WRITE (ccurv,'(I12)') icurv; ccurv = ADJUSTL(ccurv)
   datfiles(ilist) = TRIM(plotfile)//'.'//TRIM(ccurv)//'ZPdat'//TRIM(cfig)
   CALL open_file(iodat(ilist),datfiles(ilist),'OUT','A')
ENDDO icurv_410
ENDDO ifig_410

!
!5. Write data files
!-------------------
!

ilist_510: DO ilist=1,nolists

!  ---read data
   irec_511: DO irec=norec(ilist-1)+1,norec(ilist)
      IF (outform.NE.'N'.OR.irec.EQ.norec(ilist)) THEN
         CALL post_read(irec,ctime,rtime,.TRUE.)
      ENDIF
   ENDDO irec_511

!  ---store data
   ivar = postatts(ilist)%numvar
   ipos = postatts(ilist)%ipos; jpos = postatts(ilist)%jpos
   depmin = postatts(ilist)%depmin
   depmax = postatts(ilist)%depmax
   vert_unit = postatts(ilist)%vert_unit
   zdat = gzcoord(ipos,jpos,:)
   IF (gridded) THEN
      IF (ivar.GT.0) THEN
         profdat = outvals4d(ipos,jpos,:,ivar)
      ELSEIF (ivar.EQ.-1.AND.vcoord.EQ.2) THEN
         profdat = MERGE(gsigcoordatc,gsout(ipos,jpos,:),iopt_grid_vtype.EQ.2)
      ENDIF
   ELSE
      IF (ivar.GT.0) THEN
         profdat = outvals3d(ipos,:,ivar)
      ELSEIF (ivar.EQ.-1.AND.vcoord.EQ.2) THEN
         profdat = MERGE(gsigcoordatc,gsout(ipos,1,:),iopt_grid_vtype.EQ.2)
      ENDIF
   ENDIF
   IF (vert_unit.EQ.1) THEN
      kmin = depmin; kmax = depmax
      zmax = MERGE(depout(ipos,jpos),-zdat(kmin),kmin.EQ.1)
   ELSE
      kmin = 1; kmax = nzout
      zmax = MERGE(depout(ipos,jpos),depmax,depmax.EQ.0.0)
   ENDIF
 
!  ---write data
   IF (vert_unit.EQ.1) THEN
      k_512: DO k=kmin,kmax
         IF (ztime) THEN
            zetdat = gzcoord(ipos,jpos,k) - gzcoord_ini(ipos,jpos,k)
         ELSE
            zetdat = 0.0
         ENDIF
         WRITE (iodat(ilist),9001) profdat(k), zdat(k) + zetdat
      ENDDO k_512
   ELSEIF (vert_unit.EQ.2) THEN
      k_513: DO k=kmin,kmax
         IF (ztime) THEN
            zetdat = gzcoord(ipos,jpos,k) - gzcoord_ini(ipos,jpos,k)
         ELSE
            zetdat = 0.0
         ENDIF
         IF (ABS(zdat(k)).GE.depmin.AND.ABS(zdat(k)).LE.zmax.AND.&
          & (ABS(profdat(k)-fill_value).GT.fill_value_eps)) THEN
            WRITE (iodat(ilist),9001) profdat(k), zdat(k) + zetdat
         ENDIF
      ENDDO k_513
   ENDIF

ENDDO ilist_510

!
!6. Close files
!--------------
!

ilist_610: DO ilist=1,nolists
   CALL close_file(iodat(ilist),'A')
ENDDO ilist_610

!
!7. Write parameter file
!-----------------------
!

ifig_710: DO ifig=1,nofigs

!  ---select curves
   numcurves = nocurves(ifig)
   ALLOCATE (indlist(numcurves),STAT=errstat)
   CALL error_alloc('indlist',1,(/numcurves/),kndint)
   indlist = indcurv(ifig,1:numcurves)
   ilist = indlist(1)

!  ---open parameter file
   WRITE (cfig,'(I12)') noplots + ifig; cfig = ADJUSTL(cfig)
   paramfile = TRIM(plotfile)//'.ZPpar'//TRIM(cfig)
   CALL open_file(iunit,paramfile,'OUT','A')

!  ---axis titles
   ivar = postatts(ilist)%numvar
   plottitles(1) = TRIM(outvars(ivar)%long_name)//' ('//&
                 & TRIM(outvars(ivar)%units)//')'
   plottitles(2) = 'depth (m)'

!  ---main title
   plottitles(3) = TRIM(outfile)

!  ---subtitle
   IF (gridded) THEN
      IF (iopt_grid_sph.EQ.0.OR.postatts(ilist)%horz_unit.EQ.1) THEN
         WRITE (cipos,'(I12)') postatts(ilist)%ipos; cipos = ADJUSTL(cipos)
         WRITE (cjpos,'(I12)') postatts(ilist)%jpos; cjpos = ADJUSTL(cjpos)
         cloc = '('//TRIM(cipos)//','//TRIM(cjpos)//')'
      ELSEIF (iopt_grid_sph.EQ.1) THEN
         cloc(1:21) = convert_loc_to_char(postatts(ilist)%gxpos,&
                                        & postatts(ilist)%gypos)
      ENDIF
   ELSE
      numstat = postatts(ilist)%numstat
      cloc = TRIM(station_names(numstat))
   ENDIF
   plottitles(4) = TRIM(cloc)//'; '//postatts(ilist)%PDatetime(1:19)

!  ---axis styles and ranges
   IF (numcurves.EQ.1) THEN
      istyle = 0; irange = 0
   ELSE
      istyle = 2; irange = 1
   ENDIF
   rangemin = 0; rangemax = 0

!  ---write parameter file
   CALL write_muplot_par(iunit,numcurves,datfiles(indlist),&
                       & istyle,irange,rangemin,rangemax,&
                       & postatts(indlist)%linepsyms)

!  ---close parameter file
   CALL close_file(iunit,'A')

!  ---deallocate
   DEALLOCATE (indlist)

!  ---write to '.vis'-file
   WRITE (io_parvis,9002) 'muplot', TRIM(paramfile), 'y',&
                        & TRIM(outvars(ivar)%f90_name)

ENDDO ifig_710

noplots = noplots + nofigs

CALL log_timer_out()


RETURN

9001 FORMAT (2(G15.7,1X))
9002 FORMAT(4(A,1X))

END SUBROUTINE ZP_plot
