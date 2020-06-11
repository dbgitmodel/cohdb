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
!
!************************************************************************
!
! *Coupling Routines* Routines for coupling different models using MCT
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.11.1
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
! Reference -
!
! Subroutines - coh2wav_coupler, coh2wav_end, coh2wav_grid, coh2wav_params,
!               coh2wav_recv, coh2wav_send, coh2wav_switches, wav2coh_coupler,
!               wav2coh_end, wav2co_grid, wav2coh_params, wav2coh_recv, wav2coh_send,
!               wav2coh_switches, wave_model
!
!************************************************************************
!

!========================================================================

SUBROUTINE coh2wav_coupler
!************************************************************************
!
! *coh2wav_coupler* Initialises the coupling of COHERENS to the wave model
!                   using MCT
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90 V2.9
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - 
!
! Module calls - combine_mod_hgrid_2d, error_alloc, error_alloc_struc,
!                hrelativecoords_init, model_to_data_hcoords
!
!************************************************************************
!
USE couplingvars
USE datatypes
USE grid
USE gridpars
USE iopars
USE paralpars
USE switches
USE syspars
USE wavepars
USE wavevars
USE error_routines, ONLY: error_alloc, error_alloc_struc 
USE datatypes_init, ONLY: hrelativecoords_init
USE grid_interp, ONLY: model_to_data_hcoords
USE paral_comms, ONLY: combine_mod_hgrid_2d
USE time_routines, ONLY: log_timer_in, log_timer_out
#ifdef MCT
   USE m_AttrVect, AttrVect_init => init
   USE m_GlobalSegMap, GlobalSegMap_init  => init
   USE m_MCTWorld, MCTWorld_init => init
   USE m_Router, Router_init => init, Router_print => print
   USE m_SparseMatrix, SparseMatrix_init => init, &
                     & SparseMatrix_importGColInd => importGlobalColumnIndices,&
                     & SparseMatrix_importGRowInd => importGlobalRowIndices, &
                     & SparseMatrix_importMatrixElts => importMatrixElements
   USE m_SparseMatrixPlus, SparseMatrixPlus_init => init
#endif /*MCT*/

!
!*Local variables
!
INTEGER, PARAMETER :: maxlenc = lenname*(MaxCohtoWavvars+1), &
                    & maxlenw = lenname*(MaxWavtoCohvars+1)
CHARACTER (LEN=maxlenc) :: varlist_coh
CHARACTER (LEN=maxlenw) :: varlist_wav
INTEGER :: i, icoord, ivar, j, jcoord, l, n, nocohloc, nocohglb, noelements, &
         & novals, nowavloc, nowavglb, npcc
INTEGER, DIMENSION(1) :: lengthw, startw
INTEGER, DIMENSION(2) :: lbounds
INTEGER, DIMENSION(4) :: idat, jdat
REAL, DIMENSION(2,2) :: whts
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: lengthc, startc
INTEGER, POINTER, DIMENSION(:) :: columns, rows
REAL, POINTER, DIMENSION(:) :: weights
TYPE (HRelativeCoords), SAVE, ALLOCATABLE, DIMENSION(:,:) :: CohtoWavCoords, &
                                                           & CohtoWavCoordsglb
#ifdef MCT

procname(pglev+1) = 'coh2wav_coupler'
CALL log_timer_in(npcc)

!
!1. Initialise
!-------------
!
!1.1 Names of the exchanged data variables
!-----------------------------------------
!
!---COHERENS variables
cohtowav_names(1) = 'rmask'
cohtowav_names(2) = 'depmeanatc'
cohtowav_names(3) = 'zeta'
cohtowav_names(4) = 'umvel'
cohtowav_names(5) = 'vmvel'

!---wave variables
wavtocoh_names(1)  = 'rmask'
wavtocoh_names(2)  = 'waveheight'
wavtocoh_names(3)  = 'waveperiod'
wavtocoh_names(4)  = 'wavedir'
wavtocoh_names(5:13)  = ''
IF (iopt_waves_form.EQ.2) THEN
   wavtocoh_names(5)  = 'wavevel'
   wavtocoh_names(6)  = 'waveexcurs'
ENDIF
IF (iopt_waves_curr.EQ.1.AND.iopt_waves_form.EQ.2) THEN
   wavtocoh_names(7)  = 'umstokesatc'
   wavtocoh_names(8)  = 'vmstokesatc'
   IF (iopt_waves_pres.EQ.1) wavtocoh_names(9) = 'wavepres'
ENDIF
IF (iopt_waves_dissip.EQ.1) THEN
   wavtocoh_names(10) = 'umswdissipatc' 
   wavtocoh_names(11) = 'vmswdissipatc'
   wavtocoh_names(12) = 'umbwdissipatc'
   wavtocoh_names(13) = 'vmbwdissipatc'
ENDIF

!
!1.2 Initialise MCT World
!------------------------
!

CALL MCTWorld_init(nocoupledmodels,comm_world_MPI,comm_model,modelid)

!
!1.3 Allocate arrays
!-------------------
!
!1.3.1 Start and length of segments for COHERENS segmentation map
!----------------------------------------------------------------
!

ALLOCATE (lengthc(nrloc),STAT=errstat)
CALL error_alloc('lengthc',1,(/nrloc/),kndint)
ALLOCATE (startc(nrloc),STAT=errstat)
CALL error_alloc('startc',1,(/nrloc/),kndint)

!
!1.3.2 Arrays for storing sparse matrix elements
!-----------------------------------------------
!

IF (iopt_waves_couple.EQ.1.OR.iopt_waves_couple.EQ.3) THEN

!  --relative coordinates and weights
   ALLOCATE (CohtoWavCoordsglb(nc,nr),STAT=errstat)
   CALL error_alloc_struc('CohtoWavCoordsglb',2,(/nc,nr/),'HRelativeCoords')
   CALL hrelativecoords_init(CohtoWavCoordsglb,.TRUE.)
   ALLOCATE (CohtoWavCoords(ncloc,nrloc),STAT=errstat)
   CALL error_alloc_struc('CohtoWavCoords',2,(/ncloc,nrloc/),'HRelativeCoords')
   CALL hrelativecoords_init(CohtoWavCoords,.TRUE.)

ENDIF

!
!2. Sending data from COHERENS to the wave model
!-----------------------------------------------
!

 IF (iopt_waves_couple.GT.1) THEN

!
!2.1 Names of variables in MCT list format
!-----------------------------------------
!

   varlist_coh = ''
   ivar_210: DO ivar=1,MaxCohtoWavvars
      IF (TRIM(cohtowav_names(ivar)).NE.'') THEN
         varlist_coh = TRIM(varlist_coh)//TRIM(cohtowav_names(ivar))//':'
      ENDIF
   ENDDO ivar_210
   l = LEN_TRIM(varlist_coh)
   varlist_coh(l:l) = ' '
   
!
!2.2 Global segmentation map
!---------------------------
!
!  ---start and length of the segments
   lengthc = ncloc
   j_220: DO j=1,nrloc
      startc(j) = nc1loc + (nr1loc-1)*nc + (j-1)*nc
   ENDDO j_220

!  ---define segmentation map
   CALL GlobalSegMap_init(Gsm_CohSendWav,startc,lengthc,idmaster,comm_model,&
                        & modelid)

!
!2.3 Attribute vector
!--------------------
!

   nocohloc = ncloc*nrloc
   CALL AttrVect_init(Av_CohSendWav,rList=TRIM(varlist_coh),lsize=nocohloc)

ENDIF

!
!3. Receiving data from the wave model to COHERENS
!-------------------------------------------------
!

IF (iopt_waves_couple.EQ.1.OR.iopt_waves_couple.EQ.3) THEN

!
!3.1  Names of variables in MCT list format
!------------------------------------------
!

   varlist_wav = ''
   ivar_310: DO ivar=1,MaxWavtoCohvars
      IF (TRIM(wavtocoh_names(ivar)).NE.'') THEN
         varlist_wav = TRIM(varlist_wav)//TRIM(wavtocoh_names(ivar))//':'
      ENDIF
   ENDDO ivar_310
   l = LEN_TRIM(varlist_wav)
   varlist_wav(l:l) = ' '
      
!
!3.2 Global segmentation map
!---------------------------
!
!  ---start and length of the segments
   ncwavloc = ncwav
   IF (nprocs.GT.1) THEN
      nrwavloc = nrwav/nprocs
      startw = 1 + idloc*ncwavloc*nrwavloc
      IF (idloc.EQ.(nprocs-1)) nrwavloc = nrwav - nrwavloc*(nprocs-1)
   ELSE
      startw = 1
      nrwavloc = nrwav
   ENDIF
   lengthw = ncwavloc*nrwavloc
 
!  ---define segmentation map
   CALL GlobalSegMap_init(Gsm_CohRecvWav,startw,lengthw,idmaster,comm_model,&
                        & modelid)

!
!3.3 Attribute vector
!--------------------
!

   nowavloc = lengthw(1)
   CALL AttrVect_init(Av_CohRecvWav,rList=TRIM(varlist_wav),lsize=nowavloc)

!
!4. Interpolating wave data to the COHERENS grid
!-----------------------------------------------
!
!4.1 Global segmentation map
!---------------------------
!
!  ---start and length of the segments
   lengthc = ncloc
   j_410: DO j=1,nrloc
      startc(j) = nc1loc + (nr1loc-1)*nc + (j-1)*nc
   ENDDO j_410

!  ---define segmentation map
   CALL GlobalSegMap_init(Gsm_CohRecvCoh,startc,lengthc,idmaster,comm_model,&
                        & modelid)

!
!4.2 Attribute vector
!--------------------
!

   novals = ncloc*nrloc
   CALL AttrVect_init(Av_CohRecvCoh,rList=TRIM(varlist_wav),lsize=novals)

!
!4.3 Construct sparse matrix
!---------------------------
!
!4.3.1 Relative coordinates and weights
!--------------------------------------
!
   !  ---local
   CALL model_to_data_hcoords(surfacegrids(igrd_waves,1),CohtoWavCoords,'C  ',&
                            & ncwav,nrwav,gxcoordglbwav,gycoordglbwav,&
                            & maskdat=maskglbwav)
 
!  ---global on master process
   lbounds = 1
   CALL combine_mod_hgrid_2d(CohtoWavCoordsglb,CohtoWavCoords,lbounds,0,0,0.0,&
                           & commall=.FALSE.)
 
!
!4.3.2 Allocate
!--------------
!
!  ---number of non-zero matrix elements
   noelements = 4*COUNT(CohtoWavCoordsglb%icoord.NE.int_fill)
!  ---matrix rows and columns
   IF (master) THEN
      ALLOCATE (rows(noelements),STAT=errstat)
      CALL error_alloc('rows',1,(/noelements/),kndint)
      rows = 0
      ALLOCATE (columns(noelements),STAT=errstat)
      CALL error_alloc('columns',1,(/noelements/),kndint)
      columns = 0
      ALLOCATE (weights(noelements),STAT=errstat)
      CALL error_alloc('weights',1,(/noelements/),kndrtype)
      weights = 0.0
  
!
!4.3.3 Define matrix rows, columns and weights
!---------------------------------------------
!

      n = 0
      j_433: DO j=1,nr
      i_433: DO i=1,nc
         icoord = CohtoWavCoordsglb(i,j)%icoord
         IF (icoord.NE.int_fill) THEN
            jcoord = CohtoWavCoordsglb(i,j)%jcoord
            whts = CohtoWavCoordsglb(i,j)%weights
            idat = (/icoord,icoord+1,icoord+1,icoord/)
            jdat = (/jcoord,jcoord,jcoord+1,jcoord+1/)
            rows(n+1:n+4) = i + nc*(j-1)
            columns(n+1:n+4) = idat + ncwav*(jdat-1)
            weights(n+1:n+4) = (/whts(1,1),whts(2,1),whts(2,2),whts(1,2)/)
            n = n + 4
         ENDIF
      ENDDO i_433
      ENDDO j_433

!
!4.3.4 Construct sparse matrix
!-----------------------------
!
!     ---define
      nocohglb = nc*nr; nowavglb = ncwav*nrwav
      CALL SparseMatrix_init(Sm_WavinttoCoh,nocohglb,nowavglb,noelements)
!     ---rows
      CALL SparseMatrix_importGRowInd(Sm_WavinttoCoh,rows,noelements)
!     ---columns
      CALL SparseMatrix_importGColInd(Sm_WavinttoCoh,columns,noelements)
!     ---insert values
      CALL SparseMatrix_importMatrixElts(Sm_WavinttoCoh,weights,noelements)

   ENDIF

ENDIF

!
!5. Routers and Matrix+ definitions
!----------------------------------
!
!---sending data to the wave model
IF (iopt_waves_couple.GT.1) THEN
   CALL Router_init(modelidwav,Gsm_CohSendWav,comm_model,Rout_CohtoWav)
ENDIF

!---receiving data from the wave model
IF (iopt_waves_couple.EQ.1.OR.iopt_waves_couple.EQ.3) THEN
   CALL Router_init(modelidwav,Gsm_CohRecvWav,comm_model,Rout_WavtoCoh)
   CALL SparseMatrixPlus_init(Smp_WavinttoCoh,Sm_WavinttoCoh,Gsm_CohRecvWav,&
                            & Gsm_CohRecvCoh,Yonly,idmaster,comm_model,modelid)

ENDIF

!
!6. Deallocate
!-------------
!

DEALLOCATE (lengthc,startc)
IF (iopt_waves_couple.EQ.1.OR.iopt_waves_couple.EQ.3) THEN
   DEALLOCATE (CohtoWavCoordsglb,CohtoWavCoords)
   IF (master) DEALLOCATE (rows,columns,weights)
ENDIF

CALL log_timer_out(npcc,itm_MCT)

#endif /*MCT*/

RETURN

END SUBROUTINE coh2wav_coupler

!========================================================================

SUBROUTINE coh2wav_end
!************************************************************************
!
! *coh2wav_end* clean up of MCT types defined by COHERENS  
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.9
!
! Description -
!
! Reference -
!
! Calling program - simulation_end
!
! External calls - 
!
! Module calls - error_abort, write_error_message
!
!************************************************************************
!
USE couplingvars
USE iopars
USE switches
USE error_routines, ONLY: error_abort, write_error_message
USE time_routines, ONLY: log_timer_in, log_timer_out

#ifdef MCT
   USE m_AttrVect, ONLY: AttrVect_clean => clean
   USE m_GlobalSegMap, ONLY: GlobalSegMap_clean  => clean
   USE m_MCTWorld, ONLY: MCTWorld_clean => clean
   USE m_Router, ONLY: Router_clean => clean
   USE m_SparseMatrix, ONLY: SparseMatrix_clean => clean
   USE m_SparseMatrixPlus, ONLY: SparseMatrixPlus_clean => clean
#endif /*MCT*/

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc


#ifdef MCT

procname(pglev+1) = 'coh2wav_end'
CALL log_timer_in(npcc)

!---attribute vectors
IF (iopt_waves_couple.GT.1) THEN
   CALL AttrVect_clean(Av_CohSendWav,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT attribute vector Av_CohSendWav')
   ENDIF
ENDIF
IF (iopt_waves_couple.EQ.1.OR.iopt_waves_couple.EQ.3) THEN
   CALL AttrVect_clean(Av_CohRecvWav,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT attribute vector Av_CohRecvWav')
   ENDIF
   CALL AttrVect_clean(Av_CohRecvCoh,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT attribute vector Av_CohRecvCoh')
   ENDIF
ENDIF

!---global segmentation maps
IF (iopt_waves_couple.GT.1) THEN
   CALL GlobalSegMap_clean(Gsm_CohSendWav,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT global segmentation map Gsm_CohSendWav')
   ENDIF
ENDIF
IF (iopt_waves_couple.EQ.1.OR.iopt_waves_couple.EQ.3) THEN
   CALL GlobalSegMap_clean(Gsm_CohRecvWav,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT global segmentation map Gsm_CohRecvWav')
   ENDIF
   CALL GlobalSegMap_clean(Gsm_CohRecvCoh,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT global segmentation map Gsm_CohRecvCoh')
   ENDIF
ENDIF

!---routers
IF (iopt_waves_couple.GT.1) THEN
   CALL Router_clean(Rout_CohtoWav,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message('Unable to clean MCT router Rout_CohtoWav')
   ENDIF
ENDIF
IF (iopt_waves_couple.EQ.1.OR.iopt_waves_couple.EQ.3) THEN
   CALL Router_clean(Rout_WavtoCoh,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message('Unable to clean MCT router Rout_WavtoCoh')
   ENDIF
ENDIF

!---sparse matrix
IF (iopt_waves_couple.EQ.1.OR.iopt_waves_couple.EQ.3) THEN
   CALL SparseMatrix_clean(Sm_WavinttoCoh,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT sparse matrix Sm_WavinttoCoh')
   ENDIF
ENDIF

!---sparse matrix plus
IF (iopt_waves_couple.EQ.1.OR.iopt_waves_couple.EQ.3) THEN
   CALL SparseMatrixPlus_clean(Smp_WavinttoCoh,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT sparse matrix plus Smp_WavinttoCoh')
   ENDIF
ENDIF

!---MCT world
CALL MCTWorld_clean()

CALL error_abort('coh2wav_end',ierrno_MCT)

CALL log_timer_out(npcc,itm_MCT)

#endif /*MCT*/


RETURN

END SUBROUTINE coh2wav_end

!========================================================================

SUBROUTINE coh2wav_grid
!************************************************************************
!
! *coh2wav_grid* receive wave grid from the wave master, send COHERENS grid to
!                wave master
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.9
!
! Description - called by all processes within the COHERENS communicator
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - 
!
! Module calls - check_date, check_time_limits_arr_struc, comms_recv_char,
!                comms_recv_log, comms_recv_int, comms_recv_real,
!                comms_send_int, comms_send_log, comms_send_real, copy_chars,
!                copy_vars, error_abort, error_alloc, num_time_steps
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
USE wavepars
USE wavevars
USE error_routines, ONLY: error_alloc
USE comms_MPI, ONLY: comms_recv_log, comms_recv_real, comms_send_log, &
                   & comms_send_real
USE paral_comms, ONLY: copy_vars
USE time_routines, ONLY: log_timer_in, log_timer_out


IMPLICIT NONE

!
!*Local variables
!
INTEGER :: itag, novals, npcc


procname(pglev+1) = 'coh2wav_grid'
CALL log_timer_in(npcc)

!
!1. Allocate
!-----------
!

ALLOCATE (gxcoordglbwav(ncwav,nrwav),STAT=errstat)
CALL error_alloc('gxcoordglbwav',2,(/ncwav,nrwav/),kndrtype)
ALLOCATE (gycoordglbwav(ncwav,nrwav),STAT=errstat)
CALL error_alloc('gycoordglbwav',2,(/ncwav,nrwav/),kndrtype)
ALLOCATE (maskglbwav(ncwav,nrwav),STAT=errstat)
CALL error_alloc('maskglbwav',2,(/ncwav,nrwav/),kndlog)

!
!2. Receive wave grid from wave master
!-------------------------------------
!

IF (master) THEN
   novals = ncwav*nrwav
   itag = 1000
   CALL comms_recv_real(gxcoordglbwav,novals,idmasterwav,itag,comm_world_MPI,&
                      & iarr_gxcoordglbwav)
   itag = itag + 1
   CALL comms_recv_real(gycoordglbwav,novals,idmasterwav,itag,comm_world_MPI,&
                      & iarr_gycoordglbwav)
   itag = itag + 1
   CALL comms_recv_log(maskglbwav,novals,idmasterwav,itag,comm_world_MPI,&
                     & iarr_maskglbwav)

!
!3. Send COHERENS grid to wave master
!------------------------------------
!
   itag = itag + 1
   novals = (nc+2)*(nr+2)
   CALL comms_send_real(gxcoordglb,novals,idmasterwav,itag,&
                      & comm_world_MPI,1,iarr_gxcoordglb)
   itag = itag + 1
   CALL comms_send_real(gycoordglb,novals,idmasterwav,itag,&
                      & comm_world_MPI,1,iarr_gycoordglb)
   itag = itag + 1
   CALL comms_send_log(seapointglb,novals,idmasterwav,itag,&
                     & comm_world_MPI,1,iarr_seapointglb)

ENDIF


!
!3. Broadcast within the COHERENS communicator
!---------------------------------------------
!

IF (iopt_MPI.EQ.1) THEN

   CALL copy_vars(gxcoordglbwav,iarr_gxcoordglbwav)
   CALL copy_vars(gycoordglbwav,iarr_gycoordglbwav)
   CALL copy_vars(maskglbwav,iarr_maskglbwav)

ENDIF

CALL log_timer_out(npcc,itm_MCT)


RETURN

END SUBROUTINE coh2wav_grid

!========================================================================

SUBROUTINE coh2wav_params
!************************************************************************
!
! *coh2wav_params* wave model parameters received by the wave model
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.9
!
! Description - called by all processes within the COHERENS communicator
!
! Reference -
!
! Calling program - initialise_model 
!
! External calls - 
!
! Module calls - check_date, check_time_limits_arr_struc, comms_recv_char,
!                comms_recv_int, comms_recv_real, comms_send_int,
!                comms_send_log, comms_send_real, copy_chars, copy_vars,
!                error_abort, error_alloc, num_time_steps, write_error_message
!
!************************************************************************
!
USE gridpars
USE iopars
USE paralpars
USE switches
USE timepars
USE wavepars
USE comms_MPI, ONLY: comms_recv_int, comms_recv_real, comms_send_int, &
                   & comms_send_log, comms_send_real
USE error_routines, ONLY: check_time_limits_arr_struc, error_abort, &
                        & write_error_message
USE paral_comms, ONLY: copy_vars
USE time_routines, ONLY: check_date, convert_date, log_timer_in, &
                       & log_timer_out, num_time_steps

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: coordtype, itag, mdelt, mdeltwav, npcc
LOGICAL, DIMENSION(1) :: logsnd
INTEGER, DIMENSION(4) :: intrcv
INTEGER, DIMENSION(3) :: intsnd
REAL, DIMENSION(5) :: realrcv
REAL, DIMENSION(8) :: realsnd


procname(pglev+1) = 'coh2wav_params'
CALL log_timer_in(npcc)

!
!1. Receive data from the wave master
!------------------------------------
!

IF (master) THEN

!  ---start/end date
   itag = 1000
   CALL comms_recv_int(IStartDateTimeWav(1:6),6,idmasterwav,itag,&
                     & comm_world_MPI,0)
   IStartDateTimeWav(7) = 0
   itag = itag + 1
   CALL comms_recv_int(IEndDateTimeWav(1:6),6,idmasterwav,itag,comm_world_MPI,0)
   IEndDateTimeWav(7) = 0

!  ---integer parameters
   itag = itag + 1
   CALL comms_recv_int(intrcv,4,idmasterwav,itag,comm_world_MPI,0)

!  ---real parameters
   itag = itag + 1
   CALL comms_recv_real(realrcv,5,idmasterwav,itag,comm_world_MPI,0)

!
!2. Send model grid parameters to wave master
!--------------------------------------------
!

   intsnd(1) = nc; intsnd(2) = nr; intsnd(3) = iopt_grid_htype
   itag = itag + 1
   CALL comms_send_int(intsnd,3,idmasterwav,itag,comm_world_MPI,1,0)
   logsnd(1) = surfacegrids(igrd_model,1)%rotated   
   itag = itag + 1
   CALL comms_send_log(logsnd,1,idmasterwav,itag,comm_world_MPI,1,0)
   realsnd(1) = surfacegrids(igrd_model,1)%x0dat
   realsnd(2) = surfacegrids(igrd_model,1)%y0dat
   realsnd(3) = surfacegrids(igrd_model,1)%delxdat
   realsnd(4) = surfacegrids(igrd_model,1)%delydat
   realsnd(5) = surfacegrids(igrd_model,1)%gridangle
   realsnd(6) = surfacegrids(igrd_model,1)%longpole
   realsnd(7) = surfacegrids(igrd_model,1)%x0rot
   realsnd(8) = surfacegrids(igrd_model,1)%y0rot
   itag = itag + 1
   CALL comms_send_real(realsnd,8,idmasterwav,itag,comm_world_MPI,1,0) 
ENDIF

!
!3. Broadcast within the COHERENS communicator
!---------------------------------------------
!

IF (nprocs.GT.1) THEN

   CALL copy_vars(idmasterwav,0)
   CALL copy_vars(IStartDateTimeWav,0)
   CALL copy_vars(IEndDateTimeWav,0)
   CALL copy_vars(intrcv,0)
   CALL copy_vars(realrcv,0)

ENDIF

!
!4. Wave model parameters
!------------------------
!
!4.1 Wave grid
!-------------
!
!---grid sizes
ncwav = intrcv(1); nrwav = intrcv(2)

!---surface grid parameters
surfacegrids(igrd_waves,1)%n1dat = ncwav
surfacegrids(igrd_waves,1)%n2dat = nrwav
surfacegrids(igrd_waves,1)%nhtype = intrcv(3)
surfacegrids(igrd_waves,1)%x0dat = realrcv(2)
surfacegrids(igrd_waves,1)%y0dat = realrcv(3)
surfacegrids(igrd_waves,1)%delxdat = realrcv(4)
surfacegrids(igrd_waves,1)%delydat = realrcv(5)

!---type of coordinates
coordtype = intrcv(4)
IF (.NOT.(coordtype.EQ.0.EQV.iopt_grid_sph.EQ.0)) THEN
   CALL write_error_message('COHERENS and wave model must use the same '//&
        & 'type of coordinates (Cartesian or spherical)')
ENDIF

!
!4.2 Start/end dates
!-------------------
!

CStartDateTimeWav = convert_date(IStartDateTimeWav)
CEndDateTimeWav = convert_date(IEndDateTimeWav) 

!---check
IF (master) THEN
   CALL check_date(CStartDateTimeWav,'CStartDateTimeWav')
   CALL check_date(CEndDateTimeWav,'CEndDateTimeWav')
ENDIF

!
!4.3 Time steps
!--------------
!
!---time step
deltwav = realrcv(1)

!---start/end coupling time indices
CALL num_time_steps(CStartDateTime,CStartDateTimeWav,delt2d,&
                  & modfiles(io_wavsur,1,1)%tlims(1))
CALL num_time_steps(CStartDateTime,CEndDateTimeWav,delt2d,&
                  & modfiles(io_wavsur,1,1)%tlims(2))

!---time counter
mdelt = 1000*delt2d
mdeltwav = 1000*deltwav
modfiles(io_wavsur,1,1)%tlims(3) = mdeltwav/mdelt
IF (master) THEN
   CALL check_time_limits_arr_struc(modfiles(io_wavsur,1,1)%tlims,'modfiles',&
                                 & 'tlims',0,nstep,3,(/io_wavsur,1,1/),&
                                 & flag3d=.TRUE.)
ENDIF

CALL error_abort('coh2wav_params',ierrno_inival)

CALL log_timer_out(npcc,itm_MCT)


RETURN

END SUBROUTINE coh2wav_params

!========================================================================

SUBROUTINE coh2wav_recv
!************************************************************************
!
! *coh2wav_recv* receive data from the wave model for Coherens
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.9
!
! Description -
!
! Reference -
!
! Calling program - wave_input 
!
! External calls - 
!
! Module calls - rotate_vec
!
!************************************************************************
!
USE couplingvars
USE grid
USE gridpars
USE iopars
USE paralpars
USE switches
USE wavepars
USE wavevars
USE grid_routines, ONLY: rotate_vec
USE time_routines, ONLY: log_timer_in, log_timer_out

#ifdef MCT
   USE m_AttrVect,ONLY : AttrVect_indexR  => indexRA
   USE m_Transfer, ONLY : MCT_Recv => recv
   USE m_MatAttrVectMul, ONLY: MCT_MatVecMul => sMatAvMult
#endif /*MCT*/

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iivar, ivar, npcc
REAL, DIMENSION(ncloc,nrloc) :: array2du1, array2du2, array2dv1, array2dv2, &
                              & rmask



IF (iopt_waves_couple.EQ.2) RETURN

#ifdef MCT

procname(pglev+1) = 'coh2wav_recv'
CALL log_timer_in(npcc)

!
!1. Receive data on the wave grid
!--------------------------------
!

CALL MCT_Recv(Av_CohRecvWav,Rout_WavtoCoh,Tag=itagwavtocoh)

!
!2. Interpolate onto the Coherens grid
!-------------------------------------
!

CALL MCT_MatVecMul(Av_CohRecvWav,Smp_WavinttoCoh,Av_CohRecvCoh)

!
!3. Retrieve data from MCT attribute vector
!------------------------------------------
!
!3.1 Land mask
!-------------
!

ivar = AttrVect_indexR(Av_CohRecvCoh,'rmask')
rmask  =  RESHAPE(Av_CohRecvCoh%rATTr(ivar,:),(/ncloc,nrloc/))
WHERE (rmask.GT.0.0) 
   rmask = 1.0/rmask
ELSEWHERE
   rmask = 0.0
END WHERE

!
!3.2 Stokes velocities and wave induced pressure
!-----------------------------------------------
!

IF (iopt_waves_curr.EQ.1) THEN

!  ---retrieve data
   ivar_320: DO ivar=1,MaxWavtoCohvars
      SELECT CASE (TRIM(wavtocoh_names(ivar)))
         CASE ('umstokesatc')
            iivar = AttrVect_indexR(Av_CohRecvCoh,'umstokesatc')
            array2du1 = rmask*RESHAPE(Av_CohRecvCoh%rATTr(iivar,:),&
                                   & (/ncloc,nrloc/))
         CASE ('vmstokesatc')
            iivar = AttrVect_indexR(Av_CohRecvCoh,'vmstokesatc')
            array2dv1 = rmask*RESHAPE(Av_CohRecvCoh%rATTr(iivar,:),&
                                   & (/ncloc,nrloc/))

         CASE ('wavepres')
            iivar = AttrVect_indexR(Av_CohRecvCoh,'wavepres')
            wavepres(1:ncloc,1:nrloc) = rmask*&
                         RESHAPE(Av_CohRecvCoh%rATTr(iivar,:),(/ncloc,nrloc/))
      END SELECT
   ENDDO ivar_320

!  ---rotate data if needed
   IF (rotate_gvecs) THEN
      CALL rotate_vec(array2du1,array2dv1,umstokesatc(1:ncloc,1:nrloc),&
                    & vmstokesatc(1:ncloc,1:nrloc),gangleatc(1:ncloc,1:nrloc),1)
   ELSE
      umstokesatc(1:ncloc,1:nrloc) = array2du1
      vmstokesatc(1:ncloc,1:nrloc) = array2dv1
   ENDIF

ENDIF

!
!3.3 Wave dissipation forces
!---------------------------
!

IF (iopt_waves_dissip.EQ.1) THEN

!  ---retrieve surface data
   ivar_331: DO ivar=1,MaxWavtoCohvars
      SELECT CASE (TRIM(wavtocoh_names(ivar)))
         CASE ('umswdissipatc')
            iivar = AttrVect_indexR(Av_CohRecvCoh,'umswdissipatc')
            array2du1 = rmask*RESHAPE(Av_CohRecvCoh%rATTr(iivar,:),&
                                   & (/ncloc,nrloc/))
         CASE ('vmswdissipatc')
            iivar = AttrVect_indexR(Av_CohRecvCoh,'vmswdissipatc')
            array2dv1 = rmask*RESHAPE(Av_CohRecvCoh%rATTr(iivar,:),&
                                   & (/ncloc,nrloc/))
      END SELECT
   ENDDO ivar_331

!  ---retrieve bottom data
   ivar_332: DO ivar=1,MaxWavtoCohvars
      SELECT CASE (TRIM(wavtocoh_names(ivar)))
         CASE ('umbwdissipatc')
            iivar = AttrVect_indexR(Av_CohRecvCoh,'umbwdissipatc')
            array2du2 = rmask*RESHAPE(Av_CohRecvCoh%rATTr(iivar,:),&
                                   & (/ncloc,nrloc/))
         CASE ('vmbwdissipatc')
            iivar = AttrVect_indexR(Av_CohRecvCoh,'vmbwdissipatc')
            array2dv2 = rmask*RESHAPE(Av_CohRecvCoh%rATTr(iivar,:),&
                                   & (/ncloc,nrloc/))
      END SELECT
   ENDDO ivar_332

!  ---rotate data if needed
   IF (rotate_gvecs) THEN
      CALL rotate_vec(array2du1,array2dv1,umswdissipatc(1:ncloc,:),&
                    & vmswdissipatc(:,1:nrloc),gangleatc(1:ncloc,1:nrloc),1)
      CALL rotate_vec(array2du2,array2dv2,umbwdissipatc(1:ncloc,:),&
                    & vmbwdissipatc(:,1:nrloc),gangleatc(1:ncloc,1:nrloc),1)
   ELSE
      umswdissipatc(1:ncloc,:) = array2du1; vmswdissipatc(:,1:nrloc) = array2dv1
      umbwdissipatc(1:ncloc,:) = array2du2; vmbwdissipatc(:,1:nrloc) = array2dv2
   ENDIF
   
ENDIF

!
!3.4 Other spectrally averaged quantities
!----------------------------------------
!

ivar_350: DO ivar=1,MaxWavtoCohvars

   SELECT CASE (TRIM(wavtocoh_names(ivar)))
      CASE ('waveheight')
         iivar = AttrVect_indexR(Av_CohRecvCoh,'waveheight')
         waveheight(1:ncloc,1:nrloc) = rmask*&
                       & RESHAPE(Av_CohRecvCoh%rATTr(iivar,:),(/ncloc,nrloc/))
      CASE ('waveperiod')
         iivar = AttrVect_indexR(Av_CohRecvCoh,'waveperiod')
         waveperiod(1:ncloc,1:nrloc) = rmask*&
                       & RESHAPE(Av_CohRecvCoh%rATTr(iivar,:),(/ncloc,nrloc/))
      CASE ('wavedir')
         iivar = AttrVect_indexR(Av_CohRecvCoh,'wavedir')
         wavedir(1:ncloc,1:nrloc) = rmask*&
                       & RESHAPE(Av_CohRecvCoh%rATTr(iivar,:),(/ncloc,nrloc/))
      CASE ('wavevel')
         iivar = AttrVect_indexR(Av_CohRecvCoh,'wavevel')
         wavevel(1:ncloc,1:nrloc) = rmask*&
                       & RESHAPE(Av_CohRecvCoh%rATTr(iivar,:),(/ncloc,nrloc/))
      CASE ('waveexcurs')
         iivar = AttrVect_indexR(Av_CohRecvCoh,'waveexcurs')
         waveexcurs(1:ncloc,1:nrloc) = rmask*&
                       & RESHAPE(Av_CohRecvCoh%rATTr(iivar,:),(/ncloc,nrloc/))
   END SELECT

ENDDO ivar_350

CALL log_timer_out(npcc,itm_MCT)

#endif /*MCT*/


RETURN

END SUBROUTINE coh2wav_recv

!========================================================================

SUBROUTINE coh2wav_send
!************************************************************************
!
! *coh2wav_send* send data from Coherens to the wave model
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.9
!
! Description -
!
! Reference -
!
! Calling program - coherens_main, initialise_model
!
! External calls - 
!
! Module calls - Uarr_at_C, Varr_at_C
!
!************************************************************************
!
USE couplingvars
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE switches
USE syspars
USE wavepars
USE array_interp, ONLY:Uarr_at_C, Varr_at_C
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: rotate_vec
USE time_routines, ONLY: log_timer_in, log_timer_out

#ifdef MCT
   USE m_Transfer,ONLY : MCT_Send => send
   USE m_AttrVect,ONLY : AttrVect_indexR  => indexRA
#endif /*MCT*/

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, iivar, ivar, j, n, npcc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: rmask, umvelatc, umveltmp, &
                                       & umvelwav, vmvelatc, vmveltmp, vmvelwav 


IF (iopt_waves_couple.EQ.1) RETURN

#ifdef MCT

procname(pglev+1) = 'coh2wav_send'
CALL log_timer_in(npcc)

!
!1. Initialise
!-------------
!
!---allocate
ALLOCATE(rmask(ncloc,nrloc), STAT=errstat)
CALL error_alloc('rmask',2,(/ncloc,nrloc/),kndrtype)
rmask = MERGE(1.0,0.0,maskatc_int)

ALLOCATE(umvelatc(ncloc,nrloc), STAT=errstat)
CALL error_alloc('umvelatc',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE(umveltmp(ncloc,nrloc), STAT=errstat)
CALL error_alloc('umveltmp',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE(umvelwav(ncloc,nrloc), STAT=errstat)
CALL error_alloc('umvelwav',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE(vmvelatc(ncloc,nrloc), STAT=errstat)
CALL error_alloc('vmvelatc',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE(vmveltmp(ncloc,nrloc), STAT=errstat)
CALL error_alloc('vmveltmp',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE(vmvelwav(ncloc,nrloc), STAT=errstat)
CALL error_alloc('vmvelwav',2,(/ncloc,nrloc/),kndrtype)

!---2-D currents
CALL Uarr_at_C(umvel(1:ncloc+1,1:nrloc),umveltmp,1,1,&
            & (/1,1,1/),(/ncloc+1,nrloc,1/),1,iarr_umvel,.TRUE.)
CALL Varr_at_C(vmvel(1:ncloc,1:nrloc+1),vmveltmp,1,1,&
            & (/1,1,1/),(/ncloc,nrloc+1,1/),1,iarr_vmvel,.TRUE.)

!---rotate if necessary
IF (rotate_gvecs) THEN
   CALL rotate_vec(umveltmp,vmveltmp,umvelwav,vmvelwav,&
                 & gangleatc(1:ncloc,1:nrloc),-1)
ELSE
   umvelwav = umveltmp
   vmvelwav = vmveltmp
ENDIF

!
!2. Add data to MCT attribute vector
!-----------------------------------
!

ivar_210: DO ivar=1,MaxCohtoWavvars

!  ---array index within MCT attribute vector 
   SELECT CASE (TRIM(cohtowav_names(ivar)))
      CASE ('rmask')
         iivar = AttrVect_indexR(Av_CohSendWav,'rmask')
      CASE ('depmeanatc')
         iivar = AttrVect_indexR(Av_CohSendWav,'depmeanatc')
      CASE ('zeta')
         iivar = AttrVect_indexR(Av_CohSendWav,'zeta')
      CASE ('umvel')
         iivar = AttrVect_indexR(Av_CohSendWav,'umvel')
      CASE ('vmvel')
         iivar = AttrVect_indexR(Av_CohSendWav,'vmvel')
   END SELECT

!  ---store values
   n = 0
   j_211: DO j=1,nrloc
   i_211: DO i=1,ncloc
      n = n + 1
      SELECT CASE (TRIM(cohtowav_names(ivar)))
         CASE ('rmask')
            Av_CohSendWav%rATTr(iivar,n) = rmask(i,j)
         CASE ('depmeanatc')
            Av_CohSendWav%rATTr(iivar,n) = rmask(i,j)*depmeanatc(i,j)
         CASE ('zeta')
            Av_CohSendWav%rATTr(iivar,n) = rmask(i,j)*zeta(i,j)
         CASE ('umvel')
            Av_CohSendWav%rATTr(iivar,n) = rmask(i,j)*umvelwav(i,j)
         CASE ('vmvel')
            Av_CohSendWav%rATTr(iivar,n) = rmask(i,j)*vmvelwav(i,j)
      END SELECT
   ENDDO i_211
   ENDDO j_211
ENDDO ivar_210

!
!3. Send data
!------------
!

CALL MCT_Send(Av_CohSendWav,Rout_CohtoWav,Tag=itagcohtowav)

!
!4. Deallocate
!-------------
!

DEALLOCATE(rmask,umvelatc,umveltmp,umvelwav,vmvelatc,vmveltmp,vmvelwav)

CALL log_timer_out(npcc,itm_MCT)

#endif /*MCT*/


RETURN

END SUBROUTINE coh2wav_send

!========================================================================

SUBROUTINE coh2wav_switches
!************************************************************************
!
! *coh2wav_switches* define which set of variables are used for coupling
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.11.1
!
! Description - send model switches to the wave model
!
! Reference -
!
! Calling program - initialise_model 
!
! External calls - 
!
! Module calls - comms_send_int
!
!************************************************************************
!
USE iopars
USE paralpars
USE switches
USE comms_MPI, ONLY: comms_send_int
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: itag, npcc
INTEGER, DIMENSION(6) :: intdat


procname(pglev+1) = 'coh2wav_switches'
CALL log_timer_in(npcc)

!
!1. Define
!---------
!

intdat = (/iopt_waves_couple,iopt_bstres_waves_bfric,iopt_waves_curr,&
         & iopt_waves_dissip,iopt_waves_form,iopt_waves_pres/)

!
!2. Send switches to wave master
!-------------------------------
!

IF (master) THEN
   itag = 1000
   CALL comms_send_int(intdat,6,idmasterwav,itag,comm_world_MPI,1,0)
ENDIF

CALL log_timer_out(npcc,itm_MCT)


RETURN

END SUBROUTINE coh2wav_switches

!========================================================================

SUBROUTINE wav2coh_coupler(nxlocext,nylocext,nx1loc,ny1loc,nhdims)
!************************************************************************
!
! *wav2coh_coupler* initialise the coupling of the wave model to COHERENS
!                   using MCT
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90 V2.9
!
! Description - routine is called by all wave processes
!
! Reference -
!
! Calling program - wave model
!
! External calls - 
!
! Module calls - data_to_model_hcoords, error_alloc, error_alloc_struc,
!                hrelativecoords_init, monitor_files
!
!************************************************************************
!
USE couplingvars
USE datatypes
USE gridpars
USE iopars
USE paralpars
USE switches
USE syspars
USE timepars
USE wavepars
USE wavevars
USE error_routines, ONLY: error_alloc, error_alloc_struc 
USE datatypes_init, ONLY: hrelativecoords_init
USE inout_routines, ONLY: monitor_files
USE grid_interp, ONLY: data_to_model_hcoords
USE time_routines, ONLY: log_timer_in, log_timer_out
#ifdef MCT
   USE m_AttrVect, AttrVect_init => init
   USE m_GlobalSegMap, GlobalSegMap_init  => init
   USE m_MCTWorld, MCTWorld_init => init
   USE m_Router, Router_init => init, Router_print => print
   USE m_SparseMatrix, SparseMatrix_init => init, &
                     & SparseMatrix_importGColInd => importGlobalColumnIndices,&
                     & SparseMatrix_importGRowInd => importGlobalRowIndices, &
                     & SparseMatrix_importMatrixElts => importMatrixElements
   USE m_SparseMatrixPlus, SparseMatrixPlus_init => init
#endif /*MCT*/

!
!*Arguments
!
INTEGER, INTENT(IN) :: nxlocext, nx1loc, nylocext,ny1loc
INTEGER, DIMENSION(4), INTENT(IN) :: nhdims

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*nxlocext*  INTEGER  X-dimension of the local wave domain including halos
!*nylocext*  INTEGER  Y-dimension of the local wave domain including halos
!*nx1loc*    INTEGER  Global index of the first column (including halo) on the
!                     local wave domain grid
!*ny1loc*    INTEGER  Global index of the first row (including halo) on the
!                     local wave domain grid
!*nhdims*    INTEGER  Halo sizes of the local wave domain (WESN)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER, PARAMETER :: maxlenc = lenname*(MaxCohtowavvars+1), &
                    & maxlenw = lenname*(MaxWavtoCohvars+1)
CHARACTER (LEN=maxlenc) :: varlist_coh
CHARACTER (LEN=maxlenw) :: varlist_wav
INTEGER :: i, icoord, ivar, j, jcoord, l, n, nocohloc, nocohglb, noelements, &
         & novals, nowavloc, nowavglb, npcc
INTEGER, DIMENSION(1) :: lengthc, startc
INTEGER, DIMENSION(4) :: idat, jdat
REAL, DIMENSION(2,2) :: whts
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: lengthw, startw
INTEGER, POINTER, DIMENSION(:) :: columns, rows
REAL, POINTER, DIMENSION(:) :: weights
TYPE (HRelativeCoords), SAVE , ALLOCATABLE, DIMENSION(:,:) :: WavtoCohCoordsglb


#ifdef MCT

procname(pglev+1) = 'wav2coh_coupler'
CALL log_timer_in(npcc)

!
!1. Initialise
!-------------
!
!1.1 Grid parameters
!-------------------
!

ncwavloc_ext = nxlocext; nrwavloc_ext = nylocext
haloW = nhdims(1); haloE = nhdims(2)
haloS = nhdims(3); haloN = nhdims(4)
nc1wavloc = nx1loc + haloW; nr1wavloc = ny1loc + haloS
ncwavloc = ncwavloc_ext - haloW - haloE 
nrwavloc = nrwavloc_ext - haloS - haloN
nc2wavloc = nc1wavloc + ncwavloc - 1 
nr2wavloc = nr1wavloc + nrwavloc - 1 

!
!1.2 Names of the exchanged data variables
!-----------------------------------------
!
!---COHERENS variable
cohtowav_names(1) = 'rmask'
cohtowav_names(2) = 'depmeanatc'
cohtowav_names(3) = 'zeta'
cohtowav_names(4) = 'umvel'
cohtowav_names(5) = 'vmvel'

!---wave variables
wavtocoh_names(1)  = 'rmask'
wavtocoh_names(2)  = 'waveheight'
wavtocoh_names(3)  = 'waveperiod'
wavtocoh_names(4)  = 'wavedir'
wavtocoh_names(5:13)  = ''
IF (iopt_waves_form.EQ.2) THEN
   wavtocoh_names(5)  = 'wavevel'
   wavtocoh_names(6)  = 'waveexcurs'
ENDIF
IF (iopt_waves_curr.EQ.1.AND.iopt_waves_form.EQ.2) THEN
   wavtocoh_names(7)  = 'umstokesatc'
   wavtocoh_names(8)  = 'vmstokesatc'
   IF (iopt_waves_pres.EQ.1) wavtocoh_names(9)  = 'wavepres'
ENDIF
IF (iopt_waves_dissip.EQ.1) THEN
   wavtocoh_names(10) = 'umswdissipatc' 
   wavtocoh_names(11) = 'vmswdissipatc'
   wavtocoh_names(12) = 'umbwdissipatc'
   wavtocoh_names(13) = 'vmbwdissipatc'
ENDIF

!
!1.3 Initialise MCT World
!------------------------
!

CALL MCTWorld_init(nocoupledmodels,comm_world_MPI,comm_model,modelid)

!
!2. Receiving data from COHERENS to the wave model
!-------------------------------------------------
!

IF (iopt_waves_couple.GT.1) THEN
 
  IF (mastermod) THEN
      ALLOCATE (rows(noelements),STAT=errstat)
      CALL error_alloc('rows',1,(/noelements/),kndint)
      rows = 0
      ALLOCATE (columns(noelements),STAT=errstat)
      CALL error_alloc('columns',1,(/noelements/),kndint)
      columns = 0
      ALLOCATE (weights(noelements),STAT=errstat)
      CALL error_alloc('weights',1,(/noelements/),kndrtype)
      weights = 0.0
   ENDIF

!
!2.1  Names of variables in MCT list format
!------------------------------------------
!

   varlist_coh = ''
   ivar_210: DO ivar=1,MaxCohtoWavvars
      IF (TRIM(cohtowav_names(ivar)).NE.'') THEN
         varlist_coh = TRIM(varlist_coh)//TRIM(cohtowav_names(ivar))//':'
      ENDIF
   ENDDO ivar_210
   l = LEN_TRIM(varlist_coh)
   varlist_coh(l:l) = ' '

!
!2.2 Global segmentation map
!---------------------------
!
!  ---allocate
   ALLOCATE (lengthw(nrwavloc),STAT=errstat)
   CALL error_alloc('lengthw',1,(/nrwavloc/),kndint)
   ALLOCATE (startw(nrwavloc),STAT=errstat)
   CALL error_alloc('startw',1,(/nrwavloc/),kndint)

!  ---start and length of the segments
   ncloc = nc
   IF (nprocswav.GT.1) THEN
      nrloc = nr/nprocswav
      startc = 1 + idloc*ncloc*nrloc
      IF (idloc.EQ.(nprocswav-1)) nrloc = nr - nrloc*(nprocswav-1)
   ELSE
      startc = 1
      nrloc = nr
   ENDIF
   lengthc = ncloc*nrloc
   
!  ---define segmentation map
   CALL GlobalSegMap_init(Gsm_WavRecvCoh,startc,lengthc,idmastermod,comm_model,&
                        & modelid)

!
!2.3 Attribute vector
!--------------------
!

   nocohloc = lengthc(1)
   CALL AttrVect_init(Av_WavRecvCoh,rList=TRIM(varlist_coh),lsize=nocohloc)

!
!2.4 Deallocate
!--------------
!

   DEALLOCATE (lengthw,startw)
   
ENDIF
   
!
!3. Interpolating COHERENS data to the wave grid
!-----------------------------------------------
!

IF (iopt_waves_couple.GT.1) THEN

!
!3.1 Global segmentation map
!---------------------------
!
!  ---allocate
   ALLOCATE (lengthw(nrwavloc_ext),STAT=errstat)
   CALL error_alloc('lengthw',1,(/nrwavloc_ext/),kndint)
   ALLOCATE (startw(nrwavloc_ext),STAT=errstat)
   CALL error_alloc('startw',1,(/nrwavloc_ext/),kndint)

!  ---start and length of the segments
   lengthw = ncwavloc_ext
   j_310: DO j=1,nrwavloc_ext
      startw(j) = nc1wavloc - haloW + (nr1wavloc-1)*ncwav + (j-1-haloS)*ncwav
   ENDDO j_310

   !  ---define segmentation map
   nowavglb = ncwav*nrwav
   CALL GlobalSegMap_init(Gsm_WavRecvWav,startw,lengthw,idmastermod,comm_model,&
                        & modelid,gsize=nowavglb)

!
!3.2 Attribute vector
!--------------------
!

   novals = ncwavloc_ext*nrwavloc_ext
   CALL AttrVect_init(Av_WavRecvWav,rList=TRIM(varlist_coh),lsize=novals)

!
!3.3 Construct sparse matrix
!---------------------------
!
!3.3.1 Relative coordinates and weights
!--------------------------------------
!
!  ---allocate
   ALLOCATE (WavtoCohCoordsglb(ncwav,nrwav),STAT=errstat)
   CALL error_alloc_struc('WavtoCohCoordsglb',2,(/ncwav,nrwav/),& 
                        & 'HRelativeCoords')
   CALL hrelativecoords_init(WavtoCohCoordsglb,.TRUE.)

!  ---global on all procces in wave communicator
   CALL data_to_model_hcoords(WavtoCohCoordsglb,'C  ',ncwav,nrwav,&
        & gxcoordglbwav,gycoordglbwav,maskdat=maskglbwav,&
        & extrapol=iopt_waves_extrapol.EQ.1,land=.TRUE.)

!
!3.3.2 Allocate
!--------------
!
!  ---number of non-zero elements
   noelements = 4*COUNT(maskglbwav.AND.WavtoCohCoordsglb%icoord.NE.int_fill)
!  ---allocate matrix arrays
   IF (mastermod) THEN
      ALLOCATE (rows(noelements),STAT=errstat)
      CALL error_alloc('rows',1,(/noelements/),kndint)
      rows = 0
      ALLOCATE (columns(noelements),STAT=errstat)
      CALL error_alloc('columns',1,(/noelements/),kndint)
      columns = 0
      ALLOCATE (weights(noelements),STAT=errstat)
      CALL error_alloc('weights',1,(/noelements/),kndrtype)
      weights = 0.0
   ENDIF

!
!3.3.3 Define matrix rows, columns and weigths
!---------------------------------------------
!

   IF (mastermod) THEN

      n = 1
      j_333: DO j=1,nrwav
      i_333: DO i=1,ncwav
         icoord = WavtoCohCoordsglb(i,j)%icoord
         IF (maskglbwav(i,j).AND.icoord.NE.int_fill) THEN
            jcoord = WavtoCohCoordsglb(i,j)%jcoord
            whts = WavtoCohCoordsglb(i,j)%weights
            idat = (/icoord,icoord+1,icoord+1,icoord/)
            jdat = (/jcoord,jcoord,jcoord+1,jcoord+1/)
            rows(n:n+3) = i + ncwav*(j-1)
            columns(n:n+3) = idat + nc*(jdat-1)
            weights(n:n+3) = (/whts(1,1),whts(2,1),whts(2,2),whts(1,2)/)
            n = n + 4
         ENDIF
      ENDDO i_333
      ENDDO j_333

!
!3.3.4 Construct sparse matrix
!-----------------------------
!
!     ---define
      nocohglb = nc*nr
      CALL SparseMatrix_init(Sm_CohinttoWav,nowavglb,nocohglb,noelements)
!     ---rows
      CALL SparseMatrix_importGRowInd(Sm_CohinttoWav,rows,noelements)
!     ---columns
      CALL SparseMatrix_importGColInd(Sm_CohinttoWav,columns,noelements)
!     ---insert values
      CALL SparseMatrix_importMatrixElts(Sm_CohinttoWav,weights,noelements)
 
   ENDIF

!
!3.5 Router and Matrix+ definitions
!----------------------------------
!

   CALL Router_init(modelidcoh,Gsm_WavRecvCoh,comm_model,Rout_CohtoWav)
    
   CALL SparseMatrixPlus_init(Smp_CohinttoWav,Sm_CohinttoWav,Gsm_WavRecvCoh,&
                            & Gsm_WavRecvWav,Yonly,idmaster,comm_model,modelid)

!
!3.6 Deallocate
!--------------
!

   DEALLOCATE (lengthw,startw)
   DEALLOCATE (WavtoCohCoordsglb)
   IF (mastermod) DEALLOCATE (rows,columns,weights)
   
ENDIF

!
!4. Sending data from the wave model to COHERENS
!-----------------------------------------------
!

IF (iopt_waves_couple.EQ.1.OR.iopt_waves_couple.EQ.3) THEN

!
!4.1 Names of variables in MCT list format
!-----------------------------------------
!

   varlist_wav = ''
   ivar_410: DO ivar=1,MaxWavtoCohvars
      IF (TRIM(wavtocoh_names(ivar)).NE.'') THEN
         varlist_wav = TRIM(varlist_wav)//TRIM(wavtocoh_names(ivar))//':'
      ENDIF
   ENDDO ivar_410
   l = LEN_TRIM(varlist_wav)
   varlist_wav(l:l) = ' '

!
!4.2 Global segmentation map
!---------------------------
!
!  ---allocate
   ALLOCATE (lengthw(nrwavloc),STAT=errstat)
   CALL error_alloc('lengthw',1,(/nrwavloc/),kndint)
   ALLOCATE (startw(nrwavloc),STAT=errstat)
   CALL error_alloc('startw',1,(/nrwavloc/),kndint)

!  ---start and length of the segments
   lengthw = ncwavloc
   j_420: DO j=1,nrwavloc
      startw(j) = nc1wavloc + (nr1wavloc-1)*ncwav + (j-1)*ncwav
   ENDDO j_420
   
!  ---define segmentation map
   CALL GlobalSegMap_init(Gsm_WavSendCoh,startw,lengthw,idmastermod,comm_model,&
                        & modelid)

!
!4.3 Attribute vector
!--------------------
!

   nowavloc = ncwavloc*nrwavloc
   CALL AttrVect_init(Av_WavSendCoh,rList=TRIM(varlist_wav),lsize=nowavloc)

!
!4.4 Routers and Matrix+ definitions
!----------------------------------
!
!  ---receiving data from COHERENS
   CALL Router_init(modelidcoh,Gsm_WavSendCoh,comm_model,Rout_WavtoCoh)

!
!4.5 Deallocate
!--------------
!

   DEALLOCATE (lengthw,startw)
   
ENDIF

!
!5. Restart monitoring
!---------------------
!

CALL log_timer_out(npcc,itm_MCT)

nt = 1
CALL monitor_files

#endif /*MCT*/

RETURN

END SUBROUTINE wav2coh_coupler

!========================================================================

SUBROUTINE wav2coh_end
!************************************************************************
!
! *wav2coh_end* Finalise wave simulation 
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.9
!
! Description - routine is called by all wave processes
!
! Reference -
!
! Calling program - wave model
!
! External calls - 
!
! Module calls - error_abort, write_error_message
!
!************************************************************************
!
USE couplingvars
USE iopars
USE switches
USE error_routines, ONLY: error_abort, write_error_message
USE time_routines, ONLY: log_timer_in, log_timer_out

#ifdef MCT
   USE m_AttrVect, ONLY: AttrVect_clean => clean
   USE m_GlobalSegMap, ONLY: GlobalSegMap_clean  => clean
   USE m_MCTWorld, ONLY: MCTWorld_clean => clean
   USE m_Router, ONLY: Router_clean => clean
   USE m_SparseMatrix, ONLY: SparseMatrix_clean => clean
   USE m_SparseMatrixPlus, ONLY: SparseMatrixPlus_clean => clean
#endif /*MCT*/

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc


#ifdef MCT

procname(pglev+1) = 'wav2coh_end'
CALL log_timer_in(npcc)

!
!1. Clean MCT types defined by the wave model
!--------------------------------------------
!
!---attribute vectors
IF (iopt_waves_couple.EQ.1.OR.iopt_waves_couple.EQ.3) THEN
   CALL AttrVect_clean(Av_WavSendCoh,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT attribute vector Av_WavSendCoh')
   ENDIF
ENDIF
IF (iopt_waves_couple.GT.1) THEN
   CALL AttrVect_clean(Av_WavRecvCoh,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT attribute vector Av_WavRecvCoh')
   ENDIF
   CALL AttrVect_clean(Av_WavRecvWav,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT attribute vector Av_WavRecvWav')
   ENDIF
ENDIF

!---global segmentation maps
IF (iopt_waves_couple.EQ.1.OR.iopt_waves_couple.EQ.3) THEN
   CALL GlobalSegMap_clean(Gsm_WavSendCoh,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT global segmentation map Gsm_WavSendCoh')
   ENDIF
ENDIF
IF (iopt_waves_couple.GT.1) THEN
   CALL GlobalSegMap_clean(Gsm_WavRecvCoh,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT global segmentation map Gsm_WavRecvCoh')
   ENDIF
   CALL GlobalSegMap_clean(Gsm_WavRecvWav,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT global segmentation map Gsm_WavRecvWav')
   ENDIF
ENDIF

!---routers
IF (iopt_waves_couple.EQ.1.OR.iopt_waves_couple.EQ.3) THEN
   CALL Router_clean(Rout_WavtoCoh,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message('Unable to clean MCT router Rout_WavtoCoh')
   ENDIF
ENDIF
IF (iopt_waves_couple.GT.1) THEN
   CALL Router_clean(Rout_CohtoWav,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message('Unable to clean MCT router Rout_CohtoWav')
   ENDIF
ENDIF

!---sparse matrix
IF (iopt_waves_couple.GT.1) THEN
   CALL SparseMatrix_clean(Sm_CohinttoWav,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT sparse matrix Sm_CohinttoWav')
   ENDIF
ENDIF

!---sparse matrix plus
IF (iopt_waves_couple.GT.1) THEN
   CALL SparseMatrixPlus_clean(Smp_CohinttoWav,errstat)
   IF (errstat.NE.0) THEN
      CALL write_error_message(&
                 & 'Unable to clean MCT sparse matrix plus Smp_CohinttoWav')
   ENDIF
ENDIF

!---MCT world
CALL MCTWorld_clean()

CALL error_abort('wav2coh_end',ierrno_MCT)

CALL log_timer_out(npcc,itm_MCT)

#endif /*MCT*/


RETURN

END SUBROUTINE wav2coh_end

!========================================================================

SUBROUTINE wav2coh_grid(gxglbwav,gyglbwav,maskglb)
!************************************************************************
!
! *wav2coh_grid* send wave grid to the COHERENS master, receive COHERENS grid
!                 from COHERENS master
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.9
!
! Description - routine is called by all wave processes
!
! Reference -
!
! Calling program - wave model
!
! External calls - 
!
! Module calls - comms_recv_log, comms_recv_real, comms_send_log,
!                comms_send_real, copy_vars, error_alloc
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
USE wavepars
USE wavevars
USE comms_MPI, ONLY: comms_recv_log, comms_recv_real, comms_send_log, &
     & comms_send_real
USE paral_comms, ONLY: copy_vars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
LOGICAL, INTENT(IN), DIMENSION(ncwav,nrwav) :: maskglb
REAL, INTENT(IN), DIMENSION(ncwav,nrwav) :: gxglbwav, gyglbwav

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*gxglbwav*  REAL     X-coordinates of the wave grid    [m or degrees longitude]
!*gyglbwav*  REAL     Y-coordinates of the wave grid    [m or degrees latitude]
!*maskglb*   LOGICAL  Land mask for the wave grid
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: itag, novals, npcc


procname(pglev+1) = 'wav2coh_grid'
CALL log_timer_in(npcc)

!
!1. Allocate
!-----------
!
!---COHERENS grid
ALLOCATE (gxcoordglb(0:nc+1,0:nr+1),STAT=errstat)
CALL error_alloc('gxcoordglb',2,(/nc+2,nr+2/),kndrtype)
ALLOCATE (gycoordglb(0:nc+1,0:nr+1),STAT=errstat)
CALL error_alloc('gycoordglb',2,(/nc+2,nr+2/),kndrtype)
ALLOCATE (seapointglb(0:nc+1,0:nr+1),STAT=errstat)
CALL error_alloc('seapointglb',2,(/nc+2,nr+2/),kndlog)

!---wave grid
ALLOCATE (gxcoordglbwav(ncwav,nrwav),STAT=errstat)
CALL error_alloc('gxcoordglbwav',2,(/ncwav,nrwav/),kndrtype)
ALLOCATE (gycoordglbwav(ncwav,nrwav),STAT=errstat)
CALL error_alloc('gycoordglbwav',2,(/ncwav,nrwav/),kndrtype)
ALLOCATE (maskglbwav(ncwav,nrwav),STAT=errstat)
CALL error_alloc('maskglbwav',2,(/ncwav,nrwav/),kndlog)

!
!2. Send wave grid to COHERENS master
!------------------------------------
!

IF (mastermod) THEN

   novals = ncwav*nrwav
   gxcoordglbwav = gxglbwav
   itag = 1000
   CALL comms_send_real(gxcoordglbwav,novals,idmaster,itag,comm_world_MPI,1,&
                      & iarr_gxcoordglbwav)
   gycoordglbwav = gyglbwav
   itag = itag + 1
   CALL comms_send_real(gycoordglbwav,novals,idmaster,itag,comm_world_MPI,1,&
                      & iarr_gycoordglbwav)
   maskglbwav = maskglb
   itag = itag + 1
   CALL comms_send_log(maskglbwav,novals,idmaster,itag,comm_world_MPI,1,&
                     & iarr_maskglbwav)

!
!3. Receive COHERENS grid from COHERENS master
!---------------------------------------------
!

   novals = (nc+2)*(nr+2)
   itag = itag + 1
   CALL comms_recv_real(gxcoordglb,novals,idmaster,itag,comm_world_MPI,&
                      & iarr_gxcoordglb)
   itag = itag + 1
   CALL comms_recv_real(gycoordglb,novals,idmaster,itag,comm_world_MPI,&
                      & iarr_gycoordglb)
   itag = itag + 1
   CALL comms_recv_log(seapointglb,novals,idmaster,itag,comm_world_MPI,&
                      & iarr_seapointglb)

ENDIF

!
!4. Broadcast within the wave communicator
!-----------------------------------------
!

IF (iopt_MPI.EQ.1) THEN

   CALL copy_vars(gxcoordglbwav,iarr_gxcoordglbwav,idroot=idmastermod,&
                & comm=comm_model,commtype=5)
   CALL copy_vars(gycoordglbwav,iarr_gycoordglbwav,idroot=idmastermod,&
                & comm=comm_model,commtype=5)
   CALL copy_vars(maskglbwav,iarr_maskglbwav,idroot=idmastermod,&
                & comm=comm_model,commtype=5)
   CALL copy_vars(gxcoordglb,iarr_gxcoordglb,idroot=idmastermod,&
                & comm=comm_model,commtype=5)
   CALL copy_vars(gycoordglb,iarr_gycoordglb,idroot=idmastermod,&
                & comm=comm_model,commtype=5)
   CALL copy_vars(seapointglb,iarr_seapointglb,idroot=idmastermod,&
                & comm=comm_model,commtype=5)

ENDIF

CALL log_timer_out(npcc,itm_MCT)


RETURN

END SUBROUTINE wav2coh_grid

!========================================================================

SUBROUTINE wav2coh_params(ncw,nrw,nhtype,coordtype,intdate1,intdate2,timestep,&
                        & realpars)
!************************************************************************
!
! *wav2coh_params* wave model parameters send to COHERENS
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.9
!
! Description - called only by the master process within the wave communicator 
!
! Reference - routine is called by all wave processes
!           - all arguments (except the first two) only need to be defined at
!             the wave master process
!
! Calling program - wave model
!
! External calls - 
!
! Module calls - comms_recv_int, comms_recv_log, comms_recv_real,
!                comms_send_int, comms_send_real, convert_date, copy_vars,
!                error_alloc, num_time_steps
!
!************************************************************************
!
USE gridpars
USE iopars
USE paralpars
USE switches
USE timepars
USE wavepars
USE comms_MPI, ONLY: comms_recv_int, comms_recv_log, &
                   & comms_recv_real, comms_send_int, comms_send_real
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: copy_vars 
USE time_routines, ONLY: convert_date, log_timer_in, log_timer_out, &
                       & num_time_steps

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: coordtype, ncw, nhtype, nrw
INTEGER, INTENT(IN), DIMENSION(6) :: intdate1, intdate2
REAL, INTENT(IN) :: timestep
REAL, INTENT(IN), DIMENSION(4) :: realpars

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*ncw*       INTEGER  X-dimension of the global wave grid
!*nrw*       INTEGER  Y-dimension of the global wave grid
!*nhtype*    INTEGER  Type of horizontail grid of the wave model
!                = 1 => uniform rectangular
!                = 2 => non-uniform rectangular
!                = 3 => irregular (non-structured)
!*coordtype* INTEGER  Type of coordinates used by the wave model
!                = 0 => Cartesian
!                = 1 => spherical
!*intdate1*  INTEGER  Start date/time for wave/current coupling
!                     (year,month,day,hour,min,secs)
!*intdate2*  INTEGER  End date/time for wave/current coupling
!                     (year,month,day,hour,min,secs)
!*realpars*  REAL     Coordinates of the lower left corner and grid spacings of
!                     the wave grid (used only in case of a uniform rectangular
!                     grid)
!*gxglbwav*  REAL     X-coordinates of the wave grid    [m or degrees longitude]
!*gyglbwav*  REAL     Y-coordinates of the wave grid     [m or degrees latitude]
!*maskglb*   LOGICAL  Land mask for the wave grid
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: iproc, itag, npcc
LOGICAL, DIMENSION(1) :: logrcv
INTEGER, DIMENSION(3) :: intrcv 
INTEGER, DIMENSION(4) :: intsnd
REAL, DIMENSION(8) :: realrcv
REAL, DIMENSION(5) :: realsnd


pglev = 2
procname(pglev+1) = 'wav2coh_params'
CALL log_timer_in(npcc)

!
!1. Parallel settings
!--------------------
!
!---parameters for combine and collect
ALLOCATE(comprocs(2*nprocs),STAT=errstat)
CALL error_alloc('comprocs',1,(/2*nprocs/),kndint)
comprocs = (/0,(1,iproc=1,nprocs-1),-1,(2,iproc=1,nprocs-1)/)
comprocs = CSHIFT(comprocs,SHIFT=-idloc)

!---MPI switch
iopt_MPI = MERGE(0,1,nprocswav.EQ.1)

!---error checking
errchk = mastermod

!
!2. Wave model parameters
!------------------
!
!---grid sizes
ncwav = ncw; nrwav = nrw

!---surface grid parameters
surfacegrids(igrd_waves,1)%n1dat = ncwav
surfacegrids(igrd_waves,1)%n2dat = nrwav
surfacegrids(igrd_waves,1)%nhtype = nhtype
surfacegrids(igrd_waves,1)%x0dat = realpars(1)
surfacegrids(igrd_waves,1)%y0dat = realpars(2)
surfacegrids(igrd_waves,1)%delxdat = realpars(3)
surfacegrids(igrd_waves,1)%delydat = realpars(4)

!---type of coordinates
iopt_grid_sph = coordtype

!---start/end dates
IStartDateTimeWav(1:6) = intdate1; IStartDateTimeWav(7) = 0
CStartDateTimeWav = convert_date(IStartDateTimeWav)
IEndDateTimeWav(1:6) = intdate2; IEndDateTimeWav(7) = 0
CEndDateTimeWav = convert_date(IEndDateTimeWav)

!---time step
deltwav = timestep
CALL num_time_steps(CStartDateTimeWav,CEndDateTimeWav,deltwav,nstep)
runlog_count = nstep

!
!3. Send wave model data to the COHERENS master process
!------------------------------------------------------
!

IF (mastermod) THEN

!  ---start/end date
   itag = 1000
   CALL comms_send_int(intdate1,6,idmaster,itag,comm_world_MPI,1,0)
   itag = itag + 1
   CALL comms_send_int(intdate2,6,idmaster,itag,comm_world_MPI,1,0)

!  ---integer parameters
   intsnd = (/ncw,nrw,nhtype,coordtype/)
   itag = itag + 1
   CALL comms_send_int(intsnd,4,idmaster,itag,comm_world_MPI,1,0)

!  ---real parameters
   realsnd(1) = timestep; realsnd(2:5) = realpars
   itag = itag + 1
   CALL comms_send_real(realsnd,5,idmaster,itag,comm_world_MPI,1,0)

!
!4. Receive COHERENS grid parameters from COHERENS master
!--------------------------------------------------------
!

   itag = itag + 1
   CALL comms_recv_int(intrcv,3,idmaster,itag,comm_world_MPI,0)
   itag = itag + 1
   CALL comms_recv_log(logrcv,1,idmaster,itag,comm_world_MPI,0)
   itag = itag + 1
   CALL comms_recv_real(realrcv,8,idmaster,itag,comm_world_MPI,0)

ENDIF

!
!5. Broadcast within the wave communicator
!-----------------------------------------
!

IF (nprocswav.GT.1) THEN
   CALL copy_vars(intrcv,0,idroot=idmastermod,comm=comm_model,commtype=5)
   CALL copy_vars(logrcv,0,idroot=idmastermod,comm=comm_model,commtype=5)
   CALL copy_vars(realrcv,0,idroot=idmastermod,comm=comm_model,commtype=5)
ENDIF

!
!6. COHERENS grid parameters
!---------------------------
!

nc = intrcv(1); nr = intrcv(2)
surfacegrids(igrd_model,1)%n1dat = nc; surfacegrids(igrd_model,1)%n2dat = nr;
iopt_grid_htype = intrcv(3)
surfacegrids(igrd_model,1)%nhtype = iopt_grid_htype
surfacegrids(igrd_model,1)%rotated = logrcv(1)
surfacegrids(igrd_model,1)%x0dat = realrcv(1)
surfacegrids(igrd_model,1)%y0dat = realrcv(2)
surfacegrids(igrd_model,1)%delxdat = realrcv(3)
surfacegrids(igrd_model,1)%delydat = realrcv(4)
surfacegrids(igrd_model,1)%gridangle = realrcv(5)
surfacegrids(igrd_model,1)%longpole = realrcv(6)
surfacegrids(igrd_model,1)%x0rot = realrcv(7)
surfacegrids(igrd_model,1)%y0rot = realrcv(8)

CALL log_timer_out(npcc,itm_MCT)


RETURN

END SUBROUTINE wav2coh_params

!========================================================================

SUBROUTINE wav2coh_recv(depmeanwav,zetawav,umvelwav,vmvelwav)
!************************************************************************
!
! *wav2coh_recv* receive data from Coherens for the wave model
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.9
!
! Description - routine is called by all wave processes
!
! Reference -
!
! Calling program - wave model
!
! External calls - 
!
! Module calls - 
!
!************************************************************************
!
USE couplingvars
USE iopars
USE paralpars
USE switches
USE wavepars
USE time_routines, ONLY: log_timer_in, log_timer_out

#ifdef MCT
   USE m_AttrVect,ONLY : AttrVect_indexR  => indexRA
   USE m_Transfer, ONLY : MCT_Recv => recv
   USE m_MatAttrVectMul, ONLY: MCT_MatVecMul => sMatAvMult
#endif /*MCT*/

IMPLICIT NONE

!
!*Arguments
!
REAL, INTENT(OUT), DIMENSION(ncwavloc_ext,nrwavloc_ext) :: depmeanwav, &
                                                   & umvelwav, vmvelwav, zetawav
!
! Name        Type  Purpose
!------------------------------------------------------------------------------
!*depmeanwav* REAL  COHERENS bathymetry interpolated on the wave grid [m]
!*zetawav*    REAL  COHERENS elevations interpolated on the wave grid [m]
!*umvelwav*   REAL  COHERENS X-component of the depth mean current interpolated
!                   on the wave grid [m/s]
!*vmvelwav*   REAL  COHERENS Y-component of the depth mean current interpolated
!                   on the wave grid [m/s]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: iivar, ivar, npcc
REAL, DIMENSION(ncwavloc_ext,nrwavloc_ext) :: rmask


IF (iopt_waves_couple.EQ.1) RETURN

#ifdef MCT

procname(pglev+1) = 'wav2coh_recv'
CALL log_timer_in(npcc)

!
!1. Receive data on the COHERENS grid
!------------------------------------
!

CALL MCT_Recv(Av_WavRecvCoh,Rout_CohtoWav,Tag=itagcohtowav)

!
!2. Interpolate onto the wave grid
!---------------------------------
!

CALL MCT_MatVecMul(Av_WavRecvCoh,Smp_CohinttoWav,Av_WavRecvWav)

!
!3. Retrieve data from MCT attribute vector
!------------------------------------------
!
!3.1 Land mask
!-------------
!

ivar = AttrVect_indexR(Av_WavRecvWav,'rmask')
rmask = RESHAPE(Av_WavRecvWav%rATTr(ivar,:),(/ncwavloc_ext,nrwavloc_ext/))
WHERE (rmask.GT.0.0)
   rmask = 1.0/rmask
ELSEWHERE
   rmask = 0.0
END WHERE

!
!3.2 Other data
!--------------
!

ivar_320: DO ivar=1,MaxCohtoWavvars

!  ---array index within MCT attribute vector
   SELECT CASE (TRIM(cohtowav_names(ivar)))
      CASE ('depmeanatc')
         iivar = AttrVect_indexR(Av_WavRecvWav,'depmeanatc')
         depmeanwav = rmask*RESHAPE(Av_WavRecvWav%rATTr(iivar,:),&
                                 & (/ncwavloc_ext,nrwavloc_ext/))
      CASE ('zeta')
         iivar = AttrVect_indexR(Av_WavRecvWav,'zeta')
         zetawav = rmask*RESHAPE(Av_WavRecvWav%rATTr(iivar,:),&
                              & (/ncwavloc_ext,nrwavloc_ext/))
      CASE ('umvel')
         iivar = AttrVect_indexR(Av_WavRecvWav,'umvel')
         umvelwav = rmask*RESHAPE(Av_WavRecvWav%rATTr(iivar,:),&
                               & (/ncwavloc_ext,nrwavloc_ext/))
      CASE ('vmvel')
         iivar = AttrVect_indexR(Av_WavRecvWav,'vmvel')
         vmvelwav = rmask*RESHAPE(Av_WavRecvWav%rATTr(iivar,:),&
                               & (/ncwavloc_ext,nrwavloc_ext/))
   END SELECT

ENDDO ivar_320

CALL log_timer_out(npcc,itm_MCT)

#endif /*MCT*/


RETURN

END SUBROUTINE wav2coh_recv

!========================================================================

SUBROUTINE wav2coh_send(rmaskw,umstokesw,vmstokesw,umswdissipw,vmswdissipw,&
                      & umbwdissipw,vmbwdissipw,wpres,wheight,wperiod,&
                      & wdir,wcur,wexcurs)
!************************************************************************
!
! *wav2coh_send* send data from the wave model to Coherens
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.9
!
! Description - routine is called by all wave processes
!
! Reference -
!
! Calling program - wave model
!
! External calls - 
!
! Module calls - 
!
!************************************************************************
!
USE couplingvars
USE iopars
USE paralpars
USE switches
USE wavepars
USE time_routines, ONLY: log_timer_in, log_timer_out

#ifdef MCT
   USE m_Transfer,ONLY : MCT_Send => send
   USE m_AttrVect,ONLY : AttrVect_indexR  => indexRA
#endif /*MCT*/

IMPLICIT NONE

!
!*Arguments
!
REAL, INTENT(IN), DIMENSION(ncwavloc_ext,nrwavloc_ext) :: rmaskw, umbwdissipw, &
    & umstokesw, umswdissipw, vmbwdissipw, vmstokesw, vmswdissipw, wcur, wdir, &
    & wexcurs, wheight, wperiod, wpres

! Name          Type  Purpose
!------------------------------------------------------------------------------
!*rmaskw*       REAL  Local mask array set to 1 for land points on the wave
!                     model grid
!*umstokesw*    REAL  X-component of the depth-mean Stokes drift           [m/s]
!*vmstokesw*    REAL  Y-component of the depth-mean Stokes drift           [m/s]
!*umswdissipw*  REAL  X-component of the surface wave dissipation force   [m/s2]
!*vmswdissipw*  REAL  Y-component of the surface wave dissipation force   [m/s2]
!*umbwdissipw*  REAL  X-component of the bed wave dissipation force       [m/s2]
!*vmbwdissipw*  REAL  Y-component of the bed wave dissipation force       [m/s2]
!*wpres*        REAL  Wave-induced pressure                                 [Pa]
!*wheight*      REAL  Significant wave height                                [m]
!*wperiod*      REAL  Peak wave period                                       [s]
!*wdir*         REAL  Mean wave direction                                  [deg]
!*wcur*         REAL  Near-bottom wave orbital velocity                    [m/s]
!*wexcurs*      REAL  Near-bottom wave excursion amplitude                   [m]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, iivar, ivar, j, n, npcc


IF (iopt_waves_couple.EQ.2) RETURN

#ifdef MCT

procname(pglev+1) = 'wav2coh_send'
CALL log_timer_in(npcc)

!
!1. Store values in the MCT attribute vector
!-------------------------------------------
!

ivar_110: DO ivar=1,MaxWavtoCohvars

!  ---array index within MCT attribute vector 
   SELECT CASE (TRIM(wavtocoh_names(ivar)))
      CASE ('rmask')
         iivar = AttrVect_indexR(Av_WavSendCoh,'rmask')
      CASE ('umstokesatc')
         iivar = AttrVect_indexR(Av_WavSendCoh,'umstokesatc')
      CASE ('vmstokesatc')
         iivar = AttrVect_indexR(Av_WavSendCoh,'vmstokesatc')
      CASE ('umswdissipatc')
         iivar = AttrVect_indexR(Av_WavSendCoh,'umswdissipatc')
      CASE ('vmswdissipatc')
         iivar = AttrVect_indexR(Av_WavSendCoh,'vmswdissipatc')
      CASE ('umbwdissipatc')
         iivar = AttrVect_indexR(Av_WavSendCoh,'umbwdissipatc')
      CASE ('vmbwdissipatc')
         iivar = AttrVect_indexR(Av_WavSendCoh,'vmbwdissipatc')
      CASE ('wavepres')
         iivar = AttrVect_indexR(Av_WavSendCoh,'wavepres')
      CASE ('waveheight')
         iivar = AttrVect_indexR(Av_WavSendCoh,'waveheight')
      CASE ('waveperiod')
         iivar = AttrVect_indexR(Av_WavSendCoh,'waveperiod')
      CASE ('wavedir')
         iivar = AttrVect_indexR(Av_WavSendCoh,'wavedir')
      CASE ('wavevel')
         iivar = AttrVect_indexR(Av_WavSendCoh,'wavevel')
      CASE ('waveexcurs')
         iivar = AttrVect_indexR(Av_WavSendCoh,'waveexcurs')
   END SELECT

!  ---store values
   n = 0
   j_111: DO j=1+haloS,nrwavloc+haloS
   i_111: DO i=1+haloW,ncwavloc+haloW
      n = n + 1
      SELECT CASE (TRIM(wavtocoh_names(ivar)))
         CASE ('rmask')
            Av_WavSendCoh%RattR(iivar,n) = rmaskw(i,j)
         CASE ('umstokesatc')
            Av_WavSendCoh%RattR(iivar,n) = rmaskw(i,j)*umstokesw(i,j)
         CASE ('vmstokesatc')
            Av_WavSendCoh%RattR(iivar,n) = rmaskw(i,j)*vmstokesw(i,j)
         CASE ('umswdissipatc')
            Av_WavSendCoh%RattR(iivar,n) = rmaskw(i,j)*umswdissipw(i,j)
         CASE ('vmswdissipatc')
            Av_WavSendCoh%RattR(iivar,n) = rmaskw(i,j)*vmswdissipw(i,j)
         CASE ('umbwdissipatc')
            Av_WavSendCoh%RattR(iivar,n) = rmaskw(i,j)*umbwdissipw(i,j)
         CASE ('vmbwdissipatc')
            Av_WavSendCoh%RattR(iivar,n) = rmaskw(i,j)*vmbwdissipw(i,j)
         CASE ('wavepres')
            Av_WavSendCoh%RattR(iivar,n) = rmaskw(i,j)*wpres(i,j)
         CASE ('waveheight')
            Av_WavSendCoh%RattR(iivar,n) = rmaskw(i,j)*wheight(i,j)
         CASE ('waveperiod')
            Av_WavSendCoh%RattR(iivar,n) = rmaskw(i,j)*wperiod(i,j)
         CASE ('wavedir')
            Av_WavSendCoh%RattR(iivar,n) = rmaskw(i,j)*wdir(i,j)
         CASE ('wavevel')
            Av_WavSendCoh%RattR(iivar,n) = rmaskw(i,j)*wcur(i,j)
         CASE ('waveexcurs')
            Av_WavSendCoh%RattR(iivar,n) = rmaskw(i,j)*wexcurs(i,j)
      END SELECT
   ENDDO i_111
   ENDDO j_111

ENDDO ivar_110

!
!2. Send data
!------------
!

CALL MCT_Send(Av_WavSendCoh,Rout_WavtoCoh,Tag=itagwavtocoh)

CALL log_timer_out(npcc,itm_MCT)

#endif /*MCT*/


RETURN

END SUBROUTINE wav2coh_send

!========================================================================

SUBROUTINE wav2coh_switches(logswitches)
!************************************************************************
!
! *wav2coh_switches* define which set of variables are used for coupling
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.11.1
!
! Description - 
!
! Reference - routine is called by all wave processes
!
! Calling program - wave model
!
! External calls - 
!
! Module calls - comms_recv_int, copy_vars
!
!************************************************************************
!
USE iopars
USE paralpars
USE switches
USE comms_MPI, ONLY: comms_recv_int
USE paral_comms, ONLY: copy_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
LOGICAL, INTENT(OUT), DIMENSION(8) :: logswitches

! Name            Type    Purpose
!------------------------------------------------------------------------------
!*logswitches*    LOGICAL switches to determine which variables are sent to or
!                         received from COHERENS
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: itag, npcc
INTEGER, DIMENSION(6) :: intdat


procname(pglev+1) = 'wav2coh_switches'
CALL log_timer_in(npcc)

!
!1. Receive from COHERENS master
!-------------------------------
!

IF (mastermod) THEN
   itag = 1000
   CALL comms_recv_int(intdat,6,idmaster,itag,comm_world_MPI,0)
ENDIF

!
!2. Broadcast within the wave communicator
!-----------------------------------------
!

IF (nprocswav.GT.1) THEN
   CALL copy_vars(intdat,0,idroot=idmastermod,comm=comm_model,commtype=5)
ENDIF

!
!3. Model switches
!-----------------
!

iopt_waves_couple = intdat(1)
iopt_bstres_waves_bfric = intdat(2)
iopt_waves_curr = intdat(3)
iopt_waves_dissip = intdat(4)
iopt_waves_form = intdat(5)
iopt_waves_pres = intdat(6)

!
!4. Logical switches for variable selection
!------------------------------------------ 
!

logswitches(1:4) = .TRUE.
logswitches(5) = MERGE(.TRUE.,.FALSE.,iopt_waves_form.EQ.2)
logswitches(6) = MERGE(.TRUE.,.FALSE.,&
                     & iopt_waves_curr.EQ.1.AND.iopt_waves_form.EQ.2)
logswitches(7) = MERGE(.TRUE.,.FALSE.,&
                     & iopt_waves_pres.EQ.1.AND.iopt_waves_form.EQ.2)
logswitches(8) = MERGE(.TRUE.,.FALSE.,iopt_waves_dissip.EQ.1)

CALL log_timer_out(npcc,itm_MCT)


RETURN

END SUBROUTINE wav2coh_switches

!========================================================================

SUBROUTINE wave_model
!************************************************************************
!
! *wave_model* run external wave model
!
! Author - Alexander Breugem and Patrick Luyten (IMDC)
!
! Version - @(COHERENS)Wave_Coupler.F90  V2.9
!
! Description -
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - swan
!
! Module calls - error_abort, error_ubound_var
!
!************************************************************************
!
USE iopars
USE paralpars
USE switches
USE error_routines, ONLY: error_abort, error_ubound_var
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
INTEGER :: npcc


procname(pglev+1) = 'wave_model'
CALL log_timer_in(npcc)

IF (cold_start) RETURN

IF (iopt_waves_model.EQ.1) THEN

# ifdef SWAN
  IF (nosimul.GT.1) THEN
     CALL error_ubound_var(nosimul,'nosimul',1,.TRUE.)
     WRITE (ioerr,'(A)') 'Only one simulation is allowed with the SWAN model'
     CALL error_abort('wave_model',ierrno_waves)
  ENDIF
  CALL swan(comm_model)
#endif /*SWAN*/

ENDIF

CALL log_timer_out(npcc,itm_wavemod)


RETURN

END SUBROUTINE wave_model
