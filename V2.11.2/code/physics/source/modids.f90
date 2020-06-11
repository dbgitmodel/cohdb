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

MODULE modids
!************************************************************************
!
! *modids* Model variable key ids
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modids.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description -
!
! Reference -
!
!************************************************************************
!

USE syspars

IMPLICIT NONE

!---model grid
INTEGER, PARAMETER :: &
& iarr_alphatc_fld = 1, iarr_alphatu_fld = 2, iarr_alphatv_fld = 3, &
& iarr_coriolatu = 4, iarr_coriolatv = 5, iarr_delxatc = 6, iarr_delxatu = 7, &
& iarr_delxatuv = 8, iarr_delxatv = 9, iarr_delyatc = 10, iarr_delyatu = 11, &
& iarr_delyatuv = 12, iarr_delyatv = 13, iarr_delzatc = 14, iarr_delzatu = 15, &
& iarr_delzatuv = 16, iarr_delzatuw = 17, iarr_delzatv = 18, &
& iarr_delzatvw = 19, iarr_delzatw = 20, iarr_dryfac = 21, iarr_gaccatc = 22, &
& iarr_gaccatu = 23, iarr_gaccatv = 24, iarr_gangleatc = 25, iarr_garea = 26, &
& iarr_gdelxglb = 27, iarr_gdelyglb = 28, iarr_gscoordatc = 29, &
& iarr_gscoordatu = 30, iarr_gscoordatuvw = 31, iarr_gscoordatuw = 32, &
& iarr_gscoordatv = 33, iarr_gscoordatvw = 34, iarr_gscoordatw = 35, &
& iarr_gscoordglb = 36, iarr_gsigcoordatc = 37, iarr_gsigcoordatw = 38, &
& iarr_gxcoord = 39, iarr_gxcoordglb = 40, iarr_gxlon = 41, iarr_gycoord = 42, &
& iarr_gycoordglb = 43, iarr_gylat = 44, iarr_indexobu = 45, &
& iarr_indexobuprocs = 46, iarr_indexobv = 47, iarr_indexobvprocs = 48, &
& iarr_indexobx = 49, iarr_indexobxprocs = 50, iarr_indexoby = 51, &
& iarr_indexobyprocs = 52, iarr_iobu = 53, iarr_iobuloc = 54, iarr_iobv = 55, &
& iarr_iobvloc = 56, iarr_iobx = 57, iarr_iobxloc = 58, iarr_ioby = 59, &
& iarr_iobyloc = 60, iarr_jobu = 61, iarr_jobuloc = 62, iarr_jobv = 63, &
& iarr_jobvloc = 64, iarr_jobx = 65, iarr_jobxloc = 66, iarr_joby = 67, &
& iarr_jobyloc = 68, iarr_maskatc_int = 69, iarr_mg_nc1procs = 70, &
& iarr_mg_nc2procs = 71, iarr_mg_nr1procs = 72, iarr_mg_nr2procs = 73, &
& iarr_ncprocs = 74, iarr_nobuprocs = 75, iarr_nobvprocs = 76, &
& iarr_nobxprocs = 77, iarr_nobyprocs = 78, iarr_nodeatc = 79, &
& iarr_nodeatu = 80, iarr_nodeatuv = 81, iarr_nodeatuw = 82, &
& iarr_nodeatv = 83, iarr_nodeatvw = 84, iarr_node2du = 85, &
& iarr_node2duv = 86, iarr_node2dv = 87, iarr_nosbuprocs = 88, &
& iarr_nosbvprocs = 89, iarr_nrprocs = 90, iarr_nrvbuprocs = 91, &
& iarr_nrvbvprocs = 92, iarr_rlxobcatu = 93, iarr_rlxobcatv = 94, &
& iarr_seapoint = 95, iarr_seapointglb = 96, iarr_soutobv = 97, &
& iarr_soutoby = 98, iarr_westobu = 99, iarr_westobx = 100

!---depths
INTEGER, PARAMETER :: &
& iarr_depmeanatc = 101, iarr_depmeanatu = 102, iarr_depmeanatuv = 103, &
& iarr_depmeanatv = 104, iarr_depmeanglb = 105, iarr_deptotatc = 106, &
& iarr_deptotatc_err = 107, iarr_deptotatc_old = 108, iarr_deptotatu = 109, &
& iarr_deptotatu_old = 110, iarr_deptotatuv = 111, iarr_deptotatv = 112, &
& iarr_deptotatv_old = 113, iarr_dzeta = 114, iarr_zeta = 115, &
& iarr_zeta_old = 116

!---currents
INTEGER, PARAMETER :: &
& iarr_hdvelmag = 117, iarr_hmvelmag = 118, iarr_hmvelmag_hadv_cfl = 119, &
& iarr_hvelmag = 120, iarr_hvelmag_hadv_cfl = 121, iarr_p2dbcgradatu = 122, &
& iarr_p2dbcgradatv = 123, iarr_p3dbcgradatu = 124, iarr_p3dbcgradatv = 125, &
& iarr_udevint = 126, iarr_udfvel = 127, iarr_udvel = 128, &
& iarr_udvel_old = 129, iarr_ufvel = 130, iarr_umpred = 131, iarr_umvel = 132, &
& iarr_umvel_hadv_cfl = 133, iarr_umvel_old = 134, iarr_uvel = 135, &
& iarr_uvel_hadv_cfl = 136, iarr_uvel_old = 137, iarr_vdevint = 138, &
& iarr_vdfvel = 139, iarr_vdvel = 140, iarr_vdvel_old = 141, iarr_vel2d = 142, &
& iarr_vel3d = 143, iarr_vfvel = 144, iarr_vmpred = 145, iarr_vmvel = 146, &
& iarr_vmvel_hadv_cfl = 147, iarr_vmvel_old = 148, iarr_vvel = 149, &
& iarr_vvel_hadv_cfl = 150, iarr_vvel_old = 151, iarr_wfvel = 152, &
& iarr_wphys = 153, iarr_wvel = 154, iarr_wvel_vadv_cfl = 155

!---density
INTEGER, PARAMETER :: &
& iarr_beta_sal = 156, iarr_beta_temp = 157, iarr_dens = 158, iarr_sal = 159, &
& iarr_temp = 160

!---diffusion coefficients
INTEGER, PARAMETER :: &
& iarr_hdifcoef2datc = 161, iarr_hdifcoef2d_mom = 162, &
& iarr_hdifcoef2d_scal = 163, iarr_hdifcoef2datuv = 164, &
& iarr_hdifcoef3datc = 165, iarr_hdifcoef3d_mom = 166, &
& iarr_hdifcoef3d_scal = 167, iarr_hdifcoef3datu = 168, &
& iarr_hdifcoef3datuv = 169, iarr_hdifcoef3datv = 170, &
& iarr_hdifcoef3datw = 171, iarr_kinvisc = 172, iarr_mom_vdif_cfl = 173, &
& iarr_scal_vdif_cfl = 174, iarr_xslopeatu_geo = 175, &
& iarr_xslopeatu_ziso = 176, iarr_xslopeatw_geo = 177, &
& iarr_yslopeatv_geo = 178, iarr_yslopeatv_ziso = 179, &
& iarr_yslopeatw_geo = 180, iarr_vdifcoefmom = 181, iarr_vdifcoefscal = 182, &
& iarr_vdifcoefscal_norot = 183, iarr_vdifcoefscal_rot = 184, &
& iarr_vdifcoeftke = 185, iarr_2D_hdif_cfl = 186

!---turbulence
INTEGER, PARAMETER :: &
& iarr_buofreq2 = 187, iarr_dissip = 188, iarr_ricnum = 189, &
& iarr_shearfreq2 = 190, iarr_tke = 191, iarr_tke_old = 192, iarr_tkezl = 193, &
& iarr_zlmix = 194

!---tides
INTEGER, PARAMETER :: &
& iarr_astro_earth = 195, iarr_fnode_astro = 196, iarr_fnode_obc = 197, &
& iarr_fxastro = 198, iarr_fyastro = 199, iarr_index_astro = 200, &
& iarr_index_obc = 201, iarr_ispec_tides = 202, iarr_phase_astro = 203, &
& iarr_phase_obc = 204, iarr_tidal_spectrum = 205

!---meteo
INTEGER, PARAMETER :: &
& iarr_airtemp = 206, iarr_atmpres = 207, iarr_cloud_cover = 208, &
& iarr_evapminprec = 209, iarr_evaporation = 210, iarr_precipitation = 211, &
& iarr_relhum = 212, iarr_sst = 213, iarr_uwindatc = 214, &
& iarr_uwindatc_old = 215, iarr_vappres_air = 216, iarr_vwindatc = 217, &
& iarr_vwindatc_old = 218, iarr_windatc = 219

!---waves
INTEGER, PARAMETER :: &
& iarr_gxcoordglbwav = 220, iarr_gycoordglbwav = 221, iarr_hbwdissipmag = 222, &
& iarr_hmbwdissipmag = 223, iarr_hmswdissipmag = 224, iarr_hswdissipmag = 225, &
& iarr_maskglbwav = 226, iarr_ubwdissipatc = 227, iarr_umbwdissipatc = 228, &
& iarr_umswdissipatc = 229, iarr_uswdissipatc = 230, iarr_vbwdissipatc = 231, &
& iarr_vmbwdissipatc = 232, iarr_vmswdissipatc = 233, iarr_vswdissipatc = 234, &
& iarr_wavedir = 235, iarr_waveexcurs = 236, iarr_wavefreq = 237, &
& iarr_waveheight = 238, iarr_wavenum = 239, iarr_waveperiod = 240, &
& iarr_wavepres = 241, iarr_wavevel = 242

!---Stokes velocities
INTEGER, PARAMETER :: &
& iarr_hmstokesmag = 243, iarr_hmveltotmag = 244, iarr_hstokesmag = 245, &
& iarr_hveltotmag = 246, iarr_stokessource2du = 247, &
& iarr_stokessource2dv = 248, iarr_udstokesatu = 249, iarr_umstokesatc = 250, &
& iarr_umstokesatu = 251, iarr_umveltot = 252, iarr_ustokesatc = 253, &
& iarr_ustokesatu = 254, iarr_uveltot = 255, iarr_vdstokesatv = 256, &
& iarr_vmstokesatc = 257, iarr_vmstokesatv = 258, iarr_vmveltot = 259, &
& iarr_vstokesatc = 260, iarr_vstokesatv = 261, iarr_vveltot = 262, &
& iarr_wstokesatw = 263

!---optical arrays
INTEGER, PARAMETER :: &
& iarr_optattcoef2 = 264, iarr_qrad = 265, iarr_radiance = 266

!---bottom/surface fluxes
INTEGER, PARAMETER :: &
& iarr_bdragcoefatc = 267, iarr_bfricatu = 268, iarr_bfricatv = 269, &
& iarr_bfricvel = 270, iarr_bfricvel_max = 271, iarr_bfricvel_wav = 272, &
& iarr_bstresatc = 273, iarr_bstresatc_max = 274, iarr_bstresatc_wav = 275, &
& iarr_bstresatu = 276, iarr_bstresatv = 277, iarr_cds = 278, iarr_ces = 279, &
& iarr_chs = 280, iarr_fwave = 281, iarr_qlatent = 282, iarr_qlwave = 283, &
& iarr_qnonsol = 284, iarr_qsensible = 285, iarr_qtot = 286, &
& iarr_sfricatc = 287, iarr_ssalflux = 288, iarr_sstresatc = 289, &
& iarr_ubstresatc = 290, iarr_ubstresatu = 291, iarr_usstresatc = 292, &
& iarr_usstresatu = 293, iarr_vbstresatc = 294, iarr_vbstresatv = 295, &
& iarr_vsstresatc = 296, iarr_vsstresatv = 297, iarr_wavethickatc = 298, &
& iarr_zaroughatc = 299, iarr_zroughatc = 300

!---open boundary forcing
INTEGER, PARAMETER :: &
& iarr_floutobu = 301, iarr_floutobv = 302, iarr_gxslope = 303, &
& iarr_gxslope_amp = 304, iarr_gxslope_pha = 305, iarr_gyslope = 306, &
& iarr_gyslope_amp = 307, iarr_gyslope_pha = 308, iarr_iloczobu = 309, &
& iarr_iloczobv = 310, iarr_indexprof = 311, iarr_indexvar = 312, &
& iarr_index2dobuv = 313, iarr_index2dobxy = 314, iarr_iobc2dtype = 315, &
& iarr_iprofobu = 316, iarr_iprofobv = 317, iarr_iprofobx = 318, &
& iarr_iprofoby = 319, iarr_iqsecobu = 320, iarr_isur1dtype = 321, &
& iarr_itypobu = 322, iarr_itypobv = 323, iarr_itypobx = 324, &
& iarr_itypoby = 325, iarr_ityp2dobu = 326, iarr_ityp2dobv = 327, &
& iarr_ityp2dobx = 328, iarr_ityp2doby = 329, iarr_jqsecobv = 330, &
& iarr_noprofsd = 331, iarr_no2dobuv = 332, iarr_no2dobxy = 333, &
& iarr_obcbioatu = 334, iarr_obcbioatv = 335, iarr_obcsalatu = 336, &
& iarr_obcsalatv = 337, iarr_obctmpatu = 338, iarr_obctmpatv = 339, &
& iarr_obc2uvatu = 340, iarr_obc2uvatu_old = 341, iarr_obc2uvatv = 342, &
& iarr_obc2uvatv_old = 343, iarr_obc3uvatu = 344, iarr_obc3uvatv = 345, &
& iarr_profbio = 346, iarr_profsal = 347, iarr_profsed = 348, &
& iarr_proftmp = 349, iarr_profvel = 350, iarr_return_time = 351, &
& iarr_udatobu = 352, iarr_udatobu_amp = 353, iarr_udatobu_pha = 354, &
& iarr_udatoby = 355, iarr_udatoby_amp = 356, iarr_udatoby_pha = 357, &
& iarr_vdatobv = 358, iarr_vdatobv_amp = 359, iarr_vdatobv_pha = 360, &
& iarr_vdatobx = 361, iarr_vdatobx_amp = 362, iarr_vdatobx_pha = 363, &
& iarr_vel2dobc = 364, iarr_zdatobu = 365, iarr_zdatobu_amp = 366, &
& iarr_zdatobu_pha = 367, iarr_zdatobv = 368, iarr_zdatobv_amp = 369, &
& iarr_zdatobv_pha = 370, iarr_zetaobc = 371, iarr_zetasur = 372, &
& iarr_zeta_amp = 373, iarr_zeta_pha = 374

!---structures
INTEGER, PARAMETER :: &
& iarr_idry = 375, iarr_indexwbaru = 376, iarr_indexwbaruprocs = 377, &
& iarr_indexwbarv = 378, iarr_indexwbarvprocs = 379, iarr_ithinu = 380, &
& iarr_ithinuloc = 381, iarr_ithinv = 382, iarr_ithinvloc = 383, &
& iarr_iwbaru = 384, iarr_iwbaruloc = 385, iarr_iwbarv = 386, &
& iarr_iwbarvloc = 387, iarr_jdry = 388, iarr_jthinu = 389, &
& iarr_jthinuloc = 390, iarr_jthinv = 391, iarr_jthinvloc = 392, &
& iarr_jwbaru = 393, iarr_jwbaruloc = 394, iarr_jwbarv = 395, &
& iarr_jwbarvloc = 396, iarr_nowbaruprocs = 397, iarr_nowbarvprocs = 398, &
& iarr_oricoefu = 399, iarr_oricoefv = 400, iarr_oriheightu = 401, &
& iarr_oriheightv = 402, iarr_orisillu = 403, iarr_orisillv = 404, &
& iarr_wbarcoefu = 405, iarr_wbarcoefv = 406, iarr_wbarcrestu = 407, &
& iarr_wbarcrestv = 408, iarr_wbarmodlu = 409, iarr_wbarmodlv = 410, &
& iarr_wbarelossu = 411, iarr_wbarelossv = 412

!---discharges
INTEGER, PARAMETER :: &
& iarr_disarea = 413, iarr_disdir = 414, iarr_disflag = 415, &
& iarr_dissal = 416, iarr_disspeed = 417, iarr_distmp = 418, &
& iarr_disvol = 419, iarr_idis = 420, iarr_idisloc = 421, &
& iarr_indexdisloc = 422, iarr_indexdisprocs = 423, iarr_jdis = 424, &
& iarr_jdisloc = 425, iarr_kdis = 426, iarr_kdistype = 427, &
& iarr_nodisprocs = 428, iarr_xdiscoord = 429, iarr_ydiscoord = 430, &
& iarr_zdiscoord = 431

!---parameters for parallel mode
INTEGER, PARAMETER :: &
& iarr_idprocs = 432

!---energy equation, enstrophy, vorticity
INTEGER, PARAMETER :: &
& iarr_edens0d = 433, iarr_edens2d = 434, iarr_edens3d = 435, &
& iarr_edissip0d = 436, iarr_edissip2d = 437, iarr_edissip3d = 438, &
& iarr_eflux2du = 439, iarr_eflux2dv = 440, iarr_eflux3du = 441, &
& iarr_eflux3dv = 442, iarr_eflux3dw = 443, iarr_ekin0d = 444, &
& iarr_ekin2d = 445, iarr_ekin3d = 446, iarr_enstr0d = 447, iarr_epot0d = 448, &
& iarr_epot2d = 449, iarr_etot0d = 450, iarr_etot2d = 451, iarr_etot3d = 452, &
& iarr_vortic2d = 453, iarr_vortic3d = 454

!---nesting
INTEGER, PARAMETER :: &
& iarr_instbio = 455, iarr_instsed = 456, iarr_inst2dtype = 457, &
& iarr_lbhnstatc = 458, iarr_lbhnstatu = 459, iarr_lbhnstatv = 460, &
& iarr_lbhnstatx = 461, iarr_lbhnstaty = 462, iarr_nestcoords = 463, &
& iarr_nobionst = 464, iarr_nohnstatc = 465, iarr_nohnstatu = 466, &
& iarr_nohnstatv = 467, iarr_nohnstatx = 468, iarr_nohnstaty = 469, &
& iarr_nohnstcprocs = 470, iarr_nohnstglbc = 471, iarr_nohnstglbu = 472, &
& iarr_nohnstglbv = 473, iarr_nohnstglbx = 474, iarr_nohnstglby = 475, &
& iarr_nohnstuvprocs = 476, iarr_nohnstxyprocs = 477, iarr_nosednst = 478, &
& iarr_novnst = 479

!---relaxation zones
INTEGER, PARAMETER :: &
& iarr_idirrlx = 480, iarr_indexrlxatc = 481, iarr_indexrlxatuv = 482, &
& iarr_inodesrlx = 483, iarr_iposrlx = 484, iarr_iprofrlx = 485, &
& iarr_ityprlx = 486, iarr_jposrlx = 487, iarr_ncrlx = 488, iarr_nrrlx = 489, &
& iarr_rlxwghtatc = 490, iarr_rlxwghtatuv = 491

!---elliptic parameters
INTEGER, PARAMETER :: &
& iarr_ellac2d = 492, iarr_ellac3d = 493, iarr_ellcc2d = 494, &
& iarr_ellcc3d = 495, iarr_ellinc2d = 496, iarr_ellinc3d = 497, &
& iarr_ellip2d = 498, iarr_ellip3d = 499, iarr_ellmaj2d = 500, &
& iarr_ellmaj3d = 501, iarr_ellmin2d = 502, iarr_ellmin3d = 503, &
& iarr_ellpha2d = 504, iarr_ellpha3d = 505

!---relative coordinate arrays
INTEGER, PARAMETER :: &
& iarr_icoordC = 506, iarr_icoordCC = 507, iarr_icoordCU = 508, &
& iarr_icoordCV = 509, iarr_icoordU = 510, iarr_icoordUU = 511, &
& iarr_icoordUY = 512, iarr_icoordV = 513, iarr_icoordVV = 514, &
& iarr_icoordVX = 515, iarr_jcoordC = 516, iarr_jcoordCC = 517, &
& iarr_jcoordCU = 518, iarr_jcoordCV = 519, iarr_jcoordU = 520, &
& iarr_jcoordUU = 521, iarr_jcoordUY = 522, iarr_jcoordV = 523, &
& iarr_jcoordVV = 524, iarr_jcoordVX = 525, iarr_weightsC = 526, &
& iarr_weightsCC = 527, iarr_weightsCU = 528, iarr_weightsCV = 529, &
& iarr_weightsU = 530, iarr_weightsUU = 531, iarr_weightsUY = 532, &
& iarr_weightsV = 533, iarr_weightsVV = 534, iarr_weightsVX = 535

!---(absolute) data coordinate arrays
INTEGER, PARAMETER :: &
& iarr_depout = 536, iarr_levout = 537, iarr_time = 538, iarr_xcoord = 539, &
& iarr_xcoordatc = 540, iarr_xcoordatu = 541, iarr_xcoordatv = 542, &
& iarr_xcoordatx = 543, iarr_xcoordaty = 544, iarr_xout = 545, &
& iarr_ycoord = 546, iarr_ycoordatc = 547, iarr_ycoordatu = 548, &
& iarr_ycoordatv = 549, iarr_ycoordatx = 550, iarr_ycoordaty = 551, &
& iarr_yout = 552, iarr_zcmean = 553, iarr_zcoord = 554, iarr_zcoordatc = 555, &
& iarr_zcoordatu = 556, iarr_zcoordatv = 557, iarr_zcoordatx = 558, &
& iarr_zcoordaty = 559, iarr_zetout = 560

!---model parameters
INTEGER, PARAMETER :: &
& iarr_density_ref = 561, iarr_gacc_mean = 562


END MODULE modids
