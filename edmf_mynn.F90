module edmf_mynn_mod

!=======================================================================
! GFDL driver for the EDMF-MYNN (Eddy-Diffusivity/Mass-Flux Mellor–Yamada–Nakanishi–Niino) parameterization.
!
! Author: Yi-Hsuan Chen
!
! For a scientific description, refer to the following publications:

!    Olson, J. B., J. S. Kenyon, W. A. Angevine, J. M. Brown, M. Pagowski, and K. Su?selj, 2019: A description of the MYNN-EDMF scheme and the coupling to other components in WRF–ARW. NOAA Tech. Memo. OAR GSD-61, 42 pp., https://doi.org/10.25923/n9wm-be49.

!    Wu, E., H. Yang, J. Kleissl, K. Suselj, M. J. Kurowski, and J. Teixeira, 2020: On the parameterization of convective downdrafts for marine stratocumulus clouds. Mon. Weather Rev., 148, 1931–1950, https://doi.org/10.1175/MWR-D-19-0292.1.
!
!
! History
!   2020/12/15  ver 1.0
!   2020/12/20  ver 2.0, the mynn and edmf_* programs are taken from 
!                        xanadu_CM4.edmf_mynn/offline-edmf_mynn_v1_loop.F90, 939a7f5a9fe1f13c009b25925a0e4207593155d5
!   2021/01/04  ver 2.1, the mynn and edmf_* programs are taken from 
!                        xanadu_CM4.edmf_mynn/offline-edmf_mynn_v1_loop.F90, 44a5ce610f3ba13f36ca421d4052c2543668df95
!   2021/01/05  ver 2.2, the mynn and edmf_* programs are taken from 
!                        xanadu_CM4.edmf_mynn/offline-edmf_mynn_v1_loop.F90, fe4d5aa337747eb04ec6cefbc0cf0cbb738fdfed
!=======================================================================

use           mpp_mod, only: input_nml_file

use           fms_mod, only: file_exist, open_namelist_file,       &
                             error_mesg, FATAL, close_file, note,  &
                             check_nml_error, mpp_pe, mpp_root_pe, &
                             write_version_number, stdlog, stdout, &
                             mpp_chksum

use     constants_mod, only: vonkarm, rdgas, rvgas, kappa, grav, &
                             cp_air, dens_h2o, hlv, hlf, hls, pi  !, con_cliq, con_csol. con_cliq, con_csol are not accessble?

use  diag_manager_mod, only: register_diag_field, send_data

use  time_manager_mod, only: time_type

use  sat_vapor_pres_mod,    only: compute_qs

use random_numbers_mod,    only:  randomNumberStream,   &
                                  getRandomNumbers

use random_number_streams_mod, only: random_number_streams_init, &
                                     get_random_number_streams, &
                                     random_number_streams_end

use physics_types_mod,       only: alloc_physics_tendency_type, &
                                   physics_tendency_type, &
                                   phys_mp_exch_type, &
                                   phys2cosp_type, precip_flux_type, &
                                   physics_tendency_block_type, &
                                   physics_type, &
                                   physics_control_type, &
                                   physics_input_block_type, &
                                   dealloc_physics_tendency_type

use   field_manager_mod, only: MODEL_ATMOS

use  tracer_manager_mod, only: get_number_tracers, get_tracer_index,           &
                               get_tracer_names
use moist_proc_utils_mod, only : rh_calc

use random_numbers_mod,    only:  randomNumberStream,   &
                                  getRandomNumbers

use random_number_streams_mod, only: random_number_streams_init, &
                                     get_random_number_streams, &
                                     random_number_streams_end
 
!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------

public :: edmf_mynn_init, edmf_mynn_end, &
          edmf_input_type, edmf_output_type, edmf_ls_mp_type, &
          edmf_mynn_driver

!---------------------------------------------------------------------

!---------------------------------------------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
logical            :: module_is_initialized = .false.
 
!---------------------------------------------------------------------
!  derived type definition for mf_tendency_type variable
!---------------------------------------------------------------------

!==================
type edmf_input_type
  integer,              allocatable ::   &
    IDS,IDE,JDS,JDE,KDS,KDE,             &
    IMS,IME,JMS,JME,KMS,KME,             &
    ITS,ITE,JTS,JTE,KTS,KTE
  real,                 allocatable   :: &
    delt, dx
  real, dimension(:,:), allocatable   :: &   ! INPUT, DIMENSION(IMS:IME,JMS:JME)
    xland, ust, ch, rmol, ts, qsfc, qcg, ps, hfx, qfx, wspd, uoce, voce, vdfg, znt

  real, dimension(:,:,:), allocatable :: &   ! INPUT, DIMENSION(IMS:IME,KMS:KME,JMS:JME)
    dz, u, v, w, th, qv, p, exner, rho, T3D, qa, qc, ql, qi, qnc, qni

  real, dimension(:,:), allocatable   :: &   ! DIMENSION(IMS:IME,JMS:JME), diagnostic purpose
    u_star_star, shflx_star, lhflx_star, w1_thv1_surf_star, w1_th1_surf_star, w1_q1_surf_star, Obukhov_length_star,  &
    u_star_updated, shflx_updated, lhflx_updated, w1_thv1_surf_updated, w1_th1_surf_updated, w1_q1_surf_updated, Obukhov_length_updated

  real, dimension(:,:,:), allocatable :: &   ! INPUT, DIMENSION(IMS:IME,KMS:KME,JMS:JME)
    Qke, Sh3D, el_pbl, cldfra_bl, qc_bl      ! semi-prognostic variables


end type edmf_input_type

!==================
type edmf_output_type
  real, dimension(:,:),   allocatable :: &   ! INPUT/OUTPUT, DIMENSION(IMS:IME,JMS:JME)
    pblh,wstar,delta

  integer, dimension(:,:),   allocatable :: &   ! INPUT/OUTPUT, DIMENSION(IMS:IME,JMS:JME)
    KPBL, nupdraft, ktop_shallow
  
  real, dimension(:,:,:), allocatable :: &   ! INPUT/OUTPUT, DIMENSION(IMS:IME,KMS:KME,JMS:JME)
    Qke, Tsq, Qsq, Cov,         &
    RUBLTEN, RVBLTEN, RTHBLTEN, RQVBLTEN, RQLBLTEN, RQIBLTEN, RQNIBLTEN, RTHRATEN, &
    RCCBLTEN, RTHLBLTEN, RQTBLTEN,   &
    edmf_a, edmf_w, edmf_qt, edmf_thl, edmf_ent, edmf_det, edmf_qc, edmf_a_dd,edmf_w_dd,edmf_qt_dd,edmf_thl_dd,edmf_ent_dd,edmf_qc_dd, &
    edmf_debug1,edmf_debug2,edmf_debug3,edmf_debug4, &
    Q_ql, Q_qi, Q_qa, &
    Q_ql_adv, Q_qi_adv, Q_qa_adv, &
    Q_ql_ent, Q_qi_ent, Q_qa_ent, &
    Q_ql_eddy, Q_qi_eddy, Q_qa_eddy, &
    Q_ql_det, Q_qi_det, Q_qa_det, &
    Q_ql_sub, Q_qi_sub, Q_qa_sub, &
    mynn_ql, qc_bl, cldfra_bl, el_pbl, Sh3D

  real, dimension(:,:),   allocatable :: &   ! OUTPUT, DIMENSION(IMS:IME,JMS:JME)
    maxmf

  real, dimension(:,:,:), allocatable :: &   ! OUTPUT, DIMENSION(IMS:IME,KMS:KME,JMS:JME)
    exch_h, exch_m, &
    qWT, qSHEAR, qBUOY, qDISS, dqke

  real, dimension(:,:,:), allocatable :: &   ! OUTPUT, DIMENSION(IMS:IME,KMS:KME,JMS:JME)
    a_moist_half, mf_moist_half, qv_moist_half, a_moist_full, mf_moist_full, qv_moist_full, &
    a_dry_half, mf_dry_half, qv_dry_half, a_dry_full, mf_dry_full, qv_dry_full, &
    mf_all_half, mf_all_full

  real, dimension(:,:,:), allocatable :: &   ! diagnostic purpose, not used by mynn 
    t_before_mix, q_before_mix, qa_before_mix, ql_before_mix, qi_before_mix, thl_before_mix, qt_before_mix, rh_before_mix, th_before_mix, &
    qa_before_pdf, ql_before_pdf, qi_before_pdf, &
    t_after_mix, q_after_mix, qa_after_mix, ql_after_mix, qi_after_mix, thl_after_mix, qt_after_mix, rh_after_mix, th_after_mix

end type edmf_output_type

!==================
type am4_edmf_output_type

  real, dimension(:,:,:),   allocatable :: &   ! OUTPUT, DIMENSION(nlon, nlat, nfull)
    thl_edmf, qt_edmf,  & ! diagnostic purpose
    tke, Tsq, Cov_thl_qt, udt_edmf, vdt_edmf, tdt_edmf, qdt_edmf, qidt_edmf, qldt_edmf, qadt_edmf, qndt_edmf, &  ! outputs from EDMF-MYNN scheme
    qtdt_edmf, thldt_edmf, &
    diff_t_edmf, diff_m_edmf, &
    cldfra_bl, qc_bl, el_edmf, &
    Q_ql, Q_qi, Q_qa, &
    edmf_a, edmf_w, edmf_qt, edmf_thl, edmf_ent, edmf_qc                                   !

  real, dimension(:,:,:), allocatable :: &   ! OUTPUT, DIMENSION(nlon, nlat, nfull/nhalf)
    a_moist_half, mf_moist_half, qv_moist_half, a_moist_full, mf_moist_full, qv_moist_full, &
    a_dry_half, mf_dry_half, qv_dry_half, a_dry_full, mf_dry_full, qv_dry_full, &
    mf_all_half, mf_all_full

  real, dimension(:,:),     allocatable :: &   ! OUTPUT, DIMENSION(nlon, nlat)
    pbltop

  integer, dimension(:,:),     allocatable :: &   ! OUTPUT, DIMENSION(nlon, nlat)
    kpbl_edmf

  real, dimension(:,:,:), allocatable :: &   ! diagnostic purpose, not used by mynn 
    t_input, q_input, qa_input, ql_input, qi_input, thl_input, qt_input, rh_input, th_input, &
    t_before_mix, q_before_mix, qa_before_mix, ql_before_mix, qi_before_mix, thl_before_mix, qt_before_mix, rh_before_mix, th_before_mix, &
    qa_before_pdf, ql_before_pdf, qi_before_pdf, &
    t_after_mix, q_after_mix, qa_after_mix, ql_after_mix, qi_after_mix, thl_after_mix, qt_after_mix, rh_after_mix, th_after_mix, &
    rh    ! relative humidity

end type am4_edmf_output_type

!==================
type edmf_ls_mp_type

  real, dimension(:,:,:),   allocatable :: &   ! OUTPUT, DIMENSION(nlon, nlat, nfull)
    qadt_edmf, qldt_edmf, qidt_edmf, &         ! qadt_edmf: cloud fraction tendency from EDMF (1/s)
                                               ! qldt_edmf: cloud liquid specific humidity tendency from EDMF (kg/kg/s)
                                               ! qidt_edmf: cloud ice    specific humidity tendency from EDMF (kg/kg/s)
    dqa_edmf,  dql_edmf, dqi_edmf              ! incrementation changes in a time step
                                               !   i.e. dqa_edmf = qadt_edmf * dt and so on 

  integer, allocatable :: &
    option_edmf2ls_mp

end type edmf_ls_mp_type
!==================

!---------------------------------------------------------------------
! --- Constants
!---------------------------------------------------------------------

   real, parameter :: p00    = 1000.0e2     ! 1000hPa
   real, parameter :: p00inv = 1./p00
   real, parameter :: d622   = rdgas/rvgas
   real, parameter :: d378   = 1.-d622
   real, parameter :: d608   = d378/d622

!   Description of physical constants in WRF:
!     http://math.ucdenver.edu/~farguella/wrf-browsers/WRF-Fire-merge/fuel-moisture-model/html_code/share/module_model_constants.F.html

   !REAL    , PARAMETER :: karman       = 0.4        ! von Karman constant
   !REAL    , PARAMETER :: r_d          = 287.       ! gas constant of dry air (J deg^-1 kg^-1)
   !REAL    , PARAMETER :: g            = 9.81       ! acceleration due to gravity (m {s}^-2)
   !REAL    , PARAMETER :: r_v          = 461.6      ! gas constant for water vapor (J deg^-1 kg^-1)
   !REAL    , PARAMETER :: cliq         = 4190.      ! specific heat of liquid water at 0^oC
   !REAL    , PARAMETER :: Cice         = 2106.      !       ! specific heat of ice at 0^oC
   !REAL    , PARAMETER :: XLS          = 2.85E6     ! latent heat of sublimation of water at 0^oC (J kg^-1) 
   !REAL    , PARAMETER :: XLV          = 2.5E6      ! latent heat of vaporization of water at 0^oC (J kg^-1)
   !REAL    , PARAMETER :: XLF          = 3.50E5     ! latent heat of fusion of water at 0^oC (J kg^-1)
   REAL    , PARAMETER :: karman       = vonkarm        ! von Karman constant
   REAL    , PARAMETER :: g            = grav       ! acceleration due to gravity (m {s}^-2)
   REAL    , PARAMETER :: r_d          = rdgas       ! gas constant of dry air (J deg^-1 kg^-1)
   REAL    , PARAMETER :: cp           = 7.*r_d/2.
   REAL    , PARAMETER :: r_v          = rvgas      ! gas constant for water vapor (J deg^-1 kg^-1)
   REAL    , PARAMETER :: cpv          = 4.*r_v
   !REAL    , PARAMETER :: cliq         = con_cliq     ! specific heat of liquid water at 0^oC
   !REAL    , PARAMETER :: Cice         = con_csol      !       ! specific heat of ice at 0^oC
   REAL    , PARAMETER :: cliq         = 4.1855e+3     ! specific heat of liquid water at 0^oC
   REAL    , PARAMETER :: Cice         = 2.1060e+3      !       ! specific heat of ice at 0^oC
   REAL    , PARAMETER :: rcp          = r_d/cp
   REAL    , PARAMETER :: XLS          = hls      ! latent heat of sublimation of water at 0^oC (J kg^-1) 
   REAL    , PARAMETER :: XLV          = hlv      ! latent heat of vaporization of water at 0^oC (J kg^-1)
   REAL    , PARAMETER :: XLF          = hlf      ! latent heat of fusion of water at 0^oC (J kg^-1)
   REAL    , PARAMETER :: p1000mb      = 100000.
   REAL    , PARAMETER :: rvovrd       = r_v/r_d
   REAL    , PARAMETER :: SVP1         = 0.6112     ! constant for saturation vapor pressure calculation (dimensionless)
   REAL    , PARAMETER :: SVP2         = 17.67      ! constant for saturation vapor pressure calculation (dimensionless)
   REAL    , PARAMETER :: SVP3         = 29.65      ! constant for saturation vapor pressure calculation (K)
   REAL    , PARAMETER :: SVPT0        = 273.15     ! constant for saturation vapor pressure calculation (K)
   REAL    , PARAMETER :: EP_1         = R_v/R_d-1.
   REAL    , PARAMETER :: EP_2         = R_d/R_v
   REAL    , PARAMETER :: rhoair0      = 1.28        ! density of dry air at 0^oC and 1000mb pressure (kg m^-3)

! The following depends on the microphysics scheme used:
   INTEGER , PARAMETER :: param_first_scalar = 1
   INTEGER , PARAMETER :: p_qc = 2
   INTEGER , PARAMETER :: p_qr = 0
   INTEGER , PARAMETER :: p_qi = 2
   INTEGER , PARAMETER :: p_qs = 0
   INTEGER , PARAMETER :: p_qg = 0
   INTEGER , PARAMETER :: p_qnc= 0
   INTEGER , PARAMETER :: p_qni= 0

!<--- yhc note: these are pocied from the WRF MYNN-EDMF module
! The parameters below depend on stability functions of module_sf_mynn.
  REAL, PARAMETER :: cphm_st=5.0, cphm_unst=16.0, &
                     cphh_st=5.0, cphh_unst=16.0

  REAL, PARAMETER :: xlvcp=xlv/cp, xlscp=(xlv+xlf)/cp, ev=xlv, rd=r_d, &
       &rk=cp/rd, svp11=svp1*1.e3, p608=ep_1, ep_3=1.-ep_2

  REAL, PARAMETER :: tref=300.0     ! reference temperature (K)
  REAL, PARAMETER :: TKmin=253.0    ! for total water conversion, Tripoli and Cotton (1981)
  REAL, PARAMETER :: tv0=p608*tref, tv1=(1.+p608)*tref, gtr=g/tref

! Closure constants
  REAL, PARAMETER :: &
       &vk  = karman, &
       &pr  =  0.74,  &
       &g1  =  0.235, &  ! NN2009 = 0.235
       &b1  = 24.0, &
       &b2  = 15.0, &    ! CKmod     NN2009
       &c2  =  0.729, &  ! 0.729, & !0.75, &
       &c3  =  0.340, &  ! 0.340, & !0.352, &
       &c4  =  0.0, &
       &c5  =  0.2, &
       &a1  = b1*( 1.0-3.0*g1 )/6.0, &
!       &c1  = g1 -1.0/( 3.0*a1*b1**(1.0/3.0) ), &
       &c1  = g1 -1.0/( 3.0*a1*2.88449914061481660), &
       &a2  = a1*( g1-c1 )/( g1*pr ), &
       &g2  = b2/b1*( 1.0-c3 ) +2.0*a1/b1*( 3.0-2.0*c2 )

  REAL, PARAMETER :: &
       &cc2 =  1.0-c2, &
       &cc3 =  1.0-c3, &
       &e1c =  3.0*a2*b2*cc3, &
       &e2c =  9.0*a1*a2*cc2, &
       &e3c =  9.0*a2*a2*cc2*( 1.0-c5 ), &
       &e4c = 12.0*a1*a2*cc2, &
       &e5c =  6.0*a1*a1

! Constants for min tke in elt integration (qmin), max z/L in els (zmax), 
! and factor for eddy viscosity for TKE (Kq = Sqfac*Km):
  REAL, PARAMETER :: qmin=0.0, zmax=1.0, Sqfac=3.0
! Note that the following mixing-length constants are now specified in mym_length
!      &cns=3.5, alp1=0.23, alp2=0.3, alp3=3.0, alp4=10.0, alp5=0.4

! Constants for gravitational settling
!  REAL, PARAMETER :: gno=1.e6/(1.e8)**(2./3.), gpw=5./3., qcgmin=1.e-8
  REAL, PARAMETER :: gno=1.0  !original value seems too agressive: 4.64158883361278196
  REAL, PARAMETER :: gpw=5./3., qcgmin=1.e-8, qkemin=1.e-12

! Constants for cloud PDF (mym_condensation)
  REAL, PARAMETER :: rr2=0.7071068, rrp=0.3989423

! 'parameters' for Poisson distribution (StEM EDMF scheme)
  REAL, PARAMETER  :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0

  !Use Canuto/Kitamura mod (remove Ric and negative TKE) (1:yes, 0:no)
  !For more info, see Canuto et al. (2008 JAS) and Kitamura (Journal of the 
  !Meteorological Society of Japan, Vol. 88, No. 5, pp. 857-864, 2010).
  !Note that this change required further modification of other parameters
  !above (c2, c3). If you want to remove this option, set c2 and c3 constants 
  !(above) back to NN2009 values (see commented out lines next to the
  !parameters above). This only removes the negative TKE problem
  !but does not necessarily improve performance - neutral impact.
  REAL, PARAMETER :: CKmod=1.

  !Use Ito et al. (2015, BLM) scale-aware (0: no, 1: yes). Note that this also has impacts
  !on the cloud PDF and mass-flux scheme, using Honnert et al. (2011) similarity function
  !for TKE in the upper PBL/cloud layer.
  !REAL, PARAMETER :: scaleaware=0.


  !Adding top-down diffusion driven by cloud-top radiative cooling
  INTEGER, PARAMETER :: bl_mynn_topdown = 0

  !option to print out more stuff for debugging purposes
  LOGICAL, PARAMETER :: debug_code = .false.

! JAYMES-
! Constants used for empirical calculations of saturation
! vapor pressures (in function "esat") and saturation mixing ratios
! (in function "qsat"), reproduced from module_mp_thompson.F, 
! v3.6 
  REAL, PARAMETER:: J0= .611583699E03
  REAL, PARAMETER:: J1= .444606896E02
  REAL, PARAMETER:: J2= .143177157E01
  REAL, PARAMETER:: J3= .264224321E-1
  REAL, PARAMETER:: J4= .299291081E-3
  REAL, PARAMETER:: J5= .203154182E-5
  REAL, PARAMETER:: J6= .702620698E-8
  REAL, PARAMETER:: J7= .379534310E-11
  REAL, PARAMETER:: J8=-.321582393E-13

  REAL, PARAMETER:: K0= .609868993E03
  REAL, PARAMETER:: K1= .499320233E02
  REAL, PARAMETER:: K2= .184672631E01
  REAL, PARAMETER:: K3= .402737184E-1
  REAL, PARAMETER:: K4= .565392987E-3
  REAL, PARAMETER:: K5= .521693933E-5
  REAL, PARAMETER:: K6= .307839583E-7
  REAL, PARAMETER:: K7= .105785160E-9
  REAL, PARAMETER:: K8= .161444444E-12
! end-
!---> yhc note: these are pocied from the WRF MYNN-EDMF module

!---------------------------------------------------------------------
! --- Namelist
!---------------------------------------------------------------------
  CHARACTER*128 :: mynn_message
  INTEGER, PARAMETER :: kdebug=27

  INTEGER :: mynn_level        = 2         ! level <>3 to do the level 2.5 MYNN model; if 3 then will do level 3
  INTEGER :: grav_settling     = 0         ! 1 - the cloud/fog will experience gravitational settling, 0 - it does not
  INTEGER :: bl_mynn_tkebudget = 0         ! if 1 the budget terms in the TKE equation are allocated for output (WRF), if 0 then not
  INTEGER :: bl_mynn_cloudpdf  = 1         ! define the type of the subgrid PDF for cloud computation, 1= Nakanishi & Niino 2004
  INTEGER :: bl_mynn_mixlength = 2         ! defines the ED mixing length formulation
  INTEGER :: bl_mynn_edmf      = 0         ! controls  the version of the EDMF to be called 3=JPLedmf documented by Wu et al., 2020
                                           ! set “bl_mynn_edmf=0” and “bl_mynn_edmf_dd<>1” then the scheme will be MYNN only.
  INTEGER :: bl_mynn_edmf_dd   = 0         ! 0 - no downdrafts, 1 - Wu et al., 2020 downdrafts 
  REAL    :: bl_mynn_edmf_Lent = 0.        ! dummy argument
  LOGICAL :: bl_mynn_tkeadvect = .false.   ! 1 - advecting TKE, 0 - not advecting
  INTEGER :: bl_mynn_edmf_mom  = 1         ! 1- mixing of momentum by full EDMF, 0 - momenum is mixed only by ED
  INTEGER :: bl_mynn_edmf_tke  = 0         ! does not effect JPLedmf
  INTEGER :: bl_mynn_edmf_part = 1         ! partition area into updrafts and remaining environment
  INTEGER :: bl_mynn_cloudmix  = 1         ! 1 - cloud species are mixed separately, 0 - not
  INTEGER :: bl_mynn_mixqt     = 2         ! will mix moist conserved variables, after mixing invokes PDF cloud scheme to convert moist variables to dry
  INTEGER :: bl_mynn_stabfunc  = 0       !  option for stability function. 0 - MYNN, 1 - Suselj et al. (2019)
  INTEGER :: icloud_bl         = 1         ! 1, cloud cover and liquid water from the PDF scheme will be on the cloud output
  INTEGER :: spp_pbl           = 0         ! 1 stochastic perturbation to condensation  

  integer :: initflag = 1                 ! (when 1 it initializes TKE using level 2 MYNN scheme)
  logical :: FLAG_QI  = .true.            ! (flags for whether cloud and ice mixing rations and number concentrations are mixed separately)
  logical :: FLAG_QNI = .false.           ! all false except FLAG_QI that Kay Suselj said "bl_mynn_cloudmix=1 and FLAG_QI=.true." on 12/11/2020
  logical :: FLAG_QC  = .false.
  logical :: FLAG_QNC = .false.

  logical :: expmf = .true.              ! .true.  ... explicit mass-flux, .false. .... implicit mass-flux
  real    :: upwind = 1.                 ! upwind=1. ... use upwind approximation for mass-flux calculation
                                         ! upwind=0.5 ... use centered difference for mass-flux calculation
                                         ! explicit mass-flux can use either upwind or centered-difference
                                         ! implicit mass-flux must uses the centered differencing method.  
  real    :: L0  = 100.                  ! entrainemnt rate parameter
  integer :: NUP = 100                   ! the number of updrafts
  real    :: UPSTAB = 1.                 ! stability parameter for massflux, (mass flux is limited so that dt/dz*a_i*w_i<UPSTAB)

  integer :: edmf_type=0                 ! =0, the standard MYNN code, in which the PDF cloud scheme before mixing and after the mixing and compute the tendencies of liquid and cloud properties from the differences between these two.
                                         ! =1, tendencies of moist variables from the PDF scheme after mixing and from the input values (from Tiedtke, presumably)
                                         ! =2, approximate ql and qi tendencies from AM4 ED and MF terms, and then recover T and qv tendencies 

  logical :: do_qdt_same_as_qtdt = .false.  ! = .true., for testing purposes, evaporate/condensate the liquid and ice water that is produced during mixing. 

  integer :: ii_write = -999                    ! i index for column written out. Set to 0 if you want to write out in SCM
  integer :: jj_write = -999                    ! j index for column written out. Set to 0 if you want to write out in SCM
  real    :: lat_write = -999.99                ! latitude  (radian) for column written out
  real    :: lon_write = -999.99                ! longitude (radian) for column written out
  real    :: lat_range = 0.000001
  real    :: lon_range = 0.000001
  logical :: do_writeout_column_nml = .false.   ! control whether writing out the column
  logical :: do_check_consrv = .false.          ! control whether writing out the water/heat conservation
  logical :: do_check_realizability = .false.   ! whether enables the realizability limiter that forces the tracers concentration
                                                ! (q_t, q_v, q_l, q_i, q_a) not become negative after one time step 
  logical :: do_check_ent_det = .false.         ! control whether writing out entrainment and detrainment rate for debugging
  logical :: do_stop_run = .false.              ! whether to stop the simulation
  real    :: qke_min = -1.                      ! qke=2*tke. If qke < qke_min, set all EDMF tendencies to zeros
                                                !   set qke_min>0 may remove energy/water away and cause water mass is not conserved
  real    :: tracer_min = 1.E-10                ! make sure tracer value is not smaller than tracer_min
                                                ! 1.E-10 is same as qmin in lscloud_driver
  real    :: qdt_min    = 1.E-15                ! minimum value for the tendencies used by the realizability limiter
                                                ! This is to avoid nearly zero tracer concentration with extremely small negative tendency,
                                                ! which leads to extremely small rescaled ratio
  integer :: do_option_edmf2ls_mp = 0           ! option to include EDMF cloud tendencies terms into Tiedtke
                                                ! =0, add EDMF terms to accumulated ql and qi tendnecies
                                                ! =1, add EDMF terms to Tiedtke and keep Tiedtke terms except turbulence heating
                                                ! =2, add EDMF terms to Tiedtke and remove large-scale and coud erosion terms
                                                ! =3, evaporate ql and qi and put these water back to vapor and change temperature accordingly
  integer :: Qx_MF=1                            ! =1, using eddy-divergence/source form to compute MF tendencies on grid-scale cloud
                                                ! =2, using detrainment/subsidence form to compute MF tendencies on grid-scale cloud

  integer :: Qx_numerics=1                      ! =1, use upwind approximation in computing MF eddy-flux/source and subsidence/detrainment forms

  integer :: option_up_area=1                   !=1, constant updraft area assumption, i.e. the updraft area remains the same with height
                                                !=2, updraft area is changed with height to maintain mass conservation in plumes
                                                !=3, if detrainment rate>=0, assume updraft area is constant.
                                                !    otherwise, vary the area to satisfy the mass continuity equation
  logical :: do_use_tau = .true.                ! .true.  : use the T,q at the current step
                                                ! .false. : use the updated T,q
  logical :: do_simple =.false.                 ! do_simple = switch to turn on alternative definition of specific
                                                !             humidity. When true, specific humidity =
                                                !             (rdgas/rvgas)*esat/pressure
                                                ! same setting as module moist_processes_mod

  integer :: option_ent=1                       ! =1, use original stochastic entrainment formula
                                                ! =2, separate entrainment into stochastic and deterministic parts
                                                !     controlling by alpha_st (1 completely stochastic, 0 completely deterministic)
  real    :: alpha_st=0.5
  real    :: rc_MF = 10.e-6                     ! assumed cloud droplet radius in plumes (meters)

  character*20 :: do_debug_option = ""          ! debug purpose

  !character*20 :: option_surface_flux = "star"      ! surface fluxes are determined by "star" quantities, i.e. u_star, q_star, and b_star
  character*20 :: option_surface_flux = "updated"  ! surface fluxes are determined by "updated" quantities, i.e. u_flux, v_flux, shflx, and lh flx
  real    :: tdt_max     = 2000. ! K/day
  logical :: do_limit_tdt = .false.
  real    :: tdt_limit =   200. ! K/day

  logical :: do_pblh_constant = .false.    ! fix PBL depth for testing
  real    :: fixed_pblh       = 2500. 

  real    :: sgm_factor = 100.                  ! factor in computing sigma_s in MYNN

  character*20 :: option_stoch_entrain = "Poisson_knuth"  ! option for random number generator (RNG)
                                                          !   Poisson: using the Fortran intrinsic RNG
                                                          !   Poisson_knuth: using AM4 RNG
  integer :: option_rng = 1  ! in Poisson_knuth, 0 - using Fortran intrisic random_number, 1 - using AM4 RNG

  integer :: option_pblh_MF = 0   ! 0: whenever w'thv'>0, call MF
                                  ! 1: only when w'thv'>0 and PBL height > z0 (=50m), call MF. 
                                  !    This is to prevent activating MF in very stable conditions.

namelist / edmf_mynn_nml /  mynn_level, bl_mynn_edmf, bl_mynn_edmf_dd, expmf, upwind, do_qdt_same_as_qtdt, bl_mynn_mixlength, bl_mynn_stabfunc, &
                            L0, NUP, UPSTAB, edmf_type, qke_min, &
                            option_surface_flux, &
                            tdt_max, do_limit_tdt, tdt_limit, do_pblh_constant, fixed_pblh, sgm_factor, rc_MF, &  ! for testing, no need any more 2021-08-04
                            option_stoch_entrain, option_rng, option_pblh_MF, &
                            do_option_edmf2ls_mp, do_use_tau, Qx_MF, Qx_numerics, option_up_area, option_ent, alpha_st, &
                            do_debug_option, do_stop_run, do_writeout_column_nml, do_check_consrv, do_check_realizability, &
                            ii_write, jj_write, lat_write, lon_write

!---------------------------------------------------------------------
!--- Diagnostic fields       
!---------------------------------------------------------------------

character(len=10) :: mod_name = 'edmf_mynn'
real              :: missing_value = -999.

integer :: nsphum, nql, nqi, nqa, nqn, nqni  ! tracer indices for stratiform clouds
integer :: nQke, nSh3D, nel_pbl, ncldfra_bl, nqc_bl ! tracer index for EDMF-MYNN tracers
integer :: ntp          ! number of prognostic tracers

integer :: id_u_flux, id_v_flux, id_u_star_updated, id_shflx_star, id_lhflx_star, id_w1_thv1_surf_star, id_w1_thv1_surf_updated, id_Obukhov_length_star, id_Obukhov_length_updated, id_tke_edmf, id_Tsq, id_Cov_thl_qt, id_udt_edmf, id_vdt_edmf, id_tdt_edmf, id_qdt_edmf, id_qidt_edmf, id_qadt_edmf, id_qldt_edmf, id_qndt_edmf, id_edmf_a, id_edmf_w, id_edmf_qt, id_edmf_thl, id_edmf_ent, id_edmf_det, id_edmf_qc, id_thl_edmf, id_qt_edmf, id_cldfra_bl, id_qc_bl, id_z_pbl, id_z_pbl_edmf, id_qtdt_edmf, id_thldt_edmf, id_diff_t_edmf, id_diff_m_edmf, id_el_edmf

integer :: id_t_input, id_q_input, id_qa_input, id_ql_input, id_qi_input, id_thl_input, id_qt_input, id_rh_input, id_th_input, id_t_before_mix, id_q_before_mix, id_qa_before_mix, id_ql_before_mix, id_qi_before_mix, id_thl_before_mix, id_qt_before_mix, id_rh_before_mix, id_th_before_mix, id_t_after_mix, id_q_after_mix, id_qa_after_mix, id_ql_after_mix, id_qi_after_mix, id_thl_after_mix, id_qt_after_mix, id_rh_after_mix, id_th_after_mix, id_qa_before_pdf, id_ql_before_pdf, id_qi_before_pdf 

integer :: id_tdt_edmf_orig, id_qdt_edmf_orig, id_qadt_edmf_orig, id_qidt_edmf_orig, id_qldt_edmf_orig

integer :: id_qldt_edmf_ED, id_qldt_edmf_MF, id_qidt_edmf_ED, id_qidt_edmf_MF, id_qadt_edmf_ED, id_qadt_edmf_MF

integer :: id_qldt_edmf_MF_adv, id_qldt_edmf_MF_eddy, id_qldt_edmf_MF_ent, id_qidt_edmf_MF_adv, id_qidt_edmf_MF_eddy, id_qidt_edmf_MF_ent, id_qadt_edmf_MF_adv, id_qadt_edmf_MF_eddy, id_qadt_edmf_MF_ent, id_qadt_edmf_MF_det, id_qldt_edmf_MF_det, id_qidt_edmf_MF_det, id_qadt_edmf_MF_sub, id_qldt_edmf_MF_sub, id_qidt_edmf_MF_sub

integer :: id_a_moist_half, id_a_moist_full, id_mf_moist_half, id_mf_moist_full, id_qv_moist_half, id_qv_moist_full, id_a_dry_half, id_a_dry_full, id_mf_dry_half, id_mf_dry_full, id_qv_dry_half, id_qv_dry_full, id_mf_all_half, id_mf_all_full

integer :: id_num_updraft, id_num_det, id_num_ndet_zent, id_num_ndet_pent

integer :: id_rlz_ratio, id_rlz_tracer
!---------------------------------------------------------------------

  contains

!#######################################################################

!subroutine edmf_mynn_init(lonb, latb, axes, time, id, jd, kd, edmf2ls_mp)
subroutine edmf_mynn_init(lonb, latb, axes, time, id, jd, kd)

!-----------------------------------------------------------------------
!  (Intent in)
!-----------------------------------------------------------------------
!   latb, lonb  - latitudes and longitudes at grid box corners
!   axes, time  - variables needed for netcdf diagnostics
!   id, jd, kd  - size of the first 3 dimensions

 integer,              intent(in) :: id, jd, kd, axes(4)
 type(time_type),      intent(in) :: time
 real, dimension(:,:), intent(in) :: lonb, latb

!-----------------------------------------------------------------------
!  (Intent out)
!-----------------------------------------------------------------------
!  type(edmf_ls_mp_type), intent(out) :: edmf2ls_mp

!-----------------------------------------------------------------------
!  (Intent local)
!-----------------------------------------------------------------------
 integer, dimension(3) :: full = (/1,2,3/), half = (/1,2,4/)
 integer :: unit, io, ierr, logunit

!=======================================================================

  if ( module_is_initialized ) return

!-----------------------------------------------------------------------
! --- Read namelist
!-----------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=edmf_mynn_nml, iostat=io)
   ierr = check_nml_error(io,'edmf_mynn_nml')
#else   
  if( file_exist( 'input.nml' ) ) then
! -------------------------------------
   unit = open_namelist_file( )
   ierr = 1
   do while( ierr .ne. 0 )
   READ ( unit,  nml = edmf_mynn_nml, iostat = io, end = 10 ) 
   ierr = check_nml_error (io, 'edmf_mynn_nml')
   end do
10 continue
   call close_file( unit )
! -------------------------------------
  end if
#endif

!-----------------------------------------------------------------------
! --- Check namelist
!-----------------------------------------------------------------------

  if ( trim(option_surface_flux) /= 'star' .and. &
       trim(option_surface_flux) /= 'updated' )       then 
    call error_mesg( ' edmf_mynn',     &
                     ' option_surface_flux must be "star" or "updated" in the edmf_mynn_nml',&
                     FATAL )
  end if

  if (.not.expmf .and. upwind.ne.0.5) then
    call error_mesg( ' edmf_mynn',     &
                     ' when expmf=.false., upwind must be 0.5 (implicit MF with centered-diff)',&
                     FATAL )
  endif

  if ( edmf_type .gt. 2 ) then
    call error_mesg( ' edmf_mynn',     &
                     ' edmf_type must be either 0, 1, or 2',&
                     FATAL )
  endif

  if ( alpha_st.gt.1. .or. alpha_st.lt.0. ) then
    call error_mesg( ' edmf_mynn',     &
                     ' alpha_st  must be between 0 to 1',&
                     FATAL )
  endif

  if ( option_ent.gt.2 .or. option_ent.lt.0) then
    call error_mesg( ' edmf_mynn',     &
                     ' option_ent must be 1 or 2',&
                     FATAL )
  endif

  if (trim(option_stoch_entrain) /= 'Poisson' .and. trim(option_stoch_entrain) /= 'Poisson_knuth') then
    call error_mesg( ' edmf_mynn',     &
                     '  option_stoch_entrain must be Poisson or Poisson_knuth',&
                     FATAL )
  endif

!  if (     do_option_edmf2ls_mp.ne.0    &
!     .and. do_option_edmf2ls_mp.ne.1    &
!     .and. do_option_edmf2ls_mp.ne.2    &
!     .and. do_option_edmf2ls_mp.ne.3    &
!     ) then
!    call error_mesg( ' edmf_mynn',     &
!                     ' do_option_edmf2ls_mp must be 0,1,2, or 3',&
!                     FATAL )
!  endif

!  if ( trim(option_solver) /= 'explicit' .and. &
!       trim(option_solver) /= 'implicit' )       then 
!    call error_mesg( ' edmf_mynn',     &
!                     ' option_solver must be "explicit" or "implicit" in the edmf_mynn_nml',&
!                     FATAL )
!  end if

!-----------------------------------------------------------------------
! --- Output version
!-----------------------------------------------------------------------

  if ( mpp_pe() == mpp_root_pe() ) then
       call write_version_number(version, tagname)
       logunit = stdlog()
       WRITE( logunit, nml = edmf_mynn_nml ) 
  endif

!-----------------------------------------------------------------------
! --- Get tracer indices 
!-----------------------------------------------------------------------
  call get_number_tracers (MODEL_ATMOS, num_prog=ntp)

  nsphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
  nql    = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
  nqi    = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
  nqa    = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
  nqn    = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
  nqni   = get_tracer_index ( MODEL_ATMOS, 'ice_num' )

  nQke        = get_tracer_index ( MODEL_ATMOS, 'Qke' )
  nSh3D       = get_tracer_index ( MODEL_ATMOS, 'Sh3D' )
  nel_pbl     = get_tracer_index ( MODEL_ATMOS, 'el_pbl' )
  ncldfra_bl  = get_tracer_index ( MODEL_ATMOS, 'cldfra_bl' )
  nqc_bl      = get_tracer_index ( MODEL_ATMOS, 'qc_bl' )

  if (nQke <= ntp) call error_mesg ('edmf_mynn_driver_mod', &
                    'Qke can not be a prognostic tracer', FATAL)
  if (nSh3D <= ntp) call error_mesg ('edmf_mynn_driver_mod', &
                    'Sh3D can not be a prognostic tracer', FATAL)
  if (nel_pbl <= ntp) call error_mesg ('edmf_mynn_driver_mod', &
                    'el_pbl can not be a prognostic tracer', FATAL)
  if (ncldfra_bl <= ntp) call error_mesg ('edmf_mynn_driver_mod', &
                    'cldfra_bl can not be a prognostic tracer', FATAL)
  if (nqc_bl <= ntp) call error_mesg ('edmf_mynn_driver_mod', &
                    'qc_bl can not be a prognostic tracer', FATAL)

!-----------------------------------------------------------------------
!--- Register diagnostic fields       
!-----------------------------------------------------------------------

  id_u_flux = register_diag_field (mod_name, 'u_flux', axes(1:2), Time, &
                 'zonal wind stress', 'kg/m/s2' , &
                 missing_value=missing_value )

  id_v_flux = register_diag_field (mod_name, 'v_flux', axes(1:2), Time, &
                 'meridional wind stress', 'kg/m/s2' , &
                 missing_value=missing_value )

  id_u_star_updated = register_diag_field (mod_name, 'u_star_updated', axes(1:2), Time, &
                 'u_star from u_flux and v_flux', 'm/s' , &
                 missing_value=missing_value )

  id_shflx_star = register_diag_field (mod_name, 'shflx_star', axes(1:2), Time, &
                 'sensible heat flux from star', 'W/m2' , &
                 missing_value=missing_value )

  id_lhflx_star = register_diag_field (mod_name, 'lhflx_star', axes(1:2), Time, &
                 'evaporation flux from star', 'kg/m2/s' , &
                 missing_value=missing_value )

  id_w1_thv1_surf_star = register_diag_field (mod_name, 'w1_thv1_surf_star', axes(1:2), Time, &
                 'w1_theta_v1 from star', 'K m/s' , &
                 missing_value=missing_value )

  id_w1_thv1_surf_updated = register_diag_field (mod_name, 'w1_thv1_surf_updated', axes(1:2), Time, &
                 'w1_theta_v1 from updated fluxes', 'K m/s' , &
                 missing_value=missing_value )

  id_Obukhov_length_star = register_diag_field (mod_name, 'Obukhov_length_star', axes(1:2), Time, &
                 'Obukhov length from star', 'm' , &
                 missing_value=missing_value )

  id_Obukhov_length_updated = register_diag_field (mod_name, 'Obukhov_length_updated', axes(1:2), Time, &
                 'Obukhov length from updated fluxes', 'm' , &
                 missing_value=missing_value )

  id_tke_edmf = register_diag_field (mod_name, 'tke_edmf', axes(full), Time, &
                 'turbulent kinetic energy in edmf_mynn', 'm2/s2' , &
                 missing_value=missing_value )

  id_Tsq = register_diag_field (mod_name, 'Tsq', axes(full), Time, &
                 'variance of theta_l', 'K^2' , &
                 missing_value=missing_value )

  id_Cov_thl_qt = register_diag_field (mod_name, 'Cov_thl_qt', axes(full), Time, &
                 'covariance of theta_l and q_t', 'none' , &
                 missing_value=missing_value )

  id_udt_edmf = register_diag_field (mod_name, 'udt_edmf', axes(full), Time, &
                 'u tendency from edmf_mynn', 'm/s2' , &
                 missing_value=missing_value )

  id_vdt_edmf = register_diag_field (mod_name, 'vdt_edmf', axes(full), Time, &
                 'v tendency from edmf_mynn', 'm/s2' , &
                 missing_value=missing_value )

  id_tdt_edmf = register_diag_field (mod_name, 'tdt_edmf', axes(full), Time, &
                 't tendency from edmf_mynn', 'K/s' , &
                 missing_value=missing_value )

  id_qdt_edmf = register_diag_field (mod_name, 'qdt_edmf', axes(full), Time, &
                 'q tendency from edmf_mynn', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qidt_edmf = register_diag_field (mod_name, 'qidt_edmf', axes(full), Time, &
                 'qi tendency from edmf_mynn', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qadt_edmf = register_diag_field (mod_name, 'qadt_edmf', axes(full), Time, &
                 'cldfra tendency from edmf_mynn', '1/s' , &
                 missing_value=missing_value )

  id_qldt_edmf = register_diag_field (mod_name, 'qldt_edmf', axes(full), Time, &
                 'ql tendency from edmf_mynn', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qndt_edmf = register_diag_field (mod_name, 'qndt_edmf', axes(full), Time, &
                 'qn tendency from edmf_mynn', '1/kg/s' , &
                 missing_value=missing_value )

  id_edmf_a = register_diag_field (mod_name, 'edmf_a', axes(half), Time, &
                 'updraft area', 'none' , &
                 missing_value=missing_value )

  id_edmf_w = register_diag_field (mod_name, 'edmf_w', axes(half), Time, &
                 'vertical velocity of updrafts', 'm/s' , &
                 missing_value=missing_value )

  id_edmf_qt = register_diag_field (mod_name, 'edmf_qt', axes(half), Time, &
                 'qt diff between updrafts and grid-mean', 'kg/kg' , &
                 missing_value=missing_value )

  id_edmf_thl = register_diag_field (mod_name, 'edmf_thl', axes(half), Time, &
                 'thl diff between updrafts and grid-mean', 'K' , &
                 missing_value=missing_value )

  id_edmf_ent = register_diag_field (mod_name, 'edmf_ent', axes(full), Time, &
                 'entrainment in updrafts', '1/m' , &
                 missing_value=missing_value )

  id_edmf_det = register_diag_field (mod_name, 'edmf_det', axes(full), Time, &
                 'dentrainment in updrafts', '1/m' , &
                 missing_value=missing_value )

  id_edmf_qc = register_diag_field (mod_name, 'edmf_qc', axes(half), Time, &
                 'qc in updrafts', 'kg/kg' , &
                 missing_value=missing_value )

  id_thl_edmf = register_diag_field (mod_name, 'thl_edmf', axes(full), Time, &
                 'grid-scale theta_li in edmf_mynn', 'K' , &
                 missing_value=missing_value )

  id_qt_edmf = register_diag_field (mod_name, 'qt_edmf', axes(full), Time, &
                 'grid-scale qt in edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_cldfra_bl = register_diag_field (mod_name, 'cldfra_bl', axes(full), Time, &
                 'cloud fraction in edmf_mynn', 'none' , &
                 missing_value=missing_value )

  id_qc_bl = register_diag_field (mod_name, 'qc_bl', axes(full), Time, &
                 'liquid water mixing ratio in edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_z_pbl = register_diag_field (mod_name, 'z_pbl', axes(1:2), Time, &
                 'depth of planetary boundary layer', 'm' , &
                 missing_value=missing_value )

  id_z_pbl_edmf = register_diag_field (mod_name, 'z_pbl_edmf', axes(1:2), Time, &
                 'PBL depth from edmf_mynn', 'm' , &
                 missing_value=missing_value )

  id_qtdt_edmf = register_diag_field (mod_name, 'qtdt_edmf', axes(full), Time, &
                 'qt tendency from edmf_mynn', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_thldt_edmf = register_diag_field (mod_name, 'thldt_edmf', axes(full), Time, &
                 'thl tendency from edmf_mynn', 'K/s' , &
                 missing_value=missing_value )

  id_diff_t_edmf = register_diag_field (mod_name, 'diff_t_edmf', axes(half), Time, &
                 'heat diff coeffs from edmf_mynn', 'K/m/s' , &
                 missing_value=missing_value )

  id_diff_m_edmf = register_diag_field (mod_name, 'diff_m_edmf', axes(half), Time, &
                 'momentum diff coeffs from edmf_mynn', 'm2/s' , &
                 missing_value=missing_value )

  id_el_edmf = register_diag_field (mod_name, 'el_edmf', axes(full), Time, &
                 'mixing length in edmf_mynn', 'm' , &
                 missing_value=missing_value )

  id_t_input = register_diag_field (mod_name, 't_input', axes(full), Time, &
                 'T input to edmf_mynn', 'K' , &
                 missing_value=missing_value )

  id_q_input = register_diag_field (mod_name, 'q_input', axes(full), Time, &
                 'q input to edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_qa_input = register_diag_field (mod_name, 'qa_input', axes(full), Time, &
                 'qa input to edmf_mynn', 'none' , &
                 missing_value=missing_value )

  id_ql_input = register_diag_field (mod_name, 'ql_input', axes(full), Time, &
                 'ql input to edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_qi_input = register_diag_field (mod_name, 'qi_input', axes(full), Time, &
                 'qi input to edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_thl_input = register_diag_field (mod_name, 'thl_input', axes(full), Time, &
                 'thl input to edmf_mynn', 'K' , &
                 missing_value=missing_value )

  id_qt_input = register_diag_field (mod_name, 'qt_input', axes(full), Time, &
                 'qt input to edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_rh_input = register_diag_field (mod_name, 'rh_input', axes(full), Time, &
                 'rh input to edmf_mynn', 'percent' , &
                 missing_value=missing_value )

  id_th_input = register_diag_field (mod_name, 'th_input', axes(full), Time, &
                 'theta input to edmf_mynn', 'K' , &
                 missing_value=missing_value )

  id_t_before_mix = register_diag_field (mod_name, 't_before_mix', axes(full), Time, &
                 'T before_mix to edmf_mynn', 'K' , &
                 missing_value=missing_value )

  id_q_before_mix = register_diag_field (mod_name, 'q_before_mix', axes(full), Time, &
                 'q before_mix to edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_qa_before_mix = register_diag_field (mod_name, 'qa_before_mix', axes(full), Time, &
                 'qa before_mix to edmf_mynn', 'none' , &
                 missing_value=missing_value )

  id_ql_before_mix = register_diag_field (mod_name, 'ql_before_mix', axes(full), Time, &
                 'ql before_mix to edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_qi_before_mix = register_diag_field (mod_name, 'qi_before_mix', axes(full), Time, &
                 'qi before_mix to edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_thl_before_mix = register_diag_field (mod_name, 'thl_before_mix', axes(full), Time, &
                 'thl before_mix to edmf_mynn', 'K' , &
                 missing_value=missing_value )

  id_qt_before_mix = register_diag_field (mod_name, 'qt_before_mix', axes(full), Time, &
                 'qt before_mix to edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_rh_before_mix = register_diag_field (mod_name, 'rh_before_mix', axes(full), Time, &
                 'rh before_mix to edmf_mynn', 'percent' , &
                 missing_value=missing_value )

  id_th_before_mix = register_diag_field (mod_name, 'th_before_mix', axes(full), Time, &
                 'theta before_mix to edmf_mynn', 'K' , &
                 missing_value=missing_value )

  id_t_after_mix = register_diag_field (mod_name, 't_after_mix', axes(full), Time, &
                 'T after_mix to edmf_mynn', 'K' , &
                 missing_value=missing_value )

  id_q_after_mix = register_diag_field (mod_name, 'q_after_mix', axes(full), Time, &
                 'q after_mix to edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_qa_after_mix = register_diag_field (mod_name, 'qa_after_mix', axes(full), Time, &
                 'qa after_mix to edmf_mynn', 'none' , &
                 missing_value=missing_value )

  id_ql_after_mix = register_diag_field (mod_name, 'ql_after_mix', axes(full), Time, &
                 'ql after_mix to edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_qi_after_mix = register_diag_field (mod_name, 'qi_after_mix', axes(full), Time, &
                 'qi after_mix to edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_thl_after_mix = register_diag_field (mod_name, 'thl_after_mix', axes(full), Time, &
                 'thl after_mix to edmf_mynn', 'K' , &
                 missing_value=missing_value )

  id_qt_after_mix = register_diag_field (mod_name, 'qt_after_mix', axes(full), Time, &
                 'qt after_mix to edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_rh_after_mix = register_diag_field (mod_name, 'rh_after_mix', axes(full), Time, &
                 'rh after_mix to edmf_mynn', 'percent' , &
                 missing_value=missing_value )

  id_th_after_mix = register_diag_field (mod_name, 'th_after_mix', axes(full), Time, &
                 'theta after_mix to edmf_mynn', 'K' , &
                 missing_value=missing_value )

  id_qa_before_pdf = register_diag_field (mod_name, 'qa_before_pdf', axes(full), Time, &
                 'qa diagnosed by edmf_mynn', 'none' , &
                 missing_value=missing_value )

  id_ql_before_pdf = register_diag_field (mod_name, 'ql_before_pdf', axes(full), Time, &
                 'ql diagnosed by edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_qi_before_pdf = register_diag_field (mod_name, 'qi_before_pdf', axes(full), Time, &
                 'qi diagnosed by edmf_mynn', 'kg/kg' , &
                 missing_value=missing_value )

  id_tdt_edmf_orig = register_diag_field (mod_name, 'tdt_edmf_orig', axes(full), Time, &
                 't tendency from edmf_mynn original', 'K/s' , &
                 missing_value=missing_value )

  id_qdt_edmf_orig = register_diag_field (mod_name, 'qdt_edmf_orig', axes(full), Time, &
                 'q tendency from edmf_mynn original', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qadt_edmf_orig = register_diag_field (mod_name, 'qadt_edmf_orig', axes(full), Time, &
                 'cldfra tendency from edmf_mynn original', '1/s' , &
                 missing_value=missing_value )

  id_qidt_edmf_orig = register_diag_field (mod_name, 'qidt_edmf_orig', axes(full), Time, &
                 'qi tendency from edmf_mynn original', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qldt_edmf_orig = register_diag_field (mod_name, 'qldt_edmf_orig', axes(full), Time, &
                 'ql tendency from edmf_mynn original', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qldt_edmf_ED = register_diag_field (mod_name, 'qldt_edmf_ED', axes(full), Time, &
                 'ql tendency from edmf_mynn, ED', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qldt_edmf_MF = register_diag_field (mod_name, 'qldt_edmf_MF', axes(full), Time, &
                 'ql tendency from edmf_mynn, MF', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qidt_edmf_ED = register_diag_field (mod_name, 'qidt_edmf_ED', axes(full), Time, &
                 'qi tendency from edmf_mynn, ED', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qidt_edmf_MF = register_diag_field (mod_name, 'qidt_edmf_MF', axes(full), Time, &
                 'qi tendency from edmf_mynn, MF', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qadt_edmf_ED = register_diag_field (mod_name, 'qadt_edmf_ED', axes(full), Time, &
                 'qa tendency from edmf_mynn, ED', '1/s' , &
                 missing_value=missing_value )

  id_qadt_edmf_MF = register_diag_field (mod_name, 'qadt_edmf_MF', axes(full), Time, &
                 'qa tendency from edmf_mynn, MF', '1/s' , &
                 missing_value=missing_value )

  id_qldt_edmf_MF_adv = register_diag_field (mod_name, 'qldt_edmf_MF_adv', axes(full), Time, &
                 'ql tendency from edmf_mynn, MF_adv', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qldt_edmf_MF_eddy = register_diag_field (mod_name, 'qldt_edmf_MF_eddy', axes(full), Time, &
                 'ql tendency from edmf_mynn, MF_eddy', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qidt_edmf_MF_adv = register_diag_field (mod_name, 'qidt_edmf_MF_adv', axes(full), Time, &
                 'qi tendency from edmf_mynn, MF_adv', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qidt_edmf_MF_eddy = register_diag_field (mod_name, 'qidt_edmf_MF_eddy', axes(full), Time, &
                 'qi tendency from edmf_mynn, MF_eddy', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qadt_edmf_MF_adv = register_diag_field (mod_name, 'qadt_edmf_MF_adv', axes(full), Time, &
                 'qa tendency from edmf_mynn, MF_adv', '1/s' , &
                 missing_value=missing_value )

  id_qadt_edmf_MF_eddy = register_diag_field (mod_name, 'qadt_edmf_MF_eddy', axes(full), Time, &
                 'qa tendency from edmf_mynn, MF_eddy', '1/s' , &
                 missing_value=missing_value )

  id_qldt_edmf_MF_ent = register_diag_field (mod_name, 'qldt_edmf_MF_ent', axes(full), Time, &
                 'ql tendency from edmf_mynn, MF_ent', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qidt_edmf_MF_ent = register_diag_field (mod_name, 'qidt_edmf_MF_ent', axes(full), Time, &
                 'qi tendency from edmf_mynn, MF_ent', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qadt_edmf_MF_ent = register_diag_field (mod_name, 'qadt_edmf_MF_ent', axes(full), Time, &
                 'qa tendency from edmf_mynn, MF_ent', '1/s' , &
                 missing_value=missing_value )

  id_qadt_edmf_MF_det = register_diag_field (mod_name, 'qadt_edmf_MF_det', axes(full), Time, &
                 'qa tendency from edmf_mynn, MF_det', '1/s' , &
                 missing_value=missing_value )

  id_qldt_edmf_MF_det = register_diag_field (mod_name, 'qldt_edmf_MF_det', axes(full), Time, &
                 'ql tendency from edmf_mynn, MF_det', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qidt_edmf_MF_det = register_diag_field (mod_name, 'qidt_edmf_MF_det', axes(full), Time, &
                 'qi tendency from edmf_mynn, MF_det', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qadt_edmf_MF_sub = register_diag_field (mod_name, 'qadt_edmf_MF_sub', axes(full), Time, &
                 'qa tendency from edmf_mynn, MF_sub', '1/s' , &
                 missing_value=missing_value )

  id_qldt_edmf_MF_sub = register_diag_field (mod_name, 'qldt_edmf_MF_sub', axes(full), Time, &
                 'ql tendency from edmf_mynn, MF_sub', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_qidt_edmf_MF_sub = register_diag_field (mod_name, 'qidt_edmf_MF_sub', axes(full), Time, &
                 'qi tendency from edmf_mynn, MF_sub', 'kg/kg/s' , &
                 missing_value=missing_value )

  id_a_moist_half = register_diag_field (mod_name, 'a_moist_half', axes(half), Time, &
                 'moist updraft area on phalf', 'none' , &
                 missing_value=missing_value )

  id_a_moist_full = register_diag_field (mod_name, 'a_moist_full', axes(full), Time, &
                 'moist updraft area on pfull', 'none' , &
                 missing_value=missing_value )

  id_mf_moist_half = register_diag_field (mod_name, 'mf_moist_half', axes(half), Time, &
                 'moist updraft mass flux on phalf', 'kg/m2/s' , &
                 missing_value=missing_value )

  id_mf_moist_full = register_diag_field (mod_name, 'mf_moist_full', axes(full), Time, &
                 'moist updraft mass flux on pfull', 'kg/m2/s' , &
                 missing_value=missing_value )

  id_qv_moist_half = register_diag_field (mod_name, 'qv_moist_half', axes(half), Time, &
                 'spec humid of moist updraft on phalf', 'kg/kg' , &
                 missing_value=missing_value )

  id_qv_moist_full = register_diag_field (mod_name, 'qv_moist_full', axes(full), Time, &
                 'spec humid of moist updraft on pfull', 'kg/kg' , &
                 missing_value=missing_value )

  id_a_dry_half = register_diag_field (mod_name, 'a_dry_half', axes(half), Time, &
                 'dry updraft area on phalf', 'none' , &
                 missing_value=missing_value )

  id_a_dry_full = register_diag_field (mod_name, 'a_dry_full', axes(full), Time, &
                 'dry updraft area on pfull', 'none' , &
                 missing_value=missing_value )

  id_mf_dry_half = register_diag_field (mod_name, 'mf_dry_half', axes(half), Time, &
                 'dry updraft mass flux on phalf', 'kg/m2/s' , &
                 missing_value=missing_value )

  id_mf_dry_full = register_diag_field (mod_name, 'mf_dry_full', axes(full), Time, &
                 'dry updraft mass flux on pfull', 'kg/m2/s' , &
                 missing_value=missing_value )

  id_qv_dry_half = register_diag_field (mod_name, 'qv_dry_half', axes(half), Time, &
                 'spec humid of dry updraft on phalf', 'kg/kg' , &
                 missing_value=missing_value )

  id_qv_dry_full = register_diag_field (mod_name, 'qv_dry_full', axes(full), Time, &
                 'spec humid of dry updraft on pfull', 'kg/kg' , &
                 missing_value=missing_value )

  id_mf_all_half = register_diag_field (mod_name, 'mf_all_half', axes(half), Time, &
                 'updraft mass flux on phalf', 'kg/m2/s' , &
                 missing_value=missing_value )

  id_mf_all_full = register_diag_field (mod_name, 'mf_all_full', axes(full), Time, &
                 'updraft mass flux on pfull', 'kg/m2/s' , &
                 missing_value=missing_value )

  id_num_updraft = register_diag_field (mod_name, 'num_updraft', axes(half), Time, &
                 'number of edmf updrafts', 'none' , &
                 missing_value=missing_value )

  id_num_det = register_diag_field (mod_name, 'num_det', axes(full), Time, &
                 'number of dentrainment in edmf updrafts', 'none' , &
                 missing_value=missing_value )

  id_num_ndet_zent = register_diag_field (mod_name, 'num_ndet_zent', axes(full), Time, &
                 'frequency of negative dentrainment with zero entrainment in edmf updrafts', 'none' , &
                 missing_value=missing_value )

  id_num_ndet_pent = register_diag_field (mod_name, 'num_ndet_pent', axes(full), Time, &
                 'frequency of negative dentrainment with positive entrainment in edmf updrafts', 'none' , &
                 missing_value=missing_value )

  id_rlz_ratio = register_diag_field (mod_name, 'rlz_ratio', axes(1:2), Time, &
                 'ratio for realizability limiter', 'none' , &
                 missing_value=missing_value )

  id_rlz_tracer = register_diag_field (mod_name, 'rlz_tracer', axes(1:2), Time, &
                 'tracer name for realizability limiter (qv,l,i,a,t)', 'none' , &
                 missing_value=missing_value )

!-----------------------------------------------------------------------
!--- allocate edmf2ls_mp variables  
!-----------------------------------------------------------------------
!  allocate(edmf2ls_mp%option_edmf2ls_mp)       ; edmf2ls_mp%option_edmf2ls_mp = .false.
!
!  allocate(edmf2ls_mp%qadt_edmf(id,jd,kd)) ; edmf2ls_mp%qadt_edmf = 0. 
!  allocate(edmf2ls_mp%qldt_edmf(id,jd,kd)) ; edmf2ls_mp%qldt_edmf = 0. 
!  allocate(edmf2ls_mp%qidt_edmf(id,jd,kd)) ; edmf2ls_mp%qidt_edmf = 0. 
!
!  allocate(edmf2ls_mp%dqa_edmf(id,jd,kd)) ; edmf2ls_mp%dqa_edmf = 0. 
!  allocate(edmf2ls_mp%dql_edmf(id,jd,kd)) ; edmf2ls_mp%dql_edmf = 0. 
!  allocate(edmf2ls_mp%dqi_edmf(id,jd,kd)) ; edmf2ls_mp%dqi_edmf = 0. 

!-----------------------------------------------------------------------
!--- Done with initialization
!-----------------------------------------------------------------------

  module_is_initialized = .true.

!-----------------------------------------------------------------------

end subroutine edmf_mynn_init

!#######################################################################

subroutine edmf_mynn_end

!-----------------------------------------------------------------------

  if (.not.module_is_initialized) return

  module_is_initialized = .false.

!-----------------------------------------------------------------------

end subroutine edmf_mynn_end

!#######################################################################

!--------------------------------
! DESCRIPTION OF THE MYNN-PBL module
!
! module_bl_mynn
!  (1) sets key constants
!  (2) initializes variables in  mynn_bl_init_driver - do this only once at the start of simulation
!  
!  
!  (3) during time integration runs mynn_bl_driver which is the EDMF parameterization build on the WRF MYNN ED (v.4.0)
!
!
! The mynn_bl_driver should be called within the main program with the following form
!  (when compiling do not invoke WRF_CHEM in the preprocessor)
!
!   mynn_bl_driver(            &
!       &initflag,grav_settling,         &
!       &delt,dz,dx,znt,                 &
!       &u,v,w,th,qv,qc,qi,qni,qnc,      &
!       &p,exner,rho,T3D,                &
!       &xland,ts,qsfc,qcg,ps,           &
!       &ust,ch,hfx,qfx,rmol,wspd,       &
!       &uoce,voce,                      & 
!       &vdfg,                           & 
!       &Qke,tke_pbl,                    &
!       &qke_adv,bl_mynn_tkeadvect,      & 
!       &Tsq,Qsq,Cov,                    &
!       &RUBLTEN,RVBLTEN,RTHBLTEN,       &
!       &RQVBLTEN,RQCBLTEN,RQIBLTEN,     &
!       &RQNIBLTEN,                      &
!       &exch_h,exch_m,                  &
!       &Pblh,kpbl,                      & 
!       &el_pbl,                         &
!       &dqke,qWT,qSHEAR,qBUOY,qDISS,    & 
!       &wstar,delta,                    & 
!       &bl_mynn_tkebudget,              &
!       &bl_mynn_cloudpdf,Sh3D,          &
!       &bl_mynn_mixlength,              &
!       &icloud_bl,qc_bl,cldfra_bl,      &
!       &bl_mynn_edmf,                   &
!       &bl_mynn_edmf_dd,                   &
!       &bl_mynn_edmf_mom,bl_mynn_edmf_tke, &
!       &bl_mynn_edmf_part,bl_mynn_edmf_Lent,&
!       &bl_mynn_cloudmix,bl_mynn_mixqt, &
!       &edmf_a,edmf_w,edmf_qt,          &
!       &edmf_thl,edmf_ent,edmf_qc,      &
!       &edmf_debug1,edmf_debug2,        &
!       &edmf_debug3,edmf_debug4,        &
!       &edmf_a_dd,edmf_w_dd,edmf_qt_dd,    &
!       &edmf_thl_dd,edmf_ent_dd,edmf_qc_dd,&
!       &mynn_ql,                        &
!       &nupdraft,maxMF,ktop_shallow,    &
!       &spp_pbl,pattern_spp_pbl,        &
!       &RTHRATEN,                       &
!       &FLAG_QI,FLAG_QNI,FLAG_QC,FLAG_QNC &
!       &,IDS,IDE,JDS,JDE,KDS,KDE        &
!       &,IMS,IME,JMS,JME,KMS,KME        &
!       &,ITS,ITE,JTS,JTE,KTS,KTE)
!    
!
! INPUTS/FLAGS
!
! INTEGER :: level <>3 set in "mynn_bl_init_driver"
!             (to do the level 2.5 MYNN model; if 3 then will do level 3)
!  INTEGER :: initflag .... set to 1 at the start of simulations and 0 otherwise 
!                          (when 1 it initializes TKE using level 2 MYNN scheme)
!  INTEGER :: grav_settling ...  0 
!                         (1 - the cloud/fog will experience gravitational settling, 0 - it does not)
!    INTEGER::   bl_mynn_tkebudget .... 0/1
!                          (if 1 the budget terms in the TKE equation are allocated for output, if 0 then not)
!    INTEGER::  bl_mynn_cloudpdf .... 1
!                         (define the type of the subgrid PDF for cloud computation, 1= Nakanishi & Niino 2004) 
!    INTEGER  :: bl_mynn_mixlength ... 2 
!                        (defines the ED mixing length formulation)
!    INTEGER :: bl_mynn_edmf  ... 3 
!                          (controls  the version of the EDMF to be called 3=JPLedmf documented by Wu et al., 2020)
!    INTEGER :: bl_mynn_edmf_dd  .... 0
!                            (0 - no downdrafts, 1 - Wu et al., 2020 downdrafts)
!    REAL  :: bl_mynn_edmf_Lent .... 0.
!                               (dummy argument)
!    LOGICAL  :: bl_mynn_tkeadvect ... 0
!                                (1 - advecting TKE, 0 - not advecting )
!                                (if 1 then qke is stored in qke_adv)
!    INTEGER :: bl_mynn_edmf_mom .... 1 
!                              (1- mixing of momentum by full EDMF, 0 - momenum is mixed only by ED)
!    INTEGER :: bl_mynn_edmf_tke ... 0 
!                       (does not effect JPLedmf)
!    INTEGER :: bl_mynn_edmf_part .... 1
!                       (partition area into updrafts and remaining environment)    
!    INTEGER :: bl_mynn_cloudmix ..... 1
!                               (1 - cloud species for which FLAG_QI=1 are mixed separately, 0 - not)     
!    INTEGER :: bl_mynn_mixqt .... 2 
!                          (will mix moist conserved variables, after mixing invokes PDF cloud scheme to convert moist variables to dry)
!    INTEGER :: icloud_bl .... 1
!         (if 1, cloud cover and liquid water from the PDF scheme will be on the cloud output)       
!
!    LOGICAL :: FLAG_QI .... true
!               FLAG_QNI,FLAG_QC,FLAG_QNC ... all false
!        (flags for whether cloud and ice mixing rations and number concentrations are mixed separately)
!    INTEGER   :: spp_pbl .... 0 
!   (1 stochastic perturbation to condensation)
! 
! INPUTS/DOMAIN DEFINITION
!  INTEGER :: IDS,IDE,JDS,JDE,KDS,KDE,IMS,IME,JMS,JME,KMS,KME,ITS,ITE,JTS,JTE,KTS,KTE
!  ...  I,J,K stand for the three spacial directions (K is vertical)
!  .... D for domain, M ... memory, T ... tile
!  .... S start, D end
!
!
! DEFINITION OF MODEL LEVELS
!   KMS (smallest number) is the bottom level and KME (largest
! number) is the top level.  
! ---------
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
! ----------
!

! INPUTS/0D fields
!
! REAL:: 
! delt ... model time step [s]
! dx ...  horizontal length of grid [m] ... this is used in the scale aware option
!
! INPUTS/2D surface fields
!
!  REAL, DIMENSION(IMS:IME,JMS:JME) :: 
! xland ....   land mask (1 for land, 2 for water)
! ust ...  u* in similarity theory (m/s)
! ch .... drag coefficient for heat/moisture       
! rmol ....   inverse Monin-Obukhov length (1/m)
! ts  ... surface temperature (K) 
! qsfc ...   specific humidity at lower boundary (kg/kg)
! qcg   ... saturated mixing ratio at the surface (kg/kg)
! ps .... surface pressure [Pa]
! hfx ... upward heat flux at the surface (W/m^2)
! qfx  ... upward moisture flux at the surface (kg/m^2/s)
! wspd ... wind speed at lowest model level (m/s)
! uoce  ...  sea surface zonal ocean currents (m s-1)
! voce ... sea surface meridional ocean currents (m s-1)
! vdfg  ... deposition velocity of fog (m/s); set to zero
! znt ...   thermal roughness length (m)
!
!
!
! INPUTS/3D thermodynamic fields on half levels
!
!   REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(in) :: 
! dz ... dz between full levels [m], on full levels
! u  ... zonal wind speed [m s^-1], on half levels
! v ... meridional wind speed [m s^-1], on half levels
! w ... vertical wind speed [m s^-1], on full levels
! th  ... potential temperature [K], on half levels
! qv .... water vapor mixing ratio [kg kg^-1], half level
! p .... pressure [Pa], half levels
! exner .... exner function, half levels 
! rho ... density [kg m^-3], half levels
! T3D ... temperature [K], half levels
! qc ... cloud water mixing ratio [kg kg^-1], half levels 
! qi ... ice water mixing ratio [kg kg^-1], half levels 
! qnc,qni ....  cloud liq and ice number concentration (#/kg), optional, do not set them
!
!  
!
!INPUT-OUTPUT/2D fields
!
!   DIMENSION( ims:ime, kms:kme, jms:jme ),OPTIONAL  ::pattern_spp_pbl
!         (pattern for stochastic condensation .... not needed if spp_pbl == 0) 
!
!
!    REAL, DIMENSION(IMS:IME,JMS:JME) 
!         pblh,wstar,delta,tke_pbl 
! (these are essentially output variables, for coupling with one of the convection scheme, you can treat them as dummy) 
!
!
!
!INPUT-OUTPUT/3D fields
!
!    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(inout) :: &
!         prognostic variable:   Qke .... twice the turbulent kinetic energy [m^2 s^-2]
!         dummy for level 2.5:   Tsq  ...variance of theta_l (K^2)   
!         dummy for level 2.5:   Qsq ... variance of q_t (-) 
!         dummy for level 2.5: Cov ... covariance of theta_l and q_t  (K) 
!         dummy for tke_advect=0:  qke_adv   .... twice the turbulent kinetic energy after advected with resolved scale flow (advection is done outside the routine) [m^2 s^-2] 
!
!
! INPUT-OUTPUT/3D fields (tendencies, should be set to 0 on input; on full levels)
!    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: &
!         RUBLTEN ... U tendency due to EDMF parameterization (m/s^2), set to 0. on the input
!         RVBLTEN ... V tendency due to EDMF parameterization (m/s^2), set to 0. on the input
!         RTHBLTEN ... Theta tendency due to EDMF parameterization (K/s), set to 0. on the input
!         RQVBLTEN ....  Qv tendency due to EDMF parameterization (kg/kg/s), set to 0. on the input
!         RQCBLTEN ....  Qc tendency due to EDMF parameterization (kg/kg/s), set to 0. on the input
!         RQIBLTEN ...   Qi tendency due to  EDMF parameterization (kg/kg/s), set to 0. on the input
!         RQNIBLTEN ....  Qni tendency due to EDMF parameterization (#/kg/s) , dummy argument (not activated)
!         input:  RTHRATEN .... tendency from radiation [K s^-1]
!
! INPUT/OUTPUT, OPTIONAL/3D fields (updraft/downdraft properties) on full levels
!   REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME):: 
!         edmf_a ... updraft area [-] 
!         edmf_w ... vertical velocity of updrafts [m s^-1]
!         edmf_qt .... qt diff between updrafts and grid-mean  [kg kg^-1]
!         edmf_thl .... thl diff between updrafts and grid-mean  [K]
!         edmf_ent ... entrainment in updrafts [m^-1]
!         edmf_qc ... qc in updrafts [kg kg^-1]
!         edmf_a_dd,edmf_w_dd,edmf_qt_dd,edmf_thl_dd,edmf_ent_dd,edmf_qc_dd ... dummy outputs (downdraft poperties)
!         edmf_debug1,edmf_debug2,edmf_debug3,edmf_debug4 .... dummy outputs (for debugging)
!
!
! INPUT/OUTPUT, 3D cloud and turbulence fields
!
! REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME)
!         mynn_ql & qc_bl ... liquid water mixing ratio from PDF cloud scheme [kg kg-1]
!         cldfra_bl ... cloud fraction [-]
!      semi prognostic:   el_pbl ..... turbulent mixing length [m]
!
!
!  INPUT/OUTPUT 
!    INTEGER,DIMENSION(IMS:IME,JMS:JME) 
!       KPBL .... model level of boundary layer height
!       nupdraft ... diagnostic output, number of updrafts (this is only for stem MF and not for JPLedmf)
!     ktop_shallow ... diagnostic output essentially, the level of the highest updraft

   
!
! OUTPUTS, 3D mixing coefficients on full levels
!    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: 
!         exch_h ... ED mixing coefficient for heat, K m/s 
!         exch_m .... ED mixing coefficient for momentum, m^2/s
!
! OUTPUT, 2D
! REAL, DIMENSION(IMS:IME,JMS:JME) ::
!        maxmf   .... maximum mass-flux (does not work with the JPL-EDMF, set to dummy variable)
!
!
! OUTPUTS, 3D fields, budget terms for tke on half level  
! (bl_mynn_tkebudget must be 1, otherwise they are not allocated)
!  ... these fields are for diagnostic purposes only ....
!    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME) ::
!         qWT ... vertical transport [m^2 s^-3]
!         qSHEAR ... shear production [m^2 s^-3]
!         qBUOY  ... buoyancy production [m^2 s^-3]
!         qDISS .... dissipation [m^2 s^-3]
!         dqke ..... delta TKE (within model time step) [m^2 s^-2]
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!WRF:MODEL_LAYER:PHYSICS
!
! translated from NN f77 to F90 and put into WRF by Mariusz Pagowski
! NOAA/GSD & CIRA/CSU, Feb 2008
! changes to original code:
! 1. code is 1D (in z)
! 2. no advection of TKE, covariances and variances 
! 3. Cranck-Nicholson replaced with the implicit scheme
! 4. removed terrain dependent grid since input in WRF in actual
!    distances in z[m]
! 5. cosmetic changes to adhere to WRF standard (remove common blocks, 
!            intent etc)
!-------------------------------------------------------------------
!Modifications primarily from Joseph Olson and Jaymes Kenyon NOAA/GSD/MDB - CU/CIRES
!
! Departures from original MYNN (Nakanish & Niino 2009)
! 1. Addition of BouLac mixing length in the free atmosphere.
! 2. Changed the turbulent mixing length to be integrated from the
!    surface to the top of the BL + a transition layer depth.
! v3.4.1:    Option to use Kitamura/Canuto modification which removes 
!            the critical Richardson number and negative TKE (default).
!            Hybrid PBL height diagnostic, which blends a theta-v-based
!            definition in neutral/convective BL and a TKE-based definition
!            in stable conditions.
!            TKE budget output option (bl_mynn_tkebudget)
! v3.5.0:    TKE advection option (bl_mynn_tkeadvect)
! v3.5.1:    Fog deposition related changes.
! v3.6.0:    Removed fog deposition from the calculation of tendencies
!            Added mixing of qc, qi, qni
!            Added output for wstar, delta, TKE_PBL, & KPBL for correct 
!                   coupling to shcu schemes  
! v3.8.0:    Added subgrid scale cloud output for coupling to radiation
!            schemes (activated by setting icloud_bl =1 in phys namelist).
!            Added WRF_DEBUG prints (at level 3000)
!            Added Tripoli and Cotton (1981) correction.
!            Added namelist option bl_mynn_cloudmix to test effect of mixing
!                cloud species (default = 1: on). 
!            Added mass-flux option (bl_mynn_edmf, = 1 for StEM, 2 for TEMF).
!                This option is off by default (=0).
!                Related (hidden) options: 
!                 bl_mynn_edmf_mom = 1 : activate momentum transport in MF scheme
!                 bl_mynn_edmf_tke = 1 : activate TKE transport in MF scheme
!                 bl_mynn_edmf_part= 1 : activate areal partitioning of ED & MF
!            Added mixing length option (bl_mynn_mixlength, see notes below)
!            Added more sophisticated saturation checks, following Thompson scheme
!            Added new cloud PDF option (bl_mynn_cloudpdf = 2) from Chaboureau
!                and Bechtold (2002, JAS, with mods) 
!            Added capability to mix chemical species when env variable
!                WRF_CHEM = 1, thanks to Wayne Angevine.
!            Added scale-aware mixing length, following Junshi Ito's work
!                Ito et al. (2015, BLM).
! v3.9.0    Improvement to the mass-flux scheme (dynamic number of plumes,
!                better plume/cloud depth, significant speed up, better cloud
!                fraction). 
!            Added Stochastic Parameter Perturbation (SPP) implementation.
!            Many miscellaneous tweaks to the mixing lengths and stratus
!                component of the subgrid clouds.
! v.4.0      Removed or added alternatives for WRF-specific functions/modules
!                the sake of portability to other models.
!            Further refinement of mass-flux scheme from SCM experiments with
!                Wayne Angevine: switch to linear entrainment and back to
!                Simpson and Wiggert-type w-equation.
!            Addition of TKE production due to radiation cooling at top of 
!                clouds (proto-version); not activated by default.
!            Some code rewrites to move if-thens out of loops in an attempt to
!                improve computational efficiency.
!            New tridiagonal solver, which is supposedly 14% faster and more
!                conservative. Impact seems very small.
!            Many miscellaneous tweaks to the mixing lengths and stratus
!                component of the subgrid clouds.
! 
!-------------------------------------------------------------------

!For WRF:
!  USE module_model_constants, only: &
!       &karman, g, p1000mb, &
!       &cp, r_d, r_v, rcp, xlv, xlf, xls, &
!       &svp1, svp2, svp3, svpt0, ep_1, ep_2, rvovrd, &
!       &cpv, cliq, cice, rhoair0
!
!  USE module_state_description, only: param_first_scalar, &
!       &p_qc, p_qr, p_qi, p_qs, p_qg, p_qnc, p_qni 

! **********************************************************************
! *   An improved Mellor-Yamada turbulence closure model               *
! *                                                                    *
! *                                   Aug/2005  M. Nakanishi (N.D.A)   *
! *                        Modified:  Dec/2005  M. Nakanishi (N.D.A)   *
! *                                             naka@nda.ac.jp         *
! *                                                                    *
! *   Contents:                                                        *
! *     1. mym_initialize  (to be called once initially)               *
! *        gives the closure constants and initializes the turbulent   *
! *        quantities.                                                 *
! *    (2) mym_level2      (called in the other subroutines)           *
! *        calculates the stability functions at Level 2.              *
! *    (3) mym_length      (called in the other subroutines)           *
! *        calculates the master length scale.                         *
! *     4. mym_turbulence                                              *
! *        calculates the vertical diffusivity coefficients and the    *
! *        production terms for the turbulent quantities.              *
! *     5. mym_predict                                                 *
! *        predicts the turbulent quantities at the next step.         *
! *     6. mym_condensation                                            *
! *        determines the liquid water content and the cloud fraction  *
! *        diagnostically.                                             *
! *                                                                    *
! *             call mym_initialize                                    *
! *                  |                                                 *
! *                  |<----------------+                               *
! *                  |                 |                               *
! *             call mym_condensation  |                               *
! *             call mym_turbulence    |                               *
! *             call mym_predict       |                               *
! *                  |                 |                               *
! *                  |-----------------+                               *
! *                  |                                                 *
! *                 end                                                *
! *                                                                    *
! *   Variables worthy of special mention:                             *
! *     tref   : Reference temperature                                 *
! *     thl     : Liquid water potential temperature                   *
! *     qw     : Total water (water vapor+liquid water) content        *
! *     ql     : Liquid water content                                  *
! *     vt, vq : Functions for computing the buoyancy flux             *
! *                                                                    *
! *     If the water contents are unnecessary, e.g., in the case of    *
! *     ocean models, thl is the potential temperature and qw, ql, vt  *
! *     and vq are all zero.                                           *
! *                                                                    *
! *   Grid arrangement:                                                *
! *             k+1 +---------+                                        *
! *                 |         |     i = 1 - nx                         *
! *             (k) |    *    |     j = 1 - ny                         *
! *                 |         |     k = 1 - nz                         *
! *              k  +---------+                                        *
! *                 i   (i)  i+1                                       *
! *                                                                    *
! *     All the predicted variables are defined at the center (*) of   *
! *     the grid boxes. The diffusivity coefficients are, however,     *
! *     defined on the walls of the grid boxes.                        *
! *     # Upper boundary values are given at k=nz.                     *
! *                                                                    *
! *   References:                                                      *
! *     1. Nakanishi, M., 2001:                                        *
! *        Boundary-Layer Meteor., 99, 349-378.                        *
! *     2. Nakanishi, M. and H. Niino, 2004:                           *
! *        Boundary-Layer Meteor., 112, 1-31.                          *
! *     3. Nakanishi, M. and H. Niino, 2006:                           *
! *        Boundary-Layer Meteor., (in press).                         *
! *     4. Nakanishi, M. and H. Niino, 2009:                           *
! *        Jour. Meteor. Soc. Japan, 87, 895-912.                      *
! **********************************************************************
!
!     SUBROUTINE  mym_initialize:
!
!     Input variables:
!       iniflag         : <>0; turbulent quantities will be initialized
!                         = 0; turbulent quantities have been already
!                              given, i.e., they will not be initialized
!       nx, ny, nz      : Dimension sizes of the
!                         x, y and z directions, respectively
!       tref            : Reference temperature                      (K)
!       dz(nz)        : Vertical grid spacings                     (m)
!                         # dz(nz)=dz(nz-1)
!       zw(nz+1)        : Heights of the walls of the grid boxes     (m)
!                         # zw(1)=0.0 and zw(k)=zw(k-1)+dz(k-1)
!       h(nx,ny)        : G^(1/2) in the terrain-following coordinate
!                         # h=1-zg/zt, where zg is the height of the
!                           terrain and zt the top of the model domain
!       pi0(nx,my,nz) : Exner function at zw*h+zg             (J/kg K)
!                         defined by c_p*( p_basic/1000hPa )^kappa
!                         This is usually computed by integrating
!                         d(pi0)/dz = -h*g/tref.
!       rmo(nx,ny)      : Inverse of the Obukhov length         (m^(-1))
!       flt, flq(nx,ny) : Turbulent fluxes of sensible and latent heat,
!                         respectively, e.g., flt=-u_*Theta_* (K m/s)
!! flt - liquid water potential temperature surface flux
!! flq - total water flux surface flux
!       ust(nx,ny)      : Friction velocity                        (m/s)
!       pmz(nx,ny)      : phi_m-zeta at z1*h+z0, where z1 (=0.5*dz(1))
!                         is the first grid point above the surafce, z0
!                         the roughness length and zeta=(z1*h+z0)*rmo
!       phh(nx,ny)      : phi_h at z1*h+z0
!       u, v(nx,nz,ny): Components of the horizontal wind        (m/s)
!       thl(nx,nz,ny)  : Liquid water potential temperature
!                                                                    (K)
!       qw(nx,nz,ny)  : Total water content Q_w                (kg/kg)
!
!     Output variables:
!       ql(nx,nz,ny)  : Liquid water content                   (kg/kg)
!       v?(nx,nz,ny)  : Functions for computing the buoyancy flux
!       qke(nx,nz,ny) : Twice the turbulent kinetic energy q^2
!                                                              (m^2/s^2)
!       tsq(nx,nz,ny) : Variance of Theta_l                      (K^2)
!       qsq(nx,nz,ny) : Variance of Q_w
!       cov(nx,nz,ny) : Covariance of Theta_l and Q_w              (K)
!       el(nx,nz,ny)  : Master length scale L                      (m)
!                         defined on the walls of the grid boxes
!
!     Work arrays:        see subroutine mym_level2
!       pd?(nx,nz,ny) : Half of the production terms at Level 2
!                         defined on the walls of the grid boxes
!       qkw(nx,nz,ny) : q on the walls of the grid boxes         (m/s)
!
!     # As to dtl, ...gh, see subroutine mym_turbulence.
!
!-------------------------------------------------------------------
  SUBROUTINE  mym_initialize (                                & 
       &            kts,kte,                                  &
       &            dz, zw,                                   &
       &            u, v, thl, qw,                            &
!       &            ust, rmo, pmz, phh, flt, flq,             &
       &            zi, theta, sh,                            &
       &            ust, rmo, el,                             &
       &            Qke, Tsq, Qsq, Cov, Psig_bl, cldfra_bl1D, &
       &            bl_mynn_mixlength,                        &
       &            edmf_w1,edmf_a1,edmf_qc1,bl_mynn_edmf,    &
       &            edmf_w_dd1,edmf_a_dd1,edmf_qc_dd1,bl_mynn_edmf_dd,&
       &            rstoch_col,kpbl)
!
!-------------------------------------------------------------------
    
    INTEGER, INTENT(IN)   :: kts,kte,kpbl
    INTEGER, INTENT(IN)   :: bl_mynn_mixlength,bl_mynn_edmf,bl_mynn_edmf_dd
!    REAL, INTENT(IN)   :: ust, rmo, pmz, phh, flt, flq
    REAL, INTENT(IN)   :: ust, rmo, Psig_bl
    REAL, DIMENSION(kts:kte), INTENT(in) :: dz
    REAL, DIMENSION(kts:kte+1), INTENT(in) :: zw
    REAL, DIMENSION(kts:kte), INTENT(in) :: u,v,thl,qw,cldfra_bl1D,&
                                          edmf_w1,edmf_a1,edmf_qc1,&
                                          edmf_w_dd1,edmf_a_dd1,edmf_qc_dd1
    REAL, DIMENSION(kts:kte), INTENT(out) :: tsq,qsq,cov
    REAL, DIMENSION(kts:kte), INTENT(inout) :: el,qke

    REAL, DIMENSION(kts:kte) :: &
         &ql,pdk,pdt,pdq,pdc,dtl,dqw,dtv,&
         &gm,gh,sm,sh,qkw,vt,vq
    INTEGER :: k,l,lmax
    REAL :: phm,vkz,elq,elv,b1l,b2l,pmz=1.,phh=1.,flt=0.,flq=0.,tmpq
    REAL :: zi
    REAL, DIMENSION(kts:kte) :: theta

    REAL, DIMENSION(kts:kte) :: rstoch_col

!   **  At first ql, vt and vq are set to zero.  **
    DO k = kts,kte
       ql(k) = 0.0
       vt(k) = 0.0
       vq(k) = 0.0
    END DO
!
    CALL mym_level2 ( kts,kte,&
         &            dz,  &
         &            u, v, thl, qw, &
         &            ql, vt, vq, &
         &            dtl, dqw, dtv, gm, gh, sm, sh, kpbl )
!
!   **  Preliminary setting  **

    el (kts) = 0.0
    qke(kts) = ust**2 * ( b1*pmz )**(2.0/3.0)
!
    phm      = phh*b2 / ( b1*pmz )**(1.0/3.0)
    tsq(kts) = phm*( flt/ust )**2
    qsq(kts) = phm*( flq/ust )**2
    cov(kts) = phm*( flt/ust )*( flq/ust )
!
    DO k = kts+1,kte
       vkz = vk*zw(k)
       el (k) = vkz/( 1.0 + vkz/100.0 )
       qke(k) = 0.0
!
       tsq(k) = 0.0
       qsq(k) = 0.0
       cov(k) = 0.0
    END DO
!
!   **  Initialization with an iterative manner          **
!   **  lmax is the iteration count. This is arbitrary.  **
    lmax = 5 
!
    DO l = 1,lmax
!
       CALL mym_length (                     &
            &            kts,kte,            &
            &            dz, zw,             &
            &            rmo, flt, flq,      &
            &            vt, vq,             &
            &            qke,                &
            &            dtv,                &
            &            el,                 &
            &            zi,theta,           &
            &            qkw,Psig_bl,cldfra_bl1D,bl_mynn_mixlength,&
            &            edmf_w1,edmf_a1,edmf_qc1,bl_mynn_edmf,&
            &            edmf_w_dd1,edmf_a_dd1,edmf_qc_dd1,bl_mynn_edmf_dd,&
            &            rstoch_col)
!
       DO k = kts+1,kte
          elq = el(k)*qkw(k)
          pdk(k) = elq*( sm(k)*gm (k)+&
               &sh(k)*gh (k) )
          pdt(k) = elq*  sh(k)*dtl(k)**2
          pdq(k) = elq*  sh(k)*dqw(k)**2
          pdc(k) = elq*  sh(k)*dtl(k)*dqw(k)
       END DO
!
!   **  Strictly, vkz*h(i,j) -> vk*( 0.5*dz(1)*h(i,j)+z0 )  **
       vkz = vk*0.5*dz(kts)
!
       elv = 0.5*( el(kts+1)+el(kts) ) /  vkz 
       qke(kts) = ust**2 * ( b1*pmz*elv    )**(2.0/3.0)
!
       phm      = phh*b2 / ( b1*pmz/elv**2 )**(1.0/3.0)
       tsq(kts) = phm*( flt/ust )**2
       qsq(kts) = phm*( flq/ust )**2
       cov(kts) = phm*( flt/ust )*( flq/ust )
!
       DO k = kts+1,kte-1
          b1l = b1*0.25*( el(k+1)+el(k) )
          tmpq=MAX(b1l*( pdk(k+1)+pdk(k) ),qkemin)
!          PRINT *,'tmpqqqqq',tmpq,pdk(k+1),pdk(k)
          qke(k) = tmpq**(2.0/3.0)

!
          IF ( qke(k) .LE. 0.0 ) THEN
             b2l = 0.0
          ELSE
             b2l = b2*( b1l/b1 ) / SQRT( qke(k) )
          END IF
!
          tsq(k) = b2l*( pdt(k+1)+pdt(k) )
          qsq(k) = b2l*( pdq(k+1)+pdq(k) )
          cov(k) = b2l*( pdc(k+1)+pdc(k) )
       END DO

!
    END DO

    qke(kte)=qke(kte-1)
    tsq(kte)=tsq(kte-1)
    qsq(kte)=qsq(kte-1)
    cov(kte)=cov(kte-1)

  END SUBROUTINE mym_initialize
  
!
! ==================================================================
!     SUBROUTINE  mym_level2:
!
!     Input variables:    see subroutine mym_initialize
!
!     Output variables:
!       dtl(nx,nz,ny) : Vertical gradient of Theta_l             (K/m)
!       dqw(nx,nz,ny) : Vertical gradient of Q_w
!       dtv(nx,nz,ny) : Vertical gradient of Theta_V             (K/m)
!       gm (nx,nz,ny) : G_M divided by L^2/q^2                (s^(-2))
!       gh (nx,nz,ny) : G_H divided by L^2/q^2                (s^(-2))
!       sm (nx,nz,ny) : Stability function for momentum, at Level 2
!       sh (nx,nz,ny) : Stability function for heat, at Level 2
!
!       These are defined on the walls of the grid boxes.
!
  SUBROUTINE  mym_level2 (kts,kte,&
       &            dz, &
       &            u, v, thl, qw, &
       &            ql, vt, vq, &
       &            dtl, dqw, dtv, gm, gh, sm, sh, kpbl )
!
!-------------------------------------------------------------------

    INTEGER, INTENT(IN)   :: kts,kte, kpbl



    REAL, DIMENSION(kts:kte), INTENT(in) :: dz
    REAL, DIMENSION(kts:kte), INTENT(in) :: u,v,thl,qw,ql,vt,vq
    REAL, DIMENSION(kts:kte), INTENT(out) :: &
         &dtl,dqw,dtv,gm,gh,sm,sh
    INTEGER :: k

    REAL :: rfc,f1,f2,rf1,rf2,smc,shc,&
         &ri1,ri2,ri3,ri4,duz,dtz,dqz,vtt,vqq,dtq,dzk,afk,abk,ri,rf

!JOE-Canuto/Kitamura mod
    REAL ::   a2den
!JOE-end
!    ev  = 2.5e6
!    tv0 = 0.61*tref
!    tv1 = 1.61*tref
!    gtr = 9.81/tref
!
    rfc = g1/( g1+g2 )
    f1  = b1*( g1-c1 ) +3.0*a2*( 1.0    -c2 )*( 1.0-c5 ) &
    &                   +2.0*a1*( 3.0-2.0*c2 )
    f2  = b1*( g1+g2 ) -3.0*a1*( 1.0    -c2 )
    rf1 = b1*( g1-c1 )/f1
    rf2 = b1*  g1     /f2
    smc = a1 /a2*  f1/f2
    shc = 3.0*a2*( g1+g2 )
!
    ri1 = 0.5/smc
    ri2 = rf1*smc
    ri3 = 4.0*rf2*smc -2.0*ri2
    ri4 = ri2**2

    DO k = kts+1,kte
       dzk = 0.5  *( dz(k)+dz(k-1) )
       afk = dz(k)/( dz(k)+dz(k-1) )
       abk = 1.0 -afk
       duz = ( u(k)-u(k-1) )**2 +( v(k)-v(k-1) )**2
       duz =   duz                    /dzk**2
       dtz = ( thl(k)-thl(k-1) )/( dzk )
       dqz = ( qw(k)-qw(k-1) )/( dzk )

       vtt =  1.0 +vt(k)*abk +vt(k-1)*afk  ! Beta-theta in NN09, Eq. 39
       vqq =  tv0 +vq(k)*abk +vq(k-1)*afk  ! Beta-q
       dtq =  vtt*dtz +vqq*dqz
!
       dtl(k) =  dtz
       dqw(k) =  dqz
       dtv(k) =  dtq
!?      dtv(i,j,k) =  dtz +tv0*dqz
!?   :              +( ev/pi0(i,j,k)-tv1 )
!?   :              *( ql(i,j,k)-ql(i,j,k-1) )/( dzk*h(i,j) )
!
       gm (k) =  duz
       gh (k) = -dtq*gtr
!
!   **  Gradient Richardson number  **
       ri = -gh(k)/MAX( duz, 1.0e-10 )

!JOE-Canuto/Kitamura mod
    IF (CKmod .eq. 1) THEN
       a2den = 1. + MAX(ri,0.0)
    ELSE
       a2den = 1. + 0.0
    ENDIF

       rfc = g1/( g1+g2 )
       f1  = b1*( g1-c1 ) +3.0*(a2/a2den)*( 1.0    -c2 )*( 1.0-c5 ) &
    &                     +2.0*a1*( 3.0-2.0*c2 )
       f2  = b1*( g1+g2 ) -3.0*a1*( 1.0    -c2 )
       rf1 = b1*( g1-c1 )/f1
       rf2 = b1*  g1     /f2
       smc = a1 /(a2/a2den)*  f1/f2
       shc = 3.0*(a2/a2den)*( g1+g2 )

       ri1 = 0.5/smc
       ri2 = rf1*smc
       ri3 = 4.0*rf2*smc -2.0*ri2
       ri4 = ri2**2
!JOE-end

!   **  Flux Richardson number  **
       rf = MIN( ri1*( ri+ri2-SQRT(ri**2-ri3*ri+ri4) ), rfc )
!
       if (bl_mynn_stabfunc .eq. 0)  then   ! yhc_0614, use MYNN's stability function
         sh (k) = shc*( rfc-rf )/( 1.0-rf )
         sm (k) = smc*( rf1-rf )/( rf2-rf ) * sh(k)
! replace stability functions with alphas from Suselj et al. 2019        

       else if (bl_mynn_stabfunc .eq. 1) then  ! yhc_0614, use Suselj et al (2019)'s stability function
         if (ri .le. 0.) then
         ! unstable/neutral layers
          sh(k)=1.4
          sm(k)=1.
          else
          ! stable layers
          sh(k)=(1.4-0.001*ri+1.29*ri**2)/(1.+2.3*ri+19.81*ri**2)
          sm(k)=(1.+8.*ri**2)/(1.+2.3*ri+35.*ri**2)
         endif
       endif

    END DO
!
!    RETURN


  END SUBROUTINE mym_level2

! ==================================================================
!     SUBROUTINE  mym_length:
!
!     Input variables:    see subroutine mym_initialize
!
!     Output variables:   see subroutine mym_initialize
!
!     Work arrays:
!       elt(nx,ny)      : Length scale depending on the PBL depth    (m)
!       vsc(nx,ny)      : Velocity scale q_c                       (m/s)
!                         at first, used for computing elt
!
!     NOTE: the mixing lengths are meant to be calculated at the full-
!           sigmal levels (or interfaces beween the model layers).
!
  SUBROUTINE  mym_length (                     & 
    &            kts,kte,                      &
    &            dz, zw,                       &
    &            rmo, flt, flq,                &
    &            vt, vq,                       &
    &            qke,                          &
    &            dtv,                          &
    &            el,                           &
    &            zi,theta,                     &
    &            qkw,Psig_bl,cldfra_bl1D,bl_mynn_mixlength,&
    &            edmf_w1,edmf_a1,edmf_qc1,bl_mynn_edmf,&
    &            edmf_w_dd1,edmf_a_dd1,edmf_qc_dd1,bl_mynn_edmf_dd,&
    &            rstoch_col)
    
!-------------------------------------------------------------------

    INTEGER, INTENT(IN)   :: kts,kte

    INTEGER, INTENT(IN)   :: bl_mynn_mixlength,bl_mynn_edmf,bl_mynn_edmf_dd
    REAL, DIMENSION(kts:kte), INTENT(in)   :: dz
    REAL, DIMENSION(kts:kte+1), INTENT(in) :: zw
    REAL, INTENT(in) :: rmo,flt,flq,Psig_bl
    REAL, DIMENSION(kts:kte), INTENT(IN)   :: qke,vt,vq,cldfra_bl1D,&
                                          edmf_w1,edmf_a1,edmf_qc1,&
                                          edmf_w_dd1,edmf_a_dd1,edmf_qc_dd1
    REAL, DIMENSION(kts:kte), INTENT(out)  :: qkw, el
    REAL, DIMENSION(kts:kte), INTENT(in)   :: dtv
    REAL, DIMENSION(kts:kte) :: dtv_smoothed
    REAL :: elt,vsc

    REAL, DIMENSION(kts:kte), INTENT(IN) :: theta
    REAL, DIMENSION(kts:kte) :: qtke,elBLmin,elBLavg,thetaw
    REAL :: wt,wt2,zi,zi2,h1,h2,hs,elBLmin0,elBLavg0

    ! THE FOLLOWING CONSTANTS ARE IMPORTANT FOR REGULATING THE
    ! MIXING LENGTHS:
    REAL :: cns,   &   ! for surface layer (els) in stable conditions
            alp1,  &   ! for turbulent length scale (elt)
            alp2,  &   ! for buoyancy length scale (elb)
            alp3,  &   ! for buoyancy enhancement factor of elb
            alp4,  &   ! for surface layer (els) in unstable conditions
            alp5       ! for BouLac mixing length

    !THE FOLLOWING LIMITS DO NOT DIRECTLY AFFECT THE ACTUAL PBLH.
    !THEY ONLY IMPOSE LIMITS ON THE CALCULATION OF THE MIXING LENGTH 
    !SCALES SO THAT THE BOULAC MIXING LENGTH (IN FREE ATMOS) DOES
    !NOT ENCROACH UPON THE BOUNDARY LAYER MIXING LENGTH (els, elb & elt).
    REAL, PARAMETER :: minzi = 300.  !min mixed-layer height
    REAL, PARAMETER :: maxdz = 750.  !max (half) transition layer depth
                                     !=0.3*2500 m PBLH, so the transition
                                     !layer stops growing for PBLHs > 2.5 km.
    REAL, PARAMETER :: mindz = 300.  !300  !min (half) transition layer depth

    !SURFACE LAYER LENGTH SCALE MODS TO REDUCE IMPACT IN UPPER BOUNDARY LAYER
    REAL, PARAMETER :: ZSLH = 100. ! Max height correlated to surface conditions (m)
    REAL, PARAMETER :: CSL = 2.    ! CSL = constant of proportionality to L O(1)
    REAL :: z_m


    INTEGER :: i,j,k
    REAL :: afk,abk,zwk,zwk1,dzk,qdz,vflx,bv,tau_cloud,elb,els,els1,elf, &
            & el_stab,el_unstab,el_mf,el_stab_mf,elb_mf,PBLH_PLUS_ENT,el_les

    REAL, DIMENSION(kts:kte), INTENT(in)   :: rstoch_col

!    tv0 = 0.61*tref
!    gtr = 9.81/tref

    SELECT CASE(bl_mynn_mixlength)

      CASE (0) ! ORIGINAL MYNN MIXING LENGTH

        cns  = 2.7
        alp1 = 0.23
        alp2 = 1.0
        alp3 = 5.0
        alp4 = 100.
        alp5 = 0.4

        ! Impose limits on the height integration for elt and the transition layer depth
        zi2  = MIN(10000.,zw(kte-2))  !originally integrated to model top, not just 10 km.
        h1=MAX(0.3*zi2,mindz)
        h1=MIN(h1,maxdz)         ! 1/2 transition layer depth
        h2=h1/2.0                ! 1/4 transition layer depth

        qkw(kts) = SQRT(MAX(qke(kts),1.0e-10))
        DO k = kts+1,kte
           afk = dz(k)/( dz(k)+dz(k-1) )
           abk = 1.0 -afk
           qkw(k) = SQRT(MAX(qke(k)*abk+qke(k-1)*afk,1.0e-3))
        END DO

        elt = 1.0e-5
        vsc = 1.0e-5        

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        k = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. zi2+h1)
           dzk = 0.5*( dz(k)+dz(k-1) )
           qdz = MAX( qkw(k)-qmin, 0.03 )*dzk
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        elt =  alp1*elt/vsc
        vflx = ( vt(kts)+1.0 )*flt +( vq(kts)+tv0 )*flq
        vsc = ( gtr*elt*MAX( vflx, 0.0 ) )**(1.0/3.0)

        !   **  Strictly, el(i,j,1) is not zero.  **
        el(kts) = 0.0
        zwk1    = zw(kts+1)
        DO k = kts+1,kte
           zwk = zw(k)              !full-sigma levels

           !   **  Length scale limited by the buoyancy effect  **
           IF ( dtv(k) .GT. 0.0 ) THEN
              bv  = SQRT( gtr*dtv(k) )
              elb = alp2*qkw(k) / bv &
                  &       *( 1.0 + alp3/alp2*&
                  &SQRT( vsc/( bv*elt ) ) )
              elf = alp2 * qkw(k)/bv

           ELSE
              elb = 1.0e10
              elf = elb
           ENDIF


           z_m = MAX(0.,zwk - 4.)

           !   **  Length scale in the surface layer  **
           IF ( rmo .GT. 0.0 ) THEN
              els  = vk*zwk/(1.0+cns*MIN( zwk*rmo, zmax ))
              els1 = vk*z_m/(1.0+cns*MIN( zwk*rmo, zmax ))
           ELSE
              els  =  vk*zwk*( 1.0 - alp4* zwk*rmo )**0.2
              els1 =  vk*z_m*( 1.0 - alp4* zwk*rmo )**0.2
           END IF

           !   ** HARMONC AVERGING OF MIXING LENGTH SCALES:
           !       el(k) =      MIN(elb/( elb/elt+elb/els+1.0 ),elf)
           !       el(k) =      elb/( elb/elt+elb/els+1.0 )

           wt=.5*TANH((zwk - (zi2+h1))/h2) + .5

           el(k) = MIN(elb/( elb/elt+elb/els+1.0 ),elf)

        END DO

      CASE (1) !OPERATIONAL FORM OF MIXING LENGTH

        cns  = 2.3
        alp1 = 0.23
        alp2 = 0.65
        alp3 = 3.0
        alp4 = 20.
        alp5 = 0.4

        ! Impose limits on the height integration for elt and the transition layer depth
        zi2=MAX(zi,minzi)
        h1=MAX(0.3*zi2,mindz)
        h1=MIN(h1,maxdz)         ! 1/2 transition layer depth
        h2=h1/2.0                ! 1/4 transition layer depth

        qtke(kts)=MAX(qke(kts)/2.,0.01) !tke at full sigma levels
        thetaw(kts)=theta(kts)          !theta at full-sigma levels
        qkw(kts) = SQRT(MAX(qke(kts),1.0e-10))

        DO k = kts+1,kte
           afk = dz(k)/( dz(k)+dz(k-1) )
           abk = 1.0 -afk
           qkw(k) = SQRT(MAX(qke(k)*abk+qke(k-1)*afk,1.0e-3))
           qtke(k) = (qkw(k)**2.)/2.    ! q -> TKE
           thetaw(k)= theta(k)*abk + theta(k-1)*afk
        END DO

        elt = 1.0e-5
        vsc = 1.0e-5

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        k = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. zi2+h1)
           dzk = 0.5*( dz(k)+dz(k-1) )
           qdz = MAX( qkw(k)-qmin, 0.03 )*dzk
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        elt =  alp1*elt/vsc
        vflx = ( vt(kts)+1.0 )*flt +( vq(kts)+tv0 )*flq
        vsc = ( gtr*elt*MAX( vflx, 0.0 ) )**(1.0/3.0)

        !   **  Strictly, el(i,j,1) is not zero.  **
        el(kts) = 0.0
        zwk1    = zw(kts+1)              !full-sigma levels

        ! COMPUTE BouLac mixing length
        CALL boulac_length(kts,kte,zw,dz,qtke,thetaw,elBLmin,elBLavg)

        DO k = kts+1,kte
           zwk = zw(k)              !full-sigma levels

           !   **  Length scale limited by the buoyancy effect  **
           IF ( dtv(k) .GT. 0.0 ) THEN
              bv  = SQRT( gtr*dtv(k) ) 
              elb = alp2*qkw(k) / bv &                ! formulation,
                  &       *( 1.0 + alp3/alp2*&       ! except keep
                  &SQRT( vsc/( bv*elt ) ) )          ! elb bounded by
              elb = MIN(elb, zwk)                     ! zwk
              elf = alp2 * qkw(k)/bv
           ELSE
              elb = 1.0e10
              elf = elb
           ENDIF

           z_m = MAX(0.,zwk - 4.)

           !   **  Length scale in the surface layer  **
           IF ( rmo .GT. 0.0 ) THEN
              els  = vk*zwk/(1.0+cns*MIN( zwk*rmo, zmax ))
              els1 = vk*z_m/(1.0+cns*MIN( zwk*rmo, zmax ))
           ELSE
              els  =  vk*zwk*( 1.0 - alp4* zwk*rmo )**0.2
              els1 =  vk*z_m*( 1.0 - alp4* zwk*rmo )**0.2
           END IF

           !   ** NOW BLEND THE MIXING LENGTH SCALES:
           wt=.5*TANH((zwk - (zi2+h1))/h2) + .5

           !add blending to use BouLac mixing length in free atmos;
           !defined relative to the PBLH (zi) + transition layer (h1)
           el(k) = MIN(elb/( elb/elt+elb/els+1.0 ),elf)
           el(k) = el(k)*(1.-wt) + alp5*elBLmin(k)*wt

           ! include scale-awareness, except for original MYNN
           el(k) = el(k)*Psig_bl

         END DO

      CASE (2) !Experimental mixing length formulation

        cns  = 3.5
        alp1 = 0.23
        alp2 = 0.3
        alp3 = 2.0
        alp4 = 10.
        alp5 = 0.3 !obsolete?

        ! Impose limits on the height integration for elt and the transition layer depth
        zi2=MAX(zi,minzi)
        h1=MAX(0.3*zi2,mindz)
        h1=MIN(h1,maxdz)         ! 1/2 transition layer depth
        h2=h1/2.0                ! 1/4 transition layer depth

        qtke(kts)=MAX(qke(kts)/2.,0.01) !tke at full sigma levels
        qkw(kts) = SQRT(MAX(qke(kts),1.0e-10))

        DO k = kts+1,kte
           afk = dz(k)/( dz(k)+dz(k-1) )
           abk = 1.0 -afk
           qkw(k) = SQRT(MAX(qke(k)*abk+qke(k-1)*afk,1.0e-3))
           qtke(k) = (qkw(k)**2.)/2.    ! q -> TKE
        END DO

        elt = 1.0e-5
        vsc = 1.0e-5

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        PBLH_PLUS_ENT = MAX(zi, 100.)
        k = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. PBLH_PLUS_ENT)
           dzk = 0.5*( dz(k)+dz(k-1) )
           qdz = MAX( qkw(k)-qmin, 0.03 )*dzk  !consider reducing 0.3
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        elt =  MAX(alp1*elt/vsc, 10.)
        vflx = ( vt(kts)+1.0 )*flt +( vq(kts)+tv0 )*flq
        vsc = ( gtr*elt*MAX( vflx, 0.0 ) )**(1.0/3.0)

        !   **  Strictly, el(i,j,1) is not zero.  **
        el(kts) = 0.0
        zwk1    = zw(kts+1)

        DO k = kts+1,kte
           zwk = zw(k)              !full-sigma levels

           !   **  Length scale limited by the buoyancy effect  **
           IF ( dtv(k) .GT. 0.0 ) THEN
              bv  = SQRT( gtr*dtv(k) )
              !elb_mf = alp2*qkw(k) / bv  &
              elb_mf = alp2*MAX(qkw(k),edmf_a1(k)*edmf_w1(k)) / bv  &
                  &       *( 1.0 + alp3/alp2*&
                  &SQRT( vsc/( bv*elt ) ) )
              elb = MIN(alp2*qkw(k)/bv, zwk)
              elf = elb/(1. + (elb/600.))  !bound free-atmos mixing length to < 600 m.
              !IF (zwk > zi .AND. elf > 400.) THEN
              !   ! COMPUTE BouLac mixing length
              !   !CALL boulac_length0(k,kts,kte,zw,dz,qtke,thetaw,elBLmin0,elBLavg0)
              !   !elf = alp5*elBLavg0
              !   elf = MIN(MAX(50.*SQRT(qtke(k)), 400.), zwk)
              !ENDIF

           ELSE
              ! use version in development for RAP/HRRR 2016
              ! JAYMES-
              ! tau_cloud is an eddy turnover timescale;
              ! see Teixeira and Cheinet (2004), Eq. 1, and
              ! Cheinet and Teixeira (2003), Eq. 7.  The
              ! coefficient 0.5 is tuneable. Expression in
              ! denominator is identical to vsc (a convective
              ! velocity scale), except that elt is relpaced
              ! by zi, and zero is replaced by 1.0e-4 to
              ! prevent division by zero.
              tau_cloud = MIN(MAX(0.5*zi/((gtr*zi*MAX(flt,1.0e-4))**(1.0/3.0)),10.),100.)
              !minimize influence of surface heat flux on tau far away from the PBLH.
              wt=.5*TANH((zwk - (zi2+h1))/h2) + .5
              tau_cloud = tau_cloud*(1.-wt) + 50.*wt

              elb = MIN(tau_cloud*SQRT(MIN(qtke(k),50.)), zwk)
              elf = elb
              elb_mf = elb
         END IF

         z_m = MAX(0.,zwk - 4.)

         !   **  Length scale in the surface layer  **
         IF ( rmo .GT. 0.0 ) THEN
            els  = vk*zwk/(1.0+cns*MIN( zwk*rmo, zmax ))
            els1 = vk*z_m/(1.0+cns*MIN( zwk*rmo, zmax ))
         ELSE
            els  =  vk*zwk*( 1.0 - alp4* zwk*rmo )**0.2
            els1 =  vk*z_m*( 1.0 - alp4* zwk*rmo )**0.2
         END IF

         !   ** NOW BLEND THE MIXING LENGTH SCALES:
         wt=.5*TANH((zwk - (zi2+h1))/h2) + .5

         ! "el_unstab" = blended els-elt
         el_unstab = els/(1. + (els1/elt))
         el(k) = MIN(el_unstab, elb_mf)
         el(k) = el(k)*(1.-wt) + elf*wt

         ! include scale-awareness. For now, use simple asymptotic kz -> 12 m.
         el_les= MIN(els/(1. + (els1/12.)), elb_mf)
         el(k) = el(k)*Psig_bl + (1.-Psig_bl)*el_les

       END DO

    END SELECT


!   Stochastic perturbations of turbulent mixing length
!    if (spp_pbl==1) then
!       DO k = kts+1,kte
!         if (k.lt.25) then
!            zwk = zw(k)
!            el(k)= el(k) + el(k)* rstoch_col(k) * 1.5 * MAX(exp(-MAX(zwk-3000.,0.0)/2000.),0.01)
!         endif
!       END DO
!    endif


  END SUBROUTINE mym_length

!JOE- BouLac Code Start -

! ==================================================================
  SUBROUTINE boulac_length0(k,kts,kte,zw,dz,qtke,theta,lb1,lb2)
!
!    NOTE: This subroutine was taken from the BouLac scheme in WRF-ARW
!          and modified for integration into the MYNN PBL scheme.
!          WHILE loops were added to reduce the computational expense.
!          This subroutine computes the length scales up and down
!          and then computes the min, average of the up/down
!          length scales, and also considers the distance to the
!          surface.
!
!      dlu = the distance a parcel can be lifted upwards give a finite
!            amount of TKE.
!      dld = the distance a parcel can be displaced downwards given a
!            finite amount of TKE.
!      lb1 = the minimum of the length up and length down
!      lb2 = the average of the length up and length down
!-------------------------------------------------------------------

     INTEGER, INTENT(IN) :: k,kts,kte
     REAL, DIMENSION(kts:kte), INTENT(IN) :: qtke,dz,theta
     REAL, INTENT(OUT) :: lb1,lb2
     REAL, DIMENSION(kts:kte+1), INTENT(IN) :: zw

     !LOCAL VARS
     INTEGER :: izz, found
     REAL :: dlu,dld
     REAL :: dzt, zup, beta, zup_inf, bbb, tl, zdo, zdo_sup, zzz


     !----------------------------------
     ! FIND DISTANCE UPWARD             
     !----------------------------------
     zup=0.
     dlu=zw(kte+1)-zw(k)-dz(k)/2.
     zzz=0.
     zup_inf=0.
     beta=g/theta(k)           !Buoyancy coefficient

     !print*,"FINDING Dup, k=",k," zw=",zw(k)

     if (k .lt. kte) then      !cant integrate upwards from highest level
        found = 0
        izz=k
        DO WHILE (found .EQ. 0)

           if (izz .lt. kte) then
              dzt=dz(izz)                    ! layer depth above
              zup=zup-beta*theta(k)*dzt     ! initial PE the parcel has at k
              !print*,"  ",k,izz,theta(izz),dz(izz)
              zup=zup+beta*(theta(izz+1)+theta(izz))*dzt/2. ! PE gained by lifting a parcel to izz+1
              zzz=zzz+dzt                   ! depth of layer k to izz+1
              !print*,"  PE=",zup," TKE=",qtke(k)," z=",zw(izz)
              if (qtke(k).lt.zup .and. qtke(k).ge.zup_inf) then
                 bbb=(theta(izz+1)-theta(izz))/dzt
                 if (bbb .ne. 0.) then
                    !fractional distance up into the layer where TKE becomes < PE
                    tl=(-beta*(theta(izz)-theta(k)) + &
                      & sqrt( max(0.,(beta*(theta(izz)-theta(k)))**2. + &
                      &       2.*bbb*beta*(qtke(k)-zup_inf))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(k))then
                       tl=(qtke(k)-zup_inf)/(beta*(theta(izz)-theta(k)))
                    else
                       tl=0.
                    endif
                 endif
                 dlu=zzz-dzt+tl
                 !print*,"  FOUND Dup:",dlu," z=",zw(izz)," tl=",tl
                 found =1
              endif
              zup_inf=zup
              izz=izz+1
           ELSE
              found = 1
           ENDIF

        ENDDO

     endif

     !----------------------------------
     ! FIND DISTANCE DOWN               
     !----------------------------------
     zdo=0.
     zdo_sup=0.
     dld=zw(k)
     zzz=0.

     !print*,"FINDING Ddown, k=",k," zwk=",zw(k)
     if (k .gt. kts) then  !cant integrate downwards from lowest level

        found = 0
        izz=k
        DO WHILE (found .EQ. 0)

           if (izz .gt. kts) then
              dzt=dz(izz-1)
              zdo=zdo+beta*theta(k)*dzt
              !print*,"  ",k,izz,theta(izz),dz(izz-1)
              zdo=zdo-beta*(theta(izz-1)+theta(izz))*dzt/2.
              zzz=zzz+dzt
              !print*,"  PE=",zdo," TKE=",qtke(k)," z=",zw(izz)
              if (qtke(k).lt.zdo .and. qtke(k).ge.zdo_sup) then
                 bbb=(theta(izz)-theta(izz-1))/dzt
                 if (bbb .ne. 0.) then
                    tl=(beta*(theta(izz)-theta(k))+ &
                      & sqrt( max(0.,(beta*(theta(izz)-theta(k)))**2. + &
                      &       2.*bbb*beta*(qtke(k)-zdo_sup))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(k)) then
                       tl=(qtke(k)-zdo_sup)/(beta*(theta(izz)-theta(k)))
                    else
                       tl=0.
                    endif
                 endif
                 dld=zzz-dzt+tl
                 !print*,"  FOUND Ddown:",dld," z=",zw(izz)," tl=",tl
                 found = 1
              endif
              zdo_sup=zdo
              izz=izz-1
           ELSE
              found = 1
           ENDIF
        ENDDO

     endif

     !----------------------------------
     ! GET MINIMUM (OR AVERAGE)         
     !----------------------------------
     !The surface layer length scale can exceed z for large z/L,
     !so keep maximum distance down > z.
     dld = min(dld,zw(k+1))!not used in PBL anyway, only free atmos
     lb1 = min(dlu,dld)     !minimum
     !JOE-fight floating point errors
     dlu=MAX(0.1,MIN(dlu,1000.))
     dld=MAX(0.1,MIN(dld,1000.))
     lb2 = sqrt(dlu*dld)    !average - biased towards smallest
     !lb2 = 0.5*(dlu+dld)   !average

     if (k .eq. kte) then
        lb1 = 0.
        lb2 = 0.
     endif
     !print*,"IN MYNN-BouLac",k,lb1
     !print*,"IN MYNN-BouLac",k,dld,dlu

  END SUBROUTINE boulac_length0

! ==================================================================
  SUBROUTINE boulac_length(kts,kte,zw,dz,qtke,theta,lb1,lb2)
!
!    NOTE: This subroutine was taken from the BouLac scheme in WRF-ARW
!          and modified for integration into the MYNN PBL scheme.
!          WHILE loops were added to reduce the computational expense.
!          This subroutine computes the length scales up and down
!          and then computes the min, average of the up/down
!          length scales, and also considers the distance to the
!          surface.
!
!      dlu = the distance a parcel can be lifted upwards give a finite 
!            amount of TKE.
!      dld = the distance a parcel can be displaced downwards given a
!            finite amount of TKE.
!      lb1 = the minimum of the length up and length down
!      lb2 = the average of the length up and length down
!-------------------------------------------------------------------

     INTEGER, INTENT(IN) :: kts,kte
     REAL, DIMENSION(kts:kte), INTENT(IN) :: qtke,dz,theta
     REAL, DIMENSION(kts:kte), INTENT(OUT) :: lb1,lb2
     REAL, DIMENSION(kts:kte+1), INTENT(IN) :: zw

     !LOCAL VARS
     INTEGER :: iz, izz, found
     REAL, DIMENSION(kts:kte) :: dlu,dld
     REAL, PARAMETER :: Lmax=2000.  !soft limit
     REAL :: dzt, zup, beta, zup_inf, bbb, tl, zdo, zdo_sup, zzz

     !print*,"IN MYNN-BouLac",kts, kte

     do iz=kts,kte

        !----------------------------------
        ! FIND DISTANCE UPWARD
        !----------------------------------
        zup=0.
        dlu(iz)=zw(kte+1)-zw(iz)-dz(iz)/2.
        zzz=0.
        zup_inf=0.
        beta=g/theta(iz)           !Buoyancy coefficient

        !print*,"FINDING Dup, k=",iz," zw=",zw(iz)

        if (iz .lt. kte) then      !cant integrate upwards from highest level

          found = 0
          izz=iz       
          DO WHILE (found .EQ. 0) 

            if (izz .lt. kte) then
              dzt=dz(izz)                    ! layer depth above 
              zup=zup-beta*theta(iz)*dzt     ! initial PE the parcel has at iz
              !print*,"  ",iz,izz,theta(izz),dz(izz)
              zup=zup+beta*(theta(izz+1)+theta(izz))*dzt/2. ! PE gained by lifting a parcel to izz+1
              zzz=zzz+dzt                   ! depth of layer iz to izz+1
              !print*,"  PE=",zup," TKE=",qtke(iz)," z=",zw(izz)
              if (qtke(iz).lt.zup .and. qtke(iz).ge.zup_inf) then
                 bbb=(theta(izz+1)-theta(izz))/dzt
                 if (bbb .ne. 0.) then
                    !fractional distance up into the layer where TKE becomes < PE
                    tl=(-beta*(theta(izz)-theta(iz)) + &
                      & sqrt( max(0.,(beta*(theta(izz)-theta(iz)))**2. + &
                      &       2.*bbb*beta*(qtke(iz)-zup_inf))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(iz))then
                       tl=(qtke(iz)-zup_inf)/(beta*(theta(izz)-theta(iz)))
                    else
                       tl=0.
                    endif
                 endif            
                 dlu(iz)=zzz-dzt+tl
                 !print*,"  FOUND Dup:",dlu(iz)," z=",zw(izz)," tl=",tl
                 found =1
              endif
              zup_inf=zup
              izz=izz+1
             ELSE
              found = 1
            ENDIF

          ENDDO

        endif
                   
        !----------------------------------
        ! FIND DISTANCE DOWN
        !----------------------------------
        zdo=0.
        zdo_sup=0.
        dld(iz)=zw(iz)
        zzz=0.

        !print*,"FINDING Ddown, k=",iz," zwk=",zw(iz)
        if (iz .gt. kts) then  !cant integrate downwards from lowest level

          found = 0
          izz=iz       
          DO WHILE (found .EQ. 0) 

            if (izz .gt. kts) then
              dzt=dz(izz-1)
              zdo=zdo+beta*theta(iz)*dzt
              !print*,"  ",iz,izz,theta(izz),dz(izz-1)
              zdo=zdo-beta*(theta(izz-1)+theta(izz))*dzt/2.
              zzz=zzz+dzt
              !print*,"  PE=",zdo," TKE=",qtke(iz)," z=",zw(izz)
              if (qtke(iz).lt.zdo .and. qtke(iz).ge.zdo_sup) then
                 bbb=(theta(izz)-theta(izz-1))/dzt
                 if (bbb .ne. 0.) then
                    tl=(beta*(theta(izz)-theta(iz))+ &
                      & sqrt( max(0.,(beta*(theta(izz)-theta(iz)))**2. + &
                      &       2.*bbb*beta*(qtke(iz)-zdo_sup))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(iz)) then
                       tl=(qtke(iz)-zdo_sup)/(beta*(theta(izz)-theta(iz)))
                    else
                       tl=0.
                    endif
                 endif            
                 dld(iz)=zzz-dzt+tl
                 !print*,"  FOUND Ddown:",dld(iz)," z=",zw(izz)," tl=",tl
                 found = 1
              endif
              zdo_sup=zdo
              izz=izz-1
            ELSE
              found = 1
            ENDIF
          ENDDO

        endif

        !----------------------------------
        ! GET MINIMUM (OR AVERAGE)
        !----------------------------------
        !The surface layer length scale can exceed z for large z/L,
        !so keep maximum distance down > z.
        dld(iz) = min(dld(iz),zw(iz+1))!not used in PBL anyway, only free atmos
        lb1(iz) = min(dlu(iz),dld(iz))     !minimum
        !JOE-fight floating point errors
        dlu(iz)=MAX(0.1,MIN(dlu(iz),1000.))
        dld(iz)=MAX(0.1,MIN(dld(iz),1000.))
        lb2(iz) = sqrt(dlu(iz)*dld(iz))    !average - biased towards smallest
        !lb2(iz) = 0.5*(dlu(iz)+dld(iz))   !average

        !Apply soft limit (only impacts very large lb; lb=100 by 5%, lb=500 by 20%).
        lb1(iz) = lb1(iz)/(1. + (lb1(iz)/Lmax))
        lb2(iz) = lb2(iz)/(1. + (lb2(iz)/Lmax))
 
        if (iz .eq. kte) then
           lb1(kte) = lb1(kte-1)
           lb2(kte) = lb2(kte-1)
        endif
        !print*,"IN MYNN-BouLac",kts, kte,lb1(iz)
        !print*,"IN MYNN-BouLac",iz,dld(iz),dlu(iz)

     ENDDO
                   
  END SUBROUTINE boulac_length
!
!JOE-END BOULAC CODE

! ==================================================================
!     SUBROUTINE  mym_turbulence:
!
!     Input variables:    see subroutine mym_initialize
!       levflag         : <>3;  Level 2.5
!                         = 3;  Level 3
!
!     # ql, vt, vq, qke, tsq, qsq and cov are changed to input variables.
!
!     Output variables:   see subroutine mym_initialize
!       dfm(nx,nz,ny) : Diffusivity coefficient for momentum,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       dfh(nx,nz,ny) : Diffusivity coefficient for heat,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       dfq(nx,nz,ny) : Diffusivity coefficient for q^2,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       tcd(nx,nz,ny)   : Countergradient diffusion term for Theta_l
!                                                                  (K/s)
!       qcd(nx,nz,ny)   : Countergradient diffusion term for Q_w
!                                                              (kg/kg s)
!       pd?(nx,nz,ny) : Half of the production terms
!
!       Only tcd and qcd are defined at the center of the grid boxes
!
!     # DO NOT forget that tcd and qcd are added on the right-hand side
!       of the equations for Theta_l and Q_w, respectively.
!
!     Work arrays:        see subroutine mym_initialize and level2
!
!     # dtl, dqw, dtv, gm and gh are allowed to share storage units with
!       dfm, dfh, dfq, tcd and qcd, respectively, for saving memory.
!
  SUBROUTINE  mym_turbulence (                                &
    &            kts,kte,                                     &
    &            levflag,                                     &
    &            dz, zw,                                      &
    &            u, v, thl, ql, qw,                           &
    &            qke, tsq, qsq, cov,                          &
    &            vt, vq,                                      &
    &            rmo, flt, flq,                               &
    &            zi,theta,                                    &
    &            sh,                                          &
    &            El,                                          &
    &            Dfm, Dfh, Dfq, Tcd, Qcd, Pdk, Pdt, Pdq, Pdc, &
    &        qWT1D,qSHEAR1D,qBUOY1D,qDISS1D,              &
    &            bl_mynn_tkebudget,                           &
    &            Psig_bl,Psig_shcu,cldfra_bl1D,bl_mynn_mixlength,&
    &            edmf_w1,edmf_a1,edmf_qc1,bl_mynn_edmf,       &
    &            edmf_w_dd1,edmf_a_dd1,edmf_qc_dd1,bl_mynn_edmf_dd,&
    &            TKEprodTD,                                   &
    &            rstoch_col,kpbl,KHtopdown)

!-------------------------------------------------------------------
!
    INTEGER, INTENT(IN)   :: kts,kte

    INTEGER, INTENT(IN)   :: levflag,bl_mynn_mixlength,bl_mynn_edmf,bl_mynn_edmf_dd,kpbl
    REAL, DIMENSION(kts:kte), INTENT(in) :: dz
    REAL, DIMENSION(kts:kte+1), INTENT(in) :: zw
    REAL, INTENT(in) :: rmo,flt,flq,Psig_bl,Psig_shcu
    REAL, DIMENSION(kts:kte), INTENT(in) :: u,v,thl,qw,& 
         &ql,vt,vq,qke,tsq,qsq,cov,cldfra_bl1D,edmf_w1,edmf_a1,edmf_qc1,&
         &TKEprodTD,edmf_w_dd1,edmf_a_dd1,edmf_qc_dd1, KHtopdown

    REAL, DIMENSION(kts:kte), INTENT(out) :: dfm,dfh,dfq,&
         &pdk,pdt,pdq,pdc,tcd,qcd,el

    REAL, DIMENSION(kts:kte), INTENT(inout) :: &
         qWT1D,qSHEAR1D,qBUOY1D,qDISS1D
    REAL :: q3sq_old,dlsq1,qWTP_old,qWTP_new
    REAL :: dudz,dvdz,dTdz,&
            upwp,vpwp,Tpwp
    INTEGER, INTENT(in) :: bl_mynn_tkebudget

    REAL, DIMENSION(kts:kte) :: qkw,dtl,dqw,dtv,gm,gh,sm,sh,sh_temp,gh_temp

    INTEGER :: k
!    REAL :: cc2,cc3,e1c,e2c,e3c,e4c,e5c
    REAL :: e6c,dzk,afk,abk,vtt,vqq,&
         &cw25,clow,cupp,gamt,gamq,smd,gamv,elq,elh

    REAL :: zi
    REAL, DIMENSION(kts:kte), INTENT(in) :: theta

    REAL ::  a2den, duz, ri, HLmod  !JOE-Canuto/Kitamura mod
!JOE-stability criteria for cw
    REAL:: auh,aum,adh,adm,aeh,aem,Req,Rsl,Rsl2
!JOE-end

    DOUBLE PRECISION  q2sq, t2sq, r2sq, c2sq, elsq, gmel, ghel
    DOUBLE PRECISION  q3sq, t3sq, r3sq, c3sq, dlsq, qdiv
    DOUBLE PRECISION  e1, e2, e3, e4, enum, eden, wden

!   Stochastic

    REAL, DIMENSION(KTS:KTE)                      ::    rstoch_col
    REAL :: prlimit


!
!    tv0 = 0.61*tref
!    gtr = 9.81/tref
!
!    cc2 =  1.0-c2
!    cc3 =  1.0-c3
!    e1c =  3.0*a2*b2*cc3
!    e2c =  9.0*a1*a2*cc2
!    e3c =  9.0*a2*a2*cc2*( 1.0-c5 )
!    e4c = 12.0*a1*a2*cc2
!    e5c =  6.0*a1*a1
!

    CALL mym_level2 (kts,kte,&
    &            dz, &
    &            u, v, thl, qw, &
    &            ql, vt, vq, &
    &            dtl, dqw, dtv, gm, gh, sm, sh, kpbl)
!
    CALL mym_length (                           &
    &            kts,kte,                       &
    &            dz, zw,                        &
    &            rmo, flt, flq,                 &
    &            vt, vq,                        &
    &            qke,                           &
    &            dtv,                           &
    &            el,                            &
    &            zi,theta,                      &
    &            qkw,Psig_bl,cldfra_bl1D,bl_mynn_mixlength, &
    &            edmf_w1,edmf_a1,edmf_qc1,bl_mynn_edmf, &
    &            edmf_w_dd1,edmf_a_dd1,edmf_qc_dd1,bl_mynn_edmf_dd, &
    &            rstoch_col)
    
    DO k = kts+1,kte
       dzk = 0.5  *( dz(k)+dz(k-1) )
       afk = dz(k)/( dz(k)+dz(k-1) )
       abk = 1.0 -afk
       elsq = el (k)**2
       q2sq = b1*elsq*( sm(k)*gm(k)+sh(k)*gh(k) )
       q3sq = qkw(k)**2

!JOE-Canuto/Kitamura mod
       duz = ( u(k)-u(k-1) )**2 +( v(k)-v(k-1) )**2
       duz =   duz                    /dzk**2
       !   **  Gradient Richardson number  **
       ri = -gh(k)/MAX( duz, 1.0e-10 )
       IF (CKmod .eq. 1) THEN
          a2den = 1. + MAX(ri,0.0)
       ELSE
          a2den = 1. + 0.0
       ENDIF
!JOE-end
!
!  Modified: Dec/22/2005, from here, (dlsq -> elsq)
       gmel = gm (k)*elsq
       ghel = gh (k)*elsq
!  Modified: Dec/22/2005, up to here

       ! Level 2.0 debug prints
       IF ( debug_code ) THEN
         IF (sh(k)<0.0 .OR. sm(k)<0.0) THEN
           print*,"MYNN; mym_turbulence2.0; sh=",sh(k)," k=",k
           print*," gm=",gm(k)," gh=",gh(k)," sm=",sm(k)
           print*," q2sq=",q2sq," q3sq=",q3sq," q3/q2=",q3sq/q2sq
           print*," qke=",qke(k)," el=",el(k)," ri=",ri
           print*," PBLH=",zi," u=",u(k)," v=",v(k)
         ENDIF
       ENDIF

!JOE-Apply Helfand & Labraga stability check for all Ric
!      when CKmod == 1. (currently not forced below)
       IF (CKmod .eq. 1) THEN
          HLmod = q2sq -1.
       ELSE
          HLmod = q3sq
       ENDIF

!     **  Since qkw is set to more than 0.0, q3sq > 0.0.  **

!JOE-test new stability criteria in level 2.5 (as well as level 3) - little/no impact
!     **  Limitation on q, instead of L/q  **
          dlsq =  elsq
          IF ( q3sq/dlsq .LT. -gh(k) ) q3sq = -dlsq*gh(k)
!JOE-end

       IF ( q3sq .LT. q2sq ) THEN
       !IF ( HLmod .LT. q2sq ) THEN
          !Apply Helfand & Labraga mod
          qdiv = SQRT( q3sq/q2sq )   !HL89: (1-alfa)
          sm(k) = sm(k) * qdiv
          sh(k) = sh(k) * qdiv
!
          !JOE-Canuto/Kitamura mod
          !e1   = q3sq - e1c*ghel * qdiv**2
          !e2   = q3sq - e2c*ghel * qdiv**2
          !e3   = e1   + e3c*ghel * qdiv**2
          !e4   = e1   - e4c*ghel * qdiv**2
          e1   = q3sq - e1c*ghel/a2den * qdiv**2
          e2   = q3sq - e2c*ghel/a2den * qdiv**2
          e3   = e1   + e3c*ghel/(a2den**2) * qdiv**2
          e4   = e1   - e4c*ghel/a2den * qdiv**2
          eden = e2*e4 + e3*e5c*gmel * qdiv**2
          eden = MAX( eden, 1.0d-20 )
       ELSE
          !JOE-Canuto/Kitamura mod
          !e1   = q3sq - e1c*ghel
          !e2   = q3sq - e2c*ghel
          !e3   = e1   + e3c*ghel
          !e4   = e1   - e4c*ghel
          e1   = q3sq - e1c*ghel/a2den
          e2   = q3sq - e2c*ghel/a2den
          e3   = e1   + e3c*ghel/(a2den**2)
          e4   = e1   - e4c*ghel/a2den
          eden = e2*e4 + e3*e5c*gmel
          eden = MAX( eden, 1.0d-20 )

          qdiv = 1.0
          sm(k) = q3sq*a1*( e3-3.0*c1*e4       )/eden
          !JOE-Canuto/Kitamura mod
          !sh(k) = q3sq*a2*( e2+3.0*c1*e5c*gmel )/eden
          sh(k) = q3sq*(a2/a2den)*( e2+3.0*c1*e5c*gmel )/eden
       END IF !end Helfand & Labraga check

       !JOE: Level 2.5 debug prints
       ! HL88 , lev2.5 criteria from eqs. 3.17, 3.19, & 3.20
       IF ( debug_code ) THEN
         IF (sh(k)<0.0 .OR. sm(k)<0.0 .OR. &
           sh(k) > 0.76*b2 .or. (sm(k)**2*gm(k) .gt. .44**2)) THEN
           print*,"MYNN; mym_turbulence2.5; sh=",sh(k)," k=",k
           print*," gm=",gm(k)," gh=",gh(k)," sm=",sm(k)
           print*," q2sq=",q2sq," q3sq=",q3sq," q3/q2=",q3sq/q2sq
           print*," qke=",qke(k)," el=",el(k)," ri=",ri
           print*," PBLH=",zi," u=",u(k)," v=",v(k)
         ENDIF
       ENDIF

!   **  Level 3 : start  **
       IF ( levflag .EQ. 3 ) THEN
          t2sq = qdiv*b2*elsq*sh(k)*dtl(k)**2
          r2sq = qdiv*b2*elsq*sh(k)*dqw(k)**2
          c2sq = qdiv*b2*elsq*sh(k)*dtl(k)*dqw(k)
          t3sq = MAX( tsq(k)*abk+tsq(k-1)*afk, 0.0 )
          r3sq = MAX( qsq(k)*abk+qsq(k-1)*afk, 0.0 )
          c3sq =      cov(k)*abk+cov(k-1)*afk

!  Modified: Dec/22/2005, from here
          c3sq = SIGN( MIN( ABS(c3sq), SQRT(t3sq*r3sq) ), c3sq )
!
          vtt  = 1.0 +vt(k)*abk +vt(k-1)*afk
          vqq  = tv0 +vq(k)*abk +vq(k-1)*afk
          t2sq = vtt*t2sq +vqq*c2sq
          r2sq = vtt*c2sq +vqq*r2sq
          c2sq = MAX( vtt*t2sq+vqq*r2sq, 0.0d0 )
          t3sq = vtt*t3sq +vqq*c3sq
          r3sq = vtt*c3sq +vqq*r3sq
          c3sq = MAX( vtt*t3sq+vqq*r3sq, 0.0d0 )
!
          cw25 = e1*( e2 + 3.0*c1*e5c*gmel*qdiv**2 )/( 3.0*eden )
!
!     **  Limitation on q, instead of L/q  **
          dlsq =  elsq
          IF ( q3sq/dlsq .LT. -gh(k) ) q3sq = -dlsq*gh(k)
!
!     **  Limitation on c3sq (0.12 =< cw =< 0.76) **
          !JOE: use Janjic's (2001; p 13-17) methodology (eqs 4.11-414 and 5.7-5.10)
          ! to calculate an exact limit for c3sq:
          auh = 27.*a1*((a2/a2den)**2)*b2*(g/tref)**2
          aum = 54.*(a1**2)*(a2/a2den)*b2*c1*(g/tref)
          adh = 9.*a1*((a2/a2den)**2)*(12.*a1 + 3.*b2)*(g/tref)**2
          adm = 18.*(a1**2)*(a2/a2den)*(b2 - 3.*(a2/a2den))*(g/tref)

          aeh = (9.*a1*((a2/a2den)**2)*b1 +9.*a1*((a2/a2den)**2)* &
                (12.*a1 + 3.*b2))*(g/tref)
          aem = 3.*a1*(a2/a2den)*b1*(3.*(a2/a2den) + 3.*b2*c1 + &
                (18.*a1*c1 - b2)) + &
                (18.)*(a1**2)*(a2/a2den)*(b2 - 3.*(a2/a2den))

          Req = -aeh/aem
          Rsl = (auh + aum*Req)/(3.*adh + 3.*adm*Req)
          !For now, use default values, since tests showed little/no sensitivity
          Rsl = .12             !lower limit
          Rsl2= 1.0 - 2.*Rsl    !upper limit
          !IF (k==2)print*,"Dynamic limit RSL=",Rsl
          !IF (Rsl < 0.10 .OR. Rsl > 0.18) THEN
          !   wrf_err_message = '--- ERROR: MYNN: Dynamic Cw '// &
          !        'limit exceeds reasonable limits'
          !   CALL wrf_message ( wrf_err_message )
          !   WRITE ( mynn_message , FMT='(A,F8.3)' ) &
          !   " MYNN: Dynamic Cw limit needs attention=",Rsl
          !   CALL wrf_debug ( 0 , mynn_message )
          !ENDIF

          !JOE-Canuto/Kitamura mod
          !e2   = q3sq - e2c*ghel * qdiv**2
          !e3   = q3sq + e3c*ghel * qdiv**2
          !e4   = q3sq - e4c*ghel * qdiv**2
          e2   = q3sq - e2c*ghel/a2den * qdiv**2
          e3   = q3sq + e3c*ghel/(a2den**2) * qdiv**2
          e4   = q3sq - e4c*ghel/a2den * qdiv**2
          eden = e2*e4  + e3 *e5c*gmel * qdiv**2

          !JOE-Canuto/Kitamura mod
          !wden = cc3*gtr**2 * dlsq**2/elsq * qdiv**2 &
          !     &        *( e2*e4c - e3c*e5c*gmel * qdiv**2 )
          wden = cc3*gtr**2 * dlsq**2/elsq * qdiv**2 &
               &        *( e2*e4c/a2den - e3c*e5c*gmel/(a2den**2) * qdiv**2 )

          IF ( wden .NE. 0.0 ) THEN
             !JOE: test dynamic limits
             !clow = q3sq*( 0.12-cw25 )*eden/wden
             !cupp = q3sq*( 0.76-cw25 )*eden/wden
             clow = q3sq*( Rsl -cw25 )*eden/wden
             cupp = q3sq*( Rsl2-cw25 )*eden/wden
!
             IF ( wden .GT. 0.0 ) THEN
                c3sq  = MIN( MAX( c3sq, c2sq+clow ), c2sq+cupp )
             ELSE
                c3sq  = MAX( MIN( c3sq, c2sq+clow ), c2sq+cupp )
             END IF
          END IF
!
          e1   = e2 + e5c*gmel * qdiv**2
          eden = MAX( eden, 1.0d-20 )
!  Modified: Dec/22/2005, up to here

          !JOE-Canuto/Kitamura mod
          !e6c  = 3.0*a2*cc3*gtr * dlsq/elsq
          e6c  = 3.0*(a2/a2den)*cc3*gtr * dlsq/elsq

          !============================
          !     **  for Gamma_theta  **
          !!          enum = qdiv*e6c*( t3sq-t2sq )
          IF ( t2sq .GE. 0.0 ) THEN
             enum = MAX( qdiv*e6c*( t3sq-t2sq ), 0.0d0 )
          ELSE
             enum = MIN( qdiv*e6c*( t3sq-t2sq ), 0.0d0 )
          ENDIF
          gamt =-e1  *enum    /eden

          !============================
          !     **  for Gamma_q  **
          !!          enum = qdiv*e6c*( r3sq-r2sq )
          IF ( r2sq .GE. 0.0 ) THEN
             enum = MAX( qdiv*e6c*( r3sq-r2sq ), 0.0d0 )
          ELSE
             enum = MIN( qdiv*e6c*( r3sq-r2sq ), 0.0d0 )
          ENDIF
          gamq =-e1  *enum    /eden

          !============================
          !     **  for Sm' and Sh'd(Theta_V)/dz  **
          !!          enum = qdiv*e6c*( c3sq-c2sq )
          enum = MAX( qdiv*e6c*( c3sq-c2sq ), 0.0d0)

          !JOE-Canuto/Kitamura mod
          !smd  = dlsq*enum*gtr/eden * qdiv**2 * (e3c+e4c)*a1/a2
          smd  = dlsq*enum*gtr/eden * qdiv**2 * (e3c/(a2den**2) + &
               & e4c/a2den)*a1/(a2/a2den)

          gamv = e1  *enum*gtr/eden
          sm(k) = sm(k) +smd

          !============================
          !     **  For elh (see below), qdiv at Level 3 is reset to 1.0.  **
          qdiv = 1.0

          ! Level 3 debug prints
          IF ( debug_code ) THEN
            IF (sh(k)<-0.3 .OR. sm(k)<-0.3 .OR. &
              qke(k) < -0.1 .or. ABS(smd) .gt. 2.0) THEN
              print*," MYNN; mym_turbulence3.0; sh=",sh(k)," k=",k
              print*," gm=",gm(k)," gh=",gh(k)," sm=",sm(k)
              print*," q2sq=",q2sq," q3sq=",q3sq," q3/q2=",q3sq/q2sq
              print*," qke=",qke(k)," el=",el(k)," ri=",ri
              print*," PBLH=",zi," u=",u(k)," v=",v(k)
            ENDIF
          ENDIF

!   **  Level 3 : end  **

       ELSE
!     **  At Level 2.5, qdiv is not reset.  **
          gamt = 0.0
          gamq = 0.0
          gamv = 0.0
       END IF
!
!      Add stochastic perturbation of prandtl number limit
!       if (spp_pbl==1) then
!          prlimit = MIN(MAX(1.,2.5 + 5.0*rstoch_col(k)), 10.)
!          IF(sm(k) > sh(k)*Prlimit) THEN
!             sm(k) = sh(k)*Prlimit
!          ENDIF
!       ENDIF
!
       elq = el(k)*qkw(k)
       elh = elq*qdiv

       ! Production of TKE (pdk), T-variance (pdt),
       ! q-variance (pdq), and covariance (pdc)
       pdk(k) = elq*( sm(k)*gm(k) &
            &                    +sh(k)*gh(k)+gamv ) + & ! JAYMES TKE
            &   TKEprodTD(k)                             ! JOE-top-down
       pdt(k) = elh*( sh(k)*dtl(k)+gamt )*dtl(k)
       pdq(k) = elh*( sh(k)*dqw(k)+gamq )*dqw(k)
       pdc(k) = elh*( sh(k)*dtl(k)+gamt )&
            &*dqw(k)*0.5 &
                  &+elh*( sh(k)*dqw(k)+gamq )*dtl(k)*0.5

       ! Contergradient terms
       tcd(k) = elq*gamt
       qcd(k) = elq*gamq

       ! Eddy Diffusivity/Viscosity divided by dz
       dfm(k) = elq*sm(k) / dzk
       dfh(k) = elq*sh(k) / dzk + KHtopdown(k) / dzk !EW: tagged on topdown ED
!  Modified: Dec/22/2005, from here
!   **  In sub.mym_predict, dfq for the TKE and scalar variance **
!   **  are set to 3.0*dfm and 1.0*dfm, respectively. (Sqfac)   **
       dfq(k) =     dfm(k)
!  Modified: Dec/22/2005, up to here

   IF ( bl_mynn_tkebudget == 1) THEN
       !TKE BUDGET
       dudz = ( u(k)-u(k-1) )/dzk
       dvdz = ( v(k)-v(k-1) )/dzk
       dTdz = ( thl(k)-thl(k-1) )/dzk

       upwp = -elq*sm(k)*dudz
       vpwp = -elq*sm(k)*dvdz
       Tpwp = -elq*sh(k)*dTdz
       Tpwp = SIGN(MAX(ABS(Tpwp),1.E-6),Tpwp)

       IF ( k .EQ. kts+1 ) THEN
          qWT1D(kts)=0.
          q3sq_old =0.
          qWTP_old =0.
          !**  Limitation on q, instead of L/q  **
          dlsq1 = MAX(el(kts)**2,1.0)
          IF ( q3sq_old/dlsq1 .LT. -gh(k) ) q3sq_old = -dlsq1*gh(k)
       ENDIF

       !!!Vertical Transport Term
       qWTP_new = elq*Sqfac*sm(k)*(q3sq - q3sq_old)/dzk
       qWT1D(k) = 0.5*(qWTP_new - qWTP_old)/dzk
       qWTP_old = elq*Sqfac*sm(k)*(q3sq - q3sq_old)/dzk
       q3sq_old = q3sq

       !!!Shear Term
       !!!qSHEAR1D(k)=-(upwp*dudz + vpwp*dvdz)
       qSHEAR1D(k) = elq*sm(k)*gm(k)

       !!!Buoyancy Term    
       !!!qBUOY1D(k)=g*Tpwp/thl(k)
       !qBUOY1D(k)= elq*(sh(k)*gh(k) + gamv)
       qBUOY1D(k) = elq*(sh(k)*(-dTdz*g/thl(k)) + gamv)

       !!!Dissipation Term
       qDISS1D(k) = (q3sq**(3./2.))/(b1*MAX(el(k),1.))
    ENDIF

    END DO
!

    dfm(kts) = 0.0
    dfh(kts) = 0.0
    dfq(kts) = 0.0
    tcd(kts) = 0.0
    qcd(kts) = 0.0

    tcd(kte) = 0.0
    qcd(kte) = 0.0

!
    DO k = kts,kte-1
       dzk = dz(k)
       tcd(k) = ( tcd(k+1)-tcd(k) )/( dzk )
       qcd(k) = ( qcd(k+1)-qcd(k) )/( dzk )
    END DO
!

   IF ( bl_mynn_tkebudget == 1) THEN
      !JOE-TKE BUDGET
      qWT1D(kts)=0.
      qSHEAR1D(kts)=qSHEAR1D(kts+1)
      qBUOY1D(kts)=qBUOY1D(kts+1)
      qDISS1D(kts)=qDISS1D(kts+1)
   ENDIF

!    RETURN


  END SUBROUTINE mym_turbulence

! ==================================================================
!     SUBROUTINE  mym_predict:
!
!     Input variables:    see subroutine mym_initialize and turbulence
!       qke(nx,nz,ny) : qke at (n)th time level
!       tsq, ...cov     : ditto
!
!     Output variables:
!       qke(nx,nz,ny) : qke at (n+1)th time level
!       tsq, ...cov     : ditto
!
!     Work arrays:
!       qkw(nx,nz,ny)   : q at the center of the grid boxes        (m/s)
!       bp (nx,nz,ny)   : = 1/2*F,     see below
!       rp (nx,nz,ny)   : = P-1/2*F*Q, see below
!
!     # The equation for a turbulent quantity Q can be expressed as
!          dQ/dt + Ah + Av = Dh + Dv + P - F*Q,                      (1)
!       where A is the advection, D the diffusion, P the production,
!       F*Q the dissipation and h and v denote horizontal and vertical,
!       respectively. If Q is q^2, F is 2q/B_1L.
!       Using the Crank-Nicholson scheme for Av, Dv and F*Q, a finite
!       difference equation is written as
!          Q{n+1} - Q{n} = dt  *( Dh{n}   - Ah{n}   + P{n} )
!                        + dt/2*( Dv{n}   - Av{n}   - F*Q{n}   )
!                        + dt/2*( Dv{n+1} - Av{n+1} - F*Q{n+1} ),    (2)
!       where n denotes the time level.
!       When the advection and diffusion terms are discretized as
!          dt/2*( Dv - Av ) = a(k)Q(k+1) - b(k)Q(k) + c(k)Q(k-1),    (3)
!       Eq.(2) can be rewritten as
!          - a(k)Q(k+1) + [ 1 + b(k) + dt/2*F ]Q(k) - c(k)Q(k-1)
!                 = Q{n} + dt  *( Dh{n}   - Ah{n}   + P{n} )
!                        + dt/2*( Dv{n}   - Av{n}   - F*Q{n}   ),    (4)
!       where Q on the left-hand side is at (n+1)th time level.
!
!       In this subroutine, a(k), b(k) and c(k) are obtained from
!       subprogram coefvu and are passed to subprogram tinteg via
!       common. 1/2*F and P-1/2*F*Q are stored in bp and rp,
!       respectively. Subprogram tinteg solves Eq.(4).
!
!       Modify this subroutine according to your numerical integration
!       scheme (program).
!
!-------------------------------------------------------------------
  SUBROUTINE  mym_predict (kts,kte,&
       &            levflag,  &
       &            delt,&
       &            dz, &
       &            ust, flt, flq, pmz, phh, &
       &            el, dfq, &
       &            pdk, pdt, pdq, pdc,&
       &            qke, tsq, qsq, cov, &
       &            s_aw,s_awqke,sd_aw,sd_awqke,bl_mynn_edmf_tke &
       &)

!-------------------------------------------------------------------
    INTEGER, INTENT(IN) :: kts,kte    



    INTEGER, INTENT(IN) :: levflag
    INTEGER, INTENT(IN) :: bl_mynn_edmf_tke
    REAL, INTENT(IN)    :: delt
    REAL, DIMENSION(kts:kte), INTENT(IN) :: dz, dfq,el
    REAL, DIMENSION(kts:kte), INTENT(INOUT) :: pdk, pdt, pdq, pdc
    REAL, INTENT(IN)    ::  flt, flq, ust, pmz, phh
    REAL, DIMENSION(kts:kte), INTENT(INOUT) :: qke,tsq, qsq, cov
! WA 8/3/15
    REAL, DIMENSION(kts:kte+1), INTENT(INOUT) :: s_awqke,s_aw,sd_aw,sd_awqke

    INTEGER :: k,nz
    REAL, DIMENSION(kts:kte) :: qkw, bp, rp, df3q
    REAL :: vkz,pdk1,phm,pdt1,pdq1,pdc1,b1l,b2l,onoff
    REAL, DIMENSION(kts:kte) :: dtz
    REAL, DIMENSION(kts:kte) :: a,b,c,d,x

    nz=kte

    ! REGULATE THE MOMENTUM MIXING FROM THE MASS-FLUX SCHEME (on or off)
    IF (bl_mynn_edmf_tke == 0) THEN
       onoff=0.0
    ELSE
       onoff=1.0
    ENDIF

!   **  Strictly, vkz*h(i,j) -> vk*( 0.5*dz(1)*h(i,j)+z0 )  **
    vkz = vk*0.5*dz(kts)
!
!   **  dfq for the TKE is 3.0*dfm.  **
!
    DO k = kts,kte
!!       qke(k) = MAX(qke(k), 0.0)
       qkw(k) = SQRT( MAX( qke(k), 0.0 ) )
       df3q(k)=Sqfac*dfq(k)
       dtz(k)=delt/dz(k)
    END DO
!
    pdk1 = 2.0*ust**3*pmz/( vkz )
    phm  = 2.0/ust   *phh/( vkz )
    pdt1 = phm*flt**2
    pdq1 = phm*flq**2
    pdc1 = phm*flt*flq
!
!   **  pdk(i,j,1)+pdk(i,j,2) corresponds to pdk1.  **
    pdk(kts) = pdk1 -pdk(kts+1)

!!    pdt(kts) = pdt1 -pdt(kts+1)
!!    pdq(kts) = pdq1 -pdq(kts+1)
!!    pdc(kts) = pdc1 -pdc(kts+1)
    pdt(kts) = pdt(kts+1)
    pdq(kts) = pdq(kts+1)
    pdc(kts) = pdc(kts+1)
!
!   **  Prediction of twice the turbulent kinetic energy  **
!!    DO k = kts+1,kte-1
    DO k = kts,kte-1
       b1l = b1*0.5*( el(k+1)+el(k) )
       bp(k) = 2.*qkw(k) / b1l
       rp(k) = pdk(k+1) + pdk(k)
    END DO

!!    a(1)=0.
!!    b(1)=1.
!!    c(1)=-1.
!!    d(1)=0.

! Since df3q(kts)=0.0, a(1)=0.0 and b(1)=1.+dtz(k)*df3q(k+1)+bp(k)*delt.
    DO k=kts,kte-1
!       a(k-kts+1)=-dtz(k)*df3q(k)
!       b(k-kts+1)=1.+dtz(k)*(df3q(k)+df3q(k+1))+bp(k)*delt
!       c(k-kts+1)=-dtz(k)*df3q(k+1)
!       d(k-kts+1)=rp(k)*delt + qke(k)
! WA 8/3/15 add EDMF contribution
       a(k-kts+1)=-dtz(k)*df3q(k) + 0.5*dtz(k)*s_aw(k)*onoff
       b(k-kts+1)=1. + dtz(k)*(df3q(k)+df3q(k+1)) &
                     + 0.5*dtz(k)*(s_aw(k)-s_aw(k+1))*onoff + bp(k)*delt
       c(k-kts+1)=-dtz(k)*df3q(k+1) - 0.5*dtz(k)*s_aw(k+1)*onoff
       d(k-kts+1)=rp(k)*delt + qke(k) + dtz(k)*(s_awqke(k)-s_awqke(k+1))*onoff
    ENDDO

!!    DO k=kts+1,kte-1
!!       a(k-kts+1)=-dtz(k)*df3q(k)
!!       b(k-kts+1)=1.+dtz(k)*(df3q(k)+df3q(k+1))
!!       c(k-kts+1)=-dtz(k)*df3q(k+1)
!!       d(k-kts+1)=rp(k)*delt + qke(k) - qke(k)*bp(k)*delt
!!    ENDDO

    a(nz)=-1. !0.
    b(nz)=1.
    c(nz)=0.
    d(nz)=0.

!    CALL tridiag(nz,a,b,c,d)
    CALL tridiag2(nz,a,b,c,d,x)

    DO k=kts,kte
!       qke(k)=max(d(k-kts+1), 1.e-4)
       qke(k)=max(x(k), 1.e-4)
    ENDDO
      
!      
! level 3     
!
    IF ( levflag .EQ. 3 ) THEN
!
!  Modified: Dec/22/2005, from here
!   **  dfq for the scalar variance is 1.0*dfm.  **
!       CALL coefvu ( dfq, 1.0 ) make change here 
!  Modified: Dec/22/2005, up to here
!
!   **  Prediction of the temperature variance  **
!!       DO k = kts+1,kte-1
       DO k = kts,kte-1
          b2l = b2*0.5*( el(k+1)+el(k) )
          bp(k) = 2.*qkw(k) / b2l
          rp(k) = pdt(k+1) + pdt(k) 
       END DO
       
!zero gradient for tsq at bottom and top
       
!!       a(1)=0.
!!       b(1)=1.
!!       c(1)=-1.
!!       d(1)=0.

! Since dfq(kts)=0.0, a(1)=0.0 and b(1)=1.+dtz(k)*dfq(k+1)+bp(k)*delt.
       DO k=kts,kte-1
          a(k-kts+1)=-dtz(k)*dfq(k)
          b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))+bp(k)*delt
          c(k-kts+1)=-dtz(k)*dfq(k+1)
          d(k-kts+1)=rp(k)*delt + tsq(k)
       ENDDO

!!       DO k=kts+1,kte-1
!!          a(k-kts+1)=-dtz(k)*dfq(k)
!!          b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))
!!          c(k-kts+1)=-dtz(k)*dfq(k+1)
!!          d(k-kts+1)=rp(k)*delt + tsq(k) - tsq(k)*bp(k)*delt
!!       ENDDO

       a(nz)=-1. !0.
       b(nz)=1.
       c(nz)=0.
       d(nz)=0.

!       CALL tridiag(nz,a,b,c,d)
    CALL tridiag2(nz,a,b,c,d,x)
       
       DO k=kts,kte
!          tsq(k)=d(k-kts+1)
           tsq(k)=x(k)
       ENDDO
       
!   **  Prediction of the moisture variance  **
!!       DO k = kts+1,kte-1
       DO k = kts,kte-1
          b2l = b2*0.5*( el(k+1)+el(k) )
          bp(k) = 2.*qkw(k) / b2l
          rp(k) = pdq(k+1) +pdq(k) 
       END DO
       
!zero gradient for qsq at bottom and top
       
!!       a(1)=0.
!!       b(1)=1.
!!       c(1)=-1.
!!       d(1)=0.

! Since dfq(kts)=0.0, a(1)=0.0 and b(1)=1.+dtz(k)*dfq(k+1)+bp(k)*delt.
       DO k=kts,kte-1
          a(k-kts+1)=-dtz(k)*dfq(k)
          b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))+bp(k)*delt
          c(k-kts+1)=-dtz(k)*dfq(k+1)
          d(k-kts+1)=rp(k)*delt + qsq(k)
       ENDDO

!!       DO k=kts+1,kte-1
!!          a(k-kts+1)=-dtz(k)*dfq(k)
!!          b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))
!!          c(k-kts+1)=-dtz(k)*dfq(k+1)
!!          d(k-kts+1)=rp(k)*delt + qsq(k) -qsq(k)*bp(k)*delt
!!       ENDDO

       a(nz)=-1. !0.
       b(nz)=1.
       c(nz)=0.
       d(nz)=0.
       
!       CALL tridiag(nz,a,b,c,d)
       CALL tridiag2(nz,a,b,c,d,x)

       DO k=kts,kte
!          qsq(k)=d(k-kts+1)
           qsq(k)=x(k)
       ENDDO
       
!   **  Prediction of the temperature-moisture covariance  **
!!       DO k = kts+1,kte-1
       DO k = kts,kte-1
          b2l = b2*0.5*( el(k+1)+el(k) )
          bp(k) = 2.*qkw(k) / b2l
          rp(k) = pdc(k+1) + pdc(k) 
       END DO
       
!zero gradient for tqcov at bottom and top
       
!!       a(1)=0.
!!       b(1)=1.
!!       c(1)=-1.
!!       d(1)=0.

! Since dfq(kts)=0.0, a(1)=0.0 and b(1)=1.+dtz(k)*dfq(k+1)+bp(k)*delt.
       DO k=kts,kte-1
          a(k-kts+1)=-dtz(k)*dfq(k)
          b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))+bp(k)*delt
          c(k-kts+1)=-dtz(k)*dfq(k+1)
          d(k-kts+1)=rp(k)*delt + cov(k)
       ENDDO

!!       DO k=kts+1,kte-1
!!          a(k-kts+1)=-dtz(k)*dfq(k)
!!          b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))
!!          c(k-kts+1)=-dtz(k)*dfq(k+1)
!!          d(k-kts+1)=rp(k)*delt + cov(k) - cov(k)*bp(k)*delt
!!       ENDDO

       a(nz)=-1. !0.
       b(nz)=1.
       c(nz)=0.
       d(nz)=0.

!       CALL tridiag(nz,a,b,c,d)
    CALL tridiag2(nz,a,b,c,d,x)
       
       DO k=kts,kte
!          cov(k)=d(k-kts+1)
          cov(k)=x(k)
       ENDDO
       
    ELSE
!!       DO k = kts+1,kte-1
       DO k = kts,kte-1
          IF ( qkw(k) .LE. 0.0 ) THEN
             b2l = 0.0
          ELSE
             b2l = b2*0.25*( el(k+1)+el(k) )/qkw(k)
          END IF
!
          tsq(k) = b2l*( pdt(k+1)+pdt(k) )
          qsq(k) = b2l*( pdq(k+1)+pdq(k) )
          cov(k) = b2l*( pdc(k+1)+pdc(k) )
       END DO
       
!!       tsq(kts)=tsq(kts+1)
!!       qsq(kts)=qsq(kts+1)
!!       cov(kts)=cov(kts+1)

       tsq(kte)=tsq(kte-1)
       qsq(kte)=qsq(kte-1)
       cov(kte)=cov(kte-1)
      
    END IF

!
! end of level 3
!

  END SUBROUTINE mym_predict
  
  
  ! computes vt and vq for MYNN buoyancy equation from ql and cc

  SUBROUTINE  mym_vtvq (kts,kte,        &
    &            th,ql,thl,qw,cld,       &
    &            p,exner,liquid_frac,   &
    &            Vt, Vq)

!-------------------------------------------------------------------

    INTEGER, INTENT(IN)   :: kts,kte

    REAL, DIMENSION(kts:kte), INTENT(IN) :: p,exner, thl, qw, &
         & th,liquid_frac,ql,cld


    REAL, DIMENSION(kts:kte), INTENT(INOUT) :: vt,vq

    REAL :: qt,alp,bet,rac,cctq1,sgm,t,esat,lhb,qsl,dqsl,q2p,eq1,cct,q1
    INTEGER :: k

 
DO k = kts,kte
      qt   = 1.0 +p608*qw(k) -(1.+p608)*ql(k) ! Eq. B8, first three terms
   

    IF (cld(k) .eq. 0.) THEN
    ! simplest case  
      bet=0.
      rac=0.
   ELSE
  
    ! cc=1 is not compatible with the PDF approach
     cct=min(cld(k),0.98)
  
     ! compute q1 and sigma from CC and ql
      ! erfi is not fortran intrinsic function (we define it below) 
     q1=erfi(2.*cct-1.)/rr2
     sgm=0.5*ql(k)/(cct*q1+rrp*exp(-0.5*q1**2))
!    print *,'q1,sgm',q1,sgm
     sgm = max(min(sgm,1.0e-3),1.0e-6) 
   
   
     t  = th(k)*exner(k)
    !SATURATED VAPOR PRESSURE
    esat = esatLF_blend(t,liquid_frac(k))
    lhb=xlLF_blend(t,liquid_frac(k))
    !SATURATED SPECIFIC HUMIDITY
    qsl=ep_2*esat/(p(k)-ep_3*esat)
    !dqw/dT: Clausius-Clapeyron
    dqsl = qsl*ep_2*ev/( rd*t**2 )

    alp = 1.0/( 1.0+dqsl*lhb/cp ) ! a; Eq. B5
    bet = dqsl*exner(k)           ! b; Eq. B5

   q2p=lhb/(cp*exner(k))

    eq1  = rrp*EXP( -0.5*q1**2)  ! rrp=1/(2*pi)
    rac=(cld(k)-0.5*ql(k)/sgm*eq1)*alp*( q2p*qt-(1.+p608)*th(k) ) ! ~R*a*c
  
  endif
  
  vt(k) =   qt-1.0 -rac*bet ! beta_theta-1: Eq. B8
  vq(k) = p608*th(k)-tv0 +rac ! Eq. beta_qt-tv; Eq. B9


!  print *,'cc,vt,vq',cld(k),vt(k)+1,vq(k)+tv0
  
ENDDO  
  
  

END SUBROUTINE mym_vtvq

 real function erfi(x)
    real, intent(in) :: x
! Taylor expansion   
  erfi=0.8862*(x+0.2618*x**3+0.1439*x**5)
! higher order expnsion  
! erfi=0.8862*(x+0.2618*x**3+0.1439*x**5+0.0977*x**7+0.0733*x**9)  
 end function erfi



! ==================================================================
!     SUBROUTINE  mym_condensation:
!
!     Input variables:    see subroutine mym_initialize and turbulence
!       exner(nz)    : Perturbation of the Exner function    (J/kg K)
!                         defined on the walls of the grid boxes
!                         This is usually computed by integrating
!                         d(pi)/dz = h*g*tv/tref**2
!                         from the upper boundary, where tv is the
!                         virtual potential temperature minus tref.
!
!     Output variables:   see subroutine mym_initialize
!       cld(nx,nz,ny)   : Cloud fraction
!
!     Work arrays:
!       qmq(nx,nz,ny)   : Q_w-Q_{sl}, where Q_{sl} is the saturation
!                         specific humidity at T=Tl
!       alp(nx,nz,ny)   : Functions in the condensation process
!       bet(nx,nz,ny)   : ditto
!       sgm(nx,nz,ny)   : Combined standard deviation sigma_s
!                         multiplied by 2/alp
!
!     # qmq, alp, bet and sgm are allowed to share storage units with
!       any four of other work arrays for saving memory.
!
!     # Results are sensitive particularly to values of cp and rd.
!       Set these values to those adopted by you.
!
!-------------------------------------------------------------------



  SUBROUTINE  mym_condensation (kts,kte,  &
    &            dx, dz,                  &
    &            thl, qw,                 &
    &            p,exner,                 &
    &            tsq, qsq, cov,           &
    &            Sh, el, liquid_frac,     &
    &            ql, cld,    &
    &            PBLH1,HFX1,              &
    &            Vt, Vq, th, sgm)

!-------------------------------------------------------------------

    INTEGER, INTENT(IN)   :: kts,kte


    REAL, INTENT(IN)      :: dx,PBLH1,HFX1
    REAL, DIMENSION(kts:kte), INTENT(IN) :: dz
    REAL, DIMENSION(kts:kte), INTENT(IN) :: p,exner, thl, qw, &
         &tsq, qsq, cov, th,liquid_frac
    REAL, DIMENSION(kts:kte), INTENT(IN) :: Sh,el
    REAL, DIMENSION(kts:kte), INTENT(OUT) :: ql,cld

    REAL, DIMENSION(kts:kte), INTENT(INOUT) :: vt,vq,sgm

!    REAL, DIMENSION(kts:kte) :: qmq,a,b
!    DOUBLE PRECISION :: t3sq, r3sq, c3sq
    REAL :: alp,bet,qi,lhb,q1


    REAL :: qsl,esat,qsat,tlk,qsat_tl,dqsl,eq1,qll,&
         &q2p,pt,rac,qt,t,xl,rsl,cpm,cdhdz,Fng,qww,alpha,beta,bb,ls_min,ls,wt
    INTEGER :: i,j,k

    ! REAL :: erf,lhb

    !JOE: NEW VARIABLES FOR ALTERNATE SIGMA
    REAL::dth,dtl,dqw,dzk


    !JOE: variables for BL clouds
   ! REAL::zagl,cld9,damp,edown,RHcrit,RHmean,RHsum,RHnum,Hshcu,PBLH2,ql_limit
   ! REAL, PARAMETER :: Hfac = 3.0     !cloud depth factor for HFX (m^3/W)
   ! REAL, PARAMETER :: HFXmin = 50.0  !min W/m^2 for BL clouds
   ! REAL            :: RH_00L, RH_00O, phi_dz, lfac
   ! REAL, PARAMETER :: cdz = 2.0
   ! REAL, PARAMETER :: mdz = 1.5

! ORIGINAL MYNN PARTIAL-CONDENSATION SCHEME
! OR KUWANO ET AL.


!    zagl = 0.


DO k = kts,kte-1
   t  = th(k)*exner(k)
   !SATURATED VAPOR PRESSURE
   esat = esatLF_blend(t,liquid_frac(k))
   lhb=xlLF_blend(t,liquid_frac(k))
   !SATURATED SPECIFIC HUMIDITY
   qsl=ep_2*esat/(p(k)-ep_3*esat)
   !dqw/dT: Clausius-Clapeyron
   dqsl = qsl*ep_2*ev/( rd*t**2 )
   !RH (0 to 1.0)
   !RH(k)=MAX(MIN(1.0,qw(k)/MAX(1.E-8,qsl)),0.001)

   alp = 1.0/( 1.0+dqsl*lhb/cp ) ! a; Eq. B5
   bet = dqsl*exner(k)           ! b; Eq. B5


! sigma_s: EQ. B6
!   if (k .eq. kts) then 
!   	 dzk = dz(k)
!   else
!	 dzk = 0.5*( dz(k) + dz(k-1) )
!   end if
!   dth = 0.5*(thl(k+1)+thl(k)) - 0.5*(thl(k)+thl(MAX(k-1,kts)))
!   dqw = 0.5*(qw(k+1) + qw(k)) - 0.5*(qw(k) + qw(MAX(k-1,kts)))
    
   if (k .eq. kts) then
     dzk=dz(k)
     dth=thl(k+1)-thl(k)
     dqw=qw(k+1)-qw(k)
   else
     dzk=dz(k-1)
     dth=thl(k)-thl(k-1)
     dqw=qw(k)-qw(k-1)     
   endif


   sgm(k) = SQRT( MAX( (alp**2 * MAX(el(k)**2,0.1) * &
					 b2 * MAX(Sh(k),0.03))/4. * &
			  (dqw/dzk - bet*(dth/dzk ))**2 , 1.0e-12) ) 
  ! sgm(k) = min(sgm(k),1.0e-3) 
   
  ! sgm(k)=1.e-5  ! yhc_mynn add, 2021-04-12

  !sgm(k)=100.*alp*abs(dqw/dzk)
  sgm(k)=sgm_factor*alp*abs(dqw/dzk)    ! yhc_mynn, mkae sgm_factor
  sgm(k)=max(min(sgm(k),1.e-3),1.e-6)  
 
   
   q1   = alp*(qw(k)-qsl) / (2.*sgm(k)) ! Q1; Eq. B4
   
   cld(k) = 0.5*( 1.0+erf( q1*rr2 ) ) ! Eq. B2
   
   eq1  = rrp*EXP( -0.5*q1**2)  ! rrp=1/(2*pi)
   ql (k) = 2.*sgm(k)*MAX( cld(k)*q1 + eq1, 0.0 ) ! Eq. B1
 
! sanity check

   if (ql(k) > qw(k)) then
      ql(k)=qw(k)
    endif

   if ( ql(k) < 1.e-8 ) then
      ql(k)=0.
      cld(k)=0.
   endif
  

   if (cld(k)< 0.01) then
     cld(k)=0.
     ql(k)=0.
   endif


   if ((cld(k)==0.) .or. (ql(k)==0.)) then
    cld(k)=0.
    ql(k)=0.
  endif 



 ! qll=MAX( cld(k)*q1 + eq1, 0.0 )
 
   q2p = lhb/(cp*exner(k)) ! L/(cp*Pi)
   pt = thl(k) +q2p*ql(k) ! potential temp

   !qt is a THETA-V CONVERSION FOR TOTAL WATER (i.e., THETA-V = qt*THETA)
   qt   = 1.0 +p608*qw(k) -(1.+p608)*ql(k) ! Eq. B8, first three terms

   rac=(cld(k)-0.5*ql(k)/sgm(k)*eq1)*alp*( q2p*qt-(1.+p608)*pt ) ! ~R*a*c
  ! print *,'qt,pt,r,a,c',qt,pt,(cld(k)-0.5*ql(k)/sgm(k)*eq1),alp,( q2p*qt-(1.+p608)*pt )
  ! print *,'a',alp
  ! print *,'c',( q2p*qt-(1.+p608)*pt )
   !BUOYANCY FACTORS: wherever vt and vq are used, there is a
   !"+1" and "+tv0", respectively, so these are subtracted out here.
   !vt is unitless and vq has units of K.
   vt(k) =      qt-1.0 -rac*bet ! beta_theta-1: Eq. B8
   vq(k) = p608*pt-tv0 +rac ! Eq. beta_qt-tv;

      

!print*,'k,ql,cld',k,ql(k),cld(k)  ! yhc_mynn

    
END DO
    
    cld(kte) = cld(kte-1)
    ql(kte) = ql(kte-1)
    vt(kte) = vt(kte-1)
    vq(kte) = vq(kte-1)

!  print *,'beta_theta',vt+1.
!  print *,'beta_qt',vq+tv0
!  print *,'cld',cld
!  print *,'ql',ql
!  print *,'rac',rac

END SUBROUTINE mym_condensation

! ==================================================================
  SUBROUTINE mynn_tendencies(kts,kte,      &
       &levflag,grav_settling,rho,rhoh,    &
       &delt,dz,                           &
       &u,v,th,tk,qv,qc,qi,qni,qnc,        &
       &p,exner,                           &
       &thl,sqv,sqc,sqi,sqw,               &
       &ust,flt,flq,flqv,flqc,wspd,qcg,    &
       &uoce,voce,                         &
       &tsq,qsq,cov,                       &
       &tcd,qcd,                           &
       &dfm1,dfh1,dfq,                       &
       &Du,Dv,Dthl,Dsqw,Dqni,        &!Dqnc,   &
       &vdfg1,                             &
       &s_aw1,s_awthl1,s_awqt1, &
       &s_awu1,s_awv1,                       &
       &sd_aw1,sd_awthl1,sd_awqt1, &
       &sd_awu1,sd_awv1,                       &
       &ztop_shallow,ktop_shallow,         &
       &bl_mynn_cloudmix,                  &
       &bl_mynn_mixqt,                     &
       &bl_mynn_edmf,                      &
       &bl_mynn_edmf_dd,                   &
       &bl_mynn_edmf_mom,                  &
       &liquid_frac,                       &
       &dx,PBLH,HFX,Sh,el,Vt,Vq)

!-------------------------------------------------------------------
    INTEGER, INTENT(in) :: kts,kte



    INTEGER, INTENT(in) :: grav_settling,levflag
    INTEGER, INTENT(in) :: bl_mynn_cloudmix,bl_mynn_mixqt,&
                           bl_mynn_edmf,bl_mynn_edmf_mom,bl_mynn_edmf_dd
 !   LOGICAL, INTENT(IN) :: FLAG_QI,FLAG_QNI,FLAG_QC,FLAG_QNC

!! grav_settling = 1 or 2 for gravitational settling of droplets
!! grav_settling = 0 otherwise
! thl - liquid water potential temperature
! qw - total water
! dfm,dfh,dfq - as above
! flt - surface flux of thl
! flq - surface flux of qw

   REAL,DIMENSION(kts:kte), INTENT(IN) :: liquid_frac

    REAL,DIMENSION(kts:kte+1), INTENT(in) :: rhoh,s_aw1,s_awthl1,s_awqt1,&
                             s_awu1,s_awv1,&
                             sd_aw1,sd_awthl1,sd_awqt1,&
                             sd_awu1,sd_awv1
    REAL, DIMENSION(kts:kte), INTENT(in) :: u,v,th,tk,qv,qc,qi,qni,qnc,&
         &p,exner,dfq,dz,tsq,qsq,cov,tcd,qcd,rho
    REAL, DIMENSION(kts:kte), INTENT(inout) :: thl,sqw
    REAL,DIMENSION(kts:kte), INTENT(in) :: sqv,sqc,sqi,&
         &dfm1,dfh1
    REAL, DIMENSION(kts:kte), INTENT(inout) :: du,dv,dthl,dsqw,dqni !,dqnc
    REAL, INTENT(IN) :: delt,ust,flt,flq,flqv,flqc,wspd,uoce,voce,qcg,&
         ztop_shallow
    INTEGER, INTENT(IN) :: ktop_shallow

!    REAL, INTENT(IN) :: delt,ust,flt,flq,qcg,&
!         &gradu_top,gradv_top,gradth_top,gradqv_top

!local vars

    REAL, DIMENSION(kts:kte) :: dtz,vt,vq,dfhc,dfmc !Kh for clouds (Pr < 2)
    REAL, DIMENSION(kts:kte) :: sqv2,sqc2,sqi2,sqw2,qni2 !,qnc2 !AFTER MIXING
    REAL, DIMENSION(kts:kte) :: zfac,plumeKh,th_temp
    REAl, DIMENSION(kts:kte) :: upcont,dncont ! updraft/downdraft contribution to fluxes for explicit calculation
    REAL, DIMENSION(1:kte-kts+1) :: a,b,c,d,x

    REAL :: rhs,gfluxm,gfluxp,dztop,maxdfh,mindfh,maxcf,maxKh,zw
    REAL :: grav_settling2,vdfg1    !Katata-fogdes
    REAL :: t,esat,qsl,onoff,ustovwsp
    INTEGER :: k,kk,nz,itr
    REAL, INTENT(IN)      :: dx,PBLH,HFX
!    REAL, DIMENSION(kts:kte), INTENT(INOUT) :: qc_bl1D,cldfra_bl1d,sgm
    REAL, DIMENSION(kts:kte), INTENT(IN) :: Sh,el
!yhc move to namelist    REAL upwind 
!yhc move to namelist    LOGICAl expmf 
!yhc move to namelist    
!yhc move to namelist    ! upwind=1. ... use upwind approximation for mass-flux calculation
!yhc move to namelist    ! upwind=0.5 ... use centered difference for mass-flux calculation
!yhc move to namelist    upwind=1. 
!yhc move to namelist
!yhc move to namelist   ! expmf=.true.  ... explicit mass-flux
!yhc move to namelist   ! expmf =.false. .... implicit mass-flux
!yhc move to namelist    expmf=.true. 


    REAL,DIMENSION(kts:kte+1) :: s_aw,s_awthl,s_awqt,&
                             s_awu,s_awv,&
                             sd_aw,sd_awthl,sd_awqt,&
                             sd_awu,sd_awv
                             
                             
   REAL,DIMENSION(kts:kte) :: dfm,dfh 
   REAL(kind=selected_real_kind(14)) :: tt                          

   REAL :: test

    nz=kte-kts+1

    dztop=.5*(dz(kte)+dz(kte-1))


   ! USTAR/WSPD make sure it does not blow up when WSPD = 0.
   IF (wspd .le. 0.) THEN
     ustovwsp=0.
   ELSE 
       ustovwsp=ust/wspd
   ENDIF

    ! REGULATE THE MOMENTUM MIXING FROM THE MASS-FLUX SCHEME (on or off)
    ! Note that s_awu and s_awv already come in as 0.0 if bl_mynn_edmf_mom == 0, so
    ! we only need to zero-out the MF term
    IF (bl_mynn_edmf_mom == 0) THEN
       onoff=0.0
    ELSE
       onoff=1.0
    ENDIF

    !set up values for background diffusivity when MF scheme is active
!    maxdfh=maxval(dfh(1:14))
!    maxcf=maxval(cldfra_bl1D(kts:MAX(ktop_shallow,14)))
    !allow maxKh to vary according to cloud fraction in lowest ~2 km
!    maxKh = 1.*(1.-MIN(MAX(maxcf-0.5,0.0)/0.25, 0.9))
!    mindfh=min(maxKh,maxdfh*0.01)

!    zw=0.


!
! we are solving equation d phi/dt=-1/rho*d/dz(rho <w'phi'>)
! instead of the one without density
! to do so we need to divide dtz with density and multiply K, a_i*w_i and a_i*w_i*psi_i with density 
! and also modify d(1)'s
   dtz(kts:kte)=delt/(dz(kts:kte)*rho(kts:kte))

   s_aw=rhoh*s_aw1
   s_awthl=rhoh*s_awthl1
   s_awqt=rhoh*s_awqt1
   s_awu=rhoh*s_awu1
   s_awv=rhoh*s_awv1
   
   sd_aw=rhoh*sd_aw1
   sd_awthl=rhoh*sd_awthl1
   sd_awqt=rhoh*sd_awqt1
   sd_awu=rhoh*sd_awu1
   sd_awv=rhoh*sd_awv1

   dfm(kts:kte)=dfm1(kts:kte)*rhoh(kts:kte)
   dfh(kts:kte)=dfh1(kts:kte)*rhoh(kts:kte)

!    DO k=kts,kte
!       dtz(k)=delt/(dz(k)*rho(k))
       !IF (dfm(k) > dfh(k)) THEN 
       !  !in stable regime only, limit Prandtl number to < 2 within clouds
       !  IF (qc(k) > 1.e-6 .OR. &
       !      qi(k) > 1.e-6 .OR. &
       !      cldfra_bl1D(k) > 0.05 ) THEN
       !      dfh(k)= MAX(dfh(k),dfm(k)*0.5)
       !  ENDIF
       !ENDIF
       !Add small minimum Km & Kh in MF updrafts is no stratus is in the column. 
       !Note that maxval of plumeKh is mindfh*0.15, with max at about 0.75*ztop_shallow
!        IF (ktop_shallow > 0) THEN
!           zfac(k) = min( max(1.-(zw/ztop_shallow), 0.01), 1.)
!           plumeKh(k)=mindfh*max((ztop_shallow-zw)/ztop_shallow,0.0)*(1.-zfac(k))**2
!           dfh(k)=MAX(mindfh,dfh(k))
!           dfm(k)=MAX(mindfh,dfm(k))
!       ENDIF
!       zw=zw+dz(k)
!    ENDDO

!!============================================
!! u
!!============================================

 IF(expmf) THEN


   DO k=kts+1,kte-1
    upcont(k)=onoff*(s_awu(k)-s_aw(k)*(u(k)*upwind+u(k-1)*(1.-upwind)))
    dncont(k)=onoff*(sd_awu(k)-sd_aw(k)*(u(k)*upwind+u(k-1)*(1.-upwind)))
   ENDDO
! no flux at the top of the atmosphere
    upcont(kte)=0. 
    dncont(kte)=0.  
! upcont(1) and dncont(1) are not used so they don't need to be set

    k=kts
    a(1)=0.
    b(1)=1. + dtz(k)*(dfm(k+1)+ust*ustovwsp*rhoh(1)) 
    c(1)=-dtz(k)*dfm(k+1) 
    d(1)=u(k) + dtz(k)*uoce*ust*ustovwsp*rhoh(1) -dtz(k)*(upcont(k+1)+dncont(k+1))

    DO k=kts+1,kte-1
       a(k)=   - dtz(k)*dfm(k)          
       b(k)=1. + dtz(k)*(dfm(k)+dfm(k+1)) 
       c(k)=   - dtz(k)*dfm(k+1)          
       d(k)=u(k) -dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
    ENDDO

  ELSE
   
    k=kts
    a(1)=0.
    b(1)=1. + dtz(k)*(dfm(k+1)+ust*ustovwsp*rhoh(1)) - 0.5*dtz(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*sd_aw(k+1)*onoff
    c(1)=-dtz(k)*dfm(k+1) - 0.5*dtz(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*sd_aw(k+1)*onoff
    d(1)=u(k) + dtz(k)*uoce*ust*ustovwsp*rhoh(1) - dtz(k)*s_awu(k+1)*onoff - dtz(k)*sd_awu(k+1)*onoff 

    DO k=kts+1,kte-1
       a(k)=   - dtz(k)*dfm(k)            + 0.5*dtz(k)*s_aw(k)*onoff + 0.5*dtz(k)*sd_aw(k)*onoff
       b(k)=1. + dtz(k)*(dfm(k)+dfm(k+1)) + 0.5*dtz(k)*(s_aw(k)-s_aw(k+1))*onoff + 0.5*dtz(k)*(sd_aw(k)-sd_aw(k+1))*onoff
       c(k)=   - dtz(k)*dfm(k+1)          - 0.5*dtz(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*sd_aw(k+1)*onoff
       d(k)=u(k) + dtz(k)*(s_awu(k)-s_awu(k+1))*onoff + dtz(k)*(sd_awu(k)-sd_awu(k+1))*onoff
    ENDDO
  
 ENDIF 


!! no flux at the top
!    a(nz)=-1.
!    b(nz)=1.
!    c(nz)=0.
!    d(nz)=0.

!! specified gradient at the top 
!    a(nz)=-1.
!    b(nz)=1.
!    c(nz)=0.
!    d(nz)=gradu_top*dztop

!! prescribed value
    a(nz)=0
    b(nz)=1.
    c(nz)=0.
    d(nz)=u(kte)

!    CALL tridiag(nz,a,b,c,d)
    CALL tridiag2(nz,a,b,c,d,x)

    DO k=kts,kte
!       du(k)=(d(k-kts+1)-u(k))/delt
       du(k)=(x(k)-u(k))/delt
    ENDDO

!!============================================
!! v
!!============================================

 IF(expmf) THEN

   DO k=kts+1,kte-1
    upcont(k)=onoff*(s_awv(k)-s_aw(k)*(v(k)*upwind+v(k-1)*(1.-upwind)))
    dncont(k)=onoff*(sd_awv(k)-sd_aw(k)*(v(k)*upwind+v(k-1)*(1.-upwind)))
   ENDDO
    upcont(kte)=0. 
    dncont(kte)=0.

    k=kts
    a(1)=0.
    b(1)=1. + dtz(k)*(dfm(k+1)+ust*ustovwsp*rhoh(1)) 
    c(1)=   - dtz(k)*dfm(k+1)               
    d(1)=v(k) + dtz(k)*voce*ust*ustovwsp*rhoh(1) -dtz(k)*(upcont(k+1)+dncont(k+1))

    DO k=kts+1,kte-1
       a(k)=   - dtz(k)*dfm(k)           
       b(k)=1. + dtz(k)*(dfm(k)+dfm(k+1)) 
       c(k)=   - dtz(k)*dfm(k+1)         
       d(k)=v(k) -dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
    ENDDO

 ELSE 

    k=kts
    a(1)=0.
    b(1)=1. + dtz(k)*(dfm(k+1)+ust*ustovwsp*rhoh(1)) - 0.5*dtz(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*sd_aw(k+1)*onoff
    c(1)=   - dtz(k)*dfm(k+1)               - 0.5*dtz(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*sd_aw(k+1)*onoff
    d(1)=v(k) + dtz(k)*voce*ust*ustovwsp*rhoh(1) - dtz(k)*s_awv(k+1)*onoff

    DO k=kts+1,kte-1
       a(k)=   - dtz(k)*dfm(k)            + 0.5*dtz(k)*s_aw(k)*onoff + 0.5*dtz(k)*sd_aw(k)*onoff 
       b(k)=1. + dtz(k)*(dfm(k)+dfm(k+1)) + 0.5*dtz(k)*(s_aw(k)-s_aw(k+1))*onoff + 0.5*dtz(k)*(sd_aw(k)-sd_aw(k+1))*onoff
       c(k)=   - dtz(k)*dfm(k+1)          - 0.5*dtz(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*sd_aw(k+1)*onoff
       d(k)=v(k) + dtz(k)*(s_awv(k)-s_awv(k+1))*onoff + dtz(k)*(sd_awv(k)-sd_awv(k+1))*onoff
    ENDDO

  ENDIF


!! prescribed value
    a(nz)=0
    b(nz)=1.
    c(nz)=0.
    d(nz)=v(kte)

!    CALL tridiag(nz,a,b,c,d)
    CALL tridiag2(nz,a,b,c,d,x)

    DO k=kts,kte
!       dv(k)=(d(k-kts+1)-v(k))/delt
       dv(k)=(x(k)-v(k))/delt
    ENDDO

!!============================================
!! thl tendency
!! NOTE: currently, gravitational settling is removed
!!============================================
 
 
 ! explicit mass-flux
 
 !  upcont,dncont=sum a_i*w_i*(thl_i-<thl>) for updrafts and downdrafts, respectively 
 ! upwind=1. upwind approximation, upwind=0.5 centered difference	
 ! (note: we do not need surface conditions)
 !
 
 
  IF(expmf) THEN
  
   DO k=kts+1,kte-1
   upcont(k)=s_awthl(k)-s_aw(k)*(thl(k)*upwind+thl(k-1)*(1.-upwind))
   dncont(k)=sd_awthl(k)-sd_aw(k)*(thl(k)*upwind+thl(k-1)*(1.-upwind))
   ENDDO
   upcont(kte)=0. 
   dncont(kte)=0.

    k=kts

    a(k)=0.
    b(k)=1.+dtz(k)*dfh(k+1) 
    c(k)=  -dtz(k)*dfh(k+1) 
    d(k)=thl(k) + dtz(k)*flt*rhoh(1) -dtz(k)*(upcont(k+1)+dncont(k+1))

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*dfh(k)            
       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) 
       c(k)=  -dtz(k)*dfh(k+1)          
       d(k)=thl(k) -dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
    ENDDO

 ELSE
 
    k=kts
    a(k)=0.
    b(k)=1.+dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)
    c(k)=  -dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)
    d(k)=thl(k) + dtz(k)*flt*rhoh(1) -dtz(k)*s_awthl(kts+1) -dtz(k)*sd_awthl(kts+1)

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*dfh(k)            + 0.5*dtz(k)*s_aw(k) + 0.5*dtz(k)*sd_aw(k)
       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) + 0.5*dtz(k)*(s_aw(k)-s_aw(k+1)) + 0.5*dtz(k)*(sd_aw(k)-sd_aw(k+1))
       c(k)=  -dtz(k)*dfh(k+1)          - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)
       d(k)=thl(k) + dtz(k)*(s_awthl(k)-s_awthl(k+1)) + dtz(k)*(sd_awthl(k)-sd_awthl(k+1))
    ENDDO
 
  ENDIF



!! prescribed value
    a(nz)=0.
    b(nz)=1.
    c(nz)=0.
    d(nz)=thl(kte)

!    CALL tridiag(nz,a,b,c,d)
    CALL tridiag2(nz,a,b,c,d,x)

 !   IF(expmf) THEN
 !   x(kts)=X(kts)-dtz(k)*upcont(k+1)
 !   DO k=kts+1,kte-1
 !   x(k)=x(K) -dtz(k)*(upcont(k+1)-upcont(k))
 !   ENDDO
 !   ENDIF

    DO k=kts,kte
        dthl(k)=(x(k)-thl(k))/delt
        thl(k)=x(k)
    ENDDO


 !============================================
 ! MIX total water (sqw = sqc + sqv + sqi)
 ! NOTE: no total water tendency is output; instead, we must calculate
 !       the saturation specific humidity and then 
 !       subtract out the moisture excess (sqc & sqi)
 !============================================

 IF(expmf) THEN
   DO k=kts+1,kte-1
     upcont(k)=s_awqt(k)-s_aw(k)*(sqw(k)*upwind+sqw(k-1)*(1.-upwind))
     dncont(k)=sd_awqt(k)-sd_aw(k)*(sqw(k)*upwind+sqw(k-1)*(1.-upwind))
   ENDDO
    upcont(kte)=0. 
    dncont(kte)=0.
 
    k=kts
    a(k)=0.
    b(k)=1.+dtz(k)*dfh(k+1) 
    c(k)=  -dtz(k)*dfh(k+1)
    d(k)=sqw(k) + dtz(k)*flq*rhoh(1) - dtz(k)*(upcont(k+1)+dncont(k+1))


    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*dfh(k)           
       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) 
       c(k)=  -dtz(k)*dfh(k+1) 
       d(k)=sqw(k) -dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
    ENDDO
    
  ELSE  
    k=kts

    a(k)=0.
    b(k)=1.+dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)
    c(k)=  -dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)

    !rhs= qcd(k) !+ (gfluxp - gfluxm)/dz(k)& 

    d(k)=sqw(k) + dtz(k)*flq*rhoh(1)  - (dtz(k)*s_awqt(k+1) + dtz(k)*sd_awqt(k+1))

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*dfh(k)            + 0.5*dtz(k)*s_aw(k) + 0.5*dtz(k)*sd_aw(k)
       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) + 0.5*dtz(k)*(s_aw(k)-s_aw(k+1)) + 0.5*dtz(k)*(sd_aw(k)-sd_aw(k+1))
       c(k)=  -dtz(k)*dfh(k+1)          - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)

       d(k)=sqw(k) + dtz(k)*(s_awqt(k)-s_awqt(k+1)) + dtz(k)*(sd_awqt(k)-sd_awqt(k+1))
    ENDDO
  ENDIF   
     
!! prescribed value
    a(nz)=0.
    b(nz)=1.
    c(nz)=0.
    d(nz)=sqw(kte)

!    CALL tridiag(nz,a,b,c,d)
    CALL tridiag2(nz,a,b,c,d,x)
 !   IF(expmf) THEN
 !   x(kts)=X(kts)-dtz(k)*upcont(k+1)
 !   DO k=kts+1,kte-1
 !   x(k)=x(K) -dtz(k)*(upcont(k+1)-upcont(k))
 !   ENDDO
 !  ENDIF

    DO k=kts,kte
        dsqw(k)=(x(k)-sqw(k))/delt
         sqw(k)=x(k) 
   ENDDO

!<--- yhc, Kay add check total water conservation 
! rho <w'q_t'>|surf=int_z (rho * d q_t/dz) dz
!test=sum(dsqw*dz*rho) ! test is the RHS in the upper equation
!print *,'surface flux qt',flq*rhoh(1)
!print *,'zero q_t',test-flq*rhoh(1) ! difference between the LHS and RHS of  upper equation
!-->

!!============================================
!! cloud ice number concentration (qni)
!!============================================
! diasbled this since scalar_pblmix option can be invoked instead
!IF (bl_mynn_cloudmix > 0 .AND. FLAG_QNI) THEN

    qni2=qni

  END SUBROUTINE mynn_tendencies

! ==================================================================


! ==================================================================
  SUBROUTINE retrieve_exchange_coeffs(kts,kte,&
       &dfm,dfh,dfq,dz,&
       &K_m,K_h,K_q)

!-------------------------------------------------------------------

    INTEGER , INTENT(in) :: kts,kte

    REAL, DIMENSION(KtS:KtE), INTENT(in) :: dz,dfm,dfh,dfq

    REAL, DIMENSION(KtS:KtE), INTENT(out) :: &
         &K_m, K_h, K_q


    INTEGER :: k
    REAL :: dzk

    K_m(kts)=0.
    K_h(kts)=0.
    K_q(kts)=0.

    DO k=kts+1,kte
       dzk = 0.5  *( dz(k)+dz(k-1) )
       K_m(k)=dfm(k)*dzk
       K_h(k)=dfh(k)*dzk
       K_q(k)=Sqfac*dfq(k)*dzk
    ENDDO

  END SUBROUTINE retrieve_exchange_coeffs

! ==================================================================
  SUBROUTINE tridiag(n,a,b,c,d)

!! to solve system of linear eqs on tridiagonal matrix n times n
!! after Peaceman and Rachford, 1955
!! a,b,c,d - are vectors of order n 
!! a,b,c - are coefficients on the LHS
!! d - is initially RHS on the output becomes a solution vector
    
!-------------------------------------------------------------------

    INTEGER, INTENT(in):: n
    REAL, DIMENSION(n), INTENT(in) :: a,b
    REAL, DIMENSION(n), INTENT(inout) :: c,d
    
    INTEGER :: i
    REAL :: p
    REAL, DIMENSION(n) :: q
    
    c(n)=0.
    q(1)=-c(1)/b(1)
    d(1)=d(1)/b(1)
    
    DO i=2,n
       p=1./(b(i)+a(i)*q(i-1))
       q(i)=-c(i)*p
       d(i)=(d(i)-a(i)*d(i-1))*p
    ENDDO
    
    DO i=n-1,1,-1
       d(i)=d(i)+q(i)*d(i+1)
    ENDDO

  END SUBROUTINE tridiag

! ==================================================================
      subroutine tridiag2(n,a,b,c,d,x)
      implicit none
!      a - sub-diagonal (means it is the diagonal below the main diagonal)
!      b - the main diagonal
!      c - sup-diagonal (means it is the diagonal above the main diagonal)
!      d - right part
!      x - the answer
!      n - number of unknowns (levels)

        integer,intent(in) :: n
        real, dimension(n),intent(in) :: a,b,c,d
        real ,dimension(n),intent(out) :: x
        real ,dimension(n) :: cp,dp
        real :: m
        integer :: i

        ! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
        ! solve for vectors c-prime and d-prime
        do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
        enddo
        ! initialize x
        x(n) = dp(n)
        ! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
           x(i) = dp(i)-cp(i)*x(i+1)
        end do

    end subroutine tridiag2
! ==================================================================
  SUBROUTINE mynn_bl_driver(            &
       &initflag,grav_settling,         &
       &delt,dz,dx,znt,                 &
       &u,v,w,th,qv,ql,qi,cc,qni,qnc,      &
       &p,exner,rho,T3D,                &
       &xland,ts,qsfc,qcg,ps,           &
       &ust,ch,hfx,qfx,rmol,wspd,       &
       &uoce,voce,                      & !ocean current
       &vdfg,                           & !Katata-added for fog dep
       &Qke,                    &
       &Tsq,Qsq,Cov,                    &
       &RUBLTEN,RVBLTEN,RTHBLTEN,       &
       &RQVBLTEN,RQLBLTEN,RQIBLTEN,     &
       &RQNIBLTEN,                      &
       &RCCBLTEN, RTHLBLTEN, RQTBLTEN,  & ! yhc_mynn add
       &qa_before_mix, ql_before_mix, qi_before_mix, thl_before_mix, qt_before_mix, th_before_mix, &  ! yhc_mynn add
       &qa_after_mix, ql_after_mix, qi_after_mix, thl_after_mix, qt_after_mix, th_after_mix,       &  ! yhc_mynn add
       &qa_before_pdf, ql_before_pdf, qi_before_pdf,                                               &  ! yhc_mynn add
       &Q_ql,Q_qi,Q_a,                                                                             &  ! yhc_mynn add
       &Q_ql_adv,Q_qi_adv,Q_a_adv, Q_ql_eddy,Q_qi_eddy,Q_a_eddy, Q_ql_ent,Q_qi_ent,Q_a_ent, Q_ql_det,Q_qi_det,Q_a_det, Q_ql_sub,Q_qi_sub,Q_a_sub,       &  ! yhc_mynn_add
       &a_moist_half, mf_moist_half, qv_moist_half, a_moist_full, mf_moist_full, qv_moist_full, &  ! yhc 2021-09-08
       &a_dry_half, mf_dry_half, qv_dry_half, a_dry_full, mf_dry_full, qv_dry_full, &            ! yhc 2021-09-08
       &mf_all_half, mf_all_full, &                                                     ! yhc 2021-09-08
       &num_updraft, num_DET, num_nDET_pENT, num_nDET_zENT, &                                        ! yhc 2021-09-08
       &streams, &  ! yhc 2021-11-18
       &exch_h,exch_m,                  &
       &Pblh,kpbl,                      & 
       &el_pbl,                         &
       &dqke,qWT,qSHEAR,qBUOY,qDISS,    & !JOE-TKE BUDGET
       &bl_mynn_tkebudget,              &
       &bl_mynn_cloudpdf,Sh3D,          &
       &bl_mynn_mixlength,              &
       &icloud_bl,qc_bl,cldfra_bl,      &
       &bl_mynn_edmf,                   &
       &bl_mynn_edmf_dd,                   &
       &bl_mynn_edmf_mom,bl_mynn_edmf_tke, &
       &bl_mynn_edmf_part,bl_mynn_edmf_Lent,&
       &bl_mynn_cloudmix,bl_mynn_mixqt, &
       &edmf_a,edmf_w,edmf_qt,          &
       &edmf_thl,edmf_ent,edmf_det,edmf_qc,      &
       &edmf_debug1,edmf_debug2,        &
       &edmf_debug3,edmf_debug4,        &
       &edmf_a_dd,edmf_w_dd,edmf_qt_dd,    &
       &edmf_thl_dd,edmf_ent_dd,edmf_qc_dd,&
       &mynn_ql,                        &
       &ktop_shallow,    &
       &RTHRATEN,                       &
       &FLAG_QI,FLAG_QNI,FLAG_QC,FLAG_QNC &
       &,IDS,IDE,JDS,JDE,KDS,KDE        &
       &,IMS,IME,JMS,JME,KMS,KME        &
       &,ITS,ITE,JTS,JTE,KTS,KTE)
    
!-------------------------------------------------------------------

    INTEGER, INTENT(in) :: initflag
    !INPUT NAMELIST OPTIONS:
    INTEGER, INTENT(in) :: grav_settling
    INTEGER, INTENT(in) :: bl_mynn_tkebudget
    INTEGER, INTENT(in) :: bl_mynn_cloudpdf
    INTEGER, INTENT(in) :: bl_mynn_mixlength
    INTEGER, INTENT(in) :: bl_mynn_edmf
    INTEGER, INTENT(in) :: bl_mynn_edmf_dd
    REAL,    INTENT(in) :: bl_mynn_edmf_Lent
    INTEGER, INTENT(in) :: bl_mynn_edmf_mom
    INTEGER, INTENT(in) :: bl_mynn_edmf_tke
    INTEGER, INTENT(in) :: bl_mynn_edmf_part
    INTEGER, INTENT(in) :: bl_mynn_cloudmix
    INTEGER, INTENT(in) :: bl_mynn_mixqt
    INTEGER, INTENT(in) :: icloud_bl

    LOGICAL, INTENT(IN) :: FLAG_QI,FLAG_QNI,FLAG_QC,FLAG_QNC
    
    INTEGER,INTENT(IN) :: &
         & IDS,IDE,JDS,JDE,KDS,KDE &
         &,IMS,IME,JMS,JME,KMS,KME &
         &,ITS,ITE,JTS,JTE,KTS,KTE



! initflag > 0  for TRUE
! else        for FALSE
!       levflag         : <>3;  Level 2.5
!                         = 3;  Level 3
! grav_settling = 1 when gravitational settling accounted for
! grav_settling = 0 when gravitational settling NOT accounted for
    
    REAL, INTENT(in) :: delt,dx
    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(in) :: dz,&
         &u,v,w,th,qv,cc,p,exner,rho,T3D
    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), OPTIONAL, INTENT(in)::&
         &ql,qi,qni,qnc
    REAL, DIMENSION(IMS:IME,JMS:JME), INTENT(in) :: xland,ust,&
         &ch,rmol,ts,qsfc,qcg,ps,hfx,qfx, wspd,uoce,voce, vdfg,znt

    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(inout) :: &
         &Qke,Tsq,Qsq,Cov
       
          ! terms to couple EDMF with Tiedtke   
          ! Yi-Hsuan ... add INTENT(out) and output these terms 
        !REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: Q_ql,Q_qi,Q_a
        REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(out) :: Q_ql,Q_qi,Q_a, &
           Q_ql_det,Q_qi_det,Q_a_det, Q_ql_sub,Q_qi_sub,Q_a_sub, &
           Q_ql_adv,Q_qi_adv,Q_a_adv, Q_ql_eddy,Q_qi_eddy,Q_a_eddy, Q_ql_ent,Q_qi_ent,Q_a_ent

    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(inout) :: &
         &RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,RQLBLTEN,&
         &RQIBLTEN,RQNIBLTEN,RTHRATEN !,RQNCBLTEN

    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(out) :: &
         &exch_h,exch_m

   REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), OPTIONAL, INTENT(inout) :: &
         & edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_det,edmf_qc, &
         & edmf_a_dd,edmf_w_dd,edmf_qt_dd,edmf_thl_dd,edmf_ent_dd,edmf_qc_dd,& 
         & mynn_ql,edmf_debug1,edmf_debug2,edmf_debug3,edmf_debug4

    REAL, DIMENSION(IMS:IME,JMS:JME), INTENT(inout) :: Pblh

    REAL, DIMENSION(IMS:IME,JMS:JME) :: &
         &Psig_bl,Psig_shcu

    INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: & 
         &KPBL,ktop_shallow


    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(inout) :: &
         &el_pbl

    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(out) :: &
         &qWT,qSHEAR,qBUOY,qDISS,dqke
    ! 3D budget arrays are not allocated when bl_mynn_tkebudget == 0.
    ! 1D (local) budget arrays are used for passing between subroutines.
    REAL, DIMENSION(KTS:KTE) :: qWT1,qSHEAR1,qBUOY1,qDISS1,dqke1

    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: K_q,Sh3D

    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(inout) :: &
         &qc_bl,cldfra_bl
    REAL, DIMENSION(KTS:KTE) :: qc_bm,cldfra_bm,&
                             qc_am,cldfra_am
    REAl, DIMENSION(KTS:KTE) :: liquid_frac                        
  
    REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(out) :: RCCBLTEN,RTHLBLTEN,RQTBLTEN
  
!local vars
!  INTEGER :: ITF,JTF,KTF, IMD,JMD
    INTEGER :: i,j,k
    REAL, DIMENSION(KTS:KTE) :: thl,thvl,tl,sqv,sqc,sqi,sqw,sql,&
         &El, Dfm, Dfh, Dfq, Tcd, Qcd, Pdk, Pdt, Pdq, Pdc, &
         &Vt, Vq, sgm

    REAL, DIMENSION(KTS:KTE) :: thetav,sh,u1,v1,w1,p1,ex1,dz1,th1,tk1,rho1,&
           & qke1,tsq1,qsq1,cov1,qv1,qi1,ql1,qc1,du1,dv1,dth1,dqv1,dsqw1,dqi1,dql1,dthl1, &
           & k_m1,k_h1,k_q1,qni1,dqni1,qnc1,cc1,dcc1 !,dqnc1

!JOE: mass-flux variables
    REAL, DIMENSION(KTS:KTE) :: dth1mf,dqv1mf,dqc1mf,du1mf,dv1mf
    REAL, DIMENSION(KTS:KTE) :: edmf_a1,edmf_w1,edmf_qt1,edmf_thl1,&
                                edmf_ent1,edmf_qc1,&
                                edmf_a_dd1,edmf_w_dd1,edmf_qt_dd1,edmf_thl_dd1,&
                                edmf_ent_dd1,edmf_qc_dd1,&
                                edmf_det1, & ! yhc
                                edmf_debug11,edmf_debug21,edmf_debug31,edmf_debug41
    REAL,DIMENSION(KTS:KTE+1) :: s_aw1,s_awthl1,s_awqt1,&
                  s_awqv1,s_awqc1,s_awu1,s_awv1,s_awqke1,&
                  sd_aw1,sd_awthl1,sd_awqt1,&
                  sd_awqv1,sd_awqc1,sd_awu1,sd_awv1,sd_awqke1

    REAL, DIMENSION(KTS:KTE+1) :: zw,rhoh1
    REAL :: cpm,sqcg,flt,flq,flqv,flqc,pmz,phh,exnerg,zet,& 
              &afk,abk,ts_decay,th_sfc,ztop_shallow
    REAL :: dqcTT,lfTT,lvT


    REAL,DIMENSION(KTS:KTE) :: Q_ql1,Q_qi1,Q_a1, &
      Q_ql1_adv,Q_qi1_adv,Q_a1_adv, Q_ql1_eddy,Q_qi1_eddy,Q_a1_eddy, Q_ql1_ent,Q_qi1_ent,Q_a1_ent, &
      Q_ql1_sub,Q_qi1_sub,Q_a1_sub, Q_ql1_det,Q_qi1_det,Q_a1_det

    !<--- yhc 2021-09-08 

    !--- input argument
    type(randomNumberStream), dimension(IMS:IME,JMS:JME), INTENT(INOUT) :: &     ! yhc 2021-11-18
      streams

    !--- output argument
    REAL, DIMENSION(IMS:IME,KMS:KME+1,JMS:JME), INTENT(out) :: &
      num_updraft  , &
      a_moist_half , &
      a_dry_half   , &
      mf_moist_half, &
      mf_dry_half  , &
      mf_all_half  , &
      qv_moist_half, &
      qv_dry_half

    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(out) :: &
      num_DET, num_nDET_pENT, num_nDET_zENT,  &
      a_moist_full , &
      a_dry_full , &
      mf_moist_full,  &
      mf_dry_full  ,  &
      mf_all_full  , &
      qv_moist_full,  &
      qv_dry_full

    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(out) :: &
    !REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: &
       qa_before_mix, ql_before_mix, qi_before_mix, thl_before_mix, qt_before_mix, th_before_mix, &
       qa_before_pdf, ql_before_pdf, qi_before_pdf, &
       qa_after_mix, ql_after_mix, qi_after_mix, thl_after_mix, qt_after_mix, th_after_mix

    !REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: &
    !   qa_before_pdf, ql_before_pdf, qi_before_pdf

    !--- local argument
    REAL, DIMENSION(KMS:KME) :: &
       dum_1D

    REAL,DIMENSION(KTS:KTE+1) :: &
      num_updraft1  , &
      a_moist_half1 , &
      a_dry_half1   , &
      mf_moist_half1, &
      mf_dry_half1  , &
      mf_all_half1  , &
      qv_moist_half1, &
      qv_dry_half1

    REAL,DIMENSION(KTS:KTE) :: &
      num_DET1, num_nDET_pENT1, num_nDET_zENT1,  &
      a_moist_full1 ,  &
      a_dry_full1 ,  &
      mf_moist_full1,  &
      mf_dry_full1  ,  &
      mf_all_full1  , &
      qv_moist_full1,  &
      qv_dry_full1

    REAL, DIMENSION(IMS:IME,JMS:JME,KMS:KME) :: &
      diag_full

    !---> yhc 2021-09-08

! 0 ... default thing 
! 1 ... new solver


!JOE-add GRIMS parameters & variables
   real,parameter    ::  d1 = 0.02, d2 = 0.05, d3 = 0.001
   real,parameter    ::  h1 = 0.33333335, h2 = 0.6666667
   REAL :: govrth, sflux, bfx0, wstar3, wm2, wm3, delb
!JOE-end GRIMS
!JOE-top-down diffusion
   REAL, DIMENSION(ITS:ITE,JTS:JTE) :: maxKHtopdown
   REAL,DIMENSION(KTS:KTE) :: KHtopdown,zfac,wscalek2,&
                             zfacent,TKEprodTD
   REAL :: bfxpbl,dthvx,tmp1,temps,templ,zl1,wstar3_2
   real :: ent_eff,radsum,radflux,we,rcldb,rvls,&
           minrad,zminrad
   real, parameter :: pfac =2.0, zfmin = 0.01, phifac=8.0
   integer :: kk,kminrad
   logical :: cloudflg
!JOE-end top down

    INTEGER, SAVE :: levflag

! Stochastic fields 
     REAL, DIMENSION(KTS:KTE)                         ::    rstoch_col

    IF ( debug_code ) THEN
       print*,'in MYNN driver; at beginning'
    ENDIF
  
   levflag=mynn_level

!
! initialize updraft and downdraft variables
!
    IF (bl_mynn_edmf > 0) THEN
      edmf_a(its:ite,kts:kte,jts:jte)=0.
      edmf_w(its:ite,kts:kte,jts:jte)=0.
      edmf_qt(its:ite,kts:kte,jts:jte)=0.
      edmf_thl(its:ite,kts:kte,jts:jte)=0.
      edmf_ent(its:ite,kts:kte,jts:jte)=0.
      edmf_qc(its:ite,kts:kte,jts:jte)=0.
      ktop_shallow(its:ite,jts:jte)=0 !int
      edmf_debug1(its:ite,kts:kte,jts:jte)=0.
      edmf_debug2(its:ite,kts:kte,jts:jte)=0.
      edmf_debug3(its:ite,kts:kte,jts:jte)=0.
      edmf_debug4(its:ite,kts:kte,jts:jte)=0.
    ENDIF

    IF (bl_mynn_edmf_dd > 0) THEN
      edmf_a_dd(its:ite,kts:kte,jts:jte)=0.
      edmf_w_dd(its:ite,kts:kte,jts:jte)=0.
      edmf_qt_dd(its:ite,kts:kte,jts:jte)=0.
      edmf_thl_dd(its:ite,kts:kte,jts:jte)=0.
      edmf_ent_dd(its:ite,kts:kte,jts:jte)=0.
      edmf_qc_dd(its:ite,kts:kte,jts:jte)=0.
    ENDIF
    
       
    maxKHtopdown(its:ite,jts:jte)=0.


!
! initialization of turbulence fields when EDMF is run for the first time. Computes:
!  (a) PBLH & KPBL
!  (b) turbulent properties, using MYNN-level2: EL,SH,QKE,TSQ,QSQ,COV
!   


    IF (initflag > 0) THEN
 
       Sh3D(its:ite,kts:kte,jts:jte)=0.
       el_pbl(its:ite,kts:kte,jts:jte)=0.
       tsq(its:ite,kts:kte,jts:jte)=0.
       qsq(its:ite,kts:kte,jts:jte)=0.
       cov(its:ite,kts:kte,jts:jte)=0.
       dqi1(kts:kte)=0.0
       dqni1(kts:kte)=0.0
       qc_bm(kts:kte)=0.0
       cldfra_bm(kts:kte)=0.0
    !   qc_bl1D_old(kts:kte)=0.0
    !   cldfra_bl1D_old(kts:kte)=0.0
       sgm(kts:kte)=0.0
       vt(kts:kte)=0.0
       vq(kts:kte)=0.0

       DO j=JTS,JTE
          DO i=ITS,ITE
             DO k=KTS,KTE !KTF
                dz1(k)=dz(i,k,j)
                u1(k) = u(i,k,j)
                v1(k) = v(i,k,j)
                w1(k) = w(i,k,j)
                cc1(k)=cc(i,k,j)
                qc1(k)=ql(i,k,j)+qi(i,k,j)
                th1(k)=th(i,k,j)
                tk1(k)=T3D(i,k,j)
                rho1(k)=rho(i,k,j)
                sql(k)=ql(i,k,j)/(1.+ql(i,k,j))
                sqv(k)=qv(i,k,j)/(1.+qv(i,k,j))
                thetav(k)=th(i,k,j)*(1.+0.61*sqv(k))
                sqi(k)=qi(i,k,j)/(1.+qi(i,k,j))
                sqc(k)=sqi(k)+sql(k)
                sqw(k)=sqv(k)+sqc(k)
                thl(k)=th(i,k,j)- xlvcp/exner(i,k,j)*sql(k) &
                       &           - xlscp/exner(i,k,j)*sqi(k)
                IF (k==kts) THEN
                   zw(k)=0.
                ELSE
                   zw(k)=zw(k-1)+dz(i,k-1,j)
                ENDIF
                thvl(k)=thl(k)*(1.+0.61*sqv(k))
                exch_m(i,k,j)=0.
                exch_h(i,k,j)=0.
                K_q(i,k,j)=0.
                qke(i,k,j)=0.1-MIN(zw(k)*0.001, 0.0) !for initial PBLH calc only
                qke1(k)=qke(i,k,j)
                el(k)=el_pbl(i,k,j)
                sh(k)=Sh3D(i,k,j)
                tsq1(k)=tsq(i,k,j)
                qsq1(k)=qsq(i,k,j)
                cov1(k)=cov(i,k,j)
                rstoch_col(k)=0.0
                IF ( bl_mynn_tkebudget == 1) THEN
                   !TKE BUDGET VARIABLES
                   qWT(i,k,j)=0.
                   qSHEAR(i,k,j)=0.
                   qBUOY(i,k,j)=0.
                   qDISS(i,k,j)=0.
                   dqke(i,k,j)=0.
                ENDIF
             ENDDO

             zw(kte+1)=zw(kte)+dz(i,kte,j)

!         output: PBLH (boundary layer height) and KPBL(vertical level corresponding to PBLH) 
!                 QKE is taken from the prescribed value above
             CALL GET_PBLH(KTS,KTE,PBLH(i,j),thvl,&
               &  Qke1,zw,dz1,xland(i,j),KPBL(i,j))
             
!             IF (scaleaware > 0.) THEN
!                CALL SCALE_AWARE(dx,PBLH(i,j),Psig_bl(i,j),Psig_shcu(i,j))
!             ELSE
                Psig_bl(i,j)=1.0
                Psig_shcu(i,j)=1.0
!             ENDIF


             CALL mym_initialize (             & 
                  &kts,kte,                    &
                  &dz1, zw, u1, v1, thl, sqv,  &
                  &PBLH(i,j), th1, sh,         &
                  &ust(i,j), rmol(i,j),        &
                  &el, Qke1, Tsq1, Qsq1, Cov1, &
                  &Psig_bl(i,j), cldfra_bm,  &
                  &bl_mynn_mixlength,          &
                  &edmf_w1,edmf_a1,edmf_qc1,bl_mynn_edmf,&
                  &edmf_w_dd1,edmf_a_dd1,edmf_qc_dd1,bl_mynn_edmf_dd,&
                  &rstoch_col,KPBL(i,j) )


             !UPDATE 3D VARIABLE
             DO k=KTS,KTE !KTF
                el_pbl(i,k,j)=el(k)
                sh3d(i,k,j)=sh(k)
                qke(i,k,j)=qke1(k)
                tsq(i,k,j)=tsq1(k)
                qsq(i,k,j)=qsq1(k)
                cov(i,k,j)=cov1(k)
             ENDDO

          ENDDO
       ENDDO

    ENDIF 

!
! end of initialization
!




    DO j=JTS,JTE
       DO i=ITS,ITE
          DO k=KTS,KTE !KTF
            !JOE-TKE BUDGET
             IF ( bl_mynn_tkebudget == 1) THEN
                dqke(i,k,j)=qke(i,k,j)
             END IF
             dz1(k)= dz(i,k,j)
             u1(k) = u(i,k,j)
             v1(k) = v(i,k,j)
             w1(k) = w(i,k,j)
             th1(k)= th(i,k,j)
             tk1(k)=T3D(i,k,j)
             rho1(k)=rho(i,k,j)
             qv1(k)= qv(i,k,j)
             ql1(k)= ql(i,k,j)
             qi1(k)= qi(i,k,j)
             qc1(k)=ql1(k)+qi1(k)
             cc1(k)=cc(i,k,j)
             sqv(k)= qv(i,k,j)/(1.+qv(i,k,j))
             sql(k)= ql(i,k,j)/(1.+ql(i,k,j))
             sqi(k)= qi(i,k,j)/(1.+qi(i,k,j))
             sqc(k)=sqi(k)+sql(k)
             sqw(k)= sqv(k)+sqc(k)
             thl(k)= th(i,k,j) - xlvcp/exner(i,k,j)*sql(k) &
                     &            - xlscp/exner(i,k,j)*sqi(k)

            
  !           IF(icloud_bl > 0)cldfra_bl1D_old(k)=cldfra_bl(i,k,j)
  !           IF(icloud_bl > 0)qc_bl1D_old(k)=qc_bl(i,k,j)
          
             dql1(k)=0.0
             dqi1(k)=0.0
             dqni1(k)=0.0
      
    
             
          


             IF (PRESENT(qni) .AND. FLAG_QNI ) THEN
                qni1(k)=qni(i,k,j)
             ELSE
                qni1(k)=0.0
             ENDIF
             
             IF (PRESENT(qnc) .AND. FLAG_QNC ) THEN
                qnc1(k)=qnc(i,k,j)
             ELSE
                qnc1(k)=0.0
             ENDIF
             
             thetav(k)=th(i,k,j)*(1.+0.608*sqv(k))
             thvl(k)=thl(k)*(1.+0.61*sqv(k))
             p1(k) = p(i,k,j)
             ex1(k)= exner(i,k,j)
             el(k) = el_pbl(i,k,j)
             qke1(k)=qke(i,k,j)
             sh(k) = sh3d(i,k,j)
             tsq1(k)=tsq(i,k,j)
             qsq1(k)=qsq(i,k,j)
             cov1(k)=cov(i,k,j)
             
!if (spp_pbl==1) then
 !               rstoch_col(k)=pattern_spp_pbl(i,k,j)
 !            else
                rstoch_col(k)=0.0
 !            endif

             !edmf
             edmf_a1(k)=0.0
             edmf_w1(k)=0.0
             edmf_qc1(k)=0.0
             s_aw1(k)=0.
             s_awthl1(k)=0.
             s_awqt1(k)=0.
             s_awqv1(k)=0.
             s_awqc1(k)=0.
             s_awu1(k)=0.
             s_awv1(k)=0.
             s_awqke1(k)=0.
             ![EWDD]
             edmf_a_dd1(k)=0.0
             edmf_w_dd1(k)=0.0
             edmf_qc_dd1(k)=0.0
             sd_aw1(k)=0.
             sd_awthl1(k)=0.
             sd_awqt1(k)=0.
             sd_awqv1(k)=0.
             sd_awqc1(k)=0.
             sd_awu1(k)=0.
             sd_awv1(k)=0.
             sd_awqke1(k)=0.

             IF (k==kts) THEN
                zw(k)=0.
             ELSE
                zw(k)=zw(k-1)+dz(i,k-1,j)
             ENDIF
          ENDDO ! end k

! compute liquid fraction of possible condensate (1.= all liquid, 0. = all ice)
! which is a function of temperature and assumed to stay constant within the mixing time step
!
          call ComputeLiquidFrac(kts,kte,th1,ex1,liquid_frac)

          zw(kte+1)=zw(kte)+dz(i,kte,j)
          !EDMF
          s_aw1(kte+1)=0.
          s_awthl1(kte+1)=0.
          s_awqt1(kte+1)=0.
          s_awqv1(kte+1)=0.
          s_awqc1(kte+1)=0.
          s_awu1(kte+1)=0.
          s_awv1(kte+1)=0.
          s_awqke1(kte+1)=0.
          ! [EWDD]
          sd_aw1(kte+1)=0.
          sd_awthl1(kte+1)=0.
          sd_awqt1(kte+1)=0.
          sd_awqv1(kte+1)=0.
          sd_awqc1(kte+1)=0.
          sd_awu1(kte+1)=0.
          sd_awv1(kte+1)=0.
          sd_awqke1(kte+1)=0.

!         output: PBLH (boundary layer height) and KPBL(vertical level corresponding to PBLH) 
!                 QKE is taken from the previous time step
!         
          CALL GET_PBLH(KTS,KTE,PBLH(i,j),thvl,&
          & Qke1,zw,dz1,xland(i,j),KPBL(i,j))



        !  IF (scaleaware > 0.) THEN
        !     CALL SCALE_AWARE(dx,PBLH(i,j),Psig_bl(i,j),Psig_shcu(i,j))
        !  ELSE
             Psig_bl(i,j)=1.0
             Psig_shcu(i,j)=1.0
        !  ENDIF

          sqcg= 0.0   !JOE, it was: qcg(i,j)/(1.+qcg(i,j))
          cpm=cp*(1.+0.84*qv(i,kts,j))
          exnerg=(ps(i,j)/p1000mb)**rcp

          !-----------------------------------------------------
          !ORIGINAL CODE
          !flt = hfx(i,j)/( rho(i,kts,j)*cpm ) &
          ! +xlvcp*ch(i,j)*(sqc(kts)/exner(i,kts,j) -sqcg/exnerg)
          !flq = qfx(i,j)/  rho(i,kts,j)       &
          !    -ch(i,j)*(sqc(kts)   -sqcg )
          !-----------------------------------------------------
          ! Katata-added - The deposition velocity of cloud (fog)
          ! water is used instead of CH.
          flt = hfx(i,j)/( rho(i,kts,j)*cpm ) &
            & +xlvcp*vdfg(i,j)*(sqc(kts)/exner(i,kts,j)- sqcg/exnerg)
          flq = qfx(i,j)/  rho(i,kts,j)       &
            & -vdfg(i,j)*(sqc(kts) - sqcg )
!JOE-test- should this be after the call to mym_condensation?-using old vt & vq
!same as original form
!         flt = flt + xlvcp*ch(i,j)*(sqc(kts)/exner(i,kts,j) -sqcg/exnerg)
          flqv = qfx(i,j)/rho(i,kts,j)
          flqc = -vdfg(i,j)*(sqc(kts) - sqcg )
          th_sfc = ts(i,j)/ex1(kts)

          zet = 0.5*dz(i,kts,j)*rmol(i,j)
          if ( zet >= 0.0 ) then
            pmz = 1.0 + (cphm_st-1.0) * zet
            phh = 1.0 +  cphh_st      * zet
          else
            pmz = 1.0/    (1.0-cphm_unst*zet)**0.25 - zet
            phh = 1.0/SQRT(1.0-cphh_unst*zet)
          end if

!
! qc_bm and cldfra_bm are values of qc and cloud fraction before mixing
!
! condensation uses sh and el from the previous time step 
! 

    if (edmf_type .eq. 0) then
!     
! run the PDF scheme to compute qc and cloud fraction before mixing     
!     
          CALL  mym_condensation ( kts,kte,      &
               &dx,dz1,thl,sqw,p1,ex1,           &
               &tsq1, qsq1, cov1,                &
               &Sh,el,liquid_frac,          &
               &qc_bm,cldfra_bm,             &
               &PBLH(i,j),HFX(i,j),              &
               &Vt, Vq, th1, sgm )

          !<--- yhc, save cloud fraction from the PDF cloud scheme
          qa_before_pdf(i,:,j) = cldfra_bm(:)
          ql_before_pdf(i,:,j) = liquid_frac(:)      * qc_bm(:)   ! cloud liquid water content (kg/kg)
          qi_before_pdf(i,:,j) = (1.-liquid_frac(:)) * qc_bm(:)   ! cloud ice    water content (kg/kg)
          !---> yhc 

    !elseif (edmf_type .eq. 1) then
    elseif (edmf_type .eq. 1 .or. edmf_type .eq. 2) then  ! yhc add edmf_type=2, 2021-08-31

          !<--- yhc, save cloud fraction from the PDF cloud scheme
          CALL  mym_condensation ( kts,kte,      &
               &dx,dz1,thl,sqw,p1,ex1,           &
               &tsq1, qsq1, cov1,                &
               &Sh,el,liquid_frac,          &
               &qc_bm,cldfra_bm,             &
               &PBLH(i,j),HFX(i,j),              &
               &Vt, Vq, th1, sgm )

          qa_before_pdf(i,:,j) = cldfra_bm(:)
          ql_before_pdf(i,:,j) = liquid_frac(:)      * qc_bm(:)   ! cloud liquid water content (kg/kg)
          qi_before_pdf(i,:,j) = (1.-liquid_frac(:)) * qc_bm(:)   ! cloud ice    water content (kg/kg)
          !---> yhc
!   
! input values of qc and cc are taken as the values befor mixing
!   
           qc_bm=qc1
           cldfra_bm=cc1
!
! compute vt and vq for TKE buoyancy calculation   
!
         call  mym_vtvq (kts,kte,  &
          &th1,qc1,thl,sqw,cc1,    &
          &p1,ex1,liquid_frac,     &
          &Vt, Vq)
    else
          print *,'Error ... wrong edmf_type!!!'    
    endif
        
    !<--- yhc_mynn, output before-mixing states, 2021-04-02
    qa_before_mix(i,:,j) = cldfra_bm(:)                     ! cloud fraction (unit: none)
    ql_before_mix(i,:,j) = liquid_frac(:)      * qc_bm(:)   ! cloud liquid water content (kg/kg)
    qi_before_mix(i,:,j) = (1.-liquid_frac(:)) * qc_bm(:)   ! cloud ice    water content (kg/kg)
    th_before_mix(i,:,j) = th1(:)                           ! potential temperature (K)

    thl_before_mix (i,:,j)  = thl(:)                        ! ice-liquid potential temperature (K)
  
    !sqw =sqw +Dsqw1*delt                                   ! total water content (vapor+liquid+ice) (kg/kg)
    dum_1D =sqw                                             ! total water content (vapor+liquid+ice) (kg/kg)
    qt_before_mix   (i,:,j)  = dum_1D(:)/(1.-dum_1D(:))                       !

!print*,'qa_before_mix',qa_before_mix
!print*,'qa_before_pdf',qa_before_pdf

    !---> yhc_mynn
                           
 
!          !ADD TKE source driven by cloud top cooling
!          IF (bl_mynn_topdown.eq.1)then
!             cloudflg=.false.
!             minrad=100.
!             kminrad=kpbl(i,j)
!             zminrad=PBLH(i,j)
!             KHtopdown(kts:kte)=0.0
!             TKEprodTD(kts:kte)=0.0
!             maxKHtopdown(i,j)=0.0
!             !CHECK FOR STRATOCUMULUS-TOPPED BOUNDARY LAYERS
!             DO kk = MAX(1,kpbl(i,j)-2),kpbl(i,j)+3
!                if(sqc(kk).gt. 1.e-6 .OR. sqi(kk).gt. 1.e-6 .OR. &
!                   cldfra_bl1D(kk).gt.0.5) then
!                   cloudflg=.true.
!                endif
!                if(rthraten(i,kk,j) < minrad)then
!                   minrad=rthraten(i,kk,j)
!                   kminrad=kk
!                   zminrad=zw(kk) + 0.5*dz1(kk)
!                endif
!             ENDDO
!             IF (cloudflg) THEN
!                zl1 = dz1(kts)
!                k = MAX(kpbl(i,j)-1, kminrad-1)
!                !Best estimate of height of TKE source (top of downdrafts):
!                !zminrad = 0.5*pblh(i,j) + 0.5*zminrad
!
!                templ=thl(k)*ex1(k)
!                !rvls is ws at full level
!                rvls=100.*6.112*EXP(17.67*(templ-273.16)/(templ-29.65))*(ep_2/p1(k+1))
!                temps=templ + (sqw(k)-rvls)/(cp/xlv  +  ep_2*xlv*rvls/(rd*templ**2))
!                rvls=100.*6.112*EXP(17.67*(temps-273.15)/(temps-29.65))*(ep_2/p1(k+1))
!                rcldb=max(sqw(k)-rvls,0.)
!
!                !entrainment efficiency
!                dthvx     = (thl(k+2) + th1(k+2)*ep_1*sqw(k+2)) &
!                          - (thl(k)   + th1(k)  *ep_1*sqw(k))
!                dthvx     = max(dthvx,0.1)
!                tmp1      = xlvcp * rcldb/(ex1(k)*dthvx)
!                !Originally from Nichols and Turton (1986), where a2 = 60, but lowered
!                !here to 8, as in Grenier and Bretherton (2001).
!                ent_eff   = 0.2 + 0.2*8.*tmp1
!
!                radsum=0.
!                DO kk = MAX(1,kpbl(i,j)-3),kpbl(i,j)+3
!                   radflux=rthraten(i,kk,j)*ex1(kk)         !converts theta/s to temp/s
!                   radflux=radflux*cp/g*(p1(kk)-p1(kk+1)) ! converts temp/s to W/m^2
!                   if (radflux < 0.0 ) radsum=abs(radflux)+radsum
!                ENDDO
!                radsum=MIN(radsum,60.0)
!
!                !entrainment from PBL top thermals
!                bfx0 = max(radsum/rho1(k)/cp - max(sflux,0.0),0.)
!                !bfx0 = max(radsum/rho1(k)/cp,0.)
!                wm3    = g/thetav(k)*bfx0*MIN(pblh(i,j),1500.) ! this is wstar3(i)
!                wm2    = wm2 + wm3**h2
!                bfxpbl = - ent_eff * bfx0
!                dthvx  = max(thetav(k+1)-thetav(k),0.1)
!                we     = max(bfxpbl/dthvx,-sqrt(wm3**h2))
!
!                DO kk = kts,kpbl(i,j)+3
!                   !Analytic vertical profile
!                   zfac(kk) = min(max((1.-(zw(kk+1)-zl1)/(zminrad-zl1)),zfmin),1.)
!                   zfacent(kk) = 10.*MAX((zminrad-zw(kk+1))/zminrad,0.0)*(1.-zfac(kk))**3
!
!                   !Calculate an eddy diffusivity profile (not used at the moment)
!                   wscalek2(kk) = (phifac*karman*wm3*(zfac(kk)))**h1
!                   !Modify shape of KH to be similar to Lock et al (2000): use pfac = 3.0
!                   KHtopdown(kk) = wscalek2(kk)*karman*(zminrad-zw(kk+1))*(1.-zfac(kk))**3 !pfac
!                   KHtopdown(kk) = MAX(KHtopdown(kk),0.0)
!                   !Do not include xkzm at kpbl-1 since it changes entrainment
!                   !if (kk.eq.kpbl(i,j)-1 .and. cloudflg .and. we.lt.0.0) then
!                   !   KHtopdown(kk) = 0.0
!                   !endif
!                   
!                   !Calculate TKE production = 2(g/TH)(w'TH'), where w'TH' = A(TH/g)wstar^3/PBLH,
!                   !A = ent_eff, and wstar is associated with the radiative cooling at top of PBL.
!                   !An analytic profile controls the magnitude of this TKE prod in the vertical. 
!                   TKEprodTD(kk)=2.*ent_eff*wm3/MAX(pblh(i,j),100.)*zfacent(kk)
!                   TKEprodTD(kk)= MAX(TKEprodTD(kk),0.0)
!                ENDDO
!             ENDIF !end cloud check
!             maxKHtopdown(i,j)=MAXVAL(KHtopdown(:))
!          ELSE
             maxKHtopdown(i,j)=0.0
             KHtopdown(kts:kte) = 0.0
             TKEprodTD(kts:kte)=0.0
!          ENDIF 
!          !end top-down check


          IF (bl_mynn_edmf > 0) THEN
            CALL edmf_JPL(                            &
               &kts,kte,delt,zw,p1,                   &
               !&u1,v1,th1,thl,thetav,tk1,sqw,sqv,sqc, &  yhc comment out
               &u1,v1,th1,thl,thetav,tk1,sqw,sqv,sqc,sql,sqi, &  ! yhc add sql and sqi
               &ex1,rho1,                            &
               &ust(i,j),ps(i,j),flt,flq,PBLH(i,j),       &
               & liquid_frac,                     &
               & edmf_a1,edmf_w1,edmf_qt1,        &
               & edmf_thl1,edmf_ent1,edmf_qc1,    &
               & edmf_debug11,edmf_debug21,       &
               & edmf_debug31,edmf_debug41,       &
               & s_aw1,s_awthl1,s_awqt1,          &
               & s_awqv1,s_awqc1,s_awu1,s_awv1,   &
               & s_awqke1,                        &
               & qc_bm,cldfra_bm,             &
               & a_moist_half1, mf_moist_half1, qv_moist_half1, a_moist_full1, mf_moist_full1, qv_moist_full1, &  ! yhc 2021-09-08
               & a_dry_half1, mf_dry_half1, qv_dry_half1, a_dry_full1, mf_dry_full1, qv_dry_full1, &            ! yhc 2021-09-08
               & mf_all_half1, mf_all_full1, edmf_det1, &                                                     ! yhc 2021-09-08
               & num_updraft1, num_DET1, num_nDET_zENT1, num_nDET_pENT1, &                                               ! yhc 2021-09-08
               &ktop_shallow(i,j),ztop_shallow,   &
               & KPBL(i,j),                        &
               & Q_ql1,Q_qi1,Q_a1,                 &
               & streams(i,j),   & ! ych 2021-11-18
               & Q_ql1_adv,Q_qi1_adv,Q_a1_adv, Q_ql1_eddy,Q_qi1_eddy,Q_a1_eddy, Q_ql1_ent,Q_qi1_ent,Q_a1_ent, Q_ql1_det,Q_qi1_det,Q_a1_det, Q_ql1_sub,Q_qi1_sub,Q_a1_sub  &
             )

          ENDIF

!          IF (bl_mynn_edmf_dd == 1) THEN
!            ! DO k = kts, KtE
!            !     print *, "[EWDEBUG] k = ",k, " sqc = ", sqc(k), " cldfra = ", cldfra_bl1d(k)
!            ! ENDDO
!            CALL DDMF_JPL(kts,kte,delt,zw,p1,        &
!              &u1,v1,th1,thl,thetav,tk1,sqw,sqv,sqc, &
!              &ex1,                                  &
!              &ust(i,j),flt,flq,PBLH(i,j),KPBL(i,j), &
!              &edmf_a_dd1,edmf_w_dd1, edmf_qt_dd1,   &
!              &edmf_thl_dd1,edmf_ent_dd1,edmf_qc_dd1,&
!              & edmf_debug31,edmf_debug41,           &
!              &sd_aw1,sd_awthl1,sd_awqt1,            &
!              &sd_awqv1,sd_awqc1,sd_awu1,sd_awv1,    &
!              &sd_awqke1,                            &
!              &qc_bl1d,cldfra_bl1d,                  &
!              &rthraten(i,:,j)&
!              )
!          ENDIF

          CALL mym_turbulence (                  & 
               &kts,kte,levflag,                 &
               &dz1, zw, u1, v1, thl, sqc, sqw,  &
               &qke1, tsq1, qsq1, cov1,          &
               &vt, vq,                          &
               &rmol(i,j), flt, flq,             &
               &PBLH(i,j),th1,                   &
               &Sh,el,                           &
               &Dfm,Dfh,Dfq,                     &
               &Tcd,Qcd,Pdk,                     &
               &Pdt,Pdq,Pdc,                     &
               &qWT1,qSHEAR1,qBUOY1,qDISS1,      &
               &bl_mynn_tkebudget,               &
               &Psig_bl(i,j),Psig_shcu(i,j),     &     
               &cldfra_bm,bl_mynn_mixlength,   &
               &edmf_w1,edmf_a1,edmf_qc1,bl_mynn_edmf,   &
               &edmf_w_dd1,edmf_a_dd1,edmf_qc_dd1,bl_mynn_edmf_dd,   &
               &TKEprodTD,                       &
               &rstoch_col,KPBL(i,j),KHtopdown)




          IF (bl_mynn_edmf_part > 0 .AND. bl_mynn_edmf > 0) THEN
             !Partition the fluxes from each component (ed & mf).
             !Assume overlap of 50%: Reduce eddy diffusivities by 50% of the estimated
             !area fraction of mass-flux scheme's updraft.
             DO k=kts,kte
                dfm(k)=dfm(k) * (1. - 0.5*edmf_a1(k) - 0.5*edmf_a_dd1(k))
                dfh(k)=dfh(k) * (1. - 0.5*edmf_a1(k) - 0.5*edmf_a_dd1(k))
                dfq(k)=dfq(k) * (1. - 0.5*edmf_a1(k) - 0.5*edmf_a_dd1(k))
             ENDDO
          ENDIF
              
          CALL mym_predict (kts,kte,levflag,     &
               &delt, dz1,                       &
               &ust(i,j), flt, flq, pmz, phh,    &
               &el, dfq, pdk, pdt, pdq, pdc,     &
               &Qke1, Tsq1, Qsq1, Cov1,          &
               &s_aw1, s_awqke1, sd_aw1, sd_awqke1,&
               bl_mynn_edmf_tke)


          rhoh1(2:kte)=0.5*(rho1(1:kte-1)+rho1(2:kte))
          !rhoh1(kte+1)=rho1(kte)
          rhoh1(1)=ps(i,j)/(ts(i,j)*rdgas)
          rhoh1(1)=rhoh1(2)      

   
          CALL mynn_tendencies(kts,kte,          &
               &levflag,grav_settling,rho1,rhoh1,&
               &delt, dz1,                       &
               &u1, v1, th1, tk1, qv1, qc1, qi1, &
               &qni1,qnc1,                       &
               &p1, ex1, thl, sqv, sqc, sqi, sqw,&
               &ust(i,j),flt,flq,flqv,flqc,      &
               &wspd(i,j),qcg(i,j),              &
               &uoce(i,j),voce(i,j),             &
               &tsq1, qsq1, cov1,                &
               &tcd, qcd,                        &
               &dfm, dfh, dfq,                   &
               &Du1, Dv1, Dthl1, Dsqw1,            &
               &Dqni1,               & !Dqnc1,        &
               &vdfg(i,j),                       & !JOE/Katata- fog deposition
               ! mass flux components
               &s_aw1,s_awthl1,s_awqt1,          &
               &s_awu1,s_awv1,   &
               &sd_aw1,sd_awthl1,sd_awqt1,          &
               &sd_awu1,sd_awv1,   &
               &ztop_shallow,ktop_shallow(i,j),  &
               &bl_mynn_cloudmix,                &
               &bl_mynn_mixqt,                   &
               &bl_mynn_edmf,                    &
               &bl_mynn_edmf_dd,                 &
               &bl_mynn_edmf_mom,                &
               &liquid_frac,                     &
               &dx,PBLH(i,j),HFX(i,j),Sh,el,Vt,Vq)


! run cloud scheme to compute cloudiness after mixing

          CALL  mym_condensation ( kts,kte,      &
               &dx,dz1,thl,sqw,p1,ex1,           &
               &tsq1, qsq1, cov1,                &
               &Sh,el,liquid_frac,          &
               &qc_am,cldfra_am,             &
               &PBLH(i,j),HFX(i,j),              &
               &Vt, Vq, th1, sgm )

 
 !
 ! compute tendencies of dry variables
 !  d theta/dt=d thetal/dt + L_blend/(cp*pi) dql/dt
 !  d qv/dt=d qt/dt - d qc /dt
 !  d qi/dt =(1.-liquid_frac) dqc/dt; d ql/dt =liquid_frac dqc/dt
 
 
         DO k=kts,kte
           dqcTT=(qc_am(k)-qc_bm(k))/delt
           lfTT=liquid_frac(k)
           lvT=xlLF_blend(th1(k)*ex1(k),liquid_frac(k))
           
           dqv1(k)=dsqw1(k)/(1.-sqw(k))-dqcTT
           dql1(k)=dqcTT*lfTT
           dqi1(k)=(1.-lfTT)*dqcTT
           dth1(k)=dthl1(k)+lvT/(ex1(k)*cp)*dqcTT/(1.+qc1(k))
           dcc1(k)=(cldfra_am(k)-cldfra_bm(k))/delt 
         ENDDO
         
    !<--- yhc_mynn, output after-mixing states, 2021-04-02
    qa_after_mix  (i,:,j)  = cldfra_am(:)                     ! cloud fraction (unit: none)
    ql_after_mix  (i,:,j)  = liquid_frac(:)      * qc_am(:)   ! cloud liquid water content (kg/kg)
    qi_after_mix  (i,:,j)  = (1.-liquid_frac(:)) * qc_am(:)   ! cloud ice    water content (kg/kg)
    th_after_mix  (i,:,j)  = th1(:)+dth1*delt                 ! potential temperature (K)
    thl_after_mix (i,:,j)  = thl(:)                           ! ice-liquid potential temperature (K)
    qt_after_mix  (i,:,j)  = sqw(:)/(1.-sqw(:))               ! total water content (vapor+liquid+ice) (kg/kg)
    !---> yhc_mynn
 
          CALL retrieve_exchange_coeffs(kts,kte,&
               &dfm, dfh, dfq, dz1,&
               &K_m1, K_h1, K_q1)

          !UPDATE 3D ARRAYS
          DO k=KTS,KTE !KTF
             exch_m(i,k,j)=K_m1(k)
             exch_h(i,k,j)=K_h1(k)
             K_q(i,k,j)=K_q1(k)
             RUBLTEN(i,k,j)=du1(k)
             RVBLTEN(i,k,j)=dv1(k)
             RTHBLTEN(i,k,j)=dth1(k)
             RQVBLTEN(i,k,j)=dqv1(k)
             RQLBLTEN(i,k,j)=dql1(k)
             RQIBLTEN(i,k,j)=dqi1(k)
             RCCBLTEN(i,k,j)=dcc1(k)
             RTHLBLTEN(i,k,j)=dthl1(k)
             RQTBLTEN(i,k,j)=dsqw1(k)/(1.-sqw(k))
             IF (PRESENT(qni) .AND. FLAG_QNI) RQNIBLTEN(i,k,j)=dqni1(k)
             qc_bl(i,k,j)=qc_am(k)
             cldfra_bl(i,k,j)=cldfra_am(k) 
             
             
              
  !           IF(bl_mynn_cloudmix > 0)THEN
  !             IF (PRESENT(qc) .AND. FLAG_QC) RQCBLTEN(i,k,j)=dqc1(k)
  !             IF (PRESENT(qi) .AND. FLAG_QI) RQIBLTEN(i,k,j)=dqi1(k)
  !             !IF (PRESENT(qnc) .AND. FLAG_QNC) RQNCBLTEN(i,k,j)=dqnc1(k)
  !             IF (PRESENT(qni) .AND. FLAG_QNI) RQNIBLTEN(i,k,j)=dqni1(k)
  !           ELSE
  !             IF (PRESENT(qc) .AND. FLAG_QC) RQCBLTEN(i,k,j)=0.
  !             IF (PRESENT(qi) .AND. FLAG_QI) RQIBLTEN(i,k,j)=0.
  !             !IF (PRESENT(qnc) .AND. FLAG_QNC) RQNCBLTEN(i,k,j)=0.
  !             IF (PRESENT(qni) .AND. FLAG_QNI) RQNIBLTEN(i,k,j)=0.
  !           ENDIF

  !           IF(icloud_bl > 0)THEN
  !             !make BL clouds scale aware - may already be done in mym_condensation
  !             qc_bl(i,k,j)=qc_bl1D(k) !*Psig_shcu(i,j)
  !             cldfra_bl(i,k,j)=cldfra_bl1D(k) !*Psig_shcu(i,j)

               !Stochastic perturbations to cldfra_bl and qc_bl
              ! if (spp_pbl==1) then
              !    cldfra_bl(i,k,j)= cldfra_bl(i,k,j)*(1.0-rstoch_col(k))
              !    IF ((cldfra_bl(i,k,j) > 1.0) .or. (cldfra_bl(i,k,j) < 0.0)) then
              !       cldfra_bl(i,k,j)=MAX(MIN(cldfra_bl(i,k,j),1.0),0.0)
              !    ENDIF
              ! ELSE
                  !DIAGNOSTIC-DECAY FOR SUBGRID-SCALE CLOUDS
  !                IF (CLDFRA_BL(i,k,j) > cldfra_bl1D_old(k)) THEN
  !                   !KEEP UPDATED CLOUD FRACTION & MIXING RATIO
  !                ELSE
  !                   !DECAY TIMESCALE FOR CALM CONDITION IS THE EDDY TURNOVER TIMESCALE, 
  !                   !BUT FOR WINDY CONDITIONS, IT IS THE ADVECTIVE TIMESCALE.
  !                   !USE THE MINIMUM OF THE TWO.
  !                   ts_decay = MIN( 1800., 3.*dx/MAX(SQRT(u1(k)**2 + v1(k)**2), 1.0) )
  !                   cldfra_bl(i,k,j)= MAX(cldfra_bl1D(k), cldfra_bl1D_old(k)-(0.25*delt/ts_decay))
  !                   IF (cldfra_bl(i,k,j) > 0.01) THEN
  !                      IF (QC_BL(i,k,j) < 1E-5)QC_BL(i,k,j)= MAX(qc_bl1D_old(k), 1E-5)
  !                   ELSE
  !                      CLDFRA_BL(i,k,j)= 0.
  !                      QC_BL(i,k,j)    = 0.
  !                   ENDIF
  !                ENDIF
               !ENDIF

               !Reapply checks on cldfra_bl and qc_bl to avoid FPEs in radiation driver
               ! when these two quantities are multiplied by eachother (they may have changed
               ! in the MF scheme:
 !              IF (icloud_bl > 0) THEN
 !                IF (QC_BL(i,k,j) < 1E-6 .AND. ABS(CLDFRA_BL(i,k,j)) > 0.1)QC_BL(i,k,j)= 1E-6
 !                IF (CLDFRA_BL(i,k,j) < 1E-2)CLDFRA_BL(i,k,j)= 0.
 !              ENDIF
 !            ENDIF
             mynn_ql(i,k,j) =QC_BL(i,k,j)
             el_pbl(i,k,j)=el(k)
             qke(i,k,j)=qke1(k)
             tsq(i,k,j)=tsq1(k)
             qsq(i,k,j)=qsq1(k)
             cov(i,k,j)=cov1(k)
             sh3d(i,k,j)=sh(k)

             IF ( bl_mynn_tkebudget == 1) THEN
                dqke(i,k,j)  = (qke1(k)-dqke(i,k,j))*0.5  !qke->tke
                qWT(i,k,j)   = qWT1(k)*delt
                qSHEAR(i,k,j)= qSHEAR1(k)*delt
                qBUOY(i,k,j) = qBUOY1(k)*delt
                qDISS(i,k,j) = qDISS1(k)*delt
             ENDIF

             !update updraft properties
             IF (bl_mynn_edmf > 0) THEN
               edmf_a(i,k,j)=edmf_a1(k)
               edmf_w(i,k,j)=edmf_w1(k)
               edmf_qt(i,k,j)=edmf_qt1(k)
               edmf_thl(i,k,j)=edmf_thl1(k)
               edmf_ent(i,k,j)=edmf_ent1(k)
               edmf_qc(i,k,j)=edmf_qc1(k)
               edmf_debug1(i,k,j)=edmf_debug11(k)
               edmf_debug2(i,k,j)=edmf_debug21(k)
               edmf_debug3(i,k,j)=edmf_debug31(k)
               edmf_debug4(i,k,j)=edmf_debug41(k)
               
               
               Q_ql(i,k,j)=Q_ql1(k)
               Q_qi(i,k,j)=Q_qi1(k)
               Q_a(i,k,j)=Q_a1(k)

               Q_ql_adv(i,k,j)=Q_ql1_adv(k)
               Q_qi_adv(i,k,j)=Q_qi1_adv(k)
               Q_a_adv(i,k,j)=Q_a1_adv(k)

               Q_ql_eddy(i,k,j)=Q_ql1_eddy(k)
               Q_qi_eddy(i,k,j)=Q_qi1_eddy(k)
               Q_a_eddy(i,k,j)=Q_a1_eddy(k)

               Q_ql_ent(i,k,j)=Q_ql1_ent(k)
               Q_qi_ent(i,k,j)=Q_qi1_ent(k)
               Q_a_ent(i,k,j)=Q_a1_ent(k)

               Q_ql_det(i,k,j)=Q_ql1_det(k)
               Q_qi_det(i,k,j)=Q_qi1_det(k)
               Q_a_det(i,k,j)=Q_a1_det(k)

               Q_ql_sub(i,k,j)=Q_ql1_sub(k)
               Q_qi_sub(i,k,j)=Q_qi1_sub(k)
               Q_a_sub(i,k,j)=Q_a1_sub(k)

               !<--- yhc 2021-09-08
               a_moist_half  (i,k,j) = a_moist_half1(k)    
               a_moist_full  (i,k,j) = a_moist_full1(k)    
               mf_moist_half (i,k,j) = mf_moist_half1(k)    
               mf_moist_full (i,k,j) = mf_moist_full1(k)
               qv_moist_half (i,k,j) = qv_moist_half1(k)    
               qv_moist_full (i,k,j) = qv_moist_full1(k)
             
               a_dry_half    (i,k,j) = a_dry_half1(k)    
               a_dry_full    (i,k,j) = a_dry_full1(k)    
               mf_dry_half   (i,k,j) = mf_dry_half1(k)
               mf_dry_full   (i,k,j) = mf_dry_full1(k)
               qv_dry_half   (i,k,j) = qv_dry_half1(k)
               qv_dry_full   (i,k,j) = qv_dry_full1(k)
             
               mf_all_half   (i,k,j) = mf_all_half1(k)
               mf_all_full   (i,k,j) = mf_all_full1(k)

               edmf_det(i,k,j)=edmf_det1(k)

               num_updraft  (i,k,j) = num_updraft1  (k)
               num_DET      (i,k,j) = num_DET1      (k)
               num_nDET_zENT(i,k,j) = num_nDET_zENT1(k)
               num_nDET_pENT(i,k,j) = num_nDET_pENT1(k)

               !---> yhc 2021-09-08

               ELSE
               
                Q_ql(i,k,j)=0.
                Q_qi(i,k,j)=0.
                Q_a(i,k,j)=0.

                Q_ql_adv(i,k,j)=0.
                Q_qi_adv(i,k,j)=0.
                Q_a_adv(i,k,j)=0.
 
                Q_ql_eddy(i,k,j)=0.
                Q_qi_eddy(i,k,j)=0.
                Q_a_eddy(i,k,j)=0.

                Q_ql_ent(i,k,j)=0.
                Q_qi_ent(i,k,j)=0.
                Q_a_ent(i,k,j)=0.

                Q_ql_det(i,k,j)=0.
                Q_qi_det(i,k,j)=0.
                Q_a_det(i,k,j)=0.

                Q_ql_sub(i,k,j)=0.
                Q_qi_sub(i,k,j)=0.
                Q_a_sub(i,k,j)=0.
               
                ENDIF
               
               
               
               if (bl_mynn_edmf_dd > 0) THEN
                   !update downdraft properties
                   edmf_a_dd(i,k,j)=edmf_a_dd1(k)
                   edmf_w_dd(i,k,j)=edmf_w_dd1(k)
                   edmf_qt_dd(i,k,j)=edmf_qt_dd1(k)
                   edmf_thl_dd(i,k,j)=edmf_thl_dd1(k)
                   edmf_ent_dd(i,k,j)=edmf_ent_dd1(k)
                   edmf_qc_dd(i,k,j)=edmf_qc_dd1(k)
               ENDIF

            

             !***  Begin debug prints
             IF ( debug_code ) THEN
               IF ( sh(k) < 0. .OR. sh(k)> 200.)print*,&
                  "SUSPICIOUS VALUES AT: i,j,k=",i,j,k," sh=",sh(k)
               IF ( qke(i,k,j) < -1. .OR. qke(i,k,j)> 200.)print*,&
                  "SUSPICIOUS VALUES AT: i,j,k=",i,j,k," qke=",qke(i,k,j)
               IF ( el_pbl(i,k,j) < 0. .OR. el_pbl(i,k,j)> 2000.)print*,&
                  "SUSPICIOUS VALUES AT: i,j,k=",i,j,k," el_pbl=",el_pbl(i,k,j)
               IF ( ABS(vt(k)) > 0.8 )print*,&
                  "SUSPICIOUS VALUES AT: i,j,k=",i,j,k," vt=",vt(k)
               IF ( ABS(vq(k)) > 6000.)print*,&
                  "SUSPICIOUS VALUES AT: i,j,k=",i,j,k," vq=",vq(k) 
               IF ( exch_m(i,k,j) < 0. .OR. exch_m(i,k,j)> 2000.)print*,&
                  "SUSPICIOUS VALUES AT: i,j,k=",i,j,k," exxch_m=",exch_m(i,k,j)
               IF ( vdfg(i,j) < 0. .OR. vdfg(i,j)>5. )print*,&
                  "SUSPICIOUS VALUES AT: i,j,k=",i,j,k," vdfg=",vdfg(i,j)
               IF ( ABS(QFX(i,j))>.001)print*,&
                  "SUSPICIOUS VALUES AT: i,j=",i,j," QFX=",QFX(i,j)
               IF ( ABS(HFX(i,j))>1000.)print*,&
                  "SUSPICIOUS VALUES AT: i,j=",i,j," HFX=",HFX(i,j)
               IF (icloud_bl > 0) then
                  IF( cldfra_bl(i,k,j) < 0.0 .OR. cldfra_bl(i,k,j)> 1.)THEN
                  PRINT*,"SUSPICIOUS VALUES: CLDFRA_BL=",cldfra_bl(i,k,j)," qc_bl=",QC_BL(i,k,j)
                  ENDIF
               ENDIF
             ENDIF
             !***  End debug prints
             
             
          ENDDO

       ENDDO
    ENDDO




     !print *,'Q_ql',Q_ql
     !print *,'Q_qi',Q_qi
     !print *,'Q_a',Q_a


!print *,'qkeEND',qke



  END SUBROUTINE mynn_bl_driver

subroutine  ComputeLiquidFrac(kts,kte,th,ex,liquid_frac)
     integer, intent(in) :: kts,kte
     real,dimension(kts:kte), intent(in) :: th,ex
     real,dimension(kts:kte), intent(out) :: liquid_frac

     integer k
     real t    
 
     
     DO k=kts,kte
       t=ex(k)*th(k)
     
      IF (t .GE. 273.16) THEN
            liquid_frac(k)=1.
      ELSE IF (t .LE. 253.) THEN
           liquid_frac(k)=0.
      ELSE
          liquid_frac(k)  = 1.-(273.16-t)/20.16
      END IF
     
    ENDDO

end subroutine ComputeLiquidFrac



! ==================================================================
  SUBROUTINE mynn_bl_init_driver(                   &
       &RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,          &
       &RQLBLTEN,RQIBLTEN & !,RQNIBLTEN,RQNCBLTEN       &
       &,QKE,EXCH_H                         &
!       &,icloud_bl,qc_bl,cldfra_bl                 & !JOE-subgrid bl clouds 
       &,RESTART,ALLOWED_TO_READ,LEVEL              &
       &,IDS,IDE,JDS,JDE,KDS,KDE                    &
       &,IMS,IME,JMS,JME,KMS,KME                    &
       &,ITS,ITE,JTS,JTE,KTS,KTE)

    !---------------------------------------------------------------
    LOGICAL,INTENT(IN) :: ALLOWED_TO_READ,RESTART
    INTEGER,INTENT(IN) :: LEVEL !,icloud_bl

    INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE,                    &
         &                IMS,IME,JMS,JME,KMS,KME,                    &
         &                ITS,ITE,JTS,JTE,KTS,KTE
    
    
    REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(INOUT) :: &
         &RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,                 &
         &RQLBLTEN,RQIBLTEN,& !RQNIBLTEN,RQNCBLTEN       &
         &QKE,EXCH_H

!    REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(INOUT) :: &
!         &qc_bl,cldfra_bl

    INTEGER :: I,J,K
    
 
    
    IF(.NOT.RESTART)THEN
       DO J=JTS,JTE
          DO K=KTS,KTE
             DO I=ITS,ITE
                RUBLTEN(i,k,j)=0.
                RVBLTEN(i,k,j)=0.
                RTHBLTEN(i,k,j)=0.
                RQVBLTEN(i,k,j)=0.
                if( p_qc >= param_first_scalar ) RQLBLTEN(i,k,j)=0.
                if( p_qi >= param_first_scalar ) RQIBLTEN(i,k,j)=0.
                !if( p_qnc >= param_first_scalar ) RQNCBLTEN(i,k,j)=0.
                !if( p_qni >= param_first_scalar ) RQNIBLTEN(i,k,j)=0.
                !QKE(i,k,j)=0.
                EXCH_H(i,k,j)=0.
!                if(icloud_bl > 0) qc_bl(i,k,j)=0.
!                if(icloud_bl > 0) cldfra_bl(i,k,j)=0.
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    mynn_level=level

  END SUBROUTINE mynn_bl_init_driver

! ==================================================================

  SUBROUTINE GET_PBLH(KTS,KTE,zi,thetav1D,qke1D,zw1D,dz1D,landsea,kzi)

    !---------------------------------------------------------------
    !             NOTES ON THE PBLH FORMULATION
    !
    !The 1.5-theta-increase method defines PBL heights as the level at 
    !which the potential temperature first exceeds the minimum potential 
    !temperature within the boundary layer by 1.5 K. When applied to 
    !observed temperatures, this method has been shown to produce PBL-
    !height estimates that are unbiased relative to profiler-based 
    !estimates (Nielsen-Gammon et al. 2008). However, their study did not
    !include LLJs. Banta and Pichugina (2008) show that a TKE-based 
    !threshold is a good estimate of the PBL height in LLJs. Therefore,
    !a hybrid definition is implemented that uses both methods, weighting
    !the TKE-method more during stable conditions (PBLH < 400 m).
    !A variable tke threshold (TKEeps) is used since no hard-wired
    !value could be found to work best in all conditions.
    !---------------------------------------------------------------

    INTEGER,INTENT(IN) :: KTS,KTE



    REAL, INTENT(OUT) :: zi
    REAL, INTENT(IN) :: landsea
    REAL, DIMENSION(KTS:KTE), INTENT(IN) :: thetav1D, qke1D, dz1D
    REAL, DIMENSION(KTS:KTE+1), INTENT(IN) :: zw1D
    !LOCAL VARS
    REAL ::  PBLH_TKE,qtke,qtkem1,wt,maxqke,TKEeps,minthv
    REAL :: delt_thv   !delta theta-v; dependent on land/sea point
    REAL, PARAMETER :: sbl_lim  = 200. !upper limit of stable BL height (m).
    REAL, PARAMETER :: sbl_damp = 400. !transition length for blending (m).
    INTEGER :: I,J,K,kthv,ktke,kzi,kzi2

    !ADD KPBL (kzi)
    !KZI2 is the TKE-based part of the hybrid KPBL
    kzi = 2
    kzi2= 2

    !FIND MIN THETAV IN THE LOWEST 200 M AGL
    k = kts+1
    kthv = 1
    minthv = 9.E9
    DO WHILE (zw1D(k) .LE. 200.)
    !DO k=kts+1,kte-1
       IF (minthv > thetav1D(k)) then
           minthv = thetav1D(k)
           kthv = k
       ENDIF
       k = k+1
       !IF (zw1D(k) .GT. sbl_lim) exit
    ENDDO

    !FIND THETAV-BASED PBLH (BEST FOR DAYTIME).
    zi=0.
    k = kthv+1
    IF((landsea-1.5).GE.0)THEN
        ! WATER
        delt_thv = 0.75
    ELSE
        ! LAND
        delt_thv = 1.25
    ENDIF

    zi=0.
    k = kthv+1
!    DO WHILE (zi .EQ. 0.) 
    DO k=kts+1,kte-1
       IF (thetav1D(k) .GE. (minthv + delt_thv))THEN
          !kzi = MAX(k-1,1)
          zi = zw1D(k) - dz1D(k-1)* &
             & MIN((thetav1D(k)-(minthv + delt_thv))/ &
             & MAX(thetav1D(k)-thetav1D(k-1),1E-6),1.0)
          kzi= MAX(k-1,1) + NINT((zi-zw1D(k-1))/dz1D(k-1))
       ENDIF
       !k = k+1
       IF (k .EQ. kte-1) zi = zw1D(kts+1) !EXIT SAFEGUARD
       IF (zi .NE. 0.0) exit
    ENDDO
    !print*,"IN GET_PBLH:",thsfc,zi

    !FOR STABLE BOUNDARY LAYERS, USE TKE METHOD TO COMPLEMENT THE
    !THETAV-BASED DEFINITION (WHEN THE THETA-V BASED PBLH IS BELOW ~0.5 KM).
    !THE TANH WEIGHTING FUNCTION WILL MAKE THE TKE-BASED DEFINITION NEGLIGIBLE 
    !WHEN THE THETA-V-BASED DEFINITION IS ABOVE ~1 KM.
    ktke = 1
    maxqke = MAX(Qke1D(kts),0.)
    !Use 5% of tke max (Kosovic and Curry, 2000; JAS)
    !TKEeps = maxtke/20. = maxqke/40.
    TKEeps = maxqke/40.
    TKEeps = MAX(TKEeps,0.02) !0.025) 
    PBLH_TKE=0.

    k = ktke+1
!    DO WHILE (PBLH_TKE .EQ. 0.) 
    DO k=kts+1,kte-1
       !QKE CAN BE NEGATIVE (IF CKmod == 0)... MAKE TKE NON-NEGATIVE.
       qtke  =MAX(Qke1D(k)/2.,0.)      ! maximum TKE
       qtkem1=MAX(Qke1D(k-1)/2.,0.)
       IF (qtke .LE. TKEeps) THEN
           !kzi2 = MAX(k-1,1)
           PBLH_TKE = zw1D(k) - dz1D(k-1)* &
             & MIN((TKEeps-qtke)/MAX(qtkem1-qtke, 1E-6), 1.0)
           !IN CASE OF NEAR ZERO TKE, SET PBLH = LOWEST LEVEL.
           PBLH_TKE = MAX(PBLH_TKE,zw1D(kts+1))
           kzi2 = MAX(k-1,1) + NINT((PBLH_TKE-zw1D(k-1))/dz1D(k-1))
           !print *,"PBLH_TKE:",i,j,PBLH_TKE, Qke1D(k)/2., zw1D(kts+1)
       ENDIF
       !k = k+1
       IF (k .EQ. kte-1) PBLH_TKE = zw1D(kts+1) !EXIT SAFEGUARD
       IF (PBLH_TKE .NE. 0.) exit
    ENDDO

    !With TKE advection turned on, the TKE-based PBLH can be very large 
    !in grid points with convective precipitation (> 8 km!),
    !so an artificial limit is imposed to not let PBLH_TKE exceed the
    !theta_v-based PBL height +/- 350 m.
    !This has no impact on 98-99% of the domain, but is the simplest patch
    !that adequately addresses these extremely large PBLHs.
    PBLH_TKE = MIN(PBLH_TKE,zi+350.)
    PBLH_TKE = MAX(PBLH_TKE,MAX(zi-350.,10.))

    wt=.5*TANH((zi - sbl_lim)/sbl_damp) + .5
    IF (maxqke <= 0.05) THEN
       !Cold pool situation - default to theta_v-based def
    ELSE
       !BLEND THE TWO PBLH TYPES HERE: 
       zi=PBLH_TKE*(1.-wt) + zi*wt
    ENDIF

    !ADD KPBL (kzi) for coupling to some Cu schemes
     kzi = MAX(INT(kzi2*(1.-wt) + kzi*wt),1)



  END SUBROUTINE GET_PBLH
  
  
subroutine Poisson(istart,iend,jstart,jend,mu,POI,seed)
implicit none
integer, intent(in) :: istart,iend,jstart,jend
real,dimension(istart:iend,jstart:jend),intent(in) :: MU
integer, dimension(istart:iend,jstart:jend), intent(out) :: POI
integer,dimension(2),intent(in) :: seed
integer :: seed_len,i,j
integer,allocatable:: the_seed(:)

!if (seed .le. 0) then seed=max(-seed,1)


call random_seed(SIZE=seed_len)
allocate(the_seed(seed_len))
the_seed(1:2)=seed
! Gfortran uses longer seeds, so fill the rest with zero
if (seed_len > 2) the_seed(3:) = seed(2)
 
 
call random_seed(put=the_seed)


do i=istart,iend
 do j=jstart,jend
    poi(i,j)=poidev(mu(i,j))
    
enddo
 enddo

end subroutine Poisson

FUNCTION poidev(xm)
!USE nrtype
!USE nr, ONLY : gammln,ran1
IMPLICIT NONE
INTEGER, PARAMETER :: SP = KIND(1.0)
REAL(SP), INTENT(IN) :: xm
REAL(SP) :: poidev
REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
!Returns as a floating-point number an integer value that is a random deviate drawn from a
!Poisson distribution of mean xm, using ran1 as a source of uniform random deviates.
REAL(SP) :: em,harvest,t,y
REAL(SP), SAVE :: alxm,g,oldm=-1.0_sp,sq
!oldm is a flag for whether xm has changed since last call.
if (xm < 12.0) then !Use direct method.
if (xm /= oldm) then
oldm=xm
g=exp(-xm) !If xm is new, compute the exponential.
end if
em=-1
t=1.0
do
em=em+1.0_sp     !Instead of adding exponential deviates it is
                 !equivalent to multiply uniform deviates.
                 !We never actually have to take the log;
                 !merely compare to the pre-computed exponential.
call random_number(harvest)
t=t*harvest
if (t <= g) exit
end do
else      !    Use rejection method.
if (xm /= oldm) then  !If xm has changed since the last call, then precompute
                       !some functions that occur below.
oldm=xm
sq=sqrt(2.0_sp*xm)
alxm=log(xm)
g=xm*alxm-gammln_s(xm+1.0_sp) ! The function gammln is the natural log of the
end if                      ! gamma function, as given in §6.1.
do
do
call random_number(harvest)  !y is a deviate from a Lorentzian comparison
y=tan(PI*harvest)   !function.
em=sq*y+xm          !em is y, shifted and scaled.
if (em >= 0.0) exit !Reject if in regime of zero probability.
end do
em=int(em)          ! The trick for integer-valued distributions.
t=0.9_sp*(1.0_sp+y**2)*exp(em*alxm-gammln_s(em+1.0_sp)-g)
!The ratio of the desired distribution to the comparison function; we accept or reject
!by comparing it to another uniform deviate. The factor 0.9 is chosen so that t never
!exceeds 1.
call random_number(harvest)
if (harvest <= t) exit
end do
end if
poidev=em
END FUNCTION poidev
        
FUNCTION arth_d(first,increment,n)
implicit none
INTEGER, PARAMETER :: SP = KIND(1.0)
INTEGER, PARAMETER :: DP = KIND(1.0D0)
INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
REAL(DP), INTENT(IN) :: first,increment
INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
INTEGER(I4B), INTENT(IN) :: n
REAL(DP), DIMENSION(n) :: arth_d
INTEGER(I4B) :: k,k2
REAL(DP) :: temp
if (n > 0) arth_d(1)=first
if (n <= NPAR_ARTH) then
do k=2,n
arth_d(k)=arth_d(k-1)+increment
end do
else
do k=2,NPAR2_ARTH
arth_d(k)=arth_d(k-1)+increment
end do
temp=increment*NPAR2_ARTH
k=NPAR2_ARTH
do
if (k >= n) exit
k2=k+k
arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
temp=temp+temp
k=k2
end do
end if
END FUNCTION arth_d
      
FUNCTION gammln_s(xx)
IMPLICIT NONE
INTEGER, PARAMETER :: SP = KIND(1.0)
INTEGER, PARAMETER :: DP = KIND(1.0D0)
REAL(SP), INTENT(IN) :: xx
REAL(SP) :: gammln_s
!Returns the value ln[Γ(xx)] for xx > 0.
REAL(DP) :: tmp,x
!Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
!accuracy is good enough.
REAL(DP) :: stp = 2.5066282746310005_dp
REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
-86.50532032941677_dp,24.01409824083091_dp,&
-1.231739572450155_dp,0.1208650973866179e-2_dp,&
-0.5395239384953e-5_dp/)
!call assert(xx > 0.0, ’gammln_s arg’)
if (xx .le. 0.) print *,'gammaln fails'
x=xx
tmp=x+5.5_dp
tmp=(x+0.5_dp)*log(tmp)-tmp
gammln_s=tmp+log(stp*(1.000000000190015_dp+&
sum(coef(:)/arth_d(x+1.0_dp,1.0_dp,size(coef))))/x)
END FUNCTION gammln_s
  

subroutine condensation_edmf(QT,THL,P,lf,THV,QC)
!
real,intent(in)   :: QT,THL,P,lf
real,intent(out)  :: THV
real,intent(inout):: QC

integer :: niter,i
real :: diff,exn,t,th,qs,qcold,xlvv

! constants used from module_model_constants.F
! p1000mb
! rcp ... Rd/cp
! xlv ... latent heat for water (2.5e6)
! cp
! rvord .. rv/rd  (1.6) 


! number of iterations
  niter=50
! minimum difference
  diff=2.e-5


  EXN=(P/p1000mb)**rcp

  xlvv=xlLF_blend(THL*EXN,lf)

  !QC=0.  !better first guess QC is incoming from lower level, do not set to zero

  do i=1,NITER
     T=EXN*THL + xlvv/cp*QC        
     QS=qsatLF_blend(T,P,lf)
     QCOLD=QC
     QC=0.5*QC + 0.5*MAX((QT-QS),0.)
     if (abs(QC-QCOLD)<Diff) exit
  enddo

  T=EXN*THL + xlvv/cp*QC
  QS=qsatLF_blend(T,P,lf)
  QC=max(QT-QS,0.)

  THV=(THL+xlvv/cp*QC)*(1.+QT*(rvovrd-1.)-rvovrd*QC)

end subroutine condensation_edmf



subroutine condensation_edmfA(THV,QT,P,LF,THL,QC)
!
! zero or one condensation for edmf: calculates QC and THL from THV and QT
!


real,intent(in) :: THV,QT,P,LF
real,intent(out):: THL,QC


integer :: niter,i
real :: exn,t,qs,qcold,xlvv
 
! max number of iterations (we dont need that many because th~thv)
niter=2
! minimum difference
!diff=2.e-5

EXN=(P/p1000mb)**rcp
QC=0. 


do i=1,NITER
   T=EXN*THV/(1.+QT*(rvovrd-1.)-rvovrd*QC)
   QS=qsatLF_blend(T,P,lf)
   QCOLD=QC
   QC=max(0.5*QC+0.5*(QT-QS),0.)
!if (abs(QC-QCOLD)<Diff) exit
enddo

 xlvv=xlLF_blend(THL*EXN,lf)

 THL=(T-xlvv/cp*QC)/EXN
 

end subroutine condensation_edmfA


!===============================================================

! subroutine condensation_edmf_r(QT,THL,P,zagl,THV,QC)
! !
! ! zero or one condensation for edmf: calculates THL and QC
! ! similar to condensation_edmf but with different inputs
! !
! real,intent(in)   :: QT,THV,P,zagl
! real,intent(out)  :: THL, QC
! 
! integer :: niter,i
! real :: diff,exn,t,th,qs,qcold
! 
! ! number of iterations
!   niter=50
! ! minimum difference
!   diff=2.e-5
! 
!   EXN=(P/p1000mb)**rcp
!   ! assume first that th = thv
!   T = THV*EXN
!   !QS = qsat_blend(T,P)
!   !QC = QS - QT
!   
!   QC=0.
! 
!   do i=1,NITER    
!      QCOLD = QC
!      T = EXN*THV/(1.+QT*(rvovrd-1.)-rvovrd*QC)   
!      QS=qsat_blend(T,P)
!      QC= MAX((QT-QS),0.)
!      if (abs(QC-QCOLD)<Diff) exit
!   enddo
!   THL = (T - xlv/cp*QC)/EXN
! 
! end subroutine condensation_edmf_r
! =====================================================================

  FUNCTION esatLF_blend(t,lf) 
!
! This calculates saturation vapor pressure.  Separate ice and liquid functions 
! are used (identical to those in module_mp_thompson.F, v3.6).  Then, the 
! final returned value is a temperature-dependant "blend".  Because the final 
! value is "phase-aware", this formulation may be preferred for use throughout 
! the module (replacing "svp").
!
      IMPLICIT NONE
      
      REAL, INTENT(IN):: t,lf
      REAL :: esatLF_blend,XC,ESL,ESI,chi

      XC=MAX(-80.,t-273.16)

! For 253 < t < 273.16 K, the vapor pressures are "blended" as a function of temperature, 
! using the approach of Chaboureau and Bechtold (2002), JAS, p. 2363.  The resulting 
! values are returned from the function.
      IF (lf .GE. 1.) THEN
          esatLF_blend = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8))))))) 
      ELSE IF (lf .LE. 0.) THEN
          esatLF_blend = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
      ELSE
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          esatLF_blend = lf*ESL  + (1.-lf)*ESI
      END IF

  END FUNCTION esatLF_blend



!   FUNCTION esat_blend(t) 
! ! JAYMES- added 22 Apr 2015
! ! 
! ! This calculates saturation vapor pressure.  Separate ice and liquid functions 
! ! are used (identical to those in module_mp_thompson.F, v3.6).  Then, the 
! ! final returned value is a temperature-dependant "blend".  Because the final 
! ! value is "phase-aware", this formulation may be preferred for use throughout 
! ! the module (replacing "svp").
! 
!       IMPLICIT NONE
!       
!       REAL, INTENT(IN):: t
!       REAL :: esat_blend,XC,ESL,ESI,chi
! 
!       XC=MAX(-80.,t-273.16)
! 
! ! For 253 < t < 273.16 K, the vapor pressures are "blended" as a function of temperature, 
! ! using the approach of Chaboureau and Bechtold (2002), JAS, p. 2363.  The resulting 
! ! values are returned from the function.
!       IF (t .GE. 273.16) THEN
!           esat_blend = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8))))))) 
!       ELSE IF (t .LE. 253.) THEN
!           esat_blend = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
!       ELSE
!           ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
!           ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
!           chi  = (273.16-t)/20.16
!           esat_blend = (1.-chi)*ESL  + chi*ESI
!       END IF
! 
!   END FUNCTION esat_blend
! 
! ! ====================================================================

!   FUNCTION qsat_blend(t, P)
! ! JAYMES- this function extends function "esat" and returns a "blended"
! ! saturation mixing ratio.
! 
!       IMPLICIT NONE
! 
!       REAL, INTENT(IN):: t, P
!       REAL :: qsat_blend,XC,ESL,ESI,RSLF,RSIF,chi
! 
!       XC=MAX(-80.,t-273.16)
! 
!       IF (t .GE. 273.16) THEN
!           ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8))))))) 
!           qsat_blend = 0.622*ESL/(P-ESL) 
!       ELSE IF (t .LE. 253.) THEN
!           ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
!           qsat_blend = 0.622*ESI/(P-ESI)
!       ELSE
!           ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
!           ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
!           RSLF = 0.622*ESL/(P-ESL)
!           RSIF = 0.622*ESI/(P-ESI)
!           chi  = (273.16-t)/20.16
!           qsat_blend = (1.-chi)*RSLF + chi*RSIF
!       END IF
! 
!   END FUNCTION qsat_blend

  FUNCTION qsatLF_blend(t, P,lf)
! JAYMES- this function extends function "esat" and returns a "blended"
! saturation mixing ratio.

      IMPLICIT NONE

      REAL, INTENT(IN):: t, P,lf
      REAL :: qsatLF_blend,XC,ESL,ESI,RSLF,RSIF,chi

      XC=MAX(-80.,t-273.16)

      IF (lf .GE. 1.) THEN
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8))))))) 
          qsatLF_blend = 0.622*ESL/(P-ESL) 
      ELSE IF (lf .LE. 0.) THEN
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          qsatLF_blend = 0.622*ESI/(P-ESI)
      ELSE
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          RSLF = 0.622*ESL/(P-ESL)
          RSIF = 0.622*ESI/(P-ESI)
          !chi  = (273.16-t)/20.16
          qsatLF_blend = lf*RSLF + (1.-lf)*RSIF
      END IF

  END FUNCTION qsatLF_blend

! ===================================================================

!   FUNCTION xl_blend(t)
! ! JAYMES- this function interpolates the latent heats of vaporization and
! ! sublimation into a single, temperature-dependant, "blended" value, following
! ! Chaboureau and Bechtold (2002), Appendix.
! 
!       IMPLICIT NONE
! 
!       REAL, INTENT(IN):: t
!       REAL :: xl_blend,xlvt,xlst,chi
! 
!       IF (t .GE. 273.16) THEN
!           xl_blend = xlv + (cpv-cliq)*(t-273.16)  !vaporization/condensation
!       ELSE IF (t .LE. 253.) THEN
!           xl_blend = xls + (cpv-cice)*(t-273.16)  !sublimation/deposition
!       ELSE
!           xlvt = xlv + (cpv-cliq)*(t-273.16)  !vaporization/condensation
!           xlst = xls + (cpv-cice)*(t-273.16)  !sublimation/deposition
!           chi  = (273.16-t)/20.16
!           xl_blend = (1.-chi)*xlvt + chi*xlst     !blended
!       END IF
! 
!   END FUNCTION xl_blend


  FUNCTION xlLF_blend(t,lf)
! JAYMES- this function interpolates the latent heats of vaporization and
! sublimation into a single, temperature-dependant, "blended" value, following
! Chaboureau and Bechtold (2002), Appendix.

      IMPLICIT NONE

      REAL, INTENT(IN):: t,lf
      REAL :: xlLF_blend,xlvt,xlst,chi

      IF (lf .GE. 1.) THEN
          xlLF_blend = xlv + (cpv-cliq)*(t-273.16)  !vaporization/condensation
      ELSE IF (lf .LE. 0.) THEN
          xlLF_blend = xls + (cpv-cice)*(t-273.16)  !sublimation/deposition
      ELSE
          xlvt = xlv + (cpv-cliq)*(t-273.16)  !vaporization/condensation
          xlst = xls + (cpv-cice)*(t-273.16)  !sublimation/deposition
          xlLF_blend = lf*xlvt + (1.-lf)*xlst     !blended
      END IF

  END FUNCTION xlLF_blend

! ===================================================================
! ===================================================================
! This is the original mass flux scheme from Kay Suselj of NASA-JPL,


SUBROUTINE edmf_JPL(kts,kte,dt,zw,p,         &
              !&u,v,th,thl,thv,tk,qt,qv,qc,   &  yhc comment out
              &u,v,th,thl,thv,tk,qt,qv,qc,ql,qi,   &  ! yhc add ql and qi
              &exner,rho,                        &
              &ust,ps,wthl,wqt,pblh,            &
              & liquid_frac,                 &
            ! outputs - updraft properties : new implementation
              & edmf_a,edmf_w, edmf_qt,      &
              & edmf_thl,edmf_ent,edmf_qc,   &
              & edmf_debug1,edmf_debug2,     &
              & edmf_debug3,edmf_debug4,     &
            ! outputs - variables needed for solver 
              & s_aw,s_awthl,s_awqt,         &
              & s_awqv,s_awqc,s_awu,s_awv,   &
              & s_awqke,                     &
            ! in/outputs - subgrid scale clouds
              & qc_bl1d,cldfra_bl1d,         &
            ! output info
              & a_moist_half, mf_moist_half, qv_moist_half, a_moist_full, mf_moist_full, qv_moist_full, &  ! yhc 2021-09-08
              & a_dry_half, mf_dry_half, qv_dry_half, a_dry_full, mf_dry_full, qv_dry_full, &            ! yhc 2021-09-08
              & mf_all_half, mf_all_full, edmf_det, &                                                  ! yhc 2021-09-08
              & num_updraft, num_DET, num_nDET_zENT, num_nDET_pENT, &                                               ! yhc 2021-09-08
              &ktop,ztop,kpbl,Qql,Qqi,Qa, &
              & streams1, &    ! yhc 2021-11-18
              & Qql_adv,Qqi_adv,Qa_adv, Qql_eddy,Qqi_eddy,Qa_eddy, Qql_ent,Qqi_ent,Qa_ent, Qql_det,Qqi_det,Qa_det, Qql_sub,Qqi_sub,Qa_sub &
              ) ! yhc 2021-09-08



        INTEGER, INTENT(IN) :: KTS,KTE, kpbl
        REAL,DIMENSION(KTS:KTE), INTENT(IN) :: U,V,TH,THL,TK,QT,QV,QC, QL, QI  ! yhc add QL, QI
        REAL,DIMENSION(KTS:KTE), INTENT(IN) :: THV,P,exner,rho,liquid_frac
        ! zw .. heights of the updraft levels (edges of boxes)
        REAL,DIMENSION(KTS:KTE+1), INTENT(IN) :: ZW
        REAL, INTENT(IN) :: DT,UST,PS,WTHL,WQT,PBLH


   !     LOGICAL :: cloudflg
    ! outputs
    !     REAL,DIMENSION(KTS:KTE), INTENT(OUT) :: DTH,DQV,DQC,DU,DV

        ! REAL,DIMENSION(KTS:KTE) :: dry_a, moist_a,dry_w,moist_w, &
        !        dry_qt,moist_qt,dry_thl,moist_thl, dry_ent, moist_ent, moist_qc
    ! outputs - updraft properties
        REAL,DIMENSION(KTS:KTE), INTENT(OUT) :: edmf_a,edmf_w,        &
                      & edmf_qt,edmf_thl, edmf_ent,edmf_qc, &
                      & edmf_debug1, edmf_debug2, edmf_debug3, edmf_debug4
    ! output
        INTEGER, INTENT(OUT) :: ktop
        REAL, INTENT(OUT) :: ztop

    ! outputs - variables needed for solver (s_aw - sum ai*wi, s_awphi - sum ai*wi*phii)
        REAL,DIMENSION(KTS:KTE+1) :: s_aw, s_awthl, s_awqt, s_awu, s_awv, s_awqc, s_awqv, s_awqke, s_aw2
        REAL,DIMENSION(KTS:KTE), INTENT(IN) :: qc_bl1d, cldfra_bl1d

      REAL,DIMENSION(KTS:KTE), INTENT(OUT) :: Qql,Qqi,Qa, &
        Qql_det,Qqi_det,Qa_det, Qql_sub,Qqi_sub,Qa_sub, &
        Qql_adv,Qqi_adv,Qa_adv, Qql_eddy,Qqi_eddy,Qa_eddy, Qql_ent,Qqi_ent,Qa_ent

        !INTEGER, PARAMETER :: NUP=100, debug_mf=0 !fixing number of plumes to 10
        INTEGER, PARAMETER :: debug_mf=0 !fixing number of plumes to 10, yhc - move NUP to namelist parameter

    ! updraft properties
        REAL,DIMENSION(KTS:KTE+1,1:NUP) :: UPW,UPTHL,UPQT,UPQC,UPA,UPU,UPV,UPTHV,UPRHO

    ! entrainment variables
        REAl,DIMENSION(KTS:KTE,1:NUP) :: ENT,ENTf
        REAL,DIMENSION(KTS:KTE) :: L0s
        INTEGER,DIMENSION(KTS:KTE,1:NUP) :: ENTi

        REAl,DIMENSION(KTS:KTE,1:NUP) :: DET
        REAL,DIMENSION(KTS:KTE), INTENT(out) :: edmf_det

    ! internal variables
        INTEGER :: K,I, qlTop
        REAL :: wthv,wstar,qstar,thstar,sigmaW,sigmaQT,sigmaTH,z0, &
            pwmin,pwmax,wmin,wmax,wlv,wtv,dthv_dz
        REAL :: B,QTn,THLn,THVn,UPAn,QCn,Un,Vn,Wn2,EntEXP,EntW, deltaZ,EntExp_M, Beta_un, Z00, Z_i

    ! VARIABLES FOR CHABOUREAU-BECHTOLD CLOUD FRACTION
        ! REAL,DIMENSION(KTS:KTE), INTENT(INOUT) :: vt, vq, sgm
        REAL :: sigq,xl,tlk,qsat_tl,rsl,cpm,a,qmq,mf_cf,diffqt,&
               Fng,qww,alpha,beta,bb,f,pt,t,q2p,b9,satvp,rhgrid
        
        REAL,DIMENSION(kts:kte) :: ph,lfh 
        REAL :: &
          F1_ql, F1_qi, F2_ql, F2_qi, F3_ql, F3_qi, &
          THVsrfF,QTsrfF,maxS,stabF,dz,F0,F1,F2,F3,CCp1,CCp0,mf,mfp1  
        INTEGER, DIMENSION(2) :: seedmf
    ! w parameters
        REAL,PARAMETER :: &
            &Wa=1., &
            &Wb=1.5

    ! entrainment parameters
        REAL,PARAMETER :: &
        !& L0=100.,&    yhc move L0 to namelist paramter 
        & ENT0=0.2
       ! maximum wind speed in updraft
       REAL,PARAMETER :: MAXW=2.  

    !<--- yhc 2021-09-08
    !  naming convention: 
    !     a: fractional area (none)
    !    mf: mass flux (kg/m2/s)
    !    qv: specific humifity (kg/kg)
    !    moist/dry/all: moist updrafts, dry updrafts, or combined (moist+dry) updrafts
    !    full/half: on full levels (where T,q are defined) to half levels (those levels in between full levels)
       REAL,DIMENSION(kts:kte+1), INTENT(OUT) :: & 
         num_updraft  , &
         a_moist_half , &   
         a_dry_half   , &        
         mf_moist_half, &     
         mf_dry_half  , &
         mf_all_half  , & 
         qv_moist_half, &     
         qv_dry_half  

       REAL,DIMENSION(kts:kte), INTENT(OUT) :: & 
         num_DET, num_nDET_zENT, num_nDET_pENT,  &
         a_moist_full ,  &
         a_dry_full ,  &
         mf_moist_full,  &
         mf_dry_full  ,  & 
         mf_all_full  , & 
         qv_moist_full,  & 
         qv_dry_full   
  
       REAL,DIMENSION(KTS:KTE) :: &
         ZFULL, ice_frac

       REAL :: &
         qcp0, qcp1, qlp0, qlp1, qip0, qip1, cldp1, cldp0, det_temp

       REAl,DIMENSION(KTS:KTE,1:NUP) :: &
        Qql_det_i, Qqi_det_i, Qa_det_i, Qql_sub_i , Qqi_sub_i , Qa_sub_i, &
        Qql_adv_i, Qqi_adv_i, Qa_adv_i, Qql_eddy_i, Qqi_eddy_i, Qa_eddy_i, Qql_ent_i, Qqi_ent_i, Qa_ent_i

       integer, parameter :: rx = 500  ! AM4 random number generator
       real :: rr(rx)

       logical :: do_MF

       type(randomNumberStream), INTENT(INOUT) :: streams1
    !---> yhc 2021-09-08 

! stability parameter for massflux
! (mass flux is limited so that dt/dz*a_i*w_i<UPSTAB)
!       REAL,PARAMETER :: UPSTAB=1.   ! yhc - move UPSTAB to namelist parameter

!Initialize values:
ktop = 0
ztop = 0.0

UPW=0.
UPRHO=0.
UPTHL=0.
UPTHV=0.
UPQT=0.
UPA=0.
UPU=0.
UPV=0.
UPQC=0.
ENT=0.
L0s=100.

wthv=wthl+svp1*wqt

! Initialize mean updraft properties
edmf_a  =0.
edmf_w  =0.
edmf_qt =0.
edmf_thl=0.
edmf_ent=0.
edmf_qc =0.

edmf_debug1 = 0.
edmf_debug2 = 0.
edmf_debug3 = 0.
edmf_debug4 = 0.

s_aw=0.
s_awthl=0.
s_awqt=0.
s_awqv=0.
s_awqc=0.
s_awu=0.
s_awv=0.
s_awqke=0.


 Qql=0.
 Qqi=0. 
 Qa=0. 

 Qql_adv  = 0.
 Qqi_adv  = 0.
 Qa_adv  = 0.

 Qql_eddy = 0.
 Qqi_eddy = 0.
 Qa_eddy = 0.

 Qql_ent  = 0.
 Qqi_ent  = 0.
 Qa_ent  = 0.

 Qql_det  = 0.
 Qqi_det  = 0.
 Qa_det  = 0.

 Qql_sub  = 0.
 Qqi_sub  = 0.
 Qa_sub  = 0.

 Qql_adv_i = 0.
 Qqi_adv_i = 0.
 Qa_adv_i = 0.

 Qql_eddy_i = 0.
 Qqi_eddy_i = 0.
 Qa_eddy_i = 0.

 Qql_ent_i = 0.
 Qqi_ent_i = 0.
 Qa_ent_i = 0.

 Qql_det_i = 0.
 Qqi_det_i = 0.
 Qa_det_i = 0.

 Qql_sub_i = 0.
 Qqi_sub_i = 0.
 Qa_sub_i = 0.

 edmf_det = 0.
 DET=0.

 ice_frac(:) = 1. - liquid_frac(:)

! This part is specific for Stratocumulus
! cloudflg = .false.
!DO k = MAX(1,kpbl-2),kpbl+5
!    IF(qc(k).gt. 1.e-6 .AND. cldfra_bl1D(k).gt.0.5) THEN
!       cloudflg=.true. !found Sc cloud
!       qlTop = k ! index for Sc cloud top
!    ENDIF
!ENDDO
!Z00 = 1000.
!IF (cloudflg) THEN ! if Sc exist, slows down updraft near inversion
!    Z_i = ZW(qlTop+1) ! let updraft go above cloud top a little
!ELSE
!    Z_i = ZW(kte) ! otherwise set to model top so it does not affect updraft w
!ENDIF
! print*, "[EDMF_JPL]: Sc found at z = ", Z_i
! if surface buoyancy is positive we do integration otherwise not

    !<-- yhc 2021-11-26
    z0=50.
    do_MF = .true.
    if (option_pblh_MF.eq.1 .and. pblh < z0) then  ! only when w'thv'>0 and PBL height > z0 (=50m), call MF. 
                                                   ! This is to prevent activating MF in very stable conditions.
      do_MF = .false.
    endif
    !--> yhc 2021-11-26

IF ( wthv >= 0.0 .and. do_MF) then  ! yhc 2021-11-26 add
!IF ( wthv >= 0.0 ) then  ! yhc 2021-11-26 comment out

! get the pressure and liquid_fraction on updraft levels

  lfH(KTS)=liquid_frac(KTS)
  pH(KTS)=ps 
  DO K=KTS+1,KTE
   lfH(K)=0.5*(liquid_frac(K)+liquid_frac(K-1))
    ph(K)=0.5*(p(K)+p(K-1)) 
  ENDDO




 ! set initial conditions for updrafts
 !   z0=50.  ! yhc 2021-11-26 comment out. z0 is set outside of the wthv if-statement
    pwmin=1.
    pwmax=3.

    wstar=max(0.,(g/thv(KTS)*wthv*pblh)**(1./3.))
    qstar=wqt/wstar
    thstar=wthl/wstar
    sigmaW=1.34*wstar*(min(z0/pblh,1.))**(1./3.)*(1-0.8*min(z0/pblh,1.))  ! force z0/pblh<1. to avoid negative sigmaW and negative w
    sigmaQT=1.34*qstar*(min(z0/pblh,1.))**(-1./3.)
    sigmaTH=1.34*thstar*(min(z0/pblh,1.))**(-1./3.)
    !sigmaW=1.34*wstar*(z0/pblh)**(1./3.)*(1-0.8*z0/pblh)
    !sigmaQT=1.34*qstar*(z0/pblh)**(-1./3.)
    !sigmaTH=1.34*thstar*(z0/pblh)**(-1./3.)


    wmin=sigmaW*pwmin
    wmax=sigmaW*pwmax

    ! set entrainment L0 to be a function of dthv/dz, height, and w_star
    ! do k=kts+1,kte
    !     dthv_dz = 1./(ZW(k)-ZW(k-1))*(thv(k)-thv(k-1))
    !     if (dthv_dz>0.) then
    !         L0s(k) = L0s(k) - wstar**2*300./dthv_dz/g/(Z_i-ZW(k))
    !     endif
    !     if (L0s(k)<20.) then
    !         L0s(k) = 20.
    !     elseif (L0s(k)>100.) then
    !         L0s(k) = 100.
    !     endif
    !     print *, "EDMF_JPL L0s at k = ", k, " z = ", ZW(k), " L0 = ", L0s(k)
    ! enddo

!
! get entrainment coefficient
!

! get dz/L0
    do i=1,Nup
        do k=kts,kte
           ENTf(k,i)=(ZW(k+1)-ZW(k))/L0! s(k)
! Elynn's modification of entrainment rate
!          if (ZW(k) < Z_i) then
!            ENTf(k,i) = (ZW(k)-ZW(k-1)) * (3.5 * exp((-(Z_i-ZW(k)))/70.) + 1./L0)
!          else
!            ENTf(k,i) = (ZW(k)-ZW(k-1)) * (3.5 * exp(1./70.) + 1./L0)
!          endif
         enddo
    enddo


! create seed for Poisson
! create seed from last digits of temperature   
seedmf(1) = 1000000 * ( 100*thl(1) - INT(100*thl(1)))
seedmf(2) = 1000000 * ( 100*thl(2) - INT(100*thl(2))) 

! get Poisson P(dz/L0)
!    call Poisson(kts,kte,1,Nup,ENTf,ENTi,seedmf)

  !<-- yhc 2021-11-18, use AM4 random number geneartor
  if (option_stoch_entrain.eq."Poisson") then  
    call Poisson(kts,kte,1,Nup,ENTf,ENTi,seedmf)   ! the original Possion function

  elseif (option_stoch_entrain.eq."Poisson_knuth") then   ! use AM4 random number geneartor
    if (option_rng == 0) call random_number(rr)
    if (option_rng == 1) call getRandomNumbers (streams1,rr)
    call Poisson_knuth (kte-kts+1, Nup, rx, rr, ENTf, ENTi)

  endif    ! end if of option_stoch_entrain.eq
  !-->


 ! entrainent: Ent=Ent0/dz*P(dz/L0)
    do i=1,Nup
        do k=kts,kte
          if (option_ent.eq.1) then   ! original stochstic entrainment. But when it produce zero entrainment with constant
                                      ! area assumption, this can violate mass continuity in plumes. 
            ENT(k,i)=real(ENTi(k,i))*Ent0/(ZW(k+1)-ZW(k))

          elseif (option_ent.eq.2) then   ! separate Ent0 to stochastic part and deterministic part to avoid zero entrainment problems
            ENT(k,i)=alpha_st*real(ENTi(k,i))*Ent0/(ZW(k+1)-ZW(k)) + (1.-alpha_st)*Ent0/L0

          endif  ! end if of option_ent
        enddo
    enddo


    DO I=1,NUP
        wlv=wmin+(wmax-wmin)/NUP*(i-1)
        wtv=wmin+(wmax-wmin)/NUP*i

        UPW(KTS,I)=min(0.5*(wlv+wtv),maxw)
        UPA(KTS,I)=0.5*ERF(wtv/(sqrt(2.)*sigmaW))-0.5*ERF(wlv/(sqrt(2.)*sigmaW))

        UPU(KTS,I)=U(KTS)
        UPV(KTS,I)=V(KTS)

       ! UPQC(KTS,I)=0 ! computed with condensation routine
        UPQT(KTS,I)=QT(KTS)+0.32*UPW(KTS,I)*sigmaQT/sigmaW 
        UPTHV(KTS,I)=THV(KTS)+0.58*UPW(KTS,I)*sigmaTH/sigmaW
       ! UPTHL(KTS,I)=UPTHV(KTS,I)/(1.+svp1*UPQT(KTS,I)) ! with condensation routine
    ENDDO



!
! make sure that the surface flux of THV and QT does not exceed given flux
!

   QTsrfF=0.
   THVsrfF=0.
   
   DO I=1,NUP
     QTsrfF=QTsrfF+UPW(KTS,I)*UPA(KTS,I)*(UPQT(KTS,I)-QT(KTS))
     THVsrfF=THVsrfF+UPW(KTS,I)*UPA(KTS,I)*(UPTHV(KTS,I)-THV(KTS))   
   ENDDO
  
   
   IF (THVsrfF .gt. wthv) THEN
   ! change surface THV so that the fluxes from the mass flux equal prescribed values
       DO I=1,NUP
        UPTHV(KTS,I)=(UPTHV(KTS,I)-THV(KTS))*wthv/THVsrfF+THV(KTS)
       ENDDO
  ENDIF     
      
   IF ( QTsrfF .gt. wqt)  THEN
   ! change surface QT so that the fluxes from the mass flux equal prescribed values
       DO I=1,NUP
        UPQT(KTS,I)=(UPQT(KTS,I)-QT(KTS))*wqt/QTsrfF+QT(KTS)
        ENDDO
    ENDIF     


! surface condensation (compute THL and QC)
  DO I=1,NUP  
  
        UPRHO(KTS,I)=Ph(KTS)/(r_d*UPTHV(KTS,I)*(Ph(kts)/p1000mb)**rcp)
        call condensation_edmfA(UPTHV(KTS,I),UPQT(KTS,I),ph(KTS),lfh(KTS),UPTHL(KTS,I),UPQC(KTS,i))
 ENDDO
   
   
   



  ! do integration  updraft
    DO I=1,NUP
        DO k=KTS+1,KTE
            deltaZ = ZW(k)-ZW(k-1)
            EntExp=exp(-ENT(K-1,I)*deltaZ)
            EntExp_M=exp(-ENT(K-1,I)/3.*deltaZ)

            QTn=QT(K-1)*(1-EntExp)+UPQT(K-1,I)*EntExp
            THLn=THL(K-1)*(1-EntExp)+UPTHL(K-1,I)*EntExp
            Un=U(K-1)*(1-EntExp)+UPU(K-1,I)*EntExp_M
            Vn=V(K-1)*(1-EntExp)+UPV(K-1,I)*EntExp_M
            ! get thvn,qcn

            if (K.eq.KTS+1) QCn = 0.  ! yhc, first guess of QCn. If not set, Qcn may be randomly set (e.g. 0.7) which wouldcause problem
            call condensation_edmf(QTn,THLn,Ph(K),lfh(k),THVn,QCn)

            B=g*(0.5*(THVn+UPTHV(k-1,I))/THV(k-1)-1)
 !           IF (ZW(k) <= Z_i) THEN
 !               Beta_un = Wb*ENT(K,I) + 0.5/(Z_i-ZW(k)+deltaZ) * max(1. - exp((Z_i-ZW(k)+deltaZ)/Z00-1.) , 0.)
 !           ELSE
 !               Beta_un = Wb*ENT(K,I) + 0.5/deltaZ * (1. - exp(deltaZ/Z00-1.))
 !           END IF
            Beta_un=wb*ENT(K-1,I)   
          
            IF (Beta_un>0) THEN
                EntW = exp(-2.*Wb*ENT(K,I)*deltaZ)
                Wn2=UPW(K-1,I)**2*EntW+(1.-EntW)*Wa*B/Beta_un
            ELSE
                Wn2=UPW(K-1,I)**2+2*Wa*B*deltaZ
            END IF
            ! print *, "Plume number = ", I, " k = ", k, " upthv = ", THVn, " thv = ", thv(k)
            ! print *, "upthl = ", THLn, " thl = ", thl(k)
            ! print *, "Beta = ", Beta_un, " Wn2 = ", Wn2, " Bouy = ", B
         
         
             UPRHO(K,I)=Ph(k)/(r_d*THVn*(Ph(k)/p1000mb)**rcp)
         
            IF (Wn2 >0.) THEN     
              UPW(K,I)=sqrt(Wn2)

              !yhc, compute detrainment rate, assuming updraft area does not change
              mf  =UPRHO(K-1,I)*UPW(K-1,I)*UPA(K-1,I)
              mfp1=UPRHO(K,I)  *UPW(K,I)  *UPA(K-1,I)   
              det_temp=0.
              if (mf.gt.0) then
                det_temp = ENT(K-1,I) - (mfp1-mf)/mf/deltaZ
              endif
         
             ! 1/M dM/dz=eps ... conservation equation for entraining plumes with 0 detrainment
             ! integrate conservation equation from k-1 to k to get
             ! M(k)=M(k-1)*exp(int_(k-1)^k  eps*dz)
         
               !UPAn=UPA(K-1,I)*UPRHO(K-1,I)*UPW(k-1,I)/(UPRHO(K,I)*UPW(K,I)*EntExp)
         
               if (option_up_area.eq.1) then   ! updraft area is constant
                 UPAn=UPA(K-1,I)

               elseif (option_up_area.eq.2) then   ! updraft area varies with height to satisfy the mass continuity equation
                 UPAn=UPA(K-1,I)*UPRHO(K-1,I)*UPW(k-1,I)/(UPRHO(K,I)*UPW(K,I)*EntExp)

               elseif (option_up_area.eq.3) then   ! if detrainment rate>=0, assume updraft area is constant.
                                                   ! otherwise, vary the area to satisfy the mass continuity equation
                 if (det_temp.ge.0.) then
                   UPAn=UPA(K-1,I)
                 else
                   UPAn=UPA(K-1,I)*UPRHO(K-1,I)*UPW(k-1,I)/(UPRHO(K,I)*UPW(K,I)*EntExp)
                 endif
               endif

!print*,'i,k-1,upa',i,k-1,UPA(K-1,I)

               ! exit vertical loop
              !  - this should never happen, because of the exponential MF form
         !      IF (UPAn<= 0.) THEN
         !      UPRHO(K,I)=0.
         !        UPW(K,I)=0.
         !        exit
         !      ENDIF
         
                
                
           
                UPTHV(K,I)=THVn
                UPTHL(K,I)=THLn
                UPQT(K,I)=QTn
                UPQC(K,I)=QCn
                UPU(K,I)=Un
                UPV(K,I)=Vn
                UPA(K,I)=UPAn

                ktop = MAX(ktop,k)
            ELSE
            ! exit vertical loop
                  exit
            END IF
        ENDDO
    ENDDO

END IF

ktop=MIN(ktop,KTE-1)  !  Just to be safe...
IF (ktop == 0) THEN
    ztop = 0.0
ELSE
    ztop=zw(ktop)
ENDIF

!
!  get updraft properties, for saving
!

! dry_a=0.
! moist_a=0.

! dry_w=0.
! moist_w=0.

! dry_qt=0.
! moist_qt=0.

! dry_thl=0.
! moist_thl=0.

! dry_ent=0.
! moist_ent=0.

! moist_qc=0.


  ! writing updraft properties in their variable
  ! all variables, except Areas are now multipled by the area
  ! to confirm with WRF grid setup we do not save the first and the last row

! Splitting moist and dry plumes, commented to match v4 fields
! DO k=KTS,KTE-1
!     DO I=1,NUP
!         IF (UPQC(K+1,I)>0.) THEN
!             moist_a(K)=moist_a(K)+UPA(K+1,I)
!             moist_w(K)=moist_w(K)+UPA(K+1,I)*UPW(K+1,I)
!             moist_qt(K)=moist_qt(K)+UPA(K+1,I)*UPQT(K+1,I)
!             moist_thl(K)=moist_thl(K)+UPA(K+1,I)*UPTHL(K+1,I)
!             moist_ent(K)=moist_ent(K)+UPA(K+1,I)*ENT(K+1,I)
!             moist_qc(K)=moist_qc(K)+UPA(K+1,I)*UPQC(K+1,I)
!         ELSE
!             dry_a(K)=dry_a(K)+UPA(K+1,I)
!             dry_w(K)=dry_w(K)+UPA(K+1,I)*UPW(K+1,I)
!             dry_qt(K)=dry_qt(K)+UPA(K+1,I)*UPQT(K+1,I)
!             dry_thl(K)=dry_thl(K)+UPA(K+1,I)*UPTHL(K+1,I)
!             dry_ent(K)=dry_ent(K)+UPA(K+1,I)*ENT(K+1,I)
!         ENDIF
!    ! lop for n
!     ENDDO

!     IF (dry_a(k)>0.) THEN
!         dry_w(k)=dry_w(k)/dry_a(k)
!         dry_qt(k)=dry_qt(k)/dry_a(k)
!         dry_thl(k)=dry_thl(k)/dry_a(k)
!         dry_ent(k)=dry_ent(k)/dry_a(k)
!     ENDIF

!     IF (moist_a(k)>0.) THEN
!         moist_w(k)=moist_w(k)/moist_a(k)
!         moist_qt(k)=moist_qt(k)/moist_a(k)
!         moist_thl(k)=moist_thl(k)/moist_a(k)
!         moist_ent(k)=moist_ent(k)/moist_a(k)
!         moist_qc(k)=moist_qc(k)/moist_a(k)
!     ENDIF
! ENDDO

! Combine both moist and dry plume, write as one averaged plume
DO k=KTS,KTE-1
    IF(k > KTOP) exit

    DO I=1,NUP
        IF(I > NUP) exit
        ! yhc: because edmf_X variables are written to half levels in the history files, UP* variables should use K
        edmf_a(K)=edmf_a(K)+UPA(K,I)
        edmf_w(K)=edmf_w(K)+UPA(K,I)*UPW(K,I)
        edmf_qt(K)=edmf_qt(K)+UPA(K,I) * (UPQT(K,I)-QT(K))
        edmf_thl(K)=edmf_thl(K)+UPA(K,I)* (UPTHL(K,I)-THL(K))
        edmf_ent(K)=edmf_ent(K)+UPA(K,I)*ENT(K,I)
        edmf_qc(K)=edmf_qc(K)+UPA(K,I)*UPQC(K,I)
        edmf_debug1(K)=edmf_debug1(K) + UPA(K,I) * (UPTHV(K,I)-THV(K))
        edmf_debug2(K)=QT(K) !debug2 is the mean qt to compare with actual output
        ! edmf_debug3(K)=edmf_debug3(K)+UPA(K+1,I)* UPTHL(K+1,I) !debug3 is the actual ud thl
        ! edmf_debug4(K)=THL(K+1) !debug4 is the mean thl to compare with actual output

!        edmf_a(K)=edmf_a(K)+UPA(K+1,I)
!        edmf_w(K)=edmf_w(K)+UPA(K+1,I)*UPW(K+1,I)
!        edmf_qt(K)=edmf_qt(K)+UPA(K+1,I) * (UPQT(K+1,I)-QT(K+1))
!        edmf_thl(K)=edmf_thl(K)+UPA(K+1,I)* (UPTHL(K+1,I)-THL(K+1))
!        edmf_ent(K)=edmf_ent(K)+UPA(K+1,I)*ENT(K+1,I)
!        edmf_qc(K)=edmf_qc(K)+UPA(K+1,I)*UPQC(K+1,I)
!        edmf_debug1(K)=edmf_debug1(K) + UPA(K+1,I) * (UPTHV(K+1,I)-THV(K+1))
!        edmf_debug2(K)=QT(K+1) !debug2 is the mean qt to compare with actual output
!        ! edmf_debug3(K)=edmf_debug3(K)+UPA(K+1,I)* UPTHL(K+1,I) !debug3 is the actual ud thl
!        ! edmf_debug4(K)=THL(K+1) !debug4 is the mean thl to compare with actual output
    ENDDO

    IF (edmf_a(k)>0.) THEN
        edmf_w(k)=edmf_w(k)/edmf_a(k)
        edmf_qt(k)=edmf_qt(k)/edmf_a(k)
        edmf_thl(k)=edmf_thl(k)/edmf_a(k)
        edmf_ent(k)=edmf_ent(k)/edmf_a(k)
        edmf_qc(k)=edmf_qc(k)/edmf_a(k)
        edmf_debug1(k)=edmf_debug1(k)/edmf_a(k)
    ENDIF
ENDDO

  !
  ! computing variables needed for solver
  !

DO k=KTS,KTE+1
    DO I=1,NUP
        s_aw(k)=s_aw(k)+UPA(k,I)*UPW(k,I)
        s_awthl(k)=s_awthl(k)+UPA(k,i)*UPW(k,I)*UPTHL(k,I)
        s_awqt(k)=s_awqt(k)+UPA(k,i)*UPW(k,I)*UPQT(k,I)
        s_awqc(k)=s_awqc(k)+UPA(k,i)*UPW(k,I)*UPQC(k,I)
        s_awu(k)=s_awu(k)+UPA(k,i)*UPW(k,I)*UPU(k,I)
        s_awv(k)=s_awv(k)+UPA(k,i)*UPW(k,I)*UPV(k,I)
    ENDDO
    s_awqv(k) = s_awqt(k)  - s_awqc(k)
ENDDO


!<--- yhc 2021-09-08

!--- compute z at full levels (where T,q are defined)
DO K=KTS,KTE
  ZFULL(K) = 0.5 * (ZW(K)+ZW(K+1))
ENDDO

!--- diagnose detrainment rate
DO K=KTS,KTE-1
   !IF(k > KTOP) exit

  DO I=1,NUP
    IF(I > NUP) exit

    dz=ZW(k+1)-ZW(k)

    mfp1=UPRHO(K+1,I)*UPW(K+1,I)*UPA(K+1,I)
    mf  =UPRHO(K,I)  *UPW(K,I)  *UPA(K,I)
    !if (k.eq.1) then
    !  mf=UPRHO(1,I)*UPW(K,I)*UPA(K,I)
    !else
    !  mf=UPRHO(K-1,I)*UPW(K,I)*UPA(K,I)
    !endif

    if (option_up_area.eq.2) then         ! for lateral entraining plume with detrainment that occurs only at plume top
      if (mf.gt.0 .and. mfp1.le.0.) then
        dz=ZW(K)-ZW(K-1)
        DET(K,I) = 1./dz 
      endif 
    else       ! compute detrainment given entrainment and mass flux
      if (mf.gt.0) then
        DET(K,I) = ENT(K,I) - (mfp1-mf)/mf/dz   
      endif 

      if (option_up_area.eq.3) then       ! hybrid approach
        if (DET(K,I).lt.0.) DET(K,I)=0.   ! in this approach, detrainment should >=0. 
                                          ! but becasue of numeric trancation, DET may<0 (e.g. -1.e-9). In this case, reset DET to 0  
      endif
    endif 

    edmf_det(K)=edmf_det(K)+UPA(K,I)*DET(K,I)
  ENDDO

  IF (edmf_a(k)>0.) THEN
    edmf_det(k)=edmf_det(k)/edmf_a(k)
  ENDIF

ENDDO

!--- write out
if (do_writeout_column_nml) then
  do I=1,NUP
    !where (DET(:,I) .lt. 0) 
    !  print*,'gg,i',I
    !end where
    F1=0  
    DO K=KTS,KTE-1
      if (DET(K,I) .lt. 0) then
        F1=1  
        exit  
      endif 
    ENDDO 
  
    if (F1.eq.1) then
      write(6,3001) 'ZW  = (/',ZW(:)
      write(6,3002) 'UPRHO = (/',UPRHO(:,I)
      write(6,3002) 'UPA = (/',UPA(:,I)
      write(6,3002) 'UPW = (/',UPW(:,I)
      write(6,3002) 'ENT = (/',ENT(:,I)
      write(6,3002) 'DET = (/',DET(:,I)
      write(6,*) '---------------------------'
    endif 
  enddo 
endif

!--- compute variables for coupling with Tiedtke
DO K=KTS,KTE-1
  DO I=1,NUP
  !DO I=1,1
     dz=ZW(k+1)-ZW(k)
   
     !--- mass flux 
     mfp1= UPRHO(K+1,I) * UPA(K+1,I) * UPW(K+1,I)   ! mass flux at k+1/2 level
     mf  = UPRHO(K,I)   * UPA(K,I)   * UPW(K,I)     ! mass fkyx at k-1/2 level

     qcp1 = 0.5 * (qc(K)+qc(K+1))                   ! grid-scale qc at k+1/2 level
     if (K.eq.1) then
       qcp0 = qc(K)                                 ! grid-scale qc at k-1/2 level
     else
       qcp0 = 0.5 * (qc(K)+qc(K-1))                 ! grid-scale qc at k-1/2 level
     endif

     qlp1 = 0.5 * (ql(K)+ql(K+1))                   ! grid-scale ql at k+1/2 level
     if (K.eq.1) then
       qlp0 = ql(K)                                 ! grid-scale ql at k-1/2 level
     else
       qlp0 = 0.5 * (ql(K)+ql(K-1))                 ! grid-scale ql at k-1/2 level
     endif

     qip1 = 0.5 * (qi(K)+qi(K+1))                   ! grid-scale qi at k+1/2 level
     if (K.eq.1) then
       qip0 = qi(K)                                 ! grid-scale qi at k-1/2 level
     else
       qip0 = 0.5 * (qi(K)+qi(K-1))                 ! grid-scale qi at k-1/2 level
     endif

     cldp1 = 0.5 * (cldfra_bl1d(K)+cldfra_bl1d(K+1))          ! grid-scale cloud fraction at k+1/2 level
     if (K.eq.1) then
       cldp0 = cldfra_bl1d(K)                                 ! grid-scale cloud fraction at k-1/2 level
     else
       cldp0 = 0.5 * (cldfra_bl1d(K)+cldfra_bl1d(K-1))        ! grid-scale cloud fraction at k-1/2 level
     endif

     !--- saturated fraction
     IF (UPQC(K+1,I) > 0.) THEN    ! saturated fraction at k+1/2 level
       CCp1=1.
     ELSE
       CCp1=0.
     ENDIF
   
     IF (UPQC(K,I) > 0.) THEN      ! saturated fraction at k-1/2 level
       CCp0=1.
     ELSE
       CCp0=0.
     ENDIF

     !************************
     !************************
     !
     !  eddy-flux/source form
     !
     !************************
     !************************

       !-----------------------------------------------
       !--- cloud condensate, eddy-flux/source form ---
       !-----------------------------------------------
       F1=0. ; F1_ql=0. ; F1_qi=0.
       F2=0. ; F2_ql=0. ; F2_qi=0.
       F3=0. ; F3_ql=0. ; F3_qi=0.

       if (Qx_numerics.eq.1) then   ! upwind approximation
         !--- F1: eddy-flux convergence term for cloud condensate
         !if (K.eq.1) then
           !F1 = -1./rho(k) * ( mfp1*(UPQC(K+1,I)-qcp1) - mf*(UPQC(K,I)-0.) )  / dz     ! qc(K-1)=qc(0)=0.
         !  F1_ql = -1./rho(k) * ( mfp1 * (liquid_frac(k+1)*UPQC(K+1,I)-qlp1) - mf * (liquid_frac(k)*UPQC(K,I)-0.) )  / dz     ! qc(K-1)=qc(0)=0.
         !  F1_qi = -1./rho(k) * ( mfp1 * (ice_frac   (k+1)*UPQC(K+1,I)-qip1) - mf * (ice_frac   (k)*UPQC(K,I)-0.) )  / dz     ! qc(K-1)=qc(0)=0.
         !else
           !F1 = -1./rho(k) * ( mfp1*(UPQC(K+1,I)-qcp1) - mf*(UPQC(K,I)-qcp0)) / dz
           F1_ql = -1./rho(k) * ( mfp1 * (liquid_frac(k+1)*UPQC(K+1,I)-qlp1) - mf * (liquid_frac(k)*UPQC(K,I)-qlp0) )  / dz     ! qc(K-1)=qc(0)=0.
           F1_qi = -1./rho(k) * ( mfp1 * (ice_frac   (k+1)*UPQC(K+1,I)-qip1) - mf * (ice_frac   (k)*UPQC(K,I)-qip0) )  / dz     ! qc(K-1)=qc(0)=0.
         !endif
  
         if (UPW(K,I).gt.0. .and. UPW(K+1,I).gt.0.) then  ! source terms are only present in the plume
           !--- F2: source/vertical advection term for cloud condensate
           !F2 = UPA(K,I) * UPW(K,I) * (UPQC(K+1,I)-UPQC(K,I)) / dz  !+ENT(K,I)*mf*(UPQC(K,I)-qc(K))
           F2_ql = UPA(K,I) * UPW(K,I) * (liquid_frac(K+1)*UPQC(K+1,I)-liquid_frac(K)*UPQC(K,I)) / dz  !+ENT(K,I)*mf*(UPQC(K,I)-qc(K))
           F2_qi = UPA(K,I) * UPW(K,I) * (ice_frac   (K+1)*UPQC(K+1,I)-ice_frac   (K)*UPQC(K,I)) / dz  !+ENT(K,I)*mf*(UPQC(K,I)-qc(K))
    
           !--- F3: source/entrainment term for cloud condensate
           F3_ql = UPA(K,I) * UPW(K,I) * ENT(K,I) * (liquid_frac(K)*UPQC(K,I)-ql(K))
           F3_qi = UPA(K,I) * UPW(K,I) * ENT(K,I) * (ice_frac   (K)*UPQC(K,I)-qi(K))
         endif
       endif  ! end if of Qx_numerics   

         Qql_eddy_i (K,I) = F1_ql !liquid_frac(k)*F1
         Qql_adv_i  (K,I) = F2_ql !liquid_frac(k)*F2
         Qql_ent_i  (K,I) = F3_ql !liquid_frac(k)*F3
    
         Qqi_eddy_i (K,I) = F1_qi !(1.-liquid_frac(k))*F1
         Qqi_adv_i  (K,I) = F2_qi !(1.-liquid_frac(k))*F2
         Qqi_ent_i  (K,I) = F3_qi !(1.-liquid_frac(k))*F3

         !--- sum over the i-th plume
         Qql_eddy(k) = Qql_eddy(k) + Qql_eddy_i (K,I)
         Qql_adv (k) = Qql_adv (k) + Qql_adv_i  (K,I)
         Qql_ent (k) = Qql_ent (k) + Qql_ent_i  (K,I)

         Qqi_eddy(k) = Qqi_eddy(k) + Qqi_eddy_i (K,I)
         Qqi_adv (k) = Qqi_adv (k) + Qqi_adv_i  (K,I)
         Qqi_ent (k) = Qqi_ent (k) + Qqi_ent_i  (K,I)

       !---------------------------------------------
       !--- cloud fraction, eddy-flux/source form ---
       !---------------------------------------------
       F1=0. ; F1_ql=0. ; F1_qi=0.
       F2=0. ; F2_ql=0. ; F2_qi=0.
       F3=0. ; F3_ql=0. ; F3_qi=0.

       if (Qx_numerics.eq.1) then   ! upwind approximation
         !--- F1: eddy-flux convergence term for cloud fraction
         !if (K.eq.1) then
         !  Qa_eddy_i (K,I) = -1./rho(k) * ( mfp1*(CCp1-cldp1) - mf*(CCp0-0.)    ) / dz     ! qc(K-1)=qc(0)=0.
         !else
           Qa_eddy_i (K,I) = -1./rho(k) * ( mfp1*(CCp1-cldp1) - mf*(CCp0-cldp0) ) / dz
         !endif

         if (UPW(K,I).gt.0. .and. UPW(K+1,I).gt.0.) then  ! source terms are only present in the plume
           !--- F2: vertical advection term
           Qa_adv_i (K,I) = UPA(K,I) * UPW(K,I) * (CCp1-CCp0) / dz  !+ENT(K,I)*mf*(UPQC(K,I)-qc(K))
    
           !--- F3: entrainment term
           Qa_ent_i (K,I) = UPA(K,I) * UPW(K,I) * ENT(K,I) * (CCp0-cldfra_bl1d(K))
         endif
       endif  ! end if of Qx_numerics   
  
         Qa_eddy(k) = Qa_eddy(k) + Qa_eddy_i (K,I)
         Qa_adv (k) = Qa_adv (k) + Qa_adv_i  (K,I)
         Qa_ent (k) = Qa_ent (k) + Qa_ent_i  (K,I)

     !************************
     !************************
     !
     !  detrainment/subsidence form
     !
     !************************
     !************************

       !-----------------------------------------------------
       !--- cloud condensate, detrainment/subsidence form ---
       !-----------------------------------------------------
       F1=0. ; F1_ql=0. ; F1_qi=0.
       F2=0. ; F2_ql=0. ; F2_qi=0.
       F3=0. ; F3_ql=0. ; F3_qi=0.

       if (Qx_numerics.eq.1) then   ! upwind approximation
         ! F1: subsidence term for cloud condensate
         !F1 = mfp1/rho(k) * (qc(K+1)-qc(K))/(ZFULL(K+1)-ZFULL(K))  
         F1_ql = mfp1/rho(k) * (ql(K+1)-ql(K))/(ZFULL(K+1)-ZFULL(K))  
         F1_qi = mfp1/rho(k) * (qi(K+1)-qi(K))/(ZFULL(K+1)-ZFULL(K))  

         ! F2: detrainment term for cloud condensate
         !F2 = mf*DET(K,I)/rho(k) * (UPQC(K,I)-qc(K))
         F2_ql = mf*DET(K,I)/rho(k) * (liquid_frac(K)*UPQC(K,I)-ql(K))
         F2_qi = mf*DET(K,I)/rho(k) * (ice_frac   (K)*UPQC(K,I)-qi(K))
       endif  ! end if of Qx_numerics   

         Qql_sub_i (K,I) = F1_ql !liquid_frac(k)*F1
         Qql_det_i (K,I) = F2_ql !liquid_frac(k)*F2

         Qqi_sub_i (K,I) = F1_qi !(1.-liquid_frac(k))*F1
         Qqi_det_i (K,I) = F2_qi !(1.-liquid_frac(k))*F2

       ! sum over all plumes
       Qql_sub(k) = Qql_sub(k) + Qql_sub_i (K,I)
       Qqi_sub(k) = Qqi_sub(k) + Qqi_sub_i (K,I)

       Qql_det(k) = Qql_det(k) + Qql_det_i (K,I)
       Qqi_det(k) = Qqi_det(k) + Qqi_det_i (K,I)

       !---------------------------------------------------
       !--- cloud fraction, detrainment/subsidence form ---
       !---------------------------------------------------
       F1=0. ; F1_ql=0. ; F1_qi=0.
       F2=0. ; F2_ql=0. ; F2_qi=0.
       F3=0. ; F3_ql=0. ; F3_qi=0.

       if (Qx_numerics.eq.1) then   ! upwind approximation
         ! F1: subsidence term for cloud condensate
         Qa_sub_i (K,I) = mfp1/rho(k) * (cldfra_bl1d(K+1)-cldfra_bl1d(K))/(ZFULL(K+1)-ZFULL(K))  

         ! F2: detrainment term for cloud condensate
         Qa_det_i (K,I) = mf*DET(K,I)/rho(k) * (CCp0-cldfra_bl1d(K))
       endif  ! end if of Qx_numerics   

       ! sum over all plumes
       Qa_sub(k) = Qa_sub(k) + Qa_sub_i (K,I)
       Qa_det(k) = Qa_det(k) + Qa_det_i (K,I)

  ENDDO
ENDDO  

!--- return MF cloud tendencies  
if (Qx_MF.eq.1) then   ! eddy-divergence/source form 
  Qql(:) = Qql_eddy(:) + Qql_adv(:) +Qql_ent(:) 
  Qqi(:) = Qqi_eddy(:) + Qqi_adv(:) +Qqi_ent(:) 
  Qa (:) = Qa_eddy (:) + Qa_adv (:) +Qa_ent (:) 

elseif (Qx_MF.eq.2) then ! detrainment/subsidence form
  Qql(:) = Qql_det(:) + Qql_sub(:)
  Qqi(:) = Qqi_det(:) + Qqi_sub(:)
  Qa (:) = Qa_det (:) + Qa_sub (:)

else
  print*,'ERROR: Qx_MF must be 1 or 2'
endif

  !print*,'Qql',Qql
  !print*,'Qqi',Qqi
  !print*,'Qa',Qa

!--- obtain (1) mass flux, (2) fraction, and (3) plume-averaged specific humidiry for moist and dry updrafts

  !--- initial variables
  a_moist_half  = 0.    
  a_moist_full  = 0.
  mf_moist_half = 0.    
  mf_moist_full = 0.
  qv_moist_half = 0.    
  qv_moist_full = 0.

  a_dry_half    = 0.    
  a_dry_full    = 0.    
  mf_dry_half   = 0.
  mf_dry_full   = 0.
  qv_dry_half   = 0.
  qv_dry_full   = 0.

  mf_all_half   = 0.
  mf_all_full   = 0.

!--- obtain (1) mass flux, (2) fraction, and (3) plume-averaged specific humidiry for moist and dry updrafts
DO K=KTS,KTE-1
  !IF(k > KTOP) exit

  DO I=1,NUP
    IF(I > NUP) exit

    if (UPQC(K,I) > 0.) then   ! sum of individual moist updrafts
      a_moist_half  (K) = a_moist_half  (K) + UPA(K,I) 
      mf_moist_half (K) = mf_moist_half (K) + UPA(K,I)*UPRHO(K,I)*UPW(K,I)
      qv_moist_half (K) = qv_moist_half (K) + UPA(K,I)*(UPQT(K,I)-UPQC(K,I))  

    else                       ! sum of individual dry updrafts
      a_dry_half    (K) = a_dry_half  (K) + UPA(K,I) 
      mf_dry_half   (K) = mf_dry_half (K) + UPA(K,I)*UPRHO(K,I)*UPW(K,I)
      qv_dry_half   (K) = qv_dry_half (K) + UPA(K,I)*(UPQT(K,I)-UPQC(K,I))        
    endif
  ENDDO  ! end loop of i

  if (a_moist_half(K) > 0.) then
    qv_moist_half (K) = qv_moist_half (K) / a_moist_half (K)
  end if

  if (a_dry_half  (K) > 0.) then
    qv_dry_half   (K) = qv_dry_half   (K) / a_dry_half (K)
  end if
ENDDO    ! end loop of k

!--- interpolate a, mf, qv from half levels to full levels
!      interpolate from the updraft base. Right above the updraft top (KTOP, mf>0), using upwind approximation

  !--- updraft top
  !a_moist_full  (KTOP) = a_moist_half  (KTOP)
  !mf_moist_full (KTOP) = mf_moist_half (KTOP)
  !qv_moist_full (KTOP) = qv_moist_half (KTOP)

  !a_dry_full    (KTOP) = a_dry_half    (KTOP)
  !mf_dry_full   (KTOP) = mf_dry_half   (KTOP)
  !qv_dry_full   (KTOP) = qv_dry_half   (KTOP)

DO K=KTS,KTE-1
  if (mf_moist_half (K) > 0.) then
    a_moist_full  (K) = 0.5 * (a_moist_half  (K)+a_moist_half  (K+1))
    mf_moist_full (K) = 0.5 * (mf_moist_half (K)+mf_moist_half (K+1))
    qv_moist_full (K) = 0.5 * (qv_moist_half (K)+qv_moist_half (K+1))
  endif

  if (mf_dry_half (K) > 0.) then
    a_dry_full  (K) = 0.5 * (a_dry_half  (K)+a_dry_half  (K+1))
    mf_dry_full (K) = 0.5 * (mf_dry_half (K)+mf_dry_half (K+1))
    qv_dry_full (K) = 0.5 * (qv_dry_half (K)+qv_dry_half (K+1))
  endif
ENDDO


!print*,'a_dry_full',a_dry_full
!print*,'a_dry_half',a_dry_half
!print*,'mf_dry_full',mf_dry_full
!print*,'mf_dry_half',mf_dry_half
!print*,'qv_dry_full',qv_dry_full
!print*,'qv_dry_half',qv_dry_half
!
!print*,'a_moist_full',a_moist_full
!print*,'a_moist_half',a_moist_half
!print*,'mf_moist_full',mf_moist_full
!print*,'mf_moist_half',mf_moist_half
!print*,'qv_moist_full',qv_moist_full
!print*,'qv_moist_half',qv_moist_half

! moist+dry mass flux
  mf_all_half (:) = mf_moist_half (:) + mf_dry_half (:)
  mf_all_full (:) = mf_moist_full (:) + mf_dry_full (:)

!--- check part
if (do_check_ent_det) then
  DO I=1,NUP
    DO K=KTS,KTE-1
      IF(I > NUP) exit
        dz=ZW(k+1)-ZW(k)
        mfp1=UPRHO(K+1,I)*UPW(K+1,I)*UPA(K+1,I)
        mf  =UPRHO(K,I)  *UPW(K,I)  *UPA(K,I)    
  
        F1=0.
        F2=0.
        F3=0.
  
        if (mf.gt.0.) then
          F1 = 1./mf * (mfp1-mf)/dz
          F2 = ENT(K,I) - DET(K,I)  
          if (F1.gt.0.) F3 = (F1-F2)/F1 * 100.
  
          !print*,'I,K, ENT, DET, 1/m*dm/dz, ENT-DET, %diff',I,K, ENT(K,I),DET(K,I),F1,F2,F3
          print*,'I,K, 1/m*dm/dz, ENT, DET, ENT-DET, %diff',I,K,F1,ENT(K,I),DET(K,I),F2,F3
        endif
    ENDDO  ! end loop of K
  ENDDO    ! end loop of I
endif      ! end if of do_check_ent_det


!--- compute the frequency of negative DET
num_updraft=0.
num_DET=0.
num_nDET_zENT=0.
num_nDET_pENT=0.

DO K=KTS,KTE-1
  DO I=1,NUP
    mf = UPRHO(K,I)*UPW(K,I)*UPA(K,I)

    if (mf.gt.0) then
      num_updraft(K) = num_updraft(K)+1.  ! count # of updrafts at each level

      if (DET(K,I).ne.0.) then      ! count # of detrainment at each level
        num_DET (K) = num_DET(K)+1.
      endif

      if (DET(K,I).lt.0. .and. ENT(K,I).eq.0.) then  ! count # of negative detrainment with zero entrainment
        num_nDET_zENT (K) = num_nDET_zENT (K)+1. 
      endif

      if (DET(K,I).lt.0. .and. ENT(K,I).gt.0.) then  ! count # of negative detrainment with positive entrainment
        num_nDET_pENT (K) = num_nDET_pENT (K)+1. 
      endif
    endif  ! end if of mf>0
  ENDDO

  if (num_DET (K).gt.0.) then  ! compute frequency
    num_nDET_zENT (K) = num_nDET_zENT (K) / num_DET (K)
    num_nDET_pENT (K) = num_nDET_pENT (K) / num_DET (K)
  endif
ENDDO
!---> yhc 2021-09-08


!
! make sure that the mass-flux does not exceed CFL criteria, and if it does scale it back
! (see Beljaars et al., 2018)
!

maxS=0.
DO k=KTS,KTE
  maxS=max(maxS,0.5*(s_aw(k)+s_aw(k+1))*dt/(zw(k+1)-zw(k)))
ENDDO


IF (maxS .gt. upstab) THEN
! if stability exceeded, scale the fluxes
! (note: the updraft properties are not modified)
  stabF=upstab/maxS
!  print *,'maxS',maxS
  DO k=kts,KTE
        s_aw(k)=s_aw(k)*stabF
        s_awthl(k)=s_awthl(k)*stabF
        s_awqt(k)=s_awqt(k)*stabF
        s_awqc(k)=s_awqc(k)*stabF
        s_awu(k)=s_awu(k)*stabF
        s_awv(k)=s_awv(k)*stabF
        s_awqv(k)=s_awqv(k)*stabF
  ENDDO
ENDIF

3001 format (A35,2X,34(F10.3,2X,','))
3002 format (A35,2X,34(E12.4,2X,','))

END SUBROUTINE edmf_JPL

! ===================================================================
! ===================================================================
! This is the downdraft mass flux scheme - analogus to edmf_JPL but 
! flipped updraft to downdraft. This scheme is currently only tested
! for Stratocumulus cloud conditions. For a detailed desctiption of the
! model, see paper.

! SUBROUTINE DDMF_JPL(kts,kte,dt,zw,p,                 &
!               &u,v,th,thl,thv,tk,qt,qv,qc,           &
!               &exner,                                &
!               &ust,wthl,wqt,pblh,kpbl,               &
!               &edmf_a_dd,edmf_w_dd, edmf_qt_dd,      &
!               &edmf_thl_dd,edmf_ent_dd,edmf_qc_dd,   &
!               &edmf_debug3,edmf_debug4,              &
!               &sd_aw,sd_awthl,sd_awqt,               &
!               &sd_awqv,sd_awqc,sd_awu,sd_awv,        &
!               &sd_awqke,                             &
!               &qc_bl1d,cldfra_bl1d,                  &
!               &rthraten                              &
!               )
! 
!         INTEGER, INTENT(IN) :: KTS,KTE,KPBL
!         REAL,DIMENSION(KTS:KTE), INTENT(IN) :: U,V,TH,THL,TK,QT,QV,QC,THV,P,exner,rthraten
!         ! zw .. heights of the downdraft levels (edges of boxes)
!         REAL,DIMENSION(KTS:KTE+1), INTENT(IN) :: ZW
!         REAL, INTENT(IN) :: DT,UST,WTHL,WQT,PBLH
! 
!     ! outputs - downdraft properties
!         REAL,DIMENSION(KTS:KTE), INTENT(OUT) :: edmf_a_dd,edmf_w_dd,        &
!                       & edmf_qt_dd,edmf_thl_dd, edmf_ent_dd,edmf_qc_dd,edmf_debug3,edmf_debug4
! 
!     ! outputs - variables needed for solver (sd_aw - sum ai*wi, sd_awphi - sum ai*wi*phii)
!         REAL,DIMENSION(KTS:KTE+1) :: sd_aw, sd_awthl, sd_awqt, sd_awu, sd_awv, &
!                                      sd_awqc, sd_awqv, sd_awqke, sd_aw2
! 
!         REAL,DIMENSION(KTS:KTE), INTENT(IN) :: qc_bl1d, cldfra_bl1d
! 
!         INTEGER, PARAMETER :: NDOWN=10, debug_mf=0 !fixing number of plumes to 10
!     ! draw downdraft starting height randomly between cloud base and cloud top
!         INTEGER, DIMENSION(1:NDOWN) :: DD_initK
!         REAL   , DIMENSION(1:NDOWN) :: randNum
!     ! downdraft properties
!         REAL,DIMENSION(KTS:KTE+1,1:NDOWN) :: DOWNW,DOWNTHL,DOWNQT,DOWNQC,DOWNA,DOWNU,DOWNV,DOWNTHV
! 
!     ! entrainment variables
!         REAl,DIMENSION(KTS+1:KTE+1,1:NDOWN) :: ENT,ENTf
!         INTEGER,DIMENSION(KTS+1:KTE+1,1:NDOWN) :: ENTi
! 
!     ! internal variables
!         INTEGER :: K,I, kminrad, qlTop, p700_ind, qlBase
!         REAL :: wthv,wstar,qstar,thstar,sigmaW,sigmaQT,sigmaTH,z0, &
!             pwmin,pwmax,wmin,wmax,wlv,wtv,went
!         REAL :: B,QTn,THLn,THVn,QCn,Un,Vn,Wn2,EntEXP,EntW, Beta_dm, EntExp_M
!         REAL :: jump_thetav, jump_qt, jump_thetal, refTHL, refTHV, refQT
!     ! DD specific internal variables
!         REAL :: minrad,zminrad, radflux, F0, wst_rad, wst_dd, deltaZ
!         logical :: cloudflg
! 
!     ! VARIABLES FOR CHABOUREAU-BECHTOLD CLOUD FRACTION
!         ! REAL,DIMENSION(KTS:KTE), INTENT(INOUT) :: vt, vq, sgm
!         REAL :: sigq,xl,tlk,qsat_tl,rsl,cpm,a,qmq,mf_cf,diffqt,&
!                Fng,qww,alpha,beta,bb,f,pt,t,q2p,b9,satvp,rhgrid
! 
!         INTEGER, DIMENSION(2) :: seedmf
! 
!     ! w parameters
!         REAL,PARAMETER :: &
!             &Wa=1., &
!             &Wb=1.5,&
!             &Z00=100.
!     ! entrainment parameters
!         REAL,PARAMETER :: &
!         & L0=80,&
!         & ENT0=0.2
! 
!         pwmin=-3. ! drawing from the neagtive tail -3sigma to -1sigma
!         pwmax=-1.
! 
!     ! initialize downdraft properties
!         DOWNW=0.
!         DOWNTHL=0.
!         DOWNTHV=0.
!         DOWNQT=0.
!         DOWNA=0.
!         DOWNU=0.
!         DOWNV=0.
!         DOWNQC=0.
!         ENT=0.
!         DD_initK=0
! 
!         edmf_a_dd  =0.
!         edmf_w_dd  =0.
!         edmf_qt_dd =0.
!         edmf_thl_dd=0.
!         edmf_ent_dd=0.
!         edmf_qc_dd =0.
!         edmf_debug3=0.
!         edmf_debug4=0.
! 
!         sd_aw=0.
!         sd_awthl=0.
!         sd_awqt=0.
!         sd_awqv=0.
!         sd_awqc=0.
!         sd_awu=0.
!         sd_awv=0.
!         sd_awqke=0.
! 
!     ! FIRST, CHECK FOR STRATOCUMULUS-TOPPED BOUNDARY LAYERS 
!         cloudflg=.false.
!         minrad=100.
!         kminrad=kpbl
!         zminrad=PBLH
!         qlTop = 1 !initialize at 0
!         qlBase = 1
!         wthv=wthl+svp1*wqt
!         do k = MAX(1,kpbl-2),kpbl+3
!             if(qc(k).gt. 1.e-6 .AND. cldfra_bl1D(k).gt.0.5) then
!                cloudflg=.true. !found Sc cloud
!                qlTop = k ! index for Sc cloud top
!             endif
!         enddo
! 
!         do k = qlTop, kts, -1
!             if(qc(k).gt. 1E-6) then
!                 qlBase = k ! index for Sc cloud base
!             endif
!         enddo
!         qlBase = (qlTop+qlBase)/2 ! changed base to half way through the cloud
! 
! !        call init_random_seed_1()
!         call RANDOM_NUMBER(randNum)
!         do i=1,NDOWN
!             ! downdraft starts somewhere between cloud base to cloud top
!             ! the probability is equally distributed
!             DD_initK(i) = qlTop ! nint(randNum(i)*REAL(qlTop-qlBase)) + qlBase
!         enddo
! 
!         ! LOOP RADFLUX
!         F0 = 0.
!         do k = 1, qlTop ! Snippet from YSU, YSU loops until qlTop - 1
!            radflux = rthraten(k) * exner(k) ! Converts theta/s to temperature/s
!            radflux = radflux * cp / g * ( p(k) - p(k+1) ) ! Converts temperature/s to W/m^2
!            if ( radflux < 0.0 ) F0 = abs(radflux) + F0
!         enddo
!         F0 = max(F0, 1.0)
!         !found Sc cloud and cloud not at surface, trigger downdraft
!         if (cloudflg) then 
! 
!             !get entrainent coefficient
!             do i=1,NDOWN
!                 do k=kts+1,kte
!                   ENTf(k,i)=(ZW(k+1)-ZW(k))/L0
!                  enddo
!             enddo
! 
! 
! ! create seed for Poisson
! ! create seed from last digits of temperature   
! seedmf(1) = 1000000 * ( 100*thl(1) - INT(100*thl(1)))
! seedmf(2) = 1000000 * ( 100*thl(2) - INT(100*thl(2))) 
! 
! 
!             ! get Poisson P(dz/L0)
!             call Poisson(1,NDOWN,kts+1,kte,ENTf,ENTi,seedmf)
! 
! 
!             ! entrainent: Ent=Ent0/dz*P(dz/L0)
!             do i=1,NDOWN
!                 do k=kts+1,kte
!                   ENT(k,i)=real(ENTi(k,i))*Ent0/(ZW(k+1)-ZW(k))
!                 enddo
!             enddo            
! 
!             !!![EW: INVJUMP] find 700mb height then subtract trpospheric lapse rate!!!
!             p700_ind = MINLOC(ABS(p-70000),1)!p1D is 70000
!             jump_thetav = thv(p700_ind) - thv(1) - (thv(p700_ind)-thv(qlTop+3))/(ZW(p700_ind)-ZW(qlTop+3))*(ZW(p700_ind)-ZW(qlTop))
!             jump_qt = qc(p700_ind) + qv(p700_ind) - qc(1) - qv(1)
!             jump_thetal = thl(p700_ind) - thl(1) - (thl(p700_ind)-thl(qlTop+3))/(ZW(p700_ind)-ZW(qlTop+3))*(ZW(p700_ind)-ZW(qlTop))
!             
!             refTHL = thl(qlTop) !sum(thl(1:qlTop)) / (qlTop) ! avg over BL for now or just at qlTop
!             refTHV = thv(qlTop) !sum(thv(1:qlTop)) / (qlTop)
!             refQT = qt(qlTop) !sum(qt(1:qlTop)) / (qlTop)            
!             wst_rad = ( g * zw(qlTop) * F0 / (refTHL * rhoair0 * cp) ) ** (0.333) ! wstar_rad, following Lock and MacVean (1999a)
!             wst_rad = max(wst_rad, 0.1)
!             wstar = max(0.,(g/thv(1)*wthv*pblh)**(1./3.))
!             went = thv(1) / ( g * jump_thetav * zw(qlTop) ) * (0.15 * (wstar**3 + 5*ust**3) + 0.35 * wst_rad**3 )
!             qstar = abs(went*jump_qt/wst_rad)
!             thstar = F0/rhoair0/cp/wst_rad - went*jump_thetav/wst_rad
!             wst_dd = (0.15 * (wstar**3 + 5*ust**3) + 0.35 * wst_rad**3 ) ** (0.333) !wstar_dd - mix rad + surface wst
!             sigmaW = 0.2*wst_dd! 0.8*wst_dd!wst_rad !tuning parameter ! 0.5 was good
!             sigmaQT = 40  * qstar ! 50 was good
!             sigmaTH = 1.0 * thstar! 0.5 was good
! 
!             wmin=sigmaW*pwmin
!             wmax=sigmaW*pwmax
! 
!             do I=1,NDOWN !downdraft now starts at different height
!                 wlv=wmin+(wmax-wmin)/NDOWN*(i-1)
!                 wtv=wmin+(wmax-wmin)/NDOWN*i
! 
!                 DOWNW(DD_initK(I),I)=0.5*(wlv+wtv)
!                 DOWNA(DD_initK(I),I)=0.5*ERF(wtv/(sqrt(2.)*sigmaW))-0.5*ERF(wlv/(sqrt(2.)*sigmaW))
! 
!                 DOWNU(DD_initK(I),I)=U(DD_initK(I))
!                 DOWNV(DD_initK(I),I)=V(DD_initK(I))
! 
!                 refTHL = 0.5 * (thl(DD_initK(I))+thl(DD_initK(I)-1)) !reference now depends on where dd starts
!                 refTHV = 0.5 * (thv(DD_initK(I))+thv(DD_initK(I)-1))
!                 refQT  = 0.5 * (qt(DD_initK(I))+qt(DD_initK(I)-1))
! 
!                 DOWNQC(DD_initK(I),I) = 0.
!                 DOWNQT(DD_initK(I),I) = refQT  + 0.5  *DOWNW(DD_initK(I),I)*sigmaQT/sigmaW
!                 DOWNTHV(DD_initK(I),I)= refTHV + 0.01 *DOWNW(DD_initK(I),I)*sigmaTH/sigmaW ! close to no difference
! 
!                 call condensation_edmf_r(DOWNQT(DD_initK(I),I),DOWNTHL(DD_initK(I),I),(P(DD_initK(I))+P(DD_initK(I)-1))/2.,ZW(DD_initK(I)),DOWNTHV(DD_initK(I),I),DOWNQC(DD_initK(I),I))
! 
!                 ! DOWNTHL(DD_initK(I),I)=DOWNTHV(DD_initK(I),I)/(1.+svp1*DOWNQT(DD_initK(I),I))
!                 ! print *, "Plume ", I, " DOWNQT = ", DOWNQT(DD_initK(I),I), " refQT = ", refQT, &
!                 !          " DOWNTHV = ", DOWNTHV(DD_initK(I),I), " refTHV = ", refTHV, &
!                 !          " DOWNTHL = ", DOWNTHL(DD_initK(I),I), " refTHL = ", refTHL
!             enddo
!              ! print *, "RefTHV = ", refTHV, " sigmaQT = ", sigmaQT, " sigmaTH = ", sigmaTH, " sigmaW = ", sigmaW, " total = ", -0.2*DOWNW(qlTop,1)*sigmaTH/sigmaW
!             print*, "[EWDD] wst = ", wstar, " wst_rad = ", wst_rad, " wst_dd = ", wst_dd, " went = ", went
!             ! do integration downdraft
!             DO I=1,NDOWN
!                 DO k=DD_initK(I)-1,KTS+1,-1 !starting one point below cloud top
!                     deltaZ = ZW(k+1)-ZW(k)
!                     EntExp=exp(-ENT(K,I)*deltaZ)
!                     EntExp_M=exp(-ENT(K,I)/3.*deltaZ)
! 
!                     QTn=DOWNQT(K+1,I)+(QT(K)-DOWNQT(K+1,I))*(1.-EntExp)
!                     THLn=DOWNTHL(K+1,I)+(THL(K)-DOWNTHL(K+1,I))*(1.-EntExp)
!                     Un=DOWNU(K+1,I)+(U(K)-DOWNU(K+1,I))*(1.-EntExp_M)
!                     Vn=DOWNV(K+1,I)+(V(K)-DOWNV(K+1,I))*(1.-EntExp_M)
!                     ! get thvn,qcn
!                     call condensation_edmf(QTn,THLn,(P(K)+P(K-1))/2.,ZW(k),THVn,QCn)
!                     B=g*(0.5*(THVn+DOWNTHV(k+1,I))/THV(k)-1.)
! 
!                     Beta_dm = 2*Wb*ENT(K,I) + 0.5/(ZW(k)-deltaZ) * max(1. - exp((ZW(k) -deltaZ)/Z00 - 1. ) , 0.)
!                     EntW=exp(-Beta_dm * deltaZ)
!                     if (Beta_dm >0) then
!                         Wn2=DOWNW(K+1,I)**2*EntW - Wa*B/Beta_dm * (1. - EntW)
!                     else
!                         Wn2=DOWNW(K+1,I)**2 - 2.*Wa*B*deltaZ
!                     end if
!                      print *, "Plume number = ", I, " k = ", k, " z = ", ZW(k), " entw = ", ENT(K,I)
!                      print *, "downthv = ", THVn, " thv = ", thv(k)
!                      print *, "downthl = ", THLn, " thl = ", thl(k)
!                      print *, "downqt  = ", QTn , " qt  = ", qt(k)
!                      print *, "Beta = ", Beta_dm, " Wn2 = ", Wn2, " Bouy = ", B
! 
!                     IF (Wn2 .gt. 0.) THEN !terminate when velocity is too small
!                         DOWNW(K,I)=-sqrt(Wn2)
!                         DOWNTHV(K,I)=THVn
!                         DOWNTHL(K,I)=THLn
!                         DOWNQT(K,I)=QTn
!                         DOWNQC(K,I)=QCn
!                         DOWNU(K,I)=Un
!                         DOWNV(K,I)=Vn
!                         DOWNA(K,I)=DOWNA(K+1,I)
!                     ELSE
!                           exit
!                     ENDIF
!                 ENDDO
!             ENDDO
!         endif ! end cloud flag
! 
!         DOWNW(1,:) = 0. !make sure downdraft does not go to the surface
!         DOWNA(1,:) = 0.
!         ! do I=1,NDOWN
!         !     DO k=qlTop,1,-1
!         !         if (DOWNA(k,I)>0.) then
!         !             print *, "[EWDD] Plume number = ", I, " k = ", k, " Downdraft w = ", DOWNW(k,I), " area = ", DOWNA(k,I)
!         !             print *, "[EWDD] qt = ", DOWNQT(k,I), " thv = ", DOWNTHV(k,I), " thl = ", DOWNTHL(k,I)
!         !             print *, "[EWDD] refqt = ", qt(k), " refthv = ", thv(k), " refthl = ", thl(k)
!         !         endif
!         !     ENDDO
!         ! enddo
!         ! print*, "[EWDD] wst = ", wstar, " wst_rad = ", wst_rad, " wst_dd = ", wst_dd, " went = ", went
! 
!         ! Combine both moist and dry plume, write as one averaged plume
!         ! Even though downdraft starts at different height, average all up to qlTop
!         DO k=qlTop,KTS,-1
!             DO I=1,NDOWN
!                 IF(I > NDOWN) exit
!                 edmf_a_dd(K)=edmf_a_dd(K)+DOWNA(K-1,I)
!                 edmf_w_dd(K)=edmf_w_dd(K)+DOWNA(K-1,I)*DOWNW(K-1,I)
!                 edmf_qt_dd(K)=edmf_qt_dd(K)+DOWNA(K-1,I)* (DOWNQT(K-1,I)-QT(K-1))
!                 edmf_thl_dd(K)=edmf_thl_dd(K)+DOWNA(K-1,I)* (DOWNTHL(K-1,I)-THL(K-1))
!                 edmf_ent_dd(K)=edmf_ent_dd(K)+DOWNA(K-1,I)*ENT(K-1,I)
!                 edmf_qc_dd(K)=edmf_qc_dd(K)+DOWNA(K-1,I)*DOWNQC(K-1,I)
!                 edmf_debug3(K)=edmf_debug3(K)+DOWNA(K-1,I)* (DOWNTHV(K-1,I)-THV(K-1))
!                 edmf_debug4(K)=edmf_debug4(K)+THL(k)
!             ENDDO
! 
!             IF (edmf_a_dd(k)>0.) THEN
!                 edmf_w_dd(k)=edmf_w_dd(k)/edmf_a_dd(k)
!                 edmf_qt_dd(k)=edmf_qt_dd(k)/edmf_a_dd(k)
!                 edmf_thl_dd(k)=edmf_thl_dd(k)/edmf_a_dd(k)
!                 edmf_ent_dd(k)=edmf_ent_dd(k)/edmf_a_dd(k)
!                 edmf_qc_dd(k)=edmf_qc_dd(k)/edmf_a_dd(k)
!                 edmf_debug3(k)=edmf_debug3(k)/edmf_a_dd(k)
!             ENDIF
!         ENDDO        
! 
!           !
!           ! computing variables needed for solver
!           !
! 
!         DO k=KTS,qlTop
!             DO I=1,NDOWN
!                 sd_aw(k)=sd_aw(k)+DOWNA(k,I)*DOWNW(k,I)
!                 sd_awthl(k)=sd_awthl(k)+DOWNA(k,i)*DOWNW(k,I)*DOWNTHL(k,I)
!                 sd_awqt(k)=sd_awqt(k)+DOWNA(k,i)*DOWNW(k,I)*DOWNQT(k,I)
!                 sd_awqc(k)=sd_awqc(k)+DOWNA(k,i)*DOWNW(k,I)*DOWNQC(k,I)
!                 sd_awu(k)=sd_awu(k)+DOWNA(k,i)*DOWNW(k,I)*DOWNU(k,I)
!                 sd_awv(k)=sd_awv(k)+DOWNA(k,i)*DOWNW(k,I)*DOWNV(k,I)
!             ENDDO
!             sd_awqv(k) = sd_awqt(k)  - sd_awqc(k)
!         ENDDO
! 
!         ! ! This last piece is from STEM_MF
!         ! DO K=KTS,qlTop
!         !     IF(edmf_qc_dd(k)>0.0)THEN
!         !         satvp = 3.80*exp(17.27*(th(k)-273.)/ &
!         !            (th(k)-36.))/(.01*p(k))
!         !         rhgrid = max(.01,MIN( 1., qv(k) /satvp))
! 
!         !         !COMPUTE CLDFRA & QC_BL FROM MASS-FLUX SCHEME and recompute vt & vq
! 
!         !         xl = xl_blend(tk(k))                ! obtain blended heat capacity 
!         !         tlk = thl(k)*(p(k)/p1000mb)**rcp    ! recover liquid temp (tl) from thl
!         !         qsat_tl = qsat_blend(tlk,p(k))      ! get saturation water vapor mixing ratio
!         !                                         !   at tl and p
!         !         rsl = xl*qsat_tl / (r_v*tlk**2)     ! slope of C-C curve at t = tl
!         !                                         ! CB02, Eqn. 4
!         !         cpm = cp + qt(k)*cpv                ! CB02, sec. 2, para. 1
!         !         a   = 1./(1. + xl*rsl/cpm)          ! CB02 variable "a"
!         !         b9  = a*rsl                         ! CB02 variable "b" 
! 
!         !         q2p  = xlvcp/exner(k)
!         !         pt = thl(k) +q2p*edmf_qc_dd(k) ! potential temp
!         !         bb = b9*tk(k)/pt ! bb is "b9" in BCMT95.  Their "b9" differs from
!         !                    ! "b9" in CB02 by a factor
!         !                    ! of T/theta.  Strictly, b9 above is formulated in
!         !                    ! terms of sat. mixing ratio, but bb in BCMT95 is
!         !                    ! cast in terms of sat. specific humidity.  The
!         !                    ! conversion is neglected here.
!         !         qww   = 1.+0.61*qt(k)
!         !         alpha = 0.61*pt
!         !         t     = th(k)*exner(k)
!         !         beta  = pt*xl/(t*cp) - 1.61*pt
!         !         !Buoyancy flux terms have been moved to the end of this section...
! 
!         !         !Now calculate convective component of the cloud fraction:
!         !         if (a > 0.0) then
!         !             f = MIN(1.0/a, 4.0)              ! f is vertical profile scaling function (CB2005)
!         !         else
!         !             f = 1.0
!         !         endif
!         !         sigq = 9.E-3 * edmf_a_dd(k) * edmf_w_dd(k) * f ! convective component of sigma (CB2005)
!         !         !sigq = MAX(sigq, 1.0E-4)         
!         !         sigq = SQRT(sigq**2 + sgm(k)**2)    ! combined conv + stratus components
! 
!         !         qmq = a * (qt(k) - qsat_tl)         ! saturation deficit/excess;
!         !                                         !   the numerator of Q1
!         !         mf_cf = min(max(0.5 + 0.36 * atan(1.55*(qmq/sigq)),0.02),0.6)
!         !         IF ( debug_code ) THEN
!         !             print*,"In MYNN, EDMF JPL"
!         !             print*,"  CB: qt=",qt(k)," qsat=",qsat_tl," satdef=",qt(k) - qsat_tl
!         !             print*,"  CB: sigq=",sigq," qmq=",qmq," tlk=",tlk
!         !             print*,"  CB: mf_cf=",mf_cf," cldfra_bl=",cldfra_bl1d(k)," edmf_a_dd=",edmf_a_dd(k)
!         !         ENDIF
! 
!         !         IF (rhgrid >= .93) THEN
!         !             !IN high RH, defer to stratus component if > convective component
!         !             cldfra_bl1d(k) = MAX(mf_cf, cldfra_bl1d(k))
!         !             IF (cldfra_bl1d(k) > edmf_a_dd(k)) THEN
!         !                 qc_bl1d(k) = edmf_qc_dd(k)*edmf_a_dd(k)/cldfra_bl1d(k)
!         !             ELSE
!         !                 cldfra_bl1d(k)=edmf_a_dd(k)
!         !                 qc_bl1d(k) = edmf_qc_dd(k)
!         !             ENDIF
!         !         ELSE
!         !             IF (mf_cf > edmf_a_dd(k)) THEN
!         !                 cldfra_bl1d(k) = mf_cf
!         !                 qc_bl1d(k) = edmf_qc_dd(k)*edmf_a_dd(k)/mf_cf
!         !             ELSE
!         !                 cldfra_bl1d(k)=edmf_a_dd(k)
!         !                 qc_bl1d(k) = edmf_qc_dd(k)
!         !             ENDIF
!         !         ENDIF
!         !         !Now recalculate the terms for the buoyancy flux for mass-flux clouds:
!         !         !See mym_condensation for details on these formulations.  The
!         !         !cloud-fraction bounding was added to improve cloud retention,
!         !         !following RAP and HRRR testing.
!         !         Fng = 2.05 ! the non-Gaussian transport factor (assumed constant)
!         !         vt(k) = qww   - MIN(0.20,cldfra_bl1D(k))*beta*bb*Fng - 1.
!         !         vq(k) = alpha + MIN(0.20,cldfra_bl1D(k))*beta*a*Fng  - tv0
!         !     ENDIF
!         ! ENDDO
! END SUBROUTINE DDMF_JPL


!################################
!################################
!
!  GFDL model interface, yhc
!    Mellor-Yamada
!
!################################
!################################

subroutine edmf_mynn_driver ( &
              is, ie, js, je, npz, Time_next, dt, lon, lat, frac_land, area, u_star,  &
              b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, &
              rdt_mynn_ed_am4, &
              do_edmf_mynn_diagnostic, do_return_edmf_mynn_diff_only, do_edmf_mynn_in_physics, do_tracers_selective, &
              option_edmf2ls_mp, qadt_edmf, qldt_edmf, qidt_edmf, dqa_edmf,  dql_edmf, dqi_edmf, diff_t_edmf, diff_m_edmf, kpbl_edmf, &
              edmf_mc_full, edmf_mc_half, edmf_moist_area, edmf_dry_area, edmf_moist_humidity, edmf_dry_humidity, &
              pbltop, udt, vdt, tdt, rdt, rdiag)

!---------------------------------------------------------------------
! Arguments (Intent in)  
!
!   Descriptions of the input arguments are in subroutine edmf_alloc
!
!     do_edmf_mynn_diagnostic       = .true.,  edmf_mynn is only for diagnostic purpose, i.e. no changes to the model
!                                              so the model should be exactly the same as the standard one
!                                     .false., edmf_mynn will affect the model simulation
!
!     do_return_edmf_mynn_diff_only = .true.,  return edmf_mynn diffusion coefficients 
!                                              and let vert_diff do the diffusion rather than in edmf_mynn
!                                     .false., tendencies are updated in edmf_mynn
!
!     do_edmf_mynn_in_physics       = "down",  edmf_mynn_driver is called in physics_driver_down
!                                     "up"  ,  edmf_mynn_driver is called in physics_driver_up      
!
!     rdt_mynn_ed_am4
!       multiple tracer tendencies [ unit / unit / sec ] only for eddy-diffusion
!       This is used to recover dry variable tendencies when edmf_type=2
!---------------------------------------------------------------------
  integer, intent(in)                   :: is, ie, js, je, npz
  type(time_type), intent(in)           :: Time_next
  real,    intent(in)                   :: dt
  real,    intent(in), dimension(:,:)   :: &
   lon, lat, &  ! longitude and latitude in radians
   frac_land, area, u_star, b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux

  real, intent(in), dimension(:,:,:,:) :: &
    rdt_mynn_ed_am4

  type(physics_input_block_type)        :: Physics_input_block
  logical, intent(in)                   :: do_edmf_mynn_diagnostic
  logical, intent(in)                   :: do_return_edmf_mynn_diff_only
  character*5, intent(in)               :: do_edmf_mynn_in_physics
  integer, intent(in)                   :: do_tracers_selective

!---------------------------------------------------------------------
! Arguments (Intent inout)  
!
!      Physics_tendency_block derived type variable containing:
!          1) udt           zonal wind tendency [ m / s**2 ]
!          2) vdt           meridional wind tendency [ m / s**2 ]
!          3) tdt           temperature tendency [ deg k / sec ]
!          4) rdt           multiple tracer tendencies 
!                            (index 1 = specific humidity) 
!                            [ unit / unit / sec ]
!          5) rdiag          multiple 3d diagnostic tracer fields 
!                            [ unit / unit ]
!
!      pbltop - PBL height (m)
!---------------------------------------------------------------------
  real, intent(inout), dimension(:,:,:) :: &
    udt, vdt, tdt
  real, intent(inout), dimension(:,:,:,:) :: &
    rdt
  !real, intent(inout), dimension(:,:,:,:) :: &   ! Mellor-Yamada, use this in offline mode
  real, intent(inout), dimension(:,:,:,ntp+1:) :: &
    rdiag

  real, intent(inout), dimension(:,:) :: &
    pbltop

!---------------------------------------------------------------------
! Arguments (Intent out)
!
!   pbltop 		- PBL depth (m)								, dimension (nlon, nlat) 
!   kpbl_edmf 		- vertical index of PBL top (none)					, dimension (nlon, nlat) 
!   qadt_edmf 		- cloud fraction tendency from edmf_mynn (1/s), 		  	, dimension (nlon,nlat,nlay)
!   qldt_edmf 		- cloud liquid specific humidity tendency from edmf_mynn (kg/kg/s)	, dimension (nlon,nlat,nlay)
!   qidt_edmf 		- cloud ice    specific humidity tendency from edmf_mynn (kg/kg/s)	, dimension (nlon,nlat,nlay)
!   diff_t_edmf		- eddy diffusion coefficient for heat (K m/s)			  	, dimension (nlon,nlat,nlay)
!   diff_m_edmf		- eddy diffusion coefficient for momentum (m2/s)		  	, dimension (nlon,nlat,nlay)
!   option_edmf2ls_mp 	- options for linkage to Tiedtke, same as "do_option_edmf2ls_mp" in the namelist
!
!   Note
!     1. diff_t_edmf is at half levels, although AM4 use dimension of nlay instead of nlay+1
!        k=1 is at 1/2 level and k=N is at N-1/2 level 
!---------------------------------------------------------------------

  integer, intent(out), dimension(:,:) :: &  ! (lon, lat)
    kpbl_edmf
 
  real, intent(out), dimension(:,:,:) :: &   ! (lon, lat, nlay)
    qadt_edmf, qldt_edmf, qidt_edmf, &   
    dqa_edmf,  dql_edmf, dqi_edmf,  &
    diff_t_edmf, diff_m_edmf

  integer, intent(out) :: &
    option_edmf2ls_mp

  real, intent(out), dimension(:,:,:) :: &  ! (lon, lat, nlay+1)
    edmf_mc_half

  real, intent(out), dimension(:,:,:) :: &  ! (lon, lat, nlay)
    edmf_mc_full, edmf_moist_area, edmf_dry_area, edmf_moist_humidity, edmf_dry_humidity

!---------------------------------------------------------------------
! local variables  
!---------------------------------------------------------------------
  type(edmf_input_type)  :: Input_edmf            ! Input variable for mynn
  type(edmf_output_type) :: Output_edmf           ! Output variable from mynn
  type(am4_edmf_output_type) :: am4_Output_edmf   ! convert Output_edmf to AM4 needed variables and dimension

  logical used
  logical do_writeout_column
  real    :: lat_lower, lat_upper, lon_lower, lon_upper, lat_temp, lon_temp
  real    :: tt1, tt2
  integer :: i,j,k
  integer :: ix,jx,kx

  real, dimension (size(Physics_input_block%t,1), &
                   size(Physics_input_block%t,2), &
                   size(Physics_input_block%t,3)) :: &
          tmp_3d, diag_full

  real, dimension (size(Physics_input_block%t,1), &
                   size(Physics_input_block%t,2), &
                   size(Physics_input_block%t,3)+1) :: &
          diag_half

  logical, dimension (size(Physics_input_block%t,1), &
                      size(Physics_input_block%t,2), &
                      size(Physics_input_block%t,3)+1) :: &
          lmask_half 

  real, dimension (size(Physics_input_block%t,1), &
                   size(Physics_input_block%t,3), &
                   size(Physics_input_block%t,2)) :: &   ! (i,k,j)
          num_DET, num_nDET_pENT, num_nDET_zENT

  real, dimension (size(Physics_input_block%t,1), &
                   size(Physics_input_block%t,3)+1, &
                   size(Physics_input_block%t,2)) :: &   ! (i,k,j)
          num_updraft

  type(randomNumberStream), dimension(size(Physics_input_block%t,1), &
                                      size(Physics_input_block%t,2)) :: &     ! dimension (nlon,nlat)
    streams

  real, dimension (size(Physics_input_block%t,1), &
                   size(Physics_input_block%t,2)) :: &  ! dimension (nlon,nlat)
    rlz_ratio, rlz_tracer

  real :: randomNumbers

!-------------------------
  ix = size(Physics_input_block%t,1)
  jx = size(Physics_input_block%t,2)
  kx = size(Physics_input_block%t,3)  

  if (option_stoch_entrain.eq."Poisson_knuth") then
    call get_random_number_streams ( is, js, Time_next, Physics_input_block%t(:,:,kx), streams)
  endif

!-----------------
! check namelist
!----------------
  if (do_edmf_mynn_in_physics.eq."down" .and. option_surface_flux.ne."star") then
    call error_mesg( ' edmf_mynn',     &
                     ' when do_edmf_mynn_in_physics=down, option_surface_flux must be "star" in the edmf_mynn_nml',&
                     FATAL )
  endif

  !--- MYNN EDonly configuration. 
  !    MYNN only returns ED coefficients and all diffusion is handled by AM4 vertical diffusion module, similar to Lock
  if (do_return_edmf_mynn_diff_only) then

    !!!!!!!!!!!!!!!!!
    if (do_edmf_mynn_in_physics.eq."down" .and. option_surface_flux.eq."star") then
      tt1=0.
    else
      call error_mesg( ' edmf_mynn',     &
                       ' when do_return_edmf_mynn_diff_only=T, do_edmf_mynn_in_physics must be "down" and option_surface_flux must be "star" in the edmf_mynn_nml',&
                       FATAL )
    endif

    !!!!!!!!!!!!!!!!!
    if (edmf_type.eq.1 .and. bl_mynn_edmf.eq.0) then
      tt1=0.
    else
      call error_mesg( ' edmf_mynn',     &
                       ' when do_return_edmf_mynn_diff_only=T, edmf_type must be 1 and bl_mynn_edmf must be 0',&
                       FATAL )
    endif
  end if   ! end if of do_return_edmf_mynn_diff_only

  !--- Method 3, EDonly configuration. 
  !    approximate ql and qi tendencies from AM4 ED and MF terms, and then recover T and qv tendencies
  if (edmf_type.eq.2 .and. bl_mynn_edmf.eq.0) then
    !!!!!!!!!!!!!!!!!
    if (do_tracers_selective.ne.0 .or. do_option_edmf2ls_mp.ne.99) then
      call error_mesg( ' edmf_mynn',     &
                       ' For M3_EDonly (edmf_type=2 .and. bl_mynn_edmf=0), do_tracers_selective must be 0 and do_option_edmf2ls_mp must be 99',&
                       FATAL )
    endif
  end if   ! end if of do_return_edmf_mynn_diff_only
  
  !--- Method 3, EDMF configuration. 
  !    approximate ql and qi tendencies from AM4 ED and MF terms, and then recover T and qv tendencies
  if (edmf_type.eq.2 .and. bl_mynn_edmf.eq.3) then
    !!!!!!!!!!!!!!!!!
    if (do_tracers_selective.ne.6 .or. do_option_edmf2ls_mp.ne.4) then
      call error_mesg( ' edmf_mynn',     &
                       ' For M3_EDMF (edmf_type=2 .and. bl_mynn_edmf=3), do_tracers_selective must be 6 and do_option_edmf2ls_mp must be 4',&
                       FATAL )
    endif
  end if   ! end if of do_return_edmf_mynn_diff_only


!! debug01
!write(6,*) 'edmf_mynn, beginning'
!!write(6,*) 'Physics_input_block%omega',Physics_input_block%omega
!write(6,*) 'initflag,',initflag
!write(6,*) 'nQke, rdiag(:,:,:,nQke)',nQke, rdiag(:,:,:,nQke)
!write(6,*) 'rdiag(:,:,:,nel_pbl)',rdiag(:,:,:,nel_pbl)
!write(6,*) 'rdiag(:,:,:,ncldfra_bl)',rdiag(:,:,:,nqc_bl)
!write(6,*) 'rdiag(:,:,:,nqc_bl)',rdiag(:,:,:,nqc_bl)
!write(6,*) 'rdiag(:,:,:,nSh3D)',rdiag(:,:,:,nSh3D)
!!write(6,*) 'ntp,nQke, nSh3D, nel_pbl, ncldfra_bl, nqc_bl',ntp,nQke, nSh3D, nel_pbl, ncldfra_bl, nqc_bl

!-------------------------------------------------------------------------
!  determine whether writing out the selected column
!-------------------------------------------------------------------------
  do_writeout_column = .false.
  if (do_writeout_column_nml) then

    !--- for global simulations
    if (ii_write.ne.-999 .and. jj_write.ne.-999) then
      do_writeout_column = .true.

      if (lat_write.ne.-999.99 .and. lon_write.ne.-999.99) then

        lat_lower = lat_write - lat_range
        lat_upper = lat_write + lat_range
        lon_lower = lon_write - lon_range
        lon_upper = lon_write + lon_range

        if (lat_lower.gt.lat_upper) then
          lat_temp  = lat_upper
          lat_upper = lat_lower
          lat_lower = lat_temp
        endif

        if (lon_lower.gt.lon_upper) then
          lon_temp  = lon_upper
          lon_upper = lon_lower
          lon_lower = lon_temp
        endif

        if (lat (ii_write,jj_write).gt.lat_lower .and. lat (ii_write,jj_write).lt.lat_upper .and. &
            lon (ii_write,jj_write).gt.lon_lower .and. lon (ii_write,jj_write).lt.lon_upper ) then
          do_writeout_column = .true.
        else
          do_writeout_column = .false.
        endif
      endif
    endif

    !--- SCM
    if (ii_write.eq.0 .and. jj_write.eq.0) then
      do_writeout_column = .true.
      ii_write = 1
      jj_write = 1
    endif
  endif  ! end if of do_writeout_column_nml

!if (do_writeout_column) then
!  write(6,*) 'initflag, in',initflag
!  write(6,*) 'rdiag(:,:,:,nQke), in',rdiag(ii_write,jj_write,:,nQke)
!endif

!---------------------------------------------------------------------
! allocate input and output variables for the EDMF-MYNN program
!---------------------------------------------------------------------
  call edmf_alloc ( &
              is, ie, js, je, npz, Time_next, dt, frac_land, area, u_star,  &
              b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, rdiag, udt, vdt, tdt, rdt, &
              rdiag(:,:,:,nQke), rdiag(:,:,:,nel_pbl), rdiag(:,:,:,ncldfra_bl), rdiag(:,:,:,nqc_bl), rdiag(:,:,:,nSh3D), &
              Input_edmf, Output_edmf, am4_Output_edmf)

!---------------------------------------------------------------------
! run the EDMF-MYNN program
!---------------------------------------------------------------------

   call mynn_bl_driver(            &
       &initflag=initflag,grav_settling=grav_settling,         &
       &delt=Input_edmf%delt,dz=Input_edmf%dz,dx=Input_edmf%dx,znt=Input_edmf%znt,                 &
       &u=Input_edmf%u,v=Input_edmf%v,w=Input_edmf%w,th=Input_edmf%th,qv=Input_edmf%qv,             &
       &ql=Input_edmf%ql,qi=Input_edmf%qi,cc=Input_edmf%qa,qni=Input_edmf%qni,qnc=Input_edmf%qnc,                    &
       &p=Input_edmf%p,exner=Input_edmf%exner,rho=Input_edmf%rho,T3D=Input_edmf%T3D,                &
       &xland=Input_edmf%xland,ts=Input_edmf%ts,qsfc=Input_edmf%qsfc,qcg=Input_edmf%qcg,ps=Input_edmf%ps,           &
       &ust=Input_edmf%ust,ch=Input_edmf%ch,hfx=Input_edmf%hfx,qfx=Input_edmf%qfx,rmol=Input_edmf%rmol,wspd=Input_edmf%wspd,       &
       &uoce=Input_edmf%uoce,voce=Input_edmf%voce,                      & 
       &vdfg=Input_edmf%vdfg,                           & 
       &qke=Output_edmf%Qke,                    &
       &Tsq=Output_edmf%Tsq,Qsq=Output_edmf%Qsq,Cov=Output_edmf%Cov,                    &
       &RUBLTEN=Output_edmf%RUBLTEN,RVBLTEN=Output_edmf%RVBLTEN,RTHBLTEN=Output_edmf%RTHBLTEN,       &
       &RQVBLTEN=Output_edmf%RQVBLTEN,RQLBLTEN=Output_edmf%RQLBLTEN,RQIBLTEN=Output_edmf%RQIBLTEN,     &
       &RQNIBLTEN=Output_edmf%RQNIBLTEN,                      &
       &RCCBLTEN=Output_edmf%RCCBLTEN, RTHLBLTEN=Output_edmf%RTHLBLTEN, RQTBLTEN=Output_edmf%RQTBLTEN, &  ! yhc_mynn add
       &qa_before_mix=Output_edmf%qa_before_mix, ql_before_mix=Output_edmf%ql_before_mix, qi_before_mix=Output_edmf%qi_before_mix, thl_before_mix=Output_edmf%thl_before_mix, qt_before_mix=Output_edmf%qt_before_mix, th_before_mix=Output_edmf%th_before_mix, &      ! yhc_mynn add
       &qa_after_mix=Output_edmf%qa_after_mix, ql_after_mix=Output_edmf%ql_after_mix, qi_after_mix=Output_edmf%qi_after_mix, thl_after_mix=Output_edmf%thl_after_mix, qt_after_mix=Output_edmf%qt_after_mix, th_after_mix=Output_edmf%th_after_mix,        &      ! yhc_mynn add
        &qa_before_pdf=Output_edmf%qa_before_pdf, ql_before_pdf=Output_edmf%ql_before_pdf, qi_before_pdf=Output_edmf%qi_before_pdf, & ! yhc_mynn add
       &Q_ql=Output_edmf%Q_ql, Q_qi=Output_edmf%Q_qi, Q_a=Output_edmf%Q_qa,   &  ! yhc_mynn add
       &Q_ql_adv=Output_edmf%Q_ql_adv, Q_qi_adv=Output_edmf%Q_qi_adv, Q_a_adv=Output_edmf%Q_qa_adv, Q_ql_eddy=Output_edmf%Q_ql_eddy, Q_qi_eddy=Output_edmf%Q_qi_eddy, Q_a_eddy=Output_edmf%Q_qa_eddy, Q_ql_ent=Output_edmf%Q_ql_ent, Q_qi_ent=Output_edmf%Q_qi_ent, Q_a_ent=Output_edmf%Q_qa_ent, Q_ql_det=Output_edmf%Q_ql_det, Q_qi_det=Output_edmf%Q_qi_det, Q_a_det=Output_edmf%Q_qa_det, Q_ql_sub=Output_edmf%Q_ql_sub, Q_qi_sub=Output_edmf%Q_qi_sub, Q_a_sub=Output_edmf%Q_qa_sub, &  ! yhc_mynn add
       &a_moist_half=Output_edmf%a_moist_half, mf_moist_half=Output_edmf%mf_moist_half, qv_moist_half=Output_edmf%qv_moist_half, a_moist_full=Output_edmf%a_moist_full, mf_moist_full=Output_edmf%mf_moist_full, qv_moist_full=Output_edmf%qv_moist_full, &  ! yhc 2021-09-08
       &a_dry_half=Output_edmf%a_dry_half, mf_dry_half=Output_edmf%mf_dry_half, qv_dry_half=Output_edmf%qv_dry_half, a_dry_full=Output_edmf%a_dry_full, mf_dry_full=Output_edmf%mf_dry_full, qv_dry_full=Output_edmf%qv_dry_full, &            ! yhc 2021-09-08
       &mf_all_half=Output_edmf%mf_all_half, mf_all_full=Output_edmf%mf_all_full, &      ! yhc 2021-09-08
       &num_updraft=num_updraft, num_DET=num_DET, num_nDET_pENT=num_nDET_pENT, num_nDET_zENT=num_nDET_zENT, &            ! yhc 2021-09-08
       &streams=streams, &  ! yhc 2021-11-18
       &exch_h=Output_edmf%exch_h,exch_m=Output_edmf%exch_m,                  &
       &pblh=Output_edmf%Pblh,kpbl=Output_edmf%kpbl,                      & 
       &el_pbl=Output_edmf%el_pbl,                         &
       &dqke=Output_edmf%dqke,qwt=Output_edmf%qWT,qshear=Output_edmf%qSHEAR,qbuoy=Output_edmf%qBUOY,qdiss=Output_edmf%qDISS,    &
       &bl_mynn_tkebudget=bl_mynn_tkebudget,              &
       &bl_mynn_cloudpdf=bl_mynn_cloudpdf,sh3D=Output_edmf%Sh3D,          &
       &bl_mynn_mixlength=bl_mynn_mixlength,              &
       &icloud_bl=icloud_bl,qc_bl=Output_edmf%qc_bl,cldfra_bl=Output_edmf%cldfra_bl,      &
       &bl_mynn_edmf=bl_mynn_edmf,                   &
       &bl_mynn_edmf_dd=bl_mynn_edmf_dd,                   &
       &bl_mynn_edmf_mom=bl_mynn_edmf_mom,bl_mynn_edmf_tke=bl_mynn_edmf_tke, &
       &bl_mynn_edmf_part=bl_mynn_edmf_part,bl_mynn_edmf_Lent=bl_mynn_edmf_Lent,&
       &bl_mynn_cloudmix=bl_mynn_cloudmix,bl_mynn_mixqt=bl_mynn_mixqt, &
       &edmf_a=Output_edmf%edmf_a,edmf_w=Output_edmf%edmf_w,edmf_qt=Output_edmf%edmf_qt,          &
       &edmf_thl=Output_edmf%edmf_thl,edmf_ent=Output_edmf%edmf_ent, edmf_det=Output_edmf%edmf_det, edmf_qc=Output_edmf%edmf_qc,      &
       &edmf_debug1=Output_edmf%edmf_debug1,edmf_debug2=Output_edmf%edmf_debug2,        &
       &edmf_debug3=Output_edmf%edmf_debug3,edmf_debug4=Output_edmf%edmf_debug4,        &
       &edmf_a_dd=Output_edmf%edmf_a_dd,edmf_w_dd=Output_edmf%edmf_w_dd,edmf_qt_dd=Output_edmf%edmf_qt_dd,    &
       &edmf_thl_dd=Output_edmf%edmf_thl_dd,edmf_ent_dd=Output_edmf%edmf_ent_dd,edmf_qc_dd=Output_edmf%edmf_qc_dd,&
       &mynn_ql=Output_edmf%mynn_ql,                        &
       &ktop_shallow=Output_edmf%ktop_shallow,    &
       &RTHRATEN=Output_edmf%RTHRATEN,                       &
       &FLAG_QI=FLAG_QI,FLAG_QNI=FLAG_QNI,FLAG_QC=FLAG_QC,FLAG_QNC=FLAG_QNC &
       &,IDS=Input_edmf%IDS,IDE=Input_edmf%IDE,JDS=Input_edmf%JDS,JDE=Input_edmf%JDE,KDS=Input_edmf%KDS,KDE=Input_edmf%KDE        &
       &,IMS=Input_edmf%IMS,IME=Input_edmf%IME,JMS=Input_edmf%JMS,JME=Input_edmf%JME,KMS=Input_edmf%KMS,KME=Input_edmf%KME        &
       &,ITS=Input_edmf%ITS,ITE=Input_edmf%ITE,JTS=Input_edmf%JTS,JTE=Input_edmf%JTE,KTS=Input_edmf%KTS,KTE=Input_edmf%KTE)

  !--- SCM, set initflag
  if (initflag == 1) then
    initflag = 0          ! no initialization
  endif 

!---------------------------------------------------------------------
! recover dry variable tendencies from mynn_edmf
!---------------------------------------------------------------------
  call modify_mynn_edmf_tendencies( is, ie, js, je, Time_next, dt,  &
                                    do_writeout_column,         &
                                    Physics_input_block, Input_edmf, rdt_mynn_ed_am4, &
                                    size(Physics_input_block%t,1), size(Physics_input_block%t,2), size(Physics_input_block%t,3), &
                                    Output_edmf, rlz_ratio, rlz_tracer)

!---------------------------------------------------------------------
! process the outputs from the EDMF-MYNN program
!---------------------------------------------------------------------

  !--- convert Output_edmf to am4_Output_edmf
  call convert_edmf_to_am4_array (Physics_input_block, size(Physics_input_block%t,1), size(Physics_input_block%t,2), size(Physics_input_block%t,3), &
                                  Input_edmf, Output_edmf, am4_Output_edmf, rdiag, &
                                  rdiag(:,:,:,nQke), rdiag(:,:,:,nel_pbl), rdiag(:,:,:,ncldfra_bl), rdiag(:,:,:,nqc_bl), rdiag(:,:,:,nSh3D) )

!! debug01
!write(6,*) 'edmf_mynn, after mynn'
!!write(6,*) 'Physics_input_block%omega',Physics_input_block%omega
!write(6,*) 'initflag,',initflag
!write(6,*) 'nQke, rdiag(:,:,:,nQke)',nQke, rdiag(:,:,:,nQke)
!write(6,*) 'rdiag(:,:,:,nel_pbl)',rdiag(:,:,:,nel_pbl)
!write(6,*) 'rdiag(:,:,:,ncldfra_bl)',rdiag(:,:,:,nqc_bl)
!write(6,*) 'rdiag(:,:,:,nqc_bl)',rdiag(:,:,:,nqc_bl)
!write(6,*) 'rdiag(:,:,:,nSh3D)',rdiag(:,:,:,nSh3D)

!---------------------------------------------------------------------
! return EDMF-MYNN terms to GFDL model
!---------------------------------------------------------------------

  !--- initialize return variables 
  option_edmf2ls_mp = 0
  qadt_edmf   = 0.
  qldt_edmf   = 0.
  qidt_edmf   = 0.
  dqa_edmf    = 0.
  dql_edmf    = 0.
  dqi_edmf    = 0.
  diff_t_edmf = 0.
  diff_m_edmf = 0.
  edmf_mc_full         = 0.
  edmf_mc_half         = 0.
  edmf_moist_area      = 0.
  edmf_moist_humidity  = 0.
  edmf_dry_area        = 0.
  edmf_dry_humidity    = 0.

  if (.not.do_edmf_mynn_diagnostic) then

    ! return edmf_mynn diffusion coefficients and let vert_diff do the diffusion rather than in edmf_mynn
    if (do_return_edmf_mynn_diff_only) then
      diff_t_edmf (:,:,:) = am4_Output_edmf%diff_t_edmf (:,:,:)      ! diffusion coefficient for heat (K m/s)
      diff_m_edmf (:,:,:) = am4_Output_edmf%diff_m_edmf (:,:,:)      ! diffusion coefficient for heat (m2/s)
      pbltop      (:,:)   = am4_Output_edmf%pbltop      (:,:)        ! PBL depth (m)
      kpbl_edmf   (:,:)   = am4_Output_edmf%kpbl_edmf   (:,:)        ! index of PBL top (none)

    ! update tendencies
    else
      !--- updated tendencies
      udt(:,:,:) = udt(:,:,:) + am4_Output_edmf%udt_edmf(:,:,:)
      vdt(:,:,:) = vdt(:,:,:) + am4_Output_edmf%vdt_edmf(:,:,:)
      tdt(:,:,:) = tdt(:,:,:) + am4_Output_edmf%tdt_edmf(:,:,:)

      rdt(:,:,:,nsphum) = rdt(:,:,:,nsphum) + am4_Output_edmf%qdt_edmf(:,:,:)

      !--- return PBL depth, index of PBL top, and diffusion coefficient for heat
      pbltop      (:,:)   = am4_Output_edmf%pbltop      (:,:)        ! PBL depth (m)
      kpbl_edmf   (:,:)   = am4_Output_edmf%kpbl_edmf   (:,:)        ! index of PBL top (none)
      diff_t_edmf (:,:,:) = am4_Output_edmf%diff_t_edmf (:,:,:)      ! diffusion coefficient for heat (K m/s)
      diff_m_edmf (:,:,:) = am4_Output_edmf%diff_m_edmf (:,:,:)      ! diffusion coefficient for heat (m2/s)
    endif  ! end if do_return_edmf_mynn_diff_only

    !--- set edmf to ls_mp
    option_edmf2ls_mp = do_option_edmf2ls_mp

    ! accumulate EDMF cloud tendencies to the model tendencies. No variables are passed to Tiedtke.
    if (option_edmf2ls_mp.eq.0) then
      rdt(:,:,:,nqa)  = rdt(:,:,:,nqa) + am4_Output_edmf%qadt_edmf(:,:,:)  
      rdt(:,:,:,nql)  = rdt(:,:,:,nql) + am4_Output_edmf%qldt_edmf(:,:,:)  
      rdt(:,:,:,nqi)  = rdt(:,:,:,nqi) + am4_Output_edmf%qidt_edmf(:,:,:)

    ! save EDMF cloud tendencies that would be included in Tiedtke
    elseif (option_edmf2ls_mp.eq.1 .or. option_edmf2ls_mp.eq.2) then
      qadt_edmf     (:,:,:) = am4_Output_edmf%qadt_edmf(:,:,:)
      qldt_edmf     (:,:,:) = am4_Output_edmf%qldt_edmf(:,:,:)
      qidt_edmf     (:,:,:) = am4_Output_edmf%qidt_edmf(:,:,:)

    ! set EDMF cloud tendencies to zeros as the changes of qdt_edmf already include these changes
    ! the qdt_edmf is modified in convert_edmf_to_am4_array
    elseif (option_edmf2ls_mp.eq.3) then   ! the changes of qdt_edmf is handled in convert_edmf_to_am4_array
      qadt_edmf   = 0.
      qldt_edmf   = 0.
      qidt_edmf   = 0.
    
    ! accumulate EDMF cloud tendencies to the model tendencies. Pass MF mass flux to Tiedtke
    elseif (option_edmf2ls_mp.eq.4) then 
      rdt(:,:,:,nqa)  = rdt(:,:,:,nqa) + am4_Output_edmf%qadt_edmf(:,:,:)  
      rdt(:,:,:,nql)  = rdt(:,:,:,nql) + am4_Output_edmf%qldt_edmf(:,:,:)  
      rdt(:,:,:,nqi)  = rdt(:,:,:,nqi) + am4_Output_edmf%qidt_edmf(:,:,:)

      ! for testing, assume cloud droplet radius (rc_MF) and then compute cloud liquid number tendency
      am4_Output_edmf%qndt_edmf(:,:,:) = am4_Output_edmf%qldt_edmf(:,:,:) / dens_h2o / (4./3.*pi*rc_MF**3)  
      rdt(:,:,:,nqn)  = rdt(:,:,:,nqn) + am4_Output_edmf%qndt_edmf(:,:,:)

      edmf_mc_half (:,:,:) = am4_Output_edmf%mf_all_half (:,:,:)
      edmf_mc_full (:,:,:) = am4_Output_edmf%mf_all_full (:,:,:)

    ! accumulate EDMF cloud tendencies to the model tendencies. Pass MF mass flux, moist area and humidity to Tiedtke
    elseif (option_edmf2ls_mp.eq.5) then 
      rdt(:,:,:,nqa)  = rdt(:,:,:,nqa) + am4_Output_edmf%qadt_edmf(:,:,:)  
      rdt(:,:,:,nql)  = rdt(:,:,:,nql) + am4_Output_edmf%qldt_edmf(:,:,:)  
      rdt(:,:,:,nqi)  = rdt(:,:,:,nqi) + am4_Output_edmf%qidt_edmf(:,:,:)

      edmf_mc_full        (:,:,:) = am4_Output_edmf%mf_all_half  (:,:,:)
      edmf_mc_half        (:,:,:) = am4_Output_edmf%mf_all_full  (:,:,:)
      edmf_moist_area     (:,:,:) = am4_Output_edmf%a_moist_full (:,:,:)
      edmf_moist_humidity (:,:,:) = am4_Output_edmf%qv_moist_full(:,:,:)
      edmf_dry_area       (:,:,:) = am4_Output_edmf%a_dry_full   (:,:,:)
      edmf_dry_humidity   (:,:,:) = am4_Output_edmf%qv_dry_full  (:,:,:)

    ! do not change cloud tendencies 
    elseif (option_edmf2ls_mp.eq.99) then 
      qadt_edmf   = 0.
      qldt_edmf   = 0.
      qidt_edmf   = 0.

    ! set to zeros if option_edmf2ls_mp is not supported
    else
      call error_mesg( ' edmf_mynn',     &
                       ' do_option_edmf2ls_mp is not valid',&
                       FATAL )

    endif  ! end if of option_edmf2ls_mp

    dqa_edmf      (:,:,:) = qadt_edmf(:,:,:) * dt
    dql_edmf      (:,:,:) = qldt_edmf(:,:,:) * dt
    dqi_edmf      (:,:,:) = qidt_edmf(:,:,:) * dt
    
  end if  ! end if of do_edmf_mynn_diagnostic

  !--- write out EDMF-MYNN input and output fields for debugging purpose
  call edmf_writeout_column ( &
              do_writeout_column, &
              is, ie, js, je, npz, Time_next, dt, lon, lat, frac_land, area, u_star,  &
              b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, rdt_mynn_ed_am4,  &
              rdiag(:,:,:,nQke), rdiag(:,:,:,nel_pbl), rdiag(:,:,:,ncldfra_bl), rdiag(:,:,:,nqc_bl), rdiag(:,:,:,nSh3D), &
              Input_edmf, Output_edmf, am4_Output_edmf, rdiag)

!---------------------------------------------------------------------
! write out fields to history files
!---------------------------------------------------------------------

      !--- set up local mask for fields without surface data
      lmask_half(:,:,1:kx) = .true.
      lmask_half(:,:,kx+1) = .false.

!send_data
!------- zonal wind stress (units: kg/m/s2) at one level -------
      if ( id_u_flux > 0) then
        used = send_data (id_u_flux, u_flux, Time_next, is, js )
      endif

!------- meridional wind stress (units: kg/m/s2) at one level -------
      if ( id_v_flux > 0) then
        used = send_data (id_v_flux, v_flux, Time_next, is, js )
      endif

!------- u_star from u_flux and v_flux (units: m/s) at one level -------
      if ( id_u_star_updated > 0) then
        used = send_data (id_u_star_updated, Input_edmf%u_star_updated, Time_next, is, js )
      endif

!------- sensible heat flux from star (units: W/m2) at one level -------
      if ( id_shflx_star > 0) then
        used = send_data (id_shflx_star, Input_edmf%shflx_star, Time_next, is, js )
      endif

!------- evaporation flux from star (units: kg/m2/s) at one level -------
      if ( id_lhflx_star > 0) then
        used = send_data (id_lhflx_star, Input_edmf%lhflx_star, Time_next, is, js )
      endif

!------- kinematic virtual temperature flux from star (units: K m/s) at one level -------
      if ( id_w1_thv1_surf_star > 0) then
        used = send_data (id_w1_thv1_surf_star, Input_edmf%w1_thv1_surf_star, Time_next, is, js )
      endif

!------- kinematic virtual temperature flux from updated fluxes (units: K m/s) at one level -------
      if ( id_w1_thv1_surf_updated > 0) then
        used = send_data (id_w1_thv1_surf_updated, Input_edmf%w1_th1_surf_updated, Time_next, is, js )
      endif

!------- Obukhov length from star (units: m) at one level -------
      if ( id_Obukhov_length_star > 0) then
        used = send_data (id_Obukhov_length_star, Input_edmf%Obukhov_length_star, Time_next, is, js )
      endif

!------- Obukhov length from updated fluxes (units: m) at one level -------
      if ( id_Obukhov_length_updated > 0) then
        used = send_data (id_Obukhov_length_updated, Input_edmf%Obukhov_length_updated, Time_next, is, js )
      endif

!------- turbulent kinetic energy (units: m2/s2) at full level -------
      if ( id_tke_edmf > 0) then
        used = send_data (id_tke_edmf, am4_Output_edmf%tke, Time_next, is, js, 1 )
      endif

!------- variance of theta_l (units: K^2) at full level -------
      if ( id_Tsq > 0) then
        used = send_data (id_Tsq, am4_Output_edmf%Tsq, Time_next, is, js, 1 )
      endif

!------- covariance of theta_l and q_t (units: none) at full level -------
      if ( id_Cov_thl_qt > 0) then
        used = send_data (id_Cov_thl_qt, am4_Output_edmf%Cov_thl_qt, Time_next, is, js, 1 )
      endif

!------- u tendency from edmf_mynn (units: m/s2) at full level -------
      if ( id_udt_edmf > 0) then
        used = send_data (id_udt_edmf, am4_Output_edmf%udt_edmf, Time_next, is, js, 1 )
      endif

!------- v tendency from edmf_mynn (units: m/s2) at full level -------
      if ( id_vdt_edmf > 0) then
        used = send_data (id_vdt_edmf, am4_Output_edmf%vdt_edmf, Time_next, is, js, 1 )
      endif

!------- t tendency from edmf_mynn (units: K/s) at full level -------
      if ( id_tdt_edmf > 0) then
        used = send_data (id_tdt_edmf, am4_Output_edmf%tdt_edmf, Time_next, is, js, 1 )
      endif

!------- q tendency from edmf_mynn (units: kg/kg/s) at full level -------
      if ( id_qdt_edmf > 0) then
        used = send_data (id_qdt_edmf, am4_Output_edmf%qdt_edmf, Time_next, is, js, 1 )
      endif

!------- qi tendency from edmf_mynn (units: kg/kg/s) at full level -------
      if ( id_qidt_edmf > 0) then
        used = send_data (id_qidt_edmf, am4_Output_edmf%qidt_edmf, Time_next, is, js, 1 )
      endif

!------- qc tendency from edmf_mynn (units: kg/kg/s) at full level -------
      if ( id_qldt_edmf > 0) then
        used = send_data (id_qldt_edmf, am4_Output_edmf%qldt_edmf, Time_next, is, js, 1 )
      endif

!------- cldfrac tendency from edmf_mynn (units: 1/s) at full level -------
      if ( id_qadt_edmf > 0) then
        used = send_data (id_qadt_edmf, am4_Output_edmf%qadt_edmf, Time_next, is, js, 1 )
      endif

!------- cloud droplet number tendency from edmf_mynn (units: 1/kg/s) at full level -------
      if ( id_qndt_edmf > 0) then
        used = send_data (id_qndt_edmf, am4_Output_edmf%qndt_edmf, Time_next, is, js, 1 )
      endif

!------- updraft area (units: none) at half level -------
      if ( id_edmf_a > 0) then
        call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_a, diag_half)
        used = send_data (id_edmf_a, diag_half, Time_next, is, js, 1 )
      endif

!------- vertical velocity of updrafts (units: m/s) at half level -------
      if ( id_edmf_w > 0) then
        call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_w, diag_half)
        used = send_data (id_edmf_w, diag_half, Time_next, is, js, 1 )
      endif

!------- qt in updrafts (units: kg/kg) at half level -------
      if ( id_edmf_qt > 0) then
        call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_qt, diag_half)
        used = send_data (id_edmf_qt, diag_half, Time_next, is, js, 1 )
      endif

!------- thl in updrafts (units: K) at half level -------
      if ( id_edmf_thl > 0) then
        call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_thl, diag_half)
        used = send_data (id_edmf_thl, diag_half, Time_next, is, js, 1 )
      endif

!------- entrainment in updrafts (units: 1/m) at full level -------
      if ( id_edmf_ent > 0) then
        call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%edmf_ent, diag_full)
        used = send_data (id_edmf_ent, diag_full, Time_next, is, js, 1 )
      endif

!------- dentrainment in updrafts (units: 1/m) at full level -------
      if ( id_edmf_det > 0) then
        call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%edmf_det, diag_full)
        used = send_data (id_edmf_det, diag_full, Time_next, is, js, 1 )
      endif

!------- qc in updrafts (units: kg/kg) at half level -------
      if ( id_edmf_qc > 0) then
        call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_qc, diag_half)
        used = send_data (id_edmf_qc, diag_half, Time_next, is, js, 1 )
      endif

!------- theta_li in edmf_mynn (units: K) at full level -------
      if ( id_thl_edmf > 0) then
        used = send_data (id_thl_edmf, am4_Output_edmf%thl_edmf, Time_next, is, js, 1 )
      endif

!------- qt in edmf_mynn (units: kg/kg) at full level -------
      if ( id_qt_edmf > 0) then
        used = send_data (id_qt_edmf, am4_Output_edmf%qt_edmf, Time_next, is, js, 1 )
      endif

!------- cloud fraction in edmf_mynn (units: fraction) at full level -------
      if ( id_cldfra_bl > 0) then
        used = send_data (id_cldfra_bl, am4_Output_edmf%cldfra_bl, Time_next, is, js, 1 )
      endif

!------- liquid water mixing ratio in edmf_mynn (units: mynn) at full level -------
      if ( id_qc_bl > 0) then
        used = send_data (id_qc_bl, am4_Output_edmf%qc_bl, Time_next, is, js, 1 )
      endif

!------- depth of planetary boundary layer (units: m) at one level -------
      if ( id_z_pbl > 0) then
        used = send_data (id_z_pbl, pbltop, Time_next, is, js )
      endif

!------- PBL depth from edmf_mynn (units: m) at one level -------
      if ( id_z_pbl_edmf > 0) then
        used = send_data (id_z_pbl_edmf, am4_Output_edmf%pbltop, Time_next, is, js )
      endif

!------- qt tendency from edmf_mynn (units: kg/kg/s) at full level -------
      if ( id_qtdt_edmf > 0) then
        used = send_data (id_qtdt_edmf, am4_Output_edmf%qtdt_edmf, Time_next, is, js, 1 )
      endif

!------- thl tendency from edmf_mynn (units: K/s) at full level -------
      if ( id_thldt_edmf > 0) then
        used = send_data (id_thldt_edmf, am4_Output_edmf%thldt_edmf, Time_next, is, js, 1 )
      endif

!------- heat diff coeffs from edmf_mynn (units: K/m/s) at half level -------
      if ( id_diff_t_edmf > 0) then
        ! the dimension size of am4_Output_edmf%diff_t_edmf is kx, but it is on half level. Write out on half levels
        diag_half(:,:,kx+1) = 0.
        diag_half(:,:,1:kx) = am4_Output_edmf%diff_t_edmf(:,:,1:kx)
        used = send_data (id_diff_t_edmf, diag_half, Time_next, is, js, 1, mask=lmask_half )
        !used = send_data (id_diff_t_edmf, am4_Output_edmf%diff_t_edmf, Time_next, is, js, 1 )
      endif

!------- momentum diff coeffs from edmf_mynn (units: m2/s) at half level -------
      if ( id_diff_m_edmf > 0) then
        ! the dimension size of am4_Output_edmf%diff_m_edmf is kx, but it is on half level. Write out on half levels
        diag_half(:,:,kx+1) = 0.
        diag_half(:,:,1:kx) = am4_Output_edmf%diff_m_edmf(:,:,1:kx)
        used = send_data (id_diff_m_edmf, diag_half, Time_next, is, js, 1, mask=lmask_half )
        !used = send_data (id_diff_m_edmf, am4_Output_edmf%diff_m_edmf, Time_next, is, js, 1 )
      endif

!------- mixing length in edmf_mynn (units: m) at full level -------
      if ( id_el_edmf > 0) then
        used = send_data (id_el_edmf, am4_Output_edmf%el_edmf, Time_next, is, js, 1 )
      endif

!------- T input to edmf_mynn (units: K) at full level -------
      if ( id_t_input > 0) then
        used = send_data (id_t_input, am4_Output_edmf%t_input, Time_next, is, js, 1 )
      endif

!------- q input to edmf_mynn (units: kg/kg) at full level -------
      if ( id_q_input > 0) then
        used = send_data (id_q_input, am4_Output_edmf%q_input, Time_next, is, js, 1 )
      endif

!------- qa input to edmf_mynn (units: none) at full level -------
      if ( id_qa_input > 0) then
        used = send_data (id_qa_input, am4_Output_edmf%qa_input, Time_next, is, js, 1 )
      endif

!------- ql input to edmf_mynn (units: kg/kg) at full level -------
      if ( id_ql_input > 0) then
        used = send_data (id_ql_input, am4_Output_edmf%ql_input, Time_next, is, js, 1 )
      endif

!------- qi input to edmf_mynn (units: kg/kg) at full level -------
      if ( id_qi_input > 0) then
        used = send_data (id_qi_input, am4_Output_edmf%qi_input, Time_next, is, js, 1 )
      endif

!------- thl input to edmf_mynn (units: K) at full level -------
      if ( id_thl_input > 0) then
        used = send_data (id_thl_input, am4_Output_edmf%thl_input, Time_next, is, js, 1 )
      endif

!------- qt input to edmf_mynn (units: kg/kg) at full level -------
      if ( id_qt_input > 0) then
        used = send_data (id_qt_input, am4_Output_edmf%qt_input, Time_next, is, js, 1 )
      endif

!------- rh input to edmf_mynn (units: percent) at full level -------
      if ( id_rh_input > 0) then
        used = send_data (id_rh_input, am4_Output_edmf%rh_input, Time_next, is, js, 1 )
      endif

!------- theta input to edmf_mynn (units: K) at full level -------
      if ( id_th_input > 0) then
        used = send_data (id_th_input, am4_Output_edmf%th_input, Time_next, is, js, 1 )
      endif

!------- T before_mix to edmf_mynn (units: K) at full level -------
      if ( id_t_before_mix > 0) then
        used = send_data (id_t_before_mix, am4_Output_edmf%t_before_mix, Time_next, is, js, 1 )
      endif

!------- q before_mix to edmf_mynn (units: kg/kg) at full level -------
      if ( id_q_before_mix > 0) then
        used = send_data (id_q_before_mix, am4_Output_edmf%q_before_mix, Time_next, is, js, 1 )
      endif

!------- qa before_mix to edmf_mynn (units: none) at full level -------
      if ( id_qa_before_mix > 0) then
        used = send_data (id_qa_before_mix, am4_Output_edmf%qa_before_mix, Time_next, is, js, 1 )
      endif

!------- ql before_mix to edmf_mynn (units: kg/kg) at full level -------
      if ( id_ql_before_mix > 0) then
        used = send_data (id_ql_before_mix, am4_Output_edmf%ql_before_mix, Time_next, is, js, 1 )
      endif

!------- qi before_mix to edmf_mynn (units: kg/kg) at full level -------
      if ( id_qi_before_mix > 0) then
        used = send_data (id_qi_before_mix, am4_Output_edmf%qi_before_mix, Time_next, is, js, 1 )
      endif

!------- thl before_mix to edmf_mynn (units: K) at full level -------
      if ( id_thl_before_mix > 0) then
        used = send_data (id_thl_before_mix, am4_Output_edmf%thl_before_mix, Time_next, is, js, 1 )
      endif

!------- qt before_mix to edmf_mynn (units: kg/kg) at full level -------
      if ( id_qt_before_mix > 0) then
        used = send_data (id_qt_before_mix, am4_Output_edmf%qt_before_mix, Time_next, is, js, 1 )
      endif

!------- rh before_mix to edmf_mynn (units: percent) at full level -------
      if ( id_rh_before_mix > 0) then
        used = send_data (id_rh_before_mix, am4_Output_edmf%rh_before_mix, Time_next, is, js, 1 )
      endif

!------- theta before_mix to edmf_mynn (units: K) at full level -------
      if ( id_th_before_mix > 0) then
        used = send_data (id_th_before_mix, am4_Output_edmf%th_before_mix, Time_next, is, js, 1 )
      endif

!------- T after_mix to edmf_mynn (units: K) at full level -------
      if ( id_t_after_mix > 0) then
        used = send_data (id_t_after_mix, am4_Output_edmf%t_after_mix, Time_next, is, js, 1 )
      endif

!------- q after_mix to edmf_mynn (units: kg/kg) at full level -------
      if ( id_q_after_mix > 0) then
        used = send_data (id_q_after_mix, am4_Output_edmf%q_after_mix, Time_next, is, js, 1 )
      endif

!------- qa after_mix to edmf_mynn (units: none) at full level -------
      if ( id_qa_after_mix > 0) then
        used = send_data (id_qa_after_mix, am4_Output_edmf%qa_after_mix, Time_next, is, js, 1 )
      endif

!------- ql after_mix to edmf_mynn (units: kg/kg) at full level -------
      if ( id_ql_after_mix > 0) then
        used = send_data (id_ql_after_mix, am4_Output_edmf%ql_after_mix, Time_next, is, js, 1 )
      endif

!------- qi after_mix to edmf_mynn (units: kg/kg) at full level -------
      if ( id_qi_after_mix > 0) then
        used = send_data (id_qi_after_mix, am4_Output_edmf%qi_after_mix, Time_next, is, js, 1 )
      endif

!------- thl after_mix to edmf_mynn (units: K) at full level -------
      if ( id_thl_after_mix > 0) then
        used = send_data (id_thl_after_mix, am4_Output_edmf%thl_after_mix, Time_next, is, js, 1 )
      endif

!------- qt after_mix to edmf_mynn (units: kg/kg) at full level -------
      if ( id_qt_after_mix > 0) then
        used = send_data (id_qt_after_mix, am4_Output_edmf%qt_after_mix, Time_next, is, js, 1 )
      endif

!------- rh after_mix to edmf_mynn (units: percent) at full level -------
      if ( id_rh_after_mix > 0) then
        used = send_data (id_rh_after_mix, am4_Output_edmf%rh_after_mix, Time_next, is, js, 1 )
      endif

!------- theta after_mix to edmf_mynn (units: K) at full level -------
      if ( id_th_after_mix > 0) then
        used = send_data (id_th_after_mix, am4_Output_edmf%th_after_mix, Time_next, is, js, 1 )
      endif

!------- qa before_pdf to edmf_mynn (units: none) at full level -------
      if ( id_qa_before_pdf > 0) then
        used = send_data (id_qa_before_pdf, am4_Output_edmf%qa_before_pdf, Time_next, is, js, 1 )
      endif

!------- ql before_pdf to edmf_mynn (units: kg/kg) at full level -------
      if ( id_ql_before_pdf > 0) then
        used = send_data (id_ql_before_pdf, am4_Output_edmf%ql_before_pdf, Time_next, is, js, 1 )
      endif

!------- qi before_pdf to edmf_mynn (units: kg/kg) at full level -------
      if ( id_qi_before_pdf > 0) then
        used = send_data (id_qi_before_pdf, am4_Output_edmf%qi_before_pdf, Time_next, is, js, 1 )
      endif

!------- moist updraft area on phalf (units: none) at half level -------
      if ( id_a_moist_half > 0) then
        used = send_data (id_a_moist_half, am4_Output_edmf%a_moist_half, Time_next, is, js, 1 )
      endif

!------- moist updraft area on pfull (units: none) at full level -------
      if ( id_a_moist_full > 0) then
        used = send_data (id_a_moist_full, am4_Output_edmf%a_moist_full, Time_next, is, js, 1 )
      endif

!------- moist updraft mass flux on phalf (units: kg/m2/s) at half level -------
      if ( id_mf_moist_half > 0) then
        used = send_data (id_mf_moist_half, am4_Output_edmf%mf_moist_half, Time_next, is, js, 1 )
      endif

!------- moist updraft mass flux on pfull (units: kg/m2/s) at full level -------
      if ( id_mf_moist_full > 0) then
        used = send_data (id_mf_moist_full, am4_Output_edmf%mf_moist_full, Time_next, is, js, 1 )
      endif

!------- spec humid of moist updraft on phalf (units: kg/kg) at half level -------
      if ( id_qv_moist_half > 0) then
        used = send_data (id_qv_moist_half, am4_Output_edmf%qv_moist_half, Time_next, is, js, 1 )
      endif

!------- spec humid of moist updraft on pfull (units: kg/kg) at full level -------
      if ( id_qv_moist_full > 0) then
        used = send_data (id_qv_moist_full, am4_Output_edmf%qv_moist_full, Time_next, is, js, 1 )
      endif

!------- dry updraft area on phalf (units: none) at half level -------
      if ( id_a_dry_half > 0) then
        used = send_data (id_a_dry_half, am4_Output_edmf%a_dry_half, Time_next, is, js, 1 )
      endif

!------- dry updraft area on pfull (units: none) at full level -------
      if ( id_a_dry_full > 0) then
        used = send_data (id_a_dry_full, am4_Output_edmf%a_dry_full, Time_next, is, js, 1 )
      endif

!------- dry updraft mass flux on phalf (units: kg/m2/s) at half level -------
      if ( id_mf_dry_half > 0) then
        used = send_data (id_mf_dry_half, am4_Output_edmf%mf_dry_half, Time_next, is, js, 1 )
      endif

!------- dry updraft mass flux on pfull (units: kg/m2/s) at full level -------
      if ( id_mf_dry_full > 0) then
        used = send_data (id_mf_dry_full, am4_Output_edmf%mf_dry_full, Time_next, is, js, 1 )
      endif

!------- spec humid of dry updraft on phalf (units: kg/kg) at half level -------
      if ( id_qv_dry_half > 0) then
        used = send_data (id_qv_dry_half, am4_Output_edmf%qv_dry_half, Time_next, is, js, 1 )
      endif

!------- spec humid of dry updraft on pfull (units: kg/kg) at full level -------
      if ( id_qv_dry_full > 0) then
        used = send_data (id_qv_dry_full, am4_Output_edmf%qv_dry_full, Time_next, is, js, 1 )
      endif

!------- updraft mass flux on phalf (units: kg/m2/s) at half level -------
      if ( id_mf_all_half > 0) then
        used = send_data (id_mf_all_half, am4_Output_edmf%mf_all_half, Time_next, is, js, 1 )
      endif

!------- updraft mass flux on pfull (units: kg/m2/s) at full level -------
      if ( id_mf_all_full > 0) then
        used = send_data (id_mf_all_full, am4_Output_edmf%mf_all_full, Time_next, is, js, 1 )
      endif

!------- number of edmf updrafts (units: none) at half level -------
      if ( id_num_updraft > 0) then
        call reshape_mynn_array_to_am4_half(ix, jx, kx, num_updraft, diag_half)
        used = send_data (id_num_updraft, diag_half, Time_next, is, js, 1 )
      endif

!------- number of dentrainment in edmf updrafts (units: none) at full level -------
      if ( id_num_det > 0) then
        call reshape_mynn_array_to_am4(ix, jx, kx, num_DET, diag_full)
        used = send_data (id_num_det, diag_full, Time_next, is, js, 1 )
      endif

!------- number of negative dentrainment with zero entrainment in edmf updrafts (units: none) at full level -------
      if ( id_num_ndet_zent > 0) then
        call reshape_mynn_array_to_am4(ix, jx, kx, num_nDET_zENT, diag_full)
        used = send_data (id_num_ndet_zent, diag_full, Time_next, is, js, 1 )
      endif

!------- number of negative dentrainment with positive entrainment in edmf updrafts (units: none) at full level -------
      if ( id_num_ndet_pent > 0) then
        call reshape_mynn_array_to_am4(ix, jx, kx, num_nDET_pENT, diag_full)
        used = send_data (id_num_ndet_pent, diag_full, Time_next, is, js, 1 )
      endif

!------- ratio for realizability limiter (units: none) at one level -------
      if ( id_rlz_ratio > 0) then
        used = send_data (id_rlz_ratio, rlz_ratio, Time_next, is, js )
      endif

!------- tracer name for realizability limiter (units: none) at one level -------
      if ( id_rlz_tracer > 0) then
        used = send_data (id_rlz_tracer, rlz_tracer, Time_next, is, js )
      endif
!send_data

!----------
! debug
!----------
  
  !--- check whether edmf_a>1
  call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_a, diag_half)
  if (do_debug_option.eq."edmf_a" .or. do_debug_option.eq."all") then
    do i=1,ix
    do j=1,jx
    do k=1,kx+1
      if (diag_half(i,j,k).gt.1.) then
        print*,'gg01, i,j,k,lon,lat,edmf_a', i,j,k,lon(i,j),lat(i,j),diag_half(i,j,k)
      endif
    enddo
    enddo
    enddo
  endif

  !--- check whether tdt is too large
  if (do_debug_option.eq."tdt_check" .or. do_debug_option.eq."all") then
    tt1 = tdt_max / 86400. ! change unit from K/day to K/s
    do i=1,ix
    do j=1,jx
    do k=1,kx
      if (abs(am4_Output_edmf%tdt_edmf(i,j,k)).gt.tt1) then
        print*,'gg02, i,j,k,lon,lat,tdt_edmf', i,j,k,lon(i,j),lat(i,j),am4_Output_edmf%tdt_edmf(i,j,k)
      endif
    enddo
    enddo
    enddo
  end if

  !--- check whether tracers become negative
  if (do_debug_option.eq."check_rlz" .or. do_debug_option.eq."all") then
    print*,'new qv',Physics_input_block%q(:,:,:,nsphum) + dt * am4_Output_edmf%qdt_edmf(:,:,:)
    print*,'new ql',Physics_input_block%q(:,:,:,nql)    + dt * am4_Output_edmf%qldt_edmf(:,:,:)
    print*,'new qi',Physics_input_block%q(:,:,:,nqi)    + dt * am4_Output_edmf%qidt_edmf(:,:,:)
    print*,'new qa',Physics_input_block%q(:,:,:,nqa)    + dt * am4_Output_edmf%qadt_edmf(:,:,:)
  endif

!---------------------------------------------------------------------
! deallocate EDMF-MYNN input and output variables 
!---------------------------------------------------------------------

  call edmf_dealloc (Input_edmf, Output_edmf, am4_Output_edmf)

  !--- stop the model if prefered
  if (do_stop_run) then
    call error_mesg('edmf_mynn_driver',  &
       'stop by yihsuan', FATAL)
  endif


end subroutine edmf_mynn_driver

!###########################################################

subroutine edmf_alloc ( &
              is, ie, js, je, npz, Time_next, dt, frac_land, area, u_star,  &
              b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, rdiag, udt, vdt, tdt, rdt, &
              Qke, el_pbl, cldfra_bl, qc_bl, Sh3D, &
              Input_edmf, Output_edmf, am4_Output_edmf )

!---------------------------------------------------------------------
! Arguments (Intent in)  
!    is,ie               -  starting and ending i indices for window
!    js,je               -  starting and ending j indices for window
!    npz                 -  number of model layers
!    Time_next           -  variable needed for netcdf output diagnostics
!    dt                  -  Time step               (sec)
!    frac_land           -  fraction of surface covered by land, dimension (nlon, nlat)
!    area                -  grid box area           (m^2)   , dimension (nlon, nlat) 
!    u_star              -  friction velocity       (m/s)   , dimension (nlon, nlat)        ,see note 1 below
!    b_star              -  buoyancy scale          (m/s^2) , dimension (nlon, nlat)        ,see note 1 below
!    q_star              -  moisture scale          (kg/kg) , dimension (nlon, nlat)        ,see note 1 below
!    shflx               -  sensible heat flux      (W/m^2) , dimension (nlon, nlat)
!    lhflx               -  surface moisture flux   (kg/m^2/s) , dimension (nlon, nlat)
!    Physics_input_block -  derived type variable containing
!                             1) p_half         pressure at half levels (offset from t,q,u,v,r)
!                                              [ Pa ]
!                             2) p_full         pressure at full levels [ Pa }
!                             3) z_half         height at half levels [ m ]
!                             4) z_full         height at full levels [ m ]
!                             5) u              zonal wind at current time step [ m / s ]
!                             6) v              meridional wind at current time step [ m / s ]
!                             7) t              temperature at current time step [ deg k ]
!                             9) q              multiple 3d tracer fields at current time step
!                            10) ...
!
!    rdiag               - multiple 3d diagnostic tracer fields [ unit / unit ]
!
! Note:
!    1. u_star and b_star are from monin_obukhov program, src/atmos_param/monin_obukhov/monin_obukhov_kernel.F90
!       q_star            is  from surface flux  program, src/FMScoupler/surface_flux.F90    
!
!            The magnitude of the wind stress is 
!                 density*(u_star**2)
!            The buoyancy flux, (g/theta_v) * <w'theta_v'>|surface, is
!                 u_star*b_star
!            The evaporation rate), kg vapor/m^2/s, is
!                 density*u_star*q_star
!                 So, u_star*q_star is surface moisture flux, <w'q'>, unit: (m/s * kg vapor/kg air)
!
!    2. The u,v,t,q are at full levels. The vertical indexing is counted downward, i.e. k=1 is at the top of the model.
!
!         --------- 1   (top of the atmospheric model)
!           * 1    
!         --------- 2
!           * 2          -->  full levels where u,v,t,q are located
!         --------- k-1
!           * k-1
!         --------- k    
!           * k               
!         --------- k+1
!           * k+1
!         --------- ...
!           ....
!         --------- kx
!           * kx
!         --------- kxp=kx+1, surface
!
!    3. Note that diff_t and diff_m are no half levels, but GFDL physics_driver use diff_t(nlay) instead of nlay+1.
!       So, although diff_t_edmf and diff_m_edmf vertical dimension is "nlay", but they are located on half levels. 
!---------------------------------------------------------------------
  integer, intent(in)                   :: is, ie, js, je, npz
  type(time_type), intent(in)           :: Time_next
  real,    intent(in)                   :: dt
  real,    intent(in), dimension(:,:)   :: &
    frac_land, area, u_star, b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux
  type(physics_input_block_type)        :: Physics_input_block
  real, intent(in), dimension(:,:,:,:)  :: rdiag
  real, intent(in), dimension(:,:,:)    :: &
    Qke, el_pbl, cldfra_bl, qc_bl, Sh3D
  real, intent(in), dimension(:,:,:)    :: &  ! (nlon,nlat,nlay)
    udt, vdt, tdt  
  real, intent(in), dimension(:,:,:,:)  :: &  ! (nlon,nlat,nlay,ntracers)
    rdt

!---------------------------------------------------------------------
! Arguments (Intent out)
!
!***  vertical indexing
!   KMS (smallest number) is the bottom level and KME (largest
! number) is the top level.  
! ---------
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
! ----------
!
!*****************
!*****************
!  Input_edmf%
!*****************
!*****************
!
!--- INPUTS/0D fields
!
! REAL:: 
! delt ... model time step [s]
! dx ...  horizontal length of grid [m] ... this is used in the scale aware option
!
!--- INPUTS/2D surface fields
!
!  REAL, DIMENSION(IMS:IME,JMS:JME) :: 
! xland ....   land mask (1 for land, 2 for water)
! ust ...  u* in similarity theory (m/s)
! ch .... drag coefficient for heat/moisture       
! rmol ....   inverse Monin-Obukhov length (1/m)
! ts  ... surface temperature (K) 
! qsfc ...   specific humidity at lower boundary (kg/kg)
! qcg   ... saturated mixing ratio at the surface (kg/kg)
! ps .... surface pressure [Pa]
! hfx ... upward heat flux at the surface (W/m^2)
! qfx  ... upward moisture flux at the surface (kg/m^2/s)
! wspd ... wind speed at lowest model level (m/s)
! uoce  ...  sea surface zonal ocean currents (m s-1)
! voce ... sea surface meridional ocean currents (m s-1)
! vdfg  ... deposition velocity of fog (m/s); set to zero
! znt ...   thermal roughness length (m)
!
!--- INPUTS/3D thermodynamic fields on half levels
!
!   REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(in) :: 
! dz ... dz between full levels [m], on full levels
! u  ... zonal wind speed [m s^-1], on half levels
! v ... meridional wind speed [m s^-1], on half levels
! w ... vertical wind speed [m s^-1], on full levels
! th  ... potential temperature [K], on half levels
! qv .... water vapor mixing ratio [kg kg^-1], half level
! p .... pressure [Pa], half levels
! exner .... exner function, half levels 
! rho ... density [kg m^-3], half levels
! T3D ... temperature [K], half levels
! ql ... cloud liquid water mixing ratio [kg kg^-1], half levels 
! qi ... cloud ice water mixing ratio [kg kg^-1], half levels 
! qnc,qni ....  cloud liq and ice number concentration (#/kg), optional, do not set them
!
!*****************
!*****************
!  Output_edmf
!*****************
!*****************
!
!--- INPUT-OUTPUT/2D fields
!
!   DIMENSION( ims:ime, kms:kme, jms:jme ),OPTIONAL  ::pattern_spp_pbl
!         (pattern for stochastic condensation .... not needed if spp_pbl == 0) 
!
!    REAL, DIMENSION(IMS:IME,JMS:JME) 
!         pblh,wstar,delta,tke_pbl 
! (these are essentially output variables, for coupling with one of the convection scheme, you can treat them as dummy) 
!
!
!--- INPUT-OUTPUT/3D fields
!
!    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(inout) :: &
!         Qke .... twice the turbulent kinetic energy [m^2 s^-2]
!         Tsq  ...variance of theta_l (K^2)
!         Qsq ... variance of q_t (-)
!         Cov ... covariance of theta_l and q_t  (K)
!         qke_adv   .... Twice the turbulent kinetic energy after advected with resolved scale flow (advection is done outside the routine) [m^2 s^-2]
!
!
!--- INPUT-OUTPUT/3D fields (tendencies, should be set to 0 on input; on full levels)
!    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: &
!         RUBLTEN ... U tendency due to EDMF parameterization (m/s^2), set to 0. on the input
!         RVBLTEN ... V tendency due to EDMF parameterization (m/s^2), set to 0. on the input
!         RTHBLTEN ... Theta tendency due to EDMF parameterization (K/s), set to 0. on the input
!         RQVBLTEN ....  Qv tendency due to EDMF parameterization (kg/kg/s), set to 0. on the input
!         RQLBLTEN ....  Qc tendency due to EDMF parameterization (kg/kg/s), set to 0. on the input
!         RQIBLTEN ...   Qi tendency due to  EDMF parameterization (kg/kg/s), set to 0. on the input
!         RQNIBLTEN ....  Qni tendency due to EDMF parameterization (#/kg/s) , dummy argument (not activated)
!         RTHRATEN .... tendency from radiation [K s^-1]
!
!--- INPUT/OUTPUT, OPTIONAL/3D fields (updraft/downdraft properties) on full levels
!   REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME):: 
!         edmf_a ... updraft area [-] 
!         edmf_w ... vertical velocity of updrafts [m s^-1]
!         edmf_qt .... qt in updrafts  [kg kg^-1]
!         edmf_thl .... thl in updrafts  [K]
!         edmf_ent ... entrainment in updrafts [m^-1]
!         edmf_qc ... qc in updrafts [kg kg^-1]
!         edmf_a_dd,edmf_w_dd,edmf_qt_dd,edmf_thl_dd,edmf_ent_dd,edmf_qc_dd ... dummy outputs (downdraft poperties)
!         edmf_debug1,edmf_debug2,edmf_debug3,edmf_debug4 .... dummy outputs (for debugging)
!
!
!--- INPUT/OUTPUT, 3D cloud and turbulence fields
!
! REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME)
!         mynn_ql & qc_bl ... liquid water mixing ratio from PDF cloud scheme [kg kg-1]
!         cldfra_bl ... cloud fraction [-]
!         el_pbl ..... turbulent mixing length [m]
!
!
!---  INPUT/OUTPUT 
!    INTEGER,DIMENSION(IMS:IME,JMS:JME) 
!       KPBL .... model level of boundary layer height
!       nupdraft ... diagnostic output, number of updrafts (this is only for stem MF and not for JPLedmf)
!     ktop_shallow ... diagnostic output essentially, the level of the highest updraft
!
!--- OUTPUTS, 3D mixing coefficients on full levels
!    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: 
!         exch_h ... ED mixing coefficient for heat, K m/s 
!         exch_m .... ED mixing coefficient for momentum, m^2/s
!
!--- OUTPUT, 2D
! REAL, DIMENSION(IMS:IME,JMS:JME) ::
!        maxmf   .... maximum mass-flux (does not work with the JPL-EDMF, set to dummy variable)
!---------------------------------------------------------------------

  type(edmf_input_type) , intent(out) :: Input_edmf
  type(edmf_output_type), intent(out) :: Output_edmf
  type(am4_edmf_output_type), intent(out) :: am4_Output_edmf

!---------------------------------------------------------------------
! local variables  
!---------------------------------------------------------------------
  integer :: ix, jx, kx, nt, kxp
  integer :: &
    IDS,IDE,JDS,JDE,KDS,KDE,        &
    IMS,IME,JMS,JME,KMS,KME,        &
    ITS,ITE,JTS,JTE,KTS,KTE

  real, dimension (size(Physics_input_block%t,1), &
                   size(Physics_input_block%t,2), &
                   size(Physics_input_block%t,3)) :: &
    dz_host, u_host, v_host, w_host, th_host, qv_host, p_host, exner_host, rho_host, T3D_host, qa_host, qc_host, ql_host, qi_host, qnc_host, qni_host, rh_host, &
    thl_host, qt_host, &
    tv_host, thv_host, omega_host

  real, dimension (size(Physics_input_block%t,1), &
                   size(Physics_input_block%t,2)) :: &
    u_star_host, shflx_host, lhflx_host, w1_thv1_surf_host, w1_th1_surf_host, w1_q1_surf_host, Obukhov_length_host, &
    u_star_star, shflx_star, lhflx_star, w1_thv1_surf_star, w1_th1_surf_star, w1_t1_surf_star, w1_q1_surf_star, Obukhov_length_star,  &
    u_star_updated, shflx_updated, lhflx_updated, w1_thv1_surf_updated, w1_th1_surf_updated, w1_t1_surf_updated, w1_q1_surf_updated, Obukhov_length_updated, &
    tv_ref, thv_ref, rho_ref

  integer :: i, j, k, kk

!-------------------------------------------------------------------------
!  define input array sizes.
!-------------------------------------------------------------------------
  ix = size(Physics_input_block%t,1)
  jx = size(Physics_input_block%t,2)
  kx = size(Physics_input_block%t,3)
  nt = size(Physics_input_block%q,4)
  kxp = kx + 1

  ! 0 means unused
  !IMS =  1  ; IME =  ix ; JMS =  1  ; JME =  jx ; KMS  = 1  ; KME = kx
  !IDS =  0  ; IDE =   0 ; JDS =  0  ; JDE =   0 ; KDS =  0  ; KDE =  0
  !ITS =  0  ; ITE =   0 ; JTS =  0  ; JTE =   0 ; KTS =  1  ; KTE = kx
  IMS =  1  ; IME =  ix  ; JMS =  1  ; JME =  jx ; KMS  = 1  ; KME = kx
  IDS =  1  ; IDE =  ix  ; JDS =  1  ; JDE =  jx ; KDS =  1  ; KDE = kx
  ITS =  1  ; ITE =  ix  ; JTS =  1  ; JTE =  jx ; KTS =  1  ; KTE = kx

!-------------------------------------------------------------------------
! allocate and initialize variables 
!-------------------------------------------------------------------------

!******************
!--- Input_edmf
!******************
 
  ! diagnostic purpose
  allocate (Input_edmf%u_star_star(IMS:IME,JMS:JME))          ; Input_edmf%u_star_star = 0.
  allocate (Input_edmf%shflx_star(IMS:IME,JMS:JME))           ; Input_edmf%shflx_star = 0.
  allocate (Input_edmf%lhflx_star(IMS:IME,JMS:JME))           ; Input_edmf%lhflx_star = 0.
  allocate (Input_edmf%w1_thv1_surf_star(IMS:IME,JMS:JME))    ; Input_edmf%w1_thv1_surf_star = 0.
  allocate (Input_edmf%w1_th1_surf_star(IMS:IME,JMS:JME))     ; Input_edmf%w1_th1_surf_star = 0.
  allocate (Input_edmf%w1_q1_surf_star(IMS:IME,JMS:JME))      ; Input_edmf%w1_q1_surf_star = 0.
  allocate (Input_edmf%Obukhov_length_star(IMS:IME,JMS:JME))  ; Input_edmf%Obukhov_length_star = 0.

  allocate (Input_edmf%u_star_updated(IMS:IME,JMS:JME))          ; Input_edmf%u_star_updated = 0.
  allocate (Input_edmf%shflx_updated(IMS:IME,JMS:JME))           ; Input_edmf%shflx_updated = 0.
  allocate (Input_edmf%lhflx_updated(IMS:IME,JMS:JME))           ; Input_edmf%lhflx_updated = 0.
  allocate (Input_edmf%w1_thv1_surf_updated(IMS:IME,JMS:JME))    ; Input_edmf%w1_thv1_surf_updated = 0.
  allocate (Input_edmf%w1_th1_surf_updated(IMS:IME,JMS:JME))     ; Input_edmf%w1_th1_surf_updated = 0.
  allocate (Input_edmf%w1_q1_surf_updated(IMS:IME,JMS:JME))      ; Input_edmf%w1_q1_surf_updated = 0.
  allocate (Input_edmf%Obukhov_length_updated(IMS:IME,JMS:JME))  ; Input_edmf%Obukhov_length_updated = 0.
  !allocate (Input_edmf%(IMS:IME,JMS:JME))  ; Input_edmf% = 0.
 
  ! index variable
  allocate (Input_edmf%IMS) ; Input_edmf%IMS = IMS
  allocate (Input_edmf%IME) ; Input_edmf%IME = IME
  allocate (Input_edmf%JMS) ; Input_edmf%JMS = JMS
  allocate (Input_edmf%JME) ; Input_edmf%JME = JME
  allocate (Input_edmf%KMS) ; Input_edmf%KMS = KMS
  allocate (Input_edmf%KME) ; Input_edmf%KME = KME

  allocate (Input_edmf%IDS) ; Input_edmf%IDS = IDS
  allocate (Input_edmf%IDE) ; Input_edmf%IDE = IDE
  allocate (Input_edmf%JDS) ; Input_edmf%JDS = JDS
  allocate (Input_edmf%JDE) ; Input_edmf%JDE = JDE
  allocate (Input_edmf%KDS) ; Input_edmf%KDS = KDS
  allocate (Input_edmf%KDE) ; Input_edmf%KDE = KDE

  allocate (Input_edmf%ITS) ; Input_edmf%ITS = ITS
  allocate (Input_edmf%ITE) ; Input_edmf%ITE = ITE
  allocate (Input_edmf%JTS) ; Input_edmf%JTS = JTS
  allocate (Input_edmf%JTE) ; Input_edmf%JTE = JTE
  allocate (Input_edmf%KTS) ; Input_edmf%KTS = KTS
  allocate (Input_edmf%KTE) ; Input_edmf%KTE = KTE

  ! 0-D variable
  allocate (Input_edmf%delt) ; Input_edmf%delt = 0. 
  allocate (Input_edmf%dx)   ; Input_edmf%dx   = 0. 
  !allocate (Input_edmf%)   ; Input_edmf% = 

  ! 2-D variable
  allocate (Input_edmf%xland (IMS:IME,JMS:JME))  ; Input_edmf%xland = 0.
  allocate (Input_edmf%ust   (IMS:IME,JMS:JME))  ; Input_edmf%ust   = 0.
  allocate (Input_edmf%ch    (IMS:IME,JMS:JME))  ; Input_edmf%ch    = 0.
  allocate (Input_edmf%rmol  (IMS:IME,JMS:JME))  ; Input_edmf%rmol  = 0.
  allocate (Input_edmf%ts    (IMS:IME,JMS:JME))  ; Input_edmf%ts    = 0.
  allocate (Input_edmf%qsfc  (IMS:IME,JMS:JME))  ; Input_edmf%qsfc  = 0.
  allocate (Input_edmf%qcg   (IMS:IME,JMS:JME))  ; Input_edmf%qcg   = 0.
  allocate (Input_edmf%ps    (IMS:IME,JMS:JME))  ; Input_edmf%ps    = 0.
  allocate (Input_edmf%hfx   (IMS:IME,JMS:JME))  ; Input_edmf%hfx   = 0.
  allocate (Input_edmf%qfx   (IMS:IME,JMS:JME))  ; Input_edmf%qfx   = 0.
  allocate (Input_edmf%wspd  (IMS:IME,JMS:JME))  ; Input_edmf%wspd  = 0.
  allocate (Input_edmf%uoce  (IMS:IME,JMS:JME))  ; Input_edmf%uoce  = 0.
  allocate (Input_edmf%voce  (IMS:IME,JMS:JME))  ; Input_edmf%voce  = 0.
  allocate (Input_edmf%vdfg  (IMS:IME,JMS:JME))  ; Input_edmf%vdfg  = 0.
  allocate (Input_edmf%znt   (IMS:IME,JMS:JME))  ; Input_edmf%znt   = 0.
  !allocate (Input_edmf%(IMS:IME,JMS:JME))  ; Input_edmf% = 0.

  ! 3-D variable
  allocate (Input_edmf%dz    (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%dz    = 0.
  allocate (Input_edmf%u     (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%u     = 0.
  allocate (Input_edmf%v     (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%v     = 0.
  allocate (Input_edmf%w     (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%w     = 0.
  allocate (Input_edmf%th    (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%th    = 0.
  allocate (Input_edmf%qv    (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%qv    = 0.
  allocate (Input_edmf%p     (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%p     = 0.
  allocate (Input_edmf%exner (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%exner = 0.
  allocate (Input_edmf%rho   (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%rho   = 0.
  allocate (Input_edmf%T3D   (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%T3D   = 0.
  allocate (Input_edmf%qc    (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%qc    = 0.
  allocate (Input_edmf%ql    (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%ql    = 0.
  allocate (Input_edmf%qi    (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%qi    = 0.
  allocate (Input_edmf%qa    (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%qa    = 0.
  allocate (Input_edmf%qnc   (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%qnc   = 0.
  allocate (Input_edmf%qni   (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%qni   = 0.
  !allocate (Input_edmf%(IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf% = 0.

  ! semi-prognostic variables
  allocate (Input_edmf%Qke         (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%Qke         = 0.
  allocate (Input_edmf%Sh3D        (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%Sh3D        = 0.
  allocate (Input_edmf%el_pbl      (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%el_pbl      = 0.
  allocate (Input_edmf%cldfra_bl   (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%cldfra_bl   = 0.
  allocate (Input_edmf%qc_bl       (IMS:IME,KMS:KME,JMS:JME))  ; Input_edmf%qc_bl       = 0.


!******************
!--- Output_edmf
!******************

  ! 2-D variable
  allocate (Output_edmf%pblh         (IMS:IME,JMS:JME))  ; Output_edmf%pblh         = 0.
  allocate (Output_edmf%KPBL         (IMS:IME,JMS:JME))  ; Output_edmf%KPBL         = 0
  allocate (Output_edmf%nupdraft     (IMS:IME,JMS:JME))  ; Output_edmf%nupdraft     = 0
  allocate (Output_edmf%ktop_shallow (IMS:IME,JMS:JME))  ; Output_edmf%ktop_shallow = 0
  allocate (Output_edmf%maxmf        (IMS:IME,JMS:JME))  ; Output_edmf%maxmf        = 0.
  !allocate (Output_edmf%(IMS:IME,JMS:JME))  ; Output_edmf% = 0.

  ! 3-D variable
  allocate (Output_edmf%Qke         (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Qke         = 0.
!print*,'alloc, Output_edmf%Qke, ix,jx,kx',size(Output_edmf%Qke,1),size(Output_edmf%Qke,2),size(Output_edmf%Qke,3)
  allocate (Output_edmf%Tsq         (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Tsq         = 0.
  allocate (Output_edmf%Qsq         (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Qsq         = 0.
  allocate (Output_edmf%Cov         (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Cov         = 0.
  allocate (Output_edmf%RUBLTEN     (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RUBLTEN     = 0.
  allocate (Output_edmf%RVBLTEN     (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RVBLTEN     = 0.
  allocate (Output_edmf%RTHBLTEN    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RTHBLTEN    = 0.
  allocate (Output_edmf%RQVBLTEN    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RQVBLTEN    = 0.
  allocate (Output_edmf%RQLBLTEN    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RQLBLTEN    = 0.
  allocate (Output_edmf%RQIBLTEN    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RQIBLTEN    = 0.
  allocate (Output_edmf%RQNIBLTEN   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RQNIBLTEN   = 0.
  allocate (Output_edmf%RTHRATEN    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RTHRATEN    = 0.
  allocate (Output_edmf%RCCBLTEN    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RCCBLTEN    = 0.
  allocate (Output_edmf%RTHLBLTEN   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RTHLBLTEN   = 0.
  allocate (Output_edmf%RQTBLTEN    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RQTBLTEN    = 0.
  allocate (Output_edmf%edmf_a      (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_a      = 0.
  allocate (Output_edmf%edmf_w      (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_w      = 0.
  allocate (Output_edmf%edmf_qt     (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_qt     = 0.
  allocate (Output_edmf%edmf_thl    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_thl    = 0.
  allocate (Output_edmf%edmf_ent    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_ent    = 0.
  allocate (Output_edmf%edmf_det    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_det    = 0.
  allocate (Output_edmf%edmf_qc     (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_qc     = 0.
  allocate (Output_edmf%edmf_a_dd   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_a_dd   = 0.
  allocate (Output_edmf%edmf_w_dd   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_w_dd   = 0.
  allocate (Output_edmf%edmf_qt_dd  (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_qt_dd  = 0.
  allocate (Output_edmf%edmf_thl_dd (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_thl_dd = 0.
  allocate (Output_edmf%edmf_ent_dd (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_ent_dd = 0.
  allocate (Output_edmf%edmf_qc_dd  (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_qc_dd  = 0.
  allocate (Output_edmf%edmf_debug1 (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_debug1 = 0.
  allocate (Output_edmf%edmf_debug2 (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_debug2 = 0.
  allocate (Output_edmf%edmf_debug3 (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_debug3 = 0.
  allocate (Output_edmf%edmf_debug4 (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_debug4 = 0.
  allocate (Output_edmf%mynn_ql     (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%mynn_ql     = 0.
  allocate (Output_edmf%qc_bl       (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%qc_bl       = 0.
  allocate (Output_edmf%cldfra_bl   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%cldfra_bl   = 0.
  allocate (Output_edmf%el_pbl      (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%el_pbl      = 0.
  allocate (Output_edmf%Sh3D        (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Sh3D        = 0.
  allocate (Output_edmf%exch_h      (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%exch_h      = 0.
  allocate (Output_edmf%exch_m      (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%exch_m = 0.
  allocate (Output_edmf%qWT         (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%qWT = 0.
  allocate (Output_edmf%qSHEAR      (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%qSHEAR = 0.
  allocate (Output_edmf%qBUOY       (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%qBUOY = 0.
  allocate (Output_edmf%qDISS       (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%qDISS = 0.
  allocate (Output_edmf%dqke        (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%dqke = 0.

  allocate (Output_edmf%t_before_mix    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%t_before_mix    = 0.
  allocate (Output_edmf%q_before_mix    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%q_before_mix    = 0.
  allocate (Output_edmf%qa_before_mix   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%qa_before_mix   = 0.
  allocate (Output_edmf%ql_before_mix   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%ql_before_mix   = 0.
  allocate (Output_edmf%qi_before_mix   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%qi_before_mix   = 0.
  allocate (Output_edmf%thl_before_mix  (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%thl_before_mix  = 0.
  allocate (Output_edmf%qt_before_mix   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%qt_before_mix   = 0.
  allocate (Output_edmf%rh_before_mix   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%rh_before_mix   = 0.
  allocate (Output_edmf%th_before_mix   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%th_before_mix   = 0.
  allocate (Output_edmf%t_after_mix     (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%t_after_mix     = 0.
  allocate (Output_edmf%q_after_mix     (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%q_after_mix     = 0.
  allocate (Output_edmf%qa_after_mix    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%qa_after_mix    = 0.
  allocate (Output_edmf%ql_after_mix    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%ql_after_mix    = 0.
  allocate (Output_edmf%qi_after_mix    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%qi_after_mix    = 0.
  allocate (Output_edmf%thl_after_mix   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%thl_after_mix   = 0.
  allocate (Output_edmf%qt_after_mix    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%qt_after_mix    = 0.
  allocate (Output_edmf%rh_after_mix    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%rh_after_mix    = 0.
  allocate (Output_edmf%th_after_mix    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%th_after_mix    = 0.
  allocate (Output_edmf%qa_before_pdf   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%qa_before_pdf   = 0.
  allocate (Output_edmf%ql_before_pdf   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%ql_before_pdf   = 0.
  allocate (Output_edmf%qi_before_pdf   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%qi_before_pdf   = 0.
  allocate (Output_edmf%Q_ql            (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_ql = 0.
  allocate (Output_edmf%Q_qi            (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_qi = 0.
  allocate (Output_edmf%Q_qa            (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_qa = 0.
  allocate (Output_edmf%Q_ql_adv        (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_ql_adv = 0.
  allocate (Output_edmf%Q_qi_adv        (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_qi_adv = 0.
  allocate (Output_edmf%Q_qa_adv        (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_qa_adv = 0.
  allocate (Output_edmf%Q_ql_eddy       (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_ql_eddy = 0.
  allocate (Output_edmf%Q_qi_eddy       (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_qi_eddy = 0.
  allocate (Output_edmf%Q_qa_eddy       (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_qa_eddy = 0.
  allocate (Output_edmf%Q_ql_ent        (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_ql_ent = 0.
  allocate (Output_edmf%Q_qi_ent        (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_qi_ent = 0.
  allocate (Output_edmf%Q_qa_ent        (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_qa_ent = 0.
  allocate (Output_edmf%Q_ql_det        (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_ql_det = 0.
  allocate (Output_edmf%Q_qi_det        (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_qi_det = 0.
  allocate (Output_edmf%Q_qa_det        (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_qa_det = 0.
  allocate (Output_edmf%Q_ql_sub        (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_ql_sub = 0.
  allocate (Output_edmf%Q_qi_sub        (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_qi_sub = 0.
  allocate (Output_edmf%Q_qa_sub        (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%Q_qa_sub = 0.

  allocate (Output_edmf%a_moist_half    (IMS:IME,KMS:KME+1,JMS:JME))  ; Output_edmf%a_moist_half  = 0.
  allocate (Output_edmf%mf_moist_half   (IMS:IME,KMS:KME+1,JMS:JME))  ; Output_edmf%mf_moist_half = 0.
  allocate (Output_edmf%qv_moist_half   (IMS:IME,KMS:KME+1,JMS:JME))  ; Output_edmf%qv_moist_half = 0.
  allocate (Output_edmf%a_moist_full    (IMS:IME,KMS:KME  ,JMS:JME))  ; Output_edmf%a_moist_full  = 0.
  allocate (Output_edmf%mf_moist_full   (IMS:IME,KMS:KME  ,JMS:JME))  ; Output_edmf%mf_moist_full = 0.
  allocate (Output_edmf%qv_moist_full   (IMS:IME,KMS:KME  ,JMS:JME))  ; Output_edmf%qv_moist_full = 0.

  allocate (Output_edmf%a_dry_half      (IMS:IME,KMS:KME+1,JMS:JME))  ; Output_edmf%a_dry_half    = 0.
  allocate (Output_edmf%a_dry_full      (IMS:IME,KMS:KME,  JMS:JME))  ; Output_edmf%a_dry_full    = 0.
  allocate (Output_edmf%mf_dry_half     (IMS:IME,KMS:KME+1,JMS:JME))  ; Output_edmf%mf_dry_half   = 0.
  allocate (Output_edmf%mf_dry_full     (IMS:IME,KMS:KME  ,JMS:JME))  ; Output_edmf%mf_dry_full   = 0.
  allocate (Output_edmf%qv_dry_half     (IMS:IME,KMS:KME+1,JMS:JME))  ; Output_edmf%qv_dry_half   = 0.
  allocate (Output_edmf%qv_dry_full     (IMS:IME,KMS:KME  ,JMS:JME))  ; Output_edmf%qv_dry_full   = 0.

  allocate (Output_edmf%mf_all_half     (IMS:IME,KMS:KME+1,JMS:JME))  ; Output_edmf%mf_all_half   = 0.
  allocate (Output_edmf%mf_all_full     (IMS:IME,KMS:KME  ,JMS:JME))  ; Output_edmf%mf_all_full   = 0.
  !allocate (Output_edmf%(IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf% = 0.

!**********************
!--- am4 Output_edmf
!**********************

  !print*,'alloc am4_Output_edmf, ix,jx,kx',ix,jx,kx
  allocate (am4_Output_edmf%tke         (ix,jx,kx))  ; am4_Output_edmf%tke         = 0.
  allocate (am4_Output_edmf%Tsq         (ix,jx,kx))  ; am4_Output_edmf%Tsq         = 0.
  allocate (am4_Output_edmf%Cov_thl_qt  (ix,jx,kx))  ; am4_Output_edmf%Cov_thl_qt  = 0.
  allocate (am4_Output_edmf%udt_edmf    (ix,jx,kx))  ; am4_Output_edmf%udt_edmf    = 0.
  allocate (am4_Output_edmf%vdt_edmf    (ix,jx,kx))  ; am4_Output_edmf%vdt_edmf    = 0.
  allocate (am4_Output_edmf%tdt_edmf    (ix,jx,kx))  ; am4_Output_edmf%tdt_edmf    = 0.
  allocate (am4_Output_edmf%qdt_edmf    (ix,jx,kx))  ; am4_Output_edmf%qdt_edmf    = 0.
  allocate (am4_Output_edmf%qidt_edmf   (ix,jx,kx))  ; am4_Output_edmf%qidt_edmf   = 0.
  allocate (am4_Output_edmf%qldt_edmf   (ix,jx,kx))  ; am4_Output_edmf%qldt_edmf   = 0.
  allocate (am4_Output_edmf%qadt_edmf   (ix,jx,kx))  ; am4_Output_edmf%qadt_edmf   = 0.
  allocate (am4_Output_edmf%qndt_edmf   (ix,jx,kx))  ; am4_Output_edmf%qndt_edmf   = 0.
  allocate (am4_Output_edmf%qtdt_edmf   (ix,jx,kx))  ; am4_Output_edmf%qtdt_edmf   = 0.
  allocate (am4_Output_edmf%thldt_edmf  (ix,jx,kx))  ; am4_Output_edmf%thldt_edmf  = 0.
  allocate (am4_Output_edmf%edmf_a      (ix,jx,kx))  ; am4_Output_edmf%edmf_a      = 0.
  allocate (am4_Output_edmf%edmf_w      (ix,jx,kx))  ; am4_Output_edmf%edmf_w      = 0.
  allocate (am4_Output_edmf%edmf_qt     (ix,jx,kx))  ; am4_Output_edmf%edmf_qt     = 0.
  allocate (am4_Output_edmf%edmf_thl    (ix,jx,kx))  ; am4_Output_edmf%edmf_thl    = 0.
  allocate (am4_Output_edmf%edmf_ent    (ix,jx,kx))  ; am4_Output_edmf%edmf_ent    = 0.
  allocate (am4_Output_edmf%edmf_qc     (ix,jx,kx))  ; am4_Output_edmf%edmf_qc     = 0.
  allocate (am4_Output_edmf%thl_edmf    (ix,jx,kx))  ; am4_Output_edmf%thl_edmf    = 0.
  allocate (am4_Output_edmf%qt_edmf     (ix,jx,kx))  ; am4_Output_edmf%qt_edmf     = 0.
  allocate (am4_Output_edmf%cldfra_bl   (ix,jx,kx))  ; am4_Output_edmf%cldfra_bl   = 0.
  allocate (am4_Output_edmf%qc_bl       (ix,jx,kx))  ; am4_Output_edmf%qc_bl       = 0.
  allocate (am4_Output_edmf%pbltop      (ix,jx))     ; am4_Output_edmf%pbltop      = 0.
  allocate (am4_Output_edmf%diff_t_edmf (ix,jx,kx))  ; am4_Output_edmf%diff_t_edmf = 0.
  allocate (am4_Output_edmf%diff_m_edmf (ix,jx,kx))  ; am4_Output_edmf%diff_m_edmf = 0.
  allocate (am4_Output_edmf%kpbl_edmf   (ix,jx))     ; am4_Output_edmf%kpbl_edmf   = 0
  allocate (am4_Output_edmf%el_edmf     (ix,jx,kx))  ; am4_Output_edmf%el_edmf     = 0.

  !--- diagnostic purpose
  allocate (am4_Output_edmf%t_input         (ix,jx,kx))  ; am4_Output_edmf%t_input         = 0.
  allocate (am4_Output_edmf%q_input         (ix,jx,kx))  ; am4_Output_edmf%q_input         = 0.
  allocate (am4_Output_edmf%qa_input        (ix,jx,kx))  ; am4_Output_edmf%qa_input        = 0.
  allocate (am4_Output_edmf%ql_input        (ix,jx,kx))  ; am4_Output_edmf%ql_input        = 0.
  allocate (am4_Output_edmf%qi_input        (ix,jx,kx))  ; am4_Output_edmf%qi_input        = 0.
  allocate (am4_Output_edmf%thl_input       (ix,jx,kx))  ; am4_Output_edmf%thl_input       = 0.
  allocate (am4_Output_edmf%qt_input        (ix,jx,kx))  ; am4_Output_edmf%qt_input        = 0.
  allocate (am4_Output_edmf%rh_input        (ix,jx,kx))  ; am4_Output_edmf%rh_input        = 0.
  allocate (am4_Output_edmf%th_input        (ix,jx,kx))  ; am4_Output_edmf%th_input        = 0.
  allocate (am4_Output_edmf%t_before_mix    (ix,jx,kx))  ; am4_Output_edmf%t_before_mix    = 0.
  allocate (am4_Output_edmf%q_before_mix    (ix,jx,kx))  ; am4_Output_edmf%q_before_mix    = 0.
  allocate (am4_Output_edmf%qa_before_mix   (ix,jx,kx))  ; am4_Output_edmf%qa_before_mix   = 0.
  allocate (am4_Output_edmf%ql_before_mix   (ix,jx,kx))  ; am4_Output_edmf%ql_before_mix   = 0.
  allocate (am4_Output_edmf%qi_before_mix   (ix,jx,kx))  ; am4_Output_edmf%qi_before_mix   = 0.
  allocate (am4_Output_edmf%thl_before_mix  (ix,jx,kx))  ; am4_Output_edmf%thl_before_mix  = 0.
  allocate (am4_Output_edmf%qt_before_mix   (ix,jx,kx))  ; am4_Output_edmf%qt_before_mix   = 0.
  allocate (am4_Output_edmf%rh_before_mix   (ix,jx,kx))  ; am4_Output_edmf%rh_before_mix   = 0.
  allocate (am4_Output_edmf%th_before_mix   (ix,jx,kx))  ; am4_Output_edmf%th_before_mix   = 0.
  allocate (am4_Output_edmf%t_after_mix     (ix,jx,kx))  ; am4_Output_edmf%t_after_mix     = 0.
  allocate (am4_Output_edmf%q_after_mix     (ix,jx,kx))  ; am4_Output_edmf%q_after_mix     = 0.
  allocate (am4_Output_edmf%qa_after_mix    (ix,jx,kx))  ; am4_Output_edmf%qa_after_mix    = 0.
  allocate (am4_Output_edmf%ql_after_mix    (ix,jx,kx))  ; am4_Output_edmf%ql_after_mix    = 0.
  allocate (am4_Output_edmf%qi_after_mix    (ix,jx,kx))  ; am4_Output_edmf%qi_after_mix    = 0.
  allocate (am4_Output_edmf%thl_after_mix   (ix,jx,kx))  ; am4_Output_edmf%thl_after_mix   = 0.
  allocate (am4_Output_edmf%qt_after_mix    (ix,jx,kx))  ; am4_Output_edmf%qt_after_mix    = 0.
  allocate (am4_Output_edmf%rh_after_mix    (ix,jx,kx))  ; am4_Output_edmf%rh_after_mix    = 0.
  allocate (am4_Output_edmf%th_after_mix    (ix,jx,kx))  ; am4_Output_edmf%th_after_mix    = 0.
  allocate (am4_Output_edmf%qa_before_pdf   (ix,jx,kx))  ; am4_Output_edmf%qa_before_pdf   = 0.
  allocate (am4_Output_edmf%ql_before_pdf   (ix,jx,kx))  ; am4_Output_edmf%ql_before_pdf   = 0.
  allocate (am4_Output_edmf%qi_before_pdf   (ix,jx,kx))  ; am4_Output_edmf%qi_before_pdf   = 0.
  allocate (am4_Output_edmf%Q_ql            (ix,jx,kx))  ; am4_Output_edmf%Q_ql            = 0.
  allocate (am4_Output_edmf%Q_qi            (ix,jx,kx))  ; am4_Output_edmf%Q_qi            = 0.
  allocate (am4_Output_edmf%Q_qa            (ix,jx,kx))  ; am4_Output_edmf%Q_qa            = 0.
  allocate (am4_Output_edmf%a_moist_half    (ix,jx,kx+1))  ; am4_Output_edmf%a_moist_half  = 0.
  allocate (am4_Output_edmf%mf_moist_half   (ix,jx,kx+1))  ; am4_Output_edmf%mf_moist_half = 0.
  allocate (am4_Output_edmf%qv_moist_half   (ix,jx,kx+1))  ; am4_Output_edmf%qv_moist_half = 0.
  allocate (am4_Output_edmf%a_moist_full    (ix,jx,kx))    ; am4_Output_edmf%a_moist_full  = 0.
  allocate (am4_Output_edmf%mf_moist_full   (ix,jx,kx))    ; am4_Output_edmf%mf_moist_full = 0.
  allocate (am4_Output_edmf%qv_moist_full   (ix,jx,kx))    ; am4_Output_edmf%qv_moist_full = 0.

  allocate (am4_Output_edmf%a_dry_half      (ix,jx,kx+1))  ; am4_Output_edmf%a_dry_half    = 0.
  allocate (am4_Output_edmf%a_dry_full      (ix,jx,kx))    ; am4_Output_edmf%a_dry_full    = 0.
  allocate (am4_Output_edmf%mf_dry_half     (ix,jx,kx+1))  ; am4_Output_edmf%mf_dry_half   = 0.
  allocate (am4_Output_edmf%mf_dry_full     (ix,jx,kx))    ; am4_Output_edmf%mf_dry_full   = 0.
  allocate (am4_Output_edmf%qv_dry_half     (ix,jx,kx+1))  ; am4_Output_edmf%qv_dry_half   = 0.
  allocate (am4_Output_edmf%qv_dry_full     (ix,jx,kx))    ; am4_Output_edmf%qv_dry_full   = 0.

  allocate (am4_Output_edmf%mf_all_half     (ix,jx,kx+1))  ; am4_Output_edmf%mf_all_half   = 0.
  allocate (am4_Output_edmf%mf_all_full     (ix,jx,kx))    ; am4_Output_edmf%mf_all_full   = 0.
  !allocate (am4_Output_edmf%         (ix,jx,kx))  ; am4_Output_edmf%         = 0.

!-------------------------------------------------------------------------
! set values for Input_edmf 
!-------------------------------------------------------------------------

  !--- 0D variables
  Input_edmf%delt = dt
  Input_edmf%dx   = sqrt(area(1,1))

  !--- 3-D variable
  omega_host (:,:,:) = Physics_input_block%omega
  p_host     (:,:,:) = Physics_input_block%p_full

  if (do_use_tau) then  ! use T,q,u,v,clouds at the current step
    u_host     (:,:,:) = Physics_input_block%u
    v_host     (:,:,:) = Physics_input_block%v
    T3D_host   (:,:,:) = Physics_input_block%t
    qv_host    (:,:,:) = Physics_input_block%q(:,:,:,nsphum)
    ql_host    (:,:,:) = Physics_input_block%q(:,:,:,nql)
    qi_host    (:,:,:) = Physics_input_block%q(:,:,:,nqi)
    qa_host    (:,:,:) = Physics_input_block%q(:,:,:,nqa)

  else                  ! use updated T,q,u,v,clouds
    u_host     (:,:,:) = Physics_input_block%u(:,:,:) + udt(:,:,:)*dt
    v_host     (:,:,:) = Physics_input_block%v(:,:,:) + vdt(:,:,:)*dt
    T3D_host   (:,:,:) = Physics_input_block%t(:,:,:) + tdt(:,:,:)*dt

    qv_host    (:,:,:) = Physics_input_block%q(:,:,:,nsphum) +    &
                             rdt(:,:,:,nsphum)*dt
    qa_host    (:,:,:) = Physics_input_block%q(:,:,:,nqa)    +    &
                             rdt(:,:,:,nqa)   *dt
    ql_host    (:,:,:) = Physics_input_block%q(:,:,:,nql)    +    &
                             rdt(:,:,:,nql)   *dt
    qi_host    (:,:,:) = Physics_input_block%q(:,:,:,nqi)    +    &
                             rdt(:,:,:,nqi)   *dt
  endif

  qc_host    (:,:,:) = ql_host(:,:,:) + qi_host(:,:,:) 
  exner_host(:,:,:) = (p_host(:,:,:)*p00inv)**(kappa)
  th_host   (:,:,:) = T3D_host(:,:,:) / exner_host(:,:,:)
  tv_host   (:,:,:) = T3D_host(:,:,:)*(qv_host(:,:,:)*d608+1.0)
  thv_host  (:,:,:) = tv_host(:,:,:) / exner_host(:,:,:)
  rho_host  (:,:,:) = p_host(:,:,:) / rdgas / tv_host(:,:,:)
  w_host    (:,:,:) = -1.*omega_host(:,:,:) / rho_host(:,:,:) / g 

  thl_host   (:,:,:) = th_host(:,:,:) - ( hlv*ql_host(:,:,:)+hls*qi_host(:,:,:) ) / cp_air / exner_host(:,:,:)
  qt_host    (:,:,:) = qv_host(:,:,:) + ql_host(:,:,:) + qi_host(:,:,:)

  qnc_host = 0.  ! not used, set to zero
  qni_host = 0.  ! not used, set to zero

  do k=1,kx
    dz_host(:,:,k) = Physics_input_block%z_half(:,:,k) - Physics_input_block%z_half(:,:,k+1)
  enddo

  !--- 2-D variables
  do i=1,ix
  do j=1,jx
    if (frac_land(i,j) .ge. 0.5) then
      Input_edmf%xland(i,j) = 1.   ! land
    else
      Input_edmf%xland(i,j) = 0.   ! water
    endif
  enddo
  enddo

  ! surface air virtual temperature (K)
  tv_ref (:,:) = t_ref(:,:)*(1.+d608*q_ref(:,:))
  thv_ref(:,:) = tv_ref(:,:) / (Physics_input_block%p_half(:,:,kxp)*p00inv)**(kappa)
  rho_ref(:,:) = Physics_input_block%p_half(:,:,kxp) / rdgas / tv_ref (:,:)

  ! compute surface fluxes and Obukhov length from star quantities
  u_star_star           (:,:) = u_star(:,:)
  Obukhov_length_star   (:,:) = -1.*u_star_star(:,:)**2 / vonkarm / b_star(:,:)
  w1_q1_surf_star       (:,:) = u_star(:,:) * q_star(:,:)
  w1_thv1_surf_star     (:,:) = u_star(:,:) * b_star(:,:) * thv_host(:,:,kx)/g
  w1_th1_surf_star      (:,:) = (w1_thv1_surf_star(:,:) - d608*th_host(:,:,kx)*w1_q1_surf_star(:,:)) &
                                 / (1.+d608*qv_host(:,:,kx))
  w1_t1_surf_star       (:,:) = w1_th1_surf_star * exner_host(:,:,kx)
  lhflx_star            (:,:) = w1_q1_surf_star(:,:) * rho_host(:,:,kx)
  shflx_star            (:,:) = w1_t1_surf_star(:,:) * rho_host(:,:,kx) * cp_air

  ! compute surface fluxes and Obukhov length from quantities that after surface-atmosphere coupling
  u_star_updated           (:,:) = sqrt( sqrt( u_flux(:,:)**2+v_flux(:,:)**2 ) / rho_host(:,:,kx) )
  shflx_updated            (:,:) = shflx(:,:)
  lhflx_updated            (:,:) = lhflx(:,:)
  w1_t1_surf_updated       (:,:) = shflx_updated(:,:) / rho_host(:,:,kx) / cp_air
  w1_th1_surf_updated      (:,:) = w1_t1_surf_updated(:,:) / exner_host(:,:,kx)
  w1_q1_surf_updated       (:,:) = lhflx(:,:) / rho_host(:,:,kx)
  w1_thv1_surf_updated     (:,:) = (1.+d608*qv_host(:,:,kx)) * w1_th1_surf_updated(:,:) + d608*th_host(:,:,kx)*w1_q1_surf_updated(:,:)
  Obukhov_length_updated   (:,:) = -1.*u_star_updated(:,:)**3*thv_host(:,:,kx) / vonkarm / g / w1_thv1_surf_updated(:,:)

  ! assign to Input_edmf
  if (option_surface_flux.eq."star") then
    u_star_host           (:,:) = u_star_star         (:,:)
    lhflx_host            (:,:) = lhflx_star          (:,:)
    shflx_host            (:,:) = shflx_star          (:,:)
    Obukhov_length_host   (:,:) = Obukhov_length_star (:,:)

  elseif (option_surface_flux.eq."updated") then
    u_star_host           (:,:) = u_star_updated         (:,:)
    lhflx_host            (:,:) = lhflx_updated          (:,:)
    shflx_host            (:,:) = shflx_updated          (:,:)
    Obukhov_length_host   (:,:) = Obukhov_length_updated (:,:)
  else
    call error_mesg('edmf_mynn_driver',  &
       'option_surface_flux must be star or updated', FATAL)
  endif

  ! set values in Input_edmf
  Input_edmf%ust   = u_star_host(:,:)
  Input_edmf%ch    = 0.   ! no need in JPL EDMF scheme 
  Input_edmf%rmol  = 1./Obukhov_length_host(:,:) 
  Input_edmf%ts    = t_ref   
  Input_edmf%qsfc  = 0.   ! no need in JPL EDMF scheme   
  Input_edmf%qcg   = 0.   ! no need in JPL EDMF scheme
  Input_edmf%ps    = Physics_input_block%p_half(:,:,kx+1)
  Input_edmf%hfx   = shflx_host(:,:)
  Input_edmf%qfx   = lhflx_host(:,:)
  Input_edmf%wspd  = sqrt(Physics_input_block%u(:,:,kx)**2 + Physics_input_block%v(:,:,kx)**2)
  Input_edmf%uoce  = 0.   ! set to zero for now. This is needed when specifying the surface condition for u and v wind components.
  Input_edmf%voce  = 0.   ! set to zero for now. This is needed when specifying the surface condition for u and v wind components.
  Input_edmf%vdfg  = 0.   
  Input_edmf%znt   = 0.   ! no need in JPL EDMF scheme   

  do i=1,ix    
  do j=1,jx    
  do k=1,kx 
    kk=kx-k+1   
    Input_edmf%dz    (i,kk,j) = dz_host    (i,j,k)
    Input_edmf%u     (i,kk,j) = u_host     (i,j,k)
    Input_edmf%v     (i,kk,j) = v_host     (i,j,k)
    Input_edmf%w     (i,kk,j) = w_host     (i,j,k)
    Input_edmf%th    (i,kk,j) = th_host    (i,j,k)
    Input_edmf%qv    (i,kk,j) = qv_host    (i,j,k)
    Input_edmf%p     (i,kk,j) = p_host     (i,j,k)
    Input_edmf%exner (i,kk,j) = exner_host (i,j,k)
    Input_edmf%rho   (i,kk,j) = rho_host   (i,j,k)
    Input_edmf%T3D   (i,kk,j) = T3D_host   (i,j,k)
    Input_edmf%qc    (i,kk,j) = qc_host    (i,j,k)
    Input_edmf%ql    (i,kk,j) = ql_host    (i,j,k)
    Input_edmf%qi    (i,kk,j) = qi_host    (i,j,k)
    Input_edmf%qa    (i,kk,j) = qa_host    (i,j,k)
    Input_edmf%qnc   (i,kk,j) = qnc_host   (i,j,k)
    Input_edmf%qni   (i,kk,j) = qni_host   (i,j,k)
  enddo  ! end loop of k
  enddo  ! end loop of j
  enddo  ! end loop of i

  ! make sure qc,ql,qi,qa,qnc,qni value is not smaller than tracer_min
  where (Input_edmf%qc.lt.tracer_min)
    Input_edmf%qc = 0.
  endwhere
  where (Input_edmf%ql.lt.tracer_min)
    Input_edmf%ql = 0.
  endwhere
  where (Input_edmf%qi.lt.tracer_min)
    Input_edmf%qi = 0.
  endwhere
  where (Input_edmf%qa.lt.tracer_min)
    Input_edmf%qa = 0.
  endwhere
  where (Input_edmf%qnc.lt.tracer_min)
    Input_edmf%qnc = 0.
  endwhere
  where (Input_edmf%qni.lt.tracer_min)
    Input_edmf%qni = 0.
  endwhere

  ! diagnostic purpose
  Input_edmf%u_star_star         = u_star_star (:,:)
  Input_edmf%shflx_star          = shflx_star (:,:)
  Input_edmf%lhflx_star          = lhflx_star (:,:)
  Input_edmf%w1_thv1_surf_star   = w1_thv1_surf_star (:,:)
  Input_edmf%w1_th1_surf_star    = w1_th1_surf_star (:,:)
  Input_edmf%w1_q1_surf_star     = w1_q1_surf_star (:,:)
  Input_edmf%Obukhov_length_star = Obukhov_length_star (:,:)

  Input_edmf%u_star_updated         = u_star_updated (:,:)
  Input_edmf%shflx_updated          = shflx_updated (:,:)
  Input_edmf%lhflx_updated          = lhflx_updated (:,:)
  Input_edmf%w1_thv1_surf_updated   = w1_thv1_surf_updated (:,:)
  Input_edmf%w1_th1_surf_updated    = w1_th1_surf_updated (:,:)
  Input_edmf%w1_q1_surf_updated     = w1_q1_surf_updated (:,:)
  Input_edmf%Obukhov_length_updated = Obukhov_length_updated (:,:)

!-------------------------------------------------------------------------
! set values for Output_edmf 
!-------------------------------------------------------------------------

  ! semi-prognostic variables
  do i=1,ix    
  do j=1,jx    
  do k=1,kx 
    kk=kx-k+1   
    Input_edmf%Qke        (i,k,j) = Qke       (i,j,k)   ! for write out purpose so the vertical dimension is same as AM4
    Input_edmf%el_pbl     (i,k,j) = el_pbl    (i,j,k)
    Input_edmf%cldfra_bl  (i,k,j) = cldfra_bl (i,j,k)
    Input_edmf%qc_bl      (i,k,j) = qc_bl     (i,j,k)
    Input_edmf%Sh3D       (i,k,j) = Sh3D      (i,j,k)

    Output_edmf%Qke       (i,kk,j) = Qke       (i,j,k) 
    Output_edmf%el_pbl    (i,kk,j) = el_pbl    (i,j,k)
    Output_edmf%cldfra_bl (i,kk,j) = cldfra_bl (i,j,k)
    Output_edmf%qc_bl     (i,kk,j) = qc_bl     (i,j,k)
    Output_edmf%Sh3D      (i,kk,j) = Sh3D      (i,j,k)
  enddo  ! end loop of k
  enddo  ! end loop of j
  enddo  ! end loop of i

!-------------------------------------------------------------------------
! set values for am4_Output_edmf 
!-------------------------------------------------------------------------
  am4_Output_edmf%thl_edmf (:,:,:) = thl_host (:,:,:)
  am4_Output_edmf%qt_edmf  (:,:,:) = qt_host  (:,:,:)

!-------------------------------------------------------------------------
! save the input values for diagnostic purpose
!-------------------------------------------------------------------------
  am4_Output_edmf%t_input   (:,:,:) = T3D_host   (:,:,:)
  am4_Output_edmf%q_input   (:,:,:) = qv_host    (:,:,:)
  am4_Output_edmf%qa_input  (:,:,:) = qa_host    (:,:,:)
  am4_Output_edmf%ql_input  (:,:,:) = ql_host    (:,:,:)
  am4_Output_edmf%qi_input  (:,:,:) = qi_host    (:,:,:)
  am4_Output_edmf%thl_input (:,:,:) = thl_host   (:,:,:) 
  am4_Output_edmf%qt_input  (:,:,:) = qt_host    (:,:,:)
  am4_Output_edmf%th_input  (:,:,:) = th_host    (:,:,:)

  ! amip run this will fail, comment out for a moment. 2021-05-03
  !call rh_calc (Physics_input_block%p_full(:,:,:), am4_Output_edmf%t_input(:,:,:),  &
  !              am4_Output_edmf%q_input(:,:,:), rh_host(:,:,:), do_simple )
  !am4_Output_edmf%rh_input  (:,:,:) = rh_host(:,:,:)*100.

!-------------------------------

end subroutine edmf_alloc

!###################################

subroutine edmf_dealloc (Input_edmf, Output_edmf, am4_Output_edmf)

!---------------------------------------------------------------------
! Arguments (Intent inout)
!
  type(edmf_input_type)     , intent(inout) :: Input_edmf
  type(edmf_output_type)    , intent(inout) :: Output_edmf
  type(am4_edmf_output_type), intent(inout) :: am4_Output_edmf
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    deallocate the components of the aerosol_type variable.
!---------------------------------------------------------------------

!******************
!--- Input_edmf
!******************
 
  ! diagnostic purpose
  deallocate (Input_edmf%u_star_star)          
  deallocate (Input_edmf%shflx_star)           
  deallocate (Input_edmf%lhflx_star)           
  deallocate (Input_edmf%w1_thv1_surf_star)    
  deallocate (Input_edmf%w1_th1_surf_star)     
  deallocate (Input_edmf%w1_q1_surf_star)      
  deallocate (Input_edmf%Obukhov_length_star)  

  deallocate (Input_edmf%u_star_updated)          
  deallocate (Input_edmf%shflx_updated)           
  deallocate (Input_edmf%lhflx_updated)           
  deallocate (Input_edmf%w1_thv1_surf_updated)    
  deallocate (Input_edmf%w1_th1_surf_updated)     
  deallocate (Input_edmf%w1_q1_surf_updated)      
  deallocate (Input_edmf%Obukhov_length_updated)  
  !deallocate (Input_edmf%)  
 
  ! index variable
  deallocate (Input_edmf%IMS) 
  deallocate (Input_edmf%IME) 
  deallocate (Input_edmf%JMS) 
  deallocate (Input_edmf%JME) 
  deallocate (Input_edmf%KMS) 
  deallocate (Input_edmf%KME) 

  deallocate (Input_edmf%IDS) 
  deallocate (Input_edmf%IDE) 
  deallocate (Input_edmf%JDS) 
  deallocate (Input_edmf%JDE) 
  deallocate (Input_edmf%KDS) 
  deallocate (Input_edmf%KDE) 

  deallocate (Input_edmf%ITS) 
  deallocate (Input_edmf%ITE) 
  deallocate (Input_edmf%JTS) 
  deallocate (Input_edmf%JTE) 
  deallocate (Input_edmf%KTS) 
  deallocate (Input_edmf%KTE) 

  ! 0-D variable
  deallocate (Input_edmf%delt) 
  deallocate (Input_edmf%dx)   
  !deallocate (Input_edmf%)   

  ! 2-D variable
  deallocate (Input_edmf%xland )  
  deallocate (Input_edmf%ust   )  
  deallocate (Input_edmf%ch    )  
  deallocate (Input_edmf%rmol  )  
  deallocate (Input_edmf%ts    )  
  deallocate (Input_edmf%qsfc  )  
  deallocate (Input_edmf%qcg   )  
  deallocate (Input_edmf%ps    )  
  deallocate (Input_edmf%hfx   )  
  deallocate (Input_edmf%qfx   )  
  deallocate (Input_edmf%wspd  )  
  deallocate (Input_edmf%uoce  )  
  deallocate (Input_edmf%voce  )  
  deallocate (Input_edmf%vdfg  )  
  deallocate (Input_edmf%znt   )  
  !deallocate (Input_edmf%)  

  ! 3-D variable
  deallocate (Input_edmf%dz    )  
  deallocate (Input_edmf%u     )  
  deallocate (Input_edmf%v     )  
  deallocate (Input_edmf%w     )  
  deallocate (Input_edmf%th    )  
  deallocate (Input_edmf%qv    )  
  deallocate (Input_edmf%p     )  
  deallocate (Input_edmf%exner )  
  deallocate (Input_edmf%rho   )  
  deallocate (Input_edmf%T3D   )  
  deallocate (Input_edmf%qc    )  
  deallocate (Input_edmf%ql    )  
  deallocate (Input_edmf%qi    )  
  deallocate (Input_edmf%qa    )  
  deallocate (Input_edmf%qnc   )  
  deallocate (Input_edmf%qni   )  
  !deallocate (Input_edmf%)  

  deallocate (Input_edmf%Qke         )
  deallocate (Input_edmf%Sh3D        )
  deallocate (Input_edmf%el_pbl      )
  deallocate (Input_edmf%cldfra_bl   )
  deallocate (Input_edmf%qc_bl       )


!******************
!--- Output_edmf
!******************

  ! 2-D variable
  deallocate (Output_edmf%pblh         )  
  deallocate (Output_edmf%KPBL         )  
  deallocate (Output_edmf%nupdraft     )  
  deallocate (Output_edmf%ktop_shallow )  
  deallocate (Output_edmf%maxmf        )  
  !deallocate (Output_edmf%)  

  ! 3-D variable
  deallocate (Output_edmf%Qke         )  
!print*,'alloc, Output_edmf%Qke, ix,jx,kx',size(Output_edmf%Qke,1),size(Output_edmf%Qke,2),size(Output_edmf%Qke,3)
  deallocate (Output_edmf%Tsq         )  
  deallocate (Output_edmf%Qsq         )  
  deallocate (Output_edmf%Cov         )  
  deallocate (Output_edmf%RUBLTEN     )  
  deallocate (Output_edmf%RVBLTEN     )  
  deallocate (Output_edmf%RTHBLTEN    )  
  deallocate (Output_edmf%RQVBLTEN    )  
  deallocate (Output_edmf%RQLBLTEN    )  
  deallocate (Output_edmf%RQIBLTEN    )  
  deallocate (Output_edmf%RQNIBLTEN   )  
  deallocate (Output_edmf%RTHRATEN    )  
  deallocate (Output_edmf%RCCBLTEN    )
  deallocate (Output_edmf%RTHLBLTEN   )
  deallocate (Output_edmf%RQTBLTEN    )
  deallocate (Output_edmf%edmf_a      )  
  deallocate (Output_edmf%edmf_w      )  
  deallocate (Output_edmf%edmf_qt     )  
  deallocate (Output_edmf%edmf_thl    )  
  deallocate (Output_edmf%edmf_ent    )  
  deallocate (Output_edmf%edmf_det    )  
  deallocate (Output_edmf%edmf_qc     )  
  deallocate (Output_edmf%edmf_a_dd   )  
  deallocate (Output_edmf%edmf_w_dd   )  
  deallocate (Output_edmf%edmf_qt_dd  )  
  deallocate (Output_edmf%edmf_thl_dd )  
  deallocate (Output_edmf%edmf_ent_dd )  
  deallocate (Output_edmf%edmf_qc_dd  )  
  deallocate (Output_edmf%edmf_debug1 )  
  deallocate (Output_edmf%edmf_debug2 )  
  deallocate (Output_edmf%edmf_debug3 )  
  deallocate (Output_edmf%edmf_debug4 )   
  deallocate (Output_edmf%mynn_ql     )  
  deallocate (Output_edmf%qc_bl       )  
  deallocate (Output_edmf%cldfra_bl   )  
  deallocate (Output_edmf%el_pbl      )  
  deallocate (Output_edmf%Sh3D        )  
  deallocate (Output_edmf%exch_h      )  
  deallocate (Output_edmf%exch_m      )  
  deallocate (Output_edmf%qWT         )  
  deallocate (Output_edmf%qSHEAR      )  
  deallocate (Output_edmf%qBUOY       )  
  deallocate (Output_edmf%qDISS       )  
  deallocate (Output_edmf%dqke        )  
  !deallocate (Output_edmf%)  

  deallocate (Output_edmf%t_before_mix    )
  deallocate (Output_edmf%q_before_mix    )
  deallocate (Output_edmf%qa_before_mix   )
  deallocate (Output_edmf%ql_before_mix   )
  deallocate (Output_edmf%qi_before_mix   )
  deallocate (Output_edmf%thl_before_mix  )
  deallocate (Output_edmf%qt_before_mix   )
  deallocate (Output_edmf%rh_before_mix   )
  deallocate (Output_edmf%th_before_mix   )
  deallocate (Output_edmf%t_after_mix     )
  deallocate (Output_edmf%q_after_mix     )
  deallocate (Output_edmf%qa_after_mix    )
  deallocate (Output_edmf%ql_after_mix    )
  deallocate (Output_edmf%qi_after_mix    )
  deallocate (Output_edmf%thl_after_mix   )
  deallocate (Output_edmf%qt_after_mix    )
  deallocate (Output_edmf%rh_after_mix    )
  deallocate (Output_edmf%th_after_mix    )
  deallocate (Output_edmf%qa_before_pdf   )
  deallocate (Output_edmf%ql_before_pdf   )
  deallocate (Output_edmf%qi_before_pdf   )
  deallocate (Output_edmf%Q_ql            )
  deallocate (Output_edmf%Q_qi            )
  deallocate (Output_edmf%Q_qa            )
  deallocate (Output_edmf%Q_ql_adv        )
  deallocate (Output_edmf%Q_qi_adv        )
  deallocate (Output_edmf%Q_qa_adv        )
  deallocate (Output_edmf%Q_ql_eddy       )
  deallocate (Output_edmf%Q_qi_eddy       )
  deallocate (Output_edmf%Q_qa_eddy       )
  deallocate (Output_edmf%Q_ql_ent        )
  deallocate (Output_edmf%Q_qi_ent        )
  deallocate (Output_edmf%Q_qa_ent        )
  deallocate (Output_edmf%Q_ql_det        )
  deallocate (Output_edmf%Q_qi_det        )
  deallocate (Output_edmf%Q_qa_det        )
  deallocate (Output_edmf%Q_ql_sub        )
  deallocate (Output_edmf%Q_qi_sub        )
  deallocate (Output_edmf%Q_qa_sub        )
  deallocate (Output_edmf%a_moist_half    )
  deallocate (Output_edmf%mf_moist_half   )
  deallocate (Output_edmf%qv_moist_half   )
  deallocate (Output_edmf%a_moist_full    )
  deallocate (Output_edmf%mf_moist_full   )
  deallocate (Output_edmf%qv_moist_full   )
  deallocate (Output_edmf%a_dry_half      )
  deallocate (Output_edmf%a_dry_full      )
  deallocate (Output_edmf%mf_dry_half     )
  deallocate (Output_edmf%mf_dry_full     )
  deallocate (Output_edmf%qv_dry_half     )
  deallocate (Output_edmf%qv_dry_full     )
  deallocate (Output_edmf%mf_all_half     )
  deallocate (Output_edmf%mf_all_full     )

!**********************
!--- am4 Output_edmf
!**********************

  !print*,'alloc am4_Output_edmf, ix,jx,kx',ix,jx,kx
  deallocate (am4_Output_edmf%tke         )  
  deallocate (am4_Output_edmf%Tsq         )  
  deallocate (am4_Output_edmf%Cov_thl_qt  )  
  deallocate (am4_Output_edmf%udt_edmf    )  
  deallocate (am4_Output_edmf%vdt_edmf    )  
  deallocate (am4_Output_edmf%tdt_edmf    )  
  deallocate (am4_Output_edmf%qdt_edmf    )  
  deallocate (am4_Output_edmf%qidt_edmf   )  
  deallocate (am4_Output_edmf%qldt_edmf   )  
  deallocate (am4_Output_edmf%qadt_edmf   )  
  deallocate (am4_Output_edmf%qndt_edmf   )  
  deallocate (am4_Output_edmf%qtdt_edmf   )  
  deallocate (am4_Output_edmf%thldt_edmf  )  
  deallocate (am4_Output_edmf%edmf_a      )
  deallocate (am4_Output_edmf%edmf_w      )
  deallocate (am4_Output_edmf%edmf_qt     )
  deallocate (am4_Output_edmf%edmf_thl    )
  deallocate (am4_Output_edmf%edmf_ent    )
  deallocate (am4_Output_edmf%edmf_qc     )
  deallocate (am4_Output_edmf%thl_edmf    )
  deallocate (am4_Output_edmf%qt_edmf     )
  deallocate (am4_Output_edmf%cldfra_bl   )
  deallocate (am4_Output_edmf%qc_bl       )
  deallocate (am4_Output_edmf%pbltop      )  
  deallocate (am4_Output_edmf%diff_t_edmf )  
  deallocate (am4_Output_edmf%diff_m_edmf )  
  deallocate (am4_Output_edmf%kpbl_edmf   )  
  deallocate (am4_Output_edmf%el_edmf     )  
  deallocate (am4_Output_edmf%Q_ql        )  
  deallocate (am4_Output_edmf%Q_qi        )  
  deallocate (am4_Output_edmf%Q_qa        )  
  !deallocate (am4_Output_edmf%         )  

  deallocate (am4_Output_edmf%t_input         )
  deallocate (am4_Output_edmf%q_input         )
  deallocate (am4_Output_edmf%qa_input        )
  deallocate (am4_Output_edmf%ql_input        )
  deallocate (am4_Output_edmf%qi_input        )
  deallocate (am4_Output_edmf%thl_input       )
  deallocate (am4_Output_edmf%qt_input        )
  deallocate (am4_Output_edmf%rh_input        )
  deallocate (am4_Output_edmf%th_input        )
  deallocate (am4_Output_edmf%t_before_mix    )
  deallocate (am4_Output_edmf%q_before_mix    )
  deallocate (am4_Output_edmf%qa_before_mix   )
  deallocate (am4_Output_edmf%ql_before_mix   )
  deallocate (am4_Output_edmf%qi_before_mix   )
  deallocate (am4_Output_edmf%thl_before_mix  )
  deallocate (am4_Output_edmf%qt_before_mix   )
  deallocate (am4_Output_edmf%rh_before_mix   )
  deallocate (am4_Output_edmf%th_before_mix   )
  deallocate (am4_Output_edmf%t_after_mix     )
  deallocate (am4_Output_edmf%q_after_mix     )
  deallocate (am4_Output_edmf%qa_after_mix    )
  deallocate (am4_Output_edmf%ql_after_mix    )
  deallocate (am4_Output_edmf%qi_after_mix    )
  deallocate (am4_Output_edmf%thl_after_mix   )
  deallocate (am4_Output_edmf%qt_after_mix    )
  deallocate (am4_Output_edmf%rh_after_mix    )
  deallocate (am4_Output_edmf%th_after_mix    )
  deallocate (am4_Output_edmf%qa_before_pdf   )
  deallocate (am4_Output_edmf%ql_before_pdf   )
  deallocate (am4_Output_edmf%qi_before_pdf   )
  deallocate (am4_Output_edmf%a_moist_half    )
  deallocate (am4_Output_edmf%mf_moist_half   )
  deallocate (am4_Output_edmf%qv_moist_half   )
  deallocate (am4_Output_edmf%a_moist_full    )
  deallocate (am4_Output_edmf%mf_moist_full   )
  deallocate (am4_Output_edmf%qv_moist_full   )
  deallocate (am4_Output_edmf%a_dry_half      )
  deallocate (am4_Output_edmf%a_dry_full      )
  deallocate (am4_Output_edmf%mf_dry_half     )
  deallocate (am4_Output_edmf%mf_dry_full     )
  deallocate (am4_Output_edmf%qv_dry_half     )
  deallocate (am4_Output_edmf%qv_dry_full     )
  deallocate (am4_Output_edmf%mf_all_half     )
  deallocate (am4_Output_edmf%mf_all_full     )

!--------------------
!---  vi command  ---
!--------------------
!  create deallocate	: [s/allocate/deallocate/g]
!  remove ";.*":	: [s/;.*//g]
!  remove "(IM.*ME)"	: [s/(IM.*ME)//g]
!  remove "(IM.*ME)"	: [s/(ix.*kx)//g]

end subroutine edmf_dealloc

!###################################

subroutine edmf_writeout_column ( &
              do_writeout_column, &
              is, ie, js, je, npz, Time_next, dt, lon, lat, frac_land, area, u_star,  &
              b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, rdt_mynn_ed_am4, &
              Qke, el_pbl, cldfra_bl, qc_bl, Sh3D, &
              Input_edmf, Output_edmf, am4_Output_edmf, rdiag)
!---------------------------------------------------------------------
! Arguments (Intent in)
!---------------------------------------------------------------------

  integer, intent(in)                   :: is, ie, js, je, npz
  type(time_type), intent(in)           :: Time_next
  real,    intent(in)                   :: dt
  real,    intent(in), dimension(:,:)   :: &
    frac_land, area, u_star, b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux
  type(physics_input_block_type)        :: Physics_input_block

  real,    intent(in), dimension(:,:)   :: &  ! dimension(nlon,nlat)
    lon,  &     ! longitude in radians
    lat         ! latitude  in radians

  real, intent(in), dimension(:,:,:)    :: &
    Qke, el_pbl, cldfra_bl, qc_bl, Sh3D

  type(edmf_input_type) , intent(in) :: Input_edmf
  type(edmf_output_type), intent(in) :: Output_edmf

  type(am4_edmf_output_type), intent(in) :: am4_Output_edmf
  real, intent(in), dimension(:,:,:,:) :: &   ! Mellor-Yamada, use this in offline mode
  !real, intent(in), dimension(:,:,:,ntp+1:) :: &
    rdiag

  real, intent(in), dimension(:,:,:,:) :: &
    rdt_mynn_ed_am4

  logical, intent(in) :: do_writeout_column

!---------------------------------------------------------------------
! local variables  
!---------------------------------------------------------------------
  integer :: ix, jx, kx, nt, kxp
  integer :: i, j, k, kk
  real    :: lat_lower, lat_upper, lon_lower, lon_upper, lat_temp, lon_temp

  real, dimension(size(Physics_input_block%t,3)) ::  &
    var_temp1

  real, dimension (size(Physics_input_block%t,1), &
                   size(Physics_input_block%t,2), &
                   size(Physics_input_block%t,3)+1) :: &
    diag_half

  real, dimension (size(Physics_input_block%t,1), &
                   size(Physics_input_block%t,2), &
                   size(Physics_input_block%t,3)) :: &
    rh, qsat, diag_full

  real ::  &
    tk, qtk, qtdtk, rhok, dzk, &
    tt1, tt2, tt3

  real :: tt_temp, pp_temp, esat_edmf, qq_gfdl, esat_gfdl, liquid_frac

!-------------------------------------------------------------------------
!  define input array sizes.
!-------------------------------------------------------------------------
  ix = size(Physics_input_block%t,1)
  jx = size(Physics_input_block%t,2)
  kx = size(Physics_input_block%t,3)
  nt = size(Physics_input_block%q,4)

  kxp = kx + 1

!-------------------------------------------------------------------------
!  diagnose fields 
!-------------------------------------------------------------------------

  !--- compute rh
  !call rh_calc (Physics_input_block%p_full(:,:,:), Physics_input_block%t(:,:,:),  &
  !              Physics_input_block%q(:,:,:,nsphum), rh(:,:,:), do_simple )
  !rh(:,:,:) = 100. * rh(:,:,:)

  !call compute_qs(Physics_input_block%t, Physics_input_block%p_full, qsat )
  !rh(:,:,:) = 100. * Physics_input_block%q(:,:,:,nsphum) / qsat(:,:,:)

!-------------------------------------------------------------------------
! check tracer concentration 
!   if needed, enabling codes in convert_edmf_to_am4_array to modify cloud tendencies to prevent negative tendencies
!
!   03-01-2022 Note: In the offline program with SCM_am4p0_RF01_02, the updated cloud fraction become negative (-1.e-8)
!                    Although I set am4_Output_edmf%qadt_edmf = -1*Physics_input_block%q(i,j,k,nqa)/dt, the problem is still
!                    I thought there is some numeric problems in the offline program.
!-------------------------------------------------------------------------
!  if (do_check_realizability) then
!    i=1
!    j=1
!
!    !--- qa
!    do k=1,kx
!      kk=kx-k+1
!      tt1 = Physics_input_block%q(i,j,k,nqa)
!      tt2 = Physics_input_block%q(i,j,k,nqa) + am4_Output_edmf%qadt_edmf(i,j,k) * dt
!      tt3 = Output_edmf%cldfra_bl(i,kk,j)
!      if (tt2.ne.tt3) then
!        print*,'k,qa_new,qa_bl,qa_old',k,tt2,tt3,tt1
!      endif
!    enddo
!
!    !--- ql
!    write(6,*) ''
!    do k=1,kx
!      kk=kx-k+1
!      tt1 = Physics_input_block%q(i,j,k,nql)
!      tt2 = Physics_input_block%q(i,j,k,nql) + am4_Output_edmf%qldt_edmf(i,j,k) * dt
!      tt3 = Output_edmf%qc_bl(i,kk,j)
!      if (tt2.ne.tt3) then
!        print*,'k,ql_new,qc_bl,ql_old',k,tt2,tt3,tt1
!      endif
!    enddo
!
!    write(6,*) ''
!  endif

!-------------------------------------------------------------------------
! check water and energy conservation 
!-------------------------------------------------------------------------

  if (do_check_consrv) then
    tt1 = 0.
    tt2 = 0.
    tt3 = 0.

    !i=ii_write
    !j=jj_write
    i=1
    j=1

    do k=1,kx
      kk=kx-k+1
      rhok    = Input_edmf%rho(i,kk,j)
      tk      = Physics_input_block%t(i,j,k)
      qtk     =   Physics_input_block%q(i,j,k,nsphum)  &
                + Physics_input_block%q(i,j,k,nql)     &
                + Physics_input_block%q(i,j,k,nqi)
      qtdtk   =   am4_Output_edmf%qdt_edmf(i,j,k)     &
                + am4_Output_edmf%qldt_edmf(i,j,k)    &
                + am4_Output_edmf%qidt_edmf(i,j,k)

      !qtk     =   Physics_input_block%q(i,j,k,nsphum)
      dzk  = Physics_input_block%z_half(i,j,k) - Physics_input_block%z_half(i,j,k+1)

      tt1 = tt1 + qtk  * rhok * dzk
      tt2 = tt2 + qtdtk* rhok * dzk
      tt3 = tt3 + tk   * rhok * dzk

      !print*,'k,rhok,qtk,dzk',k,rhok,qtk,dzk
      !print*,'k,qdt,qldt,qidt',k,am4_Output_edmf%qdt_edmf(i,j,k),am4_Output_edmf%qldt_edmf(i,j,k),am4_Output_edmf%qidt_edmf(i,j,k)
      !print*,'k,rhok,tk,',k,rhok,tk
    enddo

    print*,'column-integrated total water (kg/m2) = ',tt1
    print*,'column-integrated total water tendency (kg/m2/s) = ',tt2
    print*,'moisture flux, qfx (kg/m2/s) = ',Input_edmf%qfx(i,j)
    print*,'diff column-integrated total water tendency - qfx (kg/m2/s)= ',tt2 - Input_edmf%qfx(i,j)
    print*,'column-integrated rho*T*dz = ',tt3

  endif  ! end if of do_check_consrv

!-------------------------------------------------------------------------
! writing out the selected column
!-------------------------------------------------------------------------

! debug purpose 
!           do k=1,kx
!             kk=kx-k+1
!             var_temp1(k) = Input_edmf%th(1,kk,1)
!           enddo
!        write(6,3001) '  th = (/'    ,var_temp1(:)

  if (do_writeout_column) then
        write(6,*)    ';@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        write(6,*)    '; i,j,',ii_write,jj_write
        write(6,*)    '; lat',lat (ii_write,jj_write)
        write(6,*)    '; lon',lon (ii_write,jj_write)
        write(6,*)    ';@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        write(6,*)    ''
        write(6,*)    ';=========================='
        write(6,*)    ';=========================='
        write(6,*)    ''
        write(6,*)    ';   Physics_input_block '
        write(6,*)    ''
        write(6,*)    ';=========================='
        write(6,*)    ';=========================='
        write(6,*)    ''
        write(6,*)    '; pressure at half level (Pa)'
        write(6,3001) '  p_half  = (/',Physics_input_block%p_half(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; pressure at full level (Pa)'
        write(6,3001) '  p_full  = (/',Physics_input_block%p_full(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; actual height at half level (m)'
        write(6,3001) '  z_half  = (/',Physics_input_block%z_half(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; actual height at full level (m)'
        write(6,3001) '  z_full  = (/',Physics_input_block%z_full(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; height at half level above the surface (m)'
        write(6,3001) '  z_half_surf0  = (/',Physics_input_block%z_half(ii_write,jj_write,:) - Physics_input_block%z_half(ii_write,jj_write,kxp)
        write(6,*)    ''
        write(6,*)    '; height at full level above the surface (m)'
        write(6,3001) '  z_full_surf0  = (/',Physics_input_block%z_full(ii_write,jj_write,:) - Physics_input_block%z_half(ii_write,jj_write,kxp)
        write(6,*)    ''
        write(6,*)    '; zonal wind velocity at full levels (m/s)'
        !write(6,3001) '  uu  = (/'    ,Physics_input_block%u(ii_write,jj_write,:)
           do k=1,kx
             kk=kx-k+1
             var_temp1(k) = Input_edmf%u(ii_write,kk,jj_write)
           enddo
        write(6,3001) '  uu  = (/'    ,var_temp1(:)
        write(6,*)    ''
        write(6,*)    '; meridional wind velocity at full levels (m/s)'
        !write(6,3001) '  vv  = (/'    ,Physics_input_block%v(ii_write,jj_write,:)
           do k=1,kx
             kk=kx-k+1
             var_temp1(k) = Input_edmf%v(ii_write,kk,jj_write)
           enddo
        write(6,3001) '  vv  = (/'    ,var_temp1(:)
        write(6,*)    ''
        write(6,*)    '; vertical velocity (Pa/s)'
        write(6,3002) '  omega  = (/'    ,Physics_input_block%omega(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; temperatur at full levels (K)'
        !write(6,3001) '  tt  = (/'    ,Physics_input_block%t(ii_write,jj_write,:)
           do k=1,kx
             kk=kx-k+1
             var_temp1(k) = Input_edmf%T3D(ii_write,kk,jj_write)
           enddo
        write(6,3001) '  tt  = (/'    ,var_temp1(:)
        write(6,*)    ''
        write(6,*)    '; potential temperatur at full levels (K)'
           do k=1,kx
             kk=kx-k+1
             var_temp1(k) = Input_edmf%th(ii_write,kk,jj_write)
           enddo
        write(6,3001) '  th  = (/'    ,var_temp1(:)
        write(6,*)    ''
        write(6,*)    '; ice-liquid water potential temperatur at full levels (K)'
        write(6,3001) '  thl  = (/'    , am4_Output_edmf%thl_edmf(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; input relative humidity (%)'
        write(6,3001) '  rh  = (/'    ,am4_Output_edmf%rh_input(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; specific humidity at full levels (kg/kg)'
        !write(6,3002) '  qq  = (/'    ,Physics_input_block%q(ii_write,jj_write,:,nsphum)
           do k=1,kx
             kk=kx-k+1
             var_temp1(k) = Input_edmf%qv(ii_write,kk,jj_write)
           enddo
        write(6,3002) '  qq  = (/'    ,var_temp1(:)
        write(6,*)    ''
        write(6,*)    '; cloud fraction (none)'
        !write(6,3002) '  qa  = (/'    ,Physics_input_block%q(ii_write,jj_write,:,nqa)
           do k=1,kx
             kk=kx-k+1
             var_temp1(k) = Input_edmf%qa(ii_write,kk,jj_write)
           enddo
        write(6,3002) '  qa  = (/'    ,var_temp1(:)
        write(6,*)    ''
        write(6,*)    '; cloud liquid water mixing ratio at full levels (kg/kg)'
        !write(6,3002) '  ql  = (/'    ,Physics_input_block%q(ii_write,jj_write,:,nql)
           do k=1,kx
             kk=kx-k+1
             var_temp1(k) = Input_edmf%ql(ii_write,kk,jj_write)
           enddo
        write(6,3002) '  ql  = (/'    ,var_temp1(:)
        write(6,*)    ''
        write(6,*)    '; cloud ice water mixing ratio at full levels (kg/kg)'
        !write(6,3002) '  qi  = (/'    ,Physics_input_block%q(ii_write,jj_write,:,nqi)
           do k=1,kx
             kk=kx-k+1
             var_temp1(k) = Input_edmf%qi(ii_write,kk,jj_write)
           enddo
        write(6,3002) '  qi  = (/'    ,var_temp1(:)
        write(6,*)    ''
        write(6,*)    '; total water mixing ratio (qv+ql+qi) at full levels (kg/kg)'
        write(6,3002) '  qt  = (/'    ,am4_Output_edmf%qt_edmf(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; surface air temperature (K)'
        write(6,3000) '  t_ref= (/'    , t_ref (ii_write,jj_write)
        write(6,*)    ''
        write(6,*)    '; surface air specific humidity (kg/kg)'
        write(6,3003) '  q_ref = (/'    , q_ref (ii_write,jj_write)
        write(6,*)    ''
        write(6,*)    '; zonal wind stress'
        write(6,3002) '  u_flux = (/'    ,u_flux (ii_write,jj_write)
        write(6,*)    ''
        write(6,*)    '; meridional wind stress'
        write(6,3002) '  v_flux = (/'    ,v_flux (ii_write,jj_write)
        write(6,*)    ''
        write(6,*)    '; '
        write(6,3003) '  u_star = (/',u_star(ii_write,jj_write)
        write(6,3003) '  b_star = (/',b_star(ii_write,jj_write)
        write(6,3003) '  q_star = (/',q_star(ii_write,jj_write)
        write(6,3003) '  shflx  = (/',shflx(ii_write,jj_write)
        write(6,3003) '  lhflx  = (/',lhflx(ii_write,jj_write)
        write(6,*)    ' '
        write(6,*)    '; rdiag(1,1,:,nQke), input'
        write(6,3002) '  rdiag(1,1,:,nQke) = (/'    ,Input_edmf%Qke(ii_write,:,jj_write)
        write(6,*)    ' '
        write(6,*)    '; rdiag(1,1,:,nSh3D), input'
        write(6,3002) '  rdiag(1,1,:,nSh3D) = (/'    ,Input_edmf%Sh3D(ii_write,:,jj_write)
        write(6,*)    ' '
        write(6,*)    '; rdiag(1,1,:,nel_pbl), input'
        write(6,3002) '  rdiag(1,1,:,nel_pbl) = (/'    ,Input_edmf%el_pbl(ii_write,:,jj_write)
        write(6,*)    ' '
        write(6,*)    '; rdiag(1,1,:,ncldfra_bl), input'
        write(6,3002) '  rdiag(1,1,:,ncldfra_bl) = (/'    ,Input_edmf%cldfra_bl(ii_write,:,jj_write)
        write(6,*)    ' '
        write(6,*)    '; rdiag(1,1,:,nqc_bl), input'
        write(6,3002) '  rdiag(1,1,:,nqc_bl) = (/'    ,Input_edmf%qc_bl(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; rdt_mynn_ed_am4(1,1,:,nql)'
        write(6,3002) '  rdt_mynn_ed_am4(1,1,:,nql) = (/'    ,rdt_mynn_ed_am4(ii_write,jj_write,:,nql)
        write(6,*)    ''
        write(6,*)    '; rdt_mynn_ed_am4(1,1,:,nqa)'
        write(6,3002) '  rdt_mynn_ed_am4(1,1,:,nqa) = (/'    ,rdt_mynn_ed_am4(ii_write,jj_write,:,nqa)
        write(6,*)    ''
        write(6,*)    '; rdt_mynn_ed_am4(1,1,:,nqi)'
        write(6,3002) '  rdt_mynn_ed_am4(1,1,:,nqi) = (/'    ,rdt_mynn_ed_am4(ii_write,jj_write,:,nqi)
        write(6,*)    ''
        write(6,*)    ';----- end of fieles needed by the offline program ---'
        write(6,*)    ''
        write(6,*)    '; friction velocity (m/s)'
        write(6,3003) '  u_star_star    = (/',Input_edmf%u_star_star(ii_write,jj_write)
        write(6,3003) '  u_star_updated = (/',Input_edmf%u_star_updated(ii_write,jj_write)
        write(6,*)    ''
        write(6,*)    '; sensible heat flux (W/m2)'
        write(6,3003) '  shflx_star    = (/',Input_edmf%shflx_star(ii_write,jj_write)
        write(6,3003) '  shflx_updated = (/',Input_edmf%shflx_updated(ii_write,jj_write)
        write(6,*)    ''
        write(6,*)    '; evaporation flux (kg water/m2/s)'
        write(6,3003) '  lhflx_star    = (/',Input_edmf%lhflx_star(ii_write,jj_write)
        write(6,3003) '  lhflx_updated = (/',Input_edmf%lhflx_updated(ii_write,jj_write)
        write(6,*)    ''
        write(6,*)    '; surface heat flux (K m/s)'
        write(6,3003) '  w1_th1_surf_star    = (/',Input_edmf%w1_th1_surf_star(ii_write,jj_write)
        write(6,3003) '  w1_th1_surf_updated = (/',Input_edmf%w1_th1_surf_updated(ii_write,jj_write)
        write(6,*)    ''
        write(6,*)    '; surface moisture flux (kg water/kg air * m/s)'
        write(6,3003) '  w1_q1_surf_star    = (/',Input_edmf%w1_q1_surf_star(ii_write,jj_write)
        write(6,3003) '  w1_q1_surf_updated = (/',Input_edmf%w1_q1_surf_updated(ii_write,jj_write)
        write(6,*)    ''
        write(6,*)    '; kinematic virtual temperature flux (K m/s)'
        write(6,3003) '  w1_thv1_surf_star    = (/',Input_edmf%w1_thv1_surf_star(ii_write,jj_write)
        write(6,3003) '  w1_thv1_surf_updated = (/',Input_edmf%w1_thv1_surf_updated(ii_write,jj_write)
        write(6,*)    ''
        write(6,*)    '; Obukhov length (m)'
        write(6,3000) '  Obukhov_length_star    = (/',Input_edmf%Obukhov_length_star(ii_write,jj_write)
        write(6,3000) '  Obukhov_length_updated = (/',Input_edmf%Obukhov_length_updated(ii_write,jj_write)
        !write(6,*)    ' '
        !write(6,*)    '; '
        !write(6,3002) '  = (/'    ,
        write(6,*)    ''
        write(6,*)    ';==========================='
        write(6,*)    ';==========================='
        write(6,*)    ''
        write(6,*)    ';   am4_Output_edmf'
        write(6,*)    ''
        write(6,*)    ';==========================='
        write(6,*)    ';==========================='
        write(6,*)    ''
        write(6,*)    '; turbulent kinetic energy (m2/s2)'
        write(6,3002) ' tke = (/'    ,am4_Output_edmf%tke(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; variance of theta_l (K^2)'
        write(6,3002) ' Tsq = (/'    ,am4_Output_edmf%Tsq (ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; covariance of theta_l and q_t (none)'
        write(6,3002) ' Cov_thl_qt = (/'    ,am4_Output_edmf%Cov_thl_qt (ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; EDMF diffusion coefficients for heat (K m/s)'
        write(6,3002) ' diff_t_edmf = (/'    ,am4_Output_edmf%diff_t_edmf (ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; EDMF diffusion coefficients for momentum (m2/s)'
        write(6,3002) ' diff_m_edmf = (/'    ,am4_Output_edmf%diff_m_edmf (ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; vertical velocity of updrafts [m s^-1]'
                      call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_w, diag_half)
        write(6,3002) ' edmf_w = (/'    ,diag_half(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; updraft area [-]'
                      call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_a, diag_half)
        write(6,3002) ' edmf_a = (/'    ,diag_half(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; ensemble updraft mass flux [kg m^-2 s^-1]'
                      call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%mf_all_half, diag_half)
        write(6,3002) ' mf_all_half = (/'    ,diag_half(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; moist updraft mass flux [kg m^-2 s^-1]'
                      call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%mf_moist_half, diag_half)
        write(6,3002) ' mf_moist_half = (/'    ,diag_half(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; dry updraft mass flux [kg m^-2 s^-1]'
                      call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%mf_dry_half, diag_half)
        write(6,3002) ' mf_dry_half = (/'    ,diag_half(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; dry updrafts area at half levels []'
                      call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%a_dry_half, diag_half)
        write(6,3002) ' a_dry_half = (/'    ,diag_half(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; moist updrafts area at half levels []'
                      call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%a_moist_half, diag_half)
        write(6,3002) ' a_moist_half = (/'    ,diag_half(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; ensemble-mean thl in updrafts [kg/kg]'
                      call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_thl, diag_half)
        write(6,3002) ' edmf_thl = (/'    ,diag_half(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; ensemble-mean qt in updrafts [kg/kg]'
                      call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_qt, diag_half)
        write(6,3002) ' edmf_qt = (/'    ,diag_half(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; ensemble-mean qc in updrafts [kg/kg]'
                      call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_qc, diag_half)
        write(6,3002) ' edmf_qc = (/'    ,diag_half(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; dry updrafts area at full levels []'
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%a_dry_full, diag_full)
        write(6,3002) ' a_dry_full = (/'    ,diag_full(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; moist updrafts area at full levels []'
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%a_moist_full, diag_full)
        write(6,3002) ' a_moist_full = (/'    ,diag_full(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; qa before mix []'
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%qa_before_mix, diag_full)
        write(6,3002) ' qa_before_mix = (/'    ,diag_full(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; entrainment in updrafts [m^-1]'
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%edmf_ent, diag_full)
        write(6,3002) ' edmf_ent = (/'  ,diag_full(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; detrainment in updrafts [m^-1]'
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%edmf_det, diag_full)
        write(6,3002) ' edmf_det = (/'  ,diag_full(ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; u tendency from edmf_mynn (m/s2)'
        write(6,3002) ' udt_edmf = (/'    ,am4_Output_edmf%udt_edmf (ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; v tendency from edmf_mynn (m/s2)'
        write(6,3002) ' vdt_edmf = (/'    ,am4_Output_edmf%vdt_edmf (ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; t tendency from edmf_mynn (K/s)'
        write(6,3002) ' tdt_edmf = (/'    ,am4_Output_edmf%tdt_edmf (ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; q tendency from edmf_mynn (kg/kg/s)'
        write(6,3002) ' qdt_edmf = (/'    ,am4_Output_edmf%qdt_edmf (ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; qi tendency from edmf_mynn (kg/kg/s)'
        write(6,3002) ' qidt_edmf = (/'    ,am4_Output_edmf%qidt_edmf (ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; qc tendency from edmf_mynn (kg/kg/s)'
        write(6,3002) ' qldt_edmf = (/'    ,am4_Output_edmf%qldt_edmf (ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; qt tendency from edmf_mynn (kg/kg/s)'
        write(6,3002) ' qtdt_edmf = (/'    ,am4_Output_edmf%qtdt_edmf (ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; qa tendency from edmf_mynn (1/s)'
        write(6,3002) ' qadt_edmf = (/'    ,am4_Output_edmf%qadt_edmf (ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; thl tendency from edmf_mynn (K/s)'
        write(6,3002) ' thldt_edmf = (/'    ,am4_Output_edmf%thldt_edmf (ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; cloud fraction in edmf_mynn'
        write(6,3002) ' cldfra_bl = (/'    ,am4_Output_edmf%cldfra_bl (ii_write,jj_write,:)
        write(6,*)    ' '
        write(6,*)    '; liquid water mixing ratio in edmf_mynn'
        write(6,3002) ' qc_bl = (/'    ,am4_Output_edmf%qc_bl (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa, diag_full)
        write(6,*)    '; Q_qa'
        write(6,3002) '  Q_qa = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_adv, diag_full)
        write(6,*)    '; Q_qa_adv'
        write(6,3002) '  Q_qa_adv = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_eddy, diag_full)
        write(6,*)    '; Q_qa_eddy'
        write(6,3002) '  Q_qa_eddy = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_ent, diag_full)
        write(6,*)    '; Q_qa_ent'
        write(6,3002) '  Q_qa_ent = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_det, diag_full)
        write(6,*)    '; Q_qa_det'
        write(6,3002) '  Q_qa_det = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_sub, diag_full)
        write(6,*)    '; Q_qa_sub'
        write(6,3002) '  Q_qa_sub = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ''
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql, diag_full)
        write(6,*)    '; Q_ql'
        write(6,3002) '  Q_ql = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_adv, diag_full)
        write(6,*)    '; Q_ql_adv'
        write(6,3002) '  Q_ql_adv = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_eddy, diag_full)
        write(6,*)    '; Q_ql_eddy'
        write(6,3002) '  Q_ql_eddy = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_ent, diag_full)
        write(6,*)    '; Q_ql_ent'
        write(6,3002) '  Q_ql_ent = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_det, diag_full)
        write(6,*)    '; Q_ql_det'
        write(6,3002) '  Q_ql_det = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_sub, diag_full)
        write(6,*)    '; Q_ql_sub'
        write(6,3002) '  Q_ql_sub = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ''
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi, diag_full)
        write(6,*)    '; Q_qi'
        write(6,3002) '  Q_qi = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_adv, diag_full)
        write(6,*)    '; Q_qi_adv'
        write(6,3002) '  Q_qi_adv = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_eddy, diag_full)
        write(6,*)    '; Q_qi_eddy'
        write(6,3002) '  Q_qi_eddy = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_ent, diag_full)
        write(6,*)    '; Q_qi_ent'
        write(6,3002) '  Q_qi_ent = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_det, diag_full)
        write(6,*)    '; Q_qi_det'
        write(6,3002) '  Q_qi_det = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ' '
                      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_sub, diag_full)
        write(6,*)    '; Q_qi_sub'
        write(6,3002) '  Q_qi_sub = (/'    ,diag_full (ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; k_pbl, ',am4_Output_edmf%kpbl_edmf(ii_write,jj_write)
        write(6,*)    ''
        write(6,*)    '; pbl depth,',am4_Output_edmf%pbltop(ii_write,jj_write)
        write(6,*)    ''
        write(6,*)    ';=============='
        write(6,*)    ';=============='
        write(6,*)    ''
        write(6,*)    ';   rdiag after mynn'
        write(6,*)    ''
        write(6,*)    ';=============='
        write(6,*)    ';=============='
        write(6,*)    ''
        write(6,*)    'rdiag_Qke',Qke(ii_write,jj_write,:)
        write(6,*)    'rdiag_Sh3D',Sh3D(ii_write,jj_write,:)
        write(6,*)    'rdiag_el_pbl',el_pbl(ii_write,jj_write,:)
        write(6,*)    'rdiag_cldfra_bl',cldfra_bl(ii_write,jj_write,:)
        write(6,*)    'rdiag_qc_bl',qc_bl(ii_write,jj_write,:)
        !write(6,*)    'rdiag_',rdiag(ii_write,jj_write,:,n)
        write(6,*)    ''
        write(6,*)    ';======================'
        write(6,*)    ';======================'
        write(6,*)    ''
        write(6,*)    ';   Input_edmf'
        write(6,*)    ''
        write(6,*)    ';======================'
        write(6,*)    ';======================'
        write(6,*)    ''
        write(6,*)    ';--- 0-D fields'
        write(6,*)    'Input_edmf%delt', Input_edmf%delt
        write(6,*)    'Input_edmf%dx', Input_edmf%dx
        write(6,*)    ''
        write(6,*)    ';--- 2-D fields'
        write(6,*)    'Input_edmf%ust', Input_edmf%ust(ii_write,jj_write)
        write(6,*)    'Input_edmf%rmol,', Input_edmf%rmol(ii_write,jj_write)
        write(6,*)    'Input_edmf%ps',Input_edmf%ps(ii_write,jj_write)
        write(6,*)    'Input_edmf%hfx',Input_edmf%hfx(ii_write,jj_write)
        write(6,*)    'Input_edmf%qfx',Input_edmf%qfx(ii_write,jj_write)
        write(6,*)    'Input_edmf%wspd',Input_edmf%wspd(ii_write,jj_write)
        write(6,*)    ''
        write(6,*)    ';--- 3-D fields'
        write(6,*)    'Input_edmf%dz',Input_edmf%dz(ii_write,:,jj_write)
        write(6,*)    'Input_edmf%u ',Input_edmf%u (ii_write,:,jj_write)
        write(6,*)    'Input_edmf%v',Input_edmf%v(ii_write,:,jj_write)
        write(6,*)    'Input_edmf%w',Input_edmf%w(ii_write,:,jj_write)
        write(6,*)    'Input_edmf%th',Input_edmf%th(ii_write,:,jj_write)
        write(6,*)    'Input_edmf%qv',Input_edmf%qv(ii_write,:,jj_write)
        write(6,*)    'Input_edmf%p',Input_edmf%p(ii_write,:,jj_write)
        write(6,*)    'Input_edmf%exner',Input_edmf%exner(ii_write,:,jj_write)
        write(6,*)    'Input_edmf%rho',Input_edmf%rho(ii_write,:,jj_write)
        write(6,*)    'Input_edmf%T3D',Input_edmf%T3D (ii_write,:,jj_write)
        write(6,*)    'Input_edmf%qc',Input_edmf%qc(ii_write,:,jj_write)
        write(6,*)    'Input_edmf%qi',Input_edmf%qi(ii_write,:,jj_write)
        write(6,*)    'Input_edmf%qa',Input_edmf%qa(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    ';======================'
        write(6,*)    ';======================'
        write(6,*)    ''
        write(6,*)    ';   Output_edmf'
        write(6,*)    ''
        write(6,*)    ';======================'
        write(6,*)    ';======================'
        write(6,*)    ''
        write(6,*)    ''
        write(6,*)    '; twice the turbulent kinetic energy [m^2 s^-2]'
        !write(*,"(A35,2X,34(E12.4,2X,','))"),' Qke = (/',Output_edmf%Qke(ii_write,:,jj_write)   ! don't know why, format 3002 doesn't work
        write(6,3002) ' Qke = (/'    ,Output_edmf%Qke(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; variance of theta_l (K^2)'
        write(6,3002) ' Tsq = (/'    ,Output_edmf%Tsq(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; variance of q_t (-)'
        write(6,3002) ' Qsq = (/'    ,Output_edmf%Qsq(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; covariance of theta_l and q_t  (K)'
        write(6,3002) ' Cov = (/'    ,Output_edmf%Cov(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; U tendency due to EDMF parameterization (m/s^2)' 
        write(6,3002) ' RUBLTEN = (/'    ,Output_edmf%RUBLTEN(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; V tendency due to EDMF parameterization (m/s^2)'
        write(6,3002) ' RVBLTEN = (/'    ,Output_edmf%RVBLTEN(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; Theta tendency due to EDMF parameterization (K/s)'
        write(6,3002) ' RTHBLTEN = (/'    ,Output_edmf%RTHBLTEN(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; Qv tendency due to EDMF parameterization (kg/kg/s)'
        write(6,3002) ' RQVBLTEN = (/'    ,Output_edmf%RQVBLTEN(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; Ql tendency due to EDMF parameterization (kg/kg/s)'
        write(6,3002) ' RQLBLTEN = (/'    ,Output_edmf%RQLBLTEN(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; Qi tendency due to  EDMF parameterization (kg/kg/s)'
        write(6,3002) ' RQIBLTEN = (/'    ,Output_edmf%RQIBLTEN(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; Qni tendency due to EDMF parameterization (#/kg/s)'
        write(6,3002) ' RQNIBLTEN = (/'    ,Output_edmf%RQNIBLTEN(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; tendency from radiation [K s^-1]'
        write(6,3002) ' RTHRATEN = (/'    ,Output_edmf%RTHRATEN(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; updraft area [-]'
        write(6,3002) ' edmf_a = (/'    ,Output_edmf%edmf_a(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; vertical velocity of updrafts [m s^-1]'
        write(6,3002) ' edmf_w = (/'    ,Output_edmf%edmf_w(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; qt in updrafts  [kg kg^-1]'
        write(6,3002) ' edmf_qt = (/'    ,Output_edmf%edmf_qt(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; thl in updrafts  [K]'
        write(6,3002) ' edmf_thl = (/'    ,Output_edmf%edmf_thl(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; entrainment in updrafts [m^-1]'
        write(6,3002) ' edmf_ent = (/'    ,Output_edmf%edmf_ent(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; qc in updrafts [kg kg^-1]'
        write(6,3002) ' edmf_qc = (/'    ,Output_edmf%edmf_qc(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; qt_before_mix [kg kg^-1]'
        write(6,3002) '  qt_before_mix = (/'    ,Output_edmf%qt_before_mix(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; qt_after_mix [kg kg^-1]'
        write(6,3002) '  qt_after_mix = (/'    ,Output_edmf%qt_after_mix(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; qa_before_pdf'
        write(6,3002) '  qa_before_pdf = (/'    ,Output_edmf%qa_before_pdf(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; qa_before_mix'
        write(6,3002) '  qa_before_mix = (/'    ,Output_edmf%qa_before_mix(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; qa_after_mix'
        write(6,3002) '  qa_after_mix = (/'    ,Output_edmf%qa_after_mix(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; ql_before_pdf [kg kg^-1]'
        write(6,3002) '  ql_before_pdf = (/'    ,Output_edmf%ql_before_pdf(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; ql_before_mix [kg kg^-1]'
        write(6,3002) '  ql_before_mix = (/'    ,Output_edmf%ql_before_mix(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; ql_after_mix [kg kg^-1]'
        write(6,3002) '  ql_after_mix = (/'    ,Output_edmf%ql_after_mix(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; Output_edmf%exch_h'
        write(6,3002) ' exch_h = (/', Output_edmf%exch_h(ii_write,:,jj_write)
        write(6,*)    ''
        write(6,*)    '; Output_edmf%exch_m'
        write(6,3002) ' exch_m = (/', Output_edmf%exch_m(ii_write,:,jj_write)
        write(6,*)    ''
        !--------------------
        !write(6,*)    'Output_edmf%a_dry_half,',Output_edmf%a_dry_half(ii_write,:,jj_write)
        !write(6,*)    'Output_edmf%a_moist_half,',Output_edmf%a_moist_half(ii_write,:,jj_write)
        !write(6,*)    'Output_edmf%mf_dry_half,',Output_edmf%mf_dry_half(ii_write,:,jj_write)
        !write(6,*)    'Output_edmf%mf_moist_half,',Output_edmf%mf_moist_half(ii_write,:,jj_write)
        !--------------------
        write(6,*)    ';-----------------------------'
        write(6,*)    ';  Some vi commands'
        write(6,*)    ';-----------------------------'
        write(6,*)    ''
        write(6,*)    ';*** NCL '
        write(6,*)    '; remove the space at the beginning	: [1,$s/^ */  /g]'
        write(6,*)    '; replace ",$" to "/)" at the end	: [1,$s/,$/\/)/g]'
        write(6,*)    ''
        write(6,*)    ';*** Fortran'
        write(6,*)    '; remove the space at the beginning	: [1,$s/^ */  /g]'
        write(6,*)    '; replace ";" to "!"			: [1,$s/;/!/g]'
        write(6,*)    '; replace ",$" to "/)" at the end	: [1,$s/,$/\/)/g]'
        write(6,*)    ''
        write(6,*)    '; for 3-D variable'
        write(6,*)    '; replace "=" to "(1,1,:)"		: [1,Ns/ = /(1,1,:) = /g]'
        write(6,*)    ''
        write(6,*)    '; for 2-D variable'
        write(6,*)    '; replace "(/" to ""                	: [N,$s/(\///g]' 
        write(6,*)    '; replace "/)" to "" 			: [N,$s/\/)//g]'
        write(6,*)    ''
  endif  ! end if of do_writeout_column

!******************************************* 
!  --- check saturation vapor pressure from EDMF and GFDL
!******************************************* 
!
!  write(6,*) 'T(K), esat_edmf (Pa), esat_gfdl (Pa)'
!do i=195,350
!  tt_temp = float(i)  ! temperature in K
!
!  !--- EDMF formula
!  IF (tt_temp .GE. 273.16) THEN
!    liquid_frac=1.
!  ELSE IF (tt_temp .LE. 253.) THEN
!    liquid_frac=0.
!  ELSE
!    liquid_frac  = 1.-(273.16-tt_temp)/20.16
!  END IF
!  
!  esat_edmf = esatLF_blend(tt_temp, liquid_frac)
!
!  !--- GFDL formula
!  pp_temp = 1013. * 100. ! pressure in Pa
!  call compute_qs(tt_temp, pp_temp, qq_gfdl, esat = esat_gfdl )
!
!  write(6,*) tt_temp, esat_edmf, esat_gfdl
!enddo
!**************************** end check saturation vapor pressure

3000 format (A35,2X,F10.3,',')
3001 format (A35,2X,34(F10.3,2X,','))
3002 format (A35,2X,34(E12.4,2X,','))
3003 format (A35,2X,E12.4,',')
3004 format (A35,2X,33(F10.3,2X,','),A5)

end subroutine edmf_writeout_column

!########################

subroutine convert_edmf_to_am4_array (Physics_input_block, ix, jx, kx, &
                                      Input_edmf, Output_edmf, am4_Output_edmf, rdiag, &
                                      Qke, el_pbl, cldfra_bl, qc_bl, Sh3D )

!--- input arguments
  type(physics_input_block_type), intent(in)  :: Physics_input_block 
  type(edmf_input_type)     , intent(in)  :: Input_edmf
  type(edmf_output_type)    , intent(in)  :: Output_edmf
  integer                   , intent(in)  :: ix, jx, kx

!--- output arguments
  type(am4_edmf_output_type), intent(inout) :: am4_Output_edmf
  real, dimension (:,:,:,:) , intent(inout) :: rdiag
  real, intent(inout), dimension(:,:,:)    :: &
    Qke, el_pbl, cldfra_bl, qc_bl, Sh3D

!--- local variable
  integer i,j,k,kk
  real :: &
    dum, dum1, qa1, qc1, qi1, qt1
!------------------------------------------

!---------------
! 2D variables
!---------------
  do i=1,ix
  do j=1,jx
    kk=kx-Output_edmf%kpbl (i,j)+1 

    am4_Output_edmf%pbltop   (i,j) = Physics_input_block%z_full(i,j,kk) - Physics_input_block%z_half(i,j,kx+1)
    am4_Output_edmf%kpbl_edmf(i,j) = kk
  enddo  ! end loop of j
  enddo  ! end loop of 1

!  !--- fixed pbl top for testing
!  if (do_pblh_constant) then
!    do i=1,ix
!    do j=1,jx
!      dum1 = 1.e+10
!      kk=1
!      do k=1,kx    
!        dum  = abs(Physics_input_block%z_full(i,j,k) - fixed_pblh)
!        if (dum.lt.dum1) then
!          dum1 = dum
!          kk = k
!        endif
!      enddo
!
!      am4_Output_edmf%pbltop   (i,j) = Physics_input_block%z_full(i,j,kk) - Physics_input_block%z_half(i,j,kx+1)
!      am4_Output_edmf%kpbl_edmf(i,j) = kk      
!    enddo
!    enddo
!  endif  ! end if of do_pblh_constant

!---------------
! 3D variables
!---------------

  do i=1,ix
  do j=1,jx

    !======================
    !  diff_t_edmf and diff_m_edmf are on half levels but the dimension is nlay, instead of nlay+1
    do k=2,kx      ! k index for half levels. Skip k=1 because k=1 in MYNN is right at the surface
      kk=kx-k+2

      !--- edmf_* variables are written out to history files, no need to convert to am4_Output_edmf variables
      !    also note that although edmf_* dimension is nlay, but they are on half levels (nlay+1)
      !am4_Output_edmf%edmf_a      (i,j,kk) = Output_edmf%edmf_a    (i,k,j)
      !am4_Output_edmf%edmf_w      (i,j,kk) = Output_edmf%edmf_w    (i,k,j)
      !am4_Output_edmf%edmf_qt     (i,j,kk) = Output_edmf%edmf_qt   (i,k,j)
      !am4_Output_edmf%edmf_thl    (i,j,kk) = Output_edmf%edmf_thl  (i,k,j)
      !am4_Output_edmf%edmf_ent    (i,j,kk) = Output_edmf%edmf_ent  (i,k,j)
      !am4_Output_edmf%edmf_qt     (i,j,kk) = Output_edmf%edmf_qt   (i,k,j)
      !am4_Output_edmf%edmf_qc     (i,j,kk) = Output_edmf%edmf_qc   (i,k,j)
  
      am4_Output_edmf%diff_t_edmf (i,j,kk) = Output_edmf%exch_h    (i,k,j)
      am4_Output_edmf%diff_m_edmf (i,j,kk) = Output_edmf%exch_m    (i,k,j)
    enddo          ! end loop of k, half levels

    !======================
    do k=1,kx      ! k index for full levels
      kk=kx-k+1
  
      !--- set am4_Output_edmf, on GFDL full levels
      am4_Output_edmf%tke         (i,j,kk) = 0.5 * Output_edmf%Qke (i,k,j)
      am4_Output_edmf%Tsq         (i,j,kk) = Output_edmf%Tsq       (i,k,j)
      am4_Output_edmf%Cov_thl_qt  (i,j,kk) = Output_edmf%Cov       (i,k,j)
      am4_Output_edmf%udt_edmf    (i,j,kk) = Output_edmf%RUBLTEN   (i,k,j) 
      am4_Output_edmf%vdt_edmf    (i,j,kk) = Output_edmf%RVBLTEN   (i,k,j)
      am4_Output_edmf%tdt_edmf    (i,j,kk) = Output_edmf%RTHBLTEN  (i,k,j) * Input_edmf%exner (i,k,j)
      am4_Output_edmf%qdt_edmf    (i,j,kk) = Output_edmf%RQVBLTEN  (i,k,j)
      am4_Output_edmf%qidt_edmf   (i,j,kk) = Output_edmf%RQIBLTEN  (i,k,j)
      am4_Output_edmf%qldt_edmf   (i,j,kk) = Output_edmf%RQLBLTEN  (i,k,j)
      am4_Output_edmf%qadt_edmf   (i,j,kk) = Output_edmf%RCCBLTEN  (i,k,j)
      am4_Output_edmf%qtdt_edmf   (i,j,kk) = Output_edmf%RQTBLTEN  (i,k,j)
      am4_Output_edmf%thldt_edmf  (i,j,kk) = Output_edmf%RTHLBLTEN (i,k,j)
      am4_Output_edmf%cldfra_bl   (i,j,kk) = Output_edmf%cldfra_bl (i,k,j)
      am4_Output_edmf%qc_bl       (i,j,kk) = Output_edmf%qc_bl     (i,k,j)
      am4_Output_edmf%el_edmf     (i,j,kk) = Output_edmf%el_pbl    (i,k,j)
      am4_Output_edmf%Q_ql        (i,j,kk) = Output_edmf%Q_ql      (i,k,j)
      am4_Output_edmf%Q_qi        (i,j,kk) = Output_edmf%Q_qi      (i,k,j)
      am4_Output_edmf%Q_qa        (i,j,kk) = Output_edmf%Q_qa      (i,k,j)

      am4_Output_edmf%a_moist_full   (i,j,kk) = Output_edmf%a_moist_full   (i,k,j)
      am4_Output_edmf%mf_moist_full  (i,j,kk) = Output_edmf%mf_moist_full  (i,k,j)
      am4_Output_edmf%qv_moist_full  (i,j,kk) = Output_edmf%qv_moist_full  (i,k,j)

      am4_Output_edmf%a_dry_full     (i,j,kk) = Output_edmf%a_dry_full   (i,k,j)
      am4_Output_edmf%mf_dry_full    (i,j,kk) = Output_edmf%mf_dry_full  (i,k,j)
      am4_Output_edmf%qv_dry_full    (i,j,kk) = Output_edmf%qv_dry_full  (i,k,j)
 
      am4_Output_edmf%mf_all_full    (i,j,kk) = Output_edmf%mf_all_full  (i,k,j)
      !!! am4_Output_edmf% (i,j,kk) = Output_edmf% (i,k,j)
  
      !--- if needed, modify am4_Output_edmf tendencies to make sure the updated qa, ql, qc, qi, qnd qt are larger than zero
      !qa1 = Physics_input_block%q(i,j,kk,nqa) + am4_Output_edmf%qadt_edmf(i,j,kk) * dt
      !if (qa1 .lt. 0.) then
      !  am4_Output_edmf%qadt_edmf(i,j,kk) = -1.*Physics_input_block%q(i,j,kk,nqa) / dt
      !  print*,'no,k,am4,mynn',kk,am4_Output_edmf%qadt_edmf(i,j,kk),Output_edmf%RCCBLTEN  (i,k,j), &
      !                           am4_Output_edmf%qadt_edmf(i,j,kk)-Output_edmf%RCCBLTEN  (i,k,j)
      !endif
  
  
      !--- change rdiag
      Qke       (i,j,kk) = Output_edmf%Qke       (i,k,j)
      el_pbl    (i,j,kk) = Output_edmf%el_pbl    (i,k,j)
      cldfra_bl (i,j,kk) = Output_edmf%cldfra_bl (i,k,j)
      qc_bl     (i,j,kk) = Output_edmf%qc_bl     (i,k,j)
      Sh3D      (i,j,kk) = Output_edmf%Sh3D      (i,k,j)
      !!! rdiag(i,j,kk,n)      = Output_edmf%      (i,k,j)
 
      !--- To avoid MYNN producing weird tendencies when turbulent mixing is small, 
      !    Kay Suselj suggested to set tendencies to zeros when TKE is small (e.g. <0.02 m2/s2)
      !    However, When MF is included, it is possible that updrafts can exist where TKE is very small, 
      !    which means MYNN-EDMF tendencies are not zeros. In such cases, the limiter would reset these tendencies to zeros, 
      !    removing some water/energy away and causing the conservation problem. 
      !    I think I added this limiter because if using edmf_type=1 (tendencies are computed based on the difference 
      !    between the variables from the PDF scheme and the input variables from the GFDL model). Even though there was no mixing, 
      !    it could have spurious tendencies due to the incompatibility between Tiedtke and the PDF cloud scheme. 
      !    I guess that is why I added the limiter to remove spurious tendencies when TKE is small. 
      !    Diabling this limiter to avoid a=confusion.
      !if (Qke(i,j,kk) .lt. qke_min) then
      !  am4_Output_edmf%udt_edmf    (i,j,kk) = 0.
      !  am4_Output_edmf%vdt_edmf    (i,j,kk) = 0.
      !  am4_Output_edmf%tdt_edmf    (i,j,kk) = 0.
      !  am4_Output_edmf%qdt_edmf    (i,j,kk) = 0.
      !  am4_Output_edmf%qidt_edmf   (i,j,kk) = 0.
      !  am4_Output_edmf%qldt_edmf   (i,j,kk) = 0.
      !  am4_Output_edmf%qadt_edmf   (i,j,kk) = 0.
      !  am4_Output_edmf%qtdt_edmf   (i,j,kk) = 0.
      !  am4_Output_edmf%thldt_edmf  (i,j,kk) = 0.
      !endif
  
    enddo  ! end loop of k, full levels

    !======================
    do k=1,kx+1      ! k index for half levels
      kk=kx-k+2
      am4_Output_edmf%a_moist_half   (i,j,kk) = Output_edmf%a_moist_half   (i,k,j)
      am4_Output_edmf%mf_moist_half  (i,j,kk) = Output_edmf%mf_moist_half  (i,k,j)
      am4_Output_edmf%qv_moist_half  (i,j,kk) = Output_edmf%qv_moist_half  (i,k,j)

      am4_Output_edmf%a_dry_half     (i,j,kk) = Output_edmf%a_dry_half   (i,k,j)
      am4_Output_edmf%mf_dry_half    (i,j,kk) = Output_edmf%mf_dry_half  (i,k,j)
      am4_Output_edmf%qv_dry_half    (i,j,kk) = Output_edmf%qv_dry_half  (i,k,j)
 
      am4_Output_edmf%mf_all_half    (i,j,kk) = Output_edmf%mf_all_half  (i,k,j)
    enddo  ! end loop of k, half levels

  enddo  ! end loop of j
  enddo  ! end loop of 1

!----------------------------
! variables for diagnostics
!----------------------------
  do i=1,ix
  do j=1,jx
    do k=1,kx      ! k index for full levels
      kk=kx-k+1

      !--- before mixing
      am4_Output_edmf%qa_before_mix   (i,j,kk) = Output_edmf%qa_before_mix   (i,k,j)
      am4_Output_edmf%ql_before_mix   (i,j,kk) = Output_edmf%ql_before_mix   (i,k,j)
      am4_Output_edmf%qi_before_mix   (i,j,kk) = Output_edmf%qi_before_mix   (i,k,j)
      am4_Output_edmf%thl_before_mix  (i,j,kk) = Output_edmf%thl_before_mix  (i,k,j)
      am4_Output_edmf%qt_before_mix   (i,j,kk) = Output_edmf%qt_before_mix   (i,k,j)
      am4_Output_edmf%th_before_mix   (i,j,kk) = Output_edmf%th_before_mix   (i,k,j)
      am4_Output_edmf%t_before_mix    (i,j,kk) = Output_edmf%th_before_mix   (i,k,j) * Input_edmf%exner (i,k,j)

      !--- after mixing
      am4_Output_edmf%qa_after_mix   (i,j,kk) = Output_edmf%qa_after_mix   (i,k,j)
      am4_Output_edmf%ql_after_mix   (i,j,kk) = Output_edmf%ql_after_mix   (i,k,j)
      am4_Output_edmf%qi_after_mix   (i,j,kk) = Output_edmf%qi_after_mix   (i,k,j)
      am4_Output_edmf%thl_after_mix  (i,j,kk) = Output_edmf%thl_after_mix  (i,k,j)
      am4_Output_edmf%qt_after_mix   (i,j,kk) = Output_edmf%qt_after_mix   (i,k,j)
      am4_Output_edmf%th_after_mix   (i,j,kk) = Output_edmf%th_after_mix   (i,k,j)
      am4_Output_edmf%t_after_mix    (i,j,kk) = Output_edmf%th_after_mix   (i,k,j) * Input_edmf%exner (i,k,j)

      !--- before mixing, from the PDF cloud scheme
      am4_Output_edmf%qa_before_pdf   (i,j,kk) = Output_edmf%qa_before_pdf   (i,k,j)
      am4_Output_edmf%ql_before_pdf   (i,j,kk) = Output_edmf%ql_before_pdf   (i,k,j)
      am4_Output_edmf%qi_before_pdf   (i,j,kk) = Output_edmf%qi_before_pdf   (i,k,j)

    enddo  ! end loop of k, full levels
  enddo  ! end loop of j
  enddo  ! end loop of 1

  am4_Output_edmf%q_before_mix (:,:,:) =    am4_Output_edmf%qt_before_mix(:,:,:)   &
                                         -  am4_Output_edmf%ql_before_mix(:,:,:)   &
                                         -  am4_Output_edmf%qi_before_mix(:,:,:)   
  am4_Output_edmf%q_after_mix (:,:,:)  =    am4_Output_edmf%qt_after_mix(:,:,:)    &
                                         -  am4_Output_edmf%ql_after_mix(:,:,:)    &
                                         -  am4_Output_edmf%qi_after_mix(:,:,:)    


  ! amip run this will fail, comment out for a moment. 2021-04-20
  !call rh_calc (Physics_input_block%p_full(:,:,:), am4_Output_edmf%t_before_mix(:,:,:),  &
  !              am4_Output_edmf%q_before_mix(:,:,:), am4_Output_edmf%rh_before_mix(:,:,:), do_simple )
  !am4_Output_edmf%rh_before_mix  (:,:,:) = am4_Output_edmf%rh_before_mix(:,:,:)*100.

  !call rh_calc (Physics_input_block%p_full(:,:,:), am4_Output_edmf%t_after_mix(:,:,:),  &
  !              am4_Output_edmf%q_after_mix(:,:,:), am4_Output_edmf%rh_after_mix(:,:,:), do_simple )
  !am4_Output_edmf%rh_after_mix  (:,:,:) = am4_Output_edmf%rh_after_mix(:,:,:)*100.

end subroutine convert_edmf_to_am4_array

!###################################

subroutine modify_mynn_edmf_tendencies (is, ie, js, je, Time_next, dt,     &
                                        do_writeout_column, &
                                        Physics_input_block, Input_edmf, rdt_mynn_ed_am4, &
                                        ix, jx, kx,  &
                                        Output_edmf, rlz_ratio, rlz_tracer)

!--- input arguments
  integer, intent(in)                   :: is, ie, js, je
  real,    intent(in)                   :: dt
  type(time_type), intent(in)           :: Time_next
  type(physics_input_block_type), intent(in)  :: Physics_input_block
  type(edmf_input_type)     , intent(in)  :: Input_edmf
  integer                   , intent(in)  :: ix, jx, kx
  real, dimension(:,:,:,:)  , intent(in)  :: &
    rdt_mynn_ed_am4
  logical, intent(in) :: do_writeout_column

!--- input/output arguments
  type(edmf_output_type)    , intent(inout)  :: Output_edmf
  real, dimension(ix,jx)    , intent(out)    :: rlz_ratio  ! ratio for scaling the tendencies 
  real, dimension(ix,jx)    , intent(out)    :: rlz_tracer    ! tracer name for the rlz_ratio

!--- local variable
  integer, parameter :: nx = 5     ! number of tracers for realizability check
  real, dimension(ix,jx,kx,nx) ::  &  ! tracers for realizability check
    tracers, tracers_tend
  real, dimension(ix,jx,nx)   ::  &  ! realizability ratio for each tracer 
    tends_ratio, &   ! realizability ratio for each tracer 
    tends_ratio_k    ! level of the realizability ratio
  real, dimension(ix,jx,kx) :: tmp_3d
  logical used
  integer i,j,k,kk
!------------------------------------------

!--- initialze return variables
  rlz_ratio  = 1.
  rlz_tracer = 0.

!----------------------------
! save the original mynn_edmf tendencies
!----------------------------

!send_data
   !------- t tendency from edmf_mynn original (units: K/s) at full level -------
   if ( id_tdt_edmf_orig > 0) then
     call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RTHBLTEN(:,:,:)*Input_edmf%exner(:,:,:), tmp_3d)
     used = send_data (id_tdt_edmf_orig, tmp_3d, Time_next, is, js, 1 )
   endif

   !------- q tendency from edmf_mynn original (units: kg/kg/s) at full level -------
   if ( id_qdt_edmf_orig > 0) then
     call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RQVBLTEN(:,:,:), tmp_3d)
     used = send_data (id_qdt_edmf_orig, tmp_3d, Time_next, is, js, 1 )
   endif

   !------- cldfra tendency from edmf_mynn original (units: 1/s) at full level -------
   if ( id_qadt_edmf_orig > 0) then
     call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RCCBLTEN(:,:,:), tmp_3d)
     used = send_data (id_qadt_edmf_orig, tmp_3d, Time_next, is, js, 1 )
   endif

   !------- qi tendency from edmf_mynn original (units: kg/kg/s) at full level -------
   if ( id_qidt_edmf_orig > 0) then
     call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RQIBLTEN(:,:,:), tmp_3d)
     used = send_data (id_qidt_edmf_orig, tmp_3d, Time_next, is, js, 1 )
   endif

   !------- ql tendency from edmf_mynn original (units: kg/kg/s) at full level -------
   if ( id_qldt_edmf_orig > 0) then
     call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RQLBLTEN(:,:,:), tmp_3d)
     used = send_data (id_qldt_edmf_orig, tmp_3d, Time_next, is, js, 1 )
   endif

!------- ql tendency from edmf_mynn, ED (units: kg/kg/s) at full level -------
    if ( id_qldt_edmf_ED > 0) then
      used = send_data (id_qldt_edmf_ED, rdt_mynn_ed_am4(:,:,:,nql), Time_next, is, js, 1 )
    endif

!------- ql tendency from edmf_mynn, MF (units: kg/kg/s) at full level -------
    if ( id_qldt_edmf_MF > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql(:,:,:), tmp_3d)
      used = send_data (id_qldt_edmf_MF, tmp_3d, Time_next, is, js, 1 )
    endif

!------- qi tendency from edmf_mynn, ED (units: kg/kg/s) at full level -------
    if ( id_qidt_edmf_ED > 0) then
      used = send_data (id_qidt_edmf_ED, rdt_mynn_ed_am4(:,:,:,nqi), Time_next, is, js, 1 )
    endif

!------- qi tendency from edmf_mynn, MF (units: kg/kg/s) at full level -------
    if ( id_qidt_edmf_MF > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi(:,:,:), tmp_3d)
      used = send_data (id_qidt_edmf_MF, tmp_3d, Time_next, is, js, 1 )
    endif

!------- qa tendency from edmf_mynn, ED (units: 1/s) at full level -------
    if ( id_qadt_edmf_ED > 0) then
      used = send_data (id_qadt_edmf_ED, rdt_mynn_ed_am4(:,:,:,nqa), Time_next, is, js, 1 )
    endif

!------- qa tendency from edmf_mynn, MF (units: 1/s) at full level -------
    if ( id_qadt_edmf_MF > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa(:,:,:), tmp_3d)
      used = send_data (id_qadt_edmf_MF, tmp_3d, Time_next, is, js, 1 )
    endif

!------- qa tendency from edmf_mynn, MF_adv (units: 1/s) at full level -------
    if ( id_qadt_edmf_MF_adv > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_adv(:,:,:), tmp_3d)
      used = send_data (id_qadt_edmf_MF_adv, tmp_3d, Time_next, is, js, 1 )
    endif

!------- ql tendency from edmf_mynn, MF_adv (units: kg/kg/s) at full level -------
    if ( id_qldt_edmf_MF_adv > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_adv(:,:,:), tmp_3d)
      used = send_data (id_qldt_edmf_MF_adv, tmp_3d, Time_next, is, js, 1 )
    endif

!------- qi tendency from edmf_mynn, MF_adv (units: kg/kg/s) at full level -------
    if ( id_qidt_edmf_MF_adv > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_adv(:,:,:), tmp_3d)
      used = send_data (id_qidt_edmf_MF_adv, tmp_3d, Time_next, is, js, 1 )
    endif

!------- qa tendency from edmf_mynn, MF_eddy (units: 1/s) at full level -------
    if ( id_qadt_edmf_MF_eddy > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_eddy(:,:,:), tmp_3d)
      used = send_data (id_qadt_edmf_MF_eddy, tmp_3d, Time_next, is, js, 1 )
    endif

!------- ql tendency from edmf_mynn, MF_eddy (units: kg/kg/s) at full level -------
    if ( id_qldt_edmf_MF_eddy > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_eddy(:,:,:), tmp_3d)
      used = send_data (id_qldt_edmf_MF_eddy, tmp_3d, Time_next, is, js, 1 )
    endif

!------- qi tendency from edmf_mynn, MF_eddy (units: kg/kg/s) at full level -------
    if ( id_qidt_edmf_MF_eddy > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_eddy(:,:,:), tmp_3d)
      used = send_data (id_qidt_edmf_MF_eddy, tmp_3d, Time_next, is, js, 1 )
    endif

!------- qa tendency from edmf_mynn, MF_ent (units: 1/s) at full level -------
    if ( id_qadt_edmf_MF_ent > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_ent(:,:,:), tmp_3d)
      used = send_data (id_qadt_edmf_MF_ent, tmp_3d, Time_next, is, js, 1 )
    endif

!------- ql tendency from edmf_mynn, MF_ent (units: kg/kg/s) at full level -------
    if ( id_qldt_edmf_MF_ent > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_ent(:,:,:), tmp_3d)
      used = send_data (id_qldt_edmf_MF_ent, tmp_3d, Time_next, is, js, 1 )
    endif

!------- qi tendency from edmf_mynn, MF_ent (units: kg/kg/s) at full level -------
    if ( id_qidt_edmf_MF_ent > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_ent(:,:,:), tmp_3d)
      used = send_data (id_qidt_edmf_MF_ent, tmp_3d, Time_next, is, js, 1 )
    endif

!------- qa tendency from edmf_mynn, MF_det (units: 1/s) at full level -------
    if ( id_qadt_edmf_MF_det > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_det(:,:,:), tmp_3d)
      used = send_data (id_qadt_edmf_MF_det, tmp_3d, Time_next, is, js, 1 )
    endif

!------- ql tendency from edmf_mynn, MF_det (units: kg/kg/s) at full level -------
    if ( id_qldt_edmf_MF_det > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_det(:,:,:), tmp_3d)
      used = send_data (id_qldt_edmf_MF_det, tmp_3d, Time_next, is, js, 1 )
    endif

!------- qi tendency from edmf_mynn, MF_det (units: kg/kg/s) at full level -------
    if ( id_qidt_edmf_MF_det > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_det(:,:,:), tmp_3d)
      used = send_data (id_qidt_edmf_MF_det, tmp_3d, Time_next, is, js, 1 )
    endif

!------- qa tendency from edmf_mynn, MF_sub (units: 1/s) at full level -------
    if ( id_qadt_edmf_MF_sub > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_sub(:,:,:), tmp_3d)
      used = send_data (id_qadt_edmf_MF_sub, tmp_3d, Time_next, is, js, 1 )
    endif

!------- ql tendency from edmf_mynn, MF_sub (units: kg/kg/s) at full level -------
    if ( id_qldt_edmf_MF_sub > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_sub(:,:,:), tmp_3d)
      used = send_data (id_qldt_edmf_MF_sub, tmp_3d, Time_next, is, js, 1 )
    endif

!------- qi tendency from edmf_mynn, MF_sub (units: kg/kg/s) at full level -------
    if ( id_qidt_edmf_MF_sub > 0) then
      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_sub(:,:,:), tmp_3d)
      used = send_data (id_qidt_edmf_MF_sub, tmp_3d, Time_next, is, js, 1 )
    endif
!send_data

!----------------------------
! modify mynn_edmf tendencies 
!----------------------------

   !******************************
   !---  do_option_edmf2ls_mp=1 or 2, 
   !      “evaporate/condensate” the liquid and ice water that is produced during mixing  
   !       above PBL to prevent EDMF produce weird cloud tendencies (e.g. EDMF sometime produces ~0.5 cloud fraction at ~200 hPa)
   !******************************
   if ( do_option_edmf2ls_mp.eq.1 .or. do_option_edmf2ls_mp.eq.2 ) then
     do i=1,ix
     do j=1,jx
       k = Output_edmf%kpbl(i,j) 
       Output_edmf%RQLBLTEN (:,k+1:kx,:) =  0.
       Output_edmf%RQIBLTEN (:,k+1:kx,:) =  0.
       Output_edmf%RQVBLTEN (:,k+1:kx,:) =  Output_edmf%RQTBLTEN(:,k+1:kx,:)
       Output_edmf%RTHBLTEN (:,k+1:kx,:) =  Output_edmf%RTHLBLTEN (:,k+1:kx,:)  &
                                         + (hlv*Output_edmf%RQLBLTEN (:,k+1:kx,:)+hls*Output_edmf%RQIBLTEN(:,k+1:kx,:)) / cp_air / Input_edmf%exner(:,k+1:kx,:)
      enddo  ! end loop of i
      enddo  ! end loop of j
   endif

   !******************************
   !---  do_option_edmf2ls_mp=3 
   !      “evaporate/condensate” the liquid and ice water that is produced during mixing  
   !******************************
   if ( do_option_edmf2ls_mp.eq.3 ) then
     Output_edmf%RQLBLTEN (:,:,:) =  0.
     Output_edmf%RQIBLTEN (:,:,:) =  0.
     Output_edmf%RQVBLTEN (:,:,:) =  Output_edmf%RQTBLTEN(:,:,:)
     Output_edmf%RTHBLTEN (:,:,:) =  Output_edmf%RTHLBLTEN (:,:,:)  &
                                   + (hlv*Output_edmf%RQLBLTEN (:,:,:)+hls*Output_edmf%RQIBLTEN(:,:,:)) / cp_air / Input_edmf%exner(:,:,:)
   endif

   !******************************
   !--- edmf_type=2, 
   !      recover dry variable tendencies by approximating cloud liquid/ice tendencies
   !******************************
   if (edmf_type .eq. 2) then
     do i=1,ix
     do j=1,jx
     do k=1,kx      ! k index for full levels
       kk=kx-k+1
 
       Output_edmf%RQIBLTEN  (i,k,j) = rdt_mynn_ed_am4(i,j,kk,nqi) + Output_edmf%Q_qi(i,k,j)  ! modify qi tendency
       Output_edmf%RQLBLTEN  (i,k,j) = rdt_mynn_ed_am4(i,j,kk,nql) + Output_edmf%Q_ql(i,k,j)  ! modify ql tendency
       Output_edmf%RCCBLTEN  (i,k,j) = rdt_mynn_ed_am4(i,j,kk,nqa) + Output_edmf%Q_qa(i,k,j)  ! modify qa tendency
 
     enddo  ! end loop of i
     enddo  ! end loop of j
     enddo  ! end loop of k
 
     ! modify qv tendecy, qvdt = qtdt - modified qldt & qidt
     Output_edmf%RQVBLTEN  (:,:,:) =   Output_edmf%RQTBLTEN  (:,:,:)  &
                                     - Output_edmf%RQLBLTEN  (:,:,:) - Output_edmf%RQIBLTEN  (:,:,:)
 
     ! modify theta tendency accordingly, keep theta_li tendency unchanged
     Output_edmf%RTHBLTEN  (:,:,:) =   Output_edmf%RTHLBLTEN (:,:,:)  &
                                     + (hlv*Output_edmf%RQLBLTEN (:,:,:)+hls*Output_edmf%RQIBLTEN(:,:,:)) / cp_air / Input_edmf%exner(:,:,:)
 
   end if  ! end if of edmf_type=2

!====================================
!
! Check for tracer realizability. If MF tendencies would
!  produce negative tracer mixing ratios, scale down tracer tendency
!  terms uniformly for this tracer throughout convective column.
!
!====================================

if (do_check_realizability) then

   !--- initialzie tracers varibles for realizability check
   tracers      = 0.
   tracers_tend = 0.  

   !--- assign tracers values, qv, ql, qi, qa, and qt(=qv+ql+qi)
   tracers(:,:,:,1) = Physics_input_block%q(:,:,:,nsphum) 
   tracers(:,:,:,2) = Physics_input_block%q(:,:,:,nql) 
   tracers(:,:,:,3) = Physics_input_block%q(:,:,:,nqi) 
   tracers(:,:,:,4) = Physics_input_block%q(:,:,:,nqa)
   tmp_3d (:,:,:)   = Physics_input_block%q(:,:,:,nsphum)+Physics_input_block%q(:,:,:,nql)+Physics_input_block%q(:,:,:,nqi)
   tracers(:,:,:,5) = tmp_3d(:,:,:)

   !--- assign tracers tendencies values
   call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RQVBLTEN(:,:,:), tracers_tend(:,:,:,1))
   call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RQLBLTEN(:,:,:), tracers_tend(:,:,:,2))
   call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RQIBLTEN(:,:,:), tracers_tend(:,:,:,3))
   call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RCCBLTEN(:,:,:), tracers_tend(:,:,:,4))
   call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RQTBLTEN(:,:,:), tracers_tend(:,:,:,5))

   call check_trc_rlzbility (dt, tracers, tracers_tend, &
                             tends_ratio, tends_ratio_k, rlz_ratio, rlz_tracer)

   !--- ratio must be between 0 and 1
   if (any(tends_ratio.lt.0.) .or. any(tends_ratio.gt.1.)) then
     call error_mesg( ' edmf_mynn',     &
                      ' the ratio from the realizability check must be between 0 to 1',&
                      FATAL )
   endif

   !--- printout for debugging purpose
   if (do_debug_option.eq."check_rlz" .or. do_debug_option.eq."all") then
     print*,'**********************'
     print*,''
     print*,'do_check_realizability = ',do_check_realizability
     print*,'index: qv, ql, qi, qa, qt'
     print*,'tends_ratio',tends_ratio
     print*,'tends_ratio_k',tends_ratio_k
     print*,'rlz_ratio',rlz_ratio
     print*,'rlz_tracer',rlz_tracer
     print*,'------------------------------'
     print*,'original tendencies'
     print*,''
     print*,'RUBLTEN',Output_edmf%RUBLTEN
     print*,'RVBLTEN',Output_edmf%RVBLTEN
     print*,'RQTBLTEN',Output_edmf%RQTBLTEN
     print*,'RTHLBLTEN',Output_edmf%RTHLBLTEN
     print*,'RQVBLTEN',Output_edmf%RQVBLTEN
     print*,'RQLBLTEN',Output_edmf%RQLBLTEN
     print*,'RQIBLTEN',Output_edmf%RQIBLTEN
     print*,'RCCBLTEN',Output_edmf%RCCBLTEN
     print*,''
     print*,'**********************'
   end if

   !--- rescale the tendencies
   do i=1,ix
   do j=1,jx
     ! u and v tendencies
     Output_edmf%RUBLTEN (i,:,j) = Output_edmf%RUBLTEN (i,:,j) * rlz_ratio(i,j)
     Output_edmf%RVBLTEN (i,:,j) = Output_edmf%RVBLTEN (i,:,j) * rlz_ratio(i,j)

     ! theta_l and q_t tendencies
     Output_edmf%RQTBLTEN (i,:,j) = Output_edmf%RQTBLTEN (i,:,j) * rlz_ratio(i,j)
     Output_edmf%RTHLBLTEN(i,:,j) = Output_edmf%RTHLBLTEN(i,:,j) * rlz_ratio(i,j)

     ! qv, ql, qi, qa tendencies
     Output_edmf%RQVBLTEN(i,:,j) = Output_edmf%RQVBLTEN(i,:,j) * rlz_ratio(i,j)
     Output_edmf%RQLBLTEN(i,:,j) = Output_edmf%RQLBLTEN(i,:,j) * rlz_ratio(i,j)
     Output_edmf%RQIBLTEN(i,:,j) = Output_edmf%RQIBLTEN(i,:,j) * rlz_ratio(i,j)
     Output_edmf%RCCBLTEN(i,:,j) = Output_edmf%RCCBLTEN(i,:,j) * rlz_ratio(i,j)
   enddo  ! end loop of j
   enddo  ! end loop of i

   !--- recompute the potential temperature tendencies using the rescaled theta_li, q_l, and q_i tendencies
   Output_edmf%RTHBLTEN  (:,:,:) =   Output_edmf%RTHLBLTEN (:,:,:)  &
                                   + (hlv*Output_edmf%RQLBLTEN (:,:,:)+hls*Output_edmf%RQIBLTEN(:,:,:)) / cp_air / Input_edmf%exner(:,:,:)

endif  ! end if of do_check_realizability


!----------------------------
! printout statements 
!----------------------------
!   if (do_writeout_column) then
!     write(6,*) '; i,j,',ii_write,jj_write
!     !--- qa
!     write(6,*) 'rdt_mynn_ed_am4, qa',rdt_mynn_ed_am4(ii_write,jj_write,:,nqa)
!     write(6,*) 'Output_edmf%Q_qa',Output_edmf%Q_qa(ii_write,:,jj_write)
!     write(6,*) 'Output_edmf%RCCBLTEN',Output_edmf%RCCBLTEN(ii_write,:,jj_write)
!
!     !--- ql
!     write(6,*) 'rdt_mynn_ed_am4, ql',rdt_mynn_ed_am4(ii_write,jj_write,:,nql)
!     write(6,*) 'Output_edmf%Q_ql',Output_edmf%Q_ql(ii_write,:,jj_write)
!     write(6,*) 'Output_edmf%RQLBLTEN',Output_edmf%RQLBLTEN(ii_write,:,jj_write)
!
!     !--- qi
!     write(6,*) 'rdt_mynn_ed_am4, qi',rdt_mynn_ed_am4(ii_write,jj_write,:,nqi)
!     write(6,*) 'Output_edmf%Q_qi',Output_edmf%Q_qi(ii_write,:,jj_write)
!     write(6,*) 'Output_edmf%RQIBLTEN',Output_edmf%RQIBLTEN(ii_write,:,jj_write)
!   end if

!
!send_data
!------- tracer name for realizability limiter (units: none) at one level -------
! ! yhc 2021-12-15, fail to compile, "error #6284: There is no matching specific function for this generic function reference.   [SEND_DATA]"
! !                 I have no idea why...
!      if ( id_rlz_tracer > 0) then
!        used = send_data (id_rlz_tracer, rlz_tracer, Time_next, is, js )
!      endif
!send_data

end subroutine modify_mynn_edmf_tendencies

!###################################

subroutine reshape_mynn_array_to_am4 (ix, jx, kx, mynn_array_3d, am4_array_3d)

!--- input arguments
  integer                   , intent(in)  :: ix, jx, kx
  real, dimension(ix, kx, jx), intent(in) :: &
    mynn_array_3d

!--- output arguments
  real, dimension(ix, jx, kx), intent(out) :: &
    am4_array_3d

!--- local variable
  integer i,j,k,kk
!----------------------------------

  am4_array_3d = 0.

  do i=1,ix
  do j=1,jx
    do k=1,kx      ! k index for full levels
      kk=kx-k+1
      am4_array_3d   (i,j,kk) = mynn_array_3d (i,k,j)
    enddo  ! end loop of k, full levels
  enddo  ! end loop of j
  enddo  ! end loop of 1

end subroutine reshape_mynn_array_to_am4

!###################################

subroutine reshape_mynn_array_to_am4_half (ix_dum, jx_dum, kx_dum, mynn_array_3d, am4_array_3d)

!--- input arguments 
  integer                   , intent(in)  :: ix_dum, jx_dum, kx_dum  ! not used
  real, dimension(:, :, :), intent(in) :: &    ! (i,k,j)
    mynn_array_3d

!--- output arguments
  real, dimension(:, :, :), intent(out) :: &   ! (i,j,k)
    am4_array_3d

!--- local variable
  integer i,j,k,kk
  integer ix_am4, jx_am4, kx_am4, ix_mynn, jx_mynn, kx_mynn
!----------------------------------

  ix_am4 = size(am4_array_3d,1)
  jx_am4 = size(am4_array_3d,2)
  kx_am4 = size(am4_array_3d,3)

  ix_mynn = size(mynn_array_3d,1)
  jx_mynn = size(mynn_array_3d,3)
  kx_mynn = size(mynn_array_3d,2)

  if ( ix_am4.ne.ix_mynn .or. jx_am4.ne.jx_mynn) then
    call error_mesg( ' edmf_mynn',     &
                     ' in sub reshape_mynn_array_to_am4_half, ix and jx must be the same',&
                     FATAL )
  endif

  am4_array_3d = 0.

  do i=1,ix_am4
  do j=1,jx_am4

    !--- both am4 and mynn has k+1 levels
    if (kx_am4.eq.kx_mynn) then
      do k=1,kx_am4      ! k index for half levels
        kk=kx_am4+1-k
        am4_array_3d   (i,j,kk) = mynn_array_3d (i,k,j)
      enddo  ! end loop of k, half levels

    !--- am4 has kx+1 levels but mynn has kx levels, but both indicate half levels
    elseif (kx_am4.eq.kx_mynn+1) then
      do k=1,kx_mynn      ! k index for half levels
        kk=kx_mynn+2-k
        am4_array_3d   (i,j,kk) = mynn_array_3d (i,k,j)
      enddo  ! end loop of k, half levels

      ! k at the topmost level
      am4_array_3d (i,j,1) = mynn_array_3d(i,kx_mynn,j)
    endif
  enddo  ! end loop of j
  enddo  ! end loop of 1

end subroutine reshape_mynn_array_to_am4_half

!###################################

subroutine Poisson_knuth (kx, nx, rx, rr, ENTf, ENTi)
!----------
! Desecription:
!   a Poisson random number generator attributed to Donald Knuth:
!
! Reference:
!   Wiki: https://en.wikipedia.org/wiki/Poisson_distribution#Generating_Poisson-distributed_random_variables
!         https://www.johndcook.com/blog/2010/06/14/generating-poisson-random-values/
! 
!
!---------
! algorithm - Knuth_Junhao method
!init:
!        Let λLeft ← λ, k ← 0 and p ← 1.
!    do:
!        k ← k + 1.
!        Generate uniform random number u in (0,1) and let p ← p × u.
!        while p < 1 and λLeft > 0:
!            if λLeft > STEP:
!                p ← p × eSTEP
!                λLeft ← λLeft − STEP
!            else:
!                p ← p × eλLeft
!                λLeft ← 0
!    while p > 1.
!    return k − 1.
!
!---------
! algorithm - Knuth method
!
!   init:
!           Let L ← exp(−λ), k ← 0 and p ← 1.
!      do:
!           k ← k + 1.
!           Generate uniform random number u in [0,1] and let p ← p × u.
!      while p > L.
!      return k − 1.
!
! Author: Yi-Hsuan Chen
!----------

  !--- input/output arguments
  integer, intent(in)  :: kx, nx, rx   ! dimension of input/output variables
  real   , intent(in)  :: rr   (rx)
  real   , intent(in)  :: ENTf (kx,nx)  ! Poisson parameter

  integer, intent(out) :: ENTi(kx,nx)  ! a random integer number drawn from the Poisson distribution

  !--- local variables
  integer, parameter :: itermax = 500   ! maximum of iteration loop
  real,    parameter :: step    = 70

  real ::        &
    LL,          &    ! L, lambda, p, and u in the equations in reference
    lambda,      &
    lambda_left, &
    pp,          &
    uu
  integer ::     &
    kk             ! k in the equations in reference

  integer i,j,n,k

  character*40 method
!-------------------------------
  !method = "Knuth"
  method = "Knuth_Junhao"

  !--- initialize 
  ENTi = 0
  j = 1

!--- loop over each element
do n=1,nx
do k=1,kx

!======================================
if (method == "Knuth_Junhao") then
!======================================

  !--- initialize values
  lambda = ENTf(k,n)
  ENTi(k,n) = int(lambda)
  lambda_left = lambda
  kk = 0
  pp = 1.

  !************
  do i=1,itermax
  !************


    !--- restart the random number array
    if (j > rx) then
!print*,'restart the index in rr'
      j=1
    endif

    kk = kk + 1
    uu = rr(j)
    pp = pp * uu

!print*,'kk,uu,pp,lambda_left',kk,uu,pp,lambda_left

    if (pp < 1 .and. lambda_left > 0.) then
      if (lambda_left > step) then
        pp = pp * exp(step)
        lambda_left = lambda_left - step
      else
        pp = pp * exp(lambda_left)
        lambda_left = 0.
      endif
    endif

    if (pp < 1) then
      ENTi(k,n) = kk - 1
      j=j+1
!print*,'yaya'
      exit
    endif
!
!
    !--- advance to the next index in the random number array
    j=j+1
!
  !************
  enddo  ! end loop of i
  !************

!============
end if   ! end if of method = "Knuth_Junhao"
!============
  if (i >= itermax) then
    ENTi(k,n) = int(lambda)
!    print*,'qq,itermax,',lambda
  endif

enddo  ! end loop of k
enddo  ! end loop of n

!print*,'ENTi',ENTi

end subroutine Poisson_knuth

!###################################

subroutine check_trc_rlzbility (dt, tracers, tracers_tend, &
                                tends_ratio, tends_ratio_k, rlz_ratio, rlz_tracer)

!---------------------------------------------------------------------
!  Check for tracer realizability. If tracer tendencies would
!  produce negative tracer mixing ratios, scale down tracer tendency
!  terms uniformly for this tracer throughout convective column. 
!
!  Reference: subroutine don_d_check_trc_rlzbility, src/atmos_param/donner_deep/donner_deep_k.F90
!
!  In the Tiedtke scheme, it calls a subroutine "impose_realizability" to force that there are no negative or extremely small qa, ql, or qi values. Any values of the prognostic variables which are less than qmin are reset to zero, while conserving total moisture.
!  (src/atmos_param/lscloud_driver/lscloud_driver.F90)
!
!  The minimum permissible value, qmin, is 1.e-10.
!  (src/atmos_param/physics_driver/physics_driver.F90)
!
!  Given that Tiedtke will do the realizability check for the cloud tracers, the realizability limiter in the MYNN-EDMF does the following:
!
!  tracer0: the input tracer concentration at a level (qa, ql, qi, qv, or qt)
!  tracer1: the updated tracer concentration at a level by the MYNN-EDMF
!  qmin: the minimum permissible tracer value in Tiedtke, 1.e-10 as default
!
!                         | tracer0<0   |   0<tracer0<qmin (e.g. tracer0=1.e-50)   |   tracer0 > qmin
!  ----------------------------------------------------------------------------------------------------
!  tracer1<-qmin          | STOP the    |   Limiter. The rescaled ratio            |    Limiter 
!  (e.g. tracer1=1.e-4)   | model       |   should be very very small              |
!  ------------------------             ---------------------------------------------------------------
!  -qmin < tracer1 < 0    | sth is wrong|   No Limiter. Let Tiedtke do the         |    Limiter 
!  (e.g. tracer1=-1.e-30) |             |   realizability check                    |
!  ----------------------------------------------------------------------------------------------------
!  tracer1 > 0            |             |   No Limiter                             !    No Limiter
!  ----------------------------------------------------------------------------------------------------
!
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! Arguments (Intent in)
!     dt             physics time step               , [ sec ]
!     tracers         tracer mixing ratios            , [ kg(tracer) / kg (dry air) ]
!     tracers_tend    tendency of tracer mixing ratios, [ kg(tracer) / kg (dry air) / sec ]
!---------------------------------------------------------------------
  real,    intent(in)                     :: dt
  real,    intent(in), dimension(:,:,:,:) :: tracers, tracers_tend  ! dimension (nlon, nlat, nlev, ntracer)

!---------------------------------------------------------------------
! Arguments (Intent out)
!     tends_ratio    ratio by which tracer tendencies need to 
!                    be reduced to permit realizability (i.e., to prevent
!                    negative tracer mixing ratios)
!     rlz_ratio      the smallest ratio of all tracers. This would be used to scale the all tendencies (Temperature, tracerts, etc)
!---------------------------------------------------------------------
  real   , intent(out), dimension(:,:,:)     :: tends_ratio          ! dimension (nlon, nlat, ntracer)

  real   , intent(out), dimension(:,:,:)     :: tends_ratio_k        ! dimension (nlon, nlat, ntracer)
 
  real   , intent(out), dimension(:,:)       :: rlz_ratio         ! dimension (nlon, nlat)

  real   , intent(out), dimension(:,:)       :: rlz_tracer           ! dimension (nlon, nlat)
 
!---------------------------------------------------------------------
! Arguments (Intent local)
!     tracer0        column tracer mixing ratios before MF
!     tracer1        column tracer mixing ratios after  MF transport only
!     trtend         column tracer mixing ratio tendencies due to convective transport [ (tracer units) / s ]
!     tracer_min     minimum of tracer0
!     tracer_max     maximum of tracer0
!---------------------------------------------------------------------

  real, dimension(size(tracers, 3)) ::    &  ! dimension (nlay)
    tracer0, trtend, tracer1

  real :: &
    tracer0_min, tracer0_max, ratio, ratio0

  character*40 cause

  !--- index variables & dimension
  integer i,j,k,n
  integer ix,jx,kx,nx

!--------------------------

!--- set dimensions
  ix  = size( tracers, 1 )
  jx  = size( tracers, 2 )
  kx  = size( tracers, 3 )
  nx  = size( tracers, 4 )

!--- initialze return variable
  tends_ratio   = 1.
  tends_ratio_k = 0.
  rlz_ratio     = 1.
  rlz_tracer    = 0.

!---------------------
! compute tend_ratio
!   Updated tracer concentation must be 
!   (1) not negative, 
!   (2) in the range of max/min of tracer0
!---------------------

  !do n=1,nx
  do n=1,1
  do i=1,ix
  do j=1,jx

    !--- set column tracer concentration
    tracer0(:)  = tracers(i,j,:,n)
    trtend (:)  = tracers_tend(i,j,:,n)
    tracer1(:)  = tracer0(:) + dt * trtend(:)
    !write(6,*),'tracer0',tracer0
    !write(6,*),'tracer1',tracer1
    !write(6,*),'trtend',trtend

    !--- testing values 
    !tracer1(28) = 1.
    !tracer1(28) = -1.E-25
    !tracer0(28) = 1.E-30
    !tracer0(28) = 1.E-8 

    if (any(tracer0.lt.0.)) then
      write(6,*) 'i,j,#tracer',i,j,n
      write(6,*) 'tracer0',tracer0
      call error_mesg( ' edmf_mynn',     &
                       ' check_trc_rlzbility, tracers < 0.',&
                       FATAL )
    endif

    !--- get max/min of tracer0
    tracer0_min = 1.e20
    tracer0_max = -1.e20

    do k = 1,kx
       if (trtend(k) /= 0.) then
          tracer0_max = max(tracer0(k),tracer0_max)
          tracer0_min = min(tracer0(k),tracer0_min)
       end if
    end do
    !print*,'tracer_max',tracer_max
    !print*,'tracer_min',tracer_min

    !--- compute ratio
    ratio = 1.
    do k = 1,kx

       !--- call realizability limiter in these conditions
       if (tracer0(k) > 0. .and. tracer1(k) < -qmin) then
          ratio0 = ratio
          ratio = MIN( ratio,tracer0(k)/(-trtend(k)*dt) )
          if (ratio.ne.ratio0) tends_ratio_k(i,j,n) = k          

          !write(6,*),'-------'
          !write(6,*),'aaa1, n,k,ratio',n,k,ratio
          !write(6,*),'  tracer0(k), tracer1(k), trtend(k)',tracer0(k), tracer1(k), trtend(k)

       elseif (tracer0(k) > qmin .and. tracer1(k) < 0. .and. tracer1(k) > -qmin) then  
          ratio0 = ratio
          ratio = MIN( ratio,tracer0(k)/(-trtend(k)*dt) )
          if (ratio.ne.ratio0) tends_ratio_k(i,j,n) = k          

          !write(6,*),'-------'
          !write(6,*),'aaa2, n,k,ratio',n,k,ratio
          !write(6,*),'  tracer0(k), tracer1(k), trtend(k)',tracer0(k), tracer1(k), trtend(k)

       endif

       !--- if tracer1 is less than zero
       !if (tracer0(k) > 0. .and. tracer1(k)<0.) then
!       if (tracer0(k)>0. .and. tracer1(k)<0. .and. abs(trtend(k))>qdt_min ) then   ! to avoid tracer0=1.E-100, trtend=-1.E-50. tracer1<0/ 
!
!          if (tracer0(k) .lt. tracer_min) then  ! if tracer0 is too small, set ratio to zero
!            ratio0 = 0. 
!            ratio  = 0.
!          else
!            ratio0 = ratio
!            ratio = MIN( ratio,tracer0(k)/(-trtend(k)*dt) )
!          endif
!
!            !cause = "tracer1 is less than zero"
!            !write(6,*),'-------'
!            !write(6,*),'aa1, n,k,ratio',n,k,ratio
!            !write(6,*),'  tracer0(k), tracer1(k), trtend(k)',tracer0(k), tracer1(k), trtend(k)
!            if (ratio.ne.ratio0) tends_ratio_k(i,j,n) = k          
!       end if

       !--- yhc 2021-12-13, comment this out because MF can produce tracer values less than the original minimum 
       !--- if tracer1 is less than tracer_min
       !if (tracer1(k)<tracer_min .and. trtend(k) /= 0.0 ) then
       !   ratio = MIN( ratio,(tracer0(k)-tracer_min)/(-trtend(k)*dt) )
       !   cause = "tracer1 is less than tracer_min"
       !   !write(6,*),'-------'
       !   write(6,*),'aa2, less than min, k,ratio',k,ratio
       !   write(6,*),'  tracer1(k), tracer_min, ',tracer1(k), tracer_min
       !end if

       !--- yhc 2021-12-13, comment this out because when detrainment occurs, new clouds can form so the ql>tracer_max
       !--- if tracer1 is larger than tracer_max
       !if (tracer1(k)>tracer_max  .and. trtend(k) /= 0.0 ) then
       !   ratio = MIN( ratio,(tracer_max-tracer0(k))/(trtend(k)*dt) )
       !   cause = "tracer1 is larger than tracer_max"
          !write(6,*),'-------'
          !write(6,*),'aa3, larger than max, k,ratio',k,ratio
          !write(6,*),'  tracer1(k), tracer_max, ',tracer1(k), tracer_max
       !end if

    end do

    !--- make sure 1 > ratio > 0
    ratio = MAX(0.,MIN(1.,ratio))

    tends_ratio(i,j,n) = ratio
  enddo   ! end do of j
  enddo   ! end do of i
  enddo   ! end do of n 

!----------------------------------------------
! find out the smallest ratio in all tracers, 
! and then return this ratio
!----------------------------------------------
  do i=1,ix
  do j=1,jx
  do n=1,nx
    rlz_ratio(i,j) = MAX(0.,MIN( rlz_ratio(i,j),tends_ratio(i,j,n) ) )
    if (rlz_ratio(i,j).ne.1. .and. rlz_ratio(i,j).eq.tends_ratio(i,j,n)) rlz_tracer(i,j) = float(n)   ! save the tracer name
  enddo   ! end do of n 
  enddo   ! end do of j
  enddo   ! end do of i

end subroutine check_trc_rlzbility

!#############################
! Mellor-Yamada
!#############################

end module edmf_mynn_mod
