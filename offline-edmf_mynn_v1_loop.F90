MODULE module_bl_mynn
!   took from module_bl_mynn_v1.F90

!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------

!############################
!############################

!  Don't copy these codes except edmf_input_type and edmf_output_type
!
!############################
!############################
  !character*50 :: input_profile = "SCM_am4p0_DCBL_C1_begin"
  !character*50 :: input_profile = "SCM_am4p0_DCBL_C1_01"
  !character*50 :: input_profile = "SCM_am4p0_DCBL_C1_02_u,vdt_NaN"
  !character*50 :: input_profile = "SCM_am4p0_BOMEX_01"
  !character*50 :: input_profile = "AMIP_i27_j01_IndOcn"
  !character*50 :: input_profile = "SCM_am4p0_BOMEX_02"
  !character*50 :: input_profile = "SCM_am4p0_RF01_01"
  ! character*50 :: input_profile = "SCM_am4p0_RF01_02"
  ! character*50 :: input_profile = "SCM_am4p0_RF01_03_cloudy"
  !character*50 :: input_profile = "SCM_RF01_modQDT-Gmy_aTnTtT_a3_5.0h"
  !character*50 :: input_profile = "xxx"
  !character*50 :: input_profile = "SCM_BOMEX_MYNN_ED_mixleng3"
  !character*50 :: input_profile = "SCM_RF01_mynn_EDMFexpUP_Gmy_ADD_0.5h"
  !character*50 :: input_profile = "SCM_RF01_rfo76a-M3_EDMFexpUP_NOsm01"
  character*50 :: input_profile = "SCM_RF01_rfo76a-M3_EDMFexpUP_NOsm02"

  integer, parameter :: loop_times = 1
 ! integer, parameter :: loop_times = 24 
 ! integer, parameter :: loop_times = 100
 ! integer, parameter :: loop_times = 60 

  integer, parameter :: ni = 1
  integer, parameter :: nj = 1
  integer, parameter :: npz = 33
  integer, parameter :: nfull = npz
  integer, parameter :: nhalf = nfull+1
  integer, parameter :: nfull2 = nfull+1
  integer, parameter :: nhalf2 = nfull2+1
  integer, parameter :: tr = 5 

  integer, parameter :: ntracer = 32 
  integer, parameter :: ntp = 26
  integer, parameter :: nQke = 28 
  integer, parameter :: nSh3D = 29
  integer, parameter :: nel_pbl = 30 
  integer, parameter :: ncldfra_bl = 31
  integer, parameter :: nqc_bl = ntracer

  integer, parameter :: nsphum = 1
  integer, parameter :: nql = 2
  integer, parameter :: nqi = 3
  integer, parameter :: nqa = 4
  integer, parameter :: nqn = 5

  real, parameter :: dt = 1800.
!  real, parameter :: dt = 300.
  integer :: FATAL
  logical :: do_stop_run = .false.
  !character*20 :: option_surface_flux = "star"
  character*20 :: option_surface_flux = "updated"

type time_type
  integer :: tt1, tt2
end type time_type

type physics_input_block_type
  real, dimension(ni,nj,nfull) :: &
    t,u,v,omega,p_full,z_full
  real, dimension(ni,nj,nhalf) :: &
    p_half, z_half
  real, dimension(ni,nj,nfull,tr) :: &
    q
end type physics_input_block_type

!type physics_tendency_block_type
!    real, dimension(ni,nj,nfull), target :: &
!      u_dt, v_dt, t_dt
!    real, dimension(ni,nj,nfull,ntracer), target ::  &
!      q_dt, qdiag
!end type physics_tendency_block_type

real, public, parameter :: VONKARM     = 0.40       !< Von Karman constant [dimensionless]
real, public, parameter :: rdgas  = 287.04
real, public, parameter :: rvgas = 461.50
real, public, parameter :: grav   = 9.80
real, public, parameter :: hlv      = 2.5e6     !< Latent heat of evaporation [J/kg]
real, public, parameter :: hlf      = 3.3358e5  !< Latent heat of fusion [J/kg]
real, public, parameter :: HLS = HLV + HLF          !< Latent heat of sublimation [J/kg]
real, public, parameter :: con_cliq = 4.1855e+3 !< spec heat H2O liq [J/kg/K]
real, public, parameter :: con_csol = 2.1060e+3 !< spec heat H2O ice [J/kg/K]
real, public, parameter :: kappa  = 2./7.
real, public, parameter :: cp_air   = 1004.6      !< Specific heat capacity of dry air at constant pressure [J/kg/deg]

real, public, parameter :: dens_h2o = 1000.     ! Density of liquid water [kg/m^3]
real, public, parameter :: pi = 3.14159265358979323846  ! Ratio of circle circumference to diameter [N/A]   

   real, parameter :: p00    = 1000.0e2     ! 1000hPa
   real, parameter :: p00inv = 1./p00
   real, parameter :: d622   = rdgas/rvgas
   real, parameter :: d378   = 1.-d622
   real, parameter :: d608   = d378/d622

   INTEGER :: mynn_level        = 2         ! level <>3 to do the level 2.5 MYNN model; if 3 then will do level 3
   INTEGER :: grav_settling     = 0         ! 1 - the cloud/fog will experience gravitational settling, 0 - it does not
   INTEGER :: bl_mynn_tkebudget = 0         ! if 1 the budget terms in the TKE equation are allocated for output (WRF), if 0 then not
   INTEGER :: bl_mynn_cloudpdf  = 1         ! define the type of the subgrid PDF for cloud computation, 1= Nakanishi & Niino 2004
   INTEGER :: bl_mynn_mixlength = 2         ! defines the ED mixing length formulation
   INTEGER :: bl_mynn_edmf      = 3         ! controls  the version of the EDMF to be called 3=JPLedmf documented by Wu et al., 2020
                                            ! set “bl_mynn_edmf=0” and “bl_mynn_edmf_dd<>1” then the scheme will be MYNN only.
   INTEGER :: bl_mynn_edmf_dd   = 0         ! 0 - no downdrafts, 1 - Wu et al., 2020 downdrafts 
   REAL    :: bl_mynn_edmf_Lent = 0.        ! dummy argument
   INTEGER :: bl_mynn_edmf_mom  = 1         ! 1- mixing of momentum by full EDMF, 0 - momenum is mixed only by ED
   INTEGER :: bl_mynn_edmf_tke  = 0         ! does not effect JPLedmf
   INTEGER :: bl_mynn_edmf_part = 1         ! partition area into updrafts and remaining environment
   INTEGER :: bl_mynn_cloudmix  = 1         ! 1 - cloud species are mixed separately, 0 - not
   INTEGER :: bl_mynn_mixqt     = 2         ! will mix moist conserved variables, after mixing invokes PDF cloud scheme to convert moist variables to dry
   INTEGER :: icloud_bl         = 1         ! 1, cloud cover and liquid water from the PDF scheme will be on the cloud output

   !integer :: initflag = 1                 ! (when 1 it initializes TKE using level 2 MYNN scheme) 
   integer :: initflag = 0 
   logical :: FLAG_QI  = .true.             ! (flags for whether cloud and ice mixing rations and number concentrations are mixed separately)
   logical :: FLAG_QNI = .false.            ! all false
   logical :: FLAG_QC  = .false.
   logical :: FLAG_QNC = .false.

  !logical :: expmf  = .false.
  !real    :: upwind = 0.5

  logical :: expmf = .true.              ! .true.  ... explicit mass-flux, .false. .... implicit mass-flux
  real    :: upwind = 1.                 ! upwind=1. ... use upwind approximation for mass-flux calculation
                                         ! upwind=0.5 ... use centered difference for mass-flux calculation
  real :: L0 = 100.                      ! entrainemnt rate parameter
  integer :: NUP=10                      ! the number of updrafts
  REAL :: UPSTAB=1.            ! stability parameter for massflux, (mass flux is limited so that dt/dz*a_i*w_i<UPSTAB)
  INTEGER :: bl_mynn_stabfunc  = 1       !  option for stability function. 0 - MYNN, 1 - Suselj et al. (2019)

!############################
!############################

!  copy these codes to AM4 
!
!############################
!############################

! namelist parameters
   integer :: ii_write = 0       ! i index for column written out. Set to 0 if you want to write out in SCM
   integer :: jj_write = 0       ! j index for column written out. Set to 0 if you want to write out in SCM
   !integer :: ii_write = -999       ! i index for column written out. Set to 0 if you want to write out in SCM
   !integer :: jj_write = -999       ! j index for column written out. Set to 0 if you want to write out in SCM
   real    :: lat_write = -999.99   ! latitude  (radian) for column written out
   real    :: lon_write = -999.99   ! longitude (radian) for column written out
   real    :: lat_range = 0.001
   real    :: lon_range = 0.001
   !logical :: do_writeout_column_nml = .true.
   logical :: do_writeout_column_nml = .false.
   !logical :: do_edmf_mynn_diagnostic = .true.
   logical :: do_edmf_mynn_diagnostic = .false.
   logical :: do_edmf2ls_mp = .true.
   logical :: do_qdt_same_as_qtdt = .false.

   real    :: tdt_max     = 500. ! K/day
   logical :: do_limit_tdt = .false.
   real    :: tdt_limit =   200. ! K/day

   !logical :: do_check_consrv = .true.
   logical :: do_check_consrv = .false.

   !logical :: do_check_realizability = .true.
   logical :: do_check_realizability = .false.

  integer :: edmf_type=2                        ! =0, the standard MYNN code, in which the PDF cloud scheme before mixing and after the mixing and compute the tendencies of liquid and cloud properties from the differences between these two.
                                                ! =1, tendencies of moist variables from the PDF scheme after mixing and from the input values (from Tiedtke, presumably)
  real    :: qke_min = -1.                      ! qke=2*tke. If qke < qke_min, set all EDMF tendencies to zeros
                                                !   set qke_min>0 may remove energy/water away and cause water mass is not conserved
  real    :: tracer_min = 1.E-10                ! make sure tracer value is not smaller than tracer_min
                                                ! 1.E-10 is same as qmin in lscloud_driver
  integer :: do_option_edmf2ls_mp = 0           ! option to include EDMF cloud tendencies terms into Tiedtke
                                                ! =0, add EDMF terms to accumulated ql and qi tendnecies
                                                ! =1, add EDMF terms to Tiedtke and keep Tiedtke terms except turbulence heating
                                                ! =2, add EDMF terms to Tiedtke and remove large-scale and coud erosion terms
                                                ! =3, evaporate ql and qi and put these water back to vapor and change temperature accordingly
  logical :: do_use_tau = .true.                ! .true.  : use the T,q at the current step
                                                ! .false. : use the updated T,q
  logical :: do_simple =.false.                 ! do_simple = switch to turn on alternative definition of specific
                                                !             humidity. When true, specific humidity =
                                                !             (rdgas/rvgas)*esat/pressure
                                                ! same setting as module moist_processes_mod
  character*20 :: do_debug_option = ""          ! debug purpose

logical :: do_return_edmf_mynn_diff_only = .false. ! .true.,  return edmf_mynn diffusion coefficients 
                                                   !          and let vert_diff do the diffusion rather than in edmf_mynn
                                                   ! .false., tendencies are updated in edmf_mynn
character*5 :: do_edmf_mynn_in_physics = "up"     ! where to call edmf_mynn. "up" in physics_driver_up; "down" in physics_driver_down

  logical :: do_pblh_constant = .false.    ! fix PBL depth for testing
  real    :: fixed_pblh       = 2500.

  real    :: sgm_factor = 100.                  ! factor in computing sigma_s in MYNN

  integer :: Qx_MF = 1 

  real    :: rc_MF = 10.e-6                     ! assumed cloud droplet radius in plumes (meters)

  integer :: do_tracers_selective=6

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
    dz, u, v, w, th, qv, p, exner, rho, T3D, qc, ql, qi, qnc, qni, qa

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
    edmf_a, edmf_w, edmf_qt, edmf_thl, edmf_ent, edmf_qc 

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


!############################
!############################


!############################
!############################
!
! module_by_mynn from Kay
!   module_bl_mynn_v1.F90
!############################
!############################


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
!         edmf_qt .... qt in updrafts  [kg kg^-1]
!         edmf_thl .... thl in updrafts  [K]
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

!-------------------------------------------------------------------
!  IMPLICIT NONE
!-------------------------------------------------------------------
!For non-WRF
   REAL    , PARAMETER :: karman       = 0.4
   REAL    , PARAMETER :: g            = 9.81
   REAL    , PARAMETER :: r_d          = 287.
   REAL    , PARAMETER :: cp           = 7.*r_d/2.
   REAL    , PARAMETER :: r_v          = 461.6
   REAL    , PARAMETER :: cpv          = 4.*r_v
   REAL    , PARAMETER :: cliq         = 4190.
   REAL    , PARAMETER :: Cice         = 2106.
   REAL    , PARAMETER :: rcp          = r_d/cp
   REAL    , PARAMETER :: XLS          = 2.85E6
   REAL    , PARAMETER :: XLV          = 2.5E6
   REAL    , PARAMETER :: XLF          = 3.50E5
   REAL    , PARAMETER :: p1000mb      = 100000.
   REAL    , PARAMETER :: rvovrd       = r_v/r_d
   REAL    , PARAMETER :: SVP1         = 0.6112
   REAL    , PARAMETER :: SVP2         = 17.67
   REAL    , PARAMETER :: SVP3         = 29.65
   REAL    , PARAMETER :: SVPT0        = 273.15
   REAL    , PARAMETER :: EP_1         = R_v/R_d-1.
   REAL    , PARAMETER :: EP_2         = R_d/R_v
! The following depends on the microphysics scheme used:
   INTEGER , PARAMETER :: param_first_scalar = 1
   INTEGER , PARAMETER :: p_qc = 2
   INTEGER , PARAMETER :: p_qr = 0
   INTEGER , PARAMETER :: p_qi = 2
   INTEGER , PARAMETER :: p_qs = 0
   INTEGER , PARAMETER :: p_qg = 0
   INTEGER , PARAMETER :: p_qnc= 0
   INTEGER , PARAMETER :: p_qni= 0
   REAL    , PARAMETER :: rhoair0      = 1.28        ! density of dry air at 0^oC and 1000mb pressure (kg m^-3)

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

! Constants for cloud PDF (mym_condensation; 1/sqrt(2.) and  sqrt(1./(2*pi)))
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

!JOE & JAYMES'S mods
!
! Mixing Length Options 
!   specifed through namelist:  bl_mynn_mixlength
!   added:  16 Apr 2015
!
! 0: Uses original MYNN mixing length formulation (except elt is calculated from 
!    a 10-km vertical integration).  No scale-awareness is applied to the master
!    mixing length (el), regardless of "scaleaware" setting. 
!
! 1 (*DEFAULT*): Instead of (0), uses BouLac mixing length in free atmosphere.  
!    This helps remove excessively large mixing in unstable layers aloft.  Scale-
!    awareness in dx is available via the "scaleaware" setting.  As of Apr 2015, 
!    this mixing length formulation option is used in the ESRL RAP/HRRR configuration.
!
! 2: As in (1), but elb is lengthened using separate cloud mixing length functions 
!    for statically stable and unstable regimes.  This elb adjustment is only 
!    possible for nonzero cloud fractions, such that cloud-free cells are treated 
!    as in (1), but BouLac calculation is used more sparingly - when elb > 500 m. 
!    This is to reduce the computational expense that comes with the BouLac calculation.
!    Also, This option is  scale-aware in dx if "scaleaware" = 1. (Following Ito et al. 2015). 
!
!JOE & JAYMES- end



!  INTEGER :: mynn_level

  CHARACTER*128 :: mynn_message

  INTEGER, PARAMETER :: kdebug=27

CONTAINS

!###################################
!###################################
!
! yhc, copy codes below to AM4 
!
!###################################
!###################################

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
    REAL, DIMENSION(IMS:IME,KMS:KME+1,JMS:JME), INTENT(out) :: &
      a_moist_half , &
      a_dry_half   , &
      mf_moist_half, &
      mf_dry_half  , &
      mf_all_half  , &
      qv_moist_half, &
      qv_dry_half

    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(out) :: &
      a_moist_full , &
      a_dry_full , &
      mf_moist_full,  &
      mf_dry_full  ,  &
      mf_all_full  , &
      qv_moist_full,  &
      qv_dry_full

    REAL,DIMENSION(KTS:KTE+1) :: &
      a_moist_half1 , &
      a_dry_half1   , &
      mf_moist_half1, &
      mf_dry_half1  , &
      mf_all_half1  , &
      qv_moist_half1, &
      qv_dry_half1

    REAL,DIMENSION(KTS:KTE) :: &
      a_moist_full1 ,  &
      a_dry_full1 ,  &
      mf_moist_full1,  &
      mf_dry_full1  ,  &
      mf_all_full1  , &
      qv_moist_full1,  &
      qv_dry_full1
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


    !<--- yhc_mynn, add new output variables, 2021-04-02
    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(out) :: &
    !REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: &
       qa_before_mix, ql_before_mix, qi_before_mix, thl_before_mix, qt_before_mix, th_before_mix, &
       qa_before_pdf, ql_before_pdf, qi_before_pdf, &
       qa_after_mix, ql_after_mix, qi_after_mix, thl_after_mix, qt_after_mix, th_after_mix

    !REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: &
    !   qa_before_pdf, ql_before_pdf, qi_before_pdf

    REAL, DIMENSION(KMS:KME) :: &
       dum_1D
    !---> yhc_mynn

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
               &u1,v1,th1,thl,thetav,tk1,sqw,sqv,sqc, &
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
               &ktop_shallow(i,j),ztop_shallow,   &
               & KPBL(i,j),                        &
               & Q_ql1,Q_qi1,Q_a1,                 &
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
              &u,v,th,thl,thv,tk,qt,qv,qc,   &
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
              &ktop,ztop,kpbl,Qql,Qqi,Qa, &
              & Qql_adv,Qqi_adv,Qa_adv, Qql_eddy,Qqi_eddy,Qa_eddy, Qql_ent,Qqi_ent,Qa_ent, Qql_det,Qqi_det,Qa_det, Qql_sub,Qqi_sub,Qa_sub &
              ) ! yhc 2021-09-08



        INTEGER, INTENT(IN) :: KTS,KTE, kpbl
        REAL,DIMENSION(KTS:KTE), INTENT(IN) :: U,V,TH,THL,TK,QT,QV,QC
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
        REAL,DIMENSION(KTS:KTE+1,1:NUP) :: UPW,UPTHL,UPQT,UPQC,UPA,UPU,UPV,UPTHV

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
        REAL :: B,QTn,THLn,THVn,QCn,Un,Vn,Wn2,EntEXP,EntW, deltaZ,EntExp_M, Beta_un, Z00, Z_i

    ! VARIABLES FOR CHABOUREAU-BECHTOLD CLOUD FRACTION
        ! REAL,DIMENSION(KTS:KTE), INTENT(INOUT) :: vt, vq, sgm
        REAL :: sigq,xl,tlk,qsat_tl,rsl,cpm,a,qmq,mf_cf,diffqt,&
               Fng,qww,alpha,beta,bb,f,pt,t,q2p,b9,satvp,rhgrid
        
        REAL,DIMENSION(kts:kte) :: ph,lfh 
        REAL :: THVsrfF,QTsrfF,maxS,stabF,dz,F0,F1,F2,F3,CCp1,CCp0,mf,mfp1  

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
         a_moist_half , &   
         a_dry_half   , &        
         mf_moist_half, &     
         mf_dry_half  , &
         mf_all_half  , & 
         qv_moist_half, &     
         qv_dry_half  

       REAL,DIMENSION(kts:kte), INTENT(OUT) :: & 
         a_moist_full ,  &
         a_dry_full ,  &
         mf_moist_full,  &
         mf_dry_full  ,  & 
         mf_all_full  , & 
         qv_moist_full,  & 
         qv_dry_full   
  
       REAL,DIMENSION(KTS:KTE) :: ZFULL
    !---> yhc 2021-09-08 

! stability parameter for massflux
! (mass flux is limited so that dt/dz*a_i*w_i<UPSTAB)
!       REAL,PARAMETER :: UPSTAB=1.   ! yhc - move UPSTAB to namelist parameter

!Initialize values:
ktop = 0
ztop = 0.0

UPW=0.
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

 edmf_det = 0.
 DET=0.

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

IF ( wthv >= 0.0 ) then

! get the pressure and liquid_fraction on updraft levels

  lfH(KTS)=liquid_frac(KTS)
  pH(KTS)=ps 
  DO K=KTS+1,KTE
   lfH(K)=0.5*(liquid_frac(K)+liquid_frac(K-1))
    ph(K)=0.5*(p(K)+p(K-1)) 
  ENDDO




 ! set initial conditions for updrafts
    z0=50.
    pwmin=1.
    pwmax=3.

    wstar=max(0.,(g/thv(KTS)*wthv*pblh)**(1./3.))
    qstar=wqt/wstar
    thstar=wthl/wstar
    sigmaW=1.34*wstar*(z0/pblh)**(1./3.)*(1-0.8*z0/pblh)
    sigmaQT=1.34*qstar*(z0/pblh)**(-1./3.)
    sigmaTH=1.34*thstar*(z0/pblh)**(-1./3.)


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
    call Poisson(kts,kte,1,Nup,ENTf,ENTi,seedmf)


 ! entrainent: Ent=Ent0/dz*P(dz/L0)
    do i=1,Nup
        do k=kts,kte
          ENT(k,i)=real(ENTi(k,i))*Ent0/(ZW(k+1)-ZW(k))
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
            IF (Wn2 >0.) THEN
                UPW(K,I)=sqrt(Wn2)
                UPTHV(K,I)=THVn
                UPTHL(K,I)=THLn
                UPQT(K,I)=QTn
                UPQC(K,I)=QCn
                UPU(K,I)=Un
                UPV(K,I)=Vn
                UPA(K,I)=UPA(K-1,I)

                ktop = MAX(ktop,k)
            ELSE
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

    mfp1=rho(K)*UPW(K+1,I)
    if (k.eq.1) then
      mf=rho(1)*UPW(K,I)
    else
      mf=rho(K-1)*UPW(K,I)
    endif

    if (mf.gt.0) then
      DET(K,I) = ENT(K,I) - (mfp1-mf)/mf/dz   
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
     dz=ZW(k+1)-ZW(k)
   
     mfp1=UPA(K+1,I)*UPW(K+1,I)
     mf=UPA(K,I)*UPW(K,I)
   
     ! for liquid and ice water
     !F1=-(mfp1*(UPQC(K+1,I)-qc(K+1))-mf*(UPQC(K,I)-qc(K)))/dz
     if (K.eq.1) then
       F1=-(mfp1*(UPQC(K+1,I)-qc(K))-mf*(UPQC(K,I)-0.))/dz     ! qc(K-1)=qc(0)=0.
     else
       F1=-(mfp1*(UPQC(K+1,I)-qc(K))-mf*(UPQC(K,I)-qc(K-1)))/dz
     endif

     F2=mf*(UPQC(K+1,I)-UPQC(K,I))/dz  !+ENT(K,I)*mf*(UPQC(K,I)-qc(K))
     F3=ENT(K,I)*mf*(UPQC(K,I)-qc(K))
       
     ! yhc: F1 and F2 do not need to be divided by rho 
     !Qql(k)=Qql(k)+liquid_frac(k)*(F1+F2+F3)
     !Qqi(k)=Qqi(k)+(1.-liquid_frac(k))*(F1+F2+F3)

     !--- eddy-flux divergence/source form
     Qql_eddy(k)=Qql_eddy(k)+liquid_frac(k)*F1
     Qql_adv (k)=Qql_adv(k)+liquid_frac(k)*F2
     Qql_ent (k)=Qql_ent(k)+liquid_frac(k)*F3

     Qqi_eddy(k)=Qqi_eddy(k)+(1.-liquid_frac(k))*F1
     Qqi_adv (k)=Qqi_adv(k)+(1.-liquid_frac(k))*F2
     Qqi_ent (k)=Qqi_ent(k)+(1.-liquid_frac(k))*F3

     !--- detrainment/subsidence form
     if (k.eq.1) then
       mf=rho(1)*UPW(K,I)*UPA(K,I)
     else
       mf=rho(K-1)*UPW(K,I)*UPA(K,I)
     endif

     mfp1=rho(K)*UPA(K+1,I)*UPW(K+1,I)

     F1=mf*DET(K,I)*(UPQC(K,I)-qc(K))
     F2=mfp1*(qc(K+1)-qc(K))/(ZFULL(K+1)-ZFULL(K))

     Qql_det(k) = Qql_det(k)+liquid_frac(k)*F1/rho(K)
     Qqi_det(k) = Qqi_det(k)+(1.-liquid_frac(k))*F1/rho(K)

     Qql_sub(k) = Qql_sub(k)+liquid_frac(k)*F2/rho(K)
     Qqi_sub(k) = Qqi_sub(k)+(1.-liquid_frac(k))*F2/rho(K)

     !Qql(k)=Qql(k)+liquid_frac(k)*(F1+F2)/rho(k)
     !Qqi(k)=Qqi(k)+(1.-liquid_frac(k))*(F1+F2)/rho(k)

!<--- yhc 2021-09-08
!<--- yhc: this Qa recovery code is not correct  
!  CCp1=0.
!  CCp0=0.
!  IF (UPQC(K+1,I) > 0.) CCp1=UPA(K+1,I)
!    IF (UPQC(K,I) > 0.) CCp0=UPA(K,I)
!  
!    ! for area fraction
!     F1=-(mfp1*(CCp1-cldfra_bl1d(K+1))-mf*(CCp0-cldfra_bl1d(K)) )/dz
!     F2=mf*(CCp1-CCp0)/dz+ENT(K,I)*mf*(CCp0-cldfra_bl1d(K))
!  
!     Qa(k)=Qa(k)+(F1+F2)/rho(k)
!--->  
  
  CCp1=0.
  CCp0=0.
  IF (UPQC(K+1,I) > 0.) THEN
    CCp1=1.
  ELSE
    CCp1=0.
  ENDIF

  IF (UPQC(K,I) > 0.) THEN
    CCp0=1.
  ELSE
    CCp0=0.
  ENDIF

  IF (UPQC(K,I) <= 0. .and. UPQC(K+1,I) > 0.) then
    F2 = UPA(K,I)*UPW(K,I)/dz
  ELSEIF (UPQC(K,I) > 0. .and. UPQC(K+1,I) <= 0.) then
    F2 = -1.*UPA(K,I)*UPW(K,I)/dz
  ELSE
    F2 = 0.
  ENDIF
  
    ! for area fraction
     if (K.eq.1) then
       F1=-(mfp1*(CCp1-cldfra_bl1d(K))-mf*(CCp0-0.) )/dz   ! cldfra_bl1d(0) = 0.
     else
       F1=-(mfp1*(CCp1-cldfra_bl1d(K))-mf*(CCp0-cldfra_bl1d(K-1)) )/dz
     endif

     !F2=F0 !+ ENT(K,I)*mf*(CCp0-cldfra_bl1d(K))
     F3=ENT(K,I)*mf*(CCp0-cldfra_bl1d(K))
  
     !Qa(k)=Qa(k)+(F1+F2)/rho(k)
     !Qa(k)=Qa(k)+(F1+F2+F3)

     !--- eddy-flux divergence/source form
     Qa_eddy(k)=Qa_eddy(k)+F1
     Qa_adv (k)=Qa_adv(k)+F2
     Qa_ent (k)=Qa_ent(k)+F3

     !--- detrainment/subsidence form
     if (k.eq.1) then
       mf=rho(1)*UPW(K,I)*UPA(K,I)
     else
       mf=rho(K-1)*UPW(K,I)*UPA(K,I)
     endif

     mfp1=rho(K)*UPA(K+1,I)*UPW(K+1,I)

     F1=mf*DET(K,I)*(CCp0-cldfra_bl1d(K))
     F2=mfp1*(cldfra_bl1d(K+1)-cldfra_bl1d(K))/(ZFULL(K+1)-ZFULL(K))

     Qa_det(k) = Qa_det(k)+F1/rho(K)
     Qa_sub(k) = Qa_sub(k)+F2/rho(K)
 
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

!print*,'Qx_MF',Qx_MF
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
  IF(k > KTOP) exit

  DO I=1,NUP
    IF(I > NUP) exit

    if (UPQC(K,I) > 0.) then   ! sum of individual moist updrafts
      a_moist_half  (K) = a_moist_half  (K) + UPA(K,I) 
      mf_moist_half (K) = mf_moist_half (K) + UPA(K,I)*rho(K)*UPW(K,I)
      qv_moist_half (K) = qv_moist_half (K) + UPA(K,I)*(UPQT(K,I)-UPQC(K,I))  

    else                       ! sum of individual dry updrafts
      a_dry_half    (K) = a_dry_half  (K) + UPA(K,I) 
      mf_dry_half   (K) = mf_dry_half (K) + UPA(K,I)*rho(K)*UPW(K,I)
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
  a_moist_full  (KTOP) = a_moist_half  (KTOP)
  mf_moist_full (KTOP) = mf_moist_half (KTOP)
  qv_moist_full (KTOP) = qv_moist_half (KTOP)

  a_dry_full    (KTOP) = a_dry_half    (KTOP)
  mf_dry_full   (KTOP) = mf_dry_half   (KTOP)
  qv_dry_full   (KTOP) = qv_dry_half   (KTOP)

DO K=KTS,KTOP-1
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
  real, intent(inout), dimension(:,:,:,:) :: &   ! Mellor-Yamada, use this in offline mode
  !real, intent(inout), dimension(:,:,:,ntp+1:) :: &
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
  real    :: tt1
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
!-------------------------
  ix = size(Physics_input_block%t,1)
  jx = size(Physics_input_block%t,2)
  kx = size(Physics_input_block%t,3)  

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
  call modify_mynn_edmf_tendencies( is, ie, js, je, Time_next,  &
                                    do_writeout_column,         &
                                    Physics_input_block, Input_edmf, rdt_mynn_ed_am4, &
                                    size(Physics_input_block%t,1), size(Physics_input_block%t,2), size(Physics_input_block%t,3), &
                                    Output_edmf)

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

!!send_data
!!------- zonal wind stress (units: kg/m/s2) at one level -------
!      if ( id_u_flux > 0) then
!        used = send_data (id_u_flux, u_flux, Time_next, is, js )
!      endif
!
!!------- meridional wind stress (units: kg/m/s2) at one level -------
!      if ( id_v_flux > 0) then
!        used = send_data (id_v_flux, v_flux, Time_next, is, js )
!      endif
!
!!------- u_star from u_flux and v_flux (units: m/s) at one level -------
!      if ( id_u_star_updated > 0) then
!        used = send_data (id_u_star_updated, Input_edmf%u_star_updated, Time_next, is, js )
!      endif
!
!!------- sensible heat flux from star (units: W/m2) at one level -------
!      if ( id_shflx_star > 0) then
!        used = send_data (id_shflx_star, Input_edmf%shflx_star, Time_next, is, js )
!      endif
!
!!------- evaporation flux from star (units: kg/m2/s) at one level -------
!      if ( id_lhflx_star > 0) then
!        used = send_data (id_lhflx_star, Input_edmf%lhflx_star, Time_next, is, js )
!      endif
!
!!------- kinematic virtual temperature flux from star (units: K m/s) at one level -------
!      if ( id_w1_thv1_surf_star > 0) then
!        used = send_data (id_w1_thv1_surf_star, Input_edmf%w1_thv1_surf_star, Time_next, is, js )
!      endif
!
!!------- kinematic virtual temperature flux from updated fluxes (units: K m/s) at one level -------
!      if ( id_w1_thv1_surf_updated > 0) then
!        used = send_data (id_w1_thv1_surf_updated, Input_edmf%w1_th1_surf_updated, Time_next, is, js )
!      endif
!
!!------- Obukhov length from star (units: m) at one level -------
!      if ( id_Obukhov_length_star > 0) then
!        used = send_data (id_Obukhov_length_star, Input_edmf%Obukhov_length_star, Time_next, is, js )
!      endif
!
!!------- Obukhov length from updated fluxes (units: m) at one level -------
!      if ( id_Obukhov_length_updated > 0) then
!        used = send_data (id_Obukhov_length_updated, Input_edmf%Obukhov_length_updated, Time_next, is, js )
!      endif
!
!!------- turbulent kinetic energy (units: m2/s2) at full level -------
!      if ( id_tke_edmf > 0) then
!        used = send_data (id_tke_edmf, am4_Output_edmf%tke, Time_next, is, js, 1 )
!      endif
!
!!------- variance of theta_l (units: K^2) at full level -------
!      if ( id_Tsq > 0) then
!        used = send_data (id_Tsq, am4_Output_edmf%Tsq, Time_next, is, js, 1 )
!      endif
!
!!------- covariance of theta_l and q_t (units: none) at full level -------
!      if ( id_Cov_thl_qt > 0) then
!        used = send_data (id_Cov_thl_qt, am4_Output_edmf%Cov_thl_qt, Time_next, is, js, 1 )
!      endif
!
!!------- u tendency from edmf_mynn (units: m/s2) at full level -------
!      if ( id_udt_edmf > 0) then
!        used = send_data (id_udt_edmf, am4_Output_edmf%udt_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- v tendency from edmf_mynn (units: m/s2) at full level -------
!      if ( id_vdt_edmf > 0) then
!        used = send_data (id_vdt_edmf, am4_Output_edmf%vdt_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- t tendency from edmf_mynn (units: K/s) at full level -------
!      if ( id_tdt_edmf > 0) then
!        used = send_data (id_tdt_edmf, am4_Output_edmf%tdt_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- q tendency from edmf_mynn (units: kg/kg/s) at full level -------
!      if ( id_qdt_edmf > 0) then
!        used = send_data (id_qdt_edmf, am4_Output_edmf%qdt_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- qi tendency from edmf_mynn (units: kg/kg/s) at full level -------
!      if ( id_qidt_edmf > 0) then
!        used = send_data (id_qidt_edmf, am4_Output_edmf%qidt_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- qc tendency from edmf_mynn (units: kg/kg/s) at full level -------
!      if ( id_qldt_edmf > 0) then
!        used = send_data (id_qldt_edmf, am4_Output_edmf%qldt_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- cldfrac tendency from edmf_mynn (units: 1/s) at full level -------
!      if ( id_qadt_edmf > 0) then
!        used = send_data (id_qadt_edmf, am4_Output_edmf%qadt_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- cloud droplet number tendency from edmf_mynn (units: 1/kg/s) at full level -------
!      if ( id_qndt_edmf > 0) then
!        used = send_data (id_qndt_edmf, am4_Output_edmf%qndt_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- updraft area (units: none) at half level -------
!      if ( id_edmf_a > 0) then
!        call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_a, diag_half)
!        used = send_data (id_edmf_a, diag_half, Time_next, is, js, 1 )
!      endif
!
!!------- vertical velocity of updrafts (units: m/s) at half level -------
!      if ( id_edmf_w > 0) then
!        call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_w, diag_half)
!        used = send_data (id_edmf_w, diag_half, Time_next, is, js, 1 )
!      endif
!
!!------- qt in updrafts (units: kg/kg) at half level -------
!      if ( id_edmf_qt > 0) then
!        call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_qt, diag_half)
!        used = send_data (id_edmf_qt, diag_half, Time_next, is, js, 1 )
!      endif
!
!!------- thl in updrafts (units: K) at half level -------
!      if ( id_edmf_thl > 0) then
!        call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_thl, diag_half)
!        used = send_data (id_edmf_thl, diag_half, Time_next, is, js, 1 )
!      endif
!
!!------- entrainment in updrafts (units: 1/m) at full level -------
!      if ( id_edmf_ent > 0) then
!        call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%edmf_ent, diag_full)
!        used = send_data (id_edmf_ent, diag_full, Time_next, is, js, 1 )
!      endif
!
!!------- dentrainment in updrafts (units: 1/m) at full level -------
!      if ( id_edmf_det > 0) then
!        call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%edmf_det, diag_full)
!        used = send_data (id_edmf_det, diag_full, Time_next, is, js, 1 )
!      endif
!
!!------- qc in updrafts (units: kg/kg) at half level -------
!      if ( id_edmf_qc > 0) then
!        call reshape_mynn_array_to_am4_half(ix, jx, kx, Output_edmf%edmf_qc, diag_half)
!        used = send_data (id_edmf_qc, diag_half, Time_next, is, js, 1 )
!      endif
!
!!------- theta_li in edmf_mynn (units: K) at full level -------
!      if ( id_thl_edmf > 0) then
!        used = send_data (id_thl_edmf, am4_Output_edmf%thl_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- qt in edmf_mynn (units: kg/kg) at full level -------
!      if ( id_qt_edmf > 0) then
!        used = send_data (id_qt_edmf, am4_Output_edmf%qt_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- cloud fraction in edmf_mynn (units: fraction) at full level -------
!      if ( id_cldfra_bl > 0) then
!        used = send_data (id_cldfra_bl, am4_Output_edmf%cldfra_bl, Time_next, is, js, 1 )
!      endif
!
!!------- liquid water mixing ratio in edmf_mynn (units: mynn) at full level -------
!      if ( id_qc_bl > 0) then
!        used = send_data (id_qc_bl, am4_Output_edmf%qc_bl, Time_next, is, js, 1 )
!      endif
!
!!------- depth of planetary boundary layer (units: m) at one level -------
!      if ( id_z_pbl > 0) then
!        used = send_data (id_z_pbl, pbltop, Time_next, is, js )
!      endif
!
!!------- PBL depth from edmf_mynn (units: m) at one level -------
!      if ( id_z_pbl_edmf > 0) then
!        used = send_data (id_z_pbl_edmf, am4_Output_edmf%pbltop, Time_next, is, js )
!      endif
!
!!------- qt tendency from edmf_mynn (units: kg/kg/s) at full level -------
!      if ( id_qtdt_edmf > 0) then
!        used = send_data (id_qtdt_edmf, am4_Output_edmf%qtdt_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- thl tendency from edmf_mynn (units: K/s) at full level -------
!      if ( id_thldt_edmf > 0) then
!        used = send_data (id_thldt_edmf, am4_Output_edmf%thldt_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- heat diff coeffs from edmf_mynn (units: K/m/s) at half level -------
!      if ( id_diff_t_edmf > 0) then
!        ! the dimension size of am4_Output_edmf%diff_t_edmf is kx, but it is on half level. Write out on half levels
!        diag_half(:,:,kx+1) = 0.
!        diag_half(:,:,1:kx) = am4_Output_edmf%diff_t_edmf(:,:,1:kx)
!        used = send_data (id_diff_t_edmf, diag_half, Time_next, is, js, 1, mask=lmask_half )
!        !used = send_data (id_diff_t_edmf, am4_Output_edmf%diff_t_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- momentum diff coeffs from edmf_mynn (units: m2/s) at half level -------
!      if ( id_diff_m_edmf > 0) then
!        ! the dimension size of am4_Output_edmf%diff_m_edmf is kx, but it is on half level. Write out on half levels
!        diag_half(:,:,kx+1) = 0.
!        diag_half(:,:,1:kx) = am4_Output_edmf%diff_m_edmf(:,:,1:kx)
!        used = send_data (id_diff_m_edmf, diag_half, Time_next, is, js, 1, mask=lmask_half )
!        !used = send_data (id_diff_m_edmf, am4_Output_edmf%diff_m_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- mixing length in edmf_mynn (units: m) at full level -------
!      if ( id_el_edmf > 0) then
!        used = send_data (id_el_edmf, am4_Output_edmf%el_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- T input to edmf_mynn (units: K) at full level -------
!      if ( id_t_input > 0) then
!        used = send_data (id_t_input, am4_Output_edmf%t_input, Time_next, is, js, 1 )
!      endif
!
!!------- q input to edmf_mynn (units: kg/kg) at full level -------
!      if ( id_q_input > 0) then
!        used = send_data (id_q_input, am4_Output_edmf%q_input, Time_next, is, js, 1 )
!      endif
!
!!------- qa input to edmf_mynn (units: none) at full level -------
!      if ( id_qa_input > 0) then
!        used = send_data (id_qa_input, am4_Output_edmf%qa_input, Time_next, is, js, 1 )
!      endif
!
!!------- ql input to edmf_mynn (units: kg/kg) at full level -------
!      if ( id_ql_input > 0) then
!        used = send_data (id_ql_input, am4_Output_edmf%ql_input, Time_next, is, js, 1 )
!      endif
!
!!------- qi input to edmf_mynn (units: kg/kg) at full level -------
!      if ( id_qi_input > 0) then
!        used = send_data (id_qi_input, am4_Output_edmf%qi_input, Time_next, is, js, 1 )
!      endif
!
!!------- thl input to edmf_mynn (units: K) at full level -------
!      if ( id_thl_input > 0) then
!        used = send_data (id_thl_input, am4_Output_edmf%thl_input, Time_next, is, js, 1 )
!      endif
!
!!------- qt input to edmf_mynn (units: kg/kg) at full level -------
!      if ( id_qt_input > 0) then
!        used = send_data (id_qt_input, am4_Output_edmf%qt_input, Time_next, is, js, 1 )
!      endif
!
!!------- rh input to edmf_mynn (units: percent) at full level -------
!      if ( id_rh_input > 0) then
!        used = send_data (id_rh_input, am4_Output_edmf%rh_input, Time_next, is, js, 1 )
!      endif
!
!!------- theta input to edmf_mynn (units: K) at full level -------
!      if ( id_th_input > 0) then
!        used = send_data (id_th_input, am4_Output_edmf%th_input, Time_next, is, js, 1 )
!      endif
!
!!------- T before_mix to edmf_mynn (units: K) at full level -------
!      if ( id_t_before_mix > 0) then
!        used = send_data (id_t_before_mix, am4_Output_edmf%t_before_mix, Time_next, is, js, 1 )
!      endif
!
!!------- q before_mix to edmf_mynn (units: kg/kg) at full level -------
!      if ( id_q_before_mix > 0) then
!        used = send_data (id_q_before_mix, am4_Output_edmf%q_before_mix, Time_next, is, js, 1 )
!      endif
!
!!------- qa before_mix to edmf_mynn (units: none) at full level -------
!      if ( id_qa_before_mix > 0) then
!        used = send_data (id_qa_before_mix, am4_Output_edmf%qa_before_mix, Time_next, is, js, 1 )
!      endif
!
!!------- ql before_mix to edmf_mynn (units: kg/kg) at full level -------
!      if ( id_ql_before_mix > 0) then
!        used = send_data (id_ql_before_mix, am4_Output_edmf%ql_before_mix, Time_next, is, js, 1 )
!      endif
!
!!------- qi before_mix to edmf_mynn (units: kg/kg) at full level -------
!      if ( id_qi_before_mix > 0) then
!        used = send_data (id_qi_before_mix, am4_Output_edmf%qi_before_mix, Time_next, is, js, 1 )
!      endif
!
!!------- thl before_mix to edmf_mynn (units: K) at full level -------
!      if ( id_thl_before_mix > 0) then
!        used = send_data (id_thl_before_mix, am4_Output_edmf%thl_before_mix, Time_next, is, js, 1 )
!      endif
!
!!------- qt before_mix to edmf_mynn (units: kg/kg) at full level -------
!      if ( id_qt_before_mix > 0) then
!        used = send_data (id_qt_before_mix, am4_Output_edmf%qt_before_mix, Time_next, is, js, 1 )
!      endif
!
!!------- rh before_mix to edmf_mynn (units: percent) at full level -------
!      if ( id_rh_before_mix > 0) then
!        used = send_data (id_rh_before_mix, am4_Output_edmf%rh_before_mix, Time_next, is, js, 1 )
!      endif
!
!!------- theta before_mix to edmf_mynn (units: K) at full level -------
!      if ( id_th_before_mix > 0) then
!        used = send_data (id_th_before_mix, am4_Output_edmf%th_before_mix, Time_next, is, js, 1 )
!      endif
!
!!------- T after_mix to edmf_mynn (units: K) at full level -------
!      if ( id_t_after_mix > 0) then
!        used = send_data (id_t_after_mix, am4_Output_edmf%t_after_mix, Time_next, is, js, 1 )
!      endif
!
!!------- q after_mix to edmf_mynn (units: kg/kg) at full level -------
!      if ( id_q_after_mix > 0) then
!        used = send_data (id_q_after_mix, am4_Output_edmf%q_after_mix, Time_next, is, js, 1 )
!      endif
!
!!------- qa after_mix to edmf_mynn (units: none) at full level -------
!      if ( id_qa_after_mix > 0) then
!        used = send_data (id_qa_after_mix, am4_Output_edmf%qa_after_mix, Time_next, is, js, 1 )
!      endif
!
!!------- ql after_mix to edmf_mynn (units: kg/kg) at full level -------
!      if ( id_ql_after_mix > 0) then
!        used = send_data (id_ql_after_mix, am4_Output_edmf%ql_after_mix, Time_next, is, js, 1 )
!      endif
!
!!------- qi after_mix to edmf_mynn (units: kg/kg) at full level -------
!      if ( id_qi_after_mix > 0) then
!        used = send_data (id_qi_after_mix, am4_Output_edmf%qi_after_mix, Time_next, is, js, 1 )
!      endif
!
!!------- thl after_mix to edmf_mynn (units: K) at full level -------
!      if ( id_thl_after_mix > 0) then
!        used = send_data (id_thl_after_mix, am4_Output_edmf%thl_after_mix, Time_next, is, js, 1 )
!      endif
!
!!------- qt after_mix to edmf_mynn (units: kg/kg) at full level -------
!      if ( id_qt_after_mix > 0) then
!        used = send_data (id_qt_after_mix, am4_Output_edmf%qt_after_mix, Time_next, is, js, 1 )
!      endif
!
!!------- rh after_mix to edmf_mynn (units: percent) at full level -------
!      if ( id_rh_after_mix > 0) then
!        used = send_data (id_rh_after_mix, am4_Output_edmf%rh_after_mix, Time_next, is, js, 1 )
!      endif
!
!!------- theta after_mix to edmf_mynn (units: K) at full level -------
!      if ( id_th_after_mix > 0) then
!        used = send_data (id_th_after_mix, am4_Output_edmf%th_after_mix, Time_next, is, js, 1 )
!      endif
!
!!------- qa before_pdf to edmf_mynn (units: none) at full level -------
!      if ( id_qa_before_pdf > 0) then
!        used = send_data (id_qa_before_pdf, am4_Output_edmf%qa_before_pdf, Time_next, is, js, 1 )
!      endif
!
!!------- ql before_pdf to edmf_mynn (units: kg/kg) at full level -------
!      if ( id_ql_before_pdf > 0) then
!        used = send_data (id_ql_before_pdf, am4_Output_edmf%ql_before_pdf, Time_next, is, js, 1 )
!      endif
!
!!------- qi before_pdf to edmf_mynn (units: kg/kg) at full level -------
!      if ( id_qi_before_pdf > 0) then
!        used = send_data (id_qi_before_pdf, am4_Output_edmf%qi_before_pdf, Time_next, is, js, 1 )
!      endif
!
!!------- moist updraft area on phalf (units: none) at half level -------
!      if ( id_a_moist_half > 0) then
!        used = send_data (id_a_moist_half, am4_Output_edmf%a_moist_half, Time_next, is, js, 1 )
!      endif
!
!!------- moist updraft area on pfull (units: none) at full level -------
!      if ( id_a_moist_full > 0) then
!        used = send_data (id_a_moist_full, am4_Output_edmf%a_moist_full, Time_next, is, js, 1 )
!      endif
!
!!------- moist updraft mass flux on phalf (units: kg/m2/s) at half level -------
!      if ( id_mf_moist_half > 0) then
!        used = send_data (id_mf_moist_half, am4_Output_edmf%mf_moist_half, Time_next, is, js, 1 )
!      endif
!
!!------- moist updraft mass flux on pfull (units: kg/m2/s) at full level -------
!      if ( id_mf_moist_full > 0) then
!        used = send_data (id_mf_moist_full, am4_Output_edmf%mf_moist_full, Time_next, is, js, 1 )
!      endif
!
!!------- spec humid of moist updraft on phalf (units: kg/kg) at half level -------
!      if ( id_qv_moist_half > 0) then
!        used = send_data (id_qv_moist_half, am4_Output_edmf%qv_moist_half, Time_next, is, js, 1 )
!      endif
!
!!------- spec humid of moist updraft on pfull (units: kg/kg) at full level -------
!      if ( id_qv_moist_full > 0) then
!        used = send_data (id_qv_moist_full, am4_Output_edmf%qv_moist_full, Time_next, is, js, 1 )
!      endif
!
!!------- dry updraft area on phalf (units: none) at half level -------
!      if ( id_a_dry_half > 0) then
!        used = send_data (id_a_dry_half, am4_Output_edmf%a_dry_half, Time_next, is, js, 1 )
!      endif
!
!!------- dry updraft area on pfull (units: none) at full level -------
!      if ( id_a_dry_full > 0) then
!        used = send_data (id_a_dry_full, am4_Output_edmf%a_dry_full, Time_next, is, js, 1 )
!      endif
!
!!------- dry updraft mass flux on phalf (units: kg/m2/s) at half level -------
!      if ( id_mf_dry_half > 0) then
!        used = send_data (id_mf_dry_half, am4_Output_edmf%mf_dry_half, Time_next, is, js, 1 )
!      endif
!
!!------- dry updraft mass flux on pfull (units: kg/m2/s) at full level -------
!      if ( id_mf_dry_full > 0) then
!        used = send_data (id_mf_dry_full, am4_Output_edmf%mf_dry_full, Time_next, is, js, 1 )
!      endif
!
!!------- spec humid of dry updraft on phalf (units: kg/kg) at half level -------
!      if ( id_qv_dry_half > 0) then
!        used = send_data (id_qv_dry_half, am4_Output_edmf%qv_dry_half, Time_next, is, js, 1 )
!      endif
!
!!------- spec humid of dry updraft on pfull (units: kg/kg) at full level -------
!      if ( id_qv_dry_full > 0) then
!        used = send_data (id_qv_dry_full, am4_Output_edmf%qv_dry_full, Time_next, is, js, 1 )
!      endif
!
!!------- updraft mass flux on phalf (units: kg/m2/s) at half level -------
!      if ( id_mf_all_half > 0) then
!        used = send_data (id_mf_all_half, am4_Output_edmf%mf_all_half, Time_next, is, js, 1 )
!      endif
!
!!------- updraft mass flux on pfull (units: kg/m2/s) at full level -------
!      if ( id_mf_all_full > 0) then
!        used = send_data (id_mf_all_full, am4_Output_edmf%mf_all_full, Time_next, is, js, 1 )
!      endif
!
!
!!send_data

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
                   size(Physics_input_block%t,3)) :: &
    rh, qsat

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
  if (do_check_realizability) then
    i=1
    j=1

    !--- qa
    do k=1,kx
      kk=kx-k+1
      tt1 = Physics_input_block%q(i,j,k,nqa)
      tt2 = Physics_input_block%q(i,j,k,nqa) + am4_Output_edmf%qadt_edmf(i,j,k) * dt
      tt3 = Output_edmf%cldfra_bl(i,kk,j)
      if (tt2.ne.tt3) then
        print*,'k,qa_new,qa_bl,qa_old',k,tt2,tt3,tt1
      endif
    enddo

    !--- ql
    write(6,*) ''
    do k=1,kx
      kk=kx-k+1
      tt1 = Physics_input_block%q(i,j,k,nql)
      tt2 = Physics_input_block%q(i,j,k,nql) + am4_Output_edmf%qldt_edmf(i,j,k) * dt
      tt3 = Output_edmf%qc_bl(i,kk,j)
      if (tt2.ne.tt3) then
        print*,'k,ql_new,qc_bl,ql_old',k,tt2,tt3,tt1
      endif
    enddo

    write(6,*) ''
  endif

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
        write(6,*)    '----- end of fieles needed by the offline program ---'
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
        write(6,*)    ''
        write(6,*)    '; k_pbl, ',am4_Output_edmf%kpbl_edmf(ii_write,jj_write)
        write(6,*)    ''
        write(6,*)    ''
        write(6,*)    '; pbl depth,',am4_Output_edmf%pbltop(ii_write,jj_write)
        write(6,*)    ''
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

3000 format (A35,2X,F8.2,',')
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

subroutine modify_mynn_edmf_tendencies (is, ie, js, je, Time_next,      &
                                        do_writeout_column, &
                                        Physics_input_block, Input_edmf, rdt_mynn_ed_am4, &
                                        ix, jx, kx,  &
                                        Output_edmf)

!--- input arguments
  integer, intent(in)                   :: is, ie, js, je
  type(time_type), intent(in)           :: Time_next
  type(physics_input_block_type), intent(in)  :: Physics_input_block
  type(edmf_input_type)     , intent(in)  :: Input_edmf
  integer                   , intent(in)  :: ix, jx, kx
  real, dimension(:,:,:,:)  , intent(in)  :: &
    rdt_mynn_ed_am4
  logical, intent(in) :: do_writeout_column

!--- input/output arguments
  type(edmf_output_type)    , intent(inout)  :: Output_edmf

!--- local variable
  real, dimension(ix,jx,kx) :: tmp_3d
  logical used
  integer i,j,k,kk
!------------------------------------------

!----------------------------
! save the original mynn_edmf tendencies
!----------------------------

!!send_data
!   !------- t tendency from edmf_mynn original (units: K/s) at full level -------
!   if ( id_tdt_edmf_orig > 0) then
!     call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RTHBLTEN(:,:,:)*Input_edmf%exner(:,:,:), tmp_3d)
!     used = send_data (id_tdt_edmf_orig, tmp_3d, Time_next, is, js, 1 )
!   endif
!
!   !------- q tendency from edmf_mynn original (units: kg/kg/s) at full level -------
!   if ( id_qdt_edmf_orig > 0) then
!     call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RQVBLTEN(:,:,:), tmp_3d)
!     used = send_data (id_qdt_edmf_orig, tmp_3d, Time_next, is, js, 1 )
!   endif
!
!   !------- cldfra tendency from edmf_mynn original (units: 1/s) at full level -------
!   if ( id_qadt_edmf_orig > 0) then
!     call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RCCBLTEN(:,:,:), tmp_3d)
!     used = send_data (id_qadt_edmf_orig, tmp_3d, Time_next, is, js, 1 )
!   endif
!
!   !------- qi tendency from edmf_mynn original (units: kg/kg/s) at full level -------
!   if ( id_qidt_edmf_orig > 0) then
!     call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RQIBLTEN(:,:,:), tmp_3d)
!     used = send_data (id_qidt_edmf_orig, tmp_3d, Time_next, is, js, 1 )
!   endif
!
!   !------- ql tendency from edmf_mynn original (units: kg/kg/s) at full level -------
!   if ( id_qldt_edmf_orig > 0) then
!     call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%RQLBLTEN(:,:,:), tmp_3d)
!     used = send_data (id_qldt_edmf_orig, tmp_3d, Time_next, is, js, 1 )
!   endif
!
!!------- ql tendency from edmf_mynn, ED (units: kg/kg/s) at full level -------
!    if ( id_qldt_edmf_ED > 0) then
!      used = send_data (id_qldt_edmf_ED, rdt_mynn_ed_am4(:,:,:,nql), Time_next, is, js, 1 )
!    endif
!
!!------- ql tendency from edmf_mynn, MF (units: kg/kg/s) at full level -------
!    if ( id_qldt_edmf_MF > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql(:,:,:), tmp_3d)
!      used = send_data (id_qldt_edmf_MF, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- qi tendency from edmf_mynn, ED (units: kg/kg/s) at full level -------
!    if ( id_qidt_edmf_ED > 0) then
!      used = send_data (id_qidt_edmf_ED, rdt_mynn_ed_am4(:,:,:,nqi), Time_next, is, js, 1 )
!    endif
!
!!------- qi tendency from edmf_mynn, MF (units: kg/kg/s) at full level -------
!    if ( id_qidt_edmf_MF > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi(:,:,:), tmp_3d)
!      used = send_data (id_qidt_edmf_MF, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- qa tendency from edmf_mynn, ED (units: 1/s) at full level -------
!    if ( id_qadt_edmf_ED > 0) then
!      used = send_data (id_qadt_edmf_ED, rdt_mynn_ed_am4(:,:,:,nqa), Time_next, is, js, 1 )
!    endif
!
!!------- qa tendency from edmf_mynn, MF (units: 1/s) at full level -------
!    if ( id_qadt_edmf_MF > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa(:,:,:), tmp_3d)
!      used = send_data (id_qadt_edmf_MF, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- qa tendency from edmf_mynn, MF_adv (units: 1/s) at full level -------
!    if ( id_qadt_edmf_MF_adv > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_adv(:,:,:), tmp_3d)
!      used = send_data (id_qadt_edmf_MF_adv, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- ql tendency from edmf_mynn, MF_adv (units: kg/kg/s) at full level -------
!    if ( id_qldt_edmf_MF_adv > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_adv(:,:,:), tmp_3d)
!      used = send_data (id_qldt_edmf_MF_adv, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- qi tendency from edmf_mynn, MF_adv (units: kg/kg/s) at full level -------
!    if ( id_qidt_edmf_MF_adv > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_adv(:,:,:), tmp_3d)
!      used = send_data (id_qidt_edmf_MF_adv, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- qa tendency from edmf_mynn, MF_eddy (units: 1/s) at full level -------
!    if ( id_qadt_edmf_MF_eddy > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_eddy(:,:,:), tmp_3d)
!      used = send_data (id_qadt_edmf_MF_eddy, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- ql tendency from edmf_mynn, MF_eddy (units: kg/kg/s) at full level -------
!    if ( id_qldt_edmf_MF_eddy > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_eddy(:,:,:), tmp_3d)
!      used = send_data (id_qldt_edmf_MF_eddy, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- qi tendency from edmf_mynn, MF_eddy (units: kg/kg/s) at full level -------
!    if ( id_qidt_edmf_MF_eddy > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_eddy(:,:,:), tmp_3d)
!      used = send_data (id_qidt_edmf_MF_eddy, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- qa tendency from edmf_mynn, MF_ent (units: 1/s) at full level -------
!    if ( id_qadt_edmf_MF_ent > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_ent(:,:,:), tmp_3d)
!      used = send_data (id_qadt_edmf_MF_ent, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- ql tendency from edmf_mynn, MF_ent (units: kg/kg/s) at full level -------
!    if ( id_qldt_edmf_MF_ent > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_ent(:,:,:), tmp_3d)
!      used = send_data (id_qldt_edmf_MF_ent, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- qi tendency from edmf_mynn, MF_ent (units: kg/kg/s) at full level -------
!    if ( id_qidt_edmf_MF_ent > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_ent(:,:,:), tmp_3d)
!      used = send_data (id_qidt_edmf_MF_ent, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- qa tendency from edmf_mynn, MF_det (units: 1/s) at full level -------
!    if ( id_qadt_edmf_MF_det > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_det(:,:,:), tmp_3d)
!      used = send_data (id_qadt_edmf_MF_det, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- ql tendency from edmf_mynn, MF_det (units: kg/kg/s) at full level -------
!    if ( id_qldt_edmf_MF_det > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_det(:,:,:), tmp_3d)
!      used = send_data (id_qldt_edmf_MF_det, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- qi tendency from edmf_mynn, MF_det (units: kg/kg/s) at full level -------
!    if ( id_qidt_edmf_MF_det > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_det(:,:,:), tmp_3d)
!      used = send_data (id_qidt_edmf_MF_det, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- qa tendency from edmf_mynn, MF_sub (units: 1/s) at full level -------
!    if ( id_qadt_edmf_MF_sub > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qa_sub(:,:,:), tmp_3d)
!      used = send_data (id_qadt_edmf_MF_sub, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- ql tendency from edmf_mynn, MF_sub (units: kg/kg/s) at full level -------
!    if ( id_qldt_edmf_MF_sub > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_ql_sub(:,:,:), tmp_3d)
!      used = send_data (id_qldt_edmf_MF_sub, tmp_3d, Time_next, is, js, 1 )
!    endif
!
!!------- qi tendency from edmf_mynn, MF_sub (units: kg/kg/s) at full level -------
!    if ( id_qidt_edmf_MF_sub > 0) then
!      call reshape_mynn_array_to_am4(ix, jx, kx, Output_edmf%Q_qi_sub(:,:,:), tmp_3d)
!      used = send_data (id_qidt_edmf_MF_sub, tmp_3d, Time_next, is, js, 1 )
!    endif
!!send_data

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

!----------------------------
! printout statements 
!----------------------------
   if (do_writeout_column) then
     write(6,*) '; i,j,',ii_write,jj_write
     !--- qa
     write(6,*) 'rdt_mynn_ed_am4, qa',rdt_mynn_ed_am4(ii_write,jj_write,:,nqa)
     write(6,*) 'Output_edmf%Q_qa',Output_edmf%Q_qa(ii_write,:,jj_write)
     write(6,*) 'Output_edmf%RCCBLTEN',Output_edmf%RCCBLTEN(ii_write,:,jj_write)

     !--- ql
     write(6,*) 'rdt_mynn_ed_am4, ql',rdt_mynn_ed_am4(ii_write,jj_write,:,nql)
     write(6,*) 'Output_edmf%Q_ql',Output_edmf%Q_ql(ii_write,:,jj_write)
     write(6,*) 'Output_edmf%RQLBLTEN',Output_edmf%RQLBLTEN(ii_write,:,jj_write)

     !--- qi
     write(6,*) 'rdt_mynn_ed_am4, qi',rdt_mynn_ed_am4(ii_write,jj_write,:,nqi)
     write(6,*) 'Output_edmf%Q_qi',Output_edmf%Q_qi(ii_write,:,jj_write)
     write(6,*) 'Output_edmf%RQIBLTEN',Output_edmf%RQIBLTEN(ii_write,:,jj_write)
   end if

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

subroutine reshape_mynn_array_to_am4_half (ix, jx, kx, mynn_array_3d, am4_array_3d)

!--- input arguments
  integer                   , intent(in)  :: ix, jx, kx
  real, dimension(ix, kx, jx), intent(in) :: &
    mynn_array_3d

!--- output arguments
  real, dimension(ix, jx, kx+1), intent(out) :: &
    am4_array_3d

!--- local variable
  integer i,j,k,kk
!----------------------------------

  am4_array_3d = 0.

  do i=1,ix
  do j=1,jx
    do k=1,kx      ! k index for full levels
      kk=kx+2-k
      am4_array_3d   (i,j,kk) = mynn_array_3d (i,k,j)
    enddo  ! end loop of k, full levels
  enddo  ! end loop of j
  enddo  ! end loop of 1

end subroutine reshape_mynn_array_to_am4_half


!#############################
! Mellor-Yamada
!#############################


!###################################
!  do not copy these subroutiens

subroutine error_mesg(character1, character2, fatal0)
  character character1
  character character2
  integer fatal0
end subroutine error_mesg

!------------------
subroutine compute_qs_one(t, p, qs)
  real, intent(in) :: t, p    ! t in K; p in Pa
  real, intent(out) :: qs     ! kg/kg
  real esl_hpa, esi_hpa, esl_pa, esi_pa

  esl_hpa = exp (53.67957 - 6743.769  /t - 4.8451 *log(t))
  esi_hpa = exp (23.33086 - 6111.72784/t + 0.15215*log(t))

  esl_pa = esl_hpa*100.  ! change unit to Pa
  esi_pa = esi_hpa*100.  ! change unit to Pa

  qs = esl_pa / p * (18./28.8) ! change unit to kg/kg
end subroutine compute_qs_one

!------------------
subroutine compute_qs(t, p, qs)
  real, intent(in), dimension(:,:,:) :: t, p    ! t in K; p in Pa
  real, intent(out), dimension(:,:,:) :: qs     ! kg/kg
  integer i,j,k,ix, jx, kx

  ix = size(t,1)
  jx = size(t,2)
  kx = size(t,3)  

  do i=1,ix
  do j=1,jx
  do k=1,kx
    call compute_qs_one(t(i,j,k),p(i,j,k),qs(i,j,k))
  enddo
  enddo
  enddo

end subroutine compute_qs

!----------------
subroutine rh_calc(p, t, q, rh, do_simple)
  real, dimension(:,:,:) :: &
    p, t, q, rh
  logical :: do_simple

  rh = -99.99

end subroutine rh_calc

!###################################
!###################################

END MODULE module_bl_mynn

!###################################
!###################################
!###################################

program test111
  use module_bl_mynn

  implicit none
  type(physics_input_block_type) :: Physics_input_block
  !type(physics_tendency_block_type) :: Physics_tendency_block
  type(time_type)           :: Time_next
  integer                   :: is, ie, js, je
  !real                      :: dt
  real,    dimension(ni,nj,nhalf) :: p_half, z_half, z_half_actual, z_half_surf0
  real,    dimension(ni,nj,nfull) :: p_full, z_full, z_full_actual, z_full_surf0, uu, vv, tt, qq, qa, ql, qi, thv, ww, omega, th, thl, qt, rh
  !real,    dimension(nhalf) :: p_half, z_half, z_half_actual
  !real,    dimension(nfull) :: p_full, z_full, z_full_actual, uu, vv, tt, qq, thv, ww, ql, qi
  real,    dimension(ni,nj,nfull) :: diff_t, diff_m, diff_t_edmf, diff_m_edmf
  real,    dimension(ni,nj,nfull) :: udt_mf, vdt_mf, tdt_mf, qdt_mf, thvdt_mf, qtdt_mf, thlidt_mf
  real,    dimension(ni,nj,nfull) :: qadt_edmf, qldt_edmf, qidt_edmf, dqa_edmf,  dql_edmf, dqi_edmf
  real,    dimension(ni,nj,nfull) :: qldt_vdif, qidt_vdif

  !real,    dimension(ni,nj)   :: cov_w_thv, cov_w_qt
  real,    dimension(ni,nj)   :: z_pbl, pbltop
  integer, dimension(ni,nj)   :: kpbl_edmf
  real,    dimension(ni,nj)   :: lat, lon, u_star, b_star, q_star, shflx, lhflx, frac_land, area, t_ref, q_ref, u_flux, v_flux
  real,    dimension(ni,nj)   :: buoy_flux, w1_thv1_surf, w1_qt1_surf

  real,    dimension(ni,nj,nfull,ntracer) :: rdiag
  real,    dimension(ni,nj,nfull) :: udt, vdt, tdt
  real,    dimension(ni,nj,nfull,tr) :: rdt 
  real,    dimension(ni,nj,nfull,tr) :: rdt_mynn_ed_am4 

  real,    dimension(ni,nj,nfull) ::  &
    edmf_mc_full, edmf_moist_area, edmf_dry_area, edmf_moist_humidity, edmf_dry_humidity

  real,    dimension(ni,nj,nhalf) ::  &
    edmf_mc_half
  
  integer mm

  integer option_edmf2ls_mp

!==============================
!==============================
!  Input profiles
!==============================
!==============================

is = 1

! initialize
ql = 0.
qi = 0.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if (input_profile == "SCM_am4p0_BOMEX_01") then

  ! pressure at half level (Pa)
  p_half(1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7427.886  ,  9734.321  , 12561.150  , 15967.110  , 20012.316  , 24748.842  , 30221.828  , 36436.993  , 43290.018  , 50486.037  , 57600.283  , 64271.380  , 70313.783  , 75674.609  , 80362.382  , 84411.590  , 87874.281  , 90812.218  , 93288.617  , 95356.558  , 97068.423  , 98466.016  , 99582.586  ,100451.620  ,101083.625  ,101500.000  /)
  
  ! height at half level above the surface (m)
  z_half(1,1,:) = (/   43740.235  , 35431.493  , 31195.300  , 28119.396  , 25663.132  , 23571.964  , 21671.263  , 19892.089  , 18212.585  , 16627.370  , 15132.718  , 13708.833  , 12298.552  , 10907.458  ,  9540.579  ,  8208.839  ,  6945.259  ,  5778.065  ,  4746.878  ,  3868.656  ,  3134.036  ,  2523.258  ,  2017.634  ,  1600.194  ,  1255.884  ,   972.382  ,   739.298  ,   548.501  ,   393.074  ,   267.776  ,   168.570  ,    91.895  ,    36.421  ,     0.000  /)
  
  ! pressure at full level (Pa)
  p_full(1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6459.524  ,  8529.192  , 11087.742  , 14196.098  , 17913.655  , 22296.793  , 27394.277  , 33232.603  , 39765.134  , 46795.850  , 53965.027  , 60874.922  , 67247.343  , 72961.375  , 77995.018  , 82370.399  , 86131.335  , 89335.198  , 92044.865  , 94318.809  , 96209.952  , 97765.555  , 99023.252  ,100016.474  ,100767.292  ,101291.670  /)
  
  ! height at full level above the surface (m)
  z_full(1,1,:) = (/   38655.403  , 33062.725  , 29524.293  , 26806.207  , 24555.793  , 22570.519  , 20736.856  , 19012.337  , 17384.299  , 15848.323  , 14392.334  , 12977.176  , 11578.398  , 10201.277  ,  8853.966  ,  7558.911  ,  6346.711  ,  5251.147  ,  4299.749  ,  3495.846  ,  2824.908  ,  2267.914  ,  1807.204  ,  1426.886  ,  1113.356  ,   855.318  ,   643.551  ,   470.557  ,   330.276  ,   218.080  ,   130.177  ,    64.129  ,    18.198  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu(1,1,:) = (/      -0.000  ,    -0.000  ,     0.000  ,    -0.000  ,    -0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,    -0.114  ,    -0.695  ,    -2.229  ,    -3.681  ,    -4.920  ,    -5.967  ,    -6.777  ,    -7.356  ,    -7.656  ,    -7.600  ,    -7.237  ,    -6.739  ,    -6.358  ,    -5.989  ,    -5.715  ,    -5.453  ,    -5.093  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv(1,1,:) = (/       0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.148  ,     0.489  ,     0.138  ,     0.050  ,     0.036  ,    -0.069  ,    -0.360  ,    -0.867  ,    -1.539  ,    -2.191  ,    -2.660  ,    -2.887  ,    -2.925  ,    -2.886  ,    -2.817  ,    -2.723  ,    -2.566  /)
  
  ! vertical velocity (Pa/s)
  omega(1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.3089E-01  ,  0.6249E-01  ,  0.5029E-01  ,  0.3965E-01  ,  0.3046E-01  ,  0.2265E-01  ,  0.1612E-01  ,  0.1074E-01  ,  0.6458E-02  ,  0.3198E-02  ,  0.9108E-03  /)
  
  ! temperatur at full levels (K)
  tt(1,1,:) = (/     204.627  ,   201.960  ,   201.399  ,   201.227  ,   201.046  ,   200.818  ,   200.632  ,   200.358  ,   200.142  ,   200.155  ,   202.623  ,   213.222  ,   223.575  ,   233.587  ,   243.117  ,   250.252  ,   258.995  ,   266.880  ,   273.320  ,   278.625  ,   283.077  ,   286.294  ,   288.768  ,   290.970  ,   292.590  ,   293.773  ,   294.831  ,   295.662  ,   296.360  ,   297.411  ,   298.273  ,   298.932  ,   299.409  /)
  
  ! specific humidity at full levels (kg/kg)
  qq(1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.4911E-03  ,  0.9510E-03  ,  0.1091E-02  ,  0.1733E-02  ,  0.3003E-02  ,  0.4261E-02  ,  0.5313E-02  ,  0.6562E-02  ,  0.8099E-02  ,  0.9730E-02  ,  0.1126E-01  ,  0.1271E-01  ,  0.1432E-01  ,  0.1604E-01  ,  0.1643E-01  ,  0.1661E-01  ,  0.1676E-01  ,  0.1697E-01  /)
  
  ! surface air temperature (K)
  t_ref=     299.67
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.1779E-01
  
  ! zonal wind stress
  u_flux =     0.8168E-01  
  
  ! meridional wind stress
  v_flux =     0.4115E-01  
  
  !
  u_star =     0.2800E+00
  b_star =     0.2034E-02
  q_star =     0.1857E-03
  shflx  =     0.9416E+01
  lhflx  =     0.6066E-04

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "SCM_am4p0_DCBL_C1_01") then

  ! pressure at half level (Pa)
  p_half(1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7420.191  ,  9704.786  , 12496.665  , 15854.955  , 19839.696  , 24502.722  , 29888.858  , 36004.018  , 42745.803  , 49824.357  , 56822.053  , 63383.605  , 69326.633  , 74599.199  , 79209.737  , 83192.195  , 86597.781  , 89487.253  , 91922.792  , 93956.593  , 95640.213  , 97014.766  , 98113.066  , 98968.000  , 99590.000  ,100000.000  /)
  
  ! height at half level above the surface (m)
  z_half(1,1,:) = (/   42707.378  , 34586.522  , 30391.441  , 27336.905  , 24895.622  , 22815.334  , 20922.374  , 19148.808  , 17478.379  , 15906.018  , 14424.876  , 13030.563  , 11717.193  , 10437.076  ,  9169.515  ,  7924.721  ,  6725.957  ,  5613.088  ,  4625.656  ,  3780.827  ,  3071.551  ,  2480.319  ,  1989.529  ,  1581.806  ,  1244.040  ,   964.809  ,   734.394  ,   545.295  ,   390.947  ,   266.359  ,   167.708  ,    91.455  ,    36.263  ,     0.000  /)
  
  ! pressure at full level (Pa)
  p_full(1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6456.018  ,  8511.448  , 11041.963  , 14109.261  , 17772.939  , 22089.240  , 27106.662  , 32851.634  , 39278.527  , 46194.726  , 53246.591  , 60043.087  , 66310.738  , 71930.712  , 76881.428  , 81184.687  , 84883.602  , 88034.614  , 90699.572  , 92935.984  , 94795.911  , 96325.855  , 97562.886  , 98539.915  , 99278.675  , 99794.860  /)
  
  ! height at full level above the surface (m)
  z_full(1,1,:) = (/   37737.529  , 32240.742  , 28732.042  , 26031.725  , 23794.044  , 21817.968  , 19990.912  , 18273.953  , 16657.070  , 15134.272  , 13700.090  , 12349.360  , 11054.632  ,  9782.320  ,  8527.820  ,  7308.201  ,  6155.317  ,  5108.561  ,  4195.549  ,  3420.892  ,  2772.323  ,  2232.471  ,  1784.001  ,  1411.794  ,  1103.661  ,   849.086  ,   639.499  ,   467.892  ,   328.505  ,   216.941  ,   129.526  ,    63.830  ,    18.119  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu(1,1,:) = (/      -0.000  ,    -0.000  ,    -0.000  ,    -0.000  ,    -0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv(1,1,:) = (/       0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  /)
  
  ! vertical velocity (Pa/s)
  omega(1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! temperatur at full levels (K)
  tt(1,1,:) = (/     200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   207.036  ,   217.797  ,   228.313  ,   238.450  ,   247.957  ,   256.524  ,   263.942  ,   270.194  ,   275.381  ,   279.414  ,   283.773  ,   287.430  ,   290.457  ,   292.958  ,   295.018  ,   296.708  ,   298.085  ,   299.192  ,   300.070  ,   300.761  ,   301.350  /)
  
  ! specific humidity at full levels (kg/kg)
  qq(1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  /)
  
  
  ! surface air temperature (K)
  t_ref=     301.36
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.2218E-02
  
  ! zonal wind stress
  u_flux =     0.0000E+00  
  
  ! meridional wind stress
  v_flux =     0.0000E+00  
  
  !
  u_star =     0.2800E+00
  b_star =     0.6965E-02
  q_star =     0.0000E+00
  shflx  =     0.6954E+02
  lhflx  =     0.0000E+00


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "SCM_am4p0_DCBL_C1_02_u,vdt_NaN") then
  ! pressure at half level (Pa)
  p_half(1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7420.191  ,  9704.786  , 12496.665  , 15854.955  , 19839.696  , 24502.722  , 29888.858  , 36004.018  , 42745.803  , 49824.357  , 56822.053  , 63383.605  , 69326.633  , 74599.199  , 79209.737  , 83192.195  , 86597.781  , 89487.253  , 91922.792  , 93956.593  , 95640.213  , 97014.766  , 98113.066  , 98968.000  , 99590.000  ,100000.000  /)
  
  ! height at half level above the surface (m)
  z_half(1,1,:) = (/   42693.885  , 34573.030  , 30377.948  , 27323.412  , 24882.129  , 22801.841  , 20908.882  , 19135.316  , 17464.886  , 15892.525  , 14411.384  , 13017.071  , 11703.700  , 10423.583  ,  9156.022  ,  7911.228  ,  6712.465  ,  5599.595  ,  4612.163  ,  3767.335  ,  3058.058  ,  2466.793  ,  1975.490  ,  1568.491  ,  1232.161  ,   954.689  ,   726.114  ,   538.788  ,   386.066  ,   262.908  ,   165.465  ,    90.195  ,    35.747  ,     0.000  /)
  
  ! pressure at full level (Pa)
  p_full(1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6456.018  ,  8511.448  , 11041.963  , 14109.261  , 17772.939  , 22089.240  , 27106.662  , 32851.634  , 39278.527  , 46194.726  , 53246.591  , 60043.087  , 66310.738  , 71930.712  , 76881.428  , 81184.687  , 84883.602  , 88034.614  , 90699.572  , 92935.984  , 94795.911  , 96325.855  , 97562.886  , 98539.915  , 99278.675  , 99794.860  /)
  
  ! height at full level above the surface (m)
  z_full(1,1,:) = (/   37724.037  , 32227.250  , 28718.550  , 26018.232  , 23780.551  , 21804.475  , 19977.420  , 18260.460  , 16643.578  , 15120.780  , 13686.597  , 12335.868  , 11041.139  ,  9768.828  ,  8514.327  ,  7294.708  ,  6141.825  ,  5095.068  ,  4182.057  ,  3407.400  ,  2758.814  ,  2218.686  ,  1770.327  ,  1399.202  ,  1092.666  ,   839.890  ,   632.109  ,   462.201  ,   324.340  ,   214.095  ,   127.775  ,    62.942  ,    17.861  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu(1,1,:) = (/       0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv(1,1,:) = (/       0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  /)
  
  ! vertical velocity (Pa/s)
  omega(1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! temperatur at full levels (K)
  tt(1,1,:) = (/     200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   207.036  ,   217.797  ,   228.313  ,   238.450  ,   247.957  ,   256.524  ,   263.942  ,   270.194  ,   275.396  ,   279.707  ,   283.269  ,   286.208  ,   288.628  ,   290.618  ,   292.251  ,   293.584  ,   294.664  ,   295.527  ,   296.202  ,   296.708  ,   297.060  /)
  
  ! specific humidity at full levels (kg/kg)
  qq(1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! surface air temperature (K)
  t_ref=     297.50
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.2372E-02
  
  ! zonal wind stress
  u_flux =     0.0000E+00  
  
  ! meridional wind stress
  v_flux =     0.0000E+00  

  u_star =     0.2800E+00
  b_star =     0.7065E-02
  q_star =     0.0000E+00
  shflx  =     0.7055E+02
  lhflx  =     0.0000E+00

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "AMIP_i27_j01_IndOcn") then

  ! pressure at half level (Pa)
  p_half (1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7428.381  ,  9736.221  , 12565.298  , 15974.325  , 20023.420  , 24764.674  , 30243.246  , 36464.844  , 43325.024  , 50528.599  , 57650.343  , 64328.485  , 70377.281  , 75743.784  , 80436.526  , 84490.026  , 87956.392  , 90897.445  , 93376.473  , 95446.610  , 97160.292  , 98559.367  , 99677.113  ,100547.053  ,101179.702  ,101596.487  /)
  
  ! height at half level above the surface (m)
  z_half (1,1,:) = (/   47252.234  , 36890.495  , 31857.423  , 28366.625  , 25696.291  , 23612.758  , 21896.202  , 20256.037  , 18623.669  , 17016.207  , 15454.237  , 13939.326  , 12463.676  , 11022.190  ,  9613.907  ,  8250.307  ,  6955.784  ,  5770.165  ,  4734.007  ,  3855.928  ,  3123.765  ,  2516.300  ,  2012.134  ,  1595.198  ,  1251.494  ,   968.583  ,   736.215  ,   546.193  ,   391.380  ,   266.548  ,   167.761  ,    91.439  ,    36.237  ,     0.000  /)
  
  ! pressure at full level (Pa)
  p_full (1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6459.750  ,  8530.333  , 11090.686  , 14201.684  , 17922.706  , 22310.144  , 27412.777  , 33257.109  , 39796.435  , 46834.517  , 54011.240  , 60928.429  , 67307.590  , 73027.672  , 78066.649  , 82446.669  , 86211.595  , 89418.857  , 92131.400  , 94407.759  , 96300.910  , 97858.163  , 99117.189  ,100111.453  ,100863.047  ,101387.951  /)
  
  ! height at full level above the surface (m)
  z_full (1,1,:) = (/   42071.364  , 34373.959  , 30112.024  , 27031.458  , 24654.524  , 22754.480  , 21076.119  , 19439.853  , 17819.938  , 16235.222  , 14696.781  , 13201.501  , 11742.933  , 10318.049  ,  8932.107  ,  7603.045  ,  6362.974  ,  5252.086  ,  4294.967  ,  3489.846  ,  2820.033  ,  2264.217  ,  1803.666  ,  1423.346  ,  1110.039  ,   852.399  ,   641.204  ,   468.787  ,   328.964  ,   217.155  ,   129.600  ,    63.838  ,    18.119  /)
  
  ! actual height at half level (m)
  z_half_actual (1,1,:) = (/   47252.234  , 36890.495  , 31857.423  , 28366.625  , 25696.291  , 23612.758  , 21896.202  , 20256.037  , 18623.669  , 17016.207  , 15454.237  , 13939.326  , 12463.676  , 11022.190  ,  9613.907  ,  8250.307  ,  6955.784  ,  5770.165  ,  4734.007  ,  3855.928  ,  3123.765  ,  2516.300  ,  2012.134  ,  1595.198  ,  1251.494  ,   968.583  ,   736.215  ,   546.193  ,   391.380  ,   266.548  ,   167.761  ,    91.439  ,    36.237  ,     0.000  /)
  
  ! actual height at full level (m)
  z_full_actual (1,1,:) = (/   42071.364  , 34373.959  , 30112.024  , 27031.458  , 24654.524  , 22754.480  , 21076.119  , 19439.853  , 17819.938  , 16235.222  , 14696.781  , 13201.501  , 11742.933  , 10318.049  ,  8932.107  ,  7603.045  ,  6362.974  ,  5252.086  ,  4294.967  ,  3489.846  ,  2820.033  ,  2264.217  ,  1803.666  ,  1423.346  ,  1110.039  ,   852.399  ,   641.204  ,   468.787  ,   328.964  ,   217.155  ,   129.600  ,    63.838  ,    18.119  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu (1,1,:) = (/     -23.416  ,    -0.006  ,     1.985  ,    -4.653  ,   -26.218  ,   -44.980  ,   -20.413  ,   -22.460  ,   -20.748  ,   -13.546  ,    -5.120  ,     3.187  ,    10.548  ,    15.150  ,    17.409  ,    18.019  ,    17.664  ,    17.114  ,    16.444  ,    15.754  ,    15.000  ,    14.146  ,    13.209  ,    12.186  ,    11.100  ,    10.003  ,     8.902  ,     7.946  ,     7.019  ,     6.277  ,     5.750  ,     5.291  ,     4.929  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv (1,1,:) = (/      -7.918  ,    -4.954  ,    -3.367  ,     3.134  ,    21.217  ,    23.349  ,    12.231  ,     1.335  ,    -1.325  ,    -2.582  ,    -2.298  ,    -0.840  ,     0.736  ,     0.281  ,    -2.190  ,    -4.683  ,    -5.294  ,    -3.918  ,    -2.374  ,    -1.683  ,    -1.277  ,    -0.883  ,    -0.682  ,    -0.554  ,    -0.402  ,    -0.313  ,    -0.342  ,    -0.448  ,    -0.674  ,    -0.915  ,    -1.117  ,    -1.342  ,    -1.542  /)
  
  ! vertical velocity (Pa/s)
  omega (1,1,:) = (/    0.3477E-01  ,  0.1060E+00  ,  0.1797E+00  ,  0.2362E+00  ,  0.2499E+00  , -0.3274E+00  , -0.1115E+01  , -0.1778E+01  , -0.2472E+01  , -0.3111E+01  , -0.3626E+01  , -0.3989E+01  , -0.4059E+01  , -0.3897E+01  , -0.3458E+01  , -0.2692E+01  , -0.1816E+01  , -0.9075E+00  , -0.8331E-01  ,  0.5620E+00  ,  0.1033E+01  ,  0.1366E+01  ,  0.1599E+01  ,  0.1758E+01  ,  0.1859E+01  ,  0.1916E+01  ,  0.1939E+01  ,  0.1941E+01  ,  0.1931E+01  ,  0.1917E+01  ,  0.1903E+01  ,  0.1891E+01  ,  0.1881E+01  /)
  
  ! temperatur at full levels (K)
  tt (1,1,:) = (/     255.188  ,   239.951  ,   228.565  ,   218.765  ,   200.312  ,   181.351  ,   184.919  ,   194.615  ,   202.737  ,   208.888  ,   215.204  ,   222.660  ,   231.171  ,   240.073  ,   248.226  ,   255.692  ,   262.458  ,   267.404  ,   272.354  ,   276.815  ,   280.706  ,   284.559  ,   287.469  ,   289.562  ,   291.267  ,   292.335  ,   293.254  ,   294.280  ,   295.217  ,   295.978  ,   296.698  ,   297.256  ,   297.691  /)
  
  ! potential temperatur at full levels (K)
  th (1,1,:) = (/    1473.035  ,  1042.743  ,   834.681  ,   699.500  ,   573.950  ,   471.757  ,   439.961  ,   425.710  ,   409.611  ,   391.545  ,   375.871  ,   363.878  ,   354.876  ,   347.479  ,   339.978  ,   332.696  ,   325.974  ,   318.861  ,   313.772  ,   309.966  ,   307.083  ,   305.420  ,   303.768  ,   302.101  ,   300.724  ,   299.261  ,   298.116  ,   297.466  ,   297.049  ,   296.729  ,   296.603  ,   296.528  ,   296.521  /)
  
  ! ice-liquid water potential temperatur at full levels (K)
  thl (1,1,:) = (/    1473.035  ,  1042.743  ,   834.681  ,   699.500  ,   573.950  ,   308.587  ,   289.640  ,   296.487  ,   302.232  ,   314.212  ,   324.487  ,   318.724  ,   324.485  ,   323.274  ,   320.678  ,   317.918  ,   316.088  ,   315.234  ,   313.281  ,   309.924  ,   307.062  ,   305.409  ,   303.764  ,   302.100  ,   300.724  ,   299.261  ,   298.116  ,   297.466  ,   297.049  ,   296.729  ,   296.603  ,   296.528  ,   296.521  /)
  
  ! specific humidity at full levels (kg/kg)
  qq (1,1,:) = (/    0.1736E-05  ,  0.1760E-05  ,  0.1754E-05  ,  0.1473E-05  ,  0.1953E-05  ,  0.1037E-03  ,  0.3360E-03  ,  0.6352E-03  ,  0.9513E-03  ,  0.1354E-02  ,  0.2010E-02  ,  0.2539E-02  ,  0.2917E-02  ,  0.3472E-02  ,  0.4198E-02  ,  0.4496E-02  ,  0.4505E-02  ,  0.5463E-02  ,  0.7016E-02  ,  0.7963E-02  ,  0.8925E-02  ,  0.1035E-01  ,  0.1180E-01  ,  0.1299E-01  ,  0.1358E-01  ,  0.1409E-01  ,  0.1465E-01  ,  0.1534E-01  ,  0.1608E-01  ,  0.1725E-01  ,  0.1757E-01  ,  0.1782E-01  ,  0.1808E-01  /)
  
  ! cloud liquid water mixing ratio at full levels (kg/kg)
  ql (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.4486E-10  ,  0.1613E-10  ,  0.0000E+00  ,  0.1190E-07  ,  0.8196E-06  ,  0.4413E-02  ,  0.5882E-02  ,  0.2861E-02  ,  0.2010E-02  ,  0.1671E-02  ,  0.8721E-03  ,  0.7515E-04  ,  0.1745E-03  ,  0.1415E-03  ,  0.1484E-04  ,  0.7707E-05  ,  0.4014E-05  ,  0.1743E-05  ,  0.3188E-06  ,  0.8335E-08  ,  0.5524E-08  ,  0.4088E-08  ,  0.2301E-08  ,  0.9814E-08  ,  0.5044E-08  ,  0.3409E-09  ,  0.1802E-14  ,  0.2397E-30  /)
  
  ! cloud ice water mixing ratio at full levels (kg/kg)
  !qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.2224E-01  ,  0.2240E-01  ,  0.2094E-01  ,  0.1884E-01  ,  0.1462E-01  ,  0.6536E-02  ,  0.4606E-02  ,  0.4494E-02  ,  0.4155E-02  ,  0.3521E-02  ,  0.3257E-02  ,  0.2755E-02  ,  0.9241E-03  ,  0.2604E-04  ,  0.4037E-11  ,  0.6534E-28  ,  0.3008E-56  ,  0.7458E-73  ,  0.3929-196  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1204-213  ,  0.3481-215  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.2224E-01  ,  0.2240E-01  ,  0.2094E-01  ,  0.1884E-01  ,  0.1462E-01  ,  0.6536E-02  ,  0.4606E-02  ,  0.4494E-02  ,  0.4155E-02  ,  0.3521E-02  ,  0.3257E-02  ,  0.2755E-02  ,  0.9241E-03  ,  0.2604E-04  ,  0.4037E-11  ,  0.6534E-28  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! total water mixing ratio (qv+ql+qi) at full levels (kg/kg)
  qt (1,1,:) = (/    0.1736E-05  ,  0.1760E-05  ,  0.1754E-05  ,  0.1473E-05  ,  0.1953E-05  ,  0.2234E-01  ,  0.2273E-01  ,  0.2158E-01  ,  0.1979E-01  ,  0.1598E-01  ,  0.1296E-01  ,  0.1303E-01  ,  0.1027E-01  ,  0.9637E-02  ,  0.9390E-02  ,  0.8625E-02  ,  0.7335E-02  ,  0.6562E-02  ,  0.7183E-02  ,  0.7978E-02  ,  0.8932E-02  ,  0.1036E-01  ,  0.1180E-01  ,  0.1299E-01  ,  0.1358E-01  ,  0.1409E-01  ,  0.1465E-01  ,  0.1534E-01  ,  0.1608E-01  ,  0.1725E-01  ,  0.1757E-01  ,  0.1782E-01  ,  0.1808E-01  /)
  
  ! surface air temperature (K)
  t_ref=     297.95
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.1858E-01
  
  ! zonal wind stress
  u_flux =    -0.2945E-01  
  
  ! meridional wind stress
  v_flux =    -0.2334E-01  
  
  ! 
  u_star =     0.2337E+00
  b_star =     0.4722E-02
  q_star =     0.1663E-03
  shflx  =     0.2273E+02
  lhflx  =     0.3312E-04
  
  ! rdiag(1,1,:,nQke), input
  rdiag(1,1,:,nQke) = (/    0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1850E-03  ,  0.2976E-02  ,  0.4362E-02  ,  0.8159E-03  ,  0.6526E-03  ,  0.6317E-03  ,  0.4102E-03  ,  0.1351E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.2492E-03  ,  0.3872E-03  ,  0.7568E-03  ,  0.8774E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1186E-03  ,  0.8824E-02  ,  0.6692E-01  ,  0.1804E+00  ,  0.4197E+00  ,  0.7561E+00  ,  0.8824E+00  ,  0.8172E+00  ,  0.7611E+00  /)
  
  ! rdiag(1,1,:,nSh3D), input
  rdiag(1,1,:,nSh3D) = (/    0.3157E-02  ,  0.7487E-03  ,  0.3744E-04  ,  0.2541E-02  ,  0.1389E-01  ,  0.1850E-01  ,  0.9224E-01  ,  0.2169E+00  ,  0.1313E+00  ,  0.6939E-01  ,  0.3555E-01  ,  0.1016E+00  ,  0.8470E-01  ,  0.1817E-01  ,  0.7454E-03  ,  0.7781E-03  ,  0.7433E-02  ,  0.6084E-01  ,  0.4820E-01  ,  0.8846E+00  ,  0.1466E+01  ,  0.2477E-01  ,  0.2554E-01  ,  0.2559E-01  ,  0.2425E-01  ,  0.1591E-01  ,  0.1672E-01  ,  0.1665E-01  ,  0.2204E-01  ,  0.4805E+00  ,  0.5401E+00  ,  0.6148E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nel_pbl), input
  rdiag(1,1,:,nel_pbl) = (/    0.2225E+00  ,  0.2454E+00  ,  0.2632E+00  ,  0.2817E+00  ,  0.2512E+00  ,  0.2384E+00  ,  0.2407E+00  ,  0.1118E+01  ,  0.1212E+01  ,  0.1630E+01  ,  0.1642E+01  ,  0.1895E+01  ,  0.1323E+01  ,  0.9574E+00  ,  0.9654E+00  ,  0.1009E+01  ,  0.1034E+01  ,  0.1912E+01  ,  0.1803E+01  ,  0.1118E+01  ,  0.1118E+01  ,  0.1171E+01  ,  0.8944E+00  ,  0.7647E+00  ,  0.1517E+01  ,  0.1573E+02  ,  0.5120E+02  ,  0.4966E+02  ,  0.4643E+02  ,  0.4099E+02  ,  0.3166E+02  ,  0.1676E+02  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,ncldfra_bl), input
  rdiag(1,1,:,ncldfra_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.2500E+00  ,  0.1000E+01  ,  0.9883E+00  ,  0.7442E+00  ,  0.7387E+00  ,  0.6922E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nqc_bl), input
  rdiag(1,1,:,nqc_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.7090E-03  ,  0.2156E-03  ,  0.8703E-04  ,  0.2856E-04  ,  0.7529E-04  ,  0.5062E-04  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "AMIP_i24_j02_GreenLand_edge") then

  ! pressure at half level (Pa)
  p_half (1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7367.101  ,  9501.016  , 12051.763  , 15081.164  , 18648.739  , 22804.667  , 27591.598  , 33016.794  , 38991.100  , 45259.229  , 51452.812  , 57258.579  , 62515.990  , 67179.623  , 71257.295  , 74779.224  , 77790.826  , 80345.923  , 82499.557  , 84297.816  , 85786.565  , 87002.158  , 87974.408  , 88732.062  , 89285.034  , 89651.051  /)
  
  ! height at half level above the surface (m)
  z_half_actual (1,1,:) = (/   42854.187  , 34090.199  , 29573.783  , 26375.541  , 23883.747  , 21788.197  , 19883.733  , 18133.894  , 16638.747  , 15142.021  , 13749.349  , 12423.072  , 11127.544  ,  9861.366  ,  8623.532  ,  7414.303  ,  6260.664  ,  5209.289  ,  4276.248  ,  3486.350  ,  2830.945  ,  2279.788  ,  1822.710  ,  1445.225  ,  1134.175  ,   878.428  ,   668.157  ,   495.588  ,   355.025  ,   241.867  ,   152.351  ,    83.197  ,    33.049  ,     0.000  /)
  
  ! pressure at full level (Pa)
  p_full (1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6431.799  ,  8388.872  , 10725.887  , 13509.903  , 16801.873  , 20657.073  , 25122.167  , 30223.085  , 35921.183  , 42047.326  , 48289.840  , 54303.980  , 59848.803  , 64819.847  , 69198.436  , 73004.101  , 76275.116  , 79061.493  , 81417.992  , 83395.455  , 85040.019  , 86392.936  , 87487.383  , 88352.694  , 89008.262  , 89467.918  /)
  
  ! height at full level above the surface (m)
  z_full_actual (1,1,:) = (/   38472.193  , 31831.991  , 27974.662  , 25129.644  , 22835.972  , 20835.965  , 19008.813  , 17386.321  , 15890.384  , 14445.685  , 13086.211  , 11775.308  , 10494.455  ,  9242.449  ,  8018.918  ,  6837.484  ,  5734.977  ,  4742.769  ,  3881.299  ,  3158.648  ,  2555.367  ,  2051.249  ,  1633.967  ,  1289.700  ,  1006.302  ,   773.293  ,   581.873  ,   425.306  ,   298.446  ,   197.109  ,   117.774  ,    58.123  ,    16.525  /)
  
  ! actual height at half level (m)
  z_half (1,1,:) = (/   43323.489  , 34559.501  , 30043.085  , 26844.843  , 24353.049  , 22257.499  , 20353.035  , 18603.197  , 17108.050  , 15611.323  , 14218.652  , 12892.375  , 11596.846  , 10330.668  ,  9092.834  ,  7883.606  ,  6729.966  ,  5678.592  ,  4745.551  ,  3955.652  ,  3300.248  ,  2749.090  ,  2292.012  ,  1914.528  ,  1603.478  ,  1347.730  ,  1137.460  ,   964.891  ,   824.327  ,   711.170  ,   621.653  ,   552.500  ,   502.352  ,   469.302  /)
  
  ! actual height at full level (m)
  z_full (1,1,:) = (/   38941.495  , 32301.293  , 28443.964  , 25598.946  , 23305.274  , 21305.267  , 19478.116  , 17855.623  , 16359.686  , 14914.987  , 13555.513  , 12244.610  , 10963.757  ,  9711.751  ,  8488.220  ,  7306.786  ,  6204.279  ,  5212.071  ,  4350.602  ,  3627.950  ,  3024.669  ,  2520.551  ,  2103.270  ,  1759.003  ,  1475.604  ,  1242.595  ,  1051.175  ,   894.609  ,   767.748  ,   666.412  ,   587.077  ,   527.426  ,   485.827  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu (1,1,:) = (/      58.199  ,    39.909  ,    28.622  ,    21.619  ,    20.128  ,    19.520  ,    15.577  ,    -6.039  ,   -19.920  ,   -13.812  ,    -7.243  ,   -12.567  ,   -13.829  ,   -18.276  ,   -21.589  ,   -21.180  ,   -21.653  ,   -22.721  ,   -21.293  ,   -19.214  ,   -15.823  ,   -11.758  ,    -9.357  ,    -6.196  ,    -1.994  ,     1.367  ,     3.627  ,     5.160  ,     6.305  ,     7.790  ,     9.260  ,    13.626  ,    24.728  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv (1,1,:) = (/      22.217  ,    16.256  ,    10.028  ,     6.631  ,     5.016  ,     5.919  ,     7.868  ,     3.534  ,    -5.813  ,    11.002  ,    23.220  ,    17.090  ,     9.450  ,    -1.169  ,    -6.904  ,   -11.565  ,   -17.407  ,   -22.856  ,   -26.247  ,   -32.022  ,   -35.422  ,   -38.666  ,   -42.043  ,   -44.820  ,   -46.854  ,   -48.798  ,   -50.492  ,   -51.665  ,   -53.430  ,   -54.230  ,   -55.628  ,   -55.723  ,   -48.272  /)
  
  ! vertical velocity (Pa/s)
  omega (1,1,:) = (/    0.1451E-01  ,  0.5539E-01  ,  0.1161E+00  ,  0.1908E+00  ,  0.2937E+00  ,  0.4420E+00  ,  0.6230E+00  ,  0.6303E+00  , -0.9104E+00  , -0.2376E+01  , -0.3775E+01  , -0.5115E+01  , -0.5909E+01  , -0.6357E+01  , -0.6483E+01  , -0.6308E+01  , -0.5807E+01  , -0.5017E+01  , -0.3932E+01  , -0.2796E+01  , -0.1690E+01  , -0.6529E+00  ,  0.2977E+00  ,  0.1126E+01  ,  0.1877E+01  ,  0.2513E+01  ,  0.3049E+01  ,  0.3495E+01  ,  0.3870E+01  ,  0.4169E+01  ,  0.4345E+01  ,  0.4326E+01  ,  0.4187E+01  /)
  
  ! temperatur at full levels (K)
  tt (1,1,:) = (/     215.839  ,   215.319  ,   209.409  ,   204.138  ,   201.467  ,   201.215  ,   197.324  ,   183.635  ,   200.832  ,   199.821  ,   201.721  ,   207.999  ,   214.507  ,   221.349  ,   229.408  ,   236.212  ,   240.263  ,   247.804  ,   251.693  ,   254.244  ,   260.946  ,   264.223  ,   266.563  ,   268.405  ,   269.644  ,   270.874  ,   272.688  ,   273.584  ,   274.036  ,   274.486  ,   274.802  ,   275.054  ,   275.229  /)
  
  ! potential temperatur at full levels (K)
  th (1,1,:) = (/    1245.897  ,   935.703  ,   764.728  ,   652.731  ,   577.261  ,   523.432  ,   469.475  ,   402.190  ,   407.705  ,   378.146  ,   357.385  ,   346.248  ,   336.618  ,   328.466  ,   322.912  ,   316.478  ,   307.744  ,   305.095  ,   299.663  ,   294.408  ,   295.358  ,   293.533  ,   291.638  ,   289.998  ,   288.366  ,   287.260  ,   287.208  ,   286.549  ,   285.731  ,   285.172  ,   284.698  ,   284.358  ,   284.121  /)
  
  ! ice-liquid water potential temperatur at full levels (K)
  thl (1,1,:) = (/    1245.897  ,   935.703  ,   764.728  ,   652.731  ,   577.261  ,   523.432  ,   469.454  ,   348.263  ,   245.210  ,   250.843  ,   253.681  ,   258.102  ,   265.988  ,   268.945  ,   272.659  ,   272.209  ,   270.954  ,   272.551  ,   274.063  ,   273.409  ,   273.885  ,   273.950  ,   274.026  ,   273.771  ,   273.509  ,   273.044  ,   272.459  ,   272.516  ,   272.801  ,   272.437  ,   272.255  ,   272.213  ,   272.314  /)
  
  ! specific humidity at full levels (kg/kg)
  qq (1,1,:) = (/    0.1769E-05  ,  0.1753E-05  ,  0.1747E-05  ,  0.1746E-05  ,  0.1747E-05  ,  0.1787E-05  ,  0.1651E-05  ,  0.2710E-04  ,  0.4532E-03  ,  0.9755E-03  ,  0.1752E-02  ,  0.2482E-02  ,  0.2795E-02  ,  0.3294E-02  ,  0.4196E-02  ,  0.4233E-02  ,  0.3624E-02  ,  0.3758E-02  ,  0.3626E-02  ,  0.3131E-02  ,  0.3764E-02  ,  0.3748E-02  ,  0.3602E-02  ,  0.3448E-02  ,  0.3269E-02  ,  0.3199E-02  ,  0.3300E-02  ,  0.3288E-02  ,  0.3217E-02  ,  0.3162E-02  ,  0.3143E-02  ,  0.3216E-02  ,  0.3469E-02  /)
  
  ! cloud liquid water mixing ratio at full levels (kg/kg)
  ql (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1176E-01  ,  0.1123E-01  ,  0.8319E-02  ,  0.7413E-02  ,  0.6025E-02  ,  0.5581E-02  ,  0.4592E-02  ,  0.4197E-02  ,  0.2671E-02  ,  0.1783E-02  ,  0.2461E-02  ,  0.2245E-02  ,  0.1991E-02  ,  0.1842E-02  ,  0.1659E-02  ,  0.1590E-02  ,  0.1730E-02  ,  0.1705E-02  ,  0.1679E-02  ,  0.1691E-02  ,  0.1699E-02  ,  0.1587E-02  ,  0.1379E-02  /)
  
  ! cloud ice water mixing ratio at full levels (kg/kg)
  qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.3162E-05  ,  0.8728E-02  ,  0.2838E-01  ,  0.2385E-01  ,  0.1038E-01  ,  0.8867E-02  ,  0.8617E-02  ,  0.7679E-02  ,  0.7341E-02  ,  0.6790E-02  ,  0.6131E-02  ,  0.5667E-02  ,  0.5266E-02  ,  0.4856E-02  ,  0.4555E-02  ,  0.4269E-02  ,  0.3950E-02  ,  0.3699E-02  ,  0.3461E-02  ,  0.3350E-02  ,  0.3439E-02  ,  0.3246E-02  ,  0.2915E-02  ,  0.2854E-02  ,  0.2759E-02  ,  0.2764E-02  ,  0.2838E-02  /)
  
  ! total water mixing ratio (qv+ql+qi) at full levels (kg/kg)
  qt (1,1,:) = (/    0.1769E-05  ,  0.1753E-05  ,  0.1747E-05  ,  0.1746E-05  ,  0.1747E-05  ,  0.1787E-05  ,  0.4813E-05  ,  0.8756E-02  ,  0.2883E-01  ,  0.2482E-01  ,  0.2389E-01  ,  0.2258E-01  ,  0.1973E-01  ,  0.1839E-01  ,  0.1756E-01  ,  0.1660E-01  ,  0.1435E-01  ,  0.1362E-01  ,  0.1156E-01  ,  0.9769E-02  ,  0.1078E-01  ,  0.1026E-01  ,  0.9543E-02  ,  0.8989E-02  ,  0.8389E-02  ,  0.8139E-02  ,  0.8469E-02  ,  0.8238E-02  ,  0.7811E-02  ,  0.7706E-02  ,  0.7601E-02  ,  0.7568E-02  ,  0.7686E-02  /)
  
  ! surface air temperature (K)
  t_ref=     275.35
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.5075E-02
  
  ! zonal wind stress
  u_flux =    -0.2682E+00  
  
  ! meridional wind stress
  v_flux =     0.4064E+00  
  
  ! 
  u_star =     0.5343E+01
  b_star =    -0.1020E-01
  q_star =    -0.4899E-04
  shflx  =    -0.8335E+02
  lhflx  =    -0.3377E-05
  
  ! rdiag(1,1,:,nQke), input
  rdiag(1,1,:,nQke) = (/    0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.2407E-03  ,  0.2349E-03  ,  0.1104E-03  ,  0.3740E-01  ,  0.6355E-01  ,  0.1453E+00  ,  0.1655E+00  ,  0.1547E+00  ,  0.1405E+00  ,  0.1748E+00  ,  0.2522E+00  ,  0.1683E+00  ,  0.4309E-01  ,  0.1280E+00  ,  0.1672E+00  ,  0.3794E-01  ,  0.4031E-01  ,  0.4862E-01  ,  0.7614E-01  ,  0.2653E+00  ,  0.6045E+00  ,  0.1192E+01  ,  0.2411E+01  ,  0.4585E+01  ,  0.7122E+01  ,  0.2941E+01  /)
  
  ! rdiag(1,1,:,nSh3D), input
  rdiag(1,1,:,nSh3D) = (/    0.3136E-02  ,  0.2902E-02  ,  0.2621E-02  ,  0.6468E-03  ,  0.2006E-02  ,  0.8455E-02  ,  0.4207E-01  ,  0.1029E+00  ,  0.6893E-01  ,  0.6744E-02  ,  0.1333E+00  ,  0.1368E+00  ,  0.1207E+01  ,  0.1390E+00  ,  0.1389E+00  ,  0.1063E+00  ,  0.1192E+00  ,  0.1343E+00  ,  0.6641E-01  ,  0.2302E-01  ,  0.6790E+00  ,  0.1306E+00  ,  0.1282E+00  ,  0.1386E+00  ,  0.1293E+00  ,  0.1314E+00  ,  0.1315E+00  ,  0.1719E+00  ,  0.1778E+00  ,  0.5953E+00  ,  0.3347E+00  ,  0.1549E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nel_pbl), input
  rdiag(1,1,:,nel_pbl) = (/    0.2433E+00  ,  0.2536E+00  ,  0.2676E+00  ,  0.2889E+00  ,  0.3140E+00  ,  0.3114E+00  ,  0.2595E+00  ,  0.3481E+00  ,  0.5943E+00  ,  0.7262E+00  ,  0.1387E+02  ,  0.1973E+02  ,  0.1370E+02  ,  0.2213E+02  ,  0.2144E+02  ,  0.2308E+02  ,  0.2357E+02  ,  0.3029E+02  ,  0.1617E+02  ,  0.7782E+01  ,  0.1446E+02  ,  0.1107E+02  ,  0.5992E+01  ,  0.1210E+02  ,  0.6364E+01  ,  0.1203E+02  ,  0.2536E+02  ,  0.2321E+02  ,  0.2158E+02  ,  0.1911E+02  ,  0.1521E+02  ,  0.9043E+01  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,ncldfra_bl), input
  rdiag(1,1,:,ncldfra_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nqc_bl), input
  rdiag(1,1,:,nqc_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "SCM_am4p0_BOMEX_02") then

  ! pressure at half level (Pa)
  p_half (1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7427.886  ,  9734.321  , 12561.150  , 15967.110  , 20012.316  , 24748.842  , 30221.828  , 36436.993  , 43290.018  , 50486.037  , 57600.283  , 64271.380  , 70313.783  , 75674.609  , 80362.382  , 84411.590  , 87874.281  , 90812.218  , 93288.617  , 95356.558  , 97068.423  , 98466.016  , 99582.586  ,100451.620  ,101083.625  ,101500.000  /)
  
  ! height at half level above the surface (m)
  z_half (1,1,:) = (/   43468.893  , 35355.218  , 31161.668  , 28107.563  , 25666.089  , 23585.306  , 21691.799  , 19917.736  , 18240.972  , 16656.822  , 15163.233  , 13740.654  , 12330.753  , 10939.593  ,  9572.336  ,  8239.522  ,  6965.192  ,  5790.498  ,  4757.148  ,  3876.377  ,  3139.488  ,  2526.977  ,  2019.490  ,  1601.285  ,  1256.716  ,   973.158  ,   739.848  ,   548.762  ,   392.970  ,   267.692  ,   168.510  ,    91.860  ,    36.407  ,     0.000  /)
  
  ! pressure at full level (Pa)
  p_full (1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6459.524  ,  8529.192  , 11087.742  , 14196.098  , 17913.655  , 22296.793  , 27394.277  , 33232.603  , 39765.134  , 46795.850  , 53965.027  , 60874.922  , 67247.343  , 72961.375  , 77995.018  , 82370.399  , 86131.335  , 89335.198  , 92044.865  , 94318.809  , 96209.952  , 97765.555  , 99023.252  ,100016.474  ,100767.292  ,101291.670  /)
  
  ! height at full level above the surface (m)
  z_full (1,1,:) = (/   38503.439  , 33010.294  , 29502.503  , 26802.281  , 24564.249  , 22587.652  , 20760.076  , 19039.419  , 17413.243  , 15878.330  , 14423.528  , 13009.194  , 11610.564  , 10233.216  ,  8885.169  ,  7584.065  ,  6362.798  ,  5262.474  ,  4308.721  ,  3502.415  ,  2829.482  ,  2270.692  ,  1808.675  ,  1427.846  ,  1114.160  ,   855.980  ,   643.956  ,   470.635  ,   330.182  ,   218.008  ,   130.130  ,    64.105  ,    18.191  /)
  
  ! actual height at half level (m)
  z_half_actual (1,1,:) = (/   43468.893  , 35355.218  , 31161.668  , 28107.563  , 25666.089  , 23585.306  , 21691.799  , 19917.736  , 18240.972  , 16656.822  , 15163.233  , 13740.654  , 12330.753  , 10939.593  ,  9572.336  ,  8239.522  ,  6965.192  ,  5790.498  ,  4757.148  ,  3876.377  ,  3139.488  ,  2526.977  ,  2019.490  ,  1601.285  ,  1256.716  ,   973.158  ,   739.848  ,   548.762  ,   392.970  ,   267.692  ,   168.510  ,    91.860  ,    36.407  ,     0.000  /)
  
  ! actual height at full level (m)
  z_full_actual (1,1,:) = (/   38503.439  , 33010.294  , 29502.503  , 26802.281  , 24564.249  , 22587.652  , 20760.076  , 19039.419  , 17413.243  , 15878.330  , 14423.528  , 13009.194  , 11610.564  , 10233.216  ,  8885.169  ,  7584.065  ,  6362.798  ,  5262.474  ,  4308.721  ,  3502.415  ,  2829.482  ,  2270.692  ,  1808.675  ,  1427.846  ,  1114.160  ,   855.980  ,   643.956  ,   470.635  ,   330.182  ,   218.008  ,   130.130  ,    64.105  ,    18.191  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu (1,1,:) = (/      -0.000  ,    -0.000  ,     0.000  ,    -0.000  ,    -0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,    -0.629  ,    -2.248  ,    -3.713  ,    -4.949  ,    -5.951  ,    -6.748  ,    -7.458  ,    -7.977  ,    -8.117  ,    -8.131  ,    -8.002  ,    -7.771  ,    -7.503  ,    -7.251  ,    -6.972  ,    -6.549  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv (1,1,:) = (/      -0.000  ,    -0.000  ,     0.000  ,    -0.000  ,    -0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.017  ,     0.003  ,     0.004  ,     0.007  ,     0.006  ,     0.001  ,    -0.007  ,    -0.041  ,    -0.096  ,    -0.170  ,    -0.250  ,    -0.315  ,    -0.358  ,    -0.378  ,    -0.386  ,    -0.381  /)
  
  ! vertical velocity (Pa/s)
  omega (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.3069E-01  ,  0.6254E-01  ,  0.5037E-01  ,  0.3967E-01  ,  0.3045E-01  ,  0.2261E-01  ,  0.1612E-01  ,  0.1074E-01  ,  0.6457E-02  ,  0.3198E-02  ,  0.9107E-03  /)
  
  ! temperatur at full levels (K)
  tt (1,1,:) = (/     199.823  ,   199.927  ,   199.972  ,   200.016  ,   200.047  ,   200.058  ,   200.056  ,   200.031  ,   200.008  ,   200.013  ,   202.437  ,   213.165  ,   223.585  ,   233.652  ,   243.313  ,   252.455  ,   260.700  ,   267.520  ,   274.142  ,   279.519  ,   283.975  ,   287.498  ,   289.196  ,   290.917  ,   292.376  ,   293.811  ,   295.073  ,   296.242  ,   296.335  ,   297.385  ,   298.223  ,   298.865  ,   299.339  /)
  
  ! potential temperatur at full levels (K)
  th (1,1,:) = (/    1153.449  ,   868.813  ,   730.265  ,   639.550  ,   573.193  ,   520.421  ,   475.975  ,   437.563  ,   404.112  ,   374.939  ,   353.612  ,   348.411  ,   343.289  ,   338.250  ,   333.320  ,   328.557  ,   323.867  ,   319.077  ,   315.911  ,   313.074  ,   310.739  ,   308.655  ,   305.674  ,   303.595  ,   301.950  ,   300.853  ,   300.045  ,   299.530  ,   298.255  ,   298.220  ,   298.209  ,   298.213  ,   298.244  /)
  
  ! ice-liquid water potential temperatur at full levels (K)
  thl (1,1,:) = (/    1153.449  ,   868.813  ,   730.265  ,   639.550  ,   573.193  ,   520.421  ,   475.975  ,   437.563  ,   404.112  ,   374.939  ,   353.612  ,   348.411  ,   343.289  ,   338.250  ,   333.320  ,   328.557  ,   323.861  ,   319.076  ,   315.911  ,   313.074  ,   310.739  ,   308.655  ,   305.674  ,   303.595  ,   301.950  ,   300.852  ,   300.045  ,   299.530  ,   298.255  ,   298.220  ,   298.209  ,   298.213  ,   298.244  /)
  
  ! specific humidity at full levels (kg/kg)
  qq (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.8716E-05  ,  0.6939E-03  ,  0.6001E-03  ,  0.1561E-02  ,  0.2804E-02  ,  0.3713E-02  ,  0.4451E-02  ,  0.7145E-02  ,  0.9640E-02  ,  0.1127E-01  ,  0.1266E-01  ,  0.1386E-01  ,  0.1497E-01  ,  0.1592E-01  ,  0.1617E-01  ,  0.1634E-01  ,  0.1650E-01  ,  0.1674E-01  /)
  
  ! cloud liquid water mixing ratio at full levels (kg/kg)
  ql (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.3672E-08  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.4867E-07  ,  0.7777E-07  ,  0.6465E-07  ,  0.7452E-07  ,  0.2496E-07  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  , -0.5680E-26  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! cloud ice water mixing ratio at full levels (kg/kg)
  qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.2592E-09  ,  0.1876E-05  ,  0.3932E-06  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! total water mixing ratio (qv+ql+qi) at full levels (kg/kg)
  qt (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.8716E-05  ,  0.6958E-03  ,  0.6005E-03  ,  0.1561E-02  ,  0.2804E-02  ,  0.3713E-02  ,  0.4451E-02  ,  0.7145E-02  ,  0.9640E-02  ,  0.1127E-01  ,  0.1266E-01  ,  0.1386E-01  ,  0.1497E-01  ,  0.1592E-01  ,  0.1617E-01  ,  0.1634E-01  ,  0.1650E-01  ,  0.1674E-01  /)
  
  ! surface air temperature (K)
  t_ref=     299.61
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.1759E-01
  
  ! zonal wind stress
  u_flux =     0.9134E-01  
  
  ! meridional wind stress
  v_flux =     0.5315E-02  
  
  ! 
  u_star =     0.2800E+00
  b_star =     0.2034E-02
  q_star =     0.1857E-03
  shflx  =     0.9419E+01
  lhflx  =     0.6068E-04
  
  ! rdiag(1,1,:,nQke), input
  rdiag(1,1,:,nQke) = (/    0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.2374E-02  ,  0.5934E-01  ,  0.1542E+00  ,  0.2672E+00  ,  0.4408E+00  ,  0.5838E+00  ,  0.6532E+00  ,  0.6743E+00  ,  0.6588E+00  /)
  
  ! rdiag(1,1,:,nSh3D), input
  rdiag(1,1,:,nSh3D) = (/    0.4405E-07  ,  0.5336E-07  ,  0.6184E-07  ,  0.6953E-07  ,  0.7680E-07  ,  0.8411E-07  ,  0.9156E-07  ,  0.9931E-07  ,  0.1074E-06  ,  0.1392E-06  ,  0.5557E-06  ,  0.5581E-06  ,  0.5579E-06  ,  0.5573E-06  ,  0.5570E-06  ,  0.5472E-06  ,  0.1476E-02  ,  0.1746E-01  ,  0.1837E-01  ,  0.1927E-01  ,  0.1699E-01  ,  0.9064E-02  ,  0.1561E-01  ,  0.1321E-01  ,  0.2473E-02  ,  0.7787E-03  ,  0.1925E-02  ,  0.7796E-02  ,  0.8841E+00  ,  0.8834E+00  ,  0.9474E+00  ,  0.1037E+01  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nel_pbl), input
  rdiag(1,1,:,nel_pbl) = (/    0.2442E+00  ,  0.2687E+00  ,  0.2893E+00  ,  0.3067E+00  ,  0.3223E+00  ,  0.3373E+00  ,  0.3519E+00  ,  0.3665E+00  ,  0.3812E+00  ,  0.4339E+00  ,  0.8663E+00  ,  0.8681E+00  ,  0.8680E+00  ,  0.8675E+00  ,  0.8673E+00  ,  0.8596E+00  ,  0.8008E+00  ,  0.9455E+00  ,  0.9155E+00  ,  0.9307E+00  ,  0.8860E+00  ,  0.6532E+00  ,  0.8161E+00  ,  0.8620E+00  ,  0.9783E+01  ,  0.4581E+02  ,  0.5771E+02  ,  0.5487E+02  ,  0.4770E+02  ,  0.4298E+02  ,  0.3158E+02  ,  0.1566E+02  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,ncldfra_bl), input
  rdiag(1,1,:,ncldfra_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nqc_bl), input
  rdiag(1,1,:,nqc_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "SCM_am4p0_DCBL_C1_begin") then

  ! pressure at half level (Pa)
  p_half (1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7420.191  ,  9704.786  , 12496.665  , 15854.955  , 19839.696  , 24502.722  , 29888.858  , 36004.018  , 42745.803  , 49824.357  , 56822.053  , 63383.605  , 69326.633  , 74599.199  , 79209.737  , 83192.195  , 86597.781  , 89487.253  , 91922.792  , 93956.593  , 95640.213  , 97014.766  , 98113.066  , 98968.000  , 99590.000  ,100000.000  /)
  
  ! pressure at full level (Pa)
  p_full (1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6456.018  ,  8511.448  , 11041.963  , 14109.261  , 17772.939  , 22089.240  , 27106.662  , 32851.634  , 39278.527  , 46194.726  , 53246.591  , 60043.087  , 66310.738  , 71930.712  , 76881.428  , 81184.687  , 84883.602  , 88034.614  , 90699.572  , 92935.984  , 94795.911  , 96325.855  , 97562.886  , 98539.915  , 99278.675  , 99794.860  /)
  
  ! actual height at half level (m)
  z_half (1,1,:) = (/   42693.885  , 34573.030  , 30377.948  , 27323.412  , 24882.129  , 22801.841  , 20908.882  , 19135.316  , 17464.886  , 15892.525  , 14411.384  , 13017.071  , 11703.700  , 10423.583  ,  9156.022  ,  7911.228  ,  6712.465  ,  5599.595  ,  4612.163  ,  3767.335  ,  3058.058  ,  2466.793  ,  1975.490  ,  1568.491  ,  1232.161  ,   954.689  ,   726.114  ,   538.788  ,   386.066  ,   262.908  ,   165.465  ,    90.195  ,    35.747  ,     0.000  /)
  
  ! actual height at full level (m)
  z_full (1,1,:) = (/   37724.037  , 32227.250  , 28718.550  , 26018.232  , 23780.551  , 21804.475  , 19977.420  , 18260.460  , 16643.578  , 15120.780  , 13686.597  , 12335.868  , 11041.139  ,  9768.828  ,  8514.327  ,  7294.708  ,  6141.825  ,  5095.068  ,  4182.057  ,  3407.400  ,  2758.814  ,  2218.686  ,  1770.327  ,  1399.202  ,  1092.666  ,   839.890  ,   632.109  ,   462.201  ,   324.340  ,   214.095  ,   127.775  ,    62.942  ,    17.861  /)
  
  ! height at half level above the surface (m)
  z_half_surf0 (1,1,:) = (/   42693.885  , 34573.030  , 30377.948  , 27323.412  , 24882.129  , 22801.841  , 20908.882  , 19135.316  , 17464.886  , 15892.525  , 14411.384  , 13017.071  , 11703.700  , 10423.583  ,  9156.022  ,  7911.228  ,  6712.465  ,  5599.595  ,  4612.163  ,  3767.335  ,  3058.058  ,  2466.793  ,  1975.490  ,  1568.491  ,  1232.161  ,   954.689  ,   726.114  ,   538.788  ,   386.066  ,   262.908  ,   165.465  ,    90.195  ,    35.747  ,     0.000  /)
  
  ! height at full level above the surface (m)
  z_full_surf0 (1,1,:) = (/   37724.037  , 32227.250  , 28718.550  , 26018.232  , 23780.551  , 21804.475  , 19977.420  , 18260.460  , 16643.578  , 15120.780  , 13686.597  , 12335.868  , 11041.139  ,  9768.828  ,  8514.327  ,  7294.708  ,  6141.825  ,  5095.068  ,  4182.057  ,  3407.400  ,  2758.814  ,  2218.686  ,  1770.327  ,  1399.202  ,  1092.666  ,   839.890  ,   632.109  ,   462.201  ,   324.340  ,   214.095  ,   127.775  ,    62.942  ,    17.861  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu (1,1,:) = (/       0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv (1,1,:) = (/       0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  /)
  
  ! vertical velocity (Pa/s)
  omega (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! temperatur at full levels (K)
  tt (1,1,:) = (/     200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   207.036  ,   217.797  ,   228.313  ,   238.450  ,   247.957  ,   256.524  ,   263.942  ,   270.194  ,   275.396  ,   279.707  ,   283.269  ,   286.208  ,   288.628  ,   290.618  ,   292.251  ,   293.584  ,   294.664  ,   295.527  ,   296.202  ,   296.708  ,   297.060  /)
  
  ! potential temperatur at full levels (K)
  th (1,1,:) = (/    1154.470  ,   869.131  ,   730.368  ,   639.500  ,   573.057  ,   520.270  ,   475.841  ,   437.562  ,   404.337  ,   375.358  ,   349.968  ,   327.631  ,   318.730  ,   316.249  ,   313.803  ,   311.425  ,   309.177  ,   307.135  ,   305.355  ,   303.844  ,   302.580  ,   301.526  ,   300.652  ,   299.928  ,   299.331  ,   298.838  ,   298.433  ,   298.101  ,   297.832  ,   297.617  ,   297.449  ,   297.323  ,   297.235  /)
  
  ! ice-liquid water potential temperatur at full levels (K)
  thl (1,1,:) = (/    1154.470  ,   869.131  ,   730.368  ,   639.500  ,   573.057  ,   520.270  ,   475.841  ,   437.562  ,   404.337  ,   375.358  ,   349.968  ,   327.631  ,   318.730  ,   316.249  ,   313.803  ,   311.425  ,   309.177  ,   307.135  ,   305.355  ,   303.844  ,   302.580  ,   301.526  ,   300.652  ,   299.928  ,   299.331  ,   298.838  ,   298.433  ,   298.101  ,   297.832  ,   297.617  ,   297.449  ,   297.323  ,   297.235  /)
  
  ! specific humidity at full levels (kg/kg)
  qq (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! cloud liquid water mixing ratio at full levels (kg/kg)
  ql (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! cloud ice water mixing ratio at full levels (kg/kg)
  qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! total water mixing ratio (qv+ql+qi) at full levels (kg/kg)
  qt (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! surface air temperature (K)
  t_ref=     297.50
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.2372E-02
  
  ! zonal wind stress
  u_flux =     0.0000E+00  
  
  ! meridional wind stress
  v_flux =     0.0000E+00  
  
  ! 
  u_star =     0.2800E+00
  b_star =     0.7065E-02
  q_star =     0.0000E+00
  shflx  =     0.7055E+02
  lhflx  =     0.0000E+00
  
  ! rdiag(1,1,:,nQke), input
  rdiag(1,1,:,nQke) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nSh3D), input
  rdiag(1,1,:,nSh3D) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nel_pbl), input
  rdiag(1,1,:,nel_pbl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,ncldfra_bl), input
  rdiag(1,1,:,ncldfra_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nqc_bl), input
  rdiag(1,1,:,nqc_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "SCM_am4p0_RF01_01") then

  ! pressure at half level (Pa)
  p_half (1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7429.322  ,  9739.835  , 12573.187  , 15988.046  , 20044.538  , 24794.785  , 30283.982  , 36517.815  , 43391.604  , 50609.551  , 57745.553  , 64437.098  , 70498.051  , 75875.352  , 80577.543  , 84639.210  , 88112.561  , 91059.544  , 93543.571  , 95617.885  , 97335.022  , 98736.916  , 99856.897  ,100728.562  ,101362.435  ,101780.000  /)
  
  ! pressure at full level (Pa)
  p_full (1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6460.179  ,  8532.503  , 11096.287  , 14212.307  , 17939.921  , 22335.536  , 27447.964  , 33303.717  , 39855.967  , 46908.059  , 54099.135  , 61030.197  , 67422.176  , 73153.765  , 78202.888  , 82591.732  , 86364.245  , 89577.974  , 92295.987  , 94576.937  , 96473.907  , 98034.299  , 99295.854  ,100292.098  ,101045.167  ,101571.074  /)
  
  ! actual height at half level (m)
  z_half (1,1,:) = (/   44249.972  , 36128.820  , 31933.701  , 28879.152  , 26437.866  , 24357.574  , 22464.613  , 20691.047  , 19013.412  , 17419.087  , 15836.495  , 14276.114  , 12744.186  , 11247.032  ,  9789.651  ,  8382.452  ,  7049.474  ,  5831.230  ,  4765.610  ,  3865.194  ,  3117.048  ,  2498.664  ,  1988.179  ,  1568.166  ,  1222.375  ,   940.873  ,   715.847  ,   531.453  ,   380.969  ,   259.511  ,   163.355  ,    89.045  ,    35.289  ,     0.000  /)
  
  ! actual height at full level (m)
  z_full (1,1,:) = (/   39279.942  , 33783.019  , 30274.296  , 27573.970  , 25336.286  , 23360.207  , 21533.151  , 19812.247  , 18180.316  , 16594.153  , 15025.091  , 13481.308  , 11969.095  , 10494.070  ,  9064.114  ,  7696.815  ,  6424.737  ,  5286.710  ,  4307.177  ,  3485.517  ,  2804.069  ,  2240.864  ,  1776.451  ,  1394.112  ,  1080.852  ,   827.855  ,   623.313  ,   455.988  ,   320.095  ,   211.343  ,   126.146  ,    62.139  ,    17.633  /)
  
  ! height at half level above the surface (m)
  z_half_surf0 (1,1,:) = (/   44249.972  , 36128.820  , 31933.701  , 28879.152  , 26437.866  , 24357.574  , 22464.613  , 20691.047  , 19013.412  , 17419.087  , 15836.495  , 14276.114  , 12744.186  , 11247.032  ,  9789.651  ,  8382.452  ,  7049.474  ,  5831.230  ,  4765.610  ,  3865.194  ,  3117.048  ,  2498.664  ,  1988.179  ,  1568.166  ,  1222.375  ,   940.873  ,   715.847  ,   531.453  ,   380.969  ,   259.511  ,   163.355  ,    89.045  ,    35.289  ,     0.000  /)
  
  ! height at full level above the surface (m)
  z_full_surf0 (1,1,:) = (/   39279.942  , 33783.019  , 30274.296  , 27573.970  , 25336.286  , 23360.207  , 21533.151  , 19812.247  , 18180.316  , 16594.153  , 15025.091  , 13481.308  , 11969.095  , 10494.070  ,  9064.114  ,  7696.815  ,  6424.737  ,  5286.710  ,  4307.177  ,  3485.517  ,  2804.069  ,  2240.864  ,  1776.451  ,  1394.112  ,  1080.852  ,   827.855  ,   623.313  ,   455.988  ,   320.095  ,   211.343  ,   126.146  ,    62.139  ,    17.633  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu (1,1,:) = (/       7.420  ,     7.106  ,     7.031  ,     7.013  ,     7.012  ,     7.007  ,     7.014  ,     7.015  ,     7.004  ,     7.001  ,     6.994  ,     6.995  ,     7.000  ,     7.003  ,     7.003  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.001  ,     6.993  ,     7.068  ,     8.138  ,     8.119  ,     8.081  ,     8.025  ,     7.951  ,     7.853  ,     7.712  ,     7.366  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv (1,1,:) = (/      -5.145  ,    -5.420  ,    -5.442  ,    -5.482  ,    -5.475  ,    -5.493  ,    -5.514  ,    -5.499  ,    -5.490  ,    -5.498  ,    -5.499  ,    -5.493  ,    -5.500  ,    -5.504  ,    -5.504  ,    -5.497  ,    -5.497  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.499  ,    -5.497  ,    -4.267  ,    -4.239  ,    -4.202  ,    -4.157  ,    -4.104  ,    -4.040  ,    -3.956  ,    -3.764  /)
  
  ! vertical velocity (Pa/s)
  omega (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1100E-01  ,  0.4527E-01  ,  0.5256E-01  ,  0.4250E-01  ,  0.3445E-01  ,  0.2644E-01  ,  0.1962E-01  ,  0.1393E-01  ,  0.9283E-02  ,  0.5580E-02  ,  0.2763E-02  ,  0.7869E-03  /)
  
  ! temperatur at full levels (K)
  tt (1,1,:) = (/     200.007  ,   200.002  ,   200.001  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   201.015  ,   211.607  ,   221.721  ,   231.309  ,   240.342  ,   248.804  ,   256.671  ,   263.878  ,   270.304  ,   275.817  ,   280.379  ,   284.068  ,   287.028  ,   289.597  ,   291.328  ,   293.283  ,   291.681  ,   283.947  ,   285.474  ,   287.070  ,   288.382  ,   289.443  ,   290.283  ,   290.924  ,   291.403  /)
  
  ! potential temperatur at full levels (K)
  th (1,1,:) = (/    1154.512  ,   869.138  ,   730.371  ,   639.501  ,   573.058  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   322.863  ,   317.933  ,   313.844  ,   310.672  ,   307.691  ,   305.828  ,   300.999  ,   290.526  ,   290.058  ,   290.029  ,   290.023  ,   290.028  ,   290.041  ,   290.061  ,   290.108  /)
  
  ! ice-liquid water potential temperatur at full levels (K)
  thl (1,1,:) = (/    1154.512  ,   869.138  ,   730.371  ,   639.501  ,   573.058  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   322.863  ,   317.933  ,   313.844  ,   310.672  ,   307.691  ,   305.828  ,   300.999  ,   290.021  ,   289.991  ,   290.001  ,   290.010  ,   290.023  ,   290.039  ,   290.060  ,   290.108  /)
  
  ! specific humidity at full levels (kg/kg)
  qq (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.4155E-03  ,  0.1097E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.2582E-02  ,  0.8763E-02  ,  0.9016E-02  ,  0.9083E-02  ,  0.9135E-02  ,  0.9180E-02  ,  0.9224E-02  ,  0.9277E-02  ,  0.9419E-02  /)

  ! cloud fraction
  qa (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.05  ,  0.9  ,  0.9  ,  0.9  ,  0.45  ,  0.45  ,  0.3  ,  0.3  ,  0.3  /)
  
  ! cloud liquid water mixing ratio at full levels (kg/kg)
  ql (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1144E-07  ,  0.1981E-03  ,  0.2658E-04  ,  0.1141E-04  ,  0.4852E-05  ,  0.2023E-05  ,  0.8174E-06  ,  0.3228E-06  ,  0.1737E-06  /)
  
  ! cloud ice water mixing ratio at full levels (kg/kg)
  qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! total water mixing ratio (qv+ql+qi) at full levels (kg/kg)
  qt (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.4155E-03  ,  0.1097E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.2582E-02  ,  0.8961E-02  ,  0.9043E-02  ,  0.9095E-02  ,  0.9140E-02  ,  0.9182E-02  ,  0.9225E-02  ,  0.9277E-02  ,  0.9419E-02  /)

  
  ! surface air temperature (K)
  t_ref=     291.37
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.9769E-02
  
  ! zonal wind stress
  u_flux =    -0.6720E-01  
  
  ! meridional wind stress
  v_flux =     0.3434E-01  
  
  ! 
  u_star =     0.2500E+00
  b_star =     0.2566E-02
  q_star =     0.1524E-03
  shflx  =     0.1500E+02
  lhflx  =     0.4600E-04
  
  ! rdiag(1,1,:,nQke), input
  rdiag(1,1,:,nQke) = (/    0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.2294E-02  ,  0.2765E-01  ,  0.1298E+00  ,  0.2535E+00  ,  0.3633E+00  ,  0.4417E+00  ,  0.4975E+00  ,  0.5332E+00  ,  0.5455E+00  /)
  
  ! rdiag(1,1,:,nSh3D), input
  rdiag(1,1,:,nSh3D) = (/    0.2036E-05  ,  0.2471E-06  ,  0.1564E-06  ,  0.6949E-07  ,  0.7680E-07  ,  0.1202E-06  ,  0.9169E-07  ,  0.1062E-06  ,  0.3404E-06  ,  0.3406E-06  ,  0.3400E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3397E-06  ,  0.3400E-06  ,  0.3402E-06  ,  0.3402E-06  ,  0.3457E-06  ,  0.3510E-06  ,  0.3694E-06  ,  0.3185E-06  ,  0.1272E-05  ,  0.7059E-04  ,  0.1424E-01  ,  0.7208E+00  ,  0.1086E+01  ,  0.1225E+01  ,  0.1402E+01  ,  0.1091E+01  ,  0.1092E+01  ,  0.8883E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nel_pbl), input
  rdiag(1,1,:,nel_pbl) = (/    0.2440E+00  ,  0.2686E+00  ,  0.2892E+00  ,  0.3066E+00  ,  0.3223E+00  ,  0.3374E+00  ,  0.3522E+00  ,  0.3791E+00  ,  0.6782E+00  ,  0.6784E+00  ,  0.6778E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6776E+00  ,  0.6778E+00  ,  0.6780E+00  ,  0.6781E+00  ,  0.6835E+00  ,  0.6887E+00  ,  0.7065E+00  ,  0.6561E+00  ,  0.7510E+00  ,  0.7272E+00  ,  0.4215E+01  ,  0.1959E+02  ,  0.3053E+02  ,  0.3898E+02  ,  0.4464E+02  ,  0.4827E+02  ,  0.3531E+02  ,  0.1679E+02  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,ncldfra_bl), input
  rdiag(1,1,:,ncldfra_bl) = (/    0.4191E-01  ,  0.4191E-01  ,  0.1735E+00  ,  0.2774E+00  ,  0.3437E+00  ,  0.3871E+00  ,  0.4169E+00  ,  0.4379E+00  ,  0.4450E+00  ,  0.3118E+00  ,  0.7450E-01  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.8919E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nqc_bl), input
  rdiag(1,1,:,nqc_bl) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.5257E-04  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "SCM_am4p0_RF01_02") then

  ! pressure at half level (Pa)
  p_half (1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7429.322  ,  9739.835  , 12573.187  , 15988.046  , 20044.538  , 24794.785  , 30283.982  , 36517.815  , 43391.604  , 50609.551  , 57745.553  , 64437.098  , 70498.051  , 75875.352  , 80577.543  , 84639.210  , 88112.561  , 91059.544  , 93543.571  , 95617.885  , 97335.022  , 98736.916  , 99856.897  ,100728.562  ,101362.435  ,101780.000  /)
  
  ! pressure at full level (Pa)
  p_full (1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6460.179  ,  8532.503  , 11096.287  , 14212.307  , 17939.921  , 22335.536  , 27447.964  , 33303.717  , 39855.967  , 46908.059  , 54099.135  , 61030.197  , 67422.176  , 73153.765  , 78202.888  , 82591.732  , 86364.245  , 89577.974  , 92295.987  , 94576.937  , 96473.907  , 98034.299  , 99295.854  ,100292.098  ,101045.167  ,101571.074  /)
  
  ! actual height at half level (m)
  z_half (1,1,:) = (/   44484.676  , 36363.560  , 32168.446  , 29113.899  , 26672.613  , 24592.321  , 22699.361  , 20925.795  , 19248.160  , 17653.835  , 16071.243  , 14510.861  , 12978.934  , 11481.780  , 10024.399  ,  8617.200  ,  7284.222  ,  6065.978  ,  5000.358  ,  4092.232  ,  3326.179  ,  2684.974  ,  2150.316  ,  1706.064  ,  1335.285  ,  1028.802  ,   775.899  ,   568.341  ,   398.926  ,   262.184  ,   164.964  ,    89.922  ,    35.638  ,     0.000  /)
  
  ! actual height at full level (m)
  z_full (1,1,:) = (/   39514.668  , 34017.762  , 30509.041  , 27808.717  , 25571.033  , 23594.955  , 21767.899  , 20046.994  , 18415.064  , 16828.900  , 15259.839  , 13716.056  , 12203.843  , 10728.818  ,  9298.862  ,  7931.563  ,  6659.485  ,  5521.458  ,  4537.999  ,  3703.468  ,  3001.649  ,  2414.966  ,  1926.370  ,  1519.432  ,  1181.203  ,   901.783  ,   671.740  ,   483.382  ,   330.392  ,   213.483  ,   127.389  ,    62.752  ,    17.807  /)
  
  ! height at half level above the surface (m)
  z_half_surf0 (1,1,:) = (/   44484.676  , 36363.560  , 32168.446  , 29113.899  , 26672.613  , 24592.321  , 22699.361  , 20925.795  , 19248.160  , 17653.835  , 16071.243  , 14510.861  , 12978.934  , 11481.780  , 10024.399  ,  8617.200  ,  7284.222  ,  6065.978  ,  5000.358  ,  4092.232  ,  3326.179  ,  2684.974  ,  2150.316  ,  1706.064  ,  1335.285  ,  1028.802  ,   775.899  ,   568.341  ,   398.926  ,   262.184  ,   164.964  ,    89.922  ,    35.638  ,     0.000  /)
  
  ! height at full level above the surface (m)
  z_full_surf0 (1,1,:) = (/   39514.668  , 34017.762  , 30509.041  , 27808.717  , 25571.033  , 23594.955  , 21767.899  , 20046.994  , 18415.064  , 16828.900  , 15259.839  , 13716.056  , 12203.843  , 10728.818  ,  9298.862  ,  7931.563  ,  6659.485  ,  5521.458  ,  4537.999  ,  3703.468  ,  3001.649  ,  2414.966  ,  1926.370  ,  1519.432  ,  1181.203  ,   901.783  ,   671.740  ,   483.382  ,   330.392  ,   213.483  ,   127.389  ,    62.752  ,    17.807  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu (1,1,:) = (/       7.326  ,     7.084  ,     7.017  ,     7.009  ,     7.005  ,     7.005  ,     7.017  ,     7.015  ,     7.002  ,     7.001  ,     6.994  ,     6.994  ,     7.000  ,     7.004  ,     7.004  ,     7.000  ,     7.000  ,     7.000  ,     7.006  ,     7.006  ,     7.006  ,     7.006  ,     7.006  ,     7.006  ,     7.006  ,     7.006  ,     7.006  ,     7.006  ,     7.005  ,     8.255  ,     8.034  ,     7.604  ,     7.060  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv (1,1,:) = (/      -5.063  ,    -5.400  ,    -5.437  ,    -5.480  ,    -5.473  ,    -5.492  ,    -5.510  ,    -5.495  ,    -5.489  ,    -5.497  ,    -5.501  ,    -5.494  ,    -5.500  ,    -5.503  ,    -5.503  ,    -5.497  ,    -5.497  ,    -5.500  ,    -5.509  ,    -5.509  ,    -5.509  ,    -5.509  ,    -5.509  ,    -5.509  ,    -5.509  ,    -5.510  ,    -5.510  ,    -5.511  ,    -5.511  ,    -0.695  ,    -0.635  ,    -0.559  ,    -0.498  /)
  
  ! vertical velocity (Pa/s)
  omega (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.3248E-01  ,  0.5338E-01  ,  0.4259E-01  ,  0.3322E-01  ,  0.2517E-01  ,  0.1837E-01  ,  0.1270E-01  ,  0.9296E-02  ,  0.5594E-02  ,  0.2771E-02  ,  0.7890E-03  /)
  
  ! temperatur at full levels (K)
  tt (1,1,:) = (/     200.006  ,   200.002  ,   200.001  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   201.015  ,   211.607  ,   221.721  ,   231.309  ,   240.342  ,   248.804  ,   256.671  ,   263.878  ,   270.304  ,   275.817  ,   282.778  ,   290.939  ,   297.818  ,   303.585  ,   308.420  ,   314.764  ,   318.063  ,   320.822  ,   323.098  ,   324.967  ,   326.476  ,   291.946  ,   292.397  ,   293.001  ,   293.478  /)
  
  ! potential temperatur at full levels (K)
  th (1,1,:) = (/    1154.507  ,   869.137  ,   730.370  ,   639.501  ,   573.058  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   325.626  ,   325.624  ,   325.642  ,   325.678  ,   325.743  ,   328.227  ,   328.224  ,   328.256  ,   328.286  ,   328.317  ,   328.333  ,   292.536  ,   292.153  ,   292.132  ,   292.174  /)
  
  ! ice-liquid water potential temperatur at full levels (K)
  thl (1,1,:) = (/    1154.507  ,   869.137  ,   730.370  ,   639.501  ,   573.058  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   325.626  ,   325.624  ,   325.642  ,   325.678  ,   325.743  ,   328.227  ,   328.224  ,   328.256  ,   328.286  ,   328.317  ,   328.333  ,   292.431  ,   292.153  ,   292.131  ,   292.174  /)
  
  ! specific humidity at full levels (kg/kg)
  qq (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.7304E-05  ,  0.7192E-05  ,  0.7057E-05  ,  0.6898E-05  ,  0.6684E-05  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5001E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.1313E-01  ,  0.1345E-01  ,  0.1370E-01  ,  0.1395E-01  /)
  
  ! cloud fraction (none)
  qa (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  , 0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  , 0.0000E+00  , 0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  , 0.0000E+00  ,  0.8917E+00  ,  0.8957E+00  ,  0.8977E+00  ,  0.8984E+00  /)
  
  ! cloud liquid water mixing ratio at full levels (kg/kg)
  ql (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  , 0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  , 0.0000E+00  ,  0.4233E-04  ,  0.9195E-07  ,  0.1722E-07  ,  0.5562E-08  /)
  
  ! cloud ice water mixing ratio at full levels (kg/kg)
  qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  , 0.0000E+00  , 0.0000E+00  , 0.0000E+00  , 0.0000E+00  , 0.0000E+00  , 0.0000E+00  , 0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! total water mixing ratio (qv+ql+qi) at full levels (kg/kg)
  qt (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.7304E-05  ,  0.7192E-05  ,  0.7057E-05  ,  0.6898E-05  ,  0.6684E-05  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5001E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.1318E-01  ,  0.1345E-01  ,  0.1370E-01  ,  0.1395E-01  /)
  
  ! surface air temperature (K)
  t_ref=     293.20
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.1377E-01
  
  ! zonal wind stress
  u_flux =    -0.7454E-01  
  
  ! meridional wind stress
  v_flux =     0.5262E-02  
  
  ! 
  u_star =     0.2500E+00
  b_star =     0.2577E-02
  q_star =     0.1539E-03
  shflx  =     0.1500E+02
  lhflx  =     0.4600E-04
  
  ! rdiag(1,1,:,nQke), input
  rdiag(1,1,:,nQke) = (/    0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1883E-03  ,  0.5164E-02  ,  0.2318E+00  ,  0.1930E+01  ,  0.5731E+01  ,  0.7132E+01  ,  0.6354E+01  ,  0.2530E+01  ,  0.1138E+01  ,  0.6344E+00  ,  0.4754E+00  ,  0.3501E+00  ,  0.3735E+00  ,  0.4450E+00  ,  0.4742E+00  /)
  
  ! rdiag(1,1,:,nSh3D), input
  rdiag(1,1,:,nSh3D) = (/    0.1944E-05  ,  0.2360E-06  ,  0.1494E-06  ,  0.6949E-07  ,  0.7680E-07  ,  0.1148E-06  ,  0.9169E-07  ,  0.1062E-06  ,  0.3404E-06  ,  0.3406E-06  ,  0.3400E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3397E-06  ,  0.3400E-06  ,  0.3402E-06  ,  0.5093E-06  ,  0.6739E+00  ,  0.6856E+00  ,  0.7098E+00  ,  0.7617E+00  ,  0.5159E+00  ,  0.8657E-03  ,  0.7487E+00  ,  0.7988E+00  ,  0.9492E+00  ,  0.8318E+00  ,  0.2180E-01  ,  0.2777E-01  ,  0.7089E+00  ,  0.7399E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nel_pbl), input
  rdiag(1,1,:,nel_pbl) = (/    0.2440E+00  ,  0.2686E+00  ,  0.2892E+00  ,  0.3066E+00  ,  0.3223E+00  ,  0.3374E+00  ,  0.3522E+00  ,  0.3791E+00  ,  0.6782E+00  ,  0.6784E+00  ,  0.6778E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6776E+00  ,  0.6778E+00  ,  0.6780E+00  ,  0.8294E+00  ,  0.2041E+01  ,  0.1029E+02  ,  0.3454E+02  ,  0.6726E+02  ,  0.7938E+02  ,  0.4137E+03  ,  0.4956E+02  ,  0.4242E+02  ,  0.3894E+02  ,  0.3488E+02  ,  0.5402E+01  ,  0.2968E+02  ,  0.2418E+02  ,  0.1403E+02  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,ncldfra_bl), input
  rdiag(1,1,:,ncldfra_bl) = (/    0.4191E-01  ,  0.4191E-01  ,  0.1735E+00  ,  0.2774E+00  ,  0.3437E+00  ,  0.3871E+00  ,  0.4169E+00  ,  0.4379E+00  ,  0.4450E+00  ,  0.3118E+00  ,  0.7450E-01  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1225E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  

  ! rdiag(1,1,:,nqc_bl), input
  rdiag(1,1,:,nqc_bl) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1936E-05  ,  0.0000E+00  ,  0.0000E+00  /)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "SCM_am4p0_RF01_03_cloudy") then
  
  ! pressure at half level (Pa)
  p_half (1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7429.322  ,  9739.835  , 12573.187  , 15988.046  , 20044.538  , 24794.785  , 30283.982  , 36517.815  , 43391.604  , 50609.551  , 57745.553  , 64437.098  , 70498.051  , 75875.352  , 80577.543  , 84639.210  , 88112.561  , 91059.544  , 93543.571  , 95617.885  , 97335.022  , 98736.916  , 99856.897  ,100728.562  ,101362.435  ,101780.000  /)
  
  ! pressure at full level (Pa)
  p_full (1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6460.179  ,  8532.503  , 11096.287  , 14212.307  , 17939.921  , 22335.536  , 27447.964  , 33303.717  , 39855.967  , 46908.059  , 54099.135  , 61030.197  , 67422.176  , 73153.765  , 78202.888  , 82591.732  , 86364.245  , 89577.974  , 92295.987  , 94576.937  , 96473.907  , 98034.299  , 99295.854  ,100292.098  ,101045.167  ,101571.074  /)
  
  ! actual height at half level (m)
  z_half (1,1,:) = (/   44248.847  , 36127.907  , 31932.816  , 28878.276  , 26436.992  , 24356.703  , 22463.743  , 20690.177  , 19012.543  , 17418.219  , 15835.627  , 14275.245  , 12743.318  , 11246.163  ,  9788.783  ,  8381.584  ,  7048.606  ,  5830.362  ,  4764.742  ,  3864.326  ,  3116.180  ,  2497.796  ,  1987.557  ,  1567.513  ,  1221.949  ,   938.762  ,   714.071  ,   530.095  ,   379.977  ,   258.827  ,   162.921  ,    88.807  ,    35.195  ,     0.000  /)
  
  ! actual height at full level (m)
  z_full (1,1,:) = (/   39278.947  , 33782.122  , 30273.415  , 27573.096  , 25335.414  , 23359.337  , 21532.281  , 19811.377  , 18179.448  , 16593.284  , 15024.223  , 13480.440  , 11968.226  , 10493.201  ,  9063.246  ,  7695.947  ,  6423.869  ,  5285.842  ,  4306.309  ,  3484.649  ,  2803.201  ,  2240.120  ,  1775.814  ,  1393.573  ,  1079.579  ,   825.913  ,   621.747  ,   454.813  ,   319.257  ,   210.784  ,   125.810  ,    61.973  ,    17.585  /)
  
  ! height at half level above the surface (m)
  z_half_surf0 (1,1,:) = (/   44248.847  , 36127.907  , 31932.816  , 28878.276  , 26436.992  , 24356.703  , 22463.743  , 20690.177  , 19012.543  , 17418.219  , 15835.627  , 14275.245  , 12743.318  , 11246.163  ,  9788.783  ,  8381.584  ,  7048.606  ,  5830.362  ,  4764.742  ,  3864.326  ,  3116.180  ,  2497.796  ,  1987.557  ,  1567.513  ,  1221.949  ,   938.762  ,   714.071  ,   530.095  ,   379.977  ,   258.827  ,   162.921  ,    88.807  ,    35.195  ,     0.000  /)
  
  ! height at full level above the surface (m)
  z_full_surf0 (1,1,:) = (/   39278.947  , 33782.122  , 30273.415  , 27573.096  , 25335.414  , 23359.337  , 21532.281  , 19811.377  , 18179.448  , 16593.284  , 15024.223  , 13480.440  , 11968.226  , 10493.201  ,  9063.246  ,  7695.947  ,  6423.869  ,  5285.842  ,  4306.309  ,  3484.649  ,  2803.201  ,  2240.120  ,  1775.814  ,  1393.573  ,  1079.579  ,   825.913  ,   621.747  ,   454.813  ,   319.257  ,   210.784  ,   125.810  ,    61.973  ,    17.585  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu (1,1,:) = (/       6.901  ,     6.979  ,     6.980  ,     6.994  ,     6.991  ,     6.998  ,     7.007  ,     7.001  ,     6.996  ,     6.999  ,     6.999  ,     6.997  ,     7.000  ,     7.002  ,     7.002  ,     6.999  ,     6.999  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     6.999  ,     6.730  ,     6.701  ,     6.655  ,     6.592  ,     6.513  ,     6.410  ,     6.268  ,     5.929  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv (1,1,:) = (/      -5.288  ,    -5.448  ,    -5.481  ,    -5.493  ,    -5.492  ,    -5.497  ,    -5.496  ,    -5.494  ,    -5.497  ,    -5.499  ,    -5.502  ,    -5.501  ,    -5.500  ,    -5.499  ,    -5.499  ,    -5.499  ,    -5.499  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.495  ,    -4.922  ,    -4.884  ,    -4.827  ,    -4.756  ,    -4.673  ,    -4.574  ,    -4.450  ,    -4.180  /)
  
  ! vertical velocity (Pa/s)
  omega (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1106E-01  ,  0.4531E-01  ,  0.5257E-01  ,  0.4218E-01  ,  0.3442E-01  ,  0.2643E-01  ,  0.1962E-01  ,  0.1393E-01  ,  0.9282E-02  ,  0.5580E-02  ,  0.2763E-02  ,  0.7869E-03  /)
  
  ! temperatur at full levels (K)
  tt (1,1,:) = (/     200.002  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   201.015  ,   211.607  ,   221.721  ,   231.309  ,   240.342  ,   248.804  ,   256.671  ,   263.878  ,   270.304  ,   275.817  ,   280.379  ,   284.068  ,   287.028  ,   289.457  ,   291.350  ,   293.091  ,   293.561  ,   283.562  ,   284.847  ,   286.386  ,   287.660  ,   288.696  ,   289.519  ,   290.149  ,   290.623  /)
  
  ! potential temperatur at full levels (K)
  th (1,1,:) = (/    1154.482  ,   869.133  ,   730.369  ,   639.501  ,   573.057  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   322.863  ,   317.933  ,   313.844  ,   310.521  ,   307.714  ,   305.628  ,   302.938  ,   290.132  ,   289.421  ,   289.338  ,   289.297  ,   289.280  ,   289.278  ,   289.289  ,   289.331  /)
  
  ! ice-liquid water potential temperatur at full levels (K)
  thl (1,1,:) = (/    1154.482  ,   869.133  ,   730.369  ,   639.501  ,   573.057  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   322.863  ,   317.933  ,   313.844  ,   310.521  ,   307.714  ,   305.628  ,   302.938  ,   289.464  ,   289.316  ,   289.290  ,   289.273  ,   289.267  ,   289.271  ,   289.285  ,   289.329  /)
  
  ! relative humidity (%)
  rh (1,1,:) = (/       0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     3.429  ,     8.078  ,    10.092  ,     9.449  ,     8.858  ,    10.877  ,   100.000  ,    98.052  ,    91.388  ,    86.178  ,    82.196  ,    79.250  ,    77.229  ,    76.606  /)
  
  ! specific humidity at full levels (kg/kg)
  qq (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.4155E-03  ,  0.1097E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1828E-02  ,  0.8540E-02  ,  0.8903E-02  ,  0.9004E-02  ,  0.9081E-02  ,  0.9144E-02  ,  0.9203E-02  ,  0.9268E-02  ,  0.9425E-02  /)
  
  ! cloud fraction (none)
  qa (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1454E-04  ,  0.1000E+01  ,  0.3963E+00  ,  0.2320E+00  ,  0.1425E+00  ,  0.9425E-01  ,  0.6758E-01  ,  0.5277E-01  ,  0.4508E-01  /)
  
  ! cloud liquid water mixing ratio at full levels (kg/kg)
  ql (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1111E-07  ,  0.2624E-03  ,  0.4171E-04  ,  0.1931E-04  ,  0.9515E-05  ,  0.4856E-05  ,  0.2582E-05  ,  0.1449E-05  ,  0.1059E-05  /)
  
  ! cloud ice water mixing ratio at full levels (kg/kg)
  qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! total water mixing ratio (qv+ql+qi) at full levels (kg/kg)
  qt (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.4155E-03  ,  0.1097E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1828E-02  ,  0.8802E-02  ,  0.8944E-02  ,  0.9023E-02  ,  0.9090E-02  ,  0.9149E-02  ,  0.9206E-02  ,  0.9269E-02  ,  0.9426E-02  /)
  
  ! surface air temperature (K)
  t_ref=     290.68
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.9773E-02
  
  ! zonal wind stress
  u_flux =    -0.6184E-01  
  
  ! meridional wind stress
  v_flux =     0.4360E-01  
  
  ! 
  u_star =     0.2500E+00
  b_star =     0.2564E-02
  q_star =     0.1520E-03
  shflx  =     0.1500E+02
  lhflx  =     0.4600E-04
  
  ! rdiag(1,1,:,nQke), input
  rdiag(1,1,:,nQke) = (/    0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1295E-01  ,  0.1306E+00  ,  0.1879E+00  ,  0.2795E+00  ,  0.3817E+00  ,  0.4675E+00  ,  0.5271E+00  ,  0.5646E+00  ,  0.5601E+00  /)
  
  ! rdiag(1,1,:,nSh3D), input
  rdiag(1,1,:,nSh3D) = (/    0.2706E-06  ,  0.5314E-07  ,  0.6152E-07  ,  0.6922E-07  ,  0.7656E-07  ,  0.8394E-07  ,  0.9151E-07  ,  0.1061E-06  ,  0.3395E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3397E-06  ,  0.3400E-06  ,  0.3402E-06  ,  0.3402E-06  ,  0.3457E-06  ,  0.3510E-06  ,  0.3524E-06  ,  0.3371E-06  ,  0.3747E-06  ,  0.3646E-06  ,  0.2609E-02  ,  0.1159E-01  ,  0.8333E+00  ,  0.8241E+00  ,  0.1036E+01  ,  0.1146E+01  ,  0.1071E+01  ,  0.8293E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nel_pbl), input
  rdiag(1,1,:,nel_pbl) = (/    0.2437E+00  ,  0.2682E+00  ,  0.2885E+00  ,  0.3060E+00  ,  0.3218E+00  ,  0.3370E+00  ,  0.3518E+00  ,  0.3788E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6776E+00  ,  0.6778E+00  ,  0.6780E+00  ,  0.6781E+00  ,  0.6835E+00  ,  0.6887E+00  ,  0.6901E+00  ,  0.6750E+00  ,  0.7140E+00  ,  0.7652E+00  ,  0.9613E+01  ,  0.7719E+02  ,  0.4137E+02  ,  0.4770E+02  ,  0.5193E+02  ,  0.5192E+02  ,  0.3652E+02  ,  0.1701E+02  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,ncldfra_bl), input
  rdiag(1,1,:,ncldfra_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.5399E+00  ,  0.3065E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nqc_bl), input
  rdiag(1,1,:,nqc_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.3541E-03  ,  0.1449E-04  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "SCM_RF01_modQDT-Gmy_aTnTtT_a3_5.0h") then

  ! pressure at half level (Pa)
  p_half (1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7429.322  ,  9739.835  , 12573.187  , 15988.046  , 20044.538  , 24794.785  , 30283.982  , 36517.815  , 43391.604  , 50609.551  , 57745.553  , 64437.098  , 70498.051  , 75875.352  , 80577.543  , 84639.210  , 88112.561  , 91059.544  , 93543.571  , 95617.885  , 97335.022  , 98736.916  , 99856.897  ,100728.562  ,101362.435  ,101780.000  /)
  
  ! pressure at full level (Pa)
  p_full (1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6460.179  ,  8532.503  , 11096.287  , 14212.307  , 17939.921  , 22335.536  , 27447.964  , 33303.717  , 39855.967  , 46908.059  , 54099.135  , 61030.197  , 67422.176  , 73153.765  , 78202.888  , 82591.732  , 86364.245  , 89577.974  , 92295.987  , 94576.937  , 96473.907  , 98034.299  , 99295.854  ,100292.098  ,101045.167  ,101571.074  /)
  
  ! actual height at half level (m)
  z_half (1,1,:) = (/   44249.004  , 36128.030  , 31932.933  , 28878.392  , 26437.108  , 24356.818  , 22463.858  , 20690.292  , 19012.658  , 17418.334  , 15835.742  , 14275.360  , 12743.433  , 11246.278  ,  9788.898  ,  8381.698  ,  7048.721  ,  5830.477  ,  4764.857  ,  3864.441  ,  3116.295  ,  2497.911  ,  1987.634  ,  1567.550  ,  1222.072  ,   938.239  ,   713.865  ,   529.895  ,   379.815  ,   258.713  ,   162.848  ,    88.766  ,    35.178  ,     0.000  /)
  
  ! actual height at full level (m)
  z_full (1,1,:) = (/   39279.083  , 33782.241  , 30273.532  , 27573.211  , 25335.529  , 23359.452  , 21532.396  , 19811.492  , 18179.563  , 16593.399  , 15024.338  , 13480.555  , 11968.341  , 10493.316  ,  9063.361  ,  7696.062  ,  6423.984  ,  5285.957  ,  4306.424  ,  3484.764  ,  2803.315  ,  2240.216  ,  1775.871  ,  1393.653  ,  1079.377  ,   825.549  ,   621.544  ,   454.633  ,   319.119  ,   210.690  ,   125.753  ,    61.944  ,    17.577  /)
  
  ! height at half level above the surface (m)
  z_half_surf0 (1,1,:) = (/   44249.004  , 36128.030  , 31932.933  , 28878.392  , 26437.108  , 24356.818  , 22463.858  , 20690.292  , 19012.658  , 17418.334  , 15835.742  , 14275.360  , 12743.433  , 11246.278  ,  9788.898  ,  8381.698  ,  7048.721  ,  5830.477  ,  4764.857  ,  3864.441  ,  3116.295  ,  2497.911  ,  1987.634  ,  1567.550  ,  1222.072  ,   938.239  ,   713.865  ,   529.895  ,   379.815  ,   258.713  ,   162.848  ,    88.766  ,    35.178  ,     0.000  /)
  
  ! height at full level above the surface (m)
  z_full_surf0 (1,1,:) = (/   39279.083  , 33782.241  , 30273.532  , 27573.211  , 25335.529  , 23359.452  , 21532.396  , 19811.492  , 18179.563  , 16593.399  , 15024.338  , 13480.555  , 11968.341  , 10493.316  ,  9063.361  ,  7696.062  ,  6423.984  ,  5285.957  ,  4306.424  ,  3484.764  ,  2803.315  ,  2240.216  ,  1775.871  ,  1393.653  ,  1079.377  ,   825.549  ,   621.544  ,   454.633  ,   319.119  ,   210.690  ,   125.753  ,    61.944  ,    17.577  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu (1,1,:) = (/       6.927  ,     6.987  ,     6.978  ,     6.994  ,     6.990  ,     6.998  ,     7.011  ,     7.003  ,     6.996  ,     6.999  ,     6.998  ,     6.995  ,     7.000  ,     7.003  ,     7.003  ,     6.998  ,     6.998  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.002  ,     6.875  ,     6.802  ,     6.731  ,     6.692  ,     6.652  ,     6.605  ,     6.550  ,     6.477  ,     6.323  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv (1,1,:) = (/      -5.183  ,    -5.423  ,    -5.468  ,    -5.488  ,    -5.487  ,    -5.495  ,    -5.496  ,    -5.492  ,    -5.495  ,    -5.499  ,    -5.503  ,    -5.501  ,    -5.500  ,    -5.499  ,    -5.499  ,    -5.499  ,    -5.499  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.504  ,    -5.143  ,    -4.863  ,    -4.622  ,    -4.528  ,    -4.450  ,    -4.378  ,    -4.306  ,    -4.224  ,    -4.072  /)
  
  ! vertical velocity (Pa/s)
  omega (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1105E-01  ,  0.4530E-01  ,  0.5259E-01  ,  0.4206E-01  ,  0.3445E-01  ,  0.2642E-01  ,  0.1962E-01  ,  0.1393E-01  ,  0.9284E-02  ,  0.5581E-02  ,  0.2764E-02  ,  0.7870E-03  /)
  
  ! temperatur at full levels (K)
  tt (1,1,:) = (/     200.003  ,   200.001  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   201.015  ,   211.607  ,   221.721  ,   231.309  ,   240.342  ,   248.804  ,   256.671  ,   263.878  ,   270.304  ,   275.817  ,   280.379  ,   284.068  ,   287.028  ,   289.478  ,   291.378  ,   293.018  ,   294.287  ,   283.198  ,   284.826  ,   286.253  ,   287.487  ,   288.516  ,   289.339  ,   289.969  ,   290.434  /)
  
  ! potential temperatur at full levels (K)
  th (1,1,:) = (/    1154.487  ,   869.134  ,   730.369  ,   639.501  ,   573.057  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   322.863  ,   317.933  ,   313.844  ,   310.545  ,   307.744  ,   305.551  ,   303.688  ,   289.759  ,   289.400  ,   289.204  ,   289.122  ,   289.099  ,   289.098  ,   289.109  ,   289.144  /)
  
  ! ice-liquid water potential temperatur at full levels (K)
  thl (1,1,:) = (/    1154.487  ,   869.134  ,   730.369  ,   639.501  ,   573.057  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   322.863  ,   317.933  ,   313.844  ,   310.545  ,   307.744  ,   305.551  ,   303.688  ,   289.209  ,   289.249  ,   289.124  ,   289.088  ,   289.083  ,   289.090  ,   289.104  ,   289.141  /)
  
  ! relative humidity (%)
  rh (1,1,:) = (/       0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     3.453  ,     8.139  ,    10.160  ,     9.515  ,     8.982  ,     8.684  ,   100.000  ,    98.877  ,    95.821  ,    90.503  ,    86.226  ,    83.011  ,    80.765  ,    79.729  /)
  
  ! specific humidity at full levels (kg/kg)
  qq (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.4155E-03  ,  0.1097E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1512E-02  ,  0.8333E-02  ,  0.8964E-02  ,  0.9356E-02  ,  0.9424E-02  ,  0.9473E-02  ,  0.9518E-02  ,  0.9568E-02  ,  0.9678E-02  /)
  
  ! cloud fraction (none)
  qa (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1854E-03  ,  0.1000E+01  ,  0.5393E+00  ,  0.3126E+00  ,  0.1518E+00  ,  0.8631E-01  ,  0.5596E-01  ,  0.4049E-01  ,  0.3338E-01  /)
  
  ! cloud liquid water mixing ratio at full levels (kg/kg)
  ql (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.2777E-07  ,  0.2160E-03  ,  0.5953E-04  ,  0.3172E-04  ,  0.1362E-04  ,  0.6559E-05  ,  0.3357E-05  ,  0.1819E-05  ,  0.1206E-05  /)
  
  ! cloud ice water mixing ratio at full levels (kg/kg)
  qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! total water mixing ratio (qv+ql+qi) at full levels (kg/kg)
  qt (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.4155E-03  ,  0.1097E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1512E-02  ,  0.8549E-02  ,  0.9024E-02  ,  0.9388E-02  ,  0.9437E-02  ,  0.9480E-02  ,  0.9521E-02  ,  0.9569E-02  ,  0.9679E-02  /)
  
  ! surface air temperature (K)
  t_ref=     290.51
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.9996E-02
  
  ! zonal wind stress
  u_flux =    -0.6365E-01  
  
  ! meridional wind stress
  v_flux =     0.4099E-01  
  
  ! 
  u_star =     0.2500E+00
  b_star =     0.2563E-02
  q_star =     0.1519E-03
  shflx  =     0.1500E+02
  lhflx  =     0.4600E-04
  
  ! rdiag(1,1,:,nQke), input
  rdiag(1,1,:,nQke) = (/    0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1414E-03  ,  0.3569E-01  ,  0.4693E+00  ,  0.1693E+01  ,  0.1702E+01  ,  0.1385E+01  ,  0.1232E+01  ,  0.1124E+01  ,  0.1031E+01  ,  0.9362E+00  ,  0.8125E+00  /)
  
  ! rdiag(1,1,:,nSh3D), input
  rdiag(1,1,:,nSh3D) = (/    0.5835E-06  ,  0.7072E-07  ,  0.6152E-07  ,  0.6922E-07  ,  0.7656E-07  ,  0.8394E-07  ,  0.9151E-07  ,  0.1061E-06  ,  0.3395E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3397E-06  ,  0.3400E-06  ,  0.3402E-06  ,  0.3402E-06  ,  0.3457E-06  ,  0.3510E-06  ,  0.3550E-06  ,  0.3375E-06  ,  0.3589E-06  ,  0.2028E-02  ,  0.1033E-02  ,  0.4244E+00  ,  0.1141E+00  ,  0.1038E+01  ,  0.1218E+01  ,  0.1257E+01  ,  0.1257E+01  ,  0.9340E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nel_pbl), input
  rdiag(1,1,:,nel_pbl) = (/    0.2437E+00  ,  0.2682E+00  ,  0.2885E+00  ,  0.3060E+00  ,  0.3218E+00  ,  0.3370E+00  ,  0.3518E+00  ,  0.3788E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6776E+00  ,  0.6778E+00  ,  0.6780E+00  ,  0.6781E+00  ,  0.6835E+00  ,  0.6887E+00  ,  0.6926E+00  ,  0.6753E+00  ,  0.2232E+01  ,  0.1435E+02  ,  0.3172E+02  ,  0.8956E+02  ,  0.8530E+02  ,  0.7906E+02  ,  0.7006E+02  ,  0.5702E+02  ,  0.3891E+02  ,  0.1747E+02  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,ncldfra_bl), input
  rdiag(1,1,:,ncldfra_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.9983E+00  ,  0.3789E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nqc_bl), input
  rdiag(1,1,:,nqc_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1933E-03  ,  0.4525E-04  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "SCM_BOMEX_MYNN_ED_mixleng3") then

  ! pressure at half level (Pa)
  p_half (1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7427.886  ,  9734.321  , 12561.150  , 15967.110  , 20012.316  , 24748.842  , 30221.828  , 36436.993  , 43290.018  , 50486.037  , 57600.283  , 64271.380  , 70313.783  , 75674.609  , 80362.382  , 84411.590  , 87874.281  , 90812.218  , 93288.617  , 95356.558  , 97068.423  , 98466.016  , 99582.586  ,100451.620  ,101083.625  ,101500.000  /)
  
  ! pressure at full level (Pa)
  p_full (1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6459.524  ,  8529.192  , 11087.742  , 14196.098  , 17913.655  , 22296.793  , 27394.277  , 33232.603  , 39765.134  , 46795.850  , 53965.027  , 60874.922  , 67247.343  , 72961.375  , 77995.018  , 82370.399  , 86131.335  , 89335.198  , 92044.865  , 94318.809  , 96209.952  , 97765.555  , 99023.252  ,100016.474  ,100767.292  ,101291.670  /)
  
  ! actual height at half level (m)
  z_half (1,1,:) = (/   43483.296  , 35362.441  , 31167.359  , 28112.823  , 25671.540  , 23591.252  , 21698.293  , 19924.727  , 18248.226  , 16664.136  , 15170.644  , 13748.231  , 12338.338  , 10947.118  ,  9579.755  ,  8246.798  ,  6972.288  ,  5797.075  ,  4760.537  ,  3877.825  ,  3139.682  ,  2526.403  ,  2018.250  ,  1600.427  ,  1256.450  ,   973.021  ,   739.821  ,   548.971  ,   393.559  ,   268.120  ,   168.815  ,    92.069  ,    36.553  ,     0.000  /)
  
  ! actual height at full level (m)
  z_full (1,1,:) = (/   38513.448  , 33016.661  , 29507.961  , 26807.643  , 24569.962  , 22593.886  , 20766.831  , 19046.547  , 17420.527  , 15885.694  , 14431.026  , 13016.775  , 11618.118  , 10240.687  ,  8892.515  ,  7591.249  ,  6369.628  ,  5267.422  ,  4311.121  ,  3503.227  ,  2829.288  ,  2269.782  ,  1807.627  ,  1427.286  ,  1113.959  ,   855.898  ,   644.047  ,   471.034  ,   330.690  ,   218.374  ,   130.386  ,    64.282  ,    18.264  /)
  
  ! height at half level above the surface (m)
  z_half_surf0 (1,1,:) = (/   43483.296  , 35362.441  , 31167.359  , 28112.823  , 25671.540  , 23591.252  , 21698.293  , 19924.727  , 18248.226  , 16664.136  , 15170.644  , 13748.231  , 12338.338  , 10947.118  ,  9579.755  ,  8246.798  ,  6972.288  ,  5797.075  ,  4760.537  ,  3877.825  ,  3139.682  ,  2526.403  ,  2018.250  ,  1600.427  ,  1256.450  ,   973.021  ,   739.821  ,   548.971  ,   393.559  ,   268.120  ,   168.815  ,    92.069  ,    36.553  ,     0.000  /)
  
  ! height at full level above the surface (m)
  z_full_surf0 (1,1,:) = (/   38513.448  , 33016.661  , 29507.961  , 26807.643  , 24569.962  , 22593.886  , 20766.831  , 19046.547  , 17420.527  , 15885.694  , 14431.026  , 13016.775  , 11618.118  , 10240.687  ,  8892.515  ,  7591.249  ,  6369.628  ,  5267.422  ,  4311.121  ,  3503.227  ,  2829.288  ,  2269.782  ,  1807.627  ,  1427.286  ,  1113.959  ,   855.898  ,   644.047  ,   471.034  ,   330.690  ,   218.374  ,   130.386  ,    64.282  ,    18.264  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu (1,1,:) = (/      -0.000  ,    -0.000  ,     0.000  ,    -0.000  ,    -0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,    -0.530  ,    -2.251  ,    -3.705  ,    -4.918  ,    -5.923  ,    -6.739  ,    -7.407  ,    -7.973  ,    -8.447  ,    -8.742  ,    -8.752  ,    -8.753  ,    -8.754  ,    -8.755  ,    -8.765  ,    -4.901  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv (1,1,:) = (/      -0.000  ,    -0.000  ,     0.000  ,    -0.000  ,    -0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.001  ,     0.001  ,     0.001  ,     0.001  ,     0.001  ,     0.001  ,     0.000  ,     0.000  ,     0.001  ,    -0.012  ,    -0.054  ,    -0.091  ,    -0.112  ,    -0.140  ,    -0.150  ,    -0.427  /)
  
  ! vertical velocity (Pa/s)
  omega (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.3084E-01  ,  0.6273E-01  ,  0.5040E-01  ,  0.3971E-01  ,  0.3051E-01  ,  0.2270E-01  ,  0.1612E-01  ,  0.1074E-01  ,  0.6461E-02  ,  0.3203E-02  ,  0.9119E-03  /)
  
  ! temperatur at full levels (K)
  tt (1,1,:) = (/     200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   202.413  ,   213.164  ,   223.595  ,   233.670  ,   243.339  ,   252.492  ,   260.925  ,   268.398  ,   274.769  ,   280.064  ,   284.418  ,   287.966  ,   288.890  ,   289.962  ,   292.128  ,   293.509  ,   294.508  ,   295.345  ,   296.660  ,   297.732  ,   298.566  ,   299.195  ,   300.100  /)
  
  ! potential temperatur at full levels (K)
  th (1,1,:) = (/    1154.470  ,   869.131  ,   730.368  ,   639.500  ,   573.057  ,   520.270  ,   475.841  ,   437.495  ,   404.096  ,   374.914  ,   353.571  ,   348.409  ,   343.304  ,   338.276  ,   333.355  ,   328.606  ,   324.147  ,   320.124  ,   316.633  ,   313.684  ,   311.225  ,   309.157  ,   305.350  ,   302.599  ,   301.694  ,   300.544  ,   299.471  ,   298.623  ,   298.582  ,   298.568  ,   298.552  ,   298.543  ,   299.001  /)
  
  ! ice-liquid water potential temperatur at full levels (K)
  thl (1,1,:) = (/    1154.470  ,   869.131  ,   730.368  ,   639.500  ,   573.057  ,   520.270  ,   475.841  ,   437.495  ,   404.096  ,   374.914  ,   353.571  ,   348.409  ,   343.304  ,   338.276  ,   333.355  ,   328.606  ,   324.147  ,   320.124  ,   316.633  ,   313.684  ,   311.225  ,   309.157  ,   305.350  ,   302.599  ,   301.694  ,   300.544  ,   299.471  ,   298.623  ,   298.582  ,   298.568  ,   298.552  ,   298.543  ,   299.001  /)
  
  ! input relative humidity (%)
  rh (1,1,:) = (/       0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  /)
  
  ! specific humidity at full levels (kg/kg)
  qq (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.2799E-03  ,  0.1427E-02  ,  0.2397E-02  ,  0.3206E-02  ,  0.3871E-02  ,  0.7398E-02  ,  0.1229E-01  ,  0.1186E-01  ,  0.1345E-01  ,  0.1488E-01  ,  0.1575E-01  ,  0.1597E-01  ,  0.1604E-01  ,  0.1625E-01  ,  0.1627E-01  ,  0.1890E-01  /)
  
  ! cloud fraction (none)
  qa (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! cloud liquid water mixing ratio at full levels (kg/kg)
  ql (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! cloud ice water mixing ratio at full levels (kg/kg)
  qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! total water mixing ratio (qv+ql+qi) at full levels (kg/kg)
  qt (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.2799E-03  ,  0.1427E-02  ,  0.2397E-02  ,  0.3206E-02  ,  0.3871E-02  ,  0.7398E-02  ,  0.1229E-01  ,  0.1186E-01  ,  0.1345E-01  ,  0.1488E-01  ,  0.1575E-01  ,  0.1597E-01  ,  0.1604E-01  ,  0.1625E-01  ,  0.1627E-01  ,  0.1890E-01  /)
  
  ! surface air temperature (K)
  t_ref=     300.31
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.1952E-01
  
  ! zonal wind stress
  u_flux =     0.0000E+00  
  
  ! meridional wind stress
  v_flux =     0.0000E+00  
  
  ! 
  u_star =     0.2800E+00
  b_star =     0.2030E-02
  q_star =     0.1857E-03
  shflx  =     0.0000E+00
  lhflx  =     0.0000E+00
  
  ! rdiag(1,1,:,nQke), input
  rdiag(1,1,:,nQke) = (/    0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.2308E-03  ,  0.1972E+00  ,  0.1970E+00  ,  0.2697E+00  ,  0.2697E+00  ,  0.5365E+00  /)
  
  ! rdiag(1,1,:,nSh3D), input
  rdiag(1,1,:,nSh3D) = (/    0.4390E-07  ,  0.5314E-07  ,  0.6152E-07  ,  0.6922E-07  ,  0.7656E-07  ,  0.8394E-07  ,  0.9151E-07  ,  0.9930E-07  ,  0.1073E-06  ,  0.1390E-06  ,  0.5581E-06  ,  0.5581E-06  ,  0.5581E-06  ,  0.5581E-06  ,  0.5584E-06  ,  0.5588E-06  ,  0.1233E-02  ,  0.1030E-01  ,  0.1030E-01  ,  0.1029E-01  ,  0.1024E-01  ,  0.7597E-02  ,  0.5962E-02  ,  0.1074E-01  ,  0.1159E-01  ,  0.8246E-02  ,  0.1521E-03  ,  0.0000E+00  ,  0.7128E+00  ,  0.0000E+00  ,  0.7601E+00  ,  0.3113E-04  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nel_pbl), input
  rdiag(1,1,:,nel_pbl) = (/    0.9930E+02  ,  0.9920E+02  ,  0.9912E+02  ,  0.9904E+02  ,  0.9895E+02  ,  0.9886E+02  ,  0.9876E+02  ,  0.9865E+02  ,  0.9852E+02  ,  0.9838E+02  ,  0.9821E+02  ,  0.9801E+02  ,  0.9777E+02  ,  0.9746E+02  ,  0.9706E+02  ,  0.9654E+02  ,  0.9587E+02  ,  0.9501E+02  ,  0.9394E+02  ,  0.9262E+02  ,  0.9099E+02  ,  0.8898E+02  ,  0.8648E+02  ,  0.8340E+02  ,  0.7956E+02  ,  0.7475E+02  ,  0.6872E+02  ,  0.6116E+02  ,  0.5176E+02  ,  0.4031E+02  ,  0.2691E+02  ,  0.1273E+02  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,ncldfra_bl), input
  rdiag(1,1,:,ncldfra_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1339E-01  ,  0.4794E-01  ,  0.9487E-01  ,  0.9678E-01  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.2860E+00  ,  0.3851E+00  /)
  
  ! rdiag(1,1,:,nqc_bl), input
  rdiag(1,1,:,nqc_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1910E-05  ,  0.7611E-05  ,  0.1644E-04  ,  0.1291E-04  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.3568E-03  ,  0.5396E-03  /)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "SCM_RF01_mynn_EDMFexpUP_Gmy_ADD_0.5h") then

  ! pressure at half level (Pa)
  p_half (1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7429.322  ,  9739.835  , 12573.187  , 15988.046  , 20044.538  , 24794.785  , 30283.982  , 36517.815  , 43391.604  , 50609.551  , 57745.553  , 64437.098  , 70498.051  , 75875.352  , 80577.543  , 84639.210  , 88112.561  , 91059.544  , 93543.571  , 95617.885  , 97335.022  , 98736.916  , 99856.897  ,100728.562  ,101362.435  ,101780.000  /)
  
  ! pressure at full level (Pa)
  p_full (1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6460.179  ,  8532.503  , 11096.287  , 14212.307  , 17939.921  , 22335.536  , 27447.964  , 33303.717  , 39855.967  , 46908.059  , 54099.135  , 61030.197  , 67422.176  , 73153.765  , 78202.888  , 82591.732  , 86364.245  , 89577.974  , 92295.987  , 94576.937  , 96473.907  , 98034.299  , 99295.854  ,100292.098  ,101045.167  ,101571.074  /)
  
  ! actual height at half level (m)
  z_half (1,1,:) = (/   44248.731  , 36127.863  , 31932.780  , 28878.244  , 26436.960  , 24356.672  , 22463.713  , 20690.146  , 19012.513  , 17418.189  , 15835.597  , 14275.215  , 12743.288  , 11246.134  ,  9788.753  ,  8381.554  ,  7048.576  ,  5830.333  ,  4764.712  ,  3864.296  ,  3116.150  ,  2497.766  ,  1987.606  ,  1567.601  ,  1221.962  ,   938.115  ,   713.426  ,   529.622  ,   379.664  ,   258.643  ,   162.838  ,    88.805  ,    35.194  ,     0.000  /)
  
  ! actual height at full level (m)
  z_full (1,1,:) = (/   39278.875  , 33782.082  , 30273.381  , 27573.063  , 25335.382  , 23359.306  , 21532.250  , 19811.347  , 18179.418  , 16593.254  , 15024.193  , 13480.410  , 11968.197  , 10493.171  ,  9063.216  ,  7695.917  ,  6423.839  ,  5285.812  ,  4306.279  ,  3484.619  ,  2803.171  ,  2240.130  ,  1775.882  ,  1393.623  ,  1079.260  ,   825.266  ,   621.188  ,   454.421  ,   319.009  ,   210.651  ,   125.768  ,    61.971  ,    17.585  /)
  
  ! height at half level above the surface (m)
  z_half_surf0 (1,1,:) = (/   44248.731  , 36127.863  , 31932.780  , 28878.244  , 26436.960  , 24356.672  , 22463.713  , 20690.146  , 19012.513  , 17418.189  , 15835.597  , 14275.215  , 12743.288  , 11246.134  ,  9788.753  ,  8381.554  ,  7048.576  ,  5830.333  ,  4764.712  ,  3864.296  ,  3116.150  ,  2497.766  ,  1987.606  ,  1567.601  ,  1221.962  ,   938.115  ,   713.426  ,   529.622  ,   379.664  ,   258.643  ,   162.838  ,    88.805  ,    35.194  ,     0.000  /)
  
  ! height at full level above the surface (m)
  z_full_surf0 (1,1,:) = (/   39278.875  , 33782.082  , 30273.381  , 27573.063  , 25335.382  , 23359.306  , 21532.250  , 19811.347  , 18179.418  , 16593.254  , 15024.193  , 13480.410  , 11968.197  , 10493.171  ,  9063.216  ,  7695.917  ,  6423.839  ,  5285.812  ,  4306.279  ,  3484.619  ,  2803.171  ,  2240.130  ,  1775.882  ,  1393.623  ,  1079.260  ,   825.266  ,   621.188  ,   454.421  ,   319.009  ,   210.651  ,   125.768  ,    61.971  ,    17.585  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu (1,1,:) = (/       6.952  ,     6.989  ,     6.993  ,     6.998  ,     6.997  ,     6.999  ,     7.002  ,     7.000  ,     6.999  ,     7.000  ,     7.000  ,     6.999  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.002  ,     6.828  ,     7.213  ,     7.054  ,     6.828  ,     6.998  ,     7.004  ,     6.133  ,     5.993  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv (1,1,:) = (/      -5.452  ,    -5.488  ,    -5.497  ,    -5.499  ,    -5.499  ,    -5.499  ,    -5.498  ,    -5.498  ,    -5.500  ,    -5.500  ,    -5.501  ,    -5.501  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.502  ,    -5.322  ,    -5.722  ,    -5.555  ,    -5.322  ,    -5.499  ,    -5.503  ,    -4.601  ,    -4.458  /)
  
  ! vertical velocity (Pa/s)
  omega (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1106E-01  ,  0.4531E-01  ,  0.5257E-01  ,  0.4206E-01  ,  0.3438E-01  ,  0.2643E-01  ,  0.1962E-01  ,  0.1393E-01  ,  0.9285E-02  ,  0.5583E-02  ,  0.2764E-02  ,  0.7871E-03  /)
  
  ! temperatur at full levels (K)
  tt (1,1,:) = (/     200.001  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   201.015  ,   211.607  ,   221.721  ,   231.309  ,   240.342  ,   248.804  ,   256.671  ,   263.878  ,   270.304  ,   275.817  ,   280.379  ,   284.068  ,   287.028  ,   289.420  ,   291.327  ,   293.139  ,   294.304  ,   283.476  ,   284.606  ,   286.081  ,   287.367  ,   288.417  ,   289.239  ,   290.071  ,   290.546  /)
  
  ! potential temperatur at full levels (K)
  th (1,1,:) = (/    1154.473  ,   869.131  ,   730.368  ,   639.500  ,   573.057  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   322.863  ,   317.933  ,   313.844  ,   310.482  ,   307.690  ,   305.678  ,   303.706  ,   290.044  ,   289.177  ,   289.030  ,   289.002  ,   289.000  ,   288.998  ,   289.210  ,   289.255  /)
  
  ! ice-liquid water potential temperatur at full levels (K)
  thl (1,1,:) = (/    1154.473  ,   869.131  ,   730.368  ,   639.500  ,   573.057  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   322.863  ,   317.933  ,   313.844  ,   310.482  ,   307.690  ,   305.678  ,   303.706  ,   289.880  ,   289.036  ,   289.028  ,   289.002  ,   289.000  ,   288.998  ,   289.210  ,   289.255  /)
  
  ! input relative humidity (%)
  rh (1,1,:) = (/       0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  /)
  
  ! specific humidity at full levels (kg/kg)
  qq (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.4155E-03  ,  0.1097E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.7936E-02  ,  0.8912E-02  ,  0.9007E-02  ,  0.9005E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9665E-02  ,  0.9811E-02  /)
  
  ! cloud fraction (none)
  qa (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.8243E+00  ,  0.9994E+00  ,  0.2382E-01  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! cloud liquid water mixing ratio at full levels (kg/kg)
  ql (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.6456E-04  ,  0.5562E-04  ,  0.8756E-06  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! cloud ice water mixing ratio at full levels (kg/kg)
  qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! total water mixing ratio (qv+ql+qi) at full levels (kg/kg)
  qt (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.4155E-03  ,  0.1097E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.8001E-02  ,  0.8968E-02  ,  0.9008E-02  ,  0.9005E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9665E-02  ,  0.9811E-02  /)
  
  ! surface air temperature (K)
  t_ref=     290.61
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.1011E-01
  
  ! zonal wind stress
  u_flux =    -0.5950E-01  
  
  ! meridional wind stress
  v_flux =     0.4675E-01  
  
  ! 
  u_star =     0.2500E+00
  b_star =     0.2564E-02
  q_star =     0.1520E-03
  shflx  =     0.1500E+02
  lhflx  =     0.4600E-04
  
  ! rdiag(1,1,:,nQke), input
  rdiag(1,1,:,nQke) =     (/ 0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.6897E-02  ,  0.1141E+00  ,  0.1784E+00  ,  0.1277E+00  ,  0.4536E-01  ,  0.1116E-02  ,  0.1718E-01  ,  0.4556E+00  ,  0.5002E+00 /) 
  
  ! rdiag(1,1,:,nSh3D), input
  rdiag(1,1,:,nSh3D) =   (/  0.4390E-07  ,  0.5314E-07  ,  0.6152E-07  ,  0.6922E-07  ,  0.7656E-07  ,  0.8394E-07  ,  0.9151E-07  ,  0.1061E-06  ,  0.3395E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3397E-06  ,  0.3400E-06  ,  0.3402E-06  ,  0.3402E-06  ,  0.3457E-06  ,  0.3510E-06  ,  0.3482E-06  ,  0.3389E-06  ,  0.3907E-06  ,  0.3217E-06  ,  0.6623E-07  ,  0.4382E+00  ,  0.3688E-05  ,  0.1376E-02  ,  0.6698E+00  ,  0.6656E+00  ,  0.6656E+00  ,  0.6647E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nel_pbl), input
  rdiag(1,1,:,nel_pbl) =   (/  0.2437E+00  ,  0.2682E+00  ,  0.2885E+00  ,  0.3060E+00  ,  0.3218E+00  ,  0.3370E+00  ,  0.3518E+00  ,  0.3788E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6776E+00  ,  0.6778E+00  ,  0.6780E+00  ,  0.6781E+00  ,  0.6835E+00  ,  0.6887E+00  ,  0.6859E+00  ,  0.6767E+00  ,  0.7282E+00  ,  0.7945E+00  ,  0.5609E+01  ,  0.2310E+02  ,  0.9460E+02  ,  0.8713E+02  ,  0.2236E+01  ,  0.2236E+01  ,  0.2236E+01  ,  0.1778E+02  ,  0.0000E+00 /) 
  
  ! rdiag(1,1,:,ncldfra_bl), input
  rdiag(1,1,:,ncldfra_bl) =  (/   0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.3532E+00  ,  0.7872E-01  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1538E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nqc_bl), input
  rdiag(1,1,:,nqc_bl) =  (/   0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.7400E-04  ,  0.2089E-05  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.6950E-04  ,  0.0000E+00  ,  0.0000E+00 /) 
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "WRONG.SCM_RF01_mynn_EDMFexpUP_Gmy_ADD_0.5h") then

  ! pressure at half level (Pa)
  p_half (1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7429.322  ,  9739.835  , 12573.187  , 15988.046  , 20044.538  , 24794.785  , 30283.982  , 36517.815  , 43391.604  , 50609.551  , 57745.553  , 64437.098  , 70498.051  , 75875.352  , 80577.543  , 84639.210  , 88112.561  , 91059.544  , 93543.571  , 95617.885  , 97335.022  , 98736.916  , 99856.897  ,100728.562  ,101362.435  ,101780.000  /)
  
  ! pressure at full level (Pa)
  p_full (1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6460.179  ,  8532.503  , 11096.287  , 14212.307  , 17939.921  , 22335.536  , 27447.964  , 33303.717  , 39855.967  , 46908.059  , 54099.135  , 61030.197  , 67422.176  , 73153.765  , 78202.888  , 82591.732  , 86364.245  , 89577.974  , 92295.987  , 94576.937  , 96473.907  , 98034.299  , 99295.854  ,100292.098  ,101045.167  ,101571.074  /)
  
  ! actual height at half level (m)
  z_half (1,1,:) = (/   44248.632  , 36127.777  , 31932.695  , 28878.159  , 26436.876  , 24356.588  , 22463.629  , 20690.063  , 19012.429  , 17418.105  , 15835.513  , 14275.131  , 12743.204  , 11246.050  ,  9788.669  ,  8381.470  ,  7048.493  ,  5830.249  ,  4764.629  ,  3864.212  ,  3116.067  ,  2497.683  ,  1987.535  ,  1567.536  ,  1221.878  ,   938.029  ,   713.317  ,   529.496  ,   379.554  ,   258.534  ,   162.729  ,    88.696  ,    35.146  ,     0.000  /)
  
  ! actual height at full level (m)
  z_full (1,1,:) = (/   39278.784  , 33781.997  , 30273.297  , 27572.979  , 25335.298  , 23359.222  , 21532.167  , 19811.263  , 18179.334  , 16593.171  , 15024.109  , 13480.326  , 11968.113  , 10493.088  ,  9063.132  ,  7695.833  ,  6423.756  ,  5285.729  ,  4306.195  ,  3484.536  ,  2803.087  ,  2240.053  ,  1775.814  ,  1393.549  ,  1079.175  ,   825.169  ,   621.071  ,   454.302  ,   318.900  ,   210.541  ,   125.659  ,    61.893  ,    17.561  /)
  
  ! height at half level above the surface (m)
  z_half_surf0 (1,1,:) = (/   44248.632  , 36127.777  , 31932.695  , 28878.159  , 26436.876  , 24356.588  , 22463.629  , 20690.063  , 19012.429  , 17418.105  , 15835.513  , 14275.131  , 12743.204  , 11246.050  ,  9788.669  ,  8381.470  ,  7048.493  ,  5830.249  ,  4764.629  ,  3864.212  ,  3116.067  ,  2497.683  ,  1987.535  ,  1567.536  ,  1221.878  ,   938.029  ,   713.317  ,   529.496  ,   379.554  ,   258.534  ,   162.729  ,    88.696  ,    35.146  ,     0.000  /)
  
  ! height at full level above the surface (m)
  z_full_surf0 (1,1,:) = (/   39278.784  , 33781.997  , 30273.297  , 27572.979  , 25335.298  , 23359.222  , 21532.167  , 19811.263  , 18179.334  , 16593.171  , 15024.109  , 13480.326  , 11968.113  , 10493.088  ,  9063.132  ,  7695.833  ,  6423.756  ,  5285.729  ,  4306.195  ,  3484.536  ,  2803.087  ,  2240.053  ,  1775.814  ,  1393.549  ,  1079.175  ,   825.169  ,   621.071  ,   454.302  ,   318.900  ,   210.541  ,   125.659  ,    61.893  ,    17.561  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu (1,1,:) = (/       6.975  ,     6.994  ,     6.996  ,     6.999  ,     6.998  ,     7.000  ,     7.001  ,     7.000  ,     6.999  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv (1,1,:) = (/      -5.478  ,    -5.494  ,    -5.499  ,    -5.499  ,    -5.500  ,    -5.500  ,    -5.499  ,    -5.499  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  /)
  
  ! vertical velocity (Pa/s)
  omega (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1107E-01  ,  0.4531E-01  ,  0.5256E-01  ,  0.4205E-01  ,  0.3438E-01  ,  0.2642E-01  ,  0.1962E-01  ,  0.1393E-01  ,  0.9280E-02  ,  0.5578E-02  ,  0.2762E-02  ,  0.7867E-03  /)
  
  ! temperatur at full levels (K)
  tt (1,1,:) = (/     200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   201.015  ,   211.607  ,   221.721  ,   231.309  ,   240.342  ,   248.804  ,   256.671  ,   263.878  ,   270.304  ,   275.817  ,   280.379  ,   284.068  ,   287.028  ,   289.413  ,   291.323  ,   293.154  ,   294.303  ,   283.320  ,   284.734  ,   286.051  ,   287.365  ,   288.417  ,   289.241  ,   289.860  ,   290.290  /)
  
  ! potential temperatur at full levels (K)
  th (1,1,:) = (/    1154.472  ,   869.131  ,   730.368  ,   639.500  ,   573.057  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   322.863  ,   317.933  ,   313.844  ,   310.474  ,   307.685  ,   305.694  ,   303.705  ,   289.884  ,   289.306  ,   289.000  ,   289.000  ,   289.000  ,   289.000  ,   289.000  ,   289.000  /)
  
  ! ice-liquid water potential temperatur at full levels (K)
  thl (1,1,:) = (/    1154.472  ,   869.131  ,   730.368  ,   639.500  ,   573.057  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   322.863  ,   317.933  ,   313.844  ,   310.474  ,   307.685  ,   305.694  ,   303.705  ,   288.769  ,   289.107  ,   289.000  ,   289.000  ,   289.000  ,   289.000  ,   289.000  ,   289.000  /)
  
  ! input relative humidity (%)
  rh (1,1,:) = (/       0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  /)
  
  ! specific humidity at full levels (kg/kg)
  qq (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.4155E-03  ,  0.1097E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.8350E-02  ,  0.8929E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  /)
  
  ! cloud fraction (none)
  qa (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.9726E+00  ,  0.1000E+01  ,  0.2381E-01  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! cloud liquid water mixing ratio at full levels (kg/kg)
  ql (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.4379E-03  ,  0.7878E-04  ,  0.2390E-06  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! cloud ice water mixing ratio at full levels (kg/kg)
  qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! total water mixing ratio (qv+ql+qi) at full levels (kg/kg)
  qt (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.4155E-03  ,  0.1097E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.8788E-02  ,  0.9008E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  /)
  
  ! surface air temperature (K)
  t_ref=     290.39
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.9417E-02
  
  ! zonal wind stress
  u_flux =    -0.5958E-01  
  
  ! meridional wind stress
  v_flux =     0.4681E-01  
  
  ! 
  u_star =     0.2500E+00
  b_star =     0.2562E-02
  q_star =     0.1518E-03
  shflx  =     0.1500E+02
  lhflx  =     0.4600E-04
  
  ! rdiag(1,1,:,nQke), input
  rdiag(1,1,:,nQke) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nSh3D), input
  rdiag(1,1,:,nSh3D) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nel_pbl), input
  rdiag(1,1,:,nel_pbl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,ncldfra_bl), input
  rdiag(1,1,:,ncldfra_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nqc_bl), input
  rdiag(1,1,:,nqc_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "SCM_RF01_rfo76a-M3_EDMFexpUP_NOsm01") then

  ! pressure at half level (Pa)
  p_half (1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7429.322  ,  9739.835  , 12573.187  , 15988.046  , 20044.538  , 24794.785  , 30283.982  , 36517.815  , 43391.604  , 50609.551  , 57745.553  , 64437.098  , 70498.051  , 75875.352  , 80577.543  , 84639.210  , 88112.561  , 91059.544  , 93543.571  , 95617.885  , 97335.022  , 98736.916  , 99856.897  ,100728.562  ,101362.435  ,101780.000  /)
  
  ! pressure at full level (Pa)
  p_full (1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6460.179  ,  8532.503  , 11096.287  , 14212.307  , 17939.921  , 22335.536  , 27447.964  , 33303.717  , 39855.967  , 46908.059  , 54099.135  , 61030.197  , 67422.176  , 73153.765  , 78202.888  , 82591.732  , 86364.245  , 89577.974  , 92295.987  , 94576.937  , 96473.907  , 98034.299  , 99295.854  ,100292.098  ,101045.167  ,101571.074  /)
  
  ! actual height at half level (m)
  z_half (1,1,:) = (/   44248.632  , 36127.777  , 31932.695  , 28878.159  , 26436.876  , 24356.588  , 22463.629  , 20690.063  , 19012.429  , 17418.105  , 15835.513  , 14275.131  , 12743.204  , 11246.050  ,  9788.669  ,  8381.470  ,  7048.493  ,  5830.249  ,  4764.629  ,  3864.212  ,  3116.067  ,  2497.683  ,  1987.535  ,  1567.536  ,  1221.878  ,   938.029  ,   713.317  ,   529.496  ,   379.554  ,   258.534  ,   162.729  ,    88.696  ,    35.146  ,     0.000  /)
  
  ! actual height at full level (m)
  z_full (1,1,:) = (/   39278.784  , 33781.997  , 30273.297  , 27572.979  , 25335.298  , 23359.222  , 21532.167  , 19811.263  , 18179.334  , 16593.171  , 15024.109  , 13480.326  , 11968.113  , 10493.088  ,  9063.132  ,  7695.833  ,  6423.756  ,  5285.729  ,  4306.195  ,  3484.536  ,  2803.087  ,  2240.053  ,  1775.814  ,  1393.549  ,  1079.175  ,   825.169  ,   621.071  ,   454.302  ,   318.900  ,   210.541  ,   125.659  ,    61.893  ,    17.561  /)
  
  ! height at half level above the surface (m)
  z_half_surf0 (1,1,:) = (/   44248.632  , 36127.777  , 31932.695  , 28878.159  , 26436.876  , 24356.588  , 22463.629  , 20690.063  , 19012.429  , 17418.105  , 15835.513  , 14275.131  , 12743.204  , 11246.050  ,  9788.669  ,  8381.470  ,  7048.493  ,  5830.249  ,  4764.629  ,  3864.212  ,  3116.067  ,  2497.683  ,  1987.535  ,  1567.536  ,  1221.878  ,   938.029  ,   713.317  ,   529.496  ,   379.554  ,   258.534  ,   162.729  ,    88.696  ,    35.146  ,     0.000  /)
  
  ! height at full level above the surface (m)
  z_full_surf0 (1,1,:) = (/   39278.784  , 33781.997  , 30273.297  , 27572.979  , 25335.298  , 23359.222  , 21532.167  , 19811.263  , 18179.334  , 16593.171  , 15024.109  , 13480.326  , 11968.113  , 10493.088  ,  9063.132  ,  7695.833  ,  6423.756  ,  5285.729  ,  4306.195  ,  3484.536  ,  2803.087  ,  2240.053  ,  1775.814  ,  1393.549  ,  1079.175  ,   825.169  ,   621.071  ,   454.302  ,   318.900  ,   210.541  ,   125.659  ,    61.893  ,    17.561  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu (1,1,:) = (/       6.975  ,     6.994  ,     6.996  ,     6.999  ,     6.998  ,     7.000  ,     7.001  ,     7.000  ,     6.999  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv (1,1,:) = (/      -5.478  ,    -5.494  ,    -5.499  ,    -5.499  ,    -5.500  ,    -5.500  ,    -5.499  ,    -5.499  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  /)
  
  ! vertical velocity (Pa/s)
  omega (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1107E-01  ,  0.4531E-01  ,  0.5256E-01  ,  0.4205E-01  ,  0.3438E-01  ,  0.2642E-01  ,  0.1962E-01  ,  0.1393E-01  ,  0.9280E-02  ,  0.5578E-02  ,  0.2762E-02  ,  0.7867E-03  /)
  
  ! temperatur at full levels (K)
  tt (1,1,:) = (/     200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   201.015  ,   211.607  ,   221.721  ,   231.309  ,   240.342  ,   248.804  ,   256.671  ,   263.878  ,   270.304  ,   275.817  ,   280.379  ,   284.068  ,   287.028  ,   289.413  ,   291.323  ,   293.154  ,   294.303  ,   283.320  ,   284.734  ,   286.051  ,   287.365  ,   288.417  ,   289.241  ,   289.860  ,   290.290  /)
  
  ! potential temperatur at full levels (K)
  th (1,1,:) = (/    1154.472  ,   869.131  ,   730.368  ,   639.500  ,   573.057  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   322.863  ,   317.933  ,   313.844  ,   310.474  ,   307.685  ,   305.694  ,   303.705  ,   289.884  ,   289.306  ,   289.000  ,   289.000  ,   289.000  ,   289.000  ,   289.000  ,   289.000  /)
  
  ! ice-liquid water potential temperatur at full levels (K)
  thl (1,1,:) = (/    1154.472  ,   869.131  ,   730.368  ,   639.500  ,   573.057  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   322.863  ,   317.933  ,   313.844  ,   310.474  ,   307.685  ,   305.694  ,   303.705  ,   288.769  ,   289.107  ,   289.000  ,   289.000  ,   289.000  ,   289.000  ,   289.000  ,   289.000  /)
  
  ! input relative humidity (%)
  rh (1,1,:) = (/       0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  /)
  
  ! specific humidity at full levels (kg/kg)
  qq (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.4155E-03  ,  0.1097E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.8350E-02  ,  0.8929E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  /)
  
  ! cloud fraction (none)
  qa (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.9726E+00  ,  0.1000E+01  ,  0.2381E-01  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! cloud liquid water mixing ratio at full levels (kg/kg)
  ql (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.4379E-03  ,  0.7878E-04  ,  0.2390E-06  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! cloud ice water mixing ratio at full levels (kg/kg)
  qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! total water mixing ratio (qv+ql+qi) at full levels (kg/kg)
  qt (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.4155E-03  ,  0.1097E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.8788E-02  ,  0.9008E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  ,  0.9000E-02  /)
  
  ! surface air temperature (K)
  t_ref=     290.39
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.9417E-02
  
  ! zonal wind stress
  u_flux =    -0.5958E-01  
  
  ! meridional wind stress
  v_flux =     0.4681E-01  
  
  ! 
  u_star =     0.2500E+00
  b_star =     0.2562E-02
  q_star =     0.1518E-03
  shflx  =     0.1500E+02
  lhflx  =     0.4600E-04
  
  ! rdiag(1,1,:,nQke), input
  rdiag(1,1,:,nQke) =    (/ 0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nSh3D), input
  rdiag(1,1,:,nSh3D) =   (/  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nel_pbl), input
  rdiag(1,1,:,nel_pbl) =  (/   0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00 /) 
  
  ! rdiag(1,1,:,ncldfra_bl), input
  rdiag(1,1,:,ncldfra_bl) =  (/   0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nqc_bl), input
  rdiag(1,1,:,nqc_bl) =   (/  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00 /) 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif (input_profile == "SCM_RF01_rfo76a-M3_EDMFexpUP_NOsm02") then

  ! pressure at half level (Pa)
  p_half (1,1,:) = (/     100.000  ,   400.000  ,   818.602  ,  1378.887  ,  2091.795  ,  2983.641  ,  4121.790  ,  5579.221  ,  7429.322  ,  9739.835  , 12573.187  , 15988.046  , 20044.538  , 24794.785  , 30283.982  , 36517.815  , 43391.604  , 50609.551  , 57745.553  , 64437.098  , 70498.051  , 75875.352  , 80577.543  , 84639.210  , 88112.561  , 91059.544  , 93543.571  , 95617.885  , 97335.022  , 98736.916  , 99856.897  ,100728.562  ,101362.435  ,101780.000  /)
  
  ! pressure at full level (Pa)
  p_full (1,1,:) = (/     216.404  ,   584.531  ,  1074.508  ,  1710.654  ,  2511.381  ,  3522.120  ,  4813.790  ,  6460.179  ,  8532.503  , 11096.287  , 14212.307  , 17939.921  , 22335.536  , 27447.964  , 33303.717  , 39855.967  , 46908.059  , 54099.135  , 61030.197  , 67422.176  , 73153.765  , 78202.888  , 82591.732  , 86364.245  , 89577.974  , 92295.987  , 94576.937  , 96473.907  , 98034.299  , 99295.854  ,100292.098  ,101045.167  ,101571.074  /)
  
  ! actual height at half level (m)
  z_half (1,1,:) = (/   44250.701  , 36129.821  , 31934.737  , 28880.200  , 26438.917  , 24358.628  , 22465.669  , 20692.103  , 19014.469  , 17420.145  , 15837.553  , 14277.171  , 12745.244  , 11248.090  ,  9790.709  ,  8383.510  ,  7050.533  ,  5832.289  ,  4766.669  ,  3866.252  ,  3118.106  ,  2499.723  ,  1989.549  ,  1569.540  ,  1223.917  ,   941.031  ,   714.963  ,   529.891  ,   379.693  ,   258.611  ,   162.780  ,    88.775  ,    35.183  ,     0.000  /)
  
  ! actual height at full level (m)
  z_full (1,1,:) = (/   39280.838  , 33784.040  , 30275.338  , 27575.020  , 25337.339  , 23361.262  , 21534.207  , 19813.303  , 18181.374  , 16595.211  , 15026.149  , 13482.366  , 11970.153  , 10495.128  ,  9065.172  ,  7697.873  ,  6425.795  ,  5287.769  ,  4308.235  ,  3486.576  ,  2805.127  ,  2242.079  ,  1777.823  ,  1395.570  ,  1081.699  ,   827.490  ,   622.089  ,   454.569  ,   319.008  ,   210.606  ,   125.724  ,    61.951  ,    17.580  /)
  
  ! height at half level above the surface (m)
  z_half_surf0 (1,1,:) = (/   44250.701  , 36129.821  , 31934.737  , 28880.200  , 26438.917  , 24358.628  , 22465.669  , 20692.103  , 19014.469  , 17420.145  , 15837.553  , 14277.171  , 12745.244  , 11248.090  ,  9790.709  ,  8383.510  ,  7050.533  ,  5832.289  ,  4766.669  ,  3866.252  ,  3118.106  ,  2499.723  ,  1989.549  ,  1569.540  ,  1223.917  ,   941.031  ,   714.963  ,   529.891  ,   379.693  ,   258.611  ,   162.780  ,    88.775  ,    35.183  ,     0.000  /)
  
  ! height at full level above the surface (m)
  z_full_surf0 (1,1,:) = (/   39280.838  , 33784.040  , 30275.338  , 27575.020  , 25337.339  , 23361.262  , 21534.207  , 19813.303  , 18181.374  , 16595.211  , 15026.149  , 13482.366  , 11970.153  , 10495.128  ,  9065.172  ,  7697.873  ,  6425.795  ,  5287.769  ,  4308.235  ,  3486.576  ,  2805.127  ,  2242.079  ,  1777.823  ,  1395.570  ,  1081.699  ,   827.490  ,   622.089  ,   454.569  ,   319.008  ,   210.606  ,   125.724  ,    61.951  ,    17.580  /)
  
  ! zonal wind velocity at full levels (m/s)
  uu (1,1,:) = (/       6.934  ,     6.985  ,     6.989  ,     6.997  ,     6.995  ,     6.999  ,     7.003  ,     7.000  ,     6.998  ,     7.000  ,     7.000  ,     6.999  ,     7.000  ,     7.001  ,     7.001  ,     6.999  ,     6.999  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     7.000  ,     6.982  ,     8.648  ,     6.721  ,     6.596  ,     6.595  ,     6.419  ,     6.279  ,     6.192  ,     5.492  ,     5.225  /)
  
  ! meridional wind velocity at full levels (m/s)
  vv (1,1,:) = (/      -5.423  ,    -5.481  ,    -5.494  ,    -5.498  ,    -5.498  ,    -5.499  ,    -5.497  ,    -5.497  ,    -5.499  ,    -5.500  ,    -5.501  ,    -5.501  ,    -5.500  ,    -5.499  ,    -5.499  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.500  ,    -5.488  ,    -7.168  ,    -5.163  ,    -5.091  ,    -5.082  ,    -4.895  ,    -4.704  ,    -4.751  ,    -3.906  ,    -3.630  /)
  
  ! vertical velocity (Pa/s)
  omega (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1093E-01  ,  0.4517E-01  ,  0.5264E-01  ,  0.4232E-01  ,  0.3424E-01  ,  0.2627E-01  ,  0.1959E-01  ,  0.1393E-01  ,  0.9281E-02  ,  0.5583E-02  ,  0.2764E-02  ,  0.7870E-03  /)
  
  ! temperatur at full levels (K)
  tt (1,1,:) = (/     200.001  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   200.000  ,   201.015  ,   211.607  ,   221.721  ,   231.309  ,   240.342  ,   248.804  ,   256.671  ,   263.878  ,   270.304  ,   275.817  ,   280.379  ,   284.068  ,   287.028  ,   289.427  ,   291.323  ,   293.158  ,   293.131  ,   285.770  ,   286.641  ,   286.544  ,   287.506  ,   288.483  ,   289.147  ,   289.981  ,   290.463  /)
  
  ! potential temperatur at full levels (K)
  th (1,1,:) = (/    1154.475  ,   869.132  ,   730.368  ,   639.501  ,   573.057  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   322.863  ,   317.933  ,   313.844  ,   310.490  ,   307.686  ,   305.698  ,   302.495  ,   292.391  ,   291.244  ,   289.499  ,   289.141  ,   289.066  ,   288.907  ,   289.121  ,   289.172  /)
  
  ! ice-liquid water potential temperatur at full levels (K)
  thl (1,1,:) = (/    1154.475  ,   869.132  ,   730.368  ,   639.501  ,   573.057  ,   520.270  ,   475.841  ,   437.482  ,   406.102  ,   396.585  ,   387.170  ,   377.908  ,   368.834  ,   359.984  ,   351.405  ,   343.201  ,   335.568  ,   328.740  ,   322.863  ,   317.933  ,   313.844  ,   310.490  ,   307.686  ,   305.698  ,   302.494  ,   292.391  ,   291.244  ,   289.498  ,   289.141  ,   289.066  ,   288.907  ,   289.121  ,   289.172  /)
  
  ! input relative humidity (%)
  rh (1,1,:) = (/       0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  ,     0.000  /)
  
  ! specific humidity at full levels (kg/kg)
  qq (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.1573E-13  ,  0.4155E-03  ,  0.1097E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1504E-02  ,  0.2547E-02  ,  0.6953E-02  ,  0.7849E-02  ,  0.8801E-02  ,  0.8995E-02  ,  0.9091E-02  ,  0.8897E-02  ,  0.9592E-02  ,  0.9776E-02  /)
  
  ! cloud fraction (none)
  qa (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.5811E-02  ,  0.1644E-03  ,  0.0000E+00  ,  0.6072E-01  ,  0.1316E-02  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! cloud liquid water mixing ratio at full levels (kg/kg)
  ql (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.6791E-06  ,  0.1921E-07  ,  0.0000E+00  ,  0.2225E-07  ,  0.4823E-09  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! cloud ice water mixing ratio at full levels (kg/kg)
  qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! total water mixing ratio (qv+ql+qi) at full levels (kg/kg)
  qt (1,1,:) = (/    0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.5000E-15  ,  0.1573E-13  ,  0.4155E-03  ,  0.1097E-02  ,  0.1500E-02  ,  0.1500E-02  ,  0.1504E-02  ,  0.2548E-02  ,  0.6953E-02  ,  0.7849E-02  ,  0.8801E-02  ,  0.8995E-02  ,  0.9091E-02  ,  0.8897E-02  ,  0.9592E-02  ,  0.9776E-02  /)
  
  ! surface air temperature (K)
  t_ref=     290.54
  
  ! surface air specific humidity (kg/kg)
  q_ref =     0.1008E-01
  
  ! zonal wind stress
  u_flux =    -0.5969E-01  
  
  ! meridional wind stress
  v_flux =     0.4654E-01  
  
  ! 
  u_star =     0.2500E+00
  b_star =     0.2563E-02
  q_star =     0.1519E-03
  shflx  =     0.1500E+02
  lhflx  =     0.4600E-04
  
  ! rdiag(1,1,:,nQke), input
  rdiag(1,1,:,nQke) = (/    0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.1000E-03  ,  0.6316E-02  ,  0.1780E+00  ,  0.2442E+00  ,  0.1602E+00  ,  0.1151E+00  ,  0.2901E+00  ,  0.1206E+01  ,  0.9146E+00  ,  0.7850E+00  /)
  
  ! rdiag(1,1,:,nSh3D), input
  rdiag(1,1,:,nSh3D) = (/    0.4390E-07  ,  0.5314E-07  ,  0.6152E-07  ,  0.6922E-07  ,  0.7656E-07  ,  0.8394E-07  ,  0.9151E-07  ,  0.1061E-06  ,  0.3395E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3396E-06  ,  0.3397E-06  ,  0.3400E-06  ,  0.3402E-06  ,  0.3402E-06  ,  0.3457E-06  ,  0.3510E-06  ,  0.3490E-06  ,  0.3383E-06  ,  0.6319E-06  ,  0.4471E-02  ,  0.3989E-02  ,  0.3827E+00  ,  0.2752E-01  ,  0.7252E-01  ,  0.1029E+01  ,  0.1672E-01  ,  0.1545E+00  ,  0.8839E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nel_pbl), input
  rdiag(1,1,:,nel_pbl) = (/    0.2437E+00  ,  0.2682E+00  ,  0.2885E+00  ,  0.3060E+00  ,  0.3218E+00  ,  0.3370E+00  ,  0.3518E+00  ,  0.3788E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6774E+00  ,  0.6776E+00  ,  0.6778E+00  ,  0.6780E+00  ,  0.6781E+00  ,  0.6835E+00  ,  0.6887E+00  ,  0.6868E+00  ,  0.6762E+00  ,  0.7267E+00  ,  0.1394E+01  ,  0.7366E+01  ,  0.2707E+02  ,  0.6214E+02  ,  0.5889E+02  ,  0.1017E+02  ,  0.4586E+02  ,  0.3348E+02  ,  0.1638E+02  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,ncldfra_bl), input
  rdiag(1,1,:,ncldfra_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1653E+00  ,  0.1551E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1582E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
  ! rdiag(1,1,:,nqc_bl), input
  rdiag(1,1,:,nqc_bl) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.4883E-04  ,  0.2462E-04  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.7686E-04  ,  0.0000E+00  ,  0.0000E+00  /)
 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!elseif (input_profile == "") then

else
  print*,'ERROR: unsupported input_profile, ', input_profile
  stop
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
endif ! end if of input profile  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
!==============================
!==============================
! run edmf_mynn
!==============================
!==============================

  print*,'------------------------'
  print*,'hello world!'
  print*,'input profile   ,', input_profile
  print*,'bl_mynn_edmf    ,', bl_mynn_edmf
  print*,'bl_mynn_stabfunc,', bl_mynn_stabfunc
  print*,'expmf           ,', expmf
  print*,'upwind          ,', upwind
  print*,'------------------------'
  print*,''

  area = 1.e+10
  Physics_input_block%p_full = p_full
  Physics_input_block%p_half = p_half
  Physics_input_block%z_full = z_full
  Physics_input_block%z_half = z_half
  Physics_input_block%t = tt
  Physics_input_block%u = uu
  Physics_input_block%v = vv
  Physics_input_block%omega = omega
  Physics_input_block%q(:,:,:,nsphum) = qq
  Physics_input_block%q(:,:,:,nql) = ql
  Physics_input_block%q(:,:,:,nqi) = qi
  Physics_input_block%q(:,:,:,nqa) = qa

  !Physics_input_block%q(:,:,:,nql) = 1.e-4
  !Physics_input_block%q(:,:,:,nqi) = 0.
  !Physics_input_block%q(:,:,:,nqa) = 1.

  !Physics_input_block%q(:,:,:,nql) = 0.
  !Physics_input_block%q(:,:,:,nqi) = 0.
  !Physics_input_block%q(:,:,:,nqa) = 0.

  do mm = 1,loop_times
    ! set tendencies to zeros
    udt = 0. ; vdt = 0. ; tdt = 0.; rdt = 0.   

!print*,'tt',Physics_input_block%t
!print *,'qq',qq

    print*,''
    print*,'----------------------'
    print*, 'loop times=',mm
    print*,'----------------------'
!    call edmf_mynn_driver ( &
!              is, ie, js, je, npz, Time_next, dt, lon, lat, frac_land, area, u_star,  &
!              b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, & 
!              do_edmf_mynn_diagnostic, &
!              option_edmf2ls_mp, qadt_edmf, qldt_edmf, qidt_edmf, dqa_edmf,  dql_edmf, dqi_edmf, &
!              pbltop, udt, vdt, tdt, rdt, rdiag(:,:,:,ntp+1:))
!              !udt, vdt, tdt, rdt, rdiag)

!write(6,*) 'nQke, rdiag(:,:,:,nQke)',nQke, rdiag(:,:,:,nQke)
!write(6,*) 'rdiag(:,:,:,nel_pbl)',rdiag(:,:,:,nel_pbl)
!write(6,*) 'rdiag(:,:,:,ncldfra_bl)',rdiag(:,:,:,nqc_bl)
!write(6,*) 'rdiag(:,:,:,nqc_bl)',rdiag(:,:,:,nqc_bl)
!write(6,*) 'rdiag(:,:,:,nSh3D)',rdiag(:,:,:,nSh3D)
!write(6,*) 'ntp,nQke, nSh3D, nel_pbl, ncldfra_bl, nqc_bl',ntp,nQke, nSh3D, nel_pbl, ncldfra_bl, nqc_bl

  !rdt_mynn_ed_am4(:,:,:,nql) = 1.e-9
  !print*,'rdt_mynn_ed_am4',rdt_mynn_ed_am4(:,:,:,nql)
  !rdt_mynn_ed_am4(:,:,:,:) = 0.
  !rdt_mynn_ed_am4(:,:,:,nql) = 1.
  !rdt_mynn_ed_am4(:,:,:,nqa) = 2.
  !rdt_mynn_ed_am4(:,:,:,nqi) = 3.

  !call edmf_mynn_driver ( &
  !            is, ie, js, je, npz, Time_next, dt, lon, lat, frac_land, area, u_star,  &
  !            b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, &
  !            rdt_mynn_ed_am4,  & 
  !            do_edmf_mynn_diagnostic, do_return_edmf_mynn_diff_only, do_edmf_mynn_in_physics, &
  !            option_edmf2ls_mp, qadt_edmf, qldt_edmf, qidt_edmf, dqa_edmf,  dql_edmf, dqi_edmf, diff_t_edmf, diff_m_edmf, kpbl_edmf, &
  !            edmf_mc_full, edmf_mc_half, edmf_moist_area, edmf_dry_area, edmf_moist_humidity, edmf_dry_humidity, &
  !           pbltop, udt, vdt, tdt, rdt, rdiag)

  call edmf_mynn_driver ( &
              is, ie, js, je, npz, Time_next, dt, lon, lat, frac_land, area, u_star,  &
              b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, &
              rdt_mynn_ed_am4, &
              do_edmf_mynn_diagnostic, do_return_edmf_mynn_diff_only, do_edmf_mynn_in_physics, do_tracers_selective, &
              option_edmf2ls_mp, qadt_edmf, qldt_edmf, qidt_edmf, dqa_edmf,  dql_edmf, dqi_edmf, diff_t_edmf, diff_m_edmf, kpbl_edmf, &
              edmf_mc_full, edmf_mc_half, edmf_moist_area, edmf_dry_area, edmf_moist_humidity, edmf_dry_humidity, &
              pbltop, udt, vdt, tdt, rdt, rdiag)

              !is, ie, js, je, npz, Time_next, dt, lon, lat, frac_land, area, u_star,  &
              !b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, &
              !do_edmf_mynn_diagnostic, &
              !option_edmf2ls_mp, qadt_edmf, qldt_edmf, qidt_edmf, dqa_edmf,  dql_edmf, dqi_edmf, diff_t_edmf, diff_m_edmf, kpbl_edmf, &
              !pbltop, udt, vdt, tdt, rdt, rdiag(:,:,:,ntp+1:))

    ! update fields
    Physics_input_block%u = Physics_input_block%u + udt(:,:,:)*dt
    Physics_input_block%v = Physics_input_block%v + vdt(:,:,:)*dt
    Physics_input_block%t = Physics_input_block%t + tdt(:,:,:)*dt
    Physics_input_block%q(:,:,:,nsphum) = Physics_input_block%q(:,:,:,1) + rdt(:,:,:,nsphum)*dt



  enddo


  area = 0.
  area = 1/area

end program test111


