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
  !character*50 :: input_profile = "SCM_am4p0_DCBL_C1_01"
  !character*50 :: input_profile = "SCM_am4p0_DCBL_C1_02_u,vdt_NaN"
  !character*50 :: input_profile = "SCM_am4p0_BOMEX_01"
  character*50 :: input_profile = "AMIP_i27_j01_IndOcn"
  !character*50 :: input_profile = "xxx"

  integer, parameter :: loop_times = 1
 ! integer, parameter :: loop_times = 100
 ! integer, parameter :: loop_times = 60 

  integer, parameter :: ni = 1
  integer, parameter :: nj = 1
  integer, parameter :: npz = 33
  integer, parameter :: nfull = npz
  integer, parameter :: nhalf = nfull+1
  integer, parameter :: nfull2 = nfull+1
  integer, parameter :: nhalf2 = nfull2+1
  integer, parameter :: tr = 4

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

  real, parameter :: dt = 1800.
  !real, parameter :: dt = 300.
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
   INTEGER :: bl_mynn_edmf      = 3        ! controls  the version of the EDMF to be called 3=JPLedmf documented by Wu et al., 2020
                                            ! set “bl_mynn_edmf=0” and “bl_mynn_edmf_dd<>1” then the scheme will be MYNN only.
   INTEGER :: bl_mynn_edmf_dd   = 0         ! 0 - no downdrafts, 1 - Wu et al., 2020 downdrafts 
   REAL    :: bl_mynn_edmf_Lent = 0.        ! dummy argument
   INTEGER :: bl_mynn_edmf_mom  = 1         ! 1- mixing of momentum by full EDMF, 0 - momenum is mixed only by ED
   INTEGER :: bl_mynn_edmf_tke  = 0         ! does not effect JPLedmf
   INTEGER :: bl_mynn_edmf_part = 1         ! partition area into updrafts and remaining environment
   INTEGER :: bl_mynn_cloudmix  = 1         ! 1 - cloud species are mixed separately, 0 - not
   INTEGER :: bl_mynn_mixqt     = 2         ! will mix moist conserved variables, after mixing invokes PDF cloud scheme to convert moist variables to dry
   INTEGER :: icloud_bl         = 1         ! 1, cloud cover and liquid water from the PDF scheme will be on the cloud output

   integer :: initflag = 1 
   logical :: FLAG_QI  = .false.             ! (flags for whether cloud and ice mixing rations and number concentrations are mixed separately)
   logical :: FLAG_QNI = .false.            ! all false
   logical :: FLAG_QC  = .false.
   logical :: FLAG_QNC = .false.


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
   logical :: do_writeout_column_nml = .true.
   !logical :: do_writeout_column_nml = .false.
   !logical :: do_edmf_mynn_diagnostic = .true.
   logical :: do_edmf_mynn_diagnostic = .false.

   real    :: tdt_max     = 500. ! K/day
   logical :: do_limit_tdt = .false.
   real    :: tdt_limit =   200. ! K/day

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
    dz, u, v, w, th, qv, p, exner, rho, T3D, qc, ql, qi, qnc, qni

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
    RUBLTEN, RVBLTEN, RTHBLTEN, RQVBLTEN, RQCBLTEN, RQIBLTEN, RQNIBLTEN, RTHRATEN, &
    edmf_a, edmf_w, edmf_qt, edmf_thl, edmf_ent, edmf_qc, edmf_a_dd,edmf_w_dd,edmf_qt_dd,edmf_thl_dd,edmf_ent_dd,edmf_qc_dd, &
    edmf_debug1,edmf_debug2,edmf_debug3,edmf_debug4, &
    mynn_ql, qc_bl, cldfra_bl, el_pbl, Sh3D

  real, dimension(:,:),   allocatable :: &   ! OUTPUT, DIMENSION(IMS:IME,JMS:JME)
    maxmf

  real, dimension(:,:,:), allocatable :: &   ! OUTPUT, DIMENSION(IMS:IME,KMS:KME,JMS:JME)
    exch_h, exch_m, &
    qWT, qSHEAR, qBUOY, qDISS, dqke

end type edmf_output_type

!==================
type am4_edmf_output_type

  real, dimension(:,:,:),   allocatable :: &   ! OUTPUT, DIMENSION(nlon, nlat, nfull)
    thl_edmf, qt_edmf,  & ! diagnostic purpose
    tke, Tsq, Cov_thl_qt, udt_edmf, vdt_edmf, tdt_edmf, qdt_edmf, qidt_edmf, qcdt_edmf, &  ! outputs from EDMF-MYNN scheme
    edmf_a, edmf_w, edmf_qt, edmf_thl, edmf_ent, edmf_qc  

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

!!    qke(kts)=qke(kts+1)
!!    tsq(kts)=tsq(kts+1)
!!    qsq(kts)=qsq(kts+1)
!!    cov(kts)=cov(kts+1)

    qke(kte)=qke(kte-1)
    tsq(kte)=tsq(kte-1)
    qsq(kte)=qsq(kte-1)
    cov(kte)=cov(kte-1)

!
!    RETURN

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
       sh (k) = shc*( rfc-rf )/( 1.0-rf )
       sm (k) = smc*( rf1-rf )/( rf2-rf ) * sh(k)
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
    &            Sh, el, bl_mynn_cloudpdf,&
    &            qc_bl1D, cldfra_bl1D,    &
    &            PBLH1,HFX1,              &
    &            Vt, Vq, th, sgm)

!-------------------------------------------------------------------

    INTEGER, INTENT(IN)   :: kts,kte, bl_mynn_cloudpdf


    REAL, INTENT(IN)      :: dx,PBLH1,HFX1
    REAL, DIMENSION(kts:kte), INTENT(IN) :: dz
    REAL, DIMENSION(kts:kte), INTENT(IN) :: p,exner, thl, qw, &
         &tsq, qsq, cov, th

    REAL, DIMENSION(kts:kte), INTENT(INOUT) :: vt,vq,sgm

    REAL, DIMENSION(kts:kte) :: qmq,alp,a,bet,b,ql,q1,cld,RH
    REAL, DIMENSION(kts:kte), INTENT(OUT) :: qc_bl1D,cldfra_bl1D
    DOUBLE PRECISION :: t3sq, r3sq, c3sq

    REAL :: qsl,esat,qsat,tlk,qsat_tl,dqsl,cld0,q1k,eq1,qll,&
         &q2p,pt,rac,qt,t,xl,rsl,cpm,cdhdz,Fng,qww,alpha,beta,bb,ls_min,ls,wt
    INTEGER :: i,j,k

    REAL :: erf

    !JOE: NEW VARIABLES FOR ALTERNATE SIGMA
    REAL::dth,dtl,dqw,dzk
    REAL, DIMENSION(kts:kte), INTENT(IN) :: Sh,el

    !JOE: variables for BL clouds
    REAL::zagl,cld9,damp,edown,RHcrit,RHmean,RHsum,RHnum,Hshcu,PBLH2,ql_limit
    REAL, PARAMETER :: Hfac = 3.0     !cloud depth factor for HFX (m^3/W)
    REAL, PARAMETER :: HFXmin = 50.0  !min W/m^2 for BL clouds
    REAL            :: RH_00L, RH_00O, phi_dz, lfac
    REAL, PARAMETER :: cdz = 2.0
    REAL, PARAMETER :: mdz = 1.5

    !JAYMES:  variables for tropopause-height estimation
    REAL            :: theta1, theta2, ht1, ht2
    INTEGER         :: k_tropo

! First, obtain an estimate for the tropopause height (k), using the method employed in the
! Thompson subgrid-cloud scheme.  This height will be a consideration later when determining 
! the "final" subgrid-cloud properties.
! JAYMES:  added 3 Nov 2016, adapted from G. Thompson

    DO k = kte-3, kts, -1
       theta1 = th(k)
       theta2 = th(k+2)
       ht1 = 44307.692 * (1.0 - (p(k)/101325.)**0.190)
       ht2 = 44307.692 * (1.0 - (p(k+2)/101325.)**0.190)
       if ( (((theta2-theta1)/(ht2-ht1)) .lt. 10./1500. ) .AND.       &
     &                       (ht1.lt.19000.) .and. (ht1.gt.4000.) ) then 
          goto 86
       endif
    ENDDO
 86   continue
    k_tropo = MAX(kts+2, k+2)

    zagl = 0.

    SELECT CASE(bl_mynn_cloudpdf)

      CASE (0) ! ORIGINAL MYNN PARTIAL-CONDENSATION SCHEME

        DO k = kts,kte-1
           t  = th(k)*exner(k)

           !SATURATED VAPOR PRESSURE
           esat = esat_blend(t)
           !SATURATED SPECIFIC HUMIDITY
           qsl=ep_2*esat/(p(k)-ep_3*esat)
           !dqw/dT: Clausius-Clapeyron
           dqsl = qsl*ep_2*ev/( rd*t**2 )
           !RH (0 to 1.0)
           RH(k)=MAX(MIN(1.0,qw(k)/MAX(1.E-8,qsl)),0.001)

           alp(k) = 1.0/( 1.0+dqsl*xlvcp )
           bet(k) = dqsl*exner(k)

           !NOTE: negative bl_mynn_cloudpdf will zero-out the stratus subgrid clouds
           !      at the end of this subroutine. 
           !Sommeria and Deardorff (1977) scheme, as implemented
           !in Nakanishi and Niino (2009), Appendix B
           t3sq = MAX( tsq(k), 0.0 )
           r3sq = MAX( qsq(k), 0.0 )
           c3sq =      cov(k)
           c3sq = SIGN( MIN( ABS(c3sq), SQRT(t3sq*r3sq) ), c3sq )
           r3sq = r3sq +bet(k)**2*t3sq -2.0*bet(k)*c3sq
           !DEFICIT/EXCESS WATER CONTENT
           qmq(k) = qw(k) -qsl
           !ORIGINAL STANDARD DEVIATION: limit e-6 produces ~10% more BL clouds
           !than e-10
           sgm(k) = SQRT( MAX( r3sq, 1.0d-10 ))
           !NORMALIZED DEPARTURE FROM SATURATION
           q1(k)   = qmq(k) / sgm(k)
           !CLOUD FRACTION. rr2 = 1/SQRT(2) = 0.707
           cld(k) = 0.5*( 1.0+erf( q1(k)*rr2 ) )

        END DO

      CASE (1, -1) !ALTERNATIVE FORM (Nakanishi & Niino 2004 BLM, eq. B6, and
                       !Kuwano-Yoshida et al. 2010 QJRMS, eq. 7):
        DO k = kts,kte-1
           t  = th(k)*exner(k)
           !SATURATED VAPOR PRESSURE
           esat = esat_blend(t)
           !SATURATED SPECIFIC HUMIDITY
           qsl=ep_2*esat/(p(k)-ep_3*esat)
           !dqw/dT: Clausius-Clapeyron
           dqsl = qsl*ep_2*ev/( rd*t**2 )
           !RH (0 to 1.0)
           RH(k)=MAX(MIN(1.0,qw(k)/MAX(1.E-8,qsl)),0.001)

           alp(k) = 1.0/( 1.0+dqsl*xlvcp )
           bet(k) = dqsl*exner(k)

           if (k .eq. kts) then 
             dzk = 0.5*dz(k)
           else
             dzk = 0.5*( dz(k) + dz(k-1) )
           end if
           dth = 0.5*(thl(k+1)+thl(k)) - 0.5*(thl(k)+thl(MAX(k-1,kts)))
           dqw = 0.5*(qw(k+1) + qw(k)) - 0.5*(qw(k) + qw(MAX(k-1,kts)))
           sgm(k) = SQRT( MAX( (alp(k)**2 * MAX(el(k)**2,0.1) * &
                             b2 * MAX(Sh(k),0.03))/4. * &
                      (dqw/dzk - bet(k)*(dth/dzk ))**2 , 1.0e-8) ) ! set the limit smaller
           sgm(k) = min(sgm(k),1.0e-4) ! set upper limit of sigma
           qmq(k) = qw(k) -qsl
           q1(k)   = qmq(k) / sgm(k)
           cld(k) = 0.5*( 1.0+erf( q1(k)*rr2 ) )
        END DO

      CASE (2, -2)
          !Diagnostic statistical scheme of Chaboureau and Bechtold (2002), JAS
          !JAYMES- this added 27 Apr 2015
        DO k = kts,kte-1
           t  = th(k)*exner(k)
           !SATURATED VAPOR PRESSURE
           esat = esat_blend(t)
           !SATURATED SPECIFIC HUMIDITY
           qsl=ep_2*esat/(p(k)-ep_3*esat)
           !dqw/dT: Clausius-Clapeyron
           dqsl = qsl*ep_2*ev/( rd*t**2 )
           !RH (0 to 1.0)
           RH(k)=MAX(MIN(1.0,qw(k)/MAX(1.E-8,qsl)),0.001)

           alp(k) = 1.0/( 1.0+dqsl*xlvcp )
           bet(k) = dqsl*exner(k)

           xl = xl_blend(t)                    ! obtain latent heat

           tlk = thl(k)*(p(k)/p1000mb)**rcp    ! recover liquid temp (tl) from thl

           qsat_tl = qsat_blend(tlk,p(k))      ! get saturation water vapor mixing ratio
                                               !   at tl and p

           rsl = xl*qsat_tl / (r_v*tlk**2)     ! slope of C-C curve at t = tl
                                               ! CB02, Eqn. 4
 
           cpm = cp + qw(k)*cpv                ! CB02, sec. 2, para. 1
     
           a(k) = 1./(1. + xl*rsl/cpm)         ! CB02 variable "a"

           qmq(k) = a(k) * (qw(k) - qsat_tl) ! saturation deficit/excess;
                                               !   the numerator of Q1

           b(k) = a(k)*rsl                     ! CB02 variable "b"

           dtl =    0.5*(thl(k+1)*(p(k+1)/p1000mb)**rcp + tlk) &
               & - 0.5*(tlk + thl(MAX(k-1,kts))*(p(MAX(k-1,kts))/p1000mb)**rcp)

           dqw = 0.5*(qw(k+1) + qw(k)) - 0.5*(qw(k) + qw(MAX(k-1,kts)))

           if (k .eq. kts) then
             dzk = 0.5*dz(k)
           else
             dzk = 0.5*( dz(k) + dz(k-1) )
           end if

           cdhdz = dtl/dzk + (g/cpm)*(1.+qw(k))  ! expression below Eq. 9
                                                 ! in CB02

           zagl = zagl + dz(k)

           ls_min = 400. + MIN(3.*MAX(HFX1,0.),500.)
           ls_min = MIN(MAX(zagl,25.),ls_min) ! Let this be the minimum possible length scale:
           if (zagl > PBLH1+2000.) ls_min = MAX(ls_min + 0.5*(PBLH1+2000.-zagl),400.)
                                        !   25 m < ls_min(=zagl) < 300 m
           lfac=MIN(4.25+dx/4000.,6.)   ! A dx-dependent multiplier for the master length scale:
                                        !   lfac(750 m) = 4.4
                                        !   lfac(3 km)  = 5.0
                                        !   lfac(13 km) = 6.0

           ls = MAX(MIN(lfac*el(k),900.),ls_min)  ! Bounded:  ls_min < ls < 900 m
                   ! Note: CB02 use 900 m as a constant free-atmosphere length scale. 

                   ! Above 300 m AGL, ls_min remains 300 m.  For dx = 3 km, the 
                   ! MYNN master length scale (el) must exceed 60 m before ls
                   ! becomes responsive to el, otherwise ls = ls_min = 300 m.

           sgm(k) = MAX(1.e-10, 0.225*ls*SQRT(MAX(0., & ! Eq. 9 in CB02:
                   & (a(k)*dqw/dzk)**2              & ! < 1st term in brackets,
                   & -2*a(k)*b(k)*cdhdz*dqw/dzk     & ! < 2nd term,
                   & +b(k)**2 * cdhdz**2)))           ! < 3rd term
                   ! CB02 use a multiplier of 0.2, but 0.225 is chosen
                   ! based on tests

           q1(k) = qmq(k) / sgm(k)  ! Q1, the normalized saturation

           cld(k) = MAX(0., MIN(1., 0.5+0.36*ATAN(1.55*q1(k)))) ! Eq. 7 in CB02

         END DO

    END SELECT

    zagl = 0.
    RHsum=0.
    RHnum=0.
    RHmean=0.1 !initialize with small value for small PBLH cases
    damp =0
    PBLH2=MAX(10.,PBLH1)

    SELECT CASE(bl_mynn_cloudpdf)

      CASE (-1 : 1) ! ORIGINAL MYNN PARTIAL-CONDENSATION SCHEME
                    ! OR KUWANO ET AL.
        DO k = kts,kte-1
           t    = th(k)*exner(k)
           q1k  = q1(k)
           zagl = zagl + dz(k)
           !q1=0.
           !cld(k)=0.

           !COMPUTE MEAN RH IN PBL (NOT PRESSURE WEIGHTED).
           IF (zagl < PBLH2 .AND. PBLH2 > 400.) THEN
              RHsum=RHsum+RH(k)
              RHnum=RHnum+1.0
              RHmean=RHsum/RHnum
           ENDIF

           RHcrit = 1. - 0.35*(1.0 - (MAX(250.- MAX(HFX1,HFXmin),0.0)/200.)**2)
           if (HFX1 > HFXmin) then
              cld9=MIN(MAX(0., (rh(k)-RHcrit)/(1.1-RHcrit)), 1.)**2
           else
              cld9=0.0
           endif
       
           edown=PBLH2*.1
           !Vary BL cloud depth (Hshcu) by mean RH in PBL and HFX 
           !(somewhat following results from Zhang and Klein (2013, JAS))
           Hshcu=200. + (RHmean+0.5)**1.5*MAX(HFX1,0.)*Hfac
           if (zagl < PBLH2-edown) then
              damp=MIN(1.0,exp(-ABS(((PBLH2-edown)-zagl)/edown)))
           elseif(zagl >= PBLH2-edown .AND. zagl < PBLH2+Hshcu)then
              damp=1.
           elseif (zagl >= PBLH2+Hshcu)then
              damp=MIN(1.0,exp(-ABS((zagl-(PBLH2+Hshcu))/500.)))
           endif
           ! cldfra_bl1D(k)=cld9*damp
           cldfra_bl1D(k)=cld(k) ! JAYMES: use this form to retain the Sommeria-Deardorff value
       
           !use alternate cloud fraction to estimate qc for use in BL clouds-radiation
           eq1  = rrp*EXP( -0.5*q1k*q1k )
           qll  = MAX( cldfra_bl1D(k)*q1k + eq1, 0.0 )
           !ESTIMATED LIQUID WATER CONTENT (UNNORMALIZED)
           ql (k) = alp(k)*sgm(k)*qll
           ! print*, "[MYNN condensation]: k = , ", k, " alp = ", alp(k), " sgm = ", sgm(k), " qll = ", qll, " ql = ", ql(k)
           if(cldfra_bl1D(k)>0.01 .and. ql(k)<1.E-6)ql(k)=1.E-6
           ! qc_bl1D(k)=ql(k)*damp
           qc_bl1D(k)=ql(k) ! JAYMES: use this form to retain the Sommeria-Deardorff value
       
           !now recompute estimated lwc for PBL scheme's use
           !qll IS THE NORMALIZED LIQUID WATER CONTENT (Sommeria and
           !Deardorff (1977, eq 29a). rrp = 1/(sqrt(2*pi)) = 0.3989
           eq1  = rrp*EXP( -0.5*q1k*q1k )
           qll  = MAX( cld(k)*q1k + eq1, 0.0 )
           !ESTIMATED LIQUID WATER CONTENT (UNNORMALIZED)
           ql (k) = alp(k)*sgm(k)*qll
       
           q2p = xlvcp/exner(k)
           pt = thl(k) +q2p*ql(k) ! potential temp
       
           !qt is a THETA-V CONVERSION FOR TOTAL WATER (i.e., THETA-V = qt*THETA)
           qt   = 1.0 +p608*qw(k) -(1.+p608)*ql(k)
           rac  = alp(k)*( cld(k)-qll*eq1 )*( q2p*qt-(1.+p608)*pt )
       
           !BUOYANCY FACTORS: wherever vt and vq are used, there is a
           !"+1" and "+tv0", respectively, so these are subtracted out here.
           !vt is unitless and vq has units of K.
           vt(k) =      qt-1.0 -rac*bet(k)
           vq(k) = p608*pt-tv0 +rac

           !To avoid FPE in radiation driver, when these two quantities are multiplied by eachother,
           ! add limit to qc_bl and cldfra_bl:
           IF (QC_BL1D(k) < 1E-6 .AND. ABS(CLDFRA_BL1D(k)) > 0.01) QC_BL1D(k)= 1E-6
           IF (CLDFRA_BL1D(k) < 1E-2)THEN
              CLDFRA_BL1D(k)=0.
              QC_BL1D(k)=0.
           ENDIF
           IF (t<273.) THEN
              CLDFRA_BL1D(k)=0.
              QC_BL1D(k)=0.               
           ENDIF

        END DO
      CASE ( 2, -2)
        ! JAYMES- this option added 8 May 2015
        ! The cloud water formulations are taken from CB02, Eq. 8.
        ! "fng" represents the non-Gaussian contribution to the liquid
        ! water flux; these formulations are from Cuijpers and Bechtold
        ! (1995), Eq. 7.  CB95 also draws from Bechtold et al. 1995,
        ! hereafter BCMT95
        DO k = kts,kte-1
           t    = th(k)*exner(k)
           q1k  = q1(k)
           zagl = zagl + dz(k)
           IF (q1k < 0.) THEN                 
              ql (k) = sgm(k)*EXP(1.2*q1k-1)
           ELSE IF (q1k > 2.) THEN
              ql (k) = sgm(k)*q1k
           ELSE
              ql (k) = sgm(k)*(EXP(-1.) + 0.66*q1k + 0.086*q1k**2)
           ENDIF
  
           !Next, adjust our initial estimates of cldfra and ql based
           !on tropopause-height and PBLH considerations
           !JAYMES:  added 4 Nov 2016
           if ((cld(k) .gt. 0.) .or. (ql(k) .gt. 0.))  then
              if (k .le. k_tropo) then
                 !At and below tropopause: impose an upper limit on ql; assume that
                 !a maximum of 0.5 percent supersaturation in water vapor can be
                 !available for cloud production
                 ql_limit = 0.005 * qsat_blend( th(k)*exner(k), p(k) )
                 ql(k) = MIN( ql(k), ql_limit )
              else
                 !Above tropopause:  eliminate subgrid clouds from CB scheme
                 cld(k) = 0.
                 ql(k) = 0.
              endif 
           endif
       
           !Buoyancy-flux-related calculations follow...
           ! "Fng" represents the non-Gaussian transport factor
           ! (non-dimensional) from from Bechtold et al. 1995 
           ! (hereafter BCMT95), section 3(c).  Their suggested 
           ! forms for Fng (from their Eq. 20) are:
           !IF (q1k < -2.) THEN
           !  Fng = 2.-q1k
           !ELSE IF (q1k > 0.) THEN
           !  Fng = 1.
           !ELSE
           !  Fng = 1.-1.5*q1k
           !ENDIF
           ! For purposes of the buoyancy flux in stratus, we will use Fng = 1
           Fng = 1.
       
           xl    = xl_blend(t)
           bb = b(k)*t/th(k) ! bb is "b" in BCMT95.  Their "b" differs from 
                             ! "b" in CB02 (i.e., b(k) above) by a factor 
                             ! of T/theta.  Strictly, b(k) above is formulated in
                             ! terms of sat. mixing ratio, but bb in BCMT95 is
                             ! cast in terms of sat. specific humidity.  The
                             ! conversion is neglected here. 
           qww   = 1.+0.61*qw(k)
           alpha = 0.61*th(k)
           beta  = (th(k)/t)*(xl/cp) - 1.61*th(k)
       
           vt(k) = qww   - MIN(cld(k),0.5)*beta*bb*Fng   - 1.
           vq(k) = alpha + MIN(cld(k),0.5)*beta*a(k)*Fng - tv0
           ! vt and vq correspond to beta-theta and beta-q, respectively,  
           ! in NN09, Eq. B8.  They also correspond to the bracketed
           ! expressions in BCMT95, Eq. 15, since (s*ql/sigma^2) = cldfra*Fng
           ! The "-1" and "-tv0" terms are included for consistency with 
           ! the legacy vt and vq formulations (above).

           ! increase the cloud fraction estimate below PBLH+1km
           if (zagl .lt. PBLH2+1000.) cld(k) = MIN( 1., 2.0*cld(k) ) 
           ! return a cloud condensate and cloud fraction for icloud_bl option:
           cldfra_bl1D(k) = cld(k)
           qc_bl1D(k) = ql(k)

           !To avoid FPE in radiation driver, when these two quantities are multiplied by eachother,
           ! add limit to qc_bl and cldfra_bl:
           IF (QC_BL1D(k) < 1E-6 .AND. ABS(CLDFRA_BL1D(k)) > 0.01) QC_BL1D(k)= 1E-6
           IF (CLDFRA_BL1D(k) < 1E-2)THEN
              CLDFRA_BL1D(k)=0.
              QC_BL1D(k)=0.
           ENDIF

         END DO

      END SELECT !end cloudPDF option

      !FOR TESTING PURPOSES ONLY, ISOLATE ON THE MASS-CLOUDS.
      IF (bl_mynn_cloudpdf .LT. 0) THEN
         DO k = kts,kte-1
              cldfra_bl1D(k) = 0.0
              qc_bl1D(k) = 0.0
         END DO
      ENDIF
!
      cld(kte) = cld(kte-1)
      ql(kte) = ql(kte-1)
      vt(kte) = vt(kte-1)
      vq(kte) = vq(kte-1)
      qc_bl1D(kte)=0.
      cldfra_bl1D(kte)=0.

    RETURN



  END SUBROUTINE mym_condensation

! ==================================================================
  SUBROUTINE mynn_tendencies(kts,kte,      &
       &levflag,grav_settling,             &
       &delt,dz,                           &
       &u,v,th,tk,qv,qc,qi,qni,qnc,        &
       &p,exner,                           &
       &thl,sqv,sqc,sqi,sqw,               &
       &ust,flt,flq,flqv,flqc,wspd,qcg,    &
       &uoce,voce,                         &
       &tsq,qsq,cov,                       &
       &tcd,qcd,                           &
       &dfm,dfh,dfq,                       &
       &Du,Dv,Dth,Dqv,Dqc,Dqi,Dqni,        &!Dqnc,   &
       &vdfg1,                             &
       &s_aw,s_awthl,s_awqt,s_awqv,s_awqc, &
       &s_awu,s_awv,                       &
       &sd_aw,sd_awthl,sd_awqt,sd_awqv,sd_awqc, &
       &sd_awu,sd_awv,                       &
       &FLAG_QI,FLAG_QNI,FLAG_QC,FLAG_QNC, &
       &cldfra_bl1d,                       &
       &ztop_shallow,ktop_shallow,         &
       &bl_mynn_cloudmix,                  &
       &bl_mynn_mixqt,                     &
       &bl_mynn_edmf,                      &
       &bl_mynn_edmf_dd,                   &
       &bl_mynn_edmf_mom,                  &
       &bl_mynn_cloudpdf,                  &
       &dx,PBLH,HFX,qc_bl1D,Sh,el,Vt,Vq,sgm)

!-------------------------------------------------------------------
    INTEGER, INTENT(in) :: kts,kte



    INTEGER, INTENT(in) :: grav_settling,levflag
    INTEGER, INTENT(in) :: bl_mynn_cloudmix,bl_mynn_mixqt,&
                           bl_mynn_edmf,bl_mynn_edmf_mom,bl_mynn_edmf_dd,bl_mynn_cloudpdf
    LOGICAL, INTENT(IN) :: FLAG_QI,FLAG_QNI,FLAG_QC,FLAG_QNC

!! grav_settling = 1 or 2 for gravitational settling of droplets
!! grav_settling = 0 otherwise
! thl - liquid water potential temperature
! qw - total water
! dfm,dfh,dfq - as above
! flt - surface flux of thl
! flq - surface flux of qw

    REAL,DIMENSION(kts:kte+1), INTENT(in) :: s_aw,s_awthl,s_awqt,&
                             s_awqv,s_awqc,s_awu,s_awv,&
                             sd_aw,sd_awthl,sd_awqt,&
                             sd_awqv,sd_awqc,sd_awu,sd_awv
    REAL, DIMENSION(kts:kte), INTENT(in) :: u,v,th,tk,qv,qc,qi,qni,qnc,&
         &p,exner,dfq,dz,tsq,qsq,cov,tcd,qcd
    REAL, DIMENSION(kts:kte), INTENT(inout) :: thl,sqw,sqv,sqc,sqi,&
         &dfm,dfh
    REAL, DIMENSION(kts:kte), INTENT(inout) :: du,dv,dth,dqv,dqc,dqi,&
         &dqni !,dqnc
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
    REAL, DIMENSION(kts:kte), INTENT(INOUT) :: qc_bl1D,cldfra_bl1d,sgm
    REAL, DIMENSION(kts:kte), INTENT(IN) :: Sh,el
    REAL upwind 
    LOGICAl expmf 
    
    ! upwind=1. ... use upwind approximation for mass-flux calculation
    ! upwind=0.5 ... use centered difference for mass-flux calculation
    upwind=1. 

   ! expmf=.true.  ... explicit mass-flux
   ! expmf =.false. .... implicit mass-flux
    expmf=.true. 

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
    DO k=kts,kte
       dtz(k)=delt/dz(k)
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
    ENDDO

!!============================================
!! u
!!============================================

 IF(expmf) THEN

   DO k=kts+1,kte
    upcont(k)=onoff*(s_awu(k)-s_aw(k)*(u(k+1)*upwind+u(k)*(1.-upwind)))
    dncont(k)=onoff*(sd_awu(k)-sd_aw(k)*(u(k+1)*upwind+u(k)*(1.-upwind)))
   ENDDO

    k=kts
    a(1)=0.
    b(1)=1. + dtz(k)*(dfm(k+1)+ust*ustovwsp) 
    c(1)=-dtz(k)*dfm(k+1) 
    d(1)=u(k) + dtz(k)*uoce*ust*ustovwsp -dtz(k)*(upcont(k+1)+dncont(k+1))

    DO k=kts+1,kte-1
       a(k)=   - dtz(k)*dfm(k)          
       b(k)=1. + dtz(k)*(dfm(k)+dfm(k+1)) 
       c(k)=   - dtz(k)*dfm(k+1)          
       d(k)=u(k) -dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
    ENDDO

  ELSE
   
    k=kts
    a(1)=0.
    b(1)=1. + dtz(k)*(dfm(k+1)+ust**2/wspd) - 0.5*dtz(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*sd_aw(k+1)*onoff
    c(1)=-dtz(k)*dfm(k+1) - 0.5*dtz(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*sd_aw(k+1)*onoff
    d(1)=u(k) + dtz(k)*uoce*ust**2/wspd - dtz(k)*s_awu(k+1)*onoff - dtz(k)*sd_awu(k+1)*onoff 

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

   DO k=kts+1,kte
    upcont(k)=onoff*(s_awv(k)-s_aw(k)*(v(k+1)*upwind+v(k)*(1.-upwind)))
    dncont(k)=onoff*(sd_awv(k)-sd_aw(k)*(v(k+1)*upwind+v(k)*(1.-upwind)))
   ENDDO


    k=kts
    a(1)=0.
    b(1)=1. + dtz(k)*(dfm(k+1)+ust*ustovwsp) 
    c(1)=   - dtz(k)*dfm(k+1)               
    d(1)=v(k) + dtz(k)*voce*ust*ustovwsp -dtz(k)*(upcont(k+1)+dncont(k+1))

    DO k=kts+1,kte-1
       a(k)=   - dtz(k)*dfm(k)           
       b(k)=1. + dtz(k)*(dfm(k)+dfm(k+1)) 
       c(k)=   - dtz(k)*dfm(k+1)         
       d(k)=v(k) -dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
    ENDDO

 ELSE 

    k=kts
    a(1)=0.
    b(1)=1. + dtz(k)*(dfm(k+1)+ust**2/wspd) - 0.5*dtz(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*sd_aw(k+1)*onoff
    c(1)=   - dtz(k)*dfm(k+1)               - 0.5*dtz(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*sd_aw(k+1)*onoff
    d(1)=v(k) + dtz(k)*voce*ust**2/wspd - dtz(k)*s_awv(k+1)*onoff

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
  
   DO k=kts+1,kte
   upcont(k)=s_awthl(k)-s_aw(k)*(thl(k+1)*upwind+thl(k)*(1.-upwind))
   dncont(k)=sd_awthl(k)-sd_aw(k)*(thl(k+1)*upwind+thl(k)*(1.-upwind))
   ENDDO
 
    k=kts

    a(k)=0.
    b(k)=1.+dtz(k)*dfh(k+1) 
    c(k)=  -dtz(k)*dfh(k+1) 
    d(k)=thl(k) + dtz(k)*flt + tcd(k)*delt -dtz(k)*(upcont(k+1)+dncont(k+1))

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*dfh(k)            
       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) 
       c(k)=  -dtz(k)*dfh(k+1)          
       d(k)=thl(k) + tcd(k)*delt -dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
    ENDDO

 ELSE
 
    k=kts
    a(k)=0.
    b(k)=1.+dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)
    c(k)=  -dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)
    d(k)=thl(k) + dtz(k)*flt + tcd(k)*delt -dtz(k)*s_awthl(kts+1) -dtz(k)*sd_awthl(kts+1)

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*dfh(k)            + 0.5*dtz(k)*s_aw(k) + 0.5*dtz(k)*sd_aw(k)
       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) + 0.5*dtz(k)*(s_aw(k)-s_aw(k+1)) + 0.5*dtz(k)*(sd_aw(k)-sd_aw(k+1))
       c(k)=  -dtz(k)*dfh(k+1)          - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)
       d(k)=thl(k) + tcd(k)*delt + dtz(k)*(s_awthl(k)-s_awthl(k+1)) + dtz(k)*(sd_awthl(k)-sd_awthl(k+1))
    ENDDO
 
  ENDIF



!! prescribed value
    a(nz)=0.
    b(nz)=1.
    c(nz)=0.
    d(nz)=thl(kte)

!    CALL tridiag(nz,a,b,c,d)
    CALL tridiag2(nz,a,b,c,d,x)

    DO k=kts,kte
       !thl(k)=d(k-kts+1)
       thl(k)=x(k)
    ENDDO

IF (bl_mynn_mixqt > 0) THEN
 !============================================
 ! MIX total water (sqw = sqc + sqv + sqi)
 ! NOTE: no total water tendency is output; instead, we must calculate
 !       the saturation specific humidity and then 
 !       subtract out the moisture excess (sqc & sqi)
 !============================================

 IF(expmf) THEN
   DO k=kts+1,kte
     upcont(k)=s_awqt(k)-s_aw(k)*(sqw(k+1)*upwind+sqw(k)*(1.-upwind))
     dncont(k)=sd_awqt(k)-sd_aw(k)*(sqw(k+1)*upwind+sqw(k)*(1.-upwind))
   ENDDO

    k=kts
    a(k)=0.
    b(k)=1.+dtz(k)*dfh(k+1) 
    c(k)=  -dtz(k)*dfh(k+1)

    d(k)=sqw(k) + dtz(k)*flq + qcd(k)*delt -dtz(k)*(upcont(k+1)+dncont(k+1))

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*dfh(k)           
       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) 
       c(k)=  -dtz(k)*dfh(k+1) 
       d(k)=sqw(k) + qcd(k)*delt  -dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
    ENDDO
    
  ELSE  
    k=kts

    a(k)=0.
    b(k)=1.+dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)
    c(k)=  -dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)

    !rhs= qcd(k) !+ (gfluxp - gfluxm)/dz(k)& 

    d(k)=sqw(k) + dtz(k)*flq + qcd(k)*delt - dtz(k)*s_awqt(k+1) - dtz(k)*sd_awqt(k+1)

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*dfh(k)            + 0.5*dtz(k)*s_aw(k) + 0.5*dtz(k)*sd_aw(k)
       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) + 0.5*dtz(k)*(s_aw(k)-s_aw(k+1)) + 0.5*dtz(k)*(sd_aw(k)-sd_aw(k+1))
       c(k)=  -dtz(k)*dfh(k+1)          - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)

       d(k)=sqw(k) + qcd(k)*delt + dtz(k)*(s_awqt(k)-s_awqt(k+1)) + dtz(k)*(sd_awqt(k)-sd_awqt(k+1))
    ENDDO
  ENDIF   
     
!! prescribed value
    a(nz)=0.
    b(nz)=1.
    c(nz)=0.
    d(nz)=sqw(kte)

!    CALL tridiag(nz,a,b,c,d)
    CALL tridiag2(nz,a,b,c,d,sqw2)

!    DO k=kts,kte
!       sqw2(k)=d(k-kts+1)
!    ENDDO
ELSE
    sqw2=sqw
ENDIF

IF (bl_mynn_mixqt == 0) THEN
!============================================
! cloud water ( sqc ). If mixing total water (bl_mynn_mixqt > 0),
! then sqc will be backed out of saturation check (below).
!============================================
  IF (bl_mynn_cloudmix > 0 .AND. FLAG_QC) THEN


  IF(expmf) THEN

   DO k=kts+1,kte
    upcont(k)=s_awqc(k)-s_aw(k)*(sqc(k+1)*upwind+sqc(k)*(1.-upwind))
    dncont(k)=sd_awqc(k)-sd_aw(k)*(sqc(k+1)*upwind+sqc(k)*(1.-upwind))
   ENDDO
 
    k=kts
    a(k)=0.
    b(k)=1.+dtz(k)*dfh(k+1) 
    c(k)=  -dtz(k)*dfh(k+1) 

    d(k)=sqc(k) + dtz(k)*flqc + qcd(k)*delt -dtz(k)*(upcont(k+1)+dncont(k+1))

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*dfh(k)            
       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) 
       c(k)=  -dtz(k)*dfh(k+1)         

       d(k)=sqc(k) + qcd(k)*delt -dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
    ENDDO   

  ELSE
  
     k=kts

    a(k)=0.
    b(k)=1.+dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)
    c(k)=  -dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)

    d(k)=sqc(k) + dtz(k)*flqc + qcd(k)*delt -dtz(k)*s_awqc(k+1)

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*dfh(k)            + 0.5*dtz(k)*s_aw(k) + 0.5*dtz(k)*sd_aw(k)
       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) + 0.5*dtz(k)*(s_aw(k)-s_aw(k+1)) + 0.5*dtz(k)*(sd_aw(k)-sd_aw(k+1))
       c(k)=  -dtz(k)*dfh(k+1)          - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)

       d(k)=sqc(k) + qcd(k)*delt + dtz(k)*(s_awqc(k)-s_awqc(k+1)) + dtz(k)*(sd_awqc(k)-sd_awqc(k+1))
    ENDDO
  ENDIF

! prescribed value
    a(nz)=0.
    b(nz)=1.
    c(nz)=0.
    d(nz)=sqc(kte)

!    CALL tridiag(nz,a,b,c,d)
    CALL tridiag2(nz,a,b,c,d,sqc2)

!    DO k=kts,kte
!       sqc2(k)=d(k-kts+1)
!    ENDDO
  ELSE
    !If not mixing clouds, set "updated" array equal to original array
    sqc2=sqc
  ENDIF
ENDIF

IF (bl_mynn_mixqt == 0) THEN
  !============================================
  ! MIX WATER VAPOR ONLY ( sqv ). If mixing total water (bl_mynn_mixqt > 0),
  ! then sqv will be backed out of saturation check (below).
  !============================================

 IF(expmf) THEN
   DO k=kts+1,kte
    upcont(k)=s_awqv(k)-s_aw(k)*(sqv(k+1)*upwind+sqv(k)*(1.-upwind))
    dncont(k)=sd_awqv(k)-sd_aw(k)*(sqv(k+1)*upwind+sqv(k)*(1.-upwind))
   ENDDO


    k=kts
    a(k)=0.
    b(k)=1.+dtz(k)*dfh(k+1) 
    c(k)=  -dtz(k)*dfh(k+1)
    d(k)=sqv(k) + dtz(k)*flqv + qcd(k)*delt -dtz(k)*(upcont(k+1)+dncont(k+1))

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*dfh(k)           
       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) 
       c(k)=  -dtz(k)*dfh(k+1)          
       d(k)=sqv(k) + qcd(k)*delt -dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
    ENDDO

  ELSE
    k=kts

    a(k)=0.
    b(k)=1.+dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)
    c(k)=  -dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)
    d(k)=sqv(k) + dtz(k)*flqv + qcd(k)*delt - dtz(k)*s_awqv(k+1)  !note: using qt, not qv...

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*dfh(k)            + 0.5*dtz(k)*s_aw(k) + 0.5*dtz(k)*sd_aw(k)
       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) + 0.5*dtz(k)*(s_aw(k)-s_aw(k+1)) + 0.5*dtz(k)*(sd_aw(k)-sd_aw(k+1))
       c(k)=  -dtz(k)*dfh(k+1)          - 0.5*dtz(k)*s_aw(k+1) - 0.5*dtz(k)*sd_aw(k+1)
       d(k)=sqv(k) + qcd(k)*delt + dtz(k)*(s_awqv(k)-s_awqv(k+1)) + dtz(k)*(sd_awqv(k)-sd_awqv(k+1))
    ENDDO
  ENDIF


! prescribed value
    a(nz)=0.
    b(nz)=1.
    c(nz)=0.
    d(nz)=sqv(kte)

!    CALL tridiag(nz,a,b,c,d)
    CALL tridiag2(nz,a,b,c,d,sqv2)

ELSE
    sqv2=sqv
ENDIF

!============================================
! MIX CLOUD ICE ( sqi )                      
!============================================
IF (bl_mynn_cloudmix > 0 .AND. FLAG_QI) THEN

    k=kts

    a(k)=0.
    b(k)=1.+dtz(k)*dfh(k+1)
    c(k)=  -dtz(k)*dfh(k+1)
    d(k)=sqi(k) !+ qcd(k)*delt !should we have qcd for ice?

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*dfh(k)
       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1))
       c(k)=  -dtz(k)*dfh(k+1)
       d(k)=sqi(k) !+ qcd(k)*delt
    ENDDO

!! no flux at the top
!    a(nz)=-1.       
!    b(nz)=1.        
!    c(nz)=0.        
!    d(nz)=0.        

!! specified gradient at the top
!assume gradqw_top=gradqv_top
!    a(nz)=-1.
!    b(nz)=1.
!    c(nz)=0.
!    d(nz)=gradqv_top*dztop

!! prescribed value
    a(nz)=0.
    b(nz)=1.
    c(nz)=0.
    d(nz)=sqi(kte)

!    CALL tridiag(nz,a,b,c,d)
    CALL tridiag2(nz,a,b,c,d,sqi2)

!    DO k=kts,kte
!       sqi2(k)=d(k-kts+1)
!    ENDDO
ELSE
   sqi2=sqi
ENDIF

!!============================================
!! cloud ice number concentration (qni)
!!============================================
! diasbled this since scalar_pblmix option can be invoked instead
!IF (bl_mynn_cloudmix > 0 .AND. FLAG_QNI) THEN
!
!    k=kts
!
!    a(1)=0.
!    b(1)=1.+dtz(k)*dfh(k+1)
!    c(1)=-dtz(k)*dfh(k+1)
!
!    rhs = qcd(k)
!
!    d(1)=qni(k) !+ dtz(k)*flqc + rhs*delt
!
!    DO k=kts+1,kte-1
!       kk=k-kts+1
!       a(kk)=-dtz(k)*dfh(k)
!       b(kk)=1.+dtz(k)*(dfh(k)+dfh(k+1))
!       c(kk)=-dtz(k)*dfh(k+1)
!
!       rhs = qcd(k)
!       d(kk)=qni(k) !+ rhs*delt
!
!    ENDDO
!
!! prescribed value
!    a(nz)=0.
!    b(nz)=1.
!    c(nz)=0.
!    d(nz)=qni(kte)
!
!    CALL tridiag(nz,a,b,c,d)
!
!    DO k=kts,kte
!       qni2(k)=d(k-kts+1)
!    ENDDO
!ELSE
    qni2=qni
!ENDIF

!!=============================================================
!!![EW]: call condensation again to get consistent ql as before
!!=============================================================
    IF (bl_mynn_mixqt == 2) THEN
        th_temp = th
        DO itr = 1, 2 !a few iterations to get a more accurate qc
          CALL  mym_condensation ( kts,kte,      &
               &dx,dz,thl,sqw2,p,exner,          &
               &tsq, qsq, cov,                   &
               &Sh,el,bl_mynn_cloudpdf,          &
               &qc_bl1D,cldfra_bl1D,             &
               &PBLH,HFX,                        &
               &Vt, Vq, th_temp, sgm )
          th_temp= thl + xlvcp/exner*qc_bl1D
      ENDDO
    ENDIF
!!============================================
!! Compute tendencies and convert to mixing ratios for WRF.
!! Note that the momentum tendencies are calculated above.
!!============================================

    DO k=kts,kte
       IF (bl_mynn_mixqt > 0) THEN
         t  = th(k)*exner(k)
         !SATURATED VAPOR PRESSURE
         esat=esat_blend(t)
         !SATURATED SPECIFIC HUMIDITY
         qsl=ep_2*esat/(p(k)-ep_3*esat)

         !IF (qsl >= sqw2(k)) THEN !unsaturated
         !   sqv2(k) = MAX(0.0,sqw2(k))
         !   sqi2(k) = MAX(0.0,sqi2(k))
         !   sqc2(k) = MAX(0.0,sqw2(k) - sqv2(k) - sqi2(k))
         !ELSE                     !saturated
            IF (FLAG_QI) THEN
              !sqv2(k) = qsl
              sqi2(k) = MAX(0., sqi2(k))
              sqc2(k) = MAX(0., sqw2(k) - sqi2(k) - qsl)      !updated cloud water
              sqv2(k) = MAX(0., sqw2(k) - sqc2(k) - sqi2(k))  !updated water vapor
            ELSE
              !sqv2(k) = qsl
              sqi2(k) = 0.0
              sqc2(k) = MAX(0., sqw2(k) - qsl)         !updated cloud water
              sqv2(k) = MAX(0., sqw2(k) - sqc2(k))     ! updated water vapor
            ENDIF
         !ENDIF
       ENDIF
       IF (bl_mynn_mixqt == 2) THEN
            sqi2(k) = 0.0 ! no ice for now
            sqc2(k) = qc_bl1D(k)/(1.+qc_bl1D(k)) ! updated liquid water from condensation routine
            sqv2(k) = MAX(0., sqw2(k) - sqc2(k)) ! updated water vapor
       ENDIF
       ! print *, "k = ", k, " , qc_bl = ", qc_bl1D(k)/(1.+qc_bl1D(k)), " qc orig = ", sqc2(k), " sqw2 = ", sqw2(k)
       !=====================
       ! WATER VAPOR TENDENCY
       !=====================
       Dqv(k)=(sqv2(k)/(1.-sqv2(k)) - qv(k))/delt
       !IF(-Dqv(k) > qv(k)) Dqv(k)=-qv(k)

       !=====================
       ! CLOUD WATER TENDENCY
       !=====================
       !qc fog settling tendency is now computed in module_bl_fogdes.F, so
       !sqc should only be changed by eddy diffusion or mass-flux.
       !print*,"FLAG_QC:",FLAG_QC
       ! print *, "In MYNN tendencies, k = ," , k, ", qc previous = ", qc(k)
       IF (bl_mynn_cloudmix > 0 .AND. FLAG_QC) THEN
         Dqc(k)=(sqc2(k)/(1.-sqc2(k)) - qc(k))/delt
         IF(Dqc(k)*delt + qc(k) < 0.) THEN
           ! print*,'  neg qc: ',qsl,' ',sqw2(k),' ',sqi2(k),' ',sqc2(k),' ',qc(k),' ',tk(k)
           Dqc(k)=-qc(k)/delt 
         ENDIF

         !REMOVED MIXING OF QNC - PERFORMED IN THE SCALAR_PBLMIX OPTION
         !IF (FLAG_QNC) THEN
         !  IF(sqc2(k)>1.e-9)qnc2(k)=MAX(qnc2(k),1.e6)
         !  Dqnc(k) = (qnc2(k)-qnc(k))/delt
         !  IF(Dqnc(k)*delt + qnc(k) < 0.)Dqnc(k)=-qnc(k)/delt 
         !ELSE
         !  Dqnc(k) = 0.
         !ENDIF
       ELSE
         Dqc(k)=0.
         !Dqnc(k)=0.
       ENDIF

       !===================
       ! CLOUD ICE TENDENCY
       !===================
       IF (bl_mynn_cloudmix > 0 .AND. FLAG_QI) THEN
         Dqi(k)=(sqi2(k)/(1.-sqi2(k)) - qi(k))/delt
         IF(Dqi(k)*delt + qi(k) < 0.) THEN
           !print*,' neg qi; ',qsl,' ',sqw2(k),' ',sqi2(k),' ',sqc2(k),' ',qi(k),' ',tk(k)
           Dqi(k)=-qi(k)/delt
         ENDIF

         !REMOVED MIXING OF QNI - PERFORMED IN THE SCALAR_PBLMIX OPTION
         !SET qni2 = qni above, so all tendencies are zero
         IF (FLAG_QNI) THEN
           Dqni(k)=(qni2(k)-qni(k))/delt
           IF(Dqni(k)*delt + qni(k) < 0.)Dqni(k)=-qni(k)/delt
         ELSE
           Dqni(k)=0.
         ENDIF
       ELSE
         Dqi(k)=0.
         Dqni(k)=0.
       ENDIF

       !===================
       ! THETA TENDENCY
       !===================
       IF (FLAG_QI) THEN
         Dth(k)=(thl(k) + xlvcp/exner(k)*sqc(k) &
           &            + xlscp/exner(k)*sqi(k) &
           &            - th(k))/delt
         !Use form from Tripoli and Cotton (1981) with their
         !suggested min temperature to improve accuracy:
         !Dth(k)=(thl(k)*(1.+ xlvcp/MAX(tk(k),TKmin)*sqc2(k)  &
         !  &               + xlscp/MAX(tk(k),TKmin)*sqi2(k)) &
         !  &               - th(k))/delt
       ELSE
         Dth(k)=(thl(k)+xlvcp/exner(k)*sqc2(k) - th(k))/delt
         !Use form from Tripoli and Cotton (1981) with their
         !suggested min temperature to improve accuracy.
         !Dth(k)=(thl(k)*(1.+ xlvcp/MAX(tk(k),TKmin)*sqc2(k))  &
         !&               - th(k))/delt
       ENDIF

    ENDDO



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
       &u,v,w,th,qv,qc,qi,qni,qnc,      &
       &p,exner,rho,T3D,                &
       &xland,ts,qsfc,qcg,ps,           &
       &ust,ch,hfx,qfx,rmol,wspd,       &
       &uoce,voce,                      & !ocean current
       &vdfg,                           & !Katata-added for fog dep
       &Qke,                    &
       &Tsq,Qsq,Cov,                    &
       &RUBLTEN,RVBLTEN,RTHBLTEN,       &
       &RQVBLTEN,RQCBLTEN,RQIBLTEN,     &
       &RQNIBLTEN,                      &
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
       &edmf_thl,edmf_ent,edmf_qc,      &
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
         &u,v,w,th,qv,p,exner,rho,T3D
    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), OPTIONAL, INTENT(in)::&
         &qc,qi,qni,qnc
    REAL, DIMENSION(IMS:IME,JMS:JME), INTENT(in) :: xland,ust,&
         &ch,rmol,ts,qsfc,qcg,ps,hfx,qfx, wspd,uoce,voce, vdfg,znt

    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(inout) :: &
         &Qke,Tsq,Qsq,Cov
       
         

    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(inout) :: &
         &RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,RQCBLTEN,&
         &RQIBLTEN,RQNIBLTEN,RTHRATEN !,RQNCBLTEN

    REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(out) :: &
         &exch_h,exch_m

   REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), OPTIONAL, INTENT(inout) :: &
         & edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc, &
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
    REAL, DIMENSION(KTS:KTE) :: qc_bl1D,cldfra_bl1D,&
                            qc_bl1D_old,cldfra_bl1D_old

!local vars
!  INTEGER :: ITF,JTF,KTF, IMD,JMD
    INTEGER :: i,j,k
    REAL, DIMENSION(KTS:KTE) :: thl,thvl,tl,sqv,sqc,sqi,sqw,&
         &El, Dfm, Dfh, Dfq, Tcd, Qcd, Pdk, Pdt, Pdq, Pdc, &
         &Vt, Vq, sgm

    REAL, DIMENSION(KTS:KTE) :: thetav,sh,u1,v1,w1,p1,ex1,dz1,th1,tk1,rho1,&
           & qke1,tsq1,qsq1,cov1,qv1,qi1,qc1,du1,dv1,dth1,dqv1,dqc1,dqi1, &
           & k_m1,k_h1,k_q1,qni1,dqni1,qnc1 !,dqnc1

!JOE: mass-flux variables
    REAL, DIMENSION(KTS:KTE) :: dth1mf,dqv1mf,dqc1mf,du1mf,dv1mf
    REAL, DIMENSION(KTS:KTE) :: edmf_a1,edmf_w1,edmf_qt1,edmf_thl1,&
                                edmf_ent1,edmf_qc1,&
                                edmf_a_dd1,edmf_w_dd1,edmf_qt_dd1,edmf_thl_dd1,&
                                edmf_ent_dd1,edmf_qc_dd1,&
                                edmf_debug11,edmf_debug21,edmf_debug31,edmf_debug41
    REAL,DIMENSION(KTS:KTE+1) :: s_aw1,s_awthl1,s_awqt1,&
                  s_awqv1,s_awqc1,s_awu1,s_awv1,s_awqke1,&
                  sd_aw1,sd_awthl1,sd_awqt1,&
                  sd_awqv1,sd_awqc1,sd_awu1,sd_awv1,sd_awqke1

    REAL, DIMENSION(KTS:KTE+1) :: zw
    REAL :: cpm,sqcg,flt,flq,flqv,flqc,pmz,phh,exnerg,zet,& 
              &afk,abk,ts_decay,th_sfc,ztop_shallow

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
       dqc1(kts:kte)=0.0
       dqi1(kts:kte)=0.0
       dqni1(kts:kte)=0.0
       qc_bl1D(kts:kte)=0.0
       cldfra_bl1D(kts:kte)=0.0
       qc_bl1D_old(kts:kte)=0.0
       cldfra_bl1D_old(kts:kte)=0.0
!       edmf_a1(kts:kte)=0.0
!       edmf_w1(kts:kte)=0.0
!       edmf_qc1(kts:kte)=0.0
!       edmf_a_dd1(kts:kte)=0.0
!       edmf_w_dd1(kts:kte)=0.0
!       edmf_qc_dd1(kts:kte)=0.0
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
                th1(k)=th(i,k,j)
                tk1(k)=T3D(i,k,j)
                rho1(k)=rho(i,k,j)
                sqc(k)=qc(i,k,j)/(1.+qc(i,k,j))
                sqv(k)=qv(i,k,j)/(1.+qv(i,k,j))
                thetav(k)=th(i,k,j)*(1.+0.61*sqv(k))
                IF (PRESENT(qi) .AND. FLAG_QI ) THEN
                   sqi(k)=qi(i,k,j)/(1.+qi(i,k,j))
                   sqw(k)=sqv(k)+sqc(k)+sqi(k)
                   thl(k)=th(i,k,j)- xlvcp/exner(i,k,j)*sqc(k) &
                       &           - xlscp/exner(i,k,j)*sqi(k)
                   !Use form from Tripoli and Cotton (1981) with their
                   !suggested min temperature to improve accuracy.
                   !thl(k)=th(i,k,j)*(1.- xlvcp/MAX(tk1(k),TKmin)*sqc(k) &
                   !    &               - xlscp/MAX(tk1(k),TKmin)*sqi(k))
                ELSE
                   sqi(k)=0.0
                   sqw(k)=sqv(k)+sqc(k)
                   thl(k)=th(i,k,j)-xlvcp/exner(i,k,j)*sqc(k)
                   !Use form from Tripoli and Cotton (1981) with their 
                   !suggested min temperature to improve accuracy.      
                   !thl(k)=th(i,k,j)*(1.- xlvcp/MAX(tk1(k),TKmin)*sqc(k))
                ENDIF

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
               ! if (spp_pbl==1) then
               !     rstoch_col(k)=pattern_spp_pbl(i,k,j)
               ! else
                    rstoch_col(k)=0.0
               ! endif


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
                  &Psig_bl(i,j), cldfra_bl1D,  &
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
             qc1(k)= qc(i,k,j)
            
             sqv(k)= qv(i,k,j)/(1.+qv(i,k,j))
             sqc(k)= qc(i,k,j)/(1.+qc(i,k,j))
            
             IF(icloud_bl > 0)cldfra_bl1D_old(k)=cldfra_bl(i,k,j)
             IF(icloud_bl > 0)qc_bl1D_old(k)=qc_bl(i,k,j)
          
             dqc1(k)=0.0
             dqi1(k)=0.0
             dqni1(k)=0.0
             
             IF(PRESENT(qi) .AND. FLAG_QI)THEN
                qi1(k)= qi(i,k,j)
                sqi(k)= qi(i,k,j)/(1.+qi(i,k,j))
                sqw(k)= sqv(k)+sqc(k)+sqi(k)
                thl(k)= th(i,k,j) - xlvcp/exner(i,k,j)*sqc(k) &
                     &            - xlscp/exner(i,k,j)*sqi(k)
                !Use form from Tripoli and Cotton (1981) with their
                !suggested min temperature to improve accuracy.    
                !thl(k)=th(i,k,j)*(1.- xlvcp/MAX(tk1(k),TKmin)*sqc(k) &
                !    &               - xlscp/MAX(tk1(k),TKmin)*sqi(k))
             ELSE
                qi1(k)=0.0
                sqi(k)=0.0
                sqw(k)= sqv(k)+sqc(k)
                thl(k)= th(i,k,j)-xlvcp/exner(i,k,j)*sqc(k)
                !Use form from Tripoli and Cotton (1981) with their
                !suggested min temperature to improve accuracy.    
                !thl(k)=th(i,k,j)*(1.- xlvcp/MAX(tk1(k),TKmin)*sqc(k))
             ENDIF

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
! condensation uses sh and el from the previous time step 
! sh is needed when bl_mynn_pdf=1/-1
! el is needed when bl_mynn_pdf=1/-1 or 2/-2
          CALL  mym_condensation ( kts,kte,      &
               &dx,dz1,thl,sqw,p1,ex1,           &
               &tsq1, qsq1, cov1,                &
               &Sh,el,bl_mynn_cloudpdf,          &
               &qc_bl1D,cldfra_bl1D,             &
               &PBLH(i,j),HFX(i,j),              &
               &Vt, Vq, th1, sgm )
               
               
               
          !ADD TKE source driven by cloud top cooling
          IF (bl_mynn_topdown.eq.1)then
             cloudflg=.false.
             minrad=100.
             kminrad=kpbl(i,j)
             zminrad=PBLH(i,j)
             KHtopdown(kts:kte)=0.0
             TKEprodTD(kts:kte)=0.0
             maxKHtopdown(i,j)=0.0
             !CHECK FOR STRATOCUMULUS-TOPPED BOUNDARY LAYERS
             DO kk = MAX(1,kpbl(i,j)-2),kpbl(i,j)+3
                if(sqc(kk).gt. 1.e-6 .OR. sqi(kk).gt. 1.e-6 .OR. &
                   cldfra_bl1D(kk).gt.0.5) then
                   cloudflg=.true.
                endif
                if(rthraten(i,kk,j) < minrad)then
                   minrad=rthraten(i,kk,j)
                   kminrad=kk
                   zminrad=zw(kk) + 0.5*dz1(kk)
                endif
             ENDDO
             IF (cloudflg) THEN
                zl1 = dz1(kts)
                k = MAX(kpbl(i,j)-1, kminrad-1)
                !Best estimate of height of TKE source (top of downdrafts):
                !zminrad = 0.5*pblh(i,j) + 0.5*zminrad

                templ=thl(k)*ex1(k)
                !rvls is ws at full level
                rvls=100.*6.112*EXP(17.67*(templ-273.16)/(templ-29.65))*(ep_2/p1(k+1))
                temps=templ + (sqw(k)-rvls)/(cp/xlv  +  ep_2*xlv*rvls/(rd*templ**2))
                rvls=100.*6.112*EXP(17.67*(temps-273.15)/(temps-29.65))*(ep_2/p1(k+1))
                rcldb=max(sqw(k)-rvls,0.)

                !entrainment efficiency
                dthvx     = (thl(k+2) + th1(k+2)*ep_1*sqw(k+2)) &
                          - (thl(k)   + th1(k)  *ep_1*sqw(k))
                dthvx     = max(dthvx,0.1)
                tmp1      = xlvcp * rcldb/(ex1(k)*dthvx)
                !Originally from Nichols and Turton (1986), where a2 = 60, but lowered
                !here to 8, as in Grenier and Bretherton (2001).
                ent_eff   = 0.2 + 0.2*8.*tmp1

                radsum=0.
                DO kk = MAX(1,kpbl(i,j)-3),kpbl(i,j)+3
                   radflux=rthraten(i,kk,j)*ex1(kk)         !converts theta/s to temp/s
                   radflux=radflux*cp/g*(p1(kk)-p1(kk+1)) ! converts temp/s to W/m^2
                   if (radflux < 0.0 ) radsum=abs(radflux)+radsum
                ENDDO
                radsum=MIN(radsum,60.0)

                !entrainment from PBL top thermals
                bfx0 = max(radsum/rho1(k)/cp - max(sflux,0.0),0.)
                !bfx0 = max(radsum/rho1(k)/cp,0.)
                wm3    = g/thetav(k)*bfx0*MIN(pblh(i,j),1500.) ! this is wstar3(i)
                wm2    = wm2 + wm3**h2
                bfxpbl = - ent_eff * bfx0
                dthvx  = max(thetav(k+1)-thetav(k),0.1)
                we     = max(bfxpbl/dthvx,-sqrt(wm3**h2))

                DO kk = kts,kpbl(i,j)+3
                   !Analytic vertical profile
                   zfac(kk) = min(max((1.-(zw(kk+1)-zl1)/(zminrad-zl1)),zfmin),1.)
                   zfacent(kk) = 10.*MAX((zminrad-zw(kk+1))/zminrad,0.0)*(1.-zfac(kk))**3

                   !Calculate an eddy diffusivity profile (not used at the moment)
                   wscalek2(kk) = (phifac*karman*wm3*(zfac(kk)))**h1
                   !Modify shape of KH to be similar to Lock et al (2000): use pfac = 3.0
                   KHtopdown(kk) = wscalek2(kk)*karman*(zminrad-zw(kk+1))*(1.-zfac(kk))**3 !pfac
                   KHtopdown(kk) = MAX(KHtopdown(kk),0.0)
                   !Do not include xkzm at kpbl-1 since it changes entrainment
                   !if (kk.eq.kpbl(i,j)-1 .and. cloudflg .and. we.lt.0.0) then
                   !   KHtopdown(kk) = 0.0
                   !endif
                   
                   !Calculate TKE production = 2(g/TH)(w'TH'), where w'TH' = A(TH/g)wstar^3/PBLH,
                   !A = ent_eff, and wstar is associated with the radiative cooling at top of PBL.
                   !An analytic profile controls the magnitude of this TKE prod in the vertical. 
                   TKEprodTD(kk)=2.*ent_eff*wm3/MAX(pblh(i,j),100.)*zfacent(kk)
                   TKEprodTD(kk)= MAX(TKEprodTD(kk),0.0)
                ENDDO
             ENDIF !end cloud check
             maxKHtopdown(i,j)=MAXVAL(KHtopdown(:))
          ELSE
             maxKHtopdown(i,j)=0.0
             KHtopdown(kts:kte) = 0.0
             TKEprodTD(kts:kte)=0.0
          ENDIF 
          !end top-down check


          IF (bl_mynn_edmf > 0) THEN
            CALL edmf_JPL(                            &
               &kts,kte,delt,zw,p1,                   &
               &u1,v1,th1,thl,thetav,tk1,sqw,sqv,sqc, &
               &ex1,                              &
               &ust(i,j),flt,flq,PBLH(i,j),       &
               & edmf_a1,edmf_w1,edmf_qt1,        &
               & edmf_thl1,edmf_ent1,edmf_qc1,    &
               & edmf_debug11,edmf_debug21,       &
               & edmf_debug31,edmf_debug41,       &
               & s_aw1,s_awthl1,s_awqt1,          &
               & s_awqv1,s_awqc1,s_awu1,s_awv1,   &
               & s_awqke1,                        &
               & qc_bl1D,cldfra_bl1D,             &
               &FLAG_QI,FLAG_QC,                  & 
               &Psig_shcu(i,j),                   &
               &ktop_shallow(i,j),ztop_shallow,   &
               KPBL(i,j)&
            )

          ENDIF

          IF (bl_mynn_edmf_dd == 1) THEN
            ! DO k = kts, KtE
            !     print *, "[EWDEBUG] k = ",k, " sqc = ", sqc(k), " cldfra = ", cldfra_bl1d(k)
            ! ENDDO
            CALL DDMF_JPL(kts,kte,delt,zw,p1,        &
              &u1,v1,th1,thl,thetav,tk1,sqw,sqv,sqc, &
              &ex1,                                  &
              &ust(i,j),flt,flq,PBLH(i,j),KPBL(i,j), &
              &edmf_a_dd1,edmf_w_dd1, edmf_qt_dd1,   &
              &edmf_thl_dd1,edmf_ent_dd1,edmf_qc_dd1,&
              & edmf_debug31,edmf_debug41,           &
              &sd_aw1,sd_awthl1,sd_awqt1,            &
              &sd_awqv1,sd_awqc1,sd_awu1,sd_awv1,    &
              &sd_awqke1,                            &
              &qc_bl1d,cldfra_bl1d,                  &
              &rthraten(i,:,j)&
              )
          ENDIF

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
               &cldfra_bl1D,bl_mynn_mixlength,   &
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

          !!!! [EWDD] Smooth ED here !!!!!!!!!!!!
          ! DO k=kts+1, KPBL(i,j)-1
          !   dfm(k) = 0.25 * dfm(k-1) + 0.5 * dfm(k) + 0.25 * dfm(k+1)
          !   dfh(k) = 0.25 * dfm(k-1) + 0.5 * dfh(k) + 0.25 * dfh(k+1)
          !   dfq(k) = 0.25 * dfm(k-1) + 0.5 * dfq(k) + 0.25 * dfq(k+1)
          ! ENDDO    
          !!!! [EWDD] End smoothing !!!!!!!!!!!!
              
          CALL mym_predict (kts,kte,levflag,     &
               &delt, dz1,                       &
               &ust(i,j), flt, flq, pmz, phh,    &
               &el, dfq, pdk, pdt, pdq, pdc,     &
               &Qke1, Tsq1, Qsq1, Cov1,          &
               &s_aw1, s_awqke1, sd_aw1, sd_awqke1,&
               bl_mynn_edmf_tke)

          CALL mynn_tendencies(kts,kte,          &
               &levflag,grav_settling,           &
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
               &Du1, Dv1, Dth1, Dqv1,            &
               &Dqc1, Dqi1, Dqni1,               & !Dqnc1,        &
               &vdfg(i,j),                       & !JOE/Katata- fog deposition
               ! mass flux components
               &s_aw1,s_awthl1,s_awqt1,          &
               &s_awqv1,s_awqc1,s_awu1,s_awv1,   &
               &sd_aw1,sd_awthl1,sd_awqt1,          &
               &sd_awqv1,sd_awqc1,sd_awu1,sd_awv1,   &
               &FLAG_QI,FLAG_QNI,FLAG_QC,FLAG_QNC,&
               &cldfra_bl1d,                     &
               &ztop_shallow,ktop_shallow(i,j),  &
               &bl_mynn_cloudmix,                &
               &bl_mynn_mixqt,                   &
               &bl_mynn_edmf,                    &
               &bl_mynn_edmf_dd,                 &
               &bl_mynn_edmf_mom,                &
               &bl_mynn_cloudpdf,                &
               &dx,PBLH(i,j),HFX(i,j),qc_bl1D,Sh,el,Vt,Vq,sgm)


 
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
             
             IF(bl_mynn_cloudmix > 0)THEN
               IF (PRESENT(qc) .AND. FLAG_QC) RQCBLTEN(i,k,j)=dqc1(k)
               IF (PRESENT(qi) .AND. FLAG_QI) RQIBLTEN(i,k,j)=dqi1(k)
               !IF (PRESENT(qnc) .AND. FLAG_QNC) RQNCBLTEN(i,k,j)=dqnc1(k)
               IF (PRESENT(qni) .AND. FLAG_QNI) RQNIBLTEN(i,k,j)=dqni1(k)
             ELSE
               IF (PRESENT(qc) .AND. FLAG_QC) RQCBLTEN(i,k,j)=0.
               IF (PRESENT(qi) .AND. FLAG_QI) RQIBLTEN(i,k,j)=0.
               !IF (PRESENT(qnc) .AND. FLAG_QNC) RQNCBLTEN(i,k,j)=0.
               IF (PRESENT(qni) .AND. FLAG_QNI) RQNIBLTEN(i,k,j)=0.
             ENDIF

             IF(icloud_bl > 0)THEN
               !make BL clouds scale aware - may already be done in mym_condensation
               qc_bl(i,k,j)=qc_bl1D(k) !*Psig_shcu(i,j)
               cldfra_bl(i,k,j)=cldfra_bl1D(k) !*Psig_shcu(i,j)

               !Stochastic perturbations to cldfra_bl and qc_bl
              ! if (spp_pbl==1) then
              !    cldfra_bl(i,k,j)= cldfra_bl(i,k,j)*(1.0-rstoch_col(k))
              !    IF ((cldfra_bl(i,k,j) > 1.0) .or. (cldfra_bl(i,k,j) < 0.0)) then
              !       cldfra_bl(i,k,j)=MAX(MIN(cldfra_bl(i,k,j),1.0),0.0)
              !    ENDIF
              ! ELSE
                  !DIAGNOSTIC-DECAY FOR SUBGRID-SCALE CLOUDS
                  IF (CLDFRA_BL(i,k,j) > cldfra_bl1D_old(k)) THEN
                     !KEEP UPDATED CLOUD FRACTION & MIXING RATIO
                  ELSE
                     !DECAY TIMESCALE FOR CALM CONDITION IS THE EDDY TURNOVER TIMESCALE, 
                     !BUT FOR WINDY CONDITIONS, IT IS THE ADVECTIVE TIMESCALE.
                     !USE THE MINIMUM OF THE TWO.
                     ts_decay = MIN( 1800., 3.*dx/MAX(SQRT(u1(k)**2 + v1(k)**2), 1.0) )
                     cldfra_bl(i,k,j)= MAX(cldfra_bl1D(k), cldfra_bl1D_old(k)-(0.25*delt/ts_decay))
                     IF (cldfra_bl(i,k,j) > 0.01) THEN
                        IF (QC_BL(i,k,j) < 1E-5)QC_BL(i,k,j)= MAX(qc_bl1D_old(k), 1E-5)
                     ELSE
                        CLDFRA_BL(i,k,j)= 0.
                        QC_BL(i,k,j)    = 0.
                     ENDIF
                  ENDIF
               !ENDIF

               !Reapply checks on cldfra_bl and qc_bl to avoid FPEs in radiation driver
               ! when these two quantities are multiplied by eachother (they may have changed
               ! in the MF scheme:
               IF (icloud_bl > 0) THEN
                 IF (QC_BL(i,k,j) < 1E-6 .AND. ABS(CLDFRA_BL(i,k,j)) > 0.1)QC_BL(i,k,j)= 1E-6
                 IF (CLDFRA_BL(i,k,j) < 1E-2)CLDFRA_BL(i,k,j)= 0.
               ENDIF
             ENDIF
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
               if (bl_mynn_edmf_dd > 0) THEN
                   !update downdraft properties
                   edmf_a_dd(i,k,j)=edmf_a_dd1(k)
                   edmf_w_dd(i,k,j)=edmf_w_dd1(k)
                   edmf_qt_dd(i,k,j)=edmf_qt_dd1(k)
                   edmf_thl_dd(i,k,j)=edmf_thl_dd1(k)
                   edmf_ent_dd(i,k,j)=edmf_ent_dd1(k)
                   edmf_qc_dd(i,k,j)=edmf_qc_dd1(k)
               ENDIF

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



print *,'qkeEND',qke

  END SUBROUTINE mynn_bl_driver

! ==================================================================
  SUBROUTINE mynn_bl_init_driver(                   &
       &RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,          &
       &RQCBLTEN,RQIBLTEN & !,RQNIBLTEN,RQNCBLTEN       &
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
         &RQCBLTEN,RQIBLTEN,& !RQNIBLTEN,RQNCBLTEN       &
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
                if( p_qc >= param_first_scalar ) RQCBLTEN(i,k,j)=0.
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
  
!=================================================================
subroutine Poisson(istart,iend,jstart,jend,mu,POI)

  integer, intent(in) :: istart,iend,jstart,jend 
  real,dimension(istart:iend,jstart:jend),intent(in) :: MU
  integer, dimension(istart:iend,jstart:jend), intent(out) :: POI
  integer :: i,j
  !
  ! do this only once
  ! call init_random_seed

  do i=istart,iend
    do j=jstart,jend
       !print *, "Poisson mu: ", mu(i,j) 
      call   random_Poisson(mu(i,j),.true.,POI(i,j))
       !print *, "POI: ", POI(i,j)
    enddo
  enddo

end subroutine Poisson
!=================================================================  
subroutine init_random_seed()
   !JOE: PGI had problem! use iso_fortran_env, only: int64
   !JOE: PGI had problem! use ifport, only: getpid 
   implicit none
   integer, allocatable :: seed(:)
   integer :: i, n, un, istat, dt(8), pid
   !JOE: PGI had problem! integer(int64) :: t
   integer :: t

   call random_seed(size = n)
   allocate(seed(n))

   ! First try if the OS provides a random number generator
   !JOE: PGI had problem! open(newunit=un, file="/dev/urandom", access="stream", &
   un=191
   open(unit=un, file="/dev/urandom", access="stream", &
   form="unformatted", action="read", status="old", iostat=istat)

   if (istat == 0) then
      read(un) seed
      close(un)
   else
      ! Fallback to XOR:ing the current time and pid. The PID is
      ! useful in case one launches multiple instances of the same
      ! program in parallel.
      call system_clock(t)
      if (t == 0) then
         call date_and_time(values=dt)
         !t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
         !   + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
         !   + dt(3) * 24_int64 * 60 * 60 * 1000 &
         !   + dt(5) * 60 * 60 * 1000 &
         !   + dt(6) * 60 * 1000 + dt(7) * 1000 &
         !   + dt(8)
         t = dt(6) * 60 &  ! only return seconds for smaller t
           + dt(7)
      end if

      !JOE: PGI had problem!pid = getpid()
      ! for distributed memory jobs we need to fix this
      !pid=1
      pid = 666 + MOD(t,10)  !JOE: doesnt work for PG compilers: getpid()
 
      t = ieor(t, int(pid, kind(t)))
      do i = 1, n
         seed(i) = lcg(t)
      end do
   end if
   call random_seed(put=seed)

  contains

  ! Pseudo-random number generator (PRNG) 
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)

   integer :: lcg
   !JOE: PGI had problem! integer(int64) :: s
   integer :: s

   if (s == 0) then
      !s = 104729
      s = 1047
   else
      !s = mod(s, 4294967296_int64)
      s = mod(s, 71)
   end if
   !s = mod(s * 279470273_int64, 4294967291_int64)
   s = mod(s * 23, 17)
   !lcg = int(mod(s, int(huge(0), int64)), kind(0))
   lcg = int(mod(s, int(s/3.5)))

  end function lcg

  end subroutine init_random_seed

  SUBROUTINE init_random_seed_1()
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
  
    CALL SYSTEM_CLOCK(COUNT=clock)
  
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
  
    DEALLOCATE(seed)
  END SUBROUTINE init_random_seed_1

subroutine random_Poisson(mu,first,ival) 
!**********************************************************************
!     Translated to Fortran 90 by Alan Miller from:  RANLIB
!
!     Library of Fortran Routines for Random Number Generation
!
!                    Compiled and Written by:
!
!                         Barry W. Brown
!                          James Lovato
!
!             Department of Biomathematics, Box 237
!             The University of Texas, M.D. Anderson Cancer Center
!             1515 Holcombe Boulevard
!             Houston, TX      77030
!
! Generates a single random deviate from a Poisson distribution with mean mu.
! Scalar Arguments:
    REAL, INTENT(IN)    :: mu  !The mean of the Poisson distribution from which
                                   !a random deviate is to be generated.
    LOGICAL, INTENT(IN) :: first
        INTEGER             :: ival

!     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
!     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
!     SEPARATION OF CASES A AND B
!
!     .. Local Scalars ..
!JOE: since many of these scalars conflict with globally declared closure constants (above),
!     need to change XX to XX_s
!   REAL          :: b1, b2, c, c0, c1, c2, c3, del, difmuk, e, fk, fx, fy, g,  &
!                    omega, px, py, t, u, v, x, xx
    REAL          :: b1_s, b2_s, c, c0, c1_s, c2_s, c3_s, del, difmuk, e, fk, fx, fy, g_s,  &
                    omega, px, py, t, u, v, x, xx
    REAL, SAVE    :: s, d, p, q, p0
        INTEGER       :: j, k, kflag
    LOGICAL, SAVE :: full_init
        INTEGER, SAVE :: l, m
!     ..
!     .. Local Arrays ..
    REAL, SAVE    :: pp(35)
!     ..
!     .. Data statements ..
!JOE: since many of these scalars conflict with globally declared closure constants (above),
!     need to change XX to XX_s
!   REAL, PARAMETER :: a0 = -.5, a1 = .3333333, a2 = -.2500068, a3 = .2000118,  &
    REAL, PARAMETER :: a0 = -.5, a1_s = .3333333, a2_s = -.2500068, a3 = .2000118,  &
                           a4 = -.1661269, a5 = .1421878, a6 = -0.1384794,  &
                           a7 = .1250060

    REAL, PARAMETER :: fact(10) = (/ 1., 1., 2., 6., 24., 120., 720., 5040.,  &
                 40320., 362880. /)

!JOE: difmuk,fk,u errors - undefined
   difmuk = 0.
   fk = 1.0
   u = 0.

!     ..
!     .. Executable Statements ..
   IF (mu > 10.0) THEN
!     C A S E  A. (RECALCULATION OF S, D, L IF MU HAS CHANGED)

      IF (first) THEN
         s = SQRT(mu)
         d = 6.0*mu*mu

!             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
!             PROBABILITIES FK WHENEVER K >= M(MU). L=IFIX(MU-1.1484)
!             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .

         l = mu - 1.1484
         full_init = .false.
      END IF

!     STEP N. NORMAL SAMPLE - random_normal() FOR STANDARD NORMAL DEVIATE
      g_s = mu + s*random_normal()
      IF (g_s > 0.0) THEN
         ival = g_s

     !     STEP I. IMMEDIATE ACCEPTANCE IF ival IS LARGE ENOUGH
         IF (ival>=l) RETURN

     !     STEP S. SQUEEZE ACCEPTANCE - SAMPLE U
        fk = ival
        difmuk = mu - fk
        CALL RANDOM_NUMBER(u)
        IF (d*u >= difmuk*difmuk*difmuk) RETURN
      END IF

      !     STEP P. PREPARATIONS FOR STEPS Q AND H.
      !             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
      !             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
      !             THE QUANTITIES B1_S, B2_S, C3_S, C2_S, C1_S, C0 ARE FOR THE HERMITE
      !             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
      !             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.

      IF (.NOT. full_init) THEN
         omega = .3989423/s
        b1_s = .4166667E-1/mu
        b2_s = .3*b1_s*b1_s
        c3_s = .1428571*b1_s*b2_s
        c2_s = b2_s - 15.*c3_s
        c1_s = b1_s - 6.*b2_s + 45.*c3_s
        c0 = 1. - b1_s + 3.*b2_s - 15.*c3_s
        c = .1069/mu
        full_init = .true.
      END IF

      IF (g_s < 0.0) GO TO 50

    !             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)

      kflag = 0
      GO TO 70

    !     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)

      40 IF (fy-u*fy <= py*EXP(px-fx)) RETURN

    !     STEP E. EXPONENTIAL SAMPLE - random_exponential() FOR STANDARD EXPONENTIAL
    !             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
    !             (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)

      50 e = random_exponential()
      CALL RANDOM_NUMBER(u)
      u = u + u - one
      t = 1.8 + SIGN(e, u)
      IF (t <= (-.6744)) GO TO 50
      ival = mu + s*t
      fk = ival
      difmuk = mu - fk

    !             'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)

      kflag = 1
      GO TO 70

    !     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)

      60 IF (c*ABS(u) > py*EXP(px+e) - fy*EXP(fx+e)) GO TO 50
      RETURN

    !     STEP F. 'SUBROUTINE' F. CALCULATION OF PX, PY, FX, FY.
    !             CASE ival < 10 USES FACTORIALS FROM TABLE FACT

      70 IF (ival>=10) GO TO 80
      px = -mu
!JOE: had error " Subscript #1 of FACT has value -858993459"; shouldn't be < 1.
         !py = mu**ival/fact(ival+1)
      py = mu**ival/fact(MAX(ival+1,1))
      GO TO 110

    !             CASE ival >= 10 USES POLYNOMIAL APPROXIMATION
    !             A0-A7 FOR ACCURACY WHEN ADVISABLE
    !             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)

      80 del = .8333333E-1/fk
      del = del - 4.8*del*del*del
      v = difmuk/fk
      IF (ABS(v)>0.25) THEN
        px = fk*LOG(one + v) - difmuk - del
      ELSE
        px = fk*v*v* (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2_s)*v+a1_s)*v+a0) - del
      END IF
      py = .3989423/SQRT(fk)
      110 x = (half - difmuk)/s
      xx = x*x
      fx = -half*xx
      fy = omega* (((c3_s*xx + c2_s)*xx + c1_s)*xx + c0)
      IF (kflag <= 0) GO TO 40
      GO TO 60

    !---------------------------------------------------------------------------
    !     C A S E  B.    mu < 10
    !     START NEW TABLE AND CALCULATE P0 IF NECESSARY
      ELSE

      IF (first) THEN
        m = MAX(1, INT(mu))
        l = 0
                !print*,"mu=",mu
                !print*," mu=",mu," p=",EXP(-mu)
        p = EXP(-mu)
        q = p
        p0 = p
      END IF

    !     STEP U. UNIFORM SAMPLE FOR INVERSION METHOD

      DO
        CALL RANDOM_NUMBER(u)
        ival = 0
        IF (u <= p0) RETURN

    !     STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
    !             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
    !             (0.458=PP(9) FOR MU=10)

        IF (l == 0) GO TO 150
        j = 1
        IF (u > 0.458) j = MIN(l, m)
        DO k = j, l
          IF (u <= pp(k)) GO TO 180
        END DO
        IF (l == 35) CYCLE

    !     STEP C. CREATION OF NEW POISSON PROBABILITIES P
    !             AND THEIR CUMULATIVES Q=PP(K)

        150 l = l + 1
        DO k = l, 35
          p = p*mu / k
          q = q + p
          pp(k) = q
          IF (u <= q) GO TO 170
        END DO
        l = 35
      END DO

      170 l = k
      180 ival = k
      RETURN
    END IF

    RETURN
    END subroutine random_Poisson

!==================================================================

    FUNCTION random_normal() RESULT(fn_val)

    ! Adapted from the following Fortran 77 code
    !      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
    !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
    !      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

    !  The function random_normal() returns a normally distributed pseudo-random
    !  number with zero mean and unit variance.

    !  The algorithm uses the ratio of uniforms method of A.J. Kinderman
    !  and J.F. Monahan augmented with quadratic bounding curves.

    REAL :: fn_val

    !     Local variables
    REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
                r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

    !     Generate P = (u,v) uniform in rectangle enclosing acceptance region

    DO
      CALL RANDOM_NUMBER(u)
      CALL RANDOM_NUMBER(v)
      v = 1.7156 * (v - half)

    !     Evaluate the quadratic form
      x = u - s
      y = ABS(v) - t
      q = x**2 + y*(a*y - b*x)

    !     Accept P if inside inner ellipse
      IF (q < r1) EXIT
    !     Reject P if outside outer ellipse
      IF (q > r2) CYCLE
    !     Reject P if outside acceptance region
      IF (v**2 < -4.0*LOG(u)*u**2) EXIT
    END DO

    !     Return ratio of P coordinates as the normal deviate
    fn_val = v/u
    RETURN

    END FUNCTION random_normal

!===============================================================

    FUNCTION random_exponential() RESULT(fn_val)

    ! Adapted from Fortran 77 code from the book:
    !     Dagpunar, J. 'Principles of random variate generation'
    !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

    ! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
    ! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
    ! TO EXP(-random_exponential), USING INVERSION.

    REAL  :: fn_val

    !     Local variable
    REAL  :: r

    DO
      CALL RANDOM_NUMBER(r)
      IF (r > zero) EXIT
    END DO

    fn_val = -LOG(r)
    RETURN

    END FUNCTION random_exponential

!===============================================================

subroutine condensation_edmf(QT,THL,P,zagl,THV,QC)
!
! zero or one condensation for edmf: calculates THV and QC
!
real,intent(in)   :: QT,THL,P,zagl
real,intent(out)  :: THV
real,intent(inout):: QC

integer :: niter,i
real :: diff,exn,t,th,qs,qcold

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
  !QC=0.  !better first guess QC is incoming from lower level, do not set to zero
  do i=1,NITER
     T=EXN*THL + xlv/cp*QC        
     QS=qsat_blend(T,P)
     QCOLD=QC
     QC=0.5*QC + 0.5*MAX((QT-QS),0.)
     if (abs(QC-QCOLD)<Diff) exit
  enddo

  T=EXN*THL + xlv/cp*QC
  QS=qsat_blend(T,P)
  QC=max(QT-QS,0.)

  !Do not allow saturation below 100 m
  !if(zagl < 100.)QC=0.

  !THV=(THL+xlv/cp*QC).*(1+(1-rvovrd)*(QT-QC)-QC);
  THV=(THL+xlv/cp*QC)*(1.+QT*(rvovrd-1.)-rvovrd*QC)
  !THIS BASICALLY GIVE THE SAME RESULT AS THE PREVIOUS LINE
  !TH = THL + xlv/cp/EXN*QC
  !THV= TH*(1. + 0.608*QT)

  !print *,'t,p,qt,qs,qc'
  !print *,t,p,qt,qs,qc 


end subroutine condensation_edmf

!===============================================================

subroutine condensation_edmf_r(QT,THL,P,zagl,THV,QC)
!
! zero or one condensation for edmf: calculates THL and QC
! similar to condensation_edmf but with different inputs
!
real,intent(in)   :: QT,THV,P,zagl
real,intent(out)  :: THL, QC

integer :: niter,i
real :: diff,exn,t,th,qs,qcold

! number of iterations
  niter=50
! minimum difference
  diff=2.e-5

  EXN=(P/p1000mb)**rcp
  ! assume first that th = thv
  T = THV*EXN
  !QS = qsat_blend(T,P)
  !QC = QS - QT
  
  QC=0.

  do i=1,NITER    
     QCOLD = QC
     T = EXN*THV/(1.+QT*(rvovrd-1.)-rvovrd*QC)   
     QS=qsat_blend(T,P)
     QC= MAX((QT-QS),0.)
     if (abs(QC-QCOLD)<Diff) exit
  enddo
  THL = (T - xlv/cp*QC)/EXN

end subroutine condensation_edmf_r

!===============================================================
function qs_edmf(t,p)
!
! calculates saturated water pressure
! at temperature t and pressure p
real, intent(in) :: t,p
real qsl,esl
real qs_edmf  
       esl=svp11*EXP(svp2*(t-svpt0)/(t-svp3))
       !SATURATED SPECIFIC HUMIDITY
       qs_edmf=ep_2*esl/(p-ep_3*esl)
 !  print *,'esl,qsl,svp11,svp2,svpt0,svp3,ep_2esl'  
!   print *,esl,qsl,svp11,svp2,svpt0,svp3,ep_2
    

end function qs_edmf

!===============================================================

! =====================================================================

  FUNCTION esat_blend(t) 
! JAYMES- added 22 Apr 2015
! 
! This calculates saturation vapor pressure.  Separate ice and liquid functions 
! are used (identical to those in module_mp_thompson.F, v3.6).  Then, the 
! final returned value is a temperature-dependant "blend".  Because the final 
! value is "phase-aware", this formulation may be preferred for use throughout 
! the module (replacing "svp").

      IMPLICIT NONE
      
      REAL, INTENT(IN):: t
      REAL :: esat_blend,XC,ESL,ESI,chi

      XC=MAX(-80.,t-273.16)

! For 253 < t < 273.16 K, the vapor pressures are "blended" as a function of temperature, 
! using the approach of Chaboureau and Bechtold (2002), JAS, p. 2363.  The resulting 
! values are returned from the function.
      IF (t .GE. 273.16) THEN
          esat_blend = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8))))))) 
      ELSE IF (t .LE. 253.) THEN
          esat_blend = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
      ELSE
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          chi  = (273.16-t)/20.16
          esat_blend = (1.-chi)*ESL  + chi*ESI
      END IF

  END FUNCTION esat_blend

! ====================================================================

  FUNCTION qsat_blend(t, P)
! JAYMES- this function extends function "esat" and returns a "blended"
! saturation mixing ratio.

      IMPLICIT NONE

      REAL, INTENT(IN):: t, P
      REAL :: qsat_blend,XC,ESL,ESI,RSLF,RSIF,chi

      XC=MAX(-80.,t-273.16)

      IF (t .GE. 273.16) THEN
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8))))))) 
          qsat_blend = 0.622*ESL/(P-ESL) 
      ELSE IF (t .LE. 253.) THEN
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          qsat_blend = 0.622*ESI/(P-ESI)
      ELSE
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          RSLF = 0.622*ESL/(P-ESL)
          RSIF = 0.622*ESI/(P-ESI)
          chi  = (273.16-t)/20.16
          qsat_blend = (1.-chi)*RSLF + chi*RSIF
      END IF

  END FUNCTION qsat_blend

! ===================================================================

  FUNCTION xl_blend(t)
! JAYMES- this function interpolates the latent heats of vaporization and
! sublimation into a single, temperature-dependant, "blended" value, following
! Chaboureau and Bechtold (2002), Appendix.

      IMPLICIT NONE

      REAL, INTENT(IN):: t
      REAL :: xl_blend,xlvt,xlst,chi

      IF (t .GE. 273.16) THEN
          xl_blend = xlv + (cpv-cliq)*(t-273.16)  !vaporization/condensation
      ELSE IF (t .LE. 253.) THEN
          xl_blend = xls + (cpv-cice)*(t-273.16)  !sublimation/deposition
      ELSE
          xlvt = xlv + (cpv-cliq)*(t-273.16)  !vaporization/condensation
          xlst = xls + (cpv-cice)*(t-273.16)  !sublimation/deposition
          chi  = (273.16-t)/20.16
          xl_blend = (1.-chi)*xlvt + chi*xlst     !blended
      END IF

  END FUNCTION xl_blend

! ===================================================================
! ===================================================================
! This is the original mass flux scheme from Kay Suselj of NASA-JPL,


SUBROUTINE edmf_JPL(kts,kte,dt,zw,p,         &
              &u,v,th,thl,thv,tk,qt,qv,qc,   &
              &exner,                        &
              &ust,wthl,wqt,pblh,            &
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
            ! inputs - flags for moist arrays
              &F_QC,F_QI,                    &
              &Psig_shcu,                    &
            ! output info
              &ktop,ztop,kpbl)

        INTEGER, INTENT(IN) :: KTS,KTE, kpbl
        REAL,DIMENSION(KTS:KTE), INTENT(IN) :: U,V,TH,THL,TK,QT,QV,QC,THV,P,exner
        ! zw .. heights of the updraft levels (edges of boxes)
        REAL,DIMENSION(KTS:KTE+1), INTENT(IN) :: ZW
        REAL, INTENT(IN) :: DT,UST,WTHL,WQT,PBLH,Psig_shcu

        LOGICAL, OPTIONAL :: F_QC,F_QI
        LOGICAL :: cloudflg
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

        INTEGER, PARAMETER :: NUP=100, debug_mf=0 !fixing number of plumes to 10

    ! updraft properties
        REAL,DIMENSION(KTS:KTE+1,1:NUP) :: UPW,UPTHL,UPQT,UPQC,UPA,UPU,UPV,UPTHV

    ! entrainment variables
        REAl,DIMENSION(KTS+1:KTE+1,1:NUP) :: ENT,ENTf
        REAL,DIMENSION(KTS+1:KTE+1) :: L0s
        INTEGER,DIMENSION(KTS+1:KTE+1,1:NUP) :: ENTi

    ! internal variables
        INTEGER :: K,I, qlTop
        REAL :: wthv,wstar,qstar,thstar,sigmaW,sigmaQT,sigmaTH,z0, &
            pwmin,pwmax,wmin,wmax,wlv,wtv,dthv_dz
        REAL :: B,QTn,THLn,THVn,QCn,Un,Vn,Wn2,EntEXP,EntW, deltaZ,EntExp_M, Beta_un, Z00, Z_i

    ! VARIABLES FOR CHABOUREAU-BECHTOLD CLOUD FRACTION
        ! REAL,DIMENSION(KTS:KTE), INTENT(INOUT) :: vt, vq, sgm
        REAL :: sigq,xl,tlk,qsat_tl,rsl,cpm,a,qmq,mf_cf,diffqt,&
               Fng,qww,alpha,beta,bb,f,pt,t,q2p,b9,satvp,rhgrid

    ! w parameters
        REAL,PARAMETER :: &
            &Wa=1., &
            &Wb=1.5

    ! entrainment parameters
        REAL,PARAMETER :: &
        & L0=100.,&
        & ENT0=0.2



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

! This part is specific for Stratocumulus
cloudflg = .false.
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

 ! set initial conditions for updrafts
    z0=50.
    pwmin=1.
    pwmax=3.

    wstar=max(0.,(g/thv(1)*wthv*pblh)**(1./3.))
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
        do k=kts+1,kte
           ENTf(k,i)=(ZW(k)-ZW(k-1))/L0! s(k)
! Elynn's modification of entrainment rate
!          if (ZW(k) < Z_i) then
!            ENTf(k,i) = (ZW(k)-ZW(k-1)) * (3.5 * exp((-(Z_i-ZW(k)))/70.) + 1./L0)
!          else
!            ENTf(k,i) = (ZW(k)-ZW(k-1)) * (3.5 * exp(1./70.) + 1./L0)
!          endif
         enddo
    enddo

! get Poisson P(dz/L0)
    call Poisson(1,Nup,kts+1,kte,ENTf,ENTi)


 ! entrainent: Ent=Ent0/dz*P(dz/L0)
    do i=1,Nup
        do k=kts+1,kte
          ENT(k,i)=real(ENTi(k,i))*Ent0/(ZW(k)-ZW(k-1))
        enddo
    enddo

!    print *,'Entrainment',ENT


    DO I=1,NUP
        wlv=wmin+(wmax-wmin)/NUP*(i-1)
        wtv=wmin+(wmax-wmin)/NUP*i

        UPW(1,I)=0.5*(wlv+wtv)
        UPA(1,I)=0.5*ERF(wtv/(sqrt(2.)*sigmaW))-0.5*ERF(wlv/(sqrt(2.)*sigmaW))

        UPU(1,I)=U(1)
        UPV(1,I)=V(1)

        UPQC(1,I)=0
        UPQT(1,I)=QT(1)+0.58*UPW(1,I)*sigmaQT/sigmaW ! was 0.32
        UPTHV(1,I)=THV(1)+0.58*UPW(1,I)*sigmaTH/sigmaW
        UPTHL(1,I)=UPTHV(1,I)/(1+svp1*UPQT(1,I))
    ENDDO


  ! do integration  updraft
    DO I=1,NUP
        DO k=KTS+1,KTE
            deltaZ = ZW(k)-ZW(k-1)
            EntExp=exp(-ENT(K,I)*deltaZ)
            EntExp_M=exp(-ENT(K,I)/3.*deltaZ)

            QTn=QT(K-1)*(1-EntExp)+UPQT(K-1,I)*EntExp
            THLn=THL(K-1)*(1-EntExp)+UPTHL(K-1,I)*EntExp
            Un=U(K-1)*(1-EntExp)+UPU(K-1,I)*EntExp_M
            Vn=V(K-1)*(1-EntExp)+UPV(K-1,I)*EntExp_M
            ! get thvn,qcn
            call condensation_edmf(QTn,THLn,(P(K)+P(K-1))/2.,ZW(k),THVn,QCn)

            B=g*(0.5*(THVn+UPTHV(k-1,I))/THV(k-1)-1)
 !           IF (ZW(k) <= Z_i) THEN
 !               Beta_un = Wb*ENT(K,I) + 0.5/(Z_i-ZW(k)+deltaZ) * max(1. - exp((Z_i-ZW(k)+deltaZ)/Z00-1.) , 0.)
 !           ELSE
 !               Beta_un = Wb*ENT(K,I) + 0.5/deltaZ * (1. - exp(deltaZ/Z00-1.))
 !           END IF
            Beta_un=wb*ENT(K,I)   
           EntW = exp(-2.*Beta_un*deltaZ) ! exp(-2.*(Wb*ENT(K,I))*deltaZ)
            IF (Beta_un>0) THEN
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
        edmf_a(K)=edmf_a(K)+UPA(K+1,I)
        edmf_w(K)=edmf_w(K)+UPA(K+1,I)*UPW(K+1,I)
        edmf_qt(K)=edmf_qt(K)+UPA(K+1,I) * (UPQT(K+1,I)-QT(K+1))
        edmf_thl(K)=edmf_thl(K)+UPA(K+1,I)* (UPTHL(K+1,I)-THL(K+1))
        edmf_ent(K)=edmf_ent(K)+UPA(K+1,I)*ENT(K+1,I)
        edmf_qc(K)=edmf_qc(K)+UPA(K+1,I)*UPQC(K+1,I)
        edmf_debug1(K)=edmf_debug1(K) + UPA(K+1,I) * (UPTHV(K+1,I)-THV(K+1))
        edmf_debug2(K)=QT(K+1) !debug2 is the mean qt to compare with actual output
        ! edmf_debug3(K)=edmf_debug3(K)+UPA(K+1,I)* UPTHL(K+1,I) !debug3 is the actual ud thl
        ! edmf_debug4(K)=THL(K+1) !debug4 is the mean thl to compare with actual output
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


END SUBROUTINE edmf_JPL

! ===================================================================
! ===================================================================
! This is the downdraft mass flux scheme - analogus to edmf_JPL but 
! flipped updraft to downdraft. This scheme is currently only tested
! for Stratocumulus cloud conditions. For a detailed desctiption of the
! model, see paper.

SUBROUTINE DDMF_JPL(kts,kte,dt,zw,p,                 &
              &u,v,th,thl,thv,tk,qt,qv,qc,           &
              &exner,                                &
              &ust,wthl,wqt,pblh,kpbl,               &
              &edmf_a_dd,edmf_w_dd, edmf_qt_dd,      &
              &edmf_thl_dd,edmf_ent_dd,edmf_qc_dd,   &
              &edmf_debug3,edmf_debug4,              &
              &sd_aw,sd_awthl,sd_awqt,               &
              &sd_awqv,sd_awqc,sd_awu,sd_awv,        &
              &sd_awqke,                             &
              &qc_bl1d,cldfra_bl1d,                  &
              &rthraten                              &
              )

        INTEGER, INTENT(IN) :: KTS,KTE,KPBL
        REAL,DIMENSION(KTS:KTE), INTENT(IN) :: U,V,TH,THL,TK,QT,QV,QC,THV,P,exner,rthraten
        ! zw .. heights of the downdraft levels (edges of boxes)
        REAL,DIMENSION(KTS:KTE+1), INTENT(IN) :: ZW
        REAL, INTENT(IN) :: DT,UST,WTHL,WQT,PBLH

    ! outputs - downdraft properties
        REAL,DIMENSION(KTS:KTE), INTENT(OUT) :: edmf_a_dd,edmf_w_dd,        &
                      & edmf_qt_dd,edmf_thl_dd, edmf_ent_dd,edmf_qc_dd,edmf_debug3,edmf_debug4

    ! outputs - variables needed for solver (sd_aw - sum ai*wi, sd_awphi - sum ai*wi*phii)
        REAL,DIMENSION(KTS:KTE+1) :: sd_aw, sd_awthl, sd_awqt, sd_awu, sd_awv, &
                                     sd_awqc, sd_awqv, sd_awqke, sd_aw2

        REAL,DIMENSION(KTS:KTE), INTENT(IN) :: qc_bl1d, cldfra_bl1d

        INTEGER, PARAMETER :: NDOWN=10, debug_mf=0 !fixing number of plumes to 10
    ! draw downdraft starting height randomly between cloud base and cloud top
        INTEGER, DIMENSION(1:NDOWN) :: DD_initK
        REAL   , DIMENSION(1:NDOWN) :: randNum
    ! downdraft properties
        REAL,DIMENSION(KTS:KTE+1,1:NDOWN) :: DOWNW,DOWNTHL,DOWNQT,DOWNQC,DOWNA,DOWNU,DOWNV,DOWNTHV

    ! entrainment variables
        REAl,DIMENSION(KTS+1:KTE+1,1:NDOWN) :: ENT,ENTf
        INTEGER,DIMENSION(KTS+1:KTE+1,1:NDOWN) :: ENTi

    ! internal variables
        INTEGER :: K,I, kminrad, qlTop, p700_ind, qlBase
        REAL :: wthv,wstar,qstar,thstar,sigmaW,sigmaQT,sigmaTH,z0, &
            pwmin,pwmax,wmin,wmax,wlv,wtv,went
        REAL :: B,QTn,THLn,THVn,QCn,Un,Vn,Wn2,EntEXP,EntW, Beta_dm, EntExp_M
        REAL :: jump_thetav, jump_qt, jump_thetal, refTHL, refTHV, refQT
    ! DD specific internal variables
        REAL :: minrad,zminrad, radflux, F0, wst_rad, wst_dd, deltaZ
        logical :: cloudflg

    ! VARIABLES FOR CHABOUREAU-BECHTOLD CLOUD FRACTION
        ! REAL,DIMENSION(KTS:KTE), INTENT(INOUT) :: vt, vq, sgm
        REAL :: sigq,xl,tlk,qsat_tl,rsl,cpm,a,qmq,mf_cf,diffqt,&
               Fng,qww,alpha,beta,bb,f,pt,t,q2p,b9,satvp,rhgrid

    ! w parameters
        REAL,PARAMETER :: &
            &Wa=1., &
            &Wb=1.5,&
            &Z00=100.
    ! entrainment parameters
        REAL,PARAMETER :: &
        & L0=80,&
        & ENT0=0.2

        pwmin=-3. ! drawing from the neagtive tail -3sigma to -1sigma
        pwmax=-1.

    ! initialize downdraft properties
        DOWNW=0.
        DOWNTHL=0.
        DOWNTHV=0.
        DOWNQT=0.
        DOWNA=0.
        DOWNU=0.
        DOWNV=0.
        DOWNQC=0.
        ENT=0.
        DD_initK=0

        edmf_a_dd  =0.
        edmf_w_dd  =0.
        edmf_qt_dd =0.
        edmf_thl_dd=0.
        edmf_ent_dd=0.
        edmf_qc_dd =0.
        edmf_debug3=0.
        edmf_debug4=0.

        sd_aw=0.
        sd_awthl=0.
        sd_awqt=0.
        sd_awqv=0.
        sd_awqc=0.
        sd_awu=0.
        sd_awv=0.
        sd_awqke=0.

    ! FIRST, CHECK FOR STRATOCUMULUS-TOPPED BOUNDARY LAYERS 
        cloudflg=.false.
        minrad=100.
        kminrad=kpbl
        zminrad=PBLH
        qlTop = 1 !initialize at 0
        qlBase = 1
        wthv=wthl+svp1*wqt
        do k = MAX(1,kpbl-2),kpbl+3
            if(qc(k).gt. 1.e-6 .AND. cldfra_bl1D(k).gt.0.5) then
               cloudflg=.true. !found Sc cloud
               qlTop = k ! index for Sc cloud top
            endif
        enddo

        do k = qlTop, kts, -1
            if(qc(k).gt. 1E-6) then
                qlBase = k ! index for Sc cloud base
            endif
        enddo
        qlBase = (qlTop+qlBase)/2 ! changed base to half way through the cloud

        call init_random_seed_1()
        call RANDOM_NUMBER(randNum)
        do i=1,NDOWN
            ! downdraft starts somewhere between cloud base to cloud top
            ! the probability is equally distributed
            DD_initK(i) = qlTop ! nint(randNum(i)*REAL(qlTop-qlBase)) + qlBase
        enddo

        ! LOOP RADFLUX
        F0 = 0.
        do k = 1, qlTop ! Snippet from YSU, YSU loops until qlTop - 1
           radflux = rthraten(k) * exner(k) ! Converts theta/s to temperature/s
           radflux = radflux * cp / g * ( p(k) - p(k+1) ) ! Converts temperature/s to W/m^2
           if ( radflux < 0.0 ) F0 = abs(radflux) + F0
        enddo
        F0 = max(F0, 1.0)
        !found Sc cloud and cloud not at surface, trigger downdraft
        if (cloudflg) then 

            !get entrainent coefficient
            do i=1,NDOWN
                do k=kts+1,kte
                  ENTf(k,i)=(ZW(k+1)-ZW(k))/L0
                 enddo
            enddo

            ! get Poisson P(dz/L0)
            call Poisson(1,NDOWN,kts+1,kte,ENTf,ENTi)


            ! entrainent: Ent=Ent0/dz*P(dz/L0)
            do i=1,NDOWN
                do k=kts+1,kte
                  ENT(k,i)=real(ENTi(k,i))*Ent0/(ZW(k+1)-ZW(k))
                enddo
            enddo            

            !!![EW: INVJUMP] find 700mb height then subtract trpospheric lapse rate!!!
            p700_ind = MINLOC(ABS(p-70000),1)!p1D is 70000
            jump_thetav = thv(p700_ind) - thv(1) - (thv(p700_ind)-thv(qlTop+3))/(ZW(p700_ind)-ZW(qlTop+3))*(ZW(p700_ind)-ZW(qlTop))
            jump_qt = qc(p700_ind) + qv(p700_ind) - qc(1) - qv(1)
            jump_thetal = thl(p700_ind) - thl(1) - (thl(p700_ind)-thl(qlTop+3))/(ZW(p700_ind)-ZW(qlTop+3))*(ZW(p700_ind)-ZW(qlTop))
            
            refTHL = thl(qlTop) !sum(thl(1:qlTop)) / (qlTop) ! avg over BL for now or just at qlTop
            refTHV = thv(qlTop) !sum(thv(1:qlTop)) / (qlTop)
            refQT = qt(qlTop) !sum(qt(1:qlTop)) / (qlTop)            
            wst_rad = ( g * zw(qlTop) * F0 / (refTHL * rhoair0 * cp) ) ** (0.333) ! wstar_rad, following Lock and MacVean (1999a)
            wst_rad = max(wst_rad, 0.1)
            wstar = max(0.,(g/thv(1)*wthv*pblh)**(1./3.))
            went = thv(1) / ( g * jump_thetav * zw(qlTop) ) * (0.15 * (wstar**3 + 5*ust**3) + 0.35 * wst_rad**3 )
            qstar = abs(went*jump_qt/wst_rad)
            thstar = F0/rhoair0/cp/wst_rad - went*jump_thetav/wst_rad
            wst_dd = (0.15 * (wstar**3 + 5*ust**3) + 0.35 * wst_rad**3 ) ** (0.333) !wstar_dd - mix rad + surface wst
            sigmaW = 0.2*wst_dd! 0.8*wst_dd!wst_rad !tuning parameter ! 0.5 was good
            sigmaQT = 40  * qstar ! 50 was good
            sigmaTH = 1.0 * thstar! 0.5 was good

            wmin=sigmaW*pwmin
            wmax=sigmaW*pwmax

            do I=1,NDOWN !downdraft now starts at different height
                wlv=wmin+(wmax-wmin)/NDOWN*(i-1)
                wtv=wmin+(wmax-wmin)/NDOWN*i

                DOWNW(DD_initK(I),I)=0.5*(wlv+wtv)
                DOWNA(DD_initK(I),I)=0.5*ERF(wtv/(sqrt(2.)*sigmaW))-0.5*ERF(wlv/(sqrt(2.)*sigmaW))

                DOWNU(DD_initK(I),I)=U(DD_initK(I))
                DOWNV(DD_initK(I),I)=V(DD_initK(I))

                refTHL = 0.5 * (thl(DD_initK(I))+thl(DD_initK(I)-1)) !reference now depends on where dd starts
                refTHV = 0.5 * (thv(DD_initK(I))+thv(DD_initK(I)-1))
                refQT  = 0.5 * (qt(DD_initK(I))+qt(DD_initK(I)-1))

                DOWNQC(DD_initK(I),I) = 0.
                DOWNQT(DD_initK(I),I) = refQT  + 0.5  *DOWNW(DD_initK(I),I)*sigmaQT/sigmaW
                DOWNTHV(DD_initK(I),I)= refTHV + 0.01 *DOWNW(DD_initK(I),I)*sigmaTH/sigmaW ! close to no difference

                call condensation_edmf_r(DOWNQT(DD_initK(I),I),DOWNTHL(DD_initK(I),I),(P(DD_initK(I))+P(DD_initK(I)-1))/2.,ZW(DD_initK(I)),DOWNTHV(DD_initK(I),I),DOWNQC(DD_initK(I),I))

                ! DOWNTHL(DD_initK(I),I)=DOWNTHV(DD_initK(I),I)/(1.+svp1*DOWNQT(DD_initK(I),I))
                ! print *, "Plume ", I, " DOWNQT = ", DOWNQT(DD_initK(I),I), " refQT = ", refQT, &
                !          " DOWNTHV = ", DOWNTHV(DD_initK(I),I), " refTHV = ", refTHV, &
                !          " DOWNTHL = ", DOWNTHL(DD_initK(I),I), " refTHL = ", refTHL
            enddo
             ! print *, "RefTHV = ", refTHV, " sigmaQT = ", sigmaQT, " sigmaTH = ", sigmaTH, " sigmaW = ", sigmaW, " total = ", -0.2*DOWNW(qlTop,1)*sigmaTH/sigmaW
            print*, "[EWDD] wst = ", wstar, " wst_rad = ", wst_rad, " wst_dd = ", wst_dd, " went = ", went
            ! do integration downdraft
            DO I=1,NDOWN
                DO k=DD_initK(I)-1,KTS+1,-1 !starting one point below cloud top
                    deltaZ = ZW(k+1)-ZW(k)
                    EntExp=exp(-ENT(K,I)*deltaZ)
                    EntExp_M=exp(-ENT(K,I)/3.*deltaZ)

                    QTn=DOWNQT(K+1,I)+(QT(K)-DOWNQT(K+1,I))*(1.-EntExp)
                    THLn=DOWNTHL(K+1,I)+(THL(K)-DOWNTHL(K+1,I))*(1.-EntExp)
                    Un=DOWNU(K+1,I)+(U(K)-DOWNU(K+1,I))*(1.-EntExp_M)
                    Vn=DOWNV(K+1,I)+(V(K)-DOWNV(K+1,I))*(1.-EntExp_M)
                    ! get thvn,qcn
                    call condensation_edmf(QTn,THLn,(P(K)+P(K-1))/2.,ZW(k),THVn,QCn)
                    B=g*(0.5*(THVn+DOWNTHV(k+1,I))/THV(k)-1.)

                    Beta_dm = 2*Wb*ENT(K,I) + 0.5/(ZW(k)-deltaZ) * max(1. - exp((ZW(k) -deltaZ)/Z00 - 1. ) , 0.)
                    EntW=exp(-Beta_dm * deltaZ)
                    if (Beta_dm >0) then
                        Wn2=DOWNW(K+1,I)**2*EntW - Wa*B/Beta_dm * (1. - EntW)
                    else
                        Wn2=DOWNW(K+1,I)**2 - 2.*Wa*B*deltaZ
                    end if
                     print *, "Plume number = ", I, " k = ", k, " z = ", ZW(k), " entw = ", ENT(K,I)
                     print *, "downthv = ", THVn, " thv = ", thv(k)
                     print *, "downthl = ", THLn, " thl = ", thl(k)
                     print *, "downqt  = ", QTn , " qt  = ", qt(k)
                     print *, "Beta = ", Beta_dm, " Wn2 = ", Wn2, " Bouy = ", B

                    IF (Wn2 .gt. 0.) THEN !terminate when velocity is too small
                        DOWNW(K,I)=-sqrt(Wn2)
                        DOWNTHV(K,I)=THVn
                        DOWNTHL(K,I)=THLn
                        DOWNQT(K,I)=QTn
                        DOWNQC(K,I)=QCn
                        DOWNU(K,I)=Un
                        DOWNV(K,I)=Vn
                        DOWNA(K,I)=DOWNA(K+1,I)
                    ELSE
                          exit
                    ENDIF
                ENDDO
            ENDDO
        endif ! end cloud flag

        DOWNW(1,:) = 0. !make sure downdraft does not go to the surface
        DOWNA(1,:) = 0.
        ! do I=1,NDOWN
        !     DO k=qlTop,1,-1
        !         if (DOWNA(k,I)>0.) then
        !             print *, "[EWDD] Plume number = ", I, " k = ", k, " Downdraft w = ", DOWNW(k,I), " area = ", DOWNA(k,I)
        !             print *, "[EWDD] qt = ", DOWNQT(k,I), " thv = ", DOWNTHV(k,I), " thl = ", DOWNTHL(k,I)
        !             print *, "[EWDD] refqt = ", qt(k), " refthv = ", thv(k), " refthl = ", thl(k)
        !         endif
        !     ENDDO
        ! enddo
        ! print*, "[EWDD] wst = ", wstar, " wst_rad = ", wst_rad, " wst_dd = ", wst_dd, " went = ", went

        ! Combine both moist and dry plume, write as one averaged plume
        ! Even though downdraft starts at different height, average all up to qlTop
        DO k=qlTop,KTS,-1
            DO I=1,NDOWN
                IF(I > NDOWN) exit
                edmf_a_dd(K)=edmf_a_dd(K)+DOWNA(K-1,I)
                edmf_w_dd(K)=edmf_w_dd(K)+DOWNA(K-1,I)*DOWNW(K-1,I)
                edmf_qt_dd(K)=edmf_qt_dd(K)+DOWNA(K-1,I)* (DOWNQT(K-1,I)-QT(K-1))
                edmf_thl_dd(K)=edmf_thl_dd(K)+DOWNA(K-1,I)* (DOWNTHL(K-1,I)-THL(K-1))
                edmf_ent_dd(K)=edmf_ent_dd(K)+DOWNA(K-1,I)*ENT(K-1,I)
                edmf_qc_dd(K)=edmf_qc_dd(K)+DOWNA(K-1,I)*DOWNQC(K-1,I)
                edmf_debug3(K)=edmf_debug3(K)+DOWNA(K-1,I)* (DOWNTHV(K-1,I)-THV(K-1))
                edmf_debug4(K)=edmf_debug4(K)+THL(k)
            ENDDO

            IF (edmf_a_dd(k)>0.) THEN
                edmf_w_dd(k)=edmf_w_dd(k)/edmf_a_dd(k)
                edmf_qt_dd(k)=edmf_qt_dd(k)/edmf_a_dd(k)
                edmf_thl_dd(k)=edmf_thl_dd(k)/edmf_a_dd(k)
                edmf_ent_dd(k)=edmf_ent_dd(k)/edmf_a_dd(k)
                edmf_qc_dd(k)=edmf_qc_dd(k)/edmf_a_dd(k)
                edmf_debug3(k)=edmf_debug3(k)/edmf_a_dd(k)
            ENDIF
        ENDDO        

          !
          ! computing variables needed for solver
          !

        DO k=KTS,qlTop
            DO I=1,NDOWN
                sd_aw(k)=sd_aw(k)+DOWNA(k,I)*DOWNW(k,I)
                sd_awthl(k)=sd_awthl(k)+DOWNA(k,i)*DOWNW(k,I)*DOWNTHL(k,I)
                sd_awqt(k)=sd_awqt(k)+DOWNA(k,i)*DOWNW(k,I)*DOWNQT(k,I)
                sd_awqc(k)=sd_awqc(k)+DOWNA(k,i)*DOWNW(k,I)*DOWNQC(k,I)
                sd_awu(k)=sd_awu(k)+DOWNA(k,i)*DOWNW(k,I)*DOWNU(k,I)
                sd_awv(k)=sd_awv(k)+DOWNA(k,i)*DOWNW(k,I)*DOWNV(k,I)
            ENDDO
            sd_awqv(k) = sd_awqt(k)  - sd_awqc(k)
        ENDDO

        ! ! This last piece is from STEM_MF
        ! DO K=KTS,qlTop
        !     IF(edmf_qc_dd(k)>0.0)THEN
        !         satvp = 3.80*exp(17.27*(th(k)-273.)/ &
        !            (th(k)-36.))/(.01*p(k))
        !         rhgrid = max(.01,MIN( 1., qv(k) /satvp))

        !         !COMPUTE CLDFRA & QC_BL FROM MASS-FLUX SCHEME and recompute vt & vq

        !         xl = xl_blend(tk(k))                ! obtain blended heat capacity 
        !         tlk = thl(k)*(p(k)/p1000mb)**rcp    ! recover liquid temp (tl) from thl
        !         qsat_tl = qsat_blend(tlk,p(k))      ! get saturation water vapor mixing ratio
        !                                         !   at tl and p
        !         rsl = xl*qsat_tl / (r_v*tlk**2)     ! slope of C-C curve at t = tl
        !                                         ! CB02, Eqn. 4
        !         cpm = cp + qt(k)*cpv                ! CB02, sec. 2, para. 1
        !         a   = 1./(1. + xl*rsl/cpm)          ! CB02 variable "a"
        !         b9  = a*rsl                         ! CB02 variable "b" 

        !         q2p  = xlvcp/exner(k)
        !         pt = thl(k) +q2p*edmf_qc_dd(k) ! potential temp
        !         bb = b9*tk(k)/pt ! bb is "b9" in BCMT95.  Their "b9" differs from
        !                    ! "b9" in CB02 by a factor
        !                    ! of T/theta.  Strictly, b9 above is formulated in
        !                    ! terms of sat. mixing ratio, but bb in BCMT95 is
        !                    ! cast in terms of sat. specific humidity.  The
        !                    ! conversion is neglected here.
        !         qww   = 1.+0.61*qt(k)
        !         alpha = 0.61*pt
        !         t     = th(k)*exner(k)
        !         beta  = pt*xl/(t*cp) - 1.61*pt
        !         !Buoyancy flux terms have been moved to the end of this section...

        !         !Now calculate convective component of the cloud fraction:
        !         if (a > 0.0) then
        !             f = MIN(1.0/a, 4.0)              ! f is vertical profile scaling function (CB2005)
        !         else
        !             f = 1.0
        !         endif
        !         sigq = 9.E-3 * edmf_a_dd(k) * edmf_w_dd(k) * f ! convective component of sigma (CB2005)
        !         !sigq = MAX(sigq, 1.0E-4)         
        !         sigq = SQRT(sigq**2 + sgm(k)**2)    ! combined conv + stratus components

        !         qmq = a * (qt(k) - qsat_tl)         ! saturation deficit/excess;
        !                                         !   the numerator of Q1
        !         mf_cf = min(max(0.5 + 0.36 * atan(1.55*(qmq/sigq)),0.02),0.6)
        !         IF ( debug_code ) THEN
        !             print*,"In MYNN, EDMF JPL"
        !             print*,"  CB: qt=",qt(k)," qsat=",qsat_tl," satdef=",qt(k) - qsat_tl
        !             print*,"  CB: sigq=",sigq," qmq=",qmq," tlk=",tlk
        !             print*,"  CB: mf_cf=",mf_cf," cldfra_bl=",cldfra_bl1d(k)," edmf_a_dd=",edmf_a_dd(k)
        !         ENDIF

        !         IF (rhgrid >= .93) THEN
        !             !IN high RH, defer to stratus component if > convective component
        !             cldfra_bl1d(k) = MAX(mf_cf, cldfra_bl1d(k))
        !             IF (cldfra_bl1d(k) > edmf_a_dd(k)) THEN
        !                 qc_bl1d(k) = edmf_qc_dd(k)*edmf_a_dd(k)/cldfra_bl1d(k)
        !             ELSE
        !                 cldfra_bl1d(k)=edmf_a_dd(k)
        !                 qc_bl1d(k) = edmf_qc_dd(k)
        !             ENDIF
        !         ELSE
        !             IF (mf_cf > edmf_a_dd(k)) THEN
        !                 cldfra_bl1d(k) = mf_cf
        !                 qc_bl1d(k) = edmf_qc_dd(k)*edmf_a_dd(k)/mf_cf
        !             ELSE
        !                 cldfra_bl1d(k)=edmf_a_dd(k)
        !                 qc_bl1d(k) = edmf_qc_dd(k)
        !             ENDIF
        !         ENDIF
        !         !Now recalculate the terms for the buoyancy flux for mass-flux clouds:
        !         !See mym_condensation for details on these formulations.  The
        !         !cloud-fraction bounding was added to improve cloud retention,
        !         !following RAP and HRRR testing.
        !         Fng = 2.05 ! the non-Gaussian transport factor (assumed constant)
        !         vt(k) = qww   - MIN(0.20,cldfra_bl1D(k))*beta*bb*Fng - 1.
        !         vq(k) = alpha + MIN(0.20,cldfra_bl1D(k))*beta*a*Fng  - tv0
        !     ENDIF
        ! ENDDO
END SUBROUTINE DDMF_JPL


!################################
!################################
!
!  GFDL model interface, yhc
!
!################################
!################################

subroutine edmf_mynn_driver ( &
              is, ie, js, je, npz, Time_next, dt, lon, lat, frac_land, area, u_star,  &
              b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, &
              do_edmf_mynn_diagnostic, &
              udt, vdt, tdt, rdt, rdiag)

!---------------------------------------------------------------------
! Arguments (Intent in)  
!   Descriptions of the input arguments are in subroutine edmf_alloc
!---------------------------------------------------------------------
  integer, intent(in)                   :: is, ie, js, je, npz
  type(time_type), intent(in)           :: Time_next
  real,    intent(in)                   :: dt
  real,    intent(in), dimension(:,:)   :: &
   lon, lat, &  ! longitude and latitude in radians
   frac_land, area, u_star, b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux
  type(physics_input_block_type)        :: Physics_input_block
  logical, intent(in)                   :: do_edmf_mynn_diagnostic

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
!---------------------------------------------------------------------
  real, intent(inout), dimension(:,:,:) :: &
    udt, vdt, tdt
  real, intent(inout), dimension(:,:,:,:) :: &
    rdt
  !real, intent(inout), dimension(:,:,:,:) :: &
  real, intent(inout), dimension(:,:,:,ntp+1:) :: &
    rdiag

!---------------------------------------------------------------------
! local variables  
!---------------------------------------------------------------------
  type(edmf_input_type)  :: Input_edmf
  type(edmf_output_type) :: Output_edmf
  type(am4_edmf_output_type) :: am4_Output_edmf

  logical used
  logical do_writeout_column
  real    :: lat_lower, lat_upper, lon_lower, lon_upper, lat_temp, lon_temp
  real    :: tt1
  integer :: i,j,k
!-------------------------

!! debug01
!write(6,*) 'edmf_mynn, beginning'
!!write(6,*) 'Physics_input_block%omega',Physics_input_block%omega
write(6,*) 'initflag,',initflag
write(6,*) 'nQke, rdiag(:,:,:,nQke)',nQke, rdiag(:,:,:,nQke)
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
              b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, rdiag, &
              rdiag(:,:,:,nQke), rdiag(:,:,:,nel_pbl), rdiag(:,:,:,ncldfra_bl), rdiag(:,:,:,nqc_bl), rdiag(:,:,:,nSh3D), &
              Input_edmf, Output_edmf, am4_Output_edmf)

!---------------------------------------------------------------------
! run the EDMF-MYNN program
!---------------------------------------------------------------------

   call mynn_bl_driver(            &
       &initflag=initflag,grav_settling=grav_settling,         &
       &delt=Input_edmf%delt,dz=Input_edmf%dz,dx=Input_edmf%dx,znt=Input_edmf%znt,                 &
       &u=Input_edmf%u,v=Input_edmf%v,w=Input_edmf%w,th=Input_edmf%th,qv=Input_edmf%qv,             &
       &qc=Input_edmf%qc,qi=Input_edmf%qi,qni=Input_edmf%qni,qnc=Input_edmf%qnc,                    &
       &p=Input_edmf%p,exner=Input_edmf%exner,rho=Input_edmf%rho,T3D=Input_edmf%T3D,                &
       &xland=Input_edmf%xland,ts=Input_edmf%ts,qsfc=Input_edmf%qsfc,qcg=Input_edmf%qcg,ps=Input_edmf%ps,           &
       &ust=Input_edmf%ust,ch=Input_edmf%ch,hfx=Input_edmf%hfx,qfx=Input_edmf%qfx,rmol=Input_edmf%rmol,wspd=Input_edmf%wspd,       &
       &uoce=Input_edmf%uoce,voce=Input_edmf%voce,                      & 
       &vdfg=Input_edmf%vdfg,                           & 
       &qke=Output_edmf%Qke,                    &
       &Tsq=Output_edmf%Tsq,Qsq=Output_edmf%Qsq,Cov=Output_edmf%Cov,                    &
       &RUBLTEN=Output_edmf%RUBLTEN,RVBLTEN=Output_edmf%RVBLTEN,RTHBLTEN=Output_edmf%RTHBLTEN,       &
       &RQVBLTEN=Output_edmf%RQVBLTEN,RQCBLTEN=Output_edmf%RQCBLTEN,RQIBLTEN=Output_edmf%RQIBLTEN,     &
       &RQNIBLTEN=Output_edmf%RQNIBLTEN,                      &
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
       &edmf_thl=Output_edmf%edmf_thl,edmf_ent=Output_edmf%edmf_ent,edmf_qc=Output_edmf%edmf_qc,      &
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

  !--- SCM, set initfla
  if (initflag == 1) then
    initflag = 0          ! no initialization
  endif

!! debug01, check semi-prognostic variables in offline code
!  call random_number (Output_edmf%Qke)
!  call random_number (Output_edmf%el_pbl     )
!  call random_number (Output_edmf%cldfra_bl  )
!  call random_number (Output_edmf%qc_bl      )
!  call random_number (Output_edmf%Sh3D       )

!---------------------------------------------------------------------
! process the outputs from the EDMF-MYNN program
!---------------------------------------------------------------------

  !--- convert Output_edmf to am4_Output_edmf
  call convert_edmf_to_am4_array (size(Physics_input_block%t,1), size(Physics_input_block%t,2), size(Physics_input_block%t,3), &
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

!  !<-- test tendency purposes
!  am4_Output_edmf%udt_edmf = 1./dt
!  am4_Output_edmf%vdt_edmf = 2./dt
!  am4_Output_edmf%tdt_edmf = -1./dt
!  am4_Output_edmf%qdt_edmf = -2./dt
!  !--> test tendency purposes

!<-- debug 
  tt1 = tdt_max / 86400.  ! change unit from K/day to K/sec
  do i=1,size(am4_Output_edmf%tdt_edmf,1)
  do j=1,size(am4_Output_edmf%tdt_edmf,2)
  do k=1,size(am4_Output_edmf%tdt_edmf,3)
    if ( abs(am4_Output_edmf%tdt_edmf(i,j,k)) .ge. tt1 ) then
      write(6,*) 'edmf, >tdt_max,i,j,lat,lon,',tdt_max,i,j,lat(i,j),lon(i,j)

      if (do_limit_tdt) then
        if (am4_Output_edmf%tdt_edmf(i,j,k).ge.0.) then
          am4_Output_edmf%tdt_edmf(i,j,k) = tdt_limit / 86400.
        else
          am4_Output_edmf%tdt_edmf(i,j,k) = -1.*tdt_limit / 86400.
        endif
      endif
    endif
  enddo
  enddo
  enddo
!-->

  !--- updated tendencies
  if (.not.do_edmf_mynn_diagnostic) then
    udt(:,:,:) = udt(:,:,:) + am4_Output_edmf%udt_edmf(:,:,:)
    vdt(:,:,:) = vdt(:,:,:) + am4_Output_edmf%vdt_edmf(:,:,:)
    tdt(:,:,:) = tdt(:,:,:) + am4_Output_edmf%tdt_edmf(:,:,:)

    rdt(:,:,:,nsphum) = rdt(:,:,:,nsphum) + am4_Output_edmf%qdt_edmf(:,:,:)
  end if

  !--- write out EDMF-MYNN input and output fields for debugging purpose
  call edmf_writeout_column ( &
              do_writeout_column, &
              is, ie, js, je, npz, Time_next, dt, lon, lat, frac_land, area, u_star,  &
              b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, &
              rdiag(:,:,:,nQke), rdiag(:,:,:,nel_pbl), rdiag(:,:,:,ncldfra_bl), rdiag(:,:,:,nqc_bl), rdiag(:,:,:,nSh3D), &
              Input_edmf, Output_edmf, am4_Output_edmf, rdiag)

!if (do_writeout_column) then
!  write(6,*) 'initflag, out',initflag
!  write(6,*) 'rdiag(:,:,:,nQke), out',rdiag(ii_write,jj_write,:,nQke)
!endif

!!---------------------------------------------------------------------
!! write out fields to history files
!!---------------------------------------------------------------------
!
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
!        used = send_data (id_w1_thv1_surf_star, Input_edmf%w1_th1_surf_star, Time_next, is, js )
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
!      if ( id_tke > 0) then
!        used = send_data (id_tke, am4_Output_edmf%tke, Time_next, is, js, 1 )
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
!      if ( id_qcdt_edmf > 0) then
!        used = send_data (id_qcdt_edmf, am4_Output_edmf%qcdt_edmf, Time_next, is, js, 1 )
!      endif
!
!!------- updraft area (units: none) at full level -------
!      if ( id_edmf_a > 0) then
!        used = send_data (id_edmf_a, am4_Output_edmf%edmf_a, Time_next, is, js, 1 )
!      endif
!
!!------- vertical velocity of updrafts (units: m/s) at full level -------
!      if ( id_edmf_w > 0) then
!        used = send_data (id_edmf_w, am4_Output_edmf%edmf_w, Time_next, is, js, 1 )
!      endif
!
!!------- qt in updrafts (units: kg/kg) at full level -------
!      if ( id_edmf_qt > 0) then
!        used = send_data (id_edmf_qt, am4_Output_edmf%edmf_qt, Time_next, is, js, 1 )
!      endif
!
!!------- thl in updrafts (units: K) at full level -------
!      if ( id_edmf_thl > 0) then
!        used = send_data (id_edmf_thl, am4_Output_edmf%edmf_thl, Time_next, is, js, 1 )
!      endif
!
!!------- entrainment in updrafts (units: 1/m) at full level -------
!      if ( id_edmf_ent > 0) then
!        used = send_data (id_edmf_ent, am4_Output_edmf%edmf_ent, Time_next, is, js, 1 )
!      endif
!
!!------- qc in updrafts (units: kg/kg) at full level -------
!      if ( id_edmf_qc > 0) then
!        used = send_data (id_edmf_qc, am4_Output_edmf%edmf_qc, Time_next, is, js, 1 )
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
              b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, rdiag, &
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
! qc ... cloud water mixing ratio [kg kg^-1], half levels 
! qi ... ice water mixing ratio [kg kg^-1], half levels 
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
!         RQCBLTEN ....  Qc tendency due to EDMF parameterization (kg/kg/s), set to 0. on the input
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
    dz_host, u_host, v_host, w_host, th_host, qv_host, p_host, exner_host, rho_host, T3D_host, qc_host, ql_host, qi_host, qnc_host, qni_host, &
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
  allocate (Output_edmf%RQCBLTEN    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RQCBLTEN    = 0.
  allocate (Output_edmf%RQIBLTEN    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RQIBLTEN    = 0.
  allocate (Output_edmf%RQNIBLTEN   (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RQNIBLTEN   = 0.
  allocate (Output_edmf%RTHRATEN    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%RTHRATEN    = 0.
  allocate (Output_edmf%edmf_a      (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_a      = 0.
  allocate (Output_edmf%edmf_w      (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_w      = 0.
  allocate (Output_edmf%edmf_qt     (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_qt     = 0.
  allocate (Output_edmf%edmf_thl    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_thl    = 0.
  allocate (Output_edmf%edmf_ent    (IMS:IME,KMS:KME,JMS:JME))  ; Output_edmf%edmf_ent    = 0.
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
  allocate (am4_Output_edmf%qcdt_edmf   (ix,jx,kx))  ; am4_Output_edmf%qcdt_edmf   = 0.
  allocate (am4_Output_edmf%edmf_a      (ix,jx,kx))  ; am4_Output_edmf%edmf_a      = 0.
  allocate (am4_Output_edmf%edmf_w      (ix,jx,kx))  ; am4_Output_edmf%edmf_w      = 0.
  allocate (am4_Output_edmf%edmf_qt     (ix,jx,kx))  ; am4_Output_edmf%edmf_qt     = 0.
  allocate (am4_Output_edmf%edmf_thl    (ix,jx,kx))  ; am4_Output_edmf%edmf_thl    = 0.
  allocate (am4_Output_edmf%edmf_ent    (ix,jx,kx))  ; am4_Output_edmf%edmf_ent    = 0.
  allocate (am4_Output_edmf%edmf_qc     (ix,jx,kx))  ; am4_Output_edmf%edmf_qc     = 0.
  allocate (am4_Output_edmf%thl_edmf    (ix,jx,kx))  ; am4_Output_edmf%thl_edmf    = 0.
  allocate (am4_Output_edmf%qt_edmf     (ix,jx,kx))  ; am4_Output_edmf%qt_edmf     = 0.

  !allocate (am4_Output_edmf%         (ix,jx,kx))  ; am4_Output_edmf%         = 0.

!-------------------------------------------------------------------------
! set values for Input_edmf 
!-------------------------------------------------------------------------

  !--- 0D variables
  Input_edmf%delt = dt
  Input_edmf%dx   = sqrt(area(1,1))

  !--- 3-D variable
  u_host     (:,:,:) = Physics_input_block%u
  v_host     (:,:,:) = Physics_input_block%v
  omega_host (:,:,:) = Physics_input_block%omega
  qv_host    (:,:,:) = Physics_input_block%q(:,:,:,nsphum)
  qc_host    (:,:,:) = Physics_input_block%q(:,:,:,nql)
  ql_host    (:,:,:) = Physics_input_block%q(:,:,:,nql)
  qi_host    (:,:,:) = Physics_input_block%q(:,:,:,nqi)
  p_host     (:,:,:) = Physics_input_block%p_full
  T3D_host   (:,:,:) = Physics_input_block%t

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
  tv_ref (:,:) = t_ref(:,:)*(1+d608*q_ref(:,:))
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
  Input_edmf%ts    = 0.   ! no need in JPL EDMF scheme  
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
    Input_edmf%qnc   (i,kk,j) = qnc_host   (i,j,k)
    Input_edmf%qni   (i,kk,j) = qni_host   (i,j,k)
  enddo  ! end loop of k
  enddo  ! end loop of j
  enddo  ! end loop of i

  ! diagnostic purpose
  Input_edmf%u_star_star = u_star_star (:,:)
  Input_edmf%shflx_star = shflx_star (:,:)
  Input_edmf%lhflx_star = lhflx_star (:,:)
  Input_edmf%w1_thv1_surf_star = w1_thv1_surf_star (:,:)
  Input_edmf%w1_th1_surf_star = w1_th1_surf_star (:,:)
  Input_edmf%w1_q1_surf_star = w1_q1_surf_star (:,:)
  Input_edmf%Obukhov_length_star = Obukhov_length_star (:,:)

  Input_edmf%u_star_updated = u_star_updated (:,:)
  Input_edmf%shflx_updated = shflx_updated (:,:)
  Input_edmf%lhflx_updated = lhflx_updated (:,:)
  Input_edmf%w1_thv1_surf_updated = w1_thv1_surf_updated (:,:)
  Input_edmf%w1_th1_surf_updated = w1_th1_surf_updated (:,:)
  Input_edmf%w1_q1_surf_updated = w1_q1_surf_updated (:,:)
  Input_edmf%Obukhov_length_updated = Obukhov_length_updated (:,:)

!-------------------------------------------------------------------------
! set values for Output_edmf 
!-------------------------------------------------------------------------

  ! semi-prognostic variables
  do i=1,ix    
  do j=1,jx    
  do k=1,kx 
    kk=kx-k+1   
    Input_edmf%Qke        (i,kk,j) = Qke       (i,j,k) 
    Input_edmf%el_pbl     (i,kk,j) = el_pbl    (i,j,k)
    Input_edmf%cldfra_bl  (i,kk,j) = cldfra_bl (i,j,k)
    Input_edmf%qc_bl      (i,kk,j) = qc_bl     (i,j,k)
    Input_edmf%Sh3D       (i,kk,j) = Sh3D      (i,j,k)

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
  deallocate (Output_edmf%RQCBLTEN    )  
  deallocate (Output_edmf%RQIBLTEN    )  
  deallocate (Output_edmf%RQNIBLTEN   )  
  deallocate (Output_edmf%RTHRATEN    )  
  deallocate (Output_edmf%edmf_a      )  
  deallocate (Output_edmf%edmf_w      )  
  deallocate (Output_edmf%edmf_qt     )  
  deallocate (Output_edmf%edmf_thl    )  
  deallocate (Output_edmf%edmf_ent    )  
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
  deallocate (am4_Output_edmf%qcdt_edmf   )  
  deallocate (am4_Output_edmf%edmf_a      )
  deallocate (am4_Output_edmf%edmf_w      )
  deallocate (am4_Output_edmf%edmf_qt     )
  deallocate (am4_Output_edmf%edmf_thl    )
  deallocate (am4_Output_edmf%edmf_ent    )
  deallocate (am4_Output_edmf%edmf_qc     )
  deallocate (am4_Output_edmf%thl_edmf    )
  deallocate (am4_Output_edmf%qt_edmf     )
  !deallocate (am4_Output_edmf%         )  

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
              b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, &
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
  !real, dimension (:,:,:,:) , intent(in) :: rdiag
  real, intent(inout), dimension(:,:,:,ntp+1:) :: &
    rdiag

  logical, intent(in) :: do_writeout_column

!---------------------------------------------------------------------
! local variables  
!---------------------------------------------------------------------
  integer :: ix, jx, kx, nt, kxp
  integer :: i, j, k, kk
  real    :: lat_lower, lat_upper, lon_lower, lon_upper, lat_temp, lon_temp

  real, dimension(size(Physics_input_block%t,3)) ::  &
    var_temp1
!-------------------------------------------------------------------------
!  define input array sizes.
!-------------------------------------------------------------------------
  ix = size(Physics_input_block%t,1)
  jx = size(Physics_input_block%t,2)
  kx = size(Physics_input_block%t,3)
  nt = size(Physics_input_block%q,4)

  kxp = kx + 1

!!-------------------------------------------------------------------------
!!  determine whether writing out the selected column
!!-------------------------------------------------------------------------
!  do_writeout_column = .false.
!  if (do_writeout_column_nml) then
!
!    !--- for global simulations
!    if (ii_write.ne.-999 .and. jj_write.ne.-999) then
!      do_writeout_column = .true.
!
!      if (lat_write.ne.-999.99 .and. lon_write.ne.-999.99) then
!
!        lat_lower = lat_write - lat_range
!        lat_upper = lat_write + lat_range
!        lon_lower = lon_write - lon_range
!        lon_upper = lon_write + lon_range
!
!        if (lat_lower.gt.lat_upper) then
!          lat_temp  = lat_upper
!          lat_upper = lat_lower
!          lat_lower = lat_temp
!        endif
!
!        if (lon_lower.gt.lon_upper) then
!          lon_temp  = lon_upper
!          lon_upper = lon_lower
!          lon_lower = lon_temp
!        endif
!
!        if (lat (ii_write,jj_write).gt.lat_lower .and. lat (ii_write,jj_write).lt.lat_upper .and. &
!            lon (ii_write,jj_write).gt.lon_lower .and. lon (ii_write,jj_write).lt.lon_upper ) then
!          do_writeout_column = .true.
!        else
!          do_writeout_column = .false.
!        endif
!      endif
!    endif
!
!    !--- SCM
!    if (ii_write.eq.0 .and. jj_write.eq.0) then
!      do_writeout_column = .true.
!      ii_write = 1
!      jj_write = 1
!    endif
!  endif  ! end if of do_writeout_column_nml

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
        write(6,3001) '  p_half = (/',Physics_input_block%p_half(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; height at half level above the surface (m)'
        write(6,3001) '  z_half = (/',Physics_input_block%z_half(ii_write,jj_write,:) - Physics_input_block%z_half(ii_write,jj_write,kxp)
        write(6,*)    ''
        write(6,*)    '; pressure at full level (Pa)'
        write(6,3001) '  p_full = (/',Physics_input_block%p_full(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; height at full level above the surface (m)'
        write(6,3001) '  z_full = (/',Physics_input_block%z_full(ii_write,jj_write,:) - Physics_input_block%z_half(ii_write,jj_write,kxp)
        write(6,*)    ''
        write(6,*)    '; actual height at half level (m)'
        write(6,3001) '  z_half_actual = (/',Physics_input_block%z_half(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; actual height at full level (m)'
        write(6,3001) '  z_full_actual = (/',Physics_input_block%z_full(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; zonal wind velocity at full levels (m/s)'
        write(6,3001) '  uu = (/'    ,Physics_input_block%u(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; meridional wind velocity at full levels (m/s)'
        write(6,3001) '  vv = (/'    ,Physics_input_block%v(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; vertical velocity (Pa/s)'
        write(6,3002) '  omega = (/'    ,Physics_input_block%omega(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; temperatur at full levels (K)'
        write(6,3001) '  tt = (/'    ,Physics_input_block%t(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; potential temperatur at full levels (K)'
           do k=1,kx
             kk=kx-k+1
             var_temp1(k) = Input_edmf%th(ii_write,kk,jj_write)
           enddo
        write(6,3001) '  th = (/'    ,var_temp1(:)
        write(6,*)    ''
        write(6,*)    '; ice-liquid water potential temperatur at full levels (K)'
        write(6,3001) '  thl = (/'    , am4_Output_edmf%thl_edmf(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '; specific humidity at full levels (kg/kg)'
        write(6,3002) '  qq = (/'    ,Physics_input_block%q(ii_write,jj_write,:,nsphum)
        write(6,*)    ''
        write(6,*)    '; cloud liquid water mixing ratio at full levels (kg/kg)'
        write(6,3002) '  ql = (/'    ,Physics_input_block%q(ii_write,jj_write,:,nql)
        write(6,*)    ''
        write(6,*)    '; cloud ice water mixing ratio at full levels (kg/kg)'
        write(6,3002) '  qi = (/'    ,Physics_input_block%q(ii_write,jj_write,:,nqi)
        write(6,*)    ''
        write(6,*)    '; total water mixing ratio (qv+ql+qi) at full levels (kg/kg)'
        write(6,3002) '  qt = (/'    ,am4_Output_edmf%qt_edmf(ii_write,jj_write,:)
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
        write(6,3002) ' qcdt_edmf = (/'    ,am4_Output_edmf%qcdt_edmf (ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    ';=============='
        write(6,*)    ';=============='
        write(6,*)    ''
        write(6,*)    ';   rdiag after mynn'
        write(6,*)    ''
        write(6,*)    ';=============='
        write(6,*)    ';=============='
        write(6,*)    ''
        write(6,*)    'rdiag_Qke',Qke
        write(6,*)    'rdiag_Sh3D',Sh3D
        write(6,*)    'rdiag_el_pbl',el_pbl
        write(6,*)    'rdiag_cldfra_bl',cldfra_bl
        write(6,*)    'rdiag_qc_bl',qc_bl
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
        write(6,*)    '; Qc tendency due to EDMF parameterization (kg/kg/s)'
        write(6,3002) ' RQCBLTEN = (/'    ,Output_edmf%RQCBLTEN(ii_write,:,jj_write)
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
        !write(6,3002) '  = (/'    ,Output_edmf%
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

3000 format (A35,2X,F8.2,',')
3001 format (A35,2X,34(F10.3,2X,','))
3002 format (A35,2X,34(E12.4,2X,','))
3003 format (A35,2X,E12.4,',')
3004 format (A35,2X,33(F10.3,2X,','),A5)

end subroutine edmf_writeout_column

!########################

subroutine convert_edmf_to_am4_array (ix, jx, kx, &
                                      Input_edmf, Output_edmf, am4_Output_edmf, rdiag, &
                                      Qke, el_pbl, cldfra_bl, qc_bl, Sh3D )

!--- input arguments
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
!------------------------------------------

!print*,'convert, Output_edmf%Qke, ix,jx,kx',size(Output_edmf%Qke,1),size(Output_edmf%Qke,2),size(Output_edmf%Qke,3)
  !ix = size(Output_edmf%Qke,1)
  !jx = size(Output_edmf%Qke,3)
  !kx = size(Output_edmf%Qke,2)
  !print*,'convert, ix,jx,kx',ix,jx,kx

  do i=1,ix
  do j=1,jx
  do k=1,kx
    kk=kx-k+1    

    !--- set am4_Output_edmf
    am4_Output_edmf%tke         (i,j,kk) = 0.5 * Output_edmf%Qke (i,k,j)
    am4_Output_edmf%Tsq         (i,j,kk) = Output_edmf%Tsq (i,k,j)
    am4_Output_edmf%Cov_thl_qt  (i,j,kk) = Output_edmf%Cov (i,k,j)
    am4_Output_edmf%udt_edmf    (i,j,kk) = Output_edmf%RUBLTEN (i,k,j) 
    am4_Output_edmf%vdt_edmf    (i,j,kk) = Output_edmf%RVBLTEN (i,k,j)
    am4_Output_edmf%tdt_edmf    (i,j,kk) = Output_edmf%RTHBLTEN (i,k,j) / Input_edmf%exner (i,k,j)
    am4_Output_edmf%qdt_edmf    (i,j,kk) = Output_edmf%RQVBLTEN (i,k,j)
    am4_Output_edmf%qidt_edmf   (i,j,kk) = Output_edmf%RQIBLTEN (i,k,j)
    am4_Output_edmf%qcdt_edmf   (i,j,kk) = Output_edmf%RQCBLTEN (i,k,j)
    am4_Output_edmf%edmf_a      (i,j,kk) = Output_edmf%edmf_a   (i,k,j)
    am4_Output_edmf%edmf_w      (i,j,kk) = Output_edmf%edmf_w   (i,k,j)
    am4_Output_edmf%edmf_qt     (i,j,kk) = Output_edmf%edmf_qt  (i,k,j)
    am4_Output_edmf%edmf_thl    (i,j,kk) = Output_edmf%edmf_thl (i,k,j)
    am4_Output_edmf%edmf_ent    (i,j,kk) = Output_edmf%edmf_ent (i,k,j)
    am4_Output_edmf%edmf_qt     (i,j,kk) = Output_edmf%edmf_qt  (i,k,j)
    !!! am4_Output_edmf% (i,j,kk) = Output_edmf% (i,k,j)

    !--- change rdiag
    Qke       (i,j,kk) = Output_edmf%Qke       (i,k,j)
    el_pbl    (i,j,kk) = Output_edmf%el_pbl    (i,k,j)
    cldfra_bl (i,j,kk) = Output_edmf%cldfra_bl (i,k,j)
    qc_bl     (i,j,kk) = Output_edmf%qc_bl     (i,k,j)
    Sh3D      (i,j,kk) = Output_edmf%Sh3D      (i,k,j)
    !!! rdiag(i,j,kk,n)      = Output_edmf%      (i,k,j)

  enddo  ! end loop of k
  enddo  ! end loop of j
  enddo  ! end loop of 1

end subroutine convert_edmf_to_am4_array

!#############################


!###################################
!  do not copy these subroutiens

subroutine error_mesg(character1, character2, fatal0)
  character character1
  character character2
  integer fatal0
end subroutine error_mesg

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
  real,    dimension(ni,nj,nhalf) :: p_half, z_half, z_half_actual
  real,    dimension(ni,nj,nfull) :: p_full, z_full, z_full_actual, uu, vv, tt, qq, ql, qi, thv, ww, omega, th, thl, qt
  !real,    dimension(nhalf) :: p_half, z_half, z_half_actual
  !real,    dimension(nfull) :: p_full, z_full, z_full_actual, uu, vv, tt, qq, thv, ww, ql, qi
  real,    dimension(ni,nj,nfull) :: diff_t, diff_m
  real,    dimension(ni,nj,nfull) :: udt_mf, vdt_mf, tdt_mf, qdt_mf, thvdt_mf, qtdt_mf, thlidt_mf
  !real,    dimension(ni,nj)   :: cov_w_thv, cov_w_qt
  real,    dimension(ni,nj)   :: z_pbl
  real,    dimension(ni,nj)   :: lat, lon, u_star, b_star, q_star, shflx, lhflx, frac_land, area, t_ref, q_ref, u_flux, v_flux
  real,    dimension(ni,nj)   :: buoy_flux, w1_thv1_surf, w1_qt1_surf

  real,    dimension(ni,nj,nfull,ntracer) :: rdiag
  real,    dimension(ni,nj,nfull) :: udt, vdt, tdt
  real,    dimension(ni,nj,nfull,tr) :: rdt 
  
  integer mm

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
  qi (1,1,:) = (/    0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.2224E-01  ,  0.2240E-01  ,  0.2094E-01  ,  0.1884E-01  ,  0.1462E-01  ,  0.6536E-02  ,  0.4606E-02  ,  0.4494E-02  ,  0.4155E-02  ,  0.3521E-02  ,  0.3257E-02  ,  0.2755E-02  ,  0.9241E-03  ,  0.2604E-04  ,  0.4037E-11  ,  0.6534E-28  ,  0.3008E-56  ,  0.7458E-73  ,  0.3929-196  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  ,  0.1204-213  ,  0.3481-215  ,  0.0000E+00  ,  0.0000E+00  ,  0.0000E+00  /)
  
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

else
  print*,'ERROR: unsupported input_profile, ', input_profile
  stop
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!elseif (input_profile == "") then
endif ! end if of input profile  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
!==============================
!==============================
! run edmf_mynn
!==============================
!==============================

  print*,'hello world!'
  print*,'input profile, ',input_profile

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

  do mm = 1,loop_times
    ! set tendencies to zeros
    udt = 0. ; vdt = 0. ; tdt = 0.; rdt = 0.   

print*,'tt',Physics_input_block%t

    print*,''
    print*,'----------------------'
    print*, 'loop times=',mm
    print*,'----------------------'
    call edmf_mynn_driver ( &
              is, ie, js, je, npz, Time_next, dt, lon, lat, frac_land, area, u_star,  &
              b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, & 
              do_edmf_mynn_diagnostic, &
              udt, vdt, tdt, rdt, rdiag(:,:,:,ntp+1:))
              !udt, vdt, tdt, rdt, rdiag)

    ! update fields
    Physics_input_block%u = Physics_input_block%u + udt(:,:,:)*dt
    Physics_input_block%v = Physics_input_block%v + vdt(:,:,:)*dt
    Physics_input_block%t = Physics_input_block%t + tdt(:,:,:)*dt
    Physics_input_block%q(:,:,:,nsphum) = Physics_input_block%q(:,:,:,1) + rdt(:,:,:,nsphum)*dt
  enddo


  area = 0.
  area = 1/area

end program test111


