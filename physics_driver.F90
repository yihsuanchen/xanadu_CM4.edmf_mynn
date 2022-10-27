module physics_driver_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="">
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!     Provides high level interfaces for calling the entire
!     FMS atmospheric physics package.
!
!    physics_driver_mod accesses the model's physics modules and
!    obtains tendencies and boundary fluxes due to the physical
!    processes that drive atmospheric time tendencies and supply 
!    boundary forcing to the surface models.
! </OVERVIEW>
! <DESCRIPTION>
!     This version of physics_driver_mod has been designed around the implicit
!     version diffusion scheme of the GCM. It requires two routines to advance
!     the model one time step into the future. These two routines
!     correspond to the down and up sweeps of the standard tridiagonal solver.
!     Radiation, Rayleigh damping, gravity wave drag, vertical diffusion of
!     momentum and tracers, and the downward pass of vertical diffusion for
!     temperature and specific humidity are performed in the down routine.
!     The up routine finishes the vertical diffusion and computes moisture
!     related terms (convection,large-scale condensation, and precipitation).
! </DESCRIPTION>
! <DIAGFIELDS>
! </DIAGFIELDS>
! <DATASET NAME="physics_driver.res">
! native format restart file
! </DATASET>
!
! <DATASET NAME="physics_driver.res.nc">
! netcdf format restart file
! </DATASET>


! <INFO>

!   <REFERENCE>            </REFERENCE>
!   <COMPILER NAME="">     </COMPILER>
!   <PRECOMP FLAG="">      </PRECOMP>
!   <LOADER FLAG="">       </LOADER>
!   <TESTPROGRAM NAME="">  </TESTPROGRAM>
!   <BUG>                  </BUG>
!   <NOTE> 
!   </NOTE>
!   <FUTURE> Deal with conservation of total energy?              </FUTURE>

! </INFO>
!   shared modules:

use time_manager_mod,        only: time_type, get_time, operator (-), &
                                   time_manager_init, operator(*)
use field_manager_mod,       only: field_manager_init, MODEL_ATMOS
use tracer_manager_mod,      only: tracer_manager_init, &
                                   get_number_tracers, &
                                   get_tracer_names, &
                                   get_tracer_index, NO_TRACER
use block_control_mod,       only: block_control_type
use atmos_tracer_driver_mod, only: atmos_tracer_driver_init,    &
                                   atmos_tracer_driver_time_vary, &
                                   atmos_tracer_driver_endts, &
                                   atmos_tracer_driver,  &
                                   atmos_tracer_driver_end
use mpp_mod,                 only: input_nml_file
use fms_mod,                 only: mpp_clock_id, mpp_clock_begin,   &
                                   mpp_clock_end, CLOCK_MODULE_DRIVER, &
                                   fms_init,  &
                                   open_namelist_file, stdlog, stdout,  &
                                   write_version_number, field_size, &
                                   file_exist, error_mesg, FATAL,   &
                                   WARNING, NOTE, check_nml_error, &
                                   close_file, mpp_pe, mpp_root_pe, &
                                   mpp_error, mpp_chksum, string
use fms_io_mod,              only: restore_state, &
                                   register_restart_field, restart_file_type, &
                                   save_restart, get_mosaic_tile_file

use diag_manager_mod,        only: register_diag_field, send_data

! shared atmospheric package modules:

use atmos_cmip_diag_mod,     only: register_cmip_diag_field_3d, &
                                   send_cmip_data_3d, &
                                   cmip_diag_id_type, &
                                   query_cmip_diag_id

!    shared radiation package modules:

use aerosol_types_mod,       only: aerosol_type, aerosol_time_vary_type

use physics_radiation_exch_mod, only: exchange_control_type, &
                                      clouds_from_moist_type, &
                                      clouds_from_moist_block_type, &
                                      cosp_from_rad_type, &
                                      cosp_from_rad_control_type, &
                                      cosp_from_rad_block_type, &
                                      radiation_flux_control_type, & 
                                      radiation_flux_block_type, & 
                                      alloc_clouds_from_moist_type, &
                                      alloc_cloud_scheme_data_type

use physics_types_mod,       only: alloc_physics_tendency_type, &
                                   physics_tendency_type, & 
                                   phys_mp_exch_type, &
                                   phys2cosp_type, precip_flux_type, &
                                   physics_tendency_block_type, &
                                   physics_type, & 
                                   physics_control_type, & 
                                   physics_input_block_type, &
                                   dealloc_physics_tendency_type

use moist_proc_utils_mod, only:    mp_removal_type, column_diag
  
use aerosol_mod,             only: aerosol_init, aerosol_driver, &
                                   aerosol_time_vary, &
                                   aerosol_endts, &
                                   aerosol_dealloc, aerosol_end
!    component modules:

use cosp_driver_mod,         only: cosp_driver_init, cosp_driver, &
                                   cosp_driver_end, cosp_driver_time_vary, &
                                   cosp_driver_endts
use  moist_processes_mod,    only: moist_processes,    &
                                   moist_processes_init,  &
                                   set_cosp_precip_sources, &
                                   define_cosp_precip_fluxes, &
                                   moist_processes_time_vary, &
                                   moist_processes_endts, &
                                   moist_processes_restart, &
                                   moist_processes_end

use vert_turb_driver_mod,    only: vert_turb_driver,  &
                                   vert_turb_driver_init,  &
                                   vert_turb_driver_end, &
                                   vert_turb_driver_restart

use vert_diff_driver_mod,    only: vert_diff_driver_down,  &
                                   vert_diff_driver_up,    &
                                   vert_diff_driver_init,  &
                                   vert_diff_driver_end,   &
                                   surf_diff_type
 
use damping_driver_mod,      only: damping_driver,      &
                                   damping_driver_init, &
                                   damping_driver_time_vary,  &
                                   damping_driver_endts, &
                                   damping_driver_end,  &
                                   damping_driver_restart

use grey_radiation_mod,       only: grey_radiation_init, grey_radiation, &
                                    grey_radiation_end

use monin_obukhov_mod,        only: monin_obukhov_init

!<-- yhc, add edmf_mynn
use edmf_mynn_mod,            only: edmf_mynn_init, edmf_mynn_end, &
                                    edmf_input_type, edmf_output_type, edmf_ls_mp_type, &
                                    edmf_mynn_driver

use     constants_mod, only: vonkarm, rdgas, rvgas, kappa, grav, &
                             cp_air, dens_h2o, hlv, hlf, hls, pi
!--> yhc, add edmf_mynn

#ifdef SCM
! Option to add SCM radiative tendencies from forcing to lw_tendency
! and radturbten

use scm_forc_mod,            only: use_scm_rad, add_scm_tdtlw, add_scm_tdtsw

#endif

!-----------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    physics_driver_mod accesses the model's physics modules and
!    obtains tendencies and boundary fluxes due to the physical
!    processes that drive atmospheric time tendencies and supply 
!    boundary forcing to the surface models.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'


!---------------------------------------------------------------------
!-------  interfaces --------

public  physics_driver_init, physics_driver_down,   &
        physics_driver_down_time_vary, physics_driver_up_time_vary, &
        physics_driver_down_endts, physics_driver_up_endts, &
        physics_driver_up, physics_driver_end, &
        do_moist_in_phys_up, get_diff_t, &
        get_radturbten, zero_radturbten, physics_driver_restart, &
        cosp_driver_init, set_cosp_precip_sources

private          &

!  called from physics_driver_down:
         check_args, &

!  called from physics_driver_init:
         physics_driver_register_restart, &

!  called from physics_driver_restart:
         physics_driver_netcdf, &

!  called from check_args:
         check_dim

interface check_dim
     module procedure check_dim_2d, check_dim_3d, check_dim_4d
end interface


!---------------------------------------------------------------------
!------- namelist ------
 
! <NAMELIST NAME="physics_driver_nml">
!  <DATA NAME="do_radiation" UNITS="" TYPE="logical" DIM="" DEFAULT=".true.     ">
!   calculating radiative fluxes and  heating rates?
!  </DATA>
!  <DATA NAME="do_clubb" UNITS="" TYPE="integer" DIM="" DEFAULT="0">
!   do_clubb > 0 implies clubb is active in some manner
!  </DATA>
!  <DATA NAME="do_cosp" UNITS="" TYPE="logical" DIM="" DEFAULT=".false.">
!   activate COSP simulator ?
!  </DATA>
!  <DATA NAME="do_modis_yim" UNITS="" TYPE="logical" DIM="" DEFAULT=".true.">
!   activate simple modis simulator ?
!  </DATA>
!  <DATA NAME="donner_meso_is_largescale" UNITS="" TYPE="logical" DIM="" DEFAULT=".true.">
!   donner meso clouds are treated as largescale (rather than convective)
!   as far as the COSP simulator is concerned ?
!  </DATA>
!  <DATA NAME="do_moist_processes" UNITS="" TYPE="logical" DIM="" DEFAULT="     .true.">
!   call moist_processes routines ?
!  </DATA>
!  <DATA NAME="tau_diff" UNITS="" TYPE="real" DIM="" DEFAULT="3600.">
!   time scale for smoothing diffusion coefficients
!  </DATA>
!  <DATA NAME="diff_min" UNITS="" TYPE="real" DIM="" DEFAULT="1.e-3">
!   minimum value of a diffusion coefficient beneath which the 
!   coefficient is reset to zero
!  </DATA>
!  <DATA NAME="diffusion_smooth" UNITS="" TYPE="logical" DIM="" DEFAULT=".t     rue.">
!   diffusion coefficients should be smoothed in time?
!  </DATA>
!  <DATA NAME="do_grey_radiation" UNITS="" TYPE="logical" DIM="" DEFAULT=".false.">
!   do grey radiation scheme?
! rif:(09/10/09) In Grey radiation we are computing just the total   
! SW radiation. We need to divide it into 4 components
! to go through the Coupler and Ice modules. Sum[R(i)*SW] = SW  
!  </DATA>
!  <DATA NAME="R1" UNITS="" TYPE="real" DIM="" DEFAULT="0.25">
!   component  number 1 of SW radiation with grey radiation scheme
!  </DATA>
!  <DATA NAME="R2" UNITS="" TYPE="real" DIM="" DEFAULT="0.25">
!   component  number 1 of SW radiation with grey radiation scheme
!  </DATA>
!  <DATA NAME="R3" UNITS="" TYPE="real" DIM="" DEFAULT="0.25">
!   component  number 1 of SW radiation with grey radiation scheme
!  </DATA>
!  <DATA NAME="R4" UNITS="" TYPE="real" DIM="" DEFAULT="0.25">
!   component  number 1 of SW radiation with grey radiation scheme
!  </DATA>
!  <DATA NAME="override_aerosols_cloud" UNITS="" TYPE="logical" DIM="" DEFA     ULT=".false.">
!   use offline aerosols for cloud calculation
!   (via data_override in aerosol_driver)?
!  </DATA>
!  <DATA NAME="l_host_applies_sfc_fluxes" UNITS="" TYPE="logical" DIM="" DEFAULT=".true.">
!   applying surface fluxes in host-model ?
!  </DATA>
!  <DATA NAME="qmin" UNITS="kg h2o/kg air" TYPE="real"  DEFAULT="1.E-10">
!   minimum permissible value of cloud liquid, cloud ice, saturated volume
!   fraction, or rain and snow areas.
!  NOTE: qmin should be chosen such that the range of {qmin, max(qa,ql,qi)}
!   is resolved by the precision of the numbers used.
!  </DATA>
!  <DATA NAME="N_land" UNITS="1/(m*m*m)" TYPE="real" DEFAULT="250.E+06">
!   assumed number of cloud drops per unit volume in liquid clouds
!   over land when droplet number is not prdicted.
!  </DATA>
!  <DATA NAME="N_ocean" UNITS="1/(m*m*m)" TYPE="real" DEFAULT="100.E+06">
!   assumed number of cloud drops per unit volume in liquid clouds
!   over ocean when droplet number is not predicted.
!  </DATA>
!  <DATA NAME="do_liq_num" UNITS="" TYPE="logical"  DEFAULT=".false.">
!   the prognostic droplet number option is activated ?
!  </DATA>
!  <DATA NAME="do_ice_num" UNITS="" TYPE="logical"  DEFAULT=".false.">
!   the prognostic ice particle number option is activated ?
!  </DATA>
!  <DATA NAME="qcvar" UNITS="" TYPE="real"  DEFAULT="1.0">
!    1 / relative variance of sub-grid cloud water distribution
!    see morrison and gettelman, 2007, J. Climate for details
!  </DATA>
!  <DATA NAME="overlap" UNITS="" TYPE="integer"  DEFAULT="2">
!      overlap        integer variable indicating which overlap 
!                     assumption to use:
!                     overlap = 1. means condensate in adjacent levels 
!                                  is treated as part of the same cloud
!                                  i.e. maximum-random overlap
!                     overlap = 2. means condensate in adjacent levels 
!                                  is treated as different clouds
!                                  i.e. random overlap
!  </DATA>
!  <DATA NAME="N_min"  UNITS="1/(m*m*m) " TYPE="real" DEFAULT="1.0E6">
!   minimum number of droplets allowed in a grid box when predicted
!   droplet number code is activated
!  </DATA>
!  <DATA NAME="min_diam_ice" UNITS="m" TYPE="real" DEFAULT="10.E-6">
!   minimum size of ice crystals allowed in microphysics and radiation
!   calculations
!  </DATA>
!  <DATA NAME="dcs" UNITS="m" TYPE="real" DEFAULT="200.E-6">
!   ice crystal size at which autoconversion to falling ice occurs;
!   used in radiation and microphysics calculations
!  </DATA>
!  <DATA NAME="min_diam_drop" UNITS="m" TYPE="real" DEFAULT="2.E-6">
!   minimum size of droplets allowed in microphysics and radiation
!   calculations
!  </DATA>
!  <DATA NAME="max_diam_drop" UNITS="m" TYPE="real" DEFAULT="50.E-6">
!   maximum size of droplets allowed in microphysics and radiation
!   calculations
!  </DATA>
!  <DATA NAME="use_tau" UNITS="" TYPE="logical"  DEFAULT=".false.">
!    switch to determine whether current time level (tau)
!    will be used or else future time level (tau+1).
!    if use_tau = true then the input values for t,q, and r
!    are used; if use_tau = false then input values
!    tm+tdt*dt, etc. are used in moist_processes and vert_turb_driver.
!  </DATA>
!  <DATA NAME="cosp_frequency" UNITS="sec" TYPE="real" DEFAULT="10800.">
!   frequency at which the COSP simulator is to be called
!  </DATA>
! </NAMELIST>
 

logical :: do_radiation = .true.
integer :: do_clubb = 0        
logical :: do_cosp = .false.   
logical :: do_modis_yim = .true.
logical :: donner_meso_is_largescale = .true.
logical :: do_moist_processes = .true.
real    :: tau_diff = 3600.    
real    :: diff_min = 1.e-3   
logical :: diffusion_smooth = .true.
logical :: do_grey_radiation = .false.
real    :: R1 = 0.25
real    :: R2 = 0.25
real    :: R3 = 0.25
real    :: R4 = 0.25
logical :: override_aerosols_cloud = .false.
logical :: l_host_applies_sfc_fluxes = .true.
real    :: qmin = 1.0e-10
real    :: N_land = 3.e8
real    :: N_ocean = 1.e8
logical :: do_liq_num = .true.
logical :: do_ice_num = .false.
real    :: qcvar = 1.
integer :: overlap = 2
real    :: N_min = 1.E6
real    :: min_diam_ice    = 10.e-6
real    :: min_diam_drop    = 2.e-6
real    :: max_diam_drop    = 50.e-6
real    :: dcs    =  200.e-6
logical :: use_tau = .false.
real    :: cosp_frequency = 10800.

!<-- yhc, add edmf_mynn nml
logical :: do_edmf_mynn = .false.             ! switch to turn on/off edmf_mynn scheme
character*5 :: do_edmf_mynn_in_physics = "up"     ! where to call edmf_mynn. "up" in physics_driver_up; "down" in physics_driver_down
logical :: do_edmf_mynn_diagnostic = .false.  ! .true.  : diagnostic edmf_mynn, no update on tendencies
                                              ! .false. : interactive edmf_mynn, update on tendencies

logical :: do_tracers_in_edmf_mynn = .false.  ! do tracer diffusion in edmf_mynn
                                              ! .true.  : tracer diffusion is handled by edmf_mynn
                                              ! .false. : tracer diffusion is handled by vert_diff_driver_down, using 
                                              !           ED diffusion coefficients from edmf_mynn
logical :: do_return_edmf_mynn_diff_only = .false. ! .true.,  return edmf_mynn diffusion coefficients 
                                                   !          and let vert_diff do the diffusion rather than in edmf_mynn
                                                   ! .false., tendencies are updated in edmf_mynn

integer :: do_tracers_selective = 0           ! determine which tracer would be processed by vert_diff_down
                                              !   0: vert_diff would process all tracers, such as qa, ql, qn, dust, sea salt, etc
                                              !   1: vert_diff would process all tracers except qn and qni
                                              !   2: vert_diff would NOT process any tracers
                                              !   3: vert_diff would ONLY process qn and qni        
                                              !   4: process all tracers except cloud fields: qa,ql,qi,nqn,nqi
                                              !   5: process all tracers except cloud fraction and specific humidity qa,ql,qi. Cloud number are processed by vert_diff

logical :: do_edmf_mynn_diffusion_smooth = .false.  ! smooth the diffusion coefficients following AM4 "diffusion_smooth"
logical :: do_bomex_radf                 = .false.  ! add BOMEX 2K/day radiative cooling in physics_driver

integer :: ii_write = -999                    ! i index for column written out. Set to 0 if you want to write out in SCM
integer :: jj_write = -999                    ! j index for column written out. Set to 0 if you want to write out in SCM
real    :: lat_write = -999.99                ! latitude  (radian) for column written out
real    :: lon_write = -999.99                ! longitude (radian) for column written out
logical :: do_stop_run = .false.              ! whether to stop the simulation
logical :: do_writeout_column_nml = .false.   ! switch to control whether writing out the column

logical :: do_writeout_column = .false.       ! local variable
real    :: lat_lower, lat_upper, lon_lower, lon_upper
real    :: tdt_max     = 500. ! K/day
logical :: do_limit_tdt = .false.
real    :: tdt_limit =   200. ! K/day


!--- namelist for compute_vert_adv_tend_offline
logical :: do_vert_adv_tend_offline = .false.
integer :: input_profile = -1
integer :: temp_vert_advec_scheme = 4
integer :: tracer_vert_advec_scheme =4
integer :: printout_vert_adv_tend = 1
!--> yhc, add edmf_mynn nml


namelist / physics_driver_nml / do_radiation, do_clubb,  do_cosp, &
                                do_modis_yim, donner_meso_is_largescale, &
                                do_moist_processes, tau_diff,      &
                                diff_min, diffusion_smooth, &
                                do_grey_radiation, R1, R2, R3, R4,  &
                                override_aerosols_cloud,    &
                                l_host_applies_sfc_fluxes, &
                                qmin, N_land, N_ocean, do_liq_num,  &
                                do_ice_num, qcvar, overlap, N_min, &
                                min_diam_ice, dcs, min_diam_drop, &
                                do_edmf_mynn, do_edmf_mynn_diagnostic, do_tracers_in_edmf_mynn, do_tracers_selective, do_edmf_mynn_diffusion_smooth, do_return_edmf_mynn_diff_only, do_edmf_mynn_in_physics, & ! yhc add
                                do_stop_run, do_writeout_column_nml, do_bomex_radf, ii_write, jj_write, lat_write, lon_write, & ! yhc add
                                tdt_max, do_limit_tdt, tdt_limit, &  ! yhc add
                                do_vert_adv_tend_offline, input_profile, printout_vert_adv_tend, temp_vert_advec_scheme, tracer_vert_advec_scheme, & ! yhc add
                                max_diam_drop, use_tau, cosp_frequency


!---------------------------------------------------------------------
!------- public data ------
! <DATA NAME="surf_diff_type" UNITS="" TYPE="surf_diff_type" DIM="" DEFAULT="">
! Defined in vert_diff_driver_mod, republished here. See vert_diff_mod for details.
! </DATA>

public  surf_diff_type   ! defined in  vert_diff_driver_mod, republished
                         ! here

!---------------------------------------------------------------------
!------- private data ------

!--------------------------------------------------------------------
! list of restart versions readable by this module:
!
! version 1: initial implementation 1/2003, contains diffusion coef-
!            ficient contribution from cu_mo_trans_mod. This variable
!            is generated in physics_driver_up (moist_processes) and
!            used on the next step in vert_diff_down, necessitating
!            its storage.
!
! version 2: adds pbltop as generated in vert_turb_driver_mod. This 
!            variable is then used on the next timestep by topo_drag
!            (called from damping_driver_mod), necessitating its 
!            storage.
!
! version 3: adds the diffusion coefficients which are passed to 
!            vert_diff_driver.  These diffusion are saved should
!            smoothing of vertical diffusion coefficients be turned
!            on.
!
! version 4: adds a logical variable, convect, which indicates whether
!            or not the grid column is convecting. This diagnostic is
!            needed by the entrain_module in vert_turb_driver.
!
! version 5: adds radturbten when strat_cloud_mod is active, adds 
!            lw_tendency when edt_mod or entrain_mod is active.
!
! version 6: adds donner cell and meso cloud variables when donner_deep
!            is activated.

! version 7: adds shallow convection cloud variables when uw_conv
!            is activated.

! version 8: adds lsc cloud props for radiation. only readable when in
!            netcdf mode.


!---------------------------------------------------------------------
integer, dimension(8) :: restart_versions = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)

!--------------------------------------------------------------------
!    the following allocatable arrays are either used to hold physics 
!    data between timesteps when required, or hold physics data between
!    physics_down and physics_up.
!  
!    diff_cu_mo     contains contribution to difusion coefficient
!                   coming from cu_mo_trans_mod (called from 
!                   moist_processes in physics_driver_up) and then used 
!                   as input on the next time step to vert_diff_down 
!                   called in physics_driver_down.
!    diff_t         vertical diffusion coefficient for temperature
!                   which optionally may be time smoothed, meaning
!                   values must be saved between steps
!    diff_m         vertical diffusion coefficient for momentum
!                   which optionally may be time smoothed, meaning
!                   values must be saved between steps
!    radturbten     the sum of the radiational and turbulent heating,
!                   generated in both physics_driver_down (radiation)
!                   and physics_driver_up (turbulence) and then used
!                   in moist_processes
!    pbltop         top of boundary layer obtained from vert_turb_driver
!                   and then used on the next timestep in topo_drag_mod
!                   called from damping_driver_down        
!    cush
!    cbmf
!    convect        flag indicating whether convection is occurring in
!                   a grid column. generated in physics_driver_up and
!                   then used in vert_turb_driver called from 
!                   physics_driver_down on the next step.
!    temp_last
!    q_last
!    diff_t_clubb
!----------------------------------------------------------------------
real,    dimension(:,:,:), allocatable,target :: diff_cu_mo, diff_t, diff_m
real,    dimension(:,:,:), allocatable,target :: radturbten
real,    dimension(:,:)  , allocatable,target :: pbltop, cush, cbmf
real,    dimension(:,:)  , allocatable,target :: hmint, cgust, tke 
real,    dimension(:,:)  , allocatable,target :: pblhto, rkmo, taudpo
logical, dimension(:,:)  , allocatable,target :: convect
integer, dimension(:,:,:), allocatable,target :: exist_shconv, exist_dpconv
real,    dimension(:,:,:), allocatable,target :: pblht_prev, hlsrc_prev, &
                                                 qtsrc_prev, cape_prev,  &
                                                 cin_prev, tke_prev !miz
real,    dimension(:,:,:), allocatable,target ::  diff_t_clubb

real,    dimension(:,:,:), allocatable        :: temp_last, q_last

!<-- yhc
real,    dimension(:,:,:), allocatable,target :: &  ! yhc, the description of these variables is in edmf_ls_mp_type, edmf_mynn_mod
  qadt_edmf, qldt_edmf, qidt_edmf, &
  diff_t_edmf, diff_m_edmf, &
  dqa_edmf,  dql_edmf, dqi_edmf, &
  edmf_mc_full, edmf_mc_half, edmf_moist_area, edmf_dry_area, edmf_moist_humidity, edmf_dry_humidity

real,    dimension(:,:,:,:), allocatable,target :: & 
  rdt_before_vdiff_down

real,    dimension(:,:,:), allocatable,target :: & 
  tdt_before_vdiff_down

integer, dimension(:,:), allocatable,target :: kpbl_edmf  

integer,                   allocatable,target :: option_edmf2ls_mp  ! yhc
!type(edmf_ls_mp_type)                         :: edmf2ls_mp
!--> yhc

!--- for netcdf restart
type(restart_file_type), pointer, save :: Phy_restart => NULL()
type(restart_file_type), pointer, save :: Til_restart => NULL()
logical                                :: in_different_file = .false.
integer                                :: vers
integer                                :: now_doing_strat = 0
integer                                :: now_doing_entrain = 0
integer                                :: now_doing_edt = 0
real, allocatable                      :: r_convect(:,:)

type(aerosol_time_vary_type)           :: Aerosol_cld

!---------------------------------------------------------------------
!    internal timing clock variables:
!---------------------------------------------------------------------
integer :: damping_clock, turb_clock,   &
           tracer_clock, diff_up_clock, diff_down_clock, &
           edmf_mynn_clock, &  ! yhc add 
           moist_processes_clock, cosp_clock

!--------------------------------------------------------------------
!    miscellaneous control variables:
!---------------------------------------------------------------------
logical   :: do_check_args = .true.   ! argument dimensions should 
                                      ! be checked ?
logical   :: module_is_initialized = .false.
                                      ! module has been initialized ?
logical   :: doing_edt                ! edt_mod has been activated ?
logical   :: doing_entrain            ! entrain_mod has been activated ?
logical   :: doing_uw_conv            ! uw_conv shallow cu mod has been 
                                      ! activated ?
logical   :: doing_liq_num = .false.  ! Prognostic cloud droplet number has 
                                      ! been activated?
integer   :: nt                       ! total no. of tracers
integer   :: ntp                      ! total no. of prognostic tracers
!integer   :: ncol                     ! number of stochastic columns
integer   ::  nsphum                  ! index for specific humidity tracer
integer   :: nqa, nql, nqi, nqn, nqni ! yhc 

logical   :: step_to_call_cosp
logical   :: doing_prog_clouds
real      :: rad_time_step

!---------------------------------------------------------------------
!---------------------------------------------------------------------

character(len=4)     :: mod_name = 'phys'
character(len=32)    :: tracer_units, tracer_name
character(len=128)   :: diaglname
real                 :: missing_value = -999.

integer                            :: id_tdt_phys,         &
                                      id_tdt_phys_vdif_dn, &
                                      id_tdt_phys_vdif_up, &
                                      id_tdt_phys_turb,    &
                                      id_tdt_phys_moist

integer, dimension(:), allocatable :: id_tracer_phys,         &
                                      id_tracer_phys_vdif_dn, &
                                      id_tracer_phys_vdif_up, &
                                      id_tracer_phys_turb,    &
                                      id_tracer_phys_moist

type(cmip_diag_id_type) :: ID_tntmp, ID_tnhusmp, &
                           ID_pfull, ID_phalf, ID_zg !, ID_zfull

type (clouds_from_moist_block_type) :: Restart

type(precip_flux_type)              :: Precip_flux

!<-- yhc
integer :: id_diff_t_vdif, id_diff_m_vdif, id_num_updraft, id_qldt_vdif, id_qadt_vdif, id_qidt_vdif, id_qdt_vdif_test, id_tdt_vdif_test, id_tdt_radturb, id_tt_phys, id_qq_phys

integer :: id_LWP0
!--> yhc
                            contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="physics_driver_init">
!  <OVERVIEW>
!    physics_driver_init is the constructor for physics_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_init is the constructor for physics_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_init (Time, lonb, latb, axes, pref, &
!                             trs, Surf_diff, phalf )
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="pref" TYPE="real">
!   reference prssure profiles
!  </IN>
!  <IN NAME="latb" TYPE="real">
!   array of model latitudes at cell corners [radians]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   array of model longitudes at cell corners [radians]
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!   axis indices, (/x,y,pf,ph/)
!                (returned from diag axis manager)
!  </IN>
!  <INOUT NAME="trs" TYPE="real">
!   atmospheric tracer fields
!  </INOUT>
!  <INOUT NAME="Surf_diff" TYPE="surf_diff_type">
!   surface diffusion derived type
!  </INOUT>
!  <IN NAME="phalf" TYPE="real">
!   pressure at model interface levels
!  </IN>
! <ERROR MSG="physics_driver_init must be called first" STATUS="FATAL">
! </ERROR>
! </SUBROUTINE>
!
subroutine physics_driver_init (Time, lonb, latb, lon, lat, axes, &
                                Surf_diff, Exch_ctrl, Atm_block,   &
                                Moist_clouds, Physics, Physics_tendency, &
                                diffm, difft)

!---------------------------------------------------------------------
!    physics_driver_init is the constructor for physics_driver_mod.
!---------------------------------------------------------------------

type(time_type),              intent(in)    :: Time
real,    dimension(:,:),      intent(in)    :: lonb, latb
real,    dimension(:,:),      intent(in)    :: lon, lat
integer, dimension(4),        intent(in)    :: axes
type(surf_diff_type),         intent(inout) :: Surf_diff
type (exchange_control_type), intent(inout) :: Exch_ctrl
type (block_control_type),    intent(in)    :: Atm_block
type(clouds_from_moist_type), intent(inout) :: Moist_clouds(:)
type(physics_type),           intent(inout) :: Physics
type(physics_tendency_type),  intent(inout) :: Physics_tendency
real,    dimension(:,:,:),    intent(out),  optional :: diffm, difft

!---------------------------------------------------------------------
!  intent(in) variables:
!     Time       current time (time_type)
!     lonb       longitude of the grid box corners [ radians ]
!     latb       latitude of the grid box corners [ radians ]
!     lon
!     lat
!     axes       axis indices, (/x,y,pf,ph/)
!                (returned from diag axis manager)
!
!   intent(inout) variables:
!     Surf_diff  surface diffusion derived type variable
!
!   intent(out), optional variables:
!
!     diffm
!     difft
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      integer :: nb, ibs, ibe, jbs, jbe
      real, dimension (size(lonb,1)-1, size(latb,2)-1) :: sgsmtn
      integer          ::  id, jd, kd, n, k, nc
      integer          ::  ierr, io, unit, logunit, outunit
      integer          ::  ndum

      integer          ::  moist_processes_init_clock, damping_init_clock,&
                           turb_init_clock, diff_init_clock, &
                           aerosol_init_clock, &
                           grey_radiation_init_clock , &
                           edmf_mynn_init_clock, &  ! yhc add
                           tracer_init_clock
      real, dimension(:,:,:),   allocatable :: phalf
      real, dimension(:,:,:,:), allocatable :: trs

!---------------------------------------------------------------------
!  local variables:
!
!       sgsmtn        sgs orography obtained from mg_drag_mod;
!                     appears to not be currently used
!       id,jd,kd      model dimensions on the processor  
!       n             loop index
!       ierr          error code
!       io            io status returned from an io call
!       unit          unit number used for an i/ operation
!       logunit
!       outunit
!       ndum          dummy argument
!       x_clock_init  clock for timing the initialization of process x
!                     where x is moist_processes, damping, turb, diff,
!                     aerosol, grey_radiation, tracer 
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, return.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that the modules used by this module that are not called 
!    later in this subroutine have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call time_manager_init
      call tracer_manager_init
      call field_manager_init (ndum)
 
!--------------------------------------------------------------------
!    read namelist.
!--------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=physics_driver_nml, iostat=io)
      ierr = check_nml_error(io,"physics_driver_nml")
#else
      if ( file_exist('input.nml')) then
        unit = open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=physics_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'physics_driver_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!--------------------------------------------------------------------
!    consistency checks for namelist options
!--------------------------------------------------------------------
      if(do_radiation .and. do_grey_radiation) & 
        call error_mesg('physics_driver_init','do_radiation and do_grey_radiation cannot both be .true.',FATAL)
      if (do_cosp .and. .not. do_radiation) &
        call error_mesg('physics_driver_init',  &
            'do_radiation must be .true. if do_cosp is .true.',FATAL)
      if (do_cosp .and. .not. do_moist_processes) &
         call error_mesg('physics_driver_init',  &
         'do_moist_processes must be .true. if do_cosp is .true.', FATAL)

!<-- yhc check edmf_mynn
      if (do_edmf_mynn_diagnostic) then
        if (.not.do_edmf_mynn) then
           call error_mesg('physics_driver_init',  &
           'When do_edmf_mynn_diagnostic is .true., do_edmf_mynn must be true', FATAL)
        endif
      endif

      if (do_tracers_in_edmf_mynn) then
        if (.not.do_edmf_mynn) then
           call error_mesg('physics_driver_init',  &
           'When do_tracers_in_edmf_mynn .true., do_edmf_mynn must be true', FATAL)
        endif
      endif

      if ( trim(do_edmf_mynn_in_physics) /= 'up' .and. &
           trim(do_edmf_mynn_in_physics) /= 'down' )       then
        call error_mesg('physics_driver_init',  &
                        'do_edmf_mynn_in_physics must be either [up] or [down]', FATAL)
      endif
!--> yhc check edmf_mynn

!--------------------------------------------------------------------
!    write version number and namelist to log file.
!--------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
               write(logunit, nml=physics_driver_nml)
 
!---------------------------------------------------------------------
!    define the model dimensions on the local processor (id, jd, kd). 
!    retrieve the total number of tracers (nt) and prognostic 
!    tracers (ntp). Save the number of prognostic tracers in a 
!    physics_control_type for use in other modules.
!---------------------------------------------------------------------
      id = size(lonb,1)-1 
      jd = size(latb,2)-1 
      kd = Atm_block%npz
      call get_number_tracers (MODEL_ATMOS, num_tracers=nt, &
                               num_prog=ntp)
      Physics%control%num_prog_tracers = ntp

!---------------------------------------------------------------------
      cosp_clock       =       &
                mpp_clock_id( '   Physics_down: COSP',    &
                   grain=CLOCK_MODULE_DRIVER )
      damping_clock         =     &
                mpp_clock_id( '   Physics_down: Damping',    &
                  grain=CLOCK_MODULE_DRIVER )
      turb_clock            =      &
                mpp_clock_id( '   Physics_down: Vert. Turb.', &
                  grain=CLOCK_MODULE_DRIVER )
      tracer_clock          =      &
                mpp_clock_id( '   Physics_down: Tracer',    &
                 grain=CLOCK_MODULE_DRIVER )
      diff_down_clock       =     &
                mpp_clock_id( '   Physics_down: Vert. Diff.',   &
                 grain=CLOCK_MODULE_DRIVER )
      diff_up_clock         =     &
                mpp_clock_id( '   Physics_up: Vert. Diff.',     &
                grain=CLOCK_MODULE_DRIVER )
      moist_processes_clock =      &
                mpp_clock_id( '   Physics_up: Moist Processes', &
                grain=CLOCK_MODULE_DRIVER )

      moist_processes_init_clock =      &
        mpp_clock_id( '   Physics_driver_init: Moist Processes: Initialization', &
                grain=CLOCK_MODULE_DRIVER )
      damping_init_clock         =     &
        mpp_clock_id( '   Physics_driver_init: Damping: Initialization',    &
                  grain=CLOCK_MODULE_DRIVER )
      turb_init_clock            =      &
        mpp_clock_id( '   Physics_driver_init: Vert. Turb.: Initialization', &
                  grain=CLOCK_MODULE_DRIVER )
      diff_init_clock       =     &
        mpp_clock_id( '   Physics_driver_init: Vert. Diff.: Initialization',   &
                 grain=CLOCK_MODULE_DRIVER )
      aerosol_init_clock       =       &
        mpp_clock_id( '   Physics_driver_init: Aerosol: Initialization', &
                       grain=CLOCK_MODULE_DRIVER )
      grey_radiation_init_clock       =       &
        mpp_clock_id( '   Physics_driver_init: Grey Radiation: Initialization', &
                       grain=CLOCK_MODULE_DRIVER )
      tracer_init_clock          =      &
        mpp_clock_id( '   Physics_driver_init: Tracer: Initialization',    &
                 grain=CLOCK_MODULE_DRIVER )

!<-- yhc add 
      edmf_mynn_clock            =      &
                mpp_clock_id( '   Physics_down/up: EDMF-MYNN', &
                  grain=CLOCK_MODULE_DRIVER )
      edmf_mynn_init_clock            =      &
        mpp_clock_id( '   Physics_driver_init: EDMF-MYNN: Initialization', &
                  grain=CLOCK_MODULE_DRIVER )
!--> yhc add 


!-----------------------------------------------------------------------
!     dummy checks
!-----------------------------------------------------------------------
      if (do_ice_num .and. .not. do_liq_num) then
        call error_mesg ('physics_driver_mod',  &
           'do_ice_num can only be .true. if do_liq_num is .true.', FATAL)
      endif

!------------------------------------------------------------------------
!   place some control variables needed in multiple physics modules
!   into the Physics%control derived type for easy movement.
!------------------------------------------------------------------------
      Physics%control%use_tau = use_tau
      Physics%control%l_host_applies_sfc_fluxes = l_host_applies_sfc_fluxes
      Physics%control%nsphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
      Physics%control%nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
      Physics%control%nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
      Physics%control%nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
      Physics%control%nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
      Physics%control%nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )
      Physics%control%nqr = get_tracer_index (MODEL_ATMOS, 'rainwat')
      Physics%control%nqs = get_tracer_index (MODEL_ATMOS, 'snowwat')
      Physics%control%nqg = get_tracer_index (MODEL_ATMOS, 'graupel')

!-----------------------------------------------------------------------
!   allocate a logical array to define whether a tracer is a cloud tracer
!   (one of those defined above), or not.
!-----------------------------------------------------------------------
      allocate (Physics%control%cloud_tracer   &
                                   (Physics%control%num_prog_tracers))
      Physics%control%cloud_tracer = .FALSE. 

      if (Physics%control%nsphum /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nsphum) = .TRUE.
      endif
      if (Physics%control%nql    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nql   ) = .TRUE.
      endif
      if (Physics%control%nqi    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqi   ) = .TRUE.
      endif
      if (Physics%control%nqa    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqa   ) = .TRUE.
      endif
      if (Physics%control%nqn    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqn   ) = .TRUE.
      endif
      if (Physics%control%nqni   /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqni  ) = .TRUE.
      endif
      if (Physics%control%nqr    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqr   ) = .TRUE.
      endif
      if (Physics%control%nqs    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqs   ) = .TRUE.
      endif
      if (Physics%control%nqg    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqg   ) = .TRUE.
      endif

!----------------------------------------------------------------------
!   define logical variable indicating whether prognostic clouds (using 
!   tracer fields) are active. 
!----------------------------------------------------------------------
      if (min(Physics%control%nql, Physics%control%nqi,   &
                                            Physics%control%nqa) > 0) then
        doing_prog_clouds = .true.
      else
        doing_prog_clouds = .false.
      endif
      
!----------------------------------------------------------------------
!do some dummy checks on the tracer indices.
!----------------------------------------------------------------------
      if (doing_prog_clouds) then
        if (min(Physics%control%nql, Physics%control%nqi,    &
                                           Physics%control%nqa) <= 0) &
          call error_mesg ('physics_driver_init', &
                        'stratiform cloud tracer(s) not found', FATAL)
        if (Physics%control%nql == Physics%control%nqi .or.   &
            Physics%control%nqa == Physics%control%nqi .or.  &
            Physics%control%nql == Physics%control%nqa)    & 
          call error_mesg ('physics_driver_init',  &
          'tracers indices cannot be the same (i.e., nql=nqi=nqa).', FATAL)
        if (mpp_pe() == mpp_root_pe()) &
            write (logunit,'(a,3i4)')   &
              'Stratiform cloud tracer indices: nql,nqi,nqa =',  &
                  Physics%control%nql, Physics%control%nqi,   &
                                                    Physics%control%nqa
                                  
        if (Physics%control%nqn == NO_TRACER .and.   &
                              Exch_ctrl%do_liq_num ) &
           call error_mesg ('physics_driver_init', &
                     'prognostic droplet number scheme requested but&
                                           &  tracer not found', FATAL)
        if (Physics%control%nqni == NO_TRACER .and.   &
                                  Exch_ctrl%do_ice_num ) &
                call error_mesg ('physics_driver_init', &
                     'prognostic ice number scheme requested but &
                                              &tracer not found', FATAL)
      endif

!----------------------------------------------------------------------
!    place the physics_driver_nml variables that are needed by both the 
!    moist_processes and radiation codes into Exch_ctrl.
!----------------------------------------------------------------------
      Exch_ctrl%cosp_frequency = cosp_frequency
      Exch_ctrl%N_min = N_min
      Exch_ctrl%qmin = qmin
      Exch_ctrl%qcvar   = qcvar
      Exch_ctrl%N_land = N_land
      Exch_ctrl%N_ocean = N_ocean
      if (overlap.ne.1 .and. overlap.ne.2) &
          call error_mesg  ('physics_driver_init',&
                               'overlap must be either 1 or 2 ', FATAL)
      Exch_ctrl%overlap = overlap
      Exch_ctrl%do_liq_num = do_liq_num
      Exch_ctrl%do_ice_num = do_ice_num

      Exch_ctrl%min_diam_ice = min_diam_ice
      Exch_ctrl%dcs          = dcs            
      Exch_ctrl%min_diam_drop = min_diam_drop
      Exch_ctrl%max_diam_drop = max_diam_drop

      Exch_ctrl%do_clubb = do_clubb
      Exch_ctrl%do_cosp = do_cosp
      Exch_ctrl%donner_meso_is_largescale =  donner_meso_is_largescale
      Exch_ctrl%do_modis_yim              =  do_modis_yim
      Exch_ctrl%doing_prog_clouds = doing_prog_clouds

!-----------------------------------------------------------------------
!--- allocate Physics_tendency to hold the physics tendencies
!-----------------------------------------------------------------------
      call alloc_physics_tendency_type (Physics_tendency, Atm_block)

!--- define trs and p_half on the full domain 
      allocate (trs(id,jd,kd,nt), phalf(id,jd,kd+1))
      do nb = 1, Atm_block%nblks
        ibs = Atm_block%ibs(nb)-Atm_block%isc+1
        ibe = Atm_block%ibe(nb)-Atm_block%isc+1
        jbs = Atm_block%jbs(nb)-Atm_block%jsc+1
        jbe = Atm_block%jbe(nb)-Atm_block%jsc+1
        trs(ibs:ibe,jbs:jbe,:,1:ntp)    = Physics%block(nb)%q
        trs(ibs:ibe,jbs:jbe,:,ntp+1:nt) = Physics%block(nb)%tmp_4d
        phalf(ibs:ibe,jbs:jbe,:)        = Physics%block(nb)%p_half
!--- the 'temp' variable inside of Physics is no longer needed - deallocate it
        deallocate(Physics%block(nb)%tmp_4d)
      enddo

!-----------------------------------------------------------------------
!---------- initialize physics -------

      if (do_moist_processes) then
        call mpp_clock_begin ( moist_processes_init_clock )
        call moist_processes_init (id, jd, kd, lonb, latb, lon, lat,  &
                                   phalf, Physics%glbl_qty%pref(:,1),&
                                   axes, Time, Physics%control, Exch_ctrl) 

        call mpp_clock_end ( moist_processes_init_clock )
      else
        diff_cu_mo = 0.0
        convect = .false.
      endif
     
!-----------------------------------------------------------------------
!    initialize damping_driver_mod.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( damping_init_clock )
      call damping_driver_init (lonb, latb, Physics%glbl_qty%pref(:,1), &
                                axes, Time, sgsmtn)
      call mpp_clock_end ( damping_init_clock )

!-----------------------------------------------------------------------
!    initialize vert_turb_driver_mod.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( turb_init_clock )
      call vert_turb_driver_init (lonb, latb, id, jd, kd, axes, Time, &
                                  Exch_ctrl, Physics%control,  &
                                  doing_edt, doing_entrain, do_clubb)
      call mpp_clock_end ( turb_init_clock )

!-----------------------------------------------------------------------
!    initialize vert_diff_driver_mod.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( diff_init_clock )
      call vert_diff_driver_init (Surf_diff, id, jd, kd, axes, Time,   &
                                  do_clubb )
      call mpp_clock_end ( diff_init_clock )

      if (do_moist_processes) then
!-----------------------------------------------------------------------
!    initialize aerosol_mod     
!-----------------------------------------------------------------------
        call mpp_clock_begin ( aerosol_init_clock )
        call aerosol_init (lonb, latb, Aerosol_cld)
        call mpp_clock_end ( aerosol_init_clock )
      endif ! do_moist_processes

!----------------------------------------------------------------------
!    if grey_radiation is active, initialize that module.
!----------------------------------------------------------------------
      if(do_grey_radiation) then
         call mpp_clock_begin ( grey_radiation_init_clock )
         call grey_radiation_init(axes, Time) 
         call mpp_clock_end ( grey_radiation_init_clock )
      endif
        
!-----------------------------------------------------------------------
!    initialize atmos_tracer_driver_mod.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( tracer_init_clock )
      call atmos_tracer_driver_init (lonb, latb, trs, axes, Time, phalf)
      call mpp_clock_end ( tracer_init_clock )

!<-- yhc, initialize edmf_mynn module
!-----------------------------------------------------------------------
!    initialize edmf_mynn_mod.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( edmf_mynn_init_clock )
      call edmf_mynn_init  (lonb, latb, axes, Time, id, jd, kd)
      call mpp_clock_end ( edmf_mynn_init_clock )

      nqa = Physics%control%nqa      
      nql = Physics%control%nql 
      nqi = Physics%control%nqi      
      nqn = Physics%control%nqn     
      nqni = Physics%control%nqni   ! nqni is not supported by SCM yet    

!--> yhc, initialize edmf_mynn module

!---------------------------------------------------------------------
!    allocate space for the module variables and initialize them.
!---------------------------------------------------------------------
      allocate ( diff_t     (id, jd, kd) ) ; diff_t = 0.0
      allocate ( diff_m     (id, jd, kd) ) ; diff_m = 0.0
      allocate ( diff_cu_mo (id, jd, kd) ) ; diff_cu_mo = 0.0
      allocate ( pbltop     (id, jd) )     ; pbltop     = -999.0
      allocate ( cush       (id, jd) )     ; cush=-1. !miz
      allocate ( cbmf       (id, jd) )     ; cbmf=0.0 !miz
      allocate ( hmint      (id, jd) )     ; hmint=0. !miz
      allocate ( cgust      (id, jd) )     ; cgust=0.0 !miz
      allocate ( tke        (id, jd) )     ; tke  =0.0 !miz
      allocate ( pblhto     (id, jd) )     ; pblhto=0.0 !miz
      allocate ( rkmo       (id, jd) )     ; rkmo=15.0 !miz
      allocate ( taudpo     (id, jd) )     ; taudpo=28800.    !miz
      allocate ( exist_shconv(id, jd,48) ) ; exist_shconv = 0 !miz
      allocate ( exist_dpconv(id, jd,48) ) ; exist_dpconv = 0 !miz
      allocate ( pblht_prev  (id, jd,48) ) ; pblht_prev   = 0.!miz
      allocate ( hlsrc_prev  (id, jd,48) ) ; hlsrc_prev   = 0.!miz
      allocate ( qtsrc_prev  (id, jd,48) ) ; qtsrc_prev   = 0.!miz
      allocate ( cape_prev   (id, jd,48) ) ; cape_prev    = 0.!miz
      allocate ( cin_prev    (id, jd,48) ) ; cin_prev     = 0.!miz
      allocate ( tke_prev    (id, jd,48) ) ; tke_prev     = 0.!miz

      allocate ( convect    (id, jd) )     ; convect = .false.
      allocate ( radturbten (id, jd, kd))  ; radturbten = 0.0
      allocate ( r_convect  (id, jd) )     ; r_convect   = 0.0
      allocate ( diff_t_clubb(id, jd, kd) ); diff_t_clubb = 0.0

      !<--- yhc, `id` and `jd` correspond to the whole grid
      !          `is,ie,js,je` get allocated based on each thread's section of the grid.
      allocate (  option_edmf2ls_mp )          ; option_edmf2ls_mp = 0  
      allocate (  qadt_edmf (id, jd, kd))  ; qadt_edmf = 0.0  
      allocate (  qldt_edmf (id, jd, kd))  ; qldt_edmf = 0.0  
      allocate (  qidt_edmf (id, jd, kd))  ; qidt_edmf = 0.0  
      allocate (  dqa_edmf  (id, jd, kd))  ; dqa_edmf  = 0.0  
      allocate (  dql_edmf  (id, jd, kd))  ; dql_edmf  = 0.0  
      allocate (  dqi_edmf  (id, jd, kd))  ; dqi_edmf  = 0.0  
      allocate (  diff_t_edmf(id, jd, kd)) ; diff_t_edmf = 0.0  
      allocate (  diff_m_edmf(id, jd, kd)) ; diff_m_edmf = 0.0  
      allocate (  kpbl_edmf  (id, jd))  ;  ; kpbl_edmf= 0     
      allocate (  rdt_before_vdiff_down  (id, jd, kd, ntp))  ;  rdt_before_vdiff_down= 0     
      allocate (  tdt_before_vdiff_down  (id, jd, kd))  ;  tdt_before_vdiff_down= 0     
      allocate (  edmf_mc_full       (id, jd, kd))    ;  edmf_mc_full = 0.0  
      allocate (  edmf_mc_half       (id, jd, kd+1))  ;  edmf_mc_half = 0.0  
      allocate (  edmf_moist_area    (id, jd, kd))    ;  edmf_moist_area = 0.0  
      allocate (  edmf_dry_area      (id, jd, kd))    ;  edmf_dry_area = 0.0  
      allocate (  edmf_moist_humidity(id, jd, kd))    ;  edmf_moist_humidity = 0.0  
      allocate (  edmf_dry_humidity  (id, jd, kd))    ;  edmf_dry_humidity = 0.0  
      !allocate (  (id, jd, kd))  ;  = 0.0  
      !---> yhc

      if (do_cosp) then
!--------------------------------------------------------------------
!    these variables are needed to preserve values of rain fluxes, q and T
!    from the step preceding the COSP call for use in the COSP simulator 
!    on the next step.
!--------------------------------------------------------------------
        allocate ( Precip_flux%fl_lsrain  (id, jd, kd))
        allocate ( Precip_flux%fl_lssnow  (id, jd, kd))
        allocate ( Precip_flux%fl_lsgrpl  (id, jd, kd))
        allocate ( Precip_flux%fl_ccrain  (id, jd, kd))
        allocate ( Precip_flux%fl_ccsnow  (id, jd, kd))
        allocate ( Precip_flux%fl_donmca_snow  (id, jd, kd))
        allocate ( Precip_flux%fl_donmca_rain  (id, jd, kd))
        allocate ( temp_last (id, jd, kd))
        allocate ( q_last    (id, jd, kd))
        Precip_flux%fl_lsrain = 0.
        Precip_flux%fl_lssnow = 0.
        Precip_flux%fl_lsgrpl = 0.
        Precip_flux%fl_ccrain = 0.
        Precip_flux%fl_ccsnow = 0.
        Precip_flux%fl_donmca_rain = 0.
        Precip_flux%fl_donmca_snow = 0.
        temp_last = 0.
        q_last    = 0.
      endif

!-----------------------------------------------------------------------
!    allocate derived-type that stores cloud properties
!    return from moist processes
!-----------------------------------------------------------------------
     call error_mesg('physics_driver_mod', 'number of cloud schemes found = '//trim(string(Exch_ctrl%ncld)), NOTE)

     call alloc_clouds_from_moist_type(Moist_clouds, Exch_ctrl, Atm_block)

!--------------------------------------------------------------------
!    call physics_driver_read_restart to obtain initial values for the module
!    variables. Also register restart fields to be ready for intermediate 
!    restart.
!--------------------------------------------------------------------
      allocate(Restart%Cloud_data(Exch_ctrl%ncld))

      do nc = 1, Exch_ctrl%ncld
        ! restart values allocated on the full domain
        call alloc_cloud_scheme_data_type( Moist_clouds(1)%block(1)%Cloud_data(nc)%scheme_name, &
                                           id, jd, kd, Restart%Cloud_data(nc))
      enddo

      call physics_driver_register_restart (Restart)
      if(file_exist('INPUT/physics_driver.res.nc')) then
         call restore_state(Phy_restart)
         if(in_different_file) call restore_state(Til_restart)
      endif
!---------------------------------------------------------------------
!    convert the real variable (r_convect) indicating columns with 
!    convection to a logical variable (convect). this will be used in 
!    vert_turb_driver_mod.
!---------------------------------------------------------------------
      convect = .false.
      where(r_convect .GT. 0.) 
         convect = .true.
      end where
         
100 FORMAT("CHECKSUM::",A32," = ",Z20)
      outunit = stdout()
      write(outunit,*) 'BEGIN CHECKSUM(physics_driver_init):: '
      write(outunit,100) 'diff_cu_mo             ', mpp_chksum(diff_cu_mo            )
      write(outunit,100) 'pbltop                 ', mpp_chksum(pbltop                )
      write(outunit,100) 'cush                   ', mpp_chksum(cush                  )
      write(outunit,100) 'cbmf                   ', mpp_chksum(cbmf                  )
      write(outunit,100) 'hmint                  ', mpp_chksum(hmint                 )
      write(outunit,100) 'cgust                  ', mpp_chksum(cgust                 )
      write(outunit,100) 'tke                    ', mpp_chksum(tke                   )
      write(outunit,100) 'pblhto                 ', mpp_chksum(pblhto                )
      write(outunit,100) 'rkmo                   ', mpp_chksum(rkmo                  )
      write(outunit,100) 'taudpo                 ', mpp_chksum(taudpo                )
      write(outunit,100) 'exist_shconv           ', mpp_chksum(exist_shconv          )
      write(outunit,100) 'exist_dpconv           ', mpp_chksum(exist_dpconv          )
      write(outunit,100) 'pblht_prev             ', mpp_chksum(pblht_prev            )
      write(outunit,100) 'hlsrc_prev             ', mpp_chksum(hlsrc_prev            )
      write(outunit,100) 'qtsrc_prev             ', mpp_chksum(qtsrc_prev            )
      write(outunit,100) 'cape_prev              ', mpp_chksum(cape_prev             )
      write(outunit,100) 'cin_prev               ', mpp_chksum(cin_prev              )
      write(outunit,100) 'tke_prev               ', mpp_chksum(tke_prev              )
      write(outunit,100) 'diff_t                 ', mpp_chksum(diff_t                )
      write(outunit,100) 'diff_m                 ', mpp_chksum(diff_m                )
      write(outunit,100) 'r_convect              ', mpp_chksum(r_convect             )
   if (doing_prog_clouds) then
      write(outunit,100) 'radturbten             ', mpp_chksum(radturbten            )
   endif
      do nc = 1, size(Restart%Cloud_data,1)
        ! NOTE: the order of the checksums in stdout will be different
        if ( trim(Restart%Cloud_data(nc)%scheme_name).eq.'donner_cell' ) then
          write(outunit,100) 'cell_cld_frac          ', mpp_chksum(Restart%Cloud_data(nc)%cloud_area )
          write(outunit,100) 'cell_liq_amt           ', mpp_chksum(Restart%Cloud_data(nc)%liquid_amt )
          write(outunit,100) 'cell_liq_size          ', mpp_chksum(Restart%Cloud_data(nc)%liquid_size)
          write(outunit,100) 'cell_ice_amt           ', mpp_chksum(Restart%Cloud_data(nc)%ice_amt    )
          write(outunit,100) 'cell_ice_size          ', mpp_chksum(Restart%Cloud_data(nc)%ice_size   )
        endif
        if ( trim(Restart%Cloud_data(nc)%scheme_name).eq.'donner_meso' ) then
          write(outunit,100) 'meso_cld_frac          ', mpp_chksum(Restart%Cloud_data(nc)%cloud_area )
          write(outunit,100) 'meso_liq_amt           ', mpp_chksum(Restart%Cloud_data(nc)%liquid_amt )
          write(outunit,100) 'meso_liq_size          ', mpp_chksum(Restart%Cloud_data(nc)%liquid_size)
          write(outunit,100) 'meso_ice_amt           ', mpp_chksum(Restart%Cloud_data(nc)%ice_amt    )
          write(outunit,100) 'meso_ice_size          ', mpp_chksum(Restart%Cloud_data(nc)%ice_size   )
          write(outunit,100) 'nsum_out               ', mpp_chksum(Restart%Cloud_data(nc)%nsum_out   )
        endif
        if ( trim(Restart%Cloud_data(nc)%scheme_name).eq.'uw_conv' ) then
          write(outunit,100) 'shallow_cloud_area     ', mpp_chksum(Restart%Cloud_data(nc)%cloud_area    )
          write(outunit,100) 'shallow_liquid         ', mpp_chksum(Restart%Cloud_data(nc)%liquid_amt    )
          write(outunit,100) 'shallow_ice            ', mpp_chksum(Restart%Cloud_data(nc)%ice_amt       )
          write(outunit,100) 'shallow_droplet_number ', mpp_chksum(Restart%Cloud_data(nc)%droplet_number)
          write(outunit,100) 'shallow_ice_number     ', mpp_chksum(Restart%Cloud_data(nc)%ice_number    )
        endif
        if ( trim(Restart%Cloud_data(nc)%scheme_name).eq.'strat_cloud' ) then
          write(outunit,100) 'lsc_cloud_area         ', mpp_chksum(Restart%Cloud_data(nc)%cloud_area    )
          write(outunit,100) 'lsc_liquid             ', mpp_chksum(Restart%Cloud_data(nc)%liquid_amt    )
          write(outunit,100) 'lsc_ice                ', mpp_chksum(Restart%Cloud_data(nc)%ice_amt       )
          write(outunit,100) 'lsc_droplet_number     ', mpp_chksum(Restart%Cloud_data(nc)%droplet_number)
          write(outunit,100) 'lsc_ice_number         ', mpp_chksum(Restart%Cloud_data(nc)%ice_number    )
          write(outunit,100) 'lsc_snow               ', mpp_chksum(Restart%Cloud_data(nc)%snow          )
          write(outunit,100) 'lsc_rain               ', mpp_chksum(Restart%Cloud_data(nc)%rain          )
          write(outunit,100) 'lsc_snow_size          ', mpp_chksum(Restart%Cloud_data(nc)%snow_size     )
          write(outunit,100) 'lsc_rain_size          ', mpp_chksum(Restart%Cloud_data(nc)%rain_size     )
        endif
      enddo ! nc

      do nb = 1, Atm_block%nblks
        ibs = Atm_block%ibs(nb)-Atm_block%isc+1
        ibe = Atm_block%ibe(nb)-Atm_block%isc+1
        jbs = Atm_block%jbs(nb)-Atm_block%jsc+1
        jbe = Atm_block%jbe(nb)-Atm_block%jsc+1

        !-- copy cloud data from restart
        do nc = 1, size(Restart%Cloud_data,1)

          ! common to all cloud schemes
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%cloud_area     = Restart%Cloud_data(nc)%cloud_area    (ibs:ibe,jbs:jbe,:)
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%liquid_amt     = Restart%Cloud_data(nc)%liquid_amt    (ibs:ibe,jbs:jbe,:)
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_amt        = Restart%Cloud_data(nc)%ice_amt       (ibs:ibe,jbs:jbe,:)
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%droplet_number = Restart%Cloud_data(nc)%droplet_number(ibs:ibe,jbs:jbe,:)
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name    = Restart%Cloud_data(nc)%scheme_name

          ! properties specific to large-scale/stratiform clouds
          if (trim(Restart%Cloud_data(nc)%scheme_name) .eq. 'strat_cloud') then
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_number = Restart%Cloud_data(nc)%ice_number (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%rain       = Restart%Cloud_data(nc)%rain       (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%snow       = Restart%Cloud_data(nc)%snow       (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%rain_size  = Restart%Cloud_data(nc)%rain_size  (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%snow_size  = Restart%Cloud_data(nc)%snow_size  (ibs:ibe,jbs:jbe,:)
          endif
  
          ! properties specific to donner deep clouds (both cell and meso)
          if (trim(Restart%Cloud_data(nc)%scheme_name) .eq. 'donner_cell' .or. &
              trim(Restart%Cloud_data(nc)%scheme_name) .eq. 'donner_meso') then
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%liquid_size = Restart%Cloud_data(nc)%liquid_size (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_size    = Restart%Cloud_data(nc)%ice_size    (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%nsum_out    = Restart%Cloud_data(nc)%nsum_out    (ibs:ibe,jbs:jbe)
          endif

          ! properties specific to uw shallow convective clouds
          if (trim(Restart%Cloud_data(nc)%scheme_name) .eq. 'uw_conv') then
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_number = Restart%Cloud_data(nc)%ice_number (ibs:ibe,jbs:jbe,:)
          endif

        enddo

        !--- return trs to the blocked data structure
        Physics%block(nb)%q = trs(ibs:ibe,jbs:jbe,:,1:ntp)
        Physics_tendency%block(nb)%qdiag = trs(ibs:ibe,jbs:jbe,:,ntp+1:nt)
      enddo
      deallocate (trs, phalf)

      vers = restart_versions(size(restart_versions(:)))

!---------------------------------------------------------------------
!    if desired, define variables to return  the restart fields of 
!    diff_m and diff_t.
!---------------------------------------------------------------------
      if (present(difft)) then
        difft = diff_t
      endif
      if (present(diffm)) then
        diffm = diff_m
      endif

!---------------------------------------------------------------------
!    register module diagnostics
!---------------------------------------------------------------------

      id_tdt_phys_vdif_dn = register_diag_field ( mod_name,    &
         'tdt_phys_vdif_dn', axes(1:3), Time,                  &
         'temperature tendency from physics driver vdif down', &
         'K/s', missing_value=missing_value)

      id_tdt_phys_vdif_up = register_diag_field ( mod_name,    &
         'tdt_phys_vdif_up', axes(1:3), Time,                  &
         'temperature tendency from physics driver vdif up',   &
         'K/s', missing_value=missing_value)

      id_tdt_phys_turb = register_diag_field ( mod_name,       &
         'tdt_phys_turb', axes(1:3), Time,                     &
         'temperature tendency from physics driver vdif turb', &
         'K/s', missing_value=missing_value)

      id_tdt_phys_moist = register_diag_field ( mod_name,            &
         'tdt_phys_moist', axes(1:3), Time,                          &
         'temperature tendency from physics driver moist processes', &
         'K/s', missing_value=missing_value)

      id_tdt_phys = register_diag_field ( mod_name,            &
         'tdt_phys', axes(1:3), Time,                          &
         'temperature tendency from physics ', &
         'K/s', missing_value=missing_value)

      !<--- yhc
      id_diff_t_vdif = register_diag_field (mod_name, 'diff_t_vdif', axes(1:3), Time, &
                       'heat diff coeffs used by vdif', 'K/m/s' , &
                       missing_value=missing_value )

      id_diff_m_vdif = register_diag_field (mod_name, 'diff_m_vdif', axes(1:3), Time, &
                       'momentum diff coeffs used by vdif', 'm2/s' , &
                       missing_value=missing_value )

      id_qldt_vdif = register_diag_field (mod_name, 'qldt_vdif', axes(1:3), Time, &
                     'liquid water tendency from vert diff', 'kg/kg/s' , &
                     missing_value=missing_value )
    
      id_qadt_vdif = register_diag_field (mod_name, 'qadt_vdif', axes(1:3), Time, &
                     'cloud fraction tendency from vert diff', '1/s' , &
                     missing_value=missing_value )
    
      id_qidt_vdif = register_diag_field (mod_name, 'qidt_vdif', axes(1:3), Time, &
                     'ice water tendency from vert diff', 'kg/kg/s' , &
                     missing_value=missing_value )
    
      id_qdt_vdif_test = register_diag_field (mod_name, 'qdt_vdif_test', axes(1:3), Time, &
                     'spec humid tendency from vert diff (test)', 'kg/kg/s' , &
                     missing_value=missing_value )

      id_tdt_vdif_test = register_diag_field (mod_name, 'tdt_vdif_test', axes(1:3), Time, &
                     'Temperature tendency from vert diff (test)', 'K/s' , &
                     missing_value=missing_value )

      id_tdt_radturb = register_diag_field (mod_name, 'tdt_radturb', axes(1:3), Time, &
                     'tempearture tendency from rad and vert diff', 'K/s' , &
                     missing_value=missing_value )

      id_tt_phys = register_diag_field (mod_name, 'tt_phys', axes(1:3), Time, &
                    'temperature in physics_down', 'K' , &
                    missing_value=missing_value )

      id_qq_phys = register_diag_field (mod_name, 'qq_phys', axes(1:3), Time, &
                    'specific humidity in physics_down', 'kg/kg' , &
                    missing_value=missing_value )

      id_LWP0    = register_diag_field (mod_name, 'LWP0', axes(1:2), Time, &
                    'Liquid water path (phy_in)', 'kg/m2' , &
                     missing_value=missing_value )
      !---> yhc

     !-------- CMIP diagnostics --------
      ID_pfull = register_cmip_diag_field_3d ( mod_name, 'pfull', Time, &
                                     'Pressure on Model Levels', 'Pa', &
                                      standard_name = 'air_pressure')

      ID_phalf = register_cmip_diag_field_3d ( mod_name, 'phalf', Time, &
                                'Pressure on Model Half-Levels', 'Pa', &
                              standard_name='air_pressure', axis='half')

      ID_zg = register_cmip_diag_field_3d ( mod_name, 'zg', Time, &
                                      'Geopotential Height', 'm', &
                             standard_name = 'geopotential_height')

     !ID_zfull = register_cmip_diag_field_3d ( mod_name, 'zfull', Time, &
     !                            'Altitude of Model Full-Levels', 'm', &
     !                standard_name = 'height_above_reference_ellipsoid')

     !-------- CMIP diagnostics (tendencies due to physics) --------
      ID_tntmp = register_cmip_diag_field_3d ( mod_name, 'tntmp', Time, &
                  'Tendency of Air Temperature due to Model Physics', 'K s-1', & 
             standard_name='tendency_of_air_temperature_due_to_model_physics' )

      nsphum = get_tracer_index(MODEL_ATMOS,'sphum')
      if (nsphum /= NO_TRACER) then
        ID_tnhusmp = register_cmip_diag_field_3d ( mod_name, 'tnhusmp', Time, &
                    'Tendency of Specific Humidity due to Model Physics', 's-1', &
               standard_name='tendency_of_specific_humidity_due_to_model_physics' )
      endif

      allocate (id_tracer_phys(ntp))
      allocate (id_tracer_phys_vdif_dn(ntp))
      allocate (id_tracer_phys_vdif_up(ntp))
      allocate (id_tracer_phys_turb(ntp))
      allocate (id_tracer_phys_moist(ntp))

      do n = 1,ntp

        call get_tracer_names (MODEL_ATMOS, n, name = tracer_name,  &
                               units = tracer_units)
        
        diaglname = trim(tracer_name)//  &
                    ' tendency from physics'
        id_tracer_phys(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

        diaglname = trim(tracer_name)//  &
                    ' tendency from physics driver vdif down'
        id_tracer_phys_vdif_dn(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys_vdif_dn',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

        diaglname = trim(tracer_name)//  &
                    ' tendency from physics driver vdif up'
        id_tracer_phys_vdif_up(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys_vdif_up',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

        diaglname = trim(tracer_name)//  &
                    ' tendency from physics driver vert turb'
        id_tracer_phys_turb(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys_turb',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

        diaglname = trim(tracer_name)//  &
                    ' tendency from physics driver moist processes'
        id_tracer_phys_moist(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys_moist',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

      end do

      call monin_obukhov_init
!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!-----------------------------------------------------------------------

 end subroutine physics_driver_init


!######################################################################
! <SUBROUTINE NAME="physics_driver_down_time_vary">
!  <OVERVIEW>
!    physics_driver_time_vary makes sure that all time-dependent, spacially-
!    independent calculations are completed before entering window or thread
!    loops. Resultant fields are usually saved as module variables in the
!    module where needed.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_time_vary makes sure that all time-dependent, spacially-
!    independent calculations are completed before entering window or thread
!    loops. Resultant fields are usually saved as module variables in the
!    module where needed.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_down_time_vary (Time, Time_next)
!
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   time of next time step
!  </IN>
! </SUBROUTINE>
!

subroutine physics_driver_down_time_vary (Time, Time_next, dt)

!---------------------------------------------------------------------
!    physics_driver_down_time_vary makes sure that all time-dependent, 
!    spacially-independent calculations are completed before entering window 
!    or thread loops. Resultant fields are usually saved as module variables in 
!    the module where needed.
!-----------------------------------------------------------------------

type(time_type),         intent(in)             :: Time, Time_next
real,                    intent(in)             :: dt

type(time_type) :: Time_last
!---------------------------------------------------------------------      
!------------------------------------------------------------------------
!    call damping_driver_time_vary to update the counter determining when
!    the convective drag module will be again called.
!------------------------------------------------------------------------
      call damping_driver_time_vary (dt)

!------------------------------------------------------------------------
!    call atmos_tracer_driver_time_vary to obtain values from the tracer 
!    climatology and emission data at the appropriate time, if needed. 
!------------------------------------------------------------------------
      call atmos_tracer_driver_time_vary (Time)

!-------------------------------------------------------------------------      

end subroutine physics_driver_down_time_vary



!######################################################################

subroutine physics_driver_down_endts(is,js)

integer, intent(in)  :: is,js

!-----------------------------------------------------------------------
!    call the component xxx_endts routines to perform needed updates to
!    primarily flag and counter variables at the end of the time step, 
!    after all spacial-dependent calculations are completed.
!-----------------------------------------------------------------------
      call damping_driver_endts
      call atmos_tracer_driver_endts

!--------------------------------------------------------------------
!    set a flag to indicate that this check was done and need not be
!    done again.
!--------------------------------------------------------------------
      do_check_args = .false.

end subroutine physics_driver_down_endts

!###################################################################

subroutine physics_driver_up_time_vary (Time, Time_next, dt, &
                                        step_to_call_cosp_in)

!---------------------------------------------------------------------
!    physics_driver_up_time_vary makes sure that all time-dependent, 
!    spacially-independent calculations are completed before entering 
!    window or thread loops. Resultant fields are usually saved as 
!    module variables in the module where needed.
!-----------------------------------------------------------------------

type(time_type),         intent(in)             :: Time
type(time_type),         intent(in)             :: Time_next
real,                    intent(in)             :: dt
logical,                 intent(in)             :: step_to_call_cosp_in

   
!----------------------------------------------------------------------
!    save the flag indicating if this is step to call cosp.
!----------------------------------------------------------------------
    step_to_call_cosp = step_to_call_cosp_in

    if (do_moist_processes) then
!----------------------------------------------------------------------
!    call aerosol_time_vary to retrieve appropriate aerosol fields from
!    the climatology, if that source of aerosol is being used.
!----------------------------------------------------------------------
      call aerosol_time_vary (Time, Aerosol_cld)
!----------------------------------------------------------------------
!    call moist_processes_time_vary to pass needed time-dependent fields 
!    to subordinate modules.
!----------------------------------------------------------------------
      call moist_processes_time_vary (dt)
    endif
!----------------------------------------------------------------------
!    call cosp_driver_time_vary to obtain satellite location at current
!    time if orbital data is being collected.
!----------------------------------------------------------------------
    if (do_cosp) call cosp_driver_time_vary (Time_next)

!----------------------------------------------------------------------      

end subroutine physics_driver_up_time_vary


!######################################################################

subroutine physics_driver_up_endts 

!-----------------------------------------------------------------------
!    call the component xxx_endts routines to perform needed updates to
!    primarily flag and counter variables at the end of the time step, 
!    after all spacial-dependent calculations are completed.
!-----------------------------------------------------------------------
    if (do_cosp) call cosp_driver_endts
    if (do_moist_processes) then
      call moist_processes_endts 
      call aerosol_endts (Aerosol_cld)
    endif

end subroutine physics_driver_up_endts


!######################################################################


!######################################################################
! <SUBROUTINE NAME="physics_driver_down">
!  <OVERVIEW>
!    physics_driver_down calculates "first pass" physics tendencies,
!    associated with radiation, damping and turbulence, and obtains
!    the vertical diffusion tendencies to be passed to the surface and
!    used in the semi-implicit vertical diffusion calculation.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_down calculates "first pass" physics tendencies,
!    associated with radiation, damping and turbulence, and obtains
!    the vertical diffusion tendencies to be passed to the surface and
!    used in the semi-implicit vertical diffusion calculation.    
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_down (is, ie, js, je,                       &
!                                Time_prev, Time, Time_next,           &
!                                lat, lon, area,                       &
!                                p_half, p_full, z_half, z_full,       &
!                                u, v, t, q, r, um, vm, tm, qm, rm,    &
!                                frac_land, rough_mom,                 &
!                                albedo,    t_surf_rad, u_ref, v_ref,  &
!                                t_ref, q_ref,                         &
!                                u_star,    b_star, q_star,            &
!                                dtau_du,  dtau_dv,  tau_x,  tau_y,    &
!                                udt, vdt, tdt, qdt, rdt,              &
!                                flux_sw,  flux_lw,  coszen,  gust,    &
!                                Surf_diff
!  </TEMPLATE>
!  <IN NAME="Time_prev" TYPE="time_type">
!   previous time, for variable um, vm, tm, qm, rm
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   next time, used for diagnostics
!  </IN>
!  <IN NAME="lat" TYPE="real">
!   array of model latitudes at model points [radians]
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   array of model longitudes at model points [radians]
!  </IN>
!  <IN NAME="area" TYPE="real">
!   grid box area - current not used
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   pressure at model interface levels (offset from t,q,u,v,r)
!  </IN>
!  <IN NAME="p_full" TPYE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="z_half" TYPE="real">
!   height at model interface levels
!  </IN>
!  <IN NAME="z_full" TPYE="real">
!   height at full levels
!  </IN>
!  <IN NAME="u" TYPE="real">
!   zonal wind at current time step
!  </IN>
!  <IN NAME="v" TYPE="real">
!   meridional wind at current time step
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at current time step
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity at current time step
!  </IN>
!  <IN NAME="r" TPYE="real">
!   multiple 3d tracer fields at current time step
!  </IN>
!  <IN NAME="um" TYPE="real">
!   zonal wind at previous time step
!  </IN>
!  <IN NAME="vm" TYPE="real">
!   meridional wind at previous time step
!  </IN>
!  <IN NAME="tm" TYPE="real">
!   temperature at previous time step
!  </IN>
!  <IN NAME="qm" TYPE="real">
!   specific humidity at previous time step
!  </IN>
!  <IN NAME="rm" TPYE="real">
!   multiple 3d tracer fields at previous time step
!  </IN>
!  <INOUT NAME="rd" TYPE="real">
!   multiple 3d diagnostic tracer fields 
!  </INOUT>
!  <IN NAME="frac_land" TYPE="real">
!   fraction of land coverage in a model grid point
!  </IN>
!  <IN NAME="rough_mom" TYPE="real">
!   boundary layer roughness
!  </IN>
!  <IN NAME="albedo" TYPE="real">
!   surface albedo
!  </IN>
!  <IN NAME="t_surf_rad" TYPE="real">
!   surface radiative temperature
!  </IN>
!  <IN NAME="u_ref" TYPE="real">
!   10m zonal wind
!  </IN>
!  <IN NAME="v_ref" TYPE="real">
!   10m meridional wind
!  </IN>
!  <IN NAME="u_star" TYPE="real">
!   boundary layer wind speed (frictional speed)
!  </IN>
!  <IN NAME="b_star" TYPE="real">
!   ???
!  </IN>
!  <IN NAME="q_star" TYPE="real">
!   boundary layer specific humidity
!  </IN>
!  <IN NAME="dtau_du" TYPE="real">
!   derivative of zonal surface stress w.r.t zonal wind speed
!  </IN>
!  <IN NAME="dtau_dv" TYPE="real">
!   derivative of meridional surface stress w.r.t meridional wind speed
!  </IN>
!  <INOUT NAME="tau_x" TYPE="real">
!   boundary layer meridional component of wind shear
!  </INOUT>
!  <INOUT NAME="tau_y" TYPE="real">
!   boundary layer zonal component of wind shear
!  </INOUT>
!  <INOUT NAME="udt" TYPE="real">
!   zonal wind tendency
!  </INOUT>
!  <INOUT NAME="vdt" TYPE="real">
!   meridional wind tendency
!  </INOUT>
!  <INOUT NAME="tdt" TYPE="real">
!   temperature tendency
!  </INOUT>
!  <INOUT NAME="qdt" TYPE="real">
!   moisture tracer tendencies
!  </INOUT>
!  <INOUT NAME="rdt" TYPE="real">
!   multiple tracer tendencies
!  </INOUT>
!  <OUT NAME="flux_sw" TYPE="real">
!   Shortwave flux from radiation package
!  </OUT>
!  <OUT NAME="flux_lw" TYPE="real">
!   Longwave flux from radiation package
!  </OUT>
!  <OUT NAME="coszen" TYPE="real">
!   cosine of zenith angle
!  </OUT>
!  <OUT NAME="gust" TYPE="real">
!  </OUT>
!  <INOUT NAME="Surf_diff" TYPE="surface_diffusion_type">
!   Surface diffusion 
!  </INOUT>
!
! </SUBROUTINE>
!
subroutine physics_driver_down (is, ie, js, je, npz,              &
                                Time_prev, Time, Time_next,       &
                                lat, lon, area,                   &
                                Physics_input_block,              &
                                frac_land, rough_mom,             &
                                frac_open_sea,                    &
                                albedo,                           &
                                t_surf_rad, u_ref, v_ref,         & !bqx+ u_ref, v_ref
                                t_ref, q_ref,                     &
                                u_star,    b_star, q_star,        &
                                dtau_du, dtau_dv,  tau_x,  tau_y, &
                                Physics_tendency_block,           &
                                Surf_diff,                        &
                                gust,                             &
                                Rad_flux_control,                 &
                                Rad_flux_block,                   &
                                diffm, difft  )

!---------------------------------------------------------------------
!    physics_driver_down calculates "first pass" physics tendencies,
!    associated with radiation, damping and turbulence, and obtains
!    the vertical diffusion tendencies to be passed to the surface and
!    used in the semi-implicit vertical diffusion calculation.
!-----------------------------------------------------------------------

integer,                 intent(in)             :: is, ie, js, je, npz
type(time_type),         intent(in)             :: Time_prev, Time, Time_next
real,dimension(:,:),     intent(in)             :: lat, lon, area
type(physics_input_block_type), intent(in)      :: Physics_input_block
real,dimension(:,:),     intent(in)             :: frac_land,   &
                                                   rough_mom, &
                                                   albedo, t_surf_rad, &
                                                   u_ref, v_ref, & !bqx+
                                                   t_ref, q_ref, &  ! cjg: PBL depth mods
                                                   u_star, b_star, &
                                                   q_star, dtau_du,   &
                                                   dtau_dv, frac_open_sea
real,dimension(:,:),     intent(inout)          :: tau_x,  tau_y
type(physics_tendency_block_type), intent(inout):: Physics_tendency_block
real,dimension(:,:),     intent(out)            :: gust
type(surf_diff_type),    intent(inout)          :: Surf_diff
type(radiation_flux_control_type),  intent(in)  :: Rad_flux_control
type(radiation_flux_block_type),    intent(in)  :: Rad_flux_block
real,  dimension(:,:,:), intent(out)  ,optional :: diffm, difft 

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      npz            number of model levels
!      Time_prev      previous time, for variables um,vm,tm,qm,rm 
!                     (time_type)
!      Time           current time, for variables u,v,t,q,r  (time_type)
!      Time_next      next time, used for diagnostics   (time_type)
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      Physics_input_block  derived type variable containing: 
!         1) p_half         pressure at half levels (offset from t,q,u,v,r)
!                          [ Pa ]
!         2) p_full         pressure at full levels [ Pa }
!         3) z_half         height at half levels [ m ]
!         4) z_full         height at full levels [ m ]
!         5) u              zonal wind at current time step [ m / s ]
!         6) v              meridional wind at current time step [ m / s ]
!         7) t              temperature at current time step [ deg k ]
!         9) q              multiple 3d tracer fields at current time step
!        10) um,vm          zonal and meridional wind at previous time step
!        11) tm             temperature at previous time step
!        12) qm             multiple 3d tracer fields at previous time step
!      Rad_flux_control
!      Rad_flux_block
!      frac_land
!      rough_mom
!      frac_open_sea
!      albedo
!      t_surf_rad
!      t_ref
!      q_ref
!      u_star
!      b_star
!      q_star
!      dtau_du
!      dtau_dv
!
!  intent(inout) variables:
!
!      tau_x
!      tau_y
!      Physics_tendency_block derived type variable containing:
!          1) u_dt           zonal wind tendency [ m / s**2 ]
!          2) v_dt           meridional wind tendency [ m / s**2 ]
!          3) t_dt           temperature tendency [ deg k / sec ]
!          4) q_dt           multiple tracer tendencies 
!                            (index 1 = specific humidity) 
!                            [ unit / unit / sec ]
!          5) qdiag          multiple 3d diagnostic tracer fields 
!                            [ unit / unit ]
!      Surf_diff      surface_diffusion_type variable
!
!   intent(out) variables:
!
!      gust
!
!   intent(out), optional variables:
!
!      diffm
!      difft
!
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!    local variables:

      real, dimension(ie-is+1,je-js+1,npz)   :: diff_t_vert, diff_m_vert
      real, dimension(ie-is+1,je-js+1,npz)   :: tdt_rad, tdt_lw
      real, dimension(ie-is+1,je-js+1,npz+1) :: lphalf
      real, dimension(ie-is+1,je-js+1)       :: z_pbl
      integer                                :: sec, day, n, nextinct
      real                                   :: dt, alpha, dt2
      logical                                :: used

!---> h1g, 2015-08-11
      real, dimension(ie-is+1,je-js+1) :: tke_avg
!<--- h1g, 2015-08-11

!---------------------------------------------------------------------
!   local variables:
!
!      diff_t_vert     vertical diffusion coefficient for temperature
!                      calculated on the current step
!      diff_m_vert     vertical diffusion coefficient for momentum   
!                      calculated on the current step
!      z_pbl           height of planetary boundary layer
!      sec, day        second and day components of the time_type 
!                      variable
!      dt              model physics time step [ seconds ]
!      alpha           ratio of physics time step to diffusion-smoothing
!                      time scale
!
!---------------------------------------------------------------------
      real, dimension(:,:,:,:), pointer :: r, rm
      real, dimension(:,:,:), pointer :: p_full, p_half, z_full, z_half
      real, dimension(:,:,:), pointer :: udt, vdt, tdt
      real, dimension(:,:,:,:), pointer :: rdt, rdiag
      real, dimension(:,:,:), pointer :: u, v, t, um, vm, tm 

!<-- yhc
      real, dimension( size(Physics_tendency_block%u_dt,1), &
                       size(Physics_tendency_block%u_dt,2), &
                       size(Physics_tendency_block%u_dt,3)  )  :: &
        udt_before_vdiff_down, vdt_before_vdiff_down, tdt_before_vdiff_down2,  &
        tdt_dum, tdt_mf_dum, tdt_rad_only, &
        pmass, &
        rdt_dum1, rdt_dum2

      real, dimension( size(Physics_tendency_block%q_dt,1), &
                       size(Physics_tendency_block%q_dt,2), &
                       size(Physics_tendency_block%q_dt,3), &
                       size(Physics_tendency_block%q_dt,4)  )  :: &
        rdt_before_vdiff_down2, rdt_dum

      real, dimension( size(Physics_tendency_block%u_dt,1), &
                       size(Physics_tendency_block%u_dt,2)) :: &
        tau_x_before_vdiff_down, tau_y_before_vdiff_down,  &
        shflx, lhflx, u_flux, v_flux
      integer ii,jj,kk,rr, kx
!--> yhc

!---------------------------------------------------------------------
!    set up local pointers into the physics input and physics tendency
!    blocks.
!---------------------------------------------------------------------
      u => Physics_input_block%u
      v => Physics_input_block%v
      t => Physics_input_block%t
      r => Physics_input_block%q
      if (associated(Physics_input_block%um)) then
        um => Physics_input_block%um
        vm => Physics_input_block%vm
        tm => Physics_input_block%tm
        rm => Physics_input_block%qm
      else
        um => Physics_input_block%u
        vm => Physics_input_block%v
        tm => Physics_input_block%t
        rm => Physics_input_block%q
      endif
      p_full => Physics_input_block%p_full
      p_half => Physics_input_block%p_half
      z_full => Physics_input_block%z_full
      z_half => Physics_input_block%z_half
      udt => Physics_tendency_block%u_dt
      vdt => Physics_tendency_block%v_dt
      tdt => Physics_tendency_block%t_dt
      rdt => Physics_tendency_block%q_dt
      rdiag => Physics_tendency_block%qdiag

!<--- yhc, 2022-06-28
      !------- temperature in physics_down (units: K) at full level -------
      if ( id_tt_phys > 0) then
        used = send_data (id_tt_phys, t, Time_next, is, js, 1 )
      endif

      !------- specific humidity in physics_down (units: kg/kg) at full level -------
      if ( id_qq_phys > 0) then
        used = send_data (id_qq_phys, r(:,:,:,1), Time_next, is, js, 1 )
      endif

      !------- LWP0 -------
      if ( id_LWP0 > 0) then
        do kk=1,size(Physics_tendency_block%q_dt,3)
          pmass(:,:,kk) =    &
                    (Physics_input_block%p_half(:,:,kk+1) - Physics_input_block%p_half(:,:,kk))/grav
        end do
        call column_diag    &
             (id_LWP0, is, js, Time_next, Physics_input_block%q(:,:,:,nql), 1.0,   &
                                                          pmass)
      endif
!---> yhc, 2022-06-28

!<--- yhc, 2022-07-15
      !--- compute vertical advection tendencies offline and stop the model
      if (do_vert_adv_tend_offline) then
        call compute_vert_adv_tend_offline(dt, t, temp_vert_advec_scheme, tracer_vert_advec_scheme) 
      end if
!---> yhc, 2022-07-15

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
                         'module has not been initialized', FATAL)
      endif

!---------------------------------------------------------------------
!    check the size of the input arguments. this is only done on the
!    first call to physics_driver_down.
!---------------------------------------------------------------------
      if (do_check_args) call check_args  &
                   (lat, lon, area, p_half, p_full, z_half, z_full, &
                    u, v, t, r(:,:,:,1), r, um, vm, tm, r(:,:,:,1), rm, &
                    udt, vdt, tdt, rdt(:,:,:,1), rdt, rdiag)

!---------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
      call get_time (Time_next - Time_prev, sec, day)
      dt = real(sec + day*86400)

!---------------------------------------------------------------------
!     save cmip diagnostics
!---------------------------------------------------------------------
      if (query_cmip_diag_id(ID_pfull)) then
         used = send_cmip_data_3d (ID_pfull, p_full, Time_next, is, js, 1)
      endif

      if (query_cmip_diag_id(ID_phalf)) then
         used = send_cmip_data_3d (ID_phalf, p_half, Time_next, is, js, 1)
      endif

     ! only allow output of geop height on model levels
     ! output from dycore for pressure levels
      if (query_cmip_diag_id(ID_zg)) then
         used = send_cmip_data_3d (ID_zg, z_full, Time_next, is, js, 1)
      endif

!<-- yhc, determine whether writing out the selected column 
  do_writeout_column = .false.
  if (do_writeout_column_nml) then

    !--- for global simulations
    if (ii_write.ne.-999 .and. jj_write.ne.-999) then
      do_writeout_column = .true.

      if (lat_write.ne.-999.99 .and. lon_write.ne.-999.99) then

        lat_lower = lat_write - 0.000001
        lat_upper = lat_write + 0.000001
        lon_lower = lon_write - 0.000001
        lon_upper = lon_write + 0.000001

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
!--> yhc, determine whether writing out the selected column

!<-- yhc
  if (do_writeout_column) then
        write(6,*) '-------------- i,j,',ii_write,jj_write
        write(6,*) 'lat',lat (ii_write,jj_write)
        write(6,*) 'lon',lon (ii_write,jj_write)
        write(6,*) 'data t_physics_down_begin/'    ,t(ii_write,jj_write,:)
        write(6,*) 'data q_physics_down_begin/'    ,r(ii_write,jj_write,:,1)
        write(6,*) 'data udt_physics_down_begin/'    ,udt(ii_write,jj_write,:)
        write(6,*) 'data vdt_physics_down_begin/'    ,vdt(ii_write,jj_write,:)
        write(6,*) 'data tdt_physics_down_begin/'    ,tdt(ii_write,jj_write,:)
        write(6,*) 'data qdt_physics_down_begin/'    ,rdt(ii_write,jj_write,:,1)
        write(6,*) 'data qadt_physics_down_begin/'    ,rdt(ii_write,jj_write,:,nqa)
        write(6,*) 'data qldt_physics_down_begin/'    ,rdt(ii_write,jj_write,:,nql)
        write(6,*) 'data qidt_physics_down_begin/'    ,rdt(ii_write,jj_write,:,nqi)
        write(6,*) 'data qndt_physics_down_begin/'    ,rdt(ii_write,jj_write,:,nqn)
        write(6,*) 'data omega_physics_down_begin/'    ,Physics_input_block%omega(ii_write,jj_write,:)
        write(6,*) 'radturbten_physics_down_begin/', radturbten(ii_write,jj_write,:)
        write(6,*) 'Rad_flux_block%tdt_rad', Rad_flux_block%tdt_rad(ii_write,jj_write,:)
        do rr=1, size(Physics_tendency_block%q_dt,4)
          write(6,*) 'data dn_begin, rr, r'    ,rr, r(ii_write,jj_write,:,rr)
        enddo
  endif
!--> yhc

!---------------------------------------------------------------------

!rab      if(do_grey_radiation) then !rif:(09/10/09) 
!rab        call grey_radiation(is, js, Time, Time_next, lat, lon, phalfgrey, albedo, t_surf_rad, t, tdt, flux_sw, flux_lw)
!rab        coszen = 1.0
!rab        flux_sw_dir     = R1*flux_sw
!rab        flux_sw_dif     = R2*flux_sw
!rab        flux_sw_vis_dir = R3*flux_sw
!rab        flux_sw_vis_dif = R4*flux_sw
!rab      endif

      if (do_radiation) then
        radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) + Rad_flux_block%tdt_rad(:,:,:)
        surf_diff%tdt_rad(is:ie,js:je,:)=Rad_flux_block%tdt_rad(:,:,:) !miz
      endif
#ifdef SCM
! Option to add SCM radiative tendencies from forcing to Rad_flux_block%tdt_lw
! and radturbten

      if (use_scm_rad) then
        call add_scm_tdtlw( Rad_flux_block%tdt_lw )
        call add_scm_tdtlw( radturbten (is:ie,js:je,:) )
        call add_scm_tdtsw( radturbten (is:ie,js:je,:) )
      endif

#endif

      !<--- yhc
      if (do_bomex_radf) then  ! yhc add, consistent with tdt_radf in history files
        kx = size(Physics_tendency_block%q_dt,3)
        radturbten = 0.
        radturbten(is:ie,js:je,kx-10:kx)    = -2. /86400.  ! prescribe -2K/day cooling rate
        radturbten(is:ie,js:je,kx-11:kx-11) = -0.4/86400.  ! -0.4 K/day
      endif

      tdt_rad_only(:,:,:) = radturbten(is:ie,js:je,:) ! yhc save tdt_rad
      !---> yhc

!----------------------------------------------------------------------
!    call damping_driver to calculate the various model dampings that
!    are desired. 
!----------------------------------------------------------------------
      z_pbl(:,:) = pbltop(is:ie,js:je) 
      call mpp_clock_begin ( damping_clock )
      call damping_driver (is, js, lat, Time_next, dt, area,        &
                           p_full, p_half, z_full, z_half,          &
                           um, vm, tm, rm(:,:,:,1), rm(:,:,:,1:ntp),&
                           u_ref, v_ref, z_pbl,                     &  !bqx+
                           udt, vdt, tdt, rdt(:,:,:,1), rdt)
     call mpp_clock_end ( damping_clock )

!<-- yhc
  if (do_writeout_column) then
        write(6,*) '-------------- i,j,',ii_write,jj_write
        write(6,*) 'lat',lat (ii_write,jj_write)
        write(6,*) 'lon',lon (ii_write,jj_write)
        write(6,*) 'data tdt_rad_only/'    ,tdt_rad_only(ii_write,jj_write,:)
        write(6,*) 'data t_damping/'    ,t(ii_write,jj_write,:)
        write(6,*) 'data q_damping/'    ,r(ii_write,jj_write,:,1)
        write(6,*) 'data tdt_damping/'    ,tdt(ii_write,jj_write,:)
        write(6,*) 'data qdt_damping/'    ,rdt(ii_write,jj_write,:,1)
        write(6,*) 'data qadt_damping/'    ,rdt(ii_write,jj_write,:,nqa)
        write(6,*) 'data qldt_damping/'    ,rdt(ii_write,jj_write,:,nql)
        write(6,*) 'data qidt_damping/'    ,rdt(ii_write,jj_write,:,nqi)
        write(6,*) 'data qndt_damping/'    ,rdt(ii_write,jj_write,:,nqn)
  endif
!--> yhc

!---------------------------------------------------------------------
!    call vert_turb_driver to calculate diffusion coefficients. save
!    the planetary boundary layer height on return.
!---------------------------------------------------------------------

      if (id_tdt_phys_turb > 0) then
        used = send_data ( id_tdt_phys_turb, -2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_turb(n) > 0) then
          used = send_data ( id_tracer_phys_turb(n), -2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

      call mpp_clock_begin ( turb_clock )
      call vert_turb_driver (is, js, Time, Time_next, dt,            &
                             Rad_flux_block%tdt_lw, frac_land,  &
                             p_half, p_full, z_half, z_full,         &
                             t_ref, q_ref,                           &  ! cjg: PBL depth mods
                             u_star, b_star, q_star, rough_mom,      &
                             lat, convect(is:ie,js:je),              &
                             u, v, t, r(:,:,:,1), r, um, vm,                  &
                             tm, rm(:,:,:,1), rm, rdiag,                      &
                             udt, vdt, tdt, rdt(:,:,:,1), rdt,                &
                             diff_t_vert, diff_m_vert, gust, z_pbl, tke_avg = tke_avg)
     call mpp_clock_end ( turb_clock )
     pbltop(is:ie,js:je) = z_pbl(:,:)
     tke   (is:ie,js:je) = tke_avg(:,:)

      if (id_tdt_phys_turb > 0) then
        used = send_data ( id_tdt_phys_turb, +2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_turb(n) > 0) then
          used = send_data ( id_tracer_phys_turb(n), +2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

  if (do_writeout_column) then
        write(6,*) '-------------- i,j,',ii_write,jj_write
        write(6,*) 'lat',lat (ii_write,jj_write)
        write(6,*) 'lon',lon (ii_write,jj_write)
        write(6,*) 'pbltop',pbltop (ii_write,jj_write)
        write(6,*) 'data diff_t_vert/'    ,diff_t_vert(ii_write,jj_write,:)
        write(6,*) 'data diff_m_vert/'    ,diff_m_vert(ii_write,jj_write,:)
        write(6,*) 'data diff_t/'    ,diff_t(ii_write,jj_write,:)
        write(6,*) 'data diff_m/'    ,diff_m(ii_write,jj_write,:)
  endif

!---------------------------------------------------------------------
!    call edmf_mynn to to calculate tendencies from PBL and convective mass flux
!---------------------------------------------------------------------
  if (do_edmf_mynn .and. do_edmf_mynn_in_physics.eq."down") then

    !--- updated fluxes are not available
    shflx  = 0.
    lhflx  = 0.
    u_flux = 0.
    v_flux = 0.

    rdt_dum = 0.   ! set rdt_dum and tdt_dum to zeros
    tdt_dum = 0. 

    call edmf_mynn_driver ( &
               is, ie, js, je, npz, Time_next, dt, lon, lat, frac_land, area, u_star,  &
               b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, &
               tdt_dum, rdt_dum,  &
               do_edmf_mynn_diagnostic, do_return_edmf_mynn_diff_only, do_edmf_mynn_in_physics, do_tracers_selective, &
               option_edmf2ls_mp, qadt_edmf(is:ie,js:je,:), qldt_edmf(is:ie,js:je,:), qidt_edmf(is:ie,js:je,:), dqa_edmf(is:ie,js:je,:),  dql_edmf(is:ie,js:je,:), dqi_edmf(is:ie,js:je,:), diff_t_vert, diff_m_vert, kpbl_edmf(is:ie,js:je), &
               edmf_mc_full(is:ie,js:je,:), edmf_mc_half(is:ie,js:je,:), edmf_moist_area(is:ie,js:je,:), edmf_dry_area(is:ie,js:je,:), edmf_moist_humidity(is:ie,js:je,:), edmf_dry_humidity(is:ie,js:je,:), &
               z_pbl, udt, vdt, tdt, rdt, tdt_mf_dum, rdiag)

    !--- yhc note, 2021-05-03
    ! diff_t & diff_m from edmf_mynn must be (ie-is+1,je-js+1,npz) dimensions, otherwise it will have grid projection problems to diff_t and diff_m
    ! For example, 
    !   diff_t,m_vert is (ie-is+1,je-js+1,npz) but diff_t,m_edmf is (id,jd,kd) 
    !     diff_t,m_vert (:,:,:)         = diff_t,m_edmf(:,:,:), this will have trouble projecting to correct grids
    !     diff_t,m_vert (is:ie,js:je,:) = diff_t,m_edmf(:,:,:), this will have very obvious problems in projecting to correct grids
    !   The regression 30-day mean will have weird "stripes" over the Pacific and Atlantic.

    !--- replace model diffusion coefficients and z_pbl with edmf_mynn ones when edmf_mynn is not diagnostic
    if (.not.do_edmf_mynn_diagnostic) then
      pbltop(is:ie,js:je) = z_pbl(:,:)
    endif

  if (do_writeout_column) then
        write(6,*) '-------------- i,j,',ii_write,jj_write
        write(6,*) 'lat',lat (ii_write,jj_write)
        write(6,*) 'lon',lon (ii_write,jj_write)
        write(6,*) 'pbltop',pbltop (ii_write,jj_write)
        write(6,*) 'data t_edmf_mynn/'    ,t(ii_write,jj_write,:)
        write(6,*) 'data q_edmf_mynn/'    ,r(ii_write,jj_write,:,nsphum)
        write(6,*) 'data udt_edmf_mynn/'    ,udt(ii_write,jj_write,:)
        write(6,*) 'data vdt_edmf_mynn/'    ,vdt(ii_write,jj_write,:)
        write(6,*) 'data tdt_edmf_mynn/'    ,tdt(ii_write,jj_write,:)
        write(6,*) 'data qdt_edmf_mynn/'    ,rdt(ii_write,jj_write,:,nsphum)
        write(6,*) 'data qadt_edmf_mynn/'    ,rdt(ii_write,jj_write,:,nqa)
        write(6,*) 'data qldt_edmf_mynn/'    ,rdt(ii_write,jj_write,:,nql)
        write(6,*) 'data qidt_edmf_mynn/'    ,rdt(ii_write,jj_write,:,nqi)
        write(6,*) 'data qndt_edmf_mynn/'    ,rdt(ii_write,jj_write,:,nqn)
        write(6,*) 'data diff_t_vert/'    ,diff_t_vert(ii_write,jj_write,:)
        write(6,*) 'data diff_m_vert/'    ,diff_m_vert(ii_write,jj_write,:)
  endif

  endif ! end if of do_edmf_mynn

!-----------------------------------------------------------------------
!    process any tracer fields.
!-----------------------------------------------------------------------

      nextinct = get_tracer_index(MODEL_ATMOS,'Extinction')
      if (Rad_flux_control%do_rad .and. nextinct /= NO_TRACER) then
        rdiag(:,:,:,nextinct) = Rad_flux_block%extinction(:,:,:)
      endif

      call mpp_clock_begin ( tracer_clock )
      call atmos_tracer_driver (is, ie, js, je, Time, lon, lat,  &
                                area, z_pbl, rough_mom,         &
                                frac_open_sea, frac_land, &
                                p_half, p_full,  &
                                u, v, t, r(:,:,:,1), r, rm, rdt, rdiag, dt, &
                                u_star, b_star, q_star, z_half, z_full, &
                                t_surf_rad, albedo, Time_next, &
                                Rad_flux_block%flux_sw_down_vis_dir, &
                                Rad_flux_block%flux_sw_down_vis_dif)
      call mpp_clock_end ( tracer_clock )

!<-- yhc
  if (do_writeout_column) then
        write(6,*) '-------------- i,j,',ii_write,jj_write
        write(6,*) 'lat',lat (ii_write,jj_write)
        write(6,*) 'lon',lon (ii_write,jj_write)
        write(6,*) 'data t_tracer/'    ,t(ii_write,jj_write,:)
        write(6,*) 'data q_tracer/'    ,r(ii_write,jj_write,:,1)
        write(6,*) 'data udt_tracer/'    ,udt(ii_write,jj_write,:)
        write(6,*) 'data vdt_tracer/'    ,vdt(ii_write,jj_write,:)
        write(6,*) 'data tdt_tracer/'    ,tdt(ii_write,jj_write,:)
        write(6,*) 'data qdt_tracer/'    ,rdt(ii_write,jj_write,:,1)
        write(6,*) 'data qadt_tracer/'    ,rdt(ii_write,jj_write,:,nqa)
        write(6,*) 'data qldt_tracer/'    ,rdt(ii_write,jj_write,:,nql)
        write(6,*) 'data qidt_tracer/'    ,rdt(ii_write,jj_write,:,nqi)
        write(6,*) 'data qndt_tracer/'    ,rdt(ii_write,jj_write,:,nqn)
        do rr=1, size(Physics_tendency_block%q_dt,4)
          write(6,*) 'data tr, rr, rdr'    ,rr, rdt(ii_write,jj_write,:,rr)
        enddo
  endif
!--> yhc

!-----------------------------------------------------------------------
!    optionally use an implicit calculation of the vertical diffusion 
!    coefficients.
!
!    the vertical diffusion coefficients are solved using an implicit
!    solution to the following equation:
!
!    dK/dt   = - ( K - K_cur) / tau_diff
!
!    where K         = diffusion coefficient
!          K_cur     = diffusion coefficient diagnosed from current 
!                      time steps' state
!          tau_diff  = time scale for adjustment
!
!    in the code below alpha = dt / tau_diff
!---------------------------------------------------------------------
      if (diffusion_smooth) then
        call get_time (Time_next - Time, sec, day)
        dt2 = real(sec + day*86400)
        alpha = dt2/tau_diff
        diff_m(is:ie,js:je,:) = (diff_m(is:ie,js:je,:) +       &
                                 alpha*(diff_m_vert(:,:,:) +  &
                                 diff_cu_mo(is:ie,js:je,:)) )/&
                                 (1. + alpha)
        where (diff_m(is:ie,js:je,:) < diff_min)
          diff_m(is:ie,js:je,:) = 0.0
        end where
        diff_t(is:ie,js:je,:) = (diff_t(is:ie,js:je,:) +      &
                                 alpha*diff_t_vert(:,:,:) )/  &
                                 (1. + alpha)
        where (diff_t(is:ie,js:je,:) < diff_min)
          diff_t(is:ie,js:je,:) = 0.0
        end where
      else
        diff_t(is:ie,js:je,:) = diff_t_vert
        diff_m(is:ie,js:je,:) = diff_m_vert + diff_cu_mo(is:ie, js:je,:)
      end if

!-----------------------------------------------------------------------
!    call vert_diff_driver_down to calculate the first pass atmos-
!    pheric vertical diffusion.
!-----------------------------------------------------------------------

      !<-- yhc
      !--- save tendencies before vert_diff_driver_down
      udt_before_vdiff_down(:,:,:) = udt(:,:,:)
      vdt_before_vdiff_down(:,:,:) = vdt(:,:,:)
      tdt_before_vdiff_down(is:ie,js:je,:) = tdt(:,:,:)
      rdt_before_vdiff_down(is:ie,js:je,:,:) = rdt(:,:,:,:)
      tau_x_before_vdiff_down(:,:) = tau_x(:,:)
      tau_y_before_vdiff_down(:,:) = tau_y(:,:)
      !--> yhc

      if (id_tdt_phys_vdif_dn > 0) then
        used = send_data ( id_tdt_phys_vdif_dn, -2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_vdif_dn(n) > 0) then
          used = send_data ( id_tracer_phys_vdif_dn(n), -2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

      call mpp_clock_begin ( diff_down_clock )
      radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) - tdt(:,:,:)
      if (do_clubb > 0) then
        call vert_diff_driver_down (is, js, Time_next, dt, p_half,   &
                                    p_full, z_full,   &
                                    diff_m(is:ie,js:je,:),         &
                                    diff_t(is:ie,js:je,:),         &
                                    u ,v ,t ,r(:,:,:,1) ,r(:,:,:,1:ntp), &
                                    dtau_du, dtau_dv, tau_x, tau_y,  &
                                    udt, vdt, tdt, rdt(:,:,:,1), rdt,       &
                                    Surf_diff,                     &
                                    diff_t_clubb=diff_t_clubb(is:ie,js:je,:))

      !<-- yhc 
      elseif (do_edmf_mynn .and. .not.do_edmf_mynn_diagnostic) then

        !--- check tracers
        if (do_tracers_in_edmf_mynn) then     ! do all tracers in edmf_mynn so setting diff_m and diff_t to zeros
          diff_m(is:ie,js:je,:) = 0.          ! vert_diff_driver_down still process variables that would be used in surface-atmosphere coupling
          diff_t(is:ie,js:je,:) = 0.
        elseif (do_edmf_mynn_in_physics.eq."up") then   ! let vert_diff_driver_down process eddy diffusion,
                                                        ! using edmf_mynn diffusion coefficients computed in physics_up
          diff_m(is:ie,js:je,:) = diff_m_edmf(is:ie,js:je,:)
          diff_t(is:ie,js:je,:) = diff_t_edmf(is:ie,js:je,:)
        endif    ! end if of do_tracers_in_edmf_mynn

        call vert_diff_driver_down (is, js, Time_next, dt, p_half,   &
                                    p_full, z_full,   &
                                    diff_m(is:ie,js:je,:),         &  
                                    diff_t(is:ie,js:je,:),         &  
                                    u ,v ,t ,r(:,:,:,1) ,r(:,:,:,1:ntp), &
                                    dtau_du, dtau_dv, tau_x, tau_y,  &
                                    udt, vdt, tdt, rdt(:,:,:,1), rdt,        &
                                    Surf_diff)

  if (do_writeout_column) then
        write(6,*) '-------------- i,j,',ii_write,jj_write
        write(6,*) 'dt,dt2,alpha',dt,dt2,alpha
        write(6,*) 'lat',lat (ii_write,jj_write)
        write(6,*) 'lon',lon (ii_write,jj_write)
        write(6,*) 'tau_x_before_vdiff_down',tau_x_before_vdiff_down(ii_write,jj_write)
        write(6,*) 'tau_x_after_vdiff_down ',tau_x(ii_write,jj_write)
        write(6,*) 'tau_y_before_vdiff_down',tau_y_before_vdiff_down(ii_write,jj_write)
        write(6,*) 'tau_y_after_vdiff_down',tau_y(ii_write,jj_write)
        write(6,*) 'surf_diff%delta_u, after vdiff_down',surf_diff%delta_u(ii_write,jj_write)
        write(6,*) 'surf_diff%delta_v, after vdiff_down',surf_diff%delta_v(ii_write,jj_write)
        write(6,*) 'diff_t',diff_t(ii_write,jj_write,:)
        write(6,*) 'diff_m',diff_m(ii_write,jj_write,:)
        write(6,*) 'udt_before_vdiff_down',udt_before_vdiff_down(ii_write,jj_write,:)
        write(6,*) 'udt_after_vdiff_down ',udt(ii_write,jj_write,:)
        write(6,*) 'rdt_before_vdiff_down, 12',rdt_before_vdiff_down(ii_write,jj_write,:,12)
        write(6,*) 'rdt_after_vdiff_down , 12',rdt(ii_write,jj_write,:,12)
        write(6,*) 'radturbten_after_vdiff_down', radturbten(ii_write,jj_write,:)
  endif

        !--- reset tendencies only the diffusion is handled by edmf_mynn
        if (.not.do_return_edmf_mynn_diff_only) then
          !--- reset tendencies to pre-vert_diff values
          udt(:,:,:) = udt_before_vdiff_down(:,:,:)  ! yhc, even diff_m=0., udt and vdt at the lowest atm loevel are still changed.
          vdt(:,:,:) = vdt_before_vdiff_down(:,:,:)  !      make sure udt and vdt are reset to the values before vdiff_down
          tdt(:,:,:) = tdt_before_vdiff_down(is:ie,js:je,:)
          rdt(:,:,:,nsphum) = rdt_before_vdiff_down(is:ie,js:je,:,nsphum)
          tau_x(:,:) = tau_x_before_vdiff_down(:,:)  !      reset tau_x, tau_y, surf_diff%delta_u, and surf_diff%delta_v
          tau_y(:,:) = tau_y_before_vdiff_down(:,:) 
          surf_diff%delta_u(:,:) = 0.
          surf_diff%delta_v(:,:) = 0.
          radturbten(is:ie,js:je,:) = tdt_rad_only(:,:,:)    ! reset to tdt_rad
        endif

        !--- select which tracers would be updated by vert_diff_driver_down. 
        !    The tracers considered here are qn, qni, and others such as dust, sea salt, etc
        if (do_tracers_in_edmf_mynn) then     ! do all tracers in edmf_mynn so setting all rdt to pre-vert_diff values
          rdt(:,:,:,:)   = rdt_before_vdiff_down(is:ie,js:je,:,:)

        else  ! some tracers are handled by vert_diff
    
          if     (do_tracers_selective.eq.0) then      ! ! do nothing, all tracers are processed by vert_diff excpet specific humidity
            rdt_dum1 = 0.

          elseif     (do_tracers_selective.eq.1) then      ! Reset qn and qni, i.e. vert_diff process all tracers except qn,qni
            if (nqn  > 0) rdt(:,:,:,nqn)  = rdt_before_vdiff_down(is:ie,js:je,:,nqn)
            if (nqni > 0) rdt(:,:,:,nqni) = rdt_before_vdiff_down(is:ie,js:je,:,nqni)
  
          elseif (do_tracers_selective.eq.2) then      ! Reset all tracers, i.e. vert_diff does not process any tracers
            rdt(:,:,:,:)   = rdt_before_vdiff_down(is:ie,js:je,:,:)
  
          elseif (do_tracers_selective.eq.3) then      ! Reset all tracers except qn and qni, i.e. vert_diff only processes qn and qni 
            if (nqn  > 0) rdt_dum1 (:,:,:) = rdt(:,:,:,nqn)
            if (nqni > 0) rdt_dum2 (:,:,:) = rdt(:,:,:,nqni)
  
            rdt(:,:,:,:)    = rdt_before_vdiff_down(is:ie,js:je,:,:)
            if (nqn  > 0) rdt(:,:,:,nqn)  = rdt_dum1(:,:,:)
            if (nqni > 0) rdt(:,:,:,nqni) = rdt_dum2(:,:,:)

          elseif (do_tracers_selective.eq.4) then      ! Reset qa,ql,qi,nqn,nqi, i.e. vert_diff does not process these tracers 
            if (nqa  > 0) rdt(:,:,:,nqa)  = rdt_before_vdiff_down(is:ie,js:je,:,nqa)
            if (nql  > 0) rdt(:,:,:,nql)  = rdt_before_vdiff_down(is:ie,js:je,:,nql)
            if (nqi  > 0) rdt(:,:,:,nqi)  = rdt_before_vdiff_down(is:ie,js:je,:,nqi)
            if (nqn  > 0) rdt(:,:,:,nqn)  = rdt_before_vdiff_down(is:ie,js:je,:,nqn)
            if (nqni > 0) rdt(:,:,:,nqni) = rdt_before_vdiff_down(is:ie,js:je,:,nqni)

          elseif (do_tracers_selective.eq.5) then      ! Reset qa, ql,qi, i.e. vert_diff process the other tracers including qn, qni 
            if (nqa  > 0) rdt(:,:,:,nqa)  = rdt_before_vdiff_down(is:ie,js:je,:,nqa)
            if (nql  > 0) rdt(:,:,:,nql)  = rdt_before_vdiff_down(is:ie,js:je,:,nql)
            if (nqi  > 0) rdt(:,:,:,nqi)  = rdt_before_vdiff_down(is:ie,js:je,:,nqi)

          elseif (do_tracers_selective.eq.6) then      ! process all tracers and qa,ql,qi,qn from EDMF would be computed in mynn_edmf 
            rdt_dum1 = 0.
          
          else   ! do nothing, all tracers are processed by vert_diff excpet specific humidity

           call error_mesg ('physics_driver_mod',  &
                            'do_tracers_selective is not supported', FATAL)

          endif  ! end if of do_tracers_selective
        endif    ! end if of do_tracers_in_edmf_mynn 
        !<-- yhc

      else
        call vert_diff_driver_down (is, js, Time_next, dt, p_half,   &
                                    p_full, z_full,   &
                                    diff_m(is:ie,js:je,:),         &
                                    diff_t(is:ie,js:je,:),         &
                                    u ,v ,t ,r(:,:,:,1) ,r(:,:,:,1:ntp), &
                                    dtau_du, dtau_dv, tau_x, tau_y,  &
                                    udt, vdt, tdt, rdt(:,:,:,1), rdt,        &
                                    Surf_diff)

        !<--- yhc, testing purpose
        !--- select which tracers would be updated by vert_diff_driver_down. 
        !    The tracers considered here are qn, qni, and others such as dust, sea salt, etc
        if     (do_tracers_selective.eq.1) then      ! process all tracers except qn and qni
          if (nqn  > 0) rdt(:,:,:,nqn)  = rdt_before_vdiff_down(is:ie,js:je,:,nqn)
          if (nqni > 0) rdt(:,:,:,nqni) = rdt_before_vdiff_down(is:ie,js:je,:,nqni)

        elseif (do_tracers_selective.eq.2) then      ! do not process any tracers
          rdt(:,:,:,:)   = rdt_before_vdiff_down(is:ie,js:je,:,:)

        elseif (do_tracers_selective.eq.3) then      ! only process qn and qni
          if (nqn  > 0) rdt_dum1 (:,:,:) = rdt(:,:,:,nqn)
          if (nqni > 0) rdt_dum2 (:,:,:) = rdt(:,:,:,nqni)

          rdt(:,:,:,:)    = rdt_before_vdiff_down(is:ie,js:je,:,:)
          if (nqn  > 0) rdt(:,:,:,nqn)  = rdt_dum1(:,:,:)
          if (nqni > 0) rdt(:,:,:,nqni) = rdt_dum2(:,:,:)
        endif
        !---> yhc

      endif

      if (id_tdt_phys_vdif_dn > 0) then
        used = send_data ( id_tdt_phys_vdif_dn, +2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_vdif_dn(n) > 0) then
          used = send_data ( id_tracer_phys_vdif_dn(n), +2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

!<-- yhc
  if (do_writeout_column) then
        write(6,*) '-------------- i,j,',ii_write,jj_write
        write(6,*) 'lat',lat (ii_write,jj_write)
        write(6,*) 'lon',lon (ii_write,jj_write)
        write(6,*) 'data t_physics_down_end/'    ,t(ii_write,jj_write,:)
        write(6,*) 'data q_physics_down_end/'    ,r(ii_write,jj_write,:,1)
        write(6,*) 'data qn_physics_down_end/'    ,r(ii_write,jj_write,:,nqn)
        write(6,*) 'data udt_physics_down_end/'    ,udt(ii_write,jj_write,:)
        write(6,*) 'data vdt_physics_down_end/'    ,vdt(ii_write,jj_write,:)
        write(6,*) 'data tdt_physics_down_end/'    ,tdt(ii_write,jj_write,:)
        write(6,*) 'data qdt_physics_down_end/'    ,rdt(ii_write,jj_write,:,1)
        write(6,*) 'data qadt_physics_down_end/'    ,rdt(ii_write,jj_write,:,nqa)
        write(6,*) 'data qldt_physics_down_end/'    ,rdt(ii_write,jj_write,:,nql)
        write(6,*) 'data qidt_physics_down_end/'    ,rdt(ii_write,jj_write,:,nqi)
        write(6,*) 'data qndt_physics_down_end/'    ,rdt(ii_write,jj_write,:,nqn)
        write(6,*) 'data diff_t_physics_down_end/', diff_t(ii_write,jj_write,:)
        write(6,*) 'data diff_m_physics_down_end/', diff_m(ii_write,jj_write,:)
        write(6,*) 'radturbten_physics_down_end', radturbten(ii_write,jj_write,:)
        write(6,*) 'yhc, nqa, nql, nqi, nqn',nqa, nql, nqi, nqn
        do rr=1, size(Physics_tendency_block%q_dt,4)
          write(6,*) 'data dn, rr, rdr'    ,rr, rdt(ii_write,jj_write,:,rr)
        enddo
  endif
!--> yhc

!<-- yhc
!------- heat diff coeffs used by vdif (units: K/m/s) at full level -------
      if ( id_diff_t_vdif > 0) then
        used = send_data (id_diff_t_vdif, diff_t, Time_next, is, js, 1 )
      endif

!------- momentum diff coeffs used by vdif (units: m2/s) at full level -------
      if ( id_diff_m_vdif > 0) then
        used = send_data (id_diff_m_vdif, diff_m, Time_next, is, js, 1 )
      endif
!--> yhc


!---------------------------------------------------------------------
!    if desired, return diff_m and diff_t to calling routine.
!-----------------------------------------------------------------------
      if (present(difft)) then
        difft = diff_t(is:ie,js:je,:)
      endif
      if (present(diffm)) then
        diffm = diff_m(is:ie,js:je,:)
      endif

     call mpp_clock_end ( diff_down_clock )

      u => null()
      v => null()
      t => null()
      r => null()
      um => null()
      vm => null()
      tm => null()
      rm => null()
      p_full => null()
      p_half => null()
      z_full => null()
      z_half => null()
      udt => null()
      vdt => null()
      tdt => null()
      rdt => null()
      rdiag => null()

 end subroutine physics_driver_down



!#######################################################################
! <SUBROUTINE NAME="physics_driver_up">
!  <OVERVIEW>
!    physics_driver_up completes the calculation of vertical diffusion 
!    and also handles moist physical processes.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_up completes the calculation of vertical diffusion 
!    and also handles moist physical processes.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_up (is, ie, js, je,                    &
!                               Time_prev, Time, Time_next,        &
!                               lat, lon, area,                    &
!                               p_half, p_full, z_half, z_full,    & 
!                               omega,                             &
!                               u, v, t, q, r, um, vm, tm, qm, rm, &
!                               frac_land,                         &
!                               udt, vdt, tdt, qdt, rdt,           &
!                               Surf_diff,                         &
!                               lprec,   fprec, gust  )            &
!  </TEMPLATE>
!  <IN NAME="Time_prev" TYPE="time_type">
!   previous time, for variable um, vm, tm, qm, rm
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   next time, used for diagnostics
!  </IN>
!  <IN NAME="lat" TYPE="real">
!   array of model latitudes at model points [radians]
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   array of model longitudes at model points [radians]
!  </IN>
!  <IN NAME="area" TYPE="real">
!   grid box area - current not used
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   pressure at model interface levels (offset from t,q,u,v,r)
!  </IN>
!  <IN NAME="p_full" TPYE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="z_half" TYPE="real">
!   height at model interface levels
!  </IN>
!  <IN NAME="z_full" TPYE="real">
!   height at full levels
!  </IN>
!  <IN NAME="omega" TYPE="real">
!   Veritical pressure tendency
!  </IN>
!  <IN NAME="u" TYPE="real">
!   zonal wind at current time step
!  </IN>
!  <IN NAME="v" TYPE="real">
!   meridional wind at current time step
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at current time step
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity at current time step
!  </IN>
!  <IN NAME="r" TPYE="real">
!   multiple 3d tracer fields at current time step
!  </IN>
!  <IN NAME="um" TYPE="real">
!   zonal wind at previous time step
!  </IN>
!  <IN NAME="vm" TYPE="real">
!   meridional wind at previous time step
!  </IN>
!  <IN NAME="tm" TYPE="real">
!   temperature at previous time step
!  </IN>
!  <IN NAME="qm" TYPE="real">
!   specific humidity at previous time step
!  </IN>
!  <IN NAME="rm" TPYE="real">
!   multiple 3d tracer fields at previous time step
!  </IN>
!  <IN NAME="frac_land" TYPE="real">
!   fraction of land coverage in a model grid point
!  </IN>
!  <INOUT NAME="udt" TYPE="real">
!   zonal wind tendency
!  </INOUT>
!  <INOUT NAME="vdt" TYPE="real">
!   meridional wind tendency
!  </INOUT>
!  <INOUT NAME="tdt" TYPE="real">
!   temperature tendency
!  </INOUT>
!  <INOUT NAME="qdt" TYPE="real">
!   moisture tracer tendencies
!  </INOUT>
!  <INOUT NAME="rdt" TYPE="real">
!   multiple tracer tendencies
!  </INOUT>
!  <OUT NAME="lprec" TYPE="real">
!  </OUT>
!  <OUT NAME="fprec" TYPE="real">
!  </OUT>
!  <OUT NAME="gust" TYPE="real">
!  </OUT>
!  <INOUT NAME="Surf_diff" TYPE="surface_diffusion_type">
!   Surface diffusion 
!  </INOUT>
! </SUBROUTINE>
!
 subroutine physics_driver_up (is, ie, js, je, npz,        &
                               Time_prev, Time, Time_next, &
                               lat, lon, area,             &
                               Physics_input_block,        &
                               frac_land,                  &
                               u_star, b_star, q_star,     &
                               shflx, lhflx,               &!miz
                               t_ref, q_ref, u_flux, v_flux, & !yhc
                               Physics_tendency_block,     &
                               Moist_clouds_block,         &
                               Cosp_block, Surf_diff,      &
                               lprec, fprec, gust)

!----------------------------------------------------------------------
!    physics_driver_up completes the calculation of vertical diffusion 
!    and also handles moist physical processes.
!---------------------------------------------------------------------

integer,                intent(in)                :: is, ie, js, je, npz
type(time_type),        intent(in)                :: Time_prev, Time, Time_next
real,dimension(:,:),    intent(in)                :: lat, lon, area
type(physics_input_block_type), intent(inout)     :: Physics_input_block
real,dimension(:,:),    intent(in)                :: frac_land
real,dimension(:,:),    intent(in)                :: u_star, b_star, q_star, shflx, lhflx!miz
real,dimension(:,:),    intent(in)                :: t_ref, q_ref, u_flux, v_flux  !yhc 
type(physics_tendency_block_type), intent(inout)  :: Physics_tendency_block
type(clouds_from_moist_block_type), intent(inout) :: Moist_clouds_block
type(cosp_from_rad_block_type),     intent(inout) :: Cosp_block
type(surf_diff_type),   intent(inout)             :: Surf_diff
real,dimension(:,:),    intent(out)               :: lprec, fprec
real,dimension(:,:),    intent(inout)             :: gust

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      npz            number of vertical levels
!      Time_prev      previous time, for variables um,vm,tm,qm,rm 
!                     (time_type)
!      Time           current time, for variables u,v,t,q,r  (time_type)
!      Time_next      next time, used for diagnostics   (time_type)
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      frac_land
!      u_star
!      b_star
!      q_star
!
!  intent(inout) variables:
!
!      Physics_input_block  derived type variable containing: 
!         1) p_half         pressure at half levels (offset from t,q,u,v,r)
!                          [ Pa ]
!         2) p_full         pressure at full levels [ Pa }
!         3) z_half         height at half levels [ m ]
!         4) z_full         height at full levels [ m ]
!         5) u              zonal wind at current time step [ m / s ]
!         6) v              meridional wind at current time step [ m / s ]
!         7) t              temperature at current time step [ deg k ]
!         9) q              multiple 3d tracer fields at current time step
!        10) um,vm          zonal and meridional wind at previous time step
!        11) tm             temperature at previous time step
!        12) qm             multiple 3d tracer fields at previous time step
!        13) omega
!      Physics_tendency_block derived type variable containing:
!          1) u_dt           zonal wind tendency [ m / s**2 ]
!          2) v_dt           meridional wind tendency [ m / s**2 ]
!          3) t_dt           temperature tendency [ deg k / sec ]
!          4) q_dt           multiple tracer tendencies 
!                            (index 1 = specific humidity) 
!                            [ unit / unit / sec ]
!          5) qdiag          multiple 3d diagnostic tracer fields 
!                            [ unit / unit ]
!      Moist_clouds_block
!      Cosp_block
!      Surf_diff      surface_diffusion_type variable
!      gust
!
!   intent(out) variables:
!
!      lprec     
!      fprec       
!
!   intent(in), optional variables:
!
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!   local variables:

      type (precip_flux_type)          :: MP2cosp
      type (Phys2cosp_type)            :: Phys2cosp
      type (phys_mp_exch_type)         :: Phys_mp_exch
      type(aerosol_type)               :: Aerosol
      real, dimension(ie-is+1, je-js+1)          :: gust_cv
      real, dimension(ie-is+1, je-js+1, npz+1)   :: pflux, lphalf
      real, dimension(ie-is+1, je-js+1), target  :: tdt_shf,  qdt_lhf
      integer :: sec, day
      real    :: dt
      integer :: i, j , k, n
      real    :: alphb
      integer :: imax, jmax, kmax
      logical :: used

      type(MP_removal_type) :: Removal_mp
   
!---------------------------------------------------------------------
!   local variables:
!
!        MP2cosp
!        Phys2cosp
!        Phys_mp_exch
!        Aerosol          aerosol_type variable describing the aerosol
!                         fields to be seen by the moist processes routines
!        gust_cv
!        pflux
!        tdt_shf          temperature tendency from sensible heat flux
!        qdt_lhf          moisture tendency from latent heat flux
!        sec, day         second and day components of the time_type 
!                         variable
!        dt               physics time step [ seconds ]
!        i,j,k,n
!        alphb
!        imax,jmax,kmax
!        used
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local pointers to derived type components
!----------------------------------------------------------------------
      real, dimension(:,:,:),   pointer :: t                    
      real, dimension(:,:,:,:), pointer :: r
      real, dimension(:,:,:),   pointer :: p_full, p_half                 
      real, dimension(:,:,:),   pointer :: tdt
      real, dimension(:,:,:,:), pointer :: rdt  
      !<-- yhc        
      real, dimension(ie-is+1,je-js+1,npz)   :: diff_t_vert, diff_m_vert
      real, dimension(:,:,:,:), pointer :: rdiag  ! yhc 
      real, dimension(:,:,:),   pointer :: udt    ! yhc
      real, dimension(:,:,:),   pointer :: vdt    ! yhc
      real :: tt1
      integer :: rr
      real                              :: alpha, dt2
      real, dimension(ie-is+1,je-js+1)       :: z_pbl
      real, dimension( size(Physics_tendency_block%q_dt,1), &
                       size(Physics_tendency_block%q_dt,2), &
                       size(Physics_tendency_block%q_dt,3), &
                       size(Physics_tendency_block%q_dt,4)  )  :: &
         rdt_mynn_ed_am4, rdt_dum
      real, dimension( size(Physics_tendency_block%q_dt,1), &
                       size(Physics_tendency_block%q_dt,2), &
                       size(Physics_tendency_block%q_dt,3)  )  :: &
         tdt_mynn_ed_am4, tdt_mf, qdt_mynn_ed_am4, tdt_dum, qdt_dum
      !--> yhc        

      t => Physics_input_block%t
      r => Physics_input_block%q
      p_full => Physics_input_block%p_full
      p_half => Physics_input_block%p_half
      tdt => Physics_tendency_block%t_dt
      rdt => Physics_tendency_block%q_dt

      udt => Physics_tendency_block%u_dt      ! yhc
      vdt => Physics_tendency_block%v_dt      ! yhc
      rdiag => Physics_tendency_block%qdiag   ! yhc
      tdt_mf = 0. ! yhc, initialize

!<-- yhc
  if (do_writeout_column) then
    write(6,*) '-------------- i,j,',ii_write,jj_write
    write(6,*) 'lat',lat(ii_write,jj_write)  ! radians
    write(6,*) 'lon',lon(ii_write,jj_write)  ! radians
    write(6,*) 'data t_physics_up_begin/'    ,t(ii_write,jj_write,:)
    write(6,*) 'data q_physics_up_begin/'    ,r(ii_write,jj_write,:,1)
    write(6,*) 'data udt_physics_up_begin/'    ,udt(ii_write,jj_write,:)
    write(6,*) 'data vdt_physics_up_begin/'    ,vdt(ii_write,jj_write,:)
    write(6,*) 'data tdt_physics_up_begin/'    ,tdt(ii_write,jj_write,:)
    write(6,*) 'data qdt_physics_up_begin/'    ,rdt(ii_write,jj_write,:,1)
    write(6,*) 'data qadt_physics_up_begin/'    ,rdt(ii_write,jj_write,:,nqa)
    write(6,*) 'data qldt_physics_up_begin/'    ,rdt(ii_write,jj_write,:,nql)
    write(6,*) 'data qidt_physics_up_begin/'    ,rdt(ii_write,jj_write,:,nqi)
  endif
!--> yhc

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
             'module has not been initialized', FATAL)
      endif

!----------------------------------------------------------------------
!    define model spatial dimensions.
!----------------------------------------------------------------------
      imax = ie -is + 1
      jmax = je- js + 1
      kmax = npz        

!-------------------------------------------------------------------------
!    if cosp is activated and this is a step on which cosp input data is
!    to be collected, set up pointers or allocate the necessary derived 
!    type variable components.    
!-------------------------------------------------------------------------
      if (do_cosp) then
        if (step_to_call_cosp) then
          allocate (MP2cosp%fl_lsrain(imax, jmax, kmax))
          allocate (MP2cosp%fl_lssnow(imax, jmax, kmax))
          allocate (MP2cosp%fl_lsgrpl(imax, jmax, kmax))
          allocate (MP2cosp%fl_ccrain(imax, jmax, kmax))
          allocate (MP2cosp%fl_ccsnow(imax, jmax, kmax))
          allocate (MP2cosp%fl_donmca_rain(imax, jmax, kmax))
          allocate (MP2cosp%fl_donmca_snow(imax, jmax, kmax))

          allocate (Phys2cosp%temp_last(imax, jmax, kmax))
          allocate (Phys2cosp%q_last(imax, jmax, kmax))
          Phys2cosp%p_full => Physics_input_block%p_full
          Phys2cosp%z_full => Physics_input_block%z_full
          Phys2cosp%p_half => Physics_input_block%p_half
          Phys2cosp%z_half => Physics_input_block%z_half
          Phys2cosp%u  => Physics_input_block%u
          Phys2cosp%v  => Physics_input_block%v
          allocate (Phys2cosp%frac_land(imax, jmax      ))
          allocate (Phys2cosp%lat   (imax, jmax      ))
         allocate (Phys2cosp%lon   (imax, jmax      ))
        endif
      endif

!----------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
      call get_time (Time_next-Time_prev, sec, day)
      dt = real(sec+day*86400)

!-------------------------------------------------------------------------
!    save temp and moisture tendencies before calculating vertical
!    diffusion. 
!-------------------------------------------------------------------------
      if (id_tdt_phys_vdif_up > 0) then
        used = send_data ( id_tdt_phys_vdif_up, -2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_vdif_up(n) > 0) then
          used = send_data ( id_tracer_phys_vdif_up(n), -2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

!--------------------------------------------------------------------------
!    save temperature and moisture tendencies due to surface fluxes at 
!    lowest-level before calculating vertical diffusion, in the case where
!    these tendencies are not yet to be applied (ie, clubb is active).
!------------------------------------------------------------------------
      if( .not. l_host_applies_sfc_fluxes ) then
          tdt_shf(:,:) = tdt(:, :, kmax)
          qdt_lhf(:,:) = rdt(:, :, kmax, 1)
      endif

      call mpp_clock_begin ( diff_up_clock )
!------------------------------------------------------------------
!    call vert_diff_driver_up to complete the vertical diffusion
!    calculation.
!------------------------------------------------------------------

!<-- yhc
      !--- call vert_diff_driver_up is not doing edmf_mynn
      if (.not.do_edmf_mynn) then    
        call vert_diff_driver_up (is, js, Time_next, dt, p_half,   &
                                  Surf_diff, tdt, rdt(:,:,:,1), rdt )

        tdt_mynn_ed_am4(:,:,:) = tdt(:,:,:) - tdt_before_vdiff_down(is:ie,js:je,:) ! compute tdt_vdif_test

      !--- call vert_diff_driver_up if edmf_mynn is diagnostic
      elseif (do_edmf_mynn .and. do_edmf_mynn_diagnostic) then  
        call vert_diff_driver_up (is, js, Time_next, dt, p_half,   &
                                  Surf_diff, tdt, rdt(:,:,:,1), rdt )

      !--- call vert_diff_driver_up if edmf_mynn only returns diffusion coefficients and diffusion is handle by vert_diff
      elseif (do_edmf_mynn .and. do_return_edmf_mynn_diff_only) then  ! 
        call vert_diff_driver_up (is, js, Time_next, dt, p_half,   &
                                  Surf_diff, tdt, rdt(:,:,:,1), rdt )

      else   ! do edmf_mynn. estimate tdt, qdt from AM4 diffusion solver
        tdt_dum=0.
        qdt_dum=0.
        rdt_dum=0.
        call vert_diff_driver_up (is, js, Time_next, dt, p_half,   &
                                  Surf_diff, tdt_dum, qdt_dum, rdt_dum )
        tdt_mynn_ed_am4(:,:,:) = tdt_dum(:,:,:) - tdt_before_vdiff_down(is:ie,js:je,:) 
        qdt_mynn_ed_am4(:,:,:) = qdt_dum(:,:,:) - rdt_before_vdiff_down(is:ie,js:je,:,nsphum)  
     endif

      !--- get tracer tendencies from vert_diff
      rdt_mynn_ed_am4(:,:,:,:) =  rdt(:,:,:,:) - rdt_before_vdiff_down (is:ie,js:je,:,:)

      !------- liquid water tendency from vert diff (units: kg/kg/s) at full level -------
      if ( id_qldt_vdif > 0) then
        used = send_data (id_qldt_vdif, rdt_mynn_ed_am4(:,:,:,nql), Time_next, is, js, 1 )
      endif

      !------- cloud fraction tendency from vert diff (units: 1/s) at full level -------
      if ( id_qadt_vdif > 0) then
        used = send_data (id_qadt_vdif, rdt_mynn_ed_am4(:,:,:,nqa), Time_next, is, js, 1 )
      endif

      !------- ice water tendency from vert diff (units: kg/kg/s) at full level -------
      if ( id_qidt_vdif > 0) then
        used = send_data (id_qidt_vdif, rdt_mynn_ed_am4(:,:,:,nqi), Time_next, is, js, 1 )
      endif

      !------- spec humid tendency from vert diff (units: kg/kg/s) at full level -------
      if ( id_qdt_vdif_test > 0) then
        used = send_data (id_qdt_vdif_test, rdt_mynn_ed_am4(:,:,:,nsphum), Time_next, is, js, 1 )
      endif

      !------- spec humid tendency from vert diff (units: kg/kg/s) at full level -------
      if ( id_tdt_vdif_test > 0) then
        used = send_data (id_tdt_vdif_test, tdt_mynn_ed_am4(:,:,:), Time_next, is, js, 1 )
      endif

      if (do_writeout_column) then
        write(6,*) 'data qdt/        '    ,rdt(ii_write,jj_write,:,1)
        write(6,*) 'data qdt_before_ed/'  ,rdt_before_vdiff_down(ii_write,jj_write,:,1)
        write(6,*) 'data qdt_mynn_ed/'    ,rdt_mynn_ed_am4(ii_write,jj_write,:,1)
        write(6,*) 'data qadt/        '    ,rdt(ii_write,jj_write,:,nqa)
        write(6,*) 'data qadt_before_ed/'  ,rdt_before_vdiff_down(ii_write,jj_write,:,nqa)
        write(6,*) 'data qadt_mynn_ed/'    ,rdt_mynn_ed_am4(ii_write,jj_write,:,nqa)
        write(6,*) 'data qldt/        '    ,rdt(ii_write,jj_write,:,nql)
        write(6,*) 'data qldt_before_ed/'  ,rdt_before_vdiff_down(ii_write,jj_write,:,nql)
        write(6,*) 'data qldt_mynn_ed/'    ,rdt_mynn_ed_am4(ii_write,jj_write,:,nql)
        write(6,*) 'data qidt/        '    ,rdt(ii_write,jj_write,:,nqi)
        write(6,*) 'data qidt_before_ed/'  ,rdt_before_vdiff_down(ii_write,jj_write,:,nqi)
        write(6,*) 'data qidt_mynn_ed/'    ,rdt_mynn_ed_am4(ii_write,jj_write,:,nqi)
        write(6,*) 'do_tracers_selective, nqa, nql, nqi',do_tracers_selective, nqa, nql, nqi
        write(6,*) 'data tdt_before_vdiff_down/  '    ,tdt_before_vdiff_down(ii_write,jj_write,:)
        write(6,*) 'data tdt_mynn_ed_am4/        '    ,tdt_mynn_ed_am4(ii_write,jj_write,:)
        write(6,*) 'data tdt_dum        /        '    ,tdt_dum(ii_write,jj_write,:)
        write(6,*) 'data qdt_before_vdiff_down/  '    ,rdt_before_vdiff_down(:,:,:,nsphum)
        write(6,*) 'data qdt_mynn_ed_am4/        '    ,qdt_mynn_ed_am4(ii_write,jj_write,:)
        write(6,*) 'data qdt_dum        /        '    ,qdt_dum(ii_write,jj_write,:)
      endif

      if (do_tracers_selective.eq.6) then  ! because ED and MF tendencies would be computed in mynn_edmf, so resetting these
                                           ! to values prior to vert_diff_down
        if (nqa  > 0) rdt(:,:,:,nqa)  = rdt_before_vdiff_down(is:ie,js:je,:,nqa)
        if (nql  > 0) rdt(:,:,:,nql)  = rdt_before_vdiff_down(is:ie,js:je,:,nql)
        if (nqi  > 0) rdt(:,:,:,nqi)  = rdt_before_vdiff_down(is:ie,js:je,:,nqi)
        if (nqn  > 0) rdt(:,:,:,nqn)  = rdt_before_vdiff_down(is:ie,js:je,:,nqn)
        if (nqni > 0) rdt(:,:,:,nqni) = rdt_before_vdiff_down(is:ie,js:je,:,nqni)
      endif
!--> yhc

!<-- yhc
  if (do_writeout_column) then
        write(6,*) '-------------- i,j,',ii_write,jj_write
        write(6,*) 'lat',lat (ii_write,jj_write)
        write(6,*) 'lon',lon (ii_write,jj_write)
        write(6,*) 'data t_diff_up/'    ,t(ii_write,jj_write,:)
        write(6,*) 'data q_diff_up/'    ,r(ii_write,jj_write,:,1)
        write(6,*) 'data qn_diff_up/'    ,r(ii_write,jj_write,:,nqn)
        write(6,*) 'data udt_diff_up/'    ,udt(ii_write,jj_write,:)
        write(6,*) 'data vdt_diff_up/'    ,vdt(ii_write,jj_write,:)
        write(6,*) 'data tdt_diff_up/'    ,tdt(ii_write,jj_write,:)
        write(6,*) 'data qdt_diff_up/'    ,rdt(ii_write,jj_write,:,1)
        write(6,*) 'data qadt_diff_up/'    ,rdt(ii_write,jj_write,:,nqa)
        write(6,*) 'data qldt_diff_up/'    ,rdt(ii_write,jj_write,:,nql)
        write(6,*) 'data qidt_diff_up/'    ,rdt(ii_write,jj_write,:,nqi)
        write(6,*) 'data qndt_diff_up/'    ,rdt(ii_write,jj_write,:,nqn)
        write(6,*) 'yhc, nqa, nql, nqi, nqn',nqa, nql, nqi, nqn
        do rr=1, size(Physics_tendency_block%q_dt,4)
          write(6,*) 'data up, rr, rdr'    ,rr, rdt(ii_write,jj_write,:,rr)
        enddo
  endif
!--> yhc

!      call vert_diff_driver_up (is, js, Time_next, dt, p_half,   &
!                                Surf_diff, tdt, rdt(:,:,:,1), rdt )

!--------------------------------------------------------------------------
!    if the surface tendencies are not to be applied here (ie, clubb),  
!    define those values and remove them from the accumulated time 
!    tendencies. otherwise, set these tendencies to 0.0.
!------------------------------------------------------------------------
      if( .not. l_host_applies_sfc_fluxes ) then
          tdt_shf(:,:) = tdt(:, :, kmax) - tdt_shf(:,:)
          qdt_lhf(:,:) = rdt(:, :, kmax, 1) - qdt_lhf(:,:)

          tdt(:, :, kmax) = tdt(:, :, kmax) - tdt_shf(:,:)
          rdt(:, :, kmax, 1) = rdt(:, :, kmax, 1) - qdt_lhf(:,:)
      endif

      call mpp_clock_end ( diff_up_clock )

!-------------------------------------------------------------------------
!    complete calculation of vertical diffusion tendency diagnostics.
!-------------------------------------------------------------------------
      if (id_tdt_phys_vdif_up > 0) then
        used = send_data ( id_tdt_phys_vdif_up, +2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_vdif_up(n) > 0) then
          used = send_data ( id_tracer_phys_vdif_up(n), +2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

!<-- yhc
!---------------------------------------------------------------------
!    call edmf_mynn to to calculate tendencies from convective mass flux
!---------------------------------------------------------------------

  if (do_edmf_mynn .and. do_edmf_mynn_in_physics.eq."up") then
    call edmf_mynn_driver ( &
               is, ie, js, je, npz, Time_next, dt, lon, lat, frac_land, area, u_star,  &
               b_star, q_star, shflx, lhflx, t_ref, q_ref, u_flux, v_flux, Physics_input_block, &
               tdt_mynn_ed_am4, rdt_mynn_ed_am4,  &
               do_edmf_mynn_diagnostic, do_return_edmf_mynn_diff_only, do_edmf_mynn_in_physics, do_tracers_selective, &
               option_edmf2ls_mp, qadt_edmf(is:ie,js:je,:), qldt_edmf(is:ie,js:je,:), qidt_edmf(is:ie,js:je,:), dqa_edmf(is:ie,js:je,:),  dql_edmf(is:ie,js:je,:), dqi_edmf(is:ie,js:je,:), diff_t_vert, diff_m_vert, kpbl_edmf(is:ie,js:je), &
               edmf_mc_full(is:ie,js:je,:), edmf_mc_half(is:ie,js:je,:), edmf_moist_area(is:ie,js:je,:), edmf_dry_area(is:ie,js:je,:), edmf_moist_humidity(is:ie,js:je,:), edmf_dry_humidity(is:ie,js:je,:), &
               z_pbl, udt, vdt, tdt, rdt, tdt_mf, rdiag)
               !pbltop, udt, vdt, tdt, rdt, rdiag)  ! if using pbltop, amip4 run will fail and such failure happen any time and is not reproducible, yhc 2021-04-21

    !--- yhc note, 2021-05-03
    ! diff_t & diff_m from edmf_mynn must be (ie-is+1,je-js+1,npz) dimensions, otherwise it will have grid projection problems to diff_t and diff_m
    ! For example, 
    !   diff_t,m_vert is (ie-is+1,je-js+1,npz) but diff_t,m_edmf is (id,jd,kd) 
    !     diff_t,m_vert (:,:,:)         = diff_t,m_edmf(:,:,:), this will have trouble projecting to correct grids
    !     diff_t,m_vert (is:ie,js:je,:) = diff_t,m_edmf(:,:,:), this will have very obvious problems in projecting to correct grids
    !   The regression 30-day mean will have weird "stripes" over the Pacific and Atlantic.

    !--- replace pbltop with edmf_mynn ones when edmf_mynn is not diagnostic
    if (.not.do_edmf_mynn_diagnostic) then
      pbltop(is:ie,js:je) = z_pbl(:,:)
    endif

  !--- only modify diff_t and diff_m when tracers are not handled by edmf_mynn and edmf_mynn is not diagnositic purpose
  if (.not.do_tracers_in_edmf_mynn .and. .not.do_edmf_mynn_diagnostic) then

     !--- smooth diffusion coefficients following GFDL diffusion_smooth method
     if (do_edmf_mynn_diffusion_smooth) then
       call get_time (Time_next - Time, sec, day)
       dt2 = real(sec + day*86400)
       alpha = dt2/tau_diff
       diff_m(is:ie,js:je,:) = (diff_m(is:ie,js:je,:) +       &   ! does not consider diff_cu_mo term because this is supposed
                                alpha*diff_m_vert(:,:,:))/&       ! be included in MYNN-EDMF
                                (1. + alpha)
       where (diff_m(is:ie,js:je,:) < diff_min)
         diff_m(is:ie,js:je,:) = 0.0
       end where
       diff_t(is:ie,js:je,:) = (diff_t(is:ie,js:je,:) +      &
                                alpha*diff_t_vert(:,:,:) )/  &
                                (1. + alpha)
       where (diff_t(is:ie,js:je,:) < diff_min)
         diff_t(is:ie,js:je,:) = 0.0
       end where
     else
       diff_t(is:ie,js:je,:) = diff_t_vert(:,:,:)
       diff_m(is:ie,js:je,:) = diff_m_vert(:,:,:)
     end if  ! end if of do_edmf_mynn_diffusion_smooth

     !--- save updated diff_t,m into diff_t,m_edmf
     diff_t_edmf(is:ie,js:je,:) = diff_t(is:ie,js:je,:)
     diff_m_edmf(is:ie,js:je,:) = diff_m(is:ie,js:je,:)
  endif

  if (do_writeout_column) then
        write(6,*) '-------------- i,j,',ii_write,jj_write
        write(6,*) 'lat',lat (ii_write,jj_write)
        write(6,*) 'lon',lon (ii_write,jj_write)
        write(6,*) 'pbltop',pbltop (ii_write,jj_write)
        write(6,*) 'data t_edmf_mynn/'    ,t(ii_write,jj_write,:)
        write(6,*) 'data q_edmf_mynn/'    ,r(ii_write,jj_write,:,nsphum)
        write(6,*) 'data udt_edmf_mynn/'    ,udt(ii_write,jj_write,:)
        write(6,*) 'data vdt_edmf_mynn/'    ,vdt(ii_write,jj_write,:)
        write(6,*) 'data tdt_edmf_mynn/'    ,tdt(ii_write,jj_write,:)
        !write(6,*) 'data radturb_edmf_mynn/'    ,radturbten(ii_write,jj_write,:)
        write(6,*) 'data qdt_edmf_mynn/'    ,rdt(ii_write,jj_write,:,nsphum)
        write(6,*) 'data qadt_edmf_mynn/'    ,rdt(ii_write,jj_write,:,nqa)
        write(6,*) 'data qldt_edmf_mynn/'    ,rdt(ii_write,jj_write,:,nql)
        write(6,*) 'data qidt_edmf_mynn/'    ,rdt(ii_write,jj_write,:,nqi)
        write(6,*) 'data qndt_edmf_mynn/'    ,rdt(ii_write,jj_write,:,nqn)
        write(6,*) 'data diff_t/'    ,diff_t(ii_write,jj_write,:)
        write(6,*) 'data diff_t_vert/'    ,diff_t_vert(ii_write,jj_write,:)
        write(6,*) 'data diff_m_vert/'    ,diff_m_vert(ii_write,jj_write,:)
  endif

  endif  ! end if of do_edmf_mynn
!--> yhc

!-----------------------------------------------------------------------
!    add the temperature tendency due to vertical  diffusion to radturbten.
!-----------------------------------------------------------------------
      ! yhc notes: in AM4, the journey of radturbten
      !    1. In physics_driver_down, radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) + Rad_flux_block%tdt_rad(:,:,:)
      !    2. tdt is accumulated, tdt_accu = tdt_rad + tdt_others
      !    3. Before vert_diff_driver_down, radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) - tdt(:,:,:), 
      !       So, radturbten(is:ie,js:je,:) = -tdt_others
      !    4. After vert_diff_driver_up, tdt_accu = tdt_rad + tdt_others + tdt_vdif
      !    5. Here, radturbten(is:ie,js:je,:) + tdt(:,:,:) = -tdt_others + tdt_accu = tdt_rad + tdt_vdif
      !    6. if tdt includes MF contribution, remove it otherwise there might be double-counting problem 
      !       (MF appears both subsidence and turbulence heating terms)
      if (do_edmf_mynn .and. .not.do_edmf_mynn_diagnostic .and. .not.do_return_edmf_mynn_diff_only) then   ! do_edmf_mynn
        radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) + tdt_mynn_ed_am4(:,:,:)   ! tdt_rad + tdt_vdif 
      else   ! AM4 way
        radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) + tdt(:,:,:) 
      endif

      !------- tempearture tendency from rad and vert diff (units: K/s) at full level -------
      if ( id_tdt_radturb > 0) then
        used = send_data (id_tdt_radturb, radturbten(is:ie,js:je,:), Time_next, is, js, 1 )
      endif

!-----------------------------------------------------------------------
!    prepare to call moist_processes, which calculates moist physics terms,
!    including convection and processes involving condensation, if 
!    desired.
!-----------------------------------------------------------------------
      if (do_moist_processes) then

!-----------------------------------------------------------------------
!    set up diagnostics to capture tendencies due to moist_processes.
!-----------------------------------------------------------------------
        if (id_tdt_phys_moist > 0) then
          used = send_data ( id_tdt_phys_moist, -2.0*tdt(:,:,:), &
                             Time_next, is, js, 1)
        endif

        do n=1,ntp
          if (id_tracer_phys_moist(n) > 0) then
            used = send_data ( id_tracer_phys_moist(n), -2.0*rdt(:,:,:,n), &
                               Time_next, is, js, 1)
          endif
        end do

        call mpp_clock_begin ( moist_processes_clock )

!-----------------------------------------------------------------------
!    call aerosol driver to obtain aerosol data needed in condensation 
!    calculations. if using grey radiation, this data is not needed.
!-----------------------------------------------------------------------
        if (.NOT. do_grey_radiation) then
          pflux(:,:,1) = 0.0E+00
          do k=2,size(p_full,3)
            pflux(:,:,k) = 0.5E+00*(p_full(:,:,k-1) + p_full(:,:,k))
          end do
          pflux(:,:,size(p_full,3)+1) = p_full(:,:,size(p_full,3))
          call aerosol_driver (is, js, Time, r, p_half, pflux, &
                             Aerosol_cld,Aerosol, override_aerosols_cloud)
        endif

!------------------------------------------------------------------------
!   set up pointers to the module variables that are transferred between
!   physics_driver and moist_processes.
!------------------------------------------------------------------------
        Phys_mp_exch%diff_t => diff_t(is:ie,js:je,:)
        Phys_mp_exch%radturbten => radturbten(is:ie,js:je,:)
        Phys_mp_exch%cush       => cush      (is:ie,js:je  )
        Phys_mp_exch%cbmf       => cbmf      (is:ie,js:je  )
        Phys_mp_exch%pbltop     => pbltop    (is:ie,js:je  )
        Phys_mp_exch%diff_cu_mo => diff_cu_mo(is:ie,js:je,:)
        Phys_mp_exch%convect    => convect   (is:ie,js:je  )
        Phys_mp_exch%diff_t_clubb => diff_t_clubb(is:ie,js:je,:)
        Phys_mp_exch%tdt_shf    => tdt_shf 
        Phys_mp_exch%qdt_lhf    => qdt_lhf 
        Phys_mp_exch%hmint      => hmint     (is:ie,js:je  )
        Phys_mp_exch%cgust      => cgust    (is:ie,js:je  )
        Phys_mp_exch%tke        => tke       (is:ie,js:je  )
        Phys_mp_exch%pblhto     => pblhto    (is:ie,js:je  )
        Phys_mp_exch%rkmo       => rkmo      (is:ie,js:je  )
        Phys_mp_exch%taudpo     => taudpo    (is:ie,js:je  )
        Phys_mp_exch%exist_shconv  => exist_shconv (is:ie,js:je,:)
        Phys_mp_exch%exist_dpconv  => exist_dpconv (is:ie,js:je,:)
        Phys_mp_exch%pblht_prev    => pblht_prev   (is:ie,js:je,:)
        Phys_mp_exch%hlsrc_prev    => pblht_prev   (is:ie,js:je,:)
        Phys_mp_exch%qtsrc_prev    => pblht_prev   (is:ie,js:je,:)
        Phys_mp_exch%cape_prev     => pblht_prev   (is:ie,js:je,:)
        Phys_mp_exch%cin_prev      => pblht_prev   (is:ie,js:je,:)
        Phys_mp_exch%tke_prev      => pblht_prev   (is:ie,js:je,:)

        !<-- yhc
        Phys_mp_exch%option_edmf2ls_mp  =>    option_edmf2ls_mp             
        Phys_mp_exch%qadt_edmf      =>    qadt_edmf(is:ie,js:je,:)  
        Phys_mp_exch%qldt_edmf      =>    qldt_edmf(is:ie,js:je,:)  
        Phys_mp_exch%qidt_edmf      =>    qidt_edmf(is:ie,js:je,:)  
        Phys_mp_exch%dqa_edmf       =>    dqa_edmf(is:ie,js:je,:)   
        Phys_mp_exch%dql_edmf       =>    dql_edmf(is:ie,js:je,:)   
        Phys_mp_exch%dqi_edmf       =>    dqi_edmf(is:ie,js:je,:)   
        Phys_mp_exch%kpbl_edmf      =>    kpbl_edmf(is:ie,js:je)   
        Phys_mp_exch%edmf_mc_full         =>    edmf_mc_full        (is:ie,js:je,:)   
        Phys_mp_exch%edmf_mc_half         =>    edmf_mc_half        (is:ie,js:je,:)   
        Phys_mp_exch%edmf_moist_area      =>    edmf_moist_area     (is:ie,js:je,:)   
        Phys_mp_exch%edmf_moist_humidity  =>    edmf_moist_humidity (is:ie,js:je,:)   
        Phys_mp_exch%edmf_dry_area        =>    edmf_dry_area       (is:ie,js:je,:)   
        Phys_mp_exch%edmf_dry_humidity    =>    edmf_dry_humidity   (is:ie,js:je,:)   

        !if (do_edmf_mynn .and. .not.do_edmf_mynn_diagnostic) then
        !  Phys_mp_exch%diff_t => diff_t(is:ie,js:je,:)
        !endif
        !--> yhc

!-----------------------------------------------------------------------
!    call moist processes to compute moist physics, including convection 
!    and processes involving condenstion.
!-----------------------------------------------------------------------
        call moist_processes (    &
              is, ie, js, je, npz, Time_next, dt, frac_land, u_star,  &
              b_star, q_star, area, lon, lat, Physics_input_block,   &
              Moist_clouds_block, Physics_tendency_block, Phys_mp_exch, &
              Surf_diff, Removal_mp, shflx, lhflx,  &
              lprec, fprec, gust_cv, Aerosol=Aerosol)
        call mpp_clock_end ( moist_processes_clock )

!<-- yhc
  if (do_writeout_column) then
    write(6,*) '-------------- i,j,',ii_write,jj_write
    write(6,*) 'lat',lat(ii_write,jj_write)  ! radians
    write(6,*) 'lon',lon(ii_write,jj_write)  ! radians
        write(6,*) 'data t_moist_up/'    ,t(ii_write,jj_write,:)
        write(6,*) 'data q_moist_up/'    ,r(ii_write,jj_write,:,1)
        write(6,*) 'data qn_moist_up/'    ,r(ii_write,jj_write,:,nqn)
        write(6,*) 'data tdt_moist_up/'    ,tdt(ii_write,jj_write,:)
        write(6,*) 'data qdt_moist_up/'    ,rdt(ii_write,jj_write,:,nsphum)
        write(6,*) 'data qadt_moist_up/'    ,rdt(ii_write,jj_write,:,nqa)
        write(6,*) 'data qldt_moist_up/'    ,rdt(ii_write,jj_write,:,nql)
        write(6,*) 'data qidt_moist_up/'    ,rdt(ii_write,jj_write,:,nqi)
        write(6,*) 'data qndt_moist_up/'    ,rdt(ii_write,jj_write,:,nqn)
        write(6,*) 'Phys_mp_exch%diff_t', Phys_mp_exch%diff_t(ii_write,jj_write,:)
        write(6,*) 'radturbten_moist_up', radturbten(ii_write,jj_write,:)
        write(6,*) 'tdt_mynn_ed_am4', tdt_mynn_ed_am4(ii_write,jj_write,:)
  endif
!--> yhc

!-------------------------------------------------------------------------
!    save the cumulus momentum output field in a module variable for use
!    in vert diff calculation done in physics_driver_down. reinitialize
!    radturbten for use on next time step.
!-------------------------------------------------------------------------
        radturbten(is:ie,js:je,:) = 0.0

!---------------------------------------------------------------------
!    add the convective gustiness effect to that previously obtained 
!    from non-convective parameterizations.
!---------------------------------------------------------------------
        gust = sqrt( gust*gust + gust_cv*gust_cv)

!------------------------------------------------------------------------
!    complete calculation of moist processes tendency diagnostics.
!------------------------------------------------------------------------
        if (id_tdt_phys_moist > 0) then
          used = send_data ( id_tdt_phys_moist, +2.0*tdt(:,:,:), &
                             Time_next, is, js, 1)
        endif
        do n=1,ntp
          if (id_tracer_phys_moist(n) > 0) then
            used = send_data ( id_tracer_phys_moist(n), +2.0*rdt(:,:,:,n), &
                               Time_next, is, js, 1)
          endif
        end do

!------------------------------------------------------------------------
!    calculate temperature and tracer tendencies due to model physics.
!------------------------------------------------------------------------
        if (id_tdt_phys > 0) then
           used = send_data ( id_tdt_phys, tdt(:,:,:), &
                              Time_next, is, js, 1)
        endif
        do n=1,ntp
          if (id_tracer_phys(n) > 0) then
            used = send_data ( id_tracer_phys(n), rdt(:,:,:,n), &
                               Time_next, is, js, 1)
          endif
        end do

!----------------------------------------------------------------------
!    if the Aerosol derived type variable component arrays were allocated, 
!    call aerosol_dealloc to deallocate them.
!----------------------------------------------------------------------
        if (.not. do_grey_radiation) call aerosol_dealloc (Aerosol)

      !------ CMIP diagnostics (tendencies due to physics) ------
      if (query_cmip_diag_id(ID_tntmp) .or. query_cmip_diag_id(ID_tnhusmp)) then
        lphalf = log(p_half)
      endif
      if (query_cmip_diag_id(ID_tntmp)) then
        used = send_cmip_data_3d (ID_tntmp, tdt(:,:,:), Time_next, is, js, 1, phalf=lphalf)
      endif
      if (query_cmip_diag_id(ID_tnhusmp)) then
        used = send_cmip_data_3d (ID_tnhusmp, rdt(:,:,:,nsphum), Time_next, is, js, 1, phalf=lphalf)
      endif


!-----------------------------------------------------------------------
!    code needed to execute COSP
!-----------------------------------------------------------------------
        if (do_cosp) then
          call mpp_clock_begin ( cosp_clock )
          alphb = SUM(temp_last(is:ie,js:je,:))

!---------------------------------------------------------------------
!    on the first step of a job segment, the values of t,q and precip 
!    flux will not be available at the proper time level. in this case
!    denoted by temp-_last = 0.0, use values from the current step for 
!    t, q and precip flux.
!---------------------------------------------------------------------
          if (alphb == 0.) then
            call define_cosp_precip_fluxes (is, js, Precip_flux,   &
                                                               Removal_mp)
            if (step_to_call_cosp) then
              Phys2cosp%temp_last(:,:,:) = t(:,:,:) + dt*tdt(:,:,:)
              Phys2cosp%q_last(:,:,:) = r(:,:,:,1) + dt*rdt(:,:,:,1)
              MP2cosp%fl_lsrain(:,:,:) =  &
                                       Precip_flux%fl_lsrain(is:ie,js:je,:)
              MP2cosp%fl_lssnow(:,:,:) =   &
                                       Precip_flux%fl_lssnow(is:ie,js:je,:)
              MP2cosp%fl_lsgrpl(:,:,:) =    &
                                       Precip_flux%fl_lsgrpl(is:ie,js:je,:)
              MP2cosp%fl_ccrain(:,:,:) =    &
                                       Precip_flux%fl_ccrain(is:ie,js:je,:)
              MP2cosp%fl_ccsnow(:,:,:) =     &
                                       Precip_flux%fl_ccsnow(is:ie,js:je,:)
              MP2cosp%fl_donmca_rain(:,:,:) =    &
                                  Precip_flux%fl_donmca_rain(is:ie,js:je,:)
              MP2cosp%fl_donmca_snow(:,:,:) =    &
                                  Precip_flux%fl_donmca_snow(is:ie,js:je,:)
            endif
          else

!--------------------------------------------------------------------
!    on all other steps of the job on which the cosp simulator is 
!    called, define input variables needed by COSP from values computed on
!    the last step that are currently available, before calculating new
!    values for the current step.
!--------------------------------------------------------------------
            if (step_to_call_cosp) then
              MP2cosp%fl_lsrain(:,:,:) =    &
                                    Precip_flux%fl_lsrain(is:ie,js:je,:)
              MP2cosp%fl_lssnow(:,:,:) =    &
                                     Precip_flux%fl_lssnow(is:ie,js:je,:)
              MP2cosp%fl_lsgrpl(:,:,:) =    &
                                     Precip_flux%fl_lsgrpl(is:ie,js:je,:)
              MP2cosp%fl_ccrain(:,:,:) =    &
                                     Precip_flux%fl_ccrain(is:ie,js:je,:)
              MP2cosp%fl_ccsnow(:,:,:) =    &
                                     Precip_flux%fl_ccsnow(is:ie,js:je,:)
              MP2cosp%fl_donmca_rain(:,:,:) =  &
                                Precip_flux%fl_donmca_rain(is:ie,js:je,:)
              MP2cosp%fl_donmca_snow(:,:,:) =  &
                                Precip_flux%fl_donmca_snow(is:ie,js:je,:)
              Phys2cosp%temp_last(:,:,:) = temp_last(is:ie,js:je,:)
              Phys2cosp%q_last(:,:,:)    = q_last(is:ie,js:je,:)
            endif
            call define_cosp_precip_fluxes (is, js, Precip_flux,   &
                                                               Removal_mp)
          endif

          if (step_to_call_cosp) then
!----------------------------------------------------------------------
!    define the remaining input fields needed by cosp_driver.
!----------------------------------------------------------------------
            Phys2cosp%lat    = lat*180./ACOS(-1.0)
            Phys2cosp%lon    = lon*180./ACOS(-1.0)
            Phys2cosp%frac_land = frac_land

            call cosp_driver (   &
             is, ie, js, je, Time_next, MP2cosp, Phys2cosp, Cosp_block)

!-----------------------------------------------------------------------
!    deallocate arrays used to hold cosp input fields.
!-----------------------------------------------------------------------
            deallocate (MP2cosp%fl_lsrain)
            deallocate (MP2cosp%fl_lssnow)
            deallocate (MP2cosp%fl_lsgrpl)
            deallocate (MP2cosp%fl_ccrain)
            deallocate (MP2cosp%fl_ccsnow)
            deallocate (MP2cosp%fl_donmca_rain)
            deallocate (MP2cosp%fl_donmca_snow)
 
            deallocate (Phys2cosp%temp_last)
            deallocate (Phys2cosp%q_last)
            Phys2cosp%p_full => null()
            Phys2cosp%z_full => null()
            Phys2cosp%p_half => null()
            Phys2cosp%z_half => null()
            Phys2cosp%u => null()
            Phys2cosp%v => null()
            deallocate (Phys2cosp%frac_land)
            deallocate (Phys2cosp%lat   )
            deallocate (Phys2cosp%lon   )
          endif ! (step_to_call_cosp)

!--------------------------------------------------------------------
!    save t and q from end of step for use with next call to COSP
!--------------------------------------------------------------------
          temp_last(is:ie,js:je,:) = t(:,:,:) + tdt(:,:,:)*dt
          q_last   (is:ie,js:je,:) = r(:,:,:,1) + rdt(:,:,:,1)*dt
          call mpp_clock_end ( cosp_clock )
        else
            if (allocated(Removal_mp%ice_precflxh)) then
                deallocate(Removal_mp%ice_precflxh)
            endif
            if (allocated(Removal_mp%liq_precflxh)) then
                deallocate(Removal_mp%liq_precflxh)
            endif
            if (allocated(Removal_mp%frz_mesoh)) then
                deallocate(Removal_mp%frz_mesoh)
            endif
            if (allocated(Removal_mp%liq_mesoh)) then
                deallocate(Removal_mp%liq_mesoh)
            endif
            if (allocated(Removal_mp%frz_cellh)) then
                deallocate(Removal_mp%frz_cellh)
            endif
            if (allocated(Removal_mp%liq_cellh)) then
                deallocate(Removal_mp%liq_cellh)
            endif
            if (allocated(Removal_mp%mca_frzh)) then
                deallocate(Removal_mp%mca_frzh)
            endif
            if (allocated(Removal_mp%mca_liqh)) then
                deallocate(Removal_mp%mca_liqh)
            endif
            if (allocated(Removal_mp%rain3d)) then
                deallocate(Removal_mp%rain3d)
            endif
            if (allocated(Removal_mp%snowclr3d)) then
                deallocate(Removal_mp%snowclr3d)
            endif
        endif ! (do_cosp)
      endif ! do_moist_processes
       
!<-- yhc
  if (do_writeout_column) then
    write(6,*) '-------------- i,j,',ii_write,jj_write
    write(6,*) 'lat',lat(ii_write,jj_write)  ! radians
    write(6,*) 'lon',lon(ii_write,jj_write)  ! radians
    write(6,*) 'data t_physics_up_end/'    ,t(ii_write,jj_write,:)
    write(6,*) 'data q_physics_up_end/'    ,r(ii_write,jj_write,:,1)
    write(6,*) 'data udt_physics_up_end/'    ,udt(ii_write,jj_write,:)
    write(6,*) 'data vdt_physics_up_end/'    ,vdt(ii_write,jj_write,:)
    write(6,*) 'data tdt_physics_up_end/'    ,tdt(ii_write,jj_write,:)
    write(6,*) 'data qdt_physics_up_end/'    ,rdt(ii_write,jj_write,:,1)
    write(6,*) 'data qadt_physics_up_end/'    ,rdt(ii_write,jj_write,:,nqa)
    write(6,*) 'data qldt_physics_up_end/'    ,rdt(ii_write,jj_write,:,nql)
    write(6,*) 'data qidt_physics_up_end/'    ,rdt(ii_write,jj_write,:,nqi)
  endif
!--> yhc


!<-- yhc
  if (do_stop_run) then
    call error_mesg('physics_driver_up',  &
                    'stop by yihsuan', FATAL)
  endif
!--> yhc

!<-- debug 
!  tt1 = tdt_max / 86400.  ! change unit from K/day to K/sec
!  do i=1,size(tdt,1)
!  do j=1,size(tdt,2)
!  do k=1,size(tdt,3)
!    if ( abs(tdt(i,j,k)) .ge. tt1 ) then
!      write(6,*) 'phys, >tdt_max,i,j,lat,lon,',tdt_max,i,j,lat(i,j),lon(i,j)
!
!      if (do_limit_tdt) then
!        if (tdt(i,j,k).ge.0.) then
!          tdt(i,j,k) = tdt_limit / 86400.
!        else
!          tdt(i,j,k) = -1.*tdt_limit / 86400.
!        endif
!      endif
!    endif
!  enddo
!  enddo
!  enddo
!-->




!-----------------------------------------------------------------------
!    nullify all local pointers.
!-----------------------------------------------------------------------

!write(6,*) 'a1.0, nQke, rdiag(:,:,:,nQke)', rdiag(:,:,:,28)

      t => null()
      r => null()
      p_full => null()
      p_half => null()
      tdt => null()
      rdt => null()
      udt => null()     ! yhc
      vdt => null()     ! yhc
      rdiag => null()   ! yhc

      Phys_mp_exch%diff_t => null()
      Phys_mp_exch%radturbten => null()
      Phys_mp_exch%diff_cu_mo => null()
      Phys_mp_exch%diff_t_clubb => null()
      Phys_mp_exch%cush   => null()
      Phys_mp_exch%cbmf   => null()
      Phys_mp_exch%pbltop => null()
      Phys_mp_exch%convect => null()
      Phys_mp_exch%tdt_shf => null()
      Phys_mp_exch%qdt_lhf => null()
      Phys_mp_exch%hmint   => null()
      Phys_mp_exch%cgust   => null()
      Phys_mp_exch%tke     => null()
      Phys_mp_exch%pblhto  => null()
      Phys_mp_exch%rkmo    => null()
      Phys_mp_exch%taudpo  => null()
      Phys_mp_exch%exist_shconv => null()
      Phys_mp_exch%exist_dpconv => null()
      Phys_mp_exch%pblht_prev   => null()
      Phys_mp_exch%hlsrc_prev   => null()
      Phys_mp_exch%qtsrc_prev   => null()
      Phys_mp_exch%cape_prev   => null()
      Phys_mp_exch%cin_prev   => null()
      Phys_mp_exch%tke_prev   => null()

      !<--- yhc
      Phys_mp_exch%option_edmf2ls_mp    => null()  
      Phys_mp_exch%qadt_edmf            => null()  
      Phys_mp_exch%qldt_edmf            => null()  
      Phys_mp_exch%qidt_edmf            => null()  
      Phys_mp_exch%dqa_edmf             => null()  
      Phys_mp_exch%dql_edmf             => null()  
      Phys_mp_exch%dqi_edmf             => null()  
      Phys_mp_exch%kpbl_edmf            => null()  
      Phys_mp_exch%edmf_mc_full         => null()
      Phys_mp_exch%edmf_mc_half         => null()
      Phys_mp_exch%edmf_moist_area      => null()
      Phys_mp_exch%edmf_moist_humidity  => null()
      Phys_mp_exch%edmf_dry_area        => null()
      Phys_mp_exch%edmf_dry_humidity    => null()
      !---> yhc
!-----------------------------------------------------------------------



 end subroutine physics_driver_up


!#######################################################################
! <SUBROUTINE NAME="physics_driver_end">
!  <OVERVIEW>
!   physics_driver_end is the destructor for physics_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_end is the destructor for physics_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_end (Time)
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
! </SUBROUTINE>
!
subroutine physics_driver_end (Time, Physics, Moist_clouds,  &
                               Physics_tendency, Atm_block)

!---------------------------------------------------------------------
!    physics_driver_end is the destructor for physics_driver_mod.
!---------------------------------------------------------------------

type(time_type), intent(in) :: Time
type(physics_type), intent(in) :: Physics
type(clouds_from_moist_type), intent(inout) :: Moist_clouds(:)
type(physics_tendency_type),  intent(inout) :: Physics_tendency
type(block_control_type), intent(in) :: Atm_block

!--------------------------------------------------------------------
!   intent(in) variables:
! 
!      Time      current time [ time_type(days, seconds) ]
!
!--------------------------------------------------------------------
integer :: n, nb, nc, ibs, ibe, jbs, jbe
integer :: moist_processes_term_clock, damping_term_clock, turb_term_clock, &
           diff_term_clock, aerosol_term_clock, clubb_term_clock, &
           edmf_mynn_term_clock, & ! yhc
           grey_radiation_term_clock, tracer_term_clock, cosp_term_clock

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
              'module has not been initialized', FATAL)
      endif

      clubb_term_clock =      &
        mpp_clock_id( '   Phys_driver_term: clubb: Termination', &
                grain=CLOCK_MODULE_DRIVER )
      moist_processes_term_clock =      &
        mpp_clock_id( '   Phys_driver_term: MP: Termination', &
                grain=CLOCK_MODULE_DRIVER )
      damping_term_clock         =     &
        mpp_clock_id( '   Phys_driver_term: Damping: Termination',    &
                  grain=CLOCK_MODULE_DRIVER )
      turb_term_clock            =      &
        mpp_clock_id( '   Phys_driver_term: Vert. Turb.: Termination', &
                  grain=CLOCK_MODULE_DRIVER )
      diff_term_clock       =     &
        mpp_clock_id( '   Phys_driver_term: Vert. Diff.: Termination',   &
                 grain=CLOCK_MODULE_DRIVER )
      cosp_term_clock       =       &
        mpp_clock_id( '   Phys_driver_term: COSP: Termination', &
                       grain=CLOCK_MODULE_DRIVER )
      edmf_mynn_term_clock       =       &
        mpp_clock_id( '   Phys_driver_term: edmf_mynn: Termination', &
                       grain=CLOCK_MODULE_DRIVER )  ! yhc

      if (do_moist_processes) &
      aerosol_term_clock       =       &
        mpp_clock_id( '   Phys_driver_term: Aerosol: Termination', &
                       grain=CLOCK_MODULE_DRIVER )
      grey_radiation_term_clock       =       &
        mpp_clock_id( '   Phys_driver_term: Grey Radiation: Termination', &
                       grain=CLOCK_MODULE_DRIVER )
      tracer_term_clock          =      &
        mpp_clock_id( '   Phys_driver_term: Tracer: Termination',    &
                 grain=CLOCK_MODULE_DRIVER )

      do nb = 1, Atm_block%nblks
        ibs = Atm_block%ibs(nb)-Atm_block%isc+1
        ibe = Atm_block%ibe(nb)-Atm_block%isc+1
        jbs = Atm_block%jbs(nb)-Atm_block%jsc+1
        jbe = Atm_block%jbe(nb)-Atm_block%jsc+1

        do nc = 1, size(Moist_clouds(1)%block(nb)%Cloud_data,1)

          ! common to all cloud schemes
          Restart%Cloud_data(nc)%cloud_area    (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%cloud_area
          Restart%Cloud_data(nc)%liquid_amt    (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%liquid_amt
          Restart%Cloud_data(nc)%ice_amt       (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_amt
          Restart%Cloud_data(nc)%droplet_number(ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%droplet_number
          Restart%Cloud_data(nc)%scheme_name                       = Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name

          ! properties specific to large-scale/stratiform clouds
          if (trim(Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'strat_cloud') then
            Restart%Cloud_data(nc)%ice_number(ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_number
            Restart%Cloud_data(nc)%rain      (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%rain
            Restart%Cloud_data(nc)%snow      (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%snow
            Restart%Cloud_data(nc)%rain_size (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%rain_size
            Restart%Cloud_data(nc)%snow_size (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%snow_size
          endif
 
          ! properties specific to donner deep clouds (both cell and meso)
          if (trim(Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'donner_cell' .or. &
              trim(Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'donner_meso') then
            Restart%Cloud_data(nc)%liquid_size(ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%liquid_size
            Restart%Cloud_data(nc)%ice_size   (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_size
            Restart%Cloud_data(nc)%nsum_out   (ibs:ibe,jbs:jbe)   = Moist_clouds(1)%block(nb)%Cloud_data(nc)%nsum_out
          endif

          ! properties specific to uw shallow convective clouds
          if (trim(Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'uw_conv') then
            Restart%Cloud_data(nc)%ice_number(ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_number
          endif

        enddo
      enddo

!----------------------------------------------------------------------
!    call physics_driver_netcdf to process physics_driver restart file.
!----------------------------------------------------------------------
      call physics_driver_netcdf

!--------------------------------------------------------------------
!    call the destructor routines for those modules who were initial-
!    ized from this module.
!--------------------------------------------------------------------
      call mpp_clock_begin ( turb_term_clock )
      call vert_turb_driver_end
      call mpp_clock_end ( turb_term_clock )
      call mpp_clock_begin ( diff_term_clock )
      call vert_diff_driver_end
      call mpp_clock_end ( diff_term_clock )

!<-- yhc
      call mpp_clock_begin ( edmf_mynn_term_clock )
      call edmf_mynn_end  
      call mpp_clock_end ( edmf_mynn_term_clock )
!--> yhc

!--------------------------------------------------------------------
!    terminate radiation routines and data
!--------------------------------------------------------------------
      if (do_moist_processes) then
        call mpp_clock_begin ( aerosol_term_clock )
        call aerosol_end (Aerosol_cld)
        call mpp_clock_end ( aerosol_term_clock )
      endif

      call mpp_clock_begin ( grey_radiation_term_clock )
      if(do_grey_radiation) call grey_radiation_end 
      call mpp_clock_end ( grey_radiation_term_clock )

      if (do_moist_processes) then  
        call mpp_clock_begin ( moist_processes_term_clock )
        call moist_processes_end ()
        call mpp_clock_end ( moist_processes_term_clock )
      endif

      call mpp_clock_begin ( tracer_term_clock )
      call atmos_tracer_driver_end
      call mpp_clock_end ( tracer_term_clock )
      call mpp_clock_begin ( damping_term_clock )
      call damping_driver_end
      call mpp_clock_end ( damping_term_clock )
      if (do_cosp) then
        call mpp_clock_begin ( cosp_term_clock )
        call cosp_driver_end
        call mpp_clock_end ( cosp_term_clock )
      endif

!---------------------------------------------------------------------
!    deallocate the module variables.
!---------------------------------------------------------------------
      deallocate (diff_cu_mo, diff_t, diff_m, pbltop, cush, cbmf,  &
                  hmint, cgust, tke, pblhto, rkmo, taudpo, exist_shconv, &  ! h1g, 2017-01-31
                  exist_dpconv, & 
                  pblht_prev, hlsrc_prev, qtsrc_prev, cape_prev, cin_prev, tke_prev, & !h1g, 2017-01-31
                  option_edmf2ls_mp, qadt_edmf, qldt_edmf, qidt_edmf, diff_t_edmf, diff_m_edmf, kpbl_edmf, &  !yhc
                  dqa_edmf, dql_edmf, dqi_edmf, & ! yhc
                  edmf_mc_full, edmf_mc_half, edmf_moist_area, edmf_dry_area, edmf_moist_humidity, edmf_dry_humidity, & ! yhc
                  rdt_before_vdiff_down, tdt_before_vdiff_down, &  ! yhc
                  convect, radturbten, r_convect)

      if (do_cosp) then
        deallocate ( temp_last, q_last)
        deallocate (&
           Precip_flux%fl_lsrain, Precip_flux%fl_lssnow,   &
           Precip_flux%fl_lsgrpl, Precip_flux%fl_ccrain,   &
           Precip_flux%fl_ccsnow, Precip_flux%fl_donmca_rain,    &
           Precip_flux%fl_donmca_snow)
      endif
 
      deallocate ( diff_t_clubb )
      
      deallocate (id_tracer_phys_vdif_dn)
      deallocate (id_tracer_phys_vdif_up)
      deallocate (id_tracer_phys_turb)
      deallocate (id_tracer_phys_moist)

      call dealloc_physics_tendency_type (Physics_tendency)

      deallocate ( Restart%Cloud_data           )   ! h1g, 2017-02-02

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.


!-----------------------------------------------------------------------

 end subroutine physics_driver_end

!#######################################################################
! <SUBROUTINE NAME="physics_driver_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine physics_driver_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp


  if (mpp_pe() == mpp_root_pe() ) then
     call error_mesg('physics_driver_mod', 'Writing netCDF formatted restart file: RESTART/physics_driver.res.nc', NOTE)
  endif
  call physics_driver_netcdf(timestamp)
  call vert_turb_driver_restart(timestamp)

  call moist_processes_restart(timestamp)
  call damping_driver_restart(timestamp)

end subroutine physics_driver_restart
! </SUBROUTINE> NAME="physics_driver_restart"

! <SUBROUTINE NAME="physics_driver_netcdf">
!
! <DESCRIPTION>
! Write out restart file for physics driver.
! This routine is needed so that physics_driver_restart and physics_driver_end
! can call a routine which will not result in multiple copies of restart files 
! being written by the destructor routines.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine physics_driver_netcdf(timestamp)
  character(len=*), intent(in), optional :: timestamp

    r_convect = 0.
    where(convect)
       r_convect = 1.0
    end where
    call save_restart(Phy_restart, timestamp)
    if(in_different_file) call save_restart(Til_restart, timestamp)

end subroutine physics_driver_netcdf
! </SUBROUTINE> NAME="physics_driver_netcdf"

!#######################################################################
! <FUNCTION NAME="do_moist_in_phys_up">
!  <OVERVIEW>
!    do_moist_in_phys_up returns the value of do_moist_processes
!  </OVERVIEW>
!  <DESCRIPTION>
!    do_moist_in_phys_up returns the value of do_moist_processes
!  </DESCRIPTION>
!  <TEMPLATE>
!   logical = do_moist_in_phys_up()
!  </TEMPLATE>
! </FUNCTION>
!
function do_moist_in_phys_up()

!--------------------------------------------------------------------
!    do_moist_in_phys_up returns the value of do_moist_processes
!----------------------------------------------------------------------

logical :: do_moist_in_phys_up

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('do_moist_in_phys_up',  &
              'module has not been initialized', FATAL)
      endif
 
!-------------------------------------------------------------------
!    define output variable.
!-------------------------------------------------------------------
      do_moist_in_phys_up = do_moist_processes

 
end function do_moist_in_phys_up

!#####################################################################
! <FUNCTION NAME="get_diff_t">
!  <OVERVIEW>
!    returns the values of array diff_t
!  </OVERVIEW>
!  <DESCRIPTION>
!    returns the values of array diff_t
!  </DESCRIPTION>
!  <TEMPLATE>
!   diff_t(:,:,:) = get_diff_t()
!  </TEMPLATE>
! </FUNCTION>
!
!#####################################################################
function get_diff_t() result(diff_t_out)
real, dimension(size(diff_t,1),size(diff_t,2),size(diff_t,3)) :: diff_t_out

  if ( .not. module_is_initialized) then
    call error_mesg ('get_diff_t','module has not been initialized', FATAL)
  endif

  diff_t_out = diff_t

end function get_diff_t

!#####################################################################
! <FUNCTION NAME="get_radturbten">
!  <OVERVIEW>
!    returns the values of array radturbten
!  </OVERVIEW>
!  <DESCRIPTION>
!    returns the values of array radturbten
!  </DESCRIPTION>
!  <TEMPLATE>
!   radturbten(:,:,:) = get_radturbten()
!  </TEMPLATE>
! </FUNCTION>
!
!#####################################################################
function get_radturbten() result(radturbten_out)
real, dimension(size(radturbten,1),size(radturbten,2),size(radturbten,3)) :: radturbten_out

  if ( .not. module_is_initialized) then
    call error_mesg ('get_radturbten','module has not been initialized', FATAL)
  endif

  radturbten_out = radturbten

end function get_radturbten
!#####################################################################
! <SUBROUTINE NAME="zero_radturbten">
!  <OVERVIEW>
!    sets all values of array radturbten to zero
!  </OVERVIEW>
!  <DESCRIPTION>
!    sets all values of array radturbten to zero
!  </DESCRIPTION>
!  <TEMPLATE>
!   call zero_radturbten()
!  </TEMPLATE>
! </SUBROUTINE>
!
!#####################################################################
subroutine zero_radturbten()

  if ( .not. module_is_initialized) then
    call error_mesg ('zero_radturbten','module has not been initialized', FATAL)
  endif

  radturbten = 0.0

end subroutine zero_radturbten



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
               
     
!#####################################################################
! <SUBROUTINE NAME="physics_driver_register_restart">
!  <OVERVIEW>
!    physics_driver_register_restart will register restart field when do_netcdf file 
!    is true. 
!  </OVERVIEW>
subroutine physics_driver_register_restart (Restart)
  type(clouds_from_moist_block_type), intent(inout), target :: Restart
  character(len=64) :: fname, fname2
  integer           :: id_restart
  integer           :: nc
  logical           :: reproduce_ulm_restart = .true.
  integer           :: index_strat

  if (do_moist_processes) then  
    if(doing_prog_clouds) then 
       now_doing_strat = 1
    else
       now_doing_strat = 0
    endif

    if(doing_edt) then 
       now_doing_edt = 1
    else
       now_doing_edt = 0
    endif

    if(doing_entrain) then 
       now_doing_entrain = 1
    else
       now_doing_entrain = 0
    endif
  endif

  fname = 'physics_driver.res.nc'
  call get_mosaic_tile_file(fname, fname2, .false. ) 
  allocate(Phy_restart)
  if(trim(fname2) == trim(fname)) then
     Til_restart => Phy_restart
     in_different_file = .false.
  else
     in_different_file = .true.
     allocate(Til_restart)
  endif

  id_restart = register_restart_field(Phy_restart, fname, 'vers',          vers,              no_domain=.true.)
  id_restart = register_restart_field(Phy_restart, fname, 'doing_strat',   now_doing_strat,   no_domain=.true.)
  id_restart = register_restart_field(Phy_restart, fname, 'doing_edt',     now_doing_edt,     no_domain=.true.)
  id_restart = register_restart_field(Phy_restart, fname, 'doing_entrain', now_doing_entrain, no_domain=.true.)

  id_restart = register_restart_field(Til_restart, fname, 'diff_cu_mo', diff_cu_mo)
  id_restart = register_restart_field(Til_restart, fname, 'pbltop',     pbltop)
  id_restart = register_restart_field(Til_restart, fname, 'cush',       cush, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'cbmf',       cbmf, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'hmint',      hmint, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'cgust',      cgust, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'tke',        tke, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'pblhto',     pblhto, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'rkmo',       rkmo, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'taudpo',     taudpo, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'exist_shconv', exist_shconv, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'exist_dpconv', exist_dpconv, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'pblht_prev',   pblht_prev, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'hlsrc_prev',   hlsrc_prev, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'qtsrc_prev',   qtsrc_prev, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'cape_prev',    cape_prev, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'cin_prev',     cin_prev, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'tke_prev',     tke_prev, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'diff_t',     diff_t)
  id_restart = register_restart_field(Til_restart, fname, 'diff_m',     diff_m)
  id_restart = register_restart_field(Til_restart, fname, 'convect',    r_convect) 

  !<-- yhc
  id_restart = register_restart_field(Til_restart, fname, 'option_edmf2ls_mp', option_edmf2ls_mp   ) 
  id_restart = register_restart_field(Til_restart, fname, 'qadt_edmf', qadt_edmf  ) 
  id_restart = register_restart_field(Til_restart, fname, 'qldt_edmf', qldt_edmf   ) 
  id_restart = register_restart_field(Til_restart, fname, 'qidt_edmf', qidt_edmf   ) 
  id_restart = register_restart_field(Til_restart, fname, 'diff_t_edmf', diff_t_edmf   ) 
  id_restart = register_restart_field(Til_restart, fname, 'diff_m_edmf', diff_m_edmf      ) 
  id_restart = register_restart_field(Til_restart, fname, 'kpbl_edmf', kpbl_edmf      ) 
  id_restart = register_restart_field(Til_restart, fname, 'dqa_edmf', dqa_edmf      ) 
  id_restart = register_restart_field(Til_restart, fname, 'dql_edmf', dql_edmf      ) 
  id_restart = register_restart_field(Til_restart, fname, 'dqi_edmf', dqi_edmf      ) 
  id_restart = register_restart_field(Til_restart, fname, 'rdt_before_vdiff_down', rdt_before_vdiff_down      ) 
  id_restart = register_restart_field(Til_restart, fname, 'edmf_mc_full', edmf_mc_full      ) 
  id_restart = register_restart_field(Til_restart, fname, 'edmf_mc_half', edmf_mc_half      ) 
  id_restart = register_restart_field(Til_restart, fname, 'edmf_moist_area', edmf_moist_area      ) 
  id_restart = register_restart_field(Til_restart, fname, 'edmf_moist_humidity', edmf_moist_humidity      ) 
  id_restart = register_restart_field(Til_restart, fname, 'edmf_dry_area', edmf_dry_area      ) 
  id_restart = register_restart_field(Til_restart, fname, 'edmf_dry_humidity', edmf_dry_humidity      ) 
  !--> yhc

  if (do_clubb > 0) then
    id_restart = register_restart_field(Til_restart, fname, 'diff_t_clubb', diff_t_clubb, mandatory = .false.)
  end if
  if (doing_prog_clouds) then
    id_restart = register_restart_field(Til_restart, fname, 'radturbten',       radturbten)
  endif

  index_strat = 0
  do nc = 1, size(Restart%Cloud_data,1)
    if (trim(Restart%Cloud_data(nc)%scheme_name).eq.'strat_cloud' .and. .not. reproduce_ulm_restart) then
      id_restart = register_restart_field(Til_restart, fname, 'lsc_cloud_area',     Restart%Cloud_data(nc)%cloud_area,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_liquid',         Restart%Cloud_data(nc)%liquid_amt,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_ice',            Restart%Cloud_data(nc)%ice_amt,        mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_droplet_number', Restart%Cloud_data(nc)%droplet_number, mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_ice_number',     Restart%Cloud_data(nc)%ice_number,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_snow',           Restart%Cloud_data(nc)%snow,           mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_rain',           Restart%Cloud_data(nc)%rain,           mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_snow_size',      Restart%Cloud_data(nc)%snow_size,      mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_rain_size',      Restart%Cloud_data(nc)%rain_size,      mandatory = .false.)
    endif
    if (trim(Restart%Cloud_data(nc)%scheme_name).eq.'strat_cloud' .and. reproduce_ulm_restart) index_strat = nc

    if (trim(Restart%Cloud_data(nc)%scheme_name).eq.'donner_cell') then
      id_restart = register_restart_field(Til_restart, fname, 'cell_cloud_frac',  Restart%Cloud_data(nc)%cloud_area,  mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'cell_liquid_amt',  Restart%Cloud_data(nc)%liquid_amt,  mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'cell_liquid_size', Restart%Cloud_data(nc)%liquid_size, mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'cell_ice_amt',     Restart%Cloud_data(nc)%ice_amt,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'cell_ice_size',    Restart%Cloud_data(nc)%ice_size,    mandatory = .false.)
    endif

    if (trim(Restart%Cloud_data(nc)%scheme_name).eq.'donner_meso') then
      id_restart = register_restart_field(Til_restart, fname, 'meso_cloud_frac',  Restart%Cloud_data(nc)%cloud_area,  mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'meso_liquid_amt',  Restart%Cloud_data(nc)%liquid_amt,  mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'meso_liquid_size', Restart%Cloud_data(nc)%liquid_size, mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'meso_ice_amt',     Restart%Cloud_data(nc)%ice_amt,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'meso_ice_size',    Restart%Cloud_data(nc)%ice_size,    mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'nsum',             Restart%Cloud_data(nc)%nsum_out,    mandatory = .false.)
    endif

    if (trim(Restart%Cloud_data(nc)%scheme_name).eq.'uw_conv') then
      id_restart = register_restart_field(Til_restart, fname, 'shallow_cloud_area',     Restart%Cloud_data(nc)%cloud_area,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'shallow_liquid',         Restart%Cloud_data(nc)%liquid_amt,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'shallow_ice',            Restart%Cloud_data(nc)%ice_amt,        mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'shallow_droplet_number', Restart%Cloud_data(nc)%droplet_number, mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'shallow_ice_number',     Restart%Cloud_data(nc)%ice_number,     mandatory = .false.)
    endif
  enddo

    ! save large-scale clouds last to reproduce ulm code
    if (index_strat > 0) then
      nc = index_strat
      id_restart = register_restart_field(Til_restart, fname, 'lsc_cloud_area',     Restart%Cloud_data(nc)%cloud_area,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_liquid',         Restart%Cloud_data(nc)%liquid_amt,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_ice',            Restart%Cloud_data(nc)%ice_amt,        mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_droplet_number', Restart%Cloud_data(nc)%droplet_number, mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_ice_number',     Restart%Cloud_data(nc)%ice_number,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_snow',           Restart%Cloud_data(nc)%snow,           mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_rain',           Restart%Cloud_data(nc)%rain,           mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_snow_size',      Restart%Cloud_data(nc)%snow_size,      mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_rain_size',      Restart%Cloud_data(nc)%rain_size,      mandatory = .false.)
    endif

end subroutine physics_driver_register_restart
! </SUBROUTINE>    
!#####################################################################
! <SUBROUTINE NAME="check_args">
!  <OVERVIEW>
!    check_args determines if the input arrays to physics_driver_down
!    are of a consistent size.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_args determines if the input arrays to physics_driver_down
!    are of a consistent size.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call check_args (lat, lon, area, p_half, p_full, z_half, z_full,&
!                        u, v, t, q, r, um, vm, tm, qm, rm,             &
!                        udt, vdt, tdt, qdt, rdt)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   array of model latitudes at model points [radians]
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   array of model longitudes at model points [radians]
!  </IN>
!  <IN NAME="area" TYPE="real">
!   grid box area - current not used
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   pressure at model interface levels (offset from t,q,u,v,r)
!  </IN>
!  <IN NAME="p_full" TPYE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="z_half" TYPE="real">
!   height at model interface levels
!  </IN>
!  <IN NAME="z_full" TPYE="real">
!   height at full levels
!  </IN>
!  <IN NAME="u" TYPE="real">
!   zonal wind at current time step
!  </IN>
!  <IN NAME="v" TYPE="real">
!   meridional wind at current time step
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at current time step
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity at current time step
!  </IN>
!  <IN NAME="r" TPYE="real">
!   multiple 3d tracer fields at current time step
!  </IN>
!  <IN NAME="um" TYPE="real">
!   zonal wind at previous time step
!  </IN>
!  <IN NAME="vm" TYPE="real">
!   meridional wind at previous time step
!  </IN>
!  <IN NAME="tm" TYPE="real">
!   temperature at previous time step
!  </IN>
!  <IN NAME="qm" TYPE="real">
!   specific humidity at previous time step
!  </IN>
!  <IN NAME="rm" TPYE="real">
!   multiple 3d tracer fields at previous time step
!  </IN>
!  <IN NAME="udt" TYPE="real">
!   zonal wind tendency
!  </IN>
!  <IN NAME="vdt" TYPE="real">
!   meridional wind tendency
!  </IN>
!  <IN NAME="tdt" TYPE="real">
!   temperature tendency
!  </IN>
!  <IN NAME="qdt" TYPE="real">
!   moisture tracer tendencies
!  </IN>
!  <IN NAME="rdt" TYPE="real">
!   multiple tracer tendencies
!  </IN>
! </SUBROUTINE>
!
subroutine check_args (lat, lon, area, p_half, p_full, z_half, z_full,&
                        u, v, t, q, r, um, vm, tm, qm, rm,             &
                        udt, vdt, tdt, qdt, rdt, rdiag)

!----------------------------------------------------------------------
!    check_args determines if the input arrays to physics_driver_down
!    are of a consistent size.
!-----------------------------------------------------------------------

real,    dimension(:,:),    intent(in)           :: lat, lon, area
real,    dimension(:,:,:),  intent(in)           :: p_half, p_full,   &
                                                    z_half, z_full,   &
                                                    u, v, t, q, um, vm, &
                                                    tm, qm
real,    dimension(:,:,:,:),intent(in)           :: r, rm
real,    dimension(:,:,:),  intent(in)           :: udt, vdt, tdt, qdt
real,    dimension(:,:,:,:),intent(in)           :: rdt
real,    dimension(:,:,:,ntp+1:),intent(in)      :: rdiag

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      p_half         pressure at half levels (offset from t,q,u,v,r)
!                     [ Pa ]
!      p_full         pressure at full levels [ Pa }
!      z_half         height at half levels [ m ]
!      z_full         height at full levels [ m ]
!      u              zonal wind at current time step [ m / s ]
!      v              meridional wind at current time step [ m / s ]
!      t              temperature at current time step [ deg k ]
!      q              specific humidity at current time step  kg / kg ]
!      r              multiple 3d tracer fields at current time step
!      um,vm          zonal and meridional wind at previous time step
!      tm,qm          temperature and specific humidity at previous 
!                     time step
!      rm             multiple 3d tracer fields at previous time step
!      udt            zonal wind tendency [ m / s**2 ]
!      vdt            meridional wind tendency [ m / s**2 ]
!      tdt            temperature tendency [ deg k / sec ]
!      qdt            specific humidity tendency 
!                     [  kg vapor / kg air / sec ]
!      rdt            multiple tracer tendencies [ unit / unit / sec ]
!
!   intent(in), optional:
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      integer ::  id, jd, kd  ! model dimensions on the processor  
      integer ::  ierr        ! error flag

!--------------------------------------------------------------------
!    define the sizes that the arrays should be.
!--------------------------------------------------------------------
      id = size(u,1) 
      jd = size(u,2) 
      kd = size(u,3) 

!--------------------------------------------------------------------
!    check the dimensions of each input array. if they are incompat-
!    ible in size with the standard, the error flag is set to so
!    indicate.
!--------------------------------------------------------------------
      ierr = 0
      ierr = ierr + check_dim (lat, 'lat',  id,jd)
      ierr = ierr + check_dim (lon, 'lon',  id,jd)
      ierr = ierr + check_dim (area,'area', id,jd)

      ierr = ierr + check_dim (p_half,'p_half', id,jd,kd+1)
      ierr = ierr + check_dim (p_full,'p_full', id,jd,kd)
      ierr = ierr + check_dim (z_half,'z_half', id,jd,kd+1)
      ierr = ierr + check_dim (z_full,'z_full', id,jd,kd)

      ierr = ierr + check_dim (u, 'u',  id,jd,kd)
      ierr = ierr + check_dim (v, 'v',  id,jd,kd)
      ierr = ierr + check_dim (t, 't',  id,jd,kd)
      ierr = ierr + check_dim (q, 'q',  id,jd,kd)
      ierr = ierr + check_dim (um,'um', id,jd,kd)
      ierr = ierr + check_dim (vm,'vm', id,jd,kd)
      ierr = ierr + check_dim (tm,'tm', id,jd,kd)
      ierr = ierr + check_dim (qm,'qm', id,jd,kd)

      ierr = ierr + check_dim (udt,'udt', id,jd,kd)
      ierr = ierr + check_dim (vdt,'vdt', id,jd,kd)
      ierr = ierr + check_dim (tdt,'tdt', id,jd,kd)
      ierr = ierr + check_dim (qdt,'qdt', id,jd,kd)

      if (ntp > 0) then
        ierr = ierr + check_dim (r,  'r',   id,jd,kd,ntp)
        ierr = ierr + check_dim (rm, 'rm',  id,jd,kd,ntp)
      endif
      if (ntp > 0) then
        ierr = ierr + check_dim (rdt,'rdt', id,jd,kd,ntp)
      endif
      if (nt > ntp) then
        ierr = ierr + check_dim (rdiag,'rdiag', id,jd,kd,nt-ntp)
      endif

!--------------------------------------------------------------------
!    if any problems were detected, exit with an error message.
!--------------------------------------------------------------------
      if (ierr > 0) then
        call error_mesg ('physics_driver_mod', 'bad dimensions', FATAL)
      endif

!-----------------------------------------------------------------------


      end subroutine check_args


!#######################################################################
! <FUNCTION NAME="check_dim_2d">
!  <OVERVIEW>
!    check_dim_2d compares the size of two-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_dim_2d compares the size of two-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </DESCRIPTION>
!  <TEMPLATE>
!    check_dim_2d (data,name,id,jd) result (ierr)
!  </TEMPLATE>
!  <IN NAME="data" TYPE="real">
!   array of data to be checked
!  </IN>
!  <IN NAME="name" TYPE="character">
!   name associated with array to be checked
!  </IN>
!  <IN NAME="id, jd" TYPE="integer">
!   expected i and j dimensions
!  </IN>
! </FUNCTION>
!
function check_dim_2d (data,name,id,jd) result (ierr)

!--------------------------------------------------------------------
!    check_dim_2d compares the size of two-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------

real,    intent(in), dimension(:,:) :: data
character(len=*), intent(in)        :: name
integer, intent(in)                 :: id, jd
integer                             :: ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data        array to be checked
!     name        name associated with array to be checked
!     id, jd      expected i and j dimensions
!     
!  result variable:
!
!     ierr        set to 0 if ok, otherwise is a count of the number
!                 of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
             'dimension 1 of argument ' //  &
              name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
           call error_mesg ('physics_driver_mod',  &
                'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
           ierr = ierr + 1
      endif

!----------------------------------------------------------------------

      end function check_dim_2d

!#######################################################################
! <FUNCTION NAME="check_dim_3d">
!  <OVERVIEW>
!    check_dim_3d compares the size of three-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_dim_3d compares the size of three-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </DESCRIPTION>
!  <TEMPLATE>
!    check_dim_3d (data,name,id,jd, kd) result (ierr)
!  </TEMPLATE>
!  <IN NAME="data" TYPE="real">
!   array of data to be checked
!  </IN>
!  <IN NAME="name" TYPE="character">
!   name associated with array to be checked
!  </IN>
!  <IN NAME="id, jd, kd" TYPE="integer">
!   expected i, j and k dimensions
!  </IN>
! </FUNCTION>
!
function check_dim_3d (data,name,id,jd,kd) result (ierr)

!--------------------------------------------------------------------
!    check_dim_3d compares the size of thr1eedimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------

real,    intent(in), dimension(:,:,:) :: data
character(len=*), intent(in)          :: name
integer, intent(in)                   :: id, jd, kd
integer  ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data        array to be checked
!     name        name associated with array to be checked
!     id, jd,kd   expected i, j and k dimensions
!     
!  result variable:
!
!     ierr        set to 0 if ok, otherwise is a count of the number
!                 of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 1 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
        call error_mesg ('physics_driver_mod',  &
              'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,3) /= kd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 3 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif

!---------------------------------------------------------------------


      end function check_dim_3d


!#######################################################################
! <FUNCTION NAME="check_dim_4d">
!  <OVERVIEW>
!    check_dim_4d compares the size of four-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_dim_4d compares the size of four-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </DESCRIPTION>
!  <TEMPLATE>
!    check_dim_4d (data,name,id,jd, kd, nt) result (ierr)
!  </TEMPLATE>
!  <IN NAME="data" TYPE="real">
!   array of data to be checked
!  </IN>
!  <IN NAME="name" TYPE="character">
!   name associated with array to be checked
!  </IN>
!  <IN NAME="id, jd, kd, nt" TYPE="integer">
!   expected i, j, k and 4th dimensions
!  </IN>
! </FUNCTION>
!
function check_dim_4d (data,name,id,jd,kd,nt) result (ierr)

!--------------------------------------------------------------------
!    check_dim_4d compares the size of four dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------
real,    intent(in), dimension(:,:,:,:) :: data
character(len=*), intent(in)            :: name
integer, intent(in)                     :: id, jd, kd, nt
integer                                 :: ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data          array to be checked
!     name          name associated with array to be checked
!     id,jd,kd,nt   expected i, j and k dimensions
!     
!  result variable:
!
!     ierr          set to 0 if ok, otherwise is a count of the number
!                   of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 1 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,3) /= kd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 3 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,4) /= nt) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 4 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif

!---------------------------------------------------------------------


      end function check_dim_4d

!<--- yhc, 2022-07-15
!#######################################################################
  subroutine compute_vert_adv_tend_offline (dts, temp, temp_vert_advec_scheme, tracer_vert_advec_scheme)

    use vert_advection_mod,       only: vert_advection, SECOND_CENTERED, &
                                        FOURTH_CENTERED, FINITE_VOLUME_LINEAR, &
                                        FINITE_VOLUME_PARABOLIC, &
                                        SECOND_CENTERED_WTS, FOURTH_CENTERED_WTS, &
                                        ADVECTIVE_FORM

    use             fv_pack, only: nlev
    use      constants_mod, only:  rdgas, rvgas, cp_air, hlv, kappa, grav, pi, SECONDS_PER_DAY

    !------------------
    ! purpose    
    !------------------
    ! Given necessary information (set by input_profile), 
    ! use AM4 vert_advection module to compute offline vertical advection of temperaure and specific humidity
    ! This can be used to decompose full dynamical tendencies (*dt_dyn) into horizontal and vertical components

    !------------------
    ! input argument    
    !------------------
    real   , intent(in) :: dts                        ! time step (sec)
    real   , intent(in) :: temp(:,:,:)                ! temperature, for getting dimension
    integer, intent(in) :: temp_vert_advec_scheme     ! Which advection scheme should be used, =4 in SCM
    integer, intent(in) :: tracer_vert_advec_scheme   ! Which advection scheme should be used, =4 in SCM  

    !------------------
    ! local argument
    !------------------
    integer, parameter :: nprof = 10  ! number of profiles

    real, dimension(size(temp,1),size(temp,2),size(temp,3)+1) :: &
      ph,        &   !  pressure at half levels (Pa)
      omega_h        !  vertical pressure velocity at half levels (Pa/s)

    real, dimension(size(temp,1),size(temp,2),size(temp,3))   :: &
      dT_vadv,   &   !  temperature vertical advection (K/s)
      dT_adi,    &   !  temperature tendencies due to adiabatic heating (K/s)
      tdt_vadv,  &   !  sum of dT_vadv and dT_adi (K/s)
      dqv_vadv,  &   !  specific humidity vertical advection (kg/kg/s)
      dp,        &   !  pressure thickness between half levels (Pa)
      qq,        &   !  specific humidity at full levels (K)
      omega_f,   &   !  vertical pressure velocity at full levels (Pa/s)
      pf,        &   !  pressure at full levels (Pa)
      pt             !  temperature at full levels (K)

    real, dimension(nprof, size(temp,1),size(temp,2),size(temp,3)+1) :: &
      ph_np             !  pressure at half levels (Pa)

    real, dimension(nprof, size(temp,1),size(temp,2),size(temp,3))   :: &
      dT_vadv_np,   &   !  temperature vertical advection (K/s)
      dqv_vadv_np,  &   !  specific humidity vertical advection (kg/kg/s)
      dp_np,        &   !  pressure thickness between half levels (Pa)
      qq_np,        &   !  specific humidity at full levels (K)
      omega_f_np,   &   !  vertical pressure velocity at full levels (Pa/s)
      pt_np             !  temperature at full levels (K)
        
    integer i,j,k,kdim, n,ndim
    integer ii_write, jj_write
    integer :: vadvec_scheme
 
    character*100 input_profile_longname
    !----------------------

    !--- initialization
    dT_vadv=0.; dT_adi=0.; tdt_vadv=0.; dqv_vadv=0. ; dp=0.

    ii_write = 1
    jj_write = 1
    kdim = size(temp,3)
    ndim = nprof

    !--- determine input profile

    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    if (input_profile .eq. 0) then
      input_profile_longname = "test111"
      ndim=1
      ph_np      (1, ii_write, jj_write, :) = (/100.000000000000,    400.000000000000, &
   818.602110000000,    1378.88653000000,    2091.79519000000, &
   2983.64084000000,    4121.78960000000,    5579.22148000000, &
   7429.32203000000,    9739.83459000000,    12573.1869600000, &
   15988.0458500000,    20044.5383400000,    24794.7845300000, &
   30283.9819600000,    36517.8148800000,    43391.6043100000, &
   50609.5508700000,    57745.5530700000,    64437.0976500000, &
   70498.0509400000,    75875.3518300000,    80577.5427300000, &
   84639.2099500000,    88112.5613700000,    91059.5443200000, &
   93543.5709900000,    95617.8850400000,    97335.0224200000, &
   98736.9161300000,    99856.8968600000,    100728.562400000, &
   101362.435000000,    101780.000000000/)
      omega_f_np (1, ii_write, jj_write, :) = (/0.000000000000000E+000,  0.000000000000000E+000, &
  0.000000000000000E+000,  0.000000000000000E+000,  0.000000000000000E+000, &
  0.000000000000000E+000,  0.000000000000000E+000,  0.000000000000000E+000, &
  0.000000000000000E+000,  0.000000000000000E+000,  0.000000000000000E+000, &
  0.000000000000000E+000,  0.000000000000000E+000,  0.000000000000000E+000, &
  0.000000000000000E+000,  0.000000000000000E+000,  0.000000000000000E+000, &
  0.000000000000000E+000,  0.000000000000000E+000,  0.000000000000000E+000, &
  0.000000000000000E+000,  0.000000000000000E+000,  2.930938536955648E-002, &
  5.812430974580787E-002,  4.692672232055620E-002,  3.784775569538005E-002, &
  3.007107786950254E-002,  2.271821830454597E-002,  1.649746359500637E-002, &
  1.135228502210560E-002,  7.203072036325996E-003,  3.950445515424451E-003, &
  1.572365337671026E-003/)
      pt_np      (1, ii_write, jj_write, :) = (/200.000000000000,    200.000000000000, &
   200.000000000000,    200.000000000000,    200.000000000000, &
   200.000000000000,    200.000000000000,    200.000000000000, &
   201.014876573338,    211.606846860821,    221.720769700700, &
   231.308545253451,    240.342448573245,    248.804323922616, &
   256.670990366677,    263.878216874829,    270.303508236814, &
   275.817264667225,    280.378898054894,    284.067885838342, &
   287.028330245440,    289.405274346612,    291.318346199574, &
   293.170646722972,    294.305579368352,    283.586238287798, &
   284.602464182821,    286.051031894629,    287.365367815623, &
   288.417108332689,    289.240938553385,    289.859807312425, &
   290.290044380658/)
      qq_np      (1, ii_write, jj_write, :) = (/    0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01  ,  0.1000E-01/) 

    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    elseif (input_profile .eq. 1) then
      input_profile_longname = "nudgeAM4_2001Jul10_DYCOMS_6_gridpnts"
      ndim=6
     
      n=1
      ph_np      (n, ii_write, jj_write, :) = (/1.00000000e+02, 4.00000000e+02, 8.18602100e+02, 1.37888600e+03,&
       2.09179500e+03, 2.98364100e+03, 4.12179000e+03, 5.57922200e+03,&
       7.43014597e+03, 9.74300006e+03, 1.25800972e+04, 1.60000637e+04,&
       2.00630357e+04, 2.48211577e+04, 3.03196612e+04, 3.65642110e+04,&
       4.34499208e+04, 5.06804537e+04, 5.78289458e+04, 6.45322276e+04,&
       7.06038303e+04, 7.59905891e+04, 8.07010555e+04, 8.47698762e+04,&
       8.82493460e+04, 9.12015231e+04, 9.36899277e+04, 9.57679001e+04,&
       9.74880640e+04, 9.88924266e+04, 1.00014365e+05, 1.00887542e+05,&
       1.01522486e+05, 1.01940734e+05/)
      omega_f_np (n, ii_write, jj_write, :) = (/5.5864159e-05,  4.1863739e-05,  6.9493122e-05,  2.1614606e-04, &
        6.1794353e-04,  1.4780449e-03,  2.1667634e-03,  1.1605999e-03, &
        5.1287045e-03,  2.4840644e-02,  4.7111008e-02,  4.2171512e-02, &
        8.3119757e-03, -2.8849840e-02, -2.8391883e-02,  1.3964850e-02, &
        4.7824223e-02,  6.1862223e-02,  6.8485878e-02,  7.0093594e-02,&
        8.0621913e-02,  9.0339243e-02,  8.8166632e-02,  8.3732627e-02,&
        8.0987364e-02,  7.3665336e-02,  5.9194673e-02,  4.2879146e-02,&
        2.9343367e-02,  1.8671820e-02,  1.0217882e-02,  3.8292645e-03,&
       -6.0854829e-04 /)
      pt_np      (n, ii_write, jj_write, :) = (/259.45935, 239.50665, 231.75102, 227.15317, 223.42767, 219.33342,&
       214.53229, 210.01433, 205.73952, 205.2344 , 211.90053, 219.18068,&
       223.57361, 230.57608, 240.46982, 251.60287, 261.12836, 267.81866,&
       274.59195, 279.92966, 283.33636, 286.89807, 289.14258, 291.8757 ,&
       291.96075, 285.6379 , 284.8002 , 286.30807, 287.6244 , 288.68655,&
       289.52362, 290.15817, 290.61783/)
      qq_np      (n, ii_write, jj_write, :) = (/1.7761962e-06, 1.8233893e-06, 1.8778217e-06, 1.8789498e-06,&
       1.8576874e-06, 1.8248562e-06, 1.6796107e-06, 1.6547541e-06,&
       2.3708410e-06, 2.3808147e-06, 2.3180116e-06, 4.1625153e-06,&
       1.1502249e-05, 5.4811208e-05, 3.3545142e-04, 8.4202422e-04,&
       1.3463127e-03, 1.8920604e-03, 1.9157311e-03, 3.0488023e-03,&
       4.7540176e-03, 5.0785756e-03, 4.9187848e-03, 4.9800370e-03,&
       3.1521651e-03, 4.8433845e-03, 8.1914850e-03, 8.6888187e-03,&
       8.7844329e-03, 8.8426853e-03, 8.8921972e-03, 8.9459913e-03,&
       9.0630408e-03/)

      n=2
      ph_np      (n, ii_write, jj_write, :) = (/1.00000000e+02, 4.00000000e+02, 8.18602100e+02, 1.37888600e+03,&
       2.09179500e+03, 2.98364100e+03, 4.12179000e+03, 5.57922200e+03,&
       7.42988666e+03, 9.74200479e+03, 1.25779242e+04, 1.59962843e+04,&
       2.00572188e+04, 2.48128640e+04, 3.03084408e+04, 3.65496206e+04,&
       4.34315819e+04, 5.06581565e+04, 5.78027211e+04, 6.45023115e+04,&
       7.05705654e+04, 7.59543500e+04, 8.06622138e+04, 8.47287851e+04,&
       8.82063306e+04, 9.11568745e+04, 9.36439022e+04, 9.57207242e+04,&
       9.74399363e+04, 9.88435225e+04, 9.99648454e+04, 1.00837547e+05,&
       1.01472154e+05, 1.01890188e+05/)
      omega_f_np (n, ii_write, jj_write, :) = (/-3.80701604e-05, -9.42479237e-05, -6.99302182e-05,  4.66579731e-05, &
        4.45159618e-04,  8.71945522e-04,  2.20187684e-03,  1.21356396e-03,&
        3.22627323e-03,  2.34250221e-02,  4.35893722e-02,  4.04218324e-02,&
        1.70142725e-02, -4.56110388e-03, -1.21179726e-02, -1.76784862e-03,&
        3.73918228e-02,  6.48340657e-02,  6.30719066e-02,  6.82196394e-02,&
        8.88698921e-02,  1.08857207e-01,  1.12968870e-01,  1.04503594e-01,&
        9.35490355e-02,  8.29649270e-02,  6.75555691e-02,  4.81525064e-02,&
        3.17999721e-02,  1.94539521e-02,  9.97573510e-03,  2.93357461e-03,&
       -1.89023686e-03 /)
      pt_np      (n, ii_write, jj_write, :) = (/259.4465 , 239.507  , 231.76668, 227.12183, 223.35129, 219.27335,&
       214.6297 , 209.7517 , 205.4097 , 204.3045 , 210.90666, 219.23444,&
       223.9434 , 231.10272, 240.71408, 251.65454, 261.12726, 268.93414,&
       275.5061 , 280.0773 , 283.73364, 287.1842 , 289.92114, 291.61136,&
       294.117  , 288.8731 , 285.14566, 286.28885, 287.58633, 288.649  ,&
       289.4885 , 290.12668, 290.59073/)
      qq_np      (n, ii_write, jj_write, :) = (/1.7764813e-06, 1.8224730e-06, 1.8787580e-06, 1.8764346e-06,&
       1.8566434e-06, 1.8211915e-06, 1.6808081e-06, 1.6499712e-06,&
       2.4174512e-06, 2.1596550e-06, 2.4104070e-06, 3.0660274e-06,&
       1.1503719e-05, 4.9068629e-05, 3.3379259e-04, 6.3536590e-04,&
       1.1464410e-03, 1.6599769e-03, 1.7649520e-03, 3.3047206e-03,&
       4.9139624e-03, 5.6036133e-03, 5.3928299e-03, 5.3772102e-03,&
       4.5086094e-03, 3.5713578e-03, 7.6173735e-03, 8.7193158e-03,&
       8.9478139e-03, 9.0318723e-03, 9.0940222e-03, 9.1573792e-03,&
       9.2759673e-03/)

      n=3
      ph_np      (n, ii_write, jj_write, :) = (/1.00000000e+02, 4.00000000e+02, 8.18602100e+02, 1.37888600e+03,&
       2.09179500e+03, 2.98364100e+03, 4.12179000e+03, 5.57922200e+03,&
       7.42954864e+03, 9.74070741e+03, 1.25750915e+04, 1.59913577e+04,&
       2.00496361e+04, 2.48020526e+04, 3.02938144e+04, 3.65306013e+04,&
       4.34076761e+04, 5.06290908e+04, 5.77685357e+04, 6.44633141e+04,&
       7.05272028e+04, 7.59071104e+04, 8.06115814e+04, 8.46752207e+04,&
       8.81502576e+04, 9.10986727e+04, 9.35839055e+04, 9.56592278e+04,&
       9.73771992e+04, 9.87797733e+04, 9.99002937e+04, 1.00772376e+05,&
       1.01406544e+05, 1.01824297e+05/)
      omega_f_np (n, ii_write, jj_write, :) = (/-1.5646929e-04, -3.5355880e-04, -4.8761873e-04, -6.1965146e-04, &
       -5.1226310e-04, -5.6678458e-04,  8.9274318e-04,  1.0616900e-04, &
       -7.3587806e-05,  1.4965368e-02,  3.0534089e-02,  3.4318790e-02,&
        2.7774157e-02,  2.0657543e-02, -4.2138384e-03, -3.0219957e-02,&
        1.1667107e-02,  6.2012438e-02,  7.1553968e-02,  7.6540418e-02,&
        8.5921444e-02,  9.5617831e-02,  9.9485628e-02,  9.2792720e-02,&
        7.9984531e-02,  6.9873638e-02,  5.9503235e-02,  4.2936314e-02,&
        2.5582071e-02,  1.1460497e-02,  7.9563487e-04, -6.9983769e-03,&
       -1.2269462e-02 /)
      pt_np      (n, ii_write, jj_write, :) = (/259.41943, 239.47134, 231.75943, 227.06264, 223.29308, 219.17064,&
       214.7469 , 209.50075, 205.0984 , 203.46696, 210.22252, 219.17328,&
       224.35645, 231.65575, 241.06029, 251.39082, 260.98102, 269.74472,&
       276.25543, 280.3309 , 284.13727, 287.54175, 290.84753, 292.05405,&
       294.74423, 292.49033, 286.38773, 286.28198, 287.5152 , 288.57535,&
       289.41718, 290.0592 , 290.52917/)
      qq_np      (n, ii_write, jj_write, :) = (/1.7770907e-06, 1.8218230e-06, 1.8782342e-06, 1.8740881e-06,&
       1.8556129e-06, 1.8172200e-06, 1.6830505e-06, 1.6588933e-06,&
       2.5348015e-06, 2.3578232e-06, 2.6123425e-06, 2.8268344e-06,&
       1.1949348e-05, 4.0499141e-05, 3.1161102e-04, 5.2557967e-04,&
       8.1448845e-04, 1.6482240e-03, 2.3045661e-03, 4.1260906e-03,&
       5.3721131e-03, 5.6491895e-03, 5.7562813e-03, 5.3291898e-03,&
       5.1380335e-03, 3.3276873e-03, 6.5492103e-03, 8.8409791e-03,&
       9.2665097e-03, 9.3802800e-03, 9.4541032e-03, 9.5247030e-03,&
       9.6436469e-03/)

      n=4
      ph_np      (n, ii_write, jj_write, :) = (/1.00000000e+02, 4.00000000e+02, 8.18602100e+02, 1.37888600e+03,&
       2.09179500e+03, 2.98364100e+03, 4.12179000e+03, 5.57922200e+03,&
       7.42992209e+03, 9.74214078e+03, 1.25782211e+04, 1.59968007e+04,&
       2.00580135e+04, 2.48139971e+04, 3.03099739e+04, 3.65516141e+04,&
       4.34340876e+04, 5.06612030e+04, 5.78063042e+04, 6.45063989e+04,&
       7.05751104e+04, 7.59593014e+04, 8.06675208e+04, 8.47343994e+04,&
       8.82122078e+04, 9.11629749e+04, 9.36501907e+04, 9.57271699e+04,&
       9.74465120e+04, 9.88502043e+04, 9.99716113e+04, 1.00844377e+05,&
       1.01479031e+05, 1.01897094e+05/)
      omega_f_np (n, ii_write, jj_write, :) = (/4.7953494e-05, -3.6808473e-05,  3.5740519e-05,  1.9297475e-04,&
        3.1848802e-04,  1.2558178e-03,  2.9671229e-03,  2.5163901e-03,&
        3.1901787e-03,  2.0888757e-02,  4.2862691e-02,  3.5957269e-02,&
        1.0077652e-02, -1.5303524e-02, -1.6749227e-02,  1.7087562e-02,&
        2.4977865e-02,  1.8438639e-02,  3.3242777e-02,  4.7898896e-02,&
        6.0433261e-02,  6.2855504e-02,  6.3961752e-02,  7.3570587e-02,&
        8.2023971e-02,  7.4456051e-02,  5.4907698e-02,  3.6152974e-02,&
        2.3210175e-02,  1.4129854e-02,  7.4563599e-03,  2.6335313e-03,&
       -6.1447278e-04 /)
      pt_np      (n, ii_write, jj_write, :) = (/259.7153 , 239.78131, 231.94554, 227.3303 , 223.52255, 219.37563,&
       214.8545 , 210.29414, 206.33498, 206.26268, 213.11755, 219.66954,&
       222.79913, 229.94989, 240.24171, 251.54802, 261.00104, 267.46176,&
       274.22842, 279.84683, 283.5188 , 287.21613, 289.5535 , 292.4041 ,&
       290.52347, 284.68314, 284.52908, 286.03766, 287.347  , 288.40305,&
       289.23462, 289.86377, 290.31464/)
      qq_np      (n, ii_write, jj_write, :) = (/1.7752716e-06, 1.8119434e-06, 1.8585856e-06, 1.8827270e-06,&
       1.8695056e-06, 1.8330701e-06, 1.6761986e-06, 1.6601078e-06,&
       2.3278369e-06, 2.3174000e-06, 2.2752315e-06, 5.1278198e-06,&
       1.3375928e-05, 6.2444713e-05, 3.4220610e-04, 7.7448162e-04,&
       1.1322818e-03, 1.8825288e-03, 1.8853350e-03, 3.2220115e-03,&
       4.7324067e-03, 4.1785790e-03, 3.4078218e-03, 3.2965872e-03,&
       2.6125771e-03, 5.4507689e-03, 8.2147652e-03, 8.6833518e-03,&
       8.7782005e-03, 8.8373795e-03, 8.8883061e-03, 8.9439852e-03,&
       9.0650013e-03/)

      n=5
      ph_np      (n, ii_write, jj_write, :) = (/1.00000000e+02, 4.00000000e+02, 8.18602100e+02, 1.37888600e+03,&
       2.09179500e+03, 2.98364100e+03, 4.12179000e+03, 5.57922200e+03,&
       7.42968367e+03, 9.74122565e+03, 1.25762230e+04, 1.59933256e+04,&
       2.00526650e+04, 2.48063713e+04, 3.02996570e+04, 3.65381987e+04,&
       4.34172254e+04, 5.06407013e+04, 5.77821912e+04, 6.44788918e+04,&
       7.05445242e+04, 7.59259805e+04, 8.06318068e+04, 8.46966172e+04,&
       8.81726562e+04, 9.11219217e+04, 9.36078715e+04, 9.56837928e+04,&
       9.74022598e+04, 9.88052382e+04, 9.99260791e+04, 1.00798408e+05,&
       1.01432752e+05, 1.01850617e+05/)
      omega_f_np (n, ii_write, jj_write, :) = (/-6.3231586e-05, -2.7610478e-04, -3.5218231e-04, -9.3256996e-05,&
        2.1956344e-04,  2.2489096e-04,  2.1414873e-03,  2.6651162e-03,&
        8.0242194e-04,  1.9178173e-02,  4.3063503e-02,  3.7578147e-02,&
        1.6121805e-02, -6.2171775e-03, -1.5616204e-02,  2.6258023e-03,&
        1.3603573e-02,  1.1346980e-02,  2.4330884e-02,  5.1727772e-02,&
        7.7207938e-02,  8.3295152e-02,  7.9164259e-02,  7.8530915e-02,&
        8.3224148e-02,  8.2255416e-02,  6.8430178e-02,  4.7953449e-02,&
        3.0316550e-02,  1.7493878e-02,  8.0931075e-03,  1.2869849e-03,&
       -3.3162215e-03 /)
      pt_np      (n, ii_write, jj_write, :) = (/259.69287, 239.75174, 231.9556 , 227.32692, 223.46388, 219.2427 ,&
       214.97527, 210.04463, 205.91733, 205.3323 , 211.97504, 219.77632,&
       223.2343 , 230.41588, 240.5069 , 251.62495, 261.03403, 268.54645,&
       275.3341 , 280.2932 , 284.02798, 287.4824 , 290.00717, 292.0785 ,&
       293.2936 , 287.34592, 284.62305, 285.9292 , 287.21512, 288.27133,&
       289.1053 , 289.7381 , 290.1943/)
      qq_np      (n, ii_write, jj_write, :) = (/1.7756363e-06, 1.8115660e-06, 1.8606580e-06, 1.8800641e-06,&
       1.8692566e-06, 1.8331184e-06, 1.6782950e-06, 1.6621650e-06,&
       2.3519117e-06, 2.0480265e-06, 2.2392971e-06, 3.3703277e-06,&
       1.3221874e-05, 6.0117793e-05, 3.7239285e-04, 7.5924041e-04,&
       1.2746624e-03, 1.5166403e-03, 1.4585337e-03, 2.9409439e-03,&
       4.5657489e-03, 5.3046015e-03, 4.3280907e-03, 4.1191629e-03,&
       3.0206868e-03, 3.6845051e-03, 7.8127496e-03, 8.6197015e-03,&
       8.8479845e-03, 8.9301132e-03, 8.9909090e-03, 9.0530468e-03,&
       9.1692153e-03/)

      n=6
      ph_np      (n, ii_write, jj_write, :) = (/1.00000000e+02, 4.00000000e+02, 8.18602100e+02, 1.37888600e+03,&
       2.09179500e+03, 2.98364100e+03, 4.12179000e+03, 5.57922200e+03,&
       7.42939719e+03, 9.74012609e+03, 1.25738223e+04, 1.59891502e+04,&
       2.00462385e+04, 2.47972084e+04, 3.02872608e+04, 3.65220794e+04,&
       4.33969648e+04, 5.06160675e+04, 5.77532184e+04, 6.44458407e+04,&
       7.05077734e+04, 7.58859439e+04, 8.05888948e+04, 8.46512202e+04,&
       8.81251332e+04, 9.10725943e+04, 9.35570230e+04, 9.56316733e+04,&
       9.73490888e+04, 9.87512094e+04, 9.98713701e+04, 1.00743175e+05,&
       1.01377146e+05, 1.01794773e+05/)
      omega_f_np (n, ii_write, jj_write, :) = (/-0.00010322, -0.00038814, -0.00068608, -0.00061517, -0.00035002,&
       -0.00092965,  0.00078254,  0.00159943, -0.00247145,  0.01351988,&
        0.0351072 ,  0.03307324,  0.02088267,  0.0014126 , -0.02795155,&
       -0.03900209, -0.00514794,  0.03117283,  0.04580639,  0.05880513,&
        0.06726731,  0.06786502,  0.06424921,  0.05976126,  0.06019506,&
        0.06403068,  0.06187926,  0.04875997,  0.0307681 ,  0.01476835,&
        0.00226887, -0.00692206, -0.01313566/)
      pt_np      (n, ii_write, jj_write, :) = (/259.65872, 239.70558, 231.92563, 227.28397, 223.42842, 219.1079 ,&
       215.04976, 209.83934, 205.49283, 204.43167, 211.19556, 219.72845,&
       223.71588, 230.95016, 240.76646, 251.40977, 260.91693, 269.5335 ,&
       276.30588, 280.77014, 284.59546, 287.7794 , 290.74222, 292.12775,&
       294.7637 , 290.96698, 285.71646, 285.69888, 286.94574, 288.00418,&
       288.84497, 289.48608, 289.95566/)
      qq_np      (n, ii_write, jj_write, :) = (/1.7765080e-06, 1.8114630e-06, 1.8620234e-06, 1.8773757e-06,&
       1.8686152e-06, 1.8323010e-06, 1.6801076e-06, 1.6690658e-06,&
       2.4474382e-06, 2.2226275e-06, 2.3388525e-06, 2.7896328e-06,&
       1.3447439e-05, 5.0336628e-05, 3.5922619e-04, 5.8292417e-04,&
       9.6952851e-04, 1.3408856e-03, 1.4487513e-03, 2.9519931e-03,&
       4.1361265e-03, 5.3722425e-03, 5.1302067e-03, 4.2821206e-03,&
       3.8701706e-03, 3.0000824e-03, 6.3300505e-03, 8.5702883e-03,&
       8.9733787e-03, 9.1105038e-03, 9.1870734e-03, 9.2576258e-03,&
       9.3720397e-03/)

    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !elseif (input_profile .eq. 1) then
    !  input_profile_longname = ""
    !  ndim=6
    !  n=1
    !  ph_np      (n, ii_write, jj_write, :) = (//)
    !  omega_h_np (n, ii_write, jj_write, :) = (//)
    !  pt_np      (n, ii_write, jj_write, :) = (//)
    !  qq_np      (n, ii_write, jj_write, :) = (//)
    !  
    !  vi command: s/,$/,\&/gc
    !  add 0. to omega_h

    else
     call error_mesg('compute_vert_adv_tend_offline',  &
                     'input_profile is not supported', FATAL)

    end if ! end if of input_profile


!-------------------------
! loop for all profiles
!-------------------------

do n=1,ndim

  !--- assign values
  ph(:,:,:) = ph_np(n,:,:,:) 
  omega_f(:,:,:) = omega_f_np(n,:,:,:) 
  pt(:,:,:) = pt_np(n,:,:,:) 
  qq(:,:,:) = qq_np(n,:,:,:) 

    ! --- compute dp, pi_fac, theta
    do k = 2,kdim+1
      dp(:,:,k-1) = ph(:,:,k) - ph(:,:,k-1)
    enddo

    !--- get pressure at full levels
    do k = 2,kdim+1
      pf(:,:,k-1)      = 0.5 * (ph(:,:,k) + ph(:,:,k-1))
    enddo

    !--- get omega at half levels
    omega_h(:,:,kdim+1) = 0.
    omega_h(:,:,1)      = omega_f(:,:,1)
    do k = 2,kdim
      omega_h(:,:,k) = 0.5 * (omega_f(:,:,k-1) + omega_f(:,:,k))
    enddo

    !--- compute temperature tendency due to adiabatic heating
    dT_adi=rdgas*pt(:,:,:)*omega_f(:,:,:)/cp_air/pf(:,:,:)

    !--- compute temperature vertical advection
    select case (temp_vert_advec_scheme)
     case(1)
        call vert_advection(dts,omega_h,dp,pt,dT_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
     case(2)
        call vert_advection(dts,omega_h,dp,pt,dT_vadv,scheme=FOURTH_CENTERED,form=ADVECTIVE_FORM)
     case(3)
        call vert_advection(dts,omega_h,dp,pt,dT_vadv,scheme=FINITE_VOLUME_LINEAR,form=ADVECTIVE_FORM)
     case(4)
        call vert_advection(dts,omega_h,dp,pt,dT_vadv,scheme=FINITE_VOLUME_PARABOLIC,form=ADVECTIVE_FORM)
     case(5)
        call vert_advection(dts,omega_h,dp,pt,dT_vadv,scheme=SECOND_CENTERED_WTS,form=ADVECTIVE_FORM)
     case(6)
        call vert_advection(dts,omega_h,dp,pt,dT_vadv,scheme=FOURTH_CENTERED_WTS,form=ADVECTIVE_FORM)
     end select

    !---
    tdt_vadv = dT_vadv + dT_adi

    !--- compute tracer vertical advection
     SELECT CASE (tracer_vert_advec_scheme)
            CASE(1)
                 vadvec_scheme = SECOND_CENTERED
            CASE(2)
                 vadvec_scheme = FOURTH_CENTERED
            CASE(3)
                 vadvec_scheme = FINITE_VOLUME_LINEAR
            CASE(4)
                 vadvec_scheme = FINITE_VOLUME_PARABOLIC
            CASE(5)
                 vadvec_scheme = SECOND_CENTERED_WTS
            CASE(6)
                 vadvec_scheme = FOURTH_CENTERED_WTS
     END SELECT
     call vert_advection(dts,omega_h,dp,qq,dqv_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)

     !--- write out
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     if (printout_vert_adv_tend.eq.1) then   ! print out input and output fields
        write(6,*)    '#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        write(6,*)    ''
        write(6,*)    '#=========================='
        write(6,*)    '#=========================='
        write(6,*)    ''
        write(6,*)    '#   input '
        write(6,*)    '#   input_profile: ', input_profile, trim(input_profile_longname)
        write(6,*)    '#   n-th of profile : ', n
        write(6,*)    ''
        write(6,*)    '#=========================='
        write(6,*)    '#=========================='
        write(6,*)    ''
        write(6,*)    '# pressure at half level (Pa)'
        write(6,3001) '  ph  = [',ph(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '# pressure at full level (Pa)'
        write(6,3001) '  pf  = [',pf(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '# pressure thickness (Pa)'
        write(6,3001) '  dp  = [',dp(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '# vertical velocity at half levels (Pa/s)'
        write(6,3002) '  omega_h  = ['    ,omega_h(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '# vertical velocity at full levels (Pa/s)'
        write(6,3002) '  omega_f  = ['    ,omega_f(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '# temperature at full levels (K)'
        write(6,3001) '  pt  = [',pt(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '# specific humidity at full levels (kg/kg)'
        write(6,3002) '  qq  = ['    ,qq(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '#=========================='
        write(6,*)    '#=========================='
        write(6,*)    ''
        write(6,*)    '#   output '
        write(6,*)    '#   input_profile: ', input_profile, trim(input_profile_longname)
        write(6,*)    '#   n-th of profile : ', n
        write(6,*)    ''
        write(6,*)    '#=========================='
        write(6,*)    '#=========================='
        write(6,*)    ''
        write(6,*)    '# temperature vertical advection tendency (K/s)'
        write(6,3002) '  dT_vadv  = ['    ,dT_vadv(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '# temperature tendency due to adiabatic heating (K/s)'
        write(6,3002) '  dT_adi   = ['    ,dT_adi(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '# temperature tendency due to adiabatic heating and advection (K/s)'
        write(6,3002) '  tdt_vadv   = ['    ,tdt_vadv(ii_write,jj_write,:)
        write(6,*)    ''
        write(6,*)    '# specific humidity vertical advection tendency (K/s)'
        write(6,3002) '  dqv_vadv  = ['    ,dqv_vadv(ii_write,jj_write,:)
        write(6,*)    ''
     end if ! (printout_vert_adv_tend.eq.1) then

     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     if (printout_vert_adv_tend.eq.2) then   ! print out
        write(6,*)    ''
        write(6,*)    '#================================='
        write(6,*)    '#   input_profile: ', input_profile, trim(input_profile_longname)
        write(6,*)    '#   n-th of profile : ', n
        write(6,*)    '#================================='
        write(6,*)    ''

       do i=1,6
         select case (i)
          case(1)
             call vert_advection(dts,omega_h,dp,pt,dT_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
          case(2)
             call vert_advection(dts,omega_h,dp,pt,dT_vadv,scheme=FOURTH_CENTERED,form=ADVECTIVE_FORM)
          case(3)
             call vert_advection(dts,omega_h,dp,pt,dT_vadv,scheme=FINITE_VOLUME_LINEAR,form=ADVECTIVE_FORM)
          case(4)
             call vert_advection(dts,omega_h,dp,pt,dT_vadv,scheme=FINITE_VOLUME_PARABOLIC,form=ADVECTIVE_FORM)
          case(5)
             call vert_advection(dts,omega_h,dp,pt,dT_vadv,scheme=SECOND_CENTERED_WTS,form=ADVECTIVE_FORM)
          case(6)
             call vert_advection(dts,omega_h,dp,pt,dT_vadv,scheme=FOURTH_CENTERED_WTS,form=ADVECTIVE_FORM)
          end select
     
         !--- compute tracer vertical advection
          SELECT CASE (i)
                 CASE(1)
                      vadvec_scheme = SECOND_CENTERED
                 CASE(2)
                      vadvec_scheme = FOURTH_CENTERED
                 CASE(3)
                      vadvec_scheme = FINITE_VOLUME_LINEAR
                 CASE(4)
                      vadvec_scheme = FINITE_VOLUME_PARABOLIC
                 CASE(5)
                      vadvec_scheme = SECOND_CENTERED_WTS
                 CASE(6)
                      vadvec_scheme = FOURTH_CENTERED_WTS
          END SELECT
          call vert_advection(dts,omega_h,dp,qq,dqv_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)

          write(6,*)    '### vadvec_scheme: ', i
          write(6,3002) '  dT_vadv   = ['    ,dT_vadv(ii_write,jj_write,:)
          write(6,3002) '  dqv_vadv   = ['    ,dqv_vadv(ii_write,jj_write,:)

        end do  ! end loop of i
     end if ! (printout_vert_adv_tend.eq.2) then

end do   ! end loop of n

     !--- stop the model
     call error_mesg('compute_vert_adv_tend_offline',  &
                     'stop', FATAL)

3000 format (A35,2X,F10.3,',')
3001 format (A35,2X,34(F10.3,2X,','))
3002 format (A35,2X,34(E12.4,2X,','))
3003 format (A35,2X,E12.4,',')
3004 format (A35,2X,33(F10.3,2X,','),A5)

  end subroutine compute_vert_adv_tend_offline
!---> yhc, 2022-07-15


end module physics_driver_mod
