module scm_cgils_mod

   use            fms_mod, only:  check_nml_error,                      &
                                  mpp_pe, mpp_root_pe,                  &
                                  write_version_number,                 &
                                  open_file, close_file, file_exist,    &
                                  read_data, write_data,                &
                                  mpp_error,                            &
                                  error_mesg, FATAL, NOTE,              &
                                  field_exist, field_size
#ifdef INTERNAL_FILE_NML
   USE              mpp_mod, ONLY: input_nml_file
#else
   USE              fms_mod, ONLY: open_namelist_file
#endif

   use   diag_manager_mod, only:  register_diag_field, send_data

   use            mpp_mod, only:  mpp_pe, mpp_root_pe, stdlog
   use         mpp_io_mod, only:  mpp_open, mpp_close, mpp_get_times, mpp_get_info,  &
                                  MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE  ! ZNT 03/29/2021
   use   diag_manager_mod, only:  register_diag_field, send_data

   use   time_manager_mod

   use sat_vapor_pres_mod, only:  lookup_es,  compute_qs, sat_vapor_pres_init   ! ZNT 02/20/2020
   use      constants_mod, only:  rdgas, rvgas, cp_air, hlv, kappa, grav, pi, SECONDS_PER_DAY

   use      scm_utils_mod, only:  us_std_atm, locate, thetal_to_temp

   use vert_advection_mod, only:  vert_advection, SECOND_CENTERED, &
                                  FOURTH_CENTERED, FINITE_VOLUME_LINEAR, &
                                  FINITE_VOLUME_PARABOLIC, &
                                  SECOND_CENTERED_WTS, FOURTH_CENTERED_WTS, &
                                  ADVECTIVE_FORM
   use   field_manager_mod, only: MODEL_ATMOS
   use  tracer_manager_mod, only: get_tracer_index

! use       constants_mod, only: kappa
use             fv_pack, only: nlon, mlat, nlev, beglat, endlat, beglon, &
                               endlon, rlonb, rlatb,  cold_start, ncnst, &
                               pnats, consv_te, ptop, fv_init, fv_domain, &
                               fv_end, change_time, p_var, restart_format, area, &
                               ak, bk, rlon, rlat, ng_d, f_d, nt_prog, get_eta_level


   implicit none
   private

   public cgils_data_read, cgils_forc_init, cgils_forc_end, update_cgils_forc, &
          get_cgils_flx, &
          cgils_forc_diagnostic_init

   character(len=8) :: mod_name = 'scm_cgils'
   character(len=7) :: mod_name_diag = 'forcing'

!--------------------- version number ----------------------------------
!
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
integer, dimension(1) :: restart_versions = (/ 1 /)

!--- puclib variables
logical, private               :: initialized = .false.
character(len=200)             :: cgils_nc 
integer, dimension(4)          :: siz   
integer                        :: nlev_cgils
integer                        :: ntime_cgils
integer, parameter  :: nlat_cgils = 1
integer, parameter  :: nlon_cgils = 1
real, allocatable, dimension(:,:,:) :: buffer_cgils_3D
real, allocatable, dimension(:,:,:) :: buffer_cgils_3D_time
real, allocatable, dimension(:,:)   :: buffer_cgils_2D
real, allocatable, dimension(:)     :: buffer_cgils_1D_lev
real, allocatable, dimension(:)     :: &  ! dimension (nlev)
                                       T_cgils, U_cgils, V_cgils, q_cgils, pfull_cgils, &   ! initial conditions
                                       omega_cgils, divT_cgils, divq_cgils                  ! forcing
real, allocatable, dimension(:)     :: Ps_cgils   ! dimension (ntime)

real,    private               :: missing_value = -999.
integer :: nsphum, nql, nqi, nqa, nqn, nqni
integer :: vadvec_scheme

!--- namelist parameters
integer, public                :: tracer_vert_advec_scheme = 3
integer, public                :: temp_vert_advec_scheme = 3
integer, public                :: momentum_vert_advec_scheme = 3
character(len=20) , public     :: do_surface_fluxes = 'none'
character(len=20) , public     :: do_nudge_terms = 'u_v_t_q'
real, public                   :: tao_nudging = 3600.         ! nudging scale. units: second
real, public                   :: plev_nudging = 600.e+2      ! nudging above this pressure level (Pa) 
character(len=20) , public     :: cgils_case    = "none"
character(len=200), public     :: cgils_ctl_s6_nc  = '/ncrc/home2/Yi-hsuan.Chen/work/research/edmf_AM4/data/CGILS/ctl_s6.nc'
character(len=200), public     :: cgils_ctl_s11_nc = '/ncrc/home2/Yi-hsuan.Chen/work/research/edmf_AM4/data/CGILS/ctl_s11.nc'
character(len=200), public     :: cgils_ctl_s12_nc = '/ncrc/home2/Yi-hsuan.Chen/work/research/edmf_AM4/data/CGILS/ctl_s12.nc'
character(len=200), public     :: cgils_p2k_s6_nc  = '/ncrc/home2/Yi-hsuan.Chen/work/research/edmf_AM4/data/CGILS/p2k_s6.nc'
character(len=200), public     :: cgils_p2k_s11_nc = '/ncrc/home2/Yi-hsuan.Chen/work/research/edmf_AM4/data/CGILS/p2k_s11.nc'
character(len=200), public     :: cgils_p2k_s12_nc = '/ncrc/home2/Yi-hsuan.Chen/work/research/edmf_AM4/data/CGILS/p2k_s12.nc'
character(len=200), public     :: dephy_nc = '/ncrc/home2/Yi-hsuan.Chen/work/research/edmf_AM4/code/xanadu_SCM.original/BOMEX_REF_SCM_driver.nc'

integer :: do_stop = -1
integer :: do_debug_printout = -1  ! =-1, do not call
                                   ! = 1, call it in cgils_forc_init

namelist /scm_cgils_nml/ &
                        do_stop, do_debug_printout, &
                        cgils_case, do_surface_fluxes, do_nudge_terms, tao_nudging, plev_nudging, &
                        tracer_vert_advec_scheme, &
                        temp_vert_advec_scheme,   &
                        momentum_vert_advec_scheme

!--- diagnostics
integer ::  id_tdt_vadv, id_tdt_adi, id_tdt_adiPvadv, id_tdt_lf, id_tdt_nudge,                      &
            id_qvdt_vadv, id_qvdt_lf, id_qvdt_nudge,                                 &
            id_qldt_vadv, id_qidt_vadv, id_qadt_vadv,     id_qndt_vadv,     & 
            id_udt_vadv, id_udt_nudge, id_udt_lf,                                   &
            id_vdt_vadv, id_vdt_nudge, id_vdt_lf,                                      &
            id_pf_forc, id_ph_forc 

!#######################################################################

contains

!#######################################################################
! Subroutine to read case specific namelist
!   Copy from subroutine bomex_data_read, scm_bomex.F90

subroutine cgils_data_read()

implicit none

integer                                   :: i
integer                                   :: unit,ierr,io, logunit
    
integer :: year, month, day

   if (initialized) return
   initialized = .true.

!-------- read namelist --------
#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=scm_cgils_nml, iostat=io)
   ierr = check_nml_error(io, 'scm_cgils_nml')
#else
   if (file_exist('input.nml')) then
      unit = open_namelist_file()
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=scm_cgils_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'scm_cgils_nml')
      enddo
10    call close_file (unit)
   endif
#endif

!--------- write version number and namelist --------

   call write_version_number ( version, tagname )
   if(mpp_pe() == mpp_root_pe() ) then
     logunit =stdlog()
     write(logunit,nml=scm_cgils_nml)
   endif

!if (do_debug_printout) then
!    call error_mesg( ' scm_cgils',     &
!                     ' data_read successfully..',&
!                     FATAL ) 
!endif

end subroutine cgils_data_read


!#######################################################################
! Subroutine to initialize case forcings

subroutine cgils_forc_init(time_interp,As,Bs)
#include "fv_arrays.h"
#include "fv_point.inc"

!      VARIABLES
!      ------
!      INPUT:
!      ------
!      time_interp     time
!      As, Bs          A's and B's of half levels in hybrid coordinate
!                      ph(k) = A(k) + B(k) * psurf

  type(time_type)                          :: time_interp
  real,  intent (in), dimension(:)         :: As,Bs

!------------------
! local variables
!------------------

  real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: &   ! dimension(lat, lon, lev+1)
    phalf   ! pressure at half levels (Pa)

  real, dimension(size(pt,1),size(pt,2),size(pt,3)) :: &     ! dimension(lat, lon, lev)
    pfull  ! pressure at full levels (Pa)
 
  real dum1

  integer itime
  integer i,j,k,ix,jx,kx,tx

!---------------------------------

  nsphum = get_tracer_index(MODEL_ATMOS, 'sphum')
  nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
  nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
  nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
  nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )   
  nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' ) 

!------------------------------- 
! Determine CGILS input file
!------------------------------- 
  if (trim(cgils_case) .eq. "ctl_s6") then
    cgils_nc = cgils_ctl_s6_nc
  elseif (trim(cgils_case) .eq. "ctl_s11") then
    cgils_nc = cgils_ctl_s11_nc
  elseif (trim(cgils_case) .eq. "ctl_s12") then
    cgils_nc = cgils_ctl_s12_nc
  elseif (trim(cgils_case) .eq. "p2k_s6") then
    cgils_nc = cgils_p2k_s6_nc
  elseif (trim(cgils_case) .eq. "p2k_s11") then
    cgils_nc = cgils_p2k_s11_nc
  elseif (trim(cgils_case) .eq. "p2k_s12") then
    cgils_nc = cgils_p2k_s12_nc
  else
    call error_mesg( ' scm_cgils',     &
                     ' nml cgils_case is not supported.',&
                     FATAL ) 
  endif

!---------------------
! allocate variables
!---------------------

  !--- set buffer variables
  call field_size (cgils_nc, 'lev' , siz) ; nlev_cgils  = siz(1)
  call field_size (cgils_nc, 'time', siz) ; ntime_cgils = siz(1)
  
  if (allocated(buffer_cgils_3D)) deallocate (buffer_cgils_3D)
    allocate(buffer_cgils_3D(nlat_cgils, nlon_cgils, nlev_cgils)); buffer_cgils_3D = missing_value

  if (allocated(buffer_cgils_3D_time)) deallocate (buffer_cgils_3D_time)
    allocate(buffer_cgils_3D_time(nlat_cgils, nlon_cgils, ntime_cgils)); buffer_cgils_3D_time = missing_value

  if (allocated(buffer_cgils_2D)) deallocate (buffer_cgils_2D)
    allocate(buffer_cgils_2D(nlat_cgils, nlon_cgils)); buffer_cgils_2D = missing_value

  if (allocated(buffer_cgils_1D_lev)) deallocate (buffer_cgils_1D_lev)
    allocate(buffer_cgils_1D_lev(nlev_cgils)); buffer_cgils_1D_lev = missing_value

  !--- set cgils variables
  !if (allocated(T_cgils)) deallocate(T_cgils)
  !  allocate(T_cgils(nlev_cgils)); T_cgils = missing_value

!--------------------------------
!  read data from cgils file
!    subroutine read_data, $scmsrc/FMS/fms/fms_io.F90
!
!  4D variable in CGILS file is var_4D(time, lev, lat, lon). buffer_cgils_3D (lat, lon, lev). Time dimension is dropped
!  3D variable in CGILS file is var_3D(time, lat, lon).      buffer_cgils_2D (lat, lon)     . Time dimension is dropped
!--------------------------------

  !--- read CGILD 4D variables
  if (allocated(T_cgils)) deallocate(T_cgils) ; allocate(T_cgils(nlev_cgils)); T_cgils = missing_value
  call read_data(cgils_nc, 'T',   buffer_cgils_3D(:,:,:), no_domain=.true.)
       T_cgils(:) = buffer_cgils_3D(1,1,:)

  if (allocated(U_cgils)) deallocate(U_cgils) ; allocate(U_cgils(nlev_cgils)); U_cgils = missing_value
  call read_data(cgils_nc, 'U',   buffer_cgils_3D(:,:,:), no_domain=.true.)
       U_cgils(:) = buffer_cgils_3D(1,1,:)

  if (allocated(V_cgils)) deallocate(V_cgils) ; allocate(V_cgils(nlev_cgils)); V_cgils = missing_value
  call read_data(cgils_nc, 'V',   buffer_cgils_3D(:,:,:), no_domain=.true.)
       V_cgils(:) = buffer_cgils_3D(1,1,:)

  if (allocated(q_cgils)) deallocate(q_cgils) ; allocate(q_cgils(nlev_cgils)); q_cgils = missing_value
  call read_data(cgils_nc, 'q',   buffer_cgils_3D(:,:,:), no_domain=.true.)
       q_cgils(:) = buffer_cgils_3D(1,1,:)

  !--- read CGILD 3D variables & save to SCM
  !      [ps, u_srf, v_srf] are used by SCM; they are defined in /src/atmos_fv_dynamics/model/fv_point.inc
  !
  !      Using BOMEX_REF_SCM_driver.nc, sfc_sens_flx(time, lat, lon), read_data OK.
  !        call read_data(dephy_nc, 'sfc_sens_flx', buffer_cgils_2D(:,:), no_domain=.true., timelevel=itime)
  !        write(6,*) 'dephy_nc, sfc_sens_flx',buffer_cgils_2D
  !
  !      But this does not work for cgils nc. I have no idea why this happend.
  !        "FATAL: fms_io(read_data_3d_new), field Ps in file /ncrc/home2/Yi-hsuan.Chen/work/research/edmf_AM4/data/CGILS/ctl_s11.nc: field size mismatch 1". 
  !         gxsize, gysize, size(data, 3), siz_in(1), siz_in(2), siz_in(3) (1, 1, 1, 1, 1, 4)

  if (allocated(Ps_cgils)) deallocate(Ps_cgils) ; allocate(Ps_cgils(ntime_cgils)); Ps_cgils = missing_value
  call read_data(cgils_nc, 'Ps',   buffer_cgils_3D_time(:,:,:))
       Ps_cgils(:) = buffer_cgils_3D_time(1,1,:)
       ps(1,1) = Ps_cgils(1)

  call read_data(cgils_nc, 'u_srf',   buffer_cgils_3D_time(:,:,:))
       u_srf(:,:) = buffer_cgils_3D_time(:,:,1)

  call read_data(cgils_nc, 'v_srf',   buffer_cgils_3D_time(:,:,:))
       v_srf(:,:) = buffer_cgils_3D_time(:,:,1)

  !--- read CGILD 1D variables
  if (allocated(pfull_cgils)) deallocate(pfull_cgils) ; allocate(pfull_cgils(nlev_cgils)); pfull_cgils = missing_value
  call read_data(cgils_nc, 'lev',   buffer_cgils_1D_lev(:), no_domain=.true.)
       pfull_cgils(:) = buffer_cgils_1D_lev(:)

!-------------------------------------------------
! interpoate cgils data on SCM pressure levels
!   variables [ua, va, pt, q, ps] are used by SCM; they are defined in /src/atmos_fv_dynamics/model/fv_point.inc
!-------------------------------------------------

  !--- initialize
  ua=0.; va=0.; pt=0.; q=0.

  !--- get pfull and phalf
  call get_eta_level(nlev, ps(1,1), pfull(1,1,:), phalf(1,1,:))

  !--- interpolate cgils data on SCM pressure levels
  call interp_cgils_to_SCM(pfull(1,1,:), ua(1,1,:), pfull_cgils(:), U_cgils(:))
  call interp_cgils_to_SCM(pfull(1,1,:), va(1,1,:), pfull_cgils(:), V_cgils(:))
  call interp_cgils_to_SCM(pfull(1,1,:), pt(1,1,:), pfull_cgils(:), T_cgils(:))
  call interp_cgils_to_SCM(pfull(1,1,:), q(1,1,:,nsphum), pfull_cgils(:), q_cgils(:))

!-------------------------------------------------
! use subroutine p_var to set up varialbes in fv_point.inc, namely, delp, pk, pe, pkz, and peln.
!
! If these variables are not set properly, the SCM would fail in longwave_utilities.F90. The error message would be like
!
!   NOTE: lw_gases_stdtf_mod: Reading NetCDF formatted input data file: INPUT/cns_co2_600_HITRAN2012_10701200_495lyr.nc
!   NOTE: lw_gases_stdtf_mod: Reading NetCDF formatted input data file: INPUT/cns_co2_360_HITRAN2012_10701200_495lyr.nc
!   forrtl: severe (174): SIGSEGV, segmentation fault occurred
!   Image              PC                Routine            Line        Source
!   fms_SCM_am4_xanad  0000000004FA31F3  Unknown               Unknown  Unknown
!   fms_SCM_am4_xanad  0000000003096040  Unknown               Unknown  Unknown
!   fms_SCM_am4_xanad  00000000017613DC  longwave_utilitie         477  longwave_utilities.F90
!   fms_SCM_am4_xanad  00000000009A1FD8  sealw99_mod_mp_e1        4105  sealw99.F90
!   fms_SCM_am4_xanad  00000000009662CE  sealw99_mod_mp_se        1394  sealw99.F90
!   fms_SCM_am4_xanad  000000000139E58A  longwave_driver_m         481  longwave_driver.F90
!   fms_SCM_am4_xanad  00000000008CF03B  radiation_driver_        4247  radiation_driver.F90
!   fms_SCM_am4_xanad  00000000008AB286  radiation_driver_        2025  radiation_driver.F90
!   fms_SCM_am4_xanad  0000000000488954  atmos_model_mod_m         384  atmos_model.F90
!
!-------------------------------------------------

  ! --- Create delp (from hydro_eq in init_dry_atm.F90)
  do k=1,size(pt,3)
    do j=1,size(pt,2)
      do i=1,size(pt,1)
        delp(i,j,k) = As(k+1)-As(k) + ps(i,j)*(Bs(k+1)-Bs(k))
      enddo
    enddo
  enddo
  call p_var(nlon, mlat, nlev, beglat, endlat, ptop, delp, ps,     &
             pe, peln,  pk,  pkz,  kappa, q, ng_d, ncnst, .false. )

!------------
! debug
!------------

if (do_debug_printout.eq.1 .or. do_debug_printout.eq.99) then
  write(6,*) 'nlev_cgils, ntime_cgils',nlev_cgils, ntime_cgils
  write(6,*) 'siz',siz
  write(6,*) 'pfull_cgils',pfull_cgils
  write(6,*) 'pfull',pfull
  write(6,*) 'T_cgils',T_cgils
  write(6,*) 'pt',pt
  write(6,*) 'U_cgils',U_cgils
  write(6,*) 'V_cgils',V_cgils
  write(6,*) 'q_cgils',q_cgils
  write(6,*) 'Ps_cgils',Ps_cgils
  !write(6,*) 'phalf',phalf
  !write(6,*) '',

  if (do_stop.eq.1) then
    call error_mesg( ' scm_cgils. cgils_forc_init',     &
                     ' STOP.',&
                     FATAL ) 
  endif
endif

end subroutine cgils_forc_init

!#######################################################################
! Subroutine to end case forcings

subroutine cgils_forc_end ()
character*64                 :: fname_res='RESTART/cgils.res.nc'

  if (.not.initialized) return
  initialized = .false.


end subroutine cgils_forc_end


!#######################################################################
! Subroutine to update CGILS forcings

subroutine update_cgils_forc(time_interp,time_diag,dt_int)
#include "fv_arrays.h"
#include "fv_point.inc"

!      ------
!      INPUT:
!      ------
!      time_interp  time for interpolation
!      time_diag    time for diagnostics
!      dt_int       time step

type(time_type), intent(in)              :: time_interp,time_diag,dt_int

!------------------
! local variables
!------------------

  real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: &   ! dimension(lat, lon, lev+1)
    omga_half, &   ! omega at SCM half levels (Pa/s)
    phalf          ! pressure at half levels (Pa)

  real, dimension(size(pt,1),size(pt,2),size(pt,3)) :: &     ! dimension(lat, lon, lev)
    T_cgils_scm, U_cgils_scm, V_cgils_scm, q_cgils_scm,  &   ! T, q, U, V of CGILS on SCM vertical levels 
    dT_adi, dT_vadv, dT_adiPvadv, dqv_vadv, dql_vadv, dqi_vadv, dqa_vadv, dqn_vadv, du_vadv, dv_vadv, &  ! vertical tendencies 
    dT_nudge, dqv_nudge, du_nudge, dv_nudge,  &   ! nudging tendencies
    dT_lf, dqv_lf, du_lf, dv_lf,              &   ! large-scale horizontal tendencies
    pfull  ! pressure at full levels (Pa)

  real    :: dts
  integer :: dt_seconds,dt_days

  integer i,j,k,kdim

  logical                                          :: used

  real dum1

!-------------------------------------------

  kdim = size(pt,3)

!---------------------
! allocate variables
!---------------------

  !--- set buffer variables
  call field_size (cgils_nc, 'lev' , siz) ; nlev_cgils  = siz(1)
  call field_size (cgils_nc, 'time', siz) ; ntime_cgils = siz(1)
  
  if (allocated(buffer_cgils_3D)) deallocate (buffer_cgils_3D)
    allocate(buffer_cgils_3D(nlat_cgils, nlon_cgils, nlev_cgils)); buffer_cgils_3D = missing_value

!--------------------------------
!  read data from cgils file
!    subroutine read_data, $scmsrc/FMS/fms/fms_io.F90
!
!  4D variable in CGILS file is var_4D(time, lev, lat, lon). buffer_cgils_3D (lat, lon, lev). Time dimension is dropped
!--------------------------------

  if (allocated(Ps_cgils)) deallocate(Ps_cgils) ; allocate(Ps_cgils(ntime_cgils)); Ps_cgils = missing_value
  call read_data(cgils_nc, 'Ps',   buffer_cgils_3D_time(:,:,:))
       Ps_cgils(:) = buffer_cgils_3D_time(1,1,:)

  !--- read CGILD 1D variables
  if (allocated(pfull_cgils)) deallocate(pfull_cgils) ; allocate(pfull_cgils(nlev_cgils)); pfull_cgils = missing_value
  call read_data(cgils_nc, 'lev',   buffer_cgils_1D_lev(:), no_domain=.true.)
       pfull_cgils(:) = buffer_cgils_1D_lev(:)

  if (allocated(T_cgils)) deallocate(T_cgils) ; allocate(T_cgils(nlev_cgils)); T_cgils = missing_value
  call read_data(cgils_nc, 'T',   buffer_cgils_3D(:,:,:), no_domain=.true.)
       T_cgils(:) = buffer_cgils_3D(1,1,:)

  if (allocated(U_cgils)) deallocate(U_cgils) ; allocate(U_cgils(nlev_cgils)); U_cgils = missing_value
  call read_data(cgils_nc, 'U',   buffer_cgils_3D(:,:,:), no_domain=.true.)
       U_cgils(:) = buffer_cgils_3D(1,1,:)

  if (allocated(V_cgils)) deallocate(V_cgils) ; allocate(V_cgils(nlev_cgils)); V_cgils = missing_value
  call read_data(cgils_nc, 'V',   buffer_cgils_3D(:,:,:), no_domain=.true.)
       V_cgils(:) = buffer_cgils_3D(1,1,:)

  if (allocated(q_cgils)) deallocate(q_cgils) ; allocate(q_cgils(nlev_cgils)); q_cgils = missing_value
  call read_data(cgils_nc, 'q',   buffer_cgils_3D(:,:,:), no_domain=.true.)
       q_cgils(:) = buffer_cgils_3D(1,1,:)

  if (allocated(omega_cgils)) deallocate(omega_cgils) ; allocate(omega_cgils(nlev_cgils)); omega_cgils = missing_value
  call read_data(cgils_nc, 'omega',   buffer_cgils_3D(:,:,:), no_domain=.true.)
       omega_cgils(:) = buffer_cgils_3D(1,1,:)

  if (allocated(divT_cgils)) deallocate(divT_cgils) ; allocate(divT_cgils(nlev_cgils)); divT_cgils = missing_value
  call read_data(cgils_nc, 'divT',   buffer_cgils_3D(:,:,:), no_domain=.true.)
       divT_cgils(:) = buffer_cgils_3D(1,1,:)

  if (allocated(divq_cgils)) deallocate(divq_cgils) ; allocate(divq_cgils(nlev_cgils)); divq_cgils = missing_value
  call read_data(cgils_nc, 'divq',   buffer_cgils_3D(:,:,:), no_domain=.true.)
       divq_cgils(:) = buffer_cgils_3D(1,1,:)

!-------------------------------------------------
! interpoate cgils data on SCM pressure levels
!   variables [omga, delp] are used by SCM; they are defined in /src/atmos_fv_dynamics/model/fv_point.inc
!-------------------------------------------------

  !--- initialize
  omga = 0.0; omga_half=0.; dT_lf = 0.0; dqv_lf= 0.0; du_lf=0.; dv_lf=0.

  !--- get pfull and phalf, and delp
  call get_eta_level(nlev, Ps_cgils(1), pfull(1,1,:), phalf(1,1,:))

  do i=1,size(pt,1)
  do j=1,size(pt,2)
    do k=1,kdim
      delp(i,j,k) = phalf(i,j,k+1) - phalf(i,j,k)   ! pressure layer thickness (Pa)
    enddo
  enddo
  enddo

  !--- interpolate cgils data on SCM pressure levels
  call interp_cgils_to_SCM(pfull(1,1,:), T_cgils_scm (1,1,:), pfull_cgils(:), T_cgils(:))
  call interp_cgils_to_SCM(pfull(1,1,:), U_cgils_scm (1,1,:), pfull_cgils(:), U_cgils(:))
  call interp_cgils_to_SCM(pfull(1,1,:), V_cgils_scm (1,1,:), pfull_cgils(:), V_cgils(:))
  call interp_cgils_to_SCM(pfull(1,1,:), q_cgils_scm (1,1,:), pfull_cgils(:), q_cgils(:))
  call interp_cgils_to_SCM(pfull(1,1,:), omga        (1,1,:), pfull_cgils(:), omega_cgils(:))
  call interp_cgils_to_SCM(pfull(1,1,:), dT_lf       (1,1,:), pfull_cgils(:), divT_cgils(:))
  call interp_cgils_to_SCM(pfull(1,1,:), dqv_lf      (1,1,:), pfull_cgils(:), divq_cgils(:))

  call interp_cgils_to_SCM(phalf(1,1,:), omga_half(1,1,:), pfull_cgils(:), omega_cgils(:))
  omga_half(:,:,1)=0. ; omga_half(:,:,kdim+1)=0.  ! omega_half is zero at the top and the bottom levels

  ps(1,1) = Ps_cgils(1)

!-------------------------------------------------
! compute large-scale subsidence tendencies. For temperature, this consists of adiabatic warming and temperature vertical advection
!   module vert_advection_mod, src/atmos_shared/vert_advection/vert_advection.F90
!-------------------------------------------------

  !---initialization
  dT_vadv=0.0; dT_adi=0.0
  dqv_vadv=0.0; dql_vadv=0.0; dqi_vadv=0.0; dqa_vadv=0.0; dqn_vadv=0.0
  du_vadv=0.0; dv_vadv=0.0

  call get_time(dt_int,dt_seconds,dt_days)
  dts = real(dt_seconds + 86400*dt_days)

  !--- large-scale subsidence tendencies for temperature  
  dT_adi=rdgas*pt(:,:,:)*omga(:,:,:)/cp_air/phalf(:,:,:)  ! adiabatic term

  select case (temp_vert_advec_scheme)                 ! vertical advection term 
    case(1)
       call vert_advection(dts,omga_half,delp,pt,dT_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
    case(2)
       call vert_advection(dts,omga_half,delp,pt,dT_vadv,scheme=FOURTH_CENTERED,form=ADVECTIVE_FORM)
    case(3)
       call vert_advection(dts,omga_half,delp,pt,dT_vadv,scheme=FINITE_VOLUME_LINEAR,form=ADVECTIVE_FORM)
    case(4)
       call vert_advection(dts,omga_half,delp,pt,dT_vadv,scheme=FINITE_VOLUME_PARABOLIC,form=ADVECTIVE_FORM)
    case(5)
       call vert_advection(dts,omga_half,delp,pt,dT_vadv,scheme=SECOND_CENTERED_WTS,form=ADVECTIVE_FORM)
    case(6)
       call vert_advection(dts,omga_half,delp,pt,dT_vadv,scheme=FOURTH_CENTERED_WTS,form=ADVECTIVE_FORM)
  end select

  dT_adiPvadv = dT_vadv + dT_adi   ! sum of adiabatic term and vertical advection term

  !--- large-scale subsidence tendencies for tracers, namely, qv, ql, qi, qa
  select case (tracer_vert_advec_scheme)
    case(1)
         vadvec_scheme = SECOND_CENTERED
    case(2)
         vadvec_scheme = FOURTH_CENTERED
    case(3)
         vadvec_scheme = FINITE_VOLUME_LINEAR
    case(4)
         vadvec_scheme = FINITE_VOLUME_PARABOLIC
    case(5)
         vadvec_scheme = SECOND_CENTERED_WTS
    case(6)
         vadvec_scheme = FOURTH_CENTERED_WTS
  end select

  call vert_advection(dts,omga_half,delp,q(:,:,:,nsphum),dqv_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)
  call vert_advection(dts,omga_half,delp,q(:,:,:,nql),dql_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)
  call vert_advection(dts,omga_half,delp,q(:,:,:,nqi),dqi_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)
  call vert_advection(dts,omga_half,delp,q(:,:,:,nqa),dqa_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)

  if( nqn > 0 ) &
  call vert_advection(dts,omga_half,delp,q(:,:,:,nqn),dqn_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)

  !--- large-scale subsidence tendencies for momentum
  select case (momentum_vert_advec_scheme)
  case(1)
     call vert_advection(dts,omga_half,delp,ua,du_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
     call vert_advection(dts,omga_half,delp,va,dv_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
  end select

!-------------------------------------------------
! compute nudging tendencies
!   [ua, va, pt, q] are updated by SCM. Nudge these variables to CGILS profiles
!-------------------------------------------------

  du_nudge = 0.; dv_nudge = 0.; dT_nudge = 0.; dqv_nudge = 0.

  if (do_nudge_terms .eq. "none") then
    du_nudge  = 0.
    dv_nudge  = 0.
    dT_nudge  = 0.
    dqv_nudge = 0.

  elseif (do_nudge_terms .eq. "u_v_t_q") then
    do k=1, kdim
      if (pfull(1,1,k) < plev_nudging) then
        du_nudge (1, 1, k) = ( U_cgils_scm(1, 1, k) - ua(1, 1, k) ) / tao_nudging
        dv_nudge (1, 1, k) = ( V_cgils_scm(1, 1, k) - va(1, 1, k) ) / tao_nudging
        dT_nudge (1, 1, k) = ( T_cgils_scm(1, 1, k) - pt(1, 1, k) ) / tao_nudging
        dqv_nudge(1, 1, k) = ( T_cgils_scm(1, 1, k) - q(1, 1, k, nsphum) ) / tao_nudging
      endif
    enddo   

  elseif (do_nudge_terms .eq. "u_v") then
    do k=1, kdim
      if (pfull(1,1,k) < plev_nudging) then
        du_nudge(1, 1, k) = ( U_cgils_scm(1, 1, k) - ua(1, 1, k) ) / tao_nudging
        dv_nudge(1, 1, k) = ( V_cgils_scm(1, 1, k) - va(1, 1, k) ) / tao_nudging
      endif
    enddo   

  elseif (do_nudge_terms .eq. "t_q") then
    do k=1, kdim
      if (pfull(1,1,k) < plev_nudging) then
        dT_nudge (1, 1, k) = ( T_cgils_scm(1, 1, k) - pt(1, 1, k) ) / tao_nudging
        dqv_nudge(1, 1, k) = ( T_cgils_scm(1, 1, k) - q(1, 1, k, nsphum) ) / tao_nudging
      endif
    enddo   

  else
    call error_mesg( ' scm_cgils',     &
                     ' The value of do_nudge_terms is not supported. STOP',&
                     FATAL ) 
  endif

!----------------------------
! sum all tendencies
!----------------------------

  !--- wind
  u_dt = du_vadv + du_lf + du_nudge
  v_dt = dv_vadv + dv_lf + dv_nudge

  !--- temperature
  t_dt = dT_vadv + dT_adi + dT_lf

  !--- tracers
  q_dt(:,:,:,nsphum) = dqv_vadv + dqv_lf
  q_dt(:,:,:,nql) = dql_vadv
  q_dt(:,:,:,nqi) = dqi_vadv
  q_dt(:,:,:,nqa) = dqa_vadv
  if(nqn > 0) q_dt(:,:,:,nqn) = dqn_vadv

!----------------------------
! write out history files
!----------------------------

if (id_pf_forc > 0)  used = send_data( id_pf_forc,  pfull  (:,:,:), time_diag, 1, 1)

if (id_ph_forc > 0)  used = send_data( id_ph_forc,  phalf  (:,:,:), time_diag, 1, 1)

!--- subsidence
if (id_tdt_adi > 0)  used = send_data( id_tdt_adi, dT_adi (:,:,:), time_diag, 1, 1)
if (id_tdt_vadv > 0) used = send_data( id_tdt_vadv, dT_vadv (:,:,:), time_diag, 1, 1)
if (id_tdt_adiPvadv > 0)  used = send_data( id_tdt_adiPvadv, dT_adiPvadv (:,:,:), time_diag, 1, 1)
if (id_udt_vadv > 0) used = send_data( id_udt_vadv, du_vadv (:,:,:), time_diag, 1, 1)
if (id_vdt_vadv > 0) used = send_data( id_vdt_vadv, dv_vadv (:,:,:), time_diag, 1, 1)
if (id_qvdt_vadv> 0) used = send_data( id_qvdt_vadv,dqv_vadv(:,:,:), time_diag, 1, 1)
if (id_qldt_vadv> 0) used = send_data( id_qldt_vadv,dql_vadv(:,:,:), time_diag, 1, 1)
if (id_qidt_vadv> 0) used = send_data( id_qidt_vadv,dqi_vadv(:,:,:), time_diag, 1, 1)
if (id_qadt_vadv> 0) used = send_data( id_qadt_vadv,dqa_vadv(:,:,:), time_diag, 1, 1)

!--- nudging
if (id_tdt_nudge > 0) used = send_data( id_tdt_nudge, dT_nudge (:,:,:), time_diag, 1, 1)
if (id_udt_nudge > 0) used = send_data( id_udt_nudge, du_nudge (:,:,:), time_diag, 1, 1)
if (id_vdt_nudge > 0) used = send_data( id_vdt_nudge, dv_nudge (:,:,:), time_diag, 1, 1)
if (id_qvdt_nudge> 0) used = send_data( id_qvdt_nudge,dqv_nudge(:,:,:), time_diag, 1, 1)

!--- horizontal advection
if (id_udt_lf > 0) used = send_data( id_udt_lf, du_lf (:,:,:), time_diag, 1, 1)
if (id_vdt_lf > 0) used = send_data( id_vdt_lf, dv_lf (:,:,:), time_diag, 1, 1)
if (id_tdt_lf > 0) used = send_data( id_tdt_lf, dT_lf (:,:,:), time_diag, 1, 1)
if (id_qvdt_lf > 0) used = send_data( id_qvdt_lf, dqv_lf (:,:,:), time_diag, 1, 1)

!-------------
! debugging
!-------------

if (do_debug_printout.eq.2  .or. do_debug_printout.eq.99) then
  write(6,*) 'nlev_cgils, ntime_cgils',nlev_cgils, ntime_cgils
  write(6,*) '2omega_cgils',omega_cgils
  write(6,*) '2omga',omga
  write(6,*) '2omga_half',omga_half
  write(6,*) '2divT_cgils', divT_cgils
  write(6,*) '2dT_lf',dT_lf 
  write(6,*) '2divq_cgils', divq_cgils
  write(6,*) '2dqv_lf', dqv_lf 
  write(6,*) '2ua',ua
  write(6,*) '2va',va
  write(6,*) '2pt',pt
  write(6,*) '2qv',q(1,1,:,nsphum)
  write(6,*) '2T_cgils_scm',T_cgils_scm 
  write(6,*) '2U_cgils_scm',U_cgils_scm 
  write(6,*) '2V_cgils_scm',V_cgils_scm 
  write(6,*) '2q_cgils_scm',q_cgils_scm 
  write(6,*) '2dT_adi',dT_adi
  write(6,*) '2phalf',phalf

  !write(6,*) '', 

  !write(6,*) '_cgils', _cgils
  !write(6,*) '', 

  if (do_stop.eq.2) then
    call error_mesg( ' scm_cgils. update_cgils_forc',     &
                     ' STOP.',&
                     FATAL ) 
  endif

endif

end subroutine update_cgils_forc

!########################################################################
! This subroutine returns surface fluxes

subroutine get_cgils_flx( rho, u_star, flux_t, flux_q )

  !--- input argument
  real, intent(in), dimension(:) :: rho

  !--- output argument
  real, intent(out), dimension(:) :: &
     u_star,   & ! friction velocity (m/s)
     flux_t,   & ! surface sensible flux (W/m2)
     flux_q      ! surface latent heat flux (W/m2)

  !--- local variables

!-------------------------------------

  !--- intialize
  u_star=0.; flux_t=0.; flux_q=0.

  !--- get surface fluxes

  !--- prescribed surface fluxes
  if (trim(do_surface_fluxes).eq.'prescribed') then
    u_star = 0.25  ! use the same value in RF01

    call read_data(cgils_nc, 'shflx',   buffer_cgils_3D_time(:,:,:))
         flux_t(:) = buffer_cgils_3D_time(1,1,1)

    call read_data(cgils_nc, 'lhflx',   buffer_cgils_3D_time(:,:,:))
         flux_q(:) = buffer_cgils_3D_time(1,1,1)
    
  else
    call error_mesg( ' scm_cgils. get_cgils_flx. The value of do_surface_fluxes is not supported',     &
                     ' STOP.',&
                     FATAL ) 
  endif

if (do_debug_printout.eq.3  .or. do_debug_printout.eq.99) then
  write(6,*) 'gg3, flux_t',flux_t
  write(6,*) 'gg3, flux_q',flux_q
  write(6,*) 'gg3, u_star',u_star

  if (do_stop.eq.3) then
    call error_mesg( ' scm_cgils. get_cgils_flx',     &
                     ' STOP.',&
                     FATAL ) 
  endif
endif

end subroutine get_cgils_flx


!#######################################################################
! Subroutine to initialize CGILS diagnostic fields

subroutine cgils_forc_diagnostic_init(axes, Time)

implicit none

integer, dimension(3) :: half = (/1,2,4/)
integer, dimension(:),intent(in) :: axes
type(time_type), intent(in)      :: Time

!---------------------------------

! --- initialize axes -------------------------------------------------!
id_pf_forc = register_diag_field (mod_name_diag, 'pf_forc', axes(1:3), Time, &
     'Pressure at full level', 'hPa',  missing_value = missing_value)

id_ph_forc = register_diag_field (mod_name_diag, 'ph_forc', axes(half), Time, &
     'Pressure at half level', 'hPa',  missing_value = missing_value)

id_tdt_vadv = register_diag_field (mod_name_diag, 'tdt_vadv', axes(1:3), Time, &
     'Temperature tendencies due to vertical advection', 'K/s', missing_value = missing_value)

id_tdt_adi = register_diag_field (mod_name_diag, 'tdt_adi', axes(1:3), Time, &
     'Temperature tendencies due to subsidence adiabatic warming', 'K/s', missing_value = missing_value)

id_tdt_adiPvadv = register_diag_field (mod_name_diag, 'tdt_adiPvadv', axes(1:3), Time, &
     'Temperature tendencies due to adiabatic warming and vertical advection', 'K/s', missing_value = missing_value)

id_udt_vadv = register_diag_field (mod_name_diag, 'udt_vadv', axes(1:3), Time, &
     'U tendencies due to vertical advection', 'm/s2', missing_value = missing_value)

id_vdt_vadv = register_diag_field (mod_name_diag, 'vdt_vadv', axes(1:3), Time, &
     'V tendencies due to vertical advection', 'm/s2', missing_value = missing_value)

id_qvdt_vadv = register_diag_field (mod_name_diag, 'qvdt_vadv', axes(1:3), Time, &
     'Vapor tendencies due to vertical advection', 'kg/kg/s', missing_value = missing_value)

id_qldt_vadv = register_diag_field (mod_name_diag, 'qldt_vadv', axes(1:3), Time, &
     'liquid tendencies due to vertical advection', 'kg/kg/s', missing_value = missing_value)

id_qidt_vadv = register_diag_field (mod_name_diag, 'qidt_vadv', axes(1:3), Time, &
     'ice tendencies due to vertical advection', 'kg/kg/s', missing_value = missing_value)

id_qadt_vadv = register_diag_field (mod_name_diag, 'qadt_vadv', axes(1:3), Time, &
     'cloud amount tendencies due to vertical advection', '1/s', missing_value = missing_value)

id_qndt_vadv = register_diag_field (mod_name_diag, 'qndt_vadv', axes(1:3), Time, &
     'liquid droplet number concentration tendencies due to vertical advection', '1/cm3/s', missing_value = missing_value)

id_udt_nudge = register_diag_field (mod_name_diag, 'udt_nudge', axes(1:3), Time, &
     'U tendencies due to nudging', 'm/s2', missing_value = missing_value)

id_vdt_nudge = register_diag_field (mod_name_diag, 'vdt_nudge', axes(1:3), Time, &
     'V tendencies due to nudging', 'm/s2', missing_value = missing_value)

id_tdt_nudge = register_diag_field (mod_name_diag, 'tdt_nudge', axes(1:3), Time, &
     'Temperature tendencies due to nudging', 'K/s', missing_value = missing_value)

id_qvdt_nudge = register_diag_field (mod_name_diag, 'qvdt_nudge', axes(1:3), Time, &
     'Vapor tendencies due to nudging', 'kg/kg/s', missing_value = missing_value)

id_udt_lf = register_diag_field (mod_name_diag, 'udt_lf', axes(1:3), Time, &
     'U tendencies due to large-scale horizontal forcing', 'm/s2', missing_value = missing_value)

id_vdt_lf = register_diag_field (mod_name_diag, 'vdt_lf', axes(1:3), Time, &
     'V tendencies due to large-scale horizontal forcing', 'm/s2', missing_value = missing_value)

id_tdt_lf = register_diag_field (mod_name_diag, 'tdt_lf', axes(1:3), Time, &
     'Temperature tendencies due to large-scale horizontal forcing', 'K/s', missing_value = missing_value)

id_qvdt_lf = register_diag_field (mod_name_diag, 'qvdt_lf', axes(1:3), Time, &
     'Vapor tendencies due to large-scale horizontal forcing', 'kg/kg/s', missing_value = missing_value)


end subroutine cgils_forc_diagnostic_init

!#######################################################################
! Subroutine to interpolat CGILS data on SCM levels
!   reference: subroutine dephy_ini, scm_bomex.F90

subroutine interp_cgils_to_SCM (plev_SCM, var_SCM, plev_cgils, var_cgils)

  !--- input argument
  real, intent(in) , dimension(:) :: plev_SCM               ! pressure levels in SCM (Pa)  , e.g. pfull_scm (nlev_scm)
  real, intent(in) , dimension(:) :: plev_cgils, var_cgils  ! variable and pressure levels in CGILS (Pa), 
                                                            !   e.g. var_cgils(nlev_cgils), pfull_cgils(nlev_cgils)
  !--- output argument
  real, intent(out), dimension(:) :: var_SCM 

  !--- local variables
  logical :: found
  real    :: al, ah
  integer :: k_scm, k_cgils, kmax_scm, kmax_cgils ! d for DEPHY; k_dp = k_d + 1
  integer :: k, k1, k1p
!----------------------------------------

  !--- initialize output
  var_SCM = 0.

  !--- find kmax
  kmax_scm   = size(plev_SCM  ,1)
  kmax_cgils = size(plev_cgils,1)

  do k=1,kmax_scm
    !--- below CGILS data range
    if (plev_SCM(k) >= plev_cgils(kmax_cgils)) then
      var_SCM(k) = var_cgils(kmax_cgils)

    !--- above CGILS data range
    elseif (plev_SCM(k) <= plev_cgils(1)) then
      var_SCM(k) = var_cgils(1)

    !--- do interpolation within CGILS profile
    elseif (plev_SCM(k) > plev_cgils(1) .and. plev_SCM(k) < plev_cgils(kmax_cgils)) then
      found = .false.

      k1 = 1; k1p = k1 + 1
      found = .false.
      do while ( (k1p .le. kmax_cgils ) .and. (.not. found) )
        if ( plev_SCM(k) >= plev_cgils(k1) .and. plev_SCM(k) < plev_cgils(k1p) ) then
           ah = (plev_cgils(k1p) - plev_SCM(k)) / (plev_cgils(k1p) - plev_cgils(k1))
           al = 1.0 - ah
           var_SCM(k) = ah*var_cgils(k1) + al*var_cgils(k1p)
           found = .true.
        end if
        k1 = k1p ; k1p = k1 + 1 ! step up
      end do

    endif  ! end if of plev_SCM
  enddo    ! end loop of k


end subroutine interp_cgils_to_SCM



end module scm_cgils_mod

