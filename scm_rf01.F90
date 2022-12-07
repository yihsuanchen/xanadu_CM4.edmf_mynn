module scm_rf01_mod
 
   use            fms_mod, only:  check_nml_error,  &
                                  mpp_pe, mpp_root_pe,                  &
                                  write_version_number,                 &
                                  open_file, close_file, file_exist,    &
                                  read_data, write_data, mpp_error,     &
                                  error_mesg, FATAL, NOTE

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

   public rf01_data_read, rf01_forc_init, rf01_forc_end, update_rf01_forc, &
          rf01_forc_diagnostic_init, add_rf01_tdtlw, add_rf01_tdtsw,       &
          rf01_aer_init,  &  ! yhc add
          get_rf01_flx,  get_rf01_flx_online

   character(len=8) :: mod_name = 'scm_rf01'
   character(len=7) :: mod_name_diag = 'forcing'

   real, public, allocatable, dimension(:,:,:)  :: u_geos, v_geos, tdt_lw, tdt_sw

   ! arrays to hold initial sounding

   integer, parameter :: ksnd = 147
   real, dimension(ksnd) :: z_snd, p_snd, T_snd, qv_snd, ql_snd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!                       PARAMETERS OF THE MODULE
!
!
!       p00                      reference pressure (pascals)
!
!       configuration            case configuration (used only if multiple
!                                configurations share the same module)
!
!       tracer_vert_advec_scheme Which advection scheme should be used?
!
!       temp_vert_advec_scheme   Which advection scheme should be used?
!
!                                SECOND_CENTERED         1
!                                FOURTH_CENTERED         2
!                                FINITE_VOLUME_LINEAR    3 
!                                FINITE_VOLUME_PARABOLIC 4
!                                SECOND_CENTERED_WTS     5 
!                                FOURTH_CENTERED_WTS     6
!
!       vert_advec_cond          Should condensate be vertically 
!                                advected?
!
!       divf                     large scale divergence to compute omega
!
!       p_omega_zero             pressure above which omega is set to
!                                to zero
!
!       do_rad                   include radiation forcings?
!       do_geo                   include geostrophic forcings?
!       do_vadv                  include large scale subsidence?
!
!       psfc                     surface pressure (Pa)
!       zsfc                     topography height (m)
!       ustar_sfc                surface u_star (m/s)
!       flux_t_sfc               surface sensible flux (W/m2)
!       flux_q_sfc               surface latent heat flux (W/m2)

real,    private               :: p00 = 100000.

character(len=64)              :: configuration = 'base'

integer, public                :: tracer_vert_advec_scheme = 3
integer, public                :: temp_vert_advec_scheme = 3
integer, public                :: momentum_vert_advec_scheme = 3
logical, public                :: vert_advec_cond = .true.

real,    public                :: divf = 3.75e-6
real,    public                :: p_omega_zero = 70000.

logical, public                :: do_rad = .true.
logical, public                :: do_geo = .true.
logical, public                :: do_vadv = .true.

real,    private               :: psfc = 1017.8e2
real,    private               :: zsfc =    0.0
real,    private               :: ustar_sfc   =   0.25   
! ZNT 04/02/2020: For drag = (u_*)^2*u/sqrt(u^2+v^2) 
!                          = C_D*u*sqrt(u^2+v^2),
!     we have the effective u_* = sqrt[C_D*(u^2+v^2)] 
!                               = sqrt(0.0011*(7^2+5.5^2)) = 0.2953
!     but the actual difference is small, and the lowest atmos level
!     is unspecified.
real,    private               :: flux_t_sfc  =  15.0
real,    private               :: flux_q_sfc  = 115.0

real,    private               :: rho_i = 1.13  ! ZNT 05/19/2020

real,    private               :: missing_value = -999.

logical, private               :: initialized = .false.
logical                        :: do_netcdf_restart = .true.
logical                        :: do_stevens = .true.   ! ZNT 05/18/2020
logical                        :: do_aer_prof = .false. ! ZNT 03/29/2021
logical                        :: do_init_cloud_free = .false. ! yhc 05/10/2022
real                           :: zi_stevens = 840.     ! yhc 05/16/2022, the top of mixed layer top (m) 
                                                        ! Stevens et al. (2005) uses 840m
logical                        :: do_any_profiles     = .false.         ! yhc 2022-08-06, read specified profiles as initial conditions 
character*100                  :: option_any_profiles = "none"          ! yhc 2022-08-06, option for initial conditions
logical                        :: do_read_rf01_forc_any = .false.       ! yhc 2022-08-06, read specified forcing profiles
character*100                  :: option_read_rf01_forc_any = "none"    ! yhc 2022-08-06, option for specified forcing profiles
character*100                  :: option_lf_forcing = "none"            ! horizontal advective tendencies.
                                                                        !   availabe: "MERRA2_DYCOMS_3hr_10Jul1030Z"

integer :: nsphum, nql, nqi, nqa, nqn, nqni
integer :: vadvec_scheme

!--------------------- version number ----------------------------------
!
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
integer, dimension(1) :: restart_versions = (/ 1 /)
        
! namelist

namelist /scm_rf01_nml/ tracer_vert_advec_scheme,                   &
                        temp_vert_advec_scheme,                     &
                        momentum_vert_advec_scheme,                 &
                        vert_advec_cond,                            &
                        p_omega_zero, divf,                         &
                        do_rad, do_geo, do_vadv,                    &
                        ustar_sfc, flux_t_sfc, flux_q_sfc,          &
                        configuration, &
                        do_stevens, &                     ! ZNT 05/18/2020
                        do_init_cloud_free, &             ! yhc 05/10/2022
                        zi_stevens,         &             ! yhc 05/16/2022
                        do_any_profiles, option_any_profiles, do_read_rf01_forc_any, option_read_rf01_forc_any, & ! yhc 2022-08-06
                        option_lf_forcing,  & ! yhc 2022-09-21
                        do_aer_prof                       ! ZNT 03/29/2021


! diagnostics

integer ::  id_tdt_radf, id_tdt_vadv, id_tdt_lf,                &
            id_qvdt_vadv, id_qvdt_lf,                           &
            id_qldt_vadv, id_qidt_vadv, id_qadt_vadv,    id_qndt_vadv,         &
            id_udt_vadv, id_udt_geos, id_udt_lf,                &
            id_vdt_vadv, id_vdt_geos, id_vdt_lf,                &
            id_flx_radf, id_zi_forc,                            &
            id_pf_forc, id_ph_forc, id_zf_forc, id_zh_forc,     &
            id_tdt_dyn_forc, id_qvdt_dyn_forc,                    &  ! yhc, 2022-08-06
            id_u_geos, id_v_geos

! ---> h1g, 2010-09-27
integer ::    id_qvdt_forc_col
integer ::    id_qldt_vadv_col
! <--- h1g, 2010-09_27

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!#######################################################################
! Subroutine to read case specific namelist

subroutine rf01_data_read(kmax)

implicit none

integer,  intent (in)                     :: kmax
integer                                   :: unit,ierr,io, logunit
character*23                              :: tracer_ascheme,temp_ascheme
character*64                              :: fname

character*64                              :: fname_res='INPUT/rf01.res.nc'
! integer                                   :: vers  !restart file version
integer k

   if (initialized) return
   initialized = .true.
      
!-------- read namelist --------
#ifdef INTERNAL_FILE_NML
   READ (input_nml_file, nml=scm_rf01_nml, iostat=io)
   ierr = check_nml_error(io, 'scm_rf01_nml')
#else
   if (file_exist('input.nml')) then
      unit = open_namelist_file()
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=scm_rf01_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'scm_rf01_nml')
      enddo
10    call close_file (unit)
   endif
#endif

!--------- write version number and namelist --------

   call write_version_number ( version, tagname )
   if(mpp_pe() == mpp_root_pe() ) then
     logunit =stdlog()
     write(logunit,nml=scm_rf01_nml)
   endif

!--------- allocate memory ---------

   if (allocated(u_geos)) deallocate(u_geos);  allocate(u_geos(1,1,kmax))
   if (allocated(v_geos)) deallocate(v_geos);  allocate(v_geos(1,1,kmax))
   if (allocated(tdt_lw)) deallocate(tdt_lw);  allocate(tdt_lw(1,1,kmax))
   if (allocated(tdt_sw)) deallocate(tdt_sw);  allocate(tdt_sw(1,1,kmax))

! ---> h1g, 2010-10-05
   if (file_exist('INPUT/rf01.res.nc') ) then
      if(mpp_pe() == mpp_root_pe() ) call mpp_error ('scm_rf01_mod', &
         'Reading netCDF formatted restart file: INPUT/rf01.res.nc', NOTE)
      if( allocated(u_geos) ) call read_data (fname_res, 'u_geos',  u_geos)
      if( allocated(v_geos) ) call read_data (fname_res, 'v_geos',  v_geos)
   endif
! <--- h1g, 2010-10-05

       
!--------- read initial sounding from file ---------
  if (.not. do_stevens) then
   fname='INPUT/rf01_idealb.dat.new2'
   call mpp_open(unit,fname,action=MPP_RDONLY)
   do k=1,ksnd
      read(unit,'(f9.2,f8.2,f7.2,f7.2,f8.2)')  &
        z_snd(k),p_snd(k),T_snd(k),qv_snd(k),ql_snd(k)
   end do
   call close_file(unit)
  endif

end subroutine rf01_data_read

!#######################################################################
! Subroutine to initialize case forcings

subroutine rf01_forc_init(time_interp,As,Bs)
#include "fv_arrays.h"

!      VARIABLES
!      ------
!      INPUT:
!      ------
!      time_interp     time
!      As, Bs          A's and B's of half levels in hybrid coordinate
!                      ph(k) = A(k) + B(k) * psurf
!
!      -------
!      OUTPUT:
!      -------
!      elev            topography elevation
!      ps              surface pressure
!
!      -------------
!      INPUT/OUTPUT:
!      -------------
!      pt               temperature in degrees Kelvin 
!      u               zonal wind (m/s)
!      v               meridional wind (m/s)
!      qv              water vapor specific humidity (kg water vapor/     &
!                                                     kg air)
!      ql              liquid water condensate specific humidity 
!                      (kg condensate/ kg air)
!      qi              ice H2O condensate specific humidity 
!                      (kg condensate/kg dry air)
!      qa              cloud fraction (fraction)

type(time_type)                          :: time_interp
real,  intent (in), dimension(:)         :: As,Bs

!  Internal variables
!  ------------------

integer  :: kdim, k

integer, parameter :: itmax = 10

real, dimension(size(pt,3)+1) :: eta, peta
real, dimension(size(pt,3))   :: T_us_std, qv_us_std
real, dimension(size(pt,3))   :: u_rf01, v_rf01, T_rf01, qv_rf01, qa_rf01, ql_rf01 

real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph, zh, zhnew
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf, zf, zfnew

! ZNT 02/20/2020: Test
real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph0
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf0
real :: u_snd, v_snd, thetal_snd, qt_snd, T_snd, qv_snd, ql_snd  ! ZNT 05/18/2020

real maxerror

integer i
real,  dimension(size(pt,1),size(pt,2))    :: elev
integer :: j
#include "fv_point.inc"
       nsphum = get_tracer_index(MODEL_ATMOS, 'sphum')
       nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
       nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
       nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
       nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )   !h1g
       nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )   !h1g

!   Initialize surface pressure and topography

    ps = psfc
    elev = zsfc
!   ZNT 02/20/2020: initialize vapor pressure table
    call sat_vapor_pres_init
    
!   Setup temporary vertical grid: this grid is needed in order to use

    kdim = size(pt,3)
       ! --- Create delp (from hydro_eq in init_dry_atm.F90)
       do k=1,KDIM
         do j=1,size(ps,2)
           do i=1,size(ps,1)
             delp(i,j,k) = As(k+1)-As(k) + ps(i,j)*(Bs(k+1)-Bs(k))
           enddo
         enddo
       enddo
       call p_var(nlon, mlat, nlev, beglat, endlat, ptop, delp, ps,     &
              pe, peln,  pk,  pkz,  kappa, q, ng_d, ncnst, .false. )
       
       ! ZNT 02/20/2020: Note - '.false.' means no dry air mass adj,
       !     thus the calculations are consistent, given Bs(1)=0, and
       !     (ak,bk) from fv_pack = (As,Bs) as input, and
       !     i.e., ph (with get_eta level) = pe = p_half (from compute_p_z), 
       !           pf (with get_eta_level) = p_full (with compute_p_z).
       ! 
       ! p_var: pe(1) = ptop = ak(1); 
       !        pe(k) = pe(k-1) + delp(k-1)
       !              = pe(k-1) + As(k)-As(k-1) + ps*(Bs(k)-Bs(k-1))
       !              = ... = pe(1) + As(k)-As(1) + ps*(Bs(k)-Bs(1))
       !              = As(k) + Bs(k)*ps + (ak(1)-As(1)-Bs(1)*ps)
       ! compute_p_z: p_half = pe; p_full(k) = delp(k)/(peln(k+1)-peln(k))
       ! 
       ! get_eta_level: ph(1) = ak(1); ph(k) = ak(k) + bk(k)*ps
       !                pf(k) = (ph(k+1)-ph(k))/log(ph(k+1)/ph(k))
       !                (ak(1) > ptop_min = 1e-6 is always satisfied).
       
       do j=1,size(ps,2)
         do i=1,size(ps,1)
           ! call get_eta_level(nlev, ps(i,j) , pf(i,j,:), ph(i,j,:))
           call get_eta_level(nlev, ps(i,j) , pf0(i,j,:), ph0(i,j,:))  ! ZNT
         enddo
       enddo

!   Initialize pt, qv first with US standard atmosphere

    do k=1,kdim
      ! call us_std_atm( pf(1,1,k), T_us_std(k), qv_us_std(k) )
      call us_std_atm( pf0(1,1,k), T_us_std(k), qv_us_std(k) )
      pt(:,:,k)  = T_us_std(k)
      q(:,:,k,nsphum) = qv_us_std(k)
    end do
    q(:,:,:,nql) = 0.0
    ua  = 0.0
    va  = 0.0

!   Compute heights
    call compute_p_z (nlev, grav*elev, pe, peln, delp, pt, q(:,:,:,nsphum), pf, ph, zf, zh)
    ! ZNT 02/20/2020: DEBUG output
    write (*,*) 'PF0 - PF', pf0 - pf
    write (*,*) 'PH0 - PH', ph0 - ph
    write (*,*) 'ZF', zf
    write (*,*) 'ZH', zh
    ! ZNT 02/20/2020: DEBUG result: ph0 == ph, pf0 ~= pf (differ by 1e-9)

    !  call compute_height (Vgrid, grav*elev, T, qv, pf, ph, zf, zh )
    !        call error_mesg ('rf01 ',   &
    !        'need to compute height', FATAL )

!   Iteration to compute blended RF01 sounding

    i=0
    maxerror = 1.0
    do while ( maxerror > 0.001 .and. i < itmax )

     i = i + 1
     if (do_stevens) then   ! ZNT: 05/18/2020
         do k=1,kdim

           call rf01_snd_stevens( zf(1,1,k), u_snd, v_snd, thetal_snd, qt_snd )
           ! ZNT 05/15/2020: Note the following subroutine takes mixing ratio as input but outputs sphum.
           call thetal_to_temp( thetal_snd, qt_snd, pf(1,1,k), T_snd, qv_snd, ql_snd )
           pt(:,:,k)  = max ( T_snd, 200.0 )
           ua(:,:,k) = u_snd
           va(:,:,k) = v_snd
           q(:,:,k,nsphum) = qv_snd
           q(:,:,k,nql) = ql_snd

         enddo

     !--- read specified profiles, e.g. those from AM4
     elseif (do_any_profiles) then ! yhc, 2022-08-06
           call rf01_snd_any_profiles(u_rf01, v_rf01, T_rf01, qv_rf01, qa_rf01, ql_rf01)
           do k=1,kdim
             pt(:,:,k)  = T_rf01(k)
             ua(:,:,k) = u_rf01(k)
             va(:,:,k) = v_rf01(k)
             q(:,:,k,nsphum) = qv_rf01(k)
             q(:,:,k,nqa) = qa_rf01(k)
             q(:,:,k,nql) = ql_rf01(k)
           enddo

     else
         do k=1,kdim

          call rf01_snd( zf(1,1,k), u_rf01(k), v_rf01(k), T_rf01(k), qv_rf01(k), ql_rf01(k) )
          pt(:,:,k)  = T_rf01(k)
          ua(:,:,k) = u_rf01(k)
          va(:,:,k) = v_rf01(k)
          q(:,:,k,nsphum) = qv_rf01(k)
          q(:,:,k,nql) = ql_rf01(k)

         enddo
     endif
     call compute_p_z (nlev, grav*elev, pe, peln, delp, pt, q(:,:,:,nsphum), pf, ph, zfnew, zhnew)

     ! call compute_height (Vgrid, grav*elev, T, qv, pf, ph, zfnew, zhnew )
     !       call error_mesg ('rf01 ',   &
     !       'need to compute height', FATAL )
     maxerror = maxval( abs( zfnew - zf ) )
     zf = zfnew
     zh = zhnew

    enddo

    if ( i >= itmax ) then
      call error_mesg('rf01_forc_init',  &
                      'failed to converge while creating sounding', FATAL)
    endif

!   Initialize cloud amount

    !<--- yhc: set cloud-free initial condition
    if (do_init_cloud_free) then
      q(:,:,:,nql) = 0.
    end if
    !---> yhc

    if (.not.do_any_profiles) then ! yhc, 2022-11-29. Make sure the read-in cloud fraction is not overwritten here. 
      where ( q(:,:,:,nql) > 0.0 )
        q(:,:,:,nqa) = 1.0
      elsewhere
        q(:,:,:,nqa) = 0.0
      endwhere
    endif  

    if( nqn > 0) q(:,:,:,nqn) = 0.0

!   Initialize geostrophic winds

    u_geos = ua
    v_geos = va

    ! ZNT 02/23/2020 - Note: Also need to initialize u_srf and v_srf
    u_srf(:,:)=ua(:,:,KDIM)
    v_srf(:,:)=va(:,:,KDIM)

    ! ZNT 03/29/2021: Initialize aerosol profiles from AM4 AMIP July Climatology
    if (do_aer_prof) then
       call rf01_aer_init(kdim)
    end if

end subroutine rf01_forc_init

!#######################################################################

subroutine rf01_forc_end ()
character*64                 :: fname_res='RESTART/rf01.res.nc'

  if (.not.initialized) return
  initialized = .false.
  
! ---> h1g, 2010-10-05
     if( do_netcdf_restart ) then
           if (mpp_pe() == mpp_root_pe()) then
              call mpp_error ('scm_rf01_mod', 'Writing netCDF formatted restart file: RESTART/rf01.res.nc', NOTE)
           endif
           call write_data (fname_res, 'vers', restart_versions(size(restart_versions(:))) , no_domain=.true.)	   
           call write_data (fname_res,  'u_geos',  u_geos)
           call write_data (fname_res,  'v_geos',  v_geos)
      endif
! <--- h1g, 2010-10-05

  deallocate (  u_geos, v_geos, tdt_lw, tdt_sw )
 
end subroutine rf01_forc_end

!#######################################################################

subroutine rf01_forc_diagnostic_init(axes, Time)

implicit none

integer, dimension(3) :: half = (/1,2,4/)
integer, dimension(:),intent(in) :: axes
type(time_type), intent(in)      :: Time
   
! --- initialize axes -------------------------------------------------!
id_pf_forc = register_diag_field (mod_name_diag, 'pf_forc', axes(1:3), Time, &
     'Pressure at full level', 'hPa',  missing_value = missing_value)

id_ph_forc = register_diag_field (mod_name_diag, 'ph_forc', axes(half), Time, &
     'Pressure at half level', 'hPa',  missing_value = missing_value)

id_zf_forc = register_diag_field (mod_name_diag, 'zf_forc', axes(1:3), Time, &
     'Height at full level', 'm',  missing_value = missing_value)

id_zh_forc = register_diag_field (mod_name_diag, 'zh_forc', axes(half), Time, &
     'Height at half level', 'm',  missing_value = missing_value)

id_tdt_radf = register_diag_field (mod_name_diag, 'tdt_radf', axes(1:3), Time, &
     'Temperature tendencies due to radiative forcing', 'K/s',  missing_value = missing_value)

id_flx_radf = register_diag_field (mod_name_diag, 'flx_radf', axes(half), Time, &
     'Vertical radiative flux', 'W/m2', missing_value = missing_value)

id_tdt_vadv = register_diag_field (mod_name_diag, 'tdt_vadv', axes(1:3), Time, &
     'Temperature tendencies due to vertical advection', 'K/s', missing_value = missing_value)

id_udt_vadv = register_diag_field (mod_name_diag, 'udt_vadv', axes(1:3), Time, &
     'U tendencies due to vertical advection', 'm/s2', missing_value = missing_value)

id_vdt_vadv = register_diag_field (mod_name_diag, 'vdt_vadv', axes(1:3), Time, &
     'V tendencies due to vertical advection', 'm/s2', missing_value = missing_value)

id_udt_geos = register_diag_field (mod_name_diag, 'udt_geos', axes(1:3), Time, &
     'U tendencies due to geostrophic wind', 'm/s2', missing_value = missing_value)

id_vdt_geos = register_diag_field (mod_name_diag, 'vdt_geos', axes(1:3), Time, &
     'V tendencies due to geostrophic wind', 'm/s2', missing_value = missing_value)

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

id_udt_lf = register_diag_field (mod_name_diag, 'udt_lf', axes(1:3), Time, &
     'U tendencies due to large-scale forcing', 'm/s2', missing_value = missing_value)

id_vdt_lf = register_diag_field (mod_name_diag, 'vdt_lf', axes(1:3), Time, &
     'V tendencies due to large-scale forcing', 'm/s2', missing_value = missing_value)

id_tdt_lf = register_diag_field (mod_name_diag, 'tdt_lf', axes(1:3), Time, &
     'Temperature tendencies due to large-scale forcing', 'K/s', missing_value = missing_value)

id_qvdt_lf = register_diag_field (mod_name_diag, 'qvdt_lf', axes(1:3), Time, &
     'Vapor tendencies due to large-scale forcing', 'kg/kg/s', missing_value = missing_value)

id_zi_forc = register_diag_field (mod_name_diag, 'zi_forc', axes(1:2), Time, &
     'inversion height', 'm', missing_value = missing_value)

id_u_geos = register_diag_field (mod_name_diag, 'u_geos', axes(1:3), Time, &
     'U geostrophic wind', 'm/s',  missing_value = missing_value)

id_v_geos = register_diag_field (mod_name_diag, 'v_geos', axes(1:3), Time, &
     'V geostrophic wind', 'm/s',  missing_value = missing_value)

! ---> h1g, 2010-09-27
 id_qvdt_forc_col =  register_diag_field (mod_name_diag, 'qvdt_forc_col', axes(1:2), Time, &
     'column integrated vapor forcing', 'kg/m2/s',  missing_value = missing_value)
 id_qldt_vadv_col =  register_diag_field (mod_name_diag, 'qldt_vadv_col', axes(1:2), Time, &
     'column integrated cloud water vertical advection', 'kg/m2/s',  missing_value = missing_value)
! <--- h1g, 2010-09-27

!<--- yhc, 2022-08-06
id_tdt_dyn_forc = register_diag_field (mod_name_diag, 'tdt_dyn_forc', axes(1:3), Time, &
     'Temperature tendencies due to large-scale forcings', 'K/s', missing_value = missing_value)
id_qvdt_dyn_forc = register_diag_field (mod_name_diag, 'qvdt_dyn_forc', axes(1:3), Time, &
     'Vapor tendencies due to large-scale forcings', 'kg/kg/s', missing_value = missing_value)
!---> yhc, 2022-08-06

end subroutine rf01_forc_diagnostic_init

!#######################################################################
! Subroutine to apply RF01 forcings

subroutine update_rf01_forc(time_interp,time_diag,dt_int)
#include "fv_arrays.h"

!      ------
!      INPUT:
!      ------
!      time_interp  time for interpolation
!      time_diag    time for diagnostics
!      dt_int       time step
!      Hgrid        horizontal grid
!      Vgrid        vertical grid
!      elev         elevation

!      ------------
!      INPUT/OUTPUT
!      ------------
!
!      pt            temperature in degrees Kelvin 
!      u            zonal wind (m/s)
!      v            meridional wind (m/s)
!      qv           water vapor specific humidity
!                   (kg water vapor/kg air)
!      ql           liquid water condensate specific humidity 
!                   (kg condensate/ kg air)
!      qi           ice H2O condensate specific humidity
!                   (kg condensate/kg dry air)
!      qa           cloud fraction (fraction)
! 
!      ------
!      OUTPUT:
!      ------
!
!      ps               surface pressure
!      AT               temperature tendency due to large-scale forcing
!                       (K/sec)
!      AU               zonal wind tendency due to large-scale forcing
!                       (m/s*s)
!      AU               meridional wind tendency due to large-scale 
!                       forcing (m/s*s)
!      AQ               water vapor tendency due to large-scale forcing 
!                       (kg vapor/kg air/sec)
!      AL               liquid water condensate tendency due to 
!                       large-scale forcing (kg condensate/ kg air/sec)
!      AI               H2O ice condensate tendency due to large-scale
!                       forcing (kg condensate/ kg air/ sec)
!      AA               cloud fraction tendency due to large-scale 
!                       forcing (fraction / sec)
!      omega_f          omega interpolated to full model levels 
!                       (Pa/sec)
!

type(time_type), intent(in)              :: time_interp,time_diag,dt_int

integer                                          :: i,j,k,kdim
integer                                          :: dt_seconds,dt_days
logical                                          :: used
real                                             :: fcriolis, term3, dzi, Qz1, Qz2, dts
real, dimension(size(pt,3))                       :: klwp
real, dimension(size(pt,1),size(pt,2))             :: zi
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf,zf,th
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pi_fac, dp
real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph, omega_h, zh, frad
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: dT_rad, dT_adi, dT_vadv, dT_lf, &
                                                    tdt_any, qvdt_any, omega_any, tdt_dyn_forc, qvdt_dyn_forc, & ! yhc 2022-08-06
                                                    dqv_vadv, dqv_lf, dql_vadv, dqi_vadv, dqa_vadv, dqn_vadv
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: du_vadv, dv_vadv, du_geos, dv_geos, du_lf, dv_lf
real,  dimension(size(pt,1),size(pt,2))    :: elev  ! znt 20200226

! ---> h1g, 2010-09-27
real, dimension(size(pt,1),size(pt,2))  ::  qvdt_forcing_col
real, dimension(size(pt,1),size(pt,2))  ::  qldt_vadv_col
! <--- h1g, 2010-09-27
!  ------------------------------------------------------------------------------------------------------------
#include "fv_point.inc"

! --- update pf, ph, and zf, zh
elev=zsfc   ! znt 20200226
ps=psfc; ph(:,:,1)=0.;
       do j=1,size(ps,2)
         do i=1,size(ps,1)
           call get_eta_level(nlev, ps(i,j) , pf(i,j,:), ph(i,j,:))
         enddo
       enddo

       call compute_p_z (nlev, grav*elev, pe, peln, delp, pt, q(:,:,:,nsphum), pf, ph, zf, zh)

       ! call compute_height    (Vgrid, grav*elev, T, qv, pf, ph, zf, zh)
       ! call error_mesg ('rf02 ',   &
       !      'need to compute height', FATAL )


! --- compute large-scale subsidence
kdim = size(pt,3)
if (do_stevens) then     
! ZNT 05/18/2020: use linear profile for subsidence rather than a sharp cutoff
!                 following CLUBB-SCM.
   do k=1,kdim
      if ( zf(1,1,k) < 1600.0 ) then
         omga(:,:,k) = pf(1,1,k)/(rdgas*pt(1,1,k))*grav  &
                       * divf*zf(1,1,k)
      elseif ( zf(1,1,k) < 2400.0 ) then
         omga(:,:,k) = pf(1,1,k)/(rdgas*pt(1,1,k))*grav  &
                       * (divf*1600. - divf*2.*(zf(1,1,k)-1600.) )
      else
         omga(:,:,k) = 0.
      end if
   end do

   omega_h = 0.0
   do k=2,kdim
      if ( zh(1,1,k) < 1600.0 ) then
         omega_h(:,:,k) = ph(1,1,k)/(rdgas*0.5*(pt(1,1,k)+pt(1,1,k-1)))*grav &
                       * divf*zh(1,1,k)
      elseif ( zh(1,1,k) < 2400.0 ) then
         omega_h(:,:,k) = ph(1,1,k)/(rdgas*0.5*(pt(1,1,k)+pt(1,1,k-1)))*grav &
                       * (divf*1600. - divf*2.*(zh(1,1,k)-1600.))
      else
         omega_h(:,:,k) = 0.
      end if
   end do

else
   do k=1,kdim
      if (pf(1,1,k) > p_omega_zero) then
         omga(:,:,k) = pf(1,1,k)/(rdgas*pt(1,1,k))*grav * zf(1,1,k)*divf
      else
         omga(:,:,k) =0.
      end if
   end do

   omega_h=0.0
   do k=2,kdim
      if (ph(1,1,k) > p_omega_zero) then
         omega_h(:,:,k) = ph(1,1,k)/(rdgas*0.5*(pt(1,1,k)+pt(1,1,k-1))) * grav * zh(1,1,k)*divf
      else
         omega_h(:,:,k) =0.
      end if
   end do
endif

! --- compute dp, pi_fac, theta
do k = 2,kdim+1
   dp(:,:,k-1) = ph(:,:,k) - ph(:,:,k-1)
   pi_fac(:,:,k-1)= (pf(:,:,k-1)/p00)**(rdgas/cp_air)
enddo
th(:,:,:) = pt(:,:,:) / pi_fac(:,:,:)

! --- compute geostrophic tendencies
du_geos = 0.0
dv_geos = 0.0
if (do_geo) then
   fcriolis = f_d(1,1)
   do k=1, kdim
      du_geos(:,:,k) =   fcriolis * (va(1,1,k)-v_geos(1,1,k))
      dv_geos(:,:,k) = - fcriolis * (ua(1,1,k)-u_geos(1,1,k))
   end do
end if

! --- large-scale forcing
du_lf = 0.0
dv_lf = 0.0
dT_lf = 0.0
dqv_lf = 0.

if (trim(option_lf_forcing).eq."test") then
  ! the horizontal forcing is diagnosed from MERRA-2 3-hourly data *data*. The horizontal advective tendencies
  ! are diagnosed by NCL advect_variable function.

  dT_lf = -1./86400.
  dqv_lf = 1.e-3/86400.

  do k=1,kdim
    if ( zf(1,1,k) > 2400.0 ) then   ! set horizontal forcings to zero where omega is zero
         dT_lf (:,:,k) = 0.
         dqv_lf(:,:,k) = 0.
    end if
  end do

elseif (trim(option_lf_forcing).eq."none") then
  ! no horizontal forcing

  dT_lf = 0. 
  dqv_lf = 0.

elseif (trim(option_lf_forcing).eq."PBL_MERRA2_wCGILSsty_FreeTrop_0") then
  ! Similar to DYCOMS SCM setup,
  !   Z<1600m, CGILS-style horizontal advection using MERRA-2 data
  !   2400m>Z>1600m, Linear interpolation by using the values at 1600m and 2400m
  !   Z>2400m, No horizontal advection

  dT_lf(1,1,:) = (/ &
        0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, &
        0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, &
        0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, &
        0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, &
        0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, &
        0.000000e+00, -4.629755e-06, -1.835072e-05, -2.360223e-05, &
       -2.360223e-05, -2.360223e-05, -2.360223e-05, -2.360223e-05, &
       -2.360223e-05, -2.360223e-05, -2.360223e-05, -2.360223e-05, &
       -2.360223e-05/)

  dqv_lf(1,1,:) = (/ &
        0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, &
        0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, &
        0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, &
        0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, &
        0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, &
        0.000000e+00, -3.066373e-09, -1.215403e-08, -1.563220e-08, &
       -1.563220e-08, -1.563220e-08, -1.563220e-08, -1.563220e-08, &
       -1.563220e-08, -1.563220e-08, -1.563220e-08, -1.563220e-08, &
       -1.563220e-08 /)

elseif (trim(option_lf_forcing).eq."PBL_MERRA2_wCGILSsty_FreeTrop_MERRA2") then
  ! Simiar to CGILS setup, &
  ! P>900hPa, CGILS-style horizontal advection using MERRA-2 data
  ! 800hPa > P > 900hPa, Linear interpolation by using the values at 900 hPa and 800 hPa
  ! P<800hPa, MERRA-2 horizontal advection

  dT_lf(1,1,:) = (/ &
       -1.464067e-05,  4.225336e-05,  8.497543e-07, -2.187984e-05, &
       -1.101093e-05, -2.185082e-06, -7.134623e-06, -1.959490e-05, &
       -4.012173e-05, -8.633712e-05, -1.050757e-04, -1.093736e-04, &
       -3.294569e-05,  1.376715e-05,  5.901304e-06, -7.226663e-06, &
        3.141057e-06, -3.530690e-05, -4.581254e-05, -4.330122e-05, &
       -3.292953e-05, -2.428754e-05, -2.407418e-05, -2.389075e-05, &
       -2.373444e-05, -2.360223e-05, -2.360223e-05, -2.360223e-05, &
       -2.360223e-05, -2.360223e-05, -2.360223e-05, -2.360223e-05, &
       -2.360223e-05 /)

  dqv_lf(1,1,:) = (/ &
        4.233625e-13,  7.152405e-13, -4.695386e-13,  9.734660e-14, &
       -2.571354e-13, -1.187783e-13, -1.400123e-13, -5.445460e-14, &
        8.269690e-13, -1.635144e-12,  1.432642e-12,  4.472376e-11, &
        2.996606e-10, -8.855314e-11,  2.046432e-12,  2.496785e-09, &
       -3.834130e-09,  1.496060e-08,  1.833367e-08,  2.273357e-08, &
        4.736444e-09, -1.438761e-08, -1.477508e-08, -1.510822e-08, &
       -1.539208e-08, -1.563220e-08, -1.563220e-08, -1.563220e-08, &
       -1.563220e-08, -1.563220e-08, -1.563220e-08, -1.563220e-08, &
       -1.563220e-08/)

elseif (trim(option_lf_forcing).eq."CGILSstyle_hadv_RF01") then
  ! Simiar to CGILS setup, &
  ! P>900hPa, CGILS-style horizontal advection using MERRA-2 data
  ! 800hPa > P > 900hPa, Linear interpolation by using the values at 900 hPa and 800 hPa
  ! P<800hPa, 0
  !    2022-11-27 add
  dT_lf (1,1,:) = &
    (/0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, -6.11707006e-06, -1.50210283e-05, &
       -2.26061568e-05, -2.36022260e-05, -2.36022260e-05, -2.36022260e-05, &
       -2.36022260e-05, -2.36022260e-05, -2.36022260e-05, -2.36022260e-05, &
       -2.36022260e-05/)
  
  dqv_lf(1,1,:) = &
     (/0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, -4.05145059e-09, -9.94870965e-09, &
       -1.49724830e-08, -1.56321984e-08, -1.56321984e-08, -1.56321984e-08, &
       -1.56321984e-08, -1.56321984e-08, -1.56321984e-08, -1.56321984e-08, &
       -1.56321984e-08/)

else
  call error_mesg('update_rf01_forc',  &
                  'option_lf_forcing nml value is not supported', FATAL)

endif  ! end if of option_lf_forcing


!<-- yhc, 2022-08-07
if (do_read_rf01_forc_any) then
  call read_rf01_forc_any(tdt_any, qvdt_any, omega_any)
  omega_h(:,:,:) = omega_any(:,:,:)    ! omega_h is used in vert_advection
  omga   (:,:,:) = omega_any(:,:,:)    ! omga is in fv_point.inc
end if
!--> yhc, 2022-08-07

! --- longwave radiative forcing
!     F(z) = F0*exp(-Q(z,inf))
!          + F1*exp(-Q(0,z))
!          + alpha*rho_i*cp*Divf*H(z-zi)*(0.25(z-zi)4/3+zi(z-zi)1/3)
dT_rad=0.
if (do_rad) then
   do k=1, kdim
      klwp(k)=85.* q(1,1,k,nql) * dp(1,1,k)/grav
   end do
   frad(1,1,:)=0.;

      do i = kdim, 2, -1
         if ( ( q(1,1,i,nsphum)/( 1.-q(1,1,i,nsphum) )+q(1,1,i,nql)   )<0.008 ) then 

!h1g
           ! zi(1,1) = zf(1,1,i) ! original cjg
           zi(1,1) =  zf(1,1,i)*( (  q(1,1,i+1,nsphum)/( 1.- q(1,1,i+1,nsphum) )+q(1,1,i+1,nql))-0.008) &
                     +zf(1,1,i+1)*( 0.008-( q(1,1,i,nsphum)/( 1.- q(1,1,i,nsphum) ) +q(1,1,i,nql)) )

           zi(1,1) = zi(1,1) / ( ( q(1,1,i+1,nsphum)/( 1.- q(1,1,i+1,nsphum) )+q(1,1,i+1,nql)) &
                                  -( q(1,1,i,nsphum)/(1.- q(1,1,i,nsphum) )+q(1,1,i,nql)) )
          exit 
!h1g
         endif
      end do

   do k=1, kdim
      Qz1=0.; Qz2=0.;
      do i=1, k
         Qz1 = Qz1 + klwp(i)
      end do
      do i=k+1, kdim
         Qz2 = Qz2 + klwp(i)
      end do

      dzi=zf(1,1,k)-zi(1,1);
      !h1g
       !if ( (dzi > 0.) .and. (ph(1,1,k+1) > p_omega_zero)) then
       !   term3=1. * 1.12 * cp_air * divf * (0.25*dzi**(4./3.) + zi(1,1)*dzi**(1./3.))

      if (do_stevens) then
      ! ZNT 05/18/2020: use linear profile for subsidence rather than a sharp cutoff
      !                 following CLUBB-SCM.
      ! Note that although following exactly the Stevens paper, the temperature profile
      ! will still drift because the LW flux is constant in the free troposphere before
      ! rho-weighting.

         if (zi(1,1) > 1600.) then
         call error_mesg('rf01_snd',  &
                         'z inv higher than 1600m, dT/dt rad not valid', FATAL)
         else
            if (zf(1,1,k) > 2400.) then  ! No subsidence, constant term3
               term3 = cp_air * divf * (0.25*(1600.-zi(1,1))**(4./3.) + &
                                     zi(1,1)*(1600.-zi(1,1))**(1./3.))
               term3 = term3 + cp_air * divf * 6.0E-3 * &
                                   ( (3*1600.*2400. - 2400.**2.) - &
                                     (3*1600.*1600. - 1600.**2.))

            elseif (zf(1,1,k) > 1600.) then
               term3 = cp_air * divf * (0.25*(1600.-zi(1,1))**(4./3.) + &
                                     zi(1,1)*(1600.-zi(1,1))**(1./3.))
               term3 = term3 + cp_air * divf * 6.0E-3 * &
                                   ( (3*1600.*zf(1,1,k) - zf(1,1,k)**2.) - &
                                     (3*1600.*1600. - 1600.**2.)) 
            elseif (zf(1,1,k) > zi(1,1)) then
               term3= cp_air * divf * (0.25*(zf(1,1,k)-zi(1,1))**(4./3.) + &
                                    zi(1,1)*(zf(1,1,k)-zi(1,1))**(1./3.))
            else
               term3= 0.
            endif
            ! ZNT 05/19/2020: note the constant density in the Stevens (2005)
            !                 this does not exactly balance out the tdt_vadv.
            term3 = rho_i*term3
            ! term3= term3*pf(1,1,k)/(rdgas*pt(1,1,k))
         endif 
      else
         if ( dzi > 0. ) then
            term3= cp_air * divf * (0.25*dzi**(4./3.) + zi(1,1)*dzi**(1./3.))
            term3= term3*pf(1,1,k)/(rdgas*pt(1,1,k))
         else
            term3=0.;
         end if
      end if
      !h1g

      frad(1,1,k+1)=70.*exp(-Qz1) + 22.*exp(-Qz2) + term3
   end do
   frad(1,1,1)=frad(1,1,2);
   do k=1,kdim
      dT_rad(1,1,k)=grav/cp_air*(frad(1,1,k+1)-frad(1,1,k))/dp(1,1,k)  
   end do
end if
tdt_lw = dT_rad
tdt_sw = 0.0

dT_vadv=0.; dT_adi=0.; dqv_vadv=0.; dql_vadv=0.; dqi_vadv=0.;dqa_vadv=0.;dqn_vadv=0.;
du_vadv=0.; dv_vadv=0.;

if (do_vadv) then
   call get_time(dt_int,dt_seconds,dt_days)
   dts = real(dt_seconds + 86400*dt_days)

   dT_adi=rdgas*pt(:,:,:)*omga(:,:,:)/cp_air/pf(:,:,:)

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
      call vert_advection(dts,omega_h,dp,q(:,:,:,nsphum),dqv_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)
      call vert_advection(dts,omega_h,dp,q(:,:,:,nql),dql_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)
      call vert_advection(dts,omega_h,dp,q(:,:,:,nqi),dqi_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)
      call vert_advection(dts,omega_h,dp,q(:,:,:,nqa),dqa_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)
      if( nqn > 0 ) &
      call vert_advection(dts,omega_h,dp,q(:,:,:,nqn),dqn_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)


   select case (momentum_vert_advec_scheme)
   case(1)
      call vert_advection(dts,omega_h,dp,ua,du_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
      call vert_advection(dts,omega_h,dp,va,dv_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
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

end if

!--- ggg
!write(6,*) 'ggg, t_dt_in', t_dt
!write(6,*) 'ggg, dT_rad', dT_rad
!write(6,*) 'ggg, dT_vadv', dT_vadv

!<--- yhc, 2022-08-06
if (do_read_rf01_forc_any) then 
  t_dt = t_dt + dT_rad + tdt_any
  q_dt(:,:,:,nsphum) = qvdt_any

  tdt_dyn_forc=0.; qvdt_dyn_forc=0.
  tdt_dyn_forc  (:,:,:) = tdt_any
  qvdt_dyn_forc (:,:,:) = qvdt_any

  q_dt(:,:,:,nql) = dql_vadv
  q_dt(:,:,:,nqi) = dqi_vadv
  q_dt(:,:,:,nqa) = dqa_vadv
  if(nqn > 0) q_dt(:,:,:,nqn) = dqn_vadv

else  ! 
  dT_vadv= dT_vadv + dT_adi
  t_dt = t_dt +  dT_rad + dT_vadv + dT_lf  ! adding large-scale term
  q_dt(:,:,:,nsphum) = dqv_vadv   + dqv_lf ! adding large-scale term
  q_dt(:,:,:,nql) = dql_vadv
  q_dt(:,:,:,nqi) = dqi_vadv
  q_dt(:,:,:,nqa) = dqa_vadv
  if(nqn > 0) q_dt(:,:,:,nqn) = dqn_vadv

  tdt_dyn_forc=0.; qvdt_dyn_forc=0.
  tdt_dyn_forc  (:,:,:) = dT_vadv  + dT_lf
  qvdt_dyn_forc (:,:,:) = dqv_vadv + dqv_lf

end if   ! end if of do_any_forcings

!---> yhc, 2022-08-06

u_dt = du_vadv+du_geos; v_dt = dv_vadv+dv_geos;

!--- ggg
!write(6,*) 'ggg, vert_adv, dts',dts
!write(6,*) 'ggg, vert_adv, omega_h',omega_h
!write(6,*) 'ggg, vert_adv, ph',ph
!write(6,*) 'ggg, vert_adv, dp',dp
!write(6,*) 'ggg, vert_adv, pt',pt
!write(6,*) 'ggg, vert_adv, qq',q(:,:,:,nsphum)
!write(6,*) 'ggg, vert_adv, dT_vadv',dT_vadv
!write(6,*) 'ggg, vert_adv, dqv_vadv',dqv_vadv


if (id_pf_forc > 0)  used = send_data( id_pf_forc,  pf      (:,:,:), time_diag, 1, 1)

if (id_ph_forc > 0)  used = send_data( id_ph_forc,  ph      (:,:,:), time_diag, 1, 1)

if (id_zf_forc > 0)  used = send_data( id_zf_forc,  zf      (:,:,:), time_diag, 1, 1)

if (id_zh_forc > 0)  used = send_data( id_zh_forc,  zh      (:,:,:), time_diag, 1, 1)

if (id_tdt_vadv > 0) used = send_data( id_tdt_vadv, dT_vadv (:,:,:), time_diag, 1, 1)

if (id_udt_vadv > 0) used = send_data( id_udt_vadv, du_vadv (:,:,:), time_diag, 1, 1)

if (id_vdt_vadv > 0) used = send_data( id_vdt_vadv, dv_vadv (:,:,:), time_diag, 1, 1)

if (id_udt_geos > 0) used = send_data( id_udt_geos, du_geos (:,:,:), time_diag, 1, 1)

if (id_vdt_geos > 0) used = send_data( id_vdt_geos, dv_geos (:,:,:), time_diag, 1, 1)

if (id_qvdt_vadv> 0) used = send_data( id_qvdt_vadv,dqv_vadv(:,:,:), time_diag, 1, 1)

if (id_qldt_vadv> 0) used = send_data( id_qldt_vadv,dql_vadv(:,:,:), time_diag, 1, 1)

if (id_qidt_vadv> 0) used = send_data( id_qidt_vadv,dqi_vadv(:,:,:), time_diag, 1, 1)

if (id_qadt_vadv> 0) used = send_data( id_qadt_vadv,dqa_vadv(:,:,:), time_diag, 1, 1)


! ---> h1g, 2010-09-27
   qvdt_forcing_col = 0.0
   do k=1, kdim
         qvdt_forcing_col =  qvdt_forcing_col &
            + ( dqv_vadv( :,:,k )  ) * ( ph( :,:,k+1 ) -  ph( :,:,k ) ) / grav
   end do

   qldt_vadv_col = 0.0
   do k=1, kdim
         qldt_vadv_col =  qldt_vadv_col &
            + (  dql_vadv( :,:,k )  ) * ( ph( :,:,k+1 ) -  ph( :,:,k ) ) / grav
   end do
! <--- h1g, 2010-09-27


if( nqn > 0 ) then
  if (id_qndt_vadv> 0) used = send_data( id_qndt_vadv,dqn_vadv(:,:,:), time_diag, 1, 1)
endif

if (id_tdt_radf > 0) used = send_data( id_tdt_radf, dT_rad  (:,:,:), time_diag, 1, 1)

if (id_flx_radf > 0) used = send_data( id_flx_radf, frad    (:,:,:), time_diag, 1, 1)

if (id_udt_lf > 0)   used = send_data( id_udt_lf,   du_lf   (:,:,:), time_diag, 1, 1)

if (id_vdt_lf > 0)   used = send_data( id_vdt_lf,   dv_lf   (:,:,:), time_diag, 1, 1)

if (id_tdt_lf > 0)   used = send_data( id_tdt_lf,   dT_lf   (:,:,:), time_diag, 1, 1)

if (id_qvdt_lf > 0)  used = send_data( id_qvdt_lf,  dqv_lf  (:,:,:), time_diag, 1, 1)

if (id_zi_forc > 0)  used = send_data( id_zi_forc,  zi      (:,:),   time_diag, 1, 1)

if (id_u_geos > 0)   used = send_data( id_u_geos,   u_geos  (:,:,:), time_diag, 1, 1)

if (id_v_geos > 0)   used = send_data( id_v_geos,   v_geos  (:,:,:), time_diag, 1, 1)

! ---> h1g, 2010-09-27
if ( id_qvdt_forc_col > 0 )  used = send_data(  id_qvdt_forc_col, qvdt_forcing_col(:,:), time_diag, 1, 1 )
if ( id_qldt_vadv_col > 0 )  used = send_data(  id_qldt_vadv_col, qldt_vadv_col(:,:), time_diag, 1, 1 )
! <--- h1g, 2010-09-27

!<-- yhc 2022-08-06
if (id_tdt_dyn_forc > 0)  used = send_data( id_tdt_dyn_forc,  tdt_dyn_forc (:,:,:), time_diag, 1, 1)
if (id_qvdt_dyn_forc > 0) used = send_data( id_qvdt_dyn_forc, qvdt_dyn_forc (:,:,:), time_diag, 1, 1)
!--> yhc 2022-08-06


end subroutine update_rf01_forc

!########################################################################
! This subroutine adds longwave radiative heating from forcings
! to input arrays

subroutine add_rf01_tdtlw( x )

  implicit none
  real, intent(inout) :: x(:,:,:)

  if (allocated(tdt_lw)) x = x + tdt_lw

end subroutine add_rf01_tdtlw

!########################################################################
! This subroutine adds shortwave radiative heating from forcings
! to input arrays

subroutine add_rf01_tdtsw( x )

  implicit none
  real, intent(inout) :: x(:,:,:)

  if (allocated(tdt_sw)) x = x + tdt_sw

end subroutine add_rf01_tdtsw

!########################################################################
! This subroutine returns imposed surface fluxes

subroutine get_rf01_flx( ustar, flux_t, flux_q )

  implicit none
  real, intent(out), dimension(:) :: ustar, flux_t, flux_q

  ! ZNT 05/19/2020: Here we fix ustar = 0.25 which is not specified in 
  !                 Stevens (2005), but it is consistent with SCAM. 
  !                 Stevens (2005) specifies C_D = 0.0011, but this is not
  !                 well defined because the lowest air level is not given. 
  ustar = ustar_sfc
  flux_t = flux_t_sfc
  flux_q = flux_q_sfc

end subroutine get_rf01_flx
 !########################################################################


!########################################################################
subroutine get_rf01_flx_online( um_sfc,   vm_sfc,   thlm_sfc, rtm_sfc, rho, &
                                            flux_t, flux_q,  ustar)
implicit none
real, parameter ::  &
  ubmin = 1.0, &
  SST   = 292.5, &
  psfc  =   1017.e2,  &
  Cd    = 0.0011

! Input variables
real, intent(in),  dimension(:)  ::  &
  um_sfc,       & ! The upward flux of u-momentum         [(m^2 s^-2]
  vm_sfc,     & ! The Upward flux of v-momentum         [(m^2 s^-2]
  thlm_sfc,    & ! Theta_l at zt(2)      [K]
  rtm_sfc,     &         ! rt at zt(2)           [kg/kg]
  rho

! Output variables
real, intent(out),  dimension(:)  ::  &
  flux_t,  & ! w'theta_l' surface flux   [(m K)/s]
  flux_q,   &      ! w'rt' surface flux        [(m kg)/(kg s)]
  ustar

! Internal variables
  real ::    qsat

  ! Internal variables
  real, dimension(size(um_sfc))  :: &
    ubar   ! This is root (u^2 + v^2), per ATEX and RICO spec.

! ------------------------------------------------------------------------------------------
     ubar = max(ubmin, sqrt(um_sfc*um_sfc + vm_sfc*vm_sfc))
     call compute_qs( SST, psfc, qsat )
      flux_t = - Cd  * ubar * ( thlm_sfc - SST * (p00/psfc)**kappa ) * rho * cp_air
      flux_q  = - Cd  * ubar * ( rtm_sfc -  qsat ) * rho * hlv
      ustar = ustar_sfc
  return
  
end subroutine  get_rf01_flx_online
 !########################################################################



!########################################################################
! Subroutine retuns RF01 initial sounding

       subroutine rf01_snd( z, u, v, T, qv, ql )
       implicit none

       real, intent(in) :: z
       real, intent(out) :: u, v, T, qv, ql

       integer k
       real a, b

!      Constant wind profile

       u =  7.0
       v = -5.5

!      Interpolate to get T, qv, ql

       call locate( z_snd, ksnd, z, k )
       if ( k < 1 .or. k >= ksnd ) then
         call error_mesg('rf01_snd',  &
                         'z value ouf of sounding bounds', FATAL)
       endif

       b = ( z - z_snd(k) ) / ( z_snd(k+1) - z_snd(k) )
       a = 1.0 - b
       T  = a * T_snd(k)  + b * T_snd(k+1)
       qv = ( a * qv_snd(k) + b * qv_snd(k+1) ) * 0.001
       ql = ( a * ql_snd(k) + b * ql_snd(k+1) ) * 0.001

       return
       end subroutine rf01_snd


! ZNT 05/18/2020: Subroutine returns RF01 initial sounding as in Stevens et al (2005)
! The extended profile above 1600m is similar to that of CLUBB SCM

       subroutine rf01_snd_stevens( z, u, v, thetal, qt )
       implicit none

       real, intent(in) :: z
       real, intent(out) :: u, v, thetal, qt

!      Constant wind profile

       u =  7.0
       v = -5.5

!      qt (specific total water content)

       if ( z < zi_stevens ) then
         qt = 9.0
       else if ( z < 2400.0 ) then
         qt = 1.5
       else
         qt = max(0.0, 1.5 - 1.0e-3 * (z - 2400.))
       end if

!      convert units to kg/kg

       qt = 0.001 * qt

!      convert qt from specific humidity to mixing ratio

       qt = qt / ( 1.0 - qt )

!      thetal

       if ( z < zi_stevens ) then
         thetal = 289.0
       else if ( z < 1600.0 ) then
         thetal = 297.5 +  (z - zi_stevens)**(1./3.)
       else
         thetal = 297.5 + (760.)**(1./3.) +  6.0E-3 * (z - 1600.)
       end if

       return
       end subroutine rf01_snd_stevens



!########################################################################
! ZNT 03/29/2021: Initialize aerosol profiles
subroutine rf01_aer_init(kdim)

#include "fv_arrays.h"
#include "fv_point.inc"

integer, intent(in)  :: kdim
character(len=64) :: aer_nc='INPUT/tracer_rf01.nc'
integer :: unit, logunit

integer, dimension(21) :: ntracer
integer :: iz, itrac

real, allocatable, dimension(:,:,:,:) :: buffer_aer


if (allocated(buffer_aer)) deallocate (buffer_aer)
allocate(buffer_aer(1,1,kdim,21)); buffer_aer = 0.0

if( file_exist(trim(aer_nc)) ) then
    call mpp_open( unit, trim(aer_nc), action=MPP_RDONLY, &
             form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE )

    ! Variable names in the field table: 
    ! soa, dust1, dust2, dust3, dust4, dust5, simpledms, simpleso2, simpleso4, simplemsa, simpleh2o2,
    ! ssalt1, ssalt2, ssalt3, ssalt4, ssalt5, bcphob, bcphil, omphob, omphil, radon
    ntracer(1)  = get_tracer_index(MODEL_ATMOS, 'soa')
    ntracer(2)  = get_tracer_index(MODEL_ATMOS, 'dust1')
    ntracer(3)  = get_tracer_index(MODEL_ATMOS, 'dust2')
    ntracer(4)  = get_tracer_index(MODEL_ATMOS, 'dust3')
    ntracer(5)  = get_tracer_index(MODEL_ATMOS, 'dust4')
    ntracer(6)  = get_tracer_index(MODEL_ATMOS, 'dust5')
    ntracer(7)  = get_tracer_index(MODEL_ATMOS, 'simpledms')
    ntracer(8)  = get_tracer_index(MODEL_ATMOS, 'simpleso2')
    ntracer(9)  = get_tracer_index(MODEL_ATMOS, 'simpleso4')
    ntracer(10) = get_tracer_index(MODEL_ATMOS, 'simplemsa')
    ntracer(11) = get_tracer_index(MODEL_ATMOS, 'simpleh2o2')
    ntracer(12) = get_tracer_index(MODEL_ATMOS, 'ssalt1')
    ntracer(13) = get_tracer_index(MODEL_ATMOS, 'ssalt2')
    ntracer(14) = get_tracer_index(MODEL_ATMOS, 'ssalt3')
    ntracer(15) = get_tracer_index(MODEL_ATMOS, 'ssalt4')
    ntracer(16) = get_tracer_index(MODEL_ATMOS, 'ssalt5')
    ntracer(17) = get_tracer_index(MODEL_ATMOS, 'bcphob')
    ntracer(18) = get_tracer_index(MODEL_ATMOS, 'bcphil')
    ntracer(19) = get_tracer_index(MODEL_ATMOS, 'omphob')
    ntracer(20) = get_tracer_index(MODEL_ATMOS, 'omphil')
    ntracer(21) = get_tracer_index(MODEL_ATMOS, 'radon')

    ! Variable names in the input NC file:
    ! SOA, dust1, dust2, dust3, dust4, dust5, DMS, SO2, SO4, MSA, H2O2,
    ! ssalt1, ssalt2, ssalt3, ssalt4, ssalt5, bcphob, bcphil, omphob, omphil, radon
    call read_data(aer_nc, 'SOA',    buffer_aer(:,:,:,1),  no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'dust1',  buffer_aer(:,:,:,2),  no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'dust2',  buffer_aer(:,:,:,3),  no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'dust3',  buffer_aer(:,:,:,4),  no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'dust4',  buffer_aer(:,:,:,5),  no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'dust5',  buffer_aer(:,:,:,6),  no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'DMS',    buffer_aer(:,:,:,7),  no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'SO2',    buffer_aer(:,:,:,8),  no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'SO4',    buffer_aer(:,:,:,9),  no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'MSA',    buffer_aer(:,:,:,10), no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'H2O2',   buffer_aer(:,:,:,11), no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'ssalt1', buffer_aer(:,:,:,12), no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'ssalt2', buffer_aer(:,:,:,13), no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'ssalt3', buffer_aer(:,:,:,14), no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'ssalt4', buffer_aer(:,:,:,15), no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'ssalt5', buffer_aer(:,:,:,16), no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'bcphob', buffer_aer(:,:,:,17), no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'bcphil', buffer_aer(:,:,:,18), no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'omphob', buffer_aer(:,:,:,19), no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'omphil', buffer_aer(:,:,:,20), no_domain=.true., timelevel=1)
    call read_data(aer_nc, 'radon',  buffer_aer(:,:,:,21), no_domain=.true., timelevel=1)

    do iz = 1, kdim
      do itrac = 1,21
        q(:,:,iz,ntracer(itrac)) = buffer_aer(1,1,iz,itrac)
      enddo
    enddo
    ! write(*,*) 'SOA',    q(:,:,:,ntracer(1))
    ! write(*,*) 'ssalt1', q(:,:,:,ntracer(12))


    call mpp_close( unit )

    if(mpp_pe() == mpp_root_pe() ) then
        logunit =stdlog()
        write (logunit,'(a)') 'used AEROSOL profile: '// aer_nc
        call close_file(logunit)
    endif

end if

end subroutine rf01_aer_init

!########################################################################

! ZNT 02/20/2020: Note - Existing code for calculating z
!
! ------ atmos_fv_dynamics/tools/fv_diagnostics.F90, Line 653 ------
! Compute height at layer edges
!         do i=1,im
!            wz(i,j,km+1) = phis(i,j) * ginv
!         enddo
!
!         do k=km,1,-1
!            do i=1,im
! #ifdef SW_DYN
!                wz(i,j,k) = wz(i,j,k+1) + gg*(pk(i,j,k+1)-pk(i,j,k))
! #else
!                wz(i,j,k) = wz(i,j,k+1) + gg*pt(i,j,k)*(1.+zvir*q(i,j,k,1))   &
!                     *(peln(i,k+1,j)-peln(i,k,j))
! #endif
!             enddo
!          enddo
!       end do
!
! ------ atmos_scm/driver/coupled/atmosphere.F90, fv_compute_p_z ------
!  if (hydrostatic ) then
!       do k=npz,1,-1
!         do j=1,size(phis,2)
!           do i=1,size(phis,1)
!             tvm = rrg*pt(i,j,k)*(1.+zvir*q_sph(i,j,k))
!             p_full(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
!             z_full(i,j,k) = z_half(i,j,k+1) + tvm*(1.-p_half(i,j,k)/p_full(i,j,k))
!             z_half(i,j,k) = z_half(i,j,k+1) + tvm*(peln(i,k+1,j)-peln(i,k,j))
!           enddo
!         enddo
!       enddo
!  endif


subroutine compute_p_z (npz, phis, pe, peln, delp, pt, q_sph, p_full, p_half, z_full, z_half)
  ! ZNT 02/20/2020: adapted from the hydrostatic version of 'fv_compute_p_z' subroutine in atmosphere.F90
     integer, intent(in)  :: npz
     real, dimension(:,:),   intent(in)  :: phis
     real, dimension(:,:,:), intent(in)  :: pe, peln, delp, pt, q_sph
     real, dimension(:,:,:), intent(out) :: p_full, p_half, z_full, z_half
  !--- local variables
     integer i,j,k
     real    tvm
     real    :: zvir, rrg, ginv

     zvir = rvgas/rdgas - 1.
     ginv = 1./ grav
     rrg  = rdgas / grav

  !----------------------------------------------------
  ! Compute pressure and height at full and half levels
  !----------------------------------------------------
     z_half(:,:,npz+1) = phis(:,:) * ginv

     do k=1,npz+1
       do j=1,size(phis,2)
          do i=1,size(phis,1)
            p_half(i,j,k) = pe(i,k,j)
         enddo
       enddo
     enddo

     do k=npz,1,-1
       do j=1,size(phis,2)
         do i=1,size(phis,1)
           tvm = rrg*pt(i,j,k)*(1.+zvir*q_sph(i,j,k))
           p_full(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
           z_full(i,j,k) = z_half(i,j,k+1) + tvm*(1.-p_half(i,j,k)/p_full(i,j,k))
           z_half(i,j,k) = z_half(i,j,k+1) + tvm*(peln(i,k+1,j)-peln(i,k,j))
         enddo
       enddo
     enddo

end subroutine compute_p_z

!###################################

!<-- yhc, 2022-08-06
subroutine rf01_snd_any_profiles(u_rf01, v_rf01, T_rf01, qv_rf01, qa_rf01, ql_rf01)
  !=================================
  ! Description
  !    read specific profiles as RF01 initial conditions
  !=================================

  !-------------------
  ! Output argument
  !   u_rf01: zonal wind speed  (m/s)
  !   v_rf01: meridional wind speed (m/s)
  !   T_rf01: temperature (K) 
  !   qv_rf01: specific humidity (kg/kg) 
  !   qa_rf01: cloud fraction (0-1, unitless) 
  !   ql_rf01: cloud liquid specific humidity (kg/kg)  
  !-------------------
  real, intent(out), dimension(:) :: &
    u_rf01, v_rf01, T_rf01, qv_rf01, qa_rf01, ql_rf01

  !------------------
  ! local argument
  !------------------

!------------------------------------

  !--- initialize
  u_rf01=0.; v_rf01=0.; T_rf01=0.; qv_rf01=0.; qa_rf01=0.; ql_rf01=0.

  !--- set u and v following subsoutine rf01_snd_stevens
  u_rf01(:) = 7.0
  v_rf01(:) = -5.5

  !--- no cloud liquid
  ql_rf01(:) = 0.

  !--- set temperature and specific humidity
  if (trim(option_any_profiles) .eq. "test") then
    T_rf01 (:) = 290.
    qv_rf01(:) = 8.e-3

    !call error_mesg('rf01_snd_any_profiles',  &
    !                'test profile. STOP', FATAL)

  elseif (trim(option_any_profiles) .eq. "nudgeAM4_DYCOMS_3havg_10Jul1030Z") then
    !--- T & Q profiles in nudge AM4 in DYCOMS region, 3-hourly avg fileds at 1030 UTC, July 10, 2001

    T_rf01(:) = (/258.74707, 239.25804, 231.56488, 226.97307, 223.35747, 219.15889, &
        214.8654 , 209.97337, 205.57759, 204.78589, 211.66156, 219.53874, &
        223.87708, 230.67662, 240.41756, 251.41612, 261.01166, 268.3155 , &
        275.00217, 280.10532, 283.84042, 287.29593, 289.93237, 291.9652 , &
        292.96103, 287.9403 , 285.0339 , 286.06296, 287.34573, 288.40512, &
        289.24265, 289.8793 , 290.34232/)

    qv_rf01(:) = (/1.7768771e-06, 1.8155506e-06, 1.8678315e-06, 1.8763882e-06, &
        1.8619512e-06, 1.8223737e-06, 1.6822354e-06, 1.6580552e-06, &
        2.4243193e-06, 2.3362234e-06, 2.4666645e-06, 4.0264349e-06, &
        1.1465544e-05, 5.2478834e-05, 3.1549469e-04, 6.8929815e-04, &
        1.1172843e-03, 1.7424213e-03, 1.8735757e-03, 3.1218224e-03, &
        4.6328190e-03, 5.0644558e-03, 4.6861884e-03, 4.5053461e-03, &
        3.6602162e-03, 3.9170138e-03, 7.6311412e-03, 8.7521821e-03, &
        8.9964895e-03, 9.0813302e-03, 9.1430033e-03, 9.2054578e-03, &
        9.3241045e-03/)

  elseif (trim(option_any_profiles) .eq. "nudgeAM4_DYCOMS_inst_10Jul09Z") then
    !--- T & Q profiles in nudge AM4 in DYCOMS region, instantaneous fileds at 09 UTC, July 10, 2001
    T_rf01(:) = (/259.56537, 239.6206 , 231.85065, 227.21313, 223.41449, 219.25061,  &
       214.79808, 209.90749, 205.66545, 204.83875, 211.55298, 219.46045, &
       223.6038 , 230.77509, 240.62654, 251.53851, 261.03146, 268.67322, &
       275.37033, 280.208  , 283.89157, 287.3503 , 290.03568, 292.02527, &
       293.2338 , 288.3329 , 285.20035, 286.09076, 287.37228, 288.43158, &
       289.269  , 289.9053 , 290.36707/)
    qv_rf01(:) = (/1.7761973e-06, 1.8171096e-06, 1.8693468e-06, 1.8782731e-06,  &
       1.8628867e-06, 1.8269596e-06, 1.6796785e-06, 1.6591595e-06, &
       2.4083802e-06, 2.2477245e-06, 2.3656903e-06, 3.5571929e-06, &
       1.2500092e-05, 5.2879685e-05, 3.4244670e-04, 6.8660267e-04, &
       1.1139525e-03, 1.6567194e-03, 1.7963117e-03, 3.2657601e-03, &
       4.7457293e-03, 5.1978002e-03, 4.8223357e-03, 4.5640515e-03, &
       3.7170404e-03, 3.9796312e-03, 7.4526053e-03, 8.6870762e-03, &
       8.9330534e-03, 9.0221390e-03, 9.0844352e-03, 9.1471216e-03, &
       9.2648193e-03/)

  elseif (trim(option_any_profiles) .eq. "MERRA2_DYCOMS_3hr_10Jul1030Z") then
    !--- T & Q profiles in MERRA-2 in DYCOMS region, 3-hourly average fileds at 1030 UTC, July 10, 2001

   T_rf01(:) = (/256.13543586, 241.21235689, 232.25942617, 226.71262839, 222.36388726, &
 217.91831726, 213.50071576, 209.70515883, 204.97553234, 204.08718919, &
 211.45008232, 221.41861183, 225.34816312, 231.10761687, 240.60936584, &
 252.28905087, 261.34257744, 267.58192342, 274.2622736, 280.17972151, &
 284.21902098, 287.44686975, 290.03266595, 291.91502186, 291.9165224, &
 286.29239834, 284.79288802, 286.05984094, 287.36214758, 288.42562299, &
 289.27423139, 289.90759277, 289.90759277 /)
 
  qv_rf01(:) = (/3.79629153e-06, 3.43418906e-06, 3.17845803e-06, 3.00516667e-06, &
 2.92830573e-06, 2.85784785e-06, 2.75097319e-06, 2.65936623e-06, &
 2.81214308e-06, 3.39195961e-06, 3.57775655e-06, 6.08632194e-06, &
 1.45774457e-05, 4.96127500e-05, 1.58722538e-04, 4.00973008e-04, &
 5.70193068e-04, 9.72203468e-04, 1.12979584e-03, 1.34697952e-03, &
 1.75247570e-03, 2.24887987e-03, 2.84085824e-03, 3.68561834e-03, &
 4.85196994e-03, 7.32555364e-03, 8.68808131e-03, 8.96757626e-03, &
 9.08925645e-03, 9.16563727e-03, 9.25226324e-03, 9.32184514e-03, &
 9.32184514e-03 /)

  elseif (trim(option_any_profiles) .eq. "MERRA2_DYCOMS_3havg_10Jul1030Z_RF01") then
    !--- T & Q profiles in MERRA-2 in DYCOMS region, 3-hourly average fileds at 1030 UTC, July 10, 2001
    !    2022-11-27 add

    T_rf01(:) = &
      (/256.13368863, 241.21046744, 232.2577286 , 226.71260078, &
       222.36269731, 217.91733437, 213.49985416, 209.70254965, &
       204.97370929, 204.08261604, 211.44499554, 221.41427917, &
       225.35195344, 231.11066055, 240.61084672, 252.28921403, &
       261.34279018, 267.58367671, 274.26340006, 280.17848477, &
       284.21688651, 287.44532522, 290.03194201, 291.91423907, &
       291.92165924, 286.30007853, 284.79525246, 286.06155979, &
       287.3638988 , 288.4274641 , 289.27614432, 289.90971236, &
       289.90971236/)

    qv_rf01(:) = &
     (/3.79621829e-06, 3.43411190e-06, 3.17834379e-06, 3.00508828e-06, &
       2.92823285e-06, 2.85778077e-06, 2.75093258e-06, 2.65932996e-06, &
       2.81219639e-06, 3.39190113e-06, 3.57773848e-06, 6.08497312e-06, &
       1.45712937e-05, 4.96102363e-05, 1.58759977e-04, 4.01234568e-04, &
       5.70596927e-04, 9.73605775e-04, 1.13124655e-03, 1.34960075e-03, &
       1.75667966e-03, 2.25285912e-03, 2.84432266e-03, 3.69341833e-03, &
       4.86087000e-03, 7.32719245e-03, 8.68883566e-03, 8.96865425e-03, &
       9.09039778e-03, 9.16679204e-03, 9.25340033e-03, 9.32298580e-03, &
       9.32298580e-03/)

    qa_rf01(:) = &
     (/0.        , 0.        , 0.        , 0.        , 0.        , &
       0.        , 0.        , 0.        , 0.        , 0.        , &
       0.        , 0.        , 0.        , 0.        , 0.        , &
       0.        , 0.        , 0.        , 0.        , 0.        , &
       0.        , 0.        , 0.        , 0.        , 0.00048463, &
       0.18675185, 0.147899  , 0.00024432, 0.        , 0.        , &
       0.        , 0.        , 0. /)

    ql_rf01(:) = &
     (/0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       5.84131941e-07, 3.12382097e-05, 3.29937776e-05, 1.32037487e-05, &
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, &
       0.00000000e+00/)

  elseif (trim(option_any_profiles) .eq. "nudgeAM4_DYCOMS_3havg_10Jul1030Z_RF01") then
    !--- T & Q profiles in nudge AM4 in DYCOMS region, 3-hourly avg fileds at 1030 UTC, July 10, 2001 (RF01 period)
    !    2022-11-27 add

   T_rf01(:) =  (/258.74643792, 240.60091205, 231.94549015, 227.15072919, &
       223.47576725, 219.28556651, 214.98620826, 210.10304275, &
       205.68791231, 204.80240676, 211.50121592, 219.36495816, &
       223.78829474, 230.54268129, 240.23320863, 251.22000625, &
       260.85520862, 268.20920814, 274.91364207, 280.04053853, &
       283.79303711, 287.25064327, 289.89551294, 291.93276848, &
       292.94559406, 288.04685467, 285.10582577, 286.03403882, &
       287.30124507, 288.35955338, 289.19690631, 289.83314288, &
       290.29411126/)

   qv_rf01(:) = (/ &
       1.77687855e-06, 1.81290892e-06, 1.86528272e-06, 1.87604803e-06, &
       1.86239377e-06, 1.82351660e-06, 1.68620747e-06, 1.65868526e-06, &
       2.40507274e-06, 2.33851854e-06, 2.46405425e-06, 3.99077121e-06, &
       1.13068233e-05, 5.16358502e-05, 3.10407906e-04, 6.82554970e-04, &
       1.11030508e-03, 1.73348452e-03, 1.87271328e-03, 3.10708080e-03, &
       4.61512798e-03, 5.06027303e-03, 4.69440784e-03, 4.51175333e-03, &
       3.67787137e-03, 3.91129016e-03, 7.54075477e-03, 8.71966380e-03, &
       8.98804491e-03, 9.07782667e-03, 9.13978769e-03, 9.20106090e-03, &
       9.31171610e-03 /)

  else
    call error_mesg('rf01_snd_any_profiles',  &
                    'The input option_any_profiles is not supported', FATAL)

  endif

end subroutine rf01_snd_any_profiles
!--> yhc, 2022-08-06

!###################################

!<-- yhc, 2022-08-06
subroutine read_rf01_forc_any(tdt_any, qvdt_any, omega_any)
  !=================================
  ! Description
  !    read specific forcings profiles, e.g. those from AM4
  !=================================

  !-------------------
  ! Output argument
  !   tdt_any : temperature tendency (K/s)
  !   qvdt_any: specific humidity tendency (kg/kg/s)
  !-------------------
  real, intent(out), dimension(:,:,:) :: &
    tdt_any, qvdt_any, omega_any 

  !-------------------
  ! local argument
  !-------------------
  integer i,j,k, ix, jx, kx

!---------------------------
  ix = size(tdt_any,1)
  jx = size(tdt_any,2)
  kx = size(tdt_any,3)

  !--- initialize
  tdt_any=0.;  qvdt_any=0. ; omega_any=0.

  !!!!! replace ',' to ', &'.  A,Bs/,$/, \&/gc

  !--- set temperature and specific humidity tendencies
  if (trim(option_read_rf01_forc_any) .eq. "test") then
    do i=1,ix
    do j=1,jx
      tdt_any (i,j,:)  = 1./86400.       ! 1 K/day
      qvdt_any(i,j,:)  = 1.e-3 / 86400.  ! 1 g/kg/day
      omega_any(i,j,:) = 40.*100./86400. ! 40 hPa/day
    enddo
    enddo

  elseif (trim(option_read_rf01_forc_any) .eq. "nudgeAM4_DYCOMS_3havg_10Jul1030Z_tdtDyn") then
    !--- T & Q tendencies and omega in nudge AM4 in DYCOMS region, 3-hourly fileds at 1030 UTC, July 10, 2001
    !    T tendency consist of dynamics, no nudging 
    do i=1,ix
    do j=1,jx

      tdt_any (i,j,:)  = &
        (/-1.63755994e-05, -9.97397819e-06, -1.84662258e-05, &
        -2.46192940e-05,  1.49949019e-05,  8.32301066e-06, &
         7.32212948e-06,  1.85725348e-05, -7.35925505e-06, &
        -5.91633398e-06,  2.42959741e-05, -5.32737567e-06, &
         3.92294132e-05,  3.42931958e-06,  3.83025827e-06, &
        -7.09798087e-06,  1.21614057e-05, -1.20696905e-05, &
        -1.60586660e-05,  8.64245249e-06,  8.77882758e-06, &
         1.20371114e-05,  2.41232283e-06,  1.97684767e-05, &
         2.99433505e-05,  2.10677772e-05,  7.77955847e-06, &
        -2.06766417e-05, -1.49885227e-05, -1.62530632e-05, &
        -1.60085619e-05, -1.60602794e-05, -1.57481645e-05  /)

      qvdt_any (i,j,:)  = &
        (/8.28537271e-14, -1.55279164e-13, -2.23506652e-13, &
        -2.90465650e-13, -8.16109672e-14, -7.37974189e-13, &
         3.42449159e-13, -3.53970542e-13,  1.97936871e-12, &
         1.20599685e-11,  1.77701014e-11,  6.16989515e-11, &
        -9.35818159e-11, -6.50476409e-11, -4.77419038e-09, &
         2.47446286e-10,  8.17462986e-10,  1.03181685e-08, &
         1.02892033e-08, -2.58301238e-08, -2.25005081e-08, &
        -2.49675889e-08, -2.63582987e-08, -1.43745362e-08, &
        -1.09896767e-08, -1.65990475e-08, -3.70988680e-08, &
        -1.13442944e-08, -1.33867912e-08, -1.18559988e-08, &
        -1.19625403e-08, -1.20990729e-08, -1.15211094e-08 /)

      omega_any (i,j,:)  = &
         (/-8.2690676e-05, -1.3924930e-04, -3.4543127e-04, -3.9579478e-04, &
         4.7645532e-04,  4.7036645e-04,  1.0444897e-03,  3.1555819e-03, &
         2.0553563e-03,  1.5838591e-02,  4.2532876e-02,  4.7241226e-02, &
         2.2664821e-02, -7.8164358e-03, -2.2913152e-02, -6.3161091e-03, &
         2.8577944e-02,  5.7582229e-02,  7.7478357e-02,  8.7921083e-02, &
         9.2251174e-02,  8.9447342e-02,  8.4838212e-02,  8.2225561e-02, &
         8.0350965e-02,  7.4529313e-02,  6.1468210e-02,  4.4626404e-02, &
         2.9814288e-02,  1.8460358e-02,  9.8877270e-03,  3.6214052e-03, &
        -6.2030775e-04 /)
    enddo   ! end loop of j
    enddo   ! end loop of i

  elseif (trim(option_read_rf01_forc_any) .eq. "nudgeAM4_DYCOMS_3havg_10Jul1030Z_tdtDyn_RF01") then
    !--- T & Q tendencies and omega in nudge AM4 in DYCOMS region, 3-hourly fileds at 1030 UTC, July 10, 2001
    !    T tendency consist of dynamics, no nudging
    !    2022-11-27 add 
    do i=1,ix
    do j=1,jx

      tdt_any (i,j,:)  = &
        (/-1.63558040e-05, -1.04333125e-05, -1.80526060e-05, -2.43540613e-05, &
        1.36958513e-05,  8.52546729e-06,  7.38162330e-06,  1.82219474e-05, &
       -6.69272614e-06, -5.94449713e-06,  2.35992861e-05, -4.69731398e-06, &
        3.82888843e-05,  4.15344211e-06,  3.83771965e-06, -6.89712382e-06, &
        1.18232430e-05, -1.17139560e-05, -1.59631337e-05,  8.34172538e-06, &
        8.79578629e-06,  1.20086355e-05,  2.57677225e-06,  1.94612200e-05, &
        2.98035013e-05,  2.12598222e-05,  8.11304147e-06, -1.98677476e-05, &
       -1.51980099e-05, -1.62056800e-05, -1.60297178e-05, -1.60634939e-05, &
       -1.57875239e-05/)

      qvdt_any (i,j,:)  = &
        (/8.28025796e-14, -1.39134530e-13, -2.21113083e-13, -2.88064544e-13, &
       -8.85638651e-14, -7.21227637e-13,  3.12029917e-13, -3.34487806e-13, &
        1.92317944e-12,  1.17964728e-11,  1.76258101e-11,  6.06428431e-11, &
       -9.00929462e-11, -6.59704138e-11, -4.68444512e-09,  1.64237794e-10, &
        8.19167091e-10,  1.01557943e-08,  1.02638721e-08, -2.54025790e-08, &
       -2.25847643e-08, -2.49219144e-08, -2.63112398e-08, -1.45353748e-08, &
       -1.10773305e-08, -1.65244173e-08, -3.66524360e-08, -1.21091483e-08, &
       -1.33361237e-08, -1.19440339e-08, -1.19771443e-08, -1.21093130e-08, &
       -1.16021277e-08/)

      omega_any (i,j,:)  = &
        (/-8.26218031e-05, -1.35470637e-04, -3.35077260e-04, -3.92950026e-04, &
        4.48274300e-04,  4.71067753e-04,  1.02914267e-03,  3.09596788e-03, &
        2.08487332e-03,  1.55110879e-02,  4.19269430e-02,  4.71605311e-02, &
        2.31982317e-02, -7.19773224e-03, -2.26219468e-02, -6.62380858e-03, &
        2.80375004e-02,  5.72306147e-02,  7.72755045e-02,  8.78440540e-02, &
        9.22717099e-02,  8.95878576e-02,  8.50071472e-02,  8.23418004e-02, &
        8.04227699e-02,  7.46676925e-02,  6.17984352e-02,  4.51275010e-02, &
        3.03481395e-02,  1.89663193e-02,  1.03710865e-02,  4.08822655e-03, &
       -1.68477368e-04/)

    enddo   ! end loop of j
    enddo   ! end loop of i

  elseif (trim(option_read_rf01_forc_any) .eq. "nudgeAM4_DYCOMS_inst_10Jul09Z_tdtDpN") then
    !--- T & Q tendencies and omega in nudge AM4 in DYCOMS region, instantaneous fileds at 09 UTC, July 10, 2001
    !    T tendency consist of dynamics and nudging
    do i=1,ix
    do j=1,jx
      tdt_any (i,j,:)  = (/1.3144995e-06, -1.5953374e-05, -1.4663427e-05, -1.0413538e-05, &
       -4.7731737e-06, -1.0182043e-05,  2.5994472e-05,  1.2129881e-06, &
       -1.9836882e-05, -9.5418200e-06,  1.3509988e-05,  2.0307789e-05, &
        4.7310405e-05,  3.2322289e-06,  9.2622759e-07,  9.7163811e-06, &
        2.2719691e-05, -3.7756588e-05, -4.0724448e-05,  6.1348537e-06, &
        1.1783287e-05,  1.1151913e-05,  3.2755852e-06,  1.2919715e-05, &
       -8.0487262e-06, -4.1700496e-05, -3.3128190e-06, -1.8820676e-05, &
       -1.3733691e-05, -1.4773523e-05, -1.4161271e-05, -1.3770061e-05, &
       -1.3664116e-05/)
      qvdt_any(i,j,:)  = (/9.96605074e-14, -3.34396139e-13, -2.43753504e-13, -3.12022896e-13,  &
       -2.04319593e-13, -6.41629978e-13,  4.93797765e-13, -1.51837280e-14, &
        3.10330516e-12,  1.53420974e-11,  1.37990825e-11,  8.04331254e-11, &
       -2.14189999e-10, -1.09768028e-10, -3.63161701e-09,  7.19752868e-10, &
        2.83813195e-10,  1.64030869e-08,  1.36852973e-08, -1.91410354e-08, &
       -1.33786893e-08, -1.74636856e-08, -1.73407475e-08, -4.83570561e-09, &
       -7.51298668e-09, -2.25630430e-08, -4.28043450e-08, -1.95569587e-08, &
       -1.65458509e-08, -1.52714481e-08, -1.50549386e-08, -1.51398591e-08, &
       -1.43465586e-08/) 
      omega_any(i,j,:) = (/-4.2861629e-05, -1.8450012e-04, -2.4843027e-04, -1.4538245e-04, &
        1.2314529e-04,  3.8904420e-04,  1.8587558e-03,  1.5435455e-03, &
        1.6337573e-03,  1.9469643e-02,  4.0377978e-02,  3.7253466e-02, &
        1.6697088e-02, -5.4769162e-03, -1.7506778e-02, -6.2186136e-03, &
        2.1719443e-02,  4.1611195e-02,  5.1081970e-02,  6.2214240e-02, &
        7.6720290e-02,  8.4804989e-02,  8.4666051e-02,  8.2148612e-02, &
        7.9994023e-02,  7.4541010e-02,  6.1911773e-02,  4.4472393e-02, &
        2.8503373e-02,  1.5996391e-02,  6.4679310e-03, -5.3951336e-04, &
       -5.3057675e-03/)
    enddo   ! end loop of j
    enddo   ! end loop of i

  elseif (trim(option_read_rf01_forc_any) .eq. "MERRA2_DYCOMS_3havg_10Jul1030Z_RF01") then
    !--- T & Q profiles in MERRA-2 in DYCOMS region, 3-hourly average fileds at 1030 UTC, July 10, 2001
    !    2022-11-27 add
    do i=1,ix
    do j=1,jx
      tdt_any (i,j,:)  = &
        (/9.44607347e-06, -9.52326739e-06, -6.03630776e-05, -3.26233477e-05, &
        1.75147497e-05,  5.55060107e-05,  3.28256499e-05, -8.66710467e-06, &
       -1.09391627e-06,  9.84080226e-06, -1.57696264e-05,  1.60015662e-06, &
        2.99088145e-05,  2.12256133e-05, -8.54442185e-07, -1.62506290e-05, &
        2.36705046e-06, -1.19462430e-05, -2.45026529e-05, -8.43949536e-06, &
        9.19961600e-06,  1.12073810e-05,  9.93037855e-06,  8.71417910e-06, &
        8.00498779e-06,  8.87235235e-06, -1.14678007e-05, -1.48464145e-05, &
       -1.43275626e-05, -1.43990512e-05, -1.44389532e-05, -1.44389532e-05, &
       -1.44389532e-05/)

      qvdt_any (i,j,:)  = &
       (/2.86371900e-13,  8.81574542e-13, -8.39186219e-13, -3.67980619e-13, &
       -8.14417274e-14,  7.21071764e-14,  8.98667430e-14,  7.27019498e-14, &
       -3.75115548e-14, -5.55387902e-13, -1.73093676e-12, -1.88155157e-11, &
       -2.48616547e-11,  1.63579885e-10,  1.82205795e-09,  2.99042655e-09, &
       -2.92331630e-09,  1.15088388e-08,  1.55786544e-08,  1.56155047e-08, &
       -2.36106498e-09, -2.39255509e-08, -3.40012022e-08, -3.17843564e-08, &
       -3.09384324e-08, -1.79534882e-08, -1.73319111e-08, -2.02995458e-08, &
       -2.21108256e-08, -2.17014108e-08, -2.14728920e-08, -2.14728920e-08, &
       -2.14728920e-08/)

      omega_any (i,j,:)  = &
       (/-0.00012452, -0.00025997, -0.00076808, -0.00012822,  0.00066273, &
        0.0029815 ,  0.00268346,  0.00203249, -0.00386734,  0.02994247, &
        0.06928966,  0.06661236,  0.043951  , -0.00159991, -0.03658661, &
       -0.02580937,  0.0118829 ,  0.04512052,  0.08075962,  0.09527853, &
        0.09993955,  0.09256994,  0.08409546,  0.07632338,  0.0708227 , &
        0.06717046,  0.05440002,  0.03859051,  0.02534778,  0.01517401, &
        0.00742699,  0.00171494,  0.00171494/)

    enddo   ! end loop of j
    enddo   ! end loop of i

  elseif (trim(option_read_rf01_forc_any) .eq. "MERRA2_DYCOMS_3hr_10Jul1030Z") then
    !--- T & Q profiles in MERRA-2 in DYCOMS region, 3-hourly average fileds at 1030 UTC, July 10, 2001
    do i=1,ix
    do j=1,jx

tdt_any (i,j,:)  = (/9.52688447e-06, -9.53216593e-06, -6.03244882e-05, -3.26703694e-05, &
  1.74554512e-05,  5.55472071e-05,  3.28820802e-05, -8.58800158e-06, &
 -1.10264785e-06,  9.75809207e-06, -1.57884724e-05,  1.65871049e-06, &
  2.99583437e-05,  2.12050268e-05, -8.76309534e-07, -1.62874356e-05, &
  2.40745300e-06, -1.18866695e-05, -2.45735394e-05, -8.45710695e-06, &
  9.17615361e-06,  1.11620677e-05,  9.89498994e-06,  8.70885840e-06, &
  7.90669762e-06,  8.78167714e-06, -1.14120187e-05, -1.47730981e-05, &
 -1.42651468e-05, -1.43430519e-05, -1.43865354e-05, -1.43865354e-05, &
 -1.43865354e-05/)

qvdt_any (i,j,:)  = (/2.87791589e-13,  8.80895932e-13, -8.38098870e-13, -3.67820389e-13, &
 -8.15827370e-14,  7.25638179e-14,  9.08845989e-14,  7.24858672e-14, &
 -3.98414485e-14, -5.67121904e-13, -1.76491621e-12, -1.86259695e-11, &
 -2.47291576e-11,  1.62107009e-10,  1.82310461e-09,  3.00888050e-09, &
 -2.92863563e-09,  1.14272725e-08,  1.55235247e-08,  1.54998679e-08, &
 -2.49130957e-09, -2.39567700e-08, -3.39707366e-08, -3.17741164e-08, &
 -3.09201577e-08, -1.78716227e-08, -1.72616136e-08, -2.02300027e-08, &
 -2.20387823e-08, -2.16277938e-08, -2.13983959e-08, -2.13983959e-08, &
 -2.13983959e-08 /)

omega_any (i,j,:)  = (/-0.00012409, -0.00025993, -0.00076685, -0.00012839,  0.00065766,  0.00298474, &
  0.00268975,  0.00203718, -0.00386964,  0.02992577,  0.06928324,  0.06656368, &
  0.0439096,  -0.00159038, -0.03653991, -0.02570238,  0.01176625,  0.04492528, &
  0.08060706,  0.09511077,  0.09970001,  0.09225194,  0.0837952, 0.07611704, &
  0.07071545,  0.06710755,  0.05434628,  0.0385482, 0.02531743,  0.01515457, &
  0.00741664,  0.00171273,  0.00171273 /)
    enddo   ! end loop of j
    enddo   ! end loop of i

  else
    call error_mesg('read_rf01_forc_any',  &
                    'The input option_read_rf01_forc_any is not supported', FATAL)

  endif
  

end subroutine read_rf01_forc_any



end module scm_rf01_mod
 
