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

   public cgils_data_read, cgils_forc_init, cgils_forc_end, update_cgils_forc

   character(len=8) :: mod_name = 'scm_cgils'
   character(len=7) :: mod_name_diag = 'forcing'

!--------------------- version number ----------------------------------
!
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
integer, dimension(1) :: restart_versions = (/ 1 /)

!--- local variables
logical, private               :: initialized = .false.
integer                        :: nlev_cgils
integer                        :: time_cgils
integer, parameter  :: nlat_cgils = 1
integer, parameter  :: nlon_cgils = 1
real, allocatable, dimension(:,:,:) :: buffer_cgils_3D
real, allocatable, dimension(:,:)   :: buffer_cgils_2D
real, allocatable, dimension(:)     :: buffer_cgils_1D_lev
real, allocatable, dimension(:)     :: T_cgils, U_cgils, V_cgils, q_cgils, pfull_cgils   ! dimension (nlev)
real, allocatable, dimension(:)     :: Ps_cgils   ! dimension (ntime)

real,    private               :: missing_value = -999.

!--- namelist parameters
integer, public                :: tracer_vert_advec_scheme = 3
integer, public                :: temp_vert_advec_scheme = 3
integer, public                :: momentum_vert_advec_scheme = 3
character(len=20) , public     :: cgils_case    = "none"
character(len=200), public     :: cgils_ctl_s6_nc  = '/ncrc/home2/Yi-hsuan.Chen/work/research/edmf_AM4/data/CGILS/ctl_s6.nc'
character(len=200), public     :: cgils_ctl_s11_nc = '/ncrc/home2/Yi-hsuan.Chen/work/research/edmf_AM4/data/CGILS/ctl_s11.nc'
character(len=200), public     :: cgils_ctl_s12_nc = '/ncrc/home2/Yi-hsuan.Chen/work/research/edmf_AM4/data/CGILS/ctl_s12.nc'
character(len=200), public     :: cgils_p2k_s6_nc  = '/ncrc/home2/Yi-hsuan.Chen/work/research/edmf_AM4/data/CGILS/p2k_s6.nc'
character(len=200), public     :: cgils_p2k_s11_nc = '/ncrc/home2/Yi-hsuan.Chen/work/research/edmf_AM4/data/CGILS/p2k_s11.nc'
character(len=200), public     :: cgils_p2k_s12_nc = '/ncrc/home2/Yi-hsuan.Chen/work/research/edmf_AM4/data/CGILS/p2k_s12.nc'

logical :: do_debug_printout = .true.

namelist /scm_cgils_nml/ &
                        do_debug_printout, &
                        cgils_case, &
                        tracer_vert_advec_scheme, &
                        temp_vert_advec_scheme,   &
                        momentum_vert_advec_scheme

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
  character(len=200)     :: cgils_nc 
  integer, dimension(4) :: siz   

  real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: &   ! dimension(lat, lon, lev+1)
    phalf   ! pressure at half levels (Pa)

  real, dimension(size(pt,1),size(pt,2),size(pt,3)) :: &     ! dimension(lat, lon, lev)
    pfull  ! pressure at full levels (Pa)
 
  integer itime
  integer i,j,k,ix,jx,kx,tx

!---------------------------------

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
  call field_size (cgils_nc, 'time', siz) ; time_cgils = siz(1)
  
  if (allocated(buffer_cgils_3D)) deallocate (buffer_cgils_3D)
    allocate(buffer_cgils_3D(nlat_cgils, nlon_cgils, nlev_cgils)); buffer_cgils_3D = missing_value

  if (allocated(buffer_cgils_2D)) deallocate (buffer_cgils_2D)
    allocate(buffer_cgils_2D(nlat_cgils, nlon_cgils)); buffer_cgils_2D = missing_value

  if (allocated(buffer_cgils_1D_lev)) deallocate (buffer_cgils_1D_lev)
    allocate(buffer_cgils_1D_lev(nlev_cgils)); buffer_cgils_1D_lev = missing_value

  !--- set cgils variables
  !if (allocated(T_cgils)) deallocate(T_cgils)
  !  allocate(T_cgils(nlev_cgils)); T_cgils = missing_value

if (do_debug_printout) then
  write(6,*) 'nlev_cgils, time_cgils',nlev_cgils, time_cgils
  write(6,*) 'siz',siz
endif

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

  !--- read CGILD 3D variables
  itime = 1

  if (allocated(Ps_cgils)) deallocate(Ps_cgils) ; allocate(Ps_cgils(time_cgils)); Ps_cgils = missing_value
  !call read_data(cgils_nc, 'Ps',   buffer_cgils_2D(:,:), no_domain=.true., timelevel=itime)
  call read_data(cgils_nc, 'Ps',   buffer_cgils_2D(:,:), no_domain=.true.)
       Ps_cgils(:) = buffer_cgils_2D(1,1)

  !--- read CGILD 1D variables
  if (allocated(pfull_cgils)) deallocate(pfull_cgils) ; allocate(pfull_cgils(nlev_cgils)); pfull_cgils = missing_value
  call read_data(cgils_nc, 'lev',   buffer_cgils_1D_lev(:), no_domain=.true.)
       pfull_cgils(:) = buffer_cgils_1D_lev(:)

!-------------------------------------------------
! interpoate cgils data on SCM pressure levels
!   variables [ua, va, pt, q] are used by SCM; they are defined in /src/atmos_fv_dynamics/model/fv_point.inc
!-------------------------------------------------

  !--- get pfull and phalf
  call get_eta_level(nlev, Ps_cgils(1), pfull(1,1,:), phalf(1,1,:))

  !--- interpolate cgils data on SCM pressure levels
  !call interp_plev_cgils_to_SCM(pfull, pt, pfull_cgils, T_cgils)

if (do_debug_printout) then
  write(6,*) 'T_cgils',T_cgils
  write(6,*) 'U_cgils',U_cgils
  write(6,*) 'V_cgils',V_cgils
  write(6,*) 'q_cgils',q_cgils
  write(6,*) 'Ps_cgils',Ps_cgils
  write(6,*) 'pfull_cgils',pfull_cgils
  write(6,*) 'pfull',pfull
  write(6,*) 'phalf',phalf
  !write(6,*) '',
    call error_mesg( ' scm_cgils',     &
                     ' STOP.',&
                     FATAL ) 
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

end subroutine update_cgils_forc


end module scm_cgils_mod

