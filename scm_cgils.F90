module scm_cgils_mod

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

   public cgils_forc_init, cgils_forc_end, update_cgils_forc

   character(len=8) :: mod_name = 'scm_cgils'
   character(len=7) :: mod_name_diag = 'forcing'

!--------------------- version number ----------------------------------
!
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
integer, dimension(1) :: restart_versions = (/ 1 /)

!--- local variables
logical, private               :: initialized = .false.

!--- namelist parameters
integer, public                :: tracer_vert_advec_scheme = 3
integer, public                :: temp_vert_advec_scheme = 3
integer, public                :: momentum_vert_advec_scheme = 3

namelist /scm_cgils_nml/ tracer_vert_advec_scheme, &
                        temp_vert_advec_scheme,   &
                        momentum_vert_advec_scheme

!#######################################################################

contains

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
!---------------------------------

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

