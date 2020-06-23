!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Program   : m a k e _ i s m i p _ o u t p u t . F 9 0
!
!> @file
!!
!! Generating ISMIP6 output from the time-dependent NetCDF data produced by
!! SICOPOLIS.
!!
!! @section Date
!!
!! 2019-11-21
!!
!! @section Copyright
!!
!! Copyright 2010-2019 Ralf Greve, Chris Chambers
!!
!! @section License
!!
!! This file is part of SICOPOLIS.
!!
!! SICOPOLIS is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! SICOPOLIS is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with SICOPOLIS.  If not, see <http://www.gnu.org/licenses/>.
!<
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------- Settings --------

#define TIME_UNIT 1
!                     Time unit of ISMIP6 output to be generated:
!                      1 - days for time itself, seconds for other
!                          time-dependent variables (ISMIP6 default)
!                      2 - years

#define YEAR_REF 1990
!                     Reference year for the day count (for TIME_UNIT==1)

!-------- Inclusion of specification header --------

#include RUN_SPECS_HEADER

!-------------------------------------------------------------------------------
!> Declarations of kind types.
!<------------------------------------------------------------------------------
module make_ismip_output_common

integer, parameter :: i1b = selected_int_kind(2)   !< 1-byte integers
integer, parameter :: i4b = selected_int_kind(9)   ! 4-byte integers
integer, parameter :: sp  = kind(1.0)              ! single-precision reals
integer, parameter :: dp  = kind(1.0d0)            ! double-precision reals

real(dp), parameter :: no_value_large_dp = 1.0e+20_dp
real(dp), parameter :: no_value_neg_dp = -9999.0_dp
real(dp), parameter :: eps_dp = 1.0e-05_dp

end module make_ismip_output_common

!-------------------------------------------------------------------------------
!> Main program:
!! Generating ISMIP6 output from the time-dependent NetCDF data produced by
!! SICOPOLIS.
!<------------------------------------------------------------------------------
program make_ismip_output

use make_ismip_output_common

implicit none

integer(i4b)       :: i, j, n
integer(i4b)       :: n_variable_dim, n_variable_type
integer(i4b)       :: ndat, ncid
character(len=256) :: runname
character(len=  4) :: ergnum
logical            :: flag_init_output

integer(i1b)       :: mapping_val
real(dp)           :: mapping_semi_major_axis_val
real(dp)           :: mapping_inv_flattening_val
real(dp)           :: mapping_radius_of_sphere_val
real(dp)           :: mapping_latitude_origin_val
real(dp)           :: mapping_standard_parallel_val
real(dp)           :: mapping_reference_longitude_val
real(dp)           :: mapping_false_E_val
real(dp)           :: mapping_false_N_val
real(dp)           :: time_val, year_val
real(dp)           :: time_bnds_val(2)
real(dp)           :: x_val(0:IMAX)
real(dp)           :: y_val(0:JMAX)
real(dp)           :: lon_val(0:IMAX,0:JMAX)
real(dp)           :: lat_val(0:IMAX,0:JMAX)
real(dp)           :: lim_val
real(dp)           :: limnsw_val
real(dp)           :: iareagr_val
real(dp)           :: iareafl_val
real(dp)           :: tendacabf_val
real(dp)           :: tendlibmassbf_val
real(dp)           :: tendlibmassbffl_val
real(dp)           :: tendlicalvf_val
real(dp)           :: tendlifmassbf_val
real(dp)           :: tendligroundf_val
real(dp)           :: lithk_val(0:IMAX,0:JMAX)
real(dp)           :: orog_val(0:IMAX,0:JMAX)
real(dp)           :: base_val(0:IMAX,0:JMAX)
real(dp)           :: topg_val(0:IMAX,0:JMAX)
real(dp)           :: acabf_val(0:IMAX,0:JMAX)
real(dp)           :: libmassbf_val(0:IMAX,0:JMAX)
real(dp)           :: libmassbfgr_val(0:IMAX,0:JMAX)
real(dp)           :: libmassbffl_val(0:IMAX,0:JMAX)
real(dp)           :: hfgeoubed_val(0:IMAX,0:JMAX)
real(dp)           :: dlithkdt_val(0:IMAX,0:JMAX)
real(dp)           :: xvelsurf_val(0:IMAX,0:JMAX)
real(dp)           :: yvelsurf_val(0:IMAX,0:JMAX)
real(dp)           :: zvelsurf_val(0:IMAX,0:JMAX)
real(dp)           :: horvelsurf_val(0:IMAX,0:JMAX)
real(dp)           :: xvelbase_val(0:IMAX,0:JMAX)
real(dp)           :: yvelbase_val(0:IMAX,0:JMAX)
real(dp)           :: zvelbase_val(0:IMAX,0:JMAX)
real(dp)           :: horvelbase_val(0:IMAX,0:JMAX)
real(dp)           :: xvelmean_val(0:IMAX,0:JMAX)
real(dp)           :: yvelmean_val(0:IMAX,0:JMAX)
real(dp)           :: horvelmean_val(0:IMAX,0:JMAX)
real(dp)           :: litemptop_val(0:IMAX,0:JMAX)
real(dp)           :: litempbot_val(0:IMAX,0:JMAX)
real(dp)           :: litempbotgr_val(0:IMAX,0:JMAX)
real(dp)           :: litempbotfl_val(0:IMAX,0:JMAX)
real(dp)           :: strbasemag_val(0:IMAX,0:JMAX)
real(dp)           :: licalvf_val(0:IMAX,0:JMAX)
real(dp)           :: lifmassbf_val(0:IMAX,0:JMAX)
real(dp)           :: ligroundf_val(0:IMAX,0:JMAX)
real(dp)           :: sftgif_val(0:IMAX,0:JMAX)
real(dp)           :: sftgrf_val(0:IMAX,0:JMAX)
real(dp)           :: sftflf_val(0:IMAX,0:JMAX)
character(len= 64) :: mapping_grid_mapping_name_val
character(len= 64) :: mapping_ellipsoid_val

write(6,'(a)') ' '

runname = RUNNAME

#if (defined(OUTPUT_INIT))

#if (OUTPUT_INIT==0)
   flag_init_output = .false.
#elif (OUTPUT_INIT==1)
   flag_init_output = .true.
#else
   stop ' >>> make_ismip_output: OUTPUT_INIT must be either 0 or 1!'
#endif

#else
   flag_init_output = .true.   ! default for undefined parameter OUTPUT_INIT
#endif

#if (OUTPUT_FLUX_VARS != 2)
   stop ' >>> make_ismip_output: OUTPUT_FLUX_VARS==2 required!'
#endif

do n_variable_dim = 1, 2
                    ! 1 - 2d fields
                    ! 2 - scalar variables

do n_variable_type = 1, 2
                    ! 1 - state variables "ST"
                    ! 2 - flux variables "FL"

   if (n_variable_dim == 1) then

#if (OUTPUT==1 || OUTPUT==3)
      if (flag_init_output) then
         ndat = nint((TIME_END0-TIME_INIT0)/DTIME_OUT0)+1
      else
         ndat = nint((TIME_END0-TIME_INIT0)/DTIME_OUT0)
      end if
#else   /* OUTPUT==2 */
      stop ' >>> make_ismip_output: OUTPUT==2 not allowed!'
#endif

   else if (n_variable_dim == 2) then

      if (flag_init_output) then
         ndat = nint((TIME_END0-TIME_INIT0)/DTIME_SER0)+1
      else
         ndat = nint((TIME_END0-TIME_INIT0)/DTIME_SER0)
      end if

   else
      stop ' >>> make_ismip_output: n_variable_dim must be either 1 or 2!'
   end if

   if (n_variable_dim == 1) then
      if (ndat > 9999) stop ' >>> make_ismip_output: Too many records!'
   end if

   do n=1, ndat

      write(ergnum, '(i0.4)') n

      call read_nc(runname, n_variable_dim, n_variable_type, ergnum, n, &
                mapping_val, &
                mapping_grid_mapping_name_val, &
                mapping_ellipsoid_val, &
                mapping_semi_major_axis_val, &
                mapping_inv_flattening_val, &
                mapping_radius_of_sphere_val, &
                mapping_latitude_origin_val, &
                mapping_standard_parallel_val, &
                mapping_reference_longitude_val, &
                mapping_false_E_val, &
                mapping_false_N_val, &
                time_val, year_val, time_bnds_val, x_val, y_val, &
                lon_val, lat_val, &
                lim_val, limnsw_val, iareagr_val, iareafl_val, &
                tendacabf_val, &
                tendlibmassbf_val, tendlibmassbffl_val, &
                tendlicalvf_val, tendlifmassbf_val, &
                tendligroundf_val, &
                lithk_val, orog_val, base_val, topg_val, &
                acabf_val, &
                libmassbf_val, libmassbfgr_val, libmassbffl_val, &
                hfgeoubed_val, &
                dlithkdt_val, &
                xvelsurf_val, yvelsurf_val, zvelsurf_val, horvelsurf_val, &
                xvelbase_val, yvelbase_val, zvelbase_val, horvelbase_val, &
                xvelmean_val, yvelmean_val, horvelmean_val, &
                litemptop_val, &
                litempbot_val, litempbotgr_val, litempbotfl_val, &
                strbasemag_val, &
                licalvf_val, lifmassbf_val, ligroundf_val, &
                sftgif_val, sftgrf_val, sftflf_val)

      if (n==1) then   ! initialization only once
         call init_ismip_netcdf(runname, n_variable_dim, n_variable_type, &
                    mapping_grid_mapping_name_val, &
                    mapping_ellipsoid_val, &
                    mapping_semi_major_axis_val, &
                    mapping_inv_flattening_val, &
                    mapping_radius_of_sphere_val, &
                    mapping_latitude_origin_val, &
                    mapping_standard_parallel_val, &
                    mapping_reference_longitude_val, &
                    mapping_false_E_val, &
                    mapping_false_N_val, &
                    ncid)
      end if

      call write_ismip_netcdf(runname, n_variable_dim, n_variable_type, &
                    n, ncid, mapping_val, &
                    time_val, year_val, time_bnds_val, x_val, y_val, &
                    lon_val, lat_val, &
                    lim_val, limnsw_val, iareagr_val, iareafl_val, &
                    tendacabf_val, &
                    tendlibmassbf_val, tendlibmassbffl_val, &
                    tendlicalvf_val, tendlifmassbf_val, &
                    tendligroundf_val, &
                    lithk_val, orog_val, base_val, topg_val, &
                    acabf_val, &
                    libmassbf_val, libmassbfgr_val, libmassbffl_val, &
                    hfgeoubed_val, &
                    dlithkdt_val, &
                    xvelsurf_val, yvelsurf_val, zvelsurf_val, horvelsurf_val, &
                    xvelbase_val, yvelbase_val, zvelbase_val, horvelbase_val, &
                    xvelmean_val, yvelmean_val, horvelmean_val, &
                    litemptop_val, &
                    litempbot_val, litempbotgr_val, litempbotfl_val, &
                    strbasemag_val, &
                    licalvf_val, lifmassbf_val, ligroundf_val, &
                    sftgif_val, sftgrf_val, sftflf_val)

   end do

   call close_ismip_netcdf(ncid)

end do   ! n_variable_type = 1, 2

end do   ! n_variable_dim = 1, 2

end program make_ismip_output

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Reading and processing of data of the time-dependent NetCDF data.
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_nc(runname, n_variable_dim, n_variable_type, ergnum, n, &
                   mapping_r, &
                   mapping_grid_mapping_name_r, &
                   mapping_ellipsoid_r, &
                   mapping_semi_major_axis_r, &
                   mapping_inv_flattening_r, &
                   mapping_radius_of_sphere_r, &
                   mapping_latitude_origin_r, &
                   mapping_standard_parallel_r, &
                   mapping_reference_longitude_r, &
                   mapping_false_E_r, &
                   mapping_false_N_r, &
                   time_r, year_r, time_bnds_r, x_r, y_r, &
                   lon_r, lat_r, &
                   lim_r, limnsw_r, iareagr_r, iareafl_r, &
                   tendacabf_r, &
                   tendlibmassbf_r, tendlibmassbffl_r, &
                   tendlicalvf_r, tendlifmassbf_r, &
                   tendligroundf_r, &
                   lithk_r, orog_r, base_r, topg_r, &
                   acabf_r, &
                   libmassbf_r, libmassbfgr_r, libmassbffl_r, &
                   hfgeoubed_r, &
                   dlithkdt_r, &
                   xvelsurf_r, yvelsurf_r, zvelsurf_r, horvelsurf_r, &
                   xvelbase_r, yvelbase_r, zvelbase_r, horvelbase_r, &
                   xvelmean_r, yvelmean_r, horvelmean_r, &
                   litemptop_r, &
                   litempbot_r, litempbotgr_r, litempbotfl_r, &
                   strbasemag_r, &
                   licalvf_r, lifmassbf_r, ligroundf_r, &
                   sftgif_r, sftgrf_r, sftflf_r)

use make_ismip_output_common
use netcdf

implicit none

character(len=256), intent(in) :: runname
character(len=  4), intent(in) :: ergnum
integer(i4b),       intent(in) :: n_variable_dim, n_variable_type, n

integer(i1b),      intent(out) :: mapping_r
real(dp),          intent(out) :: mapping_semi_major_axis_r
real(dp),          intent(out) :: mapping_inv_flattening_r
real(dp),          intent(out) :: mapping_radius_of_sphere_r
real(dp),          intent(out) :: mapping_latitude_origin_r
real(dp),          intent(out) :: mapping_standard_parallel_r
real(dp),          intent(out) :: mapping_reference_longitude_r
real(dp),          intent(out) :: mapping_false_E_r
real(dp),          intent(out) :: mapping_false_N_r
real(dp),          intent(out) :: time_r, year_r
real(dp),          intent(out) :: time_bnds_r(2)
real(dp),          intent(out) :: x_r(0:IMAX)
real(dp),          intent(out) :: y_r(0:JMAX)
real(dp),          intent(out) :: lon_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: lat_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: lim_r
real(dp),          intent(out) :: limnsw_r
real(dp),          intent(out) :: iareagr_r
real(dp),          intent(out) :: iareafl_r
real(dp),          intent(out) :: tendacabf_r
real(dp),          intent(out) :: tendlibmassbf_r
real(dp),          intent(out) :: tendlibmassbffl_r
real(dp),          intent(out) :: tendlicalvf_r
real(dp),          intent(out) :: tendlifmassbf_r
real(dp),          intent(out) :: tendligroundf_r
real(dp),          intent(out) :: lithk_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: orog_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: base_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: topg_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: acabf_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: libmassbf_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: libmassbfgr_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: libmassbffl_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: hfgeoubed_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: dlithkdt_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: xvelsurf_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: yvelsurf_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: zvelsurf_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: horvelsurf_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: xvelbase_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: yvelbase_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: zvelbase_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: horvelbase_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: xvelmean_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: yvelmean_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: horvelmean_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: litemptop_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: litempbot_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: litempbotgr_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: litempbotfl_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: strbasemag_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: licalvf_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: lifmassbf_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: ligroundf_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: sftgif_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: sftgrf_r(0:IMAX,0:JMAX)
real(dp),          intent(out) :: sftflf_r(0:IMAX,0:JMAX)
character(len=64), intent(out) :: mapping_grid_mapping_name_r
character(len=64), intent(out) :: mapping_ellipsoid_r

integer(i4b)       :: i, j
real(dp)           :: year_to_year_or_sec
character(len=256) :: filename, filename_with_path

integer(i1b), dimension(0:IMAX,0:JMAX) :: maske_erg
real(dp) :: time_erg=0.0_dp, &
            V_tot_erg=0.0_dp, V_af_erg=0.0_dp, &
            A_grounded_erg=0.0_dp, A_floating_erg=0.0_dp, &
            Q_s_erg=0.0_dp, &
            bmb_tot_erg=0.0_dp, bmb_fl_tot_erg=0.0_dp, &
            calv_tot_erg=0.0_dp, &
            xi_erg(0:IMAX), eta_erg(0:JMAX), &
            sigma_level_c_erg(0:KCMAX), sigma_level_t_erg(0:KTMAX), &
            sigma_level_r_erg(0:KRMAX)
real(dp), dimension(0:IMAX,0:JMAX) :: lon_erg, lat_erg, &
            H_erg, zs_erg, zb_erg, zl_erg, zl0_erg, &
            as_perp_apl_erg, Q_b_apl_erg, calving_apl_erg, &
            q_geo_erg, &
            dH_dtau_erg, &
            vx_s_g_erg, vy_s_g_erg, vz_s_erg, vh_s_erg, &
            vx_b_g_erg, vy_b_g_erg, vz_b_erg, vh_b_erg, &
            vx_m_g_erg, vy_m_g_erg,           vh_m_erg, &
            temp_s_erg, temp_b_erg, temph_b_erg, &
            H_w_erg, p_b_w_erg, &
            tau_b_driving_erg, tau_b_drag_erg, &
            q_gl_g_erg

integer(i4b) :: ios
integer(i4b) :: ierr1, ierr2

integer(i4b) :: ncid, ncv, ncv_test1, ncv_test2
!     ncid:      ID of the NetCDF file
!     ncv:       Variable ID
integer(i4b) :: nc1cor(1)

real(dp) :: time_aux_1, time_aux_2, year_aux_1, year_aux_2
real(dp) :: r_aux

character(len=64) :: ch_aux

real(dp), parameter :: T0  = 273.15_dp
real(dp), parameter :: rho = 910.0_dp

real(dp), parameter :: year_to_sec = 31556926.0_dp
real(dp), parameter :: year_to_day =      360.0_dp   ! 360_day calendar used

!-------- Name of file --------

if (n_variable_dim == 1) then

#if (OUTPUT==1 && ERGDAT==0)
   filename = trim(runname)//'_2d_'//trim(ergnum)//'.nc'
#elif (OUTPUT==1 && ERGDAT==1)
   filename = trim(runname)//trim(ergnum)//'.nc'
#elif (OUTPUT==3)
   filename = trim(runname)//'_2d_'//trim(ergnum)//'.nc'
#else
   write(6,'(/a)') ' >>> read_nc: File name cannot be determined!'
   write(6, '(a)') '          (OUTPUT==2 not allowed, other settings wrong?)'
   stop
#endif

write(6,'(a)') ' Now reading '//trim(filename)//' ...'

else if (n_variable_dim == 2) then

   filename = trim(runname)//'_ser.nc'

   write(6,'(a,i0,a)') ' Now reading ' &
                        //trim(filename)//', record no. ', n, ' ...'

end if

!-------- Opening of file --------

filename_with_path = trim(OUTPATH)//'/'//trim(filename)

ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

if (ios /= nf90_noerr) then
   write(6,'(/a)') ' >>> read_nc: Error when opening the time-slice file!'
   stop
end if

!-------- Reading of data --------

mapping_r = 1_i1b                     ! initial value

mapping_grid_mapping_name_r = 'xxx'   ! initial value
mapping_ellipsoid_r         = 'xxx'   ! initial value

mapping_semi_major_axis_r     = no_value_neg_dp   ! initial value
mapping_inv_flattening_r      = no_value_neg_dp   ! initial value
mapping_radius_of_sphere_r    = no_value_neg_dp   ! initial value
mapping_latitude_origin_r     = no_value_neg_dp   ! initial value
mapping_standard_parallel_r   = no_value_neg_dp   ! initial value
mapping_reference_longitude_r = no_value_neg_dp   ! initial value
mapping_false_E_r             = no_value_neg_dp   ! initial value
mapping_false_N_r             = no_value_neg_dp   ! initial value

if (n_variable_dim == 1) then

if ( nf90_inq_varid(ncid, 'mapping', ncv_test1) == nf90_noerr ) then 
   ncv = ncv_test1
else if ( nf90_inq_varid(ncid, 'crs', ncv_test2) == nf90_noerr ) then 
   ncv = ncv_test2
end if

call check( nf90_get_var(ncid, ncv, mapping_r) )

call check( nf90_get_att(ncid, ncv, 'grid_mapping_name', &
                                     mapping_grid_mapping_name_r) )

ierr1 = nf90_get_att(ncid, ncv, 'ellipsoid', ch_aux)
if (ierr1 == nf90_noerr) mapping_ellipsoid_r = trim(ch_aux)

ierr1 = nf90_get_att(ncid, ncv, 'semi_major_axis', r_aux)

if (ierr1 == nf90_noerr) then
   mapping_semi_major_axis_r = r_aux
   call check( nf90_get_att(ncid, ncv, 'inverse_flattening', &
                                        mapping_inv_flattening_r) )
else
   ierr2 = nf90_get_att(ncid, ncv, 'radius_of_sphere', r_aux)
   if (ierr2 == nf90_noerr) mapping_radius_of_sphere_r = r_aux
end if

#if (GRID==0 || GRID==1)

call check( nf90_get_att(ncid, ncv, 'latitude_of_projection_origin', &
                                     mapping_latitude_origin_r) )
call check( nf90_get_att(ncid, ncv, 'standard_parallel', &
                                     mapping_standard_parallel_r) )
call check( nf90_get_att(ncid, ncv, 'straight_vertical_longitude_from_pole', &
                                     mapping_reference_longitude_r) )
call check( nf90_get_att(ncid, ncv, 'false_easting',  mapping_false_E_r) )
call check( nf90_get_att(ncid, ncv, 'false_northing', mapping_false_N_r) )

#endif

call check( nf90_inq_varid(ncid, 'time', ncv) )
call check( nf90_get_var(ncid, ncv, time_erg) )

call check( nf90_inq_varid(ncid, 'x', ncv) )
call check( nf90_get_var(ncid, ncv, xi_erg) )

call check( nf90_inq_varid(ncid, 'y', ncv) )
call check( nf90_get_var(ncid, ncv, eta_erg) )

call check( nf90_inq_varid(ncid, 'sigma_level_c', ncv) )
call check( nf90_get_var(ncid, ncv, sigma_level_c_erg) )

call check( nf90_inq_varid(ncid, 'sigma_level_t', ncv) )
call check( nf90_get_var(ncid, ncv, sigma_level_t_erg) )

call check( nf90_inq_varid(ncid, 'sigma_level_r', ncv) )
call check( nf90_get_var(ncid, ncv, sigma_level_r_erg) )

call check( nf90_inq_varid(ncid, 'lon', ncv) )
call check( nf90_get_var(ncid, ncv, lon_erg) )

call check( nf90_inq_varid(ncid, 'lat', ncv) )
call check( nf90_get_var(ncid, ncv, lat_erg) )

call check( nf90_inq_varid(ncid, 'H', ncv) )
call check( nf90_get_var(ncid, ncv, H_erg) )

call check( nf90_inq_varid(ncid, 'zs', ncv) )
call check( nf90_get_var(ncid, ncv, zs_erg) )

call check( nf90_inq_varid(ncid, 'zb', ncv) )
call check( nf90_get_var(ncid, ncv, zb_erg) )

call check( nf90_inq_varid(ncid, 'zl', ncv) )
call check( nf90_get_var(ncid, ncv, zl_erg) )

call check( nf90_inq_varid(ncid, 'zl0', ncv) )
call check( nf90_get_var(ncid, ncv, zl0_erg) )

call check( nf90_inq_varid(ncid, 'as_perp_apl', ncv) )
call check( nf90_get_var(ncid, ncv, as_perp_apl_erg) )

call check( nf90_inq_varid(ncid, 'Q_b_apl', ncv) )
call check( nf90_get_var(ncid, ncv, Q_b_apl_erg) )

ierr1 = nf90_inq_varid(ncid, 'calving_apl', ncv)
if (ierr1 == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, calving_apl_erg) )
else
   calving_apl_erg = no_value_large_dp
end if

call check( nf90_inq_varid(ncid, 'q_geo', ncv) )
call check( nf90_get_var(ncid, ncv, q_geo_erg) )

call check( nf90_inq_varid(ncid, 'dH_dt', ncv) )
call check( nf90_get_var(ncid, ncv, dH_dtau_erg) )

call check( nf90_inq_varid(ncid, 'vx_s_g', ncv) )
call check( nf90_get_var(ncid, ncv, vx_s_g_erg) )

call check( nf90_inq_varid(ncid, 'vy_s_g', ncv) )
call check( nf90_get_var(ncid, ncv, vy_s_g_erg) )

call check( nf90_inq_varid(ncid, 'vz_s', ncv) )
call check( nf90_get_var(ncid, ncv, vz_s_erg) )

call check( nf90_inq_varid(ncid, 'vh_s', ncv) )
call check( nf90_get_var(ncid, ncv, vh_s_erg) )

call check( nf90_inq_varid(ncid, 'vx_b_g', ncv) )
call check( nf90_get_var(ncid, ncv, vx_b_g_erg) )

call check( nf90_inq_varid(ncid, 'vy_b_g', ncv) )
call check( nf90_get_var(ncid, ncv, vy_b_g_erg) )

call check( nf90_inq_varid(ncid, 'vz_b', ncv) )
call check( nf90_get_var(ncid, ncv, vz_b_erg) )

call check( nf90_inq_varid(ncid, 'vh_b', ncv) )
call check( nf90_get_var(ncid, ncv, vh_b_erg) )

call check( nf90_inq_varid(ncid, 'vx_m_g', ncv) )
call check( nf90_get_var(ncid, ncv, vx_m_g_erg) )

call check( nf90_inq_varid(ncid, 'vy_m_g', ncv) )
call check( nf90_get_var(ncid, ncv, vy_m_g_erg) )

call check( nf90_inq_varid(ncid, 'vh_m', ncv) )
call check( nf90_get_var(ncid, ncv, vh_m_erg) )

call check( nf90_inq_varid(ncid, 'temp_s', ncv) )
call check( nf90_get_var(ncid, ncv, temp_s_erg) )

call check( nf90_inq_varid(ncid, 'temp_b', ncv) )
call check( nf90_get_var(ncid, ncv, temp_b_erg) )

call check( nf90_inq_varid(ncid, 'temph_b', ncv) )
call check( nf90_get_var(ncid, ncv, temph_b_erg) )

call check( nf90_inq_varid(ncid, 'H_w', ncv) )
call check( nf90_get_var(ncid, ncv, H_w_erg) )

call check( nf90_inq_varid(ncid, 'p_b_w', ncv) )
call check( nf90_get_var(ncid, ncv, p_b_w_erg) )

call check( nf90_inq_varid(ncid, 'tau_b_driving', ncv) )
call check( nf90_get_var(ncid, ncv, tau_b_driving_erg) )

call check( nf90_inq_varid(ncid, 'tau_b_drag', ncv) )
call check( nf90_get_var(ncid, ncv, tau_b_drag_erg) )

call check( nf90_inq_varid(ncid, 'q_gl_g', ncv) )
call check( nf90_get_var(ncid, ncv, q_gl_g_erg) )

call check( nf90_inq_varid(ncid, 'maske', ncv) )
call check( nf90_get_var(ncid, ncv, maske_erg) )

else if (n_variable_dim == 2) then

nc1cor(1) = n

call check( nf90_inq_varid(ncid, 't', ncv) )
call check( nf90_get_var(ncid, ncv, time_erg, start=nc1cor) )

call check( nf90_inq_varid(ncid, 'V_tot', ncv) )
call check( nf90_get_var(ncid, ncv, V_tot_erg, start=nc1cor) )

call check( nf90_inq_varid(ncid, 'V_af', ncv) )
call check( nf90_get_var(ncid, ncv, V_af_erg, start=nc1cor) )

call check( nf90_inq_varid(ncid, 'A_grounded', ncv) )
call check( nf90_get_var(ncid, ncv, A_grounded_erg, start=nc1cor) )

call check( nf90_inq_varid(ncid, 'A_floating', ncv) )
call check( nf90_get_var(ncid, ncv, A_floating_erg, start=nc1cor) )

call check( nf90_inq_varid(ncid, 'Q_s', ncv) )
call check( nf90_get_var(ncid, ncv, Q_s_erg, start=nc1cor) )

call check( nf90_inq_varid(ncid, 'bmb_tot', ncv) )
call check( nf90_get_var(ncid, ncv, bmb_tot_erg, start=nc1cor) )

ierr1 = nf90_inq_varid(ncid, 'bmb_fl_tot', ncv)
if (ierr1 == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, bmb_fl_tot_erg, start=nc1cor) )
else
   ierr2 = nf90_inq_varid(ncid, 'bmb_si_tot', ncv)
   if (ierr2 == nf90_noerr) then
      call check( nf90_get_var(ncid, ncv, bmb_fl_tot_erg, start=nc1cor) )
   end if
end if

call check( nf90_inq_varid(ncid, 'calv_tot', ncv) )
call check( nf90_get_var(ncid, ncv, calv_tot_erg, start=nc1cor) )

end if

!-------- Closing of file --------

call check( nf90_close(ncid) )

!-------- Conversion of read data --------

#if (TIME_UNIT==1)
   year_to_year_or_sec = year_to_sec   ! a -> s
#elif (TIME_UNIT==2)
   year_to_year_or_sec = 1.0_dp        ! a -> a
#else
   stop ' >>> read_nc: TIME_UNIT must be either 1 or 2!'
#endif

if (n_variable_type == 1) then   ! state variables

   time_aux_1 = time_erg
   time_aux_2 = time_erg

else if (n_variable_type == 2) then   ! flux variables

   if (n_variable_dim == 1) then
      time_aux_1 = time_erg-real(DTIME_OUT0,dp)
      time_aux_2 = time_erg
      time_erg   = time_erg-0.5_dp*real(DTIME_OUT0,dp)
                   ! shifted to the centre of the time interval
   else if (n_variable_dim == 2) then
      time_aux_1 = time_erg-real(DTIME_SER0,dp)
      time_aux_2 = time_erg
      time_erg   = time_erg-0.5_dp*real(DTIME_SER0,dp)
                   ! shifted to the centre of the time interval
   end if

end if

#if (!defined(OUT_TIMES) || OUT_TIMES==1)
   year_r     = time_erg   + YEAR_ZERO   ! year CE
   year_aux_1 = time_aux_1 + YEAR_ZERO   ! year CE
   year_aux_2 = time_aux_2 + YEAR_ZERO   ! year CE
#elif (OUT_TIMES==2)
   year_r     = time_erg     ! year CE
   year_aux_1 = time_aux_1   ! year CE
   year_aux_2 = time_aux_2   ! year CE
#endif

#if (TIME_UNIT==1)
   time_r         = (year_r-real(YEAR_REF,dp))     * year_to_day
                                 ! year CE -> days since reference datum
   time_bnds_r(1) = (year_aux_1-real(YEAR_REF,dp)) * year_to_day
                                 ! year CE -> days since reference datum
   time_bnds_r(2) = (year_aux_2-real(YEAR_REF,dp)) * year_to_day
                                 ! year CE -> days since reference datum
#elif (TIME_UNIT==2)
   time_r         = time_erg     ! a
   time_bnds_r(1) = time_aux_1   ! a
   time_bnds_r(2) = time_aux_2   ! a
#endif

if (n_variable_dim == 1) then

lim_r             = no_value_large_dp
limnsw_r          = no_value_large_dp
iareagr_r         = no_value_large_dp
iareafl_r         = no_value_large_dp
tendacabf_r       = no_value_large_dp
tendlibmassbf_r   = no_value_large_dp
tendlibmassbffl_r = no_value_large_dp
tendlicalvf_r     = no_value_large_dp
tendlifmassbf_r   = no_value_large_dp
tendligroundf_r   = no_value_large_dp

do i=0, IMAX
   x_r(i) = xi_erg(i)   ! m
end do

do j=0, JMAX
   y_r(j) = eta_erg(j)  ! m
end do

sftgif_r = 0.0_dp
sftgrf_r = 0.0_dp
sftflf_r = 0.0_dp

do i=0, IMAX
do j=0, JMAX

   lon_r(i,j)   = lon_erg(i,j)             ! deg E
   lat_r(i,j)   = lat_erg(i,j)             ! deg N
   lithk_r(i,j) = H_erg(i,j)               ! m
   orog_r(i,j)  = zs_erg(i,j)              ! m
   base_r(i,j)  = zb_erg(i,j)              ! m
   topg_r(i,j)  = zl_erg(i,j)              ! m

   acabf_r(i,j) = as_perp_apl_erg(i,j) * rho/year_to_year_or_sec
                                       ! m/a -> kg/(m2*a) | kg/(m2*s)

   libmassbf_r(i,j) = Q_b_apl_erg(i,j) * (-rho/year_to_year_or_sec)
                                       ! m/a -> kg/(m2*a) | kg/(m2*s),
                                       ! sign changed to negative for loss

   if (maske_erg(i,j)==0_i1b) then
      libmassbfgr_r(i,j) = Q_b_apl_erg(i,j) * (-rho/year_to_year_or_sec)
                                            ! m/a -> kg/(m2*a) | kg/(m2*s),
                                            ! sign changed to negative for loss
   else
      libmassbfgr_r(i,j) = no_value_large_dp
   end if

   if (maske_erg(i,j)==3_i1b) then
      libmassbffl_r(i,j) = Q_b_apl_erg(i,j) * (-rho/year_to_year_or_sec)
                                            ! m/a -> kg/(m2*a) | kg/(m2*s),
                                            ! sign changed to negative for loss
   else
      libmassbffl_r(i,j) = no_value_large_dp
   end if

   hfgeoubed_r(i,j)  = q_geo_erg(i,j)   ! W/m2

   dlithkdt_r(i,j)   = dH_dtau_erg(i,j) / year_to_year_or_sec   ! m/a | m/s
   xvelsurf_r(i,j)   = vx_s_g_erg(i,j)  / year_to_year_or_sec   ! m/a | m/s
   yvelsurf_r(i,j)   = vy_s_g_erg(i,j)  / year_to_year_or_sec   ! m/a | m/s
   zvelsurf_r(i,j)   = vz_s_erg(i,j)    / year_to_year_or_sec   ! m/a | m/s
   horvelsurf_r(i,j) = vh_s_erg(i,j)    / year_to_year_or_sec   ! m/a | m/s
   xvelbase_r(i,j)   = vx_b_g_erg(i,j)  / year_to_year_or_sec   ! m/a | m/s
   yvelbase_r(i,j)   = vy_b_g_erg(i,j)  / year_to_year_or_sec   ! m/a | m/s
   zvelbase_r(i,j)   = vz_b_erg(i,j)    / year_to_year_or_sec   ! m/a | m/s
   horvelbase_r(i,j) = vh_b_erg(i,j)    / year_to_year_or_sec   ! m/a | m/s
   xvelmean_r(i,j)   = vx_m_g_erg(i,j)  / year_to_year_or_sec   ! m/a | m/s
   yvelmean_r(i,j)   = vy_m_g_erg(i,j)  / year_to_year_or_sec   ! m/a | m/s
   horvelmean_r(i,j) = vh_m_erg(i,j)    / year_to_year_or_sec   ! m/a | m/s

   litemptop_r(i,j)  = temp_s_erg(i,j) + T0     ! C -> K
   litempbot_r(i,j)  = temp_b_erg(i,j) + T0     ! C -> K

   if (maske_erg(i,j)==0_i1b) then
      litempbotgr_r(i,j) = temp_b_erg(i,j) + T0     ! C -> K
   else
      litempbotgr_r(i,j) = no_value_large_dp
   end if

   if (maske_erg(i,j)==3_i1b) then
      litempbotfl_r(i,j) = temp_b_erg(i,j) + T0     ! C -> K
   else
      litempbotfl_r(i,j) = no_value_large_dp
   end if

#if (DYNAMICS==1 && MARGIN!=3)
   strbasemag_r(i,j) = tau_b_drag_erg(i,j)      ! Pa (SIA only)
#else
   strbasemag_r(i,j) = no_value_large_dp
#endif
   if (abs(calving_apl_erg(i,j)) < 0.999_dp*no_value_large_dp) then
      licalvf_r(i,j)   = calving_apl_erg(i,j) * (-rho/year_to_year_or_sec)
                                          ! m/a -> kg/(m2*a) | kg/(m2*s),
                                          ! sign changed to negative for loss
      lifmassbf_r(i,j) = calving_apl_erg(i,j) * (-rho/year_to_year_or_sec)
                                          ! m/a -> kg/(m2*a) | kg/(m2*s),
                                          ! sign changed to negative for loss
   else
      licalvf_r(i,j)   = no_value_large_dp
      lifmassbf_r(i,j) = no_value_large_dp
   end if

   ligroundf_r(i,j)  = no_value_large_dp
                       !!! q_gl_g_erg(i,j) * rho/year_to_year_or_sec
                                                ! m/a -> kg/(m2*a) | kg/(m2*s)

   if ((maske_erg(i,j)==0_i1b).or.(maske_erg(i,j)==3_i1b)) &
      sftgif_r(i,j)  = 1.0_dp   ! grounded or floating ice

   if (maske_erg(i,j)==0_i1b) &
      sftgrf_r(i,j)  = 1.0_dp   ! grounded ice

   if (maske_erg(i,j)==3_i1b) &
      sftflf_r(i,j)  = 1.0_dp   ! floating ice

end do
end do

else if (n_variable_dim == 2) then

lim_r             = V_tot_erg * rho       ! m3 -> kg
limnsw_r          = V_af_erg  * rho       ! m3 -> kg
iareagr_r         = A_grounded_erg        ! m2
iareafl_r         = A_floating_erg        ! m2
tendacabf_r       = Q_s_erg * rho/year_to_year_or_sec
                                ! m3/a -> kg/a | kg/s
tendlibmassbf_r   = bmb_tot_erg * rho/year_to_year_or_sec
                                ! m3/a -> kg/a | kg/s
tendlibmassbffl_r = bmb_fl_tot_erg * rho/year_to_year_or_sec
                                ! m3/a -> kg/a | kg/s
tendlicalvf_r     = calv_tot_erg * (-rho/year_to_year_or_sec)
                                ! m3/a -> kg/a | kg/s,
                                ! sign changed to negative for loss
tendlifmassbf_r   = calv_tot_erg * (-rho/year_to_year_or_sec)
                                ! m3/a -> kg/a | kg/s,
                                ! sign changed to negative for loss
tendligroundf_r   = no_value_large_dp
                                ! kg/a | kg/s

x_r           = no_value_large_dp
y_r           = no_value_large_dp
lon_r         = no_value_large_dp
lat_r         = no_value_large_dp
lithk_r       = no_value_large_dp
orog_r        = no_value_large_dp
base_r        = no_value_large_dp
topg_r        = no_value_large_dp
acabf_r       = no_value_large_dp
libmassbf_r   = no_value_large_dp
libmassbfgr_r = no_value_large_dp
libmassbffl_r = no_value_large_dp
hfgeoubed_r   = no_value_large_dp
dlithkdt_r    = no_value_large_dp
xvelsurf_r    = no_value_large_dp
yvelsurf_r    = no_value_large_dp
zvelsurf_r    = no_value_large_dp
horvelsurf_r  = no_value_large_dp
xvelbase_r    = no_value_large_dp
yvelbase_r    = no_value_large_dp
zvelbase_r    = no_value_large_dp
horvelbase_r  = no_value_large_dp
xvelmean_r    = no_value_large_dp
yvelmean_r    = no_value_large_dp
horvelmean_r  = no_value_large_dp
litemptop_r   = no_value_large_dp
litempbot_r   = no_value_large_dp
litempbotgr_r = no_value_large_dp
litempbotfl_r = no_value_large_dp
strbasemag_r  = no_value_large_dp
licalvf_r     = no_value_large_dp
lifmassbf_r   = no_value_large_dp
ligroundf_r   = no_value_large_dp
sftgif_r      = no_value_large_dp
sftgrf_r      = no_value_large_dp
sftflf_r      = no_value_large_dp

end if

end subroutine read_nc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Initialization of ISMIP6 NetCDF file.
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine init_ismip_netcdf(runname, n_variable_dim, n_variable_type, &
                   mapping_grid_mapping_name_val, &
                   mapping_ellipsoid_val, &
                   mapping_semi_major_axis_val, &
                   mapping_inv_flattening_val, &
                   mapping_radius_of_sphere_val, &
                   mapping_latitude_origin_val, &
                   mapping_standard_parallel_val, &
                   mapping_reference_longitude_val, &
                   mapping_false_E_val, &
                   mapping_false_N_val, &
                   ncid)

use make_ismip_output_common
use netcdf

implicit none

integer(i4b),       intent(in) :: n_variable_dim, n_variable_type
real(dp),           intent(in) :: mapping_semi_major_axis_val
real(dp),           intent(in) :: mapping_inv_flattening_val
real(dp),           intent(in) :: mapping_radius_of_sphere_val
real(dp),           intent(in) :: mapping_latitude_origin_val
real(dp),           intent(in) :: mapping_standard_parallel_val
real(dp),           intent(in) :: mapping_reference_longitude_val
real(dp),           intent(in) :: mapping_false_E_val
real(dp),           intent(in) :: mapping_false_N_val
character(len= 64), intent(in) :: mapping_grid_mapping_name_val
character(len= 64), intent(in) :: mapping_ellipsoid_val
character(len=256), intent(in) :: runname

integer(i4b),      intent(out) :: ncid   ! ID of the NetCDF file

integer(i4b) :: ios
integer(i4b) :: n
integer(i4b) :: ncv
!     ncv:       Variable ID
integer(i4b) :: ncd, nc1d, nc2d(2), nc3d(3)
!     ncd:       Dimension ID
!     nc1d:      Dimension of a 1-d array
!     nc2d:      Vector with the dimensions of a 2-d array
!     nc3d:      Vector with the dimensions of a 3-d array
integer(i4b) :: nc3flag(3), nc4flag(4)
!     nc3flag:   Vector with the 3 possible values of a flag variable
!     nc4flag:   Vector with the 4 possible values of a flag variable
integer(i4b) :: nc1cor(1), nc2cor(2), nc3cor(3)
!     nc1cor(1): Corner of a 1-d array
!     nc2cor(2): Corner of a 2-d array
!     nc3cor(3): Corner of a 3-d array
integer(i4b) :: nc1cnt(1), nc2cnt(2), nc3cnt(3)
!     nc1cnt(1): Count of a 1-d array
!     nc2cnt(2): Count of a 2-d array
!     nc3cnt(3): Count of a 3-d array
character(len=  16) :: ch_date, ch_time, ch_zone
character(len=  16) :: ch_year_ref
character(len= 256) :: filename, filename_with_path, buffer
character(len=8192) :: long_buffer
character           :: ch_time_unit

#if (TIME_UNIT==1)
ch_time_unit = 's'
#elif (TIME_UNIT==2)
ch_time_unit = 'a'
#else
stop ' >>> init_ismip_netcdf: TIME_UNIT must be either 1 or 2!'
#endif

!-------- Open NetCDF file --------

if ((n_variable_dim == 1).and.(n_variable_type == 1)) then
   filename = trim(runname)//'_ismip6_'//ch_time_unit//'_st_2d.nc'
else if ((n_variable_dim == 1).and.(n_variable_type == 2)) then
   filename = trim(runname)//'_ismip6_'//ch_time_unit//'_fl_2d.nc'
else if ((n_variable_dim == 2).and.(n_variable_type == 1)) then
   filename = trim(runname)//'_ismip6_'//ch_time_unit//'_st_scalar.nc'
else if ((n_variable_dim == 2).and.(n_variable_type == 2)) then
   filename = trim(runname)//'_ismip6_'//ch_time_unit//'_fl_scalar.nc'
end if

write(6,'(a)') ' Now creating new NetCDF file '//trim(filename)//' ...'

filename_with_path = trim(OUTPATH)//'/'//trim(filename)

ios = nf90_create(trim(filename_with_path), NF90_NOCLOBBER, ncid)

if (ios /= nf90_noerr) then
   write(6,'(/a)') ' >>> init_ismip_netcdf:'
   write(6,'(/a)') '     Error when opening the new ISMIP6 output file!'
   stop
end if

!-------- Global attributes --------

buffer = 'ISMIP6 output of simulation '//trim(runname)
call check( nf90_put_att(ncid, NF90_GLOBAL, 'title', trim(buffer)) )

buffer = 'Institute of Low Temperature Science, Hokkaido University, '// &
         'Sapporo, Japan'
call check( nf90_put_att(ncid, NF90_GLOBAL, 'institution', trim(buffer)) )

buffer = 'SICOPOLIS Version '//SICO_VERSION
call check( nf90_put_att(ncid, NF90_GLOBAL, 'source', trim(buffer)) )

call date_and_time(ch_date, ch_time, ch_zone)
buffer = ch_date(1:4)//'-'//ch_date(5:6)//'-'//ch_date(7:8)//' '// &
         ch_time(1:2)//':'//ch_time(3:4)//':'//ch_time(5:6)//' '// &
         ch_zone(1:3)//':'//ch_zone(4:5)//' - Data produced'
call check( nf90_put_att(ncid, NF90_GLOBAL, 'history', trim(buffer)) )

buffer = 'http://www.sicopolis.net/'
call check( nf90_put_att(ncid, NF90_GLOBAL, 'references', trim(buffer)) )

buffer = 'CF-1.6'
call check( nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', trim(buffer)) )

!-------- Definition of the dimensions --------

if (n_variable_dim == 1) then
   call check( nf90_def_dim(ncid, 'x', IMAX+1, ncd) )
   call check( nf90_def_dim(ncid, 'y', JMAX+1, ncd) )
end if

call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, ncd) )

if (n_variable_type == 2) &
   call check( nf90_def_dim(ncid, 'bnds', 2, ncd) )

!-------- Definition of the variables --------

!  ------ Mapping

if (n_variable_dim == 1) then

!    ---- mapping

call check( nf90_def_var(ncid, 'mapping', NF90_BYTE, ncv) )

call check( nf90_put_att(ncid, ncv, 'grid_mapping_name', &
                                     trim(mapping_grid_mapping_name_val)) )

if (trim(mapping_ellipsoid_val) /= 'xxx' ) &
   call check( nf90_put_att(ncid, ncv, 'ellipsoid', &
                                        trim(mapping_ellipsoid_val)) )

if (mapping_semi_major_axis_val > (no_value_neg_dp+eps_dp) ) &
   call check( nf90_put_att(ncid, ncv, 'semi_major_axis', &
                                        mapping_semi_major_axis_val) )

if (mapping_inv_flattening_val > (no_value_neg_dp+eps_dp) ) &
   call check( nf90_put_att(ncid, ncv, 'inverse_flattening', &
                                        mapping_inv_flattening_val) )

if (mapping_radius_of_sphere_val > (no_value_neg_dp+eps_dp) ) &
   call check( nf90_put_att(ncid, ncv, 'radius_of_sphere', &
                                        mapping_radius_of_sphere_val) )

call check( nf90_put_att(ncid, ncv, 'latitude_of_projection_origin', &
                                     mapping_latitude_origin_val) )
call check( nf90_put_att(ncid, ncv, 'standard_parallel', &
                                     mapping_standard_parallel_val) )
call check( nf90_put_att(ncid, ncv, 'straight_vertical_longitude_from_pole', &
                                     mapping_reference_longitude_val) )
call check( nf90_put_att(ncid, ncv, 'false_easting', mapping_false_E_val) )
call check( nf90_put_att(ncid, ncv, 'false_northing', mapping_false_N_val) )

end if

!  ------ Time

write(ch_year_ref, '(i0)') YEAR_REF

call check( nf90_inq_dimid(ncid, 'time', nc1d) )
call check( nf90_def_var(ncid, 'time', NF90_FLOAT, nc1d, ncv) )
if (n_variable_type == 2) &
   call check( nf90_put_att(ncid, ncv, 'bounds', 'time_bnds') )
#if (TIME_UNIT==1)
buffer = 'days since '//trim(adjustl(ch_year_ref))//'-1-1'
#elif (TIME_UNIT==2)
buffer = 'a'
#endif
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'time'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'time'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
#if (TIME_UNIT==1)
buffer = '360_day'
call check( nf90_put_att(ncid, ncv, 'calendar', trim(buffer)) )
#endif
call check( nf90_put_att(ncid, ncv, 'axis', 'T') )

!  ------ Time_bnds

if (n_variable_type == 2) then
   call check( nf90_inq_dimid(ncid, 'bnds', nc2d(1)) )
   call check( nf90_inq_dimid(ncid, 'time', nc2d(2)) )
   call check( nf90_def_var(ncid, 'time_bnds', NF90_FLOAT, nc2d, ncv) )
end if

!  ------ Time in years CE

call check( nf90_inq_dimid(ncid, 'time', nc1d) )
call check( nf90_def_var(ncid, 'year', NF90_FLOAT, nc1d, ncv) )
buffer = 'year CE'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'year'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Year'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!  ------ Projection coordinates

if (n_variable_dim == 1) then

!    ---- x

call check( nf90_inq_dimid(ncid, 'x', nc1d) )
call check( nf90_def_var(ncid, 'x', NF90_FLOAT, nc1d, ncv) )
buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'projection_x_coordinate'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'x-coordinate of the grid point i'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'axis', 'X') )

!    ---- y

call check( nf90_inq_dimid(ncid, 'y', nc1d) )
call check( nf90_def_var(ncid, 'y', NF90_FLOAT, nc1d, ncv) )
buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'projection_y_coordinate'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'y-coordinate of the grid point j'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'axis', 'Y') )

!    ---- lon

call check( nf90_inq_dimid(ncid, 'x', nc2d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc2d(2)) )
call check( nf90_def_var(ncid, 'lon', NF90_FLOAT, nc2d, ncv) )
buffer = 'degrees_E'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'longitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Geographical longitude'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!    ---- lat

call check( nf90_inq_dimid(ncid, 'x', nc2d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc2d(2)) )
call check( nf90_def_var(ncid, 'lat', NF90_FLOAT, nc2d, ncv) )
buffer = 'degrees_N'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'latitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Geographical latitude'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

end if

!  ------ Scalar output variables

if (n_variable_dim == 2) then

!    ---- State variables

if (n_variable_type == 1) then

!      -- lim

call check( nf90_inq_dimid(ncid, 'time', nc1d) )
call check( nf90_def_var(ncid, 'lim', NF90_FLOAT, nc1d, ncv) )
buffer = 'kg'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_mass'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Total ice mass'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!      -- limnsw

call check( nf90_inq_dimid(ncid, 'time', nc1d) )
call check( nf90_def_var(ncid, 'limnsw', NF90_FLOAT, nc1d, ncv) )
buffer = 'kg'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_mass_not_displacing_sea_water'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Mass above floatation'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!      -- iareagr

call check( nf90_inq_dimid(ncid, 'time', nc1d) )
call check( nf90_def_var(ncid, 'iareagr', NF90_FLOAT, nc1d, ncv) )
buffer = 'm2'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'grounded_land_ice_area'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Grounded ice area'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!      -- iareafl

call check( nf90_inq_dimid(ncid, 'time', nc1d) )
call check( nf90_def_var(ncid, 'iareafl', NF90_FLOAT, nc1d, ncv) )
buffer = 'm2'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'floating_ice_shelf_area'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Floating ice area'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!    ---- Flux variables

else if (n_variable_type == 2) then

!      -- tendacabf

call check( nf90_inq_dimid(ncid, 'time', nc1d) )
call check( nf90_def_var(ncid, 'tendacabf', NF90_FLOAT, nc1d, ncv) )
buffer = 'kg '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'tendency_of_land_ice_mass_due_to_surface_mass_balance'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Total SMB flux'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!      -- tendlibmassbf

call check( nf90_inq_dimid(ncid, 'time', nc1d) )
call check( nf90_def_var(ncid, 'tendlibmassbf', NF90_FLOAT, nc1d, ncv) )
buffer = 'kg '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'tendency_of_land_ice_mass_due_to_basal_mass_balance'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Total BMB flux'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!      -- tendlibmassbffl

call check( nf90_inq_dimid(ncid, 'time', nc1d) )
call check( nf90_def_var(ncid, 'tendlibmassbffl', NF90_FLOAT, nc1d, ncv) )
buffer = 'kg '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'tendency_of_land_ice_mass_due_to_basal_mass_balance_for_floating_ice'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Total BMB flux for floating ice'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!      -- tendlicalvf

call check( nf90_inq_dimid(ncid, 'time', nc1d) )
call check( nf90_def_var(ncid, 'tendlicalvf', NF90_FLOAT, nc1d, ncv) )
buffer = 'kg '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'tendency_of_land_ice_mass_due_to_calving'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Total calving flux'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!      -- tendlifmassbf

call check( nf90_inq_dimid(ncid, 'time', nc1d) )
call check( nf90_def_var(ncid, 'tendlifmassbf', NF90_FLOAT, nc1d, ncv) )
buffer = 'kg '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'tendency_of_land_ice_mass_due_to_calving_and_ice_front_melting'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Total calving and ice front melting flux'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!!! !      -- tendligroundf
!!! 
!!! call check( nf90_inq_dimid(ncid, 'time', nc1d) )
!!! call check( nf90_def_var(ncid, 'tendligroundf', NF90_FLOAT, nc1d, ncv) )
!!! buffer = 'kg '//ch_time_unit//'-1'
!!! call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
!!! buffer = 'tendency_of_land_ice_mass_due_to_flux_at_grounding_line'
!!! call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
!!! buffer = 'Total grounding line flux'
!!! call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

end if

end if

!  ------ 2D output variables

if (n_variable_dim == 1) then

!    ---- State variables

if (n_variable_type == 1) then

!      -- lithk

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'lithk', NF90_FLOAT, nc3d, ncv) )
buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Ice thickness'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- orog

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'orog', NF90_FLOAT, nc3d, ncv) )
buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'surface_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Surface elevation'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- base

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'base', NF90_FLOAT, nc3d, ncv) )
buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'base_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Ice base elevation'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- topg

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'topg', NF90_FLOAT, nc3d, ncv) )
buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'bedrock_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Bedrock elevation'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- xvelsurf

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'xvelsurf', NF90_FLOAT, nc3d, ncv) )
buffer = 'm '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_surface_x_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Surface velocity in x'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- yvelsurf

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'yvelsurf', NF90_FLOAT, nc3d, ncv) )
buffer = 'm '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_surface_y_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Surface velocity in y'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- zvelsurf

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'zvelsurf', NF90_FLOAT, nc3d, ncv) )
buffer = 'm '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_surface_upward_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Surface velocity in z'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- horvelsurf

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'horvelsurf', NF90_FLOAT, nc3d, ncv) )
buffer = 'm '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_surface_horizontal_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal surface velocity'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- xvelbase

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'xvelbase', NF90_FLOAT, nc3d, ncv) )
buffer = 'm '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_basal_x_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Basal velocity in x'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- yvelbase

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'yvelbase', NF90_FLOAT, nc3d, ncv) )
buffer = 'm '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_basal_y_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Basal velocity in y'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- zvelbase

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'zvelbase', NF90_FLOAT, nc3d, ncv) )
buffer = 'm '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_basal_upward_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Basal velocity in z'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- horvelbase

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'horvelbase', NF90_FLOAT, nc3d, ncv) )
buffer = 'm '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_basal_horizontal_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal basal velocity'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- xvelmean

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'xvelmean', NF90_FLOAT, nc3d, ncv) )
buffer = 'm '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_vertical_mean_x_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Mean velocity in x'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- yvelmean

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'yvelmean', NF90_FLOAT, nc3d, ncv) )
buffer = 'm '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_vertical_mean_y_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Mean velocity in y'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- horvelmean

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'horvelmean', NF90_FLOAT, nc3d, ncv) )
buffer = 'm '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_vertical_mean_horizontal_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal mean velocity'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- litemptop

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'litemptop', NF90_FLOAT, nc3d, ncv) )
buffer = 'K'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'temperature_at_top_of_ice_sheet_model'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Surface temperature'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- litempbot

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'litempbot', NF90_FLOAT, nc3d, ncv) )
buffer = 'K'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'temperature_at_base_of_ice_sheet_model'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Basal temperature'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- litempbotgr

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'litempbotgr', NF90_FLOAT, nc3d, ncv) )
buffer = 'K'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'temperature_at_base_of_ice_sheet_model_for_grounded_ice'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Basal temperature for grounded ice'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- litempbotfl

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'litempbotfl', NF90_FLOAT, nc3d, ncv) )
buffer = 'K'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'temperature_at_base_of_ice_sheet_model_for_floating_ice'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Basal temperature for floating ice'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- strbasemag

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'strbasemag', NF90_FLOAT, nc3d, ncv) )
buffer = 'Pa'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_basal_drag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Basal drag'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- sftgif

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'sftgif', NF90_FLOAT, nc3d, ncv) )
buffer = '1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_area_fraction'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Land ice area fraction'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- sftgrf

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'sftgrf', NF90_FLOAT, nc3d, ncv) )
buffer = '1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'grounded_ice_sheet_area_fraction'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Grounded ice sheet area fraction'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- sftflf

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'sftflf', NF90_FLOAT, nc3d, ncv) )
buffer = '1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'floating_ice_shelf_area_fraction'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Floating ice sheet area fraction'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- Flux variables

else if (n_variable_type == 2) then

!      -- acabf

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'acabf', NF90_FLOAT, nc3d, ncv) )
buffer = 'kg m-2 '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_surface_specific_mass_balance_flux'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Surface mass balance flux'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- libmassbf

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'libmassbf', NF90_FLOAT, nc3d, ncv) )
buffer = 'kg m-2 '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_basal_specific_mass_balance_flux'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Basal mass balance flux'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- libmassbfgr

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'libmassbfgr', NF90_FLOAT, nc3d, ncv) )
buffer = 'kg m-2 '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_basal_specific_mass_balance_flux_for_grounded_ice'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Basal mass balance flux for grounded ice'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- libmassbffl

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'libmassbffl', NF90_FLOAT, nc3d, ncv) )
buffer = 'kg m-2 '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_basal_specific_mass_balance_flux_for_floating_ice'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Basal mass balance flux for floating ice'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- hfgeoubed

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'hfgeoubed', NF90_FLOAT, nc3d, ncv) )
buffer = 'W m-2'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'upward_geothermal_heat_flux_in_land_ice'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Geothermal heat flux'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- dlithkdt

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'dlithkdt', NF90_FLOAT, nc3d, ncv) )
buffer = 'm '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'tendency_of_land_ice_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Ice thickness imbalance'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- licalvf

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'licalvf', NF90_FLOAT, nc3d, ncv) )
buffer = 'kg m-2 '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_specific_mass_flux_due_to_calving'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Calving flux'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!      -- lifmassbf

call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
call check( nf90_def_var(ncid, 'lifmassbf', NF90_FLOAT, nc3d, ncv) )
buffer = 'kg m-2 '//ch_time_unit//'-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_specific_mass_flux_due_to_calving_and_ice_front_melting'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Ice front melt and calving flux'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, '_FillValue', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'missing_value', &
                                    real(no_value_large_dp,sp)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!!! !      -- ligroundf
!!! 
!!! call check( nf90_inq_dimid(ncid, 'x', nc3d(1)) )
!!! call check( nf90_inq_dimid(ncid, 'y', nc3d(2)) )
!!! call check( nf90_inq_dimid(ncid, 'time', nc3d(3)) )
!!! call check( nf90_def_var(ncid, 'ligroundf', NF90_FLOAT, nc3d, ncv) )
!!! buffer = 'kg m-2 '//ch_time_unit//'-1'
!!! call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
!!! buffer = 'land_ice_specific_mass_flux_at_grounding_line'
!!! call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
!!! buffer = 'Grounding line flux'
!!! call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
!!! call check( nf90_put_att(ncid, ncv, '_FillValue', &
!!!                                     real(no_value_large_dp,sp)) )
!!! call check( nf90_put_att(ncid, ncv, 'missing_value', &
!!!                                     real(no_value_large_dp,sp)) )
!!! call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

end if

end if

!    ---- End of the definitions

call check( nf90_enddef(ncid) )

end subroutine init_ismip_netcdf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Writing of ISMIP6 NetCDF file.
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine write_ismip_netcdf(runname, n_variable_dim, n_variable_type, &
                    n, ncid, mapping_val, &
                    time_val, year_val, time_bnds_val, x_val, y_val, &
                    lon_val, lat_val, &
                    lim_val, limnsw_val, iareagr_val, iareafl_val, &
                    tendacabf_val, &
                    tendlibmassbf_val, tendlibmassbffl_val, &
                    tendlicalvf_val, tendlifmassbf_val, &
                    tendligroundf_val, &
                    lithk_val, orog_val, base_val, topg_val, &
                    acabf_val, &
                    libmassbf_val, libmassbfgr_val, libmassbffl_val, &
                    hfgeoubed_val, &
                    dlithkdt_val, &
                    xvelsurf_val, yvelsurf_val, zvelsurf_val, horvelsurf_val, &
                    xvelbase_val, yvelbase_val, zvelbase_val, horvelbase_val, &
                    xvelmean_val, yvelmean_val, horvelmean_val, &
                    litemptop_val, &
                    litempbot_val, litempbotgr_val, litempbotfl_val, &
                    strbasemag_val, &
                    licalvf_val, lifmassbf_val, ligroundf_val, &
                    sftgif_val, sftgrf_val, sftflf_val)

use make_ismip_output_common
use netcdf

implicit none

character(len=256), intent(in) :: runname
integer(i4b),       intent(in) :: n_variable_dim, n_variable_type, n
integer(i4b),       intent(in) :: ncid   ! ID of the NetCDF file
integer(i1b),       intent(in) :: mapping_val(1)
real(dp),           intent(in) :: time_val(1), year_val(1)
real(dp),           intent(in) :: time_bnds_val(2)
real(dp),           intent(in) :: x_val(0:IMAX)
real(dp),           intent(in) :: y_val(0:JMAX)
real(dp),           intent(in) :: lon_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: lat_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: lim_val(1)
real(dp),           intent(in) :: limnsw_val(1)
real(dp),           intent(in) :: iareagr_val(1)
real(dp),           intent(in) :: iareafl_val(1)
real(dp),           intent(in) :: tendacabf_val(1)
real(dp),           intent(in) :: tendlibmassbf_val(1)
real(dp),           intent(in) :: tendlibmassbffl_val(1)
real(dp),           intent(in) :: tendlicalvf_val(1)
real(dp),           intent(in) :: tendlifmassbf_val(1)
real(dp),           intent(in) :: tendligroundf_val(1)
real(dp),           intent(in) :: lithk_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: orog_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: base_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: topg_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: acabf_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: libmassbf_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: libmassbfgr_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: libmassbffl_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: hfgeoubed_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: dlithkdt_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: xvelsurf_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: yvelsurf_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: zvelsurf_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: horvelsurf_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: xvelbase_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: yvelbase_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: zvelbase_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: horvelbase_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: xvelmean_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: yvelmean_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: horvelmean_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: litemptop_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: litempbot_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: litempbotgr_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: litempbotfl_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: strbasemag_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: licalvf_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: lifmassbf_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: ligroundf_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: sftgif_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: sftgrf_val(0:IMAX,0:JMAX)
real(dp),           intent(in) :: sftflf_val(0:IMAX,0:JMAX)

integer(i4b) :: i, j
integer(i4b) :: ncv
!     ncv:       Variable ID
integer(i4b) :: ncd, nc1d, nc2d(2), nc3d(3)
!     ncd:       Dimension ID
!     nc1d:      Dimension of a 1-d array
!     nc2d:      Vector with the dimensions of a 2-d array
!     nc3d:      Vector with the dimensions of a 3-d array
integer(i4b) :: nc3flag(3), nc4flag(4)
!     nc3flag:   Vector with the 3 possible values of a flag variable
!     nc4flag:   Vector with the 4 possible values of a flag variable
integer(i4b) :: nc1cor(1), nc2cor(2), nc3cor(3)
!     nc1cor(1): Corner of a 1-d array
!     nc2cor(2): Corner of a 2-d array
!     nc3cor(3): Corner of a 3-d array
integer(i4b) :: nc1cnt(1), nc2cnt(2), nc3cnt(3)
!     nc1cnt(1): Count of a 1-d array
!     nc2cnt(2): Count of a 2-d array
!     nc3cnt(3): Count of a 3-d array
character(len= 16) :: ch_date, ch_time, ch_zone
character(len=256) :: filename, buffer
character          :: ch_empty

!-------- Writing of data on NetCDF file --------

write(6,'(a)') ' Now writing data in new NetCDF file ...'

!  ------ Mapping variable and projection coordinates

if (n_variable_dim == 1) then

if (n==1) then   ! output only once

   call check( nf90_inq_varid(ncid, 'mapping', ncv) )
   nc1cor(1) = 1
   nc1cnt(1) = 1
   call check( nf90_put_var(ncid, ncv, mapping_val, &
                            start=nc1cor, count=nc1cnt) )

   call check( nf90_inq_varid(ncid, 'x', ncv) )
   nc1cor(1) = 1
   nc1cnt(1) = IMAX + 1
   call check( nf90_put_var(ncid, ncv, x_val, start=nc1cor, count=nc1cnt) )

   call check( nf90_inq_varid(ncid, 'y', ncv) )
   nc1cor(1) = 1
   nc1cnt(1) = JMAX + 1
   call check( nf90_put_var(ncid, ncv, y_val, start=nc1cor, count=nc1cnt) )

   call check( nf90_inq_varid(ncid, 'lon', ncv) )
   nc2cor(1) = 1
   nc2cor(2) = 1
   nc2cnt(1) = IMAX + 1
   nc2cnt(2) = JMAX + 1
   call check( nf90_put_var(ncid, ncv, lon_val, start=nc2cor, count=nc2cnt) )

   call check( nf90_inq_varid(ncid, 'lat', ncv) )
   nc2cor(1) = 1
   nc2cor(2) = 1
   nc2cnt(1) = IMAX + 1
   nc2cnt(2) = JMAX + 1
   call check( nf90_put_var(ncid, ncv, lat_val, start=nc2cor, count=nc2cnt) )

end if

end if

!  ------ Time

call check( nf90_inq_varid(ncid, 'time', ncv) )
nc1cor(1) = n
nc1cnt(1) = 1
call check( nf90_put_var(ncid, ncv, time_val, start=nc1cor, count=nc1cnt) )

!  ------ Time_bnds

if (n_variable_type == 2) then
   call check( nf90_inq_varid(ncid, 'time_bnds', ncv) )
   nc2cor(1) = 1
   nc2cor(2) = n
   nc2cnt(1) = 2
   nc2cnt(2) = 1
   call check( nf90_put_var(ncid, ncv, time_bnds_val, start=nc2cor, count=nc2cnt) )
end if

!  ------ Time in years CE

call check( nf90_inq_varid(ncid, 'year', ncv) )
nc1cor(1) = n
nc1cnt(1) = 1
call check( nf90_put_var(ncid, ncv, year_val, start=nc1cor, count=nc1cnt) )

!  ------ Scalar output variables

if (n_variable_dim == 2) then

!    ---- State variables

if (n_variable_type == 1) then

call check( nf90_inq_varid(ncid, 'lim', ncv) )
nc1cor(1) = n
nc1cnt(1) = 1
call check( nf90_put_var(ncid, ncv, lim_val, start=nc1cor, count=nc1cnt) )

call check( nf90_inq_varid(ncid, 'limnsw', ncv) )
nc1cor(1) = n
nc1cnt(1) = 1
call check( nf90_put_var(ncid, ncv, limnsw_val, start=nc1cor, count=nc1cnt) )

call check( nf90_inq_varid(ncid, 'iareagr', ncv) )
nc1cor(1) = n
nc1cnt(1) = 1
call check( nf90_put_var(ncid, ncv, iareagr_val, start=nc1cor, count=nc1cnt) )

call check( nf90_inq_varid(ncid, 'iareafl', ncv) )
nc1cor(1) = n
nc1cnt(1) = 1
call check( nf90_put_var(ncid, ncv, iareafl_val, start=nc1cor, count=nc1cnt) )

!    ---- Flux variables

else if (n_variable_type == 2) then

call check( nf90_inq_varid(ncid, 'tendacabf', ncv) )
nc1cor(1) = n
nc1cnt(1) = 1
call check( nf90_put_var(ncid, ncv, tendacabf_val, start=nc1cor, count=nc1cnt) )

call check( nf90_inq_varid(ncid, 'tendlibmassbf', ncv) )
nc1cor(1) = n
nc1cnt(1) = 1
call check( nf90_put_var(ncid, ncv, tendlibmassbf_val, start=nc1cor, count=nc1cnt) )

call check( nf90_inq_varid(ncid, 'tendlibmassbffl', ncv) )
nc1cor(1) = n
nc1cnt(1) = 1
call check( nf90_put_var(ncid, ncv, tendlibmassbffl_val, start=nc1cor, count=nc1cnt) )

call check( nf90_inq_varid(ncid, 'tendlicalvf', ncv) )
nc1cor(1) = n
nc1cnt(1) = 1
call check( nf90_put_var(ncid, ncv, tendlicalvf_val, start=nc1cor, count=nc1cnt) )

call check( nf90_inq_varid(ncid, 'tendlifmassbf', ncv) )
nc1cor(1) = n
nc1cnt(1) = 1
call check( nf90_put_var(ncid, ncv, tendlifmassbf_val, start=nc1cor, count=nc1cnt) )

!!! call check( nf90_inq_varid(ncid, 'tendligroundf', ncv) )
!!! nc1cor(1) = n
!!! nc1cnt(1) = 1
!!! call check( nf90_put_var(ncid, ncv, tendligroundf_val, start=nc1cor, count=nc1cnt) )

end if

end if

!  ------ 2D output variables

if (n_variable_dim == 1) then

!    ---- State variables

if (n_variable_type == 1) then

call check( nf90_inq_varid(ncid, 'lithk', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, lithk_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'orog', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, orog_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'base', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, base_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'topg', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, topg_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'xvelsurf', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, xvelsurf_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'yvelsurf', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, yvelsurf_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'zvelsurf', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, zvelsurf_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'horvelsurf', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, horvelsurf_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'xvelbase', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, xvelbase_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'yvelbase', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, yvelbase_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'zvelbase', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, zvelbase_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'horvelbase', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, horvelbase_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'xvelmean', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, xvelmean_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'yvelmean', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, yvelmean_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'horvelmean', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, horvelmean_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'litemptop', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, litemptop_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'litempbot', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, litempbot_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'litempbotgr', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, litempbotgr_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'litempbotfl', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, litempbotfl_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'strbasemag', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, strbasemag_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'sftgif', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, sftgif_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'sftgrf', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, sftgrf_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'sftflf', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, sftflf_val, start=nc3cor, count=nc3cnt) )

!    ---- Flux variables

else if (n_variable_type == 2) then

call check( nf90_inq_varid(ncid, 'acabf', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, acabf_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'libmassbf', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, libmassbf_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'libmassbfgr', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, libmassbfgr_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'libmassbffl', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, libmassbffl_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'hfgeoubed', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, hfgeoubed_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'dlithkdt', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, dlithkdt_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'licalvf', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, licalvf_val, start=nc3cor, count=nc3cnt) )

call check( nf90_inq_varid(ncid, 'lifmassbf', ncv) )
nc3cor(1) = 1
nc3cor(2) = 1
nc3cor(3) = n
nc3cnt(1) = IMAX + 1
nc3cnt(2) = JMAX + 1
nc3cnt(3) = 1
call check( nf90_put_var(ncid, ncv, lifmassbf_val, start=nc3cor, count=nc3cnt) )

!!! call check( nf90_inq_varid(ncid, 'ligroundf', ncv) )
!!! nc3cor(1) = 1
!!! nc3cor(2) = 1
!!! nc3cor(3) = n
!!! nc3cnt(1) = IMAX + 1
!!! nc3cnt(2) = JMAX + 1
!!! nc3cnt(3) = 1
!!! call check( nf90_put_var(ncid, ncv, ligroundf_val, start=nc3cor, count=nc3cnt) )

end if

end if

end subroutine write_ismip_netcdf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Closing of ISMIP6 NetCDF file.
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine close_ismip_netcdf(ncid)

use make_ismip_output_common
use netcdf

implicit none

integer(i4b), intent(in) :: ncid   ! ID of the NetCDF file

write(6,'(a)') ' Now closing new NetCDF file ...'

call check( nf90_sync(ncid) )

call check( nf90_close(ncid) )

write(6,'(a/)') ' Done!'

end subroutine close_ismip_netcdf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> NetCDF error capturing.
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine check(status)

use make_ismip_output_common
use netcdf

implicit none

integer(i4b), intent (in) :: status

if (status /= nf90_noerr) then 
   write(6,'(1x,a)') trim(nf90_strerror(status))
   stop ' >>> check: Stopped due to NetCDF error!'
end if

end subroutine check  

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                        End of make_ismip_output.F90
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
