!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : s i c o _ m a i n _ l o o p _ w r a p p e r _ m 
!
!> @file
!!
!! A catch-all module for openad-related subroutines.
!!
!! @section Copyright
!!
!! Copyright 2017-2020 Liz Curry-Logan, Sri Hari Krishna Narayanan
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

module sico_main_loop_wrapper_m

  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Calls sico_main_loop_iter_m. Necessarily decoupled from time loop due to use
!! of revolve. 
!<------------------------------------------------------------------------------
  subroutine sico_main_loop_wrapper(delta_ts, glac_index, &
                      mean_accum, &
                      dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
                      time, time_init, time_end, time_output, &
                      dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                      z_sl, dzsl_dtau, z_mar, &
                      ndat2d, ndat3d, n_output, &
                      runname, &
                      itercount_max,iter_temp,iter_wss,iter_ser,&
                      iter_out,iter_output)

!$openad xxx template ad_revolve_template.F
  use sico_types_m
  use sico_variables_m
  use sico_vars_m
  use boundary_m
#if (defined(GRL) && DISC>0)
  use discharge_workers_m
#endif
#if (CALCMOD==0 || CALCMOD==1 || CALCMOD==-1)
  use calc_temp_m
#elif (CALCMOD==2 || CALCMOD==3)
  use calc_temp_enth_m
#endif

  use calc_enhance_m
  use flag_update_gf_gl_cf_m
  use calc_vxy_m
  use calc_vz_m
  use calc_dxyz_m
  use calc_gia_m
  use calc_thk_m
  use calc_temp_melt_bas_m
  use calc_bas_melt_m
  use calc_thk_water_bas_m
  use sico_main_loop_iter_m
  use ctrl_m, only: myfloor

  implicit none

  integer(i4b),       intent(inout) :: n_output
  real(dp),           intent(inout) :: mean_accum
  real(dp),           intent(inout) :: dtime, dtime_temp, dtime_wss, &
                                       dtime_out, dtime_ser
  real(dp),           intent(inout) :: time_init, time_end, time_output(100)
  real(dp),           intent(inout) :: dxi, deta, dzeta_c, dzeta_t, dzeta_r
  character(len=100), intent(inout) :: runname

  integer(i4b),       intent(inout) :: ndat2d, ndat3d
  real(dp),           intent(inout) :: delta_ts, glac_index
  real(dp),           intent(inout) :: time
  real(dp),           intent(inout) :: z_sl, dzsl_dtau, z_mar

  integer(i4b)                      :: i, j, kc, kt, kr, n
  integer(i4b)                      :: itercount
  integer(i4b),       intent(inout) :: itercount_max, iter_temp, iter_wss 
  integer(i4b),       intent(inout) :: iter_ser, iter_out, iter_output(100)
  real(dp)                          :: dtime_temp_inv
  logical                           :: flag_3d_output

  real(dp)                          :: tmp, valmod
  integer(i4b)                      :: valmodint

  main_loop : do itercount=1, itercount_max
    call sico_main_loop_iter(delta_ts, glac_index, &
                      mean_accum, &
                      dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
                      time, time_init, time_end, time_output, &
                      dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                      z_sl, dzsl_dtau, z_mar, &
                      ndat2d, ndat3d, n_output, &
                      runname, &
                      itercount,iter_temp,iter_wss,iter_ser,&
                      iter_out,iter_output)
  end do main_loop   ! End of main loop

  end subroutine sico_main_loop_wrapper

end module sico_main_loop_wrapper_m
