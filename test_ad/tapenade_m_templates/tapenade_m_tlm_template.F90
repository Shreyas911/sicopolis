!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  t a p e n a d e _ m
!
!> @file
!!
!! A catch-all module for tapenade-related subroutines. 
!!
!! @section Copyright
!!
!! Copyright 2017-2022 Shreyas Sunil Gaikwad,
!!                     Liz Curry-Logan, Sri Hari Krishna Narayanan,
!!                     Patrick Heimbach, Ralf Greve
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

!-------------------------------------------------------------------------------
!> Module for all tapenade-related subroutines 
!<------------------------------------------------------------------------------
module tapenade_m

  implicit none

  private
#ifdef ALLOW_TAPENADE
  public :: adjoint_master
#endif

contains

!-------------------------------------------------------------------------------
!> Adjoint master is the main tool by which sicopolis.F90 invokes the
!! adjoint/tlm code. Its job is to figure out what mode of the adjoint code is
!! being invoked and run the appropriate subroutine. 
!<------------------------------------------------------------------------------
#ifdef ALLOW_TAPENADE
  subroutine adjoint_master

use sico_variables_m_diff
#if (defined(GRL) && DISC>0)
  use discharge_workers_m_diff
#endif
  use ice_material_properties_m_diff
  use enth_temp_omega_m_diff
  use sico_init_m_diff
  USE COST_M_DIFF
  USE SICO_TYPES_M
  USE SICO_VARS_M
  USE SICO_MAIN_LOOP_M_DIFF
  USE SICO_END_M_DIFF

  implicit none
  integer(i4b)                               :: ndat2d, ndat3d
  integer(i4b)                               :: n_output
  real(dp)                                   :: delta_ts, glac_index
  real(dp)                                   :: mean_accum
  real(dp)                                   :: dtime, dtime_temp, &
                                                dtime_wss, dtime_out, dtime_ser
  real(dp)                                   :: time, time_init, time_end
  real(dp), dimension(100)                   :: time_output
  real(dp)                                   :: dxi, deta, dzeta_c, &
                                                dzeta_t, dzeta_r
  real(dp)                                   :: z_mar
  integer(i4b), parameter                    :: points = 5
  integer(i4b), dimension(points)            :: ipoints, jpoints
  integer(i4b)                               :: i, j, p
   !-------- Test points along spines of the ice sheets
   do p = 1, points
#if (defined(GRL))
      ipoints(p) = int(real(IMAX/2))
      jpoints(p) = int(real(JMAX/5)) + (p-1) * points
#elif (defined(ANT))
      ipoints(p) = int(real(IMAX/3)) + int(real((.85-.33)*IMAX/points)) * (p - 1)
      jpoints(p) = int(real(JMAX/2))
#endif
   end do

!@ python_automated_tlm IO begin @

   !-------- Loop over points
   do p = 1, points !@ python_automated_tlm limited_or_block_or_full @
     i = ipoints(p)
     j = jpoints(p)

!@ python_automated_tlm dep_vard set 0 @

  CALL SICO_INIT_D(delta_ts, glac_index, mean_accum, dtime, dtime_temp, &
&            dtime_wss, dtime_out, dtime_ser, time, time_init, time_end&
&            , time_output, dxi, deta, dzeta_c, dzeta_t, dzeta_r, z_mar&
&            , ndat2d, ndat3d, n_output)

!@ python_automated_tlm dep_vard set 1 @
		
!-------- Main loop --------
  CALL SICO_MAIN_LOOP_D(delta_ts, glac_index, mean_accum, dtime, &
&                 dtime_temp, dtime_wss, dtime_out, dtime_ser, time, &
&                 time_init, time_end, time_output, dxi, deta, dzeta_c, &
&                 dzeta_t, dzeta_r, z_mar, ndat2d, ndat3d, n_output)
  CALL COST_FINAL_D()
     
  CALL SICO_END()

  ! Initialize compatible fields to 0
  ! 2D fields
  q_geo        = 0.0
  c_slide_init = 0.0
  H            = 0.0 ! Only compatible with ANF_DAT==1

  ! 3D fields
  temp_c       = 0.0 ! Not compatible with TEMP_INIT==5

  !@ python_automated_tlm IO write @
   end do ! (close loop over points)

!@ python_automated_tlm IO end @
  end subroutine adjoint_master
#endif

end module tapenade_m 
!