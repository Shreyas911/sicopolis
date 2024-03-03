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
#if defined(ALLOW_GENCTRL)
    use ad_input_m
    use ad_output_m
#endif
  
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
  
#if defined(ALLOW_GENCTRL)
    call ad_input()
#endif

    CALL SICO_INIT_D(delta_ts, glac_index, mean_accum, dtime, dtime_temp, &
  &            dtime_wss, dtime_out, dtime_ser, time, time_init, time_end&
  &            , time_output, dxi, deta, dzeta_c, dzeta_t, dzeta_r, z_mar&
  &            , ndat2d, ndat3d, n_output)
          
  !-------- Main loop --------
    CALL SICO_MAIN_LOOP_D(delta_ts, glac_index, mean_accum, dtime, &
  &                 dtime_temp, dtime_wss, dtime_out, dtime_ser, time, &
  &                 time_init, time_end, time_output, dxi, deta, dzeta_c, &
  &                 dzeta_t, dzeta_r, z_mar, ndat2d, ndat3d, n_output)
    CALL COST_FINAL_D()
       
    CALL SICO_END()

#if defined(ALLOW_GENCTRL)
     call ad_output()
#endif

    end subroutine adjoint_master
#endif

  end module tapenade_m 
  !