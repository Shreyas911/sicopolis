!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Program :  t a p e n a d e _ m a i n

!! Main program of SICOPOLIS_Tapenade.
!!
!! Copyright 2009-2024 SICOPOLIS Authors
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
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with SICOPOLIS. If not, see <https://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------------------------------------------------------------
!> Main program of SICOPOLIS_Tapenade.
!-------------------------------------------------------------------------------
program tapenade_main

#if defined(ALLOW_TAPENADE)

    use sico_variables_m_diff
#if (defined(GRL) && DISC>0)
    use discharge_workers_m_diff
#endif /* GRL and DISC */
    use ice_material_properties_m_diff
    use enth_temp_omega_m_diff
    use sico_init_m_diff

#if defined(ALLOW_GENCTRL)
    use ad_input_m
    use ad_output_m
#endif /* ALLOW_GENCTRL */

#if defined(ALLOW_TAP_TLM)
    use cost_m_diff
    use sico_main_loop_m_diff
    use sico_end_m_diff
#endif /* ALLOW_TAP_TLM */

#endif /* ALLOW_TAPENADE */

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

#if (defined(ALLOW_TAPENADE) && defined(ALLOW_TAP_TLM))
    integer(i4b), parameter                    :: points = 5
    integer(i4b), dimension(points)            :: ipoints, jpoints
    integer(i4b)                               :: i, j, p
#endif /* ALLOW_TAPENADE and ALLOW_TAP_TLM */

#if defined(ALLOW_TAPENADE)

#if defined(ALLOW_GENCTRL)
    call ad_input()
#endif /* ALLOW_GENCTRL */

#if defined(ALLOW_TAP_ADJ)
    fcb = 1.
    call SICOPOLIS_TAPENADE_B(delta_ts, glac_index, mean_accum, dtime, &
    & dtime_temp, dtime_wss, dtime_out, dtime_ser, time, time_init, time_end&
    & , time_output, dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
    & z_mar, ndat2d, ndat3d, n_output)

#elif (defined(ALLOW_TAP_TLM) && !defined(ALLOW_GENCTRL))

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
        q_geo          = 0.0
        c_slide_init   = 0.0
        H              = 0.0 ! Only compatible with ANF_DAT==1
#if (ACCSURFACE==2 || ACCSURFACE==3)
        gamma_s_arr    = 0.0
#endif
#if (ABLSURFACE==1 || ABLSURFACE==2 || (ACCSURFACE<=5 && SOLID_PRECIP==3))
        s_stat_arr     = 0.0
#endif
#if (ABLSURFACE==1 || ABLSURFACE==2)
        beta1_arr_orig = 0.0
        beta2_arr_orig = 0.0
        Pmax_arr       = 0.0
        mu_arr_orig    = 0.0
#endif

        ! 3D fields
        temp_c       = 0.0 ! Not compatible with TEMP_INIT==5

        ! Reset flag_ad_sico_init for next iteration
        flag_ad_sico_init = .false.

!@ python_automated_tlm IO write @
    end do ! (close loop over points)

!@ python_automated_tlm IO end @

#elif (defined(ALLOW_TAP_TLM) && defined(ALLOW_GENCTRL))

    CALL SICO_INIT_D(delta_ts, glac_index, mean_accum, dtime, dtime_temp, &
    &            dtime_wss, dtime_out, dtime_ser, time, time_init, time_end&
    &            , time_output, dxi, deta, dzeta_c, dzeta_t, dzeta_r, z_mar&
    &            , ndat2d, ndat3d, n_output)
        
!-------- Main loop --------
    CALL SICO_MAIN_LOOP_D(delta_ts, glac_index, mean_accum, dtime, &
    &            dtime_temp, dtime_wss, dtime_out, dtime_ser, time, &
    &            time_init, time_end, time_output, dxi, deta, dzeta_c, &
    &            dzeta_t, dzeta_r, z_mar, ndat2d, ndat3d, n_output)
    CALL COST_FINAL_D()
     
    CALL SICO_END()

#endif /* ALLOW_TAP_{ADJ,TLM} */

#if defined(ALLOW_GENCTRL)
    call ad_output()
#endif /* ALLOW_GENCTRL */

#endif /* ALLOW_TAPENADE */
end program tapenade_main
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                       End of tapenade_main.F90
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
