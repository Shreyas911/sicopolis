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
    real(dp)                                   :: dtime, dtime_temp, &
                                                dtime_wss, dtime_out, dtime_ser
    real(dp)                                   :: time, time_init, time_end
    real(dp), dimension(100)                   :: time_output
    real(dp)                                   :: dxi, deta, dzeta_c, &
                                                dzeta_t, dzeta_r

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
#if !defined(ALLOW_TAP_ADJ_AT_ACTION)
    fcb = 1.
#endif
    call SICOPOLIS_TAPENADE_B(dtime, &
    & dtime_temp, dtime_wss, dtime_out, dtime_ser, time, time_init, time_end, &
    & time_output, dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
    & ndat2d, ndat3d, n_output)
#if defined(ALLOW_TAP_ADJ_PROF)
    call adstack_showpeaksize()
    call adstack_showtotaltraffic()
    call adprofileadj_showprofiles()
#endif

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
    do p = 1, points !@ python_automated_tlm limited_or_block_or_full_or_scalar @
        i = ipoints(p)
        j = jpoints(p)

!@ python_automated_tlm dep_vard set 0 @

!@ python_automated_tlm_before dep_vard set 1 @

        CALL SICO_INIT_D(dtime, dtime_temp, &
        &            dtime_wss, dtime_out, dtime_ser, time, time_init, time_end, &
        &            time_output, dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
        &            ndat2d, ndat3d, n_output)

!@ python_automated_tlm_after dep_vard set 1 @
      
!-------- Main loop --------
        CALL SICO_MAIN_LOOP_D(dtime, &
        &                 dtime_temp, dtime_wss, dtime_out, dtime_ser, time, &
        &                 time_init, time_end, time_output, dxi, deta, dzeta_c, &
        &                 dzeta_t, dzeta_r, ndat2d, ndat3d, n_output)
        CALL COST_FINAL_D()
   
        CALL SICO_END()

        ! Initialize compatible fields to 0
        ! 2D fields
        delta_tda_const = 0.0 ! Not compatible with TSURFACE>4
        q_geo           = 0.0
        c_slide_init    = 0.0
        p_weert         = 0.0
        q_weert         = 0.0
#if (ENHMOD==1 || ENHMOD==2 || ENHMOD==3)
        enh_fact_da_dummy2d_scalar = 0.0 ! Only used as a scalar ! Not compatible with ANF_DAT==3
#endif
#if (ENHMOD==2 || ENHMOD==3)
        enh_intg_da_dummy2d_scalar = 0.0 ! Only used as a scalar ! Not compatible with ANF_DAT==3
#endif
#if (FLOW_LAW==1)
        n_glen_da_dummy2d_scalar   = 0.0 ! Only used as a scalar ! Only compatible with :
                                         ! ((!defined(N_POWER_LAW_INT) OR N_POWER_LAW_INT < 0) && (defined(N_POWER_LAW_REAL)))
#endif
        H               = 0.0
        zs              = 0.0 ! Not compatible with ANF_DAT==2
        zl              = 0.0 ! Not compatible with ANF_DAT==2
        zl0             = 0.0
        zb              = 0.0 ! Not compatible with ANF_DAT==2, or if ANF_DAT==1 and ZB_PRESENT_FILE is 'none' or not defined
#if (ACCSURFACE==2 || ACCSURFACE==3)
        gamma_s         = 0.0
#endif
#if (ABLSURFACE==1 || ABLSURFACE==2)
        s_stat          = 0.0 ! Not compatible with SOLID_PRECIP==3, would need to also activate s_stat_solid_precip in that case.
        beta1           = 0.0
        beta2           = 0.0
        Pmax            = 0.0
        mu              = 0.0
#endif
#if (defined(GRL) && DISC>0)
        c_dis_da        = 0.0
#endif
#if (defined(PARAM_RHO_A))
        RHO_A           = 0.0
#endif
#if (REBOUND==1 || REBOUND==2)
        time_lag_asth   = 0.0
#endif
#if (REBOUND==2)
        flex_rig_lith   = 0.0
#endif
        ! 3D fields
        omega_c         = 0.0 ! Only compatible with ANF_DAT==3
        temp_c          = 0.0 ! Not compatible with ANF_DAT==1 and TEMP_INIT==5
        age_c           = 0.0 ! Not compatible with ANF_DAT==1 and TEMP_INIT==5
        temp_r          = 0.0 ! Only compatible with ANF_DAT==3
        delta_tda       = 0.0 ! Not compatible with TSURFACE>4

!@ python_automated_tlm IO write @
    end do ! (close loop over points)

!@ python_automated_tlm IO end @

#elif (defined(ALLOW_TAP_TLM) && defined(ALLOW_GENCTRL))

    CALL SICO_INIT_D(dtime, dtime_temp, &
    &            dtime_wss, dtime_out, dtime_ser, time, time_init, time_end, &
    &            time_output, dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
    &            ndat2d, ndat3d, n_output)
        
!-------- Main loop --------
    CALL SICO_MAIN_LOOP_D(dtime, &
    &            dtime_temp, dtime_wss, dtime_out, dtime_ser, time, &
    &            time_init, time_end, time_output, dxi, deta, dzeta_c, &
    &            dzeta_t, dzeta_r, ndat2d, ndat3d, n_output)

#if !defined(ALLOW_TAP_TLM_A_ACTION)
    CALL COST_FINAL_D()
#endif
     
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
