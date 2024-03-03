#include RUN_SPECS_HEADER

subroutine sicopolis_tapenade(delta_ts, glac_index, &
                      mean_accum, &
                      dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
                      time, time_init, time_end, time_output, &
                      dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                      z_mar, &
                      ndat2d, ndat3d, n_output)

use cost_m
use sico_types_m
use sico_variables_m
use sico_vars_m
use sico_init_m
use sico_main_loop_m
use sico_end_m

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

call sico_init(delta_ts, glac_index, &
     mean_accum, &
     dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
     time, time_init, time_end, time_output, &
     dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
     z_mar, &
     ndat2d, ndat3d, n_output)

!-------- Main loop --------

call sico_main_loop(delta_ts, glac_index, &
     mean_accum, &
     dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
     time, time_init, time_end, time_output, &
     dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
     z_mar, &
     ndat2d, ndat3d, n_output)

call cost_final()
call sico_end()

end subroutine sicopolis_tapenade

