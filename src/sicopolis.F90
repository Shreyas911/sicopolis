!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Program :  s i c o p o l i s
!             (SImulation COde for POLythermal Ice Sheets)
!
#define       MODEL_SICOPOLIS
#define       VERSION '24'
#define       DATE    '2024-02-27'
!
!! Main program of SICOPOLIS.
!!
!!##### Authors
!!
!! [SICOPOLIS Authors](https://sicopolis.readthedocs.io/en/latest/introduction.html#authorship)
!!
!! The model employs either hybrid shallow-ice-shelfy-stream dynamics
!! (Bernales et al. 2017a,b) or the shallow-ice approximation for
!! grounded ice, and the shallow-shelf approximation for floating ice
!! (e.g., Greve and Blatter 2009). It is coded in Fortran and uses finite
!! difference discretization on a staggered Arakawa C grid, the velocity
!! components being taken between grid points. A variety of different
!! thermodynamics solvers are available, namely the polythermal two-layer
!! method, two versions of the one-layer enthalpy method, the cold-ice
!! method and the isothermal method (Greve and Blatter 2016).
!!
!! The coding is based on a low-tech, ease-of-use philosophy. All
!! structures are kept as simple as possible, and advanced coding
!! techniques are only employed where it is deemed appropriate. The use
!! of external libraries is kept at an absolute minimum, which makes the
!! installation very easy and fast.
!!
!! References:
!! @li <https://sicopolis.readthedocs.io/en/latest/references.html>
!!
!! Resources:
!! @li Model website: <https://www.sicopolis.net/>
!! @li User manual: <https://sicopolis.readthedocs.io/>
!! @li GitLab repository: <https://gitlab.awi.de/sicopolis/sicopolis/>
!! @li SICOPOLIS community @ Zenodo:
!!     <https://zenodo.org/communities/sicopolis/>
!!
!! @section Copyright
!!
!! Copyright 2009-2024 SICOPOLIS Authors\n
!! (<https://sicopolis.readthedocs.io/en/latest/introduction.html#authorship>)
!!
!! @section License
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
!!
!! @file
!!
!! Main program file of SICOPOLIS.
!!
!! @section Copyright
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

#if defined(ALLOW_TAPENADE)
#if defined(ALLOW_GENCTRL)
#include "ad_specs.h"
#endif
#endif

!-------- Include run specification header --------

#include RUN_SPECS_HEADER

!-------- Include header for the Library of Iterative Solvers Lis
!                                               (only if required) --------

#if !defined(ALLOW_TAPENADE) /* Normal */
#if (CALCTHK==3 || CALCTHK==6 || MARGIN==3 || DYNAMICS==2)
#include "lisf.h"
#endif
#endif /* Normal */

!-------- Include modules --------

#if !defined(ALLOW_TAPENADE) /* Normal */
#include "subroutines/general/sico_types_m.F90"
#endif /* Normal */
#include "subroutines/general/sico_variables_m.F90"

#if (defined(ANT))
#include "subroutines/ant/sico_vars_m.F90"
#elif (defined(ASF))
#include "subroutines/asf/sico_vars_m.F90"
#elif (defined(EISMINT))
#include "subroutines/eismint/sico_vars_m.F90"
#elif (defined(GRL))
#include "subroutines/grl/sico_vars_m.F90"
#elif (defined(NHEM))
#include "subroutines/nhem/sico_vars_m.F90"
#elif (defined(NMARS))
#include "subroutines/nmars/sico_vars_m.F90"
#elif (defined(SCAND))
#include "subroutines/scand/sico_vars_m.F90"
#elif (defined(SMARS))
#include "subroutines/smars/sico_vars_m.F90"
#elif (defined(TIBET))
#include "subroutines/tibet/sico_vars_m.F90"
#elif (defined(XYZ))
#include "subroutines/xyz/sico_vars_m.F90"
#endif

#include "subroutines/general/error_m.F90"

#if (defined(ALLOW_GRDCHK) || defined(ALLOW_TAPENADE)) /* Tapenade */
#include "subroutines/tapenade/cost/cost_m.F90"
#if defined(ALLOW_GENCTRL)
#include "subroutines/tapenade/ctrl/ctrl_map_genarr_m.F90"
#include "subroutines/tapenade/ctrl/ctrl_init_genarr_m.F90"
#include "subroutines/tapenade/ctrl/ctrl_map_gentim_m.F90"
#endif
#endif /* Tapenade */

#include "subroutines/general/ice_material_properties_m.F90"
#include "subroutines/general/stereo_proj_m.F90"
#include "subroutines/general/metric_m.F90"
#include "subroutines/general/sico_maths_m.F90"

#if !defined(ALLOW_TAPENADE) /* Normal */
#include "subroutines/general/compare_float_m.F90"
#endif /* Normal */

#include "subroutines/general/nc_check_m.F90"

#include "subroutines/general/read_m.F90"

#include "subroutines/general/mask_update_sea_level_m.F90"
#include "subroutines/general/flag_update_gf_gl_cf_m.F90"
#include "subroutines/general/pdd_m.F90"

#if (MARGIN==2 || MARGIN==3)
#include "subroutines/general/calving_m.F90"
#endif

#if (defined(GRL))
#if (DISC>0)
#include "subroutines/grl/discharge_workers_m.F90"
#endif
#endif

#include "subroutines/general/calc_enhance_m.F90"

#include "subroutines/general/calc_vxy_m.F90"
#include "subroutines/general/calc_vz_m.F90"
#include "subroutines/general/calc_dxyz_m.F90"

#include "subroutines/general/calc_gia_m.F90"
#include "subroutines/general/topograd_m.F90"
#include "subroutines/general/calc_thk_m.F90"

#include "subroutines/general/enth_temp_omega_m.F90"

#if (CALCMOD==0 || CALCMOD==1 || CALCMOD==-1)
#include "subroutines/general/calc_temp_m.F90"
#elif (CALCMOD==2 || CALCMOD==3)
#include "subroutines/general/calc_temp_enth_m.F90"
#endif

#if (BASAL_HYDROLOGY==1)
#include "subroutines/general/hydro_m.F90"
#endif

#if (defined(NMARS) || defined(SMARS))
#include "subroutines/general/mars_instemp_m.f90"
#endif

#include "subroutines/general/calc_temp_melt_bas_m.F90"
#include "subroutines/general/calc_bas_melt_m.F90"
#include "subroutines/general/calc_thk_water_bas_m.F90"

#include "subroutines/general/output_m.F90"

#if (defined(ANT))
#include "subroutines/ant/boundary_m.F90"
#elif (defined(ASF))
#include "subroutines/asf/boundary_m.F90"
#elif (defined(EISMINT))
#include "subroutines/eismint/boundary_m.F90"
#elif (defined(GRL))
#include "subroutines/grl/boundary_m.F90"
#elif (defined(NHEM))
#include "subroutines/nhem/boundary_m.F90"
#elif (defined(NMARS))
#include "subroutines/nmars/boundary_m.F90"
#elif (defined(SCAND))
#include "subroutines/scand/boundary_m.F90"
#elif (defined(SMARS))
#include "subroutines/smars/boundary_m.F90"
#elif (defined(TIBET))
#include "subroutines/tibet/boundary_m.F90"
#elif (defined(XYZ))
#include "subroutines/xyz/boundary_m.F90"
#endif

#include "subroutines/general/init_temp_water_age_m.F90"
#if (defined(ANT))
#include "subroutines/ant/sico_init_m.F90"
#elif (defined(ASF))
#include "subroutines/asf/sico_init_m.F90"
#elif (defined(EISMINT))
#include "subroutines/eismint/sico_init_m.F90"
#elif (defined(GRL))
#include "subroutines/grl/sico_init_m.F90"
#elif (defined(NHEM))
#include "subroutines/nhem/sico_init_m.F90"
#elif (defined(NMARS))
#include "subroutines/nmars/sico_init_m.F90"
#elif (defined(SCAND))
#include "subroutines/scand/sico_init_m.F90"
#elif (defined(SMARS))
#include "subroutines/smars/sico_init_m.F90"
#elif (defined(TIBET))
#include "subroutines/tibet/sico_init_m.F90"
#elif (defined(XYZ))
#include "subroutines/xyz/sico_init_m.F90"
#endif

#include "subroutines/general/sico_main_loop_m.F90"
#include "subroutines/general/sico_end_m.F90"

#if (defined(ALLOW_GRDCHK) || defined(ALLOW_TAPENADE)) /* Tapenade */
#include "subroutines/tapenade/src/tapenade_m.F90"
#endif /* Tapenade */

!-------------------------------------------------------------------------------
!> Main program of SICOPOLIS.
!-------------------------------------------------------------------------------
program sicopolis

!-------- Declaration of variables --------

use sico_types_m
use sico_variables_m
use sico_vars_m

use sico_init_m
use sico_main_loop_m
use sico_end_m

#if defined(ALLOW_GRDCHK) /* Tapenade */
use tapenade_m, only: grdchk_main
#endif /* Tapenade */

#if defined(ALLOW_TAPENADE) /* Tapenade */
use tapenade_m, only: adjoint_master
#endif /* Tapenade */

implicit none

integer(i4b) :: ndat2d
   !! Counter for the time-slice files

integer(i4b) :: ndat3d
   !! Counter for the time-slice files

integer(i4b) :: n_output
   !! Number of time-slice files to be produced (for OUTPUT==2, 3)

real(dp) :: delta_ts
   !! Time-dependent surface-temperature anomaly

real(dp) :: glac_index
   !! Time-dependent glacial index

real(dp) :: mean_accum
   !! Mean present accumulation over land

real(dp) :: dtime
   !! Time step for computation of velocity and topography

real(dp) :: dtime_temp
   !! Time step for computation of temperature, water content
   !! and age of the ice

real(dp) :: dtime_wss
   !! Time step for computation of isostatic steady-state displacement
   !! of the lithosphere (ELRA model)

real(dp) :: dtime_out
   !! Time step for output of time-slice files

real(dp) :: dtime_ser
   !! Time step for writing of data in time-series file

real(dp) :: time
   !! Current time of simulation

real(dp) :: time_init
   !! Initial time of simulation

real(dp) :: time_end
   !! Final time of simulation

real(dp), dimension(100) :: time_output
   !! Specified times for output of time-slice files (for OUTPUT==2, 3)

real(dp) :: dxi
   !! Grid spacing in x-direction

real(dp) :: deta
   !! Grid spacing in y-direction

real(dp) :: dzeta_c
   !! Grid spacing in z-direction in the upper (kc) ice domain
   !! (in sigma-coordinate zeta_c)

real(dp) :: dzeta_t
   !! Grid spacing in z-direction in the lower (kt) ice domain
   !! (in sigma-coordinate zeta_t)

real(dp) :: dzeta_r
   !! Grid spacing in z-direction in the bedrock (kr) domain
   !! (in sigma-coordinate zeta_r)

real(dp) :: z_mar
   !! Minimum bedrock (sea bed) elevation allowed to be covered by marine ice

#if (!defined(ALLOW_GRDCHK) && !defined(ALLOW_TAPENADE)) /* Normal */

!-------- Initializations --------

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

!-------- Endings --------

call sico_end()

!-------- End of program --------

#else /* Tapenade */
#if (defined(ALLOW_GRDCHK))
call grdchk_main()
#endif
#if (defined(ALLOW_TAPENADE))
call adjoint_master()
#endif
#endif /* Normal vs. Tapenade */

end program sicopolis

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                       End of sicopolis.F90
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
