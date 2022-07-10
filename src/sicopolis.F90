!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Program :  s i c o p o l i s
!             (SImulation COde for POLythermal Ice Sheets)
!
#define       MODEL_SICOPOLIS
#define       VERSION '5-dev'
#define       DATE    '2022-07-09'
!
!> @mainpage
!!
!! @section Description
!!
!! SICOPOLIS (SImulation COde for POLythermal Ice Sheets) is a 3-d
!! dynamic/thermodynamic model that simulates the evolution of large ice
!! sheets and ice caps. It was originally created by Greve (1997a,b) in a
!! version for the Greenland ice sheet. Since then, SICOPOLIS has been
!! developed continuously and applied to problems of past, present and
!! future glaciation of Greenland, Antarctica, the entire northern
!! hemisphere, the polar ice caps of the planet Mars and others.
!!
!! The model is based on the shallow ice approximation for grounded ice, the
!! shallow shelf approximation for floating ice (e.g., Greve and Blatter 2009)
!! and, optionally, hybrid shallow-ice--shelfy-stream dynamics for ice streams
!! (Bernales et al. 2017a,b). It is coded in Fortran and uses finite difference
!! discretisation on a staggered (Arakawa C) grid, the velocity
!! components being taken between grid points. A variety of different
!! thermodynamics solvers are available, namely the polythermal two-layer
!! method, two versions of the one-layer enthalpy method, the cold-ice
!! method and the isothermal method (Greve and Blatter 2016).
!!
!! The coding is based on an ease-of-use, low-tech philosophy. All
!! structures are kept as simple as possible, and advanced coding
!! techniques are only employed where it is deemed appropriate. The use
!! of external libraries is kept at an absolute minimum. In fact,
!! SICOPOLIS can be run without external libraries at all, which makes
!! the installation very easy and fast.
!!
!! Required model forcing:
!! @li Surface mass balance (precipitation, evaporation, runoff).
!! @li Mean annual air temperature above the ice.
!! @li Eustatic sea level.
!! @li Geothermal heat flux.
!!
!! Main output (as functions of position and time):
!! @li Extent and thickness of the ice sheet.
!! @li Velocity field.
!! @li Temperature field.
!! @li Water-content field (temperate regions).
!! @li Age of the ice.
!! @li Isostatic displacement and temperature of the lithosphere.
!!
!! References:
!! @li Bernales, J., I. Rogozhina, R. Greve and M. Thomas. 2017a.\n
!!     Comparison of hybrid schemes for the combination of
!!     shallow approximations in numerical simulations of the
!!     Antarctic Ice Sheet.\n
!!     Cryosphere 11 (1), 247-265.
!! @li Bernales, J., I. Rogozhina and M. Thomas. 2017b.\n
!!     Melting and freezing under Antarctic ice shelves from a
!!     combination of ice-sheet modelling and observations.\n
!!     Journal of Glaciology 63 (240), 731-744.
!! @li Greve, R. 1997a.\n
!!     A continuum-mechanical formulation for shallow polythermal ice sheets.\n
!!     Phil. Trans. R. Soc. A 355 (1726), 921-974.
!! @li Greve, R. 1997b.\n
!!     Application of a polythermal three-dimensional ice sheet model to the
!!     Greenland ice sheet: Response to steady-state and transient climate
!!     scenarios.\n
!!     J. Climate 10 (5), 901-918.
!! @li Greve, R. and H. Blatter. 2009.\n
!!     Dynamics of Ice Sheets and Glaciers.\n
!!     Springer, Berlin, Germany etc., 287 pp.
!! @li Greve, R. and H. Blatter. 2016.\n
!!     Comparison of thermodynamics solvers in the polythermal ice sheet model
!!     SICOPOLIS.\n
!!     Polar Sci. 10 (1), 11-23.
!! @li SICOPOLIS website: <http://www.sicopolis.net/>
!! @li Changelog:
!!     <https://gitlab.awi.de/sicopolis/sicopolis/> -> Repository -> Commits
!!
!! @section Copyright
!!
!! Copyright 2009-2022 Ralf Greve\n
!! (with contributions by Jorge Bernales, Sebastian Beyer, Heinz Blatter,
!! Reinhard Calov, Thorben Dunse, Eduardo Flandez, Ben Galton-Fenzi,
!! Thomas Goelles, Bjoern Grieger, Philipp Hancke, Patrick Heimbach,
!! Nina Kirchner, Thomas Kleiner, Sascha Knell, Anne Le Brocq, Liz Curry Logan,
!! Sri Hari Krishna Narayanan, Alex Robinson, Fuyuki Saito, Tatsuru Sato,
!! Marius Schaefer, Matthias Scheiter, Oliver J. Stenzel, Malte Thoma,
!! Roland Warner)
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
!! along with SICOPOLIS. If not, see <http://www.gnu.org/licenses/>.
!!
!! @file
!!
!! Main program file of SICOPOLIS.
!!
!! @section Copyright
!!
!! Copyright 2009-2022 Ralf Greve
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
!! along with SICOPOLIS. If not, see <http://www.gnu.org/licenses/>.
!<
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!@ begin openad_extract @

!-------- Include run specification header --------

#include RUN_SPECS_HEADER

!-------- Include header for the Library of Iterative Solvers Lis
!                                               (only if required) --------

#if !defined(ALLOW_OPENAD) /* Normal */
#if (CALCTHK==3 || CALCTHK==6 || MARGIN==3 || DYNAMICS==2)
#include "lisf.h"
#endif
#endif /* Normal */

!-------- Include modules --------

#include "subroutines/general/sico_types_m.F90"
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

#if (defined(ALLOW_GRDCHK) || defined(ALLOW_OPENAD)) /* OpenAD */
#include "subroutines/openad/ctrl_m.F90"
#endif /* OpenAD */

#include "subroutines/general/ice_material_properties_m.F90"
#include "subroutines/general/stereo_proj_m.F90"
#include "subroutines/general/metric_m.F90"

#if (!defined(ALLOW_OPENAD) || defined(ALLOW_OPENAD_DIFFERENTIATE)) /* Normal */
#include "subroutines/general/sico_maths_m.F90"
#endif /* Normal */

#if !defined(ALLOW_OPENAD) /* Normal */
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

#if !defined(ALLOW_OPENAD) /* Normal */
#include "subroutines/general/output_m.F90"
#endif /* Normal */

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

!@ end openad_extract @

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

!@ begin openad_extract @

#if defined(ALLOW_OPENAD) /* OpenAD */
#include "subroutines/openad/sico_main_loop_iter_m.F90"
#include "subroutines/openad/sico_main_loop_wrapper_m.F90"
#endif /* OpenAD */

#include "subroutines/general/sico_main_loop_m.F90"
#include "subroutines/general/sico_end_m.F90"

!@ end openad_extract @

#if (defined(ALLOW_GRDCHK) || defined(ALLOW_OPENAD)) /* OpenAD */
#include "subroutines/openad/openad_m.F90"
#endif /* OpenAD */
!-------------------------------------------------------------------------------
!> Main program of SICOPOLIS.
!<------------------------------------------------------------------------------
program sicopolis

!-------- Description of variables declared locally --------

!  i                 : Index for coordinate xi
!  j                 : Index for coordinate eta
!  kc                : Index for coordinate zeta_c
!  kt                : Index for coordinate zeta_t
!  kr                : Index for coordinate zeta_r
!                      (in the bedrock)
!  ios               : IOSTAT variable for files
!  itercount         : Counter for the time steps
!  ndat2d/3d         : Counters for the time-slice files
!  n_output          : For OUTPUT==2 number of time-slice files
!                      to be produced
!  dH_t_smooth(j,i)  : Amount of smoothing of H_t
!  delta_ts          : Time-dependent surface-temperature variation
!  glac_index        : Time-dependent glacial index
!  forcing_flag      : 1 - forcing by delta_ts, 2 - forcing by glac_index
!  precip_mam_present(j,i) : Measured present spring precipitation
!  precip_jja_present(j,i) : Measured present summer precipitation
!  precip_son_present(j,i) : Measured present autumn precipitation
!  precip_djf_present(j,i) : Measured present winter precipitation
!  temp_mam_present(j,i)   : Present spring surface temperature
!                            from data
!  temp_jja_present(j,i)   : Present summer surface temperature
!                            from data
!  temp_son_present(j,i)   : Present autumn surface temperature
!                            from data
!  temp_djf_present(j,i)   : Present winter surface temperature
!                            from data
!  mean_accum        : Mean present accumulation over land
!  time              : Current time of simulation
!  time_init         : Initial time of simulation
!  time_end          : Final time of simulation
!  dtime             : Time step for computation of velocity and
!                      topography
!  dtime_temp        : Time step for computation of temperature,
!                      water content and age of the ice
!  dtime_wss         : Time step for computation of isostatic steady-state
!                      displacement of the lithosphere (ELRA model)
!  dtime_out         : Time step for output of time-slice files
!  dtime_ser         : Time step for writing of data in time-series
!                      file
!  time_output(n)    : For OUTPUT==2 specified times for output of
!                      time-slice files
!  (.)0              : Quantity (.) before its conversion to MKS units
!  dxi               : Grid spacing in x-direction
!  deta              : Grid spacing in y-direction
!  dzeta_c           : Grid spacing in z-direction in the upper (kc) ice domain
!                      (in sigma-coordinate zeta_c)
!  dzeta_t           : Grid spacing in z-direction in the lower (kt) ice domain
!                      (in sigma-coordinate zeta_t)
!  dzeta_r           : Grid spacing in z-direction in the bedrock (kr) domain
!                      (in sigma-coordinate zeta_r)

!-------- Declaration of variables --------

!@ begin openad_extract @
!openad begin subroutine sicopolis_openad
use sico_types_m
use sico_variables_m
use sico_vars_m
!@ end openad_extract @

use sico_init_m
!@ begin openad_extract @
use sico_main_loop_m
use sico_end_m
!@ end openad_extract @

#if defined(ALLOW_GRDCHK) /* OpenAD */
use openad_m, only: adjoint_master, grdchk_main
#endif /* OpenAD */

#if defined(ALLOW_OPENAD) /* OpenAD */
use OAD_tape
use OAD_rev
use openad_m, only: adjoint_master
#endif /* OpenAD */

implicit none

!@ begin openad_extract @
integer(i4b) :: ndat2d, ndat3d
integer(i4b) :: n_output
real(dp) :: delta_ts, glac_index
real(dp) :: mean_accum
real(dp) :: dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser
real(dp) :: time, time_init, time_end
real(dp), dimension(100) :: time_output
real(dp) :: dxi, deta, dzeta_c, dzeta_t, dzeta_r
real(dp) :: z_mar
!openad sicopolis_independents_cost
!@ end openad_extract @

#if defined(ALLOW_OPENAD) /* OpenAD */
integer(i4b) :: mode
character(len=32) :: arg
logical :: ISPLAIN, ISTAPE, ISADJOINT
#endif /* OpenAD */

#if (!defined(ALLOW_GRDCHK) && !defined(ALLOW_OPENAD)) /* Normal */
!-------- Initialisations --------

call sico_init(delta_ts, glac_index, &
     mean_accum, &
     dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
     time, time_init, time_end, time_output, &
     dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
     z_mar, &
     ndat2d, ndat3d, n_output)

!@ begin openad_extract @
!-------- Main loop --------

call sico_main_loop(delta_ts, glac_index, &
     mean_accum, &
     dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
     time, time_init, time_end, time_output, &
     dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
     z_mar, &
     ndat2d, ndat3d, n_output)

!openad end subroutine sicopolis_openad

!@ end openad_extract @
!-------- Endings --------

call sico_end()

!-------- End of program --------

#else /* OpenAD */

call adjoint_master

#endif /* Normal vs. OpenAD */

end program sicopolis

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                       End of sicopolis.F90
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
