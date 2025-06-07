!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Specification file for control variables ad_repo_grl16_bm5_ss25ka.h
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#undef ALLOW_GENCTRL
!       Flag to enable specific code for genctrl setup

#undef ALLOW_GENCTRL_BEFORE_SICO_INIT
!       Flag to enable activation of genctrl before sico_init

#undef DO_GENCTRL_PRIOR
!       Flag to enable prior regularization

#define CTRL_STRLENGTH   128
!       The fixed string length for individual entries of arrays
!       like xx_genarr2d_vars, xx_genarr2d_bounds.

#define DO_CTRL_GENARR2D
#define DO_CTRL_GENARR3D
#define DO_CTRL_GENTIM2D
!       Flags to enable specific codes for various types of genctrl

#define NUM_CTRL_GENARR2D 3
#define NUM_CTRL_GENARR3D 1
#define NUM_CTRL_GENTIM2D 1
!       Number of control variables,
!       Has to be 0 for all 3 as dummy value when not in use.

!-------- Settings for genarr2D --------

#define XX_GENARR2D_VARS_ARR         [ character(CTRL_STRLENGTH) ::\
                                         'xx_c_slide_init',\
                                         'xx_q_geo',\
                                         'xx_H' ]
!       List of 2D time-invariant control variables

!-------- Settings for genarr3D --------

#define XX_GENARR3D_VARS_ARR    [ character(CTRL_STRLENGTH) :: 'xx_temp_c' ]
!       List of 3D time-invariant control variables

!-------- Settings for gentim2D --------

#define XX_GENTIM2D_VARS_ARR [ character(CTRL_STRLENGTH) :: 'xx_delta_tda' ]
!       List of 3D time-varying control variables

!-------- Settings for genctrl I/O --------

!#define AD_INPUT_PATH '/home/shreyas/update_to_develop_sicopolis/sicopolis/src/subroutines/tapenade/ad_io'
!        Absolute input path to read xx_gen* fields

#define AD_OUTPUT_PATH '/home/shreyas/update_to_develop_sicopolis/sicopolis/src/subroutines/tapenade/ad_io'
!
!       n_glen_da_dummy2d_scalar can be added to genctrl but it is hard-coded in some places, and tuning it changes units of Arrhenius factor A, which is only known for n = 3.
!       WARNING: enh_fact_da_dummy2d_scalar, enh_intg_da_dummy2d_scalar, n_glen_da_dummy2d_scalar are special cases.
!       They are only supposed to be scalars. Illustrating examples below.
!
!       Comment out in ice_material_properties_m for n_glen_da_dummy2d_scalar.
!       n_power_law = 3.0_dp + SUM(n_glen_da_dummy2d_scalar) / SIZE(n_glen_da_dummy2d_scalar)
!       But add this line since n_glen_da_dummy2d_scalar = exp(xx_n_glen_da_dummy2d_scalar) = 3.0
!       n_power_law = SUM(n_glen_da_dummy2d_scalar) / SIZE(n_glen_da_dummy2d_scalar)
!
!       Set ENH_FACT = 0 in header files and replace the following line in calc_enhance_m.
!       enh_c = ENH_FACT
!       with this line.
!       enh_c = SUM(enh_fact_da_dummy2d_scalar) / SIZE(enh_fact_da_dummy2d_scalar) + ENH_FACT
!
!       WARNING: enh_fact_da_dummy2d_scalar, enh_intg_da_dummy2d_scalar are special cases.
!       They currently cannot be used with pickups, hence you do not see them in the list above by default.
!       They need to be accounted for correctly in read_tms_nc, since it is reading enh_c and enh_t from the spinup.
!       You can add an extreme value in read_tms_nc to enh_c and enh_t (say 1000) and verify for yourself that the simulation crashes.
!       This is unlike vx_c, vy_c, vz_c, which get completely overwritten.
!       Absolute output path to write xx_gen* fields
