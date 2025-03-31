!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Specification file for control variables ad_repo_ant64_bm3_ss25ka.h
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

#define XX_GENARR2D_VARS_ARR    [ character(CTRL_STRLENGTH) ::\
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
                                         'xx_enh_fact_da_dummy2d_scalar',\
                                         'xx_enh_intg_da_dummy2d_scalar',\

!-------- Settings for genctrl I/O --------

!#define AD_INPUT_PATH '/home/shreyas/update_to_develop_sicopolis/sicopolis/src/subroutines/tapenade/ad_io'
!        Absolute input path to read xx_gen* fields

#define AD_OUTPUT_PATH '/home/shreyas/update_to_develop_sicopolis/sicopolis/src/subroutines/tapenade/ad_io'
!       Absolute output path to write xx_gen* fields
