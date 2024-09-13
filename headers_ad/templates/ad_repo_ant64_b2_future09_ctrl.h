!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Specification file for control variables ad_repo_ant64_b2_future09_ctrl.h
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#undef ALLOW_GENCTRL
!       Flag to enable specific code for genctrl setup

#undef ALLOW_GENCTRL_BEFORE_SICO_INIT
!       Flag to enable activation of genctrl before sico_init

#define DO_CTRL_GENARR2D
#define DO_CTRL_GENARR3D
#define DO_CTRL_GENTIM2D
!       Flags to enable specific codes for various types of genctrl

#define CTRL_STRLENGTH   128
!       The fixed string length for individual entries of arrays
!       like xx_genarr2d_vars, xx_genarr2d_bounds.

!-------- Settings for genarr2D --------

#define NUM_CTRL_GENARR2D 3
!       Number of 2D time-invariant control variables

#define NUMCTRLPROCARR2D  1
!       Maximum number of preprocessing steps for ctrl variables

#define XX_GENARR2D_VARS_ARR    [ character(CTRL_STRLENGTH) ::\
                                    'xx_c_slide_init',\
                                    'xx_q_geo',\
                                    'xx_H' ]
!       List of 2D time-invariant control variables

#define XX_GENARR2D_PREPROC_ARR    [ character(CTRL_STRLENGTH) ::\
                                       'none',\
                                       'none',\
                                       'none' ]
!       Define preprocessing steps for ctrl variables

!-------- Settings for genarr3D --------

#define NUM_CTRL_GENARR3D 1
!       Number of 3D time-invariant control variables

#define NUMCTRLPROCARR3D  1
!       Maximum number of preprocessing steps for ctrl variables

#define XX_GENARR3D_VARS_ARR    [ character(CTRL_STRLENGTH) :: 'xx_temp_c' ]
!       List of 3D time-invariant control variables

#define XX_GENARR3D_PREPROC_ARR [ character(CTRL_STRLENGTH) :: 'none' ]
!       Define preprocessing steps for ctrl variables

!-------- Settings for gentim2D --------

#define NUM_CTRL_GENTIM2D 1
!       Number of 2D time-varying control variables

#define NUMCTRLPROCTIM2D  1
!       Maximum number of preprocessing steps for ctrl variables

#define XX_GENTIM2D_VARS_ARR [ character(CTRL_STRLENGTH) :: 'xx_delta_tda' ]
!       List of 3D time-varying control variables

!-------- Settings for genctrl I/O --------

#undef AD_INPUT_PATH '/home/shreyas/update_to_develop_sicopolis/sicopolis/src/subroutines/tapenade/ad_io'
!      Absolute input path to read xx_gen* fields

#define AD_OUTPUT_PATH '/home/shreyas/update_to_develop_sicopolis/sicopolis/src/subroutines/tapenade/ad_io'
!       Absolute output path to write xx_gen* fields
