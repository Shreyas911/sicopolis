!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Specification file for control variables ad_repo_grl40_bm5_paleo17a_BH0_AC.h
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#undef ALLOW_GENCTRL
!       Flag to enable specific code for genctrl setup

#undef ALLOW_GENCTRL_BEFORE_SICO_INIT
!       Flag to enable activation of genctrl before sico_init

#define CTRL_STRLENGTH   128
!       The fixed string length for individual entries of arrays
!       like xx_genarr2d_vars, xx_genarr2d_bounds.

#define DO_CTRL_GENARR2D
#define DO_CTRL_GENARR3D
#define DO_CTRL_GENTIM2D
!       Flags to enable specific codes for various types of genctrl

#define NUM_CTRL_GENARR2D 9
#define NUM_CTRL_GENARR3D 2
#define NUM_CTRL_GENTIM2D 1
!       Number of control variables,
!       Has to be 0 for all 3 as dummy value when not in use.

!-------- Settings for genarr2D --------

#define XX_GENARR2D_VARS_ARR         [ character(CTRL_STRLENGTH) ::\
                                         'xx_c_slide_init',\
                                         'xx_q_geo',\
                                         'xx_H',\
                                         'xx_gamma_s',\
                                         'xx_s_stat',\
                                         'xx_beta1',\
                                         'xx_beta2',\
                                         'xx_Pmax',\
                                         'xx_mu' ]
!       List of 2D time-invariant control variables

!#define NUMCTRLPROCARR2D 1
!!       Maximum number of preprocessing steps for ctrl variables

!#define XX_GENARR2D_PREPROC_ARR      [ character(CTRL_STRLENGTH) ::\
!                                         'log10initval',\
!                                         'none',\
!                                         'none',\
!                                         'log10initval',\
!                                         'log10initval',\
!                                         'log10initval',\
!                                         'log10initval',\
!                                         'log10initval',\
!                                         'log10initval' ]
!!       Define preprocessing steps for ctrl variables

!#define XX_GENARR2D_LOG10INITVAL_ARR [ real :: \
!                                          0.39794000867,  0.0          , 0.0,\
!                                         -1.15206968873,  0.69897000434, 0.43616264704,\
!                                          0.86213137931, -0.22184874962, 0.98746515611 ]
!!       log10initval if preproc=log10ctrl
!!       Otherwise has no effect
!!       WARNING: If for example using for c_slide_init
!!       Set C_SLIDE == 0.0 in the HEADER file

!-------- Settings for genarr3D --------

#define XX_GENARR3D_VARS_ARR    [ character(CTRL_STRLENGTH) ::\
                                    'xx_temp_c',\
                                    'xx_age_c' ]
!       List of 3D time-invariant control variables

!-------- Settings for gentim2D --------

#define XX_GENTIM2D_VARS_ARR [ character(CTRL_STRLENGTH) :: 'xx_delta_tda' ]
!       List of 3D time-varying control variables

!-------- Settings for genctrl I/O --------

!#define AD_INPUT_PATH '/home/shreyas/update_to_develop_sicopolis/sicopolis/src/subroutines/tapenade/ad_io'
!        Absolute input path to read xx_gen* fields

#define AD_OUTPUT_PATH '/home/shreyas/update_to_develop_sicopolis/sicopolis/src/subroutines/tapenade/ad_io'
!       Absolute output path to write xx_gen* fields
