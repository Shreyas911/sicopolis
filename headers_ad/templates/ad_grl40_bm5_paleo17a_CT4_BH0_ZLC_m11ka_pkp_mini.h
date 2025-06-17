!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Specification file for control variables ad_repo_grl40_bm5_paleo17a_CT4_BH0_ZLC_m11ka_pkp_mini.h
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

#define NUM_CTRL_GENARR2D 23
#define NUM_CTRL_GENARR3D 3
#define NUM_CTRL_GENTIM2D 1
!       Number of control variables,
!       Has to be 0 for all 3 as dummy value when not in use.

!-------- Settings for genarr2D --------

#define XX_GENARR2D_VARS_ARR         [ character(CTRL_STRLENGTH) ::\
                                         'xx_c_slide_init',\
                                         'xx_delta_tda_const',\
                                         'xx_c_dis_da',\
                                         'xx_q_geo',\
                                         'xx_H',\
                                         'xx_gamma_s',\
                                         'xx_s_stat',\
                                         'xx_beta1',\
                                         'xx_beta2',\
                                         'xx_Pmax',\
                                         'xx_mu',\
                                         'xx_RHO_A',\
                                         'xx_time_lag_asth',\
                                         'xx_flex_rig_lith',\
                                         'xx_p_weert',\
                                         'xx_q_weert',\
                                         'xx_enh_fact_da_dummy2d_scalar',\
                                         'xx_enh_intg_da_dummy2d_scalar',\
                                         'xx_n_glen_da_dummy2d_scalar',\
                                         'xx_zs',\
                                         'xx_zl',\
                                         'xx_zl0',\
                                         'xx_zb' ]
!       List of 2D time-invariant control variables

!#define NUMCTRLPROCARR2D 1
!!       Maximum number of preprocessing steps for ctrl variables

!#define XX_GENARR2D_PREPROC_ARR      [ character(CTRL_STRLENGTH) ::\
!                                         'log10ctrl',\
!                                         'none',\
!                                         'log10ctrl',\
!                                         'none',\
!                                         'none',\
!                                         'log10ctrl',\
!                                         'log10ctrl',\
!                                         'log10ctrl',\
!                                         'log10ctrl',\
!                                         'log10ctrl',\
!                                         'log10ctrl',\
!                                         'log10ctrl',\
!                                         'log10ctrl',\
!                                         'log10ctrl',\
!                                         'log10ctrl',\
!                                         'log10ctrl',\
!                                         'log10ctrl',\
!                                         'log10ctrl',\
!                                         'log10ctrl',\
!                                         'none',\
!                                         'none',\
!                                         'none',\
!                                         'none' ]
!!       Define preprocessing steps for ctrl variables
!!       log10ctrl only works if defined(ALLOW_GENCTRL_BEFORE_SICO_INIT)
!!       WARNING: If for example using for c_slide_init
!!       Set C_SLIDE == 0.0 in the HEADER file

!!       WARNING: If using for spatially varying fields, ensure you manually comment some lines.
!!                This is to prevent fields from getting added twice.	
!!	 
!!       WARNING: If for example using for H.
!!                Just to illustrate, assuming the right ANF_DAT flags are in use.
!!
!!       Comment out in read_tms_nc if reading from file.
!!       H(j,i)    = H(j,i) + real(H_conv(i,j),dp)
!!       Comment out in topography1 if initializing.
!!       H(j,i)   = H(j,i) + zs(j,i)-zm(j,i)
!!                
!!       WARNING: Take care of units while commenting out, showing for a 3D field age_c.
!!                Just to illustrate, assuming the right ANF_DAT flags are in use.
!!                This should be done for XX_GENARR3D_PREPROC_ARR of course, the comment is here for now.
!!
!!       Comment out in read_tms_nc if reading from file.
!!       age_c(kc,j,i)   = (age_c(kc,j,i) + real(age_c_conv(i,j,kc),dp))*year2sec
!!       But to preserve the unit change, add a line.
!!       age_c(kc,j,i) = real(age_c(kc,j,i),dp)*year2sec
!!
!!       Comment out in init_age if initializing.
!!       age_c = (age_c + 15000.0_dp)*year2sec
!!       But to preserve the unit change, add a line.
!!       age_c(kc,j,i) = real(age_c(kc,j,i),dp)*year2sec

!#define XX_GENARR2D_LOG10INITVAL_ARR [ real :: \
!                                          0.92941892571,  0.0,         , 4.19476402411,\
!                                          0.0          ,  0.0,\
!                                         -1.15206968873,  0.69897000434, 0.43616264704,\
!                                          0.86213137931, -0.22184874962, 0.98746515611,\
!                                          3.51851393988,  3.47712125472, 25.0000000000,\
!                                          0.47712125472,  0.30102999566, 0.47712125472, 0.0,\
!                                          0.0, 0.0, 0.0, 0.0 ]
!!       log10initval is used only if preproc=log10ctrl and AD_INPUT_PATH is not defined.
!!       Has no effect (not even read) if AD_INPUT_PATH is defined.
!!       WARNING: If for example using for c_slide_init
!!       Set C_SLIDE == 0.0 in the HEADER file

!-------- Settings for genarr3D --------

#define XX_GENARR3D_VARS_ARR    [ character(CTRL_STRLENGTH) ::\
                                    'xx_temp_c',\
                                    'xx_omega_c',\
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
