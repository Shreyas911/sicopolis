!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Specification file for control variables ad_repo_grl16_bm5_ss25ka.h
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------- Basic settings --------

#define CTRL_STRLENGTH		128
!         	The fixed string length for individual entries of arrays
!        	like xx_genarr2d_vars, xx_genarr2d_bounds.

#define NUM_CTRL_GENARR2D	3
!          	Number of 2D time-invariant control variables

#define NUMCTRLPROCARR2D	2
!		Maximum number of preprocessing steps for ctrl variables		

#define XX_GENARR2D_VARS_ARR	(/ 'xx_c_slide_init', \
				   				   'xx_q_geo       ', \
			  	   				   'xx_H           ' /)
!		List of 2D time-invariant control variables

#define XX_GENARR2D_PREPROC_ARR    (/ 'nnnnnnnnn,bounds', \
			              'nnnnnnnnn,nnnnnn', \
				      'bounds   ,nnnnnn' /) 
!               Define preprocessing steps for ctrl variables
!		Fill with nnnn... to ensure same length

#define XX_GENARR2D_BOUNDS_ARR  (/ ' ', \
				   ' ', \
			    	   ' ' /)
!		2D mask for the 2D time-invariant control variables
!		If empty, defaults to 0, IMAX, 0, JMAX

#define XX_GENARR2D_LOG10INITVAL_ARR (/ -6.67172541073, \
				         0.00000000000, \
					 0.00000000000 /)
!		log10initval if preproc=log10ctrl
!               Otherwise has no effect
!		WARNING: If for example using for c_slide_init
!		Set C_SLIDE == 0.0 in the HEADER file

#define XX_GENARR2D_WEIGHT_ARR       (/ 'w_c_slide_init.dat', \
				        				'w_q_geo.dat       ', \
										'w_H.dat           ' /)
!		Weight files to be used if preproc=scaling
!		Otherwise defaults to weight 1.0 for all indices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define NUM_CTRL_GENARR3D	1
!          	Number of 3D time-invariant control variables

#define NUMCTRLPROCARR3D	1
!		Maximum number of preprocessing steps for ctrl variables		

#define XX_GENARR3D_VARS_ARR	(/ 'xx_temp_c' /)
!		List of 3D time-invariant control variables

#define XX_GENARR3D_PREPROC_ARR    (/ 'bounds' /)
!               Define preprocessing steps for ctrl variables
!		Fill with nnnn... to ensure same length

#define XX_GENARR3D_BOUNDS_ARR  (/ ' ' /)
!		3D mask for the 3D time-invariant control variables
!		If empty, defaults to 0, IMAX, 0, JMAX, 0, KCMAX 

#define XX_GENARR3D_LOG10INITVAL_ARR (/ -1 /)
!		log10initval if preproc=log10ctrl
!               Otherwise has no effect
!		WARNING: If for example using for temp_c
!		Set TEMP_INIT_VALUE (or some flag) == 0.0 in the HEADER file

#define XX_GENARR3D_WEIGHT_ARR       (/ 'w_temp_c.dat' /)
!		Weight files to be used if preproc=scaling
!		Otherwise defaults to weight 1.0 for all indices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define NUM_CTRL_GENTIM2D	1
!          	Number of 2D time-varying control variables

#define NUMCTRLPROCTIM2D	1
!		Maximum number of preprocessing steps for ctrl variables		

#define XX_GENTIM2D_VARS_ARR	(/ 'xx_temp' /)
!		List of 3D time-varying control variables

#define XX_GENTIM2D_PERIOD	10.0
!		Time period for gentim2d

#define ADNMAX			4
!		(TIME_END0-TIME_INIT0)/XX_GENTIM2D_PERIOD 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#define AD_INPUT_PATH '/home/shreyas/update_to_develop_sicopolis/sicopolis/src/subroutines/tapenade/ad_io'
!		 Absolute input path to read xx_gen* fields

#define AD_OUTPUT_PATH '/home/shreyas/update_to_develop_sicopolis/sicopolis/src/subroutines/tapenade/ad_io'
!		Absolute output path to write xx_gen* fields