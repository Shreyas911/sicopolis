!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Specification file for control variables ad_repo_grl16_bm5_ss25ka.h
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------- Basic settings --------

#define CTRL_STRLENGTH		128
!         	The fixed string length for individual entries of arrays
!        	like xx_genarr2d_vars, xx_genarr2d_bounds.

#define NUM_CTRL_GENARR2D	3
!          	Number of 2D time-invariant control variables

#define NUMCTRLPROC		1
!		Maximum number of preprocessing steps for ctrl variables		

#define XX_GENARR2D_VARS_ARR	(/ 'xx_c_slide_init', \
				   'xx_q_geo       ', \
			  	   'xx_zs          ' /)
!		List of 2D time-invariant control variables

#define XX_GENARR2D_PREPROC_ARR    (/ 'log10ctrl', \
			              'nnnnnnnnn', \
				      'bounds   ' /) 
!               Define preprocessing steps for ctrl variables
!		Fill with nnnn... to ensure same length

#define XX_GENARR2D_BOUNDS_ARR  (/ ' ', \
				   ' ', \
			    	   ' ' /)
!		2D mask for the 2D time-invariant control variables
!		If empty, defaults to 1, IMAX+1, 1, JMAX+1

#define XX_GENARR2D_LOG10INITVAL_ARR (/ -6.67172541073, \
				         0.00000000000, \
					 0.00000000000 /)
!		log10initval if preproc=log10ctrl
!               Otherwise has no effect
!		WARNING: If for example using for c_slide_init
!		Set C_SLIDE == 0.0 in the HEADER file

#define XX_GENARR2D_WEIGHT_ARR       (/ 'w_c_slide_init.dat', \
				        'w_q_geo.dat       ', \
					'w_zs.dat          ' /)
!		Weight files to be used if preproc=scaling
!		Otherwise defaults to weight 1.0 for all indices
