!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Specification file for control variables CTRL_VARS.h
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------- Basic settings --------

#define CTRL_STRLENGTH		15
!         	The fixed string length for individual entries of arrays
!        	like xx_genarr2d_vars, xx_genarr2d_bounds.

#define NUM_CTRL_GENARR2D	3
!          	Number of 2D time-invariant control variables

#define XX_GENARR2D_VARS_ARR	(/ 'xx_c_slide     ', \
				   'xx_zs          ', \
			  	   'xx_q_geo       ' /)
!		List of 2D time-invariant control variables

#define XX_GENARR2D_BOUNDS_ARR  (/ '               ', \
				   '               ', \
			    	   '               ' /)
!		2D mask for the 2D time-invariant control variables
!		If empty or none, defaults to 1, IMAX+1, 1, JMAX+1
