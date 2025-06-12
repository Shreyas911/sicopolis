module ctrl_map_gentim_m

  use sico_types_m
  use sico_variables_m
  use error_m

  implicit none

#ifdef DO_CTRL_GENTIM2D
  public :: ctrl_map_ini_gentim2d, ctrl_map_gentim2d

contains

  subroutine ctrl_map_ini_gentim2d

    implicit none

    integer(i4b)                :: ctrl_index
    integer(i4b)                :: igen_delta_tda

#ifdef XX_GENTIM2D_VARS_ARR
    xx_gentim2d_vars            = XX_GENTIM2D_VARS_ARR
#endif

    igen_delta_tda = 0

#ifdef XX_GENTIM2D_VARS_ARR
    do ctrl_index = 1, NUM_CTRL_GENTIM2D
      
      if (trim(adjustl(xx_gentim2d_vars(ctrl_index))) .EQ. 'xx_delta_tda') then
        igen_delta_tda = ctrl_index
#if (TSURFACE > 4)
        errormsg = ' >>> ctrl_map_ini_gentim2d: ' &
        //'delta_tda as a control param is only compatible with TSURFACE == 1, 2, 3, 4!'
        call error(errormsg)
#endif
      else
        errormsg = ' >>> ctrl_map_ini_gentim2d: ' &
          //"This control variable is not in the gentim2d setup yet!"
        call error(errormsg)
      end if
  
    end do

#if (TSURFACE <= 4)
    if (igen_delta_tda .GT. 0) then
      call ctrl_map_gentim2d(delta_tda, igen_delta_tda)
    end if
#endif
#endif
    
  end subroutine ctrl_map_ini_gentim2d

  subroutine ctrl_map_gentim2d(fld, iarr)

    implicit none  

    real(dp), dimension(0:NTDAMAX,0:JMAX,0:IMAX) :: fld
    integer(i4b)                                 :: iarr, k2

    fld = fld + xx_gentim2d(iarr,:,:,:)

  end subroutine ctrl_map_gentim2d
#endif

end module ctrl_map_gentim_m 
