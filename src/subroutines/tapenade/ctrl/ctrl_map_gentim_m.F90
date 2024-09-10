module ctrl_map_gentim_m

  use sico_types_m
  use sico_variables_m
  use error_m

  implicit none

#ifdef DO_CTRL_GENTIM2D
  public :: ctrl_map_ini_gentim2d, ctrl_map_gentim2d
#endif

contains

#ifdef DO_CTRL_GENTIM2D
  subroutine ctrl_map_ini_gentim2d(time_init, dtime, itercount)

    implicit none

    real(dp), intent(in)        :: time_init, dtime
    integer(i4b), intent(in)    :: itercount
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
      else
        errormsg = ' >>> ctrl_map_ini_gentim2d: ' &
          //"This control variable is not in the gentim2d setup yet!"
        call error(errormsg)
      end if
  
    end do

    if (igen_temp .GT. 0) then
      call ctrl_map_gentim2d(delta_tda, igen_delta_tda)
    end if
#endif
    
  end subroutine ctrl_map_ini_gentim2d

  subroutine ctrl_map_gentim2d(fld, iarr)

    implicit none  

    real(dp), dimension(0:JMAX,0:IMAX)          :: fld, xx_gen, xx_gen_mask
    integer(i4b)                                :: iarr, k2

    fld = fld + xx_gentim2d(iarr,:,:,:)

  end subroutine ctrl_map_gentim2d
#endif

end module ctrl_map_gentim_m 
