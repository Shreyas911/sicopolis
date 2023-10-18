module ctrl_map_gentim_m

  use sico_types_m
  use sico_variables_m
  use error_m
  use globals

  implicit none

  public :: ctrl_map_ini_gentim2d, ctrl_map_gentim2d

contains

  subroutine ctrl_map_ini_gentim2d(time_init, dtime, itercount)

    implicit none

    real(dp), intent(in)        :: time_init, dtime
    integer(i4b), intent(in)    :: itercount
    integer(i4b)                :: ctrl_index
    integer(i4b)                :: igen_temp

    xx_gentim2d_vars            = XX_GENTIM2D_VARS_ARR
    xx_gentim2d_period          = XX_GENTIM2D_PERIOD

    igen_temp = 0

    do ctrl_index = 1, NUM_CTRL_GENTIM2D
      
      if (trim(xx_gentim2d_vars(ctrl_index)) .EQ. 'xx_temp') then
        igen_temp = ctrl_index
      else
        errormsg = ' >>> ctrl_map_ini_gentim2d: ' &
          //"This control variable is not in the gentim2d setup yet!"
        call error(errormsg)
      end if
  
    end do

    if (igen_temp .GT. 0) then
      call ctrl_map_gentim2d(temp_ma, igen_temp, time_init, dtime, itercount)
      call ctrl_map_gentim2d(temp_mm(:,:,7), igen_temp, time_init, dtime, itercount)
    end if
    
  end subroutine ctrl_map_ini_gentim2d

  subroutine ctrl_map_gentim2d(fld, iarr, time_init, dtime, itercount)

    implicit none  

    real(dp), intent(in)                        :: time_init, dtime
    integer(i4b), intent(in)                    :: itercount
    real(dp), dimension(0:JMAX,0:IMAX)          :: fld, xx_gen, xx_gen_mask
    integer(i4b)                                :: iarr, j, i, k2
    integer(i4b)                                :: iLow, iHigh, jLow, jHigh
    integer(i4b)                                :: adLow, adHigh
    real(dp)                                    :: time, alpha

    time = (time_init + real(itercount,dp)*dtime)*sec2year

    xx_gentim2d_mask = 1.0
   
    adLow  = floor(time/xx_gentim2d_period)
    adHigh = ceiling(time/xx_gentim2d_period)
     
    alpha  = (time - adLow*xx_gentim2d_period)/xx_gentim2d_period

    xx_gen = xx_gentim2d(iarr,adLow,:,:)  * (1-alpha) &
           + xx_gentim2d(iarr,adHigh,:,:) * alpha
    xx_gen_mask = xx_gentim2d_mask(iarr,adLow,:,:)  * (1-alpha) &
                + xx_gentim2d_mask(iarr,adHigh,:,:) * alpha

    fld = fld + xx_gen * xx_gen_mask

  end subroutine ctrl_map_gentim2d

end module ctrl_map_gentim_m 
