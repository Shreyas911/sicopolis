module ctrl_map_genarr

  implicit none

  public :: ctrl_map_genarr2d

contains

  subroutine ctrl_map_genarr2d

    use sico_types_m
    use sico_variables_m
        
    implicit none  

    integer(i4b)                 :: ctrl_index, j, i
    integer(i4b)                 :: iLow, iHigh, jLow, jHigh

    xx_genarr2d_vars    = XX_GENARR2D_VARS_ARR
    xx_genarr2d_bounds  = XX_GENARR2D_BOUNDS_ARR

    do ctrl_index = 1, NUM_CTRL_GENARR2D

      if (xx_genarr2d_bounds(ctrl_index) .NE. ' ') then
        read (unit=xx_genarr2d_bounds(ctrl_index),fmt=*) iLow, iHigh, jLow, jHigh
      else
                iLow  = 0
                iHigh = IMAX
                jLow  = 0
                jHigh = JMAX
      end if

      do j = 0, JMAX
        do i = 0, IMAX
          if ((i .GE. iLow) .AND. (i .LE. iHigh) .AND. (j .GE. jLow) .AND. (j .LE. jHigh)) then
            xx_genarr2d_mask(ctrl_index,j,i) = 1.0
          else
            xx_genarr2d_mask(ctrl_index,j,i) = 0.0
          end if
        end do
      end do 

      if(trim(xx_genarr2d_vars(ctrl_index)) .EQ. 'xx_c_slide') then
        
        c_slide = c_slide &
                + xx_genarr2d(ctrl_index,:,:) * xx_genarr2d_mask(ctrl_index,:,:)

      else if (trim(xx_genarr2d_vars(ctrl_index)) .EQ. 'xx_q_geo') then

        q_geo = q_geo + xx_genarr2d(ctrl_index,:,:) * xx_genarr2d_mask(ctrl_index,:,:)
 
      else 

        print *, "WARNING! ", xx_genarr2d_vars(ctrl_index), " is not in the genctrl setup yet!"

      endif

    end do
       
  end subroutine ctrl_map_genarr2d

end module ctrl_map_genarr 
