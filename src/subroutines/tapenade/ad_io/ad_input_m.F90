module ad_input_m

    use sico_types_m
    use sico_variables_m_diff
    use sico_vars_m
    use error_m
  
    implicit none

    public :: ad_input
  
    contains
  
    subroutine ad_input
    
        use netcdf
        use nc_check_m
    
        implicit none
    
        integer(i4b) :: ncid, ncv
        !     ncid:      ID of the input file
        !     ncv:       Variable ID
        integer(i4b) :: i, j, kc, tad
        integer(i4b) :: ctrl_index
        integer(i4b) :: ios

        real(dp), dimension(NUM_CTRL_GENARR2D,0:IMAX,0:JMAX) :: xx_genarr2d_conv
        real(dp), dimension(NUM_CTRL_GENARR3D,0:IMAX,0:JMAX,0:KCMAX) :: xx_genarr3d_conv
        real(dp), dimension(NUM_CTRL_GENTIM2D,0:IMAX,0:JMAX,0:ADNMAX) :: xx_gentim2d_conv
    
        character(len=64), parameter :: thisroutine = 'ad_input'
        character(len=256) :: filename, filename_with_path, temp_path
        character(len= 16), parameter :: filename_extension = '.nc'
    
#if defined(AD_INPUT_PATH)

        xx_genarr2d_vars            = XX_GENARR2D_VARS_ARR
        xx_genarr3d_vars            = XX_GENARR3D_VARS_ARR
        xx_gentim2d_vars            = XX_GENTIM2D_VARS_ARR

        !-------- Create file name --------
        
        temp_path = AD_INPUT_PATH
        filename = 'ad_input'//trim(filename_extension)
        filename_with_path = trim(temp_path)//'/'//trim(filename)
    
        !  ------ Open NetCDF file
        ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

        if (ios /= nf90_noerr) then
            errormsg = ' >>> sico_init: Error when opening the file' &
                    //                 end_of_line &
                    //'                for the surface-temperature and SMB climatology!'
            call error(errormsg)
        end if
    
        do ctrl_index = 1, NUM_CTRL_GENARR2D
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_genarr2d_vars(ctrl_index))), ncv) )
            call check( nf90_get_var(ncid, ncv, xx_genarr2d_conv(ctrl_index,:,:)) )
        end do

        do ctrl_index = 1, NUM_CTRL_GENARR3D
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_genarr3d_vars(ctrl_index))), ncv) )
            call check( nf90_get_var(ncid, ncv, xx_genarr3d_conv(ctrl_index,:,:,:)) )
        end do

        do ctrl_index = 1, NUM_CTRL_GENTIM2D
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_gentim2d_vars(ctrl_index))), ncv) )
            call check( nf90_get_var(ncid, ncv, xx_gentim2d_conv(ctrl_index,:,:,:)) )
        end do

        call check( nf90_close(ncid) )

        do i = 0, IMAX
        do j = 0, JMAX
        do ctrl_index = 1, NUM_CTRL_GENARR2D
            xx_genarr2d(ctrl_index,j,i) = xx_genarr2d_conv(ctrl_index,i,j)
        end do
        end do
        end do
    
        do i = 0, IMAX
        do j = 0, JMAX
        do kc = 0, KCMAX
        do ctrl_index = 1, NUM_CTRL_GENARR3D
            xx_genarr3d(ctrl_index,kc,j,i) = xx_genarr3d_conv(ctrl_index,i,j,kc)
        end do
        end do
        end do
        end do
    
        do i = 0, IMAX
        do j = 0, JMAX
        do tad = 0, ADNMAX
        do ctrl_index = 1, NUM_CTRL_GENTIM2D
            xx_gentim2d(ctrl_index,tad,j,i) = xx_gentim2d_conv(ctrl_index,i,j,tad)
        end do
        end do
        end do
        end do

#endif

    end subroutine ad_input
  
end module ad_input_m
  