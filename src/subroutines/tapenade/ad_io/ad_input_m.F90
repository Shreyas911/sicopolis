module ad_input_m

    use sico_types_m
    use sico_variables_m_diff
    use error_m
  
    implicit none

#ifdef ALLOW_GENCTRL
    public :: ad_input
#endif
  
    contains

#ifdef ALLOW_GENCTRL
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

#ifdef DO_CTRL_GENARR2D
        real(dp), dimension(NUM_CTRL_GENARR2D,0:IMAX,0:JMAX) :: xx_genarr2d_conv
#ifdef ALLOW_TAP_TLM
        real(dp), dimension(NUM_CTRL_GENARR2D,0:IMAX,0:JMAX) :: xx_genarr2dd_conv
#endif /* ALLOW_TAP_TLM */
#endif

#ifdef DO_CTRL_GENARR3D
        real(dp), dimension(NUM_CTRL_GENARR3D,0:IMAX,0:JMAX,0:KCMAX) :: xx_genarr3d_conv
#ifdef ALLOW_TAP_TLM
        real(dp), dimension(NUM_CTRL_GENARR3D,0:IMAX,0:JMAX,0:KCMAX) :: xx_genarr3dd_conv
#endif /* ALLOW_TAP_TLM */
#endif

#ifdef DO_CTRL_GENTIM2D
        real(dp), dimension(NUM_CTRL_GENTIM2D,0:IMAX,0:JMAX,0:NTDAMAX) :: xx_gentim2d_conv
#ifdef ALLOW_TAP_TLM
        real(dp), dimension(NUM_CTRL_GENTIM2D,0:IMAX,0:JMAX,0:NTDAMAX) :: xx_gentim2dd_conv
#endif /* ALLOW_TAP_TLM */
#endif

        character(len=64), parameter :: thisroutine = 'ad_input'
        character(len=256) :: filename, filename_with_path, temp_path
        character(len= 16), parameter :: filename_extension = '.nc'
    
#if defined(AD_INPUT_PATH)

        xx_genarr2d_vars            = XX_GENARR2D_VARS_ARR
        xx_genarr3d_vars            = XX_GENARR3D_VARS_ARR
        xx_gentim2d_vars            = XX_GENTIM2D_VARS_ARR

        !-------- Create file name --------
        
        temp_path = AD_INPUT_PATH
#ifdef ALLOW_TAP_TLM
        filename = 'ad_input_tlm'//trim(filename_extension)
#endif /* ALLOW_TAP_TLM */
#ifdef ALLOW_TAP_ADJ
        filename = 'ad_input_adj'//trim(filename_extension)
#endif /* ALLOW_TAP_ADJ */
        filename_with_path = trim(temp_path)//'/'//trim(filename)
    
        !  ------ Open NetCDF file
        ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

        if (ios /= nf90_noerr) then
            errormsg = ' >>> '//trim(thisroutine)//': Error when opening a' &
            //               end_of_line &
            //'              NetCDF AD-input file!'
            call error(errormsg)
        end if

#ifdef DO_CTRL_GENARR2D
        do ctrl_index = 1, NUM_CTRL_GENARR2D
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_genarr2d_vars(ctrl_index))), ncv) )
            call check( nf90_get_var(ncid, ncv, xx_genarr2d_conv(ctrl_index,:,:)) )
#ifdef ALLOW_TAP_TLM
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_genarr2d_vars(ctrl_index)))//'d', ncv) )
            call check( nf90_get_var(ncid, ncv, xx_genarr2dd_conv(ctrl_index,:,:)) )
#endif /* ALLOW_TAP_TLM */
        end do
#endif

#ifdef DO_CTRL_GENARR3D
        do ctrl_index = 1, NUM_CTRL_GENARR3D
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_genarr3d_vars(ctrl_index))), ncv) )
            call check( nf90_get_var(ncid, ncv, xx_genarr3d_conv(ctrl_index,:,:,:)) )
#ifdef ALLOW_TAP_TLM
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_genarr3d_vars(ctrl_index)))//'d', ncv) )
            call check( nf90_get_var(ncid, ncv, xx_genarr3dd_conv(ctrl_index,:,:,:)) )
#endif /* ALLOW_TAP_TLM */
        end do
#endif

#ifdef DO_CTRL_GENTIM2D
        do ctrl_index = 1, NUM_CTRL_GENTIM2D
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_gentim2d_vars(ctrl_index))), ncv) )
            call check( nf90_get_var(ncid, ncv, xx_gentim2d_conv(ctrl_index,:,:,:)) )
#ifdef ALLOW_TAP_TLM
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_gentim2d_vars(ctrl_index)))//'d', ncv) )
            call check( nf90_get_var(ncid, ncv, xx_gentim2dd_conv(ctrl_index,:,:,:)) )
#endif /* ALLOW_TAP_TLM */
        end do
#endif

        !  ------ Close NetCDF file
        call check( nf90_close(ncid) )

#ifdef DO_CTRL_GENARR2D
        do i = 0, IMAX
        do j = 0, JMAX
        do ctrl_index = 1, NUM_CTRL_GENARR2D
            xx_genarr2d(ctrl_index,j,i) = xx_genarr2d_conv(ctrl_index,i,j)
#ifdef ALLOW_TAP_TLM
            xx_genarr2dd(ctrl_index,j,i) = xx_genarr2dd_conv(ctrl_index,i,j)
#endif /* ALLOW_TAP_TLM */
        end do
        end do
        end do
#endif

#ifdef DO_CTRL_GENARR3D
        do i = 0, IMAX
        do j = 0, JMAX
        do kc = 0, KCMAX
        do ctrl_index = 1, NUM_CTRL_GENARR3D
            xx_genarr3d(ctrl_index,kc,j,i) = xx_genarr3d_conv(ctrl_index,i,j,kc)
#ifdef ALLOW_TAP_TLM
            xx_genarr3dd(ctrl_index,kc,j,i) = xx_genarr3dd_conv(ctrl_index,i,j,kc)
#endif /* ALLOW_TAP_TLM */
        end do
        end do
        end do
        end do
#endif

#ifdef DO_CTRL_GENTIM2D
        do i = 0, IMAX
        do j = 0, JMAX
        do tad = 0, NTDAMAX
        do ctrl_index = 1, NUM_CTRL_GENTIM2D
            xx_gentim2d(ctrl_index,tad,j,i) = xx_gentim2d_conv(ctrl_index,i,j,tad)
#ifdef ALLOW_TAP_TLM
            xx_gentim2dd(ctrl_index,tad,j,i) = xx_gentim2dd_conv(ctrl_index,i,j,tad)
#endif /* ALLOW_TAP_TLM */
        end do
        end do
        end do
        end do
#endif

#endif

    end subroutine ad_input
#endif
  
end module ad_input_m
