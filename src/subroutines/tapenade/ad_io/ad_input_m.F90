module ad_input_m

    use sico_types_m
#if defined(ALLOW_NODIFF)
    use sico_variables_m
#else
    use sico_variables_m_diff
#endif
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

#if (defined(ALLOW_TAP_ADJ) && defined(ALLOW_TAP_ADJ_AT_ACTION))

#if (defined(AGE_COST) || defined(FAKE_AGE_COST))
        real(dp), dimension(0:IMAX,0:JMAX,0:KCMAX) :: age_cb_conv
#endif

#if (defined(BEDMACHINE_COST) || defined(FAKE_BEDMACHINE_COST))
        real(dp), dimension(0:IMAX,0:JMAX)         :: Hb_conv
#endif

#if (defined(ZS_COST) || defined(FAKE_ZS_COST))
        real(dp), dimension(0:IMAX,0:JMAX)         :: zsb_conv
#endif

#if (defined(ZL_COST) || defined(FAKE_ZL_COST))
        real(dp), dimension(0:IMAX,0:JMAX)         :: zlb_conv
#endif

#if (defined(SURFVEL_COST) || defined(FAKE_SURFVEL_COST))
#if !defined(SURF_VXVY_COST)
        real(dp), dimension(0:JMAX,0:IMAX)         :: vs
        real(dp), dimension(0:IMAX,0:JMAX)         :: vx_s_g_final_conv, vy_s_g_final_conv, vsb_conv
#else
        real(dp), dimension(0:IMAX,0:JMAX)         :: vx_s_gb_conv, vy_s_gb_conv
#endif
#endif

#endif /* ALLOW_TAP_ADJ && ALLOW_TAP_ADJ_AT_ACTION */

        character(len=64), parameter :: thisroutine = 'ad_input'
        character(len=256) :: filename, filename_with_path, temp_path
        character(len= 16), parameter :: filename_extension = '.nc'
    
#if defined(AD_INPUT_PATH)

#ifdef DO_CTRL_GENARR2D
        xx_genarr2d_vars            = XX_GENARR2D_VARS_ARR
#endif
#ifdef DO_CTRL_GENARR3D
        xx_genarr3d_vars            = XX_GENARR3D_VARS_ARR
#endif
#ifdef DO_CTRL_GENTIM2D
        xx_gentim2d_vars            = XX_GENTIM2D_VARS_ARR
#endif

        !-------- Create file name --------
        
        temp_path = AD_INPUT_PATH
#ifdef ALLOW_TAP_TLM
#ifdef ALLOW_TAP_TLM_A_ACTION
        filename = 'ad_input_tlm_hessaction'//trim(filename_extension)
#else
        filename = 'ad_input_tlm'//trim(filename_extension)
#endif
#endif /* ALLOW_TAP_TLM */
#ifdef ALLOW_TAP_ADJ
#ifdef ALLOW_TAP_ADJ_AT_ACTION
        filename = 'ad_input_adj_hessaction'//trim(filename_extension)
#else
        filename = 'ad_input_adj'//trim(filename_extension)
#endif
#endif /* ALLOW_TAP_ADJ */
#ifdef ALLOW_NODIFF
        filename = 'ad_input_nodiff'//trim(filename_extension)
#endif

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

#if (defined(ALLOW_TAP_ADJ) && defined(ALLOW_TAP_ADJ_AT_ACTION))

#if (defined(AGE_COST) || defined(FAKE_AGE_COST))
        call check( nf90_inq_varid(ncid, 'age_cb', ncv) )
        call check( nf90_get_var(ncid, ncv, age_cb_conv) )
#endif

#if (defined(BEDMACHINE_COST) || defined(FAKE_BEDMACHINE_COST))
        call check( nf90_inq_varid(ncid, 'Hb', ncv) )
        call check( nf90_get_var(ncid, ncv, Hb_conv) )
#endif

#if (defined(ZS_COST) || defined(FAKE_ZS_COST))
        call check( nf90_inq_varid(ncid, 'zsb', ncv) )
        call check( nf90_get_var(ncid, ncv, zsb_conv) )
#endif

#if (defined(ZL_COST) || defined(FAKE_ZL_COST))
        call check( nf90_inq_varid(ncid, 'zlb', ncv) )
        call check( nf90_get_var(ncid, ncv, zlb_conv) )
#endif

#if (defined(SURFVEL_COST) || defined(FAKE_SURFVEL_COST))
#if !defined(SURF_VXVY_COST)
        call check( nf90_inq_varid(ncid, 'vsb', ncv) )
        call check( nf90_get_var(ncid, ncv, vsb_conv) ) 
        call check( nf90_inq_varid(ncid, 'vx_s_g_final', ncv) )
        call check( nf90_get_var(ncid, ncv, vx_s_g_final_conv) ) 
        call check( nf90_inq_varid(ncid, 'vy_s_g_final', ncv) )
        call check( nf90_get_var(ncid, ncv, vy_s_g_final_conv) ) 
#else
        call check( nf90_inq_varid(ncid, 'vx_s_gb', ncv) )
        call check( nf90_get_var(ncid, ncv, vx_s_gb_conv) )
        call check( nf90_inq_varid(ncid, 'vy_s_gb', ncv) )
        call check( nf90_get_var(ncid, ncv, vy_s_gb_conv) )
#endif
#endif

#endif /* ALLOW_TAP_ADJ && ALLOW_TAP_ADJ_AT_ACTION */

        !  ------ Close NetCDF file
        call check( nf90_close(ncid) )

#ifdef DO_CTRL_GENARR2D
        do i = 0, IMAX
        do j = 0, JMAX
        do ctrl_index = 1, NUM_CTRL_GENARR2D
            xx_genarr2d(ctrl_index,j,i) = xx_genarr2d_conv(ctrl_index,i,j)
            xx_genarr2d_orig(ctrl_index,j,i) = xx_genarr2d_conv(ctrl_index,i,j)
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
            xx_genarr3d_orig(ctrl_index,kc,j,i) = xx_genarr3d_conv(ctrl_index,i,j,kc)
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
            xx_gentim2d_orig(ctrl_index,tad,j,i) = xx_gentim2d_conv(ctrl_index,i,j,tad)
#ifdef ALLOW_TAP_TLM
            xx_gentim2dd(ctrl_index,tad,j,i) = xx_gentim2dd_conv(ctrl_index,i,j,tad)
#endif /* ALLOW_TAP_TLM */
        end do
        end do
        end do
        end do
#endif

#if (defined(ALLOW_TAP_ADJ) && defined(ALLOW_TAP_ADJ_AT_ACTION))

#if (defined(AGE_COST) || defined(FAKE_AGE_COST))
        do i = 0, IMAX
        do j = 0, JMAX
        do kc = 0, KCMAX
            age_cb(kc,j,i) = age_cb_conv(i,j,kc)
        end do
        end do
        end do
#endif

#if (defined(BEDMACHINE_COST) || defined(FAKE_BEDMACHINE_COST))
        do i = 0, IMAX
        do j = 0, JMAX
            Hb(j,i) = Hb_conv(i,j)
        end do
        end do
#endif

#if (defined(ZS_COST) || defined(FAKE_ZS_COST))
        do i = 0, IMAX
        do j = 0, JMAX
            zsb(j,i) = zsb_conv(i,j)
        end do
        end do
#endif

#if (defined(ZL_COST) || defined(FAKE_ZL_COST))
        do i = 0, IMAX
        do j = 0, JMAX
            zlb(j,i) = zlb_conv(i,j)
        end do
        end do
#endif

#if (defined(SURFVEL_COST) || defined(FAKE_SURFVEL_COST))

#if (defined(YEAR_SEC))
year2sec = YEAR_SEC
#else
year2sec = 3.1556925445e+07_dp
#endif

        do i = 0, IMAX
        do j = 0, JMAX
#if !defined(SURF_VXVY_COST)
! /* ALLOW_TAPENADE: guarding against non-differentiable sqrt(0) */
          vs(j,i) = sqrt(vx_s_g_final_conv(i,j)**2 + vy_s_g_final_conv(i,j)**2)*year2sec
          if (vs(j,i) > 0) then
            vx_s_gb(j,i)  = vx_s_gb(j,i) + vsb_conv(i,j)*(vx_s_g_final_conv(i,j)/vs(j,i))*year2sec**2
            vy_s_gb(j,i)  = vy_s_gb(j,i) + vsb_conv(i,j)*(vy_s_g_final_conv(i,j)/vs(j,i))*year2sec**2
            vsb_conv(i,j) = 0.0
          else
            vx_s_gb(j,i)  = 0.0
            vy_s_gb(j,i)  = 0.0
          end if
#else
          vx_s_gb(j,i) = vx_s_gb_conv(i,j)
          vy_s_gb(j,i) = vy_s_gb_conv(i,j)
#endif
        end do
        end do

#endif

#endif /* ALLOW_TAP_ADJ && ALLOW_TAP_ADJ_AT_ACTION */

#if defined(DO_GENCTRL_PRIOR)

        filename = 'ad_input_nodiff_prior'//trim(filename_extension)
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
        call check( nf90_inq_varid(ncid, "genarr2d_gamma_arr", ncv) )
        call check( nf90_get_var(ncid, ncv, genarr2d_gamma_arr) )

        call check( nf90_inq_varid(ncid, "genarr2d_delta_arr", ncv) )
        call check( nf90_get_var(ncid, ncv, genarr2d_delta_arr) )

        call check( nf90_inq_varid(ncid, "genarr2d_sigma_arr", ncv) )
        call check( nf90_get_var(ncid, ncv, genarr2d_sigma_arr) )

        do ctrl_index = 1, NUM_CTRL_GENARR2D
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_genarr2d_vars(ctrl_index))), ncv) )
            call check( nf90_get_var(ncid, ncv, xx_genarr2d_conv(ctrl_index,:,:)) )
        end do
#endif

#ifdef DO_CTRL_GENARR3D
        call check( nf90_inq_varid(ncid, "genarr3d_gamma_arr", ncv) )
        call check( nf90_get_var(ncid, ncv, genarr3d_gamma_arr) )

        call check( nf90_inq_varid(ncid, "genarr3d_delta_arr", ncv) )
        call check( nf90_get_var(ncid, ncv, genarr3d_delta_arr) )

        call check( nf90_inq_varid(ncid, "genarr3d_sigma_arr", ncv) )
        call check( nf90_get_var(ncid, ncv, genarr3d_sigma_arr) )

        do ctrl_index = 1, NUM_CTRL_GENARR3D
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_genarr3d_vars(ctrl_index))), ncv) )
            call check( nf90_get_var(ncid, ncv, xx_genarr3d_conv(ctrl_index,:,:,:)) )
        end do
#endif

#ifdef DO_CTRL_GENTIM2D
        call check( nf90_inq_varid(ncid, "gentim2d_gamma_arr", ncv) )
        call check( nf90_get_var(ncid, ncv, gentim2d_gamma_arr) )

        call check( nf90_inq_varid(ncid, "gentim2d_delta_arr", ncv) )
        call check( nf90_get_var(ncid, ncv, gentim2d_delta_arr) )

        call check( nf90_inq_varid(ncid, "gentim2d_sigma_arr", ncv) )
        call check( nf90_get_var(ncid, ncv, gentim2d_sigma_arr) )

        do ctrl_index = 1, NUM_CTRL_GENTIM2D
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_gentim2d_vars(ctrl_index))), ncv) )
            call check( nf90_get_var(ncid, ncv, xx_gentim2d_conv(ctrl_index,:,:,:)) )
        end do
#endif

        !  ------ Close NetCDF file
        call check( nf90_close(ncid) )

#ifdef DO_CTRL_GENARR2D
        do i = 0, IMAX
        do j = 0, JMAX
        do ctrl_index = 1, NUM_CTRL_GENARR2D
            xx_genarr2d_prior(ctrl_index,j,i) = xx_genarr2d_conv(ctrl_index,i,j)
        end do
        end do
        end do
#endif

#ifdef DO_CTRL_GENARR3D
        do i = 0, IMAX
        do j = 0, JMAX
        do kc = 0, KCMAX
        do ctrl_index = 1, NUM_CTRL_GENARR3D
            xx_genarr3d_prior(ctrl_index,kc,j,i) = xx_genarr3d_conv(ctrl_index,i,j,kc)
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
            xx_gentim2d_prior(ctrl_index,tad,j,i) = xx_gentim2d_conv(ctrl_index,i,j,tad)
        end do
        end do
        end do
        end do
#endif

        filename = 'ad_input_nodiff_prior_X'//trim(filename_extension)
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
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_genarr2d_vars(ctrl_index)))//'d', ncv) )
            call check( nf90_get_var(ncid, ncv, xx_genarr2d_conv(ctrl_index,:,:)) )
        end do
#endif

#ifdef DO_CTRL_GENARR3D
        do ctrl_index = 1, NUM_CTRL_GENARR3D
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_genarr3d_vars(ctrl_index)))//'d', ncv) )
            call check( nf90_get_var(ncid, ncv, xx_genarr3d_conv(ctrl_index,:,:,:)) )
        end do
#endif

#ifdef DO_CTRL_GENTIM2D
        do ctrl_index = 1, NUM_CTRL_GENTIM2D
            call check( nf90_inq_varid(ncid, trim(adjustl(xx_gentim2d_vars(ctrl_index)))//'d', ncv) )
            call check( nf90_get_var(ncid, ncv, xx_gentim2d_conv(ctrl_index,:,:,:)) )
        end do
#endif

        !  ------ Close NetCDF file
        call check( nf90_close(ncid) )

#ifdef DO_CTRL_GENARR2D
        do i = 0, IMAX
        do j = 0, JMAX
        do ctrl_index = 1, NUM_CTRL_GENARR2D
            xx_genarr2d_prior_X(ctrl_index,j,i) = xx_genarr2d_conv(ctrl_index,i,j)
        end do
        end do
        end do
#endif

#ifdef DO_CTRL_GENARR3D
        do i = 0, IMAX
        do j = 0, JMAX
        do kc = 0, KCMAX
        do ctrl_index = 1, NUM_CTRL_GENARR3D
            xx_genarr3d_prior_X(ctrl_index,kc,j,i) = xx_genarr3d_conv(ctrl_index,i,j,kc)
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
            xx_gentim2d_prior_X(ctrl_index,tad,j,i) = xx_gentim2d_conv(ctrl_index,i,j,tad)
        end do
        end do
        end do
        end do
#endif

#endif

#endif

    end subroutine ad_input
#endif
  
end module ad_input_m
