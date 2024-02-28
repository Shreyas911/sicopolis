#if defined(ALLOW_GENCTRL)
module ad_output_m

  use sico_types_m
  use sico_variables_m_diff
  use sico_vars_m
  use error_m

  implicit none

  private :: set_cmode
  public :: ad_output

  contains

  subroutine ad_output

    use netcdf
    use nc_check_m
    use sico_variables_m_diff

    implicit none

    integer(i4b) :: i, j, kc, tad
    integer(i4b) :: ctrl_index
    integer(i4b) :: ios
    integer(i4b) :: cmode
    integer(i4b) :: n_deflate_level
    logical      :: flag_shuffle
    integer(i4b) :: ncid, ncv
    !     ncid:      ID of the output file
    !     ncv:       Variable ID
    integer(i4b) :: ncd, nc1d, nc2d(2), nc3d(3)
    !     ncd:       Dimension ID
    !     nc1d:      Dimension of a 1-d array
    !     nc2d:      Vector with the dimensions of a 2-d array
    !     nc3d:      Vector with the dimensions of a 3-d array

    integer(i4b) :: nc1cor_i(1), nc1cor_j(1), &
                    nc1cor_kc(1), nc1cor_tad(1), &
                    nc2cor_ij(2), nc2cor_ijtad(3), &
                    nc3cor_ijkc(3)
    !     nc1cor(1): Corner of a 1-d array
    !     nc2cor(2): Corner of a 2-d array
    !     nc3cor(3): Corner of a 3-d array
    integer(i4b) :: nc1cnt_i(1), nc1cnt_j(1), &
                    nc1cnt_kc(1), nc1cnt_tad(1), &
                    nc2cnt_ij(2), nc2cnt_ijtad(3), &
                    nc3cnt_ijkc(3)
    !     nc1cnt(1): Count of a 1-d array
    !     nc2cnt(2): Count of a 2-d array
    !     nc3cnt(3): Count of a 3-d array

#ifdef TAP_GENCTRL_TLM
    integer(i4b) :: nc0cor_fcd(1)
    integer(i4b) :: nc0cnt_fcd(1)
    !     nc0cor_fcd: Defined specially for fcd
    !     nc0cnt_fcd: Defined specially for fcd
#endif

    real(dp) :: xi_conv(0:IMAX), eta_conv(0:JMAX), &
    sigma_level_c_conv(0:KCMAX), time_ad_conv(0:ADNMAX)
#ifdef TAP_GENCTRL_TLM
    real(dp), dimension(1) :: fcd_arr
#endif

    real(dp), dimension(NUM_CTRL_GENARR2D,0:IMAX,0:JMAX) :: xx_genarr2d_conv
    real(dp), dimension(NUM_CTRL_GENARR3D,0:IMAX,0:JMAX,0:KCMAX) :: xx_genarr3d_conv
    real(dp), dimension(NUM_CTRL_GENTIM2D,0:IMAX,0:JMAX,0:ADNMAX) :: xx_gentim2d_conv

#ifdef TAP_GENCTRL_ADJ
    real(dp), dimension(NUM_CTRL_GENARR2D,0:IMAX,0:JMAX) :: xx_genarr2db_conv
    real(dp), dimension(NUM_CTRL_GENARR2D,0:IMAX,0:JMAX,0:KCMAX) :: xx_genarr3db_conv
    real(dp), dimension(NUM_CTRL_GENARR2D,0:IMAX,0:JMAX,0:ADNMAX) :: xx_gentim2db_conv
#endif

    character(len=64), parameter :: thisroutine = 'ad_output'

    character(len=256) :: filename, filename_with_path, temp_path
    character(len= 16) :: ch_date, ch_time, ch_zone

    character(len=256) :: buffer
    character(len= 16), parameter :: filename_extension = '.nc'
    character(len= 16), allocatable :: coord_id(:)

#ifdef TAP_GENCTRL_TLM
    fcd_arr(1) = fcd
    nc0cor_fcd = (/ 1 /)
    nc0cnt_fcd = (/ 1 /)
#endif

    nc1cor_i = (/ 1 /)
    nc1cnt_i = (/ IMAX+1 /)
    
    nc1cor_j = (/ 1 /)
    nc1cnt_j = (/ JMAX+1 /)
    
    nc1cor_kc = (/ 1 /)
    nc1cnt_kc = (/ KCMAX+1 /)

    nc1cor_tad = (/ 1 /)
    nc1cnt_tad = (/ ADNMAX+1 /)
    
    nc2cor_ij = (/ 1, 1 /)
    nc2cnt_ij = (/ IMAX+1, JMAX+1 /)
    
    nc3cor_ijkc = (/ 1, 1, 1 /)
    nc3cnt_ijkc = (/ IMAX+1, JMAX+1, KCMAX+1 /)

    nc2cor_ijtad = (/ 1, 1, 1 /)
    nc2cnt_ijtad = (/ IMAX+1, JMAX+1, ADNMAX+1 /)

    !-------- Create file name --------
    
    temp_path = AD_OUTPUT_PATH
    filename = 'ad_output'//trim(filename_extension)
    filename_with_path = trim(temp_path)//'/'//trim(filename)

    !-------- File initialization --------

    if (allocated(coord_id)) deallocate(coord_id); allocate(coord_id(5))
    coord_id(1) = 'x'; coord_id(2) = 'y'
    coord_id(3) = 'zeta_c'; coord_id(4) = 'time_ad'
    coord_id(5) = 'fcd_dummy_dim'
    
    !  ------ Open NetCDF file
    
    call set_cmode(cmode, n_deflate_level, flag_shuffle)

    ios = nf90_create(trim(filename_with_path), cmode, ncid)
    if (ios /= nf90_noerr) then
        errormsg = ' >>> '//trim(thisroutine)//': Error when opening a' &
                //               end_of_line &
                //'              NetCDF AD-output file!'
        call error(errormsg)
    end if
    
    !  ------ Global attributes
    
    call date_and_time(ch_date, ch_time, ch_zone)
    buffer = ch_date(1:4)//'-'//ch_date(5:6)//'-'//ch_date(7:8)//' '// &
              ch_time(1:2)//':'//ch_time(3:4)//':'//ch_time(5:6)//' '// &
              ch_zone(1:3)//':'//ch_zone(4:5)//' - Data produced'
    call check( nf90_put_att(ncid, NF90_GLOBAL, 'history', trim(buffer)), &
                thisroutine )

    !  ------ Definition of the dimensions

    call check( nf90_def_dim(ncid, trim(coord_id(1)), IMAX+1,  ncd), thisroutine )
    call check( nf90_def_dim(ncid, trim(coord_id(2)), JMAX+1,  ncd), thisroutine )
    call check( nf90_def_dim(ncid, trim(coord_id(3)), KCMAX+1, ncd), thisroutine )
    call check( nf90_def_dim(ncid, trim(coord_id(4)), ADNMAX+1, ncd), thisroutine )
    call check( nf90_def_dim(ncid, trim(coord_id(5)), 1, ncd), thisroutine )

    !    ---- x (= xi)

    call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc1d), &
          thisroutine )

    call check( nf90_def_var(ncid, 'x', NF90_DOUBLE, nc1d, ncv), &
          thisroutine )

    buffer = 'm'
    call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
          thisroutine )
    buffer = 'projection_x_coordinate'
    call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
          thisroutine )
    buffer = 'x-coordinate of the grid point i'
    call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
          thisroutine )
    call check( nf90_put_att(ncid, ncv, 'axis', 'x'), &
          thisroutine )

    !    ---- y (= eta)

    call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc1d), &
          thisroutine )

    call check( nf90_def_var(ncid, 'y', NF90_DOUBLE, nc1d, ncv), &
          thisroutine )

    buffer = 'm'
    call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
          thisroutine )
    buffer = 'projection_y_coordinate'
    call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
          thisroutine )
    buffer = 'y-coordinate of the grid point j'
    call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
          thisroutine )
    call check( nf90_put_att(ncid, ncv, 'axis', 'y'), &
          thisroutine )

    !    ---- sigma_level_c

    call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc1d), &
          thisroutine )

    call check( nf90_def_var(ncid, 'sigma_level_c', NF90_DOUBLE, nc1d, ncv), &
          thisroutine )

    buffer = 'up'
    call check( nf90_put_att(ncid, ncv, 'positive', trim(buffer)), &
          thisroutine )
    buffer = 'land_ice_kc_layer_sigma_coordinate'
    call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
          thisroutine )
    buffer = 'sigma-coordinate of the grid point kc'
    call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
          thisroutine )

    !    ---- time_ad

    call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc1d), &
          thisroutine )

    call check( nf90_def_var(ncid, 'time_ad', NF90_DOUBLE, nc1d, ncv), &
          thisroutine )

    buffer = 'time_ad_gentim2d'
    call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
          thisroutine )
    buffer = 'Times between which linear interp. for gentim2d'
    call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
          thisroutine )

#ifdef TAP_GENCTRL_TLM
    !    ---- Define fcd
    call check( nf90_inq_dimid(ncid, trim(coord_id(5)), nc1d), &
                      thisroutine )
      
#if (NETCDF4_ENABLED==1)
    call check( nf90_def_var(ncid, 'fcd', &
                      NF90_DOUBLE, nc1d, ncv, &
                      deflate_level=n_deflate_level, shuffle=flag_shuffle), &
                      thisroutine )
#else
    call check( nf90_def_var(ncid, 'fcd', &
                      NF90_DOUBLE, nc1d, ncv), &
                      thisroutine )     
#endif
#endif     

    !    ---- Define xx_genarr2d variables
    do ctrl_index = 1, NUM_CTRL_GENARR2D

      call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
                thisroutine )
      call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
                thisroutine )

#if (NETCDF4_ENABLED==1)
      call check( nf90_def_var(ncid, trim(adjustl(xx_genarr2d_vars(ctrl_index))), &
                NF90_DOUBLE, nc2d, ncv, &
                deflate_level=n_deflate_level, shuffle=flag_shuffle), &
                thisroutine )
#else
      call check( nf90_def_var(ncid, trim(adjustl(xx_genarr2d_vars(ctrl_index))), &
                NF90_DOUBLE, nc2d, ncv), &
                thisroutine )     
#endif      

      call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
                thisroutine )
      call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
                thisroutine )

#if (NETCDF4_ENABLED==1)
      call check( nf90_def_var(ncid, trim(adjustl(xx_genarr2d_vars(ctrl_index)))//'b', &
                NF90_DOUBLE, nc2d, ncv, &
                deflate_level=n_deflate_level, shuffle=flag_shuffle), &
                thisroutine )
#else
      call check( nf90_def_var(ncid, trim(adjustl(xx_genarr2d_vars(ctrl_index)))//'b', &
                NF90_DOUBLE, nc2d, ncv), &
                thisroutine )
#endif
    end do

    !    ---- Define xx_genarr3d variables
    do ctrl_index = 1, NUM_CTRL_GENARR3D

      call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
                thisroutine )
      call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
                thisroutine )
      call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)), &
                thisroutine )

#if (NETCDF4_ENABLED==1)
      call check( nf90_def_var(ncid, trim(adjustl(xx_genarr3d_vars(ctrl_index))), &
                NF90_DOUBLE, nc3d, ncv, &
                deflate_level=n_deflate_level, shuffle=flag_shuffle), &
                thisroutine )
#else
      call check( nf90_def_var(ncid, trim(adjustl(xx_genarr3d_vars(ctrl_index))), &
                NF90_DOUBLE, nc3d, ncv), &
                thisroutine )     
#endif      

      call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
                thisroutine )
      call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
                thisroutine )
      call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)), &
                thisroutine )

#if (NETCDF4_ENABLED==1)
      call check( nf90_def_var(ncid, trim(adjustl(xx_genarr3d_vars(ctrl_index)))//'b', &
                NF90_DOUBLE, nc3d, ncv, &
                deflate_level=n_deflate_level, shuffle=flag_shuffle), &
                thisroutine )
#else
      call check( nf90_def_var(ncid, trim(adjustl(xx_genarr3d_vars(ctrl_index)))//'b', &
                NF90_DOUBLE, nc3d, ncv), &
                thisroutine )
#endif
    end do

    !    ---- Define xx_gentim2d variables
    do ctrl_index = 1, NUM_CTRL_GENTIM2D
      
      call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
                thisroutine )
      call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
                thisroutine )
      call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)), &
                thisroutine )

#if (NETCDF4_ENABLED==1)
      call check( nf90_def_var(ncid, trim(adjustl(xx_gentim2d_vars(ctrl_index))), &
                NF90_DOUBLE, nc3d, ncv, &
                deflate_level=n_deflate_level, shuffle=flag_shuffle), &
                thisroutine )
#else
      call check( nf90_def_var(ncid, trim(adjustl(xx_gentim2d_vars(ctrl_index))), &
                NF90_DOUBLE, nc3d, ncv), &
                thisroutine )     
#endif      

      call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
                thisroutine )
      call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
                thisroutine )
      call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)), &
                thisroutine )

#if (NETCDF4_ENABLED==1)
      call check( nf90_def_var(ncid, trim(adjustl(xx_gentim2d_vars(ctrl_index)))//'b', &
                NF90_DOUBLE, nc3d, ncv, &
                deflate_level=n_deflate_level, shuffle=flag_shuffle), &
                thisroutine )
#else
      call check( nf90_def_var(ncid, trim(adjustl(xx_gentim2d_vars(ctrl_index)))//'b', &
                NF90_DOUBLE, nc3d, ncv), &
                thisroutine )
#endif     

    end do
    
!    ---- End of the definitions

    call check( nf90_enddef(ncid), thisroutine )

!    ---- Put variables in

    do i=0, IMAX
      xi_conv(i) = xi(i)
    end do
   
    do j=0, JMAX
      eta_conv(j) = eta(j)
    end do
   
    do kc=0, KCMAX
      sigma_level_c_conv(kc) = eaz_c_quotient(kc)
    end do

    do tad=0, ADNMAX
      time_ad_conv(tad) = tad*XX_GENTIM2D_PERIOD
    end do

    call check( nf90_inq_varid(ncid, 'x', ncv), thisroutine )
    call check( nf90_put_var(ncid, ncv, xi_conv, &
                             start=nc1cor_i, count=nc1cnt_i), &
                thisroutine )
    
    call check( nf90_inq_varid(ncid, 'y', ncv), thisroutine )
    call check( nf90_put_var(ncid, ncv, eta_conv, &
                             start=nc1cor_j, count=nc1cnt_j), &
                thisroutine )
    
    call check( nf90_inq_varid(ncid, 'sigma_level_c', ncv), thisroutine )
    call check( nf90_put_var(ncid, ncv, sigma_level_c_conv, &
                             start=nc1cor_kc, count=nc1cnt_kc), &
                thisroutine )

    call check( nf90_inq_varid(ncid, 'time_ad', ncv), thisroutine )
    call check( nf90_put_var(ncid, ncv, time_ad_conv, &
                             start=nc1cor_tad, count=nc1cnt_tad), &
                thisroutine )

#ifdef TAP_GENCTRL_TLM
    call check( nf90_inq_varid(ncid, 'fcd', &
                ncv), &
                thisroutine )
    call check( nf90_put_var(ncid, ncv, fcd_arr, &
                             start=nc0cor_fcd, count=nc0cnt_fcd), &
                thisroutine )
#endif

    do i=0, IMAX
    do j=0, JMAX
    do ctrl_index = 1, NUM_CTRL_GENARR2D
      xx_genarr2d_conv(ctrl_index,i,j) = xx_genarr2d(ctrl_index,j,i)
#ifdef TAP_GENCTRL_ADJ
      xx_genarr2db_conv(ctrl_index,i,j) = xx_genarr2db(ctrl_index,j,i)
#endif
    end do
    end do
    end do

    do i=0, IMAX
    do j=0, JMAX
    do kc=0, KCMAX
    do ctrl_index = 1, NUM_CTRL_GENARR3D
      xx_genarr3d_conv(ctrl_index,i,j,kc) = xx_genarr3d(ctrl_index,kc,j,i)
#ifdef TAP_GENCTRL_ADJ
      xx_genarr3db_conv(ctrl_index,i,j,kc) = xx_genarr3db(ctrl_index,kc,j,i)
#endif
    end do
    end do
    end do
    end do

    do i=0, IMAX
    do j=0, JMAX
    do tad=0, ADNMAX
    do ctrl_index = 1, NUM_CTRL_GENTIM2D
      xx_gentim2d_conv(ctrl_index,i,j,tad) = xx_gentim2d(ctrl_index,tad,j,i)
#ifdef TAP_GENCTRL_ADJ
      xx_gentim2db_conv(ctrl_index,i,j,tad) = xx_gentim2db(ctrl_index,tad,j,i)
#endif
    end do
    end do
    end do
    end do

    do ctrl_index = 1, NUM_CTRL_GENARR2D

      call check( nf90_inq_varid(ncid, trim(adjustl(xx_genarr2d_vars(ctrl_index))), &
                  ncv), &
                  thisroutine )
      call check( nf90_put_var(ncid, ncv, xx_genarr2d_conv(ctrl_index,:,:), &
                               start=nc2cor_ij, count=nc2cnt_ij), &
                  thisroutine )
#ifdef TAP_GENCTRL_ADJ
      call check( nf90_inq_varid(ncid, trim(adjustl(xx_genarr2d_vars(ctrl_index)))//'b', &
                  ncv), &
                  thisroutine )
      call check( nf90_put_var(ncid, ncv, xx_genarr2db_conv(ctrl_index,:,:), &
                               start=nc2cor_ij, count=nc2cnt_ij), &
                  thisroutine )
#endif
    end do

    do ctrl_index = 1, NUM_CTRL_GENARR3D

      call check( nf90_inq_varid(ncid, trim(adjustl(xx_genarr3d_vars(ctrl_index))), &
                  ncv), &
                  thisroutine )
      call check( nf90_put_var(ncid, ncv, xx_genarr3d_conv(ctrl_index,:,:,:), &
                               start=nc3cor_ijkc, count=nc3cnt_ijkc), &
                  thisroutine )
#ifdef TAP_GENCTRL_ADJ
      call check( nf90_inq_varid(ncid, trim(adjustl(xx_genarr3d_vars(ctrl_index)))//'b', &
                  ncv), &
                  thisroutine )
      call check( nf90_put_var(ncid, ncv, xx_genarr3db_conv(ctrl_index,:,:,:), &
                               start=nc3cor_ijkc, count=nc3cnt_ijkc), &
                  thisroutine )
#endif
    end do

    do ctrl_index = 1, NUM_CTRL_GENTIM2D

      call check( nf90_inq_varid(ncid, trim(adjustl(xx_gentim2d_vars(ctrl_index))), &
                  ncv), &
                  thisroutine )
      call check( nf90_put_var(ncid, ncv, xx_gentim2d_conv(ctrl_index,:,:,:), &
                               start=nc2cor_ijtad, count=nc2cnt_ijtad), &
                  thisroutine )
#ifdef TAP_GENCTRL_ADJ
      call check( nf90_inq_varid(ncid, trim(adjustl(xx_gentim2d_vars(ctrl_index)))//'b', &
                  ncv), &
                  thisroutine )
      call check( nf90_put_var(ncid, ncv, xx_gentim2db_conv(ctrl_index,:,:,:), &
                               start=nc2cor_ijtad, count=nc2cnt_ijtad), &
                  thisroutine )
#endif
    end do

    call check( nf90_sync(ncid),  thisroutine )
    call check( nf90_close(ncid), thisroutine )

    deallocate(coord_id)

  end subroutine ad_output

!-------------------------------------------------------------------------------
!> Set the creation mode and compression type for NetCDF files.
!<------------------------------------------------------------------------------
  subroutine set_cmode(cmode, n_deflate_level, flag_shuffle)

    use netcdf
  
    implicit none
  
    integer(i4b), intent(out) :: cmode
    integer(i4b), intent(out) :: n_deflate_level
    logical,      intent(out) :: flag_shuffle
  
#if (NETCDF4_ENABLED==1)
  
    cmode           = ior(NF90_NETCDF4, NF90_CLASSIC_MODEL)
    n_deflate_level = 1
    flag_shuffle    = .true.
  
#else
  
    cmode           = NF90_NOCLOBBER
    n_deflate_level = 0
    flag_shuffle    = .false.
  
#endif
  
    end subroutine set_cmode

end module ad_output_m
#endif