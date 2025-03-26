!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module  :  c o s t _ i o _ m
!
!! Declarations of control variables for adjointing.
!!
!!##### Authors
!!
!! Liz Curry-Logan, Shreyas Sunil Gaikwad,
!! Sri Hari Krishna Narayanan
!!
!!##### License
!!
!! This file is part of SICOPOLIS.
!!
!! SICOPOLIS is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! SICOPOLIS is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with SICOPOLIS. If not, see <https://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------------------------------------------------------------
!> Declarations of control variables for adjointing.
!<------------------------------------------------------------------------------
module cost_io_m

    use sico_types_m  
    use sico_variables_m
#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
    use sico_vars_m
#endif
    use error_m
  
    implicit none
  
    public :: read_cost_data
  
    contains

subroutine read_cost_data()
    
    use netcdf
    use nc_check_m
    
    implicit none
    
    integer(i4b) :: ncid, ncv
    !     ncid:      ID of the input file
    !     ncv:       Variable ID
    integer(i4b) :: i, j, kc, tad, KDATA
    integer(i4b) :: ctrl_index
    integer(i4b) :: ios

    character(len=64), parameter :: thisroutine = 'read_cost_data'
    character(len=256) :: filename, filename_with_path, temp_path
    character(len= 16), parameter :: filename_extension = '.nc'

#if (defined(BEDMACHINE_COST) || defined(FAKE_BEDMACHINE_COST) || defined(AGE_COST) || defined(FAKE_AGE_COST))
    real(dp), dimension(0:IMAX,0:JMAX) :: H_BedMachine_data_conv
    real(dp), dimension(0:IMAX,0:JMAX) :: H_unc_BedMachine_data_conv
#endif
#if (defined(BEDMACHINE_COST) || defined(FAKE_BEDMACHINE_COST) || defined(ZS_COST) || defined(FAKE_ZS_COST) || defined(ZL_COST) || defined(FAKE_ZL_COST) || defined(SURFVEL_COST) || defined(FAKE_SURFVEL_COST))
    real(dp), dimension(0:IMAX,0:JMAX) :: zs_BedMachine_data_conv
    real(dp), dimension(0:IMAX,0:JMAX) :: zs_unc_BedMachine_data_conv
#endif
#if (defined(ZL_COST) || defined(FAKE_ZL_COST))
    real(dp), dimension(0:IMAX,0:JMAX) :: zl_BedMachine_data_conv
    real(dp), dimension(0:IMAX,0:JMAX) :: zl_unc_BedMachine_data_conv
#endif
#if (defined(AGE_COST) || defined(FAKE_AGE_COST))
    real(dp), dimension(0:IMAX,0:JMAX,0:KCMAX) :: age_data_conv
    real(dp), dimension(0:IMAX,0:JMAX,0:KCMAX) :: age_unc_data_conv
#endif
#if (defined(SURFVEL_COST) || defined(FAKE_SURFVEL_COST))
#if !defined(SURF_VXVY_COST)
    real(dp), dimension(0:IMAX,0:JMAX) :: vs_MEaSUREs_data_conv
    real(dp), dimension(0:IMAX,0:JMAX) :: vs_unc_MEaSUREs_data_conv
#else
    real(dp), dimension(0:IMAX,0:JMAX) :: vx_MEaSUREs_data_conv, vy_MEaSUREs_data_conv
    real(dp), dimension(0:IMAX,0:JMAX) :: vx_unc_MEaSUREs_data_conv, vy_unc_MEaSUREs_data_conv
#endif
#endif

    !-------- Create file name --------
    
#ifdef COST_INPUT_PATH
    temp_path = COST_INPUT_PATH
#endif

#if (defined(BEDMACHINE_COST) || defined(AGE_COST))
#if (IMAX==168)
    filename = 'bm5_data_10kms'//trim(filename_extension)
#elif (IMAX==42)
    filename = 'bm5_data_40kms'//trim(filename_extension)
#elif (IMAX==105)
    filename = 'bm5_data_16kms'//trim(filename_extension)
#else
    errormsg = ' >>> '//trim(thisroutine)//': Error when looking for a' &
    //               end_of_line &
    //'              BedMachine data file!'
    call error(errormsg)
#endif
    filename_with_path = trim(temp_path)//'/'//trim(filename)

    !  ------ Open NetCDF file
    ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

    if (ios /= nf90_noerr) then
        errormsg = ' >>> '//trim(thisroutine)//': Error when opening a' &
        //               end_of_line &
        //'              NetCDF BedMachine data file!'
        call error(errormsg)
    end if

    call check( nf90_inq_varid(ncid, 'H', ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, H_BedMachine_data_conv), thisroutine )
#ifdef ALLOW_BEDMACHINE_UNCERT
    call check( nf90_inq_varid(ncid, BEDMACHINE_UNCERT_FIELD, ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, H_unc_BedMachine_data_conv), thisroutine )
#endif
    !  ------ Close NetCDF file
    call check( nf90_close(ncid) )

    do i = 0, IMAX
        do j = 0, JMAX
            H_BedMachine_data(j,i) = H_BedMachine_data_conv(i,j)
#ifdef ALLOW_BEDMACHINE_UNCERT
            H_unc_BedMachine_data(j,i) = H_unc_BedMachine_data_conv(i,j)
#endif
        end do
    end do
#endif

#if (defined(FAKE_BEDMACHINE_COST) || defined(FAKE_AGE_COST))
#if (IMAX==168)
    filename = 'fake_bm5_data_10kms'//trim(filename_extension)
#elif (IMAX==42)
    filename = 'fake_bm5_data_40kms'//trim(filename_extension)
#elif (IMAX==105)
    filename = 'fake_bm5_data_16kms'//trim(filename_extension)
#else
    errormsg = ' >>> '//trim(thisroutine)//': Error when looking for a' &
    //               end_of_line &
    //'              Fake BedMachine data file!'
    call error(errormsg)
#endif
    filename_with_path = trim(temp_path)//'/'//trim(filename)

    !  ------ Open NetCDF file
    ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

    if (ios /= nf90_noerr) then
        errormsg = ' >>> '//trim(thisroutine)//': Error when opening a' &
        //               end_of_line &
        //'              Fake NetCDF BedMachine data file!'
        call error(errormsg)
    end if

    call check( nf90_inq_varid(ncid, 'H', ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, H_BedMachine_data_conv), thisroutine )
#ifdef ALLOW_BEDMACHINE_UNCERT
    call check( nf90_inq_varid(ncid, BEDMACHINE_UNCERT_FIELD, ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, H_unc_BedMachine_data_conv), thisroutine )
#endif
    !  ------ Close NetCDF file
    call check( nf90_close(ncid) )

    do i = 0, IMAX
        do j = 0, JMAX
            H_BedMachine_data(j,i) = H_BedMachine_data_conv(i,j)
#ifdef ALLOW_BEDMACHINE_UNCERT
            H_unc_BedMachine_data(j,i) = H_unc_BedMachine_data_conv(i,j)
#endif
        end do
    end do
#endif

#if (defined(BEDMACHINE_COST) || defined(ZS_COST) || defined(ZL_COST) || defined(SURFVEL_COST))
#if (IMAX==168)
    filename = 'bm5_data_10kms'//trim(filename_extension)
#elif (IMAX==42)
    filename = 'bm5_data_40kms'//trim(filename_extension)
#elif (IMAX==105)
    filename = 'bm5_data_16kms'//trim(filename_extension)
#else
    errormsg = ' >>> '//trim(thisroutine)//': Error when looking for a' &
    //               end_of_line &
    //'              BedMachine data file!'
    call error(errormsg)
#endif
    filename_with_path = trim(temp_path)//'/'//trim(filename)

    !  ------ Open NetCDF file
    ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

    if (ios /= nf90_noerr) then
        errormsg = ' >>> '//trim(thisroutine)//': Error when opening a' &
        //               end_of_line &
        //'              NetCDF BedMachine data file!'
        call error(errormsg)
    end if

    call check( nf90_inq_varid(ncid, 'zs', ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, zs_BedMachine_data_conv), thisroutine )
#ifdef ALLOW_ZS_UNCERT
    call check( nf90_inq_varid(ncid, ZS_UNCERT_FIELD, ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, zs_unc_BedMachine_data_conv), thisroutine )
#endif
    !  ------ Close NetCDF file
    call check( nf90_close(ncid) )

    do i = 0, IMAX
        do j = 0, JMAX
            zs_BedMachine_data(j,i) = zs_BedMachine_data_conv(i,j)
#ifdef ALLOW_ZS_UNCERT
            zs_unc_BedMachine_data(j,i) = zs_unc_BedMachine_data_conv(i,j)
#endif
        end do
    end do
#endif

#if (defined(FAKE_BEDMACHINE_COST) || defined(FAKE_ZS_COST) || defined(FAKE_ZL_COST) || defined(FAKE_SURFVEL_COST))
#if (IMAX==168)
    filename = 'fake_bm5_data_10kms'//trim(filename_extension)
#elif (IMAX==42)
    filename = 'fake_bm5_data_40kms'//trim(filename_extension)
#elif (IMAX==105)
    filename = 'fake_bm5_data_16kms'//trim(filename_extension)
#else
    errormsg = ' >>> '//trim(thisroutine)//': Error when looking for a' &
    //               end_of_line &
    //'              Fake BedMachine data file!'
    call error(errormsg)
#endif
    filename_with_path = trim(temp_path)//'/'//trim(filename)

    !  ------ Open NetCDF file
    ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

    if (ios /= nf90_noerr) then
        errormsg = ' >>> '//trim(thisroutine)//': Error when opening a' &
        //               end_of_line &
        //'              Fake NetCDF BedMachine data file!'
        call error(errormsg)
    end if

    call check( nf90_inq_varid(ncid, 'zs', ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, zs_BedMachine_data_conv), thisroutine )
#ifdef ALLOW_ZS_UNCERT
    call check( nf90_inq_varid(ncid, ZS_UNCERT_FIELD, ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, zs_unc_BedMachine_data_conv), thisroutine )
#endif
    !  ------ Close NetCDF file
    call check( nf90_close(ncid) )

    do i = 0, IMAX
        do j = 0, JMAX
            zs_BedMachine_data(j,i) = zs_BedMachine_data_conv(i,j)
#ifdef ALLOW_ZS_UNCERT
            zs_unc_BedMachine_data(j,i) = zs_unc_BedMachine_data_conv(i,j)
#endif
        end do
    end do
#endif

#if (defined(ZL_COST))
#if (IMAX==168)
    filename = 'bm5_data_10kms'//trim(filename_extension)
#elif (IMAX==42)
    filename = 'bm5_data_40kms'//trim(filename_extension)
#elif (IMAX==105)
    filename = 'bm5_data_16kms'//trim(filename_extension)
#else
    errormsg = ' >>> '//trim(thisroutine)//': Error when looking for a' &
    //               end_of_line &
    //'              BedMachine data file!'
    call error(errormsg)
#endif
    filename_with_path = trim(temp_path)//'/'//trim(filename)

    !  ------ Open NetCDF file
    ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

    if (ios /= nf90_noerr) then
        errormsg = ' >>> '//trim(thisroutine)//': Error when opening a' &
        //               end_of_line &
        //'              NetCDF BedMachine data file!'
        call error(errormsg)
    end if

    call check( nf90_inq_varid(ncid, 'zl', ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, zl_BedMachine_data_conv), thisroutine )
#ifdef ALLOW_ZL_UNCERT
    call check( nf90_inq_varid(ncid, ZL_UNCERT_FIELD, ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, zl_unc_BedMachine_data_conv), thisroutine )
#endif
    !  ------ Close NetCDF file
    call check( nf90_close(ncid) )

    do i = 0, IMAX
        do j = 0, JMAX
            zl_BedMachine_data(j,i) = zl_BedMachine_data_conv(i,j)
#ifdef ALLOW_ZL_UNCERT
            zl_unc_BedMachine_data(j,i) = zl_unc_BedMachine_data_conv(i,j)
#endif
        end do
    end do
#endif

#if (defined(FAKE_ZL_COST))
#if (IMAX==168)
    filename = 'fake_bm5_data_10kms'//trim(filename_extension)
#elif (IMAX==42)
    filename = 'fake_bm5_data_40kms'//trim(filename_extension)
#elif (IMAX==105)
    filename = 'fake_bm5_data_16kms'//trim(filename_extension)
#else
    errormsg = ' >>> '//trim(thisroutine)//': Error when looking for a' &
    //               end_of_line &
    //'              Fake BedMachine data file!'
    call error(errormsg)
#endif
    filename_with_path = trim(temp_path)//'/'//trim(filename)

    !  ------ Open NetCDF file
    ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

    if (ios /= nf90_noerr) then
        errormsg = ' >>> '//trim(thisroutine)//': Error when opening a' &
        //               end_of_line &
        //'              Fake NetCDF BedMachine data file!'
        call error(errormsg)
    end if

    call check( nf90_inq_varid(ncid, 'zl', ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, zl_BedMachine_data_conv), thisroutine )
#ifdef ALLOW_ZL_UNCERT
    call check( nf90_inq_varid(ncid, ZL_UNCERT_FIELD, ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, zl_unc_BedMachine_data_conv), thisroutine )
#endif
    !  ------ Close NetCDF file
    call check( nf90_close(ncid) )

    do i = 0, IMAX
        do j = 0, JMAX
            zl_BedMachine_data(j,i) = zl_BedMachine_data_conv(i,j)
#ifdef ALLOW_ZL_UNCERT
            zl_unc_BedMachine_data(j,i) = zl_unc_BedMachine_data_conv(i,j)
#endif
        end do
    end do
#endif

#if defined(SURFVEL_COST)
#if (IMAX==168)
    filename = 'surfvel_data_10kms'//trim(filename_extension)
#elif (IMAX==42)
    filename = 'surfvel_data_40kms'//trim(filename_extension)
#elif (IMAX==105)
    filename = 'surfvel_data_16kms'//trim(filename_extension)
#else
    errormsg = ' >>> '//trim(thisroutine)//': Error when looking for a' &
    //               end_of_line &
    //'              velocity surface MEaSUREs data file!'
    call error(errormsg)
#endif
    filename_with_path = trim(temp_path)//'/'//trim(filename)

    !  ------ Open NetCDF file
    ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

    if (ios /= nf90_noerr) then
        errormsg = ' >>> '//trim(thisroutine)//': Error when opening a' &
        //               end_of_line &
        //'              NetCDF velocity surface MEaSUREs data file!'
        call error(errormsg)
    end if

#if !defined(SURF_VXVY_COST)

    call check( nf90_inq_varid(ncid, 'vs', ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, vs_MEaSUREs_data_conv), thisroutine )
#ifdef ALLOW_SURFVEL_UNCERT
    call check( nf90_inq_varid(ncid, SURFVEL_UNCERT_FIELD, ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, vs_unc_MEaSUREs_data_conv), thisroutine )
#endif

    !  ------ Close NetCDF file
    call check( nf90_close(ncid) )

    do i = 0, IMAX
        do j = 0, JMAX
            vs_MEaSUREs_data(j,i) = vs_MEaSUREs_data_conv(i,j)
#ifdef ALLOW_SURFVEL_UNCERT
            vs_unc_MEaSUREs_data(j,i) = vs_unc_MEaSUREs_data_conv(i,j)
#endif
        end do
    end do

#else

    call check( nf90_inq_varid(ncid, 'vx', ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, vx_MEaSUREs_data_conv), thisroutine )
#ifdef ALLOW_SURFVEL_UNCERT
    call check( nf90_inq_varid(ncid, SURFVX_UNCERT_FIELD, ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, vx_unc_MEaSUREs_data_conv), thisroutine )
#endif
    call check( nf90_inq_varid(ncid, 'vy', ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, vy_MEaSUREs_data_conv), thisroutine )
#ifdef ALLOW_SURFVEL_UNCERT
    call check( nf90_inq_varid(ncid, SURFVY_UNCERT_FIELD, ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, vy_unc_MEaSUREs_data_conv), thisroutine )
#endif

    !  ------ Close NetCDF file
    call check( nf90_close(ncid) )

    do i = 0, IMAX
        do j = 0, JMAX
            vx_MEaSUREs_data(j,i) = vx_MEaSUREs_data_conv(i,j)
            vy_MEaSUREs_data(j,i) = vy_MEaSUREs_data_conv(i,j)
#ifdef ALLOW_SURFVEL_UNCERT
            vx_unc_MEaSUREs_data(j,i) = vx_unc_MEaSUREs_data_conv(i,j)
            vy_unc_MEaSUREs_data(j,i) = vy_unc_MEaSUREs_data_conv(i,j)
#endif
        end do
    end do

#endif

#endif

#if defined(FAKE_SURFVEL_COST)
#if (IMAX==168)
    filename = 'fake_surfvel_data_10kms'//trim(filename_extension)
#elif (IMAX==42)
    filename = 'fake_surfvel_data_40kms'//trim(filename_extension)
#elif (IMAX==105)
    filename = 'fake_surfvel_data_16kms'//trim(filename_extension)
#else
    errormsg = ' >>> '//trim(thisroutine)//': Error when looking for a' &
    //               end_of_line &
    //'              Fake velocity surface MEaSUREs data file!'
    call error(errormsg)
#endif
    filename_with_path = trim(temp_path)//'/'//trim(filename)

    !  ------ Open NetCDF file
    ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

    if (ios /= nf90_noerr) then
        errormsg = ' >>> '//trim(thisroutine)//': Error when opening a' &
        //               end_of_line &
        //'              Fake NetCDF velocity surface MEaSUREs data file!'
        call error(errormsg)
    end if

#if !defined(SURF_VXVY_COST)

    call check( nf90_inq_varid(ncid, 'vs', ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, vs_MEaSUREs_data_conv), thisroutine )
#ifdef ALLOW_SURFVEL_UNCERT
    call check( nf90_inq_varid(ncid, SURFVEL_UNCERT_FIELD, ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, vs_unc_MEaSUREs_data_conv), thisroutine )
#endif

    !  ------ Close NetCDF file
    call check( nf90_close(ncid) )

    do i = 0, IMAX
        do j = 0, JMAX
            vs_MEaSUREs_data(j,i) = vs_MEaSUREs_data_conv(i,j)
#ifdef ALLOW_SURFVEL_UNCERT
            vs_unc_MEaSUREs_data(j,i) = vs_unc_MEaSUREs_data_conv(i,j)
#endif
       end do
    end do

#else

    call check( nf90_inq_varid(ncid, 'vx', ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, vx_MEaSUREs_data_conv), thisroutine )
#ifdef ALLOW_SURFVEL_UNCERT
    call check( nf90_inq_varid(ncid, SURFVX_UNCERT_FIELD, ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, vx_unc_MEaSUREs_data_conv), thisroutine )
#endif
    call check( nf90_inq_varid(ncid, 'vy', ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, vy_MEaSUREs_data_conv), thisroutine )
#ifdef ALLOW_SURFVEL_UNCERT
    call check( nf90_inq_varid(ncid, SURFVY_UNCERT_FIELD, ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, vy_unc_MEaSUREs_data_conv), thisroutine )
#endif

    !  ------ Close NetCDF file
    call check( nf90_close(ncid) )

    do i = 0, IMAX
        do j = 0, JMAX
            vx_MEaSUREs_data(j,i) = vx_MEaSUREs_data_conv(i,j)
            vy_MEaSUREs_data(j,i) = vy_MEaSUREs_data_conv(i,j)
#ifdef ALLOW_SURFVEL_UNCERT
            vx_unc_MEaSUREs_data(j,i) = vx_unc_MEaSUREs_data_conv(i,j)
            vy_unc_MEaSUREs_data(j,i) = vy_unc_MEaSUREs_data_conv(i,j)
#endif
        end do
    end do

#endif

#endif

#ifdef AGE_COST

#if (CALCMOD!=1)
  KDATA = KCMAX
#else 
  KDATA = KCMAX + KTMAX
  errormsg = ' >>> '//trim(thisroutine)//': Age model-data misfit not compatible' &
  //               end_of_line &
  //'              with CALCMOD==1!'
  call error(errormsg)
#endif

#if (IMAX==168)
    filename = 'age_data_10kms'//trim(filename_extension)
#elif (IMAX==42)
    filename = 'age_data_40kms'//trim(filename_extension)
#elif (IMAX==105)
    filename = 'age_data_16kms'//trim(filename_extension)
#else
    errormsg = ' >>> '//trim(thisroutine)//': Error when looking for a' &
    //               end_of_line &
    //'              Age data file!'
    call error(errormsg)
#endif
    filename_with_path = trim(temp_path)//'/'//trim(filename)

    !  ------ Open NetCDF file
    ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

    if (ios /= nf90_noerr) then
        errormsg = ' >>> '//trim(thisroutine)//': Error when opening a' &
        //               end_of_line &
        //'              NetCDF Age data file!'
        call error(errormsg)
    end if

    call check( nf90_inq_varid(ncid, 'age_c', ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, age_data_conv), thisroutine )
#ifdef ALLOW_AGE_UNCERT
    call check( nf90_inq_varid(ncid, AGE_UNCERT_FIELD, ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, age_unc_data_conv), thisroutine )
#endif
    !  ------ Close NetCDF file
    call check( nf90_close(ncid) )

    do kc=0, KDATA
        do j=0, JMAX
            do i=0, IMAX
                age_data(kc,j,i) = age_data_conv(i,j,kc)
#ifdef ALLOW_AGE_UNCERT
                age_unc_data(kc,j,i) = age_unc_data_conv(i,j,kc)
#endif
            end do
        end do
    end do
#endif

#ifdef FAKE_AGE_COST

#if (CALCMOD!=1)
  KDATA = KCMAX
#else
  KDATA = KCMAX + KTMAX
  errormsg = ' >>> '//trim(thisroutine)//': Fake Age model-data misfit not compatible' &
  //               end_of_line &
  //'              with CALCMOD==1!'
  call error(errormsg)
#endif

#if (IMAX==168)
    filename = 'fake_age_data_10kms'//trim(filename_extension)
#elif (IMAX==42)
    filename = 'fake_age_data_40kms'//trim(filename_extension)
#elif (IMAX==105)
    filename = 'fake_age_data_16kms'//trim(filename_extension)
#else
    errormsg = ' >>> '//trim(thisroutine)//': Error when looking for a' &
    //               end_of_line &
    //'              Fake Age data file!'
    call error(errormsg)
#endif
    filename_with_path = trim(temp_path)//'/'//trim(filename)

    !  ------ Open NetCDF file
    ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

    if (ios /= nf90_noerr) then
        errormsg = ' >>> '//trim(thisroutine)//': Error when opening a' &
        //               end_of_line &
        //'              Fake NetCDF Age data file!'
        call error(errormsg)
    end if

    call check( nf90_inq_varid(ncid, 'age_c', ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, age_data_conv), thisroutine )
#ifdef ALLOW_AGE_UNCERT
    call check( nf90_inq_varid(ncid, AGE_UNCERT_FIELD, ncv), thisroutine )
    call check( nf90_get_var(ncid, ncv, age_unc_data_conv), thisroutine )
#endif
    !  ------ Close NetCDF file
    call check( nf90_close(ncid) )

    do kc=0, KDATA
        do j=0, JMAX
            do i=0, IMAX
                age_data(kc,j,i) = age_data_conv(i,j,kc)
#ifdef ALLOW_AGE_UNCERT
                age_unc_data(kc,j,i) = age_unc_data_conv(i,j,kc)
#endif
            end do
        end do
    end do
#endif

end subroutine read_cost_data
  
end module cost_io_m
 
