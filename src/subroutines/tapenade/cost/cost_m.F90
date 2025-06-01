!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module  :  c o s t _ m
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
module cost_m

  use sico_types_m  
  use sico_variables_m
  use cost_io_m
  use error_m
#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  implicit none

  public :: cost_final, laplace_smoothing_2D_reg_cost, laplace_smoothing_3D_reg_cost

contains
 
!-------------------------------------------------------------------------------
!> This is the final cost calculation.
!! The cost function structure is defined here. 
!!
!! Currently there are three options - 
!!
!! 1. AGE_COST
!! Currently is a "observed age" - modeled age summed over the entire 
!! domain. The "observed age" is a fake, generated age
!! field performed by the 125 ka run in headers.
!!
!! 2. BEDMACHINE_COST
!! Based on L2 misfit with dataset in -
!! BedMachine v3: Complete Bed Topography and Ocean Bathymetry Mapping of
!! Greenland From Multibeam Echo Sounding Combined With Mass Conservation
!! by Morlighem et. al in 2017.
!!
!! 3. Total volume - Default (Only ALLOW_COST used)
!! This just defines the cost function to be the total volume of the ice sheet.
!! This is currently the default.
!!
!! Other cost functions are certainly possible, and recommended! 
!-------------------------------------------------------------------------------
  subroutine cost_final()
  
  implicit none
  
  integer(i4b) :: i, j, k, kc, kt, tad, ios, KDATA, ctrl_index
  character(len=64), parameter :: thisroutine = 'cost_final'
  real(dp), dimension(0:JMAX,0:IMAX) :: vs

  !-------- Initialize the cost functions:
  fc = 0.0
  fc_data = 0.0
  fc_reg = 0.0
  fc_bm5 = 0.0
  fc_ac = 0.0
  fc_svc = 0.0
  fc_vxc = 0.0
  fc_vyc = 0.0
  fc_zsc = 0.0
  fc_zlc = 0.0

  !-------- Read any necessary NetCDF cost files:
  call read_cost_data()

#if (defined(AGE_COST) || defined(FAKE_AGE_COST))

#if (CALCMOD!=1)
  KDATA = KCMAX
#else 
  KDATA = KCMAX + KTMAX
  errormsg = ' >>> '//trim(thisroutine)//': Age model-data misfit not compatible' &
  //               end_of_line &
  //'              with CALCMOD==1!'
  call error(errormsg)
#endif
 
  do i=0, IMAX 
    do j=0, JMAX
      do k=0, KDATA
#ifdef ALLOW_AGE_UNCERT
        ! only counting points that are real in the data: 
        if (age_unc_data(k,j,i) .gt. 0.0 .and. age_data(k,j,i) .ge. 0.0 .and. age_data(k,j,i) .le. 134000.0 .and. H_BedMachine_data(j,i) .ge. 1500.0) then
          fc = fc &
          + 0.5*(age_data(k,j,i)*year2sec - age_c(k,j,i))**2/(age_unc_data(k,j,i)*year2sec)**2
#else
        ! only counting points that are real in the data:
        if (age_data(k,j,i) .ge. 0.0 .and. age_data(k,j,i) .le. 134000.0 .and. H_BedMachine_data(j,i) .ge. 1500.0) then
          fc = fc &
          + 0.5*(age_data(k,j,i)*year2sec - age_c(k,j,i))**2
#endif
        end if
      end do
    end do
  end do

  fc_ac = fc
  print *, 'Final age cost, fc_ac = ', fc_ac
#endif

#if (defined(BEDMACHINE_COST) || defined(FAKE_BEDMACHINE_COST))
    do i=0, IMAX
      do j=0, JMAX
        fc = fc &
#ifdef ALLOW_BEDMACHINE_UNCERT
        + 0.5*(H(j,i) - H_BedMachine_data(j,i))**2/H_unc_BedMachine_data(j,i)**2
#else
        + 0.5*(H(j,i) - H_BedMachine_data(j,i))**2
#endif
      end do
    end do

  fc_bm5 = fc - fc_ac
  print *, 'Final BM5 cost, fc_bm5 = ', fc_bm5
#endif

#if (defined(ZS_COST) || defined(FAKE_ZS_COST))
    do i=0, IMAX
      do j=0, JMAX
        fc = fc &
#ifdef ALLOW_ZS_UNCERT
        + 0.5*(zs(j,i) - zs_BedMachine_data(j,i))**2/zs_unc_BedMachine_data(j,i)**2
#else
        + 0.5*(zs(j,i) - zs_BedMachine_data(j,i))**2
#endif
      end do
    end do

  fc_zsc = fc - (fc_ac + fc_bm5)
  print *, 'Final zs cost, fc_zsc = ', fc_zsc
#endif

#if (defined(ZL_COST) || defined(FAKE_ZL_COST))
    do i=0, IMAX
      do j=0, JMAX
        fc = fc &
#ifdef ALLOW_ZL_UNCERT
        + 0.5*(zl(j,i) - zl_BedMachine_data(j,i))**2/zl_unc_BedMachine_data(j,i)**2
#else
        + 0.5*(zl(j,i) - zl_BedMachine_data(j,i))**2
#endif
      end do
    end do

  fc_zlc = fc - (fc_ac + fc_bm5 + fc_zsc)
  print *, 'Final zl cost, fc_zlc = ', fc_zlc
#endif

#if (defined(SURFVEL_COST) || defined(FAKE_SURFVEL_COST))
    do i=0, IMAX
      do j=0, JMAX
        if (zs_BedMachine_data(j,i) .ge. -50.0) then

#if !defined(SURF_VXVY_COST)

#if !defined(ALLOW_TAPENADE)
          vs(j,i) = sqrt(vx_s_g(j,i)**2 + vy_s_g(j,i)**2)*year2sec
#else /* ALLOW_TAPENADE: guarding against non-differentiable sqrt(0) */
          if ((vx_s_g(j,i)**2 + vy_s_g(j,i)**2) > 0) then
            vs(j,i) = sqrt(vx_s_g(j,i)**2 + vy_s_g(j,i)**2)*year2sec
          else
            vs(j,i) = 0.0
          end if
#endif
          fc = fc &
#ifdef ALLOW_SURFVEL_UNCERT
          + 0.5*(vs(j,i) - vs_MEaSUREs_data(j,i))**2/vs_unc_MEaSUREs_data(j,i)**2
#else
          + 0.5*(vs(j,i) - vs_MEaSUREs_data(j,i))**2
#endif

#else

          fc = fc &
#ifdef ALLOW_SURFVEL_UNCERT
          + 0.5*(vx_s_g(j,i) - vx_MEaSUREs_data(j,i)*sec2year)**2/(vx_unc_MEaSUREs_data(j,i)*sec2year)**2 &
          + 0.5*(vy_s_g(j,i) - vy_MEaSUREs_data(j,i)*sec2year)**2/(vy_unc_MEaSUREs_data(j,i)*sec2year)**2
#else
          + 0.5*(vx_s_g(j,i) - vx_MEaSUREs_data(j,i)*sec2year)**2 &
          + 0.5*(vy_s_g(j,i) - vy_MEaSUREs_data(j,i)*sec2year)**2
#endif
          fc_vxc = fc_vxc &
#ifdef ALLOW_SURFVEL_UNCERT
          + 0.5*(vx_s_g(j,i) - vx_MEaSUREs_data(j,i)*sec2year)**2/(vx_unc_MEaSUREs_data(j,i)*sec2year)**2
#else
          + 0.5*(vx_s_g(j,i) - vx_MEaSUREs_data(j,i)*sec2year)**2
#endif
          fc_vyc = fc_vyc &
#ifdef ALLOW_SURFVEL_UNCERT
          + 0.5*(vy_s_g(j,i) - vy_MEaSUREs_data(j,i)*sec2year)**2/(vy_unc_MEaSUREs_data(j,i)*sec2year)**2
#else
          + 0.5*(vy_s_g(j,i) - vy_MEaSUREs_data(j,i)*sec2year)**2
#endif

#endif
        end if
      end do
    end do

#if !defined(SURF_VXVY_COST)
  fc_svc = fc - (fc_ac + fc_bm5 + fc_zsc + fc_zlc)
  print *, 'Final vs cost, fc_svc = ', fc_svc
#else
  print *, 'Final vx cost, fc_vxc = ', fc_vxc
  print *, 'Final vy cost, fc_vyc = ', fc_vyc
#endif
#endif

#if (!defined(BEDMACHINE_COST) && !defined(AGE_COST) && !defined(SURFVEL_COST) && !defined(ZS_COST) && !defined(ZL_COST))
#if (!defined(FAKE_BEDMACHINE_COST) && !defined(FAKE_AGE_COST) && !defined(FAKE_SURFVEL_COST) && !defined(FAKE_ZS_COST) && !defined(FAKE_ZL_COST))
    do i=0, IMAX
      do j=0, JMAX
        !--- Total volume cost function:
        fc = fc + H(j,i)*cell_area(j,i)
      end do
    end do
#endif
#endif

  fc_data = fc

#if (defined(ALLOW_GENCTRL) && defined(DO_GENCTRL_PRIOR))

#if (defined(DO_CTRL_GENARR2D) && defined(XX_GENARR2D_VARS_ARR))
  do ctrl_index = 1, NUM_CTRL_GENARR2D
    call laplace_smoothing_2D_reg_cost(xx_genarr2d_orig(ctrl_index,:,:), &
                                       xx_genarr2d_prior(ctrl_index,:,:), &
                                       xx_genarr2d_prior_X(ctrl_index,:,:), &
                                       genarr2d_gamma_arr(ctrl_index), &
                                       genarr2d_delta_arr(ctrl_index), &
                                       genarr2d_sigma_arr(ctrl_index))
  end do
#endif

#if (defined(DO_CTRL_GENARR3D) && defined(XX_GENARR3D_VARS_ARR))
  do ctrl_index = 1, NUM_CTRL_GENARR3D
    call laplace_smoothing_3D_reg_cost(xx_genarr3d_orig(ctrl_index,:,:,:), &
                                       xx_genarr3d_prior(ctrl_index,:,:,:), &
                                       xx_genarr3d_prior_X(ctrl_index,:,:,:), &
                                       genarr3d_gamma_arr(ctrl_index), &
                                       genarr3d_delta_arr(ctrl_index), &
                                       genarr2d_sigma_arr(ctrl_index))
  end do
#endif

#if (defined(DO_CTRL_GENTIM2D) && defined(XX_GENTIM2D_VARS_ARR))
  do ctrl_index = 1, NUM_CTRL_GENTIM2D
    do tad = 0, NTDAMAX
      call laplace_smoothing_2D_reg_cost(xx_gentim2d_orig(ctrl_index,tad,:,:), &
                                         xx_gentim2d_prior(ctrl_index,tad,:,:), &
                                         xx_gentim2d_prior_X(ctrl_index,tad,:,:), &
                                         gentim2d_gamma_arr(ctrl_index), &
                                         gentim2d_delta_arr(ctrl_index), &
                                         gentim2d_sigma_arr(ctrl_index))
    end do
  end do
#endif

#endif

  fc = fc + fc_reg
  
  !-------- Print to screen just in case something gets
  !         crazy with the file outputting:
  print *, 'Final model-data misfit cost, fc_data = ', fc_data
  print *, 'Final prior/regularization cost, fc_reg = ', fc_reg
  print *, 'Final cost, fc = ', fc
  print *, trim(OUT_PATH)

  end subroutine cost_final

  subroutine laplace_smoothing_2D_reg_cost(field, field_prior, field_prior_X, gamm, delta, sigma)

  implicit none

  real(dp), dimension(0:JMAX,0:IMAX) :: field, field_prior, field_prior_X

  integer(i4b) :: i, j
  character(len=64), parameter :: thisroutine = 'laplace_smoothing_2D_reg_cost'
  real(dp) :: delta, gamm, sigma

  field(:,:) = field(:,:) / (field_prior_X(:,:)*sigma)
  field_prior(:,:) = field_prior(:,:) / (field_prior_X(:,:)*sigma)

  fc_reg = fc_reg + 0.5*(delta*(field(0,0)-field_prior(0,0)) &
                               - gamm*((field(0,1) + field(1,0) - 2*field(0,0)) &
                               -(field_prior(0,1) + field_prior(1,0) - 2*field_prior(0,0))) / DX**2)**2
  fc_reg = fc_reg + 0.5*(delta*(field(JMAX,0)-field_prior(JMAX,0)) &
                               - gamm*((field(JMAX,1) + field(JMAX-1,0) - 2*field(JMAX,0)) &
                               -(field_prior(JMAX,1) + field_prior(JMAX-1,0) - 2*field_prior(JMAX,0))) / DX**2)**2
  fc_reg = fc_reg + 0.5*(delta*(field(0,IMAX)-field_prior(0,IMAX)) &
                               - gamm*((field(0,IMAX-1) + field(1,IMAX) - 2*field(0,IMAX)) &
                               -(field_prior(0,IMAX-1) + field_prior(1,IMAX) - 2*field_prior(0,IMAX))) / DX**2)**2
  fc_reg = fc_reg + 0.5*(delta*(field(JMAX,IMAX)-field_prior(JMAX,IMAX)) &
                               - gamm*((field(JMAX,IMAX-1) + field(JMAX-1,IMAX) - 2*field(JMAX,IMAX)) &
                               -(field_prior(JMAX,IMAX-1) + field_prior(JMAX-1,IMAX) - 2*field_prior(JMAX,IMAX))) / DX**2)**2

  do i=1, IMAX-1
    fc_reg = fc_reg &
           + 0.5*(delta*(field(0,i)-field_prior(0,i)) &
                        - gamm*(((field(1,i) - field(0,i)) - (field_prior(1,i) - field_prior(0,i))) / DX**2 &
                        + ((field(0,i-1) - 2*field(0,i) + field(0,i+1)) - (field_prior(0,i-1) - 2*field_prior(0,i) + field_prior(0,i+1))) / DX**2))**2
    fc_reg = fc_reg &
           + 0.5*(delta*(field(JMAX,i)-field_prior(JMAX,i)) &
                        - gamm*(((field(JMAX-1,i) - field(JMAX,i)) - (field_prior(JMAX-1,i) - field_prior(JMAX,i))) / DX**2 &
                        + ((field(JMAX,i-1) - 2*field(JMAX,i) + field(JMAX,i+1)) - (field_prior(JMAX,i-1) - 2*field_prior(JMAX,i) + field_prior(JMAX,i+1))) / DX**2))**2
  end do

  do j=1, JMAX-1
    fc_reg = fc_reg &
           + 0.5*(delta*(field(j,0)-field_prior(j,0)) &
                        - gamm*(((field(j-1,0) - 2*field(j,0) + field(j+1,0)) - (field_prior(j-1,0) - 2*field_prior(j,0) + field_prior(j+1,0))) / DX**2 &
                                 + ((field(j,1) - field(j,0)) - (field_prior(j,1) - field_prior(j,0))) / DX**2))**2
    fc_reg = fc_reg &
           + 0.5*(delta*(field(j,IMAX)-field_prior(j,IMAX)) &
                        - gamm*(((field(j-1,IMAX) - 2*field(j,IMAX) + field(j+1,IMAX)) - (field_prior(j-1,IMAX) - 2*field_prior(j,IMAX) + field_prior(j+1,IMAX))) / DX**2 &
                                + ((field(j,IMAX-1) - field(j,IMAX)) - (field_prior(j,IMAX-1) - field_prior(j,IMAX))) / DX**2))**2
  end do

  do i=1, IMAX-1
    do j=1, JMAX-1
      fc_reg = fc_reg &
      + 0.5*(delta*(field(j,i)-field_prior(j,i)) &
          - gamm*((field(j,i-1) - 2*field(j,i) + field(j,i+1)) - (field_prior(j,i-1) - 2*field_prior(j,i) + field_prior(j,i+1))) / DX**2 &
          - gamm*((field(j-1,i) - 2*field(j,i) + field(j+1,i)) - (field_prior(j-1,i) - 2*field_prior(j,i) + field_prior(j+1,i))) / DX**2)**2
    end do
  end do

  end subroutine laplace_smoothing_2D_reg_cost

  subroutine laplace_smoothing_3D_reg_cost(field, field_prior, field_prior_X, gamm, delta, sigma)

  implicit none

  real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: field, field_prior, field_prior_X

  integer(i4b) :: i, j, kc
  character(len=64), parameter :: thisroutine = 'laplace_smoothing_3D_reg_cost'
  real(dp) :: delta, gamm, sigma
  real(dp), dimension(1:KCMAX) :: delta_z

  do kc=1, KCMAX
    delta_z(kc) = 1.e6*(eaz_c(kc)-eaz_c(kc-1))
  end do

  field(:,:,:) = field(:,:,:) / (field_prior_X(:,:,:)*sigma)
  field_prior(:,:,:) = field_prior(:,:,:) / (field_prior_X(:,:,:)*sigma)

  do i=0,IMAX
    do j=0,JMAX
      fc_reg = fc_reg + 0.5*(delta*(field(0,j,i)-field_prior(0,j,i)) &
                                         - gamm*((field(1,j,i) - field(0,j,i)) &
                                         -(field_prior(1,j,i) - field_prior(0,j,i))) / delta_z(1)**2)**2
      fc_reg = fc_reg + 0.5*(delta*(field(KCMAX,j,i)-field_prior(KCMAX,j,i)) &
                                         - gamm*((field(KCMAX-1,j,i) - field(KCMAX,j,i)) &
                                         -(field_prior(KCMAX-1,j,i) - field_prior(KCMAX,j,i))) / delta_z(KCMAX)**2)**2
    end do
  end do

  do i=0,IMAX
    do j=0,JMAX
      do kc=1, KCMAX-1
        fc_reg = fc_reg &
        + 0.5*(delta*(field(kc,j,i)-field_prior(kc,j,i)) &
            - gamm*(((field(kc+1,j,i)-field(kc,j,i))/delta_z(kc+1) - (field(kc,j,i)-field(kc-1,j,i))/delta_z(kc))*(2.0/(delta_z(kc) + delta_z(kc+1))) &
            -((field_prior(kc+1,j,i)-field_prior(kc,j,i))/delta_z(kc+1) - (field_prior(kc,j,i)-field_prior(kc-1,j,i))/delta_z(kc))*(2.0/(delta_z(kc) + delta_z(kc+1)))))**2
      end do
    end do
  end do

  end subroutine laplace_smoothing_3D_reg_cost

end module cost_m
