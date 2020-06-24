!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ v x y _ m
!
!> @file
!!
!! Computation of the horizontal velocity vx, vy.
!!
!! @section Copyright
!!
!! Copyright 2009-2020 Ralf Greve, Tatsuru Sato, Thomas Goelles, Jorge Bernales
!!
!! @section License
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
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with SICOPOLIS.  If not, see <http://www.gnu.org/licenses/>.
!<
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------------------------------------------------------------
!> Computation of the horizontal velocity vx, vy.
!<------------------------------------------------------------------------------
module calc_vxy_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m
  use error_m

  implicit none

  real(dp), dimension(0:JMAX,0:IMAX), save :: dzs_dx_aux, dzs_dy_aux

  private
  public :: calc_dzs_dxy_aux, &
            calc_vxy_b_sia, calc_vxy_sia, calc_vxy_static, calc_vxy_ssa

contains

!-------------------------------------------------------------------------------
!> Computation of the auxiliary surface gradients dzs_dx_aux, dzs_dy_aux
!! (optional one-sided gradients at the grounding line).
!<------------------------------------------------------------------------------
subroutine calc_dzs_dxy_aux(z_sl, dxi, deta)

implicit none

real(dp), intent(in) :: z_sl, dxi, deta

integer(i4b) :: i, j
real(dp)     :: inv_dx, inv_dy
real(dp)     :: rhosw_rho_ratio
real(dp)     :: H_mid, zl_mid, zs_mid

inv_dx          = 1.0_dp/dxi
inv_dy          = 1.0_dp/deta
rhosw_rho_ratio = RHO_SW/RHO

dzs_dx_aux = dzs_dxi
dzs_dy_aux = dzs_deta

#if (MARGIN==3)

#if (!defined(GL_SURF_GRAD) || GL_SURF_GRAD==1)

!!! continue

#elif (GL_SURF_GRAD==2)

do i=0, IMAX-1
do j=1, JMAX-1
   ! inner point on the staggered grid in x-direction

   if ( (flag_grounding_line_1(j,i).and.flag_grounding_line_2(j,i+1)) &
        .or. &
        (flag_grounding_line_2(j,i).and.flag_grounding_line_1(j,i+1)) &
      ) then
        ! one neighbour is floating ice and the other is grounded ice
        ! (grounding line)

      H_mid  = 0.5_dp*((H_c(j,i)+H_t(j,i))+(H_c(j,i+1)+H_t(j,i+1)))
      zl_mid = 0.5_dp*(zl(j,i)+zl(j,i+1))
      zs_mid = 0.5_dp*(zs(j,i)+zs(j,i+1))

      if (H_mid < (z_sl-zl_mid)*rhosw_rho_ratio) then
         ! floating condition is satisfied

         if ( &
              (flag_grounding_line_1(j,i).and.flag_grounding_line_2(j,i+1)) &
              .and. &
              (i+2 <= IMAX) &
            ) then

            if ((maske(j,i+2) == 3_i1b).or.(maske(j,i+2) == 2_i1b)) &
               dzs_dx_aux(j,i) = (0.5_dp*(zs(j,i+1)+zs(j,i+2))-zs_mid)*inv_dx
                                 ! one-sided gradient into floating ice

         else if ( &
              (flag_grounding_line_2(j,i).and.flag_grounding_line_1(j,i+1)) &
              .and. &
              (i-1 >= 0) &
            ) then

            if ((maske(j,i-1) == 3_i1b).or.(maske(j,i-1) == 2_i1b)) &
               dzs_dx_aux(j,i) = (zs_mid-0.5_dp*(zs(j,i)+zs(j,i-1)))*inv_dx
                                 ! one-sided gradient into floating ice

         end if

      else   ! H_mid >= (z_sl-zl_mid)*rhosw_rho_ratio,
             ! floating condition is not satisfied

         if ( &
              (flag_grounding_line_1(j,i).and.flag_grounding_line_2(j,i+1)) &
              .and. &
              (i-1 >= 0) &
            ) then

            if ((maske(j,i-1) == 0_i1b).or.(maske(j,i-1) == 1_i1b)) &
               dzs_dx_aux(j,i) = (zs_mid-0.5_dp*(zs(j,i)+zs(j,i-1)))*inv_dx
                                 ! one-sided gradient into grounded ice

         else if ( &
              (flag_grounding_line_2(j,i).and.flag_grounding_line_1(j,i+1)) &
              .and. &
              (i+2 <= IMAX) &
            ) then

            if ((maske(j,i+2) == 0_i1b).or.(maske(j,i+2) == 1_i1b)) &
               dzs_dx_aux(j,i) = (0.5_dp*(zs(j,i+1)+zs(j,i+2))-zs_mid)*inv_dx
                                 ! one-sided gradient into grounded ice

         end if

      end if


   end if

end do
end do

do i=1, IMAX-1
do j=0, JMAX-1
   ! inner point on the staggered grid in y-direction

   if ( (flag_grounding_line_1(j,i).and.flag_grounding_line_2(j+1,i)) &
        .or. &
        (flag_grounding_line_2(j,i).and.flag_grounding_line_1(j+1,i)) &
      ) then
        ! one neighbour is floating ice and the other is grounded ice
        ! (grounding line)

      H_mid  = 0.5_dp*((H_c(j,i)+H_t(j,i))+(H_c(j+1,i)+H_t(j+1,i)))
      zl_mid = 0.5_dp*(zl(j,i)+zl(j+1,i))
      zs_mid = 0.5_dp*(zs(j,i)+zs(j+1,i))

      if (H_mid < (z_sl-zl_mid)*rhosw_rho_ratio) then
         ! floating condition is satisfied

         if ( &
              (flag_grounding_line_1(j,i).and.flag_grounding_line_2(j+1,i)) &
              .and. &
              (j+2 <= JMAX) &
            ) then

            if ((maske(j+2,i) == 3_i1b).or.(maske(j+2,i) == 2_i1b)) &
               dzs_dy_aux(j,i) = (0.5_dp*(zs(j+1,i)+zs(j+2,i))-zs_mid)*inv_dy
                                 ! one-sided gradient into floating ice

         else if ( &
              (flag_grounding_line_2(j,i).and.flag_grounding_line_1(j+1,i)) &
              .and. &
              (j-1 >= 0) &
            ) then

            if ((maske(j-1,i) == 3_i1b).or.(maske(j-1,i) == 2_i1b)) &
               dzs_dy_aux(j,i) = (zs_mid-0.5_dp*(zs(j,i)+zs(j-1,i)))*inv_dy
                                 ! one-sided gradient into floating ice

         end if

      else   ! H_mid >= (z_sl-zl_mid)*rhosw_rho_ratio,
             ! floating condition is not satisfied

         if ( &
              (flag_grounding_line_1(j,i).and.flag_grounding_line_2(j+1,i)) &
              .and. &
              (j-1 >= 0) &
            ) then

            if ((maske(j-1,i) == 0_i1b).or.(maske(j-1,i) == 1_i1b)) &
               dzs_dy_aux(j,i) = (zs_mid-0.5_dp*(zs(j,i)+zs(j-1,i)))*inv_dy
                                 ! one-sided gradient into grounded ice

         else if ( &
              (flag_grounding_line_2(j,i).and.flag_grounding_line_1(j+1,i)) &
              .and. &
              (j+2 <= JMAX) &
            ) then

            if ((maske(j+2,i) == 0_i1b).or.(maske(j+2,i) == 1_i1b)) &
               dzs_dy_aux(j,i) = (0.5_dp*(zs(j+1,i)+zs(j+2,i))-zs_mid)*inv_dy
                                 ! one-sided gradient into grounded ice

         end if

      end if

   end if

end do
end do

#else

errormsg = ' >>> calc_dzs_dxy_aux: ' &
              //'GL_SURF_GRAD must be either 1 or 2!'
call error(errormsg)

#endif

#endif

end subroutine calc_dzs_dxy_aux

!-------------------------------------------------------------------------------
!> Computation of the basal horizontal velocity vx_b, vy_b in the shallow ice
!! approximation.
!<------------------------------------------------------------------------------
subroutine calc_vxy_b_sia(time, z_sl)

implicit none

real(dp), intent(in) :: time, z_sl

integer(i4b) :: i, j, n, n_slide_regions
#if (!defined(N_SLIDE_REGIONS) || N_SLIDE_REGIONS<=1)
integer(i4b) :: p_weert_aux(1)
integer(i4b) :: q_weert_aux(1)
real(dp) :: c_slide_aux(1)
real(dp) :: gamma_slide_aux(1)
real(dp) :: gamma_slide_inv_aux(1)
#else
integer(i4b) :: p_weert_aux(N_SLIDE_REGIONS)
integer(i4b) :: q_weert_aux(N_SLIDE_REGIONS)
real(dp) :: c_slide_aux(N_SLIDE_REGIONS)
real(dp) :: gamma_slide_aux(N_SLIDE_REGIONS)
real(dp) :: gamma_slide_inv_aux(N_SLIDE_REGIONS)
#endif
real(dp), dimension(0:JMAX,0:IMAX) :: gamma_slide_inv
real(dp), dimension(0:JMAX,0:IMAX) :: p_b, p_b_red, p_b_red_lim, tau_b
real(dp) :: cvxy1, cvxy1a, cvxy1b, ctau1, ctau1a, ctau1b
real(dp) :: temp_diff
real(dp) :: c_Hw_slide, Hw0_slide, Hw0_slide_inv, ratio_Hw_slide
real(dp) :: vh_max, vh_max_inv
real(dp) :: year_sec_inv, time_in_years
real(dp) :: ramp_up_factor
logical, dimension(0:JMAX,0:IMAX) :: sub_melt_flag

year_sec_inv  = 1.0_dp/YEAR_SEC
time_in_years = time*year_sec_inv

!-------- Sliding-law coefficients --------

#if (!defined(N_SLIDE_REGIONS) || N_SLIDE_REGIONS<=1)
n_slide_regions = 1
#else
n_slide_regions = N_SLIDE_REGIONS
#endif

p_weert_aux = P_WEERT
q_weert_aux = Q_WEERT
c_slide_aux = C_SLIDE
gamma_slide_aux = GAMMA_SLIDE

do n=1, n_slide_regions
   gamma_slide_inv_aux(n) = 1.0_dp/max(gamma_slide_aux(n), eps)
end do

do i=0, IMAX
do j=0, JMAX
   if ( (n_slide_region(j,i) >= 1) &
        .and. &
        (n_slide_region(j,i) <= n_slide_regions) ) then
      p_weert(j,i)         = p_weert_aux(n_slide_region(j,i))
      q_weert(j,i)         = q_weert_aux(n_slide_region(j,i))
      c_slide(j,i)         = c_slide_aux(n_slide_region(j,i))*year_sec_inv
      gamma_slide_inv(j,i) = gamma_slide_inv_aux(n_slide_region(j,i))
      sub_melt_flag(j,i)   = (gamma_slide_aux(n_slide_region(j,i)) >= eps)
   else
      errormsg = ' >>> calc_vxy_b_sia: ' &
                    //'Region number out of allowed range!'
      call error(errormsg)
   end if
end do
end do

do i=0, IMAX
do j=0, JMAX
   p_weert_inv(j,i) = 1.0_dp/max(real(p_weert(j,i),dp), eps)
end do
end do

!  ------ Ramping up basal sliding

ramp_up_factor = 1.0_dp

#if (defined(TIME_RAMP_UP_SLIDE))

if ( (TIME_RAMP_UP_SLIDE > eps_dp) &
     .and. &
     ((time_in_years-(TIME_INIT0)) < (TIME_RAMP_UP_SLIDE)) ) then

   ramp_up_factor = (time_in_years-(TIME_INIT0))/(TIME_RAMP_UP_SLIDE)

   ramp_up_factor = max(min(ramp_up_factor, 1.0_dp), 0.0_dp)
                                  ! constrain to interval [0,1]

   ramp_up_factor = ramp_up_factor*ramp_up_factor*ramp_up_factor &
                                  *(10.0_dp + ramp_up_factor &
                                              *(-15.0_dp+6.0_dp*ramp_up_factor))
                                  ! make transition smooth (quintic function)

   c_slide = c_slide * ramp_up_factor

end if

#endif

!  ------ Coefficients for the contribution of the basal water layer

#if (BASAL_HYDROLOGY==1)

#if (!defined(HYDRO_SLIDE_SAT_FCT))
  errormsg = ' >>> calc_vxy_b_sia: ' &
                //'HYDRO_SLIDE_SAT_FCT must be defined!'
  call error(errormsg)
#endif

#if (defined(C_HW_SLIDE))
  c_Hw_slide = C_HW_SLIDE
#else
  errormsg = ' >>> calc_vxy_b_sia: C_HW_SLIDE must be defined!'
  call error(errormsg)
#endif

#if (defined(HW0_SLIDE))
  Hw0_slide = HW0_SLIDE
#else
  errormsg = ' >>> calc_vxy_b_sia: HW0_SLIDE must be defined!'
  call error(errormsg)
#endif

  Hw0_slide_inv = 1.0_dp/max(Hw0_slide, eps_dp)

#else   /* BASAL_HYDROLOGY==0 */

  c_Hw_slide    = 0.0_dp   ! dummy value
  Hw0_slide     = 1.0_dp   ! dummy value
  Hw0_slide_inv = 1.0_dp   ! dummy value

#endif

!-------- Computation of basal stresses --------

!  ------ Basal pressure p_b, basal water pressure p_b_w,
!         reduced pressure p_b_red

do i=0, IMAX
do j=0, JMAX

   if ((maske(j,i) == 0_i1b).or.flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

      p_b(j,i)         = max(RHO*G*(H_c(j,i)+H_t(j,i)), 0.0_dp)
      p_b_w(j,i)       = RHO_SW*G*max((z_sl-zb(j,i)), 0.0_dp)
      p_b_red(j,i)     = max(p_b(j,i)-p_b_w(j,i), 0.0_dp)
      p_b_red_lim(j,i) = max(p_b_red(j,i), RED_PRES_LIMIT_FACT*p_b(j,i))
                         ! in order to avoid very small values, which may lead
                         ! to huge sliding velocities in the SIA

   else   ! maske(j,i) == 1_i1b, 2_i1b or 3_i1b away from the grounding line

      p_b(j,i)         = 0.0_dp
      p_b_w(j,i)       = 0.0_dp
      p_b_red(j,i)     = 0.0_dp
      p_b_red_lim(j,i) = 0.0_dp

   end if

end do
end do

!  ------ Absolute value of the basal shear stress, tau_b

do i=0, IMAX
do j=0, JMAX

#if !defined(ALLOW_OPENAD) /* Normal */

   tau_b(j,i) = p_b(j,i)*sqrt(dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2)

#else /* OpenAD: guarding against non-differentiable sqrt(0) */

   if ((dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2) > 0) then
      tau_b(j,i) = p_b(j,i)*sqrt(dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2)
   else
      tau_b(j,i) = 0.0_dp
   end if

#endif /* Normal vs. OpenAD */

end do
end do

!-------- Computation of d_help_b (defined on the grid points (i,j)) --------

do i=0, IMAX
do j=0, JMAX

   if ((maske(j,i) == 0_i1b).or.flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

!  ------ Abbreviations

#if (SLIDE_LAW==1)
      cvxy1 = c_slide(j,i) &
              * ( (tau_b(j,i)+eps_dp)**(p_weert(j,i)-1) &
                  /(p_b(j,i)+eps_dp)**q_weert(j,i) ) &
              * p_b(j,i)
      ctau1 = 1.0_dp/(c_slide(j,i)+eps_dp)**p_weert_inv(j,i) &
              * (p_b(j,i)+eps_dp)**(q_weert(j,i)*p_weert_inv(j,i))
              ! Basal sliding at pressure melting
#elif (SLIDE_LAW==2)
      cvxy1 = c_slide(j,i) &
              * ( (tau_b(j,i)+eps_dp)**(p_weert(j,i)-1) &
                  /(p_b_red_lim(j,i)+eps_dp)**q_weert(j,i) ) &
              * p_b(j,i)
      ctau1 = 1.0_dp/(c_slide(j,i)+eps_dp)**p_weert_inv(j,i) &
              * (p_b_red_lim(j,i)+eps_dp)**(q_weert(j,i)*p_weert_inv(j,i))
              ! Basal sliding at pressure melting
#elif (SLIDE_LAW==3)
      cvxy1 = c_slide(j,i) &
              * ( (tau_b(j,i)+eps_dp)**(p_weert(j,i)-1) &
                  /(p_b_red_lim(j,i)+eps_dp)**q_weert(j,i) ) &
              * p_b(j,i)
      ctau1 = 1.0_dp/(c_slide(j,i)+eps_dp)**p_weert_inv(j,i) &
              * (p_b_red(j,i)+eps_dp)**(q_weert(j,i)*p_weert_inv(j,i))
              ! Basal sliding at pressure melting
#else
      errormsg = ' >>> calc_vxy_b_sia: ' &
                    //'SLIDE_LAW must be 1, 2 or 3!'
      call error(errormsg)
#endif

!  ------ d_help_b, c_drag

      if (n_cts(j,i) == -1_i1b) then   ! cold ice base

         if (sub_melt_flag(j,i)) then
            temp_diff = max((temp_c_m(0,j,i)-temp_c(0,j,i)), 0.0_dp)
            cvxy1a    = exp(-gamma_slide_inv(j,i)*temp_diff)  ! sub-melt sliding
            ctau1a    = 1.0_dp/(cvxy1a+eps_dp)**p_weert_inv(j,i)
         else
            cvxy1a    = 0.0_dp   ! no sub-melt sliding
            ctau1a    = 1.0_dp/eps_dp**p_weert_inv(j,i)   ! dummy value
         end if

         d_help_b(j,i) = cvxy1*cvxy1a
         c_drag(j,i)   = ctau1*ctau1a

      else if (n_cts(j,i) == 0_i1b) then   ! temperate ice base

         d_help_b(j,i) = cvxy1   ! basal sliding
         c_drag(j,i)   = ctau1   ! (pressure-melting conditions)

      else   ! n_cts(j,i) == 1_i1b, temperate ice layer

         d_help_b(j,i) = cvxy1   ! basal sliding
         c_drag(j,i)   = ctau1   ! (pressure-melting conditions)

      end if

!    ---- Contribution of the basal water layer

#if (BASAL_HYDROLOGY==1)

#if (HYDRO_SLIDE_SAT_FCT==0)

      ratio_Hw_slide = max(H_w(j,i)*Hw0_slide_inv, 0.0_dp)
                        ! constrain to interval [0,infty)

      cvxy1b = 1.0_dp + c_Hw_slide*(1.0_dp-exp(-ratio_Hw_slide))
                        ! exponential saturation function
                        ! by Kleiner and Humbert (2014, J. Glaciol. 60)

#elif (HYDRO_SLIDE_SAT_FCT==1)

      ratio_Hw_slide = max(min(H_w(j,i)*Hw0_slide_inv, 1.0_dp), 0.0_dp)
                        ! constrain to interval [0,1]

      cvxy1b = 1.0_dp + c_Hw_slide*ratio_Hw_slide
                        ! linear saturation function

#elif (HYDRO_SLIDE_SAT_FCT==2)

      ratio_Hw_slide = max(min(H_w(j,i)*Hw0_slide_inv, 1.0_dp), 0.0_dp)
                        ! constrain to interval [0,1]

      cvxy1b = 1.0_dp + c_Hw_slide &
                           *ratio_Hw_slide*ratio_Hw_slide &
                              *(3.0_dp-2.0_dp*ratio_Hw_slide)
                        ! cubic S-shape saturation function

#elif (HYDRO_SLIDE_SAT_FCT==3)

      ratio_Hw_slide = max(min(H_w(j,i)*Hw0_slide_inv, 1.0_dp), 0.0_dp)
                        ! constrain to interval [0,1]

      cvxy1b = 1.0_dp + c_Hw_slide &
                           *ratio_Hw_slide*ratio_Hw_slide*ratio_Hw_slide &
                              *(10.0_dp + ratio_Hw_slide &
                                 *(-15.0_dp+6.0_dp*ratio_Hw_slide))
                        ! quintic S-shape saturation function

#else
      errormsg = ' >>> calc_vxy_b_sia: ' &
                    //'HYDRO_SLIDE_SAT_FCT must be 0, 1, 2 or 3!'
      call error(errormsg)
#endif

      ctau1b = 1.0_dp/(cvxy1b+eps_dp)**p_weert_inv(j,i)

      d_help_b(j,i) = d_help_b(j,i) *cvxy1b
      c_drag(j,i)   = c_drag(j,i)   *ctau1b

#endif

   else   ! maske(j,i) == 1_i1b, 2_i1b or 3_i1b away from the grounding line

      d_help_b(j,i) = 0.0_dp
      c_drag(j,i)   = 0.0_dp

   end if

end do
end do

!-------- Computation of vx_b (defined at (i+1/2,j)) --------

do i=0, IMAX-1
do j=1, JMAX-1
   vx_b(j,i) = -0.5_dp*(d_help_b(j,i)+d_help_b(j,i+1))*dzs_dx_aux(j,i)
end do
end do

!-------- Computation of vy_b (defined at (i,j+1/2)) --------

do i=1, IMAX-1
do j=0, JMAX-1
   vy_b(j,i) = -0.5_dp*(d_help_b(j,i)+d_help_b(j+1,i))*dzs_dy_aux(j,i)
end do
end do

!-------- Computation of vx_b_g and vy_b_g (defined at (i,j)) --------

do i=0, IMAX
do j=0, JMAX
   vx_b_g(j,i) = -d_help_b(j,i)*dzs_dxi_g(j,i)
   vy_b_g(j,i) = -d_help_b(j,i)*dzs_deta_g(j,i)
end do
end do

!-------- Limitation of computed vx_b, vy_b, vx_b_g, vy_b_g to the interval
!         [-VH_MAX, VH_MAX] --------

vh_max     = max(VH_MAX, eps_dp)/YEAR_SEC
vh_max_inv = 1.0_dp/vh_max

#if !defined(ALLOW_OPENAD) /* Normal */

call velocity_limiter_gradual(vx_b, vh_max, vh_max_inv)
call velocity_limiter_gradual(vy_b, vh_max, vh_max_inv)

call velocity_limiter_gradual(vx_b_g, vh_max, vh_max_inv)
call velocity_limiter_gradual(vy_b_g, vh_max, vh_max_inv)

#else /* OpenAD */

do i=0, IMAX
do j=0, JMAX

   call velocity_limiter_gradual(vx_b(j,i), vh_max, vh_max_inv)
   call velocity_limiter_gradual(vy_b(j,i), vh_max, vh_max_inv)
  
   call velocity_limiter_gradual(vx_b_g(j,i), vh_max, vh_max_inv)
   call velocity_limiter_gradual(vy_b_g(j,i), vh_max, vh_max_inv)

end do
end do

#endif /* Normal vs. OpenAD */

end subroutine calc_vxy_b_sia

!-------------------------------------------------------------------------------
!> Computation of the shear stresses txz, tyz, the effective shear stress
!! sigma, the depth-averaged fluidity flui_ave_sia, the horizontal
!! velocity vx, vy and the horizontal volume flux qx, qy in the shallow ice
!! approximation.
!<------------------------------------------------------------------------------
subroutine calc_vxy_sia(dzeta_c, dzeta_t)

use ice_material_properties_m, only : ratefac_c, ratefac_t, ratefac_c_t, creep

implicit none

real(dp), intent(in) :: dzeta_c, dzeta_t

integer(i4b) :: i, j, kc, kt
real(dp) :: avxy3(0:KCMAX), aqxy1(0:KCMAX)
real(dp) :: ctxyz1(0:KCMAX,0:JMAX,0:IMAX), &
            ctxyz2(0:KTMAX,0:JMAX,0:IMAX)
real(dp) :: flui_t(0:KTMAX), flui_c(0:KCMAX)
real(dp) :: cflui0(0:KTMAX), cflui1(0:KCMAX)
real(dp) :: cvxy2(0:KTMAX), cvxy3(0:KCMAX)
real(dp) :: cqxy0(0:KTMAX), cqxy1(0:KCMAX)
real(dp) :: vh_max, vh_max_inv
real(dp) :: flui_min, flui_max, flui_init
real(dp) :: ratio_sl_threshold

!-------- Term abbreviations --------

do kc=0, KCMAX
   if (flag_aa_nonzero) then
      avxy3(kc) = aa*eaz_c(kc)/(ea-1.0_dp)*dzeta_c
      aqxy1(kc) = aa/(ea-1.0_dp)*eaz_c(kc)*dzeta_c
   else
      avxy3(kc) = dzeta_c
      aqxy1(kc) = dzeta_c
   end if
end do

!-------- Computation of stresses --------

!  ------ Term abbreviations

do i=0, IMAX
do j=0, JMAX

   if ((maske(j,i) == 0_i1b).or.flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

      do kc=0, KCMAX
         ctxyz1(kc,j,i) = RHO*G*H_c(j,i)*(1.0_dp-eaz_c_quotient(kc))
      end do

      if (n_cts(j,i) == 1_i1b) then   ! temperate layer

         do kt=0, KTMAX
            ctxyz2(kt,j,i) = RHO*G*H_t(j,i)*(1.0_dp-zeta_t(kt))
         end do

      else   ! cold base (-1_i1b), temperate base (0_i1b)

         do kt=0, KTMAX
            ctxyz2(kt,j,i) = 0.0_dp
         end do

      end if

   else   ! maske(j,i) == 1_i1b, 2_i1b or 3_i1b away from the grounding line

      do kc=0, KCMAX
         ctxyz1(kc,j,i) = 0.0_dp
      end do

      do kt=0, KTMAX
         ctxyz2(kt,j,i) = 0.0_dp
      end do

   end if

end do
end do

!  ------ Shear stress txz (defined at (i+1/2,j,kc/t))

do i=0, IMAX-1
do j=0, JMAX

   do kc=0, KCMAX
      txz_c(kc,j,i) = -0.5_dp*(ctxyz1(kc,j,i)+ctxyz1(kc,j,i+1)) &
                      *dzs_dx_aux(j,i)
   end do

   do kt=0, KTMAX
      txz_t(kt,j,i) = txz_c(0,j,i) &
                      -0.5_dp*(ctxyz2(kt,j,i)+ctxyz2(kt,j,i+1)) &
                      *dzs_dx_aux(j,i)
   end do

end do
end do

!  ------ Shear stress tyz (defined at (i,j+1/2,kc/t))

do i=0, IMAX
do j=0, JMAX-1

   do kc=0, KCMAX
      tyz_c(kc,j,i) = -0.5_dp*(ctxyz1(kc,j,i)+ctxyz1(kc,j+1,i)) &
                      *dzs_dy_aux(j,i)
   end do

   do kt=0, KTMAX
      tyz_t(kt,j,i) = tyz_c(0,j,i) &
                      -0.5_dp*(ctxyz2(kt,j,i)+ctxyz2(kt,j+1,i)) &
                      *dzs_dy_aux(j,i)
   end do

end do
end do

!  ------ Effective shear stress sigma (defined at (i,j,kc/t))

do i=0, IMAX
do j=0, JMAX

   do kc=0, KCMAX

#if !defined(ALLOW_OPENAD) /* Normal */

      sigma_c(kc,j,i) = ctxyz1(kc,j,i) &
                        *sqrt(dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2)

#else /* OpenAD: guarding against non-differentiable sqrt(0) */

      if ( (dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2) > 0 ) then
         sigma_c(kc,j,i) = ctxyz1(kc,j,i) &
                           *sqrt(dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2)
      else
         sigma_c(kc,j,i) = 0.0_dp 
      end if

#endif /* Normal vs. OpenAD */

   end do

   do kt=0, KTMAX

#if !defined(ALLOW_OPENAD) /* Normal */

      sigma_t(kt,j,i) = sigma_c(0,j,i) &
                        + ctxyz2(kt,j,i) &
                          *sqrt(dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2)

#else /* OpenAD: guarding against non-differentiable sqrt(0) */

      if ( (dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2) > 0 ) then
         sigma_t(kt,j,i) = sigma_c(0,j,i) &
                           + ctxyz2(kt,j,i) &
                             *sqrt(dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2)
      else
         sigma_t(kt,j,i) = sigma_c(0,j,i) 
      end if 

#endif /* Normal vs. OpenAD */

   end do

end do
end do

!-------- Computation of the depth-averaged fluidity
!                 (defined on the grid points (i,j)) --------

#if (defined(VISC_MIN) && defined(VISC_MAX))
  flui_min = 1.0_dp/VISC_MAX
  flui_max = 1.0_dp/VISC_MIN
#else
  flui_min = 1.0e-25_dp   ! 1/(Pa s)
  flui_max = 1.0e-10_dp   ! 1/(Pa s)
#endif

#if (defined(VISC_INIT_SSA))
  flui_init = 1.0_dp/VISC_INIT_SSA
#else
  flui_init = 1.0e-15_dp   ! 1/(Pa s)
#endif

do i=0, IMAX
do j=0, JMAX

   if ((maske(j,i) == 0_i1b).or.flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

!  ------ Fluidity, abbreviations

      do kt=0, KTMAX
         flui_t(kt) = 2.0_dp &
                      *enh_t(kt,j,i) &
                      *ratefac_t(omega_t(kt,j,i)) &
                      *creep(sigma_t(kt,j,i))
         cflui0(kt) = H_t(j,i)*flui_t(kt)*dzeta_t
      end do

      do kc=0, KCMAX
         flui_c(kc) = 2.0_dp &
                      *enh_c(kc,j,i) &
#if (CALCMOD==0 || CALCMOD==1 || CALCMOD==-1)
                      *ratefac_c(temp_c(kc,j,i), temp_c_m(kc,j,i)) &
#elif (CALCMOD==2 || CALCMOD==3)
                      *ratefac_c_t(temp_c(kc,j,i), omega_c(kc,j,i), &
                                                   temp_c_m(kc,j,i)) &
#endif
                      *creep(sigma_c(kc,j,i))
         cflui1(kc) = aqxy1(kc)*H_c(j,i)*flui_c(kc)
      end do

!  ------ Depth average

      flui_ave_sia(j,i) = 0.0_dp

      if (n_cts(j,i) == 1_i1b) then

         do kt=0, KTMAX-1
            flui_ave_sia(j,i) = flui_ave_sia(j,i)+0.5_dp*(cflui0(kt+1)+cflui0(kt))
         end do

      end if

      do kc=0, KCMAX-1
         flui_ave_sia(j,i) = flui_ave_sia(j,i)+0.5_dp*(cflui1(kc+1)+cflui1(kc))
      end do

      flui_ave_sia(j,i) = flui_ave_sia(j,i)/max((H_c(j,i)+H_t(j,i)), eps_dp)

      flui_ave_sia(j,i) = max(min(flui_ave_sia(j,i), flui_max), flui_min)

   else   ! maske(j,i) == 1_i1b, 2_i1b or 3_i1b away from the grounding line

      flui_ave_sia(j,i) = flui_init   ! dummy value

   end if

end do
end do

!-------- Computation of d_help_c/t
!         (defined on the grid points (i,j,kc/t)) --------

do i=0, IMAX
do j=0, JMAX

   if ((maske(j,i) == 0_i1b).or.flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

!  ------ Abbreviations

      do kt=0, KTMAX
         cvxy2(kt) = 2.0_dp*H_t(j,i) &
                     *enh_t(kt,j,i) &
                     *ratefac_t(omega_t(kt,j,i)) &
                     *creep(sigma_t(kt,j,i)) &
                     *(ctxyz1(0,j,i)+ctxyz2(kt,j,i)) &
                     *dzeta_t
      end do

      do kc=0, KCMAX
         cvxy3(kc) = 2.0_dp*avxy3(kc)*H_c(j,i) &
                     *enh_c(kc,j,i) &
#if (CALCMOD==0 || CALCMOD==1 || CALCMOD==-1)
                     *ratefac_c(temp_c(kc,j,i), temp_c_m(kc,j,i)) &
#elif (CALCMOD==2 || CALCMOD==3)
                     *ratefac_c_t(temp_c(kc,j,i), omega_c(kc,j,i), &
                                                  temp_c_m(kc,j,i)) &
#endif
                     *creep(sigma_c(kc,j,i)) &
                     *ctxyz1(kc,j,i)
      end do

!  ------ d_help_c, d_help_t

      if (n_cts(j,i) == -1_i1b) then   ! cold ice base

         do kt=0, KTMAX
            d_help_t(kt,j,i) = d_help_b(j,i)
         end do

         d_help_c(0,j,i) = d_help_t(KTMAX,j,i)

         do kc=0, KCMAX-1
            d_help_c(kc+1,j,i) = d_help_c(kc,j,i) &
                                +0.5_dp*(cvxy3(kc+1)+cvxy3(kc))
         end do

      else if (n_cts(j,i) == 0_i1b) then   ! temperate ice base

         do kt=0, KTMAX
            d_help_t(kt,j,i) = d_help_b(j,i)
         end do

         d_help_c(0,j,i) = d_help_t(KTMAX,j,i)

         do kc=0, KCMAX-1
            d_help_c(kc+1,j,i) = d_help_c(kc,j,i) &
                                +0.5_dp*(cvxy3(kc+1)+cvxy3(kc))
         end do

      else   ! n_cts(j,i) == 1_i1b, temperate ice layer

         d_help_t(0,j,i) = d_help_b(j,i)

         do kt=0, KTMAX-1
            d_help_t(kt+1,j,i) = d_help_t(kt,j,i) &
                                +0.5_dp*(cvxy2(kt+1)+cvxy2(kt))
         end do

         d_help_c(0,j,i) = d_help_t(KTMAX,j,i)

         do kc=0, KCMAX-1
            d_help_c(kc+1,j,i) = d_help_c(kc,j,i) &
                                +0.5_dp*(cvxy3(kc+1)+cvxy3(kc))
         end do

      end if

   else   ! maske(j,i) == 1_i1b, 2_i1b or 3_i1b away from the grounding line

      do kt=0, KTMAX
         d_help_t(kt,j,i) = 0.0_dp
      end do

      do kc=0, KCMAX
         d_help_c(kc,j,i) = 0.0_dp
      end do

   end if

end do
end do

!-------- Computation of vx_c/t (defined at (i+1/2,j,kc/t)) --------

do i=0, IMAX-1
do j=1, JMAX-1

   do kt=0, KTMAX
      vx_t(kt,j,i) = -0.5_dp*(d_help_t(kt,j,i)+d_help_t(kt,j,i+1)) &
                     *dzs_dx_aux(j,i)
   end do

   do kc=0, KCMAX
      vx_c(kc,j,i) = -0.5_dp*(d_help_c(kc,j,i)+d_help_c(kc,j,i+1)) &
                     *dzs_dx_aux(j,i)
   end do

end do
end do

!-------- Computation of vy_c/t (defined at (i,j+1/2,kc/t)) --------

do i=1, IMAX-1
do j=0, JMAX-1

   do kt=0, KTMAX
      vy_t(kt,j,i) = -0.5_dp*(d_help_t(kt,j,i)+d_help_t(kt,j+1,i)) &
                     *dzs_dy_aux(j,i)
   end do

   do kc=0, KCMAX
      vy_c(kc,j,i) = -0.5_dp*(d_help_c(kc,j,i)+d_help_c(kc,j+1,i)) &
                     *dzs_dy_aux(j,i)
   end do

end do
end do

!-------- Computation of the surface velocities vx_s_g and vy_s_g
!                                              (defined at (i,j)) --------

do i=0, IMAX
do j=0, JMAX
   vx_s_g(j,i) = -d_help_c(KCMAX,j,i)*dzs_dxi_g(j,i)
   vy_s_g(j,i) = -d_help_c(KCMAX,j,i)*dzs_deta_g(j,i)
end do
end do

!-------- Limitation of computed vx_c/t, vy_c/t, vx_s_g, vy_s_g
!         to the interval [-VH_MAX, VH_MAX] --------

vh_max     = max(VH_MAX, eps_dp)/YEAR_SEC
vh_max_inv = 1.0_dp/vh_max

#if !defined (ALLOW_OPENAD) /* Normal */

call velocity_limiter_gradual(vx_s_g, vh_max, vh_max_inv)
call velocity_limiter_gradual(vy_s_g, vh_max, vh_max_inv)

call velocity_limiter_gradual(vx_t, vh_max, vh_max_inv)
call velocity_limiter_gradual(vy_t, vh_max, vh_max_inv)
call velocity_limiter_gradual(vx_c, vh_max, vh_max_inv)
call velocity_limiter_gradual(vy_c, vh_max, vh_max_inv)

#else /* OpenAD */

do i=0, IMAX
do j=0, JMAX

   call velocity_limiter_gradual(vx_s_g(j,i), vh_max, vh_max_inv)
   call velocity_limiter_gradual(vy_s_g(j,i), vh_max, vh_max_inv)

   do kt=0, KTMAX
      call velocity_limiter_gradual(vx_t(kt,j,i), vh_max, vh_max_inv)
      call velocity_limiter_gradual(vy_t(kt,j,i), vh_max, vh_max_inv)
   end do

   do kc=0, KCMAX
      call velocity_limiter_gradual(vx_c(kc,j,i), vh_max, vh_max_inv)
      call velocity_limiter_gradual(vy_c(kc,j,i), vh_max, vh_max_inv)
   end do

end do
end do

#endif /* Normal vs. OpenAD */

!-------- Computation of h_diff
!         (defined on the grid points (i,j)) --------

do i=0, IMAX
do j=0, JMAX

   if ((maske(j,i) == 0_i1b).or.flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

!  ------ Abbreviations

      do kt=0, KTMAX
         cqxy0(kt) = H_t(j,i)*d_help_t(kt,j,i)*dzeta_t
      end do

      do kc=0, KCMAX
         cqxy1(kc) = aqxy1(kc)*H_c(j,i)*d_help_c(kc,j,i)
      end do

!  ------ h_diff

      h_diff(j,i) = 0.0_dp

      if (n_cts(j,i) == 1_i1b) then

         do kt=0, KTMAX-1
            h_diff(j,i) = h_diff(j,i)+0.5_dp*(cqxy0(kt+1)+cqxy0(kt))
         end do

      end if

      do kc=0, KCMAX-1
         h_diff(j,i) = h_diff(j,i)+0.5_dp*(cqxy1(kc+1)+cqxy1(kc))
      end do

!  ------ Limitation of h_diff

      if (h_diff(j,i) < HD_MIN) h_diff(j,i) = 0.0_dp
      if (h_diff(j,i) > HD_MAX) h_diff(j,i) = HD_MAX

   else   ! maske(j,i) == 1_i1b, 2_i1b or 3_i1b away from the grounding line

      h_diff(j,i) = 0.0_dp

   end if

end do
end do

!-------- Computation of the horizontal volume flux
!                            and the depth-averaged velocity --------

do i=0, IMAX-1
do j=0, JMAX

   qx(j,i) = -0.5_dp*(h_diff(j,i)+h_diff(j,i+1))*dzs_dx_aux(j,i)

   if ( (maske(j,i)==0_i1b).or.(maske(j,i+1)==0_i1b) ) then
                               ! at least one neighbour point is grounded ice

      vx_m(j,i) = qx(j,i) &
                    / ( 0.5_dp &
                      * ( (H_c(j,i)+H_t(j,i)) + (H_c(j,i+1)+H_t(j,i+1)) ) )

      call velocity_limiter_gradual(vx_m(j,i), vh_max, vh_max_inv)

      ratio_sl_x(j,i) = abs(vx_t(0,j,i)) / max(abs(vx_c(KCMAX,j,i)), eps_dp)

   else 

      vx_m(j,i)       = 0.0_dp
      ratio_sl_x(j,i) = 0.0_dp

   end if

end do
end do

do i=0, IMAX
do j=0, JMAX-1

   qy(j,i) = -0.5_dp*(h_diff(j,i)+h_diff(j+1,i))*dzs_dy_aux(j,i)

   if ( (maske(j,i)==0_i1b).or.(maske(j+1,i)==0_i1b) ) then
                               ! at least one neighbour point is grounded ice

      vy_m(j,i) = qy(j,i) &
                    / ( 0.5_dp &
                      * ( (H_c(j,i)+H_t(j,i)) + (H_c(j+1,i)+H_t(j+1,i)) ) )

      call velocity_limiter_gradual(vy_m(j,i), vh_max, vh_max_inv)

      ratio_sl_y(j,i) = abs(vy_t(0,j,i)) / max(abs(vy_c(KCMAX,j,i)), eps_dp)

   else 

      vy_m(j,i)       = 0.0_dp
      ratio_sl_y(j,i) = 0.0_dp

   end if

end do
end do

!-------- Detection of shelfy stream points --------

flag_shelfy_stream_x = .false.
flag_shelfy_stream_y = .false.
flag_shelfy_stream   = .false.

#if (DYNAMICS==0 || DYNAMICS==1)

ratio_sl_threshold = 1.11e+11_dp   ! dummy value

#elif (DYNAMICS==2)

#if ( defined(RATIO_SL_THRESH) )
ratio_sl_threshold = RATIO_SL_THRESH
#else
ratio_sl_threshold = 0.5_dp   ! default value
#endif

do i=0, IMAX-1
do j=0, JMAX
   if (ratio_sl_x(j,i) > ratio_sl_threshold) flag_shelfy_stream_x(j,i) = .true.
end do
end do

do i=0, IMAX
do j=0, JMAX-1
   if (ratio_sl_y(j,i) > ratio_sl_threshold) flag_shelfy_stream_y(j,i) = .true.
end do
end do

do i=1, IMAX-1
do j=1, JMAX-1

   if ( (maske(j,i) == 0_i1b) &   ! grounded ice
        .and. &
        (     flag_shelfy_stream_x(j,i-1)   &   ! at least
          .or.flag_shelfy_stream_x(j,i)     &   ! one neighbour
          .or.flag_shelfy_stream_y(j-1,i)   &   ! on the staggered grid
          .or.flag_shelfy_stream_y(j,i)   ) &   ! is a shelfy stream point
      ) then

      flag_shelfy_stream(j,i) = .true.

   end if

end do
end do

#else

errormsg = ' >>> calc_vxy_sia: DYNAMICS must be 0, 1 or 2!'
call error(errormsg)

#endif

!-------- Save mean (depth-averaged) horizontal velocities from SIA --------

vx_m_sia = vx_m
vy_m_sia = vy_m

!-------- Initialisation of the variable q_gl_g
!         (volume flux across the grounding line, to be
!         computed in the routine calc_vxy_ssa
!         if ice shelves are present)

q_gl_g = 0.0_dp

end subroutine calc_vxy_sia

!-------------------------------------------------------------------------------
!> Computation of the horizontal velocity vx, vy, the horizontal volume flux
!> qx, qy etc. for static ice.
!<------------------------------------------------------------------------------
subroutine calc_vxy_static()

implicit none

real(dp) :: flui_init

#if (defined(VISC_INIT_SSA))
  flui_init = 1.0_dp/VISC_INIT_SSA
#else
  flui_init = 1.0e-15_dp   ! 1/(Pa s)
#endif

c_slide = 0.0_dp
p_weert = 0
q_weert = 0
p_b_w   = 0.0_dp

d_help_b = 0.0_dp
c_drag   = 0.0_dp

vx_b   = 0.0_dp
vy_b   = 0.0_dp
vx_b_g = 0.0_dp
vy_b_g = 0.0_dp

txz_c = 0.0_dp
txz_t = 0.0_dp

tyz_c = 0.0_dp
tyz_t = 0.0_dp

sigma_c = 0.0_dp
sigma_t = 0.0_dp

flui_ave_sia = flui_init
de_ssa       = 0.0_dp
vis_int_g    = 0.0_dp

d_help_t = 0.0_dp
d_help_c = 0.0_dp

vx_c = 0.0_dp
vy_c = 0.0_dp

vx_t = 0.0_dp
vy_t = 0.0_dp

vx_s_g = 0.0_dp
vy_s_g = 0.0_dp

h_diff = 0.0_dp

qx = 0.0_dp
qy = 0.0_dp

vx_m = 0.0_dp
vy_m = 0.0_dp

ratio_sl_x = 0.0_dp
ratio_sl_y = 0.0_dp

flag_shelfy_stream_x = .false.
flag_shelfy_stream_y = .false.
flag_shelfy_stream   = .false.

vx_m_sia = 0.0_dp
vy_m_sia = 0.0_dp

q_gl_g = 0.0_dp

end subroutine calc_vxy_static

!-------------------------------------------------------------------------------
!> Computation of the horizontal velocity vx, vy, the horizontal volume flux
!! qx, qy and the flux across the grounding line q_gl_g in the shallow shelf
!! approximation (SSA) or the shelfy stream approximation (SStA).
!<------------------------------------------------------------------------------
subroutine calc_vxy_ssa(z_sl, dxi, deta, dzeta_c, dzeta_t)

implicit none

real(dp), intent(in) :: z_sl, dxi, deta, dzeta_c, dzeta_t

integer(i4b) :: i, j, kc, kt, m
integer(i4b) :: iter_ssa
real(dp), dimension(0:JMAX,0:IMAX) :: vx_m_prev, vy_m_prev
real(dp) :: tol_ssa, rel_ssa
real(dp) :: res_vxy_m_ssa_1, res_vxy_m_ssa_2, res_vxy_m_ssa
real(dp) :: dxi_inv, deta_inv
real(dp) :: visc_init
real(dp) :: vh_max, vh_max_inv
real(dp) :: ratio_sl_threshold, ratio_help
real(dp), dimension(0:JMAX,0:IMAX) :: weigh_ssta_sia_x, weigh_ssta_sia_y
real(dp) :: qx_gl_g, qy_gl_g
logical, dimension(0:JMAX,0:IMAX) :: flag_calc_vxy_ssa_x, flag_calc_vxy_ssa_y

#if (MARGIN==3 || DYNAMICS==2)

!-------- Parameters for the relaxation scheme --------

#if (defined(TOL_ITER_SSA))
   tol_ssa = TOL_ITER_SSA   ! tolerance of iterations 
#else
   tol_ssa = 0.1_dp         ! default value
#endif

#if (defined(N_ITER_SSA))
   iter_ssa = max(N_ITER_SSA, 1)   ! max. number of iterations
#else
   iter_ssa = 3                    ! default value
#endif

#if (defined(RELAX_FACT_SSA))
   rel_ssa = RELAX_FACT_SSA   ! relaxation factor
#else
   rel_ssa = 0.7_dp           ! default value
#endif

#if (defined(VISC_INIT_SSA))
  visc_init = VISC_INIT_SSA
#else
  visc_init = 1.0e+15_dp   ! Pa s
#endif

vh_max     = max(VH_MAX, eps_dp)/YEAR_SEC
vh_max_inv = 1.0_dp/vh_max

write(6,'(10x,a)') 'calc_vxy_ssa:'

!-------- Iterations --------

res_vxy_m_ssa = 1.11e+11_dp   ! initial, very large value of the residual

m=0

do while ((m < iter_ssa).and.(res_vxy_m_ssa > tol_ssa))

   m = m+1

   write(6,'(13x,a,i0,a)', advance='no') 'Iter ', m, ': '

!  ------ Save velocities from previous iteration

   vx_m_prev = vx_m_ssa
   vy_m_prev = vy_m_ssa

!  ------ Depth-integrated viscosity vis_int_g

   if (m > 1) then

      call calc_vis_ssa(dxi, deta, dzeta_c, dzeta_t)

   else   ! (m == 1, first iteration)

#if (!defined(ITER_INIT_SSA) || ITER_INIT_SSA==1)
      vis_int_g = (H_c+H_t)*visc_init
                  ! constant viscosity times ice thickness
#elif (ITER_INIT_SSA==2)
      vis_int_g = (H_c+H_t)*vis_ave_g
                  ! previous depth-averaged viscosity times
                  ! ice thickness
#elif (ITER_INIT_SSA==3)
      call calc_vis_ssa(dxi, deta, dzeta_c, dzeta_t)
                  ! standard computation by subroutine
                  ! calc_vis_ssa
#else
      errormsg = ' >>> calc_vxy_ssa: ITER_INIT_SSA must be 1, 2 or 3!'
      call error(errormsg)
#endif

   end if

!  ------ Horizontal velocity vx_m_ssa, vy_m_ssa

   flag_calc_vxy_ssa_x = .false.   ! initialization
   flag_calc_vxy_ssa_y = .false.   ! initialization

   call calc_vxy_ssa_matrix(z_sl, dxi, deta, &
                            flag_calc_vxy_ssa_x, flag_calc_vxy_ssa_y)

#if !defined(ALLOW_OPENAD) /* Normal */

   call velocity_limiter_gradual(vx_m_ssa, vh_max, vh_max_inv)
   call velocity_limiter_gradual(vy_m_ssa, vh_max, vh_max_inv)

#else /* OpenAD */

   do i=0, IMAX
   do j=0, JMAX
      call velocity_limiter_gradual(vx_m_ssa(j,i), vh_max, vh_max_inv)
      call velocity_limiter_gradual(vy_m_ssa(j,i), vh_max, vh_max_inv)
   end do
   end do

#endif /* Normal vs. OpenAD */

!  ------ Relaxation scheme

   if (m > 1) then
      vx_m_ssa = rel_ssa*vx_m_ssa + (1.0_dp-rel_ssa)*vx_m_prev
      vy_m_ssa = rel_ssa*vy_m_ssa + (1.0_dp-rel_ssa)*vy_m_prev
   end if

!  ------ Residual

   res_vxy_m_ssa_1 &
      = sqrt( sum((vx_m_ssa-vx_m_prev)*(vx_m_ssa-vx_m_prev), &
                  mask=flag_calc_vxy_ssa_x) &
             +sum((vy_m_ssa-vy_m_prev)*(vy_m_ssa-vy_m_prev), &
                  mask=flag_calc_vxy_ssa_y) )
   res_vxy_m_ssa_2 &
      = sqrt( sum((vx_m_ssa+vx_m_prev)*(vx_m_ssa+vx_m_prev), &
                  mask=flag_calc_vxy_ssa_x) &
             +sum((vy_m_ssa+vy_m_prev)*(vy_m_ssa+vy_m_prev), &
                  mask=flag_calc_vxy_ssa_y) )

   res_vxy_m_ssa = 2.0_dp*res_vxy_m_ssa_1/max(res_vxy_m_ssa_2, eps_dp)

   write(6,'(a,es9.2)') 'res =', res_vxy_m_ssa

end do

!  ------ Depth-integrated viscosity vis_int_g

call calc_vis_ssa(dxi, deta, dzeta_c, dzeta_t)

!-------- 3-D velocities, basal velocities and volume flux --------

#if (DYNAMICS==0 || DYNAMICS==1)

ratio_sl_threshold = 1.11e+11_dp   ! dummy value
ratio_help         = 0.0_dp

#elif (DYNAMICS==2)

#if ( defined(RATIO_SL_THRESH) )
ratio_sl_threshold = RATIO_SL_THRESH
#else
ratio_sl_threshold = 0.5_dp   ! default value
#endif

ratio_help = 1.0_dp/(1.0_dp-ratio_sl_threshold)

#else

errormsg = ' >>> calc_vxy_ssa: DYNAMICS must be 0, 1 or 2!'
call error(errormsg)

#endif

weigh_ssta_sia_x = 0.0_dp
weigh_ssta_sia_y = 0.0_dp

!  ------ x-component

do i=0, IMAX-1
do j=0, JMAX

   if (flag_shelfy_stream_x(j,i)) then   ! shelfy stream

      weigh_ssta_sia_x(j,i) = (ratio_sl_x(j,i)-ratio_sl_threshold)*ratio_help

      weigh_ssta_sia_x(j,i) = max(min(weigh_ssta_sia_x(j,i), 1.0_dp), 0.0_dp)
                              ! constrain to interval [0,1]

#if (SSTA_SIA_WEIGH_FCT==0)

      ! stick to the linear function set above

#elif (SSTA_SIA_WEIGH_FCT==1)

      weigh_ssta_sia_x(j,i) = weigh_ssta_sia_x(j,i)*weigh_ssta_sia_x(j,i) &
                                   *(3.0_dp-2.0_dp*weigh_ssta_sia_x(j,i))
                              ! make transition smooth (cubic function)

#elif (SSTA_SIA_WEIGH_FCT==2)

      weigh_ssta_sia_x(j,i) = weigh_ssta_sia_x(j,i)*weigh_ssta_sia_x(j,i) &
                                                   *weigh_ssta_sia_x(j,i) &
                                   *(10.0_dp + weigh_ssta_sia_x(j,i) &
                                       *(-15.0_dp+6.0_dp*weigh_ssta_sia_x(j,i)))
                              ! make transition even smoother (quintic function)

#else
      errormsg = ' >>> calc_vxy_ssa: SSTA_SIA_WEIGH_FCT must be 0, 1 or 2!'
      call error(errormsg)
#endif

      do kt=0, KTMAX
         vx_t(kt,j,i) = weigh_ssta_sia_x(j,i)*vx_m_ssa(j,i) &
                        + (1.0_dp-weigh_ssta_sia_x(j,i))*vx_t(kt,j,i)
      end do

      do kc=0, KCMAX
         vx_c(kc,j,i) = weigh_ssta_sia_x(j,i)*vx_m_ssa(j,i) &
                        + (1.0_dp-weigh_ssta_sia_x(j,i))*vx_c(kc,j,i)
      end do

      vx_b(j,i) = vx_t(0,j,i)

      vx_m(j,i) = weigh_ssta_sia_x(j,i)*vx_m_ssa(j,i) &
                  + (1.0_dp-weigh_ssta_sia_x(j,i))*vx_m_sia(j,i)

      qx(j,i)   = vx_m(j,i) &
                     * 0.5_dp * ( (H_c(j,i)+H_t(j,i))+(H_c(j,i+1)+H_t(j,i+1)) )

   else if (flag_calc_vxy_ssa_x(j,i)) then   ! floating ice

      do kt=0, KTMAX
         vx_t(kt,j,i) = vx_m_ssa(j,i)
      end do

      do kc=0, KCMAX
         vx_c(kc,j,i) = vx_m_ssa(j,i)
      end do

      vx_b(j,i) = vx_m_ssa(j,i)

      vx_m(j,i) = vx_m_ssa(j,i)

      qx(j,i)   = vx_m(j,i) &
                     * 0.5_dp * ( (H_c(j,i)+H_t(j,i))+(H_c(j,i+1)+H_t(j,i+1)) )

!  else
!     In all other cases, the depth-averaged velocities vx_m_ssa(j,i) computed
!     by the SSA/SStA solver are discarded.

   end if

end do
end do

!  ------ y-component

do i=0, IMAX
do j=0, JMAX-1

   if (flag_shelfy_stream_y(j,i)) then   ! shelfy stream

      weigh_ssta_sia_y(j,i) = (ratio_sl_y(j,i)-ratio_sl_threshold)*ratio_help

      weigh_ssta_sia_y(j,i) = max(min(weigh_ssta_sia_y(j,i), 1.0_dp), 0.0_dp)
                              ! constrain to interval [0,1]

#if (SSTA_SIA_WEIGH_FCT==0)

      ! stick to the linear function set above

#elif (SSTA_SIA_WEIGH_FCT==1)

      weigh_ssta_sia_y(j,i) = weigh_ssta_sia_y(j,i)*weigh_ssta_sia_y(j,i) &
                                   *(3.0_dp-2.0_dp*weigh_ssta_sia_y(j,i))
                              ! make transition smooth (cubic function)

#elif (SSTA_SIA_WEIGH_FCT==2)

      weigh_ssta_sia_y(j,i) = weigh_ssta_sia_y(j,i)*weigh_ssta_sia_y(j,i) &
                                                   *weigh_ssta_sia_y(j,i) &
                                   *(10.0_dp + weigh_ssta_sia_y(j,i) &
                                       *(-15.0_dp+6.0_dp*weigh_ssta_sia_y(j,i)))
                              ! make transition even smoother (quintic function)

#else
      errormsg = ' >>> calc_vxy_ssa: SSTA_SIA_WEIGH_FCT must be 0, 1 or 2!'
      call error(errormsg)
#endif

      do kt=0, KTMAX
         vy_t(kt,j,i) = weigh_ssta_sia_y(j,i)*vy_m_ssa(j,i) &
                        + (1.0_dp-weigh_ssta_sia_y(j,i))*vy_t(kt,j,i)
      end do

      do kc=0, KCMAX
         vy_c(kc,j,i) = weigh_ssta_sia_y(j,i)*vy_m_ssa(j,i) &
                        + (1.0_dp-weigh_ssta_sia_y(j,i))*vy_c(kc,j,i)
      end do

      vy_b(j,i) = vy_t(0,j,i)

      vy_m(j,i) = weigh_ssta_sia_y(j,i)*vy_m_ssa(j,i) &
                  + (1.0_dp-weigh_ssta_sia_y(j,i))*vy_m_sia(j,i)

      qy(j,i)   = vy_m(j,i) &
                     * 0.5_dp * ( (H_c(j,i)+H_t(j,i))+(H_c(j+1,i)+H_t(j+1,i)) )

   else if (flag_calc_vxy_ssa_y(j,i)) then   ! floating ice

      do kt=0, KTMAX
         vy_t(kt,j,i) = vy_m_ssa(j,i)
      end do

      do kc=0, KCMAX
         vy_c(kc,j,i) = vy_m_ssa(j,i)
      end do

      vy_b(j,i) = vy_m_ssa(j,i)

      vy_m(j,i) = vy_m_ssa(j,i)

      qy(j,i)   = vy_m(j,i) &
                     * 0.5_dp * ( (H_c(j,i)+H_t(j,i))+(H_c(j+1,i)+H_t(j+1,i)) )

!  else
!     In all other cases, the depth-averaged velocities vy_m_ssa(j,i) computed
!     by the SSA/SStA solver are discarded.

   end if

end do
end do

!-------- Surface and basal velocities vx_s_g vy_s_g, vx_b_g vy_b_g
!                                                (defined at (i,j)) --------

do i=1, IMAX-1
do j=1, JMAX-1

   if (flag_shelfy_stream(j,i)) then   ! shelfy stream

      vx_s_g(j,i) = 0.5_dp*(vx_c(KCMAX,j,i-1)+vx_c(KCMAX,j,i))
      vx_b_g(j,i) = 0.5_dp*(vx_b(      j,i-1)+vx_b(      j,i))

      vy_s_g(j,i) = 0.5_dp*(vy_c(KCMAX,j-1,i)+vy_c(KCMAX,j,i))
      vy_b_g(j,i) = 0.5_dp*(vy_b(      j-1,i)+vy_b(      j,i))

   else if (maske(j,i)==3_i1b) then   ! floating ice

      vx_s_g(j,i) = 0.5_dp*(vx_m(j,i-1)+vx_m(j,i))
      vx_b_g(j,i) = vx_s_g(j,i)

      vy_s_g(j,i) = 0.5_dp*(vy_m(j-1,i)+vy_m(j,i))
      vy_b_g(j,i) = vy_s_g(j,i)

   end if

end do
end do

!-------- Computation of the flux across the grounding line q_gl_g

do i=1, IMAX-1
do j=1, JMAX-1

   if ( flag_grounding_line_1(j,i) ) then   ! grounding line

      qx_gl_g = 0.5_dp*(qx(j,i)+qx(j,i-1))
      qy_gl_g = 0.5_dp*(qy(j,i)+qy(j-1,i))

#if !defined(ALLOW_OPENAD) /* Normal */

      q_gl_g(j,i) = sqrt(qx_gl_g*qx_gl_g+qy_gl_g*qy_gl_g)

#else /* OpenAD: guarding against non-differentiable sqrt(0) */

      if ( (qx_gl_g*qx_gl_g+qy_gl_g*qy_gl_g) > 0 ) then
         q_gl_g(j,i) = sqrt(qx_gl_g*qx_gl_g+qy_gl_g*qy_gl_g)
      else
         q_gl_g(j,i) = 0.0_dp 
      end if

#endif /* Normal vs. OpenAD */

   end if

end do
end do

#else

errormsg = ' >>> calc_vxy_ssa: Only to be called for MARGIN==3 or DYNAMICS==2!'
call error(errormsg)

#endif

end subroutine calc_vxy_ssa

!-------------------------------------------------------------------------------
!> Solution of the system of linear equations for the horizontal velocities
!! vx_m_ssa, vy_m_ssa in the shallow shelf approximation.
!<------------------------------------------------------------------------------
subroutine calc_vxy_ssa_matrix(z_sl, dxi, deta, &
                               flag_calc_vxy_ssa_x, flag_calc_vxy_ssa_y)

#if (MARGIN==3 || DYNAMICS==2) 
#if defined(ALLOW_OPENAD) /* OpenAD */
use sico_maths_m
#endif /* OpenAD */
#endif

implicit none

real(dp), intent(in) :: dxi, deta, z_sl
logical, dimension(0:JMAX,0:IMAX), intent(inout) :: &
                            flag_calc_vxy_ssa_x, flag_calc_vxy_ssa_y

integer(i4b) :: i, j, k, n
integer(i4b) :: i1, j1
real(dp) :: inv_dxi, inv_deta, inv_dxi_deta, inv_dxi2, inv_deta2
real(dp) :: factor_rhs_1, factor_rhs_2, factor_rhs_3a, factor_rhs_3b
real(dp) :: rhosw_rho_ratio
real(dp) :: H_mid, zl_mid
real(dp), dimension(0:JMAX,0:IMAX) :: vis_int_sgxy, beta_drag
character(len=256) :: ch_solver_set_option

#if (MARGIN==3 || DYNAMICS==2)

#if !defined(ALLOW_OPENAD) /* Normal */

LIS_INTEGER :: ierr
LIS_INTEGER :: nc, nr
LIS_INTEGER :: lin_iter
LIS_MATRIX  :: lgs_a
LIS_VECTOR  :: lgs_b, lgs_x
LIS_SOLVER  :: solver

LIS_INTEGER, parameter                 :: nmax   =  2*(IMAX+1)*(JMAX+1)
LIS_INTEGER, parameter                 :: n_sprs = 20*(IMAX+1)*(JMAX+1)
LIS_INTEGER, allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
LIS_SCALAR,  allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value

#else /* OpenAD */

integer(i4b) :: ierr
integer(i4b) :: nc, nr
integer(i4b) :: lin_iter

integer(i4b), parameter                 :: nmax   =  2*(IMAX+1)*(JMAX+1)
integer(i4b), parameter                 :: n_sprs = 20*(IMAX+1)*(JMAX+1)
integer(i4b), dimension(nmax+1)         :: lgs_a_ptr
real(dp), dimension(n_sprs)             :: lgs_a_value
real(dp), dimension(n_sprs)             :: lgs_a_index
integer(i4b), dimension(n_sprs)         :: lgs_a_index_pass
real(dp), dimension(nmax)               :: lgs_b_value, lgs_x_value

#endif /* Normal vs. OpenAD */

!-------- Abbreviations --------

inv_dxi      = 1.0_dp/dxi
inv_deta     = 1.0_dp/deta
inv_dxi_deta = 1.0_dp/(dxi*deta)
inv_dxi2     = 1.0_dp/(dxi*dxi)
inv_deta2    = 1.0_dp/(deta*deta)

factor_rhs_1  = RHO*G
factor_rhs_2  = 0.5_dp*RHO*G*(RHO_SW-RHO)/RHO_SW
factor_rhs_3a = 0.5_dp*RHO*G
factor_rhs_3b = 0.5_dp*RHO_SW*G

rhosw_rho_ratio = RHO_SW/RHO

!-------- Depth-integrated viscosity on the staggered grid
!                                       [at (i+1/2,j+1/2)] --------

vis_int_sgxy = 0.0_dp   ! initialisation

do i=0, IMAX-1
do j=0, JMAX-1

   k=0

   if ((maske(j,i)==0_i1b).or.(maske(j,i)==3_i1b)) then
      k = k+1                              ! floating or grounded ice
      vis_int_sgxy(j,i) = vis_int_sgxy(j,i) + vis_int_g(j,i)
   end if

   if ((maske(j,i+1)==0_i1b).or.(maske(j,i+1)==3_i1b)) then
      k = k+1                                  ! floating or grounded ice
      vis_int_sgxy(j,i) = vis_int_sgxy(j,i) + vis_int_g(j,i+1)
   end if

   if ((maske(j+1,i)==0_i1b).or.(maske(j+1,i)==3_i1b)) then
      k = k+1                                  ! floating or grounded ice
      vis_int_sgxy(j,i) = vis_int_sgxy(j,i) + vis_int_g(j+1,i)
   end if

   if ((maske(j+1,i+1)==0_i1b).or.(maske(j+1,i+1)==3_i1b)) then
      k = k+1                                      ! floating or grounded ice
      vis_int_sgxy(j,i) = vis_int_sgxy(j,i) + vis_int_g(j+1,i+1)
   end if

   if (k>0) vis_int_sgxy(j,i) = vis_int_sgxy(j,i)/real(k,dp)

end do
end do

!-------- Basal drag parameter (for shelfy stream) --------

beta_drag = 0.0_dp   ! initialisation

do i=1, IMAX-1
do j=1, JMAX-1

   if (flag_shelfy_stream(j,i)) then

      beta_drag(j,i) = c_drag(j,i) &
                     / sqrt( (   (0.5_dp*(vx_m_ssa(j,i)+vx_m_ssa(j,i-1)))**2  &
                               + (0.5_dp*(vy_m_ssa(j,i)+vy_m_ssa(j-1,i)))**2 ) &
                             + eps_dp**2 ) &
                                     **(1.0_dp-p_weert_inv(j,i))

   end if

end do
end do

!-------- Assembly of the system of linear equations
!                         (matrix storage: compressed sparse row CSR) --------

#if !defined(ALLOW_OPENAD) /* Normal */
allocate(lgs_a_value(n_sprs), lgs_a_index(n_sprs), lgs_a_ptr(nmax+1))
allocate(lgs_b_value(nmax), lgs_x_value(nmax))
#endif /* Normal */

lgs_a_value = 0.0_dp
#if !defined(ALLOW_OPENAD) /* Normal */
lgs_a_index = 0
#else /* OpenAD */
lgs_a_index = 0.0_dp
#endif /* Normal vs. OpenAD */
lgs_a_ptr   = 0

lgs_b_value = 0.0_dp
lgs_x_value = 0.0_dp

lgs_a_ptr(1) = 1

k = 0

do n=1, nmax-1, 2

   i = n2i((n+1)/2)
   j = n2j((n+1)/2)

!  ------ Equations for vx_m_ssa (at (i+1/2,j))

   nr = n   ! row counter

   if ( (i /= IMAX).and.(j /= 0).and.(j /= JMAX) ) then
      ! inner point on the staggered grid in x-direction

      H_mid  = 0.5_dp*((H_c(j,i)+H_t(j,i))+(H_c(j,i+1)+H_t(j,i+1)))
      zl_mid = 0.5_dp*(zl(j,i)+zl(j,i+1))

      if ( &
           ( (maske(j,i)==3_i1b).and.(maske(j,i+1)==3_i1b) ) &
           .or. &
           ( flag_grounding_line_1(j,i).and.flag_grounding_line_2(j,i+1) &
             .and.(H_mid < (z_sl-zl_mid)*rhosw_rho_ratio) ) &
           .or. &
           ( flag_grounding_line_2(j,i).and.flag_grounding_line_1(j,i+1) &
             .and.(H_mid < (z_sl-zl_mid)*rhosw_rho_ratio) ) &
         ) then
           ! both neighbours are floating ice
           !   or
           ! one neighbour is floating ice and the other is grounded ice
           ! (grounding line)
           ! and floating conditions are satisfied;
           ! discretization of the x-component of the PDE

         flag_calc_vxy_ssa_x(j,i) = .true.

         flag_shelfy_stream_x(j,i) = .false.
                                   ! make sure not to treat as shelfy stream

         nc = 2*ij2n(j,i-1)-1
                  ! smallest nc (column counter), for vx_m(j,i-1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 4.0_dp*inv_dxi2*vis_int_g(j,i)
         lgs_a_index(k) = nc

         nc = 2*ij2n(j-1,i)-1
                  ! next nc (column counter), for vx_m(j-1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_deta2*vis_int_sgxy(j-1,i)
         lgs_a_index(k) = nc

         nc = 2*ij2n(j-1,i)
                  ! next nc (column counter), for vy_m(j-1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j-1,i))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i)-1
                  ! next nc (column counter), for vx_m(j,i)
         if (nc /= nr) then   ! (diagonal element)
            errormsg = ' >>> calc_vxy_ssa_matrix: ' &
                          //'Check for diagonal element failed!'
            call error(errormsg)
         end if
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -4.0_dp*inv_dxi2 &
                                 *(vis_int_g(j,i+1)+vis_int_g(j,i)) &
                          -inv_deta2 &
                                 *(vis_int_sgxy(j,i)+vis_int_sgxy(j-1,i))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i)
                  ! next nc (column counter), for vy_m(j,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j,i))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j+1,i)-1
                  ! next nc (column counter), for vx_m(j+1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_deta2*vis_int_sgxy(j,i)
         lgs_a_index(k) = nc

         nc = 2*ij2n(j-1,i+1)
                  ! next nc (column counter), for vy_m(j-1,i+1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j,i+1)+vis_int_sgxy(j-1,i))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i+1)-1
                  ! next nc (column counter), for vx_m(j,i+1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 4.0_dp*inv_dxi2*vis_int_g(j,i+1)
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i+1)
                  ! largest nc (column counter), for vy_m(j,i+1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j,i+1)+vis_int_sgxy(j,i))
         lgs_a_index(k) = nc

         lgs_b_value(nr) = factor_rhs_1*H_mid*dzs_dx_aux(j,i)

         lgs_x_value(nr) = vx_m_ssa(j,i)
 
      else if (flag_shelfy_stream_x(j,i)) then
           ! shelfy stream (as determined by routine calc_vxy_sia)

         flag_calc_vxy_ssa_x(j,i) = .true.

#if (!defined(BC_SSA_LTGF) || BC_SSA_LTGF==1)

         if ( &
              ( ( flag_grounded_front_b_1(j,i) &
                       .and.flag_grounded_front_b_2(j,i+1) ) &
                .or. &
                ( flag_grounded_front_b_2(j,i) &
                       .and.flag_grounded_front_b_1(j,i+1) ) ) &
              .and. &
              ( zl_mid < z_sl ) &
            ) then
            ! one neighbour is grounded ice and the other is ocean
            ! (ocean-terminating grounded front)

#elif (BC_SSA_LTGF==2)

         if ( &
              ( flag_grounded_front_b_1(j,i) &
                     .and.flag_grounded_front_b_2(j,i+1) ) &
              .or. &
              ( flag_grounded_front_b_2(j,i) &
                     .and.flag_grounded_front_b_1(j,i+1) ) &
            ) then
            ! one neighbour is grounded ice and the other is ocean
            ! (ocean-terminating grounded front)

#else
         errormsg = ' >>> calc_vxy_ssa_matrix: ' &
                       //'BC_SSA_LTGF must be either 1 or 2!'
         call error(errormsg)
#endif

            if (flag_grounded_front_b_1(j,i)) then
               i1 = i     ! grounded ice marker
            else   ! flag_grounded_front_b_1(j,i+1)==.true.
               i1 = i+1   ! grounded ice marker
            end if

            if (.not.( flag_grounded_front_b_2(j,i1-1) &
                       .and. &
                       flag_grounded_front_b_2(j,i1+1) ) ) then
               ! discretization of the x-component of the BC

               nc = 2*ij2n(j,i1-1)-1
                        ! smallest nc (column counter), for vx_m(j,i1-1)
               k  = k+1
               lgs_a_value(k) = -4.0_dp*inv_dxi*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j-1,i1)
                        ! next nc (column counter), for vy_m(j-1,i1)
               k  = k+1
               lgs_a_value(k) = -2.0_dp*inv_deta*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j,i1)-1
                        ! next nc (column counter), for vx_m(j,i1)
               k  = k+1
               lgs_a_value(k) = 4.0_dp*inv_dxi*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j,i1)
                        ! largest nc (column counter), for vy_m(j,i1)
               k  = k+1
               lgs_a_value(k) = 2.0_dp*inv_deta*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               lgs_b_value(nr) = factor_rhs_3a &
                                    *(H_c(j,i1)+H_t(j,i1))**2 &
                               - factor_rhs_3b &
                                    *(max((z_sl-zb(j,i1)), 0.0_dp))**2

               lgs_x_value(nr) = vx_m_ssa(j,i)

            else   !      (flag_grounded_front_b_2(j,i1-1)==.true.)
                   ! .and.(flag_grounded_front_b_2(j,i1+1)==.true.);
                   ! velocity assumed to be zero

               k  = k+1
               lgs_a_value(k) = 1.0_dp   ! diagonal element only
               lgs_a_index(k) = nr

               lgs_b_value(nr) = 0.0_dp

               lgs_x_value(nr) = 0.0_dp

            end if

#if (BC_SSA_LTGF==2)

         else if ( &
              ( flag_grounded_front_a_1(j,i) &
                     .and.flag_grounded_front_a_2(j,i+1) ) &
              .or. &
              ( flag_grounded_front_a_2(j,i) &
                     .and.flag_grounded_front_a_1(j,i+1) ) &
            ) then
            ! one neighbour is grounded ice and the other is ice-free land
            ! (land-terminating grounded front)

            if (flag_grounded_front_a_1(j,i)) then
               i1 = i     ! grounded ice marker
            else   ! flag_grounded_front_a_1(j,i+1)==.true.
               i1 = i+1   ! grounded ice marker
            end if

            if (.not.( flag_grounded_front_a_2(j,i1-1) &
                       .and. &
                       flag_grounded_front_a_2(j,i1+1) ) ) then
               ! discretization of the x-component of the BC

               nc = 2*ij2n(j,i1-1)-1
                        ! smallest nc (column counter), for vx_m(j,i1-1)
               k  = k+1
               lgs_a_value(k) = -4.0_dp*inv_dxi*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j-1,i1)
                        ! next nc (column counter), for vy_m(j-1,i1)
               k  = k+1
               lgs_a_value(k) = -2.0_dp*inv_deta*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j,i1)-1
                        ! next nc (column counter), for vx_m(j,i1)
               k  = k+1
               lgs_a_value(k) = 4.0_dp*inv_dxi*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j,i1)
                        ! largest nc (column counter), for vy_m(j,i1)
               k  = k+1
               lgs_a_value(k) = 2.0_dp*inv_deta*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               lgs_b_value(nr) = factor_rhs_3a &
                                    *(H_c(j,i1)+H_t(j,i1))**2

               lgs_x_value(nr) = vx_m_ssa(j,i)

            else   !      (flag_grounded_front_a_2(j,i1-1)==.true.)
                   ! .and.(flag_grounded_front_a_2(j,i1+1)==.true.);
                   ! velocity assumed to be zero

               k  = k+1
               lgs_a_value(k) = 1.0_dp   ! diagonal element only
               lgs_a_index(k) = nr

               lgs_b_value(nr) = 0.0_dp

               lgs_x_value(nr) = 0.0_dp

            end if

#endif

         else
            ! inner shelfy stream
            !   or
            ! one neighbour is floating ice and the other is grounded ice
            ! (grounding line)
            ! and floating conditions are not satisfied
#if (!defined(BC_SSA_LTGF) || BC_SSA_LTGF==1)
            !   or
            ! land-terminating grounded front
#endif

            nc = 2*ij2n(j,i-1)-1
                     ! smallest nc (column counter), for vx_m(j,i-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 4.0_dp*inv_dxi2*vis_int_g(j,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j-1,i)-1
                     ! next nc (column counter), for vx_m(j-1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_deta2*vis_int_sgxy(j-1,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j-1,i)
                     ! next nc (column counter), for vy_m(j-1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j-1,i))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i)-1
                     ! next nc (column counter), for vx_m(j,i)
            if (nc /= nr) then   ! (diagonal element)
               errormsg = ' >>> calc_vxy_ssa_matrix: ' &
                             //'Check for diagonal element failed!'
               call error(errormsg)
            end if
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -4.0_dp*inv_dxi2 &
                                    *(vis_int_g(j,i+1)+vis_int_g(j,i)) &
                             -inv_deta2 &
                                    *(vis_int_sgxy(j,i)+vis_int_sgxy(j-1,i)) &
                             -0.5_dp*(beta_drag(j,i+1)+beta_drag(j,i))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i)
                     ! next nc (column counter), for vy_m(j,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                    *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j,i))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j+1,i)-1
                     ! next nc (column counter), for vx_m(j+1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_deta2*vis_int_sgxy(j,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j-1,i+1)
                     ! next nc (column counter), for vy_m(j-1,i+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                  *(2.0_dp*vis_int_g(j,i+1)+vis_int_sgxy(j-1,i))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i+1)-1
                     ! next nc (column counter), for vx_m(j,i+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 4.0_dp*inv_dxi2*vis_int_g(j,i+1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i+1)
                     ! largest nc (column counter), for vy_m(j,i+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_dp*vis_int_g(j,i+1)+vis_int_sgxy(j,i))
            lgs_a_index(k) = nc

            lgs_b_value(nr) = factor_rhs_1*H_mid*dzs_dx_aux(j,i)

            lgs_x_value(nr) = vx_m_ssa(j,i)

         end if

      else if ( &
                ( flag_grounding_line_1(j,i).and.flag_grounding_line_2(j,i+1) &
                  .and.(H_mid >= (z_sl-zl_mid)*rhosw_rho_ratio) ) &
                .or. &
                ( flag_grounding_line_2(j,i).and.flag_grounding_line_1(j,i+1) &
                  .and.(H_mid >= (z_sl-zl_mid)*rhosw_rho_ratio) ) &
              ) then 
              ! one neighbour is floating ice and the other is grounded ice
              ! (grounding line)
              ! and floating conditions are not satisfied;
              ! velocity taken from the SIA solution for grounded ice

         flag_calc_vxy_ssa_x(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = vx_m_sia(j,i)

         lgs_x_value(nr) = vx_m_sia(j,i)

      else if ( &
                ( (maske(j,i)==3_i1b).and.(maske(j,i+1)==1_i1b) ) &
                .or. &
                ( (maske(j,i)==1_i1b).and.(maske(j,i+1)==3_i1b) ) &
              ) then
              ! one neighbour is floating ice and the other is ice-free land;
              ! velocity assumed to be zero

         flag_calc_vxy_ssa_x(j,i) = .true.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = 0.0_dp

         lgs_x_value(nr) = 0.0_dp

      else if ( &
                ( flag_calving_front_1(j,i).and.flag_calving_front_2(j,i+1) ) &
                .or. &
                ( flag_calving_front_2(j,i).and.flag_calving_front_1(j,i+1) ) &
              ) then
              ! one neighbour is floating ice and the other is ocean
              ! (calving front)

         flag_calc_vxy_ssa_x(j,i) = .true.

         if (flag_calving_front_1(j,i)) then
            i1 = i     ! floating ice marker
         else   ! flag_calving_front_1(j,i+1)==.true.
            i1 = i+1   ! floating ice marker
         end if

         if (.not.( flag_calving_front_2(j,i1-1) &
                    .and. &
                    flag_calving_front_2(j,i1+1) ) ) then
            ! discretization of the x-component of the BC

            nc = 2*ij2n(j,i1-1)-1
                     ! smallest nc (column counter), for vx_m(j,i1-1)
            k  = k+1
            lgs_a_value(k) = -4.0_dp*inv_dxi*vis_int_g(j,i1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j-1,i1)
                     ! next nc (column counter), for vy_m(j-1,i1)
            k  = k+1
            lgs_a_value(k) = -2.0_dp*inv_deta*vis_int_g(j,i1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i1)-1
                     ! next nc (column counter), for vx_m(j,i1)
            k  = k+1
            lgs_a_value(k) = 4.0_dp*inv_dxi*vis_int_g(j,i1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i1)
                     ! largest nc (column counter), for vy_m(j,i1)
            k  = k+1
            lgs_a_value(k) = 2.0_dp*inv_deta*vis_int_g(j,i1)
            lgs_a_index(k) = nc

            lgs_b_value(nr) = factor_rhs_2 &
                                 *(H_c(j,i1)+H_t(j,i1))**2

            lgs_x_value(nr) = vx_m_ssa(j,i)

         else   !      (flag_calving_front_2(j,i1-1)==.true.)
                ! .and.(flag_calving_front_2(j,i1+1)==.true.);
                ! velocity assumed to be zero

            k  = k+1
            lgs_a_value(k) = 1.0_dp   ! diagonal element only
            lgs_a_index(k) = nr

            lgs_b_value(nr) = 0.0_dp

            lgs_x_value(nr) = 0.0_dp

         end if

      else if ( (maske(j,i)==0_i1b).or.(maske(j,i+1)==0_i1b) ) then
           ! neither neighbour is floating ice, but at least one neighbour is
           ! grounded ice; velocity taken from the SIA solution for grounded ice

         flag_calc_vxy_ssa_x(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = vx_m_sia(j,i)

         lgs_x_value(nr) = vx_m_sia(j,i)

      else   ! neither neighbour is floating or grounded ice,
             ! velocity assumed to be zero

         flag_calc_vxy_ssa_x(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = 0.0_dp

         lgs_x_value(nr) = 0.0_dp

      end if

   else   ! boundary condition, velocity assumed to be zero

      flag_calc_vxy_ssa_x(j,i) = .false.

      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
      lgs_a_value(k) = 1.0_dp   ! diagonal element only
      lgs_a_index(k) = nr

      lgs_b_value(nr) = 0.0_dp

      lgs_x_value(nr) = 0.0_dp

   end if

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

!  ------ Equations for vy_m_ssa (at (i,j+1/2))

   nr = n+1   ! row counter

   if ( (j /= JMAX).and.(i /= 0).and.(i /= IMAX) ) then
      ! inner point on the staggered grid in y-direction

      H_mid  = 0.5_dp*((H_c(j,i)+H_t(j,i))+(H_c(j+1,i)+H_t(j+1,i)))
      zl_mid = 0.5_dp*(zl(j,i)+zl(j+1,i))
   
      if ( &
           ( (maske(j,i)==3_i1b).and.(maske(j+1,i)==3_i1b) ) &
           .or. &
           ( flag_grounding_line_1(j,i).and.flag_grounding_line_2(j+1,i) &
             .and.(H_mid < (z_sl-zl_mid)*rhosw_rho_ratio) ) &
           .or. &
           ( flag_grounding_line_2(j,i).and.flag_grounding_line_1(j+1,i) &
             .and.(H_mid < (z_sl-zl_mid)*rhosw_rho_ratio) ) &
         ) then
           ! both neighbours are floating ice
           !   or
           ! one neighbour is floating ice and the other is grounded ice
           ! (grounding line)
           ! and floating conditions are satisfied;
           ! discretization of the y-component of the PDE

         flag_calc_vxy_ssa_y(j,i) = .true.

         flag_shelfy_stream_y(j,i) = .false.
                                   ! make sure not to treat as shelfy stream

         nc = 2*ij2n(j,i-1)-1
                  ! smallest nc (column counter), for vx_m(j,i-1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j,i-1))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i-1)
                  ! next nc (column counter), for vy_m(j,i-1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi2*vis_int_sgxy(j,i-1)
         lgs_a_index(k) = nc

         nc = 2*ij2n(j+1,i-1)-1
                  ! next nc (column counter), for vx_m(j+1,i-1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j+1,i)+vis_int_sgxy(j,i-1))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j-1,i)
                  ! next nc (column counter), for vy_m(j-1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 4.0_dp*inv_deta2*vis_int_g(j,i)
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i)-1
                  ! next nc (column counter), for vx_m(j,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j,i))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i)
                  ! next nc (column counter), for vy_m(j,i)
         if (nc /= nr) then   ! (diagonal element)
            errormsg = ' >>> calc_vxy_ssa_matrix: ' &
                          //'Check for diagonal element failed!'
            call error(errormsg)
         end if
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -4.0_dp*inv_deta2 &
                                 *(vis_int_g(j+1,i)+vis_int_g(j,i)) &
                          -inv_dxi2 &
                                 *(vis_int_sgxy(j,i)+vis_int_sgxy(j,i-1))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j+1,i)-1
                  ! next nc (column counter), for vx_m(j+1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j+1,i)+vis_int_sgxy(j,i))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j+1,i)
                  ! next nc (column counter), for vy_m(j+1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 4.0_dp*inv_deta2*vis_int_g(j+1,i)
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i+1)
                  ! largest nc (column counter), for vy_m(j,i+1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi2*vis_int_sgxy(j,i)
         lgs_a_index(k) = nc

         lgs_b_value(nr) = factor_rhs_1*H_mid*dzs_dy_aux(j,i)

         lgs_x_value(nr) = vy_m_ssa(j,i)
 
      else if (flag_shelfy_stream_y(j,i)) then
           ! shelfy stream (as determined by routine calc_vxy_sia)

         flag_calc_vxy_ssa_y(j,i) = .true.

#if (!defined(BC_SSA_LTGF) || BC_SSA_LTGF==1)

         if ( &
              ( ( flag_grounded_front_b_1(j,i) &
                       .and.flag_grounded_front_b_2(j+1,i) ) &
                .or. &
                ( flag_grounded_front_b_2(j,i) &
                       .and.flag_grounded_front_b_1(j+1,i) ) ) &
              .and. &
              ( zl_mid < z_sl ) &
            ) then
            ! one neighbour is grounded ice and the other is ocean
            ! (ocean-terminating grounded front)

#elif (BC_SSA_LTGF==2)

         if ( &
              ( flag_grounded_front_b_1(j,i) &
                     .and.flag_grounded_front_b_2(j+1,i) ) &
              .or. &
              ( flag_grounded_front_b_2(j,i) &
                     .and.flag_grounded_front_b_1(j+1,i) ) &
            ) then
            ! one neighbour is grounded ice and the other is ocean
            ! (ocean-terminating grounded front)

#else
         errormsg = ' >>> calc_vxy_ssa_matrix: ' &
                       //'BC_SSA_LTGF must be either 1 or 2!'
         call error(errormsg)
#endif

            if (flag_grounded_front_b_1(j,i)) then
               j1 = j     ! grounded ice marker
            else   ! flag_grounded_front_b_1(j+1,i)==.true.
               j1 = j+1   ! grounded ice marker
            end if

            if (.not.( flag_grounded_front_b_2(j1-1,i) &
                       .and. &
                       flag_grounded_front_b_2(j1+1,i) ) ) then
               ! discretization of the y-component of the BC

               nc = 2*ij2n(j1,i-1)-1
                        ! smallest nc (column counter), for vx_m(j1,i-1)
               k  = k+1
               lgs_a_value(k) = -2.0_dp*inv_dxi*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j1-1,i)
                        ! next nc (column counter), for vy_m(j1-1,i)
               k  = k+1
               lgs_a_value(k) = -4.0_dp*inv_deta*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j1,i)-1
                        ! next nc (column counter), for vx_m(j1,i)
               k  = k+1
               lgs_a_value(k) = 2.0_dp*inv_dxi*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j1,i)
                        ! largest nc (column counter), for vy_m(j1,i)
               k  = k+1
               lgs_a_value(k) = 4.0_dp*inv_deta*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               lgs_b_value(nr) = factor_rhs_3a &
                                    *(H_c(j1,i)+H_t(j1,i))**2 &
                               - factor_rhs_3b &
                                    *(max((z_sl-zb(j1,i)), 0.0_dp))**2

               lgs_x_value(nr) = vy_m_ssa(j,i)

            else   !      (flag_grounded_front_b_2(j1-1,i)==.true.)
                   ! .and.(flag_grounded_front_b_2(j1+1,i)==.true.);
                   ! velocity assumed to be zero

               k  = k+1
               lgs_a_value(k) = 1.0_dp   ! diagonal element only
               lgs_a_index(k) = nr

               lgs_b_value(nr) = 0.0_dp

               lgs_x_value(nr) = 0.0_dp

            end if

#if (BC_SSA_LTGF==2)

         else if ( &
              ( flag_grounded_front_a_1(j,i) &
                     .and.flag_grounded_front_a_2(j+1,i) ) &
              .or. &
              ( flag_grounded_front_a_2(j,i) &
                     .and.flag_grounded_front_a_1(j+1,i) ) &
            ) then
            ! one neighbour is grounded ice and the other is ice-free land
            ! (land-terminating grounded front)

            if (flag_grounded_front_a_1(j,i)) then
               j1 = j     ! grounded ice marker
            else   ! flag_grounded_front_a_1(j+1,i)==.true.
               j1 = j+1   ! grounded ice marker
            end if

            if (.not.( flag_grounded_front_a_2(j1-1,i) &
                       .and. &
                       flag_grounded_front_a_2(j1+1,i) ) ) then
               ! discretization of the y-component of the BC

               nc = 2*ij2n(j1,i-1)-1
                        ! smallest nc (column counter), for vx_m(j1,i-1)
               k  = k+1
               lgs_a_value(k) = -2.0_dp*inv_dxi*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j1-1,i)
                        ! next nc (column counter), for vy_m(j1-1,i)
               k  = k+1
               lgs_a_value(k) = -4.0_dp*inv_deta*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j1,i)-1
                        ! next nc (column counter), for vx_m(j1,i)
               k  = k+1
               lgs_a_value(k) = 2.0_dp*inv_dxi*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j1,i)
                        ! largest nc (column counter), for vy_m(j1,i)
               k  = k+1
               lgs_a_value(k) = 4.0_dp*inv_deta*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               lgs_b_value(nr) = factor_rhs_3a &
                                    *(H_c(j1,i)+H_t(j1,i))**2

               lgs_x_value(nr) = vy_m_ssa(j,i)

            else   !      (flag_grounded_front_a_2(j1-1,i)==.true.)
                   ! .and.(flag_grounded_front_a_2(j1+1,i)==.true.);
                   ! velocity assumed to be zero

               k  = k+1
               lgs_a_value(k) = 1.0_dp   ! diagonal element only
               lgs_a_index(k) = nr

               lgs_b_value(nr) = 0.0_dp

               lgs_x_value(nr) = 0.0_dp

            end if

#endif

         else
            ! inner shelfy stream
            !   or
            ! one neighbour is floating ice and the other is grounded ice
            ! (grounding line)
            ! and floating conditions are not satisfied
#if (!defined(BC_SSA_LTGF) || BC_SSA_LTGF==1)
            !   or
            ! land-terminating grounded front
#endif

            nc = 2*ij2n(j,i-1)-1
                     ! smallest nc (column counter), for vx_m(j,i-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j,i-1))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i-1)
                     ! next nc (column counter), for vy_m(j,i-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi2*vis_int_sgxy(j,i-1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j+1,i-1)-1
                     ! next nc (column counter), for vx_m(j+1,i-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                  *(2.0_dp*vis_int_g(j+1,i)+vis_int_sgxy(j,i-1))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j-1,i)
                     ! next nc (column counter), for vy_m(j-1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 4.0_dp*inv_deta2*vis_int_g(j,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i)-1
                     ! next nc (column counter), for vx_m(j,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                    *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j,i))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i)
                     ! next nc (column counter), for vy_m(j,i)
            if (nc /= nr) then   ! (diagonal element)
               errormsg = ' >>> calc_vxy_ssa_matrix: ' &
                             //'Check for diagonal element failed!'
               call error(errormsg)
            end if
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -4.0_dp*inv_deta2 &
                                    *(vis_int_g(j+1,i)+vis_int_g(j,i)) &
                             -inv_dxi2 &
                                    *(vis_int_sgxy(j,i)+vis_int_sgxy(j,i-1)) &
                             -0.5_dp*(beta_drag(j+1,i)+beta_drag(j,i))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j+1,i)-1
                     ! next nc (column counter), for vx_m(j+1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_dp*vis_int_g(j+1,i)+vis_int_sgxy(j,i))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j+1,i)
                     ! next nc (column counter), for vy_m(j+1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 4.0_dp*inv_deta2*vis_int_g(j+1,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i+1)
                     ! largest nc (column counter), for vy_m(j,i+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi2*vis_int_sgxy(j,i)
            lgs_a_index(k) = nc

            lgs_b_value(nr) = factor_rhs_1*H_mid*dzs_dy_aux(j,i)

            lgs_x_value(nr) = vy_m_ssa(j,i)

         end if

      else if ( &
                ( flag_grounding_line_1(j,i).and.flag_grounding_line_2(j+1,i) &
                  .and.(H_mid >= (z_sl-zl_mid)*rhosw_rho_ratio) ) &
                .or. &
                ( flag_grounding_line_2(j,i).and.flag_grounding_line_1(j+1,i) &
                  .and.(H_mid >= (z_sl-zl_mid)*rhosw_rho_ratio) ) &
              ) then
              ! one neighbour is floating ice and the other is grounded ice
              ! (grounding line)
              ! and floating conditions are not satisfied;
              ! velocity taken from the SIA solution for grounded ice

         flag_calc_vxy_ssa_y(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = vy_m_sia(j,i)

         lgs_x_value(nr) = vy_m_sia(j,i)

      else if ( &
                ( (maske(j,i)==3_i1b).and.(maske(j+1,i)==1_i1b) ) &
                .or. &
                ( (maske(j,i)==1_i1b).and.(maske(j+1,i)==3_i1b) ) &
              ) then
           ! one neighbour is floating ice and the other is ice-free land;
           ! velocity assumed to be zero

         flag_calc_vxy_ssa_y(j,i) = .true.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = 0.0_dp

         lgs_x_value(nr) = 0.0_dp

      else if ( &
                ( flag_calving_front_1(j,i).and.flag_calving_front_2(j+1,i) ) &
                .or. &
                ( flag_calving_front_2(j,i).and.flag_calving_front_1(j+1,i) ) &
              ) then
              ! one neighbour is floating ice and the other is ocean
              ! (calving front)

         flag_calc_vxy_ssa_y(j,i) = .true.

         if (flag_calving_front_1(j,i)) then
            j1 = j     ! floating ice marker
         else   ! flag_calving_front_1(j+1,i)==.true.
            j1 = j+1   ! floating ice marker
         end if

         if (.not.( flag_calving_front_2(j1-1,i) &
                    .and. &
                    flag_calving_front_2(j1+1,i) ) ) then
            ! discretization of the y-component of the BC

            nc = 2*ij2n(j1,i-1)-1
                     ! smallest nc (column counter), for vx_m(j1,i-1)
            k  = k+1
            lgs_a_value(k) = -2.0_dp*inv_dxi*vis_int_g(j1,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j1-1,i)
                     ! next nc (column counter), for vy_m(j1-1,i)
            k  = k+1
            lgs_a_value(k) = -4.0_dp*inv_deta*vis_int_g(j1,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j1,i)-1
                     ! next nc (column counter), for vx_m(j1,i)
            k  = k+1
            lgs_a_value(k) = 2.0_dp*inv_dxi*vis_int_g(j1,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j1,i)
                     ! largest nc (column counter), for vy_m(j1,i)
            k  = k+1
            lgs_a_value(k) = 4.0_dp*inv_deta*vis_int_g(j1,i)
            lgs_a_index(k) = nc

            lgs_b_value(nr) = factor_rhs_2 &
                                 *(H_c(j1,i)+H_t(j1,i))**2

            lgs_x_value(nr) = vy_m_ssa(j,i)

         else   !      (flag_calving_front_2(j1-1,i)==.true.)
                ! .and.(flag_calving_front_2(j1+1,i)==.true.);
                ! velocity assumed to be zero

            k  = k+1
            lgs_a_value(k) = 1.0_dp   ! diagonal element only
            lgs_a_index(k) = nr

            lgs_b_value(nr) = 0.0_dp

            lgs_x_value(nr) = 0.0_dp

         end if

      else if ( (maske(j,i)==0_i1b).or.(maske(j+1,i)==0_i1b) ) then
           ! neither neighbour is floating ice, but at least one neighbour is
           ! grounded ice; velocity taken from the SIA solution for grounded ice

         flag_calc_vxy_ssa_y(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = vy_m_sia(j,i)

         lgs_x_value(nr) = vy_m_sia(j,i)

      else   ! neither neighbour is floating or grounded ice,
             ! velocity assumed to be zero

         flag_calc_vxy_ssa_y(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = 0.0_dp

         lgs_x_value(nr) = 0.0_dp

      end if

   else   ! boundary condition, velocity assumed to be zero

      flag_calc_vxy_ssa_y(j,i) = .false.

      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
      lgs_a_value(k) = 1.0_dp   ! diagonal element only
      lgs_a_index(k) = nr

      lgs_b_value(nr) = 0.0_dp

      lgs_x_value(nr) = 0.0_dp

   end if

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

end do

!-------- Settings for Lis --------

#if !defined(ALLOW_OPENAD) /* Normal */

call lis_matrix_create(LIS_COMM_WORLD, lgs_a, ierr)
call lis_vector_create(LIS_COMM_WORLD, lgs_b, ierr)
call lis_vector_create(LIS_COMM_WORLD, lgs_x, ierr)

call lis_matrix_set_size(lgs_a, 0, nmax, ierr)
call lis_vector_set_size(lgs_b, 0, nmax, ierr)
call lis_vector_set_size(lgs_x, 0, nmax, ierr)

do nr=1, nmax

   do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
      call lis_matrix_set_value(LIS_INS_VALUE, nr, lgs_a_index(nc), &
                                               lgs_a_value(nc), lgs_a, ierr)
   end do

   call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_b_value(nr), lgs_b, ierr)
   call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_x_value(nr), lgs_x, ierr)

end do

call lis_matrix_set_type(lgs_a, LIS_MATRIX_CSR, ierr)
call lis_matrix_assemble(lgs_a, ierr)

!-------- Solution of the system of linear equations with Lis --------

call lis_solver_create(solver, ierr)

#if (defined(LIS_OPTS))
    ch_solver_set_option = trim(LIS_OPTS)
#else
    ch_solver_set_option = '-i bicgsafe -p jacobi '// &
                           '-maxiter 100 -tol 1.0e-4 -initx_zeros false'
#endif

call lis_solver_set_option(trim(ch_solver_set_option), solver, ierr)
call CHKERR(ierr)

call lis_solve(lgs_a, lgs_b, lgs_x, solver, ierr)
call CHKERR(ierr)

call lis_solver_get_iter(solver, lin_iter, ierr)

write(6,'(a,i0,a)', advance='no') 'lin_iter = ', lin_iter, ', '

!!! call lis_solver_get_time(solver,solver_time,ierr)
!!! print *, 'calc_vxy_ssa_matrix: time (s) = ', solver_time

lgs_x_value = 0.0_dp
call lis_vector_gather(lgs_x, lgs_x_value, ierr)
call lis_matrix_destroy(lgs_a, ierr)
call lis_vector_destroy(lgs_b, ierr)
call lis_vector_destroy(lgs_x, ierr)
call lis_solver_destroy(solver, ierr)

#else /* OpenAD */

lgs_a_index_pass = int(lgs_a_index)
call sico_lis_solver(nmax, n_sprs, &
                           lgs_a_ptr, lgs_a_index_pass, &
                           lgs_a_value, lgs_b_value, lgs_x_value)

#endif /* Normal vs. OpenAD */

do n=1, nmax-1, 2

   i = n2i((n+1)/2)
   j = n2j((n+1)/2)

   nr = n
   vx_m_ssa(j,i) = lgs_x_value(nr)

   nr = n+1
   vy_m_ssa(j,i) = lgs_x_value(nr)

end do

#if !defined(ALLOW_OPENAD) /* Normal */
deallocate(lgs_a_value, lgs_a_index, lgs_a_ptr)
deallocate(lgs_b_value, lgs_x_value)
#endif /* Normal */

#else

errormsg = ' >>> calc_vxy_ssa_matrix: ' &
              //'Only to be called for MARGIN==3 or DYNAMICS==2!'
call error(errormsg)

#endif

end subroutine calc_vxy_ssa_matrix

!-------------------------------------------------------------------------------
!> Computation of the depth-integrated viscosity vis_int_g in the
!! shallow shelf approximation.
!<------------------------------------------------------------------------------
subroutine calc_vis_ssa(dxi, deta, dzeta_c, dzeta_t)

use ice_material_properties_m, only : viscosity

implicit none

real(dp), intent(in) :: dxi, deta, dzeta_c, dzeta_t

integer(i4b) :: i, j, kc, kt, m
integer(i4b) :: m_smooth, m_smooth_abs
real(dp) :: visc_min, visc_max, visc_init
real(dp) :: dxi_inv, deta_inv
real(dp) :: dvx_dxi, dvx_deta, dvy_dxi, dvy_deta
real(dp) :: diff_vis
real(dp) :: aqxy1(0:KCMAX)
real(dp) :: cvis0(0:KTMAX), cvis1(0:KCMAX)
real(dp), dimension(0:JMAX,0:IMAX) :: vis_ave_g_smooth

#if (MARGIN==3 || DYNAMICS==2)

!-------- Parameters, term abbreviations --------

#if (defined(VISC_MIN) && defined(VISC_MAX))
  visc_min = VISC_MIN
  visc_max = VISC_MAX
#else
  visc_min = 1.0e+10_dp   ! Pa s
  visc_max = 1.0e+25_dp   ! Pa s
#endif

#if (defined(VISC_INIT_SSA))
  visc_init = VISC_INIT_SSA
#else
  visc_init = 1.0e+15_dp   ! Pa s
#endif

dxi_inv  = 1.0_dp/dxi
deta_inv = 1.0_dp/deta

do kc=0, KCMAX
   if (flag_aa_nonzero) then
      aqxy1(kc) = aa/(ea-1.0_dp)*eaz_c(kc)*dzeta_c
   else
      aqxy1(kc) = dzeta_c
   end if
end do

!-------- Computation of the depth-integrated viscosity --------

do i=0, IMAX
do j=0, JMAX

   if ((maske(j,i)==0_i1b).and.(.not.flag_shelfy_stream(j,i))) then
                                                   ! grounded ice, but
                                                   ! not shelfy stream
      de_ssa(j,i) = 0.0_dp   ! dummy value

      vis_ave_g(j,i) = 1.0_dp/flui_ave_sia(j,i)
      vis_int_g(j,i) = (H_c(j,i)+H_t(j,i)) * vis_ave_g(j,i)

   else if ((maske(j,i)==1_i1b).or.(maske(j,i)==2_i1b)) then
                                                   ! ice-free land or ocean
      de_ssa(j,i) = 0.0_dp   ! dummy value

      vis_ave_g(j,i) = visc_init   ! dummy value
      vis_int_g(j,i) = 0.0_dp      ! dummy value

   else   ! (maske(j,i)==3_i1b).or.(flag_shelfy_stream(j,i)),
          ! floating ice or shelfy stream; 
          ! must not be at the margin of the computational domain

!  ------ Effective strain rate

      dvx_dxi  = (vx_m_ssa(j,i)-vx_m_ssa(j,i-1))*dxi_inv 
      dvy_deta = (vy_m_ssa(j,i)-vy_m_ssa(j-1,i))*deta_inv

      dvx_deta = 0.25_dp*deta_inv &
                    *( (vx_m_ssa(j+1,i)+vx_m_ssa(j+1,i-1)) &
                      -(vx_m_ssa(j-1,i)+vx_m_ssa(j-1,i-1)) )
      dvy_dxi  = 0.25_dp*dxi_inv &
                    *( (vy_m_ssa(j,i+1)+vy_m_ssa(j-1,i+1)) &
                      -(vy_m_ssa(j,i-1)+vy_m_ssa(j-1,i-1)) )

#if !defined(ALLOW_OPENAD) /* Normal */

      de_ssa(j,i) = sqrt( dvx_dxi*dvx_dxi &
                        + dvy_deta*dvy_deta &
                        + dvx_dxi*dvy_deta &
                        + 0.25_dp*(dvx_deta+dvy_dxi)*(dvx_deta+dvy_dxi) )

#else /* OpenAD: guarding against non-differentiable sqrt(0) */

      if ( ( dvx_dxi*dvx_dxi &
           + dvy_deta*dvy_deta &
           + dvx_dxi*dvy_deta &
           + 0.25_dp*(dvx_deta+dvy_dxi)*(dvx_deta+dvy_dxi) ) > 0 ) then
         de_ssa(j,i) = sqrt( dvx_dxi*dvx_dxi &
                           + dvy_deta*dvy_deta &
                           + dvx_dxi*dvy_deta &
                           + 0.25_dp*(dvx_deta+dvy_dxi)*(dvx_deta+dvy_dxi) )
      else
         de_ssa(j,i) = 0.0_dp 
      end if

#endif /* Normal vs. OpenAD */

!  ------ Term abbreviations

#if (DYNAMICS==2)
      if (.not.flag_shelfy_stream(j,i)) then
#endif

         do kc=0, KCMAX
            cvis1(kc) = aqxy1(kc)*H_c(j,i) &
                           *viscosity(de_ssa(j,i), &
                              temp_c(kc,j,i), temp_c_m(kc,j,i), &
                              0.0_dp, enh_c(kc,j,i), 0_i1b)
         end do
         ! Ice shelves (floating ice) are assumed to consist of cold ice only

#if (DYNAMICS==2)
      else   ! flag_shelfy_stream(j,i) == .true.

#if (CALCMOD==-1 || CALCMOD==0)

         do kc=0, KCMAX
            cvis1(kc) = aqxy1(kc)*H_c(j,i) &
                           *viscosity(de_ssa(j,i), &
                              temp_c(kc,j,i), temp_c_m(kc,j,i), &
                              0.0_dp, enh_c(kc,j,i), 0_i1b)
         end do

#elif (CALCMOD==1)

         do kt=0, KTMAX
            cvis0(kt) = dzeta_t*H_t(j,i) &
                           *viscosity(de_ssa(j,i), &
                              temp_t_m(kt,j,i), temp_t_m(kt,j,i), &
                              omega_t(kt,j,i), enh_t(kt,j,i), 1_i1b)
         end do

         do kc=0, KCMAX
            cvis1(kc) = aqxy1(kc)*H_c(j,i) &
                           *viscosity(de_ssa(j,i), &
                              temp_c(kc,j,i), temp_c_m(kc,j,i), &
                              0.0_dp, enh_c(kc,j,i), 0_i1b)
         end do

#elif (CALCMOD==2 || CALCMOD==3)

         do kc=0, KCMAX
            cvis1(kc) = aqxy1(kc)*H_c(j,i) &
                         *viscosity(de_ssa(j,i), &
                           temp_c(kc,j,i), temp_c_m(kc,j,i), &
                           omega_c(kc,j,i), enh_c(kc,j,i), 2_i1b)
         end do

#else
         errormsg = ' >>> calc_vis_ssa: CALCMOD must be -1, 0, 1, 2 or 3!'
         call error(errormsg)
#endif

      end if

#endif

!  ------ Depth-integrated viscosity

      vis_int_g(j,i) = 0.0_dp

#if (CALCMOD==1)
      do kt=0, KTMAX-1
         vis_int_g(j,i) = vis_int_g(j,i)+0.5_dp*(cvis0(kt+1)+cvis0(kt))
      end do
#endif

      do kc=0, KCMAX-1
         vis_int_g(j,i) = vis_int_g(j,i)+0.5_dp*(cvis1(kc+1)+cvis1(kc))
      end do

!  ------ Depth-averaged viscosity

      vis_ave_g(j,i) = vis_int_g(j,i)/max((H_c(j,i)+H_t(j,i)), eps_dp)

      vis_ave_g(j,i) = max(min(vis_ave_g(j,i), visc_max), visc_min)

   end if
  
end do
end do

!-------- Smoothing of the depth-averaged viscosity --------

#if (defined(N_VISC_SMOOTH))
  m_smooth = N_VISC_SMOOTH
#else
  m_smooth = 0
#endif

if (m_smooth /= 0) then

   m_smooth_abs = abs(m_smooth)

#if (defined(VISC_SMOOTH_DIFF))
   diff_vis = VISC_SMOOTH_DIFF
#else
   diff_vis = 0.0_dp
#endif

   if (m_smooth < 0) vis_ave_g = log(vis_ave_g)   ! logarithmic smoothing

   vis_ave_g_smooth = vis_ave_g

   do m=1, m_smooth_abs

      do i=1, IMAX-1
      do j=1, JMAX-1
         vis_ave_g_smooth(j,i) = (1.0_dp-4.0_dp*diff_vis)*vis_ave_g(j,i) &
                                    + diff_vis &
                                       *( (vis_ave_g(j,i+1)+vis_ave_g(j,i-1)) &
                                         +(vis_ave_g(j+1,i)+vis_ave_g(j-1,i)) )
      end do
      end do

      vis_ave_g = vis_ave_g_smooth

   end do

   if (m_smooth < 0) vis_ave_g = exp(vis_ave_g)   ! logarithmic smoothing

end if

!-------- Final depth-integrated viscosity --------

vis_int_g = vis_ave_g*(H_c+H_t)

#else

errormsg = ' >>> calc_vis_ssa: Only to be called for MARGIN==3 or DYNAMICS==2!'
call error(errormsg)

#endif

end subroutine calc_vis_ssa

!-------------------------------------------------------------------------------
!> Gradual limitation of computed horizontal velocities to the interval
!! [-vel_max, vel_max].
!<------------------------------------------------------------------------------
elemental subroutine velocity_limiter_gradual(velocity, vel_max, vel_max_inv)

implicit none

real(dp), intent(in)    :: vel_max, vel_max_inv
real(dp), intent(inout) :: velocity

real(dp) :: vel_abs, vel_sign, vel_scaled, vel_scaled_lim

vel_abs = abs(velocity)

if (vel_abs >= 1.1_dp*vel_max) then

   vel_sign = sign(1.0_dp, velocity)
   velocity = vel_sign * vel_max

else if (vel_abs > 0.9_dp*vel_max) then

   ! gradual limitation between 0.9*vel_max and 1.1*vel_max

   vel_sign = sign(1.0_dp, velocity)

   vel_scaled     = (vel_abs-0.9_dp*vel_max)*(10.0_dp*vel_max_inv)
                       ! between 0 and 2
   vel_scaled_lim = vel_scaled &
                       *(1.0_dp-0.25_dp*vel_scaled*vel_scaled &
                                       *(1.0_dp-0.25_dp*vel_scaled))
                       ! between 0 and 1

   velocity = vel_sign * vel_max * (0.9_dp + 0.1_dp*vel_scaled_lim)

end if

end subroutine velocity_limiter_gradual

!-------------------------------------------------------------------------------

end module calc_vxy_m
!
