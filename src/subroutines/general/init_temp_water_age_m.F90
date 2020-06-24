!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  i n i t _ t e m p _ w a t e r _ a g e _ m
!
!> @file
!!
!! Initial temperature, water content and age.
!!
!! @section Copyright
!!
!! Copyright 2009-2020 Ralf Greve, Thorben Dunse
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
!> Initial temperature, water content and age.
!<------------------------------------------------------------------------------
module init_temp_water_age_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m

  implicit none

  private
  public :: init_temp_water_age_1_1, init_temp_water_age_1_2
  public :: init_temp_water_age_1_3, init_temp_water_age_1_4
  public :: init_temp_water_age_1_5
  public :: init_temp_water_age_2

contains

!-------------------------------------------------------------------------------
!> Initial temperature, water content and age
!! (case ANF_DAT==1, TEMP_INIT==1:
!! present-day initial topography, isothermal conditions).
!<------------------------------------------------------------------------------
  subroutine init_temp_water_age_1_1()

  implicit none

  integer(i4b) :: i, j, kc

!-------- Initial ice temperature --------

  do i=0, IMAX
  do j=0, JMAX

#if (defined(NMARS) || defined(SMARS))   /* Polar caps of Mars */

     do kc=0, KCMAX
        temp_c(kc,j,i) = -100.0_dp
     end do

#else   /* all other domains */

     do kc=0, KCMAX
        temp_c(kc,j,i) =  -10.0_dp
     end do

#endif

  end do
  end do

!-------- Initial lithosphere temperature, water content and age --------

  call init_temp_r()
  call init_water()
  call init_age()

  end subroutine init_temp_water_age_1_1

!-------------------------------------------------------------------------------
!> Initial temperature, water content and age
!! (case ANF_DAT==1, TEMP_INIT==2:
!! present-day initial topography,
!! ice temperature equal to local surface temperature).
!<------------------------------------------------------------------------------
  subroutine init_temp_water_age_1_2()

  implicit none

  integer(i4b) :: i, j, kc

!-------- Initial ice temperature --------

  do i=0, IMAX
  do j=0, JMAX

     do kc=0, KCMAX
        temp_c(kc,j,i) = temp_s(j,i)
     end do

  end do
  end do

!-------- Initial lithosphere temperature, water content and age --------

  call init_temp_r()
  call init_water()
  call init_age()

  end subroutine init_temp_water_age_1_2

!-------------------------------------------------------------------------------
!> Initial temperature, water content and age
!! (case ANF_DAT==1, TEMP_INIT==3:
!! present-day initial topography,
!! ice temperature linearly increasing with depth).
!<------------------------------------------------------------------------------
  subroutine init_temp_water_age_1_3()

  use ice_material_properties_m, only : kappa_val

  implicit none

  integer(i4b) :: i, j, kc
  real(dp)     :: kappa_const_val
  real(dp)     :: temp_ice_base

!-------- Initial ice temperature --------

#if (defined(NMARS) || defined(SMARS))   /* Polar caps of Mars */
  kappa_const_val = kappa_val(-100.0_dp)
#else   /* all other domains */
  kappa_const_val = kappa_val(-10.0_dp)
#endif

  do i=0, IMAX
  do j=0, JMAX

     if (maske(j,i)<=2_i1b) then

        do kc=0, KCMAX

           temp_c(kc,j,i) = temp_s(j,i) &
                            + (q_geo(j,i)/kappa_const_val) &
                              *H_c(j,i)*(1.0_dp-eaz_c_quotient(kc))
                            ! linear temperature distribution according to the
                            ! geothermal heat flux
        end do

        if (temp_c(0,j,i) >= -BETA*H_c(j,i)) then

           temp_ice_base = -BETA*H_c(j,i)

           do kc=0, KCMAX
              temp_c(kc,j,i) = temp_s(j,i) &
                               + (temp_ice_base-temp_s(j,i)) &
                                 *(1.0_dp-eaz_c_quotient(kc))
           end do

        end if

     else   ! maske(j,i)==3_i1b, floating ice

        temp_ice_base = -BETA*H_c(j,i) - DELTA_TM_SW

        do kc=0, KCMAX
           temp_c(kc,j,i) = temp_s(j,i) &
                            + (temp_ice_base-temp_s(j,i)) &
                              *(1.0_dp-eaz_c_quotient(kc))
        end do

     end if

  end do
  end do

!-------- Initial lithosphere temperature, water content and age --------

  call init_temp_r()
  call init_water()
  call init_age()

  end subroutine init_temp_water_age_1_3

!-------------------------------------------------------------------------------
!> Initial temperature, water content and age
!! (case ANF_DAT==1, TEMP_INIT==4:
!! present-day initial topography, ice temperature from Robin (1955) solution).
!<------------------------------------------------------------------------------
  subroutine init_temp_water_age_1_4()

  use ice_material_properties_m, only : kappa_val, c_val

  implicit none

  integer(i4b) :: i, j, kc
  real(dp)     :: kappa_const_val, c_const_val
  real(dp)     :: as_val, H_val, qgeo_val, K, z_above_base
  real(dp)     :: erf_val_1, erf_val_2
  real(dp)     :: temp_ice_base, temp_scale_factor

!-------- Initial ice temperature --------

#if (defined(NMARS) || defined(SMARS))   /* Polar caps of Mars */
  kappa_const_val = kappa_val(-100.0_dp)
      c_const_val =     c_val(-100.0_dp)
#else   /* all other domains */
  kappa_const_val = kappa_val(-10.0_dp)
      c_const_val =     c_val(-10.0_dp)
#endif

  do i=0, IMAX
  do j=0, JMAX

     if (maske(j,i)<=2_i1b) then
        as_val = max(as_perp(j,i), epsi)
     else   ! maske(j,i)==3_i1b, floating ice
        as_val = epsi   ! this will produce an almost linear temperature profile
     end if

     H_val    = max(H_c(j,i)  , eps)
     qgeo_val = max(q_geo(j,i), eps)

     K = sqrt( (kappa_const_val/(RHO*c_const_val)) * (H_val/as_val) )

     erf_val_1 = erf(H_c(j,i)/(sqrt(2.0_dp)*K))

     do kc=0, KCMAX
        z_above_base   = H_c(j,i)*eaz_c_quotient(kc)
        erf_val_2      = erf(z_above_base/(sqrt(2.0_dp)*K))
        temp_c(kc,j,i) = temp_s(j,i) &
                          + (qgeo_val/kappa_const_val) &
                            * sqrt(0.5_dp*pi)*K*(erf_val_1-erf_val_2)
     end do

     if ( (maske(j,i) <= 2_i1b).and.(temp_c(0,j,i) >= -BETA*H_c(j,i)) ) then
        temp_ice_base     = -BETA*H_c(j,i)
        temp_scale_factor = (temp_ice_base-temp_s(j,i)) &
                                  /(temp_c(0,j,i)-temp_s(j,i))
     else if (maske(j,i) == 3_i1b) then
        temp_ice_base     = -BETA*H_c(j,i)-DELTA_TM_SW
        temp_scale_factor = (temp_ice_base-temp_s(j,i)) &
                                /(temp_c(0,j,i)-temp_s(j,i))
     else
        temp_scale_factor = 1.0_dp
     end if

     do kc=0, KCMAX
        temp_c(kc,j,i) = temp_s(j,i) &
                         + temp_scale_factor*(temp_c(kc,j,i)-temp_s(j,i))
     end do

  end do
  end do

!-------- Initial lithosphere temperature, water content and age --------

  call init_temp_r()
  call init_water()
  call init_age()

  end subroutine init_temp_water_age_1_4

!-------------------------------------------------------------------------------
!> Initial temperature, water content and age
!! (case ANF_DAT==1, TEMP_INIT==5:
!! present-day initial topography,
!! ice temperature, water content and age from previous simulation).
!<------------------------------------------------------------------------------
  subroutine init_temp_water_age_1_5(z_sl, filename)

  use read_m, only : read_erg_nc

  implicit none

  character(len=100), intent(in)    :: filename
  real(dp),           intent(inout) :: z_sl

  integer(i4b) :: i, j, kc, kt, kr
  real(dp)     :: temp_ice_base, temp_scale_factor

  integer(i1b), dimension(0:JMAX,0:IMAX)     :: maske_aux, n_cts_aux
  integer(i4b), dimension(0:JMAX,0:IMAX)     :: kc_cts_aux
  real(dp), dimension(0:JMAX,0:IMAX)         :: H
  real(dp), dimension(0:JMAX,0:IMAX)         :: H_cold_aux, H_temp_aux, H_aux
  real(dp), dimension(0:KRMAX,0:JMAX,0:IMAX) :: temp_r_aux
  real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: omega_t_aux, age_t_aux
  real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: temp_c_aux, age_c_aux
  real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: omega_c_aux

  H = H_c + H_t

  call init_temp_water_age_1_4()

  call read_erg_nc(z_sl, filename, &
                   opt_maske   = maske_aux   , &
                   opt_n_cts   = n_cts_aux   , &
                   opt_kc_cts  = kc_cts_aux  , &
                   opt_H_cold  = H_cold_aux  , &
                   opt_H_temp  = H_temp_aux  , &
                   opt_H       = H_aux       , &
                   opt_temp_r  = temp_r_aux  , &
                   opt_omega_t = omega_t_aux , &
                   opt_age_t   = age_t_aux   , &
                   opt_temp_c  = temp_c_aux  , &
                   opt_age_c   = age_c_aux   , &
                   opt_omega_c = omega_c_aux , &
                   opt_flag_temp_age_only = .true.)

  do i=0, IMAX
  do j=0, JMAX

     if ( ((maske(j,i)==0_i1b).or.(maske(j,i)==3_i1b)) &
          .and. &
          ((maske_aux(j,i)==0_i1b).or.(maske_aux(j,i)==3_i1b)) ) then

        n_cts(j,i)  = n_cts_aux(j,i)
        kc_cts(j,i) = kc_cts_aux(j,i)

#if (CALCMOD==1)
        H_t(j,i) = H(j,i) &
                   * min( &
                          max(H_temp_aux(j,i),0.0_dp)/max(H_aux(j,i),eps_dp), &
                          1.0_dp &
                        )
        H_c(j,i) = H(j,i)  - H_t(j,i)
        zm(j,i)  = zb(j,i) + H_t(j,i)
#endif

        temp_r(:,j,i)  = temp_r_aux(:,j,i)

        omega_t(:,j,i) = omega_t_aux(:,j,i)
        age_t(:,j,i)   = age_t_aux(:,j,i)

        temp_c(:,j,i)  = temp_c_aux(:,j,i)
        age_c(:,j,i)   = age_c_aux(:,j,i)
        omega_c(:,j,i) = omega_c_aux(:,j,i)

     end if

     if ( (maske(j,i)==3_i1b).and.(maske_aux(j,i)==0_i1b) ) then
        ! correction for ice shelves

        n_cts(j,i)  = -1_i1b
        kc_cts(j,i) = 0

        H_t(j,i) = 0.0_dp
        H_c(j,i) = H(j,i)
        zm(j,i)  = zb(j,i)

        omega_t(:,j,i) = 0.0_dp
        omega_c(:,j,i) = 0.0_dp

        temp_ice_base     = -BETA*H_c(j,i)-DELTA_TM_SW
        temp_scale_factor = (temp_ice_base-temp_s(j,i)) &
                                /(temp_c(0,j,i)-temp_s(j,i))

        do kc=0, KCMAX
           temp_c(kc,j,i) = temp_s(j,i) &
                            + temp_scale_factor*(temp_c(kc,j,i)-temp_s(j,i))
        end do


     end if

  end do
  end do

  end subroutine init_temp_water_age_1_5

!-------------------------------------------------------------------------------
!> Initial temperature, water content and age
!! (case ANF_DAT==2: ice-free conditions with relaxed bedrock).
!<------------------------------------------------------------------------------
  subroutine init_temp_water_age_2()

  implicit none

  integer(i4b) :: i, j, kc

!-------- Initial ice temperature --------

  do i=0, IMAX
  do j=0, JMAX

     do kc=0, KCMAX
        temp_c(kc,j,i) = temp_s(j,i)
     end do

  end do
  end do

!-------- Initial lithosphere temperature, water content and age --------

  call init_temp_r()
  call init_water()
  call init_age()

  end subroutine init_temp_water_age_2

!-------------------------------------------------------------------------------
!> Initial lithosphere temperature.
!<------------------------------------------------------------------------------
  subroutine init_temp_r()

  implicit none

  integer(i4b) :: i, j, kr

  do i=0, IMAX
  do j=0, JMAX

     do kr=0, KRMAX
        temp_r(kr,j,i) = temp_c(0,j,i) &
                         + (q_geo(j,i)/KAPPA_R) &
                           *H_R*(1.0_dp-zeta_r(kr))
              ! linear temperature distribution according to the
              ! geothermal heat flux
     end do

  end do
  end do

  end subroutine init_temp_r

!-------------------------------------------------------------------------------
!> Initial water content.
!<------------------------------------------------------------------------------
  subroutine init_water()

  implicit none

  omega_c = 0.0_dp   ! only required for the enthalpy method
  omega_t = 0.0_dp

  end subroutine init_water

!-------------------------------------------------------------------------------
!> Initial age.
!<------------------------------------------------------------------------------
  subroutine init_age()

  implicit none

#if (defined(ASF))   /* Austfonna */

  age_c =  3500.0_dp*year2sec
  age_t =  3500.0_dp*year2sec

#elif (defined(NMARS) || defined(SMARS))   /* Polar caps of Mars */

  age_c =  1.0e+06_dp*year2sec
  age_t =  1.0e+06_dp*year2sec

#else   /* all other domains */

  age_c = 15000.0_dp*year2sec
  age_t = 15000.0_dp*year2sec

#endif

  end subroutine init_age

!-------------------------------------------------------------------------------

end module init_temp_water_age_m
!
