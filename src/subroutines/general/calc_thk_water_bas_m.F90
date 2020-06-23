!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ t h k _ w a t e r _ b a s _ m
!
!> @file
!!
!! Computation of the thickness of the water column under the ice base.
!!
!! @section Copyright
!!
!! Copyright 2009-2019 Ralf Greve
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
!> Computation of the thickness of the water column under the ice base.
!<------------------------------------------------------------------------------
module calc_thk_water_bas_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m

#if (BASAL_HYDROLOGY==1)
  use hydro_m
#endif

  implicit none

  private
  public :: calc_thk_water_bas

contains

!-------------------------------------------------------------------------------
!> Main subroutine of calc_thk_water_bas_m:
!! Computation of the thickness of the water column under the ice base.
!<------------------------------------------------------------------------------
  subroutine calc_thk_water_bas(z_sl)

  implicit none

  real(dp), intent(in) :: z_sl

#if defined(ALLOW_OPENAD) /* OpenAD */
  integer(i4b) :: i, j
#endif /* OpenAD */

  logical, save :: firstcall = .true.

#if (BASAL_HYDROLOGY==1)
  real(dp), save                     :: rho_rho_w_ratio
  integer , dimension(0:IMAX,0:JMAX) :: hydro_icemask
  real(dp), dimension(0:IMAX,0:JMAX) :: hydro_topg, hydro_thk, &
                                        hydro_temppabase, hydro_supply, &
                                        hydro_sflux, hydro_vflux, &
                                        hydro_vfluxX, hydro_vfluxY, &
                                        hydro_bwat

#if defined(ALLOW_OPENAD) /* OpenAD */
  integer(i1b), dimension(0:IMAX,0:JMAX) :: t_maske
  real(dp)    , dimension(0:IMAX,0:JMAX) :: t_H_c, t_H_t, t_Q_b_tot
#endif /* OpenAD */

  type(hydro_t), save :: hydro
                         !!! RG: Does this need a save attribute?
#endif

!-------- Water column --------

#if !defined(ALLOW_OPENAD) /* Normal */

#if (BASAL_HYDROLOGY==1)

  if (firstcall) then

     rho_rho_w_ratio = RHO/RHO_W

     call hydro_init(hydro, xi, eta)
     call hydro_gen_conf(hydro, &
          & method='quinn', &
          & avoid_frz=.false., &
          & filter_len=0.0_dp, &
          & rho_seawater=RHO_SW, &
          & rho_freshwater=RHO_W, &
          & rho_ice=RHO)

  end if

  hydro_topg       = transpose(zl)-z_sl
  hydro_temppabase = transpose(temph_b)

  where (transpose(maske)==0_i1b)   ! grounded ice
     hydro_icemask = 1
     hydro_thk     = transpose(H_c+H_t)
     hydro_supply  = rho_rho_w_ratio*transpose(Q_b_tot)
  elsewhere
     hydro_icemask = 0
     hydro_thk     = 0.0_dp
     hydro_supply  = 0.0_dp
  end where

  call hydro_set_topg(hydro, hydro_topg)
  call hydro_set_thk(hydro, hydro_thk)
  call hydro_set_temppabase(hydro, hydro_temppabase)
  call hydro_set_supply(hydro, hydro_supply)
  call hydro_set_mask(hydro, hydro_icemask)

  call hydro_update(hydro)

  call hydro_get_sflux(hydro, hydro_sflux)
  call hydro_get_vflux(hydro, hydro_vflux, hydro_vfluxX, hydro_vfluxY)
  call hydro_get_bwat(hydro, hydro_bwat)

  q_w   = transpose(hydro_vflux)
  q_w_x = transpose(hydro_vfluxX)
  q_w_y = transpose(hydro_vfluxY)
  H_w   = transpose(hydro_bwat)

#else

  where (maske==0_i1b)   ! grounded ice
     q_w   = 0.0_dp
     q_w_x = 0.0_dp
     q_w_y = 0.0_dp
     H_w   = 0.0_dp
  end where

#endif

  where (maske==2_i1b)   ! ocean
     q_w   = 0.0_dp
     q_w_x = 0.0_dp
     q_w_y = 0.0_dp
     H_w   = z_sl-zl
  elsewhere (maske==3_i1b)   ! floating ice
     q_w   = 0.0_dp
     q_w_x = 0.0_dp
     q_w_y = 0.0_dp
     H_w   = zb-zl
  elsewhere (maske==1_i1b)   ! ice-free land
     q_w   = 0.0_dp
     q_w_x = 0.0_dp
     q_w_y = 0.0_dp
     H_w   = 0.0_dp
  end where

  if (firstcall) firstcall = .false.

#else /* OpenAD */

#if (BASAL_HYDROLOGY==1)

  if (firstcall) then

     rho_rho_w_ratio = RHO/RHO_W

     call hydro_init(hydro, xi, eta)
     call hydro_gen_conf(hydro, &
          & method='quinn', &
          & avoid_frz=.false., &
          & filter_len=0.0_dp, &
          & rho_seawater=RHO_SW, &
          & rho_freshwater=RHO_W, &
          & rho_ice=RHO)

  end if

  ! hard-coding transpose
  do j=0,JMAX
  do i=0,IMAX
     hydro_topg(i,j) = zl(j,i) - z_sl
     hydro_temppabase(i,j) = temph_b(j,i)
     ! transpose these arrays for easy searching below
     t_maske(i,j) = maske(j,i)
     t_H_c(i,j) = H_c(j,i)
     t_H_t(i,j) = H_t(j,i)
     t_Q_b_tot(i,j) = Q_b_tot(j,i)
  end do
  end do

  do j=0,JMAX
  do i=0,IMAX
     if (t_maske(i,j)==0_i1b) then
        hydro_icemask(i,j) = 1
        hydro_thk(i,j)     = t_H_c(i,j) + t_H_t(i,j)
        hydro_supply(i,j)  = rho_rho_w_ratio*t_Q_b_tot(i,j)
     else
        hydro_icemask(i,j) = 0
        hydro_thk(i,j)     = 0.0_dp
        hydro_supply(i,j)  = 0.0_dp
     end if
  end do
  end do

  call hydro_set_topg(hydro, hydro_topg)
  call hydro_set_thk(hydro, hydro_thk)
  call hydro_set_temppabase(hydro, hydro_temppabase)
  call hydro_set_supply(hydro, hydro_supply)
  call hydro_set_mask(hydro, hydro_icemask)

  call hydro_update(hydro)

  call hydro_get_sflux(hydro, hydro_sflux)
  call hydro_get_vflux(hydro, hydro_vflux, hydro_vfluxX, hydro_vfluxY)
  call hydro_get_bwat(hydro, hydro_bwat)

  do i=0,IMAX
  do j=0,JMAX
     q_w(j,i)   = hydro_vflux(i,j)
     q_w_x(j,i) = hydro_vfluxX(i,j)
     q_w_y(j,i) = hydro_vfluxY(i,j)
     H_w(j,i)   = hydro_bwat(i,j)
  end do
  end do

#else

  do i=0,IMAX
  do j=0,JMAX
     if ( maske(j,i)==0_i1b ) then   ! grounded ice
        q_w(j,i)   = 0.0_dp
        q_w_x(j,i) = 0.0_dp
        q_w_y(j,i) = 0.0_dp
        H_w(j,i)   = 0.0_dp
     end if
  end do
  end do

#endif

  do i=0,IMAX
  do j=0,JMAX

     if ( maske(j,i)==2_i1b ) then
        q_w(j,i)   = 0.0_dp
        q_w_x(j,i) = 0.0_dp
        q_w_y(j,i) = 0.0_dp
        H_w(j,i)   = z_sl-zl(j,i)
     else if ( maske(j,i)==3_i1b ) then
        q_w(j,i)   = 0.0_dp
        q_w_x(j,i) = 0.0_dp
        q_w_y(j,i) = 0.0_dp
        H_w(j,i)   = zb(j,i)-zl(j,i)
     else if ( maske(j,i)==1_i1b ) then
        q_w(j,i)   = 0.0_dp
        q_w_x(j,i) = 0.0_dp
        q_w_y(j,i) = 0.0_dp
        H_w(j,i)   = 0.0_dp
     end if

  end do
  end do

  if (firstcall) firstcall = .false.

#endif /* Normal vs. OpenAD */

  end subroutine calc_thk_water_bas

!-------------------------------------------------------------------------------

end module calc_thk_water_bas_m
!
