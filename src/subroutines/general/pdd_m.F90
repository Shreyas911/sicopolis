!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  p d d _ m
!
!> @file
!!
!! Computation of the positive degree days (PDD) with statistical temperature
!! fluctuations; based on semi-analytical solution by Calov and Greve (2005).
!!
!! @section Copyright
!!
!! Copyright 2009-2019 Reinhard Calov, Ralf Greve
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
!> Computation of the positive degree days (PDD) with statistical temperature
!! fluctuations; based on semi-analytical solution by Calov and Greve (2005).
!<------------------------------------------------------------------------------
module pdd_m

  use sico_types_m
  use sico_variables_m, only : pi

  implicit none

  private
  public :: pdd

contains

!-------------------------------------------------------------------------------
!> Main subroutine of pdd_m:
!! Computation of the positive degree days (PDD) with statistical temperature
!! fluctuations; based on semi-analytical solution by Calov and Greve (2005).
!! Note that this routine uses years as time unit
!! (as opposed to the seconds which are usually used by SICOPOLIS).
!<------------------------------------------------------------------------------
  subroutine pdd(temp_mm, s_stat, ET)

#if defined(ALLOW_OPENAD) /* OpenAD */
  use sico_maths_m
#endif /* OpenAD */
  
  implicit none

  real(dp), dimension(12), intent(in) :: temp_mm
  real(dp), intent(in) :: s_stat

  real(dp), intent(out) :: ET

  integer(i4b) :: n
  real(dp) :: inv_sqrt2pi, inv_s_stat, inv_sqrt2
  real(dp) :: pdd_sum

  real(dp), parameter :: time_year     = 1.0_dp, &        ! period 1 year
                         time_year_inv = 1.0_dp/time_year, &
                         d_time        = 1.0_dp/12.0_dp   ! time-step 1 month

#if defined(ALLOW_OPENAD) /* OpenAD */
  real(dp)  :: my_erfc_retval
#endif /* OpenAD */

  inv_sqrt2pi = 1.0_dp/sqrt(2.0_dp*pi)
  inv_s_stat  = 1.0_dp/s_stat
  inv_sqrt2   = 1.0_dp/sqrt(2.0_dp)

  pdd_sum = 0.0_dp

  do n=1, 12   ! month counter

#if !defined(ALLOW_OPENAD) /* Normal */
     pdd_sum = pdd_sum &
               + ( s_stat*inv_sqrt2pi*exp(-0.5_dp*(temp_mm(n)*inv_s_stat)**2) &
                   + 0.5_dp*temp_mm(n) &
                           *erfc(-temp_mm(n)*inv_s_stat*inv_sqrt2) ) &
                 *d_time   ! positive degree days (in a * deg C)
#else /* OpenAD */
     call my_erfc(-temp_mm(n)*inv_s_stat*inv_sqrt2, my_erfc_retval)
     pdd_sum = pdd_sum &
               + ( s_stat*inv_sqrt2pi*exp(-0.5_dp*(temp_mm(n)*inv_s_stat)**2) &
                   + 0.5_dp*temp_mm(n) &
                           *my_erfc_retval ) &
                 *d_time   ! positive degree days (in a * deg C)
#endif /* Normal vs. OpenAD */

  end do

  ET   = pdd_sum*time_year_inv   ! temperature excess   (deg C)

  end subroutine pdd

!-------------------------------------------------------------------------------

end module pdd_m
!
