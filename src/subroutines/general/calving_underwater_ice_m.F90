!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l v i n g _ u n d e r w a t e r _ i c e _ m
!
!> @file
!!
!! Calving of "underwater ice".
!!
!! @section Copyright
!!
!! Copyright 2009-2021 Ralf Greve, Thorben Dunse
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
!> Calving of "underwater ice".
!<------------------------------------------------------------------------------
module calving_underwater_ice_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m
  use error_m

  implicit none

  real(dp), dimension(0:JMAX,0:IMAX) :: calv_uw_ice

  public

contains

!-------------------------------------------------------------------------------
!> Main routine: Calving of "underwater ice".
!<------------------------------------------------------------------------------
  subroutine calving_underwater_ice(z_sl)

  implicit none

  real(dp), intent(in) :: z_sl

  real(dp)                           :: year_sec_inv
  real(dp)                           :: rhosw_rho_ratio
  real(dp)                           :: calv_uw_coeff, r1_calv_uw, r2_calv_uw
  real(dp), dimension(0:JMAX,0:IMAX) :: H, H_sea
  integer(i4b)                       :: i, j

!-------- Term abbreviations --------

  year_sec_inv = 1.0_dp/year2sec

  rhosw_rho_ratio = RHO_SW/RHO

!-------- Setting of parameters --------

#if (defined(CALV_UW_COEFF))
  calv_uw_coeff = CALV_UW_COEFF * year_sec_inv
#else
  errormsg = ' >>> calving_underwater_ice: CALV_UW_COEFF undefined!'
  call error(errormsg)
#endif

#if (defined(R1_CALV_UW))
  r1_calv_uw = R1_CALV_UW
#else
  errormsg = ' >>> calving_underwater_ice: R1_CALV_UW undefined!'
  call error(errormsg)
#endif

#if (defined(R2_CALV_UW))
  r2_calv_uw = R2_CALV_UW
#else
  errormsg = ' >>> calving_underwater_ice: R2_CALV_UW undefined!'
  call error(errormsg)
#endif

!-------- Ice thickness and sea depth --------

  H     = max(H_c + H_t, 0.0_dp)   ! ice thickness
  H_sea = max(z_sl - zl, 0.0_dp)   ! sea depth

!-------- Calving of "underwater ice" --------

#if !defined(ALLOW_OPENAD) /* Normal */

  where ( (mask == 0_i1b).and.(H < rhosw_rho_ratio*H_sea) )
     calv_uw_ice = calv_uw_coeff * H**r1_calv_uw * H_sea**r2_calv_uw
  elsewhere
     calv_uw_ice = 0.0_dp
  end where

#else /* OpenAD */

  do i=0, IMAX
  do j=0, JMAX
     if ( (mask(j,i) == 0_i1b) .and. (H(j,i) < rhosw_rho_ratio*H_sea(j,i)) ) then
        calv_uw_ice(j,i) = calv_uw_coeff * H(j,i)**r1_calv_uw * H_sea(j,i)**r2_calv_uw
     else
        calv_uw_ice(j,i) = 0.0_dp
     end if
  end do
  end do

#endif /* Normal vs. OpenAD */

  end subroutine calving_underwater_ice

!-------------------------------------------------------------------------------

end module calving_underwater_ice_m
!
