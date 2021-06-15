!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  m a s k _ u p d a t e _ s e a _ l e v e l _ m
!
!> @file
!!
!! Update of the ice-land-ocean mask due to changes of the sea level.
!!
!! @section Copyright
!!
!! Copyright 2009-2021 Ralf Greve
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
!> Update of the ice-land-ocean mask due to changes of the sea level.
!<------------------------------------------------------------------------------
module mask_update_sea_level_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m

  implicit none

  private
  public :: mask_update_sea_level

contains

!-------------------------------------------------------------------------------
!> Main function of mask_update_m:
!! Update of the ice-land-ocean mask due to changes of the sea level.
!<------------------------------------------------------------------------------
  function mask_update_sea_level(z_sl, i, j)

  implicit none

#if !defined(ALLOW_OPENAD)
  integer(i4b), intent(in) :: i, j
#else
  integer(i4b), intent(inout) :: i, j
#endif
  real(dp),     intent(in) :: z_sl

  integer(i1b) :: mask_update_sea_level
  real(dp)     :: rhosw_rho_ratio, H_ice, H_sea

  rhosw_rho_ratio = RHO_SW/RHO

!-------- Previously ice-free land point or sea point --------

  if ( (mask(j,i) == 1_i1b).or.(mask(j,i) == 2_i1b) ) then

     if (zl(j,i) > z_sl) then
        mask_update_sea_level = 1_i1b   ! now ice-free land
        return
     else
        mask_update_sea_level = 2_i1b   ! now sea point
        return
     end if

!-------- Previously grounded-ice or floating-ice point --------

  else   ! (mask(j,i) == 0_i1b, 3_i1b)

     if (zl(j,i) > z_sl) then

        mask_update_sea_level = 0_i1b   ! now grounded ice
        return

     else

        H_ice = zs(j,i)-zb(j,i)   ! ice thickness
        H_sea = z_sl   -zl(j,i)   ! sea depth

        if ( H_ice < (rhosw_rho_ratio*H_sea) ) then

#if (MARGIN==1 || (MARGIN==2 && MARINE_ICE_FORMATION==1))
           mask_update_sea_level = 2_i1b     ! ice becomes floating, therefore
                                             ! now sea point (ice cut off)
#elif (MARGIN==2 && MARINE_ICE_FORMATION==2)
           mask_update_sea_level = 0_i1b     ! now "underwater ice"
#elif (MARGIN==3)
           mask_update_sea_level = 3_i1b     ! now floating ice
#endif
           return

        else

           mask_update_sea_level = 0_i1b     ! now grounded ice
           return

        end if

     end if

  end if

  end function mask_update_sea_level

!-------------------------------------------------------------------------------

end module mask_update_sea_level_m
!
