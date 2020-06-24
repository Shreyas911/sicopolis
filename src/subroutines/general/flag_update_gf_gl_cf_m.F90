!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  f l a g _ u p d a t e _ g f _ g l _ c f _ m
!
!> @file
!!
!! Update of the flags for the land-terminating grounded front,
!! marine-terminating grounded front, grounding line and calving front.
!!
!! @section Copyright
!!
!! Copyright 2018-2020 Ralf Greve
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
!> Update of the flags for the land-terminating grounded front,
!! marine-terminating grounded front, grounding line and calving front.
!<------------------------------------------------------------------------------
module flag_update_gf_gl_cf_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m

  implicit none

  private
  public :: flag_update_gf_gl_cf

contains

!-------------------------------------------------------------------------------
!> Main subroutine of flag_update_gf_gl_cf_m:
!! Update of the flags for the land-terminating grounded front,
!! marine-terminating grounded front, grounding line and calving front.
!<------------------------------------------------------------------------------
  subroutine flag_update_gf_gl_cf()

  implicit none

  integer(i4b) :: i, j

  flag_grounded_front_a_1 = .false.
  flag_grounded_front_a_2 = .false.

  flag_grounded_front_b_1 = .false.
  flag_grounded_front_b_2 = .false.

  flag_grounding_line_1 = .false.
  flag_grounding_line_2 = .false.

  flag_calving_front_1 = .false.
  flag_calving_front_2 = .false.

  do i=1, IMAX-1
  do j=1, JMAX-1

     if ( (maske(j,i)==0_i1b) &   ! grounded ice point
          .and. &
            (    (maske(j,i+1)==1_i1b)   &   ! with
             .or.(maske(j,i-1)==1_i1b)   &   ! one
             .or.(maske(j+1,i)==1_i1b)   &   ! neighbouring
             .or.(maske(j-1,i)==1_i1b) ) &   ! ice-free land point
        ) &
     flag_grounded_front_a_1(j,i) = .true.

     if ( (maske(j,i)==1_i1b) &   ! ice-free land point
          .and. &
            (    (maske(j,i+1)==0_i1b)   &   ! with
             .or.(maske(j,i-1)==0_i1b)   &   ! one
             .or.(maske(j+1,i)==0_i1b)   &   ! neighbouring
             .or.(maske(j-1,i)==0_i1b) ) &   ! grounded ice point
        ) &
     flag_grounded_front_a_2(j,i) = .true.

     if ( (maske(j,i)==0_i1b) &   ! grounded ice point
          .and. &
            (    (maske(j,i+1)==2_i1b)   &   ! with
             .or.(maske(j,i-1)==2_i1b)   &   ! one
             .or.(maske(j+1,i)==2_i1b)   &   ! neighbouring
             .or.(maske(j-1,i)==2_i1b) ) &   ! ocean point
        ) &
     flag_grounded_front_b_1(j,i) = .true.

     if ( (maske(j,i)==2_i1b) &   ! ocean point
          .and. &
            (    (maske(j,i+1)==0_i1b)   &   ! with
             .or.(maske(j,i-1)==0_i1b)   &   ! one
             .or.(maske(j+1,i)==0_i1b)   &   ! neighbouring
             .or.(maske(j-1,i)==0_i1b) ) &   ! grounded ice point
        ) &
     flag_grounded_front_b_2(j,i) = .true.

#if (MARGIN==3)

     if ( (maske(j,i)==0_i1b) &   ! grounded ice point
          .and. &
            (    (maske(j,i+1)==3_i1b)   &   ! with
             .or.(maske(j,i-1)==3_i1b)   &   ! one
             .or.(maske(j+1,i)==3_i1b)   &   ! neighbouring
             .or.(maske(j-1,i)==3_i1b) ) &   ! floating ice point
        ) &
     flag_grounding_line_1(j,i) = .true.

     if ( (maske(j,i)==3_i1b) &   ! floating ice point
          .and. &
            (    (maske(j,i+1)==0_i1b)   &   ! with
             .or.(maske(j,i-1)==0_i1b)   &   ! one
             .or.(maske(j+1,i)==0_i1b)   &   ! neighbouring
             .or.(maske(j-1,i)==0_i1b) ) &   ! grounded ice point
        ) &
     flag_grounding_line_2(j,i) = .true.

     if ( (maske(j,i)==3_i1b) &   ! floating ice point
          .and. &
            (    (maske(j,i+1)==2_i1b)   &   ! with
             .or.(maske(j,i-1)==2_i1b)   &   ! one
             .or.(maske(j+1,i)==2_i1b)   &   ! neighbouring
             .or.(maske(j-1,i)==2_i1b) ) &   ! ocean point
        ) &
     flag_calving_front_1(j,i) = .true.

     if ( (maske(j,i)==2_i1b) &   ! ocean point
          .and. &
            (    (maske(j,i+1)==3_i1b)   &   ! with
             .or.(maske(j,i-1)==3_i1b)   &   ! one
             .or.(maske(j+1,i)==3_i1b)   &   ! neighbouring
             .or.(maske(j-1,i)==3_i1b) ) &   ! floating ice point
        ) &
     flag_calving_front_2(j,i) = .true.

#endif   /* MARGIN==3 */

  end do
  end do

  do i=1, IMAX-1

     j=0

     if ( (maske(j,i)==1_i1b) &   ! ice-free land point
          .and. (maske(j+1,i)==0_i1b) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_a_2(j,i) = .true.

     if ( (maske(j,i)==2_i1b) &   ! ocean point
          .and. (maske(j+1,i)==0_i1b) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_b_2(j,i) = .true.

#if (MARGIN==3)

     if ( (maske(j,i)==2_i1b) &   ! ocean point
          .and. (maske(j+1,i)==3_i1b) &   ! with one neighbouring
        ) &                               ! floating ice point
     flag_calving_front_2(j,i) = .true.

#endif   /* MARGIN==3 */

     j=JMAX

     if ( (maske(j,i)==1_i1b) &   ! ice-free land point
          .and. (maske(j-1,i)==0_i1b) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_a_2(j,i) = .true.

     if ( (maske(j,i)==2_i1b) &   ! ocean point
          .and. (maske(j-1,i)==0_i1b) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_b_2(j,i) = .true.

#if (MARGIN==3)

     if ( (maske(j,i)==2_i1b) &   ! ocean point
          .and. (maske(j-1,i)==3_i1b) &   ! with one neighbouring
        ) &                               ! floating ice point
     flag_calving_front_2(j,i) = .true.

#endif   /* MARGIN==3 */

  end do

  do j=1, JMAX-1

     i=0

     if ( (maske(j,i)==1_i1b) &   ! ice-free land point
          .and. (maske(j,i+1)==0_i1b) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_a_2(j,i) = .true.

     if ( (maske(j,i)==2_i1b) &   ! ocean point
          .and. (maske(j,i+1)==0_i1b) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_b_2(j,i) = .true.

#if (MARGIN==3)

     if ( (maske(j,i)==2_i1b) &   ! ocean point
          .and. (maske(j,i+1)==3_i1b) &   ! with one neighbouring
        ) &                               ! floating ice point
     flag_calving_front_2(j,i) = .true.

#endif   /* MARGIN==3 */

     i=IMAX

     if ( (maske(j,i)==1_i1b) &   ! ice-free land point
          .and. (maske(j,i-1)==0_i1b) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_a_2(j,i) = .true.

     if ( (maske(j,i)==2_i1b) &   ! ocean point
          .and. (maske(j,i-1)==0_i1b) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_b_2(j,i) = .true.

#if (MARGIN==3)

     if ( (maske(j,i)==2_i1b) &   ! ocean point
          .and. (maske(j,i-1)==3_i1b) &   ! with one neighbouring
        ) &                               ! floating ice point
     flag_calving_front_2(j,i) = .true.

#endif   /* MARGIN==3 */

  end do

  end subroutine flag_update_gf_gl_cf

!-------------------------------------------------------------------------------

end module flag_update_gf_gl_cf_m
!
