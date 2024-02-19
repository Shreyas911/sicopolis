!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  g l o b a l s
!
!! TAPENADE does not like the save statement and instead prefers the global
!! variables to be declared in a separate file.
!!
!!##### Authors
!!
!! Ralf Greve, Shreyas Sunil Gaikwad, Liz Curry-Logan
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
!> TAPENADE does not like the save statement and instead prefers the global
!! variables to be declared in a separate file.
!-++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module globals
        
  use sico_types_m

  implicit none

  integer(i4b) :: just_a_dummy_variable

end module globals
!
