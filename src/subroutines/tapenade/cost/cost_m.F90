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

  public :: cost_final 

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
  
  integer(i4b) :: i, j, k, kc, kt, ios, KDATA
  character(len=64), parameter :: thisroutine = 'cost_final'

  !-------- Calculate the difference between the modeled and 'observed' ages:
fc = 0.0

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
        ! only counting points that are real in the data: 
        if (  age_data(k,j,i) .ge. -0.5 .and. age_data(k,j,i) .le. 60000.0 .and. H_BedMachine_data(j,i) .ge. 2000.0) then
          fc = fc &
#ifdef ALLOW_AGE_UNCERT
          + (age_data(k,j,i) - age_c(k,j,i)/year2sec)**2/age_unc_data(k,j,i)**2
#else
          + (age_data(k,j,i) - age_c(k,j,i)/year2sec)**2
#endif
        end if
      end do
    end do
  end do

#endif

#if (defined(BEDMACHINE_COST) || defined(FAKE_BEDMACHINE_COST))
    do i=0, IMAX
      do j=0, JMAX
        fc = fc &
#ifdef ALLOW_BEDMACHINE_UNCERT
        + (H(j,i) - H_BedMachine_data(j,i))**2/H_unc_BedMachine_data(j,i)**2
#else
        + (H(j,i) - H_BedMachine_data(j,i))**2
#endif
      end do
    end do
#endif

#if (!defined(BEDMACHINE_COST) && !defined(AGE_COST) && !defined(FAKE_BEDMACHINE_COST) && !defined(FAKE_AGE_COST))
    do i=0, IMAX
      do j=0, JMAX
        !--- Other cost functions:
        fc = fc + (H_c(j,i) + H_t(j,i))*cell_area(j,i)
      end do
    end do
#endif
  
  !-------- Print to screen just in case something gets
  !         crazy with the file outputting:
  print *, 'Final cost, fc = ', fc
  print *, trim(OUT_PATH)

  end subroutine cost_final

end module cost_m
