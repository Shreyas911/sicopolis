!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
!   Module  :  c t r l
  
!   Purpose :  Declarations of control variables for adjointing

!   Copyright 2017-2022 Liz Curry-Logan, Shreyas Sunil Gaikwad,

!                       Sri Hari Krishna Narayanan
  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
module ctrl_m

  use sico_types_m  
  use sico_variables_m
  use sico_vars_m

  implicit none

  public :: cost_final 
#if (defined(ALLOW_GRDCHK) || defined(ALLOW_TAPENADE))
  public :: myceiling, myfloor 
#endif

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
!<------------------------------------------------------------------------------
  subroutine cost_final()
  
  implicit none
  
  integer(i4b) :: i, j, k, kc, kt, ios, KDATA 
  
  !-------- Calculate the difference between the modeled and 'observed' ages:
fc = 0.0

#ifdef AGE_COST

#if (CALCMOD!=1)
  KDATA = KCMAX
#else 
  KDATA = KCMAX + KTMAX
print *, '>>> error: CALCMOD == 1 but final cost not properly working for '
print *, '           AGE_COST simulations'
#endif
 
  do k=0, KDATA 
    do i=0, IMAX
      do j=0, JMAX

        ! only counting points that are real in the data: 
        if (  age_data(k,j,i).lt. 0.0       &
         .and.age_unc(k,j,i) .lt. 0.0    ) then

          fc = fc + &
                      sqrt(((age_data(k,j,i) - age_c(k,j,i)))**2 &
                      / ((age_unc(k,j,i))**2))

        end if
      end do
    end do
  end do

#elif defined(BEDMACHINE_COST)
    do i=0, IMAX
      do j=0, JMAX
        !--- Other cost functions:
        fc = fc &
        + (H(j,i) - H_BedMachine_data(j,i))**2/H_unc_BedMachine_data(j,i)**2
      end do
    end do

#else
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

#if (defined(ALLOW_GRDCHK) || defined(ALLOW_TAPENADE))
  subroutine myceiling(num, output)
  
  implicit none
  
  integer(i4b)             :: inum
  real(dp), intent(in) :: num
  integer(i4b), intent(out) :: output
  
  inum = int(num);
  
  inum = int(num);
  if (abs(num-real(inum,dp)) <= &
    (abs(num+real(inum,dp))*myepsilon_dp) ) then
    output = inum;
  else if (num>0) then
    output = inum + 1;
  else if (num<0) then
    output = inum;
  end if
  end subroutine myceiling
  
  subroutine myfloor(num, output)
  
  implicit none
  
  integer(i4b)             :: inum
  real(dp), intent(in) :: num
  integer(i4b), intent(out) :: output
  
  inum = int(num);
  
  if (abs(num-real(inum,dp)) <= &
    (abs(num+real(inum,dp))*myepsilon_dp) ) then
    output = inum;
  else if (num>0) then
    output = inum;
  else if (num<0) then
    output = inum - 1;
  end if
  
  end subroutine myfloor
#endif

end module ctrl_m
