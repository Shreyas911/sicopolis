!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ s l e _ s o l v e r s
!
!> @file
!!
!! Solvers for systems of linear equations used by SICOPOLIS.
!!
!! @section Copyright
!!
!! Copyright 2009-2020 Ralf Greve, Tatsuru Sato
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
!> Solvers for systems of linear equations used by SICOPOLIS.
!<------------------------------------------------------------------------------
module sico_maths_m_stub

use sico_types_m

#if (CALCTHK!=3  && CALCTHK!=4 && CALCTHK!=6 )
  public :: sor_sprs, tri_sle, my_erfc
  !public :: sor_sprs, tri_sle, bilinint, my_erfc
#else  
  public :: tri_sle, my_erfc, sico_lis_solver
  !public :: tri_sle, bilinint, my_erfc, sico_lis_solver
#endif
#if (CALCTHK!=3 && CALCTHK!=4 && CALCTHK!=6 )
  interface sor_sprs
     module procedure sor_sprs_stub
  end interface
#endif
  interface tri_sle
     module procedure tri_sle_stub
  end interface

!  interface bilinint 
!     module procedure bilinint_stub
!  end interface

  interface my_erfc
     module procedure my_erfc_stub
  end interface

#if (CALCTHK==3 || CALCTHK==6 || MARGIN==3 || DYNAMICS==2)
!#if (CALCTHK==3 || CALCTHK==4 || CALCTHK==6 || MARGIN==3 || DYNAMICS==2)
  interface sico_lis_solver
     module procedure sico_lis_solver_stub
  end interface
#endif

contains

!-------------------------------------------------------------------------------
!> SOR solver for a system of linear equations lgs_a*lgs_x=lgs_b
!! [matrix storage: compressed sparse row CSR,
!! represented by arrays lgs_a_value(values), lgs_a_index (indices)
!! and lgs_a_ptr (pointers)].
!<------------------------------------------------------------------------------
#if (CALCTHK!=3 && CALCTHK!=4 && CALCTHK!=6 )
subroutine sor_sprs_stub(lgs_a_value, lgs_a_index, lgs_a_diag_index,&
                    lgs_a_ptr, &
                    lgs_b_value, &
                    nnz, nmax,&
#ifdef ALLOW_OPENAD 
                    n_sprs,&
#endif 
                    omega, eps_sor, lgs_x_value, ierr)
    !$openad xxx template oad_template_sor_sprs.f90
implicit none

integer(i4b),                     intent(in) :: nnz, nmax, n_sprs
real(dp),                         intent(in) :: omega, eps_sor
integer(i4b), dimension(nmax+1),  intent(inout) :: lgs_a_ptr
integer(i4b), dimension(nnz),     intent(inout) :: lgs_a_index
integer(i4b), dimension(nmax),    intent(inout) :: lgs_a_diag_index
real(dp),     dimension(nnz),     intent(inout) :: lgs_a_value
real(dp),     dimension(nmax),    intent(inout) :: lgs_b_value

integer(i4b),                    intent(out) :: ierr
real(dp),     dimension(nmax), intent(inout) :: lgs_x_value

integer(i4b) :: iter
integer(i4b) :: nr, k
real(dp), allocatable, dimension(:) :: lgs_x_value_prev
real(dp)     :: b_nr
logical      :: flag_convergence

lgs_x_value_prev (0) = lgs_x_value(0) + omega +  eps_sor + lgs_a_value(0) + lgs_b_value(0)
b_nr = lgs_x_value_prev(0)
lgs_x_value(0) = b_nr

end subroutine sor_sprs_stub
#endif

!-------------------------------------------------------------------------------
!> Solution of a system of linear equations Ax=b with tridiagonal matrix A.
!! @param[in]  a0       a0(j) is element A_(j,j-1) of Matrix A
!! @param[in]  a1       a1(j) is element A_(j,j)   of Matrix A
!! @param[in]  a2       a2(j) is element A_(j,j+1) of Matrix A
!! @param[in]  b        inhomogeneity vector
!! @param[in]  nrows    size of matrix A (indices run from 0 (!!!) to nrows)
!! @param[out] x        Solution vector.
!<------------------------------------------------------------------------------
subroutine tri_sle_stub(a0, a1, a2, x, b, nrows)

implicit none

integer(i4b),             intent(in)    :: nrows

real(dp), dimension(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), intent(in)    :: a0, a2
real(dp), dimension(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), intent(inout) :: a1, b

real(dp), dimension(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), intent(out)   :: x
real(dp), dimension(0:nrows) :: help_x
integer(i4b) :: n

!--------  Generate an upper triangular matrix
!                      ('obere Dreiecksmatrix') --------
!x(1) = a0(0) + a1(0) + a2(0) + b(0) 
!x(0) = x(1)
!b(0) = x(1)
!a1(0) = x(1)

do n=1, nrows
   a1(n)   = a1(n) - a0(n)/a1(n-1)*a2(n-1)
end do

do n=1, nrows
   b(n)    = b(n) - a0(n)/a1(n-1)*b(n-1)
!         a0(n)  = 0.0_dp , not needed in the following, therefore
!                           not set
end do

!-------- Iterative solution of the new system --------

!      x(nrows) = b(nrows)/a1(nrows)

!      do n=nrows-1, 0, -1
!         x(n) = (b(n)-a2(n)*x(n+1))/a1(n)
!      end do

help_x(0) = b(nrows)/a1(nrows)

do n=1, nrows
   help_x(n) = b(nrows-n)/a1(nrows-n) &
              -a2(nrows-n)/a1(nrows-n)*help_x(n-1)
end do

do n=0, nrows
   x(n) = help_x(nrows-n)
end do

end subroutine tri_sle_stub

!-------------------------------------------------------------------------------
!> Bilinear interpolation.
!<------------------------------------------------------------------------------
!  subroutine bilinint_stub(x1, x2, y1, y2, z11, z12, z21, z22, x, y, retval)
!
!  implicit none
!
!  real(dp), intent(in) :: x1, x2, y1, y2, z11, z12, z21, z22, x, y
!
!  real(dp) :: t, u
!  real(dp), intent(out) :: retval
!
!  real(dp), parameter :: I = 1.0_dp
!
!  t = (x-x1)/(x2-x1)
!  u = (y-y1)/(y2-y1)
!
!  retval = (I-t)*(I-u)*z11 + (I-t)*u*z12 + t*(I-u)*z21 + t*u*z22
!
!  end subroutine bilinint_stub

!-------------------------------------------------------------------------------
!> Computation of the complementary error function erfc(x) = 1-erf(x)
!! with a fractional error everywhere less than 1.2 x 10^(-7)
!! (formula by Press et al., 'Numerical Recipes in Fortran 77').
!<------------------------------------------------------------------------------
  subroutine my_erfc_stub(x, retval)

  implicit none

  real(dp), intent(in)  :: x
  real(dp), intent(out) :: retval

  real(dp) :: t, z

  z = abs(x)
  t = 1.0_dp/(1.0_dp+0.5_dp*z)

  retval = t * exp( -z*z     -1.26551223_dp &
                     + t  * (  1.00002368_dp &
                     + t  * (  0.37409196_dp &
                     + t  * (  0.09678418_dp &
                     + t  * ( -0.18628806_dp &
                     + t  * (  0.27886807_dp &
                     + t  * ( -1.13520398_dp &
                     + t  * (  1.48851587_dp &
                     + t  * ( -0.82215223_dp &
                     + t  *    0.17087277_dp ) ) ) ) ) ) ) ) )

  if (x < 0.0_dp) retval = 2.0_dp-retval

  end subroutine my_erfc_stub

!-------------------------------------------------------------------------------
#ifdef ALLOW_OPENAD   /* OAD VERSION */
  subroutine my_erfc(x, retval)

  implicit none

  real(dp), intent(in)  :: x
  real(dp), intent(out) :: retval

  real(dp) :: t, z

  z = abs(x)
  t = 1.0_dp/(1.0_dp+0.5_dp*z)

  retval = t * exp( -z*z     -1.26551223_dp &
                     + t  * (  1.00002368_dp &
                     + t  * (  0.37409196_dp &
                     + t  * (  0.09678418_dp &
                     + t  * ( -0.18628806_dp &
                     + t  * (  0.27886807_dp &
                     + t  * ( -1.13520398_dp &
                     + t  * (  1.48851587_dp &
                     + t  * ( -0.82215223_dp &
                     + t  *    0.17087277_dp ) ) ) ) ) ) ) ) )

  if (x < 0.0_dp) retval = 2.0_dp-retval

  end subroutine my_erfc
#endif /* OAD VERSION */

#if (CALCTHK==3 || CALCTHK==6 || MARGIN==3 || DYNAMICS==2)
subroutine sico_lis_solver_stub(nmax, nnz, &
                           lgs_a_ptr, lgs_a_index, &
                           lgs_a_value, lgs_b_value, lgs_x_value)

    !$openad xxx template oad_template_sico_lis_solver.f90
implicit none

integer(i4b),                   intent(in) :: nmax
integer(i4b),                   intent(in) :: nnz
integer(i4b), dimension(nmax+1),intent(inout) :: lgs_a_ptr
integer(i4b), dimension(nnz),   intent(inout) :: lgs_a_index
real(dp),     dimension(nnz),   intent(inout) :: lgs_a_value
real(dp),     dimension(nmax),  intent(inout) :: lgs_b_value
real(dp),     dimension(nmax),  intent(inout) :: lgs_x_value

!  ------ Settings for Lis

lgs_x_value(1) =  lgs_b_value(1) / lgs_a_value(1) 

end subroutine sico_lis_solver_stub
#endif

end module sico_maths_m_stub
!
