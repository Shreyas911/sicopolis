!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ m a t h s _ m
!
!> @file
!!
!! Several mathematical tools used by SICOPOLIS.
!!
!! @section Copyright
!!
!! Copyright 2009-2022 Ralf Greve,
!!                     Liz Curry-Logan, Sri Hari Krishna Narayanan
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
!> Several mathematical tools used by SICOPOLIS.
!<------------------------------------------------------------------------------
module sico_maths_m

  use sico_types_m

  implicit none

#if !defined(ALLOW_TAPENADE)
  public 
#else

  public :: sor_sprs, tri_sle, bilinint, my_erfc

#if (CALCTHK==3 || CALCTHK==6 || MARGIN==3 || DYNAMICS==2)
  public :: sico_lis_solver
#endif

  interface sor_sprs
     module procedure sor_sprs_local
  end interface

  interface tri_sle
     module procedure tri_sle_local
  end interface

#if (CALCTHK==3 || CALCTHK==6 || MARGIN==3 || DYNAMICS==2)
  interface sico_lis_solver 
     module procedure sico_lis_solver_local
  end interface
#endif

#endif

contains

!-------------------------------------------------------------------------------
!> SOR solver for a system of linear equations lgs_a*lgs_x=lgs_b
!! [matrix storage: compressed sparse row CSR,
!! represented by arrays lgs_a_value(values), lgs_a_index (indices)
!! and lgs_a_ptr (pointers)].
!<------------------------------------------------------------------------------
#if !defined(ALLOW_TAPENADE)
  subroutine sor_sprs(lgs_a_value, lgs_a_index, lgs_a_diag_index, lgs_a_ptr, &
#else
  subroutine sor_sprs_local(lgs_a_value, lgs_a_index, lgs_a_diag_index, lgs_a_ptr, &
#endif
                      lgs_b_value, &
                      nnz, nmax, &
#if defined(ALLOW_TAPENADE)
                      n_sprs, &
#endif
                      omega, eps_sor, lgs_x_value, ierr)

  implicit none

#if defined(ALLOW_TAPENADE)
  integer(i4b),                     intent(in) :: n_sprs
#endif
  integer(i4b),                     intent(in) :: nnz, nmax
  real(dp),                         intent(in) :: omega, eps_sor
  integer(i4b), dimension(nmax+1),  intent(in) :: lgs_a_ptr
  integer(i4b), dimension(nnz),     intent(in) :: lgs_a_index
  integer(i4b), dimension(nmax),    intent(in) :: lgs_a_diag_index
  real(dp),     dimension(nnz),     intent(in) :: lgs_a_value
  real(dp),     dimension(nmax),    intent(in) :: lgs_b_value

  integer(i4b),                    intent(out) :: ierr
  real(dp),     dimension(nmax), intent(inout) :: lgs_x_value

  integer(i4b) :: iter
  integer(i4b) :: iter_max
  integer(i4b) :: nr, k
#if !defined(ALLOW_TAPENADE)
  real(dp), allocatable, dimension(:) :: lgs_x_value_prev
#else
  real(dp),           dimension(nmax) :: lgs_x_value_prev
  real(dp)     :: temp1, temp2
  logical      :: isnanflag1, isnanflag2, isnanflag3
#endif
  real(dp)     :: b_nr
  logical      :: flag_convergence

#if (ITER_MAX_SOR > 0)
  iter_max = ITER_MAX_SOR
#else
  iter_max = 1000   ! default value
#endif

#if !defined(ALLOW_TAPENADE)
  allocate(lgs_x_value_prev(nmax))
#endif

  iter_loop : do iter=1, iter_max

     lgs_x_value_prev = lgs_x_value

     do nr=1, nmax

        b_nr = 0.0_dp 

        do k=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
           b_nr = b_nr + lgs_a_value(k)*lgs_x_value(lgs_a_index(k))
        end do

        lgs_x_value(nr) = lgs_x_value(nr) &
                          -omega*(b_nr-lgs_b_value(nr)) &
                                /lgs_a_value(lgs_a_diag_index(nr))

     end do

     flag_convergence = .true.
     do nr=1, nmax
        if (abs(lgs_x_value(nr)-lgs_x_value_prev(nr)) > eps_sor) then
           flag_convergence = .false.
           exit
        end if
     end do

     if (flag_convergence) then
        write(6,'(10x,a,i0)') 'sor_sprs: iter = ', iter
        ierr = 0   ! convergence criterion fulfilled
#if !defined(ALLOW_TAPENADE)
        deallocate(lgs_x_value_prev)
#endif
        return
     end if

  end do iter_loop

  write(6,'(10x,a,i0)') 'sor_sprs: iter = ', iter
  ierr = -1   ! convergence criterion not fulfilled
#if !defined(ALLOW_TAPENADE)
  deallocate(lgs_x_value_prev)
#endif

#if !defined(ALLOW_TAPENADE)
  end subroutine sor_sprs
#else
  end subroutine sor_sprs_local
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
#if !defined(ALLOW_TAPENADE)
  subroutine tri_sle(a0, a1, a2, x, b, nrows)
#else
  subroutine tri_sle_local(a0, a1, a2, x, b, nrows)
#endif

  implicit none

  integer(i4b),                 intent(in) :: nrows
  real(dp), dimension(0:nrows), intent(in) :: a0, a1, a2, b

  real(dp), dimension(0:nrows), intent(out) :: x

  integer(i4b) :: n
  real(dp), dimension(0:nrows) :: a0_aux, a1_aux, a2_aux, b_aux, x_aux

!--------  Define local variables --------

  a0_aux = a0
  a1_aux = a1
  a2_aux = a2
  b_aux  = b
  x_aux  = 0.0_dp   ! initialization

!--------  Generate an upper triangular matrix --------

  do n=1, nrows
     a1_aux(n) = a1_aux(n) - a0_aux(n)/a1_aux(n-1)*a2_aux(n-1)
  end do

  do n=1, nrows
     b_aux(n) = b_aux(n) - a0_aux(n)/a1_aux(n-1)*b_aux(n-1)
     ! a0_aux(n) = 0.0_dp , not needed in the following, therefore not set
  end do

!-------- Iterative solution of the new system --------

  x_aux(0) = b_aux(nrows)/a1_aux(nrows)

  do n=1, nrows
     x_aux(n) = b_aux(nrows-n)/a1_aux(nrows-n) &
                -a2_aux(nrows-n)/a1_aux(nrows-n)*x_aux(n-1)
  end do

  do n=0, nrows
     x(n) = x_aux(nrows-n)
  end do

  !  WARNING: Subroutine does not check for elements of the main
  !           diagonal becoming zero. In this case it crashes even
  !           though the system may be solvable. Otherwise ok.

#if !defined(ALLOW_TAPENADE)
  end subroutine tri_sle
#else
  end subroutine tri_sle_local
#endif

!-------------------------------------------------------------------------------
!> Bilinear interpolation.
!<------------------------------------------------------------------------------
  function bilinint(x1, x2, y1, y2, z11, z12, z21, z22, x, y)

  implicit none

  real(dp), intent(in) :: x1, x2, y1, y2, z11, z12, z21, z22, x, y

  real(dp) :: t, u
  real(dp) :: bilinint

  real(dp), parameter :: I = 1.0_dp

  t = (x-x1)/(x2-x1)
  u = (y-y1)/(y2-y1)

  bilinint = (I-t)*(I-u)*z11 + (I-t)*u*z12 + t*(I-u)*z21 + t*u*z22

  end function bilinint

!-------------------------------------------------------------------------------
!> Computation of the complementary error function erfc(x) = 1-erf(x)
!! with a fractional error everywhere less than 1.2 x 10^(-7)
!! (formula by Press et al., 'Numerical Recipes in Fortran 77').
!<------------------------------------------------------------------------------
#if defined(ALLOW_TAPENADE)
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
#endif

!-------------------------------------------------------------------------------
!> Tapenade needs a template to help differentiate through the LIS solver when  
!! it is used. This is substituted in for adjoint modes in Antarctica with 
!! ice shelves.
!<------------------------------------------------------------------------------
#if (CALCTHK==3 || CALCTHK==6 || MARGIN==3 || DYNAMICS==2)
#if defined(ALLOW_TAPENADE)
#include "lisf.h"
subroutine sico_lis_solver_local(nmax, nnz, &
#else
subroutine sico_lis_solver(nmax, nnz, &
#endif
                           lgs_a_ptr, lgs_a_index, &
                           lgs_a_value, lgs_b_value, lgs_x_value)

implicit none

integer(i4b)                                 :: ierr
integer(i4b)                                 :: iter
integer(i4b)                                 :: nc, nr
integer(i4b),                     intent(in) :: nmax
integer(i4b),                     intent(in) :: nnz
integer(i4b), dimension(nmax+1),  intent(in) :: lgs_a_ptr
integer(i4b), dimension(nnz),  intent(in)    :: lgs_a_index

LIS_MATRIX                                   :: lgs_a
LIS_VECTOR                                   :: lgs_b
LIS_VECTOR                                   :: lgs_x
LIS_SOLVER                                   :: solver

real(dp),     dimension(nnz),  intent(in)    :: lgs_a_value
real(dp),     dimension(nmax),    intent(in) :: lgs_b_value
real(dp),     dimension(nmax), intent(inout) :: lgs_x_value

character(len=256)                           :: ch_solver_set_option

!  ------ Settings for Lis

call lis_matrix_create(LIS_COMM_WORLD, lgs_a, ierr)
call lis_vector_create(LIS_COMM_WORLD, lgs_b, ierr)
call lis_vector_create(LIS_COMM_WORLD, lgs_x, ierr)

call lis_matrix_set_size(lgs_a, 0, nmax, ierr)
call lis_vector_set_size(lgs_b, 0, nmax, ierr)
call lis_vector_set_size(lgs_x, 0, nmax, ierr)

do nr=1, nmax

   do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
      call lis_matrix_set_value(LIS_INS_VALUE, nr, lgs_a_index(nc), &
                                               lgs_a_value(nc), lgs_a, ierr)
   end do

   call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_b_value(nr), lgs_b, ierr)
   call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_x_value(nr), lgs_x, ierr)

end do

call lis_matrix_set_type(lgs_a, LIS_MATRIX_CSR, ierr)
call lis_matrix_assemble(lgs_a, ierr)

!  ------ Solution with Lis

call lis_solver_create(solver, ierr)

#if (defined(LIS_OPTS))
    ch_solver_set_option = trim(LIS_OPTS)
#else
    ch_solver_set_option = '-i bicgsafe -p jacobi '// &
                           '-maxiter 2000 -tol 1.0e-18 -initx_zeros false'
#endif

call lis_solver_set_option(trim(ch_solver_set_option), solver, ierr)
call CHKERR(ierr)

call lis_solve(lgs_a, lgs_b, lgs_x, solver, ierr)
call CHKERR(ierr)

call lis_solver_get_iter(solver, iter, ierr)
write(6,'(10x,a,i0)') 'calc_thk_sia_impl: iter = ', iter

lgs_x_value = 0.0_dp
call lis_vector_gather(lgs_x, lgs_x_value, ierr)
call lis_matrix_destroy(lgs_a, ierr)
call lis_vector_destroy(lgs_b, ierr)
call lis_vector_destroy(lgs_x, ierr)
call lis_solver_destroy(solver, ierr)

#if defined(ALLOW_TAPENADE)
end subroutine sico_lis_solver_local
#else
end subroutine sico_lis_solver
#endif
#endif
!-------------------------------------------------------------------------------

end module sico_maths_m
!
