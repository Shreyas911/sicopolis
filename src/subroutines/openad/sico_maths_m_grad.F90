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
module sico_maths_m_grad

use sico_types_m

  private:: transpose_csr
  public :: sor_spr_grad, tri_sle_grad
#ifdef BUILD_LIS
  public :: sico_lis_solver_grad
#endif

  interface sor_sprs_grad
     module procedure sor_sprs_grad_local
  end interface

  interface tri_sle_grad
     module procedure tri_sle_grad_local
  end interface

#ifdef BUILD_LIS
  interface sico_lis_solver_grad
     module procedure sico_lis_solver_grad_local
  end interface
#endif
contains

subroutine transpose_csr(a_value, a_index, a_diag_index, a_ptr, &
                    nnz, nmax, &
                    b_value, b_index, b_diag_index, b_ptr)

implicit none

integer(i4b),                     intent(in)  :: nnz, nmax
integer(i4b), dimension(nmax+1),  intent(in)  :: a_ptr
integer(i4b), dimension(nnz),     intent(in)  :: a_index
integer(i4b), dimension(nmax),    intent(in)  :: a_diag_index
real(dp),     dimension(nnz),     intent(in)  :: a_value
integer(i4b), dimension(nmax+1),  intent(out) :: b_ptr
integer(i4b), dimension(nnz),     intent(out) :: b_index
integer(i4b), dimension(nmax),    intent(out) :: b_diag_index
real(dp),     dimension(nnz),     intent(out) :: b_value

integer(i4b) :: nr, k
integer(i4b), dimension(nmax)                 :: b_ptr_plus

   b_ptr = 0
   do nr=1, nmax
      do k=a_ptr(nr), a_ptr(nr+1)-1
         b_ptr(a_index(k)+1) = b_ptr(a_index(k)+1) + 1
      end do
   end do
   b_ptr(1) = 1
   do nr=1, nmax
     b_ptr(nr+1) = b_ptr(nr) + b_ptr(nr+1)
   end do
   b_ptr_plus = 0
   do nr=1, nmax
      do k=a_ptr(nr), a_ptr(nr+1)-1
         b_index(b_ptr(a_index(k)) + b_ptr_plus(a_index(k))) = nr
         if(a_index(k).eq.nr)then
           b_diag_index(nr) = b_ptr(a_index(k)) + b_ptr_plus(a_index(k))
         end if
         b_value(b_ptr(a_index(k)) + b_ptr_plus(a_index(k))) = a_value(k)
         b_ptr_plus(a_index(k)) = b_ptr_plus(a_index(k)) + 1
      end do
   end do
end subroutine transpose_csr

subroutine transpose_csr_no_diagonal(a_value, a_index, a_ptr, &
                    nnz, nmax, &
                    b_value, b_index, b_ptr)

implicit none

integer(i4b),                     intent(in)  :: nnz, nmax
integer(i4b), dimension(nmax+1),  intent(in)  :: a_ptr
integer(i4b), dimension(nnz),     intent(in)  :: a_index
real(dp),     dimension(nnz),     intent(in)  :: a_value
integer(i4b), dimension(nmax+1),  intent(out) :: b_ptr
integer(i4b), dimension(nnz),     intent(out) :: b_index
real(dp),     dimension(nnz),     intent(out) :: b_value

integer(i4b) :: nr, k
integer(i4b), dimension(nmax)                 :: b_ptr_plus

   b_ptr = 0
   do nr=1, nmax
      do k=a_ptr(nr), a_ptr(nr+1)-1
         b_ptr(a_index(k)+1) = b_ptr(a_index(k)+1) + 1
      end do
   end do
   b_ptr(1) = 1
   do nr=1, nmax
     b_ptr(nr+1) = b_ptr(nr) + b_ptr(nr+1)
   end do
   b_ptr_plus = 0
   do nr=1, nmax
      do k=a_ptr(nr), a_ptr(nr+1)-1
         b_index(b_ptr(a_index(k)) + b_ptr_plus(a_index(k))) = nr
         b_value(b_ptr(a_index(k)) + b_ptr_plus(a_index(k))) = a_value(k)
         b_ptr_plus(a_index(k)) = b_ptr_plus(a_index(k)) + 1
      end do
   end do
end subroutine transpose_csr_no_diagonal
!-------------------------------------------------------------------------------
!> SOR solver for a system of linear equations lgs_a*lgs_x=lgs_b
!! [matrix storage: compressed sparse row CSR,
!! represented by arrays lgs_a_value(values), lgs_a_index (indices)
!! and lgs_a_ptr (pointers)].
!<------------------------------------------------------------------------------
subroutine sor_sprs_grad_local(lgs_a_value, lgs_a_value_b, lgs_a_index, lgs_a_diag_index, lgs_a_ptr, &
                    lgs_b_value, lgs_b_value_b, &
                    nnz, nmax, n_sprs, omega, eps_sor, lgs_x_value, lgs_x_value_b, ierr)

use sico_maths_m, only : sor_sprs
!use sico_maths_m, only : sor_sprs_local
implicit none

integer(i4b),                     intent(in) :: nnz, nmax, n_sprs
real(dp),                         intent(in) :: omega, eps_sor
integer(i4b), dimension(nmax+1),  intent(in) :: lgs_a_ptr
integer(i4b), dimension(nnz),     intent(in) :: lgs_a_index
integer(i4b), dimension(nmax),    intent(in) :: lgs_a_diag_index
real(dp),     dimension(nnz),     intent(in) :: lgs_a_value
real(dp),     dimension(nmax),    intent(in) :: lgs_b_value

integer(i4b),                    intent(out) :: ierr
real(dp),     dimension(nmax), intent(inout) :: lgs_x_value

real(dp),     dimension(nnz),  intent(inout) :: lgs_a_value_b
real(dp),     dimension(nmax), intent(inout) :: lgs_b_value_b
real(dp),     dimension(nmax), intent(inout) :: lgs_x_value_b

integer(i4b), dimension(nmax+1)              :: lgs_aT_ptr
integer(i4b), dimension(nnz)                 :: lgs_aT_index
integer(i4b), dimension(nmax)                :: lgs_aT_diag_index
real(dp),     dimension(nnz)                 :: lgs_aT_value
real(dp),     dimension(nmax)                :: incrbb
integer(i4b) :: iter
integer(i4b) :: nr, k
real(dp), allocatable, dimension(:) :: lgs_x_value_prev
real(dp)     :: b_nr
logical      :: flag_convergence
logical      :: inif, infor

  incrbb = 0.0
  call transpose_csr(lgs_a_value, lgs_a_index, lgs_a_diag_index, lgs_a_ptr, &
                    nnz, nmax, &
                    lgs_aT_value, lgs_aT_index, lgs_aT_diag_index, lgs_aT_ptr)

  call sor_sprs(lgs_aT_value, lgs_aT_index, lgs_aT_diag_index, lgs_aT_ptr, &
                    lgs_x_value_b, &
                    nnz, nmax, n_sprs, omega, eps_sor, incrbb, ierr)

  DO nr=1,nmax
    lgs_b_value_b(nr) = lgs_b_value_b(nr) + incrbb(nr)
  ENDDO

  call sor_sprs(lgs_a_value, lgs_a_index, lgs_a_diag_index, lgs_a_ptr, &
                    lgs_b_value, &
                    nnz, nmax, n_sprs, omega, eps_sor, lgs_x_value, ierr)
  do nr=1, nmax
    inif = .false. 
    infor = .false. 
    if(nr .eq. lgs_a_diag_index(nr)) then
      inif = .true.
      lgs_a_value_b(lgs_a_diag_index(nr)) = lgs_a_value_b(lgs_a_diag_index(nr)) -&
             lgs_x_value(lgs_a_index(nr)) * incrbb(nr)
    end if
    do k=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
      infor = .true.
      lgs_a_value_b(k) = lgs_a_value_b(k) - lgs_x_value(lgs_a_index(k)) * incrbb(nr)
    end do
  end do
  lgs_x_value_b = 0.0

end subroutine sor_sprs_grad_local

!-------------------------------------------------------------------------------
!> Solution of a system of linear equations Ax=b with tridiagonal matrix A.
!! @param[in]  a0       a0(j) is element A_(j,j-1) of Matrix A
!! @param[in]  a1       a1(j) is element A_(j,j)   of Matrix A
!! @param[in]  a2       a2(j) is element A_(j,j+1) of Matrix A
!! @param[in]  b        inhomogeneity vector
!! @param[in]  nrows    size of matrix A (indices run from 0 (!!!) to nrows)
!! @param[out] x        Solution vector.
!<------------------------------------------------------------------------------
subroutine tri_sle_grad_local(a0, a0_b, a1, a1_b, a2, a2_b, x, x_b, b, b_b, nrows)

use sico_maths_m, only : tri_sle
!use sico_maths_m, only : tri_sle_local
implicit none

integer(i4b),             intent(in)    :: nrows
real(dp), dimension(0:*), intent(in)    :: a0, a2
real(dp), dimension(0:*), intent(inout) :: a1, b

real(dp), dimension(0:*), intent(out)   :: x


real(dp), dimension(0:*), intent(inout) :: a0_b, a2_b
real(dp), dimension(0:*), intent(inout) :: a1_b, b_b
real(dp), dimension(0:*), intent(inout) :: x_b
real(dp), dimension(0:nrows)            :: a0T, a1T, a2T
real(dp), dimension(0:nrows)            :: incrbb
integer(i4b) :: i
        a0T(0) = 0.0   
        a0T(1:nrows) = a2(0:nrows-1)
        a1T = a1(0:nrows)
        a2T(0:nrows-1) = a0(1:nrows)
        a2T(nrows) = 0.0
        call tri_sle(a0T, a1T, a2T, incrbb, x_b, nrows)
        DO i=0,nrows
          b_b(i) = b_b(i) + incrbb(i)
        ENDDO
        call tri_sle(a0, a1, a2, x, b, nrows)
        DO i=0,nrows
          a0_b(i) = a0_b(i) - x(i-1)*incrbb(i)
        ENDDO
        DO i=0,nrows
          a1_b(i) = a1_b(i) - x(i)*incrbb(i)
        ENDDO
        DO i=0,nrows
          a2_b(i) = a2_b(i) - x(i+1)*incrbb(i)
        ENDDO
        DO i=0,nrows
          x_b(i) = 0
        ENDDO
end subroutine tri_sle_grad_local

!-------------------------------------------------------------------------------
#ifdef BUILD_LIS
subroutine sico_lis_solver_grad_local(nmax, nnz, &
                           lgs_a_ptr, lgs_a_index, &
                           lgs_a_value, lgs_a_value_b, lgs_b_value,& 
                           lgs_b_value_b, lgs_x_value, lgs_x_value_b)

use sico_maths_m, only :sico_lis_solver
implicit none

integer(i4b)                                 :: nc, nr
integer(i4b),                     intent(in) :: nmax
integer(i4b),                     intent(in) :: nnz
integer(i4b), dimension(nmax+1),  intent(in) :: lgs_a_ptr
integer(i4b), dimension(nnz),  intent(in) :: lgs_a_index

real(dp),     dimension(nnz),  intent(in) :: lgs_a_value
real(dp),     dimension(nmax),    intent(in) :: lgs_b_value
real(dp),     dimension(nmax), intent(inout) :: lgs_x_value

real(dp),     dimension(nnz),  intent(inout) :: lgs_a_value_b
real(dp),     dimension(nmax), intent(inout) :: lgs_b_value_b
real(dp),     dimension(nmax), intent(inout) :: lgs_x_value_b

integer(i4b), dimension(nmax+1)              :: lgs_aT_ptr
integer(i4b), dimension(nnz)              :: lgs_aT_index
real(dp),     dimension(nnz)              :: lgs_aT_value
real(dp),     dimension(nmax)                :: incrbb

  call transpose_csr_no_diagonal(lgs_a_value, lgs_a_index, lgs_a_ptr, &
                    nnz, nmax, &
                    lgs_aT_value, lgs_aT_index, lgs_aT_ptr)
  call sico_lis_solver(nmax, nnz, &
                           lgs_aT_ptr, lgs_aT_index, &
                           lgs_aT_value, lgs_x_value_b, incrbb)
  DO nr=1,nmax
    lgs_b_value_b(nr) = lgs_b_value_b(nr) + incrbb(nr)
  ENDDO

  call sico_lis_solver(nmax, nnz, &
                           lgs_a_ptr, lgs_a_index, &
                           lgs_a_value, lgs_b_value, lgs_x_value)

  do nr=1, nmax
    do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
      lgs_a_value_b(nc) = lgs_a_value_b(nc) - lgs_x_value(lgs_a_index(nc)) * incrbb(nr)
    end do
  end do
  lgs_x_value_b = 0.0
end subroutine sico_lis_solver_grad_local
#endif
end module sico_maths_m_grad
!
